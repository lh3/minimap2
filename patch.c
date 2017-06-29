#include <string.h>
#include <assert.h>
#include "mmpriv.h"
#include "kalloc.h"

static int mm_squeeze_a_core(void *km, int n_regs, mm_reg1_t *regs, mm128_t *a, const uint64_t *aux)
{ // squeeze out regions in a[] that are not referenced by regs[]
	int i, as = 0;
	for (i = 0; i < n_regs; ++i) {
		mm_reg1_t *r = &regs[i];
		if (r->as != as) {
			memmove(&a[as], &a[r->as], r->cnt * 16);
			r->as = as;
		}
		as += r->cnt;
	}
	return as;
}

void mm_join_long(void *km, const mm_mapopt_t *opt, const mm_idx_t *mi, int qlen, int n_regs, mm_reg1_t *regs, mm128_t *a)
{
	int i, n_aux;
	uint64_t *aux;

	if (n_regs < 2) return; // nothing to join
	aux = (uint64_t*)kmalloc(km, n_regs * 8);
	for (i = n_aux = 0; i < n_regs; ++i)
		aux[n_aux++] = (uint64_t)regs[i].as << 32 | i;
	assert(n_aux == n_regs); // TODO: may be relaxed in future for other use cases
	radix_sort_64(aux, aux + n_aux);
	mm_squeeze_a_core(km, n_regs, regs, a, aux);

	for (i = n_aux - 1; i >= 1; --i) {
		mm_reg1_t *r0 = &regs[(int32_t)aux[i-1]], *r1 = &regs[(int32_t)aux[i]];
		mm128_t *a0e, *a1s;
		int max_gap, min_gap;

		// test
		if (r0->as + r0->cnt != r1->as) continue; // not adjacent in a[]
		if (r0->rid != r1->rid || r0->rev != r1->rev) continue; // make sure on the same target and strand
		if (r0->score < opt->min_join_flank_sc || r1->score < opt->min_join_flank_sc) continue; // require good flanking chains
		a0e = &a[r0->as + r0->cnt - 1];
		a1s = &a[r1->as];
		if (a1s->x <= a0e->x || (int32_t)a1s->y <= (int32_t)a0e->y) continue; // keep colinearity
		max_gap = min_gap = (int32_t)a1s->y - (int32_t)a0e->y;
		max_gap = max_gap > a1s->x - a0e->x? max_gap : a1s->x - a0e->x;
		min_gap = min_gap < a1s->x - a0e->x? min_gap : a1s->x - a0e->x;
		if (max_gap > opt->max_join_long || min_gap > opt->max_join_short) continue;

		// all conditions satisfied; join
		a[r1->as].y |= 1ULL<<40;
		r0->cnt += r1->cnt, r0->score += r1->score;
		mm_reg_set_coor(r0, qlen, a);
		r1->cnt = 0;
	}
	kfree(km, aux);
}
