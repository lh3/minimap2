#include <string.h>
#include <math.h>
#include "mmpriv.h"
#include "kalloc.h"

mm_reg1_t *mm_gen_regs(void *km, int qlen, int n_u, uint64_t *u, mm128_t *a) // convert chains to hits
{
	mm128_t *z, tmp;
	mm_reg1_t *r;
	int i, k;

	if (n_u == 0) return 0;

	// sort by score
	z = (mm128_t*)kmalloc(km, n_u * 16);
	for (i = k = 0; i < n_u; ++i) {
		z[i].x = u[i] >> 32;
		z[i].y = (uint64_t)k << 32 | (int32_t)u[i];
		k += (int32_t)u[i];
	}
	radix_sort_128x(z, z + n_u);
	for (i = 0; i < n_u>>1; ++i) // reverse, s.t. larger score first
		tmp = z[i], z[i] = z[n_u-1-i], z[n_u-1-i] = tmp;

	// populate r[]
	r = (mm_reg1_t*)calloc(n_u, sizeof(mm_reg1_t));
	for (i = 0; i < n_u; ++i) {
		mm_reg1_t *ri = &r[i];
		ri->id = i;
		ri->parent = MM_PARENT_UNSET;
		ri->score = z[i].x;
		ri->cnt = (int32_t)z[i].y;
		ri->as = z[i].y >> 32;
		mm_reg_set_coor(ri, qlen, a);
	}
	kfree(km, z);
	return r;
}

void mm_split_reg(mm_reg1_t *r, mm_reg1_t *r2, int n, int qlen, mm128_t *a)
{
	if (n <= 0 || n >= r->cnt) return;
	*r2 = *r;
	r2->id = -1;
	r2->sam_pri = 0;
	r2->p = 0;
	r2->cnt = r->cnt - n;
	r2->score = (int32_t)(r->score * ((float)r2->cnt / r->cnt) + .499);
	r2->as = r->as + n;
	if (r->parent == r->id) r2->parent = MM_PARENT_TMP_PRI;
	mm_reg_set_coor(r2, qlen, a);
	r->cnt -= r2->cnt;
	r->score -= r2->score;
	mm_reg_set_coor(r, qlen, a);
	r->split |= 1, r2->split |= 2;
}

void mm_set_parent(void *km, float mask_level, int n, mm_reg1_t *r) // and compute mm_reg1_t::subsc
{
	int i, j, k, *w;
	if (n <= 0) return;
	w = (int*)kmalloc(km, n * sizeof(int));
	w[0] = 0, r[0].parent = 0;
	for (i = 1, k = 1; i < n; ++i) {
		int si = r[i].qs, ei = r[i].qe;
		for (j = 0; j < k; ++j) {
			int sj = r[w[j]].qs, ej = r[w[j]].qe;
			int min = ej - sj < ei - si? ej - sj : ei - si;
			int ol = si < sj? (ei < sj? 0 : ei < ej? ei - sj : ej - sj) : (ej < si? 0 : ej < ei? ej - si : ei - si);
			if (ol > mask_level * min) {
				r[i].parent = r[w[j]].parent;
				if (r[w[j]].subsc < r[i].score)
					r[w[j]].subsc = r[i].score;
				break;
			}
		}
		if (j == k) w[k++] = i, r[i].parent = i;
	}
	kfree(km, w);
}

void mm_update_parent(void *km, float mask_level, int n, mm_reg1_t *r) // due to changes to r.{qs,qe} after DP extension
{
	int i, j, k, *w, n_pri = 0;
	if (n <= 0) return;
	for (i = 0; i < n; ++i)
		if (r[i].id == r[i].parent) ++n_pri;
	if (n_pri <= 1) return;
	w = (int*)kmalloc(km, n_pri * sizeof(int));
	for (i = j = 0; i < n; ++i) // find the first primary
		if (r[i].id == r[i].parent) break;
	for (w[0] = i, i = i + 1, k = 1; i < n; ++i) {
		int si = r[i].qs, ei = r[i].qe;
		if (r[i].id != r[i].parent && r[i].parent >= 0) {
			r[i].parent = r[r[i].parent].parent;
			continue;
		}
		for (j = 0; j < k; ++j) {
			int sj = r[w[j]].qs, ej = r[w[j]].qe;
			int min = ej - sj < ei - si? ej - sj : ei - si;
			int ol = si < sj? (ei < sj? 0 : ei < ej? ei - sj : ej - sj) : (ej < si? 0 : ej < ei? ej - si : ei - si);
			if (ol > mask_level * min) {
				r[i].parent = r[w[j]].parent;
				if (r[w[j]].subsc < r[i].score)
					r[w[j]].subsc = r[i].score;
				break;
			}
		}
		if (j == k) w[k++] = i;
	}
	kfree(km, w);
}

void mm_sync_regs(void *km, int n_regs, mm_reg1_t *regs) // keep mm_reg1_t::{id,parent} in sync; also reset id
{
	int *tmp, i, max_id = -1, n_tmp, n_pri;
	if (n_regs <= 0) return;
	for (i = 0; i < n_regs; ++i)
		max_id = max_id > regs[i].id? max_id : regs[i].id;
	n_tmp = max_id + 1;
	tmp = (int*)kmalloc(km, n_tmp * sizeof(int));
	for (i = 0; i < n_tmp; ++i) tmp[i] = -1;
	for (i = 0; i < n_regs; ++i)
		if (regs[i].id >= 0) tmp[regs[i].id] = i;
	for (i = n_pri = 0; i < n_regs; ++i) {
		mm_reg1_t *r = &regs[i];
		r->id = i;
		if (r->parent == MM_PARENT_TMP_PRI)
			r->parent = i;
		else if (r->parent >= 0 && tmp[r->parent] >= 0)
			r->parent = tmp[r->parent];
		else r->parent = MM_PARENT_UNSET;
		if (r->id == r->parent) {
			if (n_pri++ == 0) r->sam_pri = 1;
		}
	}
	kfree(km, tmp);
}

void mm_select_sub(void *km, float mask_level, float pri_ratio, int best_n, int *n_, mm_reg1_t *r)
{
	if (pri_ratio > 0.0f && *n_ > 0) {
		int i, k, n = *n_, n_2nd = 0;
		for (i = k = 0; i < n; ++i)
			if (r[i].parent == i || r[i].parent < 0) r[k++] = r[i]; // NB: r[i].parent may be -1 if its parent has been filtered
			else if (r[i].score >= r[r[i].parent].score * pri_ratio && n_2nd++ < best_n) r[k++] = r[i];
			else if (r[i].p) free(r[i].p);
		if (k != n) mm_sync_regs(km, k, r); // removing hits requires sync()
		*n_ = k;
	}
}

void mm_filter_regs(void *km, const mm_mapopt_t *opt, int *n_regs, mm_reg1_t *regs)
{ // NB: after this call, mm_reg1_t::parent can be -1 if its parent filtered out
	int i, k;
	for (i = k = 0; i < *n_regs; ++i) {
		mm_reg1_t *r = &regs[i];
		int flt = 0;
		if (r->cnt < opt->min_cnt) flt = 1;
		else {
			int blen = r->qe - r->qs < r->re - r->rs? r->qe - r->qs : r->re - r->rs;
			if (r->score < blen * opt->min_seedcov_ratio) flt = 1;
		}
		if (r->p) {
			if (r->p->blen - r->p->n_ambi - r->p->n_diff < opt->min_chain_score) flt = 1;
			else if (r->p->dp_max < opt->min_dp_max) flt = 1;
			if (flt) free(r->p);
		}
		if (!flt) {
			if (k < i) regs[k++] = regs[i];
			else ++k;
		}
	}
	*n_regs = k;
	mm_sync_regs(km, k, regs);
}

int mm_squeeze_a(void *km, int n_regs, mm_reg1_t *regs, mm128_t *a)
{ // squeeze out regions in a[] that are not referenced by regs[]
	int i, as = 0;
	uint64_t *aux;
	aux = (uint64_t*)kmalloc(km, n_regs * 8);
	for (i = 0; i < n_regs; ++i)
		aux[i] = (uint64_t)regs[i].as << 32 | i;
	radix_sort_64(aux, aux + n_regs);
	for (i = 0; i < n_regs; ++i) {
		mm_reg1_t *r = &regs[(int32_t)aux[i]];
		if (r->as != as) {
			memmove(&a[as], &a[r->as], r->cnt * 16);
			r->as = as;
		}
		as += r->cnt;
	}
	kfree(km, aux);
	return as;
}

void mm_join_long(void *km, const mm_mapopt_t *opt, int qlen, int *n_regs_, mm_reg1_t *regs, mm128_t *a)
{
	int i, n_aux, n_regs = *n_regs_, n_drop = 0;
	uint64_t *aux;

	if (n_regs < 2) return; // nothing to join
	mm_squeeze_a(km, n_regs, regs, a);

	aux = (uint64_t*)kmalloc(km, n_regs * 8);
	for (i = n_aux = 0; i < n_regs; ++i)
		if (regs[i].parent == i || regs[i].parent < 0)
			aux[n_aux++] = (uint64_t)regs[i].as << 32 | i;
	radix_sort_64(aux, aux + n_aux);

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
		r1->parent = r0->id;
		++n_drop;
	}
	kfree(km, aux);

	if (n_drop > 0) { // then fix the hits hierarchy
		for (i = 0; i < n_regs; ++i) { // adjust the mm_reg1_t::parent
			mm_reg1_t *r = &regs[i];
			if (r->parent >= 0 && r->id != r->parent) { // fix for secondary hits only
				if (regs[r->parent].parent >= 0 && regs[r->parent].parent != r->parent)
					r->parent = regs[r->parent].parent;
			}
		}
		mm_filter_regs(km, opt, n_regs_, regs);
	}
}

void mm_set_mapq(int n_regs, mm_reg1_t *regs)
{
	int i;
	for (i = 0; i < n_regs; ++i) {
		mm_reg1_t *r = &regs[i];
		if (r->parent == r->id) {
			int mapq;
			mapq = (int)(30.0 * (1. - (float)r->subsc / r->score) * logf(r->score));
			mapq = mapq > 0? mapq : 0;
			r->mapq = mapq < 60? mapq : 60;
		} else r->mapq = 0;
	}
}
