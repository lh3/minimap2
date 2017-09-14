#include <string.h>
#include <math.h>
#include "mmpriv.h"
#include "kalloc.h"

static inline void mm_cal_fuzzy_len(mm_reg1_t *r, const mm128_t *a)
{
	int i;
	r->fuzzy_mlen = r->fuzzy_blen = 0;
	if (r->cnt <= 0) return;
	r->fuzzy_mlen = r->fuzzy_blen = a[r->as].y>>32&0xff;
	for (i = r->as + 1; i < r->as + r->cnt; ++i) {
		int span = a[i].y>>32&0xff;
		int tl = (int32_t)a[i].x - (int32_t)a[i-1].x;
		int ql = (int32_t)a[i].y - (int32_t)a[i-1].y;
		r->fuzzy_blen += tl > ql? tl : ql;
		r->fuzzy_mlen += tl > span && ql > span? span : tl < ql? tl : ql;
	}
}

static inline void mm_reg_set_coor(mm_reg1_t *r, int32_t qlen, const mm128_t *a)
{ // NB: r->as and r->cnt MUST BE set correctly for this function to work
	int32_t k = r->as, q_span = (int32_t)(a[k].y>>32&0xff);
	r->rev = a[k].x>>63;
	r->rid = a[k].x<<1>>33;
	r->rs = (int32_t)a[k].x + 1 > q_span? (int32_t)a[k].x + 1 - q_span : 0; // NB: target span may be shorter, so this test is necessary
	r->re = (int32_t)a[k + r->cnt - 1].x + 1;
	if (!r->rev) {
		r->qs = (int32_t)a[k].y + 1 - q_span;
		r->qe = (int32_t)a[k + r->cnt - 1].y + 1;
	} else {
		r->qs = qlen - ((int32_t)a[k + r->cnt - 1].y + 1);
		r->qe = qlen - ((int32_t)a[k].y + 1 - q_span);
	}
	mm_cal_fuzzy_len(r, a);
}

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

void mm_set_parent(void *km, float mask_level, int n, mm_reg1_t *r, int sub_diff) // and compute mm_reg1_t::subsc
{
	int i, j, k, *w;
	if (n <= 0) return;
	for (i = 0; i < n; ++i) r[i].id = i;
	w = (int*)kmalloc(km, n * sizeof(int));
	w[0] = 0, r[0].parent = 0;
	for (i = 1, k = 1; i < n; ++i) {
		mm_reg1_t *ri = &r[i];
		int si = ri->qs, ei = ri->qe;
		for (j = 0; j < k; ++j) {
			mm_reg1_t *rp = &r[w[j]];
			int sj = rp->qs, ej = rp->qe;
			int min = ej - sj < ei - si? ej - sj : ei - si;
			int ol = si < sj? (ei < sj? 0 : ei < ej? ei - sj : ej - sj) : (ej < si? 0 : ej < ei? ej - si : ei - si);
			if (ol > mask_level * min) {
				int cnt_sub = 0;
				ri->parent = rp->parent;
				rp->subsc = rp->subsc > ri->score? rp->subsc : ri->score;
				if (ri->cnt >= rp->cnt) cnt_sub = 1;
				if (rp->p && ri->p) {
					rp->p->dp_max2 = rp->p->dp_max2 > ri->p->dp_max? rp->p->dp_max2 : ri->p->dp_max;
					if (rp->p->dp_max - ri->p->dp_max <= sub_diff) cnt_sub = 1;
				}
				if (cnt_sub) ++rp->n_sub;
				break;
			}
		}
		if (j == k) w[k++] = i, ri->parent = i, ri->n_sub = 0;
	}
	kfree(km, w);
}

void mm_hit_sort_by_dp(void *km, int *n_regs, mm_reg1_t *r)
{
	int32_t i, n_aux, n = *n_regs;
	uint64_t *aux;
	mm_reg1_t *t;

	if (n <= 1) return;
	aux = (uint64_t*)kmalloc(km, n * 8);
	t = (mm_reg1_t*)kmalloc(km, n * sizeof(mm_reg1_t));
	for (i = n_aux = 0; i < n; ++i) {
		if (r[i].inv || r[i].cnt > 0) { // squeeze out elements with cnt==0 (soft deleted)
			assert(r[i].p);
			aux[n_aux++] = (uint64_t)r[i].p->dp_max << 32 | i;
		} else if (r[i].p) {
			free(r[i].p);
			r[i].p = 0;
		}
	}
	radix_sort_64(aux, aux + n_aux);
	for (i = n_aux - 1; i >= 0; --i)
		t[n_aux - 1 - i] = r[(int32_t)aux[i]];
	memcpy(r, t, sizeof(mm_reg1_t) * n_aux);
	*n_regs = n_aux;
	kfree(km, aux);
	kfree(km, t);
}

int mm_set_sam_pri(int n, mm_reg1_t *r)
{
	int i, n_pri = 0;
	for (i = 0; i < n; ++i)
		if (r[i].id == r[i].parent) {
			++n_pri;
			r[i].sam_pri = (n_pri == 1);
		} else r[i].sam_pri = 0;
	return n_pri;
}

void mm_sync_regs(void *km, int n_regs, mm_reg1_t *regs) // keep mm_reg1_t::{id,parent} in sync; also reset id
{
	int *tmp, i, max_id = -1, n_tmp;
	if (n_regs <= 0) return;
	for (i = 0; i < n_regs; ++i) // NB: doesn't work if mm_reg1_t::id is negative
		max_id = max_id > regs[i].id? max_id : regs[i].id;
	n_tmp = max_id + 1;
	tmp = (int*)kmalloc(km, n_tmp * sizeof(int));
	for (i = 0; i < n_tmp; ++i) tmp[i] = -1;
	for (i = 0; i < n_regs; ++i)
		if (regs[i].id >= 0) tmp[regs[i].id] = i;
	for (i = 0; i < n_regs; ++i) {
		mm_reg1_t *r = &regs[i];
		r->id = i;
		if (r->parent == MM_PARENT_TMP_PRI)
			r->parent = i;
		else if (r->parent >= 0 && tmp[r->parent] >= 0)
			r->parent = tmp[r->parent];
		else r->parent = MM_PARENT_UNSET;
	}
	kfree(km, tmp);
	mm_set_sam_pri(n_regs, regs);
}

void mm_select_sub(void *km, float mask_level, float pri_ratio, int min_diff, int best_n, int *n_, mm_reg1_t *r)
{
	if (pri_ratio > 0.0f && *n_ > 0) {
		int i, k, n = *n_, n_2nd = 0;
		for (i = k = 0; i < n; ++i)
			if (r[i].parent == i) r[k++] = r[i];
			else if ((r[i].score >= r[r[i].parent].score * pri_ratio || r[i].score + min_diff >= r[r[i].parent].score) && n_2nd++ < best_n)
				r[k++] = r[i];
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
		if (!r->inv && r->cnt < opt->min_cnt) flt = 1;
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
		int max_gap, min_gap, sc_thres;

		// test
		if (r0->as + r0->cnt != r1->as) continue; // not adjacent in a[]
		if (r0->rid != r1->rid || r0->rev != r1->rev) continue; // make sure on the same target and strand
		a0e = &a[r0->as + r0->cnt - 1];
		a1s = &a[r1->as];
		if (a1s->x <= a0e->x || (int32_t)a1s->y <= (int32_t)a0e->y) continue; // keep colinearity
		max_gap = min_gap = (int32_t)a1s->y - (int32_t)a0e->y;
		max_gap = max_gap > a1s->x - a0e->x? max_gap : a1s->x - a0e->x;
		min_gap = min_gap < a1s->x - a0e->x? min_gap : a1s->x - a0e->x;
		if (max_gap > opt->max_join_long || min_gap > opt->max_join_short) continue;
		sc_thres = (int)((float)opt->min_join_flank_sc / opt->max_join_long * max_gap + .499);
		if (r0->score < sc_thres || r1->score < sc_thres) continue; // require good flanking chains
		if (r0->re - r0->rs < max_gap>>1 || r0->qe - r0->qs < max_gap>>1) continue; // require enough flanking length
		if (r1->re - r1->rs < max_gap>>1 || r1->qe - r1->qs < max_gap>>1) continue;

		// all conditions satisfied; join
		a[r1->as].y |= MM_SEED_LONG_JOIN;
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
		mm_sync_regs(km, *n_regs_, regs);
	}
}

void mm_set_mapq(int n_regs, mm_reg1_t *regs, int min_chain_sc, int match_sc, int rep_len)
{
	static const float q_coef = 40.0f;
	int i;
	for (i = 0; i < n_regs; ++i) {
		mm_reg1_t *r = &regs[i];
		if (r->inv) {
			r->mapq = 0;
		} else if (r->parent == r->id) {
			int mapq, subsc;
			float pen_s1 = r->score > 100? 1.0f : 0.01f * r->score;
			float pen_cm = r->cnt > 10? 1.0f : 0.1f * r->cnt;
			if (r->score <= 100 && rep_len > 0) {
				pen_s1 = 0.01f * (r->score - rep_len);
				pen_s1 = pen_s1 > 0.1f? pen_s1 : 0.1f;
			}
			pen_cm = pen_s1 < pen_cm? pen_s1 : pen_cm;
			subsc = r->subsc > min_chain_sc? r->subsc : min_chain_sc;
			if (r->p && r->p->dp_max2 > 0 && r->p->dp_max > 0) {
				float identity = (float)(r->p->blen - r->p->n_diff - r->p->n_ambi) / (r->p->blen - r->p->n_ambi);
				int mapq_alt = (int)(6.02f * identity * identity * (r->p->dp_max - r->p->dp_max2) / match_sc + .499f); // BWA-MEM like mapQ, mostly for short reads
				mapq = (int)(identity * pen_cm * q_coef * (1. - (float)r->p->dp_max2 * subsc / r->p->dp_max / r->score) * logf(r->score)); // more for long reads
				mapq = mapq < mapq_alt? mapq : mapq_alt; // in case the long-read heuristic fails
			} else mapq = (int)(pen_cm * q_coef * (1. - (float)subsc / r->score) * logf(r->score));
			mapq -= (int)(4.343f * logf(r->n_sub + 1) + .499f);
			mapq = mapq > 0? mapq : 0;
			r->mapq = mapq < 60? mapq : 60;
		} else r->mapq = 0;
	}
}
