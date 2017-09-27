#include <stdlib.h>
#include "mmpriv.h"
#include "kalloc.h"

void mm_select_sub_multi(void *km, float pri_ratio, float pri1, float pri2, int max_gap_ref, int min_diff, int best_n, int n_segs, const int *qlens, int *n_, mm_reg1_t *r)
{
	if (pri_ratio > 0.0f && *n_ > 0) {
		int i, k, n = *n_, n_2nd = 0;
		int max_dist = n_segs == 2? qlens[0] + qlens[1] + max_gap_ref : 0;
		for (i = k = 0; i < n; ++i) {
			int to_keep = 0;
			if (r[i].parent == i) { // primary
				to_keep = 1;
			} else if (r[i].score + min_diff >= r[r[i].parent].score) {
				to_keep = 1;
			} else {
				mm_reg1_t *p = &r[r[i].parent], *q = &r[i];
				if (p->rev == q->rev && p->rid == q->rid && q->re - p->rs < max_dist && p->re - q->rs < max_dist) { // child and parent are close on the ref
					if (q->score >= p->score * pri1)
						to_keep = 1;
				} else {
					int is_par_both = (n_segs == 2 && p->qs < qlens[0] && p->qe > qlens[0]);
					int is_chi_both = (n_segs == 2 && q->qs < qlens[0] && q->qe > qlens[0]);
					if (is_chi_both || is_chi_both == is_par_both) {
						if (q->score >= p->score * pri_ratio)
							to_keep = 1;
					} else { // the remaining case: is_chi_both == 0 && is_par_both == 1
						if (q->score >= p->score * pri2)
							to_keep = 1;
					}
				}
			}
			if (to_keep && r[i].parent != i) {
				if (n_2nd++ >= best_n) to_keep = 0; // don't keep if there are too many secondary hits
			}
			if (to_keep) r[k++] = r[i];
			else if (r[i].p) free(r[i].p);
		}
		if (k != n) mm_sync_regs(km, k, r); // removing hits requires sync()
		*n_ = k;
	}
}

#include "ksort.h"

typedef struct {
	int s, rev;
	uint64_t key;
	mm_reg1_t *r;
} pair_arr_t;

#define sort_key_pair(a) ((a).key)
KRADIX_SORT_INIT(pair, pair_arr_t, sort_key_pair, 8)

void mm_pair(void *km, int max_gap_ref, int pe_bonus, const int *qlens, int *n_regs, mm_reg1_t **regs)
{
	int i, j, s, n, last[2], dp_thres, segs = 0, n_pairs = 0, max_idx[2];
	int64_t max, max2;
	pair_arr_t *a;

	a = (pair_arr_t*)kmalloc(km, (n_regs[0] + n_regs[1]) * sizeof(pair_arr_t));
	for (s = n = 0, dp_thres = 0; s < 2; ++s) {
		int max = 0;
		for (i = 0; i < n_regs[s]; ++i) {
			a[n].s = s;
			a[n].r = &regs[s][i];
			a[n].rev = a[n].r->rev;
			a[n].key = (uint64_t)a[n].r->rid << 32 | a[n].r->rs<<1 | (s^a[n].rev);
			max = max > a[n].r->p->dp_max? max : a[n].r->p->dp_max;
			++n;
			segs |= 1<<s;
		}
		dp_thres += max;
	}
	if (segs != 3) {
		kfree(km, a); // only one end is mapped
		return;
	}
	dp_thres -= pe_bonus;
	if (dp_thres < 0) dp_thres = 0;
	radix_sort_pair(a, a + n);

	max = max2 = -1;
	max_idx[0] = max_idx[1] = -1;
	last[0] = last[1] = -1;
	for (i = 0; i < n; ++i) {
		if (a[i].key & 1) { // reverse first read or forward second read
			mm_reg1_t *q, *r;
			if (last[a[i].rev] < 0) continue;
			r = a[i].r;
			q = a[last[a[i].rev]].r;
			if (r->rid != q->rid || r->rs - q->re > max_gap_ref) continue;
			for (j = last[a[i].rev]; j >= 0; --j) {
				int sum;
				if (a[j].rev != a[i].rev || a[j].s == a[i].s) continue;
				q = a[j].r;
				if (r->rid != q->rid || r->rs - q->re > max_gap_ref) break;
				if (r->p->dp_max + q->p->dp_max < dp_thres) continue;
				++n_pairs;
				sum = r->p->dp_max + q->p->dp_max;
				if (sum > max) {
					max_idx[a[j].s] = j, max_idx[a[i].s] = i;
					max2 = max;
					max = sum;
				} else if (sum > max2)
					max2 = sum;
			}
		} else { // forward first read or reverse second read
			last[a[i].rev] = i;
		}
	}

	if (max > 0) {
		for (s = 0; s < 2; ++s) {
			mm_reg1_t *r = a[max_idx[s]].r;
			if (r->id != r->parent) {
				mm_reg1_t *p = &regs[s][r->parent];
				for (i = 0; i < n_regs[s]; ++i)
					if (regs[s][i].parent == p->id)
						regs[s][i].parent = r->id;
			}
		}
	}

	kfree(km, a);
}
