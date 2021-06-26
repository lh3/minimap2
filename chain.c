/* The MIT License

Copyright (c) 2018-     Dana-Farber Cancer Institute
              2017-2018 Broad Institute, Inc.

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
Modified Copyright (C) 2021 Intel Corporation
   Contacts: Saurabh Kalikar <saurabh.kalikar@intel.com>; 
	Vasimuddin Md <vasimuddin.md@intel.com>; Sanchit Misra <sanchit.misra@intel.com>; 
	Chirag Jain <chirag@iisc.ac.in>; Heng Li <hli@jimmy.harvard.edu>
*/
#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include "minimap.h"
#include "mmpriv.h"
#include "kalloc.h"
#include "parallel_chaining_32_bit.h"
static const char LogTable256[256] = {
#define LT(n) n, n, n, n, n, n, n, n, n, n, n, n, n, n, n, n
	-1, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3,
	LT(4), LT(5), LT(5), LT(6), LT(6), LT(6), LT(6),
	LT(7), LT(7), LT(7), LT(7), LT(7), LT(7), LT(7), LT(7)
};

static inline int ilog2_32(uint32_t v)
{
	uint32_t t, tt;
	if ((tt = v>>16)) return (t = tt>>8) ? 24 + LogTable256[t] : 16 + LogTable256[tt];
	return (t = v>>8) ? 8 + LogTable256[t] : LogTable256[v];
}


mm128_t *mm_chain_dp(int max_dist_x, int max_dist_y, int bw, int max_skip, int max_iter, int min_cnt, int min_sc, float gap_scale, int is_cdna, int n_segs, int64_t n, mm128_t *a, int *n_u_, uint64_t **_u, void *km)
{ // TODO: make sure this works when n has more than 32 bits
	int32_t k, *p, *t, *v, n_u, n_v;
	uint32_t *f;
	int64_t i, j, st = 0;
	uint64_t *u, *u2, sum_qspan = 0;
	float avg_qspan;
	mm128_t *b, *w;

	if (_u) *_u = 0, *n_u_ = 0;
	if (n == 0 || a == 0) {
		kfree(km, a);
		return 0;
	}
	f = (uint32_t*)kmalloc(km, n * 4);
	p = (int32_t*)kmalloc(km, n * 4);
	t = (int32_t*)kmalloc(km, n * 4);
	v = (int32_t*)kmalloc(km, n * 4);
	memset(t, 0, n * 4);
#if VECTORIZED_CHAINING && (__AVX512BW__ || __AVX2__)
	anchor_t* anchors = (anchor_t*)malloc(n* sizeof(anchor_t));
	for (i = 0; i < n; ++i) {
		uint64_t ri = a[i].x;
		int32_t qi = (int32_t)a[i].y, q_span = a[i].y>>32&0xff; // NB: only 8 bits of span is used!!!
		anchors[i].r = ri;
		anchors[i].q = qi;
		anchors[i].l = q_span;
	}
	num_bits_t *anchor_r, *anchor_q, *anchor_l;
	create_SoA_Anchors_32_bit(anchors, n, anchor_r, anchor_q, anchor_l);
	dp_chain obj(max_dist_x, max_dist_y, bw, max_skip, max_iter, gap_scale, is_cdna, n_segs);

	obj.mm_dp_vectorized(n, &anchors[0], anchor_r, anchor_q, anchor_l, f, p, v, max_dist_x, max_dist_y, NULL, NULL);	

	// -16 is due to extra padding at the start of arrays
	anchor_r -= 16; anchor_q -= 16; anchor_l -= 16;
	free(anchor_r); 
	free(anchor_q); 
	free(anchor_l);
	free(anchors);
#else

	for (i = 0; i < n; ++i) sum_qspan += a[i].y>>32&0xff;
	avg_qspan = (float)sum_qspan / n;

	// fill the score and backtrack arrays
	for (i = 0; i < n; ++i) {
		uint64_t ri = a[i].x;
		int64_t max_j = -1;
		int32_t qi = (int32_t)a[i].y, q_span = a[i].y>>32&0xff; // NB: only 8 bits of span is used!!!
		int32_t max_f = q_span, n_skip = 0, min_d;
		int32_t sidi = (a[i].y & MM_SEED_SEG_MASK) >> MM_SEED_SEG_SHIFT;
		while (st < i && ri > a[st].x + max_dist_x) ++st;
		if (i - st > max_iter) st = i - max_iter;
		for (j = i - 1; j >= st; --j) {
			int64_t dr = ri - a[j].x;
			int32_t dq = qi - (int32_t)a[j].y, dd, sc, log_dd, gap_cost;
			int32_t sidj = (a[j].y & MM_SEED_SEG_MASK) >> MM_SEED_SEG_SHIFT;
			if ((sidi == sidj && dr == 0) || dq <= 0) continue; // don't skip if an anchor is used by multiple segments; see below
			if ((sidi == sidj && dq > max_dist_y) || dq > max_dist_x) continue;
			dd = dr > dq? dr - dq : dq - dr;
			if (sidi == sidj && dd > bw) continue;
			if (n_segs > 1 && !is_cdna && sidi == sidj && dr > max_dist_y) continue;
			min_d = dq < dr? dq : dr;
			sc = min_d > q_span? q_span : dq < dr? dq : dr;
			log_dd = dd? ilog2_32(dd) : 0;
			gap_cost = 0;
			if (is_cdna || sidi != sidj) {
				int c_log, c_lin;
				c_lin = (int)(dd * .01 * avg_qspan);
				c_log = log_dd;
				if (sidi != sidj && dr == 0) ++sc; // possibly due to overlapping paired ends; give a minor bonus
				else if (dr > dq || sidi != sidj) gap_cost = c_lin < c_log? c_lin : c_log;
				else gap_cost = c_lin + (c_log>>1);
			} else gap_cost = (int)(dd * .01 * avg_qspan) + (log_dd>>1);
			sc -= (int)((double)gap_cost * gap_scale + .499);
			sc += f[j];
			if (sc > max_f) {
				max_f = sc, max_j = j;
				if (n_skip > 0) --n_skip;
			} else if (t[j] == i) {
				if (++n_skip > max_skip)
					break;
			}
			if (p[j] >= 0) t[p[j]] = i;
		}
		f[i] = max_f, p[i] = max_j;
		v[i] = max_j >= 0 && v[max_j] > max_f? v[max_j] : max_f; // v[] keeps the peak score up to i; f[] is the score ending at i, not always the peak
	}
#endif
	// find the ending positions of chains
	memset(t, 0, n * 4);
	for (i = 0; i < n; ++i)
		if (p[i] >= 0) t[p[i]] = 1;
	for (i = n_u = 0; i < n; ++i)
		if (t[i] == 0 && v[i] >= min_sc)
			++n_u;
	if (n_u == 0) {
		kfree(km, a); kfree(km, f); kfree(km, p); kfree(km, t); kfree(km, v);
		return 0;
	}
	u = (uint64_t*)kmalloc(km, n_u * 8);
	for (i = n_u = 0; i < n; ++i) {
		if (t[i] == 0 && v[i] >= min_sc) {
			j = i;
			while (j >= 0 && f[j] < v[j]) j = p[j]; // find the peak that maximizes f[]
			if (j < 0) j = i; // TODO: this should really be assert(j>=0)
			u[n_u++] = (uint64_t)f[j] << 32 | j;
		}
	}
	radix_sort_64(u, u + n_u);
	for (i = 0; i < n_u>>1; ++i) { // reverse, s.t. the highest scoring chain is the first
		uint64_t t = u[i];
		u[i] = u[n_u - i - 1], u[n_u - i - 1] = t;
	}

	// backtrack
	memset(t, 0, n * 4);
	for (i = n_v = k = 0; i < n_u; ++i) { // starting from the highest score
		int32_t n_v0 = n_v, k0 = k;
		j = (int32_t)u[i];
		do {
			v[n_v++] = j;
			t[j] = 1;
			j = p[j];
		} while (j >= 0 && t[j] == 0);
		if (j < 0) {
			if (n_v - n_v0 >= min_cnt) u[k++] = u[i]>>32<<32 | (n_v - n_v0);
		} else if ((int32_t)(u[i]>>32) - f[j] >= min_sc) {
			if (n_v - n_v0 >= min_cnt) u[k++] = ((u[i]>>32) - f[j]) << 32 | (n_v - n_v0);
		}
		if (k0 == k) n_v = n_v0; // no new chain added, reset
	}
	*n_u_ = n_u = k, *_u = u; // NB: note that u[] may not be sorted by score here

	// free temporary arrays
	kfree(km, f); kfree(km, p); kfree(km, t);

	// write the result to b[]
	b = (mm128_t*)kmalloc(km, n_v * sizeof(mm128_t));
	for (i = 0, k = 0; i < n_u; ++i) {
		int32_t k0 = k, ni = (int32_t)u[i];
		for (j = 0; j < ni; ++j)
			b[k] = a[v[k0 + (ni - j - 1)]], ++k;
	}
	kfree(km, v);

	// sort u[] and a[] by a[].x, such that adjacent chains may be joined (required by mm_join_long)
	w = (mm128_t*)kmalloc(km, n_u * sizeof(mm128_t));
	for (i = k = 0; i < n_u; ++i) {
		w[i].x = b[k].x, w[i].y = (uint64_t)k<<32|i;
		k += (int32_t)u[i];
	}
	radix_sort_128x(w, w + n_u);
	u2 = (uint64_t*)kmalloc(km, n_u * 8);
	for (i = k = 0; i < n_u; ++i) {
		int32_t j = (int32_t)w[i].y, n = (int32_t)u[j];
		u2[i] = u[j];
		memcpy(&a[k], &b[w[i].y>>32], n * sizeof(mm128_t));
		k += n;
	}
	if (n_u) memcpy(u, u2, n_u * 8);
	if (k) memcpy(b, a, k * sizeof(mm128_t)); // write _a_ to _b_ and deallocate _a_ because _a_ is oversized, sometimes a lot
	kfree(km, a); kfree(km, w); kfree(km, u2);
	return b;
}
