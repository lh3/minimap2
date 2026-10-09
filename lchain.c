#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include "mmpriv.h"
#include "kalloc.h"
#include "krmq.h"

static int64_t mg_chain_bk_end(int32_t max_drop, const mm128_t *z, const int32_t *f, const int64_t *p, int32_t *t, int64_t k)
{
	int64_t i = z[k].y, end_i = -1, max_i = i;
	int32_t max_s = 0;
	if (i < 0 || t[i] != 0) return i;
	do {
		int32_t s;
		t[i] = 2;
		end_i = i = p[i];
		s = i < 0? z[k].x : (int32_t)z[k].x - f[i];
		if (s > max_s) max_s = s, max_i = i;
		else if (max_s - s > max_drop) break;
	} while (i >= 0 && t[i] == 0);
	for (i = z[k].y; i >= 0 && i != end_i; i = p[i]) // reset modified t[]
		t[i] = 0;
	return max_i;
}

uint64_t *mg_chain_backtrack(void *km, int64_t n, const int32_t *f, const int64_t *p, int32_t *v, int32_t *t, int32_t min_cnt, int32_t min_sc, int32_t max_drop, int32_t *n_u_, int32_t *n_v_)
{
	mm128_t *z;
	uint64_t *u;
	int64_t i, k, n_z, n_v;
	int32_t n_u;

	*n_u_ = *n_v_ = 0;
	for (i = 0, n_z = 0; i < n; ++i) // precompute n_z
		if (f[i] >= min_sc) ++n_z;
	if (n_z == 0) return 0;
	z = Kmalloc(km, mm128_t, n_z);
	for (i = 0, k = 0; i < n; ++i) // populate z[]
		if (f[i] >= min_sc) z[k].x = f[i], z[k++].y = i;
	radix_sort_128x(z, z + n_z);

	memset(t, 0, n * 4);
	for (k = n_z - 1, n_v = n_u = 0; k >= 0; --k) { // precompute n_u
		if (t[z[k].y] == 0) {
			int64_t n_v0 = n_v, end_i;
			int32_t sc;
			end_i = mg_chain_bk_end(max_drop, z, f, p, t, k);
			for (i = z[k].y; i != end_i; i = p[i])
				++n_v, t[i] = 1;
			sc = i < 0? z[k].x : (int32_t)z[k].x - f[i];
			if (sc >= min_sc && n_v > n_v0 && n_v - n_v0 >= min_cnt)
				++n_u;
			else n_v = n_v0;
		}
	}
	u = Kmalloc(km, uint64_t, n_u);
	memset(t, 0, n * 4);
	for (k = n_z - 1, n_v = n_u = 0; k >= 0; --k) { // populate u[]
		if (t[z[k].y] == 0) {
			int64_t n_v0 = n_v, end_i;
			int32_t sc;
			end_i = mg_chain_bk_end(max_drop, z, f, p, t, k);
			for (i = z[k].y; i != end_i; i = p[i])
				v[n_v++] = i, t[i] = 1;
			sc = i < 0? z[k].x : (int32_t)z[k].x - f[i];
			if (sc >= min_sc && n_v > n_v0 && n_v - n_v0 >= min_cnt)
				u[n_u++] = (uint64_t)sc << 32 | (n_v - n_v0);
			else n_v = n_v0;
		}
	}
	kfree(km, z);
	assert(n_v < INT32_MAX);
	*n_u_ = n_u, *n_v_ = n_v;
	return u;
}

static mm128_t *compact_a(void *km, int32_t n_u, uint64_t *u, int32_t n_v, int32_t *v, mm128_t *a)
{
	mm128_t *b, *w;
	uint64_t *u2;
	int64_t i, j, k;

	// write the result to b[]
	b = Kmalloc(km, mm128_t, n_v);
	for (i = 0, k = 0; i < n_u; ++i) {
		int32_t k0 = k, ni = (int32_t)u[i];
		for (j = 0; j < ni; ++j)
			b[k++] = a[v[k0 + (ni - j - 1)]];
	}
	kfree(km, v);

	// sort u[] and a[] by the target position, such that adjacent chains may be joined
	w = Kmalloc(km, mm128_t, n_u);
	for (i = k = 0; i < n_u; ++i) {
		w[i].x = b[k].x, w[i].y = (uint64_t)k<<32|i;
		k += (int32_t)u[i];
	}
	radix_sort_128x(w, w + n_u);
	u2 = Kmalloc(km, uint64_t, n_u);
	for (i = k = 0; i < n_u; ++i) {
		int32_t j = (int32_t)w[i].y, n = (int32_t)u[j];
		u2[i] = u[j];
		memcpy(&a[k], &b[w[i].y>>32], n * sizeof(mm128_t));
		k += n;
	}
	memcpy(u, u2, n_u * 8);
	memcpy(b, a, k * sizeof(mm128_t)); // write _a_ to _b_ and deallocate _a_ because _a_ is oversized, sometimes a lot
	kfree(km, a); kfree(km, w); kfree(km, u2);
	return b;
}

static inline int32_t comput_sc(const mm128_t *ai, const mm128_t *aj, int32_t max_dist_x, int32_t max_dist_y, int32_t bw, float chn_pen_gap, float chn_pen_skip, int is_cdna, int n_seg)
{
	int32_t dq = (int32_t)ai->y - (int32_t)aj->y, dr, dd, dg, q_span, sc;
	int32_t sidi = (ai->y & MM_SEED_SEG_MASK) >> MM_SEED_SEG_SHIFT;
	int32_t sidj = (aj->y & MM_SEED_SEG_MASK) >> MM_SEED_SEG_SHIFT;
	if (dq <= 0 || dq > max_dist_x) return INT32_MIN;
	dr = (int32_t)(ai->x - aj->x);
	if (sidi == sidj && (dr == 0 || dq > max_dist_y)) return INT32_MIN;
	dd = dr > dq? dr - dq : dq - dr;
	if (sidi == sidj && dd > bw) return INT32_MIN;
	if (n_seg > 1 && !is_cdna && sidi == sidj && dr > max_dist_y) return INT32_MIN;
	dg = dr < dq? dr : dq;
	q_span = aj->y>>32&0xff;
	sc = q_span < dg? q_span : dg;
	if (dd || dg > q_span) {
		float lin_pen, log_pen;
		lin_pen = chn_pen_gap * (float)dd + chn_pen_skip * (float)dg;
		log_pen = dd >= 1? mg_log2(dd + 1) : 0.0f; // mg_log2() only works for dd>=2
		if (is_cdna || sidi != sidj) {
			if (sidi != sidj && dr == 0) ++sc; // possibly due to overlapping paired ends; give a minor bonus
			else if (dr > dq || sidi != sidj) sc -= (int)(lin_pen < log_pen? lin_pen : log_pen); // deletion or jump between paired ends
			else sc -= (int)(lin_pen + .5f * log_pen);
		} else sc -= (int)(lin_pen + .5f * log_pen);
	}
	return sc;
}

/* Input:
 *   a[].x: rev<<63 | tid<<32 | tpos
 *   a[].y: flags<<40 | q_span<<32 | q_pos
 * Output:
 *   n_u: #chains
 *   u[]: score<<32 | #anchors (sum of lower 32 bits of u[] is the returned length of a[])
 * input a[] is deallocated on return
 */
#ifdef MM_AVX2_DISPATCH
#include "ksw2.h"
#include <immintrin.h>
#define MM_CB 8
// The fill loop of mg_lchain_dp() for n_seg == 1 and !is_cdna (the mask is from the fast-hifi fork). Most predecessors fail the
// distance or band test of comput_sc() and leave the DP state unchanged. An AVX2 mask finds the valid predecessors
// 8 at a time; the scalar code then runs on them in the original order, so f, p, v and t stay the same.
__attribute__((target("avx2"))) static void mg_lchain_dp_fill_avx2(int max_dist_x, int max_dist_y, int bw, int max_skip, int max_iter, float chn_pen_gap, float chn_pen_skip,
								  int64_t n, const mm128_t *a, int32_t *f, int64_t *p, int32_t *v, int32_t *t, void *km)
{
	int64_t i, j, max_ii, st = 0;
	int32_t *x32 = Kmalloc(km, int32_t, n + MM_CB), *y32 = Kmalloc(km, int32_t, n + MM_CB);
	for (i = 0; i < n; ++i) x32[i + MM_CB] = (int32_t)a[i].x, y32[i + MM_CB] = (int32_t)a[i].y;
	for (i = 0; i < MM_CB; ++i) x32[i] = y32[i] = 0;
	const __m256i lane = _mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7), zero = _mm256_setzero_si256();
	const __m256i vmx = _mm256_set1_epi32(max_dist_x), vmy = _mm256_set1_epi32(max_dist_y), vbw = _mm256_set1_epi32(bw);
	for (i = 0, max_ii = -1; i < n; ++i) {
		int64_t max_j = -1, end_j, brk_j = -1;
		int32_t max_f = a[i].y>>32&0xff, n_skip = 0;
		int stop = 0;
		while (st < i && (a[i].x>>32 != a[st].x>>32 || a[i].x > a[st].x + max_dist_x)) ++st;
		if (i - st > max_iter) st = i - max_iter;
		const __m256i vxi = _mm256_set1_epi32((int32_t)a[i].x), vyi = _mm256_set1_epi32((int32_t)a[i].y);
		for (j = i - 1; j >= st && !stop; j -= 8) {
			const int64_t b = j - 7; // lanes 0..7 = anchors j-7..j; lanes below st are masked off
			__m256i vx = _mm256_loadu_si256((const __m256i*)&x32[b + MM_CB]), vy = _mm256_loadu_si256((const __m256i*)&y32[b + MM_CB]);
			__m256i dq = _mm256_sub_epi32(vyi, vy), dr = _mm256_sub_epi32(vxi, vx);
			__m256i dd = _mm256_abs_epi32(_mm256_sub_epi32(dr, dq));
			__m256i ok = _mm256_cmpgt_epi32(dq, zero);
			__m256i bad = _mm256_or_si256(_mm256_or_si256(_mm256_cmpgt_epi32(dq, vmx), _mm256_cmpgt_epi32(dq, vmy)),
										  _mm256_or_si256(_mm256_cmpeq_epi32(dr, zero), _mm256_cmpgt_epi32(dd, vbw)));
			unsigned m;
			if (b < st) bad = _mm256_or_si256(bad, _mm256_cmpgt_epi32(_mm256_set1_epi32((int32_t)(st - b)), lane));
			m = (unsigned)_mm256_movemask_ps(_mm256_castsi256_ps(_mm256_andnot_si256(bad, ok)));
			while (m) {
				int k = 31 - __builtin_clz(m); // the highest lane is the largest j: the scalar order
				int64_t jj = b + k;
				int32_t sc, dq_j = (int32_t)a[i].y - (int32_t)a[jj].y, dr_j = (int32_t)(a[i].x - a[jj].x);
				int32_t q_span = a[jj].y>>32&0xff, dg = dr_j < dq_j? dr_j : dq_j;
				m &= ~(1u << k);
				// comput_sc() <= min(q_span, dg) as its gap penalty is >= 0: if that bound cannot beat max_f,
				// the float score is not needed (the lane passed all tests, so comput_sc() is not INT32_MIN)
				if (f[jj] + (q_span < dg? q_span : dg) <= max_f) sc = INT32_MIN;
				else {
					sc = comput_sc(&a[i], &a[jj], max_dist_x, max_dist_y, bw, chn_pen_gap, chn_pen_skip, 0, 1);
					if (sc == INT32_MIN) continue;
					sc += f[jj];
				}
				if (sc > max_f) {
					max_f = sc, max_j = jj;
					if (n_skip > 0) --n_skip;
				} else if (t[jj] == (int32_t)i) {
					if (++n_skip > max_skip) { stop = 1, brk_j = jj; break; }
				}
				if (p[jj] >= 0) t[p[jj]] = i;
			}
		}
		end_j = stop? brk_j : st - 1;
		if (max_ii < 0 || a[i].x - a[max_ii].x > (int64_t)max_dist_x) { // argmax of f[st..i-1]; ties to the largest j
			int32_t max = INT32_MIN;
			max_ii = -1;
			for (j = i - 1; j >= st; --j)
				if (max < f[j]) max = f[j], max_ii = j;
		}
		if (max_ii >= 0 && max_ii < end_j) {
			int32_t tmp;
			tmp = comput_sc(&a[i], &a[max_ii], max_dist_x, max_dist_y, bw, chn_pen_gap, chn_pen_skip, 0, 1);
			if (tmp != INT32_MIN && max_f < tmp + f[max_ii])
				max_f = tmp + f[max_ii], max_j = max_ii;
		}
		f[i] = max_f, p[i] = max_j;
		v[i] = max_j >= 0 && v[max_j] > max_f? v[max_j] : max_f;
		if (max_ii < 0 || (a[i].x - a[max_ii].x <= (int64_t)max_dist_x && f[max_ii] < f[i]))
			max_ii = i;
	}
	kfree(km, x32); kfree(km, y32);
}
#endif

mm128_t *mg_lchain_dp(int max_dist_x, int max_dist_y, int bw, int max_skip, int max_iter, int min_cnt, int min_sc, float chn_pen_gap, float chn_pen_skip,
					  int is_cdna, int n_seg, int64_t n, mm128_t *a, int *n_u_, uint64_t **_u, void *km)
{ // TODO: make sure this works when n has more than 32 bits
	int32_t *f, *t, *v, n_u, n_v, mmax_f = 0, max_drop = bw;
	int64_t *p, i, j, max_ii, st = 0;
	uint64_t *u;

	if (_u) *_u = 0, *n_u_ = 0;
	if (n == 0 || a == 0) {
		kfree(km, a);
		return 0;
	}
	if (max_dist_x < bw) max_dist_x = bw;
	if (max_dist_y < bw && !is_cdna) max_dist_y = bw;
	if (is_cdna) max_drop = INT32_MAX;
	p = Kmalloc(km, int64_t, n);
	f = Kmalloc(km, int32_t, n);
	v = Kmalloc(km, int32_t, n);
	t = Kcalloc(km, int32_t, n);

#ifdef MM_AVX2_DISPATCH
	if (n_seg == 1 && !is_cdna && ksw_cpu_avx2()) { // the same f, p, v and t as the loop below
		mg_lchain_dp_fill_avx2(max_dist_x, max_dist_y, bw, max_skip, max_iter, chn_pen_gap, chn_pen_skip, n, a, f, p, v, t, km);
		goto fill_done;
	}
#endif
	// fill the score and backtrack arrays
	for (i = 0, max_ii = -1; i < n; ++i) {
		int64_t max_j = -1, end_j;
		int32_t max_f = a[i].y>>32&0xff, n_skip = 0;
		while (st < i && (a[i].x>>32 != a[st].x>>32 || a[i].x > a[st].x + max_dist_x)) ++st;
		if (i - st > max_iter) st = i - max_iter;
		for (j = i - 1; j >= st; --j) {
			int32_t sc;
			sc = comput_sc(&a[i], &a[j], max_dist_x, max_dist_y, bw, chn_pen_gap, chn_pen_skip, is_cdna, n_seg);
			if (sc == INT32_MIN) continue;
			sc += f[j];
			if (sc > max_f) {
				max_f = sc, max_j = j;
				if (n_skip > 0) --n_skip;
			} else if (t[j] == (int32_t)i) {
				if (++n_skip > max_skip)
					break;
			}
			if (p[j] >= 0) t[p[j]] = i;
		}
		end_j = j;
		if (max_ii < 0 || a[i].x - a[max_ii].x > (int64_t)max_dist_x) {
			int32_t max = INT32_MIN;
			max_ii = -1;
			for (j = i - 1; j >= st; --j)
				if (max < f[j]) max = f[j], max_ii = j;
		}
		if (max_ii >= 0 && max_ii < end_j) {
			int32_t tmp;
			tmp = comput_sc(&a[i], &a[max_ii], max_dist_x, max_dist_y, bw, chn_pen_gap, chn_pen_skip, is_cdna, n_seg);
			if (tmp != INT32_MIN && max_f < tmp + f[max_ii])
				max_f = tmp + f[max_ii], max_j = max_ii;
		}
		f[i] = max_f, p[i] = max_j;
		v[i] = max_j >= 0 && v[max_j] > max_f? v[max_j] : max_f; // v[] keeps the peak score up to i; f[] is the score ending at i, not always the peak
		if (max_ii < 0 || (a[i].x - a[max_ii].x <= (int64_t)max_dist_x && f[max_ii] < f[i]))
			max_ii = i;
		if (mmax_f < max_f) mmax_f = max_f;
		//fprintf(stderr, "X1\t%ld\t%ld:%d\t%ld\t%ld:%d\t%ld\t%ld\n", (long)i, (long)(a[i].x>>32), (int32_t)a[i].x, (long)max_j, max_j<0?-1L:(long)(a[max_j].x>>32), max_j<0?-1:(int32_t)a[max_j].x, (long)max_f, (long)v[i]);
	}

#ifdef MM_AVX2_DISPATCH
fill_done:
#endif
	u = mg_chain_backtrack(km, n, f, p, v, t, min_cnt, min_sc, max_drop, &n_u, &n_v);
	*n_u_ = n_u, *_u = u; // NB: note that u[] may not be sorted by score here
	kfree(km, p); kfree(km, f); kfree(km, t);
	if (n_u == 0) {
		kfree(km, a); kfree(km, v);
		return 0;
	}
	return compact_a(km, n_u, u, n_v, v, a);
}

typedef struct lc_elem_s {
	int32_t y;
	int64_t i;
	double pri;
	KRMQ_HEAD(struct lc_elem_s) head;
} lc_elem_t;

#define lc_elem_cmp(a, b) ((a)->y < (b)->y? -1 : (a)->y > (b)->y? 1 : ((a)->i > (b)->i) - ((a)->i < (b)->i))
#define lc_elem_lt2(a, b) ((a)->pri < (b)->pri)
KRMQ_INIT(lc_elem, lc_elem_t, head, lc_elem_cmp, lc_elem_lt2)

KALLOC_POOL_INIT(rmq, lc_elem_t)

static inline int32_t comput_sc_simple(const mm128_t *ai, const mm128_t *aj, float chn_pen_gap, float chn_pen_skip, int32_t *exact, int32_t *width)
{
	int32_t dq = (int32_t)ai->y - (int32_t)aj->y, dr, dd, dg, q_span, sc;
	dr = (int32_t)(ai->x - aj->x);
	*width = dd = dr > dq? dr - dq : dq - dr;
	dg = dr < dq? dr : dq;
	q_span = aj->y>>32&0xff;
	sc = q_span < dg? q_span : dg;
	if (exact) *exact = (dd == 0 && dg <= q_span);
	if (dd || dq > q_span) {
		float lin_pen, log_pen;
		lin_pen = chn_pen_gap * (float)dd + chn_pen_skip * (float)dg;
		log_pen = dd >= 1? mg_log2(dd + 1) : 0.0f; // mg_log2() only works for dd>=2
		sc -= (int)(lin_pen + .5f * log_pen);
	}
	return sc;
}

// The inner tree of mg_lchain_rmq() is only iterated in key order, so any ordered set of the same elements gives the
// same result. While the window [st_inner, i0) holds at most RMQ_INNER_HI anchors, there is no inner tree and a query
// scans the window; above that, a krmq tree is built from the window and kept until the window is below RMQ_INNER_LO.
#define RMQ_INNER_HI 256
#define RMQ_INNER_LO 64

static int64_t mm_rmq_inner_scan(const mm128_t *a, int64_t st, int64_t en, int32_t y_hi, int32_t y_lo, uint64_t *buf)
{ // the anchors of [st, en) with y_lo <= y <= y_hi in descending (y, i) order, as the inner tree iterates them
	int64_t j, n = 0, k;
	for (j = en - 1; j >= st; --j) {
		int32_t y = (int32_t)a[j].y;
		uint64_t key;
		if (y > y_hi || y < y_lo) continue;
		key = (uint64_t)(uint32_t)(y ^ INT32_MIN) << 32 | (uint64_t)j; // y as unsigned in the order of int32_t
		for (k = n; k > 0 && buf[k - 1] < key; --k) buf[k] = buf[k - 1]; // descending; nearly sorted in practice
		buf[k] = key, ++n;
	}
	return n;
}

mm128_t *mg_lchain_rmq(int max_dist, int max_dist_inner, int bw, int max_chn_skip, int cap_rmq_size, int min_cnt, int min_sc, float chn_pen_gap, float chn_pen_skip,
					   int64_t n, mm128_t *a, int *n_u_, uint64_t **_u, void *km)
{
	int32_t *f,*t, *v, n_u, n_v, mmax_f = 0, max_rmq_size = 0, max_drop = bw;
	int64_t *p, i, i0, st = 0, st_inner = 0;
	uint64_t *u;
	lc_elem_t *root = 0, *root_inner = 0;
	void *mem_mp = 0;
	kmp_rmq_t *mp;

	if (_u) *_u = 0, *n_u_ = 0;
	if (n == 0 || a == 0) {
		kfree(km, a);
		return 0;
	}
	if (max_dist < bw) max_dist = bw;
	if (max_dist_inner < 0) max_dist_inner = 0;
	if (max_dist_inner > max_dist) max_dist_inner = max_dist;
	p = Kmalloc(km, int64_t, n);
	f = Kmalloc(km, int32_t, n);
	t = Kcalloc(km, int32_t, n);
	v = Kmalloc(km, int32_t, n);
	mem_mp = km_init2(km, 0x10000);
	mp = kmp_init_rmq(mem_mp);
	int hybrid = max_dist_inner > 0 && cap_rmq_size > RMQ_INNER_HI; // the inner tree only for large windows
	int inner_tree = !hybrid;
	uint64_t *inner_buf = hybrid? Kmalloc(km, uint64_t, RMQ_INNER_HI + 1) : 0;

	// fill the score and backtrack arrays
	for (i = i0 = 0; i < n; ++i) {
		int64_t max_j = -1;
		int32_t q_span = a[i].y>>32&0xff, max_f = q_span;
		lc_elem_t s, *q, *r, lo, hi;
		// add in-range anchors
		if (i0 < i && a[i0].x != a[i].x) {
			int64_t j;
			for (j = i0; j < i; ++j) {
				q = kmp_alloc_rmq(mp);
				q->y = (int32_t)a[j].y, q->i = j, q->pri = -(f[j] + 0.5 * chn_pen_gap * ((int32_t)a[j].x + (int32_t)a[j].y));
				krmq_insert(lc_elem, &root, q, 0);
				if (max_dist_inner > 0 && inner_tree) {
					r = kmp_alloc_rmq(mp);
					*r = *q;
					krmq_insert(lc_elem, &root_inner, r, 0);
				}
			}
			i0 = i;
		}
		// get rid of active chains out of range
		while (st < i && (a[i].x>>32 != a[st].x>>32 || a[i].x > a[st].x + max_dist || krmq_size(head, root) > cap_rmq_size)) {
			s.y = (int32_t)a[st].y, s.i = st;
			if ((q = krmq_erase(lc_elem, &root, &s, 0)) != 0) // one walk: the same tree as find, then erase
				kmp_free_rmq(mp, q);
			++st;
		}
		if (max_dist_inner > 0 && !inner_tree) { // no inner tree: the window [st_inner, i0) holds all its elements
			while (st_inner < i && (a[i].x>>32 != a[st_inner].x>>32 || a[i].x > a[st_inner].x + max_dist_inner || (i0 > st_inner? i0 - st_inner : 0) > cap_rmq_size))
				++st_inner;
			if (i0 - st_inner > RMQ_INNER_HI) { // a large window: build the tree from it
				int64_t j;
				for (j = st_inner; j < i0; ++j) {
					r = kmp_alloc_rmq(mp);
					r->y = (int32_t)a[j].y, r->i = j, r->pri = -(f[j] + 0.5 * chn_pen_gap * ((int32_t)a[j].x + (int32_t)a[j].y));
					krmq_insert(lc_elem, &root_inner, r, 0);
				}
				inner_tree = 1;
			}
		} else if (max_dist_inner > 0)  { // similar to the block above, but applied to the inner tree
			while (st_inner < i && (a[i].x>>32 != a[st_inner].x>>32 || a[i].x > a[st_inner].x + max_dist_inner || krmq_size(head, root_inner) > cap_rmq_size)) {
				s.y = (int32_t)a[st_inner].y, s.i = st_inner;
				if ((q = krmq_erase(lc_elem, &root_inner, &s, 0)) != 0)
					kmp_free_rmq(mp, q);
				++st_inner;
			}
			if (hybrid && i0 - st_inner < RMQ_INNER_LO && (int64_t)krmq_size(head, root_inner) == (i0 > st_inner? i0 - st_inner : 0)) {
				while (root_inner) { // the window holds every element of the tree: drop the tree
					q = krmq_erase_first(lc_elem, &root_inner);
					kmp_free_rmq(mp, q);
				}
				inner_tree = 0;
			}
		}
		// RMQ
		lo.i = INT32_MAX, lo.y = (int32_t)a[i].y - max_dist;
		hi.i = 0, hi.y = (int32_t)a[i].y;
		if ((q = krmq_rmq(lc_elem, root, &lo, &hi)) != 0) {
			int32_t sc, exact, width, n_skip = 0;
			int64_t j = q->i;
			assert(q->y >= lo.y && q->y <= hi.y);
			sc = f[j] + comput_sc_simple(&a[i], &a[j], chn_pen_gap, chn_pen_skip, &exact, &width);
			if (width <= bw && sc > max_f) max_f = sc, max_j = j;
			if (!exact && !inner_tree && i0 > st_inner && (int32_t)a[i].y > 0) { // the inner query by a scan of the window
				int64_t k, nc = mm_rmq_inner_scan(a, st_inner, i0, (int32_t)a[i].y - 1, (int32_t)a[i].y - max_dist_inner, inner_buf);
				for (k = 0; k < nc; ++k) {
					int32_t width;
					j = (uint32_t)inner_buf[k];
					sc = f[j] + comput_sc_simple(&a[i], &a[j], chn_pen_gap, chn_pen_skip, 0, &width);
					if (width <= bw) {
						if (sc > max_f) {
							max_f = sc, max_j = j;
							if (n_skip > 0) --n_skip;
						} else if (t[j] == (int32_t)i) {
							if (++n_skip > max_chn_skip)
								break;
						}
						if (p[j] >= 0) t[p[j]] = i;
					}
				}
			} else if (!exact && root_inner && (int32_t)a[i].y > 0) {
				lc_elem_t *lo, *hi;
				s.y = (int32_t)a[i].y - 1, s.i = n;
				krmq_interval(lc_elem, root_inner, &s, &lo, &hi);
				if (lo) {
					const lc_elem_t *q;
					int32_t width;
					krmq_itr_t(lc_elem) itr;
					krmq_itr_find(lc_elem, root_inner, lo, &itr);
					while ((q = krmq_at(&itr)) != 0) {
						if (q->y < (int32_t)a[i].y - max_dist_inner) break;
						j = q->i;
						sc = f[j] + comput_sc_simple(&a[i], &a[j], chn_pen_gap, chn_pen_skip, 0, &width);
						if (width <= bw) {
							if (sc > max_f) {
								max_f = sc, max_j = j;
								if (n_skip > 0) --n_skip;
							} else if (t[j] == (int32_t)i) {
								if (++n_skip > max_chn_skip)
									break;
							}
							if (p[j] >= 0) t[p[j]] = i;
						}
						if (!krmq_itr_prev(lc_elem, &itr)) break;
					}
				}
			}
		}
		// set max
		assert(max_j < 0 || (a[max_j].x < a[i].x && (int32_t)a[max_j].y < (int32_t)a[i].y));
		f[i] = max_f, p[i] = max_j;
		v[i] = max_j >= 0 && v[max_j] > max_f? v[max_j] : max_f; // v[] keeps the peak score up to i; f[] is the score ending at i, not always the peak
		if (mmax_f < max_f) mmax_f = max_f;
		if (max_rmq_size < krmq_size(head, root)) max_rmq_size = krmq_size(head, root);
	}
	km_destroy(mem_mp);
	kfree(km, inner_buf);

	u = mg_chain_backtrack(km, n, f, p, v, t, min_cnt, min_sc, max_drop, &n_u, &n_v);
	*n_u_ = n_u, *_u = u; // NB: note that u[] may not be sorted by score here
	kfree(km, p); kfree(km, f); kfree(km, t);
	if (n_u == 0) {
		kfree(km, a); kfree(km, v);
		return 0;
	}
	return compact_a(km, n_u, u, n_v, v, a);
}
