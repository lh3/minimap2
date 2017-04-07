/* The MIT License

   Copyright (c) 2008, 2011 Attractive Chaos <attractor@live.co.uk>

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
*/

// This is a simplified version of ksort.h

#ifndef AC_KSORT_H
#define AC_KSORT_H

#include <stdlib.h>
#include <string.h>

typedef struct {
	void *left, *right;
	int depth;
} ks_isort_stack_t;

#define KSORT_SWAP(type_t, a, b) { register type_t t=(a); (a)=(b); (b)=t; }

#define KSORT_INIT(name, type_t, __sort_lt)								\
	size_t ks_lis_##name(size_t n, const type_t *a, size_t *b, size_t *_p) \
	{ /* translated from: http://www.algorithmist.com/index.php/Longest_Increasing_Subsequence.cpp */ \
		size_t i, u, v, *top = b, *p; \
		if (n == 0) return 0; \
		p = _p? _p : (size_t*)calloc(n, sizeof(size_t)); \
		*top++ = 0; \
		for (i = 1; i < n; i++) { \
			if (__sort_lt(a[*(top-1)], a[i])) { \
				p[i] = *(top-1); \
				*top++ = i; \
				continue; \
			} \
			for (u = 0, v = top - b - 1; u < v;) { \
				size_t c = (u + v) >> 1; \
				if (__sort_lt(a[b[c]], a[i])) u = c + 1; \
				else v = c; \
			} \
			if (__sort_lt(a[i], a[b[u]])) { \
				if (u > 0) p[i] = b[u-1]; \
				b[u] = i; \
			} \
		} \
		for (u = top - b, v = *(top-1); u--; v = p[v]) b[u] = v; \
		if (!_p) free(p); \
		return top - b; \
	} \
	type_t ks_ksmall_##name(size_t n, type_t arr[], size_t kk)			\
	{																	\
		type_t *low, *high, *k, *ll, *hh, *mid;							\
		low = arr; high = arr + n - 1; k = arr + kk;					\
		for (;;) {														\
			if (high <= low) return *k;									\
			if (high == low + 1) {										\
				if (__sort_lt(*high, *low)) KSORT_SWAP(type_t, *low, *high); \
				return *k;												\
			}															\
			mid = low + (high - low) / 2;								\
			if (__sort_lt(*high, *mid)) KSORT_SWAP(type_t, *mid, *high); \
			if (__sort_lt(*high, *low)) KSORT_SWAP(type_t, *low, *high); \
			if (__sort_lt(*low, *mid)) KSORT_SWAP(type_t, *mid, *low);	\
			KSORT_SWAP(type_t, *mid, *(low+1));							\
			ll = low + 1; hh = high;									\
			for (;;) {													\
				do ++ll; while (__sort_lt(*ll, *low));					\
				do --hh; while (__sort_lt(*low, *hh));					\
				if (hh < ll) break;										\
				KSORT_SWAP(type_t, *ll, *hh);							\
			}															\
			KSORT_SWAP(type_t, *low, *hh);								\
			if (hh <= k) low = ll;										\
			if (hh >= k) high = hh - 1;									\
		}																\
	}																	\

#define ks_ksmall(name, n, a, k) ks_ksmall_##name(n, a, k)

#define ks_lt_generic(a, b) ((a) < (b))
#define ks_lt_str(a, b) (strcmp((a), (b)) < 0)

typedef const char *ksstr_t;

#define KSORT_INIT_GENERIC(type_t) KSORT_INIT(type_t, type_t, ks_lt_generic)
#define KSORT_INIT_STR KSORT_INIT(str, ksstr_t, ks_lt_str)

#define RS_MIN_SIZE 64

#define KRADIX_SORT_INIT(name, rstype_t, rskey, sizeof_key) \
	typedef struct { \
		rstype_t *b, *e; \
	} rsbucket_##name##_t; \
	void rs_insertsort_##name(rstype_t *beg, rstype_t *end) \
	{ \
		rstype_t *i; \
		for (i = beg + 1; i < end; ++i) \
			if (rskey(*i) < rskey(*(i - 1))) { \
				rstype_t *j, tmp = *i; \
				for (j = i; j > beg && rskey(tmp) < rskey(*(j-1)); --j) \
					*j = *(j - 1); \
				*j = tmp; \
			} \
	} \
	void rs_sort_##name(rstype_t *beg, rstype_t *end, int n_bits, int s) \
	{ \
		rstype_t *i; \
		int size = 1<<n_bits, m = size - 1; \
		rsbucket_##name##_t *k, b[size], *be = b + size; \
		for (k = b; k != be; ++k) k->b = k->e = beg; \
		for (i = beg; i != end; ++i) ++b[rskey(*i)>>s&m].e; \
		for (k = b + 1; k != be; ++k) \
			k->e += (k-1)->e - beg, k->b = (k-1)->e; \
		for (k = b; k != be;) { \
			if (k->b != k->e) { \
				rsbucket_##name##_t *l; \
				if ((l = b + (rskey(*k->b)>>s&m)) != k) { \
					rstype_t tmp = *k->b, swap; \
					do { \
						swap = tmp; tmp = *l->b; *l->b++ = swap; \
						l = b + (rskey(tmp)>>s&m); \
					} while (l != k); \
					*k->b++ = tmp; \
				} else ++k->b; \
			} else ++k; \
		} \
		for (b->b = beg, k = b + 1; k != be; ++k) k->b = (k-1)->e; \
		if (s) { \
			s = s > n_bits? s - n_bits : 0; \
			for (k = b; k != be; ++k) \
				if (k->e - k->b > RS_MIN_SIZE) rs_sort_##name(k->b, k->e, n_bits, s); \
				else if (k->e - k->b > 1) rs_insertsort_##name(k->b, k->e); \
		} \
	} \
	void radix_sort_##name(rstype_t *beg, rstype_t *end) \
	{ \
		if (end - beg <= RS_MIN_SIZE) rs_insertsort_##name(beg, end); \
		else rs_sort_##name(beg, end, 8, sizeof_key * 8 - 8); \
	}

#endif
