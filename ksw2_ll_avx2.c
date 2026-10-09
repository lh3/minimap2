// ksw2_ll_avx2.c -- the 16-bit striped local alignment of ksw2_ll_sse.c (Farrar) with 16 lanes, score only.
//
// mm_test_zdrop() uses only the score of ksw_ll_i16(). The cell values are the same as in the SSE version (same
// saturating operations; padding cells after the query end score -1, as in ksw_ll_qinit() of 2.31, and never feed
// a cell of the query), so the score is the same. Without the end position, the copy of the H row at each new
// maximum is not necessary.
#include <stdint.h>
#include <string.h>
#include "kalloc.h"
#include "ksw2.h"

#ifdef __AVX2__
#include <immintrin.h>

// Shift the 16 lanes up by one 16-bit element (lane 0 gets 0), as _mm_slli_si128(x, 2) does for 8 lanes.
static inline __m256i mm_ll_shift1(__m256i x)
{
	return _mm256_alignr_epi8(x, _mm256_permute2x128_si256(x, x, 0x08), 14);
}

int ksw_ll_i16_score_avx2(void *km, int qlen, const uint8_t *query, int m, const int8_t *mat, int tlen, const uint8_t *target, int gapo, int gape)
{
	const int slen = (qlen + 15) / 16;
	int a, i, j, k, gmax = 0;
	uint8_t *mem = (uint8_t*)kmalloc(km, 32 * (size_t)slen * (m + 3) + 32);
	__m256i *qp = (__m256i*)(((size_t)mem + 31) >> 5 << 5), *H0 = qp + slen * m, *H1 = H0 + slen, *E = H1 + slen, *S;
	const __m256i zero = _mm256_setzero_si256(), gapoe = _mm256_set1_epi16(gapo + gape), vgape = _mm256_set1_epi16(gape);
	int16_t *t = (int16_t*)qp;

	for (a = 0; a < m; ++a) { // query profile: lane l of vector i is query position l * slen + i
		const int8_t *ma = mat + a * m;
		for (i = 0; i < slen; ++i)
			for (k = i; k < slen * 16; k += slen)
				*t++ = k >= qlen? -1 : ma[query[k]];
	}
	memset(E, 0, slen * 32);
	memset(H0, 0, slen * 32);
	for (i = 0; i < tlen; ++i) {
		__m256i e, h, f = zero, max = zero;
		S = qp + target[i] * slen;
		h = mm_ll_shift1(_mm256_load_si256(H0 + slen - 1));
		for (j = 0; j < slen; ++j) {
			h = _mm256_adds_epi16(h, _mm256_load_si256(S + j));
			e = _mm256_load_si256(E + j);
			h = _mm256_max_epi16(h, e);
			h = _mm256_max_epi16(h, f);
			max = _mm256_max_epi16(max, h);
			_mm256_store_si256(H1 + j, h);
			h = _mm256_subs_epu16(h, gapoe);
			e = _mm256_subs_epu16(e, vgape);
			e = _mm256_max_epi16(e, h);
			_mm256_store_si256(E + j, e);
			f = _mm256_subs_epu16(f, vgape);
			f = _mm256_max_epi16(f, h);
			h = _mm256_load_si256(H0 + j);
		}
		for (k = 0; k < 16; ++k) { // lazy F loop
			f = mm_ll_shift1(f);
			for (j = 0; j < slen; ++j) {
				h = _mm256_load_si256(H1 + j);
				h = _mm256_max_epi16(h, f);
				_mm256_store_si256(H1 + j, h);
				h = _mm256_subs_epu16(h, gapoe);
				f = _mm256_subs_epu16(f, vgape);
				if (!_mm256_movemask_epi8(_mm256_cmpgt_epi16(f, h))) goto end_loop;
			}
		}
end_loop:
		{
			__m128i x = _mm_max_epi16(_mm256_castsi256_si128(max), _mm256_extracti128_si256(max, 1));
			x = _mm_max_epi16(x, _mm_srli_si128(x, 8));
			x = _mm_max_epi16(x, _mm_srli_si128(x, 4));
			x = _mm_max_epi16(x, _mm_srli_si128(x, 2));
			const int imax = (uint16_t)_mm_extract_epi16(x, 0);
			if (imax > gmax) gmax = imax;
		}
		S = H1, H1 = H0, H0 = S;
	}
	kfree(km, mem);
	return gmax;
}
#endif
