#include <string.h>
#include <stdio.h>
#include <assert.h>
#include "ksw2.h"

#ifdef __SSE2__

#if defined(__AVX2__)
#include <immintrin.h>
#define SIMD_INT __m256i
#define SIMD_SHIFT 5
#define simd_func(func) _mm256_##func
#define simd_funcw(func) _mm256_##func##_si256

#elif defined(__SSE2__)
#include <emmintrin.h>
#define SIMD_INT __m128i
#define SIMD_SHIFT 4
#define simd_func(func) _mm_##func
#define simd_funcw(func) _mm_##func##_si128

#ifdef KSW_SSE2_ONLY
#undef __SSE4_1__
#endif

#ifdef __SSE4_1__
#include <smmintrin.h>
#endif
#endif // defined(__SSE2__)

#define SIMD_WIDTH (1<<SIMD_SHIFT)


#if defined(__AVX2__)
static inline __m256i simd_slli_1(__m256i x)
{
	return _mm256_insert_epi8(_mm256_slli_si256(x, 1), _mm256_extract_epi8(x, 15), 16);
}
static inline __m256i simd_srli_last(__m256i x)
{
	return _mm256_insert_epi8(_mm256_setzero_si256(), _mm256_extract_epi8(x, 31), 0);
}
#elif defined(__SSE2__)
static inline __m128i simd_slli_1(__m128i x) { return _mm_slli_si128(x, 1); }
static inline __m128i simd_srli_last(__m128i x) { return _mm_srli_si128(x, 15); }
#endif


#ifdef KSW_CPU_DISPATCH
#if defined(__AVX2__)
void ksw_extd2_avx2(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat,
				   int8_t q, int8_t e, int8_t q2, int8_t e2, int w, int zdrop, int end_bonus, int flag, ksw_extz_t *ez)
#elif defined(__SSE4_1__)
void ksw_extd2_sse41(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat,
				   int8_t q, int8_t e, int8_t q2, int8_t e2, int w, int zdrop, int end_bonus, int flag, ksw_extz_t *ez)
#elif defined(__SSE2__)
void ksw_extd2_sse2(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat,
				   int8_t q, int8_t e, int8_t q2, int8_t e2, int w, int zdrop, int end_bonus, int flag, ksw_extz_t *ez)
#endif
#else
void ksw_extd2_sse(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat,
				   int8_t q, int8_t e, int8_t q2, int8_t e2, int w, int zdrop, int end_bonus, int flag, ksw_extz_t *ez)
#endif // ~KSW_CPU_DISPATCH
{
#define __dp_code_block1 \
	z = simd_funcw(load)(&s[t]); \
	xt1 = simd_funcw(load)(&x[t]);                        /* xt1 <- x[r-1][t..t+15] */ \
	tmp = simd_srli_last(xt1);                            /* tmp <- x[r-1][t+15] */ \
	xt1 = simd_funcw(or)(simd_slli_1(xt1), x1_);          /* xt1 <- x[r-1][t-1..t+14] */ \
	x1_ = tmp; \
	vt1 = simd_funcw(load)(&v[t]);                        /* vt1 <- v[r-1][t..t+15] */ \
	tmp = simd_srli_last(vt1);                            /* tmp <- v[r-1][t+15] */ \
	vt1 = simd_funcw(or)(simd_slli_1(vt1), v1_);          /* vt1 <- v[r-1][t-1..t+14] */ \
	v1_ = tmp; \
	a = simd_func(add_epi8)(xt1, vt1);                    /* a <- x[r-1][t-1..t+14] + v[r-1][t-1..t+14] */ \
	ut = simd_funcw(load)(&u[t]);                         /* ut <- u[t..t+15] */ \
	b = simd_func(add_epi8)(simd_funcw(load)(&y[t]), ut); /* b <- y[r-1][t..t+15] + u[r-1][t..t+15] */ \
	x2t1= simd_funcw(load)(&x2[t]); \
	tmp = simd_srli_last(x2t1); \
	x2t1= simd_funcw(or)(simd_slli_1(x2t1), x21_); \
	x21_= tmp; \
	a2= simd_func(add_epi8)(x2t1, vt1); \
	b2= simd_func(add_epi8)(simd_funcw(load)(&y2[t]), ut);

#define __dp_code_block2 \
	simd_funcw(store)(&u[t], simd_func(sub_epi8)(z, vt1));/* u[r][t..t+15] <- z - v[r-1][t-1..t+14] */ \
	simd_funcw(store)(&v[t], simd_func(sub_epi8)(z, ut)); /* v[r][t..t+15] <- z - u[r-1][t..t+15] */ \
	tmp = simd_func(sub_epi8)(z, q_); \
	a = simd_func(sub_epi8)(a, tmp); \
	b = simd_func(sub_epi8)(b, tmp); \
	tmp = simd_func(sub_epi8)(z, q2_); \
	a2= simd_func(sub_epi8)(a2, tmp); \
	b2= simd_func(sub_epi8)(b2, tmp);

	int r, t, qe = q + e, n_col_, *off = 0, *off_end = 0, tlen_, qlen_, last_st, last_en, wl, wr, max_sc, min_sc, long_thres, long_diff;
	int with_cigar = !(flag&KSW_EZ_SCORE_ONLY), approx_max = !!(flag&KSW_EZ_APPROX_MAX);
	int32_t *H = 0, H0 = 0, last_H0_t = 0;
	uint8_t *qr, *sf, *mem, *mem2 = 0;
	SIMD_INT q_, q2_, qe_, qe2_, zero_, sc_mch_, sc_mis_, m1_, sc_N_, mask1_;
	SIMD_INT *u, *v, *x, *y, *x2, *y2, *s, *p = 0;

	ksw_reset_extz(ez);
	if (m <= 1 || qlen <= 0 || tlen <= 0) return;

	if (q2 + e2 < q + e) t = q, q = q2, q2 = t, t = e, e = e2, e2 = t; // make sure q+e no larger than q2+e2

	zero_   = simd_func(set1_epi8)(0);
	q_      = simd_func(set1_epi8)(q);
	q2_     = simd_func(set1_epi8)(q2);
	qe_     = simd_func(set1_epi8)(q + e);
	qe2_    = simd_func(set1_epi8)(q2 + e2);
	sc_mch_ = simd_func(set1_epi8)(mat[0]);
	sc_mis_ = simd_func(set1_epi8)(mat[1]);
	sc_N_   = mat[m*m-1] == 0? simd_func(set1_epi8)(-e2) : simd_func(set1_epi8)(mat[m*m-1]);
	m1_     = simd_func(set1_epi8)(m - 1); // wildcard

#if defined(__AVX2__)
	mask1_ = _mm256_setr_epi32(0xff, 0, 0, 0, 0, 0, 0, 0);
#elif defined(__SSE2__)
	mask1_ = _mm_setr_epi32(0xff, 0, 0, 0);
#endif

	if (w < 0) w = tlen > qlen? tlen : qlen;
	wl = wr = w;
	tlen_ = (tlen + SIMD_WIDTH - 1) / SIMD_WIDTH;
	n_col_ = qlen < tlen? qlen : tlen;
	n_col_ = ((n_col_ < w + 1? n_col_ : w + 1) + SIMD_WIDTH - 1) / SIMD_WIDTH + 1;
	qlen_ = (qlen + SIMD_WIDTH - 1) / SIMD_WIDTH;
	for (t = 1, max_sc = mat[0], min_sc = mat[1]; t < m * m; ++t) {
		max_sc = max_sc > mat[t]? max_sc : mat[t];
		min_sc = min_sc < mat[t]? min_sc : mat[t];
	}
	if (-min_sc > 2 * (q + e)) return; // otherwise, we won't see any mismatches

	long_thres = e != e2? (q2 - q) / (e - e2) - 1 : 0;
	if (q2 + e2 + long_thres * e2 > q + e + long_thres * e)
		++long_thres;
	long_diff = long_thres * (e - e2) - (q2 - q) - e2;

	mem = (uint8_t*)kcalloc(km, tlen_ * 8 + qlen_ + 1, SIMD_WIDTH);
	u = (SIMD_INT*)(((size_t)mem + SIMD_WIDTH - 1) >> SIMD_SHIFT << SIMD_SHIFT); // 16-byte aligned
	v = u + tlen_, x = v + tlen_, y = x + tlen_, x2 = y + tlen_, y2 = x2 + tlen_;
	s = y2 + tlen_, sf = (uint8_t*)(s + tlen_), qr = sf + tlen_ * SIMD_WIDTH;
	memset(u,  -q  - e,  tlen_ * SIMD_WIDTH);
	memset(v,  -q  - e,  tlen_ * SIMD_WIDTH);
	memset(x,  -q  - e,  tlen_ * SIMD_WIDTH);
	memset(y,  -q  - e,  tlen_ * SIMD_WIDTH);
	memset(x2, -q2 - e2, tlen_ * SIMD_WIDTH);
	memset(y2, -q2 - e2, tlen_ * SIMD_WIDTH);
	if (!approx_max) {
		H = (int32_t*)kmalloc(km, tlen_ * SIMD_WIDTH * 4);
		for (t = 0; t < tlen_ * SIMD_WIDTH; ++t) H[t] = KSW_NEG_INF;
	}
	if (with_cigar) {
		mem2 = (uint8_t*)kmalloc(km, ((size_t)(qlen + tlen - 1) * n_col_ + 1) * SIMD_WIDTH);
		p = (SIMD_INT*)(((size_t)mem2 + SIMD_WIDTH - 1) >> SIMD_SHIFT << SIMD_SHIFT);
		off = (int*)kmalloc(km, (qlen + tlen - 1) * sizeof(int) * 2);
		off_end = off + qlen + tlen - 1;
	}

	for (t = 0; t < qlen; ++t) qr[t] = query[qlen - 1 - t];
	memcpy(sf, target, tlen);

	for (r = 0, last_st = last_en = -1; r < qlen + tlen - 1; ++r) {
		int st = 0, en = tlen - 1, st0, en0, st_, en_;
		int8_t x1, x21, v1;
		uint8_t *qrr = qr + (qlen - 1 - r);
		int8_t *u8 = (int8_t*)u, *v8 = (int8_t*)v, *x8 = (int8_t*)x, *x28 = (int8_t*)x2;
		SIMD_INT x1_, x21_, v1_;
		// find the boundaries
		if (st < r - qlen + 1) st = r - qlen + 1;
		if (en > r) en = r;
		if (st < (r-wr+1)>>1) st = (r-wr+1)>>1; // take the ceil
		if (en > (r+wl)>>1) en = (r+wl)>>1; // take the floor
		if (st > en) {
			ez->zdropped = 1;
			break;
		}
		st0 = st, en0 = en;
		st = st / SIMD_WIDTH * SIMD_WIDTH, en = (en + SIMD_WIDTH) / SIMD_WIDTH * SIMD_WIDTH - 1;
		// set boundary conditions
		if (st > 0) {
			if (st - 1 >= last_st && st - 1 <= last_en) {
				x1 = x8[st - 1], x21 = x28[st - 1], v1 = v8[st - 1]; // (r-1,s-1) calculated in the last round
			} else {
				x1 = -q - e, x21 = -q2 - e2;
				v1 = -q - e;
			}
		} else {
			x1 = -q - e, x21 = -q2 - e2;
			v1 = r == 0? -q - e : r < long_thres? -e : r == long_thres? long_diff : -e2;
		}
		if (en >= r) {
			((int8_t*)y)[r] = -q - e, ((int8_t*)y2)[r] = -q2 - e2;
			u8[r] = r == 0? -q - e : r < long_thres? -e : r == long_thres? long_diff : -e2;
		}
		// loop fission: set scores first
		if (!(flag & KSW_EZ_GENERIC_SC)) {
			for (t = st0; t <= en0; t += SIMD_WIDTH) {
				SIMD_INT sq, st, tmp, mask;
				sq = simd_funcw(loadu)((SIMD_INT*)&sf[t]);
				st = simd_funcw(loadu)((SIMD_INT*)&qrr[t]);
				mask = simd_funcw(or)(simd_func(cmpeq_epi8)(sq, m1_), simd_func(cmpeq_epi8)(st, m1_));
				tmp = simd_func(cmpeq_epi8)(sq, st);
#if defined(__SSE4_1__) || defined(__AVX2__)
				tmp = simd_func(blendv_epi8)(sc_mis_, sc_mch_, tmp);
				tmp = simd_func(blendv_epi8)(tmp,     sc_N_,   mask);
#elif defined(__SSE2__) // emulate blendv
				tmp = _mm_or_si128(_mm_andnot_si128(tmp,  sc_mis_), _mm_and_si128(tmp,  sc_mch_));
				tmp = _mm_or_si128(_mm_andnot_si128(mask, tmp),     _mm_and_si128(mask, sc_N_));
#endif
				simd_funcw(storeu)((SIMD_INT*)((int8_t*)s + t), tmp);
			}
		} else {
			for (t = st0; t <= en0; ++t)
				((uint8_t*)s)[t] = mat[sf[t] * m + qrr[t]];
		}
		// core loop
		x1_  = simd_funcw(and)(simd_func(set1_epi8)((uint8_t)x1),  mask1_);
		x21_ = simd_funcw(and)(simd_func(set1_epi8)((uint8_t)x21), mask1_);
		v1_  = simd_funcw(and)(simd_func(set1_epi8)((uint8_t)v1),  mask1_);
		st_ = st / SIMD_WIDTH, en_ = en / SIMD_WIDTH;
		assert(en_ - st_ + 1 <= n_col_);
		if (!with_cigar) { // score only
			for (t = st_; t <= en_; ++t) {
				SIMD_INT z, a, b, a2, b2, xt1, x2t1, vt1, ut, tmp;
				__dp_code_block1;
#if defined(__SSE4_1__) || defined(__AVX2__)
				z = simd_func(max_epi8)(z, a);
				z = simd_func(max_epi8)(z, b);
				z = simd_func(max_epi8)(z, a2);
				z = simd_func(max_epi8)(z, b2);
				z = simd_func(min_epi8)(z, sc_mch_);
				__dp_code_block2; // save u[] and v[]; update a, b, a2 and b2
				simd_funcw(store)(&x[t],  simd_func(sub_epi8)(simd_func(max_epi8)(a,  zero_), qe_));
				simd_funcw(store)(&y[t],  simd_func(sub_epi8)(simd_func(max_epi8)(b,  zero_), qe_));
				simd_funcw(store)(&x2[t], simd_func(sub_epi8)(simd_func(max_epi8)(a2, zero_), qe2_));
				simd_funcw(store)(&y2[t], simd_func(sub_epi8)(simd_func(max_epi8)(b2, zero_), qe2_));
#elif defined(__SSE2__)
				tmp = _mm_cmpgt_epi8(a,  z);
				z = _mm_or_si128(_mm_andnot_si128(tmp, z), _mm_and_si128(tmp, a));
				tmp = _mm_cmpgt_epi8(b,  z);
				z = _mm_or_si128(_mm_andnot_si128(tmp, z), _mm_and_si128(tmp, b));
				tmp = _mm_cmpgt_epi8(a2, z);
				z = _mm_or_si128(_mm_andnot_si128(tmp, z), _mm_and_si128(tmp, a2));
				tmp = _mm_cmpgt_epi8(b2, z);
				z = _mm_or_si128(_mm_andnot_si128(tmp, z), _mm_and_si128(tmp, b2));
				tmp = _mm_cmplt_epi8(sc_mch_, z);
				z = _mm_or_si128(_mm_and_si128(tmp, sc_mch_), _mm_andnot_si128(tmp, z));
				__dp_code_block2;
				tmp = _mm_cmpgt_epi8(a, zero_);
				_mm_store_si128(&x[t],  _mm_sub_epi8(_mm_and_si128(tmp, a),  qe_));
				tmp = _mm_cmpgt_epi8(b, zero_);
				_mm_store_si128(&y[t],  _mm_sub_epi8(_mm_and_si128(tmp, b),  qe_));
				tmp = _mm_cmpgt_epi8(a2, zero_);
				_mm_store_si128(&x2[t], _mm_sub_epi8(_mm_and_si128(tmp, a2), qe2_));
				tmp = _mm_cmpgt_epi8(b2, zero_);
				_mm_store_si128(&y2[t], _mm_sub_epi8(_mm_and_si128(tmp, b2), qe2_));
#endif
			}
		} else if (!(flag&KSW_EZ_RIGHT)) { // gap left-alignment
			SIMD_INT *pr = p + (size_t)r * n_col_ - st_;
			off[r] = st, off_end[r] = en;
			for (t = st_; t <= en_; ++t) {
				SIMD_INT d, z, a, b, a2, b2, xt1, x2t1, vt1, ut, tmp;
				__dp_code_block1;
#if defined(__SSE4_1__) || defined(__AVX2__)
				d = simd_funcw(and)(simd_func(cmpgt_epi8)(a, z), simd_func(set1_epi8)(1));       // d = a  > z? 1 : 0
				z = simd_func(max_epi8)(z, a);
				d = simd_func(blendv_epi8)(d, simd_func(set1_epi8)(2), simd_func(cmpgt_epi8)(b,  z)); // d = b  > z? 2 : d
				z = simd_func(max_epi8)(z, b);
				d = simd_func(blendv_epi8)(d, simd_func(set1_epi8)(3), simd_func(cmpgt_epi8)(a2, z)); // d = a2 > z? 3 : d
				z = simd_func(max_epi8)(z, a2);
				d = simd_func(blendv_epi8)(d, simd_func(set1_epi8)(4), simd_func(cmpgt_epi8)(b2, z)); // d = a2 > z? 3 : d
				z = simd_func(max_epi8)(z, b2);
				z = simd_func(min_epi8)(z, sc_mch_);
#elif defined(__SSE2__) // emulate SSE4.1 intrinsics _mm_max_epi8() and _mm_blendv_epi8()
				tmp = _mm_cmpgt_epi8(a,  z);
				d = _mm_and_si128(tmp, _mm_set1_epi8(1));
				z = _mm_or_si128(_mm_andnot_si128(tmp, z), _mm_and_si128(tmp, a));
				tmp = _mm_cmpgt_epi8(b,  z);
				d = _mm_or_si128(_mm_andnot_si128(tmp, d), _mm_and_si128(tmp, _mm_set1_epi8(2)));
				z = _mm_or_si128(_mm_andnot_si128(tmp, z), _mm_and_si128(tmp, b));
				tmp = _mm_cmpgt_epi8(a2, z);
				d = _mm_or_si128(_mm_andnot_si128(tmp, d), _mm_and_si128(tmp, _mm_set1_epi8(3)));
				z = _mm_or_si128(_mm_andnot_si128(tmp, z), _mm_and_si128(tmp, a2));
				tmp = _mm_cmpgt_epi8(b2, z);
				d = _mm_or_si128(_mm_andnot_si128(tmp, d), _mm_and_si128(tmp, _mm_set1_epi8(4)));
				z = _mm_or_si128(_mm_andnot_si128(tmp, z), _mm_and_si128(tmp, b2));
				tmp = _mm_cmplt_epi8(sc_mch_, z);
				z = _mm_or_si128(_mm_and_si128(tmp, sc_mch_), _mm_andnot_si128(tmp, z));
#endif
				__dp_code_block2;
				tmp = simd_func(cmpgt_epi8)(a, zero_);
				simd_funcw(store)(&x[t],  simd_func(sub_epi8)(simd_funcw(and)(tmp, a),  qe_));
				d = simd_funcw(or)(d, simd_funcw(and)(tmp, simd_func(set1_epi8)(0x08))); // d = a > 0? 1<<3 : 0
				tmp = simd_func(cmpgt_epi8)(b, zero_);
				simd_funcw(store)(&y[t],  simd_func(sub_epi8)(simd_funcw(and)(tmp, b),  qe_));
				d = simd_funcw(or)(d, simd_funcw(and)(tmp, simd_func(set1_epi8)(0x10))); // d = b > 0? 1<<4 : 0
				tmp = simd_func(cmpgt_epi8)(a2, zero_);
				simd_funcw(store)(&x2[t], simd_func(sub_epi8)(simd_funcw(and)(tmp, a2), qe2_));
				d = simd_funcw(or)(d, simd_funcw(and)(tmp, simd_func(set1_epi8)(0x20))); // d = a > 0? 1<<5 : 0
				tmp = simd_func(cmpgt_epi8)(b2, zero_);
				simd_funcw(store)(&y2[t], simd_func(sub_epi8)(simd_funcw(and)(tmp, b2), qe2_));
				d = simd_funcw(or)(d, simd_funcw(and)(tmp, simd_func(set1_epi8)(0x40))); // d = b > 0? 1<<6 : 0
				simd_funcw(store)(&pr[t], d);
			}
		} else { // gap right-alignment
			SIMD_INT *pr = p + (size_t)r * n_col_ - st_;
			off[r] = st, off_end[r] = en;
			for (t = st_; t <= en_; ++t) {
				SIMD_INT d, z, a, b, a2, b2, xt1, x2t1, vt1, ut, tmp;
				__dp_code_block1;
#if defined(__SSE4_1__) || defined(__AVX2__)
				d = simd_funcw(andnot)(simd_func(cmpgt_epi8)(z, a), simd_func(set1_epi8)(1));    // d = z > a?  0 : 1
				z = simd_func(max_epi8)(z, a);
				d = simd_func(blendv_epi8)(simd_func(set1_epi8)(2), d, simd_func(cmpgt_epi8)(z, b));  // d = z > b?  d : 2
				z = simd_func(max_epi8)(z, b);
				d = simd_func(blendv_epi8)(simd_func(set1_epi8)(3), d, simd_func(cmpgt_epi8)(z, a2)); // d = z > a2? d : 3
				z = simd_func(max_epi8)(z, a2);
				d = simd_func(blendv_epi8)(simd_func(set1_epi8)(4), d, simd_func(cmpgt_epi8)(z, b2)); // d = z > b2? d : 4
				z = simd_func(max_epi8)(z, b2);
				z = simd_func(min_epi8)(z, sc_mch_);
#elif defined(__SSE2__)
				tmp = _mm_cmpgt_epi8(z, a);
				d = _mm_andnot_si128(tmp, _mm_set1_epi8(1));
				z = _mm_or_si128(_mm_and_si128(tmp, z), _mm_andnot_si128(tmp, a));
				tmp = _mm_cmpgt_epi8(z, b);
				d = _mm_or_si128(_mm_and_si128(tmp, d), _mm_andnot_si128(tmp, _mm_set1_epi8(2)));
				z = _mm_or_si128(_mm_and_si128(tmp, z), _mm_andnot_si128(tmp, b));
				tmp = _mm_cmpgt_epi8(z, a2);
				d = _mm_or_si128(_mm_and_si128(tmp, d), _mm_andnot_si128(tmp, _mm_set1_epi8(3)));
				z = _mm_or_si128(_mm_and_si128(tmp, z), _mm_andnot_si128(tmp, a2));
				tmp = _mm_cmpgt_epi8(z, b2);
				d = _mm_or_si128(_mm_and_si128(tmp, d), _mm_andnot_si128(tmp, _mm_set1_epi8(4)));
				z = _mm_or_si128(_mm_and_si128(tmp, z), _mm_andnot_si128(tmp, b2));
				tmp = _mm_cmplt_epi8(sc_mch_, z);
				z = _mm_or_si128(_mm_and_si128(tmp, sc_mch_), _mm_andnot_si128(tmp, z));
#endif
				__dp_code_block2;
				tmp = simd_func(cmpgt_epi8)(zero_, a);
				simd_funcw(store)(&x[t],  simd_func(sub_epi8)(simd_funcw(andnot)(tmp, a),  qe_));
				d = simd_funcw(or)(d, simd_funcw(andnot)(tmp, simd_func(set1_epi8)(0x08))); // d = a > 0? 1<<3 : 0
				tmp = simd_func(cmpgt_epi8)(zero_, b);
				simd_funcw(store)(&y[t],  simd_func(sub_epi8)(simd_funcw(andnot)(tmp, b),  qe_));
				d = simd_funcw(or)(d, simd_funcw(andnot)(tmp, simd_func(set1_epi8)(0x10))); // d = b > 0? 1<<4 : 0
				tmp = simd_func(cmpgt_epi8)(zero_, a2);
				simd_funcw(store)(&x2[t], simd_func(sub_epi8)(simd_funcw(andnot)(tmp, a2), qe2_));
				d = simd_funcw(or)(d, simd_funcw(andnot)(tmp, simd_func(set1_epi8)(0x20))); // d = a > 0? 1<<5 : 0
				tmp = simd_func(cmpgt_epi8)(zero_, b2);
				simd_funcw(store)(&y2[t], simd_func(sub_epi8)(simd_funcw(andnot)(tmp, b2), qe2_));
				d = simd_funcw(or)(d, simd_funcw(andnot)(tmp, simd_func(set1_epi8)(0x40))); // d = b > 0? 1<<6 : 0
				simd_funcw(store)(&pr[t], d);
			}
		}
		if (!approx_max) { // find the exact max with a 32-bit score array
			int32_t max_H, max_t;
			// compute H[], max_H and max_t
			if (r > 0) {
				int32_t HH[SIMD_WIDTH/4], tt[SIMD_WIDTH/4], en1 = st0 + (en0 - st0) / (SIMD_WIDTH/4) * (SIMD_WIDTH/4), i;
				SIMD_INT max_H_, max_t_;
				max_H = H[en0] = en0 > 0? H[en0-1] + u8[en0] : H[en0] + v8[en0]; // special casing the last element
				max_t = en0;
				max_H_ = simd_func(set1_epi32)(max_H);
				max_t_ = simd_func(set1_epi32)(max_t);
				for (t = st0; t < en1; t += SIMD_WIDTH/4) { // this implements: H[t]+=v8[t]-qe; if(H[t]>max_H) max_H=H[t],max_t=t;
					SIMD_INT H1, tmp, t_;
					H1 = simd_funcw(loadu)((SIMD_INT*)&H[t]);
#if defined(__AVX2__)
					t_ = _mm256_setr_epi32(v8[t], v8[t+1], v8[t+2], v8[t+3], v8[t+4], v8[t+5], v8[t+6], v8[t+7]);
#elif defined(__SSE2__)
					t_ = _mm_setr_epi32(v8[t], v8[t+1], v8[t+2], v8[t+3]);
#endif
					H1 = simd_func(add_epi32)(H1, t_);
					simd_funcw(storeu)((SIMD_INT*)&H[t], H1);
					t_ = simd_func(set1_epi32)(t);
					tmp = simd_func(cmpgt_epi32)(H1, max_H_);
#if defined(__SSE4_1__) || defined(__AVX2__)
					max_H_ = simd_func(blendv_epi8)(max_H_, H1, tmp);
					max_t_ = simd_func(blendv_epi8)(max_t_, t_, tmp);
#elif defined(__SSE2__)
					max_H_ = simd_funcw(or)(simd_funcw(and)(tmp, H1), simd_funcw(andnot)(tmp, max_H_));
					max_t_ = simd_funcw(or)(simd_funcw(and)(tmp, t_), simd_funcw(andnot)(tmp, max_t_));
#endif
				}
				simd_funcw(storeu)((SIMD_INT*)HH, max_H_);
				simd_funcw(storeu)((SIMD_INT*)tt, max_t_);
				for (i = 0; i < SIMD_WIDTH/4; ++i)
					if (max_H < HH[i]) max_H = HH[i], max_t = tt[i] + i;
				for (; t < en0; ++t) { // for the rest of values that haven't been computed with SSE
					H[t] += (int32_t)v8[t];
					if (H[t] > max_H)
						max_H = H[t], max_t = t;
				}
			} else H[0] = v8[0] - qe, max_H = H[0], max_t = 0; // special casing r==0
			// update ez
			if (en0 == tlen - 1 && H[en0] > ez->mte)
				ez->mte = H[en0], ez->mte_q = r - en;
			if (r - st0 == qlen - 1 && H[st0] > ez->mqe)
				ez->mqe = H[st0], ez->mqe_t = st0;
			if (ksw_apply_zdrop(ez, 1, max_H, r, max_t, zdrop, e2)) break;
			if (r == qlen + tlen - 2 && en0 == tlen - 1)
				ez->score = H[tlen - 1];
		} else { // find approximate max; Z-drop might be inaccurate, too.
			if (r > 0) {
				if (last_H0_t >= st0 && last_H0_t <= en0 && last_H0_t + 1 >= st0 && last_H0_t + 1 <= en0) {
					int32_t d0 = v8[last_H0_t];
					int32_t d1 = u8[last_H0_t + 1];
					if (d0 > d1) H0 += d0;
					else H0 += d1, ++last_H0_t;
				} else if (last_H0_t >= st0 && last_H0_t <= en0) {
					H0 += v8[last_H0_t];
				} else {
					++last_H0_t, H0 += u8[last_H0_t];
				}
			} else H0 = v8[0] - qe, last_H0_t = 0;
			if ((flag & KSW_EZ_APPROX_DROP) && ksw_apply_zdrop(ez, 1, H0, r, last_H0_t, zdrop, e2)) break;
			if (r == qlen + tlen - 2 && en0 == tlen - 1)
				ez->score = H0;
		}
		last_st = st, last_en = en;
		//for (t = st0; t <= en0; ++t) printf("(%d,%d)\t(%d,%d,%d,%d)\t%d\n", r, t, ((int8_t*)u)[t], ((int8_t*)v)[t], ((int8_t*)x)[t], ((int8_t*)y)[t], H[t]); // for debugging
	}
	kfree(km, mem);
	if (!approx_max) kfree(km, H);
	if (with_cigar) { // backtrack
		int rev_cigar = !!(flag & KSW_EZ_REV_CIGAR);
		if (!ez->zdropped && !(flag&KSW_EZ_EXTZ_ONLY)) {
			ksw_backtrack(km, 1, rev_cigar, 0, (uint8_t*)p, off, off_end, n_col_*SIMD_WIDTH, tlen-1, qlen-1, &ez->m_cigar, &ez->n_cigar, &ez->cigar);
		} else if (!ez->zdropped && (flag&KSW_EZ_EXTZ_ONLY) && ez->mqe + end_bonus > (int)ez->max) {
			ez->reach_end = 1;
			ksw_backtrack(km, 1, rev_cigar, 0, (uint8_t*)p, off, off_end, n_col_*SIMD_WIDTH, ez->mqe_t, qlen-1, &ez->m_cigar, &ez->n_cigar, &ez->cigar);
		} else if (ez->max_t >= 0 && ez->max_q >= 0) {
			ksw_backtrack(km, 1, rev_cigar, 0, (uint8_t*)p, off, off_end, n_col_*SIMD_WIDTH, ez->max_t, ez->max_q, &ez->m_cigar, &ez->n_cigar, &ez->cigar);
		}
		kfree(km, mem2); kfree(km, off);
	}
}
#endif // __SSE2__
