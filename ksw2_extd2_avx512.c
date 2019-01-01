#include <string.h>
#include <stdio.h>
#include <assert.h>
#include "ksw2.h"

#ifdef __AVX512BW__
#include <immintrin.h>

void ksw_extd2_avx512(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat,
				   int8_t q, int8_t e, int8_t q2, int8_t e2, int w, int zdrop, int end_bonus, int flag, ksw_extz_t *ez)
{
#define __dp_code_block1 \
	z = _mm512_load_si512(&s[t]); \
	tmp = _mm512_loadu_si512((uint8_t*)&x[t] - 1); \
	xt1 = _mm512_mask_blend_epi8(1, tmp, x1_); \
	x1_ = _mm512_maskz_set1_epi8(1, *((uint8_t*)&x[t] + 63)); \
	tmp = _mm512_loadu_si512((uint8_t*)&v[t] - 1); \
	vt1 = _mm512_mask_blend_epi8(1, tmp, v1_); \
	v1_ = _mm512_maskz_set1_epi8(1, *((uint8_t*)&v[t] + 63)); \
	a = _mm512_add_epi8(xt1, vt1);                     /* a <- x[r-1][t-1..t+62] + v[r-1][t-1..t+62] */ \
	ut = _mm512_load_si512(&u[t]);                     /* ut <- u[t..t+63] */ \
	b = _mm512_add_epi8(_mm512_load_si512(&y[t]), ut); /* b <- y[r-1][t..t+63] + u[r-1][t..t+63] */ \
	tmp = _mm512_loadu_si512((uint8_t*)&x2[t] - 1); \
	x2t1 = _mm512_mask_blend_epi8(1, tmp, x21_); \
	x21_ = _mm512_maskz_set1_epi8(1, *((uint8_t*)&x2[t] + 63)); \
	a2= _mm512_add_epi8(x2t1, vt1); \
	b2= _mm512_add_epi8(_mm512_load_si512(&y2[t]), ut);

#define __dp_code_block2 \
	_mm512_store_si512(&u[t], _mm512_sub_epi8(z, vt1));/* u[r][t..t+15] <- z - v[r-1][t-1..t+14] */ \
	_mm512_store_si512(&v[t], _mm512_sub_epi8(z, ut)); /* v[r][t..t+15] <- z - u[r-1][t..t+15] */ \
	tmp = _mm512_sub_epi8(z, q_); \
	a = _mm512_sub_epi8(a, tmp); \
	b = _mm512_sub_epi8(b, tmp); \
	tmp = _mm512_sub_epi8(z, q2_); \
	a2= _mm512_sub_epi8(a2, tmp); \
	b2= _mm512_sub_epi8(b2, tmp);

	int r, t, qe = q + e, n_col_, *off = 0, *off_end = 0, tlen_, qlen_, last_st, last_en, wl, wr, max_sc, min_sc, long_thres, long_diff;
	int with_cigar = !(flag&KSW_EZ_SCORE_ONLY), approx_max = !!(flag&KSW_EZ_APPROX_MAX);
	int32_t *H = 0, H0 = 0, last_H0_t = 0;
	uint8_t *qr, *sf, *mem, *mem2 = 0;
	__m512i q_, q2_, qe_, qe2_, zero_, sc_mch_, sc_mis_, m1_, sc_N_;
	__m512i *u, *v, *x, *y, *x2, *y2, *s, *p = 0;

	ksw_reset_extz(ez);
	if (m <= 1 || qlen <= 0 || tlen <= 0) return;

	if (q2 + e2 < q + e) t = q, q = q2, q2 = t, t = e, e = e2, e2 = t; // make sure q+e no larger than q2+e2

	zero_   = _mm512_set1_epi8(0);
	q_      = _mm512_set1_epi8(q);
	q2_     = _mm512_set1_epi8(q2);
	qe_     = _mm512_set1_epi8(q + e);
	qe2_    = _mm512_set1_epi8(q2 + e2);
	sc_mch_ = _mm512_set1_epi8(mat[0]);
	sc_mis_ = _mm512_set1_epi8(mat[1]);
	sc_N_   = mat[m*m-1] == 0? _mm512_set1_epi8(-e2) : _mm512_set1_epi8(mat[m*m-1]);
	m1_     = _mm512_set1_epi8(m - 1); // wildcard

	if (w < 0) w = tlen > qlen? tlen : qlen;
	wl = wr = w;
	tlen_ = (tlen + 64 - 1) / 64;
	n_col_ = qlen < tlen? qlen : tlen;
	n_col_ = ((n_col_ < w + 1? n_col_ : w + 1) + 64 - 1) / 64 + 1;
	qlen_ = (qlen + 64 - 1) / 64;
	for (t = 1, max_sc = mat[0], min_sc = mat[1]; t < m * m; ++t) {
		max_sc = max_sc > mat[t]? max_sc : mat[t];
		min_sc = min_sc < mat[t]? min_sc : mat[t];
	}
	if (-min_sc > 2 * (q + e)) return; // otherwise, we won't see any mismatches

	long_thres = e != e2? (q2 - q) / (e - e2) - 1 : 0;
	if (q2 + e2 + long_thres * e2 > q + e + long_thres * e)
		++long_thres;
	long_diff = long_thres * (e - e2) - (q2 - q) - e2;

	mem = (uint8_t*)kcalloc(km, tlen_ * 8 + qlen_ + 1, 64);
	u = (__m512i*)(((size_t)mem + 64 - 1) >> 6 << 6); // 16-byte aligned
	v = u + tlen_, x = v + tlen_, y = x + tlen_, x2 = y + tlen_, y2 = x2 + tlen_;
	s = y2 + tlen_, sf = (uint8_t*)(s + tlen_), qr = sf + tlen_ * 64;
	memset(u,  -q  - e,  tlen_ * 64);
	memset(v,  -q  - e,  tlen_ * 64);
	memset(x,  -q  - e,  tlen_ * 64);
	memset(y,  -q  - e,  tlen_ * 64);
	memset(x2, -q2 - e2, tlen_ * 64);
	memset(y2, -q2 - e2, tlen_ * 64);
	if (!approx_max) {
		H = (int32_t*)kmalloc(km, tlen_ * 64 * 4);
		for (t = 0; t < tlen_ * 64; ++t) H[t] = KSW_NEG_INF;
	}
	if (with_cigar) {
		mem2 = (uint8_t*)kmalloc(km, ((size_t)(qlen + tlen - 1) * n_col_ + 1) * 64);
		p = (__m512i*)(((size_t)mem2 + 64 - 1) >> 6 << 6);
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
		__m512i x1_, x21_, v1_;
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
		st = st / 64 * 64, en = (en + 64) / 64 * 64 - 1;
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
			for (t = st0; t <= en0; t += 64) {
				__m512i sq, st, tmp;
				__mmask64 mask;
				sq = _mm512_loadu_si512((__m512i*)&sf[t]);
				st = _mm512_loadu_si512((__m512i*)&qrr[t]);
				mask = _mm512_cmpeq_epi8_mask(sq, m1_) | _mm512_cmpeq_epi8_mask(st, m1_);
				tmp = _mm512_mask_blend_epi8(_mm512_cmpeq_epi8_mask(sq, st), sc_mis_, sc_mch_);
				tmp = _mm512_mask_blend_epi8(mask, tmp, sc_N_);
				_mm512_storeu_si512((__m512i*)((int8_t*)s + t), tmp);
			}
		} else {
			for (t = st0; t <= en0; ++t)
				((uint8_t*)s)[t] = mat[sf[t] * m + qrr[t]];
		}
		// core loop
		x1_  = _mm512_maskz_set1_epi8(1, x1);
		x21_ = _mm512_maskz_set1_epi8(1, x21);
		v1_  = _mm512_maskz_set1_epi8(1, v1);
		st_ = st / 64, en_ = en / 64;
		assert(en_ - st_ + 1 <= n_col_);
		if (!with_cigar) { // score only
			for (t = st_; t <= en_; ++t) {
				__m512i z, a, b, a2, b2, xt1, x2t1, vt1, ut, tmp;
				__dp_code_block1;
				z = _mm512_max_epi8(z, a);
				z = _mm512_max_epi8(z, b);
				z = _mm512_max_epi8(z, a2);
				z = _mm512_max_epi8(z, b2);
				z = _mm512_min_epi8(z, sc_mch_);
				__dp_code_block2; // save u[] and v[]; update a, b, a2 and b2
				_mm512_store_si512(&x[t],  _mm512_sub_epi8(_mm512_max_epi8(a,  zero_), qe_));
				_mm512_store_si512(&y[t],  _mm512_sub_epi8(_mm512_max_epi8(b,  zero_), qe_));
				_mm512_store_si512(&x2[t], _mm512_sub_epi8(_mm512_max_epi8(a2, zero_), qe2_));
				_mm512_store_si512(&y2[t], _mm512_sub_epi8(_mm512_max_epi8(b2, zero_), qe2_));
			}
		} else if (!(flag&KSW_EZ_RIGHT)) { // gap left-alignment
			__m512i *pr = p + (size_t)r * n_col_ - st_;
			off[r] = st, off_end[r] = en;
			for (t = st_; t <= en_; ++t) {
				__m512i d, z, a, b, a2, b2, xt1, x2t1, vt1, ut, tmp;
				__dp_code_block1;
				d = _mm512_maskz_set1_epi8(_mm512_cmpgt_epi8_mask(a,  z), 1);
				z = _mm512_max_epi8(z, a);
				d = _mm512_mask_blend_epi8(_mm512_cmpgt_epi8_mask(b,  z), d, _mm512_set1_epi8(2));
				z = _mm512_max_epi8(z, b);
				d = _mm512_mask_blend_epi8(_mm512_cmpgt_epi8_mask(a2, z), d, _mm512_set1_epi8(3));
				z = _mm512_max_epi8(z, a2);
				d = _mm512_mask_blend_epi8(_mm512_cmpgt_epi8_mask(b2, z), d, _mm512_set1_epi8(4));
				z = _mm512_max_epi8(z, b2);
				z = _mm512_min_epi8(z, sc_mch_);
				__dp_code_block2;
				d = _mm512_or_si512(d, _mm512_maskz_set1_epi8(_mm512_cmpgt_epi8_mask(a,  zero_), 0x08)); // d = a  > 0? 1<<3 : 0
				_mm512_store_si512(&x[t],  _mm512_sub_epi8(_mm512_max_epi8(a,  zero_), qe_));
				d = _mm512_or_si512(d, _mm512_maskz_set1_epi8(_mm512_cmpgt_epi8_mask(b,  zero_), 0x10)); // d = b  > 0? 1<<4 : 0
				_mm512_store_si512(&y[t],  _mm512_sub_epi8(_mm512_max_epi8(b,  zero_), qe_));
				d = _mm512_or_si512(d, _mm512_maskz_set1_epi8(_mm512_cmpgt_epi8_mask(a2, zero_), 0x20)); // d = a2 > 0? 1<<5 : 0
				_mm512_store_si512(&x2[t], _mm512_sub_epi8(_mm512_max_epi8(a2, zero_), qe2_));
				d = _mm512_or_si512(d, _mm512_maskz_set1_epi8(_mm512_cmpgt_epi8_mask(b2, zero_), 0x40)); // d = b2 > 0? 1<<6 : 0
				_mm512_store_si512(&y2[t], _mm512_sub_epi8(_mm512_max_epi8(b2, zero_), qe2_));
				_mm512_store_si512(&pr[t], d);
			}
		} else { // gap right-alignment
			__m512i *pr = p + (size_t)r * n_col_ - st_;
			off[r] = st, off_end[r] = en;
			for (t = st_; t <= en_; ++t) {
				__m512i d, z, a, b, a2, b2, xt1, x2t1, vt1, ut, tmp;
				__dp_code_block1;
				d = _mm512_maskz_set1_epi8(_mm512_cmpge_epi8_mask(a,  z), 1);
				z = _mm512_max_epi8(z, a);
				d = _mm512_mask_blend_epi8(_mm512_cmpge_epi8_mask(b,  z), d, _mm512_set1_epi8(2));
				z = _mm512_max_epi8(z, b);
				d = _mm512_mask_blend_epi8(_mm512_cmpge_epi8_mask(a2, z), d, _mm512_set1_epi8(3));
				z = _mm512_max_epi8(z, a2);
				d = _mm512_mask_blend_epi8(_mm512_cmpge_epi8_mask(b2, z), d, _mm512_set1_epi8(4));
				z = _mm512_max_epi8(z, b2);
				z = _mm512_min_epi8(z, sc_mch_);
				__dp_code_block2;
				d = _mm512_or_si512(d, _mm512_maskz_set1_epi8(_mm512_cmpge_epi8_mask(a,  zero_), 0x08)); // d = a  >= 0? 1<<3 : 0
				_mm512_store_si512(&x[t],  _mm512_sub_epi8(_mm512_max_epi8(a,  zero_), qe_));
				d = _mm512_or_si512(d, _mm512_maskz_set1_epi8(_mm512_cmpge_epi8_mask(b,  zero_), 0x10)); // d = b  >= 0? 1<<4 : 0
				_mm512_store_si512(&y[t],  _mm512_sub_epi8(_mm512_max_epi8(b,  zero_), qe_));
				d = _mm512_or_si512(d, _mm512_maskz_set1_epi8(_mm512_cmpge_epi8_mask(a2, zero_), 0x20)); // d = a2 >= 0? 1<<5 : 0
				_mm512_store_si512(&x2[t], _mm512_sub_epi8(_mm512_max_epi8(a2, zero_), qe2_));
				d = _mm512_or_si512(d, _mm512_maskz_set1_epi8(_mm512_cmpge_epi8_mask(b2, zero_), 0x40)); // d = b2 >= 0? 1<<6 : 0
				_mm512_store_si512(&y2[t], _mm512_sub_epi8(_mm512_max_epi8(b2, zero_), qe2_));
				_mm512_store_si512(&pr[t], d);
			}
		}
		if (!approx_max) { // find the exact max with a 32-bit score array
			int32_t max_H, max_t;
			// compute H[], max_H and max_t
			if (r > 0) {
				int32_t HH[16], tt[16], en1 = st0 + (en0 - st0) / 16 * 16, i;
				__m512i max_H_, max_t_;
				max_H = H[en0] = en0 > 0? H[en0-1] + u8[en0] : H[en0] + v8[en0]; // special casing the last element
				max_t = en0;
				max_H_ = _mm512_set1_epi32(max_H);
				max_t_ = _mm512_set1_epi32(max_t);
				for (t = st0; t < en1; t += 16) { // this implements: H[t]+=v8[t]; if(H[t]>max_H) max_H=H[t],max_t=t;
					__m512i H1, t_;
					__mmask64 tmp2;
					H1 = _mm512_loadu_si512((__m512i*)&H[t]);
					t_ = _mm512_cvtepi8_epi32(_mm_loadu_si128((__m128i*)&v8[t]));
					H1 = _mm512_add_epi32(H1, t_);
					_mm512_storeu_si512((__m512i*)&H[t], H1);
					t_ = _mm512_set1_epi32(t);
					tmp2 = _mm512_cmpgt_epi32_mask(H1, max_H_);
					max_H_ = _mm512_mask_blend_epi32(tmp2, max_H_, H1);
					max_t_ = _mm512_mask_blend_epi32(tmp2, max_t_, t_);
				}
				_mm512_storeu_si512((__m512i*)HH, max_H_);
				_mm512_storeu_si512((__m512i*)tt, max_t_);
				for (i = 0; i < 16; ++i)
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
			ksw_backtrack(km, 1, rev_cigar, 0, (uint8_t*)p, off, off_end, n_col_*64, tlen-1, qlen-1, &ez->m_cigar, &ez->n_cigar, &ez->cigar);
		} else if (!ez->zdropped && (flag&KSW_EZ_EXTZ_ONLY) && ez->mqe + end_bonus > (int)ez->max) {
			ez->reach_end = 1;
			ksw_backtrack(km, 1, rev_cigar, 0, (uint8_t*)p, off, off_end, n_col_*64, ez->mqe_t, qlen-1, &ez->m_cigar, &ez->n_cigar, &ez->cigar);
		} else if (ez->max_t >= 0 && ez->max_q >= 0) {
			ksw_backtrack(km, 1, rev_cigar, 0, (uint8_t*)p, off, off_end, n_col_*64, ez->max_t, ez->max_q, &ez->m_cigar, &ez->n_cigar, &ez->cigar);
		}
		kfree(km, mem2); kfree(km, off);
	}
}
#endif // __AVX512BW__
