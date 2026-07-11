/* Native ARM NEON DP chaining kernel (minimap2 mm2-fast port, hand-written 4-wide).
 *
 * This is a faithful port of the mm2-fast AVX2 vectorized chaining kernel
 * (contrib/parallel_chaining_v2_22.h) to native 128-bit NEON (4x int32), avoiding
 * the 8->4 splitting overhead that SIMDe incurs when lowering the AVX2 kernel.
 * It exposes the same C entry point as mm2fast_chain.cpp (mm_dp_vectorized_fill),
 * so it is a drop-in alternative selected by `make neonchain=1`.
 *
 * The scoring is byte-identical to scalar comput_sc() in the routed regime
 * (single-segment, non-cDNA, equal ref/query gap limits); the caller (lchain.c)
 * enforces that. Float ops mirror the AVX2 sequence exactly; build with
 * -ffp-contract=off so no FMA contraction changes rounding.
 */
#include <stdint.h>
#include <stdlib.h>
#include <limits.h>
#include <arm_neon.h>
#include "mmpriv.h"   /* mm128_t, mg_log2() */

// ---- module-level scoring parameters (set per call) ----
static int32_t g_max_dist_x, g_max_dist_y, g_bw, g_max_iter;
static float   g_chn_pen_gap, g_chn_pen_skip;

// mg_log2 approximation, vectorized to match mm2-fast's mg_log2_avx2 (NB: x>=2).
static inline float32x4_t mg_log2_neon(int32x4_t dd_v)
{
	const int32x4_t v255   = vdupq_n_s32(255);
	const int32x4_t v128   = vdupq_n_s32(128);
	const int32x4_t shift1 = vdupq_n_s32(~(255 << 23));
	const int32x4_t shift2 = vdupq_n_s32(127 << 23);
	const float32x4_t fc1 = vdupq_n_f32(-0.34484843f);
	const float32x4_t fc2 = vdupq_n_f32(2.02466578f);
	const float32x4_t fc3 = vdupq_n_f32(0.67487759f);

	float32x4_t dd_v_f = vcvtq_f32_s32(dd_v);
	int32x4_t   dd_v_i = vreinterpretq_s32_f32(dd_v_f);

	int32x4_t log2_v_i = vsubq_s32(
		vandq_s32(vreinterpretq_s32_u32(vshrq_n_u32(vreinterpretq_u32_s32(dd_v_i), 23)), v255),
		v128);

	dd_v_i = vandq_s32(dd_v_i, shift1);
	dd_v_i = vaddq_s32(dd_v_i, shift2);
	dd_v_f = vreinterpretq_f32_s32(dd_v_i);

	float32x4_t t1 = vaddq_f32(vmulq_f32(fc1, dd_v_f), fc2);
	float32x4_t t2 = vsubq_f32(vmulq_f32(t1, dd_v_f), fc3);
	return vaddq_f32(vcvtq_f32_s32(log2_v_i), t2);
}

// Vectorized comput_sc for 4 predecessors; mirrors comput_sc_vectorized_avx2.
static inline int32x4_t comput_sc_vec_neon(int32x4_t ai_x, int32x4_t ai_y,
	int32x4_t aj_x, int32x4_t aj_y, int32x4_t q_span)
{
	const int32x4_t zero = vdupq_n_s32(0);
	const int32x4_t one  = vdupq_n_s32(1);
	const int32x4_t int_min = vdupq_n_s32(INT32_MIN);
	const int32x4_t mdx = vdupq_n_s32(g_max_dist_x);
	const int32x4_t mdy = vdupq_n_s32(g_max_dist_y);
	const int32x4_t bw  = vdupq_n_s32(g_bw);
	const float32x4_t pen_gap = vdupq_n_f32(g_chn_pen_gap);

	int32x4_t dq = vsubq_s32(ai_y, aj_y);
	uint32x4_t m_leq1 = vcltq_s32(dq, zero);          // dq < 0
	uint32x4_t m_leq2 = vceqq_s32(dq, zero);          // dq == 0
	uint32x4_t m_gt2  = vcgtq_s32(dq, mdx);           // dq > max_dist_x

	int32x4_t dr = vsubq_s32(ai_x, aj_x);
	uint32x4_t m_eq  = vceqq_s32(dr, zero);           // dr == 0
	uint32x4_t m_gt1 = vcgtq_s32(dq, mdy);            // dq > max_dist_y

	int32x4_t dd = vabsq_s32(vsubq_s32(dr, dq));
	uint32x4_t bw_gt = vcgtq_s32(dd, bw);             // dd > bw

	int32x4_t dg = vminq_s32(dr, dq);

	uint32x4_t cont = vorrq_u32(vorrq_u32(vorrq_u32(bw_gt, m_eq),
	                                      vorrq_u32(m_leq1, m_leq2)),
	                            vorrq_u32(m_gt1, m_gt2));

	// sc = cont ? INT_MIN : min(q_span, dg)
	int32x4_t sc = vbslq_s32(cont, int_min, vminq_s32(q_span, dg));

	uint32x4_t m_dd_gt0 = vcgtq_s32(dd, zero);
	uint32x4_t m_dg_gt_q = vcgtq_s32(dg, q_span);
	// penalty applies where !cont && (dd>0 || dg>q_span)
	uint32x4_t pen_mask = vbicq_u32(vorrq_u32(m_dd_gt0, m_dg_gt_q), cont);

	float32x4_t lin_pen = vmulq_f32(pen_gap, vcvtq_f32_s32(dd));
	float32x4_t log2v = mg_log2_neon(vaddq_s32(dd, one));
	uint32x4_t m_dd_ge1 = vorrq_u32(vcgtq_s32(dd, one), vceqq_s32(dd, one));
	float32x4_t log_pen = vbslq_f32(m_dd_ge1, log2v, vdupq_n_f32(0.0f));

	float32x4_t half_log = vmulq_f32(vdupq_n_f32(0.5f), log_pen);
	float32x4_t sum = vaddq_f32(lin_pen, half_log);
	int32x4_t add_t2 = vcvtq_s32_f32(vrndmq_f32(sum));   // floor then trunc
	add_t2 = vandq_s32(vreinterpretq_s32_u32(pen_mask), add_t2);
	sc = vsubq_s32(sc, add_t2);
	return sc;
}

// Scalar comput_sc (single-segment), byte-identical to minimap2's comput_sc.
static inline int32_t comput_sc_scalar(uint64_t ai_x, int32_t ai_q,
	uint64_t aj_x, int32_t aj_q, int32_t q_span)
{
	int32_t dq = ai_q - aj_q, dr, dd, dg, sc;
	if (dq <= 0 || dq > g_max_dist_x) return INT32_MIN;
	dr = (int32_t)(ai_x - aj_x);
	if (dr == 0 || dq > g_max_dist_y) return INT32_MIN;
	dd = dr > dq? dr - dq : dq - dr;
	if (dd > g_bw) return INT32_MIN;
	dg = dr < dq? dr : dq;
	sc = q_span < dg? q_span : dg;
	if (dd || dg > q_span) {
		float lin_pen = g_chn_pen_gap * (float)dd + g_chn_pen_skip * (float)dg;
		float log_pen = dd >= 1? mg_log2(dd + 1) : 0.0f;
		sc -= (int)(lin_pen + .5f * log_pen);
	}
	return sc;
}

void mm_dp_vectorized_fill(
	int64_t n, const mm128_t *a,
	int max_dist_x, int max_dist_y, int bw, int max_skip, int max_iter,
	int min_cnt, int min_sc, float chn_pen_gap, float chn_pen_skip,
	int is_cdna, int n_seg,
	int32_t *f, int64_t *p, int32_t *v)
{
	(void)max_skip; (void)min_cnt; (void)min_sc; (void)is_cdna; (void)n_seg;
	g_max_dist_x = max_dist_x; g_max_dist_y = max_dist_y; g_bw = bw;
	g_max_iter = max_iter; g_chn_pen_gap = chn_pen_gap; g_chn_pen_skip = chn_pen_skip;

	// 16-padded uint32 SoA for aligned/safe vector loads below the window start.
	const int PAD = 16;
	uint32_t *ar = (uint32_t*)malloc((PAD + n) * 4);
	uint32_t *aq = (uint32_t*)malloc((PAD + n) * 4);
	uint32_t *al = (uint32_t*)malloc((PAD + n) * 4);
	ar += PAD; aq += PAD; al += PAD;
	for (int64_t i = 0; i < n; ++i) {
		ar[i] = (uint32_t)a[i].x;
		aq[i] = (uint32_t)a[i].y;
		al[i] = a[i].y >> 32 & 0xff;
	}

	int64_t st = 0, max_ii = -1;
	for (int64_t i = 0; i < n; ++i) {
		int32_t max_f = al[i];
		int64_t max_j = -1, end_j;
		while (st < i && (a[i].x >> 32 != a[st].x >> 32 || a[i].x > a[st].x + max_dist_x)) ++st;
		if (i - st > max_iter) st = i - max_iter;

		int64_t j = i - 1;
		if (!(j - st <= 3)) { // vectorized: at least one full 4-lane batch
			int32x4_t ri = vdupq_n_s32((int32_t)ar[i]);
			int32x4_t qi = vdupq_n_s32((int32_t)aq[i]);
			int32x4_t maxf_v = vdupq_n_s32((int32_t)al[i]);
			int32x4_t maxj_v = vdupq_n_s32(-1);
			const int32_t idx_init[4] = {0, 1, 2, 3};
			const int32x4_t idx_base = vld1q_s32(idx_init);

			for (j = i - 1; (j - 3) >= st; j -= 4) {
				int32x4_t rj = vreinterpretq_s32_u32(vld1q_u32(&ar[j-3]));
				int32x4_t qj = vreinterpretq_s32_u32(vld1q_u32(&aq[j-3]));
				int32x4_t lj = vreinterpretq_s32_u32(vld1q_u32(&al[j-3]));
				int32x4_t fj = vld1q_s32(&f[j-3]);
				int32x4_t sc = vaddq_s32(fj, comput_sc_vec_neon(ri, qi, rj, qj, lj));
				uint32x4_t gt = vcgtq_s32(sc, maxf_v);
				int32x4_t jidx = vaddq_s32(idx_base, vdupq_n_s32((int32_t)(j - 3)));
				maxf_v = vmaxq_s32(sc, maxf_v);
				maxj_v = vbslq_s32(gt, jidx, maxj_v);
			}
			if (j >= st) { // masked tail batch (lanes below st zeroed out)
				int32x4_t rj = vreinterpretq_s32_u32(vld1q_u32(&ar[j-3]));
				int32x4_t qj = vreinterpretq_s32_u32(vld1q_u32(&aq[j-3]));
				int32x4_t lj = vreinterpretq_s32_u32(vld1q_u32(&al[j-3]));
				int32x4_t fj = vld1q_s32(&f[j-3]);
				int32x4_t sc = vaddq_s32(fj, comput_sc_vec_neon(ri, qi, rj, qj, lj));
				int shift = (int)(st - (j - 3));
				int32_t msk[4];
				for (int it = 0; it < 4; ++it) msk[it] = it < shift? -1 : 0;
				uint32x4_t contmsk = vld1q_u32((uint32_t*)msk);
				sc = vbicq_s32(sc, vreinterpretq_s32_u32(contmsk)); // zero masked lanes
				uint32x4_t gt = vcgtq_s32(sc, maxf_v);
				int32x4_t jidx = vaddq_s32(idx_base, vdupq_n_s32((int32_t)(j - 3)));
				maxf_v = vmaxq_s32(sc, maxf_v);
				maxj_v = vbslq_s32(gt, jidx, maxj_v);
			}
			int32_t mfv[4], mjv[4];
			vst1q_s32(mfv, maxf_v); vst1q_s32(mjv, maxj_v);
			for (int it = 3; it >= 0; --it) {
				if (mfv[it] > max_f) { max_f = mfv[it]; max_j = mjv[it]; }
				else if (mfv[it] == max_f) {
					if (mjv[it] > max_j) max_j = mjv[it];
					if (max_f == (int32_t)al[i]) max_j = -1;
				}
			}
			j = st - 1;
		} else { // scalar predecessor loop
			for (; j >= st; --j) {
				int32_t sc = comput_sc_scalar(a[i].x, (int32_t)a[i].y, a[j].x, (int32_t)a[j].y, al[j]);
				if (sc == INT32_MIN) continue;
				sc += f[j];
				if (sc > max_f) { max_f = sc; max_j = j; }
			}
		}
		end_j = j;

		if (max_ii < 0 || a[i].x - a[max_ii].x > (int64_t)max_dist_x) {
			int32_t mx = INT32_MIN; max_ii = -1;
			for (j = i - 1; j >= st; --j)
				if (mx < f[j]) mx = f[j], max_ii = j;
		}
		if (max_ii >= 0 && max_ii < end_j) {
			int32_t tmp = comput_sc_scalar(a[i].x, (int32_t)a[i].y, a[max_ii].x, (int32_t)a[max_ii].y, al[max_ii]);
			if (tmp != INT32_MIN && max_f < tmp + f[max_ii])
				max_f = tmp + f[max_ii], max_j = max_ii;
		}
		f[i] = max_f; p[i] = max_j;
		v[i] = max_j >= 0 && v[max_j] > max_f? v[max_j] : max_f;
		if (max_ii < 0 || (a[i].x - a[max_ii].x <= (int64_t)max_dist_x && f[max_ii] < f[i]))
			max_ii = i;
	}
	free(ar - PAD); free(aq - PAD); free(al - PAD);
}
