/* Vectorized DP chaining wrapper (mm2-fast port).
 *
 * This C++ translation unit isolates the C++ vectorized chaining kernel from
 * the Trans-Omics Acceleration Library (Kalikar et al., mm2-fast) and exposes a
 * plain C entry point, mm_dp_vectorized_fill(), callable from lchain.c.
 *
 * The kernel fills the f[] (score), p[] (predecessor) and v[] (peak score)
 * arrays exactly as minimap2's scalar DP-fill loop does, so the caller can then
 * reuse the unchanged mg_chain_backtrack()/compact_a() path. Byte-identical
 * output to scalar minimap2 is the design goal; the caller is responsible for
 * only routing anchor sets where the vectorized formula matches scalar
 * comput_sc() (single-segment, non-cDNA, max_dist_x == max_dist_y).
 *
 * On x86 this compiles to native AVX2/AVX-512; on other ISAs (e.g. Apple ARM
 * NEON) the intrinsics are lowered by SIMDe (see contrib/parallel_chaining_v2_22.h).
 */

#include <cstdint>
#include <cstdlib>

extern "C" {
#include "mmpriv.h"   /* mm128_t, mg_log2() used by the kernel's scalar fallback */
}

#include "parallel_chaining_v2_22.h"

extern "C" void mm_dp_vectorized_fill(
	int64_t n, const mm128_t *a,
	int max_dist_x, int max_dist_y, int bw, int max_skip, int max_iter,
	int min_cnt, int min_sc, float chn_pen_gap, float chn_pen_skip,
	int is_cdna, int n_seg,
	int32_t *f, int64_t *p, int32_t *v)
{
	// Array-of-structs -> the kernel's own Structure-of-Arrays layout.
	anchor_t *anchors = (anchor_t*)malloc(n * sizeof(anchor_t));
	for (int64_t i = 0; i < n; ++i) {
		anchors[i].r = a[i].x;
		anchors[i].q = (int32_t)a[i].y;
		anchors[i].l = a[i].y >> 32 & 0xff; // only 8 bits of q_span are used
	}
	num_bits_t *anchor_r, *anchor_q, *anchor_l;
	create_SoA_Anchors_32_bit(anchors, (num_bits_t)n, anchor_r, anchor_q, anchor_l);

	dp_chain obj(max_dist_x, max_dist_y, bw, max_skip, max_iter, min_cnt, min_sc,
	             chn_pen_gap, chn_pen_skip, is_cdna, n_seg);

	// Kernel writes into pre-allocated arrays (f as uint32_t, p/v as int32_t).
	uint32_t *f1 = (uint32_t*)malloc(n * sizeof(uint32_t));
	int32_t  *p1 = (int32_t*) malloc(n * sizeof(int32_t));
	int32_t  *v1 = (int32_t*) malloc(n * sizeof(int32_t));

	obj.mm_dp_vectorized(n, anchors, anchor_r, anchor_q, anchor_l,
	                     f1, p1, v1, max_dist_x, max_dist_y, NULL, NULL);

	for (int64_t i = 0; i < n; ++i) {
		f[i] = (int32_t)f1[i];
		p[i] = p1[i];
		v[i] = v1[i];
	}

	// create_SoA_Anchors_32_bit() returns pointers offset by 16 for padding.
	anchor_r -= 16; anchor_q -= 16; anchor_l -= 16;
	free(anchor_r); free(anchor_q); free(anchor_l);
	free(anchors); free(f1); free(p1); free(v1);
}
