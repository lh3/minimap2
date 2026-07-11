/*************************************************************************************
MIT License

Copyright (c) 2020 Intel Labs

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

Authors: Saurabh Kalikar <saurabh.kalikar@intel.com>; Vasimuddin Md <vasimuddin.md@intel.com>; Sanchit Misra <sanchit.misra@intel.com>; Chirag Jain
*****************************************************************************************/

#ifdef USE_SIMDE
// Portable SIMD: translate AVX2/AVX-512 intrinsics to the host ISA (e.g. ARM NEON)
// via the SIMD-everywhere (SIMDe) library. Vendored patch for minimap2 (mm2-fast port).
#define SIMDE_ENABLE_NATIVE_ALIASES
#include <simde/x86/avx2.h>
#include <simde/x86/avx512.h>
// SIMDe exposes the rounding-mode constants under its SIMDE_ prefix but does not
// always provide the bare x86 aliases the AVX code below relies on.
#ifndef _MM_FROUND_TO_NEG_INF
#define _MM_FROUND_TO_NEG_INF SIMDE_MM_FROUND_TO_NEG_INF
#endif
#ifndef _MM_FROUND_NO_EXC
#define _MM_FROUND_NO_EXC SIMDE_MM_FROUND_NO_EXC
#endif
#else
#include <immintrin.h>
#endif
#include <assert.h>
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <algorithm>
#include <climits>
#include <vector>
#include <map>
using namespace std;


static const char LogTable256_dp_lib[256] = {
#define LT(n) n, n, n, n, n, n, n, n, n, n, n, n, n, n, n, n
	-1, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3,
	LT(4), LT(5), LT(5), LT(6), LT(6), LT(6), LT(6),
	LT(7), LT(7), LT(7), LT(7), LT(7), LT(7), LT(7), LT(7)
};

static inline int ilog2_32_dp_lib(uint32_t v)
{
	uint32_t t, tt;
	if ((tt = v>>16)) return (t = tt>>8) ? 24 + LogTable256_dp_lib[t] : 16 + LogTable256_dp_lib[tt];
	return (t = v>>8) ? 8 + LogTable256_dp_lib[t] : LogTable256_dp_lib[v];
}



class anchor_t {
	public:
		uint64_t r;
		int32_t q;
		int32_t l;
		anchor_t(){}

		anchor_t(uint64_t x, int32_t y, int32_t length){
			r = x; q = y; l = length;
		}	
};

typedef uint32_t num_bits_t;
void create_SoA_Anchors_32_bit(anchor_t* anc, num_bits_t anc_size, num_bits_t* &anchor_r, num_bits_t* &anchor_q, num_bits_t* &anchor_l){
	//num_bits_t anc_size = anc.size();

	anchor_r = (num_bits_t*) malloc((16+anc_size)*sizeof(num_bits_t));
	anchor_q = (num_bits_t*) malloc((16+anc_size)*sizeof(num_bits_t));
	anchor_l = (num_bits_t*) malloc((16+anc_size)*sizeof(num_bits_t));
	// Padding
	anchor_r = &anchor_r[16];
	anchor_q = &anchor_q[16];
	anchor_l = &anchor_l[16];

	for(uint32_t i = 0; i < anc_size; i++){
		anchor_r[i] = (num_bits_t)anc[i].r;
		anchor_q[i] = (num_bits_t)anc[i].q;
		anchor_l[i] = (num_bits_t)anc[i].l;
	}


}







class dp_chain {
	public:
	//Tunable parameters 
	int max_dist_x, max_dist_y, bw, max_skip, max_iter, is_cdna, n_segs;
	float gap_scale;
	int min_cnt;
	int min_sc;
	float chn_pen_gap;
	float chn_pen_skip;

	int n_seg;
 



#ifdef __AVX512BW__
	__m512i zero_v;// = _mm512_setzero_si512();
#elif __AVX2__
	__m256i zero_avx2_v;// = _mm512_setzero_si512();
#endif

	dp_chain(){
		this->max_dist_x = 0;
		this->max_dist_y = 0;
		this->bw = 0;
		this->max_skip = 0;
		this->max_iter = 0;
		this->min_cnt = 0;
		this->min_sc = 0;
		this->chn_pen_gap = 0.0;
		this->chn_pen_skip = 0.0;

		this->is_cdna = 0;
		this->n_seg = 0;
		this->n_segs = 0;
		this->gap_scale = 0;
#ifdef __AVX512BW__
		zero_v = _mm512_setzero_si512();
#elif __AVX2__
		zero_avx2_v = _mm256_setzero_si256();
#endif

	}

	void test(){
		printf("hyper-parameters: %d %d %d %d %d %f %d %d\n",max_dist_x, max_dist_y, bw, max_skip, max_iter, gap_scale, is_cdna, n_segs);
	}


	dp_chain(int max_dist_x, int max_dist_y, int bw, int max_skip, int max_iter,
	int min_cnt,
	int min_sc,
	float chn_pen_gap,
	float chn_pen_skip,
	int is_cdna,
	int n_seg){

		this->max_dist_x = max_dist_x;
		this->max_dist_y = max_dist_y;
		this->bw = bw;
		this->max_skip = max_skip;
		this->max_iter = max_iter;
		this->min_cnt = min_cnt;
		this->min_sc = min_sc;
		this->chn_pen_gap = chn_pen_gap;
		this->chn_pen_skip = chn_pen_skip;

		this->is_cdna = is_cdna;
		this->n_seg = n_seg;
#ifdef __AVX512BW__
		zero_v = _mm512_setzero_si512();
#elif __AVX2__
		zero_avx2_v = _mm256_setzero_si256();
#endif
		this->n_segs = 0;
		this->gap_scale = 0;

	}



	dp_chain(int max_dist_x, int max_dist_y, int bw, int max_skip, int max_iter, float gap_scale, int is_cdna, int n_segs){
 		this->max_dist_x = max_dist_x;
		this->max_dist_y = max_dist_y;
		this->bw = bw;
		this->max_skip = max_skip;
		this->max_iter = max_iter;
		this->gap_scale = gap_scale;
		this->is_cdna = is_cdna;
		this->n_segs = n_segs;
#ifdef __AVX512BW__
		zero_v = _mm512_setzero_si512();
#elif __AVX2__
		zero_avx2_v = _mm256_setzero_si256();
#endif

	}

#ifdef __AVX2__
__m256 mg_log2_avx2(__m256i dd_v) // NB: this doesn't work when x<2
{
	// -------- constant vectors --------------------
	__m256i v255 = _mm256_set1_epi32(255);
	__m256i v128 = _mm256_set1_epi32(128);
	__m256i shift1 =  _mm256_set1_epi32(~(255 << 23));
	__m256i shift2 =  _mm256_set1_epi32(127 << 23);
        __m256 fc1 = _mm256_set1_ps(-0.34484843f);
        __m256 fc2 = _mm256_set1_ps(2.02466578f);
        __m256 fc3 = _mm256_set1_ps(0.67487759f);
	// ---------------------------------------------
	__m256 dd_v_f =_mm256_cvtepi32_ps(dd_v);	
        __m256i dd_v_i = _mm256_castps_si256(dd_v_f);	

	__m256i log2_v_i = _mm256_sub_epi32 (_mm256_and_si256( _mm256_srli_epi32(dd_v_i, 23), v255) , v128);  
		
	dd_v_i = _mm256_and_si256(dd_v_i, shift1);
	dd_v_i = _mm256_add_epi32(dd_v_i, shift2);

	dd_v_f = _mm256_castsi256_ps(dd_v_i); 

	__m256 t1 =_mm256_add_ps (_mm256_mul_ps(fc1, dd_v_f), fc2); 
	__m256 t2 = _mm256_sub_ps(_mm256_mul_ps(t1, dd_v_f), fc3);

	__m256 log2_v_f = _mm256_add_ps(_mm256_cvtepi32_ps(log2_v_i), t2);	

	//print(log2_v_f);

	return log2_v_f;
}

// AVX2
//int32_t comput_sc_vectorized_avx2(uint64_t ai_x, uint64_t ai_y, uint64_t aj_x, uint64_t aj_y){//const mm128_t *ai, const mm128_t *aj) {
__m256i comput_sc_vectorized_avx2(__m256i ai_x_v, __m256i ai_y_v, __m256i aj_x_v, __m256i aj_y_v, __m256i q_span_v, bool flag = false) {

//	uint64_t ai_x, ai_y, aj_x, aj_y;
//	ai_x = 17264860985; ai_y = 64424509728; aj_x = 17264860960; aj_y = 64424509720; // input

/*	int32_t max_dist_x = 5000;
 	int32_t max_dist_y = 5000;
	int32_t bw = 500;
	float chn_pen_gap = 0.120000; 
	float chn_pen_skip = 0.000;  
	int is_cdna = 0; 
	int n_seg = 1;
*/
	int32_t a[16];


	__m256i zero_v = _mm256_set1_epi32(0);
	__m256i max_dist_x_v = _mm256_set1_epi32(max_dist_x);
	__m256i max_dist_y_v = _mm256_set1_epi32(max_dist_y);
	__m256i bw_v = _mm256_set1_epi32(bw);
	__m256i int_min_v = _mm256_set1_epi32(INT32_MIN);
	__m256 chn_pen_gap_v = _mm256_set1_ps(chn_pen_gap);	
	__m256 chn_pen_skip_v = _mm256_set1_ps(chn_pen_skip);	
	

	
	//int32_t dq = (int32_t)ai_y - (int32_t)aj_y, dr, dd, dg, q_span, sc;
	 __m256i dq_v = _mm256_sub_epi32(ai_y_v, aj_y_v);



	//int32_t sidi = 1; //(ai_y & MM_SEED_SEG_MASK) >> MM_SEED_SEG_SHIFT;
	//int32_t sidj = 1; //(aj_y & MM_SEED_SEG_MASK) >> MM_SEED_SEG_SHIFT;

	//if (dq <= 0 || dq > max_dist_x) return INT32_MIN;
//	__mmask16 mask_leq = _mm512_cmple_epi32_mask(dq_v, zero_v);
	__m256i mask_leq1 = _mm256_cmpgt_epi32(zero_v, dq_v);						
	__m256i mask_leq2 = _mm256_cmpeq_epi32(dq_v, zero_v);

//	__mmask16 mask_gt2 = _mm512_cmpgt_epi32_mask(dq_v, max_dist_x_v);
	__m256i mask_gt2 = _mm256_cmpgt_epi32(dq_v, max_dist_x_v);


	//dr = (int32_t)(ai_x - aj_x);
	__m256i dr_v = _mm256_sub_epi32(ai_x_v, aj_x_v);


	//if (/*sidi == sidj &&*/ (dr == 0 || dq > max_dist_y)) return INT32_MIN;
	//__mmask16 mask_eq = _mm512_cmpeq_epi32_mask(dr_v, zero_v);
	__m256i mask_eq = _mm256_cmpeq_epi32(dr_v, zero_v);

	//__mmask16 mask_gt1 = _mm512_cmpgt_epi32_mask(dq_v, max_dist_y_v);	
	__m256i mask_gt1 = _mm256_cmpgt_epi32(dq_v, max_dist_y_v);

	//dd = dr > dq? dr - dq : dq - dr;
	//__m512i dd_v = _mm512_abs_epi32(_mm512_sub_epi32(dr_v, dq_v)); 	
	__m256i dd_v = _mm256_abs_epi32(_mm256_sub_epi32(dr_v, dq_v));

	//if (/*sidi == sidj &&*/ dd > bw) return INT32_MIN;
	//__mmask16 bw_gt = _mm512_cmpgt_epi32_mask(dd_v, bw_v);	
	__m256i bw_gt = _mm256_cmpgt_epi32(dd_v, bw_v);	


	// n_seg == 1 , so always false
	//if (/*n_seg > 1 && !is_cdna && sidi == sidj &&*/ dr > max_dist_y) return INT32_MIN;
	

	//dg = dr < dq? dr : dq;
	//__m512i dg_v = _mm512_min_epi32(dr_v, dq_v);	
	__m256i dg_v = _mm256_min_epi32(dr_v, dq_v);	

	//q_span = aj_y>>32&0xff;
		

	//sc = q_span < dg? q_span : dg;

	//__mmask16 continue_mask = ~(mask_leq | mask_gt2 | mask_eq | mask_gt1 | bw_gt);

	__m256i tmp1 = _mm256_or_si256(bw_gt, mask_eq); 
	__m256i tmp2 = _mm256_or_si256(mask_leq1, mask_leq2); 
	__m256i tmp3 = _mm256_or_si256(mask_gt1, mask_gt2); 

	__m256i continue_mask = _mm256_or_si256(_mm256_or_si256 (tmp1, tmp2), tmp3);	





	//__m512i sc_v = _mm512_maskz_min_epi32(continue_mask, q_span_v, dg_v); //elements are zeroed ut when the corresponding mask bit is "NOT SET"

	//__m256i sc_v = _mm256_andnot_si256(continue_mask,  _mm256_min_epi32(q_span_v, dg_v));
	__m256i sc_v = _mm256_blendv_epi8(_mm256_min_epi32(q_span_v, dg_v), int_min_v, continue_mask);
	
	


	//__mmask16 mask_dd_gt_zero = _mm512_cmpgt_epi32_mask(dd_v, zero_v);	
	//__mmask16 mask_dg_gt_qspan = _mm512_cmpgt_epi32_mask(dg_v, q_span_v);	
	__m256i mask_dd_gt_zero = _mm256_cmpgt_epi32(dd_v, zero_v);	
	__m256i mask_dg_gt_qspan = _mm256_cmpgt_epi32(dg_v, q_span_v);	


	//cout<<continue_mask<<endl;
	continue_mask = (_mm256_andnot_si256(continue_mask, _mm256_or_si256 (mask_dd_gt_zero, mask_dg_gt_qspan)));
	//print(continue_mask);
	//if (dd || dg > q_span)
	{
	//	cout<<"entered into the condition\n";
	//	float lin_pen, log_pen;
	//	lin_pen = chn_pen_gap * (float)dd + chn_pen_skip * (float)dg;
		__m256 lin_pen_v = _mm256_mul_ps(chn_pen_gap_v, _mm256_cvtepi32_ps(dd_v));		

#if 0 //log2 function correctness assertion
		if(mg_log2_v1(dd + 1) != mg_log2(dd + 1)){
			cout<<"Mismatch.\n";
		}
#endif
	//	log_pen = dd >= 1? mg_log2_v1(dd + 1) : 0.0f; // mg_log2() only works for dd>=2
		__m256 log2_v = mg_log2_avx2(_mm256_add_epi32(dd_v , _mm256_set1_epi32(1)));
	

	//	__mmask16 mask_dd_ge_one = _mm512_cmpge_epi32_mask(dd_v, _mm512_set1_epi32(1));
		__m256i mask_dd_ge_one_1 = _mm256_cmpgt_epi32(dd_v, _mm256_set1_epi32(1));
		__m256i mask_dd_ge_one_2 = _mm256_cmpeq_epi32(dd_v, _mm256_set1_epi32(1));
		//__m256i mask_dd_ge_one = _mm256_or_epi32(mask_dd_ge_one_1, mask_dd_ge_one_2);
		__m256i mask_dd_ge_one = _mm256_or_si256(mask_dd_ge_one_1, mask_dd_ge_one_2);
		

		//__m512 log_pen_v = _mm512_mask_blend_ps(mask_dd_ge_one, _mm512_set1_ps(0.0f) ,log2_v);	
		__m256 log_pen_v = _mm256_blendv_ps(_mm256_set1_ps(0.0f) ,log2_v, _mm256_cvtepi32_ps(mask_dd_ge_one));	
#if 0
		if (is_cdna || sidi != sidj) {
			if (sidi != sidj && dr == 0) ++sc; // possibly due to overlapping paired ends; give a minor bonus
			else if (dr > dq || sidi != sidj) sc -= (int)(lin_pen < log_pen? lin_pen : log_pen); // deletion or jump between paired ends
			else sc -= (int)(lin_pen + .5f * log_pen);
		}
		else
#endif
		{
	//		sc -= (int)(lin_pen + .5f * log_pen);
			__m256 mul_t1 = _mm256_mul_ps(_mm256_set1_ps(0.5f), log_pen_v);
	//		print(mul_t1);	

			//__m512i add_t2 = _mm512_cvt_roundps_epi32(_mm256_and_ps(continue_mask, _mm256_add_ps(lin_pen_v, mul_t1)),_MM_FROUND_TO_NEG_INF |_MM_FROUND_NO_EXC );
			//__m256i add_t2 = _mm256_cvtps_epi32(_mm256_round_ps(_mm256_and_ps(_mm256_cvtepi32_ps(continue_mask), _mm256_add_ps(lin_pen_v, mul_t1)),_MM_FROUND_TO_NEG_INF |_MM_FROUND_NO_EXC ));
			__m256i add_t2 = _mm256_cvtps_epi32(_mm256_round_ps(/*_mm256_and_ps(_mm256_cvtepi32_ps(continue_mask),*/ _mm256_add_ps(lin_pen_v, mul_t1)/*)*/,_MM_FROUND_TO_NEG_INF |_MM_FROUND_NO_EXC ));

			add_t2 = _mm256_and_si256(continue_mask, add_t2);
			sc_v = _mm256_sub_epi32(sc_v, add_t2);
		}
	}
	_mm256_store_si256((__m256i *)a, sc_v);
	//sc = a[0];
	//cout<<a[0]<<endl;
//	return a[0];
	return sc_v;	
}
#endif
	


#ifdef __AVX2__

int32_t comput_sc_vectorized_avx2_caller(uint64_t ai_x, uint64_t ai_y, uint64_t aj_x, uint64_t aj_y, int32_t qspan = 0, bool flag = false){//const mm128_t *ai, const mm128_t *aj) {
	#if 0
	int32_t max_dist_x = 5000;
 	int32_t max_dist_y = 5000;
	int32_t bw = 500;
	float chn_pen_gap = 0.120000; 
	float chn_pen_skip = 0.000;  
	int is_cdna = 0; 
	int n_seg = 1;
	#endif
	int32_t a[16];
	__m256i ai_x_v = _mm256_set1_epi32((int32_t)ai_x);
	__m256i ai_y_v = _mm256_set1_epi32((int32_t)ai_y);
	__m256i aj_x_v = _mm256_set1_epi32((int32_t)aj_x);
	__m256i aj_y_v = _mm256_set1_epi32((int32_t)aj_y);
	__m256i q_span_v = _mm256_set1_epi32(qspan);

	__m256i sc_v = comput_sc_vectorized_avx2(ai_x_v, ai_y_v, aj_x_v, aj_y_v, q_span_v, flag);	
	_mm256_store_si256((__m256i *)a, sc_v);
	//sc = a[0];
	//cout<<a[0]<<endl;
	return a[0];
}
#endif



#ifdef __AVX512BW__
void print(__m512i dd_v){
	int32_t a[16];
	_mm512_store_epi32(a, dd_v);

	for(int i = 0; i < 16; i++){
		fprintf(stderr,"%d ", a[i]);// << " ";
	}
		fprintf(stderr,"\n ");// << " ";
	
}
void print(__m512 dd_v){
	float a[16];
	_mm512_store_ps(a, dd_v);

	for(int i = 0; i < 16; i++){
		fprintf(stderr,"%d ", a[i]);// << " ";
	}
		fprintf(stderr,"\n ");// << " ";
}
__m512 mg_log2_v2(__m512i dd_v) // NB: this doesn't work when x<2
{
	// -------- constant vectors --------------------
	__m512i v255 = _mm512_set1_epi32(255);
	__m512i v128 = _mm512_set1_epi32(128);
	__m512i shift1 =  _mm512_set1_epi32(~(255 << 23));
	__m512i shift2 =  _mm512_set1_epi32(127 << 23);
        __m512 fc1 = _mm512_set1_ps(-0.34484843f);
        __m512 fc2 = _mm512_set1_ps(2.02466578f);
        __m512 fc3 = _mm512_set1_ps(0.67487759f);
	// ---------------------------------------------

	__m512 dd_v_f =_mm512_cvtepi32_ps(dd_v);	
        __m512i dd_v_i = _mm512_castps_si512(dd_v_f);	

	__m512i log2_v_i = _mm512_sub_epi32 (_mm512_and_epi32( _mm512_srli_epi32(dd_v_i, 23), v255) , v128);  
	
	dd_v_i = _mm512_and_epi32(dd_v_i, shift1);
	dd_v_i = _mm512_add_epi32(dd_v_i, shift2);

	dd_v_f = _mm512_castsi512_ps(dd_v_i); 

	__m512 t1 =_mm512_add_ps (_mm512_mul_ps(fc1, dd_v_f), fc2); 
	__m512 t2 = _mm512_sub_ps(_mm512_mul_ps(t1, dd_v_f), fc3);

	__m512 log2_v_f = _mm512_add_ps(_mm512_cvtepi32_ps(log2_v_i), t2);	

	//print(log2_v_f);

	return log2_v_f;
}
float mg_log2_v1(float x) // NB: this doesn't work when x<2
{
	__m512 v = mg_log2_v2(_mm512_set1_epi32((int) x));	
	float a[16];
	_mm512_store_ps(a, v);
//	cout<<"log2: "<<a[0]<<endl;
	return a[0];	
}


__m512i comput_sc_vectorized(__m512i ai_x_v, __m512i ai_y_v, __m512i aj_x_v, __m512i aj_y_v, __m512i q_span_v, bool flag = false){

#if 0

	int32_t max_dist_x = 5000;
 	int32_t max_dist_y = 5000;
	int32_t bw = 500;
	float chn_pen_gap = 0.120000; 
	float chn_pen_skip = 0.000;  
	int is_cdna = 0; 
	int n_seg = 1;

//	__m512i zero_v = _mm512_set1_epi32(0);
	__m512i max_dist_x_v = _mm512_set1_epi32(5000);
	__m512i max_dist_y_v = _mm512_set1_epi32(5000);
	__m512i bw_v = _mm512_set1_epi32(500);
	__m512i int_min_v = _mm512_set1_epi32(INT32_MIN);
	__m512 chn_pen_gap_v = _mm512_set1_ps(0.120000);	
	__m512 chn_pen_skip_v = _mm512_set1_ps(0.000);	
#endif
	__m512i max_dist_x_v = _mm512_set1_epi32(max_dist_x);
	__m512i max_dist_y_v = _mm512_set1_epi32(max_dist_y);
	__m512i bw_v = _mm512_set1_epi32(bw);
	__m512i int_min_v = _mm512_set1_epi32(INT32_MIN);
	__m512 chn_pen_gap_v = _mm512_set1_ps(chn_pen_gap);	
	__m512 chn_pen_skip_v = _mm512_set1_ps(chn_pen_skip);	
	
	//int32_t dq = (int32_t)ai_y - (int32_t)aj_y, dr, dd, dg, q_span, sc;
	 __m512i dq_v = _mm512_sub_epi32(ai_y_v, aj_y_v);



	//int32_t sidi = 1; //(ai_y & MM_SEED_SEG_MASK) >> MM_SEED_SEG_SHIFT;
	//int32_t sidj = 1; //(aj_y & MM_SEED_SEG_MASK) >> MM_SEED_SEG_SHIFT;

	//if (dq <= 0 || dq > max_dist_x) return INT32_MIN;
	__mmask16 mask_leq = _mm512_cmple_epi32_mask(dq_v, zero_v);
	__mmask16 mask_gt2 = _mm512_cmpgt_epi32_mask(dq_v, max_dist_x_v);



	//dr = (int32_t)(ai_x - aj_x);
	__m512i dr_v = _mm512_sub_epi32(ai_x_v, aj_x_v);


	//if (/*sidi == sidj &&*/ (dr == 0 || dq > max_dist_y)) return INT32_MIN;
	__mmask16 mask_eq = _mm512_cmpeq_epi32_mask(dr_v, zero_v);
	__mmask16 mask_gt1 = _mm512_cmpgt_epi32_mask(dq_v, dq_v);	


	//dd = dr > dq? dr - dq : dq - dr;
	__m512i dd_v = _mm512_abs_epi32(_mm512_sub_epi32(dr_v, dq_v)); 	

	//if (/*sidi == sidj &&*/ dd > bw) return INT32_MIN;
	__mmask16 bw_gt = _mm512_cmpgt_epi32_mask(dd_v, bw_v);	

	// n_seg == 1 , so always false
	//if (/*n_seg > 1 && !is_cdna && sidi == sidj &&*/ dr > max_dist_y) return INT32_MIN;
	

	//dg = dr < dq? dr : dq;
	__m512i dg_v = _mm512_min_epi32(dr_v, dq_v);	

		

	//sc = q_span < dg? q_span : dg;
	__mmask16 continue_mask = ~(mask_leq | mask_gt2 | mask_eq | mask_gt1 | bw_gt);
	
	__m512i sc_v = _mm512_min_epi32(q_span_v, dg_v);
	sc_v = _mm512_mask_blend_epi32(continue_mask ,int_min_v, sc_v);

	__mmask16 mask_dd_gt_zero = _mm512_cmpgt_epi32_mask(dd_v, zero_v);	
	__mmask16 mask_dg_gt_qspan = _mm512_cmpgt_epi32_mask(dg_v, q_span_v);	

	continue_mask = (continue_mask & (mask_dd_gt_zero | mask_dg_gt_qspan));
	//if (dd || dg > q_span)
	{
		//float lin_pen, log_pen;
		//lin_pen = chn_pen_gap * (float)dd + chn_pen_skip * (float)dg;
		__m512 lin_pen_v = _mm512_mul_ps(chn_pen_gap_v, _mm512_cvtepi32_ps(dd_v));		

#if 0 //log2 function correctness assertion
		if(mg_log2_v1(dd + 1) != mg_log2(dd + 1)){
			cout<<"Mismatch.\n";
		}
#endif
		//log_pen = dd >= 1? mg_log2_v1(dd + 1) : 0.0f; // mg_log2() only works for dd>=2
		__m512 log2_v = mg_log2_v2(_mm512_add_epi32(dd_v , _mm512_set1_epi32(1)));
		__mmask16 mask_dd_ge_one = _mm512_cmpge_epi32_mask(dd_v, _mm512_set1_epi32(1));
		__m512 log_pen_v = _mm512_mask_blend_ps(mask_dd_ge_one, _mm512_set1_ps(0.0f) ,log2_v);	
#if 0
		if (is_cdna || sidi != sidj) {
			if (sidi != sidj && dr == 0) ++sc; // possibly due to overlapping paired ends; give a minor bonus
			else if (dr > dq || sidi != sidj) sc -= (int)(lin_pen < log_pen? lin_pen : log_pen); // deletion or jump between paired ends
			else sc -= (int)(lin_pen + .5f * log_pen);
		}
		else
#endif
		{
			//sc -= (int)(lin_pen + .5f * log_pen);
			__m512 mul_t1 = _mm512_mul_ps(_mm512_set1_ps(0.5f), log_pen_v);

			__m512i add_t2 = _mm512_cvt_roundps_epi32(_mm512_maskz_add_ps(continue_mask, lin_pen_v, mul_t1),_MM_FROUND_TO_NEG_INF |_MM_FROUND_NO_EXC );

		
			sc_v = _mm512_sub_epi32(sc_v, add_t2);
		}
	}
			if(flag) print(sc_v);
	return sc_v;
}
#endif

#ifdef __AVX512BW__
int32_t comput_sc_vectorized_caller(uint64_t ai_x, uint64_t ai_y, uint64_t aj_x, uint64_t aj_y, int32_t qspan = 0, bool flag = false){//const mm128_t *ai, const mm128_t *aj) {
	int32_t max_dist_x = 5000;
 	int32_t max_dist_y = 5000;
	int32_t bw = 500;
	float chn_pen_gap = 0.120000; 
	float chn_pen_skip = 0.000;  
	int is_cdna = 0; 
	int n_seg = 1;
	int32_t a[16];
	__m512i ai_x_v = _mm512_set1_epi32((int32_t)ai_x);
	__m512i ai_y_v = _mm512_set1_epi32((int32_t)ai_y);
	__m512i aj_x_v = _mm512_set1_epi32((int32_t)aj_x);
	__m512i aj_y_v = _mm512_set1_epi32((int32_t)aj_y);
	__m512i q_span_v = _mm512_set1_epi32(qspan);

	__m512i sc_v = comput_sc_vectorized(ai_x_v, ai_y_v, aj_x_v, aj_y_v, q_span_v, flag);	
	_mm512_store_epi32(a, sc_v);
	//sc = a[0];
	//cout<<a[0]<<endl;
	return a[0];
}
#endif

int32_t comput_sc(uint64_t ai_x, uint64_t ai_y, uint64_t aj_x, uint64_t aj_y, int32_t q_sp){//const mm128_t *ai, const mm128_t *aj) {

	//uint64_t ai_x, ai_y, aj_x, aj_y;
	//ai_x = 17264860985; ai_y = 64424509728; aj_x = 17264860960; aj_y = 64424509720;
//	int32_t max_dist_x = 5000;
// 	int32_t max_dist_y = 5000;
//	int32_t bw = 500;
//	float chn_pen_gap = 0.120000; 
//	float chn_pen_skip = 0.000;  
//	int is_cdna = 0; 
//	int n_seg = 1;
	int32_t sidi = 1; //(ai_y & MM_SEED_SEG_MASK) >> MM_SEED_SEG_SHIFT;
	int32_t sidj = 1; //(aj_y & MM_SEED_SEG_MASK) >> MM_SEED_SEG_SHIFT;
	
	int32_t dq = (int32_t)ai_y - (int32_t)aj_y, dr, dd, dg, q_span, sc;
	dr = (int32_t)(ai_x - aj_x);
	dd = dr > dq? dr - dq : dq - dr;
	

	if (/*sidi == sidj &&*/ dd > bw) return INT32_MIN;
	if (/*sidi == sidj &&*/ (dr == 0 || dq > max_dist_y)) return INT32_MIN;
	if (dq <= 0 || dq > max_dist_x) return INT32_MIN;
	if (/*n_seg > 1 && !is_cdna && sidi == sidj &&*/ dr > max_dist_y) return INT32_MIN;
	

	dg = dr < dq? dr : dq;
	q_span = q_sp;//aj_y>>32&0xff;
	sc = q_span < dg? q_span : dg;
	if (dd || dg > q_span) {
		float lin_pen, log_pen;
		lin_pen = chn_pen_gap * (float)dd + chn_pen_skip * (float)dg;
	//	if(mg_log2_v1(dd + 1) != mg_log2(dd + 1)){
	//		cout<<"Mismatch.\n";
	//	}
		log_pen = dd >= 1? mg_log2(dd + 1) : 0.0f; // mg_log2() only works for dd>=2
#if 0
		if (is_cdna || sidi != sidj) {
			if (sidi != sidj && dr == 0) ++sc; // possibly due to overlapping paired ends; give a minor bonus
			else if (dr > dq || sidi != sidj) sc -= (int)(lin_pen < log_pen? lin_pen : log_pen); // deletion or jump between paired ends
			else sc -= (int)(lin_pen + .5f * log_pen);
		}
		else
#endif
		{
			sc -= (int)(lin_pen + .5f * log_pen);
		}
	}
	return sc;
}



	

#ifdef __AVX512BW__
	//Vector code with SoA function parameters 32-bit number representation - avx512
	void mm_dp_vectorized(int64_t n, anchor_t *anchors, uint32_t* anchor_r, uint32_t* anchor_q, uint32_t* anchor_l, uint32_t* &f, int32_t* &p, int32_t* &v, int32_t dr, int32_t dq, int (*gap_cost)(anchor_t a, anchor_t b, void* meta_data), void* meta_data)
	{ 
		uint64_t sum_qspan = 0;	
		float avg_qspan;
//		for (int i = 0; i < n; ++i) sum_qspan += anchor_l[i];
//		avg_qspan = (float)sum_qspan / n;
		int st = 0;	
		
		__m512i dr_v = _mm512_set1_epi32((int64_t)dr);
		__m512i dq_v = _mm512_set1_epi32((int64_t)dq);
		__m512i j_idx_base = _mm512_set_epi32(15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0);		
		int32_t maxfVector_v[16];
		int32_t maxjVector_v[16];
		__m512i neg_one_v = _mm512_set1_epi32((int32_t)-1);

		int64_t max_ii = -1;
		for(int i = 0; i < n; i++){
			int64_t end_j;

			int32_t max_j = -1;
			int32_t max_f = anchor_l[i];

			

			//uint64_t ri = anchors[i].r;
			while (st < i && !(anchors[i].r - anchors[st].r <= (uint32_t)dr)) ++st; //predecessor's position is too far

			// TODO: Minimap specific max_iter parameter
			if (i - st > max_iter) st = i - max_iter; //predecessor's index is too far

			
			int j = i - 1;
			if( !(j - st <= 5))
			{
				//broadcast ri and qi
				__m512i ri_v = _mm512_set1_epi32(anchor_r[i]); 
				__m512i qi_v = _mm512_set1_epi32(anchor_q[i]);
			

				__m512i maxj_v = neg_one_v;
				__m512i maxf_v = _mm512_set1_epi32((int32_t)anchor_l[i]);
				__m512i li_v = maxf_v;

			// 16-way vectorized
			//_mm_prefetch((const char *)(&anchor_r[j - 30]), _MM_HINT_T2);	
			//_mm_prefetch((const char *)(&anchor_q[j - 30]), _MM_HINT_T2);	

			for(j = i - 1; (j - 15) >= st; j = j - 16){
		
				_mm_prefetch((const char *)(&anchor_r[j-60]), _MM_HINT_T0);	
				_mm_prefetch((const char *)(&anchor_q[j-60]), _MM_HINT_T0);	
	
				uint32_t *rj, *qj, *lj; 
				rj = &anchor_r[j-15];	qj = &anchor_q[j-15]; lj = &anchor_l[j-15];

 
				// Load rj and qj
    				__m512i rj_v = _mm512_loadu_si512(rj);
    				__m512i qj_v = _mm512_loadu_si512(qj);
    				__m512i lj_v = _mm512_loadu_si512(lj);
			 

    				__m512i fj_v = _mm512_loadu_si512(&f[j-15]);
				__m512i sc_v = comput_sc_vectorized(ri_v, qi_v, rj_v, qj_v, lj_v);
		//		if (i == 177) {
					//fprintf(stderr, "score vector loop\n");
					//print(sc_v);
					//fprintf(stderr, "score single scalar: %d\n", comput_sc_vectorized_caller(12639904, 27952, 12638743, 48485, 15, false));
					//print(comput_sc_vectorized(ri_v, qi_v, rj_v, qj_v, lj_v,false));
					//print(comput_sc_vectorized(ri_v, qi_v, rj_v, qj_v, lj_v,false));
					 //comput_sc_vectorized_caller(12638743, 29416, 12637582, 27952, 15, true);
		//		}
				sc_v = _mm512_add_epi32(fj_v, sc_v);
				

				// Update Maxf and Maxj
				__mmask16 mask_max_sc = _mm512_cmpgt_epi32_mask(sc_v, maxf_v);
				__m512i j_idx_v = _mm512_add_epi32(j_idx_base, _mm512_set1_epi32(j - 15));

				maxf_v = _mm512_max_epi32(sc_v, maxf_v);
				maxj_v = _mm512_mask_blend_epi32(mask_max_sc, maxj_v, j_idx_v);
			}
			

			if(j >= st)
			{
				uint32_t *rj, *qj, *lj; 
				rj = &anchor_r[j-15];	qj = &anchor_q[j-15]; lj = &anchor_l[j-15];

 
				// Load rj and qj
    				__m512i rj_v = _mm512_loadu_si512(rj);
    				__m512i qj_v = _mm512_loadu_si512(qj);
    				__m512i lj_v = _mm512_loadu_si512(lj);
    				
				// padding 
				__m512i fj_v = _mm512_loadu_si512(&f[j-15]);
				__m512i sc_v = comput_sc_vectorized(ri_v, qi_v, rj_v, qj_v, lj_v);
				
  				int shift = st - (j-15);
				__mmask16 loopContinueMask = 0xFFFF;
				loopContinueMask = loopContinueMask>>(shift);
				
				loopContinueMask = loopContinueMask<<(shift);
				sc_v = _mm512_maskz_add_epi32(loopContinueMask, fj_v, sc_v);
				//if (i == 177) {
					//fprintf(stderr, "score single scalar: %d\n", comput_sc_vectorized_caller(12639904, 27952, 12638743, 48485, 15, false));
					 //comput_sc_vectorized_caller(12638743, 29416, 12637582, 27952, 15, true);
				//	print(sc_v);
					//comput_sc_vectorized(ri_v, qi_v, rj_v, qj_v, lj_v,true);
				//	print(ri_v);
				//	print(qi_v);
				//	print(rj_v);
				//	print(qj_v);
				//	print(lj_v);
				//}

				// Update Maxf and Maxj
				__mmask16 mask_max_sc = _mm512_cmpgt_epi32_mask(sc_v, maxf_v);
				__m512i j_idx_v = _mm512_add_epi32( j_idx_base, _mm512_set1_epi32(j - 15));

				maxf_v = _mm512_max_epi32(sc_v, maxf_v);
				maxj_v = _mm512_mask_blend_epi32(mask_max_sc, maxj_v, j_idx_v);
			 	
					
			}


				_mm512_store_epi32(maxfVector_v, maxf_v);
				_mm512_store_epi32(maxjVector_v, maxj_v);
			
				for(int iter = 15; iter >=0; iter--){
					if(maxfVector_v[iter] > max_f) {
						max_f = maxfVector_v[iter];
						max_j = maxjVector_v[iter];
					}
					else if (maxfVector_v[iter] == max_f){
						max_j = max (max_j, maxjVector_v[iter]);
						if((uint32_t)max_f == anchor_l[i]) max_j = -1;
					}
				}
		
				j = st - 1;	
			}
			else{
			int32_t ri = anchor_r[i], qi = anchor_q[i];
			for(; j >=st; j--){

				int32_t ddr, ddq;

				int32_t rj = anchor_r[j], 
					qj = anchor_q[j]; 
			
				//int32_t score = comput_sc_vectorized_caller(anchors[i].r, anchors[i].q, anchors[j].r, anchors[j].q, anchors[j].l);
				int32_t score = comput_sc(anchors[i].r, anchors[i].q, anchors[j].r, anchors[j].q, anchors[j].l);
				//if (i == 177) {
					//fprintf(stderr, "score scalar %d\n", score);
				//}
				if (score == INT32_MIN) continue;
				score += f[j];//anchor_l[i]; 

				if(score > max_f){
					max_f = score;
					max_j = j;
				}

			}

			}
			end_j = j;	
			int debug_iter = 2057329;
			//if (i == debug_iter) fprintf(stderr, "lisa -- endj: %d max_ii: %d max_f: %d \n", end_j, max_ii, max_f);	
#if 1	
		if (max_ii < 0 || anchors[i].r - anchors[max_ii].r > (int64_t)dr) {
			int32_t max = INT32_MIN;
			max_ii = -1;
			for (j = i - 1; j >= st; --j) {
				if (max < (int32_t)f[j]) max = f[j], max_ii = j;
			}
		}
#endif			
#if 1
			if (max_ii >= 0 && max_ii < end_j) {
				int32_t tmp;
				//tmp = comput_sc_vectorized_caller(anchors[i].r, anchors[i].q, anchors[max_ii].r, anchors[max_ii].q, anchors[max_ii].l);
				tmp = comput_sc(anchors[i].r, anchors[i].q, anchors[max_ii].r, anchors[max_ii].q, anchors[max_ii].l);

				if (/*tmp > 0 &&*/ max_f < int32_t(tmp + f[max_ii])){
					max_f = tmp + f[max_ii], max_j = max_ii;
				}
			}
#endif

			f[i] = max_f; 
			p[i] = max_j;
			v[i] = max_j >= 0 && v[max_j] > max_f? v[max_j] : max_f; // v[] keeps the peak score up to i; f[] is the score ending at i, not always the peak

#if 1

		if (max_ii < 0 || (anchors[i].r - anchors[max_ii].r <= (int64_t)dr && f[max_ii] < f[i]))
			max_ii = i;
		//if (mmax_f < max_f) mmax_f = max_f;
#endif
		}
	}
#elif __AVX2__
	void mm_dp_vectorized(int64_t n, anchor_t *anchors, uint32_t* anchor_r, uint32_t* anchor_q, uint32_t* anchor_l, uint32_t* &f, int32_t* &p, int32_t* &v, int32_t dr, int32_t dq, int (*gap_cost)(anchor_t a, anchor_t b, void* meta_data), void* meta_data)
	{ 
		uint64_t sum_qspan = 0;	
		float avg_qspan;
//		for (int i = 0; i < n; ++i) sum_qspan += anchor_l[i];
//		avg_qspan = (float)sum_qspan / n;
		int st = 0;	
		
		__m256i dr_v = _mm256_set1_epi32((int64_t)dr);
		__m256i dq_v = _mm256_set1_epi32((int64_t)dq);
		//__m512i j_idx_base = _mm512_set_epi32(15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0);		
		
		__m256i j_idx_base = _mm256_set_epi32(7, 6, 5, 4, 3, 2, 1, 0);
		int32_t maxfVector_v[8];
		int32_t maxjVector_v[8];
		__m256i neg_one_v = _mm256_set1_epi32((int32_t)-1);

		int64_t max_ii = -1;
		for(int i = 0; i < n; i++){
			int64_t end_j;

			int32_t max_j = -1;
			int32_t max_f = anchor_l[i];

			

			//uint64_t ri = anchors[i].r;
			while (st < i && !(anchors[i].r - anchors[st].r <= (uint32_t)dr)) ++st; //predecessor's position is too far

			// TODO: Minimap specific max_iter parameter
			if (i - st > max_iter) st = i - max_iter; //predecessor's index is too far

			
			int j = i - 1;
			if(!(j - st <= 5))
			{
			//	fprintf(stderr, "False\n");
				//broadcast ri and qi
				__m256i ri_v = _mm256_set1_epi32(anchor_r[i]); 
				__m256i qi_v = _mm256_set1_epi32(anchor_q[i]);
			

				__m256i maxj_v = neg_one_v;
				__m256i maxf_v = _mm256_set1_epi32((int32_t)anchor_l[i]);
				__m256i li_v = maxf_v;

			// 16-way vectorized
			//_mm_prefetch((const char *)(&anchor_r[j - 30]), _MM_HINT_T2);	
			//_mm_prefetch((const char *)(&anchor_q[j - 30]), _MM_HINT_T2);	

			for(j = i - 1; (j - 7) >= st; j = j - 8){
		
				_mm_prefetch((const char *)(&anchor_r[j-60]), _MM_HINT_T0);	
				_mm_prefetch((const char *)(&anchor_q[j-60]), _MM_HINT_T0);	
	
				uint32_t *rj, *qj, *lj; 
				rj = &anchor_r[j-7];	qj = &anchor_q[j-7]; lj = &anchor_l[j-7];

 
				// Load rj and qj
    				__m256i rj_v = _mm256_loadu_si256((__m256i*)rj);
    				__m256i qj_v = _mm256_loadu_si256((__m256i*)qj);
    				__m256i lj_v = _mm256_loadu_si256((__m256i*)lj);
			 

    				__m256i fj_v = _mm256_loadu_si256((__m256i*)&f[j-7]);
				__m256i sc_v = comput_sc_vectorized_avx2(ri_v, qi_v, rj_v, qj_v, lj_v);
		//		if (i == 177) {
					//fprintf(stderr, "score vector loop\n");
					//print(sc_v);
					//fprintf(stderr, "score single scalar: %d\n", comput_sc_vectorized_caller(12639904, 27952, 12638743, 48485, 15, false));
					//print(comput_sc_vectorized(ri_v, qi_v, rj_v, qj_v, lj_v,false));
					//print(comput_sc_vectorized(ri_v, qi_v, rj_v, qj_v, lj_v,false));
					 //comput_sc_vectorized_caller(12638743, 29416, 12637582, 27952, 15, true);
		//		}
				sc_v = _mm256_add_epi32(fj_v, sc_v);
				

				// Update Maxf and Maxj
				__m256i mask_max_sc = _mm256_cmpgt_epi32(sc_v, maxf_v);
				__m256i j_idx_v = _mm256_add_epi32(j_idx_base, _mm256_set1_epi32(j - 7));

				maxf_v = _mm256_max_epi32(sc_v, maxf_v);
				maxj_v = _mm256_or_si256(_mm256_andnot_si256(mask_max_sc, maxj_v), _mm256_and_si256(mask_max_sc, j_idx_v));

#if 0
				__mmask16 mask_max_sc = _mm512_cmpgt_epi32_mask(sc_v, maxf_v);
				__m512i j_idx_v = _mm512_add_epi32(j_idx_base, _mm512_set1_epi32(j - 15));

				maxf_v = _mm512_max_epi32(sc_v, maxf_v);
				maxj_v = _mm512_mask_blend_epi32(mask_max_sc, maxj_v, j_idx_v);
#endif
			}
			

			if(j >= st)
			{
				uint32_t *rj, *qj, *lj; 
				rj = &anchor_r[j-7];	qj = &anchor_q[j-7]; lj = &anchor_l[j-7];

 
				// Load rj and qj
    				__m256i rj_v = _mm256_loadu_si256((__m256i*)rj);
    				__m256i qj_v = _mm256_loadu_si256((__m256i*)qj);
    				__m256i lj_v = _mm256_loadu_si256((__m256i*)lj);
    				
				// padding 
				__m256i fj_v = _mm256_loadu_si256((__m256i*)&f[j-7]);
				__m256i sc_v = comput_sc_vectorized_avx2(ri_v, qi_v, rj_v, qj_v, lj_v);
				
  				int shift = st - (j-7);

				int32_t msk_ar[8];
				for(int it = 0; it < 8; it++){
					msk_ar[it] = (it < (shift))?0xFFFFFFFF:0;
				}

				//__mmask16 loopContinueMask = 0xFFFF;
				//loopContinueMask = loopContinueMask>>(shift);
				
				//loopContinueMask = loopContinueMask<<(shift);
				__m256i loopContinueMask = _mm256_loadu_si256((__m256i*)msk_ar);

				//sc_v = _mm512_maskz_add_epi32(loopContinueMask, fj_v, sc_v);
				sc_v = _mm256_andnot_si256(loopContinueMask, _mm256_add_epi32(fj_v, sc_v));

				//if (i == 177) {
					//fprintf(stderr, "score single scalar: %d\n", comput_sc_vectorized_caller(12639904, 27952, 12638743, 48485, 15, false));
					 //comput_sc_vectorized_caller(12638743, 29416, 12637582, 27952, 15, true);
				//	print(sc_v);
					//comput_sc_vectorized(ri_v, qi_v, rj_v, qj_v, lj_v,true);
				//	print(ri_v);
				//	print(qi_v);
				//	print(rj_v);
				//	print(qj_v);
				//	print(lj_v);
				//}

				// Update Maxf and Maxj
				
#if 0
				__mmask16 mask_max_sc = _mm512_cmpgt_epi32_mask(sc_v, maxf_v);
				__m512i j_idx_v = _mm512_add_epi32( j_idx_base, _mm512_set1_epi32(j - 15));

				maxf_v = _mm512_max_epi32(sc_v, maxf_v);
				maxj_v = _mm512_mask_blend_epi32(mask_max_sc, maxj_v, j_idx_v);
#endif			 	
				__m256i mask_max_sc = _mm256_cmpgt_epi32(sc_v, maxf_v);
				__m256i j_idx_v = _mm256_add_epi32(j_idx_base, _mm256_set1_epi32(j - 7));

				maxf_v = _mm256_max_epi32(sc_v, maxf_v);
				maxj_v = _mm256_or_si256(_mm256_andnot_si256(mask_max_sc, maxj_v), _mm256_and_si256(mask_max_sc, j_idx_v));
					
			}

			_mm256_store_si256((__m256i*) maxfVector_v, maxf_v);
			_mm256_store_si256((__m256i*) maxjVector_v, maxj_v);
		
			for(int iter = 7; iter >=0; iter--){
				if(maxfVector_v[iter] > max_f) {
					max_f = maxfVector_v[iter];
					max_j = maxjVector_v[iter];
				}
				else if (maxfVector_v[iter] == max_f){
					max_j = max (max_j, maxjVector_v[iter]);
					if(max_f == anchor_l[i]) max_j = -1;
				}
			}


		
				j = st - 1;	
			}
			else{
			int32_t ri = anchor_r[i], qi = anchor_q[i];
			for(; j >=st; j--){

				int32_t ddr, ddq;

				int32_t rj = anchor_r[j], 
					qj = anchor_q[j]; 
			
				//int32_t score = comput_sc_vectorized_avx2_caller(anchors[i].r, anchors[i].q, anchors[j].r, anchors[j].q, anchors[j].l);
				int32_t score = comput_sc(anchors[i].r, anchors[i].q, anchors[j].r, anchors[j].q, anchors[j].l);
				//if (score != score1) {
				//	fprintf(stderr, "score scalar %d %d\n", score, score1);
				//}
				if (score == INT32_MIN) continue;
				 score += f[j];//anchor_l[i]; 

				if(score > max_f){
					max_f = score;
					max_j = j;
				}

			}

			}
			end_j = j;	
			int debug_iter = 2057329;
			//if (i == debug_iter) fprintf(stderr, "lisa -- endj: %d max_ii: %d max_f: %d \n", end_j, max_ii, max_f);	
#if 1	
		if (max_ii < 0 || anchors[i].r - anchors[max_ii].r > (int64_t)dr) {
			int32_t max = INT32_MIN;
			max_ii = -1;
			for (j = i - 1; j >= st; --j) {
				if (max < (int32_t)f[j]) max = f[j], max_ii = j;
			}
		}
#endif			
#if 1
			if (max_ii >= 0 && max_ii < end_j) {
				int32_t tmp;
				//tmp = comput_sc(&a[i], &a[max_ii], max_dist_x, max_dist_y, bw, chn_pen_gap, chn_pen_skip, is_cdna, n_seg);
				//tmp = comput_sc_vectorized_caller(anchors[i].r, anchors[i].q, anchors[max_ii].r, anchors[max_ii].q, anchors[max_ii].l);
				tmp = comput_sc(anchors[i].r, anchors[i].q, anchors[max_ii].r, anchors[max_ii].q, anchors[max_ii].l);
			//	if (i == debug_iter) fprintf(stderr, "lisa: endj : %d max_ii: %d max_f: %d tmp_score: %d\n", end_j, max_ii, max_f, tmp);	

			//	if (i == debug_iter) fprintf(stderr, "tmp : %d max_f: %d tmp: %d f_max_ii: %d\n", tmp, max_f, tmp, f[max_ii]);	
				if (/*tmp > 0 &&*/ max_f < int32_t(tmp + f[max_ii])){
			//	if (i == debug_iter) fprintf(stderr, "lisa: endj : %d max_ii: %d max_f: %d tmp_score: %d\n", end_j, max_ii, max_f, tmp);	

					max_f = tmp + f[max_ii], max_j = max_ii;
			//	if (i == debug_iter) fprintf(stderr, "lisa: endj : %ld max_ii: %ld max_f: %ld tmp_score: %ld\n", end_j, max_ii, max_f, tmp);	

				}
			}
#endif

			f[i] = max_f; 
			p[i] = max_j;
			v[i] = max_j >= 0 && v[max_j] > max_f? v[max_j] : max_f; // v[] keeps the peak score up to i; f[] is the score ending at i, not always the peak

#if 1

		if (max_ii < 0 || (anchors[i].r - anchors[max_ii].r <= (int64_t)dr && f[max_ii] < f[i]))
			max_ii = i;
		//if (mmax_f < max_f) mmax_f = max_f;
#endif
		}
	}
#endif

};
