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
#include "ksw2_extd2_avx.h"

#ifdef __AVX512BW__ 
    #define __dp_code_block1_pcl_512_mm2_fast                                        \
        /*__mmask64 mskc_ = mskc;  */                                       \
        /*if (t == st_)                */                                   \
        __mmask64 mskc_ = (t == st_) ? mskc_ar[(st0 - t*64)/16]:mskc;                    \
        z = _mm512_load_si512(&s[t]);                                   \
        xt1 = _mm512_load_si512(&x[t]);                     /* xt1 <- x[r-1][t..t+15] */ \
        /*  tmp = _mm_srli_si128(xt1, 15);   */                /* tmp <- x[r-1][t+15] */ \
        tmp = _mm512_set1_epi8(((int8_t*)x)[t*64 + 63]); \
	/*  xt1 = _mm_or_si128(_mm_slli_si128(xt1, 1), x1_);*/ /* xt1 <- x[r-1][t-1..t+14] */ \
        shf512a = _mm512_shuffle_epi8(xt1, ind512_slli);                \
        shf512b = _mm512_shuffle_i32x4(shf512a, shf512a, mska);         \
        slli512 = _mm512_mask_blend_epi8(mskb, shf512a, shf512b);       \
        xt1 = _mm512_mask_blend_epi8(mskc_, slli512, x1_);               \
        x1_ = tmp;                                                      \
        vt1 = _mm512_load_si512(&v[t]);                     /* vt1 <- v[r-1][t..t+15] */ \
        /*    tmp = _mm_srli_si128(vt1, 15);             */      /* tmp <- v[r-1][t+15] */ \
        tmp = _mm512_set1_epi8(((int8_t*)v)[t*64 + 63]); \
        /*    vt1 = _mm_or_si128(_mm_slli_si128(vt1, 1), v1_); *//* vt1 <- v[r-1][t-1..t+14] */ \
        shf512a = _mm512_shuffle_epi8(vt1, ind512_slli);                \
        shf512b = _mm512_shuffle_i32x4(shf512a, shf512a, mska);         \
        slli512 = _mm512_mask_blend_epi8(mskb, shf512a, shf512b);       \
        vt1 = _mm512_mask_blend_epi8(mskc_, slli512, v1_);               \
        v1_ = tmp;                                                      \
        a = _mm512_add_epi8(xt1, vt1);                      /* a <- x[r-1][t-1..t+14] + v[r-1][t-1..t+14] */ \
        ut = _mm512_load_si512(&u[t]);                      /* ut <- u[t..t+15] */ \
        b = _mm512_add_epi8(_mm512_load_si512(&y[t]), ut);     /* b <- y[r-1][t..t+15] + u[r-1][t..t+15] */ \
        x2t1= _mm512_load_si512(&x2[t]);                                \
        /*  tmp = _mm_srli_si128(x2t1, 15);*/                           \
        tmp = _mm512_set1_epi8(((int8_t*)x2)[t*64 + 63]); \
        /*  x2t1= _mm_or_si128(_mm_slli_si128(x2t1, 1), x21_); */       \
        shf512a = _mm512_shuffle_epi8(x2t1, ind512_slli);               \
        shf512b = _mm512_shuffle_i32x4(shf512a, shf512a, mska);         \
        slli512 = _mm512_mask_blend_epi8(mskb, shf512a, shf512b);       \
        x2t1 = _mm512_mask_blend_epi8(mskc_, slli512, x21_);             \
        x21_= tmp;                                                      \
        a2= _mm512_add_epi8(x2t1, vt1);                                 \
        b2= _mm512_add_epi8(_mm512_load_si512(&y2[t]), ut);


void ksw_extd2_avx512(void *km, int qlen, const uint8_t *query, int tlen,
                      const uint8_t *target, int8_t m, const int8_t *mat,
                      int8_t q, int8_t e, int8_t q2, int8_t e2, int w,
                      int zdrop, int end_bonus, int flag, ksw_extz_t *ez)
{
    
    // printf("AVX512 call\n");

    __m512i ind512_srli = _mm512_set1_epi8(15);
    __mmask8 msk_srli = 0x0003;
    

    int8_t index[64] __attribute((aligned(64)));

    for (int i=0; i<64; i++) 
        index[i] = i%16 - 1;

    index[0]  = 15;
    index[16] = 31;
    index[32] = 47;
    index[48] = 63;
    
    __m512i shf512a, shf512b, slli512;
    __m512i ind512_slli = _mm512_load_si512((__m512i*) index);
    __mmask8 mska = 0x90;
    __mmask64 mskb = 0x0001000100010000;
    __mmask64 mskc = 0x1;
    __mmask64 mskc_ar[4] = {0x1, 0x10000, 0x100000000, 0x1000000000000};
    
    

    
    __mmask64 msk_ar[5] = {0xFFFF, 0xFFFFFFFF, 0xFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF};
    __mmask64 msk_ar2[5] = {0xFFFF, 0xFFFF, 0xFFFFFFFF, 0xFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF};
    
    int r, t, qe = q + e, n_col_, *off = 0, *off_end = 0, tlen_, qlen_, last_st, last_en, wl, wr, max_sc, min_sc, long_thres, long_diff;
    int with_cigar = !(flag&KSW_EZ_SCORE_ONLY), approx_max = !!(flag&KSW_EZ_APPROX_MAX);
    int32_t *H = 0, H0 = 0, last_H0_t = 0;
    uint8_t *qr, *sf, *mem, *mem2 = 0;
    
    __m512i q_, q2_, qe_, qe2_, zero_, sc_mch_, sc_mis_, m1_, sc_N_;
    __m512i *u, *v, *x, *y, *x2, *y2, *s, *p = 0;
    __m512i one_, two_, three_, four_, s1_, s2_, s3_, s4_;
        
    ksw_reset_extz(ez);
    if (m <= 1 || qlen <= 0 || tlen <= 0) return;

    if (q2 + e2 < q + e) t = q, q = q2, q2 = t, t = e, e = e2, e2 = t; // make sure q+e no larger than q2+e2
    s1_   = _mm512_set1_epi8(0x08);
    s2_   = _mm512_set1_epi8(0x10);
    s3_   = _mm512_set1_epi8(0x20);
    s4_   = _mm512_set1_epi8(0x40);
    
    one_   = _mm512_set1_epi8(1);
    two_   = _mm512_set1_epi8(2);
    three_   = _mm512_set1_epi8(3);
    four_   = _mm512_set1_epi8(4);
    __m512i four32_   = _mm512_set1_epi32(4);
    __m512i bt32_     = _mm512_setr_epi32(0,0,0,0,4,4,4,4,8,8,8,8,12,12,12,12);
    
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
    tlen_ = (tlen + 63) / 64;
    n_col_ = qlen < tlen? qlen : tlen;
    n_col_ = ((n_col_ < w + 1? n_col_ : w + 1) + 63) / 64 + 1;
    qlen_ = (qlen + 63) / 64;
    for (t = 1, max_sc = mat[0], min_sc = mat[1]; t < m * m; ++t) {
        max_sc = max_sc > mat[t]? max_sc : mat[t];
        min_sc = min_sc < mat[t]? min_sc : mat[t];
    }
    if (-min_sc > 2 * (q + e)) return; // otherwise, we won't see any mismatches

    long_thres = e != e2? (q2 - q) / (e - e2) - 1 : 0;
    if (q2 + e2 + long_thres * e2 > q + e + long_thres * e)
        ++long_thres;
    long_diff = long_thres * (e - e2) - (q2 - q) - e2;

    mem = (uint8_t*)kcalloc(km, tlen_ * 8 + qlen_ + 1 + 63, 64); //TODO: saurabh - 63
    u = (__m512i*)(((size_t)mem + 63) >> 6 << 6); // 16-byte aligned //TODO: saurabh mm_malloc
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
        // printf("qlen: %dm tlen: %dm ncol: %d\n", qlen, tlen, n_col_);
        mem2 = (uint8_t*)kmalloc(km, ((size_t)(qlen + tlen - 1) * n_col_ + 1) * 64);
        p = (__m512i*)(((size_t)mem2 + 63) >> 6 << 6);//TODO: saurabh
        off = (int*)kmalloc(km, (qlen + tlen - 1) * sizeof(int) * 2);
        off_end = off + qlen + tlen - 1;
    }

    
    #define  ORIG 0             // for B1 block optimization
    for (t = 0; t < qlen; ++t) {
        #if ORIG
        qr[t] = query[qlen - 1 - t];
        #else
        if (query[qlen - 1 - t] == 4)
            qr[t] = 8;
        else
            qr[t] = query[qlen - 1 - t];
        #endif
        //assert(qr[t] < 5);
    }
    memcpy(sf, target, tlen);
    
    
    #if TIMER_O
    uint64_t tic = _rdtsc();
    #endif

    #if 1     // changes due to B1 modification
    int8_t cN   = mat[m*m-1] == 0? -e2 : mat[m*m-1];
    int8_t pmat[64] __attribute((aligned(64))) = {0};
    pmat[0] = mat[0];
    pmat[1] = pmat[2] = pmat[3] = mat[1];
    pmat[4] = pmat[5] = pmat[6] = pmat[7] = pmat[8] = pmat[9] = pmat[10] = pmat[11] = pmat[12] = cN;
    __m512i vpmat = _mm512_load_si512(pmat);
    vpmat = _mm512_shuffle_i32x4(vpmat, vpmat, 0x0000);
    #endif

    static int itr = 0;
    itr++;
    
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
        {
          //  _mm_prefetch(sf + st0, _MM_HINT_T0);
          //  _mm_prefetch(qrr + st0, _MM_HINT_T0);
          //  _mm_prefetch(s + st0, _MM_HINT_T0);
          //  _mm_prefetch(sf + st0 + 64, _MM_HINT_T0);
          //  _mm_prefetch(qrr + st0 + 64, _MM_HINT_T0);
          //  _mm_prefetch(s + st0 + 64, _MM_HINT_T0);
        }
        
        int st_new = st / 16 * 16, en_new = (en + 16) / 16 * 16 - 1;//TODO: saurabh
        // int st_new = st / 64 * 64, en_new = (en + 64) / 64 * 64 - 1;//TODO: saurabh
        int stb = st, enb = en;
        st = st / 64 * 64, en = (en + 64) / 64 * 64 - 1;//TODO: saurabh
        int stn = stb / 16 * 16, enn = (enb + 16) / 16 * 16 - 1;//TODO: saurabh
        // set boundary conditions
        //if(st_new != st)
	{
        	//printf("st0 %d en0 %d st16 %d en16 %d st64 %d en64 %d\n", st0, en0, st_new, en_new, st, en);
		//processed_vect64 += ceil(((double)(en - st) + 1)/64);                                                                                             possible_vect64 += ceil(((double)(en_new - st_new) + 1)/64);
		//exit(0);
	}
        if (st_new > 0) {
            if (st_new - 1 >= last_st && st_new - 1 <= last_en) {
                x1 = x8[st_new - 1], x21 = x28[st_new - 1], v1 = v8[st_new - 1]; // (r-1,s-1) calculated in the last round
            } else {
                x1 = -q - e, x21 = -q2 - e2;
                v1 = -q - e;
            }
        } else {
            x1 = -q - e, x21 = -q2 - e2;
            v1 = r == 0? -q - e : r < long_thres? -e : r == long_thres? long_diff : -e2;
        }
        if (en_new >= r) {
            ((int8_t*)y)[r] = -q - e, ((int8_t*)y2)[r] = -q2 - e2;
            u8[r] = r == 0? -q - e : r < long_thres? -e : r == long_thres? long_diff : -e2;
        }
        #if TIMER_I
        uint64_t tic = _rdtsc();
        #endif

     
   
        // loop fission: set scores first //TODO: saurabh: introducing mask64
        if (!(flag & KSW_EZ_GENERIC_SC)) {
            #if ORIG
            for (t = st0; t <= en0; t += 64) {
                __m512i sq, st, tmp_512, mask_512;
                __mmask64 tmp, mask;
                sq = _mm512_loadu_si512((__m512i*)&sf[t]);
                st = _mm512_loadu_si512((__m512i*)&qrr[t]);
                mask = (_mm512_cmpeq_epi8_mask(sq, m1_) | _mm512_cmpeq_epi8_mask(st, m1_));
                tmp = _mm512_cmpeq_epi8_mask(sq, st);
                
                #if 1
                tmp_512 = _mm512_mask_blend_epi8(tmp, sc_mis_, sc_mch_);
                tmp_512 = _mm512_mask_blend_epi8(mask, tmp_512, sc_N_);
                #else
                tmp = _mm_or_si128(_mm_andnot_si128(tmp,  sc_mis_), _mm_and_si128(tmp,  sc_mch_));//TODO: saurabh - left unchanged
                tmp = _mm_or_si128(_mm_andnot_si128(mask, tmp),     _mm_and_si128(mask, sc_N_));
                #endif
                if (t + 64 > en0)
                {
                    __mmask64 msk;
                    int ind = (en0 - t + 16) >> 4;
                    // assert(ind >= 0 && ind < 5);
                    msk = msk_ar2[ind];
                    _mm512_mask_storeu_epi8((__m512i*)((int8_t*)s + t), msk, tmp_512);
                }
                else
                    _mm512_storeu_si512((__m512i*)((int8_t*)s + t), tmp_512);
            }
            #else
            for (t = 0; t <= (en0 - st0) - 64; t += 64) {
                __m512i sq, st, tmp_512, mask_512;
                __mmask64 tmp, mask;
                // const int pdf = 4;
                //_mm_prefetch(sf + t + 64*pdf, _MM_HINT_T0);
                //_mm_prefetch(qrr + t + 64*pdf, _MM_HINT_T0);
                //_mm_prefetch(s + t + 64*pdf, _MM_HINT_T0);
                
                sq = _mm512_loadu_si512((__m512i*)&sf[st0 + t]);
                st = _mm512_loadu_si512((__m512i*)&qrr[st0 + t]);
                tmp_512 = _mm512_xor_si512(sq, st);
                tmp_512 = _mm512_shuffle_epi8(vpmat, tmp_512);                
                _mm512_storeu_si512((__m512i*)((int8_t*)s + st0 + t), tmp_512);
            }
            {
                __m512i sq, st, tmp_512, mask_512;
                __mmask64 tmp, mask;
                sq = _mm512_loadu_si512((__m512i*)&sf[st0 + t]);
                st = _mm512_loadu_si512((__m512i*)&qrr[st0 + t]);
                tmp_512 = _mm512_xor_si512(sq, st);
                tmp_512 = _mm512_shuffle_epi8(vpmat, tmp_512);
                
                if (t + 64 > en0 - st0)
                {
                    __mmask64 msk;
                    int ind = ((en0 - st0) - t + 16) >> 4;
                    msk = msk_ar2[ind];
                    _mm512_mask_storeu_epi8((__m512i*)((int8_t*)s + st0 + t), msk, tmp_512);
                }
                else
                    _mm512_storeu_si512((__m512i*)((int8_t*)s + st0 + t), tmp_512);
            }
            #endif
        }
        else
        {
        //    printf("ELSE\n");
            for (t = st0; t <= en0; ++t) {
                #if ORIG
                ((uint8_t*)s)[t] = mat[sf[t] * m + qrr[t]];
                #else
                int8_t c = qrr[t];
                if (qrr[t] == 8) {
                    c = 4;
                }
                ((uint8_t*)s)[t] = mat[sf[t] * m + c];
                #endif
            }
        }
	

        #if TIMER_I
        uint64_t toc = _rdtsc();
        tprof[B1][0] += toc - tic;
        #endif
        //continue;
        // core loop //TODO: saurabh - need more understanding of logic
        // fprintf(stderr, "- r: %d, x1: %d, x21: %d, v1: %d, en_new: %d, e: %d, q: %d\n",
        //r, x1, x21, v1, en_new, e, q);
        x1_  = _mm512_set1_epi8((uint8_t)x1);
        x21_ = _mm512_set1_epi8((uint8_t)x21);
        v1_  = _mm512_set1_epi8((uint8_t)v1);
        
        //st_ = st / 16, en_ = en / 16;     
        st_ = st / 64, en_ = en / 64;
        //assert(en_ - st_ + 1 <= n_col_);
        if (!with_cigar) { // score only
            #if TIMER_I
            uint64_t tic = _rdtsc();
            #endif
            for (t = st_; t <= en_; ++t) {
                __m512i z, a, b, a2, b2, xt1, x2t1, vt1, ut, tmp;

                //__dp_code_block1_pcl;
                __dp_code_block1_pcl_512_mm2_fast;
                z = _mm512_max_epi8(z, a);
                z = _mm512_max_epi8(z, b);
                z = _mm512_max_epi8(z, a2);
                z = _mm512_max_epi8(z, b2);
                z = _mm512_min_epi8(z, sc_mch_);
                // __dp_code_block2_pcl; // save u[] and v[]; update a, b, a2 and b2
                #if 1
                if (t == en_) {
                    __mmask64 msk;
                    int ind = (en0 - t*64 + 16)/16;
                    msk = msk_ar2[ind];
                    // fprintf(stderr, "en0: %d, t: %d, ind: %d, msk: %d\n", en0, t, ind, msk);
                    _mm512_mask_storeu_epi8(&u[t], msk, _mm512_sub_epi8(z, vt1));
                    _mm512_mask_storeu_epi8(&v[t], msk, _mm512_sub_epi8(z, ut));
                    tmp = _mm512_sub_epi8(z, q_);                   
                    a = _mm512_sub_epi8(a, tmp);                    
                    b = _mm512_sub_epi8(b, tmp);                    
                    tmp = _mm512_sub_epi8(z, q2_);                  
                    a2= _mm512_sub_epi8(a2, tmp);                   
                    b2= _mm512_sub_epi8(b2, tmp);
                    
                }
                else {
                    _mm512_store_si512(&u[t], _mm512_sub_epi8(z, vt1));
                    _mm512_store_si512(&v[t], _mm512_sub_epi8(z, ut));
                    tmp = _mm512_sub_epi8(z, q_);                   
                    a = _mm512_sub_epi8(a, tmp);                    
                    b = _mm512_sub_epi8(b, tmp);                    
                    tmp = _mm512_sub_epi8(z, q2_);                  
                    a2= _mm512_sub_epi8(a2, tmp);                   
                    b2= _mm512_sub_epi8(b2, tmp);
                }
                #endif
                
                #if 1
                if (t == en_) {
                    __mmask64 msk;
                    int ind = (en0 - t*64 + 16)/16;
                    //assert(ind >= 0);
                    msk = msk_ar2[ind];
                    // fprintf(stderr, "en0: %d, t: %d, ind: %d, msk: %d\n", en0, t, ind, msk);
                    _mm512_mask_storeu_epi8(&x[t], msk, _mm512_sub_epi8(_mm512_max_epi8(a,  zero_), qe_));
                    _mm512_mask_storeu_epi8(&y[t], msk,  _mm512_sub_epi8(_mm512_max_epi8(b,  zero_), qe_));
                    _mm512_mask_storeu_epi8(&x2[t], msk,  _mm512_sub_epi8(_mm512_max_epi8(a2, zero_), qe2_));
                    _mm512_mask_storeu_epi8(&y2[t], msk, _mm512_sub_epi8(_mm512_max_epi8(b2, zero_), qe2_));                    
                }
                else
                #endif
                {
                    _mm512_store_si512(&x[t],  _mm512_sub_epi8(_mm512_max_epi8(a,  zero_), qe_));
                    _mm512_store_si512(&y[t],  _mm512_sub_epi8(_mm512_max_epi8(b,  zero_), qe_));
                    _mm512_store_si512(&x2[t], _mm512_sub_epi8(_mm512_max_epi8(a2, zero_), qe2_));
                    _mm512_store_si512(&y2[t], _mm512_sub_epi8(_mm512_max_epi8(b2, zero_), qe2_));
                }
                // for (int l=0; l<64; l++)
                //    fprintf(stderr, "%d ", ((int8_t*)x)[l]);

            }
            #if TIMER_I
            uint64_t toc = _rdtsc();
            tprof[B2][0] += toc - tic;
            #endif
            
        } else if (!(flag&KSW_EZ_RIGHT)) { // gap left-alignment
            __m512i *pr = p + (size_t)r * n_col_ - st_;
            off[r] = st, off_end[r] = en;
            
           	//_mm_prefetch(pr, _MM_HINT_T0); 
            for (t = st_; t < en_; ++t) {
                __m512i d, z, a, b, a2, b2, xt1, x2t1, vt1, ut, tmp;
                __mmask64 tmp_mask;
                //__dp_code_block1_pcl;                
		__dp_code_block1_pcl_512_mm2_fast;
                d = _mm512_mask_blend_epi8(_mm512_cmpgt_epi8_mask(a, z),zero_,  one_);       // d = a  > z? 1 : 0
                z = _mm512_max_epi8(z, a);              
                d = _mm512_mask_blend_epi8(_mm512_cmpgt_epi8_mask(b,  z), d, two_); // d = b  > z? 2 : d
                z = _mm512_max_epi8(z, b);
                d = _mm512_mask_blend_epi8(_mm512_cmpgt_epi8_mask(a2, z), d, three_); // d = a2 > z? 3 : d
                z = _mm512_max_epi8(z, a2);
                d = _mm512_mask_blend_epi8(_mm512_cmpgt_epi8_mask(b2, z), d, four_); // d = b2 > z? 4 : d
                z = _mm512_max_epi8(z, b2);
                z = _mm512_min_epi8(z, sc_mch_);
                // __dp_code_block2_pcl;                
                    _mm512_store_si512(&u[t], _mm512_sub_epi8(z, vt1));
                    _mm512_store_si512(&v[t], _mm512_sub_epi8(z, ut));
                    tmp = _mm512_sub_epi8(z, q_);                   
                    a = _mm512_sub_epi8(a, tmp);                    
                    b = _mm512_sub_epi8(b, tmp);                    
                    tmp = _mm512_sub_epi8(z, q2_);                  
                    a2= _mm512_sub_epi8(a2, tmp);                   
                    b2= _mm512_sub_epi8(b2, tmp);
                    tmp_mask = _mm512_cmpgt_epi8_mask(a, zero_);
                    _mm512_store_si512(&x[t],  _mm512_sub_epi8(_mm512_mask_blend_epi8(tmp_mask, zero_, a),  qe_));
                    d = _mm512_or_si512(d, _mm512_mask_blend_epi8(tmp_mask, zero_, s1_)); // d = a > 0? 1<<3 : 0
                    tmp_mask = _mm512_cmpgt_epi8_mask(b, zero_);
                    _mm512_store_si512(&y[t],  _mm512_sub_epi8(_mm512_mask_blend_epi8(tmp_mask, zero_, b),  qe_));
                    d = _mm512_or_si512(d, _mm512_mask_blend_epi8(tmp_mask, zero_, s2_)); // d = b > 0? 1<<4 : 0
                    tmp_mask = _mm512_cmpgt_epi8_mask(a2, zero_);
                    _mm512_store_si512(&x2[t], _mm512_sub_epi8(_mm512_mask_blend_epi8(tmp_mask, zero_, a2), qe2_));
                    d = _mm512_or_si512(d, _mm512_mask_blend_epi8(tmp_mask, zero_, s3_)); // d = a > 0? 1<<5 : 0
                    tmp_mask = _mm512_cmpgt_epi8_mask(b2, zero_);
                    _mm512_store_si512(&y2[t], _mm512_sub_epi8(_mm512_mask_blend_epi8(tmp_mask, zero_, b2), qe2_));
                    d = _mm512_or_si512(d, _mm512_mask_blend_epi8(tmp_mask, zero_, s4_)); // d = b > 0? 1<<6 : 0             
                    _mm512_store_si512(&pr[t], d);
           	    //_mm_prefetch(&pr[t+1], _MM_HINT_T0); 
            }
            #if TIMER_I
            uint64_t tic = _rdtsc();
            #endif
            {
                __m512i d, z, a, b, a2, b2, xt1, x2t1, vt1, ut, tmp;
                __mmask64 tmp_mask;
                //__dp_code_block1_pcl;                
		__dp_code_block1_pcl_512_mm2_fast;
                d = _mm512_mask_blend_epi8(_mm512_cmpgt_epi8_mask(a, z),zero_,  one_);       // d = a  > z? 1 : 0
                z = _mm512_max_epi8(z, a);              
                d = _mm512_mask_blend_epi8(_mm512_cmpgt_epi8_mask(b,  z), d, two_); // d = b  > z? 2 : d
                z = _mm512_max_epi8(z, b);
                d = _mm512_mask_blend_epi8(_mm512_cmpgt_epi8_mask(a2, z), d, three_); // d = a2 > z? 3 : d
                z = _mm512_max_epi8(z, a2);
                d = _mm512_mask_blend_epi8(_mm512_cmpgt_epi8_mask(b2, z), d, four_); // d = b2 > z? 3 : d
                z = _mm512_max_epi8(z, b2);
                z = _mm512_min_epi8(z, sc_mch_);
                // __dp_code_block2_pcl;                
                #if 1
                {
                    __mmask64 msk;
                    int ind = (en0 - t*64 + 16)/16;
                  //  assert(ind >= 0 && ind < 5);
                    msk = msk_ar2[ind];
                    
                    _mm512_mask_storeu_epi8(&u[t], msk, _mm512_sub_epi8(z, vt1));
                    _mm512_mask_storeu_epi8(&v[t], msk, _mm512_sub_epi8(z, ut));
                    tmp = _mm512_sub_epi8(z, q_);                   
                    a = _mm512_sub_epi8(a, tmp);                    
                    b = _mm512_sub_epi8(b, tmp);                    
                    tmp = _mm512_sub_epi8(z, q2_);                  
                    a2= _mm512_sub_epi8(a2, tmp);                   
                    b2= _mm512_sub_epi8(b2, tmp);
                }
                #endif
                
                #if 1
                {
                    __mmask64 msk;
                    int ind = (en0 - t*64 + 16)/16;
                    msk = msk_ar2[ind];
                    off_end[r] -= (4-ind)*16;
                    
                    tmp_mask = _mm512_cmpgt_epi8_mask(a, zero_);
                    _mm512_mask_storeu_epi8(&x[t], msk, _mm512_sub_epi8(_mm512_mask_blend_epi8(tmp_mask, zero_, a), qe_));
                    d = _mm512_or_si512(d, _mm512_mask_blend_epi8(tmp_mask, zero_, s1_)); // d = a > 0? 1<<3 : 0
                    
                    tmp_mask = _mm512_cmpgt_epi8_mask(b, zero_);
                    _mm512_mask_storeu_epi8(&y[t], msk, _mm512_sub_epi8(_mm512_mask_blend_epi8(tmp_mask, zero_, b), qe_));
                    d = _mm512_or_si512(d, _mm512_mask_blend_epi8(tmp_mask, zero_, s2_)); // d = b > 0? 1<<4 : 0                    
                    tmp_mask = _mm512_cmpgt_epi8_mask(a2, zero_);
                    _mm512_mask_storeu_epi8(&x2[t], msk, _mm512_sub_epi8(_mm512_mask_blend_epi8(tmp_mask, zero_, a2), qe2_));
                    d = _mm512_or_si512(d, _mm512_mask_blend_epi8(tmp_mask, zero_, s3_)); // d = a > 0? 1<<5 : 0                    
                    tmp_mask = _mm512_cmpgt_epi8_mask(b2, zero_);
                    _mm512_mask_storeu_epi8(&y2[t], msk, _mm512_sub_epi8(_mm512_mask_blend_epi8(tmp_mask, zero_, b2), qe2_));
                    d = _mm512_or_si512(d, _mm512_mask_blend_epi8(tmp_mask, zero_, s4_)); // d = b > 0? 1<<6 : 0
                    //_mm512_store_si512(&pr[t], d);
                    _mm512_mask_storeu_epi8(&pr[t], msk, d);
                }
                #endif    
            }
            #if TIMER_I
            uint64_t toc = _rdtsc();
            tprof[B3][0]+= toc - tic;
            #endif
/*
            for (t = st_; t <= en_; ++t) {

                __m512i d, z, a, b, a2, b2, xt1, x2t1, vt1, ut, tmp;
                __mmask64 tmp_mask;
                __dp_code_block1_pcl;                
                
                d = _mm512_mask_blend_epi8(_mm512_cmpgt_epi8_mask(a, z), _mm512_set1_epi8(0),  _mm512_set1_epi8(1));       // d = a  > z? 1 : 0
                z = _mm512_max_epi8(z, a);              
                d = _mm512_mask_blend_epi8(_mm512_cmpgt_epi8_mask(b,  z), d, _mm512_set1_epi8(2)); // d = b  > z? 2 : d
                z = _mm512_max_epi8(z, b);
                d = _mm512_mask_blend_epi8(_mm512_cmpgt_epi8_mask(a2, z), d, _mm512_set1_epi8(3)); // d = a2 > z? 3 : d
                z = _mm512_max_epi8(z, a2);
                d = _mm512_mask_blend_epi8(_mm512_cmpgt_epi8_mask(b2, z), d, _mm512_set1_epi8(4)); // d = b2 > z? 3 : d
                z = _mm512_max_epi8(z, b2);
                z = _mm512_min_epi8(z, sc_mch_);
                // __dp_code_block2_pcl;                
                #if 1
                if (t == en_) {
                    __mmask64 msk;
                    int ind = (en0 - t*64 + 16)/16;
                    //assert(ind >= 0 && ind < 5);
                    msk = msk_ar2[ind];
                    
                    _mm512_mask_storeu_epi8(&u[t], msk, _mm512_sub_epi8(z, vt1));
                    _mm512_mask_storeu_epi8(&v[t], msk, _mm512_sub_epi8(z, ut));
                    tmp = _mm512_sub_epi8(z, q_);                   
                    a = _mm512_sub_epi8(a, tmp);                    
                    b = _mm512_sub_epi8(b, tmp);                    
                    tmp = _mm512_sub_epi8(z, q2_);                  
                    a2= _mm512_sub_epi8(a2, tmp);                   
                    b2= _mm512_sub_epi8(b2, tmp);

                }
                else {
                    _mm512_store_si512(&u[t], _mm512_sub_epi8(z, vt1));
                    _mm512_store_si512(&v[t], _mm512_sub_epi8(z, ut));
                    tmp = _mm512_sub_epi8(z, q_);                   
                    a = _mm512_sub_epi8(a, tmp);                    
                    b = _mm512_sub_epi8(b, tmp);                    
                    tmp = _mm512_sub_epi8(z, q2_);                  
                    a2= _mm512_sub_epi8(a2, tmp);                   
                    b2= _mm512_sub_epi8(b2, tmp);
                }
                #endif
                
                #if 1
                if (t == en_) {
                    __mmask64 msk;
                    int ind = (en0 - t*64 + 16) >> 4;
                    msk = msk_ar2[ind];
                    off_end[r] -= (4-ind)*16;
                    
                    tmp_mask = _mm512_cmpgt_epi8_mask(a, zero_);
                    _mm512_mask_storeu_epi8(&x[t], msk, _mm512_sub_epi8(_mm512_mask_blend_epi8(tmp_mask, zero_, a), qe_));
                    d = _mm512_or_si512(d, _mm512_mask_blend_epi8(tmp_mask, zero_, _mm512_set1_epi8(0x08))); // d = a > 0? 1<<3 : 0
                    
                    tmp_mask = _mm512_cmpgt_epi8_mask(b, zero_);
                    _mm512_mask_storeu_epi8(&y[t], msk, _mm512_sub_epi8(_mm512_mask_blend_epi8(tmp_mask, zero_, b), qe_));
                    d = _mm512_or_si512(d, _mm512_mask_blend_epi8(tmp_mask, zero_, _mm512_set1_epi8(0x10))); // d = b > 0? 1<<4 : 0                    
                    tmp_mask = _mm512_cmpgt_epi8_mask(a2, zero_);
                    _mm512_mask_storeu_epi8(&x2[t], msk, _mm512_sub_epi8(_mm512_mask_blend_epi8(tmp_mask, zero_, a2), qe2_));
                    d = _mm512_or_si512(d, _mm512_mask_blend_epi8(tmp_mask, zero_, _mm512_set1_epi8(0x20))); // d = a > 0? 1<<5 : 0                    
                    tmp_mask = _mm512_cmpgt_epi8_mask(b2, zero_);
                    _mm512_mask_storeu_epi8(&y2[t], msk, _mm512_sub_epi8(_mm512_mask_blend_epi8(tmp_mask, zero_, b2), qe2_));
                    d = _mm512_or_si512(d, _mm512_mask_blend_epi8(tmp_mask, zero_, _mm512_set1_epi8(0x40))); // d = b > 0? 1<<6 : 0
                    //_mm512_store_si512(&pr[t], d);
                    _mm512_mask_storeu_epi8(&pr[t], msk, d);
                }
                else
                #endif    
                {
                    tmp_mask = _mm512_cmpgt_epi8_mask(a, zero_);
                    _mm512_store_si512(&x[t],  _mm512_sub_epi8(_mm512_mask_blend_epi8(tmp_mask, zero_, a),  qe_));
                    d = _mm512_or_si512(d, _mm512_mask_blend_epi8(tmp_mask, zero_, _mm512_set1_epi8(0x08))); // d = a > 0? 1<<3 : 0
                    tmp_mask = _mm512_cmpgt_epi8_mask(b, zero_);
                    _mm512_store_si512(&y[t],  _mm512_sub_epi8(_mm512_mask_blend_epi8(tmp_mask, zero_, b),  qe_));
                    d = _mm512_or_si512(d, _mm512_mask_blend_epi8(tmp_mask, zero_, _mm512_set1_epi8(0x10))); // d = b > 0? 1<<4 : 0
                    tmp_mask = _mm512_cmpgt_epi8_mask(a2, zero_);
                    _mm512_store_si512(&x2[t], _mm512_sub_epi8(_mm512_mask_blend_epi8(tmp_mask, zero_, a2), qe2_));
                    d = _mm512_or_si512(d, _mm512_mask_blend_epi8(tmp_mask, zero_, _mm512_set1_epi8(0x20))); // d = a > 0? 1<<5 : 0
                    tmp_mask = _mm512_cmpgt_epi8_mask(b2, zero_);
                    _mm512_store_si512(&y2[t], _mm512_sub_epi8(_mm512_mask_blend_epi8(tmp_mask, zero_, b2), qe2_));
                    d = _mm512_or_si512(d, _mm512_mask_blend_epi8(tmp_mask, zero_, _mm512_set1_epi8(0x40))); // d = b > 0? 1<<6 : 0
                    // int64_t tic = _rdtsc();
                    _mm512_store_si512(&pr[t], d);
                    // int64_t toc = _rdtsc();
                    //fprintf(stderr, "store time: %d %ld\n", t, toc - tic);
                }
            }
*/
	    //continue;
        } else { // gap right-alignment
	//	continue;
            __m512i *pr = p + (size_t)r * n_col_ - st_;
            off[r] = st, off_end[r] = en;
            // off[r] = stn, off_end[r] = enn;
            // fprintf(stderr, "t: %d, st0: %d\n", st_, st0);
            #if TIMER_I
            uint64_t tic = _rdtsc();
            #endif
            for (t = st_; t <= en_; ++t) {
                __m512i d, z, a, b, a2, b2, xt1, x2t1, vt1, ut, tmp;
                __mmask64 tmp_mask;
                
                //__dp_code_block1_pcl;
                __dp_code_block1_pcl_512_mm2_fast;
                //d = _mm512_mask_blend_epi8(_mm512_cmpgt_epi8_mask(z, a),_mm512_set1_epi8(1),  zero_);
                d = _mm512_mask_blend_epi8(_mm512_cmpgt_epi8_mask(z, a), one_,  zero_);
                z = _mm512_max_epi8(z, a);
                // d = _mm512_mask_blend_epi8(_mm512_cmpgt_epi8_mask(z, b), _mm512_set1_epi8(2), d);  // d = z > b?  d : 2
                d = _mm512_mask_blend_epi8(_mm512_cmpgt_epi8_mask(z, b), two_, d);  // d = z > b?  d : 2
                z = _mm512_max_epi8(z, b);
                // d = _mm512_mask_blend_epi8(_mm512_cmpgt_epi8_mask(z, a2), _mm512_set1_epi8(3), d); // d = z > a2? d : 3
                d = _mm512_mask_blend_epi8(_mm512_cmpgt_epi8_mask(z, a2), three_, d); // d = z > a2? d : 3
                z = _mm512_max_epi8(z, a2);
                // d = _mm512_mask_blend_epi8(_mm512_cmpgt_epi8_mask(z, b2), _mm512_set1_epi8(4), d); // d = z > b2? d : 4
                d = _mm512_mask_blend_epi8(_mm512_cmpgt_epi8_mask(z, b2), four_, d); // d = z > b2? d : 4
                z = _mm512_max_epi8(z, b2);                
                z = _mm512_min_epi8(z, sc_mch_);

                // __dp_code_block2_pcl;
                
                #if 1
                __mmask64 msk;
                if (t == en_) {
                    // __mmask64 msk;
                    int ind = (en0 - t*64 + 16)/16;
                    msk = msk_ar2[ind];
                    off_end[r] -= (4-ind)*16;
                    
                    _mm512_mask_storeu_epi8(&u[t], msk, _mm512_sub_epi8(z, vt1));
                    _mm512_mask_storeu_epi8(&v[t], msk, _mm512_sub_epi8(z, ut));
                    tmp = _mm512_sub_epi8(z, q_);                   
                    a = _mm512_sub_epi8(a, tmp);                    
                    b = _mm512_sub_epi8(b, tmp);                    
                    tmp = _mm512_sub_epi8(z, q2_);                  
                    a2= _mm512_sub_epi8(a2, tmp);                   
                    b2= _mm512_sub_epi8(b2, tmp);                    
                }
                else {
                    _mm512_store_si512(&u[t], _mm512_sub_epi8(z, vt1));
                    _mm512_store_si512(&v[t], _mm512_sub_epi8(z, ut));
                    tmp = _mm512_sub_epi8(z, q_);                   
                    a = _mm512_sub_epi8(a, tmp);                    
                    b = _mm512_sub_epi8(b, tmp);                    
                    tmp = _mm512_sub_epi8(z, q2_);                  
                    a2= _mm512_sub_epi8(a2, tmp);                   
                    b2= _mm512_sub_epi8(b2, tmp);
                }
                #endif
                
                #if 1
                if (t == en_) {
                    //__mmask64 msk;
                    //int ind = (en0 - t*64 + 16)/16;
                    //msk = msk_ar2[ind];                                        
                    // off_end[r] -= (4-ind)*16;
                    
                    tmp_mask = _mm512_cmpgt_epi8_mask(zero_, a);
                    _mm512_mask_storeu_epi8(&x[t], msk, _mm512_sub_epi8(_mm512_mask_blend_epi8(tmp_mask, a, zero_), qe_));
                    d = _mm512_or_si512(d, _mm512_mask_blend_epi8(tmp_mask, _mm512_set1_epi8(0x08), zero_)); // d = a > 0? 1<<3 : 0
                    tmp_mask = _mm512_cmpgt_epi8_mask(zero_, b);
                    _mm512_mask_storeu_epi8(&y[t], msk, _mm512_sub_epi8(_mm512_mask_blend_epi8(tmp_mask, b, zero_), qe_));
                    d = _mm512_or_si512(d, _mm512_mask_blend_epi8(tmp_mask, _mm512_set1_epi8(0x10), zero_)); // d = b > 0? 1<<4 : 0
                    tmp_mask = _mm512_cmpgt_epi8_mask(zero_, a2);
                    _mm512_mask_storeu_epi8(&x2[t], msk, _mm512_sub_epi8(_mm512_mask_blend_epi8(tmp_mask, a2, zero_), qe2_));
                    d = _mm512_or_si512(d, _mm512_mask_blend_epi8(tmp_mask, _mm512_set1_epi8(0x20), zero_)); // d = a > 0? 1<<5 : 0
                    tmp_mask = _mm512_cmpgt_epi8_mask(zero_, b2);
                    _mm512_mask_storeu_epi8(&y2[t], msk, _mm512_sub_epi8(_mm512_mask_blend_epi8(tmp_mask, b2, zero_), qe2_));
                    d = _mm512_or_si512(d, _mm512_mask_blend_epi8(tmp_mask, _mm512_set1_epi8(0x40), zero_)); // d = b > 0? 1<<6 : 0
                    // _mm512_store_si512(&pr[t], d);
                    _mm512_mask_storeu_epi8(&pr[t], msk, d);

                }
                else
                #endif    
                {
                    tmp_mask = _mm512_cmpgt_epi8_mask(zero_, a);
                    _mm512_store_si512(&x[t],  _mm512_sub_epi8(_mm512_mask_blend_epi8(tmp_mask, a, zero_),  qe_));
                    // d = _mm512_or_si512(d, _mm512_mask_blend_epi8(tmp_mask, _mm512_set1_epi8(0x08), zero_)); // d = a > 0? 1<<3 : 0
                    d = _mm512_or_si512(d, _mm512_mask_blend_epi8(tmp_mask, s1_, zero_)); // d = a > 0? 1<<3 : 0
                    tmp_mask = _mm512_cmpgt_epi8_mask(zero_, b);
                    _mm512_store_si512(&y[t],  _mm512_sub_epi8(_mm512_mask_blend_epi8(tmp_mask, b, zero_),  qe_));
                    // d = _mm512_or_si512(d, _mm512_mask_blend_epi8(tmp_mask, _mm512_set1_epi8(0x10), zero_)); // d = b > 0? 1<<4 : 0
                    d = _mm512_or_si512(d, _mm512_mask_blend_epi8(tmp_mask, s2_, zero_)); // d = b > 0? 1<<4 : 0
                    tmp_mask = _mm512_cmpgt_epi8_mask(zero_, a2);
                    _mm512_store_si512(&x2[t], _mm512_sub_epi8(_mm512_mask_blend_epi8(tmp_mask, a2, zero_), qe2_));
                    // d = _mm512_or_si512(d, _mm512_mask_blend_epi8(tmp_mask, _mm512_set1_epi8(0x20), zero_)); // d = a > 0? 1<<5 : 0
                    d = _mm512_or_si512(d, _mm512_mask_blend_epi8(tmp_mask, s3_, zero_)); // d = a > 0? 1<<5 : 0
                    tmp_mask = _mm512_cmpgt_epi8_mask(zero_, b2);
                    _mm512_store_si512(&y2[t], _mm512_sub_epi8(_mm512_mask_blend_epi8(tmp_mask, b2, zero_), qe2_));
                    // d = _mm512_or_si512(d, _mm512_mask_blend_epi8(tmp_mask, _mm512_set1_epi8(0x40), zero_)); // d = b > 0? 1<<6 : 0
                    d = _mm512_or_si512(d, _mm512_mask_blend_epi8(tmp_mask, s4_, zero_)); // d = b > 0? 1<<6 : 0
                    _mm512_store_si512(&pr[t], d);
                    
                }
            }
            #if TIMER_I
            uint64_t toc = _rdtsc();
            tprof[B4][0] += toc - tic;
            #endif
        }
//	continue;

        #if TIMER_I
        tic = _rdtsc();
        #endif
        if (!approx_max) { // find the exact max with a 32-bit score array
            int32_t max_H, max_t;
            // compute H[], max_H and max_t
            if (r > 0) {
                int32_t HH[16], tt[16], en1 = st0 + (en0 - st0) / 16 * 16, i; //TODO: saurabh HH[4] -> HH[16]. 4*4?? 
                __m512i max_H_, max_t_;
                max_H = H[en0] = en0 > 0? H[en0-1] + u8[en0] : H[en0] + v8[en0]; // special casing the last element

                max_t = en0;
                max_H_ = _mm512_set1_epi32(max_H);
                max_t_ = _mm512_set1_epi32(max_t);
                for (t = st0; t < en1; t += /*4*/16) { // this implements: H[t]+=v8[t]-qe; if(H[t]>max_H) max_H=H[t],max_t=t;
                    __m512i H1, tmp, t_;
                    __mmask16 tmp_mask;
                    H1 = _mm512_loadu_si512((__m512i*)&H[t]);
                    __m128i t__ = _mm_load_si128((__m128i*) &v8[t]);
                    t_ = _mm512_cvtepi8_epi32(t__);
                    // t_ = _mm512_setr_epi32(v8[t], v8[t+1], v8[t+2], v8[t+3],v8[t+4], v8[t+5], v8[t+6], v8[t+7],v8[t+8], v8[t+9], v8[t+10], v8[t+11],v8[t+12], v8[t+13], v8[t+14], v8[t+15]); //TODO: need to understand well
                    H1 = _mm512_add_epi32(H1, t_);
                    _mm512_storeu_si512((__m512i*)&H[t], H1);
                    // new code by vasim: making it 4 lanes to match accuracy
                    {
                        #if 0
                        t_ = _mm512_set1_epi32(t);
                        __m512i shfH, shft, max1, max2;
                        for (int l=0; l<4; l++) {
                            //t_ = _mm512_set1_epi32(t + l*4);
                            // __m512i shfH = H1, shft = t_;
                            shfH = _mm512_shuffle_i32x4(H1, H1, l);
                            shft = _mm512_shuffle_i32x4(t_, t_, l);
                            tmp_mask = _mm512_cmpgt_epi32_mask(shfH, max_H_);
                            max_H_ = _mm512_mask_blend_epi32(tmp_mask, max_H_, shfH);
                            max_t_ = _mm512_mask_blend_epi32(tmp_mask, max_t_, shft);
                            t_ = _mm512_add_epi32(t_, four32_);
                        }
                        #else                        
                        __m512i shfH, shft, max1, max2;
                        t_ = _mm512_set1_epi32(t);
                        t_ = _mm512_add_epi32(t_, bt32_);
                        shfH = _mm512_shuffle_i32x4(H1, H1, 0x31);
                        shft = _mm512_shuffle_i32x4(t_, t_, 0x31);
                        tmp_mask = _mm512_cmpgt_epi32_mask(shfH, H1);
                        max1 = _mm512_mask_blend_epi32(tmp_mask, H1, shfH);
                        max2 = _mm512_mask_blend_epi32(tmp_mask, t_, shft);
                        shfH = _mm512_shuffle_i32x4(max1, max1, 0x2);
                        shft = _mm512_shuffle_i32x4(max2, max2, 0x2);
                        tmp_mask = _mm512_cmpgt_epi32_mask(shfH, max1);
                        max1 = _mm512_mask_blend_epi32(tmp_mask, max1, shfH);
                        max2 = _mm512_mask_blend_epi32(tmp_mask, max2, shft);
                        tmp_mask = _mm512_cmpgt_epi32_mask(max1, max_H_);
                        max_H_ = _mm512_mask_blend_epi32(tmp_mask, max_H_, max1);
                        max_t_ = _mm512_mask_blend_epi32(tmp_mask, max_t_, max2);
                        #endif
                    }                    
                }
                _mm512_storeu_si512((__m512i*)HH, max_H_);
                _mm512_storeu_si512((__m512i*)tt, max_t_);

                // fprintf(stderr, "max_H: %d, max_t: %d\n", max_H, max_t);
                int rem = (en0 - t) / 4;
                for (int l=0; l<rem; l++) {
                    int bt = t;
                    #pragma unroll
                    for (int j=0; j<4; j++) {
                        H[t] += (int32_t)v8[t];
                        if (H[t] > HH[j]) {
                            HH[j] = H[t];
                            tt[j] = bt;
                        }
                        t++;
                    }
                }
                #pragma unroll
                for (i = 0; i < 4; ++i)
                    if (max_H < HH[i]) max_H = HH[i], max_t = tt[i] + i;
                                                                                
                for (; t < en0; ++t) { // for the rest of values that haven't been computed with SSE     
                    H[t] += (int32_t)v8[t];
                    if (H[t] > max_H) {
                        max_H = H[t], max_t = t;
                    }
                }
            } else H[0] = v8[0] - qe, max_H = H[0], max_t = 0; // special casing r==0

            // update ez
            if (en0 == tlen - 1 && H[en0] > ez->mte) {
                // ez->mte = H[en0], ez->mte_q = r - en;
                ez->mte = H[en0], ez->mte_q = r - en_new;
                // fprintf(stderr, "r: %d, en: %d\n", r, en);
            }
            if (r - st0 == qlen - 1 && H[st0] > ez->mqe) {
                ez->mqe = H[st0], ez->mqe_t = st0;
            }
            
            if (ksw_apply_zdrop(ez, 1, max_H, r, max_t, zdrop, e2)) {
                // fprintf(stderr, "Dropping out, r: %d\n", r);
                break;
            }
            if (r == qlen + tlen - 2 && en0 == tlen - 1) {                
                ez->score = H[tlen - 1];
                // fprintf(stderr, "ez->score: %d\n", ez->score);
            }
        } else { // find approximate max; Z-drop might be inaccurate, too.
            if (r > 0) {
                if (last_H0_t >= st0 && last_H0_t <= en0 && last_H0_t + 1 >= st0 && last_H0_t + 1 <= en0) {
                    int32_t d0 = v8[last_H0_t];
                    int32_t d1 = u8[last_H0_t + 1];
                    if (d0 > d1) H0 += d0;
                    else H0 += d1, ++last_H0_t;
                    // fprintf(stderr, "1. r: %d, d0: %d, d1: %d H0: %d, last_H0_t: %d\n", r, d0, d1, H0, last_H0_t);
                } else if (last_H0_t >= st0 && last_H0_t <= en0) {
                    H0 += v8[last_H0_t];
                    // fprintf(stderr, "2. r: %d, H0: %d, last_H0_t: %d, v8: %d\n", r, H0, last_H0_t, v8[last_H0_t]);
                } else {
                    ++last_H0_t, H0 += u8[last_H0_t];
                    // fprintf(stderr, "3. r: %d, H0: %d, last_H0_t: %d, u8: %d\n", r, H0, last_H0_t, u8[last_H0_t]);
                }
            } else H0 = v8[0] - qe, last_H0_t = 0;
            if ((flag & KSW_EZ_APPROX_DROP) && ksw_apply_zdrop(ez, 1, H0, r, last_H0_t, zdrop, e2)) {
                // fprintf(stderr, "Dropping II\n");
                break;
            }
            if (r == qlen + tlen - 2 && en0 == tlen - 1) {
                ez->score = H0;
                // fprintf(stderr, ">>r: %d, ez->score: %d\n", r, ez->score);
            }
        }
        // last_st = st, last_en = en;
        last_st = st_new, last_en = en_new;

        #if TIMER_I
        toc = _rdtsc();
        tprof[B5][0] += toc - tic;
        #endif
    }
    #if TIMER_O
    uint64_t toc = _rdtsc();
    tprof[OUTER][0] += toc - tic;
    #endif
    
    kfree(km, mem);
     if (!approx_max) kfree(km, H);

    #if 1
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
    #endif

    //kfree_all();
}



void ksw_extd2_avx512_v1(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat,
                   int8_t q, int8_t e, int8_t q2, int8_t e2, int w, int zdrop, int end_bonus, int flag, ksw_extz_t *ez)
{
    __m512i bt32_ = _mm512_setr_epi32(0,0,0,0,4,4,4,4,8,8,8,8,12,12,12,12);
    
    int8_t index[64] __attribute((aligned(64)));

    for (int i=0; i<64; i++) 
        index[i] = i%16 - 1;

    index[0]  = 15;
    index[16] = 31;
    index[32] = 47;
    index[48] = 63;
    
    __m512i shf512a, shf512b, slli512;
    __m512i ind512_slli = _mm512_load_si512((__m512i*) index);
    __mmask8 mska = 0x90;
    __mmask64 mskb = 0x0001000100010000;
    __mmask64 mskc = 0x1;
    __mmask64 mskc_ar[4] = {0x1, 0x10000, 0x100000000, 0x1000000000000};



    #define __dp_code_block1_pcl                                        \
        /*__mmask64 mskc_ = mskc;  */                                       \
        /*if (t == st_)                */                                   \
        __mmask64 mskc_ = (t == st_) ? mskc_ar[(st0 - t*64)/16]:mskc;                    \
        z = _mm512_load_si512(&s[t]);                                   \
        xt1 = _mm512_load_si512(&x[t]);                     /* xt1 <- x[r-1][t..t+15] */ \
        /*  tmp = _mm_srli_si128(xt1, 15);   */                /* tmp <- x[r-1][t+15] */ \
        tmp = _mm512_set1_epi8(((int8_t*)x)[t*64 + 63]); \
	/*  xt1 = _mm_or_si128(_mm_slli_si128(xt1, 1), x1_);*/ /* xt1 <- x[r-1][t-1..t+14] */ \
        shf512a = _mm512_shuffle_epi8(xt1, ind512_slli);                \
        shf512b = _mm512_shuffle_i32x4(shf512a, shf512a, mska);         \
        slli512 = _mm512_mask_blend_epi8(mskb, shf512a, shf512b);       \
        xt1 = _mm512_mask_blend_epi8(mskc_, slli512, x1_);               \
        x1_ = tmp;                                                      \
        vt1 = _mm512_load_si512(&v[t]);                     /* vt1 <- v[r-1][t..t+15] */ \
        /*    tmp = _mm_srli_si128(vt1, 15);             */      /* tmp <- v[r-1][t+15] */ \
        tmp = _mm512_set1_epi8(((int8_t*)v)[t*64 + 63]); \
        /*    vt1 = _mm_or_si128(_mm_slli_si128(vt1, 1), v1_); *//* vt1 <- v[r-1][t-1..t+14] */ \
        shf512a = _mm512_shuffle_epi8(vt1, ind512_slli);                \
        shf512b = _mm512_shuffle_i32x4(shf512a, shf512a, mska);         \
        slli512 = _mm512_mask_blend_epi8(mskb, shf512a, shf512b);       \
        vt1 = _mm512_mask_blend_epi8(mskc_, slli512, v1_);               \
        v1_ = tmp;                                                      \
        a = _mm512_add_epi8(xt1, vt1);                      /* a <- x[r-1][t-1..t+14] + v[r-1][t-1..t+14] */ \
        ut = _mm512_load_si512(&u[t]);                      /* ut <- u[t..t+15] */ \
        b = _mm512_add_epi8(_mm512_load_si512(&y[t]), ut);     /* b <- y[r-1][t..t+15] + u[r-1][t..t+15] */ \
        x2t1= _mm512_load_si512(&x2[t]);                                \
        /*  tmp = _mm_srli_si128(x2t1, 15);*/                           \
        tmp = _mm512_set1_epi8(((int8_t*)x2)[t*64 + 63]); \
        /*  x2t1= _mm_or_si128(_mm_slli_si128(x2t1, 1), x21_); */       \
        shf512a = _mm512_shuffle_epi8(x2t1, ind512_slli);               \
        shf512b = _mm512_shuffle_i32x4(shf512a, shf512a, mska);         \
        slli512 = _mm512_mask_blend_epi8(mskb, shf512a, shf512b);       \
        x2t1 = _mm512_mask_blend_epi8(mskc_, slli512, x21_);             \
        x21_= tmp;                                                      \
        a2= _mm512_add_epi8(x2t1, vt1);                                 \
        b2= _mm512_add_epi8(_mm512_load_si512(&y2[t]), ut);

    
    #define __dp_code_block2_pcl                                        \
        _mm512_store_si512(&u[t], _mm512_sub_epi8(z, vt1));    /* u[r][t..t+15] <- z - v[r-1][t-1..t+14] */ \
        _mm512_store_si512(&v[t], _mm512_sub_epi8(z, ut));     /* v[r][t..t+15] <- z - u[r-1][t..t+15] */ \
        tmp = _mm512_sub_epi8(z, q_);                                   \
        a = _mm512_sub_epi8(a, tmp);                                    \
        b = _mm512_sub_epi8(b, tmp);                                    \
        tmp = _mm512_sub_epi8(z, q2_);                                  \
        a2= _mm512_sub_epi8(a2, tmp);                                   \
        b2= _mm512_sub_epi8(b2, tmp);

    __mmask64 msk_ar2[5] = {0xFFFF, 0xFFFF, 0xFFFFFFFF, 0xFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF};
    
    int r, t, qe = q + e, n_col_, *off = 0, *off_end = 0, tlen_, qlen_, last_st, last_en, wl, wr, max_sc, min_sc, long_thres, long_diff;
    int with_cigar = !(flag&KSW_EZ_SCORE_ONLY), approx_max = !!(flag&KSW_EZ_APPROX_MAX);
    int32_t *H = 0, H0 = 0, last_H0_t = 0;
    uint8_t *qr, *sf, *mem, *mem2 = 0;
    
    __m512i q_, q2_, qe_, qe2_, zero_, sc_mch_, sc_mis_, m1_, sc_N_;
    __m512i *u, *v, *x, *y, *x2, *y2, *s, *p = 0;
    __m512i one_, two_, three_, four_, s1_, s2_, s3_, s4_;
        
    ksw_reset_extz(ez);
    if (m <= 1 || qlen <= 0 || tlen <= 0) return;

    if (q2 + e2 < q + e) t = q, q = q2, q2 = t, t = e, e = e2, e2 = t; // make sure q+e no larger than q2+e2
    s1_   = _mm512_set1_epi8(0x08);
    s2_   = _mm512_set1_epi8(0x10);
    s3_   = _mm512_set1_epi8(0x20);
    s4_   = _mm512_set1_epi8(0x40);
    
    one_   = _mm512_set1_epi8(1);
    two_   = _mm512_set1_epi8(2);
    three_   = _mm512_set1_epi8(3);
    four_   = _mm512_set1_epi8(4);
    
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
    tlen_ = (tlen + 63) / 64;
    n_col_ = qlen < tlen? qlen : tlen;
    n_col_ = ((n_col_ < w + 1? n_col_ : w + 1) + 63) / 64 + 1;
    qlen_ = (qlen + 63) / 64;
    for (t = 1, max_sc = mat[0], min_sc = mat[1]; t < m * m; ++t) {
        max_sc = max_sc > mat[t]? max_sc : mat[t];
        min_sc = min_sc < mat[t]? min_sc : mat[t];
    }
    if (-min_sc > 2 * (q + e)) return; // otherwise, we won't see any mismatches

    long_thres = e != e2? (q2 - q) / (e - e2) - 1 : 0;
    if (q2 + e2 + long_thres * e2 > q + e + long_thres * e)
        ++long_thres;
    long_diff = long_thres * (e - e2) - (q2 - q) - e2;

    mem = (uint8_t*)kcalloc(km, tlen_ * 8 + qlen_ + 1 + 63, 64); 
    u = (__m512i*)(((size_t)mem + 63) >> 6 << 6); // 16-byte aligned 
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
        p = (__m512i*)(((size_t)mem2 + 63) >> 6 << 6);
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
        int st_new = st / 16 * 16, en_new = (en + 16) / 16 * 16 - 1;
        // int st_new = st / 64 * 64, en_new = (en + 64) / 64 * 64 - 1;
        //int stb = st, enb = en;
        st = st / 64 * 64, en = (en + 64) / 64 * 64 - 1;
        //int stn = stb / 16 * 16, enn = (enb + 16) / 16 * 16 - 1;
        // set boundary conditions
        if (st_new > 0) {
            if (st_new - 1 >= last_st && st_new - 1 <= last_en) {
                x1 = x8[st_new - 1], x21 = x28[st_new - 1], v1 = v8[st_new - 1]; // (r-1,s-1) calculated in the last round
            } else {
                x1 = -q - e, x21 = -q2 - e2;
                v1 = -q - e;
            }
        } else {
            x1 = -q - e, x21 = -q2 - e2;
            v1 = r == 0? -q - e : r < long_thres? -e : r == long_thres? long_diff : -e2;
        }
        if (en_new >= r) {
            ((int8_t*)y)[r] = -q - e, ((int8_t*)y2)[r] = -q2 - e2;
            u8[r] = r == 0? -q - e : r < long_thres? -e : r == long_thres? long_diff : -e2;
        }
        
        // loop fission: set scores first 
        if (!(flag & KSW_EZ_GENERIC_SC)) {
            for (t = st0; t <= en0; t += 64) {
                __m512i sq, st, tmp_512;
                __mmask64 tmp, mask;
                sq = _mm512_loadu_si512((__m512i*)&sf[t]);
                st = _mm512_loadu_si512((__m512i*)&qrr[t]);
//              mask = _mm512_or_si512(_mm_cmpeq_epi8(sq, m1_), _mm_cmpeq_epi8(st, m1_));
                mask = (_mm512_cmpeq_epi8_mask(sq, m1_) | _mm512_cmpeq_epi8_mask(st, m1_));
                tmp = _mm512_cmpeq_epi8_mask(sq, st);
                
                tmp_512 = _mm512_mask_blend_epi8(tmp, sc_mis_, sc_mch_);
                tmp_512 = _mm512_mask_blend_epi8(mask, tmp_512, sc_N_);
                if (t + 64 > en0)
                {
                    __mmask64 msk;
                    int ind = (en0 - t + 16)/16;
                    //assert(ind >= 0 && ind < 5);
                    msk = msk_ar2[ind];
                    _mm512_mask_storeu_epi8((__m512i*)((int8_t*)s + t), msk, tmp_512);
                }
                else
                    _mm512_storeu_si512((__m512i*)((int8_t*)s + t), tmp_512);
                
            }
        } else {
            for (t = st0; t <= en0; ++t)
                ((uint8_t*)s)[t] = mat[sf[t] * m + qrr[t]];
        }
        
        // core loop 
        // fprintf(stderr, "- r: %d, x1: %d, x21: %d, v1: %d, en_new: %d, e: %d, q: %d\n",
        //r, x1, x21, v1, en_new, e, q);
        x1_  = _mm512_set1_epi8((uint8_t)x1);
        x21_ = _mm512_set1_epi8((uint8_t)x21);
        v1_  = _mm512_set1_epi8((uint8_t)v1);
        
        //st_ = st / 16, en_ = en / 16;     
        st_ = st / 64, en_ = en / 64;
        //assert(en_ - st_ + 1 <= n_col_);
        if (!with_cigar) { // score only
            for (t = st_; t <= en_; ++t) {
                __m512i z, a, b, a2, b2, xt1, x2t1, vt1, ut, tmp;
                __dp_code_block1_pcl;
                
                z = _mm512_max_epi8(z, a);
                z = _mm512_max_epi8(z, b);
                z = _mm512_max_epi8(z, a2);
                z = _mm512_max_epi8(z, b2);
                z = _mm512_min_epi8(z, sc_mch_);
                // __dp_code_block2_pcl; // save u[] and v[]; update a, b, a2 and b2
                if (t == en_) {
                    __mmask64 msk;
                    int ind = (en0 - t*64 + 16)/16;
                    msk = msk_ar2[ind];
                    // fprintf(stderr, "en0: %d, t: %d, ind: %d, msk: %d\n", en0, t, ind, msk);
                    _mm512_mask_storeu_epi8(&u[t], msk, _mm512_sub_epi8(z, vt1));
                    _mm512_mask_storeu_epi8(&v[t], msk, _mm512_sub_epi8(z, ut));
                    tmp = _mm512_sub_epi8(z, q_);                   
                    a = _mm512_sub_epi8(a, tmp);                    
                    b = _mm512_sub_epi8(b, tmp);                    
                    tmp = _mm512_sub_epi8(z, q2_);                  
                    a2= _mm512_sub_epi8(a2, tmp);                   
                    b2= _mm512_sub_epi8(b2, tmp);
                    
                }
                else {
                    _mm512_store_si512(&u[t], _mm512_sub_epi8(z, vt1));
                    _mm512_store_si512(&v[t], _mm512_sub_epi8(z, ut));
                    tmp = _mm512_sub_epi8(z, q_);                   
                    a = _mm512_sub_epi8(a, tmp);                    
                    b = _mm512_sub_epi8(b, tmp);                    
                    tmp = _mm512_sub_epi8(z, q2_);                  
                    a2= _mm512_sub_epi8(a2, tmp);                   
                    b2= _mm512_sub_epi8(b2, tmp);
                }
                
                if (t == en_) {
                    __mmask64 msk;
                    int ind = (en0 - t*64 + 16)/16;
                    //assert(ind >= 0);
                    msk = msk_ar2[ind];
                    // fprintf(stderr, "en0: %d, t: %d, ind: %d, msk: %d\n", en0, t, ind, msk);
                    _mm512_mask_storeu_epi8(&x[t], msk, _mm512_sub_epi8(_mm512_max_epi8(a,  zero_), qe_));
                    _mm512_mask_storeu_epi8(&y[t], msk,  _mm512_sub_epi8(_mm512_max_epi8(b,  zero_), qe_));
                    _mm512_mask_storeu_epi8(&x2[t], msk,  _mm512_sub_epi8(_mm512_max_epi8(a2, zero_), qe2_));
                    _mm512_mask_storeu_epi8(&y2[t], msk, _mm512_sub_epi8(_mm512_max_epi8(b2, zero_), qe2_));                    
                }
                else
                {
                    _mm512_store_si512(&x[t],  _mm512_sub_epi8(_mm512_max_epi8(a,  zero_), qe_));
                    _mm512_store_si512(&y[t],  _mm512_sub_epi8(_mm512_max_epi8(b,  zero_), qe_));
                    _mm512_store_si512(&x2[t], _mm512_sub_epi8(_mm512_max_epi8(a2, zero_), qe2_));
                    _mm512_store_si512(&y2[t], _mm512_sub_epi8(_mm512_max_epi8(b2, zero_), qe2_));
                }
                // for (int l=0; l<64; l++)
                //    fprintf(stderr, "%d ", ((int8_t*)x)[l]);

            }
            
        } else if (!(flag&KSW_EZ_RIGHT)) { // gap left-alignment
            __m512i *pr = p + (size_t)r * n_col_ - st_;
            off[r] = st, off_end[r] = en;

           
            for (t = st_; t < en_; ++t) {
                __m512i d, z, a, b, a2, b2, xt1, x2t1, vt1, ut, tmp;
                __mmask64 tmp_mask;
                __dp_code_block1_pcl;                

                d = _mm512_mask_blend_epi8(_mm512_cmpgt_epi8_mask(a, z),zero_,  one_);       // d = a  > z? 1 : 0
                z = _mm512_max_epi8(z, a);              
                d = _mm512_mask_blend_epi8(_mm512_cmpgt_epi8_mask(b,  z), d, two_); // d = b  > z? 2 : d
                z = _mm512_max_epi8(z, b);
                d = _mm512_mask_blend_epi8(_mm512_cmpgt_epi8_mask(a2, z), d, three_); // d = a2 > z? 3 : d
                z = _mm512_max_epi8(z, a2);
                d = _mm512_mask_blend_epi8(_mm512_cmpgt_epi8_mask(b2, z), d, four_); // d = b2 > z? 4 : d
                z = _mm512_max_epi8(z, b2);
                z = _mm512_min_epi8(z, sc_mch_);
                // __dp_code_block2_pcl;                
                    _mm512_store_si512(&u[t], _mm512_sub_epi8(z, vt1));
                    _mm512_store_si512(&v[t], _mm512_sub_epi8(z, ut));
                    tmp = _mm512_sub_epi8(z, q_);                   
                    a = _mm512_sub_epi8(a, tmp);                    
                    b = _mm512_sub_epi8(b, tmp);                    
                    tmp = _mm512_sub_epi8(z, q2_);                  
                    a2= _mm512_sub_epi8(a2, tmp);                   
                    b2= _mm512_sub_epi8(b2, tmp);
                    tmp_mask = _mm512_cmpgt_epi8_mask(a, zero_);
                    _mm512_store_si512(&x[t],  _mm512_sub_epi8(_mm512_mask_blend_epi8(tmp_mask, zero_, a),  qe_));
                    d = _mm512_or_si512(d, _mm512_mask_blend_epi8(tmp_mask, zero_, s1_)); // d = a > 0? 1<<3 : 0
                    tmp_mask = _mm512_cmpgt_epi8_mask(b, zero_);
                    _mm512_store_si512(&y[t],  _mm512_sub_epi8(_mm512_mask_blend_epi8(tmp_mask, zero_, b),  qe_));
                    d = _mm512_or_si512(d, _mm512_mask_blend_epi8(tmp_mask, zero_, s2_)); // d = b > 0? 1<<4 : 0
                    tmp_mask = _mm512_cmpgt_epi8_mask(a2, zero_);
                    _mm512_store_si512(&x2[t], _mm512_sub_epi8(_mm512_mask_blend_epi8(tmp_mask, zero_, a2), qe2_));
                    d = _mm512_or_si512(d, _mm512_mask_blend_epi8(tmp_mask, zero_, s3_)); // d = a > 0? 1<<5 : 0
                    tmp_mask = _mm512_cmpgt_epi8_mask(b2, zero_);
                    _mm512_store_si512(&y2[t], _mm512_sub_epi8(_mm512_mask_blend_epi8(tmp_mask, zero_, b2), qe2_));
                    d = _mm512_or_si512(d, _mm512_mask_blend_epi8(tmp_mask, zero_, s4_)); // d = b > 0? 1<<6 : 0             
                    _mm512_store_si512(&pr[t], d);
            }
            {
                __m512i d, z, a, b, a2, b2, xt1, x2t1, vt1, ut, tmp;
                __mmask64 tmp_mask;
                __dp_code_block1_pcl;                

                d = _mm512_mask_blend_epi8(_mm512_cmpgt_epi8_mask(a, z),zero_,  one_);       // d = a  > z? 1 : 0
                z = _mm512_max_epi8(z, a);              
                d = _mm512_mask_blend_epi8(_mm512_cmpgt_epi8_mask(b,  z), d, two_); // d = b  > z? 2 : d
                z = _mm512_max_epi8(z, b);
                d = _mm512_mask_blend_epi8(_mm512_cmpgt_epi8_mask(a2, z), d, three_); // d = a2 > z? 3 : d
                z = _mm512_max_epi8(z, a2);
                d = _mm512_mask_blend_epi8(_mm512_cmpgt_epi8_mask(b2, z), d, four_); // d = b2 > z? 3 : d
                z = _mm512_max_epi8(z, b2);
                z = _mm512_min_epi8(z, sc_mch_);
                // __dp_code_block2_pcl;                
                {
                    __mmask64 msk;
                    int ind = (en0 - t*64 + 16)/16;
                  //  //assert(ind >= 0 && ind < 5);
                    msk = msk_ar2[ind];
                    
                    _mm512_mask_storeu_epi8(&u[t], msk, _mm512_sub_epi8(z, vt1));
                    _mm512_mask_storeu_epi8(&v[t], msk, _mm512_sub_epi8(z, ut));
                    tmp = _mm512_sub_epi8(z, q_);                   
                    a = _mm512_sub_epi8(a, tmp);                    
                    b = _mm512_sub_epi8(b, tmp);                    
                    tmp = _mm512_sub_epi8(z, q2_);                  
                    a2= _mm512_sub_epi8(a2, tmp);                   
                    b2= _mm512_sub_epi8(b2, tmp);
                }
                
                {
                    __mmask64 msk;
                    int ind = (en0 - t*64 + 16)/16;
                    msk = msk_ar2[ind];
                    off_end[r] -= (4-ind)*16;
                    
                    tmp_mask = _mm512_cmpgt_epi8_mask(a, zero_);
                    _mm512_mask_storeu_epi8(&x[t], msk, _mm512_sub_epi8(_mm512_mask_blend_epi8(tmp_mask, zero_, a), qe_));
                    d = _mm512_or_si512(d, _mm512_mask_blend_epi8(tmp_mask, zero_, s1_)); // d = a > 0? 1<<3 : 0
                    
                    tmp_mask = _mm512_cmpgt_epi8_mask(b, zero_);
                    _mm512_mask_storeu_epi8(&y[t], msk, _mm512_sub_epi8(_mm512_mask_blend_epi8(tmp_mask, zero_, b), qe_));
                    d = _mm512_or_si512(d, _mm512_mask_blend_epi8(tmp_mask, zero_, s2_)); // d = b > 0? 1<<4 : 0                    
                    tmp_mask = _mm512_cmpgt_epi8_mask(a2, zero_);
                    _mm512_mask_storeu_epi8(&x2[t], msk, _mm512_sub_epi8(_mm512_mask_blend_epi8(tmp_mask, zero_, a2), qe2_));
                    d = _mm512_or_si512(d, _mm512_mask_blend_epi8(tmp_mask, zero_, s3_)); // d = a > 0? 1<<5 : 0                    
                    tmp_mask = _mm512_cmpgt_epi8_mask(b2, zero_);
                    _mm512_mask_storeu_epi8(&y2[t], msk, _mm512_sub_epi8(_mm512_mask_blend_epi8(tmp_mask, zero_, b2), qe2_));
                    d = _mm512_or_si512(d, _mm512_mask_blend_epi8(tmp_mask, zero_, s4_)); // d = b > 0? 1<<6 : 0
                    //_mm512_store_si512(&pr[t], d);
                    _mm512_mask_storeu_epi8(&pr[t], msk, d);
                }
            }


        } else { // gap right-alignment
            __m512i *pr = p + (size_t)r * n_col_ - st_;
            off[r] = st, off_end[r] = en;
            // off[r] = stn, off_end[r] = enn;
            // fprintf(stderr, "t: %d, st0: %d\n", st_, st0);
            for (t = st_; t < en_; ++t) {
                __m512i d, z, a, b, a2, b2, xt1, x2t1, vt1, ut, tmp;
                __mmask64 tmp_mask;
                
                __dp_code_block1_pcl;
                
                d = _mm512_mask_blend_epi8(_mm512_cmpgt_epi8_mask(z, a), one_,  zero_);
                z = _mm512_max_epi8(z, a);
                // d = _mm512_mask_blend_epi8(_mm512_cmpgt_epi8_mask(z, b), _mm512_set1_epi8(2), d);  
                // d = z > b?  d : 2
                d = _mm512_mask_blend_epi8(_mm512_cmpgt_epi8_mask(z, b), two_, d);  // d = z > b?  d : 2
                z = _mm512_max_epi8(z, b);
                // d = z > a2? d : 3
                d = _mm512_mask_blend_epi8(_mm512_cmpgt_epi8_mask(z, a2), three_, d); // d = z > a2? d : 3
                z = _mm512_max_epi8(z, a2);
                // d = z > b2? d : 4
                d = _mm512_mask_blend_epi8(_mm512_cmpgt_epi8_mask(z, b2), four_, d); // d = z > b2? d : 4
                z = _mm512_max_epi8(z, b2);                
                z = _mm512_min_epi8(z, sc_mch_);

                // __dp_code_block2_pcl;
                
                {
                    _mm512_store_si512(&u[t], _mm512_sub_epi8(z, vt1));
                    _mm512_store_si512(&v[t], _mm512_sub_epi8(z, ut));
                    tmp = _mm512_sub_epi8(z, q_);                   
                    a = _mm512_sub_epi8(a, tmp);                    
                    b = _mm512_sub_epi8(b, tmp);                    
                    tmp = _mm512_sub_epi8(z, q2_);                  
                    a2= _mm512_sub_epi8(a2, tmp);                   
                    b2= _mm512_sub_epi8(b2, tmp);
                }
                
                {
                    tmp_mask = _mm512_cmpgt_epi8_mask(zero_, a);
                    _mm512_store_si512(&x[t],  _mm512_sub_epi8(_mm512_mask_blend_epi8(tmp_mask, a, zero_),  qe_));
                    // d = a > 0? 1<<3 : 0
                    d = _mm512_or_si512(d, _mm512_mask_blend_epi8(tmp_mask, s1_, zero_)); // d = a > 0? 1<<3 : 0
                    tmp_mask = _mm512_cmpgt_epi8_mask(zero_, b);
                    _mm512_store_si512(&y[t],  _mm512_sub_epi8(_mm512_mask_blend_epi8(tmp_mask, b, zero_),  qe_));
                    // d = b > 0? 1<<4 : 0
                    d = _mm512_or_si512(d, _mm512_mask_blend_epi8(tmp_mask, s2_, zero_)); // d = b > 0? 1<<4 : 0
                    tmp_mask = _mm512_cmpgt_epi8_mask(zero_, a2);
                    _mm512_store_si512(&x2[t], _mm512_sub_epi8(_mm512_mask_blend_epi8(tmp_mask, a2, zero_), qe2_));
                    // d = a > 0? 1<<5 : 0
                    d = _mm512_or_si512(d, _mm512_mask_blend_epi8(tmp_mask, s3_, zero_)); // d = a > 0? 1<<5 : 0
                    tmp_mask = _mm512_cmpgt_epi8_mask(zero_, b2);
                    _mm512_store_si512(&y2[t], _mm512_sub_epi8(_mm512_mask_blend_epi8(tmp_mask, b2, zero_), qe2_));
                    // d = b > 0? 1<<6 : 0
                    d = _mm512_or_si512(d, _mm512_mask_blend_epi8(tmp_mask, s4_, zero_)); // d = b > 0? 1<<6 : 0
                    _mm512_store_si512(&pr[t], d);
                    
                }
            }
            //for (t = st_; t <= en_; ++t)// Last iteration unrolled 
	    {
                __m512i d, z, a, b, a2, b2, xt1, x2t1, vt1, ut, tmp;
                __mmask64 tmp_mask;
                
                __dp_code_block1_pcl;
                
                d = _mm512_mask_blend_epi8(_mm512_cmpgt_epi8_mask(z, a), one_,  zero_);
                z = _mm512_max_epi8(z, a);
                // d = z > b?  d : 2
                d = _mm512_mask_blend_epi8(_mm512_cmpgt_epi8_mask(z, b), two_, d);  // d = z > b?  d : 2
                z = _mm512_max_epi8(z, b);
                // d = z > a2? d : 3
                d = _mm512_mask_blend_epi8(_mm512_cmpgt_epi8_mask(z, a2), three_, d); // d = z > a2? d : 3
                z = _mm512_max_epi8(z, a2);
                // d = z > b2? d : 4
                d = _mm512_mask_blend_epi8(_mm512_cmpgt_epi8_mask(z, b2), four_, d); // d = z > b2? d : 4
                z = _mm512_max_epi8(z, b2);                
                z = _mm512_min_epi8(z, sc_mch_);

                // __dp_code_block2_pcl;
                
                __mmask64 msk;
                {
                    // __mmask64 msk;
                    int ind = (en0 - t*64 + 16)/16;
                    msk = msk_ar2[ind];
                    off_end[r] -= (4-ind)*16;
                    
                    _mm512_mask_storeu_epi8(&u[t], msk, _mm512_sub_epi8(z, vt1));
                    _mm512_mask_storeu_epi8(&v[t], msk, _mm512_sub_epi8(z, ut));
                    tmp = _mm512_sub_epi8(z, q_);                   
                    a = _mm512_sub_epi8(a, tmp);                    
                    b = _mm512_sub_epi8(b, tmp);                    
                    tmp = _mm512_sub_epi8(z, q2_);                  
                    a2= _mm512_sub_epi8(a2, tmp);                   
                    b2= _mm512_sub_epi8(b2, tmp);                    
                }
                
                {
                    
                    tmp_mask = _mm512_cmpgt_epi8_mask(zero_, a);
                    _mm512_mask_storeu_epi8(&x[t], msk, _mm512_sub_epi8(_mm512_mask_blend_epi8(tmp_mask, a, zero_), qe_));
                    d = _mm512_or_si512(d, _mm512_mask_blend_epi8(tmp_mask, s1_, zero_)); // d = a > 0? 1<<3 : 0
                    tmp_mask = _mm512_cmpgt_epi8_mask(zero_, b);
                    _mm512_mask_storeu_epi8(&y[t], msk, _mm512_sub_epi8(_mm512_mask_blend_epi8(tmp_mask, b, zero_), qe_));
                    d = _mm512_or_si512(d, _mm512_mask_blend_epi8(tmp_mask, s2_, zero_)); // d = b > 0? 1<<4 : 0
                    tmp_mask = _mm512_cmpgt_epi8_mask(zero_, a2);
                    _mm512_mask_storeu_epi8(&x2[t], msk, _mm512_sub_epi8(_mm512_mask_blend_epi8(tmp_mask, a2, zero_), qe2_));
                    d = _mm512_or_si512(d, _mm512_mask_blend_epi8(tmp_mask, s3_, zero_)); // d = a > 0? 1<<5 : 0
                    tmp_mask = _mm512_cmpgt_epi8_mask(zero_, b2);
                    _mm512_mask_storeu_epi8(&y2[t], msk, _mm512_sub_epi8(_mm512_mask_blend_epi8(tmp_mask, b2, zero_), qe2_));
                    d = _mm512_or_si512(d, _mm512_mask_blend_epi8(tmp_mask, s4_, zero_)); // d = b > 0? 1<<6 : 0
                    // _mm512_store_si512(&pr[t], d);
                    _mm512_mask_storeu_epi8(&pr[t], msk, d);

                }
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
                for (t = st0; t < en1; t += /*4*/16) { // this implements: H[t]+=v8[t]-qe; if(H[t]>max_H) max_H=H[t],max_t=t;
                    __m512i H1, t_;
                    __mmask16 tmp_mask;
                    H1 = _mm512_loadu_si512((__m512i*)&H[t]);
                    __m128i t__ = _mm_load_si128((__m128i*) &v8[t]);
                    t_ = _mm512_cvtepi8_epi32(t__);
                    H1 = _mm512_add_epi32(H1, t_);
                    _mm512_storeu_si512((__m512i*)&H[t], H1);
                    // making it 4 lanes to match accuracy
                   
			__m512i shfH, shft, max1, max2;
			t_ = _mm512_set1_epi32(t);
			t_ = _mm512_add_epi32(t_, bt32_);
			shfH = _mm512_shuffle_i32x4(H1, H1, 0x31);
			shft = _mm512_shuffle_i32x4(t_, t_, 0x31);
			tmp_mask = _mm512_cmpgt_epi32_mask(shfH, H1);
			max1 = _mm512_mask_blend_epi32(tmp_mask, H1, shfH);
			max2 = _mm512_mask_blend_epi32(tmp_mask, t_, shft);
			shfH = _mm512_shuffle_i32x4(max1, max1, 0x2);
			shft = _mm512_shuffle_i32x4(max2, max2, 0x2);
			tmp_mask = _mm512_cmpgt_epi32_mask(shfH, max1);
			max1 = _mm512_mask_blend_epi32(tmp_mask, max1, shfH);
			max2 = _mm512_mask_blend_epi32(tmp_mask, max2, shft);
			tmp_mask = _mm512_cmpgt_epi32_mask(max1, max_H_);
			max_H_ = _mm512_mask_blend_epi32(tmp_mask, max_H_, max1);
			max_t_ = _mm512_mask_blend_epi32(tmp_mask, max_t_, max2);


                }
                _mm512_storeu_si512((__m512i*)HH, max_H_);
                _mm512_storeu_si512((__m512i*)tt, max_t_);

                int rem = (en0 - t) / 4;
                for (int l=0; l<rem; l++) {
                    int bt = t;
                    for (int j=0; j<4; j++) {
                        H[t] += (int32_t)v8[t];
                        if (H[t] > HH[j]) {
                            HH[j] = H[t];
                            tt[j] = bt;
                        }
                        t++;
                    }
                }

                for (i = 0; i < 4; ++i)
                    if (max_H < HH[i]) max_H = HH[i], max_t = tt[i] + i;
                                                                                
                for (; t < en0; ++t) { // for the rest of values that haven't been computed with SSE     
                    H[t] += (int32_t)v8[t];
                    if (H[t] > max_H) {
                        max_H = H[t], max_t = t;
                    }
                }
            } else H[0] = v8[0] - qe, max_H = H[0], max_t = 0; // special casing r==0

            // update ez
            if (en0 == tlen - 1 && H[en0] > ez->mte) {
                ez->mte = H[en0], ez->mte_q = r - en_new;
            }
            if (r - st0 == qlen - 1 && H[st0] > ez->mqe) {
                ez->mqe = H[st0], ez->mqe_t = st0;
            }
            
            if (ksw_apply_zdrop(ez, 1, max_H, r, max_t, zdrop, e2)) {
                break;
            }
            if (r == qlen + tlen - 2 && en0 == tlen - 1) {                
                ez->score = H[tlen - 1];
            }
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
            if ((flag & KSW_EZ_APPROX_DROP) && ksw_apply_zdrop(ez, 1, H0, r, last_H0_t, zdrop, e2)) {
                break;
            }
            if (r == qlen + tlen - 2 && en0 == tlen - 1) {
                ez->score = H0;
            }
        }
        // last_st = st, last_en = en;
        last_st = st_new, last_en = en_new;

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
#undef __dp_code_block1_pcl
#undef __dp_code_block2_pcl
}

#endif

#ifdef __AVX2__
/*
__m256i get_mask_store(__m256i msk, void* addr){
	return _mm256_or_si256(msk , _mm256_loadu_si256((__m256i*)addr));
}
*/
__m256i get_mask_store(__m256i msk, void* addr, __m256i store_data){
	return _mm256_blendv_epi8(_mm256_loadu_si256((__m256i*)addr), store_data, msk);
}

extern void *km1;
extern uint64_t km_size;// = 500000000; // 500 MB
extern int km_top;
extern uint64_t avg;
extern uint64_t alignment_time;
void *kcalloc_(void* km, int count, int size)
{
 
//   assert(km_top == 0);
//   assert(count*size < km_size);
    if(!(km_top == 0 && count*size < km_size)){
	fprintf(stderr, "%d %lld %lld %lld %lld", km_top, count, size, km_size, count*size);
	
	}
    assert(km_top == 0);
    assert(count*size < km_size);
    km_top += count*size + 1024;
    memset(km, 0, count * size);
        
    // printf("km_top: %d\n", km_top);
    return km;
}

void *kmalloc_(void* km, int count) {
    if(km_top + count >= km_size)
        printf("count: %d\n", count);
    assert(km_top + count < km_size);
    void *mem = (void*) ((int8_t*) km + km_top);
    km_top += count + 1024;
    // printf("km_top: %d\n", km_top);
    return mem;
}

void kfree_all() { km_top = 0;}
    
void ksw_extd2_avx2(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat,
                   int8_t q, int8_t e, int8_t q2, int8_t e2, int w, int zdrop, int end_bonus, int flag, ksw_extz_t *ez)
{


    __m256i bt32_ = _mm256_setr_epi32(0,0,0,0,4,4,4,4);//8,8,8,8,12,12,12,12);
    
    int8_t index[32] __attribute((aligned(64)));

    for (int i=0; i<32; i++) 
        index[i] = i%16 - 1;

    index[0]  = 15;
    index[16] = 31;
//    index[32] = 47;
//    index[48] = 63;
    
    __m256i shf256a, shf256b, slli256;
    __m256i ind256_slli = _mm256_load_si256((__m256i*) index);
    //__mmask8 mska = 0x00;//0x90
    //__mmask32 mskb = 0x00010000;//0000 0000 0000 0001 0000 0000 0000 0000
    //__mmask32 mskc = 0x1;
	
    __m256i mskb_v = _mm256_set_epi32(0,0,0,255,0,0,0,0);
    __m256i mskc_v = _mm256_set_epi32(0,0,0,0,0,0,0,255);

    //__mmask32 mskc_ar[2] = {0x1, 0x10000};
    __m256i mskc_ar_v[2];// = {0x1, 0x10000};
    mskc_ar_v[0] = _mm256_set_epi32(0,0,0,0,0,0,0,255); 
    mskc_ar_v[1] = _mm256_set_epi32(0,0,0,255,0,0,0,0); 



    #define __dp_code_block1_pcl                                        \
       /* __mmask32 mskc_ = (t == st_) ? mskc_ar[(st0 - t*32)/16]:mskc;     */               \
        __m256i mskc_ = (t == st_) ? mskc_ar_v[(st0 - t*32)/16]:mskc_v;                    \
        z = _mm256_load_si256(&s[t]);                                   \
        xt1 = _mm256_load_si256(&x[t]);                     /* xt1 <- x[r-1][t..t+15] */ \
        /*  tmp = _mm_srli_si128(xt1, 15);   */                /* tmp <- x[r-1][t+15] */ \
        tmp = _mm256_set1_epi8(((int8_t*)x)[t*32 + 31]); \
	/*  xt1 = _mm_or_si128(_mm_slli_si128(xt1, 1), x1_);*/ /* xt1 <- x[r-1][t-1..t+14] */ \
        shf256a = _mm256_shuffle_epi8(xt1, ind256_slli);                \
        /*shf256b = _mm256_shuffle_i32x4(shf256a, shf256a, 0x00);   */      \
        shf256b = _mm256_permute2x128_si256(shf256a, shf256a, 0);         \
/*        slli256 = _mm256_mask_blend_epi8(mskb, shf256a, shf256b);  */     \
        slli256 = _mm256_blendv_epi8(shf256a, shf256b, mskb_v);       \
        xt1 = _mm256_blendv_epi8(slli256, x1_, mskc_);               \
        x1_ = tmp;                                                      \
        vt1 = _mm256_load_si256(&v[t]);                     /* vt1 <- v[r-1][t..t+15] */ \
        /*    tmp = _mm_srli_si128(vt1, 15);             */      /* tmp <- v[r-1][t+15] */ \
        tmp = _mm256_set1_epi8(((int8_t*)v)[t*32 + 31]); \
        /*    vt1 = _mm_or_si128(_mm_slli_si128(vt1, 1), v1_); *//* vt1 <- v[r-1][t-1..t+14] */ \
        shf256a = _mm256_shuffle_epi8(vt1, ind256_slli);                \
        /*shf256b = _mm256_shuffle_i32x4(shf256a, shf256a, 0x00); */        \
        shf256b =  _mm256_permute2x128_si256(shf256a, shf256a, 0);         \
/*        slli256 = _mm256_mask_blend_epi8(mskb, shf256a, shf256b);    */   \
        slli256 = _mm256_blendv_epi8(shf256a, shf256b, mskb_v);       \
        vt1 = _mm256_blendv_epi8(slli256, v1_, mskc_);               \
        v1_ = tmp;                                                      \
        a = _mm256_add_epi8(xt1, vt1);                      /* a <- x[r-1][t-1..t+14] + v[r-1][t-1..t+14] */ \
        ut = _mm256_load_si256(&u[t]);                      /* ut <- u[t..t+15] */ \
        b = _mm256_add_epi8(_mm256_load_si256(&y[t]), ut);     /* b <- y[r-1][t..t+15] + u[r-1][t..t+15] */ \
        x2t1= _mm256_load_si256(&x2[t]);                                \
        /*  tmp = _mm_srli_si128(x2t1, 15);*/                           \
        tmp = _mm256_set1_epi8(((int8_t*)x2)[t*32 + 31]); \
        /*  x2t1= _mm_or_si128(_mm_slli_si128(x2t1, 1), x21_); */       \
        shf256a = _mm256_shuffle_epi8(x2t1, ind256_slli);               \
        shf256b =  _mm256_permute2x128_si256(shf256a, shf256a, 0);         \
        slli256 = _mm256_blendv_epi8(shf256a, shf256b, mskb_v);       \
        x2t1 = _mm256_blendv_epi8(slli256, x21_, mskc_);             \
        x21_= tmp;                                                      \
        a2= _mm256_add_epi8(x2t1, vt1);                                 \
        b2= _mm256_add_epi8(_mm256_load_si256(&y2[t]), ut);

    
    #define __dp_code_block2_pcl                                        \
        _mm256_storeu_si256(&u[t], _mm256_sub_epi8(z, vt1));    /* u[r][t..t+15] <- z - v[r-1][t-1..t+14] */ \
        _mm256_storeu_si256(&v[t], _mm256_sub_epi8(z, ut));     /* v[r][t..t+15] <- z - u[r-1][t..t+15] */ \
        tmp = _mm256_sub_epi8(z, q_);                                   \
        a = _mm256_sub_epi8(a, tmp);                                    \
        b = _mm256_sub_epi8(b, tmp);                                    \
        tmp = _mm256_sub_epi8(z, q2_);                                  \
        a2= _mm256_sub_epi8(a2, tmp);                                   \
        b2= _mm256_sub_epi8(b2, tmp);

    //__mmask32 msk_ar2[3] = {0xFFFF, 0xFFFF, 0xFFFFFFFF};
    __m256i msk_ar2_v[3];
    msk_ar2_v[0] = _mm256_set_epi32(0,0,0,0,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF); 
    msk_ar2_v[1] = _mm256_set_epi32(0,0,0,0,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF); 
    msk_ar2_v[2] = _mm256_set_epi32(0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF); 
 

    
    int r, t, qe = q + e, n_col_, *off = 0, *off_end = 0, tlen_, qlen_, last_st, last_en, wl, wr, max_sc, min_sc, long_thres, long_diff;
    int with_cigar = !(flag&KSW_EZ_SCORE_ONLY), approx_max = !!(flag&KSW_EZ_APPROX_MAX);
    int32_t *H = 0, H0 = 0, last_H0_t = 0;
    uint8_t *qr, *sf, *mem, *mem2 = 0;
    
    __m256i q_, q2_, qe_, qe2_, zero_, sc_mch_, sc_mis_, m1_, sc_N_;
    __m256i *u, *v, *x, *y, *x2, *y2, *s, *p = 0;
    __m256i one_, two_, three_, four_, s1_, s2_, s3_, s4_;
        
    ksw_reset_extz(ez);
    if (m <= 1 || qlen <= 0 || tlen <= 0) return;

    if (q2 + e2 < q + e) t = q, q = q2, q2 = t, t = e, e = e2, e2 = t; // make sure q+e no larger than q2+e2
    s1_   = _mm256_set1_epi8(0x08);
    s2_   = _mm256_set1_epi8(0x10);
    s3_   = _mm256_set1_epi8(0x20);
    s4_   = _mm256_set1_epi8(0x40);
    
    one_   = _mm256_set1_epi8(1);
    two_   = _mm256_set1_epi8(2);
    three_   = _mm256_set1_epi8(3);
    four_   = _mm256_set1_epi8(4);
    
    zero_   = _mm256_set1_epi8(0);
    q_      = _mm256_set1_epi8(q);
    q2_     = _mm256_set1_epi8(q2);
    qe_     = _mm256_set1_epi8(q + e);
    qe2_    = _mm256_set1_epi8(q2 + e2);
    sc_mch_ = _mm256_set1_epi8(mat[0]);
    sc_mis_ = _mm256_set1_epi8(mat[1]);
    sc_N_   = mat[m*m-1] == 0? _mm256_set1_epi8(-e2) : _mm256_set1_epi8(mat[m*m-1]);
    m1_     = _mm256_set1_epi8(m - 1); // wildcard

    if (w < 0) w = tlen > qlen? tlen : qlen;
    wl = wr = w;
    tlen_ = (tlen + 31) / 32;
    n_col_ = qlen < tlen? qlen : tlen;
    n_col_ = ((n_col_ < w + 1? n_col_ : w + 1) + 31) / 32 + 1;
    qlen_ = (qlen + 31) / 32;
    for (t = 1, max_sc = mat[0], min_sc = mat[1]; t < m * m; ++t) {
        max_sc = max_sc > mat[t]? max_sc : mat[t];
        min_sc = min_sc < mat[t]? min_sc : mat[t];
    }
    if (-min_sc > 2 * (q + e)) return; // otherwise, we won't see any mismatches

    long_thres = e != e2? (q2 - q) / (e - e2) - 1 : 0;
    if (q2 + e2 + long_thres * e2 > q + e + long_thres * e)
        ++long_thres;
    long_diff = long_thres * (e - e2) - (q2 - q) - e2;

    //mem = (uint8_t*)kcalloc_(km1, tlen_ * 8 + qlen_ + 1 + 63, 64); 
    mem = (uint8_t*)kcalloc(km, tlen_ * 8 + qlen_ + 1 + 63, 64); 
    u = (__m256i*)(((size_t)mem + 31) >> 5 << 5); // 16-byte aligned 
    v = u + tlen_, x = v + tlen_, y = x + tlen_, x2 = y + tlen_, y2 = x2 + tlen_;
    s = y2 + tlen_, sf = (uint8_t*)(s + tlen_), qr = sf + tlen_ * 32;
    memset(u,  -q  - e,  tlen_ * 32);
    memset(v,  -q  - e,  tlen_ * 32);
    memset(x,  -q  - e,  tlen_ * 32);
    memset(y,  -q  - e,  tlen_ * 32);
    memset(x2, -q2 - e2, tlen_ * 32);
    memset(y2, -q2 - e2, tlen_ * 32);
    if (!approx_max) {
        //H = (int32_t*)kmalloc_(km1, tlen_ * 32 * 4);
        H = (int32_t*)kmalloc(km, tlen_ * 32 * 4);
        for (t = 0; t < tlen_ * 32; ++t) H[t] = KSW_NEG_INF;
    }
    if (with_cigar) {
        //mem2 = (uint8_t*)kmalloc_(km1, ((size_t)(qlen + tlen - 1) * n_col_ + 1) * 32);
        mem2 = (uint8_t*)kmalloc(km, ((size_t)(qlen + tlen - 1) * n_col_ + 1) * 32);
        p = (__m256i*)(((size_t)mem2 + 31) >> 5 << 5);
        //off = (int*)kmalloc_(km1, (qlen + tlen - 1) * sizeof(int) * 2);
        off = (int*)kmalloc(km, (qlen + tlen - 1) * sizeof(int) * 2);
        off_end = off + qlen + tlen - 1;
    }
//	for(uint64_t itr = 0; itr < ((qlen + tlen - 1) * n_col_ + 1) * 1; itr ++){
//		avg+= mem2[itr];
//	}
#ifdef MANUAL_PROFILING
//	uint64_t align_start = __rdtsc();
#endif
	

    for (t = 0; t < qlen; ++t) qr[t] = query[qlen - 1 - t];
    memcpy(sf, target, tlen);

    
    for (r = 0, last_st = last_en = -1; r < qlen + tlen - 1; ++r) {
        int st = 0, en = tlen - 1, st0, en0, st_, en_;
        int8_t x1, x21, v1;
        uint8_t *qrr = qr + (qlen - 1 - r);
        int8_t *u8 = (int8_t*)u, *v8 = (int8_t*)v, *x8 = (int8_t*)x, *x28 = (int8_t*)x2;
        __m256i x1_, x21_, v1_;
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
        int st_new = st / 16 * 16, en_new = (en + 16) / 16 * 16 - 1;
        // int st_new = st / 64 * 64, en_new = (en + 64) / 64 * 64 - 1;
        //int stb = st, enb = en;
        st = st / 32 * 32, en = (en + 32) / 32 * 32 - 1;
        //int stn = stb / 16 * 16, enn = (enb + 16) / 16 * 16 - 1;
        // set boundary conditions
        if (st_new > 0) {
            if (st_new - 1 >= last_st && st_new - 1 <= last_en) {
                x1 = x8[st_new - 1], x21 = x28[st_new - 1], v1 = v8[st_new - 1]; // (r-1,s-1) calculated in the last round
            } else {
                x1 = -q - e, x21 = -q2 - e2;
                v1 = -q - e;
            }
        } else {
            x1 = -q - e, x21 = -q2 - e2;
            v1 = r == 0? -q - e : r < long_thres? -e : r == long_thres? long_diff : -e2;
        }
        if (en_new >= r) {
            ((int8_t*)y)[r] = -q - e, ((int8_t*)y2)[r] = -q2 - e2;
            u8[r] = r == 0? -q - e : r < long_thres? -e : r == long_thres? long_diff : -e2;
        }
        
        // loop fission: set scores first 
        if (!(flag & KSW_EZ_GENERIC_SC)) {
            for (t = st0; t <= en0; t += 32) {
                __m256i sq, st, tmp_256, mask_256;
                sq = _mm256_loadu_si256((__m256i*)&sf[t]);
                st = _mm256_loadu_si256((__m256i*)&qrr[t]);
//              mask = (_mm256_cmpeq_epi8_mask(sq, m1_) | _mm256_cmpeq_epi8_mask(st, m1_));
                mask_256 = _mm256_or_si256(_mm256_cmpeq_epi8(sq, m1_), _mm256_cmpeq_epi8(st, m1_));
                
		tmp_256 = _mm256_cmpeq_epi8(sq, st);
                
                tmp_256 = _mm256_blendv_epi8(sc_mis_, sc_mch_, tmp_256);
                tmp_256 = _mm256_blendv_epi8(tmp_256, sc_N_, mask_256);
                if (t + 32 > en0)
                {
                    int ind = (en0 - t + 16)/16;
                    //assert(ind >= 0 && ind < 3);
		    __m256i msk_v = msk_ar2_v[ind];
		    //__m256i str = (get_mask_store(msk_v, ((int8_t*)s + t), tmp_256));// msk_ar2_v[ind];
		    //__m256i str =_mm256_and_si256(get_mask_store(msk_v,(int8_t*)s + t), tmp_256);// msk_ar2_v[ind];

			
                    //_mm256_storeu_si256((__m256i*)((int8_t*)s + t), str);
		    
                    _mm256_maskstore_epi32((int32_t *)((int8_t*)s + t), msk_v,  tmp_256);
                }
                else
                    _mm256_storeu_si256((__m256i*)((int8_t*)s + t), tmp_256);
                
            }
        } else {
            for (t = st0; t <= en0; ++t)
                ((uint8_t*)s)[t] = mat[sf[t] * m + qrr[t]];
        }
        
        // core loop 
        // fprintf(stderr, "- r: %d, x1: %d, x21: %d, v1: %d, en_new: %d, e: %d, q: %d\n",
        //r, x1, x21, v1, en_new, e, q);
        x1_  = _mm256_set1_epi8((uint8_t)x1);
        x21_ = _mm256_set1_epi8((uint8_t)x21);
        v1_  = _mm256_set1_epi8((uint8_t)v1);
        
        //st_ = st / 16, en_ = en / 16;     
        st_ = st / 32, en_ = en / 32;
        //assert(en_ - st_ + 1 <= n_col_);
        if (!with_cigar) { // score only
            for (t = st_; t <= en_; ++t) {
                __m256i z, a, b, a2, b2, xt1, x2t1, vt1, ut, tmp;
                __dp_code_block1_pcl;
                
                z = _mm256_max_epi8(z, a);
                z = _mm256_max_epi8(z, b);
                z = _mm256_max_epi8(z, a2);
                z = _mm256_max_epi8(z, b2);
                z = _mm256_min_epi8(z, sc_mch_);
                // __dp_code_block2_pcl; // save u[] and v[]; update a, b, a2 and b2
                if (t == en_) {
                    int ind = (en0 - t*32 + 16)/16;//doubt
                    // fprintf(stderr, "en0: %d, t: %d, ind: %d, msk: %d\n", en0, t, ind, msk);
                    //_mm256_storeu_si256(&u[t], (get_mask_store(msk_ar2_v[ind], &u[t], _mm256_sub_epi8(z, vt1))));
                    _mm256_maskstore_epi32((int32_t *)&u[t], msk_ar2_v[ind],  _mm256_sub_epi8(z, vt1));
                    //_mm256_storeu_si256(&v[t], (get_mask_store(msk_ar2_v[ind], &v[t], _mm256_sub_epi8(z, ut))));
                    _mm256_maskstore_epi32((int32_t *)&v[t], msk_ar2_v[ind],  _mm256_sub_epi8(z, ut));
                    tmp = _mm256_sub_epi8(z, q_);                   
                    a = _mm256_sub_epi8(a, tmp);                    
                    b = _mm256_sub_epi8(b, tmp);                    
                    tmp = _mm256_sub_epi8(z, q2_);                  
                    a2= _mm256_sub_epi8(a2, tmp);                   
                    b2= _mm256_sub_epi8(b2, tmp);
                    
                }
                else {
                    _mm256_storeu_si256(&u[t], _mm256_sub_epi8(z, vt1));
                    _mm256_storeu_si256(&v[t], _mm256_sub_epi8(z, ut));
                    tmp = _mm256_sub_epi8(z, q_);                   
                    a = _mm256_sub_epi8(a, tmp);                    
                    b = _mm256_sub_epi8(b, tmp);                    
                    tmp = _mm256_sub_epi8(z, q2_);                  
                    a2= _mm256_sub_epi8(a2, tmp);                   
                    b2= _mm256_sub_epi8(b2, tmp);
                }
                
                if (t == en_) {
                    //__mmask32 msk;
                    int ind = (en0 - t*32 + 16)/16;//doubt
                    //assert(ind >= 0);
                    //msk = msk_ar2[ind];
                    // fprintf(stderr, "en0: %d, t: %d, ind: %d, msk: %d\n", en0, t, ind, msk);
                    //_mm256_storeu_si256(&x[t], (get_mask_store(msk_ar2_v[ind], &x[t], _mm256_sub_epi8(_mm256_max_epi8(a,  zero_), qe_))));
                    _mm256_maskstore_epi32((int32_t *)&x[t], msk_ar2_v[ind],  _mm256_sub_epi8(_mm256_max_epi8(a,  zero_), qe_));
                    //_mm256_storeu_si256(&y[t], (get_mask_store(msk_ar2_v[ind],&y[t],  _mm256_sub_epi8(_mm256_max_epi8(b,  zero_), qe_))));
                    _mm256_maskstore_epi32((int32_t *)&y[t], msk_ar2_v[ind],  _mm256_sub_epi8(_mm256_max_epi8(b,  zero_), qe_));
                    //_mm256_storeu_si256(&x2[t], (get_mask_store(msk_ar2_v[ind], &x2[t] ,  _mm256_sub_epi8(_mm256_max_epi8(a2, zero_), qe2_))));
                    _mm256_maskstore_epi32((int32_t *)&x2[t], msk_ar2_v[ind],  _mm256_sub_epi8(_mm256_max_epi8(a2,  zero_), qe2_));
                    //_mm256_storeu_si256(&y2[t], (get_mask_store(msk_ar2_v[ind], &y2[t], _mm256_sub_epi8(_mm256_max_epi8(b2, zero_), qe2_))));                    
                    _mm256_maskstore_epi32((int32_t *)&y2[t], msk_ar2_v[ind],  _mm256_sub_epi8(_mm256_max_epi8(b2,  zero_), qe2_));
                }
                else
                {
                    _mm256_storeu_si256(&x[t],  _mm256_sub_epi8(_mm256_max_epi8(a,  zero_), qe_));
                    _mm256_storeu_si256(&y[t],  _mm256_sub_epi8(_mm256_max_epi8(b,  zero_), qe_));
                    _mm256_storeu_si256(&x2[t], _mm256_sub_epi8(_mm256_max_epi8(a2, zero_), qe2_));
                    _mm256_storeu_si256(&y2[t], _mm256_sub_epi8(_mm256_max_epi8(b2, zero_), qe2_));
                }
                // for (int l=0; l<64; l++)
                //    fprintf(stderr, "%d ", ((int8_t*)x)[l]);

            }
            
        } else if (!(flag&KSW_EZ_RIGHT)) { // gap left-alignment
            __m256i *pr = p + (size_t)r * n_col_ - st_;
            off[r] = st, off_end[r] = en;

           
            _mm_prefetch(pr + 4, _MM_HINT_T0); 
            for (t = st_; t < en_; ++t) {
                __m256i d, z, a, b, a2, b2, xt1, x2t1, vt1, ut, tmp;
            //    __mmask32 tmp_mask;
             	__m256i tmp_mask;
		__dp_code_block1_pcl;                

                d = _mm256_blendv_epi8(zero_,  one_, _mm256_cmpgt_epi8(a, z));       // d = a  > z? 1 : 0
                z = _mm256_max_epi8(z, a);              
                d = _mm256_blendv_epi8(d, two_,_mm256_cmpgt_epi8(b,  z)); // d = b  > z? 2 : d
                z = _mm256_max_epi8(z, b);
                d = _mm256_blendv_epi8(d, three_,_mm256_cmpgt_epi8(a2, z)); // d = a2 > z? 3 : d
                z = _mm256_max_epi8(z, a2);
                d = _mm256_blendv_epi8(d, four_,_mm256_cmpgt_epi8(b2, z)); // d = b2 > z? 4 : d
                z = _mm256_max_epi8(z, b2);
                z = _mm256_min_epi8(z, sc_mch_);
                // __dp_code_block2_pcl;                
                    _mm256_storeu_si256(&u[t], _mm256_sub_epi8(z, vt1));
                    _mm256_storeu_si256(&v[t], _mm256_sub_epi8(z, ut));
                    tmp = _mm256_sub_epi8(z, q_);                   
                    a = _mm256_sub_epi8(a, tmp);                    
                    b = _mm256_sub_epi8(b, tmp);                    
                    tmp = _mm256_sub_epi8(z, q2_);                  
                    a2= _mm256_sub_epi8(a2, tmp);                   
                    b2= _mm256_sub_epi8(b2, tmp);
                    tmp_mask = _mm256_cmpgt_epi8(a, zero_);
                    _mm256_storeu_si256(&x[t],  _mm256_sub_epi8(_mm256_blendv_epi8(zero_, a, tmp_mask),  qe_));
                    d = _mm256_or_si256(d, _mm256_blendv_epi8( zero_, s1_, tmp_mask)); // d = a > 0? 1<<3 : 0
                    tmp_mask = _mm256_cmpgt_epi8(b, zero_);
                    _mm256_storeu_si256(&y[t],  _mm256_sub_epi8(_mm256_blendv_epi8( zero_, b, tmp_mask),  qe_));
                    d = _mm256_or_si256(d, _mm256_blendv_epi8(zero_, s2_, tmp_mask)); // d = b > 0? 1<<4 : 0
                    tmp_mask = _mm256_cmpgt_epi8(a2, zero_);
                    _mm256_storeu_si256(&x2[t], _mm256_sub_epi8(_mm256_blendv_epi8(zero_, a2,tmp_mask), qe2_));
                    d = _mm256_or_si256(d, _mm256_blendv_epi8( zero_, s3_,tmp_mask)); // d = a > 0? 1<<5 : 0
                    tmp_mask = _mm256_cmpgt_epi8(b2, zero_);
                    _mm256_storeu_si256(&y2[t], _mm256_sub_epi8(_mm256_blendv_epi8(zero_, b2,tmp_mask), qe2_));
                    d = _mm256_or_si256(d, _mm256_blendv_epi8(zero_, s4_,tmp_mask)); // d = b > 0? 1<<6 : 0             
                    _mm256_storeu_si256(&pr[t], d);
            	    _mm_prefetch(pr + 4, _MM_HINT_T0); 
            }
            {
                __m256i d, z, a, b, a2, b2, xt1, x2t1, vt1, ut, tmp;
                //__mmask32 tmp_mask;
             	__m256i tmp_mask;
                __dp_code_block1_pcl;                

                d = _mm256_blendv_epi8(zero_,one_, _mm256_cmpgt_epi8(a, z));       // d = a  > z? 1 : 0
                z = _mm256_max_epi8(z, a);              
                d = _mm256_blendv_epi8(d, two_, _mm256_cmpgt_epi8(b,  z)); // d = b  > z? 2 : d
                z = _mm256_max_epi8(z, b);
                d = _mm256_blendv_epi8(d, three_, _mm256_cmpgt_epi8(a2, z)); // d = a2 > z? 3 : d
                z = _mm256_max_epi8(z, a2);
                d = _mm256_blendv_epi8(d, four_, _mm256_cmpgt_epi8(b2, z)); // d = b2 > z? 3 : d
                z = _mm256_max_epi8(z, b2);
                z = _mm256_min_epi8(z, sc_mch_);
                // __dp_code_block2_pcl;                
                {
                    //__mmask32 msk;
                    int ind = (en0 - t*32 + 16)/16;//doubt
                  //  //assert(ind >= 0 && ind < 5);
                    //msk = msk_ar2[ind];
                    
                    //_mm256_storeu_si256(&u[t], (get_mask_store(msk_ar2_v[ind], &u[t], _mm256_sub_epi8(z, vt1))));
                    //_mm256_storeu_si256(&v[t], (get_mask_store(msk_ar2_v[ind],&v[t], _mm256_sub_epi8(z, ut))));
                    _mm256_maskstore_epi32((int32_t *)&u[t], msk_ar2_v[ind],  _mm256_sub_epi8(z, vt1));
                    _mm256_maskstore_epi32((int32_t *)&v[t], msk_ar2_v[ind],  _mm256_sub_epi8(z, ut));
                    

		    tmp = _mm256_sub_epi8(z, q_);                   
                    a = _mm256_sub_epi8(a, tmp);                    
                    b = _mm256_sub_epi8(b, tmp);                    
                    tmp = _mm256_sub_epi8(z, q2_);                  
                    a2= _mm256_sub_epi8(a2, tmp);                   
                    b2= _mm256_sub_epi8(b2, tmp);
                }
                
                {
                    //__mmask32 msk;
                    int ind = (en0 - t*32 + 16)/16;//doubt
                    //msk = msk_ar2[ind];
		    __m256i msk_v= msk_ar2_v[ind];
                    off_end[r] -= (2-ind)*16;//doubt
                    
                    tmp_mask = _mm256_cmpgt_epi8(a, zero_);
                    //_mm256_storeu_si256(&x[t], (get_mask_store(msk_v, &x[t], _mm256_sub_epi8(_mm256_blendv_epi8(zero_, a, tmp_mask), qe_))));
                    _mm256_maskstore_epi32((int32_t *)&x[t], msk_v,  _mm256_sub_epi8(_mm256_blendv_epi8(zero_, a, tmp_mask), qe_));
                    d = _mm256_or_si256(d, _mm256_blendv_epi8(zero_, s1_, tmp_mask)); // d = a > 0? 1<<3 : 0
                    
                    tmp_mask = _mm256_cmpgt_epi8(b, zero_);
                    //_mm256_storeu_si256(&y[t], (get_mask_store(msk_v,&y[t] , _mm256_sub_epi8(_mm256_blendv_epi8(zero_, b, tmp_mask), qe_))));
                    _mm256_maskstore_epi32((int32_t *)&y[t], msk_v,  _mm256_sub_epi8(_mm256_blendv_epi8(zero_, b, tmp_mask), qe_));
                    d = _mm256_or_si256(d, _mm256_blendv_epi8(zero_, s2_, tmp_mask)); // d = b > 0? 1<<4 : 0                    
                    tmp_mask = _mm256_cmpgt_epi8(a2, zero_);
                    //_mm256_storeu_si256(&x2[t], (get_mask_store(msk_v,&x2[t], _mm256_sub_epi8(_mm256_blendv_epi8(zero_, a2, tmp_mask), qe2_))));
                    _mm256_maskstore_epi32((int32_t *)&x2[t], msk_v,  _mm256_sub_epi8(_mm256_blendv_epi8(zero_, a2, tmp_mask), qe2_));
                    d = _mm256_or_si256(d, _mm256_blendv_epi8(zero_, s3_, tmp_mask)); // d = a > 0? 1<<5 : 0                    
                    tmp_mask = _mm256_cmpgt_epi8(b2, zero_);
                    //_mm256_storeu_si256(&y2[t], (get_mask_store(msk_v,&y2[t], _mm256_sub_epi8(_mm256_blendv_epi8(zero_, b2, tmp_mask), qe2_))));
                    _mm256_maskstore_epi32((int32_t *)&y2[t], msk_v,  _mm256_sub_epi8(_mm256_blendv_epi8(zero_, b2, tmp_mask), qe2_));
                    d = _mm256_or_si256(d, _mm256_blendv_epi8(zero_, s4_, tmp_mask)); // d = b > 0? 1<<6 : 0
                    //_mm256_storeu_si256(&pr[t], (get_mask_store(msk_v, &pr[t], d)));
                    _mm256_maskstore_epi32((int32_t *)&pr[t], msk_v, d);
                    //_mm256_mask_storeu_epi8(&pr[t], msk, d);
                }
            }


        } else { // gap right-alignment
            __m256i *pr = p + (size_t)r * n_col_ - st_;
            off[r] = st, off_end[r] = en;
            // off[r] = stn, off_end[r] = enn;
            // fprintf(stderr, "t: %d, st0: %d\n", st_, st0);
            for (t = st_; t < en_; ++t) {
                __m256i d, z, a, b, a2, b2, xt1, x2t1, vt1, ut, tmp;
                //__mmask32 tmp_mask;
                __m256i tmp_mask;
                
                __dp_code_block1_pcl;
                
                d = _mm256_blendv_epi8(one_,  zero_, _mm256_cmpgt_epi8(z, a));
                z = _mm256_max_epi8(z, a);
                // d = _mm256_blendv_epi8(_mm256_cmpgt_epi8(z, b), _mm256_set1_epi8(2), d);  
                // d = z > b?  d : 2
                d = _mm256_blendv_epi8(two_, d,_mm256_cmpgt_epi8(z, b) );  // d = z > b?  d : 2
                z = _mm256_max_epi8(z, b);
                // d = z > a2? d : 3
                d = _mm256_blendv_epi8(three_, d, _mm256_cmpgt_epi8(z, a2)); // d = z > a2? d : 3
                z = _mm256_max_epi8(z, a2);
                // d = z > b2? d : 4
                d = _mm256_blendv_epi8(four_, d, _mm256_cmpgt_epi8(z, b2)); // d = z > b2? d : 4
                z = _mm256_max_epi8(z, b2);                
                z = _mm256_min_epi8(z, sc_mch_);

                // __dp_code_block2_pcl;
                
                //__mmask32 msk;
                {
                    _mm256_storeu_si256(&u[t], _mm256_sub_epi8(z, vt1));
                    _mm256_storeu_si256(&v[t], _mm256_sub_epi8(z, ut));
                    tmp = _mm256_sub_epi8(z, q_);                   
                    a = _mm256_sub_epi8(a, tmp);                    
                    b = _mm256_sub_epi8(b, tmp);                    
                    tmp = _mm256_sub_epi8(z, q2_);                  
                    a2= _mm256_sub_epi8(a2, tmp);                   
                    b2= _mm256_sub_epi8(b2, tmp);
                }
                
                {
                    tmp_mask = _mm256_cmpgt_epi8(zero_, a);
                    _mm256_storeu_si256(&x[t],  _mm256_sub_epi8(_mm256_blendv_epi8(a, zero_,tmp_mask),  qe_));
                    // d = a > 0? 1<<3 : 0
                    d = _mm256_or_si256(d, _mm256_blendv_epi8(s1_, zero_,tmp_mask)); // d = a > 0? 1<<3 : 0
                    tmp_mask = _mm256_cmpgt_epi8(zero_, b);
                    _mm256_storeu_si256(&y[t],  _mm256_sub_epi8(_mm256_blendv_epi8(b, zero_,tmp_mask),  qe_));
                    // d = b > 0? 1<<4 : 0
                    d = _mm256_or_si256(d, _mm256_blendv_epi8(s2_, zero_,tmp_mask)); // d = b > 0? 1<<4 : 0
                    tmp_mask = _mm256_cmpgt_epi8(zero_, a2);
                    _mm256_storeu_si256(&x2[t], _mm256_sub_epi8(_mm256_blendv_epi8(a2, zero_,tmp_mask), qe2_));
                    // d = a > 0? 1<<5 : 0
                    d = _mm256_or_si256(d, _mm256_blendv_epi8(s3_, zero_,tmp_mask)); // d = a > 0? 1<<5 : 0
                    tmp_mask = _mm256_cmpgt_epi8(zero_, b2);
                    _mm256_storeu_si256(&y2[t], _mm256_sub_epi8(_mm256_blendv_epi8(b2, zero_,tmp_mask), qe2_));
                    // d = b > 0? 1<<6 : 0
                    d = _mm256_or_si256(d, _mm256_blendv_epi8(s4_, zero_,tmp_mask)); // d = b > 0? 1<<6 : 0
                    _mm256_storeu_si256(&pr[t], d);
                    
                }
            }
            //for (t = st_; t <= en_; ++t)// Last iteration unrolled 
	    {
                __m256i d, z, a, b, a2, b2, xt1, x2t1, vt1, ut, tmp;
                //__mmask32 tmp_mask;
                __m256i tmp_mask;
                
                __dp_code_block1_pcl;
                
                d = _mm256_blendv_epi8(one_,  zero_,_mm256_cmpgt_epi8(z, a) );
                z = _mm256_max_epi8(z, a);
                // d = z > b?  d : 2
                d = _mm256_blendv_epi8( two_, d, _mm256_cmpgt_epi8(z, b));  // d = z > b?  d : 2
                z = _mm256_max_epi8(z, b);
                // d = z > a2? d : 3
                d = _mm256_blendv_epi8( three_, d, _mm256_cmpgt_epi8(z, a2)); // d = z > a2? d : 3
                z = _mm256_max_epi8(z, a2);
                // d = z > b2? d : 4
                d = _mm256_blendv_epi8(four_, d, _mm256_cmpgt_epi8(z, b2)); // d = z > b2? d : 4
                z = _mm256_max_epi8(z, b2);                
                z = _mm256_min_epi8(z, sc_mch_);

                // __dp_code_block2_pcl;
                
                //__mmask32 msk;
                
                    // __mmask64 msk;
                    int ind = (en0 - t*32 + 16)/16;//doubt
                    //msk = msk_ar2[ind];
                    off_end[r] -= (2-ind)*16;//doubt
                    
                    //_mm256_storeu_si256(&u[t], (get_mask_store(msk_ar2_v[ind], &u[t], _mm256_sub_epi8(z, vt1))));
                    //_mm256_storeu_si256(&v[t], (get_mask_store(msk_ar2_v[ind], &v[t], _mm256_sub_epi8(z, ut))));
                    _mm256_maskstore_epi32((int32_t *)&u[t], msk_ar2_v[ind],  _mm256_sub_epi8(z, vt1));
                    _mm256_maskstore_epi32((int32_t *)&v[t], msk_ar2_v[ind],  _mm256_sub_epi8(z, ut));
                    

		    tmp = _mm256_sub_epi8(z, q_);                   
                    a = _mm256_sub_epi8(a, tmp);                    
                    b = _mm256_sub_epi8(b, tmp);                    
                    tmp = _mm256_sub_epi8(z, q2_);                  
                    a2= _mm256_sub_epi8(a2, tmp);                   
                    b2= _mm256_sub_epi8(b2, tmp);                    
                
                
                
                    
                    tmp_mask = _mm256_cmpgt_epi8(zero_, a);
                    //_mm256_storeu_si256(&x[t], (get_mask_store(msk_ar2_v[ind],&x[t], _mm256_sub_epi8(_mm256_blendv_epi8(a, zero_, tmp_mask), qe_))));
                    _mm256_maskstore_epi32((int32_t *)&x[t], msk_ar2_v[ind],  _mm256_sub_epi8(_mm256_blendv_epi8(a, zero_,tmp_mask), qe_));
                    d = _mm256_or_si256(d, _mm256_blendv_epi8( s1_, zero_, tmp_mask)); // d = a > 0? 1<<3 : 0
                    tmp_mask = _mm256_cmpgt_epi8(zero_, b);
                    //_mm256_storeu_si256(&y[t], (get_mask_store(msk_ar2_v[ind], &y[t], _mm256_sub_epi8(_mm256_blendv_epi8(b, zero_,tmp_mask), qe_))));
                    _mm256_maskstore_epi32((int32_t *)&y[t], msk_ar2_v[ind],  _mm256_sub_epi8(_mm256_blendv_epi8(b, zero_,tmp_mask), qe_));
                    d = _mm256_or_si256(d, _mm256_blendv_epi8(s2_, zero_,tmp_mask)); // d = b > 0? 1<<4 : 0
                    tmp_mask = _mm256_cmpgt_epi8(zero_, a2);
                    //_mm256_storeu_si256(&x2[t], (get_mask_store(msk_ar2_v[ind], &x2[t], _mm256_sub_epi8(_mm256_blendv_epi8(a2, zero_,tmp_mask), qe2_))));
                    _mm256_maskstore_epi32((int32_t *)&x2[t], msk_ar2_v[ind],  _mm256_sub_epi8(_mm256_blendv_epi8(a2, zero_,tmp_mask), qe2_));
                    d = _mm256_or_si256(d, _mm256_blendv_epi8(s3_, zero_,tmp_mask)); // d = a > 0? 1<<5 : 0
                    tmp_mask = _mm256_cmpgt_epi8(zero_, b2);
                    //_mm256_storeu_si256(&y2[t], (get_mask_store(msk_ar2_v[ind], &y2[t], _mm256_sub_epi8(_mm256_blendv_epi8(b2, zero_,tmp_mask), qe2_))));
                    _mm256_maskstore_epi32((int32_t *)&y2[t], msk_ar2_v[ind],  _mm256_sub_epi8(_mm256_blendv_epi8(b2, zero_,tmp_mask), qe2_));
              
	      		d = _mm256_or_si256(d, _mm256_blendv_epi8(s4_, zero_, tmp_mask)); // d = b > 0? 1<<6 : 0
                    // _mm256_storeu_si256(&pr[t], d);
//                    _mm256_storeu_si256(&pr[t], (get_mask_store(msk_ar2_v[ind], &pr[t], d)));
                    _mm256_maskstore_epi32((int32_t *)&pr[t], msk_ar2_v[ind], d);

                
            }
        }

        if (!approx_max) { // find the exact max with a 32-bit score array
            int32_t max_H, max_t;
            // compute H[], max_H and max_t
            if (r > 0) {
                int32_t HH[8], tt[8], en1 = st0 + (en0 - st0) / 8 * 8, i; //doubt
                __m256i max_H_, max_t_;
                max_H = H[en0] = en0 > 0? H[en0-1] + u8[en0] : H[en0] + v8[en0]; // special casing the last element

                max_t = en0;
                max_H_ = _mm256_set1_epi32(max_H);
                max_t_ = _mm256_set1_epi32(max_t);
#if 1
         //      for (t = st0; t < en1; t +=1) { // this implements: H[t]+=v8[t]-qe; if(H[t]>max_H) max_H=H[t],max_t=t;
	//		H[t] += v8[t];
	//	}
                for (t = st0; t < en1; t += /*4*/8) { // this implements: H[t]+=v8[t]-qe; if(H[t]>max_H) max_H=H[t],max_t=t;
                    __m256i H1, t_;
                    H1 = _mm256_loadu_si256((__m256i*)&H[t]);
                    __m128i t__ = _mm_loadu_si128((__m128i*) &v8[t]);
                    t_ = _mm256_cvtepi8_epi32(t__);
               //     t_ = _mm256_setr_epi32(v8[t], v8[t+1], v8[t+2], v8[t+3], v8[t+4], v8[t+5], v8[t+6], v8[t+7]);
                    H1 = _mm256_add_epi32(H1, t_);
                    _mm256_storeu_si256((__m256i*)&H[t], H1);
                    // making it 4 lanes to match accuracy
                   
			__m256i shfH, shft, max1, max2;
                        __m256i tmp_mask_v;
			t_ = _mm256_set1_epi32(t);
			t_ = _mm256_add_epi32(t_, bt32_);
			//shfH = _mm256_shuffle_i32x4(H1, H1, 0x1);//doubt
			//shft = _mm256_shuffle_i32x4(t_, t_, 0x1);//doubt
			shfH =  _mm256_permute2x128_si256(H1, H1, 1);
			shft =  _mm256_permute2x128_si256(t_, t_, 1);
			tmp_mask_v = _mm256_cmpgt_epi32(shfH, H1);
			//max1 = _mm256_or_si256(_mm256_and_si256(tmp_mask_v, shfH) , _mm256_andnot_si256(tmp_mask_v, H1));
			max1 = _mm256_blendv_epi8(H1, shfH, tmp_mask_v);
			//max2 = _mm256_or_si256(_mm256_and_si256(tmp_mask_v, shft) , _mm256_andnot_si256(tmp_mask_v, t_));
			max2 = _mm256_blendv_epi8(t_, shft, tmp_mask_v);

/*			//--shfH = _mm256_shuffle_i32x4(max1, max1, 0x2);//doubt
			//--shft = _mm256_shuffle_i32x4(max2, max2, 0x2);//doubt
			//--tmp_mask = _mm256_cmpgt_epi32_mask(shfH, max1);
			//-max1 = _mm256_mask_blend_epi32(tmp_mask, max1, shfH);
			//--max2 = _mm256_mask_blend_epi32(tmp_mask, max2, shft);

*/
//			tmp_mask = _mm256_cmpgt_epi32_mask(max1, max_H_);
//			max_H_ = _mm256_mask_blend_epi32(tmp_mask, max_H_, max1);
//			max_t_ = _mm256_mask_blend_epi32(tmp_mask, max_t_, max2);
			tmp_mask_v = _mm256_cmpgt_epi32(max1, max_H_);
			//max_H_ = _mm256_or_si256(_mm256_and_si256(tmp_mask_v, max1) , _mm256_andnot_si256(tmp_mask_v, max_H_));
			//max_t_ = _mm256_or_si256(_mm256_and_si256(tmp_mask_v, max2) , _mm256_andnot_si256(tmp_mask_v, max_t_));
			max_H_ = _mm256_blendv_epi8(max_H_, max1, tmp_mask_v);
			max_t_ = _mm256_blendv_epi8(max_t_, max2, tmp_mask_v);


                }
                _mm256_storeu_si256((__m256i*)HH, max_H_);
                _mm256_storeu_si256((__m256i*)tt, max_t_);

                int rem = (en0 - t) / 4;
                for (int l=0; l<rem; l++) {
                    int bt = t;
                    for (int j=0; j<4; j++) {
                        H[t] += (int32_t)v8[t];
                        if (H[t] > HH[j]) {
                            HH[j] = H[t];
                            tt[j] = bt;
                        }
                        t++;
                    }
                }

                for (i = 0; i < 4; ++i)
                    if (max_H < HH[i]) max_H = HH[i], max_t = tt[i] + i;
                for (; t < en0; ++t) { // for the rest of values that haven't been computed with SSE     
                    H[t] += (int32_t)v8[t];
                    if (H[t] > max_H) {
                        max_H = H[t], max_t = t;
                    }
                }
	#else
                for (t = st0 ; t < en0; ++t) { // for the rest of values that haven't been computed with SSE     
                    int32_t tmp = H[t] + v8[t];
		    H[t] = tmp;
                    if (tmp > max_H) {
                        max_H = tmp, max_t = t;
                    }
                }
       #endif                                                                         
            } else H[0] = v8[0] - qe, max_H = H[0], max_t = 0; // special casing r==0

            // update ez
            if (en0 == tlen - 1 && H[en0] > ez->mte) {
                ez->mte = H[en0], ez->mte_q = r - en_new;
            }
            if (r - st0 == qlen - 1 && H[st0] > ez->mqe) {
                ez->mqe = H[st0], ez->mqe_t = st0;
            }
            
            if (ksw_apply_zdrop(ez, 1, max_H, r, max_t, zdrop, e2)) {
                break;
            }
            if (r == qlen + tlen - 2 && en0 == tlen - 1) {                
                ez->score = H[tlen - 1];
            }
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
            if ((flag & KSW_EZ_APPROX_DROP) && ksw_apply_zdrop(ez, 1, H0, r, last_H0_t, zdrop, e2)) {
                break;
            }
            if (r == qlen + tlen - 2 && en0 == tlen - 1) {
                ez->score = H0;
            }
        }
        // last_st = st, last_en = en;
        last_st = st_new, last_en = en_new;

    }
#ifdef MANUAL_PROFILING
//	alignment_time += (__rdtsc() - align_start);
#endif
    
    kfree(km, mem);
    if (!approx_max) kfree(km, H);
    if (with_cigar) { // backtrack
        int rev_cigar = !!(flag & KSW_EZ_REV_CIGAR);
        if (!ez->zdropped && !(flag&KSW_EZ_EXTZ_ONLY)) {
            ksw_backtrack(km, 1, rev_cigar, 0, (uint8_t*)p, off, off_end, n_col_*32, tlen-1, qlen-1, &ez->m_cigar, &ez->n_cigar, &ez->cigar);
        } else if (!ez->zdropped && (flag&KSW_EZ_EXTZ_ONLY) && ez->mqe + end_bonus > (int)ez->max) {
            ez->reach_end = 1;
            ksw_backtrack(km, 1, rev_cigar, 0, (uint8_t*)p, off, off_end, n_col_*32, ez->mqe_t, qlen-1, &ez->m_cigar, &ez->n_cigar, &ez->cigar);
        } else if (ez->max_t >= 0 && ez->max_q >= 0) {
            ksw_backtrack(km, 1, rev_cigar, 0, (uint8_t*)p, off, off_end, n_col_*32, ez->max_t, ez->max_q, &ez->m_cigar, &ez->n_cigar, &ez->cigar);
        }
        kfree(km, mem2); kfree(km, off);
    }
//	kfree_all();
#undef __dp_code_block1_pcl
#undef __dp_code_block2_pcl
}

#endif

