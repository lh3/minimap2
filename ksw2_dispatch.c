#ifdef KSW_CPU_DISPATCH
#include <stdlib.h>
#include "ksw2.h"

#define SIMD_SSE     0x1
#define SIMD_SSE2    0x2
#define SIMD_SSE3    0x4
#define SIMD_SSE4_1  0x8
#define SIMD_SSE4_2  0x10
#define SIMD_AVX     0x20
#define SIMD_AVX2    0x40
#define SIMD_AVX512F 0x80

unsigned x86_simd(void)
{
	unsigned eax, ebx, ecx, edx, flag = 0;
#ifdef _MSC_VER
	int cpuid[4];
	__cpuid(cpuid, 1);
	eax = cpuid[0], ebx = cpuid[1], ecx = cpuid[2], edx = cpuid[3];
#else
	asm volatile("cpuid" : "=a" (eax), "=b" (ebx), "=c" (ecx), "=d" (edx) : "a" (1));
#endif
	if (edx>>25&1) flag |= SIMD_SSE;
	if (edx>>26&1) flag |= SIMD_SSE2;
	if (ecx>>0 &1) flag |= SIMD_SSE3;
	if (ecx>>19&1) flag |= SIMD_SSE4_1;
	if (ecx>>20&1) flag |= SIMD_SSE4_2;
	if (ecx>>28&1) flag |= SIMD_AVX;
	if (ebx>>5 &1) flag |= SIMD_AVX2;
	if (ebx>>16&1) flag |= SIMD_AVX512F;
	return flag;
}

void ksw_extz2_sse(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat, int8_t q, int8_t e, int w, int zdrop, int flag, ksw_extz_t *ez)
{
	extern void ksw_extz2_sse2(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat, int8_t q, int8_t e, int w, int zdrop, int flag, ksw_extz_t *ez);
	extern void ksw_extz2_sse41(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat, int8_t q, int8_t e, int w, int zdrop, int flag, ksw_extz_t *ez);
	unsigned simd;
	simd = x86_simd();
	if (simd & SIMD_SSE4_1)
		ksw_extz2_sse41(km, qlen, query, tlen, target, m, mat, q, e, w, zdrop, flag, ez);
	else if (simd & SIMD_SSE2)
		ksw_extz2_sse2(km, qlen, query, tlen, target, m, mat, q, e, w, zdrop, flag, ez);
	else abort();
}

void ksw_extd2_sse(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat,
				   int8_t q, int8_t e, int8_t q2, int8_t e2, int w, int zdrop, int flag, ksw_extz_t *ez)
{
	extern void ksw_extd2_sse2(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat,
				   int8_t q, int8_t e, int8_t q2, int8_t e2, int w, int zdrop, int flag, ksw_extz_t *ez);
	extern void ksw_extd2_sse41(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat,
				   int8_t q, int8_t e, int8_t q2, int8_t e2, int w, int zdrop, int flag, ksw_extz_t *ez);
	unsigned simd;
	simd = x86_simd();
	if (simd & SIMD_SSE4_1)
		ksw_extd2_sse41(km, qlen, query, tlen, target, m, mat, q, e, q2, e2, w, zdrop, flag, ez);
	else if (simd & SIMD_SSE2)
		ksw_extd2_sse2(km, qlen, query, tlen, target, m, mat, q, e, q2, e2, w, zdrop, flag, ez);
	else abort();
}

void ksw_exts2_sse(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat,
				   int8_t q, int8_t e, int8_t q2, int8_t noncan, int zdrop, int flag, ksw_extz_t *ez)
{
	extern void ksw_exts2_sse2(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat,
				   int8_t q, int8_t e, int8_t q2, int8_t noncan, int zdrop, int flag, ksw_extz_t *ez);
	extern void ksw_exts2_sse41(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat,
				   int8_t q, int8_t e, int8_t q2, int8_t noncan, int zdrop, int flag, ksw_extz_t *ez);
	unsigned simd;
	simd = x86_simd();
	if (simd & SIMD_SSE4_1)
		ksw_exts2_sse41(km, qlen, query, tlen, target, m, mat, q, e, q2, noncan, zdrop, flag, ez);
	else if (simd & SIMD_SSE2)
		ksw_exts2_sse2(km, qlen, query, tlen, target, m, mat, q, e, q2, noncan, zdrop, flag, ez);
	else abort();
}
#endif
