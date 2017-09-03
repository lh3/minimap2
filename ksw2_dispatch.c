#ifdef KSW_CPU_DISPATCH
#include <stdlib.h>
#include "ksw2.h"

void ksw_extz2_sse(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat, int8_t q, int8_t e, int w, int zdrop, int flag, ksw_extz_t *ez)
{
	extern void ksw_extz2_sse2(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat, int8_t q, int8_t e, int w, int zdrop, int flag, ksw_extz_t *ez);
	extern void ksw_extz2_sse41(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat, int8_t q, int8_t e, int w, int zdrop, int flag, ksw_extz_t *ez);
	if (__builtin_cpu_supports("sse4.1"))
		ksw_extz2_sse41(km, qlen, query, tlen, target, m, mat, q, e, w, zdrop, flag, ez);
	else if (__builtin_cpu_supports("sse2"))
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
	if (__builtin_cpu_supports("sse4.1"))
		ksw_extd2_sse41(km, qlen, query, tlen, target, m, mat, q, e, q2, e2, w, zdrop, flag, ez);
	else if (__builtin_cpu_supports("sse2"))
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
	if (__builtin_cpu_supports("sse4.1"))
		ksw_exts2_sse41(km, qlen, query, tlen, target, m, mat, q, e, q2, noncan, zdrop, flag, ez);
	else if (__builtin_cpu_supports("sse2"))
		ksw_exts2_sse2(km, qlen, query, tlen, target, m, mat, q, e, q2, noncan, zdrop, flag, ez);
	else abort();
}
#endif
