#ifndef MMPRIV2_H
#define MMPRIV2_H

#include "minimap.h"

#ifdef __cplusplus
extern "C" {
#endif

double cputime(void);
double realtime(void);

void radix_sort_128x(mm128_t *beg, mm128_t *end);
void radix_sort_64(uint64_t *beg, uint64_t *end);
uint32_t ks_ksmall_uint32_t(size_t n, uint32_t arr[], size_t kk);

int mm_chain_dp(int match_len, int max_dist, int bw, int max_skip, int min_sc, int n, mm128_t *a, uint64_t **_u, int32_t **_v, void *km);

#ifdef __cplusplus
}
#endif

#endif
