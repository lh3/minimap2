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

int mm_chain_dp(int max_dist, int bw, int max_skip, int min_sc, int n, mm128_t *a, uint64_t **_u, void *km);

void mm_align_skeleton(void *km, const mm_mapopt_t *opt, const mm_idx_t *mi, int qlen, const char *qstr, int n_regs, mm_reg1_t *regs, mm128_t *a);

#ifdef __cplusplus
}
#endif

#endif
