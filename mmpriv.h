#ifndef MMPRIV2_H
#define MMPRIV2_H

#include <assert.h>
#include "minimap.h"
#include "bseq.h"

#define MM_PARENT_UNSET   (-1)
#define MM_PARENT_TMP_PRI (-2)

#define MM_DBG_NO_KALLOC   0x1
#define MM_DBG_PRINT_QNAME 0x2

#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

#ifndef KSTRING_T
#define KSTRING_T kstring_t
typedef struct __kstring_t {
	unsigned l, m;
	char *s;
} kstring_t;
#endif

#ifdef __cplusplus
extern "C" {
#endif

double cputime(void);
double realtime(void);

void radix_sort_128x(mm128_t *beg, mm128_t *end);
void radix_sort_64(uint64_t *beg, uint64_t *end);
uint32_t ks_ksmall_uint32_t(size_t n, uint32_t arr[], size_t kk);

void mm_write_paf(kstring_t *s, const mm_idx_t *mi, const bseq1_t *t, const mm_reg1_t *r);
void mm_write_sam(kstring_t *s, const mm_idx_t *mi, const bseq1_t *t, const mm_reg1_t *r);
int mm_chain_dp(int max_dist, int bw, int max_skip, int min_cnt, int min_sc, int64_t n, mm128_t *a, uint64_t **_u, void *km);
mm_reg1_t *mm_align_skeleton(void *km, const mm_mapopt_t *opt, const mm_idx_t *mi, int qlen, const char *qstr, int *n_regs_, mm_reg1_t *regs, mm128_t *a);

mm_reg1_t *mm_gen_regs(void *km, int qlen, int n_u, uint64_t *u, mm128_t *a);
void mm_split_reg(mm_reg1_t *r, mm_reg1_t *r2, int n, int qlen, mm128_t *a);
void mm_sync_regs(void *km, int n_regs, mm_reg1_t *regs);
void mm_set_parent(void *km, float mask_level, int n, mm_reg1_t *r);
void mm_select_sub(void *km, float mask_level, float pri_ratio, int best_n, int *n_, mm_reg1_t *r);
void mm_filter_regs(void *km, const mm_mapopt_t *opt, int *n_regs, mm_reg1_t *regs);
void mm_join_long(void *km, const mm_mapopt_t *opt, int qlen, int *n_regs, mm_reg1_t *regs, mm128_t *a);
void mm_hit_sort_by_dp(void *km, int *n_regs, mm_reg1_t *r);
void mm_set_mapq(int n_regs, mm_reg1_t *regs);

#ifdef __cplusplus
}
#endif

#endif
