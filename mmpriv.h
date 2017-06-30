#ifndef MMPRIV2_H
#define MMPRIV2_H

#include <assert.h>
#include "minimap.h"
#include "bseq.h"

#define MM_PARENT_UNSET   (-1)
#define MM_PARENT_TMP_PRI (-2)

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

mm_reg1_t *mm_gen_regs(int qlen, int n_u, uint64_t *u, mm128_t *a);
void mm_split_reg(mm_reg1_t *r, mm_reg1_t *r2, int n, int qlen, mm128_t *a);
void mm_sync_regs(void *km, int n_regs, mm_reg1_t *regs);
void mm_select_sub(void *km, float mask_level, float pri_ratio, int *n_, mm_reg1_t *r);
void mm_join_long(void *km, const mm_mapopt_t *opt, int qlen, int n_regs, mm_reg1_t *regs, mm128_t *a);
void mm_set_mapq(int n_regs, mm_reg1_t *regs);

#ifdef __cplusplus
}
#endif

static inline void mm_reg_set_coor(mm_reg1_t *r, int32_t qlen, const mm128_t *a)
{
	int32_t k = r->as, q_span = (int32_t)(a[k].y>>32&0xff);
	r->rev = a[k].x>>63;
	r->rid = a[k].x<<1>>33;
	r->rs = (int32_t)a[k].x + 1 > q_span? (int32_t)a[k].x + 1 - q_span : 0; // NB: target span may be shorter, so this test is necessary
	r->re = (int32_t)a[k + r->cnt - 1].x + 1;
	if (!r->rev) {
		r->qs = (int32_t)a[k].y + 1 - q_span;
		r->qe = (int32_t)a[k + r->cnt - 1].y + 1;
	} else {
		r->qs = qlen - ((int32_t)a[k + r->cnt - 1].y + 1);
		r->qe = qlen - ((int32_t)a[k].y + 1 - q_span);
	}
}

#endif
