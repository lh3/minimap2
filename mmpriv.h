#ifndef MMPRIV2_H
#define MMPRIV2_H

#include <assert.h>
#include "minimap.h"
#include "bseq.h"
#include "kseq.h"

#define MM_PARENT_UNSET   (-1)
#define MM_PARENT_TMP_PRI (-2)

#define MM_DBG_NO_KALLOC     0x1
#define MM_DBG_PRINT_QNAME   0x2
#define MM_DBG_PRINT_SEED    0x4
#define MM_DBG_PRINT_ALN_SEQ 0x8

#define MM_SEED_LONG_JOIN  (1ULL<<40)
#define MM_SEED_IGNORE     (1ULL<<41)
#define MM_SEED_TANDEM     (1ULL<<42)
#define MM_SEED_SELF       (1ULL<<43)

#define MM_SEED_SEG_SHIFT  48
#define MM_SEED_SEG_MASK   (0xffULL<<(MM_SEED_SEG_SHIFT))

#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

#define mm_seq4_set(s, i, c) ((s)[(i)>>3] |= (uint32_t)(c) << (((i)&7)<<2))
#define mm_seq4_get(s, i)    ((s)[(i)>>3] >> (((i)&7)<<2) & 0xf)

#define MALLOC(type, len) ((type*)malloc((len) * sizeof(type)))
#define CALLOC(type, len) ((type*)calloc((len), sizeof(type)))

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
	uint32_t n;
	uint32_t q_pos;
	uint32_t q_span:31, flt:1;
	uint32_t seg_id:31, is_tandem:1;
	const uint64_t *cr;
} mm_seed_t;

typedef struct {
	int n_u, n_a;
	uint64_t *u;
	mm128_t *a;
} mm_seg_t;

double cputime(void);
double realtime(void);
long peakrss(void);

void radix_sort_128x(mm128_t *beg, mm128_t *end);
void radix_sort_64(uint64_t *beg, uint64_t *end);
uint32_t ks_ksmall_uint32_t(size_t n, uint32_t arr[], size_t kk);

void mm_sketch(void *km, const char *str, int len, int w, int k, uint32_t rid, int is_hpc, mm128_v *p);

mm_seed_t *mm_collect_matches(void *km, int *_n_m, int qlen, int max_occ, int max_max_occ, int dist, const mm_idx_t *mi, const mm128_v *mv, int64_t *n_a, int *rep_len, int *n_mini_pos, uint64_t **mini_pos);

double mm_event_identity(const mm_reg1_t *r);
int mm_write_sam_hdr(const mm_idx_t *mi, const char *rg, const char *ver, int argc, char *argv[]);
void mm_write_paf(kstring_t *s, const mm_idx_t *mi, const mm_bseq1_t *t, const mm_reg1_t *r, void *km, int64_t opt_flag);
void mm_write_paf3(kstring_t *s, const mm_idx_t *mi, const mm_bseq1_t *t, const mm_reg1_t *r, void *km, int64_t opt_flag, int rep_len);
void mm_write_sam(kstring_t *s, const mm_idx_t *mi, const mm_bseq1_t *t, const mm_reg1_t *r, int n_regs, const mm_reg1_t *regs);
void mm_write_sam2(kstring_t *s, const mm_idx_t *mi, const mm_bseq1_t *t, int seg_idx, int reg_idx, int n_seg, const int *n_regs, const mm_reg1_t *const* regs, void *km, int64_t opt_flag);
void mm_write_sam3(kstring_t *s, const mm_idx_t *mi, const mm_bseq1_t *t, int seg_idx, int reg_idx, int n_seg, const int *n_regss, const mm_reg1_t *const* regss, void *km, int64_t opt_flag, int rep_len);

void mm_idxopt_init(mm_idxopt_t *opt);
const uint64_t *mm_idx_get(const mm_idx_t *mi, uint64_t minier, int *n);
int32_t mm_idx_cal_max_occ(const mm_idx_t *mi, float f);
int mm_idx_getseq2(const mm_idx_t *mi, int is_rev, uint32_t rid, uint32_t st, uint32_t en, uint8_t *seq);
mm_reg1_t *mm_align_skeleton(void *km, const mm_mapopt_t *opt, const mm_idx_t *mi, int qlen, const char *qstr, int *n_regs_, mm_reg1_t *regs, mm128_t *a);
mm_reg1_t *mm_gen_regs(void *km, uint32_t hash, int qlen, int n_u, uint64_t *u, mm128_t *a, int is_qstrand);

mm128_t *mm_chain_dp(int max_dist_x, int max_dist_y, int bw, int max_skip, int max_iter, int min_cnt, int min_sc, float gap_scale,
					 int is_cdna, int n_segs, int64_t n, mm128_t *a, int *n_u_, uint64_t **_u, void *km);
mm128_t *mg_lchain_dp(int max_dist_x, int max_dist_y, int bw, int max_skip, int max_iter, int min_cnt, int min_sc, float chn_pen_gap, float chn_pen_skip,
					  int is_cdna, int n_segs, int64_t n, mm128_t *a, int *n_u_, uint64_t **_u, void *km);
mm128_t *mg_lchain_rmq(int max_dist, int max_dist_inner, int bw, int max_chn_skip, int cap_rmq_size, int min_cnt, int min_sc, float chn_pen_gap, float chn_pen_skip,
					   int64_t n, mm128_t *a, int *n_u_, uint64_t **_u, void *km);

void mm_mark_alt(const mm_idx_t *mi, int n, mm_reg1_t *r);
void mm_split_reg(mm_reg1_t *r, mm_reg1_t *r2, int n, int qlen, mm128_t *a, int is_qstrand);
void mm_sync_regs(void *km, int n_regs, mm_reg1_t *regs);
int mm_squeeze_a(void *km, int n_regs, mm_reg1_t *regs, mm128_t *a);
int mm_set_sam_pri(int n, mm_reg1_t *r);
void mm_set_parent(void *km, float mask_level, int mask_len, int n, mm_reg1_t *r, int sub_diff, int hard_mask_level, float alt_diff_frac);
void mm_select_sub(void *km, float pri_ratio, int min_diff, int best_n, int *n_, mm_reg1_t *r);
void mm_select_sub_multi(void *km, float pri_ratio, float pri1, float pri2, int max_gap_ref, int min_diff, int best_n, int n_segs, const int *qlens, int *n_, mm_reg1_t *r);
void mm_filter_regs(const mm_mapopt_t *opt, int qlen, int *n_regs, mm_reg1_t *regs);
void mm_hit_sort(void *km, int *n_regs, mm_reg1_t *r, float alt_diff_frac);
void mm_set_mapq(void *km, int n_regs, mm_reg1_t *regs, int min_chain_sc, int match_sc, int rep_len, int is_sr);
void mm_update_dp_max(int qlen, int n_regs, mm_reg1_t *regs, float frac, int a, int b);

void mm_est_err(const mm_idx_t *mi, int qlen, int n_regs, mm_reg1_t *regs, const mm128_t *a, int32_t n, const uint64_t *mini_pos);

mm_seg_t *mm_seg_gen(void *km, uint32_t hash, int n_segs, const int *qlens, int n_regs0, const mm_reg1_t *regs0, int *n_regs, mm_reg1_t **regs, const mm128_t *a);
void mm_seg_free(void *km, int n_segs, mm_seg_t *segs);
void mm_pair(void *km, int max_gap_ref, int dp_bonus, int sub_diff, int match_sc, const int *qlens, int *n_regs, mm_reg1_t **regs);

FILE *mm_split_init(const char *prefix, const mm_idx_t *mi);
mm_idx_t *mm_split_merge_prep(const char *prefix, int n_splits, FILE **fp, uint32_t *n_seq_part);
int mm_split_merge(int n_segs, const char **fn, const mm_mapopt_t *opt, int n_split_idx);
void mm_split_rm_tmp(const char *prefix, int n_splits);

void mm_err_puts(const char *str);
void mm_err_fwrite(const void *p, size_t size, size_t nitems, FILE *fp);
void mm_err_fread(void *p, size_t size, size_t nitems, FILE *fp);

static inline float mg_log2(float x) // NB: this doesn't work when x<2
{
	union { float f; uint32_t i; } z = { x };
	float log_2 = ((z.i >> 23) & 255) - 128;
	z.i &= ~(255 << 23);
	z.i += 127 << 23;
	log_2 += (-0.34484843f * z.f + 2.02466578f) * z.f - 0.67487759f;
	return log_2;
}

#ifdef __cplusplus
}
#endif

#endif
