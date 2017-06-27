#ifndef MINIMAP2_H
#define MINIMAP2_H

#include <stdint.h>
#include <stdio.h>
#include <sys/types.h>

#define MM_IDX_DEF_B    14

#define MM_F_NO_SELF   0x01
#define MM_F_AVA       0x02
#define MM_F_CIGAR     0x04
#define MM_F_OUT_SAM   0x08

#define MM_IDX_MAGIC   "MMI\2"

typedef struct {
 	uint64_t x, y;
} mm128_t;

typedef struct { size_t n, m; mm128_t *a; } mm128_v;
typedef struct { size_t n, m; uint64_t *a; } uint64_v;
typedef struct { size_t n, m; uint32_t *a; } uint32_v;

typedef struct {
	mm128_v a;   // (minimizer, position) array
	int32_t n;   // size of the _p_ array
	uint64_t *p; // position array for minimizers appearing >1 times
	void *h;     // hash table indexing _p_ and minimizers appearing once
} mm_idx_bucket_t;

typedef struct {
	char *name;      // name of the db sequence
	uint64_t offset; // offset in mm_idx_t::S
	uint32_t len;    // length
} mm_idx_seq_t;

typedef struct {
	int32_t b, w, k, is_hpc;
	uint32_t n_seq;     // number of reference sequences
	mm_idx_seq_t *seq;  // sequence name, length and offset
	uint32_t *S;        // 4-bit packed sequence
	mm_idx_bucket_t *B; // index
} mm_idx_t;

typedef struct {
	uint32_t capacity;
	int32_t score;
	uint32_t blen;
	uint32_t n_diff, n_ambi;
	uint32_t n_cigar;
	uint32_t cigar[];
} mm_extra_t;

typedef struct {
	uint32_t cnt:30, rev:1, split:1;
	uint32_t rid:31, rep:1;
	int32_t score;
	int32_t qs, qe, rs, re;
	int32_t parent, subsc;
	int32_t as;
	int32_t mapq, n_sub; // TODO: n_sub is not used for now
	mm_extra_t *p;
} mm_reg1_t;

typedef struct {
	float max_occ_frac;
	float mid_occ_frac;
	int sdust_thres;  // score threshold for SDUST; 0 to disable
	int flag;    // see MM_F_* macros
	int bw;  // bandwidth
	int max_gap; // break a chain if there are no minimizers in a max_gap window
	int max_skip;
	int min_cnt, min_score;
	float pri_ratio;
	float mask_level;
	int a, b, q, e; // matching score, mismatch, gap-open and gap-ext penalties
	int zdrop;

	int max_occ;
	int mid_occ;
	int min_ksw_len;
} mm_mapopt_t;

extern int mm_verbose;
extern double mm_realtime0;

struct mm_tbuf_s;
typedef struct mm_tbuf_s mm_tbuf_t;

struct bseq_file_s;

#define mm_seq4_set(s, i, c) ((s)[(i)>>3] |= (uint32_t)(c) << (((i)&7)<<2))
#define mm_seq4_get(s, i)    ((s)[(i)>>3] >> (((i)&7)<<2) & 0xf)

#ifdef __cplusplus
extern "C" {
#endif

// compute minimizers
void mm_sketch(void *km, const char *str, int len, int w, int k, uint32_t rid, int is_hpc, mm128_v *p);

// minimizer indexing
mm_idx_t *mm_idx_init(int w, int k, int b, int is_hpc);
void mm_idx_destroy(mm_idx_t *mi);
mm_idx_t *mm_idx_gen(struct bseq_file_s *fp, int w, int k, int b, int is_hpc, int mini_batch_size, int n_threads, uint64_t batch_size, int keep_name);
uint32_t mm_idx_cal_max_occ(const mm_idx_t *mi, float f);
void mm_idx_stat(const mm_idx_t *idx);
const uint64_t *mm_idx_get(const mm_idx_t *mi, uint64_t minier, int *n);
int mm_idx_getseq(const mm_idx_t *mi, uint32_t rid, uint32_t st, uint32_t en, uint8_t *seq);

mm_idx_t *mm_idx_build(const char *fn, int w, int k, int is_hpc, int n_threads);

// minimizer index I/O
void mm_idx_dump(FILE *fp, const mm_idx_t *mi);
mm_idx_t *mm_idx_load(FILE *fp);

// mapping
void mm_mapopt_init(mm_mapopt_t *opt);
void mm_mapopt_update(mm_mapopt_t *opt, const mm_idx_t *mi);
mm_tbuf_t *mm_tbuf_init(void);
void mm_tbuf_destroy(mm_tbuf_t *b);
mm_reg1_t *mm_map(const mm_idx_t *mi, int l_seq, const char *seq, int *n_regs, mm_tbuf_t *b, const mm_mapopt_t *opt, const char *name);

int mm_map_file(const mm_idx_t *idx, const char *fn, const mm_mapopt_t *opt, int n_threads, int tbatch_size);

#ifdef __cplusplus
}
#endif

#endif // MINIMAP2_H
