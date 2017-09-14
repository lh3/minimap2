#ifndef MINIMAP2_H
#define MINIMAP2_H

#include <stdint.h>
#include <stdio.h>
#include <sys/types.h>

#define MM_F_NO_SELF     0x001
#define MM_F_AVA         0x002
#define MM_F_CIGAR       0x004
#define MM_F_OUT_SAM     0x008
#define MM_F_NO_QUAL     0x010
#define MM_F_OUT_CG      0x020
#define MM_F_OUT_CS      0x040
#define MM_F_SPLICE      0x080
#define MM_F_SPLICE_FOR  0x100
#define MM_F_SPLICE_REV  0x200
#define MM_F_SPLICE_BOTH 0x400
#define MM_F_NO_SAM_SQ   0x800
#define MM_F_APPROX_EXT  0x1000

#define MM_IDX_MAGIC   "MMI\2"

#ifdef __cplusplus
extern "C" {
#endif

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
	void *km;
} mm_idx_t;

typedef struct {
	uint32_t capacity;
	int32_t dp_score, dp_max, dp_max2;
	uint32_t blen;
	uint32_t n_diff;
	uint32_t n_ambi:30, trans_strand:2;
	uint32_t n_cigar;
	uint32_t cigar[];
} mm_extra_t;

typedef struct {
	int32_t id;
	uint32_t cnt:31, rev:1;
	uint32_t rid:31, inv:1;
	int32_t score;
	int32_t qs, qe, rs, re;
	int32_t parent, subsc;
	int32_t as;
	int32_t fuzzy_mlen, fuzzy_blen;
	uint32_t mapq:8, split:2, sam_pri:1, n_sub:21; // TODO: n_sub is not used for now
	mm_extra_t *p;
} mm_reg1_t;

typedef struct {
	float max_occ_frac;
	float mid_occ_frac;
	int sdust_thres;  // score threshold for SDUST; 0 to disable
	int flag;    // see MM_F_* macros

	int bw;  // bandwidth
	int max_gap, max_gap_ref; // break a chain if there are no minimizers in a max_gap window
	int max_chain_skip;
	int min_cnt;
	int min_chain_score;

	float mask_level;
	float pri_ratio;
	int best_n;

	int max_join_long, max_join_short;
	int min_join_flank_sc;

	int a, b, q, e, q2, e2; // matching score, mismatch, gap-open and gap-ext penalties
	int noncan;
	int zdrop;
	int min_dp_max;
	int min_ksw_len;

	int mini_batch_size;
	int32_t max_occ, mid_occ;
} mm_mapopt_t;

typedef struct {
	short k, w, is_hpc, bucket_bits;
	int mini_batch_size;
	uint64_t batch_size;
} mm_idxopt_t;

struct mm_bseq_file_s;

typedef struct {
	int is_idx, n_parts;
	mm_idxopt_t opt;
	union {
		struct mm_bseq_file_s *seq;
		FILE *idx;
	} fp;
} mm_idx_reader_t;

extern int mm_verbose, mm_dbg_flag;
extern double mm_realtime0;

struct mm_tbuf_s;
typedef struct mm_tbuf_s mm_tbuf_t;

#define mm_seq4_set(s, i, c) ((s)[(i)>>3] |= (uint32_t)(c) << (((i)&7)<<2))
#define mm_seq4_get(s, i)    ((s)[(i)>>3] >> (((i)&7)<<2) & 0xf)

void mm_idxopt_init(mm_idxopt_t *opt);
void mm_mapopt_init(mm_mapopt_t *opt);
int mm_preset(const char *preset, mm_idxopt_t *io, mm_mapopt_t *mo);
void mm_mapopt_update(mm_mapopt_t *opt, const mm_idx_t *mi);

mm_idx_reader_t *mm_idx_reader_open(const char *fn, const mm_idxopt_t *opt);
int mm_idx_reader_is_idx(const mm_idx_reader_t *r);
mm_idx_t *mm_idx_reader_read(mm_idx_reader_t *r, int n_threads);
void mm_idx_reader_close(mm_idx_reader_t *r);
void mm_idx_destroy(mm_idx_t *mi);

// minimizer indexing
void mm_idx_stat(const mm_idx_t *idx);
const uint64_t *mm_idx_get(const mm_idx_t *mi, uint64_t minier, int *n);
int mm_idx_getseq(const mm_idx_t *mi, uint32_t rid, uint32_t st, uint32_t en, uint8_t *seq);

// minimizer index I/O
void mm_idx_dump(FILE *fp, const mm_idx_t *mi);
mm_idx_t *mm_idx_load(FILE *fp);

// mapping
mm_tbuf_t *mm_tbuf_init(void);
void mm_tbuf_destroy(mm_tbuf_t *b);
mm_reg1_t *mm_map(const mm_idx_t *mi, int l_seq, const char *seq, int *n_regs, mm_tbuf_t *b, const mm_mapopt_t *opt, const char *name);

int mm_map_file(const mm_idx_t *idx, const char *fn, const mm_mapopt_t *opt, int n_threads);

// obsolete APIs (for backward compatibility)
mm_idx_t *mm_idx_build(const char *fn, int w, int k, int is_hpc, int n_threads);

#ifdef __cplusplus
}
#endif

#endif // MINIMAP2_H
