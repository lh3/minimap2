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

// hidden structures
struct mm_idx_bucket_s;
struct mm_bseq_file_s;
struct mm_tbuf_s;

// emulate 128-bit integers and arrays
typedef struct { uint64_t x, y; } mm128_t;
typedef struct { size_t n, m; mm128_t *a; } mm128_v;

// minimap2 index
typedef struct {
	char *name;      // name of the db sequence
	uint64_t offset; // offset in mm_idx_t::S
	uint32_t len;    // length
} mm_idx_seq_t;

typedef struct {
	int32_t b, w, k, is_hpc;
	uint32_t n_seq;            // number of reference sequences
	mm_idx_seq_t *seq;         // sequence name, length and offset
	uint32_t *S;               // 4-bit packed sequence
	struct mm_idx_bucket_s *B; // index (hidden)
	void *km;
} mm_idx_t;

// minimap2 alignment
typedef struct {
	uint32_t capacity;                  // the capacity of cigar[]
	int32_t dp_score, dp_max, dp_max2;  // DP score; score of the max-scoring segment; score of the best alternate mappings
	uint32_t blen;                      // block length
	uint32_t n_diff;                    // number of differences, including ambiguous bases
	uint32_t n_ambi:30, trans_strand:2; // number of ambiguous bases; transcript strand: 0 for unknown, 1 for +, 2 for -
	uint32_t n_cigar;                   // number of cigar operations in cigar[]
	uint32_t cigar[];
} mm_extra_t;

typedef struct {
	int32_t id;                     // ID for internal uses (see also parent below)
	uint32_t cnt:31, rev:1;         // number of minimizers; if on the reverse strand
	uint32_t rid:31, inv:1;         // reference index; if this is an alignment from inversion rescue
	int32_t score;                  // DP alignment score
	int32_t qs, qe, rs, re;         // query start and end; reference start and end
	int32_t parent, subsc;          // parent==id if primary; best alternate mapping score
	int32_t as;                     // offset in the a[] array (for internal uses only)
	int32_t fuzzy_mlen, fuzzy_blen; // seeded exact match length; seeded alignment block length (approximate)
	uint32_t mapq:8, split:2, sam_pri:1, n_sub:21; // mapQ; split pattern; if SAM primary; number of suboptimal mappings
	mm_extra_t *p;
} mm_reg1_t;

// indexing and mapping options
typedef struct {
	short k, w, is_hpc, bucket_bits;
	int mini_batch_size;
	uint64_t batch_size;
} mm_idxopt_t;

typedef struct {
	float mid_occ_frac;
	int sdust_thres; // score threshold for SDUST; 0 to disable
	int flag;        // see MM_F_* macros

	int bw;          // bandwidth
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
	int32_t mid_occ;
} mm_mapopt_t;

// index reader
typedef struct {
	int is_idx, n_parts;
	mm_idxopt_t opt;
	FILE *fp_out;
	union {
		struct mm_bseq_file_s *seq;
		FILE *idx;
	} fp;
} mm_idx_reader_t;

// memory buffer for thread-local storage during mapping
typedef struct mm_tbuf_s mm_tbuf_t;

// global variables
extern int mm_verbose, mm_dbg_flag; // verbose level: 0 for no info, 1 for error, 2 for warning, 3 for message (default); debugging flag
extern double mm_realtime0; // wall-clock timer

int mm_set_opt(const char *preset, mm_idxopt_t *io, mm_mapopt_t *mo);
void mm_mapopt_update(mm_mapopt_t *opt, const mm_idx_t *mi);

mm_idx_reader_t *mm_idx_reader_open(const char *fn, const mm_idxopt_t *opt, const char *fn_out);
mm_idx_t *mm_idx_reader_read(mm_idx_reader_t *r, int n_threads);
void mm_idx_reader_close(mm_idx_reader_t *r);
void mm_idx_stat(const mm_idx_t *idx);
int mm_idx_getseq(const mm_idx_t *mi, uint32_t rid, uint32_t st, uint32_t en, uint8_t *seq);
void mm_idx_destroy(mm_idx_t *mi);

mm_tbuf_t *mm_tbuf_init(void);
void mm_tbuf_destroy(mm_tbuf_t *b);
mm_reg1_t *mm_map(const mm_idx_t *mi, int l_seq, const char *seq, int *n_regs, mm_tbuf_t *b, const mm_mapopt_t *opt, const char *name);
int mm_map_file(const mm_idx_t *idx, const char *fn, const mm_mapopt_t *opt, int n_threads);

// deprecated APIs for backward compatibility
void mm_mapopt_init(mm_mapopt_t *opt);
mm_idx_t *mm_idx_build(const char *fn, int w, int k, int is_hpc, int n_threads);

#ifdef __cplusplus
}
#endif

#endif // MINIMAP2_H
