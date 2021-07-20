#ifndef MINIMAP2_H
#define MINIMAP2_H

#include <stdint.h>
#include <stdio.h>
#include <sys/types.h>

#define MM_F_NO_DIAG       0x001 // no exact diagonal hit
#define MM_F_NO_DUAL       0x002 // skip pairs where query name is lexicographically larger than target name
#define MM_F_CIGAR         0x004
#define MM_F_OUT_SAM       0x008
#define MM_F_NO_QUAL       0x010
#define MM_F_OUT_CG        0x020
#define MM_F_OUT_CS        0x040
#define MM_F_SPLICE        0x080 // splice mode
#define MM_F_SPLICE_FOR    0x100 // match GT-AG
#define MM_F_SPLICE_REV    0x200 // match CT-AC, the reverse complement of GT-AG
#define MM_F_NO_LJOIN      0x400
#define MM_F_OUT_CS_LONG   0x800
#define MM_F_SR            0x1000
#define MM_F_FRAG_MODE     0x2000
#define MM_F_NO_PRINT_2ND  0x4000
#define MM_F_2_IO_THREADS  0x8000
#define MM_F_LONG_CIGAR    0x10000
#define MM_F_INDEPEND_SEG  0x20000
#define MM_F_SPLICE_FLANK  0x40000
#define MM_F_SOFTCLIP      0x80000
#define MM_F_FOR_ONLY      0x100000
#define MM_F_REV_ONLY      0x200000
#define MM_F_HEAP_SORT     0x400000
#define MM_F_ALL_CHAINS    0x800000
#define MM_F_OUT_MD        0x1000000
#define MM_F_COPY_COMMENT  0x2000000
#define MM_F_EQX           0x4000000 // use =/X instead of M
#define MM_F_PAF_NO_HIT    0x8000000 // output unmapped reads to PAF
#define MM_F_NO_END_FLT    0x10000000
#define MM_F_HARD_MLEVEL   0x20000000
#define MM_F_SAM_HIT_ONLY  0x40000000
#define MM_F_RMQ           (0x80000000LL)
#define MM_F_QSTRAND       (0x100000000LL)
#define MM_F_NO_INV        (0x200000000LL)

#define MM_I_HPC          0x1
#define MM_I_NO_SEQ       0x2
#define MM_I_NO_NAME      0x4

#define MM_IDX_MAGIC   "MMI\2"

#define MM_MAX_SEG       255

#define MM_CIGAR_MATCH      0
#define MM_CIGAR_INS        1
#define MM_CIGAR_DEL        2
#define MM_CIGAR_N_SKIP     3
#define MM_CIGAR_SOFTCLIP   4
#define MM_CIGAR_HARDCLIP   5
#define MM_CIGAR_PADDING    6
#define MM_CIGAR_EQ_MATCH   7
#define MM_CIGAR_X_MISMATCH 8

#define MM_CIGAR_STR  "MIDNSHP=XB"

#ifdef __cplusplus
extern "C" {
#endif

// emulate 128-bit integers and arrays
typedef struct { uint64_t x, y; } mm128_t;
typedef struct { size_t n, m; mm128_t *a; } mm128_v;

// minimap2 index
typedef struct {
	char *name;      // name of the db sequence
	uint64_t offset; // offset in mm_idx_t::S
	uint32_t len;    // length
	uint32_t is_alt;
} mm_idx_seq_t;

typedef struct {
	int32_t b, w, k, flag;
	uint32_t n_seq;            // number of reference sequences
	int32_t index;
	int32_t n_alt;
	mm_idx_seq_t *seq;         // sequence name, length and offset
	uint32_t *S;               // 4-bit packed sequence
	struct mm_idx_bucket_s *B; // index (hidden)
	struct mm_idx_intv_s *I;   // intervals (hidden)
	void *km, *h;
} mm_idx_t;

// minimap2 alignment
typedef struct {
	uint32_t capacity;                  // the capacity of cigar[]
	int32_t dp_score, dp_max, dp_max2;  // DP score; score of the max-scoring segment; score of the best alternate mappings
	uint32_t n_ambi:30, trans_strand:2; // number of ambiguous bases; transcript strand: 0 for unknown, 1 for +, 2 for -
	uint32_t n_cigar;                   // number of cigar operations in cigar[]
	uint32_t cigar[];
} mm_extra_t;

typedef struct {
	int32_t id;             // ID for internal uses (see also parent below)
	int32_t cnt;            // number of minimizers; if on the reverse strand
	int32_t rid;            // reference index; if this is an alignment from inversion rescue
	int32_t score;          // DP alignment score
	int32_t qs, qe, rs, re; // query start and end; reference start and end
	int32_t parent, subsc;  // parent==id if primary; best alternate mapping score
	int32_t as;             // offset in the a[] array (for internal uses only)
	int32_t mlen, blen;     // seeded exact match length; seeded alignment block length
	int32_t n_sub;          // number of suboptimal mappings
	int32_t score0;         // initial chaining score (before chain merging/spliting)
	uint32_t mapq:8, split:2, rev:1, inv:1, sam_pri:1, proper_frag:1, pe_thru:1, seg_split:1, seg_id:8, split_inv:1, is_alt:1, dummy:6;
	uint32_t hash;
	float div;
	mm_extra_t *p;
} mm_reg1_t;

// indexing and mapping options
typedef struct {
	short k, w, flag, bucket_bits;
	int64_t mini_batch_size;
	uint64_t batch_size;
} mm_idxopt_t;

typedef struct {
	int64_t flag;    // see MM_F_* macros
	int seed;
	int sdust_thres; // score threshold for SDUST; 0 to disable

	int max_qlen;    // max query length

	int bw, bw_long; // bandwidth
	int max_gap, max_gap_ref; // break a chain if there are no minimizers in a max_gap window
	int max_frag_len;
	int max_chain_skip, max_chain_iter;
	int min_cnt;         // min number of minimizers on each chain
	int min_chain_score; // min chaining score
	float chain_gap_scale;
	int rmq_size_cap, rmq_inner_dist;
	int rmq_rescue_size;
	float rmq_rescue_ratio;

	float mask_level;
	int mask_len;
	float pri_ratio;
	int best_n;      // top best_n chains are subjected to DP alignment

	float alt_drop;

	int a, b, q, e, q2, e2; // matching score, mismatch, gap-open and gap-ext penalties
	int sc_ambi; // score when one or both bases are "N"
	int noncan;      // cost of non-canonical splicing sites
	int junc_bonus;
	int zdrop, zdrop_inv;   // break alignment if alignment score drops too fast along the diagonal
	int end_bonus;
	int min_dp_max;  // drop an alignment if the score of the max scoring segment is below this threshold
	int min_ksw_len;
	int anchor_ext_len, anchor_ext_shift;
	float max_clip_ratio; // drop an alignment if BOTH ends are clipped above this ratio

	int rank_min_len;
	float rank_frac;

	int pe_ori, pe_bonus;

	float mid_occ_frac;  // only used by mm_mapopt_update(); see below
	int32_t min_mid_occ, max_mid_occ;
	int32_t mid_occ;     // ignore seeds with occurrences above this threshold
	int32_t max_occ, max_max_occ, occ_dist;
	int64_t mini_batch_size; // size of a batch of query bases to process in parallel
	int64_t max_sw_mat;
	int64_t cap_kalloc;

	const char *split_prefix;
} mm_mapopt_t;

// index reader
typedef struct {
	int is_idx, n_parts;
	int64_t idx_size;
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

/**
 * Set default or preset parameters
 *
 * @param preset     NULL to set all parameters as default; otherwise apply preset to affected parameters
 * @param io         pointer to indexing parameters
 * @param mo         pointer to mapping parameters
 *
 * @return 0 if success; -1 if _present_ unknown
 */
int mm_set_opt(const char *preset, mm_idxopt_t *io, mm_mapopt_t *mo);
int mm_check_opt(const mm_idxopt_t *io, const mm_mapopt_t *mo);

/**
 * Update mm_mapopt_t::mid_occ via mm_mapopt_t::mid_occ_frac
 *
 * If mm_mapopt_t::mid_occ is 0, this function sets it to a number such that no
 * more than mm_mapopt_t::mid_occ_frac of minimizers in the index have a higher
 * occurrence.
 *
 * @param opt        mapping parameters
 * @param mi         minimap2 index
 */
void mm_mapopt_update(mm_mapopt_t *opt, const mm_idx_t *mi);

void mm_mapopt_max_intron_len(mm_mapopt_t *opt, int max_intron_len);

/**
 * Initialize an index reader
 *
 * @param fn         index or fasta/fastq file name (this function tests the file type)
 * @param opt        indexing parameters
 * @param fn_out     if not NULL, write built index to this file
 *
 * @return an index reader on success; NULL if fail to open _fn_
 */
mm_idx_reader_t *mm_idx_reader_open(const char *fn, const mm_idxopt_t *opt, const char *fn_out);

/**
 * Read/build an index
 *
 * If the input file is an index file, this function reads one part of the
 * index and returns. If the input file is a sequence file (fasta or fastq),
 * this function constructs the index for about mm_idxopt_t::batch_size bases.
 * Importantly, for a huge collection of sequences, this function may only
 * return an index for part of sequences. It needs to be repeatedly called
 * to traverse the entire index/sequence file.
 *
 * @param r          index reader
 * @param n_threads  number of threads for constructing index
 *
 * @return an index on success; NULL if reaching the end of the input file
 */
mm_idx_t *mm_idx_reader_read(mm_idx_reader_t *r, int n_threads);

/**
 * Destroy/deallocate an index reader
 *
 * @param r          index reader
 */
void mm_idx_reader_close(mm_idx_reader_t *r);

int mm_idx_reader_eof(const mm_idx_reader_t *r);

/**
 * Check whether the file contains a minimap2 index
 *
 * @param fn         file name
 *
 * @return the file size if fn is an index file; 0 if fn is not.
 */
int64_t mm_idx_is_idx(const char *fn);

/**
 * Load a part of an index
 *
 * Given a uni-part index, this function loads the entire index into memory.
 * Given a multi-part index, it loads one part only and places the file pointer
 * at the end of that part.
 *
 * @param fp         pointer to FILE object
 *
 * @return minimap2 index read from fp
 */
mm_idx_t *mm_idx_load(FILE *fp);

/**
 * Append an index (or one part of a full index) to file
 *
 * @param fp         pointer to FILE object
 * @param mi         minimap2 index
 */
void mm_idx_dump(FILE *fp, const mm_idx_t *mi);

/**
 * Create an index from strings in memory
 *
 * @param w            minimizer window size
 * @param k            minimizer k-mer size
 * @param is_hpc       use HPC k-mer if true
 * @param bucket_bits  number of bits for the first level of the hash table
 * @param n            number of sequences
 * @param seq          sequences in A/C/G/T
 * @param name         sequence names; could be NULL
 *
 * @return minimap2 index
 */
mm_idx_t *mm_idx_str(int w, int k, int is_hpc, int bucket_bits, int n, const char **seq, const char **name);

/**
 * Print index statistics to stderr
 *
 * @param mi         minimap2 index
 */
void mm_idx_stat(const mm_idx_t *idx);

/**
 * Destroy/deallocate an index
 *
 * @param r          minimap2 index
 */
void mm_idx_destroy(mm_idx_t *mi);

/**
 * Initialize a thread-local buffer for mapping
 *
 * Each mapping thread requires a buffer specific to the thread (see mm_map()
 * below). The primary purpose of this buffer is to reduce frequent heap
 * allocations across threads. A buffer shall not be used by two or more
 * threads.
 *
 * @return pointer to a thread-local buffer
 */
mm_tbuf_t *mm_tbuf_init(void);

/**
 * Destroy/deallocate a thread-local buffer for mapping
 *
 * @param b          the buffer
 */
void mm_tbuf_destroy(mm_tbuf_t *b);

void *mm_tbuf_get_km(mm_tbuf_t *b);

/**
 * Align a query sequence against an index
 *
 * This function possibly finds multiple alignments of the query sequence.
 * The returned array and the mm_reg1_t::p field of each element are allocated
 * with malloc().
 *
 * @param mi         minimap2 index
 * @param l_seq      length of the query sequence
 * @param seq        the query sequence
 * @param n_regs     number of hits (out)
 * @param b          thread-local buffer; two mm_map() calls shall not use one buffer at the same time!
 * @param opt        mapping parameters
 * @param name       query name, used for all-vs-all overlapping and debugging
 *
 * @return an array of hits which need to be deallocated with free() together
 *         with mm_reg1_t::p of each element. The size is written to _n_regs_.
 */
mm_reg1_t *mm_map(const mm_idx_t *mi, int l_seq, const char *seq, int *n_regs, mm_tbuf_t *b, const mm_mapopt_t *opt, const char *name);

void mm_map_frag(const mm_idx_t *mi, int n_segs, const int *qlens, const char **seqs, int *n_regs, mm_reg1_t **regs, mm_tbuf_t *b, const mm_mapopt_t *opt, const char *qname);

/**
 * Align a fasta/fastq file and print alignments to stdout
 *
 * @param idx        minimap2 index
 * @param fn         fasta/fastq file name
 * @param opt        mapping parameters
 * @param n_threads  number of threads
 *
 * @return 0 on success; -1 if _fn_ can't be read
 */
int mm_map_file(const mm_idx_t *idx, const char *fn, const mm_mapopt_t *opt, int n_threads);

int mm_map_file_frag(const mm_idx_t *idx, int n_segs, const char **fn, const mm_mapopt_t *opt, int n_threads);

/**
 * Generate the cs tag (new in 2.12)
 *
 * @param km         memory blocks; set to NULL if unsure
 * @param buf        buffer to write the cs/MD tag; typicall NULL on the first call
 * @param max_len    max length of the buffer; typically set to 0 on the first call
 * @param mi         index
 * @param r          alignment
 * @param seq        query sequence
 * @param no_iden    true to use : instead of =
 *
 * @return the length of cs
 */
int mm_gen_cs(void *km, char **buf, int *max_len, const mm_idx_t *mi, const mm_reg1_t *r, const char *seq, int no_iden);
int mm_gen_MD(void *km, char **buf, int *max_len, const mm_idx_t *mi, const mm_reg1_t *r, const char *seq);

// query sequence name and sequence in the minimap2 index
int mm_idx_index_name(mm_idx_t *mi);
int mm_idx_name2id(const mm_idx_t *mi, const char *name);
int mm_idx_getseq(const mm_idx_t *mi, uint32_t rid, uint32_t st, uint32_t en, uint8_t *seq);

int mm_idx_alt_read(mm_idx_t *mi, const char *fn);
int mm_idx_bed_read(mm_idx_t *mi, const char *fn, int read_junc);
int mm_idx_bed_junc(const mm_idx_t *mi, int32_t ctg, int32_t st, int32_t en, uint8_t *s);

// deprecated APIs for backward compatibility
void mm_mapopt_init(mm_mapopt_t *opt);
mm_idx_t *mm_idx_build(const char *fn, int w, int k, int flag, int n_threads);

#ifdef __cplusplus
}
#endif

#endif // MINIMAP2_H
