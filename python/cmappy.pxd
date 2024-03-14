from libc.stdint cimport int8_t, uint8_t, int32_t, int64_t, uint32_t, uint64_t

cdef extern from "minimap.h":
	#
	# Options
	#
	ctypedef struct mm_idxopt_t:
		short k, w, flag, bucket_bits
		int64_t mini_batch_size
		uint64_t batch_size

	ctypedef struct mm_mapopt_t:
		int64_t flag
		int seed
		int sdust_thres

		int max_qlen

		int bw, bw_long
		int max_gap, max_gap_ref
		int max_frag_len
		int max_chain_skip, max_chain_iter
		int min_cnt
		int min_chain_score
		float chain_gap_scale
		float chain_skip_scale
		int rmq_size_cap, rmq_inner_dist
		int rmq_rescue_size
		float rmq_rescue_ratio

		float mask_level
		int mask_len
		float pri_ratio
		int best_n

		float alt_drop

		int a, b, q, e, q2, e2
		int transition
		int sc_ambi
		int noncan
		int junc_bonus
		int zdrop, zdrop_inv
		int end_bonus
		int min_dp_max
		int min_ksw_len
		int anchor_ext_len, anchor_ext_shift
		float max_clip_ratio

		int rank_min_len
		float rank_frac

		int pe_ori, pe_bonus

		float mid_occ_frac
		float q_occ_frac
		int32_t min_mid_occ
		int32_t mid_occ
		int32_t max_occ
		int64_t mini_batch_size
		int64_t max_sw_mat
		int64_t cap_kalloc

		const char *split_prefix

	int mm_set_opt(char *preset, mm_idxopt_t *io, mm_mapopt_t *mo)
	int mm_verbose

	#
	# Indexing
	#
	ctypedef struct mm_idx_seq_t:
		char *name
		uint64_t offset
		uint32_t len

	ctypedef struct mm_idx_bucket_t:
		pass

	ctypedef struct mm_idx_t:
		int32_t b, w, k, flag
		uint32_t n_seq
		mm_idx_seq_t *seq
		uint32_t *S
		mm_idx_bucket_t *B
		void *km
		void *h

	ctypedef struct mm_idx_reader_t:
		pass

	mm_idx_reader_t *mm_idx_reader_open(const char *fn, const mm_idxopt_t *opt, const char *fn_out)
	mm_idx_t *mm_idx_reader_read(mm_idx_reader_t *r, int n_threads)
	void mm_idx_reader_close(mm_idx_reader_t *r)
	void mm_idx_destroy(mm_idx_t *mi)
	void mm_mapopt_update(mm_mapopt_t *opt, const mm_idx_t *mi)

	int mm_idx_index_name(mm_idx_t *mi)

	#
	# Mapping (key struct defined in cmappy.h below)
	#
	ctypedef struct mm_reg1_t:
		pass

	ctypedef struct mm_tbuf_t:
		pass

	mm_tbuf_t *mm_tbuf_init()
	void mm_tbuf_destroy(mm_tbuf_t *b)
	void *mm_tbuf_get_km(mm_tbuf_t *b)
	int mm_gen_cs(void *km, char **buf, int *max_len, const mm_idx_t *mi, const mm_reg1_t *r, const char *seq, int no_iden)
	int mm_gen_MD(void *km, char **buf, int *max_len, const mm_idx_t *mi, const mm_reg1_t *r, const char *seq)

#
# Helper header (because it is hard to expose mm_reg1_t with Cython)
#
cdef extern from "cmappy.h":
	ctypedef struct mm_hitpy_t:
		const char *ctg
		int32_t ctg_start, ctg_end
		int32_t qry_start, qry_end
		int32_t blen, mlen, NM, ctg_len
		uint8_t mapq, is_primary
		int8_t strand, trans_strand
		int32_t seg_id
		int32_t n_cigar32
		uint32_t *cigar32

	void mm_reg2hitpy(const mm_idx_t *mi, mm_reg1_t *r, mm_hitpy_t *h)
	void mm_free_reg1(mm_reg1_t *r)
	mm_reg1_t *mm_map_aux(const mm_idx_t *mi, const char *seq1, const char *seq2, int *n_regs, mm_tbuf_t *b, const mm_mapopt_t *opt)
	char *mappy_fetch_seq(const mm_idx_t *mi, const char *name, int st, int en, int *l)
	mm_idx_t *mappy_idx_seq(int w, int k, int is_hpc, int bucket_bits, const char *seq, int l)

	ctypedef struct kstring_t:
		unsigned l, m
		char *s

	ctypedef struct kstream_t:
		pass

	ctypedef struct kseq_t:
		kstring_t name, comment, seq, qual
		int last_char
		kstream_t *f

	kseq_t *mm_fastx_open(const char *fn)
	void mm_fastx_close(kseq_t *ks)
	int kseq_read(kseq_t *seq)

	char *mappy_revcomp(int l, const uint8_t *seq)
	int mm_verbose_level(int v)
	void mm_reset_timer()
