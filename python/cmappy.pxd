from libc.stdint cimport int8_t, uint8_t, int32_t, int64_t, uint32_t, uint64_t

cdef extern from "minimap.h":
	#
	# Options
	#
	ctypedef struct mm_idxopt_t:
		short k, w, is_hpc, bucket_bits
		int mini_batch_size
		uint64_t batch_size

	ctypedef struct mm_mapopt_t:
		int seed
		int sdust_thres
		int flag
		int bw
		int max_gap, max_gap_ref
		int max_frag_len
		int max_chain_skip
		int min_cnt
		int min_chain_score
		float mask_level
		float pri_ratio
		int best_n
		int max_join_long, max_join_short
		int min_join_flank_sc
		int a, b, q, e, q2, e2
		int noncan
		int zdrop
		int end_bonus
		int min_dp_max
		int min_ksw_len
		int pe_ori, pe_bonus
		float mid_occ_frac
		int32_t mid_occ
		int32_t max_occ
		int mini_batch_size

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
		int32_t b, w, k, is_hpc
		uint32_t n_seq
		mm_idx_seq_t *seq
		uint32_t *S
		mm_idx_bucket_t *B
		void *km

	ctypedef struct mm_idx_reader_t:
		pass

	mm_idx_reader_t *mm_idx_reader_open(const char *fn, const mm_idxopt_t *opt, const char *fn_out)
	mm_idx_t *mm_idx_reader_read(mm_idx_reader_t *r, int n_threads)
	void mm_idx_reader_close(mm_idx_reader_t *r)
	void mm_idx_destroy(mm_idx_t *mi)
	void mm_mapopt_update(mm_mapopt_t *opt, const mm_idx_t *mi)

	#
	# Mapping (key struct defined in cmappy.h below)
	#
	ctypedef struct mm_reg1_t:
		pass

	ctypedef struct mm_tbuf_t:
		pass

	mm_tbuf_t *mm_tbuf_init()
	void mm_tbuf_destroy(mm_tbuf_t *b)
	mm_reg1_t *mm_map(const mm_idx_t *mi, int l_seq, const char *seq, int *n_regs, mm_tbuf_t *b, const mm_mapopt_t *opt, const char *name)

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
		int32_t n_cigar32
		uint32_t *cigar32

	void mm_reg2hitpy(const mm_idx_t *mi, mm_reg1_t *r, mm_hitpy_t *h)
	void mm_free_reg1(mm_reg1_t *r)

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

	int mm_verbose_level(int v)
	void mm_reset_timer()
