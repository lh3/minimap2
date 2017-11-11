from libc.stdint cimport uint8_t, int8_t
from libc.stdlib cimport free
cimport cmappy

cmappy.mm_reset_timer()

cdef class Alignment:
	cdef int _ctg_len, _r_st, _r_en
	cdef int _q_st, _q_en
	cdef int _NM, _mlen, _blen
	cdef int8_t _strand, _trans_strand
	cdef uint8_t _mapq, _is_primary
	cdef _ctg, _cigar # these are python objects

	def __cinit__(self, ctg, cl, cs, ce, strand, qs, qe, mapq, cigar, is_primary, mlen, blen, NM, trans_strand):
		self._ctg = ctg if isinstance(ctg, str) else ctg.decode()
		self._ctg_len, self._r_st, self._r_en = cl, cs, ce
		self._strand, self._q_st, self._q_en = strand, qs, qe
		self._NM, self._mlen, self._blen = NM, mlen, blen
		self._mapq = mapq
		self._cigar = cigar
		self._is_primary = is_primary
		self._trans_strand = trans_strand

	@property
	def ctg(self): return self._ctg

	@property
	def ctg_len(self): return self._ctg_len

	@property
	def r_st(self): return self._r_st

	@property
	def r_en(self): return self._r_en

	@property
	def strand(self): return self._strand

	@property
	def trans_strand(self): return self._trans_strand

	@property
	def blen(self): return self._blen

	@property
	def mlen(self): return self._mlen

	@property
	def NM(self): return self._NM

	@property
	def is_primary(self): return (self._is_primary != 0)

	@property
	def q_st(self): return self._q_st

	@property
	def q_en(self): return self._q_en

	@property
	def mapq(self): return self._mapq

	@property
	def cigar(self): return self._cigar

	@property
	def cigar_str(self):
		return "".join(map(lambda x: str(x[0]) + 'MIDNSH'[x[1]], self._cigar))

	def __str__(self):
		if self._strand > 0: strand = '+'
		elif self._strand < 0: strand = '-'
		else: strand = '?'
		if self._is_primary != 0: tp = 'tp:A:P'
		else: tp = 'tp:A:S'
		if self._trans_strand > 0: ts = 'ts:A:+'
		elif self._trans_strand < 0: ts = 'ts:A:-'
		else: ts = 'ts:A:.'
		return "\t".join([str(self._q_st), str(self._q_en), strand, self._ctg, str(self._ctg_len), str(self._r_st), str(self._r_en),
				str(self._mlen), str(self._blen), str(self._mapq), tp, ts, "cg:Z:" + self.cigar_str])

cdef class ThreadBuffer:
	cdef cmappy.mm_tbuf_t *_b

	def __cinit__(self):
		self._b = cmappy.mm_tbuf_init()

	def __dealloc__(self):
		cmappy.mm_tbuf_destroy(self._b)

cdef class Aligner:
	cdef cmappy.mm_idx_t *_idx
	cdef cmappy.mm_idxopt_t idx_opt
	cdef cmappy.mm_mapopt_t map_opt

	def __cinit__(self, fn_idx_in, preset=None, k=None, w=None, min_cnt=None, min_chain_score=None, min_dp_score=None, bw=None, best_n=None, n_threads=3, fn_idx_out=None):
		cmappy.mm_set_opt(NULL, &self.idx_opt, &self.map_opt) # set the default options
		if preset is not None:
			cmappy.mm_set_opt(str.encode(preset), &self.idx_opt, &self.map_opt) # apply preset
		self.map_opt.flag |= 4 # always perform alignment
		self.idx_opt.batch_size = 0x7fffffffffffffffL # always build a uni-part index
		if k is not None: self.idx_opt.k = k
		if w is not None: self.idx_opt.w = w
		if min_cnt is not None: self.map_opt.min_cnt = min_cnt
		if min_chain_score is not None: self.map_opt.min_chain_score = min_chain_score
		if min_dp_score is not None: self.map_opt.min_dp_max = min_dp_score
		if bw is not None: self.map_opt.bw = bw
		if best_n is not None: self.best_n = best_n

		cdef cmappy.mm_idx_reader_t *r;
		if fn_idx_out is None:
			r = cmappy.mm_idx_reader_open(str.encode(fn_idx_in), &self.idx_opt, NULL)
		else:
			r = cmappy.mm_idx_reader_open(str.encode(fn_idx_in), &self.idx_opt, fn_idx_out)
		if r is not NULL:
			self._idx = cmappy.mm_idx_reader_read(r, n_threads) # NB: ONLY read the first part
			cmappy.mm_idx_reader_close(r)
			cmappy.mm_mapopt_update(&self.map_opt, self._idx)

	def __dealloc__(self):
		if self._idx is not NULL:
			cmappy.mm_idx_destroy(self._idx)

	def __bool__(self):
		return (self._idx != NULL)

	def map(self, seq, buf=None):
		cdef cmappy.mm_reg1_t *regs
		cdef cmappy.mm_hitpy_t h
		cdef ThreadBuffer b
		cdef int n_regs

		if self._idx is NULL: return None
		if buf is None: b = ThreadBuffer()
		else: b = buf
		regs = cmappy.mm_map(self._idx, len(seq), str.encode(seq), &n_regs, b._b, &self.map_opt, NULL)

		for i in range(n_regs):
			cmappy.mm_reg2hitpy(self._idx, &regs[i], &h)
			cigar = []
			for k in range(h.n_cigar32):
				c = h.cigar32[k]
				cigar.append([c>>4, c&0xf])
			yield Alignment(h.ctg, h.ctg_len, h.ctg_start, h.ctg_end, h.strand, h.qry_start, h.qry_end, h.mapq, cigar, h.is_primary, h.mlen, h.blen, h.NM, h.trans_strand)
			cmappy.mm_free_reg1(&regs[i])
		free(regs)

def fastx_read(fn):
	cdef cmappy.kseq_t *ks
	ks = cmappy.mm_fastx_open(str.encode(fn))
	if ks is NULL: return None
	while cmappy.kseq_read(ks) >= 0:
		if ks.qual.l > 0: qual = ks.qual.s if isinstance(ks.qual.s, str) else ks.qual.s.decode()
		else: qual = None
		name = ks.name.s if isinstance(ks.name.s, str) else ks.name.s.decode()
		seq = ks.seq.s if isinstance(ks.seq.s, str) else ks.seq.s.decode()
		yield name, seq, qual
	cmappy.mm_fastx_close(ks)

def verbose(v=None):
	if v is None: v = -1
	return cmappy.mm_verbose_level(v)
