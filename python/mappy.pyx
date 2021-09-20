from libc.stdint cimport uint8_t, int8_t
from libc.stdlib cimport free
cimport cmappy
import sys

__version__ = '2.22'

cmappy.mm_reset_timer()

cdef class Alignment:
	cdef int _ctg_len, _r_st, _r_en
	cdef int _q_st, _q_en
	cdef int _NM, _mlen, _blen
	cdef int8_t _strand, _trans_strand
	cdef uint8_t _mapq, _is_primary
	cdef int _seg_id
	cdef _ctg, _cigar, _cs, _MD # these are python objects

	def __cinit__(self, ctg, cl, cs, ce, strand, qs, qe, mapq, cigar, is_primary, mlen, blen, NM, trans_strand, seg_id, cs_str, MD_str):
		self._ctg = ctg if isinstance(ctg, str) else ctg.decode()
		self._ctg_len, self._r_st, self._r_en = cl, cs, ce
		self._strand, self._q_st, self._q_en = strand, qs, qe
		self._NM, self._mlen, self._blen = NM, mlen, blen
		self._mapq = mapq
		self._cigar = cigar
		self._is_primary = is_primary
		self._trans_strand = trans_strand
		self._seg_id = seg_id
		self._cs = cs_str
		self._MD = MD_str

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
	def read_num(self): return self._seg_id + 1

	@property
	def cs(self): return self._cs

	@property
	def MD(self): return self._MD

	@property
	def cigar_str(self):
		return "".join(map(lambda x: str(x[0]) + 'MIDNSHP=XB'[x[1]], self._cigar))

	def __str__(self):
		if self._strand > 0: strand = '+'
		elif self._strand < 0: strand = '-'
		else: strand = '?'
		if self._is_primary != 0: tp = 'tp:A:P'
		else: tp = 'tp:A:S'
		if self._trans_strand > 0: ts = 'ts:A:+'
		elif self._trans_strand < 0: ts = 'ts:A:-'
		else: ts = 'ts:A:.'
		a = [str(self._q_st), str(self._q_en), strand, self._ctg, str(self._ctg_len), str(self._r_st), str(self._r_en),
			str(self._mlen), str(self._blen), str(self._mapq), tp, ts, "cg:Z:" + self.cigar_str]
		if self._cs != "": a.append("cs:Z:" + self._cs)
		return "\t".join(a)

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

	def __cinit__(self, fn_idx_in=None, preset=None, k=None, w=None, min_cnt=None, min_chain_score=None, min_dp_score=None, bw=None, best_n=None, n_threads=3, fn_idx_out=None, max_frag_len=None, extra_flags=None, seq=None, scoring=None):
		self._idx = NULL
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
		if best_n is not None: self.map_opt.best_n = best_n
		if max_frag_len is not None: self.map_opt.max_frag_len = max_frag_len
		if extra_flags is not None: self.map_opt.flag |= extra_flags
		if scoring is not None and len(scoring) >= 4:
			self.map_opt.a, self.map_opt.b = scoring[0], scoring[1]
			self.map_opt.q, self.map_opt.e = scoring[2], scoring[3]
			self.map_opt.q2, self.map_opt.e2 = self.map_opt.q, self.map_opt.e
			if len(scoring) >= 6:
				self.map_opt.q2, self.map_opt.e2 = scoring[4], scoring[5]
				if len(scoring) >= 7:
					self.map_opt.sc_ambi = scoring[6]

		cdef cmappy.mm_idx_reader_t *r;

		if seq is None:
			if fn_idx_out is None:
				r = cmappy.mm_idx_reader_open(str.encode(fn_idx_in), &self.idx_opt, NULL)
			else:
				r = cmappy.mm_idx_reader_open(str.encode(fn_idx_in), &self.idx_opt, str.encode(fn_idx_out))
			if r is not NULL:
				self._idx = cmappy.mm_idx_reader_read(r, n_threads) # NB: ONLY read the first part
				cmappy.mm_idx_reader_close(r)
				cmappy.mm_mapopt_update(&self.map_opt, self._idx)
				cmappy.mm_idx_index_name(self._idx)
		else:
			self._idx = cmappy.mappy_idx_seq(self.idx_opt.w, self.idx_opt.k, self.idx_opt.flag&1, self.idx_opt.bucket_bits, str.encode(seq), len(seq))
			cmappy.mm_mapopt_update(&self.map_opt, self._idx)
			self.map_opt.mid_occ = 1000 # don't filter high-occ seeds

	def __dealloc__(self):
		if self._idx is not NULL:
			cmappy.mm_idx_destroy(self._idx)

	def __bool__(self):
		return (self._idx != NULL)

	def map(self, seq, seq2=None, buf=None, cs=False, MD=False, max_frag_len=None, extra_flags=None):
		cdef cmappy.mm_reg1_t *regs
		cdef cmappy.mm_hitpy_t h
		cdef ThreadBuffer b
		cdef int n_regs
		cdef char *cs_str = NULL
		cdef int l_cs_str, m_cs_str = 0
		cdef void *km
		cdef cmappy.mm_mapopt_t map_opt

		if self._idx == NULL: return
		map_opt = self.map_opt
		if max_frag_len is not None: map_opt.max_frag_len = max_frag_len
		if extra_flags is not None: map_opt.flag |= extra_flags

		if self._idx is NULL: return None
		if buf is None: b = ThreadBuffer()
		else: b = buf
		km = cmappy.mm_tbuf_get_km(b._b)

		_seq = seq if isinstance(seq, bytes) else seq.encode()
		if seq2 is None:
			regs = cmappy.mm_map_aux(self._idx, _seq, NULL,  &n_regs, b._b, &map_opt)
		else:
			_seq2 = seq2 if isinstance(seq2, bytes) else seq2.encode()
			regs = cmappy.mm_map_aux(self._idx, _seq, _seq2, &n_regs, b._b, &map_opt)

		try:
			i = 0
			while i < n_regs:
				cmappy.mm_reg2hitpy(self._idx, &regs[i], &h)
				cigar, _cs, _MD = [], '', ''
				for k in range(h.n_cigar32): # convert the 32-bit CIGAR encoding to Python array
					c = h.cigar32[k]
					cigar.append([c>>4, c&0xf])
				if cs or MD: # generate the cs and/or the MD tag, if requested
					if cs:
						l_cs_str = cmappy.mm_gen_cs(km, &cs_str, &m_cs_str, self._idx, &regs[i], _seq, 1)
						_cs = cs_str[:l_cs_str] if isinstance(cs_str, str) else cs_str[:l_cs_str].decode()
					if MD:
						l_cs_str = cmappy.mm_gen_MD(km, &cs_str, &m_cs_str, self._idx, &regs[i], _seq)
						_MD = cs_str[:l_cs_str] if isinstance(cs_str, str) else cs_str[:l_cs_str].decode()
				yield Alignment(h.ctg, h.ctg_len, h.ctg_start, h.ctg_end, h.strand, h.qry_start, h.qry_end, h.mapq, cigar, h.is_primary, h.mlen, h.blen, h.NM, h.trans_strand, h.seg_id, _cs, _MD)
				cmappy.mm_free_reg1(&regs[i])
				i += 1
		finally:
			while i < n_regs:
				cmappy.mm_free_reg1(&regs[i])
				i += 1
			free(regs)
			free(cs_str)

	def seq(self, str name, int start=0, int end=0x7fffffff):
		cdef int l
		cdef char *s
		if self._idx == NULL: return
		s = cmappy.mappy_fetch_seq(self._idx, name.encode(), start, end, &l)
		if l == 0: return None
		r = s[:l] if isinstance(s, str) else s[:l].decode()
		free(s)
		return r

	@property
	def k(self): return self._idx.k

	@property
	def w(self): return self._idx.w

	@property
	def n_seq(self): return self._idx.n_seq

	@property
	def seq_names(self):
		cdef char *p
		if self._idx == NULL: return
		sn = []
		for i in range(self._idx.n_seq):
			p = self._idx.seq[i].name
			s = p if isinstance(p, str) else p.decode()
			sn.append(s)
		return sn

def fastx_read(fn, read_comment=False):
	cdef cmappy.kseq_t *ks
	ks = cmappy.mm_fastx_open(str.encode(fn))
	if ks is NULL: return None
	while cmappy.kseq_read(ks) >= 0:
		if ks.qual.l > 0: qual = ks.qual.s if isinstance(ks.qual.s, str) else ks.qual.s.decode()
		else: qual = None
		name = ks.name.s if isinstance(ks.name.s, str) else ks.name.s.decode()
		seq = ks.seq.s if isinstance(ks.seq.s, str) else ks.seq.s.decode()
		if read_comment:
			if ks.comment.l > 0: comment = ks.comment.s if isinstance(ks.comment.s, str) else ks.comment.s.decode()
			else: comment = None
			yield name, seq, qual, comment
		else:
			yield name, seq, qual
	cmappy.mm_fastx_close(ks)

def revcomp(seq):
	l = len(seq)
	bseq = seq if isinstance(seq, bytes) else seq.encode()
	cdef char *s = cmappy.mappy_revcomp(l, bseq)
	r = s[:l] if isinstance(s, str) else s[:l].decode()
	free(s)
	return r

def verbose(v=None):
	if v is None: v = -1
	return cmappy.mm_verbose_level(v)
