from libc.stdint cimport uint8_t, int8_t
from libc.stdlib cimport free
cimport cminimap2

cdef class Alignment:
	cdef const char *_ctg
	cdef int _r_st, _r_en
	cdef int _q_st, _q_en
	cdef int8_t _is_rev, _is_primary
	cdef uint8_t _mapq
	cdef _cigar

	def __cinit__(self, ctg, cs, ce, strand, qs, qe, mapq, cigar, is_primary):
		self._ctg, self._r_st, self._r_en = ctg, cs, ce
		self._is_rev, self._q_st, self._q_en = strand, qs, qe
		self._mapq = mapq
		self._cigar = cigar
		self._is_primary = is_primary

	@property
	def ctg(self):
		return self._ctg

	@property
	def r_st(self):
		return self._r_st

	@property
	def r_en(self):
		return self._r_en

	@property
	def is_rev(self):
		return self._is_rev

	@property
	def is_primary(self):
		return self._is_primary

	@property
	def q_st(self):
		return self._q_st

	@property
	def q_en(self):
		return self._q_en

	@property
	def mapq(self):
		return self._mapq

	@property
	def cigar(self):
		return self._cigar

cdef class ThreadBuffer:
	cdef cminimap2.mm_tbuf_t *_b

	def __cinit__(self):
		self._b = cminimap2.mm_tbuf_init()

	def __dealloc__(self):
		cminimap2.mm_tbuf_destroy(self._b)

cdef class Aligner:
	cdef cminimap2.mm_idx_t *_idx
	cdef public cminimap2.mm_idxopt_t idx_opt
	cdef public cminimap2.mm_mapopt_t map_opt

	def __cinit__(self, fn, preset=None):
		self.config(preset)
		cdef cminimap2.mm_idx_reader_t *r;
		r = cminimap2.mm_idx_reader_open(fn, &self.idx_opt, NULL)
		self._idx = cminimap2.mm_idx_reader_read(r, 3) # NB: ONLY read the first part
		cminimap2.mm_idx_reader_close(r)
		cminimap2.mm_mapopt_update(&self.map_opt, self._idx)
	
	def __dealloc__(self):
		if self._idx is not NULL:
			cminimap2.mm_idx_destroy(self._idx)

	def config(self, preset=None):
		cminimap2.mm_set_opt(NULL, &self.idx_opt, &self.map_opt)
		if preset is not None:
			cminimap2.mm_set_opt(preset, &self.idx_opt, &self.map_opt)
		self.map_opt.flag |= 4 # always perform alignment

	def map(self, seq, buf=None):
		cdef cminimap2.mm_reg1_t *regs
		cdef ThreadBuffer b
		cdef int n_regs

		if buf is None: b = ThreadBuffer()
		else: b = buf
		regs = cminimap2.mm_map(self._idx, len(seq), seq, &n_regs, b._b, &self.map_opt, NULL)
		hits = cminimap2.mm_reg2hitpy(self._idx, n_regs, regs)

		arr = []
		for i in range(n_regs):
			h = hits[i]
			cigar = []
			for k in range(h.n_cigar32):
				c = h.cigar32[k]
				cigar.append([c>>4, c&0xf])
			arr.append(Alignment(h.ctg, h.ctg_start, h.ctg_end, h.strand, h.qry_start, h.qry_end, h.mapq, cigar, h.is_primary))
			cminimap2.mm_free_reg1(&regs[i])
		free(regs)
		free(hits)

		return arr
