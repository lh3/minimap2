#ifndef CMAPPY_H
#define CMAPPY_H

#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include "minimap.h"
#include "kseq.h"
KSEQ_DECLARE(gzFile)

typedef struct {
	const char *ctg;
	int32_t ctg_start, ctg_end;
	int32_t qry_start, qry_end;
	int32_t blen, mlen, NM, ctg_len;
	uint8_t mapq, is_primary;
	int8_t strand, trans_strand;
	int32_t seg_id;
	int32_t n_cigar32;
	uint32_t *cigar32;
} mm_hitpy_t;

static inline void mm_reg2hitpy(const mm_idx_t *mi, mm_reg1_t *r, mm_hitpy_t *h)
{
	h->ctg = mi->seq[r->rid].name;
	h->ctg_len = mi->seq[r->rid].len;
	h->ctg_start = r->rs, h->ctg_end = r->re;
	h->qry_start = r->qs, h->qry_end = r->qe;
	h->strand = r->rev? -1 : 1;
	h->mapq = r->mapq;
	h->mlen = r->mlen;
	h->blen = r->blen;
	h->NM = r->blen - r->mlen + r->p->n_ambi;
	h->trans_strand = r->p->trans_strand == 1? 1 : r->p->trans_strand == 2? -1 : 0;
	h->is_primary = (r->id == r->parent);
	h->seg_id = r->seg_id;
	h->n_cigar32 = r->p->n_cigar;
	h->cigar32 = r->p->cigar;
}

static inline void mm_free_reg1(mm_reg1_t *r)
{
	free(r->p);
}

static inline kseq_t *mm_fastx_open(const char *fn)
{
	gzFile fp;
	fp = fn && strcmp(fn, "-") != 0? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
	return kseq_init(fp);
}

static inline void mm_fastx_close(kseq_t *ks)
{
	gzFile fp;
	fp = ks->f->f;
	kseq_destroy(ks);
	gzclose(fp);
}

static inline int mm_verbose_level(int v)
{
	if (v >= 0) mm_verbose = v;
	return mm_verbose;
}

static inline void mm_reset_timer(void)
{
	extern double realtime(void);
	mm_realtime0 = realtime();
}

extern unsigned char seq_comp_table[256];
static inline mm_reg1_t *mm_map_aux(const mm_idx_t *mi, const char *seq1, const char *seq2, int *n_regs, mm_tbuf_t *b, const mm_mapopt_t *opt)
{
	mm_reg1_t *r;

	Py_BEGIN_ALLOW_THREADS
	if (seq2 == 0) {
		r = mm_map(mi, strlen(seq1), seq1, n_regs, b, opt, NULL);
	} else {
		int _n_regs[2];
		mm_reg1_t *regs[2];
		char *seq[2];
		int i, len[2];

		len[0] = strlen(seq1);
		len[1] = strlen(seq2);
		seq[0] = (char*)seq1;
		seq[1] = strdup(seq2);
		for (i = 0; i < len[1]>>1; ++i) {
			int t = seq[1][len[1] - i - 1];
			seq[1][len[1] - i - 1] = seq_comp_table[(uint8_t)seq[1][i]];
			seq[1][i] = seq_comp_table[t];
		}
		if (len[1]&1) seq[1][len[1]>>1] = seq_comp_table[(uint8_t)seq[1][len[1]>>1]];
		mm_map_frag(mi, 2, len, (const char**)seq, _n_regs, regs, b, opt, NULL);
		for (i = 0; i < _n_regs[1]; ++i)
			regs[1][i].rev = !regs[1][i].rev;
		*n_regs = _n_regs[0] + _n_regs[1];
		regs[0] = (mm_reg1_t*)realloc(regs[0], sizeof(mm_reg1_t) * (*n_regs));
		memcpy(&regs[0][_n_regs[0]], regs[1], _n_regs[1] * sizeof(mm_reg1_t));
		free(regs[1]);
		r = regs[0];
	}
	Py_END_ALLOW_THREADS

	return r;
}

static inline char *mappy_revcomp(int len, const uint8_t *seq)
{
	int i;
	char *rev;
	rev = (char*)malloc(len + 1);
	for (i = 0; i < len; ++i)
		rev[len - i - 1] = seq_comp_table[seq[i]];
	rev[len] = 0;
	return rev;
}

static char *mappy_fetch_seq(const mm_idx_t *mi, const char *name, int st, int en, int *len)
{
	int i, rid;
	char *s;
	*len = 0;
	rid = mm_idx_name2id(mi, name);
	if (rid < 0) return 0;
	if ((uint32_t)st >= mi->seq[rid].len || st >= en) return 0;
	if (en < 0 || (uint32_t)en > mi->seq[rid].len)
		en = mi->seq[rid].len;
	s = (char*)malloc(en - st + 1);
	*len = mm_idx_getseq(mi, rid, st, en, (uint8_t*)s);
	for (i = 0; i < *len; ++i)
		s[i] = "ACGTN"[(uint8_t)s[i]];
	s[*len] = 0;
	return s;
}

static mm_idx_t *mappy_idx_seq(int w, int k, int is_hpc, int bucket_bits, const char *seq, int len)
{
	const char *fake_name = "N/A";
	char *s;
	mm_idx_t *mi;
	s = (char*)calloc(len + 1, 1);
	memcpy(s, seq, len);
	mi = mm_idx_str(w, k, is_hpc, bucket_bits, 1, (const char**)&s, (const char**)&fake_name);
	free(s);
	return mi;
}

#endif
