#ifndef CMAPPY_H
#define CMAPPY_H

#include <stdlib.h>
#include "minimap.h"

typedef struct {
	const char *ctg;
	int32_t ctg_start, ctg_end;
	int32_t qry_start, qry_end;
	int32_t blen, NM, ctg_len;
	uint8_t mapq, is_primary;
	int8_t strand, trans_strand;
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
	h->blen = r->p->blen;
	h->NM = r->p->n_diff;
	h->trans_strand = r->p->trans_strand == 1? 1 : r->p->trans_strand == 2? -1 : 0;
	h->is_primary = (r->id == r->parent);
	h->n_cigar32 = r->p->n_cigar;
	h->cigar32 = r->p->cigar;
}

static inline void mm_free_reg1(mm_reg1_t *r)
{
	free(r->p);
}

#endif
