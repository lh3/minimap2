#ifndef CMINIMAP2_H
#define CMINIMAP2_H

#include <stdlib.h>
#include "minimap.h"

typedef struct {
	const char *ctg;
	int32_t ctg_start, ctg_end;
	int32_t qry_start, qry_end;
	int32_t n_cigar32;
	uint8_t strand, mapq, is_primary;
	uint32_t *cigar32;
} mm_hitpy_t;

static inline mm_hitpy_t *mm_reg2hitpy(const mm_idx_t *mi, int n_regs, mm_reg1_t *regs)
{
	mm_hitpy_t *h;
	int i;
	h = (mm_hitpy_t*)calloc(n_regs, sizeof(mm_hitpy_t));
	for (i = 0; i < n_regs; ++i) {
		mm_reg1_t *r = &regs[i];
		h[i].ctg = mi->seq[r->rid].name;
		h[i].ctg_start = r->rs, h[i].ctg_end = r->re;
		h[i].qry_start = r->qs, h[i].qry_end = r->qe;
		h[i].strand = r->rev;
		h[i].mapq = r->mapq;
		h[i].is_primary = (r->id == r->parent);
		h[i].n_cigar32 = r->p->n_cigar;
		h[i].cigar32 = r->p->cigar;
	}
	return h;
}

static inline void mm_free_reg1(mm_reg1_t *r)
{
	free(r->p);
}

#endif
