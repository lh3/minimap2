#include <assert.h>
#include "minimap.h"
#include "ksw2.h"

void mm_align_skeleton(void *km, const mm_mapopt_t *opt, const mm_idx_t *mi, int qlen, const char *qstr, int n_regs, mm_reg1_t *regs, mm128_t *a)
{
	extern unsigned char seq_nt4_table[256];
	int i, k, reg;
	uint8_t *qseq0[2];

	qseq0[0] = (uint8_t*)kmalloc(km, qlen);
	qseq0[1] = (uint8_t*)kmalloc(km, qlen);
	for (i = 0; i < qlen; ++i) {
		qseq0[0][i] = seq_nt4_table[(uint8_t)qstr[i]];
		qseq0[1][qlen - 1 - i] = qseq0[0][i] < 4? 3 - qseq0[0][i] : 4;
	}

	for (reg = 0; reg < n_regs; ++reg) {
		mm_reg1_t *r = &regs[reg];
		int32_t rid = a[r->as].x<<1>>33, rev = a[r->as].x>>63;
		uint8_t *tseq0, *tseq, *qseq;
		int32_t rs, re, qs, qe, ret;
		tseq0 = (uint8_t*)kmalloc(km, r->re - r->rs);
		ret = mm_idx_getseq(mi, rid, r->rs, r->re, tseq0);
		assert(ret > 0);
		rs = (int32_t)a[r->as].x + 1;
		qs = (int32_t)a[r->as].y + 1;
		for (i = 1; i < r->cnt; ++i) {
			re = (int32_t)a[r->as + i].x + 1;
			qe = (int32_t)a[r->as + i].y + 1;
			if (qe - qs >= opt->min_ksw_len || re - rs >= opt->min_ksw_len) {
				qseq = &qseq0[rev][qs];
				tseq = &tseq0[rs - r->rs];
				fprintf(stderr, "===> [%d,%d] %d-%d %c (%s:%d-%d) <===\n", reg, i, qs, qe, "+-"[rev], mi->seq[rid].name, rs, re);
				for (k = 0; k < re - rs; ++k) fputc("ACGTN"[tseq[k]], stderr); fputc('\n', stderr);
				for (k = 0; k < qe - qs; ++k) fputc("ACGTN"[qseq[k]], stderr); fputc('\n', stderr);
				rs = re, qs = qe;
			}
		}
		kfree(km, tseq0);
	}

	kfree(km, qseq0[0]); kfree(km, qseq0[1]);
}
