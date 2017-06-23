#include <assert.h>
#include "minimap.h"
#include "ksw2.h"

static inline void mm_seq_rev(uint32_t len, uint8_t *seq)
{
	uint32_t i;
	uint8_t t;
	for (i = 0; i < len>>1; ++i)
		t = seq[i], seq[i] = seq[len - 1 - i], seq[len - 1 - i] = t;
}

static void mm_align1(void *km, const mm_mapopt_t *opt, const mm_idx_t *mi, int qlen, uint8_t *qseq0[2], mm_reg1_t *r, mm_reg1_t *r2, mm128_t *a)
{
	int32_t rid = a[r->as].x<<1>>33, rev = a[r->as].x>>63;
	uint8_t *tseq, *qseq;
	int32_t i, k, l, rs0, re0, qs0, qe0;
	int32_t rs, re, qs, qe, ret;
	mm_reg1_t r1;

	rs = (int32_t)a[r->as].x + 1; // NB: this is the same as r->{rs,re}
	re = (int32_t)a[r->as + r->cnt - 1].x + 1;
	qs = (int32_t)a[r->as].y + 1; // NB: this is the coordinate on the reverse strand; r->{qs,qe} are on the reverse strand
	qe = (int32_t)a[r->as + r->cnt - 1].y + 1;

	if (qs > 0 && rs > 0) {
		l = qs < opt->max_gap? qs : opt->max_gap;
		qs0 = qs - l;
		l = (l * opt->a - opt->q) / opt->e;
		l = l < opt->max_gap? l : opt->max_gap;
		l = l < rs? l : rs;
		rs0 = rs - l;
	} else rs0 = rs, qs0 = qs;

	if (qe < qlen && re < mi->seq[rid].len) {
		l = qlen - re < opt->max_gap? qlen - re : opt->max_gap;
		qe0 = qe + l;
		l = (l * opt->a - opt->q) / opt->e;
		l = l < opt->max_gap? l : opt->max_gap;
		l = l < mi->seq[rid].len - re? l : mi->seq[rid].len - re;
		re0 = re + l;
	} else re0 = re, qe0 = qe;

	tseq = (uint8_t*)kmalloc(km, re0 - rs0);

	if (qs > 0 && rs > 0) { // left extension
		uint32_t ql = qs - qs0, tl = rs - rs0;
		qseq = &qseq0[rev][qs0];
		ret = mm_idx_getseq(mi, rid, rs0, rs, tseq);
		assert(ret > 0);
		mm_seq_rev(ql, qseq);
		mm_seq_rev(tl, tseq);
		fprintf(stderr, "===> [-1] %d-%d %c (%s:%d-%d) <===\n", qs0, qs, "+-"[rev], mi->seq[rid].name, rs0, rs);
		for (k = 0; k < tl; ++k) fputc("ACGTN"[tseq[k]], stderr); fputc('\n', stderr);
		for (k = 0; k < ql; ++k) fputc("ACGTN"[qseq[k]], stderr); fputc('\n', stderr);
		mm_seq_rev(ql, qseq);
	}
/*
	for (i = 1; i < r->cnt; ++i) {
		re = (int32_t)a[r->as + i].x + 1;
		qe = (int32_t)a[r->as + i].y + 1;
		if (i == r->cnt - 1 || qe - qs >= opt->min_ksw_len || re - rs >= opt->min_ksw_len) {
			qseq = &qseq0[rev][qs];
			tseq = &tseq0[rs - rs0];
			#if 0
			fprintf(stderr, "===> [%d] %d-%d %c (%s:%d-%d) <===\n", i, qs, qe, "+-"[rev], mi->seq[rid].name, rs, re);
			for (k = 0; k < re - rs; ++k) fputc("ACGTN"[tseq[k]], stderr); fputc('\n', stderr);
			for (k = 0; k < qe - qs; ++k) fputc("ACGTN"[qseq[k]], stderr); fputc('\n', stderr);
			#endif
			rs = re, qs = qe;
		}
	}
*/
	kfree(km, tseq);
}

void mm_align_skeleton(void *km, const mm_mapopt_t *opt, const mm_idx_t *mi, int qlen, const char *qstr, int n_regs, mm_reg1_t *regs, mm128_t *a)
{
	extern unsigned char seq_nt4_table[256];
	int i, reg;
	uint8_t *qseq0[2];

	qseq0[0] = (uint8_t*)kmalloc(km, qlen);
	qseq0[1] = (uint8_t*)kmalloc(km, qlen);
	for (i = 0; i < qlen; ++i) {
		qseq0[0][i] = seq_nt4_table[(uint8_t)qstr[i]];
		qseq0[1][qlen - 1 - i] = qseq0[0][i] < 4? 3 - qseq0[0][i] : 4;
	}

	for (reg = 0; reg < n_regs; ++reg) {
		mm_reg1_t r2;
		mm_align1(km, opt, mi, qlen, qseq0, &regs[reg], &r2, a);
	}

	kfree(km, qseq0[0]); kfree(km, qseq0[1]);
}
