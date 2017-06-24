#include <assert.h>
#include <string.h>
#include "minimap.h"
#include "ksw2.h"

#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

static void ksw_gen_simple_mat(int m, int8_t *mat, int8_t a, int8_t b)
{
	int i, j;
	a = a < 0? -a : a;
	b = b > 0? -b : b;
	for (i = 0; i < m - 1; ++i) {
		for (j = 0; j < m - 1; ++j)
			mat[i * m + j] = i == j? a : b;
		mat[i * m + m - 1] = 0;
	}
	for (j = 0; j < m; ++j)
		mat[(m - 1) * m + j] = 0;
}

static inline void mm_seq_rev(uint32_t len, uint8_t *seq)
{
	uint32_t i;
	uint8_t t;
	for (i = 0; i < len>>1; ++i)
		t = seq[i], seq[i] = seq[len - 1 - i], seq[len - 1 - i] = t;
}

static void mm_append_cigar(mm_reg1_t *r, uint32_t n_cigar, uint32_t *cigar) // TODO: this calls the libc realloc()
{
	if (n_cigar == 0) return;
	if (r->cigar == 0) {
		uint32_t m_cigar = n_cigar + 2;
		kroundup32(m_cigar);
		r->cigar = (mm_cigar_t*)malloc(m_cigar * 4);
		r->cigar->n = 0, r->cigar->m = m_cigar;
	} else if (r->cigar->n + n_cigar > r->cigar->m - 2) {
		r->cigar->m = r->cigar->n + n_cigar + 2;
		kroundup32(r->cigar->m);
		r->cigar = (mm_cigar_t*)realloc(r->cigar, r->cigar->m * 4);
	}
	if (r->cigar->n > 0 && (r->cigar->cigar[r->cigar->n-1]&0xf) == (cigar[0]&0xf)) { // same CIGAR op at the boundary
		r->cigar->cigar[r->cigar->n-1] += cigar[0]>>4<<4;
		if (n_cigar > 1) memcpy(r->cigar->cigar + r->cigar->n, cigar + 1, (n_cigar - 1) * 4);
		r->cigar->n += n_cigar - 1;
	} else {
		memcpy(r->cigar->cigar + r->cigar->n, cigar, n_cigar * 4);
		r->cigar->n += n_cigar;
	}
}

static void mm_align1(void *km, const mm_mapopt_t *opt, const mm_idx_t *mi, int qlen, uint8_t *qseq0[2], mm_reg1_t *r, mm_reg1_t *r2, mm128_t *a, ksw_extz_t *ez)
{
	int32_t rid = a[r->as].x<<1>>33, rev = a[r->as].x>>63;
	uint8_t *tseq, *qseq;
	int32_t i, k, l, rs0, re0, qs0, qe0;
	int32_t rs, re, qs, qe;
	int32_t rs1, qs1;
	int8_t mat[25];

	ksw_gen_simple_mat(5, mat, opt->a, opt->b);

	rs = (int32_t)a[r->as].x + 1; // NB: this is the same as r->{rs,re}
	re = (int32_t)a[r->as + r->cnt - 1].x + 1;
	qs = (int32_t)a[r->as].y + 1; // NB: this is the coordinate on the reverse strand; r->{qs,qe} are on the reverse strand
	qe = (int32_t)a[r->as + r->cnt - 1].y + 1;

	if (qs > 0 && rs > 0) {
		l = qs < opt->max_gap? qs : opt->max_gap;
		qs0 = qs - l;
		l += (l * opt->a - opt->q) / opt->e;
		l = l < opt->max_gap? l : opt->max_gap;
		l = l < rs? l : rs;
		rs0 = rs - l;
	} else rs0 = rs, qs0 = qs;

	if (qe < qlen && re < mi->seq[rid].len) {
		l = qlen - re < opt->max_gap? qlen - re : opt->max_gap;
		qe0 = qe + l;
		l += (l * opt->a - opt->q) / opt->e;
		l = l < opt->max_gap? l : opt->max_gap;
		l = l < mi->seq[rid].len - re? l : mi->seq[rid].len - re;
		re0 = re + l;
	} else re0 = re, qe0 = qe;

	tseq = (uint8_t*)kmalloc(km, re0 - rs0);

	if (qs > 0 && rs > 0) { // left extension
		uint32_t ql = qs - qs0, tl = rs - rs0;
		qseq = &qseq0[rev][qs0];
		mm_idx_getseq(mi, rid, rs0, rs, tseq);
		mm_seq_rev(ql, qseq);
		mm_seq_rev(tl, tseq);
		#if 0
		fprintf(stderr, "===> [-1] %d-%d %c (%s:%d-%d) <===\n", qs0, qs, "+-"[rev], mi->seq[rid].name, rs0, rs);
		for (k = 0; k < tl; ++k) fputc("ACGTN"[tseq[k]], stderr); fputc('\n', stderr);
		for (k = 0; k < ql; ++k) fputc("ACGTN"[qseq[k]], stderr); fputc('\n', stderr);
		#endif
		ksw_extz2_sse(km, ql, qseq, tl, tseq, 5, mat, opt->q, opt->e, (int)(opt->bw * 1.5 + .499), opt->zdrop, KSW_EZ_EXTZ_ONLY|KSW_EZ_RIGHT|KSW_EZ_REV_CIGAR, ez);
		mm_seq_rev(ql, qseq);
		mm_append_cigar(r, ez->n_cigar, ez->cigar);
		rs1 = rs - (ez->max_t + 1);
		qs1 = qs - (ez->max_q + 1);
	} else rs1 = rs, qs1 = qs;

	for (i = 1; i < r->cnt; ++i) { // gap filling
		re = (int32_t)a[r->as + i].x + 1;
		qe = (int32_t)a[r->as + i].y + 1;
		if (i == r->cnt - 1 || qe - qs >= opt->min_ksw_len || re - rs >= opt->min_ksw_len) {
			qseq = &qseq0[rev][qs];
			mm_idx_getseq(mi, rid, rs, re, tseq);
			#if 0
			fprintf(stderr, "===> [%d] %d-%d %c (%s:%d-%d) <===\n", i, qs, qe, "+-"[rev], mi->seq[rid].name, rs, re);
			for (k = 0; k < re - rs; ++k) fputc("ACGTN"[tseq[k]], stderr); fputc('\n', stderr);
			for (k = 0; k < qe - qs; ++k) fputc("ACGTN"[qseq[k]], stderr); fputc('\n', stderr);
			#endif
			ksw_extz2_sse(km, qe-qs, qseq, re-rs, tseq, 5, mat, opt->q, opt->e, (int)(opt->bw * 1.5 + .499), opt->zdrop, KSW_EZ_DYN_BAND, ez);
			if (ez->score == KSW_NEG_INF) { // truncated by Z-drop
				abort();
			} else {
				mm_append_cigar(r, ez->n_cigar, ez->cigar);
			}
			//for (k = 0; k < r->cigar->n; ++k) fprintf(stderr, "%d%c", r->cigar->cigar[k]>>4, "MID"[r->cigar->cigar[k]&0xf]); fputc('\n', stderr);
			rs = re, qs = qe;
		}
	}

	if (i == r->cnt) { // right extension
	}

	for (i = 0; i < r->cigar->n; ++i) fprintf(stderr, "%d%c", r->cigar->cigar[i]>>4, "MID"[r->cigar->cigar[i]&0xf]); fputc('\n', stderr);
	kfree(km, tseq);
}

void mm_align_skeleton(void *km, const mm_mapopt_t *opt, const mm_idx_t *mi, int qlen, const char *qstr, int n_regs, mm_reg1_t *regs, mm128_t *a)
{
	extern unsigned char seq_nt4_table[256];
	int i, reg;
	uint8_t *qseq0[2];
	ksw_extz_t ez;

	memset(&ez, 0, sizeof(ksw_extz_t));
	qseq0[0] = (uint8_t*)kmalloc(km, qlen);
	qseq0[1] = (uint8_t*)kmalloc(km, qlen);
	for (i = 0; i < qlen; ++i) {
		qseq0[0][i] = seq_nt4_table[(uint8_t)qstr[i]];
		qseq0[1][qlen - 1 - i] = qseq0[0][i] < 4? 3 - qseq0[0][i] : 4;
	}

	for (reg = 0; reg < n_regs; ++reg) {
		mm_reg1_t r2;
		mm_align1(km, opt, mi, qlen, qseq0, &regs[reg], &r2, a, &ez);
	}

	kfree(km, qseq0[0]); kfree(km, qseq0[1]); kfree(km, ez.cigar);
}
