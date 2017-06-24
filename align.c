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

static void mm_update_extra(mm_extra_t *p, const uint8_t *qseq, const uint8_t *tseq, uint32_t n_cigar, uint32_t *cigar)
{
	uint32_t k, l, toff = 0, qoff = 0;
	for (k = 0; k < n_cigar; ++k) {
		uint32_t op = cigar[k]&0xf, len = cigar[k]>>4;
		if (op == 0) {
			for (l = 0; l < len; ++l) {
				if (tseq[toff + l] > 3 || qseq[qoff + l] > 3) ++p->n_ambi;
				else if (tseq[toff + l] != qseq[qoff + l]) ++p->n_diff;
			}
			toff += len, qoff += len, p->blen += len;
		} else if (op == 1) {
			int n_ambi = 0;
			for (l = 0; l < len; ++l)
				if (qseq[qoff + l] > 3) ++n_ambi;
			qoff += len, p->blen += len;
			p->n_ambi += n_ambi, p->n_diff += len - n_ambi;
		} else if (op == 2) {
			int n_ambi = 0;
			for (l = 0; l < len; ++l)
				if (tseq[toff + l] > 3) ++n_ambi;
			toff += len, p->blen += len;
			p->n_ambi += n_ambi, p->n_diff += len - n_ambi;
		}
	}
}

static void mm_append_cigar(mm_reg1_t *r, uint32_t n_cigar, uint32_t *cigar) // TODO: this calls the libc realloc()
{
	mm_extra_t *p;
	if (n_cigar == 0) return;
	if (r->p == 0) {
		uint32_t capacity = n_cigar + sizeof(mm_extra_t);
		kroundup32(capacity);
		r->p = (mm_extra_t*)calloc(capacity, 4);
		r->p->capacity = capacity;
	} else if (r->p->n_cigar + n_cigar + sizeof(mm_extra_t) > r->p->capacity) {
		r->p->capacity = r->p->n_cigar + n_cigar + sizeof(mm_extra_t);
		kroundup32(r->p->capacity);
		r->p = (mm_extra_t*)realloc(r->p, r->p->capacity * 4);
	}
	p = r->p;
	if (p->n_cigar > 0 && (p->cigar[p->n_cigar-1]&0xf) == (cigar[0]&0xf)) { // same CIGAR op at the boundary
		p->cigar[p->n_cigar-1] += cigar[0]>>4<<4;
		if (n_cigar > 1) memcpy(p->cigar + p->n_cigar, cigar + 1, (n_cigar - 1) * 4);
		p->n_cigar += n_cigar - 1;
	} else {
		memcpy(p->cigar + p->n_cigar, cigar, n_cigar * 4);
		p->n_cigar += n_cigar;
	}
}

static inline void mm_adjust_minier(const mm_idx_t *mi, uint8_t *qseq0[2], uint8_t *tseq0, mm128_t *a, int32_t *r, int32_t *q)
{
	if (mi->is_hpc && 1) {
		int span, rev, i, r0, c;
		rev = a->x >> 63;
		*q = (int32_t)a->y;
		for (i = *q - 1, c = qseq0[rev][*q]; i > 0; --i)
			if (qseq0[rev][i] != c) break;
		*q = i + 1;
		span = a->y >> 32;
		*r = (int32_t)a->x;
		r0 = *r + 1 >= 256? *r + 1 - 256 : 0;
		mm_idx_getseq(mi, a->x<<1>>33, r0, *r + 1, tseq0);
		for (i = *r - 1 - r0, c = tseq0[*r - r0]; i > 0; --i)
			if (tseq0[i] != c) break;
		*r = i + r0 + 1;
	} else {
		*r = (int32_t)a->x + 1;
		*q = (int32_t)a->y + 1;
	}
}

static void mm_align1(void *km, const mm_mapopt_t *opt, const mm_idx_t *mi, int qlen, uint8_t *qseq0[2], mm_reg1_t *r, mm_reg1_t *r2, mm128_t *a, ksw_extz_t *ez)
{
	int32_t rid = a[r->as].x<<1>>33, rev = a[r->as].x>>63;
	uint8_t *tseq, *qseq;
	int32_t i, l, bw, rs0, re0, qs0, qe0;
	int32_t rs, re, qs, qe;
	int32_t rs1, qs1, re1, qe1;
	int8_t mat[25];

	ksw_gen_simple_mat(5, mat, opt->a, opt->b);
	bw = (int)(opt->bw * 1.5 + 1.);

	/*
	rs = (int32_t)a[r->as].x + 1; // NB: this is the same as r->{rs,re}
	re = (int32_t)a[r->as + r->cnt - 1].x + 1;
	qs = (int32_t)a[r->as].y + 1; // NB: this is the coordinate on the reverse strand; r->{qs,qe} are on the reverse strand
	qe = (int32_t)a[r->as + r->cnt - 1].y + 1;
	*/
	tseq = (uint8_t*)kmalloc(km, 256);
	mm_adjust_minier(mi, qseq0, tseq, &a[r->as], &rs, &qs);
	mm_adjust_minier(mi, qseq0, tseq, &a[r->as + r->cnt - 1], &re, &qe);
	kfree(km, tseq);

	if (qs > 0 && rs > 0) {
		l = qs < opt->max_gap? qs : opt->max_gap;
		qs0 = qs - l;
		l += (l * opt->a - opt->q) / opt->e;
		l = l < opt->max_gap? l : opt->max_gap;
		l = l < rs? l : rs;
		rs0 = rs - l;
	} else rs0 = rs, qs0 = qs;

	if (qe < qlen && re < mi->seq[rid].len) {
		l = qlen - qe < opt->max_gap? qlen - qe : opt->max_gap;
		qe0 = qe + l;
		l += (l * opt->a - opt->q) / opt->e;
		l = l < opt->max_gap? l : opt->max_gap;
		l = l < mi->seq[rid].len - re? l : mi->seq[rid].len - re;
		re0 = re + l;
	} else re0 = re, qe0 = qe;

	tseq = (uint8_t*)kmalloc(km, re0 - rs0 > 256? re0 - rs0 : 256); // TODO: we can allocate a smaller size

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
		ksw_extz2_sse(km, ql, qseq, tl, tseq, 5, mat, opt->q, opt->e, bw, opt->zdrop, KSW_EZ_EXTZ_ONLY|KSW_EZ_RIGHT|KSW_EZ_REV_CIGAR, ez);
		mm_append_cigar(r, ez->n_cigar, ez->cigar);
		mm_update_extra(r->p, qseq, tseq, ez->n_cigar, ez->cigar);
		r->p->score += ez->max;
		rs1 = rs - (ez->max_t + 1);
		qs1 = qs - (ez->max_q + 1);
		mm_seq_rev(ql, qseq);
	} else rs1 = rs, qs1 = qs;
	assert(qs1 >= 0 && rs1 >= 0);

	for (i = 1; i < r->cnt; ++i) { // gap filling
		mm_adjust_minier(mi, qseq0, tseq, &a[r->as + i], &re, &qe);
		if (i == r->cnt - 1 || qe - qs >= opt->min_ksw_len || re - rs >= opt->min_ksw_len) {
			qseq = &qseq0[rev][qs];
			mm_idx_getseq(mi, rid, rs, re, tseq);
			#if 0
			fprintf(stderr, "===> [%d] %d-%d %c (%s:%d-%d) <===\n", i, qs, qe, "+-"[rev], mi->seq[rid].name, rs, re);
			for (k = 0; k < re - rs; ++k) fputc("ACGTN"[tseq[k]], stderr); fputc('\n', stderr);
			for (k = 0; k < qe - qs; ++k) fputc("ACGTN"[qseq[k]], stderr); fputc('\n', stderr);
			#endif
			ksw_extz2_sse(km, qe-qs, qseq, re-rs, tseq, 5, mat, opt->q, opt->e, bw, opt->zdrop, KSW_EZ_DYN_BAND, ez);
			mm_append_cigar(r, ez->n_cigar, ez->cigar);
			mm_update_extra(r->p, qseq, tseq, ez->n_cigar, ez->cigar);
			if (ez->score == KSW_NEG_INF) { // truncated by Z-drop
				r->p->score += ez->max;
				abort();
			} else {
				r->p->score += ez->score;
			}
			//for (k = 0; k < r->cigar->n; ++k) fprintf(stderr, "%d%c", r->cigar->cigar[k]>>4, "MID"[r->cigar->cigar[k]&0xf]); fputc('\n', stderr);
			rs = re, qs = qe;
		}
	}

	if (i == r->cnt && qe < qe0 && re < re0) { // right extension
		qseq = &qseq0[rev][qe];
		mm_idx_getseq(mi, rid, re, re0, tseq);
		ksw_extz2_sse(km, qe0-qe, qseq, re0-re, tseq, 5, mat, opt->q, opt->e, bw, opt->zdrop, KSW_EZ_EXTZ_ONLY, ez);
		mm_append_cigar(r, ez->n_cigar, ez->cigar);
		mm_update_extra(r->p, qseq, tseq, ez->n_cigar, ez->cigar);
		r->p->score += ez->max;
		re1 = re + (ez->max_t + 1);
		qe1 = qe + (ez->max_q + 1);
	} else re1 = re, qe1 = qe;
	assert(qe1 <= qlen);

	r->rs = rs1, r->re = re1;
	if (rev) r->qs = qlen - qe1, r->qe = qlen - qs1;
	else r->qs = qs1, r->qe = qe1;

//	for (i = 0; i < r->p->n_cigar; ++i) fprintf(stderr, "%d%c", r->p->cigar[i]>>4, "MID"[r->p->cigar[i]&0xf]); fputc('\n', stderr);
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
