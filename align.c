#include <assert.h>
#include <string.h>
#include "minimap.h"
#include "mmpriv.h"
#include "ksw2.h"

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

static void mm_reg_split(mm_reg1_t *r, mm_reg1_t *r2, int n, int qlen, mm128_t *a)
{
	if (n <= 0 || n >= r->cnt) return;
	*r2 = *r;
	r2->p = 0;
	r2->cnt = r->cnt - n;
	r2->score = (int32_t)(r->score * ((float)r2->cnt / r->cnt) + .499);
	r2->as = r->as + n;
	mm_reg_set_coor(r2, qlen, a);
	r->cnt -= r2->cnt;
	r->score -= r2->score;
	mm_reg_set_coor(r, qlen, a);
	r->split = r2->split = 1;
}

static void mm_update_extra(mm_extra_t *p, const uint8_t *qseq, const uint8_t *tseq, uint32_t n_cigar, uint32_t *cigar, int rev_cigar)
{
	uint32_t l, toff = 0, qoff = 0;
	int32_t k, step, st, en;
	if (rev_cigar) step = -1, st = n_cigar - 1, en = -1;
	else step = 1, st = 0, en = n_cigar;
	for (k = st; k != en; k += step) {
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

static inline int mm_get_hplen_back(const mm_idx_t *mi, uint32_t rid, uint32_t x)
{
	int64_t i, off0 = mi->seq[rid].offset, off = off0 + x;
	int c = mm_seq4_get(mi->S, off);
	for (i = off - 1; i >= off0; --i)
		if (mm_seq4_get(mi->S, i) != c) break;
	return (int)(off - i);
}

static inline void mm_adjust_minier(const mm_idx_t *mi, uint8_t *const qseq0[2], mm128_t *a, int32_t *r, int32_t *q)
{
	if (mi->is_hpc) {
		const uint8_t *qseq = qseq0[a->x>>63];
		int i, c;
		*q = (int32_t)a->y;
		for (i = *q - 1, c = qseq[*q]; i > 0; --i)
			if (qseq[i] != c) break;
		*q = i + 1;
		c = mm_get_hplen_back(mi, a->x<<1>>33, (int32_t)a->x);
		*r = (int32_t)a->x + 1 - c;
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

	r2->cnt = 0;
	mm_adjust_minier(mi, qseq0, &a[r->as], &rs, &qs);
	mm_adjust_minier(mi, qseq0, &a[r->as + r->cnt - 1], &re, &qe);

	// compute rs0 and qs0
	if (r->split && r->as > 0) {
		mm_adjust_minier(mi, qseq0, &a[r->as-1], &rs0, &qs0);
	} else {
		if (qs > 0 && rs > 0) { // actually this is always true
			l = qs < opt->max_gap? qs : opt->max_gap;
			qs0 = qs - l;
			l += (l * opt->a - opt->q) / opt->e;
			l = l < opt->max_gap? l : opt->max_gap;
			l = l < rs? l : rs;
			rs0 = rs - l;
		} else rs0 = rs, qs0 = qs;

	}
	// compute re0 and qe0
	if (qe < qlen && re < mi->seq[rid].len) {
		l = qlen - qe < opt->max_gap? qlen - qe : opt->max_gap;
		qe0 = qe + l;
		l += (l * opt->a - opt->q) / opt->e;
		l = l < opt->max_gap? l : opt->max_gap;
		l = l < mi->seq[rid].len - re? l : mi->seq[rid].len - re;
		re0 = re + l;
	} else re0 = re, qe0 = qe;

	tseq = (uint8_t*)kmalloc(km, re0 - rs0); // TODO: we can allocate a smaller size

	if (qs > 0 && rs > 0) { // left extension
		qseq = &qseq0[rev][qs0];
		mm_idx_getseq(mi, rid, rs0, rs, tseq);
		mm_seq_rev(qs - qs0, qseq);
		mm_seq_rev(rs - rs0, tseq);
		ksw_extz2_sse(km, qs - qs0, qseq, rs - rs0, tseq, 5, mat, opt->q, opt->e, bw, opt->zdrop, KSW_EZ_EXTZ_ONLY|KSW_EZ_RIGHT|KSW_EZ_REV_CIGAR, ez);
		if (ez->n_cigar > 0) {
			mm_append_cigar(r, ez->n_cigar, ez->cigar);
			mm_update_extra(r->p, qseq, tseq, ez->n_cigar, ez->cigar, 1);
			r->p->score += ez->max;
		}
		rs1 = rs - (ez->max_t + 1);
		qs1 = qs - (ez->max_q + 1);
		mm_seq_rev(qs - qs0, qseq);
	} else rs1 = rs, qs1 = qs;
	re1 = rs, qe1 = qs;
	assert(qs1 >= 0 && rs1 >= 0);

	for (i = 1; i < r->cnt; ++i) { // gap filling
		mm_adjust_minier(mi, qseq0, &a[r->as + i], &re, &qe);
		re1 = re, qe1 = qe;
		if (i == r->cnt - 1 || qe - qs >= opt->min_ksw_len || re - rs >= opt->min_ksw_len) {
			qseq = &qseq0[rev][qs];
			mm_idx_getseq(mi, rid, rs, re, tseq);
			ksw_extz2_sse(km, qe - qs, qseq, re - rs, tseq, 5, mat, opt->q, opt->e, bw, opt->zdrop, 0, ez);
			if (ez->n_cigar > 0) {
				mm_append_cigar(r, ez->n_cigar, ez->cigar);
				mm_update_extra(r->p, qseq, tseq, ez->n_cigar, ez->cigar, 0);
			}
			if (ez->score == KSW_NEG_INF) { // truncated by Z-drop
				int j;
				for (j = i - 1; j >= 0; --j)
					if ((int32_t)a[r->as + j].x < re + ez->max_t)
						break;
				r->p->score += ez->max;
				re1 = rs + (ez->max_t + 1);
				qe1 = qs + (ez->max_q + 1);
				mm_reg_split(r, r2, j + 1, qlen, a);
				break;
			} else { // FIXME: in rare cases, r->p can be NULL, which leads to a segfault
				r->p->score += ez->score;
			}
			rs = re, qs = qe;
		}
	}

	if (r2->cnt == 0 && qe < qe0 && re < re0) { // right extension
		qseq = &qseq0[rev][qe];
		mm_idx_getseq(mi, rid, re, re0, tseq);
		ksw_extz2_sse(km, qe0 - qe, qseq, re0 - re, tseq, 5, mat, opt->q, opt->e, bw, opt->zdrop, KSW_EZ_EXTZ_ONLY, ez);
		if (ez->n_cigar > 0) {
			mm_append_cigar(r, ez->n_cigar, ez->cigar);
			mm_update_extra(r->p, qseq, tseq, ez->n_cigar, ez->cigar, 0);
			r->p->score += ez->max;
		}
		re1 = re + (ez->max_t + 1);
		qe1 = qe + (ez->max_q + 1);
	}
	assert(qe1 <= qlen);

	r->rs = rs1, r->re = re1;
	if (rev) r->qs = qlen - qe1, r->qe = qlen - qs1;
	else r->qs = qs1, r->qe = qe1;

	kfree(km, tseq);
}

mm_reg1_t *mm_align_skeleton(void *km, const mm_mapopt_t *opt, const mm_idx_t *mi, int qlen, const char *qstr, int *n_regs_, mm_reg1_t *regs, mm128_t *a)
{
	extern unsigned char seq_nt4_table[256];
	int i, r, n_regs = *n_regs_;
	uint8_t *qseq0[2];
	ksw_extz_t ez;

	memset(&ez, 0, sizeof(ksw_extz_t));
	qseq0[0] = (uint8_t*)kmalloc(km, qlen);
	qseq0[1] = (uint8_t*)kmalloc(km, qlen);
	for (i = 0; i < qlen; ++i) {
		qseq0[0][i] = seq_nt4_table[(uint8_t)qstr[i]];
		qseq0[1][qlen - 1 - i] = qseq0[0][i] < 4? 3 - qseq0[0][i] : 4;
	}

	for (r = 0; r < n_regs; ++r) {
		mm_reg1_t r2;
		mm_align1(km, opt, mi, qlen, qseq0, &regs[r], &r2, a, &ez);
		if (r2.cnt > 0) {
			regs = (mm_reg1_t*)realloc(regs, (n_regs + 1) * sizeof(mm_reg1_t)); // this should be very rare
			if (r + 1 != n_regs)
				memmove(&regs[r + 2], &regs[r + 1], sizeof(mm_reg1_t) * (n_regs - r - 1));
			regs[r + 1] = r2;
			++n_regs;
		}
	}
	*n_regs_ = n_regs;

	kfree(km, qseq0[0]); kfree(km, qseq0[1]); kfree(km, ez.cigar);
	return regs;
}
