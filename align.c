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

static inline int test_zdrop_aux(int32_t score, int i, int j, int32_t *max, int *max_i, int *max_j, int e, int zdrop)
{
	if (score < *max) {
		int li = i - *max_i;
		int lj = j - *max_j;
		int diff = li > lj? li - lj : lj - li;
		if (*max - score > zdrop + diff * e)
			return 1;
	} else *max = score, *max_i = i, *max_j = j;
	return 0;
}

static int mm_check_zdrop(const uint8_t *qseq, const uint8_t *tseq, uint32_t n_cigar, uint32_t *cigar, const int8_t *mat, int8_t q, int8_t e, int zdrop)
{
	uint32_t k;
	int32_t score = 0, max = 0, max_i = -1, max_j = -1, i = 0, j = 0;
	for (k = 0; k < n_cigar; ++k) {
		uint32_t l, op = cigar[k]&0xf, len = cigar[k]>>4;
		if (op == 0) {
			for (l = 0; l < len; ++l) {
				score += mat[tseq[i + l] * 5 + qseq[j + l]];
				if (test_zdrop_aux(score, i+l, j+l, &max, &max_i, &max_j, e, zdrop)) return 1;
			}
			i += len, j += len;
		} else if (op == 1) {
			score -= q + e * len, j += len;
			if (test_zdrop_aux(score, i, j, &max, &max_i, &max_j, e, zdrop)) return 1;
		} else if (op == 2 || op == 3) {
			score -= q + e * len, i += len;
			if (test_zdrop_aux(score, i, j, &max, &max_i, &max_j, e, zdrop)) return 1;
		}
	}
	return 0;
}

static void mm_update_extra(mm_extra_t *p, const uint8_t *qseq, const uint8_t *tseq, const int8_t *mat, int8_t q, int8_t e)
{
	uint32_t k, l, toff = 0, qoff = 0;
	int32_t s = 0, max = 0, n_gtag = 0, n_ctac = 0;
	if (p == 0) return;
	for (k = 0; k < p->n_cigar; ++k) {
		uint32_t op = p->cigar[k]&0xf, len = p->cigar[k]>>4;
		if (op == 0) { // match/mismatch
			for (l = 0; l < len; ++l) {
				int cq = qseq[qoff + l], ct = tseq[toff + l];
				if (ct > 3 || cq > 3) ++p->n_ambi;
				else if (ct != cq) ++p->n_diff;
				s += mat[ct * 5 + cq];
				if (s < 0) s = 0;
				else max = max > s? max : s;
			}
			toff += len, qoff += len, p->blen += len;
		} else if (op == 1) { // insertion
			int n_ambi = 0;
			for (l = 0; l < len; ++l)
				if (qseq[qoff + l] > 3) ++n_ambi;
			qoff += len, p->blen += len;
			p->n_ambi += n_ambi, p->n_diff += len - n_ambi;
			s -= q + e * len;
			if (s < 0) s = 0;
		} else if (op == 2) { // deletion
			int n_ambi = 0;
			for (l = 0; l < len; ++l)
				if (tseq[toff + l] > 3) ++n_ambi;
			toff += len, p->blen += len;
			p->n_ambi += n_ambi, p->n_diff += len - n_ambi;
			s -= q + e * len;
			if (s < 0) s = 0;
		} else if (op == 3) { // intron
			uint8_t b[4];
			b[0] = tseq[toff], b[1] = tseq[toff+1];
			b[2] = tseq[toff+len-2], b[3] = tseq[toff+len-1];
			if (memcmp(b, "\2\3\0\2", 4) == 0) ++n_gtag;
			else if (memcmp(b, "\1\3\0\1", 4) == 0) ++n_ctac;
			toff += len, p->blen += len;
		}
	}
	p->dp_max = max;
	if (n_gtag > n_ctac) p->trans_strand = 1;
	else if (n_gtag < n_ctac) p->trans_strand = 2;
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

static void mm_align_pair(void *km, const mm_mapopt_t *opt, int qlen, const uint8_t *qseq, int tlen, const uint8_t *tseq, const int8_t *mat, int w, int flag, ksw_extz_t *ez)
{
	if (mm_dbg_flag & MM_DBG_PRINT_ALN_SEQ) {
		int i;
		fprintf(stderr, "===> q=(%d,%d), e=(%d,%d), bw=%d, flag=%d, zdrop=%d <===\n", opt->q, opt->q2, opt->e, opt->e2, w, flag, opt->zdrop);
		for (i = 0; i < tlen; ++i) fputc("ACGTN"[tseq[i]], stderr);
		fputc('\n', stderr);
		for (i = 0; i < qlen; ++i) fputc("ACGTN"[qseq[i]], stderr);
		fputc('\n', stderr);
	}
	if (opt->flag & MM_F_APPROX_EXT) {
		flag |= KSW_EZ_APPROX_MAX;
		if (flag & KSW_EZ_EXTZ_ONLY) flag |= KSW_EZ_APPROX_DROP;
	}
	if (opt->flag & MM_F_SPLICE)
		ksw_exts2_sse(km, qlen, qseq, tlen, tseq, 5, mat, opt->q, opt->e, opt->q2, opt->noncan, opt->zdrop, flag, ez);
	else if (opt->q == opt->q2 && opt->e == opt->e2)
		ksw_extz2_sse(km, qlen, qseq, tlen, tseq, 5, mat, opt->q, opt->e, w, opt->zdrop, flag, ez);
	else
		ksw_extd2_sse(km, qlen, qseq, tlen, tseq, 5, mat, opt->q, opt->e, opt->q2, opt->e2, w, opt->zdrop, flag, ez);
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
		*r = (int32_t)a->x - (mi->k>>1);
		*q = (int32_t)a->y - (mi->k>>1);
	}
}

static void mm_filter_bad_seeds(void *km, int as1, int cnt1, mm128_t *a, int min_gap, int diff_thres, int max_ext_len, int max_ext_cnt)
{
	int max_st, max_en, n, i, k, max, *K;
	for (i = 1, n = 0; i < cnt1; ++i) { // count the number of gaps longer than min_gap
		int gap = ((int32_t)a[as1 + i].y - a[as1 + i - 1].y) - ((int32_t)a[as1 + i].x - a[as1 + i - 1].x);
		if (gap < -min_gap || gap > min_gap) ++n;
	}
	if (n <= 1) return;
	K = (int*)kmalloc(km, n * sizeof(int));
	for (i = 1, n = 0; i < cnt1; ++i) { // store the positions of long gaps
		int gap = ((int32_t)a[as1 + i].y - a[as1 + i - 1].y) - ((int32_t)a[as1 + i].x - a[as1 + i - 1].x);
		if (gap < -min_gap || gap > min_gap)
			K[n++] = i;
	}
	max = 0, max_st = max_en = -1;
	for (k = 0;; ++k) { // traverse long gaps
		int gap, l, n_ins = 0, n_del = 0, qs, rs, max_diff = 0, max_diff_l = -1;
		if (k == n || k >= max_en) {
			if (max_en > 0)
				for (i = K[max_st]; i < K[max_en]; ++i)
					a[as1 + i].y |= MM_SEED_IGNORE;
			max = 0, max_st = max_en = -1;
			if (k == n) break;
		}
		i = K[k];
		gap = ((int32_t)a[as1 + i].y - a[as1 + i - 1].y) - ((int32_t)a[as1 + i].x - a[as1 + i - 1].x);
		if (gap > 0) n_ins += gap;
		else n_del += -gap;
		qs = (int32_t)a[as1 + i - 1].y;
		rs = (int32_t)a[as1 + i - 1].x;
		for (l = k + 1; l < n && l <= k + max_ext_cnt; ++l) {
			int j = K[l], diff;
			if ((int32_t)a[as1 + j].y - qs > max_ext_len || (int32_t)a[as1 + j].x - rs > max_ext_len) break;
			gap = ((int32_t)a[as1 + j].y - (int32_t)a[as1 + j - 1].y) - (a[as1 + j].x - a[as1 + j - 1].x);
			if (gap > 0) n_ins += gap;
			else n_del += -gap;
			diff = n_ins + n_del - abs(n_ins - n_del);
			if (max_diff < diff)
				max_diff = diff, max_diff_l = l;
		}
		if (max_diff > diff_thres && max_diff > max)
			max = max_diff, max_st = k, max_en = max_diff_l;
	}
	kfree(km, K);
}

static void mm_fix_bad_ends(const mm_reg1_t *r, const mm128_t *a, int bw, int32_t *as, int32_t *cnt)
{
	int32_t i, l;
	*as = r->as, *cnt = r->cnt;
	if (r->cnt < 3) return;
	l = a[r->as].y >> 32 & 0xff;
	for (i = r->as + 1; i < r->as + r->cnt - 1; ++i) {
		int32_t lq, lr, min, max;
		lr = (int32_t)a[i].x - (int32_t)a[i-1].x;
		lq = (int32_t)a[i].y - (int32_t)a[i-1].y;
		min = lr < lq? lr : lq;
		max = lr > lq? lr : lq;
		if (max - min > l >> 1) *as = i;
		l += min;
		if (l >= bw << 1) break;
	}
	*cnt = r->as + r->cnt - *as;
	l = a[r->as + r->cnt - 1].y >> 32 & 0xff;
	for (i = r->as + r->cnt - 2; i > *as; --i) {
		int32_t lq, lr, min, max;
		lr = (int32_t)a[i+1].x - (int32_t)a[i].x;
		lq = (int32_t)a[i+1].y - (int32_t)a[i].y;
		min = lr < lq? lr : lq;
		max = lr > lq? lr : lq;
		if (max - min > l >> 1) *cnt = i + 1 - *as;
		l += min;
		if (l >= bw) break;
	}
}

static void mm_align1(void *km, const mm_mapopt_t *opt, const mm_idx_t *mi, int qlen, uint8_t *qseq0[2], mm_reg1_t *r, mm_reg1_t *r2, mm128_t *a, ksw_extz_t *ez, int splice_flag)
{
	int32_t rid = a[r->as].x<<1>>33, rev = a[r->as].x>>63, as1, cnt1;
	uint8_t *tseq, *qseq;
	int32_t i, l, bw, dropped = 0, extra_flag = 0, rs0, re0, qs0, qe0;
	int32_t rs, re, qs, qe;
	int32_t rs1, qs1, re1, qe1;
	int8_t mat[25];

	if (r->cnt == 0) return;
	ksw_gen_simple_mat(5, mat, opt->a, opt->b);
	bw = (int)(opt->bw * 1.5 + 1.);

	r2->cnt = 0;
	if (!(opt->flag & MM_F_SPLICE))
		mm_fix_bad_ends(r, a, opt->bw, &as1, &cnt1);
	else as1 = r->as, cnt1 = r->cnt;
	mm_filter_bad_seeds(km, as1, cnt1, a, 10, 40, opt->max_gap>>1, 10);
	mm_adjust_minier(mi, qseq0, &a[as1], &rs, &qs);
	mm_adjust_minier(mi, qseq0, &a[as1 + cnt1 - 1], &re, &qe);

	if (opt->flag & MM_F_SPLICE) {
		if (splice_flag & MM_F_SPLICE_FOR) extra_flag |= rev? KSW_EZ_SPLICE_REV : KSW_EZ_SPLICE_FOR;
		if (splice_flag & MM_F_SPLICE_REV) extra_flag |= rev? KSW_EZ_SPLICE_FOR : KSW_EZ_SPLICE_REV;
		if (splice_flag & MM_F_SPLICE_BOTH) extra_flag |= KSW_EZ_SPLICE_FOR|KSW_EZ_SPLICE_REV;
	}

	// compute rs0 and qs0
	if (r->split && as1 > 0) {
		mm_adjust_minier(mi, qseq0, &a[as1-1], &rs0, &qs0);
	} else {
		if (qs > 0 && rs > 0) { // actually this is always true
			l = qs < opt->max_gap? qs : opt->max_gap;
			qs0 = qs - l;
			l += l * opt->a > opt->q? (l * opt->a - opt->q) / opt->e : 0;
			l = l < opt->max_gap? l : opt->max_gap;
			l = l < rs? l : rs;
			rs0 = rs - l;
		} else rs0 = rs, qs0 = qs;

	}
	// compute re0 and qe0
	if (qe < qlen && re < mi->seq[rid].len) {
		l = qlen - qe < opt->max_gap? qlen - qe : opt->max_gap;
		qe0 = qe + l;
		l += l * opt->a > opt->q? (l * opt->a - opt->q) / opt->e : 0;
		l = l < opt->max_gap? l : opt->max_gap;
		l = l < mi->seq[rid].len - re? l : mi->seq[rid].len - re;
		re0 = re + l;
	} else re0 = re, qe0 = qe;

	assert(re0 > rs0);
	tseq = (uint8_t*)kmalloc(km, re0 - rs0);

	if (qs > 0 && rs > 0) { // left extension
		qseq = &qseq0[rev][qs0];
		mm_idx_getseq(mi, rid, rs0, rs, tseq);
		mm_seq_rev(qs - qs0, qseq);
		mm_seq_rev(rs - rs0, tseq);
		mm_align_pair(km, opt, qs - qs0, qseq, rs - rs0, tseq, mat, bw, extra_flag|KSW_EZ_EXTZ_ONLY|KSW_EZ_RIGHT|KSW_EZ_REV_CIGAR, ez);
		if (ez->n_cigar > 0) {
			mm_append_cigar(r, ez->n_cigar, ez->cigar);
			r->p->dp_score += ez->max;
		}
		rs1 = rs - (ez->max_t + 1);
		qs1 = qs - (ez->max_q + 1);
		mm_seq_rev(qs - qs0, qseq);
	} else rs1 = rs, qs1 = qs;
	re1 = rs, qe1 = qs;
	assert(qs1 >= 0 && rs1 >= 0);

	for (i = 1; i < cnt1; ++i) { // gap filling
		if ((a[as1+i].y & (MM_SEED_IGNORE|MM_SEED_TANDEM)) && i != cnt1 - 1) continue;
		mm_adjust_minier(mi, qseq0, &a[as1 + i], &re, &qe);
		re1 = re, qe1 = qe;
		if (i == cnt1 - 1 || (a[as1+i].y&MM_SEED_LONG_JOIN) || (qe - qs >= opt->min_ksw_len && re - rs >= opt->min_ksw_len)) {
			int bw1 = bw;
			if (a[as1+i].y & MM_SEED_LONG_JOIN)
				bw1 = qe - qs > re - rs? qe - qs : re - rs;
			qseq = &qseq0[rev][qs];
			mm_idx_getseq(mi, rid, rs, re, tseq);
			mm_align_pair(km, opt, qe - qs, qseq, re - rs, tseq, mat, bw1, extra_flag|KSW_EZ_APPROX_MAX, ez);
			if (mm_check_zdrop(qseq, tseq, ez->n_cigar, ez->cigar, mat, opt->q, opt->e, opt->zdrop))
				mm_align_pair(km, opt, qe - qs, qseq, re - rs, tseq, mat, bw1, extra_flag, ez);
			if (ez->n_cigar > 0)
				mm_append_cigar(r, ez->n_cigar, ez->cigar);
			if (ez->zdropped) { // truncated by Z-drop; TODO: sometimes Z-drop kicks in because the next seed placement is wrong. This can be fixed in principle.
				int j;
				for (j = i - 1; j >= 0; --j)
					if ((int32_t)a[as1 + j].x < re + ez->max_t)
						break;
				dropped = 1;
				r->p->dp_score += ez->max;
				re1 = rs + (ez->max_t + 1);
				qe1 = qs + (ez->max_q + 1);
				if (cnt1 - (j + 1) >= opt->min_cnt)
					mm_split_reg(r, r2, j + 1, qlen, a);
				break;
			} else r->p->dp_score += ez->score;
			rs = re, qs = qe;
		}
	}

	if (!dropped && qe < qe0 && re < re0) { // right extension
		qseq = &qseq0[rev][qe];
		mm_idx_getseq(mi, rid, re, re0, tseq);
		mm_align_pair(km, opt, qe0 - qe, qseq, re0 - re, tseq, mat, bw, extra_flag|KSW_EZ_EXTZ_ONLY, ez);
		if (ez->n_cigar > 0) {
			mm_append_cigar(r, ez->n_cigar, ez->cigar);
			r->p->dp_score += ez->max;
		}
		re1 = re + (ez->max_t + 1);
		qe1 = qe + (ez->max_q + 1);
	}
	assert(qe1 <= qlen);

	r->rs = rs1, r->re = re1;
	if (rev) r->qs = qlen - qe1, r->qe = qlen - qs1;
	else r->qs = qs1, r->qe = qe1;

	assert(re1 - rs1 <= re0 - rs0);
	if (r->p) {
		mm_idx_getseq(mi, rid, rs1, re1, tseq);
		mm_update_extra(r->p, &qseq0[r->rev][qs1], tseq, mat, opt->q, opt->e);
		if (rev && r->p->trans_strand)
			r->p->trans_strand ^= 3; // flip to the read strand
	}

	kfree(km, tseq);
}

static int mm_align1_inv(void *km, const mm_mapopt_t *opt, const mm_idx_t *mi, int qlen, uint8_t *qseq0[2], const mm_reg1_t *r1, const mm_reg1_t *r2, mm_reg1_t *r_inv, ksw_extz_t *ez)
{
	int tl, ql, score, ret = 0, q_off, t_off;
	uint8_t *tseq, *qseq;
	int8_t mat[25];
	void *qp;

	memset(r_inv, 0, sizeof(mm_reg1_t));
	if (!(r1->split&1) || !(r2->split&2)) return 0;
	if (r1->id != r1->parent && r1->parent != MM_PARENT_TMP_PRI) return 0;
	if (r2->id != r2->parent && r2->parent != MM_PARENT_TMP_PRI) return 0;
	if (r1->rid != r2->rid || r1->rev != r2->rev) return 0;
	ql = r2->qs - r1->qe;
	tl = r2->rs - r1->re;
	if (ql < opt->min_chain_score || ql > opt->max_gap) return 0;
	if (tl < opt->min_chain_score || tl > opt->max_gap) return 0;

	ksw_gen_simple_mat(5, mat, opt->a, opt->b);
	tseq = (uint8_t*)kmalloc(km, tl);
	mm_idx_getseq(mi, r1->rid, r1->re, r2->rs, tseq);
	qseq = &qseq0[!r1->rev][qlen - r2->qs];

	mm_seq_rev(ql, qseq);
	mm_seq_rev(tl, tseq);
	qp = ksw_ll_qinit(km, 2, ql, qseq, 5, mat);
	score = ksw_ll_i16(qp, tl, tseq, opt->q, opt->e, &q_off, &t_off);
	kfree(km, qp);
	mm_seq_rev(ql, qseq);
	mm_seq_rev(tl, tseq);
	if (score < opt->min_dp_max) goto end_align1_inv;
	q_off = ql - (q_off + 1), t_off = tl - (t_off + 1);
	mm_align_pair(km, opt, ql - q_off, qseq + q_off, tl - t_off, tseq + t_off, mat, (int)(opt->bw * 1.5), KSW_EZ_EXTZ_ONLY, ez);
	if (ez->n_cigar == 0) goto end_align1_inv; // should never be here
	mm_append_cigar(r_inv, ez->n_cigar, ez->cigar);
	r_inv->p->dp_score = ez->max;
	mm_update_extra(r_inv->p, qseq + q_off, tseq + t_off, mat, opt->q, opt->e);
	r_inv->id = -1;
	r_inv->parent = MM_PARENT_UNSET;
	r_inv->inv = 1;
	r_inv->rev = !r1->rev;
	r_inv->qs = r1->qe + q_off, r_inv->qe = r_inv->qs + ez->max_q + 1;
	r_inv->rs = r1->re + t_off, r_inv->re = r_inv->rs + ez->max_t + 1;
	ret = 1;
end_align1_inv:
	kfree(km, tseq);
	return ret;
}

static inline mm_reg1_t *mm_insert_reg(const mm_reg1_t *r, int i, int *n_regs, mm_reg1_t *regs)
{
	regs = (mm_reg1_t*)realloc(regs, (*n_regs + 1) * sizeof(mm_reg1_t));
	if (i + 1 != *n_regs)
		memmove(&regs[i + 2], &regs[i + 1], sizeof(mm_reg1_t) * (*n_regs - i - 1));
	regs[i + 1] = *r;
	++*n_regs;
	return regs;
}

mm_reg1_t *mm_align_skeleton(void *km, const mm_mapopt_t *opt, const mm_idx_t *mi, int qlen, const char *qstr, int *n_regs_, mm_reg1_t *regs, mm128_t *a)
{
	extern unsigned char seq_nt4_table[256];
	int32_t i, n_regs = *n_regs_;
	uint8_t *qseq0[2];
	ksw_extz_t ez;

	// encode the query sequence
	qseq0[0] = (uint8_t*)kmalloc(km, qlen);
	qseq0[1] = (uint8_t*)kmalloc(km, qlen);
	for (i = 0; i < qlen; ++i) {
		qseq0[0][i] = seq_nt4_table[(uint8_t)qstr[i]];
		qseq0[1][qlen - 1 - i] = qseq0[0][i] < 4? 3 - qseq0[0][i] : 4;
	}

	// align through seed hits
	memset(&ez, 0, sizeof(ksw_extz_t));
	for (i = 0; i < n_regs; ++i) {
		mm_reg1_t r2;
		if ((opt->flag&MM_F_SPLICE) && (opt->flag&MM_F_SPLICE_FOR) && (opt->flag&MM_F_SPLICE_REV)) {
			mm_reg1_t s[2], s2[2];
			int which, trans_strand;
			s[0] = s[1] = regs[i];
			mm_align1(km, opt, mi, qlen, qseq0, &s[0], &s2[0], a, &ez, MM_F_SPLICE_FOR);
			mm_align1(km, opt, mi, qlen, qseq0, &s[1], &s2[1], a, &ez, MM_F_SPLICE_REV);
			if (s[0].p->dp_score > s[1].p->dp_score) which = 0, trans_strand = 1;
			else if (s[0].p->dp_score < s[1].p->dp_score) which = 1, trans_strand = 2;
			else trans_strand = 3, which = (qlen + s[0].p->dp_score) & 1; // randomly choose a strand, effectively
			if (which == 0) {
				regs[i] = s[0], r2 = s2[0];
				free(s[1].p);
			} else {
				regs[i] = s[1], r2 = s2[1];
				free(s[0].p);
			}
			regs[i].p->trans_strand = trans_strand;
		} else {
			mm_align1(km, opt, mi, qlen, qseq0, &regs[i], &r2, a, &ez, opt->flag);
			if ((opt->flag&MM_F_SPLICE) && !(opt->flag&MM_F_SPLICE_BOTH))
				regs[i].p->trans_strand = opt->flag&MM_F_SPLICE_FOR? 1 : 2;
		}
		if (r2.cnt > 0) regs = mm_insert_reg(&r2, i, &n_regs, regs);
		if (i > 0 && mm_align1_inv(km, opt, mi, qlen, qseq0, &regs[i-1], &regs[i], &r2, &ez)) {
			regs = mm_insert_reg(&r2, i, &n_regs, regs);
			++i; // skip the inserted INV alignment
		}
	}
	*n_regs_ = n_regs;
	kfree(km, qseq0[0]); kfree(km, qseq0[1]);
	kfree(km, ez.cigar);
	mm_filter_regs(km, opt, n_regs_, regs);
	mm_hit_sort_by_dp(km, n_regs_, regs);
	return regs;
}
