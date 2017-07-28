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
		} else if (op == 2) {
			score -= q + e * len, i += len;
			if (test_zdrop_aux(score, i, j, &max, &max_i, &max_j, e, zdrop)) return 1;
		}
	}
	return 0;
}

static void mm_update_extra(mm_extra_t *p, const uint8_t *qseq, const uint8_t *tseq, const int8_t *mat, int8_t q, int8_t e)
{
	uint32_t k, l, toff = 0, qoff = 0;
	int32_t s = 0, max = 0;
	if (p == 0) return;
	for (k = 0; k < p->n_cigar; ++k) {
		uint32_t op = p->cigar[k]&0xf, len = p->cigar[k]>>4;
		if (op == 0) {
			for (l = 0; l < len; ++l) {
				int cq = qseq[qoff + l], ct = tseq[toff + l];
				if (ct > 3 || cq > 3) ++p->n_ambi;
				else if (ct != cq) ++p->n_diff;
				s += mat[ct * 5 + cq];
				if (s < 0) s = 0;
				else max = max > s? max : s;
			}
			toff += len, qoff += len, p->blen += len;
		} else if (op == 1) {
			int n_ambi = 0;
			for (l = 0; l < len; ++l)
				if (qseq[qoff + l] > 3) ++n_ambi;
			qoff += len, p->blen += len;
			p->n_ambi += n_ambi, p->n_diff += len - n_ambi;
			s -= q + e * len;
			if (s < 0) s = 0;
		} else if (op == 2) {
			int n_ambi = 0;
			for (l = 0; l < len; ++l)
				if (tseq[toff + l] > 3) ++n_ambi;
			toff += len, p->blen += len;
			p->n_ambi += n_ambi, p->n_diff += len - n_ambi;
			s -= q + e * len;
			if (s < 0) s = 0;
		}
	}
	p->dp_max = max;
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
	if (opt->q == opt->q2 && opt->e == opt->e2)
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
		*r = (int32_t)a->x + 1;
		*q = (int32_t)a->y + 1;
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

static void mm_align1(void *km, const mm_mapopt_t *opt, const mm_idx_t *mi, int qlen, uint8_t *qseq0[2], mm_reg1_t *r, mm_reg1_t *r2, mm128_t *a, ksw_extz_t *ez)
{
	int32_t rid = a[r->as].x<<1>>33, rev = a[r->as].x>>63, as1, cnt1;
	uint8_t *tseq, *qseq;
	int32_t i, l, bw, dropped = 0, rs0, re0, qs0, qe0;
	int32_t rs, re, qs, qe;
	int32_t rs1, qs1, re1, qe1;
	int8_t mat[25];

	if (r->cnt == 0) return;
	ksw_gen_simple_mat(5, mat, opt->a, opt->b);
	bw = (int)(opt->bw * 1.5 + 1.);

	r2->cnt = 0;
	mm_fix_bad_ends(r, a, opt->bw, &as1, &cnt1);
	mm_filter_bad_seeds(km, as1, cnt1, a, 10, 40, opt->max_gap>>1, 10);
	mm_adjust_minier(mi, qseq0, &a[as1], &rs, &qs);
	mm_adjust_minier(mi, qseq0, &a[as1 + cnt1 - 1], &re, &qe);

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
		mm_align_pair(km, opt, qs - qs0, qseq, rs - rs0, tseq, mat, bw, KSW_EZ_EXTZ_ONLY|KSW_EZ_RIGHT|KSW_EZ_REV_CIGAR, ez);
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
		if (a[as1+i].y & (MM_SEED_IGNORE|MM_SEED_TANDEM)) continue;
		mm_adjust_minier(mi, qseq0, &a[as1 + i], &re, &qe);
		re1 = re, qe1 = qe;
		if (i == cnt1 - 1 || (a[as1+i].y&MM_SEED_LONG_JOIN) || (qe - qs >= opt->min_ksw_len && re - rs >= opt->min_ksw_len)) {
			int bw1 = bw;
			if (a[as1+i].y & MM_SEED_LONG_JOIN)
				bw1 = qe - qs > re - rs? qe - qs : re - rs;
			qseq = &qseq0[rev][qs];
			mm_idx_getseq(mi, rid, rs, re, tseq);
			mm_align_pair(km, opt, qe - qs, qseq, re - rs, tseq, mat, bw1, KSW_EZ_APPROX_MAX, ez);
			if (mm_check_zdrop(qseq, tseq, ez->n_cigar, ez->cigar, mat, opt->q, opt->e, opt->zdrop))
				mm_align_pair(km, opt, qe - qs, qseq, re - rs, tseq, mat, bw1, 0, ez);
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
		mm_align_pair(km, opt, qe0 - qe, qseq, re0 - re, tseq, mat, bw, KSW_EZ_EXTZ_ONLY, ez);
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
	}

	kfree(km, tseq);
}

mm_reg1_t *mm_align_skeleton(void *km, const mm_mapopt_t *opt, const mm_idx_t *mi, int qlen, const char *qstr, int *n_regs_, mm_reg1_t *regs, mm128_t *a)
{
	extern unsigned char seq_nt4_table[256];
	int32_t i, r, n_regs = *n_regs_;
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
	kfree(km, qseq0[0]); kfree(km, qseq0[1]);
	kfree(km, ez.cigar);
	mm_filter_regs(km, opt, n_regs_, regs);
	mm_hit_sort_by_dp(km, n_regs_, regs);
	return regs;
}
