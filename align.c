#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "minimap.h"
#include "mmpriv.h"
#include "ksw2.h"

static void ksw_gen_simple_mat(int m, int8_t *mat, int8_t a, int8_t b, int8_t sc_ambi)
{
	int i, j;
	a = a < 0? -a : a;
	b = b > 0? -b : b;
	sc_ambi = sc_ambi > 0? -sc_ambi : sc_ambi;
	for (i = 0; i < m - 1; ++i) {
		for (j = 0; j < m - 1; ++j)
			mat[i * m + j] = i == j? a : b;
		mat[i * m + m - 1] = sc_ambi;
	}
	for (j = 0; j < m; ++j)
		mat[(m - 1) * m + j] = sc_ambi;
}

static void ksw_gen_ts_mat(int m, int8_t *mat, int8_t a, int8_t b, int8_t transition, int8_t sc_ambi)
{
	assert(m == 5);
	ksw_gen_simple_mat(m, mat, a, b, sc_ambi);
	if (transition == 0 || transition == b) return;
	transition = transition > 0? -transition : transition;
	mat[0 * m + 2] = transition;  // A->G
	mat[1 * m + 3] = transition;  // C->T
	mat[2 * m + 0] = transition;  // G->A
	mat[3 * m + 1] = transition;  // T->C
}

static inline void mm_seq_rev(uint32_t len, uint8_t *seq)
{
	uint32_t i;
	uint8_t t;
	for (i = 0; i < len>>1; ++i)
		t = seq[i], seq[i] = seq[len - 1 - i], seq[len - 1 - i] = t;
}

static inline void update_max_zdrop(int32_t score, int i, int j, int32_t *max, int *max_i, int *max_j, int e, int *max_zdrop, int pos[2][2])
{
	if (score < *max) {
		int li = i - *max_i;
		int lj = j - *max_j;
		int diff = li > lj? li - lj : lj - li;
		int z = *max - score - diff * e;
		if (z > *max_zdrop) {
			*max_zdrop = z;
			pos[0][0] = *max_i, pos[0][1] = i;
			pos[1][0] = *max_j, pos[1][1] = j;
		}
	} else *max = score, *max_i = i, *max_j = j;
}

static int mm_test_zdrop(void *km, const mm_mapopt_t *opt, const uint8_t *qseq, const uint8_t *tseq, uint32_t n_cigar, uint32_t *cigar, const int8_t *mat)
{
	uint32_t k;
	int32_t score = 0, max = INT32_MIN, max_i = -1, max_j = -1, i = 0, j = 0, max_zdrop = 0;
	int pos[2][2] = {{-1, -1}, {-1, -1}}, q_len, t_len;

	// find the score and the region where score drops most along diagonal
	for (k = 0, score = 0; k < n_cigar; ++k) {
		uint32_t l, op = cigar[k]&0xf, len = cigar[k]>>4;
		if (op == MM_CIGAR_MATCH) {
			for (l = 0; l < len; ++l) {
				score += mat[tseq[i + l] * 5 + qseq[j + l]];
				update_max_zdrop(score, i+l, j+l, &max, &max_i, &max_j, opt->e, &max_zdrop, pos);
			}
			i += len, j += len;
		} else if (op == MM_CIGAR_INS || op == MM_CIGAR_DEL || op == MM_CIGAR_N_SKIP) {
			score -= opt->q + opt->e * len;
			if (op == MM_CIGAR_INS) j += len;
			else i += len;
			update_max_zdrop(score, i, j, &max, &max_i, &max_j, opt->e, &max_zdrop, pos);
		}
	}

	// test if there is an inversion in the most dropped region
	q_len = pos[1][1] - pos[1][0], t_len = pos[0][1] - pos[0][0];
	if (!(opt->flag&(MM_F_SPLICE|MM_F_SR|MM_F_FOR_ONLY|MM_F_REV_ONLY)) && max_zdrop > opt->zdrop_inv && q_len < opt->max_gap && t_len < opt->max_gap) {
		uint8_t *qseq2;
		void *qp;
		int q_off, t_off;
		qseq2 = (uint8_t*)kmalloc(km, q_len);
		for (i = 0; i < q_len; ++i) {
			int c = qseq[pos[1][1] - i - 1];
			qseq2[i] = c >= 4? 4 : 3 - c;
		}
		qp = ksw_ll_qinit(km, 2, q_len, qseq2, 5, mat);
		score = ksw_ll_i16(qp, t_len, tseq + pos[0][0], opt->q, opt->e, &q_off, &t_off);
		kfree(km, qseq2);
		kfree(km, qp);
		if (score >= opt->min_chain_score * opt->a && score >= opt->min_dp_max)
			return 2; // there is a potential inversion
	}
	return max_zdrop > opt->zdrop? 1 : 0;
}

static void mm_fix_cigar(mm_reg1_t *r, const uint8_t *qseq, const uint8_t *tseq, int *qshift, int *tshift)
{
	mm_extra_t *p = r->p;
	int32_t toff = 0, qoff = 0, to_shrink = 0;
	uint32_t k;
	*qshift = *tshift = 0;
	if (p->n_cigar <= 1) return;
	for (k = 0; k < p->n_cigar; ++k) { // indel left alignment
		uint32_t op = p->cigar[k]&0xf, len = p->cigar[k]>>4;
		if (len == 0) to_shrink = 1;
		if (op == MM_CIGAR_MATCH) {
			toff += len, qoff += len;
		} else if (op == MM_CIGAR_INS || op == MM_CIGAR_DEL) {
			if (k > 0 && k < p->n_cigar - 1 && (p->cigar[k-1]&0xf) == 0 && (p->cigar[k+1]&0xf) == 0) {
				int l, prev_len = p->cigar[k-1] >> 4;
				if (op == MM_CIGAR_INS) {
					for (l = 0; l < prev_len; ++l)
						if (qseq[qoff - 1 - l] != qseq[qoff + len - 1 - l])
							break;
				} else {
					for (l = 0; l < prev_len; ++l)
						if (tseq[toff - 1 - l] != tseq[toff + len - 1 - l])
							break;
				}
				if (l > 0)
					p->cigar[k-1] -= l<<4, p->cigar[k+1] += l<<4, qoff -= l, toff -= l;
				if (l == prev_len) to_shrink = 1;
			}
			if (op == MM_CIGAR_INS) qoff += len;
			else toff += len;
		} else if (op == MM_CIGAR_N_SKIP) {
			toff += len;
		}
	}
	assert(qoff == r->qe - r->qs && toff == r->re - r->rs);
	for (k = 0; k < p->n_cigar - 2; ++k) { // fix CIGAR like 5I6D7I
		if ((p->cigar[k]&0xf) > 0 && (p->cigar[k]&0xf) + (p->cigar[k+1]&0xf) == 3) {
			uint32_t l, s[3] = {0,0,0};
			for (l = k; l < p->n_cigar; ++l) { // count number of adjacent I and D
				uint32_t op = p->cigar[l]&0xf;
				if (op == MM_CIGAR_INS || op == MM_CIGAR_DEL || p->cigar[l]>>4 == 0)
					s[op] += p->cigar[l] >> 4;
				else break;
			}
			if (s[1] > 0 && s[2] > 0 && l - k > 2) { // turn to a single I and a single D
				p->cigar[k]   = s[1]<<4|MM_CIGAR_INS;
				p->cigar[k+1] = s[2]<<4|MM_CIGAR_DEL;
				for (k += 2; k < l; ++k)
					p->cigar[k] &= 0xf;
				to_shrink = 1;
			}
			k = l;
		}
	}
	if (to_shrink) { // squeeze out zero-length operations
		int32_t l = 0;
		for (k = 0; k < p->n_cigar; ++k) // squeeze out zero-length operations
			if (p->cigar[k]>>4 != 0)
				p->cigar[l++] = p->cigar[k];
		p->n_cigar = l;
		for (k = l = 0; k < p->n_cigar; ++k) // merge two adjacent operations if they are the same
			if (k == p->n_cigar - 1 || (p->cigar[k]&0xf) != (p->cigar[k+1]&0xf))
				p->cigar[l++] = p->cigar[k];
			else p->cigar[k+1] += p->cigar[k]>>4<<4; // add length to the next CIGAR operator
		p->n_cigar = l;
	}
	if ((p->cigar[0]&0xf) == MM_CIGAR_INS || (p->cigar[0]&0xf) == MM_CIGAR_DEL) { // get rid of leading I or D
		int32_t l = p->cigar[0] >> 4;
		if ((p->cigar[0]&0xf) == MM_CIGAR_INS) {
			if (r->rev) r->qe -= l;
			else r->qs += l;
			*qshift = l;
		} else r->rs += l, *tshift = l;
		--p->n_cigar;
		memmove(p->cigar, p->cigar + 1, p->n_cigar * 4);
	}
}

static void mm_update_cigar_eqx(mm_reg1_t *r, const uint8_t *qseq, const uint8_t *tseq) // written by @armintoepfer
{
	uint32_t n_EQX = 0;
	uint32_t k, l, m, cap, toff = 0, qoff = 0, n_M = 0;
	mm_extra_t *p;
	if (r->p == 0) return;
	for (k = 0; k < r->p->n_cigar; ++k) {
		uint32_t op = r->p->cigar[k]&0xf, len = r->p->cigar[k]>>4;
		if (op == MM_CIGAR_MATCH) {
			while (len > 0) {
				for (l = 0; l < len && qseq[qoff + l] == tseq[toff + l]; ++l) {} // run of "="; TODO: N<=>N is converted to "="
				if (l > 0) { ++n_EQX; len -= l; toff += l; qoff += l; }

				for (l = 0; l < len && qseq[qoff + l] != tseq[toff + l]; ++l) {} // run of "X"
				if (l > 0) { ++n_EQX; len -= l; toff += l; qoff += l; }
			}
			++n_M;
		} else if (op == MM_CIGAR_INS) {
			qoff += len;
		} else if (op == MM_CIGAR_DEL) {
			toff += len;
		} else if (op == MM_CIGAR_N_SKIP) {
			toff += len;
		}
	}
	// update in-place if we can
	if (n_EQX == n_M) {
		for (k = 0; k < r->p->n_cigar; ++k) {
			uint32_t op = r->p->cigar[k]&0xf, len = r->p->cigar[k]>>4;
			if (op == MM_CIGAR_MATCH) r->p->cigar[k] = len << 4 | MM_CIGAR_EQ_MATCH;
		}
		return;
	}
	// allocate new storage
	cap = r->p->n_cigar + (n_EQX - n_M) + sizeof(mm_extra_t);
	kroundup32(cap);
	p = (mm_extra_t*)calloc(cap, 4);
	memcpy(p, r->p, sizeof(mm_extra_t));
	p->capacity = cap;
	// update cigar while copying
	toff = qoff = m = 0;
	for (k = 0; k < r->p->n_cigar; ++k) {
		uint32_t op = r->p->cigar[k]&0xf, len = r->p->cigar[k]>>4;
		if (op == MM_CIGAR_MATCH) {
			while (len > 0) {
				// match
				for (l = 0; l < len && qseq[qoff + l] == tseq[toff + l]; ++l) {}
				if (l > 0) p->cigar[m++] = l << 4 | MM_CIGAR_EQ_MATCH;
				len -= l;
				toff += l, qoff += l;
				// mismatch
				for (l = 0; l < len && qseq[qoff + l] != tseq[toff + l]; ++l) {}
				if (l > 0) p->cigar[m++] = l << 4 | MM_CIGAR_X_MISMATCH;
				len -= l;
				toff += l, qoff += l;
			}
			continue;
		} else if (op == MM_CIGAR_INS) {
			qoff += len;
		} else if (op == MM_CIGAR_DEL) {
			toff += len;
		} else if (op == MM_CIGAR_N_SKIP) {
			toff += len;
		}
		p->cigar[m++] = r->p->cigar[k];
	}
	p->n_cigar = m;
	free(r->p);
	r->p = p;
}

static void mm_update_extra(mm_reg1_t *r, const uint8_t *qseq, const uint8_t *tseq, const int8_t *mat, int8_t q, int8_t e, int is_eqx, int log_gap)
{
	uint32_t k, l;
	int32_t qshift, tshift, toff = 0, qoff = 0;
	double s = 0.0, max = 0.0;
	mm_extra_t *p = r->p;
	if (p == 0) return;
	mm_fix_cigar(r, qseq, tseq, &qshift, &tshift);
	qseq += qshift, tseq += tshift; // qseq and tseq may be shifted due to the removal of leading I/D
	r->blen = r->mlen = 0;
	for (k = 0; k < p->n_cigar; ++k) {
		uint32_t op = p->cigar[k]&0xf, len = p->cigar[k]>>4;
		if (op == MM_CIGAR_MATCH) {
			int n_ambi = 0, n_diff = 0;
			for (l = 0; l < len; ++l) {
				int cq = qseq[qoff + l], ct = tseq[toff + l];
				if (ct > 3 || cq > 3) ++n_ambi;
				else if (ct != cq) ++n_diff;
				s += mat[ct * 5 + cq];
				if (s < 0) s = 0;
				else max = max > s? max : s;
			}
			r->blen += len - n_ambi, r->mlen += len - (n_ambi + n_diff), p->n_ambi += n_ambi;
			toff += len, qoff += len;
		} else if (op == MM_CIGAR_INS) {
			int n_ambi = 0;
			for (l = 0; l < len; ++l)
				if (qseq[qoff + l] > 3) ++n_ambi;
			r->blen += len - n_ambi, p->n_ambi += n_ambi;
			if (log_gap) s -= q + (double)e * mg_log2(1.0 + len);
			else s -= q + e;
			if (s < 0) s = 0;
			qoff += len;
		} else if (op == MM_CIGAR_DEL) {
			int n_ambi = 0;
			for (l = 0; l < len; ++l)
				if (tseq[toff + l] > 3) ++n_ambi;
			r->blen += len - n_ambi, p->n_ambi += n_ambi;
			if (log_gap) s -= q + (double)e * mg_log2(1.0 + len);
			else s -= q + e;
			if (s < 0) s = 0;
			toff += len;
		} else if (op == MM_CIGAR_N_SKIP) {
			toff += len;
		}
	}
	p->dp_max = p->dp_max0 = (int32_t)(max + .499);
	assert(qoff == r->qe - r->qs && toff == r->re - r->rs);
	if (is_eqx) mm_update_cigar_eqx(r, qseq, tseq); // NB: it has to be called here as changes to qseq and tseq are not returned
}

static void mm_append_cigar(mm_reg1_t *r, uint32_t n_cigar, uint32_t *cigar) // TODO: this calls the libc realloc()
{
	mm_extra_t *p;
	if (n_cigar == 0) return;
	if (r->p == 0) {
		uint32_t capacity = n_cigar + sizeof(mm_extra_t)/4;
		kroundup32(capacity);
		r->p = (mm_extra_t*)calloc(capacity, 4);
		r->p->capacity = capacity;
	} else if (r->p->n_cigar + n_cigar + sizeof(mm_extra_t)/4 > r->p->capacity) {
		r->p->capacity = r->p->n_cigar + n_cigar + sizeof(mm_extra_t)/4;
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

static void mm_align_pair(void *km, const mm_mapopt_t *opt, int qlen, const uint8_t *qseq, int tlen, const uint8_t *tseq, const uint8_t *junc, const int8_t *mat, int w, int end_bonus, int zdrop, int flag, ksw_extz_t *ez)
{
	if (mm_dbg_flag & MM_DBG_PRINT_ALN_SEQ) {
		int i;
		fprintf(stderr, "===> q=(%d,%d), e=(%d,%d), bw=%d, flag=%d, zdrop=%d <===\n", opt->q, opt->q2, opt->e, opt->e2, w, flag, opt->zdrop);
		for (i = 0; i < tlen; ++i) fputc("ACGTN"[tseq[i]], stderr);
		fputc('\n', stderr);
		for (i = 0; i < qlen; ++i) fputc("ACGTN"[qseq[i]], stderr);
		fputc('\n', stderr);
	}
	if (opt->transition != 0 && opt->b != opt->transition)
		flag |= KSW_EZ_GENERIC_SC;
	if (opt->max_sw_mat > 0 && (int64_t)tlen * qlen > opt->max_sw_mat) {
		ksw_reset_extz(ez);
		ez->zdropped = 1;
	} else if (opt->flag & MM_F_SPLICE) {
		int flag_tmp = flag;
		if (!(opt->flag & MM_F_SPLICE_OLD)) flag_tmp |= KSW_EZ_SPLICE_CMPLX;
		ksw_exts2_sse(km, qlen, qseq, tlen, tseq, 5, mat, opt->q, opt->e, opt->q2, opt->noncan, zdrop, opt->junc_bonus, flag_tmp, junc, ez);
	} else if (opt->q == opt->q2 && opt->e == opt->e2)
		ksw_extz2_sse(km, qlen, qseq, tlen, tseq, 5, mat, opt->q, opt->e, w, zdrop, end_bonus, flag, ez);
	else
		ksw_extd2_sse(km, qlen, qseq, tlen, tseq, 5, mat, opt->q, opt->e, opt->q2, opt->e2, w, zdrop, end_bonus, flag, ez);
	if (mm_dbg_flag & MM_DBG_PRINT_ALN_SEQ) {
		int i;
		fprintf(stderr, "score=%d, cigar=", ez->score);
		for (i = 0; i < ez->n_cigar; ++i)
			fprintf(stderr, "%d%c", ez->cigar[i]>>4, MM_CIGAR_STR[ez->cigar[i]&0xf]);
		fprintf(stderr, "\n");
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
	if (mi->flag & MM_I_HPC) {
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

static int *collect_long_gaps(void *km, int as1, int cnt1, mm128_t *a, int min_gap, int *n_)
{
	int i, n, *K;
	*n_ = 0;
	for (i = 1, n = 0; i < cnt1; ++i) { // count the number of gaps longer than min_gap
		int gap = ((int32_t)a[as1 + i].y - a[as1 + i - 1].y) - ((int32_t)a[as1 + i].x - a[as1 + i - 1].x);
		if (gap < -min_gap || gap > min_gap) ++n;
	}
	if (n <= 1) return 0;
	K = (int*)kmalloc(km, n * sizeof(int));
	for (i = 1, n = 0; i < cnt1; ++i) { // store the positions of long gaps
		int gap = ((int32_t)a[as1 + i].y - a[as1 + i - 1].y) - ((int32_t)a[as1 + i].x - a[as1 + i - 1].x);
		if (gap < -min_gap || gap > min_gap)
			K[n++] = i;
	}
	*n_ = n;
	return K;
}

static void mm_filter_bad_seeds(void *km, int as1, int cnt1, mm128_t *a, int min_gap, int diff_thres, int max_ext_len, int max_ext_cnt)
{
	int max_st, max_en, n, i, k, max, *K;
	K = collect_long_gaps(km, as1, cnt1, a, min_gap, &n);
	if (K == 0) return;
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
		gap = ((int32_t)a[as1 + i].y - (int32_t)a[as1 + i - 1].y) - (int32_t)(a[as1 + i].x - a[as1 + i - 1].x);
		if (gap > 0) n_ins += gap;
		else n_del += -gap;
		qs = (int32_t)a[as1 + i - 1].y;
		rs = (int32_t)a[as1 + i - 1].x;
		for (l = k + 1; l < n && l <= k + max_ext_cnt; ++l) {
			int j = K[l], diff;
			if ((int32_t)a[as1 + j].y - qs > max_ext_len || (int32_t)a[as1 + j].x - rs > max_ext_len) break;
			gap = ((int32_t)a[as1 + j].y - (int32_t)a[as1 + j - 1].y) - (int32_t)(a[as1 + j].x - a[as1 + j - 1].x);
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

static void mm_filter_bad_seeds_alt(void *km, int as1, int cnt1, mm128_t *a, int min_gap, int max_ext)
{
	int n, k, *K;
	K = collect_long_gaps(km, as1, cnt1, a, min_gap, &n);
	if (K == 0) return;
	for (k = 0; k < n;) {
		int i = K[k], l;
		int gap1 = ((int32_t)a[as1 + i].y - (int32_t)a[as1 + i - 1].y) - ((int32_t)a[as1 + i].x - (int32_t)a[as1 + i - 1].x);
		int re1 = (int32_t)a[as1 + i].x;
		int qe1 = (int32_t)a[as1 + i].y;
		gap1 = gap1 > 0? gap1 : -gap1;
		for (l = k + 1; l < n; ++l) {
			int j = K[l], gap2, q_span_pre, rs2, qs2, m;
			if ((int32_t)a[as1 + j].y - qe1 > max_ext || (int32_t)a[as1 + j].x - re1 > max_ext) break;
			gap2 = ((int32_t)a[as1 + j].y - (int32_t)a[as1 + j - 1].y) - (int32_t)(a[as1 + j].x - a[as1 + j - 1].x);
			q_span_pre = a[as1 + j - 1].y >> 32 & 0xff;
			rs2 = (int32_t)a[as1 + j - 1].x + q_span_pre;
			qs2 = (int32_t)a[as1 + j - 1].y + q_span_pre;
			m = rs2 - re1 < qs2 - qe1? rs2 - re1 : qs2 - qe1;
			gap2 = gap2 > 0? gap2 : -gap2;
			if (m > gap1 + gap2) break;
			re1 = (int32_t)a[as1 + j].x;
			qe1 = (int32_t)a[as1 + j].y;
			gap1 = gap2;
		}
		if (l > k + 1) {
			int j, end = K[l - 1];
			for (j = K[k]; j < end; ++j)
				a[as1 + j].y |= MM_SEED_IGNORE;
			a[as1 + end].y |= MM_SEED_LONG_JOIN;
		}
		k = l;
	}
	kfree(km, K);
}

static void mm_fix_bad_ends(const mm_reg1_t *r, const mm128_t *a, int bw, int min_match, int32_t *as, int32_t *cnt)
{
	int32_t i, l, m;
	*as = r->as, *cnt = r->cnt;
	if (r->cnt < 3) return;
	m = l = a[r->as].y >> 32 & 0xff;
	for (i = r->as + 1; i < r->as + r->cnt - 1; ++i) {
		int32_t lq, lr, min, max;
		int32_t q_span = a[i].y >> 32 & 0xff;
		if (a[i].y & MM_SEED_LONG_JOIN) break;
		lr = (int32_t)a[i].x - (int32_t)a[i-1].x;
		lq = (int32_t)a[i].y - (int32_t)a[i-1].y;
		min = lr < lq? lr : lq;
		max = lr > lq? lr : lq;
		if (max - min > l >> 1) *as = i;
		l += min;
		m += min < q_span? min : q_span;
		if (l >= bw << 1 || (m >= min_match && m >= bw) || m >= r->mlen >> 1) break;
	}
	*cnt = r->as + r->cnt - *as;
	m = l = a[r->as + r->cnt - 1].y >> 32 & 0xff;
	for (i = r->as + r->cnt - 2; i > *as; --i) {
		int32_t lq, lr, min, max;
		int32_t q_span = a[i+1].y >> 32 & 0xff;
		if (a[i+1].y & MM_SEED_LONG_JOIN) break;
		lr = (int32_t)a[i+1].x - (int32_t)a[i].x;
		lq = (int32_t)a[i+1].y - (int32_t)a[i].y;
		min = lr < lq? lr : lq;
		max = lr > lq? lr : lq;
		if (max - min > l >> 1) *cnt = i + 1 - *as;
		l += min;
		m += min < q_span? min : q_span;
		if (l >= bw << 1 || (m >= min_match && m >= bw) || m >= r->mlen >> 1) break;
	}
}

static void mm_max_stretch(const mm_reg1_t *r, const mm128_t *a, int32_t *as, int32_t *cnt)
{
	int32_t i, score, max_score, len, max_i, max_len;

	*as = r->as, *cnt = r->cnt;
	if (r->cnt < 2) return;

	max_score = -1, max_i = -1, max_len = 0;
	score = a[r->as].y >> 32 & 0xff, len = 1;
	for (i = r->as + 1; i < r->as + r->cnt; ++i) {
		int32_t lq, lr, q_span;
		q_span = a[i].y >> 32 & 0xff;
		lr = (int32_t)a[i].x - (int32_t)a[i-1].x;
		lq = (int32_t)a[i].y - (int32_t)a[i-1].y;
		if (lq == lr) {
			score += lq < q_span? lq : q_span;
			++len;
		} else {
			if (score > max_score)
				max_score = score, max_len = len, max_i = i - len;
			score = q_span, len = 1;
		}
	}
	if (score > max_score)
		max_score = score, max_len = len, max_i = i - len;
	*as = max_i, *cnt = max_len;
}

static int mm_seed_ext_score(void *km, const mm_mapopt_t *opt, const mm_idx_t *mi, const int8_t mat[25], int qlen, uint8_t *qseq0[2], const mm128_t *a)
{
	uint8_t *qseq, *tseq;
	int q_span = a->y>>32&0xff, qs, qe, rs, re, rid, score, q_off, t_off, ext_len = opt->anchor_ext_len;
	void *qp;
	rid = a->x<<1>>33;
	re = (uint32_t)a->x + 1, rs = re - q_span;
	qe = (uint32_t)a->y + 1, qs = qe - q_span;
	rs = rs - ext_len > 0? rs - ext_len : 0;
	qs = qs - ext_len > 0? qs - ext_len : 0;
	re = re + ext_len < (int32_t)mi->seq[rid].len? re + ext_len : mi->seq[rid].len;
	qe = qe + ext_len < qlen? qe + ext_len : qlen;
	tseq = (uint8_t*)kmalloc(km, re - rs);
	if (opt->flag & MM_F_QSTRAND) {
		qseq = qseq0[0] + qs;
		mm_idx_getseq2(mi, a->x>>63, rid, rs, re, tseq);
	} else {
		qseq = qseq0[a->x>>63] + qs;
		mm_idx_getseq(mi, rid, rs, re, tseq);
	}
	qp = ksw_ll_qinit(km, 2, qe - qs, qseq, 5, mat);
	score = ksw_ll_i16(qp, re - rs, tseq, opt->q, opt->e, &q_off, &t_off);
	kfree(km, tseq);
	kfree(km, qp);
	return score;
}

static void mm_fix_bad_ends_splice(void *km, const mm_mapopt_t *opt, const mm_idx_t *mi, const mm_reg1_t *r, const int8_t mat[25], int qlen, uint8_t *qseq0[2], const mm128_t *a, int *as1, int *cnt1)
{ // this assumes a very crude k-mer based mode; it is not necessary to use a good model just for filtering bounary exons
	int score;
	double log_gap;
	*as1 = r->as, *cnt1 = r->cnt;
	if (r->cnt < 3) return;
	log_gap = log((int32_t)a[r->as + 1].x - (int32_t)a[r->as].x);
	if ((a[r->as].y>>32&0xff) < log_gap + opt->anchor_ext_shift) {
		score = mm_seed_ext_score(km, opt, mi, mat, qlen, qseq0, &a[r->as]);
		if ((double)score / mat[0] < log_gap + opt->anchor_ext_shift) // a more exact format is "score < log_4(gap) + shift"
			++(*as1), --(*cnt1);
	}
	log_gap = log((int32_t)a[r->as + r->cnt - 1].x - (int32_t)a[r->as + r->cnt - 2].x);
	if ((a[r->as + r->cnt - 1].y>>32&0xff) < log_gap + opt->anchor_ext_shift) {
		score = mm_seed_ext_score(km, opt, mi, mat, qlen, qseq0, &a[r->as + r->cnt - 1]);
		if ((double)score / mat[0] < log_gap + opt->anchor_ext_shift)
			--(*cnt1);
	}
}

static void mm_align1(void *km, const mm_mapopt_t *opt, const mm_idx_t *mi, int qlen, uint8_t *qseq0[2], mm_reg1_t *r, mm_reg1_t *r2, int n_a, mm128_t *a, ksw_extz_t *ez, int splice_flag)
{
	int is_sr = !!(opt->flag & MM_F_SR), is_splice = !!(opt->flag & MM_F_SPLICE);
	int32_t rid = a[r->as].x<<1>>33, rev = a[r->as].x>>63, as1, cnt1;
	uint8_t *tseq, *qseq, *junc;
	int32_t i, l, bw, bw_long, dropped = 0, extra_flag = 0, rs0, re0, qs0, qe0;
	int32_t rs, re, qs, qe;
	int32_t rs1, qs1, re1, qe1;
	int8_t mat[25];

	if (is_sr) assert(!(mi->flag & MM_I_HPC)); // HPC won't work with SR because with HPC we can't easily tell if there is a gap

	r2->cnt = 0;
	if (r->cnt == 0) return;
	ksw_gen_ts_mat(5, mat, opt->a, opt->b, opt->transition, opt->sc_ambi);
	bw = (int)(opt->bw * 1.5 + 1.);
	bw_long = (int)(opt->bw_long * 1.5 + 1.);
	if (bw_long < bw) bw_long = bw;

	if (is_sr && !(mi->flag & MM_I_HPC)) {
		mm_max_stretch(r, a, &as1, &cnt1);
		rs = (int32_t)a[as1].x + 1 - (int32_t)(a[as1].y>>32&0xff);
		qs = (int32_t)a[as1].y + 1 - (int32_t)(a[as1].y>>32&0xff);
		re = (int32_t)a[as1+cnt1-1].x + 1;
		qe = (int32_t)a[as1+cnt1-1].y + 1;
	} else {
		if (!(opt->flag & MM_F_NO_END_FLT)) {
			if (is_splice)
				mm_fix_bad_ends_splice(km, opt, mi, r, mat, qlen, qseq0, a, &as1, &cnt1);
			else
				mm_fix_bad_ends(r, a, opt->bw, opt->min_chain_score * 2, &as1, &cnt1);
		} else as1 = r->as, cnt1 = r->cnt;
		mm_filter_bad_seeds(km, as1, cnt1, a, 10, 40, opt->max_gap>>1, 10);
		mm_filter_bad_seeds_alt(km, as1, cnt1, a, 30, opt->max_gap>>1);
		mm_adjust_minier(mi, qseq0, &a[as1], &rs, &qs);
		mm_adjust_minier(mi, qseq0, &a[as1 + cnt1 - 1], &re, &qe);
	}
	assert(cnt1 > 0);

	if (is_splice) {
		if (splice_flag & MM_F_SPLICE_FOR) extra_flag |= rev? KSW_EZ_SPLICE_REV : KSW_EZ_SPLICE_FOR;
		if (splice_flag & MM_F_SPLICE_REV) extra_flag |= rev? KSW_EZ_SPLICE_FOR : KSW_EZ_SPLICE_REV;
		if (opt->flag & MM_F_SPLICE_FLANK) extra_flag |= KSW_EZ_SPLICE_FLANK;
	}

	/* Look for the start and end of regions to perform DP. This sounds easy
	 * but is in fact tricky. Excessively small regions lead to unnecessary
	 * clippings and lose alignable sequences. Excessively large regions
	 * occasionally lead to large overlaps between two chains and may cause
	 * loss of alignments in corner cases. */
	if (is_sr) {
		qs0 = 0, qe0 = qlen;
		l = qs;
		l += l * opt->a + opt->end_bonus > opt->q? (l * opt->a + opt->end_bonus - opt->q) / opt->e : 0;
		rs0 = rs - l > 0? rs - l : 0;
		l = qlen - qe;
		l += l * opt->a + opt->end_bonus > opt->q? (l * opt->a + opt->end_bonus - opt->q) / opt->e : 0;
		re0 = re + l < (int32_t)mi->seq[rid].len? re + l : mi->seq[rid].len;
	} else {
		// compute rs0 and qs0
		rs0 = (int32_t)a[r->as].x + 1 - (int32_t)(a[r->as].y>>32&0xff);
		qs0 = (int32_t)a[r->as].y + 1 - (int32_t)(a[r->as].y>>32&0xff);
		if (rs0 < 0) rs0 = 0; // this may happen when HPC is in use
		assert(qs0 >= 0); // this should never happen, or it is logic error
		rs1 = qs1 = 0;
		for (i = r->as - 1, l = 0; i >= 0 && a[i].x>>32 == a[r->as].x>>32; --i) { // inspect nearby seeds
			int32_t x = (int32_t)a[i].x + 1 - (int32_t)(a[i].y>>32&0xff);
			int32_t y = (int32_t)a[i].y + 1 - (int32_t)(a[i].y>>32&0xff);
			if (x < rs0 && y < qs0) {
				if (++l > opt->min_cnt) {
					l = rs0 - x > qs0 - y? rs0 - x : qs0 - y;
					rs1 = rs0 - l, qs1 = qs0 - l;
					if (rs1 < 0) rs1 = 0; // not strictly necessary; better have this guard for explicit
					break;
				}
			}
		}
		if (qs > 0 && rs > 0) {
			l = qs < opt->max_gap? qs : opt->max_gap;
			qs1 = qs1 > qs - l? qs1 : qs - l;
			qs0 = qs0 < qs1? qs0 : qs1; // at least include qs0
			l += l * opt->a > opt->q? (l * opt->a - opt->q) / opt->e : 0;
			l = l < opt->max_gap? l : opt->max_gap;
			l = l < rs? l : rs;
			rs1 = rs1 > rs - l? rs1 : rs - l;
			rs0 = rs0 < rs1? rs0 : rs1;
			rs0 = rs0 < rs? rs0 : rs;
		} else rs0 = rs, qs0 = qs;
		// compute re0 and qe0
		re0 = (int32_t)a[r->as + r->cnt - 1].x + 1;
		qe0 = (int32_t)a[r->as + r->cnt - 1].y + 1;
		re1 = mi->seq[rid].len, qe1 = qlen;
		for (i = r->as + r->cnt, l = 0; i < n_a && a[i].x>>32 == a[r->as].x>>32; ++i) { // inspect nearby seeds
			int32_t x = (int32_t)a[i].x + 1;
			int32_t y = (int32_t)a[i].y + 1;
			if (x > re0 && y > qe0) {
				if (++l > opt->min_cnt) {
					l = x - re0 > y - qe0? x - re0 : y - qe0;
					re1 = re0 + l, qe1 = qe0 + l;
					break;
				}
			}
		}
		if (qe < qlen && re < (int32_t)mi->seq[rid].len) {
			l = qlen - qe < opt->max_gap? qlen - qe : opt->max_gap;
			qe1 = qe1 < qe + l? qe1 : qe + l;
			qe0 = qe0 > qe1? qe0 : qe1; // at least include qe0
			l += l * opt->a > opt->q? (l * opt->a - opt->q) / opt->e : 0;
			l = l < opt->max_gap? l : opt->max_gap;
			l = l < (int32_t)mi->seq[rid].len - re? l : mi->seq[rid].len - re;
			re1 = re1 < re + l? re1 : re + l;
			re0 = re0 > re1? re0 : re1;
		} else re0 = re, qe0 = qe;
	}
	if (a[r->as].y & MM_SEED_SELF) {
		int max_ext = r->qs > r->rs? r->qs - r->rs : r->rs - r->qs;
		if (r->rs - rs0 > max_ext) rs0 = r->rs - max_ext;
		if (r->qs - qs0 > max_ext) qs0 = r->qs - max_ext;
		max_ext = r->qe > r->re? r->qe - r->re : r->re - r->qe;
		if (re0 - r->re > max_ext) re0 = r->re + max_ext;
		if (qe0 - r->qe > max_ext) qe0 = r->qe + max_ext;
	}

	assert(re0 > rs0);
	tseq = (uint8_t*)kmalloc(km, re0 - rs0);
	junc = (uint8_t*)kmalloc(km, re0 - rs0);

	if (qs > 0 && rs > 0) { // left extension; probably the condition can be changed to "qs > qs0 && rs > rs0"
		if (opt->flag & MM_F_QSTRAND) {
			qseq = &qseq0[0][qs0];
			mm_idx_getseq2(mi, rev, rid, rs0, rs, tseq);
		} else {
			qseq = &qseq0[rev][qs0];
			mm_idx_getseq(mi, rid, rs0, rs, tseq);
		}
		mm_idx_bed_junc(mi, rid, rs0, rs, junc);
		mm_seq_rev(qs - qs0, qseq);
		mm_seq_rev(rs - rs0, tseq);
		mm_seq_rev(rs - rs0, junc);
		mm_align_pair(km, opt, qs - qs0, qseq, rs - rs0, tseq, junc, mat, bw, opt->end_bonus, r->split_inv? opt->zdrop_inv : opt->zdrop, extra_flag|KSW_EZ_EXTZ_ONLY|KSW_EZ_RIGHT|KSW_EZ_REV_CIGAR, ez);
		if (ez->n_cigar > 0) {
			mm_append_cigar(r, ez->n_cigar, ez->cigar);
			r->p->dp_score += ez->max;
		}
		rs1 = rs - (ez->reach_end? ez->mqe_t + 1 : ez->max_t + 1);
		qs1 = qs - (ez->reach_end? qs - qs0 : ez->max_q + 1);
		mm_seq_rev(qs - qs0, qseq);
	} else rs1 = rs, qs1 = qs;
	re1 = rs, qe1 = qs;
	assert(qs1 >= 0 && rs1 >= 0);

	for (i = is_sr? cnt1 - 1 : 1; i < cnt1; ++i) { // gap filling
		if ((a[as1+i].y & (MM_SEED_IGNORE|MM_SEED_TANDEM)) && i != cnt1 - 1) continue;
		if (is_sr && !(mi->flag & MM_I_HPC)) {
			re = (int32_t)a[as1 + i].x + 1;
			qe = (int32_t)a[as1 + i].y + 1;
		} else mm_adjust_minier(mi, qseq0, &a[as1 + i], &re, &qe);
		re1 = re, qe1 = qe;
		if (i == cnt1 - 1 || (a[as1+i].y&MM_SEED_LONG_JOIN) || (qe - qs >= opt->min_ksw_len && re - rs >= opt->min_ksw_len)) {
			int j, bw1 = bw_long, zdrop_code;
			if (a[as1+i].y & MM_SEED_LONG_JOIN)
				bw1 = qe - qs > re - rs? qe - qs : re - rs;
			// perform alignment
			if (opt->flag & MM_F_QSTRAND) {
				qseq = &qseq0[0][qs];
				mm_idx_getseq2(mi, rev, rid, rs, re, tseq);
			} else {
				qseq = &qseq0[rev][qs];
				mm_idx_getseq(mi, rid, rs, re, tseq);
			}
			mm_idx_bed_junc(mi, rid, rs, re, junc);
			if (is_sr) { // perform ungapped alignment
				assert(qe - qs == re - rs);
				ksw_reset_extz(ez);
				for (j = 0, ez->score = 0; j < qe - qs; ++j) {
					if (qseq[j] >= 4 || tseq[j] >= 4) ez->score += opt->e2;
					else ez->score += qseq[j] == tseq[j]? opt->a : -opt->b;
				}
				ez->cigar = ksw_push_cigar(km, &ez->n_cigar, &ez->m_cigar, ez->cigar, MM_CIGAR_MATCH, qe - qs);
			} else { // perform normal gapped alignment
				mm_align_pair(km, opt, qe - qs, qseq, re - rs, tseq, junc, mat, bw1, -1, opt->zdrop, extra_flag|KSW_EZ_APPROX_MAX, ez); // first pass: with approximate Z-drop
			}
			// test Z-drop and inversion Z-drop
			if ((zdrop_code = mm_test_zdrop(km, opt, qseq, tseq, ez->n_cigar, ez->cigar, mat)) != 0)
				mm_align_pair(km, opt, qe - qs, qseq, re - rs, tseq, junc, mat, bw1, -1, zdrop_code == 2? opt->zdrop_inv : opt->zdrop, extra_flag, ez); // second pass: lift approximate
			// update CIGAR
			if (ez->n_cigar > 0)
				mm_append_cigar(r, ez->n_cigar, ez->cigar);
			if (ez->zdropped) { // truncated by Z-drop; TODO: sometimes Z-drop kicks in because the next seed placement is wrong. This can be fixed in principle.
				if (!r->p) {
					assert(ez->n_cigar == 0);
					uint32_t capacity = sizeof(mm_extra_t)/4;
					kroundup32(capacity);
					r->p = (mm_extra_t*)calloc(capacity, 4);
					r->p->capacity = capacity;
				}
				for (j = i - 1; j >= 0; --j)
					if ((int32_t)a[as1 + j].x <= rs + ez->max_t)
						break;
				dropped = 1;
				if (j < 0) j = 0;
				r->p->dp_score += ez->max;
				re1 = rs + (ez->max_t + 1);
				qe1 = qs + (ez->max_q + 1);
				if (cnt1 - (j + 1) >= opt->min_cnt) {
					mm_split_reg(r, r2, as1 + j + 1 - r->as, qlen, a, !!(opt->flag&MM_F_QSTRAND));
					if (zdrop_code == 2) r2->split_inv = 1;
				}
				break;
			} else r->p->dp_score += ez->score;
			rs = re, qs = qe;
		}
	}

	if (!dropped && qe < qe0 && re < re0) { // right extension
		if (opt->flag & MM_F_QSTRAND) {
			qseq = &qseq0[0][qe];
			mm_idx_getseq2(mi, rev, rid, re, re0, tseq);
		} else {
			qseq = &qseq0[rev][qe];
			mm_idx_getseq(mi, rid, re, re0, tseq);
		}
		mm_idx_bed_junc(mi, rid, re, re0, junc);
		mm_align_pair(km, opt, qe0 - qe, qseq, re0 - re, tseq, junc, mat, bw, opt->end_bonus, opt->zdrop, extra_flag|KSW_EZ_EXTZ_ONLY, ez);
		if (ez->n_cigar > 0) {
			mm_append_cigar(r, ez->n_cigar, ez->cigar);
			r->p->dp_score += ez->max;
		}
		re1 = re + (ez->reach_end? ez->mqe_t + 1 : ez->max_t + 1);
		qe1 = qe + (ez->reach_end? qe0 - qe : ez->max_q + 1);
	}
	assert(qe1 <= qlen);

	r->rs = rs1, r->re = re1;
	if (!rev || (opt->flag & MM_F_QSTRAND)) r->qs = qs1, r->qe = qe1;
	else r->qs = qlen - qe1, r->qe = qlen - qs1;

	assert(re1 - rs1 <= re0 - rs0);
	if (r->p) {
		if (opt->flag & MM_F_QSTRAND) {
			mm_idx_getseq2(mi, r->rev, rid, rs1, re1, tseq);
			qseq = &qseq0[0][qs1];
		} else {
			mm_idx_getseq(mi, rid, rs1, re1, tseq);
			qseq = &qseq0[r->rev][qs1];
		}
		mm_update_extra(r, qseq, tseq, mat, opt->q, opt->e, opt->flag & MM_F_EQX, !(opt->flag & MM_F_SR));
		if (rev && r->p->trans_strand)
			r->p->trans_strand ^= 3; // flip to the read strand
	}

	kfree(km, tseq);
	kfree(km, junc);
}

static int mm_align1_inv(void *km, const mm_mapopt_t *opt, const mm_idx_t *mi, int qlen, uint8_t *qseq0[2], const mm_reg1_t *r1, const mm_reg1_t *r2, mm_reg1_t *r_inv, ksw_extz_t *ez)
{ // NB: this doesn't work with the qstrand mode
	int tl, ql, score, ret = 0, q_off, t_off;
	uint8_t *tseq, *qseq;
	int8_t mat[25];
	void *qp;

	memset(r_inv, 0, sizeof(mm_reg1_t));
	if (!(r1->split&1) || !(r2->split&2)) return 0;
	if (r1->id != r1->parent && r1->parent != MM_PARENT_TMP_PRI) return 0;
	if (r2->id != r2->parent && r2->parent != MM_PARENT_TMP_PRI) return 0;
	if (r1->rid != r2->rid || r1->rev != r2->rev) return 0;
	ql = r1->rev? r1->qs - r2->qe : r2->qs - r1->qe;
	tl = r2->rs - r1->re;
	if (ql < opt->min_chain_score || ql > opt->max_gap) return 0;
	if (tl < opt->min_chain_score || tl > opt->max_gap) return 0;

	ksw_gen_ts_mat(5, mat, opt->a, opt->b, opt->transition, opt->sc_ambi);
	tseq = (uint8_t*)kmalloc(km, tl);
	mm_idx_getseq(mi, r1->rid, r1->re, r2->rs, tseq);
	qseq = r1->rev? &qseq0[0][r2->qe] : &qseq0[1][qlen - r2->qs];

	mm_seq_rev(ql, qseq);
	mm_seq_rev(tl, tseq);
	qp = ksw_ll_qinit(km, 2, ql, qseq, 5, mat);
	score = ksw_ll_i16(qp, tl, tseq, opt->q, opt->e, &q_off, &t_off);
	kfree(km, qp);
	mm_seq_rev(ql, qseq);
	mm_seq_rev(tl, tseq);
	if (score < opt->min_dp_max) goto end_align1_inv;
	q_off = ql - (q_off + 1), t_off = tl - (t_off + 1);
	mm_align_pair(km, opt, ql - q_off, qseq + q_off, tl - t_off, tseq + t_off, 0, mat, (int)(opt->bw * 1.5), -1, opt->zdrop, KSW_EZ_EXTZ_ONLY, ez);
	if (ez->n_cigar == 0) goto end_align1_inv; // should never be here
	mm_append_cigar(r_inv, ez->n_cigar, ez->cigar);
	r_inv->p->dp_score = ez->max;
	r_inv->id = -1;
	r_inv->parent = MM_PARENT_UNSET;
	r_inv->inv = 1;
	r_inv->rev = !r1->rev;
	r_inv->rid = r1->rid;
	r_inv->div = -1.0f;
	if (r_inv->rev == 0) {
		r_inv->qs = r2->qe + q_off;
		r_inv->qe = r_inv->qs + ez->max_q + 1;
	} else {
		r_inv->qe = r2->qs - q_off;
		r_inv->qs = r_inv->qe - (ez->max_q + 1);
	}
	r_inv->rs = r1->re + t_off;
	r_inv->re = r_inv->rs + ez->max_t + 1;
	mm_update_extra(r_inv, &qseq[q_off], &tseq[t_off], mat, opt->q, opt->e, opt->flag & MM_F_EQX, !(opt->flag & MM_F_SR));
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

static inline void mm_count_gaps(const mm_reg1_t *r, int32_t *n_gap_, int32_t *n_gapo_)
{
	uint32_t i;
	int32_t n_gapo = 0, n_gap = 0;
	*n_gap_ = *n_gapo_ = -1;
	if (r->p == 0) return;
	for (i = 0; i < r->p->n_cigar; ++i) {
		int32_t op = r->p->cigar[i] & 0xf, len = r->p->cigar[i] >> 4;
		if (op == MM_CIGAR_INS || op == MM_CIGAR_DEL)
			++n_gapo, n_gap += len;
	}
	*n_gap_ = n_gap, *n_gapo_ = n_gapo;
}

double mm_event_identity(const mm_reg1_t *r)
{
	int32_t n_gap, n_gapo;
	if (r->p == 0) return -1.0f;
	mm_count_gaps(r, &n_gap, &n_gapo);
	return (double)r->mlen / (r->blen + r->p->n_ambi - n_gap + n_gapo);
}

static int32_t mm_recal_max_dp(const mm_reg1_t *r, double b2, int32_t match_sc)
{
	uint32_t i;
	int32_t n_gap = 0, n_gapo = 0, n_mis;
	double gap_cost = 0.0;
	if (r->p == 0) return -1;
	for (i = 0; i < r->p->n_cigar; ++i) {
		int32_t op = r->p->cigar[i] & 0xf, len = r->p->cigar[i] >> 4;
		if (op == MM_CIGAR_INS || op == MM_CIGAR_DEL) {
			gap_cost += b2 + (double)mg_log2(1.0 + len);
			++n_gapo, n_gap += len;
		}
	}
	n_mis = r->blen + r->p->n_ambi - r->mlen - n_gap;
	return (int32_t)(match_sc * (r->mlen - b2 * n_mis - gap_cost) + .499);
}

void mm_update_dp_max(int qlen, int n_regs, mm_reg1_t *regs, float frac, int a, int b)
{
	int32_t max = -1, max2 = -1, i, max_i = -1;
	double div, b2;
	if (n_regs < 2) return;
	for (i = 0; i < n_regs; ++i) {
		mm_reg1_t *r = &regs[i];
		if (r->p == 0) continue;
		if (r->p->dp_max > max) max2 = max, max = r->p->dp_max, max_i = i;
		else if (r->p->dp_max > max2) max2 = r->p->dp_max;
	}
	if (max_i < 0 || max < 0 || max2 < 0) return;
	if (regs[max_i].qe - regs[max_i].qs < (double)qlen * frac) return;
	if (max2 < (double)max * frac) return;
	div = 1. - mm_event_identity(&regs[max_i]);
	if (div < 0.02) div = 0.02;
	b2 = 0.5 / div; // max value: 25
	if (b2 * a < b) b2 = (double)a / b;
	for (i = 0; i < n_regs; ++i) {
		mm_reg1_t *r = &regs[i];
		if (r->p == 0) continue;
		r->p->dp_max = mm_recal_max_dp(r, b2, a);
		if (r->p->dp_max < 0) r->p->dp_max = 0;
	}
}

mm_reg1_t *mm_align_skeleton(void *km, const mm_mapopt_t *opt, const mm_idx_t *mi, int qlen, const char *qstr, int *n_regs_, mm_reg1_t *regs, mm128_t *a)
{
	extern unsigned char seq_nt4_table[256];
	int32_t i, n_regs = *n_regs_, n_a;
	uint8_t *qseq0[2];
	ksw_extz_t ez;

	// encode the query sequence
	qseq0[0] = (uint8_t*)kmalloc(km, qlen * 2);
	qseq0[1] = qseq0[0] + qlen;
	for (i = 0; i < qlen; ++i) {
		qseq0[0][i] = seq_nt4_table[(uint8_t)qstr[i]];
		qseq0[1][qlen - 1 - i] = qseq0[0][i] < 4? 3 - qseq0[0][i] : 4;
	}

	// align through seed hits
	n_a = mm_squeeze_a(km, n_regs, regs, a);
	memset(&ez, 0, sizeof(ksw_extz_t));
	for (i = 0; i < n_regs; ++i) {
		mm_reg1_t r2;
		if ((opt->flag&MM_F_SPLICE) && (opt->flag&MM_F_SPLICE_FOR) && (opt->flag&MM_F_SPLICE_REV)) { // then do two rounds of alignments for both strands
			mm_reg1_t s[2], s2[2];
			int which, trans_strand;
			s[0] = s[1] = regs[i];
			mm_align1(km, opt, mi, qlen, qseq0, &s[0], &s2[0], n_a, a, &ez, MM_F_SPLICE_FOR);
			mm_align1(km, opt, mi, qlen, qseq0, &s[1], &s2[1], n_a, a, &ez, MM_F_SPLICE_REV);
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
		} else { // one round of alignment
			mm_align1(km, opt, mi, qlen, qseq0, &regs[i], &r2, n_a, a, &ez, opt->flag);
			if (opt->flag&MM_F_SPLICE)
				regs[i].p->trans_strand = opt->flag&MM_F_SPLICE_FOR? 1 : 2;
		}
		if (r2.cnt > 0) regs = mm_insert_reg(&r2, i, &n_regs, regs);
		if (i > 0 && regs[i].split_inv && !(opt->flag & MM_F_NO_INV)) {
			if (mm_align1_inv(km, opt, mi, qlen, qseq0, &regs[i-1], &regs[i], &r2, &ez)) {
				regs = mm_insert_reg(&r2, i, &n_regs, regs);
				++i; // skip the inserted INV alignment
			}
		}
	}
	*n_regs_ = n_regs;
	kfree(km, qseq0[0]);
	kfree(km, ez.cigar);
	mm_filter_regs(opt, qlen, n_regs_, regs);
	if (!(opt->flag&MM_F_SR) && !opt->split_prefix && qlen >= opt->rank_min_len) {
		mm_update_dp_max(qlen, *n_regs_, regs, opt->rank_frac, opt->a, opt->b);
		mm_filter_regs(opt, qlen, n_regs_, regs);
	}
	mm_hit_sort(km, n_regs_, regs, opt->alt_drop);
	return regs;
}
