#include <stdio.h>
#include "mmpriv.h"
#include "kalloc.h"

#define MM_MIN_EXON_LEN 20

static int32_t mm_jump_check(void *km, const mm_idx_t *mi, int32_t qlen, const uint8_t *qseq0, const mm_reg1_t *r, int32_t ext, int32_t is_left) // TODO: check close N
{
	int32_t clip, clen, e = !r->rev ^ !is_left; // 0 for left of the alignment; 1 for right
	uint32_t cigar;
	if (!r->p || r->p->n_cigar <= 0) return -1; // only working with CIGAR
	clip = e == 0? r->qs : qlen - r->qe;
	cigar = r->p->cigar[is_left? 0 : r->p->n_cigar - 1];
	clen = (cigar&0xf) == MM_CIGAR_MATCH? cigar>>4 : 0;
	if (clen <= ext) return -1;
	if (is_left) {
		if (clip >= r->rs) return -1; // no space to jump
	} else {
		if (clip >= mi->seq[r->rid].len - r->re) return -1; // no space to jump
	}
	return 0;
}

static uint8_t *mm_jump_get_qseq_seq(void *km, int32_t qlen, const uint8_t *qseq0, const mm_reg1_t *r, int32_t is_left, int32_t ql0, uint8_t *qseq)
{
	extern unsigned char seq_nt4_table[256];
	int32_t i, k = 0;
	if (!r->rev) {
		if (is_left)
			for (i = 0; i < ql0; ++i)
				qseq[k++] = seq_nt4_table[(uint8_t)qseq0[i]];
		else
			for (i = qlen - ql0; i < qlen; ++i)
				qseq[k++] = seq_nt4_table[(uint8_t)qseq0[i]];
	} else {
		if (is_left)
			for (i = qlen - 1; i >= qlen - ql0; --i) {
				uint8_t c = seq_nt4_table[(uint8_t)qseq0[i]];
				qseq[k++] = c >= 4? c : 3 - c;
			}
		else
			for (i = ql0 - 1; i >= 0; --i) {
				uint8_t c = seq_nt4_table[(uint8_t)qseq0[i]];
				qseq[k++] = c >= 4? c : 3 - c;
			}
	}
	return qseq;
}

static void mm_jump_split_left(void *km, const mm_idx_t *mi, const mm_mapopt_t *opt, int32_t qlen, const uint8_t *qseq0, mm_reg1_t *r, int32_t ts_strand)
{
	uint8_t *tseq = 0, *qseq = 0;
	int32_t i, n, l, i0, m, mm0;
	int32_t i0_anno = -1, n_anno = 0, mm0_anno = 0, i0_misc = -1, n_misc = 0, mm0_misc = 0;
	int32_t ext = 1 + (opt->b + opt->a - 1) / opt->a + 1;
	int32_t clip = !r->rev? r->qs : qlen - r->qe;
	int32_t extt = clip < ext? clip : ext;
	const mm_idx_jjump1_t *a;

	if (mm_jump_check(km, mi, qlen, qseq0, r, ext + MM_MIN_EXON_LEN, 1) < 0) return;
	a = mm_idx_jump_get(mi, r->rid, r->rs - extt, r->rs + ext, &n);
	if (n == 0) return;

	for (i = 0; i < n; ++i) { // traverse possible jumps
		const mm_idx_jjump1_t *ai = &a[i];
		int32_t tlen, tl1, j, mm1, mm2;
		assert(ai->off >= r->rs - extt && ai->off <= r->rs + ext);
		if (ts_strand * ai->strand < 0) continue; // wrong strand
		if (ai->off2 >= ai->off) continue; // wrong direction
		if (ai->off - ai->off2 < 6) continue; // intron too small
		if (ai->off2 < clip + ext) continue; // not long enough
		if (tseq == 0) {
			tseq = Kcalloc(km, uint8_t, (clip + ext) * 2); // tseq and qseq are allocated together
			qseq = tseq + clip + ext;
			mm_jump_get_qseq_seq(km, qlen, qseq0, r, 1, clip + ext, qseq);
		}
		tl1 = clip + (ai->off - r->rs);
		tlen = mm_idx_getseq2(mi, 0, r->rid, ai->off, r->rs + ext, &tseq[tl1]);
		assert(tlen == r->rs + ext - ai->off);
		tlen = mm_idx_getseq2(mi, 0, r->rid, ai->off2 - tl1, ai->off2, tseq);
		assert(tlen == tl1);
		for (j = 0, mm1 = 0; j < tl1; ++j)
			if (qseq[j] != tseq[j] || qseq[j] > 3 || tseq[j] > 3)
				++mm1;
		for (mm2 = 0; j < clip + ext; ++j)
			if (qseq[j] != tseq[j] || qseq[j] > 3 || tseq[j] > 3)
				++mm2;
		if (mm1 == 0 && mm2 <= 1) {
			if (ai->flag & MM_JUNC_ANNO)
				i0_anno = i, mm0_anno = mm1 + mm2, ++n_anno; // i0 points to the rightmost i
			else
				i0_misc = i, mm0_misc = mm1 + mm2, ++n_misc;
		}
	}
	if (n_anno > 0) m = n_anno, i0 = i0_anno, mm0 = mm0_anno;
	else m = n_misc, i0 = i0_misc, mm0 = mm0_misc;
	kfree(km, tseq);

	l = m > 0? a[i0].off - r->rs : 0; // may be negative
	if (m == 1 && clip + l >= opt->jump_min_match) { // add one more exon
		mm_enlarge_cigar(r, 2);
		memmove(r->p->cigar + 2, r->p->cigar, r->p->n_cigar * 4);
		r->p->cigar[0] = (clip + l) << 4 | MM_CIGAR_MATCH;
		r->p->cigar[1] = (a[i0].off - a[i0].off2) << 4 | MM_CIGAR_N_SKIP;
		r->p->cigar[2] = ((r->p->cigar[2]>>4) - l) << 4 | MM_CIGAR_MATCH;
		r->p->n_cigar += 2;
		r->rs = a[i0].off2 - (clip + l);
		if (!r->rev) r->qs = 0;
		else r->qe = qlen;
		r->blen += clip, r->mlen += clip - mm0;
		r->p->dp_max0 += (clip - mm0) * opt->a - mm0 * opt->b;
		r->p->dp_max += (clip - mm0) * opt->a - mm0 * opt->b;
		if (!r->is_spliced) r->is_spliced = 1, r->p->dp_max += (opt->a + opt->b) + ((opt->a + opt->b) >> 1);
	} else if (m > 0 && a[i0].off > r->rs) { // trim by l; l is always positive
		r->p->cigar[0] -= l << 4 | MM_CIGAR_MATCH;
		r->rs += l;
		if (!r->rev) r->qs += l;
		else r->qe -= l;
	}
}

static void mm_jump_split_right(void *km, const mm_idx_t *mi, const mm_mapopt_t *opt, int32_t qlen, const uint8_t *qseq0, mm_reg1_t *r, int32_t ts_strand)
{
	uint8_t *tseq = 0, *qseq = 0;
	int32_t i, n, l, i0, m, mm0;
	int32_t i0_anno = -1, n_anno = 0, mm0_anno = 0, i0_misc = -1, n_misc = 0, mm0_misc = 0;
	int32_t ext = 1 + (opt->b + opt->a - 1) / opt->a + 1;
	int32_t clip = !r->rev? qlen - r->qe : r->qs;
	int32_t extt = clip < ext? clip : ext;
	const mm_idx_jjump1_t *a;

	if (mm_jump_check(km, mi, qlen, qseq0, r, ext + MM_MIN_EXON_LEN, 0) < 0) return;
	a = mm_idx_jump_get(mi, r->rid, r->re - ext, r->re + extt, &n);
	if (n == 0) return;

	for (i = 0; i < n; ++i) { // traverse possible jumps
		const mm_idx_jjump1_t *ai = &a[i];
		int32_t tlen, tl1, j, mm1, mm2;
		assert(ai->off >= r->re - ext && ai->off <= r->re + extt);
		if (ts_strand * ai->strand < 0) continue; // wrong strand
		if (ai->off2 <= ai->off) continue; // wrong direction
		if (ai->off2 - ai->off < 6) continue; // intron too small
		if (ai->off2 + clip + ext > mi->seq[r->rid].len) continue; // not long enough
		if (tseq == 0) {
			tseq = Kcalloc(km, uint8_t, (clip + ext) * 2); // tseq and qseq are allocated together
			qseq = tseq + clip + ext;
			mm_jump_get_qseq_seq(km, qlen, qseq0, r, 0, clip + ext, qseq);
		}
		tl1 = clip + (r->re - ai->off);
		tlen = mm_idx_getseq2(mi, 0, r->rid, r->re - ext, ai->off, tseq);
		assert(tlen == ai->off - (r->re - ext));
		tlen = mm_idx_getseq2(mi, 0, r->rid, ai->off2, ai->off2 + tl1, &tseq[clip + ext - tl1]);
		assert(tlen == tl1);
		for (j = 0, mm2 = 0; j < clip + ext - tl1; ++j)
			if (qseq[j] != tseq[j] || qseq[j] > 3 || tseq[j] > 3)
				++mm2;
		for (mm1 = 0; j < clip + ext; ++j)
			if (qseq[j] != tseq[j] || qseq[j] > 3 || tseq[j] > 3)
				++mm1;
		if (mm1 == 0 && mm2 <= 1) {
			if (ai->flag & MM_JUNC_ANNO) {
				if (i0_anno < 0) i0_anno = i, mm0_anno = mm1 + mm2;
				++n_anno;
			} else {
				if (i0_misc < 0) i0_misc = i, mm0_misc = mm1 + mm2;
				++n_misc;
			}
		}
	}
	if (n_anno > 0) m = n_anno, i0 = i0_anno, mm0 = mm0_anno;
	else m = n_misc, i0 = i0_misc, mm0 = mm0_misc;
	kfree(km, tseq);

	l = m > 0? r->re - a[i0].off : 0; // may be negative
	if (m == 1 && clip + l >= opt->jump_min_match) { // add one more exon
		mm_enlarge_cigar(r, 2);
		r->p->cigar[r->p->n_cigar - 1] = ((r->p->cigar[r->p->n_cigar - 1]>>4) - l) << 4 | MM_CIGAR_MATCH;
		r->p->cigar[r->p->n_cigar] = (a[i0].off2 - a[i0].off) << 4 | MM_CIGAR_N_SKIP;
		r->p->cigar[r->p->n_cigar + 1] = (clip + l) << 4 | MM_CIGAR_MATCH;
		r->p->n_cigar += 2;
		r->re = a[i0].off2 + (clip + l);
		if (!r->rev) r->qe = qlen;
		else r->qs = 0;
		r->blen += clip, r->mlen += clip - mm0;
		r->p->dp_max0 += (clip - mm0) * opt->a - mm0 * opt->b;
		r->p->dp_max += (clip - mm0) * opt->a - mm0 * opt->b;
		if (!r->is_spliced) r->is_spliced = 1, r->p->dp_max += (opt->a + opt->b) + ((opt->a + opt->b) >> 1);
	} else if (m > 0 && r->re > a[i0].off) { // trim by l; l is always positive
		r->p->cigar[r->p->n_cigar - 1] -= l << 4 | MM_CIGAR_MATCH;
		r->re -= l;
		if (!r->rev) r->qe -= l;
		else r->qs += l;
	}
}

void mm_jump_split(void *km, const mm_idx_t *mi, const mm_mapopt_t *opt, int32_t qlen, const uint8_t *qseq, mm_reg1_t *r, int32_t ts_strand)
{
	assert((opt->flag & MM_F_EQX) == 0);
	mm_jump_split_left(km,  mi, opt, qlen, qseq, r, ts_strand);
	mm_jump_split_right(km, mi, opt, qlen, qseq, r, ts_strand);
}
