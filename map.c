#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "kthread.h"
#include "kvec.h"
#include "kalloc.h"
#include "sdust.h"
#include "mmpriv.h"
#include "bseq.h"
#include "khash.h"

void mm_mapopt_init(mm_mapopt_t *opt)
{
	memset(opt, 0, sizeof(mm_mapopt_t));
	opt->seed = 11;
	opt->mid_occ_frac = 2e-4f;
	opt->sdust_thres = 0; // no SDUST masking

	opt->min_cnt = 3;
	opt->min_chain_score = 40;
	opt->bw = 500;
	opt->max_gap = 5000;
	opt->max_gap_ref = -1;
	opt->max_chain_skip = 25;

	opt->mask_level = 0.5f;
	opt->pri_ratio = 0.8f;
	opt->best_n = 5;

	opt->max_join_long = 20000;
	opt->max_join_short = 2000;
	opt->min_join_flank_sc = 1000;

	opt->a = 2, opt->b = 4, opt->q = 4, opt->e = 2, opt->q2 = 24, opt->e2 = 1;
	opt->zdrop = 400;
	opt->end_bonus = -1;
	opt->min_dp_max = opt->min_chain_score * opt->a;
	opt->min_ksw_len = 200;
	opt->mini_batch_size = 500000000;

	opt->pe_ori = 0; // FF
	opt->pe_bonus = 33;
}

void mm_mapopt_update(mm_mapopt_t *opt, const mm_idx_t *mi)
{
	if ((opt->flag & MM_F_SPLICE_FOR) && (opt->flag & MM_F_SPLICE_REV))
		opt->flag |= MM_F_SPLICE;
	if (opt->mid_occ <= 0)
		opt->mid_occ = mm_idx_cal_max_occ(mi, opt->mid_occ_frac);
	if (mm_verbose >= 3)
		fprintf(stderr, "[M::%s::%.3f*%.2f] mid_occ = %d\n", __func__, realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0), opt->mid_occ);
}

void mm_mapopt_max_intron_len(mm_mapopt_t *opt, int max_intron_len)
{
	if ((opt->flag & MM_F_SPLICE) && max_intron_len > 0)
		opt->max_gap_ref = opt->bw = max_intron_len;
}

int mm_set_opt(const char *preset, mm_idxopt_t *io, mm_mapopt_t *mo)
{
	if (preset == 0) {
		mm_idxopt_init(io);
		mm_mapopt_init(mo);
	} else if (strcmp(preset, "ava-ont") == 0) {
		io->is_hpc = 0, io->k = 15, io->w = 5;
		mo->flag |= MM_F_AVA | MM_F_NO_SELF;
		mo->min_chain_score = 100, mo->pri_ratio = 0.0f, mo->max_gap = 10000, mo->max_chain_skip = 25;
	} else if (strcmp(preset, "ava-pb") == 0) {
		io->is_hpc = 1, io->k = 19, io->w = 5;
		mo->flag |= MM_F_AVA | MM_F_NO_SELF;
		mo->min_chain_score = 100, mo->pri_ratio = 0.0f, mo->max_gap = 10000, mo->max_chain_skip = 25;
	} else if (strcmp(preset, "map10k") == 0 || strcmp(preset, "map-pb") == 0) {
		io->is_hpc = 1, io->k = 19;
	} else if (strcmp(preset, "map-ont") == 0) {
		io->is_hpc = 0, io->k = 15;
	} else if (strcmp(preset, "asm5") == 0) {
		io->is_hpc = 0, io->k = 19, io->w = 19;
		mo->a = 1, mo->b = 19, mo->q = 39, mo->q2 = 81, mo->e = 3, mo->e2 = 1, mo->zdrop = 200;
		mo->min_dp_max = 200;
		mo->best_n = 50;
	} else if (strcmp(preset, "asm10") == 0) {
		io->is_hpc = 0, io->k = 19, io->w = 19;
		mo->a = 1, mo->b = 9, mo->q = 16, mo->q2 = 41, mo->e = 2, mo->e2 = 1, mo->zdrop = 200;
		mo->min_dp_max = 200;
		mo->best_n = 50;
	} else if (strcmp(preset, "short") == 0 || strcmp(preset, "sr") == 0) {
		io->is_hpc = 0, io->k = 21, io->w = 11;
		mo->flag |= MM_F_SR | MM_F_FRAG_MODE | MM_F_NO_PRINT_2ND | MM_F_2_IO_THREADS;
		mo->pe_ori = 0<<1|1; // FR
		mo->a = 2, mo->b = 8, mo->q = 12, mo->e = 2, mo->q2 = 24, mo->e2 = 1;
		mo->zdrop = 100;
		mo->end_bonus = 10;
		mo->max_frag_len = 800;
		mo->max_gap = 100;
		mo->bw = 100;
		mo->pri_ratio = 0.5f;
		mo->min_cnt = 2;
		mo->min_chain_score = 25;
		mo->min_dp_max = 40;
		mo->best_n = 20;
		mo->mid_occ = 1000;
		mo->max_occ = 5000;
		mo->mini_batch_size = 50000000;
	} else if (strcmp(preset, "splice") == 0 || strcmp(preset, "cdna") == 0) {
		io->is_hpc = 0, io->k = 15, io->w = 5;
		mo->flag |= MM_F_SPLICE | MM_F_SPLICE_FOR | MM_F_SPLICE_REV | MM_F_SPLICE_FLANK;
		mo->max_gap = 2000, mo->max_gap_ref = mo->bw = 200000;
		mo->a = 1, mo->b = 2, mo->q = 2, mo->e = 1, mo->q2 = 32, mo->e2 = 0;
		mo->noncan = 9;
		mo->zdrop = 200;
	} else return -1;
	return 0;
}

typedef struct {
	uint32_t n;
	uint32_t qpos;
	uint32_t seg_id;
	const uint64_t *cr;
} mm_match_t;

struct mm_tbuf_s {
	sdust_buf_t *sdb;
	mm128_v mini;
	void *km;
};

mm_tbuf_t *mm_tbuf_init(void)
{
	mm_tbuf_t *b;
	b = (mm_tbuf_t*)calloc(1, sizeof(mm_tbuf_t));
	if (!(mm_dbg_flag & 1)) b->km = km_init();
	b->sdb = sdust_buf_init(b->km);
	return b;
}

void mm_tbuf_destroy(mm_tbuf_t *b)
{
	if (b == 0) return;
	kfree(b->km, b->mini.a);
	sdust_buf_destroy(b->sdb);
	km_destroy(b->km);
	free(b);
}

static int mm_dust_minier(int n, mm128_t *a, int l_seq, const char *seq, int sdust_thres, sdust_buf_t *sdb)
{
	int n_dreg, j, k, u = 0;
	const uint64_t *dreg;
	if (sdust_thres <= 0 || sdb == 0) return n;
	dreg = sdust_core((const uint8_t*)seq, l_seq, sdust_thres, 64, &n_dreg, sdb);
	for (j = k = 0; j < n; ++j) { // squeeze out minimizers that significantly overlap with LCRs
		int32_t qpos = (uint32_t)a[j].y>>1, span = a[j].x&0xff;
		int32_t s = qpos - (span - 1), e = s + span;
		while (u < n_dreg && (uint32_t)dreg[u] <= s) ++u;
		if (u < n_dreg && dreg[u]>>32 < e) {
			int v, l = 0;
			for (v = u; v < n_dreg && dreg[v]>>32 < e; ++v) { // iterate over LCRs overlapping this minimizer
				int ss = s > dreg[v]>>32? s : dreg[v]>>32;
				int ee = e < (uint32_t)dreg[v]? e : (uint32_t)dreg[v];
				l += ee - ss;
			}
			if (l <= span>>1) a[k++] = a[j]; // keep the minimizer if less than half of it falls in masked region
		}
	}
	return k; // the new size
}

static void collect_minimizers(const mm_mapopt_t *opt, const mm_idx_t *mi, int n_segs, const int *qlens, const char **seqs, mm_tbuf_t *b)
{
	int i, j, n, sum = 0;
	b->mini.n = 0;
	for (i = n = 0; i < n_segs; ++i) {
		mm_sketch(b->km, seqs[i], qlens[i], mi->w, mi->k, i, mi->is_hpc, &b->mini);
		for (j = n; j < b->mini.n; ++j)
			b->mini.a[j].y += sum << 1;
		if (opt->sdust_thres > 0) // mask low-complexity minimizers
			b->mini.n = n + mm_dust_minier(b->mini.n - n, b->mini.a + n, qlens[i], seqs[i], opt->sdust_thres, b->sdb);
		sum += qlens[i], n = b->mini.n;
	}
}

static mm128_t *collect_seed_hits(const mm_mapopt_t *opt, int max_occ, const mm_idx_t *mi, const char *qname, int qlen, int64_t *n_a, int *rep_len, mm_tbuf_t *b)
{
	int rep_st = 0, rep_en = 0, i;
	mm_match_t *m;
	mm128_t *a;

	m = (mm_match_t*)kmalloc(b->km, b->mini.n * sizeof(mm_match_t));
	for (i = 0; i < b->mini.n; ++i) {
		int t;
		mm128_t *p = &b->mini.a[i];
		m[i].qpos = (uint32_t)p->y;
		m[i].cr = mm_idx_get(mi, p->x>>8, &t);
		m[i].n = t;
		m[i].seg_id = p->y >> 32;
	}
	for (i = 0, *n_a = 0; i < b->mini.n; ++i) // find the length of a[]
		if (m[i].n < max_occ) *n_a += m[i].n;
	a = (mm128_t*)kmalloc(b->km, *n_a * sizeof(mm128_t));
	for (i = *rep_len = 0, *n_a = 0; i < b->mini.n; ++i) {
		mm128_t *p = &b->mini.a[i];
		mm_match_t *q = &m[i];
		const uint64_t *r = q->cr;
		int k, q_span = p->x & 0xff, is_tandem = 0;
		if (q->n >= max_occ) {
			int en = (q->qpos>>1) + 1, st = en - q_span;
			if (st > rep_en) {
				*rep_len += rep_en - rep_st;
				rep_st = st, rep_en = en;
			} else rep_en = en;
			continue;
		}
		if (i > 0 && p->x>>8 == b->mini.a[i - 1].x>>8) is_tandem = 1;
		if (i < b->mini.n - 1 && p->x>>8 == b->mini.a[i + 1].x>>8) is_tandem = 1;
		for (k = 0; k < q->n; ++k) {
			int32_t rpos = (uint32_t)r[k] >> 1;
			mm128_t *p;
			if (qname && (opt->flag&(MM_F_NO_SELF|MM_F_AVA))) {
				const char *tname = mi->seq[r[k]>>32].name;
				int cmp;
				cmp = strcmp(qname, tname);
				if ((opt->flag&MM_F_NO_SELF) && cmp == 0 && rpos == (q->qpos>>1)) // avoid the diagonal
					continue;
				if ((opt->flag&MM_F_AVA) && cmp > 0) // all-vs-all mode: map once
					continue;
			}
			p = &a[(*n_a)++];
			if ((r[k]&1) == (q->qpos&1)) { // forward strand
				p->x = (r[k]&0xffffffff00000000ULL) | rpos;
				p->y = (uint64_t)q_span << 32 | q->qpos >> 1;
			} else { // reverse strand
				p->x = 1ULL<<63 | (r[k]&0xffffffff00000000ULL) | rpos;
				p->y = (uint64_t)q_span << 32 | (qlen - ((q->qpos>>1) + 1 - q_span) - 1);
			}
			p->y |= (uint64_t)q->seg_id << MM_SEED_SEG_SHIFT;
			if (is_tandem) p->y |= MM_SEED_TANDEM;
		}
	}
	*rep_len += rep_en - rep_st;
	kfree(b->km, m);
	return a;
}

static void chain_post(const mm_mapopt_t *opt, int max_chain_gap_ref, const mm_idx_t *mi, void *km, int qlen, int n_segs, const int *qlens, int *n_regs, mm_reg1_t *regs, mm128_t *a)
{
	if (!(opt->flag & MM_F_AVA)) { // don't choose primary mapping(s) for read overlap
		mm_set_parent(km, opt->mask_level, *n_regs, regs, opt->a * 2 + opt->b);
		if (n_segs <= 1) mm_select_sub(km, opt->pri_ratio, mi->k*2, opt->best_n, n_regs, regs);
		else mm_select_sub_multi(km, opt->pri_ratio, 0.2f, 0.7f, max_chain_gap_ref, mi->k*2, opt->best_n, n_segs, qlens, n_regs, regs);
		if (!(opt->flag & MM_F_SPLICE) && !(opt->flag & MM_F_SR) && !(opt->flag & MM_F_NO_LJOIN))
			mm_join_long(km, opt, qlen, n_regs, regs, a);
	}
}

static mm_reg1_t *align_regs(const mm_mapopt_t *opt, const mm_idx_t *mi, void *km, int qlen, const char *seq, const char *qual, int *n_regs, mm_reg1_t *regs, mm128_t *a)
{
	if (!(opt->flag & MM_F_CIGAR)) return regs;
	regs = mm_align_skeleton(km, opt, mi, qlen, seq, qual, n_regs, regs, a); // this calls mm_filter_regs()
	if (!(opt->flag & MM_F_AVA)) {
		mm_set_parent(km, opt->mask_level, *n_regs, regs, opt->a * 2 + opt->b);
		mm_select_sub(km, opt->pri_ratio, mi->k*2, opt->best_n, n_regs, regs);
		mm_set_sam_pri(*n_regs, regs);
	}
	return regs;
}

void mm_map_frag(const mm_idx_t *mi, int n_segs, const int *qlens, const char **seqs, const char **quals, int *n_regs, mm_reg1_t **regs, mm_tbuf_t *b, const mm_mapopt_t *opt, const char *qname)
{
	int i, j, rep_len, qlen_sum, n_regs0;
	int max_chain_gap_qry, max_chain_gap_ref, is_splice = !!(opt->flag & MM_F_SPLICE), is_sr = !!(opt->flag & MM_F_SR);
	uint32_t hash;
	int64_t n_a;
	uint64_t *u;
	mm128_t *a;
	mm_reg1_t *regs0;

	for (i = 0, qlen_sum = 0; i < n_segs; ++i)
		qlen_sum += qlens[i], n_regs[i] = 0, regs[i] = 0;

	if (qlen_sum == 0 || n_segs <= 0 || n_segs > MM_MAX_SEG) return;

	hash  = qname? __ac_X31_hash_string(qname) : 0;
	hash ^= __ac_Wang_hash(qlen_sum) + __ac_Wang_hash(opt->seed);
	hash  = __ac_Wang_hash(hash);

	collect_minimizers(opt, mi, n_segs, qlens, seqs, b);
	a = collect_seed_hits(opt, opt->mid_occ, mi, qname, qlen_sum, &n_a, &rep_len, b);
	radix_sort_128x(a, a + n_a);

	if (mm_dbg_flag & MM_DBG_PRINT_SEED) {
		fprintf(stderr, "RS\t%d\n", rep_len);
		for (i = 0; i < n_a; ++i)
			fprintf(stderr, "SD\t%s\t%d\t%c\t%d\t%d\t%d\n", mi->seq[a[i].x<<1>>33].name, (int32_t)a[i].x, "+-"[a[i].x>>63], (int32_t)a[i].y, (int32_t)(a[i].y>>32&0xff),
					i == 0? 0 : ((int32_t)a[i].y - (int32_t)a[i-1].y) - ((int32_t)a[i].x - (int32_t)a[i-1].x));
	}

	// set max chaining gap on the query and the reference sequence
	if (is_sr)
		max_chain_gap_qry = qlen_sum > opt->max_gap? qlen_sum : opt->max_gap;
	else max_chain_gap_qry = opt->max_gap;
	if (opt->max_gap_ref > 0) {
		max_chain_gap_ref = opt->max_gap_ref; // always honor mm_mapopt_t::max_gap_ref if set
	} else if (opt->max_frag_len > 0) {
		max_chain_gap_ref = opt->max_frag_len - qlen_sum;
		if (max_chain_gap_ref < opt->max_gap) max_chain_gap_ref = opt->max_gap;
	} else max_chain_gap_ref = opt->max_gap;

	a = mm_chain_dp(max_chain_gap_ref, max_chain_gap_qry, opt->bw, opt->max_chain_skip, opt->min_cnt, opt->min_chain_score, is_splice, n_segs, n_a, a, &n_regs0, &u, b->km);

	if (opt->max_occ > opt->mid_occ && rep_len > 0) {
		int rechain = 0;
		if (n_regs0 > 0) { // test if the best chain has all the segments
			int n_chained_segs = 1, max = 0, max_i = -1, max_off = -1, off = 0;
			for (i = 0; i < n_regs0; ++i) { // find the best chain
				if (max < u[i]>>32) max = u[i]>>32, max_i = i, max_off = off;
				off += (uint32_t)u[i];
			}
			for (i = 1; i < (uint32_t)u[max_i]; ++i) // count the number of segments in the best chain
				if ((a[max_off+i].y&MM_SEED_SEG_MASK) != (a[max_off+i-1].y&MM_SEED_SEG_MASK))
					++n_chained_segs;
			if (n_chained_segs < n_segs)
				rechain = 1;
		} else rechain = 1;
		if (rechain) { // redo chaining with a higher max_occ threshold
			kfree(b->km, a);
			kfree(b->km, u);
			a = collect_seed_hits(opt, opt->max_occ, mi, qname, qlen_sum, &n_a, &rep_len, b);
			radix_sort_128x(a, a + n_a);
			a = mm_chain_dp(max_chain_gap_ref, max_chain_gap_qry, opt->bw, opt->max_chain_skip, opt->min_cnt, opt->min_chain_score, is_splice, n_segs, n_a, a, &n_regs0, &u, b->km);
		}
	}

	regs0 = mm_gen_regs(b->km, hash, qlen_sum, n_regs0, u, a);

	if (mm_dbg_flag & MM_DBG_PRINT_SEED)
		for (j = 0; j < n_regs0; ++j)
			for (i = regs0[j].as; i < regs0[j].as + regs0[j].cnt; ++i)
				fprintf(stderr, "CN\t%d\t%s\t%d\t%c\t%d\t%d\t%d\n", j, mi->seq[a[i].x<<1>>33].name, (int32_t)a[i].x, "+-"[a[i].x>>63], (int32_t)a[i].y, (int32_t)(a[i].y>>32&0xff),
						i == regs0[j].as? 0 : ((int32_t)a[i].y - (int32_t)a[i-1].y) - ((int32_t)a[i].x - (int32_t)a[i-1].x));

	chain_post(opt, max_chain_gap_ref, mi, b->km, qlen_sum, n_segs, qlens, &n_regs0, regs0, a);

	if (n_segs == 1) { // uni-segment
		regs0 = align_regs(opt, mi, b->km, qlens[0], seqs[0], quals? quals[0] : 0, &n_regs0, regs0, a);
		mm_set_mapq(n_regs0, regs0, opt->min_chain_score, opt->a, rep_len, is_sr);
		n_regs[0] = n_regs0, regs[0] = regs0;
	} else { // multi-segment
		mm_seg_t *seg;
		seg = mm_seg_gen(b->km, hash, n_segs, qlens, n_regs0, regs0, n_regs, regs, a); // split fragment chain to separate segment chains
		free(regs0);
		for (i = 0; i < n_segs; ++i) {
			mm_set_parent(b->km, opt->mask_level, n_regs[i], regs[i], opt->a * 2 + opt->b); // update mm_reg1_t::parent
			regs[i] = align_regs(opt, mi, b->km, qlens[i], seqs[i], quals? quals[i] : 0, &n_regs[i], regs[i], seg[i].a);
			mm_set_mapq(n_regs[i], regs[i], opt->min_chain_score, opt->a, rep_len, is_sr);
		}
		mm_seg_free(b->km, n_segs, seg);
		if (n_segs == 2 && opt->pe_ori >= 0 && (opt->flag&MM_F_CIGAR))
			mm_pair(b->km, max_chain_gap_ref, opt->pe_bonus, opt->a * 2 + opt->b, opt->a, qlens, n_regs, regs); // pairing
	}

	kfree(b->km, a);
	kfree(b->km, u);
}

mm_reg1_t *mm_map(const mm_idx_t *mi, int qlen, const char *seq, int *n_regs, mm_tbuf_t *b, const mm_mapopt_t *opt, const char *qname)
{
	mm_reg1_t *regs;
	mm_map_frag(mi, 1, &qlen, &seq, 0, n_regs, &regs, b, opt, qname);
	return regs;
}

/**************************
 * Multi-threaded mapping *
 **************************/

typedef struct {
	int mini_batch_size, n_processed, n_threads, n_fp;
	const mm_mapopt_t *opt;
	mm_bseq_file_t **fp;
	const mm_idx_t *mi;
	kstring_t str;
} pipeline_t;

typedef struct {
	const pipeline_t *p;
    int n_seq, n_frag;
	mm_bseq1_t *seq;
	int *n_reg, *seg_off, *n_seg;
	mm_reg1_t **reg;
	mm_tbuf_t **buf;
} step_t;

static void worker_for(void *_data, long i, int tid) // kt_for() callback
{
    step_t *s = (step_t*)_data;
	int *qlens, j, off = s->seg_off[i], pe_ori = s->p->opt->pe_ori, is_sr = !!(s->p->opt->flag & MM_F_SR);
	const char **qseqs, **quals = 0;
	mm_tbuf_t *b = s->buf[tid];
	if (mm_dbg_flag & MM_DBG_PRINT_QNAME)
		fprintf(stderr, "QR\t%s\t%d\n", s->seq[off].name, tid);
	qlens = (int*)kmalloc(b->km, s->n_seg[i] * sizeof(int));
	qseqs = (const char**)kmalloc(b->km, s->n_seg[i] * sizeof(const char**));
	quals = (const char**)kmalloc(b->km, s->n_seg[i] * sizeof(const char**));
	for (j = 0; j < s->n_seg[i]; ++j) {
		if (s->n_seg[i] == 2 && ((j == 0 && (pe_ori>>1&1)) || (j == 1 && (pe_ori&1))))
			mm_revcomp_bseq(&s->seq[off + j]);
		qlens[j] = s->seq[off + j].l_seq;
		qseqs[j] = s->seq[off + j].seq;
		quals[j] = is_sr? s->seq[off + j].qual : 0;
	}
	if (s->p->opt->flag & MM_F_INDEPEND_SEG) {
		for (j = 0; j < s->n_seg[i]; ++j)
			mm_map_frag(s->p->mi, 1, &qlens[j], &qseqs[j], &quals[j], &s->n_reg[off+j], &s->reg[off+j], b, s->p->opt, s->seq[off+j].name);
	} else {
		mm_map_frag(s->p->mi, s->n_seg[i], qlens, qseqs, quals, &s->n_reg[off], &s->reg[off], b, s->p->opt, s->seq[off].name);
	}
	for (j = 0; j < s->n_seg[i]; ++j) // flip the query strand and coordinate to the original read strand
		if (s->n_seg[i] == 2 && ((j == 0 && (pe_ori>>1&1)) || (j == 1 && (pe_ori&1)))) {
			int k, t;
			mm_revcomp_bseq(&s->seq[off + j]);
			for (k = 0; k < s->n_reg[off + j]; ++k) {
				mm_reg1_t *r = &s->reg[off + j][k];
				t = r->qs;
				r->qs = qlens[j] - r->qe;
				r->qe = qlens[j] - t;
				r->rev = !r->rev;
			}
		}
	kfree(b->km, qlens);
	kfree(b->km, qseqs);
	kfree(b->km, quals);
}

static void *worker_pipeline(void *shared, int step, void *in)
{
	int i, j, k;
    pipeline_t *p = (pipeline_t*)shared;
    if (step == 0) { // step 0: read sequences
		int with_qual = (!!(p->opt->flag & MM_F_OUT_SAM) && !(p->opt->flag & MM_F_NO_QUAL));
		int frag_mode = (p->n_fp > 1 || !!(p->opt->flag & MM_F_FRAG_MODE));
        step_t *s;
        s = (step_t*)calloc(1, sizeof(step_t));
		if (p->n_fp > 1) s->seq = mm_bseq_read_frag(p->n_fp, p->fp, p->mini_batch_size, with_qual, &s->n_seq);
		else s->seq = mm_bseq_read2(p->fp[0], p->mini_batch_size, with_qual, frag_mode, &s->n_seq);
		if (s->seq) {
			s->p = p;
			for (i = 0; i < s->n_seq; ++i)
				s->seq[i].rid = p->n_processed++;
			s->buf = (mm_tbuf_t**)calloc(p->n_threads, sizeof(mm_tbuf_t*));
			for (i = 0; i < p->n_threads; ++i)
				s->buf[i] = mm_tbuf_init();
			s->n_reg = (int*)calloc(3 * s->n_seq, sizeof(int));
			s->seg_off = s->n_reg + s->n_seq; // seg_off and n_seg are allocated together with n_reg
			s->n_seg = s->seg_off + s->n_seq;
			s->reg = (mm_reg1_t**)calloc(s->n_seq, sizeof(mm_reg1_t*));
			for (i = 1, j = 0; i <= s->n_seq; ++i)
				if (i == s->n_seq || !frag_mode || !mm_qname_same(s->seq[i-1].name, s->seq[i].name)) {
					s->n_seg[s->n_frag] = i - j;
					s->seg_off[s->n_frag++] = j;
					j = i;
				}
			return s;
		} else free(s);
    } else if (step == 1) { // step 1: map
		kt_for(p->n_threads, worker_for, in, ((step_t*)in)->n_frag);
		return in;
    } else if (step == 2) { // step 2: output
		void *km = 0;
        step_t *s = (step_t*)in;
		const mm_idx_t *mi = p->mi;
		for (i = 0; i < p->n_threads; ++i) mm_tbuf_destroy(s->buf[i]);
		free(s->buf);
		if ((p->opt->flag & MM_F_OUT_CS) && !(mm_dbg_flag & MM_DBG_NO_KALLOC)) km = km_init();
		for (k = 0; k < s->n_frag; ++k) {
			int seg_st = s->seg_off[k], seg_en = s->seg_off[k] + s->n_seg[k];
			for (i = seg_st; i < seg_en; ++i) {
				mm_bseq1_t *t = &s->seq[i];
				for (j = 0; j < s->n_reg[i]; ++j) {
					mm_reg1_t *r = &s->reg[i][j];
					assert(!r->sam_pri || r->id == r->parent);
					if ((p->opt->flag & MM_F_NO_PRINT_2ND) && r->id != r->parent)
						continue;
					if (p->opt->flag & MM_F_OUT_SAM)
						mm_write_sam2(&p->str, mi, t, i - seg_st, j, s->n_seg[k], &s->n_reg[seg_st], (const mm_reg1_t*const*)&s->reg[seg_st], km, p->opt->flag);
					else
						mm_write_paf(&p->str, mi, t, r, km, p->opt->flag);
					puts(p->str.s);
				}
				if (s->n_reg[i] == 0 && (p->opt->flag & MM_F_OUT_SAM)) {
					mm_write_sam2(&p->str, mi, t, i - seg_st, -1, s->n_seg[k], &s->n_reg[seg_st], (const mm_reg1_t*const*)&s->reg[seg_st], km, p->opt->flag);
					puts(p->str.s);
				}
			}
			for (i = seg_st; i < seg_en; ++i) {
				for (j = 0; j < s->n_reg[i]; ++j) free(s->reg[i][j].p);
				free(s->reg[i]);
				free(s->seq[i].seq); free(s->seq[i].name);
				if (s->seq[i].qual) free(s->seq[i].qual);
			}
		}
		free(s->reg); free(s->n_reg); free(s->seq); // seg_off and n_seg were allocated with reg; no memory leak here
		km_destroy(km);
		if (mm_verbose >= 3)
			fprintf(stderr, "[M::%s::%.3f*%.2f] mapped %d sequences\n", __func__, realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0), s->n_seq);
		free(s);
	}
    return 0;
}

int mm_map_file_frag(const mm_idx_t *idx, int n_segs, const char **fn, const mm_mapopt_t *opt, int n_threads)
{
	int i, j, pl_threads;
	pipeline_t pl;
	if (n_segs < 1) return -1;
	memset(&pl, 0, sizeof(pipeline_t));
	pl.n_fp = n_segs;
	pl.fp = (mm_bseq_file_t**)calloc(n_segs, sizeof(mm_bseq_file_t*));
	for (i = 0; i < n_segs; ++i) {
		pl.fp[i] = mm_bseq_open(fn[i]);
		if (pl.fp[i] == 0) {
			if (mm_verbose >= 1)
				fprintf(stderr, "ERROR: failed to open file '%s'\n", fn[i]);
			for (j = 0; j < i; ++j)
				mm_bseq_close(pl.fp[j]);
			free(pl.fp);
			return -1;
		}
	}
	pl.opt = opt, pl.mi = idx;
	pl.n_threads = n_threads > 1? n_threads : 1;
	pl.mini_batch_size = opt->mini_batch_size;
	pl_threads = n_threads == 1? 1 : (opt->flag&MM_F_2_IO_THREADS)? 3 : 2;
	kt_pipeline(pl_threads, worker_pipeline, &pl, 3);
	free(pl.str.s);
	for (i = 0; i < n_segs; ++i)
		mm_bseq_close(pl.fp[i]);
	free(pl.fp);
	return 0;
}

int mm_map_file(const mm_idx_t *idx, const char *fn, const mm_mapopt_t *opt, int n_threads)
{
	return mm_map_file_frag(idx, 1, &fn, opt, n_threads);
}
