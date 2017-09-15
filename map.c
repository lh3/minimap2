#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "kthread.h"
#include "kvec.h"
#include "kalloc.h"
#include "sdust.h"
#include "mmpriv.h"
#include "bseq.h"

void mm_mapopt_init(mm_mapopt_t *opt)
{
	memset(opt, 0, sizeof(mm_mapopt_t));
	opt->mid_occ_frac = 2e-4f;
	opt->sdust_thres = 0;

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
	opt->min_dp_max = opt->min_chain_score * opt->a;
	opt->min_ksw_len = 200;
	opt->mini_batch_size = 200000000;
}

void mm_mapopt_update(mm_mapopt_t *opt, const mm_idx_t *mi)
{
	if (opt->flag & MM_F_SPLICE_BOTH)
		opt->flag &= ~(MM_F_SPLICE_FOR|MM_F_SPLICE_REV);
	if (opt->mid_occ <= 0)
		opt->mid_occ = mm_idx_cal_max_occ(mi, opt->mid_occ_frac);
	if (mm_verbose >= 3)
		fprintf(stderr, "[M::%s::%.3f*%.2f] mid_occ = %d\n", __func__, realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0), opt->mid_occ);
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
		mo->mini_batch_size = 500000000;
	} else if (strcmp(preset, "ava-pb") == 0) {
		io->is_hpc = 1, io->k = 19, io->w = 5;
		mo->flag |= MM_F_AVA | MM_F_NO_SELF;
		mo->min_chain_score = 100, mo->pri_ratio = 0.0f, mo->max_gap = 10000, mo->max_chain_skip = 25;
		mo->mini_batch_size = 500000000;
	} else if (strcmp(preset, "map10k") == 0 || strcmp(preset, "map-pb") == 0) {
		io->is_hpc = 1, io->k = 19;
	} else if (strcmp(preset, "map-ont") == 0) {
		io->is_hpc = 0, io->k = 15;
	} else if (strcmp(preset, "asm5") == 0) {
		io->is_hpc = 0, io->k = 19, io->w = 19;
		mo->a = 1, mo->b = 19, mo->q = 39, mo->q2 = 81, mo->e = 3, mo->e2 = 1, mo->zdrop = 200;
		mo->min_dp_max = 200;
	} else if (strcmp(preset, "asm10") == 0) {
		io->is_hpc = 0, io->k = 19, io->w = 19;
		mo->a = 1, mo->b = 9, mo->q = 16, mo->q2 = 41, mo->e = 2, mo->e2 = 1, mo->zdrop = 200;
		mo->min_dp_max = 200;
	} else if (strcmp(preset, "short") == 0 || strcmp(preset, "sr") == 0) {
		io->is_hpc = 0, io->k = 21, io->w = 11;
		mo->flag |= MM_F_APPROX_EXT;
		mo->a = 2, mo->b = 8, mo->q = 12, mo->e = 2, mo->q2 = 32, mo->e2 = 1;
		mo->max_gap = 100;
		mo->pri_ratio = 0.5f;
		mo->min_cnt = 2;
		mo->min_chain_score = 20;
		mo->min_dp_max = 40;
		mo->best_n = 20;
		mo->bw = 50;
		mo->mid_occ = 1000;
		mo->mini_batch_size = 50000000;
	} else if (strcmp(preset, "splice") == 0 || strcmp(preset, "cdna") == 0) {
		io->is_hpc = 0, io->k = 15, io->w = 5;
		mo->flag |= MM_F_SPLICE | MM_F_SPLICE_FOR | MM_F_SPLICE_REV;
		mo->max_gap = 2000, mo->max_gap_ref = mo->bw = 200000;
		mo->a = 1, mo->b = 2, mo->q = 2, mo->e = 1, mo->q2 = 32, mo->e2 = 0;
		mo->noncan = 5;
		mo->zdrop = 200;
	} else return -1;
	return 0;
}

typedef struct {
	uint32_t n;
	uint32_t qpos;
	union {
		const uint64_t *cr;
		uint64_t *r;
	} x;
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

static void mm_dust_minier(mm128_v *mini, int l_seq, const char *seq, int sdust_thres, sdust_buf_t *sdb)
{
	int n_dreg, j, k, u = 0;
	const uint64_t *dreg;
	if (sdust_thres <= 0 || sdb == 0) return;
	dreg = sdust_core((const uint8_t*)seq, l_seq, sdust_thres, 64, &n_dreg, sdb);
	for (j = k = 0; j < mini->n; ++j) { // squeeze out minimizers that significantly overlap with LCRs
		int32_t qpos = (uint32_t)mini->a[j].y>>1, span = mini->a[j].x&0xff;
		int32_t s = qpos - (span - 1), e = s + span;
		while (u < n_dreg && (uint32_t)dreg[u] <= s) ++u;
		if (u < n_dreg && dreg[u]>>32 < e) {
			int v, l = 0;
			for (v = u; v < n_dreg && dreg[v]>>32 < e; ++v) { // iterate over LCRs overlapping this minimizer
				int ss = s > dreg[v]>>32? s : dreg[v]>>32;
				int ee = e < (uint32_t)dreg[v]? e : (uint32_t)dreg[v];
				l += ee - ss;
			}
			if (l <= span>>1) mini->a[k++] = mini->a[j]; // keep the minimizer if less than half of it falls in masked region
		}
	}
	mini->n = k;
}

mm_reg1_t *mm_map(const mm_idx_t *mi, int qlen, const char *seq, int *n_regs, mm_tbuf_t *b, const mm_mapopt_t *opt, const char *qname)
{
	int i, n, j, n_u, max_gap_ref, rep_st = 0, rep_en = 0, rep_len = 0;
	int64_t n_a;
	uint64_t *u;
	mm_match_t *m;
	mm128_t *a;
	mm_reg1_t *regs;

	// collect minimizers
	b->mini.n = 0;
	mm_sketch(b->km, seq, qlen, mi->w, mi->k, 0, mi->is_hpc, &b->mini);
	n = b->mini.n;

	if (opt->sdust_thres > 0)
		mm_dust_minier(&b->mini, qlen, seq, opt->sdust_thres, b->sdb);

	// convert to local representation
	m = (mm_match_t*)kmalloc(b->km, n * sizeof(mm_match_t));
	for (i = 0; i < n; ++i) {
		int t;
		mm128_t *p = &b->mini.a[i];
		m[i].qpos = (uint32_t)p->y;
		m[i].x.cr = mm_idx_get(mi, p->x>>8, &t);
		m[i].n = t;
	}

	// fill the _a_ array
	for (i = 0, n_a = 0; i < n; ++i) // find the length of a[]
		if (m[i].n < opt->mid_occ) n_a += m[i].n;
	a = (mm128_t*)kmalloc(b->km, n_a * sizeof(mm128_t));
	for (i = j = 0; i < n; ++i) {
		mm128_t *p = &b->mini.a[i];
		mm_match_t *q = &m[i];
		const uint64_t *r = q->x.cr;
		int k, q_span = p->x & 0xff, is_tandem = 0;
		if (q->n >= opt->mid_occ) {
			int en = (q->qpos>>1) + 1, st = en - q_span;
			if (st > rep_en) {
				rep_len += rep_en - rep_st;
				rep_st = st, rep_en = en;
			} else rep_en = en;
			continue;
		}
		if (i > 0 && p->x>>8 == b->mini.a[i - 1].x>>8) is_tandem = 1;
		if (i < n - 1 && p->x>>8 == b->mini.a[i + 1].x>>8) is_tandem = 1;
		for (k = 0; k < q->n; ++k) {
			int32_t rpos = (uint32_t)r[k] >> 1;
			mm128_t *p;
			if (qname && (opt->flag&(MM_F_NO_SELF|MM_F_AVA))) {
				const char *tname = mi->seq[r[k]>>32].name;
				if ((opt->flag&MM_F_NO_SELF) && strcmp(qname, tname) == 0 && rpos == (q->qpos>>1)) // avoid the diagonal
					continue;
				if ((opt->flag&MM_F_AVA) && strcmp(qname, tname) > 0) // all-vs-all mode: map once
					continue;
			}
			p = &a[j++];
			if ((r[k]&1) == (q->qpos&1)) { // forward strand
				p->x = (r[k]&0xffffffff00000000ULL) | rpos;
				p->y = (uint64_t)q_span << 32 | q->qpos >> 1;
			} else { // reverse strand
				p->x = 1ULL<<63 | (r[k]&0xffffffff00000000ULL) | rpos;
				p->y = (uint64_t)q_span << 32 | (qlen - ((q->qpos>>1) + 1 - q_span) - 1);
			}
			if (is_tandem) p->y |= MM_SEED_TANDEM;
		}
	}
	rep_len += rep_en - rep_st;
	n_a = j;
	radix_sort_128x(a, a + n_a);
	kfree(b->km, m);

	if (mm_dbg_flag & MM_DBG_PRINT_SEED) {
		fprintf(stderr, "RS\t%d\n", rep_len);
		for (i = 0; i < n_a; ++i)
			fprintf(stderr, "SD\t%s\t%d\t%c\t%d\t%d\t%d\n", mi->seq[a[i].x<<1>>33].name, (int32_t)a[i].x, "+-"[a[i].x>>63], (int32_t)a[i].y, (int32_t)(a[i].y>>32&0xff),
					i == 0? 0 : ((int32_t)a[i].y - (int32_t)a[i-1].y) - ((int32_t)a[i].x - (int32_t)a[i-1].x));
	}

	max_gap_ref = opt->max_gap_ref >= 0? opt->max_gap_ref : opt->max_gap;
	n_u = mm_chain_dp(max_gap_ref, opt->max_gap, opt->bw, opt->max_chain_skip, opt->min_cnt, opt->min_chain_score, !!(opt->flag&MM_F_SPLICE), n_a, a, &u, b->km);
	regs = mm_gen_regs(b->km, qlen, n_u, u, a);
	*n_regs = n_u;

	if (mm_dbg_flag & MM_DBG_PRINT_SEED)
		for (j = 0; j < n_u; ++j)
			for (i = regs[j].as; i < regs[j].as + regs[j].cnt; ++i)
				fprintf(stderr, "CN\t%d\t%s\t%d\t%c\t%d\t%d\t%d\n", j, mi->seq[a[i].x<<1>>33].name, (int32_t)a[i].x, "+-"[a[i].x>>63], (int32_t)a[i].y, (int32_t)(a[i].y>>32&0xff),
						i == regs[j].as? 0 : ((int32_t)a[i].y - (int32_t)a[i-1].y) - ((int32_t)a[i].x - (int32_t)a[i-1].x));

	if (!(opt->flag & MM_F_AVA)) { // don't choose primary mapping(s) for read overlap
		mm_set_parent(b->km, opt->mask_level, *n_regs, regs, opt->a * 2 + opt->b);
		mm_select_sub(b->km, opt->mask_level, opt->pri_ratio, mi->k*2, opt->best_n, n_regs, regs);
		if (!(opt->flag & MM_F_SPLICE))
			mm_join_long(b->km, opt, qlen, n_regs, regs, a); // TODO: this can be applied to all-vs-all in principle
	}
	if (opt->flag & MM_F_CIGAR) {
		regs = mm_align_skeleton(b->km, opt, mi, qlen, seq, n_regs, regs, a); // this calls mm_filter_regs()
		if (!(opt->flag & MM_F_AVA)) {
			mm_set_parent(b->km, opt->mask_level, *n_regs, regs, opt->a * 2 + opt->b);
			mm_select_sub(b->km, opt->mask_level, opt->pri_ratio, mi->k*2, opt->best_n, n_regs, regs);
			mm_set_sam_pri(*n_regs, regs);
		}
	}
	mm_set_mapq(*n_regs, regs, opt->min_chain_score, opt->a, rep_len);

	// free
	kfree(b->km, a);
	kfree(b->km, u);
	return regs;
}

/**************************
 * Multi-threaded mapping *
 **************************/

typedef struct {
	int mini_batch_size, n_processed, n_threads;
	const mm_mapopt_t *opt;
	mm_bseq_file_t *fp;
	const mm_idx_t *mi;
	kstring_t str;
} pipeline_t;

typedef struct {
	const pipeline_t *p;
    int n_seq;
	mm_bseq1_t *seq;
	int *n_reg;
	mm_reg1_t **reg;
	mm_tbuf_t **buf;
} step_t;

static void worker_for(void *_data, long i, int tid) // kt_for() callback
{
    step_t *step = (step_t*)_data;
	if (mm_dbg_flag & MM_DBG_PRINT_QNAME)
		fprintf(stderr, "QR\t%s\t%d\n", step->seq[i].name, tid);
	step->reg[i] = mm_map(step->p->mi, step->seq[i].l_seq, step->seq[i].seq, &step->n_reg[i], step->buf[tid], step->p->opt, step->seq[i].name);
}

static void *worker_pipeline(void *shared, int step, void *in)
{
	int i, j;
    pipeline_t *p = (pipeline_t*)shared;
    if (step == 0) { // step 0: read sequences
		int with_qual = (!!(p->opt->flag & MM_F_OUT_SAM) && !(p->opt->flag & MM_F_NO_QUAL));
        step_t *s;
        s = (step_t*)calloc(1, sizeof(step_t));
		s->seq = mm_bseq_read(p->fp, p->mini_batch_size, with_qual, &s->n_seq);
		if (s->seq) {
			s->p = p;
			for (i = 0; i < s->n_seq; ++i)
				s->seq[i].rid = p->n_processed++;
			s->buf = (mm_tbuf_t**)calloc(p->n_threads, sizeof(mm_tbuf_t*));
			for (i = 0; i < p->n_threads; ++i)
				s->buf[i] = mm_tbuf_init();
			s->n_reg = (int*)calloc(s->n_seq, sizeof(int));
			s->reg = (mm_reg1_t**)calloc(s->n_seq, sizeof(mm_reg1_t*));
			return s;
		} else free(s);
    } else if (step == 1) { // step 1: map
		kt_for(p->n_threads, worker_for, in, ((step_t*)in)->n_seq);
		return in;
    } else if (step == 2) { // step 2: output
		void *km = 0;
        step_t *s = (step_t*)in;
		const mm_idx_t *mi = p->mi;
		for (i = 0; i < p->n_threads; ++i) mm_tbuf_destroy(s->buf[i]);
		free(s->buf);
		if ((p->opt->flag & MM_F_OUT_CS) && !(mm_dbg_flag & MM_DBG_NO_KALLOC)) km = km_init();
		for (i = 0; i < s->n_seq; ++i) {
			mm_bseq1_t *t = &s->seq[i];
			for (j = 0; j < s->n_reg[i]; ++j) {
				mm_reg1_t *r = &s->reg[i][j];
				if (p->opt->flag & MM_F_OUT_SAM)
					mm_write_sam(&p->str, mi, t, r, s->n_reg[i], s->reg[i]);
				else
					mm_write_paf(&p->str, mi, t, r, km, p->opt->flag);
				puts(p->str.s);
			}
			if (s->n_reg[i] == 0 && (p->opt->flag & MM_F_OUT_SAM)) {
				mm_write_sam(&p->str, 0, t, 0, 0, 0);
				puts(p->str.s);
			}
			for (j = 0; j < s->n_reg[i]; ++j) free(s->reg[i][j].p);
			free(s->reg[i]);
			free(s->seq[i].seq); free(s->seq[i].name);
			if (s->seq[i].qual) free(s->seq[i].qual);
		}
		free(s->reg); free(s->n_reg); free(s->seq);
		km_destroy(km);
		if (mm_verbose >= 3)
			fprintf(stderr, "[M::%s::%.3f*%.2f] mapped %d sequences\n", __func__, realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0), s->n_seq);
		free(s);
	}
    return 0;
}

int mm_map_file(const mm_idx_t *idx, const char *fn, const mm_mapopt_t *opt, int n_threads)
{
	pipeline_t pl;
	memset(&pl, 0, sizeof(pipeline_t));
	pl.fp = mm_bseq_open(fn);
	if (pl.fp == 0) {
		if (mm_verbose >= 1)
			fprintf(stderr, "ERROR: failed to open file '%s'\n", fn);
		return -1;
	}
	pl.opt = opt, pl.mi = idx;
	pl.n_threads = n_threads, pl.mini_batch_size = opt->mini_batch_size;
	if ((opt->flag & MM_F_OUT_SAM) && !(opt->flag & MM_F_NO_SAM_SQ))
		mm_write_sam_SQ(idx);
	kt_pipeline(n_threads == 1? 1 : 2, worker_pipeline, &pl, 3);
	free(pl.str.s);
	mm_bseq_close(pl.fp);
	return 0;
}
