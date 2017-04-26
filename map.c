#include <stdlib.h>
#include <string.h>
#include "kthread.h"
#include "kvec.h"
#include "kalloc.h"
#include "sdust.h"
#include "minimap.h"
#include "bseq.h"

void mm_mapopt_init(mm_mapopt_t *opt)
{
	opt->n_frag_mini = 100;
	opt->max_occ_frac = 1e-5f;
	opt->mid_occ_frac = 1e-3f;
	opt->sdust_thres = 0;

	opt->radius = 500;
	opt->max_gap = 10000;
	opt->min_cnt = 4;
	opt->min_match = 40;
}

void mm_mapopt_update(mm_mapopt_t *opt, const mm_idx_t *mi)
{
	opt->max_occ = mm_idx_cal_max_occ(mi, opt->max_occ_frac);
	opt->mid_occ = mm_idx_cal_max_occ(mi, opt->mid_occ_frac);
	if (mm_verbose >= 3)
		fprintf(stderr, "[M::%s::%.3f*%.2f] mid_occ = %d; max_occ = %d\n", __func__, realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0),
				opt->mid_occ, opt->max_occ);
}

typedef struct {
	uint32_t n, is_alloc;
	uint32_t qpos, span;
	union {
		const uint64_t *cr;
		uint64_t *r;
	} x;
} mm_match_t;

struct mm_tbuf_s {
	sdust_buf_t *sdb;
	mm128_v mini;
	kvec_t(mm_reg1_t) reg;
	void *km, *km_fixed;
};

mm_tbuf_t *mm_tbuf_init(void)
{
	mm_tbuf_t *b;
	b = (mm_tbuf_t*)calloc(1, sizeof(mm_tbuf_t));
	b->km = km_init();
	b->km_fixed = km_init();
	b->sdb = sdust_buf_init(b->km);
	return b;
}

void mm_tbuf_destroy(mm_tbuf_t *b)
{
	if (b == 0) return;
	km_destroy(b->km_fixed);
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

void mm_pair_thin(mm_tbuf_t *b, mm_match_t *_m1, mm_match_t *_m2)
{
}

void mm_map_frag(const mm_mapopt_t *opt, const mm_idx_t *mi, mm_tbuf_t *b, uint32_t m_st, uint32_t m_en)
{
	int i, n = m_en - m_st, last = -1, last2 = -1;
	mm_match_t *m;

	// convert to local representation
	m = (mm_match_t*)kmalloc(b->km_fixed, (m_en - m_st) * sizeof(mm_match_t));
	for (i = 0; i < n; ++i) {
		int t;
		mm128_t *p = &b->mini.a[i + m_st];
		m[i].is_alloc = 0;
		m[i].qpos = (uint32_t)p->y;
		m[i].span = p->x & 0xff;
		m[i].x.cr = mm_idx_get(mi, p->x, &t);
		m[i].n = t;
	}

	// pair k-mer thinning
	for (i = 0; i < n; ++i) {
		if (m[i].n >= opt->mid_occ && m[i].n < opt->max_occ) {
			if (last2 < 0) last2 = i;
			if (last < 0 || m[last].n < m[i].n) last = i;
			if (last >= 0 && m[last].span + (m[last].qpos>>1) <= m[i].qpos>>1) {
				mm_pair_thin(b, &m[last], &m[i]);
				last2 = last = -1;
			} else if (last2 >= 0 && m[last].span + (m[last].qpos>>1) <= m[i].qpos>>1) {
				mm_pair_thin(b, &m[last2], &m[i]);
				last2 = last = -1;
			}
		}
	}
	kfree(b->km_fixed, m);
}

const mm_reg1_t *mm_map(const mm_idx_t *mi, int l_seq, const char *seq, int *n_regs, mm_tbuf_t *b, const mm_mapopt_t *opt, const char *name)
{
	uint32_t proc_mini = 0;
	b->mini.n = 0;
	mm_sketch(b->km, seq, l_seq, mi->w, mi->k, 0, mi->is_hpc, &b->mini);
	if (opt->sdust_thres > 0)
		mm_dust_minier(&b->mini, l_seq, seq, opt->sdust_thres, b->sdb);
	while (proc_mini < b->mini.n) {
		uint32_t n = b->mini.n - proc_mini < opt->n_frag_mini * 1.5f? b->mini.n - proc_mini : opt->n_frag_mini;
		mm_map_frag(opt, mi, b, proc_mini, proc_mini + n);
		proc_mini += n;
	}
	*n_regs = b->reg.n;
	return b->reg.a;
}

/**************************
 * Multi-threaded mapping *
 **************************/

typedef struct {
	int mini_batch_size, n_processed, n_threads;
	const mm_mapopt_t *opt;
	bseq_file_t *fp;
	const mm_idx_t *mi;
} pipeline_t;

typedef struct {
	const pipeline_t *p;
    int n_seq;
	bseq1_t *seq;
	int *n_reg;
	mm_reg1_t **reg;
	mm_tbuf_t **buf;
} step_t;

static void worker_for(void *_data, long i, int tid) // kt_for() callback
{
    step_t *step = (step_t*)_data;
	const mm_reg1_t *regs;
	int n_regs;

	regs = mm_map(step->p->mi, step->seq[i].l_seq, step->seq[i].seq, &n_regs, step->buf[tid], step->p->opt, step->seq[i].name);
	step->n_reg[i] = n_regs;
	if (n_regs > 0) {
		step->reg[i] = (mm_reg1_t*)malloc(n_regs * sizeof(mm_reg1_t));
		memcpy(step->reg[i], regs, n_regs * sizeof(mm_reg1_t));
	}
}

static void *worker_pipeline(void *shared, int step, void *in)
{
	int i, j;
    pipeline_t *p = (pipeline_t*)shared;
    if (step == 0) { // step 0: read sequences
        step_t *s;
        s = (step_t*)calloc(1, sizeof(step_t));
		s->seq = bseq_read(p->fp, p->mini_batch_size, &s->n_seq);
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
        step_t *s = (step_t*)in;
		const mm_idx_t *mi = p->mi;
		for (i = 0; i < p->n_threads; ++i) mm_tbuf_destroy(s->buf[i]);
		free(s->buf);
		for (i = 0; i < s->n_seq; ++i) {
			bseq1_t *t = &s->seq[i];
			for (j = 0; j < s->n_reg[i]; ++j) {
				mm_reg1_t *r = &s->reg[i][j];
				if (r->len < p->opt->min_match) continue;
				printf("%s\t%d\t%d\t%d\t%c\t", t->name, t->l_seq, r->qs, r->qe, "+-"[r->rev]);
				if (mi->seq[r->rid].name) fputs(mi->seq[r->rid].name, stdout);
				else printf("%d", r->rid + 1);
				printf("\t%d\t%d\t%d\t%d\t%d\t255\tcm:i:%d\n", mi->seq[r->rid].len, r->rs, r->re, r->len,
						r->re - r->rs > r->qe - r->qs? r->re - r->rs : r->qe - r->qs, r->cnt);
			}
			free(s->reg[i]);
			free(s->seq[i].seq); free(s->seq[i].name);
		}
		free(s->reg); free(s->n_reg); free(s->seq);
		free(s);
	}
    return 0;
}

int mm_map_file(const mm_idx_t *idx, const char *fn, const mm_mapopt_t *opt, int n_threads, int mini_batch_size)
{
	pipeline_t pl;
	memset(&pl, 0, sizeof(pipeline_t));
	pl.fp = bseq_open(fn);
	if (pl.fp == 0) return -1;
	pl.opt = opt, pl.mi = idx;
	pl.n_threads = n_threads, pl.mini_batch_size = mini_batch_size;
	kt_pipeline(n_threads == 1? 1 : 2, worker_pipeline, &pl, 3);
	bseq_close(pl.fp);
	return 0;
}
