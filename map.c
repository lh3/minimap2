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
	opt->mid_occ_frac = 2e-4f;
	opt->sdust_thres = 0;
	opt->min_cnt = 2;
	opt->radius = 500;

	opt->max_gap = 10000;
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

int mm_pair_thin_core(mm_tbuf_t *b, uint64_t x, int radius, int rel, int st0, int n, const uint64_t *z, uint64_v *a)
{
	int i, st = st0, en = n, mid = en - 1;
	while (st < en) {
		uint64_t y;
		mid = st + ((en - st) >> 1);
		y = z[mid];
		if (y < x && (x - y)>>1 > radius) st = mid + 1;
		else if (y >= x && (y - x)>>1 > radius) en = mid;
		else break;
	}
	if (st < en) {
		for (en = mid + 1; en < n; ++en)
			if (z[en] > x && (z[en] - x)>>1 > radius)
				break;
		for (st = mid - 1; st >= st0; --st)
			if (z[st] < x && (x - z[st])>>1 > radius)
				break;
		++st;
		for (i = st; i < en; ++i) {
			uint64_t y = z[i];
			if (((x ^ y) & 1) == rel) {
//				printf("* %d,%d\n", (uint32_t)x>>1, (uint32_t)y>>1);
				kv_push(uint64_t, b->km_fixed, *a, y);
			}
		}
		return en;
	} else return z[st] < x? st + 1 : en;
}

void mm_pair_thin(mm_tbuf_t *b, int radius, mm_match_t *m1, mm_match_t *m2)
{
	mm_match_t *m[2];
	const uint64_t *z[2];
	uint64_v a[2];
	int i, n[2], k[2], u = 0, rel = (m1->qpos ^ m2->qpos) & 1;

	m[0] = m1, m[1] = m2;
	for (i = 0; i < 2; ++i) {
		n[i] = m[i]->n;
		z[i] = m[i]->x.cr;
		k[i] = 0;
		kv_init(a[i]);
		kv_resize(uint64_t, b->km_fixed, a[i], 256);
	}
	while (k[0] < n[0] && k[1] < n[1]) {
		//printf("%d; %d,%d\n", u, k[0], k[1]);
		int v = u^1, dist = (int)(m[v]->qpos>>1) - (int)(m[u]->qpos>>1);
		uint64_t x = z[u][k[u]];
		int uori = (x ^ m[u]->qpos) & 1, last;
		int64_t tpos = x>>1 & 0x7fffffff;
		tpos = uori == 0? tpos + dist : tpos - dist;
		if (tpos < 0) tpos = 0;
		x = x>>32<<32 | tpos<<1 | (x&1);
		last = a[v].n;
		k[v] = mm_pair_thin_core(b, x, radius, rel, k[v], n[v], z[v], &a[v]);
		if (a[v].n > last) kv_push(uint64_t, b->km_fixed, a[u], z[u][k[u]]);
		++k[u];
		u ^= 1;
	}
	for (i = 0; i < 2; ++i)
		m[i]->n = a[i].n, m[i]->x.r = a[i].a, m[i]->is_alloc = 1;
//	printf("%d,%d; %d,%d\n", m[0]->qpos>>1, m[1]->qpos>>1, m[0]->n, m[1]->n);
}

int mm_liss(int d, int min_L, void *km, int n, mm128_t *a, int *b)
{
	int i, o = 0, L = 0, *M, *P, m = 0, off = 0;
	M = (int*)kmalloc(km, (n + 1) * sizeof(int));
	P = (int*)kmalloc(km, n * sizeof(int));
	for (i = 0; i <= n; ++i) {
		int lo, hi, newL;
		if (L > 0 || i == n) {
			int reset = 0, j = M[L]; // _j_ is the index of the last element in the longest chain
			if (i == n) reset = 1; // reaching the end
			else if (a[i].x - a[j].x > d) reset = 1; // large gap on target
			else if ((uint32_t)a[i].y > (uint32_t)a[j].y && (uint32_t)a[i].y - (uint32_t)a[j].y > d) reset = 1; // large gap on query
			if (reset) {
				if (L >= min_L) { // save the chain if long enough
					int k;
					for (k = L - 1; k >= 0; --k)
						a[o + k] = a[j], j = P[j];
					o += L;
					b[m++] = L;
				}
				off = i, L = 0; // reset
			}
			if (i == n) break;
		}
		lo = off + 1, hi = off + L; // a binary search; see wiki
		while (lo <= hi) {
			int mid = (lo + hi + 1) >> 1;
			if (a[M[mid - off]].x < a[i].x) lo = mid + 1;
			else hi = mid - 1;
		}
		newL = lo - off, P[i] = M[newL - 1], M[newL] = i;
		if (newL > L) L = newL;
	}
	kfree(km, M);
	kfree(km, P);
	return m;
}

void mm_map_frag(const mm_mapopt_t *opt, const mm_idx_t *mi, mm_tbuf_t *b, uint32_t m_st, uint32_t m_en, const char *qname, int qlen)
{
	int i, n = m_en - m_st, last = -1, last2 = -1, j, n_a, *t, n_t;
	mm_match_t *m;
	mm128_t *a;

	// convert to local representation
	m = (mm_match_t*)kmalloc(b->km_fixed, n * sizeof(mm_match_t));
	if (mm_verbose >= 5) printf("%d\t", n);
	for (i = 0; i < n; ++i) {
		int t;
		mm128_t *p = &b->mini.a[i + m_st];
		m[i].is_alloc = 0;
		m[i].qpos = (uint32_t)p->y;
		m[i].span = p->x & 0xff;
		m[i].x.cr = mm_idx_get(mi, p->x>>8, &t);
		m[i].n = t;
		if (mm_verbose >= 5) {
			if (i) printf("; ");
			printf("%d,%d", m[i].qpos>>1, t);
		}
	}
	if (mm_verbose >= 5) printf("\n");

	// pair k-mer thinning
	for (i = 0; i < n; ++i) {
		if (m[i].n >= opt->mid_occ) {
			if (last2 < 0) last2 = i;
			if (last < 0 || m[last].n < m[i].n) last = i;
			if (last >= 0 && (m[last].qpos>>1) + (m[last].span>>0) <= m[i].qpos>>1) {
				mm_pair_thin(b, opt->radius, &m[last], &m[i]);
				last2 = last = -1;
			} else if (last2 >= 0 && (m[last2].qpos>>1) + (m[last2].span>>0) <= m[i].qpos>>1) {
				mm_pair_thin(b, opt->radius, &m[last2], &m[i]);
				last2 = last = -1;
			}
		}
	}

	// find the length of _a_
	for (i = 0, n_a = 0; i < n; ++i)
		if (m[i].n < opt->mid_occ) n_a += m[i].n;

	// fill the _a_ array
	a = kmalloc(b->km_fixed, n_a * sizeof(mm128_t));
	for (i = j = 0; i < n; ++i) {
		mm_match_t *q = &m[i];
		const uint64_t *r = q->x.cr;
		int k;
		if (q->n >= opt->mid_occ) continue;
		for (k = 0; k < q->n; ++k) {
			int32_t rpos = (uint32_t)r[k] >> 1;
			mm128_t *p;
			if (qname && (opt->flag&MM_F_NO_SELF) && strcmp(qname, mi->seq[r[k]>>32].name) == 0 && rpos == (q->qpos>>1)) // avoid the diagonal
				continue;
			if (qname && (opt->flag&MM_F_AVA) && strcmp(qname, mi->seq[r[k]>>32].name) > 0) // all-vs-all mode: map once
				continue;
			p = &a[j++];
			if ((r[k]&1) == (q->qpos&1)) { // forward strand
				p->x = r[k];
				p->y = (uint64_t)i << 32 | q->qpos >> 1;
			} else { // reverse strand
				p->x = 1ULL<<63 | r[k];
				p->y = (uint64_t)i << 32 | (qlen - ((q->qpos>>1) + 1 - q->span) - 1);
			}
		}
	}
	radix_sort_128x(a, a + n_a);
//	for (i = 0; i < n_a; ++i) printf("%c\t%d\t%d\n", "+-"[a[i].x>>63], (uint32_t)a[i].x>>1, (uint32_t)a[i].y);

	t = kmalloc(b->km_fixed, (n_a + opt->min_cnt - 1) / opt->min_cnt * sizeof(int));
	n_t = mm_liss(opt->radius, opt->min_cnt, b->km_fixed, n_a, a, t);

	int n_tot = 0, n_low = 0, s_low = 0;
	for (i = 0; i < n; ++i) {
		if (m[i].n > 0) ++n_tot;
		if (m[i].n > 0 && m[i].n < opt->mid_occ) ++n_low, s_low += m[i].n;
	}
	printf("%d\t%d\t%d\t%d\t%d", n_tot, n_low, s_low, n, n_t);
	for (i = 0; i < n_t; ++i) printf("\t%d", t[i]);
	printf("\n");

	// free
	for (i = 0; i < n; ++i)
		if (m[i].is_alloc) kfree(b->km_fixed, m[i].x.r);
	kfree(b->km_fixed, m);
	kfree(b->km_fixed, a);
	kfree(b->km_fixed, t);
}

const mm_reg1_t *mm_map(const mm_idx_t *mi, int l_seq, const char *seq, int *n_regs, mm_tbuf_t *b, const mm_mapopt_t *opt, const char *qname)
{
	uint32_t proc_mini = 0;
	if (mm_verbose >= 5) printf("=====> %s <=====\n", qname);
	b->mini.n = 0;
	mm_sketch(b->km, seq, l_seq, mi->w, mi->k, 0, mi->is_hpc, &b->mini);
	if (opt->sdust_thres > 0)
		mm_dust_minier(&b->mini, l_seq, seq, opt->sdust_thres, b->sdb);
	while (proc_mini < b->mini.n) {
		uint32_t n = b->mini.n - proc_mini < opt->n_frag_mini * 1.5f? b->mini.n - proc_mini : opt->n_frag_mini;
		mm_map_frag(opt, mi, b, proc_mini, proc_mini + n, qname, l_seq);
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
