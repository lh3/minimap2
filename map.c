#include <stdlib.h>
#include <string.h>
#include "kthread.h"
#include "kvec.h"
#include "kalloc.h"
#include "sdust.h"
#include "mmpriv.h"
#include "bseq.h"

void mm_mapopt_init(mm_mapopt_t *opt)
{
	memset(opt, 0, sizeof(mm_mapopt_t));
	opt->max_occ_frac = 1e-5f;
	opt->mid_occ_frac = 2e-4f;
	opt->sdust_thres = 0;

	opt->min_score = 40;
	opt->bw = 1000;
	opt->max_gap = 10000;
	opt->max_skip = 15;

	opt->pri_ratio = 2.0f;
	opt->mask_level = 0.5f;

	opt->a = 1, opt->b = 1, opt->q = 1, opt->e = 1;
	opt->zdrop = 100;
	opt->min_ksw_len = 100;
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
	uint32_t n:31, is_alloc:1;
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
	b->km = km_init();
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
#if 0
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
				kv_push(uint64_t, b->km, *a, y);
			}
		}
		return en;
	} else return st < n && z[st] < x? st + 1 : en;
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
		kv_resize(uint64_t, b->km, a[i], 256);
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
		if (a[v].n > last) kv_push(uint64_t, b->km, a[u], z[u][k[u]]);
		++k[u];
		u ^= 1;
	}
	for (i = 0; i < 2; ++i)
		m[i]->n = a[i].n, m[i]->x.r = a[i].a, m[i]->is_alloc = 1;
//	printf("%d,%d; %d,%d\n", m[0]->qpos>>1, m[1]->qpos>>1, m[0]->n, m[1]->n);
}
#endif
mm_reg1_t *mm_gen_reg(int qlen, int n_u, uint64_t *u, mm128_t *a)
{
	mm_reg1_t *r;
	int i, k;
	r = (mm_reg1_t*)calloc(n_u, sizeof(mm_reg1_t));
	for (i = k = 0; i < n_u; ++i) {
		mm_reg1_t *ri = &r[i];
		ri->parent = -1, ri->subsc = 0;
		ri->score = u[i]>>32;
		ri->cnt = (int32_t)u[i];
		ri->as = k;
		mm_reg_set_coor(ri, qlen, a);
		k += ri->cnt;
	}
	return r;
}

static void mm_set_subsc(float mask_level, int n, mm_reg1_t *r, void *km)
{
	int i, j, k, *w;
	if (n <= 0) return;
	w = (int*)kmalloc(km, n * sizeof(int));
	w[0] = 0, r[0].parent = 0;
	for (i = 1, k = 1; i < n; ++i) {
		int si = r[i].qs, ei = r[i].qe;
		for (j = 0; j < k; ++j) {
			int sj = r[w[j]].qs, ej = r[w[j]].qe;
			int min = ej - sj < ei - si? ej - sj : ei - si;
			int ol = si < sj? (ei < sj? 0 : ei < ej? ei - sj : ej - sj) : (ej < si? 0 : ej < ei? ej - si : ei - si);
			if (ol > mask_level * min) {
				r[i].parent = r[w[j]].parent;
				if (r[w[j]].subsc < r[i].score)
					r[w[j]].subsc = r[i].score;
				break;
			}
		}
		if (j == k) w[k++] = i, r[i].parent = i;
	}
	kfree(km, w);
}

void mm_select_sub(float mask_level, float pri_ratio, int *n_, mm_reg1_t *r, void *km)
{
	if (pri_ratio > 0.0f && *n_ > 0) {
		int i, k, *tmp, n = *n_;
		mm_set_subsc(mask_level, n, r, km);
		tmp = (int*)kmalloc(km, n * sizeof(int));
		for (i = 0; i < n; ++i) tmp[i] = -1;
		for (i = k = 0; i < n; ++i)
			if (r[i].parent == i || r[i].score >= r[r[i].parent].score * pri_ratio)
				tmp[i] = k, r[k++] = r[i];
		n = k;
		for (i = 0; i < n; ++i)
			if (tmp[r[i].parent] >= 0)
				r[i].parent = tmp[r[i].parent];
		kfree(km, tmp);
		*n_ = n;
	}
}

mm_reg1_t *mm_map_frag(const mm_mapopt_t *opt, const mm_idx_t *mi, mm_tbuf_t *b, uint32_t m_st, uint32_t m_en, const char *qname, int qlen, const char *seq, int *n_regs)
{
	int i, n = m_en - m_st, j, n_a, n_u;
	uint64_t *u;
	mm_match_t *m;
	mm128_t *a;
	mm_reg1_t *regs;

	// convert to local representation
	m = (mm_match_t*)kmalloc(b->km, n * sizeof(mm_match_t));
	for (i = 0; i < n; ++i) {
		int t;
		mm128_t *p = &b->mini.a[i + m_st];
		m[i].is_alloc = 0;
		m[i].qpos = (uint32_t)p->y;
		m[i].x.cr = mm_idx_get(mi, p->x>>8, &t);
		m[i].n = t;
	}
#if 0
	int last = -1, last2 = -1;
	// pair k-mer thinning
	for (i = 0; i < n; ++i) {
		if (m[i].n >= opt->mid_occ && m[i].n < opt->max_occ) {
			if (last2 < 0) last2 = i;
			if (last < 0 || m[last].n < m[i].n) last = i;
			if (last >= 0 && (m[last].qpos>>1) + (m[last].span>>1) <= m[i].qpos>>1) {
				mm_pair_thin(b, opt->bw, &m[last], &m[i]);
				last2 = last = -1;
			} else if (last2 >= 0 && (m[last2].qpos>>1) + (m[last2].span>>1) <= m[i].qpos>>1) {
				mm_pair_thin(b, opt->bw, &m[last2], &m[i]);
				last2 = last = -1;
			}
		}
	}
#endif
	// fill the _a_ array
	for (i = 0, n_a = 0; i < n; ++i) // find the length of a[]
		if (m[i].n < opt->mid_occ) n_a += m[i].n;
	a = (mm128_t*)kmalloc(b->km, n_a * sizeof(mm128_t));
	for (i = j = 0; i < n; ++i) {
		mm128_t *p = &b->mini.a[i + m_st];
		mm_match_t *q = &m[i];
		const uint64_t *r = q->x.cr;
		int k, q_span = p->x & 0xff;
		if (q->n >= opt->mid_occ) continue;
		for (k = 0; k < q->n; ++k) {
			const char *tname = mi->seq[r[k]>>32].name;
			int32_t rpos = (uint32_t)r[k] >> 1;
			mm128_t *p;
			if (qname && (opt->flag&MM_F_NO_SELF) && strcmp(qname, tname) == 0 && rpos == (q->qpos>>1)) // avoid the diagonal
				continue;
			if (qname && (opt->flag&MM_F_AVA) && strcmp(qname, tname) > 0) // all-vs-all mode: map once
				continue;
			p = &a[j++];
			if ((r[k]&1) == (q->qpos&1)) { // forward strand
				p->x = (r[k]&0xffffffff00000000ULL) | (uint32_t)r[k]>>1;
				p->y = (uint64_t)q_span << 32 | q->qpos >> 1;
			} else { // reverse strand
				p->x = 1ULL<<63 | (r[k]&0xffffffff00000000ULL) | (uint32_t)r[k]>>1;
				p->y = (uint64_t)q_span << 32 | (qlen - ((q->qpos>>1) + 1 - q_span) - 1);
			}
		}
	}
	radix_sort_128x(a, a + n_a);
	for (i = 0; i < n; ++i)
		if (m[i].is_alloc) kfree(b->km, m[i].x.r);
	kfree(b->km, m);
	//for (i = 0; i < n_a; ++i) printf("%c\t%s\t%d\t%d\n", "+-"[a[i].x>>63], mi->seq[a[i].x<<1>>33].name, (uint32_t)a[i].x, (uint32_t)a[i].y);

	n_u = mm_chain_dp(opt->max_gap, opt->bw, opt->max_skip, opt->min_score, n_a, a, &u, b->km);
	regs = mm_gen_reg(qlen, n_u, u, a);
	*n_regs = n_u;
	mm_select_sub(opt->mask_level, opt->pri_ratio, n_regs, regs, b->km);
	if (opt->flag & MM_F_CIGAR)
		regs = mm_align_skeleton(b->km, opt, mi, qlen, seq, n_regs, regs, a);

	// free
	kfree(b->km, a);
	kfree(b->km, u);
	return regs;
}

mm_reg1_t *mm_map(const mm_idx_t *mi, int l_seq, const char *seq, int *n_regs, mm_tbuf_t *b, const mm_mapopt_t *opt, const char *qname)
{
	mm_reg1_t *regs;
	if (mm_verbose >= 10) fprintf(stderr, "===> %s <===\n", qname);
	b->mini.n = 0;
	mm_sketch(b->km, seq, l_seq, mi->w, mi->k, 0, mi->is_hpc, &b->mini);
	if (opt->sdust_thres > 0)
		mm_dust_minier(&b->mini, l_seq, seq, opt->sdust_thres, b->sdb);
	regs = mm_map_frag(opt, mi, b, 0, b->mini.n, qname, l_seq, seq, n_regs);
	return regs;
}

/**************************
 * Multi-threaded mapping *
 **************************/

typedef struct {
	int mini_batch_size, n_processed, n_threads;
	const mm_mapopt_t *opt;
	bseq_file_t *fp;
	const mm_idx_t *mi;
	kstring_t str;
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
	step->reg[i] = mm_map(step->p->mi, step->seq[i].l_seq, step->seq[i].seq, &step->n_reg[i], step->buf[tid], step->p->opt, step->seq[i].name);
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
				if (r->p && r->p->blen - r->p->n_ambi - r->p->n_diff < p->opt->min_score)
					continue;
				if (p->opt->flag & MM_F_OUT_SAM) mm_write_sam(&p->str, mi, t, j, r);
				else mm_write_paf(&p->str, mi, t, j, r);
				puts(p->str.s);
				free(r->p);
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
	if (opt->flag & MM_F_OUT_SAM) {
		uint32_t i;
		for (i = 0; i < idx->n_seq; ++i)
			printf("@SQ\tID:%s\tLN:%d\n", idx->seq[i].name, idx->seq[i].len);
	}
	kt_pipeline(n_threads == 1? 1 : 2, worker_pipeline, &pl, 3);
	free(pl.str.s);
	bseq_close(pl.fp);
	return 0;
}
