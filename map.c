#include <stdlib.h>
#include "kvec.h"
#include "kalloc.h"
#include "sdust.h"
#include "minimap.h"

void mm_mapopt_init(mm_mapopt_t *opt)
{
	opt->radius = 500;
	opt->max_gap = 10000;
	opt->min_cnt = 4;
	opt->min_match = 40;
	opt->sdust_thres = 0;
	opt->flag = MM_F_WITH_REP;
	opt->merge_frac = .5;
}

struct mm_tbuf_s {
	sdust_buf_t *sdb;
	mm128_v mini;
	kvec_t(mm_reg1_t) reg;
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
	km_destroy(b->km);
	free(b);
}

const mm_reg1_t *mm_map(const mm_idx_t *mi, int l_seq, const char *seq, int *n_regs, mm_tbuf_t *b, const mm_mapopt_t *opt, const char *name)
{
	int j;

	b->mini.n = 0;
	mm_sketch(b->km, seq, l_seq, mi->w, mi->k, 0, mi->is_hpc, &b->mini);
	if (opt->sdust_thres > 0) { // squeeze out minimizers that significantly overlap with LCRs
		int n_dreg, k, u = 0;
		const uint64_t *dreg;
		dreg = sdust_core((const uint8_t*)seq, l_seq, opt->sdust_thres, 64, &n_dreg, b->sdb);
		for (j = k = 0; j < b->mini.n; ++j) {
			int32_t qpos = (uint32_t)b->mini.a[j].y>>1, span = b->mini.a[j].x&0xff;
			int32_t s = qpos - (span - 1), e = s + span;
			while (u < n_dreg && (uint32_t)dreg[u] <= s) ++u;
			if (u < n_dreg && dreg[u]>>32 < e) {
				int v, l = 0;
				for (v = u; v < n_dreg && dreg[v]>>32 < e; ++v) { // iterate over LCRs overlapping this minimizer
					int ss = s > dreg[v]>>32? s : dreg[v]>>32;
					int ee = e < (uint32_t)dreg[v]? e : (uint32_t)dreg[v];
					l += ee - ss;
				}
				if (l <= mi->k>>1) b->mini.a[k++] = b->mini.a[j];
			}
		}
		b->mini.n = k;
	}
	*n_regs = b->reg.n;
	return b->reg.a;
}
