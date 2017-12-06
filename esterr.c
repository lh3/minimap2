#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "mmpriv.h"

static inline int32_t get_for_qpos(int32_t qlen, const mm128_t *a)
{
	int32_t x = (int32_t)a->y;
	int32_t q_span = a->y>>32 & 0xff;
	if (a->x>>63)
		x = qlen - 1 - (x + 1 - q_span); // revert the position to the forward strand of query
	return x;
}

static int get_mini_idx(int qlen, const mm128_t *a, int32_t n, const uint64_t *mini_pos)
{
	int32_t x, L = 0, R = n - 1;
	x = get_for_qpos(qlen, a);
	while (L <= R) { // binary search
		int32_t m = ((uint64_t)L + R) >> 1;
		int32_t y = (int32_t)mini_pos[m];
		if (y < x) L = m + 1;
		else if (y > x) R = m - 1;
		else return m;
	}
	return -1;
}

void mm_est_err(int qlen, int n_regs, mm_reg1_t *regs, const mm128_t *a, int32_t n, uint64_t *mini_pos)
{
	int i;
	uint64_t sum_k = 0;
	float avg_k;

	if (n == 0) return;
	for (i = 0; i < n; ++i)
		sum_k += mini_pos[i] >> 32 & 0xff;
	avg_k = (float)sum_k / n;

	for (i = 0; i < n_regs; ++i) {
		mm_reg1_t *r = &regs[i];
		int32_t st, en, j, k, n_mis = 0;
		r->div = -1.0f;
		if (r->cnt == 0) continue;
		st = en = get_mini_idx(qlen, r->rev? &a[r->as + r->cnt - 1] : &a[r->as], n, mini_pos);
		assert(st >= 0);
		for (k = 1, j = st + 1; j < n && k < r->cnt; ++j) {
			int32_t x;
			x = get_for_qpos(qlen, r->rev? &a[r->as + r->cnt - 1 - k] : &a[r->as + k]);
			if (x == (int32_t)mini_pos[j])
				mini_pos[j] |= 1ULL<<40, ++k, en = j;
			else mini_pos[j] &= ~(1ULL<<40), ++n_mis;
		}
		r->div = -logf((float)(en - st + 1 - n_mis) / (en - st + 1)) / avg_k;
	}
}
