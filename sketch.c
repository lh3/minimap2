#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#define __STDC_LIMIT_MACROS
#include "kvec.h"
#include "mmpriv.h"

unsigned char seq_nt4_table[256] = {
	0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

static inline uint64_t hash64(uint64_t key, uint64_t mask)
{
	key = (~key + (key << 21)) & mask; // key = (key << 21) - key - 1;
	key = key ^ key >> 24;
	key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
	key = key ^ key >> 14;
	key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
	key = key ^ key >> 28;
	key = (key + (key << 31)) & mask;
	return key;
}

typedef struct { // a simplified version of kdq
	int front, count;
	int a[32];
} tiny_queue_t;

static inline void tq_push(tiny_queue_t *q, int x)
{
	q->a[((q->count++) + q->front) & 0x1f] = x;
}

static inline int tq_shift(tiny_queue_t *q)
{
	int x;
	if (q->count == 0) return -1;
	x = q->a[q->front++];
	q->front &= 0x1f;
	--q->count;
	return x;
}

/**
 * Find symmetric (w,k)-minimizers on a DNA sequence
 *
 * @param km     thread-local memory pool; using NULL falls back to malloc()
 * @param str    DNA sequence
 * @param len    length of $str
 * @param w      find a minimizer for every $w consecutive k-mers
 * @param k      k-mer size
 * @param rid    reference ID; will be copied to the output $p array
 * @param is_hpc homopolymer-compressed or not
 * @param p      minimizers
 *               p->a[i].x = kMer<<8 | kmerSpan
 *               p->a[i].y = rid<<32 | lastPos<<1 | strand
 *               where lastPos is the position of the last base of the i-th minimizer,
 *               and strand indicates whether the minimizer comes from the top or the bottom strand.
 *               Callers may want to set "p->n = 0"; otherwise results are appended to p
 */
static void mm_sketch_random(void *km, const char *str, int len, int w, int k, uint32_t rid, int is_hpc, mm128_v *p)
{
	uint64_t shift1 = 2 * (k - 1), mask = (1ULL<<2*k) - 1, kmer[2] = {0,0};
	int i, j, l, buf_pos, min_pos, kmer_span = 0;
	mm128_t buf[256], min = { UINT64_MAX, UINT64_MAX };
	tiny_queue_t tq;

	assert(len > 0 && (w > 0 && w < 256) && (k > 0 && k <= 28)); // 56 bits for k-mer; could use long k-mers, but 28 enough in practice
	memset(buf, 0xff, w * 16);
	memset(&tq, 0, sizeof(tiny_queue_t));
	kv_resize(mm128_t, km, *p, p->n + len/w);

	for (i = l = buf_pos = min_pos = 0; i < len; ++i) {
		int c = seq_nt4_table[(uint8_t)str[i]];
		mm128_t info = { UINT64_MAX, UINT64_MAX };
		if (c < 4) { // not an ambiguous base
			int z;
			if (is_hpc) {
				int skip_len = 1;
				if (i + 1 < len && seq_nt4_table[(uint8_t)str[i + 1]] == c) {
					for (skip_len = 2; i + skip_len < len; ++skip_len)
						if (seq_nt4_table[(uint8_t)str[i + skip_len]] != c)
							break;
					i += skip_len - 1; // put $i at the end of the current homopolymer run
				}
				tq_push(&tq, skip_len);
				kmer_span += skip_len;
				if (tq.count > k) kmer_span -= tq_shift(&tq);
			} else kmer_span = l + 1 < k? l + 1 : k;
			kmer[0] = (kmer[0] << 2 | c) & mask;           // forward k-mer
			kmer[1] = (kmer[1] >> 2) | (3ULL^c) << shift1; // reverse k-mer
			if (kmer[0] == kmer[1]) continue; // skip "symmetric k-mers" as we don't know it strand
			z = kmer[0] < kmer[1]? 0 : 1; // strand
			++l;
			if (l >= k && kmer_span < 256) {
				info.x = hash64(kmer[z], mask) << 8 | kmer_span;
				info.y = (uint64_t)rid<<32 | (uint32_t)i<<1 | z;
			}
		} else l = 0, tq.count = tq.front = 0, kmer_span = 0;
		buf[buf_pos] = info; // need to do this here as appropriate buf_pos and buf[buf_pos] are needed below
		if (l == w + k - 1 && min.x != UINT64_MAX) { // special case for the first window - because identical k-mers are not stored yet
			for (j = buf_pos + 1; j < w; ++j)
				if (min.x == buf[j].x && buf[j].y != min.y) kv_push(mm128_t, km, *p, buf[j]);
			for (j = 0; j < buf_pos; ++j)
				if (min.x == buf[j].x && buf[j].y != min.y) kv_push(mm128_t, km, *p, buf[j]);
		}
		if (info.x <= min.x) { // a new minimum; then write the old min
			if (l >= w + k && min.x != UINT64_MAX) kv_push(mm128_t, km, *p, min);
			min = info, min_pos = buf_pos;
		} else if (buf_pos == min_pos) { // old min has moved outside the window
			if (l >= w + k - 1 && min.x != UINT64_MAX) kv_push(mm128_t, km, *p, min);
			for (j = buf_pos + 1, min.x = UINT64_MAX; j < w; ++j) // the two loops are necessary when there are identical k-mers
				if (min.x >= buf[j].x) min = buf[j], min_pos = j; // >= is important s.t. min is always the closest k-mer
			for (j = 0; j <= buf_pos; ++j)
				if (min.x >= buf[j].x) min = buf[j], min_pos = j;
			if (l >= w + k - 1 && min.x != UINT64_MAX) { // write identical k-mers
				for (j = buf_pos + 1; j < w; ++j) // these two loops make sure the output is sorted
					if (min.x == buf[j].x && min.y != buf[j].y) kv_push(mm128_t, km, *p, buf[j]);
				for (j = 0; j <= buf_pos; ++j)
					if (min.x == buf[j].x && min.y != buf[j].y) kv_push(mm128_t, km, *p, buf[j]);
			}
		}
		if (++buf_pos == w) buf_pos = 0;
	}
	if (min.x != UINT64_MAX)
		kv_push(mm128_t, km, *p, min);
}

/*
 * Mod-minimizer sampling (Groot Koerkamp & Pibiri, WABI 2024).
 * The minimum t-mer in the w+k-t t-mers spanning a window chooses the k-mer
 * at its position modulo w. HPC and t==k are routed to the stock sketcher.
 */
int mm_modmin_t(int k, int w, int r)
{
	int t;
	if (k <= r) return k;
	t = r + ((k - r) % w);
	if (t > k) t = k;
	if (t < 1) t = 1;
	return t;
}

static void mm_sketch_modmin(void *km, const char *str, int len, int w, int k, uint32_t rid, mm128_v *p, int r)
{
	int t = mm_modmin_t(k, w, r);
	int wp = w + k - t;
	uint64_t tshift = 2 * (t - 1), tmask = (1ULL<<2*t) - 1, tmer[2] = {0,0};
	uint64_t kshift = 2 * (k - 1), kmask = (1ULL<<2*k) - 1, kmer[2] = {0,0};
	uint64_t dq_h[512];
	int32_t dq_q[512];
	int dq_head = 0, dq_tail = 0;
	mm128_t kbuf[256];
	uint64_t em[4] = {0,0,0,0}; // rolling dedup bitmap over the last 256 k-mer start positions
	int i, j, b, lt = 0;

	assert(len > 0 && (w > 0 && w < 256) && (k > 0 && k <= 28));
	assert(wp > 0 && wp < 512);
	memset(kbuf, 0xff, sizeof(kbuf));
	kv_resize(mm128_t, km, *p, p->n + len/w);

	for (i = 0; i < len; ++i) {
		int c = seq_nt4_table[(uint8_t)str[i]];
		int64_t q, ws, sel;
		uint64_t th;
		int zt;
		if (c > 3) {
			lt = 0; dq_head = dq_tail = 0;
			em[0] = em[1] = em[2] = em[3] = 0;
			tmer[0] = tmer[1] = kmer[0] = kmer[1] = 0;
			continue;
		}
		tmer[0] = (tmer[0] << 2 | c) & tmask;
		tmer[1] = (tmer[1] >> 2) | (3ULL^(uint64_t)c) << tshift;
		kmer[0] = (kmer[0] << 2 | c) & kmask;
		kmer[1] = (kmer[1] >> 2) | (3ULL^(uint64_t)c) << kshift;
		++lt;
		if (lt >= k) {
			int64_t ks = (int64_t)i - k + 1;
			mm128_t info = { UINT64_MAX, UINT64_MAX };
			if (kmer[0] != kmer[1]) {
				int z = kmer[0] < kmer[1]? 0 : 1;
				info.x = hash64(kmer[z], kmask) << 8 | (uint64_t)k;
				info.y = (uint64_t)rid<<32 | (uint32_t)i<<1 | (uint64_t)z;
			}
			kbuf[ks & 255] = info;
		}
		if (lt < t) continue;
		q = (int64_t)i - t + 1;
		zt = tmer[0] < tmer[1]? 0 : 1;
		th = hash64(tmer[zt], tmask);
		while (dq_tail > dq_head && dq_h[(dq_tail-1) & 511] > th) --dq_tail; // strict: keep tied minima
		dq_h[dq_tail & 511] = th;
		dq_q[dq_tail & 511] = (int32_t)q;
		++dq_tail;
		if (lt < w + k - 1) continue;
		ws = (int64_t)i - (w + k - 1) + 1;
		while (dq_head < dq_tail && (int64_t)dq_q[dq_head & 511] < ws) ++dq_head;
		b = (int)((ws - 1) & 255); // ws-1 can no longer be selected; retire its dedup bit
		em[b>>6] &= ~(1ULL << (b&63));
		// emit the k-mers induced by ALL tied minimal t-mers; with t == k (mod w) this
		// set is mirror-symmetric under reverse complement, like the random minimizer
		// with its tie emission; a single (leftmost or rightmost) pick is not
		for (j = dq_head; j < dq_tail && dq_h[j & 511] == dq_h[dq_head & 511]; ++j) {
			sel = ws + (((int64_t)dq_q[j & 511] - ws) % w);
			b = (int)(sel & 255);
			if (!(em[b>>6] >> (b&63) & 1)) {
				mm128_t s = kbuf[sel & 255];
				em[b>>6] |= 1ULL << (b&63);
				if (s.x != UINT64_MAX) kv_push(mm128_t, km, *p, s);
			}
		}
	}
}

void mm_sketch(void *km, const char *str, int len, int w, int k, uint32_t rid, int is_hpc, mm128_v *p)
{
	mm_sketch_random(km, str, len, w, k, rid, is_hpc, p);
}

void mm_sketch_mod(void *km, const char *str, int len, int w, int k, uint32_t rid, int is_hpc, mm128_v *p)
{
	if (!is_hpc && len >= w + k - 1 && mm_modmin_t(k, w, MM_MODMIN_R) < k)
		mm_sketch_modmin(km, str, len, w, k, rid, p, MM_MODMIN_R);
	else
		mm_sketch_random(km, str, len, w, k, rid, is_hpc, p);
}
