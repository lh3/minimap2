#include <stdint.h>
#include "ksw2.h"

#ifdef HAVE_KALLOC
#include "kalloc.h"
#else
#include <stdlib.h>
#define kmalloc(km, size) malloc((size))
#define kcalloc(km, count, size) calloc((count), (size))
#define krealloc(km, ptr, size) realloc((ptr), (size))
#define kfree(km, ptr) free((ptr))
#endif

typedef struct {
	int32_t h, e;
} eh_t;

/*********************
 * One-way extension *
 *********************/

int ksw_extend(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int m, const int8_t *mat, int gapo, int gape, int w,
			   int end_bonus, int zdrop, int h0, int *_qle, int *_tle, int *_gtle, int *_gscore)
{
	eh_t *eh; // score array
	int8_t *qp; // query profile
	int i, j, k, gapoe = gapo + gape, st, en, max, max_i, max_j, max_j0, max_gap, max_ie, gscore;

	// allocate memory
	qp = kmalloc(km, (long)qlen * m);
	eh = kcalloc(km, qlen + 1, 8);

	// generate the query profile
	for (k = i = 0; k < m; ++k) {
		const int8_t *p = &mat[k * m];
		for (j = 0; j < qlen; ++j) qp[i++] = p[query[j]];
	}

	// fill the first row
	eh[0].h = h0, eh[0].e = h0 - gapoe - gapoe > 0? h0 - gapoe - gapoe : 0;
	for (j = 1; j <= qlen && j <= w; ++j) {
		eh[j].h = -(gapo + gape * j), eh[j].e = eh[j-1].e - gape;
		if (eh[j].e < 0) eh[j].e = 0;
		if (eh[j].h < 0) {
			eh[j].h = 0;
			break;
		}
	}

	// adjust $w if it is too large
	k = m * m;
	for (i = 0, max = 0; i < k; ++i) // get the max score
		max = max > mat[i]? max : mat[i];
	max_gap = (int)((double)((qlen < tlen? qlen : tlen) * max + end_bonus - gapo) / gape + 1.);
	max_gap = max_gap > 1? max_gap : 1;
	w = w < max_gap? w : max_gap;

	// DP loop
	max = h0, max_i = max_j = -1; max_ie = -1, gscore = -1; max_j0 = 0;
	st = 0, en = qlen;
	for (i = 0; i < tlen; ++i) {
		int f = 0, h1, m0 = 0;
		int8_t *q = &qp[target[i] * qlen];
		// apply the band and the constraint (if provided)
		if (st < max_j0 - w)     st = max_j0 - w;
		if (en > max_j0 + w + 1) en = max_j0 + w + 1;
		// compute the first column
		if (st == 0) {
			h1 = h0 - (gapo + gape * (i + 1));
			if (h1 < 0) h1 = 0;
		} else h1 = 0;
		for (j = st; j < en; ++j) {
			// At the beginning of the loop: eh[j] = { H(i-1,j-1), E(i,j) }, f = F(i,j) and h1 = H(i,j-1)
			// Similar to SSE2-SW, cells are computed in the following order:
			//   H(i,j)   = max{H(i-1,j-1)+S(i,j), E(i,j), F(i,j)}
			//   E(i+1,j) = max{H(i,j)-gapo, E(i,j)} - gape
			//   F(i,j+1) = max{H(i,j)-gapo, F(i,j)} - gape
			eh_t *p = &eh[j];
			int h = p->h, e = p->e; // get H(i-1,j-1) and E(i-1,j)
			p->h = h1;          // set H(i,j-1) for the next row
			h += q[j];
			h = h > e? h : e;   // e and f are guaranteed to be non-negative, so h>=0 even if h<0
			h = h > f? h : f;
			h1 = h;             // save H(i,j) to h1 for the next column
			max_j0 = m0 > h? max_j0 : j; // record the position where max score is achieved
			m0     = m0 > h? m0     : h; // m0 is stored at eh[mj+1]
			h -= gapoe;
			h = h > 0? h : 0;
			e -= gape;
			e = e > h? e : h;   // computed E(i+1,j)
			p->e = e;           // save E(i+1,j) for the next row
			f -= gape;
			f = f > h? f : h;   // computed F(i,j+1)
		}
		eh[en].h = h1; eh[en].e = 0;
		if (j == qlen) {
			max_ie = gscore > h1? max_ie : i;
			gscore = gscore > h1? gscore : h1;
		}
		if (m0 == 0) break;
		if (m0 > max) {
			max = m0, max_i = i, max_j = max_j0;
		} else if (zdrop > 0) {
			int diff = (i - max_i) - (max_j0 - max_j);
			if (max - m0 - (diff < 0? -diff : diff) * gape > zdrop) break;
		}
		// update beg and end for the next round
		for (j = st; j < en && eh[j].h == 0 && eh[j].e == 0; ++j);
		st = j;
		for (j = en; j >= st && eh[j].h == 0 && eh[j].e == 0; --j);
		en = j + 2 < qlen? j + 2 : qlen;
		//beg = 0; end = qlen; // uncomment this line for debugging
	}
	kfree(km, eh); kfree(km, qp);
	if (_qle) *_qle = max_j + 1;
	if (_tle) *_tle = max_i + 1;
	if (_gtle) *_gtle = max_ie + 1;
	if (_gscore) *_gscore = gscore;
	return max;
}

/********************
 * Global alignment *
 ********************/

#define NEG_INF -0x40000000

static inline uint32_t *push_cigar(void *km, int *n_cigar, int *m_cigar, uint32_t *cigar, int op, int len)
{
	if (*n_cigar == 0 || op != (cigar[(*n_cigar) - 1]&0xf)) {
		if (*n_cigar == *m_cigar) {
			*m_cigar = *m_cigar? (*m_cigar)<<1 : 4;
			cigar = krealloc(km, cigar, (*m_cigar) << 2);
		}
		cigar[(*n_cigar)++] = len<<4 | op;
	} else cigar[(*n_cigar)-1] += len<<4;
	return cigar;
}

int ksw_global(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int m, const int8_t *mat, int gapo, int gape, int w, int *n_cigar_, uint32_t **cigar_)
{
	eh_t *eh;
	int8_t *qp; // query profile
	int32_t i, j, k, max_j = 0, gapoe = gapo + gape, score, n_col, *off = 0, last_en = -1;
	uint8_t *z = 0; // backtrack matrix; in each cell: f<<4|e<<2|h; in principle, we can halve the memory, but backtrack will be more complex

	// allocate memory
	n_col = qlen < 2*w+1? qlen : 2*w+1; // maximum #columns of the backtrack matrix
	qp = kmalloc(km, qlen * m);
	eh = kcalloc(km, qlen + 1, 8);
	if (n_cigar_ && cigar_) {
		*n_cigar_ = 0;
		z = kmalloc(km, (size_t)n_col * tlen);
		off = kcalloc(km, tlen, 4);
	}

	// generate the query profile
	for (k = i = 0; k < m; ++k) {
		const int8_t *p = &mat[k * m];
		for (j = 0; j < qlen; ++j) qp[i++] = p[query[j]];
	}

	// fill the first row
	eh[0].h = 0, eh[0].e = -gapoe - gapoe; // -gapoe*2: one deletion followed by one insertion
	for (j = 1; j <= qlen && j <= w; ++j)
		eh[j].h = -(gapo + gape * j), eh[j].e = eh[j-1].e - gape;
	for (; j <= qlen; ++j) eh[j].h = eh[j].e = NEG_INF; // everything is -inf outside the band

	// DP loop
	for (i = 0; i < tlen; ++i) { // target sequence is in the outer loop
		int32_t f = NEG_INF, h1, st, en, max = NEG_INF;
		int8_t *q = &qp[target[i] * qlen];
		st = max_j > w? max_j - w : 0;
		en = max_j + w + 1 < qlen? max_j + w + 1 : qlen;
		h1 = st > 0? NEG_INF : -gapoe - gape * (i + 1);
		f  = st > 0? NEG_INF : -gapoe - gapoe - gape * i; // similarly, -gapoe*2: one del followed by one ins
		off[i] = st;
		if (n_cigar_ && cigar_) {
			uint8_t *zi = &z[(long)i * n_col];
			for (j = st; j < en; ++j) {
				// At the beginning of the loop: eh[j] = { H(i-1,j-1), E(i,j) }, f = F(i,j) and h1 = H(i,j-1)
				// Cells are computed in the following order:
				//   H(i,j)   = max{H(i-1,j-1) + S(i,j), E(i,j), F(i,j)}
				//   E(i+1,j) = max{H(i,j)-gapo, E(i,j)} - gape
				//   F(i,j+1) = max{H(i,j)-gapo, F(i,j)} - gape
				eh_t *p = &eh[j];
				int32_t h = p->h, e = p->e;
				uint8_t d; // direction
				p->h = h1;
				h += q[j];
				d = h >= e? 0 : 1;
				h = h >= e? h : e;
				d = h >= f? d : 2;
				h = h >= f? h : f;
				h1 = h;
				max_j = max > h? max_j : j;
				max   = max > h? max   : h;
				h -= gapoe;
				e -= gape;
				d |= e > h? 1<<2 : 0;
				e  = e > h? e    : h;
				p->e = e;
				f -= gape;
				d |= f > h? 2<<4 : 0; // if we want to halve the memory, use one bit only, instead of two
				f  = f > h? f    : h;
				zi[j - st] = d; // z[i,j] keeps h for the current cell and e/f for the next cell
			}
		} else {
			for (j = st; j < en; ++j) {
				eh_t *p = &eh[j];
				int32_t h = p->h, e = p->e;
				p->h = h1;
				h += q[j];
				h = h >= e? h : e;
				h = h >= f? h : f;
				h1 = h;
				max_j = max > h? max_j : j;
				max   = max > h? max   : h;
				h -= gapoe;
				e -= gape;
				e  = e > h? e : h;
				p->e = e;
				f -= gape;
				f  = f > h? f : h;
			}
		}
		eh[en].h = h1, eh[en].e = NEG_INF, last_en = en;
	}

	// backtrack
	score = eh[qlen].h;
	if (n_cigar_ && cigar_) {
		int n_cigar = 0, m_cigar = 0, which = 0;
		uint32_t *cigar = 0, tmp;
		i = tlen - 1, k = last_en - 1; // (i,k) points to the last cell; FIXME: with a moving band, we need to take care of last deletion/insertion!!!
		while (i >= 0 && k >= 0) {
			tmp = z[i * n_col + k - off[i]];
			which = tmp >> (which << 1) & 3;
			if (which == 0 && tmp>>6) break;
			if (which == 0) which = tmp & 3;
			if (which == 0)      cigar = push_cigar(km, &n_cigar, &m_cigar, cigar, 0, 1), --i, --k; // match
			else if (which == 1) cigar = push_cigar(km, &n_cigar, &m_cigar, cigar, 2, 1), --i;      // deletion
			else                 cigar = push_cigar(km, &n_cigar, &m_cigar, cigar, 1, 1), --k;      // insertion
		}
		if (i >= 0) cigar = push_cigar(km, &n_cigar, &m_cigar, cigar, 2, i + 1); // first deletion
		if (k >= 0) cigar = push_cigar(km, &n_cigar, &m_cigar, cigar, 1, k + 1); // first insertion
		for (i = 0; i < n_cigar>>1; ++i) // reverse CIGAR
			tmp = cigar[i], cigar[i] = cigar[n_cigar-1-i], cigar[n_cigar-1-i] = tmp;
		*n_cigar_ = n_cigar, *cigar_ = cigar;
	}

	kfree(km, qp); kfree(km, eh);
	if (n_cigar_ && cigar_) {
		kfree(km, z);
		kfree(km, off);
	}
	return score;
}
