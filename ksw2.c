#include <stdint.h>
#include <stdio.h> // for debugging only
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

void ksw_gen_simple_mat(int m, int8_t *mat, int8_t a, int8_t b)
{
	int i, j;
	a = a < 0? -a : a;
	b = b > 0? -b : b;
	for (i = 0; i < m - 1; ++i) {
		for (j = 0; j < m - 1; ++j)
			mat[i * m + j] = i == j? a : b;
		mat[i * m + m - 1] = 0;
	}
	for (j = 0; j < m; ++j)
		mat[(m - 1) * m + j] = 0;
}

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

int ksw_global(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat, int8_t gapo, int8_t gape, int w, int *n_cigar_, uint32_t **cigar_)
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
	eh[0].h = 0, eh[0].e = -gapoe - gapo;
	for (j = 1; j <= qlen && j <= w; ++j)
		eh[j].h = -(gapo + gape * j), eh[j].e = -(gapoe + gapo + gape * j);
	for (; j <= qlen; ++j) eh[j].h = eh[j].e = NEG_INF; // everything is -inf outside the band

	// DP loop
	for (i = 0; i < tlen; ++i) { // target sequence is in the outer loop
		int32_t f = NEG_INF, h1, st, en, max = NEG_INF;
		int8_t *q = &qp[target[i] * qlen];
		#if 0
		st = max_j > w? max_j - w : 0;
		en = max_j + w + 1 < qlen? max_j + w + 1 : qlen;
		#else
		st = i > w? i - w : 0;
		en = i + w + 1 < qlen? i + w + 1 : qlen;
		#endif
		h1 = st > 0? NEG_INF : -(gapo + gape * i);
		f  = st > 0? NEG_INF : -(gapoe + gapo + gape * i);
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

#ifdef __SSE2__
#include <emmintrin.h>

#ifdef __SSE4_1__
#include <smmintrin.h>
#endif

int ksw_global2_sse(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat, int8_t q, int8_t e, int w, int *n_cigar_, uint32_t **cigar_)
{
	int r, t, n_col, *off, score, tlen16;
	int8_t *u, *v, *x, *y, *s;
	uint8_t *p, *qr, *mem;
	__m128i q_, qe2_, zero_, flag1_, flag2_, flag4_, flag32_;

	zero_   = _mm_set1_epi8(0);
	q_      = _mm_set1_epi8(q);
	qe2_    = _mm_set1_epi8((q + e) * 2);
	flag1_  = _mm_set1_epi8(1<<0);
	flag2_  = _mm_set1_epi8(2<<0);
	flag4_  = _mm_set1_epi8(1<<2);
	flag32_ = _mm_set1_epi8(2<<4);

	w = (w + 1 + 15) / 16 * 16 - 1;
	tlen16 = (tlen + 15) / 16 * 16;
	n_col = w + 1 < tlen16? w + 1 : tlen16; // number of columns in the backtrack matrix

	mem = (uint8_t*)kcalloc(km, tlen16 * 5 + 15, 1);
	u = (int8_t*)(((size_t)mem + 15) >> 4 << 4); // 16-byte aligned (though not necessary)
	v = u + tlen16, x = v + tlen16, y = x + tlen16, s = y + tlen16;
	qr = (uint8_t*)kcalloc(km, qlen, 1);
	p = (uint8_t*)kcalloc(km, (qlen + tlen) * n_col, 1);
	off = (int*)kmalloc(km, (qlen + tlen) * sizeof(int));

	for (t = 0; t < qlen; ++t)
		qr[t] = query[qlen - 1 - t];

	for (r = 0; r < qlen + tlen - 1; ++r) {
		int st = 0, en = tlen - 1;
		int8_t x1, v1;
		uint8_t *pr = p + r * n_col;
		__m128i x1_, v1_;
		// find the boundaries
		if (st < r - qlen + 1) st = r - qlen + 1;
		if (en > r) en = r;
		if (st < (r-w+1)>>1) st = (r-w+1)>>1; // take the ceil
		if (en > (r+w)>>1) en = (r+w)>>1; // take the floor
		off[r] = st;
		// set boundary conditions
		if (st != 0) {
			if (r > st + st + w - 1) x1 = v1 = 0;
			else x1 = x[st-1], v1 = v[st-1]; // (r-1, st-1) in the band
		} else x1 = 0, v1 = r? q : 0;
		if (en != r) {
			if (r < en + en - w - 1) y[en] = u[en] = 0; // (r-1,en) out of the band; TODO: is this line necessary?
		} else y[r] = 0, u[r] = r? q : 0;
		// loop fission: set scores first
		for (t = st; t <= en; ++t)
			s[t] = mat[target[t] * m + qr[t + qlen - 1 - r]];
		// core loop
		x1_ = _mm_cvtsi32_si128(x1);
		v1_ = _mm_cvtsi32_si128(v1);
		for (t = st; t <= en; t += 16) {
			__m128i d, z, a, b, xt1, vt1, ut, tmp;

			z = _mm_add_epi8(_mm_loadu_si128((__m128i*)&s[t]), qe2_);

			xt1 = _mm_loadu_si128((__m128i*)&x[t]);          // xt1 <- x[r-1][t..t+15]
			tmp = _mm_srli_si128(xt1, 15);                   // tmp <- x[r-1][t+15]
			xt1 = _mm_or_si128(_mm_slli_si128(xt1, 1), x1_); // xt1 <- x[r-1][t-1..t+14]
			x1_ = tmp;
			vt1 = _mm_loadu_si128((__m128i*)&v[t]);          // vt1 <- v[r-1][t..t+15]
			tmp = _mm_srli_si128(vt1, 15);                   // tmp <- v[r-1][t+15]
			vt1 = _mm_or_si128(_mm_slli_si128(vt1, 1), v1_); // vt1 <- v[r-1][t-1..t+14]
			v1_ = tmp;
			a = _mm_add_epi8(xt1, vt1);                      // a <- x[r-1][t-1..t+14] + v[r-1][t-1..t+14]

			ut = _mm_loadu_si128((__m128i*)&u[t]);           // ut <- u[t..t+15]
			b = _mm_add_epi8(_mm_loadu_si128((__m128i*)&y[t]), ut); // b <- y[r-1][t..t+15] + u[r-1][t..t+15]

			d = _mm_and_si128(flag1_, _mm_cmplt_epi8(a, z)); // d = a < z? 1 : 0
#ifdef __SSE4_1__
			z = _mm_max_epi8(z, a);                          // z = z > a? z : a (signed)
			tmp = _mm_cmplt_epi8(b, z);
			d = _mm_blendv_epi8(d, flag2_, tmp);             // d = b < z? d : 2
#else
			z = _mm_and_si128(z, _mm_cmpgt_epi8(z, zero_));  // z = z > 0? z : 0;
			z = _mm_max_epu8(z, a);                          // z = max(z, a); this works because both are non-negative
			tmp = _mm_cmplt_epi8(b, z);
			d = _mm_or_si128(_mm_andnot_si128(tmp, d), _mm_and_si128(tmp, flag2_));
#endif
			z = _mm_max_epu8(z, b);                          // z = max(z, b); this works because both are non-negative

			_mm_storeu_si128((__m128i*)&u[t], _mm_sub_epi8(z, vt1)); // u[r][t..t+15] <- z - v[r-1][t-1..t+14]
			_mm_storeu_si128((__m128i*)&v[t], _mm_sub_epi8(z, ut));  // v[r][t..t+15] <- z - u[r-1][t..t+15]

			z = _mm_sub_epi8(z, q_);
			a = _mm_sub_epi8(a, z);
			b = _mm_sub_epi8(b, z);
			tmp = _mm_cmpgt_epi8(a, zero_);
			d = _mm_or_si128(d, _mm_and_si128(flag4_,  tmp));
			_mm_storeu_si128((__m128i*)&x[t], _mm_and_si128(a, tmp));
			tmp = _mm_cmpgt_epi8(b, zero_);
			d = _mm_or_si128(d, _mm_and_si128(flag32_, tmp));
			_mm_storeu_si128((__m128i*)&y[t], _mm_and_si128(b, tmp));
			_mm_storeu_si128((__m128i*)&pr[t - st], d);
		}
	}
	kfree(km, mem); kfree(km, qr);
	{ // backtrack
		int n_cigar = 0, m_cigar = 0, which = 0, i, j, k, l;
		uint32_t *cigar = 0, tmp;
		i = tlen - 1, j = qlen - 1;
		while (i >= 0 && j >= 0) {
			r = i + j;
			tmp = p[r * n_col + i - off[r]];
			which = tmp >> (which << 1) & 3;
			if (which == 0 && tmp>>6) break;
			if (which == 0) which = tmp & 3;
			if (which == 0)      cigar = push_cigar(km, &n_cigar, &m_cigar, cigar, 0, 1), --i, --j; // match
			else if (which == 1) cigar = push_cigar(km, &n_cigar, &m_cigar, cigar, 2, 1), --i;      // deletion
			else                 cigar = push_cigar(km, &n_cigar, &m_cigar, cigar, 1, 1), --j;      // insertion
		}
		if (i >= 0) cigar = push_cigar(km, &n_cigar, &m_cigar, cigar, 2, i + 1); // first deletion
		if (j >= 0) cigar = push_cigar(km, &n_cigar, &m_cigar, cigar, 1, j + 1); // first insertion
		for (i = 0; i < n_cigar>>1; ++i) // reverse CIGAR
			tmp = cigar[i], cigar[i] = cigar[n_cigar-1-i], cigar[n_cigar-1-i] = tmp;
		*n_cigar_ = n_cigar, *cigar_ = cigar;

		// compute score
		for (k = 0, score = 0, i = j = 0; k < n_cigar; ++k) {
			int op = cigar[k] & 0xf, len = cigar[k] >> 4;
			if (op == 0) {
				for (l = 0; l < len; ++l)
					score += mat[target[i + l] * m + query[j + l]];
				i += len, j += len;
			} else if (op == 1) {
				score -= q + len * e;
				j += len;
			} else if (op == 2) {
				score -= q + len * e;
				i += len;
			}
		}
	}
	kfree(km, p); kfree(km, off);
	return score;
}
#endif // __SSE2__

int ksw_global2(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat, int8_t q, int8_t e, int w, int *n_cigar_, uint32_t **cigar_)
{
	int qe = q + e, qe2 = qe + qe, r, t, n_col, *off, score;
	int8_t *u, *v, *x, *y, *s;
	uint8_t *p, *qr;

	u = (int8_t*)kcalloc(km, tlen + 1, 1);
	v = (int8_t*)kcalloc(km, tlen + 1, 1);
	x = (int8_t*)kcalloc(km, tlen + 1, 1);
	y = (int8_t*)kcalloc(km, tlen + 1, 1);
	s = (int8_t*)kmalloc(km, tlen);
	qr = (uint8_t*)kmalloc(km, qlen);
	n_col = w + 1 < tlen? w + 1 : tlen;
	p = (uint8_t*)kcalloc(km, (qlen + tlen) * n_col, 1);
	off = (int*)kmalloc(km, (qlen + tlen) * sizeof(int));

	for (t = 0; t < qlen; ++t)
		qr[t] = query[qlen - 1 - t];

	for (r = 0; r < qlen + tlen - 1; ++r) {
		int st = 0, en = tlen - 1;
		int8_t x1, v1;
		uint8_t *pr = p + r * n_col;
		// find the boundaries
		if (st < r - qlen + 1) st = r - qlen + 1;
		if (en > r) en = r;
		if (st < (r-w+1)>>1) st = (r-w+1)>>1; // take the ceil
		if (en > (r+w)>>1) en = (r+w)>>1; // take the floor
		off[r] = st;
		// set boundary conditions
		if (st != 0) {
			if (r > st + st + w - 1) x1 = v1 = 0;
			else x1 = x[st-1], v1 = v[st-1]; // (r-1, st-1) in the band
		} else x1 = 0, v1 = r? q : 0;
		if (en != r) {
			if (r < en + en - w - 1) y[en] = u[en] = 0; // (r-1,en) out of the band; TODO: is this line necessary?
		} else y[r] = 0, u[r] = r? q : 0;
		// loop fission: set scores first
		for (t = st; t <= en; ++t)
			s[t] = mat[target[t] * m + qr[t + qlen - 1 - r]];
		// core loop
		for (t = st; t <= en; ++t) {
			/* At the beginning of the loop, v1=v(r-1,t-1), x1=x(r-1,t-1), u[t]=u(r-1,t), v[t]=v(r-1,t), x[t]=x(r-1,t), y[t]=y(r-1,t)
			   a      = x(r-1,t-1) + v(r-1,t-1)
			   b      = y(r-1,t)   + u(r-1,t)
			   z      = max{ S(t,r-t) + 2q + 2r, a, b }
			   u(r,t) = z - v(r-1,t-1)
			   v(r,t) = z - u(r-1,t)
			   x(r,t) = max{ 0, a - z + q }
			   y(r,t) = max{ 0, b - z + q }
			*/
			uint8_t d;
			int8_t u1;
			int8_t z = s[t] + qe2;
			int8_t a = x1   + v1;
			int8_t b = y[t] + u[t];
			d = z >= a? 0 : 1;
			z = z >= a? z : a;
			d = z >= b? d : 2;
			z = z >= b? z : b;
			u1 = u[t];              // u1   = u(r-1,t) (temporary variable)
			u[t] = z - v1;          // u[t] = u(r,t)
			v1 = v[t];              // v1   = v(r-1,t) (set for the next iteration)
			v[t] = z - u1;          // v[t] = v(r,t)
			z -= q;
			a -= z;
			b -= z;
			x1 = x[t];              // x1   = x(r-1,t) (set for the next iteration)
			d   |= a > 0? 1<<2 : 0;
			x[t] = a > 0? a    : 0; // x[t] = x(r,t)
			d   |= b > 0? 2<<4 : 0;
			y[t] = b > 0? b    : 0; // y[t] = y(r,t)
			pr[t - st] = d;
		}
	}
	kfree(km, u); kfree(km, v); kfree(km, x); kfree(km, y); kfree(km, s); kfree(km, qr);
	{ // backtrack
		int n_cigar = 0, m_cigar = 0, which = 0, i, j, k, l;
		uint32_t *cigar = 0, tmp;
		i = tlen - 1, j = qlen - 1;
		while (i >= 0 && j >= 0) {
			r = i + j;
			tmp = p[r * n_col + i - off[r]];
			which = tmp >> (which << 1) & 3;
			if (which == 0 && tmp>>6) break;
			if (which == 0) which = tmp & 3;
			if (which == 0)      cigar = push_cigar(km, &n_cigar, &m_cigar, cigar, 0, 1), --i, --j; // match
			else if (which == 1) cigar = push_cigar(km, &n_cigar, &m_cigar, cigar, 2, 1), --i;      // deletion
			else                 cigar = push_cigar(km, &n_cigar, &m_cigar, cigar, 1, 1), --j;      // insertion
		}
		if (i >= 0) cigar = push_cigar(km, &n_cigar, &m_cigar, cigar, 2, i + 1); // first deletion
		if (j >= 0) cigar = push_cigar(km, &n_cigar, &m_cigar, cigar, 1, j + 1); // first insertion
		for (i = 0; i < n_cigar>>1; ++i) // reverse CIGAR
			tmp = cigar[i], cigar[i] = cigar[n_cigar-1-i], cigar[n_cigar-1-i] = tmp;
		*n_cigar_ = n_cigar, *cigar_ = cigar;

		// compute score
		for (k = 0, score = 0, i = j = 0; k < n_cigar; ++k) {
			int op = cigar[k] & 0xf, len = cigar[k] >> 4;
			if (op == 0) {
				for (l = 0; l < len; ++l)
					score += mat[target[i + l] * m + query[j + l]];
				i += len, j += len;
			} else if (op == 1) {
				score -= q + len * e;
				j += len;
			} else if (op == 2) {
				score -= q + len * e;
				i += len;
			}
		}
	}
	kfree(km, p); kfree(km, off);
	return score;
}
