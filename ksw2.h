#ifndef KSW2_H_
#define KSW2_H_

#include <stdint.h>

#define KSW_NEG_INF -0x40000000

#define KSW_EZ_SCORE_ONLY  0x01 // don't record alignment path/cigar
#define KSW_EZ_RIGHT       0x02 // right-align gaps
#define KSW_EZ_GENERIC_SC  0x04 // without this flag: match/mismatch only; last symbol is a wildcard
#define KSW_EZ_APPROX_MAX  0x08 // approximate max; this is faster with sse
#define KSW_EZ_APPROX_DROP 0x10 // approximate Z-drop; faster with sse
#define KSW_EZ_DYN_BAND    0x20 // once used, ksw_extz_t::{mqe,mte} may be wrong
#define KSW_EZ_EXTZ_ONLY   0x40 // only perform extension
#define KSW_EZ_REV_CIGAR   0x80 // reverse CIGAR in the output

typedef struct {
	uint32_t max:31, zdropped:1;
	int max_q, max_t;      // max extension coordinate
	int mqe, mqe_t;        // max score when reaching the end of query
	int mte, mte_q;        // max score when reaching the end of target
	int score;             // max score reaching both ends; may be KSW_NEG_INF
	int m_cigar, n_cigar;
	uint32_t *cigar;
} ksw_extz_t;

#ifdef __cplusplus
extern "C" {
#endif

/**
 * NW-like extension
 *
 * @param km        memory pool, when used with kalloc
 * @param qlen      query length
 * @param query     query sequence with 0 <= query[i] < m
 * @param tlen      target length
 * @param target    target sequence with 0 <= target[i] < m
 * @param m         number of residue types
 * @param mat       m*m scoring mattrix in one-dimension array
 * @param gapo      gap open penalty; a gap of length l cost "-(gapo+l*gape)"
 * @param gape      gap extension penalty
 * @param w         band width
 * @param zdrop     off-diagonal drop-off to stop extension (positive)
 * @param flag      flag (see KSW_EZ_* macros)
 * @param ez        (out) scores and cigar
 */
void ksw_extz(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat, int8_t q, int8_t e, int w, int zdrop, int flag, ksw_extz_t *ez);
void ksw_extz2_sse(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat, int8_t q, int8_t e, int w, int zdrop, int flag, ksw_extz_t *ez);

/**
 * Global alignment
 *
 * (first 10 parameters identical to ksw_extz_sse())
 * @param m_cigar   (modified) max CIGAR length; feed 0 if cigar==0
 * @param n_cigar   (out) number of CIGAR elements
 * @param cigar     (out) BAM-encoded CIGAR; caller need to deallocate with kfree(km, )
 *
 * @return          score of the alignment
 */
int ksw_gg(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat, int8_t gapo, int8_t gape, int w, int *m_cigar_, int *n_cigar_, uint32_t **cigar_);
int ksw_gg2(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat, int8_t gapo, int8_t gape, int w, int *m_cigar_, int *n_cigar_, uint32_t **cigar_);
int ksw_gg2_sse(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat, int8_t gapo, int8_t gape, int w, int *m_cigar_, int *n_cigar_, uint32_t **cigar_);

#ifdef __cplusplus
}
#endif

/************************************
 *** Private macros and functions ***
 ************************************/

#ifdef HAVE_KALLOC
#include "kalloc.h"
#else
#include <stdlib.h>
#define kmalloc(km, size) malloc((size))
#define kcalloc(km, count, size) calloc((count), (size))
#define krealloc(km, ptr, size) realloc((ptr), (size))
#define kfree(km, ptr) free((ptr))
#endif

static inline uint32_t *ksw_push_cigar(void *km, int *n_cigar, int *m_cigar, uint32_t *cigar, uint32_t op, int len)
{
	if (*n_cigar == 0 || op != (cigar[(*n_cigar) - 1]&0xf)) {
		if (*n_cigar == *m_cigar) {
			*m_cigar = *m_cigar? (*m_cigar)<<1 : 4;
			cigar = (uint32_t*)krealloc(km, cigar, (*m_cigar) << 2);
		}
		cigar[(*n_cigar)++] = len<<4 | op;
	} else cigar[(*n_cigar)-1] += len<<4;
	return cigar;
}

static inline void ksw_backtrack(void *km, int is_rot, int is_rev, const uint8_t *p, const int *off, int n_col, int i0, int j0, int *m_cigar_, int *n_cigar_, uint32_t **cigar_)
{
	int n_cigar = 0, m_cigar = *m_cigar_, which = 0, i = i0, j = j0, r;
	uint32_t *cigar = *cigar_, tmp;
	while (i >= 0 && j >= 0) {
		if (is_rot) r = i + j, tmp = p[r * n_col + i - off[r]];
		else tmp = p[i * n_col + j - off[i]];
		which = tmp >> (which << 1) & 3;
		if (which == 0) which = tmp & 3;
		if (which == 0)      cigar = ksw_push_cigar(km, &n_cigar, &m_cigar, cigar, 0, 1), --i, --j; // match
		else if (which == 1) cigar = ksw_push_cigar(km, &n_cigar, &m_cigar, cigar, 2, 1), --i;      // deletion
		else                 cigar = ksw_push_cigar(km, &n_cigar, &m_cigar, cigar, 1, 1), --j;      // insertion
	}
	if (i >= 0) cigar = ksw_push_cigar(km, &n_cigar, &m_cigar, cigar, 2, i + 1); // first deletion
	if (j >= 0) cigar = ksw_push_cigar(km, &n_cigar, &m_cigar, cigar, 1, j + 1); // first insertion
	if (!is_rev)
		for (i = 0; i < n_cigar>>1; ++i) // reverse CIGAR
			tmp = cigar[i], cigar[i] = cigar[n_cigar-1-i], cigar[n_cigar-1-i] = tmp;
	*m_cigar_ = m_cigar, *n_cigar_ = n_cigar, *cigar_ = cigar;
}

static inline void ksw_reset_extz(ksw_extz_t *ez)
{
	ez->max_q = ez->max_t = ez->mqe_t = ez->mte_q = -1;
	ez->max = 0, ez->score = ez->mqe = ez->mte = KSW_NEG_INF;
	ez->n_cigar = 0, ez->zdropped = 0;
}

#endif
