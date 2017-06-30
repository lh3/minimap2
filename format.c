#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "mmpriv.h"

static inline void str_enlarge(kstring_t *s, int l)
{
	if (s->l + l + 1 > s->m) {
		s->m = s->l + l + 1;
		kroundup32(s->m);
		s->s = (char*)realloc(s->s, s->m);
	}
}

static inline void str_copy(kstring_t *s, const char *st, const char *en)
{
	str_enlarge(s, en - st);
	memcpy(&s->s[s->l], st, en - st);
	s->l += en - st;
}

static void mm_sprintf_lite(kstring_t *s, const char *fmt, ...)
{
	char buf[16]; // for integer to string conversion
	const char *p, *q;
	va_list ap;
	va_start(ap, fmt);
	for (q = p = fmt; *p; ++p) {
		if (*p == '%') {
			if (p > q) str_copy(s, q, p);
			++p;
			if (*p == 'd') {
				int c, i, l = 0;
				unsigned int x;
				c = va_arg(ap, int);
				x = c >= 0? c : -c;
				do { buf[l++] = x%10 + '0'; x /= 10; } while (x > 0);
				if (c < 0) buf[l++] = '-';
				str_enlarge(s, l);
				for (i = l - 1; i >= 0; --i) s->s[s->l++] = buf[i];
			} else if (*p == 's') {
				char *r = va_arg(ap, char*);
				str_copy(s, r, r + strlen(r));
			} else if (*p == 'c') {
				str_enlarge(s, 1);
				s->s[s->l++] = va_arg(ap, int);
			} else abort();
			q = p + 1;
		}
	}
	if (p > q) str_copy(s, q, p);
	va_end(ap);
	s->s[s->l] = 0;
}

static inline void write_tags(kstring_t *s, const mm_reg1_t *r)
{
	mm_sprintf_lite(s, "\tcm:i:%d", r->cnt);
	if (r->p) mm_sprintf_lite(s, "\ts1:i:%d", r->score);
	if (r->parent == r->id) mm_sprintf_lite(s, "\ts2:i:%d", r->subsc);
	if (r->split) mm_sprintf_lite(s, "\tzd:i:%d", r->split);
	if (r->p) mm_sprintf_lite(s, "\tNM:i:%d\tAS:i:%d\tnn:i:%d", r->p->n_diff, r->p->score, r->p->n_ambi);
}

void mm_write_paf(kstring_t *s, const mm_idx_t *mi, const bseq1_t *t, const mm_reg1_t *r)
{
	s->l = 0;
	mm_sprintf_lite(s, "%s\t%d\t%d\t%d\t%c\t", t->name, t->l_seq, r->qs, r->qe, "+-"[r->rev]);
	if (mi->seq[r->rid].name) mm_sprintf_lite(s, "%s", mi->seq[r->rid].name);
	else mm_sprintf_lite(s, "%d", r->rid);
	mm_sprintf_lite(s, "\t%d\t%d\t%d", mi->seq[r->rid].len, r->rs, r->re);
	if (r->p) mm_sprintf_lite(s, "\t%d\t%d", r->p->blen - r->p->n_ambi - r->p->n_diff, r->p->blen);
	else mm_sprintf_lite(s, "\t%d\t%d", r->score, r->re - r->rs > r->qe - r->qs? r->re - r->rs : r->qe - r->qs);
	mm_sprintf_lite(s, "\t%d", r->mapq);
	write_tags(s, r);
	if (r->p) {
		uint32_t k;
		mm_sprintf_lite(s, "cg:Z:");
		for (k = 0; k < r->p->n_cigar; ++k)
			mm_sprintf_lite(s, "%d%c", r->p->cigar[k]>>4, "MID"[r->p->cigar[k]&0xf]);
	}
}

static char comp_tab[] = {
	  0,   1,	2,	 3,	  4,   5,	6,	 7,	  8,   9,  10,	11,	 12,  13,  14,	15,
	 16,  17,  18,	19,	 20,  21,  22,	23,	 24,  25,  26,	27,	 28,  29,  30,	31,
	 32,  33,  34,	35,	 36,  37,  38,	39,	 40,  41,  42,	43,	 44,  45,  46,	47,
	 48,  49,  50,	51,	 52,  53,  54,	55,	 56,  57,  58,	59,	 60,  61,  62,	63,
	 64, 'T', 'V', 'G', 'H', 'E', 'F', 'C', 'D', 'I', 'J', 'M', 'L', 'K', 'N', 'O',
	'P', 'Q', 'Y', 'S', 'A', 'A', 'B', 'W', 'X', 'R', 'Z',	91,	 92,  93,  94,	95,
	 64, 't', 'v', 'g', 'h', 'e', 'f', 'c', 'd', 'i', 'j', 'm', 'l', 'k', 'n', 'o',
	'p', 'q', 'y', 's', 'a', 'a', 'b', 'w', 'x', 'r', 'z', 123, 124, 125, 126, 127
};

static void sam_write_sq(kstring_t *s, char *seq, int l, int rev, int comp)
{
	if (rev) {
		int i;
		str_enlarge(s, l);
		for (i = 0; i < l; ++i) {
			int c = seq[l - 1 - i];
			s->s[s->l + i] = c < 128 && comp? comp_tab[c] : c;
		}
		s->l += l;
	} else str_copy(s, seq, seq + l);
}

void mm_write_sam(kstring_t *s, const mm_idx_t *mi, const bseq1_t *t, const mm_reg1_t *r)
{
	int flag = 0;
	s->l = 0;
	if (r->rev) flag |= 0x10;
	if (r->parent != r->id) flag |= 0x100;
	if (r->id != 0) flag |= 0x800; // TODO: make sure this is always working!
	mm_sprintf_lite(s, "%s\t%d\t%s\t%d\t%d\t", t->name, flag, mi->seq[r->rid].name, r->rs+1, r->mapq);
	if (r->p) { // TODO: using hard clippings
		uint32_t k, clip_len = r->rev? t->l_seq - r->qe : r->qs;
		int clip_char = (flag&0x800)? 'H' : 'S';
		if (clip_len) mm_sprintf_lite(s, "%d%c", clip_len, clip_char);
		for (k = 0; k < r->p->n_cigar; ++k)
			mm_sprintf_lite(s, "%d%c", r->p->cigar[k]>>4, "MID"[r->p->cigar[k]&0xf]);
		clip_len = r->rev? r->qs : t->l_seq - r->qe;
		if (clip_len) mm_sprintf_lite(s, "%d%c", clip_len, clip_char);
	} else mm_sprintf_lite(s, "*");
	mm_sprintf_lite(s, "\t*\t0\t0\t");
	if ((flag & 0x900) == 0) sam_write_sq(s, t->seq, t->l_seq, r->rev, r->rev);
	else if (flag & 0x100) mm_sprintf_lite(s, "\t*");
	else sam_write_sq(s, t->seq + r->qs, r->qe - r->qs, r->rev, r->rev);
	mm_sprintf_lite(s, "\t*"); // quality
	write_tags(s, r);
}
