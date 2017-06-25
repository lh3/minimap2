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

void mm_write_paf(kstring_t *s, const mm_idx_t *mi, bseq1_t *t, int which, mm_reg1_t *r)
{
	s->l = 0;
	mm_sprintf_lite(s, "%s\t%d\t%d\t%d\t%c\t%s\t%d\t%d\t%d", t->name, t->l_seq, r->qs, r->qe, "+-"[r->rev], mi->seq[r->rid].name, mi->seq[r->rid].len, r->rs, r->re);
	if (r->p) mm_sprintf_lite(s, "\t%d\t%d\t255", r->p->blen - r->p->n_ambi - r->p->n_diff, r->p->blen);
	else mm_sprintf_lite(s, "\t%d\t%d\t255", r->score, r->re - r->rs > r->qe - r->qs? r->re - r->rs : r->qe - r->qs);
	mm_sprintf_lite(s, "\tcm:i:%d", r->cnt);
	if (r->p) mm_sprintf_lite(s, "\ts1:i:%d", r->score);
	if (r->parent == which) mm_sprintf_lite(s, "\ts2:i:%d", r->subsc);
	if (r->p) {
		uint32_t k;
		mm_sprintf_lite(s, "\tNM:i:%d\tAS:i:%d\tnn:i:%d\tcg:Z:", r->p->n_diff, r->p->score, r->p->n_ambi);
		for (k = 0; k < r->p->n_cigar; ++k)
			mm_sprintf_lite(s, "%d%c", r->p->cigar[k]>>4, "MID"[r->p->cigar[k]&0xf]);
	}
}
