#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdio.h>
#include "kalloc.h"
#include "mmpriv.h"

static char mm_rg_id[256];

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

static char *mm_escape(char *s)
{
	char *p, *q;
	for (p = q = s; *p; ++p) {
		if (*p == '\\') {
			++p;
			if (*p == 't') *q++ = '\t';
			else if (*p == '\\') *q++ = '\\';
		} else *q++ = *p;
	}
	*q = '\0';
	return s;
}

static void sam_write_rg_line(kstring_t *str, const char *s)
{
	char *p, *q, *r, *rg_line = 0;
	memset(mm_rg_id, 0, 256);
	if (s == 0) return;
	if (strstr(s, "@RG") != s) {
		if (mm_verbose >= 1) fprintf(stderr, "[ERROR] the read group line is not started with @RG\n");
		goto err_set_rg;
	}
	if (strstr(s, "\t") != NULL) {
		if (mm_verbose >= 1) fprintf(stderr, "[ERROR] the read group line contained literal <tab> characters -- replace with escaped tabs: \\t\n");
		goto err_set_rg;
	}
	rg_line = strdup(s);
	mm_escape(rg_line);
	if ((p = strstr(rg_line, "\tID:")) == 0) {
		if (mm_verbose >= 1) fprintf(stderr, "[ERROR] no ID within the read group line\n");
		goto err_set_rg;
	}
	p += 4;
	for (q = p; *q && *q != '\t' && *q != '\n'; ++q);
	if (q - p + 1 > 256) {
		if (mm_verbose >= 1) fprintf(stderr, "[ERROR] @RG:ID is longer than 255 characters\n");
		goto err_set_rg;
	}
	for (q = p, r = mm_rg_id; *q && *q != '\t' && *q != '\n'; ++q)
		*r++ = *q;
	mm_sprintf_lite(str, "%s\n", rg_line);

err_set_rg:
	free(rg_line);
}

void mm_write_sam_hdr_no_SQ(const char *rg, const char *ver, int argc, char *argv[])
{
	kstring_t str = {0,0,0};
	sam_write_rg_line(&str, rg);
	mm_sprintf_lite(&str, "@PG\tID:minimap2\tPN:minimap2");
	if (ver) mm_sprintf_lite(&str, "\tVN:%s", ver);
	if (argc > 1) {
		int i;
		mm_sprintf_lite(&str, "\tCL:minimap2");
		for (i = 1; i < argc; ++i)
			mm_sprintf_lite(&str, " %s", argv[i]);
	}
	mm_sprintf_lite(&str, "\n");
	fputs(str.s, stdout);
	free(str.s);
}

static void write_cs(void *km, kstring_t *s, const mm_idx_t *mi, const mm_bseq1_t *t, const mm_reg1_t *r)
{
	extern unsigned char seq_nt4_table[256];
	int i, q_off, t_off;
	uint8_t *qseq, *tseq;
	char *tmp;
	if (r->p == 0) return;
	mm_sprintf_lite(s, "\tcs:Z:");
	qseq = (uint8_t*)kmalloc(km, r->qe - r->qs);
	tseq = (uint8_t*)kmalloc(km, r->re - r->rs);
	tmp = (char*)kmalloc(km, r->re - r->rs > r->qe - r->qs? r->re - r->rs + 1 : r->qe - r->qs + 1);
	mm_idx_getseq(mi, r->rid, r->rs, r->re, tseq);
	if (!r->rev) {
		for (i = r->qs; i < r->qe; ++i)
			qseq[i - r->qs] = seq_nt4_table[(uint8_t)t->seq[i]];
	} else {
		for (i = r->qs; i < r->qe; ++i) {
			uint8_t c = seq_nt4_table[(uint8_t)t->seq[i]];
			qseq[r->qe - i - 1] = c >= 4? 4 : 3 - c;
		}
	}
	for (i = q_off = t_off = 0; i < r->p->n_cigar; ++i) {
		int j, op = r->p->cigar[i]&0xf, len = r->p->cigar[i]>>4;
		assert(op >= 0 && op <= 2);
		if (op == 0) {
			int l_tmp = 0;
			for (j = 0; j < len; ++j) {
				if (qseq[q_off + j] != tseq[t_off + j]) {
					if (l_tmp > 0) {
						tmp[l_tmp] = 0;
						mm_sprintf_lite(s, "=%s", tmp);
						l_tmp = 0;
					}
					mm_sprintf_lite(s, "*%c%c", "acgtn"[tseq[t_off + j]], "acgtn"[qseq[q_off + j]]);
				} else tmp[l_tmp++] = "ACGTN"[qseq[q_off + j]];
			}
			if (l_tmp > 0) {
				tmp[l_tmp] = 0;
				mm_sprintf_lite(s, "=%s", tmp);
			}
			q_off += len, t_off += len;
		} else if (op == 1) {
			for (j = 0, tmp[len] = 0; j < len; ++j)
				tmp[j] = "acgtn"[qseq[q_off + j]];
			mm_sprintf_lite(s, "+%s", tmp);
			q_off += len;
		} else if (op == 2) {
			for (j = 0, tmp[len] = 0; j < len; ++j)
				tmp[j] = "acgtn"[tseq[t_off + j]];
			mm_sprintf_lite(s, "-%s", tmp);
			t_off += len;
		}
	}
	assert(t_off == r->re - r->rs && q_off == r->qe - r->qs);
	kfree(km, qseq); kfree(km, tseq); kfree(km, tmp);
}

static inline void write_tags(kstring_t *s, const mm_reg1_t *r)
{
	int type = r->inv? 'I' : r->id == r->parent? 'P' : 'S';
	mm_sprintf_lite(s, "\ttp:A:%c\tcm:i:%d\ts1:i:%d", type, r->cnt, r->score);
	if (r->parent == r->id) mm_sprintf_lite(s, "\ts2:i:%d", r->subsc);
	if (r->split) mm_sprintf_lite(s, "\tzd:i:%d", r->split);
	if (r->p) {
		mm_sprintf_lite(s, "\tNM:i:%d\tms:i:%d\tAS:i:%d\tnn:i:%d", r->p->n_diff, r->p->dp_max, r->p->dp_score, r->p->n_ambi);
		if (r->p->trans_strand == 1 || r->p->trans_strand == 2)
			mm_sprintf_lite(s, "\tts:A:%c", "?+-?"[r->p->trans_strand]);
	}
}

void mm_write_paf(kstring_t *s, const mm_idx_t *mi, const mm_bseq1_t *t, const mm_reg1_t *r, void *km, int opt_flag)
{
	s->l = 0;
	mm_sprintf_lite(s, "%s\t%d\t%d\t%d\t%c\t", t->name, t->l_seq, r->qs, r->qe, "+-"[r->rev]);
	if (mi->seq[r->rid].name) mm_sprintf_lite(s, "%s", mi->seq[r->rid].name);
	else mm_sprintf_lite(s, "%d", r->rid);
	mm_sprintf_lite(s, "\t%d\t%d\t%d", mi->seq[r->rid].len, r->rs, r->re);
	if (r->p) mm_sprintf_lite(s, "\t%d\t%d", r->p->blen - r->p->n_ambi - r->p->n_diff, r->p->blen);
	else mm_sprintf_lite(s, "\t%d\t%d", r->fuzzy_mlen, r->fuzzy_blen);
	mm_sprintf_lite(s, "\t%d", r->mapq);
	write_tags(s, r);
	if (r->p && (opt_flag & MM_F_OUT_CG)) {
		uint32_t k;
		mm_sprintf_lite(s, "\tcg:Z:");
		for (k = 0; k < r->p->n_cigar; ++k)
			mm_sprintf_lite(s, "%d%c", r->p->cigar[k]>>4, "MIDN"[r->p->cigar[k]&0xf]);
	}
	if (r->p && (opt_flag & MM_F_OUT_CS))
		write_cs(km, s, mi, t, r);
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

void mm_write_sam_SQ(const mm_idx_t *idx)
{
	uint32_t i;
	for (i = 0; i < idx->n_seq; ++i)
		printf("@SQ\tSN:%s\tLN:%d\n", idx->seq[i].name, idx->seq[i].len);
}

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

void mm_write_sam(kstring_t *s, const mm_idx_t *mi, const mm_bseq1_t *t, const mm_reg1_t *r, int n_regs, const mm_reg1_t *regs)
{
	int flag = 0;
	s->l = 0;
	if (r == 0) {
		mm_sprintf_lite(s, "%s\t4\t*\t0\t0\t*\t*\t0\t0\t", t->name);
		sam_write_sq(s, t->seq, t->l_seq, 0, 0);
		mm_sprintf_lite(s, "\t");
		if (t->qual) sam_write_sq(s, t->qual, t->l_seq, 0, 0);
		else mm_sprintf_lite(s, "*");
	} else {
		if (r->rev) flag |= 0x10;
		if (r->parent != r->id) flag |= 0x100;
		else if (!r->sam_pri) flag |= 0x800;
		mm_sprintf_lite(s, "%s\t%d\t%s\t%d\t%d\t", t->name, flag, mi->seq[r->rid].name, r->rs+1, r->mapq);
		if (r->p) { // actually this should always be true for SAM output
			uint32_t k, clip_len = r->rev? t->l_seq - r->qe : r->qs;
			int clip_char = (flag&0x800)? 'H' : 'S';
			if (clip_len) mm_sprintf_lite(s, "%d%c", clip_len, clip_char);
			for (k = 0; k < r->p->n_cigar; ++k)
				mm_sprintf_lite(s, "%d%c", r->p->cigar[k]>>4, "MIDN"[r->p->cigar[k]&0xf]);
			clip_len = r->rev? r->qs : t->l_seq - r->qe;
			if (clip_len) mm_sprintf_lite(s, "%d%c", clip_len, clip_char);
		} else mm_sprintf_lite(s, "*");
		mm_sprintf_lite(s, "\t*\t0\t0\t");
		if ((flag & 0x900) == 0) {
			sam_write_sq(s, t->seq, t->l_seq, r->rev, r->rev);
			mm_sprintf_lite(s, "\t");
			if (t->qual) sam_write_sq(s, t->qual, t->l_seq, r->rev, 0);
			else mm_sprintf_lite(s, "*");
		} else if (flag & 0x100) {
			mm_sprintf_lite(s, "*\t*");
		} else {
			sam_write_sq(s, t->seq + r->qs, r->qe - r->qs, r->rev, r->rev);
			mm_sprintf_lite(s, "\t");
			if (t->qual) sam_write_sq(s, t->qual + r->qs, r->qe - r->qs, r->rev, 0);
			else mm_sprintf_lite(s, "*");
		}
		write_tags(s, r);
		if (mm_rg_id[0]) mm_sprintf_lite(s, "\tRG:Z:%s", mm_rg_id);
		if (r->parent == r->id && r->p && n_regs > 1 && regs && r >= regs && r - regs < n_regs) { // supplementary aln may exist
			int i, n_sa = 0; // n_sa: number of SA fields
			for (i = 0; i < n_regs; ++i)
				if (i != r - regs && regs[i].parent == regs[i].id && regs[i].p)
					++n_sa;
			if (n_sa > 0) {
				mm_sprintf_lite(s, "\tSA:Z:");
				for (i = 0; i < n_regs; ++i) {
					const mm_reg1_t *q = &regs[i];
					int l_M, l_I = 0, l_D = 0, clip5 = 0, clip3 = 0;
					if (r == q || q->parent != q->id || q->p == 0) continue;
					if (q->qe - q->qs < q->re - q->rs) l_M = q->qe - q->qs, l_D = (q->re - q->rs) - l_M;
					else l_M = q->re - q->rs, l_I = (q->qe - q->qs) - l_M;
					clip5 = q->rev? t->l_seq - q->qe : q->qs;
					clip3 = q->rev? q->qs : t->l_seq - q->qe;
					mm_sprintf_lite(s, "%s,%d,%c,", mi->seq[q->rid].name, q->rs+1, "+-"[q->rev]);
					if (clip5) mm_sprintf_lite(s, "%dS", clip5);
					if (l_M) mm_sprintf_lite(s, "%dM", l_M);
					if (l_I) mm_sprintf_lite(s, "%dI", l_I);
					if (l_D) mm_sprintf_lite(s, "%dD", l_D);
					if (clip3) mm_sprintf_lite(s, "%dS", clip3);
					mm_sprintf_lite(s, ",%d,%d;", q->mapq, q->p->n_diff);
				}
			}
		}
	}
	s->s[s->l] = 0; // we always have room for an extra byte (see str_enlarge)
}
