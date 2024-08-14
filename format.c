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
			} else if (*p == 'u') {
				int i, l = 0;
				uint32_t x;
				x = va_arg(ap, uint32_t);
				do { buf[l++] = x%10 + '0'; x /= 10; } while (x > 0);
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

static int sam_write_rg_line(kstring_t *str, const char *s)
{
	char *p, *q, *r, *rg_line = 0;
	memset(mm_rg_id, 0, 256);
	if (s == 0) return 0;
	if (strstr(s, "@RG") != s) {
		if (mm_verbose >= 1) fprintf(stderr, "[ERROR] the read group line is not started with @RG\n");
		goto err_set_rg;
	}
	if (strstr(s, "\t") != NULL) {
		if (mm_verbose >= 1) fprintf(stderr, "[ERROR] the read group line contained literal <tab> characters -- replace with escaped tabs: \\t\n");
		goto err_set_rg;
	}
	rg_line = (char*)malloc(strlen(s) + 1);
	strcpy(rg_line, s);
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
	return 0;

err_set_rg:
	free(rg_line);
	return -1;
}

int mm_write_sam_hdr(const mm_idx_t *idx, const char *rg, const char *ver, int argc, char *argv[])
{
	kstring_t str = {0,0,0};
	int ret = 0;
	mm_sprintf_lite(&str, "@HD\tVN:1.6\tSO:unsorted\tGO:query\n");
	if (idx) {
		uint32_t i;
		for (i = 0; i < idx->n_seq; ++i)
			mm_sprintf_lite(&str, "@SQ\tSN:%s\tLN:%d\n", idx->seq[i].name, idx->seq[i].len);
	}
	if (rg) ret = sam_write_rg_line(&str, rg);
	mm_sprintf_lite(&str, "@PG\tID:minimap2\tPN:minimap2");
	if (ver) mm_sprintf_lite(&str, "\tVN:%s", ver);
	if (argc > 1) {
		int i;
		mm_sprintf_lite(&str, "\tCL:minimap2");
		for (i = 1; i < argc; ++i)
			mm_sprintf_lite(&str, " %s", argv[i]);
	}
	mm_err_puts(str.s);
	free(str.s);
	return ret;
}

static void write_indel_ds(kstring_t *str, int64_t len, const uint8_t *seq, int64_t ll, int64_t lr) // write an indel to ds; adapted from minigraph
{
	int64_t i;
	if (ll + lr >= len) {
		mm_sprintf_lite(str, "[");
		for (i = 0; i < len; ++i)
			mm_sprintf_lite(str, "%c", "acgtn"[seq[i]]);
		mm_sprintf_lite(str, "]");
	} else {
		int64_t k = 0;
		if (ll > 0) {
			mm_sprintf_lite(str, "[");
			for (i = 0; i < ll; ++i)
				mm_sprintf_lite(str, "%c", "acgtn"[seq[k+i]]);
			mm_sprintf_lite(str, "]");
			k += ll;
		}
		for (i = 0; i < len - lr - ll; ++i)
			mm_sprintf_lite(str, "%c", "acgtn"[seq[k+i]]);
		k += len - lr - ll;
		if (lr > 0) {
			mm_sprintf_lite(str, "[");
			for (i = 0; i < lr; ++i)
				mm_sprintf_lite(str, "%c", "acgtn"[seq[k+i]]);
			mm_sprintf_lite(str, "]");
		}
	}
}

static void write_cs_ds_core(kstring_t *s, const uint8_t *tseq, const uint8_t *qseq, const mm_reg1_t *r, char *tmp, int no_iden, int is_ds, int write_tag)
{
	int i, q_off, t_off, q_len = 0, t_len = 0;
	if (write_tag) mm_sprintf_lite(s, "\t%cs:Z:", is_ds? 'd' : 'c');
	for (i = 0; i < (int)r->p->n_cigar; ++i) {
		int op = r->p->cigar[i]&0xf, len = r->p->cigar[i]>>4;
		if (op == MM_CIGAR_MATCH || op == MM_CIGAR_EQ_MATCH || op == MM_CIGAR_X_MISMATCH)
			q_len += len, t_len += len;
		else if (op == MM_CIGAR_INS)
			q_len += len;
		else if (op == MM_CIGAR_DEL || op == MM_CIGAR_N_SKIP)
			t_len += len;
	}
	for (i = q_off = t_off = 0; i < (int)r->p->n_cigar; ++i) {
		int j, op = r->p->cigar[i]&0xf, len = r->p->cigar[i]>>4;
		assert((op >= MM_CIGAR_MATCH && op <= MM_CIGAR_N_SKIP) || op == MM_CIGAR_EQ_MATCH || op == MM_CIGAR_X_MISMATCH);
		if (op == MM_CIGAR_MATCH || op == MM_CIGAR_EQ_MATCH || op == MM_CIGAR_X_MISMATCH) {
			int l_tmp = 0;
			for (j = 0; j < len; ++j) {
				if (qseq[q_off + j] != tseq[t_off + j]) {
					if (l_tmp > 0) {
						if (!no_iden) {
							tmp[l_tmp] = 0;
							mm_sprintf_lite(s, "=%s", tmp);
						} else mm_sprintf_lite(s, ":%d", l_tmp);
						l_tmp = 0;
					}
					mm_sprintf_lite(s, "*%c%c", "acgtn"[tseq[t_off + j]], "acgtn"[qseq[q_off + j]]);
				} else tmp[l_tmp++] = "ACGTN"[qseq[q_off + j]];
			}
			if (l_tmp > 0) {
				if (!no_iden) {
					tmp[l_tmp] = 0;
					mm_sprintf_lite(s, "=%s", tmp);
				} else mm_sprintf_lite(s, ":%d", l_tmp);
			}
			q_off += len, t_off += len;
		} else if (op == MM_CIGAR_INS) {
			if (is_ds) {
				int z, ll, lr, y = q_off;
				for (z = 1; z <= len; ++z)
					if (y - z < 0 || qseq[y + len - z] != qseq[y - z])
						break;
				lr = z - 1;
				for (z = 0; z < len; ++z)
					if (y + len + z >= q_len || qseq[y + len + z] != qseq[y + z])
						break;
				ll = z;
				mm_sprintf_lite(s, "+");
				write_indel_ds(s, len, &qseq[y], ll, lr);
			} else {
				for (j = 0, tmp[len] = 0; j < len; ++j)
					tmp[j] = "acgtn"[qseq[q_off + j]];
				mm_sprintf_lite(s, "+%s", tmp);
			}
			q_off += len;
		} else if (op == MM_CIGAR_DEL) {
			if (is_ds) {
				int z, ll, lr, x = t_off;
				for (z = 1; z <= len; ++z)
					if (x - z < 0 || tseq[x + len - z] != tseq[x - z])
						break;
				lr = z - 1;
				for (z = 0; z < len; ++z)
					if (x + len + z >= t_len || tseq[x + z] != tseq[x + len + z])
						break;
				ll = z;
				mm_sprintf_lite(s, "-");
				write_indel_ds(s, len, &tseq[x], ll, lr);
			} else {
				for (j = 0, tmp[len] = 0; j < len; ++j)
					tmp[j] = "acgtn"[tseq[t_off + j]];
				mm_sprintf_lite(s, "-%s", tmp);
			}
			t_off += len;
		} else { // intron
			assert(len >= 2);
			mm_sprintf_lite(s, "~%c%c%d%c%c", "acgtn"[tseq[t_off]], "acgtn"[tseq[t_off+1]],
				len, "acgtn"[tseq[t_off+len-2]], "acgtn"[tseq[t_off+len-1]]);
			t_off += len;
		}
	}
	assert(t_off == r->re - r->rs && q_off == r->qe - r->qs);
}

static void write_MD_core(kstring_t *s, const uint8_t *tseq, const uint8_t *qseq, const mm_reg1_t *r, char *tmp, int write_tag)
{
	int i, q_off, t_off, l_MD = 0;
	if (write_tag) mm_sprintf_lite(s, "\tMD:Z:");
	for (i = q_off = t_off = 0; i < (int)r->p->n_cigar; ++i) {
		int j, op = r->p->cigar[i]&0xf, len = r->p->cigar[i]>>4;
		assert((op >= MM_CIGAR_MATCH && op <= MM_CIGAR_N_SKIP) || op == MM_CIGAR_EQ_MATCH || op == MM_CIGAR_X_MISMATCH);
		if (op == MM_CIGAR_MATCH || op == MM_CIGAR_EQ_MATCH || op == MM_CIGAR_X_MISMATCH) {
			for (j = 0; j < len; ++j) {
				if (qseq[q_off + j] != tseq[t_off + j]) {
					mm_sprintf_lite(s, "%d%c", l_MD, "ACGTN"[tseq[t_off + j]]);
					l_MD = 0;
				} else ++l_MD;
			}
			q_off += len, t_off += len;
		} else if (op == MM_CIGAR_INS) {
			q_off += len;
		} else if (op == MM_CIGAR_DEL) {
			for (j = 0, tmp[len] = 0; j < len; ++j)
				tmp[j] = "ACGTN"[tseq[t_off + j]];
			mm_sprintf_lite(s, "%d^%s", l_MD, tmp);
			l_MD = 0;
			t_off += len;
		} else if (op == MM_CIGAR_N_SKIP) {
			t_off += len;
		}
	}
	if (l_MD > 0) mm_sprintf_lite(s, "%d", l_MD);
	assert(t_off == r->re - r->rs && q_off == r->qe - r->qs);
}

static void write_cs_ds_or_MD(void *km, kstring_t *s, const mm_idx_t *mi, const mm_bseq1_t *t, const mm_reg1_t *r, int no_iden, int is_MD, int is_ds, int write_tag, int is_qstrand)
{
	extern unsigned char seq_nt4_table[256];
	int i;
	uint8_t *qseq, *tseq;
	char *tmp;
	if (r->p == 0) return;
	qseq = (uint8_t*)kmalloc(km, r->qe - r->qs);
	tseq = (uint8_t*)kmalloc(km, r->re - r->rs);
	tmp = (char*)kmalloc(km, r->re - r->rs > r->qe - r->qs? r->re - r->rs + 1 : r->qe - r->qs + 1);
	if (is_qstrand) {
		mm_idx_getseq2(mi, r->rev, r->rid, r->rs, r->re, tseq);
		for (i = r->qs; i < r->qe; ++i)
			qseq[i - r->qs] = seq_nt4_table[(uint8_t)t->seq[i]];
	} else {
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
	}
	if (is_MD) write_MD_core(s, tseq, qseq, r, tmp, write_tag);
	else write_cs_ds_core(s, tseq, qseq, r, tmp, no_iden, is_ds, write_tag);
	kfree(km, qseq); kfree(km, tseq); kfree(km, tmp);
}

int mm_gen_cs_or_MD(void *km, char **buf, int *max_len, const mm_idx_t *mi, const mm_reg1_t *r, const char *seq, int is_MD, int no_iden, int is_qstrand)
{
	mm_bseq1_t t;
	kstring_t str;
	str.s = *buf, str.l = 0, str.m = *max_len;
	t.l_seq = strlen(seq);
	t.seq = (char*)seq;
	write_cs_ds_or_MD(km, &str, mi, &t, r, no_iden, is_MD, 0, 0, is_qstrand);
	*max_len = str.m;
	*buf = str.s;
	return str.l;
}

int mm_gen_cs(void *km, char **buf, int *max_len, const mm_idx_t *mi, const mm_reg1_t *r, const char *seq, int no_iden)
{
	return mm_gen_cs_or_MD(km, buf, max_len, mi, r, seq, 0, no_iden, 0);
}

int mm_gen_MD(void *km, char **buf, int *max_len, const mm_idx_t *mi, const mm_reg1_t *r, const char *seq)
{
	return mm_gen_cs_or_MD(km, buf, max_len, mi, r, seq, 1, 0, 0);
}

static inline void write_tags(kstring_t *s, const mm_reg1_t *r)
{
	int type;
	if (r->id == r->parent) type = r->inv? 'I' : 'P';
	else type = r->inv? 'i' : 'S';
	if (r->p) {
		mm_sprintf_lite(s, "\tNM:i:%d\tms:i:%d\tAS:i:%d\tnn:i:%d", r->blen - r->mlen + r->p->n_ambi, r->p->dp_max0, r->p->dp_score, r->p->n_ambi);
		if (r->p->trans_strand == 1 || r->p->trans_strand == 2)
			mm_sprintf_lite(s, "\tts:A:%c", "?+-?"[r->p->trans_strand]);
	}
	mm_sprintf_lite(s, "\ttp:A:%c\tcm:i:%d\ts1:i:%d", type, r->cnt, r->score);
	if (r->parent == r->id) mm_sprintf_lite(s, "\ts2:i:%d", r->subsc);
	if (r->p) {
		char buf[16];
		double div;
		div = 1.0 - mm_event_identity(r);
		if (div == 0.0) buf[0] = '0', buf[1] = 0;
		else snprintf(buf, 16, "%.4f", 1.0 - mm_event_identity(r));
		mm_sprintf_lite(s, "\tde:f:%s", buf);
	} else if (r->div >= 0.0f && r->div <= 1.0f) {
		char buf[16];
		if (r->div == 0.0f) buf[0] = '0', buf[1] = 0;
		else snprintf(buf, 16, "%.4f", r->div);
		mm_sprintf_lite(s, "\tdv:f:%s", buf);
	}
	if (r->split) mm_sprintf_lite(s, "\tzd:i:%d", r->split);
}

void mm_write_paf3(kstring_t *s, const mm_idx_t *mi, const mm_bseq1_t *t, const mm_reg1_t *r, void *km, int64_t opt_flag, int rep_len)
{
	s->l = 0;
	if (r == 0) {
		mm_sprintf_lite(s, "%s\t%d\t0\t0\t*\t*\t0\t0\t0\t0\t0\t0", t->name, t->l_seq);
		if (rep_len >= 0) mm_sprintf_lite(s, "\trl:i:%d", rep_len);
		return;
	}
	mm_sprintf_lite(s, "%s\t%d\t%d\t%d\t%c\t", t->name, t->l_seq, r->qs, r->qe, "+-"[r->rev]);
	if (mi->seq[r->rid].name) mm_sprintf_lite(s, "%s", mi->seq[r->rid].name);
	else mm_sprintf_lite(s, "%d", r->rid);
	mm_sprintf_lite(s, "\t%d", mi->seq[r->rid].len);
	if ((opt_flag & MM_F_QSTRAND) && r->rev)
		mm_sprintf_lite(s, "\t%d\t%d", mi->seq[r->rid].len - r->re, mi->seq[r->rid].len - r->rs);
	else
		mm_sprintf_lite(s, "\t%d\t%d", r->rs, r->re);
	mm_sprintf_lite(s, "\t%d\t%d", r->mlen, r->blen);
	mm_sprintf_lite(s, "\t%d", r->mapq);
	write_tags(s, r);
	if (rep_len >= 0) mm_sprintf_lite(s, "\trl:i:%d", rep_len);
	if (r->p && (opt_flag & MM_F_OUT_CG)) {
		uint32_t k;
		mm_sprintf_lite(s, "\tcg:Z:");
		for (k = 0; k < r->p->n_cigar; ++k)
			mm_sprintf_lite(s, "%d%c", r->p->cigar[k]>>4, MM_CIGAR_STR[r->p->cigar[k]&0xf]);
	}
	if (r->p && (opt_flag & (MM_F_OUT_CS|MM_F_OUT_DS|MM_F_OUT_MD)))
		write_cs_ds_or_MD(km, s, mi, t, r, !(opt_flag&MM_F_OUT_CS_LONG), !!(opt_flag&MM_F_OUT_MD), !!(opt_flag&MM_F_OUT_DS), 1, !!(opt_flag&MM_F_QSTRAND));
	if ((opt_flag & MM_F_COPY_COMMENT) && t->comment)
		mm_sprintf_lite(s, "\t%s", t->comment);
}

void mm_write_paf(kstring_t *s, const mm_idx_t *mi, const mm_bseq1_t *t, const mm_reg1_t *r, void *km, int64_t opt_flag)
{
	mm_write_paf3(s, mi, t, r, km, opt_flag, -1);
}

static void sam_write_sq(kstring_t *s, char *seq, int l, int rev, int comp)
{
	extern unsigned char seq_comp_table[256];
	if (rev) {
		int i;
		str_enlarge(s, l);
		for (i = 0; i < l; ++i) {
			int c = seq[l - 1 - i];
			s->s[s->l + i] = c < 128 && comp? seq_comp_table[c] : c;
		}
		s->l += l;
	} else str_copy(s, seq, seq + l);
}

static inline const mm_reg1_t *get_sam_pri(int n_regs, const mm_reg1_t *regs)
{
	int i;
	for (i = 0; i < n_regs; ++i)
		if (regs[i].sam_pri)
			return &regs[i];
	assert(n_regs == 0);
	return NULL;
}

static void write_sam_cigar(kstring_t *s, int sam_flag, int in_tag, int qlen, const mm_reg1_t *r, int64_t opt_flag)
{
	if (r->p == 0) {
		mm_sprintf_lite(s, "*");
	} else {
		uint32_t k, clip_len[2];
		clip_len[0] = r->rev? qlen - r->qe : r->qs;
		clip_len[1] = r->rev? r->qs : qlen - r->qe;
		if (in_tag) {
			int clip_char = (((sam_flag&0x800) || ((sam_flag&0x100) && (opt_flag&MM_F_SECONDARY_SEQ))) &&
							 !(opt_flag&MM_F_SOFTCLIP)) ? 5 : 4;
			mm_sprintf_lite(s, "\tCG:B:I");
			if (clip_len[0]) mm_sprintf_lite(s, ",%u", clip_len[0]<<4|clip_char);
			for (k = 0; k < r->p->n_cigar; ++k)
				mm_sprintf_lite(s, ",%u", r->p->cigar[k]);
			if (clip_len[1]) mm_sprintf_lite(s, ",%u", clip_len[1]<<4|clip_char);
		} else {
			int clip_char = (((sam_flag&0x800) || ((sam_flag&0x100) && (opt_flag&MM_F_SECONDARY_SEQ))) &&
							 !(opt_flag&MM_F_SOFTCLIP)) ? 'H' : 'S';
			assert(clip_len[0] < qlen && clip_len[1] < qlen);
			if (clip_len[0]) mm_sprintf_lite(s, "%d%c", clip_len[0], clip_char);
			for (k = 0; k < r->p->n_cigar; ++k)
				mm_sprintf_lite(s, "%d%c", r->p->cigar[k]>>4, MM_CIGAR_STR[r->p->cigar[k]&0xf]);
			if (clip_len[1]) mm_sprintf_lite(s, "%d%c", clip_len[1], clip_char);
		}
	}
}

void mm_write_sam3(kstring_t *s, const mm_idx_t *mi, const mm_bseq1_t *t, int seg_idx, int reg_idx, int n_seg, const int *n_regss, const mm_reg1_t *const* regss, void *km, int64_t opt_flag, int rep_len)
{
	const int max_bam_cigar_op = 65535;
	int flag, n_regs = n_regss[seg_idx], cigar_in_tag = 0;
	int this_rid = -1, this_pos = -1;
	const mm_reg1_t *regs = regss[seg_idx], *r_prev = NULL, *r_next;
	const mm_reg1_t *r = n_regs > 0 && reg_idx < n_regs && reg_idx >= 0? &regs[reg_idx] : NULL;

	// find the primary of the previous and the next segments, if they are mapped
	if (n_seg > 1) {
		int i, next_sid = (seg_idx + 1) % n_seg;
		r_next = get_sam_pri(n_regss[next_sid], regss[next_sid]);
		if (n_seg > 2) {
			for (i = 1; i <= n_seg - 1; ++i) {
				int prev_sid = (seg_idx + n_seg - i) % n_seg;
				if (n_regss[prev_sid] > 0) {
					r_prev = get_sam_pri(n_regss[prev_sid], regss[prev_sid]);
					break;
				}
			}
		} else r_prev = r_next;
	} else r_prev = r_next = NULL;

	// write QNAME
	s->l = 0;
	mm_sprintf_lite(s, "%s", t->name);
	if (n_seg > 1) s->l = mm_qname_len(t->name); // trim the suffix like /1 or /2

	// write flag
	flag = n_seg > 1? 0x1 : 0x0;
	if (r == 0) {
		flag |= 0x4;
	} else {
		if (r->rev) flag |= 0x10;
		if (r->parent != r->id) flag |= 0x100;
		else if (!r->sam_pri) flag |= 0x800;
	}
	if (n_seg > 1) {
		if (r && r->proper_frag) flag |= 0x2; // TODO: this doesn't work when there are more than 2 segments
		if (seg_idx == 0) flag |= 0x40;
		else if (seg_idx == n_seg - 1) flag |= 0x80;
		if (r_next == NULL) flag |= 0x8;
		else if (r_next->rev) flag |= 0x20;
	}
	mm_sprintf_lite(s, "\t%d", flag);

	// write coordinate, MAPQ and CIGAR
	if (r == 0) {
		if (r_prev) {
			this_rid = r_prev->rid, this_pos = r_prev->rs;
			mm_sprintf_lite(s, "\t%s\t%d\t0\t*", mi->seq[this_rid].name, this_pos+1);
		} else mm_sprintf_lite(s, "\t*\t0\t0\t*");
	} else {
		this_rid = r->rid, this_pos = r->rs;
		mm_sprintf_lite(s, "\t%s\t%d\t%d\t", mi->seq[r->rid].name, r->rs+1, r->mapq);
		if ((opt_flag & MM_F_LONG_CIGAR) && r->p && r->p->n_cigar > max_bam_cigar_op - 2) {
			int n_cigar = r->p->n_cigar;
			if (r->qs != 0) ++n_cigar;
			if (r->qe != t->l_seq) ++n_cigar;
			if (n_cigar > max_bam_cigar_op)
				cigar_in_tag = 1;
		}
		if (cigar_in_tag) {
			int slen;
			if ((flag & 0x900) == 0 || (opt_flag & MM_F_SOFTCLIP)) slen = t->l_seq;
			else if ((flag & 0x100) && !(opt_flag & MM_F_SECONDARY_SEQ)) slen = 0;
			else slen = r->qe - r->qs;
			mm_sprintf_lite(s, "%dS%dN", slen, r->re - r->rs);
		} else write_sam_cigar(s, flag, 0, t->l_seq, r, opt_flag);
	}

	// write mate positions
	if (n_seg > 1) {
		int tlen = 0;
		if (this_rid >= 0 && r_next) {
			if (this_rid == r_next->rid) {
				if (r) {
					int this_pos5 = r->rev? r->re - 1 : this_pos;
					int next_pos5 = r_next->rev? r_next->re - 1 : r_next->rs;
					tlen = next_pos5 - this_pos5;
				}
				mm_sprintf_lite(s, "\t=\t");
			} else mm_sprintf_lite(s, "\t%s\t", mi->seq[r_next->rid].name);
			mm_sprintf_lite(s, "%d\t", r_next->rs + 1);
		} else if (r_next) { // && this_rid < 0
			mm_sprintf_lite(s, "\t%s\t%d\t", mi->seq[r_next->rid].name, r_next->rs + 1);
		} else if (this_rid >= 0) { // && r_next == NULL
			mm_sprintf_lite(s, "\t=\t%d\t", this_pos + 1); // next segment will take r's coordinate
		} else mm_sprintf_lite(s, "\t*\t0\t"); // neither has coordinates
		if (tlen > 0) ++tlen;
		else if (tlen < 0) --tlen;
		mm_sprintf_lite(s, "%d\t", tlen);
	} else mm_sprintf_lite(s, "\t*\t0\t0\t");

	// write SEQ and QUAL
	if (r == 0) {
		sam_write_sq(s, t->seq, t->l_seq, 0, 0);
		mm_sprintf_lite(s, "\t");
		if (t->qual) sam_write_sq(s, t->qual, t->l_seq, 0, 0);
		else mm_sprintf_lite(s, "*");
	} else {
		if ((flag & 0x900) == 0 || (opt_flag & MM_F_SOFTCLIP)) {
			sam_write_sq(s, t->seq, t->l_seq, r->rev, r->rev);
			mm_sprintf_lite(s, "\t");
			if (t->qual) sam_write_sq(s, t->qual, t->l_seq, r->rev, 0);
			else mm_sprintf_lite(s, "*");
		} else if ((flag & 0x100) && !(opt_flag & MM_F_SECONDARY_SEQ)){
			mm_sprintf_lite(s, "*\t*");
		} else {
			sam_write_sq(s, t->seq + r->qs, r->qe - r->qs, r->rev, r->rev);
			mm_sprintf_lite(s, "\t");
			if (t->qual) sam_write_sq(s, t->qual + r->qs, r->qe - r->qs, r->rev, 0);
			else mm_sprintf_lite(s, "*");
		}
	}

	// write tags
	if (mm_rg_id[0]) mm_sprintf_lite(s, "\tRG:Z:%s", mm_rg_id);
	if (n_seg > 2) mm_sprintf_lite(s, "\tFI:i:%d", seg_idx);
	if (r) {
		write_tags(s, r);
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
					mm_sprintf_lite(s, ",%d,%d;", q->mapq, q->blen - q->mlen + q->p->n_ambi);
				}
			}
		}
		if (r->p && (opt_flag & (MM_F_OUT_CS|MM_F_OUT_DS|MM_F_OUT_MD)))
			write_cs_ds_or_MD(km, s, mi, t, r, !(opt_flag&MM_F_OUT_CS_LONG), opt_flag&MM_F_OUT_MD, !!(opt_flag&MM_F_OUT_DS), 1, 0);
		if (cigar_in_tag)
			write_sam_cigar(s, flag, 1, t->l_seq, r, opt_flag);
	}
	if (rep_len >= 0) mm_sprintf_lite(s, "\trl:i:%d", rep_len);

	if ((opt_flag & MM_F_COPY_COMMENT) && t->comment)
		mm_sprintf_lite(s, "\t%s", t->comment);

	s->s[s->l] = 0; // we always have room for an extra byte (see str_enlarge)
}

void mm_write_sam2(kstring_t *s, const mm_idx_t *mi, const mm_bseq1_t *t, int seg_idx, int reg_idx, int n_seg, const int *n_regss, const mm_reg1_t *const* regss, void *km, int64_t opt_flag)
{
	mm_write_sam3(s, mi, t, seg_idx, reg_idx, n_seg, n_regss, regss, km, opt_flag, -1);
}

void mm_write_sam(kstring_t *s, const mm_idx_t *mi, const mm_bseq1_t *t, const mm_reg1_t *r, int n_regs, const mm_reg1_t *regs)
{
	int i;
	for (i = 0; i < n_regs; ++i)
		if (r == &regs[i]) break;
	mm_write_sam2(s, mi, t, 0, i, 1, &n_regs, &regs, NULL, 0);
}
