#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "bseq.h"
#include "kvec.h"
#include "kseq.h"
KSEQ_INIT2(, gzFile, gzread)

#define CHECK_PAIR_THRES 1000000

struct mm_bseq_file_s {
	gzFile fp;
	kseq_t *ks;
	mm_bseq1_t s;
};

mm_bseq_file_t *mm_bseq_open(const char *fn)
{
	mm_bseq_file_t *fp;
	gzFile f;
	f = fn && strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
	if (f == 0) return 0;
	fp = (mm_bseq_file_t*)calloc(1, sizeof(mm_bseq_file_t));
	fp->fp = f;
	fp->ks = kseq_init(fp->fp);
	return fp;
}

void mm_bseq_close(mm_bseq_file_t *fp)
{
	kseq_destroy(fp->ks);
	gzclose(fp->fp);
	free(fp);
}

static inline int mm_qname_same(const char *s1, const char *s2)
{
	int l1, l2;
	l1 = strlen(s1);
	l2 = strlen(s2);
	if (l1 != l2 || l1 < 3) return 0;
	if (!(s1[l1-1] >= '1' && s1[l1-1] <= '2' && s1[l1-2] == '/')) return 0;
	if (!(s2[l2-1] >= '1' && s2[l2-1] <= '2' && s2[l2-2] == '/')) return 0;
	return (strncmp(s1, s2, l1 - 2) == 0);
}

static inline void kseq2bseq(kseq_t *ks, mm_bseq1_t *s, int with_qual)
{
	s->name = strdup(ks->name.s);
	s->seq = strdup(ks->seq.s);
	s->qual = with_qual && ks->qual.l? strdup(ks->qual.s) : 0;
	s->l_seq = ks->seq.l;
}

mm_bseq1_t *mm_bseq_read2(mm_bseq_file_t *fp, int chunk_size, int with_qual, int check_name, int *n_)
{
	int64_t size = 0;
	kvec_t(mm_bseq1_t) a = {0,0,0};
	kseq_t *ks = fp->ks;
	*n_ = 0;
	if (fp->s.seq) {
		kv_resize(mm_bseq1_t, 0, a, 256);
		kv_push(mm_bseq1_t, 0, a, fp->s);
		size = fp->s.l_seq;
		memset(&fp->s, 0, sizeof(mm_bseq1_t));
	}
	while (kseq_read(ks) >= 0) {
		mm_bseq1_t *s;
		assert(ks->seq.l <= INT32_MAX);
		if (a.m == 0) kv_resize(mm_bseq1_t, 0, a, 256);
		kv_pushp(mm_bseq1_t, 0, a, &s);
		kseq2bseq(ks, s, with_qual);
		size += s->l_seq;
		if (size >= chunk_size) {
			if (check_name && a.a[a.n-1].l_seq < CHECK_PAIR_THRES) {
				while (kseq_read(ks) >= 0) {
					kseq2bseq(ks, &fp->s, with_qual);
					if (mm_qname_same(fp->s.name, a.a[a.n-1].name)) {
						kv_push(mm_bseq1_t, 0, a, fp->s);
						memset(&fp->s, 0, sizeof(mm_bseq1_t));
					} else break;
				}
			}
			break;
		}
	}
	*n_ = a.n;
	return a.a;
}

mm_bseq1_t *mm_bseq_read(mm_bseq_file_t *fp, int chunk_size, int with_qual, int *n_)
{
	return mm_bseq_read2(fp, chunk_size, with_qual, 0, n_);
}

mm_bseq1_t *mm_bseq_read_multi(int n_fp, mm_bseq_file_t **fp, int chunk_size, int with_qual, int *n_)
{
	int i;
	int64_t size = 0;
	kvec_t(mm_bseq1_t) a = {0,0,0};
	*n_ = 0;
	if (n_fp < 1) return 0;
	while (1) {
		for (i = 0; i < n_fp; ++i)
			if (kseq_read(fp[i]->ks) < 0)
				break;
		if (i != n_fp) break; // some file reaches the end
		if (a.m == 0) kv_resize(mm_bseq1_t, 0, a, 256);
		for (i = 0; i < n_fp; ++i) {
			mm_bseq1_t *s;
			kv_pushp(mm_bseq1_t, 0, a, &s);
			kseq2bseq(fp[i]->ks, s, with_qual);
			size += s->l_seq;
		}
		if (size >= chunk_size) break;
	}
	*n_ = a.n;
	return a.a;
}

int mm_bseq_eof(mm_bseq_file_t *fp)
{
	return (ks_eof(fp->ks->f) && fp->s.seq == 0);
}
