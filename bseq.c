#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#define __STDC_LIMIT_MACROS
#include "bseq.h"
#include "kvec.h"
#include "kseq.h"
KSEQ_INIT2(, gzFile, gzread)

unsigned char seq_comp_table[256] = {
	  0,   1,	2,	 3,	  4,   5,	6,	 7,	  8,   9,  10,	11,	 12,  13,  14,	15,
	 16,  17,  18,	19,	 20,  21,  22,	23,	 24,  25,  26,	27,	 28,  29,  30,	31,
	 32,  33,  34,	35,	 36,  37,  38,	39,	 40,  41,  42,	43,	 44,  45,  46,	47,
	 48,  49,  50,	51,	 52,  53,  54,	55,	 56,  57,  58,	59,	 60,  61,  62,	63,
	 64, 'T', 'V', 'G', 'H', 'E', 'F', 'C', 'D', 'I', 'J', 'M', 'L', 'K', 'N', 'O',
	'P', 'Q', 'Y', 'S', 'A', 'A', 'B', 'W', 'X', 'R', 'Z',	91,	 92,  93,  94,	95,
	 96, 't', 'v', 'g', 'h', 'e', 'f', 'c', 'd', 'i', 'j', 'm', 'l', 'k', 'n', 'o',
	'p', 'q', 'y', 's', 'a', 'a', 'b', 'w', 'x', 'r', 'z', 123, 124, 125, 126, 127,
	128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143,
	144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159,
	160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175,
	176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191,
	192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207,
	208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223,
	224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239,
	240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255
};

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
	f = fn && strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(0, "r");
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

static inline char *kstrdup(const kstring_t *s)
{
	char *t;
	t = (char*)malloc(s->l + 1);
	memcpy(t, s->s, s->l + 1);
	return t;
}

static inline void kseq2bseq(kseq_t *ks, mm_bseq1_t *s, int with_qual, int with_comment)
{
	int i;
	if (ks->name.l == 0)
		fprintf(stderr, "[WARNING]\033[1;31m empty sequence name in the input.\033[0m\n");
	s->name = kstrdup(&ks->name);
	s->seq = kstrdup(&ks->seq);
	for (i = 0; i < (int)ks->seq.l; ++i) // convert U to T
		if (s->seq[i] == 'u' || s->seq[i] == 'U')
			--s->seq[i];
	s->qual = with_qual && ks->qual.l? kstrdup(&ks->qual) : 0;
	s->comment = with_comment && ks->comment.l? kstrdup(&ks->comment) : 0;
	s->l_seq = ks->seq.l;
}

mm_bseq1_t *mm_bseq_read3(mm_bseq_file_t *fp, int64_t chunk_size, int with_qual, int with_comment, int frag_mode, int *n_)
{
	int64_t size = 0;
	int ret;
	kvec_t(mm_bseq1_t) a = {0,0,0};
	kseq_t *ks = fp->ks;
	*n_ = 0;
	if (fp->s.seq) {
		kv_resize(mm_bseq1_t, 0, a, 256);
		kv_push(mm_bseq1_t, 0, a, fp->s);
		size = fp->s.l_seq;
		memset(&fp->s, 0, sizeof(mm_bseq1_t));
	}
	while ((ret = kseq_read(ks)) >= 0) {
		mm_bseq1_t *s;
		assert(ks->seq.l <= INT32_MAX);
		if (a.m == 0) kv_resize(mm_bseq1_t, 0, a, 256);
		kv_pushp(mm_bseq1_t, 0, a, &s);
		kseq2bseq(ks, s, with_qual, with_comment);
		size += s->l_seq;
		if (size >= chunk_size) {
			if (frag_mode && a.a[a.n-1].l_seq < CHECK_PAIR_THRES) {
				while ((ret = kseq_read(ks)) >= 0) {
					kseq2bseq(ks, &fp->s, with_qual, with_comment);
					if (mm_qname_same(fp->s.name, a.a[a.n-1].name)) {
						kv_push(mm_bseq1_t, 0, a, fp->s);
						memset(&fp->s, 0, sizeof(mm_bseq1_t));
					} else break;
				}
			}
			break;
		}
	}
	if (ret < -1) {
		if (a.n) fprintf(stderr, "[WARNING]\033[1;31m failed to parse the FASTA/FASTQ record next to '%s'. Continue anyway.\033[0m\n", a.a[a.n-1].name);
		else fprintf(stderr, "[WARNING]\033[1;31m failed to parse the first FASTA/FASTQ record. Continue anyway.\033[0m\n");
	}
	*n_ = a.n;
	return a.a;
}

mm_bseq1_t *mm_bseq_read2(mm_bseq_file_t *fp, int64_t chunk_size, int with_qual, int frag_mode, int *n_)
{
	return mm_bseq_read3(fp, chunk_size, with_qual, 0, frag_mode, n_);
}

mm_bseq1_t *mm_bseq_read(mm_bseq_file_t *fp, int64_t chunk_size, int with_qual, int *n_)
{
	return mm_bseq_read2(fp, chunk_size, with_qual, 0, n_);
}

mm_bseq1_t *mm_bseq_read_frag2(int n_fp, mm_bseq_file_t **fp, int64_t chunk_size, int with_qual, int with_comment, int *n_)
{
	int i;
	int64_t size = 0;
	kvec_t(mm_bseq1_t) a = {0,0,0};
	*n_ = 0;
	if (n_fp < 1) return 0;
	while (1) {
		int n_read = 0;
		for (i = 0; i < n_fp; ++i)
			if (kseq_read(fp[i]->ks) >= 0)
				++n_read;
		if (n_read < n_fp) {
			if (n_read > 0)
				fprintf(stderr, "[W::%s]\033[1;31m query files have different number of records; extra records skipped.\033[0m\n", __func__);
			break; // some file reaches the end
		}
		if (a.m == 0) kv_resize(mm_bseq1_t, 0, a, 256);
		for (i = 0; i < n_fp; ++i) {
			mm_bseq1_t *s;
			kv_pushp(mm_bseq1_t, 0, a, &s);
			kseq2bseq(fp[i]->ks, s, with_qual, with_comment);
			size += s->l_seq;
		}
		if (size >= chunk_size) break;
	}
	*n_ = a.n;
	return a.a;
}

mm_bseq1_t *mm_bseq_read_frag(int n_fp, mm_bseq_file_t **fp, int64_t chunk_size, int with_qual, int *n_)
{
	return mm_bseq_read_frag2(n_fp, fp, chunk_size, with_qual, 0, n_);
}

int mm_bseq_eof(mm_bseq_file_t *fp)
{
	return (ks_eof(fp->ks->f) && fp->s.seq == 0);
}
