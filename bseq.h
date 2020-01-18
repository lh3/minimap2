#ifndef MM_BSEQ_H
#define MM_BSEQ_H

#include <stdint.h>
#include <string.h>

#ifdef __cplusplus
extern "C" {
#endif

struct mm_bseq_file_s;
typedef struct mm_bseq_file_s mm_bseq_file_t;

typedef struct {
	int l_seq, rid;
	char *name, *seq, *qual, *comment;
} mm_bseq1_t;

mm_bseq_file_t *mm_bseq_open(const char *fn);
void mm_bseq_close(mm_bseq_file_t *fp);
mm_bseq1_t *mm_bseq_read3(mm_bseq_file_t *fp, int64_t chunk_size, int with_qual, int with_comment, int frag_mode, int *n_);
mm_bseq1_t *mm_bseq_read2(mm_bseq_file_t *fp, int64_t chunk_size, int with_qual, int frag_mode, int *n_);
mm_bseq1_t *mm_bseq_read(mm_bseq_file_t *fp, int64_t chunk_size, int with_qual, int *n_);
mm_bseq1_t *mm_bseq_read_frag2(int n_fp, mm_bseq_file_t **fp, int64_t chunk_size, int with_qual, int with_comment, int *n_);
mm_bseq1_t *mm_bseq_read_frag(int n_fp, mm_bseq_file_t **fp, int64_t chunk_size, int with_qual, int *n_);
int mm_bseq_eof(mm_bseq_file_t *fp);

extern unsigned char seq_nt4_table[256];
extern unsigned char seq_comp_table[256];

static inline int mm_qname_len(const char *s)
{
	int l;
	l = strlen(s);
	return l >= 3 && s[l-1] >= '0' && s[l-1] <= '9' && s[l-2] == '/'? l - 2 : l;
}

static inline int mm_qname_same(const char *s1, const char *s2)
{
	int l1, l2;
	l1 = mm_qname_len(s1);
	l2 = mm_qname_len(s2);
	return (l1 == l2 && strncmp(s1, s2, l1) == 0);
}

static inline void mm_revcomp_bseq(mm_bseq1_t *s)
{
	int i, t, l = s->l_seq;
	for (i = 0; i < l>>1; ++i) {
		t = s->seq[l - i - 1];
		s->seq[l - i - 1] = seq_comp_table[(uint8_t)s->seq[i]];
		s->seq[i] = seq_comp_table[t];
	}
	if (l&1) s->seq[l>>1] = seq_comp_table[(uint8_t)s->seq[l>>1]];
	if (s->qual)
		for (i = 0; i < l>>1; ++i)
			t = s->qual[l - i - 1], s->qual[l - i - 1] = s->qual[i], s->qual[i] = t;
}

#ifdef __cplusplus
}
#endif

#endif
