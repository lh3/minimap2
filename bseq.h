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
	char *name, *seq, *qual;
} mm_bseq1_t;

mm_bseq_file_t *mm_bseq_open(const char *fn);
void mm_bseq_close(mm_bseq_file_t *fp);
mm_bseq1_t *mm_bseq_read2(mm_bseq_file_t *fp, int chunk_size, int with_qual, int check_name, int *n_);
mm_bseq1_t *mm_bseq_read(mm_bseq_file_t *fp, int chunk_size, int with_qual, int *n_);
mm_bseq1_t *mm_bseq_read_multi(int n_fp, mm_bseq_file_t **fp, int chunk_size, int with_qual, int *n_);
int mm_bseq_eof(mm_bseq_file_t *fp);

extern unsigned char seq_nt4_table[256];

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

#ifdef __cplusplus
}
#endif

#endif
