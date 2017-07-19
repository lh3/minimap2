#ifndef MM_BSEQ_H
#define MM_BSEQ_H

#include <stdint.h>

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
mm_bseq1_t *mm_bseq_read(mm_bseq_file_t *fp, int chunk_size, int with_qual, int *n_);
int mm_bseq_eof(mm_bseq_file_t *fp);

extern unsigned char seq_nt4_table[256];

#ifdef __cplusplus
}
#endif

#endif
