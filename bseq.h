#ifndef MM_BSEQ_H
#define MM_BSEQ_H

#include <stdint.h>

struct bseq_file_s;
typedef struct bseq_file_s bseq_file_t;

typedef struct {
	int l_seq, rid;
	char *name, *seq;
} bseq1_t;

bseq_file_t *bseq_open(const char *fn);
void bseq_close(bseq_file_t *fp);
bseq1_t *bseq_read(bseq_file_t *fp, int chunk_size, int *n_);
int bseq_eof(bseq_file_t *fp);

#endif
