#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include "mmpriv.h"

FILE *mm_split_init(const char *prefix, const mm_idx_t *mi)
{
	char *fn;
	FILE *fp;
	uint32_t i;
	fn = (char*)calloc(strlen(prefix) + 10, 1);
	sprintf(fn, "%s.%.4d.tmp", prefix, mi->index);
	fp = fopen(fn, "wb");
	assert(fp);
	mm_err_fwrite(&mi->n_seq, 4, 1, fp);
	for (i = 0; i < mi->n_seq; ++i) {
		uint8_t l;
		l = strlen(mi->seq[i].name);
		mm_err_fwrite(&l, 1, 1, fp);
		mm_err_fwrite(mi->seq[i].name, 1, l, fp);
		mm_err_fwrite(&mi->seq[i].len, 4, 1, fp);
	}
	free(fn);
	return fp;
}
