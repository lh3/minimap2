#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include "mmpriv.h"

FILE *mm_split_init(const char *prefix, const mm_idx_t *mi)
{
	char *fn;
	FILE *fp;
	uint32_t i, k = mi->k;
	fn = (char*)calloc(strlen(prefix) + 10, 1);
	sprintf(fn, "%s.%.4d.tmp", prefix, mi->index);
	fp = fopen(fn, "wb");
	assert(fp);
	mm_err_fwrite(&k, 4, 1, fp);
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

FILE **mm_split_merge_prep(const char *prefix, int n_splits, mm_idx_t **mi_)
{
	mm_idx_t *mi = 0;
	FILE **fp;
	char *fn;
	int i, j;
	uint32_t m_seq = 0;

	if (n_splits < 1) return 0;
	*mi_ = 0;
	fp = (FILE**)calloc(n_splits, sizeof(FILE*));
	fn = (char*)calloc(strlen(prefix) + 10, 1);
	for (i = 0; i < n_splits; ++i) {
		sprintf(fn, "%s.%.4d.tmp", prefix, i);
		if ((fp[i] = fopen(fn, "rb")) == 0) {
			if (mm_verbose >= 1)
				fprintf(stderr, "ERROR: failed to open temporary file '%s'\n", fn);
			for (j = 0; j < i; ++j)
				fclose(fp[j]);
			free(fp);
			return 0;
		}
	}
	free(fn);
	mi = (mm_idx_t*)calloc(1, sizeof(mm_idx_t));
	for (i = 0; i < n_splits; ++i) {
		uint32_t k, n;
		mm_err_fread(&k, 4, 1, fp[i]);
		mm_err_fread(&n, 4, 1, fp[i]);
		mi->k = k;
		if (mi->n_seq + n > m_seq) {
			m_seq = mi->n_seq + n;
			kroundup32(m_seq);
			mi->seq = (mm_idx_seq_t*)realloc(mi->seq, sizeof(mm_idx_seq_t) * m_seq);
		}
		for (k = 0; k < n; ++i) {
			uint8_t l;
			mm_err_fread(&l, 1, 1, fp[i]);
			mi->seq[mi->n_seq].name = (char*)calloc(l + 1, 1);
			mm_err_fread(mi->seq[mi->n_seq].name, 1, l, fp[i]);
			mm_err_fread(&mi->seq[mi->n_seq].len, 4, 1, fp[i]);
			++mi->n_seq;
		}
	}
	*mi_ = mi;
	return fp;
}
