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

mm_idx_t *mm_split_merge_prep(const char *prefix, int n_splits, FILE **fp, uint32_t *n_seq_part)
{
	mm_idx_t *mi = 0;
	char *fn;
	int i, j;

	if (n_splits < 1) return 0;
	fn = CALLOC(char, strlen(prefix) + 10);
	for (i = 0; i < n_splits; ++i) {
		sprintf(fn, "%s.%.4d.tmp", prefix, i);
		if ((fp[i] = fopen(fn, "rb")) == 0) {
			if (mm_verbose >= 1)
				fprintf(stderr, "ERROR: failed to open temporary file '%s'\n", fn);
			for (j = 0; j < i; ++j)
				fclose(fp[j]);
			free(fn);
			return 0;
		}
	}
	free(fn);

	mi = CALLOC(mm_idx_t, 1);
	for (i = 0; i < n_splits; ++i) {
		mm_err_fread(&mi->k, 4, 1, fp[i]); // TODO: check if k is all the same
		mm_err_fread(&n_seq_part[i], 4, 1, fp[i]);
		mi->n_seq += n_seq_part[i];
	}
	mi->seq = CALLOC(mm_idx_seq_t, mi->n_seq);
	for (i = j = 0; i < n_splits; ++i) {
		uint32_t k;
		for (k = 0; k < n_seq_part[i]; ++k, ++j) {
			uint8_t l;
			mm_err_fread(&l, 1, 1, fp[i]);
			mi->seq[j].name = (char*)calloc(l + 1, 1);
			mm_err_fread(mi->seq[j].name, 1, l, fp[i]);
			mm_err_fread(&mi->seq[j].len, 4, 1, fp[i]);
		}
	}
	return mi;
}

void mm_split_rm_tmp(const char *prefix, int n_splits)
{
	int i;
	char *fn;
	fn = CALLOC(char, strlen(prefix) + 10);
	for (i = 0; i < n_splits; ++i) {
		sprintf(fn, "%s.%.4d.tmp", prefix, i);
		remove(fn);
	}
	free(fn);
}
