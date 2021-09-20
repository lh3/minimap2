#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include "mmpriv.h"

FILE *mm_split_init(const char *fn, const mm_idx_t *mi)
{
	FILE *fp;
	uint32_t i, k = mi->k;
	if ((fp = fopen(fn, "wb")) == NULL) {
		if (mm_verbose >= 1)
			fprintf(stderr, "[ERROR]\033[1;31m failed to write to intermediate file '%s'\033[0m: %s\n", fn, strerror(errno));
		exit(1);
	}
	mm_err_fwrite(&k, 4, 1, fp);
	mm_err_fwrite(&mi->n_seq, 4, 1, fp);
	for (i = 0; i < mi->n_seq; ++i) {
		uint32_t l;
		l = strlen(mi->seq[i].name);
		mm_err_fwrite(&l, 1, 4, fp);
		mm_err_fwrite(mi->seq[i].name, 1, l, fp);
		mm_err_fwrite(&mi->seq[i].len, 4, 1, fp);
	}
	return fp;
}

FILE *mm_split_init_tmp(const char *prefix, const mm_idx_t *mi)
{
	FILE *fp;
	char *fn;
	fn = (char*)calloc(strlen(prefix) + 10, 1);
	sprintf(fn, "%s.%.4d.tmp", prefix, mi->index);
	fp = mm_split_init(fn, mi);
	free(fn);
	return fp;
}

mm_idx_t *mm_split_merge_prep(const char **intermediates, int n_splits, FILE **fp, uint32_t *n_seq_part)
{
	mm_idx_t *mi = 0;
	int i, j;

	if (n_splits < 1) return 0;
	for (i = 0; i < n_splits; ++i) {
		if ((fp[i] = fopen(intermediates[i], "rb")) == 0) {
			if (mm_verbose >= 1)
				fprintf(stderr, "ERROR: failed to open intermediate file '%s': %s\n", intermediates[i], strerror(errno));
			for (j = 0; j < i; ++j)
				fclose(fp[j]);
			return 0;
		}
	}

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
			uint32_t l;
			mm_err_fread(&l, 1, 4, fp[i]);
			mi->seq[j].name = (char*)calloc(l + 1, 1);
			mm_err_fread(mi->seq[j].name, 1, l, fp[i]);
			mm_err_fread(&mi->seq[j].len, 4, 1, fp[i]);
		}
	}
	return mi;
}

mm_idx_t *mm_split_merge_prep_tmp(const char *prefix, int n_splits, FILE **fp, uint32_t *n_seq_part)
{
	char **filenames;
	int i;
	mm_idx_t *mi;

	if (n_splits < 1) return 0;
	filenames = CALLOC(char*, n_splits);
	for (i = 0; i < n_splits; ++i) {
		filenames[i] = CALLOC(char, strlen(prefix) + 10);
		sprintf(filenames[i], "%s.%.4d.tmp", prefix, i);
	}

	mi = mm_split_merge_prep((const char**)filenames, n_splits, fp, n_seq_part);

	for (i = 0; i < n_splits; ++i) {
		free(filenames[i]);
	}
	free(filenames);

	return mi;
}

char** mm_split_tmp_intermediates(const char *prefix, int n_splits)
{
	char **filenames;
	int i;

	if (n_splits < 1) return 0;
	filenames = CALLOC(char*, n_splits);
	for (i = 0; i < n_splits; ++i) {
		filenames[i] = CALLOC(char, strlen(prefix) + 10);
		sprintf(filenames[i], "%s.%.4d.tmp", prefix, i);
	}

	return filenames;
}
