#ifndef MINIMAP_H
#define MINIMAP_H

#include <stdint.h>
#include <stdio.h>
#include <sys/types.h>
#include "bseq.h"

#define MM_IDX_DEF_B    14
#define MM_DEREP_Q50    5.0

#define MM_F_WITH_REP  0x1
#define MM_F_NO_SELF   0x2
#define MM_F_NO_ISO    0x4
#define MM_F_AVA       0x8

typedef struct {
 	uint64_t x, y;
} mm128_t;

typedef struct { size_t n, m; mm128_t *a; } mm128_v;
typedef struct { size_t n, m; uint64_t *a; } uint64_v;
typedef struct { size_t n, m; uint32_t *a; } uint32_v;

typedef struct {
	mm128_v a;   // (minimizer, position) array
	int32_t n;   // size of the _p_ array
	uint64_t *p; // position array for minimizers appearing >1 times
	void *h;     // hash table indexing _p_ and minimizers appearing once
} mm_idx_bucket_t;

typedef struct {
	int b, w, k;
	uint32_t n;   // number of reference sequences
	mm_idx_bucket_t *B;
	uint32_t max_occ;
	float freq_thres;
	int32_t *len; // length of each reference sequence
	char **name;  // TODO: if this uses too much RAM, switch to one concatenated string
} mm_idx_t;

typedef struct {
	uint32_t cnt:31, rev:1;
	uint32_t rid:31, rep:1;
	uint32_t len;
	int32_t qs, qe, rs, re;
} mm_reg1_t;

typedef struct {
	int radius;  // bandwidth to cluster hits
	int max_gap; // break a chain if there are no minimizers in a max_gap window
	int min_cnt; // minimum number of minimizers to start a chain
	int min_match;
	int sdust_thres;  // score threshold for SDUST; 0 to disable
	int flag;    // see MM_F_* macros
	float merge_frac; // merge two chains if merge_frac fraction of minimzers are shared between the chains
} mm_mapopt_t;

extern int mm_verbose;
extern double mm_realtime0;

struct mm_tbuf_s;
typedef struct mm_tbuf_s mm_tbuf_t;

#ifdef __cplusplus
extern "C" {
#endif

// compute minimizers
void mm_sketch(void *km, const char *str, int len, int w, int k, uint32_t rid, int is_hpc, mm128_v *p);

// minimizer indexing
mm_idx_t *mm_idx_init(int w, int k, int b);
void mm_idx_destroy(mm_idx_t *mi);
mm_idx_t *mm_idx_gen(bseq_file_t *fp, int w, int k, int b, int tbatch_size, int n_threads, uint64_t ibatch_size, int keep_name);
void mm_idx_set_max_occ(mm_idx_t *mi, float f);
const uint64_t *mm_idx_get(const mm_idx_t *mi, uint64_t minier, int *n);

mm_idx_t *mm_idx_build(const char *fn, int w, int k, int n_threads);

// minimizer index I/O
void mm_idx_dump(FILE *fp, const mm_idx_t *mi);
mm_idx_t *mm_idx_load(FILE *fp);

// mapping
void mm_mapopt_init(mm_mapopt_t *opt);
mm_tbuf_t *mm_tbuf_init(void);
void mm_tbuf_destroy(mm_tbuf_t *b);
const mm_reg1_t *mm_map(const mm_idx_t *mi, int l_seq, const char *seq, int *n_regs, mm_tbuf_t *b, const mm_mapopt_t *opt, const char *name);

int mm_map_file(const mm_idx_t *idx, const char *fn, const mm_mapopt_t *opt, int n_threads, int tbatch_size);

// private functions (may be moved to a "mmpriv.h" in future)
double cputime(void);
double realtime(void);
void radix_sort_128x(mm128_t *beg, mm128_t *end);
void radix_sort_64(uint64_t *beg, uint64_t *end);
uint32_t ks_ksmall_uint32_t(size_t n, uint32_t arr[], size_t kk);

#ifdef __cplusplus
}
#endif

#endif
