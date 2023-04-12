#ifndef _PLUTILS_H_
#define _PLUTILS_H_

#include <assert.h>
#include <limits.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "kalloc.h"
#include "minimap.h"

/* Chaining Options */

/* structure for metadata and hits */
// Sequence meta data
typedef struct {
    long i;          // read id
    int seg_id;      // seg id
    char name[200];  // name of the sequence
    uint32_t len;    // name of the sequence

    // mi data
    int n_alt;
    int is_alt;  // reference sequences only

    // sequence info
    int qlen_sum;
} mm_seq_meta_t;

typedef struct {
    int max_iter, max_dist_x, max_dist_y, max_skip, bw, min_cnt, min_score,
        is_cdna, n_seg;
    float chn_pen_gap, chn_pen_skip;
} Misc;

typedef struct {
    mm_seq_meta_t *refs;
    int n_refs;
    Misc misc;
} input_meta_t;

typedef struct {
    mm_seq_meta_t seq;

    // minimap2 input data for reads
    const char **qseqs;  // sequences for each segment          <- allocated in worker_for, freed in free_read after seeding
    int *qlens;          // query length for each segment       <- allocated in worker_for, freed in free_read after seeding
    int n_seg;           // number of segs

#ifdef DEBUG_CHECK
    // for score checking after chaining
    int32_t *f;
    int64_t *p;
#endif  // DEBUG_CHECK
    int rep_len;
    int frag_gap;

    // seeding outputs
    uint64_t *mini_pos;  // minimizer positions                 <- allocated in 
    int n_mini_pos;

    // seeding output, updated in chaining
    mm128_t *a;  // array of anchors
    int64_t n;   // number of anchors = n_a

    // chaining outputs
    uint64_t *u;      // scores for chains
    int n_u;          // number of chains formed from anchors == n_reg0

} chain_read_t;

#ifdef __cplusplus
extern "C" {
#endif // __cplusplus
/* GPU chaining methods */
// // <lchain.h> backward, original chaining methods
// void chain_backword_cpu(const input_meta_t *meta, chain_read_t *read_arr,
//                         int n_read);
// // <fchain.c> forward chaining methods
// void chain_forward_cpu(const input_meta_t *meta, chain_read_t *read_arr,
//                        int n_read);

// <plchain.cu> gpu chaining methods
// initialization and cleanup
void init_blocking_gpu(size_t *max_total_n, int *max_reads,
                       int *min_n, Misc misc);  // for blocking_gpu
void init_stream_gpu(size_t *max_total_n, int *max_reads,
                     int *min_n, Misc misc);  // for stream_gpu
void finish_stream_gpu(const mm_idx_t *mi, const mm_mapopt_t *opt, chain_read_t **batches,
                        int *num_reads, int num_batch, void *km);  // for stream_gpu
void free_stream_gpu(int n_threads); // for stream_gpu free pinned memory
// chaining method
void chain_blocking_gpu(const mm_idx_t *mi, const mm_mapopt_t *opt, chain_read_t *in_arr, int n_read, void* km);
void chain_stream_gpu(const mm_idx_t *mi, const mm_mapopt_t *opt, chain_read_t **in_arr_ptr, int *n_read_ptr, int thread_id, void* km);

/* <lchain.c> Chaining backtracking methods */
uint64_t *mg_chain_backtrack(void *km, int64_t n, const int32_t *f,
                             const int64_t *p, int32_t *v, int32_t *t,
                             int32_t min_cnt, int32_t min_sc, int32_t max_drop,
                             int32_t *n_u_, int32_t *n_v_);
mm128_t *compact_a(void *km, int32_t n_u, uint64_t *u, int32_t n_v, int32_t *v, mm128_t *a);


/* Post Chaining helpers */
Misc build_misc(const mm_idx_t *mi, const mm_mapopt_t *opt, const int64_t qlen_sum, const int n_seg);
void post_chaining_helper(const mm_idx_t *mi, const mm_mapopt_t *opt,
                          chain_read_t *read, Misc misc, void *km);

#ifdef __cplusplus
}
#endif  // __cplusplus

/////////////////////////////////////////////////////
///////////         Free Input Struct   /////////////
/////////////////////////////////////////////////////
// free input_iter pointers except a, because it is freed seperately.
static inline void free_read(chain_read_t *in, void* km) {
    if (in->qseqs) kfree(km, in->qseqs);
    if (in->qlens) kfree(km, in->qlens);
    // if (in->a) kfree(km, in->a);
    // if (in->u) kfree(km, in->u);
#ifdef DEBUG_CHECK
    if (in->f) kfree(km, in->f);
    if (in->p) kfree(km, in->p);
    in->f = 0, in->p = 0;
#endif
    in->qseqs = 0, in->qlens = 0;
    in->a = 0, in->u = 0;
}

static inline void free_meta_struct(input_meta_t *meta, void *km) {
    if (meta->refs) kfree(km, meta->refs);
}
#endif  // _PLUTILS_H_
