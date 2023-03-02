#ifndef _PLUTILS_H_ 
#define _PLUTILS_H_ 

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <limits.h>
#include <assert.h>
#include "minimap.h"

#include "kalloc.h"

#define MM_SEED_SEG_SHIFT  48
#define MM_SEED_SEG_MASK   (0xffULL<<(MM_SEED_SEG_SHIFT))

/* Range Kernel configuaration */
typedef struct range_kernel_config_t {
    int blockdim;           // number of threads in each block
    int cut_check_anchors;  // number of anchors to check around each cut
    int anchor_per_block;   // number of anchors assgined to one block = max_it * blockdim
} range_kernel_config_t;

/* Score Generation Kernel configuration */
typedef struct score_kernel_config_t{
    int short_blockdim;
    int long_blockdim;
    int cut_per_block;
} score_kernel_config_t;

/* Chaining Options */

// Chaining Flags

typedef struct {
    char method[20]; // "forward_cpu", "backward_cpu", "blocking_gpu", "stream_gpu"
    char *gpu_cfg; // path to gpu json file
    int64_t flag;  // see MM_F_** macros
    int k;

    int bw, bw_long;           // bandwidth
    int max_gap, max_gap_ref;  // break a chain if there are no minimizers in a
                               // max_gap window
    int max_frag_len;
    int max_chain_skip, max_chain_iter;
    int min_cnt;          // min number of minimizers on each chain
    int min_chain_score;  // min chaining score
    float chain_gap_scale;
    float chain_skip_scale;

    // primary chain identification parameters
    float mask_level;
    int mask_len;
    int match_sc, mismatch_sc; // a,b in minimap2
    float pri_ratio;
    int best_n;

    float alt_drop;

    // rmq chaining parameters
    int rmq_size_cap, rmq_inner_dist;
    int rmq_rescue_size;
    float rmq_rescue_ratio;

    int seed; // random seed for hash

    int n_thread;

} chain_opt_t;
extern chain_opt_t opt;

/* structure for metadata and hits */
// Sequence meta data
typedef struct {
    char name[200]; // name of the sequence
    uint32_t len; // name of the sequence
    int is_alt;
    int n_alt;
    int rep_len;
} mm_seq_meta_t;

typedef struct {
    int max_iter, max_dist_x, max_dist_y, max_skip, bw, min_cnt, min_score, is_cdna, n_seg;
    float chn_pen_gap, chn_pen_skip;
} Misc;

typedef struct {
    mm_seq_meta_t *refs;
    int n_refs;
    Misc misc;
} input_meta_t;

typedef struct {
    void *km;

    mm_seq_meta_t seq;
    mm128_t *a;     // anchors
    int64_t n;      // number of anchors
    Misc misc;
    uint64_t *u;    // chains
    int64_t n_u;    // number of chains
    mm_reg1_t *regs;    // hits
    int n_reg;          // number of hits

    // for score checking after chaining
    int32_t *f;
    int64_t *p;
} chain_read_t;


/////////////////////////////////////////////////////
///////////         Free Input Struct   /////////////
/////////////////////////////////////////////////////
// free input_iter pointers except a, because it is freed seperately. 
static inline void free_read(chain_read_t *in) {
    void *km = in->km;
    if (in->a) kfree(km, in->a);
    if (in->u) kfree(km, in->u);
    if (in->regs) kfree(km, in->regs);
    if (in->f) kfree(km, in->f);
    if (in->p) kfree(km, in->p);
    in->a = 0, in->u = 0, in->regs = 0;
    in->f = 0, in->p = 0;
}

static inline void free_meta_struct(input_meta_t *meta, void* km) {
    if (meta->refs)
        kfree(km, meta->refs);
}
#endif // _PLUTILS_H_ 

