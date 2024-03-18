#ifndef _PLCHAIN_H_
#define _PLCHAIN_H_

/* Range Kernel configuaration */
typedef struct range_kernel_config_t {
    int blockdim;           // number of threads in each block
    int cut_check_anchors;  // number of anchors to check around each cut
    int anchor_per_block;   // number of anchors assgined to one block = max_it * blockdim
} range_kernel_config_t;

/* Score Generation Kernel configuration */
typedef struct score_kernel_config_t{
    int micro_batch;
    int short_blockdim;
    int long_blockdim;
    int mid_blockdim;
    int short_griddim;
    int long_griddim;
    int mid_griddim;
    int cut_unit;
    int long_seg_cutoff;
    int mid_seg_cutoff;
} score_kernel_config_t;

#endif // _PLCHAIN_H_