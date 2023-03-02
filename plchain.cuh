#ifndef _PLCHAIN_CUH_
#define _PLCHAIN_CUH_

#include "hipify.cuh"
#include "plutils.h"
#include <assert.h>


typedef struct {
    int index;       // read index / batch index
    int griddim;     // grid for range selection kernel. 
    int size;        // number of reads in the batch
    size_t total_n;  // number of anchors in the batch
    size_t cut_num;  // number of cuts in the batch

    // array size: number of anchors in the batch
    int64_t *ax;  // a[].x
    int64_t *ay;  // a[].y
    int32_t *f;   // score
    uint16_t *p;  // predecessor

    // start index for each block in range selection
    /***** range selection block assiagnment
     * One block only gets assgined one read or part of one read.
     *  start_idx:      idx of the first anchor assigned to each block
     *  read_end_idx:   idx of the last anchor OF THE READ assigned to each
     * block if a read is devided into several blocks, all the blocks take the
     * last anchor index of the read cut_start_idx:  idx of the first cut this
     * block needs to make
     */
    // array size: grid dimension
    size_t *start_idx;
    size_t *read_end_idx;
    size_t *cut_start_idx;
} hostMemPtr;

typedef struct seg_t {
    size_t start_idx;
    size_t end_idx;
} seg_t;

typedef struct {
    int size;
    int griddim;
    size_t total_n;
    size_t num_cut;
    // device memory ptrs
    // data array
    int64_t *d_ax;
    int64_t *d_ay;
    int32_t *d_range;
    int32_t *d_f;   // score
    uint16_t *d_p;  // predecessor

    // range selection index
    size_t *d_start_idx;
    size_t *d_read_end_idx;
    size_t *d_cut_start_idx;

    // cut
    size_t *d_cut;  // cut
    unsigned int *d_long_seg_count;
    seg_t *d_long_seg;
} deviceMemPtr;

typedef struct stream_ptr_t{
    chain_read_t *reads;
    hostMemPtr host_mem;
    deviceMemPtr dev_mem;
    cudaStream_t cudastream;
    cudaEvent_t cudaevent;
    bool busy = false;
} stream_ptr_t;

typedef struct gputSetup_t {
    int num_stream;
    stream_ptr_t *streams;
    size_t max_anchors_stream, max_range_grid, max_num_cut;
} streamSetup_t;

extern streamSetup_t stream_setup;

#include "plchain.cuh"

/* memory management methods */
// initialization and cleanup
void plmem_initialize(size_t *max_total_n, size_t *max_read, size_t *min_n);
void plmem_stream_initialize(size_t *max_total_n, size_t *max_read, size_t *min_n);
void plmem_stream_cleanup();

// alloc and free
void plmem_malloc_host_mem(hostMemPtr *host_mem, size_t anchor_per_batch,
                           int range_grid_size);
void plmem_free_host_mem(hostMemPtr *host_mem);
void plmem_malloc_device_mem(deviceMemPtr *dev_mem, size_t anchor_per_batch,
                             int range_grid_size, int num_cut);
void plmem_free_device_mem(deviceMemPtr *dev_mem);

// data movement
void plmem_reorg_input_arr(chain_read_t *reads, int n_read,
                           hostMemPtr *host_mem, range_kernel_config_t config);
void plmem_async_h2d_memcpy(stream_ptr_t *stream_ptrs);
void plmem_sync_h2d_memcpy(hostMemPtr *host_mem, deviceMemPtr *dev_mem);
void plmem_async_d2h_memcpy(stream_ptr_t *stream_ptrs);
void plmem_sync_d2h_memcpy(hostMemPtr *host_mem, deviceMemPtr *dev_mem);

#ifdef __cplusplus
extern "C" {
#endif

typedef __int32_t int32_t;

/* functions declaration */
void plrange_upload_misc(Misc misc);
void plrange_async_range_selection(deviceMemPtr* device_mem_ptr, cudaStream_t* stream);
void plrange_sync_range_selection(deviceMemPtr* dev_mem, Misc misc
#ifdef DEBUG_CHECK
                                  ,
                                  chain_read_t* reads
#endif // DEBUG_CHECK
);

extern range_kernel_config_t range_kernel_config;

void plscore_upload_misc(Misc misc);
void plscore_async_naive_forward_dp(deviceMemPtr* dev_mem, cudaStream_t* stream);
void plscore_async_long_short_forward_dp(deviceMemPtr* dev_mem,cudaStream_t* stream);
void plscore_sync_long_short_forward_dp(deviceMemPtr* dev_mem, Misc misc_);
void plscore_sync_naive_forward_dp(deviceMemPtr* dev_mem, Misc misc_);

extern score_kernel_config_t score_kernel_config;

// helper function
void post_chaining_helper(chain_read_t *read,
                          const mm_seq_meta_t *target,
                          const mm_seq_meta_t *query, int max_dist_t, void *km);

#ifdef __cplusplus
}
#endif


#endif // _PLCHAIN_CUH_
