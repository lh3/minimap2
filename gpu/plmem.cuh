#ifndef _PLMEM_CUH_
#define _PLMEM_CUH_
#include "hipify.cuh"
#include "plchain.h"
#include "plutils.h"

#ifndef MAX_MICRO_BATCH
#define MAX_MICRO_BATCH 8
#endif // MAX_MICRO_BATCH

#define OneK 1024
#define OneM (OneK*1024)
#define OneG (OneM*1024)


typedef struct {
    int index;       // read index / batch index
    int griddim;     // grid for range selection kernel. 
    int size;        // number of reads in the batch
    size_t total_n;  // number of anchors in the batch
    size_t cut_num;  // number of cuts in the batch

    // array size: number of anchors in the batch
    int32_t *ax;  // (int32_t) a[].x
    int32_t *ay;  // (int32_t) a[].y
    int8_t* sid;  // a[].y >> 40 & 0xff
    int32_t *xrev; // a[].x >> 32
    // outputs
    int32_t *f;   // score
    uint16_t *p;  // predecessor

    // array size: number of cuts in the batch / long_seg_cut
    // total long segs number till this batch
    unsigned int *long_segs_num;

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

typedef struct {
    // array size: number of cuts in the batch / long_seg_cut
    seg_t *long_segs_og_idx;                   // start & end idx of long segs in the original micro batch
    unsigned int *total_long_segs_num; // sum of mini batch long_segs_num
    size_t *total_long_segs_n; // number of anchors in all the long segs
    int32_t *f_long;   // score for long segs
    uint16_t *p_long;  // predecessor for long segs
} longMemPtr;

typedef struct {
    int size;
    int griddim;
    size_t total_n;
    size_t num_cut;
    // device memory ptrs
    // data array
    int32_t *d_ax;
    int32_t *d_ay;
    int8_t *d_sid;  // a[].y >> 40 & 0xff
    int32_t *d_xrev; // a[].x >> 32
    int32_t *d_range;
    int32_t *d_f;   // score
    uint16_t *d_p;  // predecessor

    // range selection index
    size_t *d_start_idx;
    size_t *d_read_end_idx;
    size_t *d_cut_start_idx;

    // cut
    size_t *d_cut;  // cut
    unsigned int *d_long_seg_count; // total number of long seg (aggregated accross micro batches)
    seg_t *d_long_seg;              // start & end idx of long segs in the long seg buffer (aggregated across micro batches)
    seg_t *d_long_seg_og;           // start & end idx of long seg in the micro batch. (aggregated accross micro batches)
    unsigned int *d_mid_seg_count;  // private to micro batch
    seg_t *d_mid_seg;               // private to micro batch

    // long segement buffer
    unsigned *d_map;
    int32_t *d_ax_long, *d_ay_long;
    int8_t *d_sid_long;
    int32_t *d_range_long;
    size_t *d_total_n_long;
    size_t buffer_size_long;
    int32_t *d_f_long;  // score, size: buffer_size_long * sizeof(int32_t)
    uint16_t *d_p_long;  // predecessor, size: buffer_size_long * sizeof(uint16_t)
} deviceMemPtr;

typedef struct stream_ptr_t{
    chain_read_t *reads;
    size_t n_read;
    hostMemPtr host_mems[MAX_MICRO_BATCH];
    longMemPtr long_mem;
    deviceMemPtr dev_mem;
    cudaStream_t cudastream;
    cudaEvent_t stopevent, startevent, long_kernel_event;
    cudaEvent_t short_kernel_start_event[MAX_MICRO_BATCH];
    cudaEvent_t short_kernel_stop_event[MAX_MICRO_BATCH];
    bool busy = false;
} stream_ptr_t;

typedef struct gputSetup_t {
    int num_stream;
    stream_ptr_t *streams;
    size_t max_anchors_stream, max_num_cut, long_seg_buffer_size_stream;
    int max_range_grid;
} streamSetup_t;

extern streamSetup_t stream_setup;

/* memory management methods */
// initialization and cleanup
void plmem_initialize(size_t *max_total_n, int *max_read, int *min_n);
void plmem_stream_initialize(size_t *max_total_n, int *max_read, int *min_n, char* gpu_config_file);
void plmem_stream_cleanup();

// alloc and free
void plmem_malloc_host_mem(hostMemPtr *host_mem, size_t anchor_per_batch,
                           int range_grid_size, size_t buffer_size_long);
void plmem_malloc_long_mem(longMemPtr *long_mem, size_t buffer_size_long);
void plmem_free_host_mem(hostMemPtr *host_mem);
void plmem_free_long_mem(longMemPtr *long_mem);
void plmem_malloc_device_mem(deviceMemPtr *dev_mem, size_t anchor_per_batch,
                             int range_grid_size, int num_cut);
void plmem_free_device_mem(deviceMemPtr *dev_mem);

// data movement
void plmem_reorg_input_arr(chain_read_t *reads, int n_read,
                           hostMemPtr *host_mem, range_kernel_config_t config);
void plmem_async_h2d_memcpy(stream_ptr_t *stream_ptrs);
void plmem_async_h2d_short_memcpy(stream_ptr_t *stream_ptrs, size_t uid);
void plmem_sync_h2d_memcpy(hostMemPtr *host_mem, deviceMemPtr *dev_mem);
void plmem_async_d2h_memcpy(stream_ptr_t *stream_ptrs);
void plmem_async_d2h_short_memcpy(stream_ptr_t *stream_ptrs, size_t uid);
void plmem_async_d2h_long_memcpy(stream_ptr_t *stream_ptrs);
void plmem_sync_d2h_memcpy(hostMemPtr *host_mem, deviceMemPtr *dev_mem);
#endif  // _PLMEM_CUH_