#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include "plrange.cuh"
#include "hipify.cuh"


/* 

CUDA/HIP kernel for range selection using forward chaining

*/

/* kernels begin */
__constant__ int d_max_dist_x;
__constant__ int d_max_iter;
__constant__ int d_cut_check_anchors;

inline __device__ int64_t range_binary_search(const int32_t* ax, const int32_t* rev, int64_t i, int64_t st_end){
    int64_t st_high = st_end, st_low=i;
    while (st_high != st_low) {
        int64_t mid = (st_high + st_low -1) / 2+1;
        if (rev[i] != rev[mid] || ax[mid] > ax[i] + d_max_dist_x) {
            st_high = mid -1;
        } else {
            st_low = mid;
        }
    }
    return st_high;
}


/**
 * Forward Range Selection Kernel using global memory and binary range search. 
 * cut reads into segements where successor range = 0. 
*/
__global__ void range_selection_kernel_binary(const int32_t* ax, const int32_t* rev, size_t *start_idx_arr, size_t *read_end_idx_arr, 
    int32_t *range, size_t* cut, size_t* cut_start_idx, size_t total_n, range_kernel_config_t config){
    int tid = threadIdx.x;
    int bid = blockIdx.x;

    size_t start_idx = start_idx_arr[bid];
    size_t read_end_idx = read_end_idx_arr[bid];
    size_t end_idx = start_idx + config.anchor_per_block;
    end_idx = end_idx > read_end_idx ? read_end_idx : end_idx;
    size_t cut_idx = cut_start_idx[bid];
    if(tid == 0 && (bid == 0 || read_end_idx_arr[bid-1] != read_end_idx)){
        cut[cut_idx] = start_idx;
    }
    cut_idx++;
    int range_op[3] = {16, 512, 5000};  // Range Options
    range_op[2] = d_max_iter;
    for (size_t i = start_idx + tid; i < end_idx; i += blockDim.x) {
        size_t st_max = i + d_max_iter;
        st_max = st_max < read_end_idx ? st_max : read_end_idx -1;
        size_t st;
        for (int j=0; j<3; ++j){
            st = i + range_op[j];
            st = st <= st_max ? st : st_max;
            assert(st < total_n);
            assert(i < total_n);
            if (st > i && (rev[st] != rev[i] || ax[st] > ax[i] + d_max_dist_x)){
                break;
            }
        }
        st = range_binary_search(ax, rev, i, st);
        range[i] = st - i;

        if (tid >= blockDim.x - d_cut_check_anchors &&
            blockDim.x - tid + i <= end_idx) {
            if (st == i) cut[cut_idx] = i+1;
        }
        cut_idx++;
    }
}

/**
 * Forward Range Selection Kernel using global memory and linear range search.
 * cut reads into segements where successor range = 0.
 */
__global__ void range_selection_kernel_naive(const int32_t* ax, const int32_t* rev, size_t *start_idx_arr, size_t *read_end_idx_arr, 
    int32_t *range, size_t* cut, size_t* cut_start_idx, size_t total_n, range_kernel_config_t config){
    int tid = threadIdx.x;
    int bid = blockIdx.x;

    size_t start_idx = start_idx_arr[bid];
    size_t read_end_idx = read_end_idx_arr[bid];
    size_t end_idx = start_idx + config.anchor_per_block;
    end_idx = end_idx > read_end_idx ? read_end_idx : end_idx;
    assert(end_idx == (bid +1 < gridDim.x) ? start_idx_arr[bid+1]: total_n);
    // if(end_idx_ref != end_idx){
    //     if (tid == 0){
    //         int grimdim = gridDim.x;
    //     printf("start idx %d anchor_per_block %d read_end_idx %d, next start idx %d gridDim %d bid %d\n", 
    //     start_idx, config.anchor_per_block, read_end_idx, start_idx_arr[bid+1], grimdim, bid);
    //     }
    // }
    // __syncthreads();
    
    size_t cut_idx = cut_start_idx[bid];
    if(tid == 0 && (bid == 0 || read_end_idx_arr[bid-1] != read_end_idx)){
        cut[cut_idx] = start_idx;
    }
    cut_idx++;
    for (size_t i = start_idx + tid; i < end_idx; i += blockDim.x){
        size_t st = i + d_max_iter;
        st = i + d_max_iter < read_end_idx ? st : read_end_idx -1;
        assert(st < total_n);
        assert(i < total_n);
        while (st > i && 
                (rev[i] != rev[st] // NOTE: different prefix cannot become predecessor 
                || ax[st] > ax[i] + d_max_dist_x)) { // NOTE: same prefix compare the value
            --st;
        }
        range[i] = st - i;

        if (tid >= blockDim.x - d_cut_check_anchors && blockDim.x - tid + i <= end_idx) {
            if (st == i) cut[cut_idx] = i+1;
        }
        cut_idx++;
    }
}

// __global__ void range_selection_kernel(const int64_t* ax, size_t *start_idx_arr, size_t *read_end_idx_arr, int32_t *range){
//     int tid = threadIdx.x;
//     int bid = blockIdx.x;

//     size_t start_idx = start_idx_arr[bid];
//     size_t read_end_idx = read_end_idx_arr[bid];
//     size_t end_idx = start_idx + MAX_ANCHOR_PER_BLOCK;
//     end_idx = end_idx > read_end_idx ? read_end_idx : end_idx;

//     size_t load_anchor_idx = 100;
//     size_t load_smem_idx;
//     size_t cal_idx = start_idx + threadIdx.x;
//     int32_t cal_smem = tid;
//     __shared__ int64_t smem[NUM_ANCHOR_IN_SMEM];

//     /* prefetch anchors */
//     load_smem_idx = tid;
//     load_anchor_idx = start_idx + tid;
//     // if (tid == 20) printf("load_smem_idx %d, load_anchor_idx %lu\n", load_smem_idx, load_anchor_idx);
//     for (int i = 0; i < PREFETCH_ANCHORS_RANGE/NUM_THREADS_RANGE && load_anchor_idx < read_end_idx; ++i){
//         // if (tid == 20) printf("load_smem_idx %d, load_anchor_idx %lu\n", load_smem_idx, load_anchor_idx);
//         smem[load_smem_idx] = ax[load_anchor_idx];
//         load_smem_idx += NUM_THREADS_RANGE;
//         load_anchor_idx += NUM_THREADS_RANGE;
//     }

//     int iter = (NUM_ANCHOR_IN_SMEM - PREFETCH_ANCHORS_RANGE)/NUM_THREADS_RANGE; // iterations before another load is needed
//     while (cal_idx < end_idx) { // tail threads may skip this loop
//         /* load anchors */
//         load_smem_idx = load_smem_idx >= NUM_ANCHOR_IN_SMEM ? load_smem_idx - NUM_ANCHOR_IN_SMEM : load_smem_idx;
//         for (int i = 0; i < iter && load_anchor_idx < end_idx + PREFETCH_ANCHORS_RANGE; ++i){
//             // if (tid == 20) printf("load it load_smem_idx %d, load_anchor_idx %lu\n", load_smem_idx, load_anchor_idx);
//             smem[load_smem_idx] = ax[load_anchor_idx];
//             load_smem_idx += NUM_THREADS_RANGE;
//             load_anchor_idx += NUM_THREADS_RANGE;
//             load_smem_idx = load_smem_idx >= NUM_ANCHOR_IN_SMEM ? load_smem_idx - NUM_ANCHOR_IN_SMEM : load_smem_idx;
//         }

//         __syncthreads();
        
//         /* calculate sucessor range */
//         for (int i = 0; i < iter && cal_idx < end_idx; ++i){
//             int64_t anchor = smem[cal_smem];

//             size_t st = cal_idx + PREFETCH_ANCHORS_RANGE < read_end_idx ? cal_idx + PREFETCH_ANCHORS_RANGE : read_end_idx-1;
//             int32_t st_smem = cal_smem + st - cal_idx;
//             st_smem = st_smem >= NUM_ANCHOR_IN_SMEM ? st_smem - NUM_ANCHOR_IN_SMEM : st_smem;
//             // if (tid == 20) printf("cal idx %lu, cal_mem %d, st %lu, st_smem %d\n", cal_idx, cal_smem, st,st_smem);

//             // if (tid == 20) printf("anchor.x %d, smem[st_smem] %d, anchor.x+MAX_DIST_X%d\n", anchor, smem[st_smem], anchor+MAX_DIST_X);

//             while (st > cal_idx && 
//                         (anchor>> 32 != smem[st_smem] >> 32 ||
//                             smem[st_smem] > anchor + d_max_dist_x
//                         )
//                     ){
//                 // if (bid == 25)
//                 // printf("while 0 bid %d tid %d cal_idx %d\n", bid, tid, cal_idx);
//                 --st;
//                 if (st_smem == 0) st_smem = NUM_ANCHOR_IN_SMEM-1;
//                 else --st_smem;
//             }
            
//             /* NOTE: fallback: succussor is not prefetched */
//             if (st >= PREFETCH_ANCHORS_RANGE + cal_idx){
//                 st = cal_idx + MAX_ITER < read_end_idx ? i + MAX_ITER : read_end_idx-1;
//                 while(
//                     anchor >> 32 != ax[st] >> 32 || 
//                     ax[st] > anchor + d_max_dist_x // check from global memory
//                 ){
//                     --st;
//                     // if (bid == 25)
//                     // printf("while 1 bid %d tid %d\n", bid, tid);
//                 }

//             }
//             range[cal_idx] = st - cal_idx;
//             cal_smem += NUM_THREADS_RANGE;
//             cal_smem = cal_smem >= NUM_ANCHOR_IN_SMEM ? cal_smem - NUM_ANCHOR_IN_SMEM : cal_smem;
//             cal_idx += NUM_THREADS_RANGE;
//             // if (bid == 25)
//             // printf("for loop i %d bid %d tid %d\n", i, bid, tid);
//         }
//         // if (bid == 25)
//         // printf("outer while bid %d tid %d\n", bid, tid);
//         __syncthreads();

//     }
    
// }

/* kernels end */

#ifdef __cplusplus
extern "C" {
#endif

/* host functions begin */
range_kernel_config_t range_kernel_config;

void plrange_upload_misc(Misc misc){
#ifdef USEHIP
    hipMemcpyToSymbol(HIP_SYMBOL(d_max_dist_x), &misc.max_dist_x, sizeof(int));
    hipMemcpyToSymbol(HIP_SYMBOL(d_max_iter), &misc.max_iter, sizeof(int));
    hipMemcpyToSymbol(HIP_SYMBOL(d_cut_check_anchors),
                      &range_kernel_config.cut_check_anchors, sizeof(int));
#else
    cudaCheck();
    cudaMemcpyToSymbol(d_max_dist_x, &misc.max_dist_x, sizeof(int));
    cudaMemcpyToSymbol(d_max_iter, &misc.max_iter, sizeof(int));
    cudaMemcpyToSymbol(d_cut_check_anchors,
                       &range_kernel_config.cut_check_anchors, sizeof(int));
#endif  // USEHIP
    cudaCheck();
}

void plrange_async_range_selection(deviceMemPtr* dev_mem, cudaStream_t* stream) {
    size_t total_n = dev_mem->total_n, cut_num = dev_mem->num_cut;
    int griddim = dev_mem->griddim;
    dim3 DimBlock(range_kernel_config.blockdim, 1, 1);
    dim3 DimGrid(griddim, 1, 1);

    // Run kernel
    range_selection_kernel_binary<<<DimGrid, DimBlock, 0, *stream>>>(
        dev_mem->d_ax, dev_mem->d_xrev, dev_mem->d_start_idx, dev_mem->d_read_end_idx,
        dev_mem->d_range, dev_mem->d_cut, dev_mem->d_cut_start_idx, total_n, range_kernel_config);
    cudaCheck();
#ifdef DEBUG_PRINT
    // fprintf(stderr, "[Info] %s (%s:%d): Batch total_n %lu, Range Kernel Launched, grid %d cut %d\n", __func__, __FILE__, __LINE__, total_n, DimGrid.x, cut_num);
#endif
}

void plrange_sync_range_selection(deviceMemPtr *dev_mem, Misc misc) {
    size_t total_n = dev_mem->total_n, cut_num = dev_mem->num_cut;
    int griddim = dev_mem->griddim;
    dim3 DimBlock(range_kernel_config.blockdim, 1, 1);
    dim3 DimGrid(griddim,1,1);

    plrange_upload_misc(misc);

    // Run kernel
#ifdef DEBUG_PRINT
        fprintf(stderr, "[Info] %s (%s:%d): Grim Dim: %d Cut: %zu Anchors: %zu\n", __func__, __FILE__, __LINE__, DimGrid.x,
                cut_num, total_n);
#endif
    range_selection_kernel_binary<<<DimGrid, DimBlock>>>(
        dev_mem->d_ax, dev_mem->d_xrev, dev_mem->d_start_idx, dev_mem->d_read_end_idx,
        dev_mem->d_range, dev_mem->d_cut, dev_mem->d_cut_start_idx, total_n, range_kernel_config);
    cudaCheck();
    cudaDeviceSynchronize();
    cudaCheck();
#ifdef DEBUG_PRINT
    fprintf(stderr, "[Info] %s: range calculation success\n", __func__);
#endif
}

#ifdef __cplusplus
}
#endif

/* host functions end */
