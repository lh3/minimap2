#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <time.h>


#include "mmpriv.h"
#include "plmem.cuh"
#include "plrange.cuh"
#include "plscore.cuh"
#include "plchain.h"

/**
 * translate relative predecessor index to abs index 
 * Input
 *  rel[]   relative predecessor index
 * Output
 *  p[]     absolute predecessor index (of each read)
 */
void p_rel2idx(const uint16_t* rel, int64_t* p, size_t n) {
    for (int i = 0; i < n; ++i) {
        if (rel[i] == 0)
            p[i] = -1;
        else
            p[i] = i - rel[i];
    }
}

//////////////////////////////////////////////////////////////////////////
///////////         Backtracking    //////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

/**
 * @brief start from end index of the chain, find the location of min score on
 * the chain until anchor has no predecessor OR anchor is in another chain
 *
 * @param max_drop
 * @param z [in] {sc, anchor idx}, sorted by sc
 * @param f [in] score
 * @param p [in] predecessor
 * @param k [in] chain end index
 * @param t [update] 0 for unchained anchor, 1 for chained anchor
 * @return min_i minmute score location in the chain
 */

static int64_t mg_chain_bk_end(int32_t max_drop, const mm128_t *z,
                               const int32_t *f, const int64_t *p, int32_t *t,
                               int64_t k) {
    int64_t i = z[k].y, end_i = -1, max_i = i;
    int32_t max_s = 0;
    if (i < 0 || t[i] != 0) return i;
    do {
        int32_t s;
        t[i] = 2;
        end_i = i = p[i];
        s = i < 0 ? z[k].x : (int32_t)z[k].x - f[i];
        if (s > max_s)
            max_s = s, max_i = i;
        else if (max_s - s > max_drop)
            break;
    } while (i >= 0 && t[i] == 0);
    for (i = z[k].y; i >= 0 && i != end_i; i = p[i])  // reset modified t[]
        t[i] = 0;
    return max_i;
}

void plchain_backtracking(hostMemPtr *host_mem, chain_read_t *reads, Misc misc, void* km){
    int max_drop = misc.bw;
    if (misc.max_dist_x < misc.bw) misc.max_dist_x = misc.bw;
    if (misc.max_dist_y < misc.bw && !misc.is_cdna) misc.max_dist_y = misc.bw;
    if (misc.is_cdna) max_drop = INT32_MAX;

    size_t n_read = host_mem->size;

    uint16_t* p_hostmem = host_mem->p;
    int32_t* f = host_mem->f;
    for (int i = 0; i < n_read; i++) {
        int64_t* p;
        KMALLOC(km, p, reads[i].n);
        p_rel2idx(p_hostmem, p, reads[i].n);
#ifdef DEBUG_VERBOSE
        debug_print_score(p, f, reads[i].n);
#endif
#ifdef DEBUG_CHECK
        debug_check_score(p, f, reads[i].p, reads[i].f, reads[i].n);
#endif

        /* Backtracking */
        uint64_t* u;
        int32_t *v, *t;
        KMALLOC(km, v, reads[i].n);
        KCALLOC(km, t, reads[i].n);
        int32_t n_u, n_v;
        u = mg_chain_backtrack(km, reads[i].n, f, p, v, t, misc.min_cnt, misc.min_score, max_drop, &n_u, &n_v);
        reads[i].u = u;
        reads[i].n_u = n_u;
        kfree(km, p);
        // here f is not managed by km memory pool
        kfree(km, t);
        if (n_u == 0) {
            kfree(km, reads[i].a);
            kfree(km, v);
            reads[i].a = 0;

            f += reads[i].n;
            p_hostmem += reads[i].n;
            continue;
        }

        mm128_t* new_a = compact_a(km, n_u, u, n_v, v, reads[i].a);
        reads[i].a = new_a;

        f += reads[i].n;
        p_hostmem += reads[i].n;
    }
}

void plchain_cal_score_sync(chain_read_t *reads, int n_read, Misc misc, void* km) { 
    hostMemPtr host_mem;
    deviceMemPtr dev_mem;

    size_t total_n = 0, cut_num = 0;
    int griddim = 0;
    for (int i = 0; i < n_read; i++) {
        total_n += reads[i].n;
        int an_p_block = range_kernel_config.anchor_per_block;
        int an_p_cut = range_kernel_config.blockdim;
        int block_num = (reads[i].n - 1) / an_p_block + 1;
        griddim += block_num;
        cut_num += (reads[i].n - 1) / an_p_cut + 1;
    }

    plmem_malloc_host_mem(&host_mem, total_n, griddim);
    plmem_malloc_device_mem(&dev_mem, total_n, griddim, cut_num);
    plmem_reorg_input_arr(reads, n_read, &host_mem, range_kernel_config);
    // sanity check
    assert(host_mem.griddim == griddim);
    assert(host_mem.cut_num == cut_num);
    assert(host_mem.total_n == total_n);

    plmem_sync_h2d_memcpy(&host_mem, &dev_mem);
    plrange_sync_range_selection(&dev_mem, misc
#ifdef DEBUG_CHECK
    , reads
#endif
    );
    // plscore_sync_long_short_forward_dp(&dev_mem, misc);
    plscore_sync_naive_forward_dp(&dev_mem, misc);
    plmem_sync_d2h_memcpy(&host_mem, &dev_mem);

    plchain_backtracking(&host_mem, reads, misc, km);

    plmem_free_host_mem(&host_mem);
    plmem_free_device_mem(&dev_mem);
}

/**
 * Wait and find a free stream
 *  stream_setup:    [in]
 *  batchid:         [in]
 *  stream_id:      [out] stream_id to schedule to
 * RETURN
 *  true if need to cleanup current stream
*/
int plchain_schedule_stream(const streamSetup_t stream_setup, const int batchid){
    /* Haven't fill all the streams*/
    if (batchid < stream_setup.num_stream) {
        return batchid;
    }
    
    // wait until one stream is free
    int streamid = -1;
    while(streamid == -1){
        for (int t = 0; t < stream_setup.num_stream; t++){
            if (!cudaEventQuery(stream_setup.streams[t].cudaevent)) {
                streamid = t;
                // FIXME: unnecessary recreate?
                cudaEventDestroy(stream_setup.streams[t].cudaevent);
                cudaEventCreate(&stream_setup.streams[t].cudaevent);
                cudaCheck();
                break;
            }
            // cudaCheck();
        }
    }
    return streamid;
}

void plchain_cal_score_launch(chain_read_t **reads_, int *n_read_, Misc misc, streamSetup_t stream_setup, int batchid, void* km){
    chain_read_t* reads = *reads_;
    *reads_ = NULL;
    int n_read = *n_read_;
    *n_read_ = 0;
    /* stream scheduler */
    int stream_id = plchain_schedule_stream(stream_setup, batchid);
    if (stream_setup.streams[stream_id].busy) {
        // cleanup previous batch in the stream
        plchain_backtracking(&stream_setup.streams[stream_id].host_mem,
                             stream_setup.streams[stream_id].reads, misc, km);
        *reads_ = stream_setup.streams[stream_id].reads;
        *n_read_ = stream_setup.streams[stream_id].host_mem.size;
        stream_setup.streams[stream_id].busy = false;
    }

    // size sanity check
    size_t total_n = 0, cut_num = 0;
    int griddim = 0;
    for (int i = 0; i < n_read; i++) {
        total_n += reads[i].n;
        int an_p_block = range_kernel_config.anchor_per_block;
        int an_p_cut = range_kernel_config.blockdim;
        int block_num = (reads[i].n - 1) / an_p_block + 1;
        griddim += block_num;
        cut_num += (reads[i].n - 1) / an_p_cut + 1;
    }
    if (stream_setup.max_anchors_stream < total_n){
        fprintf(stderr, "max_anchors_stream %lu total_n %lu n_read %d\n",
                stream_setup.max_anchors_stream, total_n, n_read);
    }

    if (stream_setup.max_range_grid < griddim) {
        fprintf(stderr, "max_range_grid %d griddim %d, total_n %lu n_read %d\n",
                stream_setup.max_range_grid, griddim, total_n, n_read);
    }

    assert(stream_setup.max_anchors_stream >= total_n);
    assert(stream_setup.max_range_grid >= griddim);
    assert(stream_setup.max_num_cut >= cut_num);

    plmem_reorg_input_arr(reads, n_read,
                          &stream_setup.streams[stream_id].host_mem,
                          range_kernel_config);

    plmem_async_h2d_memcpy(&stream_setup.streams[stream_id]);
    plrange_async_range_selection(&stream_setup.streams[stream_id].dev_mem,
                                  &stream_setup.streams[stream_id].cudastream);
    // plscore_async_naive_forward_dp(&stream_setup.streams[stream_id].dev_mem,
    //                                &stream_setup.streams[stream_id].cudastream);
    plscore_async_long_short_forward_dp(&stream_setup.streams[stream_id].dev_mem,
                                   &stream_setup.streams[stream_id].cudastream);
    plmem_async_d2h_memcpy(&stream_setup.streams[stream_id]);
    cudaEventRecord(stream_setup.streams[stream_id].cudaevent,
                    stream_setup.streams[stream_id].cudastream);
    stream_setup.streams[stream_id].busy = true;
    stream_setup.streams[stream_id].reads = reads;
    cudaCheck();
}


void plchain_cal_score_async(chain_read_t **reads_, int *n_read_, Misc misc, streamSetup_t stream_setup, int thread_id, void* km){
    chain_read_t* reads = *reads_;
    *reads_ = NULL;
    int n_read = *n_read_;
    *n_read_ = 0;
    /* sync stream and process previous batch */
    int stream_id = thread_id;
    if (stream_setup.streams[stream_id].busy) {
        cudaStreamSynchronize(stream_setup.streams[stream_id].cudastream);
#ifdef DEBUG_CHECK
        size_t total_n = stream_setup.streams[stream_id].host_mem.total_n;
        chain_read_t* reads = stream_setup.streams[stream_id].reads;
        deviceMemPtr* dev_mem = &stream_setup.streams[stream_id].dev_mem;
        size_t cut_num = stream_setup.streams[stream_id].host_mem.cut_num;
        // check range
        int32_t* range = (int32_t*)malloc(sizeof(int32_t) * total_n);
        cudaMemcpy(range, dev_mem->d_range, sizeof(int32_t) * total_n,
                   cudaMemcpyDeviceToHost);

        size_t* cut = (size_t*)malloc(sizeof(size_t) * cut_num);
        cudaMemcpy(cut, dev_mem->d_cut, sizeof(size_t) * cut_num,
                   cudaMemcpyDeviceToHost);
        for (int readid = 0, cid = 0, idx = 0; readid < dev_mem->size;
             readid++) {
#ifdef DEBUG_VERBOSE
            debug_print_cut(cut + cid, cut_num - cid, reads[readid].n, idx);
#endif
            cid += debug_check_cut(cut + cid, range, cut_num - cid,
                                   reads[readid].n, idx);
            idx += reads[readid].n;
        }
        int64_t read_start = 0;
        for (int i = 0; i < dev_mem->size; i++) {
#ifdef DEBUG_VERBOSE
            debug_print_successor_range(range + read_start, reads[i].n);
#endif
            // debug_check_range(range + read_start, input_arr[i].range,
            // input_arr[i].n);
            read_start += reads[i].n;
        }
        free(range);
        free(cut);
#endif
        // cleanup previous batch in the stream
        plchain_backtracking(&stream_setup.streams[stream_id].host_mem,
                                stream_setup.streams[stream_id].reads, misc, km);
        *reads_ = stream_setup.streams[stream_id].reads;
        *n_read_ = stream_setup.streams[stream_id].host_mem.size;
        stream_setup.streams[stream_id].busy = false;
    }

    // size sanity check
    size_t total_n = 0, cut_num = 0;
    int griddim = 0;
    for (int i = 0; i < n_read; i++) {
        total_n += reads[i].n;
        int an_p_block = range_kernel_config.anchor_per_block;
        int an_p_cut = range_kernel_config.blockdim;
        int block_num = (reads[i].n - 1) / an_p_block + 1;
        griddim += block_num;
        cut_num += (reads[i].n - 1) / an_p_cut + 1;
    }
    if (stream_setup.max_anchors_stream < total_n){
        fprintf(stderr, "max_anchors_stream %lu total_n %lu n_read %d\n",
                stream_setup.max_anchors_stream, total_n, n_read);
    }

    if (stream_setup.max_range_grid < griddim) {
        fprintf(stderr, "max_range_grid %d griddim %d, total_n %lu n_read %d\n",
                stream_setup.max_range_grid, griddim, total_n, n_read);
    }

    assert(stream_setup.max_anchors_stream >= total_n);
    assert(stream_setup.max_range_grid >= griddim);
    assert(stream_setup.max_num_cut >= cut_num);

    plmem_reorg_input_arr(reads, n_read,
                          &stream_setup.streams[stream_id].host_mem,
                          range_kernel_config);

    plmem_async_h2d_memcpy(&stream_setup.streams[stream_id]);
    plrange_async_range_selection(&stream_setup.streams[stream_id].dev_mem,
                                  &stream_setup.streams[stream_id].cudastream);
    // plscore_async_naive_forward_dp(&stream_setup.streams[stream_id].dev_mem,
    //                                &stream_setup.streams[stream_id].cudastream);
    plscore_async_long_short_forward_dp(&stream_setup.streams[stream_id].dev_mem,
                                   &stream_setup.streams[stream_id].cudastream);
    plmem_async_d2h_memcpy(&stream_setup.streams[stream_id]);
    cudaEventRecord(stream_setup.streams[stream_id].cudaevent,
                    stream_setup.streams[stream_id].cudastream);
    stream_setup.streams[stream_id].busy = true;
    stream_setup.streams[stream_id].reads = reads;
    cudaCheck();
}

#ifdef __cplusplus
extern "C" {
#endif  // __cplusplus

void init_blocking_gpu(size_t* total_n, int* max_reads, int *min_n, Misc misc) {
    plmem_initialize(total_n, max_reads, min_n);
}

void init_stream_gpu(size_t* total_n, int* max_reads, int *min_n, Misc misc) {
    fprintf(stderr, "[M::%s] gpu initialized for chaining\n", __func__);
    plmem_stream_initialize(total_n, max_reads, min_n);
    plrange_upload_misc(misc);
    plscore_upload_misc(misc);
}

/**
 * worker for forward chaining on cpu (blocking)
 * use KMALLOC and kfree for cpu memory management
 */
void chain_blocking_gpu(const mm_idx_t *mi, const mm_mapopt_t *opt, chain_read_t *in_arr, int n_read, void* km) {
    // assume only one seg. and qlen_sum desn't matter
    assert(opt->max_frag_len <= 0);
    assert(!!(opt->flag & MM_F_SR));
    assert(in_arr[0].n_seg == 1);
    Misc misc = build_misc(mi, opt, 0, 1);
    plchain_cal_score_sync(in_arr, n_read, misc, km);
    for (int i = 0; i < n_read; i++){
        post_chaining_helper(mi, opt,  &in_arr[i], misc, km);
    }
}

/**
 * worker for launching forward chaining on gpu (streaming)
 * use KMALLOC and kfree for cpu memory management
 * [in/out] in_arr_: ptr to array of reads, updated to a batch launched in previous run 
 *                  (NULL if no finishing batch)
 * [in/out] n_read_: ptr to num of reads in array, updated to a batch launched in previous run
 *                  (NULL if no finishing batch)
*/

// void chain_stream_gpu(const input_meta_t* meta, chain_read_t**in_arr_, int *n_read_) {
//     static int batchid = 0;
//     Misc misc = build_misc(INT64_MAX);
//     plchain_cal_score_launch(in_arr_, n_read_, misc, stream_setup, batchid);
//     batchid++;
//     if (in_arr_){
//         int n_read = *n_read_;
//         chain_read_t* out_arr = *in_arr_;
//         for (int i = 0; i < n_read; i++){
//             post_chaining_helper(&out_arr[i], meta->refs, &out_arr[i].seq,
//                                  misc.max_dist_x, out_arr[i].km);
//         }
//     }
// }

void chain_stream_gpu(const mm_idx_t *mi, const mm_mapopt_t *opt, chain_read_t **in_arr_, int *n_read_,
                      int thread_id, void* km) {
    // assume only one seg. and qlen_sum desn't matter
    assert(opt->max_frag_len <= 0);
    assert(!!(opt->flag & MM_F_SR));
    assert(in_arr[0].n_seg == 1);
    Misc misc = build_misc(mi, opt, 0, 1);
    plchain_cal_score_async(in_arr_, n_read_, misc, stream_setup, thread_id, km);
    if (in_arr_) {
        int n_read = *n_read_;
        chain_read_t* out_arr = *in_arr_;
        for (int i = 0; i < n_read; i++) {
            post_chaining_helper(mi, opt, &out_arr[i], misc, km);
        }
    }
}

// /** 
//  * worker for finish all forward chaining kernenls on gpu
//  * use KMALLOC and kfree for cpu memory management
//  * [out] batches:   array of batches
//  * [out] num_reads: array of number of reads in each batch
//  */
// void finish_stream_gpu(const input_meta_t* meta, chain_read_t** batches,
//                        int* num_reads, int num_batch) {
//     int batch_handled = 0;
    
//     /* Sanity Check */
// #ifdef DEBUG_VERBOSE
//     fprintf(stderr, "[Info] Finishing up %d pending batches\n", num_batch);
// #endif
// #ifdef DEBUG_CHECK
//         for (int t = 0; t < stream_setup.num_stream; t++) {
//         if (stream_setup.streams[t].busy) batch_handled++;
//     }
//     if (batch_handled != num_batch){
//         fprintf(stderr, "[Error] busy streams %d num_batch %d\n", batch_handled,
//                 num_batch);
//         exit(1);
//     }
//     batch_handled = 0;
// #endif  // DEBUG_CHECK

//     Misc misc = build_misc(INT64_MAX);
//     /* Sync all the pending batches + backtracking */
//     while(batch_handled < num_batch){
//         for (int t = 0; t < stream_setup.num_stream; t++) {
//             if (!stream_setup.streams[t].busy) continue;
//             if (!cudaEventQuery(stream_setup.streams[t].cudaevent)) {
//                 cudaCheck();
//                 assert(batch_handled < num_batch);
//                 plchain_backtracking(&stream_setup.streams[t].host_mem,
//                                         stream_setup.streams[t].reads, misc);
//                 batches[batch_handled] = stream_setup.streams[t].reads;
//                 num_reads[batch_handled] = stream_setup.streams[t].host_mem.size;
//                 for (int i = 0; i < num_reads[batch_handled]; i++) {
//                     post_chaining_helper(&batches[batch_handled][i], meta->refs,
//                                         &batches[batch_handled][i].seq,
//                                         misc.max_dist_x,
//                                         batches[batch_handled][i].km);
//                 }
//                 batch_handled++;
//                 stream_setup.streams[t].busy = false;
//             }
//         }
//     }

//     assert(num_batch == batch_handled);
//     for (int t = 0; t < stream_setup.num_stream; t++){
//         plmem_free_host_mem(&stream_setup.streams[t].host_mem);
//         plmem_free_device_mem(&stream_setup.streams[t].dev_mem);
//     }
// }

/**
 * worker for finish all forward chaining kernenls on gpu
 * use KMALLOC and kfree for cpu memory management
 * [out] batches:   array of batches
 * [out] num_reads: array of number of reads in each batch
 */
void finish_stream_gpu(const mm_idx_t *mi, const mm_mapopt_t *opt, chain_read_t** reads_,
                       int* n_read_, int t, void* km) {
    // assume only one seg. and qlen_sum desn't matter
    assert(opt->max_frag_len <= 0);
    assert(!!(opt->flag & MM_F_SR));
    assert(in_arr[0].n_seg == 1);
    Misc misc = build_misc(mi, opt, 0, 1);
    /* Sync all the pending batches + backtracking */
    if (!stream_setup.streams[t].busy) {
        *reads_ = NULL;
        *n_read_ = 0;
        return;
    }

    chain_read_t* reads;
    int n_read;
    cudaStreamSynchronize(stream_setup.streams[t].cudastream);
    cudaCheck();
    plchain_backtracking(&stream_setup.streams[t].host_mem,
                         stream_setup.streams[t].reads, misc, km);
    reads = stream_setup.streams[t].reads;
    n_read = stream_setup.streams[t].host_mem.size;
    for (int i = 0; i < n_read; i++) {
        post_chaining_helper(mi, opt, &reads[i], misc, km);
    }
    stream_setup.streams[t].busy = false;

    *reads_ = reads;
    *n_read_ = n_read;

}


void free_stream_gpu(int n_threads){
    for (int t = 0; t < n_threads; t++){
        plmem_free_host_mem(&stream_setup.streams[t].host_mem);
        plmem_free_device_mem(&stream_setup.streams[t].dev_mem);
    }
    fprintf(stderr, "[M::%s] gpu free memory\n", __func__);
}

#ifdef __cplusplus
} // extern "C"
#endif  // __cplusplus