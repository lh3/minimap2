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
#include <utility>
#include <algorithm>

#ifdef DEBUG_CHECK
#include "planalyze.cuh"
#include "debug.h"
#endif // DEBUG_CHECK

// utils functions
struct
{
    bool operator()(std::pair<size_t,unsigned> a, std::pair<size_t,unsigned> b) const {
        return a.first > b.first;
    }
}
comp;

void pairsort(seg_t *sdata, unsigned *map, unsigned N){
    std::pair<size_t,unsigned> elements[N];
    for (unsigned i = 0; i < N; ++i){
      elements[i].second = map[i];
      elements[i].first = sdata[i].end_idx - sdata[i].start_idx;
    }

    std::sort(elements, elements+N, comp); // descent order

    for (unsigned i = 0; i < N; ++i){
      map[i] = elements[i].second;
    }
}


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
// print scores
#if defined(DEBUG_VERBOSE) && 0
        debug_print_score(p, f, reads[i].n);
#endif
// Check score w.r.t to input (MAKE SURE INPUT SCORE EXISTS: search for SCORE CHECK) 
#if defined(DEBUG_CHECK) && 0
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


//////////////////////////////////////////////////////////////////////////
///////////         Stream Management    /////////////////////////////////
//////////////////////////////////////////////////////////////////////////

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
            if (!cudaEventQuery(stream_setup.streams[t].stopevent)) {
                streamid = t;
                // FIXME: unnecessary recreate?
                cudaEventDestroy(stream_setup.streams[t].stopevent);
                cudaEventCreate(&stream_setup.streams[t].stopevent);
                cudaCheck();
                break;
            }
            // cudaCheck();
        }
    }
    return streamid;
}


// Global variable for debug prints. Throughput, runtime & mem usage
#ifdef DEBUG_PRINT
    float kernel_mem_usage[MAX_MICRO_BATCH + 1] = {0};
    float kernel_throughput[MAX_MICRO_BATCH + 1] = {0};
#endif // DEBUG_PRINT

/*
 * Accepts a stream that has already been synced, and finished processing a batch
 * Finish and cleanup the stream, save primary chain results to unpinned CPU memory.  
 * RETURN: number of reads in last batch 
*/
int plchain_post_gpu_helper(streamSetup_t stream_setup, int stream_id, Misc misc, void* km){
    int n_reads = 0; // Number of reads in the batch

#if defined(DEBUG_CHECK)
    planalyze_long_kernel(stream_setup.streams[stream_id], kernel_throughput);
#endif  // DEBUG_CHECK

#ifdef DEBUG_PRINT
    kernel_mem_usage[score_kernel_config.micro_batch] = (float)(*stream_setup.streams[stream_id].long_mem.total_long_segs_n)/stream_setup.long_seg_buffer_size_stream*100;

    float kernel_runtime_ms[MAX_MICRO_BATCH + 1] = {0};
    float kernel_throughput_anchors[MAX_MICRO_BATCH + 1] = {0};
    cudaEventElapsedTime(&kernel_runtime_ms[score_kernel_config.micro_batch], stream_setup.streams[stream_id].long_kernel_event, stream_setup.streams[stream_id].stopevent);
    kernel_throughput_anchors[score_kernel_config.micro_batch] = *stream_setup.streams[stream_id].long_mem.total_long_segs_n / kernel_runtime_ms[score_kernel_config.micro_batch] / (float)1000;
#endif


    seg_t* long_segs = stream_setup.streams[stream_id].long_mem.long_segs_og_idx;
    size_t long_seg_idx = 0;
    size_t long_i = 0;
    for (int uid = 0; uid < score_kernel_config.micro_batch; uid++) {
        if (stream_setup.streams[stream_id].host_mems[uid].size == 0) continue;
        // regorg long to each host mem ptr
        // NOTE: this is the number of long segs till this microbatch
        unsigned int long_segs_num = stream_setup.streams[stream_id].host_mems[uid].long_segs_num[0];
#ifdef DEBUG_VERBOSE
    fprintf(stderr, "[Debug] %s (%s:%d) MICROBATCH$%d FINISHED: long seg %lu - %u", 
                __func__, __FILE__, __LINE__, uid, long_seg_idx, long_segs_num);
#endif // DEBUG_VERBOSE
        size_t total_n_long_segs = 0;
        for (; long_seg_idx < long_segs_num; long_seg_idx++) {
            for (size_t i = long_segs[long_seg_idx].start_idx;
                 i < long_segs[long_seg_idx].end_idx; i++, long_i++) {
                stream_setup.streams[stream_id].host_mems[uid].f[i] =
                    stream_setup.streams[stream_id].long_mem.f_long[long_i];
                stream_setup.streams[stream_id].host_mems[uid].p[i] =
                    stream_setup.streams[stream_id].long_mem.p_long[long_i];
            }
            total_n_long_segs += long_segs[long_seg_idx].end_idx - long_segs[long_seg_idx].start_idx;
        }

        // backtrack after p/f is copied
        plchain_backtracking(&stream_setup.streams[stream_id].host_mems[uid],
                            stream_setup.streams[stream_id].reads + n_reads, misc, km);
        // accumulate n_reads
        n_reads += stream_setup.streams[stream_id].host_mems[uid].size;
#ifdef DEBUG_PRINT
        cudaEventElapsedTime(&kernel_runtime_ms[uid], stream_setup.streams[stream_id].short_kernel_start_event[uid], stream_setup.streams[stream_id].short_kernel_stop_event[uid]);
        kernel_throughput_anchors[uid] =
            (stream_setup.streams[stream_id].host_mems[uid].total_n - total_n_long_segs) /
            kernel_runtime_ms[uid] / (float)1000;
#ifdef DEBUG_VERBOSE
        fprintf(stderr, ", %.2f%% anchors are in long segs. \n", (float)total_n_long_segs /  stream_setup.streams[stream_id].host_mems[uid].total_n * 100);
#endif // DEBUG_VERBOSE
#endif // DEBUG_PRINT
    }

#ifdef DEBUG_PRINT
    fprintf(stderr, "----------------------------------------------------------------------------\n                  ");
    for (int uid = 0; uid < score_kernel_config.micro_batch; uid++) fprintf(stderr, "Short%d     ", uid);
    fprintf(stderr, "Long\n");
    fprintf(stderr, "----------------------------------------------------------------------------\n");
    fprintf(stderr, "Mem Usage  = ");
    for (int uid = 0; uid < score_kernel_config.micro_batch; uid++) fprintf(stderr, " %9.2f%%", kernel_mem_usage[uid]);
    fprintf(stderr, " %9.2f %%\n", kernel_mem_usage[score_kernel_config.micro_batch]);
    fprintf(stderr, "Runtime(s) = ");
    for (int uid = 0; uid < score_kernel_config.micro_batch; uid++) fprintf(stderr, "%11.2f", kernel_runtime_ms[uid] / 1000);
    fprintf(stderr, " %11.2f\n", kernel_runtime_ms[score_kernel_config.micro_batch] / 1000);
    fprintf(stderr, "BW (Ma/s)  = ");
    for (int uid = 0; uid < score_kernel_config.micro_batch; uid++) fprintf(stderr, "%11.2f", kernel_throughput_anchors[uid]);
    fprintf(stderr, " %11.2f\n", kernel_throughput_anchors[score_kernel_config.micro_batch]);
    fprintf(stderr, "BW(Mpair/s)= ");
        for (int uid = 0; uid < score_kernel_config.micro_batch; uid++) fprintf(stderr, "%11.2f", kernel_throughput[uid]);
    fprintf(stderr, " %11.2f\n", kernel_throughput[score_kernel_config.micro_batch]);
    fprintf(stderr, "----------------------------------------------------------------------------\n");
    if (kernel_mem_usage[score_kernel_config.micro_batch] > 99){
        fprintf(stderr,
                "[WARNING] long segment buffer is full. Consider increase "
                "long_seg_buffer_size to improve performance.\n");
    }
#endif 

    return n_reads;
} 


/* 
 * 1. synchronize stream and process previous batch. cleanup stream
 * 2. launch kernels (asynchornizely) for the input batch 
 */

void plchain_cal_score_async(chain_read_t **reads_, int *n_read_, Misc misc, streamSetup_t stream_setup, int thread_id, void* km){
    chain_read_t* reads = *reads_;
    *reads_ = NULL;
    int n_read = *n_read_;
    *n_read_ = 0;

    /* sync stream and process previous batch */
    int stream_id = thread_id;
    if (stream_setup.streams[stream_id].busy) {
        cudaStreamSynchronize(stream_setup.streams[stream_id].cudastream);
        *n_read_ = plchain_post_gpu_helper(stream_setup, stream_id, misc, km);
        *reads_ = stream_setup.streams[stream_id].reads;
        stream_setup.streams[stream_id].busy = false;
    }

#ifdef DEBUG_PRINT
    for(int uid = 0; uid < score_kernel_config.micro_batch + 1; uid++) {
        kernel_mem_usage[uid] = 0;
        kernel_throughput[uid] = 0;
    }
#endif // DEBUG_PRINT


    cudaEventRecord(stream_setup.streams[stream_id].startevent,
                    stream_setup.streams[stream_id].cudastream);
    size_t total_n = 0;
    for (int i = 0; i < n_read; i++) {
        total_n += reads[i].n;
    } // compute total_n first
#ifdef DEBUG_PRINT
    fprintf(stderr, "[Info] %s (%s:%d) Launching Batch: n_read %d, total anchors %lu (mem usage: %.2f%%)\n", __func__, __FILE__, __LINE__, n_read, total_n, (float)total_n/stream_setup.max_anchors_stream*100);
#endif // DEBUG_PRINT

    // reset long seg counters
    cudaMemsetAsync(stream_setup.streams[stream_id].dev_mem.d_long_seg_count, 0, sizeof(unsigned int),
                    stream_setup.streams[stream_id].cudastream);
    cudaMemsetAsync(stream_setup.streams[stream_id].dev_mem.d_total_n_long, 0, sizeof(size_t),
                    stream_setup.streams[stream_id].cudastream);
    stream_setup.streams[stream_id].long_mem.total_long_segs_num[0] = 0;
    stream_setup.streams[stream_id].long_mem.total_long_segs_n[0] = 0;
    for(int uid = 0; uid < score_kernel_config.micro_batch; uid++) {
        stream_setup.streams[stream_id].host_mems[uid].long_segs_num[0] = 0;
        stream_setup.streams[stream_id].host_mems[uid].index = uid;
        stream_setup.streams[stream_id].host_mems[uid].griddim = 0;
        stream_setup.streams[stream_id].host_mems[uid].size = 0;
        stream_setup.streams[stream_id].host_mems[uid].total_n = 0;
        stream_setup.streams[stream_id].host_mems[uid].cut_num = 0;
    }

    stream_setup.streams[stream_id].reads = reads;
    stream_setup.streams[stream_id].n_read = n_read;
    int read_start = 0;

    for (int uid = 0; uid < score_kernel_config.micro_batch; uid++) {
        if (read_start == n_read) continue;
#ifdef USEHIP
        roctxRangePushA("microbatch");
#endif
        // decide the size of micro batch
        size_t batch_n = 0;
        int read_end = 0;
        size_t cut_num = 0;
        int griddim = 0;
        for (read_end = read_start; read_end < n_read; read_end++) {
            if (batch_n + reads[read_end].n > stream_setup.max_anchors_stream) {
                break;
            }
            batch_n += reads[read_end].n;
            int an_p_block = range_kernel_config.anchor_per_block;
            int an_p_cut = range_kernel_config.blockdim;
            int block_num = (reads[read_end].n - 1) / an_p_block + 1;
            griddim += block_num;
            cut_num += (reads[read_end].n - 1) / an_p_cut + 1;
        }
#ifdef DEBUG_VERBOSE
        fprintf(stderr, "[Debug] %s (%s:%d) MICROBATCH#%d batch_n %lu, read_start %d, read_end %d usage %.2f %%\n", __func__, __FILE__, __LINE__, uid, batch_n, read_start, read_end, (float)batch_n/stream_setup.max_anchors_stream*100);
#endif // DEBUG_VERBOS 

#ifdef DEBUG_PRINT
        kernel_mem_usage[uid] = (float)batch_n/stream_setup.max_anchors_stream*100;
#endif // DEBUG_PRINT

        // sanity check
        assert(stream_setup.max_anchors_stream >= batch_n);
        assert(stream_setup.max_range_grid >= griddim);
        assert(stream_setup.max_num_cut >= cut_num);
        // work on micro batch
#ifdef USEHIP
        roctxRangePushA("reorg");
#endif
        // step1: reorg input
        plmem_reorg_input_arr(reads + read_start, read_end - read_start,
                          &stream_setup.streams[stream_id].host_mems[uid],
                          range_kernel_config);
        // step2: copy to device
#ifdef USEHIP
        roctxRangePop();
#endif
        plmem_async_h2d_short_memcpy(&stream_setup.streams[stream_id], uid);
        // step3: range selection
#ifdef DEBUG_PRINT
        cudaEventRecord(stream_setup.streams[stream_id].short_kernel_start_event[uid],
                    stream_setup.streams[stream_id].cudastream);
#endif // DEBUG_PRINT
        plrange_async_range_selection(&stream_setup.streams[stream_id].dev_mem,
                                    &stream_setup.streams[stream_id].cudastream);
        // step4: score generation for short and mid segs
        plscore_async_short_mid_forward_dp(&stream_setup.streams[stream_id].dev_mem,
                                    &stream_setup.streams[stream_id].cudastream);
#ifdef DEBUG_PRINT
        cudaEventRecord(stream_setup.streams[stream_id].short_kernel_stop_event[uid],
                    stream_setup.streams[stream_id].cudastream);
#endif // DEBUG_PRINT
        // step5: copy short and mid results back
        plmem_async_d2h_short_memcpy(&stream_setup.streams[stream_id], uid);
        // update index
        read_start = read_end;

#ifdef USEHIP
        roctxRangePop();
#endif

#ifdef DEBUG_CHECK
        planalyze_short_kernel(stream_setup.streams[stream_id], uid, kernel_throughput);
#endif // DEBUG_CHECK
    }

// FIXME: temporary solution for microbatching
    if (read_start < n_read) {
        fprintf(stderr, "[WARNING] Unable to fit reads %d - %d into a microbatch. Fall back to cpu chaining\n", read_start, n_read-1);
    }

    // step6: copy back long_segs_og
    cudaStreamSynchronize(stream_setup.streams[stream_id].cudastream);
    unsigned int num_long_seg;
    cudaMemcpy(&num_long_seg, stream_setup.streams[stream_id].dev_mem.d_long_seg_count, sizeof(unsigned int),
                cudaMemcpyDeviceToHost);
    seg_t* long_segs_og = (seg_t*)malloc(sizeof(seg_t) * num_long_seg);
    cudaMemcpy(long_segs_og, stream_setup.streams[stream_id].dev_mem.d_long_seg_og, sizeof(seg_t) * num_long_seg,
                cudaMemcpyDeviceToHost);

    // step7: sort long segs in descent order
    unsigned *map = new unsigned[num_long_seg];
    for (unsigned i = 0; i < num_long_seg; i++) {
        map[i] = i;
    }
    pairsort(long_segs_og, map, num_long_seg);
    #ifdef DEBUG_VERBOSE
    auto last_length = long_segs_og[map[0]].end_idx - long_segs_og[map[0]].start_idx;
    for (int i = 1; i < num_long_seg; i++){
        auto this_length = long_segs_og[map[i]].end_idx - long_segs_og[map[i]].start_idx;
        if (this_length > last_length)
            fprintf(stderr, "Failed sort at: %d - %u\n", i, map[i]);
    }
    #endif // DEBUG_VERBOSE
    free(long_segs_og);

    // step8: copy map to device
    cudaMalloc(&stream_setup.streams[stream_id].dev_mem.d_map, sizeof(unsigned) * num_long_seg);
    cudaMemcpy(stream_setup.streams[stream_id].dev_mem.d_map, map, sizeof(unsigned) * num_long_seg, cudaMemcpyHostToDevice);
    free(map);

    cudaEventRecord(stream_setup.streams[stream_id].long_kernel_event,
                    stream_setup.streams[stream_id].cudastream);
    plscore_async_long_forward_dp(&stream_setup.streams[stream_id].dev_mem,
                                   &stream_setup.streams[stream_id].cudastream);
    cudaEventRecord(stream_setup.streams[stream_id].stopevent,
                    stream_setup.streams[stream_id].cudastream);
    plmem_async_d2h_long_memcpy(&stream_setup.streams[stream_id]);
    stream_setup.streams[stream_id].busy = true;
    cudaCheck();
}

#ifdef __cplusplus
extern "C" {
#endif  // __cplusplus

void init_stream_gpu(size_t* total_n, int* max_reads, int *min_n, char gpu_config_file[], Misc misc) {
    plmem_stream_initialize(total_n, max_reads, min_n, gpu_config_file);
    plrange_upload_misc(misc);
    plscore_upload_misc(misc);
#ifdef DEBUG_PRINT
    fprintf(stderr, "[Info::%s] gpu initialized for chaining with config %s\n", __func__, gpu_config_file);
    fprintf(stderr, "[Info::%s] Compile time config: \n", __func__);
#ifdef USEHIP
    fprintf(stderr, "\t\t USE HIP\n");
#else
    fprintf(stderr, "\t\t USE CUDA\n");
#endif // USEHIP
#ifdef MAX_MICRO_BATCH
    fprintf(stderr, "\t\t MAX MICRO BATCH       \t%d\n", MAX_MICRO_BATCH);
#endif // MAX_MICRO_BATCH
#endif  // DEBUG_PRINT
}

/**
 * worker for launching forward chaining on gpu (streaming)
 * use KMALLOC and kfree for cpu memory management
 * [in/out] in_arr_: ptr to array of reads, updated to a batch launched in previous run 
 *                  (NULL if no finishing batch)
 * [in/out] n_read_: ptr to num of reads in array, updated to a batch launched in previous run
 *                  (NULL if no finishing batch)
*/
void chain_stream_gpu(const mm_idx_t *mi, const mm_mapopt_t *opt, chain_read_t **in_arr_, int *n_read_,
                      int thread_id, void* km) {
    // assume only one seg. and qlen_sum desn't matter
    assert(opt->max_frag_len <= 0);
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
    Misc misc = build_misc(mi, opt, 0, 1);
    /* Sync all the pending batches + backtracking */
    if (!stream_setup.streams[t].busy) {
        *reads_ = NULL;
        *n_read_ = 0;
        return;
    }

    chain_read_t* reads;
    int n_read = 0;
    cudaStreamSynchronize(stream_setup.streams[t].cudastream);
    cudaCheck();

    n_read = plchain_post_gpu_helper(stream_setup, t, misc, km);
    reads = stream_setup.streams[t].reads;
    stream_setup.streams[t].busy = false;

    for (int i = 0; i < n_read; i++) {
        post_chaining_helper(mi, opt, &reads[i], misc, km);
    }

    *reads_ = reads;
    *n_read_ = n_read;

}


void free_stream_gpu(int n_threads){
    plmem_stream_cleanup();

    size_t gpu_free_mem, gpu_total_mem;
    cudaMemGetInfo(&gpu_free_mem, &gpu_total_mem);
#ifdef DEBUG_PRINT
        fprintf(stderr, "[Info] GPU free mem: %f GB, total mem: %f GB (after cleanup) \n", (float)gpu_free_mem / OneG, (float)gpu_total_mem / OneG);
#endif
}

#ifdef __cplusplus
} // extern "C"
#endif  // __cplusplus