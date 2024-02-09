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

#ifdef DEBUG_CHECK
#include "debug.h"
#endif // DEBUG_CHECK

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
    // FIXME: DISABLED BACKTRACK, REMOVE THE RETURN HERE
    return;
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

    plmem_malloc_host_mem(&host_mem, total_n, griddim, 0);
    plmem_malloc_device_mem(&dev_mem, total_n, griddim, cut_num);
    plmem_reorg_input_arr(reads, n_read, &host_mem, range_kernel_config);
    // sanity check
    assert(host_mem.griddim == griddim);
    assert(host_mem.cut_num == cut_num);
    assert(host_mem.total_n == total_n);

    plmem_sync_h2d_memcpy(&host_mem, &dev_mem);
    plrange_sync_range_selection(&dev_mem, misc);
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

void plchain_cal_score_launch(chain_read_t **reads_, int *n_read_, Misc misc, streamSetup_t stream_setup, int batchid, void* km){
    chain_read_t* reads = *reads_;
    *reads_ = NULL;
    int n_read = *n_read_;
    *n_read_ = 0;
    /* stream scheduler */
    int stream_id = plchain_schedule_stream(stream_setup, batchid);
    if (stream_setup.streams[stream_id].busy) {
#ifdef DEBUG_PRINT
    fprintf(stderr, "[Info] %s (%s:%d) stream %d sync, total_n %lu\n", __func__, __FILE__, __LINE__, stream_id, stream_setup.streams[stream_id].host_mems[0].total_n);
#endif // DEBUG_PRINT 
        // cleanup previous batch in the stream
        plchain_backtracking(&stream_setup.streams[stream_id].host_mems[0],
                             stream_setup.streams[stream_id].reads, misc, km);
        
        *reads_ = stream_setup.streams[stream_id].reads;
        *n_read_ = stream_setup.streams[stream_id].host_mems[0].size;
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
                          &stream_setup.streams[stream_id].host_mems[0],
                          range_kernel_config);

    plmem_async_h2d_memcpy(&stream_setup.streams[stream_id]);
    plrange_async_range_selection(&stream_setup.streams[stream_id].dev_mem,
                                  &stream_setup.streams[stream_id].cudastream);
    // plscore_async_naive_forward_dp(&stream_setup.streams[stream_id].dev_mem,
    //                                &stream_setup.streams[stream_id].cudastream);
    plscore_async_long_short_forward_dp(&stream_setup.streams[stream_id].dev_mem,
                                   &stream_setup.streams[stream_id].cudastream);
    plmem_async_d2h_memcpy(&stream_setup.streams[stream_id]);
    cudaEventRecord(stream_setup.streams[stream_id].stopevent,
                    stream_setup.streams[stream_id].cudastream);
    stream_setup.streams[stream_id].busy = true;
    stream_setup.streams[stream_id].reads = reads;
    cudaCheck();
}
#ifdef DEBUG_VERBOSE

// find long seg range distribution
void plchain_cal_long_seg_range_dis(size_t total_n, deviceMemPtr* dev_mem){
    static uint64_t long_seg_range_dis[5001] = {0};
    static unsigned int long_seg_total = 0;
    static uint64_t long_seg_anchors_total = 0;
    static FILE* fp = NULL;
    // check range
    int32_t* range = (int32_t*)malloc(sizeof(int32_t) * total_n);
    cudaMemcpy(range, dev_mem->d_range, sizeof(int32_t) * total_n,
                cudaMemcpyDeviceToHost);

    unsigned int long_seg_count;
    cudaMemcpy(&long_seg_count, dev_mem->d_long_seg_count, sizeof(unsigned int),
                cudaMemcpyDeviceToHost);
    fprintf(stderr, "[verbose] %u long segs generated\n", long_seg_count);
    seg_t* long_seg = (seg_t*)malloc(sizeof(seg_t) * long_seg_count);
    cudaMemcpy(long_seg, dev_mem->d_long_seg, sizeof(seg_t) * long_seg_count,
                cudaMemcpyDeviceToHost);

    static FILE* fp_all_data = NULL;
    if (!fp_all_data){
        fp_all_data = fopen("long_seg_range_dis_all.csv", "w+");
        fprintf(fp_all_data, "seg_length");
        for (int i = 0; i < 5000; i++) fprintf(fp_all_data, ",%d", i);
        fprintf(fp_all_data, "\n");
    }

    for (int sid = 0; sid < long_seg_count; sid++) {
        uint64_t long_seg_range_dis_all[5001] = {0};
        for (size_t i = long_seg[sid].start_idx; i < long_seg[sid].end_idx; i++){
            assert(range[i] <= 5000);
            long_seg_range_dis[range[i]]++;
            long_seg_range_dis_all[range[i]]++;
        }
        if (long_seg[sid].end_idx - long_seg[sid].start_idx < 530){
            fprintf(stderr, "Found seg of length %lu, segid %lu - %lu \n",
                    long_seg[sid].end_idx - long_seg[sid].start_idx,
                    long_seg[sid].start_segid, long_seg[sid].end_segid);
        }

        fprintf(fp_all_data, "%lu", long_seg[sid].end_idx - long_seg[sid].start_idx);
        for (int i = 0; i < 5000; i++){
            fprintf(fp_all_data, ",%lu", long_seg_range_dis_all[i]);
        }
        fprintf(fp_all_data, "\n");
        long_seg_anchors_total += long_seg[sid].end_idx - long_seg[sid].start_idx;
        ++long_seg_total;
    }   
    if (!fp) {
        fp = fopen("long_seg_range_dis.csv", "w+");
        fprintf(fp, "num_segs,num_anchors");
        for (int i = 0; i < 5000; i++) fprintf(fp, ",%d", i);
        fprintf(fp, "\n");
    }
    fprintf(fp, "%u segs,%lu anchors", long_seg_total, long_seg_anchors_total);
    for (int i = 0; i < 5000; i++){
        fprintf(fp, ",%lu", long_seg_range_dis[i]);
    }
    fprintf(fp, "\n");
    free(range);
    free(long_seg);
}

// range distribution
void plchain_cal_range_dis(size_t total_n, size_t num_cut, deviceMemPtr* dev_mem){
    static uint64_t range_dis[5001] = {0};
    static size_t seg_total = 0;
    static uint64_t anchors_total = 0;
    static FILE* fp = NULL;
    // check range
    int32_t* range = (int32_t*)malloc(sizeof(int32_t) * total_n);
    cudaMemcpy(range, dev_mem->d_range, sizeof(int32_t) * total_n,
                cudaMemcpyDeviceToHost);

    fprintf(stderr, "[verbose] %lu cuts generated\n", num_cut);
    for (size_t i = 0; i < total_n; i++){
        assert(range[i] <= 5000);
        range_dis[range[i]]++;
    }
    anchors_total += total_n;
    seg_total += num_cut;
    if (!fp) {
        fp = fopen("range_dis.csv", "w+");
        fprintf(fp, "num_segs,num_anchors");
        for (int i = 0; i < 5000; i++) fprintf(fp, ",%d", i);
        fprintf(fp, "\n");
    }
    fprintf(fp, "%lusegs,%luanchors", seg_total, anchors_total);
    for (int i = 0; i < 5000; i++){
        fprintf(fp, ",%lu", range_dis[i]);
    }
    fprintf(fp, "\n");
    free(range);
}

// sc pair vs. seg length
void plchain_cal_sc_pair_density(size_t total_n, size_t num_cut, deviceMemPtr* dev_mem){
    // bin width: 10 cuts, max 5000 cuts
    static uint64_t sc_pair_dis[501] = {0}; // number of sc pairs for each seg length
    static uint64_t anchors_dis[501] = {0};
    size_t* cut = (size_t*)malloc(sizeof(size_t) * num_cut);
    cudaMemcpy(cut, dev_mem->d_cut, sizeof(size_t) * num_cut, cudaMemcpyDeviceToHost);

    int32_t* range = (int32_t*)malloc(sizeof(int32_t) * total_n);
    cudaMemcpy(range, dev_mem->d_range, sizeof(int32_t) * total_n, cudaMemcpyDeviceToHost);

    uint64_t start_idx = 0, cut_size = 0;
    for (int cid = 0; cid < num_cut; cid++) {
        if (cut[cid] != SIZE_MAX) {
            uint64_t sc_pair_num = 0;
            for (uint64_t i = start_idx; i < cut[cid]; i++){
                sc_pair_num += range[i];
            }
            if (cut_size/10 < 500) {
                sc_pair_dis[cut_size/10] += sc_pair_num;
                anchors_dis[cut_size/10] += cut[cid] - start_idx;
            } else {
                sc_pair_dis[500] += sc_pair_num;
                anchors_dis[500] += cut[cid] - start_idx;
            }
            cut_size = 0;
            start_idx = cut[cid];
        } else {
            ++cut_size;
        }
    }
    free(range);
    free(cut);

    static FILE* f_sc_pair_dis = NULL;
    if (!f_sc_pair_dis){
        f_sc_pair_dis = fopen("sc_pair_dis.csv", "w+");
        fprintf(f_sc_pair_dis, "seg_len,");
        for (int i = 0; i < 500; i++){
            fprintf(f_sc_pair_dis, "%d,", i*10);
        }
        fprintf(f_sc_pair_dis, "5000\n");
    }
    
    fprintf(f_sc_pair_dis, "sc_pairs,");
    for (int i = 0; i <= 500; i++){
        fprintf(f_sc_pair_dis, "%lu,", sc_pair_dis[i]);
    }
    fprintf(f_sc_pair_dis, "\n");
    fprintf(f_sc_pair_dis, "anchors,");
    for (int i = 0; i <= 500; i++){
        fprintf(f_sc_pair_dis, "%lu,", anchors_dis[i]);
    }
    fprintf(f_sc_pair_dis, "\n");
    fflush(f_sc_pair_dis);
}


#endif // DEBUG_CHECK




#ifdef DEBUG_CHECK

void plchain_debug_analysis(stream_ptr_t stream){
    // TODO: analysis multiple or current host mems
    // TODO: this needs to be recalculated
    size_t uid = 0;
    size_t total_n = stream.host_mems[uid].total_n;
    chain_read_t* reads = stream.reads;
    deviceMemPtr* dev_mem = &stream.dev_mem;
    hostMemPtr* host_mem = &stream.host_mems[uid];
    longMemPtr* long_mem = &stream.long_mem;
    size_t cut_num = stream.host_mems[uid].cut_num;

    unsigned int num_mid_seg, num_long_seg;
    cudaMemcpy(&num_mid_seg, dev_mem->d_mid_seg_count, sizeof(unsigned int),
                cudaMemcpyDeviceToHost);
    cudaMemcpy(&num_long_seg, dev_mem->d_long_seg_count, sizeof(unsigned int),
                cudaMemcpyDeviceToHost);

    fprintf(stderr, "[DEBUG] total segs: %lu, short:%lu mid: %u long: %u\n", cut_num, cut_num - num_mid_seg - num_long_seg, num_mid_seg, num_long_seg);


    int32_t* range = (int32_t*)malloc(sizeof(int32_t) * total_n);
    cudaMemcpy(range, dev_mem->d_range, sizeof(int32_t) * total_n,
                cudaMemcpyDeviceToHost);
    size_t* cut = (size_t*)malloc(sizeof(size_t) * cut_num);
    cudaMemcpy(cut, dev_mem->d_cut, sizeof(size_t) * cut_num,
                cudaMemcpyDeviceToHost);
    seg_t* long_segs = (seg_t*)malloc(sizeof(seg_t) * num_long_seg);
    cudaMemcpy(long_segs, dev_mem->d_long_seg, sizeof(seg_t) * num_long_seg,
                cudaMemcpyDeviceToHost);
    int32_t* long_range = (int32_t*)malloc(sizeof(int32_t) * *(long_mem->total_long_segs_n));
    cudaMemcpy(long_range, dev_mem->d_range_long, sizeof(int32_t) * *(long_mem->total_long_segs_n), cudaMemcpyDeviceToHost);

// Calculate total workload (sc pairs)
    size_t total_sc_pairs = 0;
    for (size_t i = 0; i < total_n; i++) {
        total_sc_pairs += range[i];
    }
    size_t long_seg_sc_pairs = 0;
    for (int long_seg_id = 0; long_seg_id < num_long_seg; long_seg_id++) {
        for(size_t i = long_segs[long_seg_id].start_idx; i < long_segs[long_seg_id].end_idx; i++){
            long_seg_sc_pairs += long_range[i];
        }
    }
    fprintf(stderr, "[DEBUG] Total workload (sc pairs) in batch: %lu, in long segs %lu\n", total_sc_pairs, long_seg_sc_pairs);


    // calculate kernel throughput
    float long_kernel_runtime_ms = 0;
    cudaEventElapsedTime(&long_kernel_runtime_ms, stream.long_kernel_event, stream.stopevent);
    float long_kernel_througput = long_seg_sc_pairs / long_kernel_runtime_ms / (float)1000;
    fprintf(stderr, "[DEBUG] Long Seg kernel throughput: %.2f Mpairs/s\n", long_kernel_througput);

// Check range w.r.t input (MAKE SURE INPUT RANGE EXISTS)
#if 0
    int64_t read_start = 0;
    for (int i = 0; i < dev_mem->size; i++) {
// DEBUG: print range
#if defined(DEBUG_VERBOSE) && 0
        debug_print_successor_range(range + read_start, reads[i].n);
#endif
        debug_check_range(range + read_start, input_arr[i].range, input_arr[i].n);
        read_start += reads[i].n;
    }
#endif

// DEBUG: Check voilation of cut
#if defined(DEBUG_CHECK) && 0
    for (int readid = 0, cid = 0, idx = 0; readid < dev_mem->size; readid++) {
// DEBUG: Print cuts
#if defined(DEBUG_VERBOSE) && 0
    debug_print_cut(cut + cid, cut_num - cid, reads[readid].n, idx, reads[readid].seq.name);
#endif
    cid += debug_check_cut(cut + cid, range, cut_num - cid, reads[readid].n, idx);
    idx += reads[readid].n;
    }
#endif

#if defined(DEBUG_VERBOSE) && 0
    int32_t* ax = (int32_t*) malloc(sizeof(int32_t) * dev_mem->buffer_size_long);
    cudaMemcpy(ax, dev_mem->d_ax_long, sizeof(int32_t) * dev_mem->buffer_size_long, cudaMemcpyDeviceToHost);
    debug_print_segs(host_mem->long_segs, reads, host_mem->long_segs_num[0], stream.host_mems[uid].size);
    debug_check_anchors(host_mem->long_segs, host_mem->long_segs_num[0], ax, host_mem->ax);
#endif

    free(cut);
    free(range);

// DEBUG: Calculate workload distribution
#if defined(DEBUG_VERBOSE) && 0
    plchain_cal_sc_pair_density(total_n, cut_num, dev_mem);
#endif // DEBUG_VERBOSE

//DEBUG: Calculate range distribution (MAKE SURE start_segid & end_segid EXISTS: search for LONG_SEG_RANGE_DIS)
#if defined(DEBUG_VERBOSE) && 0
        plchain_cal_long_seg_range_dis(total_n, dev_mem);
        plchain_cal_range_dis(total_n, cut_num, dev_mem);
#endif // DEBUG_VERBOSE

}

#endif // DEBUG_CHECK

/*
Accepts a stream that has already been synced. Finish and cleanup the stream, save to unpinned CPU memory.  
RETURN: number of reads in the batch 
*/
int plchain_finish_batch(streamSetup_t stream_setup, int stream_id, Misc misc, void* km){
    int n_reads = 0; // Number of reads in the batch
#ifdef DEBUG_PRINT
        float milliseconds = 0;
        cudaEventElapsedTime(&milliseconds, stream_setup.streams[stream_id].startevent, stream_setup.streams[stream_id].stopevent);
        // fprintf(stderr, "[Info] %s (%s:%d) last launch runtime: %f ms\n", __func__, __FILE__, __LINE__, milliseconds);
#endif

#ifdef DEBUG_PRINT
        float long_kernel_runtime_ms = 0;
        float short_kernel_runtime_ms[MICRO_BATCH];
        cudaEventElapsedTime(&long_kernel_runtime_ms, stream_setup.streams[stream_id].long_kernel_event, stream_setup.streams[stream_id].stopevent);
        float long_kernel_througput = *stream_setup.streams[stream_id].long_mem.total_long_segs_n / long_kernel_runtime_ms / (float)1000;
        for (int uid = 0; uid < MICRO_BATCH; uid++){
            cudaEventElapsedTime(&short_kernel_runtime_ms[uid], stream_setup.streams[stream_id].short_kernel_event[uid], 
            uid + 1 < MICRO_BATCH ? stream_setup.streams[stream_id].short_kernel_event[uid+1] : stream_setup.streams[stream_id].long_kernel_event);
        }
        fprintf(stderr, " %9.2f %%\n", (float)(*stream_setup.streams[stream_id].long_mem.total_long_segs_n)/stream_setup.long_seg_buffer_size_stream*100);
        fprintf(stderr, "Runtime(s) = ");
        for (int uid = 0; uid < MICRO_BATCH; uid++) fprintf(stderr, " %11.2f", short_kernel_runtime_ms[uid] / 1000);
        fprintf(stderr, " %9.2f %%\n", long_kernel_runtime_ms / 1000);
        fprintf(stderr, "Throutput (Ma/s) = ");
        for (int uid = 0; uid < MICRO_BATCH; uid++) fprintf(stderr, "           ");
        fprintf(stderr, " %11.2f\n", long_kernel_througput);
        fprintf(stderr, "[Info] %s finish (%s:%d) n_read %d long seg buffer usage %.2f%%. Long seg kernel throughput %.2f Manchors/s\n", __func__, __FILE__, __LINE__, n_reads, (float)(*stream_setup.streams[stream_id].long_mem.total_long_segs_n)/stream_setup.long_seg_buffer_size_stream*100, long_kernel_througput);
#endif // DEBUG_PRINT

// DEBUG: debug analysis that involves sychronizating the whole device 
#if defined(DEBUG_CHECK)
        plchain_debug_analysis(stream_setup.streams[stream_id]);
#endif  // DEBUG_VERBOSE
        seg_t* long_segs = stream_setup.streams[stream_id].long_mem.long_segs;
        size_t long_seg_idx = 0;
        for (int uid = 0; uid < MICRO_BATCH; uid++) {
            // regorg long to each host mem ptr
            // NOTE: this is the number of long segs till this microbatch
            size_t long_segs_num = stream_setup.streams[stream_id].host_mems[uid].long_segs_num[0];
#ifdef DEBUG_PRINT
            fprintf(stderr, "[Info] %s (%s:%d) MICROBATCH$%d FINISHED: long seg %lu - %lu\n", __func__, __FILE__, __LINE__, uid, long_seg_idx, long_segs_num);
#endif
            for (; long_seg_idx < long_segs_num; long_seg_idx++) {
                // TODO: write long_segs + long_seg_idx to f/p
            }

            // backtrack after p/f is copied
            plchain_backtracking(&stream_setup.streams[stream_id].host_mems[uid],
                                stream_setup.streams[stream_id].reads + n_reads, misc, km);
            // accumulate n_reads
            n_reads += stream_setup.streams[stream_id].host_mems[uid].size;
        }

    return n_reads;
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
        *n_read_ = plchain_finish_batch(stream_setup, stream_id, misc, km);
        *reads_ = stream_setup.streams[stream_id].reads;
        stream_setup.streams[stream_id].busy = false;
    }

    cudaEventRecord(stream_setup.streams[stream_id].startevent,
                    stream_setup.streams[stream_id].cudastream);
    size_t total_n = 0;
    for (int i = 0; i < n_read; i++) {
        total_n += reads[i].n;
    } // compute total_n first

    fprintf(stderr, "[Info] %s (%s:%d) Launching Batch: n_read %d, total anchors %lu (usage: %.2f%%)\n", __func__, __FILE__, __LINE__, n_read, total_n, (float)total_n/stream_setup.max_anchors_stream*100);


    // reset long seg counters
    cudaMemsetAsync(stream_setup.streams[stream_id].dev_mem.d_long_seg_count, 0, sizeof(unsigned int),
                    stream_setup.streams[stream_id].cudastream);
    cudaMemsetAsync(stream_setup.streams[stream_id].dev_mem.d_total_n_long, 0, sizeof(size_t),
                    stream_setup.streams[stream_id].cudastream);

    stream_setup.streams[stream_id].reads = reads;
    int read_start = 0;
#ifdef DEBUG_PRINT
    fprintf(stderr, "----------------------------------------------------------------\n           ");
    for (int uid = 0; uid < MICRO_BATCH; uid++){
        fprintf(stderr, "       Short%d", uid);
    }
    fprintf(stderr, "       Long\n");
	fprintf(stderr, "----------------------------------------------------------------\n");
    fprintf(stderr, "Mem Usage  = ");
#endif
    for (int uid = 0; uid < MICRO_BATCH; uid++) {
#ifdef USEHIP
        roctxRangePushA("microbatch");
#endif
        // decide the size of micro batch
        size_t batch_n = 0;
        int read_end = 0;
        size_t cut_num = 0;
        int griddim = 0;
        for (read_end = read_start; read_end < n_read; read_end++) {
            if (batch_n > (total_n - 1) / MICRO_BATCH + 1) {
                break;
            }
            batch_n += reads[read_end].n;
            int an_p_block = range_kernel_config.anchor_per_block;
            int an_p_cut = range_kernel_config.blockdim;
            int block_num = (reads[read_end].n - 1) / an_p_block + 1;
            griddim += block_num;
            cut_num += (reads[read_end].n - 1) / an_p_cut + 1;
        }
#ifdef DEBUG_PRINT
        // fprintf(stderr, "[Info] %s (%s:%d) MICROBATCH#%d batch_n %lu, read_start %d, read_end %d usage %.2f %%\n", __func__, __FILE__, __LINE__, uid, batch_n, read_start, read_end, (float)batch_n/stream_setup.max_anchors_stream*100);

        fprintf(stderr, " %9.2f %%", (float)batch_n/stream_setup.max_anchors_stream*100);
#endif // DEBUG_PRINT
        // // sanity check
        // if (stream_setup.max_anchors_stream < total_n){
        //     fprintf(stderr, "max_anchors_stream %lu total_n %lu n_read %d\n",
        //             stream_setup.max_anchors_stream, total_n, n_read);
        // }

        // if (stream_setup.max_range_grid < griddim) {
        //     fprintf(stderr, "max_range_grid %d griddim %d, total_n %lu n_read %d\n",
        //             stream_setup.max_range_grid, griddim, total_n, n_read);
        // }

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
        
        cudaEventRecord(stream_setup.streams[stream_id].short_kernel_event[uid],
                    stream_setup.streams[stream_id].cudastream);
        plmem_async_h2d_short_memcpy(&stream_setup.streams[stream_id], uid);
        // step3: range selection
        plrange_async_range_selection(&stream_setup.streams[stream_id].dev_mem,
                                    &stream_setup.streams[stream_id].cudastream);
        // step4: score generation for short and mid segs
        plscore_async_short_mid_forward_dp(&stream_setup.streams[stream_id].dev_mem,
                                    &stream_setup.streams[stream_id].cudastream);
        // step5: copy short and mid results back
        plmem_async_d2h_short_memcpy(&stream_setup.streams[stream_id], uid);
        // update index
        read_start = read_end;
#ifdef USEHIP
        roctxRangePop();
#endif
    }

    // step6: copy back long_segs_og
    cudaStreamSynchronize(stream_setup.streams[stream_id].cudastream);
    unsigned int num_long_seg;
    cudaMemcpy(&num_long_seg, stream_setup.streams[stream_id].dev_mem.d_long_seg_count, sizeof(unsigned int),
                cudaMemcpyDeviceToHost);
    seg_t* long_segs_og = (seg_t*)malloc(sizeof(seg_t) * num_long_seg);
    cudaMemcpy(long_segs_og, stream_setup.streams[stream_id].dev_mem.d_long_seg_og, sizeof(seg_t) * num_long_seg,
                cudaMemcpyDeviceToHost);
    #ifdef DEBUG_PRINT
    for (int i = 0; i < num_long_seg; i++){
        fprintf(stderr, "long seg %d: %lu - %lu\n", i, long_segs_og[i].start_idx, long_segs_og[i].end_idx);
    }
    #endif // DEBUG_PRINT

    // step7: sort long segs in descent order
    unsigned *map = new unsigned[num_long_seg];
    for (unsigned i = 0; i < num_long_seg; i++) {
        map[i] = i;
    }
    pairsort(long_segs_og, map, num_long_seg);
    #ifdef DEBUG_PRINT
    for (int i = 0; i < num_long_seg; i++){
        fprintf(stderr, "sorted index: %u, length: %zu\n", map[i], long_segs_og[map[i]].end_idx - long_segs_og[map[i]].start_idx);
    }
    #endif // DEBUG_PRINT
    free(long_segs_og);

    // step8: copy map to device
    cudaMalloc(&stream_setup.streams[stream_id].dev_mem.d_map, sizeof(unsigned) * num_long_seg);
    cudaMemcpy(stream_setup.streams[stream_id].dev_mem.d_map, map, sizeof(unsigned) * num_long_seg, cudaMemcpyHostToDevice);
    free(map);

    cudaEventRecord(stream_setup.streams[stream_id].long_kernel_event,
                    stream_setup.streams[stream_id].cudastream);
    plscore_async_long_forward_dp(&stream_setup.streams[stream_id].dev_mem,
                                   &stream_setup.streams[stream_id].cudastream);
    plmem_async_d2h_long_memcpy(&stream_setup.streams[stream_id]);
    cudaEventRecord(stream_setup.streams[stream_id].stopevent,
                    stream_setup.streams[stream_id].cudastream);
    stream_setup.streams[stream_id].busy = true;
    cudaCheck();
}

#ifdef __cplusplus
extern "C" {
#endif  // __cplusplus

void init_blocking_gpu(size_t* total_n, int* max_reads, int *min_n, Misc misc) {
    plmem_initialize(total_n, max_reads, min_n);
}

void init_stream_gpu(size_t* total_n, int* max_reads, int *min_n, Misc misc) {
#ifdef DEBUG_PRINT
    fprintf(stderr, "[INFO::%s] gpu initialized for chaining\n", __func__);
#endif // DEBUG_PRINT
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

    n_read = plchain_finish_batch(stream_setup, t, misc, km);
    reads = stream_setup.streams[t].reads;
#ifdef DEBUG_PRINT
    fprintf(stderr, "[Debug] %s finish (%s:%d) n_read %d\n", __func__, __FILE__, __LINE__, n_read);
#endif // DEBUG_PRINT
    stream_setup.streams[t].busy = false;

    for (int i = 0; i < n_read; i++) {
        post_chaining_helper(mi, opt, &reads[i], misc, km);
    }
    // FIXME: return an array of reads
    *reads_ = reads;
    *n_read_ = n_read;

}


void free_stream_gpu(int n_threads){
    // for (int t = 0; t < n_threads; t++){
    //     for (int i = 0; i < MICRO_BATCH; i++){
    //         plmem_free_host_mem(&stream_setup.streams[t].host_mems[i]);
    //     }
    //     plmem_free_long_mem(&stream_setup.streams[t].long_mem);
    //     plmem_free_device_mem(&stream_setup.streams[t].dev_mem);
    // }
    plmem_stream_cleanup();
#ifdef DEBUG_PRINT
    fprintf(stderr, "[Info::%s] gpu free memory\n", __func__);
#endif // DEBUG_PRINT
}

#ifdef __cplusplus
} // extern "C"
#endif  // __cplusplus