/* Implement kernel performance analysis that requires extra device
 * synchornization. disabled unless DEBUG_LEVEL is set to analyze.
 * Enable individual verbose prints in planalyze.cu 
 */
#include "planalyze.cuh"

#ifdef DEBUG_CHECK
void planalyze_short_kernel(stream_ptr_t stream, int uid, float throughput[]){
    cudaStreamSynchronize(stream.cudastream);
    size_t total_n = stream.host_mems[uid].total_n;
    chain_read_t* reads = stream.reads;
    deviceMemPtr* dev_mem = &stream.dev_mem;
    hostMemPtr* host_mem = &stream.host_mems[uid];
    size_t cut_num = stream.host_mems[uid].cut_num;
    unsigned int num_mid_seg, num_long_seg;
    cudaMemcpy(&num_mid_seg, dev_mem->d_mid_seg_count, sizeof(unsigned int),
        cudaMemcpyDeviceToHost);
    num_long_seg = host_mem->long_segs_num[0] - (uid>0 ? stream.host_mems[uid-1].long_segs_num[0] : 0);
    cudaMemcpy(&num_long_seg, dev_mem->d_long_seg_count, sizeof(unsigned int),
        cudaMemcpyDeviceToHost);
#ifdef DEBUG_VERBOSE
    fprintf(stderr, "[DEBUG](MICROBATCH# %d) total segs: %lu, short:%lu mid: %u \n", uid, cut_num, cut_num - num_mid_seg - num_long_seg, num_mid_seg);
#endif // DEBUG_VERBOSE

    int32_t* range = (int32_t*)malloc(sizeof(int32_t) * total_n);
    cudaMemcpy(range, dev_mem->d_range, sizeof(int32_t) * total_n,
                cudaMemcpyDeviceToHost);
    size_t* cut = (size_t*)malloc(sizeof(size_t) * cut_num);
    cudaMemcpy(cut, dev_mem->d_cut, sizeof(size_t) * cut_num,
                cudaMemcpyDeviceToHost);

    seg_t* mid_segs = (seg_t*)malloc(sizeof(seg_t) * num_mid_seg);
    cudaMemcpy(mid_segs, dev_mem->d_mid_seg, sizeof(seg_t) * num_mid_seg,
                cudaMemcpyDeviceToHost);


    longMemPtr* long_mem = &stream.long_mem;
    unsigned int num_aggregated_long_segs;
    cudaMemcpy(&num_aggregated_long_segs, dev_mem->d_long_seg_count, sizeof(unsigned int),
                cudaMemcpyDeviceToHost);
#ifdef DEBUG_VERBOSE
    fprintf(stderr,
            "[DEBUG] aggreagated num of long segs %u, %u-%u belongs to this "
            "minibatch\n",
            num_aggregated_long_segs,
            uid > 0 ? stream.host_mems[uid-1].long_segs_num[0] : 0,
            num_aggregated_long_segs);
#endif // DEBUG_VERBOSE

    seg_t* long_segs = (seg_t*)malloc(sizeof(seg_t) * num_aggregated_long_segs);
    cudaMemcpy(long_segs, dev_mem->d_long_seg, sizeof(seg_t) * num_aggregated_long_segs,
                cudaMemcpyDeviceToHost);

    size_t long_segs_total_n;
    cudaMemcpy(&long_segs_total_n, dev_mem->d_total_n_long, sizeof(size_t), cudaMemcpyDeviceToHost);
    int32_t* long_range = (int32_t*)malloc(sizeof(int32_t) * long_segs_total_n);
    cudaMemcpy(long_range, dev_mem->d_range_long, sizeof(int32_t) * long_segs_total_n, cudaMemcpyDeviceToHost);

// Calculate long segs total workload (sc pairs)
    size_t long_seg_sc_pairs = 0;
    for(unsigned int segid = (uid>0 ? stream.host_mems[uid-1].long_segs_num[0] : 0); segid < num_aggregated_long_segs; segid++){
        for (size_t i = long_segs[segid].start_idx; i < long_segs[segid].end_idx; i++)
            long_seg_sc_pairs += long_range[i];
    }
#ifdef DEBUG_VERBOSE
    fprintf(stderr, "[DEBUG] workload (sc pairs) in long segs: %lu\n", long_seg_sc_pairs);
#endif // DEBUG_VERBOSE

// Calculate total workload (sc pairs)
    size_t total_sc_pairs = 0;
    for (size_t i = 0; i < total_n; i++) {
        total_sc_pairs += range[i];
    }

#ifdef DEBUG_VERBOSE
    fprintf(stderr, "[DEBUG] Total workload (sc pairs) in batch: %lu. %.2f%% work are in long segs.\n", total_sc_pairs, (float)long_seg_sc_pairs/total_sc_pairs*100);
#endif // DEBUG_VERBOSE
    assert(long_seg_sc_pairs <= total_sc_pairs); 

    // calculate short kernel throughput
    float short_kernel_runtime_ms = 0;
    cudaEventElapsedTime(&short_kernel_runtime_ms, stream.short_kernel_start_event[uid], stream.short_kernel_stop_event[uid]);
    throughput[uid] = (total_sc_pairs - long_seg_sc_pairs) / short_kernel_runtime_ms / (float)1000;
#ifdef DEBUG_VERBOSE
    fprintf(stderr, "[DEBUG] Short Seg kernel #%d throughput: %.2f Mpairs/s\n", uid, throughput[uid]);
#endif // DEBUG_VERBOSE

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

//DEBUG: Calculate range distribution
#if defined(DEBUG_VERBOSE) && 0
        debug_cal_range_dis(total_n, cut_num, range);
#endif // DEBUG_VERBOSE

// Calculate range distribution for mid segs
#if defined(DEBUG_VERBOSE) && 0
    for (int seg_id = 0; seg_id < num_mid_seg; seg_id++){
        debug_cal_mid_range_dis(mid_segs[seg_id].end_idx - mid_segs[seg_id].start_idx, 1, range + mid_segs[seg_id].start_idx);
    }
#endif // DEBUG_VERBOSE

// DEBUG: Calculate workload distribution
#if defined(DEBUG_VERBOSE) && 0
    debug_cal_sc_pair_density(total_n, cut_num, cut, range);
#endif // DEBUG_VERBOSE

    free(cut);
    free(range);

}
#endif 




#ifdef DEBUG_CHECK

void planalyze_long_kernel(stream_ptr_t stream, float* throughput){
    deviceMemPtr* dev_mem = &stream.dev_mem;
    longMemPtr* long_mem = &stream.long_mem;

    unsigned int num_long_seg = long_mem->total_long_segs_num[0];
#ifdef DEBUG_VERBOSE
    fprintf(stderr, "[DEBUG]LONG Kernel: num of long segs %u\n", num_long_seg);
#endif // DEBUG_VERBOSE


    seg_t* long_segs = (seg_t*)malloc(sizeof(seg_t) * num_long_seg);
    cudaMemcpy(long_segs, dev_mem->d_long_seg, sizeof(seg_t) * num_long_seg,
                cudaMemcpyDeviceToHost);
    int32_t* long_range = (int32_t*)malloc(sizeof(int32_t) * *(long_mem->total_long_segs_n));
    cudaMemcpy(long_range, dev_mem->d_range_long, sizeof(int32_t) * *(long_mem->total_long_segs_n), cudaMemcpyDeviceToHost);
#ifdef DEBUG_VERBOSE
    fprintf(stderr, "[DEBUG] Total n of anchors in long segs %lu\n", *long_mem->total_long_segs_n);
#endif // DEBUG_VERBOSE

// Calculate total workload (sc pairs)
    size_t long_seg_sc_pairs = 0;
    for(size_t i = 0; i < *long_mem->total_long_segs_n; i++){
        long_seg_sc_pairs += long_range[i];
    }
#ifdef DEBUG_VERBOSE
    fprintf(stderr, "[DEBUG] workload (sc pairs) in long kernel: %lu\n", long_seg_sc_pairs);
#endif // DEBUG_VERBOSE


    // calculate long kernel throughput
    float long_kernel_runtime_ms = 0;
    cudaEventElapsedTime(&long_kernel_runtime_ms, stream.long_kernel_event, stream.stopevent);
    float long_kernel_througput = long_seg_sc_pairs / long_kernel_runtime_ms / (float)1000;
#ifdef DEBUG_VERBOSE
    fprintf(stderr, "[DEBUG] Long Seg kernel throughput: %.2f Mpairs/s\n", long_kernel_througput);
#endif // DEBUG_VERBOSE
    throughput[score_kernel_config.micro_batch] = long_kernel_througput;

// Check range w.r.t input (MAKE SURE INPUT RANGE EXISTS)
#if defined(DEBUG_VERBOSE) && 0
    debug_print_successor_range(long_range, *long_mem->total_long_segs_n);
#endif 

//DEBUG: Calculate range distribution
#if defined(DEBUG_VERBOSE) && 0
    debug_cal_long_seg_range_dis(*long_mem->total_long_segs_n, num_long_seg, long_range);
#endif // DEBUG_VERBOSE

    free(long_segs);
    free(long_range);
    
}

#endif // DEBUG_CHECK