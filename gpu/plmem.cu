#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include "plmem.cuh"
#include "plrange.cuh"
#include "plscore.cuh"
#include <time.h>

void plmem_malloc_host_mem(hostMemPtr *host_mem, size_t anchor_per_batch,
                           int range_grid_size) {
    // data array
    cudaMallocHost((void**)&host_mem->ax, anchor_per_batch * sizeof(int32_t));
    cudaMallocHost((void**)&host_mem->ay, anchor_per_batch * sizeof(int32_t));
    cudaMallocHost((void**)&host_mem->sid, anchor_per_batch * sizeof(int8_t));
    cudaMallocHost((void**)&host_mem->xrev, anchor_per_batch * sizeof(int32_t));
    cudaMallocHost((void**)&host_mem->f, anchor_per_batch * sizeof(int32_t));
    cudaMallocHost((void**)&host_mem->p, anchor_per_batch * sizeof(uint16_t));

    //index
    cudaMallocHost((void**)&host_mem->start_idx, range_grid_size * sizeof(size_t));
    cudaMallocHost((void**)&host_mem->read_end_idx, range_grid_size * sizeof(size_t));
    cudaMallocHost((void**)&host_mem->cut_start_idx, range_grid_size * sizeof(size_t));
    cudaCheck();
}

void plmem_free_host_mem(hostMemPtr *host_mem) { 
    cudaFreeHost(host_mem->ax);
    cudaFreeHost(host_mem->ay);
    cudaFreeHost(host_mem->f);
    cudaFreeHost(host_mem->p);
    cudaFreeHost(host_mem->start_idx);
    cudaFreeHost(host_mem->read_end_idx);
    cudaFreeHost(host_mem->cut_start_idx);
    cudaCheck();
}

void plmem_malloc_device_mem(deviceMemPtr *dev_mem, size_t anchor_per_batch, int range_grid_size, int num_cut){
    // data array
    cudaMalloc(&dev_mem->d_ax, anchor_per_batch * sizeof(int32_t));
    cudaMalloc(&dev_mem->d_ay, anchor_per_batch * sizeof(int32_t));
    cudaMalloc(&dev_mem->d_sid, anchor_per_batch * sizeof(int8_t));
    cudaMalloc(&dev_mem->d_xrev, anchor_per_batch * sizeof(int32_t));
    cudaMalloc(&dev_mem->d_range, anchor_per_batch * sizeof(int32_t));
    cudaMalloc(&dev_mem->d_f, anchor_per_batch * sizeof(int32_t));
    cudaMalloc(&dev_mem->d_p, anchor_per_batch * sizeof(uint16_t));

    //index
    cudaMalloc(&dev_mem->d_start_idx, range_grid_size * sizeof(size_t));
    cudaMalloc(&dev_mem->d_read_end_idx, range_grid_size * sizeof(size_t));
    cudaMalloc(&dev_mem->d_cut_start_idx, range_grid_size * sizeof(size_t));

    // cut
    cudaMalloc(&dev_mem->d_cut, num_cut * sizeof(size_t));
    cudaMalloc(&dev_mem->d_long_seg_count, sizeof(unsigned int));
    cudaMalloc(&dev_mem->d_long_seg, num_cut/(MM_LONG_SEG_CUTOFF + 1) * sizeof(seg_t));
    cudaMalloc(&dev_mem->d_mid_seg_count, sizeof(unsigned int));
    cudaMalloc(&dev_mem->d_mid_seg, num_cut/(MM_MID_SEG_CUTOFF + 1) * sizeof(seg_t));
    cudaCheck();
}

void plmem_free_device_mem(deviceMemPtr *dev_mem) { 
    cudaFree(dev_mem->d_ax);
    cudaFree(dev_mem->d_ay);
    cudaFree(dev_mem->d_sid);
    cudaFree(dev_mem->d_xrev);
    cudaFree(dev_mem->d_range);
    cudaFree(dev_mem->d_f);
    cudaFree(dev_mem->d_p);

    cudaFree(dev_mem->d_start_idx);
    cudaFree(dev_mem->d_read_end_idx);
    cudaFree(dev_mem->d_cut_start_idx);

    cudaFree(dev_mem->d_cut);
    cudaFree(dev_mem->d_long_seg);
    cudaFree(dev_mem->d_long_seg_count);
    cudaFree(dev_mem->d_mid_seg);
    cudaFree(dev_mem->d_mid_seg_count);
    cudaCheck();
}


/**
 * Input
 *  reads[]:    array
 *  n_reads
 *  config:     range kernel configuartions
 * Output
 *  *host_mem   populate host_mem
*/
void plmem_reorg_input_arr(chain_read_t *reads, int n_read,
                           hostMemPtr *host_mem, range_kernel_config_t config) {
    size_t total_n = 0, cut_num = 0;
    size_t griddim = 0;

    host_mem->size = n_read;
    for (int i = 0; i < n_read; i++) {
        total_n += reads[i].n;
    }
    host_mem->total_n = total_n;

    size_t idx = 0;
    for (int i = 0; i < n_read; i++) {
        int n = reads[i].n;
        int block_num = (n - 1) / config.anchor_per_block + 1;

        host_mem->start_idx[griddim] = idx;
        size_t end_idx = idx + config.anchor_per_block;
        host_mem->read_end_idx[griddim] = idx + n;
        host_mem->cut_start_idx[griddim] = cut_num;
        for (int j = 1; j < block_num; j++) {
            cut_num += (config.anchor_per_block / config.blockdim);
            host_mem->start_idx[griddim + j] = end_idx;
            end_idx =
                host_mem->start_idx[griddim + j] + config.anchor_per_block;
            host_mem->read_end_idx[griddim + j] = idx + n;
            host_mem->cut_start_idx[griddim + j] = cut_num;
        }
        cut_num += (n - (block_num - 1) * config.anchor_per_block - 1) /
                       config.blockdim;
        end_idx = idx + n;

        griddim += block_num;

        for (int j = 0; j < n; j++) {
            host_mem->ax[idx] = (int32_t)reads[i].a[j].x;
            host_mem->ay[idx] = (int32_t)reads[i].a[j].y;
            host_mem->sid[idx] = (reads[i].a[j].y & MM_SEED_SEG_MASK) >> MM_SEED_SEG_SHIFT;
            host_mem->xrev[idx] = reads[i].a[j].x >> 32;
            ++idx;
        }
    }
    host_mem->cut_num = cut_num;
    host_mem->griddim = griddim;
}

void plmem_async_h2d_memcpy(stream_ptr_t* stream_ptrs) {
    hostMemPtr *host_mem = &stream_ptrs->host_mem;
    deviceMemPtr *dev_mem = &stream_ptrs->dev_mem;
    cudaStream_t *stream = &stream_ptrs->cudastream;
    cudaMemcpyAsync(dev_mem->d_ax, host_mem->ax,
                    sizeof(int32_t) * host_mem->total_n, cudaMemcpyHostToDevice,
                    *stream);
    cudaMemcpyAsync(dev_mem->d_ay, host_mem->ay,
                    sizeof(int32_t) * host_mem->total_n, cudaMemcpyHostToDevice,
                    *stream);
    cudaMemcpyAsync(dev_mem->d_sid, host_mem->sid,
                    sizeof(int8_t) * host_mem->total_n, cudaMemcpyHostToDevice,
                    *stream);
    cudaMemcpyAsync(dev_mem->d_xrev, host_mem->xrev,
                    sizeof(int32_t) * host_mem->total_n, cudaMemcpyHostToDevice,
                    *stream);
    cudaMemcpyAsync(dev_mem->d_start_idx, host_mem->start_idx,
                    sizeof(size_t) * host_mem->griddim, cudaMemcpyHostToDevice,
                    *stream);
    cudaMemcpyAsync(dev_mem->d_read_end_idx, host_mem->read_end_idx,
                    sizeof(size_t) * host_mem->griddim, cudaMemcpyHostToDevice,
                    *stream);
    cudaMemcpyAsync(dev_mem->d_cut_start_idx, host_mem->cut_start_idx,
                    sizeof(size_t) * host_mem->griddim, cudaMemcpyHostToDevice,
                    *stream);
    cudaMemsetAsync(dev_mem->d_cut, 0xff,
                    sizeof(size_t) * host_mem->cut_num, *stream);
    cudaMemsetAsync(dev_mem->d_f, 0, sizeof(int32_t) * host_mem->total_n,
                    *stream);
    cudaMemsetAsync(dev_mem->d_p, 0, sizeof(uint16_t) * host_mem->total_n,
                    *stream);
    cudaCheck();
    dev_mem->total_n = host_mem->total_n;
    dev_mem->num_cut = host_mem->cut_num;
    dev_mem->size = host_mem->size;
    dev_mem->griddim = host_mem->griddim;
}

void plmem_sync_h2d_memcpy(hostMemPtr *host_mem, deviceMemPtr *dev_mem) {
    cudaMemcpy(dev_mem->d_ax, host_mem->ax, sizeof(int32_t) * host_mem->total_n,
               cudaMemcpyHostToDevice);
    cudaMemcpy(dev_mem->d_ay, host_mem->ay, sizeof(int32_t) * host_mem->total_n,
               cudaMemcpyHostToDevice);
    cudaMemcpy(dev_mem->d_sid, host_mem->sid, sizeof(int8_t) * host_mem->total_n,
               cudaMemcpyHostToDevice);
    cudaMemcpy(dev_mem->d_xrev, host_mem->xrev,
               sizeof(int32_t) * host_mem->total_n, cudaMemcpyHostToDevice);
    cudaMemcpy(dev_mem->d_start_idx, host_mem->start_idx,
               sizeof(size_t) * host_mem->griddim, cudaMemcpyHostToDevice);
    cudaMemcpy(dev_mem->d_read_end_idx, host_mem->read_end_idx,
               sizeof(size_t) * host_mem->griddim, cudaMemcpyHostToDevice);
    cudaMemcpy(dev_mem->d_cut_start_idx, host_mem->cut_start_idx,
               sizeof(size_t) * host_mem->griddim, cudaMemcpyHostToDevice);
    cudaMemset(dev_mem->d_cut, 0xff, sizeof(size_t) * host_mem->cut_num);
    dev_mem->total_n = host_mem->total_n;
    dev_mem->num_cut = host_mem->cut_num;
    dev_mem->size = host_mem->size;
    dev_mem->griddim = host_mem->griddim;
    cudaCheck();
}

void plmem_async_d2h_memcpy(stream_ptr_t *stream_ptrs) {
    hostMemPtr *host_mem = &stream_ptrs->host_mem;
    deviceMemPtr *dev_mem = &stream_ptrs->dev_mem;
    cudaStream_t *stream = &stream_ptrs->cudastream;
    cudaMemcpyAsync(host_mem->f, dev_mem->d_f,
                    sizeof(int32_t) * host_mem->total_n, cudaMemcpyDeviceToHost,
                    *stream);
    cudaMemcpyAsync(host_mem->p, dev_mem->d_p,
                    sizeof(uint16_t) * host_mem->total_n,
                    cudaMemcpyDeviceToHost, *stream);
    cudaCheck();
}

void plmem_sync_d2h_memcpy(hostMemPtr *host_mem, deviceMemPtr *dev_mem){
    cudaMemcpy(host_mem->f, dev_mem->d_f, sizeof(int32_t) * host_mem->total_n,
               cudaMemcpyDeviceToHost);
    cudaMemcpy(host_mem->p, dev_mem->d_p, sizeof(uint16_t) * host_mem->total_n,
               cudaMemcpyDeviceToHost);
    cudaCheck();
}

//////////////////// Initialization and Cleanup <mmpriv.h>////////////////////////
streamSetup_t stream_setup;

#include "cJSON.h"
cJSON *plmem_parse_gpu_config(const char filename[]){
    // read json file to cstring
    char *buffer = 0;
    long length;
    FILE *f = fopen(filename, "rb");

    if (f) {
        fseek(f, 0, SEEK_END);
        length = ftell(f);
        fseek(f, 0, SEEK_SET);
        buffer = (char*)malloc(length);
        if (buffer) {
            fread(buffer, 1, length, f);
        }
        fclose(f);
    }

    if (!buffer) {
        fprintf(stderr, "[Error] fail to open gpu config file %s\n", filename);
        exit(1);
    }

    cJSON *json = cJSON_Parse(buffer);
    if (!json) {
        const char *error_ptr = cJSON_GetErrorPtr();
        if (error_ptr != NULL) {
            fprintf(stderr, "[Error] cJSON error before %s\n", error_ptr);
        }
        exit(1);
    }

    return json;
}

int get_json_int(cJSON *json, const char name[]) {
    cJSON *elt = cJSON_GetObjectItem(json, name);
    if (!cJSON_IsNumber(elt)) {
        fprintf(stderr, "[Error] cJSON error failed to get field %s\n", name);
        exit(1);
    }
    return elt->valueint;
}

void plmem_config_kernels(cJSON *json) {
    cJSON *range_config_json = cJSON_GetObjectItem(json, "range_kernel");
    range_kernel_config.blockdim = get_json_int(range_config_json, "blockdim");
    range_kernel_config.cut_check_anchors =
        get_json_int(range_config_json, "cut_check_anchors");
    range_kernel_config.anchor_per_block =
        get_json_int(range_config_json, "anchor_per_block");

    cJSON *score_config_json = cJSON_GetObjectItem(json, "score_kernel");
    score_kernel_config.short_blockdim =
        get_json_int(score_config_json, "short_blockdim");
    score_kernel_config.long_blockdim =
        get_json_int(score_config_json, "long_blockdim");
    score_kernel_config.mid_blockdim =
        get_json_int(score_config_json, "mid_blockdim");
    score_kernel_config.short_griddim =
        get_json_int(score_config_json, "short_griddim");
    score_kernel_config.long_griddim =
        get_json_int(score_config_json, "long_griddim");
    score_kernel_config.mid_griddim =
        get_json_int(score_config_json, "mid_griddim");
}

void plmem_config_stream(size_t *max_range_grid_, size_t *max_num_cut_, size_t max_total_n, size_t max_read, size_t min_n){
    size_t max_range_grid, max_num_cut;
    max_range_grid =
        (max_total_n - 1) / range_kernel_config.anchor_per_block + 1 + max_read;
    max_num_cut = (max_total_n - 1) / range_kernel_config.blockdim + 1 + max_read;
    *max_range_grid_ = max_range_grid;
    *max_num_cut_ = max_num_cut;

    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, 0);
    cudaCheck();

    if (*max_range_grid_ > prop.maxGridSize[0]) {
        fprintf(stderr, "Invalid memory config!\n");
        exit(1);
    }

    fprintf(stderr,
            "[M: %s] max_grid: %zu, max_anchors_per_stream: %zu, "
            "max_num_cut_per_stream: %zu \n",
            __func__, *max_range_grid_, max_total_n, *max_num_cut_);
}


template <bool is_blooking = false>
void plmem_config_batch(cJSON *json, int *num_stream_,
                         int *min_n_, size_t *max_total_n_,
                         int *max_read_) {
    if (is_blooking) 
        *num_stream_ = 16;
    else
        *num_stream_ = get_json_int(json, "num_streams");

    size_t min_anchors = get_json_int(json, "min_n");
    *min_n_ = min_anchors;

    /* If Use define max_total_n & max_read */
    cJSON *max_total_n_json = cJSON_GetObjectItem(json, "max_total_n");
    cJSON *max_read_json = cJSON_GetObjectItem(json, "max_read");
    if (max_total_n_json && max_read_json){
        *max_total_n_ = max_total_n_json->valueint;
        *max_read_ = max_read_json->valueint;
        return;
    }

    /* Determine configuration smartly */
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, 0);

    size_t avail_mem_per_stream = (prop.totalGlobalMem / *num_stream_ ) * 0.8;

    fprintf(stderr, "[M: %s] mem_per_stream: %zu x %d streams \n", __func__,
            avail_mem_per_stream, *num_stream_);

    // memory per anchor = (ax + ay + range + f + p) + (start_idx + read_end_idx
    // + cut_start_idx) + cut + long_seg size: F1 = ax + ay + range + f + p; F2
    // = start_idx + read_end_idx + cut_start_idx; F3 = cut; F4 = long_seg
    int F1 = 8 + 8 + 4 + 4 + 2, F2 = 8 + 8 + 8, F3 = 8, F4 = 16;
    // TODO: define these data types

    // max iteration of each block, must be an integer
    // int max_it =
    //     range_kernel_config.anchor_per_block / range_kernel_config.blockdim;
    // int blockdim = range_kernel_config.blockdim;

    // g = max_grid_size
    // g * F2 + g * blockdim * max_it * F1 + max_cut * F3 + max_cut/2 * F4 <
    // mem_per_stream max_cut = g * max_it
    /*
    size_t cost_per_anchor = F1;
    size_t cost_per_grid = F2;
    size_t cost_per_cut = F3 + F4 / 2;
    */
    size_t avg_read_n = get_json_int(json, "avg_read_n");

    /**
     * Assume max_total_n = max_read * avg_read_n
     * max_grid = max_total_n / range.anchor_per_block + max_read
     *          = max_read *( avg_read_n / anchor_per_block + 1)
     * max_cut  = max_total_n / range.blockdim + max_read
     *          = max_read * ( avg_read_n / blockdim + 1)
     * total_mem = max_grid * cost_per_grid + max_cut * cost_per_cut + max_total_n * cost_per_anchor
     */

    float grid_cost_per_read =
        (avg_read_n / (float)range_kernel_config.anchor_per_block + 1) * F2;
    float cut_cost_per_read =
        (avg_read_n / (float)range_kernel_config.blockdim + 1) * (F3 + F4 / 2);
    *max_read_ = floor(avail_mem_per_stream /
                 (grid_cost_per_read + cut_cost_per_read + F1 * avg_read_n));
    *max_total_n_ = *max_read_ * avg_read_n;
}

// intialize and config kernels for gpu blocking setup
void plmem_initialize(size_t *max_total_n_, int *max_read_,
                      int *min_anchors_) {
    cJSON *json = plmem_parse_gpu_config("gpu_config.json");
    plmem_config_kernels(json);
    int num_streams;
    plmem_config_batch<true>(json, &num_streams, min_anchors_,
                              max_total_n_, max_read_);
}

// initialize global variable stream_setup
void plmem_stream_initialize(size_t *max_total_n_,
                             int *max_read_, int *min_anchors_) {

    int num_stream;
    size_t max_anchors_stream, max_range_grid, max_num_cut;
#ifndef GPU_CONFIG
    cJSON *json = plmem_parse_gpu_config("gpu_config.json");
#else 
    cJSON *json = plmem_parse_gpu_config(GPU_CONFIG);
#endif
    plmem_config_kernels(json);
    plmem_config_batch<false>(json, &num_stream, min_anchors_, &max_anchors_stream,
                              max_read_);
    plmem_config_stream(&max_range_grid, &max_num_cut, max_anchors_stream,
                        *max_read_, *min_anchors_);

    stream_setup.num_stream = num_stream;
    // assert(num_stream > 1);

    stream_setup.streams = new stream_ptr_t[num_stream];
#ifdef DEBUG_CHECK
    fprintf(
        stderr,
        "[Info] max anchors per stream: %zu, max range grid %zu max_num_cut %zu\n", max_anchors_stream, max_range_grid, max_num_cut);
#endif  // DEBUG_CHECK

    for (int i = 0; i < num_stream; i++) {
        stream_setup.streams[i].busy = false;
        cudaStreamCreate(&stream_setup.streams[i].cudastream);
        cudaEventCreate(&stream_setup.streams[i].cudaevent);
        cudaEventCreate(&stream_setup.streams[i].startevent);
        cudaCheck();
        plmem_malloc_host_mem(&stream_setup.streams[i].host_mem, max_anchors_stream,
                              max_range_grid);
        plmem_malloc_device_mem(&stream_setup.streams[i].dev_mem, max_anchors_stream,
                                max_range_grid, max_num_cut);
    }

    *max_total_n_ = max_anchors_stream;

    stream_setup.max_anchors_stream = max_anchors_stream;
    stream_setup.max_range_grid = max_range_grid;
    stream_setup.max_num_cut = max_num_cut;
}

void plmem_stream_cleanup() {
    for (int i = 0; i < stream_setup.num_stream; i++) {
        cudaStreamDestroy(stream_setup.streams[i].cudastream);
        cudaEventDestroy(stream_setup.streams[i].cudaevent);
        cudaCheck();
        plmem_free_host_mem(&stream_setup.streams[i].host_mem);
        plmem_free_device_mem(&stream_setup.streams[i].dev_mem);
    }
    delete[] stream_setup.streams;
}
