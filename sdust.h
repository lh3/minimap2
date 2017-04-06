#ifndef SDUST_H
#define SDUST_H

struct sdust_buf_s;
typedef struct sdust_buf_s sdust_buf_t;

#ifdef __cplusplus
extern "C" {
#endif

// the simple interface
uint64_t *sdust(const uint8_t *seq, int l_seq, int T, int W, int *n);

// the following interface dramatically reduce heap allocations when sdust is frequently called.
sdust_buf_t *sdust_buf_init(void);
void sdust_buf_destroy(sdust_buf_t *buf);
const uint64_t *sdust_core(const uint8_t *seq, int l_seq, int T, int W, int *n, sdust_buf_t *buf);

#ifdef __cplusplus
}
#endif

#endif
