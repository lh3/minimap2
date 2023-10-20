#ifndef _PLSCORE_CUH_
#define _PLSCORE_CUH_

#include "plmem.cuh"
#include "mmpriv.h"

#ifdef __cplusplus
extern "C"{
#endif

#ifndef MM_MID_SEG_CUTOFF
#define MM_MID_SEG_CUTOFF 1
#endif

#ifndef MM_LONG_SEG_CUTOFF
#define MM_LONG_SEG_CUTOFF 10
#endif

#ifndef MM_CUT_SIZE
#define MM_CUT_SIZE 512
#endif

#define MM_QSPAN 15 

void plscore_upload_misc(Misc misc);
void plscore_async_naive_forward_dp(deviceMemPtr* dev_mem, cudaStream_t* stream);
void plscore_async_long_short_forward_dp(deviceMemPtr* dev_mem,cudaStream_t* stream);
void plscore_sync_long_short_forward_dp(deviceMemPtr* dev_mem, Misc misc_);
void plscore_sync_naive_forward_dp(deviceMemPtr* dev_mem, Misc misc_);

extern score_kernel_config_t score_kernel_config;

#ifdef __cplusplus
}
#endif

#endif // _PLSCORE_CUH_