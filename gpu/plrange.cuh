#ifndef _PLRANGE_CUH_
#define _PLRANGE_CUH_

#include "plmem.cuh"
#include <assert.h>

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


#ifdef __cplusplus
}
#endif

#endif  // _PLRANGE_CUH_
