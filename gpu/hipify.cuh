#ifndef __HIPIFY_CUH__
#define __HIPIFY_CUH__

#ifdef USEHIP
#include "hip/hip_runtime.h"
#define     cudaDeviceProp                  hipDeviceProp_t
#define     cudaGetDeviceProperties         hipGetDeviceProperties
#define     cudaMalloc                      hipMalloc
#define     cudaMallocAsync                 hipMallocAsync
#define     cudaMemcpy                      hipMemcpy
#define     cudaMemcpyAsync                 hipMemcpyAsync
#define     cudaMemcpyToSymbolAsync         hipMemcpyToSymbolAsync
#define     cudaMemcpyHostToDevice          hipMemcpyHostToDevice
#define     cudaMemcpyDeviceToHost          hipMemcpyDeviceToHost
#define     cudaDeviceSynchronize           hipDeviceSynchronize
#define     cudaFree                        hipFree
#define     cudaFreeAsync                   hipFreeAsync
#define     cudaMemcpyToSymbol              hipMemcpyToSymbol
#define     cudaMemset                      hipMemset
#define     cudaMemsetAsync                 hipMemsetAsync
#define     cudaStream_t                    hipStream_t
#define     cudaStreamCreate                hipStreamCreate
#define     cudaStreamSynchronize           hipStreamSynchronize
#define     cudaStreamDestroy               hipStreamDestroy
#define     cudaMallocHost                  hipHostMalloc
#define     cudaFreeHost                    hipHostFree
#define     cudaEvent_t                     hipEvent_t
#define     cudaEventCreate                 hipEventCreate
#define     cudaEventRecord                 hipEventRecord
#define     cudaEventQuery                  hipEventQuery
#define     cudaEventDestroy                hipEventDestroy
#define cudaCheck() {                                                       \
    hipError_t err = hipGetLastError();                                     \
    if (hipSuccess != err) {                                                \
        fprintf(stderr, "Error in %s:%i %s(): %s.\n", __FILE__, __LINE__,   \
                __func__, hipGetErrorString(err));                          \
        fflush(stderr);                                                     \
        exit(EXIT_FAILURE);                                                 \
    }                                                                       \
}
#else
#define cudaCheck() {                                                   \
    cudaError_t err = cudaGetLastError();                               \
    if (cudaSuccess != err) {                                           \
        fprintf(stderr, "Error in %s:%i %s(): %s.\n", __FILE__, __LINE__,\
                __func__, cudaGetErrorString(err));                     \
        fflush(stderr);                                                 \
        exit(EXIT_FAILURE);                                             \
    }                                                                   \
}
#include <cuda.h>

#endif


#endif // __HIPIFY_CUH__
