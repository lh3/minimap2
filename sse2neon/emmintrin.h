#ifndef SSE2NEON_EMMINTRIN_SHIM_H
#define SSE2NEON_EMMINTRIN_SHIM_H
/* minimap2 ARM/NEON build: route <emmintrin.h> to the modern DLTcollab sse2neon
 * translation layer (SSE->NEON). Replaces the legacy jratcliff shim, which is
 * kept as emmintrin_legacy.h for reference. */
#include "sse2neon.h"
#endif
