#ifndef _KALLOC_H_
#define _KALLOC_H_

#include <stdlib.h>

#define km_size(x) (*(((size_t*)(x))-1) * sizeof(size_t))

#ifdef __cplusplus
extern "C" {
#endif

void *kmalloc(void *km, size_t size);
void *krealloc(void *km, void *ptr, size_t size);
void *kcalloc(void *km, size_t count, size_t size);
void kfree(void *km, void *ptr);

void *km_init(void);
void km_destroy(void *km);

void km_stat(const void *km); // TODO: return numbers instead of print to stderr

#ifdef __cplusplus
}
#endif

#endif
