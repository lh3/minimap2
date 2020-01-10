#ifndef _KALLOC_H_
#define _KALLOC_H_

#include <stddef.h> /* for size_t */

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
	size_t capacity, available, n_blocks, n_cores, largest;
} km_stat_t;

void *kmalloc(void *km, size_t size);
void *krealloc(void *km, void *ptr, size_t size);
void *kcalloc(void *km, size_t count, size_t size);
void kfree(void *km, void *ptr);

void *km_init(void);
void *km_init2(void *km_par, size_t min_core_size);
void km_destroy(void *km);
void km_stat(const void *_km, km_stat_t *s);

#ifdef __cplusplus
}
#endif

#define KMALLOC(km, ptr, len) ((ptr) = (__typeof__(ptr))kmalloc((km), (len) * sizeof(*(ptr))))
#define KCALLOC(km, ptr, len) ((ptr) = (__typeof__(ptr))kcalloc((km), (len), sizeof(*(ptr))))
#define KREALLOC(km, ptr, len) ((ptr) = (__typeof__(ptr))krealloc((km), (ptr), (len) * sizeof(*(ptr))))

#define KEXPAND(km, a, m) do { \
		(m) = (m) >= 4? (m) + ((m)>>1) : 16; \
		KREALLOC((km), (a), (m)); \
	} while (0)

#endif
