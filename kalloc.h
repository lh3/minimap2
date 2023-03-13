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

#define KMALLOC(__type, km, ptr, len) ((ptr) = (__type*)kmalloc((km), (len) * sizeof(*(ptr))))
#define KCALLOC(__type, km, ptr, len) ((ptr) = (__type*)kcalloc((km), (len), sizeof(*(ptr))))
#define KREALLOC(__type, km, ptr, len) ((ptr) = (__type*)krealloc((km), (ptr), (len) * sizeof(*(ptr))))

#define KEXPAND(__type, km, a, m) do { \
		(m) = (m) >= 4? (m) + ((m)>>1) : 16; \
		KREALLOC(__type, (km), (a), (m)); \
	} while (0)

#ifndef klib_unused
#if (defined __clang__ && __clang_major__ >= 3) || (defined __GNUC__ && __GNUC__ >= 3)
#define klib_unused __attribute__ ((__unused__))
#else
#define klib_unused
#endif
#endif /* klib_unused */

#define KALLOC_POOL_INIT2(SCOPE, name, kmptype_t) \
	typedef struct { \
		size_t cnt, n, max; \
		kmptype_t **buf; \
		void *km; \
	} kmp_##name##_t; \
	SCOPE kmp_##name##_t *kmp_init_##name(void *km) { \
		kmp_##name##_t *mp; \
		KCALLOC(kmp_##name##_t, km, mp, 1); \
		mp->km = km; \
		return mp; \
	} \
	SCOPE void kmp_destroy_##name(kmp_##name##_t *mp) { \
		size_t k; \
		for (k = 0; k < mp->n; ++k) kfree(mp->km, mp->buf[k]); \
		kfree(mp->km, mp->buf); kfree(mp->km, mp); \
	} \
	SCOPE kmptype_t *kmp_alloc_##name(kmp_##name##_t *mp) { \
		++mp->cnt; \
		if (mp->n == 0) return (kmptype_t*)kcalloc(mp->km, 1, sizeof(kmptype_t)); \
		return mp->buf[--mp->n]; \
	} \
	SCOPE void kmp_free_##name(kmp_##name##_t *mp, kmptype_t *p) { \
		--mp->cnt; \
		if (mp->n == mp->max) KEXPAND(kmptype_t*, mp->km, mp->buf, mp->max); \
		mp->buf[mp->n++] = p; \
	}

#define KALLOC_POOL_INIT(name, kmptype_t) \
	KALLOC_POOL_INIT2(static inline klib_unused, name, kmptype_t)

#endif
