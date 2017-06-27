#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include "kalloc.h"

/* The whole thing is: ("@" for the kheader_t of the block, "-" for free
 * memory, and "+" for allocated memory. One char for one unit.)
 *                        
 *           This region is core 1.                             This region is core 2.
 *
 *   @-------@++++++@++++++++++++@------------           @----------@++++++++++++@+++++++@------------
 *   |                           |                       |                               |
 *   p=p->ptr->ptr->ptr->ptr     p->ptr             p->ptr->ptr                p->ptr->ptr->ptr
 */

#define PTR(p) ((size_t*)((size_t*)p)[1])

typedef struct _allocated_t {
	struct _allocated_t *next;
	size_t *ptr;
} allocated_t;

typedef struct {
	size_t base[2], *loop_head;
	allocated_t list_head, *list_tail;
	size_t total_allocated;
} kmem_t;

void *km_init()
{
	return calloc(1, sizeof(kmem_t));
}

static void kerror(const char *s)
{
	fprintf(stderr, "%s\n", s);
	exit(1);
}

static size_t *morecore(kmem_t *km, size_t nu)
{
	size_t rnu, *up;

	rnu = (nu + 0xfffff) & (~(size_t)0xfffff);
	up = (size_t*)malloc(rnu * sizeof(size_t));
	if (!up) { /* fail to allocate memory */
		km_stat(km);
		fprintf(stderr, "[morecore] %lu bytes requested but not available.\n", rnu * sizeof(size_t));
		exit(1);
	}
	/* put the pointer in km->list_head */
	if (km->list_tail == 0) km->list_tail = &km->list_head;
	km->list_tail->ptr = up;
	km->list_tail->next = (allocated_t*)calloc(1, sizeof(allocated_t));
	km->list_tail = km->list_tail->next;

	km->total_allocated += rnu * sizeof(size_t);
	*up = rnu; /* the size of the current block, and in this case the block is the same as the new core */
	kfree(km, up + 1); /* initialize the new "core" */
	return km->loop_head;
}

void km_destroy(void *_km)
{
	kmem_t *km = (kmem_t*)_km;
	allocated_t *p, *q;
	if (km == 0) return;
	p = &km->list_head;
	do {
		q = p->next;
		free(p->ptr);
		if (p != &km->list_head) free(p);
		p = q;
	} while (p && p->next);
	if (p != &km->list_head) free(p);
	free(km);
}

void kfree(void *_km, void *ap)
{
	size_t *p, *q;
	kmem_t *km = (kmem_t*)_km;
	
	if (!ap) return;
	if (km == 0) {
		free(ap);
		return;
	}
	p = (size_t*)ap - 1; /* *p is the size of the current block */
	/* Find the pointer that points to the block to be freed. The following loop can stop on two conditions:
	 *
	 * a) "p>q && p<q->ptr": @------@++++++++@+++++++@-------    @---------------@+++++++@-------
	 *    (can also be in    |      |                |        -> |                       |
	 *     two cores)        q      p           q->ptr           q                  q->ptr
	 *
	 *                       @--------    @+++++++++@--------    @--------    @------------------
	 *                       |            |         |         -> |            |
	 *                       q            p    q->ptr            q       q->ptr
	 *
	 * b) "q>=q->ptr && (p>q || p<q->ptr)":  @-------@+++++   @--------@+++++++     @-------@+++++   @----------------
	 *                                       |                |        |         -> |                |
	 *                                  q->ptr                q        p       q->ptr                q
	 *
	 *                                       @+++++++@-----   @++++++++@-------     @-------------   @++++++++@-------
	 *                                       |       |                 |         -> |                         |
	 *                                       p  q->ptr                 q       q->ptr                         q
	 */
	for (q = km->loop_head; !(p > q && p < PTR(q)); q = PTR(q))
		if (q >= PTR(q) && (p > q || p < PTR(q))) break;
	if (p + (*p) == PTR(q)) { /* two adjacent blocks, merge p and q->ptr (the 2nd and 4th cases) */
		*p += *PTR(q); /* this is the new q->ptr size */
		p[1] = (size_t)PTR(PTR(q)); /* this is the new q->ptr->ptr */
		/* p is actually the new q->ptr. The actual change happens a few lines below. */
	} else if (p + (*p) > PTR(q) && PTR(q) >= p) { /* the end of the allocated block is in the next free block */
		kerror("[kfree] The end of the allocated block enters a free block.");
	} else p[1] = (size_t)PTR(q); /* backup q->ptr */

	if (q + (*q) == p) { /* two adjacent blocks, merge q and p (the other two cases) */
		*q += *p;
		q[1] = (size_t)PTR(p);
		km->loop_head = q;
	} else if (q + (*q) > p && p >= q) { /* the end of a free block in the allocated block */
		kerror("[kfree] The end of a free block enters the allocated block.");
	} else km->loop_head = p, q[1] = (size_t)p; /* in two cores, cannot be merged */
}

void *krealloc(void *_km, void *ap, size_t n_bytes)
{
	kmem_t *km = (kmem_t*)_km;
	size_t n_units, *p, *q;

	if (n_bytes == 0) {
		kfree(km, ap); return 0;
	}
	if (km == 0) return realloc(ap, n_bytes);
	if (!ap) return kmalloc(km, n_bytes);
	n_units = 1 + (n_bytes + sizeof(size_t) - 1) / sizeof(size_t);
	p = (size_t*)ap - 1;
	if (*p >= n_units) return ap; /* TODO: this prevents shrinking */
	q = (size_t*)kmalloc(km, n_bytes);
	memcpy(q, ap, (*p - 1) * sizeof(size_t));
	kfree(km, ap);
	return q;
}

void *kmalloc(void *_km, size_t n_bytes)
{
	kmem_t *km = (kmem_t*)_km;
	size_t n_units, *p, *q;

	if (n_bytes == 0) return 0;
	if (km == 0) return malloc(n_bytes);
	/* "n_units" means the number of units. The size of one unit equals to sizeof(kheader_t).
	 * "1" is the kheader_t of a block, which is always required. */
	n_units = 1 + (n_bytes + sizeof(size_t) - 1) / sizeof(size_t);
	if (n_units&1) ++n_units; /* make n_units an even number, or it will segfault if only one unit remains */

	if (!(q = km->loop_head)) { /* the first time when kmalloc() is called, intialization */
		km->base[1] = (size_t)(km->loop_head = q = km->base); *q = 0;
	}
	for (p = PTR(q);; q = p, p = PTR(p)) { /* search for a suitable block */
		if (*p >= n_units) { /* p->size if the size of current block. This line means the current block is large enough. */
			if (*p == n_units) q[1] = (size_t)PTR(p); /* no need to split the block */
			else { /* split the block */
				/* memory is allocated at the end of the block */
				*p -= n_units; /* reduce the size of the free block */
				p += *p; /* skip to the kheader_t of the allocated block */
				*p = n_units; /* set the size */
			}
			km->loop_head = q; /* set the end of chain */
			return p + 1; /* skip the kheader_t */
		}
		if (p == km->loop_head) { /* then ask for more "cores" */
			if ((p = morecore(km, n_units)) == 0) return 0;
		}
	}
}

void *kcalloc(void *_km, size_t count, size_t size)
{
	kmem_t *km = (kmem_t*)_km;
	void *p;
	if (size == 0 || count == 0) return 0;
	if (km == 0) return calloc(count, size);
	p = kmalloc(km, count * size);
	memset(p, 0, count * size);
	return p;
}

void km_stat(const void *_km)
{
	kmem_t *km = (kmem_t*)_km;
	unsigned n_blocks, n_units;
	size_t max_block = 0, *p, *q;
	float frag;

	if (km == 0 || !(p = km->loop_head)) return;
	n_blocks = n_units = 0;
	do {
		q = PTR(p);
		if (*p > max_block) max_block = *p;
		n_units += *p;
		if (p + (*p) > q && q > p)
			kerror("[kr_stat] The end of a free block enters another free block.");
		p = q;
		++n_blocks;
	} while (p != km->loop_head);
	
	--n_blocks;
	frag = 1.0/1024.0 * n_units * sizeof(size_t) / n_blocks;
	fprintf(stderr, "[kr_stat] tot=%lu, free=%lu, n_block=%u, max_block=%lu, frag_len=%.3fK\n",
			km->total_allocated, n_units * sizeof(size_t), n_blocks, max_block * sizeof(size_t), frag);
}
