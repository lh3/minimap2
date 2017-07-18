// To compile:
//   gcc -g -O2 example.c libminimap2.a -lz

#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <zlib.h>
#include "minimap.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

int main(int argc, char *argv[])
{
	mm_verbose = 2; // disable message output to stderr

	if (argc < 3) {
		fprintf(stderr, "Usage: minimap2-lite <target.fa> <query.fa>\n");
		return 1;
	}
	
	// open query file for reading; you may use your favorite FASTA/Q parser
	gzFile f = gzopen(argv[2], "r");
	assert(f);
	kseq_t *ks = kseq_init(f);

	// create index for target; we are creating one index for all target sequence
	int n_threads = 4, w = 10, k = 15, is_hpc = 0;
	mm_idx_t *mi = mm_idx_build(argv[1], w, k, is_hpc, n_threads);
	assert(mi);

	// mapping
	mm_mapopt_t opt;
	mm_mapopt_init(&opt); // initialize mapping parameters
	mm_mapopt_update(&opt, mi); // this sets the maximum minimizer occurrence; TODO: set a better default in mm_mapopt_init()!
	opt.flag |= MM_F_CIGAR; // perform alignment
	mm_tbuf_t *tbuf = mm_tbuf_init(); // thread buffer; for multi-threading, allocate one tbuf for each thread
	while (kseq_read(ks) >= 0) { // each kseq_read() call reads one query sequence
		const mm_reg1_t *reg;
		int j, i, n_reg;
		// get all hits for the query
		reg = mm_map(mi, ks->seq.l, ks->seq.s, &n_reg, tbuf, &opt, 0);
		// traverse hits and print them out
		for (j = 0; j < n_reg; ++j) {
			const mm_reg1_t *r = &reg[j];
			assert(r->p); // with MM_F_CIGAR, this should not be NULL
			printf("%s\t%d\t%d\t%d\t%c\t", ks->name.s, ks->seq.l, r->qs, r->qe, "+-"[r->rev]);
			printf("%s\t%d\t%d\t%d\t%d\t%d\t%d\tcg:Z:", mi->seq[r->rid].name, mi->seq[r->rid].len, r->rs, r->re,
				   r->p->blen - r->p->n_ambi - r->p->n_diff, r->p->blen, r->mapq);
			for (i = 0; i < r->p->n_cigar; ++i) // IMPORTANT: this gives the CIGAR in the aligned regions. NO soft/hard clippings!
				printf("%d%c", r->p->cigar[i]>>4, "MIDSHN"[r->p->cigar[i]&0xf]);
			putchar('\n');
		}
	}
	mm_tbuf_destroy(tbuf);

	// deallocate index and close the query file
	mm_idx_destroy(mi);
	kseq_destroy(ks);
	gzclose(f);
	return 0;
}
