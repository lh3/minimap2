#include <getopt.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/resource.h>
#include <sys/time.h>
#include <fcntl.h>
#include "bseq.h"
#include "minimap.h"
#include "mmpriv.h"

#define MM_VERSION "2.0-r88-pre"

void liftrlimit()
{
#ifdef __linux__
	struct rlimit r;
	getrlimit(RLIMIT_AS, &r);
	r.rlim_cur = r.rlim_max;
	setrlimit(RLIMIT_AS, &r);
#endif
}

static int test_idx(const char *fn)
{
	int fd, is_idx = 0;
	off_t ret;
	char magic[4];

	if (strcmp(fn, "-") == 0) return 0; // read from pipe; not an index
	fd = open(fn, O_RDONLY);
	if (fd < 0) return -1; // error
	if ((ret = lseek(fd, 0, SEEK_END)) >= 4) {
		lseek(fd, 0, SEEK_SET);
		read(fd, magic, 4);
		if (strncmp(magic, MM_IDX_MAGIC, 4) == 0)
			is_idx = 1;
	}
	close(fd);
	return is_idx;
}

static struct option long_options[] = {
	{ "bucket-bits",    required_argument, 0, 0 },
	{ "mb-size",        required_argument, 0, 0 },
	{ "int-rname",      no_argument,       0, 0 },
	{ "version",        no_argument,       0, 'V' },
	{ 0, 0, 0, 0}
};

int main(int argc, char *argv[])
{
	mm_mapopt_t opt;
	int i, c, k = 17, w = -1, bucket_bits = MM_IDX_DEF_B, n_threads = 3, keep_name = 1, is_idx, is_hpc = 0, long_idx;
	int minibatch_size = 100000000;
	uint64_t batch_size = 4000000000ULL;
	bseq_file_t *fp = 0;
	char *fnw = 0;
	FILE *fpr = 0, *fpw = 0;

	liftrlimit();
	mm_realtime0 = realtime();
	mm_mapopt_init(&opt);

	while ((c = getopt_long(argc, argv, "w:k:t:r:f:Vv:g:I:d:ST:s:x:Hcp:m:z:F:", long_options, &long_idx)) >= 0) {
		if (c == 'w') w = atoi(optarg);
		else if (c == 'k') k = atoi(optarg);
		else if (c == 'H') is_hpc = 1;
		else if (c == 'd') fnw = optarg; // the above are indexing related options, except -I
		else if (c == 'r') opt.bw = atoi(optarg);
		else if (c == 'f') opt.mid_occ_frac = atof(optarg);
		else if (c == 't') n_threads = atoi(optarg);
		else if (c == 'v') mm_verbose = atoi(optarg);
		else if (c == 'g') opt.max_gap = atoi(optarg);
		else if (c == 'p') opt.pri_ratio = atof(optarg);
		else if (c == 'm') opt.mask_level = atof(optarg);
		else if (c == 'c') opt.flag |= MM_F_CIGAR;
		else if (c == 'S') opt.flag |= MM_F_AVA | MM_F_NO_SELF;
		else if (c == 'T') opt.sdust_thres = atoi(optarg);
		else if (c == 's') opt.min_score = atoi(optarg);
		else if (c == 'z') opt.zdrop = atoi(optarg);
		else if (c == 0 && long_idx == 0) bucket_bits = atoi(optarg); // bucket-bits
		else if (c == 0 && long_idx == 2) keep_name = 0; // int-rname
		else if (c == 'V') {
			puts(MM_VERSION);
			return 0;
		} else if (c == 'I' || (c == 0 && long_idx == 1)) {
			double x;
			char *p;
			x = strtod(optarg, &p);
			if (*p == 'G' || *p == 'g') x *= 1e9;
			else if (*p == 'M' || *p == 'm') x *= 1e6;
			else if (*p == 'K' || *p == 'k') x *= 1e3;
			if (c == 'I') batch_size = (uint64_t)(x + .499);
			else minibatch_size = (uint64_t)(x + .499);
		} else if (c == 'F') {
			if (strcmp(optarg, "sam") == 0) opt.flag |= MM_F_OUT_SAM | MM_F_CIGAR;
			else if (strcmp(optarg, "paf") == 0) opt.flag &= ~MM_F_OUT_SAM;
			else {
				fprintf(stderr, "[E::%s] unknown output format '%s'\n", __func__, optarg);
				return 0;
			}
		} else if (c == 'x') {
			if (strcmp(optarg, "ava10k") == 0) {
				opt.flag |= MM_F_AVA | MM_F_NO_SELF;
				opt.min_score = 100;
				w = 5;
			}
		}
	}
	if (w < 0) w = (int)(.6666667 * k + .499);

	if (argc == optind) {
		fprintf(stderr, "Usage: minimap2 [options] <target.fa>|<target.idx> [query.fa] [...]\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  Indexing:\n");
		fprintf(stderr, "    -H         use homopolymer-compressed k-mer\n");
		fprintf(stderr, "    -k INT     k-mer size (no larger than 28) [%d]\n", k);
		fprintf(stderr, "    -w INT     minizer window size [{-k}*2/3]\n");
		fprintf(stderr, "    -I NUM     split index for every ~NUM input bases [4G]\n");
		fprintf(stderr, "    -d FILE    dump index to FILE []\n");
		fprintf(stderr, "  Mapping:\n");
		fprintf(stderr, "    -f FLOAT   filter out top FLOAT fraction of repetitive minimizers [%.3f]\n", opt.mid_occ_frac);
		fprintf(stderr, "    -r INT     bandwidth [%d]\n", opt.bw);
		fprintf(stderr, "    -s INT     min score [%d]\n", opt.min_score);
		fprintf(stderr, "    -g INT     split a mapping if there is a gap longer than INT [%d]\n", opt.max_gap);
		fprintf(stderr, "    -T INT     SDUST threshold; 0 to disable SDUST [%d]\n", opt.sdust_thres);
		fprintf(stderr, "    -S         skip self and dual mappings\n");
		fprintf(stderr, "    -p FLOAT   threshold to output a mapping [%g]\n", opt.pri_ratio);
		fprintf(stderr, "    -z INT     Z-drop score [%d]\n", opt.zdrop);
		fprintf(stderr, "    -x STR     preset (recommended to be applied before other options) []\n");
		fprintf(stderr, "               ava10k: -Sw5 -L100 -m0 (PacBio/ONT all-vs-all read mapping)\n");
		fprintf(stderr, "  Input/Output:\n");
		fprintf(stderr, "    -F STR     output format: sam or paf [paf]\n");
		fprintf(stderr, "    -c         output CIGAR in PAF\n");
		fprintf(stderr, "    -t INT     number of threads [%d]\n", n_threads);
//		fprintf(stderr, "    -v INT     verbose level [%d]\n", mm_verbose);
		fprintf(stderr, "    -V         show version number\n");
		fprintf(stderr, "\nSee minimap2.1 for detailed description of the command-line options.\n");
		return 1;
	}

	is_idx = test_idx(argv[optind]);
	if (is_idx < 0) {
		fprintf(stderr, "[E::%s] failed to open file '%s'\n", __func__, argv[optind]);
		return 1;
	}
	if (is_idx) fpr = fopen(argv[optind], "rb");
	else fp = bseq_open(argv[optind]);
	if (fnw) fpw = fopen(fnw, "wb");
	for (;;) {
		mm_idx_t *mi = 0;
		if (fpr) mi = mm_idx_load(fpr);
		else if (!bseq_eof(fp))
			mi = mm_idx_gen(fp, w, k, bucket_bits, is_hpc, minibatch_size, n_threads, batch_size, keep_name);
		if (mi == 0) break;
		if (mm_verbose >= 3)
			fprintf(stderr, "[M::%s::%.3f*%.2f] loaded/built the index for %d target sequence(s)\n",
					__func__, realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0), mi->n_seq);
		if (fpw) {
			mm_idx_dump(fpw, mi);
			if (mm_verbose >= 3)
				fprintf(stderr, "[M::%s::%.3f*%.2f] dumpped the (partial) index to disk\n", __func__, realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0));
		}
		if (argc != optind + 1) mm_mapopt_update(&opt, mi);
		if (mm_verbose >= 3) mm_idx_stat(mi);
		for (i = optind + 1; i < argc; ++i)
			mm_map_file(mi, argv[i], &opt, n_threads, minibatch_size);
		mm_idx_destroy(mi);
	}
	if (fpw) fclose(fpw);
	if (fpr) fclose(fpr);
	if (fp)  bseq_close(fp);

	fprintf(stderr, "[M::%s] Version: %s\n", __func__, MM_VERSION);
	fprintf(stderr, "[M::%s] CMD:", __func__);
	for (i = 0; i < argc; ++i)
		fprintf(stderr, " %s", argv[i]);
	fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec\n", __func__, realtime() - mm_realtime0, cputime());
	return 0;
}
