#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "bseq.h"
#include "minimap.h"
#include "mmpriv.h"
#include "getopt.h"

#define MM_VERSION "2.1-r337-dirty"

#ifdef __linux__
#include <sys/resource.h>
#include <sys/time.h>
void liftrlimit()
{
	struct rlimit r;
	getrlimit(RLIMIT_AS, &r);
	r.rlim_cur = r.rlim_max;
	setrlimit(RLIMIT_AS, &r);
}
#else
void liftrlimit() {}
#endif

static struct option long_options[] = {
	{ "bucket-bits",    required_argument, 0, 0 },
	{ "mb-size",        required_argument, 0, 'K' },
	{ "int-rname",      no_argument,       0, 0 },
	{ "no-kalloc",      no_argument,       0, 0 },
	{ "print-qname",    no_argument,       0, 0 },
	{ "no-self",        no_argument,       0, 0 },
	{ "print-seed",     no_argument,       0, 0 },
	{ "max-chain-skip", required_argument, 0, 0 },
	{ "min-dp-len",     required_argument, 0, 0 },
	{ "print-aln-seq",  no_argument,       0, 0 },
	{ "splice",         no_argument,       0, 0 },
	{ "cost-non-gt-ag", required_argument, 0, 0 },
	{ "no-sam-sq",      no_argument,       0, 0 },
	{ "help",           no_argument,       0, 'h' },
	{ "max-intron-len", required_argument, 0, 'G' },
	{ "version",        no_argument,       0, 'V' },
	{ "min-count",      required_argument, 0, 'n' },
	{ "min-chain-score",required_argument, 0, 'm' },
	{ "mask-level",     required_argument, 0, 'M' },
	{ "min-dp-score",   required_argument, 0, 's' },
	{ "sam",            no_argument,       0, 'a' },
	{ 0, 0, 0, 0}
};

static inline int64_t mm_parse_num(const char *str)
{
	double x;
	char *p;
	x = strtod(optarg, &p);
	if (*p == 'G' || *p == 'g') x *= 1e9;
	else if (*p == 'M' || *p == 'm') x *= 1e6;
	else if (*p == 'K' || *p == 'k') x *= 1e3;
	return (int64_t)(x + .499);
}

int main(int argc, char *argv[])
{
	mm_mapopt_t opt;
	int i, c, k = 15, w = -1, bucket_bits = MM_IDX_DEF_B, n_threads = 3, keep_name = 1, is_idx, is_hpc = 0, long_idx, idx_par_set = 0, max_intron_len = 0, n_idx_part = 0;
	int minibatch_size = 200000000;
	uint64_t batch_size = 4000000000ULL;
	mm_bseq_file_t *fp = 0;
	char *fnw = 0, *rg = 0, *s;
	FILE *fpr = 0, *fpw = 0, *fp_help = stderr;

	liftrlimit();
	mm_realtime0 = realtime();
	mm_mapopt_init(&opt);

	while ((c = getopt_long(argc, argv, "aSw:k:K:t:r:f:Vv:g:G:I:d:XT:s:x:Hcp:M:n:z:A:B:O:E:m:N:Qu:R:h", long_options, &long_idx)) >= 0) {
		if (c == 'w') w = atoi(optarg), idx_par_set = 1;
		else if (c == 'k') k = atoi(optarg), idx_par_set = 1;
		else if (c == 'H') is_hpc = 1, idx_par_set = 1;
		else if (c == 'd') fnw = optarg; // the above are indexing related options, except -I
		else if (c == 'r') opt.bw = (int)mm_parse_num(optarg);
		else if (c == 'f') opt.mid_occ_frac = atof(optarg);
		else if (c == 't') n_threads = atoi(optarg);
		else if (c == 'v') mm_verbose = atoi(optarg);
		else if (c == 'g') opt.max_gap = (int)mm_parse_num(optarg);
		else if (c == 'G') max_intron_len = (int)mm_parse_num(optarg);
		else if (c == 'N') opt.best_n = atoi(optarg);
		else if (c == 'p') opt.pri_ratio = atof(optarg);
		else if (c == 'M') opt.mask_level = atof(optarg);
		else if (c == 'c') opt.flag |= MM_F_OUT_CG | MM_F_CIGAR;
		else if (c == 'S') opt.flag |= MM_F_OUT_CS | MM_F_CIGAR;
		else if (c == 'X') opt.flag |= MM_F_AVA | MM_F_NO_SELF;
		else if (c == 'a') opt.flag |= MM_F_OUT_SAM | MM_F_CIGAR;
		else if (c == 'Q') opt.flag |= MM_F_NO_QUAL;
		else if (c == 'T') opt.sdust_thres = atoi(optarg);
		else if (c == 'n') opt.min_cnt = atoi(optarg);
		else if (c == 'm') opt.min_chain_score = atoi(optarg);
		else if (c == 'A') opt.a = atoi(optarg);
		else if (c == 'B') opt.b = atoi(optarg);
		else if (c == 'z') opt.zdrop = atoi(optarg);
		else if (c == 's') opt.min_dp_max = atoi(optarg);
		else if (c == 'I') batch_size = mm_parse_num(optarg);
		else if (c == 'K') minibatch_size = (int)mm_parse_num(optarg);
		else if (c == 'R') rg = optarg;
		else if (c == 'h') fp_help = stdout;
		else if (c == 0 && long_idx == 0) bucket_bits = atoi(optarg); // --bucket-bits
		else if (c == 0 && long_idx == 2) keep_name = 0; // --int-rname
		else if (c == 0 && long_idx == 3) mm_dbg_flag |= MM_DBG_NO_KALLOC; // --no-kalloc
		else if (c == 0 && long_idx == 4) mm_dbg_flag |= MM_DBG_PRINT_QNAME; // --print-qname
		else if (c == 0 && long_idx == 5) opt.flag |= MM_F_NO_SELF; // --no-self
		else if (c == 0 && long_idx == 6) mm_dbg_flag |= MM_DBG_PRINT_QNAME | MM_DBG_PRINT_SEED; // --print-seed
		else if (c == 0 && long_idx == 7) opt.max_chain_skip = atoi(optarg); // --max-chain-skip
		else if (c == 0 && long_idx == 8) opt.min_ksw_len = atoi(optarg); // --min-dp-len
		else if (c == 0 && long_idx == 9) mm_dbg_flag |= MM_DBG_PRINT_QNAME | MM_DBG_PRINT_ALN_SEQ; // --print-aln-seq
		else if (c == 0 && long_idx ==10) opt.flag |= MM_F_SPLICE; // --splice
		else if (c == 0 && long_idx ==11) opt.noncan = atoi(optarg); // --cost-non-gt-ag
		else if (c == 0 && long_idx ==12) opt.flag |= MM_F_NO_SAM_SQ; // --no-sam-sq
		else if (c == 'V') {
			puts(MM_VERSION);
			return 0;
		} else if (c == 'u') {
			if (*optarg == 'b') opt.flag |= MM_F_SPLICE_FOR|MM_F_SPLICE_REV;
			else if (*optarg == 'B') opt.flag |= MM_F_SPLICE_BOTH;
			else if (*optarg == 'f') opt.flag |= MM_F_SPLICE_FOR, opt.flag &= ~MM_F_SPLICE_REV;
			else if (*optarg == 'r') opt.flag |= MM_F_SPLICE_REV, opt.flag &= ~MM_F_SPLICE_FOR;
			else if (*optarg == 'n') opt.flag &= ~(MM_F_SPLICE_FOR|MM_F_SPLICE_REV);
			else {
				fprintf(stderr, "[E::%s] unrecognized cDNA direction\n", __func__);
				return 1;
			}
		} else if (c == 'O') {
			opt.q = opt.q2 = strtol(optarg, &s, 10);
			if (*s == ',') opt.q2 = strtol(s + 1, &s, 10);
		} else if (c == 'E') {
			opt.e = opt.e2 = strtol(optarg, &s, 10);
			if (*s == ',') opt.e2 = strtol(s + 1, &s, 10);
		} else if (c == 'x') {
			if (strcmp(optarg, "ava-ont") == 0) {
				opt.flag |= MM_F_AVA | MM_F_NO_SELF;
				opt.min_chain_score = 100, opt.pri_ratio = 0.0f, opt.max_gap = 10000, opt.max_chain_skip = 25;
				minibatch_size = 500000000;
				k = 15, w = 5;
			} else if (strcmp(optarg, "ava-pb") == 0) {
				opt.flag |= MM_F_AVA | MM_F_NO_SELF;
				opt.min_chain_score = 100, opt.pri_ratio = 0.0f, opt.max_gap = 10000, opt.max_chain_skip = 25;
				minibatch_size = 500000000;
				is_hpc = 1, k = 19, w = 5;
			} else if (strcmp(optarg, "map10k") == 0 || strcmp(optarg, "map-pb") == 0) {
				is_hpc = 1, k = 19;
			} else if (strcmp(optarg, "map-ont") == 0) {
				is_hpc = 0, k = 15;
			} else if (strcmp(optarg, "asm5") == 0) {
				k = 19, w = 19;
				opt.a = 1, opt.b = 19, opt.q = 39, opt.q2 = 81, opt.e = 3, opt.e2 = 1, opt.zdrop = 200;
				opt.min_dp_max = 200;
			} else if (strcmp(optarg, "asm10") == 0) {
				k = 19, w = 19;
				opt.a = 1, opt.b = 9, opt.q = 16, opt.q2 = 41, opt.e = 2, opt.e2 = 1, opt.zdrop = 200;
				opt.min_dp_max = 200;
			} else if (strcmp(optarg, "splice") == 0 || strcmp(optarg, "cdna") == 0) {
				k = 15, w = 5;
				opt.flag |= MM_F_SPLICE | MM_F_SPLICE_FOR | MM_F_SPLICE_REV;
				opt.max_gap = 2000, opt.max_gap_ref = opt.bw = 200000;
				opt.a = 1, opt.b = 2, opt.q = 2, opt.e = 1, opt.q2 = 32, opt.e2 = 0;
				opt.noncan = 5;
				opt.zdrop = 200;
			} else {
				fprintf(stderr, "[E::%s] unknown preset '%s'\n", __func__, optarg);
				return 1;
			}
		}
	}
	if (w < 0) w = (int)(.6666667 * k + .499);
	if ((opt.flag & MM_F_SPLICE) && max_intron_len > 0)
		opt.max_gap_ref = opt.bw = max_intron_len;

	if (argc == optind || fp_help == stdout) {
		fprintf(fp_help, "Usage: minimap2 [options] <target.fa>|<target.idx> [query.fa] [...]\n");
		fprintf(fp_help, "Options:\n");
		fprintf(fp_help, "  Indexing:\n");
		fprintf(fp_help, "    -H           use homopolymer-compressed k-mer\n");
		fprintf(fp_help, "    -k INT       k-mer size (no larger than 28) [%d]\n", k);
		fprintf(fp_help, "    -w INT       minizer window size [{-k}*2/3]\n");
		fprintf(fp_help, "    -I NUM       split index for every ~NUM input bases [4G]\n");
		fprintf(fp_help, "    -d FILE      dump index to FILE []\n");
		fprintf(fp_help, "  Mapping:\n");
		fprintf(fp_help, "    -f FLOAT     filter out top FLOAT fraction of repetitive minimizers [%g]\n", opt.mid_occ_frac);
		fprintf(fp_help, "    -g INT       stop chain enlongation if there are no minimizers in INT-bp [%d]\n", opt.max_gap);
		fprintf(fp_help, "    -r INT       bandwidth used in chaining and DP-based alignment [%d]\n", opt.bw);
		fprintf(fp_help, "    -n INT       minimal number of minimizers on a chain [%d]\n", opt.min_cnt);
		fprintf(fp_help, "    -m INT       minimal chaining score (matching bases minus log gap penalty) [%d]\n", opt.min_chain_score);
//		fprintf(fp_help, "    -T INT       SDUST threshold; 0 to disable SDUST [%d]\n", opt.sdust_thres); // TODO: this option is never used; might be buggy
		fprintf(fp_help, "    -X           skip self and dual mappings (for the all-vs-all mode)\n");
		fprintf(fp_help, "    -p FLOAT     min secondary-to-primary score ratio [%g]\n", opt.pri_ratio);
		fprintf(fp_help, "    -N INT       retain at most INT secondary alignments [%d]\n", opt.best_n);
		fprintf(fp_help, "    -G NUM       max intron length (only effective following -x splice) [200k]\n");
		fprintf(fp_help, "  Alignment:\n");
		fprintf(fp_help, "    -A INT       matching score [%d]\n", opt.a);
		fprintf(fp_help, "    -B INT       mismatch penalty [%d]\n", opt.b);
		fprintf(fp_help, "    -O INT[,INT] gap open penalty [%d,%d]\n", opt.q, opt.q2);
		fprintf(fp_help, "    -E INT[,INT] gap extension penalty; a k-long gap costs min{O1+k*E1,O2+k*E2} [%d,%d]\n", opt.e, opt.e2);
		fprintf(fp_help, "    -z INT       Z-drop score [%d]\n", opt.zdrop);
		fprintf(fp_help, "    -s INT       minimal peak DP alignment score [%d]\n", opt.min_dp_max);
		fprintf(fp_help, "    -u CHAR      how to find GT-AG. f:transcript strand, b:both strands, n:don't match GT-AG [n]\n");
		fprintf(fp_help, "  Input/Output:\n");
		fprintf(fp_help, "    -a           output in the SAM format (PAF by default)\n");
		fprintf(fp_help, "    -Q           don't output base quality in SAM\n");
		fprintf(fp_help, "    -R STR       SAM read group line in a format like '@RG\\tID:foo\\tSM:bar' []\n");
		fprintf(fp_help, "    -c           output CIGAR in PAF\n");
		fprintf(fp_help, "    -S           output the cs tag in PAF (cs encodes both query and ref sequences)\n");
		fprintf(fp_help, "    -t INT       number of threads [%d]\n", n_threads);
		fprintf(fp_help, "    -K NUM       minibatch size [200M]\n");
//		fprintf(fp_help, "    -v INT       verbose level [%d]\n", mm_verbose);
		fprintf(fp_help, "    --version    show version number\n");
		fprintf(fp_help, "  Preset:\n");
		fprintf(fp_help, "    -x STR       preset (recommended to be applied before other options) []\n");
		fprintf(fp_help, "                 map10k/map-pb: -Hk19 (PacBio/ONT vs reference mapping)\n");
		fprintf(fp_help, "                 map-ont: -k15 (slightly more sensitive than 'map10k' for ONT vs reference)\n");
		fprintf(fp_help, "                 asm5: -k19 -w19 -A1 -B19 -O39,81 -E3,1 -s200 -z200 (asm to ref mapping; break at 5%% div.)\n");
		fprintf(fp_help, "                 asm10: -k19 -w19 -A1 -B9 -O16,41 -E2,1 -s200 -z200 (asm to ref mapping; break at 10%% div.)\n");
		fprintf(fp_help, "                 ava-pb: -Hk19 -w5 -Xp0 -m100 -g10000 -K500m --max-chain-skip 25 (PacBio read overlap)\n");
		fprintf(fp_help, "                 ava-ont: -k15 -w5 -Xp0 -m100 -g10000 -K500m --max-chain-skip 25 (ONT read overlap)\n");
		fprintf(fp_help, "                 splice: long-read spliced alignment (see minimap2.1 for details)\n");
		fprintf(fp_help, "\nSee `man ./minimap2.1' for detailed description of command-line options.\n");
		return fp_help == stdout? 0 : 1;
	}

	is_idx = mm_idx_is_idx(argv[optind]);
	if (is_idx < 0) {
		fprintf(stderr, "[ERROR] failed to open file '%s'\n", argv[optind]);
		return 1;
	}
	if (!is_idx && fnw == 0 && argc - optind < 2) {
		fprintf(stderr, "[ERROR] missing input: please specify a query file to map or option -d to keep the index\n");
		return 1;
	}
	if (is_idx) fpr = fopen(argv[optind], "rb");
	else fp = mm_bseq_open(argv[optind]);
	if (fnw) fpw = fopen(fnw, "wb");
	if (opt.flag & MM_F_OUT_SAM)
		mm_write_sam_hdr_no_SQ(rg, MM_VERSION, argc, argv);
	for (;;) {
		mm_idx_t *mi;
		if (fpr) {
			mi = mm_idx_load(fpr);
			if (mi == 0) break;
			if (idx_par_set && mm_verbose >= 2 && (mi->k != k || mi->w != w || mi->is_hpc != is_hpc))
				fprintf(stderr, "[WARNING] \033[1;31mIndexing parameters on the command line (-k/-w/-H) overridden by parameters in the prebuilt index.\033[0m\n");
		} else {
			mi = mm_idx_gen(fp, w, k, bucket_bits, is_hpc, minibatch_size, n_threads, batch_size, keep_name);
		}
		if (mi == 0) break;
		++n_idx_part;
		if (mm_verbose >= 2 && n_idx_part > 1 && (opt.flag&MM_F_OUT_SAM) && !(opt.flag&MM_F_NO_SAM_SQ))
			fprintf(stderr, "[WARNING] \033[1;31mSAM output is malformated due to internal @SQ lines. Please add option --no-sam-sq or filter afterwards.\033[0m\n");
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
	if (fp)  mm_bseq_close(fp);

	fprintf(stderr, "[M::%s] Version: %s\n", __func__, MM_VERSION);
	fprintf(stderr, "[M::%s] CMD:", __func__);
	for (i = 0; i < argc; ++i)
		fprintf(stderr, " %s", argv[i]);
	fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec\n", __func__, realtime() - mm_realtime0, cputime());
	return 0;
}
