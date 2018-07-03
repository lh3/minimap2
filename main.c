#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "bseq.h"
#include "minimap.h"
#include "mmpriv.h"
#ifdef HAVE_GETOPT
#include <getopt.h>
#else
#include "getopt.h"
#endif

#define MM_VERSION "2.11-r797"

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
	{ "seed",           required_argument, 0, 0 },
	{ "no-kalloc",      no_argument,       0, 0 },
	{ "print-qname",    no_argument,       0, 0 },
	{ "no-self",        no_argument,       0, 'D' },
	{ "print-seeds",    no_argument,       0, 0 },
	{ "max-chain-skip", required_argument, 0, 0 },
	{ "min-dp-len",     required_argument, 0, 0 },
	{ "print-aln-seq",  no_argument,       0, 0 },
	{ "splice",         no_argument,       0, 0 },
	{ "cost-non-gt-ag", required_argument, 0, 'C' },
	{ "no-long-join",   no_argument,       0, 0 },
	{ "sr",             no_argument,       0, 0 },
	{ "frag",           required_argument, 0, 0 },
	{ "secondary",      required_argument, 0, 0 },
	{ "cs",             optional_argument, 0, 0 },
	{ "end-bonus",      required_argument, 0, 0 },
	{ "no-pairing",     no_argument,       0, 0 },
	{ "splice-flank",   required_argument, 0, 0 },
	{ "idx-no-seq",     no_argument,       0, 0 },
	{ "end-seed-pen",   required_argument, 0, 0 },   // 21
	{ "for-only",       no_argument,       0, 0 },   // 22
	{ "rev-only",       no_argument,       0, 0 },   // 23
	{ "heap-sort",      required_argument, 0, 0 },   // 24
	{ "all-chain",      no_argument,       0, 'P' },
	{ "dual",           required_argument, 0, 0 },   // 26
	{ "max-clip-ratio", required_argument, 0, 0 },   // 27
	{ "min-occ-floor",  required_argument, 0, 0 },   // 28
	{ "MD",             no_argument,       0, 0 },   // 29
	{ "lj-min-ratio",   required_argument, 0, 0 },   // 30
	{ "score-N",        required_argument, 0, 0 },   // 31
	{ "eqx",            no_argument,       0, 0 },   // 32
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
	x = strtod(str, &p);
	if (*p == 'G' || *p == 'g') x *= 1e9;
	else if (*p == 'M' || *p == 'm') x *= 1e6;
	else if (*p == 'K' || *p == 'k') x *= 1e3;
	return (int64_t)(x + .499);
}

static inline void yes_or_no(mm_mapopt_t *opt, int flag, int long_idx, const char *arg, int yes_to_set)
{
	if (yes_to_set) {
		if (strcmp(arg, "yes") == 0 || strcmp(arg, "y") == 0) opt->flag |= flag;
		else if (strcmp(arg, "no") == 0 || strcmp(arg, "n") == 0) opt->flag &= ~flag;
		else fprintf(stderr, "[WARNING]\033[1;31m option '--%s' only accepts 'yes' or 'no'.\033[0m\n", long_options[long_idx].name);
	} else {
		if (strcmp(arg, "yes") == 0 || strcmp(arg, "y") == 0) opt->flag &= ~flag;
		else if (strcmp(arg, "no") == 0 || strcmp(arg, "n") == 0) opt->flag |= flag;
		else fprintf(stderr, "[WARNING]\033[1;31m option '--%s' only accepts 'yes' or 'no'.\033[0m\n", long_options[long_idx].name);
	}
}

int main(int argc, char *argv[])
{
	const char *opt_str = "2aSDw:k:K:t:r:f:Vv:g:G:I:d:XT:s:x:Hcp:M:n:z:A:B:O:E:m:N:Qu:R:hF:LC:yY";
	mm_mapopt_t opt;
	mm_idxopt_t ipt;
	int i, c, n_threads = 3, long_idx;
	char *fnw = 0, *rg = 0, *s;
	FILE *fp_help = stderr;
	mm_idx_reader_t *idx_rdr;
	mm_idx_t *mi;

	mm_verbose = 3;
	liftrlimit();
	mm_realtime0 = realtime();
	mm_set_opt(0, &ipt, &opt);

	while ((c = getopt_long(argc, argv, opt_str, long_options, &long_idx)) >= 0) // apply option -x/preset first
		if (c == 'x') {
			if (mm_set_opt(optarg, &ipt, &opt) < 0) {
				fprintf(stderr, "[ERROR] unknown preset '%s'\n", optarg);
				return 1;
			}
			break;
		}
	optind = 0; // for musl getopt, optind=0 has the same effect as optreset=1; older libc doesn't have optreset

	while ((c = getopt_long(argc, argv, opt_str, long_options, &long_idx)) >= 0) {
		if (c == 'w') ipt.w = atoi(optarg);
		else if (c == 'k') ipt.k = atoi(optarg);
		else if (c == 'H') ipt.flag |= MM_I_HPC;
		else if (c == 'd') fnw = optarg; // the above are indexing related options, except -I
		else if (c == 'r') opt.bw = (int)mm_parse_num(optarg);
		else if (c == 't') n_threads = atoi(optarg);
		else if (c == 'v') mm_verbose = atoi(optarg);
		else if (c == 'g') opt.max_gap = (int)mm_parse_num(optarg);
		else if (c == 'G') mm_mapopt_max_intron_len(&opt, (int)mm_parse_num(optarg));
		else if (c == 'F') opt.max_frag_len = (int)mm_parse_num(optarg);
		else if (c == 'N') opt.best_n = atoi(optarg);
		else if (c == 'p') opt.pri_ratio = atof(optarg);
		else if (c == 'M') opt.mask_level = atof(optarg);
		else if (c == 'c') opt.flag |= MM_F_OUT_CG | MM_F_CIGAR;
		else if (c == 'D') opt.flag |= MM_F_NO_DIAG;
		else if (c == 'P') opt.flag |= MM_F_ALL_CHAINS;
		else if (c == 'X') opt.flag |= MM_F_ALL_CHAINS | MM_F_NO_DIAG | MM_F_NO_DUAL | MM_F_NO_LJOIN; // -D -P --no-long-join --dual=no
		else if (c == 'a') opt.flag |= MM_F_OUT_SAM | MM_F_CIGAR;
		else if (c == 'Q') opt.flag |= MM_F_NO_QUAL;
		else if (c == 'Y') opt.flag |= MM_F_SOFTCLIP;
		else if (c == 'L') opt.flag |= MM_F_LONG_CIGAR;
		else if (c == 'y') opt.flag |= MM_F_COPY_COMMENT;
		else if (c == 'T') opt.sdust_thres = atoi(optarg);
		else if (c == 'n') opt.min_cnt = atoi(optarg);
		else if (c == 'm') opt.min_chain_score = atoi(optarg);
		else if (c == 'A') opt.a = atoi(optarg);
		else if (c == 'B') opt.b = atoi(optarg);
		else if (c == 's') opt.min_dp_max = atoi(optarg);
		else if (c == 'C') opt.noncan = atoi(optarg);
		else if (c == 'I') ipt.batch_size = mm_parse_num(optarg);
		else if (c == 'K') opt.mini_batch_size = (int)mm_parse_num(optarg);
		else if (c == 'R') rg = optarg;
		else if (c == 'h') fp_help = stdout;
		else if (c == '2') opt.flag |= MM_F_2_IO_THREADS;
		else if (c == 0 && long_idx == 0) ipt.bucket_bits = atoi(optarg); // --bucket-bits
		else if (c == 0 && long_idx == 2) opt.seed = atoi(optarg); // --seed
		else if (c == 0 && long_idx == 3) mm_dbg_flag |= MM_DBG_NO_KALLOC; // --no-kalloc
		else if (c == 0 && long_idx == 4) mm_dbg_flag |= MM_DBG_PRINT_QNAME; // --print-qname
		else if (c == 0 && long_idx == 6) mm_dbg_flag |= MM_DBG_PRINT_QNAME | MM_DBG_PRINT_SEED, n_threads = 1; // --print-seed
		else if (c == 0 && long_idx == 7) opt.max_chain_skip = atoi(optarg); // --max-chain-skip
		else if (c == 0 && long_idx == 8) opt.min_ksw_len = atoi(optarg); // --min-dp-len
		else if (c == 0 && long_idx == 9) mm_dbg_flag |= MM_DBG_PRINT_QNAME | MM_DBG_PRINT_ALN_SEQ, n_threads = 1; // --print-aln-seq
		else if (c == 0 && long_idx ==10) opt.flag |= MM_F_SPLICE; // --splice
		else if (c == 0 && long_idx ==12) opt.flag |= MM_F_NO_LJOIN; // --no-long-join
		else if (c == 0 && long_idx ==13) opt.flag |= MM_F_SR; // --sr
		else if (c == 0 && long_idx ==17) opt.end_bonus = atoi(optarg); // --end-bonus
		else if (c == 0 && long_idx ==18) opt.flag |= MM_F_INDEPEND_SEG; // --no-pairing
		else if (c == 0 && long_idx ==20) ipt.flag |= MM_I_NO_SEQ; // --idx-no-seq
		else if (c == 0 && long_idx ==21) opt.anchor_ext_shift = atoi(optarg); // --end-seed-pen
		else if (c == 0 && long_idx ==22) opt.flag |= MM_F_FOR_ONLY; // --for-only
		else if (c == 0 && long_idx ==23) opt.flag |= MM_F_REV_ONLY; // --rev-only
		else if (c == 0 && long_idx ==27) opt.max_clip_ratio = atof(optarg); // --max-clip-ratio
		else if (c == 0 && long_idx ==28) opt.min_mid_occ = atoi(optarg); // --min-occ-floor
		else if (c == 0 && long_idx ==29) opt.flag |= MM_F_OUT_MD; // --MD
		else if (c == 0 && long_idx ==30) opt.min_join_flank_ratio = atof(optarg); // --lj-min-ratio
		else if (c == 0 && long_idx ==31) opt.sc_ambi = atoi(optarg); // --score-N
		else if (c == 0 && long_idx ==32) opt.flag |= MM_F_EQX; // --eqx
		else if (c == 0 && long_idx == 14) { // --frag
			yes_or_no(&opt, MM_F_FRAG_MODE, long_idx, optarg, 1);
		} else if (c == 0 && long_idx == 15) { // --secondary
			yes_or_no(&opt, MM_F_NO_PRINT_2ND, long_idx, optarg, 0);
		} else if (c == 0 && long_idx == 16) { // --cs
			opt.flag |= MM_F_OUT_CS | MM_F_CIGAR;
			if (optarg == 0 || strcmp(optarg, "short") == 0) {
				opt.flag &= ~MM_F_OUT_CS_LONG;
			} else if (strcmp(optarg, "long") == 0) {
				opt.flag |= MM_F_OUT_CS_LONG;
			} else if (strcmp(optarg, "none") == 0) {
				opt.flag &= ~MM_F_OUT_CS;
			} else if (mm_verbose >= 2) {
				fprintf(stderr, "[WARNING]\033[1;31m --cs only takes 'short' or 'long'. Invalid values are assumed to be 'short'.\033[0m\n");
			}
		} else if (c == 0 && long_idx == 19) { // --splice-flank
			yes_or_no(&opt, MM_F_SPLICE_FLANK, long_idx, optarg, 1);
		} else if (c == 0 && long_idx == 24) { // --heap-sort
			yes_or_no(&opt, MM_F_HEAP_SORT, long_idx, optarg, 1);
		} else if (c == 0 && long_idx == 26) { // --dual
			yes_or_no(&opt, MM_F_NO_DUAL, long_idx, optarg, 0);
		} else if (c == 'S') {
			opt.flag |= MM_F_OUT_CS | MM_F_CIGAR | MM_F_OUT_CS_LONG;
			if (mm_verbose >= 2)
				fprintf(stderr, "[WARNING]\033[1;31m option -S is deprecated and may be removed in future. Please use --cs=long instead.\033[0m\n");
		} else if (c == 'V') {
			puts(MM_VERSION);
			return 0;
		} else if (c == 'f') {
			double x;
			char *p;
			x = strtod(optarg, &p);
			if (x < 1.0) opt.mid_occ_frac = x, opt.mid_occ = 0;
			else opt.mid_occ = (int)(x + .499);
			if (*p == ',') opt.max_occ = (int)(strtod(p+1, &p) + .499);
		} else if (c == 'u') {
			if (*optarg == 'b') opt.flag |= MM_F_SPLICE_FOR|MM_F_SPLICE_REV; // both strands
			else if (*optarg == 'f') opt.flag |= MM_F_SPLICE_FOR, opt.flag &= ~MM_F_SPLICE_REV; // match GT-AG
			else if (*optarg == 'r') opt.flag |= MM_F_SPLICE_REV, opt.flag &= ~MM_F_SPLICE_FOR; // match CT-AC (reverse complement of GT-AG)
			else if (*optarg == 'n') opt.flag &= ~(MM_F_SPLICE_FOR|MM_F_SPLICE_REV); // don't try to match the GT-AG signal
			else {
				fprintf(stderr, "[ERROR]\033[1;31m unrecognized cDNA direction\033[0m\n");
				return 1;
			}
		} else if (c == 'z') {
			opt.zdrop = opt.zdrop_inv = strtol(optarg, &s, 10);
			if (*s == ',') opt.zdrop_inv = strtol(s + 1, &s, 10);
		} else if (c == 'O') {
			opt.q = opt.q2 = strtol(optarg, &s, 10);
			if (*s == ',') opt.q2 = strtol(s + 1, &s, 10);
		} else if (c == 'E') {
			opt.e = opt.e2 = strtol(optarg, &s, 10);
			if (*s == ',') opt.e2 = strtol(s + 1, &s, 10);
		}
	}
	if ((opt.flag & MM_F_SPLICE) && (opt.flag & MM_F_FRAG_MODE)) {
		fprintf(stderr, "[ERROR]\033[1;31m --splice and --frag should not be specified at the same time.\033[0m\n");
		return 1;
	}
	if (!fnw && !(opt.flag&MM_F_CIGAR))
		ipt.flag |= MM_I_NO_SEQ;
	if (mm_check_opt(&ipt, &opt) < 0)
		return 1;

	if (argc == optind || fp_help == stdout) {
		fprintf(fp_help, "Usage: minimap2 [options] <target.fa>|<target.idx> [query.fa] [...]\n");
		fprintf(fp_help, "Options:\n");
		fprintf(fp_help, "  Indexing:\n");
		fprintf(fp_help, "    -H           use homopolymer-compressed k-mer (preferrable for PacBio)\n");
		fprintf(fp_help, "    -k INT       k-mer size (no larger than 28) [%d]\n", ipt.k);
		fprintf(fp_help, "    -w INT       minizer window size [%d]\n", ipt.w);
		fprintf(fp_help, "    -I NUM       split index for every ~NUM input bases [4G]\n");
		fprintf(fp_help, "    -d FILE      dump index to FILE []\n");
		fprintf(fp_help, "  Mapping:\n");
		fprintf(fp_help, "    -f FLOAT     filter out top FLOAT fraction of repetitive minimizers [%g]\n", opt.mid_occ_frac);
		fprintf(fp_help, "    -g NUM       stop chain enlongation if there are no minimizers in INT-bp [%d]\n", opt.max_gap);
		fprintf(fp_help, "    -G NUM       max intron length (effective with -xsplice; changing -r) [200k]\n");
		fprintf(fp_help, "    -F NUM       max fragment length (effective with -xsr or in the fragment mode) [800]\n");
		fprintf(fp_help, "    -r NUM       bandwidth used in chaining and DP-based alignment [%d]\n", opt.bw);
		fprintf(fp_help, "    -n INT       minimal number of minimizers on a chain [%d]\n", opt.min_cnt);
		fprintf(fp_help, "    -m INT       minimal chaining score (matching bases minus log gap penalty) [%d]\n", opt.min_chain_score);
//		fprintf(fp_help, "    -T INT       SDUST threshold; 0 to disable SDUST [%d]\n", opt.sdust_thres); // TODO: this option is never used; might be buggy
		fprintf(fp_help, "    -X           skip self and dual mappings (for the all-vs-all mode)\n");
		fprintf(fp_help, "    -p FLOAT     min secondary-to-primary score ratio [%g]\n", opt.pri_ratio);
		fprintf(fp_help, "    -N INT       retain at most INT secondary alignments [%d]\n", opt.best_n);
		fprintf(fp_help, "  Alignment:\n");
		fprintf(fp_help, "    -A INT       matching score [%d]\n", opt.a);
		fprintf(fp_help, "    -B INT       mismatch penalty [%d]\n", opt.b);
		fprintf(fp_help, "    -O INT[,INT] gap open penalty [%d,%d]\n", opt.q, opt.q2);
		fprintf(fp_help, "    -E INT[,INT] gap extension penalty; a k-long gap costs min{O1+k*E1,O2+k*E2} [%d,%d]\n", opt.e, opt.e2);
		fprintf(fp_help, "    -z INT[,INT] Z-drop score and inversion Z-drop score [%d,%d]\n", opt.zdrop, opt.zdrop_inv);
		fprintf(fp_help, "    -s INT       minimal peak DP alignment score [%d]\n", opt.min_dp_max);
		fprintf(fp_help, "    -u CHAR      how to find GT-AG. f:transcript strand, b:both strands, n:don't match GT-AG [n]\n");
		fprintf(fp_help, "  Input/Output:\n");
		fprintf(fp_help, "    -a           output in the SAM format (PAF by default)\n");
		fprintf(fp_help, "    -Q           don't output base quality in SAM\n");
		fprintf(fp_help, "    -L           write CIGAR with >65535 ops at the CG tag\n");
		fprintf(fp_help, "    -R STR       SAM read group line in a format like '@RG\\tID:foo\\tSM:bar' []\n");
		fprintf(fp_help, "    -c           output CIGAR in PAF\n");
		fprintf(fp_help, "    --cs[=STR]   output the cs tag; STR is 'short' (if absent) or 'long' [none]\n");
		fprintf(fp_help, "    --MD         output the MD tag\n");
		fprintf(fp_help, "    --eqx        write =/X CIGAR operators\n");
		fprintf(fp_help, "    -Y           use soft clipping for supplementary alignments\n");
		fprintf(fp_help, "    -t INT       number of threads [%d]\n", n_threads);
		fprintf(fp_help, "    -K NUM       minibatch size for mapping [500M]\n");
//		fprintf(fp_help, "    -v INT       verbose level [%d]\n", mm_verbose);
		fprintf(fp_help, "    --version    show version number\n");
		fprintf(fp_help, "  Preset:\n");
		fprintf(fp_help, "    -x STR       preset (always applied before other options; see minimap2.1 for details) []\n");
		fprintf(fp_help, "                 - map-pb/map-ont: PacBio/Nanopore vs reference mapping\n");
		fprintf(fp_help, "                 - ava-pb/ava-ont: PacBio/Nanopore read overlap\n");
		fprintf(fp_help, "                 - asm5/asm10/asm20: asm-to-ref mapping, for ~0.1/1/5%% sequence divergence\n");
		fprintf(fp_help, "                 - splice: long-read spliced alignment\n");
		fprintf(fp_help, "                 - sr: genomic short-read mapping\n");
		fprintf(fp_help, "\nSee `man ./minimap2.1' for detailed description of command-line options.\n");
		return fp_help == stdout? 0 : 1;
	}

	if ((opt.flag & MM_F_SR) && argc - optind > 3) {
		fprintf(stderr, "[ERROR] incorrect input: in the sr mode, please specify no more than two query files.\n");
		return 1;
	}
	idx_rdr = mm_idx_reader_open(argv[optind], &ipt, fnw);
	if (idx_rdr == 0) {
		fprintf(stderr, "[ERROR] failed to open file '%s'\n", argv[optind]);
		return 1;
	}
	if (!idx_rdr->is_idx && fnw == 0 && argc - optind < 2) {
		fprintf(stderr, "[ERROR] missing input: please specify a query file to map or option -d to keep the index\n");
		mm_idx_reader_close(idx_rdr);
		return 1;
	}
	if (opt.best_n == 0 && (opt.flag&MM_F_CIGAR) && mm_verbose >= 2)
		fprintf(stderr, "[WARNING]\033[1;31m `-N 0' reduces alignment accuracy. Please use --secondary=no to suppress secondary alignments.\033[0m\n");
	while ((mi = mm_idx_reader_read(idx_rdr, n_threads)) != 0) {
		if ((opt.flag & MM_F_CIGAR) && (mi->flag & MM_I_NO_SEQ)) {
			fprintf(stderr, "[ERROR] the prebuilt index doesn't contain sequences.\n");
			mm_idx_destroy(mi);
			mm_idx_reader_close(idx_rdr);
			return 1;
		}
		if ((opt.flag & MM_F_OUT_SAM) && idx_rdr->n_parts == 1) {
			if (mm_idx_reader_eof(idx_rdr)) {
				mm_write_sam_hdr(mi, rg, MM_VERSION, argc, argv);
			} else {
				mm_write_sam_hdr(0, rg, MM_VERSION, argc, argv);
				if (mm_verbose >= 2)
					fprintf(stderr, "[WARNING]\033[1;31m For a multi-part index, no @SQ lines will be outputted.\033[0m\n");
			}
		}
		if (mm_verbose >= 3)
			fprintf(stderr, "[M::%s::%.3f*%.2f] loaded/built the index for %d target sequence(s)\n",
					__func__, realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0), mi->n_seq);
		if (argc != optind + 1) mm_mapopt_update(&opt, mi);
		if (mm_verbose >= 3) mm_idx_stat(mi);
		if (!(opt.flag & MM_F_FRAG_MODE)) {
			for (i = optind + 1; i < argc; ++i)
				mm_map_file(mi, argv[i], &opt, n_threads);
		} else {
			mm_map_file_frag(mi, argc - (optind + 1), (const char**)&argv[optind + 1], &opt, n_threads);
		}
		mm_idx_destroy(mi);
	}
	mm_idx_reader_close(idx_rdr);

	if (fflush(stdout) == EOF) {
		fprintf(stderr, "[ERROR] failed to write the results\n");
		exit(EXIT_FAILURE);
	}

	if (mm_verbose >= 3) {
		fprintf(stderr, "[M::%s] Version: %s\n", __func__, MM_VERSION);
		fprintf(stderr, "[M::%s] CMD:", __func__);
		for (i = 0; i < argc; ++i)
			fprintf(stderr, " %s", argv[i]);
		fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec\n", __func__, realtime() - mm_realtime0, cputime());
	}
	return 0;
}
