#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include "bseq.h"
#include "minimap.h"
#include "mmpriv.h"
#include "ketopt.h"

#define MM_VERSION "2.24-r1155-dirty"

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

static ko_longopt_t long_options[] = {
	{ "bucket-bits",    ko_required_argument, 300 },
	{ "mb-size",        ko_required_argument, 'K' },
	{ "seed",           ko_required_argument, 302 },
	{ "no-kalloc",      ko_no_argument,       303 },
	{ "print-qname",    ko_no_argument,       304 },
	{ "no-self",        ko_no_argument,       'D' },
	{ "print-seeds",    ko_no_argument,       306 },
	{ "max-chain-skip", ko_required_argument, 307 },
	{ "min-dp-len",     ko_required_argument, 308 },
	{ "print-aln-seq",  ko_no_argument,       309 },
	{ "splice",         ko_no_argument,       310 },
	{ "cost-non-gt-ag", ko_required_argument, 'C' },
	{ "no-long-join",   ko_no_argument,       312 },
	{ "sr",             ko_no_argument,       313 },
	{ "frag",           ko_required_argument, 314 },
	{ "secondary",      ko_required_argument, 315 },
	{ "cs",             ko_optional_argument, 316 },
	{ "end-bonus",      ko_required_argument, 317 },
	{ "no-pairing",     ko_no_argument,       318 },
	{ "splice-flank",   ko_required_argument, 319 },
	{ "idx-no-seq",     ko_no_argument,       320 },
	{ "end-seed-pen",   ko_required_argument, 321 },
	{ "for-only",       ko_no_argument,       322 },
	{ "rev-only",       ko_no_argument,       323 },
	{ "heap-sort",      ko_required_argument, 324 },
	{ "all-chain",      ko_no_argument,       'P' },
	{ "dual",           ko_required_argument, 326 },
	{ "max-clip-ratio", ko_required_argument, 327 },
	{ "min-occ-floor",  ko_required_argument, 328 },
	{ "MD",             ko_no_argument,       329 },
	{ "lj-min-ratio",   ko_required_argument, 330 },
	{ "score-N",        ko_required_argument, 331 },
	{ "eqx",            ko_no_argument,       332 },
	{ "paf-no-hit",     ko_no_argument,       333 },
	{ "split-prefix",   ko_required_argument, 334 },
	{ "no-end-flt",     ko_no_argument,       335 },
	{ "hard-mask-level",ko_no_argument,       336 },
	{ "cap-sw-mem",     ko_required_argument, 337 },
	{ "max-qlen",       ko_required_argument, 338 },
	{ "max-chain-iter", ko_required_argument, 339 },
	{ "junc-bed",       ko_required_argument, 340 },
	{ "junc-bonus",     ko_required_argument, 341 },
	{ "sam-hit-only",   ko_no_argument,       342 },
	{ "chain-gap-scale",ko_required_argument, 343 },
	{ "alt",            ko_required_argument, 344 },
	{ "alt-drop",       ko_required_argument, 345 },
	{ "mask-len",       ko_required_argument, 346 },
	{ "rmq",            ko_optional_argument, 347 },
	{ "qstrand",        ko_no_argument,       348 },
	{ "cap-kalloc",     ko_required_argument, 349 },
	{ "q-occ-frac",     ko_required_argument, 350 },
	{ "chain-skip-scale",ko_required_argument,351 },
	{ "print-chains",   ko_no_argument,       352 },
	{ "no-hash-name",   ko_no_argument,       353 },
	{ "secondary-seq",  ko_no_argument,       354 },
	{ "help",           ko_no_argument,       'h' },
	{ "max-intron-len", ko_required_argument, 'G' },
	{ "version",        ko_no_argument,       'V' },
	{ "min-count",      ko_required_argument, 'n' },
	{ "min-chain-score",ko_required_argument, 'm' },
	{ "mask-level",     ko_required_argument, 'M' },
	{ "min-dp-score",   ko_required_argument, 's' },
	{ "sam",            ko_no_argument,       'a' },
	{ 0, 0, 0 }
};

static inline int64_t mm_parse_num2(const char *str, char **q)
{
	double x;
	char *p;
	x = strtod(str, &p);
	if (*p == 'G' || *p == 'g') x *= 1e9, ++p;
	else if (*p == 'M' || *p == 'm') x *= 1e6, ++p;
	else if (*p == 'K' || *p == 'k') x *= 1e3, ++p;
	if (q) *q = p;
	return (int64_t)(x + .499);
}

static inline int64_t mm_parse_num(const char *str)
{
	return mm_parse_num2(str, 0);
}

static inline void yes_or_no(mm_mapopt_t *opt, int64_t flag, int long_idx, const char *arg, int yes_to_set)
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
	const char *opt_str = "2aSDw:k:K:t:r:f:Vv:g:G:I:d:XT:s:x:Hcp:M:n:z:A:B:O:E:m:N:Qu:R:hF:LC:yYPo:e:U:j:";
	ketopt_t o = KETOPT_INIT;
	mm_mapopt_t opt;
	mm_idxopt_t ipt;
	int i, c, n_threads = 3, n_parts, old_best_n = -1;
	char *fnw = 0, *rg = 0, *junc_bed = 0, *s, *alt_list = 0;
	FILE *fp_help = stderr;
	mm_idx_reader_t *idx_rdr;
	mm_idx_t *mi;

	mm_verbose = 3;
	liftrlimit();
	mm_realtime0 = realtime();
	mm_set_opt(0, &ipt, &opt);

	while ((c = ketopt(&o, argc, argv, 1, opt_str, long_options)) >= 0) { // test command line options and apply option -x/preset first
		if (c == 'x') {
			if (mm_set_opt(o.arg, &ipt, &opt) < 0) {
				fprintf(stderr, "[ERROR] unknown preset '%s'\n", o.arg);
				return 1;
			}
		} else if (c == ':') {
			fprintf(stderr, "[ERROR] missing option argument\n");
			return 1;
		} else if (c == '?') {
			fprintf(stderr, "[ERROR] unknown option in \"%s\"\n", argv[o.i - 1]);
			return 1;
		}
	}
	o = KETOPT_INIT;

	while ((c = ketopt(&o, argc, argv, 1, opt_str, long_options)) >= 0) {
		if (c == 'w') ipt.w = atoi(o.arg), ipt.flag &= ~MM_I_SYNCMER;
		else if (c == 'j') ipt.w = atoi(o.arg), ipt.flag |= MM_I_SYNCMER;
		else if (c == 'k') ipt.k = atoi(o.arg);
		else if (c == 'H') ipt.flag |= MM_I_HPC;
		else if (c == 'd') fnw = o.arg; // the above are indexing related options, except -I
		else if (c == 't') n_threads = atoi(o.arg);
		else if (c == 'v') mm_verbose = atoi(o.arg);
		else if (c == 'g') opt.max_gap = (int)mm_parse_num(o.arg);
		else if (c == 'G') mm_mapopt_max_intron_len(&opt, (int)mm_parse_num(o.arg));
		else if (c == 'F') opt.max_frag_len = (int)mm_parse_num(o.arg);
		else if (c == 'N') old_best_n = opt.best_n, opt.best_n = atoi(o.arg);
		else if (c == 'p') opt.pri_ratio = atof(o.arg);
		else if (c == 'M') opt.mask_level = atof(o.arg);
		else if (c == 'c') opt.flag |= MM_F_OUT_CG | MM_F_CIGAR;
		else if (c == 'D') opt.flag |= MM_F_NO_DIAG;
		else if (c == 'P') opt.flag |= MM_F_ALL_CHAINS;
		else if (c == 'X') opt.flag |= MM_F_ALL_CHAINS | MM_F_NO_DIAG | MM_F_NO_DUAL | MM_F_NO_LJOIN; // -D -P --no-long-join --dual=no
		else if (c == 'a') opt.flag |= MM_F_OUT_SAM | MM_F_CIGAR;
		else if (c == 'Q') opt.flag |= MM_F_NO_QUAL;
		else if (c == 'Y') opt.flag |= MM_F_SOFTCLIP;
		else if (c == 'L') opt.flag |= MM_F_LONG_CIGAR;
		else if (c == 'y') opt.flag |= MM_F_COPY_COMMENT;
		else if (c == 'T') opt.sdust_thres = atoi(o.arg);
		else if (c == 'n') opt.min_cnt = atoi(o.arg);
		else if (c == 'm') opt.min_chain_score = atoi(o.arg);
		else if (c == 'A') opt.a = atoi(o.arg);
		else if (c == 'B') opt.b = atoi(o.arg);
		else if (c == 's') opt.min_dp_max = atoi(o.arg);
		else if (c == 'C') opt.noncan = atoi(o.arg);
		else if (c == 'I') ipt.batch_size = mm_parse_num(o.arg);
		else if (c == 'K') opt.mini_batch_size = mm_parse_num(o.arg);
		else if (c == 'e') opt.occ_dist = mm_parse_num(o.arg);
		else if (c == 'R') rg = o.arg;
		else if (c == 'h') fp_help = stdout;
		else if (c == '2') opt.flag |= MM_F_2_IO_THREADS;
		else if (c == 'o') {
			if (strcmp(o.arg, "-") != 0) {
				if (freopen(o.arg, "wb", stdout) == NULL) {
					fprintf(stderr, "[ERROR]\033[1;31m failed to write the output to file '%s'\033[0m: %s\n", o.arg, strerror(errno));
					exit(1);
				}
			}
		}
		else if (c == 300) ipt.bucket_bits = atoi(o.arg); // --bucket-bits
		else if (c == 302) opt.seed = atoi(o.arg); // --seed
		else if (c == 303) mm_dbg_flag |= MM_DBG_NO_KALLOC; // --no-kalloc
		else if (c == 304) mm_dbg_flag |= MM_DBG_PRINT_QNAME; // --print-qname
		else if (c == 306) mm_dbg_flag |= MM_DBG_PRINT_QNAME | MM_DBG_PRINT_SEED, n_threads = 1; // --print-seed
		else if (c == 307) opt.max_chain_skip = atoi(o.arg); // --max-chain-skip
		else if (c == 339) opt.max_chain_iter = atoi(o.arg); // --max-chain-iter
		else if (c == 308) opt.min_ksw_len = atoi(o.arg); // --min-dp-len
		else if (c == 309) mm_dbg_flag |= MM_DBG_PRINT_QNAME | MM_DBG_PRINT_ALN_SEQ, n_threads = 1; // --print-aln-seq
		else if (c == 310) opt.flag |= MM_F_SPLICE; // --splice
		else if (c == 312) opt.flag |= MM_F_NO_LJOIN; // --no-long-join
		else if (c == 313) opt.flag |= MM_F_SR; // --sr
		else if (c == 317) opt.end_bonus = atoi(o.arg); // --end-bonus
		else if (c == 318) opt.flag |= MM_F_INDEPEND_SEG; // --no-pairing
		else if (c == 320) ipt.flag |= MM_I_NO_SEQ; // --idx-no-seq
		else if (c == 321) opt.anchor_ext_shift = atoi(o.arg); // --end-seed-pen
		else if (c == 322) opt.flag |= MM_F_FOR_ONLY; // --for-only
		else if (c == 323) opt.flag |= MM_F_REV_ONLY; // --rev-only
		else if (c == 327) opt.max_clip_ratio = atof(o.arg); // --max-clip-ratio
		else if (c == 328) opt.min_mid_occ = atoi(o.arg); // --min-occ-floor
		else if (c == 329) opt.flag |= MM_F_OUT_MD; // --MD
		else if (c == 331) opt.sc_ambi = atoi(o.arg); // --score-N
		else if (c == 332) opt.flag |= MM_F_EQX; // --eqx
		else if (c == 333) opt.flag |= MM_F_PAF_NO_HIT; // --paf-no-hit
		else if (c == 334) opt.split_prefix = o.arg; // --split-prefix
		else if (c == 335) opt.flag |= MM_F_NO_END_FLT; // --no-end-flt
		else if (c == 336) opt.flag |= MM_F_HARD_MLEVEL; // --hard-mask-level
		else if (c == 337) opt.max_sw_mat = mm_parse_num(o.arg); // --cap-sw-mat
		else if (c == 338) opt.max_qlen = mm_parse_num(o.arg); // --max-qlen
		else if (c == 340) junc_bed = o.arg; // --junc-bed
		else if (c == 341) opt.junc_bonus = atoi(o.arg); // --junc-bonus
		else if (c == 342) opt.flag |= MM_F_SAM_HIT_ONLY; // --sam-hit-only
		else if (c == 343) opt.chain_gap_scale = atof(o.arg); // --chain-gap-scale
		else if (c == 351) opt.chain_skip_scale = atof(o.arg); // --chain-skip-scale
		else if (c == 344) alt_list = o.arg; // --alt
		else if (c == 345) opt.alt_drop = atof(o.arg); // --alt-drop
		else if (c == 346) opt.mask_len = mm_parse_num(o.arg); // --mask-len
		else if (c == 348) opt.flag |= MM_F_QSTRAND | MM_F_NO_INV; // --qstrand
		else if (c == 349) opt.cap_kalloc = mm_parse_num(o.arg); // --cap-kalloc
		else if (c == 350) opt.q_occ_frac = atof(o.arg); // --q-occ-frac
		else if (c == 352) mm_dbg_flag |= MM_DBG_PRINT_CHAIN; // --print-chains
		else if (c == 353) opt.flag |= MM_F_NO_HASH_NAME; // --no-hash-name
		else if (c == 347) opt.flag |= MM_F_SECONDARY_SEQ; // --secondary-seq
		else if (c == 330) {
			fprintf(stderr, "[WARNING] \033[1;31m --lj-min-ratio has been deprecated.\033[0m\n");
		} else if (c == 314) { // --frag
			yes_or_no(&opt, MM_F_FRAG_MODE, o.longidx, o.arg, 1);
		} else if (c == 315) { // --secondary
			yes_or_no(&opt, MM_F_NO_PRINT_2ND, o.longidx, o.arg, 0);
		} else if (c == 316) { // --cs
			opt.flag |= MM_F_OUT_CS | MM_F_CIGAR;
			if (o.arg == 0 || strcmp(o.arg, "short") == 0) {
				opt.flag &= ~MM_F_OUT_CS_LONG;
			} else if (strcmp(o.arg, "long") == 0) {
				opt.flag |= MM_F_OUT_CS_LONG;
			} else if (strcmp(o.arg, "none") == 0) {
				opt.flag &= ~MM_F_OUT_CS;
			} else if (mm_verbose >= 2) {
				fprintf(stderr, "[WARNING]\033[1;31m --cs only takes 'short' or 'long'. Invalid values are assumed to be 'short'.\033[0m\n");
			}
		} else if (c == 319) { // --splice-flank
			yes_or_no(&opt, MM_F_SPLICE_FLANK, o.longidx, o.arg, 1);
		} else if (c == 324) { // --heap-sort
			yes_or_no(&opt, MM_F_HEAP_SORT, o.longidx, o.arg, 1);
		} else if (c == 326) { // --dual
			yes_or_no(&opt, MM_F_NO_DUAL, o.longidx, o.arg, 0);
		} else if (c == 347) { // --rmq
			if (o.arg) yes_or_no(&opt, MM_F_RMQ, o.longidx, o.arg, 1);
			else opt.flag |= MM_F_RMQ;
		} else if (c == 'S') {
			opt.flag |= MM_F_OUT_CS | MM_F_CIGAR | MM_F_OUT_CS_LONG;
			if (mm_verbose >= 2)
				fprintf(stderr, "[WARNING]\033[1;31m option -S is deprecated and may be removed in future. Please use --cs=long instead.\033[0m\n");
		} else if (c == 'V') {
			puts(MM_VERSION);
			return 0;
		} else if (c == 'r') {
			opt.bw = (int)mm_parse_num2(o.arg, &s);
			if (*s == ',') opt.bw_long = (int)mm_parse_num2(s + 1, &s);
		} else if (c == 'U') {
			opt.min_mid_occ = strtol(o.arg, &s, 10);
			if (*s == ',') opt.max_mid_occ = strtol(s + 1, &s, 10);
		} else if (c == 'f') {
			double x;
			char *p;
			x = strtod(o.arg, &p);
			if (x < 1.0) opt.mid_occ_frac = x, opt.mid_occ = 0;
			else opt.mid_occ = (int)(x + .499);
			if (*p == ',') opt.max_occ = (int)(strtod(p+1, &p) + .499);
		} else if (c == 'u') {
			if (*o.arg == 'b') opt.flag |= MM_F_SPLICE_FOR|MM_F_SPLICE_REV; // both strands
			else if (*o.arg == 'f') opt.flag |= MM_F_SPLICE_FOR, opt.flag &= ~MM_F_SPLICE_REV; // match GT-AG
			else if (*o.arg == 'r') opt.flag |= MM_F_SPLICE_REV, opt.flag &= ~MM_F_SPLICE_FOR; // match CT-AC (reverse complement of GT-AG)
			else if (*o.arg == 'n') opt.flag &= ~(MM_F_SPLICE_FOR|MM_F_SPLICE_REV); // don't try to match the GT-AG signal
			else {
				fprintf(stderr, "[ERROR]\033[1;31m unrecognized cDNA direction\033[0m\n");
				return 1;
			}
		} else if (c == 'z') {
			opt.zdrop = opt.zdrop_inv = strtol(o.arg, &s, 10);
			if (*s == ',') opt.zdrop_inv = strtol(s + 1, &s, 10);
		} else if (c == 'O') {
			opt.q = opt.q2 = strtol(o.arg, &s, 10);
			if (*s == ',') opt.q2 = strtol(s + 1, &s, 10);
		} else if (c == 'E') {
			opt.e = opt.e2 = strtol(o.arg, &s, 10);
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
	if (opt.best_n == 0) {
		fprintf(stderr, "[WARNING]\033[1;31m changed '-N 0' to '-N %d --secondary=no'.\033[0m\n", old_best_n);
		opt.best_n = old_best_n, opt.flag |= MM_F_NO_PRINT_2ND;
	}

	if (argc == o.ind || fp_help == stdout) {
		fprintf(fp_help, "Usage: minimap2 [options] <target.fa>|<target.idx> [query.fa] [...]\n");
		fprintf(fp_help, "Options:\n");
		fprintf(fp_help, "  Indexing:\n");
		fprintf(fp_help, "    -H           use homopolymer-compressed k-mer (preferrable for PacBio)\n");
		fprintf(fp_help, "    -k INT       k-mer size (no larger than 28) [%d]\n", ipt.k);
		fprintf(fp_help, "    -w INT       minimizer window size [%d]\n", ipt.w);
		fprintf(fp_help, "    -j INT       syncmer submer size (overriding -w) []\n");
		fprintf(fp_help, "    -I NUM       split index for every ~NUM input bases [4G]\n");
		fprintf(fp_help, "    -d FILE      dump index to FILE []\n");
		fprintf(fp_help, "  Mapping:\n");
		fprintf(fp_help, "    -f FLOAT     filter out top FLOAT fraction of repetitive minimizers [%g]\n", opt.mid_occ_frac);
		fprintf(fp_help, "    -g NUM       stop chain enlongation if there are no minimizers in INT-bp [%d]\n", opt.max_gap);
		fprintf(fp_help, "    -G NUM       max intron length (effective with -xsplice; changing -r) [200k]\n");
		fprintf(fp_help, "    -F NUM       max fragment length (effective with -xsr or in the fragment mode) [800]\n");
		fprintf(fp_help, "    -r NUM[,NUM] chaining/alignment bandwidth and long-join bandwidth [%d,%d]\n", opt.bw, opt.bw_long);
		fprintf(fp_help, "    -n INT       minimal number of minimizers on a chain [%d]\n", opt.min_cnt);
		fprintf(fp_help, "    -m INT       minimal chaining score (matching bases minus log gap penalty) [%d]\n", opt.min_chain_score);
//		fprintf(fp_help, "    -T INT       SDUST threshold; 0 to disable SDUST [%d]\n", opt.sdust_thres); // TODO: this option is never used; might be buggy
		fprintf(fp_help, "    -X           skip self and dual mappings (for the all-vs-all mode)\n");
		fprintf(fp_help, "    -p FLOAT     min secondary-to-primary score ratio [%g]\n", opt.pri_ratio);
		fprintf(fp_help, "    -N INT       retain at most INT secondary alignments [%d]\n", opt.best_n);
		fprintf(fp_help, "  Alignment:\n");
		fprintf(fp_help, "    -A INT       matching score [%d]\n", opt.a);
		fprintf(fp_help, "    -B INT       mismatch penalty (larger value for lower divergence) [%d]\n", opt.b);
		fprintf(fp_help, "    -O INT[,INT] gap open penalty [%d,%d]\n", opt.q, opt.q2);
		fprintf(fp_help, "    -E INT[,INT] gap extension penalty; a k-long gap costs min{O1+k*E1,O2+k*E2} [%d,%d]\n", opt.e, opt.e2);
		fprintf(fp_help, "    -z INT[,INT] Z-drop score and inversion Z-drop score [%d,%d]\n", opt.zdrop, opt.zdrop_inv);
		fprintf(fp_help, "    -s INT       minimal peak DP alignment score [%d]\n", opt.min_dp_max);
		fprintf(fp_help, "    -u CHAR      how to find GT-AG. f:transcript strand, b:both strands, n:don't match GT-AG [n]\n");
		fprintf(fp_help, "  Input/Output:\n");
		fprintf(fp_help, "    -a           output in the SAM format (PAF by default)\n");
		fprintf(fp_help, "    -o FILE      output alignments to FILE [stdout]\n");
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
		fprintf(fp_help, "                 - map-pb/map-ont - PacBio CLR/Nanopore vs reference mapping\n");
		fprintf(fp_help, "                 - map-hifi - PacBio HiFi reads vs reference mapping\n");
		fprintf(fp_help, "                 - ava-pb/ava-ont - PacBio/Nanopore read overlap\n");
		fprintf(fp_help, "                 - asm5/asm10/asm20 - asm-to-ref mapping, for ~0.1/1/5%% sequence divergence\n");
		fprintf(fp_help, "                 - splice/splice:hq - long-read/Pacbio-CCS spliced alignment\n");
		fprintf(fp_help, "                 - sr - genomic short-read mapping\n");
		fprintf(fp_help, "\nSee `man ./minimap2.1' for detailed description of these and other advanced command-line options.\n");
		return fp_help == stdout? 0 : 1;
	}

	if ((opt.flag & MM_F_SR) && argc - o.ind > 3) {
		fprintf(stderr, "[ERROR] incorrect input: in the sr mode, please specify no more than two query files.\n");
		return 1;
	}
	idx_rdr = mm_idx_reader_open(argv[o.ind], &ipt, fnw);
	if (idx_rdr == 0) {
		fprintf(stderr, "[ERROR] failed to open file '%s': %s\n", argv[o.ind], strerror(errno));
		return 1;
	}
	if (!idx_rdr->is_idx && fnw == 0 && argc - o.ind < 2) {
		fprintf(stderr, "[ERROR] missing input: please specify a query file to map or option -d to keep the index\n");
		mm_idx_reader_close(idx_rdr);
		return 1;
	}
	if (opt.best_n == 0 && (opt.flag&MM_F_CIGAR) && mm_verbose >= 2)
		fprintf(stderr, "[WARNING]\033[1;31m `-N 0' reduces alignment accuracy. Please use --secondary=no to suppress secondary alignments.\033[0m\n");
	while ((mi = mm_idx_reader_read(idx_rdr, n_threads)) != 0) {
		int ret;
		if ((opt.flag & MM_F_CIGAR) && (mi->flag & MM_I_NO_SEQ)) {
			fprintf(stderr, "[ERROR] the prebuilt index doesn't contain sequences.\n");
			mm_idx_destroy(mi);
			mm_idx_reader_close(idx_rdr);
			return 1;
		}
		if ((opt.flag & MM_F_OUT_SAM) && idx_rdr->n_parts == 1) {
			if (mm_idx_reader_eof(idx_rdr)) {
				if (opt.split_prefix == 0)
					ret = mm_write_sam_hdr(mi, rg, MM_VERSION, argc, argv);
				else
					ret = mm_write_sam_hdr(0, rg, MM_VERSION, argc, argv);
			} else {
				ret = mm_write_sam_hdr(0, rg, MM_VERSION, argc, argv);
				if (opt.split_prefix == 0 && mm_verbose >= 2)
					fprintf(stderr, "[WARNING]\033[1;31m For a multi-part index, no @SQ lines will be outputted. Please use --split-prefix.\033[0m\n");
			}
			if (ret != 0) {
				mm_idx_destroy(mi);
				mm_idx_reader_close(idx_rdr);
				return 1;
			}
		}
		if (mm_verbose >= 3)
			fprintf(stderr, "[M::%s::%.3f*%.2f] loaded/built the index for %d target sequence(s)\n",
					__func__, realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0), mi->n_seq);
		if (argc != o.ind + 1) mm_mapopt_update(&opt, mi);
		if (mm_verbose >= 3) mm_idx_stat(mi);
		if (junc_bed) mm_idx_bed_read(mi, junc_bed, 1);
		if (alt_list) mm_idx_alt_read(mi, alt_list);
		if (argc - (o.ind + 1) == 0) {
			mm_idx_destroy(mi);
			continue; // no query files
		}
		ret = 0;
		if (!(opt.flag & MM_F_FRAG_MODE)) {
			for (i = o.ind + 1; i < argc; ++i) {
				ret = mm_map_file(mi, argv[i], &opt, n_threads);
				if (ret < 0) break;
			}
		} else {
			ret = mm_map_file_frag(mi, argc - (o.ind + 1), (const char**)&argv[o.ind + 1], &opt, n_threads);
		}
		mm_idx_destroy(mi);
		if (ret < 0) {
			fprintf(stderr, "ERROR: failed to map the query file\n");
			exit(EXIT_FAILURE);
		}
	}
	n_parts = idx_rdr->n_parts;
	mm_idx_reader_close(idx_rdr);

	if (opt.split_prefix)
		mm_split_merge(argc - (o.ind + 1), (const char**)&argv[o.ind + 1], &opt, n_parts);

	if (fflush(stdout) == EOF) {
		perror("[ERROR] failed to write the results");
		exit(EXIT_FAILURE);
	}

	if (mm_verbose >= 3) {
		fprintf(stderr, "[M::%s] Version: %s\n", __func__, MM_VERSION);
		fprintf(stderr, "[M::%s] CMD:", __func__);
		for (i = 0; i < argc; ++i)
			fprintf(stderr, " %s", argv[i]);
		fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec; Peak RSS: %.3f GB\n", __func__, realtime() - mm_realtime0, cputime(), peakrss() / 1024.0 / 1024.0 / 1024.0);
	}
	return 0;
}
