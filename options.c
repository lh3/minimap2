#include <stdio.h>
#include "mmpriv.h"

void mm_idxopt_init(mm_idxopt_t *opt)
{
	memset(opt, 0, sizeof(mm_idxopt_t));
	opt->k = 15, opt->w = 10, opt->flag = 0;
	opt->bucket_bits = 14;
	opt->mini_batch_size = 50000000;
	opt->batch_size = 4000000000ULL;
}

void mm_mapopt_init(mm_mapopt_t *opt)
{
	memset(opt, 0, sizeof(mm_mapopt_t));
	opt->seed = 11;
	opt->mid_occ_frac = 2e-4f;
	opt->sdust_thres = 0; // no SDUST masking

	opt->min_cnt = 3;
	opt->min_chain_score = 40;
	opt->bw = 500;
	opt->max_gap = 5000;
	opt->max_gap_ref = -1;
	opt->max_chain_skip = 25;

	opt->mask_level = 0.5f;
	opt->pri_ratio = 0.8f;
	opt->best_n = 5;

	opt->max_join_long = 20000;
	opt->max_join_short = 2000;
	opt->min_join_flank_sc = 1000;
	opt->min_join_flank_ratio = 0.5f;

	opt->a = 2, opt->b = 4, opt->q = 4, opt->e = 2, opt->q2 = 24, opt->e2 = 1;
	opt->sc_ambi = 1;
	opt->zdrop = 400, opt->zdrop_inv = 200;
	opt->end_bonus = -1;
	opt->min_dp_max = opt->min_chain_score * opt->a;
	opt->min_ksw_len = 200;
	opt->anchor_ext_len = 20, opt->anchor_ext_shift = 6;
	opt->max_clip_ratio = 1.0f;
	opt->mini_batch_size = 500000000;

	opt->pe_ori = 0; // FF
	opt->pe_bonus = 33;
}

void mm_mapopt_update(mm_mapopt_t *opt, const mm_idx_t *mi)
{
	if ((opt->flag & MM_F_SPLICE_FOR) && (opt->flag & MM_F_SPLICE_REV))
		opt->flag |= MM_F_SPLICE;
	if (opt->mid_occ <= 0)
		opt->mid_occ = mm_idx_cal_max_occ(mi, opt->mid_occ_frac);
	if (opt->mid_occ < opt->min_mid_occ)
		opt->mid_occ = opt->min_mid_occ;
	if (mm_verbose >= 3)
		fprintf(stderr, "[M::%s::%.3f*%.2f] mid_occ = %d\n", __func__, realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0), opt->mid_occ);
}

void mm_mapopt_max_intron_len(mm_mapopt_t *opt, int max_intron_len)
{
	if ((opt->flag & MM_F_SPLICE) && max_intron_len > 0)
		opt->max_gap_ref = opt->bw = max_intron_len;
}

int mm_set_opt(const char *preset, mm_idxopt_t *io, mm_mapopt_t *mo)
{
	if (preset == 0) {
		mm_idxopt_init(io);
		mm_mapopt_init(mo);
	} else if (strcmp(preset, "ava-ont") == 0) {
		io->flag = 0, io->k = 15, io->w = 5;
		mo->flag |= MM_F_ALL_CHAINS | MM_F_NO_DIAG | MM_F_NO_DUAL | MM_F_NO_LJOIN;
		mo->min_chain_score = 100, mo->pri_ratio = 0.0f, mo->max_gap = 10000, mo->max_chain_skip = 25;
		mo->bw = 2000;
	} else if (strcmp(preset, "ava-pb") == 0) {
		io->flag |= MM_I_HPC, io->k = 19, io->w = 5;
		mo->flag |= MM_F_ALL_CHAINS | MM_F_NO_DIAG | MM_F_NO_DUAL | MM_F_NO_LJOIN;
		mo->min_chain_score = 100, mo->pri_ratio = 0.0f, mo->max_gap = 10000, mo->max_chain_skip = 25;
	} else if (strcmp(preset, "map10k") == 0 || strcmp(preset, "map-pb") == 0) {
		io->flag |= MM_I_HPC, io->k = 19;
	} else if (strcmp(preset, "map-ont") == 0) {
		io->flag = 0, io->k = 15;
	} else if (strcmp(preset, "asm5") == 0) {
		io->flag = 0, io->k = 19, io->w = 19;
		mo->a = 1, mo->b = 19, mo->q = 39, mo->q2 = 81, mo->e = 3, mo->e2 = 1, mo->zdrop = mo->zdrop_inv = 200;
		mo->min_mid_occ = 100;
		mo->min_dp_max = 200;
		mo->best_n = 50;
	} else if (strcmp(preset, "asm10") == 0) {
		io->flag = 0, io->k = 19, io->w = 19;
		mo->a = 1, mo->b = 9, mo->q = 16, mo->q2 = 41, mo->e = 2, mo->e2 = 1, mo->zdrop = mo->zdrop_inv = 200;
		mo->min_mid_occ = 100;
		mo->min_dp_max = 200;
		mo->best_n = 50;
	} else if (strcmp(preset, "asm20") == 0) {
		io->flag = 0, io->k = 19, io->w = 10;
		mo->a = 1, mo->b = 4, mo->q = 6, mo->q2 = 26, mo->e = 2, mo->e2 = 1, mo->zdrop = mo->zdrop_inv = 200;
		mo->min_mid_occ = 100;
		mo->min_dp_max = 200;
		mo->best_n = 50;
	} else if (strcmp(preset, "short") == 0 || strcmp(preset, "sr") == 0) {
		io->flag = 0, io->k = 21, io->w = 11;
		mo->flag |= MM_F_SR | MM_F_FRAG_MODE | MM_F_NO_PRINT_2ND | MM_F_2_IO_THREADS | MM_F_HEAP_SORT;
		mo->pe_ori = 0<<1|1; // FR
		mo->a = 2, mo->b = 8, mo->q = 12, mo->e = 2, mo->q2 = 24, mo->e2 = 1;
		mo->zdrop = mo->zdrop_inv = 100;
		mo->end_bonus = 10;
		mo->max_frag_len = 800;
		mo->max_gap = 100;
		mo->bw = 100;
		mo->pri_ratio = 0.5f;
		mo->min_cnt = 2;
		mo->min_chain_score = 25;
		mo->min_dp_max = 40;
		mo->best_n = 20;
		mo->mid_occ = 1000;
		mo->max_occ = 5000;
		mo->mini_batch_size = 50000000;
	} else if (strcmp(preset, "splice") == 0 || strcmp(preset, "cdna") == 0) {
		io->flag = 0, io->k = 15, io->w = 5;
		mo->flag |= MM_F_SPLICE | MM_F_SPLICE_FOR | MM_F_SPLICE_REV | MM_F_SPLICE_FLANK;
		mo->max_gap = 2000, mo->max_gap_ref = mo->bw = 200000;
		mo->a = 1, mo->b = 2, mo->q = 2, mo->e = 1, mo->q2 = 32, mo->e2 = 0;
		mo->noncan = 9;
		mo->zdrop = 200, mo->zdrop_inv = 100; // because mo->a is halved
	} else return -1;
	return 0;
}

int mm_check_opt(const mm_idxopt_t *io, const mm_mapopt_t *mo)
{
	if (io->k <= 0 || io->w <= 0) {
		if (mm_verbose >= 1)
			fprintf(stderr, "[ERROR]\033[1;31m -k and -w must be positive\033[0m\n");
		return -5;
	}
	if (mo->best_n < 0) {
		if (mm_verbose >= 1)
			fprintf(stderr, "[ERROR]\033[1;31m -N must be no less than 0\033[0m\n");
		return -4;
	}
	if (mo->best_n == 0 && mm_verbose >= 2)
		fprintf(stderr, "[WARNING]\033[1;31m '-N 0' reduces mapping accuracy. Please use '--secondary=no' instead.\033[0m\n");
	if (mo->pri_ratio < 0.0f || mo->pri_ratio > 1.0f) {
		if (mm_verbose >= 1)
			fprintf(stderr, "[ERROR]\033[1;31m -p must be within 0 and 1 (including 0 and 1)\033[0m\n");
		return -4;
	}
	if ((mo->flag & MM_F_FOR_ONLY) && (mo->flag & MM_F_REV_ONLY)) {
		if (mm_verbose >= 1)
			fprintf(stderr, "[ERROR]\033[1;31m --for-only and --rev-only can't be applied at the same time\033[0m\n");
		return -3;
	}
	if ((mo->q != mo->q2 || mo->e != mo->e2) && !(mo->e > mo->e2 && mo->q + mo->e < mo->q2 + mo->e2)) {
		if (mm_verbose >= 1)
			fprintf(stderr, "[ERROR]\033[1;31m dual gap penalties violating E1>E2 and O1+E1<O2+E2\033[0m\n");
		return -2;
	}
	if ((mo->q + mo->e) + (mo->q2 + mo->e2) > 127) {
		if (mm_verbose >= 1)
			fprintf(stderr, "[ERROR]\033[1;31m scoring system violating ({-O}+{-E})+({-O2}+{-E2}) <= 127\033[0m\n");
		return -1;
	}
	if (mo->zdrop < mo->zdrop_inv) {
		if (mm_verbose >= 1)
			fprintf(stderr, "[ERROR]\033[1;31m Z-drop should not be less than inversion-Z-drop\033[0m\n");
		return -5;
	}
	return 0;
}
