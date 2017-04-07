#include "minimap.h"

void mm_mapopt_init(mm_mapopt_t *opt)
{
	opt->radius = 500;
	opt->max_gap = 10000;
	opt->min_cnt = 4;
	opt->min_match = 40;
	opt->sdust_thres = 0;
	opt->flag = MM_F_WITH_REP;
	opt->merge_frac = .5;
}
