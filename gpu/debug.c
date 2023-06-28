#include "debug.h"

#include <time.h>

#include "plchain.h"

// FILE *f_anchors = NULL;
// FILE *range_input = NULL;
// FILE *binary = NULL;
// bool use_binary_input = true;
// FILE *range = NULL;

#ifdef DEBUG_VERBOSE

/////////////////////////////////////////////////////////////////////
///////////        Print Input Files            /////////////////////
/////////////////////////////////////////////////////////////////////
void debug_output_anchors(const char debug_folder[], chain_read_t *in) {
    static FILE *f_anchors = NULL;
    if (!f_anchors) {
        char anchors_filename[50];
        strcpy(anchors_filename, debug_folder);
        strcat(anchors_filename, ".anchor.out");
        if ((f_anchors = fopen(anchors_filename, "w+")) == NULL) {
            fprintf(stderr, "[Error] Cannot create output file %s\n",
                    anchors_filename);
            exit(1);
        }
        fprintf(stderr, "[Info] Writing anchors to file %s\n",
                anchors_filename);
        fprintf(f_anchors, "@@@<qname\tqlen\n");
    }

    /* Write Sequence Name and Length, rep_len*/
    fprintf(f_anchors, "<%s\t%d\n", in->seq.name, in->seq.len);
    fprintf(f_anchors, "*%d\n", in->seq.rep_len);

    /* Read Number of Anchors */
    fprintf(f_anchors, "#%d\n", in->n);

    /* Read Anchors */
    for (int i = 0; i < in->n; i++) {
        fprintf(f_anchors, "%lx,%lx\t", in->a[i].x, in->a[i].y);
    }
    fprintf(f_anchors, "\n");
}

void debug_output_score(const char debug_folder[], chain_read_t *in) {
    static FILE *f_score = NULL;
    if (!f_score) {
        char score_filename[50];
        strcpy(score_filename, debug_folder);
        strcat(score_filename, ".goldscore.out");
        if ((f_score = fopen(score_filename, "w+")) == NULL) {
            fprintf(stderr, "[Error] Cannot crea input file %s\n",
                    score_filename);
            exit(1);
        }
        fprintf(stderr, "[Info] Writing gold score to file %s\n",
                score_filename);
        fprintf(f_score, "@@@<qname\tqlen\n");
    }

    /* Write Sequence Name and Length */
    fprintf(f_score, "<%s %d\n", in->seq.name, in->seq.len);

    /* Write Score and Predecesser */
    fprintf(f_score, "#%ld\n", in->n);

    /* Write score */
    for (int i = 0; i < in->n; i++) {
        fprintf(f_score, "%d,%ld\t", in->f[i], in->p[i]);
    }
    fprintf(f_score, "\n");
}

void debug_output_meta(const char debug_folder[], input_meta_t *meta) {
    static FILE *f_metaout = NULL;
    if (!f_metaout) {
        char *buf = NULL;
        size_t len = 0;
        char meta_filename[50];
        strcpy(meta_filename, debug_folder);
        strcat(meta_filename, ".meta.out");
        if ((f_metaout = fopen(meta_filename, "w+")) == NULL) {
            fprintf(stderr, "[Error] Cannot create input file %s\n",
                    meta_filename);
            exit(1);
        }
        fprintf(stderr, "[Info] Writing meta data to file %s\n", meta_filename);
        fprintf(f_metaout, "@@@>tname\ttlen\n");
    }

    /* Write Number of reference Sequences */
    fprintf(f_metaout, "#%d\n", meta->n_refs);

    /* Write Reference Seq metadata */
    for (int i = 0; i < meta->n_refs; i++) {
        fprintf(f_metaout, ">%s\t%d\n", meta->refs[i].name, meta->refs[i].len);
    }
}

void debug_print_successor_range(int32_t *range, int64_t n) {
    static FILE *fout_range = NULL;
    static int read_idx = 0;
    if (fout_range == NULL) {
        char fout_range_filename[50];
        strcpy(fout_range_filename, debug_folder);
        strcat(fout_range_filename, ".range.out");
        if ((fout_range = fopen(fout_range_filename, "w+")) == NULL) {
            fprintf(stderr, "[Error]: Cannot create range output file: %s \n",
                   fout_range_filename);
            exit(1);
        }
        fprintf(stderr, "[Info] Writing successor range to file %s\n",
                fout_range_filename);
    }
    fprintf(fout_range, "> %ld, len: %ld ", read_idx, n);
    for (int64_t i = 0; i < n; ++i) {
        fprintf(fout_range, "#%ld: %d, ", i, range[i]);
    }
    fprintf(fout_range, "\n");
    read_idx++;
}

int debug_print_cut(const size_t *cut, size_t max_cut, size_t n,
                    size_t offset, char* qname) {
    static FILE *fout_cut = NULL;
    static int read_idx = 0;
    if (fout_cut == NULL) {
        char fout_cut_filename[50];
        strcpy(fout_cut_filename, debug_folder);
        strcat(fout_cut_filename, ".cut.out");
        if ((fout_cut = fopen(fout_cut_filename, "w+")) == NULL) {
            fprintf(stderr, "[Error]: Cannot create cut output file: %s \n",
                   fout_cut_filename);
            exit(1);
        }
        fprintf(stderr, "[Info] Writing cut to file %s\n", fout_cut_filename);
    }
    fprintf(fout_cut, "> %s, len: %ld offset %ld ", qname == NULL ? "--" : qname, n, offset);
    size_t cid = 0;
    for (; cid < max_cut && (cut[cid] < n + offset || cut[cid] == SIZE_MAX);
         cid++) {
        if (cut[cid] != SIZE_MAX)
            fprintf(fout_cut, "%zu(%zu)\t", cut[cid] - offset, cut[cid]);
        else
            fprintf(fout_cut, "x\t");
    }
    fprintf(fout_cut, "\n");
    read_idx++;
    return cid;
}

void debug_print_score(const int64_t *p, const int32_t *score, int64_t n) {
    static FILE *fout_score = NULL;
    static int read_idx = 0;
    if (fout_score == NULL) {
        char fout_score_filename[50];
        strcpy(fout_score_filename, debug_folder);
        strcat(fout_score_filename, ".score.out");
        if ((fout_score = fopen(fout_score_filename, "w+")) == NULL) {
            fprintf(stderr, "[Error]: Cannot create score output file: %s \n",
                   fout_score_filename);
            exit(1);
        }
        fprintf(stderr, "[Info] Writing score to file %s\n",
                fout_score_filename);
        fprintf(fout_score, "@@@<qname\tqlen\n");
    }
    fprintf(fout_score, "<%ld\t\n", read_idx);
    fprintf(fout_score, "#%ld\n", n);
    for (int i = 0; i < n; ++i) {
        fprintf(fout_score, "%d,%ld\t", score[i], p[i]);
    }
    fprintf(fout_score, "\n");
    read_idx++;
}

void debug_print_chain(mm128_t *a, uint64_t *u, int32_t n_u, char* qname) {
    static FILE *fout_chain = NULL;
    if (fout_chain == NULL) {
        char fout_chain_filename[50];
        strcpy(fout_chain_filename, debug_folder);
        strcat(fout_chain_filename, ".chain.out");
        if ((fout_chain = fopen(fout_chain_filename, "w+")) == NULL) {
            fprintf(stderr, "[Error]: Cannot create chain output file: %s \n",
                fout_chain_filename);
            exit(1);
        }
        fprintf(stderr, "[Info] Writing chain to file %s\n",
                fout_chain_filename);
        fprintf(fout_chain, "[score] anchors\n");
    }
    fprintf(fout_chain, "<%s\n", qname);
    for (int i = 0, j = 0; i < n_u; i++) {
        fprintf(fout_chain, "[%d] #%d: ", u[i] >> 32, (uint32_t)u[i]);
        for (int new_j = j + (uint32_t)u[i]; j < new_j; j++) {
            fprintf(fout_chain, "%lx,%lx ", a[j].x, a[j].y);
        }
        fprintf(fout_chain, "\n");
    }
}

void debug_print_regs(mm_reg1_t* regs, int n_u, char* qname){
    static FILE *fout_regs = NULL;
    if (fout_regs == NULL){
        char fout_regs_filename[50];
        strcpy(fout_regs_filename, debug_folder);
        strcat(fout_regs_filename, ".regs.out");
        if ((fout_regs = fopen(fout_regs_filename, "w+")) == NULL) {
            fprintf(stderr, "[Error]: Cannot create print output file: %s \n",
                    fout_regs_filename);
            exit(1);
        }
        fprintf(stderr, "[Info] Writing regs to file %s\n", fout_regs_filename);
        fprintf(fout_regs, "[regs] \n");
    }
    fprintf(fout_regs, "<%s\n", qname);
    for (int i = 0; i < n_u; i++){
        fprintf(fout_regs,
                "[%d] cnt %d rid %d score %d qs %d qe %d rs %d re %d parent %d "
                "subsc %d as %d mlen %d blen %d n_sub %d score0 %d\n", regs[i].id,
                regs[i].cnt, regs[i].rid, regs[i].score, regs[i].qs, regs[i].qe,
                regs[i].rs, regs[i].re, regs[i].parent, regs[i].subsc,
                regs[i].as, regs[i].mlen, regs[i].blen, regs[i].n_sub,
                regs[i].score0);
    }
}

#endif // DEBUG_VERBOSE

///////////////////////////////////////////////////////////////////////////
/////////////           check functions     ///////////////////////////////
///////////////////////////////////////////////////////////////////////////
#ifdef DEBUG_CHECK
/**
 * Read Plaintxt input file for Chaining scores from <debug_folder>.score
 * Allocate and Populate chain_read_t.f, chain_read_t.p
 * Return number of anchors if success, -1 if failed
 */
int debug_skip_score_check = 0;
int debug_no_score_file = 0;
int debug_read_score(const char input_filename[], chain_read_t *in, void *km) {
    static FILE *f_score = NULL;
    if (debug_no_score_file) return -1;
    if (!f_score) {
        char *buf = NULL;
        size_t len = 0;
        char score_filename[50];
        strcpy(score_filename, input_filename);
        strcat(score_filename, ".score");
        if ((f_score = fopen(score_filename, "r")) == NULL) {
            in->f = NULL;
            in->p = NULL;
            debug_skip_score_check = 1;
            debug_no_score_file = 1;
            fprintf(stderr, "[Warning] Cannot open score file %s, skip score checking! \n",
                    score_filename);
            return -1;
        }
        fprintf(stderr, "[Info] Reading gold score from file %s\n",
                score_filename);
        if (getline(&buf, &len, f_score) <= 0) {
            // discard header line.
            fprintf(stderr, "[Error] wrong format: %s %s\n", score_filename, buf);
            return -1;
        }
        if (strlen(buf) < 3 || strncmp(buf, "@@@", 3)) {
            fprintf(stderr, "[Error] wrong format %s %s\n", score_filename, buf);
            return -1;
        };
        kfree(km, buf);
    }

    /* Read Sequence Name and Length */
    char buf[100];
    int seqlen = -1;
    int num_anchor = -1;
    if (fscanf(f_score, "<%s %d\n", buf, &seqlen) != 2) {
        return -1;
    }

    if (strcmp(buf, in->seq.name) || seqlen != in->seq.len) {
        fprintf(stderr, "[Error] query sequence mismatch: %s %d\n", buf,
                seqlen);
        return -1;
    }

    /* Read Score and Predecesser */
    if (fscanf(f_score, "#%ld\n", &num_anchor) != 1) return -1;
    assert(num_anchor == in->n);

    /* Read Anchors */
    KMALLOC(km, in->p, in->n);
    KMALLOC(km, in->f, in->n);
    for (int i = 0; i < in->n; i++) {
        if (fscanf(f_score, "%d,%ld", &in->f[i], &in->p[i]) != 2) return -1;
    }
    fscanf(f_score, "\n");
#ifdef DEBUG_VERBOSE
    debug_output_score(debug_folder, in);
#endif  // DEBUG_VERBOSE
    return in->n;
}

#ifdef DEBUG_CHECK_FORCE

/**
 * Build ground truth score using backward_cpu 
 * Allocate and Populate chain_read_t.f, chain_read_t.p
 * Return number of anchors
 */
int debug_build_score(chain_read_t *in, void *km) {
    if (debug_skip_score_check){
        fprintf(stderr,
                "[Info] force score checking against backward_cpu! \n");
    }
    debug_skip_score_check = 0;

    /* Read Anchors */
    KMALLOC(km, in->p, in->n);
    KMALLOC(km, in->f, in->n);

    Misc misc = build_misc(INT64_MAX);
    mg_lchain_dp(misc.max_dist_x, misc.max_dist_y, misc.bw, misc.max_skip,
                 misc.max_iter, misc.min_cnt, misc.min_score, 0.12,
                 misc.chn_pen_skip, misc.is_cdna, misc.n_seg, in->n, in->a,
                 &in->n_u, &in->u, in->km, in, in->f, in->p);
#ifdef DEBUG_VERBOSE
    debug_output_score(debug_folder, in);
#endif  // DEBUG_VERBOSE
    return in->n;
}

#endif // DEBUG_CHECK_FORCE 

/**
 * Check p[], f[] array against gold standard. 
 * Print if there is mismatch. 
*/
int debug_check_score(const int64_t *p, const int32_t *f, const int64_t *p_gold,
                      const int32_t *f_gold, int64_t n, char* qname) {
    for (int64_t i = 0; i < n; ++i) {
        if (f[i] == 0) {
#ifdef DEBUG_VERBOSE
            fprintf(stderr,
                    "[Debug] Score Mismatch: %s Anchor %ld, score: %d (gold "
                    "x), previous: %ld (gold x)\n",
                    qname == 0 ? "--" : qname, i, f[i], p[i]);
#endif
        }
    }
    if (debug_skip_score_check) return -1;
    static int readid = 0;
    size_t score_mismatches = 0;
    int rt = 0;
    for (int64_t i = 0; i < n; ++i) {
        if (p[i] != p_gold[i] || f[i] != f_gold[i]) {
// #ifdef DEBUG_VERBOSE
            fprintf(stderr,
                    "[Debug] Score Mismatch: %s Anchor %ld, score: %d (gold %d), "
                    "previous: %ld (gold %ld)\n",
                    qname == 0? "--": qname, readid, i, f[i], f_gold[i], p[i], p_gold[i]);
// #endif
            rt = 1;
            score_mismatches++;
        }
    }
    if (rt == 1)
        fprintf(stderr, "[Debug] Score Mismatch: %s %d mismatches\n",
                qname == 0 ? "--" : qname, score_mismatches);
    readid++;
    return rt;
}


void debug_check_range(const int32_t* range, size_t n){
    static int read_idx = 0;
    for (size_t i = 1; i < n; i++){
        if (range[i] < range[i-1] - 1)
        fprintf(stderr, "[debug]No realistic range sequence read #%d i %d %d %d\n", read_idx, i, range[i-1], range[i]);
    }
    read_idx++;
}

int debug_check_cut(const size_t *cut, const int32_t *range, size_t max_cut,
                    size_t n, size_t offset) {

    static int read_idx = 0;
    size_t cid = 0;
    for (; cid < max_cut && (cut[cid] < n + offset || cut[cid] == SIZE_MAX);
         cid++) {
        if (cut[cid] != SIZE_MAX) {
            if (cut[cid] != 0 && range[cut[cid] - 1] != 0)
                fprintf(
                    stderr,
                    "[debug] Cut Error: > %ld len %d, Cut at %zu %zu (%d)\n",
                    read_idx, cut[cid], offset, range[cut[cid] - 1]);
        }
        if (cid > 0 && cut[cid] != SIZE_MAX){
            static size_t prev_cut = 0;
            int cut_issue = 0;
            for (size_t i = prev_cut; i < cut[cid]; i++) {
                if (range[i] + i >= cut[cid]){
                    fprintf(stderr, "[debug] Cut Error: > %ld cid %d , Cut %zu - %zu, i %zu, range %zu\n", 
                    read_idx, cid, prev_cut, cut[cid], i, range[i]);
                    cut_issue = 1;
                }
            }
            if (cut_issue){
                for (int i = prev_cut; i < cut[cid]; i++){
                    fprintf(stderr, "%zu[%d]\t", i, range[i]);
                }
                fprintf(stderr, "\n");
            }

            prev_cut = cut[cid];
        }
    }
    read_idx++;
    return cid;
}

#endif  // DEBUG_CHECK