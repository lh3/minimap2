#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include "bseq.h"
#include "minimap.h"
#include "mmpriv.h"

#define MM_VERSION "2.10-r763-dirty"
#define INITIAL_NUM_REGS 256

void multipart_write(FILE* fp, void *buf, size_t element_size, size_t num_elements){
	size_t ret=fwrite(buf,element_size,num_elements,fp);
	if(ret!=num_elements){
		fprintf(stderr,"Writing error has occured :%s\n",strerror(errno));
		exit(EXIT_FAILURE);
	}
}

void multipart_read(FILE* fp, void *buf, size_t element_size, size_t num_elements){
	size_t ret=fread(buf,element_size,num_elements,fp);
	if(ret!=num_elements){
		fprintf(stderr,"Reading error has occured :%s\n",strerror(errno));
		exit(EXIT_FAILURE);
	}
}

FILE* multipart_init(const mm_mapopt_t *opt, mm_idx_t *mi){

	int i=0;
	char filename[256];
	sprintf(filename,"%s%d.tmp",opt->multi_prefix, mi->idx_id);
	fprintf(stderr,"filename %s\n",filename);
	//int fd=open(filename,O_WRONLY|O_CREAT , S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH);
	FILE *fp=fopen(filename,"wb");
	if (fp==NULL){
		if (mm_verbose >= 1)
			fprintf(stderr, "ERROR: failed to open file '%s'\n for writing", opt->multi_prefix);
		exit(EXIT_FAILURE);
	}

	//write sequence names
	multipart_write(fp,&(mi->n_seq), sizeof(mi->n_seq), 1);
	for (i = 0; i < mi->n_seq; ++i) {
		uint8_t l;
		l = strlen(mi->seq[i].name);
		multipart_write(fp,&l, 1, 1);
		multipart_write(fp,mi->seq[i].name, 1, l);
		multipart_write(fp,&(mi->seq[i].len), 4, 1);
	}	

	return fp;		
}

void multipart_close(FILE* fp){
	int ret=fclose(fp);
	if (ret==-1){
		fprintf(stderr, "Cannot close file");
	}
}



static inline mm_reg1_t *merge_regs(const mm_mapopt_t *opt, const mm_idx_t *mi, void *km, int qlen, int *n_regs, mm_reg1_t *regs)
{
	if (!(opt->flag & MM_F_CIGAR)) return regs;
	mm_filter_regs(km, opt, qlen, n_regs, regs);
	mm_hit_sort_by_dp(km, n_regs, regs);

	if (!(opt->flag & MM_F_ALL_CHAINS)) { // don't choose primary mapping(s)
		mm_set_parent(km, opt->mask_level, *n_regs, regs, opt->a * 2 + opt->b);
		mm_select_sub(km, opt->pri_ratio, mi->k*2, opt->best_n, n_regs, regs);
		mm_set_sam_pri(*n_regs, regs);
	}

	return regs;
}


static inline int maximum(int *replens,int num_idx_parts){
	assert(num_idx_parts>0);
	int max_rep=replens[0];
	int i;
	for(i=1;i<num_idx_parts;i++){
		assert(replens[i]>=0);
		if(replens[i]>max_rep) max_rep=replens[i];
	}
	return max_rep;
}


void merge(mm_mapopt_t *opt, int num_idx_parts, const char **fn, int argc, char** argv,const char *rg){

	int i,j;
	int n_reg;
	int n_regs_sum;
	uint32_t n_seq=0;
	uint32_t current_seq=0;
	uint32_t current_num_regs=INITIAL_NUM_REGS;

	mm_reg1_t *reg = (mm_reg1_t *)malloc(sizeof(mm_reg1_t)*current_num_regs); 
	FILE** file=(FILE**)malloc(sizeof(FILE*)*num_idx_parts);
	mm_idx_t* mi = (mm_idx_t*)calloc(1, sizeof(mm_idx_t));
	char *filename=(char *)malloc(strlen(opt->multi_prefix)+20);
	uint32_t *cum_n_seq=(uint32_t *)malloc(sizeof(uint32_t)*num_idx_parts);
	int *replens=(int *)malloc(sizeof(int)*num_idx_parts);

	//go through each index
	for(i=0;i<num_idx_parts;i++){
		cum_n_seq[i]=mi->n_seq;		
		sprintf(filename,"%s%d.tmp",opt->multi_prefix, i);
		//fprintf(stderr,"filename %s\n",filename);		
		file[i]=fopen(filename,"rb");
		if (file[i]==NULL)	fprintf(stderr,"Cannot open file %s\n",filename);
		
		//read the number of sequences and then the actuall number of sequences
		multipart_read(file[i],&n_seq, sizeof(uint32_t), 1);
		//fprintf(stderr,"We have %d sequences\n",n_seq);
		mi->n_seq += n_seq;

		mi->seq = (mm_idx_seq_t*)realloc(mi->seq, mi->n_seq*sizeof(mm_idx_seq_t));
		for (j = 0; j < n_seq; ++j) {
			uint8_t l;
			mm_idx_seq_t *s = &mi->seq[current_seq]; current_seq++;
			multipart_read(file[i],&l, 1, 1); //read name length
			s->name = (char*)malloc(l + 1);
			multipart_read(file[i],s->name, 1, l);
			s->name[l] = 0;
			multipart_read(file[i],&(s->len), 4, 1);
			//fprintf(stderr,"sequence name : %s\n",s->name);
		}
	}

	if (opt->flag & MM_F_OUT_SAM){
		mm_write_sam_hdr(mi, rg, MM_VERSION, argc, argv);
	}

	// fprintf(stderr,"We have %d sequences\n",mi->n_seq);
	// for (j = 0; j < mi->n_seq; ++j) {
	// 	mm_idx_seq_t *s = &mi->seq[j];
	// 	fprintf(stderr,"sequence name : %s\n",s->name);
	// }



	fprintf(stderr,"opening : %s\n",fn[0]);
	mm_bseq_file_t *fastq=mm_bseq_open(fn[0]);
	if (fastq == 0) {
		if (mm_verbose >= 1)
			fprintf(stderr, "ERROR: failed to open file '%s'\n", fn[0]);
		exit(EXIT_FAILURE);
	}

	int nseq=1;
	void *km=NULL;
	kstring_t *st;
	st = (kstring_t*)calloc(1, sizeof(kstring_t));
	st->s=NULL; st->l=0; st->m=0;

	while(nseq>0){
		mm_bseq1_t *seq = mm_bseq_read3(fastq, 1, 1, 1, 0, &nseq);
		if(nseq!=1 || seq==NULL){
			//fprintf(stderr,"Read %d sequenced when 1 is expected\n",nseq);
			break;
		}
		//fprintf(stderr,"Name %s\n",seq->name);
		int qlen=seq->l_seq;

		n_regs_sum=0;
		for(i=0;i<num_idx_parts;i++){
			int ret=fread(&n_reg,sizeof(int),1,file[i]);
			fprintf(stderr,"n regs %d\n",n_reg);
			if(ret!=1){
				fprintf(stderr,"Reading error for dump file for index part %d\n",i);
			}
			multipart_read(file[i],&(replens[i]),sizeof(int),1);

			if(current_num_regs<n_regs_sum+n_reg){
				current_num_regs=current_num_regs+n_regs_sum;
				reg = (mm_reg1_t *)realloc(reg,sizeof(mm_reg1_t)*current_num_regs); 
			}

			for (j = 0; j < n_reg; ++j){
				mm_reg1_t *r=&reg[n_regs_sum+j];
				multipart_read(file[i],r,sizeof(mm_reg1_t),1);
				int capacity=0;
				multipart_read(file[i],&capacity,sizeof(uint32_t),1);
				r->p = (mm_extra_t *)malloc(sizeof(mm_extra_t)+sizeof(uint32_t)*capacity);
				multipart_read(file[i],r->p,sizeof(mm_extra_t)+sizeof(uint32_t)*capacity,1);

				if(capacity!=r->p->capacity){
					fprintf(stderr,"Corrupted files?\n");
					exit(EXIT_FAILURE);
				}

				//fix the reference index
				r->rid = r->rid+cum_n_seq[i];
			}
			n_regs_sum+=n_reg;
		}


		fprintf(stderr,"Replens : ");
		for(i=0;i<num_idx_parts;i++)	fprintf(stderr,"%d\t",replens[i]);

		reg=merge_regs(opt, mi, km, qlen, &n_regs_sum, reg);

		int is_sr = !!(opt->flag & MM_F_SR);
		int rep_len=maximum(replens,num_idx_parts);
		fprintf(stderr,"\nMax replen %d\n",rep_len);
		mm_set_mapq(km, n_regs_sum, reg, opt->min_chain_score, opt->a, rep_len, is_sr);

		for (j = 0; j < n_regs_sum; ++j) {
			mm_reg1_t *r = &reg[j];
			fprintf(stderr,"sizeof mm_reg1_t is %ld\t id %d\thash %d\tdiv %f\n",sizeof(mm_reg1_t),r->id,r->hash,r->div);
			
			assert(!r->sam_pri || r->id == r->parent);
			if ((opt->flag & MM_F_NO_PRINT_2ND) && r->id != r->parent)
				continue;
			if (opt->flag & MM_F_OUT_SAM)
				//void mm_write_sam2(kstring_t *s, const mm_idx_t *mi, const mm_bseq1_t *t, int seg_idx, int reg_idx, int n_seg, const int *n_regss, const mm_reg1_t *const* regss, void *km, int opt_flag)
				mm_write_sam2(st, mi, seq, 0, j, 1, &n_regs_sum, (const mm_reg1_t*const*)&reg, km, opt->flag);
				//fprintf(stderr, "%d\t%s\t%d\t%d\t",r->rid, mi->seq[r->rid].name, r->rs+1, r->mapq);
			else
				mm_write_paf(st, mi, seq, r, km, opt->flag);
			mm_err_puts(st->s);		
		}

		if (n_regs_sum == 0 && (opt->flag & MM_F_OUT_SAM)) { // write an unmapped record
			mm_write_sam2(st, mi, seq, 0, -1, 1, &n_regs_sum, (const mm_reg1_t*const*)&reg, km, opt->flag);
			mm_err_puts(st->s);
		}

	}

	//close open files
	mm_bseq_close(fastq);
	for(i=0;i<num_idx_parts;i++) fclose(file[i]);


	// free
	free(st->s); free(st);
	free(reg);
	free(file);
	free(mi);
	free(filename);
	
	return;

}