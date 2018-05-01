#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include "bseq.h"
#include "minimap.h"
#include "mmpriv.h"
#include "kalloc.h"

#define MM_VERSION "merge_v1.0"
#define INITIAL_NUM_REGS 256

//See : https://github.com/lh3/minimap2/issues/141


//make this inline for performance reasons
void multipart_write(FILE* fp, void *buf, size_t element_size, size_t num_elements){
	size_t ret=fwrite(buf,element_size,num_elements,fp);
	if(ret!=num_elements){
		fprintf(stderr,"Writing error has occured :%s\n",strerror(errno));
		exit(EXIT_FAILURE);
	}
}

static inline void multipart_read(FILE* fp, void *buf, size_t element_size, size_t num_elements){
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


static inline void mm_hit_sort_by_score(void *km, int *n_regs, mm_reg1_t *r)
{
	int32_t i, n_aux, n = *n_regs;
	mm128_t *aux;
	mm_reg1_t *t;

	if (n <= 1) return;
	aux = (mm128_t*)kmalloc(km, n * 16);
	t = (mm_reg1_t*)kmalloc(km, n * sizeof(mm_reg1_t));

	for (i = n_aux = 0; i < n; ++i) {
		aux[n_aux].x = (uint64_t)r[i].score << 32 | r[i].hash;
		aux[n_aux++].y = i;
	}
	radix_sort_128x(aux, aux + n_aux);
	for (i = n_aux - 1; i >= 0; --i)
		t[n_aux - 1 - i] = r[aux[i].y];

	memcpy(r, t, sizeof(mm_reg1_t) * n_aux);
	*n_regs = n_aux;
	kfree(km, aux);
	kfree(km, t);
}

static inline mm_reg1_t *merge_regs(const mm_mapopt_t *opt, const mm_idx_t *mi, void *km, int qlen, int *n_regs, mm_reg1_t *regs)
{
	// if (!(opt->flag & MM_F_ALL_CHAINS)){
	// 	int backup=*n_regs;
	// 	mm_set_parent(km, opt->mask_level, *n_regs, regs, opt->a * 2 + opt->b);
	// 	mm_select_sub(km, opt->pri_ratio, mi->k*2, opt->best_n, n_regs, regs);
	// 	fprintf(stderr,"Filter check %d %d\n",backup,*n_regs);
	// 	mm_set_sam_pri(*n_regs, regs);
	// }

	//if (!(opt->flag & MM_F_CIGAR)) return regs;
	//int backup=*n_regs;
	//mm_filter_regs(km, opt, qlen, n_regs, regs);
	//fprintf(stderr,"Filter check %d %d\n",backup,*n_regs);

	if (opt->flag & MM_F_CIGAR) {
		mm_hit_sort_by_dp(km, n_regs, regs);
	}
	else{
		mm_hit_sort_by_score(km, n_regs, regs);
	}

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


//This function can be parallelised later by a kt-pipeline later if performance is a problem
void merge(mm_mapopt_t *opt, mm_idxopt_t *ipt, int num_idx_parts, const char **fn, int argc, char** argv,const char *rg){

	int i,j,n_reg,n_regs_sum;
	uint32_t n_seq=0,current_seq=0,current_num_regs=INITIAL_NUM_REGS;

	mm_reg1_t *reg = (mm_reg1_t *)malloc(sizeof(mm_reg1_t)*current_num_regs); 
	FILE** file=(FILE**)malloc(sizeof(FILE*)*num_idx_parts);
	mm_idx_t* mi = (mm_idx_t*)calloc(1, sizeof(mm_idx_t));
	char *filename=(char *)malloc(strlen(opt->multi_prefix)+20);
	uint32_t *cum_n_seq=(uint32_t *)malloc(sizeof(uint32_t)*num_idx_parts);
	int *replens=(int *)malloc(sizeof(int)*num_idx_parts);

	//go through each multi-part dump and grab the reference sequence information
	//At the end, mi will have information about all the reference sequences that were in multi-part indices. aka an emulated unipartindex
	for(i=0;i<num_idx_parts;i++){

		//the cumulative number of reference sequences
		//mi->n_seq is at the beginning 0 (calloc used for mi)
		cum_n_seq[i]=mi->n_seq;		

		//generate file name for the multi-part dump and open it
		sprintf(filename,"%s%d.tmp",opt->multi_prefix, i);
		//fprintf(stderr,"filename %s\n",filename);		
		file[i]=fopen(filename,"rb");
		if (file[i]==NULL){
			fprintf(stderr,"ERROR: Cannot open file %s\n",filename);
			exit(EXIT_FAILURE);
		}
		
		//read the number of reference sequnces in the multi-part index
		multipart_read(file[i],&n_seq, sizeof(uint32_t), 1);
		//fprintf(stderr,"We have %d sequences\n",n_seq);

		//update the total number of reference sequnces read so far
		mi->n_seq += n_seq;

		//allocate space for reference sequnces
		mi->seq = (mm_idx_seq_t*)realloc(mi->seq, mi->n_seq*sizeof(mm_idx_seq_t));

		//load each reference sequnce information
		for (j = 0; j < n_seq; ++j) {
			uint8_t l;
			mm_idx_seq_t *s = &mi->seq[current_seq]; current_seq++;
			multipart_read(file[i],&l, 1, 1); //get the read name length
			s->name = (char*)malloc(l + 1);
			multipart_read(file[i],s->name, 1, l);	//get the read name
			s->name[l] = 0;
			multipart_read(file[i],&(s->len), 4, 1);	//reference sequnce length
			//fprintf(stderr,"sequence name : %s\n",s->name);
		}
	}

	//k-mer size
	mi->k=ipt->k;
	//fprintf(stderr,"Kmer size : %d\n",mi->k);

	if (opt->flag & MM_F_OUT_SAM){
		mm_write_sam_hdr(mi, rg, MM_VERSION, argc, argv);
	}

	// fprintf(stderr,"We have %d sequences\n",mi->n_seq);
	// for (j = 0; j < mi->n_seq; ++j) {
	// 	mm_idx_seq_t *s = &mi->seq[j];
	// 	fprintf(stderr,"sequence name : %s\n",s->name);
	// }


	//open the query sequnce/fastq file
	//fprintf(stderr,"opening : %s\n",fn[0]);
	mm_bseq_file_t *fastq=mm_bseq_open(fn[0]);
	if (fastq == 0) {
		if (mm_verbose >= 1)
			fprintf(stderr, "ERROR: failed to open file '%s'\n", fn[0]);
		exit(EXIT_FAILURE);
	}

	int nseq=1; //number of sequences read
	void *km=km_init();	//kmemory
	kstring_t *st;	//kstring
	st = (kstring_t*)calloc(1, sizeof(kstring_t));
	st->s=NULL; st->l=0; st->m=0;

	//read until all the query sequences are read
	while(nseq>0){

		mm_bseq1_t *seq = mm_bseq_read3(fastq, 1, 1, 0, 0, &nseq);	//read a query sequence
		if(nseq!=1 || seq==NULL){	
			break;	//if no query sequences were read thats the end of the file. 
		}
		//fprintf(stderr,"Name %s\n",seq->name);
		int qlen=seq->l_seq; //length of a query sequence

		n_regs_sum=0;	//total number of regs for the current query

		//read the internal state (dumped in the multi-part dump files) for each query seqeunce while going through each part of the index
		for(i=0;i<num_idx_parts;i++){

			multipart_read(file[i],&n_reg,sizeof(int),1); //number of reg (mm_reg1_t)
			//fprintf(stderr,"n regs %d\n",n_reg);

			multipart_read(file[i],&(replens[i]),sizeof(int),1); // the replen : replen is calculated by collect_matches() in map.c. It is the sum of length of regions covered by highly repetitive k-mers.

			//allocate space if we have run out of space for reg
			if(current_num_regs<n_regs_sum+n_reg){
				current_num_regs=current_num_regs+n_regs_sum;
				reg = (mm_reg1_t *)realloc(reg,sizeof(mm_reg1_t)*current_num_regs); 
			}

			//now read the regs one by one
			for (j = 0; j < n_reg; ++j){
				mm_reg1_t *r=&reg[n_regs_sum+j];
				multipart_read(file[i],r,sizeof(mm_reg1_t),1);	//read the mm_reg1_t

				if(opt->flag & MM_F_CIGAR){
					int capacity=0;
					multipart_read(file[i],&capacity,sizeof(uint32_t),1);	//read the capacity of cigar[] under mm_extra_t *p in mm_reg1_t
					
					//this is freed later. not the most optimal way, but for now.
					r->p = (mm_extra_t *)malloc(sizeof(mm_extra_t)+sizeof(uint32_t)*capacity);		//cigar[] is a flexible array member  : allocate memmory as at : https://wiki.sei.cmu.edu/confluence/display/c/MEM33-C.++Allocate+and+copy+structures+containing+a+flexible+array+member+dynamically
					multipart_read(file[i],r->p,sizeof(mm_extra_t)+sizeof(uint32_t)*capacity,1);	//read the mm_extra_t

					if(capacity!=r->p->capacity){
						fprintf(stderr,"WARNING : Multi-part dump files are corrupted\n");
					}
				}
				else{
					r->p=NULL;
				}

				//fix the reference index in the emulated uni-part index by using the cumulative sum of reference sequnces in each part of multi-part index
				r->rid = r->rid+cum_n_seq[i];
			}

			//total number of regs for the current query
			n_regs_sum+=n_reg;
		}


		//fprintf(stderr,"Replens : ");
		//for(i=0;i<num_idx_parts;i++)	fprintf(stderr,"%d\t",replens[i]);

		//perform the merging of regs
		reg=merge_regs(opt, mi, km, qlen, &n_regs_sum, reg);

		int is_sr = !!(opt->flag & MM_F_SR);
		//int rep_len=(int)(maximum(replens,num_idx_parts)*1.1);
		int rep_len=maximum(replens,num_idx_parts);
		//fprintf(stderr,"\nMax replen %d\n",rep_len);
		mm_set_mapq(km, n_regs_sum, reg, opt->min_chain_score, opt->a, rep_len, is_sr);	//set the mapq

		//go through each  reg and print them
		for (j = 0; j < n_regs_sum; ++j) {
			mm_reg1_t *r = &reg[j];
			//fprintf(stderr,"id %d\thash %d\tdiv %f\n",r->id,r->hash,r->div);
			
			assert(!r->sam_pri || r->id == r->parent);
			if ((opt->flag & MM_F_NO_PRINT_2ND) && r->id != r->parent){	//don't print secondary mappings if the respective option has been set
				continue;
			}
			if (opt->flag & MM_F_OUT_SAM){	//sam output
				mm_write_sam2(st, mi, seq, 0, j, 1, &n_regs_sum, (const mm_reg1_t*const*)&reg, km, opt->flag);
				//fprintf(stderr, "%d\t%s\t%d\t%d\t",r->rid, mi->seq[r->rid].name, r->rs+1, r->mapq);
			}
			else{	//paf output
				mm_write_paf(st, mi, seq, r, km, opt->flag);
			}
			mm_err_puts(st->s);	

		}

		if (n_regs_sum == 0 && (opt->flag & MM_F_OUT_SAM)) { // write an unmapped record
			mm_write_sam2(st, mi, seq, 0, -1, 1, &n_regs_sum, (const mm_reg1_t*const*)&reg, km, opt->flag);
			mm_err_puts(st->s);
		}

		for (j = 0; j < n_regs_sum; ++j) {
			mm_reg1_t *r = &reg[j];
			//free the registers p pointer, an optimisation can be done here
			if (r->p!=NULL){ free(r->p); };
		}

		//free the sequence
		free(seq->seq); free(seq->name);
		if (seq->qual) free(seq->qual);
		free(seq);

	}

	//close fastqs
	mm_bseq_close(fastq);

	//check if all multi-part dump files have been read completely
	for(i=0;i<num_idx_parts;i++) {
		if(fread(reg,1,1,file[i])>0){
			fprintf(stderr,"WARNING: Multi-part dump files were not fully read\n");
		}
		else {
			if(feof(file[i])==0){
				fprintf(stderr,"WARNING: Multi-part dump files were not fully read\n");
			}
		}
		fclose(file[i]);
	}	

	// free
	free(st->s); free(st); //free the kstring 

	free(file);
	free(filename);
	free(cum_n_seq);
	free(replens);

	//free the index
	for (j = 0; j < mi->n_seq; ++j) {
			free(mi->seq[j].name);
	}
	free(mi->seq);
	free(mi);	

	free(reg);
	km_destroy(km);

	return;

}
