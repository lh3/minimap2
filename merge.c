#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include "bseq.h"
#include "minimap.h"
#include "mmpriv.h"

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

FILE* multipart_init(const mm_mapopt_t *opt, const mm_idx_t *mi){

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
	multipart_write(fp,&mi->n_seq, sizeof(mi->n_seq), 1);
	for (i = 0; i < mi->n_seq; ++i) {
		uint8_t l;
		l = strlen(mi->seq[i].name);
		multipart_write(fp,&l, 1, 1);
		multipart_write(fp,mi->seq[i].name, 1, l);
	}	

	return fp;		
}

void multipart_close(FILE* fp){
	int ret=fclose(fp);
	if (ret==-1){
		fprintf(stderr, "Cannot close file");
	}
}




void merge(mm_mapopt_t *opt, int num_idx_parts, const char **fn){


	int i,j;
	int n_reg;
	int n_regs_sum;
	uint32_t n_seq=0;
	uint32_t current_seq=0;
	char filename[256];

	mm_reg1_t *reg = (mm_reg1_t *)malloc(sizeof(mm_reg1_t)*1000); //free this
	FILE** file=(FILE**)malloc(sizeof(FILE*)*num_idx_parts);
	mm_idx_t* mi = (mm_idx_t*)calloc(1, sizeof(mm_idx_t));


	//go through each index
	for(i=0;i<num_idx_parts;i++){		
		sprintf(filename,"%s%d.tmp",opt->multi_prefix, i);
		//fprintf(stderr,"filename %s\n",filename);		
		file[i]=fopen(filename,"rb");
		if (file[i]==NULL)	fprintf(stderr,"Cannot open file %s\n",filename);
		
		//read the number of sequences and then the actuall number of sequences
		if (fread(&n_seq, sizeof(uint32_t), 1, file[i]) != 1) fprintf(stderr,"Reading failed for %s\n",filename);
		//fprintf(stderr,"We have %d sequences\n",n_seq);
		mi->n_seq += n_seq;

		mi->seq = (mm_idx_seq_t*)realloc(mi->seq, mi->n_seq*sizeof(mm_idx_seq_t));
		for (j = 0; j < n_seq; ++j) {
			uint8_t l;
			mm_idx_seq_t *s = &mi->seq[current_seq]; current_seq++;
			if (fread(&l, 1, 1, file[i]) != 1) fprintf(stderr,"Reading failed for %s\n",filename); //read name length
			s->name = (char*)malloc(l + 1);
			if(fread(s->name, 1, l, file[i])!=l) fprintf(stderr,"Reading failed for %s\n",filename);;
			s->name[l] = 0;
			//fprintf(stderr,"sequence name : %s\n",s->name);
		}
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

	
	// mi->n_seq=100;
	// mi->seq = (mm_idx_seq_t*)calloc(mi->n_seq, sizeof(mm_idx_seq_t));
	// for (i = 0; i < mi->n_seq; ++i) {
	// 	uint8_t l=10;
	// 	mi->seq[i].name = (char*)malloc(l + 1);
	// 	mi->seq[i].name[l] = 0;
	// 	sprintf(mi->seq[i].name, "chr%d",i);
	// }

	while(nseq>0){
		mm_bseq1_t *seq = mm_bseq_read3(fastq, 1, 1, 1, 0, &nseq);

		if(nseq!=1 || seq==NULL){
			//fprintf(stderr,"Read %d sequenced when 1 is expected\n",nseq);
			break;
		}
		//fprintf(stderr,"Name %s\n",seq->name);


		n_regs_sum=0;
		for(i=0;i<num_idx_parts;i++){
			int ret=fread(&n_reg,sizeof(int),1,file[i]);
			fprintf(stderr,"n regs %d\n",n_reg);
			if(ret!=1){
				fprintf(stderr,"Reading error");
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

			}
			n_regs_sum+=n_reg;
		}

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
			//else
				//mm_write_paf(st, mi, t, r, km, opt->flag);
			mm_err_puts(st->s);		
		}

		if (n_regs_sum == 0 && (opt->flag & MM_F_OUT_SAM)) { // write an unmapped record
			mm_write_sam2(st, mi, seq, 0, -1, 1, &n_regs_sum, (const mm_reg1_t*const*)&reg, km, opt->flag);
			mm_err_puts(st->s);
		}

		

	}

	mm_bseq_close(fastq);

	for(i=0;i<num_idx_parts;i++){	
		int ret=fclose(file[i]);
		if (ret==-1){
			fprintf(stderr, "Cannot close file descripter");
		}
	}

	// 	//mm_set_mapq(km, n_reg, reg, opt->min_chain_score, opt->a, rep_len, is_sr);
	// 	//multi segments are not supported
	// 	int n_segs=1;
	// 	if (argc - (optind + 1) > n_segs) {
	// 		fprintf(stderr,"Multi segments are not yet supported for merging %s\n");
	// 		return -1;
	// 	}
	// 	pipeline_t pl;
	// 	pipeline_t *p=&pl;

	// 	kstring_t *st;
	// 	st = (kstring_t*)calloc(1, sizeof(kstring_t));

	// 	memset(p, 0, sizeof(pipeline_t));			
	// 	pl.fp = (mm_bseq_file_t**)calloc(n_segs, sizeof(mm_bseq_file_t*));
	// 	for (i = 0; i < n_segs; ++i) {
	// 		pl.fp[i] = mm_bseq_open(argv[optind + 1]);
	// 		if (pl.fp[i] == 0) {
	// 			if (mm_verbose >= 1)
	// 				fprintf(stderr, "ERROR: failed to open file '%s'\n", argv[optind + 1]);
	// 			for (j = 0; j < i; ++j)
	// 				mm_bseq_close(pl.fp[j]);
	// 			free(pl.fp);
	// 			return -1;
	// 		}
	// 	}

//       	step_t *s;
//       	s = (step_t*)calloc(1, sizeof(step_t));
	// 	s->seq = mm_bseq_read3(p->fp[0], 100, 1, 1, 1, &s->n_seq);			
	// 	mm_bseq1_t *t = &s->seq[i];

	// 	for (j = 0; j < n_reg; ++j) {
			
	// 		mm_reg1_t *r = &reg[j];
	// 		assert(!r->sam_pri || r->id == r->parent);
	// 		if ((opt->flag & MM_F_NO_PRINT_2ND) && r->id != r->parent)
	// 			continue;
	// 		if (opt->flag & MM_F_OUT_SAM)
	// 			mm_write_sam2(st, mi, t, 0, j, 1, &n_reg, (const mm_reg1_t*const*)&reg, km, opt->flag);
	// 		else
	// 			mm_write_paf(st, mi, t, r, km, opt->flag);
	// 		mm_err_puts(st->s);	

	// 	}
	// 	if (s->n_reg[i] == 0 && (p->opt->flag & MM_F_OUT_SAM)) { // write an unmapped record
	// 		mm_write_sam2(st, mi, t, 0, -1, 1, &n_reg, (const mm_reg1_t*const*)&reg, km, opt->flag);
	// 		mm_err_puts(st->s);
	// 	}

	// 	free(st->s);
	// 	for (i = 0; i < n_segs; ++i)
	// 		mm_bseq_close(pl.fp[i]);						
				
	// }

	
	// free
	free(st->s); free(st);
	free(reg);
	free(file);
	
	return;

}