
FILE* multipart_init(const mm_mapopt_t *opt, const mm_idx_t *idx);
void multipart_close(FILE* fd);
void multipart_write(FILE* fd, void *buf, size_t element_size, size_t num_elements);
void merge(mm_mapopt_t *opt, int num_idx_parts, const char **fn, int argc, char** argv,const char *rg);