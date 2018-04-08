
int multipart_open(const mm_mapopt_t *opt, const mm_idx_t *idx);
void multipart_close(int fd);
void multipart_write(int fd, const void *buf, size_t count);
void merge(mm_mapopt_t *opt, int num_idx_parts, const char **fn);