extern double* dalloc(double* array, size_t xmax);
extern double** dalloc2(double** array, size_t xmax, size_t ymax);
extern void free2(double** array, size_t xmax);
extern double*** dalloc3(double*** array, size_t xmax, size_t ymax, size_t zmax);
extern void free3(double*** array, size_t xmax, size_t ymax);