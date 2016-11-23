extern int* ialloc(int* array, size_t xmax);
extern int** ialloc2(int** array, size_t xmax, size_t ymax);
extern void ifree2(int** array, size_t xmax);
extern int*** ialloc3(int*** array, size_t xmax, size_t ymax, size_t zmax);
extern void ifree3(int*** array, size_t xmax, size_t ymax);