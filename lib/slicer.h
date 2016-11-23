//
// Created by georg on 22/11/16.
//

extern int* islice1(int* array, int* slice, int nxmin, int nxmax);
extern float* fslice1(float* array, float* slice, int nxmin, int nxmax);
extern double* dslice1(double* array, double* slice, int nxmin, int nxmax);

extern int** islice2(int** array, int** slice, int nxmin, int nxmax, int nymin, int nymax);
extern float** fslice2(float** array, float** slice, int nxmin, int nxmax, int nymin, int nymax);
extern double** dslice2(double** array, double** slice, int nxmin, int nxmax, int nymin, int nymax);

extern int*** islice3(int*** array, int*** slice, int nxmin, int nxmax, int nymin, int nymax, int nzmin, int nzmax);
extern float*** fslice3(float*** array, float*** slice, int nxmin, int nxmax, int nymin, int nymax, int nzmin, int nzmax);
extern double*** dslice3(double*** array, double*** slice, int nxmin, int nxmax, int nymin, int nymax, int nzmin, int nzmax);