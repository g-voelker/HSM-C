//
// Created by Georg Sebastian Voelker on 18/11/16.
//

#include <stdlib.h>

// wrappers for 1-, 2- and 3-dimensional double precision allocations / frees
double* dalloc(double* array, size_t xmax) {
  array = malloc(sizeof(double) * xmax);
  return(array);
}

double** dalloc2(double** array, size_t xmax, size_t ymax){
  int nn=0;
  array = (double **) malloc(sizeof(double) * xmax);
  for (nn=0; nn<xmax; nn++){
    array[nn] = (double *) malloc(sizeof(double) * ymax);
  }
  return(array);
}

void dfree2(double** array, size_t xmax){
  int nn=0;
  for (nn=0; nn<xmax; nn++){
    free(array[nn]);
  }
  free(array);
}

double*** dalloc3(double*** array, size_t xmax, size_t ymax, size_t zmax){
  int nn=0, mm=0;
  array = (double ***) malloc(sizeof(double) * xmax);
  for (nn=0; nn<xmax; nn++){
    array[nn] = (double **) malloc(sizeof(double) * ymax);
    for (mm=0; mm<ymax; mm++){
      array[nn][mm] = (double *) malloc(sizeof(double) * zmax);
    }
  }
  return(array);
}

void dfree3(double*** array, size_t xmax, size_t ymax){
  int nn=0, mm=0;
  for (nn=0; nn<xmax; nn++){
    for (mm=0; mm<ymax; mm++){
      free(array[nn][mm]);
    }
    free(array[nn]);
  }
  free(array);
}