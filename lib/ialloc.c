//
// Created by Georg Sebastian Voelker on 18/11/16.
//

#include <stdlib.h>

// wrappers for 1-, 2- and 3-dimensional integer allocations / frees
int* ialloc(int* array, size_t xmax) {
  array = malloc(sizeof(int) * xmax);
  return(array);
}

int** ialloc2(int** array, size_t xmax, size_t ymax){
  int nn=0;
  array = (int **) malloc(sizeof(int) * xmax);
  for (nn=0; nn<xmax; nn++){
    array[nn] = (int *) malloc(sizeof(int) * ymax);
  }
  return(array);
}

void ifree2(int** array, size_t xmax){
  int nn=0;
  for (nn=0; nn<xmax; nn++){
    free(array[nn]);
  }
  free(array);
}

int*** ialloc3(int*** array, size_t xmax, size_t ymax, size_t zmax){
  int nn=0, mm=0;
  array = (int ***) malloc(sizeof(int) * xmax);
  for (nn=0; nn<xmax; nn++){
    array[nn] = (int **) malloc(sizeof(int) * ymax);
    for (mm=0; mm<zmax; mm++){
      array[nn][mm] = (int *) malloc(sizeof(int) * zmax);
    }
  }
  return(array);
}

void ifree3(int*** array, size_t xmax, size_t ymax){
  int nn=0, mm=0;
  for (nn=0; nn<xmax; nn++){
    for (mm=0; mm<ymax; mm++){
      free(array[nn][mm]);
    }
    free(array[nn]);
  }
  free(array);
}