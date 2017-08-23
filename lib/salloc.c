//
// Created by Georg Sebastian Voelker on 8/23/2017.
//

#include <stdlib.h>

// wrapper for allocation / free of array of strings
char** salloc2(char** array, size_t nmax, size_t smax){
  int nn=0;
  array = (char **) malloc(sizeof(char *) * nmax);
  for (nn=0; nn<nmax; nn++){
    array[nn] = (char *) malloc(sizeof(char) * smax);
  }
  return(array);
}

void sfree2(char** array, size_t nmax){
  int nn=0;
  for (nn=0; nn<nmax; nn++){
    free(array[nn]);
  }
  free(array);
}