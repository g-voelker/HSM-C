//
// Created by georg on 22/11/16.
//

// 1-dimensional slicer
int* islice1(int* array, int* slice, int nxmin, int nxmax){
  int nn;
  for (nn=0; nn<(nxmax - nxmin); nn++){
    slice[nn] = array[nxmin + nn];
  }
  return(slice);
}

float* fslice1(float* array, float* slice, int nxmin, int nxmax){
  int nn;
  for (nn=0; nn<(nxmax - nxmin); nn++){
    slice[nn] = array[nxmin + nn];
  }
  return(slice);
}

double* dslice1(double* array, double* slice, int nxmin, int nxmax){
  int nn;
  for (nn=0; nn<(nxmax - nxmin); nn++){
    slice[nn] = array[nxmin + nn];
  }
  return(slice);
}

// 2-dimensional slicer
int** islice2(int** array, int** slice, int nxmin, int nxmax, int nymin, int nymax){
  int nn, mm;
  for (nn=0; nn<(nxmax - nxmin); nn++){
    for (mm=0; mm<(nymax - nymin); mm++){
      slice[nn][mm] = array[nxmin + nn][nymin + mm];
    }
  }
  return(slice);
}

float** fslice2(float** array, float** slice, int nxmin, int nxmax, int nymin, int nymax){
  int nn, mm;
  for (nn=0; nn<(nxmax - nxmin); nn++){
    for (mm=0; mm<(nymax - nymin); mm++){
      slice[nn][mm] = array[nxmin + nn][nymin + mm];
    }
  }
  return(slice);
}

double** dslice2(double** array, double** slice, int nxmin, int nxmax, int nymin, int nymax){
  int nn, mm;
  for (nn=0; nn<(nxmax - nxmin); nn++){
    for (mm=0; mm<(nymax - nymin); mm++){
      slice[nn][mm] = array[nxmin + nn][nymin + mm];
    }
  }
  return(slice);
}

// 3-dimensional slicer
int*** islice3(int*** array, int*** slice, int nxmin, int nxmax, int nymin, int nymax, int nzmin, int nzmax){
  int ll, nn, mm;
  for (ll=0; ll<(nxmax - nxmin); ll++){
    for (mm=0; mm<(nymax - nymin); mm++){
      for (nn=0; nn<(nzmax - nzmin); nn++){
        slice[ll][nn][mm] = array[nxmin + ll][nymin + mm][nzmin + nn];
      }
    }
  }
  return(slice);
}

float*** fslice3(float*** array, float*** slice, int nxmin, int nxmax, int nymin, int nymax, int nzmin, int nzmax){
  int ll, nn, mm;
  for (ll=0; ll<(nxmax - nxmin); ll++){
    for (mm=0; mm<(nymax - nymin); mm++){
      for (nn=0; nn<(nzmax - nzmin); nn++){
        slice[ll][nn][mm] = array[nxmin + ll][nymin + mm][nzmin + nn];
      }
    }
  }
  return(slice);
}

double*** dslice3(double*** array, double*** slice, int nxmin, int nxmax, int nymin, int nymax, int nzmin, int nzmax){
  int ll, nn, mm;
  for (ll=0; ll<(nxmax - nxmin); ll++){
    for (mm=0; mm<(nymax - nymin); mm++){
      for (nn=0; nn<(nzmax - nzmin); nn++){
        slice[ll][nn][mm] = array[nxmin + ll][nymin + mm][nzmin + nn];
      }
    }
  }
  return(slice);
}