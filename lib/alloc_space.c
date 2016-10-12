/* Routines to allocate space for vectors and matrices dynamically */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <malloc.h>
#include "complex.h"

void nrerror(error_text)
   char error_text[];
{
   printf("Numerical Recipes run-time error...\n");
   printf("%s\n",error_text);
   printf("...now exiting to system...\n");
   
   exit(1);
}

short *svector(nl,nh)
   int nl,nh;
{
   short *v;

   v=(short *)malloc((nh-nl+1)*sizeof(short));
   if (!v) nrerror("allocation failure in ivector()");
   return v-nl;
}


void free_svector(v, nl,nh)
   int nl,nh;
   short *v;
{
   free((short*) (v+nl));
}

int *ivector(nl,nh)
   int nl,nh;
{
   int *v;

   v=(int *)malloc((nh-nl+1)*sizeof(int));
   if (!v) nrerror("allocation failure in ivector()");
   return v-nl;
}


void free_ivector(v, nl,nh)
   int nl,nh;
   int *v;
{
   free((int*) (v+nl));
}

float *fvector(nl,nh)
   int nl,nh;
{
   float *v;

   v=(float *)malloc((nh-nl+1)*sizeof(float));
   if (!v) nrerror("allocation failure in fvector()");
   return v-nl;
}

void free_fvector(v, nl,nh)
   int nl,nh;
   float *v; 
{
   free((float*) (v+nl));
}

double *dvector(nl,nh)
   int nl,nh;
{
   double *v;

   v=(double *)malloc((nh-nl+1)*sizeof(double));
   if (!v) 
      nrerror("allocation failure in dvector()");
   return v-nl;
}

void free_dvector(v, nl,nh)
   int nl,nh;
   double *v;
{
   free((double*) (v+nl));
}

dcomplex *dcvector(nl,nh)
   int nl,nh;
{
   dcomplex *v;

   v=(dcomplex *)malloc((nh-nl+1)*sizeof(dcomplex));
   if (!v) 
      nrerror("allocation failure in dcvector()");
   return(v-nl);
}

void free_dcvector(v, nl,nh)
   int nl,nh;
   dcomplex *v;
{
   free((dcomplex*) (v+nl));
}

zomplex *zvector(nl,nh)
   int nl,nh;
{
   zomplex *v;

   v=(zomplex *)malloc((nh-nl+1)*sizeof(zomplex));
   if (!v) 
      nrerror("allocation failure in zvector()");
   return(v-nl);
}

void free_zvector(v, nl,nh)
   int nl,nh;
   zomplex *v;
{
   free((zomplex*) (v+nl));
}

short **smatrix2(nrl,nrh,ncl,nch)
   int nrl,nrh,ncl,nch;
{
   int i;
   short **m;

   m=(short **) malloc((nrh-nrl+1)*sizeof(short*));
   if (!m) nrerror("allocation failure 1 in smatrix2()");
   m -= nrl;

   for(i=nrl;i<=nrh;i++) {
       m[i]=(short *) malloc((nch-ncl+1)*sizeof(short));
       if (!m[i]) nrerror("allocation failure 2 in smatrix2()");
       m[i] -= ncl;
   }
   return m;
}

void free_smatrix2(m, nrl,nrh,ncl,nch)
   int nrl,nrh,ncl,nch;
   short **m;
{
   int i;

   for(i=nrh;i>=nrl;i--) free((short*) (m[i]+ncl));
   free((short*) (m+nrl));
}

int **imatrix2(nrl,nrh,ncl,nch)
   int nrl,nrh,ncl,nch;
{
   int i;
   int **m;

   m=(int **) malloc((nrh-nrl+1)*sizeof(int*));
   if (!m) nrerror("allocation failure 1 in imatrix2()");
   m -= nrl;

   for(i=nrl;i<=nrh;i++) {
       m[i]=(int *) malloc((nch-ncl+1)*sizeof(int));
       if (!m[i]) nrerror("allocation failure 2 in imatrix2()");
       m[i] -= ncl;
   }
   return m;
}

void free_imatrix2(m, nrl,nrh,ncl,nch)
   int nrl,nrh,ncl,nch;
   int **m;
{
   int i;

   for(i=nrh;i>=nrl;i--) free((int*) (m[i]+ncl));
   free((int*) (m+nrl));
}

float **fmatrix2(nrl,nrh,ncl,nch)
   int nrl,nrh,ncl,nch;
{
   int i;
   float **m;

   m=(float **) malloc((nrh-nrl+1)*sizeof(float*));
   if (!m) nrerror("allocation failure 1 in fmatrix2()");
   m -= nrl;

   for(i=nrl;i<=nrh;i++) {
       m[i]=(float *) malloc((nch-ncl+1)*sizeof(float));
       if (!m[i]) nrerror("allocation failure 2 in fmatrix2()");
       m[i] -= ncl;
   }
   return m;
}


void free_fmatrix2(m, nrl,nrh,ncl,nch)
   int nrl,nrh,ncl,nch;
   float **m;
{
   int i;

   for(i=nrh;i>=nrl;i--) free((float*) (m[i]+ncl));
   free((float*) (m+nrl));
}

double **dmatrix2(nrl,nrh,ncl,nch)
   int nrl,nrh,ncl,nch;
{
   int i;
   double **m;

   m=(double **) malloc((nrh-nrl+1)*sizeof(double*));
   if (!m) nrerror("allocation failure 1 in dmatrix2()");
   m -= nrl;

   for(i=nrl;i<=nrh;i++) {
       m[i]=(double *) malloc((nch-ncl+1)*sizeof(double));
       if (!m[i]) nrerror("allocation failure 2 in dmatrix2()");
       m[i] -= ncl;
   }
   return m;
}


void free_dmatrix2(m, nrl,nrh,ncl,nch)
   int nrl,nrh,ncl,nch;
   double **m;
{
   int i;

   for(i=nrh;i>=nrl;i--) free((double*) (m[i]+ncl));
   free((double*) (m+nrl));
}

dcomplex **dcmatrix2(nrl,nrh,ncl,nch)
   int nrl,nrh,ncl,nch;
{
   int i;
   dcomplex **m;

   m=(dcomplex **) malloc((nrh-nrl+1)*sizeof(dcomplex*));
   if (!m) nrerror("allocation failure 1 in dcmatrix2()");
   m -= nrl;

   for(i=nrl;i<=nrh;i++) {
       m[i]=(dcomplex *) malloc((nch-ncl+1)*sizeof(dcomplex));
       if (!m[i]) nrerror("allocation failure 2 in dcmatrix2()");
       m[i] -= ncl;
   }
   return m;
}


void free_dcmatrix2(m, nrl,nrh,ncl,nch)
   int nrl,nrh,ncl,nch;
   dcomplex **m;
{
   int i;

   for(i=nrh;i>=nrl;i--) free((dcomplex*) (m[i]+ncl));
   free((dcomplex*) (m+nrl));
}

zomplex **zmatrix2(nrl,nrh,ncl,nch)
   int nrl,nrh,ncl,nch;
{
   int i;
   zomplex **m;

   m=(zomplex **) malloc((nrh-nrl+1)*sizeof(zomplex*));
   if (!m) nrerror("allocation failure 1 in zmatrix2()");
   m -= nrl;

   for(i=nrl;i<=nrh;i++) {
       m[i]=(zomplex *) malloc((nch-ncl+1)*sizeof(zomplex));
       if (!m[i]) nrerror("allocation failure 2 in zmatrix2()");
       m[i] -= ncl;
   }
   return m;
}

void free_zmatrix2(m, nrl,nrh,ncl,nch)
   int nrl,nrh,ncl,nch;
   zomplex **m;
{
   int i;

   for(i=nrh;i>=nrl;i--) free((zomplex*) (m[i]+ncl));
   free((zomplex*) (m+nrl));
}



int ***imatrix3(nrl,nrh,ncl,nch,ndl,ndh)
   int nrl,nrh,ncl,nch,ndl,ndh;
{
   int i,j;
   int ***m;

   m=(int ***) malloc((nrh-nrl+1)*sizeof(int**));
   if (!m) nrerror("allocation failure 1 in imatrix3()");
   m -= nrl;

   for(i=nrl;i<=nrh;i++) {
       m[i]=(int **) malloc((nch-ncl+1)*sizeof(int*));
       if (!m[i]) nrerror("allocation failure 2 in imatrix3()");
       m[i] -= ncl;

       for(j=ncl;j<=nch;j++) {
          m[i][j]=(int *) malloc((ndh-ndl+1)*sizeof(int));
          if (!m[i][j]) nrerror("allocation failure 3 in imatrix3()");
          m[i][j] -= ndl;
       }
      
   }
   return m;
}


void free_imatrix3(m, nrl,nrh,ncl,nch,ndl,ndh)
   int nrl,nrh,ncl,nch,ndl,ndh;
   int ***m;
{
   int i,j;

   for(i=nrh;i>=nrl;i--) 
   for(j=nch;j>=ncl;j--)
      free((char*) (m[i][j]+ndl));
   for(i=nrh;i>=nrl;i--) 
      free((char*) (m[i]+ncl));
   free((char*) (m+nrl));
}

short ***smatrix3(nrl,nrh,ncl,nch,ndl,ndh)
   int nrl,nrh,ncl,nch,ndl,ndh;
{
   int i,j;
   short ***m;

   m=(short ***) malloc((nrh-nrl+1)*sizeof(short**));
   if (!m) nrerror("allocation failure 1 in imatrix3()");
   m -= nrl;

   for(i=nrl;i<=nrh;i++) {
       m[i]=(short **) malloc((nch-ncl+1)*sizeof(short*));
       if (!m[i]) nrerror("allocation failure 2 in imatrix3()");
       m[i] -= ncl;

       for(j=ncl;j<=nch;j++) {
          m[i][j]=(short *) malloc((ndh-ndl+1)*sizeof(short));
          if (!m[i][j]) nrerror("allocation failure 3 in imatrix3()");
          m[i][j] -= ndl;
       }
      
   }
   return m;
}


void free_smatrix3(m, nrl,nrh,ncl,nch,ndl,ndh)
   int nrl,nrh,ncl,nch,ndl,ndh;
   short ***m;
{
   int i,j;

   for(i=nrh;i>=nrl;i--) 
   for(j=nch;j>=ncl;j--)
      free((short*) (m[i][j]+ndl));
   for(i=nrh;i>=nrl;i--) 
      free((short*) (m[i]+ncl));
   free((short*) (m+nrl));
}

float ***fmatrix3(nrl,nrh,ncl,nch,ndl,ndh)
   int nrl,nrh,ncl,nch,ndl,ndh;
{
   int i,j;
   float ***m;

   m=(float ***) malloc((nrh-nrl+1)*sizeof(float**));
   if (!m) nrerror("allocation failure 1 in fmatrix3()");
   m -= nrl;

   for(i=nrl;i<=nrh;i++) {
       m[i]=(float **) malloc((nch-ncl+1)*sizeof(float*));
       if (!m[i]) nrerror("allocation failure 2 in fmatrix3()");
       m[i] -= ncl;

       for(j=ncl;j<=nch;j++) {
          m[i][j]=(float *) malloc((ndh-ndl+1)*sizeof(float));
          if (!m[i][j]) nrerror("allocation failure 3 in fmatrix3()");
          m[i][j] -= ndl;
       }
      
   }
   return m;
}

void free_fmatrix3(m, nrl,nrh,ncl,nch,ndl,ndh)
   int nrl,nrh,ncl,nch,ndl,ndh;
   float ***m;
{
   int i,j;

   for(i=nrh;i>=nrl;i--) 
   for(j=nch;j>=ncl;j--)
      free((float*) (m[i][j]+ndl));
   for(i=nrh;i>=nrl;i--) 
      free((float*) (m[i]+ncl));
   free((float*) (m+nrl));
}

double ***dmatrix3(nrl,nrh,ncl,nch,ndl,ndh)
   int nrl,nrh,ncl,nch,ndl,ndh;
{
   int i,j;
   double ***m;

   m=(double ***) malloc((nrh-nrl+1)*sizeof(float**));
   if (!m) nrerror("allocation failure 1 in dmatrix3()");
   m -= nrl;

   for(i=nrl;i<=nrh;i++) {
       m[i]=(double **) malloc((nch-ncl+1)*sizeof(double*));
       if (!m[i]) nrerror("allocation failure 2 in dmatrix3()");
       m[i] -= ncl;

       for(j=ncl;j<=nch;j++) {
          m[i][j]=(double *) malloc((ndh-ndl+1)*sizeof(double));
          if (!m[i][j]) nrerror("allocation failure 3 in dmatrix3()");
          m[i][j] -= ndl;
       }
      
   }
   return m;
}

void free_dmatrix3(m, nrl,nrh,ncl,nch,ndl,ndh)
   int nrl,nrh,ncl,nch,ndl,ndh;
   double ***m;
{
   int i,j;

   for(i=nrh;i>=nrl;i--) 
   for(j=nch;j>=ncl;j--)
      free((double*) (m[i][j]+ndl));
   for(i=nrh;i>=nrl;i--) 
      free((double*) (m[i]+ncl));
   free((double*) (m+nrl));
}

dcomplex ***dcmatrix3(nrl,nrh,ncl,nch,ndl,ndh)
   int nrl,nrh,ncl,nch,ndl,ndh;
{
   int i,j;
   dcomplex ***m;

   m=(dcomplex ***) malloc((nrh-nrl+1)*sizeof(float**));
   if (!m) nrerror("allocation failure 1 in dcmatrix3()");
   m -= nrl;

   for(i=nrl;i<=nrh;i++) {
       m[i]=(dcomplex **) malloc((nch-ncl+1)*sizeof(dcomplex*));
       if (!m[i]) nrerror("allocation failure 2 in dcmatrix3()");
       m[i] -= ncl;

       for(j=ncl;j<=nch;j++) {
          m[i][j]=(dcomplex *) malloc((ndh-ndl+1)*sizeof(dcomplex));
          if (!m[i][j]) nrerror("allocation failure 3 in dcmatrix3()");
          m[i][j] -= ndl;
       }
      
   }
   return m;
}

void free_dcmatrix3(m, nrl,nrh,ncl,nch,ndl,ndh)
   int nrl,nrh,ncl,nch,ndl,ndh;
   dcomplex ***m;
{
   int i,j;

   for(i=nrh;i>=nrl;i--) 
   for(j=nch;j>=ncl;j--)
      free((dcomplex*) (m[i][j]+ndl));
   for(i=nrh;i>=nrl;i--) 
      free((dcomplex*) (m[i]+ncl));
   free((dcomplex*) (m+nrl));
}

zomplex ***zmatrix3(nrl,nrh,ncl,nch,ndl,ndh)
   int nrl,nrh,ncl,nch,ndl,ndh;
{
   int i,j;
   zomplex ***m;

   m=(zomplex ***) malloc((nrh-nrl+1)*sizeof(float**));
   if (!m) nrerror("allocation failure 1 in zmatrix3()");
   m -= nrl;

   for(i=nrl;i<=nrh;i++) {
       m[i]=(zomplex **) malloc((nch-ncl+1)*sizeof(zomplex*));
       if (!m[i]) nrerror("allocation failure 2 in zmatrix3()");
       m[i] -= ncl;

       for(j=ncl;j<=nch;j++) {
          m[i][j]=(zomplex *) malloc((ndh-ndl+1)*sizeof(zomplex));
          if (!m[i][j]) nrerror("allocation failure 3 in zmatrix3()");
          m[i][j] -= ndl;
       }
      
   }
   return m;
}

void free_zmatrix3(m, nrl,nrh,ncl,nch,ndl,ndh)
   int nrl,nrh,ncl,nch,ndl,ndh;
   zomplex ***m;
{
   int i,j;

   for(i=nrh;i>=nrl;i--) 
   for(j=nch;j>=ncl;j--)
      free((zomplex*) (m[i][j]+ndl));
   for(i=nrh;i>=nrl;i--) 
      free((zomplex*) (m[i]+ncl));
   free((zomplex*) (m+nrl));
}

float ****fmatrix4(nl1,nh1,nl2,nh2,nl3,nh3,nl4,nh4)
   int nl1,nh1,nl2,nh2,nl3,nh3,nl4,nh4;
{
   int i,j,k;
   float ****m;

   m=(float ****) malloc((nh1-nl1+1)*sizeof(float***));
   if (!m) nrerror("allocation failure 1 in fmatrix4()");
   m -= nl1;

   for(i=nl1;i<=nh1;i++) {
      m[i]=(float ***) malloc((nh2-nl2+1)*sizeof(float**));
      if (!m[i]) nrerror("allocation failure 2 in fmatrix4()");
      m[i] -= nl2;

      for(j=nl2;j<=nh2;j++) {
         m[i][j]=(float **) malloc((nh3-nl3+1)*sizeof(float*));
         if (!m[i][j]) nrerror("allocation failure 3 in fmatrix4()");
         m[i][j] -= nl3;

         for(k=nl3;k<=nh3;k++) {
            m[i][j][k]=(float *) malloc((nh4-nl4+1)*sizeof(float));
            if (!m[i][j][k]) nrerror("allocation failure 4 in fmatrix4()");
            m[i][j][k] -= nl4;
         }
      }
   }
   return m;
}

void free_fmatrix4(m, nl1,nh1,nl2,nh2,nl3,nh3,nl4,nh4)
   int nl1,nh1,nl2,nh2,nl3,nh3,nl4,nh4;
   float ****m;
{
   int i,j,k;

   for(i=nh1;i>=nl1;i--) 
   for(j=nh2;j>=nl2;j--)
   for(k=nh3;k>=nl3;k--)
      free((float*) (m[i][j][k]+nl4));
   for(i=nh1;i>=nl1;i--) 
   for(j=nh2;j>=nl2;j--)
      free((float*) (m[i][j]+nl3));
   for(i=nh1;i>=nl1;i--) 
      free((float*) (m[i]+nl2));
   free((float*) (m+nl1));
}

double ****dmatrix4(nl1,nh1,nl2,nh2,nl3,nh3,nl4,nh4)
   int nl1,nh1,nl2,nh2,nl3,nh3,nl4,nh4;
{
   int i,j,k;
   double ****m;

   m=(double ****) malloc((nh1-nl1+1)*sizeof(double***));
   if (!m) nrerror("allocation failure 1 in fmatrix4()");
   m -= nl1;

   for(i=nl1;i<=nh1;i++) {
      m[i]=(double ***) malloc((nh2-nl2+1)*sizeof(double**));
      if (!m[i]) nrerror("allocation failure 2 in fmatrix4()");
      m[i] -= nl2;

      for(j=nl2;j<=nh2;j++) {
         m[i][j]=(double **) malloc((nh3-nl3+1)*sizeof(double*));
         if (!m[i][j]) nrerror("allocation failure 3 in fmatrix4()");
         m[i][j] -= nl3;

         for(k=nl3;k<=nh3;k++) {
            m[i][j][k]=(double *) malloc((nh4-nl4+1)*sizeof(double));
            if (!m[i][j][k]) nrerror("allocation failure 4 in fmatrix4()");
            m[i][j][k] -= nl4;
         }
      }
   }
   return m;
}

void free_dmatrix4(m, nl1,nh1,nl2,nh2,nl3,nh3,nl4,nh4)
   int nl1,nh1,nl2,nh2,nl3,nh3,nl4,nh4;
   double ****m;
{
   int i,j,k;

   for(i=nh1;i>=nl1;i--) 
   for(j=nh2;j>=nl2;j--)
   for(k=nh3;k>=nl3;k--)
      free((double*) (m[i][j][k]+nl4));
   for(i=nh1;i>=nl1;i--) 
   for(j=nh2;j>=nl2;j--)
      free((double*) (m[i][j]+nl3));
   for(i=nh1;i>=nl1;i--) 
      free((double*) (m[i]+nl2));
   free((double*) (m+nl1));
}

dcomplex ****dcmatrix4(nl1,nh1,nl2,nh2,nl3,nh3,nl4,nh4)
   int nl1,nh1,nl2,nh2,nl3,nh3,nl4,nh4;
{
   int i,j,k;
   dcomplex ****m;

   m=(dcomplex ****) malloc((nh1-nl1+1)*sizeof(dcomplex***));
   if (!m) nrerror("allocation failure 1 in fmatrix4()");
   m -= nl1;

   for(i=nl1;i<=nh1;i++) {
      m[i]=(dcomplex ***) malloc((nh2-nl2+1)*sizeof(dcomplex**));
      if (!m[i]) nrerror("allocation failure 2 in fmatrix4()");
      m[i] -= nl2;

      for(j=nl2;j<=nh2;j++) {
         m[i][j]=(dcomplex **) malloc((nh3-nl3+1)*sizeof(dcomplex*));
         if (!m[i][j]) nrerror("allocation failure 3 in fmatrix4()");
         m[i][j] -= nl3;

         for(k=nl3;k<=nh3;k++) {
            m[i][j][k]=(dcomplex *) malloc((nh4-nl4+1)*sizeof(dcomplex));
            if (!m[i][j][k]) nrerror("allocation failure 4 in fmatrix4()");
            m[i][j][k] -= nl4;
         }
      }
   }
   return m;
}

void free_dcmatrix4(m, nl1,nh1,nl2,nh2,nl3,nh3,nl4,nh4)
   int nl1,nh1,nl2,nh2,nl3,nh3,nl4,nh4;
   dcomplex ****m;
{
   int i,j,k;

   for(i=nh1;i>=nl1;i--) 
   for(j=nh2;j>=nl2;j--)
   for(k=nh3;k>=nl3;k--)
      free((dcomplex*) (m[i][j][k]+nl4));
   for(i=nh1;i>=nl1;i--) 
   for(j=nh2;j>=nl2;j--)
      free((dcomplex*) (m[i][j]+nl3));
   for(i=nh1;i>=nl1;i--) 
      free((dcomplex*) (m[i]+nl2));
   free((dcomplex*) (m+nl1));
}

zomplex ****zmatrix4(nl1,nh1,nl2,nh2,nl3,nh3,nl4,nh4)
   int nl1,nh1,nl2,nh2,nl3,nh3,nl4,nh4;
{
   int i,j,k;
   zomplex ****m;

   m=(zomplex ****) malloc((nh1-nl1+1)*sizeof(zomplex***));
   if (!m) nrerror("allocation failure 1 in fmatrix4()");
   m -= nl1;

   for(i=nl1;i<=nh1;i++) {
      m[i]=(zomplex ***) malloc((nh2-nl2+1)*sizeof(zomplex**));
      if (!m[i]) nrerror("allocation failure 2 in fmatrix4()");
      m[i] -= nl2;

      for(j=nl2;j<=nh2;j++) {
         m[i][j]=(zomplex **) malloc((nh3-nl3+1)*sizeof(zomplex*));
         if (!m[i][j]) nrerror("allocation failure 3 in fmatrix4()");
         m[i][j] -= nl3;

         for(k=nl3;k<=nh3;k++) {
            m[i][j][k]=(zomplex *) malloc((nh4-nl4+1)*sizeof(zomplex));
            if (!m[i][j][k]) nrerror("allocation failure 4 in fmatrix4()");
            m[i][j][k] -= nl4;
         }
      }
   }
   return m;
}

void free_zmatrix4(m, nl1,nh1,nl2,nh2,nl3,nh3,nl4,nh4)
   int nl1,nh1,nl2,nh2,nl3,nh3,nl4,nh4;
   zomplex ****m;
{
   int i,j,k;

   for(i=nh1;i>=nl1;i--) 
   for(j=nh2;j>=nl2;j--)
   for(k=nh3;k>=nl3;k--)
      free((zomplex*) (m[i][j][k]+nl4));
   for(i=nh1;i>=nl1;i--) 
   for(j=nh2;j>=nl2;j--)
      free((zomplex*) (m[i][j]+nl3));
   for(i=nh1;i>=nl1;i--) 
      free((zomplex*) (m[i]+nl2));
   free((zomplex*) (m+nl1));
}
