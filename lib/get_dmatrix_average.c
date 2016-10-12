/* =========================================================================
   GET_DMATRIX_AVERAGE: Jan 16/95

   Average over row or column of a matrix of double values
   =========================================================================
*/
#include "myerror.h"

void get_dmatrix_average(
   double **mat, int i0, int iN, int j0, int jN,
   int rcflg,
   double *x)
{
   int i,j;
   double sum;
   
   if (rcflg==1) { /* average over first index */
      for (j=j0;j<=jN;j++) {
         sum=0.0;
         for (i=i0;i<=iN;i++) sum += mat[i][j];
         x[j]=sum/((double) (iN-i0+1));
      }
   } else if (rcflg==2) { /* average over second index */
      for (i=i0;i<=iN;i++) {
         sum=0.0;
         for (j=j0;j<=jN;j++) sum += mat[i][j];
         x[i]=sum/((double) (jN-j0+1));
      }
   } else {
      printf("rcflg=%d\n",rcflg);
      myerror("GET_DMATRIX_AVERAGE(ERROR): Improper index passed");
   }
}
