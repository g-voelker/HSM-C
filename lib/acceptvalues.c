/* =================================================================
   Routines to accept input with some checking
   =================================================================
*/
#include <stdio.h>
#include <stdlib.h>

/* ================================================================
   AcceptInteger: (July 9, 1996)

   Program to take prompt for an integer with specified bounds
   ================================================================
*/

void AcceptInteger(int *ivalue, int imin, int imax)
{
  int pass;
  char dum[10];

  *ivalue = imin-1;
  pass=0;
  while ((*ivalue)<imin || (*ivalue)>imax) {
    if (pass++ ==0) {
       printf("Input integer: ");
    } else {
       printf("INVALID CHOICE ... Input integer: ");
    }
    scanf("%d",ivalue);    
  }
  /* gets(dum); */
}
