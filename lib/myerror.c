#include <stdio.h>
#include <stdlib.h>

void myerror(error_text)
   char *error_text;
{
   printf("\nCrash and burn dude ...\n");
   printf("%s\n",error_text);
   printf("...now exiting to system...\n");
   
   exit(1);
}
