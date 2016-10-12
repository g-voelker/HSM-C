/* ===========================================================================
   Routines to store data in Postscript, RGB, etc format
   Last Update: June 16/98
   ===========================================================================
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/* #include <gl/image.h> */
#include "complex.h"
#include "alloc_space.h"
#include "readPlot.h"
#include "constants.h"
#include "macros.h"
#include "jpeglib.h"

typedef struct RGBCLR {short r, g, b;} rgbclr;

#define CHARDIMEN 8
#define MAXCHARLEN 30
#define MAXCOLMAP 33
#define COLMAPSIZE 256
#define RGBTYPE 1   /* =1 for 1 byte arrays (0..255); =2 for 2 byte arrays */
#define BWRGB 3     /* =1 for black and white; =3 for rgb */
#define BLK 0  
#define WHT 255

#define JPGQUALITY 95

/* ===========================================================================
   CHAR2IMG: MAY 03

   Define 8x8 matrices of letters, numbers and special symbols.
   NOTE: This routine assumes COLMAPSIZE=256
   ===========================================================================
*/
void char2img(c,cm)
   short **cm;
   char c;
{
   int idx;
   short i,j;
   short capletter[26][CHARDIMEN][CHARDIMEN] = 
                  {{{WHT, WHT, WHT, BLK, WHT, WHT, WHT, WHT},
                    {WHT, WHT, BLK, WHT, BLK, WHT, WHT, WHT},
                    {WHT, BLK, WHT, WHT, WHT, BLK, WHT, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, BLK, WHT},
                    {BLK, BLK, BLK, BLK, BLK, BLK, BLK, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, BLK, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, BLK, WHT},
                    {WHT, WHT, WHT, WHT, WHT, WHT, WHT, WHT}},
                   {{BLK, BLK, BLK, BLK, BLK, BLK, WHT, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, BLK, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, BLK, WHT},
                    {BLK, BLK, BLK, BLK, BLK, BLK, WHT, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, BLK, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, BLK, WHT},
                    {BLK, BLK, BLK, BLK, BLK, BLK, WHT, WHT},
                    {WHT, WHT, WHT, WHT, WHT, WHT, WHT, WHT}},
                   {{WHT, BLK, BLK, BLK, BLK, BLK, WHT, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, BLK, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, WHT, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, WHT, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, WHT, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, BLK, WHT},
                    {WHT, BLK, BLK, BLK, BLK, BLK, WHT, WHT},
                    {WHT, WHT, WHT, WHT, WHT, WHT, WHT, WHT}},
                   {{BLK, BLK, BLK, BLK, BLK, BLK, WHT, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, BLK, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, BLK, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, BLK, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, BLK, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, BLK, WHT},
                    {BLK, BLK, BLK, BLK, BLK, BLK, WHT, WHT},
                    {WHT, WHT, WHT, WHT, WHT, WHT, WHT, WHT}},
                   {{BLK, BLK, BLK, BLK, BLK, BLK, BLK, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, WHT, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, WHT, WHT},
                    {BLK, BLK, BLK, BLK, BLK, WHT, WHT, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, WHT, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, WHT, WHT},
                    {BLK, BLK, BLK, BLK, BLK, BLK, BLK, WHT},
                    {WHT, WHT, WHT, WHT, WHT, WHT, WHT, WHT}},
                   {{BLK, BLK, BLK, BLK, BLK, BLK, BLK, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, WHT, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, WHT, WHT},
                    {BLK, BLK, BLK, BLK, BLK, WHT, WHT, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, WHT, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, WHT, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, WHT, WHT},
                    {WHT, WHT, WHT, WHT, WHT, WHT, WHT, WHT}},
                   {{WHT, BLK, BLK, BLK, BLK, BLK, WHT, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, BLK, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, WHT, WHT},
                    {BLK, WHT, WHT, BLK, BLK, BLK, BLK, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, BLK, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, BLK, WHT},
                    {WHT, BLK, BLK, BLK, BLK, BLK, WHT, WHT},
                    {WHT, WHT, WHT, WHT, WHT, WHT, WHT, WHT}},
                   {{BLK, WHT, WHT, WHT, WHT, WHT, BLK, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, BLK, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, BLK, WHT},
                    {BLK, BLK, BLK, BLK, BLK, BLK, BLK, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, BLK, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, BLK, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, BLK, WHT},
                    {WHT, WHT, WHT, WHT, WHT, WHT, WHT, WHT}},
                   {{WHT, WHT, BLK, BLK, BLK, WHT, WHT, WHT},
                    {WHT, WHT, WHT, BLK, WHT, WHT, WHT, WHT},
                    {WHT, WHT, WHT, BLK, WHT, WHT, WHT, WHT},
                    {WHT, WHT, WHT, BLK, WHT, WHT, WHT, WHT},
                    {WHT, WHT, WHT, BLK, WHT, WHT, WHT, WHT},
                    {WHT, WHT, WHT, BLK, WHT, WHT, WHT, WHT},
                    {WHT, WHT, BLK, BLK, BLK, WHT, WHT, WHT},
                    {WHT, WHT, WHT, WHT, WHT, WHT, WHT, WHT}},
                   {{WHT, WHT, WHT, WHT, WHT, WHT, BLK, WHT},
                    {WHT, WHT, WHT, WHT, WHT, WHT, BLK, WHT},
                    {WHT, WHT, WHT, WHT, WHT, WHT, BLK, WHT},
                    {WHT, WHT, WHT, WHT, WHT, WHT, BLK, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, BLK, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, BLK, WHT},
                    {WHT, BLK, BLK, BLK, BLK, BLK, WHT, WHT},
                    {WHT, WHT, WHT, WHT, WHT, WHT, WHT, WHT}},
                   {{BLK, WHT, WHT, WHT, WHT, BLK, WHT, WHT},
                    {BLK, WHT, WHT, WHT, BLK, WHT, WHT, WHT},
                    {BLK, WHT, WHT, BLK, WHT, WHT, WHT, WHT},
                    {BLK, BLK, BLK, WHT, WHT, WHT, WHT, WHT},
                    {BLK, WHT, WHT, BLK, WHT, WHT, WHT, WHT},
                    {BLK, WHT, WHT, WHT, BLK, WHT, WHT, WHT},
                    {BLK, WHT, WHT, WHT, WHT, BLK, WHT, WHT},
                    {WHT, WHT, WHT, WHT, WHT, WHT, WHT, WHT}},
                   {{BLK, WHT, WHT, WHT, WHT, WHT, WHT, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, WHT, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, WHT, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, WHT, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, WHT, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, WHT, WHT},
                    {BLK, BLK, BLK, BLK, BLK, BLK, BLK, WHT},
                    {WHT, WHT, WHT, WHT, WHT, WHT, WHT, WHT}},
                   {{BLK, WHT, WHT, WHT, WHT, WHT, BLK, WHT},
                    {BLK, BLK, WHT, WHT, WHT, BLK, BLK, WHT},
                    {BLK, WHT, BLK, WHT, BLK, WHT, BLK, WHT},
                    {BLK, WHT, WHT, BLK, WHT, WHT, BLK, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, BLK, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, BLK, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, BLK, WHT},
                    {WHT, WHT, WHT, WHT, WHT, WHT, WHT, WHT}},
                   {{BLK, WHT, WHT, WHT, WHT, WHT, BLK, WHT},
                    {BLK, BLK, WHT, WHT, WHT, WHT, BLK, WHT},
                    {BLK, WHT, BLK, WHT, WHT, WHT, BLK, WHT},
                    {BLK, WHT, WHT, BLK, WHT, WHT, BLK, WHT},
                    {BLK, WHT, WHT, WHT, BLK, WHT, BLK, WHT},
                    {BLK, WHT, WHT, WHT, WHT, BLK, BLK, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, BLK, WHT},
                    {WHT, WHT, WHT, WHT, WHT, WHT, WHT, WHT}},
                   {{BLK, BLK, BLK, BLK, BLK, BLK, BLK, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, BLK, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, BLK, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, BLK, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, BLK, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, BLK, WHT},
                    {BLK, BLK, BLK, BLK, BLK, BLK, BLK, WHT},
                    {WHT, WHT, WHT, WHT, WHT, WHT, WHT, WHT}},
                   {{BLK, BLK, BLK, BLK, BLK, BLK, WHT, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, BLK, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, BLK, WHT},
                    {BLK, BLK, BLK, BLK, BLK, BLK, WHT, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, WHT, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, WHT, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, WHT, WHT},
                    {WHT, WHT, WHT, WHT, WHT, WHT, WHT, WHT}},
                   {{WHT, BLK, BLK, BLK, BLK, BLK, WHT, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, BLK, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, BLK, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, BLK, WHT},
                    {BLK, WHT, WHT, WHT, BLK, WHT, BLK, WHT},
                    {BLK, WHT, WHT, WHT, WHT, BLK, WHT, WHT},
                    {WHT, BLK, BLK, BLK, BLK, WHT, BLK, WHT},
                    {WHT, WHT, WHT, WHT, WHT, WHT, WHT, WHT}},
                   {{BLK, BLK, BLK, BLK, BLK, BLK, WHT, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, BLK, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, BLK, WHT},
                    {BLK, BLK, BLK, BLK, BLK, BLK, WHT, WHT},
                    {BLK, WHT, WHT, WHT, BLK, WHT, WHT, WHT},
                    {BLK, WHT, WHT, WHT, WHT, BLK, WHT, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, BLK, WHT},
                    {WHT, WHT, WHT, WHT, WHT, WHT, WHT, WHT}},
                   {{WHT, BLK, BLK, BLK, BLK, BLK, WHT, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, BLK, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, WHT, WHT},
                    {WHT, BLK, BLK, BLK, BLK, BLK, WHT, WHT},
                    {WHT, WHT, WHT, WHT, WHT, WHT, BLK, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, BLK, WHT},
                    {WHT, BLK, BLK, BLK, BLK, BLK, WHT, WHT},
                    {WHT, WHT, WHT, WHT, WHT, WHT, WHT, WHT}},
                   {{BLK, BLK, BLK, BLK, BLK, BLK, BLK, WHT},
                    {WHT, WHT, WHT, BLK, WHT, WHT, WHT, WHT},
                    {WHT, WHT, WHT, BLK, WHT, WHT, WHT, WHT},
                    {WHT, WHT, WHT, BLK, WHT, WHT, WHT, WHT},
                    {WHT, WHT, WHT, BLK, WHT, WHT, WHT, WHT},
                    {WHT, WHT, WHT, BLK, WHT, WHT, WHT, WHT},
                    {WHT, WHT, WHT, BLK, WHT, WHT, WHT, WHT},
                    {WHT, WHT, WHT, WHT, WHT, WHT, WHT, WHT}},
                   {{BLK, WHT, WHT, WHT, WHT, WHT, BLK, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, BLK, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, BLK, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, BLK, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, BLK, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, BLK, WHT},
                    {WHT, BLK, BLK, BLK, BLK, BLK, WHT, WHT},
                    {WHT, WHT, WHT, WHT, WHT, WHT, WHT, WHT}},
                   {{BLK, WHT, WHT, WHT, WHT, WHT, BLK, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, BLK, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, BLK, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, BLK, WHT},
                    {WHT, BLK, WHT, WHT, WHT, BLK, WHT, WHT},
                    {WHT, WHT, BLK, WHT, BLK, WHT, WHT, WHT},
                    {WHT, WHT, WHT, BLK, WHT, WHT, WHT, WHT},
                    {WHT, WHT, WHT, WHT, WHT, WHT, WHT, WHT}},
                   {{BLK, WHT, WHT, WHT, WHT, WHT, BLK, WHT},
                    {BLK, WHT, WHT, BLK, WHT, WHT, BLK, WHT},
                    {BLK, WHT, WHT, BLK, WHT, WHT, BLK, WHT},
                    {BLK, WHT, WHT, BLK, WHT, WHT, BLK, WHT},
                    {BLK, WHT, WHT, BLK, WHT, WHT, BLK, WHT},
                    {BLK, WHT, WHT, BLK, WHT, WHT, BLK, WHT},
                    {WHT, BLK, BLK, WHT, BLK, BLK, WHT, WHT},
                    {WHT, WHT, WHT, WHT, WHT, WHT, WHT, WHT}},
                   {{BLK, WHT, WHT, WHT, WHT, WHT, BLK, WHT},
                    {WHT, BLK, WHT, WHT, WHT, BLK, WHT, WHT},
                    {WHT, WHT, BLK, WHT, BLK, WHT, WHT, WHT},
                    {WHT, WHT, WHT, BLK, WHT, WHT, WHT, WHT},
                    {WHT, WHT, BLK, WHT, BLK, WHT, WHT, WHT},
                    {WHT, BLK, WHT, WHT, WHT, BLK, WHT, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, BLK, WHT},
                    {WHT, WHT, WHT, WHT, WHT, WHT, WHT, WHT}},
                   {{BLK, WHT, WHT, WHT, WHT, WHT, BLK, WHT},
                    {WHT, BLK, WHT, WHT, WHT, BLK, WHT, WHT},
                    {WHT, WHT, BLK, WHT, BLK, WHT, WHT, WHT},
                    {WHT, WHT, WHT, BLK, WHT, WHT, WHT, WHT},
                    {WHT, WHT, WHT, BLK, WHT, WHT, WHT, WHT},
                    {WHT, WHT, WHT, BLK, WHT, WHT, WHT, WHT},
                    {WHT, WHT, WHT, BLK, WHT, WHT, WHT, WHT},
                    {WHT, WHT, WHT, WHT, WHT, WHT, WHT, WHT}},
                   {{BLK, BLK, BLK, BLK, BLK, BLK, BLK, WHT},
                    {WHT, WHT, WHT, WHT, WHT, BLK, WHT, WHT},
                    {WHT, WHT, WHT, WHT, BLK, WHT, WHT, WHT},
                    {WHT, WHT, WHT, BLK, WHT, WHT, WHT, WHT},
                    {WHT, WHT, BLK, WHT, WHT, WHT, WHT, WHT},
                    {WHT, BLK, WHT, WHT, WHT, WHT, WHT, WHT},
                    {BLK, BLK, BLK, BLK, BLK, BLK, BLK, WHT},
                    {WHT, WHT, WHT, WHT, WHT, WHT, WHT, WHT}}};
   short number[10][CHARDIMEN][CHARDIMEN] = 
                  {{{WHT, WHT, BLK, BLK, BLK, WHT, WHT, WHT},
                    {WHT, BLK, WHT, WHT, WHT, BLK, WHT, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, BLK, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, BLK, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, BLK, WHT},
                    {WHT, BLK, WHT, WHT, WHT, BLK, WHT, WHT},
                    {WHT, WHT, BLK, BLK, BLK, WHT, WHT, WHT},
                    {WHT, WHT, WHT, WHT, WHT, WHT, WHT, WHT}},
                   {{WHT, WHT, WHT, BLK, WHT, WHT, WHT, WHT},
                    {WHT, WHT, BLK, BLK, WHT, WHT, WHT, WHT},
                    {WHT, BLK, WHT, BLK, WHT, WHT, WHT, WHT},
                    {WHT, WHT, WHT, BLK, WHT, WHT, WHT, WHT},
                    {WHT, WHT, WHT, BLK, WHT, WHT, WHT, WHT},
                    {WHT, WHT, WHT, BLK, WHT, WHT, WHT, WHT},
                    {WHT, BLK, BLK, BLK, BLK, BLK, WHT, WHT},
                    {WHT, WHT, WHT, WHT, WHT, WHT, WHT, WHT}},
                   {{WHT, BLK, BLK, BLK, BLK, BLK, WHT, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, BLK, WHT},
                    {WHT, WHT, WHT, WHT, WHT, WHT, BLK, WHT},
                    {WHT, BLK, BLK, BLK, BLK, BLK, WHT, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, WHT, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, WHT, WHT},
                    {BLK, BLK, BLK, BLK, BLK, BLK, BLK, WHT},
                    {WHT, WHT, WHT, WHT, WHT, WHT, WHT, WHT}},
                   {{WHT, BLK, BLK, BLK, BLK, BLK, WHT, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, BLK, WHT},
                    {WHT, WHT, WHT, WHT, WHT, WHT, BLK, WHT},
                    {WHT, BLK, BLK, BLK, BLK, BLK, WHT, WHT},
                    {WHT, WHT, WHT, WHT, WHT, WHT, BLK, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, BLK, WHT},
                    {WHT, BLK, BLK, BLK, BLK, BLK, WHT, WHT},
                    {WHT, WHT, WHT, WHT, WHT, WHT, WHT, WHT}},
                   {{BLK, WHT, WHT, WHT, WHT, WHT, WHT, WHT},
                    {BLK, WHT, WHT, WHT, WHT, BLK, WHT, WHT},
                    {BLK, WHT, WHT, WHT, WHT, BLK, WHT, WHT},
                    {BLK, WHT, WHT, WHT, WHT, BLK, WHT, WHT},
                    {BLK, BLK, BLK, BLK, BLK, BLK, BLK, WHT},
                    {WHT, WHT, WHT, WHT, WHT, BLK, WHT, WHT},
                    {WHT, WHT, WHT, WHT, WHT, BLK, WHT, WHT},
                    {WHT, WHT, WHT, WHT, WHT, WHT, WHT, WHT}},
                   {{BLK, BLK, BLK, BLK, BLK, BLK, BLK, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, WHT, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, WHT, WHT},
                    {BLK, BLK, BLK, BLK, BLK, BLK, WHT, WHT},
                    {WHT, WHT, WHT, WHT, WHT, WHT, BLK, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, BLK, WHT},
                    {WHT, BLK, BLK, BLK, BLK, BLK, WHT, WHT},
                    {WHT, WHT, WHT, WHT, WHT, WHT, WHT, WHT}},
                   {{WHT, BLK, BLK, BLK, BLK, BLK, WHT, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, BLK, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, WHT, WHT},
                    {BLK, BLK, BLK, BLK, BLK, BLK, WHT, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, BLK, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, BLK, WHT},
                    {WHT, BLK, BLK, BLK, BLK, BLK, WHT, WHT},
                    {WHT, WHT, WHT, WHT, WHT, WHT, WHT, WHT}},
                   {{BLK, BLK, BLK, BLK, BLK, BLK, BLK, WHT},
                    {BLK, WHT, WHT, WHT, WHT, BLK, WHT, WHT},
                    {WHT, WHT, WHT, WHT, BLK, WHT, WHT, WHT},
                    {WHT, WHT, WHT, BLK, WHT, WHT, WHT, WHT},
                    {WHT, WHT, BLK, WHT, WHT, WHT, WHT, WHT},
                    {WHT, WHT, BLK, WHT, WHT, WHT, WHT, WHT},
                    {WHT, WHT, BLK, WHT, WHT, WHT, WHT, WHT},
                    {WHT, WHT, WHT, WHT, WHT, WHT, WHT, WHT}},
                   {{WHT, BLK, BLK, BLK, BLK, BLK, WHT, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, BLK, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, BLK, WHT},
                    {WHT, BLK, BLK, BLK, BLK, BLK, WHT, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, BLK, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, BLK, WHT},
                    {WHT, BLK, BLK, BLK, BLK, BLK, WHT, WHT},
                    {WHT, WHT, WHT, WHT, WHT, WHT, WHT, WHT}},
                   {{WHT, BLK, BLK, BLK, BLK, BLK, WHT, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, BLK, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, BLK, WHT},
                    {WHT, BLK, BLK, BLK, BLK, BLK, BLK, WHT},
                    {WHT, WHT, WHT, WHT, WHT, WHT, BLK, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, BLK, WHT},
                    {WHT, BLK, BLK, BLK, BLK, BLK, WHT, WHT},
                    {WHT, WHT, WHT, WHT, WHT, WHT, WHT, WHT}}};
   short space[CHARDIMEN][CHARDIMEN] = 
                   {{WHT, WHT, WHT, WHT, WHT, WHT, WHT, WHT},
                    {WHT, WHT, WHT, WHT, WHT, WHT, WHT, WHT},
                    {WHT, WHT, WHT, WHT, WHT, WHT, WHT, WHT},
                    {WHT, WHT, WHT, WHT, WHT, WHT, WHT, WHT},
                    {WHT, WHT, WHT, WHT, WHT, WHT, WHT, WHT},
                    {WHT, WHT, WHT, WHT, WHT, WHT, WHT, WHT},
                    {WHT, WHT, WHT, WHT, WHT, WHT, WHT, WHT},
                    {WHT, WHT, WHT, WHT, WHT, WHT, WHT, WHT}};
   short block[CHARDIMEN][CHARDIMEN] = 
                   {{BLK, BLK, BLK, BLK, BLK, BLK, BLK, BLK},
                    {BLK, BLK, BLK, BLK, BLK, BLK, BLK, BLK},
                    {BLK, BLK, BLK, BLK, BLK, BLK, BLK, BLK},
                    {BLK, BLK, BLK, BLK, BLK, BLK, BLK, BLK},
                    {BLK, BLK, BLK, BLK, BLK, BLK, BLK, BLK},
                    {BLK, BLK, BLK, BLK, BLK, BLK, BLK, BLK},
                    {BLK, BLK, BLK, BLK, BLK, BLK, BLK, BLK},
                    {BLK, BLK, BLK, BLK, BLK, BLK, BLK, BLK}};
   short plus[CHARDIMEN][CHARDIMEN] = 
                   {{WHT, WHT, WHT, WHT, WHT, WHT, WHT, WHT},
                    {WHT, WHT, WHT, BLK, WHT, WHT, WHT, WHT},
                    {WHT, WHT, WHT, BLK, WHT, WHT, WHT, WHT},
                    {WHT, BLK, BLK, BLK, BLK, BLK, WHT, WHT},
                    {WHT, WHT, WHT, BLK, WHT, WHT, WHT, WHT},
                    {WHT, WHT, WHT, BLK, WHT, WHT, WHT, WHT},
                    {WHT, WHT, WHT, WHT, WHT, WHT, WHT, WHT},
                    {WHT, WHT, WHT, WHT, WHT, WHT, WHT, WHT}};
   short minus[CHARDIMEN][CHARDIMEN] = 
                   {{WHT, WHT, WHT, WHT, WHT, WHT, WHT, WHT},
                    {WHT, WHT, WHT, WHT, WHT, WHT, WHT, WHT},
                    {WHT, WHT, WHT, WHT, WHT, WHT, WHT, WHT},
                    {WHT, BLK, BLK, BLK, BLK, BLK, WHT, WHT},
                    {WHT, WHT, WHT, WHT, WHT, WHT, WHT, WHT},
                    {WHT, WHT, WHT, WHT, WHT, WHT, WHT, WHT},
                    {WHT, WHT, WHT, WHT, WHT, WHT, WHT, WHT},
                    {WHT, WHT, WHT, WHT, WHT, WHT, WHT, WHT}};
   short times[CHARDIMEN][CHARDIMEN] = 
                   {{WHT, WHT, WHT, WHT, WHT, WHT, WHT, WHT},
                    {WHT, BLK, WHT, WHT, WHT, BLK, WHT, WHT},
                    {WHT, WHT, BLK, WHT, BLK, WHT, WHT, WHT},
                    {BLK, BLK, BLK, BLK, BLK, BLK, BLK, WHT},
                    {WHT, WHT, BLK, WHT, BLK, WHT, WHT, WHT},
                    {WHT, BLK, WHT, WHT, WHT, BLK, WHT, WHT},
                    {WHT, WHT, WHT, WHT, WHT, WHT, WHT, WHT},
                    {WHT, WHT, WHT, WHT, WHT, WHT, WHT, WHT}};
   short div[CHARDIMEN][CHARDIMEN] = 
                   {{WHT, WHT, WHT, WHT, WHT, WHT, BLK, WHT},
                    {WHT, WHT, WHT, WHT, WHT, BLK, WHT, WHT},
                    {WHT, WHT, WHT, WHT, BLK, WHT, WHT, WHT},
                    {WHT, WHT, WHT, BLK, WHT, WHT, WHT, WHT},
                    {WHT, WHT, BLK, WHT, WHT, WHT, WHT, WHT},
                    {WHT, BLK, WHT, WHT, WHT, WHT, WHT, WHT},
                    {BLK, WHT, WHT, WHT, WHT, WHT, WHT, WHT},
                    {WHT, WHT, WHT, WHT, WHT, WHT, WHT, WHT}};
   short equals[CHARDIMEN][CHARDIMEN] = 
                   {{WHT, WHT, WHT, WHT, WHT, WHT, WHT, WHT},
                    {WHT, WHT, WHT, WHT, WHT, WHT, WHT, WHT},
                    {WHT, BLK, BLK, BLK, BLK, BLK, WHT, WHT},
                    {WHT, WHT, WHT, WHT, WHT, WHT, WHT, WHT},
                    {WHT, BLK, BLK, BLK, BLK, BLK, WHT, WHT},
                    {WHT, WHT, WHT, WHT, WHT, WHT, WHT, WHT},
                    {WHT, WHT, WHT, WHT, WHT, WHT, WHT, WHT},
                    {WHT, WHT, WHT, WHT, WHT, WHT, WHT, WHT}};
   short openbr[CHARDIMEN][CHARDIMEN] = 
                   {{WHT, WHT, WHT, BLK, BLK, WHT, WHT, WHT},
                    {WHT, WHT, BLK, WHT, WHT, WHT, WHT, WHT},
                    {WHT, BLK, WHT, WHT, WHT, WHT, WHT, WHT},
                    {WHT, BLK, WHT, WHT, WHT, WHT, WHT, WHT},
                    {WHT, BLK, WHT, WHT, WHT, WHT, WHT, WHT},
                    {WHT, WHT, BLK, WHT, WHT, WHT, WHT, WHT},
                    {WHT, WHT, WHT, BLK, BLK, WHT, WHT, WHT},
                    {WHT, WHT, WHT, WHT, WHT, WHT, WHT, WHT}};
   short closebr[CHARDIMEN][CHARDIMEN] = 
                   {{WHT, WHT, BLK, BLK, WHT, WHT, WHT, WHT},
                    {WHT, WHT, WHT, WHT, BLK, WHT, WHT, WHT},
                    {WHT, WHT, WHT, WHT, WHT, BLK, WHT, WHT},
                    {WHT, WHT, WHT, WHT, WHT, BLK, WHT, WHT},
                    {WHT, WHT, WHT, WHT, WHT, BLK, WHT, WHT},
                    {WHT, WHT, WHT, WHT, BLK, WHT, WHT, WHT},
                    {WHT, WHT, BLK, BLK, WHT, WHT, WHT, WHT},
                    {WHT, WHT, WHT, WHT, WHT, WHT, WHT, WHT}};
   short lessthan[CHARDIMEN][CHARDIMEN] = 
                   {{WHT, WHT, WHT, WHT, WHT, BLK, WHT, WHT},
                    {WHT, WHT, WHT, WHT, BLK, WHT, WHT, WHT},
                    {WHT, WHT, WHT, BLK, WHT, WHT, WHT, WHT},
                    {WHT, WHT, BLK, WHT, WHT, WHT, WHT, WHT},
                    {WHT, WHT, WHT, BLK, WHT, WHT, WHT, WHT},
                    {WHT, WHT, WHT, WHT, BLK, WHT, WHT, WHT},
                    {WHT, WHT, WHT, WHT, WHT, BLK, WHT, WHT},
                    {WHT, WHT, WHT, WHT, WHT, WHT, WHT, WHT}};
   short gtrthan[CHARDIMEN][CHARDIMEN] = 
                   {{WHT, WHT, BLK, WHT, WHT, WHT, WHT, WHT},
                    {WHT, WHT, WHT, BLK, WHT, WHT, WHT, WHT},
                    {WHT, WHT, WHT, WHT, BLK, WHT, WHT, WHT},
                    {WHT, WHT, WHT, WHT, WHT, BLK, WHT, WHT},
                    {WHT, WHT, WHT, WHT, BLK, WHT, WHT, WHT},
                    {WHT, WHT, WHT, BLK, WHT, WHT, WHT, WHT},
                    {WHT, WHT, BLK, WHT, WHT, WHT, WHT, WHT},
                    {WHT, WHT, WHT, WHT, WHT, WHT, WHT, WHT}};
   short period[CHARDIMEN][CHARDIMEN] = 
                   {{WHT, WHT, WHT, WHT, WHT, WHT, WHT, WHT},
                    {WHT, WHT, WHT, WHT, WHT, WHT, WHT, WHT},
                    {WHT, WHT, WHT, WHT, WHT, WHT, WHT, WHT},
                    {WHT, WHT, WHT, WHT, WHT, WHT, WHT, WHT},
                    {WHT, WHT, BLK, BLK, BLK, WHT, WHT, WHT},
                    {WHT, WHT, BLK, BLK, BLK, WHT, WHT, WHT},
                    {WHT, WHT, BLK, BLK, BLK, WHT, WHT, WHT},
                    {WHT, WHT, WHT, WHT, WHT, WHT, WHT, WHT}};
   short comma[CHARDIMEN][CHARDIMEN] = 
                   {{WHT, WHT, WHT, WHT, WHT, WHT, WHT, WHT},
                    {WHT, WHT, WHT, WHT, WHT, WHT, WHT, WHT},
                    {WHT, WHT, WHT, WHT, WHT, WHT, WHT, WHT},
                    {WHT, WHT, BLK, BLK, BLK, WHT, WHT, WHT},
                    {WHT, WHT, BLK, BLK, BLK, WHT, WHT, WHT},
                    {WHT, WHT, WHT, BLK, WHT, WHT, WHT, WHT},
                    {WHT, WHT, BLK, WHT, WHT, WHT, WHT, WHT},
                    {WHT, WHT, WHT, WHT, WHT, WHT, WHT, WHT}};
   short semicln[CHARDIMEN][CHARDIMEN] = 
                   {{WHT, WHT, BLK, BLK, BLK, WHT, WHT, WHT},
                    {WHT, WHT, BLK, BLK, BLK, WHT, WHT, WHT},
                    {WHT, WHT, WHT, WHT, WHT, WHT, WHT, WHT},
                    {WHT, WHT, BLK, BLK, BLK, WHT, WHT, WHT},
                    {WHT, WHT, BLK, BLK, BLK, WHT, WHT, WHT},
                    {WHT, WHT, WHT, BLK, WHT, WHT, WHT, WHT},
                    {WHT, WHT, BLK, WHT, WHT, WHT, WHT, WHT},
                    {WHT, WHT, WHT, WHT, WHT, WHT, WHT, WHT}};
   short colon[CHARDIMEN][CHARDIMEN] = 
                   {{WHT, WHT, WHT, BLK, WHT, WHT, WHT, WHT},
                    {WHT, WHT, BLK, BLK, BLK, WHT, WHT, WHT},
                    {WHT, WHT, WHT, BLK, WHT, WHT, WHT, WHT},
                    {WHT, WHT, WHT, WHT, WHT, WHT, WHT, WHT},
                    {WHT, WHT, WHT, BLK, WHT, WHT, WHT, WHT},
                    {WHT, WHT, BLK, BLK, BLK, WHT, WHT, WHT},
                    {WHT, WHT, WHT, BLK, WHT, WHT, WHT, WHT},
                    {WHT, WHT, WHT, WHT, WHT, WHT, WHT, WHT}};

   /* ========================================================================
      put selected letter in cm array 
      ========================================================================
   */

   if ( (c>='a' && c<='z') ) {  /* lower case letter (still draw upper case) */
      idx = c-'a';
      for (i=0;i<=7;i++) for (j=0;j<=7;j++) cm[i][j] = capletter[idx][i][j];
   } else if ( (c>='A' && c<='Z') ) { /* upper case letter */
      idx = c-'A'; 
      for (i=0;i<=7;i++) for (j=0;j<=7;j++) cm[i][j] = capletter[idx][i][j];
   } else if ( (c>='0' && c<='9') ) { /* number */
      idx = c-'0';
      for (i=0;i<=7;i++) for (j=0;j<=7;j++) cm[i][j] = number[idx][i][j];
   } else {
      switch (c) {
       case '_': for (i=0;i<=7;i++) for (j=0;j<=7;j++) cm[i][j]=space[i][j];
                 break;
       case ' ': for (i=0;i<=7;i++) for (j=0;j<=7;j++) cm[i][j]=space[i][j];
                 break;
       case '+': for (i=0;i<=7;i++) for (j=0;j<=7;j++) cm[i][j]=plus[i][j];
                 break;
       case '-': for (i=0;i<=7;i++) for (j=0;j<=7;j++) cm[i][j]=minus[i][j];
                 break;
       case '*': for (i=0;i<=7;i++) for (j=0;j<=7;j++) cm[i][j]=times[i][j];
                 break;
       case '/': for (i=0;i<=7;i++) for (j=0;j<=7;j++) cm[i][j]=div[i][j];
                 break;
       case '=': for (i=0;i<=7;i++) for (j=0;j<=7;j++) cm[i][j]=equals[i][j];
                 break;
       case '(': for (i=0;i<=7;i++) for (j=0;j<=7;j++) cm[i][j]=openbr[i][j];
                 break;
       case ')': for (i=0;i<=7;i++) for (j=0;j<=7;j++) cm[i][j]=closebr[i][j];
                 break;
       case '<': for (i=0;i<=7;i++) for (j=0;j<=7;j++) cm[i][j]=lessthan[i][j];
                 break;
       case '>': for (i=0;i<=7;i++) for (j=0;j<=7;j++) cm[i][j]=gtrthan[i][j];
                 break;
       case '.': for (i=0;i<=7;i++) for (j=0;j<=7;j++) cm[i][j]=period[i][j];
                 break;
       case ',': for (i=0;i<=7;i++) for (j=0;j<=7;j++) cm[i][j]=comma[i][j];
                 break;
       case ';': for (i=0;i<=7;i++) for (j=0;j<=7;j++) cm[i][j]=semicln[i][j];
                 break;
       case ':': for (i=0;i<=7;i++) for (j=0;j<=7;j++) cm[i][j]=colon[i][j];
                 break;
       default:  for (i=0;i<=7;i++) for (j=0;j<=7;j++) cm[i][j]=block[i][j];
                 break;
      }
   }
}

/* ===========================================================================
   DEFCOLMAP: MAY 03

   Define RGB colour map from predefined data.
   NOTE: This routine assumes COLMAPSIZE=256
   ===========================================================================
*/
void defcolmap(icmval,colmap)
   int colmap;
   double *icmval;
{
   short i,j;
   short blank[256] = 
        {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 };
   short gry0[256] = 
        {  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 
          16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 
          32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 
          48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 
          64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 
          80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 
          96, 97, 98, 99,100,101,102,103,104,105,106,107,108,109,110,111,
         112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,
         128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,
         144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,
         160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,
         176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,
         192,193,194,195,196,197,198,199,200,201,202,203,204,205,206,207,
         208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,
         224,225,226,227,228,229,230,231,232,233,234,235,236,237,238,239,
         240,241,242,243,244,245,246,247,248,249,250,251,252,253,254,255 };
   short gry1[256] = 
        {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 48, 96,144,  
         192,152,112, 72, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 
          32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 80,128,176,  
         224,184,144,104, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 
          64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64,112,160,208,  
         255,216,176,136, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 
          96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 
          96,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,  
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,  
         160,160,160,160,160,160,160,160,160,160,160,160,160,160,160,160,
         160,160,160,160,160,160,160,160,160,160,160,160,160,120, 80, 40,
           0, 48, 96,144,192,192,192,192,192,192,192,192,192,192,192,192,  
         192,192,192,192,192,192,192,192,192,192,192,192,192,152,112, 72,
          32, 80,128,176,224,224,224,224,224,224,224,224,224,224,224,224,
         224,224,224,224,224,224,224,224,224,224,224,224,224,184,144,104,
          64,112,160,208,255,255,255,255,255,255,255,255,255,255,255,255 };
   short gry2[256] = 
        {160,160,160,160,160,160,160,160,160,160,160,160,160,160,160,160,
         160,160,160,160,160,160,160,160,160,160,160,160,160,160,160,160,
         160,160,160,160,160,160,160,160,160,160,160,160,160,160,160,160,
         160,160,160,160,160,160,160,160,160,160,160,160,160,160,160,160,
         160,160,160,160,160,160,160,160,160,160,160,160,160,160,160,160,
         160,160,160,160,160,160,160,160,160,160,160,160,160,160,160,160,
         160,160,160,160,160,160,160,160,160,160,160,160,160,160,160,160,
         160,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         224,224,224,224,224,224,224,224,224,224,224,224,224,224,224,224,
         224,224,224,224,224,224,224,224,224,224,224,224,224,224,224,224,
         224,224,224,224,224,224,224,224,224,224,224,224,224,224,224,224,
         224,224,224,224,224,224,224,224,224,224,224,224,224,224,224,224,
         224,224,224,224,224,224,224,224,224,224,224,224,224,224,224,224,
         224,224,224,224,224,224,224,224,224,224,224,224,224,224,224,224,
         224,224,224,224,224,224,224,224,224,224,224,224,224,224,224,224 };
   short gry3[256] = 
        {  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 
          16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 
          32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 
          48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 
          64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 
          80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 
          96, 97, 98, 99,100,101,102,103,104,105,106,107,108,109,110,111,
         112,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,
         160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,
         176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,
         192,193,194,195,196,197,198,199,200,201,202,203,204,205,206,207,
         208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,
         224,225,226,227,228,229,230,231,232,233,234,235,236,237,238,239,
         240,241,242,243,244,245,246,247,248,249,250,251,252,253,254,255 };
   short gry4[256] = 
        {255,254,253,252,251,250,249,248,247,246,245,244,243,242,241,240,
         239,238,237,236,235,234,233,232,231,230,229,228,227,226,225,224,
         223,222,221,220,219,218,217,216,215,214,213,212,211,210,209,208,
         207,206,205,204,203,202,201,200,199,198,197,196,195,194,193,192,
         191,190,189,188,187,186,185,184,183,182,181,180,179,178,177,176,
         175,174,173,172,171,170,169,168,167,166,165,164,163,162,161,160,
         159,158,157,156,155,154,153,152,151,150,149,148,147,146,145,144,
         143,142,141,140,139,138,137,136,135,134,133,132,131,130,129,128,
         127,126,125,124,123,122,121,120,119,118,117,116,115,114,113,112,
         111,110,109,108,107,106,105,104,103,102,101,100, 99, 98, 97, 96,
         95, 94, 93, 92, 91, 90, 89, 88, 87, 86, 85, 84, 83, 82, 81,  80,
         79, 78, 77, 76, 75, 74, 73,  72,71, 70, 69, 68, 67, 66, 65,  64,
         63, 62, 61, 60, 59, 58, 57, 56, 55, 54, 53, 52, 51, 50, 49,  48,
         47, 46, 45, 44, 43, 42, 41, 40, 39, 38, 37, 36, 35, 34, 33,  32,
         31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17,  16,
         15, 14, 13, 12, 11, 10, 9,  8,  7,  6,  5,  4,  3,  2,  1,   0  };
   short gry5[256] = 
        {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 48, 96,144,
         192,152,112, 72, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32,
          32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 80,128,176,
         224,184,144,104, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64,
          64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64,112,160,208,
         255,216,176,136, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96,
          96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96,136,176,216,
         255,224,192,160,128,128,128,128,128,128,128,128,128,128,128,128,
         128,128,128,128,128,128,128,128,128,128,128,128,128, 96, 64, 32,
           0, 40, 80,120,160,160,160,160,160,160,160,160,160,160,160,160,
         160,160,160,160,160,160,160,160,160,160,160,160,160,120, 80, 40, 
           0, 48, 96,144,192,192,192,192,192,192,192,192,192,192,192,192,
         192,192,192,192,192,192,192,192,192,192,192,192,192,152,112, 72,
          32, 80,128,176,224,224,224,224,224,224,224,224,224,224,224,224,
         224,224,224,224,224,224,224,224,224,224,224,224,224,184,144,104,
          64,112,160,208,255,255,255,255,255,255,255,255,255,255,255,255};
   short gry6[256] = 
        {255,255,255,255,255,255,255,255,255,255,255,255,255,208,160,112,
          64,104,144,184,224,224,224,224,224,224,224,224,224,224,224,224,
         224,224,224,224,224,224,224,224,224,224,224,224,224,176,128, 80,
          32, 72,112,152,192,192,192,192,192,192,192,192,192,192,192,192,
         192,192,192,192,192,192,192,192,192,192,192,192,192,144, 96, 48,
           0, 40, 80,120,160,160,160,160,160,160,160,160,160,160,160,160,
         160,160,160,160,160,160,160,160,160,160,160,160,160,120, 80, 40, 
           0, 32, 64, 96,128,128,128,128,128,128,128,128,128,128,128,128,
         128,128,128,128,128,128,128,128,128,128,128,128,128,160,192,224,
         255,216,176,136, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96,
          96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96,136,176,216,
         255,208,160,112, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64,
          64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64,104,144,184,
         224,176,128, 80, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32,
          32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 72,112,152,
         192,144, 96, 48,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 };
   short gry7[256] = 
        { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64,
          64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 80, 96,
         112,128,128,128,128,128,128,128,128,128,128,128,128,128,128,128,
         128,128,128,128,128,128,128,128,128,128,128,128,128,128,128,128,
         130,132,134,136,138,140,142,144,146,148,150,152,154,156,158,160,
         162,164,166,168,170,172,174,176,178,180,182,184,186,188,190,192,
         194,196,198,200,202,204,206,208,210,212,214,216,218,220,222,224,
         226,228,230,232,234,236,238,240,242,244,246,248,250,252,254,255,
         255,255,255,255,254,254,253,253,252,252,251,251,250,250,249,249,
         248,248,247,247,246,246,245,245,244,244,243,243,242,242,241,241,
         240,240,239,239,238,238,237,237,236,236,235,235,234,234,233,233,
         232,232,231,231,230,230,229,229,228,228,227,227,226,226,225,225,
         224,224,224,224,224,224,224,224,224,224,224,224,224,224,224,224,
         224,224,224,224,224,224,224,224,224,224,224,224,224,224,228,232,
         236,240,240,240,240,240,240,240,240,240,240,240,240,240,240,240,
         240,240,240,240,240,240,240,240,240,240,240,240,240,240,240,240};
   short gry8[256] = 
        {  0, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96,
          96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96,
          96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96,
          96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96,
          96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96,
          96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96,
          96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96,
          96,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         160,160,160,160,160,160,160,160,160,160,160,160,160,160,160,160,
         160,160,160,160,160,160,160,160,160,160,160,160,160,160,160,160,
         160,160,160,160,160,160,160,160,160,160,160,160,160,160,160,160,
         160,160,160,160,160,160,160,160,160,160,160,160,160,160,160,160,
         160,160,160,160,160,160,160,160,160,160,160,160,160,160,160,160,
         160,160,160,160,160,160,160,160,160,160,160,160,160,160,160,160,
         160,160,160,160,160,160,160,160,160,160,160,160,160,160,160,255};
   short gry9[256] = 
        {128,128,128,128,128,128,128,128,128,128,128,128,128,128,128,128,
         128,128,128,128,128,128,128,128,128,128,128,128,128,128,128,128,
         128,128,128,128,128,128,128,128,128,128,128,128,128,128,128,128,
         128,128,128,128,128,128,128,128,128,128,128,128,128,128,128,128,
         128,128,128,128,128,128,128,128,128,128,128,128,128,128,128,128,
         128,128,128,128,128,128,128,128,128,128,128,128,128,128,128,128,
         128,128,128,128,128,128,128,128,128,128,128,128,128,128,128,128,
         128,128,128,128,128,128,128,128,128,128,128,128,128,128,128,128,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255};
   short gry10[256] = 
        {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0, 75, 75, 75, 75, 75, 75, 75, 75, 
          75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 
          75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75,108,
         108,108,108,108,108,108,108,108,108,108,108,108,108,108,108,108,
         108,108,108,108,108,108,108,108,108,108,108,108,108,108,108,108,
         108,108,108,108,108,108,122,136,150,164,178,192,206,220,234,248,
         255,238,228,218,208,198,188,178,168,158,148,148,148,148,148,148,
         148,148,148,148,148,148,148,148,148,148,148,148,148,148,148,148,
         148,148,148,148,148,148,148,148,148,148,148,148,148,148,148,148,
         148,180,180,180,180,180,180,180,180,180,180,180,180,180,180,180,
         180,180,180,180,180,180,180,180,180,180,180,180,180,180,180,180,
         180,180,180,180,180,180,180,180,224,224,224,224,224,224,224,224,
         224,224,224,224,224,224,224,224,224,224,224,224,224,224,224,224,
         224,224,224,224,224,224,224,224,224,224,224,224,224,224,224,255 };
   short gry11[256] = 
        {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 
          36, 36, 36, 36, 36, 36, 36, 36, 75, 75, 75, 75, 75, 75, 75, 75, 
          75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 90, 90, 90, 90, 90, 
          90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90,108,
         108,108,108,108,108,108,108,108,108,108,108,108,108,108,108,108,
         108,108,124,124,124,124,124,124,124,124,124,124,124,124,124,124,
         124,124,124,124,124,124,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,130,130,130,130,130,130,
         130,130,130,130,130,130,130,130,130,130,130,130,130,130,148,148,
         148,148,148,148,148,148,148,148,148,148,148,148,148,148,148,148,
         148,165,165,165,165,165,165,165,165,165,165,165,165,165,165,165,
         165,165,165,165,165,180,180,180,180,180,180,180,180,180,180,180,
         180,180,180,180,180,180,180,180,200,200,200,200,200,200,200,200,
         200,200,200,200,200,200,200,200,200,200,200,200,224,224,224,224,
         224,224,224,224,224,224,224,224,224,224,224,224,224,224,224,255 };
   short gry12[256] = 
        {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255 };
   short gry13[256] = 
        {255,128,  0,  0, 64, 64, 64, 64,112,112,112,112,112,112,112,112,
         160,160,160,160,160,160,160,160,160,160,160,160,160,160,160,160,
         192,192,192,192,192,192,192,192,192,192,192,192,192,192,192,192,
         192,192,192,192,192,192,192,192,192,192,192,192,192,192,192,192,
         240,240,240,240,240,240,240,240,240,240,240,240,240,240,240,240,
         240,240,240,240,240,240,240,240,240,240,240,240,240,240,240,240,
         240,240,240,240,240,240,240,240,240,240,240,240,240,240,240,240,
         240,240,240,240,240,240,240,240,240,240,240,240,240,240,240,240,
         224,224,224,224,224,224,224,224,224,224,224,224,224,224,224,224,
         224,224,224,224,224,224,224,224,224,224,224,224,224,224,224,224,
         224,224,224,224,224,224,224,224,224,224,224,224,224,224,224,224,
         224,224,224,224,224,224,224,224,224,224,224,224,224,224,224,224,
         224,224,224,224,224,224,224,224,224,224,224,224,224,224,224,224,
         224,224,224,224,224,224,224,224,224,224,224,224,224,224,224,224,
         224,224,224,224,224,224,224,224,224,224,224,224,224,224,224,224,
         224,224,224,224,224,224,224,224,224,224,224,224,224,224,224,224 };
   short gry14[256] = 
        {255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         236,216,196,176,156,136,116, 96, 96, 96, 96, 96, 96, 96, 96, 96,
          96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96,
          96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96,
          96, 96, 96,104,112,120,128,136,144,152,160,160,160,160,160,160,
         160,160,160,160,160,160,160,160,160,160,160,160,160,160,160,160,
         160,160,160,160,160,160,160,160,160,160,160,160,160,160,160,160,
         160,160,160,160,160,160,168,176,184,192,200,208,216,224,224,224,
         224,224,224,224,224,224,224,224,224,224,224,224,224,224,224,224,
         224,224,224,224,224,224,224,224,224,224,224,224,224,224,224,224,
         224,224,224,224,224,224,224,224,224,226,228,230,232,234,236,238,
         240,240,240,240,240,240,240,240,240,240,240,240,240,240,240,240,
         240,240,240,240,240,240,240,240,240,240,240,240,240,240,240,240,
         240,240,240,240,240,240,240,240,240,240,240,240,240,240,240,240 };
   short gry15[256] = 
        {255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         244,232,220,208,196,184,172,160,160,160,160,160,160,160,160,160,
         160,160,160,160,160,160,160,160,160,160,160,160,160,160,160,160,
         160,160,160,160,160,160,160,160,160,160,160,160,160,160,160,160,
         160,160,160,168,176,184,192,200,208,216,224,224,224,224,224,224,
         224,224,224,224,224,224,224,224,224,224,224,224,224,224,224,224,
         224,224,224,224,224,224,224,224,224,224,224,224,224,224,224,224,
         224,224,224,224,224,224,212,200,188,176,164,152,140,128,128,128,
         128,128,128,128,128,128,128,128,128,128,128,128,128,128,128,128,
         128,128,128,128,128,128,128,128,128,128,128,128,128,128,128,128,
         128,128,128,128,128,128,128,128,128,112, 96, 80, 64, 48, 32, 16,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 };
   short gry16[256] = 
        {255,253,251,249,247,245,243,241,239,237,235,233,231,229,227,225,
         223,221,219,217,215,213,211,209,207,205,203,201,199,197,195,193,
         191,189,187,185,183,181,179,177,175,173,171,169,167,165,163,161,
         159,157,155,153,151,149,147,145,143,141,139,137,135,133,131,129,
         127,125,123,121,119,117,115,113,111,109,107,105,103,101, 99, 97,
          95, 93, 91, 89, 87, 85, 83, 81, 79, 77, 75, 73, 71, 69, 67, 65,
          63, 61, 59, 57, 55, 53, 51, 49, 47, 45, 43, 41, 39, 37, 35, 33,  
          31, 29, 27, 25, 23, 21, 19, 17, 15, 13, 11,  9,  7,  5,  3,  1,
           0,  2,  4,  6,  8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 
          32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60, 62, 
          64, 66, 68, 70, 72, 74, 76, 78, 80, 82, 84, 86, 88, 90, 92, 94, 
          96, 98,100,102,104,106,108,110,112,114,116,118,120,122,124,126,
         128,130,132,134,136,138,140,142,144,146,148,150,152,154,156,158,
         160,162,164,166,168,170,172,174,176,178,180,182,184,186,188,190,
         192,194,196,198,200,202,204,206,208,210,212,214,216,218,220,222,
         224,226,228,230,232,234,236,238,240,242,244,246,248,250,252,254 };
   short gry17[256] = 
        {255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         192,128, 64,  0,24,  48, 72, 96, 96, 96, 96, 96, 96, 96, 96, 96,
          96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96,
          96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96,
          96, 96, 96, 72, 48,24,   0,40,  80,120,160,160,160,160,160,160,
         160,160,160,160,160,160,160,160,160,160,160,160,160,160,160,160,
         160,160,160,160,160,160,160,160,160,160,160,160,160,160,160,160,
         160,160,160,160,160,160,120, 80,40,   0,56, 112,178,224,224,224,
         224,224,224,224,224,224,224,224,224,224,224,224,224,224,224,224,
         224,224,224,224,224,224,224,224,224,224,224,224,224,224,224,224,
         224,224,224,224,224,224,224,224,224,178,112,56,   0, 64,128,192,
         240,240,240,240,240,240,240,240,240,240,240,240,240,240,240,240,
         240,240,240,240,240,240,240,240,240,240,240,240,240,240,240,240,
         240,240,240,240,240,240,240,240,240,240,240,240,240,240,240,240 };
   short gry18[256] = 
        {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
          40, 80,120,160,144,128,112, 96, 96, 96, 96, 96, 96, 96, 96, 96,
          96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96,
          96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96,
          96, 96, 96,136,176,216,255,232,208,184,160,160,160,160,160,160,
         160,160,160,160,160,160,160,160,160,160,160,160,160,160,160,160,
         160,160,160,160,160,160,160,160,160,160,160,160,160,160,160,160,
         160,160,160,160,160,160,120, 80,40,   0,56, 112,178,224,224,224,
         224,224,224,224,224,224,224,224,224,224,224,224,224,224,224,224,
         224,224,224,224,224,224,224,224,224,224,224,224,224,224,224,224,
         224,224,224,224,224,224,224,224,224,184,144,104, 64,112,160,208,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255 };
   short gry19[256] = 
        {255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         208,160,112, 64,104,144,184,224,224,224,224,224,224,224,224,224,
         224,224,224,224,224,224,224,224,224,224,224,224,224,224,224,224,
         224,224,224,224,224,224,224,224,224,224,224,224,224,224,224,224,
         224,224,224,184,144,104, 64, 88,112,136,160,160,160,160,160,160,
         160,160,160,160,160,160,160,160,160,160,160,160,160,160,160,160,
         160,160,160,160,160,160,160,160,160,160,160,160,160,160,160,160,
         160,160,160,160,160,160,184,208,232,255,216,176,136, 96, 96, 96,
          96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96,
          96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96,
          96, 96, 96, 96, 96, 96, 96, 96, 96,112,128,144,160,120, 80, 40,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 };
   short col0[3][256] = 
       {{  0,224,224,224,224,224,224,224,224,224,224,224,224,184,144,104,
          64,84, 104,124,144,144,144,144,144,144,144,144,144,144,144,144,
         144,144,144,144,144,144,144,144,144,144,144,144,144,124,104,84, 
          64, 48, 32, 16,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 16, 32, 48,
          64, 48, 32, 16,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,208,160,112,
          64,112,160,208,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,208,160,112,
          64,112,160,208,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,208,160,112,
          64,112,160,208,255,255,255,255,255,255,255,255,255,255,255,255 },
        {  0,184,184,184,184,184,184,184,184,184,184,184,184,154,124,94, 
          64,112,160,208,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,208,160,112,
          64,112,160,208,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,208,160,112,
          64,112,160,208,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         144,144,144,144,144,144,144,144,144,144,144,144,144,144,144,144,
         144,144,144,144,144,144,144,144,144,144,144,144,144,124,104,84, 
          64,66, 68, 70,  72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72,
          72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72,70, 68, 66, 
          64, 48, 32, 16,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 16, 32, 48,
          64, 48, 32, 16,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,255 },
        {  0,255,255,255,255,255,255,255,255,255,255,255,255,208,160,112,
          64,112,160,208,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,208,160,112,
          64,112,160,208,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,208,160,112,
          64, 48, 32, 16,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 16, 32, 48,
          64, 48, 32, 16,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 16, 32, 48,
          64, 48, 32, 16,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 16, 32, 48,
          64,112,160,208,255,255,255,255,255,255,255,255,255,255,255,255 }};
   short col1[3][256] = 
       {{224,224,224,224,224,224,224,224,224,224,224,224,224,184,144,104,
          64,84, 104,124,144,144,144,144,144,144,144,144,144,144,144,144,
         144,144,144,144,144,144,144,144,144,144,144,144,144,124,104,84, 
          64, 48, 32, 16,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 16, 32, 48,
          64, 48, 32, 16,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 16, 32, 48,
          64,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,208,160,112,
          64,112,160,208,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,208,160,112,
          64,112,160,208,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,208,160,112,
          64,112,160,208,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,208,160,112,
          64,112,160,208,255,255,255,255,255,255,255,255,255,255,255,255 },
        {184,184,184,184,184,184,184,184,184,184,184,184,184,154,124,94, 
          64,112,160,208,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,208,160,112,
          64,112,160,208,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,208,160,112,
          64,112,160,208,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,208,160,112,
          64,112,160,208,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,208,160,112,
          64,84, 104,124,144,144,144,144,144,144,144,144,144,144,144,144,
         144,144,144,144,144,144,144,144,144,144,144,144,144,124,104,84, 
          64,66, 68, 70,  72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72,
          72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72,70, 68, 66, 
          64, 48, 32, 16,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 16, 32, 48,
          64, 48, 32, 16,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,128 },
        {255,255,255,255,255,255,255,255,255,255,255,255,255,208,160,112,
          64,112,160,208,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,208,160,112,
          64,112,160,208,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,208,160,112,
          64, 48, 32, 16,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 16, 32, 48,
          64, 48, 32, 16,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 16, 32, 48,
          64, 48, 32, 16,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 16, 32, 48,
          64, 48, 32, 16,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 16, 32, 48,
          64, 48, 32, 16,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 16, 32, 48,
          64,112,160,208,255,255,255,255,255,255,255,255,255,255,255,255 }};
   short col2[3][256] = 
       {{  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255 },
        {255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         144,144,144,144,144,144,144,144,144,144,144,144,144,144,144,144,
         144,144,144,144,144,144,144,144,144,144,144,144,144,144,144,144,
         144,144,144,144,144,144,144,144,144,144,144,144,144,144,144,144,
         144,144,144,144,144,144,144,144,144,144,144,144,144,144,144,144,
         144,144,144,144,144,144,144,144,144,144,144,144,144,144,144,144,
         144,144,144,144,144,144,144,144,144,144,144,144,144,144,144,144,
         144,144,144,144,144,144,144,144,144,144,144,144,144,144,144,144 },
        {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,0  }};
   short col3[3][256] = 
       {{255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,192,
         128, 64,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 64,
         128,192,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255 },
        {255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,192,
         128, 64,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 64,
         128,192,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,228,
         200,172,144,144,144,144,144,144,144,144,144,144,144,144,144,144,
         144,144,144,144,144,144,144,144,144,144,144,144,144,144,144,126,
         108,90,  72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72,
          72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72,54, 
         36, 18,   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,128 },
        {255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,192,
         128, 64,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 64,
         128,192,255,255,255,255,255,255,255,255,255,255,255,255,255,255 }};
   short col4[3][256] = 
       {{255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
         128,128,128,128,128,128,128,128,128,128,128,128,128,128,128,128,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
         192,192,192,192,192,192,192,192,192,192,192,192,192,192,192,192,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
         192,192,192,192,192,192,192,192,192,192,192,192,192,192,192,192,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255 },
        {255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
         128,128,128,128,128,128,128,128,128,128,128,128,128,128,128,128,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         144,144,144,144,144,144,144,144,144,144,144,144,144,144,144,144,
         240,240,240,240,240,240,240,240,240,240,240,240,240,240,240,240,
          72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72,
         208,208,208,208,208,208,208,208,208,208,208,208,208,208,208,208, 
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
         192,192,192,192,192,192,192,192,192,192,192,192,192,192,192,192,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
        {255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
         192,192,192,192,192,192,192,192,192,192,192,192,192,192,192,192,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
         224,224,224,224,224,224,224,224,224,224,224,224,224,224,224,224,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
         224,224,224,224,224,224,224,224,224,224,224,224,224,224,224,224,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
         192,192,192,192,192,192,192,192,192,192,192,192,192,192,192,192,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
         192,192,192,192,192,192,192,192,192,192,192,192,192,192,192,192,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255 }};
   short col5[3][256] = 
       {{  0,224,224,224,224,224,224,224,224,224,224,224,224,184,144,104,
          64,84, 104,124,144,144,144,144,144,144,144,144,144,144,144,144,
         144,144,144,144,144,144,144,144,144,144,144,144,144,124,104,84, 
          64, 48, 32, 16,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 16, 32, 48,
          64, 48, 32, 16,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,208,160,112,
          64,112,160,208,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,208,160,112,
          64,112,160,208,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,208,160,112,
          64,112,160,208,255,255,255,255,255,255,255,255,255,255,255,255 },
        {  0,184,184,184,184,184,184,184,184,184,184,184,184,154,124, 94, 
          64,112,160,208,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,208,160,112,
          64,112,160,208,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,208,160,112,
          64,112,160,208,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         144,144,144,144,144,144,144,144,144,144,144,144,144,144,144,144,
         144,144,144,144,144,144,144,144,144,144,144,144,144,124,104, 84, 
          64, 66, 68, 70, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72,
          72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 70, 68, 66, 
          64, 48, 32, 16,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 16, 32, 48,
          64, 48, 32, 16,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,255 },
        {  0,255,255,255,255,255,255,255,255,255,255,255,255,208,160,112,
          64,112,160,208,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,208,160,112,
          64,112,160,208,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,208,160,112,
          64, 48, 32, 16,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,255,255,255,255,255,255,255,255,255,255,255,255,255,255,  0,
           0,  0,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 16, 32, 48,
          64, 48, 32, 16,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 16, 32, 48,
          64, 48, 32, 16,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 16, 32, 48,
          64,112,160,208,255,255,255,255,255,255,255,255,255,255,255,255 }};
   short col6[3][256] = 
       {{  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  2,  4,  6,  8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30,
          32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60, 62,
          64, 66, 68, 70, 72, 74, 76, 78, 80, 82, 84, 86, 88, 90, 92, 94,
          96, 98,100,102,104,106,108,110,112,114,116,118,120,122,124,126,
         128,130,132,134,136,138,140,142,144,146,148,150,152,154,156,158,
         160,162,164,166,168,170,172,174,176,178,180,182,184,186,188,190,
         192,194,196,198,200,202,204,206,208,210,212,214,216,218,220,222,
         224,226,228,230,232,234,236,238,240,242,244,246,248,250,252,254 },
        {  0,  2,  4,  6,  8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30,
          32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60, 62,
          64, 66, 68, 70, 72, 74, 76, 78, 80, 82, 84, 86, 88, 90, 92, 94,
          96, 98,100,102,104,106,108,110,112,114,116,118,120,122,124,126,
         128,130,132,134,136,138,140,142,144,146,148,150,152,154,156,158,
         160,162,164,166,168,170,172,174,176,178,180,182,184,186,188,190,
         192,194,196,198,200,202,204,206,208,210,212,214,216,218,220,222,
         224,226,228,230,232,234,236,238,240,242,244,246,248,250,252,254,
         254,252,250,248,246,244,242,240,238,236,234,232,230,228,226,224,
         222,220,218,216,214,212,210,208,206,204,202,200,198,196,194,192,
         190,188,186,184,182,180,178,176,174,172,170,168,166,164,162,160,
         158,156,154,152,150,148,146,144,142,140,138,136,134,132,130,128,
         126,124,122,120,118,116,114,112,110,108,106,104,102,100, 98, 96,
          94, 92, 90, 88, 86, 84, 82, 80, 78, 76, 74, 72, 70, 68, 66, 64,
          62, 60, 58, 56, 54, 52, 50, 48, 46, 44, 42, 40, 38, 36, 34, 32,
          30, 28, 26, 24, 22, 20, 18, 16, 14, 12, 10,  8,  6,  4,  2,  0 },
        {254,252,250,248,246,244,242,240,238,236,234,232,230,228,226,224,
         222,220,218,216,214,212,210,208,206,204,202,200,198,196,194,192,
         190,188,186,184,182,180,178,176,174,172,170,168,166,164,162,160,
         158,156,154,152,150,148,146,144,142,140,138,136,134,132,130,128,
         126,124,122,120,118,116,114,112,110,108,106,104,102,100, 98, 96,
          94, 92, 90, 88, 86, 84, 82, 80, 78, 76, 74, 72, 70, 68, 66, 64,
          62, 60, 58, 56, 54, 52, 50, 48, 46, 44, 42, 40, 38, 36, 34, 32,
          30, 28, 26, 24, 22, 20, 18, 16, 14, 12, 10,  8,  6,  4,  2,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 }};
   short col7[3][256] = 
       {{  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 
          16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 
          32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 
          48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 
          64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 
          80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 
          96, 97, 98, 99,100,101,102,103,104,105,106,107,108,109,110,111,
         112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,
         128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,
         144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,
         160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,
         176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,
         192,193,194,195,196,197,198,199,200,201,202,203,204,205,206,207,
         208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,
         224,225,226,227,228,229,230,231,232,233,234,235,236,237,238,239,
         240,241,242,243,244,245,246,247,248,249,250,251,252,253,254,255 },
        {  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 
          16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 
          32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 
          48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 
          64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 
          80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 
          96, 97, 98, 99,100,101,102,103,104,105,106,107,108,109,110,111,
         112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,
         127,126,125,124,123,122,121,120,119,118,117,116,115,114,113,112,
         111,110,109,108,107,106,105,104,103,102,101,100, 99, 98, 97, 96,
         95, 94, 93, 92, 91, 90, 89, 88, 87, 86, 85, 84, 83, 82, 81,  80,
         79, 78, 77, 76, 75, 74, 73,  72,71, 70, 69, 68, 67, 66, 65,  64,
         63, 62, 61, 60, 59, 58, 57, 56, 55, 54, 53, 52, 51, 50, 49,  48,
         47, 46, 45, 44, 43, 42, 41, 40, 39, 38, 37, 36, 35, 34, 33,  32,
         31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17,  16,
         15, 14, 13, 12, 11, 10, 9,  8,  7,  6,  5,  4,  3,  2,  1,   0  },
        {255,254,253,252,251,250,249,248,247,246,245,244,243,242,241,240,
         239,238,237,236,235,234,233,232,231,230,229,228,227,226,225,224,
         223,222,221,220,219,218,217,216,215,214,213,212,211,210,209,208,
         207,206,205,204,203,202,201,200,199,198,197,196,195,194,193,192,
         191,190,189,188,187,186,185,184,183,182,181,180,179,178,177,176,
         175,174,173,172,171,170,169,168,167,166,165,164,163,162,161,160,
         159,158,157,156,155,154,153,152,151,150,149,148,147,146,145,144,
         143,142,141,140,139,138,137,136,135,134,133,132,131,130,129,128,
         127,126,125,124,123,122,121,120,119,118,117,116,115,114,113,112,
         111,110,109,108,107,106,105,104,103,102,101,100, 99, 98, 97, 96,
         95, 94, 93, 92, 91, 90, 89, 88, 87, 86, 85, 84, 83, 82, 81,  80,
         79, 78, 77, 76, 75, 74, 73,  72,71, 70, 69, 68, 67, 66, 65,  64,
         63, 62, 61, 60, 59, 58, 57, 56, 55, 54, 53, 52, 51, 50, 49,  48,
         47, 46, 45, 44, 43, 42, 41, 40, 39, 38, 37, 36, 35, 34, 33,  32,
         31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17,  16,
         15, 14, 13, 12, 11, 10, 9,  8,  7,  6,  5,  4,  3,  2,  1,   0  }};
   short col8[3][256] = 
       {{  0,  2,  4,  6,  8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30,
          32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60, 62,
          64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 
          80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 
          96, 97, 98, 99,100,101,102,103,104,105,106,107,108,109,110,111,
         112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,
         128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,
         144,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         192,192,193,193,194,194,195,195,196,196,197,197,198,198,199,199,
         200,200,201,201,202,202,203,203,204,204,205,205,206,206,207,207,
         208,208,209,209,210,210,211,211,212,212,213,213,214,214,215,215,
         216,216,217,217,218,218,219,219,220,220,221,221,222,222,223,223,
         224,224,225,225,226,226,227,227,228,228,229,229,230,230,231,231,
         232,232,233,233,234,234,235,235,236,236,237,237,238,238,239,239,
         240,241,242,243,244,245,246,247,248,249,250,251,252,253,254,255 },
        {  0,  2,  4,  6,  8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30,
          32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60, 62,
          64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 
          80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 
          96, 97, 98, 99,100,101,102,103,104,105,106,107,108,109,110,111,
         112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,
         128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,
         144,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         191,190,189,188,187,186,185,184,183,182,181,180,179,178,177,176,
         175,174,173,172,171,170,169,168,167,166,165,164,163,162,161,160,
         158,156,154,152,150,148,146,144,142,140,138,136,134,132,130,128,
         126,124,122,120,118,116,114,112,110,108,106,104,102,100, 98, 96,
          94, 92, 90, 88, 86, 84, 82, 80, 78, 76, 74, 72, 70, 68, 66, 64,
          62, 60, 58, 56, 54, 52, 50, 48, 46, 44, 42, 40, 38, 36, 34, 32,
          30, 28, 26, 24, 22, 20, 18, 16, 14, 12, 10,  8,  6,  4,  2,  0 },
        {255,254,253,252,251,250,249,248,247,246,245,244,243,242,241,240,
         239,238,237,236,235,234,233,232,231,230,229,228,227,226,225,224,
         223,222,221,220,219,218,217,216,215,214,213,212,211,210,209,208,
         207,206,205,204,203,202,201,200,199,198,197,196,195,194,193,192,
         191,190,189,188,187,186,185,184,183,182,181,180,179,178,177,176,
         175,174,173,172,171,170,169,168,167,166,165,164,163,162,161,160,
         159,158,157,156,155,154,153,152,151,150,149,148,147,146,145,144,
         143,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         191,190,189,188,187,186,185,184,183,182,181,180,179,178,177,176,
         175,174,173,172,171,170,169,168,167,166,165,164,163,162,161,160,
         158,156,154,152,150,148,146,144,142,140,138,136,134,132,130,128,
         126,124,122,120,118,116,114,112,110,108,106,104,102,100, 98, 96,
          94, 92, 90, 88, 86, 84, 82, 80, 78, 76, 74, 72, 70, 68, 66, 64,
          62, 60, 58, 56, 54, 52, 50, 48, 46, 44, 42, 40, 38, 36, 34, 32,
          30, 28, 26, 24, 22, 20, 18, 16, 14, 12, 10,  8,  6,  4,  2,  0 }};
   short col9[3][256] = 
       {{  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 
          16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 
          32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 
          48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 
          64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 
          80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 
          96, 97, 98, 99,100,101,102,103,104,105,106,107,108,109,110,111,
         112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,
         128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,
         144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,
         160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,
         176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,
         192,193,194,195,196,197,198,199,200,201,202,203,204,205,206,207,
         208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,
         224,225,226,227,228,229,230,231,232,233,234,235,236,237,238,239,
         240,241,242,243,244,245,246,247,248,249,250,251,252,253,254,255 },
        {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
        {255,254,253,252,251,250,249,248,247,246,245,244,243,242,241,240,
         239,238,237,236,235,234,233,232,231,230,229,228,227,226,225,224,
         223,222,221,220,219,218,217,216,215,214,213,212,211,210,209,208,
         207,206,205,204,203,202,201,200,199,198,197,196,195,194,193,192,
         191,190,189,188,187,186,185,184,183,182,181,180,179,178,177,176,
         175,174,173,172,171,170,169,168,167,166,165,164,163,162,161,160,
         159,158,157,156,155,154,153,152,151,150,149,148,147,146,145,144,
         143,142,141,140,139,138,137,136,135,134,133,132,131,130,129,128,
         127,126,125,124,123,122,121,120,119,118,117,116,115,114,113,112,
         111,110,109,108,107,106,105,104,103,102,101,100, 99, 98, 97, 96,
         95, 94, 93, 92, 91, 90, 89, 88, 87, 86, 85, 84, 83, 82, 81,  80,
         79, 78, 77, 76, 75, 74, 73,  72,71, 70, 69, 68, 67, 66, 65,  64,
         63, 62, 61, 60, 59, 58, 57, 56, 55, 54, 53, 52, 51, 50, 49,  48,
         47, 46, 45, 44, 43, 42, 41, 40, 39, 38, 37, 36, 35, 34, 33,  32,
         31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17,  16,
         15, 14, 13, 12, 11, 10, 9,  8,  7,  6,  5,  4,  3,  2,  1,   0  }};
   short col10[3][256] =  
       {{  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
           0,  4,  8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60,
          64, 68, 72, 76, 80, 84, 88, 92, 96,100,104,108,112,116,120,124,
         128,132,136,140,144,148,152,156,160,164,168,172,176,180,184,188,
         192,196,200,204,208,212,216,220,224,228,232,236,240,244,248,252,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,252,248,244,240,236,232,228,224,220,216,212,208,204,200,196,
         192,188,184,180,176,172,168,164,160,156,152,148,144,140,136,132 },
        {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  4,  8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60,
          64, 68, 72, 76, 80, 84, 88, 92, 96,100,104,108,112,116,120,124,
         128,132,136,140,144,148,152,156,160,164,168,172,176,180,184,188,
         192,196,200,204,208,212,216,220,224,228,232,236,240,244,248,252,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,252,248,244,240,236,232,228,224,220,216,212,208,204,200,196,
         192,188,184,180,176,172,168,164,160,156,152,148,144,140,136,132,
         128,124,120,116,112,108,104,100, 96, 92, 88, 84, 80, 76, 72, 68,
          64, 60, 56, 52, 48, 44, 40, 36, 32, 28, 24, 20, 16, 12,  8,  4,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
        {128,132,136,140,144,148,152,156,160,164,168,172,176,180,184,188,
         192,196,200,204,208,212,216,220,224,228,232,236,240,244,248,252,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,252,248,244,240,236,232,228,224,220,216,212,208,204,200,196,
         192,188,184,180,176,172,168,164,160,156,152,148,144,140,136,132,
         128,124,120,116,112,108,104,100, 96, 92, 88, 84, 80, 76, 72, 68,
          64, 60, 56, 52, 48, 44, 40, 36, 32, 28, 24, 20, 16, 12,  8,  4,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0  }};
   short col11[3][256] =  
       {{255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         224,192,160,128, 96, 64, 32,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
           0,  0,  0,  0,  0,  0, 32, 64, 96,128,160,192,224,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255 },
        {255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         224,192,160,128, 96, 64, 32,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
           0,  0,  0, 32, 64, 96,128,160,192,224,255,255,255,255,255,255, 
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,224,192,160,128, 96, 64, 32,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0  },
        {255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,224,192,160,128, 96, 64, 32,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 }};
   short col12[3][256] = 
       {{  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
         112,112,112,112,112,112,112,112,112,112,112,112,112,112,112,112,
         112,112,112,112,112,112,112,112,112,112,112,112,112,112,112,112, 
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255 },
        {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
          92, 92, 92, 92, 92, 92, 92, 92, 92, 92, 92, 92, 92, 92, 92, 92,
          92, 92, 92, 92, 92, 92, 92, 92, 92, 92, 92, 92, 92, 92, 92, 92,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         196,196,196,196,196,196,196,196,196,196,196,196,196,196,196,196,
         196,196,196,196,196,196,196,196,196,196,196,196,196,196,196,196, 
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
         160,160,160,160,160,160,160,160,160,160,160,160,160,160,160,160,
         160,160,160,160,160,160,160,160,160,160,160,160,160,160,160,160,
         240,240,240,240,240,240,240,240,240,240,240,240,240,240,240,240 },
        {128,128,128,128,128,128,128,128,128,128,128,128,128,128,128,128,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255 }};
   short col13[3][256] = 
       {{  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 14, 28, 42,
          56, 70, 84, 98,112,112,112,112,112,112,112,112,112,112,112,112,
         112,112,112,112,112,112,112,112,112,112,112,112,112, 98, 84, 70, 
          56, 42, 28, 14,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 16, 32, 48,
          64, 48, 32, 16,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 16, 32, 48,
          64,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,208,160,112,
          64,112,160,208,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,208,160,112,
          64,112,160,208,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,208,160,112,
          64,112,160,208,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,208,160,112,
          64,112,160,208,255,255,255,255,255,255,255,255,255,255,255,255 },
        {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 23, 46, 69,
          92,115,138,161, 92, 92, 92, 92, 92, 92, 92, 92, 92, 92, 92, 92,
          92, 92, 92, 92, 92, 92, 92, 92, 92, 92, 92, 92, 92, 85, 78, 71,
          64,112,160,208,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,208,160,112,
          64,112,160,208,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,208,160,112,
          64,112,160,208,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,208,160,112,
          64,84, 104,124,144,144,144,144,144,144,144,144,144,144,144,144,
         144,144,144,144,144,144,144,144,144,144,144,144,144,124,104,84, 
          64,66, 68, 70,  72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72,
          72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72,70, 68, 66, 
          64, 48, 32, 16,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 16, 32, 48,
          64, 48, 32, 16,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
        {128,128,128,128,128,128,128,128,128,128,128,128,128,112, 96, 80,
          64,112,160,208,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,208,160,112,
          64,112,160,208,255,255,255,255,255,255,255,255,255,255,255,255,
         255,255,255,255,255,255,255,255,255,255,255,255,255,208,160,112,
          64, 48, 32, 16,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 16, 32, 48,
          64, 48, 32, 16,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 16, 32, 48,
          64, 48, 32, 16,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 16, 32, 48,
          64, 48, 32, 16,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 16, 32, 48,
          64, 48, 32, 16,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 16, 32, 48,
          64,112,160,208,255,255,255,255,255,255,255,255,255,255,255,255 }};
       

   switch (colmap) {
     default: 
       for (i=0;i<=2;i++) for (j=0;j<=255;j++)  
         icmval[j+256*i] = (double) blank[j];
       break;
     case 0: 
       for (i=0;i<=2;i++) for (j=0;j<=255;j++)  
         icmval[j+256*i] = (double) gry0[j];
       break;
     case 1: 
       for (i=0;i<=2;i++) for (j=0;j<=255;j++)  
         icmval[j+256*i] = (double) gry1[j];
       break;
     case 2: 
       for (i=0;i<=2;i++) for (j=0;j<=255;j++)  
         icmval[j+256*i] = (double) gry2[j];
       break;
     case 3: 
       for (i=0;i<=2;i++) for (j=0;j<=255;j++)  
         icmval[j+256*i] = (double) gry3[j];
       break;
     case 4: 
       for (i=0;i<=2;i++) for (j=0;j<=255;j++)  
         icmval[j+256*i] = (double) gry4[j];
       break;
     case 5: 
       for (i=0;i<=2;i++) for (j=0;j<=255;j++)  
         icmval[j+256*i] = (double) gry5[j];
       break;
     case 6: 
       for (i=0;i<=2;i++) for (j=0;j<=255;j++)  
         icmval[j+256*i] = (double) gry6[j];
       break;
     case 7: 
       for (i=0;i<=2;i++) for (j=0;j<=255;j++)  
         icmval[j+256*i] = (double) gry7[j];
       break;
     case 8: 
       for (i=0;i<=2;i++) for (j=0;j<=255;j++)  
         icmval[j+256*i] = (double) gry8[j];
       break;
     case 9: 
       for (i=0;i<=2;i++) for (j=0;j<=255;j++)  
         icmval[j+256*i] = (double) gry9[j];
       break;
     case 10: 
       for (i=0;i<=2;i++) for (j=0;j<=255;j++)  
         icmval[j+256*i] = (double) col0[i][j];
       break;
     case 11: 
       for (i=0;i<=2;i++) for (j=0;j<=255;j++)  
         icmval[j+256*i] = (double) col1[i][j];
       break;
     case 12: 
       for (i=0;i<=2;i++) for (j=0;j<=255;j++)  
         icmval[j+256*i] = (double) col2[i][j];
       break;
     case 13: 
       for (i=0;i<=2;i++) for (j=0;j<=255;j++)  
         icmval[j+256*i] = (double) col3[i][j];
       break;
     case 14: 
       for (i=0;i<=2;i++) for (j=0;j<=255;j++)  
         icmval[j+256*i] = (double) col4[i][j];
       break;
     case 15: 
       for (i=0;i<=2;i++) for (j=0;j<=255;j++)  
         icmval[j+256*i] = (double) col5[i][j];
       break;
     case 16: 
       for (i=0;i<=2;i++) for (j=0;j<=255;j++)  
         icmval[j+256*i] = (double) col6[i][j];
       break;
     case 17: 
       for (i=0;i<=2;i++) for (j=0;j<=255;j++)  
         icmval[j+256*i] = (double) col7[i][j];
       break;
     case 18: 
       for (i=0;i<=2;i++) for (j=0;j<=255;j++)  
         icmval[j+256*i] = (double) col8[i][j];
       break;
     case 19: 
       for (i=0;i<=2;i++) for (j=0;j<=255;j++)  
         icmval[j+256*i] = (double) col9[i][j];
       break;
     case 20: 
       for (i=0;i<=2;i++) for (j=0;j<=255;j++)  
         icmval[j+256*i] = (double) gry10[j];
       break;
     case 21: 
       for (i=0;i<=2;i++) for (j=0;j<=255;j++)  
         icmval[j+256*i] = (double) gry11[j];
       break;
     case 22: 
       for (i=0;i<=2;i++) for (j=0;j<=255;j++)  
         icmval[j+256*i] = (double) gry12[j];
       break;
     case 23: 
       for (i=0;i<=2;i++) for (j=0;j<=255;j++)  
         icmval[j+256*i] = (double) gry13[j];
       break;
     case 24: 
       for (i=0;i<=2;i++) for (j=0;j<=255;j++)  
         icmval[j+256*i] = (double) gry14[j];
       break;
     case 25: 
       for (i=0;i<=2;i++) for (j=0;j<=255;j++)  
         icmval[j+256*i] = (double) gry15[j];
       break;
     case 26: 
       for (i=0;i<=2;i++) for (j=0;j<=255;j++)  
         icmval[j+256*i] = (double) gry16[j];
       break;
     case 27: 
       for (i=0;i<=2;i++) for (j=0;j<=255;j++)  
         icmval[j+256*i] = (double) gry17[j];
       break;
     case 28: 
       for (i=0;i<=2;i++) for (j=0;j<=255;j++)  
         icmval[j+256*i] = (double) gry18[j];
       break;
     case 29: 
       for (i=0;i<=2;i++) for (j=0;j<=255;j++)  
         icmval[j+256*i] = (double) gry19[j];
       break;
     case 30: 
       for (i=0;i<=2;i++) for (j=0;j<=255;j++)  
         icmval[j+256*i] = (double) col10[i][j];
       break;
     case 31: 
       for (i=0;i<=2;i++) for (j=0;j<=255;j++)  
         icmval[j+256*i] = (double) col11[i][j];
       break;
     case 32: 
       for (i=0;i<=2;i++) for (j=0;j<=255;j++)  
         icmval[j+256*i] = (double) col12[i][j];
       break;
   }
}

/* =============================
   Get colour map
   =============================
*/
int getcolmapnumber()
{
  int colmap;
  char dum[MAXCHARLEN];

  printf("Choose one of the following colour maps:\n");
  printf("\t(-1) Read `colmap.xpt` having 3 curves of 256 points each\n");
  printf("\t(0) Gray scale - linear black to white\n");
  printf("\t(1) Gray scale - 9 levels: black to white, white middle\n");
  printf("\t(2) Gray scale - linear black to white with white middle\n");
  printf("\t(3) Gray scale - linear white to black with white middle\n");
  printf("\t(4) Gray scale - linear white to black\n");
  printf("\t(5) Gray scale - 9 levels: black to white, gray middle\n");
  printf("\t(6) Gray scale - 9 levels: white to black, gray middle\n");
  printf("\t(7) Gray scale - 5 levels: black,drk,white,light, v. light\n");
  printf("\t(8) Gray scale - 3 levels: drk,white,light\n");
  printf("\t(9) Gray scale - 2 levels: half gray, half white\n");
  printf("\t(10) Colour scale - 9 levels: white centre\n");
  printf("\t(11) Colour scale - 9 levels: yellow centre\n");
  printf("\t(12) Colour scale - 2 levels: G-O/white centre\n");
  printf("\t(13) Colour scale - 9 levels: white bottom\n");
  printf("\t(14) Colour scale - 9 levels: NCAR rainbow\n");
  printf("\t(15) Colour scale - 9 levels: yellow zero\n");
  printf("\t(16) Colour scale - linear blue->green->red\n");
  printf("\t(17) Colour scale - linear blue->gray->red\n");
  printf("\t(18) Colour scale - linear blue->red, white centre\n");
  printf("\t(19) Colour scale - linear blue->magenta->red\n");
  printf("\t(20) Gray scale - 7 levels: white centre\n");
  printf("\t(21) Gray scale - 13 levels: white centre\n");
  printf("\t(22) Gray scale - 2 levels: half black, half white\n");
  printf("\t(23) Gray scale - 5 levels: dark gray centre\n");
  printf("\t(24) Gray scale - 5 levels: white,drk,gray,lt.gray,white\n");
  printf("\t(25) Gray scale - 5 levels: white,md lt gry,lt gry,gry,blk\n");
  printf("\t(26) Gray scale - 3 levels: white->blk->white\n");
  printf("\t(27) Gray scale - 5 levels: white,drk,gry,lt.gray,ltlt gray\n");
  printf("\t(28) Gray scale - 5 levels: black,drk.gry,gry,lt.gray,white\n");
  printf("\t(29) Gray scale - 5 levels: white,lt.gry,gry,drk.gry,black\n");
  printf("\t(30) Colour scale - matlab rainbow: prp,bl,cyn,grn,ylw,rd,mgta\n");
  printf("\t(31) Colour scale - 5 levels: white,bl,cyn,ylw,rd\n");
  printf("\t(32) Colour scale - 9 levels: bl,prp,cy,gn,yl,or,orrd,rd,mg\n");
  scanf("%d",&colmap);
  if (colmap<-1 ||colmap>MAXCOLMAP) 
    myerror("GetColourMap: colour map chosen does not exist");

  return(colmap);
}

/* ===========================================================================
   GETCOLMAP: APR 98

   Define RGB colour map using XPlot file or predefined data
   ===========================================================================
*/
void getcolmap(colmap,cm)
   int colmap;
   short **cm;
{
   int j;
   int nc,*npts;
   double *icm,*icmval,*xzrng;
   char sfpr[50];

   npts = ivector(1,3);
   xzrng = dvector(1,4);
   icm = dvector(0,3*COLMAPSIZE-1);
   icmval = dvector(0,3*COLMAPSIZE-1);

   if (colmap<0) {
     printf("Reading colourmap from XPlot file: colmap.xpt\n");
     printf("This should contain 3 curves with 256 points in each.\n");

     sprintf(sfpr,"colmap.xpt");
   
     readXPlot_noprompt(sfpr, icm-1, icmval-1, &nc, npts, xzrng, 3, 768);
     if (nc!=3) 
        myerror("GETCOLMAP: not enough curves in colmap file");
     if (npts[1]!=COLMAPSIZE || npts[2]!=COLMAPSIZE || npts[3]!=COLMAPSIZE) 
        myerror("GETCOLMAP: not enough points in colmap file");
     if ( xzrng[1]<0.0 || xzrng[2]>((double) (COLMAPSIZE-1)) ||
          xzrng[3]<0.0 || xzrng[4]>((double) (COLMAPSIZE-1)) )
        myerror("GETCOLMAP: values in colmap file are out of range");
   } else {
     defcolmap(icmval,colmap);
   }

   for (j=0;j<=255;j++) {
      cm[0][j] = (short) icmval[j];
      cm[1][j] = (short) icmval[j+COLMAPSIZE];
      cm[2][j] = (short) icmval[j+2*COLMAPSIZE];
   }

   free_ivector(npts,1,3);
   free_dvector(xzrng,1,4);
   free_dvector(icm,0,3*COLMAPSIZE-1);
   free_dvector(icmval,0,3*COLMAPSIZE-1);
}

/* =============================
   Get colour from user input 
   =============================
*/
rgbclr getcolour()
{
  int cflg;
  char dum[MAXCHARLEN];
  rgbclr clr;

  printf("Enter colour: \n");
  printf("\t (-1) define (red, green, blue) values explicitly\n");
  printf("\t (0) black\n");
  printf("\t (1) white\n");
  printf("\t (2) red\n");
  printf("\t (3) green\n");
  printf("\t (4) blue\n");
  printf("\t (5) cyan\n");
  printf("\t (6) magenta\n");
  printf("\t (7) yellow\n");
  scanf("%d",&cflg); gets(dum);

  switch (cflg) {
  default : /* do same thing as case -1 */
  case -1 : printf("red?   [255] \n"); scanf("%hi", &clr.r); gets(dum);
            printf("green? [255] \n"); scanf("%hi", &clr.g); gets(dum);
            printf("blue?  [255] \n"); scanf("%hi", &clr.b); gets(dum);
            break;
  case 0  : clr.r =   0; clr.g =   0; clr.b =   0;
            break;
  case 1  : clr.r = 255; clr.g = 255; clr.b = 255;
            break;
  case 2  : clr.r = 255; clr.g =   0; clr.b =   0;
            break;
  case 3  : clr.r =   0; clr.g = 255; clr.b =   0;
            break;
  case 4  : clr.r =   0; clr.g =   0; clr.b =   0;
            break;
  case 5  : clr.r =   0; clr.g = 255; clr.b = 255;
            break;
  case 6  : clr.r = 255; clr.g =   0; clr.b = 255;
            break;
  case 7  : clr.r = 255; clr.g = 255; clr.b =   0;
            break;
  }
  return(clr);
}

/* ===========================================================================
   Img2RGB: June 16/98

   Write matrix of intensity data to RGB file using colour map.
   Integer image data range from 0 to 255.
   ===========================================================================
*/
void Img2RGB(imat,nx,ny,colmap,rgbmat,ni,nj)
   int nx,ny,colmap,ni,nj;
   short **imat;
   short ***rgbmat;
{
   int ix,iy;
   short ic;
   short **cm;

   if (ni<ny || nj<nx) myerror("IMG2RGB: rgb matrix too small");
   if (ni!=ny || nj!=nx) 
     printf("IMG2RGB (WARNING): rgbmat should be transpose of imat");

   cm = smatrix2(0,2,0,COLMAPSIZE-1);

   /* Define RGB colour map */
   getcolmap(colmap,cm);

   /* Write data by mapping "intensity" to rgb values in hex 
      NOTE: Image is flipped so that (ix,iy) origin at bottom-left
            becomes array origin (i=iy (col), j=ix (row)) at top-left
   */
   for (iy=0;iy<ny;iy++) for (ix=0;ix<nx;ix++) {
         ic = imat[ix][iy];
         rgbmat[0][ny-iy-1][ix]= cm[0][ic];
         rgbmat[1][ny-iy-1][ix]= cm[1][ic];
         rgbmat[2][ny-iy-1][ix]= cm[2][ic];
   }

   free_smatrix2(cm,0,2,0,COLMAPSIZE-1);
}

/* ===========================================================================
   WRITEGRYPSHEADER: APR 98

   Write Grayscale Postscript header
   ===========================================================================
*/
void writeGryPSHeader(fpw,sfpw,ni,nj)
   FILE *fpw;
   char *sfpw;
   int ni,nj;
{
   fprintf(fpw,"%%!PS-Adobe-2.0 EPSF-2.0\n");
   fprintf(fpw,"%%%%Title: %s\n", sfpw);
   fprintf(fpw,"%%%%Creator: WriteGryPS 01/01  -  by Bruce Sutherland\n");
   fprintf(fpw,"%%%%BoundingBox: 0 0 %d %d\n",nj,ni);
   fprintf(fpw,"%%%%Pages: 1\n");
   fprintf(fpw,"%%%%DocumentFonts:\n");
   fprintf(fpw,"%%%%EndComments\n");
   fprintf(fpw,"%%%%EndProlog\n");
   fprintf(fpw,"\n");
   fprintf(fpw,"%%%%Page: 1 1\n");
   fprintf(fpw,"\n");
   fprintf(fpw,"%% remember original state\n");
   fprintf(fpw,"/origstate save def\n");
   fprintf(fpw,"\n");
   fprintf(fpw,"%% build a temporary dictionary\n");
   fprintf(fpw,"20 dict begin\n");
   fprintf(fpw,"\n");
   fprintf(fpw,"%% define string to hold a scanline's worth of data\n");
   fprintf(fpw,"/pix %d string def\n",nj);
   fprintf(fpw,"\n");
   fprintf(fpw,"%% define space for color conversions\n");
   fprintf(fpw,"/grays %d string def  %% space for gray scale line\n",nj);
   fprintf(fpw,"/npixls 0 def\n");
   fprintf(fpw,"/rgbindx 0 def\n");
   fprintf(fpw,"\n");
   fprintf(fpw,"%% lower left corner\n");
   fprintf(fpw,"0 0 translate\n");
   fprintf(fpw,"\n");
   fprintf(fpw,"%% size of image (on paper, in 1/72inch coords)\n");
   fprintf(fpw,"%9.4f %9.4f scale\n", (double) nj, (double) ni);
   fprintf(fpw,"\n");
   fprintf(fpw,"%d %d 8                       %% dimensions of data\n",nj,ni);
   fprintf(fpw,"[%d 0 0 -%d 0 %d]            %% mapping matrix\n",nj,ni,ni);
   fprintf(fpw,"{currentfile pix readhexstring pop}\n");
   fprintf(fpw,"image\n");
   fprintf(fpw,"\n");
}

/* ===========================================================================
   WRITECOLPSHEADER: APR 98

   Write Colour Postscript header
   ===========================================================================
*/
void writeColPSHeader(fpw,sfpw,ni,nj)
   FILE *fpw;
   char *sfpw;
   int ni,nj;
{
   fprintf(fpw,"%%!PS-Adobe-2.0 EPSF-2.0\n");
   fprintf(fpw,"%%%%Title: %s\n", sfpw);
   fprintf(fpw,"%%%%Creator: WriteColPS 04/98  -  by Bruce Sutherland\n");
   fprintf(fpw,"%%%%BoundingBox: 0 0 %d %d\n",nj,ni);
   fprintf(fpw,"%%%%Pages: 1\n");
   fprintf(fpw,"%%%%DocumentFonts:\n");
   fprintf(fpw,"%%%%EndComments\n");
   fprintf(fpw,"%%%%EndProlog\n");
   fprintf(fpw,"\n");
   fprintf(fpw,"%%%%Page: 1 1\n");
   fprintf(fpw,"\n");
   fprintf(fpw,"%% remember original state\n");
   fprintf(fpw,"/origstate save def\n");
   fprintf(fpw,"\n");
   fprintf(fpw,"%% build a temporary dictionary\n");
   fprintf(fpw,"20 dict begin\n");
   fprintf(fpw,"\n");
   fprintf(fpw,"%% define string to hold a scanline's worth of data\n");
   fprintf(fpw,"/pix %d string def\n",3*nj);
   fprintf(fpw,"\n");
   fprintf(fpw,"%% define space for color conversions\n");
   fprintf(fpw,"/grays %d string def  %% space for gray scale line\n",nj);
   fprintf(fpw,"/npixls 0 def\n");
   fprintf(fpw,"/rgbindx 0 def\n");
   fprintf(fpw,"\n");
   fprintf(fpw,"%% lower left corner\n");
   fprintf(fpw,"0 0 translate\n");
   fprintf(fpw,"\n");
   fprintf(fpw,"%% size of image (on paper, in 1/72inch coords)\n");
   fprintf(fpw,"%9.4f %9.4f scale\n", (double) nj, (double) ni);
   fprintf(fpw,"\n");
   fprintf(fpw,"%% define 'colorimage' if it isn't defined\n");
   fprintf(fpw,"%%   ('colortogray' and 'mergeprocs' come from xwd2ps\n");
   fprintf(fpw,"%%     via xgrab)\n");
   fprintf(fpw,"/colorimage where   %% do we know about 'colorimage'?\n");
   fprintf(fpw,"{ pop }           %% yes: pop off the 'dict' returned\n");
   fprintf(fpw,"{                 %% no:  define one\n");
   fprintf(fpw,"/colortogray {  %% define an RGB->I function\n");
   fprintf(fpw,"/rgbdata exch store    %% call input 'rgbdata'\n");
   fprintf(fpw,"rgbdata length 3 idiv\n");
   fprintf(fpw,"/npixls exch store\n");
   fprintf(fpw,"/rgbindx 0 store\n");
   fprintf(fpw,"0 1 npixls 1 sub {\n");
   fprintf(fpw,"grays exch\n");
   fprintf(fpw,"rgbdata rgbindx       get 20 mul    %% Red\n");
   fprintf(fpw,"rgbdata rgbindx 1 add get 32 mul    %% Green\n");
   fprintf(fpw,"rgbdata rgbindx 2 add get 12 mul    %% Blue\n");
   fprintf(fpw,"add add 64 idiv      %% I = .5G + .31R + .18B\n");
   fprintf(fpw,"put\n");
   fprintf(fpw,"/rgbindx rgbindx 3 add store\n");
   fprintf(fpw,"} for\n");
   fprintf(fpw,"grays 0 npixls getinterval\n");
   fprintf(fpw,"} bind def\n");
   fprintf(fpw,"\n");
   fprintf(fpw,"%% Utility procedure for colorimage operator.\n");
   fprintf(fpw,"%% This procedure takes two procedures off the\n");
   fprintf(fpw,"%% stack and merges them into a single procedure.\n");
   fprintf(fpw,"\n");
   fprintf(fpw,"/mergeprocs { %% def\n");
   fprintf(fpw,"dup length\n");
   fprintf(fpw,"3 -1 roll\n");
   fprintf(fpw,"dup\n");
   fprintf(fpw,"length\n");
   fprintf(fpw,"dup\n");
   fprintf(fpw,"5 1 roll\n");
   fprintf(fpw,"3 -1 roll\n");
   fprintf(fpw,"add\n");
   fprintf(fpw,"array cvx\n");
   fprintf(fpw,"dup\n");
   fprintf(fpw,"3 -1 roll\n");
   fprintf(fpw,"0 exch\n");
   fprintf(fpw,"putinterval\n");
   fprintf(fpw,"dup\n");
   fprintf(fpw,"4 2 roll\n");
   fprintf(fpw,"putinterval\n");
   fprintf(fpw,"} bind def\n");
   fprintf(fpw,"\n");
   fprintf(fpw,"/colorimage { %% def\n");
   fprintf(fpw,"pop pop     %% remove 'false 3' operands\n");
   fprintf(fpw,"{colortogray} mergeprocs\n");
   fprintf(fpw,"image\n");
   fprintf(fpw,"} bind def\n");
   fprintf(fpw,"} ifelse          %% end of 'false' case\n");
   fprintf(fpw,"\n");
   fprintf(fpw,"\n");
   fprintf(fpw,"\n");
   fprintf(fpw,"%d %d 8                       %% dimensions of data\n",nj,ni);
   fprintf(fpw,"[%d 0 0 -%d 0 %d]            %% mapping matrix\n",nj,ni,ni);
   fprintf(fpw,"{currentfile pix readhexstring pop}\n");
   fprintf(fpw,"false 3 colorimage\n");
   fprintf(fpw,"\n");
}

/* ===========================================================================
   WRITECOLPSFOOT: APR 98

   Write Colour Postscript foot
   ===========================================================================
*/
void writePSFoot(fpw)
   FILE *fpw;
{
   fprintf(fpw,"\n");
   fprintf(fpw,"showpage\n");
   fprintf(fpw,"\n");
   fprintf(fpw,"%% stop using temporary dictionary\n");
   fprintf(fpw,"end\n");
   fprintf(fpw,"\n");
   fprintf(fpw,"%% restore original state\n");
   fprintf(fpw,"origstate restore\n");
   fprintf(fpw,"\n");
   fprintf(fpw,"%%%%Trailer\n");
}

/* ===========================================================================
   WRITECOLPS: Apr 12/98

   Write matrix of data to colour Postscript file.
   Integer data range from 0 to 255.
   ===========================================================================
*/
void writeColPS(sfpw,rgbmat,ni,nj)
   int ni,nj;
   short ***rgbmat;
   char *sfpw;
{
   FILE *fpw;
   int i,j,col;
   short irgb;

   fpw = fopen( sfpw, "w" );

   /* Write header information */
   writeColPSHeader(fpw,sfpw,ni,nj);
   
   /* Write data by mapping "intensity" to rgb values in hex */
   for (i=0;i<ni;i++) {
      col=0;
      for (j=0;j<nj;j++) {
         for (irgb=0;irgb<=2;irgb++) {
           if (rgbmat[irgb][i][j]<16) {
             fprintf(fpw,"0%x", rgbmat[irgb][i][j]);
           } else {
             fprintf(fpw,"%2x", rgbmat[irgb][i][j]);
           }
         }
         if ( (col += 6) >= 72 ) { fprintf(fpw,"\n"); col=0; }
      }
      fprintf(fpw,"\n"); fflush(fpw);
   }

   /* Write foot */
   writePSFoot(fpw);
   fclose(fpw);
}

/* ===========================================================================
   WRITEGRYPS: Jan 23/01

   Write matrix of data to Grayscale Postscript file.
   Integer data range from 0 to 255.
   ===========================================================================
*/
void writeGryPS(sfpw,rgbmat,ni,nj)
   int ni,nj;
   short ***rgbmat;
   char *sfpw;
{
   FILE *fpw;
   int i,j,col;
   short irgb;

   fpw = fopen( sfpw, "w" );

   /* Write header information */
   writeGryPSHeader(fpw,sfpw,ni,nj);
   
   /* Write data by mapping "intensity" to rgb values in hex */
   for (i=0;i<ni;i++) {
      col=0;
      for (j=0;j<nj;j++) {
         for (irgb=0;irgb<=0;irgb++) { /* only use 1st of 3 sets of data */
           if (rgbmat[irgb][i][j]<16) {
             fprintf(fpw,"0%x", rgbmat[irgb][i][j]);
           } else {
             fprintf(fpw,"%2x", rgbmat[irgb][i][j]);
           }
         }
         if ( (col += 2) >= 72 ) { fprintf(fpw,"\n"); col=0; }
      }
      fprintf(fpw,"\n"); fflush(fpw);
   }

   /* Write foot */
   writePSFoot(fpw);
   fclose(fpw);
}

/* ===========================================================================
   WRITERGB: March 19/99

   Write 3D matrix of intensity data to RGB file.
   Integer image data range from 0 to 255.

   NOTE: IMAGE files are like X,Y files (not matrix I,J files).  That is
         "iopen" expects nx==nj,ny==ni as 5th and 6th argument, respectively.
         Also, 3rd argument to putrow is row number from bottom up. 

         Because, I'm treating rgbmat like a matrix, I reflip the vertical.
   ===========================================================================
*/
void writeRGB(sfpw,rgbmat,ni,nj)
   int ni,nj;
   short ***rgbmat;
   char *sfpw;
{
   /* IMAGE *img; */
   int i;

   /* Open RGB file for writing */
   /* img = iopen(sfpw,"w",RLE(RGBTYPE),3,nj,ni,BWRGB);  */
   printf("WRITERGB (ERROR): cannot work with RGB files on non-sgi machine\n");
   exit(1);

   /* Write data by mapping "intensity" to rgb values in hex */
/*
   for (i=0;i<ni;i++) {
      putrow(img,rgbmat[0][i],ni-i-1,0);
      putrow(img,rgbmat[1][i],ni-i-1,1);
      putrow(img,rgbmat[2][i],ni-i-1,2);
   }

   iclose(img);
*/
}

/* ===========================================================================
   WRITEJPG: July 1/06

   Write 3D matrix of intensity data to JPG file.
   Integer image data range from 0 to 255.

   NOTE: IMAGE files are like X,Y files (not matrix I,J files).  That is
         "iopen" expects nx==nj,ny==ni as 5th and 6th argument, respectively.
         Also, 3rd argument to putrow is row number from bottom up. 

         Because, I'm treating rgbmat like a matrix, I reflip the vertical.
   ===========================================================================
*/
void writeJPG(sfpw,rgbmat,ni,nj)
   int ni,nj;
   short ***rgbmat;
   char *sfpw;
{

	//Needed Jpeg Structures
	struct jpeg_compress_struct cinfo;
	struct jpeg_error_mgr jerr;
	
	//Initialization
	cinfo.err = jpeg_std_error(&jerr);
	jpeg_create_compress(&cinfo);

	//Open the file and set as output destination
	FILE * outfile;

	if ((outfile = fopen(sfpw, "wb")) == NULL) {
		fprintf(stderr, "can't open %s\n", sfpw);
		exit(1);
	}
	jpeg_stdio_dest(&cinfo, outfile);

	//Set information (size, image type)
	cinfo.image_width =nj;
	cinfo.image_height = ni;
	cinfo.input_components = 3;
	cinfo.in_color_space = JCS_RGB;

	jpeg_set_defaults(&cinfo);
	jpeg_set_quality(&cinfo, JPGQUALITY, TRUE);

	//Start the compress
	jpeg_start_compress(&cinfo, TRUE);
	

	//Write code to change the format to a scaneline format
	JSAMPLE * pScanLines;
	pScanLines = (JSAMPLE *) malloc(ni * nj * 3 *sizeof(JSAMPLE));

	int i = 0,j = 0;	
	for(i = 0; i < ni; i++)
		for(j = 0; j < nj; j++)
		{
			pScanLines[(i*nj*3) + 3*j] = rgbmat[0][i][j];
			pScanLines[(i*nj*3) + 3*j + 1] = rgbmat[1][i][j];
			pScanLines[(i*nj*3) + 3*j + 2] = rgbmat[2][i][j];
		}
	
	JSAMPROW row_pointer[1];
	int row_stride = nj*3;

	//Loop to write the scanlines
	while (cinfo.next_scanline < cinfo.image_height) {
		row_pointer[0] = &pScanLines[cinfo.next_scanline * row_stride];
		jpeg_write_scanlines(&cinfo, row_pointer, 1);
	}
	
	//Finish up
	jpeg_finish_compress(&cinfo);
	fclose(outfile);

	jpeg_destroy_compress(&cinfo);
	free(pScanLines);
}


#undef CHARDIMEN
#undef MAXCHARLEN
#undef MAXCOLMAP
#undef COLMAPSIZE
#undef RGBTYPE
#undef BWRGB
#undef JPGQUALITY
