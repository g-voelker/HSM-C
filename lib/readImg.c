/* ===========================================================================
   Routines to read data in RGB, etc format
   Last Update: Mar 18/99
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

#define MAXCHARLEN 40  
#define RGBTYPE 1   /* =1 for 1 byte arrays (0..255); =2 for 2 byte arrays */
#define BWRGB 3     /* =1 for black and white; =3 for rgb */

/* ===========================================================================
   READRGB: Mar 18/99

   Read RGB file and store data in 3D matrix of short ints.
   Integer image data range from 0 to 255.

   NOTE: IMAGE files are like X,Y files (not matrix I,J files).  That is
         "iopen" expects nx==nj,ny==ni as 5th and 6th argument, respectively.
         Also, 3rd argument to getrow is row number from bottom up.

         Because, I'm treating rgbmat like a matrix, I reflip the vertical.
   ===========================================================================
*/
void readRGB(sfpr,rgbmat,ni,nj)
   int *ni,*nj;
   short ***rgbmat;
   char *sfpr;
{
/*   IMAGE *img; */
   int i,nx,ny;
   char filename[MAXCHARLEN];

   /* GET .rgb filee name */
   printf("READRGB (ERROR): Cannot work with RGB files except on SGI\n");
   exit(1);
/*
   printf("Enter RGB file name: ");
   scanf("%s",filename); sprintf(sfpr,"%s",filename);
   img = iopen(sfpr,"r");

   nx = img->xsize;
   ny = img->ysize;
   *ni=ny; *nj=nx;
*/

   /* Write data by mapping "intensity" to rgb values in hex */
   /* origin at bottom left-hand */
/*
   for (i=0;i<ny;i++) { 
      getrow(img,rgbmat[0][i],ny-i-1,0); 
      getrow(img,rgbmat[1][i],ny-i-1,1);
      getrow(img,rgbmat[2][i],ny-i-1,2);
   }

   iclose(img);
*/
}



/* ===========================================================================
   READJPG: July 11/06

   Read JPG file and store data in 3D matrix of short ints.
   Integer image data range from 0 to 255.
   ===========================================================================
*/
void readJPG(sfpr,rgbmat,ni,nj)
   int *ni,*nj;
   short ***rgbmat;
   char *sfpr;
{

  int nx,ny;
  char filename[MAXCHARLEN];

  //Allocate/Initialize the decompression object
  struct jpeg_decompress_struct cinfo;
  struct jpeg_error_mgr jerr;

  cinfo.err = jpeg_std_error(&jerr);
  jpeg_create_decompress(&cinfo);

  /* GET .jpg file name */

  printf("Enter JPG file name: ");
  scanf("%s",filename);
  sprintf(sfpr,"%s",filename);
  printf("\n");

  //Open the file
  FILE * infile;
  if((infile = fopen(filename, "rb")) == NULL ) {
    fprintf(stderr, "Can't open %s\n", filename);
    exit(1);
  }
  jpeg_stdio_src(&cinfo, infile);

  //Read the jpeg headers
  jpeg_read_header(&cinfo, TRUE);

  //Start Decompression
  jpeg_start_decompress(&cinfo);

  //Get the width and height
  nx = cinfo.output_width;
  ny = cinfo.output_height;
  *ni=ny;
  *nj=nx; 
  
  JSAMPARRAY buffer;
  buffer =(*cinfo.mem->alloc_sarray)
    ((j_common_ptr) &cinfo, JPOOL_IMAGE, cinfo.output_width * cinfo.output_components, 1);
  
  //Loop to read the scanlines into this buffer
  while(cinfo.output_scanline < cinfo.output_height) {
    jpeg_read_scanlines(&cinfo,buffer,1);
    int z;
    for(z=0;z<nx;z++) {
      rgbmat[0][cinfo.output_scanline-1][z] = buffer[0][z*3];
      rgbmat[1][cinfo.output_scanline-1][z] = buffer[0][z*3 + 1];
      rgbmat[2][cinfo.output_scanline-1][z] = buffer[0][z*3 + 2];
    }
  }

  //Complete Decompression
  jpeg_finish_decompress(&cinfo);
  jpeg_destroy_decompress(&cinfo);
  
  //Close the file
  fclose(infile);
}



/* ===========================================================================
   READJPGGRY: June 7/07

   Read JPG file and store data in 3D matrix of short ints.
   Integer image data range from 0 to 255.
   ===========================================================================
*/
void readJPGGRY(sfpr,grymat,ni,nj)
   int *ni,*nj;
   short **grymat;
   char *sfpr;
{

  int nx,ny;
  char filename[MAXCHARLEN];

  //Allocate/Initialize the decompression object
  struct jpeg_decompress_struct cinfo;
  struct jpeg_error_mgr jerr;

  cinfo.err = jpeg_std_error(&jerr);
  jpeg_create_decompress(&cinfo);

  /* GET .jpg file name */

  printf("Enter JPG file name: ");
  scanf("%s",filename);
  sprintf(sfpr,"%s",filename);
  printf("\n");

  //Open the file
  FILE * infile;
  if((infile = fopen(filename, "rb")) == NULL ) {
    fprintf(stderr, "Can't open %s\n", filename);
    exit(1);
  }
  jpeg_stdio_src(&cinfo, infile);

  //Read the jpeg headers
  jpeg_read_header(&cinfo, TRUE);

  //Start Decompression
  jpeg_start_decompress(&cinfo);

  //Get the width and height
  nx = cinfo.output_width;
  ny = cinfo.output_height;
  *ni=ny;
  *nj=nx; 
  
  JSAMPARRAY buffer;
  buffer =(*cinfo.mem->alloc_sarray)
    ((j_common_ptr) &cinfo, JPOOL_IMAGE, cinfo.output_width * cinfo.output_components, 1);
  
  //Loop to read the scanlines into this buffer
  while(cinfo.output_scanline < cinfo.output_height) {
    jpeg_read_scanlines(&cinfo,buffer,1);
    int z;
    for(z=0;z<nx;z++) {
      grymat[cinfo.output_scanline-1][z] = buffer[0][z];
    }
  }

  //Complete Decompression
  jpeg_finish_decompress(&cinfo);
  jpeg_destroy_decompress(&cinfo);
  
  //Close the file
  fclose(infile);
}

