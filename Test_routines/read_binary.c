/* this test the reading of unformatted binary files */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define FL 631855161

#define MAXCHUNK  536500000 /*largest number of 4-byte objects fread can handle */


int main(int argc, char **argv){
  char filename1[] = "/Users/frazin/tomography/bindata/ieuvi.AB.195.cr2077.bf4.1.50-1.76.Nr38";
  char filename2[] = "/Users/frazin/tomography/bindata/veuvi.AB.195.cr2077.bf4.1.50-1.76.Nr38";
  FILE *fp1, *fp2;
  int *d;
  float *h;
  int k,sz, chunk, filelength, nchunks, lastchunk;

  d = (int *) malloc(FL * sizeof(int) );
  if (d == NULL){
    fprintf(stdout,"null pointer (d).\n");
    exit(2);
  }

  h = (float *) malloc(FL * sizeof(int) );
  if (h == NULL){
    fprintf(stdout,"null pointer (h).\n");
    exit(2);
  }

  fp1 = fopen(filename1,"rb");
  if (fp1 == NULL){
    fprintf(stderr,"File not found: %s",filename1);
    exit(1);
  }

  fp2 = fopen(filename2,"rb");
  if (fp2 == NULL){
    fprintf(stderr,"File not found: %s",filename2);
    exit(1);
  }

  filelength = FL;

  nchunks = filelength / MAXCHUNK;
  lastchunk = filelength % MAXCHUNK;

  fprintf(stdout,"\n nchuncks = %d, lastchunk = %d\n",nchunks,lastchunk);
  fflush(stdout);

  chunk = 0;
  while (chunk < nchunks){
    fread(d + chunk*MAXCHUNK, sizeof (int)  , MAXCHUNK, fp1);
    fread(h + chunk*MAXCHUNK, sizeof (float), MAXCHUNK, fp2);
    chunk++;
  }
  fread(d + chunk*MAXCHUNK, sizeof (int)  , lastchunk, fp1);
  fread(h + chunk*MAXCHUNK, sizeof (float), lastchunk, fp2);
  fclose(fp1); fclose(fp2);


  fprintf(stdout,"\n");
  for (k = 0; k < 99; k++){
    fprintf(stdout,"(%d: %d %g) ",k,*(d+k),*(h+k));
      fflush(stdout);
  }
  k = MAXCHUNK;
  fprintf(stdout,"(%d: %d %g) ",k,*(d+k),*(h+k));
  k = FL - 1;
  fprintf(stdout,"(%d: %d %g) ",k,*(d+k),*(h+k));
  fprintf(stdout,"\n");

  return(0);
}



  


