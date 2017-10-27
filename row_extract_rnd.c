/*#include <malloc.h>*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "headers.h"

/* this code calculates for the matrix a 
 * c  =  d * a    and
 * cc = dd * a .
 * d is diagonal, containing only zeros or ones.  
 * dd = I - d.   d is chosen randomly.
 * the fraction of nonzero d is about 2^(-random_factor)
 *
 * c and cc are written to disk. 
 * cc contains most of a.
 * this also write the delta and y vectors to disk
 *
 * to compile:
 * gcc -o row_extract_rnd row_extract_rnd.c sparse.c -lm
 *
 */

int
main(int argc, char **argv)
{
  char filename_a[MAXPATH], filename_c[MAXPATH], filename_cc[MAXPATH];
  int random_factor;
  FILE *fid;
  struct matrix a;
  struct sparse aa, d, dd;
  int *diag_d, *diag_dd;
  float *y_c, *delta_c, *y_cc, *delta_cc;
  char fstr0[MAXPATH];
  int i, j, k;
  int nc3 = NBINS;


  fprintf(stdout,"fraction of matrix to be removed = 2^(-q), enter q: ");
  gets(filename_a);
  sscanf(filename_a,"%d",&random_factor);
  fprintf(stdout,"filename extension for input matrix: ");
  gets(filename_a);
  fprintf(stdout,"filename extension for large output matrix: ");
  gets(filename_cc);
  fprintf(stdout,"filename extension for small output matrix: ");
  gets(filename_c);

  strcpy( a.file_id, filename_a );

  /* load A matrix */ 

  strcpy(fstr0,BINDIR);
  strcat(fstr0,"n");
  strcat(fstr0,a.file_id);
  a.fid_n = fopen(fstr0,"rb");    
  if (a.fid_n == NULL) {
    fprintf(stderr, "%s\n",fstr0);
    fprintf(stderr, "Input file(s) not found\n");
    exit(2);
  }	    	    	    
  fprintf(stderr,"ROW_EXTRACT: reading %s ... ",fstr0);
  if (fread( a.n, sizeof (int), nc3+1, a.fid_n) != nc3+1){
    fprintf(stderr,"row_extract: problem reading file\n");
    exit(96);
  }
  fprintf(stderr,"read.\n");
  fclose(a.fid_n);
  
  a.iB = (int   *)malloc(a.n[nc3] * sizeof (int  ));
  a.vB = (float *)malloc(a.n[nc3] * sizeof (float));
  if ( a.iB == NULL || a.vB == NULL){
    fprintf(stderr,"malloc failed!\n");
    exit(3);
  }

  strcpy(fstr0,BINDIR);
  strcat(fstr0,"i");
  strcat(fstr0,a.file_id);
  a.fid_ind = fopen(fstr0,"rb");
  if (a.fid_ind == NULL) {
    fprintf(stderr, "%s\n",fstr0);
    fprintf(stderr, "Input file(s) not found\n");
    exit(2);
  }	    	    	    
  fprintf(stderr,"ROW_EXTRACT: reading %s ... ",fstr0);
  if (fread( a.iB, sizeof (int), a.n[nc3]
	     , a.fid_ind) != a.n[nc3]){
    fprintf(stderr,"ROW_EXTRACT: problem reading file\n");
    exit(96);
  }
  fprintf(stderr,"read.\n");
  fclose(a.fid_ind);

  strcpy(fstr0,BINDIR);
  strcat(fstr0,"v");
  strcat(fstr0,a.file_id);
  a.fid_val = fopen(fstr0,"rb");
  if (a.fid_val == NULL) {
    fprintf(stderr, "%s\n",fstr0);
    fprintf(stderr, "Input file(s) not found\n");
    exit(2);
  }	    	    	    
  fprintf(stderr,"ROW_EXTRACT: reading %s ... ",fstr0);
  if (fread( a.vB, sizeof (float), a.n[nc3]
	     , a.fid_val) != a.n[nc3]){
    fprintf(stderr,"ROW_EXTRACT: problem reading v file\n");
    exit(96);
  }
  fprintf(stderr,"read.\n");
  fclose(a.fid_val);

  j = 0;
  for (k = 0; k < a.n[nc3]; k++){
    if ( *(a.iB + k) > j  ) j = *(a.iB + k); 
  }
  a.nf = j + 1;

  a.y     = (float *)malloc( a.nf * sizeof (float) );
  a.delta = (float *)malloc( a.nf * sizeof (float) );
  if ( a.y == NULL || a.delta == NULL){
    fprintf(stderr,"malloc failed!\n");
    exit(3);
  }

  strcpy(fstr0,BINDIR);
  strcat(fstr0,"y");
  strcat(fstr0,a.file_id);
  a.fid_y = fopen(fstr0,"rb");    
  if (a.fid_y == NULL) {
    fprintf(stderr, "%s\n",fstr0);
    fprintf(stderr, "Input file(s) not found\n");
    exit(2);
  }	    	    	    
  fprintf(stderr,"ROW_EXTRACT: reading %s ... ",fstr0);
  if (fread( a.y, sizeof (float), a.nf, 
	     a.fid_y) != a.nf){
    fprintf(stderr,"row_extract: problem reading file\n");
    exit(96);
  }
  fprintf(stderr,"read.\n");
  fclose(a.fid_y);

  strcpy(fstr0,BINDIR);
  strcat(fstr0,"delta_");
  strcat(fstr0,a.file_id);
  a.fid_delta = fopen(fstr0,"rb");    
  if (a.fid_delta == NULL) {
    fprintf(stderr, "%s\n",fstr0);
    fprintf(stderr, "Input file(s) not found\n");
    exit(2);
  }	    	    	    
  fprintf(stderr,"ROW_EXTRACT: reading %s ... ",fstr0);
  if (fread( a.delta, sizeof (float), a.nf, 
	     a.fid_delta) != a.nf){
    fprintf(stderr,"row_extract: problem reading file\n");
    exit(96);
  }
  fprintf(stderr,"read.\n");
  fclose(a.fid_delta);


  /* put a into aa sparse structure */

  aa.nf = a.nf;
  aa.ncol = nc3;
  aa.n = a.n;
  aa.iB = a.iB;
  aa.vB = a.vB;

  /*  set up diag vectors */

  diag_d  = (int *)malloc(a.nf * sizeof (int));
  diag_dd = (int *)malloc(a.nf * sizeof (int));
  if ( diag_d == NULL || diag_dd == NULL){
    fprintf(stderr,"malloc failed!\n");
    exit(3);
  }

  for (k = 0; k < a.nf; k++){
    *(diag_d  + k) = 1;
  }
 
  srandom( (unsigned int) time(0) );
  for (j = 0; j < random_factor; j++){
    for (k = 0; k < a.nf; k++){
      *(diag_d  + k) *= random()%2;
    }
  }

  for (k = 0; k < a.nf; k++){
    *(diag_dd + k) = 1 - *(diag_d + k);
  }

 /* set up d and dd in sparse structure */

  d.ncol  =  a.nf; 
  dd.ncol =  a.nf;
  d.nf    =  a.nf;
  dd.nf   =  a.nf;

  d.n  = malloc( (d.ncol + 1) * sizeof (int));
  dd.n = malloc( (dd.ncol+ 1) * sizeof (int));
  if ( dd.n == NULL || d.n == NULL){
    fprintf(stderr,"malloc failed!\n");
    exit(3);
  }

  i = 0;
  j = 0;
  for (k = 0; k < d.ncol; k++){
    if ( *(diag_d + k) == 1) {
      i++;
    } else {
      j++;
    }
  }

  d.iB  = malloc( i * sizeof (int));
  dd.iB = malloc( j * sizeof (int));
  if ( dd.iB == NULL || d.iB == NULL){
    fprintf(stderr,"malloc failed!\n");
    exit(3);
  }

  d.vB  = malloc( i * sizeof (float));
  dd.vB = malloc( j * sizeof (float));
  if ( dd.vB == NULL || d.vB == NULL){
    fprintf(stderr,"malloc failed!\n");
    exit(3);
  }

  *dd.n = 0;
  *d.n  = 0;

  i = 0;
  j = 0;
  for (k = 0; k < d.ncol; k++){
    if ( *(diag_d + k) == 1) {
      *(dd.n + k + 1) = *(dd.n + k);
      *(d.n  + k + 1) = *(d.n  + k) + 1;
      *(d.iB + j) = k;
      *(d.vB + j) = 1.0;
      j++;
    } else {
      *(d.n  + k + 1) = *(d.n  + k);
      *(dd.n + k + 1) = *(dd.n + k) + 1;
      *(dd.iB + i) = k;
      *(dd.vB + i) = 1.0;
      i++;
    }
  }

  /* make space for y, delta vectors */

  y_c      = malloc( a.nf * sizeof (float));
  delta_c  = malloc( a.nf * sizeof(float));
  y_cc     = malloc( a.nf * sizeof (float));
  delta_cc = malloc( a.nf * sizeof(float));
  if ( y_c == NULL  || delta_c  == NULL ||
       y_cc == NULL || delta_cc == NULL   ){
    fprintf(stderr,"malloc failed!\n");
    exit(3);
  }

  /* multiply the matrices and write to disk */

  matmultSparse( d  , aa , filename_c  );
  matmultSparse( dd , aa , filename_cc );

  vmultSparse( d  , a.y     , y_c      , 1.0 );
  vmultSparse( dd , a.y     , y_cc     , 1.0 );
  vmultSparse( d  , a.delta , delta_c  , 1.0 );
  vmultSparse( dd , a.delta , delta_cc , 1.0 );

  strcpy(fstr0,BINDIR);
  strcat(fstr0,"y");
  strcat(fstr0,filename_c);
  fid = fopen(fstr0,"wb");    
  if (fid == NULL) {
    fprintf(stderr, "%s\n",fstr0);
    fprintf(stderr, "Unable to open file!\n");
    exit(2);
  }	    	    	    
  if ( fwrite( y_c, sizeof (int), a.nf, fid) != a.nf){
    fprintf(stderr,"row_extract: problem writing file\n%s\n",fstr0);
    exit(96);
  }
  fclose(fid);

  strcpy(fstr0,BINDIR);
  strcat(fstr0,"y");
  strcat(fstr0,filename_cc);
  fid = fopen(fstr0,"wb");    
  if (fid == NULL) {
    fprintf(stderr, "%s\n",fstr0);
    fprintf(stderr, "Unable to open file!\n");
    exit(2);
  }	    	    	    
  if ( fwrite( y_cc, sizeof (int), a.nf, fid) != a.nf){
    fprintf(stderr,"row_extract: problem writing file\n%s\n",fstr0);
    exit(96);
  }
  fclose(fid);

  strcpy(fstr0,BINDIR);
  strcat(fstr0,"delta_");
  strcat(fstr0,filename_c);
  fid = fopen(fstr0,"wb");    
  if (fid == NULL) {
    fprintf(stderr, "%s\n",fstr0);
    fprintf(stderr, "Unable to open file!\n");
    exit(2);
  }	    	    	    
  if ( fwrite( delta_c, sizeof (int), a.nf, fid) != a.nf){
    fprintf(stderr,"row_extract: problem writing file\n%s\n",fstr0);
    exit(96);
  }
  fclose(fid);

  strcpy(fstr0,BINDIR);
  strcat(fstr0,"delta_");
  strcat(fstr0,filename_cc);
  fid = fopen(fstr0,"wb");    
  if (fid == NULL) {
    fprintf(stderr, "%s\n",fstr0);
    fprintf(stderr, "Unable to open file!\n");
    exit(2);
  }	    	    	    
  if ( fwrite( delta_cc, sizeof (int), a.nf, fid) != a.nf){
    fprintf(stderr,"row_extract: problem writing file\n%s\n",fstr0);
    exit(96);
  }
  fclose(fid);

  free(a.iB); free(a.vB); free(a.y); free(a.delta);
  free(diag_d); free(diag_dd); 
  free(d.n);  free(d.iB);  free(d.vB);
  free(dd.n); free(dd.iB); free(dd.vB);
  free(y_c);  free(delta_c);
  free(y_cc); free(delta_cc);

 

  return(0);
}

