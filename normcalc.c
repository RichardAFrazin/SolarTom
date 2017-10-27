/* normcalc.c calculates || Ax - y ||^2  */
/*
 * nf is the number rows in A
 *
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "headers.h"

float normcalc(fstr2, x_infile)
char *fstr2, *x_infile;
{
  FILE *fid_x;
  float x[NBINS];
  float norm, tmp;
  struct matrix A;
  struct sparse sparmat;
  char filestring[MAXPATH];
  int i, j, k, nc3;
  /* char fstr2[MAXPATH], x_infile[MAXPATH]; */

  nc3 = NBINS;

  /*
     fprintf(stdout,"input filename for x vector: ");
     gets(x_infile);
     fprintf(stdout,"input filename extension: ");
     gets(fstr2);
   */

  /* load stuff and malloc more stuff */

  strcpy(filestring, BINDIR);
  strncat(filestring, x_infile, MAXPATH);
  fid_x = fopen(filestring, "rb");
  if (fid_x == NULL) {
    fprintf(stderr, "Input file(s) not found\n");
    fprintf(stderr, "%s\n", filestring);
    exit(2);
  }
  if (fread(x, sizeof(float), nc3, fid_x) != nc3)
    fprintf(stderr, "error reading file: %s\n", filestring);
  fclose(fid_x);


  strcpy(filestring, BINDIR);
  strcat(filestring, "n");
  strcat(filestring, fstr2);
  A.fid_n = fopen(filestring, "rb");
  if (A.fid_n == NULL) {
    fprintf(stderr, "file not found! file = %s\n", filestring);
    exit(1);
  }

  if (fread(A.n, sizeof(int), nc3 + 1, A.fid_n) != (nc3 + 1)) {
    fprintf(stderr, "problem with fread! file = %s\n", filestring);
    exit(2);
  }
  fclose(A.fid_n);

  A.iB = (int *) malloc(A.n[nc3] * sizeof(int));
  A.vB = (float *) malloc(A.n[nc3] * sizeof(float));
  if (A.iB == NULL || A.vB == NULL) {
    fprintf(stderr, "malloc error in (iB and/or vB)!\n");
    exit(1);
  }

  strcpy(filestring, BINDIR);
  strcat(filestring, "i");
  strcat(filestring, fstr2);
  A.fid_ind = fopen(filestring, "rb");
  if (A.fid_ind == NULL) {
    fprintf(stderr, "Input file(s) not found\n");
    fprintf(stderr, "%s\n", filestring);
    exit(2);
  }
  if (fread(A.iB, sizeof(int), A.n[nc3], A.fid_ind) != A.n[nc3]) {
    fprintf(stderr, "problem reading! file = %s\n", filestring);
    exit(96);
  }
  fclose(A.fid_ind);

  strcpy(filestring, BINDIR);
  strcat(filestring, "v");
  strcat(filestring, fstr2);
  A.fid_val = fopen(filestring, "rb");
  if (A.fid_val == NULL) {
    fprintf(stderr, "Input file(s) not found\n");
    fprintf(stderr, "%s\n", filestring);
    exit(2);
  }
  if (fread(A.vB, sizeof(float), A.n[nc3], A.fid_val) != A.n[nc3]) {
    fprintf(stderr, "problem reading! file = %s\n", filestring);
    exit(96);
  }
  fclose(A.fid_val);

  j = 0;
  for (k = 0; k < A.n[nc3]; k++) {
    if (*(A.iB + k) > j)
      j = *(A.iB + k);
  }
  A.nf = j + 1;

  A.huber = 0;
  A.r = (float *) malloc(A.nf * sizeof(float));
  A.y = (float *) malloc(A.nf * sizeof(float));
  A.delta = NULL;
  if (A.r == NULL || A.y == NULL) {
    fprintf(stderr, "malloc error in A.r or A.y! \n");
    exit(0);
  }

  strcpy(filestring, BINDIR);
  strcat(filestring, "y");
  strcat(filestring, fstr2);
  A.fid_y = fopen(filestring, "rb");
  if (A.fid_y == NULL) {
    fprintf(stderr, "file not found! file = %s\n", filestring);
    exit(1);
  }

  if (fread(A.y, sizeof(float), A.nf, A.fid_y) != A.nf) {
    fprintf(stderr, "problem with fread! file = %s\n", filestring);
    exit(2);
  }
  fclose(A.fid_y);

  /* set up sparmat */

  sparmat.ncol = nc3;
  sparmat.nf = A.nf;
  sparmat.n = A.n;
  sparmat.iB = A.iB;
  sparmat.vB = A.vB;

  /* set A.r = Ax */

  for (i = 0; i < A.nf; i++)
    A.r[i] = 0.0;

  vmultSparse(sparmat, x, A.r, 1.0);

  /* calculate norm */

  norm = 0.0;

  for (i = 0; i < A.nf; i++) {
    tmp = A.r[i] - A.y[i];
    norm += tmp * tmp;
  }

  free(A.iB);
  free(A.vB);
  free(A.r);
  free(A.y);


  /* fprintf(stdout,"norm = %g\n",norm); */

  return (norm);

}
