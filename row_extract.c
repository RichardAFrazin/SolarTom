#include <malloc.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <assert.h>

#include "headers.h"

/* this code calculates for the matrix a
 * c  =  d * a    and
 * cc = dd * a .
 * d is diagonal, containing only zeros or ones.  
 * dd = I - d.
 *
 * d is such that rows rmin to rmax will be extracted
 *
 * c and cc are written to disk. 
 * cc contains most of a.
 * this also writes the delta and y vectors to disk
 */

int main(int argc, char **argv) {
  char filename_a[MAXPATH], filename_d[MAXPATH], filename_dd[MAXPATH];
  struct matrix a;
  struct sparse aa;
  char fstr0[MAXPATH];
  int j, k;
  int nc3 = NBINS;
  int rmin, rmax;

  printf("rmin: ");
  scanf("%d", &rmin);
  printf("rmax: ");
  scanf("%d", &rmax);
  fprintf(stdout, "filename extension for input matrix: ");
  scanf("%s", filename_a);
  fprintf(stdout, "filename extension for d matrix: ");
  scanf("%s", filename_d);
  fprintf(stdout, "filename extension for dd matrix: ");
  scanf("%s", filename_dd);
  
  strcpy(a.file_id, filename_a);

  /* load A matrix */

  strcpy(fstr0, BINDIR);
  strcat(fstr0, "n");
  strcat(fstr0, a.file_id);
  a.fid_n = fopen(fstr0, "rb");
  if (a.fid_n == NULL) {
    fprintf(stderr, "%s\n", fstr0);
    fprintf(stderr, "Input file(s) not found\n");
    exit(2);
  }
  fprintf(stderr, "row_extract: reading %s ... ", fstr0);
  if (fread(a.n, sizeof(int), nc3 + 1, a.fid_n) != nc3 + 1) {
    fprintf(stderr, "row_extract: problem reading file\n");
    exit(96);
  }
  fprintf(stderr, "read.\n");
  fclose(a.fid_n);

  a.iB = (int *) malloc(a.n[nc3] * sizeof(int));
  a.vB = (float *) malloc(a.n[nc3] * sizeof(float));
  if (a.iB == NULL || a.vB == NULL) {
    fprintf(stderr, "malloc failed!\n");
    exit(3);
  }

  strcpy(fstr0, BINDIR);
  strcat(fstr0, "i");
  strcat(fstr0, a.file_id);
  a.fid_ind = fopen(fstr0, "rb");
  if (a.fid_ind == NULL) {
    fprintf(stderr, "%s\n", fstr0);
    fprintf(stderr, "Input file(s) not found\n");
    exit(2);
  }
  fprintf(stderr, "row_extract: reading %s ... ", fstr0);
  if (fread(a.iB, sizeof(int), a.n[nc3], a.fid_ind) != a.n[nc3]) {
    fprintf(stderr, "row_extract: problem reading file\n");
    exit(96);
  }
  fprintf(stderr, "read.\n");
  fclose(a.fid_ind);

  strcpy(fstr0, BINDIR);
  strcat(fstr0, "v");
  strcat(fstr0, a.file_id);
  a.fid_val = fopen(fstr0, "rb");
  if (a.fid_val == NULL) {
    fprintf(stderr, "%s\n", fstr0);
    fprintf(stderr, "Input file(s) not found\n");
    exit(2);
  }
  fprintf(stderr, "row_extract: reading %s ... ", fstr0);
  if (fread(a.vB, sizeof(float), a.n[nc3], a.fid_val) != a.n[nc3]) {
    fprintf(stderr, "row_extract: problem reading v file\n");
    exit(96);
  }
  fprintf(stderr, "read.\n");
  fclose(a.fid_val);

  j = 0;
  for (k = 0; k < a.n[nc3]; k++) {
    if (*(a.iB + k) > j)
      j = *(a.iB + k);
  }
  a.nf = j + 1;

  a.y = (float *) malloc(a.nf * sizeof(float));
  assert(a.y != NULL);

  a.delta = (float *) malloc(a.nf * sizeof(float));
  assert(a.delta != NULL);
  
  strcpy(fstr0, BINDIR);
  strcat(fstr0, "y");
  strcat(fstr0, a.file_id);
  a.fid_y = fopen(fstr0, "rb");
  if (a.fid_y == NULL) {
    fprintf(stderr, "%s\n", fstr0);
    fprintf(stderr, "Input file(s) not found\n");
    exit(2);
  }
  fprintf(stderr, "row_extract: reading %s ... ", fstr0);
  if (fread(a.y, sizeof(float), a.nf, a.fid_y) != a.nf) {
    fprintf(stderr, "row_extract: problem reading file\n");
    exit(96);
  }
  fprintf(stderr, "read.\n");
  fclose(a.fid_y);

  strcpy(fstr0, BINDIR);
  strcat(fstr0, "delta_");
  strcat(fstr0, a.file_id);
  a.fid_delta = fopen(fstr0, "rb");
  if (a.fid_delta == NULL) {
    fprintf(stderr, "%s\n", fstr0);
    fprintf(stderr, "Input file(s) not found\n");
    exit(2);
  }
  fprintf(stderr, "row_extract: reading %s ... ", fstr0);
  if (fread(a.delta, sizeof(float), a.nf, a.fid_delta) != a.nf) {
    fprintf(stderr, "row_extract: problem reading file\n");
    exit(96);
  }
  fprintf(stderr, "read.\n");
  fclose(a.fid_delta);


  /* put a into aa sparse structure */
  aa.nf = a.nf;
  aa.ncol = nc3;
  aa.n = a.n;
  aa.iB = a.iB;
  aa.vB = a.vB;

  /* multiply the matrices and write to disk */
  fprintf(stderr, "row_extract: multiplying matricies\n");

  row_extract(aa, a.y, a.delta, filename_d, filename_dd, rmin, rmax);
  
  free(a.iB);
  free(a.vB);
  free(a.y);
  free(a.delta);

  return (0);
}
