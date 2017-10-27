/* this code calculates for the matrix A and vector x
 * y  =  A * x    and
 *
 *   and writes the new y vector to disk.
 *   the new y vector replaces the old one in memory
 */


/*#include <malloc.h>*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>

#include "headers.h"


void usage(char *arg0) {
  printf("usage: %s <-h> [<matrix_name> <y_name> <x_name>]\n", arg0);
}

int main(int argc, char **argv) {
  char filename_a[MAXPATH], filename_x[MAXPATH], filename_y[MAXPATH];
  FILE *fid;
  struct matrix a;
  struct sparse aa;
  float x[NBINS];
  float *y_new;
  char fstr0[MAXPATH];
  int j, k;
  int nc3 = NBINS;
  
  int opt;
  char optstring[] = "h";

  /* process command line arguments */
  while ((opt = getopt(argc, argv, optstring)) != -1) {
    switch (opt) {
    case 'h':
      usage(argv[0]);
      return 0;
      break;
    default:
      usage(argv[0]);
      return 1;
    }
  }

  if (argc - optind == 0) {
    fprintf(stdout, "filename extension for input matrix: ");
    scanf("%s", filename_a);
    fprintf(stdout, "filename extension for new y vector: ");
    scanf("%s", filename_y);
    fprintf(stdout, "filename of input x vector: ");
    scanf("%s", filename_x);
    
    strcpy(a.file_id, filename_a);
  }
  else if (argc - optind == 3) {
    strcpy(filename_a, argv[optind]);
    strcpy(filename_y, argv[optind+1]);
    strcpy(filename_x, argv[optind+2]);

    strcpy(a.file_id, filename_a);
  }
  else {
    usage(argv[0]);
    return 1;
  }

  /* load x vector */

  strcpy(fstr0, BINDIR);
  strcat(fstr0, filename_x);
  fid = fopen(fstr0, "rb");
  if (fid == NULL) {
    fprintf(stderr, "%s\n", fstr0);
    fprintf(stderr, "Input file(s) not found\n");
    exit(2);
  }
  fprintf(stderr, "CALCULATE_Y: reading %s ... ", fstr0);
  if (fread(x, sizeof(float), nc3, fid) != nc3) {
    fprintf(stderr, "CALCULATE_Y: problem reading file\n");
    exit(96);
  }
  fprintf(stderr, "read.\n");
  fclose(fid);

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
  fprintf(stderr, "CALCULATE_Y: reading %s ... ", fstr0);
  if (fread(a.n, sizeof(int), nc3 + 1, a.fid_n) != nc3 + 1) {
    fprintf(stderr, "CALCULATE_Y: problem reading file\n");
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
  fprintf(stderr, "CALCULATE_Y: reading %s ... ", fstr0);
  if (fread(a.iB, sizeof(int), a.n[nc3], a.fid_ind) != a.n[nc3]) {
    fprintf(stderr, "CALCULATE_Y: problem reading file\n");
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
  fprintf(stderr, "CALCULATE_Y: reading %s ... ", fstr0);
  if (fread(a.vB, sizeof(float), a.n[nc3], a.fid_val) != a.n[nc3]) {
    fprintf(stderr, "CALCULATE_Y: problem reading v file\n");
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
  a.delta = (float *) malloc(a.nf * sizeof(float));
  if (a.y == NULL || a.delta == NULL) {
    fprintf(stderr, "malloc failed!\n");
    exit(3);
  }

  strcpy(fstr0, BINDIR);
  strcat(fstr0, "y");
  strcat(fstr0, a.file_id);
  a.fid_y = fopen(fstr0, "rb");
  if (a.fid_y == NULL) {
    fprintf(stderr, "%s\n", fstr0);
    fprintf(stderr, "Input file(s) not found\n");
    exit(2);
  }
  fprintf(stderr, "CALCULATE_Y: reading %s ... ", fstr0);
  if (fread(a.y, sizeof(float), a.nf, a.fid_y) != a.nf) {
    fprintf(stderr, "CALCULATE_Y: problem reading file\n");
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
  fprintf(stderr, "CALCULATE_Y: reading %s ... ", fstr0);
  if (fread(a.delta, sizeof(float), a.nf, a.fid_delta) != a.nf) {
    fprintf(stderr, "CALCULATE_Y: problem reading file\n");
    exit(96);
  }
  fprintf(stderr, "read.\n");
  fclose(a.fid_delta);

  /* assign y pointer */

  y_new = a.y;


  /* put a into aa sparse structure */

  aa.nf = a.nf;
  aa.ncol = nc3;
  aa.n = a.n;
  aa.iB = a.iB;
  aa.vB = a.vB;



  /* multiply and write to disk */

  vmultSparse(aa, x, y_new, 1.0);


  strcpy(fstr0, BINDIR);
  strcat(fstr0, "y");
  strcat(fstr0, filename_y);
  fid = fopen(fstr0, "wb");
  if (fid == NULL) {
    fprintf(stderr, "%s\n", fstr0);
    fprintf(stderr, "Unable to open file!\n");
    exit(2);
  }
  if (fwrite(y_new, sizeof(float), a.nf, fid) != a.nf) {
    fprintf(stderr, "CALCULATE_Y: problem writing file\n%s\n", fstr0);
    exit(96);
  }
  fclose(fid);


  return (0);
}
