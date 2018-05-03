/*             row_to_col.c by Richard Frazin 1999
 * This takes the old semi-sparse row format (e.g., the output of 
 *   derivs_hollowsph.m) and converts it into column sparse format.
 *
 * to compile just type "make row_to_col"
 */


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include "headers.h"


void usage(char *command_name) {
  printf("usage: %s <directory name> <file name> <rows> <columns>\n", command_name);
}

int main(int argc, char **argv) {
  FILE *fid_nA;
  int numcols, numrows;
  char bindir[MAXPATH], filestring[MAXPATH];
  char nfilestring[MAXPATH], ifilestring[MAXPATH], vfilestring[MAXPATH];
  char fstr1[MAXPATH], fstr2[MAXPATH];
  int *q, *count;
  int *n_colmat, *n_rowmat;
  int *i_colmat, *i_rowmat;
  float *v_colmat, *v_rowmat;
  int j, k, nnz, place, col;
  char optstring[] = "h";
  char opt;


  while((opt = getopt(argc, argv, optstring)) != -1) {
    switch(opt) {
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
    fprintf(stdout, "WARNING: THIS OVERWRITES ORIGINAL FILES!\n");
    fprintf(stdout, "directory name (include final slash) ");
    scanf("%s", bindir);
    fprintf(stdout, "filename extension: ");
    scanf("%s", fstr1);
    fprintf(stdout, "How many rows in (full) matrix? ");
    scanf("%s", fstr2);
    sscanf(fstr2, "%d", &numrows);
    fprintf(stdout, "How many columns in (full) matrix? ");
    scanf("%s", fstr2);
    sscanf(fstr2, "%d", &numcols);
  } else if (argc - optind == 4) {
    strcpy(bindir, argv[optind]);
    strcpy(fstr1, argv[optind+1]);
    numrows = strtol(argv[optind+2], NULL, 10);
    numcols = strtol(argv[optind+3], NULL, 10);
  } else {
    usage(argv[0]);
    return(1);
  }

  n_colmat = (int *) malloc((numcols + 1) * sizeof(int));
  n_rowmat = (int *) malloc((numrows + 1) * sizeof(int));
  if (n_colmat == NULL || n_rowmat == NULL) {
    fprintf(stderr,
            "ROW_TO_COL: malloc error (n_rowmat or n_colmat)!");
    exit(2);
  }


  strcpy(nfilestring, bindir);
  strcat(nfilestring, "n");
  strcat(nfilestring, fstr1);

  strcpy(ifilestring, bindir);
  strcat(ifilestring, "i");
  strcat(ifilestring, fstr1);

  strcpy(vfilestring, bindir);
  strcat(vfilestring, "v");
  strcat(vfilestring, fstr1);


  strcpy(filestring, nfilestring);
  fid_nA = fopen(filestring, "rb");
  if (fid_nA == NULL) {
    fprintf(stderr, "Cannot open input file! Filename:\n");
    fprintf(stderr, "%s\n", filestring);
    exit(2);
  }
  fprintf(stderr, "ROW_TO_COL: reading %s ... ", filestring);
  if (fread(n_rowmat, sizeof(int), numrows + 1, fid_nA) != numrows + 1) {
    fprintf(stderr, "row_to_col: problem reading file\n");
    exit(96);
  }
  fclose(fid_nA);
  fprintf(stderr, "read.\n");


  nnz = n_rowmat[numrows];


  i_rowmat = (int *) malloc(nnz * sizeof(int));
  i_colmat = (int *) malloc(nnz * sizeof(int));
  v_rowmat = (float *) malloc(nnz * sizeof(float));
  v_colmat = (float *) malloc(nnz * sizeof(float));
  if (i_colmat == NULL || i_rowmat == NULL ||
      v_colmat == NULL || v_rowmat == NULL) {
    fprintf(stderr,
            "ROW_TO_COL: malloc error ((n,v)_rowmat or (n,v)_colmat)!");
    exit(2);
  }


  strcpy(filestring, ifilestring);
  fid_nA = fopen(filestring, "rb");
  if (fid_nA == NULL) {
    fprintf(stderr, "Cannot open input file! Filename:\n");
    fprintf(stderr, "%s\n", filestring);
    exit(2);
  }
  fprintf(stderr, ": reading %s ... ", filestring);
  if (fread(i_rowmat, sizeof(int), nnz, fid_nA) != nnz) {
    fprintf(stderr, "ROW_TO_COL: problem reading file\n");
    exit(96);
  }
  fclose(fid_nA);
  fprintf(stderr, "read.\n");


  strcpy(filestring, vfilestring);
  fid_nA = fopen(filestring, "rb");
  if (fid_nA == NULL) {
    fprintf(stderr, "Cannot open input file! Filename:\n");
    fprintf(stderr, "%s\n", filestring);
    exit(2);
  }
  fprintf(stderr, "ROW_TO_COL: reading %s ... ", filestring);
  if (fread(v_rowmat, sizeof(float), nnz, fid_nA) != nnz) {
    fprintf(stderr, "row_to_col: problem reading file\n");
    exit(96);
  }
  fclose(fid_nA);
  fprintf(stderr, "read.\n");


  /* print input for debugging */

  /*
     fprintf(stdout,"nnz = %d\n",nnz);
     fprintf(stdout,"n: ");
     for (j = 0; j <= numrows; j++)
     fprintf(stdout,"%d ",n_rowmat[j]);
     fprintf(stdout,"\ni: ");
     for (j = 0; j < nnz; j++)
     fprintf(stdout,"%d ",i_rowmat[j]);
     fprintf(stdout,"\nv: ");
     for (j = 0; j < nnz; j++)
     fprintf(stdout,"%g ",v_rowmat[j]);
     fprintf(stdout,"\n\n");
   */

  /*build q vector which contains the number of elements in each column
   *      count is a handy indexing vector */

  q = (int *) malloc(numcols * sizeof(int));
  count = (int *) malloc(numcols * sizeof(int));
  if (q == NULL) {
    fprintf(stderr, "ROW_TO_COL: malloc error (q or count)!");
    exit(2);
  }

  for (j = 0; j < numcols; j++) {
    q[j] = 0;
    count[j] = 0;
  }

  for (j = 0; j < nnz; j++)
    q[i_rowmat[j]]++;

  /* fill out n_colmat vector */

  n_colmat[0] = 0;
  for (j = 1; j <= numcols; j++)
    n_colmat[j] = n_colmat[j - 1] + q[j - 1];


  /*fill out i_colmat and v_colmat vectors */

  for (j = 0; j < numrows; j++) {
    for (k = n_rowmat[j]; k < n_rowmat[j + 1]; k++) {

      col = i_rowmat[k];
      place = n_colmat[col] + count[col];

      i_colmat[place] = j;
      v_colmat[place] = v_rowmat[k];

      count[col]++;

    }
  }

  /* check for consistency */
  for (j = 0; j < numcols; j++) {
    if (count[j] != q[j]) {
      fprintf(stderr, "count != q, column = %d\n", j);
      exit(3);
    }
  }

  /* print output for debugging */

  /*
     fprintf(stdout,"\nq: ");
     for (j = 0; j < numcols; j++)
     fprintf(stdout,"%d ",q[j]);
     fprintf(stdout,"\nn: ");
     for (j = 0; j <= numcols; j++)
     fprintf(stdout,"%d ",n_colmat[j]);
     fprintf(stdout,"\ni: ");
     for (j = 0; j < nnz; j++)
     fprintf(stdout,"%d ",i_colmat[j]);
     fprintf(stdout,"\nv: ");
     for (j = 0; j < nnz; j++)
     fprintf(stdout,"%g ",v_colmat[j]);
     fprintf(stdout,"\n\n");
   */

  /* write out column formatted data */

  strcpy(filestring, nfilestring);
  fid_nA = fopen(filestring, "wb");
  if (fid_nA == NULL) {
    fprintf(stderr, "Cannot open output file! Filename:\n");
    fprintf(stderr, "%s\n", filestring);
    exit(2);
  }
  fprintf(stderr, "ROW_TO_COL: writing %s ... ", filestring);
  if (fwrite(n_colmat, sizeof(int), numcols + 1, fid_nA) != numcols + 1) {
    fprintf(stderr, "row_to_col: problem reading file\n");
    exit(96);
  }
  fclose(fid_nA);
  fprintf(stderr, "written.\n");

  strcpy(filestring, ifilestring);
  fid_nA = fopen(filestring, "wb");
  if (fid_nA == NULL) {
    fprintf(stderr, "Cannot open output file! Filename:\n");
    fprintf(stderr, "%s\n", filestring);
    exit(2);
  }
  fprintf(stderr, "ROW_TO_COL: writing %s ... ", filestring);
  if (fwrite(i_colmat, sizeof(int), nnz, fid_nA) != nnz) {
    fprintf(stderr, "row_to_col: problem reading file\n");
    exit(96);
  }
  fclose(fid_nA);
  fprintf(stderr, "written.\n");

  strcpy(filestring, vfilestring);
  fid_nA = fopen(filestring, "wb");
  if (fid_nA == NULL) {
    fprintf(stderr, "Cannot open output file! Filename:\n");
    fprintf(stderr, "%s\n", filestring);
    exit(2);
  }
  fprintf(stderr, "ROW_TO_COL: writing %s ... ", filestring);
  if (fwrite(v_colmat, sizeof(float), nnz, fid_nA) != nnz) {
    fprintf(stderr, "row_to_col: problem reading file\n");
    exit(96);
  }
  fclose(fid_nA);
  fprintf(stderr, "written.\n");


  return (0);
}
