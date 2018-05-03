/*       call_rcs_2_ccs.c by Richard Frazin May 2008 
 *  This converts row sparse format files to column sparse files.  This is the 
 *   (w,j,m) to (v,i,n) conversion.
 *
 *    Compile instructions:
 * gcc -O3 -o call_rcs_2_ccs call_rcs_2_ccs.c rcs_llist.c llist.c -lm
 *
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <assert.h>
#include <getopt.h>
#include "headers.h"

static void usage(const char *arg0);
/* matrix format converter routine */
static void convert_rcs_2_ccs(const char *rw, const char *rj, const char *rm,
   const char *cv, const char *ci, const char *cn, int n_cols, int n_rows);
static void usage(const char *arg0) {
  printf("usage: %s <-h> <rcs_filename_suffix> <n_col> \n", arg0);
}


int main(int argc, char **argv){
  char rw[MAXPATH], rj[MAXPATH], rm[MAXPATH], cv[MAXPATH], ci[MAXPATH], cn[MAXPATH];
  char fn[MAXPATH];
  int n_rows, n_cols, length, k;
  FILE *fid_j, *fid_m;
  int opt;
  char optstring[] = "h"; /* command line options*/


  fprintf(stdout,"For some reason this doesn't work and trips up inside convert_rcs_2_ccs ... compare to builda.c, which works.");
  assert(0);


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
  if (argc != 3){
    usage(argv[0]);
    return(0);
  }

  fprintf(stdout,"filename suffix is: %s\n",argv[1]);
  n_cols = atoi(argv[2]);

  strcpy(rw,BINDIR);
  strcat(rw,"w");
  strcat(rw,argv[1]);
  strcpy(rj,BINDIR);
  strcat(rj,"j");
  strcat(rj,argv[1]);
  strcpy(rm,BINDIR);
  strcat(rm,"m");
  strcat(rm,argv[1]);
  strcpy(cv,BINDIR);
  strcat(cv,"v");
  strcat(cv,argv[1]);
  strcpy(ci,BINDIR);
  strcat(ci,"i");
  strcat(ci,argv[1]);
  strcpy(cn,BINDIR);
  strcat(cn,"n");
  strcat(cn,argv[1]);

  fid_j = fopen(rj,"rb");
  if (fid_j == NULL){
    printf("can't open file: %s\n",rj);
    return(1);
  }
  fid_m = fopen(rm,"rb");
  if (fid_m == NULL){
    printf("can't open file: %s\n",rm);
    return(1);
  }

  /* The number of rows is the length of the m file - 1. 
   * The last number in the m file is the size 
   *   of the j and w files. */
  assert(1 + fseek(fid_m, -1*sizeof(int), SEEK_END));
  fread(&length, sizeof(int), 1, fid_m);
  /*
  printf("lastnumber = %d, filelength = %d\n",length, ftell(fid_m)/sizeof(int));
  */ 
  n_rows = ftell(fid_m)/sizeof (int) - 1;
  fclose(fid_m);

  /*scan the j file to get the number of columns.
   * This is disabled because the last P columns may not be 
   * present in the j file, but they need to be accounted
   * for in the n file.
   */

/*
  ja = (int *) malloc( length * sizeof(int) );
  if (ja == NULL){
    printf("malloc error: ja");
    exit(1);
  }
  assert(fread(ja, sizeof(int), length, fid_j) == length);

  n_cols = ja[0];
  for (k = 1; k < length; k++){
    if (ja[k] > n_cols)
      n_cols = ja[k];
  }
  free(ja);
*/
  fclose(fid_j);
  
  fprintf(stdout,"n_rows = %d, n_cols = %d\n",n_rows, n_cols); fflush(stdout);
  fprintf(stderr,"\nStarting row sparse to col sparse conversion...");

  convert_rcs_2_ccs(rw, rj, rm, cv, ci, cn, n_cols, n_rows);

  fprintf(stderr,"..done\n");

  return(0);
}

/* this code is taken from builda.c */
static void convert_rcs_2_ccs(const char *rcs_v,
			      const char *rcs_j,
                              const char *rcs_r,
                              const char *ccs_v,
                              const char *ccs_i,
                              const char *ccs_n,
                              int n_cols, int n_rows) {
  FILE *rcs_v_fid, *rcs_j_fid, *rcs_r_fid;
  FILE *ccs_v_fid, *ccs_i_fid, *ccs_n_fid;

  int *n;
  int col, k, row_start, next_row_start, n_elem, row_index, col_index, count;
  float val;
  
  float **v;
  int **i;
  
  assert(rcs_v && rcs_j && rcs_r);
  assert(ccs_v && ccs_i && ccs_n);
  assert(n_cols > 0);
  assert(n_rows > 0);

  assert(rcs_v_fid = fopen(rcs_v, "r"));
  assert(rcs_j_fid = fopen(rcs_j, "r"));
  assert(rcs_r_fid = fopen(rcs_r, "r"));

  n = (int *) calloc(n_cols, sizeof(int));
  assert(n);

  /* Determine the # of elements in each column */
  while(fread(&col, sizeof(int), 1, rcs_j_fid) == 1) {
    assert(col >= 0 && col < n_cols);
    n[col]++;
  }
  rewind(rcs_j_fid);


  /* Allocate memory */
  v = (float **) malloc(n_cols * sizeof(float *));
  assert(v);
  
  i = (int   **) malloc(n_cols * sizeof(int *));
  assert(i);
  
  for (col_index = 0; col_index < n_cols; col_index++) {
    /* Assertion fails if there are no lines of sight that pass
       through the col_index voxel */
    /*assert(n[col_index] > 0);*/

    if (n[col_index] > 0) {
      v[col_index] = (float *) malloc(n[col_index] * sizeof(float));
      assert(v[col_index]);
      
      i[col_index] = (int *)   malloc(n[col_index] * sizeof(int));
      assert(i[col_index]);
    }
  }

  /* Clear n - it now stores the index to write to for the column
     format in memory */
  memset(n, 0, n_cols * sizeof(int));
  
  /* Store rcs elements columnwise in memory */
  assert(fread(&row_start, sizeof(int), 1, rcs_r_fid) == 1);
  assert(row_start == 0);
  
  for (row_index = 0; row_index < n_rows; row_index++) {
    assert(fread(&next_row_start, sizeof(int), 1, rcs_r_fid) == 1);
    /* Assertion fails if there is an empty row */
    /*fprintf(stdout,"(%d %d)",next_row_start,row_start);fflush(stdout);*/
    assert(next_row_start > row_start);
    
    n_elem = next_row_start - row_start;

    for (k = 0; k < n_elem; k++) {
      assert(fread(&val, sizeof(float), 1, rcs_v_fid) == 1);
      assert(fread(&col, sizeof(int),   1, rcs_j_fid) == 1);

      v[col][n[col]] = val;
      i[col][n[col]] = row_index;

      n[col]++;
    }
    
    row_start = next_row_start;
  }

  fclose(rcs_v_fid);
  fclose(rcs_j_fid);
  fclose(rcs_r_fid);
  
  assert(ccs_v_fid = fopen(ccs_v, "w"));
  assert(ccs_i_fid = fopen(ccs_i, "w"));
  assert(ccs_n_fid = fopen(ccs_n, "w"));

  /* Write columnwise matrix from memory to file */
  for (col_index = 0; col_index < n_cols; col_index++) {
    assert(fwrite(v[col_index], sizeof(float), n[col_index], ccs_v_fid) ==
	   n[col_index]);

    assert(fwrite(i[col_index], sizeof(int), n[col_index], ccs_i_fid) ==
	   n[col_index]);
  }

  count = 0;
  assert(fwrite(&count, sizeof(int), 1, ccs_n_fid) == 1);
  
  for (col_index = 0; col_index < n_cols; col_index++) {
    count += n[col_index];

    assert(fwrite(&count, sizeof(int), 1, ccs_n_fid) == 1);
  }
  
  fclose(ccs_v_fid);
  fclose(ccs_i_fid);
  fclose(ccs_n_fid);

  /* Cleanup memory */
  for (col_index = 0; col_index < n_cols; col_index++) {
    if (n[col_index] > 0) {
      free(v[col_index]);
      free(i[col_index]);
    }
  }

  free(v);
  free(i);

  free(n);
}
