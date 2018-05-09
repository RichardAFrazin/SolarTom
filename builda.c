/*#include <malloc.h>*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <assert.h>
#include "headers.h"

static void usage(const char *arg0);
/* matrix format converter routine */
static void convert_rcs_2_ccs(const char *rcs_v,
                              const char *rcs_j,
                              const char *rcs_r,
                              const char *ccs_v,
                              const char *ccs_i,
                              const char *ccs_n,
			      int n_cols, int n_rows);
  
static void usage(const char *arg0) {
  printf("usage: %s <-h> [<config file> <outfile suffix>]\n", arg0);
}

int main(int argc, char **argv) {
  char confstring[] = CONFSTRING;
  char a_outfile[MAXPATH];
  FILE *fid_y, *fid_delta, *fid_w, *fid_j, *fid_m;
  FILE *fid_conf, *fid_log, *fid_date, *fid_block;
  int nfiles, i, er;
  char idstring[MAXPATH], filestring[MAXPATH];
  char filename_v[MAXPATH], filename_i[MAXPATH], filename_n[MAXPATH];
  char fn_w[MAXPATH], fn_m[MAXPATH], fn_j[MAXPATH];
  int nc3, n_elem_exported;
  float *yy, *dd;
  int len_yy, len_y, opt;
  rcs_llist *rcs;
  char optstring[] = "h";

  strcpy(a_outfile, A_OUTFILE);

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
    /* nothing special for this case, we proceed as before */
    ;
  } else if (argc - optind == 2) {
    strcpy(confstring, argv[optind]);
    strcpy(a_outfile, argv[optind+1]);
  } else {
    usage(argv[0]);
    return 1;
  }

  nc3 = NBINS;
  len_y = 0;

  /* get the string containing the time and remove the spaces 
   *    the resulting string time_string will be used to make
   *    unique scratch file names. 
   *
   *      to do this, one needs the variable declaration:
   *   time_t current_time;  and  #include <time.h>
   *
   * This is no longer needed since the row sparse format
   *   files (w,j,m) are treated just like the column format
   *   files (v,i,n).  
   */

  /*
  current_time = time(NULL);
  time_str = ctime(&current_time);
  k = 0; i = 0;
  while ( *(time_str + k) != '\0'){
    if (*(time_str + k) != 32){
      strncpy(time_string + i , time_str + k,1);
      i++;
    }
    k++;
  }
  sprintf(time_string + i-1,"%s","\0");
  fprintf(stdout,"The current time is: %s\n",time_str);
  fflush(stdout);
  */

  /* open output files */
  strcpy(filestring, BINDIR);
  strcat(filestring, "log_");
  strcat(filestring, a_outfile);
  fid_log = fopen(filestring, "wb");

  strcpy(filestring, BINDIR);  /* starting row number of subA */
  strcat(filestring, "block_");
  strcat(filestring, a_outfile);
  fid_block = fopen(filestring, "wb");
  
  strcpy(filestring, BINDIR);  /* modified Julian date of subA */
  strcat(filestring, "date_");
  strcat(filestring, a_outfile);
  fid_date = fopen(filestring, "wb");
  
  strcpy(filestring, BINDIR);
  strcat(filestring, "y");
  strcat(filestring, a_outfile);
  fid_y = fopen(filestring, "wb");

  strcpy(filestring, BINDIR);
  strcat(filestring, "delta_");
  strcat(filestring, a_outfile);
  fid_delta = fopen(filestring, "wb");


  if (fid_y    == NULL || fid_delta == NULL || 
      fid_block == NULL || fid_log   == NULL ||
      fid_date == NULL) {
    fprintf(stderr, "build_A: Cannot open '%s' \n", filestring);
    exit(2);
  }

  if ((fid_conf = fopen(confstring, "r")) == NULL) {
    fprintf(stderr, "Main: Input file '%s' not found\n", confstring);
    exit(1);
  }
  fgets(idstring, MAXPATH, fid_conf);
  sscanf(idstring, "%d", &nfiles);
  fprintf(stderr, "There are %d files.\n", nfiles);

  /* these are for the row sparse format files */
  strcpy(fn_w,BINDIR);
  strcat(fn_w,"w");
  strcat(fn_w,A_OUTFILE);
  fid_w = fopen(fn_w, "w");
  assert(fid_w);

  strcpy(fn_j,BINDIR);
  strcat(fn_j,"j");
  strcat(fn_j,A_OUTFILE);
  fid_j = fopen(fn_j, "w");
  assert(fid_j);

  strcpy(fn_m,BINDIR);
  strcat(fn_m,"m");
  strcat(fn_m,A_OUTFILE);
  fid_m = fopen(fn_m, "w");
  assert(fid_m);

  n_elem_exported = 0;
  
  /* make the first number in the block file be number of images */
  assert(fwrite(&nfiles, sizeof(int), 1, fid_block) == 1);

  /* get the submatrix (nAA, vAA, yy, dd) for an individual image
   * and add it to the large matrix by writing it to disk 
   */
  
  for (i = 0; i < nfiles; i++) {
    fgets(idstring, MAXPATH, fid_conf);
    if (idstring[strlen(idstring) - 1] == '\n') {
      idstring[strlen(idstring) - 1] = '\0';
    }
    fprintf(stderr, "\ncurrent file: %s, %d of %d files\n", idstring,
            i + 1, nfiles);
    fprintf(fid_log, "%s ", idstring);

    /* output the start row of the data associated with the i'th image
       in A */
    assert(fwrite(&len_y, sizeof(int), 1, fid_block) == 1);

    rcs = rcs_llist_create();

/* Albert's test printouts  
    fprintf(stderr,"\n%p\n",rcs);
    fprintf(stderr,"\n%p\n",yy);
    fprintf(stderr,"\n%p\n",dd);
    fprintf(stderr,"\nSo far, so good...\n\n");
    fprintf(stderr,"\n%p\n",fid_log);
    fprintf(stderr,"\n%p\n",fid_date);
    */
    
    build_subA(idstring, rcs, &yy, &dd, &len_yy, fid_log, fid_date);

    /* Need to step back one position to compensate for the extra end
     * element written at the last step - this element is equal to the
     * number of elements written so far */
    if (i > 0) {
      assert(fseek(fid_m, -1*sizeof(int), SEEK_CUR) != -1);
    }
    
    rcs_llist_export_bin(rcs,fid_w, fid_j, fid_m, &n_elem_exported);
    rcs_llist_destroy(&rcs);
    assert(len_yy != 0);
    len_y += len_yy;

    if ((er = fwrite(yy, sizeof(float), len_yy, fid_y)) != len_yy) {
      fprintf(stderr, "error in 3rd fwrite, %d out of %d\n", er,
              len_yy);
      exit(3);
    }
    free(yy);
    if ((er = fwrite(dd, sizeof(float), len_yy, fid_delta)) != len_yy) {
      fprintf(stderr, "error in 4th fwrite, %d out of %d\n", er,
              len_yy);
      exit(3);
    }
    free(dd);
  }				/*end i loop (over files) */

  fclose(fid_w);
  fclose(fid_j);
  fclose(fid_m);
  
  
  /* output the total number of rows in A */
  assert(fwrite(&len_y, sizeof(int), 1, fid_block) == 1);
    
  /* print some housekeeping data and close files */
  fprintf(stdout, "\nThe matrices have %d rows and %d columns\n", len_y, nc3);
  fprintf(stdout, "\nFilenames: %s*%s\n", BINDIR,a_outfile);
  fprintf(fid_log, "\nThe matrices have %d rows and %d columns\n", len_y, nc3);

  fclose(fid_conf);
  fclose(fid_y);
  fclose(fid_delta);
  fclose(fid_block);
  fclose(fid_date);

  strcpy(filename_v, BINDIR);
  strcat(filename_v, "v");
  strcat(filename_v, a_outfile);

  strcpy(filename_i, BINDIR);
  strcat(filename_i, "i");
  strcat(filename_i, a_outfile);

  strcpy(filename_n, BINDIR);
  strcat(filename_n, "n");
  strcat(filename_n, a_outfile);

  fprintf(stderr,"\nStarting row sparse to col sparse conversion...");

  /* convert A matrix from RCS forma to CCS format  */
  convert_rcs_2_ccs(fn_w, fn_j, fn_m,
                    filename_v, filename_i, filename_n, nc3, len_y);

  fprintf(stderr,"..done\n");


  /* delete scratch files */
  /*
  sprintf(command,"rm -f %s\n",fn_w);
  assert(system(command) == 0);
  sprintf(command,"rm -f %s\n",fn_j);
  assert(system(command) == 0);
  sprintf(command,"rm -f %s\n",fn_m);
  assert(system(command) == 0);
  */
  
  return (0);
}

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

 // Note for Rich: In next three lines I surrounded the whole argument within extra (),
 // as suggested by previous warnings.
  assert((rcs_v_fid = fopen(rcs_v, "r")));
  assert((rcs_j_fid = fopen(rcs_j, "r")));
  assert((rcs_r_fid = fopen(rcs_r, "r")));

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

 // Note for Rich: In next three lines I surrounded the whole argument within extra (),
 // as suggested by previous warnings. 
  assert((ccs_v_fid = fopen(ccs_v, "w")));
  assert((ccs_i_fid = fopen(ccs_i, "w")));
  assert((ccs_n_fid = fopen(ccs_n, "w")));

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
