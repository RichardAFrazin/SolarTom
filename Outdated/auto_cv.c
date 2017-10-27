/* auto_cv.c
 *
 * by: Mark D. Butala  9/04, modified R. Frazin 6/07
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <signal.h>
#include <time.h>
#include <getopt.h>

#include "headers.h"
#include "amoeba.h"


static char auto_cv_norm_matrix[MAXPATH];
static char auto_cv_system_matrix[MAXPATH];
static char x_outfile[MAXPATH];

static void load_A_y_delta(struct sparse *A, float **y, float **delta, char *matrix_name);
  static void initialize(int *n_blocks, int *block_start, char *full_matrix_name);

float cvcalc_wrapper(float *function_parameters) {
  return cvcalc(function_parameters, auto_cv_norm_matrix, auto_cv_system_matrix, x_outfile, x_outfile);
}


static void usage(char *arg0) {
  printf("usage: %s [<-d>] [-jN] [<-m A_OUTFILE>] [<-o CV_X_OUTFILE>] [<-s SUFFIX>][k]\n", arg0);
  printf("-m: Specify matrix suffix\n");
  printf("-o: Specify output file name\n");
  printf("-s: Specify suffix.  Then A_OUTFILE=suffix and  CV_X_OUTFILE=x_suffix.\n");
  printf("-j: Only evaluate every Nth image block (i.e., 0, N, 2*N, ...).  Do not combine with 'k' option\n");
  printf("[k]: Compute the k'th auto_cv (first is k=0)\n");
  printf("-d: Enabled distributed mode.  Two more parameters are then mandatory: [cpu_id] and [num_cpu].\n");
}

int main(int argc, char **argv) {
  float p0[NDIM + 1][NDIM] = AMOEBA_LAMBDA;
  float p[NDIM + 1][NDIM];
  float x[NDIM], yy[NDIM + 1];
  int i, j, k, jump;
  int status, iter;
  float (*func) (float *);
  int n_blocks;
  int block_start[MAX_NUMBER_OF_IMAGES+1];
  int row_start, row_end;
  struct sparse A;
  char full_matrix_name[MAXPATH] = A_OUTFILE;
  char infile_name[MAXPATH] = CV_X_INFILE;
  char outfile_name[MAXPATH] = CV_X_OUTFILE;
  float *y, *delta;
  FILE *auto_cv_info_fid = NULL;
  time_t current_time;
  char command[MAXPATH], buffer[MAXPATH];
  int start, stop;
  int cpu_id = 0, num_cpu = 0;
  
  int opt;
  char distributed_mode = 0;
  char optstring[] = "dhj:m:o:s:"; /* set valid command line options */

  jump = 1;

  /* process command line arguments */
  while ((opt = getopt(argc, argv, optstring)) != -1) {
    switch (opt) {
    case 'h':
      usage(argv[0]);
      return 0;
      break;
    case 'd':
      distributed_mode = 1;
      break;
    case 'j':
      jump = atoi(optarg);
      if (jump < 1){
	printf("%s: N >= 1 only!\n",argv[0]);
	usage(argv[0]);
	return 1;
      }
      break;
    case 'm':
      strncpy(full_matrix_name, optarg, MAXPATH);
      break;
    case 'o':
      strncpy(outfile_name, optarg, MAXPATH);
      break;
    case 's':
      strncpy(full_matrix_name, optarg, MAXPATH);
      strcpy(outfile_name, "x_");
      strncat(outfile_name, optarg, MAXPATH-2);
      break;
    default:
      usage(argv[0]);
      return 1;
    }
  }

  initialize(&n_blocks, block_start, full_matrix_name);

  if (argc - optind == 0) {
    start = 0;
    stop = n_blocks;

    strcpy(auto_cv_norm_matrix, AUTO_CV_NORM_MATRIX);
    strcpy(auto_cv_system_matrix, AUTO_CV_SYSTEM_MATRIX);
    
    /* Open the auto_cv info file */
    sprintf(buffer, "%sinfo_%s_auto_cv", BINDIR, full_matrix_name);
    auto_cv_info_fid = fopen(buffer, "w");
    assert(auto_cv_info_fid != NULL);
  }
  else if (argc - optind == 1) {
    if (jump != 1){
      fprintf(stderr,"%s: Do not use -j option with the k option\n",argv[0]);
      fflush(stderr);
      assert(0);
    }
    start = strtol(argv[optind], (char **)NULL, 10);
    stop = start + 1;
    
    if (start < 0) {
      fprintf(stderr, "k >= 0!\n");
      return 2;
    }
    
    if (stop >= n_blocks) {
      fprintf(stderr, "k < n_blocks!\n");
      return 2;
    }
    

    sprintf(auto_cv_norm_matrix, "%s_%d", AUTO_CV_NORM_MATRIX, start);
    sprintf(auto_cv_system_matrix, "%s_%d", AUTO_CV_SYSTEM_MATRIX, start);
    
    /* Open the auto_cv info file */
    sprintf(buffer, "%sinfo_%s_auto_cv_%d", BINDIR, full_matrix_name, start);
    auto_cv_info_fid = fopen(buffer, "w");
    assert(auto_cv_info_fid != NULL);
  }
  else if (argc - optind == 2) {
    if (distributed_mode) {
      start = 0;
      stop = n_blocks;
      
      cpu_id = strtol(argv[optind], (char **)NULL, 10);
      num_cpu = strtol(argv[optind+1], (char **)NULL, 10);
      
      if (cpu_id < 0) {
        fprintf(stderr, "cpu_id >= 0!\n");
        return 2;
      }

      if (cpu_id >= num_cpu) {
        fprintf(stderr, "cpu_id < num_cpu!\n");
      }
      
      if (num_cpu < 0) {
        fprintf(stderr, "num_cpu >= 0!\n");
        return 2;
      }

      sprintf(auto_cv_norm_matrix, "%s_%d", AUTO_CV_NORM_MATRIX, cpu_id);
      sprintf(auto_cv_system_matrix, "%s_%d", AUTO_CV_SYSTEM_MATRIX, cpu_id);
      
      /* Open the auto_cv info file */
      sprintf(buffer, "%sinfo_%s_auto_cv_%d", BINDIR, full_matrix_name, cpu_id);
      auto_cv_info_fid = fopen(buffer, "w");
      assert(auto_cv_info_fid != NULL);
    }
    else {
      usage(argv[0]);
      return 1;
    }      
  }
  else {
    usage(argv[0]);
    return 1;
  }
  
 
  printf("n_blocks: %d\n", n_blocks);

  if (!distributed_mode) {
    for (k = start; k < stop + 1; k = k + jump ) {
      printf("%d %d\n", k,block_start[k]);
    }
  }
 
  for (k = start; k < stop; k = k + jump) {
    if (distributed_mode) {
      if (k % num_cpu != cpu_id) {
        continue;
      }
    }

    printf("auto_cv: begin iter %d\n", k);
    current_time = time(NULL);
    printf("%s", ctime(&current_time));
    fflush(stdout);

    /* initialize the p array, i.e. the array of lambdas */
    for (i = 0; i < NDIM + 1; i++) {
      for (j = 0; j < NDIM; j++) {
        p[i][j] = p0[i][j];
      }
    }
    
    /* sprintf(x_outfile, "%s_auto_cv_%d", outfile_name, k); */
    sprintf(x_outfile, "%s_auto_cv", outfile_name);

    /* first delete x_outfile and then copy infile_name to
     * x_outfile */
    printf("Removing: %s\n", x_outfile);
    fflush(stdout);
    sprintf(command, "rm -f %s%s", BINDIR, x_outfile);
    assert(system(command) == 0);
    printf("Copying: %s to %s\n", infile_name, x_outfile);
    fflush(stdout);
    sprintf(command, "cp %s%s %s%s", BINDIR, infile_name, BINDIR, x_outfile);
    assert(system(command) == 0);
    
    row_start = block_start[k];
    row_end = block_start[k+1]-1;

    load_A_y_delta(&A, &y, &delta, full_matrix_name);

    printf("full matrix: %s\n", full_matrix_name);
    printf("norm matrix: %s \t cv system matrix: %s\n",
           auto_cv_norm_matrix, auto_cv_system_matrix);
    printf("Extracting row %d-%d ... ", row_start, row_end);
    fflush(stdout);
    row_extract(A, y, delta, auto_cv_norm_matrix, auto_cv_system_matrix, row_start, row_end);
    printf("done\n");

    free_sparse(&A);
    free_vf(&y);
    free_vf(&delta);

    printf("Beginning amoeba optimization\n");
    fflush(stdout);

    func = FUNCTION_NAME;

    for (i = 0; i < NDIM + 1; i++) {
      for (j = 0; j < NDIM; j++)
        x[j] = p[i][j];
      yy[i] = func(x);
      printf("initial amoeba iteration %d: %g\n", i, yy[i]);
      printf("x: ");
      for (j=0; j < NDIM; j++)
        printf(" %g", x[j]);
      printf("\n");
      fflush(stdout);
    }
    
    status = amoeba(p, yy, func, &iter);

    printf("amoeba: status = %d, iter =  %d\n", status, iter);
    printf("function value =");
    fflush(stdout);

    for (j = 0; j < NDIM; j++) {
      x[j] = p[0][j];
      printf(" %g", x[j]);
      fprintf(auto_cv_info_fid, " %g", x[j]);
    }
    printf("\n\n");
    fflush(stdout);
    fprintf(auto_cv_info_fid, "\n");
    fflush(auto_cv_info_fid);
  }

  /* cleanup */
  fclose(auto_cv_info_fid);
  
  return (0);
}

void load_A_y_delta(struct sparse *A, float **y, float **delta, char *matrix_name) {
  char buffer[MAXPATH];
  
  /* Load A matrix */
  printf("Loading A (%s) ... ", matrix_name);
  fflush(stdout);
  load_sparse(A, matrix_name, BINDIR);
  printf("done\n");
  
  /* Load y and delta vectors */
  printf("Loading y (y%s) ... ", matrix_name);
  fflush(stdout);
  sprintf(buffer, "%sy%s", BINDIR, matrix_name);
  load_vf(y, A->nf, buffer);
  printf("done\n");
  printf("Loadind delta (delta_%s) ... ", matrix_name);
  fflush(stdout);
  sprintf(buffer, "%sdelta_%s", BINDIR, matrix_name);
  load_vf(delta, A->nf, buffer);
  printf("done\n");
}

void initialize(int *n_blocks, int *block_start, char *full_matrix_name) {
  FILE *fid_info;
  char buffer[MAXPATH];
  int i, info_file_data, found_newline;
  
  strcat(strcat(strcpy(buffer, BINDIR),"info_"), full_matrix_name);
  fid_info = fopen(buffer, "r");
  assert(fid_info != NULL);

  found_newline = 0;
  
  for (i = 0; i < MAX_NUMBER_OF_IMAGES+1; i++) {
    assert(fread(&info_file_data, sizeof(int), 1, fid_info) == 1);

    if ((char) info_file_data == '\n') {
      found_newline = 1;
      break;
    }
    else {
      block_start[i] = info_file_data;
    }
  }

  if (!found_newline) {
    assert(fread(&info_file_data, sizeof(int), 1, fid_info) == 1);

    if ((char) info_file_data != '\n') {
      fprintf(stderr, "More than MAX_NUMBER_OF_IMAGES (%d) number of blocks in file %s\n!", MAX_NUMBER_OF_IMAGES, buffer);
      exit(1);
    }
  }

  fclose(fid_info);

  *n_blocks = i - 1;
}
