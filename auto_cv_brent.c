/* auto_cv_brent.c
 *
 * this uses the brent algorithm to do the cross validation
 *AUTO_CV_SYSTEM_MATRIX
 *  R. Frazin 6/07 (based on Butala's auto_cv.c)
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


static char auto_cv_norm_matrix[MAXPATH];
static char auto_cv_sys_matrix[MAXPATH];
static char x_outfile[MAXPATH];  

static void load_A_y_delta(struct sparse *A, float **y, float **delta, char *matrix_name);
static void initialize(int *n_blocks, int *block_start, char *full_matrix_name);
double fminbr(double, double, double (*)( ), double);

double dbl_cvcalc_wrapper(double lambda_value) {
  return (double) cvcalc(lambda_value, auto_cv_norm_matrix, x_outfile, x_outfile, 1, auto_cv_sys_matrix);
}

double f2(x)   /*test function */
double x;
{
  fprintf(stdout,"f2: x = %g\n",x);
  fflush(stdout);
  return pow( (pow(x,2)-2.0)*x - 5.0, 2 );
}

static void usage(char *arg0) {
  printf("usage: %s [-jN] [-rM] [-wP] [<-m matrix_suffix >] [<-i x_infile>] [<-o x_outfile>] [k]\n", arg0);
  printf("-h: Print this helpful information\n");
  printf("-j: Only evaluate every Nth image block (i.e., 0, N, 2*N, ...).\n    Do not combine with '-r', '-w' or 'k' options\n");
  printf("-r: Remove M images for each CV evaluation.\n    Thus, the first evaluation removes images 0-(M-1), the seccond M-(2M-1), etc.\n    Must be used with the '-w' option.  Do not combine with '-j' or 'k' options.\n");
  printf("-w: Use the central P images of the M images removed for the CV cost calculation.\n    Thus, if this routine is called with the '-r7 -w3', for the first evaluation\n     images 0-6 will be removed and images 2-4 will be used for the CV cost calculation.\n     Must be used with the '-r' option. Do not combine with '-j' or 'k' options.\n");
  printf("-m: Specify matrix suffix\n");
  printf("-o: Specify output file name\n");
  printf("[k]: Compute the k'th auto_cv (first is k=0)\n");
  /*printf("-d: distributed mode (disabled).  Two more parameters are then mandatory: [cpu_id] and [num_cpu].\n");*/
}

int main(int argc, char **argv) {
  double brent_upper_bound, brent_lower_bound, brent_min, brent_tol; /* brent parameters */
  int k, i, jump;
  double (*func) (double *);
  int n_blocks;
  int block_start[MAX_NUMBER_OF_IMAGES+1];
  int row_start, row_end, row_cv_start, row_cv_end;
  struct sparse A;
  char full_matrix_name[MAXPATH] = FILESTR0;
  char infile_name[MAXPATH] = CV_X_INFILE, outfile_name[MAXPATH] = CV_X_OUTFILE;
  char *time_str, time_string[33];
  float *y, *delta;
  FILE *auto_cv_info_fid = NULL;
  time_t current_time;
  char command[MAXPATH], buffer[MAXPATH];
  int start, stop;
  int remove_images, cv_images;
  int cpu_id = 0, num_cpu = 0;
  int opt;
  char distributed_mode = 0;
  char optstring[] = "hj:r:w:m:o:i:"; /* set valid command line options */

  func =  dbl_cvcalc_wrapper;
  /*func = f2;*/    /*test function, see above */

  /* set up scratch file names for CV matrices */
  current_time = time(NULL);
  time_str = ctime(&current_time);
  fprintf(stdout,"The current time is: %s\n",time_str);
  fflush(stdout);
  k = 0; i = 0;
  while ( *(time_str + k) != '\0'){
    /* remove spaces from string */
    if (*(time_str + k) != 32){
      strncpy(time_string + i , time_str + k,1);
      i++;
    }
    k++;
  }
  sprintf(time_string + i-1,"%s","\0");
  strcpy(auto_cv_norm_matrix,"cv_norm_");
  strcpy(auto_cv_sys_matrix,"cv_system_");
  strcat(auto_cv_norm_matrix,time_string);
  strcat(auto_cv_sys_matrix,time_string);


  jump = 1; remove_images = 1; cv_images = 1;

  /* process command line arguments */
  while ((opt = getopt(argc, argv, optstring)) != -1) {
    switch (opt) {
    case 'h':
      usage(argv[0]);
      return 0;
      break;
    case 'j':
      jump = atoi(optarg);
      if (jump < 1){
	printf("%s: N >= 1 only!\n",argv[0]);
	usage(argv[0]);
	return 1;
      }
      break;
    case 'r':
      remove_images = atoi(optarg);
      if ( (remove_images > 1) && (jump > 1) ){
	fprintf(stderr,"Do not use -j and -r options together!\n");
	usage(argv[0]);
	return 1;
      }
      break;
    case 'w':
      cv_images = atoi(optarg);
      if ( (cv_images > 1) && (jump > 1) ){
	fprintf(stderr,"Do not use -j and -w options together!\n");
	usage(argv[0]);
	return 1;
      }
      if ( remove_images < cv_images ){
	fprintf(stderr,"M must be greater than P!\n");
	usage(argv[0]);
	return 1;
      }
      break;
    case 'm':
      strncpy(full_matrix_name, optarg, MAXPATH);
      break;
    case 'i':
      strncpy(infile_name, optarg, MAXPATH);
      break;
    case 'o':
      strncpy(outfile_name, optarg, MAXPATH);
      break;
    default:
      usage(argv[0]);
      return 1;
    }
  }

  if (remove_images > 1)
    jump = remove_images;

  initialize(&n_blocks, block_start, full_matrix_name);

  if (argc - optind == 0) {
    start = 0;
    stop = n_blocks;

    /* Open the auto_cv info file */
    sprintf(buffer, "%scv_info_%s", BINDIR, full_matrix_name);
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
    
    /*
    sprintf(auto_cv_norm_matrix, "%s_%d", AUTO_CV_NORM_MATRIX, start);
    sprintf(auto_cv_sys_matrix, "%s_%d", AUTO_CV_SYSTEM_MATRIX, start);
    */

    /* Open the auto_cv info file */
    sprintf(buffer, "%sinfo_%s_auto_cv_%d", BINDIR, full_matrix_name, start);
    auto_cv_info_fid = fopen(buffer, "w");
    assert(auto_cv_info_fid != NULL);
  }
  else if (argc - optind == 2) {
    if (distributed_mode) { /* disabled */
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
  
 
  printf("There are %d images.\n", n_blocks);

  for (k = start; k < stop + 1; k++) {
    printf("%d %d\n", k,block_start[k]);
  }
 
  for (k = start; k < stop; k = k + jump) {
 
    printf("auto_cv: begin file %d\n", k);
    fflush(stdout);
    
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
    
    load_A_y_delta(&A, &y, &delta, full_matrix_name);

    printf("full matrix: %s\n", full_matrix_name);
    printf("norm matrix: %s \ncv system matrix: %s\n",
           auto_cv_norm_matrix, auto_cv_sys_matrix);
   
    row_start = block_start[k];
    row_end = block_start[k + remove_images] - 1;
    if ( (cv_images == 1) && (remove_images == 1) ){
      row_cv_start = row_start;
      row_cv_end   = row_end;
    } else {
      row_cv_start = block_start[k + (int) lround( ((float) remove_images - (float) cv_images)/2. ) ];
      row_cv_end   = block_start[k + (int) lround( ((float) remove_images - (float) cv_images)/2. ) + cv_images] - 1;
    }

    fprintf(stdout,"Removing rows %d-%d, Writing rows %d-%d ...", row_start, row_end,row_cv_start, row_cv_end);
    fflush(stdout);

    /*the row_extract function is in sparse.c */
    row_extract_mod(A, y, delta, auto_cv_norm_matrix, auto_cv_sys_matrix, row_start, row_end, row_cv_start, row_cv_end);
    printf("done\n");

    free_sparse(&A);
    free_vf(&y);
    free_vf(&delta);

    fprintf(stdout,"Beginning brent optimization\n");
    fflush(stdout);

    brent_lower_bound = MIN_LAMBDA;
    brent_upper_bound = MAX_LAMBDA;
    brent_tol = BRENT_TOL;
 
    brent_min = fminbr(brent_lower_bound, brent_upper_bound, func,brent_tol);
    
    fprintf(stdout,"function value = %g\n\n",brent_min);
    fflush(stdout);

    fprintf(auto_cv_info_fid, " %g\n", brent_min);
    fflush(auto_cv_info_fid);

    /* remove the scratch matrices */
    printf("Removing: %s and %s\n", auto_cv_sys_matrix,auto_cv_norm_matrix);
    fflush(stdout);
    sprintf(command, "rm -f %s*%s", BINDIR, auto_cv_sys_matrix);
    assert(system(command) == 0);
    sprintf(command, "rm -f %s*%s", BINDIR, auto_cv_norm_matrix);
    assert(system(command) == 0);
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
  
  strcat(strcat(strcpy(buffer, BINDIR),"block_"), full_matrix_name);
  fid_info = fopen(buffer, "r");
  assert(fid_info != NULL);

  assert(fread(&info_file_data, sizeof(int), 1, fid_info) == 1);
  *n_blocks = info_file_data;
  fprintf(stderr,"Entry #1: %d\n",info_file_data);
  for (i = 0; i < *n_blocks+1; i++) {
    assert(fread(&info_file_data, sizeof(int), 1, fid_info) == 1);
    block_start[i] = info_file_data;
    fprintf(stderr,"Entry # %d: %d\n",i+2,info_file_data);
  }
  fclose(fid_info);  
  //  exit(0);
}
