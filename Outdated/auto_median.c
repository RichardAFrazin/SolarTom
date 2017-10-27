/*
 * auto_median.c
 *
 * by: Mark D. Butala 10/13/04
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <getopt.h>

#include "auto_median.h"
#include "headers.h"



/* For some reason the prototype is not in math.h on my system. */
double round(double x);

static int compare_d(const void *x, const void *y) {
  double *a, *b;

  a = (double *) x;
  b = (double *) y;
  
  return (*a < *b);
}

static void usage(char *arg0) {
  printf("usage: %s <-m> [k]\n", arg0);
  printf("-m: Compute medians but do not compute reconstructions\n");
  printf("[k]: Compute the k'th median reconstruction (first is k=0)\n");
}

static int process_conf(char base_names[][MAXPATH]) {
  FILE *conf_file;
  char buffer[MAXPATH];
  int num_lines;
  int n;
  
  sprintf(buffer, "%s%s", BINDIR, AUTO_MEDIAN_CONF_FILE);
  conf_file = fopen(buffer, "r");
  assert(conf_file != NULL);

  num_lines = 0;

  n = 0;
  
  do {
    if (n >= AUTO_MEDIAN_MAX_BASE_NAMES) {
      fprintf(stderr, "Increase AUTO_MEDIAN_MAX_BASE_NAMES (%d)!\n", AUTO_MEDIAN_MAX_BASE_NAMES);
      exit(1);
    }
    
    fscanf(conf_file, "%s\n", base_names[n]);
    n++;
  } while (!feof(conf_file));
    
  fclose(conf_file);

  return n;
}


int main(int argc, char **argv) {
  FILE *info_file;
  char base_names[AUTO_MEDIAN_MAX_BASE_NAMES][MAXPATH];
  char buffer[MAXPATH];
  int i, j, k, num_lambda;
  double lambda[3][MAX_NUMBER_OF_IMAGES], median[3][MAX_NUMBER_OF_IMAGES];
  int f, r;
  int n;
  int start = 0, stop = 0;
  
  int opt;
  char optstring[] = "mh";
  int median_only_flag = 0;

  /* process command line arguments */
  while ((opt = getopt(argc, argv, optstring)) != -1) {
    switch (opt) {
    case 'h':
      usage(argv[0]);
      return 0;
      break;
    case 'm':
      median_only_flag = 1;
      break;
    default:
      usage(argv[0]);
      return 1;
    }
  }

  n = process_conf(base_names);

  if (argc - optind == 0) {
    start = 0;
    stop = n - 1;
  }
  else if (argc - optind == 1) {
    start = strtol(argv[optind], (char **)NULL, 10);
    stop = start;

    if (start < 0) {
      fprintf(stderr, "k >= 0!\n");
      return 2;
    }

    if (stop >= n) {
      fprintf(stderr, "k < n!\n");
      return 2;
    }
  }
  else {
    usage(argv[0]);
    return 1;
  }
  
  for (i = 0; i < n; i++) {
    sprintf(buffer, "%sinfo_%s_auto_cv", BINDIR, base_names[i]);
    
    info_file = fopen(buffer, "r");
    printf("%s\n", buffer);
    assert(info_file != NULL);

    j = 0;
    
    do {
      for (k = 0; k < 3; k++) {
        fscanf(info_file, "%lg", &lambda[k][j]);
      }

      fscanf(info_file, "\n");
      j = j + 1;
    } while (!feof(info_file));

    num_lambda = j;
    
    fclose(info_file);
    
    for (j = 0; j < 3; j++) {
      qsort(lambda[j],num_lambda,sizeof(double),&compare_d);
    }

    f = floor(num_lambda/2);
    r = round(num_lambda/2);

    if (f == r) {
      for (j = 0; j < 3; j++) {
        median[j][i] = (lambda[j][f-1] + lambda[j][f])/2;
      }
    }
    else {
      for (j = 0; j < 3; j++) {
        median[j][i] = lambda[j][f];
      }
    }

    if (!median_only_flag && i >= start && i <= stop) {
      printf("Processing: %s\n", base_names[i]);
      fflush(stdout);
      
      sprintf(buffer, "./callsolve_cg -v %s x_%s%s %s %g %g %g", MAIN_X_INFILE, base_names[i], AUTO_MEDIAN_OUTFILE_SUFFIX, base_names[i], median[0][i], median[1][i], median[2][i]);
      
      printf("Executing: %s\n", buffer);
      fflush(stdout);
      assert(system(buffer) == 0);
      
      sprintf(buffer, "./callsolve_fess -v x_%s%s x_%s%s %s %g %g %g", base_names[i], AUTO_MEDIAN_OUTFILE_SUFFIX, base_names[i], AUTO_MEDIAN_OUTFILE_SUFFIX, base_names[i], median[0][i], median[1][i], median[2][i]);

      printf("Executing: %s\n", buffer);
      fflush(stdout);
      assert(system(buffer) == 0);
      
      printf("\n");
    }
  }

  for (i = start; i <= stop; i++) {
    printf("(%d)\t%s:\t\t", i, base_names[i]);
    for (j = 0; j < 3; j++) {
      printf("%g\t", median[j][i]);
    }
    printf("\n");
  }

  return 0;
}
