/* this file calls the solve routine to find the solution for a given
 * set of regularization parameters */


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <assert.h>

#include "headers.h"


#define CHARLEN (256)

static void usage(char *arg0) {
    printf("usage: %s <-h> [<infile> <outfile> <matrix suffix> <lamda1> [<lambda2>] [<lambda3>]]\n", arg0);
}

int main(int argc, char **argv) {
  float lambda[] = LAMBDA;
  char x_outfile[CHARLEN], x_infile[CHARLEN], matrix_name[CHARLEN];
  int huber_flag[] = HUBER_FLAG;
  float final_normx, final_normd, final_normt;
  int status;
  int verbose_flag = 1;
  int i;
  
  int opt, n_optarg;  /*n_optarg is for solve.c, nothing to do with command line */
  char optstring[] = "h";

    

  /* - type "man 3 getopt"
   * - optind is defined in getopt.h
   * - note the C 'switch' statement is more like a "goto".  Thus, if there is 
   *   no 'break' statement the commands below are executed even though they
   *   not satisfy the switch condition.  Thus, this code will work with fewer
   *   than 3 regularization parameters.
   * - make sure LAMBDA has the correct number of components in solve_cv_params.h
   */
  
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

  n_optarg = 0;
  if (argc - optind == 0) {
      strcpy(x_infile, MAIN_X_INFILE);
      strcpy(x_outfile, MAIN_X_OUTFILE);
  } else if ((argc - optind) >= 1) {
      switch (argc - optind) {
      case 6:
        lambda[3] = strtod(argv[optind+5], (char **)NULL);
      case 5:
        lambda[2] = strtod(argv[optind+4], (char **)NULL);
      case 4:
        lambda[1] = strtod(argv[optind+3], (char **)NULL);
        lambda[0] = 1;
      case 3:
        strcpy(matrix_name, argv[optind+2]);
	n_optarg = 1;
      case 2:
        strcpy(x_outfile, argv[optind+1]);
      case 1:
        strcpy(x_infile, argv[optind]);
        break;
      default:
        assert(0);
        break;
    }
  }  else {
    usage(argv[0]);
    return 1;
  }

  if (verbose_flag) {
    printf("lambda = { ");
    for (i = 0; i < NMATS; i++) {
      printf((i == NMATS - 1) ? "%g" : "%g, ", lambda[i]);
    }
    printf(" }\n");

    printf("huber_flag = { ");
    for (i = 0; i < NMATS; i++) {
      printf((i == NMATS - 1) ? "%d" : "%d, ", huber_flag[i]);
    }
    printf(" }\n");
  }

  if (n_optarg == 0){
    status = solve(lambda, huber_flag, &final_normx, &final_normd, &final_normt, 
		   x_infile, x_outfile,n_optarg);
  } else if (n_optarg == 1){
    status = solve(lambda, huber_flag, &final_normx, &final_normd, &final_normt, 
		   x_infile, x_outfile,n_optarg,matrix_name);
  }

  printf("status = %d\n", status);

  if (verbose_flag) {
    printf("final_normt = %g\n", final_normt);
  }
  
  return (0);
}
