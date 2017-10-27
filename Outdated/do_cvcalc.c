#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>

#include "headers.h"

void usage(char *arg0) {
  /*
  printf("usage: %s <lambda1> <lambda2> <lambda3> <norm_matrix> <solve_matrix> <x_infile> <x_outfile>\n", arg0);
  */
  printf("usage: %s <lambda> <norm_matrix> <solve_matrix> <x_infile> <x_outfile>\n", arg0);
}

int main(int argc, char **argv) {
  float normval;
  /*
  float lambda[3];
  */
  float lambda;
  char *norm_matrix_name, *solve_matrix_name, *x_infile, *x_outfile;
  
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

  /*  if ((argc - optind) == 7) { */

  if ((argc - optind) == 5) {
    lambda = strtod(argv[optind], (char **)NULL);
    norm_matrix_name = argv[optind+1];
    solve_matrix_name = argv[optind+2];
    x_infile = argv[optind+3];
    x_outfile = argv[optind+4];
    /*
    lambda[0] = strtod(argv[optind], (char **)NULL);
    lambda[1] = strtod(argv[optind+1], (char **)NULL);
    lambda[2] = strtod(argv[optind+2], (char **)NULL);
    norm_matrix_name = argv[optind+3];
    solve_matrix_name = argv[optind+4];
    x_infile = argv[optind+5];
    x_outfile = argv[optind+6];
    */
  }
  else {
    usage(argv[0]);
    return 1;
  }
  
  normval = cvcalc(lambda,
                   norm_matrix_name,
                   solve_matrix_name,
                   x_infile, x_outfile);

  printf("normval: %f\n", normval);
  
  return 0;
}
