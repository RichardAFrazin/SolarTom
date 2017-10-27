/*       cv_brent_fixed.c - by Richard Frazin, May 2008
 *
 * This code does the single parameter cross validation to 
 *  determine the regularization parameter.  Unlike 
 *  auto_cv_brent.c, which extracts submatrices to be used
 *  for the cross validation, this code uses a fixed matrix
 *  for both the cross validation (cv_martrix) and the 
 *  reconstruction (sys_matrix).  This code is command 
 *  line driven.
 */

#include <string.h>
#include <signal.h>
#include <time.h>
#include <getopt.h>
#include "headers.h"

static char x_outfile[MAXPATH], x_infile[MAXPATH];
char cv_matrix[MAXPATH], sys_matrix[MAXPATH];

double fminbr(double, double, double (*)( ), double);
double dbl_cvcalc_wrapper(double lambda_value) {
  return (double) cvcalc(lambda_value, cv_matrix, x_infile, x_outfile,1,sys_matrix);
}

double f2(x)   /*test function */
double x;
{
  fprintf(stdout,"f2: x = %g\n",x);
  fflush(stdout);
  return pow( (pow(x,2)-2.0)*x - 5.0, 2 );
}

static void usage(char *arg0) {
  printf("usage: %s [-h] x_infile x_outfile system_matrix cv_matrix \n", arg0);
  printf("-h: Print this helpful information\n");
  printf(" x_infile is the initial guess of x and the x_outfile is the output x.\n sys_matrix is the one used for reconstruction.  cv_matrix is the matrix used for the cross-validation cost\n"); 
}

int main(int argc, char **argv) {
  double brent_upper_bound, brent_lower_bound, brent_min, brent_tol; /* brent parameters */
  double (*func) (double *);
  char optstring[] = "h"; /* set command line options*/
  int opt;

  while ((opt = getopt(argc, argv, optstring)) != -1) {
    switch (opt) {
    case 'h':
      usage(argv[0]);
      return(0);
      break;
    default:
      usage(argv[0]);
      return(1);
    }
  }
  if ( argc != 5){
    usage(argv[0]);
    return(0);
  }

  strcpy(x_infile, argv[1]);
  strcpy(x_outfile,argv[2]);
  strcpy(sys_matrix,argv[3]);
  strcpy(cv_matrix,argv[4]);

  func = dbl_cvcalc_wrapper;

  brent_lower_bound = MIN_LAMBDA;
  brent_upper_bound = MAX_LAMBDA;
  brent_tol = BRENT_TOL;
 
  brent_min = fminbr(brent_lower_bound, brent_upper_bound, func,brent_tol);
    
  fprintf(stdout,"function value = %g\n\n",brent_min);
  fflush(stdout);

  return(0);
}

