/*
 *  cvcalc.c by Richard Frazin 9/99 
 *
 *  this function takes some input parameters
 *   and calls solve.c to create an object x.
 *   the code then calls normcalc.c to evaluate
 *   the norm of x.
 *
 *   this has been modified to deal with a scalar lambda.  Note how the lambda array
 *     has been mofidied.  Also note how the first argument is not longer a float pointer,
 *     but instead is a single float.
 *
 * n_optarg is the number of optional arguments.  For now, the value must be
 *   0 or 1.  If n_optarg is 0, there is no argument following it in the 
 *   calling sequence and the 0th matrix name comes from FILESTR0  in 
 *   buildA_params.h.  If n_optarg is 1, the next argument is char * pointing
 *   to the name of the 0th matrix suffix.  See also "man stdarg", "solve.c", 
 *   "./Test_routines/stdarg_test.c"
 *
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdarg.h>

#include "headers.h"

#define CHARLEN (256)

#ifdef PRESERVE_SOL
static int i_count = 1;
#endif

float cvcalc(double lambda_value, char *norm_matrix_name, char *x_infile,
	     char *x_outfile, int n_optarg, ... )
{
  float lambda[NMATS];
  int huber_flag[] = HUBER_FLAG;
  float normval, waste1, waste2, waste3;
  int status, out_of_bounds;
#ifdef PRESERVE_SOL
  char buffer[MAXPATH];
#endif
  va_list pt_optarg;
  char *sys_matrix;
  
  if (n_optarg == 0){
    sys_matrix = NULL;
  } else if (n_optarg == 1) {
    va_start(pt_optarg, n_optarg);
    sys_matrix = va_arg(pt_optarg, char *);
    va_end(pt_optarg);
  } else {
      fprintf(stderr,"cv_calc: n_noptarg must be 0 or 1."); fflush(stderr);
      assert(0);
  }

  out_of_bounds = 0;
  lambda[0] = 1.0;
  lambda[1] = (float) lambda_value;

  fprintf(stdout, "\ncvcalc: x = %g\n", lambda_value);
  fflush(stdout);
  /*  fprintf(stdout, "\ncvcalc: x = %g %g %g \n", function_parameters[0], function_parameters[1], function_parameters[2]); */

  if (out_of_bounds < 1) {
    if (n_optarg == 0){ 
      status = solve(lambda, huber_flag, &waste1, &waste2, &waste3, x_infile, x_outfile,0);
    } else if (n_optarg == 1) {
      status = solve(lambda, huber_flag, &waste1, &waste2, &waste3, x_infile, x_outfile,1,sys_matrix);
    } else {
      assert(0);
    }
    normval = normcalc(norm_matrix_name, x_outfile);
    fprintf(stdout, "cvcalc: norm value = %g\n", normval);
    fflush(stdout);
  } else {
    normval = 1.0e32;
  }

#ifdef PRESERVE_SOL
  /* copy the solution at the current iteration so that it will not be
   * overwritten on the next iteration */
  sprintf(buffer, "cp %s%s %s%s_%d", BINDIR, x_outfile, BINDIR, x_outfile, i_count);
  i_count++;

  assert(system(buffer) == 0);
#endif
  
  return (normval);
}
