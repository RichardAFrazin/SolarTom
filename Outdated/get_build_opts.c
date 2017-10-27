/* Return various build options to matlab for use with the functions
 * that build the auxillary matricies.
 *
 * Mark D. Butala
 * Jul. 2004
 *
 */

#include "mex.h"
#include "matrix.h"

#include "headers.h"


#ifdef CARTESIAN
#error This code can only handle cylindrical reconstructions
#endif

#ifdef C2BUILD
#define BUILD       "C2BUILD"
#define BUILD_RMIN  C2_RMIN
#define BUILD_RMAX  C2_RMAX
#define IMSIZE      C2_IMSIZE
#elif defined MK4BUILD
#define BUILD       "MK4BUILD"
#define BUILD_RMIN  MK4_RMIN
#define BUILD_RMAX  MK4_RMAX
#define IMSIZE      MK4_IMSIZE
#elif defined MK4_POL_BUILD
#define BUILD       "MK4_CYL_BUILD"
#define BUILD_RMIN  MK4_RMIN
#define BUILD_RMAX  MK4_RMAX
#define IMSIZE      0
#else
#error Build not recognized
#endif

mxArray *myCreateDoubleScaler(double value) {
  mxArray *array;

  array = mxCreateDoubleMatrix(1, 1, mxREAL);
  *mxGetPr(array) = value;

  return array;
}

/* The gateway routine */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]) {
  char *buf;
  
  /* Check for proper number of arguments. */
  /* NOTE: You do not need an else statement when using
     mexErrMsgTxt within an if statement. It will never
     get to the else statement if mexErrMsgTxt is executed.
     (mexErrMsgTxt breaks you out of the MEX-file.) 
  */
  if (nrhs != 0)
    mexErrMsgTxt("Zero inputs required.");

  if (nlhs != 13)
    mexErrMsgTxt("12 outputs required.");


  plhs[0] = myCreateDoubleScaler(NRAD);
  plhs[1] = myCreateDoubleScaler(NPHI);
  plhs[2] = myCreateDoubleScaler(NZ);

  plhs[3] = myCreateDoubleScaler(BUILD_RMIN);
  plhs[4] = myCreateDoubleScaler(BUILD_RMAX);
  plhs[5] = myCreateDoubleScaler(RMAX);

  buf = mxCalloc(MAXPATH, sizeof(char));
  strcpy(buf, BINDIR);
  plhs[6] = mxCreateString(buf);

  strcpy(buf, MAIN_X_INFILE);
  plhs[7] = mxCreateString(buf);

  strcpy(buf, CV_X_INFILE);
  plhs[8] = mxCreateString(buf);

  plhs[9] = myCreateDoubleScaler(IMSIZE);

  strcpy(buf, BUILD);
  plhs[10] = mxCreateString(buf);

#ifdef CARTISIAN
  strcpy(buf, "CARTISIAN");
  plhs[11] = mxCreateString(buf);
#elif defined CYLINDRICAL
  strcpy(buf, "CYLINDRICAL");
  plhs[11] = mxCreateString(buf);
#else
#error Computation grid not recognized
#endif

  strcpy(buf, MAIN_X_INFILE);
  plhs[12] = mxCreateString(buf);
}
