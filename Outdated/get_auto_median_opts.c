/* Return all defines that are relevant to auto_median
 *
 * Mark D. Butala
 * Oct. 2004
 *
 */

#include "mex.h"
#include "matrix.h"

#include "headers.h"
#include "auto_median.h"


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

  if (nlhs != 2)
    mexErrMsgTxt("2 outputs required.");


  buf = mxCalloc(MAXPATH, sizeof(char));
  strcpy(buf, AUTO_MEDIAN_CONF_FILE);
  plhs[0] = mxCreateString(buf);

  strcpy(buf, AUTO_MEDIAN_OUTFILE_SUFFIX);
  plhs[1] = mxCreateString(buf);
}
