#include <sys/types.h>
#include <sys/uio.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "headers.h"


/* FESS_HU:
 *   x[nc3] object   vector 
 *   nsubit  number of iterations done before return
 *   matrices -- structure with lots of stuff 
 */
void fess_hu(float *x, float *xw, float *lambda, struct matrix *matrices,
             int *boundary, float *constraint)
{

  const int nsubit = 1;

  int nc3 = NBINS;		/* number of pixels */
  int ig, k, l, fil;
  int subcount, bindx;
  float numer, denom, rr, dt, vvB, tmp;

  /* randomizer */
  /*
     for (k = 0; k < nc3; k++)
     *(index + k) = k ;

     srandom( (unsigned int) time(0) );
     for (k = 0; k < nc3; k++){
     l = random() % nc3;
     ig = *(index + l);
     *(index + l) = *(index + k);
     *(index + k) = ig;
     }
   */

  dt = 0;

  /* calculate updates of x and r.
   *  k gives the column number
   *  fil   and l give the row    number
   */
  for (k = 0; k < nc3; k++) {
    /* ig is the pixel number */
    ig = k;

    /*      ig = *(index + k); */
    xw[ig] = x[ig];
    for (subcount = 0; subcount < nsubit; subcount++) {
      numer = 0.0;
      denom = 0.0;
      for (fil = 0; fil < NMATS; fil++) {
        for (l = matrices[fil].n[ig]; l < matrices[fil].n[ig + 1];
             l++) {
          bindx = *(matrices[fil].iB + l);
#ifdef DEBUG

          if ((bindx < 0) || (bindx >= matrices[fil].nf)) {
            fprintf(stderr,
                    "\nFESS_HU: `bindx' out of range!\n");
            fprintf(stderr,
                    "bindx %d, k %d, pixel %d, fil %d, row %d\n",
                    bindx, k, ig, fil, l);
            exit(6);
          }
#endif
          rr = matrices[fil].r[bindx];
          vvB = *(matrices[fil].vB + l);
          vvB *= lambda[fil];

          if (matrices[fil].huber == 1) {	/* huber fcn */
            dt = matrices[fil].delta[bindx];
            if (rr > dt) {

              numer += vvB * dt;
              denom += vvB * vvB * dt / rr;
              /*above line should be replace by that below if nsubit > 1 */
              /*   denom += vvB*vvB*dt/(vvB*(xw[ig]-x[ig]) + rr ); */

            } else if (rr < -dt) {

              numer -= vvB * dt;
              denom -= vvB * vvB * dt / rr;
              /*above line should be replace by that below if nsubit > 1 */
              /*   denom -= vvB*vvB*dt/(vvB*(xw[ig]-x[ig]) + rr); */

            } else {

              numer += vvB * rr;
              denom += vvB * vvB;
              /* above line should be replace by that below if nsubit > 1 */
              /* numer += vvB*(vvB*(xw[ig]-x[ig]) + rr); */

            }
          } else {	/* quadratic fcn */

            numer += vvB * rr;
            denom += vvB * vvB;
            /* above line should be replace by that below if nsubit > 1 */
            /* numer += vvB*(vvB*(xw[ig]-x[ig]) + rr); */

          }
#ifdef DEBUG
          if (isnan((double) (rr)) || isinf((double) (rr))) {
            fprintf(stderr, "FESS_HU: error in 1st loop\n");
            fprintf(stderr, "(%d,%d):[%g, %g,%g,%g] ", fil,
                    bindx, rr, vvB, numer, denom);
            exit(7);
          }
#endif

        }		/* end l loop (over rows in file) */
      }			/* end fil files (over files) */
      if (denom != 0.0) {
        
#ifdef MINIMUM_VALUE
        if ((xw[ig] -= numer / denom) < MINIMUM_VALUE)
          xw[ig] = MINIMUM_VALUE;
#else
        xw[ig] -= numer / denom;
#endif
#ifdef CONSTRAINT_FILE
        if (xw[ig] < constraint[ig]) {
          xw[ig] = constraint[ig];
        }
#endif
#ifdef BOUNDARY_FILE
        if (boundary[ig] == 0) {
          xw[ig] = 0;
        }
#endif
      } else {

        fprintf(stderr, "FESS_HU: ZERO DENOMINATOR! pixel %d\n",
                ig);
        fprintf(stderr, "%d: %g \n", l, dt);

      }

    }			/* end subcount loop (over sub-iterations) */

    /* update r, requires separate loop over rows */
    tmp = xw[ig] - x[ig];
#ifdef DEBUG

    if (isnan((double) (tmp)) || isinf((double) (tmp))) {
      fprintf(stderr, "\ntmp is inf or NaN\n");
      fprintf(stderr,
              "k %d,pixel %d: numer %g, denom %g, xw %g, x %g, tmp %g \n",
              k, ig, numer, denom, xw[ig], x[ig], tmp);
      exit(0);
    }
#endif
    for (fil = 0; fil < NMATS; fil++) {
      for (l = matrices[fil].n[ig]; l < matrices[fil].n[ig + 1]; l++) {
        bindx = *(matrices[fil].iB + l);
        matrices[fil].r[bindx] +=
          *(matrices[fil].vB + l) * tmp * lambda[fil];
#ifdef DEBUG

        rr = matrices[fil].r[bindx] * lambda[fil];
        if (isnan((double) (rr)) || isinf((double) (rr))) {
          fprintf(stderr, "\nerror in 2nd loop\n");
          fprintf(stderr, "(%d,%d):[%g, %g,%g] ", fil, bindx, rr,
                  *(matrices[fil].vB + l), tmp);
          exit(0);
        }
#endif

      }			/* end l loop (over rows) */
    }			/* end fil files (over files) */

    x[ig] = xw[ig];

  }				/* end k loop (over pixels within block) */


  return;
}
