/* amoeba.c, amotry.c, testfunc.c
 *
 *    by Richard Frazin 9/99
 *
 * numerical recipes optimization algorithm
 *
 * testfunc is a function to test the code
 */

#include <sys/types.h>
#include <sys/uio.h>
#include <malloc.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "headers.h"
#include "amoeba.h"


static float amotry(float[NDIM + 1][NDIM], float *, float *,
                    float (*)(float *),
                    int, float);


int amoeba(p, y, func, iter)
float p[NDIM + 1][NDIM], y[NDIM + 1];
float (*func) (float *);
int *iter;
{
  const int itmax = AMOEBA_ITMAX, ndim = NDIM;
  const float ftol = FTOL;
  int i, ihi, ilo, inhi, j, m, n;
  float rtol, sum, swap, ysave, ytry, psum[NDIM];
  int status;

  status = 0;
  *iter = 0;

one:
  for (n = 0; n < ndim; n++) {
    sum = 0.;
    for (m = 0; m < ndim + 1; m++)
      sum += p[m][n];
    psum[n] = sum;
  }

two:
  ilo = 0;
  if (y[0] > y[1]) {
    ihi = 0;
    inhi = 1;
  } else {
    ihi = 1;
    inhi = 0;
  }

  for (i = 0; i < ndim + 1; i++) {
    if (y[i] <= y[ilo])
      ilo = i;
    if (y[i] > y[ihi]) {
      inhi = ihi;
      ihi = i;
    } else if (y[i] > y[inhi]) {
      if (i != ihi)
        inhi = i;
    }
  }									  /* end for loop */

  rtol = 2. * fabs(y[ihi] - y[ilo]) / (fabs(y[ihi]) + fabs(y[ilo]));

  if (rtol < ftol) {
    swap = y[0];
    y[0] = y[ilo];
    y[ilo] = swap;
    for (n = 0; n < ndim; n++) {
      swap = p[0][n];
      p[0][n] = p[ilo][n];
      p[ilo][n] = swap;
    }
    fprintf(stderr, "AMOEBA: stopping criterion met\n");
    return (status);
  }

  if (*iter >= itmax) {
    fprintf(stderr, "AMOEBA: max numerber of iterations exceeded!!\n");
    status = 1;
    return (status);
  }

  *iter += 2;
  ytry = amotry(p, y, psum, func, ihi, -1.0);

  if (ytry <= y[ilo]) {
    ytry = amotry(p, y, psum, func, ihi, EXPAN_FAC);
  } else if (ytry >= y[inhi]) {
    ysave = y[ihi];
    ytry = amotry(p, y, psum, func, ihi, 1.0 / EXPAN_FAC);
    if (ytry >= ysave) {
      for (i = 0; i < ndim + 1; i++) {
        if (i != ilo) {
          for (j = 0; j < ndim; j++) {
            psum[j] = .5 * (p[i][j] + p[ilo][j]);
            p[i][j] = psum[j];
          }
          y[i] = func(psum);
        }
      }						  /* endfor (i loop) */
      *iter += ndim;
      goto one;
    }							  /* end if ytry >= ysave */
  } else {
    *iter -= 1;
  }									  /* end big if */
  goto two;
}

float amotry(p, y, psum, func, ihi, fac)
float p[NDIM + 1][NDIM], *y, *psum;
float (*func) (float *);
int ihi;
float fac;
{
  int j;
  float fac1, fac2, ytry, ptry[NDIM];
  const int ndim = NDIM;

  fac1 = (1. - fac) / ndim;
  fac2 = fac1 - fac;

  for (j = 0; j < ndim; j++)
    ptry[j] = psum[j] * fac1 - p[ihi][j] * fac2;

  ytry = func(ptry);

  if (ytry < y[ihi]) {
    y[ihi] = ytry;
    for (j = 0; j < ndim; j++) {
      psum[j] = psum[j] - p[ihi][j] + ptry[j];
      p[ihi][j] = ptry[j];
    }
  }

  return (ytry);
}



float testfunc(a)
float *a;
{
  float z;
  float o0, o1, o2, cw, cx, cy;
  float w, x, y, theta_w, theta_y;

  theta_w = 3.1415927 / 4.0;
  theta_y = -3.1415927 / 4.0;

  cw = 3.0;
  cx = 15.0;
  cy = 1.0;
  o0 = -2.0;
  o1 = 3.0;
  o2 = -1.0;


  x = cos(theta_w) * (a[0] - o1) - sin(theta_w) * (a[1] - o2);
  y = sin(theta_w) * (a[0] - o1) + cos(theta_w) * (a[1] - o2);

  x = cos(theta_y) * x - sin(theta_y) * (a[2] - o0);
  w = sin(theta_y) * x + cos(theta_y) * (a[2] - o0);

  w *= cw;
  x *= cx;
  y *= cy;

  z = fabs(w * w * w) + x * x * x * x + y * y * y * y;

  return (z);
}
