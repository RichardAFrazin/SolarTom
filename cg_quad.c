#include <sys/types.h>
#include <sys/uio.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
/*#include "/opt/SUNWspro/SC4.2/include/cc/sunmath.h"*/
#include "headers.h"

/* CG_QUAD: - conjugate gradient least-squares minimizer
 *   uses the algorithm in section 10.3 of 
 *      Golub and Van Loan "Matrix Computations" 
 *   x[nc3] object   vector 
 *   lambda[NMATS] vector of regularization params
 *   nsubit  number of iterations done before return
 *   matrices -- structure with lots of stuff 
 *
 *      see cgtest2.m for an easy example
 *
 *    by Richard Frazin, July 2000
 */
int cg_quad(x, lambda, nsubit, matrices)
float *x, *lambda;
int nsubit;
struct matrix matrices[NMATS];
{
  struct sparse sparmat;
  int nc3 = NBINS;		/* number of pixels */
  int i, j, k, fil;
  static float r[NBINS], w[NBINS], p[NBINS];
  float tmp, alpha, beta;
  static float *v1[NMATS], *v2[NMATS];
  float *rho;
  
  for (i = 0; i < nc3; i++) {
    r[i] = 0.0;
    w[i] = 0.0;
    p[i] = 0.0;
  }

  rho = (float *) malloc((nsubit + 1) * sizeof(float));
  if (rho == NULL) {
    fprintf(stderr, "cg_quad: malloc error\n");
    exit(1);
  }
  for (i = 0; i <= nsubit; i++)
    rho[i] = 0.0;

  for (fil = 0; fil < NMATS; fil++) {
    v1[fil] = (float *) malloc(matrices[fil].nf * sizeof(float));
    v2[fil] = (float *) malloc(matrices[fil].nf * sizeof(float));
    if (v1[fil] == NULL || v2[fil] == NULL) {
      fprintf(stderr, "cg_quad: malloc error\n");
      exit(1);
    }
  }


  /* initialize r = B'(y - Bx)
   *     v1 = Bx
   *     v2 = y - v1
   */

  for (j = 0; j < 2; j++) {
    for (fil = 0; fil < NMATS; fil++) {

      sparmat.ncol = nc3;
      sparmat.nf = matrices[fil].nf;
      sparmat.n = matrices[fil].n;
      sparmat.iB = matrices[fil].iB;
      sparmat.vB = matrices[fil].vB;

      if (j == 0) {
        vmultSparse(sparmat, x, v1[fil], lambda[fil]);
      }


      if (j == 1) {
        for (i = 0; i < sparmat.nf; i++) {
          *(v2[fil] + i) =
            *(matrices[fil].y + i) - *(v1[fil] + i);
        }

        vmult_transposeSparse(sparmat, v2[fil], r, lambda[fil], 0);
      }
    }			/* fil loop */
  }				/* j loop */

  rho[0] = 0.0;
  for (i = 0; i < nc3; i++) {
    rho[0] += r[i] * r[i];
  }

  /* initialization complete, start the iterations */

  for (k = 1; k <= nsubit; k++) {

    if (k == 1) {
      for (i = 0; i < nc3; i++)
        p[i] = r[i];
    } else {
      beta = rho[k - 1] / rho[k - 2];
      for (i = 0; i < nc3; i++)
        p[i] = r[i] + beta * p[i];
    }

    /* w = B'Bp , set v1 = Bp */

    for (i = 0; i < nc3; i++)
      w[i] = 0.0;

    for (j = 0; j < 2; j++) {
      for (fil = 0; fil < NMATS; fil++) {

        sparmat.ncol = nc3;
        sparmat.nf = matrices[fil].nf;
        sparmat.n = matrices[fil].n;
        sparmat.iB = matrices[fil].iB;
        sparmat.vB = matrices[fil].vB;

        if (j == 0) {
          vmultSparse(sparmat, p, v1[fil], lambda[fil]);
        }

        if (j == 1) {
          vmult_transposeSparse(sparmat, v1[fil], w, lambda[fil],
                                0);
        }
      }			/* fil loop */
    }			/* j loop */

    tmp = 0.0;
    for (i = 0; i < nc3; i++) {
      tmp += p[i] * w[i];
    }
    alpha = rho[k - 1] / tmp;
    for (i = 0; i < nc3; i++) {
      x[i] += alpha * p[i];
      r[i] -= alpha * w[i];
      rho[k] += r[i] * r[i];
    }



  }				/* k loop */


  free(rho);
  for (fil = 0; fil < NMATS; fil++) {
    free(v1[fil]);
    free(v2[fil]);
  }

  return (0);

}
