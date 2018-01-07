#include <stdio.h>
#include <string.h>
#include <signal.h>
#include <unistd.h>
#include <assert.h>
#include <stdlib.h>
#include <stdarg.h>

#include "headers.h"

#if defined BOUNDARY_FILE && defined FESSMIN
static int boundary[NBINS];
#elif defined FESSMIN
static int *boundary = NULL;
#endif
#if defined CONSTRAINT_FILE && defined FESSMIN
static float constraint[NBINS];
#elif defined FESSMIN
static float *constraint = NULL;
#endif

#define MAXCHUNK  536500000 /*largest number of 4-byte objects fread can handle */

/* solve.c  by Richard Frazin 7/99 */

/* n_optarg is the number of optional arguments.  For now, the value must be
 *   0 or 1.  If n_optarg is 0, there is no argument following it in the 
 *   calling sequence and the 0th matrix name comes from FILESTR0  in 
 *   buildA_params.h.  If n_optarg is 1, the next argument is char * pointing
 *   to the name of the 0th matrix suffix.  See also "man stdarg" and 
 *   "./Test_routines/stdarg_test.c"
 */

int solve(float *regularization_parameter, const int *huber_flag,
      float *final_normx, float *final_normd, float *final_normt, 
	  char *x_infile, char *x_outfile, int n_optarg, ...)
{
  FILE *fid_xoutput, *fid_xinput;
  float *x, *x_o;//  float x[NBINS], x_o[NBINS];
#ifdef FESSMIN
  float *x_work;// float x_work[NBINS];
#endif
  float lambda[NMATS];
  int i, jj, k, fil, itcount, status,number_read;
  int chunk, nchunks, lastchunk; /* used for breaking up the fread statements */ 
  float tmp, xoi, xni, changetol;
  float *rr;
  struct matrix matrices[NMATS];
  struct sparse sparmat;
  char filestring[MAXPATH];
  float normd[ITMAX + 2], normx[ITMAX + 2], normt[ITMAX + 2];
  const int nc3 = NBINS;
#ifdef CONJGRAD
  const int nsubit = NSUBIT;
#endif
  va_list pt_optarg; /*optional argument pointer*/
  char *matrix_name;

  x      = (float *)malloc(NBINS*sizeof(float));
  x_o    = (float *)malloc(NBINS*sizeof(float));
#ifdef FESSMIN  
  x_work = (float *)malloc(NBINS*sizeof(float));
#endif

  if (n_optarg == 0){
    strcpy(matrices[0].file_id, FILESTR0);
  } else if (n_optarg == 1) {  
    va_start(pt_optarg, n_optarg);
    matrix_name = va_arg(pt_optarg, char *);
    strcpy(matrices[0].file_id,matrix_name);
    va_end(pt_optarg);
  } else if ((n_optarg > 1) || (n_optarg < 0)) {
    fprintf(stderr,"solve: n_optarg must be 0 or 1.  If 1, the suffix of the 0th matrix must also be passed.");
    fflush(stderr);
    exit(1);
  }

  if (NMATS > 1)
#ifdef FILESTR1
       strcpy(matrices[1].file_id, FILESTR1);
#endif
  if (NMATS > 2)
#ifdef FILESTR2
       strcpy(matrices[2].file_id, FILESTR2);
#endif
  if (NMATS > 3)
#ifdef FILESTR3
       strcpy(matrices[3].file_id, FILESTR3);
#endif
#ifdef FILESTR4
       strcpy(matrices[4].file_id, FILESTR4);
#endif
  if (NMATS > 5){
       fprintf(stderr,"SOLVE: NMATS = %d!\n",NMATS);
	   fflush(stderr);
       exit(1);
  }
  	    
  for (i = 0; i < ITMAX + 2; i++) {
    normd[i] = 0.0;
    normx[i] = 0.0;
    normt[i] = 0.0;
  }

  for (i = 0; i < NMATS; i++) {

    lambda[i] = regularization_parameter[i];
    matrices[i].huber = huber_flag[i];

#ifdef PRINT_FILE_INFO
    fprintf(stdout,"reading *%s\n",matrices[i].file_id);
    fflush(stdout);
#endif

    /* put n, ro (data), delta arrays in memory */
    strcpy(filestring, BINDIR);
    strcat(filestring, "n");
    strcat(filestring, matrices[i].file_id);
    matrices[i].fid_n = fopen(filestring, "rb");
    if (matrices[i].fid_n == NULL) {
      fprintf(stderr, "%s\n", filestring);
      fprintf(stderr, "Input file(s) not found\n");
      exit(2);
    }
    /*
#ifdef PRINT_FILE_INFO
    fprintf(stderr, "SOLVE: reading %s ... ", filestring);
#endif
    */
    if (fread(matrices[i].n, sizeof(int), nc3 + 1, matrices[i].fid_n)
        != nc3 + 1) {
      fprintf(stderr, "solve: problem reading file %s\n",filestring);
      exit(96);
    }
    /*
#ifdef PRINT_FILE_INFO
    fprintf(stderr, "read.\n");
#endif
    */
    fclose(matrices[i].fid_n);

    /* allocate space for i and v matrices */
    matrices[i].iB = (int *) malloc(matrices[i].n[nc3] * sizeof(int));
    matrices[i].vB =
      (float *) malloc(matrices[i].n[nc3] * sizeof(float));
    if (matrices[i].iB == NULL || matrices[i].vB == NULL) {
      fprintf(stderr,
              "solve: malloc error in (iB and/or vB), file = %d\n",
              i);
      fprintf(stderr, "       2 x %d bytes requested\n",
              matrices[i].n[nc3] * sizeof(float));
      /*exit(1); -- DUMP CORE INSTEAD */
      kill(getpid(), SIGBUS);
    }

    /* load i and v matrices */
    strcpy(filestring, BINDIR);
    strcat(filestring, "i");
    strcat(filestring, matrices[i].file_id);
    matrices[i].fid_ind = fopen(filestring, "rb");
    if (matrices[i].fid_ind == NULL) {
      fprintf(stderr, "%s\n", filestring);
      fprintf(stderr, "Input file(s) not found\n");
      exit(2);
    }

    /* this is how it should work, but fread needs smaller chunks
     * b/c is doesn't seem to be 32-bit compliant (as of 10/22/09)
    if (fread
        (matrices[i].iB, sizeof(int), matrices[i].n[nc3],
         matrices[i].fid_ind) != matrices[i].n[nc3]) {
      fprintf(stderr, "SOLVE: problem reading file %s\n",filestring);
      exit(96);
    }
    */

    nchunks   = matrices[i].n[nc3] / MAXCHUNK; /*number of file chunks */
    lastchunk = matrices[i].n[nc3] % MAXCHUNK; /* length of last chunk */

    chunk = 0;
    while (chunk < nchunks){
      if ( fread(matrices[i].iB + chunk*MAXCHUNK, sizeof(int), MAXCHUNK,
	   matrices[i].fid_ind) != MAXCHUNK) {
	      fprintf(stderr, "SOLVE: problem reading file %s\n",filestring);
	      exit(96);
      }
      chunk++;
    }
    if ( fread(matrices[i].iB + chunk*MAXCHUNK, sizeof(int), lastchunk,
	   matrices[i].fid_ind) != lastchunk) {
	      fprintf(stderr, "SOLVE: problem reading file %s\n",filestring);
	      exit(96);
    }
    fclose(matrices[i].fid_ind);

    strcpy(filestring, BINDIR);
    strcat(filestring, "v");
    strcat(filestring, matrices[i].file_id);
    matrices[i].fid_val = fopen(filestring, "rb");
    if (matrices[i].fid_val == NULL) {
      fprintf(stderr, "%s\n", filestring);
      fprintf(stderr, "Input file(s) not found\n");
      exit(2);
    }

    chunk = 0;
    while (chunk < nchunks){
      if ( fread(matrices[i].vB + chunk*MAXCHUNK, sizeof(float), MAXCHUNK,
         matrices[i].fid_val) != MAXCHUNK) {
      fprintf(stderr, "SOLVE: problem reading v file %s\n",filestring);
      exit(96);
      }
      chunk++;
    }
    if ( fread(matrices[i].vB + chunk*MAXCHUNK, sizeof(float), lastchunk,
         matrices[i].fid_val) != lastchunk) {
      fprintf(stderr, "SOLVE: problem reading v file %s\n",filestring);
      exit(96);
    }
    fclose(matrices[i].fid_val);

    /* determine matrices[i].nf */
    jj = 0;
    for (k = 0; k < matrices[i].n[nc3]; k++) {
      if (*(matrices[i].iB + k) > jj)
        jj = *(matrices[i].iB + k);
    }
    matrices[i].nf = jj + 1;


	  
    /* allocate space for y, r, delta vectors */
    matrices[i].r = (float *) malloc(matrices[i].nf * sizeof(float));
    matrices[i].y = (float *) malloc(matrices[i].nf * sizeof(float));
    matrices[i].delta =
      (float *) malloc(matrices[i].nf * sizeof(float));
    if (matrices[i].r == NULL || matrices[i].y == NULL
        || matrices[i].delta == NULL) {
      fprintf(stderr,
              "solve: calloc error in (r, ro and/or delta), i = %d\n",
              i);
      exit(1);
    }

    /* load y and delta vectors */
    strcpy(filestring, BINDIR);
    strcat(filestring, "y");
    strcat(filestring, matrices[i].file_id);
    matrices[i].fid_y = fopen(filestring, "rb");
    if (matrices[i].fid_y == NULL) {
      fprintf(stderr, "%s\n", filestring);
      fprintf(stderr, "Input file(s) not found\n");
      exit(2);
    }

    if (fread(matrices[i].y, sizeof(float), matrices[i].nf,
              matrices[i].fid_y) != matrices[i].nf) {
      fprintf(stderr, "solve: problem reading file %s\n",filestring);
      exit(96);
    }
    fclose(matrices[i].fid_y);
    

    strcpy(filestring, BINDIR);
    strcat(filestring, "delta_");
    strcat(filestring, matrices[i].file_id);
    matrices[i].fid_delta = fopen(filestring, "rb");
    if (matrices[i].fid_delta == NULL) {
      fprintf(stderr, "%s\n", filestring);
      fprintf(stderr, "Input file(s) not found\n");
      exit(2);
    }
    if (fread(matrices[i].delta, sizeof(float), matrices[i].nf,
              matrices[i].fid_delta) != matrices[i].nf) {
      fprintf(stderr, "solve: problem reading file %s\n",filestring);
      exit(96);
    }
    fclose(matrices[i].fid_delta);


    /* rescale delta vector */
    for (jj = 0; jj < matrices[i].nf; jj++)
      matrices[i].delta[jj] *= lambda[i];

  }				/* end of i loop over matrices */


  /* initialize x, x_o */
  strcpy(filestring, BINDIR);
  strncat(filestring, x_infile, MAXPATH);
  fid_xinput = fopen(filestring, "rb");
  if (fid_xinput == NULL) {
    fprintf(stderr, "%s\n", filestring);
    fprintf(stderr, "Input file(s) not found\n");
    exit(2);
  }
 
#ifdef PRINT_FILE_INFO
  fprintf(stderr, "reading %s\n", filestring);
#endif
 
  number_read = fread(x, sizeof(float), nc3, fid_xinput);
  if (number_read != nc3){
    fprintf(stderr, "\n       error reading: %s\n %d elements read %d expected.\n", filestring,number_read,nc3);
    fflush(stderr);
    fclose(fid_xinput);
    assert(0);
  }
  /*
#ifdef PRINT_FILE_INFO
  fprintf(stderr, "read.\n");
#endif
  */
  for (i = 0; i < nc3; i++)
    x_o[i] = x[i];

  /* initialize boundary */

#if defined BOUNDARY_FILE && defined FESSMIN

  strcpy(filestring, BINDIR);
  strcat(filestring, BOUNDARY_FILE);
  fid_xinput = fopen(filestring, "rb");
  if (fid_xinput == NULL) {
    fprintf(stderr, "%s\n", filestring);
    fprintf(stderr, "Input file(s) not found\n");
    exit(2);
  }
#ifdef PRINT_FILE_INFO
  fprintf(stderr, "SOLVE: reading %s ... ", filestring);
#endif

  number_read = fread(boundary, sizeof(int), nc3, fid_xinput);
  if (number_read != nc3){
    fprintf(stderr, "\nerror reading boundary file: %s\n  only read %d elements, expected %d", filestring,number_read,nc3);
    fflush(stdout);
    fclose(fid_xinput);
    assert(0);
  }
  /*
#ifdef PRINT_FILE_INFO
  fprintf(stderr, "read.\n");
#endif
  */
#endif
  
  /* initialize constraint */

#if defined CONSTRAINT_FILE && defined FESSMIN

  strcpy(filestring, BINDIR);
  strcat(filestring, CONSTRAINT_FILE);
  fid_xinput = fopen(filestring, "rb");
  if (fid_xinput == NULL) {
    fprintf(stderr, "%s\n", filestring);
    fprintf(stderr, "Input file(s) not found\n");
    exit(2);
  }
#ifdef PRINT_FILE_INFO
  fprintf(stderr, "SOLVE: reading %s ... ", filestring);
#endif
  assert(constraint != NULL);
  if (fread(constraint, sizeof(float), nc3, fid_xinput) != nc3) {
    fprintf(stderr, "error reading xinputfile: %s\n", filestring);
    exit(2);
  }
  fclose(fid_xinput);
  /*
#ifdef PRINT_FILE_INFO
  fprintf(stderr, "read.\n");
#endif
  */
#endif
  
  /* Calculate r = Bx - y */

  /* zero remainder vector and initialize rr */
  jj = 0;
  for (fil = 0; fil < NMATS; fil++) {
    if (matrices[fil].nf > jj)
      jj = matrices[fil].nf;
    for (i = 0; i < matrices[fil].nf; i++)
      matrices[fil].r[i] = 0.0;
  }

  rr = (float *) malloc(jj * sizeof(float));
  if (rr == NULL) {
    fprintf(stderr, "SOLVE: malloc error (variable rr)!\n");
    exit(0);
  }
  for (i = 0; i < jj; i++)
    rr[i] = 0.0;

  for (fil = 0; fil < NMATS; fil++) {

    sparmat.ncol = nc3;
    sparmat.nf = matrices[fil].nf;
    sparmat.n = matrices[fil].n;
    sparmat.iB = matrices[fil].iB;
    sparmat.vB = matrices[fil].vB;

    vmultSparse(sparmat, x, rr, lambda[fil]);

    for (i = 0; i < matrices[fil].nf; i++)
      matrices[fil].r[i] += (*(rr + i));


  }				/*end fil loop */

  /*
        * matrices[*].y has the data 
        * matrices[*].r has lambda*Ax,
        * rescale by matrices[*].y by lambda and  
        * put lambda*(Ax - y) into matrices[*].r 
        */

  for (fil = 0; fil < NMATS; fil++) {
    for (jj = 0; jj < matrices[fil].nf; jj++) {
      matrices[fil].y[jj] *= lambda[fil];
      matrices[fil].r[jj] -= matrices[fil].y[jj];
    }
  }

  /* calculate normx, normd, normt before iterations */
  for (i = NUMBER_OF_DATA_MATRICES; i < NMATS; i++) {
    for (jj = 0; jj < matrices[i].nf; jj++) {
      tmp = matrices[i].r[jj] / lambda[i];
      normx[0] += tmp * tmp;
      normt[0] += matrices[i].r[jj] * matrices[i].r[jj];
    }
  }
  for (i = 0; i < NUMBER_OF_DATA_MATRICES; i++) {
    for (jj = 0; jj < matrices[i].nf; jj++) {
      tmp = matrices[i].r[jj] / lambda[i];
      normd[0] += tmp * tmp;
      normt[0] += matrices[i].r[jj] * matrices[i].r[jj];
    }
  }

  fprintf(stderr, "initial guess: normx = %g, normd = %g, normt = %g\n",
          normx[0], normd[0], normt[0]);


  itcount = 0;
  status = 0;

  normx[0] = 0.0;
  normd[0] = 0.0;
  normt[0] = 0.0;

  xoi = 0.0;
  xni = 0.0;

  changetol = START_TOL;

  while (status == 0) {

    /* do the optimization */
#ifdef FESSMIN
    fess_hu(x, x_work, lambda, matrices, boundary, constraint);
#elif defined (CONJGRAD)

    cg_quad(x, lambda, nsubit, matrices);
#endif

    /* write to disc every 100 iterations for FESSMIN
     *               every iteration          CONJGRAD
     */
#ifdef FESSMIN

    if (itcount % 100 == 0) {
#elif defined (CONJGRAD)
    if (1) {
#endif
      strcpy(filestring, BINDIR);
      strncat(filestring, x_outfile, MAXPATH);
      fid_xoutput = fopen(filestring, "wb");
      if (fwrite(x, sizeof(float), nc3, fid_xoutput) != nc3) {
        fprintf(stderr, "SOLVE: error in writing output (x)\n");
        exit(5);
      }
      fclose(fid_xoutput);
    }
#ifdef FESSMIN
    if (itcount % 20 == 0) {
#elif defined (CONJGRAD)
    if (1) {
#endif
      /* recalculate  r = Bx - y */
      for (fil = 0; fil < NMATS; fil++) {
        for (i = 0; i < matrices[fil].nf; i++)
          matrices[fil].r[i] = 0.0;
      }
      for (fil = 0; fil < NMATS; fil++) {

        sparmat.ncol = nc3;
        sparmat.nf = matrices[fil].nf;
        sparmat.n = matrices[fil].n;
        sparmat.iB = matrices[fil].iB;
        sparmat.vB = matrices[fil].vB;

        vmultSparse(sparmat, x, rr, lambda[fil]);

        for (i = 0; i < matrices[fil].nf; i++)
          matrices[fil].r[i] += (*(rr + i));


      }			/* end fil loop */

      /* matrices[*].y has the data (scaled by lambda)
       * matrices[*].r has lambda*Ax,
       * put lambda*(Ax - y) into matrices[*].r 
       */

      for (fil = 0; fil < NMATS; fil++) {
        for (jj = 0; jj < matrices[fil].nf; jj++) {
          matrices[fil].r[jj] -= matrices[fil].y[jj];
        }
      }
    }


    /* if loop */
    /* calculate normd, normx, normt */
    for (i = 0; i < NUMBER_OF_DATA_MATRICES; i++) {
      for (jj = 0; jj < matrices[i].nf; jj++) {
        tmp = matrices[i].r[jj] / lambda[i];
        normd[itcount] += tmp * tmp;
        normt[itcount] += matrices[i].r[jj] * matrices[i].r[jj];
      }
    }
    for (i = NUMBER_OF_DATA_MATRICES; i < NMATS; i++) {
      for (jj = 0; jj < matrices[i].nf; jj++) {
        tmp = matrices[i].r[jj] / lambda[i];
        normx[itcount] += tmp * tmp;
        normt[itcount] += matrices[i].r[jj] * matrices[i].r[jj];
      }
    }
#ifdef FESSMIN
    if (itcount % 50 == 0) {
#elif defined (CONJGRAD)
    if (1) {
#endif
      fprintf(stderr,
              "iteration %d: normx = %g, normd = %g, normt = %g\n",
              itcount, normx[itcount], normd[itcount],
              normt[itcount]);
    }


    /* when the fractional change in norm falls below changtol
     *  decrease changetol, calucate mean( abs(x - x_o)/x ).
     *  it it's small enough, exit.  if not, decrease changetol
     *  store x in x_o and start again.
     */


    if (itcount > 1) {

      if ((normt[itcount - 1] - normt[itcount]) / normt[itcount]
          < changetol) {

        xni = 0.0;
        k = 0;
        for (i = 0; i < nc3; i++) {
          if (x[i] > MINIMUM_VALUE) {
            k += 1;
            xni += fabs((x[i] - x_o[i]) / x[i]);
          }
        }
        xni /= (float) k;

        fprintf(stderr,
                "**** mean change %g, changetol %g, fraction counted: %g\n",
                xni, changetol, ((float) k) / ((float) nc3));

        if (xni <= FRACTIONAL_CHANGE_TOL) {
          status = 1;
        } else {
          changetol /= CHANGETOL_FACTOR;
        }

        for (i = 0; i < nc3; i++)
          x_o[i] = x[i];

      }
    }


    if (itcount > ITMAX) {
      status = 2;
      fprintf(stderr,
              "SOLVE: maximum number of iterations exceeded!!!\n");
    }

    itcount++;
    }				/* end while loop */

  *final_normx = normx[itcount - 1];
  *final_normd = normd[itcount - 1];
  *final_normt = normt[itcount - 1];
  
  for (i = 0; i < NMATS; i++) {
    free(matrices[i].r);
    free(matrices[i].y);
    free(matrices[i].delta);
    free(matrices[i].iB);
    free(matrices[i].vB);
  }
  free(rr);

  /* output */
  strcpy(filestring, BINDIR);
  strncat(filestring, x_outfile, MAXPATH);
  fid_xoutput = fopen(filestring, "wb");
  if (fwrite(x, sizeof(float), nc3, fid_xoutput) != nc3) {
    fprintf(stderr, "SOLVE: error in writing output (x)\n");
    exit(5);
  }
  fprintf(stderr, "wrote %s\n", filestring);
  fclose(fid_xoutput);

  free(x);
  free(x_o);
#ifdef FESSMIN
  free(x_work);
#endif
  
  return (status);
}
