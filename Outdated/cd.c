
/* cd.c, linmin.c, testfunc.c 
 *
 *    by Richard Frazin 9/99
 *
 * cd.c is a primative coordinate descent algorithm.
 * it calls linmin.c to do the line minimization
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


#define FUNCTION_NAME cvcalc
#define MINIMUM_VALUE (1.e-10)  /*minimum parameter value */

#define NDIMS (3)             /* # of dimensions*/
#define NCYCLES (30)           /* # of passes thru coords */
#define UPDATE_FRAC (.005)     /* continuing delta for search (see code) */
#define STOP_TOL (.001)        /* stopping criterion on objective function */
#define EXPAN_FAC (4.0)       /* step expansion factor */
#define MAXHIST (20)         /* # max number of points in line search */


float FUNCTION_NAME(float *);
void linmin( float *, float *, float *, float *, float *, float *, 
	    float (*)(float *));

void
main( int argc, char **argv)
{
float initial_guess[] = {4.e-6,4.e-6,2.e-5} ;
float x_o[NDIMS], x_n[NDIMS], x_c[NDIMS],  x_d[NDIMS], v[NDIMS];
float f_o, f_n, f_tmp, q; 
int   ncycles;
int   i, j, k;
float xhist[NDIMS][NCYCLES+1], fhist[NCYCLES+1];
float (*func)(float *);

 ncycles = NCYCLES;
 func = FUNCTION_NAME; 

 for (i = 0; i < NDIMS; i++){
   x_o[i] = initial_guess[i];
   x_n[i] = x_o[i];
   x_c[i] = 0.0;
   x_d[i] = 0.0;
 }
 f_o = func(x_o);
 f_n = f_o;

 for (i = 0; i < NDIMS; i++){
   xhist[i][0] = initial_guess[i];
 }
 fhist[0] = f_n;

 i = -1;
 k = - 1;
 fprintf(stdout,"\n[%d,%d]: x_n = [",i,k);
 for (j = 0; j < NDIMS-1; j++)
   fprintf(stdout,"%g,",x_n[j]);
 fprintf(stdout,"%g] f_n = %g\n",x_n[NDIMS-1],f_n);

 for (i = 1; i <= ncycles; i++){

   /* check for stopping criterion */
   if (i > 0){
     q = 0.0;
     for (j = 0; j < NDIMS; j++){
       f_tmp = fabs( (xhist[j][i-1]-xhist[j][i])/xhist[j][i]  );
      if (f_tmp > q)
	q = f_tmp;
     }
     if (q <= STOP_TOL){
       fprintf(stdout,"final result:\n");
       for (j = 0; j < NDIMS; j++)
	 fprintf(stdout,"%d: %g\n",j,xhist[j][i]);
       exit(0);
     }
   }

   for (j = 0; j < NDIMS; j++){
     x_c[j] = x_o[j];
   }

   for (k = 0; k < NDIMS; k++){

     /* set direction vector v */

     for (j = 0; j < NDIMS; j++)
       v[j] = 0.0;

     v[k] = 1.0;

     /* do the line search along the coordinate
      *   and try to get close to the minimum
      */ 

     linmin(x_d, x_o, &f_o, v, x_n, &f_n, func);

     for (j = 0; j < NDIMS; j++){
       x_o[j] = x_n[j];
       f_o = f_n;
     }

     fprintf(stdout,"[%d,%d]: x_n = [",i,k);
     for (j = 0; j < NDIMS-1; j++)
       fprintf(stdout,"%g,",x_n[j]);
     fprintf(stdout,"%g] f_n = %g\n",x_n[NDIMS-1],f_n);
    
   } /* loop over k */

   for (j = 0; j < NDIMS; j++){
     x_d[j] = x_c[j];
     x_c[j] = x_o[j];
     xhist[j][i+1] = x_o[j];
   }
   fhist[i+1] = f_o;


     
 } /* loop over cycles */

}

void
linmin( x_comp, x_o, f_o, v, x_n, f_n, func)
     float *x_comp, *x_o, *x_n, *v; /* these are vectors */
     float *f_o, *f_n; /* these are scalars */
     float (*func)(float *); /* this is the objective function */
{
  int   i, j, histcount;
  float len, tmp, f_tmp;
  float x_tmp[NDIMS];
  float valhist[MAXHIST], lenhist[MAXHIST];
  float len1, len2, len3, val1, val2, val3;
  float alpha, beta, gamma, f_pred; /* f_pred = alpha*(x - beta)^2 + gamma */
  int status;
#ifdef MINIMUM_VALUE
  float len_min, len_max;
#endif

#ifdef MINIMUM_VALUE
  /* check to see if initial guess is in bounds */
  status = 0;
  for (i = 0; i < NDIMS; i++){
    if (x_o[i] < MINIMUM_VALUE){
      status += 1;
      fprintf(stderr,"LINMIN: invalid input vector!  ");
      fprintf(stderr,"component %d: value %g\n",i,x_o[i]);
    }
  }
  if (status > 1)
    exit(status);

  /* find min,max values of len allowed 
   *   len_min < 0
   *   len_max > 0
   */
 
  len_min =  1.0;
  len_max = -1.0;

  for (i = 0; i < NDIMS; i++){
    if ( fabs(v[i]) > 1.e-8/NDIMS ){
      tmp = (MINIMUM_VALUE - x_o[i])/v[i];
      if ( tmp < 0 ){
	if ( tmp < len_min)
	  len_min = tmp;
      }	else {
	if (tmp > len_max)
	  len_max = tmp;
      }	
    }
  }
  if (len_min > 0.0)
    len_min = - 1.0e33;
  if (len_max < 0.0)
    len_max =   1.0e33;
#endif

  histcount = 0;
  valhist[0] = *f_o;
  lenhist[0] = 0.0;

  tmp = 0.0;
  for (i = 0; i < NDIMS; i++)
    tmp += (x_comp[i] - x_o[i])*v[i];

  if (tmp == 0.0){
    fprintf(stderr,"LINMIN: warning: tmp = 0.  Resetting.");
    for (i = 0; i < NDIMS; i++)
      tmp += x_o[i]*x_o[i]*v[i]*v[i];
    tmp = .1*sqrt(tmp)*UPDATE_FRAC;
  }
    
  len = tmp*UPDATE_FRAC;

#ifdef MINIMUM_VALUE
  if (len < len_min)
    len = len_min*.9999;
  if (len > len_max)
    len = len_max*.9999;
#endif

  for (i = 0; i < NDIMS; i++)
    x_tmp[i] = x_o[i] + len*v[i] ;

  f_tmp = func(x_tmp);
  
  histcount = 1;
  valhist[1] = f_tmp;
  lenhist[1] = len;

  /* except for the last pt., valhist should be
   *  a decreasing sequence
   */

  if (valhist[1] > valhist[0]){

    len = - EXPAN_FAC*len;

    len1 = lenhist[0];
    len2 = lenhist[1];
    val1 = valhist[0];
    val2 = valhist[1];
    lenhist[1] = len1;
    lenhist[0] = len2;
    valhist[1] = val1;
    valhist[0] = val2;

  } else {
    len = EXPAN_FAC*len;
  }

     /* the idea is to find 3 pts s.t.
      * f(b) < f(a) && f(b) < f(c)
      *   for a < b < c or a > b > c
      * valhist should be a monotonically
      *   decreasing sequence except for
      *   the last point.  Thus, the last
      *   three points are used for parabolic
      *   interpolation.
      */

  status = 1; /* exponential stepping */
  while ( (valhist[histcount] < valhist[histcount-1]) &&
	  (status > 0) && (histcount < MAXHIST - 1)){

    histcount++;
    if (status == 1) /* exponential stepping */
      len *= EXPAN_FAC;
    if (status == 2) /* linear stepping */
      len = 2.0*lenhist[histcount-1] - lenhist[histcount-2];

#ifdef MINIMUM_VALUE
    if (len < len_min){
      len = len_min*.9999;
      status = 0;
    }
    if (len > len_max){
      len = len_max*.9999;
      status = 0;
    }
#endif

    for (i = 0; i < NDIMS; i++)
      x_tmp[i] = x_o[i] + len*v[i] ;

    f_tmp = func(x_tmp);

    /* check to see if you stepped too far
     * if so, go into linear stepping mode
     */
    if ( (i > 2) && (f_tmp > valhist[histcount-1]) ){

      /* fit parabola to previous 3 data points */

      len1 = lenhist[histcount-3];
      len2 = lenhist[histcount-2];
      len3 = lenhist[histcount-1];
      val1 = valhist[histcount-3];
      val2 = valhist[histcount-2];
      val3 = valhist[histcount-1];

      beta = val1*(len3*len3 - len2*len2) + val2*(len1*len1 - len3*len3) 
	+ val3*(len2*len2 - len1*len1);
      beta /= val1*(len3 - len2) + val2*(len1-len3) + val3*(len2-len1); 
      beta /= 2.0;

      alpha = (val2 - val3);
      alpha /= (len2 - beta)*(len2 - beta) - (len3 - beta)*(len3 - beta);

      gamma = val3 - alpha*(len3 - beta)*(len3 - beta);

      f_pred = alpha*(len - beta)*(len - beta) + gamma; 

      if ( f_tmp > f_pred ){ 

	status = 2; /* go to linear stepping */

	len = 2.0*lenhist[histcount-1] - lenhist[histcount-2];
	for (i = 0; i < NDIMS; i++)
	  x_tmp[i] = x_o[i] + len*v[i] ;
	f_tmp = func(x_tmp);
      }
    }

    lenhist[histcount] = len;
    valhist[histcount] = f_tmp;

  }  /* while loop */

  /* parabolic interpolation */

  len1 = lenhist[histcount-2];
  len2 = lenhist[histcount-1];
  len3 = lenhist[histcount];

  val1 = valhist[histcount-2];
  val2 = valhist[histcount-1];
  val3 = valhist[histcount];

  len = val1*(len3*len3 - len2*len2) + val2*(len1*len1 - len3*len3) 
    + val3*(len2*len2 - len1*len1);
  len /= val1*(len3 - len2) + val2*(len1-len3) + val3*(len2-len1); 
  len /= 2.0;

#ifdef MINIMUM_VALUE
  if (len < len_min)
       len = len_min*.9999;
  if (len > len_max)
       len = len_max*.9999;
#endif

  for (i = 0; i < NDIMS; i++)
    x_tmp[i] = x_o[i] + len*v[i] ;

  histcount++;
  f_tmp = func(x_tmp);
  lenhist[histcount] = len;
  valhist[histcount] = f_tmp;

  j = histcount;
  for (i = 0; i <= histcount; i++){
    if (valhist[i] < f_tmp){
      f_tmp = valhist[i];
      len    = lenhist[i];
      j = i;
    }
  }

  if (j != histcount)
    fprintf(stderr,"LINMIN: warning: parabolic interpolation failed!\n");

  /* return the best point so far */
  for (i = 0; i < NDIMS; i++)
    x_n[i] = x_o[i] + len*v[i] ;
  
  *f_n = f_tmp;

}


float
testfunc( a )
float *a;
{
float z;
float o0,o1, o2,cw,cx,cy;
float w,x,y, theta_w, theta_y;

theta_w =   3.1415927/4.0;
theta_y = - 3.1415927/4.0;

cw =  3.0;
cx =  9.0;
cy =  1.0;
o0 = -2.0;
o1 =  3.0;
o2 = -1.0;


x = cos(theta_w)*(a[0]-o1) - sin(theta_w)*(a[1]-o2);
y = sin(theta_w)*(a[0]-o1) + cos(theta_w)*(a[1]-o2);

x = cos(theta_y)*x - sin(theta_y)*(a[2]-o0);
w = sin(theta_y)*x + cos(theta_y)*(a[2]-o0);

w *= cw;
x *= cx;
y *= cy;

z = fabs(w*w*w) + x*x*x*x + y*y*y*y;

return(z);
}
