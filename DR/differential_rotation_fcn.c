/* differential_rotation_fcn_2.c 
 *
 *   by Richard Frazin and Alberto Vasquez Oct. 2012
 *
 * Given radius (Rs), latitude (radians), time (days) and 
 *    a reference time (days), this function outputs
 *    the longitudinal displacement due to differential
 *    rotation.  Output is in radians.
 *
 * This version (2) mimics the SOLID (r=1) and DOTTED (r=1.6) 
 * curves in Fig-8 of Salvatore Mancuso and Silvio Giordano
 * The Astrophysical Journal, 729:79 (8pp), 2011.
 *
 * The parameter w is the "mixing scale height" between the
 * r=1 and r=1.6 coundary conditions.
 *
 */

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>

double differential_rotation_fcn(double r, double theta, double mjd,
				   double mjd_ref, double *w){
  double delta_phi;
  double A, B, C, B1, C1, Tau_sid, Omega;

  if (r < 1.){
    fprintf(stderr,"input radius (r) is incorrect, r must be >=1.");
    exit(1);
  }

  A = +14.713 ; /*[deg/day]*/
  B = - 2.396;
  C = - 1.787;

  B1 = -0.5;    /*[deg/day]*/
  C1 = -0.5;

  Tau_sid = 25.38 ; /* days*/

  /* Differential rotation angular velocity in [rad/day] */
  Omega = ( A + B  * pow(sin(theta),2) + C  * pow(sin(theta),4) ) * (M_PI/180.) * (   exp(-(r-1)/ w[0]) ) + ( A + B1 * pow(sin(theta),2) + C1 * pow(sin(theta),4) ) * (M_PI/180.) * ( 1-exp(-(r-1)/w[0] ) ) ;

  delta_phi = (Omega-2.*M_PI/Tau_sid) * (mjd-mjd_ref);

  /*
  fprintf(stderr,"DRF: r = %f, th = %f, w = %f, delta_phi = %f     ",
	  r,theta*180./M_PI,w[0],delta_phi*180./M_PI);
  fprintf(stderr,"mjd = %f, mjd_ref = %f, deltaT (d) = %f\n",
	  mjd,mjd_ref,mjd-mjd_ref); 
  fflush(stderr);
  */

  return delta_phi;
}
