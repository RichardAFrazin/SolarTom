/* differential_rotation_fcn.c 
 *
 *   by Richard Frazin and Alberto Vasquez Oct. 2012
 *
 * Given radius (Rs), latitude (radians), time (days) and 
 *    a reference time (days), this function outputs
 *    the longitudinal displacement due to differential
 *    rotation.   
 */

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>

double differential_rotation_fcn(double r, double theta, double mjd,
				 double mjd_ref, double *w){
  double delta_phi;
  double A, B, C, w1, w2, Tau_sid, Omega;

  if (r < 1.){
    fprintf(stderr,"input radius (r) is incorrect, r must be >=1.");
    exit(1);
  }

  w1 = w[0] + w[1] * abs(theta) + w[2] * theta * theta;
  w2 = w[3] + w[4] * abs(theta) + w[5] * theta * theta;

  A = +14.713 ; /*[deg/day]*/
  B = - 2.396;
  C = - 1.787;

  Tau_sid = 25.38 ; /* days*/

  /* Differential rotation angular velocity in [rad/day] */
  Omega = ( A + B * pow(sin(theta),2) + C * pow(sin(theta),4) ) * (M_PI/180.) + w1 * (r-1) + w2 * pow((r-1),2);

  delta_phi = (Omega-2.*M_PI/Tau_sid) * (mjd-mjd_ref);


  return delta_phi;
}
