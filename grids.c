/*
 * grids.c
 *
 * Functions to deal with user-provide non-uniform grids
 *
 * A.M.Vasquez, R.A.Frazin: CLASP Fall-2017.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "buildA_params.h"


#ifdef NONUNIFORMRAD

void rad_bin_boundaries(int bin, double *s){ // s is bin boundary array = [outer, inner]
  int j;
  double q1, q2, p, drmin = 0.25, drmax = 4.0;
  if ((bin < 0) || (bin >= NRAD)){
    fprintf(stderr, "rad_boundaries: invalid bin value.\n");
      exit(-1);
  }
  p = (drmax - drmin)/(NRAD-1);
  q1 = RMIN;
  q2 = RMIN + drmin;
  for (j=1; j <= bin; j++){ // loop not entered if bin=0
    q1 = q2;
    q2 = q1 + drmin + p*j;
  }
  s[0] = q2; //outer
  s[1] = q1; //inner
  return;
}

int rad_bin_number(double dist){  // Returns the radial bin index of given distance.
  int bin = -1;
  double ibd, obd, s[2];
  // fprintf(stderr,"Test inside rad_bin_number. dist = %g, NRAD = %d.\n",dist,NRAD);
  for (int b=0; b < NRAD; b++){
    rad_bin_boundaries(b, s);
    ibd = *(s + 1); // inner
    obd = *s; // outer
    if ((dist > ibd) && (dist <= obd)){
	bin = b;
	break;
    }
  }
  if (dist == RMIN)
    bin = 0;
  if (bin == -1){
    fprintf(stderr, "rad_bin_number: invalid bin.\n");
    exit(-1);
  }
  return(bin);
}  

#else // if UNIFORM-RADIAL-GRID then:

void rad_bin_boundaries(int bin, double *s){  // s is bin boundary array = [outer, inner]
  if ((bin < 0) || (bin >= NRAD)){
    fprintf(stderr, "rad_boundaries: invalid bin value.\n");
      exit(-1);
  }
  double dr;
  dr   = (RMAX - ((double) RMIN)) / ((double) NRAD);
  s[0] = ((double) RMIN) + dr*(bin+1); //outer
  s[1] = ((double) RMIN) + dr*(bin  ); //inner
  //fprintf(stderr,"Test inside rad_bin_boundaries. Rmin = %g, Rmax = %g, NRAD = %d.\n",RMIN,RMAX,NRAD);
  //fprintf(stderr,"Test inside rad_bin_boundaries. bin = %d, r = %g, dr = %g.\n",bin,(s[0]+s[1])/2.,dr);
  return;
}

int rad_bin_number(double dist){
  int bin = -1;
  // fprintf(stderr,"Test inside rad_bin_number. dist = %g, Rmin = %g, Rmax = %g, NRAD = %d.\n",dist,RMIN,RMAX,NRAD);
  bin = floor( (dist - (double) RMIN) * ((double) NRAD)/(RMAX-(double) RMIN) );
  return(bin);  
}

#endif

void print_grid(void){
  char gridfile[256];
  int k;
  double s[2], center, size;
  FILE *fp;
  strcpy(gridfile, GRID_FILENAME);
  fp = fopen(gridfile, "w");
  for (k=0; k<NRAD; k++){
    rad_bin_boundaries(k, s);
    center = (*s + *(s+1))/2.;
    size = *s - *(s+1);
    fprintf(fp, "%d  %g  %g\n", k, center, size);    
  }
  fclose(fp);
}

double GridDivision(double num, double denom){
  double test_ratio, max_time = 1000.;
  if (denom == 0.)
    return(copysign(max_time,num));
  if (fabs(denom) < fabs(num)){
       test_ratio = fabs(denom)/fabs(num);
       if (test_ratio < 1./max_time){
	 return(copysign(max_time,num*denom));
       }
  } 
  return(num/denom);  
}
