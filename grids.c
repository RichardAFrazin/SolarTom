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
double* rad_bin_boundaries(int bin){
static  double s[2];  // [outer, inner]
  int j;
  double q1, q2, p, drmin = 0.25, drmax = 4.;
  if ((bin < 0) || (bin >= NRAD)){
    fprintf(stderr, "rad_boundaries: invalid bin value.\n");
      exit(-1);
  }
  p = (drmax - drmin)/NRAD;
  q1 = RMIN + drmin;
  q2 = RMIN;
  for (j=1; j <= bin; j++){ // loop not entered if bin=0
    q2 = q1;
    q1 += p*j;
  }
  s[0] = q1; //outer
  s[1] = q2; //inner
  return(s);
}

int rad_bin_number(double dist){  // returns radial bin index given distance
  int bin = -1;
  double ibd, obd, *s;
  for (int b=0; b < NRAD; b++){
    s = rad_bin_boundaries(b);
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

#else

double* rad_bin_boundaries(int bin){
double dr, s[2];  // [outer, inner]
  dr   = (RMAX - ((double) RMIN)) / ((double) NRAD);
  s[1] = ((double) RMIN) + dr*(bin  );
  s[0] = ((double) RMIN) + dr*(bin+1);
  return(s);
}

int rad_bin_number(double dist){
int bin = -1;
  bin = floor( (dist - (double) RMIN)*((double) NRAD) / (RMAX - (double) RMIN));
  return(bin);  
}

#endif

/* make a function that prints the grid values */
void print_grid(void){
  char *gridfile;
  int k;
  double *s, center, size;
  FILE *fp;
  strcpy(gridfile, GRID_FILENAME);
  fp = fopen(gridfile, "w");
  for (k=0; k<NRAD; k++){
    s = rad_bin_boundaries(k);
    center = (*s + *(s+1))/2.;
    size = *s - *(s+1);
    fprintf(fp, "%d  %g  %g\n", k, center, size);    
  }
  fclose(fp);
}
