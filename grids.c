/*
 * grids.c
 *
 * Functions to deal with user-provide non-uniform grids
 *
 * A.M.Vasquez, R.A.Frazin: CLASP Fall-2017.
 *
 */


#include <stdlib.h>
#include <math.h>
#include "headers.h"

#ifdef NONUNIFORMRAD
double* rad_boundaries(int bin){
  double s[2];  // [outer, inner]
  int j;
  double q1, q2, p, drmin = 0.25, drmax = 4.;
  if ((bin < 0) || (bin >= NRAD)){
    fprintf(stderr, "rad_boundaries: invalid bin value.\n")
      exit(-1)
  }
  p = (drmax - drmin)/NRAD;
  q1 = RMIN + drmin;
  q2 = RMIN;
  for (j=1; j <= bin; j++){ // loop not entered if bin=0
    q2 = q1;
    q1 += p*j;
  }
  s[0] = q1;
  s[1] = q2;
  return(s);
}

int rad_bin_number(double dist){  // returns radial bin index given distance
  int bin = -1;
  double ibd, obd, *s;
  for (int b=0; b < NRAD; b++){
    s = rad_bin_boundaries(b);
    ibd = *(s + 1); // inner
    obd = *s; // outer
    if (dist > ibd) && (dist <= obd){
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
double* rad_boundaries(int bin){
}

int rad_bin_number(double dist){

  bin = floor( (distance - (double) RMIN)*((double) NRAD) /
	      (rmax - (double) RMIN));
  
}
#endif


/* make a function that prints the grid values */
void print_grid(void){
  char *gridfile;
  int k;
  double *s, center, size;
  strcpy(gridfile, GRID_FILENAME);
  fp = fopen(gridfile, "w");
  for (k=0; k<NRAD; k++){
    s = rad_boundaries(k);
    center = (*s + *(s+1))/2.;
    size = *s - *(s+1);
    strcpy(gridfile, "%d  %g  %g\n", k, center, size);
    fprintf(fp, gridfile);
  }
  fclose(fp);
}

