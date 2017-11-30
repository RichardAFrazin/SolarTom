/*
 * grids.c
 *
 * Functions to deal with user-provide non-uniform grids
 *
 * A.M.Vasquez: CLASP Fall-2017.
 *
 */

#include "headers.h"

// Return the radial grid bin correspoding to given a distance.
double radbin(distance)
double distance;
int    bin;
{
  readradialgrid(....);
  
  return bin;
}

// Return radius of radial grid cell boundary given an index.
double radbound(index)
int index;
double bound;
{
  return bound;
}

// Read RAD array from file.
// In this way it will be more flexible, we could handle any radial grid.
void readradialgrid(grid_filename,NRAD,rad,drad)
double rad[NRAD],drad[NRAD];
{
  // read grid_filename to load arrays: rad,drad.
}
