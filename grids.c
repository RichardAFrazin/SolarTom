/*
 * grids.c
 *
 * Functions to deal with user-provide non-uniform grids
 *
 * A.M.Vasquez: CLASP Fall-2017.
 *
 */

#include "headers.h"

// Compute radialbin from distance.
double radbin(distance)
double distance;
{
  return floor(distance);
}

// Determine radius of cell boundary from index.
double radbound(index)
int index;
{
  return (index*1.0);
}

// Read RAD array from file.
// In this way it will be more flexible, we could handle any radial grid.
double readradialgrid()
int index;
{
  return ();
}
