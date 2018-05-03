/*	$Id: rots.c,v 1.1.1.1 2008/05/02 18:13:10 rfrazin Exp $	*/
/*
 * rots.c : do rotations
 *
 * replaces matlab files:
 *	rotx.m
 *	roty.m
 *	rotz.m
 */

#ifdef OPENBSD
#include <err.h>
#endif
/*#include <malloc.h>*/
#include <stdlib.h>
#include <math.h>
#include "headers.h"

Rot *rotx(ax)
double ax;
{
  Rot *foo;

  if ((foo = (Rot *) malloc(sizeof(Rot))) == NULL) {
    printf("malloc in rotx");
    exit(1);
  }
  (*foo)[0][0] = 1.0;
  (*foo)[0][1] = 0.0;
  (*foo)[0][2] = 0.0;
  (*foo)[1][0] = 0.0;
  (*foo)[1][1] = cos(ax);
  (*foo)[1][2] = -sin(ax);
  (*foo)[2][0] = 0.0;
  (*foo)[2][1] = sin(ax);
  (*foo)[2][2] = cos(ax);
  return (foo);
}

Rot *roty(ay)
double ay;
{
  Rot *foo;

  if ((foo = (Rot *) malloc(sizeof(Rot))) == NULL) {
    printf("malloc in roty");
    exit(1);
  }
  (*foo)[0][0] = cos(ay);
  (*foo)[0][1] = 0.0;
  (*foo)[0][2] = -sin(ay);
  (*foo)[1][0] = 0.0;
  (*foo)[1][1] = 1.0;
  (*foo)[1][2] = 0.0;
  (*foo)[2][0] = sin(ay);
  (*foo)[2][1] = 0.0;
  (*foo)[2][2] = cos(ay);
  return (foo);
}

Rot *rotz(az)
double az;
{
  Rot *foo;

  if ((foo = (Rot *) malloc(sizeof(Rot))) == NULL) {
    printf("malloc in rotz");
    exit(1);
  }
  (*foo)[0][0] = cos(az);
  (*foo)[0][1] = -sin(az);
  (*foo)[0][2] = 0.0;
  (*foo)[1][0] = sin(az);
  (*foo)[1][1] = cos(az);
  (*foo)[1][2] = 0.0;
  (*foo)[2][0] = 0.0;
  (*foo)[2][1] = 0.0;
  (*foo)[2][2] = 1.0;
  return (foo);
}


/* Multiply rotation matrices: a = b * c */
void rotmul(a, b, c)
Rot *a, *b, *c;
{
  int i, j, k;
  double tmp;

  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      tmp = 0.0;
      for (k = 0; k < 3; k++) {
        tmp += (*b)[i][k] * (*c)[k][j];
      }
      (*a)[i][j] = tmp;
    }
  }
}

/* Multiply rotation matrix times vector: a = b * c, where it is assumed
   a and c are 3-d double arrays */
void rotvmul(a, b, c)
double *a;
Rot *b;
double *c;
{
  int i, k;
  double tmp;

  for (i = 0; i < 3; i++) {
    tmp = 0.0;
    for (k = 0; k < 3; k++) {
      tmp += (*b)[i][k] * (*(c + k));
    }
    *(a + i) = tmp;
  }
}
