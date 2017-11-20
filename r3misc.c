/*	$Id: r3misc.c,v 1.1.1.1 2008/05/02 18:13:10 rfrazin Exp $	*/

/*
 * r3misc.c
 *
 * Just some convenience functions
 *
 * A.M.Vasquez did a few cosmetic edits. CLASP Fall-2017.
 *
 */

#include "headers.h"

/* Normalize a 3-d vector. */
void r3norm(foo)
double *foo;
{
  double tmp;

  tmp = r3dot(foo, foo);
  tmp = 1.0 / sqrt(tmp);
  r3scalmul(foo, tmp);
}

/* More convenience */
double r3dot(foo, bar)
double *foo, *bar;
{
return (*(foo) * *(bar) + *(foo + 1) * *(bar + 1) + *(foo + 2) * *(bar + 2));
}

void r3scalmul(foo, scal)
double *foo, scal;
{
  *(foo    ) *= scal;
  *(foo + 1) *= scal;
  *(foo + 2) *= scal;
}

/* a = b x c */
void r3cross(a, b, c)
double *a, *b, *c;
{
  *(a    ) = *(b + 1) * *(c + 2) - *(b + 2) * *(c + 1);
  *(a + 1) = *(b + 2) * *(c    ) - *(b    ) * *(c + 2);
  *(a + 2) = *(b    ) * *(c + 1) - *(b + 1) * *(c    );
}

/* a = b */
void r3eq(a, b)
double *a, *b;
{
  *(a    ) = *(b    );
  *(a + 1) = *(b + 1);
  *(a + 2) = *(b + 2);
}

/* a = b + c */
void r3add(a, b, c)
double *a, *b, *c;
{
  *(a    ) = *(b    ) + *(c    );
  *(a + 1) = *(b + 1) + *(c + 1);
  *(a + 2) = *(b + 2) + *(c + 2);
}

int doublecompare(const void *x, const void *y) {
  if ( *((double *) x) > *((double *) y) ) {
    return(1);
  } else if ( *((double *) x) < *((double*) y)) {
    return(-1);
  } else {
    return(0);
  }
}
