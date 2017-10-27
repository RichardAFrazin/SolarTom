#include <sys/types.h>
#include <sys/uio.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
int
main( int argc, char **argv)
{
  double x, y[2];
  int k;
  char s[3];

  k = 14;
  sprintf(s,"%02d",k);
  fprintf(stderr,"%s\n",s);

  x = 31.41592653589793;
  y[0] = x;  y[1] = 1./9.;
  /*
  fprintf(stderr,"%g\n%f\n%G\n\n",x,x,x);
  fprintf(stderr,"%16.14g\n%2.14f\n%20.14G\n",x,x,x);
  fprintf(stderr,"%16.14g\n%2.14f\n%20.14G\n",x,x,x);
  fprintf(stderr,"%2.49f, %4.20f\n", y[0], y[1]);
  */

  return;
}
