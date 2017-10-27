#include <stdio.h>
#include <string.h>
#include "headers.h"

int main(int argc, char **argv) {
  float vA[] = {1, 4, 7, 10, 2, 5, 8, 11, 3, 6, 9, 12};
  int   iA[] = {0, 1, 2,  3, 0, 1, 2,  3, 0, 1, 2,  3};
  int   nA[] = {0, 4, 8, 12};

  float  y[] = {0, 1, 2,  3};
  float delta[] = {0, 0, 0, 0};

  
  struct sparse A;

  A.nf = 4;
  A.ncol = 3;
  A.n = nA;
  A.iB = iA;
  A.vB = vA;
  
  row_extract(A, y, delta, "tre_d", "tre_dd", 1, 2);
  
  return 0;
}
