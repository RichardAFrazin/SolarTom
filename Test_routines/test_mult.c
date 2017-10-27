#include "headers.h"

int main(int argc, char **argv) {
  float vA[] = {1, 4, 2, 5, 3, 6};
  int   iA[] = {0, 1, 0, 1, 0, 1};
  int   nA[] = {0, 2, 4, 6};

  float vB[] = {1, 3, 5, 2, 4, 6};
  int   iB[] = {0, 1, 2, 0, 1, 2};
  int   nB[] = {0, 3, 6};

  struct sparse A, B;

  A.nf = 2;
  A.ncol = 3;
  A.n = nA;
  A.iB = iA;
  A.vB = vA;

  B.nf = 3;
  B.ncol = 2;
  B.n = nB;
  B.iB = iB;
  B.vB = vB;

  matmultSparse(A,B,"tm");
}
