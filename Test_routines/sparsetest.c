#include <string.h>
#include "headers.h"

int main(int arc, char **argv) {
  char fs[] = "crap";
  struct sparse a;
  struct sparse b;
  int na[] = { 0, 2, 3, 4, 7 };
  int ia[] = { 0, 2, 2, 1, 0, 1, 3 };
  float va[] = { 1., 2., 1., 1., 2., 3., 1. };
  int nb[] = { 0, 2, 4, 6, 7 };
  int ib[] = { 0, 2, 1, 3, 0, 2, 3 };
  float vb[] = { 3., 1., 2., 1., 2., 1., 1. };

  a.nf = 4;
  b.nf = 4;
  a.ncol = 4;
  b.ncol = 4;

  a.n = na;
  a.iB = ia;
  a.vB = va;
  b.n = nb;
  b.iB = ib;
  b.vB = vb;

  matmultSparse(a, b, fs);

  return (0);
}
