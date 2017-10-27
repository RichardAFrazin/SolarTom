#include <stdio.h>
#include <assert.h>

#include "array.h"


void initialize(const ARRAY *A, ARRAY_TYPE *v, int r, int c) {
  *v = r * (*A).cols + c;
}


int main(int argc, char **argv) {
  ARRAY *A;
  
  A = array_new_s(4, 4);

  array_doall(A, &initialize);

  array_doall(A, &array_print);
  
  array_put(&A, 10, 10, 69);
  printf("\n\n");
  array_doall(A, &array_print);

  array_put(&A, 17, 17, 69);
  printf("\n\n");
  array_doall(A, &array_print);
  
  array_delete(&A);

  assert(A == NULL);
    
  return 0;
}
