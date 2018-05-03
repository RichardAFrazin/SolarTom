#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>

#include "array.h"

#define MIN(X,Y) X < Y ? X : Y
#define MAX(X,Y) X > Y ? X : Y

#define LOG2(x) log(x)/log(2)

ARRAY *array_new() {
  ARRAY *A;

  A = (ARRAY *) malloc(sizeof(ARRAY));

  (*A).a = (ARRAY_TYPE *) calloc(INITAL_ARRAY_ROWS * INITAL_ARRAY_COLS, sizeof(ARRAY_TYPE));

  (*A).rows = INITAL_ARRAY_ROWS;
  (*A).cols = INITAL_ARRAY_COLS;

  (*A).rows_a = INITAL_ARRAY_ROWS;
  (*A).cols_a = INITAL_ARRAY_COLS;

  return A;
}

ARRAY *array_new_s(const unsigned int r, const unsigned int c) {
  ARRAY *A;
  unsigned int l_r, l_c;
  unsigned int p_r, p_c;
  
  A = (ARRAY *) malloc(sizeof(ARRAY));

  l_r = ceil(LOG2(r));
  l_c = ceil(LOG2(c));

  p_r = floor(pow(2, l_r + 1));
  p_c = floor(pow(2, l_c + 1));
  
  (*A).a = (ARRAY_TYPE *) calloc(p_r * p_c, sizeof(ARRAY_TYPE));

  (*A).rows = r;
  (*A).cols = c;

  (*A).rows_a = p_r;
  (*A).cols_a = p_c;
  
  return A;
}

void array_delete(ARRAY **A) {
  free((**A).a);
  free(*A);

  *A = NULL;
}

ARRAY_TYPE *array_get(const ARRAY *A, const unsigned int r, const unsigned int c) {
  assert(r <= (*A).rows);
  assert(c <= (*A).cols);

  return ((*A).a + r * (*A).cols_a + c);
}

ARRAY_TYPE *array_get_row_ptr(const ARRAY *A, const unsigned int r) {
  return ((*A).a + r * (*A).cols_a);
}

void array_put(ARRAY **A, const unsigned int r, const unsigned int c, const ARRAY_TYPE v) {
  ARRAY *B = NULL;
  unsigned int min_r, min_c, max_r, max_c;
  int i;

  if (r > (**A).rows || c > (**A).cols) {
    if (r <= (**A).rows_a && c <= (**A).cols_a) {
      if (r > (**A).rows)
        (**A).rows = r;

      if (c > (**A).cols)
          (**A).cols = c;
    }
    else {
      min_r = MIN(r,(**A).rows);
      min_c = MIN(c,(**A).cols);
      
      max_r = MAX(r,(**A).rows);
      max_c = MAX(c,(**A).cols);

      B = array_new_s(max_r, max_c);
      
      for (i = 0; i < min_r; i++) {
        assert(memcpy(array_get_row_ptr(B, i), array_get_row_ptr(*A, i), sizeof(ARRAY_TYPE) * min_c) != NULL);
      }
      
      array_delete(A);
      
      *A = B;
    }
  }

  printf("%d %d %d %d\n", r, c, (**A).rows, (**A).cols);
  *array_get(*A, r, c) = v;
}

void array_doall(const ARRAY *A, const DOALL_FPTR fptr) {
  ARRAY_TYPE *a_ptr;
  int i, j;

  printf("%d %d\n", (*A).rows, (*A).cols);
  printf("%d %d\n", (*A).rows_a, (*A).cols_a);
  
  for (i = 0; i <= (*A).rows; i++) {
    for (j = 0; j <= (*A).cols; j++) {
      a_ptr = array_get(A, i, j);
      (*fptr)(A, a_ptr, i, j);
    }
  }
}

void array_print(const ARRAY *A, ARRAY_TYPE *v, int r, int c) {
  printf("%f ", *v);

  if(c == (*A).cols) {
    printf("\n");
  }
}
