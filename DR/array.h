#ifndef ARRAY_H
#define ARRAY_H

#define INITAL_ARRAY_ROWS 64
#define INITAL_ARRAY_COLS 64

typedef float ARRAY_TYPE;

typedef struct {
  unsigned int rows;
  unsigned int cols;
  unsigned int rows_a;
  unsigned int cols_a;
  ARRAY_TYPE *a;
}
ARRAY;


typedef void (*DOALL_FPTR)(const ARRAY *A, ARRAY_TYPE *v, int r, int c);

ARRAY *array_new();
ARRAY *array_new_s(const unsigned int r, const unsigned int c);
void array_delete(ARRAY **A);
ARRAY_TYPE *array_get(const ARRAY *A, const unsigned int r, const unsigned int c);
void array_put(ARRAY **A, const unsigned int r, const unsigned int c, const ARRAY_TYPE v);
void array_doall(const ARRAY *A, const DOALL_FPTR fptr);

void array_print(const ARRAY *A, ARRAY_TYPE *v, int r, int c);

#endif
