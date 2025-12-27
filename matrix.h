#include "stdlib.h"

typedef struct {
  size_t dim;
  double * mat;
} Matrix;

double mat_get(Matrix mat, size_t row, size_t col);

void mat_set(Matrix mat, size_t row, size_t col, double val);

void debug_matrix(Matrix mat);
void swap_rows(Matrix mat, size_t row_i, size_t row_j);

/* Equivalent to Row Operation Ri <- Ri + cRj*/
void row_addition(Matrix mat, size_t row_i, size_t row_j, double c);


