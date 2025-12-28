#include "stdlib.h"
#include "complex_polynomial.h"

typedef union {
  Polynomial polynomial;
  double     scalar;
} pmat_item;

typedef struct {
  size_t dim;
  double * mat;
} Matrix;

typedef struct {
  size_t dim;
  pmat_item * mat;
} PolyMatrix;

double mat_get(Matrix mat, size_t row, size_t col);
void mat_set(Matrix mat, size_t row, size_t col, double val);

pmat_item pmat_get(PolyMatrix mat, size_t row, size_t col);
void pmat_set(PolyMatrix mat, size_t row, size_t col, pmat_item val);

void debug_matrix(Matrix mat);
void swap_rows(Matrix mat, size_t row_i, size_t row_j);

/* Equivalent to Row Operation Ri <- Ri + cRj*/
void row_addition(Matrix mat, size_t row_i, size_t row_j, double c);


