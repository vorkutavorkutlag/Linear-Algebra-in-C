#include "matrix.h"
#include "stdio.h"

double mat_get(Matrix mat, size_t row, size_t col) {
  return mat.mat[mat.dim * row + col];
}

void mat_set(Matrix mat, size_t row, size_t col, double val) {
  mat.mat[mat.dim * row + col] = val;
}

void debug_matrix(Matrix mat) {
  for (size_t row = 0; row < mat.dim; row++) {
    if (row) printf("\n");
      for (size_t col = 0; col < mat.dim; col++) {
        printf("%g ", mat_get(mat, row, col));
      }
  }
  printf("\n");
}

void swap_rows(Matrix mat, size_t row_i, size_t row_j) {
  for (size_t index = 0; index < mat.dim; index++) {
    double temp = mat_get(mat, row_i, index);
    mat_set(mat, row_i, index, mat_get(mat, row_j, index));
    mat_set(mat, row_j, index, temp);
  }
}

/* Equivalent to Row Operation Ri <- Ri + cRj*/
void row_addition(Matrix mat, size_t row_i, size_t row_j, double c) {
  double sum;
  for (size_t col = 0; col < mat.dim; col++) {
    sum = mat_get(mat, row_i, col);
    sum += c * mat_get(mat, row_j, col);
    mat_set(mat, row_i, col, sum);
  }
}


