#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "matrix.h"

/* given an upper triangular matrix, returns its determinant */
/* that being, the product of its diagonal entries */
double __det(Matrix mat) {
  double det = 1;
  for (size_t row_col = 0; row_col < mat.dim; row_col++)
    det *= mat_get(mat, row_col, row_col);
  return det;
}

/* searches and swaps given top row with a lower row that has a more further
 * pivot */
/* returns further most pivot's index if it has more pivots */
/* returns -1 if it's a zero matrix starting from init_row */
__ssize_t __prep_REF(Matrix mat, size_t init_row, bool *swap) {
  for (size_t col = 0; col < mat.dim; col++) {
    for (size_t row = init_row; row < mat.dim; row++) {
      if (dbl_eq(mat_get(mat, row, col), 0))
        continue;
      if (*swap = row != init_row)
        swap_rows(mat, row, init_row);
      return col;
    }
  }
  return -1;
}

/* given matrix, rewrites it in REF form and returns its determinant */
double GEM(Matrix mat) {
  int det_coefficient = 1;
  double det = 1;
  __ssize_t p_col;
  double scalar;
  bool swap;

  for (size_t top_row = 0; top_row < mat.dim; top_row++) {
    if ((p_col = __prep_REF(mat, top_row, &swap)) == -1)
      return 0.0;
    det_coefficient *= swap ? -1 : 1; // multiply by -1 if swapped
    det *= mat_get(mat, top_row, top_row);

    printf("DEBUG\n");
    debug_matrix(mat);
    printf("\n");

    for (size_t bot_row = top_row + 1; bot_row < mat.dim; bot_row++) {
      scalar = mat_get(mat, bot_row, p_col) / mat_get(mat, top_row, p_col);

      row_addition(mat, bot_row, top_row, -1.0 * scalar);
    }
  }

  return det_coefficient * det;
}

int main() {
  char route;
  printf(
      "Square matrix determinant, or solve characteristic polynomial? [M/P]\n");
  if (scanf("%c", &route) != 1) {
    fprintf(stderr, "Invalid input. Choose M or P.");
    return EXIT_FAILURE;
  }

  if (route == 'M') {
    size_t dimension;
    printf("Enter dimension of real square matrix: ");
    scanf("%zu", &dimension);

    double matrix[dimension][dimension];

    printf("Enter Matrix (same as visually):\n");
    for (size_t col = 0; col < dimension; col++) {
      for (size_t element_index = 0; element_index < dimension;
           element_index++) {
        scanf("%lf", &matrix[col][element_index]);
      }
    }

    Matrix mat;
    mat.mat = matrix[0];
    mat.dim = dimension;

    double det = GEM(mat);

    printf("REF FORM\n");

    debug_matrix(mat);

    printf("Determinant: %g\n", det);
    return EXIT_SUCCESS;
  } else if (route == 'P') {
    size_t degree;

    printf("Enter polynomial degree:\n");
    scanf("%zu", &degree);

    Complex coefficients[degree + 1];
    printf("Enter real polynomial coefficients in rising degree order:\n");

    for (size_t index = 0; index <= degree; ++index) {
      if (scanf("%lf", coefficients + index) != 1) {
        fprintf(stderr, "Invalid input. Example: '2 0 1' == 'x^2 + 2'\n");
        return EXIT_FAILURE;
      }
    }

    Polynomial poly = {degree, coefficients};
    Complex *roots = polynomial_all_roots(poly);

    for (size_t index = 0; index < poly.degree; ++index) {
      printf("Root %zu: %lf + %lfi\n", index, roots[index].real,
             roots[index].im);
    }

    free(roots);

    return EXIT_SUCCESS;
  }

  fprintf(stderr, "Unkown option.\n");
  return EXIT_FAILURE;
}
