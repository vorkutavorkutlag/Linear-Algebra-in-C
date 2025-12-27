#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <stdbool.h>
#include <math.h>

#define EPS_ABS 1e-12
#define EPS_REL 1e-9

int dbl_eq(double a, double b) {
  double diff = fabs(a - b);
  return (diff < EPS_ABS) || (diff <= EPS_REL * fmax(fabs(a), fabs(b)));
}

typedef struct {
  size_t dim;
  double * mat;
} Matrix;

typedef struct {
  double real;
  double im;
} Complex;

typedef struct {
	size_t degree;
  Complex * coefficients; // of size `degree`. array of the complex coefficients of every degree
} Polynomial;

Complex complex_sum(Complex comp1, Complex comp2) {
  return (Complex) {comp1.real + comp2.real, comp1.im + comp2.im};
}

Complex complex_product(Complex comp1, Complex comp2) {
  Complex c1, c2;
  c1 = (Complex) {comp1.real * comp2.real, comp1.real * comp2.im};
  c2 = (Complex) {-1 * comp1.im * comp2.im, comp1.im * comp1.real};
  return complex_sum(c1, c2);
}

 Polynomial polynomial_product(Polynomial poly1, Polynomial poly2) {
  size_t prod_degree = poly1.degree + poly2.degree;
  Complex * prod_coefficients = (Complex *) calloc(prod_degree, sizeof(Complex));

  for (size_t degree1 = 0; degree1 <= poly1.degree; degree1++) {
    for (size_t degree2 = 0; degree2 <= poly2.degree; degree2++) {
      prod_coefficients[degree1+degree2] = complex_sum(prod_coefficients[degree1+degree2], 
                                              complex_product(poly1.coefficients[degree1], 
                                                              poly2.coefficients[degree2]));
    }
  }

  return (Polynomial) {prod_degree, prod_coefficients};
}

void polynomial_free(Polynomial poly) {
  free(poly.coefficients);
}

void polynomial_scalar_mult(Polynomial * poly, Complex scalar) {
  for (size_t deg = 0; deg <= poly->degree; deg++)
    poly->coefficients[deg] = complex_product(poly->coefficients[deg], scalar);
}

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


/* given an upper triangular matrix, returns its determinant */
/* that being, the product of its diagonal entries */
double __det(Matrix mat) {
  double det = 1;
  for (size_t row_col = 0; row_col < mat.dim; row_col++) 
    det *= mat_get(mat, row_col, row_col);
  return det;
}

/* searches and swaps given top row with a lower row that has a more further pivot */
/* returns further most pivot's index if it has more pivots */
/* returns -1 if it's a zero matrix starting from init_row */
__ssize_t __prep_REF(Matrix mat, size_t init_row, bool * swap) {
  for (size_t col = 0; col < mat.dim; col++){
    for (size_t row = init_row; row < mat.dim; row++) {
      if (dbl_eq(mat_get(mat, row, col), 0)) continue;
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
      if ((p_col = __prep_REF(mat, top_row, &swap)) == -1) return 0.0;
      det_coefficient *= swap? -1: 1; // multiply by -1 if swapped
      det *= mat_get(mat, top_row, top_row);

      printf("DEBUG\n");
      debug_matrix(mat);
      printf("\n");

      for (size_t bot_row = top_row+1; bot_row < mat.dim; bot_row++) {
        scalar = mat_get(mat, bot_row, p_col) 
                / mat_get(mat, top_row, p_col); 
            
        row_addition(mat, bot_row, top_row, -1.0 * scalar);
      }
  }

  return det_coefficient * det;  
}

int main() {

  Complex poly_c1[] = { (Complex) {1, 0}, (Complex) {3, 0}, (Complex) {4, 0} };
  Complex poly_c2[] = { (Complex) {0,0}, (Complex) {1, 0}, (Complex) {2, 0} };

  Polynomial poly1 = (Polynomial) {2, poly_c1};
  Polynomial poly2 = (Polynomial) {2, poly_c2};
  Polynomial prod = polynomial_product(poly1, poly2);

  for (size_t deg = 0; deg <= prod.degree; deg++)
    printf("%.2lf\n", prod.coefficients[prod.degree-deg].real);
  printf("\n");
  return 0;

  size_t dimension;
  printf("SQUARE MATRIX IN THE COMPLEX PLANE\n");
  printf("Enter Dimension: ");
  scanf("%zu", &dimension);
   
  double matrix[dimension][dimension];
    
  printf("Enter Matrix:\n");
  for (size_t col = 0; col < dimension; col++) {
    for (size_t element_index = 0; element_index < dimension; element_index++) {
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
}
