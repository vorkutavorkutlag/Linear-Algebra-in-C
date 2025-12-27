#include <stdlib.h>

typedef struct {
  double real;
  double im;
} Complex;

typedef struct {
	size_t degree;
  Complex * coefficients; // of size `degree`. array of the complex coefficients of every degree
} Polynomial;


Complex complex_sum(Complex comp1, Complex comp2);
Complex complex_product(Complex comp1, Complex comp2);

Polynomial polynomial_product(Polynomial poly1, Polynomial poly2);
void polynomial_scalar_mult(Polynomial * poly, Complex scalar);
void polynomial_free(Polynomial poly);

