#include <stdlib.h>

typedef struct {
  double real;
  double im;
} Complex;

typedef struct {
	size_t degree;
  Complex * coefficients; // of size `degree`. array of the complex coefficients of every degree
} Polynomial;

int dbl_eq(double d1, double d2);

Complex complex_sum(Complex comp1, Complex comp2);
Complex complex_product(Complex comp1, Complex comp2);
Complex complex_division(Complex comp1, Complex comp2);
Complex complex_neg(Complex comp);
Complex complex_conjugate(Complex comp);
int complex_nan(Complex comp);
double complex_abs(Complex comp);
double complex_diff(Complex comp1, Complex comp2);

Polynomial polynomial_product(Polynomial poly1, Polynomial poly2);
Polynomial polynomial_division(Polynomial poly1, Polynomial poly2);
Polynomial polynomial_sum(Polynomial poly1, Polynomial poly2);
Polynomial polynomial_derivative(Polynomial poly);

int dummy_assignment(Complex * c1, Complex c2);
Complex cauchy_nr_root(Polynomial poly);
Complex * polynomial_all_roots(Polynomial poly);



void polynomial_scalar_mult(Polynomial * poly, Complex scalar);
void polynomial_free(Polynomial poly);

