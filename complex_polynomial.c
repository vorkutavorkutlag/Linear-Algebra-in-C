#include "complex_polynomial.h"

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


