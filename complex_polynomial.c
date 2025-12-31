#include <stddef.h>
#include <stdio.h>
#include "complex_polynomial.h"

#define C_ZERO (Complex) {0, 0}
#define C_ONE (Complex) {1, 0}

Complex complex_sum(Complex comp1, Complex comp2) {
  return (Complex) {comp1.real + comp2.real, comp1.im + comp2.im};
}

Complex complex_conjugate(Complex comp) {
  return (Complex) {comp.real, -1.0 * comp.im};
}

Complex complex_neg(Complex comp) {
  return (Complex) {-1.0 * comp.real, -1.0 * comp.im};
}

Complex complex_product(Complex comp1, Complex comp2) {
  Complex c1, c2;
  c1 = (Complex) {comp1.real * comp2.real, comp1.real * comp2.im};
  c2 = (Complex) {-1 * comp1.im * comp2.im, comp1.im * comp1.real};
  return complex_sum(c1, c2);
}

/* Divides comp1 by comp2*/
/* Multiplied by conjugate, then divded. */
Complex complex_division(Complex comp1, Complex comp2) {
  Complex nom = complex_product(comp1, complex_conjugate(comp2));
  double denom = comp2.real * comp2.real + comp2.im + comp2.im;
  return (Complex) {nom.real / denom, nom.im / denom};
}

Complex polynomial_evaluate(Polynomial poly, Complex val) {
  Complex sum = C_ZERO;
  for (size_t cur_degree = 0; cur_degree <= poly.degree; cur_degree++) {
    Complex mult = C_ONE;
    
    // for every degree, mult = 1, mult = x, mult = x^2
    for (size_t mult_degree = 0; mult_degree < cur_degree; mult_degree++)
      mult = complex_product(mult, val);
    
    mult = complex_product(mult, poly.coefficients[cur_degree]);
    sum = complex_sum(sum, mult);
  }

  return sum;
}

/* Creates new polynomial on the heap. */
Polynomial polynomial_sum(Polynomial poly1, Polynomial poly2) {
  size_t sum_degree = poly1.degree >= poly2.degree ? poly1.degree : poly2.degree;
  Complex * sum_coefficients = (Complex *) calloc(sum_degree+1, sizeof(Complex));
  for (size_t index = 0; index <= sum_degree; index++) {
    sum_coefficients[index] = complex_sum(sum_coefficients[index],
                                          poly1.degree >= sum_degree ? poly1.coefficients[index] : C_ZERO);
    sum_coefficients[index] = complex_sum(sum_coefficients[index],
                                          poly2.degree >= sum_degree ? poly2.coefficients[index] : C_ZERO);
  }

  return (Polynomial) {sum_degree, sum_coefficients};
}

/* Creates a new polynomial on the heap. */
Polynomial polynomial_derivative(Polynomial poly) {
  size_t der_degree = poly.degree - 1;
  Complex * der_coefficients = (Complex *) calloc(der_degree, sizeof(Complex));

  for (size_t index = 1; index <= poly.degree; index++) {
    der_coefficients[index-1] = complex_product(poly.coefficients[index], (Complex) { (double) index, 0});
  }
  
  return (Polynomial) {der_degree, der_coefficients};
}

/* Creates new polynomial on the heap. */
Polynomial polynomial_product(Polynomial poly1, Polynomial poly2) {
  size_t prod_degree = poly1.degree + poly2.degree;
  Complex * prod_coefficients = (Complex *) calloc(prod_degree+1, sizeof(Complex));

  for (size_t degree1 = 0; degree1 <= poly1.degree; degree1++) {
    for (size_t degree2 = 0; degree2 <= poly2.degree; degree2++) {
      prod_coefficients[degree1+degree2] = complex_sum(prod_coefficients[degree1+degree2], 
                                              complex_product(poly1.coefficients[degree1], 
                                                              poly2.coefficients[degree2]));
    }
  }

  return (Polynomial) {prod_degree, prod_coefficients};
}

/* Assumes poly2 is single-degree. Uses synthetic division
 * Creates new polynomial on the heap. */
Polynomial polynomial_division(Polynomial poly1, Polynomial poly2) {
  // zx+y = 0 => x=-y/z
  Complex root = complex_neg(complex_division(poly2.coefficients[0], poly2.coefficients[1]));
  
  Polynomial quotient;
  quotient.degree = poly2.degree - 1;
  quotient.coefficients = (Complex *) calloc (poly2.degree, sizeof(Complex)); 

  Complex cur_coefficient = poly2.coefficients[poly2.degree];

  for (size_t cur_deg = poly1.degree; cur_deg >= 1; cur_deg--) {
    quotient.coefficients[cur_deg - 1] = cur_coefficient;
    cur_coefficient = complex_sum(complex_product(root, cur_coefficient), poly1.coefficients[cur_deg - 1]);
  }

  return quotient;
}

 

void polynomial_free(Polynomial poly) {
  free(poly.coefficients);
}

/* Modifies polynomial directly. */
void polynomial_scalar_mult(Polynomial * poly, Complex scalar) {
  for (size_t deg = 0; deg <= poly->degree; deg++)
    poly->coefficients[deg] = complex_product(poly->coefficients[deg], scalar);
}


