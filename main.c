#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <stdbool.h>
#include <math.h>

#define EPS_ABS 1e-12
#define EPS_REL 1e-9

int dbl_eq(double a, double b) {
    double diff = fabs(a - b);
    
    if (diff < EPS_ABS) return 1;
    return diff <= EPS_REL * fmax(fabs(a), fabs(b));
}

typedef struct {
    size_t dimension;
    double * mat;
} Matrix;

double mat_get(Matrix matrix, size_t dimension, size_t row, size_t col) {
    return matrix.mat[dimension * row + col];
}

void mat_set(Matrix matrix, size_t dimension, size_t row, size_t col, double val) {
    matrix.mat[dimension * row + col] = val;
}

void debug_matrix(Matrix mat, size_t dimension) {
    for (size_t row = 0; row < dimension; row++) {
        printf("\n");
        for (size_t col = 0; col < dimension; col++) {
            printf("%lf ", mat_get(mat, dimension, row, col));
        }
    }
    printf("\n");
}

void swap_rows(Matrix mat, size_t dimension, size_t row_i, size_t row_j) {
    for (size_t index = 0; index < dimension; index++) {
        double temp = mat_get(mat, dimension, row_i, index);
        mat_set(mat, dimension, row_i, index, matget(mat, dimension, row_j, index));
        mat_set(mat, dimension, row_j, index, temp);
    }
}

void prep_REF(Matrix mat, size_t dimension) {
    for (size_t col = 0; col < dimension; col++) {
        for (size_t row = 0; row < dimension; row++) {
            if (!dbl_eq(mat_get(mat, dimension, col, row), 0)) {
                swap_rows(mat, dimension, row, 0);
                return;
            }
        }
    }
}

void REF(double * matrix, size_t dimension) {
    
    prep_REF(matrix, dimension);
    
    double c;
    for (size_t j = 0; j < dimension; j++) {
        for (size_t i = 0; i < dimension; i++) {
            if (i>j) {
                c=matrix[i * dimension + j] / matrix[j * dimension + j];
                
                for (size_t k = 0; k < dimension; k++) 
                    matrix[i * dimension + k] -= c * matrix[j * dimension + k];
            }
        }
    }
}

int main() {
    size_t dimension;
    printf("SQUARE MATRIX IN THE COMPLEX PLANE\n");
    printf("Enter Dimension: ");
    scanf("%zu", &dimension);
    
    double mat[dimension][dimension];
    
    printf("Enter Matrix:\n");
    for (size_t col = 0; col < dimension; col++) {
        for (size_t element_index = 0; element_index < dimension; element_index++) {
            scanf("%lf", &mat[col][element_index]);
        }
    }
    

    Matrix matrix;
    matrix.mat = mat[0];
    printf("Top Right: %lf\n", mat_get(matrix, dimension, 0, 0));
    return 0;

    // REF(matrix[0], dimension);
    REF(mat[0], dimension);
    debug_matrix(mat[0], dimension);
}