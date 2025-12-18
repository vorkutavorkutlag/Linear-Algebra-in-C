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
        if (row) printf("\n");
        for (size_t col = 0; col < dimension; col++) {
            printf("%.2lf ", mat_get(mat, dimension, row, col));
        }
    }
    printf("\n");
}

void swap_rows(Matrix mat, size_t dimension, size_t row_i, size_t row_j) {
    if (row_i == row_j) return;
    
    for (size_t index = 0; index < dimension; index++) {
        double temp = mat_get(mat, dimension, row_i, index);
        mat_set(mat, dimension, row_i, index, mat_get(mat, dimension, row_j, index));
        mat_set(mat, dimension, row_j, index, temp);
    }
}

void row_addition(Matrix mat, size_t dimension, size_t row_i, size_t row_j, double c) {
    for (size_t index = 0; index < dimension; index++) {
        double new;
        new = mat_get(mat, dimension, row_i, index);
        new += c * mat_get(mat, dimension, row_j, index);
        mat_set(mat, dimension, row_i, index, new);
    }
}

/* searches and swaps top row with a row that has a more further pivot*/
/* returns true if it has more pivots, false if it's a zero matrix starting from init_row*/
bool prep_REF(Matrix mat, size_t dimension, size_t init_row) {
    printf("PREP\n");

    for (size_t col = 0; col < dimension; col++){
        for (size_t row = init_row; row < dimension; row++) {
            if (!dbl_eq(mat_get(mat, dimension, row, col), 0)) {
                printf("found pivot %lf at %zu, %zu\n", mat_get(mat, dimension, row, col), row+1, col+1);
                swap_rows(mat, dimension, row, init_row);
                return true;
            }
        }
    }

    return false;
}

void REF(Matrix mat, size_t dimension) {  
    double c;
    for (size_t outer = 0; outer < dimension; outer++) {
        
        // printf("ITER PRE PREP %zu\n", outer);
        // debug_matrix(mat, dimension);
        
        if (!prep_REF(mat, dimension, outer)) return;

        // printf("ITER POST PREP %zu\n", outer);
        // debug_matrix(mat, dimension);

        for (size_t inner = 0; inner < dimension; inner++) {
            if (inner<=outer) 
                continue;
                
            c = mat_get(mat, dimension, inner, outer) / mat_get(mat, dimension, outer, outer); 
            
            row_addition(mat, dimension, inner, outer, -1 * c);

        }
    }
}

int main() {
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

    REF(mat, dimension);
    
    printf("REF FORM\n");
    
    debug_matrix(mat, dimension);
}