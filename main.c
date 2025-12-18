#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <stdbool.h>
#include <math.h>

#define EPS_ABS 1e-12
#define EPS_REL 1e-9

int dbl_eq(double a, double b) {
    double diff = fabs(a - b);
    
    // if (diff < EPS_ABS) return 1;
    return (diff < EPS_ABS) || (diff <= EPS_REL * fmax(fabs(a), fabs(b)));
}

typedef struct {
    size_t dimension;
    double * mat;
} Matrix;

double mat_get(Matrix mat, size_t row, size_t col) {
    return mat.mat[mat.dimension * row + col];
}

void mat_set(Matrix mat, size_t row, size_t col, double val) {
    mat.mat[mat.dimension * row + col] = val;
}

void debug_matrix(Matrix mat) {
    for (size_t row = 0; row < mat.dimension; row++) {
        if (row) printf("\n");
        for (size_t col = 0; col < mat.dimension; col++) {
            printf("%.2lf ", mat_get(mat, row, col));
        }
    }
    printf("\n");
}

void swap_rows(Matrix mat, size_t row_i, size_t row_j) {
    if (row_i == row_j) return;
    
    for (size_t index = 0; index < mat.dimension; index++) {
        double temp = mat_get(mat, row_i, index);
        mat_set(mat, row_i, index, mat_get(mat, row_j, index));
        mat_set(mat, row_j, index, temp);
    }
}

void row_addition(Matrix mat, size_t row_i, size_t row_j, double c) {
    for (size_t col = 0; col < mat.dimension; col++) {
        double sum;
        sum = mat_get(mat, row_i, col);
        sum += c * mat_get(mat, row_j, col);
        mat_set(mat, row_i, col, sum);
    }
}

/* searches and swaps given top row with a lower row that has a more further pivot*/
/* returns true if it has more pivots, false if it's a zero matrix starting from init_row*/
bool prep_REF(Matrix mat, size_t init_row) {
    for (size_t col = 0; col < mat.dimension; col++){
        for (size_t row = init_row; row < mat.dimension; row++) {
            if (dbl_eq(mat_get(mat, row, col), 0)) continue; 
            
            swap_rows(mat, row, init_row);
            return true;
        }
    }
    return false;
}

void REF(Matrix mat) {  
    double c;
    
    for (size_t outer = 0; outer < mat.dimension; outer++) {

        printf("PRE PREP ITER %zu\n", outer);
        debug_matrix(mat);

        if (!prep_REF(mat, outer)) return;

        printf("POST PREP ITER %zu\n", outer);
        debug_matrix(mat);

        for (size_t inner = 0; inner < mat.dimension; inner++) {
            if (inner<=outer) 
                continue;
                
            c = mat_get(mat, inner, outer) 
                / mat_get(mat, outer, outer); 
            
            row_addition(mat, inner, outer, -1 * c);
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
    mat.dimension = dimension;
    
    REF(mat);
    
    printf("REF FORM\n");
    
    debug_matrix(mat);
}