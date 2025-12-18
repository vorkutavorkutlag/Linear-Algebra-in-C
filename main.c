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
/* returns further most pivot's index if it has more pivots*/
/*returns -1 if it's a zero matrix starting from init_row*/
__ssize_t __prep_REF(Matrix mat, size_t init_row) {
    for (size_t col = 0; col < mat.dimension; col++){
        for (size_t row = init_row; row < mat.dimension; row++) {
            if (dbl_eq(mat_get(mat, row, col), 0)) continue; 
            swap_rows(mat, row, init_row);
            return col;
        }
    }
    return -1;
}

void REF(Matrix mat) {  
    double scalar;
    __ssize_t p_col;
    
    for (size_t bot_row = 0; bot_row < mat.dimension; bot_row++) {
        
        if ((p_col = __prep_REF(mat, bot_row)) == -1) return;

        for (size_t top_row = bot_row+1; top_row < mat.dimension; top_row++) {

            scalar = mat_get(mat, top_row, p_col) 
                / mat_get(mat, bot_row, p_col); 
            
            row_addition(mat, top_row, bot_row, -1 * scalar);
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