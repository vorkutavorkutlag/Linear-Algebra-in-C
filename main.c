#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <stdbool.h>
#include <math.h>

#define EPS_ABS 1e-12
#define EPS_REL 1e-9

void debug_matrix(double * matrix, size_t dimension) {
    for (size_t col = 0; col < dimension; col++) {
        printf("\n");
        for (size_t element_index = 0; element_index < dimension; element_index++) {
            printf("%lf ", matrix[col * dimension + element_index]);
        }
    }

    printf("\n");
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

    debug_matrix(matrix[0], dimension);
}