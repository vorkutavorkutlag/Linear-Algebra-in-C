#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>

int main() {
    size_t dimension;
    printf("SQUARE MATRIX IN THE COMPLEX PLANE\n");
    printf("Enter Dimension: ");
    scanf("%zu", &dimension);
    
    double matrix[dimension][dimension];
    
    printf("Enter Matrix by Columns:\n");
    for (size_t col = 0; col < dimension; col++) {
        for (size_t element_index = 0; element_index < dimension; element_index++) {
            scanf("%lf", &matrix[col][element_index]);
        }
    }

    printf("Attemping to print.");

    for (size_t col = 0; col < dimension; col++) {
        printf("\n");
        for (size_t element_index = 0; element_index < dimension; element_index++) {
            printf("%lf ", matrix[col][element_index]);
        }
    }
    printf("\n");
}