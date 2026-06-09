#include "matrix_core.h"
#include "matrix_lu.h"
#include <stdio.h>

int main(void)
{
    Matrix A, L, U, b, x;
    MatrixInit(&A); MatrixInit(&L); MatrixInit(&U); MatrixInit(&b); MatrixInit(&x);
    MatrixCreate(&A, 3, 3);
    MatrixCreate(&L, 3, 3);
    MatrixCreate(&U, 3, 3);
    MatrixCreate(&b, 3, 1);
    MatrixCreate(&x, 3, 1);

    REAL a[] = {
        2, 1, 1,
        4, -6, 0,
        -2, 7, 2
    };
    REAL rhs[] = {5, -2, 9};
    for (int i = 0; i < 9; ++i) A.data[i] = a[i];
    for (int i = 0; i < 3; ++i) b.data[i] = rhs[i];

    MatrixError error = LUDecomposeNoPivot(&A, &L, &U, 1e-12);
    printf("LU status: %s\n", MatrixErrorMessage(error));
    MatrixPrint(&A, "A");
    MatrixPrint(&L, "L");
    MatrixPrint(&U, "U");

    error = LUSolve(&L, &U, &b, &x, 1e-12);
    printf("LUSolve status: %s\n", MatrixErrorMessage(error));
    MatrixPrint(&x, "x");

    REAL det = 0.0;
    LUDeterminant(&U, &det);
    printf("det(A) = %.6f\n", det);

    MatrixFree(&A); MatrixFree(&L); MatrixFree(&U); MatrixFree(&b); MatrixFree(&x);
    return 0;
}
