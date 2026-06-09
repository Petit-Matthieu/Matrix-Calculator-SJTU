#include "matrix_core.h"
#include "matrix_lu.h"
#include <stdio.h>

int main(void)
{
    Matrix A, L, U, b, x;
    MatrixInit(&A); MatrixInit(&L); MatrixInit(&U); MatrixInit(&b); MatrixInit(&x);
    MatrixCreate(&A, 2, 2);
    MatrixCreate(&L, 2, 2);
    MatrixCreate(&U, 2, 2);
    MatrixCreate(&b, 2, 1);
    MatrixCreate(&x, 2, 1);

    MatrixSet(&A, 0, 0, 4.0);
    MatrixSet(&A, 0, 1, 3.0);
    MatrixSet(&A, 1, 0, 6.0);
    MatrixSet(&A, 1, 1, 3.0);
    MatrixSet(&b, 0, 0, 10.0);
    MatrixSet(&b, 1, 0, 12.0);

    MatrixError error = LUDecomposeNoPivot(&A, &L, &U, 1e-12);
    printf("LU status: %s\n", MatrixErrorMessage(error));
    MatrixPrint(&L, "L");
    MatrixPrint(&U, "U");

    error = LUSolve(&L, &U, &b, &x, 1e-12);
    printf("LUSolve status: %s\n", MatrixErrorMessage(error));
    MatrixPrint(&x, "x");

    MatrixFree(&A); MatrixFree(&L); MatrixFree(&U); MatrixFree(&b); MatrixFree(&x);
    return 0;
}
