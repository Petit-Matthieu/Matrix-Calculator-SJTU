#ifndef MATRIX_OPS_H
#define MATRIX_OPS_H

#include "matrix_core.h"

/*
 * Matrix operation module
 *
 * The first four functions are the required HW3 operations. Additional
 * multiplication versions are retained for the performance comparison.
 */

/* Elementwise operations and transpose. */
MatrixError MatrixAdd(const Matrix *A, const Matrix *B, Matrix *C);
MatrixError MatrixSub(const Matrix *A, const Matrix *B, Matrix *C);
MatrixError MatrixScale(REAL alpha, const Matrix *A, Matrix *B);
MatrixError MatrixTranspose(const Matrix *A, Matrix *AT);
MatrixError MatrixNormFrobenius(const Matrix *A, REAL *norm_value);

/* Production multiplication and retained performance-comparison versions. */
MatrixError MatrixMultiply(const Matrix *A, const Matrix *B, Matrix *C);
MatrixError MatrixMultiplyOriginal(const Matrix *A, const Matrix *B, Matrix *C);
MatrixError MatrixMultiplyIKJ(const Matrix *A, const Matrix *B, Matrix *C);
MatrixError MatrixMultiplyBlocked(const Matrix *A, const Matrix *B, Matrix *C);
MatrixError MatrixMultiplyOptimized(const Matrix *A, const Matrix *B, Matrix *C);

#endif
