#include "matrix_ops.h"

#include <math.h>

#define MATRIX_BLOCK_SIZE 64
#define MATRIX_PORTABLE_BLOCK_SIZE 64

/*
 * Addition and subtraction share the same dimension requirement:
 * A, B and the output C must all have identical shapes.
 */
static MatrixError CheckBinaryElementwise(const Matrix *A, const Matrix *B, Matrix *C)
{
    if (!MatrixIsValid(A) || !MatrixIsValid(B) || !MatrixIsValid(C)) {
        return MATRIX_ERROR_NULL_POINTER;
    }
    if (A->row != B->row || A->column != B->column) {
        return MATRIX_ERROR_SIZE_MISMATCH;
    }
    if (C->row != A->row || C->column != A->column) {
        return MATRIX_ERROR_SIZE_MISMATCH;
    }
    return MATRIX_SUCCESS;
}

static MatrixError CheckMultiply(const Matrix *A, const Matrix *B, Matrix *C)
{
    if (!MatrixIsValid(A) || !MatrixIsValid(B) || !MatrixIsValid(C)) {
        return MATRIX_ERROR_NULL_POINTER;
    }
    if (A->column != B->row || C->row != A->row || C->column != B->column) {
        return MATRIX_ERROR_SIZE_MISMATCH;
    }
    return MATRIX_SUCCESS;
}

static int MinInt(int a, int b)
{
    return a < b ? a : b;
}

MatrixError MatrixAdd(const Matrix *A, const Matrix *B, Matrix *C)
{
    MatrixError error = CheckBinaryElementwise(A, B, C);
    if (error != MATRIX_SUCCESS) {
        return error;
    }
    int total = A->row * A->column;
    for (int k = 0; k < total; ++k) {
        C->data[k] = A->data[k] + B->data[k];
    }
    return MATRIX_SUCCESS;
}

MatrixError MatrixSub(const Matrix *A, const Matrix *B, Matrix *C)
{
    MatrixError error = CheckBinaryElementwise(A, B, C);
    if (error != MATRIX_SUCCESS) {
        return error;
    }
    int total = A->row * A->column;
    for (int k = 0; k < total; ++k) {
        C->data[k] = A->data[k] - B->data[k];
    }
    return MATRIX_SUCCESS;
}

MatrixError MatrixScale(REAL alpha, const Matrix *A, Matrix *B)
{
    if (!MatrixIsValid(A) || !MatrixIsValid(B)) {
        return MATRIX_ERROR_NULL_POINTER;
    }
    if (A->row != B->row || A->column != B->column) {
        return MATRIX_ERROR_SIZE_MISMATCH;
    }
    int total = A->row * A->column;
    for (int k = 0; k < total; ++k) {
        B->data[k] = alpha * A->data[k];
    }
    return MATRIX_SUCCESS;
}

MatrixError MatrixTranspose(const Matrix *A, Matrix *AT)
{
    if (!MatrixIsValid(A) || !MatrixIsValid(AT)) {
        return MATRIX_ERROR_NULL_POINTER;
    }
    if (AT->row != A->column || AT->column != A->row) {
        return MATRIX_ERROR_SIZE_MISMATCH;
    }
    for (int i = 0; i < A->row; ++i) {
        for (int j = 0; j < A->column; ++j) {
            AT->data[MatrixIndex(AT, j, i)] = A->data[MatrixIndex(A, i, j)];
        }
    }
    return MATRIX_SUCCESS;
}

MatrixError MatrixNormFrobenius(const Matrix *A, REAL *norm_value)
{
    if (!MatrixIsValid(A) || norm_value == NULL) {
        return MATRIX_ERROR_NULL_POINTER;
    }
    REAL sum = 0.0;
    int total = A->row * A->column;
    for (int k = 0; k < total; ++k) {
        sum += A->data[k] * A->data[k];
    }
    *norm_value = sqrt(sum);
    return MATRIX_SUCCESS;
}

MatrixError MatrixMultiplyOriginal(const Matrix *A, const Matrix *B, Matrix *C)
{
    MatrixError error = CheckMultiply(A, B, C);
    if (error != MATRIX_SUCCESS) {
        return error;
    }
    MatrixFillZero(C);
    for (int i = 0; i < A->row; ++i) {
        for (int j = 0; j < B->column; ++j) {
            REAL sum = 0.0;
            for (int k = 0; k < A->column; ++k) {
                sum += A->data[MatrixIndex(A, i, k)] * B->data[MatrixIndex(B, k, j)];
            }
            C->data[MatrixIndex(C, i, j)] = sum;
        }
    }
    return MATRIX_SUCCESS;
}

MatrixError MatrixMultiplyIKJ(const Matrix *A, const Matrix *B, Matrix *C)
{
    MatrixError error = CheckMultiply(A, B, C);
    if (error != MATRIX_SUCCESS) {
        return error;
    }
    MatrixFillZero(C);
    for (int i = 0; i < A->row; ++i) {
        for (int k = 0; k < A->column; ++k) {
            REAL aik = A->data[MatrixIndex(A, i, k)];
            for (int j = 0; j < B->column; ++j) {
                C->data[MatrixIndex(C, i, j)] += aik * B->data[MatrixIndex(B, k, j)];
            }
        }
    }
    return MATRIX_SUCCESS;
}

MatrixError MatrixMultiplyBlocked(const Matrix *A, const Matrix *B, Matrix *C)
{
    MatrixError error = CheckMultiply(A, B, C);
    if (error != MATRIX_SUCCESS) {
        return error;
    }
    MatrixFillZero(C);

    for (int ii = 0; ii < A->row; ii += MATRIX_BLOCK_SIZE) {
        int i_end = MinInt(ii + MATRIX_BLOCK_SIZE, A->row);
        for (int kk = 0; kk < A->column; kk += MATRIX_BLOCK_SIZE) {
            int k_end = MinInt(kk + MATRIX_BLOCK_SIZE, A->column);
            for (int jj = 0; jj < B->column; jj += MATRIX_BLOCK_SIZE) {
                int j_end = MinInt(jj + MATRIX_BLOCK_SIZE, B->column);

                for (int i = ii; i < i_end; ++i) {
                    const REAL *a_row = A->data + i * A->column;
                    REAL *c_row = C->data + i * C->column;
                    for (int k = kk; k < k_end; ++k) {
                        REAL aik = a_row[k];
                        const REAL *b_row = B->data + k * B->column;
                        int j = jj;

                        for (; j + 3 < j_end; j += 4) {
                            c_row[j] += aik * b_row[j];
                            c_row[j + 1] += aik * b_row[j + 1];
                            c_row[j + 2] += aik * b_row[j + 2];
                            c_row[j + 3] += aik * b_row[j + 3];
                        }
                        for (; j < j_end; ++j) {
                            c_row[j] += aik * b_row[j];
                        }
                    }
                }
            }
        }
    }
    return MATRIX_SUCCESS;
}

static void Multiply4x4(const REAL * restrict a, int a_stride,
                        const REAL * restrict b, int b_stride,
                        REAL * restrict c, int c_stride, int inner)
{
    REAL c00 = 0.0, c01 = 0.0, c02 = 0.0, c03 = 0.0;
    REAL c10 = 0.0, c11 = 0.0, c12 = 0.0, c13 = 0.0;
    REAL c20 = 0.0, c21 = 0.0, c22 = 0.0, c23 = 0.0;
    REAL c30 = 0.0, c31 = 0.0, c32 = 0.0, c33 = 0.0;

    const REAL *a0 = a;
    const REAL *a1 = a + a_stride;
    const REAL *a2 = a + 2 * a_stride;
    const REAL *a3 = a + 3 * a_stride;
    for (int k = 0; k < inner; ++k) {
        const REAL *b_row = b + k * b_stride;
        REAL b0 = b_row[0];
        REAL b1 = b_row[1];
        REAL b2 = b_row[2];
        REAL b3 = b_row[3];
        REAL a0k = a0[k];
        REAL a1k = a1[k];
        REAL a2k = a2[k];
        REAL a3k = a3[k];

        c00 += a0k * b0; c01 += a0k * b1; c02 += a0k * b2; c03 += a0k * b3;
        c10 += a1k * b0; c11 += a1k * b1; c12 += a1k * b2; c13 += a1k * b3;
        c20 += a2k * b0; c21 += a2k * b1; c22 += a2k * b2; c23 += a2k * b3;
        c30 += a3k * b0; c31 += a3k * b1; c32 += a3k * b2; c33 += a3k * b3;
    }

    c[0] = c00; c[1] = c01; c[2] = c02; c[3] = c03;
    c += c_stride;
    c[0] = c10; c[1] = c11; c[2] = c12; c[3] = c13;
    c += c_stride;
    c[0] = c20; c[1] = c21; c[2] = c22; c[3] = c23;
    c += c_stride;
    c[0] = c30; c[1] = c31; c[2] = c32; c[3] = c33;
}

static void MultiplyTail(const Matrix *A, const Matrix *B, Matrix *C,
                         int i_begin, int i_end, int j_begin, int j_end)
{
    for (int i = i_begin; i < i_end; ++i) {
        const REAL *a_row = A->data + i * A->column;
        REAL *c_row = C->data + i * C->column;
        for (int j = j_begin; j < j_end; ++j) {
            REAL sum = 0.0;
            for (int k = 0; k < A->column; ++k) {
                sum += a_row[k] * B->data[k * B->column + j];
            }
            c_row[j] = sum;
        }
    }
}

MatrixError MatrixMultiplyOptimized(const Matrix *A, const Matrix *B, Matrix *C)
{
    MatrixError error = CheckMultiply(A, B, C);
    if (error != MATRIX_SUCCESS) {
        return error;
    }

    /*
     * Standard-C 4 x 4 microkernel. The column block keeps the active slice
     * of B cache-friendly while local accumulators reduce repeated C traffic.
     */
    for (int jj = 0; jj < B->column; jj += MATRIX_PORTABLE_BLOCK_SIZE) {
        int j_end = MinInt(jj + MATRIX_PORTABLE_BLOCK_SIZE, B->column);
        for (int ii = 0; ii < A->row; ii += MATRIX_PORTABLE_BLOCK_SIZE) {
            int i_end = MinInt(ii + MATRIX_PORTABLE_BLOCK_SIZE, A->row);
            int i = ii;
            for (; i + 3 < i_end; i += 4) {
                int j = jj;
                for (; j + 3 < j_end; j += 4) {
                    Multiply4x4(A->data + i * A->column, A->column,
                                B->data + j, B->column,
                                C->data + i * C->column + j, C->column,
                                A->column);
                }
                if (j < j_end) {
                    MultiplyTail(A, B, C, i, i + 4, j, j_end);
                }
            }
            if (i < i_end) {
                MultiplyTail(A, B, C, i, i_end, jj, j_end);
            }
        }
    }
    return MATRIX_SUCCESS;
}

MatrixError MatrixMultiply(const Matrix *A, const Matrix *B, Matrix *C)
{
    return MatrixMultiplyOptimized(A, B, C);
}
