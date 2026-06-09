#include "matrix_core.h"
#include "matrix_ops.h"

#include <math.h>
#include <stdint.h>
#include <stdio.h>

static void FillPattern(Matrix *A, uint32_t seed)
{
    int count = A->row * A->column;
    for (int k = 0; k < count; ++k) {
        uint32_t hash = (uint32_t)k * 0x9e3779b9u + seed;
        hash ^= hash >> 16;
        hash *= 0x85ebca6bu;
        hash ^= hash >> 13;
        A->data[k] = ((double)(hash & 0xffffu) / 32767.5) - 1.0;
    }
}

static int MatricesMatch(const Matrix *actual, const Matrix *expected)
{
    int count = actual->row * actual->column;
    double scale = 1.0;
    double max_error = 0.0;
    for (int k = 0; k < count; ++k) {
        double expected_abs = fabs(expected->data[k]);
        double error = fabs(actual->data[k] - expected->data[k]);
        if (expected_abs > scale) {
            scale = expected_abs;
        }
        if (error > max_error) {
            max_error = error;
        }
    }
    return max_error / scale < 1e-11;
}

static int RunMultiplyCase(int rows, int inner, int columns)
{
    Matrix A, B, expected, actual;
    MatrixInit(&A);
    MatrixInit(&B);
    MatrixInit(&expected);
    MatrixInit(&actual);
    if (MatrixCreate(&A, rows, inner) != MATRIX_SUCCESS
            || MatrixCreate(&B, inner, columns) != MATRIX_SUCCESS
            || MatrixCreate(&expected, rows, columns) != MATRIX_SUCCESS
            || MatrixCreate(&actual, rows, columns) != MATRIX_SUCCESS) {
        fprintf(stderr, "allocation failed for %d x %d times %d x %d\n",
                rows, inner, inner, columns);
        return 0;
    }

    FillPattern(&A, 0x12345678u);
    FillPattern(&B, 0x87654321u);
    if (MatrixMultiplyOriginal(&A, &B, &expected) != MATRIX_SUCCESS
            || MatrixMultiply(&A, &B, &actual) != MATRIX_SUCCESS
            || !MatricesMatch(&actual, &expected)) {
        fprintf(stderr, "multiply mismatch for %d x %d times %d x %d\n",
                rows, inner, inner, columns);
        return 0;
    }

    MatrixFree(&A);
    MatrixFree(&B);
    MatrixFree(&expected);
    MatrixFree(&actual);
    return 1;
}

int main(void)
{
    static const int cases[][3] = {
        {1, 1, 1},
        {3, 5, 7},
        {7, 11, 9},
        {13, 17, 15},
        {31, 33, 129},
        {35, 127, 137},
        {193, 211, 197}
    };

    for (size_t k = 0; k < sizeof(cases) / sizeof(cases[0]); ++k) {
        if (!RunMultiplyCase(cases[k][0], cases[k][1], cases[k][2])) {
            return 1;
        }
    }
    printf("portable optimized multiplication regression cases passed\n");
    return 0;
}
