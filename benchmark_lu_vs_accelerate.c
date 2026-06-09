#include "matrix_core.h"
#include "matrix_lu.h"

#include <Accelerate/Accelerate.h>

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static double NowSeconds(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

static int CompareDouble(const void *left, const void *right)
{
    double a = *(const double *)left;
    double b = *(const double *)right;
    return (a > b) - (a < b);
}

static double Median(double *values, int count)
{
    qsort(values, (size_t)count, sizeof(*values), CompareDouble);
    return values[count / 2];
}

static double InputValue(int i, int j, int n)
{
    if (i == j) {
        return (double)n + 1.0;
    }
    uint32_t hash = (uint32_t)(i + 1) * 0x9e3779b9u;
    hash ^= (uint32_t)(j + 1) * 0x85ebca6bu;
    hash ^= hash >> 16;
    return ((double)(hash & 0xffffu) / 65535.0) - 0.5;
}

static void FillInputs(double *row_major, double *column_major, int n)
{
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            double value = InputValue(i, j, n);
            row_major[(size_t)i * (size_t)n + (size_t)j] = value;
            column_major[(size_t)i + (size_t)j * (size_t)n] = value;
        }
    }
}

static double InputScale(const double *row_major, int n)
{
    double scale = 1.0;
    size_t count = (size_t)n * (size_t)n;
    for (size_t k = 0; k < count; ++k) {
        double value = fabs(row_major[k]);
        if (value > scale) {
            scale = value;
        }
    }
    return scale;
}

static double CustomResidual(const double *input, const Matrix *L, const Matrix *U)
{
    int n = L->row;
    double max_error = 0.0;
    double scale = InputScale(input, n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            double actual = 0.0;
            int end = i < j ? i : j;
            for (int k = 0; k <= end; ++k) {
                actual += L->data[(size_t)i * (size_t)n + (size_t)k]
                        * U->data[(size_t)k * (size_t)n + (size_t)j];
            }
            double error = fabs(actual - input[(size_t)i * (size_t)n + (size_t)j]);
            if (error > max_error) {
                max_error = error;
            }
        }
    }
    return max_error / scale;
}

static void SwapRows(double *matrix, int n, int first, int second)
{
    if (first == second) {
        return;
    }
    for (int j = 0; j < n; ++j) {
        size_t a = (size_t)first * (size_t)n + (size_t)j;
        size_t b = (size_t)second * (size_t)n + (size_t)j;
        double temporary = matrix[a];
        matrix[a] = matrix[b];
        matrix[b] = temporary;
    }
}

static double LapackResidual(const double *input, double *permuted_input,
                             const double *factors, const __LAPACK_int *pivots, int n)
{
    memcpy(permuted_input, input, (size_t)n * (size_t)n * sizeof(*input));
    for (int k = 0; k < n; ++k) {
        SwapRows(permuted_input, n, k, (int)pivots[k] - 1);
    }

    double max_error = 0.0;
    double scale = InputScale(input, n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            double actual = 0.0;
            int end = i < j ? i : j;
            for (int k = 0; k <= end; ++k) {
                double l = i == k ? 1.0 : factors[(size_t)i + (size_t)k * (size_t)n];
                double u = factors[(size_t)k + (size_t)j * (size_t)n];
                actual += l * u;
            }
            double error = fabs(actual - permuted_input[(size_t)i * (size_t)n + (size_t)j]);
            if (error > max_error) {
                max_error = error;
            }
        }
    }
    return max_error / scale;
}

static int RunsForSize(int n)
{
    if (n <= 128) {
        return 9;
    }
    if (n <= 256) {
        return 7;
    }
    if (n <= 512) {
        return 5;
    }
    return 3;
}

static int BenchmarkSize(int n)
{
    size_t count = (size_t)n * (size_t)n;
    size_t bytes = count * sizeof(double);
    double *input_row_major = (double *)malloc(bytes);
    double *input_column_major = (double *)malloc(bytes);
    double *lapack_work = (double *)malloc(bytes);
    double *permuted_input = (double *)malloc(bytes);
    __LAPACK_int *pivots = (__LAPACK_int *)malloc((size_t)n * sizeof(*pivots));
    if (input_row_major == NULL || input_column_major == NULL || lapack_work == NULL
            || permuted_input == NULL || pivots == NULL) {
        fprintf(stderr, "allocation failed for n=%d\n", n);
        free(input_row_major);
        free(input_column_major);
        free(lapack_work);
        free(permuted_input);
        free(pivots);
        return 0;
    }

    Matrix A, L, U;
    MatrixInit(&A);
    MatrixInit(&L);
    MatrixInit(&U);
    if (MatrixCreate(&A, n, n) != MATRIX_SUCCESS
            || MatrixCreate(&L, n, n) != MATRIX_SUCCESS
            || MatrixCreate(&U, n, n) != MATRIX_SUCCESS) {
        fprintf(stderr, "MatrixCreate failed for n=%d\n", n);
        MatrixFree(&A);
        MatrixFree(&L);
        MatrixFree(&U);
        free(input_row_major);
        free(input_column_major);
        free(lapack_work);
        free(permuted_input);
        free(pivots);
        return 0;
    }

    FillInputs(input_row_major, input_column_major, n);

    memcpy(A.data, input_row_major, bytes);
    MatrixError custom_error = LUDecomposeNoPivot(&A, &L, &U, 1e-12);

    __LAPACK_int size = (__LAPACK_int)n;
    __LAPACK_int info = 0;
    memcpy(lapack_work, input_column_major, bytes);
    dgetrf_(&size, &size, lapack_work, &size, pivots, &info);
    if (custom_error != MATRIX_SUCCESS || info != 0) {
        fprintf(stderr, "warm-up factorization failed for n=%d: custom=%s, lapack=%d\n",
                n, MatrixErrorMessage(custom_error), (int)info);
        MatrixFree(&A);
        MatrixFree(&L);
        MatrixFree(&U);
        free(input_row_major);
        free(input_column_major);
        free(lapack_work);
        free(permuted_input);
        free(pivots);
        return 0;
    }

    int runs = RunsForSize(n);
    double custom_times[9];
    double lapack_times[9];
    for (int run = 0; run < runs; ++run) {
        memcpy(A.data, input_row_major, bytes);
        double start = NowSeconds();
        custom_error = LUDecomposeNoPivot(&A, &L, &U, 1e-12);
        custom_times[run] = NowSeconds() - start;
        if (custom_error != MATRIX_SUCCESS) {
            fprintf(stderr, "custom factorization failed for n=%d: %s\n",
                    n, MatrixErrorMessage(custom_error));
            return 0;
        }

        memcpy(lapack_work, input_column_major, bytes);
        start = NowSeconds();
        dgetrf_(&size, &size, lapack_work, &size, pivots, &info);
        lapack_times[run] = NowSeconds() - start;
        if (info != 0) {
            fprintf(stderr, "LAPACK factorization failed for n=%d: info=%d\n", n, (int)info);
            return 0;
        }
    }

    double custom_residual = CustomResidual(input_row_major, &L, &U);
    double lapack_residual = LapackResidual(input_row_major, permuted_input, lapack_work, pivots, n);
    double custom_seconds = Median(custom_times, runs);
    double lapack_seconds = Median(lapack_times, runs);
    double flops = (2.0 / 3.0) * (double)n * (double)n * (double)n;

    printf("%d,%d,%.6f,%.6f,%.2f,%.3f,%.3f,%.3e,%.3e\n",
           n, runs, custom_seconds * 1e3, lapack_seconds * 1e3,
           custom_seconds / lapack_seconds,
           flops / custom_seconds / 1e9, flops / lapack_seconds / 1e9,
           custom_residual, lapack_residual);

    MatrixFree(&A);
    MatrixFree(&L);
    MatrixFree(&U);
    free(input_row_major);
    free(input_column_major);
    free(lapack_work);
    free(permuted_input);
    free(pivots);
    return custom_residual < 1e-10 && lapack_residual < 1e-10;
}

int main(void)
{
    static const int sizes[] = {32, 64, 128, 256, 512, 768, 1024};
    int success = 1;

    printf("LU benchmark: repository LUDecomposeNoPivot vs Apple Accelerate LAPACK dgetrf\n");
    printf("Input copying and allocation are excluded from timed regions.\n");
    printf("n,runs,custom_ms,lapack_ms,lapack_speedup,custom_gflops,lapack_gflops,custom_residual,lapack_residual\n");
    for (size_t k = 0; k < sizeof(sizes) / sizeof(sizes[0]); ++k) {
        if (!BenchmarkSize(sizes[k])) {
            success = 0;
            break;
        }
    }
    return success ? 0 : 1;
}
