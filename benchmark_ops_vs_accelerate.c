#include "matrix_core.h"
#include "matrix_ops.h"

#include <Accelerate/Accelerate.h>

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

typedef int (*BenchmarkFunction)(void *context);

typedef struct {
    const Matrix *A;
    const Matrix *B;
    Matrix *C;
    REAL alpha;
    REAL norm;
} OperationContext;

static volatile double benchmark_guard;

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

static size_t MatrixElementCount(const Matrix *A)
{
    return (size_t)A->row * (size_t)A->column;
}

static void FillPattern(Matrix *A, uint32_t seed)
{
    size_t count = MatrixElementCount(A);
    for (size_t k = 0; k < count; ++k) {
        uint32_t hash = (uint32_t)k * 0x9e3779b9u + seed;
        hash ^= hash >> 16;
        hash *= 0x85ebca6bu;
        hash ^= hash >> 13;
        A->data[k] = ((double)(hash & 0xffffu) / 32767.5) - 1.0;
    }
}

static double MatrixChecksum(const Matrix *A)
{
    size_t count = MatrixElementCount(A);
    return A->data[0] + A->data[count / 2] + A->data[count - 1];
}

static double RelativeDifference(const Matrix *actual, const Matrix *expected)
{
    size_t count = MatrixElementCount(actual);
    double max_error = 0.0;
    double scale = 1.0;
    for (size_t k = 0; k < count; ++k) {
        double expected_abs = fabs(expected->data[k]);
        double error = fabs(actual->data[k] - expected->data[k]);
        if (expected_abs > scale) {
            scale = expected_abs;
        }
        if (error > max_error) {
            max_error = error;
        }
    }
    return max_error / scale;
}

static double TimeMedian(BenchmarkFunction function, void *context, int trials)
{
    const double target_seconds = 0.015;
    const int max_repeats = 100000;
    double samples[7];

    if (!function(context)) {
        return -1.0;
    }

    double start = NowSeconds();
    if (!function(context)) {
        return -1.0;
    }
    double elapsed = NowSeconds() - start;
    int repeats = elapsed > 0.0 ? (int)ceil(target_seconds / elapsed) : max_repeats;
    if (repeats < 1) {
        repeats = 1;
    } else if (repeats > max_repeats) {
        repeats = max_repeats;
    }

    for (int trial = 0; trial < trials; ++trial) {
        start = NowSeconds();
        for (int repeat = 0; repeat < repeats; ++repeat) {
            if (!function(context)) {
                return -1.0;
            }
        }
        samples[trial] = (NowSeconds() - start) / (double)repeats;
    }
    return Median(samples, trials);
}

static int CustomAdd(void *context)
{
    OperationContext *op = (OperationContext *)context;
    return MatrixAdd(op->A, op->B, op->C) == MATRIX_SUCCESS;
}

static int AccelerateAdd(void *context)
{
    OperationContext *op = (OperationContext *)context;
    vDSP_vaddD(op->A->data, 1, op->B->data, 1, op->C->data, 1,
               (vDSP_Length)MatrixElementCount(op->A));
    return 1;
}

static int CustomSub(void *context)
{
    OperationContext *op = (OperationContext *)context;
    return MatrixSub(op->A, op->B, op->C) == MATRIX_SUCCESS;
}

static int AccelerateSub(void *context)
{
    OperationContext *op = (OperationContext *)context;
    vDSP_vsubD(op->B->data, 1, op->A->data, 1, op->C->data, 1,
               (vDSP_Length)MatrixElementCount(op->A));
    return 1;
}

static int CustomScale(void *context)
{
    OperationContext *op = (OperationContext *)context;
    return MatrixScale(op->alpha, op->A, op->C) == MATRIX_SUCCESS;
}

static int AccelerateScale(void *context)
{
    OperationContext *op = (OperationContext *)context;
    vDSP_vsmulD(op->A->data, 1, &op->alpha, op->C->data, 1,
                (vDSP_Length)MatrixElementCount(op->A));
    return 1;
}

static int CustomTranspose(void *context)
{
    OperationContext *op = (OperationContext *)context;
    return MatrixTranspose(op->A, op->C) == MATRIX_SUCCESS;
}

static int AccelerateTranspose(void *context)
{
    OperationContext *op = (OperationContext *)context;
    vDSP_mtransD(op->A->data, 1, op->C->data, 1,
                 (vDSP_Length)op->A->column, (vDSP_Length)op->A->row);
    return 1;
}

static int CustomNorm(void *context)
{
    OperationContext *op = (OperationContext *)context;
    return MatrixNormFrobenius(op->A, &op->norm) == MATRIX_SUCCESS;
}

static int AccelerateNorm(void *context)
{
    OperationContext *op = (OperationContext *)context;
    op->norm = cblas_dnrm2((int)MatrixElementCount(op->A), op->A->data, 1);
    return 1;
}

static int CustomMultiplyIJK(void *context)
{
    OperationContext *op = (OperationContext *)context;
    return MatrixMultiplyOriginal(op->A, op->B, op->C) == MATRIX_SUCCESS;
}

static int CustomMultiplyIKJ(void *context)
{
    OperationContext *op = (OperationContext *)context;
    return MatrixMultiplyIKJ(op->A, op->B, op->C) == MATRIX_SUCCESS;
}

static int CustomMultiplyBlocked(void *context)
{
    OperationContext *op = (OperationContext *)context;
    return MatrixMultiplyBlocked(op->A, op->B, op->C) == MATRIX_SUCCESS;
}

static int CustomMultiplyOptimized(void *context)
{
    OperationContext *op = (OperationContext *)context;
    return MatrixMultiply(op->A, op->B, op->C) == MATRIX_SUCCESS;
}

static int AccelerateMultiply(void *context)
{
    OperationContext *op = (OperationContext *)context;
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                op->A->row, op->B->column, op->A->column,
                1.0, op->A->data, op->A->column,
                op->B->data, op->B->column,
                0.0, op->C->data, op->C->column);
    return 1;
}

static void FreeMatrices(Matrix *A, Matrix *B, Matrix *C, Matrix *D)
{
    MatrixFree(A);
    MatrixFree(B);
    MatrixFree(C);
    MatrixFree(D);
}

static int CreateMatrices(Matrix *A, Matrix *B, Matrix *C, Matrix *D, int n)
{
    MatrixInit(A);
    MatrixInit(B);
    MatrixInit(C);
    MatrixInit(D);
    if (MatrixCreate(A, n, n) != MATRIX_SUCCESS
            || MatrixCreate(B, n, n) != MATRIX_SUCCESS
            || MatrixCreate(C, n, n) != MATRIX_SUCCESS
            || MatrixCreate(D, n, n) != MATRIX_SUCCESS) {
        FreeMatrices(A, B, C, D);
        return 0;
    }
    return 1;
}

static int ValidateVectorOperation(const char *name, BenchmarkFunction custom,
                                   BenchmarkFunction accelerate,
                                   OperationContext *custom_context,
                                   OperationContext *accelerate_context)
{
    if (!custom(custom_context) || !accelerate(accelerate_context)) {
        fprintf(stderr, "%s failed during validation\n", name);
        return 0;
    }
    double difference = RelativeDifference(custom_context->C, accelerate_context->C);
    if (difference > 1e-12) {
        fprintf(stderr, "%s validation failed: relative difference %.3e\n", name, difference);
        return 0;
    }
    return 1;
}

static int BenchmarkVectorSize(int n)
{
    Matrix A, B, custom_output, accelerate_output;
    if (!CreateMatrices(&A, &B, &custom_output, &accelerate_output, n)) {
        fprintf(stderr, "allocation failed for vector benchmark n=%d\n", n);
        return 0;
    }

    FillPattern(&A, 0x12345678u);
    FillPattern(&B, 0x87654321u);
    OperationContext custom = {&A, &B, &custom_output, 2.5, 0.0};
    OperationContext accelerate = {&A, &B, &accelerate_output, 2.5, 0.0};

    if (!ValidateVectorOperation("add", CustomAdd, AccelerateAdd, &custom, &accelerate)
            || !ValidateVectorOperation("subtract", CustomSub, AccelerateSub, &custom, &accelerate)
            || !ValidateVectorOperation("scale", CustomScale, AccelerateScale, &custom, &accelerate)
            || !ValidateVectorOperation("transpose", CustomTranspose, AccelerateTranspose,
                                        &custom, &accelerate)
            || !CustomNorm(&custom) || !AccelerateNorm(&accelerate)) {
        FreeMatrices(&A, &B, &custom_output, &accelerate_output);
        return 0;
    }
    double norm_difference = fabs(custom.norm - accelerate.norm) / fmax(1.0, fabs(accelerate.norm));
    if (norm_difference > 1e-12) {
        fprintf(stderr, "norm validation failed: relative difference %.3e\n", norm_difference);
        FreeMatrices(&A, &B, &custom_output, &accelerate_output);
        return 0;
    }

    double custom_add = TimeMedian(CustomAdd, &custom, 7);
    double accelerate_add = TimeMedian(AccelerateAdd, &accelerate, 7);
    double custom_sub = TimeMedian(CustomSub, &custom, 7);
    double accelerate_sub = TimeMedian(AccelerateSub, &accelerate, 7);
    double custom_scale = TimeMedian(CustomScale, &custom, 7);
    double accelerate_scale = TimeMedian(AccelerateScale, &accelerate, 7);
    double custom_transpose = TimeMedian(CustomTranspose, &custom, 7);
    double accelerate_transpose = TimeMedian(AccelerateTranspose, &accelerate, 7);
    double custom_norm = TimeMedian(CustomNorm, &custom, 7);
    double accelerate_norm = TimeMedian(AccelerateNorm, &accelerate, 7);

    benchmark_guard += MatrixChecksum(&custom_output) + MatrixChecksum(&accelerate_output)
                     + custom.norm + accelerate.norm;
    printf("%d,%.6f,%.6f,%.2f,%.6f,%.6f,%.2f,%.6f,%.6f,%.2f,"
           "%.6f,%.6f,%.2f,%.6f,%.6f,%.2f\n",
           n,
           custom_add * 1e3, accelerate_add * 1e3, custom_add / accelerate_add,
           custom_sub * 1e3, accelerate_sub * 1e3, custom_sub / accelerate_sub,
           custom_scale * 1e3, accelerate_scale * 1e3, custom_scale / accelerate_scale,
           custom_transpose * 1e3, accelerate_transpose * 1e3,
           custom_transpose / accelerate_transpose,
           custom_norm * 1e3, accelerate_norm * 1e3, custom_norm / accelerate_norm);

    FreeMatrices(&A, &B, &custom_output, &accelerate_output);
    return 1;
}

static int ValidateMultiply(const char *name, BenchmarkFunction custom,
                            OperationContext *custom_context,
                            OperationContext *accelerate_context)
{
    if (!custom(custom_context)) {
        fprintf(stderr, "%s failed during validation\n", name);
        return 0;
    }
    double difference = RelativeDifference(custom_context->C, accelerate_context->C);
    if (difference > 1e-11) {
        fprintf(stderr, "%s validation failed: relative difference %.3e\n", name, difference);
        return 0;
    }
    return 1;
}

static int BenchmarkMultiplySize(int n)
{
    Matrix A, B, custom_output, accelerate_output;
    if (!CreateMatrices(&A, &B, &custom_output, &accelerate_output, n)) {
        fprintf(stderr, "allocation failed for multiplication benchmark n=%d\n", n);
        return 0;
    }

    FillPattern(&A, 0x12345678u);
    FillPattern(&B, 0x87654321u);
    OperationContext custom = {&A, &B, &custom_output, 0.0, 0.0};
    OperationContext accelerate = {&A, &B, &accelerate_output, 0.0, 0.0};

    if (!AccelerateMultiply(&accelerate)
            || !ValidateMultiply("multiply i-j-k", CustomMultiplyIJK, &custom, &accelerate)
            || !ValidateMultiply("multiply i-k-j", CustomMultiplyIKJ, &custom, &accelerate)
            || !ValidateMultiply("multiply blocked", CustomMultiplyBlocked, &custom, &accelerate)
            || !ValidateMultiply("multiply optimized", CustomMultiplyOptimized,
                                 &custom, &accelerate)) {
        FreeMatrices(&A, &B, &custom_output, &accelerate_output);
        return 0;
    }

    double ijk = TimeMedian(CustomMultiplyIJK, &custom, 5);
    double ikj = TimeMedian(CustomMultiplyIKJ, &custom, 5);
    double blocked = TimeMedian(CustomMultiplyBlocked, &custom, 5);
    double optimized = TimeMedian(CustomMultiplyOptimized, &custom, 5);
    double blas = TimeMedian(AccelerateMultiply, &accelerate, 5);
    double flops = 2.0 * (double)n * (double)n * (double)n;

    benchmark_guard += MatrixChecksum(&custom_output) + MatrixChecksum(&accelerate_output);
    printf("%d,%.6f,%.6f,%.6f,%.6f,%.6f,%.2f,%.2f,%.2f,%.3f,%.3f,%.3f\n",
           n, ijk * 1e3, ikj * 1e3, blocked * 1e3, optimized * 1e3, blas * 1e3,
           optimized / blas, blocked / optimized, ijk / optimized,
           flops / blocked / 1e9, flops / optimized / 1e9, flops / blas / 1e9);

    FreeMatrices(&A, &B, &custom_output, &accelerate_output);
    return 1;
}

int main(void)
{
    static const int vector_sizes[] = {256, 512, 1024, 2048};
    static const int multiply_sizes[] = {64, 128, 256, 512, 768};

    printf("Matrix operations benchmark: repository matrix_ops vs Apple Accelerate\n");
    printf("Initialization is excluded from timed regions. Values are medians in milliseconds.\n");
    printf("\nElementwise operations: Accelerate uses vDSP; norm uses BLAS cblas_dnrm2.\n");
    printf("n,add_custom_ms,add_accelerate_ms,add_speedup,"
           "sub_custom_ms,sub_accelerate_ms,sub_speedup,"
           "scale_custom_ms,scale_accelerate_ms,scale_speedup,"
           "transpose_custom_ms,transpose_accelerate_ms,transpose_speedup,"
           "norm_custom_ms,norm_blas_ms,norm_speedup\n");
    for (size_t k = 0; k < sizeof(vector_sizes) / sizeof(vector_sizes[0]); ++k) {
        if (!BenchmarkVectorSize(vector_sizes[k])) {
            return 1;
        }
    }

    printf("\nMatrix multiplication: Accelerate uses BLAS cblas_dgemm.\n");
    printf("n,ijk_ms,ikj_ms,blocked_ms,optimized_ms,blas_ms,blas_vs_optimized_speedup,"
           "optimized_vs_blocked_speedup,optimized_vs_ijk_speedup,"
           "blocked_gflops,optimized_gflops,blas_gflops\n");
    for (size_t k = 0; k < sizeof(multiply_sizes) / sizeof(multiply_sizes[0]); ++k) {
        if (!BenchmarkMultiplySize(multiply_sizes[k])) {
            return 1;
        }
    }

    printf("\nvalidation_guard=%.6f\n", benchmark_guard);
    return 0;
}
