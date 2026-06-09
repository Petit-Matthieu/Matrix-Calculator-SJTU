#include "matrix_core.h"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

void MatrixInit(Matrix *A)
{
    if (A == NULL) {
        return;
    }
    A->row = 0;
    A->column = 0;
    A->data = NULL;
}

int MatrixIsValid(const Matrix *A)
{
    return (A != NULL && A->row > 0 && A->column > 0 && A->data != NULL);
}

int MatrixHasShape(const Matrix *A, int row, int column)
{
    return MatrixIsValid(A) && A->row == row && A->column == column;
}

int MatrixIndex(const Matrix *A, int i, int j)
{
    return i * A->column + j;
}

MatrixError MatrixCreate(Matrix *A, int row, int column)
{
    if (A == NULL) {
        return MATRIX_ERROR_NULL_POINTER;
    }
    if (A->data != NULL) {
        return MATRIX_ERROR_ALREADY_ALLOCATED;
    }
    if (row <= 0 || column <= 0) {
        return MATRIX_ERROR_INVALID_SIZE;
    }

    size_t nrow = (size_t)row;
    size_t ncol = (size_t)column;
    if (nrow > SIZE_MAX / ncol || nrow * ncol > SIZE_MAX / sizeof(REAL)) {
        return MATRIX_ERROR_SIZE_OVERFLOW;
    }

    A->data = (REAL *)malloc(nrow * ncol * sizeof(REAL));
    if (A->data == NULL) {
        MatrixInit(A);
        return MATRIX_ERROR_ALLOC_FAILED;
    }
    A->row = row;
    A->column = column;
    return MATRIX_SUCCESS;
}

void MatrixFree(Matrix *A)
{
    if (A == NULL) {
        return;
    }
    free(A->data);
    A->data = NULL;
    A->row = 0;
    A->column = 0;
}

MatrixError MatrixSet(Matrix *A, int i, int j, REAL value)
{
    if (!MatrixIsValid(A)) {
        return MATRIX_ERROR_NULL_POINTER;
    }
    if (i < 0 || i >= A->row || j < 0 || j >= A->column) {
        return MATRIX_ERROR_INDEX_OUT_OF_RANGE;
    }
    A->data[MatrixIndex(A, i, j)] = value;
    return MATRIX_SUCCESS;
}

MatrixError MatrixGet(const Matrix *A, int i, int j, REAL *value)
{
    if (!MatrixIsValid(A) || value == NULL) {
        return MATRIX_ERROR_NULL_POINTER;
    }
    if (i < 0 || i >= A->row || j < 0 || j >= A->column) {
        return MATRIX_ERROR_INDEX_OUT_OF_RANGE;
    }
    *value = A->data[MatrixIndex(A, i, j)];
    return MATRIX_SUCCESS;
}

MatrixError MatrixFillZero(Matrix *A)
{
    if (!MatrixIsValid(A)) {
        return MATRIX_ERROR_NULL_POINTER;
    }
    int total = A->row * A->column;
    for (int k = 0; k < total; ++k) {
        A->data[k] = 0.0;
    }
    return MATRIX_SUCCESS;
}

MatrixError MatrixFillSequence(Matrix *A, REAL start, REAL step)
{
    if (!MatrixIsValid(A)) {
        return MATRIX_ERROR_NULL_POINTER;
    }
    int total = A->row * A->column;
    for (int k = 0; k < total; ++k) {
        A->data[k] = start + step * (REAL)k;
    }
    return MATRIX_SUCCESS;
}

MatrixError MatrixCopy(const Matrix *src, Matrix *dst)
{
    if (!MatrixIsValid(src) || !MatrixIsValid(dst)) {
        return MATRIX_ERROR_NULL_POINTER;
    }
    if (src->row != dst->row || src->column != dst->column) {
        return MATRIX_ERROR_SIZE_MISMATCH;
    }
    int total = src->row * src->column;
    for (int k = 0; k < total; ++k) {
        dst->data[k] = src->data[k];
    }
    return MATRIX_SUCCESS;
}

void MatrixPrint(const Matrix *A, const char *name)
{
    if (name == NULL) {
        name = "Matrix";
    }
    if (!MatrixIsValid(A)) {
        printf("%s is an empty or invalid matrix.\n", name);
        return;
    }
    printf("%s = (%d x %d)\n", name, A->row, A->column);
    for (int i = 0; i < A->row; ++i) {
        printf("  ");
        for (int j = 0; j < A->column; ++j) {
            printf("%10.4f ", A->data[MatrixIndex(A, i, j)]);
        }
        printf("\n");
    }
}

void MatrixPrintIndexMap(const Matrix *A, const char *name)
{
    if (name == NULL) {
        name = "Matrix";
    }
    if (!MatrixIsValid(A)) {
        printf("%s is an empty or invalid matrix.\n", name);
        return;
    }
    printf("Index map for %s:\n", name);
    for (int i = 0; i < A->row; ++i) {
        for (int j = 0; j < A->column; ++j) {
            printf("  %s(%d,%d) -> data[%d]\n", name, i, j, MatrixIndex(A, i, j));
        }
    }
}

const char *MatrixErrorMessage(MatrixError error)
{
    switch (error) {
        case MATRIX_SUCCESS: return "success";
        case MATRIX_ERROR_NULL_POINTER: return "null pointer or invalid matrix";
        case MATRIX_ERROR_INVALID_SIZE: return "invalid matrix size";
        case MATRIX_ERROR_SIZE_OVERFLOW: return "matrix size overflow";
        case MATRIX_ERROR_ALLOC_FAILED: return "memory allocation failed";
        case MATRIX_ERROR_INDEX_OUT_OF_RANGE: return "index out of range";
        case MATRIX_ERROR_SIZE_MISMATCH: return "matrix size mismatch";
        case MATRIX_ERROR_NOT_SQUARE: return "matrix must be square";
        case MATRIX_ERROR_SINGULAR: return "matrix is singular or nearly singular";
        case MATRIX_ERROR_ALREADY_ALLOCATED: return "matrix data already allocated; free it first";
        default: return "unknown error";
    }
}
