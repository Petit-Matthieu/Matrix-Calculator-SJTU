#ifndef MATRIX_CORE_H
#define MATRIX_CORE_H

#include <stddef.h>

typedef double REAL;

typedef enum {
    MATRIX_SUCCESS = 0,
    MATRIX_ERROR_NULL_POINTER = 1,
    MATRIX_ERROR_INVALID_SIZE = 2,
    MATRIX_ERROR_SIZE_OVERFLOW = 3,
    MATRIX_ERROR_ALLOC_FAILED = 4,
    MATRIX_ERROR_INDEX_OUT_OF_RANGE = 5,
    MATRIX_ERROR_SIZE_MISMATCH = 6,
    MATRIX_ERROR_NOT_SQUARE = 7,
    MATRIX_ERROR_SINGULAR = 8,
    MATRIX_ERROR_ALREADY_ALLOCATED = 9
} MatrixError;

typedef struct {
    int row;
    int column;
    REAL *data;
} Matrix;

void MatrixInit(Matrix *A);
MatrixError MatrixCreate(Matrix *A, int row, int column);
void MatrixFree(Matrix *A);
int MatrixIsValid(const Matrix *A);
int MatrixHasShape(const Matrix *A, int row, int column);
int MatrixIndex(const Matrix *A, int i, int j);
MatrixError MatrixSet(Matrix *A, int i, int j, REAL value);
MatrixError MatrixGet(const Matrix *A, int i, int j, REAL *value);
MatrixError MatrixFillZero(Matrix *A);
MatrixError MatrixFillSequence(Matrix *A, REAL start, REAL step);
MatrixError MatrixCopy(const Matrix *src, Matrix *dst);
void MatrixPrint(const Matrix *A, const char *name);
void MatrixPrintIndexMap(const Matrix *A, const char *name);
const char *MatrixErrorMessage(MatrixError error);

#endif
