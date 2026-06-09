# Methods Explained

[English](#english) | [中文](#中文)

This document explains the algorithms and implementation methods currently used
in the Matrix Calculator SJTU repository.

## English

### 1. Matrix Representation

The C library stores every matrix in the `Matrix` structure:

```c
typedef double REAL;

typedef struct {
    int row;
    int column;
    REAL *data;
} Matrix;
```

The matrix data is stored in row-major order. For an element `A(i, j)`, the
linear index is:

```c
i * A->column + j
```

This storage method is simple, matches normal C array layout, and makes each row
contiguous in memory. Row-contiguous storage is important for the optimized
matrix multiplication routines because reading or writing consecutive memory
locations is much faster than jumping through memory with a large stride.

### 2. Matrix Lifecycle and Error Handling

The basic lifecycle is:

1. Call `MatrixInit()` to put the matrix into an empty safe state.
2. Call `MatrixCreate()` to allocate `row * column` elements.
3. Use the matrix through `MatrixSet()`, `MatrixGet()`, or direct `data` access.
4. Call `MatrixFree()` to release memory.

Most public functions return a `MatrixError` value. This avoids silent failure
and makes invalid input easier to diagnose. The implementation checks common
errors such as:

- Null or invalid matrix pointers
- Invalid matrix sizes
- Memory allocation failure
- Index out of range
- Shape mismatch
- Non-square input for LU decomposition
- Singular or nearly singular pivots

### 3. Basic Matrix Operations

The basic operations are implemented with direct loops over the row-major data.

#### Addition and Subtraction

`MatrixAdd()` and `MatrixSub()` require `A`, `B`, and `C` to have the same
shape. They then scan all elements linearly:

```c
C[k] = A[k] + B[k]
C[k] = A[k] - B[k]
```

This is cache-friendly because all three arrays are read or written in
contiguous order.

#### Scalar Multiplication

`MatrixScale()` multiplies every element by one scalar:

```c
B[k] = alpha * A[k]
```

It also scans memory linearly, so its performance is mainly limited by memory
bandwidth.

#### Transpose

`MatrixTranspose()` writes:

```c
AT(j, i) = A(i, j)
```

The input matrix is read row by row, but the output matrix is written with a
stride. This is why transpose is usually slower than addition or scalar
multiplication.

#### Frobenius Norm

`MatrixNormFrobenius()` computes:

```text
sqrt(sum(A(i, j)^2))
```

The implementation performs one linear pass over the matrix, accumulates the
sum of squares, and then calls `sqrt()` from the standard math library.

### 4. Matrix Multiplication Methods

Matrix multiplication dominates the performance work in this repository. For
matrices:

- `A` with shape `m x k`
- `B` with shape `k x n`
- `C` with shape `m x n`

the mathematical definition is:

```text
C(i, j) = sum(A(i, p) * B(p, j)), p = 0 ... k - 1
```

The repository keeps several implementations so their performance can be
compared.

#### 4.1 Original `i-j-k` Method

`MatrixMultiplyOriginal()` uses the most direct triple loop:

```text
for i
  for j
    for k
      C(i, j) += A(i, k) * B(k, j)
```

This is easy to understand and directly follows the mathematical formula.
However, it has poor cache behavior for row-major matrices. `A(i, k)` is read
contiguously, but `B(k, j)` jumps by one full row each time `k` changes.

This method is retained as a correctness reference and baseline.

#### 4.2 Reordered `i-k-j` Method

`MatrixMultiplyIKJ()` changes the loop order:

```text
for i
  for k
    for j
      C(i, j) += A(i, k) * B(k, j)
```

This improves memory access for row-major storage:

- `B(k, j)` is read across one contiguous row.
- `C(i, j)` is updated across one contiguous row.
- `A(i, k)` is loaded once and reused across many columns.

The arithmetic count is unchanged, but memory access is more efficient.

#### 4.3 Blocked Method

`MatrixMultiplyBlocked()` divides the work into `64 x 64` blocks. Blocking is a
cache optimization. Instead of streaming through entire large matrices at once,
the code works on smaller regions that are more likely to stay in cache.

The blocked version still uses normal scalar C loops. Inside each block, the
inner column loop is manually unrolled by four elements:

```c
c_row[j]     += aik * b_row[j];
c_row[j + 1] += aik * b_row[j + 1];
c_row[j + 2] += aik * b_row[j + 2];
c_row[j + 3] += aik * b_row[j + 3];
```

This reduces loop overhead and gives the compiler a clearer pattern to
optimize.

#### 4.4 Current Optimized Portable-C Method

`MatrixMultiply()` currently delegates to `MatrixMultiplyOptimized()`.

This implementation follows the user's current constraint: it avoids NEON,
external libraries, platform-specific threading APIs, and packed heap buffers.
It uses standard C methods only.

The optimized method uses:

- `64 x 64` output tiling
- A scalar `4 x 4` microkernel
- Local accumulator variables
- C99 `restrict` pointers
- Generic tail handling for rectangular and uneven matrix sizes

The microkernel computes a `4 x 4` output tile at a time. Instead of repeatedly
loading and storing `C` inside the innermost loop, it keeps sixteen partial
sums in local variables:

```c
REAL c00, c01, c02, c03;
REAL c10, c11, c12, c13;
REAL c20, c21, c22, c23;
REAL c30, c31, c32, c33;
```

For every inner index `k`, it loads four values from one row of `B` and one
value from each of four rows of `A`. It updates all sixteen accumulators, then
writes the final result to `C` only once.

This improves performance because:

- The output values stay in registers during accumulation.
- Fewer writes are made to memory.
- Consecutive values of `B` are read from contiguous memory.
- The compiler has a small, regular loop body to optimize.
- The implementation remains portable C.

If the matrix dimensions are not multiples of four, `MultiplyTail()` handles the
remaining rows or columns with a generic scalar loop. This is why the optimized
method works for both square and rectangular matrices.

### 5. LU Decomposition Method

`LUDecomposeNoPivot()` factors a square matrix `A` into:

```text
A = L * U
```

where:

- `L` is lower triangular with ones on the diagonal.
- `U` is upper triangular.

The implementation uses the Doolittle-style no-pivot algorithm:

1. Initialize `L` to the identity matrix.
2. Compute each row of `U`.
3. Check the current pivot `U(k, k)`.
4. Compute the entries below the pivot in `L`.

If the absolute value of a pivot is smaller than the given tolerance, the
function returns `MATRIX_ERROR_SINGULAR`.

This no-pivot method is useful for learning and controlled input data. It is
not as numerically robust as a professional LAPACK factorization with partial
pivoting. Some matrices that are mathematically invertible can still fail or
produce inaccurate results without pivoting.

### 6. Solving Linear Systems

After LU decomposition, the system:

```text
A * x = b
```

is solved in two stages:

1. Forward substitution solves `L * y = b`.
2. Back substitution solves `U * x = y`.

This avoids directly computing the inverse of `A`, which is usually less stable
and less efficient.

### 7. Determinant Method

After no-pivot LU decomposition:

```text
det(A) = det(L) * det(U)
```

Because `L` has ones on the diagonal, `det(L) = 1`. Therefore:

```text
det(A) = product of diagonal elements of U
```

This is implemented by `LUDeterminant()`.

### 8. Testing Method

`test_matrix_ops.c` tests the optimized multiplication implementation by
comparing it with the original triple-loop version.

The test cases include:

- Very small matrices
- Rectangular matrices
- Sizes that are not multiples of four
- Larger irregular dimensions

This is important because the optimized `4 x 4` microkernel has separate tail
handling. The tests verify that the tail path agrees with the simple reference
implementation.

Run the tests with:

```sh
make run-test
```

### 9. Benchmark Method

The repository contains two macOS benchmark programs:

- `benchmark_ops_vs_accelerate.c`
- `benchmark_lu_vs_accelerate.c`

The matrix-operation benchmark compares the repository implementation with
Apple Accelerate vDSP and BLAS. The LU benchmark compares the repository
no-pivot LU implementation with Accelerate LAPACK `dgetrf`.

Important benchmark rules:

- Matrix initialization is excluded from timed regions.
- Each result is reported as a median time.
- Multiplication performance is also reported as GFLOP/s.
- Validation checks are used to prevent the compiler from removing work.

Run benchmarks with:

```sh
make run-benchmark
```

Generate the chart with:

```sh
make plot-benchmark-preview
```

### 10. Why BLAS Is Still Faster

The optimized repository multiplication is much faster than the original C
triple loop, but Apple Accelerate BLAS remains significantly faster.

That is expected. Professional BLAS implementations usually use:

- Architecture-specific SIMD instructions
- Hand-tuned assembly or compiler intrinsics
- Packed matrix panels
- Multithreading
- CPU-cache-size-specific kernels

The current repository implementation intentionally avoids those techniques so
that the production code stays portable and based on standard C.

## 中文

### 1. 矩阵表示方法

C 语言矩阵库使用 `Matrix` 结构体保存矩阵：

```c
typedef double REAL;

typedef struct {
    int row;
    int column;
    REAL *data;
} Matrix;
```

矩阵采用行优先存储。元素 `A(i, j)` 在线性数组中的位置是：

```c
i * A->column + j
```

这种存储方式简单，符合普通 C 数组的使用习惯，并且同一行的数据在内存中
是连续的。连续内存访问对矩阵乘法优化很重要，因为顺序读写内存通常比
跨大步长访问内存更快。

### 2. 矩阵生命周期和错误处理

基本使用流程是：

1. 调用 `MatrixInit()`，把矩阵初始化为空的安全状态。
2. 调用 `MatrixCreate()`，分配 `row * column` 个元素。
3. 通过 `MatrixSet()`、`MatrixGet()` 或直接访问 `data` 使用矩阵。
4. 调用 `MatrixFree()` 释放内存。

大多数公开函数返回 `MatrixError`。这样可以避免静默失败，并且更容易定位
输入错误。当前实现会检查：

- 空指针或无效矩阵
- 非法矩阵大小
- 内存分配失败
- 索引越界
- 矩阵形状不匹配
- LU 分解输入不是方阵
- 主元为零或接近零

### 3. 基本矩阵运算

基本运算直接遍历行优先数组实现。

#### 加法和减法

`MatrixAdd()` 和 `MatrixSub()` 要求 `A`、`B`、`C` 形状完全一致，然后线性
扫描所有元素：

```c
C[k] = A[k] + B[k]
C[k] = A[k] - B[k]
```

由于三个数组都是连续读写，因此缓存表现较好。

#### 标量乘法

`MatrixScale()` 将矩阵中的每个元素乘以同一个标量：

```c
B[k] = alpha * A[k]
```

它也是线性扫描内存，性能主要受内存带宽影响。

#### 转置

`MatrixTranspose()` 执行：

```c
AT(j, i) = A(i, j)
```

输入矩阵按行读取，但输出矩阵会跨步写入。因此转置通常比加法或标量乘法
慢。

#### Frobenius 范数

`MatrixNormFrobenius()` 计算：

```text
sqrt(sum(A(i, j)^2))
```

实现方式是线性遍历矩阵，累计平方和，然后调用标准数学库中的 `sqrt()`。

### 4. 矩阵乘法方法

矩阵乘法是本项目性能优化的重点。对于：

- `A` 的形状为 `m x k`
- `B` 的形状为 `k x n`
- `C` 的形状为 `m x n`

数学定义为：

```text
C(i, j) = sum(A(i, p) * B(p, j)), p = 0 ... k - 1
```

仓库中保留了多种实现，方便对比性能。

#### 4.1 原始 `i-j-k` 方法

`MatrixMultiplyOriginal()` 使用最直接的三重循环：

```text
for i
  for j
    for k
      C(i, j) += A(i, k) * B(k, j)
```

这种方法最容易理解，和数学公式一致。但是对于行优先矩阵，它的缓存表现
较差。`A(i, k)` 是连续读取的，但 `B(k, j)` 每次改变 `k` 时都会跨过一整行。

该方法现在主要作为正确性参考和性能基线。

#### 4.2 循环重排 `i-k-j` 方法

`MatrixMultiplyIKJ()` 改变循环顺序：

```text
for i
  for k
    for j
      C(i, j) += A(i, k) * B(k, j)
```

这种顺序更适合行优先存储：

- `B(k, j)` 按一整行连续读取。
- `C(i, j)` 按一整行连续更新。
- `A(i, k)` 只读取一次，然后用于多个列。

计算量不变，但内存访问更高效。

#### 4.3 分块方法

`MatrixMultiplyBlocked()` 将计算划分为 `64 x 64` 的小块。分块是一种缓存
优化：程序不一次性遍历整个大矩阵，而是反复处理更小的区域，使这些数据
更可能保留在缓存中。

分块版本仍然使用普通标量 C 循环。在块内部，列方向的内层循环手动展开
四个元素：

```c
c_row[j]     += aik * b_row[j];
c_row[j + 1] += aik * b_row[j + 1];
c_row[j + 2] += aik * b_row[j + 2];
c_row[j + 3] += aik * b_row[j + 3];
```

这样可以减少循环开销，也让编译器更容易优化。

#### 4.4 当前优化后的可移植 C 方法

`MatrixMultiply()` 当前调用 `MatrixMultiplyOptimized()`。

这个实现遵守当前限制：不使用 NEON、不使用外部库、不使用平台专用线程 API、
也不使用堆内存打包缓冲区。它只使用标准 C 方法。

优化版本使用：

- `64 x 64` 输出分块
- 标量 `4 x 4` 微内核
- 局部累加变量
- C99 `restrict` 指针
- 针对非方阵和非 4 倍数尺寸的通用尾部处理

微内核每次计算一个 `4 x 4` 的输出块。它不会在最内层循环中反复读写 `C`，
而是把 16 个部分和保存在局部变量中：

```c
REAL c00, c01, c02, c03;
REAL c10, c11, c12, c13;
REAL c20, c21, c22, c23;
REAL c30, c31, c32, c33;
```

对于每一个内部下标 `k`，它读取 `B` 中连续的 4 个元素，并读取 `A` 中 4 行
各 1 个元素，然后更新 16 个累加器，最后只把最终结果写入 `C` 一次。

这种方法更快的原因是：

- 输出值在累加过程中更容易保存在寄存器中。
- 对内存的写入次数减少。
- `B` 的连续元素可以顺序读取。
- 循环体小而规则，便于编译器优化。
- 实现仍然保持可移植 C。

如果矩阵维度不是 4 的倍数，`MultiplyTail()` 会用通用标量循环处理剩余的
行或列。因此优化版本可以正确处理方阵和非方阵。

### 5. LU 分解方法

`LUDecomposeNoPivot()` 将方阵 `A` 分解为：

```text
A = L * U
```

其中：

- `L` 是对角线为 1 的下三角矩阵。
- `U` 是上三角矩阵。

当前实现使用 Doolittle 风格的无主元 LU 分解：

1. 将 `L` 初始化为单位矩阵。
2. 计算 `U` 的当前行。
3. 检查当前主元 `U(k, k)`。
4. 计算 `L` 中主元下方的元素。

如果某个主元绝对值小于给定容差，函数返回 `MATRIX_ERROR_SINGULAR`。

无主元 LU 分解适合教学和输入受控的场景，但数值稳定性不如带部分主元
选取的专业 LAPACK 分解。有些数学上可逆的矩阵，在无主元方法下仍可能
失败或得到不准确的结果。

### 6. 线性方程组求解

完成 LU 分解后，线性方程组：

```text
A * x = b
```

分两步求解：

1. 前向代入求解 `L * y = b`。
2. 回代求解 `U * x = y`。

这种方法避免直接计算 `A` 的逆矩阵，通常更稳定也更高效。

### 7. 行列式方法

无主元 LU 分解后：

```text
det(A) = det(L) * det(U)
```

由于 `L` 的对角线全为 1，所以 `det(L) = 1`。因此：

```text
det(A) = U 对角线元素的乘积
```

这由 `LUDeterminant()` 实现。

### 8. 测试方法

`test_matrix_ops.c` 通过和原始三重循环版本对比，测试优化后的矩阵乘法。

测试用例包括：

- 非常小的矩阵
- 非方阵
- 维度不是 4 的倍数的矩阵
- 较大的不规则尺寸

这很重要，因为优化后的 `4 x 4` 微内核有单独的尾部处理逻辑。测试可以确认
尾部处理路径和简单参考实现一致。

运行测试：

```sh
make run-test
```

### 9. 性能基准测试方法

仓库中有两个 macOS 基准测试程序：

- `benchmark_ops_vs_accelerate.c`
- `benchmark_lu_vs_accelerate.c`

矩阵运算基准测试将本项目实现与 Apple Accelerate vDSP 和 BLAS 对比。
LU 基准测试将本项目无主元 LU 与 Accelerate LAPACK `dgetrf` 对比。

基准测试规则：

- 矩阵初始化不计入计时区域。
- 每个结果报告中位数时间。
- 矩阵乘法额外报告 GFLOP/s。
- 使用校验值避免编译器删除实际计算。

运行基准测试：

```sh
make run-benchmark
```

生成图表：

```sh
make plot-benchmark-preview
```

### 10. 为什么 BLAS 仍然更快

优化后的本项目矩阵乘法已经比原始 C 三重循环快很多，但 Apple Accelerate
BLAS 仍然明显更快。

这是正常结果。专业 BLAS 通常使用：

- 架构专用 SIMD 指令
- 手写汇编或编译器内置指令
- 矩阵分块打包
- 多线程
- 针对具体 CPU 缓存大小调优的内核

当前仓库实现有意避免这些技术，使生产代码保持可移植，并且基于标准 C。
