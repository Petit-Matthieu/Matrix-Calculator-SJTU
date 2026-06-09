# 矩阵乘法优化原理详解

本项目实现了 4 种矩阵乘法，按优化程度递进。以下逐一介绍每种方法的算法原理、性能瓶颈和优化思路。

---

## 1. 原始三重循环（MatrixMultiplyOriginal）

### 算法

矩阵乘法的数学定义：

```
C[i][j] = Σ A[i][k] × B[k][j]   (k = 0 .. n-1)
```

直接翻译为三层嵌套循环：

```c
for (int i = 0; i < A->row; ++i) {
    for (int j = 0; j < B->column; ++j) {
        REAL sum = 0.0;
        for (int k = 0; k < A->column; ++k) {
            sum += A->data[MatrixIndex(A, i, k)] * B->data[MatrixIndex(B, k, j)];
        }
        C->data[MatrixIndex(C, i, j)] = sum;
    }
}
```

循环顺序为 **i → j → k**。

### 存储方式

本项目使用**行优先（row-major）**存储，即矩阵的同一行元素在内存中连续排列。对于 `Matrix[A][B]`，元素 `(i, j)` 的地址为：

```
data[i * column + j]
```

### 性能瓶颈

内层循环沿 `k` 方向迭代时：

- `A[i][k]`：沿行方向移动，地址连续，**空间局部性好**
- `B[k][j]`：沿列方向移动，每次跳过一整行（例如 768 个 double = 6144 字节），**空间局部性极差**

现代 CPU 的缓存行（cache line）通常为 64 字节（8 个 double）。当访问 `B[k][j]` 时，CPU 会把包含该元素的整条缓存行加载到 L1 缓存，但内层循环下一步访问的是 `B[k+1][j]`——它距离 `B[k][j]` 恰好一整行，几乎不可能在同一条缓存行内。

结果是：**每一次内层迭代都触发一次 L1 缓存未命中**，数据被迫从更远的 L2、L3 甚至主存加载，延迟从 ~4 个时钟周期暴涨到 ~100-300 个周期。

### 复杂度

- 时间复杂度：O(n³)
- 空间复杂度：O(1)（不计输出矩阵）
- 实测 768×768：**564.9ms**

---

## 2. 循环重排 i-k-j（MatrixMultiplyIKJ）

### 核心思想

将循环顺序从 `i-j-k` 改为 `i-k-j`，让内层循环沿 `j` 方向（行方向）遍历 B 矩阵。

```c
for (int i = 0; i < A->row; ++i) {
    for (int k = 0; k < A->column; ++k) {
        REAL aik = A->data[MatrixIndex(A, i, k)];  // 提取为标量
        for (int j = 0; j < B->column; ++j) {
            C->data[MatrixIndex(C, i, j)] += aik * B->data[MatrixIndex(B, k, j)];
        }
    }
}
```

### 为什么更快

**1. B 矩阵访问模式改善**

| 版本 | 内层循环方向 | B 的访问模式 | 缓存行为 |
|------|-------------|-------------|---------|
| i-j-k | k | `B[0][j], B[1][j], B[2][j]...`（列方向跳跃） | 每次 cache miss |
| i-k-j | j | `B[k][0], B[k][1], B[k][2]...`（行方向连续） | 几乎全部 cache hit |

内层循环沿 j 方向移动时，`B[k][j]` 和 `B[k][j+1]` 在内存中相邻（行优先存储），地址相差仅 8 字节，落在同一条 64 字节缓存行内。一个缓存行可以服务 8 次连续访问。

**2. A 矩阵的标量复用**

`aik` 被提取为局部变量后，编译器可以将它保留在寄存器中。内层 j 循环的每次迭代只需要一次乘法和一次加法访问 C 矩阵，而不需要每次都索引 A 矩阵。

**3. C 矩阵的访问模式改善**

`C[i][j]` 在内层 j 循环中也是按行连续访问的，写回同样具备空间局部性。

### 复杂度

- 时间复杂度：O(n³)（不变）
- 实测 768×768：**~215ms**（约 2.6× 加速）

### 局限

当矩阵尺寸大到 C 矩阵本身也超出缓存容量时，`C[i][j]` 的反复读-改-写仍然会导致缓存颠簸（cache thrashing）。需要进一步引入分块策略。

---

## 3. 分块乘法（MatrixMultiplyBlocked）

### 核心思想

将大矩阵切分为若干个小块（tile），保证每个小块的工作数据能完全放入 L1 缓存，在块内完成计算后再移动到下一块。

```c
#define BLOCK_SIZE 64

for (int ii = 0; ii < A->row; ii += BLOCK_SIZE) {
    for (int kk = 0; kk < A->column; kk += BLOCK_SIZE) {
        for (int jj = 0; jj < B->column; jj += BLOCK_SIZE) {
            // 在 64×64 的小块内做乘法
            for (int i = ii; i < i_end; ++i) {
                for (int k = kk; k < k_end; ++k) {
                    REAL aik = a_row[k];
                    for (int j = jj; j + 3 < j_end; j += 4) {  // 4 路展开
                        c_row[j]   += aik * b_row[j];
                        c_row[j+1] += aik * b_row[j+1];
                        c_row[j+2] += aik * b_row[j+2];
                        c_row[j+3] += aik * b_row[j+3];
                    }
                }
            }
        }
    }
}
```

### 分块原理

假设使用 64×64 的分块：

- 每个 A 的行块：64 行 × 64 列 × 8 字节 = 32 KB
- 每个 B 的列块：64 行 × 64 列 × 8 字节 = 32 KB
- C 的输出块：64 × 64 × 8 字节 = 32 KB

三个块合计约 96 KB，而典型的 L1 数据缓存为 32-48 KB。通过分块，A 和 B 的工作集可以分时复用 L1，C 的输出块也能常驻。

### 4 路循环展开

```c
for (; j + 3 < j_end; j += 4) {
    c_row[j]   += aik * b_row[j];
    c_row[j+1] += aik * b_row[j+1];
    c_row[j+2] += aik * b_row[j+2];
    c_row[j+3] += aik * b_row[j+3];
}
```

循环展开（loop unrolling）的作用：

1. **减少循环开销**：每次迭代处理 4 个元素，循环分支减少 75%
2. **指令级并行**：4 次乘加相互独立，CPU 的乱序执行引擎可以同时发射
3. **减少指针运算**：`j` 的自增和比较次数减少 75%

尾部（维度不是 4 的倍数）由 `j < j_end` 的收尾循环处理，保证正确性。

### 复杂度

- 时间复杂度：O(n³)（不变）
- 实测 768×768：**~53.5ms**（约 10.5× 加速）

### 局限

分块版本虽然解决了缓存容量问题，但内层的 `j` 循环仍然对 C 矩阵进行读-改-写。当分块内的 C 子矩阵在计算过程中被反复访问时，每次迭代都要经历"读取 C[i][j] → 加上乘积 → 写回 C[i][j]"的完整缓存访问周期。

---

## 4. 4×4 微内核优化（MatrixMultiplyOptimized）

这是本项目的最终优化版本，也是 `MatrixMultiply()` 的默认实现。

### 整体结构

采用 **三层循环 + 微内核** 的架构：

```c
// 外层：按 64×64 分块（与分块版本相同）
for (jj = 0; jj < B->column; jj += BLOCK_SIZE) {
    for (ii = 0; ii < A->row; ii += BLOCK_SIZE) {
        // 中层：每次取 A 的 4 行
        for (i = ii; i + 3 < i_end; i += 4) {
            // 每次取 B 的 4 列
            for (j = jj; j + 3 < j_end; j += 4) {
                Multiply4x4(A + i*A->col, A->col,
                            B + j,       B->col,
                            C + i*C->col + j, C->col,
                            A->col);
            }
            // 尾部处理
            if (j < j_end) MultiplyTail(...);
        }
        if (i < i_end) MultiplyTail(...);
    }
}
```

### Multiply4x4 微内核详解

这是整个优化的核心。它计算 C 矩阵的一个 4×4 子块：

```c
static void Multiply4x4(const REAL * restrict a, int a_stride,
                        const REAL * restrict b, int b_stride,
                        REAL * restrict c, int c_stride, int inner)
{
    // ① 16 个局部累加器——映射到 CPU 寄存器
    REAL c00=0, c01=0, c02=0, c03=0;
    REAL c10=0, c11=0, c12=0, c13=0;
    REAL c20=0, c21=0, c22=0, c23=0;
    REAL c30=0, c31=0, c32=0, c33=0;

    // ② A 的 4 行指针
    const REAL *a0 = a;
    const REAL *a1 = a + a_stride;
    const REAL *a2 = a + 2 * a_stride;
    const REAL *a3 = a + 3 * a_stride;

    // ③ 沿内维 k 累加
    for (int k = 0; k < inner; ++k) {
        // B 的一行 4 个元素，被 A 的 4 行共享
        const REAL *b_row = b + k * b_stride;
        REAL b0 = b_row[0], b1 = b_row[1],
             b2 = b_row[2], b3 = b_row[3];

        // A 的 4 行各取 1 个元素
        REAL a0k = a0[k], a1k = a1[k],
             a2k = a2[k], a3k = a3[k];

        // ④ 4×4 = 16 次乘加（FMA），全部在寄存器中完成
        c00 += a0k*b0; c01 += a0k*b1; c02 += a0k*b2; c03 += a0k*b3;
        c10 += a1k*b0; c11 += a1k*b1; c12 += a1k*b2; c13 += a1k*b3;
        c20 += a2k*b0; c21 += a2k*b1; c22 += a2k*b2; c23 += a2k*b3;
        c30 += a3k*b0; c31 += a3k*b1; c32 += a3k*b2; c33 += a3k*b3;
    }

    // ⑤ 循环结束后一次性写回 C 矩阵
    c[0]=c00; c[1]=c01; c[2]=c02; c[3]=c03; c += c_stride;
    c[0]=c10; c[1]=c11; c[2]=c12; c[3]=c13; c += c_stride;
    c[0]=c20; c[1]=c21; c[2]=c22; c[3]=c23; c += c_stride;
    c[0]=c30; c[1]=c31; c[2]=c32; c[3]=c33;
}
```

### 5 个关键优化点

#### ① 16 个局部累加器 → 零内存访问

16 个 `c**` 变量会被编译器分配到 16 个物理寄存器（现代 x86-64 CPU 有 16 个通用寄存器，ARM64 有 31 个）。内层 `k` 循环的 16 次乘加操作完全在寄存器中完成，**不需要在每次迭代中读取和写回 C 矩阵**。

对比分块版本：

| 版本 | 内层循环中 C 的访问 |
|------|-------------------|
| 分块版本 | 每次迭代：读 C[i][j] → 加 → 写回 C[i][j]（2 次内存操作） |
| 微内核版本 | 整个 k 循环结束后：写回 16 个 C 元素（仅 16 次内存写） |

对于 768×768 的矩阵，内维 k 循环执行 768 次，微内核将 768×16×2 = 24576 次内存操作减少到 16 次。

#### ② B 行 4 列复用 → 数据复用率 ×4

在 k 循环的每次迭代中，`b_row[0..3]` 被加载到 4 个寄存器变量 `b0-b3` 中。这 4 个值分别与 `a0k, a1k, a2k, a3k` 相乘——也就是说，B 矩阵的一行数据被 A 的 4 行**共享复用**。

数据复用分析（每次 k 迭代）：
- B 矩阵：读取 4 个元素（`b0-b3`）
- A 矩阵：读取 4 个元素（`a0k-a3k`）
- 产出：16 个乘加结果

内存读取效率 = 16 次有效计算 / 8 次内存读取 = **2 FLOP/byte**，是朴素循环的 2 倍。

#### ③ C99 restrict 指针 → 消除别名假设

```c
const REAL * restrict a, ...
const REAL * restrict b, ...
REAL * restrict c, ...
```

C 语言编译器默认假设所有指针可能指向同一块内存（别名问题）。这意味着写入 `c[0]` 可能影响 `a[1]` 的值，编译器被迫在每次写入后重新从内存读取 `a` 和 `b` 的值。

`restrict` 关键字向编译器保证：**这三块内存区域互不重叠**。编译器因此可以：
- 将 `a0k, b0` 等变量保留在寄存器中，跨越 k 循环迭代而不需要重载
- 更激进地调度指令流水线
- 在某些架构上启用 SIMD 向量化

#### ④ 尾部处理 → 支持任意尺寸

```c
static void MultiplyTail(const Matrix *A, const Matrix *B, Matrix *C,
                         int i_begin, int i_end, int j_begin, int j_end)
{
    for (int i = i_begin; i < i_end; ++i) {
        for (int j = j_begin; j < j_end; ++j) {
            REAL sum = 0.0;
            for (int k = 0; k < A->column; ++k)
                sum += a_row[k] * B->data[k * B->column + j];
            c_row[j] = sum;
        }
    }
}
```

当矩阵维度不是 4 的倍数（例如 100×100）时，4×4 微内核无法覆盖边缘的行列。`MultiplyTail` 用朴素循环处理这些"尾部"，保证任意尺寸矩阵的计算正确性。

#### ⑤ 分块外层 → 工作集常驻 L1

外层的 64×64 分块保证了微内核处理的数据子集适合 L1 缓存：

```
A 的 4 行 × 64 列 = 4 × 64 × 8 = 2048 字节
B 的 64 行 × 4 列 = 64 × 4 × 8 = 2048 字节
C 的 4 行 × 4 列  = 4 × 4 × 8   = 128 字节
```

每次微内核调用的工作集仅 ~4KB，远小于 L1 缓存（32-48KB），保证了数据命中率。

### 复杂度

- 时间复杂度：O(n³)（不变）
- 实测 768×768：**30.3ms**
- 吞吐量：**29.924 GFLOP/s**
- 对比原始版本：**18.66× 加速**

---

## 总结对比

| 版本 | 关键技术 | 768×768 耗时 | 加速比 | 核心改善 |
|------|---------|-------------|--------|---------|
| 原始 i-j-k | 三重循环 | 564.9ms | 1× | 基准 |
| i-k-j 重排 | 循环重排 | ~215ms | ~2.6× | B 按行访问，空间局部性 |
| 分块 | 64×64 分块 + 4 路展开 | ~53.5ms | ~10.5× | 工作集适配 L1 缓存 |
| 4×4 微内核 | 寄存器累加 + restrict | **30.3ms** | **18.66×** | 零内存读写 + 数据复用 |

四种实现全部保留在代码中，可通过函数名直接调用对比。`MatrixMultiply()` 默认调用优化版本。
