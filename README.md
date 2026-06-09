# Matrix Calculator SJTU

矩阵工具箱 — C99 矩阵运算库 + 浏览器交互式计算器 + 乘法性能对比

## 浏览器工具箱

直接打开 [`matrix-toolbox.html`](matrix-toolbox.html) 即可使用，无需服务器或安装依赖。

### 📊 矩阵计算

支持 8 种运算，自动同步维度，方阵运算自动锁定行列：

| 运算 | 说明 |
|------|------|
| A + B / A − B | 加法、减法，AB 维度自动同步 |
| A × B | 乘法，A 列自动匹配 B 行 |
| k · A | 标量乘法 |
| Aᵀ | 转置 |
| det(A) | 行列式（方阵） |
| LU 分解 | 输出 L、U 矩阵和分解步骤 |
| Ax = b | 线性方程组求解，支持快捷输入向量 b |

附加功能：粘贴导入矩阵（支持 Excel 复制）、单位阵/随机/清空一键填充。

### ⚡ 乘法性能对比

在浏览器中对比 4 种矩阵乘法实现的耗时：

| 版本 | 关键技术 |
|------|----------|
| IJK | i-j-k 三重循环（基准） |
| IKJ | 循环重排 i-k-j |
| Blocked | 64×64 分块 + 4 路循环展开 |
| Optimized | 4×4 微内核 + 分块 |

支持 100 ~ 1500 规模选择，含进度条、内存估算和随机矩阵预览。

## C 语言矩阵库

### 头文件

| 头文件 | 用途 |
|--------|------|
| [`matrix_core.h`](matrix_core.h) | 矩阵生命周期、存储、索引、复制、填充、打印和错误处理 |
| [`matrix_ops.h`](matrix_ops.h) | 算术运算、转置、范数和多种矩阵乘法实现 |
| [`matrix_lu.h`](matrix_lu.h) | LU 分解、代入求解和行列式计算 |

### 快速开始

```sh
make            # 编译
make run        # 运行 LU 演示
make run-test   # 运行乘法回归测试
make clean      # 清理
```

### 乘法优化策略

`MatrixMultiply()` 默认调用优化的可移植 C 实现：

- 适配缓存大小的 64×64 输出分块
- 标量 4×4 微内核，局部累加变量减少写回
- C99 `restrict` 指针提供别名信息
- 通用尾部处理，支持非方阵和非 4 倍数维度

768×768 实测：优化版比原始 ijk 快 **18.66×**，比基础分块快 **1.77×**。

### 仓库结构

| 路径 | 说明 |
|------|------|
| [`matrix_core.c`](matrix_core.c) | 矩阵存储和生命周期管理 |
| [`matrix_ops.c`](matrix_ops.c) | 算术运算和可移植矩阵乘法内核 |
| [`matrix_lu.c`](matrix_lu.c) | 无主元 LU 分解及相关运算 |
| [`main.c`](main.c) | LU 分解演示程序 |
| [`student_practice.c`](student_practice.c) | 学生练习程序 |
| [`test_matrix_ops.c`](test_matrix_ops.c) | 矩阵乘法回归测试 |
| [`benchmark_ops_vs_accelerate.c`](benchmark_ops_vs_accelerate.c) | 与 Apple Accelerate 对比的矩阵运算基准测试 |
| [`benchmark_lu_vs_accelerate.c`](benchmark_lu_vs_accelerate.c) | 与 LAPACK 对比的 LU 基准测试 |
| [`matrix-toolbox.html`](matrix-toolbox.html) | 浏览器交互式矩阵工具箱 |

### 数值说明

`LUDecomposeNoPivot()` 实现为不选主元的 LU 分解，适用于教学和输入受控场景。
