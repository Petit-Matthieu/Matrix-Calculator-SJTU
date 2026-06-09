# 指定 C 编译器为 gcc
CC = gcc

# 默认构建使用优化版本；需要逐步调试时可运行 make debug
COMMON_CFLAGS = -std=c99 -Wall -Wextra -pedantic
RELEASE_CFLAGS = $(COMMON_CFLAGS) -O3 -DNDEBUG
DEBUG_CFLAGS = $(COMMON_CFLAGS) -O0 -g
CFLAGS = $(RELEASE_CFLAGS)

# 链接数学库，因为 LU 模块中使用 fabs
LDFLAGS = -lm

# 教师演示程序和学生练习程序名称
TARGET = lab6_demo
PRACTICE = lab6_practice
TEST = test_matrix_ops
LU_BENCHMARK = benchmark_lu_vs_accelerate
OPS_BENCHMARK = benchmark_ops_vs_accelerate
BENCHMARKS = $(OPS_BENCHMARK) $(LU_BENCHMARK)
BENCHMARK_CHART = benchmark_comparison.svg
BENCHMARK_CHART_PREVIEW = benchmark_comparison.png

# 公共目标文件：矩阵核心模块 + 基本运算模块 + LU 模块
OBJS = matrix_core.o matrix_ops.o matrix_lu.o

all: $(TARGET) $(PRACTICE)

debug:
	$(MAKE) clean
	$(MAKE) CFLAGS="$(DEBUG_CFLAGS)" all

benchmark: $(BENCHMARKS)

test: $(TEST)

$(TARGET): main.o $(OBJS)
	$(CC) $(CFLAGS) -o $(TARGET) main.o $(OBJS) $(LDFLAGS)

$(PRACTICE): student_practice.o $(OBJS)
	$(CC) $(CFLAGS) -o $(PRACTICE) student_practice.o $(OBJS) $(LDFLAGS)

$(TEST): test_matrix_ops.c matrix_core.c matrix_core.h matrix_ops.c matrix_ops.h
	$(CC) $(RELEASE_CFLAGS) \
		-o $(TEST) test_matrix_ops.c matrix_core.c matrix_ops.c -lm

main.o: main.c matrix_core.h matrix_lu.h
	$(CC) $(CFLAGS) -c main.c

student_practice.o: student_practice.c matrix_core.h matrix_lu.h
	$(CC) $(CFLAGS) -c student_practice.c

matrix_core.o: matrix_core.c matrix_core.h
	$(CC) $(CFLAGS) -c matrix_core.c

matrix_ops.o: matrix_ops.c matrix_ops.h matrix_core.h
	$(CC) $(CFLAGS) -c matrix_ops.c

matrix_lu.o: matrix_lu.c matrix_lu.h matrix_core.h
	$(CC) $(CFLAGS) -c matrix_lu.c

$(LU_BENCHMARK): benchmark_lu_vs_accelerate.c matrix_core.c matrix_core.h matrix_lu.c matrix_lu.h
	$(CC) $(RELEASE_CFLAGS) -DACCELERATE_NEW_LAPACK \
		-o $(LU_BENCHMARK) benchmark_lu_vs_accelerate.c matrix_core.c matrix_lu.c \
		-framework Accelerate -lm

$(OPS_BENCHMARK): benchmark_ops_vs_accelerate.c matrix_core.c matrix_core.h matrix_ops.c matrix_ops.h
	$(CC) $(RELEASE_CFLAGS) -DACCELERATE_NEW_LAPACK \
		-o $(OPS_BENCHMARK) benchmark_ops_vs_accelerate.c matrix_core.c matrix_ops.c \
		-framework Accelerate -lm

run: $(TARGET)
	./$(TARGET)

run-practice: $(PRACTICE)
	./$(PRACTICE)

run-test: $(TEST)
	./$(TEST)

run-benchmark-ops: $(OPS_BENCHMARK)
	./$(OPS_BENCHMARK)

run-benchmark-lu: $(LU_BENCHMARK)
	./$(LU_BENCHMARK)

run-benchmark: $(BENCHMARKS)
	./$(OPS_BENCHMARK)
	./$(LU_BENCHMARK)

plot-benchmark: $(BENCHMARKS) plot_benchmark_comparison.py
	python3 plot_benchmark_comparison.py

plot-benchmark-preview: plot-benchmark
	qlmanage -t -s 1800 -o /tmp $(BENCHMARK_CHART)
	cp /tmp/$(BENCHMARK_CHART).png $(BENCHMARK_CHART_PREVIEW)

clean:
	rm -f $(TARGET) $(PRACTICE) $(TEST) $(BENCHMARKS) $(BENCHMARK_CHART) $(BENCHMARK_CHART_PREVIEW) *.o
