#!/bin/bash
set -euo pipefail

# 仓库根目录和串行测试输出目录
code_path="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
result_path="${code_path}/../test/serial"
work_path="${result_path}/../code/serial"

# 串行运行参数
export OMP_STACKSIZE=2000M
export OMP_NUM_THREADS=1
export FOR_COARRAY_NUM_IMAGES=1

ulimit -s unlimited
ulimit -v unlimited

rm -rf "${result_path}"
mkdir -p "${work_path}"
mkdir -p "${result_path}/neutrinos" "${result_path}/image1"

# 复制源码，排除测试目录，避免递归复制测试输出
rsync -a --exclude='test' --exclude='test/baseline/serial' "${code_path}/" "${work_path}/"
cp -v $code_path/test/serial/parameters.f90 $work_path/parameters.f90

# 串行测试只使用 baseline/noise/image1 的输入
cp -v $code_path/test/baseline/noise/image1/noise_1.bin $result_path/image1/noise_1.bin
cp -v $code_path/test/baseline/noise/image1/seed_1.bin $result_path/image1/seed_1.bin
cp -vr $code_path/test/baseline/mnu000/neutrinos/* $result_path/neutrinos/

cd $work_path
sed -i "s#.*sim%cur_checkpoint=.* !set checkpoint number#   sim%cur_checkpoint=1 !set checkpoint number#" "main.f90"

# 编译主程序及工具。具体编译器和库环境由调用环境提供。
make clean
make
cd utilities
make clean
make
cd ..

IC=./utilities/ic.x
MAIN=./main.x
CIC=./utilities/cicpower.x

# 运行初始条件、主程序和 CIC 功率谱计算。
"${IC}"
"${MAIN}"
"${CIC}"


# 设置新的检查点
cd $work_path
sed -i "s#.*sim%cur_checkpoint=.* !set checkpoint number#   sim%cur_checkpoint=5 !set checkpoint number#" "main.f90"

# 运行主程序-测试断点续传
make
$MAIN

echo "串行 serial 测试完成：${result_path}/serial"
