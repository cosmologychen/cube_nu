#!/bin/bash
#SBATCH -p aarch64
#SBATCH -N 1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=16
#SBATCH -J cube_test
#SBATCH -o test_slurm_%j.out
#SBATCH -e test_slurm_%j.err

set -eo pipefail

# 1. 环境准备
if [ -f "/share/users/twtest/miniconda3/etc/profile.d/conda.sh" ]; then
    source "/share/users/twtest/miniconda3/etc/profile.d/conda.sh"
    conda activate base
fi
source ~/.bashrc

cd /share/users/twtest/cube_nu
source env.sh

code_path=$(pwd)
export OMP_STACKSIZE=2000M
export OMP_NUM_THREADS=16
ulimit -s unlimited

echo "=========================================================="
echo " 阶段 1: 开始串行测试 (Serial Test)"
echo "=========================================================="
export OMP_NUM_THREADS=1
bash test/serial/run.sh
echo "串行测试完成！"

echo "=========================================================="
echo " 阶段 2: 开始并行测试 (Parallel Test - mnu000, 8 Images)"
echo "=========================================================="
export OMP_NUM_THREADS=16
PROCS=8
work_path="${code_path}/../test/code/mnu000"
result_path="${code_path}/../test/mnu000"

rm -rf "${work_path}" "${result_path}"
mkdir -p "${work_path}" "${result_path}/neutrinos"

# 复制代码与测试输入
rsync -a --exclude='test' "${code_path}/" "${work_path}/"
cp -vr test/baseline/mnu000/neutrinos/* "${result_path}/neutrinos/"
cp -v test/baseline/mnu000/parameters.f90 "${work_path}/parameters.f90"

# 复制 noise 输入文件
cp -r "${code_path}/test/baseline/noise/"* "${result_path}/"

# 编译主程序与工具
cd "${work_path}"
sed -i "s#.*sim%cur_checkpoint=.* !set checkpoint number#   sim%cur_checkpoint=1 !set checkpoint number#" "main.f90"
make clean && make
cd utilities && make clean && make && cd ..

# 运行并行作业
cafrun -np ${PROCS} ./utilities/ic.x
cafrun -np ${PROCS} ./main.x
cafrun -np ${PROCS} ./utilities/cicpower.x

# 测试断点续传
sed -i "s#.*sim%cur_checkpoint=.* !set checkpoint number#   sim%cur_checkpoint=5 !set checkpoint number#" "main.f90"
make
cafrun -np ${PROCS} ./main.x
echo "并行模拟完成！"

echo "=========================================================="
echo " 阶段 3: 执行并行结果对比 (compare.py)"
echo "=========================================================="
cd "${code_path}"
python test/compare.py

echo "=========================================================="
echo " 所有测试运行完成！"
echo "=========================================================="
