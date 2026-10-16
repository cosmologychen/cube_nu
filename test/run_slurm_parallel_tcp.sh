#!/bin/bash    废弃脚本
#SBATCH -p aarch64
#SBATCH -N 2
#SBATCH --ntasks=8
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=32
#SBATCH -J cube_parallel_TCP
#SBATCH -o ./log/test_parallel_TCP_%j.out
#SBATCH -e ./log/test_parallel_TCP_%j.err

set -eo pipefail

# 1. 环境准备
if [ -f "/share/users/twtest/miniconda3/etc/profile.d/conda.sh" ]; then
    source "/share/users/twtest/miniconda3/etc/profile.d/conda.sh"
    conda activate base
fi
source ~/.bashrc

cd /share/users/twtest/cube_nu
source env.sh

# ---------------------------------------------------------------------------
export LD_PRELOAD=/share/users/twtest/lib/libforcefunneled.so
export CAF_FORCE_THREAD_LEVEL=funneled

# 约束 OpenMPI 跨节点走 172.16.12.0/24 (bond0, 50G)
export OMPI_MCA_pml=ob1
export OMPI_MCA_btl=self,vader,tcp
export OMPI_MCA_btl_tcp_if_include="172.16.12.0/24"

code_path=$(pwd)
export OMP_STACKSIZE=2000M
export OMP_NUM_THREADS=8
ulimit -s unlimited

echo "=========================================================="
echo " 开始并行测试 (Parallel Test - mu000_TCP, 8 Images, 2 Nodes, TCP)"
echo " 分配节点: $(scontrol show hostnames $SLURM_JOB_NODELIST | tr '\n' ' ')"
echo "=========================================================="
PROCS=8
work_path="${code_path}/../test/code/mu000_TCP"
result_path="${code_path}/../test/mu000_TCP"

rm -rf "${work_path}" "${result_path}"
mkdir -p "${work_path}" "${result_path}/neutrinos"

# 复制代码与测试输入
rsync -a --exclude='test' --exclude='log/' --exclude='.*/' "${code_path}/" "${work_path}/"
cp -vr test/baseline/mnu000/neutrinos/* "${result_path}/neutrinos/"
cp -v test/baseline/mnu000/parameters.f90 "${work_path}/parameters.f90"
sed -i "s#.*character(\*),parameter :: opath=.*#   character(*),parameter :: opath='../../mu000_TCP/'#" "${work_path}/parameters.f90"

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
echo " 并行测试全部完成！"
echo "=========================================================="
