#!/bin/bash
#SBATCH -p aarch64
#SBATCH -N 2
#SBATCH --ntasks=8
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=32
#SBATCH -J cube_parallel
#SBATCH -o ./log/test_parallel_%j.out
#SBATCH -e ./log/test_parallel_%j.err

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
export PATH=/share/softwares/aarch64/libs/openmpi/4.1.8-ucx/bin:$PATH
export LD_LIBRARY_PATH=/share/softwares/aarch64/libs/openmpi/4.1.8-ucx/lib:$LD_LIBRARY_PATH
export UCX_NET_DEVICES="bond0,rocep150s0f0:1"
export UCX_SOCKADDR_TLS_PRIORITY=tcp
export OMPI_MCA_pml=ucx
export OMPI_MCA_osc=pt2pt
export LD_PRELOAD=/share/users/twtest/lib/libforcefunneled.so
export CAF_FORCE_THREAD_LEVEL=funneled

# 覆盖 cafrun（OpenCoarrays 的 cafrun 硬编码了旧 mpiexec 路径，这里改用 openmpi-ucx 的 mpiexec）
cafrun() {
    /share/softwares/aarch64/libs/openmpi/4.1.8-ucx/bin/mpiexec "$@"
}
export -f cafrun

code_path=$(pwd)
export OMP_STACKSIZE=2000M
export OMP_NUM_THREADS=8
ulimit -s unlimited

echo "=========================================================="
echo " 开始并行测试 (Parallel Test - mnu000, 8 Images, 2 Nodes, RoCE)"
echo " 分配节点: $(scontrol show hostnames $SLURM_JOB_NODELIST | tr '\n' ' ')"
echo "=========================================================="
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
echo " 执行并行结果对比 (compare.py)"
echo "=========================================================="
cd "${code_path}"
python test/compare.py

echo "=========================================================="
echo " 并行测试全部完成！"
echo "=========================================================="
