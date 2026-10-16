#!/bin/bash
#SBATCH -p aarch64
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH -J cube_serial
#SBATCH -o ./log/test_serial_%j.out
#SBATCH -e ./log/test_serial_%j.err

set -eo pipefail

# 1. 环境准备
source ~/.bashrc
if [ -f "/share/users/twtest/miniconda3/etc/profile.d/conda.sh" ]; then
    source "/share/users/twtest/miniconda3/etc/profile.d/conda.sh"
    conda activate camb
fi

cd /share/users/twtest/cube_nu
source env.sh

export LD_PRELOAD=/share/users/twtest/lib/libforcefunneled.so
export CAF_FORCE_THREAD_LEVEL=funneled
export OMPI_MCA_btl_tcp_if_include="172.16.12.0/24"

code_path=$(pwd)
export OMP_STACKSIZE=2000M
export OMP_NUM_THREADS=1
ulimit -s unlimited

echo "=========================================================="
echo " 开始串行测试 (Serial Test)"
echo "=========================================================="
bash test/serial/run.sh
echo "串行测试完成！"
