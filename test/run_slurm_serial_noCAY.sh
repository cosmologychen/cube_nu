#!/bin/bash
#SBATCH -p aarch64
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH -J serial_noCAY
#SBATCH -o /share/users/twtest/cube_nu/log/test_serial_noCAY_%j.out
#SBATCH -e /share/users/twtest/cube_nu/log/test_serial_noCAY_%j.err

set -eo pipefail
source ~/.bashrc
source /share/users/twtest/cube_nu/env.sh
unset LD_PRELOAD CAF_FORCE_THREAD_LEVEL

code_path=/share/users/twtest/cube_nu
test_path="${code_path}/.test_ai/serial_noCAY"
work_path="${test_path}/code/serial_noCAY"
result_path="${test_path}/serial_noCAY"

export FC=gfortran
export OFLAG="-O3 -cpp -fopenmp -mcmodel=large -ffree-line-length-none -fcoarray=single -c"
export XFLAG="-O3 -fopenmp -mcmodel=large -fcoarray=single"
export OMP_STACKSIZE=2000M
export OMP_NUM_THREADS=1
ulimit -s unlimited

mkdir -p "${work_path}" "${result_path}/neutrinos" "${result_path}/image1"
rsync -a --exclude='test' --exclude='log/' --exclude='.*/' "${code_path}/" "${work_path}/"
cp "${code_path}/test/serial/parameters.f90" "${work_path}/parameters.f90"
sed -i "s#opath='../../serial/'#opath='../../serial_noCAY/'#" "${work_path}/parameters.f90"
cp "${code_path}/test/baseline/noise/image1/"{noise_1,seed_1}.bin "${result_path}/image1/"
cp -r "${code_path}/test/baseline/mnu000/neutrinos/". "${result_path}/neutrinos/"

cd "${work_path}"
sed -i "s#.*sim%cur_checkpoint=.* !set checkpoint number#   sim%cur_checkpoint=1 !set checkpoint number#" main.f90
make clean
make
make -C utilities clean
make -C utilities

echo "检查可执行文件依赖：应无 MPI/OpenCoarrays"
for exe in ./utilities/ic.x ./main.x ./utilities/cicpower.x; do
    dependencies=$(ldd "${exe}")
    printf '%s\n%s\n' "${exe}" "${dependencies}"
    if [[ "${dependencies}" == *libmpi* || "${dependencies}" == *libcaf* ]]; then
        echo "错误：仍存在 MPI/OpenCoarrays 动态依赖"
        exit 1
    fi
done

./utilities/ic.x
./main.x
./utilities/cicpower.x

sed -i "s#.*sim%cur_checkpoint=.* !set checkpoint number#   sim%cur_checkpoint=5 !set checkpoint number#" main.f90
make
./main.x
echo "serial_noCAY 模拟与断点续传完成：${result_path}"
