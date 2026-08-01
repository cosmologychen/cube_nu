#!/bin/bash
set -euo pipefail

code_path="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
source_path="${code_path}/../test/serial"
baseline_path="${code_path}/test/baseline/serial"

if [[ ! -d "${source_path}" ]]; then
    echo "错误：找不到串行模拟结果目录：${source_path}" >&2
    exit 1
fi

required_files=(
    "image1/0.000_delta_c_1.bin"
    "image1/0.000_proj_xy_1.bin"
    "image1/0.000_cicpower_1.bin"
)

for file in "${required_files[@]}"; do
    if [[ ! -f "${source_path}/${file}" ]]; then
        echo "错误：缺少基准文件：${source_path}/${file}" >&2
        exit 1
    fi
done

rm -rf "${baseline_path}"
mkdir -p "${baseline_path}"

cp -va "${source_path}/image1/0.000_delta_c_1.bin" "${baseline_path}/"
cp -va "${source_path}/image1/0.000_proj_xy_1.bin" "${baseline_path}/"
cp -va "${source_path}/image1/0.000_cicpower_1.bin" "${baseline_path}/"

if [[ -f "${source_path}/run.log" ]]; then
    cp -a "./run.log" "${baseline_path}/"
fi

echo "串行基准已固化到：${baseline_path}"
