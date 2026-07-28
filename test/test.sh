#!/bin/bash
# 确定运行目录
code_path=/home/cbh/18T/cube_nus/cube_nu
cd $code_path
tar -zxvf test/baseline.tar.gz

sed -i "s#.*code_path=.*#code_path=${code_path}/#" "test/baseline/mnu000/run.qsub"
sed -i "s#.*code_path=.*#code_path=${code_path}/#" "test/baseline/mnu006/run.qsub"
sed -i "s#.*code_path=.*#code_path='${code_path}/'#" "test/compare.py"
# # 运行模拟参数1
# qsub test/baseline/mnu000/run.qsub
# # 运行模拟参数2
# qsub test/baseline/mnu006/run.qsub

# #运行完成后用test/compare.py对比运行结果