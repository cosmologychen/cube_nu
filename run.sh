
exe=int2real
exe=fof
echo "编译${exe}工具..."
cd utilities/
if ! make ${exe}.x; then
    echo "错误：${exe}编译失败！"
    exit 1
fi
cd ..
# 检查编译结果
if [ ! -f "utilities/${exe}.x" ]; then
    echo "错误：${exe}可执行文件未生成！"
    exit 1
fi
echo "${exe}编译成功，提交作业..."

qsub run.qsub
# tail -f run.log