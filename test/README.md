# CUBE 测试目录说明

## 1. 当前测试流程

本目录提供 CUBE 的测试脚本、基准结果压缩包和 Python 结果对比脚本。运行测试脚本时，会自动解压基准结果：

- `mnu000`：无中微子质量，对应 `m_nu = [0.0, 0.0, 0.0]`。
- `mnu006`：一个有质量中微子，对应 `m_nu = [0.06, 0.0, 0.0]`。

当前实际流程由 [`test/test.sh`](test.sh) 驱动。使用前请先将该脚本中的 `code_path` 修改为当前仓库根目录：

1. 解压 [`test/baseline.tar.gz`](baseline.tar.gz)，生成 `test/baseline/`。
2. 在仓库根目录提交 [`test/baseline/mnu000/run.qsub`](baseline/mnu000/run.qsub)。
3. 在仓库根目录提交 [`test/baseline/mnu006/run.qsub`](baseline/mnu006/run.qsub)。
4. 等待两个 PBS 作业完成。
5. 使用 [`test/compare.py`](compare.py) 对基准结果和测试结果进行可视化比较。

[`test/test.sh`](test.sh) 只提交作业，不等待作业结束，也不自动调用 [`test/compare.py`](compare.py)。

用户获得的测试相关文件如下；运行 [`test/test.sh`](test.sh) 并解压基准压缩包后，会生成 `baseline/`：

```text
test/
├── README.md
├── baseline.tar.gz
├── compare.py
└── test.sh
```

解压后，`baseline/` 的主要结构如下：

```text
test/
├── README.md
├── baseline/
│   ├── mnu000/
│   ├── mnu006/
│   └── noise/
├── compare.py
└── test.sh
```

### `baseline/`

[`test/baseline/`](baseline/) 保存解压后的两套基准模拟结果。每套结果目前包含：

- `image1/` 到 `image8/` 下的密度场文件 `0.000_delta_c_*.bin`。
- `image1/` 下的二维投影文件 `0.000_proj_xy_1.bin`。
- `image1/` 下的功率谱二进制文件 `0.000_cicpower_1.bin`。
- `neutrinos/` 下的功率谱和相关输入文件。
- `parameters.f90`、`run.qsub` 和 `run.log` 等运行记录。
- `fig/` 下的已有结果图。

两套基准参数设置是 `box = 600`、`ng = 128`、`nn = 2`，不要将本目录说明中的参数与其他测试方案混用。

`baseline/noise/` 下保存 `image1` 到 `image8` 对应的噪声输入。根目录下的 `power_11.jpg` 以及各模拟目录中的 `fig/` 保存已有结果图。

### [`test/compare.py`](compare.py)

该脚本使用 NumPy 读取二进制密度场和功率谱，使用 Matplotlib 生成比较图。脚本中的路径由 `code_path` 计算得到，运行 [`test/test.sh`](test.sh) 时会自动更新：

- 基准结果：`<code_path>/test/baseline/`。
- 待比较结果：由 `compare.py` 中的 `result_path` 确定，运行前需确认其指向实际测试结果目录。
- 图像输出：`<code_path>/test/fig/`。

脚本比较以下内容：

- `mnu000` 和 `mnu006` 的 `z = 0` 二维投影图。
- 两套模拟的功率谱与基准功率谱。
- 功率谱相对差异。
- `image1` 到 `image8` 的三维密度场差异及相对差异。

脚本包含需要在交互式 Python 环境中执行的单元格标记，并会调用 `plt.show()` 显示图像。它不是自动化的 PASS/FAIL 比较器，也不会返回基于数值容差的回归测试结论。

## 3. PBS 作业脚本

[`test/baseline/mnu000/run.qsub`](baseline/mnu000/run.qsub) 和 [`test/baseline/mnu006/run.qsub`](baseline/mnu006/run.qsub) 的作业配置基本一致：

- PBS 资源：1 个节点、64 个 CPU、30 GB 内存、`normal` 队列。
- Coarray image 数：`FOR_COARRAY_NUM_IMAGES=8`。
- OpenMP 线程数：`OMP_NUM_THREADS=32`。
- 并行进程数：`PROCS=8`。
- 编译环境：Intel 2018、MPICH、MKL、FFTW，以及脚本中列出的 HDF5、Boost 和 Python 模块。
- 编译器：`/opt/intel2018update3/bin/ifort`。

作业脚本执行的主要步骤是：

1. 清理对应的 `/home/cossim/cube_nu/baseline_build/work-build/mnu000/` 或 `mnu006/`。
2. 将当前仓库源码复制到该隔离工作目录。
3. 复制对应的 `parameters.f90`、neutrino 输入和 noise 文件。
4. 编译主程序、初始条件程序和功率谱程序。
5. 依次运行初始条件程序、主程序和 `cicpower` 程序。
6. `mnu000` 作业还会修改主程序中的检查点编号并追加一次断点续传运行。

作业脚本使用固定的绝对路径 `/home/cossim/cube_nu/...`，因此只能在相应目录布局和平台环境下直接使用。提交前应检查 PBS 环境、模块名称、绝对路径、编译器和资源配置。

## 4. 运行方法

在仓库根目录执行：

```bash
# 先修改 test.sh 中的 code_path，再执行：
bash test/test.sh
```

该命令等价于依次提交：

```bash
qsub test/baseline/mnu000/run.qsub
qsub test/baseline/mnu006/run.qsub
```

作业完成后，确认 [`test/compare.py`](compare.py) 中的 `result_path` 指向实际待比较结果目录，再在配置了 Python、NumPy 和 Matplotlib 的环境中运行：

```bash
python test/compare.py
```

当前 `compare.py` 中的 `result_path` 默认指向仓库同级目录下的 `test/`。实际运行时，应根据 PBS 作业的输出位置确认该路径，例如：

```text
<仓库父目录>/test/
```

如果作业输出位于其他目录，需要先修改 [`test/compare.py`](compare.py) 中的路径变量；该修改属于测试结果路径配置，不应覆盖解压得到的 [`test/baseline/`](baseline/) 基准数据。

