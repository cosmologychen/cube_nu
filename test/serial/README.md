# CUBE 串行测试

本目录提供 `mnu000` 的单 image、单线程测试入口。测试不会修改现有 [`test/baseline/`](../baseline/)；运行结果写入仓库同级的 `test/serial/`。

## 测试参数

参数文件 [`parameters.f90`](parameters.f90) 基于 [`test/baseline/mnu000/parameters.f90`](../baseline/mnu000/parameters.f90)，仅将并行规模设置为：

- `box = 600`
- `ng = 128`
- `nn = 1`
- `ncore = 1`
- `nteam = 1`
- `nnest = 1`
- `m_nu = [0.0, 0.0, 0.0]`

程序仍然使用 Coarray Fortran，但运行时只启动 1 个 image；OpenMP 线程数设置为 1。

## 输入

串行测试使用 [`test/baseline/noise/image1/`](../baseline/noise/image1/) 中的：

- `noise_1.bin`
- `seed_1.bin`

中微子相关输入从 [`test/baseline/mnu000/neutrinos/`](../baseline/mnu000/neutrinos/) 复制。

## 运行

在仓库根目录执行：

```bash
bash test/serial/run.sh
```

脚本将：

1. 清理并重新创建仓库同级的 `test/serial/`；
2. 将当前源码复制到 `test/baseline/serial/code/serial/`；
3. 覆盖复制串行参数文件；
4. 复制 `image1` 的 noise 和 seed 输入；
5. 编译主程序及工具程序；
6. 依次运行初始条件程序、主程序和 CIC 功率谱程序。

输出目录为：

```text
test/baseline/serial/
├── code/serial/
└── mnu000/
```

## 编译环境

编译环境需要与项目现有测试一致，包括支持 Coarray Fortran 的 Intel Fortran 编译器、FFTW 和 MKL。脚本不加载模块，也不覆盖当前环境中的 `FC`、`XFLAG`、`OFLAG` 和 `FFTFLAG`；这些变量需要由调用环境提供。

当前脚本只执行单 image、单 OpenMP 线程配置，不包含 PBS 提交，也不执行结果自动比较。
