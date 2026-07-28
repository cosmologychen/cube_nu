# CUBE 中微子宇宙学模拟程序

## 概述

CUBE 是一个用于中微子宇宙学模拟的并行程序，基于 Fortran 2008 编写，支持 OpenMP 和 MPI 并行计算，采用 Coarray 分布式内存编程模型。该程序主要用于模拟含有中微子的宇宙结构形成过程。

## 系统要求

### 编译环境
- **Fortran 编译器**: Intel Fortran Compiler 2018u3
- **MPI 实现**: MPICH 3.2.1
- **数学库**: Intel MKL
- **FFT 库**: FFTW 3.3.8
- **HDF5**: 1.12.1 (with C++ and Fortran 支持)
- **Boost**: 1.73.0
- **Python**: 3.8.5

### 硬件配置
- PBS 作业调度系统
- 推荐内存: 30GB+
- 支持多节点并行

## 编译

```bash
# 清理并编译主程序
make clean
make

# 编译工具程序 (位于 utilities 目录)
cd utilities
make ic.x        # 初始条件生成程序
make cicpower.x  # CIC 功率谱计算程序
make fof.x       # 晕查找程序
make dsp.x       # 位移场计算程序
make cdm2nu.x    # CDM 到中微子转换程序
make gadget.x    # Gadget 格式输出程序
make int2real.x  # 整型转实型程序
cd ..
```

## 运行

### 基本用法

```bash
./RUN [选项]
```

### 主要运行选项

#### 运行模式控制
| 选项 | 说明 | 默认值 |
|------|------|--------|
| `-IC <true/false>` | 运行初始条件生成 | true |
| `-Main <true/false>` | 运行主模拟进程 | true |
| `-Halo <true/false>` | 运行晕查找 | false |
| `-Power <true/false>` | 计算 CIC 功率谱 | false |
| `-DSP <true/false>` | 计算位移场 | false |
| `-C2u <true/false>` | CDM2Nu 转换 | false |
| `-Gad <true/false>` | Gadget 输出 | true |

#### 路径配置
| 选项 | 说明 | 默认值 |
|------|------|--------|
| `-c <路径>` | 代码路径 | /home/cossim/cube_nu/test/output/code |
| `-o <路径>` | 输出路径 | /home/cossim/cube_nu/test/output |
| `-z <路径>` | 红移检查点文件 | 无 |
| `-h <路径>` | 晕查找检查点文件 | 无 |
| `-s <路径>` | 随机种子路径 | 无 |
| `-n <路径>` | 噪声图路径 | 无 |
| `-r <路径>` | 密度场路径 | 无 |

#### 模拟参数
| 选项 | 说明 | 默认值 |
|------|------|--------|
| `-box <num>` | 盒子大小 (Mpc/h) | 1200 |
| `-ng <num>` | 每维粒子数 | 1024 |
| `-nn <num>` | 每维图像数 | 1 |
| `-ncore <num>` | 每节点核心数 | 64 |
| `-nnest <num>` | 嵌套数 | 4 |

#### 宇宙学参数
| 选项 | 说明 | 默认值 |
|------|------|--------|
| `-h0 <num>` | Hubble 参数 h0 | 0.6817 |
| `-omc <num>` | 冷暗物质密度 Ω_c | 0.253209016 |
| `-omb <num>` | 重子物质密度 Ω_b | 0.048072486 |
| `-ns <num>` | 谱指数 n_s | 0.9693 |
| `-As <num>` | 功率谱振幅 A_s | 2.122e-09 |
| `-mnus <num>` | 三个中微子质量 (eV) | 0.0 0.0 0.0 |
| `-Mnu <num>` | 中微子总质量 (eV) | 0.0 |
| `-Hnu <num>` | 中微子层级 (0=DG, 1=IH) | 0 |
| `-Cpk <num>` | 功率谱算法 | 2 |

#### 功率谱算法 (-Cpk)
- `-1`: 完全使用 CAMB
- `0`: 部分使用 CAMB
- `1`: 使用拼接功率谱
- `2`: 使用全局功率谱

#### 力计算选项
| 选项 | 说明 | 默认值 |
|------|------|--------|
| `-PM3 <true/false>` | 启用 PM3 力计算 | true |
| `-PP <true/false>` | 启用 PP 力计算 | true |

#### 其他选项
| 选项 | 说明 | 默认值 |
|------|------|--------|
| `-ZIPX <0/1>` | 0=真实坐标, 1=压缩坐标 | 1 |
| `-z1 <int>` | 晕/功率计算起始点 | 1 |
| `-z2 <int>` | 晕/功率计算终止点 | n_checkpoint |
| `-flag <name>` | 运行标识 | 空 |
| `-log <name>` | 日志文件名 | run.log |

## 运行示例

### 完整模拟流程
```bash
./RUN -flag myrun -box 1200 -ng 1024 -nn 1
```

### 仅生成初始条件
```bash
./RUN -flag -RunIC only
```

### 从检查点恢复运行
```bash
./RUN -flag myrun -checkpoint 5
```

### 设置中微子质量
```bash
./RUN -Mnu 0.06 -mnus 0.02 0.02 0.02 -Cpk 2
```

### 改变宇宙学参数
```bash
./RUN -h0 0.67 -omc 0.25 -omb 0.05 -ns 0.96
```

## 程序模块

### 初始条件 (IC)
生成宇宙学模拟的初始条件，包括:
- 粒子位置和速度
- 随机种子读取/生成
- 噪声图读取/生成
- 密度场读取/生成

### 主模拟 (Main)
主要的 N体/粒子网格模拟程序，包含:
- PM3 (Particle-Mesh 3D) 力计算
- PP (Particle-Particle) 力计算
- 时间积分 (kick-drift-kick)

### CIC 功率谱
使用 CIC (Cloud-in-Cell) 方法计算功率谱:
- `cicpower_global.f90`: 全局功率谱计算
- `cicpower_segment.f90`: 分段功率谱计算

### 晕查找 (FOF)
使用 Friends-of-Friends 算法进行暗晕识别

### Gadget 输出
将模拟结果输出为 Gadget 格式

## 目录结构

```
cube_nu_new/
├── RUN              # 主运行脚本
├── run.qsub         # PBS 作业提交脚本
├── run.sh           # 快速运行脚本
├── Makefile         # 主程序编译配置
├── main.f90         # 主程序入口
├── initialize.f90   # 初始化模块
├── kick.f90         # 时间积分模块
├── drift.f90        # 漂移模块
├── checkpoint.f90   # 检查点模块
├── parameters.f90   # 模拟参数定义
├── variables.f90    # 变量定义
├── basic_functions.f08  # 基础函数库
├── buffer_grid.f90  # 网格缓冲区
├── buffer_particle.f90   # 粒子缓冲区
├── pencil_fft.f90   # Pencil FFT
├── pencil_fft_global.f90  # 全局 Pencil FFT
├── cubefft.f90      # Cube FFT
├── Green/           # Green 函数相关
├── neutrinos/       # 中微子相关模块
├── utilities/       # 工具程序
│   ├── ic.f90       # 初始条件生成
│   ├── cicpower.f90 # 功率谱计算
│   ├── fof.f90      # 晕查找
│   ├── dsp.x        # 位移场
│   ├── cdm2nu.f90   # CDM2Nu 转换
│   ├── gadget.x    # Gadget 输出
│   └── ...
└── z_checkpoint.txt # 红移检查点文件
```

## 环境变量

运行前需设置以下环境变量:
```bash
export FC='/opt/intel2018update3/bin/ifort'
export XFLAG_NO_OMP='-O3 -fpp -q -coarray=distributed -mcmodel=large'
export XFLAG='-O3 -fpp -qopenmp -coarray=distributed -mcmodel=large'
export FFTFLAG='-I/opt/intel2018update3/mkl/include/fftw -mkl'
export LD_LIBRARY_PATH="/opt/intel2018update3/mkl/lib/intel64:/opt/fftw/3.3.8/lib:..."
export OMP_STACKSIZE=2G
export OMP_NUM_THREADS=64
export FOR_COARRAY_NUM_IMAGES=8
```

## PBS 作业配置

默认 PBS 配置:
- 队列: normal
- 每节点任务数: 2
- 每节点核心数: 64
- 内存: 300GB
- 嵌套数: 4

## 输出

运行后将在输出目录生成:
- `run.log`: 运行日志
- `image*/`: 各图像数据目录
  - `*_xp_*.bin`: 粒子位置文件
  - `*_vk_*.bin`: 粒子速度文件
- `checkpoint/`: 检查点文件
- `powerSpectrum/`: 功率谱数据

## 注意事项

1. 确保所有依赖库已正确加载 (通过 `module load`)
2. 检查点文件 `z_checkpoint.txt` 定义了模拟的红移点
3. 中微子质量可通过 `-mnus` 或 `-Mnu` 设置
4. 使用 `-Y` 选项可跳过所有确认提示
5. 运行时需确保输出路径有足够磁盘空间