"""
analysis_run.py - 宇宙学模拟分析运行脚本

本脚本将 Jupyter Notebook 中的分析步骤转换为可独立运行的 Python 脚本。
使用 # %% 分隔符，可以在 IDE 中像 notebook 一样运行单元格。

使用方法:
    在 IDE 中: 点击单元格上方的"运行单元格"按钮
    命令行: python analysis_run.py
"""
# %%  
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import zoom

from analysis_core import (
    loadsim, loadz, loadfield2d, loadfield3d, 
    loadpower, loadtcat, cambpower,
    plt_proj, plt_power, plt_tcat, plt_delta_field, plt_phi_field
)


# %% 单元格 0: 加载模拟参数和红移数据
# 单元格 0: 加载模拟参数和红移数据
print("=" * 60)
print("单元格 0: 加载模拟参数和红移数据")
print("=" * 60)

Path = '/home/cbh/cube_nu/'
Dir = 'output/canp600_512_1_0.0_1/image1/'

sim = loadsim(Path + Dir + 'parameters.txt')
z = loadz(Path + Dir + 'z_checkpoint.txt')

print(f"模拟参数: {sim}")
print(f"红移数据点数: {len(z)}")


# %% 单元格 1: 绘制投影图
# 单元格 1: 绘制投影图
print("\n" + "=" * 60)
print("单元格 1: 绘制投影图")
print("=" * 60)

ng = sim['ng']
boxsize = sim['boxsize']
xgrid = [0, boxsize]

Redshift_i = ['0.000']
for Redshift in Redshift_i:
    plt_proj(Redshift, Path, Dir, xgrid)


# %% 单元格 2: 绘制功率谱
# 单元格 2: 绘制功率谱
print("\n" + "=" * 60)
print("单元格 2: 绘制功率谱")
print("=" * 60)

plt_power('200.000', Path, Dir, sim)


# %% 单元格 3: 加载3D场数据
# 单元格 3: 加载3D场数据
print("\n" + "=" * 60)
print("单元格 3: 加载3D场数据")
print("=" * 60)

ng = 320
phi, n = loadfield3d(
    '/home/cossim/cube_nu/test/200_1024/sigle_nu_200_512_2_0.06_2/image4/1.157_error_phi_4.bin',
    ng=ng
)

plt.figure()
plt.imshow(np.mean(phi[:, :, :], axis=1).reshape(ng, ng).T, cmap='gray')
plt.colorbar()
plt.show()


# %% 单元格 4: 绘制密度场
# 单元格 4: 绘制密度场
print("\n" + "=" * 60)
print("单元格 4: 绘制密度场")
print("=" * 60)

ng = sim['ng']
Redshift = '200.000'

delta_L, n = loadfield3d(Path + Dir + Redshift + '_delta_c_1.bin')

plt.figure()
plt.imshow(np.mean(delta_L[:, :, :], axis=2).reshape(n, n).T, cmap='gray')
plt.colorbar()
plt.clim(-3, 3)
plt.title('$\\delta_c z=200.0$')
plt.show()


# %% 单元格 5: 绘制初始场
# 单元格 5: 绘制初始场
print("\n" + "=" * 60)
print("单元格 5: 绘制初始场")
print("=" * 60)

ng = sim['ng']
delta_L, n = loadfield3d(Path + Dir + 'delta_L_1.bin')
delta_L = delta_L / 0.0063

plt.figure()
plt.imshow(np.mean(delta_L[:, :, :], axis=2).reshape(n, n).T, cmap='gray')
plt.colorbar()
plt.clim(-3, 3)
plt.title('$\\delta_L$')
plt.show()


# %% 单元格 6: 绘制势能场
# 单元格 6: 绘制势能场
print("\n" + "=" * 60)
print("单元格 6: 绘制势能场")
print("=" * 60)

ng = sim['ng']
Redshift = '200.000'

plt_phi_field(Redshift, Path, Dir, ng)


# %% 单元格 7: 绘制时间统计
# 单元格 7: 绘制时间统计
print("\n" + "=" * 60)
print("单元格 7: 绘制时间统计")
print("=" * 60)

nsnap_read = [100, 98, 96, 95, 92, 84, 81, 72, 66, 38, 1]
plt_tcat(z, nsnap_read, Path, Dir)


# %% 单元格 8: 比较不同分辨率的场
# 单元格 8: 比较不同分辨率的场
print("\n" + "=" * 60)
print("单元格 8: 比较不同分辨率的场")
print("=" * 60)

delta_L, ng = loadfield3d('../output/canp600_512_1_0.0_1/image1/25.000_delta_c_1.bin')
delta_C, nw = loadfield3d('../output/canp600_512_1_0.0_1/image1/25.000_delta_coarse_1.bin')

delta_L_resampled = zoom(delta_L, nw/ng, order=1)

plt.figure()
plt.imshow(np.mean(delta_L[:, :, :], axis=2).reshape(ng, ng).T, cmap='gray')
plt.colorbar()
plt.clim(-3, 3)
plt.title('$\\delta_{cicpower} z=25.0$')
plt.show()

plt.figure()
plt.imshow(np.mean(delta_L_resampled[:, :, :], axis=2).reshape(nw, nw).T, cmap='gray')
plt.colorbar()
plt.clim(-3, 3)
plt.title('$\\delta_{cicpower_zoom} z=25.0$')
plt.show()

plt.figure()
plt.imshow(np.mean(delta_C[:, :, :], axis=2).reshape(nw, nw).T, cmap='gray')
plt.colorbar()
plt.clim(-3, 3)
plt.title('$\\delta_{coarse} z=25.0$')
plt.show()


# %% 单元格 9: 比较场差异
# 单元格 9: 比较场差异
print("\n" + "=" * 60)
print("单元格 9: 比较场差异")
print("=" * 60)

delta_L, ng = loadfield3d('../output/canp600_512_1_0.0_1/image1/25.000_delta_c_1.bin')
delta_C, nw = loadfield3d('../output/canp600_512_1_0.0_1/image1/25.000_delta_coarse_1.bin')

delta_L_resampled = zoom(delta_L, nw/ng, order=1)
difference = delta_L_resampled - delta_C

plt.figure()
plt.imshow(np.mean(difference[:, :, :], axis=2).reshape(nw, nw).T, cmap='gray')
plt.colorbar()
plt.title('$difference z=25.0$')
plt.show()

plt.figure()
plt.imshow(np.mean(difference[:, :, :] / delta_C - 1, axis=2).reshape(nw, nw).T, cmap='gray')
plt.colorbar()
plt.clim(-1, 1)
plt.title('$diff z=25.0$')
plt.show()


# %% 单元格 10: 加载外部场数据
# 单元格 10: 加载外部场数据
print("\n" + "=" * 60)
print("单元格 10: 加载外部场数据")
print("=" * 60)

Redshift = '200.000'
delta_L, nc = loadfield3d('/home/cbh/visualization/plane/rhoc.bin')

plt.figure()
plt.imshow(np.mean(delta_L[:, :, :], axis=2).reshape(nc, nc).T, cmap='gray')
plt.colorbar()
plt.title('$\\delta_{NFW}$')
plt.show()


# %% 单元格 11: 比较全局delta场
# 单元格 11: 比较全局delta场
print("\n" + "=" * 60)
print("单元格 11: 比较全局delta场")
print("=" * 60)

delta_1, nc = loadfield3d('../output/gpower_600_256_1_0.1_2/image1/0.020_delta_global1_1.bin')
delta_2, nc = loadfield3d('../output/gpower_600_256_1_0.1_2/image1/0.020_delta_global2_1.bin')

plt.figure()
plt.imshow(np.mean(delta_1[:, :, :], axis=2).reshape(nc, nc).T, cmap='gray')
plt.colorbar()
plt.title('$\\delta_1$')
plt.show()

plt.figure()
plt.imshow(np.mean(delta_2[:, :, :], axis=2).reshape(nc, nc).T, cmap='gray')
plt.colorbar()
plt.title('$\\delta_2$')
plt.show()
