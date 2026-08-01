# %% 导入模块
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


code_path = Path(__file__).resolve().parents[2]
baseline_path = code_path / "test" / "baseline" / "serial"
result_path = code_path.parent / "test" / "serial/image1/"
fig_path = code_path / "test"/ "serial" / "fig"

sim = "serial"
boxsize = 600

fig_path.mkdir(parents=True, exist_ok=True)


def require_file(filename):
    path = Path(filename)
    if not path.is_file():
        raise FileNotFoundError(f"找不到文件: {path}")
    return path


def loadfield2d(filename):
    """读取方形 float32 二维场。"""
    filename = require_file(filename)
    field = np.fromfile(filename, dtype=np.float32)
    n = round(len(field) ** 0.5)
    if n * n != len(field):
        raise ValueError(f"二维场数据长度不是完全平方数: {filename}")
    field = field.reshape(n, n)
    print(f"{filename}: n={n}")
    return field


def loadfield3d(filename, ng=None):
    """读取 float32 三维密度场。"""
    filename = require_file(filename)
    field = np.fromfile(filename, dtype=np.float32)
    if ng is None:
        n = round(len(field) ** (1.0 / 3.0))
        if n * n * n != len(field):
            raise ValueError(f"三维场数据长度不是完全立方数: {filename}")
        return field.reshape(n, n, n), n

    expected = ng * ng * (ng + 2)
    if len(field) != expected:
        raise ValueError(f"三维场长度不匹配: {filename}, 实际 {len(field)}, 期望 {expected}")
    field = field.reshape(ng, ng, ng + 2)
    return field[:, :, :ng], ng


def loadpower(filename):
    """读取每行 10 个 float32 数值的功率谱文件。"""
    filename = require_file(filename)
    values = np.fromfile(filename, dtype=np.float32)
    n_row = 10
    if len(values) % n_row != 0:
        raise ValueError(f"功率谱数据长度不能按 10 列重塑: {filename}")
    values = values.reshape(-1, n_row)
    return values[:, 1], values


# %% 投影可视化
xgrid = [0, boxsize]
proj_xy = loadfield2d(result_path / "0.000_proj_xy_1.bin")
    
plt.figure(figsize=(6,5),dpi=150)
plt.imshow(proj_xy.T, cmap='gray', extent=[xgrid[0], xgrid[1], xgrid[0], xgrid[1]])

ax = plt.gca()
ax.spines['top'].set_visible(True)
ax.spines['right'].set_visible(True)

plt.colorbar()
plt.clim(-1, 3)
plt.title('$z=0.0\\,\\,\\,\\delta_c$')
plt.savefig(fig_path / 'proj_xy.jpg', format='jpeg')

plt.show()


# %% 功率谱比较
baseline_power = loadpower(
    baseline_path / "0.000_cicpower_1.bin"
)
result_power = loadpower(
    result_path /  "0.000_cicpower_1.bin"
)

if len(baseline_power[0]) != len(result_power[0]):
    raise ValueError("基准和串行结果的功率谱长度不一致")

Path0=baseline_path/ "../mnu000"
k = np.loadtxt(Path0 /'neutrinos/k_values.txt')
pk0 = np.loadtxt(Path0 /'neutrinos/Pk_cb_0.0000.txt')
ksim_baseline, xi_baseline = baseline_power
ksim_result, xi_result = result_power

fig, axes = plt.subplots(2, 1, figsize=(7, 9), dpi=150, constrained_layout=True)
axes[0].loglog(k, pk0, '--', linewidth=5, color=[.7, .7, .7],label="CAMB nonlinear")
axes[0].loglog(ksim_baseline, xi_baseline[:, 4], label="serial baseline")
axes[0].loglog(ksim_result, xi_result[:, 4], "--", label="serial result")
axes[0].set_xlabel(r"$k\,[h/{\rm Mpc}]$")
axes[0].set_ylabel(r"$P(k)$")
axes[0].set_title(r"$z=0.0\quad P_{\rm matter}$")
axes[0].grid(True, which="both", linewidth=0.3)
axes[0].legend()

relative_power_difference = np.abs(xi_result[:, 4] / xi_baseline[:, 4] - 1.0)
axes[1].loglog(ksim_baseline, relative_power_difference)
axes[1].set_xlabel(r"$k\,[h/{\rm Mpc}]$")
axes[1].set_ylabel(r"$|\Delta P/P|$")
axes[1].set_title("power spectrum relative difference")
axes[1].grid(True, which="both", linewidth=0.3)

fig.savefig(fig_path / "serial_power.jpg", format="jpeg")
plt.show()


# %% 密度场绝对差值
baseline_field, n_baseline = loadfield3d(
    baseline_path / "0.000_delta_c_1.bin"
)
result_field, n_result = loadfield3d(
    result_path / "0.000_delta_c_1.bin"
)

if n_baseline != n_result or baseline_field.shape != result_field.shape:
    raise ValueError("基准和串行结果的密度场形状不一致")

difference = np.abs(baseline_field - result_field)
print("serial density field absolute difference:")
print(f"  max  = {np.max(difference):.8e}")
print(f"  mean = {np.mean(difference):.8e}")
print(f"  std  = {np.std(difference):.8e}")

fig, axis = plt.subplots(figsize=(6, 5), dpi=150, constrained_layout=True)
difference_projection = np.mean(difference, axis=2)
image_plot = axis.imshow(
    difference_projection.T,
    cmap="gray",
    extent=[0, boxsize, 0, boxsize],
    origin="lower",
)
axis.set_title("serial baseline-result absolute difference")
axis.set_xlabel("x [Mpc/h]")
axis.set_ylabel("y [Mpc/h]")
fig.colorbar(image_plot, ax=axis, label="absolute difference")
fig.savefig(fig_path / "serial_density_absolute_difference.jpg", format="jpeg")
plt.show()
