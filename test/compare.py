# %% 导入模块
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


code_path = Path(__file__).resolve().parents[1]
baseline_path = code_path / "test" / "baseline"
result_path = code_path.parent / "test"
fig_path = code_path / "test" / "fig"

sim0 = "mnu000"
sim6 = "mnu006"
has_sim6 = (result_path / sim6 / "image1" / "0.000_proj_xy_1.bin").exists()
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
    field = field.reshape(n, n)
    print(f"{filename}: n={n}")
    return field


def loadpower(filename):
    """读取每行 10 个 float32 数值的功率谱文件。"""
    filename = require_file(filename)
    values = np.fromfile(filename, dtype=np.float32)
    n_row = 10
    values = values.reshape(-1, n_row)
    return values[:, 1], values


# %% 投影可视化
def plot_xy(sim):
    xgrid = [0, boxsize]
    proj_xy = loadfield2d(
        result_path / sim / "image1" / "0.000_proj_xy_1.bin"
    )

    fig, axis = plt.subplots(figsize=(6, 5), dpi=150, constrained_layout=True)
    axis.imshow(
        proj_xy.T,
        cmap="gray",
        extent=[xgrid[0], xgrid[1], xgrid[0], xgrid[1]],
    )
    axis.spines["top"].set_visible(True)
    axis.spines["right"].set_visible(True)
    fig.colorbar(axis.images[0], ax=axis)
    axis.images[0].set_clim(-1, 3)
    axis.set_title(r"$z=0.0\,\,\,\delta_c$")
    fig.savefig(fig_path / f"{sim}_proj_xy.jpg", format="jpeg")
    plt.close(fig)


plot_xy(sim0)
if has_sim6:
    plot_xy(sim6)


# %% 功率谱比较
baseline_path0 = baseline_path / sim0
result_path0 = result_path / sim0

k = np.loadtxt(baseline_path0 / "neutrinos" / "k_values.txt")
pk0 = np.loadtxt(baseline_path0 / "neutrinos" / "Pk_cb_0.0000.txt")

ksim0, xi0 = loadpower(
    baseline_path0 / "image1" / "0.000_cicpower_1.bin"
)
ksim0_test, xi0_test = loadpower(
    result_path0 / "image1" / "0.000_cicpower_1.bin"
)

if has_sim6:
    baseline_path6 = baseline_path / sim6
    result_path6 = result_path / sim6
    pk6 = np.loadtxt(baseline_path6 / "neutrinos" / "Pk_cb_0.0000.txt")
    ksim6, xi6 = loadpower(
        baseline_path6 / "image1" / "0.000_cicpower_1.bin"
    )
    ksim6_test, xi6_test = loadpower(
        result_path6 / "image1" / "0.000_cicpower_1.bin"
    )

fig, axes = plt.subplots(2, 1, figsize=(6, 8), dpi=150, constrained_layout=True)
axes[0].loglog(k, pk0, "--", linewidth=5, color=[0.7, 0.7, 0.7], label="CAMB nonlinear")
axes[0].loglog(ksim0, xi0[:, 4], label=r"$P_{000,baseline}(k)$")
if has_sim6:
    axes[0].loglog(ksim6, xi6[:, 4], label=r"$P_{006,baseline}(k)$")
axes[0].loglog(ksim0_test, xi0_test[:, 4], "--", label=r"$P_{000,test}(k)$")
if has_sim6:
    axes[0].loglog(ksim6_test, xi6_test[:, 4], "--", label=r"$P_{006,test}(k)$")
axes[0].grid(True, which="both", linewidth=0.3)
axes[0].set_ylim(1e-4, 1e5)
axes[0].set_xlabel(r"$k\,[h/{\rm Mpc}]$")
axes[0].set_ylabel(r"$P(k)$")
axes[0].legend(loc="lower left")
axes[0].set_title(r"$z=0.0\quad power_{matter}$")

axes[1].loglog(ksim0, np.abs(xi0_test[:, 4] / xi0[:, 4] - 1), label=r"$|\Delta P_{000}/P|$")
if has_sim6:
    axes[1].loglog(ksim6, np.abs(xi6_test[:, 4] / xi6[:, 4] - 1), label=r"$|\Delta P_{006}/P|$")
axes[1].set_ylim(1e-5, 1e-1)
axes[1].grid(True, which="both", linewidth=0.5)
axes[1].set_xlabel(r"$k\,[h/{\rm Mpc}]$")
axes[1].set_ylabel(r"$\Delta P(k)$")
axes[1].spines["top"].set_visible(True)
axes[1].spines["right"].set_visible(True)
if has_sim6:
    axes[1].legend()

fig.savefig(fig_path / "power.jpg", format="jpeg")
plt.close(fig)


# %% 密度场误差

def loadfield3d(filename, ng=None):
    """读取 float32 三维密度场。"""
    filename = require_file(filename)
    field = np.fromfile(filename, dtype=np.float32)

    if ng is None:
        n = round(len(field) ** (1.0 / 3.0))
        field = field.reshape(n, n, n)
    else:
        field = field.reshape(ng, ng, ng + 2)
        field = field[:ng, :ng, :ng]
        n = ng

    return field, n


def compare_field(sim):
    for image_number in range(1, 9):
        print(f"image{image_number}")
        baseline_field, n_baseline = loadfield3d(
            baseline_path
            / sim
            / f"image{image_number}"
            / f"0.000_delta_c_{image_number}.bin"
        )
        result_field, n_result = loadfield3d(
            result_path
            / sim
            / f"image{image_number}"
            / f"0.000_delta_c_{image_number}.bin"
        )

        difference = np.abs(baseline_field - result_field)
        maximum = np.max(difference)
        mean = np.mean(difference)
        standard_deviation = np.std(difference)
        print(
            " baseline-test:\n"
            f"    max:{maximum:.6f} "
            f"mean:{mean:.6f} std:{standard_deviation:.6f}"
        )

        if mean > 3 or standard_deviation > 15:
            print("!" * 20 + f"\n image{image_number} too large error\n" + "!" * 20)

        fig, axis = plt.subplots(figsize=(6, 5), dpi=150, constrained_layout=True)
        image_plot = axis.imshow(
            np.mean(difference, axis=2).reshape(n_baseline, n_baseline).T,
            cmap="gray",
        )
        image_plot.set_clim(-0.1, 0.1)
        axis.set_title(f"baseline-test image{image_number}")
        fig.colorbar(image_plot, ax=axis)
        fig.savefig(fig_path / f"{sim}_d-d.jpg", format="jpeg")
        plt.close(fig)


compare_field(sim0)
if has_sim6:
    compare_field(sim6)


# %% 密度场视觉检验
image_number = 8
sim = sim0
baseline_field, n = loadfield3d(
    baseline_path
    / sim
    / f"image{image_number}"
    / f"0.000_delta_c_{image_number}.bin"
)
result_field, n = loadfield3d(
    result_path
    / sim
    / f"image{image_number}"
    / f"0.000_delta_c_{image_number}.bin"
)

fig, axis = plt.subplots(figsize=(6, 5), dpi=150, constrained_layout=True)
image_plot = axis.imshow(
    np.mean(baseline_field, axis=2).reshape(n, n).T,
    cmap="gray",
)
image_plot.set_clim(-1, 3)
axis.set_title(f"baseline image{image_number}")
fig.colorbar(image_plot, ax=axis)
fig.savefig(fig_path / f"{sim}_d.jpg", format="jpeg")
plt.close(fig)

fig, axis = plt.subplots(figsize=(6, 5), dpi=150, constrained_layout=True)
image_plot = axis.imshow(
    np.mean(result_field, axis=2).reshape(n, n).T,
    cmap="gray",
)
image_plot.set_clim(-1, 3)
axis.set_title(f"test image{image_number}")
fig.colorbar(image_plot, ax=axis)
fig.savefig(fig_path / f"{sim}_d.jpg", format="jpeg")
plt.close(fig)


# %%
