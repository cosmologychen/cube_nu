# %% 导入模块
import os
import numpy as np
import matplotlib.pyplot as plt
def loadfield2d(fn):
    """
    加载2D场数据
    
    Parameters:
        fn: filename
        
    Returns:
        a: 2D field
    """
    fid = open(fn, 'rb')
    p1 = np.fromfile(fid, dtype=np.float32)
    fid.close()
    n = round(len(p1) ** 0.5)
    a = np.reshape(p1, (n, n))
    print(fn, n)
    return a

def loadpower(filename):
    """
    加载功率谱数据
    
    Parameters:
        filename: 功率谱文件名
        
    Returns:
        ksim: k值数组
        xi: 功率谱数据
    """
    print(filename)
    n_row_xi = 10
    fid = open(filename, 'rb')
    xi = np.fromfile(fid, dtype='float32')
    fid.close()
    xi = np.reshape(xi, (int(len(xi) / n_row_xi), n_row_xi))
    ksim = xi[:, 1]
    return ksim, xi


code_path='/home/cbh/18T/cube_nus/cube_nu/'
baseline_path=code_path+'/test/baseline/'
result_path=code_path+'../test/'
fig_path=code_path+'/test/fig/'
sim0='mnu000'
sim6='mnu006'

os.makedirs(fig_path,exist_ok=True)

# %% 单元格 1: 绘制投影图
# 单元格 1: 绘制投影图
print("\n" + "=" * 60)
print("单元格 1: 绘制投影图")
print("=" * 60)

ng = 256
boxsize = 600
def plot_xy(sim):
    xgrid = [0, boxsize]
    proj_xy = loadfield2d(result_path+sim +'/image1/0.000_proj_xy_1.bin')
        
    plt.figure(figsize=(6,5),dpi=150)
    plt.imshow(proj_xy.T, cmap='gray', extent=[xgrid[0], xgrid[1], xgrid[0], xgrid[1]])

    ax = plt.gca()
    ax.spines['top'].set_visible(True)
    ax.spines['right'].set_visible(True)

    plt.colorbar()
    plt.clim(-1, 3)
    plt.title('$z=0.0\\,\\,\\,\\delta_c$')
    plt.savefig(fig_path +sim+ '_proj_xy.jpg', format='jpeg')

    plt.show()

plot_xy(sim0)
plot_xy(sim6)

# %% 单元格子 2: 绘制功率谱
Path0=baseline_path+sim0
Path6=baseline_path+sim6
Path0_test=result_path+sim0
Path6_test=result_path+sim6
k = np.loadtxt(Path0 +'/neutrinos/k_values.txt')
pk0 = np.loadtxt(Path0 +'/neutrinos/Pk_cb_0.0000.txt')
pk6 = np.loadtxt(Path6 +'/neutrinos/Pk_cb_0.0000.txt')
ksim0, xi0 = loadpower(Path0 + '/image1/0.000_cicpower_1.bin')
ksim6, xi6 = loadpower(Path6 + '/image1/0.000_cicpower_1.bin')

ksim0_test, xi0_test = loadpower(Path0_test + '/image1/0.000_cicpower_1.bin')
ksim6_test, xi6_test = loadpower(Path6_test + '/image1/0.000_cicpower_1.bin')
# print(len(k), len(ksim), len(pk[0]), len(xi))

plt.figure(figsize=(6,8),dpi=150)
plt.subplot(2, 1, 1)
plt.loglog(k, pk0, '--', linewidth=5, color=[.7, .7, .7])
plt.loglog(ksim0, xi0[:, 4], ksim6, xi6[:, 4], ksim0_test, xi0_test[:, 4], ksim6_test, xi6_test[:, 4])
plt.grid(True, which='both', linewidth=0.3)
plt.ylim([1e-4, 1e5])
plt.xlabel('$k\\,[h/{ Mpc}]$')
plt.ylabel('$P(k)$')
plt.legend(['CAMB nonlinear','$P_{000,baseline}(k)$','$P_{006}(k,baseline)$','$P_{000,test}(k)$','$P_{006,test}(k)$' ], loc='lower left')
plt.title('$z=0.0\\,\\,\\,\\,power_{mater}$')
plt.subplot(2, 1, 2)
plt.loglog(ksim0, np.abs(xi0_test[:, 4] / xi0[:, 4]-1))
plt.loglog(ksim6, np.abs(xi6_test[:, 4] / xi6[:, 4]-1))
plt.ylim([1e-5, 1e-1])
plt.grid(True, which='both', linewidth=0.5)
plt.xlabel('$k\\,[h/{ Mpc}]$')
plt.ylabel('$/delta P(k)$')

ax = plt.gca()
ax.spines['top'].set_visible(True)
ax.spines['right'].set_visible(True)

plt.savefig(fig_path+'/power.jpg', format='jpeg')

plt.show()
# %%密度场误差

def loadfield3d(fn, ng=None):
    """
    加载3D场数据
    
    Parameters:
        fn: filename
        ng: 网格数 (可选)
        
    Returns:
        a: 3D field
        n: 网格尺寸
    """
    fid = open(fn, 'rb')
    p1 = np.fromfile(fid, dtype=np.float32)
    fid.close()
    
    if ng is None:
        n = round(len(p1) ** (1/3))
        a = np.reshape(p1, (n, n, n))
    else:
        a = np.reshape(p1, (ng, ng, ng+2))
        a = a[:ng, :ng, :ng]
        n = ng
    return a, n
def compare_field(sim):
    for im in range(1,9):
        print(f'image{im}')
        delta_baseline, n = loadfield3d(baseline_path+sim + f'/image{im}/0.000_delta_c_{im}.bin')
        delta_test, n = loadfield3d(result_path+sim + f'/image{im}/0.000_delta_c_{im}.bin')
        #对比两个密度场差值
        dd=abs(delta_baseline-delta_test)
        max=np.max(dd)
        mean=np.mean(dd)
        std=np.std(dd)
        print(f' baseline-test:\n    max:{max:.6f} mean:{mean:.6f} std:{std:.6f}')
        if mean>3 or std > 15: # 实际基准测试中 mean ~= 0.3 std ~= 2.5
            print('!'*20+f'\n image{im} too large error\n'+'!'*20)
        plt.imshow(np.mean(dd[:, :, :], axis=2).reshape(n, n).T, cmap='gray')
        plt.colorbar()
        plt.clim(-0.1, .1)
        plt.title(f'baseline-test image{im}')
        plt.savefig(fig_path +sim+ '_d-d.jpg', format='jpeg')
        plt.show()

compare_field(sim0)
compare_field(sim6)


# %% 密度场视觉禁言
im=8
sim=sim0
delta_baseline, n = loadfield3d(baseline_path+sim + f'/image{im}/0.000_delta_c_{im}.bin')
delta_test, n = loadfield3d(result_path+sim + f'/image{im}/0.000_delta_c_{im}.bin')
plt.imshow(np.mean(delta_baseline[:, :, :], axis=2).reshape(n, n).T, cmap='gray')
plt.colorbar()
plt.clim(-1,3)
plt.title(f'baseline image{im}')
plt.savefig(fig_path +sim+ '_d.jpg', format='jpeg')
plt.show()
plt.imshow(np.mean(delta_test[:, :, :], axis=2).reshape(n, n).T, cmap='gray')
plt.colorbar()
plt.clim(-1, 3)
plt.title(f'test image{im}')
plt.savefig(fig_path +sim+ '_d.jpg', format='jpeg')
plt.show()

# %%
