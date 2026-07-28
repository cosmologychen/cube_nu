import camb
import os
import struct
import math
import numpy as np
import matplotlib.pyplot as plt


def loadsim(filename):
    """
    加载模拟参数
    """
    sim = {}
    with open(filename, 'r') as f:
        for line in f:
            if '=' in line and not line.strip().startswith('#'):
                key, value = line.strip().split('=')
                key = key.strip()
                value = value.strip()
                try:
                    if '.' in value:
                        sim[key] = float(value)
                    else:
                        sim[key] = int(value)
                except ValueError:
                    sim[key] = value
    return sim


def loadz(filename):
    """
    加载红移数据
    """
    z = np.loadtxt(filename)
    return z


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
    print(fn, n)
    return a, n


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


def loadtcat(filename):
    """
    加载时间统计数据
    
    Parameters:
        filename: 时间统计文件名
        
    Returns:
        nstep: 步数
        tcat0: 初始时间统计
        tcat: 时间统计数组
    """
    with open(filename, 'rb') as fid:
        nstep = np.fromfile(fid, dtype=np.int32, count=1)[0]
        tcat0 = np.fromfile(fid, dtype=np.float32, count=100)
        tcat = np.fromfile(fid, dtype=np.float32, count=100 * nstep)
        tcat = tcat.reshape(nstep, 100).T
    
    return nstep, tcat0, tcat


def cambpower(z, sim):
    """
    使用CAMB计算功率谱
    
    Parameters:
        z: 红移
        sim: 模拟参数字典
        
    Returns:
        k: k值数组
        pk: 功率谱
        s_8: sigma8值
    """
    h0 = sim['h0']
    h02 = h0**2
    H0 = h0 * 100
    ombh2 = (0.04630) * h02
    omch2 = (sim['omega_m'] - sim['omega_nu'] - 0.04630) * h02
    omk = 0.0
    mnu = sim['Mass_nu']
    neutrino_hierarchy = 'degenerate'
    num_massive_neutrinos = 3
    nnu = 3.044
    standard_neutrino_neff = 3.044
    kmax = 200.0

    pars = camb.CAMBparams()
    pars.set_cosmology(
        H0=h0*100, 
        ombh2=ombh2, 
        omch2=omch2, 
        omk=omk, 
        neutrino_hierarchy=neutrino_hierarchy, 
        num_massive_neutrinos=num_massive_neutrinos, 
        mnu=mnu, 
        nnu=nnu, 
        standard_neutrino_neff=standard_neutrino_neff
    )
    pars.omch2 = omch2 - pars.omnuh2
    pars.InitPower.set_params(As=2.6025e-09, ns=0.972)
    pars.set_matter_power(redshifts=[z], kmax=kmax, nonlinear=True)
    pars.NonLinear = camb.model.NonLinear_both

    results = camb.get_results(pars)
    k, zout, pk = results.get_matter_power_spectrum(
        minkh=1e-4, 
        maxkh=kmax, 
        npoints=100,
        var1='delta_nonu',
        var2='delta_nonu'
    )
    s_8 = results.get_sigma8()
    return k, pk, s_8


def plt_proj(Redshift, Path, Dir, xgrid):
    """
    绘制投影图
    
    Parameters:
        Redshift: 红移字符串
        Path: 数据路径
        Dir: 数据目录
        xgrid: 网格范围
    """
    proj_xy = loadfield2d(Path + Dir + Redshift + '_proj_xy_1.bin')
    
    plt.figure()
    plt.imshow(proj_xy.T, cmap='gray', extent=[xgrid[0], xgrid[1], xgrid[0], xgrid[1]])
    
    ax = plt.gca()
    ax.spines['top'].set_visible(True)
    ax.spines['right'].set_visible(True)
    
    plt.colorbar()
    plt.clim(-1, 2)
    plt.title('$z=' + str(Redshift) + '\\,\\,\\,\\,\\delta_c$')
    
    savename = Path + '/fig/' + 'proj_' + str(Redshift) + '_xy.jpg'
    plt.savefig(savename, format='jpeg')
    plt.show()


def plt_power(Redshift, Path, Dir, sim):
    """
    绘制功率谱
    
    Parameters:
        Redshift: 红移字符串
        Path: 数据路径
        Dir: 数据目录
        sim: 模拟参数字典
    """
    k, pk, s_8 = cambpower(float(Redshift), sim)
    ksim, xi = loadpower(Path + Dir + Redshift + '_cicpower_1.bin')
    print(len(k), len(ksim), len(pk[0]), len(xi))

    plt.figure()
    plt.subplot(2, 1, 1)
    plt.loglog(k, pk[0], '--', linewidth=2.5, color=[.7, .7, .7])
    print(max(ksim), min(ksim), len(ksim))
    plt.loglog(ksim, xi[:, 2], ksim, xi[:, 2] - xi[:, 3], ksim, xi[:, 3], '--', ksim, xi[:, 4])
    plt.grid(True, which='both', linewidth=0.3)
    plt.xlim([min(ksim), max(ksim)])
    plt.ylim([1e-4, 1e5])
    plt.xlabel('$k\\,[h/{ Mpc}]$')
    plt.ylabel('$P(k)$')
    plt.legend(['CAMB nonlinear', '$P_{ raw}(k)$', '$D^2(k)=C_1(k)/N$', '$P_{ raw}(k)-D^2(k)$', '$P_0(k)$'], loc='upper right')
    plt.title('$z=' + str(Redshift) + '\\,\\,\\,\\,power_{mater}$')
    
    ax = plt.gca()
    ax.spines['top'].set_visible(True)
    ax.spines['right'].set_visible(True)

    plt.subplot(2, 1, 2)
    plt.loglog(ksim, xi[:, 3] / xi[:, 4])
    plt.grid(True, which='both', linewidth=0.5)
    plt.xlabel('$k\\,[h/{ Mpc}]$')
    plt.ylabel('$C_2(k)$')
    
    ax = plt.gca()
    ax.spines['top'].set_visible(True)
    ax.spines['right'].set_visible(True)

    savename = Path + '/fig/power_11' + str(Redshift) + '.jpg'
    plt.savefig(savename, format='jpeg')
    plt.show()


def plt_tcat(z, nsnap_read, Path, Dir):
    """
    绘制时间统计图
    
    Parameters:
        z: 红移数组
        nsnap_read: 快照索引列表
        Path: 数据路径
        Dir: 数据目录
    """
    for isnap in range(len(nsnap_read)):
        Redshift = '{:.3f}'.format(z[int(nsnap_read[isnap])])
        nstep1, tcat0, tcat1 = loadtcat(Path + Dir + Redshift + '_tcat_1.bin')
        for i in range(1, nstep1-1):
            for j in range(100):
                if tcat1[j, i] < 0:
                    tcat1[j, i] = np.mean(tcat1[j, [i-1, i+1]])
        
        if isnap == 0:
            tcat = tcat1
            nstep = nstep1
        else:
            n = tcat1.size // 100
            tcat[:, :n-1] = tcat1[:, :n-1]

    for i in range(1, nstep-1):
        for j in range(100):
            if tcat[j, i] == 0:
                tcat[j, i] = np.mean(tcat[j, [i-1, i+1]])
                if j == 99:
                    tcat[j, i] = tcat[j, i] + 180

    istep = np.arange(1, nstep-2).T
    tcat[99, nstep-3] = tcat[99, nstep-3] + 180
    tcat[50:52, :] = tcat[50:52, :] - 1
    tcat = tcat.T
    print('compute time: %d days' % ((sum(tcat[:, 99]))/3600/24))

    plt.figure(figsize=(10, 8))

    plt.subplot(2, 1, 1)
    plt.plot(istep, tcat[:nstep-3, 42], label='$scale\\,factor\\,a+{ d}a$', linewidth=0.5)
    plt.plot(istep, tcat[:nstep-3, 50], label='$memory\\, overhead\\,\\delta_{ node}$', linewidth=0.5)
    plt.plot(istep, tcat[:nstep-3, 51], label='$memory\\,overhead\\,\\delta_{ tile}$', linewidth=0.5)
    plt.grid(True, which='both', linestyle='dotted', linewidth=0.3)
    plt.xlabel('timestep')
    plt.ylabel('')
    plt.legend(loc='upper left')
    plt.xlim(0, nstep-3)
    plt.title('memory')
    ax = plt.gca()
    ax.spines['top'].set_visible(True)
    ax.spines['right'].set_visible(True)

    plt.subplot(2, 1, 2)
    plt.plot(istep, tcat[:nstep-3, 10], label='PM1', linewidth=0.5)
    plt.plot(istep, tcat[:nstep-3, 11], label='PM2', linewidth=0.5)
    plt.plot(istep, tcat[:nstep-3, 12], label='PM3', linewidth=0.5)
    plt.plot(istep, tcat[:nstep-3, 13], label='PM1', linewidth=0.5)
    plt.plot(istep, tcat[:nstep-3, 1], label='drift', linewidth=0.5)
    plt.plot(istep, tcat[:nstep-3, 99], label='$\\sum$', linewidth=0.5)
    plt.grid(True, which='both', linestyle='dotted', linewidth=0.3)
    plt.xlabel('timestep')
    plt.ylabel('time/sec')
    plt.legend(loc='upper left')
    plt.xlim(0, nstep-3)
    plt.title('time')
    ax = plt.gca()
    ax.spines['top'].set_visible(True)
    ax.spines['right'].set_visible(True)

    plt.tight_layout()
    savename = Path + '/fig/tcat_%.3f.jpg' % z[int(max(nsnap_read))]
    plt.savefig(savename, format='jpeg')
    plt.show()


def plt_delta_field(Redshift, Path, Dir, ng, axis=2, clim=(-3, 3)):
    """
    绘制密度场
    
    Parameters:
        Redshift: 红移字符串
        Path: 数据路径
        Dir: 数据目录
        ng: 网格数
        axis: 投影轴 (0, 1, 2)
        clim: 颜色范围
    """
    delta_L, n = loadfield3d(Path + Dir + Redshift + '_delta_c_1.bin')
    
    plt.figure()
    if axis == 0:
        img = np.mean(delta_L[:, :, :], axis=0).reshape(n, n).T
    elif axis == 1:
        img = np.mean(delta_L[:, :, :], axis=1).reshape(n, n).T
    else:
        img = np.mean(delta_L[:, :, :], axis=2).reshape(n, n).T
    
    plt.imshow(img, cmap='gray')
    plt.colorbar()
    plt.clim(clim[0], clim[1])
    plt.title('$\\delta_c z=' + str(Redshift) + '$')
    plt.show()


def plt_phi_field(Redshift, Path, Dir, ng, axis=2):
    """
    绘制势能场
    
    Parameters:
        Redshift: 红移字符串
        Path: 数据路径
        Dir: 数据目录
        ng: 网格数
        axis: 投影轴 (0, 1, 2)
    """
    phi, n = loadfield3d(Path + Dir + Redshift + '_phi1_1.bin')
    
    plt.figure()
    if axis == 0:
        img = np.mean(phi[:, :, :], axis=0).reshape(n, n).T
    elif axis == 1:
        img = np.mean(phi[:, :, :], axis=1).reshape(n, n).T
    else:
        img = np.mean(phi[:, :, :], axis=2).reshape(n, n).T
    
    plt.imshow(img, cmap='gray')
    plt.colorbar()
    plt.title('$\\phi$')
    plt.show()
