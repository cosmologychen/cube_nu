import camb,sys,os,re,time
import numpy as np
import scipy.integrate as integrate
from scipy.interpolate import interp1d

test = 0

n=50
ni=2

def match_para(para):
    # 读取Fortran文件
    file_path = './parameters.f90'
    with open(file_path, 'r') as file:
        content = file.read()

    # 使用正则表达式匹配模式
    pattern = r'parameter\s*::\s*%s\s*=\s*([^\s]+)'%para
    match = re.search(pattern, content)

    if match:
        variable_value = match.group(1).strip()
        return float(variable_value)
    else:
        print(f"Pattern({para}) not found in the file.")
        sys.exit()

def get_z(n): #get array a
    z = [0]*n
    a = [0]*n
    z[0] = z_max
    a = np.linspace(1/(1+z[0]),1,n+1)
    z=1/a-1
    return z[::-1]
    
def interp_pk(z,s,k_new):
    Pk0 = np.loadtxt('%s/IC/IC_P%s_rescaled_z%.4f.txt'%(nupath,s,z))[:,1]
    interp_00 = interp1d(k,Pk0, kind='linear', bounds_error=False, fill_value="extrapolate")
    new_Pk0 =interp_00(k_new)
    new_Pk0[new_Pk0 <= 0] = 1e-10
    return new_Pk0

def get_f_nr(z):
    y = Mass_nu/3/(sigma_nu*N_eff*K_b*T_gama*(1+z))
    Fy=integrate.quad(lambda u: (u**2*np.sqrt(u**2+y**2))/(1+np.exp(u)), 0,300)[0]
    f1=integrate.quad(lambda u: (u**2)/((1+np.exp(u))*np.sqrt(u**2+y**2)), 0,300)[0]
    f_nr = y**2*f1/Fy#*f_nu
    return f_nr
    
def get_Pk_nonlin_CDM(n):
    Pk_nonlin_CDM = [0]*len(z_nonlin)
    Pk_nonlin_nu = [0]*len(z_nonlin)
    nz = int(np.ceil(len(z_nonlin)/n_PK))
    for i in range(0,nz):
        z_n=z_nonlin[i*n_PK:(i+1)*n_PK]
        pars = camb.CAMBparams()
        pars.set_cosmology(H0=H0, ombh2=ombh2, omch2=omch2, omk=omk, neutrino_hierarchy=neutrino_hierarchy, num_massive_neutrinos=num_massive_neutrinos, mnu=mnu, nnu=nnu, standard_neutrino_neff=standard_neutrino_neff,tau=0.0611,YHe=0.245383,Alens=1.0026)
        pars.omch2=omch2-pars.omnuh2
        pars.InitPower.set_params(As=As,ns=ns)
        print('i = %d/%d;z_max = %.3f;len = %d'%(i+1,nz,z_n[-1],len(z_n)))
        pars.set_matter_power(redshifts=z_n, kmax=k_ic_max, nonlinear=True)
        pars.NonLinear = camb.model.NonLinear_both
        pars.NonLinearModel.set_params(halofit_version='takahashi')
        # print(pars)
        results = camb.get_results(pars)
        kh_nonlin, z_nonlin[i*n_PK:i*n_PK+len(z_n)],Pk_nonlin_CDM[i*n_PK:i*n_PK+len(z_n)]= results.get_matter_power_spectrum(minkh=k_ic_min, maxkh=k_ic_max, npoints = npbin,var1='delta_nonu',var2='delta_nonu')
        _, _                                        ,Pk_nonlin_nu[i*n_PK:i*n_PK+len(z_n)]= results.get_matter_power_spectrum(minkh=k_ic_min, maxkh=k_ic_max, npoints = npbin,var1='delta_nu'  ,var2='delta_nu')
        pars = 0
    return Pk_nonlin_CDM[::-1],Pk_nonlin_nu[::-1],kh_nonlin,z_nonlin[::-1]

#get cube parm
if (1):
    H0=match_para('h0')*100
    omega_bar=match_para('omega_bar')
    omega_cdm=match_para('omega_cdm')
    omega_r=match_para('omega_r')
    omk=0.0

    file_path = './parameters.f90'
    with open(file_path, 'r') as file:
        content = file.read()

    # 使用正则表达式匹配模式
    pattern = r'parameter\s*::\s*%s\s*=*([^\s]+)'%'m_nu'
    match = re.search(pattern, content)

    if match:
        m_nus = match.group(1).strip()
        str_cleaned = m_nus.replace("[", "").replace("]", "").replace(",", " ")
        m_nus = np.array(str_cleaned.split(), dtype=float)
    else:
        print("Pattern (m_nu) not found in the file.")
        sys.exit()

    mnu=m_nus.sum()
    if (mnu>0) :
        # if (np.all(m_nus == m_nus[0])):
        #     neutrino_hierarchy = 'degenerate'
        #     num_massive_neutrinos=3
        # elif (m_nus[0] < m_nus[2]):
        #     neutrino_hierarchy = 'normal'
        #     num_massive_neutrinos=3
        # elif (m_nus[0] > m_nus[2]):
        #     neutrino_hierarchy = 'inverted'
        #     num_massive_neutrinos=2
        neutrino_hierarchy = 'degenerate'
        num_massive_neutrinos=1
    else:
        neutrino_hierarchy = 'degenerate'
        num_massive_neutrinos=0
    nnu=3.044
    standard_neutrino_neff=match_para('N_eff')
    ns=match_para('n_s')
    As=match_para('A_s')
    z_max = np.loadtxt('./z_checkpoint.txt')[0]

    pi=np.pi
    ratio_cs=match_para('ratio_cs')
    ratio_sf=np.array([1,2,4,6,8,12,16])
    ng = match_para('ng')
    box = match_para('box')
    nns = match_para('nns')
    nnt = match_para('nnt')
    nn = match_para('nn')
    ngb = match_para('ngb')
    istep_max = match_para('istep_max')
    nc =  ng/nnt
    ng_global = ng*nn
    k_smooth = max(10.0,ng_global*pi/box)
    nfg = min(np.ceil((k_smooth*box/pi/nn) / nc)*nc,ng*1.)
    nfg_global=nfg*nn
    npbin=int(nfg_global/2*np.sqrt(3.))+22
    kmin = 2*pi/box
    kmax = pi/box*nfg_global*np.sqrt(3)
    ombh2 = omega_bar*(H0/100)**2
    omch2 = omega_cdm*(H0/100)**2

    sigma_nu = (4/11)**(1/3)
    N_nu = 3
    N_eff = standard_neutrino_neff
    K_b = 8.617342e-5
    T_gama = 2.7255
    Mass_nu = mnu
    omega_m = omega_bar+omega_cdm
    omega_cb = omega_bar+omega_cdm-Mass_nu/93.14/((H0/100)**2)
    omega_l = 1-omega_bar-omega_cdm-omega_r
    f_nu = Mass_nu/93.14/((H0/100)**2)/omega_m
    try:
        a_nu = match_para('a_nu')
    except:
        a_nu = 1./(595./5.47*(Mass_nu/3)/0.1+0.01) 
    z_nu_i = 1/a_nu-1
    omHsq = 2/3*np.sqrt(1/omega_m)/H0


    file_path = './parameters.f90'
    with open(file_path, 'r') as file:
        content = file.read()

    # 使用正则表达式匹配模式
    pattern = r'parameter\s*::\s*opath\s*=\s*([^\s]+)'
    match = re.search(pattern, content)
    opath = os.path.expanduser(match.group(1).strip()[1:-1])
    print(f"Variable value of opath : {opath}")
    nupath =  opath+"neutrinos"
    print(f"Variable value of nupath : {nupath}")
    print('\n'+('+'*40+'\n')*2)
    print('Cosmology  Paras:\n\n   omega_r:   %.6f\n   omega_b:   %.6f\n   omega_c:   %.6f\n   omega_l:   %.6f\n   mass_nu:   %.3f              eV\n  mass_nus:   %.3f %.3f %.3f %s  eV\n      f_nu:   %.6f\n\n'%(omega_r,omega_bar,omega_cdm,omega_l,Mass_nu,m_nus[0],m_nus[1],m_nus[2],neutrino_hierarchy,f_nu))
    print('Simulation Paras:\n\n     nfg:   %d\n     npbin:   %d\n     kmin:   %.3f\n     kmax:   %.3f\n\n\n'%(nfg,npbin,kmin,kmax))

    kh_nonlin = np.concatenate((np.array([1.0, 1.4142135623730951, 1.7320508075688772, 2.0, 
                                        2.23606797749979, 2.449489742783178, 2.8284271247461903, 
                                        3.0, 3.1622776601683795, 3.3166247903554, 3.4641016151377544, 
                                        3.605551275463989, 3.7416573867739413, 4.0, 4.123105625617661, 
                                        4.242640687119285, 4.358898943540674, 4.47213595499958, 
                                        4.58257569495584, 4.69041575982343, 4.898979485566356, 
                                        5.0, 5.0990195135927845, 5.196152422706632, 5.385164807134504, 
                                        5.477225575051661]) * 2 * np.pi / box
                                , np.exp(np.linspace(np.log(6*kmin), np.log(kmax),  npbin - 26))))


# print(kh_nonlin.shape)
# exit()



#mkdir
try:
    os.system('mkdir -p '+nupath+'/TF')
    os.system('mkdir -p '+nupath+'/Pk_m')
    os.system('mkdir -p '+nupath+'/Pk_nu')
    os.system('mkdir -p '+nupath+'/Pk_nus')
    os.system('mkdir -p %s/IC'%(nupath))
except:
    None
    
#set n_z
if test:
    n=n*ni
else:
    n=n
n_PK = 150 # max num of Pk(z) at once
k_ic_min = 1e-4
k_ic_max=max(1e2,kmax*1.5)

z_nonlin=get_z(n)



#get Expansion History
par = camb.CAMBparams()
par.set_cosmology(H0=H0, ombh2=ombh2, omch2=omch2, omk=omk, neutrino_hierarchy=neutrino_hierarchy, num_massive_neutrinos= num_massive_neutrinos, mnu= mnu, nnu= nnu, standard_neutrino_neff= standard_neutrino_neff,tau=0.0611,YHe=0.245383,Alens=1.0026)
par.omch2=omch2-par.omnuh2
par.InitPower.set_params(As=As,ns=ns)
par.set_matter_power(redshifts=[z_max], kmax=k_ic_max, nonlinear=True)
par.NonLinear = camb.model.NonLinear_both
par.NonLinearModel.set_params(halofit_version='takahashi')
print(par)
result = camb.get_results(par)
kh_ic, _,Pk_cb_ic= result.get_matter_power_spectrum(minkh=k_ic_min, maxkh=k_ic_max, npoints = npbin,var1='delta_nonu',var2='delta_nonu')
Pk_cb_ic = Pk_cb_ic[0]
par.set_matter_power(redshifts=[1/a_nu-1], kmax=k_ic_max, nonlinear=True)
result = camb.get_results(par)
# if (neutrino_hierarchy == 'degenerate'):
#     print('neutrino_hierarchy = degenerate')
#     kh_ic, _,Pk_nu_ic= result.get_matter_power_spectrum(minkh=k_ic_min, maxkh=k_ic_max, npoints = npbin,var1='delta_nu'  ,var2='delta_nu')
#     Pk_nu_ic = [Pk_nu_ic[0],Pk_nu_ic[0],Pk_nu_ic[0]]#list(np.zeros(npbin)),list(np.zeros(npbin))]
# else :
if(1):
    Pk_nu_ic=[0,0,0]
    i=0
    for m_nu_i in m_nus:
        if (m_nu_i >0):
            print(i,'m_nu = ',m_nu_i)
            par = camb.CAMBparams()
            par.set_cosmology(H0=H0, ombh2=ombh2, omch2=omch2, omk=omk, neutrino_hierarchy='degenerate', num_massive_neutrinos=1, mnu=m_nu_i, nnu=nnu, standard_neutrino_neff=standard_neutrino_neff,tau=0.0611,YHe=0.245383,Alens=1.0026)
            par.omch2=omch2-par.omnuh2
            par.InitPower.set_params(As=As,ns=ns)
            par.set_matter_power(redshifts=[z_max], kmax=k_ic_max, nonlinear=True)
            par.NonLinear = camb.model.NonLinear_both
            par.NonLinearModel.set_params(halofit_version='takahashi')
            print(par)
            result1 = camb.get_results(par)
            kh_ic, _,Pk_nu_ic_i= result1.get_matter_power_spectrum(minkh=k_ic_min, maxkh=k_ic_max, npoints = npbin,var1='delta_nu'  ,var2='delta_nu')
            Pk_nu_ic[i]=Pk_nu_ic_i[0]
        else:
            print(i,'m_nu = 0')
            Pk_nu_ic[i]=list(np.zeros(npbin))
        i+=1
def Ha(a):
    return result.hubble_parameter(1/a-1)
def Hz(z):
    return result.hubble_parameter(z)
    
def taua(a):
    z=1/a-1
    return integrate.quad(lambda z0: 1/Hz(z0), z,10000)[0]


T1 = time.time()

print('calculating Expansion History')
n_a = int(istep_max)
dt0 = 5e-4
t = -np.arange(n_a)*dt0
H_ex = np.ones(n_a)*H0
a_ex = np.ones(n_a)
a_ex[0] = 1
for i in range(n_a-1):
    Hai = Ha(a_ex[i])
    H_ex[i] = Hai
    a_ex[i+1] = -omHsq*Hai*a_ex[i]**3 * dt0 +a_ex[i]
    dz = 1/a_ex[i]-1/a_ex[i+1]

    if (a_ex[i+1]<=1./301 or np.isnan(a_ex[i+1])):
        break
i_end = i
T2 = time.time()

H_ex[:] = np.ones(H_ex.shape)*Ha(a_nu)
tau =  result.conformal_time(1/a_ex-1)/299792.458
os.system('rm '+nupath+'/*.txt')
np.savetxt(nupath+'/s_a_tau_H.txt',np.array([t,a_ex,tau,H_ex]))
print("EH:\n\n      time:   %.2f seconds\n      step:   %d\n      save:   '%s'\n\n\n"%(T2-T1,i_end,nupath+'/s_a_tau_H.txt'))


print('get Pk')
n=int(z_nonlin.shape[0])
Pk_nl = [0]*n
Pk_nu = [0]*n
#get pk from camb
Pk_cdm_cambs,Pk_nu_cambs,kh_nl,z_nonlin = get_Pk_nonlin_CDM(n)
Pk_nu_ic_k = interp1d(kh_ic,Pk_nu_ic, kind='linear', bounds_error=False, fill_value="extrapolate")
Pk_nu_ic = Pk_nu_ic_k(kh_nonlin)
print(kh_nl.max(),kh_nl.min())
print(kh_nonlin.max(),kh_nonlin.min())
# exit()

for i in range(n):
    # 找到Pk是否有小于0的值
    if np.any(Pk_cdm_cambs[i]<0):
        print('Pk_cdm_cambs[i]<0')
        print(Pk_cdm_cambs[i])
        print('z_nonlin[i] = ',z_nonlin[i])
        exit()
    
    Pk_cdm_z = interp1d(kh_nl,Pk_cdm_cambs[i], kind='linear', bounds_error=False, fill_value="extrapolate")
    Pk_nl[i] = Pk_cdm_z(kh_nonlin)
    Pk_nu_z = interp1d(kh_nl,Pk_nu_cambs[i], kind='linear', bounds_error=False, fill_value="extrapolate")
    Pk_nu[i] = Pk_nu_z(kh_nonlin)

#write Pk to nupath
z_powerpoint=open(nupath+'/z_powerpoint.txt', 'w')
z_values=open(nupath+'/z_values.txt', 'w')
k_values=open(nupath+'/k_values.txt', 'w')

print('write Pk to nupath')
np.savetxt(nupath+'/IC/Pcb_ic.txt',np.array([kh_ic,Pk_cb_ic]).T)
np.savetxt(nupath+'/IC/Pnu_ic.txt',np.array(Pk_nu_ic))
for i in range(len(kh_nonlin)):
    k_values.write('%3.12f\n'%kh_nonlin[i])

if test:
    for i in range(n):
        z_values.write('%3.4f\n'%z_nonlin[i])
        
        if i%ni==0:
            z_powerpoint.write('%3.4f\n'%z_nonlin[i])
        np.savetxt(nupath+'/Pk_cb_%3.4f.txt'%z_nonlin[i],Pk_nl[i])
        f_nr = get_f_nr(z_nonlin[i])
        np.savetxt(nupath+'/Pk_nu_%3.4f.txt'%z_nonlin[i],Pk_nu[i])
        np.savetxt(nupath+'/Tf_nu_%3.4f.txt'%z_nonlin[i],((1-f_nu)*np.sqrt(Pk_nl[i])+f_nr*(np.sqrt(Pk_nu[i])))/np.sqrt(Pk_nl[i]))
else:
    z_str = ''
    if (Mass_nu > 0):
        for i in range(n):
            
            np.savetxt(nupath+'/Pk_cb_%3.4f.txt'%z_nonlin[i],Pk_nl[i])

            # 找到Pk是否有小于0的值
            if np.any(Pk_nl[i]<0):
                print('Pk_cb has negative values')
                print(Pk_nl[i])
                print(z_nonlin[i])
                print(i)
                print(n)
                print(nupath)
                print(nupath+'/Pk_cb_%3.4f.txt'%z_nonlin[i])
                exit()
            f_nr = get_f_nr(z_nonlin[i])
            np.savetxt(nupath+'/Pk_nu_%3.4f.txt'%z_nonlin[i],Pk_nu[i][:])
            # 找到Pk是否有小于0的值
            if np.any(Pk_nu[i]<0):
                print('Pk_nu has negative values')
                print(Pk_nu[i])
                print(z_nonlin[i])
                print(i)
                print(n)
                print(nupath)
                print(nupath+'/Pk_nu_%3.4f.txt'%z_nonlin[i])
                exit()
            np.savetxt(nupath+'/Tf_nu_%3.4f.txt'%z_nonlin[i],((1-f_nu)*np.sqrt(Pk_nl[i])+f_nr*(np.sqrt(Pk_nu[i])))/np.sqrt(Pk_nl[i]))
            # print(nupath+'/Tf_nu_%3.4f.txt'%z_nonlin[i])
            z_str+='%3.4f\n'%z_nonlin[i]
        # for i in z_nonlin:
        # print(z_str)
        z_powerpoint.write(z_str)
    else:
        print('save '+nupath+'/Pk_cb_0.0000.txt'%z_nonlin[-1])
        np.savetxt(nupath+'/Pk_cb_0.0000.txt'%z_nonlin[-1],Pk_nl[-1])


z_powerpoint.close()
z_values.close()

print('******************** Pk_init done ********************\n\n\n')
