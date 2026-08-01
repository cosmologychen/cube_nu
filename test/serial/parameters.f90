module parameters
   implicit none

   ! output directory
   character(*),parameter :: opath='../../serial/'

   ! zip parameters
   integer,parameter :: ndim=3
#ifdef ZIPX
   integer,parameter :: izipx=2 ! size to store xp as
#else
   integer,parameter :: izipx=-1 ! size to store xp as
#endif
   integer,parameter :: izipv=2 ! size to store vp as
   integer, parameter :: izipi = 8 ! if pids are on, size to store as

   ! cell resolution parameters
   integer,parameter :: nn=1
   real,parameter :: box=600
   integer,parameter :: ncore=1
   integer,parameter :: nteam=1
   integer,parameter :: nnest=1
   integer(8),parameter :: ng=128

   integer,parameter :: ratio_cs=4
   integer,dimension(7),parameter :: ratio_sf=[1,2,4,6,8,12,16]
   integer,dimension(7),parameter :: ratio_cf=ratio_cs*ratio_sf ! 8,16,24,32
   integer,parameter :: nnt=4 ! number of tiles /image/dim
   integer,parameter :: nns=4 ! number of subtiles /tile/dim

   ! particle resolution parameters
   integer(8),parameter :: np_nc=ratio_cs
   integer,parameter :: nic=2   ! refined resolution for IC
   integer,parameter :: ngic=ng*nic
   logical,parameter :: body_centered_cubic=.false.

   ! force resolution parameters
   real,parameter :: apm1c=3.5            ! PM1's softening
   real,parameter :: apm1=apm1c*ratio_cs  ! PM2's range = 14
   real,parameter :: apm2=3.5             ! PM2's softening = 3.5
   real,dimension(7),parameter :: apm2f=apm2*ratio_sf  ! PM3's range =(3.5),7,14,21,28
   real,parameter :: apm3f=3.5            ! PM3's softening = 3.5
   real,dimension(7),parameter :: apm3=3.5/ratio_sf    ! PP's range
   real,dimension(7),parameter :: appr=apm3
   real,parameter :: dt_refine=1 ! rescale dt
   real,parameter :: app=0.06 ! target=0.05 !apm3f/ratio_sf ! 0.875 pp softening: apm3f/ratio_sf
   integer,parameter :: nrange=1 ! 1 actual PP computing range
   integer,parameter :: n_neighbor=((1+2*nrange)**3-1)/2 ! 13

   ! Green's function
   integer,parameter :: n_int=2
   integer,parameter :: p=2 ! order of interpolation. 1=CIC, 2=TSC
   real,parameter :: alpha=4./3
   real,dimension(4),parameter :: weight = [(alpha-1)/4,-alpha/2,alpha/2,(1-alpha)/4]

   ! derived parameters
   integer(8),parameter :: nc=ng/ratio_cs ! 192 coarse cells /image/dim, must be integer
   integer(8),parameter :: nt=nc/nnt ! 96
   integer,parameter :: ntt=nt/nns ! 24
   integer(8),parameter :: ngp=ng/nnt ! 384 physical tile
   integer,parameter :: ngb=16 ! tile buffer, floor(apm1)+2
   integer(8),parameter :: ngt=ngp+2*ngb ! 416
   integer(8),parameter :: ng_global=ng*nn
   integer(8),parameter :: nc_global=nc*nn

   integer(8),dimension(7),parameter :: nfp=(ngp/nns)*ratio_sf ! 192*[1,2,4,6,8,12,16]
   !integer(8),parameter :: nf=nfp
   integer,dimension(7),parameter :: nfb=ratio_cf
   integer(8),dimension(7),parameter :: nft=nfp+2*nfb

   ! ngrid /image/dim for pencil-fft
# ifdef FFTFINE
   integer(8),parameter :: nw=ngic ! standard grid fft, for IC, dsp, convert_zip_to_xv
# else
   integer(8),parameter :: nw=nc ! coarse grid fft, for main code
# endif
   integer(8),parameter :: npen=nw/nn ! ng /dim in shorter side of the pencil, for pencil decomposition
   integer(8),parameter :: nw_global=nw*nn
   integer(8),parameter :: nyquist=nw_global/2

   integer(8),parameter :: ncb=4 ! nc in buffer /dim, single side; 6 by default
   integer(8),parameter :: nte=nt+2*ncb ! extended nc


   real,parameter :: image_buffer=2.0
   real,parameter :: tile_buffer=8.0
   real,parameter :: vbuf=0.9

   ! derived parameters
   integer(8),parameter :: nvbin=int(2,8)**(8*izipv)
   integer(8),parameter :: ishift=-(int(2,8)**(izipx*8-1))
   real(8),parameter :: rshift=0.5-ishift

   real(8),parameter :: pi=4*atan(1d0)

   ! cosmological parameters
   real,parameter :: s8=0.8091 ! use -Dsigma_8 in initial_conditions
   integer,parameter :: zdim=2 ! the dimension being the redshift direction

   ! background parameters
   real, parameter :: h0 =0.636
   !nu
   real(8),dimension(3),parameter :: m_nu=[0.0,0.0,0.0] ! mass_nus/eV
   real(8),parameter :: Mass_nu=sum(m_nu) ! Mass_nu/eV

   ! real,dimension(3), parameter :: O_nu = m_nu/93.14/(h0**2)
   real, parameter :: omega_nu = Mass_nu/93.14/(h0**2)
   real, parameter :: omega_cdm = 0.29729389546402873 - omega_nu! cdm energy
   real, parameter :: omega_bar = 0.054907835924211855 ! baryon energy, goes into cdm
   real, parameter :: omega_mhd = 0.0 ! mhd energy, evolved separately
   real, parameter :: omega_r = 5.046734693877551e-05


   real, parameter :: omega_m = omega_cdm+omega_bar+omega_mhd+omega_nu ! total matter
   real, parameter :: omega_l = 1-omega_m-omega_r
   ! real, parameter :: omega_l = 1-omega_m
   real, parameter :: wde = -1 ! de equation of state

   ! initial conditions
   real,parameter :: f_nl=0
   real,parameter :: g_nl=0
   real,parameter :: n_s=0.9644
   real,parameter :: A_s=2.087e-09
   real,parameter :: k_o=0.05/h0

   integer(8),parameter :: istep_max=100000
   real,parameter :: ra_max=0.1
   real(8),parameter :: v_resolution=2.1/(int(2,8)**(izipv*8))

#ifdef ZIPX
   real(8),parameter :: x_resolution=1.0/(int(2,8)**(izipx*8))
#endif
   !real(8),parameter :: vdisp_boost=1.0
   real(8),parameter :: vrel_boost=2.5
   real,parameter :: b_link=0.20 ! linking length for FoF
   real,parameter :: np_halo_min=10 ! minimum number of particles to be a halo

   !! MPI image variables !!
   integer(8) image,rank,icx,icy,icz,inx,iny,inz,ipx,ipy,ipz
   integer(8) m1,m2,m3,m
   logical head
   ! checkpoint variables
   integer(8),parameter :: nmax_redshift=1000
   integer(8) n_checkpoint[*],n_halofind[*]
   real z_checkpoint(nmax_redshift)[*],z_halofind(nmax_redshift)[*]
   logical checkpoint_step[*],halofind_step[*],final_step[*]

   ! timing variables
   integer(4) istep,t1,t2,tt1,tt2,ttt1,ttt2,t_start,t_end,t_rate,tictoc(2,0:100),tnu1,tnu2
   integer(4) tp1,tp2,ttp1,ttp2,tpr,cc
   real tps(5,100)
   real tcat(100,0:10000)



   !nu
   character(*),parameter :: nupath=trim(adjustl(opath))//'neutrinos/' !path for neutrino
   real z_powerpoint(nmax_redshift)[*] ! calculate Pk
   integer(8) n_powerpoint[*]
   logical power_step[*]

   real(8),parameter :: omhsq0=2./3.*sqrt(1/omega_m)/h0/100
   real(8),parameter :: f_nu=omega_nu/omega_m ! neutrino fraction
   real(8),dimension(3),parameter :: f_nus=m_nu/max(Mass_nu,0.0000001) ! neutrino fraction
   real(8),parameter :: N_eff=3.044
   real(8),parameter :: T_gama = 2.7255 ! neutrino tempture/eV today
   real(8),parameter :: T_nu0=0.00016764 ! photon tempture/eV today
   real(8),parameter :: C=299792.458*h0 ! speed of light
   real(8),parameter :: sigma_nu = 0.71649 !the neutrino to photon temperature ratio today
   real(8),parameter :: k_b = 8.617342e-5 !the Boltzmann’s con-stant
   real(8),parameter :: a_nu=0.016667 ! 1./(595./5.47*(Mass_nu/3)/0.1+0.01) ! nu is matter in a_nu

   real(8),parameter :: k_smooth=max(10.0,ng_global*pi/box) !smooth the tf in k >k_smooth
   integer(8),parameter :: tf_smooth=100 !smooth the tf in k >tf_smooth*k_fs
   integer(8),parameter :: calculate_PK = 2
   integer(8),parameter :: nfg=min(ceiling((k_smooth*box/pi/nn) / nc)*nc*1.,ng*1.)
   integer(8),parameter :: nfg_global=nfg*nn
   integer(8),parameter :: npbin=int(nfg_global/2*sqrt(3.))+22



   type sim_header
      integer(8) nplocal,npglobal
      integer(8) izipx,izipv
      integer(8) image
      integer(8) nn,nnt,nt,ncell,ncb
      integer(8) timestep
      integer(8) cur_checkpoint
      integer(8) cur_halofind
      integer(8) cur_powerpoint !nu
      integer(8) calculate_PK !nu
      integer(8) cic_iapm !nu
      real a, t, tau
      real dt_pp, dt_pm3, dt_pm2, dt_pm1, dt_vmax
      real mass_p_cdm
      real m_nu(3)!nu
      real Mass_nu!nu
      real box
      real h0
      real omega_m
      real omega_nu!nu
      real omega_l
      real s8
      real vsim2phys
      real sigma_vres
      real sigma_vi
      real z_i
      real vz_max
   endtype
   type(sim_header) sim[*]


   type type_halo_catalog_header
      integer nhalo_tot,nhalo,ninfo
      real linking_parameter
   endtype
   type type_halo_catalog_array
      real hmass ! number of particles ! 1:1
      real xv(6)
   endtype
   integer,parameter :: ninfo=7 ! number of real numbers per halo in the halo catalog
   character(4) b_link_string


contains
   subroutine print_header(s)
      type(sim_header),intent(in) :: s
      if (this_image()==1) then
         print*,'-------------------------------- CUBE info --------------------------------'
         print*,'| np local/global =',s%nplocal,s%npglobal
         print*,'| a,t,tau         =',s%a,s%t,s%tau
         print*,'| timestep        =',s%timestep
         print*,'| dt PM123,PP     =',s%dt_pm1,s%dt_pm2,s%dt_pm3,s%dt_pp
         print*,'| dt v            =',s%dt_vmax
         print*,'| cur_checkpoint  =',int(s%cur_checkpoint,2),z_checkpoint(s%cur_checkpoint)
         print*,'| cur_powerpoint  =',int(s%cur_powerpoint,2),z_powerpoint(s%cur_powerpoint)
         print*,'| cur_halofind    =',int(s%cur_halofind,2),z_halofind(s%cur_halofind)
         print*,'| mass_p          =',s%mass_p_cdm
         print*,'| mass_nu         =',s%mass_nu,s%m_nu,s%calculate_PK
         print*,'| box             =',s%box, 'Mpc/h'
         print*,'| image           =',s%image
         print*,'| nn              =',s%nn
         print*,'| nnt             =',s%nnt
         print*,'| ncb             =',s%ncb
         print*,'| izip x,v        =',int(s%izipx,1),int(s%izipv,1)
         print*,'| h_0             =',s%h0,'*100 km/s/Mpc'
         print*,'| omega_m         =',s%omega_m
         print*,'| omega_nu         =',s%omega_nu
         print*,'| omega_l         =',s%omega_l
         print*,'| sigma_8         =',s%s8
         print*,'| vsim2phys       =',s%vsim2phys, '(km/s)/(1.0)'
         print*,'| sigma_vres      =',s%sigma_vres,'(km/s)'
         print*,'| sigma_vi        =',s%sigma_vi,'(simulation unit)'
         print*,'| z_i             =',s%z_i
         print*,'| vz_max          =',s%vz_max
         print*,'------------------------------------------------------------------------------'
      endif
      sync all
   endsubroutine

   include 'basic_functions.f08'
endmodule
