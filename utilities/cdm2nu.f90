
program cdm2nu
  use omp_lib
  use parameters
  use pencil_fft
  implicit none

  integer,parameter :: Nmnus = 3
  character(300),parameter :: top_path = '/home/cossim/cube_nu/output/2400_512_1_'
  character(10),parameter,dimension(Nmnus) :: Mnus=&
  ['0.05' &
  ,'0.1'  &
  ,'0.15' ]

  character(300) :: fn_tf
  integer inu,i,j,k,ig,jg,kg,cur_checkpoint
  real kr,kx(3)
  character(20) str_z,str_cpk
  real mnu,fnu
  real kh_lin(npbin),kh_lin_log(npbin),tf_F(npbin),tf_F_log(npbin)
  real tf1(nw_global/2+1,nw,npen)

  real(4),allocatable :: rho_c(:,:,:)[:,:,:]
  complex,allocatable :: irho_1(:,:,:)


  call geometry
  allocate(rho_c(nw,nw,nw)[nn,nn,*],irho_1(nw_global/2+1,nw,npen))
  if (head) then
    print*, 'cicpower on resolution: nw,nw_global=',nw,nw*nn
    print*, 'checkpoint at:'
    open(16,file='./z_checkpoint.txt',status='old')
    do i=1,nmax_redshift
      read(16,end=71,fmt='(f8.4)') z_checkpoint(i); print*, z_checkpoint(i)
    enddo
    71 n_checkpoint=i-1
    close(16); print*,''
  endif
  sync all
  n_checkpoint=n_checkpoint[1]; z_checkpoint(:)=z_checkpoint(:)[1]
  call create_penfft_plan

  ! do cur_checkpoint= 49,49
  do cur_checkpoint= n_checkpoint,n_checkpoint
  ! do cur_checkpoint= 1,n_checkpoint
    if (head) print*, ''
    if (head) print*, '==========================================='
    if (head) print*, '==========================================='
    if (head) print*, 'Start analyzing redshift ',z2str(z_checkpoint(cur_checkpoint))
    !print*,output_name('info')
    open(11,file=output_name('info'),access='stream'); read(11) sim; close(11)
    if (sim%izipx/=izipx .or. sim%izipv/=izipv) error stop 'zip format incompatable'
    if (head) print*, 'nplocal,npglobal =',sim%nplocal,sim%npglobal
    write(str_cpk,'(i0)') sim%calculate_PK

    if (head) print*,'Write delta_c into',output_name('delta_c')
    open(11,file=output_name('delta_c'),status='old',access='stream')
    read(11) rho1
    close(11)

    call pencil_fft_forward
    irho_1 = cxyz

    do inu = 1,Nmnus
      print*, 'Mnu = ', Mnus(inu)
      write(str_z,'(f8.4)') z_powerpoint(sim%cur_powerpoint)
      fn_tf = trim(top_path)//trim(Mnus(inu))//'_'//trim(str_cpk)//'/neutrinos/tf/Tf_nu_'//trim(adjustl(str_z))//'.txt'
      ! print*, trim(top_path),trim(Mnus(inu))
      print*,'    Reading Tf_nu form',fn_tf
      open(10,file=trim(fn_tf),status='old',access='stream'); read(10) tf_F; close(10)
      read(Mnus(inu), *) mnu
      print*, '    mnu', mnu
      fnu = mnu/93.14/(h0**2)/omega_m
      print*, '    fnu',fnu
      print*, '    tf_F',tf_F(1),tf_F(2),tf_F(npbin)
      tf_F = (tf_F-(1-fnu))/fnu
      print*, '    tf_F',tf_F(1),tf_F(2),tf_F(npbin)
      kh_lin = -1
      open(13,file=nupath//'k_values.txt',status='old')
      do i=1,npbin
      read(13,end=75,fmt='(f15.12)') kh_lin(i)
      enddo
      75 close(13)
      kh_lin_log = log(kh_lin)
      print*,'    k',kh_lin(1),kh_lin(2),kh_lin(npbin),icx,icy,icz,npen,nw,nyquist
      tf_F_log = log(tf_F)

      !$omp paralleldo default(shared) schedule(dynamic)&
      !$omp& private(kg,jg,ig,i,j,k,kx,kr)
      do k=1,npen
      do j=1,nw
      do i=1,nyquist+1
        kg=(nn*(icz-1)+icy-1)*npen+k
        jg=(icx-1)*nw+j
        ig=i
        kx=(mod([ig,jg,kg]+nyquist-1,nw_global)-nyquist)*(2*pi)/box
        if (ig==1.and.jg==1.and.kg==1) cycle ! zero frequency
        kr=sqrt(kx(1)**2+kx(2)**2+kx(3)**2)
        tf1(i,j,k) = interp_tf_F(kh_lin_log,kr,tf_F_log)
      enddo
      enddo
      enddo
      !$omp endparalleldo

      cxyz=irho_1*tf1
      call pencil_fft_backward
      rho_c = rho1
      
      print *,'    save delta nu into ',output_name('delta_nu_'//trim(Mnus(inu)))
      open(11,file=output_name('delta_nu_'//trim(Mnus(inu))),status='replace',access='stream')
      write(11) rho_c
      close(11); sync all
    enddo
  enddo

  contains

  real function interp_tf_F(k_lin_log,k_need,tf_F_log)
      implicit none
      real, intent(in) :: k_lin_log(npbin),k_need,tf_F_log(npbin)
      real a1,a2,k2,b2,c2
      integer(4) i_mid,i1,i2,i3
      real k_need_log
  
      k_need_log = log(k_need)
      
      i1=1; i2=npbin
      do while (i2-i1>1)
      i_mid=(i1+i2)/2
      if (k_need_log>k_lin_log(i_mid)) then
          i1=i_mid
      else
          i2=i_mid
      endif
      enddo
  
      if (k_need_log-k_lin_log(i1) < k_lin_log(i2)-k_need_log) then
          i3 = i2
          i2 = i1
          i1 = i1-1
      else 
          i3 = i2+1
      endif
  
      if (i3 > npbin) then
          i1 = i1-1
          i2 = i2-1
          i3 = i3-1
      endif
  
      a1 = (tf_F_log(i1)-tf_F_log(i2))/(k_lin_log(i1)-k_lin_log(i2))
      a2 = (tf_F_log(i2)-tf_F_log(i3))/(k_lin_log(i2)-k_lin_log(i3))
      k2 = (a1-a2)/(k_lin_log(i1)-k_lin_log(i3))
      b2 = a1-k2*(k_lin_log(i1)+k_lin_log(i2))
      c2 = tf_F_log(i1)-k2*k_lin_log(i1)**2-b2*k_lin_log(i1)
      interp_tf_F = exp(k2*k_need_log**2+b2*k_need_log+c2)
  
      ! print*,'interp:++++++++++++++'
      ! print*,i1,i2,i3
      ! print*,kh_lin(i1),kh_lin(i2),kh_lin(i3),k_need
      ! print*,tf_F(i1),tf_F(i2),tf_F(i3),interp_tf_F
      ! print*,'interp,done'
  
  endfunction 
  
end program
