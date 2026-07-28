!! module for power spectrum analysis
! uses pencil_fft however work for single image only
! cx1,cx2 can be memory optimized
! auto-power can be memory optimized
! check nexp frequently
#define linear_kbin
!#define pl2d

module powerspectrum
use omp_lib
use pencil_fft

#ifdef linear_kbin
  integer(8),parameter :: nbin=nint(nyquist*sqrt(3.))
#else
  integer(8),parameter :: nbin=floor(4*log(nyquist*sqrt(3.)/0.95)/log(2.))
#endif
#ifdef pl2d
  integer kp,kl
  integer nmode(nint(nyquist*sqrt(2.))+1,nyquist+1)
  real pow2d(nint(nyquist*sqrt(2.))+1,nyquist+1)
  real pow2drsd(nint(nyquist*sqrt(2.))+1,nyquist+1)
#endif
complex cx1(nw*nn/2+1,nw,npen),cx2(nw*nn/2+1,nw,npen)

contains

subroutine cross_power(xip,cube1,cube2)
  use omp_lib
  implicit none

  integer i,j,k,ig,jg,kg,ibin
  real kr,kx(3),sincx,sincy,sincz,sinc,rbin

  real cube1(nw,nw,nw),cube2(nw,nw,nw)
  real amp11,amp12,amp22,xi(10,nbin),xip(10,nbin)[*]
  complex cx1(nw*nn/2+1,nw,npen),cx2(nw*nn/2+1,nw,npen)

  real,parameter :: nexp=4.0 ! CIC kernel
  xi=0
  rho1=cube1
  call pencil_fft_forward
  cx1=cxyz/nw_global/nw_global/nw_global
  print*, ' fft 1 done'

  rho1=cube2
  call pencil_fft_forward
  cx2=cxyz/nw_global/nw_global/nw_global
  print*, ' fft 2 done'
  xi=0
  sync all
#ifdef pl2d
  print*, 'size of pow2d',nint(nyquist*sqrt(2.))+1,nyquist+1
  pow2d=0
  nmode=0
#endif

  do k=1,npen
  do j=1,nw
  do i=1,nyquist+1
    
    kg=(nn*(icz-1)+icy-1)*npen+k
    jg=(icx-1)*nw+j
    ig=i
    kx=mod((/ig,jg,kg/)+nyquist-1,nw_global)-nyquist
    if (ig==1.and.jg==1.and.kg==1) cycle ! zero frequency
    if ((ig==1.or.ig==nw*nn/2+1) .and. jg>nw*nn/2+1) cycle
    if ((ig==1.or.ig==nw*nn/2+1) .and. (jg==1.or.jg==nw*nn/2+1) .and. kg>nw*nn/2+1) cycle
    kr=sqrt(kx(1)**2+kx(2)**2+kx(3)**2)
    ! sincx=merge(1d0,sin(pi*kx(1)/nw_global)/(pi*kx(1)/nw_global),kx(1)==0.0)
    ! sincy=merge(1d0,sin(pi*kx(2)/nw_global)/(pi*kx(2)/nw_global),kx(2)==0.0)
    ! sincz=merge(1d0,sin(pi*kx(3)/nw_global)/(pi*kx(3)/nw_global),kx(3)==0.0)
    sinc=sincx*sincy*sincz
#   ifdef linear_kbin
      ibin=nint(kr)
#   else
      rbin=4.0/log(2.)*log(kr/0.95)
      ibin=merge(ceiling(rbin),floor(rbin),rbin<1)
#   endif
    xi(1,ibin)=xi(1,ibin)+1 ! number count
    xi(2,ibin)=xi(2,ibin)+kr ! k count
    amp11=real(cx1(i,j,k)*conjg(cx1(i,j,k)))!/(sinc**4.0)*4*pi*kr**3
    amp22=real(cx2(i,j,k)*conjg(cx2(i,j,k)))!/(sinc**4.0)*4*pi*kr**3
    amp12=real(cx1(i,j,k)*conjg(cx2(i,j,k)))!/(sinc**4.0)*4*pi*kr**3
#ifdef pl2d
    kp=nint(sqrt(kx(1)**2+kx(3)**2))+1
    kl=abs(kx(2))+1
    nmode(kp,kl)=nmode(kp,kl)+1
    pow2d(kp,kl)=pow2d(kp,kl)+amp11
    pow2drsd(kp,kl)=pow2drsd(kp,kl)+amp22
#endif
    xi(3,ibin)=xi(3,ibin)+amp11 ! auto power 1
    xi(4,ibin)=xi(4,ibin)+amp22 ! auto power 2
    xi(5,ibin)=xi(5,ibin)+amp12 ! cross power
    ! xi(6,ibin)=xi(6,ibin)+1/sinc**2.0 ! kernel 1
    ! xi(7,ibin)=xi(7,ibin)+1/sinc**4.0 ! kernel 2

  enddo
  enddo
  enddo
  sync all
#ifdef pl2d
  nmode=max(1,nmode)
  pow2d=pow2d/nmode
  pow2drsd=pow2drsd/nmode
  open(55,file=output_name('pow2d'),status='replace',access='stream')
  write(55) pow2d
  write(55) pow2drsd
  close(55)
#endif
sync all
xip=xi(:,1:)

  ! co_sum
  if (head) then
    do i=2,nn**3
      xip=xip+xip(:,:)[i]
    enddo
  endif
  sync all

  ! broadcast
  xip=xip(:,:)[1]
  sync all

  ! divide and normalize
  xip(2,:)=xip(2,:)/xip(1,:)*(2*pi)/box ! k_phy
  xip(3,:)=xip(3,:)/xip(1,:) ! Delta_LL
  xip(4,:)=xip(4,:)/xip(1,:) ! Delta_RR
  xip(5,:)=xip(5,:)/xip(1,:) ! Delta_LR ! cross power
  ! xip(6,:)=xip(6,:)/xip(1,:) ! kernel
  ! xip(7,:)=xip(7,:)/xip(1,:) ! kernel
  ! xip(8,:)=xip(5,:)/sqrt(xip(3,:)*xip(4,:)) ! r
  ! xip(9,:)=sqrt(xip(4,:)/xip(3,:)) ! b
  ! xip(10,:)=xip(8,:)**4/xip(9,:)**2 * xip(4,:) ! P_RR*r^4/b^2 reco power

  sync all
endsubroutine cross_power

subroutine auto_power(xi,cube1,n_particle,n_interp)
  use omp_lib
  implicit none
  
  integer(4) t_xi1,t_xi2,t_xi_rate

  integer i,j,k,ig,jg,kg,ibin,n_interp
  integer(8) n_particle
real alpha,kr,kx(3),sincx,sincy,sincz,sinc,rbin,C1k(3),Dk,amp11,cube1(nw,nw,nw),xi(10,0:nbin)[*]
  alpha=0
  xi=0
  rho1=cube1
  call pencil_fft_forward
  cxyz=cxyz/nw_global/nw_global/nw_global
  if (head) print*, 'check: min,max of rho_k = '
  sync all
  ! if (head) 
  print*, image,minval(real(cxyz)),maxval(real(cxyz))
  sync all
  call system_clock(t_xi1,t_xi_rate)
  do k=1,npen
  do j=1,nw
    do i=1,nyquist+1
    kg=(nn*(icz-1)+icy-1)*npen+k
    jg=(icx-1)*nw+j
    ig=i
    kx=mod([ig,jg,kg]+nyquist-1,nw_global)-nyquist
    if (ig==1.and.jg==1.and.kg==1) cycle ! zero frequency
    if ((ig==1.or.ig==nw*nn/2+1) .and. jg>nw*nn/2+1) cycle
    if ((ig==1.or.ig==nw*nn/2+1) .and. (jg==1.or.jg==nw*nn/2+1) .and. kg>nw*nn/2+1) cycle
    kr=sqrt(kx(1)**2+kx(2)**2+kx(3)**2)
    ibin=nint(kr)
    xi(1,ibin)=xi(1,ibin)+1 ! number count
    xi(2,ibin)=xi(2,ibin)+kr ! k count
    amp11=real(cxyz(i,j,k)*conjg(cxyz(i,j,k)))
    if (n_interp==1) then ! NGP
      C1k=1
    elseif (n_interp==2) then ! CIC
      C1k=1-(2./3.)*sin(pi*kx/nw_global)**2
    elseif (n_interp==3) then ! TSC
      C1k=1-sin(pi*kx/nw_global)**2+(2./15.)*sin(pi*kx/nw_global)**4
    endif
    Dk=(C1k(1)*C1k(2)*C1k(3))/n_particle
    xi(3,ibin)=xi(3,ibin)+amp11 ! raw power
    xi(4,ibin)=xi(4,ibin)+(amp11-Dk) ! P_r(k)
  enddo
  enddo
  enddo
  sync all
  call system_clock(t_xi2,t_xi_rate)
  ! print*,'    auto power old  elapsed time =',real(t_xi2-t_xi1)/t_xi_rate,'secs';
  if (head) then ! in head node, reduce and recover P(k)
    do i=2,nn**3
      xi=xi+xi(:,:)[i]
    enddo
    xi(2,:)=xi(2,:)/xi(1,:)
    xi(3,:)=xi(3,:)/xi(1,:) ! raw power
    xi(4,:)=xi(4,:)/xi(1,:) ! raw power - Dk
    xi(5,:)=xi(4,:)
    call pk_correction(xi,n_interp,3,alpha)
    if (.not. isnan(alpha)) call pk_correction(xi,n_interp,3,alpha)
    if (.not. isnan(alpha)) call pk_correction(xi,n_interp,3,alpha)
    if (.not. isnan(alpha)) call pk_correction(xi,n_interp,3,alpha)
    ! divide and normalize
    xi(2,:)=xi(2,:)*(2*pi)/box ! k_phys  
    xi(3,:)=xi(3,:)*(box**3) ! power_phys
    xi(4,:)=xi(4,:)*(box**3) ! power_phys
    xi(5,:)=xi(5,:)*(box**3) ! power_phys
  endif
  sync all
endsubroutine auto_power

subroutine pk_correction(xi,p,n_int,alpha)
  use omp_lib
  implicit none
  integer i,j,k,n_int,in,jn,kn,ibin,nplocal,icore,p
  real alpha,kvec(3),kmag,kmagn,kvecn(3),ks(3),Wk2Pk,Pk,cdata(0:nbin,0:ncore,3),xi(10,0:nbin)[*]

  call omp_set_num_threads(ncore)
  alpha=(log(interp1(xi(2,:),xi(5,:),real(nyquist)))-log(interp1(xi(2,:),xi(5,:),real(nyquist)/2)))/log(2.)
  print*,'  alpha =',alpha
  if (.not. isnan(alpha)) then
    print*,'pk_correction: p,n_int =',p,n_int
    print*,'  P(k_N),P(k_N/2) =',interp1(xi(2,:),xi(4,:),real(nyquist)),interp1(xi(2,:),xi(4,:),real(nyquist)/2)
    cdata=0
    call system_clock(t1,t_rate)
    !$omp paralleldo default(shared) schedule(dynamic)&
    !$omp& private(i,icore,j,k,kvec,kmag,ibin,Wk2Pk,Pk,in,jn,kn,kvecn,kmagn,ks)
    do i=1,nyquist+1
      icore=omp_get_thread_num()+1
      do j=1,i
        do k=1,j
          kvec=[i,j,k]-1.0
          kmag=norm2(kvec)
          ibin=nint(kmag)
          Wk2Pk=0
          Pk=kmag**alpha
          do in=-n_int,n_int
          do jn=-n_int,n_int
          do kn=-n_int,n_int
            kvecn=kvec+[in,jn,kn]*nw_global
            kmagn=norm2(kvecn)
            ks=pi*kvecn/nw_global
            Wk2Pk=Wk2Pk+(product(merge(1.,sin(ks)/ks,ks==0))**(2*p)) * (kmagn**alpha)
          enddo
          enddo
          enddo
          cdata(ibin,icore,:)=cdata(ibin,icore,:)+[1.,kmag,Wk2Pk/Pk]
        enddo
      enddo
    enddo
    !$omp endparalleldo
    call system_clock(t2,t_rate)
    print*, '  integration time =',real(t2-t1)/t_rate,'secs'
    cdata(1:nbin,0,:)=sum(cdata(1:nbin,1:ncore,:),dim=2)
    cdata(1:nbin,0,2)=cdata(1:nbin,0,2)/cdata(1:nbin,0,1)
    cdata(1:nbin,0,3)=cdata(1:nbin,0,3)/cdata(1:nbin,0,1)
    xi(5,1:nbin)=xi(4,1:nbin)/cdata(1:nbin,0,3)
  endif
endsubroutine
  
  real function interp1(xdata,ydata,xq)
    implicit none
    integer(4) i_mid,i1,i2
    real xdata(nbin),ydata(nbin),xq
    i1=1; i2=nbin
    do while (i2-i1>1)
      i_mid=(i1+i2)/2
      if (xq>xdata(i_mid)) then
        i1=i_mid
      else
        i2=i_mid
      endif
    enddo
    interp1=ydata(i1)+(xq-xdata(i1))/(xdata(i2)-xdata(i1))*(ydata(i2)-ydata(i1))
  endfunction

subroutine density_to_potential(cube1)
  implicit none
  integer i,j,k,ig,jg,kg
  real kr,kx(3),cube1(nw,nw,nw)
  print*,'convert density to potential'
  rho1=cube1
  call pencil_fft_forward
  cxyz=cxyz/nw_global/nw_global/nw_global
  sync all

  do k=1,npen
  do j=1,nw
  do i=1,nyquist+1
    kg=(nn*(icz-1)+icy-1)*npen+k
    jg=(icx-1)*nw+j
    ig=i
    kx=mod((/ig,jg,kg/)+nyquist-1,nw_global)-nyquist
    kx=2*sin(pi*kx/nw_global)
    kr=kx(1)**2+kx(2)**2+kx(3)**2
    kr=max(kr,1.0/nw_global**2)
    cxyz(i,j,k) = -4*pi/kr * cxyz(i,j,k)
  enddo
  enddo
  enddo
  cxyz(1,1,1)=0
  sync all

  open(11,file=output_name('phik'),status='replace',access='stream')
    write(11) cxyz
  close(11)
  sync all

  call pencil_fft_backward
  sync all

  open(11,file=output_name('phi'),status='replace',access='stream')
    write(11) rho1
  close(11)
  sync all

endsubroutine

subroutine cdm_2_nu(fn_tf)
    use omp_lib
    use parameters
    use pencil_fft
    implicit none

    character(300) :: fn_tf
    integer i,j,k,ig,jg,kg
    real kr,kx(3)
    character(20) str_z
    real kh_lin(npbin),kh_lin_log(npbin),tf_F(npbin),tf_F_log(npbin)
    real tf1(nw_global/2+1,nw,npen),k_tf1(nw_global/2+1,nw,npen)

    write(str_z,'(f8.4)') z_powerpoint(sim%cur_powerpoint)
    if (head) print*,'z=',str_z,npbin,trim(fn_tf)
    open(10,file=trim(fn_tf),status='old',access='stream'); read(10) tf_F; close(10)
    tf_F = (tf_F-(1-f_nu))/f_nu
    if (head) print*, 'tf_F',tf_F(1),tf_F(2),tf_F(npbin)
    kh_lin = -1
    open(13,file=nupath//'k_values.txt',status='old')
    read(13, '(f15.12)') kh_lin(1:npbin)
    close(13)
    kh_lin_log = log(kh_lin)
    if (head) print*,'k',kh_lin(1),kh_lin(2),kh_lin(npbin),icx,icy,icz,npen,nw,nyquist
    tf_F_log = log(tf_F)
    call pencil_fft_forward
    
    tf1 = 1
    if (Mass_nu > 0) then
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
            k_tf1(i,j,k) = kr
        enddo
        enddo
        enddo
        !$omp endparalleldo
    endif

    cxyz=cxyz*tf1
    call pencil_fft_backward

endsubroutine 

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


endmodule
