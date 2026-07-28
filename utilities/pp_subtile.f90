subroutine pp
  use omp_lib
  use variables
  implicit none
  save
  if (head) print*, 'PP',nrange,n_neighbor,app,appr,pp_range
  call system_clock(t1,t_rate)
  call tic(14)
  f2_max_pp=0
  do itz=1,nnt
  do ity=1,nnt
  do itx=1,nnt
    do ifz=1,nns
    do ify=1,nns
    do ifx=1,nns
     call PP_Force()
    enddo
    enddo
    enddo
  enddo
  enddo
  enddo
  sim%dt_pp = 1. / (sqrt(f2_max_pp)*a_mid*G_grid)
  call toc(14)
  call system_clock(t2,t_rate)
  if (head) print*, '  elapsed time =',real(t2-t1)/t_rate,'secs'
endsubroutine

subroutine PP_Force
  use omp_lib
  use variables
  implicit none
  save
  integer i,j,k,ilayer,np,l,idx(3),ip1,ip2,i_neighbor,nptile
  integer rcp, npgrid
  integer(8) ip,nzero,ip_offset
  real xvec(3),vreal(3),rvec(3),rmag,fpp(3)
  integer,allocatable :: ll(:),hoc(:,:,:),ip_local(:)
  real,allocatable :: xf(:,:),vf(:,:),af(:,:)

  rcp=ratio_sf ! = 4,6,8
  npgrid=ntt*rcp
  nl=([ifx,ify,ifz]-1)*ntt+1
  nh=[ifx,ify,ifz]*ntt
  nptile=sum(rhoc(nl(1)-1:nh(1)+1,nl(2)-1:nh(2)+1,nl(3)-1:nh(3)+1,itx,ity,itz))
  allocate(ll(nptile),hoc(1-rcp:npgrid+rcp,1-rcp:npgrid+rcp,1-rcp:npgrid+rcp))
  allocate(ip_local(nptile),xf(3,nptile),vf(3,nptile),af(3,nptile))

  call system_clock(tt1,t_rate)
  ip_offset=idx_b_l(0,0,itx,ity,itz)+sum(rhoc(:-1,0,0,itx,ity,itz))
  hoc=0; ll=0; ip1=0; af=0
  do k=nl(3)-1,nh(3)+1
  do j=nl(2)-1,nh(2)+1
  do i=nl(1)-1,nh(1)+1
    np=rhoc(i,j,k,itx,ity,itz)
    nzero=idx_b_r(j,k,itx,ity,itz)-sum(rhoc(i:,j,k,itx,ity,itz))
    do l=1,np
      ip=nzero+l; ip1=ip1+1; ip_local(ip1)=ip-ip_offset
      xvec=[i,j,k]-1+(int(xp(:,ip)+ishift,izipx)+rshift)*x_resolution
      xf(:,ip1)=xvec*ratio_cs
      vf(:,ip1)=tan(pi*real(vp(:,ip))/real(nvbin-1))/(sqrt(pi/2)/(sigma_vi*vrel_boost))
      idx=floor(rcp*(xvec-nl))+1
      ll(ip1)=hoc(idx(1),idx(2),idx(3))
      hoc(idx(1),idx(2),idx(3))=ip1
    enddo
  enddo
  enddo
  enddo
  call system_clock(tt2,t_rate)
  !print*,'  ll',real(tt2-tt1)/t_rate,'secs'

  call system_clock(tt1,t_rate)
  do ilayer=0,1
    !$omp paralleldo default(shared) schedule(dynamic)&
    !$omp& private(k,j,i,ip1,ip2,rvec,rmag,fpp,i_neighbor,idx)
    do k=ilayer,npgrid+1,2
    do j=0,npgrid+1
    do i=0,npgrid+1
      ip1=hoc(i,j,k)
      do while(ip1/=0) ! particle A
        ip2=ll(ip1)
        do while (ip2/=0) ! particle B in same cell
          rvec=xf(:,ip2)-xf(:,ip1); rmag=norm2(rvec)
          if (rmag<pp_range) then
            fpp=rvec/rmag*(F_ra(rmag,app)-F_ra(rmag,appr))
            af(:,ip1)=af(:,ip1)+fpp
            af(:,ip2)=af(:,ip2)-fpp
          endif
          ip2=ll(ip2)
        enddo
        do i_neighbor=1,n_neighbor ! neighbor cells
          idx=[i,j,k]+ijk(:,i_neighbor)
          if (maxval(idx)>ngp+ngb .or. minval(idx)<1-ngb) print*,[i,j,k],idx
          ip2=hoc(idx(1),idx(2),idx(3))
          do while (ip2/=0) ! particle B in neighbor cell
            rvec=xf(:,ip2)-xf(:,ip1); rmag=norm2(rvec)
            if (rmag<pp_range) then
              fpp=rvec/rmag*(F_ra(rmag,app)-F_ra(rmag,appr))
              af(:,ip1)=af(:,ip1)+fpp
              af(:,ip2)=af(:,ip2)-fpp
            endif  
            ip2=ll(ip2)
          enddo
        enddo
        ip1=ll(ip1)
      enddo
    enddo
    enddo
    enddo
    !$omp endparalleldo
  enddo
  call system_clock(tt2,t_rate)
  !print*,'  PP',real(tt2-tt1)/t_rate,'secs'

  call system_clock(tt1,t_rate)
  af=af*sim%mass_p_cdm
  do k=1,npgrid
  do j=1,npgrid
  do i=1,npgrid
    ip1=hoc(i,j,k)
    do while (ip1/=0)
      f2_max_pp=max(f2_max_pp,sum(af(:,ip1)**2))
      ip=ip_offset+ip_local(ip1)
      vreal=vf(:,ip1)+af(:,ip1)*a_mid*dt/6/pi
      vp(:,ip)=nint(real(nvbin-1)*atan(sqrt(pi/2)/(sigma_vi*vrel_boost)*vreal)/pi,kind=izipv)
      ip1=ll(ip1)
    enddo
  enddo
  enddo
  enddo
  call system_clock(tt2,t_rate)
  !print*,'  max',real(tt2-tt1)/t_rate,'secs'
  deallocate(ll,hoc,ip_local,xf,vf,af)
endsubroutine
