! #define Cpower_LE
! #define Cpower_LN
#define Cpower_LU
#define Cpower_NU
#define Cpower_EU
! #define Emode
! #define dspE
! #define dsp

program displacement
   use omp_lib
   use parameters
   use pencil_fft
   use powerspectrum
   use ieee_arithmetic
   implicit none
   save
   include 'fftw3.f'




   integer(4),allocatable :: rhoc(:,:,:,:,:,:),rhoc0(:,:,:,:,:,:)
   integer(izipi),allocatable :: pid(:),shift_count(:)
   integer(izipx),allocatable :: xp(:,:)

   integer :: ip,iq(3),i_dim,idx1(3),idx2(3)
   real :: dx1(3),dx2(3),xi(10,nbin)[*]
   integer(8) np,istat,nthreads,plan,iplan,nlast,pid8,irange(2,3)



#ifdef Emode
   complex,allocatable :: cdiv(:,:,:),cphi(:,:,:)
   real,allocatable :: cube1(:,:,:),cube0(:,:,:)
#endif

#ifdef dspE
   complex,allocatable :: kpsi(:,:,:)
   real(4),allocatable :: rho_grid(:,:,:)
#endif
   integer :: i,j,k,l,ix,iy,iz,ilayer,itx,ity,itz,kg,jg,ig,cur_checkpoint
   real :: pos0(3),pos1(3),dpos(3),kx(3),pdim(3)
   real,allocatable :: dsp(:,:,:,:),dsp_shift(:,:,:)

   real,allocatable :: rho_L(:,:,:),rho_c(:,:,:),rho_mu(:,:,:),rho_E(:,:,:)


   call geometry


   call create_penfft_plan


   print*, 'Displacement field analysis on resolution:'
   print*, 'ng=',ng
   print*, 'checkpoint at:'
   open(16,file='./z_checkpoint.txt',status='old')
   do i=1,nmax_redshift
      read(16,end=71,fmt='(f8.4)') z_checkpoint(i)
      print*, z_checkpoint(i)
   enddo
71 n_checkpoint=i-1
   close(16)
   print*,''
   nthreads=omp_get_max_threads()
   print*, '    omp_get_max_threads() =',nthreads
   call omp_set_num_threads(ncore)
   allocate(rhoc(nt,nt,nt,nnt,nnt,nnt),rhoc0(nt,nt,nt,nnt,nnt,nnt))
   irange(:,1) = [(icx-1)*ng+1,icx*ng]
   irange(:,2) = [(icy-1)*ng+1,icy*ng]
   irange(:,3) = [(icz-1)*ng+1,icz*ng]


   do cur_checkpoint= 1,n_checkpoint,-1
      print*, ''
      print*,'==========================================================='


      sim%cur_checkpoint=cur_checkpoint
      print*, 'Start analyzing redshift ',z2str(z_checkpoint(sim%cur_checkpoint))
      open(11,file=output_name('info'),access='stream'); read(11) sim; close(11)

#     ifdef dsp
      if (sim%npglobal /= ng_global**3) stop 'ng^3 != npglobal'
      print*, '    ng_global =',sim%npglobal,sim%nplocal

      allocate(xp(3,sim%nplocal),pid(sim%nplocal))
      open(11,file=output_name('xp'),access='stream'); read(11) xp; close(11)
      open(11,file=output_name('np'),access='stream'); read(11) rhoc; close(11)
      open(11,file=output_name('id'),access='stream'); read(11) pid; close(11)

      allocate(dsp(3,ng,ng,ng))!,dsp_shift(3,nn*3,sim%nplocal/100),shift_count(nn**3))
      print*, 'rhoc',minval(rhoc),maxval(rhoc),sum(rhoc)
      ! print*,rhoc
      ! stop

      print*,'PID range:',minval(pid),maxval(pid)
      dsp=0; nlast=0
      do iz=1,nnt
         do iy=1,nnt
            do ix=1,nnt
               ! aa
               !   !$omp paralleldo default(shared) schedule(dynamic)&
               !   !$omp& private(ilayer,k,nlast,j,i,l,np,ip,pid8,iq,pos0,pos1,dpos)
               do ilayer=0,ncore-1
                  do k=1+ilayer,nt,ncore
                     nlast = sum(rhoc(:,:,:,:,:,:iz-1))   &
                        + sum(rhoc(:,:,:,:,:iy-1,iz))  &
                        + sum(rhoc(:,:,:,:ix-1,iy,iz)) &
                        + sum(rhoc(:,:,:k-1,ix,iy,iz))
                     ! print*, '  working on',k,ix,iy,iz,nlast
                     do j=1,nt
                        do i=1,nt
                           ! print*, '  ',i,j,k,ix,iy,iz,rhoc(i,j,k,ix,iy,iz)
                           np=rhoc(i,j,k,ix,iy,iz)
                           ! if(rhoc(i,j,k,ix,iy,iz) /=0) print*, '  working on',i,j,k,ix,iy,iz,np,rhoc(i,j,k,ix,iy,iz)
                           ! if(np /=0) print*, '  working on',i,j,k,ix,iy,iz,np
                           do l=1,np
                              ip=nlast+l
                              if (ip > sim%nplocal .or. ip < 0) stop 'ip > nplocal'
                              pid8=pid(ip)-1
                              iq(3)=pid8/int(ng_global,4)**2
                              iq(2)=(pid8-iq(3)*int(ng_global,4)**2)/int(ng_global,4)
                              iq(1)=modulo(pid8,int(ng_global,4))
                              ! print*, '  iq,pid8',iq,pid8
                              pos0=iq+0.5
                              pos1=nt*([ix,iy,iz]-1) + [i,j,k]-1 + (int(xp(:,ip)+ishift,izipx)+rshift)*x_resolution
                              pos1=real(ng)*([icx,icy,icz]-1) + pos1*real(ng)/real(nc)
                              ! print*, '  pos0,pos1',pos0,pos1
                              ! stop
                              dpos=pos1-pos0
                              dpos=modulo(dpos+ng_global/2,real(ng_global))-ng_global/2
                              if (any(ieee_is_nan(dpos))) then
                                 print*, '  NaN detected in dpos',dpos
                                 print*, '  pos0,pos1',pos0,pos1
                                 print*, '  iq,ip,pid8',iq,ip,pid8
                              endif
                              if (dpos(1) == 0 .or. dpos(2) == 0 .or. dpos(3) == 0 ) then
                                 print*, '  0 detected in dpos',dpos
                                 print*, '  pos0,pos1',pos0,pos1
                                 print*, '  iq,ip,pid8',iq,ip,pid8
                              endif
                              ! if (iq(1)<irange(1,1) .or. iq(1)>irange(1,2) .or. &
                              !     iq(2)<irange(2,1) .or. iq(2)>irange(2,2) .or. &
                              !     iq(3)<irange(3,1) .or. iq(3)>irange(3,2)) then

                              !   !缺少多节点逻辑
                              !   cycle
                              ! endif
                              dsp(:,iq(1)+1,iq(2)+1,iq(3)+1)=dpos
                           enddo
                           nlast = nlast+np

                        enddo
                     enddo
                  enddo
               enddo
               ! asa
               ! !$omp endparalleldo
            enddo
         enddo
      enddo
      deallocate(xp,pid,rhoc)

      do i_dim=1,3
         print*, 'dsp: dimension',int(i_dim,1),'min,max values ='
         print*, minval(dsp(i_dim,:,:,:)), maxval(dsp(i_dim,:,:,:))
      enddo

      if (head) print*,'Write dsp into file',output_name('dsp')
      open(15,file=output_name('dsp'),status='replace',access='stream')
      write(15) dsp
      close(15)
      call get_mu_D('dsp')
#     endif




#     ifdef Emode
      print*,''
      print*,'Start computing delta_E'
      allocate(cdiv(ng*nn/2+1,ng,npen),cphi(ng*nn/2+1,ng,npen))
      cphi=0
      cdiv=0
      do i_dim=1,3
         print*,'  working on dim',int(i_dim,1)
         rho1=dsp(i_dim,1:ng,1:ng,1:ng)
         call pencil_fft_forward
         !$omp paralleldo default(shared) schedule(dynamic)&
         !$omp& private(k,j,i,kg,jg,ig,kx,pdim)
         do k=1,npen
            do j=1,ng
               do i=1,ng*nn/2+1
                  kg=(nn*(icz-1)+icy-1)*npen+k
                  jg=(icx-1)*ng+j
                  ig=i
                  kx=modulo([ig,jg,kg]+ng/2-1,ng)-ng/2 !k
                  pdim=sin(2*pi*kx/ng)
                  cphi(i,j,k)=cphi(i,j,k)+(0,1)*rho1k(i,j,k)*pdim(i_dim)/(-sum(pdim**2)) !phik
                  cdiv(i,j,k)=cdiv(i,j,k)+(0,1)*rho1k(i,j,k)*pdim(i_dim) !c means complex
               enddo
            enddo
         enddo
         !$omp endparalleldo
      enddo ! i_dim

      if (head) then
         cphi(1,1,1)=0
         cdiv(1,1,1)=0
      endif
      sync all

      cxyz=cdiv
      print*,'  btran'
      call pencil_fft_backward
      allocate(rho_E(ng,ng,ng))
      rho_E=-rho1
      print*, minval(rho_E(1:nw,1:nw,1:nw)),maxval(rho_E(1:nw,1:nw,1:nw)),sum(rho_E(1:nw,1:nw,1:nw)*1d0)/real(nw)**3
      print*,'  write delta_E into file',output_name('delta_E')
      open(15,file=output_name('delta_E'),status='replace',access='stream')
      write(15) rho_E
      close(15)
      sync all
      deallocate(rho_E,cdiv,cphi)
#     endif

#     ifdef dspE
      print*,'dspe'
      allocate(kpsi(ng*nn/2+1,ng,npen))
      kpsi=0
      do i_dim=1,3
         rho1=dsp(i_dim,1:ng,1:ng,1:ng)
         call pencil_fft_forward
         !$omp paralleldo default(shared) schedule(dynamic)&
         !$omp& private(k,j,i,kg,jg,ig,kx)
         do k=1,npen
            do j=1,ng
               do i=1,ng*nn/2+1
                  kg=(nn*(icz-1)+icy-1)*npen+k
                  jg=(icx-1)*ng+j
                  ig=i
                  kx=mod([ig,jg,kg]+ng/2-1,ng)-ng/2
                  kpsi(i,j,k)=kpsi(i,j,k)+kx(i_dim)*cxyz(i,j,k)
               enddo
            enddo
         enddo
         !$omp endparalleldo
      enddo
      print*, 'sum of kpsi',sum(kpsi)

      open(16,file=output_name('dsp_E'),status='replace',access='stream')
      print*,'compute and write into',output_name('dsp_E')
      do i_dim=1,3
         !$omp paralleldo default(shared) schedule(dynamic)&
         !$omp& private(k,j,i,kg,jg,ig,kx)
         do k=1,npen
            do j=1,ng
               do i=1,ng*nn/2+1
                  kg=(nn*(icz-1)+icy-1)*npen+k
                  jg=(icx-1)*ng+j
                  ig=i
                  kx=mod([ig,jg,kg]+ng/2-1,ng)-ng/2
                  cxyz(i,j,k)=kpsi(i,j,k)*kx(i_dim)/sum(kx**2)
               enddo
            enddo
         enddo
         !$omp endparalleldo
         if (head) then
            cxyz(1,1,1)=0
         endif
         ! print*, 'sum of cxyz',minval(abs(cxyz)),maxval(abs(cxyz))
         call pencil_fft_backward
         print*, 'rho1',minval(rho1),maxval(rho1)
         write(16) rho1
         dsp(i_dim,:,:,:) = rho1
      enddo
      close(16)
      deallocate(kpsi)

      ! CIC to get the Eulerian density field
      if (head) print*, 'CIC interpolation by dsp_E'
      allocate(rho_grid(0:ng+1,0:ng+1,0:ng+1))
      rho_grid=0
      do k=1,ng
         do j=1,ng
            do i=1,ng

               pos1=[i,j,k]-0.5+dsp(:,i,j,k)
               pos1=modulo(pos1-1,real(ng))+1
               pos1=pos1-0.5

               !pos1=250;


               idx1=floor(pos1)+1
               idx2=idx1+1
               dx1=idx1-pos1
               dx2=1-dx1

               !print*, idx1;stop

               rho_grid(idx1(1),idx1(2),idx1(3))=rho_grid(idx1(1),idx1(2),idx1(3))+dx1(1)*dx1(2)*dx1(3)
               rho_grid(idx2(1),idx1(2),idx1(3))=rho_grid(idx2(1),idx1(2),idx1(3))+dx2(1)*dx1(2)*dx1(3)
               rho_grid(idx1(1),idx2(2),idx1(3))=rho_grid(idx1(1),idx2(2),idx1(3))+dx1(1)*dx2(2)*dx1(3)
               rho_grid(idx1(1),idx1(2),idx2(3))=rho_grid(idx1(1),idx1(2),idx2(3))+dx1(1)*dx1(2)*dx2(3)
               rho_grid(idx1(1),idx2(2),idx2(3))=rho_grid(idx1(1),idx2(2),idx2(3))+dx1(1)*dx2(2)*dx2(3)
               rho_grid(idx2(1),idx1(2),idx2(3))=rho_grid(idx2(1),idx1(2),idx2(3))+dx2(1)*dx1(2)*dx2(3)
               rho_grid(idx2(1),idx2(2),idx1(3))=rho_grid(idx2(1),idx2(2),idx1(3))+dx2(1)*dx2(2)*dx1(3)
               rho_grid(idx2(1),idx2(2),idx2(3))=rho_grid(idx2(1),idx2(2),idx2(3))+dx2(1)*dx2(2)*dx2(3)

            enddo
         enddo
      enddo
      rho_grid(1,:,:)=rho_grid(1,:,:)+rho_grid(ng+1,:,:)
      rho_grid(ng,:,:)=rho_grid(ng,:,:)+rho_grid(0,:,:)
      rho_grid(:,1,:)=rho_grid(:,1,:)+rho_grid(:,ng+1,:)
      rho_grid(:,ng,:)=rho_grid(:,ng,:)+rho_grid(:,0,:)
      rho_grid(:,:,1)=rho_grid(:,:,1)+rho_grid(:,:,ng+1)
      rho_grid(:,:,ng)=rho_grid(:,:,ng)+rho_grid(:,:,0)
      rho_grid=rho_grid-1

      print*, minval(rho_grid(1:ng,1:ng,1:ng)),maxval(rho_grid(1:ng,1:ng,1:ng)),sum(rho_grid(1:ng,1:ng,1:ng)*1d0)
      open(16,file=output_name('delta_cE'),status='replace',access='stream')
      print*,'compute and write into',output_name('delta_cE')
      write(16) rho_grid(1:ng,1:ng,1:ng)
      close(16)
      deallocate(rho_grid)

      call get_mu_D  ('dsp_E')
#     endif

#ifdef Cpower_LN
      print*,'Cpower_LN'
      allocate(rho_L(nw,nw,nw))
      open(11,file=output_dir()//'delta_L'//output_suffix(),status='old',access='stream')
      read(11) rho_L

      allocate(rho_c(nw,nw,nw))
      open(11,file=output_name('delta_c'),status='old',access='stream')
      read(11) rho_c
      call cross_power(xi,rho_L,rho_c)
      print*,'   save: ',output_name('Cpower_LN')
      open(11,file=output_name('Cpower_LN'),status='replace',access='stream')
      write(11) xi
      close(11)
      deallocate(rho_c,rho_L)
#endif


#ifdef Cpower_LE
      print*,'Cpower_LE'
      allocate(rho_L(nw,nw,nw))
      open(11,file=output_dir()//'delta_L'//output_suffix(),status='old',access='stream')
      read(11) rho_L

      allocate(rho_E(nw,nw,nw))
      open(11,file=output_name('delta_cE'),status='old',access='stream')
      read(11) rho_E
      call cross_power(xi,rho_L,rho_E)
      print*,'   save: ',output_name('Cpower_LE')
      open(11,file=output_name('Cpower_LE'),status='replace',access='stream')
      write(11) xi
      close(11)
      deallocate(rho_E,rho_L)
#endif

#ifdef Cpower_LU
      print*,'Cpower_LU'   
      allocate(rho_L(nw,nw,nw))
      open(11,file=output_dir()//'delta_L'//output_suffix(),status='old',access='stream')
      read(11) rho_L

      allocate(rho_mu(nw,nw,nw))
      open(11,file=output_name('dspu_x'),status='old',access='stream')
      read(11) rho_mu
      call cross_power(xi,rho_L,rho_mu)
      print*,'   save: ',output_name('Cpower_LU')
      open(11,file=output_name('Cpower_LU'),status='replace',access='stream')
      write(11) xi
      close(11)
      WHERE (rho_L > 1)
         rho_L = 1
      END WHERE
      WHERE (rho_mu > 1)
         rho_mu = 1
      END WHERE

      call cross_power(xi,rho_L,rho_mu)
      print*,'   save: ',output_name('Cpower_LU_low')
      open(11,file=output_name('Cpower_LU_low'),status='replace',access='stream')
      write(11) xi
      close(11)
      deallocate(rho_mu,rho_L)
#endif

#ifdef Cpower_NU
      print*,'Cpower_NU'
      allocate(rho_c(nw,nw,nw))
      open(11,file=output_name('delta_c'),status='old',access='stream')
      read(11) rho_c

      allocate(rho_mu(nw,nw,nw))
      open(11,file=output_name('dspu_x'),status='old',access='stream')
      read(11) rho_mu
      call cross_power(xi,rho_c,rho_mu)
      print*,'   save: ',output_name('Cpower_NU')
      open(11,file=output_name('Cpower_NU'),status='replace',access='stream')
      write(11) xi
      close(11)
      WHERE (rho_c > 1)
         rho_c = 1
      END WHERE
      WHERE (rho_mu > 1)
         rho_mu = 1
      END WHERE

      call cross_power(xi,rho_c,rho_mu)
      print*,'   save: ',output_name('Cpower_NU_low')
      open(11,file=output_name('Cpower_NU_low'),status='replace',access='stream')
      write(11) xi
      close(11)
      deallocate(rho_mu,rho_c)
#endif


#ifdef Cpower_EU
      print*,'Cpower_EU'
      allocate(rho_E(nw,nw,nw))
      open(11,file=output_name('delta_cE'),status='old',access='stream')
      read(11) rho_E

      allocate(rho_mu(nw,nw,nw))
      open(11,file=output_name('dspu_x'),status='old',access='stream')
      read(11) rho_mu
      call cross_power(xi,rho_E,rho_mu)
      print*,'   save: ',output_name('Cpower_EU')
      open(11,file=output_name('Cpower_EU'),status='replace',access='stream')
      write(11) xi
      close(11)
      WHERE (rho_E > 1)
         rho_E = 1
      END WHERE
      WHERE (rho_mu > 1)
         rho_mu = 1
      END WHERE

      call cross_power(xi,rho_E,rho_mu)
      print*,'   save: ',output_name('Cpower_EU_low')
      open(11,file=output_name('Cpower_EU_low'),status='replace',access='stream')
      write(11) xi
      close(11)
      deallocate(rho_mu,rho_E)
#endif

   enddo


contains


   subroutine get_mu_D(namespace)
      implicit none
      character(len=*), intent(in) :: namespace
      integer :: i_0,j_0,k_0,i_n(4),j_n(4),k_n(4),dim,l1,l2,lr
      real,allocatable :: dsp(:,:,:,:)

      integer :: idx1(3),idx2(3),nwb
      real(8)  A_mesh(3,3),mu_i,det_A,pos1(3),min_r,now_r,xp_neigh(3),mu_neigh
      real(4),allocatable :: mu(:,:,:),rho_grid(:,:,:)
      integer(4),allocatable :: hoc(:,:,:),ll(:)
      integer(8) :: ip,ip_neigh
      allocate(dsp(3,ng,ng,ng))

      dsp = 0
      print*,''
      print*,'  read:'
      print*,'    ',output_name(namespace)
      open(15,file=output_name(namespace),status='old',access='stream')
      read(15) dsp
      close(15)

      if (any(ieee_is_nan(dsp))) stop 'dsp is nan'
      nwb = maxval(abs(dsp))+1

      do i_dim=1,3
         print*, 'dsp: dimension',int(i_dim,1),'min,max values ='
         print*, minval(dsp(i_dim,:,:,:)), maxval(dsp(i_dim,:,:,:))
      enddo
      ! stop

      allocate(mu(ng,ng,ng))
      mu = 0
      ! print*,'    u     : ',minval(mu),maxval(mu),sum(mu*1d0)/real(ng)**3
      !$omp paralleldo default(shared) schedule(dynamic) private(k_0,j_0,i_0,i_n,j_n,k_n,A_mesh,dim,det_A)
      do k_0=1,ng
         do j_0=1,ng
            do i_0=1,ng
               i_n=modulo(i_0+[-2,-1,1,2]+ng-1,ng)+1
               j_n=modulo(j_0+[-2,-1,1,2]+ng-1,ng)+1
               k_n=modulo(k_0+[-2,-1,1,2]+ng-1,ng)+1
               A_mesh = 0
               do dim = 1,3
                  A_mesh(1,dim)=sum(dsp(dim,i_n,j_0,k_0)*weight)
                  A_mesh(2,dim)=sum(dsp(dim,i_0,j_n,k_0)*weight)
                  A_mesh(3,dim)=sum(dsp(dim,i_0,j_0,k_n)*weight)
               enddo
               A_mesh(1,1) = A_mesh(1,1)+1
               A_mesh(2,2) = A_mesh(2,2)+1
               A_mesh(3,3) = A_mesh(3,3)+1

               if (any(ieee_is_nan(A_mesh))) then
                  print*, 'i_0,j_0,k_0',i_0,j_0,k_0
                  print*, 'A_mesh',A_mesh
                  print*, ''

                  i_n=i_0+[-2,-1,1,2]
                  j_n=j_0+[-2,-1,1,2]
                  k_n=k_0+[-2,-1,1,2]
                  do dim = 1,3
                     print*,dim , 'i',dsp(dim,i_n,j_0,k_0),i_n
                     print*,dim , 'j',dsp(dim,i_0,j_n,k_0),j_n
                     print*,dim , 'k',dsp(dim,i_0,j_0,k_n),k_n
                  enddo
                  stop 'A_mesh is nan'
               endif

               det_A = A_mesh(1,1) * (A_mesh(2,2)*A_mesh(3,3) - A_mesh(2,3)*A_mesh(3,2))  &
                  + A_mesh(1,2) * (A_mesh(2,3)*A_mesh(3,1) - A_mesh(2,1)*A_mesh(3,3))  &
                  + A_mesh(1,3) * (A_mesh(2,1)*A_mesh(3,2) - A_mesh(2,2)*A_mesh(3,1))
               if (ieee_is_nan(det_A)) stop 'det_A is nan'
               if (det_A <= 0) then
                  det_A = 1e-10
               endif
               ! if (abs(det_A) < 1)  det_A = 0
               if ( 1/det_A < 0) then
                  print*, 'i_0,j_0,k_0',i_0,j_0,k_0
                  print*, 'A_mesh',A_mesh
                  print*, 'det_A',det_A
                  stop 'det_A is negative'
               endif
               mu(i_0,j_0,k_0) = 1/det_A
               if ( mu(i_0,j_0,k_0) < 0) then
                  print*, 'i_0,j_0,k_0',i_0,j_0,k_0
                  print*, 'A_mesh',A_mesh
                  print*, 'det_A',det_A
                  stop 'det_A is negative'
               endif
            enddo
         enddo
      enddo
      !$omp endparalleldo

      print*,'    u     : ',minval(mu),maxval(mu),sum(mu*1d0)/real(ng)**3
      open(21,file=output_name(namespace//'u_q'),status='replace',access='stream')
      write(21) mu
      close(21)

      if (head) print*, trim('CIC interpolation mu by ' // namespace),ng,nw
      allocate(rho_grid(0:nw+1,0:nw+1,0:nw+1),hoc(1-nwb:ng+nwb,1-nwb:ng+nwb,1-nwb:ng+nwb),ll(ng**3))
      rho_grid=0; hoc=0; ll=0; ip = 0

      do k=1,ng
         ! print*,k
         do j=1,ng
            do i=1,ng

               pos1=[i,j,k]-0.5+dsp(:,i,j,k)
               pos1=modulo(pos1-1,real(ng))+1
               pos1=pos1*real(nw)/real(ng)-0.5

               idx1=floor(pos1)+1
               idx2=idx1+1
               dx1=idx1-pos1
               dx2=1-dx1

               mu_i = mu(i,j,k)

               rho_grid(idx1(1),idx1(2),idx1(3))=rho_grid(idx1(1),idx1(2),idx1(3))+dx1(1)*dx1(2)*dx1(3)*mu_i
               rho_grid(idx2(1),idx1(2),idx1(3))=rho_grid(idx2(1),idx1(2),idx1(3))+dx2(1)*dx1(2)*dx1(3)*mu_i
               rho_grid(idx1(1),idx2(2),idx1(3))=rho_grid(idx1(1),idx2(2),idx1(3))+dx1(1)*dx2(2)*dx1(3)*mu_i
               rho_grid(idx1(1),idx1(2),idx2(3))=rho_grid(idx1(1),idx1(2),idx2(3))+dx1(1)*dx1(2)*dx2(3)*mu_i
               rho_grid(idx1(1),idx2(2),idx2(3))=rho_grid(idx1(1),idx2(2),idx2(3))+dx1(1)*dx2(2)*dx2(3)*mu_i
               rho_grid(idx2(1),idx1(2),idx2(3))=rho_grid(idx2(1),idx1(2),idx2(3))+dx2(1)*dx1(2)*dx2(3)*mu_i
               rho_grid(idx2(1),idx2(2),idx1(3))=rho_grid(idx2(1),idx2(2),idx1(3))+dx2(1)*dx2(2)*dx1(3)*mu_i
               rho_grid(idx2(1),idx2(2),idx2(3))=rho_grid(idx2(1),idx2(2),idx2(3))+dx2(1)*dx2(2)*dx2(3)*mu_i

               ip = ip + 1
               ll(ip) = hoc(idx2(1),idx2(2),idx2(3))
               hoc(idx2(1),idx2(2),idx2(3))=ip
            enddo
         enddo
      enddo
      ! deallocate(dsp)
      rho_grid(1,:,:)=rho_grid(1,:,:)+rho_grid(nw+1,:,:)
      rho_grid(nw,:,:)=rho_grid(nw,:,:)+rho_grid(0,:,:)
      rho_grid(:,1,:)=rho_grid(:,1,:)+rho_grid(:,nw+1,:)
      rho_grid(:,nw,:)=rho_grid(:,nw,:)+rho_grid(:,0,:)
      rho_grid(:,:,1)=rho_grid(:,:,1)+rho_grid(:,:,nw+1)
      rho_grid(:,:,nw)=rho_grid(:,:,nw)+rho_grid(:,:,0)

      hoc(:1,:,:) = hoc(ng-nwb+1:ng,:,:)
      hoc(:,:1,:) = hoc(:,ng-nwb+1:ng,:)
      hoc(:,:,:1) = hoc(:,:,ng-nwb+1:ng)
      hoc(ng+1:ng+nwb,:,:) = hoc(1:nwb,:,:)
      hoc(:,ng+1:ng+nwb,:) = hoc(:,1:nwb,:)
      hoc(:,:,ng+1:ng+nwb) = hoc(:,:,1:nwb)




      print*, minval(abs(rho_grid(1:nw,1:nw,1:nw))),maxval(abs(rho_grid(1:nw,1:nw,1:nw))),sum(rho_grid(1:nw,1:nw,1:nw)*1d0)/real(nw)**3
      open(16,file=output_name(namespace//'u_x'),status='replace',access='stream')
      print*,'compute and write mu into',output_name(namespace//'u_x')
      write(16) rho_grid(1:nw,1:nw,1:nw)
      close(16)

      print*,'nwb',nwb
      ! a
      ! !$omp paralleldo num_threads(nnest) default(shared) schedule(dynamic)&
      ! !$omp& private(k,j,i,min_r,l1,k_0,j_0,i_0,ip,idx1,pos1,now_r,xp_neigh,mu_neigh,lr,l2)
      do k=1,ng
         ! print*,k
         do j=1,ng
            do i=1,ng
               if ( rho_grid(i,j,k) == 0 ) then
                  min_r = ng*2
                  do l1 = 1,nwb !l1 means layer
                     do k_0=-l1,l1
                        do j_0=-l1,l1
                           do i_0=-l1,l1
                              if (hoc(i_0+i,j_0+j,k_0+k) == 0) cycle
                              if (max(abs(k_0), abs(j_0), abs(i_0)) /= l1) cycle

                              ip = hoc(i_0+i,j_0+j,k_0+k)
                              do while (ip /= 0)
                                 idx1(3) = ((ip-1) / ng**2)+1
                                 idx1(2) = ((ip-1) - (idx1(3)-1)*ng**2)/ng+1
                                 idx1(1) = mod((ip-1),ng)+1
                                 pos1 = (idx1 + dsp(:,idx1(1),idx1(2),idx1(3))) * real(nw)/real(ng)
                                 pos1 = modulo(pos1 - [i,j,k]+ng_global/2,real(ng_global))-ng_global/2
                                 now_r = sqrt(sum((pos1)**2))
                                 if (now_r < min_r .and. mu(idx1(1),idx1(2),idx1(3)) > 0) then
                                    min_r = now_r
                                    xp_neigh = (idx1 + dsp(:,idx1(1),idx1(2),idx1(3))) * real(nw)/real(ng)
                                    mu_neigh = mu(idx1(1),idx1(2),idx1(3))
                                    ! if (mu_neigh <= 0) then
                                    !   print*, 'mu_neigh <= 0 in l1',i,j,k,l1
                                    !   print*,ip,idx1
                                    !   print*,dsp(:,idx1(1),idx1(2),idx1(3))
                                    !   print*,(idx1 + dsp(:,idx1(1),idx1(2),idx1(3))) * real(nw)/real(ng)
                                    !   print*,pos1
                                    !   print*,mu(idx1(1),idx1(2),idx1(3))
                                    !   stop 'mu_neigh <= 0'
                                    ! endif
                                 endif
                                 ip = ll(ip)
                              enddo
                           enddo
                        enddo
                     enddo
                     if (min_r < ng*2 ) then
                        lr = floor(min_r)+1

                        do l2 = l1+1,lr !l2 means layer
                           do k_0=-l2,l2
                              do j_0=-l2,l2
                                 do i_0=-l2,l2
                                    if (hoc(i_0+i,j_0+j,k_0+k) == 0) cycle
                                    if (max(abs(k_0), abs(j_0), abs(i_0)) /= l2) cycle
                                    if (sqrt(real(i_0**2+j_0**2+k_0**2)) > min_r) cycle

                                    ip = hoc(i_0+i,j_0+j,k_0+k)
                                    do while (ip /= 0)
                                       idx1(3) = ((ip-1) / ng**2)+1
                                       idx1(2) = ((ip-1) - (idx1(3)-1)*ng**2)/ng+1
                                       idx1(1) = mod((ip-1),ng)+1
                                       pos1 = (idx1 + dsp(:,idx1(1),idx1(2),idx1(3))) * real(nw)/real(ng)
                                       pos1 = modulo(pos1 - [i,j,k]+ng_global/2,real(ng_global))-ng_global/2
                                       now_r = sqrt(sum((pos1)**2))
                                       if (now_r < min_r .and. mu(idx1(1),idx1(2),idx1(3)) > 0 ) then
                                          min_r = now_r
                                          xp_neigh = (idx1 + dsp(:,idx1(1),idx1(2),idx1(3))) * real(nw)/real(ng)
                                          mu_neigh = mu(idx1(1),idx1(2),idx1(3))
                                          ! if (mu_neigh <= 0) then
                                          !   print*, 'mu_neigh <= 0 in l2',i,j,k,l2
                                          !   print*,ip,idx1
                                          !   print*,dsp(:,idx1(1),idx1(2),idx1(3))
                                          !   print*,(idx1 + dsp(:,idx1(1),idx1(2),idx1(3))) * real(nw)/real(ng)
                                          !   print*,mu(idx1(1),idx1(2),idx1(3))
                                          !   stop 'mu_neigh <= 0'
                                          ! endif
                                       endif
                                       ip = ll(ip)
                                    enddo
                                 enddo
                              enddo
                           enddo
                        enddo
                        exit

                     endif
                  enddo
                  if (min_r == ng*2 ) then
                     print*, 'no particle found in the neighbor',i,j,k


                     min_r = ng*2
                     do l1 = 1,nwb !l1 means layer
                        print*, 'layer 1',l1
                        do k_0=-l1,l1
                           do j_0=-l1,l1
                              do i_0=-l1,l1
                                 if (hoc(i_0+i,j_0+j,k_0+k) == 0) cycle
                                 if (max(abs(k_0), abs(j_0), abs(i_0)) /= l1) cycle
                                 print*, ' neigh',k_0,j_0,i_0
                                 print*, '   hoc',hoc(i_0+i,j_0+j,k_0+k)

                                 ip = hoc(i_0+i,j_0+j,k_0+k)
                                 do while (ip /= 0)
                                    idx1(3) = ((ip-1) / ng**2)+1
                                    idx1(2) = ((ip-1) - (idx1(3)-1)*ng**2)/ng+1
                                    idx1(1) = mod((ip-1),ng)+1
                                    pos1 = (idx1 + dsp(:,idx1(1),idx1(2),idx1(3))) * real(nw)/real(ng)
                                    pos1 = modulo(pos1 - [i,j,k]+ng_global/2,real(ng_global))-ng_global/2
                                    now_r = sqrt(sum((pos1)**2))
                                    print*, '   pos1',pos1
                                    print*, '   now_r',now_r,min_r

                                    if (now_r < min_r) then
                                       min_r = now_r
                                       xp_neigh = (idx1 + dsp(:,idx1(1),idx1(2),idx1(3))) * real(nw)/real(ng)
                                       mu_neigh = mu(idx1(1),idx1(2),idx1(3))
                                       if (mu_neigh <= 0) then
                                          print*, 'mu_neigh <= 0 in l2',i,j,k,l2
                                          print*,ip,idx1
                                          print*,dsp(:,idx1(1),idx1(2),idx1(3))
                                          print*,(idx1 + dsp(:,idx1(1),idx1(2),idx1(3))) * real(nw)/real(ng)
                                          print*,mu(idx1(1),idx1(2),idx1(3))
                                          stop 'mu_neigh <= 0'
                                       endif
                                    endif
                                    ip = ll(ip)
                                 enddo
                              enddo
                           enddo
                        enddo
                        if (min_r < ng*2 ) then
                           lr = floor(min_r)+1

                           do l2 = l1+1,lr !l2 means layer
                              print*, ' layer 2',l2
                              do k_0=-l2,l2
                                 do j_0=-l2,l2
                                    do i_0=-l2,l2
                                       if (hoc(i_0+i,j_0+j,k_0+k) == 0) cycle
                                       if (max(abs(k_0), abs(j_0), abs(i_0)) /= l2) cycle
                                       if (sqrt(real(i_0**2+j_0**2+k_0**2)) > min_r) cycle

                                       ip = hoc(i_0+i,j_0+j,k_0+k)
                                       do while (ip /= 0)
                                          idx1(3) = ((ip-1) / ng**2)+1
                                          idx1(2) = ((ip-1) - (idx1(3)-1)*ng**2)/ng+1
                                          idx1(1) = mod((ip-1),ng)+1
                                          pos1 = (idx1 + dsp(:,idx1(1),idx1(2),idx1(3))) * real(nw)/real(ng)
                                          pos1 = modulo(pos1 - [i,j,k]+ng_global/2,real(ng_global))-ng_global/2
                                          now_r = sqrt(sum((pos1)**2))
                                          if (now_r < min_r) then
                                             min_r = now_r
                                             xp_neigh = (idx1 + dsp(:,idx1(1),idx1(2),idx1(3))) * real(nw)/real(ng)
                                             mu_neigh = mu(idx1(1),idx1(2),idx1(3))
                                             if (mu_neigh <= 0) then
                                                print*, 'mu_neigh <= 0 in l2',i,j,k,l2
                                                print*,ip,idx1
                                                print*,dsp(:,idx1(1),idx1(2),idx1(3))
                                                print*,(idx1 + dsp(:,idx1(1),idx1(2),idx1(3))) * real(nw)/real(ng)
                                                print*,mu(idx1(1),idx1(2),idx1(3))
                                                stop 'mu_neigh <= 0'
                                             endif
                                          endif
                                          ip = ll(ip)
                                       enddo
                                    enddo
                                 enddo
                              enddo
                           enddo
                           exit

                        endif
                     enddo
                     stop 'min_r == ng*2'
                  endif
                  rho_grid(i,j,k) = mu_neigh/min_r/min_r
               else if ( rho_grid(i,j,k) < 0 )  then
                  rho_grid(i,j,k) = 1e10 !inf
               endif
            enddo
         enddo
      enddo
      ! a
      ! !$omp endparalleldo

      print*, minval(rho_grid(1:nw,1:nw,1:nw)),maxval(rho_grid(1:nw,1:nw,1:nw)),sum(rho_grid(1:nw,1:nw,1:nw)*1d0)/real(nw)**3
      open(16,file=output_name(namespace//'uf_x'),status='replace',access='stream')
      print*,'compute and write mu into',output_name(namespace//'uf_x')
      write(16) rho_grid(1:nw,1:nw,1:nw)
      close(16)
      deallocate(dsp)
      deallocate(mu,rho_grid)
   endsubroutine get_mu_D

endprogram displacement
