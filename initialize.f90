subroutine initialize
   use variables
   use cubefft
   use pencil_fft
   use power_nu
   implicit none

   include 'fftw3.f'

   logical,parameter :: read_Gks=.false.
   integer i,j,k,l

   real omega_m_zi,omega_l_zi


   istep=0; tictoc=0; tcat=0;
   sync all
   if (this_image()==1) then
      print*, ''
      print*, 'CUBE run on',int(nn**3,2),'images  x',int(ncore,1),'cores'
      print*, '  call geometry'
      print*,'  opath:',opath
   endif
   sync all
   call geometry
   call omp_set_num_threads(ncore)
   call omp_set_max_active_levels(2)
   ! call omp_set_nested(.true.)
   if (head) then
      print*,'  omp_get_max_threads() =',omp_get_max_threads()
      print*,'  omp_get_num_procs()   =',omp_get_num_procs()
      print*,'  omp_set_threads()     =',ncore
   endif
   sync all

   !t=0
   dt=0
   dt_old=0
   da=5


   ! sim%cur_checkpoint=4
   sim%cur_halofind=1
   z_checkpoint=-9999
   z_halofind=-9999
   checkpoint_step=.false.
   halofind_step=.false.
   final_step=.false.

   open(16,file='z_checkpoint.txt',status='old')
   do i=1,nmax_redshift-1
      read(16,end=71,fmt='(f8.4)') z_checkpoint(i)
   enddo
71 n_checkpoint=i-1
   close(16)
   if (n_checkpoint==0) stop 'z_checkpoint.txt empty'
# ifdef HALOFIND
   open(16,file='z_halofind.txt',status='old')
   do i=1,nmax_redshift-1
      read(16,end=81,fmt='(f8.4)') z_halofind(i)
   enddo
81 n_halofind=i-1
   close(16)
   if (n_halofind==0) stop 'z_halofind.txt empty'
   n_halofind=n_halofind[1]
   z_halofind(:)=z_halofind(:)[1]
# endif
   sim%tau=0!-3/sqrt(1./(1+z_checkpoint(sim%cur_checkpoint)))
   sync all
   n_checkpoint=n_checkpoint[1]
   z_checkpoint(:)=z_checkpoint(:)[1]
   !! create output directories
   if (head) then
      call system('mkdir -p '//opath//'/image'//image2str(image))
      call system('mkdir -p '//nupath//'tf')
      call system('mkdir -p '//nupath//'Pk')
      call system('mkdir -p '//nupath//'Pk_m')
      call system('mkdir -p '//nupath//'Pk_nu')
      call system('mkdir -p '//nupath//'Pk_nus')
   endif
   ! ! ! sync all


   open(10,file=nupath//'s_a_tau_H.txt',form='formatted')
   read(10,*) stime
   read(10,*) s2a
   read(10,*) s2tau
   read(10,*) s2H
   close(10)

   sim%a = 1./(1+z_checkpoint(1))

   ia = 1
   do while(s2a(ia+1)>=sim%a .and. ia < istep_max)
      ia = ia+1
   enddo
   H_i = s2H(ia) + ((s2H(ia+1)-s2H(ia))/(s2a(ia+1)-s2a(ia)))*(sim%a-s2a(ia))
   H_0 = h0*100
   s_fi = 0
   s_fi = sf_a(1./(1+z_checkpoint(1)))

   !nu
   s_f=0
   z_powerpoint=-9999
   power_step=.false.
   a_step = -1
   tau_step = -1
   sPk_step = 0
   Pk_cb_check = 0
   Pk_nu_check = 0
   if (Mass_nu > 0) then
      sim%cur_powerpoint=1
      if (head) then
         print*,'init nu_info'
         if (calculate_PK == 2) then
            print*,'  nfg nfg_global,npbin',nfg,nfg_global,npbin
         endif
         print*,'  z_nu_start f_nu',1/a_nu-1,f_nu
         print*,'  nupath:',nupath
      endif

      open(16,file=nupath//'z_powerpoint.txt',status='old')
      do i=1,nmax_redshift-1
         read(16,end=72,fmt='(f8.4)') z_powerpoint(i)
      enddo
72    n_powerpoint=i-1
      close(16)
      if (n_powerpoint==0) stop 'z_powerpoint.txt empty'
      n_powerpoint=n_powerpoint[1]
      z_powerpoint(:)=z_powerpoint(:)[1]

      do i=1,n_powerpoint ! Read the three earliest Pk_nu for linear interpolation
         write(str_z,'(f8.4)') z_powerpoint(i)
         open(10,file=nupath//'Pk_cb_'//trim(adjustl(str_z))//'.txt',form='formatted')
         read(10,*) Pk_cb_check(:,i)
         close(10)
         open(10,file=nupath//'Pk_nu_'//trim(adjustl(str_z))//'.txt',form='formatted')
         read(10,*) Pk_nu_check(:,i)
         close(10)
      enddo
      sPk_step(:,0)=Pk_cb_check(:,1)

      open(10,file=nupath//'IC/Pnu_ic.txt',form='formatted')
      read(10,*) Pk_nu_ic
      close(10)
      sq_Pk_nu_ic = sqrt(Pk_nu_ic)
      ! print*,'  Pk_nu_ic',Pk_nu_ic(1:3),Pk_nu_ic(npbin/3),Pk_nu_ic(npbin/3*2),Pk_nu_ic(npbin)

      kh_lin = -1

      omega_m_zi=(H_0/H_i)**2*omega_m*(1/(merge(real(a_nu,4),1./(1+z_checkpoint(1)),a_nu>0.0)))**3
      omega_l_zi=(H_0)**2/H_i**2*omega_l
      f_growth_nu=omega_m_zi**(4.0/7.0)+omega_l_zi/70.0*(1+omega_m_zi/2.0)


      ! stop

      kh_lin = -1
      open(13,file=nupath//'k_values.txt',status='old')
      read(13, '(f15.12)') kh_lin(1:npbin)
      close(13)
      kh_lin_log = log(kh_lin)
      if (head) print*,'k',kh_lin(1),kh_lin(2),kh_lin(npbin)

      do i=1,3
         k_mnu_2d(:,i)=spread(0, dim=1, ncopies=npbin)
         if (m_nu(i) /= 0) k_mnu_2d(:,i)=kh_lin/spread(m_nu(i), dim=1, ncopies=npbin)
      enddo
   else
      sim%cur_powerpoint = -1
   endif



! # ifdef HALOFIND
!    if (head) then
!       print*, ''
!       print*, 'runtime halofind information'
!       print*, '  ',z_checkpoint(1),'< CDM initial conditions'
!       do i=1,n_halofind
!          print*, '  ',z_halofind(i)
!       enddo
!    endif
!    sync all
! # endif

!    if (head .and. mass_nu>0) then
!       print*, ''
!       print*, 'powerpoint information'
!       print*, '  ',z_checkpoint(1),'< CDM initial conditions'
!       do i=2,n_powerpoint
!          print*, '  ',z_powerpoint(i)
!       enddo
!    endif
!    sync all

!    if (head) then
!       print*, ''
!       print*, 'checkpoint information'
!       print*, '  ',z_checkpoint(1),'< CDM initial conditions'
!       do i=2,n_checkpoint
!          print*, '  ',z_checkpoint(i)
!       enddo
!    endif
!    sync all

   if (head) then
      print*,'  initialize Green''s functions'
      print*,'    nc,ngt,nft,ng =',int(nc,2),int(ngt,2),int(nft,2),int(ng,2)
   endif
   allocate(Gk1(nc_global/2+1,nc,npen))
   allocate(Gk2(ngt/2+1,ngt,ngt))
   allocate(Gk3_2(nft(2)/2+1,nft(2),nft(2)))
   allocate(Gk3_4(nft(3)/2+1,nft(3),nft(3)))
   allocate(Gk3_6(nft(4)/2+1,nft(4),nft(4)))
   allocate(Gk3_8(nft(5)/2+1,nft(5),nft(5)))
   ! allocate(Gk3_12(nft(6)/2+1,nft(6),nft(6)))
   !allocate(Gk3_16(nft(7)/2+1,nft(7),nft(7)))

   !call Green_3D(Gk1,  nc,nc/2+1,  nc, nc, apm1c,    0., real(ratio_cs))
   if (read_Gks) then
      open(10,file=output_dir()//'Gk1'//output_suffix(),access='stream')
      read(10) Gk1
      close(10)
      open(10,file=output_dir()//'Gk2'//output_suffix(),access='stream')
      read(10) Gk2
      close(10)
      open(10,file=output_dir()//'Gk3_2'//output_suffix(),access='stream')
      read(10) Gk3_2
      close(10)
      open(10,file=output_dir()//'Gk3_4'//output_suffix(),access='stream')
      read(10) Gk3_4
      close(10)
      open(10,file=output_dir()//'Gk3_6'//output_suffix(),access='stream')
      read(10) Gk3_6
      close(10)
      open(10,file=output_dir()//'Gk3_8'//output_suffix(),access='stream')
      read(10) Gk3_8
      close(10)
      ! open(10,file=output_dir()//'Gk3_12'//output_suffix(),access='stream')
      ! read(10) Gk3_12
      ! close(10)
      !open(10,file=output_dir()//'Gk3_16'//output_suffix(),access='stream')
      !read(10) Gk3_16
      !close(10)
   else
      call tic(31)
      call Green_3D(Gk1,nc_global,nc_global/2+1, nc,npen,1, apm1c,   0., real(ratio_cs))
      call toc(31)
      open(10,file=output_dir()//'Gk1'//output_suffix(),status='replace',access='stream')
      write(10) Gk1
      close(10)
      call tic(32)
      call Green_3D(Gk2,      ngt,      ngt/2+1,ngt, ngt,0, apm2,  apm1,             1.)
      call toc(32)
      open(10,file=output_dir()//'Gk2'//output_suffix(),status='replace',access='stream')
      write(10) Gk2
      close(10)

      call tic(33)
      call Green_3D(Gk3_2, nft(2), nft(2)/2+1,nft(2), nft(2),0, apm3f,apm2f(2), 1./ratio_sf(2))
      call toc(33)
      open(10,file=output_dir()//'Gk3_2'//output_suffix(),status='replace',access='stream')
      write(10) Gk3_2
      close(10)
      call tic(34)
      call Green_3D(Gk3_4, nft(3), nft(3)/2+1,nft(3), nft(3),0, apm3f,apm2f(3), 1./ratio_sf(3))
      call toc(34)
      open(10,file=output_dir()//'Gk3_4'//output_suffix(),status='replace',access='stream')
      write(10) Gk3_4
      close(10)
      call tic(35)
      call Green_3D(Gk3_6, nft(4), nft(4)/2+1,nft(4), nft(4),0, apm3f,apm2f(4), 1./ratio_sf(4))
      call toc(35)
      open(10,file=output_dir()//'Gk3_6'//output_suffix(),status='replace',access='stream')
      write(10) Gk3_6
      close(10)

      call tic(36)
      call Green_3D(Gk3_8, nft(5), nft(5)/2+1,nft(5), nft(5),0, apm3f,apm2f(5), 1./ratio_sf(5))
      call toc(36)
      open(10,file=output_dir()//'Gk3_8'//output_suffix(),status='replace',access='stream')
      write(10) Gk3_8
      close(10)

      ! call tic(37)
      ! call Green_3D(Gk3_12, nft(6), nft(6)/2+1,nft(6), nft(6),0, apm3f,apm2f(6),1./ratio_sf(6))
      ! call toc(37)
      ! open(10,file=output_dir()//'Gk3_12'//output_suffix(),status='replace',access='stream')
      ! write(10) Gk3_12
      ! close(10)

      !call tic(38)
      !call Green_3D(Gk3_16, nft(7), nft(7)/2+1,nft(7), nft(7),0, apm3f,apm2f(7),1./ratio_sf(7))
      !call toc(38)
      !open(10,file=output_dir()//'Gk3_16'//output_suffix(),status='replace',access='stream')
      !write(10) Gk3_16
      !close(10)

   endif
   ! if (head) print*,'  Gk3_8 ',minval(Gk3_8 ),maxval(Gk3_8 )
   ! if (head) print*,'  Gk3_12',minval(Gk3_12),maxval(Gk3_12)
   !print*,'Gk3',Gk3(1:3,1,1)
   !allocate(Gk0(ng/2+1,ng,ng))
   !call Green_3D(Gk0,  ng,ng/2+1,  ng, ng, apm2,1,     0.,             1.)

   !sync all
   !if (this_image()==2) then
   !  print*,Gk1(1:2,1,1)
   !  print*,Gk2(1:2,1,1)
   !  print*,Gk3(1:2,1,1)
   !endif


   !nu
   if (Mass_nu > 0) then
      tf1 = (1-f_nu)*Gk1
      tf2 = (1-f_nu)*Gk2
      tf3_2  = (1-f_nu)*Gk3_2
      tf3_4  = (1-f_nu)*Gk3_4
      tf3_6  = (1-f_nu)*Gk3_6
      tf3_8  = (1-f_nu)*Gk3_8
      ! tf3_12 = (1-f_nu)*Gk3_12
   else
      tf1 = Gk1
      tf2 = Gk2
      tf3_2  = Gk3_2
      tf3_4  = Gk3_4
      tf3_6  = Gk3_6
      tf3_8  = Gk3_8
      ! tf3_12 = Gk3_12
      !deallocate(Gk1,,Gk2,Gk3_2,Gk3_4,Gk3_6,Gk3_8,Gk3_12)
      deallocate(Gk1,Gk2,Gk3_2,Gk3_4,Gk3_6,Gk3_8)!,Gk3_12)
   endif

   if (head) print*, '  call create_cubefft_plan'
   call tic(40)
   call create_cubefft_plan
   call toc(40)
   if (head) print*, '    elapsed time =',tcat(5,0),'secs'
   sync all

   if (head) print*, '  call create_penfft_plan nw = ',nw
   call tic(39)
   call create_penfft_plan
   call toc(39)
   if (head) print*, '    elapsed time =',tcat(4,0),'secs'
   sync all

   if (head) print*,'  initialize tile and subtile index'
   ipm2=0; ipm3=0
   do itz=1,nnt
      do ity=1,nnt
         do itx=1,nnt
            ipm2=ipm2+1
            ixyz2(:,ipm2)=[itx,ity,itz]
            do ifz=1,nns
               do ify=1,nns
                  do ifx=1,nns
                     ipm3=ipm3+1
                     ixyz3(:,ipm3)=[ifx,ify,ifz,itx,ity,itz]
                  enddo
               enddo
            enddo
         enddo
      enddo
   enddo

   if (head) print*,'  initialize PP neighbors'
   l=0
   do k=-nrange,-1
      do j=-nrange,nrange
         do i=-nrange,nrange
            l=l+1
            ijk(:,l)=[i,j,k]
         enddo
      enddo
   enddo
   k=0
   do j=-nrange,-1
      do i=-nrange,nrange
         l=l+1
         ijk(:,l)=[i,j,k]
      enddo
   enddo
   j=0
   do i=-nrange,-1
      l=l+1
      ijk(:,l)=[i,j,k]
   enddo

endsubroutine
