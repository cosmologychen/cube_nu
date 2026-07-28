subroutine particle_initialization
   use parameters
#ifdef PID
   use variables, only: xp,xp_new,vp,vp_new,pid,pid_new
#else
   use variables, only: xp,xp_new,vp,vp_new
#endif
   use variables, only: sigma_vi,np_image_max,np_tile_max,sPk_step,Pk_cb_check,Pk_nu_check,sq_Pk_nu_ic,a_step,tau_step,sf_step,kh_lin,s_f
! #ifdef ZIPX
   use variables, only: rhoc
! #endif
   implicit none

   character(100) fn10,fn11,fn12,fn13,fn14,fn15,fn16,fn17,fn18,fn19,str_z
   integer i
   call tic(41)
   if (head) then
      print*, '',sim%cur_checkpoint
      print*, 'particle_initialization'
      print*, '  at redshift',z_checkpoint(sim%cur_checkpoint)
   endif
   fn10=output_name('info')
   fn11=output_name('xp')
   fn12=output_name('vp')
! #ifdef ZIPX
   fn13=output_name('np')
! #endif
   fn14=output_name('vc')
   fn15=output_name('id')

   open(10,file=fn10,status='old',access='stream')
   read(10) sim
   close(10)
   sim%a=1./(1+z_checkpoint(sim%cur_checkpoint))
   sigma_vi=sim%sigma_vi


   if (sim%izipx/=izipx .or. sim%izipv/=izipv) then
      print*, '  zip format incompatable'
      close(12)
      error stop
   endif




   allocate(xp(3,np_image_max)[nn,nn,*],xp_new(3,np_tile_max))
   allocate(vp(3,np_image_max)[nn,nn,*],vp_new(3,np_tile_max))
! #ifdef ZIPX
   allocate(rhoc(1-ncb:nt+ncb,1-ncb:nt+ncb,1-ncb:nt+ncb,nnt,nnt,nnt)[nn,nn,*])
! #endif
   !allocate(vfield(3,1-ncb:nt+ncb,1-ncb:nt+ncb,1-ncb:nt+ncb,nnt,nnt,nnt)[nn,nn,*])
#ifdef PID
   allocate(pid(np_image_max)[nn,nn,*],pid_new(np_tile_max))
#endif
! open(11,file=fn11,status='old',access='stream'); read(11) xp(:,:sim%nplocal); close(11)
! open(12,file=fn12,status='old',access='stream'); read(12) vp(:,:sim%nplocal); close(12)
! open(13,file=fn13,status='old',access='stream'); read(13) rhoc(1:nt,1:nt,1:nt,:,:,:); close(13)
! # ifdef PID
! open(15,file=fn15,status='old',access='stream'); read(15) pid(:sim%nplocal); close(15)
! if (minval(pid(:sim%nplocal))<1) then
!    print*, 'warning: pid are not all positive'
!    !stop
! endif
! # endif

   !$omp parallelsections default(shared) num_threads(4)
   !$omp section
   !omp workshare
   open(11,file=fn11,status='old',access='stream'); read(11) xp(:,:sim%nplocal); close(11)
   !omp endworkshare
   !$omp section
   open(12,file=fn12,status='old',access='stream'); read(12) vp(:,:sim%nplocal); close(12)
! #ifdef ZIPX
   !$omp section
   open(13,file=fn13,status='old',access='stream'); read(13) rhoc(1:nt,1:nt,1:nt,:,:,:); close(13)
! #endif
!!  !$omp section
!!    open(14,file=fn14,status='old',access='stream'); read(14) vfield(:,1:nt,1:nt,1:nt,:,:,:); close(14)
# ifdef PID
   !$omp section
   open(15,file=fn15,status='old',access='stream'); read(15) pid(:sim%nplocal); close(15)
   !print*, '  image',this_image(),' PID range: ',minval(pid(:sim%nplocal)),maxval(pid(:sim%nplocal))
   if (minval(pid(:sim%nplocal))<1) then
      print*, 'warning: pid are not all positive'
      !stop
   endif
# endif
   if (this_image()==1 .or. this_image()==num_images()) print*,'  from image',this_image(),' read',sim%nplocal,' particles'
   !$omp endparallelsections


   !nu
   if (sim%cur_powerpoint > 1 .and. Mass_nu > 0 ) then
      if (head) print*,Mass_nu
      ! sim%cur_powerpoint = sim%cur_powerpoint-1
      open(16,file=output_name('steps'),status='old',access='stream')
      read(16) sPk_step
      read(16) Pk_cb_check
      read(16) a_step
      read(16) tau_step
      read(16) sf_step
      read(16) kh_lin
      read(16) sq_Pk_nu_ic
      close(16)

      if (calculate_PK == -1) then
         do i=4,sim%cur_powerpoint-1 ! Read the three earliest Pk_nu for linear interpolation
            write(str_z,'(f8.4)') z_powerpoint(i)
            open(10,file=nupath//'Pk_nu_'//trim(adjustl(str_z))//'.txt',form='formatted')
            read(10,*) Pk_nu_check(:,i)
            close(10)
         enddo
      endif

      if (head) then
         print*,'init '
         print*,' s_f     :',sf_step(sim%timestep-1),sim%timestep
         print*,' Pk_0    :',sPk_step(1:4,sim%timestep)
         print*,' Pk_-1   :',sPk_step(1:4,sim%timestep-1)
         print*,' Pk_check:',Pk_cb_check(1:4,sim%cur_powerpoint)
         print*,' a       :',a_step(1),a_step(sim%timestep-2:sim%timestep)
         print*,' tau     :',tau_step(1),tau_step(sim%timestep-2:sim%timestep)
      endif
   endif
   a_step(sim%timestep) = sim%a
   tau_step(sim%timestep) = sim%tau
   if(head) print*,'  a,tau init',sim%timestep,sim%a,1/sim%a-1,sim%tau

   call toc(41)
   if (head) then
      print*,'  npglobal    =', sim%npglobal
      print*,'  mass_p_cdm =', sim%mass_p_cdm
      print*,'  vsim2phys =',sim%vsim2phys, ' (km/s)/(1.0)'
      print*,'  sigma_vi    =',sigma_vi,'(simulation unit)'
      print*,'  elapsed time =',tcat(6,0),'secs'
      print*,''
   endif
   call print_header(sim)
   sync all
endsubroutine
