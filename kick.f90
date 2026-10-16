subroutine kick
   use omp_lib
   use variables
   use pencil_fft
   use cubefft
   use power_nu
   implicit none

  logical,parameter :: PM3=.true.
  logical,parameter :: PP=.true.
   logical,parameter :: PP_corr=.true.
   logical,parameter :: speedtest=.false.
   integer inode,itile,np_phy,nl(3),nh(3),i1,i2,it,i3,isort(ns3),ires(ns3)
   real,allocatable :: rho_th(:,:,:),phi_th(:,:,:),force_th(:,:,:,:)
   real(8),allocatable :: ptotal(:),ttotal(:)
   real(8) pt
   real npc_max(ns3),overden_phy

   real rho_std(ngp*nnt,ngp*nnt,ngp*nnt)
   real rho_fin(ngp*nnt*nns,ngp*nnt*nns,ngp*nnt*nns)

   !nu
   if (Mass_nu > 0) then
      call tic(99)
      call get_tf_cb2matter
      call toc(99); if (head) print*,'  real time =  ',tcat(96:99,istep),'secs';
      if (head) print*,''
      if (head) print*,''
   endif
   ! call interp_Pk_CDM

   if (head) print*,'PM1' ! =====================================================
   call tic(11)
   iapm=1; pm%pm_layer=1; pm%iapm=iapm; pm%nwork=nt; pm%nstart(:)=1; pm%nend(:)=nc
   pm%tile1(:)=1; pm%tile2(:)=nnt; pm%nloop=1; pm%nex=2; pm%nphy=nt; pm%gridsize=ratio_cs
   pm%tile_shift=1; pm%utile_shift=0; pm%m1=1; pm%m2=nt; pm%nforce=nc; pm%m1phi(:)=-2; pm%m2phi(:)=nc+3
   pm%sigv1=sigma_vi; pm%sigv2=sigma_vi; pm%nc1(:)=1; pm%nc2(:)=nt; testrho=0; pm%f2max=0
   call particle_mesh(rho1)
   sim%dt_pm1=0.5*sqrt(1./(sqrt(pm%f2max)*a_mid*G_grid))
   call toc(11)
   if (head) then
      print*,'  pm%f2max =',pm%f2max; print*,'    sum of PM1 density =',testrho
      print*,'  real time =',tcat(11,istep),'secs'; print*,''
   endif
   sync all

   if (head) print*,'PM2' ! =====================================================
   call tic(12)
   vmax=0; f2max=0; pt=0
   allocate(ptotal(nnt**3),ttotal(nnt**3)); ptotal=0
   !$omp paralleldo num_threads(nteam) schedule(dynamic,1) default(shared) &
   !$omp& private(itile,iteam,pm,rho_th,force_th,vmax_team,f2max_team) reduction(max:vmax,f2max) reduction(+:pt)
   do itile=1,nnt**3
      iteam=omp_get_thread_num()+1; iapm=1; pm%pm_layer=2; pm%iapm=iapm; pm%nwork=ngt; pm%nstart(:)=1-ngb
      pm%nend(:)=ngp+ngb; pm%nend(1)=pm%nend(1)+2; pm%tile1(:)=ixyz2(:,itile); pm%tile2(:)=ixyz2(:,itile)
      pm%nloop=ncb; pm%nex=ngb+1; pm%nphy=ngp; pm%gridsize=1; pm%tile_shift=0; pm%utile_shift=0;
      pm%m1=1-ngb; pm%m2=ngp+ngb; pm%nforce=ngp; pm%m1phi(:)=1-ngb; pm%m2phi(:)=ngp+ngb; pm%m2phi(1)=pm%m2phi(1)+2
      pm%sigv1=sigma_vi; pm%sigv2=sigma_vi_new; pm%nc1(:)=1; pm%nc2(:)=nt
      allocate(rho_th(1-pm%nex:pm%nphy+pm%nex,1-pm%nex:pm%nphy+pm%nex,1-pm%nex:pm%nphy+pm%nex))
      allocate(force_th(ndim,0:pm%nforce+1,0:pm%nforce+1,0:pm%nforce+1))
      call pm_thread(rho2(:,:,:,iteam),rho_th,force_th,pm,itile,iteam,f2max_team,vmax_team)

      deallocate(rho_th,force_th)
      f2max=max(f2max,f2max_team); vmax=max(vmax,vmax_team)
   enddo
   !$omp endparalleldo
   sigma_vi=sigma_vi_new;
   sim%vsim2phys=(1.5/sim%a)*box*h0*100.*sqrt(omega_m)/ng_global
   sim%dt_pm2=0.5*sqrt( 1. / (sqrt(f2max)*a_mid*G_grid) )
   call toc(12)
   if (head) then
      print*,'  sum of PM2 density',sum(ptotal),pt; print*,'  vmax =',vmax; print*,'  f2max =',f2max
      !print*,'  team time =',ttotal
      print*,'  total team time =',sum(ttotal),'secs'; print*,'  real time =',tcat(12,istep),'secs'; print*,''
   endif
   deallocate(ptotal,ttotal)
   sync all
   ! print*,'  node',image,'     vmax =',vmax,'     f2max =',f2max
   ! sync all


   !       print*,'+++++++++++++++++++++++++++++++++++++++++++++'
   !       print*,'+++++++++++++++++++++++++++++++++++++++++++++'
   !       print*,'+++++++++++++++++++++++++++++++++++++++++++++'
   !       print*,'+++++++++++++++++++++++++++++++++++++++++++++'
   !       print*,sim%timestep

   !       write(str_z,'(f8.1)') Mass_nu
   !       open(11,file='./delat_tile_'//trim(adjustl(str_z))//'.bin',status='replace',access='stream')
   !       do itile=1,ngp*nnt
   !         ! print*,itile
   !         write(11) rho_std(:,:,itile) ! write layer by layer to avoid bug
   !       enddo
   !       close(11)
   !       open(11,file='./tf2.bin',status='replace',access='stream')
   !       do itile=1,ngt
   !         ! print*,itile
   !         write(11) tf2(:,:,itile) ! write layer by layer to avoid bug
   !       enddo
   !       close(11)




   ! error stop


   if (PM3 .or. PP) then  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      napm=0
      do itile=1,ns3 ! count clustering
         nl=(ixyz3(1:3,itile)-1)*ntt+1; nh=ixyz3(1:3,itile)*ntt
         i1=ixyz3(4,itile); i2=ixyz3(5,itile); i3=ixyz3(6,itile);
         np_phy=sum(rhoc(nl(1):nh(1),nl(2):nh(2),nl(3):nh(3),i1,i2,i3))
         overden_phy=(real(np_phy)/(ntt*np_nc)**3)-1
         npc_max(itile)=maxval(rhoc(nl(1):nh(1),nl(2):nh(2),nl(3):nh(3),i1,i2,i3))/64.
         if (npc_max(itile)<24) then
            ires(itile)=2
         elseif (npc_max(itile)<114) then
            ires(itile)=3
         elseif (npc_max(itile)<320) then
            ires(itile)=4
            ! elseif (npc_max(itile)<1350) then
            !    ires(itile)=5
            !elseif (npc_max(itile)<2000) then
            !  ires(itile)=5
            ! else
            !    ires(itile)=6
         else
            ires(itile)=5
         endif
         !ires(itile)=4
         napm(ires(itile))=napm(ires(itile))+1
      enddo
      call indexed_sort(ns3,-npc_max,isort)
      !print*,npc_max
      ! print*,isort
      ! print*,npc_max(isort)
   endif  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   sync all

   if (speedtest) then
      allocate(ptotal(ns3),ttotal(ns3)); ptotal=0
      do iapm=2,6
         if (head) print*,'speed test: iapm =',iapm
         !$omp paralleldo num_threads(nteam) schedule(dynamic,1)&
         !$omp& default(shared) private(it,iteam,itile,nl,nh,i1,i2,i3,pm,rho_th,force_th,vmax_team,f2max_team)&
         !$omp& reduction(max:vmax,f2max) reduction(+:pt)
         do it=1,ns3
            iteam=omp_get_thread_num()+1; itile=isort(it); nl=(ixyz3(1:3,itile)-1)*ntt+1; nh=ixyz3(1:3,itile)*ntt
            i1=ixyz3(4,itile); i2=ixyz3(5,itile); i3=ixyz3(6,itile); pm%pm_layer=3; pm%iapm=iapm; pm%nwork=nft(iapm)
            pm%nstart(:)=1-nfb(iapm); pm%nend(:)=nfp(iapm)+nfb(iapm); pm%nend(1)=pm%nend(1)+2
            pm%tile1(:)=ixyz3(4:6,itile); pm%tile2(:)=ixyz3(4:6,itile); pm%nloop=1; pm%nex=nfb(iapm)+1; pm%nphy=nfp(iapm);
            pm%gridsize=1./ratio_sf(iapm); pm%tile_shift=0; pm%utile_shift=ixyz3(1:3,itile)-1; pm%m1=1-nfb(iapm);
            pm%m2=nfp(iapm)+nfb(iapm);pm%nforce=nfp(iapm); pm%m1phi(:)=1-nfb(iapm); pm%m2phi(:)=nfp(iapm)+nfb(iapm);
            pm%m2phi(1)=pm%m2phi(1)+2; pm%sigv1=sigma_vi; pm%sigv2=sigma_vi
            pm%nc1(:)=(ixyz3(1:3,itile)-1)*ntt+1; pm%nc2(:)=ixyz3(1:3,itile)*ntt
            allocate(rho_th(1-pm%nex:pm%nphy+pm%nex,1-pm%nex:pm%nphy+pm%nex,1-pm%nex:pm%nphy+pm%nex))
            allocate(force_th(ndim,0:pm%nforce+1,0:pm%nforce+1,0:pm%nforce+1))
            selectcase(iapm)
             case(2)
               call pm_thread(rho3_2(:,:,:,iteam),rho_th,force_th,pm,itile,iteam,f2max_team,vmax_team)
             case(3)
               call pm_thread(rho3_4(:,:,:,iteam),rho_th,force_th,pm,itile,iteam,f2max_team,vmax_team)
             case(4)
               call pm_thread(rho3_6(:,:,:,iteam),rho_th,force_th,pm,itile,iteam,f2max_team,vmax_team)
             case(5)
               call pm_thread(rho3_8(:,:,:,iteam),rho_th,force_th,pm,itile,iteam,f2max_team,vmax_team)
            endselect
            deallocate(rho_th,force_th)
            speed(1,iapm,itile)=npc_max(itile)
            speed(3,iapm,itile)=ttotal(itile)
         enddo
         !$omp paralleldo num_threads(nteam) schedule(dynamic,1)&
         !$omp& default(shared) private(it,iteam,itile,i1,i2,i3,pm,vmax_team,f2max_team)&
         !$omp& reduction(max:vmax,f2max) reduction(+:pt)
         do it=1,ns3
            iteam=omp_get_thread_num()+1; itile=isort(it); i1=ixyz3(4,itile); i2=ixyz3(5,itile); i3=ixyz3(6,itile);
            pm%pm_layer=3; pm%iapm=iapm; pm%nwork=nft(iapm); pm%nstart(:)=1-nfb(iapm); pm%nend(:)=nfp(iapm)+nfb(iapm); pm%nend(1)=pm%nend(1)+2
            pm%tile1(:)=ixyz3(4:6,itile); pm%tile2(:)=ixyz3(4:6,itile);pm%nloop=1; pm%nex=nfb(iapm)+1; pm%nphy=nfp(iapm);
            pm%gridsize=1./ratio_sf(iapm); pm%tile_shift=0; pm%utile_shift=ixyz3(1:3,itile)-1; pm%m1=1-nfb(iapm); pm%m2=nfp(iapm)+nfb(iapm);pm%nforce=nfp(iapm); pm%m1phi(:)=1-nfb(iapm); pm%m2phi(:)=nfp(iapm)+nfb(iapm); pm%m2phi(1)=pm%m2phi(1)+2
            pm%sigv1=sigma_vi; pm%sigv2=sigma_vi; pm%nc1(:)=(ixyz3(1:3,itile)-1)*ntt+1; pm%nc2(:)=ixyz3(1:3,itile)*ntt
            call PP_Force(pm,itile,iteam,f2max_team,vmax_team)
            speed(4,iapm,itile)=ttotal(itile)
         enddo
         !$omp endparalleldo
      enddo
      deallocate(ptotal,ttotal)
      sync all
      if (head) then
         do inode=1,nn**3
            speed_global(:,:,:,inode)=speed(:,:,:)[inode]
         enddo
         sim%cur_checkpoint=sim%cur_checkpoint-1
         open(102,file=output_name('speedtest'),status='replace',access='stream')
         write(102) speed_global
         close(102)
         print*,'speed test complete.'
      endif
      sync all
      error stop
   endif

   if (PM3) then ! =====================================================
      if (head) print*, 'PM3'
      call tic(13)
      allocate(ptotal(ns3),ttotal(ns3)); ptotal=0
      pm%f2max=0; f2_max_pp=0; vmax=0; f2max=0; pt=0;
      !$omp paralleldo num_threads(nteam) schedule(dynamic,1)&
      !$omp& default(shared) private(it,iteam,itile,iapm,nl,nh,i1,i2,i3,pm,rho_th,force_th,vmax_team,f2max_team)&
      !$omp& reduction(max:vmax,f2max) reduction(+:pt)
      do it=1,ns3
         iteam=omp_get_thread_num()+1; itile=isort(it); iapm=ires(itile); nl=(ixyz3(1:3,itile)-1)*ntt+1; nh=ixyz3(1:3,itile)*ntt
         i1=ixyz3(4,itile); i2=ixyz3(5,itile); i3=ixyz3(6,itile); pm%pm_layer=3; pm%iapm=iapm; pm%nwork=nft(iapm)
         pm%nstart(:)=1-nfb(iapm); pm%nend(:)=nfp(iapm)+nfb(iapm); pm%nend(1)=pm%nend(1)+2
         pm%tile1(:)=ixyz3(4:6,itile); pm%tile2(:)=ixyz3(4:6,itile);pm%nloop=1; pm%nex=nfb(iapm)+1; pm%nphy=nfp(iapm);
         pm%gridsize=1./ratio_sf(iapm); pm%tile_shift=0; pm%utile_shift=ixyz3(1:3,itile)-1; pm%m1=1-nfb(iapm);
         pm%m2=nfp(iapm)+nfb(iapm);pm%nforce=nfp(iapm); pm%m1phi(:)=1-nfb(iapm); pm%m2phi(:)=nfp(iapm)+nfb(iapm);
         pm%m2phi(1)=pm%m2phi(1)+2; pm%sigv1=sigma_vi; pm%sigv2=sigma_vi
         pm%nc1(:)=(ixyz3(1:3,itile)-1)*ntt+1; pm%nc2(:)=ixyz3(1:3,itile)*ntt
         allocate(rho_th(1-pm%nex:pm%nphy+pm%nex,1-pm%nex:pm%nphy+pm%nex,1-pm%nex:pm%nphy+pm%nex))
         allocate(force_th(ndim,0:pm%nforce+1,0:pm%nforce+1,0:pm%nforce+1))
         selectcase(iapm)
          case(2)
            call pm_thread(rho3_2(:,:,:,iteam),rho_th,force_th,pm,itile,iteam,f2max_team,vmax_team)
          case(3)
            call pm_thread(rho3_4(:,:,:,iteam),rho_th,force_th,pm,itile,iteam,f2max_team,vmax_team)
          case(4)
            call pm_thread(rho3_6(:,:,:,iteam),rho_th,force_th,pm,itile,iteam,f2max_team,vmax_team)
          case(5)
            call pm_thread(rho3_8(:,:,:,iteam),rho_th,force_th,pm,itile,iteam,f2max_team,vmax_team)
            !  case(6)
            !    call pm_thread(rho3_12(:,:,:,iteam),rho_th,force_th,pm,itile,iteam,f2max_team,vmax_team)
         endselect
         deallocate(rho_th,force_th)
         f2max=max(f2max,f2max_team); vmax=max(vmax,vmax_team)
         ! if (image == 6) print*, 'f2max,vmax:', f2max, vmax,iapm,npc_max(itile)
      enddo
      !$omp endparalleldo
      sim%dt_pm3=0.5*sqrt( 1. / (sqrt(f2max)*a_mid*G_grid) )
      call toc(13)
      if (head) then
         print*,'  sum of PM3 density',sum(ptotal); print*,'  vmax =',vmax; print*,'  f2max =',f2max
         print*,'  total team time =',sum(ttotal),'secs'; print*,'  real time =',tcat(13,istep),'secs'; print*,''
      endif
      deallocate(ptotal,ttotal)
      sync all
      ! print '( "PM3 node ",I0," vmax=",3(F10.3,:,1X)," f2max=",E8.1 )', image, vmax, f2max
      ! stop


      !       print*,'+++++++++++++++++++++++++++++++++++++++++++++'xx
      !       print*,'+++++++++++++++++++++++++++++++++++++++++++++'
      !       print*,'+++++++++++++++++++++++++++++++++++++++++++++'
      !       print*,'+++++++++++++++++++++++++++++++++++++++++++++'
      !       print*,sim%timestep

      !       write(str_z,'(f8.1)') Mass_nu
      !       open(11,file='./delat_tile_'//trim(adjustl(str_z))//'.bin',status='replace',access='stream')
      !       do itile=1,ngp*nnt
      !         ! print*,itile
      !         write(11) rho_std(:,:,itile) ! write layer by layer to avoid bug
      !       enddo
      !       close(11)
      !       open(11,file='./tf2.bin',status='replace',access='stream')
      !       do itile=1,ngt
      !         ! print*,itile
      !         write(11) tf2(:,:,itile) ! write layer by layer to avoid bug
      !       enddo
      !       close(11)




      ! stop
      ! if (maxval(vmax) > 100) error stop 'vmax > 100'
   endif !PM3

   sync all

   if (PP) then ! =====================================================
      if (head) print*,'PP'
      call tic(14)
      allocate(ptotal(ns3),ttotal(ns3))
      vmax=0; f2max=0; pt=0
      !$omp paralleldo num_threads(ncore) schedule(dynamic,1)&
      !$omp& default(shared) private(it,iteam,itile,iapm,i1,i2,i3,pm,vmax_team,f2max_team)&
      !$omp& reduction(max:vmax,f2max) reduction(+:pt)
      do it=1,ns3
         iteam=omp_get_thread_num()+1; itile=isort(it); iapm=ires(itile);i1=ixyz3(4,itile); i2=ixyz3(5,itile); i3=ixyz3(6,itile);
         pm%pm_layer=3; pm%iapm=iapm; pm%nwork=nft(iapm); pm%nstart(:)=1-nfb(iapm); pm%nend(:)=nfp(iapm)+nfb(iapm); pm%nend(1)=pm%nend(1)+2
         pm%tile1(:)=ixyz3(4:6,itile); pm%tile2(:)=ixyz3(4:6,itile);pm%nloop=1; pm%nex=nfb(iapm)+1; pm%nphy=nfp(iapm);
         pm%gridsize=1./ratio_sf(iapm); pm%tile_shift=0; pm%utile_shift=ixyz3(1:3,itile)-1; pm%m1=1-nfb(iapm); pm%m2=nfp(iapm)+nfb(iapm);pm%nforce=nfp(iapm); pm%m1phi(:)=1-nfb(iapm); pm%m2phi(:)=nfp(iapm)+nfb(iapm); pm%m2phi(1)=pm%m2phi(1)+2
         pm%sigv1=sigma_vi; pm%sigv2=sigma_vi; pm%nc1(:)=(ixyz3(1:3,itile)-1)*ntt+1; pm%nc2(:)=ixyz3(1:3,itile)*ntt
         call PP_Force(pm,itile,iteam,f2max_team,vmax_team)
         f2max=max(f2max,f2max_team); vmax=max(vmax,vmax_team)
      enddo
      !$omp endparalleldo
      sim%dt_pp=0.05*sqrt( 1. / (sqrt(f2max)*a_mid*G_grid) )
      call toc(14)
      if (head) then
         print*,'  sum of PP count',sum(ptotal); print*,'  vmax =',vmax; print*,'  f2max =',f2max
         if(PP_corr .and. Mass_nu > 0) print*,'  PP_corr =',((1-f_nu)**2)
         print*,'  total team time =',sum(ttotal),'secs'; print*,'  real time =',tcat(14,istep),'secs'; print*,''
         ! do inode=2,nn**3
         !    print*,'  node',inode; print*,'     vmax =',vmax(:)[inode]; print*,'     f2max =',f2max[inode]
         ! enddo
      endif
      deallocate(ptotal,ttotal)
      sync all
      ! print '( "PP node ",I0," vmax=",3(F10.3,:,1X)," f2max=",E8.1 )', image, vmax, f2max
   endif ! PP
   sync all

   sim%dt_vmax=vbuf*14/maxval(vmax)
   sim%vz_max=vmax(3)
   if (head) then ! gather info from all nodes
      inode = 1
      ! print*,'  node',1
      ! print*,'    dt_vmax=',sim[inode]%dt_vmax
      ! print*,'    dt_pm1 =',sim[inode]%dt_pm1
      ! print*,'    dt_pm2 =',sim[inode]%dt_pm2
      ! print*,'    dt_pm3 =',sim[inode]%dt_pm3
      ! print*,'    dt_pp  =',sim[inode]%dt_pp
      do inode=2,nn**3
         napm(:)=napm(:)+napm(:)[inode]
         sim%dt_pm1=min(sim%dt_pm1,sim[inode]%dt_pm1)
         sim%dt_pm2=min(sim%dt_pm2,sim[inode]%dt_pm2)
         sim%dt_pm3=min(sim%dt_pm3,sim[inode]%dt_pm3)
         sim%dt_pp=min(sim%dt_pp,sim[inode]%dt_pp)
         sim%dt_vmax=min(sim%dt_vmax,sim[inode]%dt_vmax)
         ! print*,'  node',inode
         ! print*,'    dt_vmax=',sim[inode]%dt_vmax
         ! print*,'    dt_pm1 =',sim[inode]%dt_pm1
         ! print*,'    dt_pm2 =',sim[inode]%dt_pm2
         ! print*,'    dt_pm3 =',sim[inode]%dt_pm3
         ! print*,'    dt_pp  =',sim[inode]%dt_pp
      enddo
   endif
   sync all ! broadcast
   napm(:)=napm(:)[1]
   sim%dt_pm1=sim[1]%dt_pm1
   sim%dt_pm2=sim[1]%dt_pm2
   sim%dt_pm3=sim[1]%dt_pm3
   sim%dt_pp=sim[1]%dt_pp
   sim%dt_vmax=sim[1]%dt_vmax
   sync all

   if (head) then
      tcat(61:67,istep)=napm
      print*,'    napm =',napm
      print*,'  real time =',tcat(13,istep),real(t2-t1)/t_rate,'secs'
   endif
   call toc(5)

contains




   pure elemental real function F_ra_inline(r,apm)
      real,intent(in) :: r,apm
      real ep
      F_ra_inline=0
      ep=2*r/apm
      if (apm==0 .or. ep>2) then
         F_ra_inline=r**(-2)
      elseif (ep>1) then
         F_ra_inline=(1./35./apm**2)*(12*ep**(-2) - 224 + 896*ep - 840*ep**2 + 224*ep**3 + 70*ep**4 - 48*ep**5 + 7*ep**6)
      elseif (ep>0) then
         F_ra_inline=(1./35./apm**2)*(224*ep - 224*ep**3 + 70*ep**4 + 48*ep**5 - 21*ep**6)
      else
         F_ra_inline=0
      endif
   endfunction

   subroutine PP_Force(pm,itile,iteam,f2max_t,vmax_t)
      use omp_lib
      implicit none
      type(type_pm) pm
      integer(4) i_neighbor,nptile,itile,iteam,iapm,rcp,npgrid,ncount,tteam1,tteam2,tteam_rate
      integer i,j,k,np,l,idx(3),ip1,ip2,pos,cnt,i_act,s1,e1,s2,e2
      integer n_phy_cnt,n_ext_cnt,n_act,n_act_ext
      integer,allocatable :: cell_start(:,:,:),cell_count(:,:,:),cell_offset(:,:,:),p_cell(:,:)
      integer,allocatable :: active_ijk(:,:),active_ext_ijk(:,:)
      integer(8) ip,nzero
      integer(8),allocatable :: ip_tmp(:),ip_local(:)
      real xvec(3),rvec(3),rmag,r2,fpp(3),f1(3),pp_range,pp_range2,f2max_t,vmax_t(:),nbb,mass_scale,v_scale
      real,allocatable :: xf_tmp(:,:),xf(:,:),af(:,:)

      call system_clock(tteam1,tteam_rate)
      iapm=pm%iapm
      pp_range=apm3(iapm); pp_range2=pp_range**2; rcp=ratio_sf(iapm); npgrid=ntt*rcp
      f2max_t=0; vmax_t=0; nbb = npgrid+rcp

      nptile=sum(rhoc(pm%nc1(1)-1:pm%nc2(1)+1,pm%nc1(2)-1:pm%nc2(2)+1,pm%nc1(3)-1:pm%nc2(3)+1, &
                      pm%tile1(1),pm%tile1(2),pm%tile1(3)))
      if (nptile == 0) then
         ptotal(itile)=0
         call system_clock(tteam2,tteam_rate)
         ttotal(itile)=real(tteam2-tteam1)/tteam_rate
         return
      endif

      allocate(cell_start(0:npgrid+1,0:npgrid+1,0:npgrid+1))
      allocate(cell_count(0:npgrid+1,0:npgrid+1,0:npgrid+1))
      allocate(cell_offset(0:npgrid+1,0:npgrid+1,0:npgrid+1))
      allocate(p_cell(3,nptile),xf_tmp(3,nptile),ip_tmp(nptile))

      cell_count=0; ip1=0; n_phy_cnt=0; n_ext_cnt=0

      ! Pass 1: compute fine cell indices and count particles in [0..npgrid+1]^3
      do k=pm%nc1(3)-1,pm%nc2(3)+1
         do j=pm%nc1(2)-1,pm%nc2(2)+1
            do i=pm%nc1(1)-1,pm%nc2(1)+1
               np=rhoc(i,j,k,pm%tile1(1),pm%tile1(2),pm%tile1(3))
               nzero=idx_b_r(j,k,pm%tile1(1),pm%tile1(2),pm%tile1(3))-sum(rhoc(i:,j,k,pm%tile1(1),pm%tile1(2),pm%tile1(3)))
               do l=1,np
                  ip=nzero+l; ip1=ip1+1; ip_tmp(ip1)=ip
#ifdef ZIPX
                  xvec=[i,j,k]-1+(int(xp(:,ip)+ishift,izipx)+rshift)*x_resolution
#else
                  xvec = xp(:,ip)/ratio_cs - (pm%tile1(:)-1)*nt
#endif
                  xf_tmp(:,ip1)=xvec*ratio_cs
                  xvec = rcp*(xvec-pm%utile_shift*ntt)
#ifndef ZIPX
                  if ( xvec(1) == nbb ) xvec(1)= 1-rcp
                  if ( xvec(2) == nbb ) xvec(2)= 1-rcp
                  if ( xvec(3) == nbb ) xvec(3)= 1-rcp
#endif
                  idx=floor(xvec)+1
                  if (minval(idx)<1-rcp .or. maxval(idx)>npgrid+rcp) then
                     print*,ip,'out of bound',1-rcp,npgrid+rcp,npgrid,nbb
                     print*,ip,'ijk',i,j,k
                     print*,ip,'xp',xp(:,ip)
                     print*,ip,'xvec',xp(:,ip)/ratio_cs - (pm%tile1(:)-1)*nt
                     print*,ip,'xf',xf_tmp(:,ip1)
                     print*,ip,'xfg',xp(:,ip)/ratio_cs - (pm%tile1(:)-1)*nt-pm%utile_shift*ntt
                     print*,ip,'xff',xvec
                     print*,ip,'idx',idx
                     print*,ip,'tile',pm%tile1,pm%utile_shift
                     print*,ip,'nc1',pm%nc1
                     print*,ip,'nc2',pm%nc2
                     error stop
                  endif
                  if (minval(idx)>=0 .and. maxval(idx)<=npgrid+1) then
                     p_cell(:,ip1)=idx
                     cell_count(idx(1),idx(2),idx(3))=cell_count(idx(1),idx(2),idx(3))+1
                     n_ext_cnt=n_ext_cnt+1
                     if (minval(idx)>=1 .and. maxval(idx)<=npgrid) then
                        n_phy_cnt=n_phy_cnt+1
                     endif
                  else
                     p_cell(1,ip1)=-999
                  endif
               enddo
            enddo
         enddo
      enddo

      allocate(active_ijk(3,max(1,n_phy_cnt)))
      allocate(active_ext_ijk(3,max(1,n_ext_cnt)))
      allocate(xf(3,max(1,n_ext_cnt)))
      allocate(af(3,max(1,n_ext_cnt)))
      allocate(ip_local(max(1,n_phy_cnt)))
      af=0.0

      ! Pass 2: prefix sum over fine cells (physical 1..npgrid first, then ghost cells in 0..npgrid+1)
      pos=1; n_act=0; n_act_ext=0
      do k=1,npgrid
         do j=1,npgrid
            do i=1,npgrid
               cnt=cell_count(i,j,k)
               cell_start(i,j,k)=pos
               cell_offset(i,j,k)=pos
               pos=pos+cnt
               if (cnt>0) then
                  n_act=n_act+1
                  active_ijk(:,n_act)=[i,j,k]
               endif
            enddo
         enddo
      enddo

      do k=0,npgrid+1
         do j=0,npgrid+1
            do i=0,npgrid+1
               cnt=cell_count(i,j,k)
               if (cnt>0) then
                  n_act_ext=n_act_ext+1
                  active_ext_ijk(:,n_act_ext)=[i,j,k]
               endif
               if (i>=1 .and. i<=npgrid .and. j>=1 .and. j<=npgrid .and. k>=1 .and. k<=npgrid) cycle
               cell_start(i,j,k)=pos
               cell_offset(i,j,k)=pos
               pos=pos+cnt
            enddo
         enddo
      enddo

      ! Pass 3: scatter particles into contiguous CSR order
      do l=1,ip1
         if (p_cell(1,l)/=-999) then
            pos=cell_offset(p_cell(1,l),p_cell(2,l),p_cell(3,l))
            xf(:,pos)=xf_tmp(:,l)
            if (pos<=n_phy_cnt) then
               ip_local(pos)=ip_tmp(l)
            endif
            cell_offset(p_cell(1,l),p_cell(2,l),p_cell(3,l))=pos+1
         endif
      enddo

      deallocate(p_cell,cell_offset,xf_tmp,ip_tmp)

      !self pp
      do i_act=1,n_act
         i=active_ijk(1,i_act); j=active_ijk(2,i_act); k=active_ijk(3,i_act)
         s1=cell_start(i,j,k)
         e1=s1+cell_count(i,j,k)-1
         do ip1=s1,e1-1
            f1=af(:,ip1)
            do ip2=ip1+1,e1
               rvec=xf(:,ip2)-xf(:,ip1)
               r2=rvec(1)*rvec(1)+rvec(2)*rvec(2)+rvec(3)*rvec(3)
               if (0.<r2 .and. r2<pp_range2) then
                  rmag=sqrt(r2)
                  fpp=rvec/rmag*(F_ra_inline(rmag,app)-F_ra_inline(rmag,pp_range))
                  f1=f1+fpp
                  af(:,ip2)=af(:,ip2)-fpp
               endif
            enddo
            af(:,ip1)=f1
         enddo
      enddo

      !neigh pp: include active ghost cells in 0..npgrid+1
      do i_act=1,n_act_ext
         i=active_ext_ijk(1,i_act); j=active_ext_ijk(2,i_act); k=active_ext_ijk(3,i_act)
         s1=cell_start(i,j,k)
         e1=s1+cell_count(i,j,k)-1
         do i_neighbor=1,n_neighbor
            idx=[i,j,k]+ijk(:,i_neighbor)
            if (minval(idx)<0 .or. maxval(idx)>npgrid+1) cycle
            ! Neither buffer cell contributes to a physical particle update.
            if (min(i,j,k)==0 .or. max(i,j,k)==npgrid+1) then
               if (minval(idx)==0 .or. maxval(idx)==npgrid+1) cycle
            endif
            if (cell_count(idx(1),idx(2),idx(3))==0) cycle
            s2=cell_start(idx(1),idx(2),idx(3))
            e2=s2+cell_count(idx(1),idx(2),idx(3))-1
            do ip1=s1,e1
               f1=af(:,ip1)
               do ip2=s2,e2
                  rvec=xf(:,ip2)-xf(:,ip1)
                  r2=rvec(1)*rvec(1)+rvec(2)*rvec(2)+rvec(3)*rvec(3)
                  if (r2<pp_range2) then
                     rmag=sqrt(r2)
                     fpp=rvec/rmag*(F_ra_inline(rmag,app)-F_ra_inline(rmag,pp_range))
                     f1=f1+fpp
                     af(:,ip2)=af(:,ip2)-fpp
                  endif
               enddo
               af(:,ip1)=f1
            enddo
         enddo
      enddo

      !update v: physical particles occupy the contiguous prefix
      if (PP_corr .and. Mass_nu > 0) then
         mass_scale=sim%mass_p_cdm*((1-f_nu)**2)
      else
         mass_scale=sim%mass_p_cdm
      endif
      v_scale=a_mid*dt/6/pi
      ncount=n_phy_cnt
      do ip1=1,ncount
         f1=af(:,ip1)*mass_scale
         f2max_t=max(f2max_t,sum(f1**2))
         ip=ip_local(ip1)
         vp(:,ip)=vp(:,ip)+f1*v_scale
         vmax_t=max(vmax_t,abs(vp(:,ip)))
      enddo

      deallocate(cell_start,cell_count,active_ijk,active_ext_ijk,xf,af,ip_local)
      ptotal(itile)=ncount
      call system_clock(tteam2,tteam_rate)
      ttotal(itile)=real(tteam2-tteam1)/tteam_rate
   endsubroutine PP_Force



   !! call pm_thread(rho2(:,:,:,iteam),rho_th,force_th,pm,itile,iteam)
   subroutine pm_thread(density,rho_th,force_th,pm,itile,iteam,f2max_t,vmax_t)
      use omp_lib
      use ieee_arithmetic
      implicit none
      integer,parameter :: nlayer=3 ! thread save for TSC intepolation
      type(type_pm) pm
      real density(pm%nstart(1):,pm%nstart(2):,pm%nstart(3):)
      !complex density_k(pm%nwork/2+1,pm%nwork,pm%nwork)
      real rho_th(1-pm%nex:,1-pm%nex:,1-pm%nex:),force_th(1:,0:,0:,0:)
      integer(4) itile,iteam,iapm,i_0,j_0,k_0,i_n(4),j_n(4),k_n(4),tteam1,tteam2,tteam_rate
      integer i,j,k,l,np,idx0(ndim,p+1),ilayer,i1,i2,i3,ix,iy,iz,n1(ndim),n2(ndim)
      integer(8) ip,nzero
      real(8) xpos(ndim),nbb
      real dx(ndim,p+1),l3(ndim),dv(ndim),vreal(ndim),f2max_t,vmax_t(:)
      !equivalence(density,density_k)

      call system_clock(tteam1,tteam_rate)
      iapm=pm%iapm; density=0; nbb=pm%nphy+pm%nex-1
      do iz=pm%tile1(3),pm%tile2(3) ! interpolate particles to density field
         do iy=pm%tile1(2),pm%tile2(2)
            do ix=pm%tile1(1),pm%tile2(1)
               !print*,ix,iy,iz
               rho_th=0
               do ilayer=0,nlayer-1
                  !$omp paralleldo default(shared) num_threads(nnest) schedule(dynamic,1)&
                  !$omp& private(k,j,i,np,nzero,l,ip,xpos,idx0,l3,dx,i3,i2,i1)
                  do k=pm%nc1(3)-pm%nloop+ilayer,pm%nc2(3)+pm%nloop,nlayer
                     !do k=pm%nc1(3)-pm%nloop,pm%nc2(3)+pm%nloop
                     do j=pm%nc1(2)-pm%nloop,pm%nc2(2)+pm%nloop
                        do i=pm%nc1(1)-pm%nloop,pm%nc2(1)+pm%nloop
                           np=rhoc(i,j,k,ix,iy,iz)
                           nzero=idx_b_r(j,k,ix,iy,iz)-sum(rhoc(i:,j,k,ix,iy,iz))
                           do l=1,np
                              ip=nzero+l
#ifdef ZIPX
                              xpos=([i,j,k]-1)+(int(xp(:,ip)+ishift,izipx)+rshift)*x_resolution
                              xpos=xpos*ratio_cs/pm%gridsize-pm%utile_shift*nfp(iapm)
#else
                              xpos = (xp(:,ip)-([ix,iy,iz]-1)*nt*ratio_cs)   &
                                 /pm%gridsize-pm%utile_shift*nfp(iapm)

                              if ( xpos(1) == nbb ) xpos(1)= 1-pm%nex
                              if ( xpos(2) == nbb ) xpos(2)= 1-pm%nex
                              if ( xpos(3) == nbb ) xpos(3)= 1-pm%nex

                              ! xpos = xp(:,ip)-([ix,iy,iz]-1)*nt*ratio_cs
                              ! xpos = xpos/pm%gridsize-pm%utile_shift*nfp(iapm)
                              ! ! if (minval(xpos) < -pm%nex .or. maxval(xpos) > pm%nphy+pm%nex) then
                              ! !    print*,'particle out of range'
                              ! !    print*,'tile',ix,iy,iz,i,j,k,l
                              ! !    print*,'xp',xp(:,ip)
                              ! !    print*,'xpos',xpos
                              ! !    error stop
                              ! ! endif
                              ! xpos = modulo(xpos+spread(real(pm%nex),dim=1,ncopies=3),real(pm%nphy+2*pm%nex))-spread(real(pm%nex),dim=1,ncopies=3)
#endif
                              idx0(:,2)=floor(xpos)+1
                              idx0(:,1)=idx0(:,2)-1
                              idx0(:,3)=idx0(:,2)+1
                              l3=xpos-floor(xpos)
                              dx(:,1)=(1-l3)**2/2
                              dx(:,3)=l3**2/2
                              dx(:,2)=1-dx(:,1)-dx(:,3)
                              if (minval(idx0)<1-pm%nex .or. maxval(idx0)>pm%nphy+pm%nex) then
                                 print*,'mass assignment out of range'
                                 print*,ip,ix,iy,iz,i,j,k
                                 print*,ip,pm%nex,pm%nphy
                                 print*,ip,pm%utile_shift, nfp(iapm)
                                 print*,ip, xp(:,ip)
#ifdef ZIPX
                                 print*,ip,  'tile shift',([i,j,k]-1)+(int(xp(:,ip)+ishift,izipx)+rshift)*x_resolution
#endif
                                 print*,ip,'xpos',xpos
                                 print*,ip,'idx0',idx0
                                 error stop
                              endif
                              do i3=1,p+1
                                 do i2=1,p+1
                                    do i1=1,p+1
                                       rho_th(idx0(1,i1),idx0(2,i2),idx0(3,i3))=rho_th(idx0(1,i1),idx0(2,i2),idx0(3,i3))+dx(1,i1)*dx(2,i2)*dx(3,i3)
                                    enddo
                                 enddo
                              enddo
                           enddo
                        enddo
                     enddo
                  enddo
                  !$omp endparalleldo
               enddo ! ilayer
               !print*,'    from rho_ex',pm%m1,pm%m2
               n1(:)=pm%nstart+pm%tile_shift*([ix,iy,iz]-1)*pm%nphy
               n2(:)=n1(:)+pm%nwork-1
               !print*,'     to density (',n1(1),':',n2(1),',',n1(2),':',n2(2),',',n1(3),':',n2(3),')'
               density(n1(1):n2(1),n1(2):n2(2),n1(3):n2(3))=rho_th(pm%m1:pm%m2,pm%m1:pm%m2,pm%m1:pm%m2)*sim%mass_p_cdm
               ptotal(itile)=ptotal(itile)+sum(1d0*rho_th(1:pm%nphy,1:pm%nphy,1:pm%nphy))
               pt=pt+sum(1d0*rho_th(1:pm%nphy,1:pm%nphy,1:pm%nphy))
            enddo
         enddo
      enddo

      selectcase(pm%pm_layer)
       case(0)
       case(1)
       case(2)
         call sfftw_execute(plan2(iteam))
         rho2k(:,:,:,iteam)=rho2k(:,:,:,iteam)*tf2
         call sfftw_execute(iplan2(iteam))
       case(3)
         call sfftw_execute(plan3(iteam,iapm))
         selectcase(iapm)
          case(2)
            rho3k_2(:,:,:,iteam) =rho3k_2(:,:,:,iteam)*tf3_2
          case(3)
            rho3k_4(:,:,:,iteam) =rho3k_4(:,:,:,iteam)*tf3_4
          case(4)
            rho3k_6(:,:,:,iteam) =rho3k_6(:,:,:,iteam)*tf3_6
          case(5)
            rho3k_8(:,:,:,iteam) =rho3k_8(:,:,:,iteam)*tf3_8
            !  case(6)
            !    rho3k_12(:,:,:,iteam)=rho3k_12(:,:,:,iteam)*tf3_12
         endselect
         call sfftw_execute(iplan3(iteam,iapm))
      endselect
      call potential_to_force(force_th,density,pm)
      f2max_t=maxval(sum(force_th(:,1:pm%nforce,1:pm%nforce,1:pm%nforce)**2,1))
      vmax_t=0

      ! if ( abs(f2max_t) > 1.0e11 .and. ieee_is_finite(f2max_t)) then
      !    open(11,file=output_name('error_phi'),status='replace',access='stream')
      !    write(11) density
      !    close(11)

      !    call sfftw_execute(plan3(iteam,iapm))
      !    selectcase(iapm)
      !     case(2)
      !       rho3k_2(:,:,:,iteam) =rho3k_2(:,:,:,iteam)/tf3_2
      !     case(3)
      !       rho3k_4(:,:,:,iteam) =rho3k_4(:,:,:,iteam)/tf3_4
      !     case(4)
      !       rho3k_6(:,:,:,iteam) =rho3k_6(:,:,:,iteam)/tf3_6
      !     case(5)
      !       rho3k_8(:,:,:,iteam) =rho3k_8(:,:,:,iteam)/tf3_8
      !     case(6)
      !       rho3k_12(:,:,:,iteam)=rho3k_12(:,:,:,iteam)/tf3_12
      !    endselect
      !    call sfftw_execute(iplan3(iteam,iapm))

      !    open(11,file=output_name('error_rho3'),status='replace',access='stream')
      !    write(11) density
      !    close(11)
      !    error stop
      ! endif
      nbb = pm%nforce
      do iz=pm%tile1(3),pm%tile2(3) ! update velocity
         do iy=pm%tile1(2),pm%tile2(2)
            do ix=pm%tile1(1),pm%tile2(1)
               !$omp paralleldo num_threads(nnest) default(shared) schedule(dynamic)&
               !$omp& private(k,j,i,np,nzero,l,ip,xpos,idx0,l3,dx,vreal,dv,i3,i2,i1)&
               !$omp& reduction(max:vmax_t)
               do k=pm%nc1(3),pm%nc2(3)
                  do j=pm%nc1(2),pm%nc2(2)
                     do i=pm%nc1(1),pm%nc2(1)
                        np=rhoc(i,j,k,ix,iy,iz)
                        nzero=idx_b_r(j,k,ix,iy,iz)-sum(rhoc(i:,j,k,ix,iy,iz))
                        do l=1,np ! loop over particle
                           ip=nzero+l
#ifdef ZIPX
                           xpos=pm%tile_shift*([ix,iy,iz]-1)*nt+([i,j,k]-1)+(int(xp(:,ip)+ishift,izipx)+rshift)*x_resolution
                           xpos=xpos*ratio_cs/pm%gridsize-pm%utile_shift*nfp(iapm)
#else
                           xpos = (xp(:,ip)-([ix,iy,iz]-1) * nt * ratio_cs) &
                              /pm%gridsize-pm%utile_shift*nfp(iapm)

                           if ( xpos(1) == nbb ) xpos(1)= 0
                           if ( xpos(2) == nbb ) xpos(2)= 0
                           if ( xpos(3) == nbb ) xpos(3)= 0
#endif
                           idx0(:,2)=floor(xpos)+1
                           idx0(:,1)=idx0(:,2)-1
                           idx0(:,3)=idx0(:,2)+1
                           l3=xpos-floor(xpos)
                           dx(:,1)=(1-l3)**2/2
                           dx(:,3)=l3**2/2
                           dx(:,2)=1-dx(:,1)-dx(:,3)
                           !vreal=tan(pi*real(vp(:,ip))/real(nvbin-1))/(sqrt(pi/2)/(pm%sigv1*vrel_boost))
                           dv=0
                           if (minval(idx0)<0 .or. maxval(idx0)>pm%nforce+1) then
                              print*,'update velocity out of range'
                              print*,'tile',ix,iy,iz,i,j,k,l
                              print*,'nforce',pm%nforce
                              print*,'xp',xp(:,ip)
                              print*,'xpos',xpos
                              print*,'dx',xpos-floor(xpos)
                              print*,'idx',idx0
                              error stop
                           endif
                           do i3=1,p+1
                              do i2=1,p+1
                                 do i1=1,p+1
                                    dv=dv+force_th(:,idx0(1,i1),idx0(2,i2),idx0(3,i3))*dx(1,i1)*dx(2,i2)*dx(3,i3)
                                 enddo
                              enddo
                           enddo
                           !vreal=vreal+dv*a_mid*dt/6/pi
                           vp(:,ip)=vp(:,ip)+dv*a_mid*dt/6/pi
                           !vmax_t=max(vmax_t,abs(vreal+vfield(:,i,j,k,ix,iy,iz)))
                           vmax_t=max(vmax_t,abs(vp(:,ip)))
                           !vp(:,ip)=nint(real(nvbin-1)*atan(sqrt(pi/2)/(pm%sigv2*vrel_boost)*vreal)/pi,kind=izipv)

                           ! if ( maxval(vp(:,ip)) > 10000 .or. maxval(vmax_t) > 10000) then
                           !    print*,'update velocity out of range',ix,iy,iz,i,j,k,l,xpos,iapm,vp(:,ip)
                           !    print*,output_name('error_phi')
                           !    open(11,file=output_name('error_phi'),status='replace',access='stream')
                           !    write(11) density
                           !    close(11)

                           !    ! do k_0=k-1,k+1 ! interpolate particles to density field
                           !    !    do j_0=j-1,j+1
                           !    !       do i_0=i-1,i+1
                           !    !          np=rhoc(i_0,j_0,k_0,ix,iy,iz)
                           !    !          nzero=idx_b_r(j_0,k_0,ix,iy,iz)-sum(rhoc(i_0:,j_0,k_0,ix,iy,iz))
                           !    !          do ilayer=1,np ! loop over particle
                           !    !             ip=nzero+ilayer
                           !    !             xpos=([i_0,j_0,k_0]-1)+(int(xp(:,ip)+ishift,izipx)+rshift)*x_resolution
                           !    !             xpos=xpos*ratio_cs/pm%gridsize-pm%utile_shift*nfp(iapm)
                           !    !             print '(I16,"  xp ",3(F10.6),"  vp ",3(F10.6))', ip, xpos, vp(:,ip)
                           !    !          enddo
                           !    !       enddo
                           !    !    enddo
                           !    ! enddo

                           !    error stop
                           ! endif

                        enddo
                     enddo
                  enddo
               enddo
               !$omp endparalleldo
            enddo
         enddo
      enddo
      call system_clock(tteam2,tteam_rate)
      ttotal(itile)=real(tteam2-tteam1)/tteam_rate
   endsubroutine pm_thread

   subroutine potential_to_force(force,phi,pm)
      use omp_lib
      implicit none
      type(type_pm) pm
      integer i_0,j_0,k_0,i_n(4),j_n(4),k_n(4)
      real force(1:,0:,0:,0:)
      real phi(pm%m1phi(1):,pm%m1phi(2):,pm%m1phi(3):)
      phi=phi/(real(pm%nwork)**3) ! dividing factor in cubefft
      !$omp paralleldo num_threads(nnest) default(shared) schedule(dynamic)&
      !$omp& private(k_0,j_0,i_0,i_n,j_n,k_n)
      do k_0=0,pm%nforce+1
         do j_0=0,pm%nforce+1
            do i_0=0,pm%nforce+1
               i_n=i_0+[-2,-1,1,2]
               j_n=j_0+[-2,-1,1,2]
               k_n=k_0+[-2,-1,1,2]
               force(1,i_0,j_0,k_0)=sum(phi(i_n,j_0,k_0)*weight)
               force(2,i_0,j_0,k_0)=sum(phi(i_0,j_n,k_0)*weight)
               force(3,i_0,j_0,k_0)=sum(phi(i_0,j_0,k_n)*weight)
            enddo
         enddo
      enddo
      !$omp endparalleldo
   endsubroutine

   subroutine indexed_sort(N,ARRIN,INDX)
      implicit none
      integer(4) N ! number of halos to sort
      integer(4) IR
      integer(4) INDX(:),INDXT,J,L,I
      real ARRIN(N),Q

      DO 11 J=1,N
         INDX(J)=J
11    CONTINUE
      L=N/2+1
      IR=N
10    CONTINUE
      IF(L.GT.1)THEN
         L=L-1
         INDXT=INDX(L)
         Q=ARRIN(INDXT)
      ELSE
         INDXT=INDX(IR)
         Q=ARRIN(INDXT)
         INDX(IR)=INDX(1)
         IR=IR-1
         IF(IR.EQ.1)THEN
            INDX(1)=INDXT
            RETURN
         ENDIF
      ENDIF
      I=L
      J=L+L
20    IF(J.LE.IR)THEN
         IF(J.LT.IR)THEN
            IF(ARRIN(INDX(J)).LT.ARRIN(INDX(J+1)))J=J+1
         ENDIF
         IF(Q.LT.ARRIN(INDX(J)))THEN
            INDX(I)=INDX(J)
            I=J
            J=J+J
         ELSE
            J=IR+1
         ENDIF
         GO TO 20
      ENDIF
      INDX(I)=INDXT
      GO TO 10
   endsubroutine indexed_sort
endsubroutine kick
