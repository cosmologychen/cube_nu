
module cicpower_segment
   use omp_lib
   use cubefft
   use pencil_fft
   use parameters
   use ieee_arithmetic
   implicit none
   save

   logical use_CAMB[*]
   real xiglobal(6,0:ngbin)[*],xicoarse(6,0:ncbin)[*],xi_cdm(6,0:npbin)[*],rho8[*]
   real xitile(6,0:nnbin,nteam),xitile_avange(6,0:nnbin)[*]


contains


   subroutine get_power
      use variables
      use_CAMB = .false.
      if (sim%calculate_PK == 2) then
         call system_clock(tp1,tpr)
         call global_power
         sync all; call system_clock(tp2,tpr); if (head) print*,'     global_power elapsed time =',real(tp2-tp1)/tpr
         tps(2,cc) = real(tp2-tp1)/tpr
      elseif (sim%calculate_PK == 1) then
         call system_clock(tp1,tpr)
         call system_clock(ttp1,tpr)
         call coarse_power
         sync all; call system_clock(ttp2,tpr);
         tps(3,cc) = real(ttp2-ttp1)/tpr

         call system_clock(ttp1,tpr)
         if (.not. use_CAMB) call standard_power
         sync all; call system_clock(ttp2,tpr);
         tps(4,cc) = real(ttp2-ttp1)/tpr
         ! print*,'spower ',image,' done'
         
         if (head) then
            write(str_z,'(f7.3)') 1/sim%a-1
            print*,'   Write cicpower_segment into'
            print*,'     ',nupath//'Pk/'//trim(adjustl(str_z))//'*.bin'
            open(15,file=nupath//'Pk/'//trim(adjustl(str_z))//'_cicpower_sall.bin',status='replace',access='stream')
            write(15) xi_cdm(:,1:npbin)
            close(15)
         endif

         call system_clock(ttp1,tpr)
         if (.not. use_CAMB) call fine_power
         sync all; call system_clock(ttp2,tpr);
         tps(5,cc) = real(ttp2-ttp1)/tpr
         ! print*,'power ',image,' done'
         
         if (head) then
            write(str_z,'(f7.3)') 1/sim%a-1
            print*,'   Write cicpower_segment into'
            print*,'     ',nupath//'Pk/'//trim(adjustl(str_z))//'*.bin'
            open(15,file=nupath//'Pk/'//trim(adjustl(str_z))//'_cicpower_sall.bin',status='replace',access='stream')
            write(15) xi_cdm(:,1:npbin)
            close(15)
         endif

         sync all; call system_clock(tp2,tpr); if (head) print*,'     segment_power elapsed time =',real(tp2-tp1)/tpr
         tps(1,cc) = real(tp2-tp1)/tpr
         ! print*,'time ',image,' done'
      endif
   endsubroutine

   real function interp1_global(xdata,ydata,xq)
      implicit none
      integer(4) i_mid,i1,i2
      real xdata(ngbin),ydata(ngbin),xq
      i1=1; i2=ngbin
      do while (i2-i1>1)
         i_mid=(i1+i2)/2
         if (xq>xdata(i_mid)) then
            i1=i_mid
         else
            i2=i_mid
         endif
      enddo
      interp1_global=ydata(i1)+(xq-xdata(i1))/(xdata(i2)-xdata(i1))*(ydata(i2)-ydata(i1))
   endfunction

   subroutine pk_correction_global(xi_g,p,n_int,alpha)
      use omp_lib
      implicit none
      integer i,j,k,n_int,in,jn,kn,ibin,nplocal,icore,p
      real alpha,kvec(3),kmag,kmagn,kvecn(3),ks(3),Wk2Pk,Pk,cdata(0:ngbin,0:ncore,3),xi_g(6,0:ngbin)[*],kmaxs(26)

      kmaxs=[1.0, 1.4142135623730951, 1.7320508075688772, 2.0, 2.23606797749979, 2.449489742783178, 2.8284271247461903, 3.0, 3.1622776601683795, 3.3166247903554, 3.4641016151377544, 3.605551275463989, 3.7416573867739413, 4.0, 4.123105625617661, 4.242640687119285, 4.358898943540674, 4.47213595499958, 4.58257569495584, 4.69041575982343, 4.898979485566356, 5.0, 5.0990195135927845, 5.196152422706632, 5.385164807134504, 5.477225575051661]

      if (.not. isnan(alpha)) then
         call omp_set_num_threads(ncore)
         alpha=(log(interp1_global(xi_g(2,:),xi_g(5,:),real((nfg_global/2))))-log(interp1_global(xi_g(2,:),xi_g(5,:),real((nfg_global/2))/2)))/log(2.)

         cdata=0
         !$omp paralleldo default(shared) schedule(dynamic)&
         !$omp& private(i,icore,j,k,kvec,kmag,ibin,Wk2Pk,Pk,in,jn,kn,kvecn,kmagn,ks)
         do i=1,(nfg_global/2)+1
            icore=omp_get_thread_num()+1
            do j=1,i
               do k=1,j
                  kvec=[i,j,k]-1.0
                  kmag=norm2(kvec)
                  if (0 < kmag .and. kmag < 5.5) then
                     ibin = 1
                     do while (abs(kmag/kmaxs(ibin)-1) > 1e-4 .and. ibin <27)
                        ibin=ibin+1
                     enddo
                  elseif (kmag == 0) then
                     ibin=0
                  else
                     ibin=nint(kmag)+21
                  endif
                  Wk2Pk=0
                  Pk=kmag**alpha
                  do in=-n_int,n_int
                     do jn=-n_int,n_int
                        do kn=-n_int,n_int
                           kvecn=kvec+[in,jn,kn]*nfg_global
                           kmagn=norm2(kvecn)
                           ks=pi*kvecn/nfg_global
                           Wk2Pk=Wk2Pk+(product(merge(1.,sin(ks)/ks,ks==0))**(2*p)) * (kmagn**alpha)
                        enddo
                     enddo
                  enddo
                  cdata(ibin,icore,:)=cdata(ibin,icore,:)+[1.,kmag,Wk2Pk/Pk]
                  ! if (ibin == 27) then
                  !    print*,icore,kmag,Wk2Pk,Pk
                  ! endif
               enddo
            enddo
         enddo
         !$omp endparalleldo
         print*,'   alpha',alpha
         ! print*,cdata(26,1:ncore,3)
         ! print*,cdata(27,1:ncore,3)
         ! print*,cdata(28,1:ncore,3)
         cdata(1:ngbin,0,:)=sum(cdata(1:ngbin,1:ncore,:),dim=2)
         ! print*,cdata(26,0,3)
         ! print*,cdata(27,0,3)
         ! print*,cdata(28,0,3)
         ! print*,cdata(1:ngbin,0,1)
         ! print*,cdata(1:ngbin,0,3)
         ! stop
         cdata(1:ngbin,0,2)=cdata(1:ngbin,0,2)/cdata(1:ngbin,0,1)
         cdata(1:ngbin,0,3)=cdata(1:ngbin,0,3)/cdata(1:ngbin,0,1)
         ! print*,cdata(1:ngbin,0,3)
         xi_g(5,1:ngbin)=xi_g(4,1:ngbin)/cdata(1:ngbin,0,3)
         ! print*,xi_g(5,1:ngbin)
      endif
   endsubroutine

   subroutine cicpower_global(xi_global,n_particle)
      use omp_lib
      use pencil_fft_global
      implicit none
      integer i,j,k,ig,jg,kg,ibin
      integer(8) n_particle
      real a,kr,kx(3),C1k(3),Dk,amp11,xi_global(6,0:ngbin)[*],kmaxs(26)


      call create_penfft_plan_global
      kmaxs=[1.0, 1.4142135623730951, 1.7320508075688772, 2.0, 2.23606797749979, 2.449489742783178, 2.8284271247461903, 3.0, 3.1622776601683795, 3.3166247903554, 3.4641016151377544, 3.605551275463989, 3.7416573867739413, 4.0, 4.123105625617661, 4.242640687119285, 4.358898943540674, 4.47213595499958, 4.58257569495584, 4.69041575982343, 4.898979485566356, 5.0, 5.0990195135927845, 5.196152422706632, 5.385164807134504, 5.477225575051661]

      a=0
      xi_global=0
      call pencil_fft_forward_global
      if (head) print*, 'check: min,max of rho_k = '
      ! if (head)
      ! print*, image,minval(real(cxyz_global)),maxval(real(cxyz_global))
      ! print*,'k count',image
      do k=1,nfg/nn
         do j=1,nfg
            do i=1,(nfg_global/2)+1
               kg=(nn*(icz-1)+icy-1)*nfg/nn+k
               jg=(icx-1)*nfg+j
               ig=i
               kx=mod([ig,jg,kg]+(nfg_global/2)-1,nfg_global)-(nfg_global/2)
               if (ig==1.and.jg==1.and.kg==1) cycle ! zero frequency
               if ((ig==1.or.ig==nfg_global/2+1) .and. jg>nfg_global/2+1) cycle
               if ((ig==1.or.ig==nfg_global/2+1) .and. (jg==1.or.jg==nfg_global/2+1) .and. kg>nfg_global/2+1) cycle
               kr=sqrt(kx(1)**2+kx(2)**2+kx(3)**2)
               if (0 < kr .and. kr< 5.5) then
                  ibin = 1
                  do while (abs(kr/kmaxs(ibin)-1) > 1e-4 .and. ibin <27)
                     ibin=ibin+1
                  enddo
               elseif (kr == 0) then
                  ibin=0
               else
                  ibin=nint(kr)+21
               endif
               if (ibin > ngbin .and. ibin < 0 ) then
                  print*,image,'kr out of range'
                  print*,'kr',i,j,k,kg,jg,ig,kr
                  print*,'ngbin',ngbin
                  print*,'(nfg_global/2)',(nfg_global/2)
                  error stop
               endif
               xi_global(1,ibin)=xi_global(1,ibin)+1 ! number count
               xi_global(2,ibin)=xi_global(2,ibin)+kr ! k count
               amp11=real(cxyz_global(i,j,k)*conjg(cxyz_global(i,j,k)))
               C1k=1-(2./3.)*sin(pi*kx/nfg_global)**2! CIC
               Dk=(C1k(1)*C1k(2)*C1k(3))/n_particle
               xi_global(3,ibin)=xi_global(3,ibin)+amp11 ! raw power
               xi_global(4,ibin)=xi_global(4,ibin)+(amp11-Dk) ! P_r(k)
            enddo
         enddo
      enddo
      ! print*,'k count done',image
      call destroy_penfft_plan_global
      if (head) then ! in head node, reduce and recover P(k)
         do i=2,nn**3
            xi_global=xi_global+xi_global(:,:)[i]
         enddo
         ! print*,xi_global(3,ngbin)
         xi_global(2,:)=xi_global(2,:)/xi_global(1,:)
         xi_global(3,:)=xi_global(3,:)/xi_global(1,:) ! raw power
         ! print*,xi_global(3,ngbin),xi_global(1,ngbin)
         ! print*,xi_global(1,:)
         ! print*,xi_global(2,:)
         ! print*,xi_global(3,:)
         xi_global(4,:)=xi_global(4,:)/xi_global(1,:) ! raw power - Dk
         xi_global(5,:)=xi_global(4,:)
         call pk_correction_global(xi_global,2,3,a)
         call pk_correction_global(xi_global,2,3,a)
         call pk_correction_global(xi_global,2,3,a)
         call pk_correction_global(xi_global,2,3,a)
         ! divide and normalize
         xi_global(2,:)=xi_global(2,:)*(2*pi)/box ! k_phys
         xi_global(3,:)=xi_global(3,:)*(box**3) ! power_phys
         xi_global(4,:)=xi_global(4,:)*(box**3) ! power_phys
         xi_global(5,:)=xi_global(5,:)*(box**3) ! power_phys
         ! print*,xi_global(3,:)
      endif
      sync all
      ! if(any(ieee_is_nan(xiglobal(5,:)))) then
      !    if(head) then
      !       print*,''
      !       print*,''
      !       print*,''
      !       print*,'cicpower_global err '
      !       print*,'k_global',xiglobal(2,:)
      !       print*,'xi_global_raw',xi_global(3,:)
      !       print*,'xi_global_corr',xi_global(5,:)
      !    endif
      !    error stop 'cicpower_global err '
      ! endif
   endsubroutine cicpower_global

   subroutine global_power
      use variables
      use pencil_fft_global
      implicit none
      ! save

      integer i,j,k,l,ilayer,n,idx1(3),idx2(3),ix,iy,iz,i1,i2,i3,itile,cur_checkpoint
      integer(8) ip,nplocal,npglobal,nzero,np
      real pos1(3),dx1(3),dx2(3),xpos(ndim)
      real(4),allocatable :: rho_grid(:,:,:)[:,:,:]
      ! real(8),save :: rho8[*]
      real xi_log(6,0:ngbin),kmin,xi_real(6),xi_cdm_log(6,0:npbin)
      real a1(6),a2(6),k2(6),b2(6),c2(6)
      integer(4) i_mid
      real k_need_log,a,k_global_max,k_std_max


      allocate(rho1_global(nfg,nfg,nfg)[nn,nn,*],rho_grid(0:nfg+1,0:nfg+1,0:nfg+1)[nn,nn,*])
      sync all


      if (head) print*, ' global_power'
      ! if (head) print*,'size_rhoc',size(rho_grid)
      rho_grid=0; nzero=0
      do itz=1,nnt
         do ity=1,nnt
            do itx=1,nnt
               do k=1,nt
                  do j=1,nt
                     do i=1,nt
                        np=rhoc(i,j,k,itx,ity,itz)
                        ! print*,i,j,k,itx,ity,itz,'np',np
                        nzero=idx_b_r(j,k,itx,ity,itz)-sum(rhoc(i:,j,k,itx,ity,itz))
                        do l=1,np
                           ip=nzero+l
                           if (ip<1 )then
                              print*,'ip out of range in standard_power',image
                              print*,image,'tile',ix,iy,iz,nzero
                              print*,image,'ip',ip
                              print*,image,'subtile',i,j,k,nzero
                              print*,image,'rhoc',maxval(rhoc(:,:,:,ix,iy,iz)),minval(rhoc(:,:,:,ix,iy,iz))
                              error stop
                           endif
                           pos1=nt*((/itx,ity,itz/)-1)+ ((/i,j,k/)-1) + (int(xp(:,ip)+ishift,izipx)+rshift)*x_resolution
                           pos1=pos1*real(nfg)/real(nc) - 0.5
                           idx1=floor(pos1)+1; idx2=idx1+1
                           dx1=idx1-pos1;      dx2=1-dx1
                           if (minval(idx1)< 0 .or. maxval(idx2)>nfg+1) then
                              print*,'mass assignment out of range in standard_power'
                              print*,image,'tile',ix,iy,iz,i,j,k,ipm2
                              print*,image,'range',0,nfg+1
                              print*,image,'xp', xp(:,ip)
                              print*,image,'xp_real',  ([i,j,k]-1)+(int(xp(:,ip)+ishift,izipx)+rshift)*x_resolution
                              print*,image,'xpos',pos1
                              print*,image,'idx',idx1,idx2
                              error stop
                           endif
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
               enddo
            enddo
         enddo
      enddo
      sync all
      ! if (head) print*, 'sync buffer'
      sync all
      rho_grid(1,:,:)=rho_grid(1,:,:)+rho_grid(nfg+1,:,:)[inx,icy,icz]
      rho_grid(nfg,:,:)=rho_grid(nfg,:,:)+rho_grid(0,:,:)[ipx,icy,icz]; sync all
      rho_grid(:,1,:)=rho_grid(:,1,:)+rho_grid(:,nfg+1,:)[icx,iny,icz]
      rho_grid(:,nfg,:)=rho_grid(:,nfg,:)+rho_grid(:,0,:)[icx,ipy,icz]; sync all
      rho_grid(:,:,1)=rho_grid(:,:,1)+rho_grid(:,:,nfg+1)[icx,icy,inz]
      rho_grid(:,:,nfg)=rho_grid(:,:,nfg)+rho_grid(:,:,0)[icx,icy,ipz]; sync all
      do i=1,nfg
         do j=1,nfg
            do k=1,nfg
               rho1_global(k,j,i)=rho_grid(k,j,i)
            enddo
         enddo
      enddo
      deallocate(rho_grid)
      if (head) print*, 'check: min,max,sum of rho_grid = '
      ! if (head)
      print*, minval(rho1_global),maxval(rho1_global),sum(rho1_global*1d0),sum(rhoc*1d0)
      rho8=sum(rho1_global*1d0); sync all
      if (head) then
         do i=2,nn**3
            rho8=rho8+rho8[i]
         enddo
      endif
      ! print*,image,'rho8',rho8
      rho8=rho8/nfg/nfg/nfg; sync all
      ! print*,image,'ooo rho8',rho8
      rho8=rho8[1]
      ! print*,image,'rho8',rho8
      do i=1,nfg
         rho1_global(:,:,i)=rho1_global(:,:,i)/(rho8)-1
      enddo
      ! stop
      if (head) print*,'min',minval(rho1_global),'max',maxval(rho1_global),'mean',sum(rho1_global*1d0)/nfg/nfg/nfg; sync all
      ! if (head) print*,'Write delta_c into',output_name('delta_global')
      ! open(11,file=output_name('delta_global2'),status='replace',access='stream')
      ! write(11) rho1_global
      ! close(11); sync all
      ! open(11,file=output_name('rho_c2'),status='replace',access='stream')
      ! write(11) rhoc(1:nt,1:nt,1:nt,1:nnt,1:nnt,1:nnt)
      ! close(11); sync all
      ! stop

      if (head) print*,'auto_power'
      call cicpower_global(xiglobal,sim%npglobal)

      sync all
      if (head) then
         k_std_max = ngp/tile*pi
         k_global_max = nfg_global/box*pi

         write(str_z,'(f7.3)') 1/sim%a-1
         print*,'Write cicpower_global into',nupath//'Pk/'//trim(adjustl(str_z))//'_cicpower_global.bin'
         open(15,file=nupath//'Pk/'//trim(adjustl(str_z))//'_cicpower_global.bin',status='replace',access='stream')
         write(15) xiglobal(:,0:ngbin)[1]
         close(15)

         ! print*, '**********************global_segment *******************'
         j = 1
         do while(xiglobal(2,j) < k_global_max)
            j = j+1
         enddo
         xi_cdm = 0
         ! print*,'j',j,xiglobal(2,j),xiglobal(5,j)

         if(any(ieee_is_nan(xiglobal(5,1:j-1)))) xiglobal(5,1:ngbin) = xiglobal(3,1:ngbin)

         i = 1
         do while(xiglobal(2,i) < k_std_max)
            i = i+1
         enddo

         j = 1
         do while(xiglobal(2,j) < k_std_max*2.)
            j = j+1
         enddo

         if (xiglobal(5,j)*2 > xiglobal(5,i)) then
            k_global_max = k_std_max ! grid effect
         else
            k_global_max = k_global_max*1.731
         endif
         ! print*, k_global_max, k_std_max

         j = 1
         do while(xiglobal(2,j) < k_global_max)
            j = j+1
         enddo
         xi_cdm = 0
         ! print*,'j',j,xiglobal(2,j),xiglobal(5,j)
         ! print*,xiglobal(5,1:j-1)

         if(any(ieee_is_nan(xiglobal(5,1:j-1)))) xiglobal(5,1:ngbin) = xiglobal(3,1:ngbin)


         xi_log = log(xiglobal(:,:)[1])
         j = 1
         do while(xiglobal(1,j) <= 100.)
            j = j+1
         enddo
         ! print*,'j',j,xiglobal(1,j),xiglobal(2,j)

         xi_cdm(:5,:26) = xiglobal(:5,:26)

         i = 27
         do while(kh_lin(i) <= k_global_max)
            k_need_log = kh_lin_log(i)

            i1=1; i2=ngbin
            do while (i2-i1>1)
               i_mid=(i1+i2)/2
               if (k_need_log>xi_log(2,i_mid)) then
                  i1=i_mid
               else
                  i2=i_mid
               endif
            enddo

            if (k_need_log-xi_log(2,i1) < xi_log(2,i2)-k_need_log) then
               i3 = i2
               i2 = i1
               i1 = i1-1
            else
               i3 = i2+1
            endif

            if (i3 > ngbin) then
               print*,'out range of global bin'
               i1 = ngbin-2
               i2 = ngbin-1
               i3 = ngbin
            endif

            if (i1 < 1) then
               i1 = 1
               i2 = i2+1
               i3 = i3+2
            endif

            a1 = (xi_log(:,i1)-xi_log(:,i2))/(xi_log(2,i1)-xi_log(2,i2))
            a2 = (xi_log(:,i2)-xi_log(:,i3))/(xi_log(2,i2)-xi_log(2,i3))
            k2 = (a1-a2)/(xi_log(2,i1)-xi_log(2,i3))
            b2 = a1-k2*(xi_log(2,i1)+xi_log(2,i2))
            c2 = xi_log(:,i1)-k2*xi_log(2,i1)**2-b2*xi_log(2,i1)

            xi_cdm(:,i) = exp(k2*k_need_log**2+b2*k_need_log+c2)
            xi_cdm(2,i) = kh_lin(i)

            i = i+1
         enddo

         xi_cdm_log = log(xi_cdm)

         if (k_global_max == k_std_max) then
            ! print*,'grid'
            i1 = i-1
            do while(kh_lin(i1) > k_std_max .and. (isnan(xi_cdm_log(5,i1))))
               i1 = i1-1
            enddo
         else
            i1 = i-1
            do while(kh_lin(i1) > ng/box*pi .and. (isnan(xi_cdm_log(5,i1))))
               i1 = i1-1
            enddo
         endif

         i2 = i1
         do while(kh_lin(i2) > kh_lin(i1)/3*2)
            i2 = i2-1
         enddo

         i3 = i2
         do while(kh_lin(i3) > kh_lin(i1)/3)
            i3 = i3-1
         enddo

         print*,i,i1,i2,i3
         print*,k_global_max,kh_lin(i1),kh_lin(i1)/4*3,kh_lin(i1)/2
         print*,kh_lin(i),kh_lin(i1),kh_lin(i2),kh_lin(i3)
         print*,xi_cdm_log(5,i),xi_cdm_log(5,i1),xi_cdm_log(5,i2),xi_cdm_log(5,i3)
         print*,xi_cdm(5,i),xi_cdm(5,i1),xi_cdm(5,i2),xi_cdm(5,i3)
         i=i1
         a1 = (xi_cdm_log(:,i1)-xi_cdm_log(:,i2))/(xi_cdm_log(2,i1)-xi_cdm_log(2,i2))
         a2 = (xi_cdm_log(:,i2)-xi_cdm_log(:,i3))/(xi_cdm_log(2,i2)-xi_cdm_log(2,i3))
         k2 = (a1-a2)/(xi_cdm_log(2,i1)-xi_cdm_log(2,i3))
         b2 = a1-k2*(xi_cdm_log(2,i1)+xi_cdm_log(2,i2))
         c2 = xi_cdm_log(:,i1)-k2*xi_cdm_log(2,i1)**2-b2*xi_cdm_log(2,i1)
         do while (i <= npbin)
            k_need_log = kh_lin_log(i)
            xi_real = exp(k2*k_need_log**2+b2*k_need_log+c2)
            if(ieee_is_nan(xi_real(5))) then
               print*,'global xi_real err',i
               print*,a1(1:5)
               print*,a2(1:5)
               print*,k2(1:5)
               print*,b2(1:5)
               print*,c2(1:5)
               print*,k_need_log
               print*,xi_real(1:5)
               print*,k2(1:5)*k_need_log**2+b2(1:5)*k_need_log+c2(1:5)
               error stop 'global xi_real err'
            endif
            xi_cdm(:,i) = xi_real(:)
            xi_cdm(2,i) = kh_lin(i)
            i = i+1
         enddo

         write(str_i,'(i6)') image
         write(str_z,'(f7.3)') 1/sim%a-1
         print*,'   Write cicpower_segment into'
         print*,'     ',nupath//'Pk/'//trim(adjustl(str_z))//'*.bin'
         open(15,file=nupath//'Pk/'//trim(adjustl(str_z))//'_cicpower_gall.bin',status='replace',access='stream')
         write(15) xi_cdm(:,1:npbin)
         close(15)
      endif
      sync all
      xi_cdm(5,1:) = xi_cdm(5,1:)[1]
      if(head .and. any(ieee_is_nan(xi_cdm(5,1:i)))) then
         print*,''
         print*,''
         print*,''
         print*,'global_power err ',i,k_global_max
         print*,'xi_global',xiglobal(5,1:i)
         print*,'xi_cdm',xi_cdm(5,1:i)
         print*,'k_global',xiglobal(2,1:i)
         print*,'kh',kh_lin(:i)
         error stop 'global_power err '
      endif
      ! stop

   endsubroutine global_power

   real function interp1_coarse(xdata,ydata,xq)
      implicit none
      integer(4) i_mid,i1,i2
      real xdata(ncbin),ydata(ncbin),xq
      i1=1; i2=ncbin
      do while (i2-i1>1)
         i_mid=(i1+i2)/2
         if (xq>xdata(i_mid)) then
            i1=i_mid
         else
            i2=i_mid
         endif
      enddo
      interp1_coarse=ydata(i1)+(xq-xdata(i1))/(xdata(i2)-xdata(i1))*(ydata(i2)-ydata(i1))
   endfunction

   real function interp1_tile(xdata,ydata,xq)
      implicit none
      integer(4) i_mid,i1,i2
      real xdata(nnbin),ydata(nnbin),xq
      i1=1; i2=nnbin
      do while (i2-i1>1)
         i_mid=(i1+i2)/2
         if (xq>xdata(i_mid)) then
            i1=i_mid
         else
            i2=i_mid
         endif
      enddo
      interp1_tile=ydata(i1)+(xq-xdata(i1))/(xdata(i2)-xdata(i1))*(ydata(i2)-ydata(i1))
      ! print*,'xdata',xdata(i1),xq,xdata(i2)
      ! print*,'interp1_tile',ydata(i1),interp1_tile,ydata(i2)
   endfunction

   subroutine pk_correction_coarse(xi_c,p,n_int,alpha)
      use omp_lib
      implicit none
      integer i,j,k,n_int,in,jn,kn,ibin,nplocal,icore,p
      real alpha,kvec(3),kmag,kmagn,kvecn(3),ks(3),Wk2Pk,Pk,cdata(0:ncbin,0:ncore,3),xi_c(6,0:ncbin)[*],kmaxs(26)

      kmaxs=[1.0, 1.4142135623730951, 1.7320508075688772, 2.0, 2.23606797749979, 2.449489742783178, 2.8284271247461903, 3.0, 3.1622776601683795, 3.3166247903554, 3.4641016151377544, 3.605551275463989, 3.7416573867739413, 4.0, 4.123105625617661, 4.242640687119285, 4.358898943540674, 4.47213595499958, 4.58257569495584, 4.69041575982343, 4.898979485566356, 5.0, 5.0990195135927845, 5.196152422706632, 5.385164807134504, 5.477225575051661]

      if (.not. isnan(alpha)) then
         call omp_set_num_threads(ncore)
         alpha=(log(interp1_coarse(xi_c(2,:),xi_c(5,:),real(nyquist)))-log(interp1_coarse(xi_c(2,:),xi_c(5,:),real(nyquist)/2)))/log(2.)

         cdata=0
         !$omp paralleldo default(shared) schedule(dynamic)&
         !$omp& private(i,icore,j,k,kvec,kmag,ibin,Wk2Pk,Pk,in,jn,kn,kvecn,kmagn,ks)
         do i=1,nyquist+1
            icore=omp_get_thread_num()+1
            do j=1,i
               do k=1,j
                  kvec=[i,j,k]-1.0
                  kmag=norm2(kvec)
                  if (0 < kmag .and. kmag < 5.5) then
                     ibin = 1
                     do while (abs(kmag/kmaxs(ibin)-1) > 1e-4 .and. ibin <27)
                        ibin=ibin+1
                     enddo
                  elseif (kmag == 0) then
                     ibin=0
                  else
                     ibin=nint(kmag)+21
                  endif
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
                  ! if (ibin == 27) then
                  !    print*,icore,kmag,Wk2Pk,Pk
                  ! endif
               enddo
            enddo
         enddo
         !$omp endparalleldo
         print*,'   alpha',alpha
         ! print*,cdata(26,1:ncore,3)
         ! print*,cdata(27,1:ncore,3)
         ! print*,cdata(28,1:ncore,3)
         cdata(1:ncbin,0,:)=sum(cdata(1:ncbin,1:ncore,:),dim=2)
         ! print*,cdata(26,0,3)
         ! print*,cdata(27,0,3)
         ! print*,cdata(28,0,3)
         ! print*,cdata(1:ncbin,0,1)
         ! print*,cdata(1:ncbin,0,3)
         ! stop
         cdata(1:ncbin,0,2)=cdata(1:ncbin,0,2)/cdata(1:ncbin,0,1)
         cdata(1:ncbin,0,3)=cdata(1:ncbin,0,3)/cdata(1:ncbin,0,1)
         ! print*,cdata(1:ncbin,0,3)
         xi_c(5,1:ncbin)=xi_c(4,1:ncbin)/cdata(1:ncbin,0,3)
         ! print*,xi_c(5,1:ncbin)
      endif
   endsubroutine

   subroutine pk_correction_tile(xi_t,p,n_int,nyquist_tile,alpha)
      use omp_lib
      implicit none
      integer i,j,k,n_int,in,jn,kn,ibin,nplocal,icore,p
      real alpha,kvec(3),kmag,kmagn,kvecn(3),ks(3),Wk2Pk,Pk,cdata(0:nnbin,0:ncore,3),xi_t(6,0:nnbin)
      integer(8) nyquist_tile

      call omp_set_num_threads(ncore)
      alpha=(log(interp1_tile(xi_t(2,:),xi_t(5,:),real(nyquist_tile)))-log(interp1_tile(xi_t(2,:),xi_t(5,:),real(nyquist_tile)/2)))/log(2.)

      if (.not. isnan(alpha)) then
         cdata=0
         !$omp paralleldo default(shared) num_threads(ncore) schedule(dynamic,1)&
         !$omp& private(i,icore,j,k,kvec,kmag,ibin,Wk2Pk,Pk,in,jn,kn,kvecn,kmagn,ks)
         do i=1,nyquist_tile+1
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
                           kvecn=kvec+[in,jn,kn]*nyquist_tile*2
                           kmagn=norm2(kvecn)
                           ks=pi*kvecn/nyquist_tile/2
                           Wk2Pk=Wk2Pk+(product(merge(1.,sin(ks)/ks,ks==0))**(2*p)) * (kmagn**alpha)
                        enddo
                     enddo
                  enddo
                  ! print*,'k',kvec,kmag,Wk2Pk,Pk
                  cdata(ibin,icore,:)=cdata(ibin,icore,:)+[1.,kmag,Wk2Pk/Pk]
               enddo
            enddo
         enddo
         !$omp endparalleldo
         print*,'    alpha',alpha
         cdata(1:nnbin,0,:)=sum(cdata(1:nnbin,1:ncore,:),dim=2)
         cdata(1:nnbin,0,2)=cdata(1:nnbin,0,2)/cdata(1:nnbin,0,1)
         cdata(1:nnbin,0,3)=cdata(1:nnbin,0,3)/cdata(1:nnbin,0,1)
         xi_t(5,1:nnbin)=xi_t(4,1:nnbin)/cdata(1:nnbin,0,3)
      end if
   endsubroutine

   subroutine cicpower_coarse(xi_coarse,n_particle)
      use omp_lib
      use pencil_fft
      implicit none
      integer i,j,k,ig,jg,kg,ibin
      integer(8) n_particle
      real a,kr,kx(3),C1k(3),Dk,amp11,xi_coarse(6,0:ncbin)[*],kmaxs(26)

      kmaxs=[1.0, 1.4142135623730951, 1.7320508075688772, 2.0, 2.23606797749979, 2.449489742783178, 2.8284271247461903, 3.0, 3.1622776601683795, 3.3166247903554, 3.4641016151377544, 3.605551275463989, 3.7416573867739413, 4.0, 4.123105625617661, 4.242640687119285, 4.358898943540674, 4.47213595499958, 4.58257569495584, 4.69041575982343, 4.898979485566356, 5.0, 5.0990195135927845, 5.196152422706632, 5.385164807134504, 5.477225575051661]
      a=0
      xi_coarse=0
      call pencil_fft_forward
      cxyz=cxyz/nw_global/nw_global/nw_global
      sync all
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
               if (0 < kr .and. kr< 5.5) then
                  ibin = 1
                  do while (abs(kr/kmaxs(ibin)-1) > 1e-4 .and. ibin <27)
                     ibin=ibin+1
                  enddo
               elseif (kr == 0) then
                  ibin=0
               else
                  ibin=nint(kr)+21
               endif
               if (ibin > ncbin) then
                  print*,image,'kr out of range'
                  print*,'kr',i,j,k,kg,jg,ig,kr
                  print*,'ncbin',ncbin
                  print*,'nyquist',nyquist
                  error stop
               endif
               xi_coarse(1,ibin)=xi_coarse(1,ibin)+1 ! number count
               xi_coarse(2,ibin)=xi_coarse(2,ibin)+kr ! k count
               amp11=real(cxyz(i,j,k)*conjg(cxyz(i,j,k)))
               C1k=1-(2./3.)*sin(pi*kx/nw_global)**2! CIC
               Dk=(C1k(1)*C1k(2)*C1k(3))/n_particle
               xi_coarse(3,ibin)=xi_coarse(3,ibin)+amp11 ! raw power
               xi_coarse(4,ibin)=xi_coarse(4,ibin)+(amp11-Dk) ! P_r(k)
            enddo
         enddo
      enddo
      ! print*,xi_coarse(3,ncbin)
      sync all
      if (head) then ! in head node, reduce and recover P(k)
         do i=2,nn**3
            xi_coarse=xi_coarse+xi_coarse(:,:)[i]
         enddo
         ! print*,xi_coarse(3,ncbin)
         xi_coarse(2,:)=xi_coarse(2,:)/xi_coarse(1,:)
         xi_coarse(3,:)=xi_coarse(3,:)/xi_coarse(1,:) ! raw power
         ! print*,xi_coarse(3,ncbin),xi_coarse(1,ncbin)
         ! print*,xi_coarse(1,:)
         ! print*,xi_coarse(2,:)
         ! print*,xi_coarse(3,:)
         xi_coarse(4,:)=xi_coarse(4,:)/xi_coarse(1,:) ! raw power - Dk
         xi_coarse(5,:)=xi_coarse(4,:)
         call pk_correction_coarse(xi_coarse,2,3,a)
         call pk_correction_coarse(xi_coarse,2,3,a)
         call pk_correction_coarse(xi_coarse,2,3,a)
         call pk_correction_coarse(xi_coarse,2,3,a)
         ! divide and normalize
         xi_coarse(2,:)=xi_coarse(2,:)*(2*pi)/box ! k_phys
         xi_coarse(3,:)=xi_coarse(3,:)*(box**3) ! power_phys
         xi_coarse(4,:)=xi_coarse(4,:)*(box**3) ! power_phys
         xi_coarse(5,:)=xi_coarse(5,:)*(box**3) ! power_phys
      endif
      sync all
   endsubroutine cicpower_coarse

   subroutine cicpower_tile(xi_tile,tile1_k,nn_t,n_particle)
      use omp_lib
      use cubefft

      integer(8) nyquist_tile,n_particle,nn_t
      integer i,j,k,ibin,bin
      ! integer(8),parameter :: nbin_tile=nint(nyquist_tile*sqrt(3.))
      complex tile1_k(nn_t/2+1,nn_t,nn_t)
      real kr,kx(3),C1k(3),Dk,amp11,xi_tile(6,0:nnbin)

      ! if (head) print*, '********************** cic_tile *******************'
      ! print*,n_particle
      tile1_k = tile1_k/nn_t/nn_t/nn_t
      xi_tile=0
      nyquist_tile = nn_t/2
      do k=1,nn_t
         do j=1,nn_t
            do i=1,nyquist_tile+1
               kx=mod([i,j,k]+nyquist_tile-1,nn_t)-nyquist_tile
               if (i==1.and.j==1.and.k==1) cycle ! zero frequency
               kr=sqrt(kx(1)**2+kx(2)**2+kx(3)**2)
               ibin=nint(kr)
               xi_tile(1,ibin)=xi_tile(1,ibin)+1 ! number count
               xi_tile(2,ibin)=xi_tile(2,ibin)+kr ! k count
               amp11=real(tile1_k(i,j,k)*conjg(tile1_k(i,j,k)))
               C1k=1-(2./3.)*sin(pi*kx/nn_t)**2! CIC
               Dk=(C1k(1)*C1k(2)*C1k(3))/n_particle
               xi_tile(3,ibin)=xi_tile(3,ibin)+amp11 ! raw power
               xi_tile(4,ibin)=xi_tile(4,ibin)+(amp11-Dk) ! P_r(k)
            enddo
         enddo
      enddo
      ! if (head) print*, '********************** cic_init_done *******************'
   endsubroutine cicpower_tile

   subroutine coarse_power
      use variables
      use pencil_fft
      implicit none
      ! save

      integer i,j,k,l,ilayer,n,idx1(3),idx2(3),ix,iy,iz,i1,i2,i3,itile,cur_checkpoint
      integer(8) ip,nplocal,npglobal,nzero,np
      real pos1(3),dx1(3),dx2(3),xpos(ndim)
      real(4),allocatable :: rho_grid(:,:,:)[:,:,:]
      ! real(8),save :: rho8[*]
      real xi_log(6,0:ncbin),kmin,xi_real(6)
      real a1(6),a2(6),k2(6),b2(6),c2(6)
      integer(4) i_mid
      real k_need_log,a,k_coarse_max


      allocate(rho_grid(0:nw+1,0:nw+1,0:nw+1)[nn,nn,*])


      if (head) print*, ' coarse_power'
      ! if (head) print*,'size_rhoc',size(rho_grid)
      rho_grid=0; nzero=0
      do itz=1,nnt
         do ity=1,nnt
            do itx=1,nnt
               do k=1,nt
                  do j=1,nt
                     do i=1,nt
                        np=rhoc(i,j,k,itx,ity,itz)
                        ! print*,i,j,k,itx,ity,itz,'np',np
                        nzero=idx_b_r(j,k,itx,ity,itz)-sum(rhoc(i:,j,k,itx,ity,itz))
                        do l=1,np
                           ip=nzero+l
                           if (ip<1 )then
                              print*,'ip out of range in standard_power',image
                              print*,image,'tile',ix,iy,iz,nzero
                              print*,image,'ip',ip
                              print*,image,'subtile',i,j,k,nzero
                              print*,image,'rhoc',maxval(rhoc(:,:,:,ix,iy,iz)),minval(rhoc(:,:,:,ix,iy,iz))
                              error stop
                           endif
                           pos1=nt*((/itx,ity,itz/)-1)+ ((/i,j,k/)-1) + (int(xp(:,ip)+ishift,izipx)+rshift)*x_resolution
                           pos1=pos1*real(nw)/real(nc) - 0.5
                           idx1=floor(pos1)+1; idx2=idx1+1
                           dx1=idx1-pos1;      dx2=1-dx1
                           if (minval(idx1)< 0 .or. maxval(idx2)>nw+1) then
                              print*,'mass assignment out of range in standard_power'
                              print*,image,'tile',ix,iy,iz,i,j,k,ipm2
                              print*,image,'range',1-ngb,ngp+ngb,1-ncb,nt+ncb
                              print*,image,'xp', xp(:,ip)
                              print*,image,'xp_real',  ([i,j,k]-1)+(int(xp(:,ip)+ishift,izipx)+rshift)*x_resolution
                              print*,image,'xpos',pos1
                              print*,image,'idx',idx1,idx2
                              error stop
                           endif
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
               enddo
            enddo
         enddo
      enddo
      sync all
      ! if (head) print*, 'sync buffer'
      sync all
      rho_grid(1,:,:)=rho_grid(1,:,:)+rho_grid(nw+1,:,:)[inx,icy,icz]
      rho_grid(nw,:,:)=rho_grid(nw,:,:)+rho_grid(0,:,:)[ipx,icy,icz]; sync all
      rho_grid(:,1,:)=rho_grid(:,1,:)+rho_grid(:,nw+1,:)[icx,iny,icz]
      rho_grid(:,nw,:)=rho_grid(:,nw,:)+rho_grid(:,0,:)[icx,ipy,icz]; sync all
      rho_grid(:,:,1)=rho_grid(:,:,1)+rho_grid(:,:,nw+1)[icx,icy,inz]
      rho_grid(:,:,nw)=rho_grid(:,:,nw)+rho_grid(:,:,0)[icx,icy,ipz]; sync all
      do i=1,nw
         do j=1,nw
            do k=1,nw
               rho1(k,j,i)=rho_grid(k,j,i)
            enddo
         enddo
      enddo
      ! if (head) print*, 'check: min,max,sum of rho_grid = '
      ! if (head) print*, minval(rho1),maxval(rho1),sum(rho1*1d0)
      rho8=sum(rho1*1d0); sync all
      if (head) then
         do i=2,nn**3
            rho8=rho8+rho8[i]
         enddo
      endif; sync all
      rho8=rho8/nw_global/nw_global/nw_global; sync all
      ! print*,image,'ooo rho8',rho8
      rho8=rho8[1]
      ! print*,image,'rho8',rho8
      do i=1,nw
         rho1(:,:,i)=rho1(:,:,i)/(rho8)-1
      enddo
      ! if (head) print*,'min',minval(rho1),'max',maxval(rho1),'mean',sum(rho1*1d0)/nw/nw/nw; sync all
      ! if (head) print*,'Write delta_c into',output_name('delta_coarse')
      ! open(11,file=output_name('delta_coarse'),status='replace',access='stream')
      ! write(11) rho1
      ! close(11); sync all
      deallocate(rho_grid)

      ! if (head) print*,'auto_power'
      call cicpower_coarse(xicoarse,sim%npglobal)

      sync all
      if (head) then
         k_coarse_max = nw_global/box*pi
         write(str_z,'(f7.3)') 1/sim%a-1
         ! print*,'Write cicpower_coarse into',opath//'image'//trim(adjustl(str_i))//'/'//trim(adjustl(str_z))//'cicpower_coarse.bin'
         open(15,file=nupath//'Pk/'//trim(adjustl(str_z))//'_cicpower_coarse.bin',status='replace',access='stream')
         write(15) xicoarse(:,0:ncbin)[1]
         close(15)


         ! print*, '**********************coarse_segment *******************'
         j = 1
         do while(xicoarse(2,j) < k_coarse_max)
            j = j+1
         enddo
         xi_cdm = 0
         if(any(ieee_is_nan(xicoarse(5,1:j))) .or. any(xicoarse(5,1:j)<0)) then

            use_CAMB = .true.

            xicoarse(5,1:ncbin) = xicoarse(3,1:ncbin)
            ! print*,xicoarse(5,ncbin)

            write(str_z,'(f8.4)') z_powerpoint(sim%cur_powerpoint)
            if(head) print*,'nan in xi_cdm_coarse use CAMB',str_z,z_powerpoint(sim%cur_powerpoint),sim%cur_powerpoint
            open(11,file=nupath//'Pk_cb_'//trim(adjustl(str_z))//'.txt',form='formatted')
            read(11,*) xi_cdm(6,1:)
            close(11)

            xi_log = log(xicoarse(:,:)[1])

            j = 1
            do while(kh_lin(j) < k_coarse_max/8)
               j = j+1
            enddo

            xi_cdm(:5,:26) = xicoarse(:5,:26)

            i = 27
            do while(kh_lin(i) <= k_coarse_max)
               k_need_log = kh_lin_log(i)

               i1=1; i2=ncbin
               do while (i2-i1>1)
                  i_mid=(i1+i2)/2
                  if (k_need_log>xi_log(2,i_mid)) then
                     i1=i_mid
                  else
                     i2=i_mid
                  endif
               enddo

               if (k_need_log-xi_log(2,i1) < xi_log(2,i2)-k_need_log) then
                  i3 = i2
                  i2 = i1
                  i1 = i1-1
               else
                  i3 = i2+1
               endif

               if (i3 > ncbin) then
                  print*,'out range of coarse bin'
                  i1 = ncbin-2
                  i2 = ncbin-1
                  i3 = ncbin
               endif

               if (i1 < 1) then
                  i1 = 1
                  i2 = i2+1
                  i3 = i3+2
               endif
               a1 = (xi_log(:,i1)-xi_log(:,i2))/(xi_log(2,i1)-xi_log(2,i2))
               a2 = (xi_log(:,i2)-xi_log(:,i3))/(xi_log(2,i2)-xi_log(2,i3))
               k2 = (a1-a2)/(xi_log(2,i1)-xi_log(2,i3))
               b2 = a1-k2*(xi_log(2,i1)+xi_log(2,i2))
               c2 = xi_log(:,i1)-k2*xi_log(2,i1)**2-b2*xi_log(2,i1)
               xi_real = exp(k2*k_need_log**2+b2*k_need_log+c2)
               xi_cdm(5,i) = (1.15**(j-i)*xi_real(5)+1.15**(i-j)*xi_cdm(6,i))/(1.15**(j-i)+1.15**(i-j))

               if(ieee_is_nan(xi_cdm(5,i))) then
                  if(head) then
                     print*,''
                     print*,''
                     print*,i,i1,i2,i3,ncbin
                     print*,'interp',a1(5),a2(5),k2(5),b2(5),c2(5),xi_log(5,i1),xi_log(5,i2),xi_log(5,i3)
                     print*,'xi_coarse',xi_real(5),xicoarse(5,i1),xicoarse(5,i2),xicoarse(5,i3)
                     print*,'xi_cdm',xi_cdm(5,i)
                     print*,'k_coarse',xicoarse(2,i1),xicoarse(2,i2),xicoarse(2,i3)
                     print*,'kh',kh_lin(i)
                  endif
                  error stop 'coarse_power err '
               endif
               xi_cdm(4,i) = (1.15**(j-i)*xi_real(4)+1.15**(i-j)*xi_cdm(6,i))/(1.15**(j-i)+1.15**(i-j))
               xi_cdm(3,i) = (1.15**(j-i)*xi_real(3)+1.15**(i-j)*xi_cdm(6,i))/(1.15**(j-i)+1.15**(i-j))
               xi_cdm(2,i) = kh_lin(i)
               xi_cdm(1,i) = xi_cdm(1,i)+xi_real(1)

               i = i+1
            enddo

            xi_cdm(5,i-1:) = xi_cdm(6,i-1:)
            xi_cdm(4,i-1:) = xi_cdm(6,i-1:)
            xi_cdm(3,i-1:) = xi_cdm(6,i-1:)
            xi_cdm(2,i-1:) = kh_lin(i-1:)
            xi_cdm(1,i-1:) = 0
         else
            xi_log = log(xicoarse(:,:)[1])
            j = 1
            do while(xicoarse(1,j) <= 100.)
               j = j+1
            enddo
            ! print*,'j',j,xicoarse(1,j),xicoarse(2,j)
            i = 1
            do while(kh_lin(i) <= maxval(xicoarse(2,:))/1.731)
               k_need_log = kh_lin_log(i)

               i1=1; i2=ncbin
               do while (i2-i1>1)
                  i_mid=(i1+i2)/2
                  if (k_need_log>xi_log(2,i_mid)) then
                     i1=i_mid
                  else
                     i2=i_mid
                  endif
               enddo

               if (k_need_log-xi_log(2,i1) < xi_log(2,i2)-k_need_log) then
                  i3 = i2
                  i2 = i1
                  i1 = i1-1
               else
                  i3 = i2+1
               endif

               if (i3 > ncbin) then
                  print*,'out range of coarse bin'
                  i1 = ncbin-2
                  i2 = ncbin-1
                  i3 = ncbin
               endif

               if (i1 < 1) then
                  i1 = 1
                  i2 = i2+1
                  i3 = i3+2
               endif

               ! print*,"_______________",i1,i2,i3
               ! print*,xicoarse(2,i1),xicoarse(2,i2),xicoarse(2,i3),kh_lin(i)
               ! print*,xicoarse(5,i1),xicoarse(5,i2),xicoarse(5,i3)
               ! print*,"_______________"

               a1 = (xi_log(:,i1)-xi_log(:,i2))/(xi_log(2,i1)-xi_log(2,i2))
               a2 = (xi_log(:,i2)-xi_log(:,i3))/(xi_log(2,i2)-xi_log(2,i3))
               k2 = (a1-a2)/(xi_log(2,i1)-xi_log(2,i3))
               b2 = a1-k2*(xi_log(2,i1)+xi_log(2,i2))
               c2 = xi_log(:,i1)-k2*xi_log(2,i1)**2-b2*xi_log(2,i1)

               xi_cdm(:,i) = exp(k2*k_need_log**2+b2*k_need_log+c2)
               xi_cdm(2,i) = kh_lin(i)

               i = i+1
            enddo
            xi_cdm(:,i:) = 0
            ! do j =1,i
            ! print*,xi_cdm(2,j),xi_cdm(5,j)
            ! enddo
         endif
      endif
      sync all
      xi_cdm(5,1:) = xi_cdm(5,1:)[1]
      use_CAMB = use_CAMB[1]
      if(head .and. any(ieee_is_nan(xi_cdm(5,1:i)))) then
         print*,''
         print*,''
         print*,''
         print*,'coarse_power err ',i,k_coarse_max
         print*,'xi_coarse',xicoarse(5,1:i)
         print*,'xi_cdm',xi_cdm(5,1:i)
         print*,'k_coarse',xicoarse(2,1:i)
         print*,'kh',kh_lin(:i)
         error stop 'coarse_power err '
      endif

   endsubroutine coarse_power

   subroutine standard_power
      use variables
      implicit none
      ! save

      integer i,j,k,l,ilayer,n,nntc,idx1(3),idx2(3),ix,iy,iz,i1,i2,i3,itile,cur_checkpoint,tcpu
      integer(8) nlast,ip,nplocal,npglobal,nzero,np
      real pos1(3),dx1(3),dx2(3),xpos(ndim)
      real rho_th(-ngb:ngp+ngb+1,-ngb:ngp+ngb+1,-ngb:ngp+ngb+1)
      ! real,allocatable :: rho_th(:,:,:),rho_s(:,:,:)
      integer,parameter :: nlayer=3 ! thread save for CIC intepolation
      real xi_log(6,0:nnbin),xi_real(6),xi_cdm_log(6,0:npbin)
      real a1(6),a2(6),k2(6),b2(6),c2(6)
      integer(4) i_mid
      real k_need_log,a,k_std_max


      if (head) print*, ' standard_power'
      xitile_avange = 0
      nntc = 0
      sync all

      ! print*,'********************** standard_avange0 *******************'
      ! print*,'xi_standard',xitile_avange(3,:)
      ! print*,'k_standard ',xitile_avange(2,:)
      !$omp paralleldo default(shared) num_threads(nteam) schedule(dynamic,1)&
      !$omp& private(n,ipm2,ix,iy,iz,ilayer,nlast,iteam,rho_th,str_i)
      ! do n=1,2
      do n=1,nnt**3
         ! ipm2 = ntpk(n)
         ipm2 = n
         ix = ixyz2(1,ipm2)
         iy = ixyz2(2,ipm2)
         iz = ixyz2(3,ipm2)
         nlast=sum(rhoc(1:nt,1:nt,1:nt,ix,iy,iz))
         ! print*,image,'tile',ix,iy,iz,ipm2,nlast

         ! do ilayer=0,nlayer-1
         !   do k=1-ncb+ilayer,nt+ncb,nlayer
         !   do j=1-ncb,nt+ncb
         !   do i=1-ncb,nt+ncb
         !     np=rhoc(i,j,k,ix,iy,iz)
         !     nlast = nlast+np
         !   enddo
         !   enddo
         !   enddo
         ! enddo
         ! print*,ipm2,nlast,sum(rhoc(1-ncb:nt+ncb,1-ncb:nt+ncb,1-ncb:nt+ncb,ix,iy,iz))
         ! stop
         ! print*,image,'tile',ix,iy,iz,ipm2,nlast


         if (nlast<0 )then
            print*,image,'tile',ix,iy,iz,ipm2,nlast
            print*,image,'rho_c',maxval(rhoc(:,:,:,ix,iy,iz)),minval(rhoc(:,:,:,ix,iy,iz))
            error stop
         endif
         iteam=omp_get_thread_num()+1
         xitile(:,:,iteam) = 1
         ! if (head) print*, n,ipm2,iteam
         rho_th = 0
         do ilayer=0,nlayer-1
            !$omp paralleldo default(shared) num_threads(nnest) schedule(dynamic,1)&
            !$omp& private(k,j,i,l,np,nzero,ip,pos1,dx1,dx2,idx1,idx2)
            do k=0+ilayer,nt+1,nlayer
               ! print*,image,'k',k,nlast
               do j=0,nt+1
                  do i=0,nt+1
                     np=rhoc(i,j,k,ix,iy,iz)
                     ! nlast = nlast+np
                     nzero=idx_b_r(j,k,ix,iy,iz)-sum(rhoc(i:,j,k,ix,iy,iz))
                     do l=1,np
                        ip=nzero+l
                        if (ip<1 )then
                           print*,'ip out of range in standard_power',image
                           print*,image,'tile',ix,iy,iz,nlast
                           print*,image,'ip',ip
                           print*,image,'subtile',i,j,k,nzero
                           print*,image,'rhoc',maxval(rhoc(:,:,:,ix,iy,iz)),minval(rhoc(:,:,:,ix,iy,iz))
                           error stop
                        endif
                        pos1=([i,j,k]-1)+(int(xp(:,ip)+ishift,izipx)+rshift)*x_resolution
                        pos1=pos1*ratio_cs-0.5
                        idx1=floor(pos1)+1; idx2=idx1+1
                        dx1=idx1-pos1;      dx2=1-dx1
                        if (minval(idx1)< -ngb .or. maxval(idx2)>ngp+ngb+1) then
                           print*,'mass assignment out of range in standard_power'
                           print*,image,'tile',ix,iy,iz,i,j,k,ipm2
                           print*,image,'range',1-ngb,ngp+ngb
                           print*,image,'xp', xp(:,ip)
                           print*,image,'xp_real',  ([i,j,k]-1)+(int(xp(:,ip)+ishift,izipx)+rshift)*x_resolution
                           print*,image,'xpos',pos1
                           print*,image,'idx',idx1,idx2
                           error stop
                        endif
                        rho_th(idx1(1),idx1(2),idx1(3))=rho_th(idx1(1),idx1(2),idx1(3))+dx1(1)*dx1(2)*dx1(3)
                        rho_th(idx2(1),idx1(2),idx1(3))=rho_th(idx2(1),idx1(2),idx1(3))+dx2(1)*dx1(2)*dx1(3)
                        rho_th(idx1(1),idx2(2),idx1(3))=rho_th(idx1(1),idx2(2),idx1(3))+dx1(1)*dx2(2)*dx1(3)
                        rho_th(idx1(1),idx1(2),idx2(3))=rho_th(idx1(1),idx1(2),idx2(3))+dx1(1)*dx1(2)*dx2(3)
                        rho_th(idx1(1),idx2(2),idx2(3))=rho_th(idx1(1),idx2(2),idx2(3))+dx1(1)*dx2(2)*dx2(3)
                        rho_th(idx2(1),idx1(2),idx2(3))=rho_th(idx2(1),idx1(2),idx2(3))+dx2(1)*dx1(2)*dx2(3)
                        rho_th(idx2(1),idx2(2),idx1(3))=rho_th(idx2(1),idx2(2),idx1(3))+dx2(1)*dx2(2)*dx1(3)
                        rho_th(idx2(1),idx2(2),idx2(3))=rho_th(idx2(1),idx2(2),idx2(3))+dx2(1)*dx2(2)*dx2(3)
                     enddo
                  enddo
               enddo
            enddo
            !$omp endparalleldo
         enddo

         ! rho8_s=sum(rho_th*1d0)/ngt/ngt/ngt
         ! print*,image,'rho8',rho8/ratio_cs/ratio_cs/ratio_cs,sum(rho_th*1d0)/ngt/ngt/ngt
         ! print*,'np, sum_rho :',nlast,sum(rho_th(1:ngp,1:ngp,1:ngp)),nlast/sum(rho_th(1:ngp,1:ngp,1:ngp))-1
         ! print*,sum(rhoc(1:nt,1:nt,1:nt,ix,iy,iz)),sum(rho_th(0:ngp+1,0:ngp+1,0:ngp+1)),sum(rho_th(1:ngp,1:ngp,1:ngp)),sum(rhoc(1-ncb:nt+ncb,1-ncb:nt+ncb,1-ncb:nt+ncb,ix,iy,iz)),sum(rho_th)
         ! do i=-ngb,ngp+ngb+1
         !    rho_th(:,:,i)=rho_th(:,:,i)/(rho8/ratio_cs/ratio_cs/ratio_cs)-1
         ! enddo

         ! if (head) print*, '********************** standard_fft *******************' ,iteam,ipm2
         rho2p(:,:,:,iteam) = 0         
         rho2p(:ngp,:ngp,:ngp,iteam)=rho_th(1:ngp,1:ngp,1:ngp)!/(rho8/ratio_cs/ratio_cs/ratio_cs)-1
         ! print*,sum(rho2p(:ngp,:ngp,:ngp,iteam)),nlast


         ! write(str_i,'(I)') ipm2
         ! print*,str_i,ix,iy,iz,1-ngb,ngp+ngb,ngt
         ! print*,'Write delta_std into',output_name('delta_std'//trim(adjustl(str_i)))
         ! open(15+iteam,file=output_name('delta_std'//trim(adjustl(str_i))),status='replace',access='stream')
         ! write(15+iteam) rho2(1-ngb:ngp+ngb,1-ngb:ngp+ngb,1-ngb:ngp+ngb,iteam)
         ! close(15+iteam); sync all
         ! ! stop


         ! if (head) print*, '********************** standard_cicpower *******************',iteam,maxval(rho2(:,:,:,iteam))
         ! print*,'max rho_k',real(maxval(int(rho2k(:,:,:,iteam)*conjg(rho2k(:,:,:,iteam)))))/ngt/ngt/ngt
         call sfftw_execute(plan2p(iteam))
         call cicpower_tile(xitile(:,:,iteam),rho2kp(:,:,:,iteam),ngp,nlast)
         !  print*,ipm2,'k_standard ',xitile(2,:,iteam)
         ! write(str_i,'(i6)') ipm2
         ! write(str_z,'(f7.3)') z_checkpoint(sim%cur_checkpoint)
         ! print*,'Write cicpower_std into',opath//'image1/'//trim(adjustl(str_z))//'_cicpower_stdc_'//trim(adjustl(str_i))//'.bin'
         ! open(15+iteam,file=opath//'image1/'//trim(adjustl(str_z))//'_cicpower_stdc_'//trim(adjustl(str_i))//'.bin',status='replace',access='stream')
         ! write(15+iteam) xitile(:,1:nnbin,iteam)
         ! close(15+iteam)
         ! ! stop
         ! print*,image,ipm2,'xi',xitile(3,10,iteam)
         ! print*,n,'n ',xitile(3,:20,iteam)

         !$omp barrier
         if (iteam == 1) then
            do i=1,nteam
               xitile_avange(:,0:) = xitile_avange(:,0:)+xitile(:,0:,i)
               nntc = nntc+1
            enddo
         endif
         !$omp barrier
      enddo
      !$omp endparalleldo
      ! stop

      ! print*, image,'********************** standard_avange *******************'
      ! print*,'xi_standard',xitile_avange(3,:)
      ! print*,'k_standard ',xitile_avange(2,:)
      ! print*,nntc,xitile_avange(1,:60)
      xitile_avange = xitile_avange/nntc
      sync all
      if (head) then
         do i = 2,nn**3
            xitile_avange(:,0:) = xitile_avange(:,0:) + xitile_avange(:,0:)[i]
         enddo
         xitile_avange = xitile_avange/(nn**3)


         a = 0
         xitile_avange(2,:)=xitile_avange(2,:)/xitile_avange(1,:)
         xitile_avange(3,:)=xitile_avange(3,:)/xitile_avange(1,:) ! raw power
         xitile_avange(4,:)=xitile_avange(4,:)/xitile_avange(1,:) ! raw power - Dk
         xitile_avange(5,:)=xitile_avange(4,:)
         call pk_correction_tile(xitile_avange,2,3,ngp/2,a)
         call pk_correction_tile(xitile_avange,2,3,ngp/2,a)
         call pk_correction_tile(xitile_avange,2,3,ngp/2,a)
         call pk_correction_tile(xitile_avange,2,3,ngp/2,a)
         ! divide and normalize
         xitile_avange(2,:)=xitile_avange(2,:)*(2*pi)/tile ! k_phys
         xitile_avange(3,:)=xitile_avange(3,:)*(tile**3) ! power_phys
         xitile_avange(4,:)=xitile_avange(4,:)*(tile**3) ! power_phys
         xitile_avange(5,:)=xitile_avange(5,:)*(tile**3) ! power_phys

         k_std_max = ngp/tile*pi
         write(str_i,'(i6)') image
         write(str_z,'(f7.3)') 1/sim%a-1
         ! print*,'Write cicpower_std into',opath//'image'//trim(adjustl(str_i))//'/'//trim(adjustl(str_z))//'_cicpower_std.bin'
         open(15,file=nupath//'Pk/'//trim(adjustl(str_z))//'_cicpower_std.bin',status='replace',access='stream')
         write(15) xitile_avange(:,1:)
         close(15)
         ! stop

          print*, '********************** standard_segment *******************',k_std_max

         j = 1
         do while(xitile_avange(2,j) < k_std_max)
            j = j+1
         enddo

         xi_log = log(xitile_avange)
         if(any(ieee_is_nan(xi_log(5,1:j-1)))) then

            use_CAMB = .true.

            xi_log(5,1:) = xi_log(3,1:)

            write(str_z,'(f8.4)') z_powerpoint(sim%cur_powerpoint)
            open(11,file=nupath//'Pk_cb_'//trim(adjustl(str_z))//'.txt',form='formatted')
            read(11,*) xi_cdm(6,1:)
            close(11)

            i = 1
            do while(minval(xitile_avange(2,:))*1.2 > kh_lin(i))
               i = i+1
            enddo

            j = 1
            do while(kh_lin(j) < k_std_max/4)
               j = j+1
            enddo
            if(head) print*,'nan in xi_cdm_std use CAMB',str_z,k_std_max/4
            print*,'index',i,j
            print*,'k    ',kh_lin(i),kh_lin(j)

            do while(kh_lin(i) < k_std_max)
               k_need_log = kh_lin_log(i)

               i1=1; i2=nnbin
               do while (i2-i1>1)
                  i_mid=(i1+i2)/2
                  if (k_need_log>xi_log(2,i_mid)) then
                     i1=i_mid
                  else
                     i2=i_mid
                  endif
               enddo

               if (k_need_log-xi_log(2,i1) < xi_log(2,i2)-k_need_log) then
                  i3 = i2
                  i2 = i1
                  i1 = i1-1
               else
                  i3 = i2+1
               endif

               if (i3 > nnbin-1) then
                  i1 = i1-1
                  i2 = i2-1
                  i3 = i3-1
               endif

               if (i1 < 1) then
                  i1 = i1+1
                  i2 = i2+1
                  i3 = i3+1
               endif

               a1 = (xi_log(:,i1)-xi_log(:,i2))/(xi_log(2,i1)-xi_log(2,i2))
               a2 = (xi_log(:,i2)-xi_log(:,i3))/(xi_log(2,i2)-xi_log(2,i3))
               k2 = (a1-a2)/(xi_log(2,i1)-xi_log(2,i3))
               b2 = a1-k2*(xi_log(2,i1)+xi_log(2,i2))
               c2 = xi_log(:,i1)-k2*xi_log(2,i1)**2-b2*xi_log(2,i1)
               xi_real = exp(k2*k_need_log**2+b2*k_need_log+c2)
               xi_cdm(5,i) = (1.15**(j-i)*xi_real(5)+1.15**(i-j)*xi_cdm(6,i))/(1.15**(j-i)+1.15**(i-j))
               if(ieee_is_nan(xi_cdm(5,i))) then
                  if(head) then
                     print*,''
                     print*,''
                     print*,i,i1,i2,i3,nnbin
                     print*,'interp',a1(5),a2(5),k2(5),b2(5),c2(5),xi_log(5,i1),xi_log(5,i2),xi_log(5,i3)
                     print*,'xi_std',xi_real(5),xitile_avange(5,i1),xitile_avange(5,i2),xitile_avange(5,i3)
                     print*,'xi_cdm',xi_cdm(5,i)
                     print*,'k_std',xitile_avange(2,i1),xitile_avange(2,i2),xitile_avange(2,i3)
                     print*,'kh',kh_lin(i),k_std_max
                     print*,'xi_standard',xitile_avange(5,:i3+10)
                     print*,'k_standard ',xitile_avange(2,:i3+10)
                  endif
                  error stop 'standard_power err '
               endif
               xi_cdm(4,i) = (1.15**(i-j)*xi_real(4)+1.15**(j-i)*xi_cdm(6,i))/(1.15**(i-j)+1.15**(j-i))
               xi_cdm(3,i) = (1.15**(i-j)*xi_real(3)+1.15**(j-i)*xi_cdm(6,i))/(1.15**(i-j)+1.15**(j-i))
               xi_cdm(2,i) = kh_lin(i)
               xi_cdm(1,i) = xi_cdm(1,i)+xi_real(1)

               i = i+1
            enddo

            xi_cdm(5,i-1:) = xi_cdm(6,i-1:)
            xi_cdm(4,i-1:) = xi_cdm(6,i-1:)
            xi_cdm(3,i-1:) = xi_cdm(6,i-1:)
            xi_cdm(2,i-1:) = kh_lin(i-1:)
            xi_cdm(1,i-1:) = 0


         else

            i = 1
            do while(minval(xitile_avange(2,:))/1.2 > kh_lin(i))
               i = i+1
            enddo

            j = 1
            do while(kh_lin(j) < nc/box*2*pi)
               j = j+1
            enddo


            k = (j-i)/2
            if(head) print*,'   std start '
            if(head) print*,'    ',kh_lin(i),kh_lin(j),k_std_max,minval(xitile_avange(2,:))

            do while (kh_lin(i) <= k_std_max)
               k_need_log = kh_lin_log(i)

               i1=1; i2=nnbin
               do while (i2-i1>1)
                  i_mid=(i1+i2)/2
                  if (k_need_log>xi_log(2,i_mid)) then
                     i1=i_mid
                  else
                     i2=i_mid
                  endif
               enddo

               if (k_need_log-xi_log(2,i1) < xi_log(2,i2)-k_need_log) then
                  i3 = i2
                  i2 = i1
                  i1 = i1-1
               else
                  i3 = i2+1
               endif

               if (i3 > nnbin-1) then
                  i1 = i1-1
                  i2 = i2-1
                  i3 = i3-1
               endif

               if (i1 < 1) then
                  i1 = i1+1
                  i2 = i2+1
                  i3 = i3+1
               endif

               a1 = (xi_log(:,i1)-xi_log(:,i2))/(xi_log(2,i1)-xi_log(2,i2))
               a2 = (xi_log(:,i2)-xi_log(:,i3))/(xi_log(2,i2)-xi_log(2,i3))
               k2 = (a1-a2)/(xi_log(2,i1)-xi_log(2,i3))
               b2 = a1-k2*(xi_log(2,i1)+xi_log(2,i2))
               c2 = xi_log(:,i1)-k2*xi_log(2,i1)**2-b2*xi_log(2,i1)
               xi_real = exp(k2*k_need_log**2+b2*k_need_log+c2)
               if (i <= j) then
                  xi_cdm(3:,i) = (1.15**(i-j+k)*ratio_cs**3*xi_real(1)*xi_real(3:)+1.15**(j-i-k)*xi_cdm(1,i)*xi_cdm(3:,i))/(1.15**(i-j+k)*ratio_cs**3*xi_real(1)+1.15**(j-i-k)*xi_cdm(1,i))
                  ! print*,'a',xitile_avange(2,i1),i1,i2,i3
                  ! print*,kh_lin(i),xi_cdm(5,i),xitile_avange(5,i1),xi_real(5),abs(xi_cdm(5,i)/xi_real(5))
               else
                  xi_cdm(3:,i) = xi_real(3:)
                  ! print*,'i',xitile_avange(2,i1),i1,i2,i3
                  ! print*,kh_lin(i),xi_cdm(5,i),xitile_avange(5,i1),xi_real(5)
               endif
               xi_cdm(1,i) = xi_cdm(1,i)+xi_real(1)
               xi_cdm(2,i) = kh_lin(i)

               i = i+1
            enddo
            !  print*,i,i1,i2,i3

            xi_cdm_log = log(xi_cdm)

            i1 = i-1
            do while(isnan(xi_cdm_log(5,i1)))
               i1 = i1-1
            enddo
            k_std_max=kh_lin(i1)

            i2 = i1
            do while(kh_lin(i2) > k_std_max/2)
               i2 = i2-1
            enddo

            i3 = i2
            do while(kh_lin(i3) > k_std_max/3)
               i3 = i3-1
            enddo

            print*,i,i1,i2,i3,k_std_max
            print*,kh_lin(i),kh_lin(i1),kh_lin(i2),kh_lin(i3)
            print*,xi_cdm_log(5,i),xi_cdm_log(5,i1),xi_cdm_log(5,i2),xi_cdm_log(5,i3)
            print*,xi_cdm(5,i),xi_cdm(5,i1),xi_cdm(5,i2),xi_cdm(5,i3)
            ! if(head) print*,'k_std',xitile_avange(2,i3),kh_lin(i),k_std_max
            a1 = (xi_cdm_log(:,i1)-xi_cdm_log(:,i2))/(xi_cdm_log(2,i1)-xi_cdm_log(2,i2))
            a2 = (xi_cdm_log(:,i2)-xi_cdm_log(:,i3))/(xi_cdm_log(2,i2)-xi_cdm_log(2,i3))
            k2 = (a1-a2)/(xi_cdm_log(2,i1)-xi_cdm_log(2,i3))
            b2 = a1-k2*(xi_cdm_log(2,i1)+xi_cdm_log(2,i2))
            c2 = xi_cdm_log(:,i1)-k2*xi_cdm_log(2,i1)**2-b2*xi_cdm_log(2,i1)
            i=i1
            do while (i <= npbin)
               k_need_log = kh_lin_log(i)
               xi_real = exp(k2*k_need_log**2+b2*k_need_log+c2)
               xi_cdm(:,i) = xi_real(:)
               xi_cdm(2,i) = kh_lin(i)
               i = i+1
            enddo
         endif




         ! write(str_i,'(i6)') image
         ! write(str_z,'(f7.3)') 1/sim%a-1
         ! print*,'   Write cicpower_segment into'
         ! print*,'     ',nupath//'Pk/'//trim(adjustl(str_z))//'*.bin'
         ! open(15,file=nupath//'Pk/'//trim(adjustl(str_z))//'_cicpower_sall.bin',status='replace',access='stream')
         ! write(15) xi_cdm(:,1:npbin)
         ! close(15)
         ! close(15)
         ! ! open(15,file=nupath//'z_step.txt',status='old',access='stream',position='append')
         ! ! write(15) str_z
         ! ! close(15)
         ! ! stop
      endif
      if(head .and. any(ieee_is_nan(xi_cdm(5,1:)))) then
         print*,''
         print*,''
         print*,''
         print*,'standard_power err ',i,i1,i2,i3
         print*,'xi_standard',xitile_avange(5,:)
         print*,'k_standard ',xitile_avange(2,:)
         print*,'xi_cdm',xi_cdm(5,:i)
         print*,'kh',kh_lin(:i)
         error stop 'standard_power err '
      endif
      sync all
      xi_cdm(5,1:) = xi_cdm(5,1:)[1]
      use_CAMB = use_CAMB[1]
      ! print*,'std power ',image,' done'
      sync all

   endsubroutine standard_power

   subroutine fine_power
      use variables
      implicit none
      ! save

      integer i,j,k,l,ilayer,n,nnsc,idx1(3),idx2(3),ix,iy,iz,i1,i2,i3,itile,cur_checkpoint,tcpu
      integer(8) nlast,ip,nplocal,npglobal,nzero,np,inest
      real pos1(3),dx1(3),dx2(3),xpos(ndim),nc1(3),nc2(3)
      real rho_th(-nfb(cic_iapm):nfp(cic_iapm)+nfb(cic_iapm)+1,-nfb(cic_iapm):nfp(cic_iapm)+nfb(cic_iapm)+1,-nfb(cic_iapm):nfp(cic_iapm)+nfb(cic_iapm)+1)
      ! real(8) klast,dk
      integer,parameter :: nlayer=3
      real xi_log(6,0:nnbin),xi_real(6),xi_cdm_log(6,0:npbin)
      real a1(6),a2(6),k2(6),b2(6),c2(6)
      integer(4) i_mid
      real k_need_log,a,k_fine_max,k_std_max


      if (head) print*, ' fine_power'
      xitile_avange = 0
      nnsc = 0
      ! sync all



      call system_clock(tt1,t_rate)
      !$omp paralleldo default(shared) num_threads(nteam) schedule(dynamic,1)&
      !$omp& private(n,ipm3,ix,iy,iz,ilayer,nlast,iteam,rho_th,str_i,str_z,nc1,nc2)
      ! do n=22,22
      do n=1,(nns*nnt)**3
         iteam=omp_get_thread_num()+1
         xitile(:,0:,iteam) = 0
         ! ipm3 = nspk(n)
         ! if (head) print*, n,ipm3,iteam
         ipm3 = n
         ix = ixyz3(4,ipm3)
         iy = ixyz3(5,ipm3)
         iz = ixyz3(6,ipm3)
         nc1=(ixyz3(1:3,ipm3)-1)*ntt
         nc2=ixyz3(1:3,ipm3)*ntt+1

         rho_th = 0
         do ilayer=0,nlayer-1
            !$omp paralleldo default(shared) num_threads(nnest) schedule(dynamic,1)&
            !$omp& private(np,nzero,k,j,i,ip,pos1,dx1,dx2,idx1,idx2)
            do k=nc1(3)+ilayer,nc2(3),nlayer
               inest=omp_get_thread_num()+1
               do j=nc1(2),nc2(2)
                  do i=nc1(1),nc2(1)
                     np=rhoc(i,j,k,ix,iy,iz)
                     nzero=idx_b_r(j,k,ix,iy,iz)-sum(rhoc(i:,j,k,ix,iy,iz))
                     do l=1,np
                        ip=nzero+l
                        if (ip<1 )then
                           print*,'ip out of range in standard_power',image
                           print*,image,'tile',ix,iy,iz
                           print*,image,'ip',ip
                           print*,image,'subtile',i,j,k,nzero
                           print*,image,'rhoc',maxval(rhoc(:,:,:,ix,iy,iz)),minval(rhoc(:,:,:,ix,iy,iz))
                           error stop
                        endif
                        pos1=((/i,j,k/)-1) + (int(xp(:,ip)+ishift,izipx)+rshift)*x_resolution
                        pos1=pos1*ratio_cs*ratio_sf(cic_iapm) - (ixyz3(1:3,ipm3)-1)*nfp(cic_iapm)-0.5
                        ! print*,i,j,k,np,l,pos1
                        idx1=floor(pos1)+1; idx2=idx1+1
                        dx1=idx1-pos1;      dx2=1-dx1
                        if (minval(idx1)< -nfb(cic_iapm) .or. maxval(idx2)>nfp(cic_iapm)+nfb(cic_iapm)+1) then
                           print*,'mass assignment out of range in standard_power'
                           print*,image,'tile',ixyz3(:,ipm3),ipm3
                           print*,image,'range',-nfb(cic_iapm),nfp(cic_iapm)+nfb(cic_iapm)+1
                           print*,image,'xp', xp(:,ip)
                           print*,image,'xp_real',  ([i,j,k]-1)+(int(xp(:,ip)+ishift,izipx)+rshift)*x_resolution
                           print*,image,'xpos',pos1
                           print*,image,'idx',idx1,idx2
                           print*,image,'pars',ratio_cs,ratio_sf(cic_iapm),nfp(cic_iapm),(ixyz3(1:3,ipm3)-1)
                           error stop
                        endif
                        rho_th(idx1(1),idx1(2),idx1(3))=rho_th(idx1(1),idx1(2),idx1(3))+dx1(1)*dx1(2)*dx1(3)
                        rho_th(idx2(1),idx1(2),idx1(3))=rho_th(idx2(1),idx1(2),idx1(3))+dx2(1)*dx1(2)*dx1(3)
                        rho_th(idx1(1),idx2(2),idx1(3))=rho_th(idx1(1),idx2(2),idx1(3))+dx1(1)*dx2(2)*dx1(3)
                        rho_th(idx1(1),idx1(2),idx2(3))=rho_th(idx1(1),idx1(2),idx2(3))+dx1(1)*dx1(2)*dx2(3)
                        rho_th(idx1(1),idx2(2),idx2(3))=rho_th(idx1(1),idx2(2),idx2(3))+dx1(1)*dx2(2)*dx2(3)
                        rho_th(idx2(1),idx1(2),idx2(3))=rho_th(idx2(1),idx1(2),idx2(3))+dx2(1)*dx1(2)*dx2(3)
                        rho_th(idx2(1),idx2(2),idx1(3))=rho_th(idx2(1),idx2(2),idx1(3))+dx2(1)*dx2(2)*dx1(3)
                        rho_th(idx2(1),idx2(2),idx2(3))=rho_th(idx2(1),idx2(2),idx2(3))+dx2(1)*dx2(2)*dx2(3)
                     enddo
                  enddo
               enddo
            enddo
            !$omp endparalleldo
         enddo

         do i=1-nfb(cic_iapm),nfp(cic_iapm)+nfb(cic_iapm)
            rho_th(:,:,i)=rho_th(:,:,i)/(rho8/ratio_cf(cic_iapm)/ratio_cf(cic_iapm)/ratio_cf(cic_iapm))-1
         enddo
         nlast = sum(rhoc(nc1(1)+1:nc2(1)-1,nc1(2)+1:nc2(2)-1,nc1(3)+1:nc2(3)-1,ix,iy,iz))
         ! print*,iapm,nlast,sum(rho_th(1:nfp(cic_iapm),1:nfp(cic_iapm),1:nfp(cic_iapm)))

         rho3p(:,:,:,iteam) = 0
         rho3p(1:nfp(cic_iapm),1:nfp(cic_iapm),1:nfp(cic_iapm),iteam) = rho_th(1:nfp(cic_iapm),1:nfp(cic_iapm),1:nfp(cic_iapm))
         call sfftw_execute(plan3p(iteam))
         call cicpower_tile(xitile(:,0:,iteam),rho3kp(:,:,:,iteam),nfp(cic_iapm),nlast)
         ! selectcase(cic_iapm)
         !  case(2)
         !    rho3_2(:,:,:,iteam) = 0
         !    rho3_2(1:nfp(cic_iapm),1:nfp(cic_iapm),1:nfp(cic_iapm),iteam) = rho_th(1:nfp(cic_iapm),1:nfp(cic_iapm),1:nfp(cic_iapm))
         !    call sfftw_execute(plan3p(iteam,cic_iapm))
         !    call cicpower_tile(xitile(:,0:,iteam),rho3k_2(1:nfp(cic_iapm)/2+1,1:nfp(cic_iapm),1:nfp(cic_iapm),iteam),nfp(cic_iapm),nlast)
         !  case(3)
         !    rho3_4(:,:,:,iteam) = 0
         !    rho3_4(1:nfp(cic_iapm),1:nfp(cic_iapm),1:nfp(cic_iapm),iteam) = rho_th(1:nfp(cic_iapm),1:nfp(cic_iapm),1:nfp(cic_iapm))
         !    call sfftw_execute(plan3p(iteam,cic_iapm))
         !    call cicpower_tile(xitile(:,0:,iteam),rho3k_4(1:nfp(cic_iapm)/2+1,1:nfp(cic_iapm),1:nfp(cic_iapm),iteam),nfp(cic_iapm),nlast)
         !  case(4)
         !    rho3_6(:,:,:,iteam) = 0
         !    rho3_6(1:nfp(cic_iapm),1:nfp(cic_iapm),1:nfp(cic_iapm),iteam) = rho_th(1:nfp(cic_iapm),1:nfp(cic_iapm),1:nfp(cic_iapm))
         !    call sfftw_execute(plan3p(iteam,cic_iapm))
         !    call cicpower_tile(xitile(:,0:,iteam),rho3k_6(1:nfp(cic_iapm)/2+1,1:nfp(cic_iapm),1:nfp(cic_iapm),iteam),nfp(cic_iapm),nlast)
         !  case(5)
         !    rho3_8(:,:,:,iteam) = 0
         !    rho3_8(1:nfp(cic_iapm),1:nfp(cic_iapm),1:nfp(cic_iapm),iteam) = rho_th(1:nfp(cic_iapm),1:nfp(cic_iapm),1:nfp(cic_iapm))
         !    call sfftw_execute(plan3p(iteam,cic_iapm))
         !    call cicpower_tile(xitile(:,0:,iteam),rho3k_8(1:nfp(cic_iapm)/2+1,1:nfp(cic_iapm),1:nfp(cic_iapm),iteam),nfp(cic_iapm),nlast)
         ! endselect

         !$omp barrier
         if (iteam == 1) then
            do i=1,nteam
               xitile_avange(:,0:) = xitile_avange(:,0:)+xitile(:,0:,i)
               nnsc = nnsc +1
            enddo
         endif
         !$omp barrier
         ! stop
      enddo
      !$omp endparalleldo


      xitile_avange = xitile_avange/nnsc

      sync all
      if (head) then
         do i = 2,nn**3
            xitile_avange(:,0:)[1] = xitile_avange(:,0:)[1] + xitile_avange(:,0:)[i]
         enddo
         xitile_avange = xitile_avange/(nn**3)

         a = 0
         xitile_avange(2,:)=xitile_avange(2,:)/xitile_avange(1,:)
         xitile_avange(3,:)=xitile_avange(3,:)/xitile_avange(1,:) ! raw power
         xitile_avange(4,:)=xitile_avange(4,:)/xitile_avange(1,:) ! raw power - Dk
         xitile_avange(5,:)=xitile_avange(4,:)
         ! print*,'xi_avg',xitile_avange(5,:)
         ! print*,'xi_avg',xitile_avange(3,:)
         call pk_correction_tile(xitile_avange,2,3,nfp(cic_iapm)/2,a)
         call pk_correction_tile(xitile_avange,2,3,nfp(cic_iapm)/2,a)
         call pk_correction_tile(xitile_avange,2,3,nfp(cic_iapm)/2,a)
         call pk_correction_tile(xitile_avange,2,3,nfp(cic_iapm)/2,a)
         ! divide and normalize
         xitile_avange(2,:)=xitile_avange(2,:)*(2*pi)/subtile ! k_phys
         xitile_avange(3,:)=xitile_avange(3,:)*(subtile**3) ! power_phys
         xitile_avange(4,:)=xitile_avange(4,:)*(subtile**3) ! power_phys
         xitile_avange(5,:)=xitile_avange(5,:)*(subtile**3) ! power_phys

         ! print*,'xi_corr',xitile_avange(5,:)
         ! print*,'xi_corr',xitile_avange(3,:)

         k_fine_max = nfp(cic_iapm)/subtile*pi
         k_std_max = ngp/tile*pi
         write(str_i,'(i6)') image
         write(str_z,'(f7.3)') 1/sim%a-1
         print*,'Write cicpower_fine into',opath//'image'//trim(adjustl(str_i))//'/'//trim(adjustl(str_z))//'_cicpower_finec.bin'
         open(15,file=nupath//'Pk/'//trim(adjustl(str_z))//'_cicpower_fine.bin',status='replace',access='stream')
         write(15) xitile_avange(:,1:)
         close(15)

         ! print*, '********************** fine_segment *******************'

         j = 1
         do while(xitile_avange(2,j) < k_std_max*2.)
            j = j+1
         enddo
         print*,'k_fine_index',j


         xi_log = log(xitile_avange)
         if(any(ieee_is_nan(xi_log(5,1:j-1)))) then

            use_CAMB = .true.

            xi_log(5,1:) = xi_log(3,1:)
            if(head) print*,'nan in xi_cdm_fine use CAMB'

            write(str_z,'(f8.4)') z_powerpoint(sim%cur_powerpoint)
            open(11,file=nupath//'Pk_cb_'//trim(adjustl(str_z))//'.txt',form='formatted')
            read(11,*) xi_cdm(6,1:)
            close(11)
            
            i = 1
            do while(minval(xitile_avange(2,:))*1.2 > kh_lin(i))
               i = i+1
            enddo

            j = 1
            do while(kh_lin(j) < k_fine_max/4)
               j = j+1
            enddo

            do while(kh_lin(i) <= k_fine_max)
               k_need_log = kh_lin_log(i)

               i1=1; i2=nnbin
               do while (i2-i1>1)
                  i_mid=(i1+i2)/2
                  if (k_need_log>xi_log(2,i_mid)) then
                     i1=i_mid
                  else
                     i2=i_mid
                  endif
               enddo

               if (k_need_log-xi_log(2,i1) < xi_log(2,i2)-k_need_log) then
                  i3 = i2
                  i2 = i1
                  i1 = i1-1
               else
                  i3 = i2+1
               endif

               if (i3 > nnbin-1) then
                  i1 = i1-1
                  i2 = i2-1
                  i3 = i3-1
               endif

               if (i1 < 1) then
                  i1 = i1+1
                  i2 = i2+1
                  i3 = i3+1
               endif

               a1 = (xi_log(:,i1)-xi_log(:,i2))/(xi_log(2,i1)-xi_log(2,i2))
               a2 = (xi_log(:,i2)-xi_log(:,i3))/(xi_log(2,i2)-xi_log(2,i3))
               k2 = (a1-a2)/(xi_log(2,i1)-xi_log(2,i3))
               b2 = a1-k2*(xi_log(2,i1)+xi_log(2,i2))
               c2 = xi_log(:,i1)-k2*xi_log(2,i1)**2-b2*xi_log(2,i1)
               xi_real = exp(k2*k_need_log**2+b2*k_need_log+c2)
               xi_cdm(5,i) = (1.15**(j-i)*xi_real(5)+1.15**(i-j)*xi_cdm(6,i))/(1.15**(j-i)+1.15**(i-j))
               xi_cdm(4,i) = (1.15**(j-i)*xi_real(4)+1.15**(i-j)*xi_cdm(6,i))/(1.15**(j-i)+1.15**(i-j))
               xi_cdm(3,i) = (1.15**(j-i)*xi_real(3)+1.15**(i-j)*xi_cdm(6,i))/(1.15**(j-i)+1.15**(i-j))
               xi_cdm(2,i) = kh_lin(i)
               xi_cdm(1,i) = xi_cdm(1,i)+xi_real(1)

               i = i+1
            enddo

            xi_cdm(5,i-1:) = xi_cdm(6,i-1:)
            xi_cdm(4,i-1:) = xi_cdm(6,i-1:)
            xi_cdm(3,i-1:) = xi_cdm(6,i-1:)
            xi_cdm(2,i-1:) = kh_lin(i-1:)
            xi_cdm(1,i-1:) = 0


         else
            i = 1
            do while(xitile_avange(2,i) < k_std_max)
               i = i+1
            enddo

            j = 1
            do while(xitile_avange(2,j) < k_std_max*2.)
               j = j+1
            enddo

            ! print*,k_std_max,k_fine_max,k_fine_max*1.731
            if (xitile_avange(5,j)*2 > xitile_avange(5,i)) k_fine_max = k_std_max ! grid effect


            k = 1
            do while(xitile_avange(1,k) <= 50/(nn**3)/nnsc)
               k = k+1
            enddo


            i = 1
            do while(xitile_avange(2,k) > kh_lin(i))
               i = i+1
            enddo

            ! print*,'xi_se',xitile_avange(2,:)
            ! l = 1
            ! do while(l < nnbin .and. (isnan(xitile_avange(2,l))))
            !    l = l+1
            ! enddo

            k = (j-i)/2

            print*,'   fine start '
            print*,'     ',kh_lin(i),kh_lin(j),k_fine_max,k

            do while (kh_lin(i) <= k_fine_max)
               k_need_log = kh_lin_log(i)

               i1=1; i2=nnbin
               do while (i2-i1>1)
                  i_mid=(i1+i2)/2
                  if (k_need_log>xi_log(2,i_mid)) then
                     i1=i_mid
                  else
                     i2=i_mid
                  endif
               enddo

               if (k_need_log-xi_log(2,i1) < xi_log(2,i2)-k_need_log) then
                  i3 = i2
                  i2 = i1
                  i1 = i1-1
               else
                  i3 = i2+1
               endif

               if (i3 > nnbin-1) then
                  i1 = i1-1
                  i2 = i2-1
                  i3 = i3-1
               endif

               if (i1 < 1) then
                  i1 = i1+1
                  i2 = i2+1
                  i3 = i3+1
               endif

               a1 = (xi_log(:,i1)-xi_log(:,i2))/(xi_log(2,i1)-xi_log(2,i2))
               a2 = (xi_log(:,i2)-xi_log(:,i3))/(xi_log(2,i2)-xi_log(2,i3))
               k2 = (a1-a2)/(xi_log(2,i1)-xi_log(2,i3))
               b2 = a1-k2*(xi_log(2,i1)+xi_log(2,i2))
               c2 = xi_log(:,i1)-k2*xi_log(2,i1)**2-b2*xi_log(2,i1)
               xi_real = exp(k2*k_need_log**2+b2*k_need_log+c2)
               if (i <= j) then
                  xi_cdm(3:,i) = (1.15**(i-j+k)*ratio_cs**3*xi_real(1)*xi_real(3:)+1.15**(j-i-k)*xi_cdm(1,i)*xi_cdm(3:,i))/(1.15**(i-j+k)*ratio_cs**3*xi_real(1)+1.15**(j-i-k)*xi_cdm(1,i))
                  ! print*,kh_lin(i),'a',1.15**(i-j+k),1.15**(j-i-k),1.15**(i-j+k)*ratio_cs**3*xi_real(1)/(1.15**(i-j+k)*ratio_cs**3*xi_real(1)+1.15**(j-i-k)*xi_cdm(1,i))
               else
                  xi_cdm(3:,i) = xi_real(3:)
                  ! print*,'i',xitile_avange(2,i1),i1,i2,i3
                  ! print*,kh_lin(i),xi_cdm(5,i),xitile_avange(5,i1),xi_real(5)
               endif
               if (ieee_is_nan(xi_cdm(5,i))) then
                  print*,''
                  print*,''
                  print*,''
                  print*,'fine_xi_real err ',i,i1,i2,i3
                  print*,'k',xitile_avange(2,i),xitile_avange(2,i1),xitile_avange(2,i2),xitile_avange(2,i3)
                  print*,'pk',xitile_avange(5,i),xitile_avange(5,i1),xitile_avange(5,i2),xitile_avange(5,i3)
                  print*,'interp:',a1(5),a2(5),k2(5),b2(5),c2(5)
                  print*,kh_lin(i),xi_cdm(5,i),xitile_avange(5,i1),xi_real(5)
                  error stop 'fine_xi_real err '
               endif
               xi_cdm(1,i) = xi_cdm(1,i)+xi_real(1)
               xi_cdm(2,i) = kh_lin(i)

               i = i+1
            enddo
            !  print*,i,i1,i2,i3

            xi_cdm_log = log(xi_cdm)

            if (k_fine_max == k_std_max) then
               ! print*,'grid'
               i1 = i-1
               do while(kh_lin(i1) > k_std_max .and. (isnan(xi_cdm_log(5,i1))))
                  i1 = i1-1
               enddo
            else
               i1 = i-1
               do while(kh_lin(i1) > ngp*nnt/box*pi .and. (isnan(xi_cdm_log(5,i1))))
                  i1 = i1-1
               enddo
            endif

            i2 = i1
            do while(kh_lin(i2) > kh_lin(i1)/3*2)
               i2 = i2-1
            enddo

            i3 = i2
            do while(kh_lin(i3) > kh_lin(i1)/4)
               i3 = i3-1
            enddo

            print*,i,i1,i2,i3,k_fine_max
            print*,kh_lin(i),kh_lin(i1),kh_lin(i2),kh_lin(i3)
            print*,xi_cdm_log(5,i),xi_cdm_log(5,i1),xi_cdm_log(5,i2),xi_cdm_log(5,i3)
            print*,xi_cdm(5,i),xi_cdm(5,i1),xi_cdm(5,i2),xi_cdm(5,i3)
            i=i1
            a1 = (xi_cdm_log(:,i1)-xi_cdm_log(:,i2))/(xi_cdm_log(2,i1)-xi_cdm_log(2,i2))
            a2 = (xi_cdm_log(:,i2)-xi_cdm_log(:,i3))/(xi_cdm_log(2,i2)-xi_cdm_log(2,i3))
            k2 = (a1-a2)/(xi_cdm_log(2,i1)-xi_cdm_log(2,i3))
            b2 = a1-k2*(xi_cdm_log(2,i1)+xi_cdm_log(2,i2))
            c2 = xi_cdm_log(:,i1)-k2*xi_cdm_log(2,i1)**2-b2*xi_cdm_log(2,i1)
            do while (i <= npbin)
               k_need_log = kh_lin_log(i)
               xi_real = exp(k2*k_need_log**2+b2*k_need_log+c2)
               xi_cdm(:,i) = xi_real(:)
               xi_cdm(2,i) = kh_lin(i)
               i = i+1
            enddo
         endif
      endif
      if(head .and. any(ieee_is_nan(xi_cdm(5,1:)))) then
         print*,''
         print*,''
         print*,''
         print*,'fine_power err ',i,i1,i2,i3
         print*,'xi_fine',xitile_avange(5,:)
         print*,'k_fine ',xitile_avange(2,:)
         print*,'xi_cdm',xi_cdm(5,:i)
         print*,'kh',kh_lin(:i)
         error stop 'fine_power err '
      endif
      sync all
      xi_cdm(5,1:) = xi_cdm(5,1:)[1]
      use_CAMB = use_CAMB[1]
      sync all
   endsubroutine fine_power

endmodule
