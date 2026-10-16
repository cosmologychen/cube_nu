
module cicpower_global
   use omp_lib
   use cubefft
   use pencil_fft
   use parameters
   use ieee_arithmetic
   implicit none
   save

   logical use_CAMB[*]
   real xiglobal(6,0:npbin)[*],xi_cdm(6,0:npbin)[*],rho8[*]


contains


   real function interp1_global(xdata,ydata,xq)
      implicit none
      integer(4) i_mid,i1,i2
      real xdata(npbin),ydata(npbin),xq
      i1=1; i2=npbin
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
      real alpha,kvec(3),kmag,kmagn,kvecn(3),ks(3),Wk2Pk,Pk,cdata(0:npbin,0:ncore,3),xi_g(6,0:npbin)[*],kmaxs(26)

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
         cdata(1:npbin,0,:)=sum(cdata(1:npbin,1:ncore,:),dim=2)
         cdata(1:npbin,0,2)=cdata(1:npbin,0,2)/cdata(1:npbin,0,1)
         cdata(1:npbin,0,3)=cdata(1:npbin,0,3)/cdata(1:npbin,0,1)
         xi_g(5,1:npbin)=xi_g(4,1:npbin)/cdata(1:npbin,0,3)
      endif
   endsubroutine

   subroutine power_global(xi_global,n_particle)
      use omp_lib
      use pencil_fft_global
      implicit none
      integer i,j,k,ig,jg,kg,ibin
      integer(8) n_particle
      real a,kr,kx(3),C1k(3),Dk,amp11,xi_global(6,0:npbin)[*],kmaxs(26)


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
               if (ibin > npbin .and. ibin < 0 ) then
                  print*,image,'kr out of range'
                  print*,'kr',i,j,k,kg,jg,ig,kr
                  print*,'npbin',npbin
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
         ! print*,xi_global(3,npbin)
         xi_global(2,:)=xi_global(2,:)/xi_global(1,:)
         xi_global(3,:)=xi_global(3,:)/xi_global(1,:) ! raw power
         ! print*,xi_global(3,npbin),xi_global(1,npbin)
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
      !       print*,'power_global err '
      !       print*,'k_global',xiglobal(2,:)
      !       print*,'xi_global_raw',xi_global(3,:)
      !       print*,'xi_global_corr',xi_global(5,:)
      !    endif
      !    error stop 'power_global err '
      ! endif
   endsubroutine power_global

   subroutine global_power
      use variables
      use pencil_fft_global
      implicit none
      ! save

      integer i,j,k,l,ilayer,n,idx1(3),idx2(3),ix,iy,iz,i1,i2,i3,itile,cur_checkpoint
      integer(8) ip,nplocal,npglobal,nzero,np
      real pos1(3),dx1(3),dx2(3),xpos(ndim),nbb
      real(4),allocatable :: rho_grid(:,:,:)[:,:,:]
      ! real(8),save :: rho8[*]
      real xi_log(6,0:npbin),kmin,xi_real(6),xi_cdm_log(6,0:npbin)
      real a1(6),a2(6),k2(6),b2(6),c2(6)
      integer(4) i_mid
      real k_need_log,a,k_global_max,k_std_max


      allocate(rho1_global(nfg,nfg,nfg)[nn,nn,*],rho_grid(0:nfg+1,0:nfg+1,0:nfg+1)[nn,nn,*])
      ! print*, ' global_power',image
      sync all


      if (head) print*, ' global_power'
      ! if (head) print*,'size_rhoc',size(rho_grid)
      rho_grid=0; nzero=0; nbb = nfg
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
#ifdef ZIPX
                           pos1=nt*((/itx,ity,itz/)-1)+ ((/i,j,k/)-1) + (int(xp(:,ip)+ishift,izipx)+rshift)*x_resolution
                           pos1=pos1*real(nfg)/real(nc) - 0.5
#else
                           pos1=xp(:,ip)*real(nfg)/real(ng) - 0.5
                           if ( pos1(1) == nbb ) pos1(1)= 0
                           if ( pos1(2) == nbb ) pos1(2)= 0
                           if ( pos1(3) == nbb ) pos1(3)= 0
#endif
                           idx1=floor(pos1)+1; idx2=idx1+1
                           dx1=idx1-pos1;      dx2=1-dx1
                           if (minval(idx1)< 0 .or. maxval(idx2)>nfg+1) then
                              print*,'mass assignment out of range in standard_power'
                              print*,image,'tile',ix,iy,iz,i,j,k,ipm2
                              print*,image,'range',0,nfg+1
                              print*,image,'xp', xp(:,ip)
#ifdef ZIPX
                              print*,image,'xp_real',  ([i,j,k]-1)+(int(xp(:,ip)+ishift,izipx)+rshift)*x_resolution
#endif
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
      block
         real,allocatable :: sbuf_x(:,:)[:,:,:]
         real,allocatable :: sbuf_y(:,:)[:,:,:]
         real,allocatable :: sbuf_z(:,:)[:,:,:]
         real :: rbuf_x(nfg+2,nfg+2)
         real :: rbuf_y(nfg+2,nfg+2)
         real :: rbuf_z(nfg+2,nfg+2)
         allocate(sbuf_x(nfg+2,nfg+2)[nn,nn,*])
         sync all
         sbuf_x = rho_grid(nfg+1,:,:)
         sync all
         rbuf_x = sbuf_x(:,:)[inx,icy,icz]
         rho_grid(1,:,:) = rho_grid(1,:,:) + rbuf_x
         sync all
         sbuf_x = rho_grid(0,:,:)
         sync all
         rbuf_x = sbuf_x(:,:)[ipx,icy,icz]
         rho_grid(nfg,:,:) = rho_grid(nfg,:,:) + rbuf_x
         sync all
         deallocate(sbuf_x)

         allocate(sbuf_y(nfg+2,nfg+2)[nn,nn,*])
         sync all
         sbuf_y = rho_grid(:,nfg+1,:)
         sync all
         rbuf_y = sbuf_y(:,:)[icx,iny,icz]
         rho_grid(:,1,:) = rho_grid(:,1,:) + rbuf_y
         sync all
         sbuf_y = rho_grid(:,0,:)
         sync all
         rbuf_y = sbuf_y(:,:)[icx,ipy,icz]
         rho_grid(:,nfg,:) = rho_grid(:,nfg,:) + rbuf_y
         sync all
         deallocate(sbuf_y)

         allocate(sbuf_z(nfg+2,nfg+2)[nn,nn,*])
         sync all
         sbuf_z = rho_grid(:,:,nfg+1)
         sync all
         rbuf_z = sbuf_z(:,:)[icx,icy,inz]
         rho_grid(:,:,1) = rho_grid(:,:,1) + rbuf_z
         sync all
         sbuf_z = rho_grid(:,:,0)
         sync all
         rbuf_z = sbuf_z(:,:)[icx,icy,ipz]
         rho_grid(:,:,nfg) = rho_grid(:,:,nfg) + rbuf_z
         sync all
         deallocate(sbuf_z)
      end block
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
      call power_global(xiglobal,sim%npglobal)

      sync all
      if (head) then
         xi_cdm = 0
         k_std_max = ng_global/box*pi
         k_global_max = nfg_global/box*pi

         write(str_z,'(f7.3)') z_powerpoint(sim%cur_powerpoint)
         print*,'Write cicpower_global into',nupath//'Pk/'//trim(adjustl(str_z))//'_cicpower_global.bin'
         open(15,file=nupath//'Pk/'//trim(adjustl(str_z))//'_cicpower_global.bin',status='replace',access='stream')
         write(15) xiglobal(:,0:npbin)[1]
         close(15)

         ! print*, '**********************global_segment *******************'
         j = 1
         do while(xiglobal(2,j) < k_std_max)
            j = j+1
         enddo

         if(any(ieee_is_nan(xiglobal(5,1:j-1))) .or. xiglobal(5,j-1) < xiglobal(3,j-1)) xiglobal(5,1:npbin) = xiglobal(3,1:npbin)

         i = 1
         do while(xiglobal(2,i) < k_std_max)
            i = i+1
         enddo

         j = 1
         do while(xiglobal(2,j) < k_std_max*2.)
            j = j+1
         enddo

         if (xiglobal(5,j)*2 > xiglobal(5,i)) then
            print*,'grid effect'
            k_std_max = k_std_max/2 ! grid effect
         else
            k_std_max = k_global_max
         endif

         xi_log = log(xiglobal)
         xi_cdm(:5,:26) = xiglobal(:5,:26)

         i = 26
         do while(kh_lin(i) <= k_std_max .and. i <= npbin)

            i = i+1
            k_need_log = kh_lin_log(i)

            i1=1; i2=npbin
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

            if (i3 > npbin) then
               print*,'out range of global bin'
               i1 = npbin-2
               i2 = npbin-1
               i3 = npbin
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
         enddo

         xi_cdm_log = log(xi_cdm)

         i1 = i-1

         i2 = i1
         do while(kh_lin(i2) > kh_lin(i1)/3)
            i2 = i2-1
         enddo

         i3 = i2
         do while(kh_lin(i3) > 8*pi/box)
            i3 = i3-1
         enddo

         ! print*,i,i1,i2,i3
         ! print*,k_global_max,kh_lin(i1),kh_lin(i1)/4*3,kh_lin(i1)/2
         ! print*,kh_lin(i),kh_lin(i1),kh_lin(i2),kh_lin(i3)
         ! print*,xi_cdm_log(5,i),xi_cdm_log(5,i1),xi_cdm_log(5,i2),xi_cdm_log(5,i3)
         ! print*,xi_cdm(5,i),xi_cdm(5,i1),xi_cdm(5,i2),xi_cdm(5,i3)
         i=i1
         a1 = (xi_cdm_log(:,i1)-xi_cdm_log(:,i2))/(xi_cdm_log(2,i1)-xi_cdm_log(2,i2))
         a2 = (xi_cdm_log(:,i2)-xi_cdm_log(:,i3))/(xi_cdm_log(2,i2)-xi_cdm_log(2,i3))
         k2 = (a1-a2)/(xi_cdm_log(2,i1)-xi_cdm_log(2,i3))
         b2 = a1-k2*(xi_cdm_log(2,i1)+xi_cdm_log(2,i2))
         c2 = xi_cdm_log(:,i1)-k2*xi_cdm_log(2,i1)**2-b2*xi_cdm_log(2,i1)
         do while (i <= npbin)
            i = i+1
            k_need_log = kh_lin_log(i)
            xi_real = exp(k2*k_need_log**2+b2*k_need_log+c2)
            if(ieee_is_nan(xi_real(5))) then
               print*,'global xi_real err',i,i1,i2,i3
               print*,xi_cdm(1:5,i1)
               print*,xi_cdm(1:5,i2)
               print*,xi_cdm(1:5,i3)
               print*,xi_cdm_log(1:5,i1)
               print*,xi_cdm_log(1:5,i2)
               print*,xi_cdm_log(1:5,i3)
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
         enddo

         write(str_z,'(f7.3)') z_powerpoint(sim%cur_powerpoint)
         print*,'   Write cicpower_global into'
         print*,'     ',nupath//'Pk/'//trim(adjustl(str_z))//'*.bin'
         print*,'  cicpower_global_step',xi_cdm(5,1),xi_cdm(5,npbin/3),xi_cdm(5,npbin/3*2),xi_cdm(5,npbin)
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

   endsubroutine global_power

endmodule
