!#define gadget
!#define RSD
#define merge_projection
! #define density_nu
! #define density_matter
! #define power_matter


program cicpower
   use omp_lib
   use parameters
   use pencil_fft
   use powerspectrum
   implicit none
   save

   integer i,j,k,l,np,idx(3),idx1(3),idx2(3),itx,ity,itz,cur_checkpoint,ix,iy,iz
   integer(8) nlast,ip,nplocal,npglobal
   character(20) str_z,str_i
   real pos1(3),dx1(3),dx2(3),xi(10,0:nbin)[*],nbb
   real(8) rho8[*]
   integer(4),allocatable :: rhoc(:,:,:,:,:,:)

#ifdef ZIPX
   integer(izipx),allocatable :: xp(:,:)
#else
   real(4),allocatable :: xp(:,:)
#endif

   real(4),allocatable :: xv(:,:),rho_grid(:,:,:)[:,:,:],rho_c(:,:,:)[:,:,:],rho_etc(:,:,:),proj_xy(:,:)
   real kh_lin(npbin),kh_lin_log(npbin),tf_F(npbin),tf_F_log(npbin),Pk_nu(npbin),Pk_nu_log(npbin)




   call geometry
   allocate(rhoc(nt,nt,nt,nnt,nnt,nnt),rho_grid(0:nw+1,0:nw+1,0:nw+1)[nn,nn,*],rho_c(nw,nw,nw)[nn,nn,*])
   if (head) then
      print*, 'cicpower on resolution: nw,nw_global=',nw,nw*nn
      print*, 'checkpoint at:'
      open(16,file='./z_checkpoint.txt',status='old')
      do i=1,nmax_redshift
         read(16,end=71,fmt='(f8.4)') z_checkpoint(i); print*, z_checkpoint(i)
      enddo
71    n_checkpoint=i-1
      close(16); print*,''
   endif
   sync all
   n_checkpoint=n_checkpoint[1]; z_checkpoint(:)=z_checkpoint(:)[1]
   call create_penfft_plan

   do cur_checkpoint= n_checkpoint,n_checkpoint
   ! do cur_checkpoint= 1,n_checkpoint
      sim%cur_checkpoint=cur_checkpoint
      if (head) print*, ''
      if (head) print*, '==========================================='
      if (head) print*, '==========================================='
      if (head) print*, 'Start analyzing redshift ',z2str(z_checkpoint(cur_checkpoint))
      !print*,output_name('info')
      open(11,file=output_name('info'),access='stream'); read(11) sim; close(11)
      if (sim%izipx/=izipx .or. sim%izipv/=izipv) error stop 'zip format incompatable'
      nplocal=sim%nplocal; npglobal=sim%npglobal
      if (head) print*, 'nplocal,npglobal =',nplocal,npglobal
      allocate(xp(3,nplocal))
      open(11,file=output_name('xp'),access='stream'); read(11) xp; close(11)
      open(11,file=output_name('np'),access='stream'); read(11) rhoc; close(11)
      rho_grid=0; nlast=0; nbb = nw
      do itz=1,nnt
         do ity=1,nnt
            do itx=1,nnt
               do k=1,nt
                  do j=1,nt
                     do i=1,nt
                        np=rhoc(i,j,k,itx,ity,itz)
                        do l=1,np
                           ip=nlast+l
#ifdef ZIPX
                           pos1=nt*((/itx,ity,itz/)-1)+ ((/i,j,k/)-1) + (int(xp(:,ip)+ishift,izipx)+rshift)*x_resolution
                           pos1=pos1*real(nw)/real(nc) - 0.5
#else
                           pos1 = xp(:,ip)*nic - 0.5
                           if ( pos1(1) == nbb ) pos1(1)= 0
                           if ( pos1(2) == nbb ) pos1(2)= 0
                           if ( pos1(3) == nbb ) pos1(3)= 0
#endif
                           idx1=floor(pos1)+1; idx2=idx1+1
                           dx1=idx1-pos1;      dx2=1-dx1
                           rho_grid(idx1(1),idx1(2),idx1(3))=rho_grid(idx1(1),idx1(2),idx1(3))+dx1(1)*dx1(2)*dx1(3)
                           rho_grid(idx2(1),idx1(2),idx1(3))=rho_grid(idx2(1),idx1(2),idx1(3))+dx2(1)*dx1(2)*dx1(3)
                           rho_grid(idx1(1),idx2(2),idx1(3))=rho_grid(idx1(1),idx2(2),idx1(3))+dx1(1)*dx2(2)*dx1(3)
                           rho_grid(idx1(1),idx1(2),idx2(3))=rho_grid(idx1(1),idx1(2),idx2(3))+dx1(1)*dx1(2)*dx2(3)
                           rho_grid(idx1(1),idx2(2),idx2(3))=rho_grid(idx1(1),idx2(2),idx2(3))+dx1(1)*dx2(2)*dx2(3)
                           rho_grid(idx2(1),idx1(2),idx2(3))=rho_grid(idx2(1),idx1(2),idx2(3))+dx2(1)*dx1(2)*dx2(3)
                           rho_grid(idx2(1),idx2(2),idx1(3))=rho_grid(idx2(1),idx2(2),idx1(3))+dx2(1)*dx2(2)*dx1(3)
                           rho_grid(idx2(1),idx2(2),idx2(3))=rho_grid(idx2(1),idx2(2),idx2(3))+dx2(1)*dx2(2)*dx2(3)
                        enddo
                        nlast=nlast+np
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
      sync all
      deallocate(xp)
      if (head) print*, 'sync buffer'
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
               rho_c(k,j,i)=rho_grid(k,j,i)
            enddo
         enddo
      enddo


      if (head) print*, 'check: min,max,sum of rho_grid = '
      ! if (head)
      print*, image,minval(rho_c),maxval(rho_c),sum(rho_c*1d0)
      rho8=sum(rho_c*1d0); sync all
      if (head) then
         do i=2,nn**3
            rho8=rho8+rho8[i]
         enddo
      endif; sync all
      rho8=rho8[1]; sync all
      print*,rho8/nw_global/nw_global/nw_global
      do i=1,nw
         rho_c(:,:,i)=rho_c(:,:,i)/(rho8/nw_global/nw_global/nw_global)-1
      enddo
      ! if (head)
      print*,image,'min',minval(rho_c),'max',maxval(rho_c),'mean',sum(rho_c*1d0)/nw/nw/nw; sync all
      if (nic < 3)then
         if (head) print*,'Write delta_c into',output_name('delta_c')
         open(11,file=output_name('delta_c'),status='replace',access='stream')
         write(11) rho_c
         close(11)
      endif
      sync all
#ifdef merge_projection
      allocate(rho_etc(nw,nw,nw))
      write(str_i,'(i6)') image
      write(str_z,'(f7.3)') z_checkpoint(cur_checkpoint)
      sync all
      if (head) then
         print*,'create merged density projection'
         allocate(proj_xy(nw_global,nw_global))
         proj_xy=0
         do iz=1,nn
            do iy=1,nn
               do ix=1,nn
                  idx=[ix,iy,iz]*nw ! upper index of each dim
                  rho_etc=rho_c(:,:,:)[ix,iy,iz]
                  proj_xy(idx(1)-nw+1:idx(1),idx(2)-nw+1:idx(2))&
                     =proj_xy(idx(1)-nw+1:idx(1),idx(2)-nw+1:idx(2))+sum(rho_etc*1d0,dim=3)/nw
               enddo
            enddo
         enddo
         proj_xy=proj_xy/nn
         open(11,file=output_name('proj_xy'),status='replace',access='stream')
         write(11) proj_xy
         close(11)
         deallocate(proj_xy)
      endif
      deallocate(rho_etc)
      sync all
#endif

      if (head) print*,'auto_power'
      call auto_power(xi,rho_c,npglobal,2)

      if (head) then
#ifdef power_matter
         print*,''
         print*,'matter_power'
         if( Mass_nu > 0 ) then
            write(str_z,'(f8.4)') z_checkpoint(cur_checkpoint)
            open(11,file=nupath//'tf/Tf_nu_'//trim(adjustl(str_z))//'.txt',status='old',access='stream'); read(11) tf_F;close(11)
            tf_F = (tf_F-(1-f_nu))/f_nu
            kh_lin = -1
            open(13,file=nupath//'k_values.txt',status='old')
            read(13, '(f15.12)') kh_lin(1:npbin)
            close(13)
            kh_lin_log = log(kh_lin)
            tf_F_log = log(tf_F**2)
            do i=1,nbin
               xi(6,i) = xi(5,i)*interp_tf_F(kh_lin_log,xi(2,i),tf_F_log) !Pk_nu
            enddo
            xi(7,:) = (f_nu*sqrt(xi(6,:))+(1-f_nu)*sqrt(xi(5,:)))**2  !Pk_tot
         endif
#endif
         print*,'Write cicpower',output_name('cicpower')
         open(15,file=output_name('cicpower'),status='replace',access='stream')
         write(15) xi(:,1:nbin)
         close(15)
      endif
      sync all

#ifdef density_nu
      if (head) print*,'create delta_nu'

      if (Mass_nu > 0) then
         rho1 = rho_c
         write(str_z,'(f8.4)') z_powerpoint(sim%cur_powerpoint)
         call cdm_2_nu(nupath//'tf/Tf_nu_'//trim(adjustl(str_z))//'.txt')
      else
         rho1 = 0
      endif


      if (head) print*, 'check: min,max,mean of rho_nu = '
      if (head) print*,'min',minval(rho_c),'max',maxval(rho_c),'mean',sum(rho_c*1d0)/nw/nw/nw; sync all
      if (head) print*,'Write delta_nu into',output_name('delta_nu')
      open(11,file=output_name('delta_nu'),status='replace',access='stream')
      write(11) rho1
      close(11); sync all

#endif

#ifdef density_matter
      if (head) print*,'create delta_matter'

      if (Mass_nu > 0) then
#ifndef density_nu
         rho1 = rho_c
         write(str_z,'(f8.4)') z_powerpoint(sim%cur_powerpoint)
         call cdm_2_nu(nupath//'tf/Tf_nu_'//trim(str_z)//'.txt')
#endif
         rho_1 = f_nu*rho1+(1-f_nu)*rho_c
      else
         rho1 = rho_c
      endif

      if (head) print*, 'check: min,max,mean of rho_matter = '
      if (head) print*,'min',minval(rho_c),'max',maxval(rho_c),'mean',sum(rho_c*1d0)/nw/nw/nw; sync all
      if (head) print*,'Write delta_matter into',output_name('delta_matter')
      open(11,file=output_name('delta_matter'),status='replace',access='stream')
      write(11) rho1
      close(11); sync all

#endif

#ifdef gadget
      ! read Gadget output
      allocate(xv(6,nplocal))
      open(11,file='../output/yuyu/0.000_gadget_xv.bin',access='stream')
      read(11) xv
      close(11)
      rho_grid=0;
      do ip=1,nplocal
         pos1=xv(1:3,ip) - 0.5
         idx1=floor(pos1)+1; idx2=idx1+1
         dx1=idx1-pos1;      dx2=1-dx1
         rho_grid(idx1(1),idx1(2),idx1(3))=rho_grid(idx1(1),idx1(2),idx1(3))+dx1(1)*dx1(2)*dx1(3)
         rho_grid(idx2(1),idx1(2),idx1(3))=rho_grid(idx2(1),idx1(2),idx1(3))+dx2(1)*dx1(2)*dx1(3)
         rho_grid(idx1(1),idx2(2),idx1(3))=rho_grid(idx1(1),idx2(2),idx1(3))+dx1(1)*dx2(2)*dx1(3)
         rho_grid(idx1(1),idx1(2),idx2(3))=rho_grid(idx1(1),idx1(2),idx2(3))+dx1(1)*dx1(2)*dx2(3)
         rho_grid(idx1(1),idx2(2),idx2(3))=rho_grid(idx1(1),idx2(2),idx2(3))+dx1(1)*dx2(2)*dx2(3)
         rho_grid(idx2(1),idx1(2),idx2(3))=rho_grid(idx2(1),idx1(2),idx2(3))+dx2(1)*dx1(2)*dx2(3)
         rho_grid(idx2(1),idx2(2),idx1(3))=rho_grid(idx2(1),idx2(2),idx1(3))+dx2(1)*dx2(2)*dx1(3)
         rho_grid(idx2(1),idx2(2),idx2(3))=rho_grid(idx2(1),idx2(2),idx2(3))+dx2(1)*dx2(2)*dx2(3)
      enddo
      sync all
      if (head) print*, 'sync buffer'
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
               rho_c(k,j,i)=rho_grid(k,j,i)
            enddo
         enddo
      enddo
      print*, 'check: min,max,sum of rho_gadget = '
      print*, minval(rho_c),maxval(rho_c),sum(rho_c*1d0)
      rho8=sum(rho_c*1d0); sync all
      if (head) then
         do i=2,nn**3
            rho8=rho8+rho8[i]
         enddo
      endif; sync all
      rho8=rho8[1]; sync all
      do i=1,nw
         rho_c(:,:,i)=rho_c(:,:,i)/(rho8/nw_global/nw_global/nw_global)-1
      enddo
      print*,'min',minval(rho_c),'max',maxval(rho_c),'mean',sum(rho_c*1d0)/nw/nw/nw; sync all
      if (head) print*,'Write delta_gadget into',output_name('delta_gadget')
      open(11,file=output_name('delta_gadget'),status='replace',access='stream')
      write(11) rho_c
      close(11); sync all
      open(11,file=output_name('delta_gadget_proj'),status='replace',access='stream')
      write(11) sum(rho_c(:,:,:),dim=3)/nw
      close(11); sync all
      write(str_i,'(i6)') image
      write(str_z,'(f7.3)') z_checkpoint(cur_checkpoint)

      print*,'auto_power'
      call auto_power(sim%a,xi,rho_c,npglobal,2)
      if (head) then
         open(15,file=output_name('gadgetpower'),status='replace',access='stream')
         write(15) xi(:,1:nbin)
         close(15)
      endif
      deallocate(xv)
      sync all
#endif
      if (head) print*,''
   enddo !z done

   deallocate(rho_c,rhoc)
   call destroy_penfft_plan
   sync all
   if (head) print*,'cicpower done'
   sync all
end
