module pencil_fft_global
   use omp_lib
   use parameters
   implicit none
   save

   integer(4),parameter :: NULL_global=0
   integer(8) planx_global,plany_global,planz_global,iplanx_global,iplany_global,iplanz_global
   real,allocatable ::         rho1_global(:,:,:)[:,:,:]
   complex,allocatable ::      rho1k_global(:,:,:)
   complex,allocatable ::      cxyz_global(:,:,:)
   complex,allocatable ::      cyyxz_global(:,:,:,:)[:,:,:]
   complex,allocatable ::      czzzxy_global(:,:,:,:,:)[:,:,:]


contains

   subroutine pencil_fft_forward_global
      use omp_lib
      implicit none
      save

      sync all
      call system_clock(ttt1,t_rate)
      call system_clock(tt1,t_rate)
      call c2x_global
      sync all; call system_clock(tt2,t_rate); if (head) print*,'     c2x_global',real(tt2-tt1)/t_rate
      !call flush(6)
      sync all

      call system_clock(tt1,t_rate)
      call sfftw_execute(planx_global)
      sync all; call system_clock(tt2,t_rate); if (head) print*,'     xtran_global',real(tt2-tt1)/t_rate
      !call flush(6)
      sync all

      call system_clock(tt1,t_rate)
      call x2y_global
      sync all; call system_clock(tt2,t_rate); if (head) print*,'     x2y_global',real(tt2-tt1)/t_rate
      !call flush(6)
      sync all

      call system_clock(tt1,t_rate)
      call sfftw_execute(plany_global)
      sync all; call system_clock(tt2,t_rate); if (head) print*,'     ytran_global',real(tt2-tt1)/t_rate
      !call flush(6)
      sync all

      call system_clock(tt1,t_rate)
      call y2z_global
      sync all; call system_clock(tt2,t_rate); if (head) print*,'     y2z_global',real(tt2-tt1)/t_rate
      !call flush(6)
      sync all

      call system_clock(tt1,t_rate)
      call sfftw_execute(planz_global)
      sync all; call system_clock(tt2,t_rate); if (head) print*,'     ztran_global',real(tt2-tt1)/t_rate
      !call flush(6)
      sync all

      call system_clock(tt1,t_rate)
      call z2y_global
      sync all; call system_clock(tt2,t_rate); if (head) print*,'     z2y_global',real(tt2-tt1)/t_rate
      !call flush(6)
      sync all

      call system_clock(tt1,t_rate)
      call y2x_global
      sync all; call system_clock(tt2,t_rate); if (head) print*,'     y2x_global',real(tt2-tt1)/t_rate
      !call flush(6)
      sync all
      
      cxyz_global=cxyz_global/nfg_global/nfg_global/nfg_global
      sync all; call system_clock(ttt2,t_rate); if (head) print*,'     pen_forward_global time',real(ttt2-ttt1)/t_rate
      !call flush(6)
      sync all
      ! stop
   endsubroutine

   ! subroutine pencil_fft_backward
   !   use omp_lib
   !   implicit none
   !   save
   !   sync all
   !   call x2y
   !   call y2z
   !   call sfftw_execute(iplanz)
   !   call z2y
   !   call sfftw_execute(iplany)
   !   call y2x
   !   call sfftw_execute(iplanx)
   !   call x2c
   !   rho1=rho1/nfg_global/nfg_global/nfg_global
   !   sync all
   ! endsubroutine

   subroutine c2x_global
      use omp_lib
      implicit none
      save
      integer(8) i1,islab
      complex,allocatable :: ctransfer(:,:,:)[:,:,:]

     ! print*,image
      !call flush(6)
      sync all

      rho1k_global=cmplx(rho1_global(1::2,:,:),rho1_global(2::2,:,:))
      allocate(ctransfer(nfg/2,nfg,nn)[nn,nn,*])
      sync all
      ! print*,rho1k_global(10,11,12)
      ! print*,rho1_global(19:20,11,12)
      ! stop
      do islab=1,(nfg/nn) ! loop over cells in z, extract slabs
         ctransfer=rho1k_global(:,:,islab::(nfg/nn)) ! nn slabs of rho1k copied to ctransfer1
         sync all
         do i1=1,nn ! loop over parts in x, get slabs from each y node
            cxyz_global(nfg*(i1-1)/2+1:nfg*i1/2,:,islab)=ctransfer(:,:,m2)[i1,m1,m3]
         enddo
         sync all
      enddo
      deallocate(ctransfer)
      !print*,'a',image
      !call flush(6)
      sync all
   endsubroutine

   subroutine x2y_global
      use omp_lib
      implicit none
      save
      integer(8) i1,islab
      complex,allocatable :: ctransfer(:,:,:)[:,:,:]

     ! print*,image
      !call flush(6)
      sync all

      allocate(ctransfer(nfg,nfg/2+1,nn)[nn,nn,*])
      sync all
      do islab=1,(nfg/nn) ! loop over z
         do i1=1,nn ! loop over squares in x direction
            ctransfer(:,:,i1)=transpose(cxyz_global(nfg/2*(i1-1)+1:nfg/2*i1+1,:,islab))
         enddo
         sync all
         do i1=1,nn
            cyyxz_global(:,i1,:,islab)=ctransfer(:,:,m1)[i1,m2,m3]
         enddo
         sync all
      enddo
      deallocate(ctransfer)
      sync all
      !print*,'a',image
      !call flush(6)
      sync all
   endsubroutine

   subroutine y2z_global
      use omp_lib
      implicit none
      save
      integer(8) i1,i2,islab
      complex,allocatable :: ctransfer(:,:,:,:)[:,:,:]

     ! print*,image
      !call flush(6)
      sync all

      allocate(ctransfer((nfg/nn),(nfg/nn),nn,nn)[nn,nn,*])
      sync all
      do islab=1,nfg/2+1 ! loop over slices in x direction
         do i2=1,nn
            do i1=1,nn
               ctransfer(:,:,i1,i2)=transpose(cyyxz_global((nfg/nn)*(i1-1)+1:(nfg/nn)*i1,i2,islab,:))
            enddo
         enddo
         sync all
         do i2=1,nn
            do i1=1,nn
               czzzxy_global(:,i1,i2,islab,:)=ctransfer(:,:,m2,m3)[m1,i1,i2]
            enddo
         enddo
         sync all
      enddo
      deallocate(ctransfer)
      sync all
      !print*,'a',image
      !call flush(6)
      sync all
   endsubroutine

   subroutine z2y_global
      implicit none
      save
      integer(8) i1,i2,islab
      complex,allocatable :: ctransfer(:,:,:,:)[:,:,:]

     ! print*,image
      !call flush(6)
      sync all

      allocate(ctransfer((nfg/nn),(nfg/nn),nn,nn)[nn,nn,*])
      sync all
      do islab=1,nfg/2+1 ! loop over slices in x direction
         do i2=1,nn
            do i1=1,nn
               ctransfer(:,:,i1,i2)=transpose(czzzxy_global(:,i1,i2,islab,:))
            enddo
         enddo
         sync all
         do i2=1,nn
            do i1=1,nn
               !cyyyxz(:,i1,i2,islab,:)=ctransfer(:,:,m2,m3)[m1,i1,i2]
               cyyxz_global((nfg/nn)*(i1-1)+1:(nfg/nn)*i1,i2,islab,:)=ctransfer(:,:,m2,m3)[m1,i1,i2]
            enddo
         enddo
         sync all
      enddo
      deallocate(ctransfer)
      sync all
      !print*,'a',image
      !call flush(6)
      sync all
   endsubroutine

   subroutine y2x_global
      implicit none
      save
      integer(8) i1,islab
      complex,allocatable :: ctransfer(:,:,:)[:,:,:]

     ! print*,image
      !call flush(6)
      sync all

      allocate(ctransfer(nfg/2+1,nfg,nn)[nn,nn,*])
      sync all
      do islab=1,(nfg/nn) ! loop over z
         do i1=1,nn ! loop over squares in x direction
            ctransfer(:,:,i1)=transpose(cyyxz_global(:,i1,:,islab))
         enddo
         sync all
         do i1=1,nn
            cxyz_global(nfg/2*(i1-1)+1:nfg/2*i1+1,:,islab)=ctransfer(:,:,m1)[i1,m2,m3]
         enddo
         sync all
      enddo
      deallocate(ctransfer)
      sync all
      !print*,'a',image
      !call flush(6)
      sync all
   endsubroutine

   ! subroutine x2c
   !   implicit none
   !   save
   !   integer(8) i1,islab
   !   real,allocatable :: rtransfer(:,:,:)[:,:,:], rxyz_send(:,:,:)[:,:,:]
   !   if (.true.) then
   !   allocate(rtransfer(nfg,nfg,nn)[nn,nn,*])
   !   do islab=1,(nfg/nn)
   !     do i1=1,nn
   !       rtransfer(:,:,i1)=rxyz(nfg*(i1-1)+1:nfg*i1,:,islab)
   !     enddo
   !     sync all
   !     do i1=1,nn
   !       rho1(:,:,islab+(i1-1)*(nfg/nn))=rtransfer(:,:,m1)[m2,i1,m3]
   !     enddo
   !     sync all
   !   enddo
   !   deallocate(rtransfer)
   !   endif

   !   if (.false.) then
   !     allocate(rxyz_send(nfg_global+2  ,nfg,(nfg/nn))[nn,nn,*])
   !     rxyz_send=rxyz; sync all
   !     do i1=1,nn
   !       rho1(:,:,(i1-1)*(nfg/nn)+1:i1*(nfg/nn))=rxyz_send(nfg*(m1-1)+1:nfg*m1,:,:)[m2,i1,m3]
   !     enddo
   !     deallocate(rxyz_send)
   !   endif
   !   sync all
   ! endsubroutine

   subroutine create_penfft_plan_global
      implicit none
      save
      include 'fftw3.f'
      sync all
      allocate(rho1k_global(nfg/2,nfg,nfg),&
               cxyz_global(nfg_global/2+1,nfg,(nfg/nn)),&
               cyyxz_global(nfg,nn,nfg/2+1,(nfg/nn))[nn,nn,*],&
               czzzxy_global((nfg/nn),nn,nn,nfg/2+1,(nfg/nn))[nn,nn,*])

      call sfftw_plan_many_dft_r2c(planx_global,1,nfg_global,nfg*(nfg/nn),cxyz_global,NULL_global,1,nfg_global+2,cxyz_global,NULL_global,1,nfg_global/2+1,FFTW_MEASURE)
      call sfftw_plan_many_dft(plany_global,1,nfg_global,(nfg/2+1)*(nfg/nn),cyyxz_global,NULL_global,1,nfg_global,cyyxz_global,NULL_global,1,nfg_global,FFTW_FORWARD,FFTW_MEASURE)
      call sfftw_plan_many_dft(planz_global,1,nfg_global,(nfg/2+1)*(nfg/nn),czzzxy_global,NULL_global,1,nfg_global,czzzxy_global,NULL_global,1,nfg_global,FFTW_FORWARD,FFTW_MEASURE)
      sync all
   endsubroutine

   subroutine destroy_penfft_plan_global
      implicit none
      save
      include 'fftw3.f'
      sync all
      deallocate(rho1_global,rho1k_global,cxyz_global,cyyxz_global,czzzxy_global)

      call sfftw_destroy_plan(planx_global)
      call sfftw_destroy_plan(plany_global)
      call sfftw_destroy_plan(planz_global)
      sync all
   endsubroutine
endmodule
