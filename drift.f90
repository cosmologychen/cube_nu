subroutine drift
   use variables
   implicit none


   integer i,j,k,l,np,ilayer,nlayer
   integer(8) idx,np_prev,ig(3),npcheck,ip,nzero,idx_ex_r(1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb)
   integer(4),allocatable :: rhoce(:,:,:),rholocal(:,:,:)
   integer(8),dimension(nt,nt) :: pp_l,pp_r,ppe_l,ppe_r
   real(8) xq(ndim),deltax(ndim)!,vreal(ndim)
   real(8) sum_c,sum_r,sum_s
   !real,allocatable :: vfield_new(:,:,:,:)

   if (head) print*,'drift'
   call tic(2)
   call system_clock(t1,t_rate)
   !allocate(vfield_new(3,1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb))
   allocate(rhoce(1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb))
   allocate(rholocal(1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb))
   dt_mid=(dt_old+dt)/2
   if (head) print*,'  dt_mid =',dt_mid
   overhead_tile=0
   np_prev=0

   nlayer=2*ceiling(dt_mid*sim%vz_max/ratio_cs)+2
   if (head) print*,'  nlayer =',nlayer

   do itz=1,nnt ! loop over tile
      do ity=1,nnt
         do itx=1,nnt
            sync all
            !if (head) print*,'itx,ity,itz'; sync all
            rhoce=0
            rholocal=0
            ! for empty coarse cells, use previous-step vfield
            !vfield_new=0
            !vfield_new(:,1-ncb:nt+ncb,1-ncb:nt+ncb,1-ncb:nt+ncb)=vfield(:,:,:,:,itx,ity,itz)*weight_v
            xp_new=0
            vp_new=0
#   ifdef PID
            pid_new=0
#   endif
            sync all
            !if (head) print*,'    density loop'; sync all
            do ilayer=0,nlayer-1
               !$omp paralleldo default(shared) schedule(dynamic)&
               !$omp& private(k,j,i,np,nzero,l,ip,xq,deltax,ig)
               do k=1-ncb+ilayer,nt+ncb,nlayer
                  !do k=1-ncb,nt+ncb
                  do j=1-ncb,nt+ncb
                     do i=1-ncb,nt+ncb
                        np=rhoc(i,j,k,itx,ity,itz)
                        nzero=idx_b_r(j,k,itx,ity,itz)-sum(rhoc(i:,j,k,itx,ity,itz))
                        do l=1,np
                           ip=nzero+l
                           deltax=(dt_mid*vp(:,ip))/ratio_cs
#ifdef ZIPX
                           xq=((/i,j,k/)-1d0) + (int(xp(:,ip)+ishift,izipx)+rshift)*x_resolution
#else
                           xq=xp(:,ip)/ratio_cs - ([itx,ity,itz]-1)*nt 
#endif
                           ig=ceiling(xq + deltax)
                           if (minval(ig) < 1-2*ncb .or. maxval(ig) > nt+2*ncb) then
                              print*, '  error: particle out of range',1-2*ncb,nt+2*ncb
                              print*, '  ',ig
                              print*, '  ',xp(:,ip)/ratio_cs
                              print*, '  ',deltax
                              print*, '  ',([itx,ity,itz]-1)*nt*ratio_cs
                              print*, '  on',this_image(), itx,ity,itz
                              print*,ip
                              error stop
                           endif
                           rhoce(ig(1),ig(2),ig(3))=rhoce(ig(1),ig(2),ig(3))+1
                        enddo
                     enddo
                  enddo
               enddo ! k
               !$omp endparalleldo
            enddo ! ilayer
            sync all
            

            call spine_tile(rhoce,idx_ex_r,pp_l,pp_r,ppe_l,ppe_r)
            overhead_tile=max(overhead_tile,idx_ex_r(nt+2*ncb,nt+2*ncb)/real(np_tile_max))

            if (idx_ex_r(nt+2*ncb,nt+2*ncb)>np_tile_max) then
               print*, '  error: too many particles in this tile+buffer'
               print*, '  ',idx_ex_r(nt+2*ncb,nt+2*ncb),'>',np_tile_max
               print*, '  on',this_image(), itx,ity,itz
               print*, '  please set tile_buffer larger'
               error stop
            endif
            ! create x_new v_new for this local tile
            !if (head) print*,'    xvnew loop'
            do ilayer=0,nlayer-1
               !$omp paralleldo default(shared) schedule(dynamic)&
               !$omp& private(k,j,i,np,nzero,l,ip,xq,deltax,ig,idx)
               do k=1-ncb+ilayer,nt+ncb,nlayer
                  !do k=1-ncb,nt+ncb
                  do j=1-ncb,nt+ncb
                     do i=1-ncb,nt+ncb
                        np=rhoc(i,j,k,itx,ity,itz)
                        nzero=idx_b_r(j,k,itx,ity,itz)-sum(rhoc(i:,j,k,itx,ity,itz))
                        do l=1,np
                           ip=nzero+l
                           deltax=(dt_mid*vp(:,ip))/ratio_cs
#ifdef ZIPX
                           xq=((/i,j,k/)-1d0) + (int(xp(:,ip)+ishift,izipx)+rshift)*x_resolution
#else
                           xq=xp(:,ip)/ratio_cs - ([itx,ity,itz]-1)*nt 
#endif
                           ig=ceiling(xq + deltax)
                           rholocal(ig(1),ig(2),ig(3))=rholocal(ig(1),ig(2),ig(3))+1
                           idx=idx_ex_r(ig(2),ig(3))-sum(rhoce(ig(1):,ig(2),ig(3)))+rholocal(ig(1),ig(2),ig(3))
#ifdef ZIPX       
                           xp_new(:,idx)=xp(:,ip)+nint(dt_mid*vp(:,ip)/(x_resolution*ratio_cs))
#else
                           xp_new(:,idx)=xp(:,ip)+deltax*ratio_cs
#endif
                           vp_new(:,idx)=vp(:,ip)
#         ifdef PID
                           pid_new(idx)=pid(ip)
                           if (pid(ip)==0) then
                              print*,'found pid=0'
                              print*,itx,ity,itz,i,j,k,np
                              print*,ip,xq
                              error stop
                           endif
#         endif
                        enddo
                     enddo
                  enddo
               enddo ! k
               !$omp endparalleldo
            enddo ! ilayer
            sync all

            if (sum(abs(rholocal-rhoce))/=0) then
               print*,'  warning: density different'
               print*,maxloc(abs(rholocal-rhoce))
               error stop
            endif
            if (minval(pid_new(:sum(rholocal)))==0) then
               print*,'  found 0 in pid_new'
               error stop
            endif


            ! delete particles
            sync all
            !if (head) print*,'    delete_particle loop'; sync all
            !$omp paralleldo default(shared) schedule(dynamic)&
            !$omp& private(k,j)
            do k=1,nt
               do j=1,nt
                  xp(:,np_prev+pp_l(j,k):np_prev+pp_r(j,k))=xp_new(:,ppe_l(j,k):ppe_r(j,k))
                  vp(:,np_prev+pp_l(j,k):np_prev+pp_r(j,k))=vp_new(:,ppe_l(j,k):ppe_r(j,k))
#     ifdef PID
                  pid(np_prev+pp_l(j,k):np_prev+pp_r(j,k))=pid_new(ppe_l(j,k):ppe_r(j,k))
                  if (minval(pid_new(ppe_l(j,k):ppe_r(j,k)))==0) then
                     print*,'found 0 in delete'
                     print*,itx,ity,itz
                     print*,k,j
                     print*,ppe_l(j,k),ppe_r(j,k)
                     print*,pid_new(ppe_l(j,k):ppe_r(j,k))
                     error stop
                  endif
#     endif
               enddo
            enddo
            !$omp endparalleldo
            np_prev=np_prev+pp_r(nt,nt)
            ! update rhoc
            rhoc(:,:,:,itx,ity,itz)=0
            rhoc(1:nt,1:nt,1:nt,itx,ity,itz)=rhoce(1:nt,1:nt,1:nt)
            !vfield(:,:,:,:,itx,ity,itz)=0
            !vfield(:,1:nt,1:nt,1:nt,itx,ity,itz)=vfield_new(:,1:nt,1:nt,1:nt)
         enddo
      enddo
   enddo ! end looping over tiles
   sim%nplocal=np_prev
   xp(:,sim%nplocal+1:)=0
   vp(:,sim%nplocal+1:)=0
# ifdef PID
   pid(sim%nplocal+1:)=0
   !print*, '  check PID range:',minval(pid(:sim%nplocal)),maxval(pid(:sim%nplocal))
   if (minval(pid(:sim%nplocal))<1) then
      print*, '  pid are not all positive'
      print*, '  number of 0s =',count(pid(:sim%nplocal)==0)
      print*, minloc(pid(:sim%nplocal))
      error stop
   endif
# endif
   !sim%nplocal=sum(rhoc(1:nt,1:nt,1:nt,:,:,:))
   !print*, iright,sum(rhoc(1:nt,1:nt,1:nt,:,:,:))
   !error stop
   sync all

   ! calculate std of the velocity field
   !cum=cumsum6(rhoc)
   call spine_image(rhoc,idx_b_l,idx_b_r,ppl0,pplr,pprl,ppr0,ppl,ppr)
   
   if (head) then
      do i=1,nn**3
         overhead_tile=max(overhead_tile,overhead_tile[i])
      enddo
   endif
   sync all
   overhead_tile=overhead_tile[1]
   sync all

   if (head) then
      print*, '  tile overhead',overhead_tile*100,'% full'
      print*, '  consumed tile_buffer =',overhead_tile*tile_buffer,'/',tile_buffer
      tcat(52,istep)=overhead_tile*tile_buffer
      print*, '  clean buffer of rhoc'
   endif

   if (head) then
      npcheck=0
      do i=1,nn**3
         npcheck=npcheck+sim[i]%nplocal
      enddo
      print*, '  npcheck,npglobal=', npcheck,sim%npglobal
      call system_clock(t2,t_rate)
      print*, '  elapsed time =',real(t2-t1)/t_rate,'secs'
      print*, ''
   endif
   deallocate(rhoce,rholocal)
   call toc(2)
   sync all
endsubroutine
