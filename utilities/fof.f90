program CUBE_FoF
   ! use omp_lib
   use variables
   implicit none
   character(len = 4) str_refine
   integer,parameter:: fofcore = ncore
   integer(8),parameter:: fof_buffer = ceiling(20.0/box*nc*nn)
   integer(8),parameter:: nfof = (nc+2*fof_buffer)*ratio_cs
   real(8),parameter:: n_refine = nfof*1d0/(nc + fof_buffer*2 + b_link/ratio_cs)
   real(8),parameter:: L_b    = fof_buffer
   real(8),parameter:: L_bL   = L_b  + nc
   real(8),parameter:: L_bLb  = L_bL + fof_buffer
   real(8),parameter:: L_bLb2 = L_bLb/2
   real(8),parameter:: rp2=(b_link/ratio_cs)**2

   integer :: my_id          ! 当前 image 的 ID
   integer :: log_unit       ! 文件单元号
   character(len=512) :: log_filename ! 文件名字符串

   integer total_images,layer_image,image_now,i,j,k,l,nlayer,layer_core,cur_checkpoint,np,n1,n2,idx(3),nh[*],im,idx1(3),idx2(3),ft(6,nn*3),ftr(nn*3)
   integer(8) iq1,iq2,iq3,i_neighbor,jq(3),ijk_neighbor(3,3),neighbor_b(3,3)
   integer(8) nhalo_all,halo_images(nn**3)[*]
   integer(4) nlast,ip,jp,np_iso(fofcore),np_head(fofcore+1),np_max,np_neighbors(3,3,3)[nn,nn,*],np_need(3,3,3),max_nei[nn,nn,*],offset_nei(3,3,3),offset_team(fofcore)
   integer(4),allocatable :: np_halo_all(:),np_halo(:),np_halo_team(:,:), rhoc_local(:,:,:,:,:,:),offset_map(:,:,:,:,:,:)
   integer(4),allocatable :: hoc(:,:,:),ll(:),llgp(:),hcgp(:),ecgp(:),iph_halo_all(:),iph_halo(:)!,h2(:,:,:,:),l2(:,:)
   real rsq,pos1(3),dx1(3),dx2(3),shift_xv(3)
   real(8) rho8,dxv(6)
   real,allocatable :: xv(:,:),xv_mean(:,:),xv_mean_team(:,:,:),rho_grid(:,:,:)
   type(type_halo_catalog_header) halo_header
   type(type_halo_catalog_array),allocatable :: hcat(:)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   total_images = num_images()
   layer_image = nn**3/total_images
   if (total_images*layer_image /= nn**3) then
      if(head) print*,total_images,layer_image,nn**3
      stop 'total_images*layer_image /= nn*3'
   endif
   ft = 0

   nlayer = ceiling(fof_buffer*n_refine*2)+1
   do while (mod(nfof,nlayer) /= 0)
      nlayer = nlayer + 1
   enddo

   head=(this_image()==1)
   call omp_set_num_threads(fofcore)
   layer_core = min(fofcore,nfof/nlayer)

   if(head) call system("rm "//opath//"/fof/run_output_*.log")
   if(head) call system("mkdir -p "//opath//'/fof')

   write(str_refine,'(i4)')  nfof
   if(head) print*,nfof,nlayer,layer_core,total_images,layer_image

   if (head) then
      open(16,file='z_checkpoint.txt',status='old')
      do i=1,nmax_redshift-1
         read(16,end=71,fmt='(f8.4)') z_checkpoint(i)
      enddo
71    n_checkpoint=i-1
      close(16)
      if (n_checkpoint==0) stop 'z_checkpoint.txt empty'
   endif
   sync all
   n_checkpoint=n_checkpoint[1]
   z_checkpoint(:)=z_checkpoint(:)[1]
   sync all

   if (head) print*,'  initialize FoF cell neighbors',fof_buffer
   l=0
   do j=-1,1
      do i=-1,1
         l=l+1; ijk(:,l)=[i,j,-1]
      enddo
   enddo
   do i=-1,1
      l=l+1; ijk(:,l)=[i,-1,0]
   enddo
   l=l+1; ijk(:,l)=[-1,0,0]

   neighbor_b = reshape([int(0,kind=8),fof_buffer,nc  ,int(0,kind=8),nc,int(0,kind=8) ,nc-fof_buffer,nc,-nc ],[3,3])
   if (head) print*,'linking,L _b,L_bL,L_bLb',b_link/ratio_cs,L_b,L_bL,L_bLb

   do cur_checkpoint=n_checkpoint,n_checkpoint
      ! do cur_checkpoint=1,1
      sim%cur_checkpoint=cur_checkpoint
      if (head) print*,  ''
      if (head) print*,  ''
      if (head) print*,  'FoF at redshift ',z2str(z_checkpoint(cur_checkpoint))

      nhalo_all = 0
      do image_now = this_image(),nn**3,total_images


         write(log_filename, "(A,'/fof/run_output_', I4.4, '.log')") opath,image_now
         open(newunit=log_unit, file=trim(log_filename), status='replace', action='write')
         write(*,'(A, I0, A)') '  Fof at image',image_now,'->'//trim(log_filename)
         write(log_unit,'(A, I0)') '  Fof at image',image_now


         ! init image
         call geometry_images
         shift_xv  = [(icx-1)*nc*1d0-fof_buffer,(icy-1)*nc*1d0-fof_buffer,(icz-1)*nc*1d0-fof_buffer]

         ! count  neighbor particle
         call system_clock(ft(1,image_now),ftr(image_now))
         write(log_unit, *)' count  neighbor particle '
         allocate(rhoc_local(nt,nt,nt,nnt,nnt,nnt))
         np_neighbors = 0
         do iq1 = 1, 3
            do iq2 = 1, 3
               do iq3 = 1, 3
                  ijk_neighbor(:, 1) = neighbor_b(:, 4-iq1)
                  ijk_neighbor(:, 2) = neighbor_b(:, 4-iq2)
                  ijk_neighbor(:, 3) = neighbor_b(:, 4-iq3)

                  inx=modulo(icx+1-iq1,nn)+1
                  iny=modulo(icy+1-iq2,nn)
                  inz=modulo(icz+1-iq3,nn)
                  image = inz*nn**2+iny*nn+inx

                  open(11,file=output_name('np'),access='stream'); read(11) rhoc_local; close(11)
                  ! 循环遍历各个维度
                  do itz = floor((ijk_neighbor(1, 3))*1d0/nt)+1, floor((ijk_neighbor(2, 3)-1)*1d0/nt)+1
                     do ity = floor((ijk_neighbor(1, 2))*1d0/nt)+1, floor((ijk_neighbor(2, 2)-1)*1d0/nt)+1
                        do itx = floor((ijk_neighbor(1, 1))*1d0/nt)+1, floor((ijk_neighbor(2, 1)-1)*1d0/nt)+1
                           do k = merge(1_8, mod(ijk_neighbor(1, 3)+1, nt), itz*nt-nt >= ijk_neighbor(1, 3)),merge(nt, mod(ijk_neighbor(2, 3)+1, nt)-1, itz*nt <= ijk_neighbor(2, 3))
                              do j = merge(1_8, mod(ijk_neighbor(1, 2)+1, nt), ity*nt-nt >= ijk_neighbor(1, 2)),merge(nt, mod(ijk_neighbor(2, 2)+1, nt)-1, ity*nt <= ijk_neighbor(2, 2))
                                 do i = merge(1_8, mod(ijk_neighbor(1, 1)+1, nt), itx*nt-nt >= ijk_neighbor(1, 1)),merge(nt, mod(ijk_neighbor(2, 1)+1, nt)-1, itx*nt <= ijk_neighbor(2, 1))
                                    np_neighbors(iq1,iq2,iq3) = np_neighbors(iq1,iq2,iq3) + rhoc_local(i,j,k,itx,ity,itz)
                                 enddo
                              enddo
                           enddo
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
         deallocate(rhoc_local)
         np_max = sum(np_neighbors)





         !  initialize particles
         write(log_unit, *)' initialize particles', np_max

         allocate(xv(6,np_max))
         xv = 0
         nlast = 0
         do iq1 = 1, 3;  do iq2 = 1, 3; do iq3 = 1, 3
                  offset_nei(iq1,iq2,iq3) =  nlast
                  nlast = nlast + np_neighbors(iq1,iq2,iq3)
               enddo; enddo; enddo

         do iq1 = 1, 3
            do iq2 = 1, 3
               do iq3 = 1, 3
                  ijk_neighbor(:, 1) = neighbor_b(:, 4-iq1)
                  ijk_neighbor(:, 2) = neighbor_b(:, 4-iq2)
                  ijk_neighbor(:, 3) = neighbor_b(:, 4-iq3)

                  inx=modulo(icx+1-iq1,nn)+1
                  iny=modulo(icy+1-iq2,nn)
                  inz=modulo(icz+1-iq3,nn)
                  image = inz*nn**2+iny*nn+inx

                  iteam=omp_get_thread_num()+111
                  open(iteam,file=output_name('info'),access='stream'); read(iteam) sim; close(iteam)
                  sim%cur_checkpoint=cur_checkpoint
                  allocate(xp_new(3,sim%nplocal),vp_new(3,sim%nplocal),rhoc_local(nt,nt,nt,nnt,nnt,nnt),offset_map(nt,nt,nt,nnt,nnt,nnt))
                  open(iteam,file=output_name('xp'),access='stream'); read(iteam) xp_new; close(iteam)
                  open(iteam,file=output_name('vp'),access='stream'); read(iteam) vp_new; close(iteam)
                  open(iteam,file=output_name('np'),access='stream'); read(iteam) rhoc_local; close(iteam)

                  nlast = 0
                  do itz=1,nnt; do ity=1,nnt; do itx=1,nnt; do k=1,nt; do j=1,nt; do i=1,nt
                                    np = rhoc_local(i,j,k,itx,ity,itz)
                                    if (np < 0 ) then
                                       write(log_unit, '(7I4,I10)')  i,j,k,itx,ity,itz,image,np
                                       error stop 'particle index error'
                                    endif
                                    offset_map(i,j,k,itx,ity,itz) = nlast
                                    nlast = nlast + np
                                 enddo; enddo; enddo; enddo; enddo; enddo

                  jp = offset_nei(iq1,iq2,iq3)
                  n1 = jp+1
                  ! 循环遍历各个维度
                  do itz = floor((ijk_neighbor(1, 3))*1d0/nt)+1, floor((ijk_neighbor(2, 3)-1)*1d0/nt)+1
                     do ity = floor((ijk_neighbor(1, 2))*1d0/nt)+1, floor((ijk_neighbor(2, 2)-1)*1d0/nt)+1
                        do itx = floor((ijk_neighbor(1, 1))*1d0/nt)+1, floor((ijk_neighbor(2, 1)-1)*1d0/nt)+1
                           do k = merge(1_8, mod(ijk_neighbor(1, 3)+1, nt), itz*nt-nt >= ijk_neighbor(1, 3)),merge(nt, mod(ijk_neighbor(2, 3)+1, nt)-1, itz*nt <= ijk_neighbor(2, 3))
                              do j = merge(1_8, mod(ijk_neighbor(1, 2)+1, nt), ity*nt-nt >= ijk_neighbor(1, 2)),merge(nt, mod(ijk_neighbor(2, 2)+1, nt)-1, ity*nt <= ijk_neighbor(2, 2))
                                 do i = merge(1_8, mod(ijk_neighbor(1, 1)+1, nt), itx*nt-nt >= ijk_neighbor(1, 1)),merge(nt, mod(ijk_neighbor(2, 1)+1, nt)-1, itx*nt <= ijk_neighbor(2, 1))
                                    np=rhoc_local(i,j,k,itx,ity,itz)
                                    nlast = offset_map(i,j,k,itx,ity,itz)

                                    ! if (jp+np > np_max .or. nlast+np > sim%nplocal ) then
                                    !   write(log_unit, *)  iq1-2,iq2-2,iq3-2,image
                                    !   write(log_unit, *)  i,j,k,itx,ity,itz
                                    !   write(log_unit, *)  nlast,np,np_neighbors(iq1,iq2,iq3),jp
                                    !   write(log_unit, *)  'particle index error'
                                    !   error stop 'particle index error'
                                    ! endif
#ifdef ZIPX
                                    xv(1:3,jp+1:jp+np)=(int(xp_new(:,nlast+1:nlast+np)+ishift,izipx)+rshift)*x_resolution &
                                       +spread(nt*((/itx,ity,itz/)-1)+((/i,j,k/)-1)+[fof_buffer,fof_buffer,fof_buffer]+ijk_neighbor(3,:),dim=2,ncopies=np)
#else
                                    xv(1:3,jp+1:jp+np)=xp_new(:,nlast+1:nlast+np)+spread([fof_buffer,fof_buffer,fof_buffer]+ijk_neighbor(3,:),dim=2,ncopies=np)
#endif
                                    xv(4:6,jp+1:jp+np)=vp_new(:,nlast+1:nlast+np)
                                    jp = jp+np
                                 enddo
                              enddo
                           enddo
                        enddo
                     enddo
                  enddo
                  deallocate(xp_new,vp_new,rhoc_local,offset_map)
                  n2 = jp-n1+1
                  if (n2 /= np_neighbors(iq1,iq2,iq3) ) then
                     write(log_unit, *)  iq1-2,iq2-2,iq3-2,image
                     write(log_unit, *) np_neighbors(iq1,iq2,iq3),n2,jp,n1
                     write(log_unit, *) sum(xv(1,n1:jp))/n2,minval(xv(1,n1:jp)),maxval(xv(1,n1:jp)),ijk_neighbor(3,1)
                     write(log_unit, *) sum(xv(2,n1:jp))/n2,minval(xv(2,n1:jp)),maxval(xv(2,n1:jp)),ijk_neighbor(3,2)
                     write(log_unit, *) sum(xv(3,n1:jp))/n2,minval(xv(3,n1:jp)),maxval(xv(3,n1:jp)),ijk_neighbor(3,3)
                     write(log_unit, *)  'particle len error'
                     print*,image
                     stop 'particle len error'
                  endif
                  ! write(log_unit,'(5I4,I10)')  iq3+(iq2-1)*3+(iq1-1)*9,iq1-2,iq2-2,iq3-2,image,sim%nplocal
               enddo
            enddo
         enddo
         image = image_now


         call system_clock(ft(2,image_now),ftr(image_now))
         write(log_unit, '(A, F10.3, A)') '    real time =',real(ft(2,image_now)-ft(1,image_now))/ftr(image_now),'secs';write(log_unit, *) ''
         ! cycle



         ! create hoc ll
         write(log_unit, *)' create hoc ll'
         allocate(hoc(nfof,nfof,nfof),ll(np_max),llgp(np_max),hcgp(np_max),ecgp(np_max))
         hoc=0; ll=0
         do ip=1,np_max
            idx=floor(xv(1:3,ip)*n_refine)+1 ! index of the grid
            ll(ip)=hoc(idx(1),idx(2),idx(3)) ! linked list
            hoc(idx(1),idx(2),idx(3))=ip ! head of chain
            hcgp(ip)=ip ! initialize hcgp(ip)=ip for isolated particles
         enddo
         llgp=0; ecgp=0; ! initialize group link list


         call system_clock(ft(3,image_now),ftr(image_now))
         write(log_unit, '(A, F10.3, A)') '    real time =',real(ft(3,image_now)-ft(2,image_now))/ftr(image_now),'secs';write(log_unit, *) ''

         ! rho
         if (.false.) then
            write(log_unit, *)' create rho'
            allocate(rho_grid(0:nfof+1,0:nfof+1,0:nfof+1))
            rho_grid=0
            do l = 0,nlayer-1
               call system_clock(ft(5,image_now),ftr(image_now))
               !$omp parallel do num_threads(layer_core) schedule(dynamic,1) default(shared) &
               !$omp private(iq1,iq2,iq3,ip,pos1,idx1,idx2,dx1,dx2)
               do iq3=1+l,nfof,nlayer
                  do iq2=1,nfof
                     do iq1=1,nfof
                        ip=hoc(iq1,iq2,iq3)
                        do while (ip/=0)
                           pos1 = xv(1:3,ip)
                           if ( pos1(1) == L_bLb ) pos1(1)= 0
                           if ( pos1(2) == L_bLb ) pos1(2)= 0
                           if ( pos1(3) == L_bLb ) pos1(3)= 0
                           pos1 = pos1*n_refine
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
                           ip=ll(ip) ! find next particle in the chain
                        enddo
                     enddo
                  enddo
               enddo
               !$omp endparalleldo
            enddo
            rho8 = sum(rho_grid(1:nfof,1:nfof,1:nfof))/(nfof**3)
            rho_grid = rho_grid/rho8-1
            write(log_unit, *)'min',minval(rho_grid(1:nfof,1:nfof,1:nfof)),'max',maxval(rho_grid(1:nfof,1:nfof,1:nfof)),'mean',sum(rho_grid(1:nfof,1:nfof,1:nfof)*1d0)/nfof/nfof/nfof
            write(log_unit, *)'Write delta_fof into',output_name_halo('delta_fof_'//trim(adjustl(str_refine)))
            open(11,file=output_name_halo('delta_fof_'//trim(adjustl(str_refine))),status='replace',access='stream')
            write(11) rho_grid(1:nfof,1:nfof,1:nfof)
            close(11)
            deallocate(rho_grid)


            call system_clock(ft(4,image_now),ftr(image_now))
            write(log_unit, *) '    real time =',real(ft(4,image_now)-ft(3,image_now))/ftr(image_now),'secs';write(log_unit, *) ''
         else
            call system_clock(ft(4,image_now),ftr(image_now))
         endif


         ! loop over fof cells
         write(log_unit, *) 'Loop over fof cells'

         do l = 0,nlayer-1
            call system_clock(ft(2,image_now),ftr(image_now))
            !$omp parallel do num_threads(layer_core) schedule(dynamic,1) default(shared) &
            !$omp private(iq1,iq2,iq3,ip,jp,jq,rsq,i_neighbor)
            do iq3=1+l,nfof,nlayer
               do iq2=1,nfof
                  do iq1=1,nfof
                     ip=hoc(iq1,iq2,iq3)
                     do while (ip/=0)
                        jp=ll(ip)
                        do while (jp/=0)
                           rsq=sum((xv(1:3,ip)-xv(1:3,jp))**2)
                           if (rsq<=rp2) call merge_chain(ip,jp)
                           jp=ll(jp)
                        enddo
                        do i_neighbor=1,13
                           jq = [iq1,iq2,iq3]+ijk(:,i_neighbor)
                           if (maxval(jq)>nfof .or. minval(jq)<1) cycle
                           jp=hoc(jq(1),jq(2),jq(3))
                           do while (jp/=0)
                              rsq=sum((xv(1:3,ip)-xv(1:3,jp))**2)
                              if (rsq<=rp2) call merge_chain(ip,jp)
                              jp=ll(jp)
                           enddo
                        enddo
                        ip=ll(ip) ! find next particle in the chain
                     enddo
                  enddo
               enddo
            enddo
            !$omp endparalleldo
            call system_clock(ft(3,image_now),ftr(image_now))
            write(log_unit, '(A, I4, A, F10.3, A)') 'layer = ',l,'    real time =',real(ft(3,image_now)-ft(2,image_now))/ftr(image_now),'secs'
         enddo
         deallocate(hoc,ll,ecgp)

         call system_clock(ft(5,image_now),ftr(image_now))
         write(log_unit, '(A, F10.3, A)') '    real time =',real(ft(5,image_now)-ft(4,image_now))/ftr(image_now),'secs';write(log_unit, *) ''


         ! count halos
         nlast = np_max/np_halo_min/2
         allocate(xv_mean_team(6, nlast, fofcore),np_halo_team(nlast, fofcore))
         write(log_unit, *)  'Count halos'
         np_iso = 0; np_head = 0; np_halo_team=0; xv_mean_team=0
         !$omp parallel do private(i, iteam, ip, jp, np, dxv) default(shared)
         do i=1,np_max
            iteam = omp_get_thread_num()+1
            if (hcgp(i)==i) then
               np_iso(iteam)=np_iso(iteam)+1
            elseif (hcgp(i)/=0) then
               ip = hcgp(i); jp = ip; np = 0; dxv = 0
               do while (jp /= 0)
                  np = np + 1
                  dxv(1:3) = dxv(1:3) + modulo(xv(1:3,jp) - xv(1:3,ip) + L_bLb2, L_bLb*1d0)-L_bLb2
                  dxv(4:6) = dxv(4:6) + xv(4:6,jp)
                  jp = llgp(jp)
               enddo
               dxv(1:3) = modulo(dxv(1:3)/np + xv(1:3,ip), L_bLb*1d0)
               dxv(4:6) = dxv(4:6)/np
               if (np > np_halo_min .and. minval(dxv(1:3)) > L_b .and. maxval(dxv(1:3)) < L_bL) then
                  np_head(iteam)=np_head(iteam)+1
                  np_halo_team(np_head(iteam),iteam) = np
                  xv_mean_team(1:3,np_head(iteam),iteam) =  (dxv(1:3) + shift_xv)*box/nn/nc
                  xv_mean_team(4:6,np_head(iteam),iteam) =  dxv(4:6)
               endif
            endif
         enddo
         !$omp end parallel do
         deallocate(xv,hcgp,llgp)

         offset_team = 0; nlast = 0
         do i = 1,fofcore
            np_iso(1)  = np_iso(1) + np_iso(i)
            offset_team(i) = nlast
            nlast = nlast + np_head(i)
         enddo

         allocate(np_halo(nlast))
         do i = 1,fofcore
            np_halo(offset_team(i)+1:offset_team(i)+np_head(i)) = np_halo_team(1:np_head(i),i)
         enddo
         deallocate(np_halo_team)

         allocate(xv_mean(6,nlast))
         do i = 1,fofcore
            xv_mean(:, offset_team(i)+1:offset_team(i)+np_head(i)) = xv_mean_team(:, 1:np_head(i),i)
         enddo
         deallocate(xv_mean_team)


         np_head(1) = nlast
         np_head(2) = sum(np_halo)
         write(log_unit,'(A, I10, A, I10, I10)') '        np_iso: ',np_iso(1),'  np_head: ',np_head(1),np_head(2)
         call system_clock(ft(6,image_now),ftr(image_now))
         write(log_unit, '(A, F10.3, A)') '    real time =',real(ft(6,image_now)-ft(5,image_now))/ftr(image_now),'secs';write(log_unit, *) ''


         nh=np_head(1)
         halo_header%nhalo_tot=0
         halo_header%nhalo=nh; halo_header%ninfo=ninfo; halo_header%linking_parameter=b_link

         allocate(hcat(nh))
         hcat%hmass = np_halo
         do i=1,6
            hcat%xv(i)=xv_mean(i,1:nh)
         enddo
         write(log_unit, *)'halo_header',halo_header,maxval(hcat%hmass)


         write(log_unit, *) output_name_halo(trim('halo_'//trim(adjustl(str_refine))))
         open(21,file=output_name_halo(trim('halo_'//trim(adjustl(str_refine)))),status='replace',access='stream')
         write(21) halo_header
         write(21) hcat
         close(21)

         ! write(log_unit, *) output_name_halo(trim('halo_mass_'//trim(adjustl(str_refine))))
         ! open(21,file=output_name_halo(trim('halo_mass_'//trim(adjustl(str_refine)))),status='replace',access='stream')
         ! write(21) nh
         ! write(21) np_halo(1:nh)
         ! close(21)
         deallocate(np_halo,xv_mean)

         deallocate(hcat)
         halo_images(image_now) = np_head(1)

         write(log_unit, *) np_head(1),'all time =',real(ft(6,image_now)-ft(1,image_now))/ftr(image_now),'secs'

         write(*,'(A, I0, A, A, I0, A, I10, I15, A, F10.3, A)') '  Fof at image',image_now,'->'//trim(log_filename),' |  np_max :', np_max,' |  find halo :',np_head(1),np_head(2),' |  used time =',real(ft(6,image_now)-ft(1,image_now))/ftr(image_now),'secs'

         close(log_unit)
      enddo

      ! write halo mass
      sync all
      if (head) then
         write(*, *) ''
         write(*, *) ''
         write(*, *) ''
         do image=2,total_images
            halo_images = halo_images + halo_images(:)[image]
         enddo
         nhalo_all = sum(halo_images)
         write(*, *) nhalo_all
         do image=1,nn**3
            write(*, *) image,halo_images(image)
         enddo
      endif
      sync all
      if (head) then
         image = 0
         halo_header%nhalo_tot=nhalo_all
         write(*, *) output_name_halo(trim('halo_all_'//trim(adjustl(str_refine))))
         open(111,file=output_name_halo(trim('halo_all_'//trim(adjustl(str_refine)))),status='replace',access='stream')
         write(111) halo_header
         do image=1,nn**3
            open(21,file=output_name_halo(trim('halo_'//trim(adjustl(str_refine)))),status='old',access='stream')
            read(21) halo_header; allocate(hcat(halo_header%nhalo)); read(21) hcat; close(21)
            write(*, *) output_name_halo(trim('halo_'//trim(adjustl(str_refine)))),sum(hcat%hmass),halo_header%nhalo
            print*,hcat(1)
            print*,'================'
            print*,hcat(2)
            print*,'================'
            print*,hcat(10)
            write(111) hcat
            deallocate(hcat)
         enddo
      endif




   enddo ! cur_checkpoint



contains

   subroutine merge_chain(ii,jj)
      ! llgp is a linked list: llgp(ip) means ip->llgp
      ! ip1->ip2->...->ipn
      ! ip1 is the head of chain (hoc)
      ! ipn is the end of chain (eoc)
      integer(4) ii,jj,ihead,jhead,iend,jend,ipart
      jend=merge(jj,ecgp(jj),ecgp(jj)==0)
      iend=merge(ii,ecgp(ii),ecgp(ii)==0)
      if (iend==jend) return ! same chain
      ihead=max(hcgp(ii),hcgp(iend))
      jhead=max(hcgp(jj),hcgp(jend))
      ipart=jhead ! change eoc of the chain-j
      do while (ipart/=0)
         ecgp(ipart)=iend ! set chain-j's eoc to be iend
         ipart=llgp(ipart) ! next particle
      enddo
      llgp(jend)=ihead ! link j group to i group
      ecgp(ii)=iend
      hcgp(iend)=jhead ! change hoc
      hcgp(jend)=0 ! set jend as a member
   endsubroutine

   subroutine geometry_images
      write(log_unit, *) 'geometry_images'
      rank=image_now-1               ! MPI_rank
      icz=rank/(nn**2)+1             ! image_z
      icy=(rank-nn**2*(icz-1))/nn+1  ! image_y
      icx=mod(rank,nn)+1             ! image_x
   endsubroutine

end
