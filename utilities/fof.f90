module buffered_output
   ! 完整复制自 fof_new/clean/buffered_output.f90；双缓冲和原格式输出不变。
   use iso_fortran_env, only: int8, int32, int64, real32
   use omp_lib
   implicit none

   private
   public :: write_members_buffered

   integer(int8), parameter :: mold_4(4) = 0_int8
   integer(int8), parameter :: mold_8(8) = 0_int8
   integer(int8), parameter :: mold_24(24) = 0_int8

   ! 段元数据：记录切片归属的 Halo 及其在双缓冲区内的相对偏移
   type halo_segment
      integer(int32) :: local_offset
      integer(int32) :: halo_idx
      integer(int32) :: member_start
      integer(int32) :: member_count
      logical :: write_header
   end type halo_segment

contains

   ! 流式并行双缓冲写入 Halo 成员与相空间数据
   ! 二进制布局：每个 Halo 写入 28 字节头部 (np_halo 4字节 + xv_mean 24字节)，
   ! 随后连续写入 np_halo 个 32 字节成员记录 (pid 8字节 + xv 24字节)。
   subroutine write_members_buffered(unit, nh, np_halo, xv_mean, halo_heads, &
                                    group_ptr, group_members, pid, xv, threads, &
                                    buffer_bytes)
      integer, intent(in) :: unit
      integer(int32), intent(in) :: nh
      integer(int32), intent(in) :: np_halo(nh)
      real(real32), intent(in) :: xv_mean(6, nh)
      integer(int32), intent(in) :: halo_heads(nh)
      integer(int32), intent(in) :: group_ptr(:)
      integer(int32), intent(in) :: group_members(:)
      integer(int64), intent(in) :: pid(:)
      real(real32), intent(in) :: xv(:, :)
      integer, intent(in) :: threads
      integer(int64), intent(in), optional :: buffer_bytes

      integer(int64) :: buf_cap
      integer :: n_batches, n_segs, b, s
      integer :: slot_dep(2), write_order
      integer :: ntasks

      type(halo_segment), allocatable :: segs(:)
      integer(int32), allocatable :: batch_bytes(:)
      integer(int32), allocatable :: batch_seg_start(:), batch_seg_end(:)
      integer(int8), allocatable :: buf(:, :)

      if (present(buffer_bytes)) then
         buf_cap = buffer_bytes
      else
         buf_cap = 1073741824_int64
      end if

      ! 第一阶段：规划批次与段边界
      b = 1
      n_segs = 0
      call partition_stream(1, buf_cap, nh, np_halo, b, n_segs, segs, &
                           batch_bytes, batch_seg_start, batch_seg_end)
      n_batches = b

      ! 分配元数据与双缓冲
      allocate(segs(n_segs))
      allocate(batch_bytes(n_batches), batch_seg_start(n_batches), batch_seg_end(n_batches))

      ! 第二阶段：记录段详细切分
      call partition_stream(2, buf_cap, nh, np_halo, b, n_segs, segs, &
                           batch_bytes, batch_seg_start, batch_seg_end)

      allocate(buf(buf_cap, 2))

      slot_dep = 0
      write_order = 0
      ntasks = max(1, threads * 4)

      !$omp parallel num_threads(threads) default(shared) private(b, s)
      !$omp single
      do b = 1, n_batches
         s = mod(b - 1, 2) + 1

         ! 1. 填充任务：依赖同 slot 上一轮写入释放
         !$omp task depend(inout: slot_dep(s)) firstprivate(b, s, ntasks)
         block
            integer :: iseg, n_items
            integer(int32) :: ih, ip, igroup, gpos, member
            integer(int32) :: pos, np, member_start

            n_items = max(1, min(ntasks, batch_seg_end(b) - batch_seg_start(b) + 1))

            !$omp taskloop num_tasks(n_items) default(shared) &
            !$omp private(iseg, ih, ip, igroup, gpos, member, pos, np, member_start)
            do iseg = batch_seg_start(b), batch_seg_end(b)
               ih = segs(iseg)%halo_idx
               pos = segs(iseg)%local_offset

               if (segs(iseg)%write_header) then
                  buf(pos:pos + 3, s) = transfer(np_halo(ih), mold_4)
                  buf(pos + 4:pos + 27, s) = transfer(xv_mean(:, ih), mold_24)
                  pos = pos + 28
               end if

               np = segs(iseg)%member_count
               if (np > 0) then
                  igroup = halo_heads(ih)
                  gpos = group_ptr(igroup)
                  member_start = segs(iseg)%member_start

                  do member = 0, np - 1
                     ip = group_members(gpos + (member_start - 1) + member)
                     buf(pos:pos + 7, s) = transfer(pid(ip), mold_8)
                     buf(pos + 8:pos + 31, s) = transfer(xv(:, ip), mold_24)
                     pos = pos + 32
                  end do
               end if
            end do
            !$omp end taskloop
         end block
         !$omp end task

         ! 2. 写入任务：依赖本批填充与严格顺序写入
         !$omp task depend(inout: slot_dep(s)) depend(inout: write_order) firstprivate(b, s)
         block
            integer(int32) :: cur_len
            cur_len = batch_bytes(b)

            write(unit) buf(1:cur_len, s)
         end block
         !$omp end task
      end do
      !$omp end single
      !$omp end parallel

      deallocate(buf)
      deallocate(segs)
      deallocate(batch_bytes, batch_seg_start, batch_seg_end)

   contains

      subroutine partition_stream(pass_num, cap, n_h, np_h, out_b, out_s, &
                                 s_arr, b_bytes, b_s_start, b_s_end)
         integer, intent(in) :: pass_num
         integer(int64), intent(in) :: cap
         integer(int32), intent(in) :: n_h
         integer(int32), intent(in) :: np_h(n_h)
         integer, intent(inout) :: out_b, out_s
         type(halo_segment), intent(inout), optional :: s_arr(:)
         integer(int32), intent(inout), optional :: b_bytes(:), b_s_start(:), b_s_end(:)

         integer(int32) :: ih, cur_b, cur_s, cur_b_bytes, rem_mem, cur_m_start, k
         integer(int32) :: space_left, space_for_mem, max_fit_mem, header_bytes, seg_bytes
         integer(int32), parameter :: max_chunk = 16384_int32
         logical :: need_head, write_head

         cur_b = 1
         cur_s = 0
         cur_b_bytes = 0

         if (pass_num == 2) then
            b_s_start(cur_b) = 1
         end if

         do ih = 1, n_h
            need_head = .true.
            rem_mem = np_h(ih)
            cur_m_start = 1

            do while (need_head .or. (rem_mem > 0))
               space_left = int(cap - int(cur_b_bytes, int64), int32)

               if (need_head) then
                  if (space_left < 28) then
                     if (pass_num == 2) then
                        b_bytes(cur_b) = cur_b_bytes
                        b_s_end(cur_b) = cur_s
                     end if
                     cur_b = cur_b + 1
                     cur_b_bytes = 0
                     space_left = int(cap, int32)
                     if (pass_num == 2) then
                        b_s_start(cur_b) = cur_s + 1
                     end if
                  end if
                  write_head = .true.
                  header_bytes = 28
                  need_head = .false.
               else
                  write_head = .false.
                  header_bytes = 0
               end if

               space_for_mem = space_left - header_bytes
               max_fit_mem = space_for_mem / 32

               if (rem_mem > 0 .and. max_fit_mem == 0 .and. (.not. write_head)) then
                  if (pass_num == 2) then
                     b_bytes(cur_b) = cur_b_bytes
                     b_s_end(cur_b) = cur_s
                  end if
                  cur_b = cur_b + 1
                  cur_b_bytes = 0
                  space_left = int(cap, int32)
                  if (pass_num == 2) then
                     b_s_start(cur_b) = cur_s + 1
                  end if
                  space_for_mem = space_left
                  max_fit_mem = space_for_mem / 32
               end if

               k = min(rem_mem, max_fit_mem)
               k = min(k, max_chunk)

               seg_bytes = header_bytes + k * 32
               cur_s = cur_s + 1

               if (pass_num == 2) then
                  s_arr(cur_s)%halo_idx = ih
                  s_arr(cur_s)%write_header = write_head
                  s_arr(cur_s)%member_start = cur_m_start
                  s_arr(cur_s)%member_count = k
                  s_arr(cur_s)%local_offset = cur_b_bytes + 1
               end if

               cur_b_bytes = cur_b_bytes + seg_bytes
               cur_m_start = cur_m_start + k
               rem_mem = rem_mem - k
            end do
         end do

         if (pass_num == 2) then
            b_bytes(cur_b) = cur_b_bytes
            b_s_end(cur_b) = cur_s
         end if

         out_b = cur_b
         out_s = cur_s
      end subroutine partition_stream

   end subroutine write_members_buffered

end module buffered_output

program CUBE_FoF
   use omp_lib
   ! 显式导入避免 variables 中的 pid 等名字与本入口的 CRS 所有权冲突。
   use variables, only: ncore,nn,nnt,ratio_cs,nc,nt,box,b_link,np_halo_min,ninfo, &
      izipx,izipi,ishift,rshift,x_resolution,nmax_redshift,n_checkpoint,z_checkpoint, &
      z_halofind,sim,head,image,rank,icx,icy,icz,inx,iny,inz,itx,ity,itz,iteam,ijk, &
      xp_new,vp_new,pid_new,opath,output_name,output_name_halo,z2str, &
      type_halo_catalog_header,type_halo_catalog_array
   use iso_c_binding, only: c_int
   use iso_fortran_env, only: file_storage_size
   use buffered_output, only: write_members_buffered
   use parameters, only: h0,omega_m,omega_l
   use bound_test, only: bound_check_batch,bound_release,a_bound,b_bound, &
      bound_xv=>xv,bound_pid=>pid,bound_member_table,bound_member_offset,bound_member_count
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

   ! 计时按 image_now 索引，容量必须与 nn**3 对应，复制自已测 clean 入口。
   integer total_images,layer_image,image_now,i,j,k,l,nlayer,layer_core,cur_checkpoint,np,n1,n2,idx(3),nh[*],im,idx1(3),idx2(3),ft(6,nn**3),ftr(nn**3)
   integer(8) iq1,iq2,iq3,i_neighbor,jq(3),ijk_neighbor(3,3),neighbor_b(3,3)
   integer(8) nhalo_all,halo_images(nn**3)[*]
   integer(8) nhalo_bound_all,bound_images(nn**3)[*]
   real(8) :: image_start,bound_start,bound_elapsed,output_start
   integer(4) nlast,ip,jp,np_iso(fofcore),np_head(fofcore),np_max,np_neighbors(3,3,3)[nn,nn,*],np_need(3,3,3),max_nei[nn,nn,*],offset_nei(3,3,3),offset_team(fofcore)
   integer(4),allocatable :: np_halo_all(:),np_halo(:),np_halo_team(:,:), rhoc_local(:,:,:,:,:,:),offset_map(:,:,:,:,:,:)
   integer(4),allocatable :: hoc(:,:,:),ll(:),llgp(:),hcgp(:),ecgp(:),iph_halo_all(:),iph_halo(:)!,h2(:,:,:,:),l2(:,:)
   real rsq,pos1(3),dx1(3),dx2(3),shift_xv(3)
   real(8) rho8,dxv(6)
   real,allocatable :: xv(:,:),xv_mean(:,:),xv_mean_team(:,:,:),rho_grid(:,:,:)
   type(type_halo_catalog_header) halo_header
   type(type_halo_catalog_array),allocatable :: hcat(:)

   ! CRS、PID 和所有权声明复制自 fof_new/clean/fof.f90，已通过300_1536全量回归。
   integer,parameter :: x_units=3*storage_size(xp_new)/file_storage_size
   integer,parameter :: v_units=3*storage_size(vp_new)/file_storage_size
   integer,parameter :: id_units=storage_size(0_izipi)/file_storage_size
   integer :: ux,uv,uid,u_members,ilo,ihi,row_start,row_count,row_capacity
   integer :: pass,icell,ip_source,ip_target,ithread,my_start,step,p,root,g,ngroup,gpos
   integer(8) :: file_length,expected_length
   logical :: files_checked=.false.
   integer,allocatable :: cell_ptr(:),cell_cursor(:,:,:),parent_size(:)
   integer,allocatable :: group_ptr(:),group_cursor(:),group_members(:),head_team(:,:),halo_heads(:)
   integer(8),allocatable :: pid(:)
   ! 标记已认领的层；并行区结束后才整理成员。
   logical(1) :: layer_done(nfof),claimed
   interface
      subroutine fof_atomic_union(parent,ii,jj) bind(C,name='fof_atomic_union')
         import c_int
         integer(c_int) :: parent(*)
         integer(c_int),value :: ii,jj
      end subroutine
   end interface

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
      ! 复制自 clean：路径函数通过 halofind 索引取本次快照红移。
      z_halofind(cur_checkpoint)=z_checkpoint(cur_checkpoint)
      if (head) print*,  ''
      if (head) print*,  ''
      if (head) print*,  'FoF at redshift ',z2str(z_checkpoint(cur_checkpoint))

      nhalo_all = 0
      halo_images=0; bound_images=0
      do image_now = this_image(),nn**3,total_images


         write(log_filename, "(A,'/fof/run_output_', I4.4, '.log')") opath,image_now
         open(newunit=log_unit, file=trim(log_filename), status='replace', action='write')
         write(*,'(A, I0, A)') '  Fof at image',image_now,'->'//trim(log_filename)
         write(log_unit,'(A, I0)') '  Fof at image',image_now


         ! init image
         call geometry_images
         image_start=omp_get_wtime()
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
                           do k = merge(1, mod(ijk_neighbor(1, 3)+1, nt), itz*nt-nt >= ijk_neighbor(1, 3)),merge(nt, mod(ijk_neighbor(2, 3)+1, nt)-1, itz*nt <= ijk_neighbor(2, 3))
                              do j = merge(1, mod(ijk_neighbor(1, 2)+1, nt), ity*nt-nt >= ijk_neighbor(1, 2)),merge(nt, mod(ijk_neighbor(2, 2)+1, nt)-1, ity*nt <= ijk_neighbor(2, 2))
                                 do i = merge(1, mod(ijk_neighbor(1, 1)+1, nt), itx*nt-nt >= ijk_neighbor(1, 1)),merge(nt, mod(ijk_neighbor(2, 1)+1, nt)-1, itx*nt <= ijk_neighbor(2, 1))
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

         ! 按cell两遍直接装配空间CRS，保持PID不变，不常驻第二份粒子数组。
         ! 复制自 fof_new/clean/fof.f90 的部分读取/CRS段，已通过300_1536全量回归。
         allocate(cell_ptr(int(nfof)**3+1),cell_cursor(nfof,nfof,nfof))
         cell_cursor=0
         nlast = 0
         do iq1 = 1, 3;  do iq2 = 1, 3; do iq3 = 1, 3
                  offset_nei(iq1,iq2,iq3) =  nlast
                  nlast = nlast + np_neighbors(iq1,iq2,iq3)
               enddo; enddo; enddo

         do pass=1,2
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
                     sim%cur_halofind=cur_checkpoint
                     if (.not.files_checked) then
                        inquire(file=output_name('xp'),size=file_length)
                        expected_length=sim%nplocal*x_units
                        if (file_length/=expected_length) error stop 'xp file size disagrees with configured layout'
                        inquire(file=output_name('vp'),size=file_length)
                        expected_length=sim%nplocal*v_units
                        if (file_length/=expected_length) error stop 'vp file size disagrees with configured layout'
                        files_checked=.true.
                     endif
                     allocate(rhoc_local(nt,nt,nt,nnt,nnt,nnt),offset_map(nt,nt,nt,nnt,nnt,nnt))
                     open(iteam,file=output_name('np'),access='stream'); read(iteam) rhoc_local; close(iteam)
                     row_capacity=maxval(sum(rhoc_local,dim=1))
                     allocate(xp_new(3,row_capacity),vp_new(3,row_capacity),pid_new(row_capacity))
                     open(newunit=ux,file=output_name('xp'),access='stream',status='old',action='read')
                     open(newunit=uv,file=output_name('vp'),access='stream',status='old',action='read')
                     open(newunit=uid,file=output_name('id'),access='stream',status='old',action='read')

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
                              do k = merge(1, mod(ijk_neighbor(1, 3)+1, nt), itz*nt-nt >= ijk_neighbor(1, 3)),merge(nt, mod(ijk_neighbor(2, 3)+1, nt)-1, itz*nt <= ijk_neighbor(2, 3))
                                 do j = merge(1, mod(ijk_neighbor(1, 2)+1, nt), ity*nt-nt >= ijk_neighbor(1, 2)),merge(nt, mod(ijk_neighbor(2, 2)+1, nt)-1, ity*nt <= ijk_neighbor(2, 2))
                                    ! 同一源行的连续 i cell 合并读取，POS 是64位文件单位偏移。
                                    ilo=merge(1,mod(ijk_neighbor(1,1)+1,nt),itx*nt-nt>=ijk_neighbor(1,1))
                                    ihi=merge(nt,mod(ijk_neighbor(2,1)+1,nt)-1,itx*nt<=ijk_neighbor(2,1))
                                    row_start=offset_map(ilo,j,k,itx,ity,itz)
                                    row_count=sum(rhoc_local(ilo:ihi,j,k,itx,ity,itz))
                                    if (row_count>0) then
                                       read(ux,pos=1_8+int(row_start,8)*x_units) xp_new(:,:row_count)
                                       if (pass==2) then
                                          read(uv,pos=1_8+int(row_start,8)*v_units) vp_new(:,:row_count)
                                          read(uid,pos=1_8+int(row_start,8)*id_units) pid_new(:row_count)
                                       endif
                                    endif
                                    do i = merge(1, mod(ijk_neighbor(1, 1)+1, nt), itx*nt-nt >= ijk_neighbor(1, 1)),merge(nt, mod(ijk_neighbor(2, 1)+1, nt)-1, itx*nt <= ijk_neighbor(2, 1))
                                       np=rhoc_local(i,j,k,itx,ity,itz)
                                       nlast = offset_map(i,j,k,itx,ity,itz)-row_start

                                       ! if (jp+np > np_max .or. nlast+np > sim%nplocal ) then
                                       !   write(log_unit, *)  iq1-2,iq2-2,iq3-2,image
                                       !   write(log_unit, *)  i,j,k,itx,ity,itz
                                       !   write(log_unit, *)  nlast,np,np_neighbors(iq1,iq2,iq3),jp
                                       !   write(log_unit, *)  'particle index error'
                                       !   error stop 'particle index error'
                                       ! endif
                                       do ip_source=nlast+1,nlast+np
#ifdef ZIPX
                                          pos1=(int(xp_new(:,ip_source)+ishift,izipx)+rshift)*x_resolution &
                                             +nt*([itx,ity,itz]-1)+[i,j,k]-1+[fof_buffer,fof_buffer,fof_buffer]+ijk_neighbor(3,:)
#else
                                          pos1=xp_new(:,ip_source)+[fof_buffer,fof_buffer,fof_buffer]+ijk_neighbor(3,:)
#endif
                                          idx=floor(pos1*n_refine)+1
                                          if (pass==1) then
                                             cell_cursor(idx(1),idx(2),idx(3))=cell_cursor(idx(1),idx(2),idx(3))+1
                                          else
                                             ip_target=cell_cursor(idx(1),idx(2),idx(3))
                                             xv(1:3,ip_target)=pos1
                                             xv(4:6,ip_target)=vp_new(:,ip_source)
                                             pid(ip_target)=pid_new(ip_source)
                                             cell_cursor(idx(1),idx(2),idx(3))=ip_target+1
                                          endif
                                       enddo
                                       jp = jp+np
                                    enddo
                                 enddo
                              enddo
                           enddo
                        enddo
                     enddo
                     close(ux); close(uv); close(uid)
                     deallocate(xp_new,vp_new,pid_new,rhoc_local,offset_map)
                     n2 = jp-n1+1
                     if (n2 /= np_neighbors(iq1,iq2,iq3) ) then
                        write(log_unit, *)  iq1-2,iq2-2,iq3-2,image
                        write(log_unit, *) np_neighbors(iq1,iq2,iq3),n2,jp,n1
                        ! 第一遍尚无 xv，且第二遍已改为 CRS 存储顺序。
                        write(log_unit, *)  'particle len error'
                        print*,image
                        stop 'particle len error'
                     endif
                     ! write(log_unit,'(5I4,I10)')  iq3+(iq2-1)*3+(iq1-1)*9,iq1-2,iq2-2,iq3-2,image,sim%nplocal
                  enddo
               enddo
            enddo
            if (pass==1) then
               cell_ptr(1)=1
               icell=0
               do k=1,int(nfof)
                  do j=1,int(nfof)
                     do i=1,int(nfof)
                        icell=icell+1
                        cell_ptr(icell+1)=cell_ptr(icell)+cell_cursor(i,j,k)
                        cell_cursor(i,j,k)=cell_ptr(icell)
                     enddo
                  enddo
               enddo
               allocate(xv(6,np_max),pid(np_max))
            endif
         enddo
         deallocate(cell_cursor)
         image = image_now


         call system_clock(ft(2,image_now),ftr(image_now))
         write(log_unit, '(A, F10.3, A)') '    real time =',real(ft(2,image_now)-ft(1,image_now))/ftr(image_now),'secs';write(log_unit, *) ''
         ! cycle



         ! 空间CRS已在读取时建立；原子DSU找链段复制自 fof_new/clean/fof.f90。
         write(log_unit, *)' initialize DSU'
         allocate(parent_size(np_max))
         parent_size=-1


         call system_clock(ft(3,image_now),ftr(image_now))
         write(log_unit, '(A, F10.3, A)') '    real time =',real(ft(3,image_now)-ft(2,image_now))/ftr(image_now),'secs';write(log_unit, *) ''

         call system_clock(ft(4,image_now),ftr(image_now))


         ! loop over fof cells
         write(log_unit, *) 'Loop over fof cells'

         ! 厚切片领完后环扫互助；同cell及13半邻域配对顺序不变。
         layer_done=.false.
         !$omp parallel num_threads(fofcore) default(shared) &
         !$omp private(ithread,my_start,step,iq1,iq2,iq3,claimed)
         ithread=omp_get_thread_num()+1
         my_start=int(int(nfof,8)*(ithread-1)/omp_get_num_threads())+1
         do step=0,int(nfof)-1
            iq3=mod(my_start-1+step,int(nfof))+1
            !$omp atomic capture
            claimed=layer_done(iq3)
            layer_done(iq3)=.true.
            !$omp end atomic
            if (claimed) cycle
            do iq2=1,nfof
               do iq1=1,nfof
                  call process_crs_cell(iq1,iq2,iq3)
               enddo
            enddo
         enddo
         !$omp end parallel
         deallocate(cell_ptr)

         call system_clock(ft(5,image_now),ftr(image_now))
         write(log_unit, '(A, F10.3, A)') '    real time =',real(ft(5,image_now)-ft(4,image_now))/ftr(image_now),'secs';write(log_unit, *) ''


         ! 组成员CRS复制自 clean：并发找链结束后恢复根大小，再保存原xv内部下标。
         do p=1,np_max
            root=find_root(p)
            if (root/=p) parent_size(root)=parent_size(root)-1
         enddo
         ngroup=0
         do p=1,np_max
            if (parent_size(p)<0) ngroup=ngroup+1
         enddo
         allocate(group_ptr(ngroup+1),group_cursor(ngroup))
         group_ptr(1)=1
         g=0
         do p=1,np_max
            if (parent_size(p)<0) then
               g=g+1
               group_ptr(g+1)=group_ptr(g)-parent_size(p)
            endif
         enddo
         do p=1,np_max
            root=find_root(p)
            if (p/=root) parent_size(p)=root
         enddo
         g=0
         do p=1,np_max
            if (parent_size(p)<0) then
               g=g+1
               parent_size(p)=-g
            endif
         enddo
         group_cursor=group_ptr(:ngroup)
         allocate(group_members(np_max))
         do p=1,np_max
            if (parent_size(p)<0) then
               g=-parent_size(p)
            else
               g=-parent_size(parent_size(p))
            endif
            group_members(group_cursor(g))=p
            group_cursor(g)=group_cursor(g)+1
         enddo
         deallocate(parent_size,group_cursor)

         ! count halos
         nlast = np_max/np_halo_min/2
         allocate(xv_mean_team(6, nlast, fofcore),np_halo_team(nlast, fofcore),head_team(nlast,fofcore))
         write(log_unit, *)  'Count halos'
         np_iso = 0; np_head = 0; np_halo_team=0; xv_mean_team=0
         !$omp parallel do private(i, iteam, ip, jp, np, dxv,gpos) default(shared)
         do i=1,ngroup
            iteam = omp_get_thread_num()+1
            np=group_ptr(i+1)-group_ptr(i)
            if (np==1) then
               np_iso(iteam)=np_iso(iteam)+1
            else
               ip = group_members(group_ptr(i)); dxv = 0
               do gpos=group_ptr(i),group_ptr(i+1)-1
                  jp=group_members(gpos)
                  dxv(1:3) = dxv(1:3) + modulo(xv(1:3,jp) - xv(1:3,ip) + L_bLb2, L_bLb*1d0)-L_bLb2
                  dxv(4:6) = dxv(4:6) + xv(4:6,jp)
               enddo
               dxv(1:3) = modulo(dxv(1:3)/np + xv(1:3,ip), L_bLb*1d0)
               dxv(4:6) = dxv(4:6)/np
               if (np > np_halo_min .and. minval(dxv(1:3)) > L_b .and. maxval(dxv(1:3)) < L_bL) then
                  np_head(iteam)=np_head(iteam)+1
                  np_halo_team(np_head(iteam),iteam) = np
                  head_team(np_head(iteam),iteam) = i
                  xv_mean_team(1:3,np_head(iteam),iteam) =  (dxv(1:3) + shift_xv)*box/nn/nc
                  xv_mean_team(4:6,np_head(iteam),iteam) =  dxv(4:6)
               endif
            endif
         enddo
         !$omp end parallel do

         offset_team = 0; nlast = 0
         do i = 1,fofcore
            np_iso(1)  = np_iso(1) + np_iso(i)
            offset_team(i) = nlast
            nlast = nlast + np_head(i)
         enddo

         allocate(np_halo(nlast),halo_heads(nlast))
         do i = 1,fofcore
            np_halo(offset_team(i)+1:offset_team(i)+np_head(i)) = np_halo_team(1:np_head(i),i)
            halo_heads(offset_team(i)+1:offset_team(i)+np_head(i)) = head_team(1:np_head(i),i)
         enddo
         deallocate(np_halo_team,head_team)

         allocate(xv_mean(6,nlast))
         do i = 1,fofcore
            xv_mean(:, offset_team(i)+1:offset_team(i)+np_head(i)) = xv_mean_team(:, 1:np_head(i),i)
         enddo
         deallocate(xv_mean_team)


         np_head(1) = nlast
         np_head(2) = sum(np_halo)
         ! 原格式完整成员双缓冲输出复制自 clean，不改变PID、相空间或组顺序。
         output_start=omp_get_wtime()
         open(newunit=u_members,file=output_name_halo('members'),access='stream',status='replace')
         write(u_members) nlast,np_head(2)
         call write_members_buffered(u_members,nlast,np_halo,xv_mean,halo_heads, &
            group_ptr,group_members,pid,xv,fofcore)
         close(u_members)
         write(log_unit,'(A,F14.4)') '[FOF_MEMBER_OUTPUT] seconds=',omp_get_wtime()-output_start
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
         if (nh>0) write(log_unit, *)'halo_header',halo_header,maxval(hcat%hmass)


         write(log_unit, *) output_name_halo(trim('halo_'//trim(adjustl(str_refine))))
         open(21,file=output_name_halo(trim('halo_'//trim(adjustl(str_refine)))),status='replace',access='stream')
         write(21) halo_header
         write(21) hcat
         close(21)

         ! 直接借用组成员CRS完成去束缚，随后写独立目录；接入与输出复制自 clean。
         bound_start=omp_get_wtime()
         call write_bound_catalog
         bound_elapsed=omp_get_wtime()-bound_start
         write(log_unit,'(A,F14.4)') '[BOUND_AND_OUTPUT] seconds=',bound_elapsed
         write(log_unit,'(A,F14.4)') '[IMAGE_TOTAL] seconds=',omp_get_wtime()-image_start
         write(*,'(A,F14.4)') '[IMAGE_TOTAL] seconds=',omp_get_wtime()-image_start
         deallocate(xv,pid)

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
            bound_images = bound_images + bound_images(:)[image]
         enddo
         nhalo_all = sum(halo_images)
         nhalo_bound_all = sum(bound_images)
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
            ! 与已测 clean 一致：不访问空目录或小目录中不存在的调试行。
            write(111) hcat
            deallocate(hcat)
         enddo
         close(111)
         ! 束缚目录有独立计数和header，复制自 clean，不覆盖原FoF目录。
         image=0
         halo_header%nhalo_tot=int(nhalo_bound_all)
         halo_header%nhalo=int(nhalo_bound_all)
         halo_header%ninfo=ninfo; halo_header%linking_parameter=b_link
         open(111,file=output_name_halo(trim('halo_bound_all_'//trim(adjustl(str_refine)))), &
            status='replace',access='stream')
         write(111) halo_header
         do image=1,nn**3
            open(21,file=output_name_halo(trim('halo_bound_'//trim(adjustl(str_refine)))), &
               status='old',access='stream',action='read')
            read(21) halo_header
            allocate(hcat(halo_header%nhalo)); read(21) hcat; close(21)
            write(111) hcat
            deallocate(hcat)
         enddo
         close(111)
      endif




   enddo ! cur_checkpoint



contains

   ! 下列接入、输出与DSU遍历复制自 fof_new/clean/fof.f90，仅去除重复内部断言。
   subroutine write_bound_catalog
      integer :: h,g,nbh,m,u_status,u_catalog,u_bound
      integer,allocatable :: heads(:),nb(:),status(:),final_members(:),result_ptr(:), &
         kept_heads(:),kept_counts(:)
      real,allocatable :: xm(:,:),vm(:,:),xb(:,:),vb(:,:),kept_means(:,:)
      real(8) :: t0,hubble,begin_bound,begin_output
      integer(8) :: total_members
      type(type_halo_catalog_header) :: header_bound
      type(type_halo_catalog_array),allocatable :: catalog_bound(:)

      allocate(heads(nh),bound_member_offset(nh),bound_member_count(nh), &
         nb(nh),status(nh),xm(3,nh),vm(3,nh),xb(3,nh),vb(3,nh),result_ptr(nh+1))
      total_members=sum(int(np_halo,8))
      if (total_members>=int(huge(0),8)) error stop 'FoF: bound output index overflow'
      result_ptr(1)=1
      do h=1,nh
         g=halo_heads(h)
         bound_member_offset(h)=group_ptr(g)
         bound_member_count(h)=np_halo(h)
         heads(h)=group_ptr(g)
         result_ptr(h+1)=result_ptr(h)+np_halo(h)
      enddo
      allocate(final_members(int(total_members)))
      xm=xv_mean(1:3,:)*real(nn*nc/box,4)-spread(shift_xv,2,int(nh))
      vm=xv_mean(4:6,:)
      t0=1d0/(150d0*dble(h0)*sqrt(dble(omega_m)))
      a_bound=dble(sim%a)
      hubble=100d0*dble(h0)*sqrt(dble(omega_m)/a_bound**3+dble(omega_l))
      b_bound=hubble*a_bound**2*t0
      deallocate(halo_heads,group_ptr)
      call move_alloc(xv,bound_xv)
      call move_alloc(pid,bound_pid)
      call move_alloc(group_members,bound_member_table)
      begin_bound=omp_get_wtime()
      call bound_check_batch(heads,np_halo,xm,vm,nb,xb,vb,status,configuration=8,final_members=final_members)
      write(log_unit,'(A,F14.4)') '[BOUND_BATCH] seconds=',omp_get_wtime()-begin_bound
      call bound_release
      ! 写第二套双缓冲之前释放完整组CRS；xv/PID所有权返回，数据不改写。
      deallocate(bound_member_table,bound_member_offset,bound_member_count)
      call move_alloc(bound_xv,xv)
      call move_alloc(bound_pid,pid)
      begin_output=omp_get_wtime()
      nbh=count(status==0.or.status==5)
      allocate(kept_heads(nbh),kept_counts(nbh),kept_means(6,nbh),catalog_bound(nbh))
      open(newunit=u_status,file=output_name_halo('bound_status'),access='stream',status='replace')
      ! int32候选数，逐候选五个int32；status5保留并可追溯，不宣称收敛。
      write(u_status) int(nh)
      m=0
      do h=1,nh
         g=0
         if (status(h)==0.or.status(h)==5) then
            m=m+1; g=m
            kept_heads(m)=h; kept_counts(m)=nb(h)
            kept_means(1:3,m)=(xb(:,h)+shift_xv)*box/nn/nc
            kept_means(4:6,m)=vb(:,h)
            catalog_bound(m)%hmass=real(nb(h),4)
            catalog_bound(m)%xv=kept_means(:,m)
         endif
         write(u_status) h,g,status(h),np_halo(h),nb(h)
      enddo
      close(u_status)
      header_bound%nhalo_tot=0; header_bound%nhalo=nbh
      header_bound%ninfo=ninfo; header_bound%linking_parameter=b_link
      open(newunit=u_catalog,file=output_name_halo(trim('halo_bound_'//trim(adjustl(str_refine)))), &
         access='stream',status='replace')
      write(u_catalog) header_bound
      write(u_catalog) catalog_bound
      close(u_catalog)
      total_members=sum(int(kept_counts,8))
      if (total_members>int(huge(0),8)) error stop 'FoF: bound member file count overflow'
      open(newunit=u_bound,file=output_name_halo('members_bound'),access='stream',status='replace')
      write(u_bound) nbh,int(total_members)
      call write_members_buffered(u_bound,nbh,kept_counts,kept_means,kept_heads, &
         result_ptr,final_members,pid,xv,fofcore)
      close(u_bound)
      bound_images(image_now)=nbh
      write(log_unit,'(A,3(I0,1X),A,I0)') '[BOUND_STATUS] OK/TOO_FEW/MAX_ITER=', &
         count(status==0),count(status==1),count(status==5),' written_members=',total_members
      write(log_unit,'(A,F14.4)') '[BOUND_OUTPUT] seconds=',omp_get_wtime()-begin_output
   end subroutine write_bound_catalog

   subroutine process_crs_cell(iq1,iq2,iq3)
      integer(8),intent(in) :: iq1,iq2,iq3
      integer :: icell,jcell,ip,jp,ip_begin,ip_end,jp_begin,jp_end,i_neighbor
      integer(8) :: jq(3)
      real :: rsq
      icell=int(iq1+nfof*(iq2-1+nfof*(iq3-1)))
      ip_begin=cell_ptr(icell); ip_end=cell_ptr(icell+1)-1
      do ip=ip_end,ip_begin,-1
         do jp=ip-1,ip_begin,-1
            rsq=sum((xv(1:3,ip)-xv(1:3,jp))**2)
            if (rsq<=rp2) call fof_atomic_union(parent_size,ip,jp)
         enddo
         do i_neighbor=1,13
            jq = [iq1,iq2,iq3]+ijk(:,i_neighbor)
            if (maxval(jq)>nfof .or. minval(jq)<1) cycle
            jcell=int(jq(1)+nfof*(jq(2)-1+nfof*(jq(3)-1)))
            jp_begin=cell_ptr(jcell); jp_end=cell_ptr(jcell+1)-1
            do jp=jp_end,jp_begin,-1
               rsq=sum((xv(1:3,ip)-xv(1:3,jp))**2)
               if (rsq<=rp2) call fof_atomic_union(parent_size,ip,jp)
            enddo
         enddo
      enddo
   end subroutine

   integer function find_root(ip) result(iroot)
      integer,intent(in) :: ip
      integer :: ipart,inext
      iroot=ip
      do while (parent_size(iroot)>0)
         iroot=parent_size(iroot)
      enddo
      ipart=ip
      do while (ipart/=iroot)
         inext=parent_size(ipart)
         parent_size(ipart)=iroot
         ipart=inext
      enddo
   end function

   subroutine geometry_images
      write(log_unit, *) 'geometry_images'
      rank=image_now-1               ! MPI_rank
      icz=rank/(nn**2)+1             ! image_z
      icy=(rank-nn**2*(icz-1))/nn+1  ! image_y
      icx=mod(rank,nn)+1             ! image_x
   endsubroutine

end
