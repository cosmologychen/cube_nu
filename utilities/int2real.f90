program int2real
  use omp_lib
  use variables
  implicit none
  integer,parameter:: real_images = 8
  integer,parameter:: layer_image = nn**3/real_images

  integer :: my_id          ! 当前 image 的 ID
  integer :: log_unit       ! 文件单元号
  character(len=64) :: log_filename ! 文件名字符串
  
  integer image_now,i,j,k,l,nlayer,cur_checkpoint,np,n1,n2,idx(3)
  integer(4) nlast,ip,jp
  integer(4),allocatable :: rhoc_local(:,:,:,:,:,:)
  real,allocatable :: xv(:,:)
  integer(4), allocatable :: offset_map(:,:,:,:,:,:)
  allocate(offset_map(nt,nt,nt,nnt,nnt,nnt))
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (real_images*layer_image /= nn**3) then
    if(head) print*,real_images,layer_image,nn**3
    stop 'real_images*layer_image /= nn*3'
  endif

  head=(this_image()==1)
  call omp_set_num_threads(ncore)


  if (head) then
    open(16,file='z_checkpoint.txt',status='old')
    do i=1,nmax_redshift-1
      read(16,end=71,fmt='(f8.4)') z_checkpoint(i)
    enddo
    71 n_checkpoint=i-1
    close(16)
    if (n_checkpoint==0) stop 'z_checkpoint.txt empty'
  endif
  sync all
  n_checkpoint=n_checkpoint[1]
  z_checkpoint(:)=z_checkpoint(:)[1]
  sync all

  do cur_checkpoint=n_checkpoint,n_checkpoint
  ! do cur_checkpoint=1,1
    sim%cur_checkpoint=cur_checkpoint
    if (head) print*,  ''
    if (head) print*,  ''
    if (head) print*,  'int2real at redshift ',z2str(z_checkpoint(cur_checkpoint))

    do image = this_image(),nn**3,real_images



      write(log_filename, "('int2real_error_', I4.4, '.txt')") image
      open(newunit=log_unit, file=trim(log_filename), status='replace', action='write')
      open(111,file=output_name('info'),access='stream'); read(111) sim; close(111)
      sim%cur_checkpoint=cur_checkpoint
      write(*,'(I4,I10)')  image,sim%nplocal
      allocate(xv(3,sim%nplocal),xp_new(3,sim%nplocal),rhoc_local(nt,nt,nt,nnt,nnt,nnt))
      open(111,file=output_name('xp'),access='stream'); read(111) xp_new; close(111)
      open(111,file=output_name('np'),access='stream'); read(111) rhoc_local; close(111)
      nlast = 0
      ! 循环遍历各个维度
      call system_clock(t1,t_rate)
      !$omp parallel do schedule(dynamic,1) &
      !$omp default(shared) &
      !$omp private(i,j,k,np,nlast,itx,ity,itz)
      do itz = 1,nnt
      do ity = 1,nnt
      do itx = 1,nnt
        do k = 1,nt
        do j = 1,nt
        do i = 1,nt
          np=rhoc_local(i,j,k,itx,ity,itz)
          nlast = offset_map(i,j,k,itx,ity,itz)
          
          if (np < 0 .or. nlast+np > sim%nplocal  .or. nlast < 0) then
            write(log_unit, '(I4,I10,I10)')         image,nlast,np
            write(log_unit, '(6I4)')  i,j,k,itx,ity,itz
            write(log_unit, *)  'particle index error'
            close(log_unit)
            error stop 'particle index error'
          endif
          
          xv(:,nlast+1:nlast+np)=(int(xp_new(:,nlast+1:nlast+np)+ishift,izipx)+rshift)*x_resolution &
                                +spread(nt*((/itx,ity,itz/)-1)+((/i,j,k/)-1),dim=2,ncopies=np)
  
        enddo
        enddo
        enddo
        ! write(log_unit, '(4I4,3F7.3)')  itz,itx,ity,1,sum(xv(1,nlast+1:nlast+np))/np,minval(xv(1,nlast+1:nlast+np)),maxval(xv(1,nlast+1:nlast+np))
        ! write(log_unit, '(4I4,3F7.3)')  itz,itx,ity,2,sum(xv(2,nlast+1:nlast+np))/np,minval(xv(2,nlast+1:nlast+np)),maxval(xv(2,nlast+1:nlast+np))
        ! write(log_unit, '(4I4,3F7.3)')  itz,itx,ity,3,sum(xv(3,nlast+1:nlast+np))/np,minval(xv(3,nlast+1:nlast+np)),maxval(xv(3,nlast+1:nlast+np))
      enddo
      enddo
      enddo
      !$omp end parallel do
      call system_clock(t2,t_rate)
      print*,image,'    real time =',real(t2-t1)/t_rate,'secs';print*, ''
      deallocate(xp_new,rhoc_local)
      
      close(log_unit)
      write(*, *) '  int2real at image',image,output_name('xr')
      ! write(*, *) image,sum(xv(1,:))/sim%nplocal,minval(xv(1,:)),maxval(xv(1,:))
      ! write(*, *) image,sum(xv(2,:))/sim%nplocal,minval(xv(2,:)),maxval(xv(2,:))
      ! write(*, *) image,sum(xv(3,:))/sim%nplocal,minval(xv(3,:)),maxval(xv(3,:))
      open(111,file=output_name('xr'),status='replace',access='stream'); write(111) xv; close(111)
    enddo
  enddo
end program