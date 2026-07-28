subroutine geometry
   sync all
   image=this_image()
   rank=image-1            ! MPI_rank
   icz=rank/(nn**2)+1             ! image_z
   icy=(rank-nn**2*(icz-1))/nn+1  ! image_y
   icx=mod(rank,nn)+1             ! image_x
   m1=icx ! pencil_fft convension
   m2=icy
   m3=icz
   inx=modulo(icx-2,nn)+1  ! adjacent images
   iny=modulo(icy-2,nn)+1
   inz=modulo(icz-2,nn)+1
   ipx=modulo(icx,nn)+1
   ipy=modulo(icy,nn)+1
   ipz=modulo(icz,nn)+1
   head=(this_image()==1)
endsubroutine

subroutine tic(i)
   integer i
   sync all
   if (this_image()==1) call system_clock(tictoc(1,i),t_rate)
   sync all
endsubroutine

subroutine toc(i)
   integer i
   sync all
   if (this_image()==1) then
      call system_clock(tictoc(2,i),t_rate)
      tcat(i,istep)=real(tictoc(2,i)-tictoc(1,i))/t_rate
   endif
   sync all
endsubroutine

pure function image2str(nimage)
   character(:),allocatable :: image2str
   character(20) :: str
   integer(8),intent(in) :: nimage
   write(str,'(i6)') nimage
   image2str=trim(adjustl(str))
endfunction

pure function z2str(z)
   character(:),allocatable :: z2str
   character(20) :: str
   real,intent(in) :: z
   write(str,'(f7.3)') z
   z2str=trim(adjustl(str))
endfunction

function output_dir()
   character(:),allocatable :: output_dir
   character(20) :: str_i
   write(str_i,'(i6)') image
   output_dir=opath//'image'//trim(adjustl(str_i))//'/'
endfunction

function fn_hbt(n_snap,n_rank,n_kind)
   ! snapdir_101/snapshot_101.*
   ! groups_101/group_tab_101.*;
   ! groups_101/group_ids_101.*;
   ! 其中*为子文件序号，从0开始。 -- Jiaxin
   integer(4) n_snap,n_rank,n_kind
   character(20) :: str_snap, str_rank
   character(:),allocatable :: fn_hbt
   write(str_snap,'(i6)') n_snap
   write(str_rank,'(i6)') n_rank

   if (n_kind==0) then ! snapdir_
      fn_hbt=opath//'hbt/snapdir_'//trim(adjustl(str_snap))//'/'
   elseif (n_kind==10) then ! snapdir_
      fn_hbt=opath//'hbt/snapdir_'//trim(adjustl(str_snap))//'/snapshot_'//trim(adjustl(str_snap))//'.'//trim(adjustl(str_rank))

   elseif (n_kind==1) then ! groups_
      fn_hbt=opath//'hbt/groups_'//trim(adjustl(str_snap))//'/'
   elseif (n_kind==11) then ! groups_
      fn_hbt=opath//'hbt/groups_'//trim(adjustl(str_snap))//'/group_tab_'//trim(adjustl(str_snap))//'.'//trim(adjustl(str_rank))

   elseif (n_kind==2) then
      fn_hbt=opath//'hbt/groups_'//trim(adjustl(str_snap))//'/'
   elseif (n_kind==12) then
      fn_hbt=opath//'hbt/groups_'//trim(adjustl(str_snap))//'/group_ids_'//trim(adjustl(str_snap))//'.'//trim(adjustl(str_rank))

   else
      if (head) then
         print*, 'invalid n_kind',n_kind
         stop
      endif
   endif
   ! 如果需要rank号写三位数:
   ! write(str_snap,'(i6)') cur_checkpoint+100000; write(str_rank,'(i6)') rank
   !  if (head) print*, '  fn= ', output_dir()//'snapshot_'//str_snap(4:6)//'.'//trim(adjustl(str_rank))
endfunction

function output_prefix()
   character(:),allocatable :: output_prefix
   character(20) :: str_z,str_i
   write(str_i,'(i6)') image
   write(str_z,'(f7.3)') z_checkpoint(sim%cur_checkpoint)
   output_prefix=opath//'image'//trim(adjustl(str_i))//'/'//trim(adjustl(str_z))//'_'
endfunction

function output_prefix_halo()
   character(:),allocatable :: output_prefix_halo
   character(20) :: str_z,str_i
   write(str_i,'(i6)') image
   write(str_z,'(f7.3)') z_halofind(sim%cur_halofind)
   output_prefix_halo=opath//'fof/'//trim(adjustl(str_z))//'_'
endfunction

pure function output_suffix_halo()
   character(:),allocatable :: output_suffix_halo
   character(20) :: str_i
   write(str_i,'(i6)') image
   if (image == 0 ) then
      output_suffix_halo='.bin'
   else
      output_suffix_halo='_'//trim(adjustl(str_i))//'.bin'
   endif
endfunction

pure function output_suffix()
   character(:),allocatable :: output_suffix
   character(20) :: str_i
   write(str_i,'(i6)') image
   output_suffix='_'//trim(adjustl(str_i))//'.bin'
endfunction

function output_name(zipname)
   character(*) ::  zipname
   character(:),allocatable :: output_name
   output_name=output_prefix()//zipname//output_suffix()
endfunction

function output_name_halo(zipname)
   character(*) ::  zipname
   character(:),allocatable :: output_name_halo
   output_name_halo=output_prefix_halo()//zipname//output_suffix_halo()
endfunction

elemental real function F_ra(r,apm)
   real,intent(in) :: r,apm
   real ep
   F_ra=0
   ep=2*r/apm
   if (apm==0 .or. ep>2) then
      F_ra=r**(-2)
   elseif (ep>1) then
      F_ra=(1./35./apm**2)*(12*ep**(-2) - 224 + 896*ep - 840*ep**2 + 224*ep**3 + 70*ep**4 - 48*ep**5 + 7*ep**6)
   elseif (ep>0) then
      F_ra=(1./35./apm**2)*(224*ep - 224*ep**3 + 70*ep**4 + 48*ep**5 - 21*ep**6)
   else
      F_ra=0
   endif
endfunction
