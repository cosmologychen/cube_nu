program write_gadget_snapshot_and_find_fof_halos
  use variables
  implicit none
  integer i,j,k,l,cur_checkpoint,itile,np,nptile,nfof,n_refine,n1,n2,idx(3),nh_tile,nh[*],ifile
  integer iq1,iq2,iq3,i_neighbor,jq(3),npglobal_2(2),nhlocal,nhglobal,np_offset,nhlocal_max,nid,nfile[*],nfile_global[*],nfile_start(nn**3)[*]
  integer np_iso,np_mem,np_head,np_last
  integer(8) nlast,nzero,ip,jp,npglobal_long,npglobal[*],npfile_max,p1,p2
  integer(4),allocatable :: np_halo_all(:),np_halo(:),offset(:),mass_halo(:)
  integer(izipi),allocatable :: hoc(:,:,:),ll(:),llgp(:),hcgp(:),ecgp(:),iph_halo_all(:),iph_halo(:),pid_halo(:),pid_list(:),pid_last
  real rp2,rsq,xf_hoc(3),pos1(3),box_real
  real,allocatable :: xv(:,:),dxv(:,:),xv_mean(:,:)
  character(6) str_snap,str_rank
  type(type_halo_catalog_header) halo_header
  type(type_halo_catalog_array),allocatable :: hcat(:)
  type gadget_snapshot_header
    integer(4) Npart(6)
    real(8) Massarr(6),Time,Redshift
    integer(4) Flagsfr,Flagfbk,Nall(6),FlagCooling,NumFiles
    real(8) BoxSize,Omega0,OmegaLambda,HubbleParam
    integer(4) FlagAge,FlagMetals,NallHW(6),flag_entr_ics,unused(15)
  endtype
  type(gadget_snapshot_header) io_header
  equivalence(npglobal_long,npglobal_2)
  logical first_id

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call geometry
  if (head) then
    print*,''; print*, 'Gadget format and FoF finder. Checkpoint at:'
    open(16,file='z_checkpoint.txt',status='old')
!    do i=1,nmax_redshift-1
    do i=1,nmax_redshift
      read(16,end=71,fmt='(f8.4)') z_checkpoint(i)
    enddo
!    71 n_checkpoint=i-1; close(16)
	71 n_checkpoint=i-1; close(16)
  endif
  sync all
  n_checkpoint=n_checkpoint[1]; z_checkpoint(:)=z_checkpoint(:)[1]

  if (head) print*,'  initialize FoF cell neighbors'
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
  !do i=1,13;  print*,ijk(:,i); enddo

  if (head) print*,'  initialize tile index'
  ipm2=0;
  do itz=1,nnt
  do ity=1,nnt
  do itx=1,nnt
    ipm2=ipm2+1
    ixyz2(:,ipm2)=[itx,ity,itz]
  enddo
  enddo
  enddo

  n_refine=1
  nfof=nt*n_refine*ratio_cs
  n1=1-ncb*n_refine*ratio_cs; n2=(nt+ncb)*n_refine*ratio_cs
  rp2=b_link**2

  npfile_max=512**3 ! max np in one snapshot file
  !do cur_checkpoint=n_checkpoint,n_checkpoint
  do cur_checkpoint=n_checkpoint,n_checkpoint
    ! filename test
    !print*, fn_hbt(cur_checkpoint,int(rank,4),0)
    !print*, fn_hbt(cur_checkpoint,int(rank,4),1)
    !print*, fn_hbt(cur_checkpoint,int(rank,4),2)

    sim%cur_checkpoint=cur_checkpoint
    if (head) print*, 'particle initialization at ',z2str(z_checkpoint(cur_checkpoint))
    if (head) print*, output_name('info')
    call particle_initialization
    deallocate(xp_new,vp_new,pid_new)
    
    allocate(xv(6,sim%nplocal)) ! float format xv
    nlast=0
    do itz=1,nnt
    do ity=1,nnt
    do itx=1,nnt
      do k=1,nt
      do j=1,nt
      do i=1,nt
        np=rhoc(i,j,k,itx,ity,itz)
        do l=1,np
          ip=nlast+l
          pos1=nt*([itx,ity,itz]-1)+ ([i,j,k]-1) + (int(xp(:,ip)+ishift,izipx)+rshift)*x_resolution
          xv(1:3,ip)=pos1*real(ng)/real(nc)+([m1,m2,m3]-1)*ng
          xv(4:6,ip)=vp(:,ip)
          ! for velocity of integer format
          !xv(4:6,ip)=vc(:,i,j,k,itx,ity,itz)+tan((pi*real(vp(:,ip)))/real(nvbin-1))/(sqrt(pi/2)/(sim%sigma_vi*vrel_boost))*sim%vsim2phys
        enddo
        nlast=nlast+np
      enddo
      enddo
      enddo
    enddo
    enddo
    enddo
    if (nlast/=sim%nplocal) then
      print*, 'nlast/=sim%nplocal', nlast, sim%nplocal
      stop
    endif

    if (sim%nplocal>2**31-1) then
      print*, '  warning: sim%nplocal=',sim%nplocal,' > 2^31-1'
    endif
    !npglobal=sim%npglobal

    sync all; npglobal=0
    if (head) then ! 只考虑少数节点的数据的情况
      do i=1,nn**3
        npglobal=npglobal+sim[i]%nplocal
      enddo
      print*,'sim%npglobal, npglobal =',sim%npglobal, npglobal
    endif;    sync all
    npglobal=npglobal[1]
    npglobal_long=npglobal

    if (npglobal>2**31-1 .and. head) then
      print*,'  notice: npglobal =',npglobal,' > 2^31-1'
      print*,'  npglobal LW, HW =',npglobal_2
    endif
    ! convert units
    box_real=2000 ! 原始模拟的box
    ! box_real=600 ! 原始模拟的box
    io_header%Npart=[0,int(sim%nplocal,kind=4),0,0,0,0]
    io_header%Massarr=[0.,27.75*omega_m*(box_real**3)/sim%npglobal,0.,0.,0.,0.] 
    io_header%Time=sim%a; io_header%Redshift=1./sim%a-1
    io_header%FlagSfr=0; io_header%Flagfbk=0
    io_header%Nall=[0,npglobal_2(1),0,0,0,0]
    io_header%FlagCooling=0
    io_header%NumFiles=nn**3; io_header%BoxSize=box*1000
    io_header%Omega0=omega_m; io_header%OmegaLambda=omega_l
    io_header%HubbleParam=h0
    io_header%FlagAge=0; io_header%FlagMetals=0
    io_header%NallHW=npglobal_2(2); io_header%flag_entr_ics=0; io_header%unused=0
    xv(1:3,:)=xv(1:3,:)*box*1000/real(ng*nn) ! positions in kpc/h, real(4)
    xv(4:6,:)=xv(4:6,:)/sqrt(sim%a) ! velocities in km/sec/sqrt(a), real(4)
    if (head) then
      print*, 'sim%nplocal =',sim%nplocal
      print*, 'first particle: x,y,z [kpc/h]; vx,vy,vz [km/s/sqrt(a)]'
      print*, xv(:,1)
      print*, 'last particle: x,y,z [kpc/h]; vx,vy,vz [km/s/sqrt(a)]'
      print*, xv(:,sim%nplocal)
      print*, 'first PID =', pid(1)
      print*, 'last PID  =', pid(sim%nplocal)
    endif

    ! 拆分成多个子文件
    nfile=ceiling(real(sim%nplocal)/npfile_max); nfile_global=0; sync all
    if (head) then
      do i=1,nn**3
        nfile_start(i)=nfile_global
        nfile_global=nfile_global+nfile[i]
      enddo
    endif
    sync all
    nfile_start=nfile_start(:)[1]
    nfile_global=nfile_global[1]

    if (head) print*, 'write snap files'
    call system('mkdir -p '//fn_hbt(cur_checkpoint,int(0,4),0))   

    io_header%NumFiles=nfile_global
    do ifile=1,nfile
      p1=1+(ifile-1)*npfile_max; p2=min(ifile*npfile_max,sim%nplocal)
      io_header%Npart(2)=p2-p1+1
      open(11,file=fn_hbt(cur_checkpoint,int(nfile_start(image)+ifile-1,4),10),status='replace',form='unformatted')
      write(11) io_header ! 1st block: 256 byte header
      write(11) xv(1:3,p1:p2) ! 2nd block: positions
      write(11) xv(4:6,p1:p2) ! 3rd block: velocities
      write(11) pid(p1:p2)    ! 4th block: PIDs
      close(11)
    enddo

    sync all
    if (head) print*, 'wrote snap files'
    deallocate(xv)
    if (head) print*, 'done'
    sync all
    !deallocate(xp,vp,rhoc,pid); cycle





    if (head) then
      print*, ''
      print*, 'FoF at redshift ',z2str(z_checkpoint(cur_checkpoint))
      print*, 'read checkpoint header',output_name('info')
    endif
    
    call buffer_grid
    call buffer_x
    call buffer_v

    halo_header%nhalo_tot=0; nh=0; halo_header%nhalo=0; 
    halo_header%ninfo=ninfo; halo_header%linking_parameter=b_link
    if (head) print*, output_name('halo') ! CUBE format halo catalog
    open(21,file=output_name('halo'),status='replace',access='stream')
    write(21) halo_header
    
    ! group_ids
    call system('mkdir -p '//fn_hbt(cur_checkpoint,int(rank,4),2))
    open(102,file=fn_hbt(cur_checkpoint,int(rank,4),12),status='replace',access='stream')
    write(102) nhlocal,nid,nhglobal,nfile ! will rewind and rewrite
    first_id=.true.

    nhlocal_max=sim%nplocal/100
    allocate(offset(nhlocal_max),mass_halo(nhlocal_max))
    nhlocal=0; np_offset=0; nid=0
    do itile=1,nnt**3 ! work on each tile
      nlast=0
      if (head) print*, ixyz2(:,itile)
      nptile=sum(rhoc(:,:,:,ixyz2(1,itile),ixyz2(2,itile),ixyz2(3,itile)))
      !print*,nfof,n1,n2,nptile
      allocate(xv(6,nptile),pid_list(nptile))
      do k=1-ncb,nt+ncb
      do j=1-ncb,nt+ncb
      do i=1-ncb,nt+ncb
        np=rhoc(i,j,k,ixyz2(1,itile),ixyz2(2,itile),ixyz2(3,itile))
        nzero=idx_b_r(j,k,ixyz2(1,itile),ixyz2(2,itile),ixyz2(3,itile)) &
              -sum(rhoc(i:,j,k,ixyz2(1,itile),ixyz2(2,itile),ixyz2(3,itile)))
        do l=1,np
          ip=nzero+l; jp=nlast+l
          xv(1:3,jp)=ratio_cs*([i,j,k]-1+(int(xp(:,ip)+ishift,izipx)+rshift)*x_resolution) ! coarse grid
          xv(4:6,jp)=vp(:,ip)
          pid_list(jp)=pid(ip)
        enddo
        nlast=nlast+np
      enddo
      enddo
      enddo
      !print*,minval(xv),maxval(xv),nlast;

      ! create hoc ll
      allocate(hoc(n1:n2,n1:n2,n1:n2),ll(nptile),llgp(nptile),hcgp(nptile),ecgp(nptile))
      hoc=0; ll=0
      do ip=1,nptile
        idx=ceiling(xv(1:3,ip)*n_refine) ! index of the grid
        if (minval(idx)<n1 .or. maxval(idx)>n2) then
          print*, 'idx out of range'
          print*, idx
          stop
        endif
        ll(ip)=hoc(idx(1),idx(2),idx(3)) ! linked list 
        hoc(idx(1),idx(2),idx(3))=ip ! head of chain
        hcgp(ip)=ip ! initialize hcgp(ip)=ip for isolated particles
      enddo
      llgp=0; ecgp=0; ! initialize group link list
      if (head) print*, 'percolation'

      ! loop over fof cells
      do iq3=n1,n2
      do iq2=n1,n2
      do iq1=n1,n2
        ip=hoc(iq1,iq2,iq3)
        do while (ip/=0)
          jp=ll(ip)
          do while (jp/=0)
            rsq=sum((xv(1:3,ip)-xv(1:3,jp))**2)
            if (rsq<=rp2) call merge_chain(ip,jp)
            jp=ll(jp)
          enddo
          do i_neighbor=1,13
            jq=[iq1,iq2,iq3]+ijk(:,i_neighbor)
            if (minval(jq)<n1 .or. maxval(jq)>n2) cycle
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


      ! count isolated, member and leader particles
      np_iso=0; np_mem=0; np_head=0;
      do i=1,nptile
        if (hcgp(i)==i) then
          np_iso=np_iso+1
        elseif (hcgp(i)==0) then
          np_mem=np_mem+1
        else
          np_head=np_head+1
        endif
      enddo
      !print*,'N_iso,mem,head =',np_iso,np_mem,np_head
      allocate(pid_halo(np_mem+np_head))

      !print*,'fof particles'
      allocate(xv_mean(6,np_head),iph_halo_all(np_head),dxv(6,np_head),np_halo_all(np_head))

      nh_tile=0; ! count halo number in this tile
      do ip=1,nptile
        if (hcgp(ip)/=ip .and. hcgp(ip)/=0) then
          np=0; jp=hcgp(ip) ! hoc of the group
          do while (jp/=0) ! loop over group members
            np=np+1; dxv(:,np)=xv(:,jp); pid_halo(np)=pid_list(jp)
            jp=llgp(jp) ! next particle in chain
          enddo
          if (np>=np_halo_min .and. minval(sum(dxv(1:3,:np),2)/np)>=0 .and. & 
                                    maxval(sum(dxv(1:3,:np),2)/np)<nt*ratio_cs) then
            nh_tile=nh_tile+1
            nhlocal=nhlocal+1
            !print*,'nhlocal=',nhlocal
            np_halo_all(nh_tile)=np
            iph_halo_all(nh_tile)=hcgp(ip) ! hoc of the halo "ip-header"
            xv_mean(:,nh_tile)=sum(dxv(:,:np),2)/np+(ixyz2(:,itile)-1)*ngp+([icx,icy,icz]-1)*ng
            
            mass_halo(nhlocal)=np
            offset(nhlocal)=np_offset
            !print*,'np_offset',np_offset
            np_offset=np_offset+np
            nid=nid+np
            write(102) pid_halo(:np)
            pid_last=pid_halo(np)
            np_last=np
            if (first_id .and. head) then
              print*,'first id =', pid_halo(1)
              first_id=.false.
            endif
          endif
        endif
      enddo
      
      !print*,'nh_tile',nh_tile
      if (head) print*, 'halo position range:',minval(xv_mean(:3,:nh_tile)),maxval(xv_mean(:3,:nh_tile))

      ! transfer data to smaller arrays
      allocate(hcat(nh_tile),iph_halo(nh_tile)); 
      iph_halo=iph_halo_all(:nh_tile); hcat%hmass=np_halo_all(:nh_tile)
      do i=1,6
        hcat%xv(i)=xv_mean(i,1:nh_tile)
      enddo
      nh=nh+nh_tile
      write(21) hcat
      deallocate(np_halo_all,xv_mean,dxv,iph_halo_all,iph_halo)
      deallocate(xv,hoc,ll,llgp,hcgp,ecgp,hcat)
      deallocate(pid_list,pid_halo)
    enddo ! itile

    if (head) then
      print*,'last np, last id =', np_last, pid_last
      print*,'nh,nhlocal =',nh,nhlocal
    endif

    halo_header%nhalo=nh; sync all
    if (head) then
      do i=2,nn**3
        nh=nh+nh[i]
      enddo
    endif
    sync all; nh=nh[1]; halo_header%nhalo_tot=nh; nhglobal=nh; nfile=nn**3
    rewind(21); write(21) halo_header; close(21)

    rewind(102)
    write(102) nhlocal,nid,nhglobal,nfile
    close(102) ! group_ids
    
    open(101,file=fn_hbt(cur_checkpoint,int(rank,4),11),status='replace',access='stream')
    write(101) nhlocal,nid,nhglobal,nfile,mass_halo(:nhlocal),offset(:nhlocal)
    close(101) ! group_tab

    if (head) then
      print*,'wrote ids,tab file. nhlocal,nid,nhglobal,nfile='
      print*, nhlocal,nid,nhglobal,nfile
      print*,'first and last mass_halo',mass_halo(1),mass_halo(nhlocal)
      print*,'first and last offset',offset(1),offset(nhlocal)
      print*,'halo_header',halo_header
      print*,''
    endif
    deallocate(offset,mass_halo,xp,vp,rhoc,pid)
    sync all
  enddo ! cur_checkpoint

  contains

  subroutine merge_chain(ii,jj)
    ! llgp is a linked list: llgp(ip) means ip->llgp
    ! ip1->ip2->...->ipn
    ! ip1 is the head of chain (hoc)
    ! ipn is the end of chain (eoc)
    integer(8) ii,jj,ihead,jhead,iend,jend,ipart
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
    !ecgp(i)=iend
    hcgp(iend)=jhead ! change hoc
    hcgp(jend)=0 ! set jend as a member
  endsubroutine
  
end
