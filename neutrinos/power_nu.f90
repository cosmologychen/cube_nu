module power_nu
   use cicpower_global
   use parameters
   use ieee_arithmetic
   implicit none
contains

   real function expansion(dt0)
      use variables, only: stime,s2a,ia
      implicit none
      real(8) :: a_x,adot,t_x,tdoa,a8_0
      real(4) :: dt0
      integer i1,i2

      a8_0=sim%a
      i1 = ia
      do while(s2a(i1) <= a8_0 .and. i1 > 1)
         i1 = i1-1
      enddo
      ia = i1

      tdoa = (stime(i1+1)-stime(i1))/(s2a(i1+1)-s2a(i1))
      t_x = stime(i1)+tdoa*(a8_0-s2a(i1))+dt0


      i2 = i1
      do while(stime(i2+1)<t_x .and. i2 > 1)
         i2 = i2-1
      enddo

      adot = (s2a(i2+1)-s2a(i2))/(stime(i2+1)-stime(i2))
      a_x = s2a(i2)+adot*(t_x-stime(i2))

      expansion=a_x-a8_0
   endfunction

   real function dtau_a(da0)
      use variables, only: stime,s2a,s2tau,ia
      implicit none
      real a0,da0,taudot,t_x,tdoa
      integer i1,i2

      a0 = sim%a+da0
      i1 = ia - 1
      do while(s2a(i1) <= a0 .and. i1 > 1)
         i1 = i1-1
      enddo
      i2=i1+1
      tdoa = (stime(i2)-stime(i1))/(s2a(i2)-s2a(i1))
      t_x = stime(i1)+tdoa*(a0-s2a(i1))

      taudot = (s2tau(i2)-s2tau(i1))/(stime(i2)-stime(i1))
      dtau_a = s2tau(i1)+taudot*(t_x-stime(i1))-sim%tau
   endfunction

   real function sf_a(a0)
      use variables, only: stime,s2a,ia,s_fi
      implicit none
      real a0,tdoa
      integer i1,i2

      i1 = 1
      do while(s2a(i1+1)>=a0 .and. i1 < istep_max-1)
         i1 = i1+1
      enddo
      i2=i1+1
      tdoa = (stime(i2)-stime(i1))/(s2a(i2)-s2a(i1))
      sf_a = (stime(i1)+tdoa*(a0-s2a(i1)))*omhsq0-s_fi
      if (sf_a < 0) sf_a = 0

   endfunction

   function get_f_nr(a) result(f_nr_a)
      use omp_lib
      use parameters
      implicit none

      real a, f_nr_a(3)
      real(8) :: x(3),y(3), Fy(3) , dFy(3)
      real(8) ::  k1(3), k2(3), k3(3), k4(3), h
      integer ix,n

      ! print*,'para',Mass_nu,a,sigma_nu,N_eff,k_b,T_gama

      if (a_nu == 0) then
         y = m_nu*a/(sigma_nu*N_eff*k_b*T_gama)
         ! print*,'y',y

         ! 定义步长和循环变量
         h = 1./2
         n = 300./h

         ! 初始化 Fy
         Fy = 0.0

         ! RK4 循环
         do ix=0,n
            x =  h * ix
            k1 = h * (x**2*sqrt(x**2+y**2))/(1+exp(x))
            k3 = h * ((x+h/2)**2*sqrt((x+h/2)**2+y**2))/(1+exp((x+h/2)))
            k4 = h * ((x+h)**2*sqrt((x+h)**2+y**2))/(1+exp((x+h)))
            Fy = Fy + (k1 + 4*k3 + k4)/6
         enddo


         ! 初始化 dFy
         dFy = 0.0

         ! RK4 循环
         do ix=0,n
            x =  h * ix
            k1 = h * (x**2)/((1+exp(x))*sqrt(x**2+y**2))
            k3 = h * ((x+h/2)**2)/((1+exp((x+h/2)))*sqrt((x+h/2)**2+y**2))
            k4 = h * ((x+h)**2)/((1+exp((x+h)))*sqrt((x+h)**2+y**2))
            dFy = dFy + (k1 + 4*k3 + k4)/6
         enddo

         f_nr_a = y**2*dFy/Fy
      elseif (a > a_nu) then
         f_nr_a=1
      else
         f_nr_a=0
      endif
   endfunction

   function interp_quad(a, cur_powerpoint, Pk, z_point) result(spk)
      use variables
      use parameters
      implicit none
      integer, intent(in) ::  cur_powerpoint
      real, intent(in) :: a, Pk(npbin, nmax_redshift), z_point(nmax_redshift)
      real spk(npbin)

      real a1(npbin), a2(npbin), k2(npbin), b2(npbin), c2(npbin)
      real log_a, la(3), lpk(npbin,3)

      la(1) = log(1 / (z_point(cur_powerpoint) + 1))
      la(2) = log(1 / (z_point(cur_powerpoint - 1) + 1))
      la(3) = log(1 / (z_point(cur_powerpoint - 2) + 1))

      lpk(:,1) = log(Pk(:, cur_powerpoint))
      lpk(:,2) = log(Pk(:, cur_powerpoint - 1))
      lpk(:,3) = log(Pk(:, cur_powerpoint - 2))

      a1 = (lpk(:,1) - lpk(:,2)) / (la(1) - la(2))
      a2 = (lpk(:,2) - lpk(:,3)) / (la(2) - la(3))
      k2 = (a1 - a2) / (la(1) - la(3))
      b2 = a1 - k2 * (la(1) + la(2))
      ! 强制曲率为非负：若为负，则退化为线性插值
      where (k2 < 0.0d0 .or. -b2/2/k2 > la(1))
         k2 = 0.0d0
         b2 = (lpk(:,1) - lpk(:,3)) / (la(1) - la(3))
      end where
      c2 = lpk(:,1)
      log_a = log(a)-la(1)
      spk = exp((k2 * log_a ** 2 + b2 * log_a + c2)/2)

      ! if (any((ieee_is_nan(spk(:ifs))))) then
      !    if (head) then
      !       print *, ''
      !       print *, ''
      !       print *, ''
      !       print *, 'Pk interp err ', minval(spk(:ifs)), 1. / a - 1, ifs
      !       print *, 'cur_powerpoint', cur_powerpoint, z_point(cur_powerpoint)
      !       print *, 'Pk0', Pk(:, cur_powerpoint)
      !       print *, 'Pk-1', Pk(:, cur_powerpoint - 1)
      !       print *, 'Pk-2', Pk(:, cur_powerpoint - 2)
      !       print *, 'a1', a1
      !       print *, 'a2', a2
      !       print *, 'k', k2
      !       print *, 'b', b2
      !       print *, 'c', c2
      !       print *, 'result', spk(:ifs)
      !    endif
      !    error stop 'Pk interp err '
      ! endif
   endfunction

   real function interp_tf_F(k_lin_log,k_need)
      use variables
      implicit none
      real, intent(in) :: k_lin_log(npbin),k_need
      real a1,a2,k2,b2,c2
      integer(4) i_mid,i1,i2,i3
      real k_need_log

      k_need_log = log(k_need)

      i1=1; i2=npbin
      do while (i2-i1>1)
         i_mid=(i1+i2)/2
         if (k_need_log>k_lin_log(i_mid)) then
            i1=i_mid
         else
            i2=i_mid
         endif
      enddo

      if (k_need_log-k_lin_log(i1) < k_lin_log(i2)-k_need_log) then
         i3 = i2
         i2 = i1
         i1 = i1-1
      else
         i3 = i2+1
      endif

      if (i3 > npbin) then
         i1 = i1-1
         i2 = i2-1
         i3 = i3-1
      endif

      a1 = (tf_F_log(i1)-tf_F_log(i2))/(k_lin_log(i1)-k_lin_log(i2))
      a2 = (tf_F_log(i2)-tf_F_log(i3))/(k_lin_log(i2)-k_lin_log(i3))
      k2 = (a1-a2)/(k_lin_log(i1)-k_lin_log(i3))
      b2 = a1-k2*(k_lin_log(i1)+k_lin_log(i2))
      c2 = tf_F_log(i1)-k2*k_lin_log(i1)**2-b2*k_lin_log(i1)
      interp_tf_F = exp(k2*k_need_log**2+b2*k_need_log+c2)
   endfunction

   subroutine tf_F_correction
      use variables
      implicit none
      real(8) a1,b1,d1,d2
      integer(4) i1,i2

      i1 = ifs-1
      d1 = (tf_F(i1+1)-tf_F(i1))/(kh_lin(i1+1)-kh_lin(i1))
      do while (i1 > 1 .and. d1 < -1)
         i1 = i1-1
         d1 = (tf_F(i1+1)-tf_F(i1))/(kh_lin(i1+1)-kh_lin(i1))
      enddo

      if (head) print*,'  tf_F_correction',i1,kh_lin(i1),k_fs,tf_F(i1)

      a1 = ((1-f_nu)-tf_F(i1))/(10-kh_lin(i1))
      b1 = tf_F(i1)-kh_lin(i1)*a1

      i2 = i1
      do while (i2 <= npbin .and. kh_lin(i2)<10)
         tf_F(i2) = a1*kh_lin(i2)+b1
         i2 = i2+1
      enddo
      if (i2<npbin) tf_F(i2:)=1-f_nu

   endsubroutine

   function get_sqrt_pk_nu1() result(spk)
      use variables, only: s_f,a_step,sq_Pk_nu_ic,k_mnu_2d,H_i,H_0,f_growth_nu
      use parameters
      implicit none
      real  spk(npbin,3)
      real(8) X(npbin,3)
      real omega_m_zi,omega_l_zi

      X = s_f*k_mnu_2d*T_nu0*C
      spk = ((1+0.0168*X**2+0.0407*X**4)/(1+2.1734*X**2+1.6787*X**4.1811+0.1467*X**8))*sq_Pk_nu_ic*(1+s_f*a_step(0)**2*H_i*f_growth_nu)
   endfunction

   function get_L(i_step) result(L)
      use variables, only: s_f,sf_step,sPk_step,k_mnu_2d
      use parameters
      implicit none
      integer i_step
      real si
      real(16) X(npbin,3),L(npbin,3)
      si = s_f-sf_step(i_step)
      X = si*k_mnu_2d*T_nu0*C
      L = ((1+0.0168*X**2+0.0407*X**4)/(1+2.1734*X**2+1.6787*X**4.1811+0.1467*X**8))*spread(sPk_step(:,i_step), dim=2, ncopies=3)*si
   endfunction

   function get_sqrt_pk_nu2() result(spk)
      use variables, only: tau_step
      implicit none
      real spk(npbin,3)
      real(16) dtau21,dtau10,dtau20
      real(16) spk16(npbin,3),L0(npbin,3),L1(npbin,3),L2(npbin,3)
      real(16) a1(npbin,3), a2(npbin,3), k2(npbin,3), b2(npbin,3)

      integer is

      if (istep == 1) then
         spk = (get_L(1)+get_L(0))*(tau_step(1)-tau_step(0))
         return
      endif

      L1 = get_L(0); L2 = get_L(1); L0 = get_L(2)

      dtau21 = tau_step(2) - tau_step(1); dtau10 = tau_step(1) - tau_step(0); dtau20 = tau_step(2) - tau_step(0)
      a1 = (L0 - L2) / dtau21; a2 = (L2 - L1) / dtau10; k2 = (a1 - a2) / dtau20; b2 = a2 - k2 * dtau10
      spk16 = dtau10**3/3*k2 + dtau10**2/2*b2 + L1*dtau10

      !2是is 1是is-1 0是is-2
      do is = 2,istep
         L0 = L1; L1 = L2; L2 = get_L(is)
         dtau21 = tau_step(is)-tau_step(is-1); dtau10 = tau_step(is-1)-tau_step(is-2); dtau20 = tau_step(is)-tau_step(is-2)
         a1 = (L2 - L1) / dtau21; a2 = (L1 - L0) / dtau10; k2 = (a1 - a2) / dtau20; b2 = a2 - k2 * dtau10
         spk16 =  spk16 + dtau20**3/3*k2 + dtau20**2/2*b2 + L0*dtau20
      enddo
      spk16 =  spk16 +  (k2/3 * (dtau20*(dtau20+dtau10) + dtau10**2)  + b2/2 * (dtau20+dtau10) + L0)*(dtau20-dtau10)

      spk = spk16
   endfunction

   subroutine interp_Pk_CDM
      use omp_lib
      use variables
      use parameters
      implicit none
      integer :: i,interp_powerpoint
      real a_pre

      if (power_step) then
         if (head) print*,''
         if (head) print*,'Power_step',z_powerpoint(sim%cur_powerpoint)
         if (calculate_PK > 0)then
            call system_clock(tp1,tpr)
            call global_power
            Pk_cb_check(:,sim%cur_powerpoint) = xi_cdm(5,1:)[1]
            sync all; call system_clock(tp2,tpr); if (head) print*,'     global_power elapsed time =',real(tp2-tp1)/tpr
            i = istep-1
            interp_powerpoint = merge(sim%cur_powerpoint,3,sim%cur_powerpoint > 2)
            a_pre = 1/(z_powerpoint(sim%cur_powerpoint-1)+1)! merge(1/(z_powerpoint(interp_powerpoint-1)+1),1/(z_powerpoint(1)+1),sim%cur_powerpoint == 2)
            ! if (head) print*,' ',a_pre,interp_powerpoint
            do while (a_step(i) >= a_pre .and. i >= 0)
               ! if (head) print*,'      interp',1/a_step(i)-1,i
               sPk_step(:,i) =  interp_quad(a_step(i),interp_powerpoint,Pk_cb_check,z_powerpoint)
               i = i-1
            enddo
            ! stop
         endif
         sim%cur_powerpoint=sim%cur_powerpoint+1
         power_step = .false.
      endif

      interp_powerpoint = merge(sim%cur_powerpoint-1,3,sim%cur_powerpoint > 3)
      sPk_step(:,istep) = interp_quad(a_step(istep),interp_powerpoint,Pk_cb_check,z_powerpoint)
      if (calculate_PK == -1) then
         Pk_nus(:,1) = interp_quad(a_step(istep),interp_powerpoint,Pk_nu_check,z_powerpoint)
         Pk_nus(:,2) = interp_quad(a_step(istep),interp_powerpoint,Pk_nu_check,z_powerpoint)
         Pk_nus(:,3) = interp_quad(a_step(istep),interp_powerpoint,Pk_nu_check,z_powerpoint)
      endif
      if (head) print*,'      interp',1/a_step(istep)-1,sim%cur_powerpoint,istep
   endsubroutine

   subroutine get_tf_cb2matter
      use variables
      use parameters
      implicit none

      integer is,n
      real s_pre,s_post
      real sqrt_pk_nu1(npbin,3),sqrt_pk_nu2(npbin,3)
      real(8) Lx0(npbin,3),Lx1(npbin,3),Lx2(npbin,3)

      k_fs = min((0.8*sqrt(a_step(istep)*sim%omega_m/0.3) * maxval(m_nu) * h0 )*tf_smooth,10.)
      ifs = npbin
      do while (kh_lin(ifs) > k_fs)
         ifs = ifs-1
      enddo
      f_nr = get_f_nr(a_step(istep))
      s_f  = sf_a(a_step(istep))
      sf_step(istep) = s_f

      if(head) then
         print*,''
         print*,'Get Tf_nu'
         print*,'  save Pk in z = ',1/a_step(istep)-1,istep
         print*,'  s_f    = ',s_f
         print*,'  1-f_nu = ',(1-f_nu)
         print*,'  f_nus  = ',f_nus
         print*,'  f_nr   = ',f_nr
         print*,'  ifs    = ',ifs,npbin
         print*,'  k_fs   = ',k_fs,kh_lin(ifs)
      endif
      ! print*,'Get Tf_nu',image


      call tic(96)
      call interp_Pk_CDM
      call toc(96)


      call tic(97)
      if (istep > 0) then
         if(head .and. calculate_PK > -1) then
            sqrt_pk_nu1 = get_sqrt_pk_nu1()
            sqrt_pk_nu2 = 0.75*sim%omega_m*H_0**2*get_sqrt_pk_nu2() ! 0.75 = (3/2)/2
            print*,'  sqrt_pk_nu1',sqrt_pk_nu1(1,1),sqrt_pk_nu1(npbin/3,1),sqrt_pk_nu1(npbin/3*2,1),sqrt_pk_nu1(npbin,1)
            print*,'  sqrt_pk_nu2',sqrt_pk_nu2(1,1),sqrt_pk_nu2(npbin/3,1),sqrt_pk_nu2(npbin/3*2,1),sqrt_pk_nu2(npbin,1)
            Pk_nus = (sqrt_pk_nu1+sqrt_pk_nu2)
         endif
      else if (istep == 0) then
         Pk_nus=sq_Pk_nu_ic
      else
         ! 不需要更新tf
         sPk_step(:,istep) = tf_F**2*sPk_step(:,istep) !sPk_step_cdm --> sPk_step_matter
         return
      endif

      Pk_nu = (f_nus(1)*f_nr(1)*Pk_nus(:,1)+f_nus(2)*f_nr(2)*Pk_nus(:,2)+f_nus(3)*f_nr(3)*Pk_nus(:,3))**2
      Pk_nus = Pk_nus**2
      tf_F = ((1-f_nu)*sPk_step(:,istep)+f_nu*sqrt(abs(Pk_nu)))/sPk_step(:,istep) ! sPk is sqrt(Pk)
      if (calculate_PK > -1) call tf_F_correction
      call toc(97)

      if(head) then
         write(str_z,'(f8.4)') 1/a_step(istep)-1
         print*,'  Pk_step',sPk_step(1,istep)**2,sPk_step(npbin/3,istep)**2,sPk_step(npbin/3*2,istep)**2,sPk_step(npbin,istep)**2
         print*,'  Pk_nu',Pk_nu(1),Pk_nu(npbin/3),Pk_nu(npbin/3*2),Pk_nu(npbin)
         if (m_nu(1)>0) print*,'  Pk_nu1',Pk_nus(1,1),Pk_nus(npbin/3,1),Pk_nus(npbin/3*2,1),Pk_nus(npbin,1)
         if (m_nu(2)>0) print*,'  Pk_nu2',Pk_nus(1,2),Pk_nus(npbin/3,2),Pk_nus(npbin/3*2,2),Pk_nus(npbin,2)
         if (m_nu(3)>0) print*,'  Pk_nu3',Pk_nus(1,3),Pk_nus(npbin/3,3),Pk_nus(npbin/3*2,3),Pk_nus(npbin,3)
         print*,'  tf_F',tf_F(1),tf_F(npbin/3),tf_F(npbin/3*2),tf_F(npbin)
         open(211,file=nupath//'a_step.txt',status='replace',access='stream'); write(211) a_step(:istep); close(211)
         open(311,file=nupath//'Pk_nu/Pk_nu_'//trim(adjustl(str_z))//'.txt',status='replace',access='stream'); write(311) Pk_nu; close(311)
         open(312,file=nupath//'Pk_nus/Pk_nus_'//trim(adjustl(str_z))//'.txt',status='replace',access='stream'); write(312) Pk_nus; close(312)
         open(411,file=nupath//'Pk_m/Pk_m_'//trim(adjustl(str_z))//'.txt',status='replace',access='stream'); write(411) sPk_step(:,istep)**2; close(411)
         open(511,file=nupath//'tf/Tf_nu_'//trim(adjustl(str_z))//'.txt',status='replace',access='stream'); write(511) tf_F; close(511)
         print*,'  path',nupath//'*/*_'//trim(adjustl(str_z))//'.txt'
         sPk_step(:,istep) = tf_F*sPk_step(:,istep) !sPk_step_cdm --> sPk_step_matter
      endif

      call tic(98)
      tf_F_log = log(tf_F(:)[1])
      call nu_correction_coarse()
      pm%nwork = ngt   ;call nu_correction(tf2   ,Gk2   ,box/nn/nnt    )
      pm%nwork = nft(2);call nu_correction(tf3_2 ,Gk3_2 ,box/nn/nnt/nns)
      pm%nwork = nft(3);call nu_correction(tf3_4 ,Gk3_4 ,box/nn/nnt/nns)
      pm%nwork = nft(4);call nu_correction(tf3_6 ,Gk3_6 ,box/nn/nnt/nns)
      pm%nwork = nft(5);call nu_correction(tf3_8 ,Gk3_8 ,box/nn/nnt/nns)
      call toc(98)
   endsubroutine

   subroutine nu_correction(tfk,Gk,size)
      use omp_lib
      use variables
      use parameters
      implicit none

      integer i,j,k
      real kr,kx(3),size
      real tfk(pm%nwork/2+1,pm%nwork,pm%nwork)
      ! real k_tf(pm%nwork/2+1,pm%nwork,pm%nwork)
      real Gk(pm%nwork/2+1,pm%nwork,pm%nwork)

      tfk = 1
      ! k_tf =0
      if (Mass_nu > 0) then
         !$omp paralleldo default(shared) schedule(dynamic)&
         !$omp& private(i,j,k,kx,kr)
         do k=1,pm%nwork
            do j=1,pm%nwork
               do i=1,pm%nwork/2+1
                  kx=(mod([i,j,k]+pm%nwork/2,pm%nwork)-pm%nwork/2-1)*(2*pi)/size
                  if (i==1 .and. j==1 .and. k==1) cycle ! zero frequency
                  kr=sqrt(kx(1)**2+kx(2)**2+kx(3)**2)
                  tfk(i,j,k) = interp_tf_F(kh_lin_log,kr)*Gk(i,j,k)
               enddo
            enddo
         enddo
         !$omp endparalleldo
      endif
      ! call toc_nu
   endsubroutine

   subroutine nu_correction_coarse()
      use omp_lib
      use variables
      use parameters
      use pencil_fft
      implicit none

      integer i,j,k,ig,jg,kg
      real kr,kx(3)

      if (Mass_nu > 0) then
         !$omp paralleldo default(shared) schedule(dynamic)&
         !$omp& private(kg,jg,ig,i,j,k,kx,kr)
         do k=1,npen
            do j=1,nw
               do i=1,nyquist+1
                  kg=(nn*(icz-1)+icy-1)*npen+k
                  jg=(icx-1)*nw+j
                  ig=i
                  kx=(mod([ig,jg,kg]+nyquist-1,nw_global)-nyquist)*(2*pi)/box
                  if (ig==1.and.jg==1.and.kg==1) cycle ! zero frequency
                  kr=sqrt(kx(1)**2+kx(2)**2+kx(3)**2)
                  tf1(i,j,k) = interp_tf_F(kh_lin_log,kr)*Gk1(i,j,k)
                  ! k_tf1(i,j,k) = kr
               enddo
            enddo
         enddo
         !$omp endparalleldo
      endif

      ! call toc_nu
   endsubroutine

endmodule
