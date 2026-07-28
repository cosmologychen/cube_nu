! #define new_expansion
! #define nu_expansion

subroutine timestep
   use variables
   use power_nu
   implicit none
   save
   integer ntry,j
   real ra,da_1,da_2,a_next,z_next,ai

   real(16) ::  k1, k4

   dt_old=dt
   sim%timestep=sim%timestep+1
   call tic(1)
   if (head) then
      print*, ''
      print*, '-------------------------------------------------------'
      print*, 'timestep    :',sim%timestep,istep
      dt_e=dt_max
      ntry=0
      do
         ntry=ntry+1
         da = expansion(dt_e)
         ra=da/(sim%a+da)
         if (ra>ra_max) then
            dt_e=dt_e*(ra_max/ra)
         else
            exit
         endif
         if (ntry>10) exit
      enddo
      dt=min(dt_e,sim%dt_pm1,sim%dt_pm2,sim%dt_pm3,dt_refine*sim%dt_pp,sim%dt_vmax)
      da = expansion(dt)
      ! check if checkpointing is needed
      checkpoint_step=.false.
      halofind_step=.false.
      power_step = .false.
#   ifdef HALOFIND
      z_next=max(z_checkpoint(sim%cur_checkpoint),z_halofind(sim%cur_halofind))
#   else
      z_next=z_checkpoint(sim%cur_checkpoint)
#   endif
      if (Mass_nu > 0.0) z_next = max(z_next,z_powerpoint(sim%cur_powerpoint))
      a_next=1.0/(1+z_next)
      if (da>=a_next-sim%a) then !改
         ! print*, 'z         :',z_next,1/a_nu-1
         if (z_next==z_checkpoint(sim%cur_checkpoint)) then
            checkpoint_step=.true.
            if (sim%cur_checkpoint==n_checkpoint) final_step=.true.
         endif
         if (z_next==z_powerpoint(sim%cur_powerpoint)) power_step=.true.
#   ifdef HALOFIND
         if (z_next==z_halofind(sim%cur_halofind)) halofind_step=.true.
#   endif
         do while (abs((sim%a+da)/a_next-1)>=1e-6 .or. (sim%a+da) > 1)
            dt=dt*(a_next-sim%a)/da
            da = expansion(dt)
         enddo
      endif

      ra=da/(sim%a+da)
      a_mid=sim%a+(da/2)

      tcat(41,istep)=sim%a
      tcat(42,istep)=a_mid
      tcat(43,istep)=sim%a+da

      !nu
      dtau = dtau_a(da)

      print*, 'tau       :',sim%tau,sim%tau+dtau
      print*, 'z         :',1.0/sim%a-1.0,1.0/(sim%a+da)-1.0
      print*, 'z_next    :',z_next,z_powerpoint(sim%cur_powerpoint),z_checkpoint(sim%cur_checkpoint)
      print*, 'a         :',sim%a,a_mid,sim%a+da
      print*, 'expansion :',ra
      print*, 'dt        :',dt
      print*, 'dt_e      :',dt_e
      print*, 'dt_pm1    :',sim%dt_pm1
      print*, 'dt_pm2    :',sim%dt_pm2
      print*, 'dt_pm3    :',sim%dt_pm3
      print*, 'dt_pp     :',sim%dt_pp
      print*, 'dt_vmax   :',sim%dt_vmax
      print*, 'cur_ppt   :',sim%cur_powerpoint,z_powerpoint(sim%cur_powerpoint)
      print*, ''
      sim%tau=sim%tau+dtau
      sim%t=sim%t+dt
      sim%a=sim%a+da
      tau_step(istep+1)=sim%tau
      a_step(istep+1)=sim%a
   endif
   sync all

   !a=a[1]
   sim%tau=sim[1]%tau
   sim%t=sim[1]%t
   sim%a=sim[1]%a
   a_mid=a_mid[1]
   dt=dt[1]
   checkpoint_step=checkpoint_step[1]
   power_step=power_step[1]
   tau_step=tau_step(:)[1]
   a_step=a_step(:)[1]
#ifdef HALOFIND
   halofind_step=halofind_step[1]
#endif
   final_step=final_step[1]
   call toc(1)
endsubroutine timestep

