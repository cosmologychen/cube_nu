! #define power_test
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   CUBE™ in Coarray Fortran  !
!   haoran@xmu.edu.cn         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program CUBE
   use omp_lib
   use variables
   use power_nu
   implicit none
   save

   sim%cur_checkpoint=2 !set checkpoint number
   call initialize
   call particle_initialization
   call buffer_grid
   call buffer_x
   call buffer_v
   if (sim%cur_powerpoint == 1 .and. Mass_nu >0) then
      power_step = .true.
      call interp_Pk_CDM
      sim%cur_powerpoint = 2
   endif
   sim%cur_checkpoint = sim%cur_checkpoint+1
   sim%cur_halofind = sim%cur_halofind+1

   if (head) open(77,file=output_dir()//'vinfo'//output_suffix(),access='stream',status='replace')

   if (head) print*, '---------- starting main loop ----------'
   do istep=sim%timestep,istep_max
      call system_clock(t_start,t_rate)
      call tic(100)
      call timestep
      call drift
      call buffer_grid
      call buffer_x
      call kick
      call buffer_v
      if (checkpoint_step .or. halofind_step) then
         dt_old=0
         call drift
         if (checkpoint_step) then
            call checkpoint
            sim%cur_checkpoint = sim%cur_checkpoint+1
         endif
         call buffer_grid
         call buffer_x
         call buffer_v
         ! if (halofind_step) then
         !    call halofind_FoF
         !    sim%cur_halofind = sim%cur_halofind+1
         ! endif
         call print_header(sim)
         if (final_step) exit
      endif
      call system_clock(t_end,t_rate)
      call toc(100)
      if(head) print*, 'total elapsed time =',tcat(100,istep),real(t_end-t_start)/t_rate,'secs';
   enddo
   if (head) close(77)
   call finalize

contains

   function Dgrow(scale_factor)
      implicit none
      real, parameter :: om=omega_m
      real, parameter :: ol=omega_l
      real scale_factor
      real Dgrow
      real g,ga,hsq,oma,ola
      hsq=om/scale_factor**3+(1-om-ol)/scale_factor**2+ol
      oma=om/(scale_factor**3*hsq)
      ola=ol/hsq
      g=2.5*om/(om**(4./7)-ol+(1+om/2)*(1+ol/70))
      ga=2.5*oma/(oma**(4./7)-ola+(1+oma/2)*(1+ola/70))
      Dgrow=scale_factor*ga/g
   endfunction Dgrow

endprogram

