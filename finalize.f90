subroutine finalize
  use variables
  use cubefft
  use pencil_fft
  use power_nu
  implicit none
  save
  include 'fftw3.f'
  istep = istep+1
  a_step(istep)=1
  call destroy_cubefft_plan
  call destroy_penfft_plan
  
  if (Mass_nu > 0) then
    call get_tf_cb2matter
    deallocate(Gk1,Gk2,Gk3_2,Gk3_4,Gk3_6,Gk3_8)!,Gk3_12)!,Gk3_16)
  endif
  deallocate(xp,xp_new,vp,vp_new,rhoc)
#ifdef PID
  deallocate(pid,pid_new)
#endif
  if (head) then
   open(101,file=output_name('tcat'),status='replace',access='stream')
   write(101) int(sim%timestep,kind=4), tcat(:,:sim%timestep)
   close(101)
  endif
endsubroutine
