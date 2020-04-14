program div1d_test
!
! test of reaction rate calculations
! author of code: Egbert Westerhof, DIFFER (c) 2019

   use constants
   use physics_parameters
   use numerics_parameters
   use grid_data
   use plasma_data
   use physics_routines

   implicit none
   integer, parameter :: wp = KIND(1.0D0)
   integer :: input_error, input_error_numerics, input_error_physics
   integer :: restart_error, time_step_error
   integer :: ix
   ! external right_hand_side
   
   real(wp), allocatable :: rate_cx(:), rate_ion(:), rate_exc(:), rate_rec(:), rate_imp(:)

   ! read non-default input parameters
   call read_numerics_parameters(input_error_numerics)
   call read_physics_parameters(input_error_physics)
   input_error = input_error_numerics + 10*input_error_physics
   if( input_error .ne. 0 ) call error_report(input_error, restart_error, time_step_error)
   
   allocate(rate_cx(Nx), rate_ion(Nx), rate_exc(Nx), rate_rec(Nx), rate_imp(Nx))
   allocate(temperature(Nx), density(Nx))
   ! check the atomic rate calculations between minimum_temperature and initial_T
   do ix = 1, Nx
     temperature(ix) = minimum_temperature + (ix-1) * (initial_T - minimum_temperature) / (Nx-1)
     density(ix)     = initial_n
     rate_cx(ix)  = charge_exchange(temperature(ix))
     rate_ion(ix) = ionization(density(ix),temperature(ix))
     rate_exc(ix) = excitation(density(ix),temperature(ix))
     rate_rec(ix) = recombination(density(ix),temperature(ix))
     rate_imp(ix) = impurity_radiation(temperature(ix))
     write( *, '(7(1PE15.6))' ) temperature(ix), density(ix), rate_cx(ix), rate_ion(ix), rate_exc(ix), rate_rec(ix), rate_imp(ix)
   end do

end program div1d_test


subroutine error_report(input_error, restart_error, time_step_error)

   implicit none
   integer, parameter :: wp = KIND(1.0D0)
   integer, intent(in) :: input_error, restart_error, time_step_error
   
   if( input_error .ne. 0 ) then
      write(*,*) 'input_error = ', input_error
      stop 'fatal error on input'
   elseif( restart_error .ne. 0 ) then
      write(*,*) 'restart_error = ', restart_error
      stop 'fatal error on restart'
   elseif( time_step_error .ne. 2 ) then
      write(*,*) 'time_step_error = ', time_step_error
      stop 'fatal error in time step'
   endif

   return
end subroutine error_report

