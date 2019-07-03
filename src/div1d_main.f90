program div1d
!
! program to solve the time dependant evolution on a 1D fluxtube in a divertor
! equations from Nakazawa et al. Plasma Phys. Control. Fusion 42 (2000) 401
! author of code: Egbert Westerhof, DIFFER (c) 2019

   use constants
   use physics_parameters
   use numerics_parameters
   use grid_data
   use plasma_data
   use physics_routines
   use dvode_f90_m

   implicit none
   integer, parameter :: wp = KIND(1.0D0)
   real(wp) :: start_time, end_time
   integer :: input_error, input_error_numerics, input_error_physics
   integer :: restart_error, time_step_error
   integer :: istep
   ! external right_hand_side
   
   ! the following is needed for dvode_f90
   integer :: itask, istate
   real(wp), allocatable :: abstol_vector(:)
   type (VODE_OPTS) :: options

   ! read non-default input parameters
   call read_numerics_parameters(input_error_numerics)
   call read_physics_parameters(input_error_physics)
   input_error = input_error_numerics + 10*input_error_physics
   if( input_error .ne. 0 ) call error_report(input_error, restart_error, time_step_error)

   ! initialize the grid
   call initialize_grid

   ! set the initial value of the solution vector
   call initial_values

   start_time = 0.0
   restart_error = 0
   if( restart ) call read_restart_file( restart_error )
   if( restart .and. restart_error .ne. 0 ) call error_report(input_error, restart_error, time_step_error)

   ! write the inital solution to file
   ! calculate the fluxes
   call calculate_fluxes( Nx, density, velocity, temperature, neutral, Gamma_n, pressure, q_parallel, neutral_flux )
   ! calculate the sources
   call calculate_sources( Nx, density, velocity, temperature, neutral, Source_n, Source_v, Source_Q, source_neutral )
   open( UNIT=10, FILE='div1d_output.txt' )
   call write_solution( start_time )

   ! setting the options for dvode_f90
   allocate( abstol_vector(4*Nx) )
   abstol_vector(1:Nx) = initial_n * abstol
   abstol_vector(Nx+1:2*Nx) = 1.0
   abstol_vector(2*Nx+1:3*Nx) = initial_T * abstol
   abstol_vector(3*Nx+1:4*Nx) = initial_n * abstol
   if( method .gt. 0 ) options = set_opts(RELERR=reltol, ABSERR_VECTOR=abstol_vector, METHOD_FLAG=method, MXSTEP=100000, NZSWAG=20000, MA28_ELBOW_ROOM=10)

   istate = 1
   do istep=1, ntime
      end_time = start_time+delta_t
      ! we use dvode_f90 for the integration
      ! CALL DVODE_F90(F,NEQ,Y,T,TOUT,ITASK,ISTATE,OPTIONS,J_FCN=JAC,G_FCN=GEX)
      itask = 1 ! for normal computation in dvode_f90 till end_time
      if( method .gt. 0 ) then
         call dvode_f90( right_hand_side, 4*Nx, y, start_time, end_time, itask, istate, options )
      else
         call rk4( right_hand_side, 4*Nx, y, start_time, end_time )
      endif
      ! make sure that the solution respects the minimum density and temperature
      call y2nvt( Nx, y, density, velocity, temperature, neutral )
      call nvt2y( Nx, density, velocity, temperature, neutral, y )
      time_step_error = istate
      if( istate .ne. 2 .and. method .gt. 0 ) call error_report(input_error, restart_error, time_step_error)
      if( mod( istep, nout ) .eq. 0 ) then
         ! call y2nvt( Nx, y, density, velocity, temperature, neutral )
         ! calculate the fluxes
         call calculate_fluxes( Nx, density, velocity, temperature, neutral, Gamma_n, pressure, q_parallel, neutral_flux )
         ! calculate the sources
         call calculate_sources( Nx, density, velocity, temperature, neutral, Source_n, Source_v, Source_Q, source_neutral )
         call write_solution( end_time )
      endif
      start_time = end_time
   enddo

   call write_restart_file

   stop 'div1d completed'
end program div1d


subroutine write_solution( time )

   use grid_data, only : Nx, x
   use plasma_data, only : density, velocity, temperature, neutral, Gamma_n, pressure, q_parallel, neutral_flux, Source_n, Source_v, Source_Q, source_neutral

   implicit none
   integer, parameter :: wp = KIND(1.0D0)
   integer :: i
   real(wp), intent(in) :: time
   
   write( 10, * ) 'time = ', time
   write( 10, '(A156)' ) '    X [m]     N [/m^3]    V [m/s]      T [eV]     Nn [/m^3]   Gamma_n      P [Pa]    q_parallel neutral_flux  Source_n    Source_v    Source_Q  source_neut'
   write( 10, '(13(1PE15.6))' ) ( x(i), density(i), velocity(i), temperature(i), neutral(i), &
   &                                   Gamma_n(i), pressure(i), q_parallel(i), neutral_flux(i), &
   &                                   Source_n(i), Source_v(i), Source_Q(i), source_neutral(i), i=1,Nx )
   
   return
end subroutine write_solution


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

