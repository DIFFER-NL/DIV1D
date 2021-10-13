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
   use experiments

   implicit none
   integer, parameter :: wp = KIND(1.0D0)
   real(wp) :: start_time, end_time
   integer :: input_error, input_error_numerics, input_error_physics
   integer :: restart_error, time_step_error
   integer :: istep, ix, itol, iopt, ml, mu, lrw, liw, nzswag_input
   ! external right_hand_side
   
   ! the following is needed for dvode_f90
   integer :: itask, istate
   integer :: attempt
   integer,  allocatable :: iwork(:)
   integer,  allocatable :: bounded_components(:)
   real(wp), allocatable :: lower_bounds(:), upper_bounds(:)
   real(wp), allocatable :: abstol_vector(:), rwork(:)
   type (VODE_OPTS) :: options

   ! time the program execution time exluding external system calls (ifortran, gfortran) incl external system calls (solaris f90)
   real(wp) :: T1,T2
   call cpu_time(T1)

   ! read non-default input parameters
   call read_numerics_parameters(input_error_numerics)
   call read_physics_parameters(input_error_physics)
   input_error = input_error_numerics + 10*input_error_physics
   if( input_error .ne. 0 ) call error_report(input_error, restart_error, time_step_error)

   ! initialize the grid
   call initialize_grid

   ! set the initial value of the solution vector
   if( simple_sol ) then 
      call init_simple_sol
   else
      call initial_values
   endif

   start_time = 0.0
   restart_error = 0
   if( restart ) call read_restart_file( restart_error )
   if( restart .and. restart_error .ne. 0 ) call error_report(input_error, restart_error, time_step_error)
   
   ! initialize the gas_puff
   call initialize_gas_puff(Nx)

   ! write the inital solution to file
   ! calculate the fluxes
   call calculate_fluxes( Nx, start_time,  density, velocity, temperature, neutral, Gamma_n, Gamma_mom, q_parallel, neutral_flux )
   ! calculate the sources
   call calculate_sources( Nx, start_time,  density, velocity, temperature, neutral, q_parallel, &
                          Source_n, Source_v, Source_Q, source_neutral )
   open( UNIT=10, FILE='div1d_output.txt' )
   call write_header
   call write_solution( start_time )
  
   ! setting the absolute error tolerance for dvode_f90 or dlsode
   ! note that the solution vector y is normalized so the absolute error tolarance in all elements of y can be made equal
   allocate( abstol_vector(4*Nx) )
   abstol_vector(1:4*Nx) = abstol
   ! abstol_vector(1:Nx) = initial_n * abstol
   ! abstol_vector(Nx+1:2*Nx) = 1.0
   ! abstol_vector(2*Nx+1:3*Nx) = initial_T * abstol
   ! abstol_vector(3*Nx+1:4*Nx) = initial_n * abstol
   
   ! enforce nonnegative densities and energies by lower bounds in dvode_f90
   allocate( bounded_components(Nx), lower_bounds(Nx), upper_bounds(Nx) )
   ! plasma density / only energy for now
   do ix = 1, Nx
      bounded_components(ix) = ix + 2*Nx
   enddo
   ! energy
   ! bounded_components(  Nx+1:2*Nx) = bounded_components(1:Nx) + 2 * Nx
   ! neutral density
   ! bounded_components(2*Nx+1:3*Nx) = bounded_components(1:Nx) + 3 * Nx
   ! set lower bounds to zero
   lower_bounds = 0.0d+0
   ! set artifically high upper bounds
   upper_bounds = 1.0d+50

   ! setting the options fo dvode_f90
   ! rough estimate of maximum number of nonzeros in Jacobian
   nzswag_input = 48 * Nx
   if( method .gt. 0 ) then 
      options = set_opts(RELERR=reltol, ABSERR_VECTOR=abstol_vector, &
                         CONSTRAINED=bounded_components, CLOWER=lower_bounds, CUPPER=upper_bounds,  &
                         METHOD_FLAG=method, MXSTEP=max_step, NZSWAG=nzswag_input, MA28_ELBOW_ROOM=200, MA28_RPS=.TRUE.)
      ! ! call set_jacobian_sparsity_structure
      ! call set_diagonal_jacobian
      ! call USERSETS_IAJA(IAUSER,NIAUSER,JAUSER,NJAUSER)

   elseif( method .lt. 0 ) then
      ! allocate arrays needed by dlsode
      ! RWORK :WORK   Real work array of length at least:
      !        20 + 16*NEQ                    for MF = 10,
      !        22 +  9*NEQ + NEQ**2           for MF = 21 or 22,
      !        22 + 10*NEQ + (2*ML + MU)*NEQ  for MF = 24 or 25.
      ! LRW   :IN     Declared length of RWORK (in user's DIMENSION statement).
      ! IWORK :WORK   Integer work array of length at least:
      !        20        for MF = 10,
      !        20 + NEQ  for MF = 21, 22, 24, or 25.
      ! If MF = 24 or 25, input in IWORK(1),IWORK(2) the lower and upper Jacobian half-bandwidths ML,MU.
      ! LIW   :IN     Declared length of IWORK (in user's DIMENSION statement).
      ml = 3*Nx
      mu = 3*Nx
      lrw = 22 + (20+2*ml+mu)*4*Nx
      liw = 20 + 4*Nx
      allocate( rwork(lrw), iwork(liw) )
      iwork(1) = ml
      iwork(2) = mu
      iwork(6) = max_step ! maximum number of internal iterations
   endif

   istate = 1
   do istep=1, ntime
      write(*,*) 'start integration step nr.', istep
      end_time = start_time+delta_t
      if( method .gt. 0 ) then
         ! we use dvode_f90 for the integration
         ! CALL DVODE_F90(F,NEQ,Y,T,TOUT,ITASK,ISTATE,OPTIONS,J_FCN=JAC,G_FCN=GEX)
         itask = 1 ! for normal computation in dvode_f90 till end_time
         call dvode_f90( right_hand_side, 4*Nx, y, start_time, end_time, itask, istate, options )
      elseif( method .lt. 0 ) then
         ! we use dlsode.f for the integration
         itask = 1 ! for normal computation in dvode_f90 till end_time
         itol = 2 ! the absolute error tolerance is specified in an array
         iopt = 1 ! optional inputs (upper and lower half bandwidths of the Jacobian)
         ! CALL DLSODE(F,               NEQ, Y,          T,     TOUT, ITOL,   RTOL,          ATOL, ITASK, ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, MF)
         call dlsode( right_hand_side, 4*Nx, y, start_time, end_time, itol, reltol, abstol_vector, itask, istate, iopt, rwork, lrw, iwork, liw, jac, abs(method) )
      else
         call rk4( right_hand_side, 4*Nx, y, start_time, end_time )
      endif
      ! make sure that the solution respects the minimum density and temperature
      call y2nvt( Nx, y, density, velocity, temperature, neutral )
      do ix = 1, Nx
         if( temperature(ix) .lt. minimum_temperature) temperature(ix) = minimum_temperature
         ! if( density(ix)     .lt. minimum_density)     density(ix) =     minimum_density
         ! if( neutral(ix)     .lt. minimum_density)     neutral(ix) =     minimum_density
      enddo
      call nvt2y( Nx, density, velocity, temperature, neutral, y )
      time_step_error = istate
      if( istate .ne. 2 .and. method .gt. 0 ) then
         attempt = 0
         do while( istate .ne. 2 .and. attempt .le. max_attempts )
            ! first try calling dvode again with istate = 1
            attempt = attempt + 1
            write(*,*) 'trying to continue integration. attempt nr.', attempt
            istate = 1
            call dvode_f90( right_hand_side, 4*Nx, y, start_time, end_time, itask, istate, options )
         enddo
         if( istate .ne. 2 ) then 
            call error_report(input_error, restart_error, time_step_error)
         else
            write(*,*) 'success on ', attempt, 'th attempt. continuing'
         endif
      endif
      if( mod( istep, nout ) .eq. 0 ) then
         ! call y2nvt( Nx, y, density, velocity, temperature, neutral )
         ! calculate the fluxes
         call calculate_fluxes( Nx, start_time, density, velocity, temperature, neutral, Gamma_n, Gamma_mom, q_parallel, neutral_flux )
         ! calculate the sources
         call calculate_sources( Nx, start_time, density, velocity, temperature, neutral, q_parallel, Source_n, Source_v, Source_Q, source_neutral )
         call write_solution( end_time )
      endif
      start_time = end_time
      ! reset istate to 1 every so often
      if( modulo( istep, istate_mod ) .eq. 0 ) istate = 1
   enddo

   call cpu_time(T2)


   write( 10, * ) '   cpu_time   = ', T2-T1

   call write_restart_file
   stop 'div1d completed'
end program div1d


subroutine write_header

   use numerics_parameters
   use physics_parameters
   use grid_data, only : x, B_field

   implicit none
   integer :: i
   integer, parameter :: wp = KIND(1.0D0)

   ! note these lists should still be completed (GD, complete now?)
   write( 10, * ) 'numerics parameters:'
   write( 10, * ) '   Nx         = ', Nx
   write( 10, * ) '   ntime      = ', ntime
   write( 10, * ) '   nout       = ', nout
   write( 10, * ) '   dxmin      = ', dxmin
   write( 10, * ) '   delta_t    = ', delta_t  
   write( 10, * ) '   abstol     = ', abstol  
   write( 10, * ) '   reltol     = ', reltol  
   write( 10, * ) '   viscocity  = ', viscosity
   write( 10, * ) '   method     = ', method  
   write( 10, * ) '   restart    = ', restart
   write( 10, * ) 'physics parameters:'
   write( 10, * ) '   gamma      = ', gamma
   write( 10, * ) '   L          = ', L
   write( 10, * ) '   sintheta   = ', sintheta
   write( 10, * ) '   mass       = ', mass
   write( 10, * ) '   Gamma_X    = ', Gamma_X
   write( 10, * ) '   q_parX     = ', q_parX
   write( 10, * ) '   flux_exp   = ', flux_expansion
   write( 10, * ) '   initial_n  = ', initial_n
   write( 10, * ) '   initial_v  = ', initial_v
   write( 10, * ) '   initial_T  = ', initial_T
   write( 10, * ) '   initial_a  = ', initial_a
   write( 10, * ) '   density_ramp_rate       = ', density_ramp_rate
   write( 10, * ) '   energy_loss_ion         = ', energy_loss_ion
   write( 10, * ) '   neutral_residence_time  = ', neutral_residence_time
   write( 10, * ) '   redistributed_fraction  = ', redistributed_fraction
   write( 10, * ) '   recycling               = ', recycling
   write( 10, * ) '   carbon_concentration    = ', carbon_concentration
   write( 10, * ) '   gas_puff_source         = ', gas_puff_source
   write( 10, * ) '   gas_puff_location       = ', gas_puff_location
   write( 10, * ) '   gas_puff_width          = ', gas_puff_width
   write( 10, * ) '   elm_start_time          = ', elm_start_time
   write( 10, * ) '   elm_ramp_time           = ', elm_ramp_time
   write( 10, * ) '   elm_time_between        = ', elm_time_between
   write( 10, * ) '   elm_expelled_heat       = ', elm_expelled_heat
   write( 10, * ) '   elm_expelled_particles  = ', elm_expelled_particles
   write( 10, * ) '   switch_elm_series       = ', switch_elm_series
   write( 10, * ) '   gaussian_elm            = ', gaussian_elm
   write( 10, * ) '   radial_loss_factor      = ', radial_loss_factor
   write( 10, * ) '   radial_loss_gaussian    = ', radial_loss_gaussian
   write( 10, * ) '   radial_loss_width       = ', radial_loss_width
   write( 10, * ) '   radial_loss_location    = ', radial_loss_location
   write( 10, * ) '   switch_dyn_nu           = ', switch_dyn_nu 
   write( 10, * ) '   switch_dyn_gas          = ', switch_dyn_gas 
   write( 10, * ) '   switch_dyn_rec          = ', switch_dyn_rec
   write( 10, * ) '   switch_dyn_rad_los      = ', switch_dyn_rad_los
   write( 10, * ) '   switch_car_con_prf      = ', switch_car_con_prf
   write( 10, * ) '   switch_dyn_qpar         = ', switch_dyn_qpar
   write( 10, * ) '   switch_dyn_red_frc      = ', switch_dyn_red_frc

   write( 10, '(A195)' ) ' X [m]   car_con_prf [%]    gas_puff_prf []  B_field [fraction]'  !       rad_los_prf  '
   write( 10, '(13(1PE15.6))' ) ( x(i),car_con_prf(i), gas_puff(i), B_field(i), i=1,Nx )
  !write( 10, '(13(1PE15.6))' ) ( x(i),car_con_prf(i), gas_puff(i),rad_los_prf, i=1,Nx )
   return
end subroutine write_header


subroutine write_solution( time )
   use physics_parameters, only : dyn_gas, dyn_nu, dyn_rec, dyn_rad_los, dyn_qparX, dyn_red_frc
   use numerics_parameters, only : Nx, delta_t
   use grid_data, only : x, B_field
   use plasma_data, only : density, velocity, temperature, neutral, Gamma_n, Gamma_mom, q_parallel, neutral_flux, Source_n, Source_v, Source_Q, source_neutral

   implicit none
   integer, parameter :: wp = KIND(1.0D0)
   integer :: i
   real(wp), intent(in) :: time
   integer :: itime 
   itime = time / delta_t  

   write( 10, * ) 'time        = ', time
   write( 10, * ) 'dyn_gas     = ', dyn_gas(itime)
   write( 10, * ) 'dyn_nu      = ', dyn_nu(itime)
   write( 10, * ) 'dyn_rec     = ', dyn_rec(itime)
   write( 10, * ) 'dyn_rad_los = ', dyn_rad_los(itime)
   write( 10, * ) 'dyn_qparX   = ', dyn_qparX(itime)
   write( 10, * ) 'dyn_red_frc = ', dyn_red_frc(itime)
   write( 10, '(A195)' ) '    X [m]        N [/m^3]       V [m/s]         T [eV]        Nn [/m^3]      Gamma_n    Gamma_mom [Pa]      q_parallel    neutral_flux     Source_n       Source_v       Source_Q     source_neut  '
   write( 10, '(13(1PE15.6))' ) ( x(i), density(i), velocity(i), temperature(i), neutral(i), &                 !!!! the neutral flux should not be multiplied with the B_field !!!!    
   &                                   Gamma_n(i)*B_field(i),Gamma_mom(i)*B_field(i),q_parallel(i)*B_field(i),neutral_flux(i)*B_field(i), & ! multiplied by B_field because the code uses normalized values
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

