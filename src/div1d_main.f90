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

   ! the following is used to format for Matlab integration
   real(wp), allocatable :: float_div1d(20) = 0
   integer,  allocatable :: int_div1d(20) = 0
   logical               :: log_phys = .true.
   integer, :: call_from_matlab = 0

   ! time the program execution time exluding external system calls (ifortran, gfortran) incl external system calls (solaris f90)
   real(wp) :: T1,T2
   call cpu_time(T1)

   ! Initialization from FORTRAN will look for .txt and .dat files
   call initialize_div1d(call_from_matlab, float_div1d, int_div1d, log_phys)

   start_time = 0.0

   ! only fortran writes a text file
   open( UNIT=10, FILE='div1d_output.txt' )
   call write_header
   call write_solution( start_time )
  
 
   istep = 1
   do ibigstep = 1, nbigsteps
        istate =1
     
   ! run nout time step_t in a function that can also be called by matlab
   call run_div1d(y, ydot, &
                density, velocity, temperature, neutral,&
                Gamma_n, Gamma_mom, pressure, q_parallel,&
                Source_n, Source_v, Source_Q, source_neutral,&
                x, xcb, delta_x, delta_xcb, B_field, B_field_cb,&
                e_charge,, c, K_B, amu, me, pi,&
                options, bounded_components, lower_bounds,&
                upper_bounds, abstol_vector,&
                float_div1d, int_div1d, log_phys)
   ! FORTRAN Writes the solution   
   call write_solution( end_time )
   end do
      
   ! loop is over
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
   ! NOTE that the tag is HARDCODED for backward compatibility in reading the outputs with Matlab! 
   write( 10, * ) 'git tag: v3.0.2'
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
   write( 10, * ) '   num_impurities          = ', num_impurities
   write( 10, * ) '   impurity_concentration  = ', ( impurity_concentration(i), i = 1,num_impurities )
   write( 10, * ) '   impurity_Z              = ', ( impurity_Z(i), i = 1,num_impurities )
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
!   write( 10, * ) '   switch_dyn_nu           = ', switch_dyn_nu 
!   write( 10, * ) '   switch_dyn_gas          = ', switch_dyn_gas 
!   write( 10, * ) '   switch_dyn_rec          = ', switch_dyn_rec
!   write( 10, * ) '   switch_dyn_rad_los      = ', switch_dyn_rad_los
!   write( 10, * ) '   switch_dyn_imp_con      = ', switch_dyn_imp_con
!   write( 10, * ) '   switch_dyn_qpar         = ', switch_dyn_qpar
!   write( 10, * ) '   switch_dyn_red_frc      = ', switch_dyn_red_frc

   write( 10, '(A195)' ) ' X [m]   gas_puff_prf []  B_field [fraction]'  !       rad_los_prf  '
   write( 10, '(13(1PE15.6))' ) ( x(i), gas_puff(i), B_field(i), i=1,Nx )
  !write( 10, '(13(1PE15.6))' ) ( x(i),car_con_prf(i), gas_puff(i),rad_los_prf, i=1,Nx )
   return
end subroutine write_header


subroutine write_solution( time )
   use physics_parameters, only : dyn_gas, dyn_nu, dyn_nb, dyn_rec, dyn_rad_los, dyn_qparX, dyn_red_frc, dyn_imp_con, num_impurities
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
   write( 10, * ) 'dyn_nb      = ', dyn_nu(itime)
   write( 10, * ) 'dyn_rec     = ', dyn_rec(itime)
   write( 10, * ) 'dyn_rad_los = ', dyn_rad_los(itime)
   write( 10, * ) 'dyn_qparX   = ', dyn_qparX(itime)
   write( 10, * ) 'dyn_red_frc = ', dyn_red_frc(itime)
   write( 10, * ) 'dyn_imp_con = ', ( dyn_imp_con(i,itime), i = 1,num_impurities )
   write( 10, '(A195)' ) '    X [m]        N [/m^3]       V [m/s]         T [eV]        Nn [/m^3]      Gamma_n    Gamma_mom [Pa]      q_parallel    neutral_flux     Source_n       Source_v       Source_Q     source_neut  '
   write( 10, '(13(1PE15.6))' ) ( x(i), density(i), velocity(i), temperature(i), neutral(i), &                 !!!! the neutral flux should not be multiplied with the B_field !!!!    
   &                                   Gamma_n(i)*B_field(i),Gamma_mom(i)*B_field(i),q_parallel(i)*B_field(i),neutral_flux(i), & ! multiplied by B_field because the code uses normalized values
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

