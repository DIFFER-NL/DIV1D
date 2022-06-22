
module div1d_solve

      ! module containing submodule to call DIV1D from C

      use constants  ! e_charge, c, K_B, amu, me, pi
      use physics_parameters 
      use numerics_parameters
      use grid_data
      use plasma_data ! defines plasma params
      use physics_routines
      use dvode_f90_m
      use experiments
        
      implicit none
      integer, parameter :: wp = KIND(1.0D0)
      real(wp) :: start_time, end_time
      integer :: input_error, input_error_numerics, input_error_physics, solver_setup_error
      integer :: restart_error, time_step_error
      integer :: istep

      ! solver parameters
      integer :: nzswag_input, ix, attempt                      !dvode 
      type (VODE_OPTS)      :: options                          !dvode
      integer,  allocatable :: iwork(:)                         !dvode
      integer,  allocatable :: bounded_omponents(:)             !dvode
      real(wp), allocatable :: lower_bounds(:), upper_bounds(:) !dvode
      real(wp), allocatable :: abstol_vector(:), rwork(:)       !dvode
      integer :: ml, mu, lrw, liw, itol  !dlsode
      integer :: max_step, itask, istate !dlsode and dvode
        
      ! define parameters for external call
      real(wp) floatinphys(27)
      real(wp) floatinnum(24)
      integer  intinphys(6) 
      integer  intinnum(12) 
      logical  loginnum(4)
      !logical, allocatable  :: loginphys(1) 
      integer               :: call_from_extern = 0 ! default value

  contains
        
         subroutine initialize_div1d(density, velocity, temperature, neutral, & ! plasma params (IN/OUT)
                             Gamma_n, Gamma_mom, pressure, q_parallel, &        ! plasma params (OUT)
                             Source_n, Source_v, Source_Q, source_neutral, &    ! plasma params (OUT)
                             x, xcb, delta_x, delta_xcb, B_field, B_field_cb, & ! grid data     (IN/OUT)  
                             e_charge, c, K_B, amu, me, pi,&                    ! constants     (OUT)
                             floatinnum, intinnum, loginnum,&                   ! numerics params(INPUTS)
                             floatinphys, intinphys, &    ! physics params (INPUTS)
                             call_from_extern)                                  ! used as library or on its own
                  !bind(c,name="initialize_div1d_")
            ! Returns initial-plasma values, grid data, constants, and solver options
            implicit none
            real(wp), intent(in) :: floatinnum(24), floatinphys(27)
            integer, intent(in) :: intinnum(12), intinphys(6)
            logical, intent(in) :: loginnum(4)
            integer, intent(in) :: call_from_extern
            real(wp), intent(in) :: e_charge, c, K_B, amu, me, pi ! not allowed to change constants
            !real(wp), intent(out) :: grid data
            !real(wp), intent(out) :: plasma data

            integer :: restart_error, 
            if call_from_extern .eq. 0
            ! read non-default inputs from .txt and .dat files
            call read_numerics_parameters(input_error_numerics)
            call read_physics_parameters(input_error_physics)
            else if call_from_extern .eq. 1
            ! read non-default settings from floatindiv1d and intindiv1d etc.
            call extern_read_numerics_parameters(floatinnum, intinnum, loginnum)
            call extern_read_physics_parameters(floatinphys, intinphys)              
            end if
            
            !if grid undefined (to be added!!)
            call initialize_grid
            ! else take over inputs

            ! if values undefined (to be added!!)
            call initialize_values
            ! else take over inputs

            restart_error = 0.0
            if call_from_extern .eq. 0
                if( restart )  call read_restart_file( restart_error )
                if( restart .and. restart_error .ne. 0 ) call error_report( input_error, restart_error, time_step_error)
            else        
                ! implement restart routine for library
            endif
            call initialize_gas_puff(Nx)
            call initialize_dvode_settings( options, reltol, abstol, lower_bounds, upper_bounds,&
                                            Nx, solver_setup_error )
        end subroutine initialize_div1d


        subroutine run_div1d(density, velocity, temperature, neutral, &         ! plasma params  ! State input & output
                             Gamma_n, Gamma_mom, pressure, q_parallel, &        ! plasma params  ! Output
                             Source_n, Source_v, Source_Q, source_neutral, &    ! plasma params  ! Output
                             start_time, end_time, &
                             imp_con, neu, dneu, ngb, gas, recycle, rad_los, qpar_x, red_frc ) ! external BCs and input
                             !bind(c,name="run_div1d_")
        implicit none
        real(wp), intent(in) :: imp_con, neu, dneu, ngb, gas, recycle, rad_los, qpar_x, red_frc
        real(wp), intent(in) :: end_time 
        real(wp) :: density(Nx), velocity(Nx), temperature(Nx), neutral(Nx)
        real(wp), intent(out) :: Gamma_n(Nx), Gamma_mom(Nx), pressure(Nx), q_parallel(Nx)
        real(wp), intent(out) :: Source_n, Source_v, Source_Q, source_neutral

        call normalize ! redefines normalization factors ! in  MODULE physics_routines
        istate = 1
        do istep = 1, nout
        end_time = start_time + delta_t
        call nvt2y( Nx, density, velocity, temperature, neutral, y)
        call dvode_f90( right_hand_side, 4*Nx, y, start_time, end_time,
        itask, istate, options )
        time_step_error = istate
        if (istate .ne. 2) ! restart if unsuccesfull
        attempt = 0
           do while( istate .ne. 2 .and. attempt .le. max_attempts )
                attempt = attempt + 1
                istate = 1
                call dvode_f90( right_hand_side, 4*Nx, y, start_time, end_time,
                                itask, istate, options )
           enddo
           if( istate .ne. 2 ) then
                call error_report(input_error, restart_error, time_step_error)
           else
                write(*,*) 'success on ', attempt, 'th attempt. continuing'           
           endif     
        endif
        call y2nvt( Nx, y, density, velocity, temperature, neutral )
        call calculate_fluxes( Nx, start_time, density, velocity, temeprature, neutral, Gamma_n,&
                                 Gamma_mom, q_parallel, neutral_flux )
        call calculate_sources( Nx, start_time, density, velocity, temperature, neutral, q_parallel,&
                                 Source_n, Source_v, Source_Q, source_neutral )
        start_time = end_time
        enddo

        end subroutine run_div1d
        
        !subroutine run_div1d_derivative( y, ydot, dy, &                         ! plasma params
        !                     density, velocity, temperature, neutral, &         ! plasma params 
        !                     x, xcb, delta_x, delta_xcb, B_field, B_field_cb, & ! grid data     
        !                     e_charge, c, K_B, amu, me, pi,&                    ! constants
        !                     floatindiv1d, intindiv1d, strinphys, loginphys)    ! phys+num params  ) 
                  !bind(c,name="run_div1d_derivative_") 
        ! calculates the derivative dy/dt informing other blocks in Simulink top structure
        !call nvt2y( Nx, density, velocity, temperature, neutral, y)
        !right_hand_side( neq, time, y , ydot ) 
        ! routine below needs to be build?
        !call dy2dnvt( Nx, dy ,d_density,d_velocity,d_temperature, d_neutral)
        !end subroutine run_div1d_derivative

        subroutine initialize_dvode_settings( options, reltol, abstol, lower_bounds, upper_bounds, Nx )
        ! setting options for dvode_f90
        implicit none
        
        ! setting absolute errror tolerance for dvode_f90 
        ! (note that y is normalized resulting in equal tolerances)
        allocate( abstol_vector(4*Nx) )
        abstol_vector(1:4*Nx) = abstol
        
        ! enforce nonnegative densities and energies by lower bounds in dvode_f90
        allocate ( bounded_components(Nx), lower_bounds(Nx), upper_bounds(Nx) )
        do ix = 1,Nx
                bounded_components( ix) = ix + 2*Nx
        enddo

        ! rough estimate of maximum number of nonzeros in Jacobian 
        nzswag_input = 48 * Nx
        ! build the options structure
        options = set_opts(RELERR=reltol, ABSERR_VECTOR=abstol_vector, &
                                CONSTRAINED=bounded_components, CLOWER=lower_bounds, CUPPER=upper_bounds, &
                                METHOD_FLAG=method, MXSTEP=max_step, NZSWAG=nzswag_input, MA28_ELBOW_ROOM=200, &
                                MA28_RPS=.TRUE. )
        ! call set_jacobian_sparsity_structure
        ! call set_diagnoal_jacobian
        ! call USERSETS_IAJA(IAUSER,NIAUSER,JAUSER,NJAUSER)

        end subroutine initialize_dvode_settings


end module div1d_solve
