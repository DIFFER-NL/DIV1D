
module div1d_solve

      ! module containing submodule to call DIV1D from C

      use constants
      use physics_parameters
      use numerics_parameters
      use grid_data
      use plasma_data
      use physics_routines
      use dvode_f90_m
      use experiments
        
        implicit none

  contains
        


         subroutine initialize_div1d( y, ydot, &                                 ! plasma params
                             density, velocity, temperature, neutral, &         ! plasma params 
                             Gamma_n, Gamma_mom, pressure, q_parallel, &        ! plasma params
                             Source_n, Source_v, Source_Q, source_neutral, &    ! plasma params 
                             x, xcb, delta_x, delta_xcb, B_field, B_field_cb, & ! grid data     
                             e_charge, c, K_B, amu, me, pi,&                    ! constants
                             options, bounded_components, lower_bounds, &       ! solver options (DVODE only)
                             upper_bounds, abstol_vector, &                     ! solver options (DVODE only)
                             floatindiv1d, intindiv1d, strinphys, loginphys,&   ! phys+num params(INPUTS)
                             call_from_matlab)
                  !bind(c,name="initialize_div1d_")
            ! Returns initial-plasma values, grid data, constants, and solver options
            implicit none
            integer, INTENT(IN) :: call_from_matlab
            integer :: restart_error
            if call_from_matlab .eq. 0
            ! read non-default inputs from .txt and .dat files
            call read_numerics_parameters(input_error_numerics)
            call read_physics_parameters(input_error_physics)
            else if call_from_matlab .eq. 1
            ! read non-default settings from floatindiv1d and intindiv1d etc.
            call read_numerics_par_mat(floatindiv1d, intindiv1d, strinphys, loginphys)
            call read_physics_par_mat(floatindiv1d, intindiv1d, strinphys, loginphys)              
            end if
        
            call initialize_grid
            call initialize_values

            restart_error = 0.0
            if( restart )  call read_restart_file( restart_error )
            if restart .and. restart_error .ne. 0 ) call error_report( input_error, restart_error, time_step_error)
            call initialize_gas_puff(Nx)
            call initialize_dvode_settings( solver_setup_error )
        end subroutine initialize_div1d


        subroutine run_div1d(y, &                                               ! plasma params
                             density, velocity, temperature, neutral, &         ! plasma params 
                             Gamma_n, Gamma_mom, pressure, q_parallel, &        ! plasma params
                             Source_n, Source_v, Source_Q, source_neutral, &    ! plasma params 
                             x, xcb, delta_x, delta_xcb, B_field, B_field_cb, & ! grid data     
                             e_charge, c, K_B, amu, me, pi,&                    ! constants
                             options, bounded_components, lower_bounds, &       ! solver options (DVODE only)
                             upper_bounds, abstol_vector, &                     ! solver options (DVODE only)
                             float_div1d, int_div1d, log_phys)                  ! phys+num params                              !bind(c,name="run_div1d_")
        implicit none
        integer, parameter, private :: wp = KIND(1.0D0)
        real( wp ),INTENT(IN) :: float_div1d(10)
        integer,   INTENT(IN) :: int_div1d(10)
        logical,   INTENT(IN) :: log_phys

        real( wp ), :: Nx = float_div1d(1)
        real( wp ), ::


        istate = 1
        do istep = 1, ntime
        end_time = start_time + delta_t
        call nvt2y( Nx, density, velocity, temperature, neutral, y)
        call dvode_f90( right_hand_side, 4*Nx, y, start_time, end_time,
        itask, istate, options )
        if (istate .ne. 2) ! restart if unsuccesfull
        istate = 1
        attempt = 0
        call dvode_f90( right_hand_side, 4*Nx, y, start_time, end_time,
        itask, istate, options )
        endif     
        call y2nvt( Nx, y, density, velocity, temperature, neutral )
        end_time = start_time
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
