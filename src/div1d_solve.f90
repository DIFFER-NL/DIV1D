
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
        
  contains

        subroutine initialize_div1d( y, ydot, &                                 ! plasma params
                             density, velocity, temperature, neutral, &         ! plasma params 
                             Gamma_n, Gamma_mom, pressure, q_parallel, &        ! plasma params
                             Source_n, Source_v, Source_Q, source_neutral, &    ! plasma params 
                             x, xcb, delta_x, delta_xcb, B_field, B_field_cb, & ! grid data     
                             e_charge, c, K_B, amu, me, pi,&                    ! constants
                             options, bounded_components, lower_bounds, &       ! solver options (DVODE only)
                             upper_bounds, abstol_vector, &                     ! solver options (DVODE only)
                             floatindiv1d, intindiv1d, strinphys, loginphys)    ! phys+num params(INPUTS)
                  !bind(c,name="initialize_div1d_")
            ! Returns initial-plasma values, grid data, constants, and solver options
            
        
            call initialize_grid
            call initialize_gas_puff(Nx)
            call read_restart_file( restart_error )
            call initialize_solver_settings( solver_setup_error )
            
                

        end subroutine initialize_div1d


        subroutine run_div1d(y, ydot, &                                         ! plasma params
                             density, velocity, temperature, neutral, &         ! plasma params 
                             Gamma_n, Gamma_mom, pressure, q_parallel, &        ! plasma params
                             Source_n, Source_v, Source_Q, source_neutral, &    ! plasma params 
                             x, xcb, delta_x, delta_xcb, B_field, B_field_cb, & ! grid data     
                             e_charge, c, K_B, amu, me, pi,&                    ! constants
                             options, bounded_components, lower_bounds, &       ! solver options (DVODE only)
                             upper_bounds, abstol_vector, &                     ! solver options (DVODE only)
                             float_div1d, int_div1d, log_phys)                  ! phys+num params                              !bind(c,name="run_div1d_")
       
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

end module div1d_solve
