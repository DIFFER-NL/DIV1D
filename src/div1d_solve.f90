
module div1d_solve

      ! module containing subroutines to call DIV1D externally
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
      integer :: input_error, input_error_numerics, input_error_physics !, solver_setup_error
      integer :: restart_error, time_step_error
      !integer :: istep

      ! solver parameters
      integer :: nzswag_input, ix, attempt                       
      type (VODE_OPTS)      :: options                         
      integer,  allocatable :: iwork(:)                         
      integer,  allocatable :: bounded_components(:)           
      real(wp), allocatable :: lower_bounds(:), upper_bounds(:)
      real(wp), allocatable :: abstol_vector(:), rwork(:)      
      integer :: ml, mu, lrw, liw, itol  !dlsode
      integer :: itask, istate !dlsode and dvode
      !integer :: max_step ! already available through USE statement  
      ! define parameters for external call
      real(wp), dimension(27) :: floatinphys = 0.0
      real(wp), dimension(24) :: floatinnum = 0.0
      integer, dimension(6)  :: intinphys = 0
      integer, dimension(12)  :: intinnum = 0
      logical, dimension(4)  :: loginnum = .true. 
      !logical, allocatable  :: loginphys(1) 
      integer               :: call_from_extern = 0 ! default value

 

  contains
        
        subroutine initialize_div1d_settings(floatinnum, intinnum, loginnum,&  ! numerics params(INPUTS)
                             floatinphys, intinphys, &                          ! physics params (INPUTS)
                             call_from_extern)                                  ! used as library or on its own
                  !bind(c,name="initialize_div1d_settings_")
            implicit none
            real(wp), intent(in) :: floatinnum(24), floatinphys(27)
            integer, intent(in) :: intinnum(12), intinphys(6)
            logical, intent(in) :: loginnum(4)
            integer, intent(in) :: call_from_extern

            integer :: restart_error 
            if( call_from_extern .eq. 0 ) then
            ! read non-default inputs from .txt and .dat files
            call read_numerics_parameters(input_error_numerics)
            call read_physics_parameters(input_error_physics)
            elseif( call_from_extern .eq. 1 ) then
            ! read non-default settings from floatindiv1d and intindiv1d etc.
            call extern_read_numerics_parameters(floatinnum, intinnum, loginnum)
            call extern_read_physics_parameters(floatinphys, intinphys)              
            end if
            ! use settings to allocate grid and plasma vectors
            allocate( x(Nx), xcb(Nx+1), delta_x(Nx), delta_xcb(Nx), B_field(Nx), B_field_cb(Nx+1) )
            allocate( density(Nx), velocity(Nx), temperature(Nx), neutral(Nx) )   
            allocate( E_imp_con(5,nout), E_neu(nout), E_dneu(nout), E_ngb(nout), E_gas(nout))
            allocate( E_rec(nout), E_qpar_x(nout), E_red_frc(nout) )          
    
            return  
        end subroutine initialize_div1d_settings

        subroutine initialize_div1d(E_density, E_velocity,E_temperature, E_neutral, & ! plasma params (IN)
                             E_x, E_xcb, E_delta_x, E_delta_xcb, E_B_field, E_B_field_cb, & ! grid data     (IN)  
                             density,velocity ,temperature,neutral, & ! plasma params (IN/OUT)
                             x, xcb, delta_x, delta_xcb, B_field, B_field_cb, & ! grid data     (IN/OUT)  
                             call_from_extern, Nx )  
                             !bind(c,name="initialize_div1d_")
            implicit none
            integer, intent(in) :: Nx, call_from_extern
            real(wp) :: E_x(Nx), E_xcb(Nx+1), E_delta_x(Nx), E_delta_xcb(Nx), E_B_field(Nx), E_B_field_cb(Nx+1)
            real(wp) :: x(Nx), xcb(Nx+1),delta_x(Nx), delta_xcb(Nx+1), B_field(Nx), B_field_cb(Nx+1)
            real(wp) :: E_density(Nx), E_velocity(Nx), E_temperature(Nx), E_neutral(Nx)
            real(wp) :: density(Nx), velocity(Nx), temperature(Nx), neutral(Nx)
                   
            write(*,*) 'going to initialize grid'
            if( E_x(1) .le. 0 ) then
            call initialize_grid
                E_x             = x
                E_xcb           = xcb
                E_delta_x       = delta_x
                E_delta_xcb     = delta_xcb
                E_B_field       = B_field
                E_B_field_cb    = B_field_cb
            else
                x               = E_x
                xcb             = E_xcb
                delta_x         = E_delta_x
                delta_xcb       = E_delta_xcb
                B_field         = E_B_field
                B_field_cb      = E_B_field_cb
             endif
             write(*,*) 'going to initialize values'
             if( E_density(1) .le. 0) then ! if the external density is undefined, e.g. 0 or -1 use native intialization
             if( call_from_extern .eq. 1 ) then
              write(*,*) 'assign physics arrays'
              call assign_physics_arrays
             end if
             write(*,*) 'set initial values 1'
             call initial_values
                E_density       = density
                E_velocity      = velocity
                E_temperature   = temperature
                E_neutral       = neutral
             else
             write(*,*) 'set initial values 2'
                density         = E_density
                velocity        = E_velocity
                temperature     = E_temperature
                neutral         = E_neutral
             endif 
            ! else take over inputs
            write(*, *) 'restarting'
            restart_error = 0.0
            if( call_from_extern .eq. 0 ) then
                if( restart )  call read_restart_file( restart_error )
                if( restart .and. restart_error .ne. 0 ) call error_report( input_error, restart_error, time_step_error)
            else        
                ! implement restart routine for library
            endif

            write(*,*) 'going to initilize gas puff'
            call initialize_gas_puff(Nx)
            
            write(*,*) 'going to initialize solver settings'
            ! ----------------------------------DVODE solver settings ------------------- !
            allocate( abstol_vector(4*Nx) )
            abstol_vector(1:4*Nx) = abstol
            ! enforce nonnegative densities and energies by lower bounds in dvode_f90
            allocate( bounded_components(Nx), lower_bounds(Nx), upper_bounds(Nx) )
            do ix = 1,Nx
                bounded_components( ix) = ix + 2*Nx
            enddo
            write(*,*) 'set nzswag'
            ! rough estimate of maximum number of nonzeros in Jacobian 
            nzswag_input = 48 * Nx
            ! build the options structure
            write(*,*) 'set options'
            options = set_opts(RELERR=reltol, ABSERR_VECTOR=abstol_vector, &
                                CONSTRAINED=bounded_components, CLOWER=lower_bounds, CUPPER=upper_bounds, &
                                METHOD_FLAG=method, MXSTEP=max_step, NZSWAG=nzswag_input, MA28_ELBOW_ROOM=200, &
                                MA28_RPS=.TRUE. )
            write(*,*) 'solver settings initialized'
            return
        end subroutine initialize_div1d


        subroutine run_div1d(density, velocity, temperature, neutral, &     ! plasma params  ! input & output
                             Gamma_n, Gamma_mom, q_parallel, neutral_flux, &! plasma params  ! Output
                             Source_n, Source_v, Source_Q, source_neutral, &! plasma params  ! Output
                             start_time, end_time, nout, delta_t, &
                             E_imp_con, E_neu, E_dneu, E_ngb, E_gas, E_rec,  E_qpar_x, E_red_frc) ! external BC/inputs
                            
                             !bind(c,name="run_div1d_")
        implicit none
        integer, intent(in) :: nout
        real(wp), intent(in) :: delta_t
        real(wp), dimension(5,nout), intent(in) :: E_imp_con
        real(wp), dimension(nout), intent(in) :: E_neu, E_dneu, E_ngb, E_gas, E_rec, E_qpar_x, E_red_frc      
        real(wp) :: end_time, start_time 
        real(wp) :: density(Nx), velocity(Nx), temperature(Nx), neutral(Nx) ! both input and output
        real(wp), intent(out) :: Gamma_n(Nx), Gamma_mom(Nx), q_parallel(Nx), neutral_flux(Nx)
        real(wp), intent(out) :: Source_n(Nx), Source_v(Nx), Source_Q(Nx), source_neutral(Nx)
        !integer :: itime = 10
        !integer :: i 
        !write(*,*) 'imp_con', E_imp_con
        !write(*,*) 'neu ', E_neu
        !write(*,*) 'dneu ',E_dneu 
        !write(*,*) 'ngb ', E_ngb
        !write(*,*) 'gas ',E_gas 
        !write(*,*) 'qpar', E_qpar_x
        !write(*,*) 'red frc', E_red_frc

         
        call normalize(density_norm, temperature_norm, velocity_norm, momentum_norm, energy_norm, neutral_norm, &
                        density, velocity, temperature, neutral, renormalize )
                !redefines normalization factors ! in  MODULE physics_routines
        istate = 1
        itask = 1 ! for normal computation in dvode_f90 till end_time
        do internal_istep = 1, nout
        !internal_istep = istep ! pass istep to the physics_routine module
        end_time = start_time + delta_t
        call nvt2y( Nx, density, velocity, temperature, neutral, y)
        call dvode_f90( right_hand_side, 4*Nx, y, start_time, end_time, itask, istate, options )
        time_step_error = istate
        if(istate .ne. 2) then ! restart if unsuccesfull
        attempt = 0
           do while( istate .ne. 2 .and. attempt .le. max_attempts )
                attempt = attempt + 1
                istate = 1
                call dvode_f90( right_hand_side, 4*Nx, y, start_time, end_time, &
                                itask, istate, options )
           enddo
           if( istate .ne. 2 ) then
                call error_report(input_error, restart_error, time_step_error)
           else
                write(*,*) 'success on ', attempt, 'th attempt. continuing'           
           endif     
        endif
        call y2nvt( Nx, y, density, velocity, temperature, neutral )
        start_time = end_time
        enddo
        call calculate_fluxes( Nx, start_time, density, velocity, temperature, neutral, Gamma_n,&
                                 Gamma_mom, q_parallel, neutral_flux )
        call calculate_sources( Nx, start_time, density, velocity, temperature, neutral, q_parallel,&
                                 Source_n, Source_v, Source_Q, source_neutral )
        end subroutine run_div1d
        
        !subroutine run_div1d_derivative( y, ydot, dy, &                         ! plasma params
        !                     density, velocity, temperature, neutral, &         ! plasma params 
        !                     x, xcb, delta_x, delta_xcb, B_field, B_field_cb, & ! grid data     
        !                     e_charge, c, K_B, amu, me, pi,&                    ! constants
        !                     floatindiv1d, intindiv1d, strinphys, loginphys)    ! phys+num params  ) 
                  !bind(c,name="run_div1d_derivative_") 
        ! calculates the derivative dy/dt informing other blocks in Simulink 
        !call nvt2y( Nx, density, velocity, temperature, neutral, y)
        !right_hand_side( neq, time, y , ydot ) 
        ! routine below needs to be build?
        !call dy2dnvt( Nx, dy ,d_density,d_velocity,d_temperature, d_neutral)
        !end subroutine run_div1d_derivative

    
end module div1d_solve
