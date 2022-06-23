
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
      integer :: input_error, input_error_numerics, input_error_physics !, solver_setup_error
      integer :: restart_error, time_step_error
      integer :: istep

      ! solver parameters
      integer :: nzswag_input, ix, attempt                      !dvode 
      type (VODE_OPTS)      :: options                          !dvode
      integer,  allocatable :: iwork(:)                         !dvode
      integer,  allocatable :: bounded_components(:)             !dvode
      real(wp), allocatable :: lower_bounds(:), upper_bounds(:) !dvode
      real(wp), allocatable :: abstol_vector(:), rwork(:)       !dvode
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

    !  ! parameters that are changed in time externally and need to be known by all modules
    !  real(wp), dimension(5,1) :: E_imp_con
    !  real(wp) :: E_neu, E_dneu, E_ngb, E_gas, E_rec, E_qpar_x, E_red_frc ! rad_los


  contains
        
         subroutine initialize_div1d(density, velocity, temperature, neutral, & ! plasma params (IN/OUT)
                             x, xcb, delta_x, delta_xcb, B_field, B_field_cb, & ! grid data     (IN/OUT)  
                             floatinnum, intinnum, loginnum,&                   ! numerics params(INPUTS)
                             floatinphys, intinphys, &    ! physics params (INPUTS)
                             call_from_extern)                                  ! used as library or on its own
                  !bind(c,name="initialize_div1d_")
 ! e_charge, c, K_B, amu, me, pi,&                    ! constants     (OUT)

            ! Returns initial-plasma values, grid data, constants, and solver options
!    Gamma_n, Gamma_mom, pressure, q_parallel, &        ! plasma params (OUT)
 !                            Source_n, Source_v, Source_Q, source_neutral, &    ! plasma params (OUT)

            implicit none
            real(wp), intent(in) :: floatinnum(24), floatinphys(27)
            integer, intent(in) :: intinnum(12), intinphys(6)
            logical, intent(in) :: loginnum(4)
            integer, intent(in) :: call_from_extern
            !real(wp), parameter, intent(in) :: e_charge, c, K_B, amu, me, pi ! not allowed to change constants
            real(wp) :: x(Nx), xcb(Nx+1), delta_x(Nx), delta_xcb(Nx+1), B_field(Nx), B_field_cb(Nx+1)
            real(wp) :: density(Nx), velocity(Nx), temperature(Nx), neutral(Nx)
  !          real(wp) :: Gamma_n(Nx), Gamma_mom(Nx), pressure(Nx), q_parallel(Nx)
   !         real(wp) :: Source_n(Nx), Source_v(Nx), Source_Q(Nx), source_neutral(Nx)

            !real(wp), intent(out) :: plasma data

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
            write(*,*) 'going to initialize grid'
            !if grid undefined (to be added!!)
            call initialize_grid
            ! else take over inputs
            write(*,*) 'going to initialize values'
            ! if values undefined (to be added!!)
            call initial_values
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

        end subroutine initialize_div1d


        subroutine run_div1d(density, velocity, temperature, neutral, &         ! plasma params  ! State input & output
                             Gamma_n, Gamma_mom, q_parallel, neutral_flux, &        ! plasma params  ! Output
                             Source_n, Source_v, Source_Q, source_neutral, &    ! plasma params  ! Output
                             start_time, end_time, nout, delta_t, &
                             E_imp_con, E_neu, E_dneu, E_ngb, E_gas, E_rec,  E_qpar_x, E_red_frc ) ! xternal BCs and input
                             !bind(c,name="run_div1d_")
        implicit none
        integer, intent(in) :: nout
        real(wp), intent(in) :: delta_t
        real(wp), dimension(5,1), intent(in) :: E_imp_con
        real(wp), intent(in) :: E_neu, E_dneu, E_ngb, E_gas, E_rec, E_qpar_x, E_red_frc ! rad_los
        real(wp) :: end_time, start_time 
        real(wp) :: density(Nx), velocity(Nx), temperature(Nx), neutral(Nx) ! both input and output
        real(wp), intent(out) :: Gamma_n(Nx), Gamma_mom(Nx), q_parallel(Nx), neutral_flux(Nx)
        real(wp), intent(out) :: Source_n(Nx), Source_v(Nx), Source_Q(Nx), source_neutral(Nx)
        integer :: itime = 10
        integer :: i 
        write(*,*) 'imp_con', E_imp_con
        write(*,*) 'neu ', E_neu
        write(*,*) 'dneu ',E_dneu 
        write(*,*) 'ngb ', E_ngb
        write(*,*) 'gas ',E_gas 
   write(*, * ) 'dyn_gas     = ', dyn_gas(itime)
   write(*, * ) 'dyn_nu      = ', dyn_nu(itime)
   write(*, * ) 'dyn_nb      = ', dyn_nb(itime)
   write(* , * ) 'dyn_rec     = ', dyn_rec(itime)
   write(* , * ) 'dyn_rad_los = ', dyn_rad_los(itime)
   write(* , * ) 'dyn_qparX   = ', dyn_qparX(itime)
   write(* , * ) 'dyn_red_frc = ', dyn_red_frc(itime)
   write(* , * ) 'dyn_imp_con = ', ( dyn_imp_con(i,itime), i = 1,num_impurities )
       
        call normalize(density_norm, temperature_norm, velocity_norm, momentum_norm, energy_norm, neutral_norm, &
                        density, velocity, temperature, neutral, renormalize )
                !redefines normalization factors ! in  MODULE physics_routines
        istate = 1
        itask = 1 ! for normal computation in dvode_f90 till end_time
        do istep = 1, nout
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
      !  call calculate_fluxes( Nx, start_time, density, velocity, temperature, neutral, Gamma_n,&
      !                           Gamma_mom, q_parallel, neutral_flux )
      !  call calculate_sources( Nx, start_time, density, velocity, temperature, neutral, q_parallel,&
      !                           Source_n, Source_v, Source_Q, source_neutral )
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
        ! calculates the derivative dy/dt informing other blocks in Simulink top structure
        !call nvt2y( Nx, density, velocity, temperature, neutral, y)
        !right_hand_side( neq, time, y , ydot ) 
        ! routine below needs to be build?
        !call dy2dnvt( Nx, dy ,d_density,d_velocity,d_temperature, d_neutral)
        !end subroutine run_div1d_derivative

    
end module div1d_solve
