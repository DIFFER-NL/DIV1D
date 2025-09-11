
module div1d_step

      ! module containing subroutines to call DIV1D externally
      use constants  
      use physics_parameters 
      use numerics_parameters
      use grid_data
      use plasma_data ! containts y, ys, yr, yc
      use physics_routines ! variables in local scope of physics_routines can be set from div1d_step.f90
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
      
      ! define parameters for external call
      real(wp), dimension(141) :: floatinphys = 0.0
      real(wp), dimension(90) :: floatinnum = 0.0
      integer, dimension(31)  :: intinphys = 0
      integer, dimension(30)  :: intinnum = 0
      logical, dimension(20)  :: loginnum = .true. 
      !logical, allocatable  :: loginphys(1) 
      integer               :: call_from_extern = 0 ! default value
      integer               :: init_grid_fortran = 1 ! default value
      integer               :: init_prof_fortran = 1 ! default value

  

  contains
        
        subroutine initialize_div1d_settings(floatinnum, intinnum, loginnum,&  ! numerics params(INPUTS)
                             floatinphys, intinphys, &                          ! physics params (INPUTS)
                             call_from_extern)     bind(C, name="initialize_div1d_settings_")  ! used as library or on its own
                 
            implicit none
            real(wp), intent(in) :: floatinnum(90), floatinphys(141)
            integer, intent(in) :: intinnum(30), intinphys(31)
            logical, intent(in) :: loginnum(10)
            integer, intent(in) :: call_from_extern
            integer :: counter
            integer :: restart_error 
            real(wp) :: tmpval
		
            if( call_from_extern .eq. 0 ) then
                ! read non-default inputs from .txt and .dat files
                !write(*,*) 'F reading input.txt files
                ! Write the output text file
                open( UNIT=10, FILE='div1d_output.txt' )
                
                write( 10, * ) 'git tag: v6.0.0' ! version of DIV1D
                call read_numerics_parameters(input_error_numerics)
                
                write(10, div1d_numerics)
                call read_physics_parameters(input_error_physics)
                
                write(10,div1d_physics) 
            elseif( call_from_extern .eq. 1 ) then
                open( UNIT=10, FILE='echo_div1d_inputs.txt' )
                
                write( 10, * ) 'git tag: v6.0.0' ! version of DIV1D
                ! read non-default settings from floatindiv1d and intindiv1d etc.
                call extern_read_numerics_parameters(floatinnum, intinnum, loginnum)
                
                write(10, div1d_numerics)
                call extern_read_physics_parameters(floatinphys, intinphys)    
                stop
                
                write(10, div1d_physics)       			   
            end if
		    close(10 )
            
            call allocate_grid_plasma_vectors()
            
	        !write(*,*) "F Reading floatinnum in Fortran"
            !write(*,*) "F  floatinnum ", floatinnum(:)
            !write(*,*) "F  intinnum   ", intinnum(:)
            !write(*,*) "F  floatinphys", floatinphys(:)
            !write(*,*) "F  intinphys  ", intinphys(:)
            write(*,*) "F Initialized settings in Fortran"           
        end subroutine initialize_div1d_settings
        
        subroutine allocate_grid_plasma_vectors()
        
            ! allocate grid and plasma vectors
            allocate( x(Nx), xcb(0:Nx), delta_x(Nx-1), delta_xcb(Nx), B_field(Nx), B_field_cb(0:Nx), B_trans(Nx), B_trans_cb(0:Nx)  )
	        allocate( R_cc(Nx), R_cb(0:Nx), Area_extern(Nx), sintheta_cc(Nx), sintheta_cb(0:Nx), sol_width_pol(Nx), sol_width_pol_cb(0:Nx), volumes(Nx), A_int(1:Nx) )
	        allocate( Z_cc(Nx), Z_cb(0:Nx), nr_cc(Nx), nz_cc(Nx) )
		    !allocate( Area_extern(Nx) )	    
	        allocate( gas_puff_profile(Nx) , core_source_profile_Q(Nx), core_source_profile_n(Nx) ) 
            allocate( D_imp_con(5,nout), D_neu(nout), D_dneu(nout), D_nb(5,nout), D_mb(5,nout), D_gas(nout))
            allocate( D_rec(nout), D_qpar_x(nout), D_red_frc(nout), D_Q_core(nout), D_Gamma_core(nout) ) 
	        allocate( D_neutral_puff(5,nout), D_molecule_puff(5,nout), D_core_fuelling(nout), D_core_neutral_density(nout) )
            allocate( y(6*Nx+10), ydot(6*Nx+10), ys(6*Nx) )
	        allocate( density(Nx), velocity(Nx), temperature(Nx), neutral(Nx), neutral_velocity(Nx), molecule(Nx) )
            allocate( Gamma_n(0:Nx), Gamma_mom(0:Nx), pressure(Nx), q_parallel(0:Nx), Gamma_neutral(0:Nx), Gamma_mom_neutral(0:Nx), Gamma_molecule(0:Nx) ) !, vesrz(1,3) )
	        allocate( Source_n(Nx), Source_v(Nx), Source_Q(Nx), Source_neutral(Nx), Source_vn(Nx), Source_molecule(Nx) )
	        allocate( extern2sol_flux(Nx), core2sol_flux(Nx), extern2sol_mol(Nx), core2sol_mol(Nx), sol2extern_ion_flux(Nx) )
        end subroutine allocate_grid_plasma_vectors
        
        subroutine initialize_div1d_arrays(E_density, E_velocity,E_temperature, E_neutral, E_neutral_velocity, E_molecule,  & ! plasma params (IN/OUT)
                             E_x, E_xcb, E_delta_x, E_delta_xcb, E_B_field, E_B_field_cb, E_B_trans, E_B_trans_cb, &
			     E_R_cc, E_R_cb, E_Area_extern, E_sintheta_cc, E_sintheta_cb, E_sol_width_pol, E_sol_width_pol_cb, E_volumes, &
			     E_gas_puff_profile, E_core_source_profile_Q, E_core_source_profile_n, & ! grid data     (IN/OUT)  
                             init_grid_fortran, init_prof_fortran, Nx, E_i_omp, E_i_Xpoint, E_i_baffle, E_mid_point, E_X_omp )  bind(C,name="initialize_div1d_arrays_")
            ! use allocated vectors from the DIV1D program
            use grid_data,      only :  i_omp, i_Xpoint, i_baffle, mid_point, X_omp, x, xcb, delta_x, delta_xcb, B_field, B_field_cb, &
					R_cc, R_cb, Area_extern, sintheta_cc, sintheta_cb, sol_width_pol, sol_width_pol_cb, volumes, &
					gas_puff_profile, core_source_profile_Q, core_source_profile_n
            use plasma_data,    only :  density, velocity, temperature, neutral, neutral_velocity, molecule, extern_neutral_density, extern_molecule_density
            ! this function is called in div1d_main.f90, libdiv1d.h, and test_libdiv1d.c
            implicit none
	    ! define variables used in the interface of subroutine: initialize_div1d_arrays
            integer, intent(in) :: Nx, init_grid_fortran, init_prof_fortran
	    integer :: E_i_omp, E_i_Xpoint(2), E_i_baffle(2) !, i_omp, i_Xpoint(2), mid_point
	    real(wp) :: E_X_omp, E_mid_point !, X_omp
            real(wp) :: E_x(Nx), E_xcb(0:Nx), E_delta_x(Nx-1), E_delta_xcb(Nx), E_B_field(Nx), E_B_field_cb(0:Nx)
	    real(wp) :: E_B_trans(Nx), E_B_trans_cb(0:Nx), E_R_cc(Nx), E_R_cb(0:Nx), E_Area_extern(Nx)
	    real(wp) :: E_sintheta_cc(Nx), E_sintheta_cb(0:Nx), E_sol_width_pol(Nx), E_sol_width_pol_cb(0:Nx), E_volumes(Nx)
	    real(wp) :: E_gas_puff_profile(Nx), E_core_source_profile_Q(Nx),  E_core_source_profile_n(Nx)
            real(wp) :: E_density(Nx), E_velocity(Nx), E_temperature(Nx), E_neutral(Nx), E_neutral_velocity(Nx), E_molecule(Nx)
       
            write(*,*) 'F going to initialize grid'

            if( init_grid_fortran .eq. 1 ) then
                write(*,*) 'F  using internal DIV1D grid setup'
		call initialize_grid(Nx, x, xcb, B_field, B_field_cb, B_trans, B_trans_cb, &
		            R_cc, R_cb, Area_extern, sintheta_cc, sintheta_cb, sol_width_pol, sol_width_pol_cb, volumes, &
			    gas_puff_profile, E_core_source_profile_Q, E_core_source_profile_n, i_omp, i_Xpoint, i_baffle, A_int)
		E_i_omp		= i_omp
		E_i_Xpoint	= i_Xpoint
		E_mid_point	= mid_point
		E_X_omp		= E_X_omp                
		E_x             = x
                E_xcb           = xcb
                E_delta_x       = delta_x
                E_delta_xcb     = delta_xcb
                E_B_field       = B_field
                E_B_field_cb    = B_field_cb
		E_B_trans	= B_trans	
		E_B_trans_cb	= B_trans_cb
		E_R_cc		= R_cc
		E_R_cb		= R_cb
		E_Area_extern	= Area_extern
		E_sintheta_cc	= sintheta_cc
		E_sintheta_cb	= sintheta_cb
		E_sol_width_pol	= sol_width_pol
		E_sol_width_pol_cb	= sol_width_pol_cb
		E_volumes	= volumes
		E_gas_puff_profile	= gas_puff_profile
		E_core_source_profile_Q	= core_source_profile_Q
		E_core_source_profile_n	= core_source_profile_n
            else
                write(*,*) 'F   using external grid'
		i_omp		= E_i_omp
		i_Xpoint	= E_i_Xpoint
		i_baffle	= E_i_baffle
		mid_point	= E_mid_point
		X_omp		= E_X_omp
                x               = E_x
                xcb             = E_xcb
                delta_x         = E_delta_x
                delta_xcb       = E_delta_xcb
                B_field         = E_B_field
                B_field_cb      = E_B_field_cb
		B_trans		= E_B_trans	
		B_trans_cb	= E_B_trans_cb
		R_cc		= E_R_cc
		R_cb		= E_R_cb
		Area_extern	= E_Area_extern
		sintheta_cc	= E_sintheta_cc
		sintheta_cb	= E_sintheta_cb
		sol_width_pol	= E_sol_width_pol
		sol_width_pol_cb	= E_sol_width_pol_cb
		volumes			= E_volumes
		gas_puff_profile	= E_gas_puff_profile
		core_source_profile_Q	= E_core_source_profile_Q
		core_source_profile_n	= E_core_source_profile_n
             endif

             write(*,*) 'F going to initialize profile values'
             ! if( E_density(1) .le. 0) then ! if the external density is undefined, e.g. 0 or -1 use native initialization
             if( init_prof_fortran .eq. 1 ) then
                if( restart ) then
                   write(*, *) 'F   restarting from old file'
                   restart_error = 0.0
                   call read_restart_file( restart_error )
                   if( restart .and. restart_error .ne. 0 ) call error_report( input_error, restart_error, time_step_error)
               elseif( simple_sol ) then   
		  write(*,*) 'F   using init_simple_sol from DIV1D'
		  call init_simple_sol
	       else
                  write(*,*) 'F   using initial_values from DIV1D'
                  call initial_values
                  !write(*,*) 'F density cell after initial values 30', density(30)
		  !TODO add simple sol here 
                endif
                write(*,*) 'F   pass internal profiles back to interface'
                E_density       = density
                E_velocity      = velocity
                E_temperature   = temperature
                E_neutral       = neutral
		E_neutral_velocity = neutral_velocity
		E_molecule	= molecule
		! should add to passback the other densities etc? no, just specify that in the initial conditions. 
            else
                write(*,*) 'F   using external initial profile values'
                density         = E_density
                velocity        = E_velocity
                temperature     = E_temperature
                neutral         = E_neutral
		neutral_velocity = E_neutral_velocity
		molecule	= E_molecule
                ! write(*,*) 'F density cell 30', density(30)
             endif 
        
            write(*,*) 'F going to initialize solver settings'
            ! ----------------------------------DVODE solver settings ------------------- !
            if( allocated(abstol_vector) .eqv. .false.) allocate( abstol_vector(6*Nx+10) )
            abstol_vector(1:6*Nx+10) = abstol
            ! enforce nonnegative densities and energies by lower bounds in dvode_f90
            if( allocated(bounded_components) .eqv. .false.) allocate( bounded_components(Nx), lower_bounds(Nx), upper_bounds(Nx) )
            do ix = 1,Nx
                bounded_components( ix) = ix + 2*Nx
            enddo
            ! write(*,*) 'F set nzswag'
	    lower_bounds = 0.0d+0
	    upper_bounds = 1.0d+40
	    ! possibly we can add lower and upper bounds here for other variables
	    ! might work better than max,min inside physics_routines

            ! rough estimate of maximum number of nonzeros in Jacobian 
            nzswag_input = 48 * Nx
            ! build the options structure
            ! write(*,*) 'F set options'
            options = set_opts(RELERR=reltol, ABSERR_VECTOR=abstol_vector, &
                                CONSTRAINED=bounded_components, CLOWER=lower_bounds, CUPPER=upper_bounds, &
                                METHOD_FLAG=method, MXSTEP=max_step, NZSWAG=nzswag_input, MA28_ELBOW_ROOM=200, &
                                MA28_RPS=.TRUE. )
            write(*,*) 'F initialized_div1d_arrays'
            return
        end subroutine initialize_div1d_arrays

        subroutine run_div1d(E_density, E_velocity, E_temperature, E_neutral, E_neutral_velocity, E_molecule, &     !6 plasma params  ! input & output
                             E_Gamma_n, E_Gamma_mom, E_q_parallel, E_Gamma_neutral, E_Gamma_mom_neutral, E_Gamma_molecule, &!12 plasma params  ! Output
			     E_core2sol_flux, E_core2sol_mol, E_sol2extern_ion_flux, & ! 15
      			     E_Source_n, E_Source_v, E_Source_Q, E_Source_neutral, E_Source_vn, E_Source_molecule, & !21
			     E_extern_neutral_density, & ! 22 input and output
			     E_extern_neutral_flux, E_sol2extern_flux, E_extern2sol_flux, E_tar2extern_flux, & !26 neutral fluxes related external neutral volumes ! Output
      			     E_extern2core_flux, E_sum_sol2extern_ion_flux, & ! 28
			     E_Source_extern, E_neutral_pump, & ! 30
			     E_extern_molecule_density, & ! 31
			     E_extern_molecule_flux, E_sol2extern_mol, E_extern2sol_mol, E_tar2extern_mol, & !35 neutral fluxes related external neutral volumes ! Output
      			     E_extern2core_mol, E_sum_sol2extern_ion_mol,  &! 37
                             E_Source_extern_mol, E_molecule_pump, & ! 39
			     E_core_density, E_Gamma_core2sol, E_sol2core_flux, E_sol2core_mol, E_Source_core, &! 44 plasma params  ! Output
                             E_start_time, E_end_time, E_nout, E_delta_t, & ! 48
                             E_imp_con, E_neu, E_dneu, E_nb, E_mb, E_gas, E_rec,  E_qpar_x, E_red_frc, & ! 57
			     E_Q_core, E_Gamma_core, E_core_neutral_density, & ! 60
			     E_neutral_puff, E_molecule_puff, E_core_fuelling) bind(C,name="run_div1d_") ! 63
	! E_sum_sol2extern_ion_flux 

        ! this function is called in div1d_main.f90, libdiv1d.h, and test_libdiv1d.c
        implicit none
        integer , intent(in) :: E_nout  ! should not be larger than nout given to numerics parameters.
        real(wp), intent(in) :: E_delta_t
        real(wp), intent(in), dimension(5,E_nout) :: E_imp_con, E_neutral_puff, E_molecule_puff, E_nb, E_mb
        real(wp), intent(in), dimension(E_nout) :: E_neu, E_dneu, E_gas, E_rec, E_qpar_x, E_red_frc, E_Q_core, E_Gamma_core, E_core_fuelling
	! real(wp), intent(in) :: E_pump_n, E_pump_m
        real(wp), intent(in) :: E_start_time ! start time is given to routine 
        real(wp) :: E_end_time ! end time rolls out after nout*delta_t + t_start
        real(wp) :: E_density(Nx), E_velocity(Nx), E_temperature(Nx), E_neutral(Nx), E_neutral_velocity(Nx), E_molecule(Nx) ! both input and output
        real(wp) :: E_Gamma_n(0:Nx), E_Gamma_mom(0:Nx), E_q_parallel(0:Nx), E_Gamma_neutral(0:Nx), E_Gamma_mom_neutral(0:Nx), E_Gamma_molecule(0:Nx)
	real(wp) :: E_core2sol_flux(Nx), E_core2sol_mol(Nx), E_sol2extern_ion_flux(Nx)
        real(wp) :: E_Source_n(Nx), E_Source_v(Nx), E_Source_Q(Nx), E_Source_neutral(Nx), E_Source_vn(Nx), E_Source_molecule(Nx)
	real(wp) :: elm_heat_load, elm_density_change, elm_core_particle_source !, leakage_fluxes_n(5), leakage_fluxes_m(5)
	real(wp) :: atom_association_sink(5), molecule_association_source(5)
  
   
	real(wp) :: E_extern_neutral_density(5)  
        real(wp) :: E_extern_neutral_flux(3), E_sol2extern_flux(5), E_extern2sol_flux(Nx), E_tar2extern_flux(5) ! neutral fluxes external neutral volumes
        real(wp) :: E_extern2core_flux(5),   E_sum_sol2extern_ion_flux  
	real(wp) :: E_neutral_pump(5)
	real(wp) :: E_Source_extern(5)

	real(wp) :: E_extern_molecule_density(5)
        real(wp) :: E_extern_molecule_flux(3), E_sol2extern_mol(5), E_extern2sol_mol(Nx), E_tar2extern_mol(5) ! mol fluxes external molecule volumes
        real(wp) :: E_extern2core_mol(5), E_sum_sol2extern_ion_mol   ! fluxes molecules to the core ? add ? E_sol2extern_ion_mol
	real(wp) :: E_molecule_pump(5)
	real(wp) :: E_Source_extern_mol(5)

	real(wp) :: E_core_density
	real(wp) :: E_core_neutral_density(E_nout) !(can be input and output for now)
	real(wp) :: E_Gamma_core2sol, E_sol2core_flux, E_sol2core_mol, E_Source_core   
	
	! solution vectors
	real(wp) :: yc, ys(6*Nx), yr(10) 

	! account for core ionization fraction
	real(wp) :: extern2core_shinethrough(Nx)
	real(wp) :: extern2core_shinethrough_mol(Nx)
	real(wp) :: zeroNx(Nx), zero3(3)
        extern2core_shinethrough(Nx) = 0.0d+0
        extern2core_shinethrough_mol(Nx) = 0.0d+0

	zeroNx = 0.0d+0
	zero3 = 0.0d+0
	E_Source_core = 0.0d+0

	! checks to be added
 	! (E_end_time - E_start_time)/E_delta_t = E_nout
	! E_nout <= nout (otherwise we have to allocate more memory)
	
        !write(*,*)  'F imp_con: sz, 11 21', size(E_imp_con), E_imp_con(1,1), E_imp_con(2,1)
        !write(*,*) 'F neu ', E_neu
        !write(*,*) 'F dneu ',E_dneu 
        !write(*,*) 'F ngb ', E_nb
        !write(*,*) 'F gas ',E_gas 
        !write(*,*) 'F qpar', E_qpar_x
        !write(*,*) 'F red frc', E_red_frc
        
 
        write(*, *) 'F main is calling: overwriting internal variables' 
        nout            = E_nout  ! set the size nout in D_vars, dynamic variables.
        start_time      = E_start_time
	! plasma
        density         = E_density
        velocity        = E_velocity
        temperature     = E_temperature
        neutral         = E_neutral
	neutral_velocity = E_neutral_velocity
	molecule 	= E_molecule
	! atoms
	extern_neutral_density = E_extern_neutral_density	
	! molecules
	extern_molecule_density = E_extern_molecule_density
	! core 
	core_density = E_core_density
	!write(*,*) 'if values above where not initialized outsided this routine and have no refrence memory adress this will give a segmentation fault'

	! if some of these fluxes have to be fixed from outside, this is the point to prescribe them       
	!extern_neutral_flux = E_extern_neutral_flux
	!sol2extern_flux = E_sol2extern_flux
	!extern2sol_flux = E_extern2sol_flux
	!tar2extern_flux = E_tar2extern_flux      
	!Source_extern	= E_Source_extern

	!extern_molecule_flux = E_extern_molecule_flux
	!sol2extern_mol = E_sol2extern_mol
	!extern2sol_mol = E_extern2sol_mol
	!tar2extern_mol = E_tar2extern_mol
	!Source_extern_mol   = E_Source_extern_mol
	
	!extern2core_flux = E_extern2core_flux(5)
	!extern2core_mol = E_extern2core_mol(5)
	!Gamma_core2sol = E_Gamma_core2sol
	!Source_core = E_Source_core 	
	
	!write(*,*) 'F overwrite dynamic variables used by right-hand-side'
	D_start_time = E_start_time
	D_delta_t = E_delta_t
	D_nout = E_nout
	
	D_imp_con(1:5,1:E_nout) = E_imp_con
        D_neu(1:E_nout) = E_neu
	D_dneu(1:E_nout) = E_dneu
	D_nb(1:5,1:E_nout) = E_nb
	D_mb(1:5,1:E_nout) = E_mb
	D_gas(1:E_nout) = E_gas 
	D_rec(1:E_nout) = E_rec
	D_qpar_x(1:E_nout) = E_qpar_x
	D_red_frc(1:E_nout) = E_red_frc
	D_Q_core(1:E_nout) = E_Q_core
	D_Gamma_core(1:E_nout) = E_Gamma_core     
	D_neutral_puff(1:5,1:E_nout) = E_neutral_puff
	D_molecule_puff(1:5,1:E_nout) = E_molecule_puff
	D_core_fuelling(1:E_nout) = E_core_fuelling
	D_core_neutral_density(1:E_nout) = E_core_neutral_density
	
	! add pumps as variable inputs? 

	!write(*,*) 'D_Q(:)', D_Q_core
	!write(*,*) 'D_G(:)', D_Gamma_core
	!wite(*,*) 'D_rec', D_rec
	!write(*,*), 'E_rec', E_rec
        istate = 1
        itask = 1 ! for normal computation in dvode_f90 till end_time
	



        write(*,*) 'F starting internal stepping loop'
	start_time = E_start_time        
	do internal_istep = 1, E_nout
        !internal_istep = istep ! pass istep to the physics_routine module
        !global_time = (internal_istep-1)*nout*delta_t
        end_time = start_time + delta_t	
        !write(*,*) 'F call nvt2y on internal_istep = ', internal_istep
	!write(*,*) 'velocity', velocity
        call nvt2ys( Nx, density, velocity, temperature, neutral, neutral_velocity, molecule, ys)
	call nr2yr( extern_neutral_density, extern_molecule_density, yr) 
	call nc2yc( core_density,yc)

 	if( detect_nan ) then	! check if input for integration is NaN  (this is for debugging purposes)
	do ix = 1,Nx
	if (isnan(density(ix))) stop 'step: "density" is a NaN'
	if (isnan(velocity(ix))) stop 'step: "velocity" is a NaN'
	if (isnan(temperature(ix))) stop 'step: "temperature" is a NaN'
	if (isnan(neutral(ix))) stop 'step: "neutral" is a NaN'
	if (isnan(neutral_velocity(ix))) stop 'step: "neutral_velocity" is a NaN'
	if (isnan(molecule(ix))) stop 'step: "molecule" is a NaN'
	enddo
	do ix = 1,5
	if (isnan(extern_neutral_density(ix))) stop 'step: "extern_neutral_density" is a NaN'
	if (isnan(extern_molecule_density(ix))) stop 'step: "extern_molecule_density" is a NaN'
	enddo
	if (isnan(core_density)) stop 'step: "core_density" is a NaN'
	endif
	
        !write(*,*) 'F step internal_istep = ',internal_istep
	!write(*,*) 'F step t0=', start_time, ' t1=', end_time
        !write(*,*) 'F step in T(10) = ', temperature(10)
	!write(*,*) 'nb initialized = ', exp(y(601:605))*1.0d19
	!write(*,*) 'y = ', y
	!write(*,*) 'y(201:300) = ', y(2*Nx+1:3*Nx) 
	!write(*,*) 'T = ', temperature	
	call step_div1d( ys, yr, yc, start_time, end_time, itask, istate, options )     

	! unpack 
        call ys2nvt(Nx, ys, density, velocity, temperature, neutral, neutral_velocity, molecule )
	call yr2nr( yr, extern_neutral_density, extern_molecule_density ) 
	call yc2nc( yc, core_density )
        !write(*,*) 'F step out T(10) =', temperature(10);
	if( detect_nan ) then ! check if outcome of integration is NaN  (this is for debugging purposes)
	do ix = 1,Nx
	if (isnan(density(ix))) stop 'step: "density" is a NaN'
	if (isnan(velocity(ix))) stop 'step: "velocity" is a NaN'
	if (isnan(temperature(ix))) stop 'step: "temperature" is a NaN'
	if (isnan(neutral(ix))) stop 'step: "neutral" is a NaN'
	if (isnan(neutral_velocity(ix))) stop 'step: "neutral_velocity" is a NaN'
	if (isnan(molecule(ix))) stop 'step: "molecule" is a NaN'
	enddo
	do ix = 1,5
	if (isnan(extern_neutral_density(ix))) stop 'step: "extern_neutral_density" is a NaN'
	if (isnan(extern_molecule_density(ix))) stop 'step: "extern_molecule_density" is a NaN'
	enddo
	if (isnan(core_density)) stop 'step: "core_density" is a NaN'
	endif

        start_time = end_time
        enddo
	E_end_time = end_time
	
	!write(*,*) 'F call ELMS'
	call simulate_elm(elm_heat_load, elm_density_change, elm_core_particle_source, end_time)
        !write(*,*) 'F calculate fluxes'
	if( evolve_core .eq. 1 ) then   
	call calculate_core2sol_ion_flux( core_density, Gamma_core2sol) ! this we only do if called from the core 
        else
	Gamma_core2sol = D_Gamma_core(E_nout) ! if the core is not used, use the input for upcoming calculations
	endif
        
	call calculate_extern2core_fluxes(extern2core_flux,extern_neutral_density,core_ext_neutral_pump) ! these should be passed to the core model
	call calculate_extern2core_fluxes(extern2core_mol,extern_molecule_density,core_ext_molecule_pump) ! these should be passed to the core model
	! intermediate calculation to get the neutral shinethrough fluxes
	call calculate_sol_extern_neutral_fluxes(extern2sol_flux, sol2extern_flux, zero3, & !Out
					extern_neutral_density, neutral, neutral_residence_time, extern_neutral_ex )
	call calculate_core_ionization_neutral_flux(extern2core_shinethrough, extern2sol_flux, core_ionization_fraction )
        call calculate_core_ionization_neutral_flux(extern2core_shinethrough_mol, extern2sol_mol, core_ionization_fraction_mol )
 
	! shinethrough flux is considered in the core2sol_flux
	call calculate_sol2core_neutral_flux(core2sol_flux, sol2core_flux, D_core_neutral_density(E_nout), neutral,  core_sol_neutral_ex, extern2core_shinethrough )
	call calculate_sol2core_neutral_flux(core2sol_mol , sol2core_mol , D_core_neutral_density(E_nout), molecule, core_sol_molecule_ex, extern2core_shinethrough_mol)
	
	! molecules do not have a shinethrough factor
   	!call calculate_sol2core_neutral_flux(core2sol_mol , sol2core_mol , D_core_neutral_density(E_nout), molecule, core_sol_molecule_ex, zeroNx)
 	call calculate_sol2extern_ion_flux( sol2extern_ion_flux, sum_sol2extern_ion_flux, density , Gamma_core2sol) ! is depends on the core settings

	call calculate_sol_fluxes( Nx, density, velocity, temperature, neutral, neutral_velocity, molecule, &
				extern_neutral_density,  extern_molecule_density,  &
				Gamma_n, Gamma_mom, q_parallel, Gamma_neutral, Gamma_molecule, Gamma_mom_neutral,& 
				elm_heat_load, & 
                                extern_neutral_flux, extern_molecule_flux, &
				sol2extern_flux, sol2extern_mol, & 
				extern2sol_flux, extern2sol_mol, &
				tar2extern_flux, tar2extern_mol)!, &
!	write(*,*) 'step: neutral_flux =',extern_neutral_flux
!        write(*,*) 'step: molecule =',extern_molecule_flux
        ! correct for the core ionization fluxes
        !call calculate_extern_sol_core_shinethrough_neutral_fluxes(extern2core_st,extern2sol_flux,extern_neutral_density, neutral, neutral_residence_time, core_ionization_fraction)
        !sol2core_flux = sol2core_flux - sum(extern2core_st)

  	!write(*,*) 'F calculate sources'
	call calculate_sol_sources( Nx, density, velocity, temperature, neutral, neutral_velocity, molecule, q_parallel, extern2sol_flux, extern2sol_mol, & 
                                 Source_n, Source_v, Source_Q, Source_neutral, Source_vn, Source_molecule, elm_heat_load, elm_core_particle_source, &
				 Gamma_core2sol, core2sol_flux, core2sol_mol, sol2extern_ion_flux )

	! atoms
	E_sum_sol2extern_ion_flux  = (1.0d+0-mol_rec)*sum_sol2extern_ion_flux

	call calculate_wall_association(atom_association_sink, molecule_association_source, extern_neutral_wall_area, extern_neutral_density,neutral_residence_time,wall_association_probability)
	!call calculate_leakage_fluxes(Gamma_neutral, leakage_fluxes_n)        
	call calculate_extern_sources(sol2extern_flux, tar2extern_flux, extern_neutral_flux, extern2core_flux, E_sum_sol2extern_ion_flux, &
				    Source_extern,atom_association_sink)
	!write(*,*) 'F calc_sour : molecule reservoirs'
        ! molecules
	E_sum_sol2extern_ion_mol = 0.5d+0*mol_rec*sum_sol2extern_ion_flux
	!write(*,*) 'E_sum_sol2extern_ion_mol', E_sum_sol2extern_ion_mol
	!write(*,*) 'sol2extern_mol', sol2extern_mol
        !write(*,*) 'tar2extern_mol', tar2extern_mol
  	!write(*,*) 'extern_molecule_flux', extern_molecule_flux
	!write(*,*) 'extern2core_mol', extern2core_mol
	!call calculate_leakage_fluxes(Gamma_molecule, leakage_fluxes_m)        
        call calculate_extern_sources(sol2extern_mol, tar2extern_mol, extern_molecule_flux, extern2core_mol,E_sum_sol2extern_ion_mol, &
				    Source_extern_mol,molecule_association_source)

	!write(*,*) 'F calc_sour : Pumps'
	! atoms
        call calculate_extern_pump(neutral_pump,extern_neutral_density,pump_rate_n)
	! molecules
	call calculate_extern_pump(molecule_pump,extern_molecule_density,pump_rate_m)
        !write(*,*) 'F calc_sour : core'
        call calculate_core_source(extern2core_flux, extern2core_mol, Gamma_core2sol, sol2core_flux, sol2core_mol, Source_core)
        ! pass variables back through the interface
 	!write(*,*) 'F passback states, fluxes, sources'
	!write(*,*) 'F passback : sol'
	! SOL vars
        E_density       = density
        E_velocity      = velocity
        E_temperature   = temperature
        E_neutral       = neutral
        E_Gamma_n       = Gamma_n
        E_Gamma_mom     = Gamma_mom
        E_q_parallel    = q_parallel
        E_Gamma_neutral  = Gamma_neutral
	E_Gamma_molecule  = Gamma_molecule
        E_Source_n      = Source_n
        E_Source_v      = Source_v
        E_Source_Q      = Source_Q
        E_Source_neutral = Source_neutral
	E_Source_vn 	 = Source_vn
	E_Source_molecule = Source_molecule
	!write(*,*)  'F passback : neutral reservoirs'
	! extern neutral vars
	E_extern_neutral_density = extern_neutral_density
	E_extern_neutral_flux = extern_neutral_flux
	E_sol2extern_flux = sol2extern_flux
	E_extern2sol_flux = extern2sol_flux
	E_tar2extern_flux = tar2extern_flux
	E_Source_extern = Source_extern
	E_neutral_pump = neutral_pump
 	!write(*,*) 'F passback : molecule reservoirs'
	! extern molecule vars
	E_extern_molecule_density = extern_molecule_density
	E_extern_molecule_flux = extern_molecule_flux
	E_sol2extern_mol = sol2extern_mol
	E_extern2sol_mol = extern2sol_mol
	E_tar2extern_mol = tar2extern_mol
	E_Source_extern_mol = Source_extern_mol
	E_molecule_pump = molecule_pump
	!write(*,*) 'F passback : core'
	!write(*,*) 'run-div1d: extern2core_flux', extern2core_flux
	! core vars ! if solved 
	E_core_density = core_density
	E_extern2core_flux = extern2core_flux
	E_extern2core_mol = extern2core_mol
	E_Gamma_core2sol = Gamma_core2sol
	E_sol2core_flux = sol2core_flux
	E_sol2core_mol = sol2core_mol
        E_Source_core = Source_core

        !E_end_time      = end_time
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


end module div1d_step
