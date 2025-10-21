module physics_parameters
! module defining the physics parameters and their default values
! when the div1d_physics namelist is read the normalizations are defined

   use numerics_parameters, only : Nx, density_norm, temperature_norm, velocity_norm, momentum_norm, energy_norm, ntime, delta_t, nout !, noutsteps
   use constants, only : e_charge

   implicit none
	! integer is INTEGER*4 by default and turns into signed int for C
   integer, parameter, private :: wp = KIND(1.0D0) ! double in C
   real( wp ) :: gamma                  = 6.5d+0      ! sheath heat transmission factor [-]
   real( wp ) :: L                      = 5.0d+1      ! total lenght along flux tube (from midpoint to X- and) from X-point to target/sheath [m] (value from Stangeby problem 5.1)
   real( wp ) :: L_core_SOL             = 0.0d+0      ! lenght of core-SOL boundary: i.e. lenght along flux tube from midpoint to X-point [m] (default 0.0 => no core-SOL)
   real( wp ) :: X_core_SOL             = 0.0d+0      ! position of first core-SOL boundary: (default 0.0 => symmetric core SOL, divertor case)
   real( wp ) :: L_baffle               = 0.0d+0      ! Length along flux-tube from Xpoint to baffle position: (default 0.0 => Baffle is located at Xpoint)
   real( wp ) :: sintheta               = 0.1d+0      ! sinus of angle theta between B-field and divertor target plate [-]
   real( wp ) :: mass                   = 3.3436d-27  ! mass of the dominant ion species (default value representing Deuterium) [kg]
   real( wp ) :: Gamma_X           	= 1.0d+23     ! particle flux entering the flux tube at the X-point [/m^2s]
   real( wp ) :: q_parX                 = 1.0d+8      ! parallel heat flux entering the flux tube at the X-point [W/m^2] (value from Stangeby problem 5.1)
   real( wp ) :: Gamma_core		= 1.0d+21     ! particle flux coming from the core boundary when L_core_SOL .gt. 0 (note: Gamma_X is neglected in this case)
   real( wp ) :: Q_core	                = 1.0d+6      ! heat flux coming from the core boundary when L_core_SOL .gt. 0 (note: qparX is neglected in this case)
   real( wp ) :: alpha_core_profile_Q     = 1.0d+0    ! parameter describing the ballooning of the the core losses ~ (1 - (x/L_core_SOL)^2)^alpha
   real( wp ) :: alpha_core_profile_n     = 1.0d+0    ! parameter describing the ballooning of the the core losses ~ (1 - (x/L_core_SOL)^2)^alpha
   logical     :: case_single_null = .false.   ! double null, or single null parameter describing if the core source will have a OMP biased distribution of core source
   real( wp ) :: flux_expansion         = 1.0d+0      ! the flux expansion factor between X-point and target = B_X / B_target = R_target / R_X (set to -1 to load B_field.dat)
   real( wp ) :: flux_expansion_left   = 0.0d+0       ! the flux expansion factor at the left target ( when set to 0, this will mirror on flux_expansion)
   real( wp ) :: trans_expansion 	= 1.0d+0      ! additional transport-induced expansion of heat flux channel between X-point and target ^ ^ (set to -1 to load B_trans.dat)
   real( wp ) :: trans_expansion_left 	= 0.0d+0      ! transport-induced expansion at left target (when set to 0, this will mirror on trans_expansion)
   real( wp ) :: initial_n              = 1.0d+20     ! initial plasma particle density (homogeneous) and density at X-point [/m^3]
   real( wp ) :: initial_v              = 0.0d+0      ! initial plasma velocity (homogeneous) [m/s]
   real( wp ) :: initial_vn             = 0.0d+0      ! initial neutral velocity (homogeneous) [m/s]
   real( wp ) :: initial_T              = 1.0d+2      ! initial plasma temperature (homogeneous) [eV]
   real( wp ) :: initial_a              = 1.0d+14     ! initial neutral density (homogeneous) [/m^3]
   real( wp ) :: initial_m	 	= 1.0d+14     ! initial molecule density (homogeneous) [/m^3]
   real( wp ) :: initial_ncore		= 1.0d+14     ! initial core density [/m^3]
   real( wp ) :: initial_core_neutral   = 1.0d+14     ! initial core neutral density [/m^3]
   real( wp ) :: initial_nb(5) = (/1.0d+17,1.0d+17,4.60d+16,2.8d+18,2.8d+18/)! [m^-3] external neutral background densities
   real( wp ) :: initial_mb(5)= (/1.0d+17,1.0d+17,4.60d+16,2.8d+18,2.8d+18/)! [m^-3] external molecule background densities
   real( wp ) :: density_ramp_rate      = 0.0d+0      ! ramp rate of the density at the X-point [/m^3s]
   real( wp ) :: energy_loss_ion        = 3.0d+1      ! average loss of plasma energy due to ionization [eV]
   real( wp ) :: neutral_energy	= 5d+0        ! average energy of neutrals [eV]
   !real( wp ) :: atomic_energy_reservoir = 0.5d+0    ! average energy of neutrals in reservoirs [eV] 
   !real( wp ) :: molecule_energy_sol 	= 0.5d+0     ! average energy of molecules [eV]
   !real( wp ) :: molecule_energy_reservoir = 0.5d+0  ! average energy of molecules in reservoirs [eV]
   real( wp ) :: recycling              = 0.95d+0     ! fraction of recycled neutrals coming from the target [-]
   real( wp ) :: atom_recycle_energy_fraction = 0.2  ! fraction of the energy retained in the neutral atom when recycled. See stangeby 2000 section 3.1 figure 3.1.	
   real( wp ) :: mol_rec 		= 0.0d0       ! fraction of recycled ions that are released as molecules
   !real( wp ) :: redistributed_fraction = 0.0d+0      ! fraction of recycled neutrals that is evenly redistributed along the SOL [-]
   real( wp ) :: wall_association_probability = 0.0d+0! probability of atoms hitting the wall to associate
   real( wp ) :: neutral_residence_time = 1.0d+20     ! velocity with which neutrals are lost from the SOL [s/m]
   real( wp ) :: molecule_residence_time = 1.0d+20    ! velocity with which molecules are lost from the SOL [s/m]
   real( wp ) :: ato_res_asy = 1.0d+0      ! > 1 asymmetry in neutral exchange (outgoing particles are X times more energetic) -> results in net zero when the density in SOL is lower than outside	
   real( wp ) :: mol_res_asy = 1.0d+0      ! > 1 asymmetry in neutral exchange (outgoing particles are X times more energetic) -> results in net zero when the density in SOL is lower than outside	
   real( wp ) :: minimum_density        = 1.0d+4      ! densities are not allowed to become smaller than this value [/m^3]
   real( wp ) :: maximum_density        = 1.0d+25     ! maximum density allowed in some routines [/m^3]
   real( wp ) :: minimum_temperature    = 1.0d-1      ! the temperature is not allowed to drop below this value [eV]
   !real( wp ) :: rec_loss_fraction	= 0.0d+0      ! Fraction of atoms generated by recombination that is lost to external volumes	
   integer    :: num_impurities         = 5           ! number of impurities in the list (this is not yet dynamic in size)
   real( wp ) :: impurity_concentration(5) =(/0.0d+1,0.0d+1,0.0d+1,0.0d+1,0.0d+1/)  ! the concentration of impurity ions (default = 1%)
   integer    :: impurity_Z(5)             =(/6,0,0,0,0/)          ! the Z value of the impurity used (default = carbon)
   integer    :: switch_imp_distribution= 0           ! switch to use an impurity profile along the leg loaded from "prf_imp_dis.dat"
   real( wp ) :: sol_width_omp          = 3.0d-2      ! width of the core SOL [m]
   real( wp ) :: location_omp		= 0.5d+0      ! parameter describing relative position of the outer midplane on the core SOL [0,1]
   real( wp ) :: major_radius           = 0.9d+0      ! major radius of the core SOL [m] (for now this is constant)
   real( wp ) :: gas_puff_source        = 0.0d+0      ! total particle source from gas puff per flux tube width [/m^2 s]
   real( wp ) :: gas_puff_location      = 0.0d+0      ! location of gas puff along divertor leg [m]
   real( wp ) :: gas_puff_width         = 1.0d+20     ! Gaussian width of effective gas puff source [m?]

   !real( wp ) :: core_ionization_fraction = 0.0d+0    ! fraction of neutrals that enter the sol and pass on to ionize in the core (0,1).
   real( wp ) :: extern_neutral_volumes(5) = (/0.6d+0,0.1d+0,5.0d+0,1.2d+0,0.6d+0/) ![m^3] (1) inner PFR, (2) inner divertor CFR, (3) CFR, (4) outer divertor CFR, (5) outer PFR
   real( wp ) :: extern_neutral_wall_area(5) = (/0.6d+0,0.1d+0,5.0d+0,1.2d+0,0.6d+0/) ! [m^2] (1) inner PFR, (2) inner divertor CFR, (3) CFR, (4) outer divertor CFR, (5) outer PFR
   real( wp ) :: extern_neutral_ex(3) = (/0.0d+0,0.0d+0,0.0d+0 /) ! [s] exchange between neutral volume: (1) = (5->1),(2) = volume(2->3),(3) = volume(3->4) (Clockwise around core) 
   real( wp ) :: extern_molecule_ex(3) = (/0.0d+0,0.0d+0,0.0d+0 /) ! [s] exchange between molecule volume: (1) = (5->1),(2) = volume(2->3),(3) = volume(3->4) (Clockwise around core) 
   real( wp ) :: pol_target_angle(2) 	= (/90d+0,90d+0/) ! poloidal angle at target 1,2 (between 0-180) [degrees] 0 all redistributed recycling flux -> 1 and 4, while with 180 all flux -> chamber 2 and 5. 
   real( wp ) :: puff_rate_neutral(5)	= (/0.1d+0,0.0d+0,0.0d+0,0.0d+0,0.0d+0/) ! puff rate of neutral atoms in background chamber [particles per second]
   real( wp ) :: puff_rate_molecule(5)	= (/0.1d+0,0.0d+0,1.0d+20,0.0d+0,0.0d+0/) ! puff rate of neutral molecules in background chamber [particles per second]
   logical :: switch_dyn_neutral_puff = .false.
   logical :: switch_dyn_molecule_puff = .false.
   real( wp ) :: pump_rate_n(5)	= 0.0d+0 ! pump rate of atoms in background chamber [particles per second]
   real( wp ) :: pump_rate_m(5)	= 1.0d+0 ! pump rate of molecules in background chamber [particles per second]
   real( wp ) :: core_confinement_time  = 0.2d+0      ! confinement time of particles in the core (s)
   real( wp ) :: core_ext_neutral_pump(5) = 0.0d+0    ! The core, pumping external neutrals directly
   real( wp ) :: core_ext_molecule_pump(5) = 0.0d+0    ! The core, pumping external molecules directly
   real( wp ) :: core_far_sol_ion_loss  = 0.0d+0      ! Loss of ions to the far scrape-off layer in the core-sol (where effective flux expansion is not used yet)
   real( wp ) :: core_far_sol_feedthrough = 0.0d+0    ! Direct feedtrhough of core flux to the main chamber (this is not proportional to the electron density in the SOL!!)
   real( wp ) :: core_ionization_fraction = 0.0d+0    ! fraction of particles from external reservoir that ionize in the core.
   real( wp ) :: core_ionization_fraction_mol = 0.0d+0    ! fraction of molecules from external reservoir that ionize in the core
   real( wp ) :: core_sol_neutral_ex   = 0.0d+0      ! core neutral ex between sol and core [m/s]
   real( wp ) :: core_sol_molecule_ex  = 0.0d+0      ! core molecule ex between sol and core [m/s]
   real( wp ) :: core_fuelling		= 1.0d+22     ! core fuelling rate (e.g. pellet injection or NBI) [/s] 
   real( wp ) :: core_volume		= 5	      ! volume of core plasma [m^3]
   integer    :: wide_PFR		= 0	      ! switch to extend PFR neutral volume to to X-point (0 means from target to baffle, 1 from target to X-point
   real(wp)   :: sigma_nb 		= 1.0d-10     ! Spread of sigmoid function for neutral background profile (default value <<1, so nb follows a step function)

   integer    :: elm_start_time         = 0           ! time step (outer step) at which the ELM starts
   integer    :: elm_ramp_time          = 0           ! time (outer step) over which the ELM ramps up
   integer    :: elm_time_between       = 200000000   ! time (outer step) between two ELMs
   real( wp ) :: elm_expelled_heat      = 0.0d+0      ! total heat flux owing to elm, integrated over time (if L_core_SOL==0 with units [J/m^2] otherwise in units [J])
   real( wp ) :: elm_expelled_particles = 0.0d+0      ! total particle flux owing to elm, integrated over time (if L_core_SOL==0 with units [/m^2s] otherwise in units [/s])

   integer    :: switch_elm_heat_flux   = 0           ! turns off (0, default) or on (1) the elm contribution to the heat flux
   integer    :: switch_elm_density     = 0           ! turns off (0, default) or on (1) the elm contribution to the particle flux
   integer    :: switch_elm_series      = 0           ! turns off (0, default) or on (1) the multi-elm sequence
   integer    :: gaussian_elm           = 0           ! switch between gaussian ELM (1, default) and triangular ELM (anything but 1) 
   logical    :: case_AMJUEL            = .true.      ! use collision rates from AMJUEL data base
   integer    :: load_dynamic_impurities = 0	      ! load in time-dependent impurity concentrations
   character*10 :: charge_exchange_model= "AMJUEL"    ! use charge exchange reaction rates from "AMJUEL" data base, "Havlickova", or "Freeman" and Jones
   character*10 :: ionization_model     = "AMJUEL"    ! use ionization rates from "AMJUEL" data base, "Havlickova", or "Freeman" and Jones
   character*10 :: recombination_model  = "AMJUEL"    ! use recombination rates from "AMJUEL" data base, or "Nakazawa" (combining radiative rec. from Gordeev with 3 body rec. from Hinnov et al)

  ! parameters that are changed in time externally and need to be known by all modules
  real(wp), allocatable, dimension(:,:) :: E_imp_con
  real(wp), allocatable, dimension(:) :: E_neu, E_dneu, E_ngb, E_gas, E_rec, E_qpar_x, E_red_frc, E_Q_core, E_Gamma_core
 ! parameters used in fortran simulations of multiple chunks lenght will be ntime=nout*nout_steps
  real( wp ), allocatable :: dyn_nu(:)  ! density of boundary condition 		
  real( wp ), allocatable :: dyn_dnu(:) ! derivative for ODE solver
  real( wp ), allocatable :: dyn_nb(:,:)  ! neutral background instead of initial_a
  real( wp ), allocatable :: dyn_mb(:,:)  ! molecule background instead of initial_a  
  real( wp ), allocatable :: dyn_gas(:) ! gas source [1/m2] [0,->)
  real( wp ), allocatable :: dyn_rec(:)   ! recycling coefficient [0-1]
  real( wp ), allocatable :: dyn_red_frc(:) ! redistributed fraction [0-1]
  real( wp ), allocatable :: dyn_imp_con(:,:) ! impurity concentration [0-1]
  real( wp ), allocatable :: dyn_qparX(:) ! parallel heat flux [0,->)
  real( wp ), allocatable :: dyn_Gamma_core(:) ! option for dynamic core particle flux [#/s]
  real( wp ), allocatable :: dyn_Q_core(:) ! option for dynamic heat flux [J/s]
  real( wp ), allocatable :: prf_imp_con(:,:) ! profile for impurity distribution along the leg (0 -> L/target)
  real( wp ), allocatable :: dyn_neutral_puff(:,:) ! neutral puff rates for external volumes
  real( wp ), allocatable :: dyn_molecule_puff(:,:) ! neutral puff rates for external volumes
  real( wp ), allocatable :: dyn_core_fuelling(:) ! vector for dynamic core fuelling
  real( wp ), allocatable :: dyn_core_neutral_density(:) ! vector holding values for the core neutral density

    ! dynamic arrays used in right_hand_side in physics_routines to make sure no out of bound requests are made these will have length nout
      integer :: internal_istep, D_nout
      real(wp) :: D_start_time, D_delta_t ! start time dynamically given to run_div1d_
      real(wp), allocatable ::  D_imp_con(:,:)
      real(wp), allocatable ::  D_neu(:) 
      real(wp), allocatable ::	D_dneu(:) 
      real(wp), allocatable ::  D_nb(:,:) 
      real(wp), allocatable ::  D_mb(:,:) 
      real(wp), allocatable ::	D_gas(:)  
      real(wp), allocatable ::	D_rec(:)
      real(wp), allocatable ::	D_qpar_x(:)
      real(wp), allocatable ::	D_red_frc(:) 
      real(wp), allocatable ::	D_Q_core(:) 
      real(wp), allocatable ::	D_Gamma_core(:) 
      real(wp), allocatable ::  D_neutral_puff(:,:) 
      real(wp), allocatable ::  D_molecule_puff(:,:) 
      real(wp), allocatable ::  D_core_fuelling(:) 	
      real(wp), allocatable ::  D_core_neutral_density(:)
   
  logical, private :: exists

namelist /div1d_physics/ gamma, L, sintheta, mass, Gamma_X, q_parX, Q_core, Gamma_core, flux_expansion, flux_expansion_left, trans_expansion, trans_expansion_left, initial_n, initial_v, initial_vn, initial_T, initial_a, initial_m, initial_nb, initial_mb, density_ramp_rate , &
                               L_core_SOL, X_core_SOL, alpha_core_profile_Q, alpha_core_profile_n,  & 
                               energy_loss_ion, neutral_energy,  ato_res_asy, mol_res_asy, neutral_residence_time,  molecule_residence_time, recycling, wall_association_probability, atom_recycle_energy_fraction, impurity_concentration, impurity_Z, &
                               case_AMJUEL, charge_exchange_model, ionization_model, recombination_model, &
                               minimum_temperature, minimum_density, maximum_density, gas_puff_source, gas_puff_location, gas_puff_width, &
                               elm_start_time, elm_ramp_time, elm_time_between, elm_expelled_heat, elm_expelled_particles, &
                               switch_elm_density, switch_elm_heat_flux, switch_elm_series, gaussian_elm, pump_rate_n, pump_rate_m, &
			       extern_neutral_volumes, extern_neutral_wall_area, extern_neutral_ex, extern_molecule_ex, pol_target_angle, puff_rate_neutral, puff_rate_molecule, &
			       sol_width_omp, location_omp, major_radius,  switch_imp_distribution, L_baffle, mol_rec, &	
			        core_fuelling, initial_ncore, core_confinement_time, initial_core_neutral, &
			      core_sol_neutral_ex, core_sol_molecule_ex, core_volume, wide_PFR, sigma_nb, core_far_sol_feedthrough, core_far_sol_ion_loss, core_ionization_fraction , core_ionization_fraction_mol,& 
				core_ext_neutral_pump, core_ext_molecule_pump, case_single_null, switch_dyn_neutral_puff,  switch_dyn_molecule_puff

contains

   subroutine read_physics_parameters( error )
      implicit none
      integer :: error, i, z
      integer :: load_dynamic_impurities = 0

      error = 0
      inquire(file = "input.txt", exist=exists)
      if(exists) then 
        open(unit = 1, file = "input.txt")
        read(1, div1d_physics, IOSTAT = error)
      else
        write(*,*)"Could not find or open: INPUT.TXT"
	    stop
      end if
    write(*,*) 'physics read error =', error
	if( error .ne. 0 ) then
	stop 'div1d physics not read'
	endif
	!write(*,*) 'call assign+physic arra'
	!write(*,*) 'ntime', ntime
	! constrain some of the inputs
	wall_association_probability = min(max(0.0d0,wall_association_probability),1.0d0)

        call assign_physics_arrays
    end subroutine read_physics_parameters
 
    subroutine assign_physics_arrays
      implicit none
      integer :: i, z, j
      ! assign physics arrays that are dynamic in time(for main fortran program)
      !num_impurities = 5 ! size(impurity_concentration)
      integer :: load_dynamic_impurities = 0
	
      allocate( dyn_imp_con(5,0:ntime) )
      allocate( dyn_nu(0:ntime) )
      allocate( dyn_dnu(0:ntime) )
      allocate( dyn_nb(5,0:ntime) ) ! should be 5: ntime
      allocate( dyn_mb(5,0:ntime) ) ! should be 5: ntime
      allocate( dyn_gas(0:ntime) )
      allocate( dyn_rec(0:ntime) ) 
      allocate( dyn_qparX(0:ntime) )
      allocate( dyn_Gamma_core(0:ntime) )
      allocate( dyn_Q_core(0:ntime) )
      allocate( dyn_red_frc(0:ntime) ) ! obsolete
      allocate( prf_imp_con(num_impurities,Nx) )
      allocate( dyn_neutral_puff(5,0:ntime) )
      allocate( dyn_molecule_puff(5,0:ntime) )
      allocate( dyn_core_fuelling(0:ntime) )
      allocate( dyn_core_neutral_density(0:ntime) ) 
      ! %%%%%%%%%% read position dependent parameters %%%%%%% !
	!write(*,*) 'arrays allocated'
	if( switch_imp_distribution .eq. 1 ) then
	 open(1, file = 'prf_imp_con.dat', status = 'old')
	 do z = 1,num_impurities	
	  do i = 1,Nx
		read(1,*) prf_imp_con(z,i)
      		prf_imp_con(z,i) = max(prf_imp_con(z,i),0d+0) ! no negative fractions allowed
		!write(*,*) 'prf_imp_con(z,i)', prf_imp_con(z,i)
	  enddo
	 enddo
	 close(1)
	else
	 do z = 1,num_impurities	
	  do i = 1,Nx
		prf_imp_con(z,i) = 1.0d+0
	  enddo
	 enddo
	endif
	!write(*,*) 'prf_imp_con',prf_imp_con
	
      ! %%%%%%%%%%  read time dependent parameters %%%%%%%%%% !
      ! -------- impurity concentration -------!
	do z = 1,num_impurities
		if ( impurity_concentration(z) .eq. -1 ) then 
		load_dynamic_impurities  = 1
		endif
	enddo      

       if ( load_dynamic_impurities .eq. 1 ) then
	open(1, file = 'dyn_imp_con.dat', status = 'old')
      	 do z = 1,num_impurities
          do i = 0,ntime
           read(1,*) dyn_imp_con(z,i) ! (row, column)
	   if ( impurity_concentration(z) .eq. -1 ) then 
           dyn_imp_con(z,i) = min(max(dyn_imp_con(z,i),0.0d+0),1.0d+0)
	   else  ! otherwise do read, but take the fixed value
		dyn_imp_con(z,i) = min(max(impurity_concentration(z),0.0d+0),1.0d+0)
		!write(*,*) 'dyn_imp_con(z,i)', dyn_imp_con(z,i)
		endif
       !  write(*,*) 'read dyn_imp_con.dat =0'
          enddo
	 enddo
       else
	do z = 1,num_impurities
         do i = 0,ntime
            dyn_imp_con(z,i) = min(max(impurity_concentration(z),0.0d+0),1.0d+0)
         enddo
        enddo
       ! write(*,*) 'dimpdt=0'
      endif   
      close(1)

      ! ------- upstream density ------- !
      if (initial_n .lt. 0.0d+0) then
      open(1, file = 'dyn_nu.dat', status = 'old')
       do i =  0,ntime
        read(1,*) dyn_nu(i)
       end do
       close(1)
       !initial_n = dyn_nu(1) ! overwrite initial_n -> this turns it to  absolute input
       ! derivatives for ODE solver
       do i = 1,ntime-1 ! forward difference
       dyn_dnu(i) = min( max( ( dyn_nu(i + 1) - dyn_nu(i) ) / delta_t , -1.0d+34) ,1.0d+34) 
       ! note that prescribing the derivative can result in
       ! integration errors in long simulations
       ! slope going from i to i +  1  nu_t =10^19, delta_t = 10^-6  real is single
       ! limit value to 1d34, below 32 bits of precision 1.7d+38 (wp = kind(1.0))
       end do
       dyn_dnu(ntime) = 0
       write(*,*) "nu.dat read test.", dyn_nu(1)*1.0d-20, dyn_nu(100)*1.0d-20
       !write(*,*) "dnu.dat read test", dnu_t(100)
      else
        do i = 0,ntime
               dyn_nu(i) = initial_n
               dyn_dnu(i) = 0.0d+0
        end do
      !  write(*,*) "dndt=0"  
      endif

    ! ------- neutral background ------- !
      if (initial_nb(1) .lt. 0.0d+0) then
	write(*,*) 'reading dyn_nb.dat'
      open(1, file = 'dyn_nb.dat', status = 'old')

       do i =  0,ntime
        read(1,*) (dyn_nb(j,i), j=1,5)
       end do
       close(1)
      else
        do i = 0,ntime
               dyn_nb(1:5,i) = initial_nb(1:5)
        end do
      !  write(*,*) "dndt=0"  
      endif
    ! ------- molecule background ------- !
!    do z = 1,5
      if (initial_mb(1) .lt. 0.0d+0) then
	write(*,*) 'reading dyn_mb.dat'
      open(1, file = 'dyn_mb.dat', status = 'old')
       do i =  0,ntime
        read(1,*) (dyn_mb(j,i), j=1,5)
       end do
       close(1)
      else
        do i = 0,ntime
               dyn_mb(1:5,i) = initial_mb(1:5)
        end do
      !  write(*,*) "dndt=0"  
      endif

      ! ------- core neutral density ------- !
      if (initial_core_neutral .lt. 0.0d+0) then
      open(1, file = 'dyn_core_neutral_density.dat', status = 'old')
       do i =  0,ntime
        read(1,*) dyn_core_neutral_density(i)
       end do
       close(1)
      else
        do i = 0,ntime
               dyn_core_neutral_density(i) = initial_core_neutral
        end do
      !  write(*,*) "dndt=0"  
      endif

    ! ------- core fuelling ------- !
      if (core_fuelling .eq. 0.0d+0) then
      open(1, file = 'dyn_core_fuelling.dat', status = 'old')
       do i =  0,ntime
        read(1,*) dyn_core_fuelling(i)
       end do
       close(1)
      else
        do i = 0,ntime
               dyn_core_fuelling(i) = core_fuelling
        end do
      !  write(*,*) "dndt=0"  
      endif

    ! ------- neutral puff rate ------- !
      if (puff_rate_neutral(1) .eq. -1.0d+0) then
	write(*,*) 'reading neutral_puff.dat'
      open(1, file = 'neutral_puff.dat', status = 'old')
       do i =  0,ntime
         read(1,*) (dyn_neutral_puff(j,i), j=1,5)
       end do
       close(1)
      else
	write(*,*) 'input: puff rate neutral, ', puff_rate_neutral
        do i = 0,ntime
               dyn_neutral_puff(1:5,i) = puff_rate_neutral
        end do
	write(*,*) 'input: dyn puff rate neutra;(1:5,1) =',dyn_neutral_puff(1:5,1) 
      !  write(*,*) "dndt=0"  
      endif

    ! ------- molecule puff rate ------- !
      if (puff_rate_molecule(1) .eq. -1.0d+0) then
	write(*,*) 'reading molecule_puff.dat'
      open(1, file = 'molecule_puff.dat', status = 'old')
  	!write(*,*)'ntime', ntime
       do i =  0,ntime-1 ! is this strange maybe? 
         read(1,*) (dyn_molecule_puff(j,i), j=1,5)
       end do
       close(1)
      else
        do i = 0,ntime
               dyn_molecule_puff(1:5,i) = puff_rate_molecule
	       !write(*,*) 'dyn-mol-puff:', dyn_molecule_puff(1:5,i)
        end do
      !  write(*,*) "dndt=0"  
      endif

     ! -------- upstream heat flux -----------!
      if (q_parX .lt. 0.0d+0) then
      open(1, file = 'dyn_qpar.dat', status = 'old')

       do i =  0,ntime
        read(1,*) dyn_qparX(i)
       end do
       close(1)
       !q_parX = dyn_qparX(1) ! overwrite q_parX -> turns it to  absolute input
      else
       do i= 0,ntime
        dyn_qparX(i) = q_parX
       end do 
      ! write(*,*) 'dqdt=0'
      endif

      ! -------- core particle flux ------------!
     if (Gamma_core .lt. 0.0d+0) then
	write(*,*) 'reading dyn_gamma_core.dat'
      open(1, file = 'dyn_gamma_core.dat', status = 'old')
       do i =  0,ntime
        read(1,*) dyn_Gamma_core(i)
       end do
       close(1)
       !Gamma_core = dyn_gamma_core(1) ! overwrite Gamma_core -> turns it to  absolute input
      else
       do i= 0,ntime
        dyn_Gamma_core(i) = Gamma_core
       end do 
      ! write(*,*) 'dGdt=0'
      endif
      ! -------- core heat flux ------------!
      if (Q_core .lt. 0.0d+0) then
	write(*,*) 'reading dyn_q_core.dat'  
      open(1, file = 'dyn_q_core.dat', status = 'old')
       do i =  0,ntime
	!write(*,*) i
        read(1,*) dyn_Q_core(i)
       end do
       close(1)
       !Q_core = dyn_Q_core(1) ! overwrite Q_core -> turns it to  absolute input
      else
       do i= 0,ntime
        dyn_Q_core(i) = Q_core
       end do 
      ! write(*,*) 'dGdt=0'
      endif

      ! ------- Recycling ------!
      if (recycling .lt. 0.0d+0) then
	write(*,*) 'reading dyn_rec.dat'  
        open(2, file = 'dyn_rec.dat', status = 'old')
        do i = 0,ntime
         read(2,*) dyn_rec(i)
         dyn_rec(i) = min(max(dyn_rec(i),0.0d+0),1.0d+0)
        end do
        close(2)
       else
        do i = 0,ntime
         dyn_rec(i) = min(max(recycling,0.0d+0),1.0d+0)
        end do
       ! write(*,*)  "dRdt=0"
      endif 

      ! -------- Gas puff -------!
      if (gas_puff_source .lt. 0.0d+0) then
	write(*,*) 'reading dyn_gas.dat'  
        open(3, file = 'dyn_gas.dat', status = 'old')

        do i = 0,ntime

        read(3,*) dyn_gas(i) 
        dyn_gas(i) = max(dyn_gas(i),0.0d+0)
        end do
        close(3)
        write(*,*) "dgdt=1"
      else      
        do i = 0,ntime

        dyn_gas(i) = max(gas_puff_source,0.0d+0) 
        end do
       ! write(*,*) "dgdt=0"  
      endif
 
      ! %%%%%%%%%%%% end read time dependent parameters %%%%%%%% !

      ! correct the desired normalizations
      if( density_norm .eq. 0.0d+0 ) density_norm = 1.0d+19   !initial_n
      if( temperature_norm .eq. 0.0d+0 ) temperature_norm = 1.0e+0
      if( velocity_norm .eq. 0.0d+0 ) velocity_norm = sqrt( 2.0d+0 * temperature_norm / mass )
      momentum_norm = mass * density_norm * velocity_norm
      energy_norm = density_norm * e_charge * temperature_norm
      !write(*,*) 'phys_par: energy_norm', energy_norm
      return
   end subroutine assign_physics_arrays

   subroutine extern_read_physics_parameters(floatinphys, intinphys) !, strinphys, loginphys)
        implicit none
         real(wp), intent(in) ::  floatinphys(141)
         integer, intent(in) :: intinphys(31)
           
        !logical, INTENT(IN) :: loginphys(1)

          ! floats
	  ! geometry, grid, profiles
          L                      = floatinphys(1)   
          L_core_SOL             = floatinphys(2)      
          X_core_SOL             = floatinphys(3) 
          sintheta               = floatinphys(4)
	  pol_target_angle(1:2)	= floatinphys(5:6) ! poloidal angle at target 1,2 (between 0-180) [degrees] 0 all redistributed recycling flux -> 1 and 4, while with 90 all flux -> chamber 2 and 5. 
          sol_width_omp          = floatinphys(7) !3.0d-2      ! width of the core SOL [m]
          location_omp		= floatinphys(8) !0.5d+0      ! parameter describing relative position of the outer midplane on the core SOL [0,1]
          major_radius           = floatinphys(9) !0.9d+0      ! major radius of the core SOL [m] (for now this is constant)
	  alpha_core_profile_Q   = floatinphys(10)
	  alpha_core_profile_n   = floatinphys(11)
          flux_expansion         = floatinphys(12)
	  trans_expansion         = floatinphys(13)
          gas_puff_location      = floatinphys(14)     
          gas_puff_width         = floatinphys(15) 
          sigma_nb 		=  floatinphys(16)      ! Spread of sigmoid function for neutral background profile
          L_baffle 		= floatinphys(17)
	  ! sol influx settings
          Gamma_X                = floatinphys(20)  
          q_parX                 = floatinphys(21)    
	  Gamma_core 		 = floatinphys(22)
 	  Q_core 		 = floatinphys(23)
          density_ramp_rate      = floatinphys(24)   
 	  gas_puff_source        = floatinphys(25)    

	  ! initial values
          initial_n              = floatinphys(30) 
          initial_v              = floatinphys(31)
          initial_T              = floatinphys(32)  
          initial_a              = floatinphys(33) 
	  initial_vn 	 	 = floatinphys(34) 
          initial_m 		 = floatinphys(35)
	  initial_ncore 	 = floatinphys(36)
	  initial_core_neutral   = floatinphys(37)
	  initial_nb		 = floatinphys(40:44)
	  initial_mb 		 = floatinphys(45:49)

	  ! limits
          minimum_density        = floatinphys(50) 
	  maximum_density        = floatinphys(51)
          minimum_temperature    = floatinphys(52) 

	  ! sol behavior
	  gamma                  = floatinphys(60)
 	  mass                   = floatinphys(61)
          energy_loss_ion        = floatinphys(62)  ! obsolete with amjuel?
          recycling              = floatinphys(63)  
	  mol_rec 		= floatinphys(64)     ! fraction of recycled ions that are released as molecules
          neutral_residence_time = floatinphys(65) 
	  molecule_residence_time = floatinphys(66)
          core_far_sol_ion_loss  =  floatinphys(67)      ! Loss of ions to the far scrape-off layer in the core-sol
          impurity_concentration(1:5) = floatinphys(70:74)
	 
	  ! reservoirs
	  core_ext_neutral_pump(1:5) = floatinphys(80:84) ! core pumping neutrals from external reservoirs
          core_ext_molecule_pump(1:5) = floatinphys(85:89) ! core pumping molecules from external reservoirs
          extern_neutral_volumes(1:5) = floatinphys(90:94)         	
	  extern_neutral_ex(1:3) = floatinphys(100:102)
	  extern_molecule_ex(1:3) = floatinphys(103:105)
  	  pump_rate_n(1:5)	= floatinphys(110:114)  ! pump rate of atoms in background chamber [particles per second]
   	  pump_rate_m(1:5)	= floatinphys(115:119) ! pump rate of molecules in background chamber [particles per second]

	  ! core
          core_confinement_time  = floatinphys(130)     ! confinement time of particles in the core (s)
    	  core_sol_neutral_ex   =  floatinphys(131)       ! core neutral ex_time between sol and core
   	  core_sol_molecule_ex  =  floatinphys(132)       ! core molecule ex_time between sol and core
    	  core_volume		=  floatinphys(133)    ! volume of core plasma [m^3]
 	 

          ! integers
          num_impurities         =5 !intinphys(1) ! 5 ! number of impurities
          impurity_Z(1:5)  = intinphys(2:6) 
          switch_imp_distribution= intinphys(7)
	  wide_PFR		= intinphys(8)	      ! switch to extend PFR neutral volume to to X-point (0 means from target to baffle, 1 from target to X-point
     
         ! case_AMJUEL            = loginphys(1)    ! use collision rates from AMJUEL data base
         ! matlab always uses the default AMJUEL-like settings
        call assign_physics_arrays 
   end subroutine extern_read_physics_parameters        
   
end module physics_parameters
