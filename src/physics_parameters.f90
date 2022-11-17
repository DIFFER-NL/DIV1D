module physics_parameters
! module defining the physics parameters and their default values
! when the div1d_physics namelist is read the normalizations are defined

   use numerics_parameters, only : Nx, density_norm, temperature_norm, velocity_norm, momentum_norm, energy_norm, ntime, delta_t
   use constants, only : e_charge

   implicit none

   integer, parameter, private :: wp = KIND(1.0D0)
   real( wp ) :: gamma                  = 6.5d+0      ! sheath heat transmission factor [-]
   real( wp ) :: L                      = 5.0d+1      ! total lenght along flux tube (from midpoint to X- and) from X-point to target/sheath [m] (value from Stangeby problem 5.1)
   real( wp ) :: L_core_SOL             = 0.0d+0      ! lenght of core-SOL boundary: i.e. lenght along flux tube from midpoint to X-point [m] (default 0.0 => no core-SOL)
   real( wp ) :: X_core_SOL             = 0.0d+0      ! position of first core-SOL boundary: (default 0.0 => symmetric core SOL, divertor case)
   real( wp ) :: sintheta               = 0.1d+0      ! sinus of angle theta between B-field and divertor target plate [-]
   real( wp ) :: mass                   = 3.3436d-27  ! mass of the dominant ion species (default value representing Deuterium) [kg]
   real( wp ) :: Gamma_X                = 1.0d+23     ! particle flux entering the flux tube at the X-point [/m^2s]
   real( wp ) :: q_parX                 = 1.0d+8      ! parallel heat flux entering the flux tube at the X-point [W/m^2] (value from Stangeby problem 5.1)
   real( wp ) :: alpha_core_profile     = 1.0d+0      ! parameter describing the ballooning of the the core losses ~ (1 - (x/L_core_SOL)^2)^alpha
   real( wp ) :: normalization_core_profile = 0.0d+0  ! normalization factor of the loss profile of heat and particles across the core-SOL boundary
   real( wp ) :: flux_expansion         = 1.0d+0      ! the flux expansion factor between X-point and target = B_X / B_target = R_target / R_X
   real( wp ) :: initial_n              = 1.0d+20     ! initial plasma particle density (homogeneous) and density at X-point [/m^3]
   real( wp ) :: initial_v              = 0.0d+0      ! initial plasma velocity (homogeneous) [m/s]
   real( wp ) :: initial_T              = 1.0d+2      ! initial plasma temperature (homogeneous) [eV]
   real( wp ) :: initial_a              = 0.0d+4      ! initial neutral density (homogeneous) [/m^3]
   real( wp ) :: density_ramp_rate      = 0.0d+0      ! ramp rate of the density at the X-point [/m^3s]
   real( wp ) :: energy_loss_ion        = 3.0d+1      ! average loss of plasma energy due to ionization [eV]
   real( wp ) :: recycling              = 0.95d+0     ! fraction of recycled neutrals coming from the target [-]
   real( wp ) :: redistributed_fraction = 0.8d+0      ! fraction of recycled neutrals that is evenly redistributed along the SOL [-]
   real( wp ) :: neutral_residence_time = 1.0d+20     ! time scale on which neutrals are lost from the SOL [s]
   real( wp ) :: minimum_density        = 1.0d+4      ! densities are not allowed to become smaller than this value [/m^3]
   real( wp ) :: minimum_temperature    = 1.0d-1      ! the temperature is not allowed to drop below this value [eV]
   integer    :: num_impurities         = 5           ! number of impurities in the list (this is not yet dynamic in size)
   real( wp ) :: impurity_concentration(5) =(/0.0d+1,0.0d+1,0.0d+1,0.0d+1,0.0d+1/)  ! the concentration of impurity ions (default = 1%)
   integer    :: impurity_Z(5)             =(/6,0,0,0,0/)          ! the Z value of the impurity used (default = carbon)
   real( wp ) :: gas_puff_source        = 0.0d+0      ! total particle source from gas puff per flux tube width [/m^2 s]
   real( wp ) :: gas_puff_location      = 0.0d+0      ! location of gas puff along divertor leg [m]
   real( wp ) :: gas_puff_width         = 1.0d+20     ! Gaussian width of effective gas puff source [m?]
   integer    :: elm_start_time         = 0           ! time step (outer step) at which the ELM starts
   integer    :: elm_ramp_time          = 0           ! time (outer step) over which the ELM ramps up
   integer    :: elm_time_between       = 200000000   ! time (outer step) between two ELMs
   real( wp ) :: elm_expelled_heat      = 0.0d+0      ! total heat flux owing to elm, integrated over time [J/m^2]
   real( wp ) :: elm_expelled_particles = 0.0d+0      ! total particle flux owing to elm, integrated over time [/m^2]
   integer    :: switch_elm_heat_flux   = 0           ! turns off (0, default) or on (1) the elm contribution to the heat flux
   integer    :: switch_elm_density     = 0           ! turns off (0, default) or on (1) the elm contribution to the particle flux
   integer    :: switch_elm_series      = 0           ! turns off (0, default) or on (1) the multi-elm sequence
   integer    :: gaussian_elm           = 0           ! switch between gaussian ELM (1, default) and triangular ELM (anything but 1) 
   logical    :: case_AMJUEL            = .true.      ! use collision rates from AMJUEL data base
   character*10 :: charge_exchange_model= "AMJUEL"    ! use charge exchange reaction rates from "AMJUEL" data base, "Havlickova", or "Freeman" and Jones
   character*10 :: ionization_model     = "AMJUEL"    ! use ionization rates from "AMJUEL" data base, "Havlickova", or "Freeman" and Jones
   character*10 :: recombination_model  = "AMJUEL"    ! use recombination rates from "AMJUEL" data base, or "Nakazawa" (combining radiative rec. from Gordeev with 3 body rec. from Hinnov et al)
   real( wp ) :: radial_loss_factor     = 0           ! percentage of the parallel flux that is lost radially throughout the flux tube (not exact for radial_loss_gaussian = -1)
   integer    :: radial_loss_gaussian   = 0           ! set to 0 (default) for a constant loss factor, to 1 for a gaussian distribution or to -1 for a locally dependent version 
   real( wp ) :: radial_loss_width      = 1.0d+20       ! determine width of radial loss distribution (only used for radial_loss_gaussian = 1) [m]
   real( wp ) :: radial_loss_location   = 0           ! determine peak location of radial loss distribution (only used for radial_loss_gaussian = 1) [m]


!  time dependent settings
 !  integer    :: switch_dyn_nu          = 0           ! switch now depends on initial_n value .leq. -1 !  time dependent plasma X point density requested from dyn_nu.dat (perturbation on initial_n)
 !  integer    :: switch_dyn_gas         = 0           ! time dependent gas source quested from dyn_gas.dat [/m^2 s]
 !  integer    :: switch_dyn_rec         = 0           ! time dependent recycling fraction taken from dyn_rec.dat [-]
 !  integer    :: switch_dyn_rad_los     = 0           ! time dependent radial loss factor taken from dyn_rad_loss.dat
 !  integer    :: switch_dyn_imp_con     = 0           ! time dependent impurity concentration from dyn_imp_con.dat   
 !  integer    :: switch_dyn_qpar        = 0           ! time dependent qparallel boundary condition from dyn_qpar.dat 
 !  integer    :: switch_dyn_red_frc     = 0           ! time dependent redistribution fraction from dyn_red_frc.dat


  real( wp ), allocatable :: dyn_nu(:)  ! density of boundary condition 		
  real( wp ), allocatable :: dyn_dnu(:) ! derivative for ODE solver
  real( wp ), allocatable :: dyn_nb(:)  ! neutral background instead of initial_a
  real( wp ), allocatable :: dyn_gas(:) ! gas source [1/m2] [0,->)
  real( wp ), allocatable :: dyn_rec(:)   ! recycling coefficient [0-1]
  real( wp ), allocatable :: dyn_red_frc(:) ! redistributed fraction [0-1]
  real( wp ), allocatable :: dyn_rad_los(:)  ! radial loss factor [0-1]
  real( wp ), allocatable :: dyn_imp_con(:,:) ! impurity concentration [0-1]
  real( wp ), allocatable :: gas_puff(:) ! gas puff distribution ( this is now globally accessable )
  real( wp ), allocatable :: dyn_qparX(:) ! parallel heat flux [0,->)


contains

   subroutine read_physics_parameters( error )
      implicit none
      integer :: error, i, num_impurities, z

    !  integer, parameter, private :: wp = KIND(1.0D0)
    !  real( wp ) :: tmp_imp
      num_impurities = 0

      namelist /div1d_physics/ gamma, L, sintheta, mass, Gamma_X, q_parX, flux_expansion, initial_n, initial_v, initial_T, initial_a, density_ramp_rate, &
                               L_core_SOL, X_core_SOL, alpha_core_profile, &
                               energy_loss_ion, neutral_residence_time, redistributed_fraction, recycling, num_impurities, impurity_concentration, impurity_Z, &
                               case_AMJUEL, charge_exchange_model, ionization_model, recombination_model, &
                               minimum_temperature, minimum_density, gas_puff_source, gas_puff_location, gas_puff_width, &
                               elm_start_time, elm_ramp_time, elm_time_between, elm_expelled_heat, elm_expelled_particles, &
                               switch_elm_density, switch_elm_heat_flux, switch_elm_series, gaussian_elm, &
                               radial_loss_factor, radial_loss_gaussian, radial_loss_width, radial_loss_location            
                              !density_ramp_rate
      error = 0
!      open(, file = 'input.txt', status = 'old')
      read(*,div1d_physics, IOSTAT = error)
      write(*,*) 'physics read error =', error
!      close(1) 
        
      num_impurities = 5 ! size(impurity_concentration)
      allocate( dyn_imp_con(num_impurities,ntime) )
      allocate( dyn_nu(ntime) )
      allocate( dyn_nb(ntime) )
      allocate( dyn_dnu(ntime) )
      allocate( dyn_gas(ntime) )
      allocate( dyn_rec(ntime) ) 
      allocate( dyn_rad_los(ntime) )
      allocate( dyn_qparX(ntime) )
      allocate( dyn_red_frc(ntime) )



      ! %%%%%%%%%%  read time dependent parameters %%%%%%%%%% !
      ! -------- impurity concentration -------!
      do z = 1,num_impurities
      if ( impurity_concentration(z) .eq. -1 ) then
        open(1, file = 'dyn_imp_con.dat', status = 'old')
        do i = 1,ntime
          read(1,*) dyn_imp_con(z,i) ! (row, column)
         dyn_imp_con(z,i) = min(max(dyn_imp_con(z,i),0.0d+0),1.0d+0)
       !  write(*,*) 'read dyn_imp_con.dat =0'
         end do
          close(1)
       else
        do i = 1,ntime
            dyn_imp_con(z,i) = min(max(impurity_concentration(z),0.0d+0),1.0d+0)
            end do
       ! write(*,*) 'dimpdt=0'
      endif
      enddo
      ! ------- upstream density ------- !
      if (initial_n .lt. 0.0d+0) then
      open(1, file = 'dyn_nu.dat', status = 'old')
       do i =  1,ntime
        read(1,*) dyn_nu(i)
       end do
       close(1)
       !initial_n = dyn_nu(1) ! overwrite initial_n -> this turns it to  absolute input
       ! derivatives for ODE solver
       do i = 1,ntime-1 ! forward difference
       dyn_dnu(i) = min( ( dyn_nu(i + 1) - dyn_nu(i) ) / delta_t ,1.0d+34) ! note that prescribing the derivative can result in
                                                                           ! integration errors in long simulations
       ! slope going from i to i +  1  nu_t =10^19, delta_t = 10^-6  real is single
       !  limit value to 1d34, below 32 bits of precision 1.7d+38 (wp = kind(1.0))
       end do
       dyn_dnu(ntime) = 0
       !write(*,*) "nu.dat read test.", nu_t(1010), nu_t(1011)
       !write(*,*) "dnu.dat read test", dnu_t(1010)
      else
        do i = 1,ntime
               dyn_nu(i) = initial_n
               dyn_dnu(i) = 0.0d+0
        end do
      !  write(*,*) "dndt=0"  
      endif

    ! ------- neutral background ------- !
      if (initial_a .lt. 0.0d+0) then
      open(1, file = 'dyn_nb.dat', status = 'old')
       do i =  1,ntime
        read(1,*) dyn_nb(i)
       end do
       close(1)
       !initial_a = dyn_nb(1) ! overwrite initial_a -> this turns it to  absolute input
       ! derivatives for ODE solver
      else
        do i = 1,ntime
               dyn_nb(i) = initial_a
        end do
      !  write(*,*) "dndt=0"  
      endif


     ! -------- upstream heat flux -----------!
      if (q_parX .lt. 0.0d+0) then
      open(1, file = 'dyn_qpar.dat', status = 'old')
       do i =  1,ntime
        read(1,*) dyn_qparX(i)
       end do
       close(1)
       !q_parX = dyn_qparX(1) ! overwrite q_parX -> turns it to  absolute input
      else
       do i= 1,ntime
        dyn_qparX(i) = q_parX
       end do 
      ! write(*,*) 'dqdt=0'
      endif
  
      ! ------- Recycling ------!
      if (recycling .lt. 0.0d+0) then
        open(2, file = 'dyn_rec.dat', status = 'old')
        do i = 1,ntime
         read(2,*) dyn_rec(i)
         dyn_rec(i) = min(max(dyn_rec(i),0.0d+0),1.0d+0)
        end do
        close(2)
       else
        do i = 1,ntime
         dyn_rec(i) = min(max(recycling,0.0d+0),1.0d+0)
        end do
       ! write(*,*)  "dRdt=0"
      endif 

      ! -------- Gas puff -------!
      if (gas_puff_source .lt. 0.0d+0) then
        open(3, file = 'dyn_gas.dat', status = 'old')
        do i = 1,ntime
        read(3,*) dyn_gas(i) 
        dyn_gas(i) = max(dyn_gas(i),0.0d+0)
        end do
        close(3)
        write(*,*) "dgdt=1"
      else      
        do i = 1,ntime
        dyn_gas(i) = max(gas_puff_source,0.0d+0) 
        end do
       ! write(*,*) "dgdt=0"  
      endif

      ! -------- Radial Loss fraction ----- !
      if (radial_loss_factor .lt. 0.0d+0) then
        open(4, file = 'dyn_rad_los.dat', status= 'old')
        do i = 1,ntime
        read(4,*) dyn_rad_los(i) 
        dyn_rad_los(i) = min(max(dyn_rad_los(i),0.0d+0),1.0d+0)
        end do
        close(4)
        write(*,*) "dRLdt=1"
      else
        do i = 1,ntime
        dyn_rad_los(i) = min(max(radial_loss_factor,0.0d+0),1.0d+0)
        end do
       ! write(*,*) "dRLdt=0"
      endif
      
      
      ! ------------- redistribution fraction --------------! 
      if (redistributed_fraction .lt. 0.0d+0) then
        open(4, file = 'dyn_red_frc.dat', status= 'old')
        do i = 1,ntime
        read(4,*) dyn_red_frc(i) 
        dyn_red_frc(i) = min(max(dyn_red_frc(i),0.0d+0),1.0d+0)
        end do
        close(4)
        write(*,*) "dfRLdt=1"
      else
        do i = 1,ntime
        dyn_red_frc(i) = min(max(redistributed_fraction,0.0d+0),1.0d+0)
        end do
       ! write(*,*) "dfRLdt=0"
      endif
      ! %%%%%%%%%%%% end read time dependent parameters %%%%%%%% !

      ! correct the desired normalizations
      if( density_norm .eq. 0.0d+0 ) density_norm = initial_n ! was 10**19
      if( temperature_norm .eq. 0.0d+0 ) temperature_norm = 1.0d+0
      if( velocity_norm .eq. 0.0d+0 ) velocity_norm = sqrt( 2.0d+0 * temperature_norm / mass )
      momentum_norm = mass * density_norm * velocity_norm
      energy_norm = density_norm * e_charge * temperature_norm
      
      return
   end subroutine read_physics_parameters
   
end module physics_parameters
