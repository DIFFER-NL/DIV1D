module plasma_data
! module defining and handling the plasma data

   use numerics_parameters, only : Nx
   use grid_data, only           : x
   use physics_parameters, only  : L, dyn_nu, initial_v, initial_T, dyn_nb, dyn_mb, initial_n, initial_a, initial_m, initial_vn, initial_nb, initial_mb, initial_ncore, mass, dyn_qparX, Gamma_X, gamma, density_norm, velocity_norm, temperature_norm 
   ! note dyn_nu, dyn_nb, and dyn_qpar replaced initial_n and initial_a and q_parX
   use physics_routines
   use interpolation

   implicit none
   integer, parameter, private :: wp = KIND(1.0D0)
   ! plasma
   integer ::                 neq             ! number of equation (= 3 Nx (= 4 Nx when adding the neutrals))
   real( wp ), allocatable :: y(:)            ! vector holding the solution for the plasma variables [various units] 
   real( wp ), allocatable :: ys(:)            ! vector holding the solution for the plasma variables [various units] 
   real( wp ), allocatable :: ydot(:)         ! vector holding the solution for the plasma variables [various units]
   real( wp ), allocatable :: ynb(:) 	      ! vector holding the solution for the background densities [m^-3]
   real( wp ), allocatable :: density(:)      ! vector holding the solution for the plasma density [/m^3]
   real( wp ), allocatable :: velocity(:)     ! vector holding the solution for the plasma velocity [m/s]
   real( wp ), allocatable :: temperature(:)  ! vector holding the solution for the plasma temperature [eV]
   real( wp ), allocatable :: neutral(:)      ! vector holding the solution for the neutral density [/m^3]
   real( wp ), allocatable :: molecule(:)      ! vector holding the solution for the molecule density [/m^3]
   real( wp ), allocatable :: neutral_velocity(:)      ! vector holding the solution for the neutral velocity [/m^3]
   real( wp ), allocatable :: Gamma_n(:)      ! vector holding the solution for the plasma particle flux [/m^2s]
   real( wp ), allocatable :: Gamma_mom(:)    ! vector holding the solution for the plasma momentum flux [Pa]
   real( wp ), allocatable :: Gamma_mom_neutral(:)    ! vector holding the solution for the plasma momentum flux [Pa]
   real( wp ), allocatable :: pressure(:)     ! vector holding the solution for the plasma pressure [Pa]
   real( wp ), allocatable :: q_parallel(:)   ! vector holding the solution for the plasma heat flux [W/m^2]
   real( wp ), allocatable :: Gamma_neutral(:) ! vector holding the solution for the plasma particle flux [/m^2s]
   real( wp ), allocatable :: Gamma_molecule(:) ! vector holding the solution for the molecule flux [/m^2s]
   real( wp ), allocatable :: extern2sol_flux(:) ! vector holding the fluxes from external neutral volumes into the SOL [1/s] (Nx)
   real( wp ), allocatable :: core2sol_flux(:) ! vector holding the fluxes from the core into the SOL [1/s] (Nx)
   real( wp ), allocatable :: extern2sol_mol(:) ! vector holding the fluxes from external neutral volumes into the SOL [1/s] (Nx)
   real( wp ), allocatable :: core2sol_mol(:) ! vector holding the fluxes from the core into the SOL [1/s] (Nx)
   real( wp ), allocatable :: sol2extern_ion_flux(:)
   real( wp ), allocatable :: extern2core_shinethrough_flux(:) ! vector holding the particle flux passing the SOL into the core.
   real( wp ), allocatable :: Source_n(:)     ! vector holding the solution for the plasma particle source [/m^3s]
   real( wp ), allocatable :: Source_v(:)     ! vector holding the solution for the plasma momentume source [?]
   real( wp ), allocatable :: Source_Q(:)     ! vector holding the solution for the plasma heat source [W/m^3]
   real( wp ), allocatable :: Source_neutral(:)     ! vector holding the solution for the plasma heat source [/m^3s]
   real( wp ), allocatable :: Source_molecule(:)     ! vector holding the solution for the molecule source [/m^3s]
   real( wp ), allocatable :: Source_vn(:)     ! vector holding the solution for the neutral momentum source [(kg m/s)/m^3/s]

   ! reservoirs
   real( wp ) :: yr(10) 		      ! vector for solutions of the reservoir
   ! atoms
   real( wp )  :: extern_neutral_density(5)     ! vector holding the solution for the external neutral densities [m^-3]
   real( wp )  :: extern_neutral_flux(3) = (/0.0d+0,0.0d+0,0.0d+0/) ! vector with neutral fluxes between external volumes [1/s]
   real( wp )  :: sol2extern_flux(5) = (/0.0d+0,0.0d+0,0.0d+0,0.0d+0,0.0d+0/) ! vector with neutral fluxes from SOL into external volumes [1/s]
   real( wp )  :: tar2extern_flux(5) 	= 0.0d+0 ! redistributed target recycling fluxes.
   real( wp )  :: extern2core_flux(5) = (/0.0d+0,0.0d+0,0.0d+0,0.0d+0,0.0d+0/)  !fluxes from neutral reservoirs directly to core
   real( wp )  :: sum_sol2extern_ion_flux = 0.0d+0
   real( wp )  :: neutral_pump(5) = 0.0d+0 ! neutral pumping term
   real( wp )  :: Source_extern(5) = 0.0d+0 ! the source terms for the external neutral volumes
   ! molecules
   real( wp )  :: extern_molecule_density(5)     ! vector holding the solution for the external molecule densities [m^-3]
   real( wp )  :: extern_molecule_flux(3) = (/0.0d+0,0.0d+0,0.0d+0/) ! vector with neutral fluxes between external volumes [1/s]
   real( wp )  :: sol2extern_mol(5) = (/0.0d+0,0.0d+0,0.0d+0,0.0d+0,0.0d+0/) ! molecule fluxes from SOL into external volumes [1/s]
   real( wp )  :: tar2extern_mol(5)     = 0.0d+0 ! redistributed target recycling fluxes.
   real( wp ) :: extern2core_mol(5)= (/0.0d+0,0.0d+0,0.0d+0,0.0d+0,0.0d+0/)  ! fluxes from molecule reservoirs directly to core
   real( wp ) :: sum_sol2extern_ion_mol = 0.0d+0 
   real( wp )  :: molecule_pump(5) = 0.0d+0 ! molecule pumping term
   real( wp )  :: Source_extern_mol(5)  = 0.0d+0 ! the source termsfor the external neutral volumes

   ! core
   real( wp )  :: core_density = 0.0d+0   ! Core density
   real( wp )  :: core_neutral_density = 0.0d+0 ! neutral density in the core to calculate inflow of neutrals
   real( wp )  :: Gamma_core2sol = 0.0d+0 ! Source of particles coming from core [/s]
   real( wp )  :: Source_core = 0.0d+0 ! core source

   ! legacy but might re-introduce this
   real( wp )  :: sol2core_flux = 0.0d+0 ! flux from SOL into core.  
   real( wp )  :: sol2core_mol = 0.0d+0 ! molecule flux from SOL into core. 

! make  namelists to write solutions easier and in a more general way? or maybe start making structs?
!namelist /data_sol_profiles/ density, velocity, temperature, neutral, molecule, neutral_velocity 
!namelist /data_sol_fluxes/  Gamma_n, Gamma_mom, q_parallel, Gamma_mom_neutral, neutral_flux, molecule_flux
!namelist /data_sol_sources/  Source_n, Source_v, Source_Q, Source_vn, Source_neutral, Source_molecule
!namelist /data_sol_extern/ extern2sol_flux, extern2sol_mol, core2sol_flux, core2sol_mol, Gamma_core2sol
!namelist /data_extern/ extern_neutral_density, extern_molecule_density, 
!namelist /data_sol_extern_fluxes1/sol2extern_flux, sol2extern_mol, tar2extern_flux, tar2extern_mol
!namelist /data_extern_sources/ Source_extern_atom, Source_exern_mol
!namelist /data_extern_fluxes2/ extern_neutral_flux, extern_molecule_flux
contains

 
   subroutine initial_values
   ! subroutine to initialize the plasma data and solution vetor
   ! the solution vector is defined as follows
   !    y(      1, ...,   Nx )     = density( 1, ..., Nx )
   !    y(   Nx+1, ..., 2*Nx )     = Momentum = mass * density * velocity ( 1, ..., Nx )
   !    y( 2*Nx+1, ..., 3*Nx )     = pressure = 2 * density * e_charge * temperature [eV] ( 1, ..., Nx )
   !    y( 3*Nx+1, ..., 4*Nx )     = neutral density( 1, ..., Nx )
      implicit none
      integer :: neq_tmp 
      real( wp ) :: tmp_density(Nx)      ! vector holding the solution for the plasma density [/m^3]
      real( wp ) :: tmp_velocity(Nx)     ! vector holding the solution for the plasma velocity [m/s]
      real( wp ) :: tmp_temperature(Nx)  ! vector holding the solution for the plasma temperature [eV]
      real( wp ) :: tmp_neutral(Nx)      ! vector holding the solution for the neutral density [/m^3]
      real( wp ) :: tmp_molecule(Nx)      ! vector holding the solution for the molecule density [/m^3]
      real( wp ) :: tmp_neutral_velocity(Nx)      ! vector holding the solution for the neutral velocity [/m^3]
      real( wp ) :: tmp_q_parallel(Nx)   ! vector holding the solution for the plasma heat flux [W/m^2]
      real( wp ) :: tmp_extern_neutral_density(5)     ! vector holding the solution for the external neutral densities [m^-3]
      real( wp ) :: tmp_extern_molecule_density(5)     ! vector holding the solution for the external molecule densities [m^-3]
 
      ! Arrays are allocated in div1d_step.f90 
      temperature = initial_T
      write(*,*) 'F initial value: temperature cell 30', temperature(30)
      density     = initial_n ! dyn_nu(1)
      write(*,*) 'F initial value: density cell 30', density(30)
      velocity    = initial_v
      neutral     = initial_a ! initial_a
      molecule	  = initial_m ! initial_m 
      pressure    = 2.0d+0 * density * temperature
      neutral_velocity = initial_vn
      extern_neutral_density = initial_nb
      extern_molecule_density = initial_mb
      core_density = initial_ncore
          
      ! transform to normalized solution vector y
      call nvt2ys( Nx, density, velocity, temperature, neutral, neutral_velocity, molecule, ys )
      call nr2yr( extern_neutral_density, extern_molecule_density, yr )
      ! transform back to ytmp
      call ys2nvt( Nx, ys, tmp_density, tmp_velocity, tmp_temperature, tmp_neutral, tmp_neutral_velocity, tmp_molecule )
      call yr2nr( yr, tmp_extern_neutral_density, tmp_extern_molecule_density)
      ! check if y2nvt and nvt2y commute
      !write(*,*) 'plsma data: diff n - nt', (density - tmp_density) / density_norm
      !write(*,*) 'plsma data: diff T - Tt', (temperature - tmp_temperature) / temperature_norm
      !write(*,*) 'plsma data: diff v - vt', (velocity - tmp_velocity) / velocity_norm
      !write(*,*) 'plsma data: diff a - at', (neutral - tmp_neutral) / density_norm
      !write(*,*) 'plsma data: diff m - mt', (molecule - tmp_molecule) / density_norm
      !write(*,*) 'plsma data: diff va - vat', (neutral_velocity - tmp_neutral_velocity) / velocity_norm
      !write(*,*) 'plsma data: diff ea - eat', (extern_neutral_density - tmp_extern_neutral_density) / density_norm
      !write(*,*) 'plsma data: diff em - emt', (extern_molecule_density - tmp_extern_molecule_density) / density_norm
      return
   end subroutine initial_values

   subroutine init_simple_sol
   ! subroutine to initialize the plasma data according to the simple sol model
   !  EXTEND THIS SUBROUTINE with momentum and molecules
      implicit none
      real(wp) :: density_X, temperature_X, temperature_target, temperature_target_new, kappa_0, diff
      ! Arrays are allocated in div1d_step.f90 
	

      velocity    = initial_v
      neutral     = initial_a !dyn_nb(1) !initial_a
      extern_neutral_density = initial_nb !dyn_nb(1) !initial_nb
      extern_molecule_density = initial_mb ! dyn_mb(1) !initial_mb
      core_density = initial_ncore
	!     define kappa_0
      kappa_0     = 2.0d+3
	!     set the X-point density as given on input
      density_X = dyn_nu(1) !initial_n
	!     solve 2PM by iteration starting from temperature_target = 0
      diff = 1.0d+0
      temperature_target = 0.0d+0
      write(*,*) dyn_qparX(1), L, kappa_0
      do while ( diff .gt. 0.01d+0 )
         temperature_X = (temperature_target**(7.d+0/2.d+0) + 7.0d+0 * dyn_qparX(1) * L / 2.0d+0 / kappa_0)**(2.d+0/7.d+0)
         temperature_target_new = (mass / e_charge) * 2.0d+0 * dyn_qparX(1)**2 / temperature_X**2 / (gamma * e_charge * density_X)**2
         diff = abs(temperature_target_new - temperature_target)
         ! write(*,*) temperature_target, temperature_target_new, temperature_X
         temperature_target = 0.1d+0*temperature_target_new+0.9d+0*temperature_target
      end do
      write(*,*) 'Initialization from simple 2 Point Model'
      write(*,*) 'inputs:'
      write(*,*) '        upstream density n_X =', density_X, 'm^-3, upstream heat flux q_parallel,X =', dyn_qparX(1), 'W/m^2, length of divertor leg L =', L, 'm'
      write(*,*) 'results:'
      write(*,*) '        Xpoint temperature 2PM T_X =', temperature_X, 'eV'
      write(*,*) '        target temperature 2PM T_L =', temperature_target, 'eV'

!     now set the temperature solution
      temperature = (temperature_target**(7.d+0/2.d+0) + 7.0d+0 * dyn_qparX(1) * (L-x) / 2.0d+0 / kappa_0)**(2.d+0/7.d+0)
!     set the density solution according to constant pressure (acceleration to sound speed must occur in sheath)
      density = dyn_nu(1) * temperature_X / temperature
      pressure    = 2.0d+0 * density * temperature
      extern2sol_flux = 0.0d-20
      core2sol_flux = 0.0d-20
      extern2sol_mol = 0.0d-20
      core2sol_mol = 0.0d-20
      ! transform (density, velocity, temperature) to (density, momentum, pressure) in solution vector y
      call nvt2ys( Nx, density, velocity, temperature, neutral, neutral_velocity, molecule,  ys )
      call nr2yr(extern_neutral_density, extern_molecule_density, yr)
      return
   end subroutine init_simple_sol

   
   subroutine read_restart_file( restart_error )
   ! subroutine to read a restart file in order to coninue a run
      implicit none
      integer, intent( out ) :: restart_error
      integer     :: Nx_restart
      real( wp ), allocatable  :: x_restart(:)
      real( wp ), allocatable  :: y_restart(:)
      real( wp )  :: L_restart
      real( wp )  :: mass_restart
      real( wp ) :: y(6*Nx+10), yr(10), ys(6*Nx)

        
      restart_error = 0
      write(6,*) 'entered restart_file function'
      ! open the restart file
      open( UNIT = 11, FILE = 'div1d_restart_old.txt', IOSTAT = restart_error )
      ! first read numerics parameters and check for consitency
      read(11,*, IOSTAT = restart_error) Nx_restart
      if( Nx .ne. Nx_restart ) then
         write(*,*) 'grid size in restart file: Nx_restart =', Nx_restart, '     not equal to Nx =', Nx
         write(*,*) 'interpolation is performed between old and new grid'
      endif
      if( allocated(x_restart) .eqv. .false.) allocate( x_restart(Nx_restart), y_restart(6*Nx_restart+10) )
      read(11,*) x_restart
      ! next read physics parameters that must be identical between runs
      read(11,*, IOSTAT = restart_error) L_restart, mass_restart
      if( L .ne. L_restart .or. mass .ne. mass_restart ) then
         write(*,*) 'inconsistent parameters in restart file: '
         write(*,*) 'L_restart    =', L_restart,    '     while L    =', L
         write(*,*) 'mass_restart =', mass_restart, '     while mass =', mass
         close(11)
         return
      endif
      ! next read the plasma data
      read(11,*, IOSTAT = restart_error) y_restart
      close(11)
      ! interpolate between restart grid and current grid
      call interpolate(x_restart, y_restart(             1:  Nx_restart), Nx_restart, x, y(     1:  Nx), Nx)
      call interpolate(x_restart, y_restart(  Nx_restart+1:2*Nx_restart), Nx_restart, x, y(  Nx+1:2*Nx), Nx)
      call interpolate(x_restart, y_restart(2*Nx_restart+1:3*Nx_restart), Nx_restart, x, y(2*Nx+1:3*Nx), Nx)
      call interpolate(x_restart, y_restart(3*Nx_restart+1:4*Nx_restart), Nx_restart, x, y(3*Nx+1:4*Nx), Nx)
      call interpolate(x_restart, y_restart(4*Nx_restart+1:5*Nx_restart), Nx_restart, x, y(4*Nx+1:5*Nx), Nx)
      call interpolate(x_restart, y_restart(5*Nx_restart+1:6*Nx_restart), Nx_restart, x, y(5*Nx+1:6*Nx), Nx)
      ! take over external densities
      y(6*Nx+1:6*Nx+10) = y_restart(6*Nx_restart+1:6*Nx_restart+10) 
      ys = y(1:6*Nx)
      yr = y(6*Nx+1:6*Nx+10)
      ! since the y vector is normalized all quantities are automatically rescaled in accordance with the value of initial_n / density_norm in the current input file
      ! set the secondary plasma variables
      call ys2nvt( Nx, ys, density, velocity, temperature, neutral, neutral_velocity, molecule )
      call yr2nr(yr, extern_neutral_density, extern_molecule_density )
      return
   end subroutine read_restart_file

   subroutine write_restart_file
   ! subroutine to write a restart file in order to be able to continue this run
      implicit none
      ! open the restart file
      open( UNIT = 11, FILE = 'div1d_restart_new.txt' )
      ! first write numerics parameter that cannot change between runs
      write(11,*) Nx
      write(11,*) x
      ! next write physics parameters that must be identical between runs
      write(11,*) L, mass
      ! next write the plasma data
      write(11,*) ys
      write(11,*) yr
      close(11)
      return
   end subroutine write_restart_file
end module plasma_data
