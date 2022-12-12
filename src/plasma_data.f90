module plasma_data
! module defining and handling the plasma data

   use numerics_parameters, only : Nx
   use grid_data, only           : x
   use physics_parameters, only  : L, dyn_nu, initial_v, initial_T, dyn_nb, mass, dyn_qparX, Gamma_X, gamma 
   ! note dyn_nu, dyn_nb, and dyn_qpar replaced initial_n and initial_a and q_parX
   use physics_routines
   use interpolation

   implicit none
   integer, parameter, private :: wp = KIND(1.0D0)
   integer ::                 neq             ! number of equation (= 3 Nx (= 4 Nx when adding the neutrals))
   real( wp ), allocatable :: y(:)            ! vector holding the solution for the plasma variables [various units]
   real( wp ), allocatable :: ydot(:)         ! vector holding the solution for the plasma variables [various units]
   real( wp ), allocatable :: density(:)      ! vector holding the solution for the plasma density [/m^3]
   real( wp ), allocatable :: velocity(:)     ! vector holding the solution for the plasma velocity [m/s]
   real( wp ), allocatable :: temperature(:)  ! vector holding the solution for the plasma temperature [eV]
   real( wp ), allocatable :: neutral(:)      ! vector holding the solution for the neutral density [/m^3]
!   real( wp ), allocatable :: extern_neutral_density(:) ! vector holding the neutral density of external volumes [m^-3]
   real( wp ), allocatable :: Gamma_n(:)      ! vector holding the solution for the plasma particle flux [/m^2s]
   real( wp ), allocatable :: Gamma_mom(:)    ! vector holding the solution for the plasma momentum flux [Pa]
   real( wp ), allocatable :: pressure(:)     ! vector holding the solution for the plasma pressure [Pa]
   real( wp ), allocatable :: q_parallel(:)   ! vector holding the solution for the plasma heat flux [W/m^2]
   real( wp ), allocatable :: neutral_flux(:) ! vector holding the solution for the plasma particle flux [/m^2s]
   real( wp ), allocatable :: extern2sol_flux(:) ! vector holding the fluxes from external neutral volumes into the SOL [1/s]
   real( wp ), allocatable :: Source_n(:)     ! vector holding the solution for the plasma particle source [/m^3s]
   real( wp ), allocatable :: Source_v(:)     ! vector holding the solution for the plasma momentume source [?]
   real( wp ), allocatable :: Source_Q(:)     ! vector holding the solution for the plasma heat source [W/m^3]
   real( wp ), allocatable :: source_neutral(:)     ! vector holding the solution for the plasma heat source [/m^3s]
   real( wp )  :: extern_neutral_flux(3) = (/0.0d+0,0.0d+0,0.0d+0/) ! vector with neutral fluxes between external volumes [1/s]
   real( wp )  :: sol2extern_flux(5) = (/0.0d+0,0.0d+0,0.0d+0,0.0d+0,0.0d+0/) ! vector with neutral fluxes from SOL into external volumes [1/s]
contains

   subroutine initial_values
   ! subroutine to initialize the plasma data and solution vetor
   ! the solution vector is defined as follows
   !    y(      1, ...,   Nx )     = density( 1, ..., Nx )
   !    y(   Nx+1, ..., 2*Nx )     = Momentum = mass * density * velocity ( 1, ..., Nx )
   !    y( 2*Nx+1, ..., 3*Nx )     = pressure = 2 * density * e_charge * temperature [eV] ( 1, ..., Nx )
   !    y( 3*Nx+1, ..., 4*Nx )     = neutral density( 1, ..., Nx )
      implicit none
      ! first allocate all arrays
      allocate( y(4*Nx), ydot(4*Nx), density(Nx), velocity(Nx), temperature(Nx), neutral(Nx) )
      allocate( Gamma_n(0:Nx), Gamma_mom(0:Nx), pressure(Nx), q_parallel(0:Nx), neutral_flux(0:Nx), Source_n(Nx), Source_v(Nx), Source_Q(Nx), source_neutral(Nx) )
      temperature = initial_T
      density     = dyn_nu(1) ! initial_n
      velocity    = initial_v
      neutral     = dyn_nb(1) ! initial_a
      pressure    = 2.0d+0 * density * temperature
!      if( Gamma_X .ne. 0.0 ) then
!         ! initialize velocity at the sound speed and density in accordance with Gamma_X
!         velocity = sqrt( 2.0d+0 * e_charge * temperature / mass )
!         density  = Gamma_X / velocity
!      endif
      ! transform (density, velocity, temperature) to (density, momentum, pressure) in solution vector y
      call nvt2y( Nx, density, velocity, temperature, neutral, y )
      return
   end subroutine initial_values

   subroutine init_simple_sol
   ! subroutine to initialize the plasma data according to the simple sol model
      implicit none
      real(wp) :: density_X, temperature_X, temperature_target, temperature_target_new, kappa_0, diff
      ! first allocate all arrays
      allocate( y(4*Nx), ydot(4*Nx), density(Nx), velocity(Nx), temperature(Nx), neutral(Nx) )
      allocate( Gamma_n(Nx), Gamma_mom(Nx), pressure(Nx), q_parallel(Nx), neutral_flux(Nx), Source_n(Nx), Source_v(Nx), Source_Q(Nx), source_neutral(Nx) )
!     set the velocity and neutral density arrays as defined in input
      velocity    = initial_v
      neutral     = dyn_nb(1) !initial_a
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
      ! transform (density, velocity, temperature) to (density, momentum, pressure) in solution vector y
      call nvt2y( Nx, density, velocity, temperature, neutral, y )
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
      restart_error = 0
      ! open the restart file
      open( UNIT = 10, FILE = 'div1d_restart_old.txt', IOSTAT = restart_error )
      ! first read numerics parameters and check for consitency
      read(10,*, IOSTAT = restart_error) Nx_restart
      if( Nx .ne. Nx_restart ) then
         write(*,*) 'grid size in restart file: Nx_restart =', Nx_restart, '     not equal to Nx =', Nx
         write(*,*) 'interpolation is performed between old and new grid'
      endif
      allocate( x_restart(Nx_restart), y_restart(4*Nx_restart) )
      read(10,*) x_restart
      ! next read physics parameters that must be identical between runs
      read(10,*, IOSTAT = restart_error) L_restart, mass_restart
      if( L .ne. L_restart .or. mass .ne. mass_restart ) then
         write(*,*) 'inconsistent parameters in restart file: '
         write(*,*) 'L_restart    =', L_restart,    '     while L    =', L
         write(*,*) 'mass_restart =', mass_restart, '     while mass =', mass
         close(10)
         return
      endif
      ! next read the plasma data
      read(10,*, IOSTAT = restart_error) y_restart
      close(10)
      ! interpolate between restart grid and current grid
      call interpolate(x_restart, y_restart(             1:  Nx_restart), Nx_restart, x, y(     1:  Nx), Nx)
      call interpolate(x_restart, y_restart(  Nx_restart+1:2*Nx_restart), Nx_restart, x, y(  Nx+1:2*Nx), Nx)
      call interpolate(x_restart, y_restart(2*Nx_restart+1:3*Nx_restart), Nx_restart, x, y(2*Nx+1:3*Nx), Nx)
      call interpolate(x_restart, y_restart(3*Nx_restart+1:4*Nx_restart), Nx_restart, x, y(3*Nx+1:4*Nx), Nx)
      ! since the y vector is normalized all quantities are automatically rescaled in accordance with the value of initial_n / density_norm in the current input file
      ! set the secondary plasma variables
      call y2nvt( Nx, y, density, velocity, temperature, neutral )
      return
   end subroutine read_restart_file

   subroutine write_restart_file
   ! subroutine to write a restart file in order to be able to continue this run
      implicit none
      ! open the restart file
      open( UNIT = 10, FILE = 'div1d_restart_new.txt' )
      ! first write numerics parameter that cannot change between runs
      write(10,*) Nx
      write(10,*) x
      ! next write physics parameters that must be identical between runs
      write(10,*) L, mass
      ! next write the plasma data
      write(10,*) y
      close(10)
      return
   end subroutine write_restart_file


end module plasma_data
