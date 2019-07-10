module plasma_data
! module defining and handling the plasma data

   use numerics_parameters, only : Nx
   use grid_data, only           : x
   use physics_parameters, only  : L, initial_n, initial_v, initial_T, initial_a, mass, Gamma_X
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
   real( wp ), allocatable :: Gamma_n(:)      ! vector holding the solution for the plasma particle flux [/m^2s]
   real( wp ), allocatable :: Gamma_mom(:)    ! vector holding the solution for the plasma momentum flux [Pa]
   real( wp ), allocatable :: pressure(:)     ! vector holding the solution for the plasma pressure [Pa]
   real( wp ), allocatable :: q_parallel(:)   ! vector holding the solution for the plasma heat flux [W/m^2]
   real( wp ), allocatable :: neutral_flux(:) ! vector holding the solution for the plasma particle flux [/m^2s]
   real( wp ), allocatable :: Source_n(:)     ! vector holding the solution for the plasma particle source [/m^3s]
   real( wp ), allocatable :: Source_v(:)     ! vector holding the solution for the plasma momentume source [?]
   real( wp ), allocatable :: Source_Q(:)     ! vector holding the solution for the plasma heat source [W/m^3]
   real( wp ), allocatable :: source_neutral(:)     ! vector holding the solution for the plasma heat source [/m^3s]

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
      allocate( Gamma_n(Nx), Gamma_mom(Nx), pressure(Nx), q_parallel(Nx), neutral_flux(Nx), Source_n(Nx), Source_v(Nx), Source_Q(Nx), source_neutral(Nx) )
      temperature = initial_T
      density     = initial_n
      velocity    = initial_v
      neutral     = initial_a
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
      ! reset the density at the X-point boundary to the value in the current input file
      ! this can be changed from the original run
      y(1) = initial_n
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
