module plasma_data
! module defining and handling the plasma data

   use numerics_parameters, only : Nx
   use physics_parameters, only  : initial_n, initial_v, initial_T, initial_a, mass, Gamma_X
   use physics_routines

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
   !    y( 3*Nx+1, ..., 4*Nx )     = neutral density( 1, ..., Nx ) (not yet implemented)
      implicit none
      ! first allocate all arrays
      allocate( y(3*Nx), ydot(3*Nx), density(Nx), velocity(Nx), temperature(Nx), neutral(Nx) )
      allocate( Gamma_n(Nx), pressure(Nx), q_parallel(Nx), neutral_flux(Nx), Source_n(Nx), Source_v(Nx), Source_Q(Nx), source_neutral(Nx) )
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

end module plasma_data
