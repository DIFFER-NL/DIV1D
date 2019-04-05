module plasma_data
! module defining and handling the plasma data

   use numerics_parameters, only : Nx
   use physics_parameters, only  : initial_n, initial_v, initial_T, mass, Gamma_X
   use physics_routines

   implicit none
   integer, parameter, private :: wp = KIND(1.0D0)
   integer ::                 neq             ! number of equation (= 3 Nx (= 4 Nx when adding the neutrals))
   real( wp ), allocatable :: y(:)            ! vector holding the solution for the plasma variables [various units]
   real( wp ), allocatable :: ydot(:)         ! vector holding the solution for the plasma variables [various units]
   real( wp ), allocatable :: density(:)      ! vector holding the solution for the plasma variables [/m^3]
   real( wp ), allocatable :: velocity(:)     ! vector holding the solution for the plasma variables [m/s]
   real( wp ), allocatable :: temperature(:)  ! vector holding the solution for the plasma variables [eV]
   real( wp ), allocatable :: Gamma_n(:)      ! vector holding the solution for the plasma variables [/m^2s]
   real( wp ), allocatable :: pressure(:)     ! vector holding the solution for the plasma variables [Pa]
   real( wp ), allocatable :: q_parallel(:)   ! vector holding the solution for the plasma variables [W/m^2]
   real( wp ), allocatable :: Source_n(:)     ! vector holding the solution for the plasma variables [/m^3s]
   real( wp ), allocatable :: Source_v(:)     ! vector holding the solution for the plasma variables [?]
   real( wp ), allocatable :: Source_Q(:)     ! vector holding the solution for the plasma variables [W/m^3]

contains

   subroutine initial_values
   ! subroutine to initialize the plasma data and solution vetor
   ! the solution vector is defined as follows
   !    y(      1, ...,   Nx )     = density( 1, ..., Nx )
   !    y(   Nx+1, ..., 2*Nx )     = Momentum = mass * density * velocity ( 1, ..., Nx )
   !    y( 2*Nx+1, ..., 3*Nx )     = Energy = 3 * density * e_charge * temperature [eV] ( 1, ..., Nx ) + 0.5 mass * density * velocity^2 ( 1, ..., Nx )
   !!!!!    y( 3*Nx+1, ..., 4*Nx )     = neutral density( 1, ..., Nx ) (not yet implemented) !!!!
      implicit none
      ! first allocate all arrays
      allocate( y(3*Nx), ydot(3*Nx), density(Nx), velocity(Nx), temperature(Nx) )
      allocate( Gamma_n(Nx), pressure(Nx), q_parallel(Nx), Source_n(Nx), Source_v(Nx), Source_Q(Nx) )
      density     = initial_n
      velocity    = initial_v
      if( Gamma_X .ne. 0.0 ) velocity = Gamma_X / initial_n
      temperature = initial_T
      ! transform (density, velocity, temperature) to (density, momentum, energy) in solution vector y
      call nvt2y( Nx, density, velocity, temperature, y )
      return
   end subroutine initial_values

end module plasma_data
