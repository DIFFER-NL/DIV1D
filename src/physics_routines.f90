module physics_routines
! module containing general purpose routines implementing the equations

   use grid_data, only : delta_x
   use constants, only : e_charge
   use physics_parameters, only : gamma, mass, Gamma_X, q_parX

   implicit none
   integer, parameter, private :: wp = KIND(1.0D0)

contains

   subroutine nvt2y( Nx, density, velocity, temperature, y )
   ! subroutine to transform from (density, velocity, temperature) to (density, momentum, energy) in the solution vector y
   !    y(      1, ...,   Nx )     = density( 1, ..., Nx )
   !    y(   Nx+1, ..., 2*Nx )     = Momentum = mass * density * velocity ( 1, ..., Nx )
   !    y( 2*Nx+1, ..., 3*Nx )     = Energy = 3 * density * e_charge * temperature [eV] ( 1, ..., Nx ) + 0.5 mass * density * velocity^2 ( 1, ..., Nx )
   !!!!!    y( 3*Nx+1, ..., 4*Nx )     = neutral density( 1, ..., Nx ) (not yet implemented) !!!!
      implicit none
      integer,  intent(in)  :: Nx
      real(wp), intent(in)  :: density(Nx), velocity(Nx), temperature(Nx)
      real(wp), intent(out) :: y(3*Nx)
      ! density
      y(     1:   Nx) = density(1:Nx)
      ! momentum
      y(  Nx+1: 2*Nx) = mass * density(1:Nx) * velocity(1:Nx)
      ! energy
      y(2*Nx+1: 3*Nx) = 3.0d+0 * density(1:Nx) * e_charge * temperature(1:Nx) + 0.5d+0 * mass * density(1:Nx) * velocity(1:Nx) * velocity(1:Nx)
      return
   end subroutine nvt2y


   subroutine y2nvt( Nx, y, density, velocity, temperature )
   ! subroutine to transform from the solution vector y containing (density, momentum, energy) to (density, velocity, temperature)
   !    y(      1, ...,   Nx )     = density( 1, ..., Nx )
   !    y(   Nx+1, ..., 2*Nx )     = Momentum = mass * density * velocity ( 1, ..., Nx )
   !    y( 2*Nx+1, ..., 3*Nx )     = Energy = 3 * density * e_charge * temperature [eV] ( 1, ..., Nx ) + 0.5 mass * density * velocity^2 ( 1, ..., Nx )
   !!!!!    y( 3*Nx+1, ..., 4*Nx )     = neutral density( 1, ..., Nx ) (not yet implemented) !!!!
      implicit none
      integer,  intent(in)  :: Nx
      real(wp), intent(in)  :: y(3*Nx)
      real(wp), intent(out) :: density(Nx), velocity(Nx), temperature(Nx)
      ! density
      density(1:Nx)     =  y(1:Nx)
      ! velocity
      velocity(1:Nx)    =  y(  Nx+1: 2*Nx) / mass / density(1:Nx)
      ! temperature
      temperature(1:Nx) = (y(2*Nx+1: 3*Nx) - 0.5d+0 * mass * density(1:Nx) * velocity(1:Nx) * velocity(1:Nx)) / 3.0d+0 / density(1:Nx) / e_charge
      return
   end subroutine y2nvt


   subroutine calculate_fluxes( Nx, density, velocity, temperature, Gamma_n, pressure, q_parallel )
   ! this subroutine calculates the 'fluxes' required for the calculation of ydot (the right hand side of the discretized conservation equations)
   ! the fluxes are defined at the halfway point between grid points: i.e. Flux(i) is midway between x(i) and x(i+1)
      implicit none
      integer,  intent(in)  :: Nx
      real(wp), intent(in)  :: density(Nx), velocity(Nx), temperature(Nx)
      real(wp), intent(out) :: Gamma_n(Nx), pressure(Nx), q_parallel(Nx)
      real(wp)              :: csound
      integer               :: i
      ! sheath sound velocity
      csound = sqrt( 2.0d+0 * e_charge * temperature(Nx) / mass )
      ! the particle flux = density velocity
         Gamma_n    = density * velocity !* 0.0d+0 !!!!!!!!!!!!!!!!!!!!!!!!!!
         Gamma_n(1:Nx-1) = 0.5d+0 * (Gamma_n(1:Nx-1) + Gamma_n(2:Nx))
         ! boundary condition at the sheath
         Gamma_n(Nx) = density(Nx) * csound * 0.1
      ! total pressure = mass density velocity^2 + 2 density k temperature
         pressure   = mass * density * velocity + 2.0d+0 * density * e_charge * temperature !* 0.0d+0 !!!!!!!!!!!!!!!!!!!!!!!
         pressure(1:Nx-1) = 0.5d+0 * (pressure(1:Nx-1) + pressure(2:Nx))
         ! boundary condition at the sheath
         pressure(Nx) = mass * density(Nx) * csound + 2.0d+0 * density(Nx) * e_charge * temperature(Nx)
      ! convective heat flux = 5 density k temperature + 0.5 mass density velocity^3
         q_parallel = 5.0d+0 * density * e_charge * temperature * velocity + 0.5d+0 * mass * density * velocity * velocity * velocity
      ! add the conductive heat flux in the internal region
         do i = 1, Nx-1
            q_parallel(i) = q_parallel(i) * 1.0d+0 - kappa_parallel(0.5d+0*(temperature(i)+temperature(i+1))) * (temperature(i+1)-temperature(i))/delta_x(i) !!!!!!!!!!!!!!!!!!!!!!!!!
         enddo
         ! boundary condition at the sheath: given by the sheath heat transmission
         q_parallel(Nx) = gamma * csound * (density(Nx)/1.0d+0) * e_charge * Temperature(Nx) ! we have equated the density in the sheath to 0.5 * density (Nx) because of the pressure balance, i.e. density_target = 0.5 * density(Nx)
      return
   end subroutine calculate_fluxes


   subroutine calculate_sources( Nx, density, velocity, temperature, Source_n, Source_v, Source_Q )
   ! this subroutine calculates the source terms of the discretized conservation equations
      implicit none
      integer,  intent(in)  :: Nx
      real(wp), intent(in)  :: density(Nx), velocity(Nx), temperature(Nx)
      real(wp), intent(out) :: Source_n(Nx), Source_v(Nx), Source_Q(Nx)
      Source_n = 0.0d+0
      Source_v = 0.0d+0
      Source_Q = 0.0d+0
   end subroutine calculate_sources


   subroutine right_hand_side( neq, time, y, ydot )
   ! this subroutine calculates the right hand side ydot of the discretized conservation equations
      implicit none
      integer,  intent(in)  :: neq
      real(wp), intent(in)  :: time, y(neq) !time is not used
      real(wp), intent(out) :: ydot(neq)
      integer               :: Nx
      real(wp)              :: density(neq/3), velocity(neq/3), temperature(neq/3)
      real(wp)              :: Gamma_n(neq/3), pressure(neq/3), q_parallel(neq/3)
      real(wp)              :: Source_n(neq/3), Source_v(neq/3), Source_Q(neq/3)
      ! real(wp)              :: csound, q_sheath
      Nx = neq/3
      write(*,*) 'RHS called at t =', time
      ! write(*,*) 'y =', y
      ! first tranform the solution vector to (density velcoity temperature)
      call y2nvt( Nx, y, density, velocity, temperature )
      ! write(*,*) 'density = ', density
      ! write(*,*) 'velocity =', velocity
      ! write(*,*) 'temperature =', temperature
      ! calculate the fluxes
      call calculate_fluxes( Nx, density, velocity, temperature, Gamma_n, pressure, q_parallel )
      ! calculate the sources
      call calculate_sources( Nx, density, velocity, temperature, Source_n, Source_v, Source_Q )
      ! write(*,*) 'Gamma_n =', Gamma_n
      ! write(*,*) 'pressure =', pressure
      ! write(*,*) 'q_parallel =', q_parallel
      ! ydot for the density equation
         ydot(1:Nx) = Source_n(1:Nx)
         ! add the particle flux term in the internal region
         ydot(2:Nx) = ydot(2:Nx) - 2.0d+0*(Gamma_n(2:Nx)-Gamma_n(1:Nx-1))/(delta_x(1:Nx-1)+delta_x(2:Nx))
         ! apply boundary condition at the X-point, i=1: particle flux given by Gamma(0) = Gamma_X
         ydot(1) = ydot(1) - (Gamma_n(1)-Gamma_X)/delta_x(1)
         ! ! apply boundary condition at the sheath entrance, i=Nx: for now, assume full ionization of reflected neutrals in sheath giving a net particle flux Gamma(Nx+1) = 0
         ! ydot(Nx) = ydot(Nx) + Gamma_n(Nx-1)/delta_x(Nx)
      ! ydot for the momentum equation
         ydot(Nx+1:2*Nx) = Source_v(1:Nx)
         ! add the pressure term in the internal region
         ydot(Nx+2:2*Nx) = ydot(Nx+2:2*Nx) - 2.0d+0*(pressure(2:Nx)-pressure(1:Nx-1))/(delta_x(1:Nx-1)+delta_x(2:Nx))
         ! apply boundary condition at the X-point, i=1: assume zero pressure gradient pressure(0) = pressure(1)
         ! ydot(Nx+1) = ydot(Nx+1) - (pressure(1)-pressure(0))/delta_x(1)
         ! apply boundary condition at the sheath entrance, i=Nx: for now, assume zero pressure gradient across the sheath pressure(Nx+1) = pressure(Nx)
         ydot(2*Nx) = ydot(2*Nx) - (pressure(Nx)-pressure(Nx-1))/delta_x(Nx)
      ! reset density and momentum changes to zero for now
      !   ydot(1:Nx) = 0.0d+0
      ! ydot for the energy equation
         ydot(2*Nx+1:3*Nx) = Source_Q(1:Nx)
         ! add the heat flux term in the internal region (including the sheath)
         ydot(2*Nx+2:3*Nx) = ydot(2*Nx+2:3*Nx) - 2.0d+0*(q_parallel(2:Nx)-q_parallel(1:Nx-1))/(delta_x(1:Nx-1)+delta_x(2:Nx))
         ! apply boundary condition at the X-point, i=1: energy flux given by q_parallel(0) = q_parX
         ydot(2*Nx+1) = ydot(2*Nx+1) - (q_parallel(1)-q_parX)/delta_x(1)
         ! ! apply boundary condition at the sheath entrance, i=Nx: heat flux is given by sheath heat transmission q_par(Nx+1) = gamma * csound * density * k * temperature
         ! ! sheath sound velocity
         ! csound = sqrt( 2.0d+0 * e_charge * temperature(Nx) / mass )
         ! ! sheath heat transmission
         ! q_sheath = gamma * csound * (density(Nx)/2.0d+0) * e_charge * Temperature(Nx) ! we have equated the density in the sheath to 0.5 * density (Nx) because of the pressure balance, i.e. density_target = 0.5 * density(Nx)
         ! ydot(3*Nx) = ydot(3*Nx) - (q_sheath-q_parallel(Nx-1))/delta_x(Nx)
      ! write(*,*) 'ydot =', ydot
      return
   end subroutine right_hand_side

   real(wp) function kappa_parallel( temperature )
   ! function to calculate the parallel heat conductivity
      implicit none
      real(wp) :: temperature
      ! use expression from Stangeby page 187 (Chapter 4.10.1)
      kappa_parallel = 2.0d+3 * temperature*temperature*sqrt(temperature)
      return
   end function kappa_parallel
   
end module physics_routines
