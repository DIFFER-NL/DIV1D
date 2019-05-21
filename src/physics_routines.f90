module physics_routines
! module containing general purpose routines implementing the equations

   use grid_data, only : delta_x
   use constants, only : e_charge
   use reaction_rates
   use physics_parameters, only : gamma, mass, Gamma_X, q_parX, energy_loss_ion, recycling, redistributed_fraction, L, neutral_residence_time, sintheta, minimum_density, minimum_temperature
   use numerics_parameters, only : evolve_density, evolve_momentum, evolve_energy, evolve_neutral, switch_density_source, switch_momentum_source, switch_energy_source, switch_neutral_source, viscosity

   implicit none
   integer, parameter, private :: wp = KIND(1.0D0)

contains

   subroutine nvt2y( Nx, density, velocity, temperature, neutral, y )
   ! subroutine to transform from (density, velocity, temperature) to (density, momentum, energy) in the solution vector y
   !    y(      1, ...,   Nx )     = density( 1, ..., Nx )
   !    y(   Nx+1, ..., 2*Nx )     = Momentum = mass * density * velocity ( 1, ..., Nx )
   !    y( 2*Nx+1, ..., 3*Nx )     = Energy = 3 * density * e_charge * temperature [eV] ( 1, ..., Nx ) + 0.5 mass * density * velocity^2 ( 1, ..., Nx )
   !    y( 3*Nx+1, ..., 4*Nx )     = neutral density( 1, ..., Nx )
      implicit none
      integer,  intent(in)  :: Nx
      real(wp), intent(in)  :: density(Nx), velocity(Nx), temperature(Nx), neutral(Nx)
      real(wp), intent(out) :: y(4*Nx)
      ! density
      y(     1:   Nx) = density(1:Nx)
      ! momentum
      y(  Nx+1: 2*Nx) = mass * density(1:Nx) * velocity(1:Nx)
      ! energy
      y(2*Nx+1: 3*Nx) = 3.0d+0 * density(1:Nx) * e_charge * temperature(1:Nx) + 0.5d+0 * mass * density(1:Nx) * velocity(1:Nx) * velocity(1:Nx)
      ! neutral density
      y(3*Nx+1: 4*Nx) = neutral(1:Nx)
      return
   end subroutine nvt2y


   subroutine y2nvt( Nx, y, density, velocity, temperature, neutral )
   ! subroutine to transform from the solution vector y containing (density, momentum, energy) to (density, velocity, temperature)
   !    y(      1, ...,   Nx )     = density( 1, ..., Nx )
   !    y(   Nx+1, ..., 2*Nx )     = Momentum = mass * density * velocity ( 1, ..., Nx )
   !    y( 2*Nx+1, ..., 3*Nx )     = Energy = 3 * density * e_charge * temperature [eV] ( 1, ..., Nx ) + 0.5 mass * density * velocity^2 ( 1, ..., Nx )
   !    y( 3*Nx+1, ..., 4*Nx )     = neutral density( 1, ..., Nx )
      implicit none
      integer,  intent(in)  :: Nx
      real(wp), intent(in)  :: y(4*Nx)
      real(wp), intent(out) :: density(Nx), velocity(Nx), temperature(Nx), neutral(Nx)
      ! density
      density(1:Nx)     =  max(y(1:Nx), minimum_density)
      ! velocity
      velocity(1:Nx)    =  y(  Nx+1: 2*Nx) / mass / density(1:Nx)
      ! temperature
      temperature(1:Nx) =  max((y(2*Nx+1: 3*Nx) - 0.5d+0 * mass * density(1:Nx) * velocity(1:Nx) * velocity(1:Nx)) / 3.0d+0 / density(1:Nx) / e_charge, minimum_temperature)
      ! neutral density
      neutral(1:Nx)     =  max(y(3*Nx+1: 4*Nx), minimum_density)
      return
   end subroutine y2nvt


   subroutine calculate_fluxes( Nx, density, velocity, temperature, neutral, Gamma_n, pressure, q_parallel, neutral_flux )
   ! this subroutine calculates the 'fluxes' required for the calculation of ydot (the right hand side of the discretized conservation equations)
   ! the fluxes are defined at the halfway point between grid points: i.e. Flux(i) is midway between x(i) and x(i+1)
      implicit none
      integer,  intent(in)  :: Nx
      real(wp), intent(in)  :: density(Nx), velocity(Nx), temperature(Nx), neutral(Nx)
      real(wp), intent(out) :: Gamma_n(Nx), pressure(Nx), q_parallel(Nx), neutral_flux(Nx)
      real(wp)              :: csound
      integer               :: i
      ! sheath sound velocity
      csound = sqrt( 2.0d+0 * e_charge * temperature(Nx) / mass )
      ! the particle flux = density velocity
         Gamma_n    = density * velocity !* 0.0d+0 !!!!!!!!!!!!!!!!!!!!!!!!!!
         ! Gamma_n(1:Nx-1) = 0.5d+0 * (Gamma_n(1:Nx-1) + Gamma_n(2:Nx))
         ! boundary condition at the sheath
         ! Gamma_n(Nx) = density(Nx) * csound * 0.1
      ! total pressure = mass density velocity^2 + 2 density k temperature
         pressure   = mass * density * velocity * velocity + 2.0d+0 * density * e_charge * temperature !* 0.0d+0 !!!!!!!!!!!!!!!!!!!!!!!
         ! pressure(1:Nx-1) = 0.5d+0 * (pressure(1:Nx-1) + pressure(2:Nx))
         ! boundary condition at the sheath
         ! pressure(Nx) = mass * density(Nx) * csound + 2.0d+0 * density(Nx) * e_charge * temperature(Nx)
      ! convective heat flux = 5 density k temperature + 0.5 mass density velocity^3
         q_parallel = 5.0d+0 * density * e_charge * temperature * velocity + 0.5d+0 * mass * density * velocity * velocity * velocity
      ! add the conductive heat flux in the internal region
         do i = 1, Nx-1
            q_parallel(i) = q_parallel(i) * 1.0d+0 - kappa_parallel(0.5d+0*(temperature(i)+temperature(i+1))) * (temperature(i+1)-temperature(i))/delta_x(i) !!!!!!!!!!!!!!!!!!!!!!!!!
         enddo
         ! boundary condition at the sheath: given by the sheath heat transmission
         q_parallel(Nx) = gamma * csound * (density(Nx)/1.0d+0) * e_charge * Temperature(Nx) ! we have equated the density in the sheath to 0.5 * density (Nx) because of the pressure balance, i.e. density_target = 0.5 * density(Nx)
      ! the neutral particle diffusion
         ! we do this in the right_hand_side routine itself
         neutral_flux = 0.0 ! to avoid warning messages at compile time
         ! write(*,*) 'temperature =', temperature
         ! write(*,*) 'q_parallel =', q_parallel
      return
   end subroutine calculate_fluxes


   subroutine calculate_sources( Nx, density, velocity, temperature, neutral, Source_n, Source_v, Source_Q, neutral_source )
   ! this subroutine calculates the source terms of the discretized conservation equations
      implicit none
      integer,  intent(in)  :: Nx
      real(wp), intent(in)  :: density(Nx), velocity(Nx), temperature(Nx), neutral(Nx)
      real(wp), intent(out) :: Source_n(Nx), Source_v(Nx), Source_Q(Nx), neutral_source(Nx)
      real(wp) :: rate_cx(Nx), rate_ion(Nx), rate_rec(Nx)
      integer  :: ix
      Source_n = 0.0d+0
      Source_v = 0.0d+0
      Source_Q = 0.0d+0
      neutral_source = 0.0d+0
      do ix = 1, Nx
         rate_cx(ix)  = density(ix) * neutral(ix) * charge_exchange(temperature(ix))
         rate_ion(ix) = density(ix) * neutral(ix) * ionization(temperature(ix))
         rate_rec(ix) = density(ix) * density(ix) * recombination(density(ix),temperature(ix))
      enddo
      ! the particle sources
      neutral_source = rate_rec - rate_ion
      Source_n = rate_ion - rate_rec
      ! the momentum sources
      Source_v = - mass * velocity * ( rate_cx + rate_rec )
      ! the energy sources
      Source_Q = - (0.5d+0 * mass * velocity * velocity + 1.5d+0 * e_charge * temperature) * (rate_cx + rate_rec)
      Source_Q = source_Q - rate_ion * e_charge * energy_loss_ion
      ! write(*,*) rate_ion, source_Q
      return
   end subroutine calculate_sources


   subroutine right_hand_side( neq, time, y, ydot )
   ! this subroutine calculates the right hand side ydot of the discretized conservation equations
      implicit none
      integer,  intent(in)  :: neq
      real(wp), intent(in)  :: time, y(neq) !time is not used
      real(wp), intent(out) :: ydot(neq)
      integer               :: Nx, ix
      real(wp)              :: density(neq/4), velocity(neq/4), temperature(neq/4), neutral(neq/4)
      real(wp)              :: Gamma_n(neq/4), pressure(neq/4), q_parallel(neq/4), neutral_flux(neq/4)
      real(wp)              :: Source_n(neq/4), Source_v(neq/4), Source_Q(neq/4), neutral_source(neq/4)
      real(wp)              :: Diff_neutral(neq/4)
      real(wp)              :: csound, q_sheath
      Nx = neq/4
      ! write(*,*) 'RHS called at t =', time
      ! write(*,*) 'y =', y
      ! first tranform the solution vector to (density velocity temperature neutral-density)
      call y2nvt( Nx, y, density, velocity, temperature, neutral )
      ! write(*,*) 'density = ', density
      ! write(*,*) 'velocity =', velocity
      ! write(*,*) 'temperature =', temperature
      ! write(*,*) 'neutral density =', neutral
      ! calculate the fluxes
      call calculate_fluxes( Nx, density, velocity, temperature, neutral, Gamma_n, pressure, q_parallel, neutral_flux )
      ! calculate the sources
      call calculate_sources( Nx, density, velocity, temperature, neutral, Source_n, Source_v, Source_Q, neutral_source )
      ! write(*,*) 'Gamma_n =', Gamma_n
      ! write(*,*) 'pressure =', pressure
      ! write(*,*) 'q_parallel =', q_parallel
      ! write(*,*) 'Source_n =', Source_n
      ! write(*,*) 'Source_v =', Source_v
      ! write(*,*) 'Source_Q =', Source_Q
      ! write(*,*) 'neutral_source =', neutral_source
      ! ydot for the density equation
         ydot(1:Nx) = switch_density_source * Source_n(1:Nx)
         ! add the particle flux term - (d nv / dx) in the internal region using upwind differencing
         do ix = 2, Nx-1
            ! if( velocity(ix) .ge. 0.0d+0 ) then
               ydot(ix) = ydot(ix) - density(ix) *(velocity(ix)-velocity(ix-1))/delta_x(ix)
               ydot(ix) = ydot(ix) - velocity(ix)*(density(ix) -density(ix-1)) /delta_x(ix)
            ! elseif( velocity(ix) .lt. 0.0d+0 ) then
            !    ydot(ix) = ydot(ix) - density(ix) *(velocity(ix+1)-velocity(ix))/delta_x(ix)
            !    ydot(ix) = ydot(ix) - velocity(ix)*(density(ix+1) -density(ix)) /delta_x(ix)
            ! else
            !    ydot(ix) = ydot(ix) - density(ix) *(velocity(ix+1)-velocity(ix-1))/delta_x(ix)/2.0d+0
            !    ydot(ix) = ydot(ix) - velocity(ix)*(density(ix+1) -density(ix-1)) /delta_x(ix)/2.0d+0
            ! endif
         enddo
         ! ! apply boundary condition at the X-point, i=1: particle flux given by Gamma(0) = Gamma_X
         ! ydot(1) = ydot(1) - (Gamma_n(1)-Gamma_X)/delta_x(1)
         ! apply boundary condition at the X-point, i=1: fixed density
         ydot(1) = 0.0
         ! apply boundary condition at the sheath entrance, i=Nx: Flux on target = csound * density(Nx)
         csound = sqrt( 2.0d+0 * e_charge * temperature(Nx) / mass )
         ydot(Nx) = ydot(Nx) - (csound*density(Nx)-Gamma_n(Nx))/delta_x(Nx)
         ! ydot(Nx) = ydot(Nx) - (csound*density(Nx)-0.5d+0*(Gamma_n(Nx)+Gamma_n(Nx-1)))/delta_x(Nx)
      ! write(*,*) 'ydot(density) =', ydot(0*Nx+1:1*Nx)
      ! ydot for the momentum equation
         ydot(Nx+1:2*Nx) = switch_momentum_source * Source_v(1:Nx)
         ! add the pressure term in the internal region using downwind differencing
         ydot(Nx+1) = ydot(Nx+1) - (pressure(2)-pressure(1))/delta_x(1)
         do ix = 2, Nx-1
            ! if( velocity(ix) .ge. 0.0d+0 ) then
               ydot(Nx+ix) = ydot(Nx+ix) - (pressure(ix+1)-pressure(ix))/delta_x(ix)
            ! elseif( velocity(ix) .lt. 0.0d+0 ) then
            !    ydot(Nx+ix) = ydot(Nx+ix) - (pressure(ix)-pressure(ix-1))/delta_x(ix)
            ! else
            !   ydot(Nx+ix) = ydot(Nx+ix) - (pressure(ix+1)-pressure(ix-1))/delta_x(ix)/2.0d+0
            ! endif
            ! add effect of numerical viscosity
            ydot(Nx+ix) = ydot(Nx+ix) + viscosity*(velocity(ix+1) + velocity(ix-1)-2.0d0*velocity(ix))
         enddo
         ! ! apply boundary condition at the X-point, i=1: fixed density and particle flux so also fixed momentum
         ! ydot(Nx+1) = 0.0
         ! apply boundary condition at the sheath entrance, i=Nx: for now, assume zero temperture/density gradient across the sheath: total pressure(Nx+1) = 4n(Nx)T(Nx)
         ydot(2*Nx) = ydot(2*Nx) - (4.0d+0*density(Nx)*e_charge*temperature(Nx)-pressure(Nx-1))/delta_x(Nx)/2.0d+0
         ! ydot(2*Nx) = ydot(2*Nx) - (pressure(Nx)-pressure(Nx-1))/delta_x(Nx)
         ydot(2*Nx) = ydot(2*Nx) + viscosity*(csound + velocity(Nx-1)-2.0d0*velocity(Nx))
      ! write(*,*) 'ydot(momentum) =', ydot(1*Nx+1:2*Nx)
      ! ydot for the energy equation
         ydot(2*Nx+1:3*Nx) = switch_energy_source * Source_Q(1:Nx)
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
      ! write(*,*) 'ydot(energy) =', ydot(2*Nx+1:3*Nx)
      ! ydot for the neutral density equation
         ydot(3*Nx+1:4*Nx) = switch_neutral_source * neutral_source(1:Nx)
         ! add the density diffusion in the internal reagion
         ! the neutral particle diffusion coefficient D == n_n kT / m charge_exchange_rate sin^2theta
         do ix = 1, Nx
            Diff_neutral(ix) = e_charge * temperature(ix) / mass / (density(ix) * charge_exchange(temperature(ix))) / sintheta**2
         enddo
         ! write(*,*) 'Diff_neutral =', Diff_neutral
         ydot(3*Nx+2:4*Nx-1) = ydot(3*Nx+2:4*Nx-1) + Diff_neutral(2:Nx-1) * (neutral(3:Nx)-2.0d0*neutral(2:Nx-1)+neutral(1:Nx-2))/delta_x(2:Nx-1)**2
         ydot(3*Nx+2:4*Nx-1) = ydot(3*Nx+2:4*Nx-1) + (Diff_neutral(3:Nx)-Diff_neutral(1:Nx-2))*(neutral(3:Nx)-neutral(1:Nx-2))/4.0d0/delta_x(2:Nx-1)**2
         ! boundary condition at X-point (zero gradient i.e. at i=0 every equals i=1)
         ydot(3*Nx+1) = ydot(3*Nx+1) + Diff_neutral(1) * (neutral(2)-neutral(1))/delta_x(1)**2
         ydot(3*Nx+1) = ydot(3*Nx+1) + (Diff_neutral(2)-Diff_neutral(1))*(neutral(2)-neutral(1))/4.0d0/delta_x(1)**2
         ! boundary condition at sheath: neutral flux = - density * csound * recycling * (1.0d-0 - redistributed_fraction)
         ydot(4*Nx) = ydot(4*Nx) + (density(Nx) * csound * recycling * (1.0d-0 - redistributed_fraction)- 0.5d+0*(Diff_neutral(Nx)+Diff_neutral(Nx-1))*(neutral(Nx)-neutral(Nx-1)))/delta_x(Nx)!!!! ++++ ?????
         ! finally add the neutral sources and losses from redistribution and finite residence time
         ydot(3*Nx+1:4*Nx) = ydot(3*Nx+1:4*Nx) + density(Nx) * csound * recycling * redistributed_fraction / L - neutral / neutral_residence_time
      ! write(*,*) 'ydot(neutrals) =', ydot(3*Nx+1:4*Nx)
      ! apply evolution switches
      ydot(     1:  Nx) = evolve_density  * ydot(     1:  Nx)
      ydot(  Nx+1:2*Nx) = evolve_momentum * ydot(  Nx+1:2*Nx)
      ydot(2*Nx+1:3*Nx) = evolve_energy   * ydot(2*Nx+1:3*Nx)
      ydot(3*Nx+1:4*Nx) = evolve_neutral  * ydot(3*Nx+1:4*Nx)
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
