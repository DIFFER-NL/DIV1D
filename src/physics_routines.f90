module physics_routines
! module containing general purpose routines implementing the equations

   use grid_data, only : delta_x, delta_xcb
   use constants, only : e_charge
   use reaction_rates
   use physics_parameters, only : gamma, mass, Gamma_X, q_parX, energy_loss_ion, recycling, redistributed_fraction, L, neutral_residence_time, sintheta, minimum_density, minimum_temperature, density_ramp_rate
   use numerics_parameters, only : evolve_density, evolve_momentum, evolve_energy, evolve_neutral, switch_density_source, switch_momentum_source, switch_energy_source, switch_neutral_source, &
                                   switch_convective_heat, switch_impurity_radiation, viscosity, density_norm, momentum_norm, energy_norm

   implicit none
   integer, parameter, private :: wp = KIND(1.0D0)

contains

   subroutine nvt2y( Nx, density, velocity, temperature, neutral, y )
   ! subroutine to transform from (density, velocity, temperature) to (density, momentum, energy) in the solution vector y
   !    y(      1, ...,   Nx )     = density( 1, ..., Nx )
   !    y(   Nx+1, ..., 2*Nx )     = Momentum = mass * density * velocity ( 1, ..., Nx )
   !    y( 2*Nx+1, ..., 3*Nx )     = internal energy = 3 * density * e_charge * temperature [eV] ( 1, ..., Nx )
   !    y( 3*Nx+1, ..., 4*Nx )     = neutral density( 1, ..., Nx )
      implicit none
      integer,  intent(in)  :: Nx
      real(wp), intent(in)  :: density(Nx), velocity(Nx), temperature(Nx), neutral(Nx)
      real(wp), intent(out) :: y(4*Nx)
      ! density
      y(     1:   Nx) = density(1:Nx) / density_norm
      ! momentum
      y(  Nx+1: 2*Nx) = mass * density(1:Nx) * velocity(1:Nx) / momentum_norm
      ! energy
      y(2*Nx+1: 3*Nx) = 3.0d+0 * density(1:Nx) * e_charge * temperature(1:Nx) / energy_norm
      ! neutral density
      y(3*Nx+1: 4*Nx) = neutral(1:Nx) / density_norm
      return
   end subroutine nvt2y


   subroutine y2nvt( Nx, y, density, velocity, temperature, neutral )
   ! subroutine to transform from the solution vector y containing (density, momentum, energy) to (density, velocity, temperature)
   ! note that the solution vector y is normalized
   !    y(      1, ...,   Nx )     = density( 1, ..., Nx ) / density_norm
   !    y(   Nx+1, ..., 2*Nx )     = Momentum / momentum_norm = mass * density/density_norm * velocity/velocity_norm ( 1, ..., Nx ) 
   !    y( 2*Nx+1, ..., 3*Nx )     = internal energy / energy_norm = 3 * density/density_norm * e_charge * temperature/temperature_norm ( 1, ..., Nx )
   !    y( 3*Nx+1, ..., 4*Nx )     = neutral density( 1, ..., Nx )
      implicit none
      integer,  intent(in)  :: Nx
      real(wp), intent(in)  :: y(4*Nx)
      real(wp), intent(out) :: density(Nx), velocity(Nx), temperature(Nx), neutral(Nx)
      ! density
      density(1:Nx)     =  y(1:Nx) * density_norm
      ! velocity
      velocity(1:Nx)    =  y(  Nx+1: 2*Nx) * momentum_norm / mass / density(1:Nx)
      ! temperature
      temperature(1:Nx) =  y(2*Nx+1: 3*Nx) * energy_norm / 3.0d+0 / density(1:Nx) / e_charge
      ! neutral density
      neutral(1:Nx)     =  y(3*Nx+1: 4*Nx) * density_norm
      return
   end subroutine y2nvt


   subroutine advection(Nx, variable, velocity, temperature, flux)
   ! this subroutine calculates the advected flux of a variable with a given velocity field
   ! the temperature is needed to calculate the sound velocity
   ! the variable and velcoty field are given on the (equidistant) grid centers
   ! the flux is returned at the grid boundaries i + 1/2
      implicit none
      integer,  intent(in)  :: Nx
      real(wp), intent(in)  :: variable(Nx), velocity(Nx), temperature(Nx)
      real(wp), intent(out) :: flux(Nx)
      real(wp)              :: csound, average_velocity
      integer               :: i
      flux = 0.0
      ! we follow here the discretization as put forward in B. Dudson et al. (2019) PPCF 61 065008
      do i = 1, Nx-1
         csound = sqrt( 2.0d+0 * e_charge * max(temperature(i),temperature(i+1)) / mass )
         average_velocity = 0.5d+0 * (velocity(i)+velocity(i+1))
         if( average_velocity .gt. csound ) then
            if( i .eq. 1 ) then
               flux(i) = average_velocity *  variable(i)
            else
               flux(i) = average_velocity * (variable(i) + 0.5 * MC_limit(variable(i-1),variable(i),variable(i+1)))
            endif
         elseif( average_velocity .lt. -csound ) then
            if( i .eq. Nx-1 ) then
               flux(i) = average_velocity * variable(Nx)
            else
               flux(i) = average_velocity * (variable(i+1) - 0.5 * MC_limit(variable(i),variable(i+1),variable(i+2)))
            endif
         else
            if( i .eq. 1 ) then
               flux(i) = 0.5 * ( (average_velocity+csound) *  variable(i) + &
               &                    (average_velocity-csound) * (variable(i+1) - 0.5 * MC_limit(variable(i),variable(i+1),variable(i+2))) )
            elseif( i .eq. Nx-1 ) then
               flux(i) = 0.5 * ( (average_velocity+csound) * (variable(i) + 0.5 * MC_limit(variable(i-1),variable(i),variable(i+1))) + &
               &                    (average_velocity-csound) *  variable(Nx) )
            else
               flux(i) = 0.5 * ( (average_velocity+csound) * (variable(i) + 0.5 * MC_limit(variable(i-1),variable(i),variable(i+1))) + &
               &                    (average_velocity-csound) * (variable(i+1) - 0.5 * MC_limit(variable(i),variable(i+1),variable(i+2))) )
            endif
         endif
      enddo
   end subroutine advection


   subroutine calculate_fluxes( Nx, density, velocity, temperature, neutral, Gamma_n, Gamma_mom, q_parallel, neutral_flux )
   ! this subroutine calculates the 'fluxes' required for the calculation of ydot (the right hand side of the discretized conservation equations)
   ! the fluxes are defined at the halfway point between grid points: i.e. Flux(i) is midway between x(i) and x(i+1)
      implicit none
      integer,  intent(in)  :: Nx
      real(wp), intent(in)  :: density(Nx), velocity(Nx), temperature(Nx), neutral(Nx)
      real(wp), intent(out) :: Gamma_n(Nx), Gamma_mom(Nx), q_parallel(Nx), neutral_flux(Nx)
      real(wp)              :: momentum(Nx), enthalpy(Nx)
      real(wp)              :: csound, average_velocity
      integer               :: i
      ! the particle flux = density velocity
      ! we follow here the discretization as put forward in B. Dudson et al. (2019) PPCF 61 065008
         call advection(Nx, density, velocity, temperature, Gamma_n)
         ! boundary condition at the sheath (note that velocity is allowed to be supersonic)
         csound = sqrt( 2.0d+0 * e_charge * temperature(Nx) / mass )
         Gamma_n(Nx) = (1.5*density(Nx)-0.5*density(Nx-1)) * max(velocity(Nx),csound)
         if(temperature(Nx) .le. minimum_temperature) Gamma_n(Nx) = 0.0d+0
      ! the momentum flux = momentum * velocity where momentum = density * mass * velocity
         momentum = density * mass * velocity
         call advection(Nx, momentum, velocity, temperature, Gamma_mom)
         ! boundary condition at the sheath
         Gamma_mom(Nx) = (1.5*density(Nx)-0.5*density(Nx-1)) * mass * max(velocity(Nx),csound)**2
         if(temperature(Nx) .le. minimum_temperature) Gamma_mom(Nx) = 0.0d+0
      ! convective heat flux = 5 density k temperature velocity (i.e. 5/2 pressure)
         enthalpy = 5.0d+0 * density * e_charge * temperature
         call advection(Nx, enthalpy, velocity, temperature, q_parallel) 
      ! add the conductive heat flux in the internal region
         do i = 1, Nx-1
            q_parallel(i) = q_parallel(i) * switch_convective_heat - kappa_parallel(0.5d+0*(temperature(i)+temperature(i+1))) * (temperature(i+1)-temperature(i))/delta_x(i)
         enddo
         ! boundary condition at the sheath: given by the sheath heat transmission
         q_parallel(Nx) = gamma * csound * (1.5*density(Nx)-0.5*density(Nx-1)) * e_charge * Temperature(Nx) ! we have extrapolated the density linear towards x = L, i.e. the sheath
         if(temperature(Nx) .le. minimum_temperature .or. q_parallel(Nx) .lt. 0.0d+0) q_parallel(Nx) = 0.0d+0
      ! the neutral particle diffusion !!!! switch-on in case you want this diagnostic
         ! we do this in the right_hand_side routine itself
         do i = 1, Nx-1
            neutral_flux(i) = - 0.5d+0*(D_neutral(temperature(i),density(i))+D_neutral(temperature(i+1),density(i+1))) * (neutral(i+1)-neutral(i))/delta_x(i)
         enddo
         ! boundary condition at the sheath (- flux of plasma density in case of full recycling)
         neutral_flux(Nx) = - Gamma_n(Nx) * recycling * (1.0d-0 - redistributed_fraction)
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
      real(wp) :: rate_cx(Nx), rate_ion(Nx), rate_exc(Nx), rate_rec(Nx), rate_imp(Nx)
      integer  :: ix
      Source_n = 0.0d+0
      Source_v = 0.0d+0
      Source_Q = 0.0d+0
      neutral_source = 0.0d+0
      do ix = 1, Nx
         rate_cx(ix)  = density(ix) * neutral(ix) * charge_exchange(temperature(ix))
         rate_ion(ix) = density(ix) * neutral(ix) * ionization(density(ix),temperature(ix))
         rate_exc(ix) = density(ix) * neutral(ix) * excitation(density(ix),temperature(ix))
         rate_rec(ix) = density(ix) * density(ix) * recombination(density(ix),temperature(ix))
         rate_imp(ix) = density(ix) * density(ix) * impurity_radiation(temperature(ix))
      enddo
      ! the particle sources
      neutral_source = rate_rec - rate_ion
      Source_n = rate_ion - rate_rec
      ! the momentum sources
      Source_v = - mass * velocity * ( rate_cx + rate_rec )
      ! the energy sources (only internal energy)
      Source_Q = - (1.5d+0 * e_charge * temperature) * (rate_cx + 2.0d+0*rate_rec) ! the factor 2 for the recombination rate accounts for the loss of both an ion and an electron
      if ( switch_excitation .eq. 0.0d+0 ) then
         Source_Q = Source_Q - rate_ion * e_charge * energy_loss_ion ! note energy loss per ionization is in eV
      else
         Source_Q = Source_Q - switch_excitation * rate_exc * e_charge ! note excitation rate is in eV m^3 / s
      endif
      Source_Q = Source_Q - switch_impurity_radiation * rate_imp * e_charge ! note impurity radiation loss rate is in eV m^3 / s
      ! write(*,*) rate_ion, Source_Q
      return
   end subroutine calculate_sources


   subroutine right_hand_side( neq, time, y, ydot )
   ! this subroutine calculates the right hand side ydot of the discretized conservation equations
   ! note that y and ydor are normalized, but the arrays density, velocity, temperature, and neutral are not!
      implicit none
      integer,  intent(in)  :: neq
      real(wp), intent(in)  :: time, y(neq) !time is not used
      real(wp), intent(out) :: ydot(neq)
      integer               :: Nx, ix
      real(wp)              :: density(neq/4), velocity(neq/4), temperature(neq/4), neutral(neq/4)
      real(wp)              :: Gamma_n(neq/4), Gamma_mom(neq/4), q_parallel(neq/4), neutral_flux(neq/4)
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
      call calculate_fluxes( Nx, density, velocity, temperature, neutral, Gamma_n, Gamma_mom, q_parallel, neutral_flux )
      ! calculate the sources
      call calculate_sources( Nx, density, velocity, temperature, neutral, Source_n, Source_v, Source_Q, neutral_source )
      ! write(*,*) 'Gamma_n =', Gamma_n
      ! write(*,*) 'Gamma_mom =', Gamma_mom
      ! write(*,*) 'q_parallel =', q_parallel
      ! write(*,*) 'Source_n =', Source_n
      ! write(*,*) 'Source_v =', Source_v
      ! write(*,*) 'Source_Q =', Source_Q
      ! write(*,*) 'neutral_source =', neutral_source
      ! sound velocity at the target
      csound = sqrt( 2.0d+0 * e_charge * temperature(Nx) / mass )
      ! ydot for the density equation
         ydot(1:Nx) = switch_density_source * Source_n(1:Nx)
         ! add the particle flux term using the flux as calculated in calculate_fluxes
         ! Gamma_n(i) contains the flux at i+1/2
         do ix = 2, Nx
            ydot(ix) = ydot(ix) - (Gamma_n(ix)-Gamma_n(ix-1))/delta_xcb(ix)
         enddo
         ! apply boundary condition at the X-point, i=1: fixed density with specified ramp rate
         ydot(1) = density_ramp_rate
      ! write(*,*) 'ydot(density) =', ydot(0*Nx+1:1*Nx)
      ! ydot for the momentum equation
         ydot(Nx+1:2*Nx) = switch_momentum_source * Source_v(1:Nx)
         ! add the momentum flux term using the flux as calculated in calculate_fluxes
         ! Gamma_mom(i) contains the flux at i+1/2
         do ix = 2, Nx
            ydot(Nx+ix) = ydot(Nx+ix) - (Gamma_mom(ix)-Gamma_mom(ix-1))/delta_xcb(ix)
         enddo
         ! ! apply boundary condition at the X-point, as following from the constant density v(1) = n(2) v(2) / n(1)
         ! ydot(Nx+1) = ydot(Nx+1) - ???
         ! add the pressure term in the internal region using downwind differencing: NB pressure =2/3 * y(2*Nx+1:3*Nx)
         ydot(Nx+1) = ydot(Nx+1) - (y(2*Nx+2)-y(2*Nx+1))*energy_norm/1.5d+0/delta_x(1)
         do ix = 2, Nx-1
            ydot(Nx+ix) = ydot(Nx+ix) - (y(2*Nx+ix+1)-y(2*Nx+ix))*energy_norm/1.5d+0/delta_x(ix)
            ! add effect of numerical viscosity
            ydot(Nx+ix) = ydot(Nx+ix) + viscosity*(velocity(ix+1) + velocity(ix-1)-2.0d0*velocity(ix))
         enddo
         ! apply boundary condition at the sheath entrance, i=Nx: 
         ! assume zero temperture and extrapolate density and pressure to the sheath, add numerical viscosity (linearly extrapolate velocity beyond the sheath)
         ydot(2*Nx) = ydot(2*Nx) - (y(3*Nx)-y(3*Nx-1))*energy_norm/1.5d+0/delta_x(Nx)
         ydot(2*Nx) = ydot(2*Nx) + viscosity*(2.0d0*csound + velocity(Nx-1)-3.0d0*velocity(Nx))
      ! write(*,*) 'ydot(momentum) =', ydot(1*Nx+1:2*Nx)
      ! ydot for the energy equation
         ydot(2*Nx+1:3*Nx) = switch_energy_source * Source_Q(1:Nx)
         ! add the heat flux term in the internal region (including the sheath)
         ydot(2*Nx+2:3*Nx) = ydot(2*Nx+2:3*Nx) - (q_parallel(2:Nx)-q_parallel(1:Nx-1))/delta_xcb(2:Nx)
         ! add the compression term
         ydot(2*Nx+2:3*Nx) = ydot(2*Nx+2:3*Nx) + velocity(2:Nx) * (y(2*Nx+2:3*Nx)-y(2*Nx+1:3*Nx-1))*energy_norm/1.5d+0/delta_x(1:Nx-1)
         ! apply boundary condition at the X-point, i=1: energy flux given by q_parallel(0) = q_parX
         ydot(2*Nx+1) = ydot(2*Nx+1) - (q_parallel(1)-q_parX)/delta_xcb(1)
         ! limit ydot(2*Nx+1:3*Nx) to prevent temperature below minimum_temperature (assuming a time step of 1.0d-6)
         ! do ix = 1, Nx
         !    ! ydot(2*Nx+ix) = max(ydot(2*Nx+ix), (3.0d+0*density(ix)*e_charge*minimum_temperature - y(2*Nx+ix))/1.0d-9)
         !    if( 3.0d+0*density(ix)*e_charge*minimum_temperature .gt. y(2*Nx+ix) ) then
         !       ydot(2*Nx+ix) = max(ydot(2*Nx+ix), (3.0d+0*density(ix)*e_charge*minimum_temperature - y(2*Nx+ix))/1.0d-4)
         !       write(*,*) 'temperature limited', ix, time, y(2*Nx+ix)/density(ix)/e_charge/3.0d+0
         !    endif
         ! enddo
      ! write(*,*) 'ydot(energy) =', ydot(2*Nx+1:3*Nx)
      ! ydot for the neutral density equation
         ydot(3*Nx+1:4*Nx) = switch_neutral_source * neutral_source(1:Nx)
         ! add the density diffusion in the internal reagion
         ! the neutral particle diffusion coefficient D == n_n kT / m charge_exchange_rate sin^2theta
         do ix = 1, Nx
            Diff_neutral(ix) = D_neutral( temperature(ix), density(ix) )
         enddo
         ! write(*,*) 'Diff_neutral =', Diff_neutral
         ydot(3*Nx+2:4*Nx-1) = ydot(3*Nx+2:4*Nx-1) + Diff_neutral(2:Nx-1) * ((neutral(3:Nx)-neutral(2:Nx-1))/delta_x(2:Nx-1)-(neutral(2:Nx-1)-neutral(1:Nx-2))/delta_x(1:Nx-2))/delta_xcb(2:Nx-1)
         ydot(3*Nx+2:4*Nx-1) = ydot(3*Nx+2:4*Nx-1) + (Diff_neutral(3:Nx)-Diff_neutral(1:Nx-2))*(neutral(3:Nx)-neutral(1:Nx-2))/(delta_x(2:Nx-1)+delta_x(1:Nx-2))**2
         ! ydot(3*Nx+2:4*Nx) = ydot(3*Nx+2:4*Nx) - (neutral_flux(2:Nx)-neutral_flux(1:Nx-1))/delta_xcb(2:Nx)
         ! boundary condition at X-point (zero gradient i.e. at i=0 every equals i=1)
         ydot(3*Nx+1) = ydot(3*Nx+1) + Diff_neutral(1) * (neutral(2)-neutral(1))/delta_x(1)/delta_xcb(1)
         ydot(3*Nx+1) = ydot(3*Nx+1) + (Diff_neutral(2)-Diff_neutral(1))*(neutral(2)-neutral(1))/4.0d0/delta_x(1)/delta_x(1)
         ! ! boundary condition at sheath: neutral flux = - Gamma_n(Nx) * recycling * (1.0d-0 - redistributed_fraction)
         ydot(4*Nx) = ydot(4*Nx) + (Gamma_n(Nx) * recycling * (1.0d-0-redistributed_fraction) - 0.5d+0*(Diff_neutral(Nx)+Diff_neutral(Nx-1))*(neutral(Nx)-neutral(Nx-1))/delta_x(Nx-1))/delta_xcb(Nx)
         ! finally add the neutral sources and losses from redistribution and finite residence time
         ydot(3*Nx+1:4*Nx) = ydot(3*Nx+1:4*Nx) + Gamma_n(Nx) * recycling * redistributed_fraction / L - neutral / neutral_residence_time
      ! write(*,*) 'ydot(neutrals) =', ydot(3*Nx+1:4*Nx)
      ! apply evolution switches
      ydot(     1:  Nx) = evolve_density  * ydot(     1:  Nx) / density_norm
      ydot(  Nx+1:2*Nx) = evolve_momentum * ydot(  Nx+1:2*Nx) / momentum_norm
      ydot(2*Nx+1:3*Nx) = evolve_energy   * ydot(2*Nx+1:3*Nx) / energy_norm
      ydot(3*Nx+1:4*Nx) = evolve_neutral  * ydot(3*Nx+1:4*Nx) / density_norm
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

   real(wp) function D_neutral( temperature, density )
   ! function to calculate the neutral particle diffusion coefficient (Ref. Nakazawa et al. 2000 PPCF 42 401 equation 2.14)
      implicit none
      real(wp) :: temperature, density
      ! the neutral particle diffusion coefficient D == n_n kT / m charge_exchange_rate sin^2theta
      !                                               =     kT / m density <sigma v>_cx sin^2theta
      D_neutral = e_charge * temperature / (mass * density * charge_exchange(temperature) * sintheta**2)
      return
   end function D_neutral

   real(wp) function MC_limit( fm, fc, fp )
   ! function to implement MC slope limiter (van Leer 1977) as in SD1D
   ! the arguments are the values of the variable at i-1, i, i+1, respectively
      implicit none
      real(wp) :: fm, fc, fp
      MC_limit = minmod( (fp-fc), 0.5d+0*(fp-fm), (fc-fm) )
      return
   end function MC_limit

   real(wp) function minmod( a, b, c )
   ! function returning 0 when any input has different sign, otherwise minimum
      implicit none
      real(wp) :: a, b, c
      if( a*b .le. 0.0d+0 .or. a*c .le. 0.0d+0 ) then
         minmod = 0.0d+0
      else
         minmod = sign( min(abs(a),abs(b),abs(c)), a )
      endif
      return
   end function minmod

   SUBROUTINE JAC (NEQ, T, Y, ML, MU, PD, NROWPD)
      !dummy subroutine for calculation of Jacobian (dlsode option 21 or 24)
      INTEGER  NEQ, ML, MU, NROWPD
      DOUBLE PRECISION  T, Y(*), PD(NROWPD,*)
      RETURN
   END SUBROUTINE JAC
   
end module physics_routines
