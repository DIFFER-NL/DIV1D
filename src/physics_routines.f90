module physics_routines
! module containing general purpose routines implementing the equations

   use grid_data, only : delta_x, delta_xcb, x, B_field, B_field_cb, i_Xpoint
   use constants, only : e_charge
   use reaction_rates
   use physics_parameters, only : gamma, mass, Gamma_X, q_parX, energy_loss_ion, recycling, redistributed_fraction, L, neutral_residence_time, sintheta, minimum_density, minimum_temperature, density_ramp_rate, &
                                  gas_puff_source, gas_puff_location, gas_puff_width, initial_a, &
                                  dyn_nu, dyn_nb, dyn_dnu, dyn_gas, dyn_rec, dyn_rad_los, dyn_imp_con, gas_puff, dyn_red_frc, dyn_qparX, &
                                  L_core_SOL, X_core_SOL, alpha_core_profile, normalization_core_profile
                                 ! switch_X_vel_con
   use numerics_parameters, only : evolve_density, evolve_momentum, evolve_energy, evolve_neutral, switch_density_source, switch_momentum_source, switch_energy_source, switch_neutral_source, &
                                   switch_convective_heat, switch_impurity_radiation, viscosity, central_differencing, density_norm, momentum_norm, energy_norm, filter_sources,&
			   	                   delta_t
   use experiments, only: simulate_elm, calculate_radial_losses

   implicit none
   integer, parameter, private :: wp = KIND(1.0D0)
   !real( wp ), allocatable :: gas_puff(:)     ! vector holding the gas puff source [/m^3s]

contains

   subroutine nvt2y( Nx, density, velocity, temperature, neutral, y )
   ! subroutine to transform from (density, velocity, temperature) to (density, momentum, energy) in the solution vector y
   !    y(      1, ...,   Nx )     = log( density( 1, ..., Nx ) / density_norm )
   !    y(   Nx+1, ..., 2*Nx )     = Momentum = mass * density * velocity ( 1, ..., Nx ) / momentum_norm
   !    y( 2*Nx+1, ..., 3*Nx )     = internal energy = 3 * density * e_charge * temperature [eV] ( 1, ..., Nx ) / energy_norm
   !    y( 3*Nx+1, ..., 4*Nx )     = log( neutral density( 1, ..., Nx ) / density_norm )
      implicit none
      integer,  intent(in)  :: Nx
      real(wp), intent(in)  :: density(Nx), velocity(Nx), temperature(Nx), neutral(Nx)
      real(wp), intent(out) :: y(4*Nx)
      ! density
      y(     1:   Nx) = log(density(1:Nx) / density_norm)
      ! momentum
      y(  Nx+1: 2*Nx) = mass * density(1:Nx) * velocity(1:Nx) / momentum_norm
      ! energy
      y(2*Nx+1: 3*Nx) = 3.0d+0 * density(1:Nx) * e_charge * temperature(1:Nx) / energy_norm
      ! neutral density
      y(3*Nx+1: 4*Nx) = log(neutral(1:Nx) / density_norm)
      return
   end subroutine nvt2y


   subroutine y2nvt( Nx, y, density, velocity, temperature, neutral )
   ! subroutine to transform from the solution vector y containing (density, momentum, energy) to (density, velocity, temperature)
   ! note that the solution vector y is normalized
   !    y(      1, ...,   Nx )     = log( density( 1, ..., Nx ) / density_norm )
   !    y(   Nx+1, ..., 2*Nx )     = Momentum = mass * density * velocity ( 1, ..., Nx ) / momentum_norm
   !    y( 2*Nx+1, ..., 3*Nx )     = internal energy = 3 * density * e_charge * temperature [eV] ( 1, ..., Nx ) / energy_norm
   !    y( 3*Nx+1, ..., 4*Nx )     = log( neutral density( 1, ..., Nx ) / density_norm )
      implicit none
      integer,  intent(in)  :: Nx
      real(wp), intent(in)  :: y(4*Nx)
      real(wp), intent(out) :: density(Nx), velocity(Nx), temperature(Nx), neutral(Nx)
      ! density
      density(1:Nx)     =  exp(y(1:Nx)) * density_norm
      ! velocity
      velocity(1:Nx)    =  y(  Nx+1: 2*Nx) * momentum_norm / mass / density(1:Nx)
      ! temperature
      temperature(1:Nx) =  y(2*Nx+1: 3*Nx) * energy_norm / 3.0d+0 / density(1:Nx) / e_charge
      ! neutral density
      neutral(1:Nx)     =  exp(y(3*Nx+1: 4*Nx)) * density_norm
      return
   end subroutine y2nvt


   subroutine advection(Nx, variable, velocity, temperature, flux)
   ! this subroutine calculates the advected flux of a variable with a given velocity field
   ! the temperature is needed to calculate the sound velocity
   ! the variable and velocity field are given on the (equidistant) grid centers
   ! the flux is returned at the grid boundaries i + 1/2 (internal only: i=1,...,Nx-1, BC's must be applied outside this routine)
      implicit none
      integer,  intent(in)  :: Nx
      real(wp), intent(in)  :: variable(Nx), velocity(Nx), temperature(Nx)
      real(wp), intent(out) :: flux(0:Nx)
      real(wp)              :: csound, average_velocity
      integer               :: i
      flux = 0.0
      ! we follow here the discretization as put forward in B. Dudson et al. (2019) PPCF 61 065008
      do i = 1, Nx-1
         csound = sqrt( 2.0d+0 * e_charge * max(temperature(i),temperature(i+1),minimum_temperature) / mass )
         ! csound = sqrt(          e_charge *    (temperature(i)+temperature(i+1)) / mass )
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


   subroutine calculate_fluxes( Nx, time,  density, velocity, temperature, neutral, Gamma_n, Gamma_mom, q_parallel, neutral_flux, elm_heat_load )
   ! this subroutine calculates the 'fluxes' required for the calculation of ydot (the right hand side of the discretized conservation equations)
   ! the fluxes are defined at the halfway point between grid points: i.e. Flux(i) is midway between x(i) and x(i+1) on the cell boundaries xcb
   ! ew 01-03-2021: modified to take account of flux expansion
      implicit none
      integer,  intent(in)  :: Nx
      real(wp), intent(in)  :: time
      integer               :: itime
      real(wp), intent(in)  :: density(Nx), velocity(Nx), temperature(Nx), neutral(Nx)
      real(wp), intent(out) :: Gamma_n(0:Nx), Gamma_mom(0:Nx), q_parallel(0:Nx), neutral_flux(0:Nx)
      real(wp), intent(in)  :: elm_heat_load
      real(wp)              :: momentum(Nx), enthalpy(Nx)
      real(wp)              :: csound_target(2), average_velocity
      integer               :: i
      ! the particle flux = density velocity / B 
      itime = time / delta_t

      ! the particle flux = density velocity
      ! we follow here the discretization as put forward in B. Dudson et al. (2019) PPCF 61 065008
         call advection(Nx, density/B_field, velocity, temperature, Gamma_n)
         ! boundary condition at the sheath (note that velocity is allowed to be supersonic)
         csound_target(2) = sqrt( 2.0d+0 * e_charge * max(1.5d+0*temperature(Nx)-0.5d+0*temperature(Nx-1),minimum_temperature) / mass )
         ! write(*,*) 'flux', e_charge, mass, temperature(Nx)
         ! Gamma_n(Nx) = (1.5*density(Nx)-0.5*density(Nx-1)) * max(velocity(Nx),csound_target(2)) / B_field_cb(Nx)
         Gamma_n(Nx) = min(density(Nx),(1.5*density(Nx)-0.5*density(Nx-1))) * max(velocity(Nx),csound_target(2)) / B_field_cb(Nx)
         ! boundary condition at i = 0 (not used when X-point density BC is used)
         if( X_core_SOL .eq. 0.0d+0 ) then
             ! case for core-SOL + divertorleg (i=0 is stagnation point)
             Gamma_n(0)  = 0.0d0
         else
             ! case with two targets so sheath boundary conditions at i=0 as well
             csound_target(1) = - sqrt( 2.0d+0 * e_charge * max(1.5d+0*temperature(1)-0.5d+0*temperature(2),minimum_temperature) / mass ) ! note the sign
             ! Gamma_n(0) = (1.5*density(1)-0.5*density(2)) * max(abs(velocity(Nx)),csound_target(1)) / B_field_cb(0)
             Gamma_n(0) = min(density(1),(1.5*density(1)-0.5*density(2))) * min(velocity(1),csound_target(1)) / B_field_cb(0)
         endif
         ! write(*,*) 'flux', density(Nx), density(Nx-1), velocity(Nx), csound_target
         ! if(temperature(Nx) .le. minimum_temperature) Gamma_n(Nx) = 0.0d+0

      ! the momentum flux = momentum * velocity where momentum = density * mass * velocity / B_field
         momentum = density * mass * velocity / B_field
         call advection(Nx, momentum, velocity, temperature, Gamma_mom)
         ! boundary condition at the sheath
         ! Gamma_mom(Nx) = (1.5*density(Nx)-0.5*density(Nx-1)) * mass * max(velocity(Nx),csound_target(2))**2 / B_field_cb(Nx)
         Gamma_mom(Nx) = min(density(Nx),(1.5*density(Nx)-0.5*density(Nx-1))) * mass * max(velocity(Nx),csound_target(2))**2 / B_field_cb(Nx)
         ! if(temperature(Nx) .le. minimum_temperature) Gamma_mom(Nx) = 0.0d+0
         ! boundary condition at i = 0
         if( X_core_SOL .eq. 0.0d+0 ) then
             ! case for core-SOL + divertorleg (i=0 is stagnation point)
             Gamma_mom(0)  = 0.0d0
         else
             ! case with two targets so sheath boundary conditions at i=0 as well
             Gamma_mom(0) = min(density(1),(1.5*density(1)-0.5*density(2))) * mass * min(velocity(1),csound_target(1))**2 / B_field_cb(0)
         endif

      ! convective heat flux = 5 density k temperature velocity / B_field (i.e. 5/2 pressure)
         enthalpy = 5.0d+0 * density * e_charge * temperature / B_field
         call advection(Nx, enthalpy, velocity, temperature, q_parallel) 
         ! write(*,*) 'advection', q_parallel
      ! add the conductive heat flux in the internal region
         do i = 1, Nx-1
            q_parallel(i) = q_parallel(i) * switch_convective_heat - kappa_parallel(0.5d+0*(temperature(i)+temperature(i+1))) * (temperature(i+1)-temperature(i))/delta_x(i) / B_field_cb(i)
         enddo
         ! boundary condition at the sheath: given by the sheath heat transmission
         q_parallel(Nx) = gamma * csound_target(2) * min(density(Nx),(1.5d+0*density(Nx)-0.5d+0*density(Nx-1))) * e_charge * max(1.5d+0*temperature(Nx)-0.5d+0*temperature(Nx-1),minimum_temperature) / B_field_cb(Nx) ! we have extrapolated the density linear towards x = L, i.e. the sheath
         ! if(temperature(Nx) .le. minimum_temperature .or. q_parallel(Nx) .lt. 0.0d+0) q_parallel(Nx) = 0.0d+0
         ! boundary condition at i = 0
         if( L_core_SOL .eq. 0.0 ) then
             ! case with prescribed heat flux at i = 0 (X-point)
             q_parallel(0)  = dyn_qparX(itime)+elm_heat_load
         elseif( X_core_SOL .eq. 0.0 ) then
             ! case for core-SOL + divertorleg (i=0 is stagnation point)
             q_parallel(0)  = 0.0d0
         else
             ! case with two targets so sheath boundary conditions at i=0 as well
             q_parallel(0) = gamma * csound_target(1) * min(density(1),(1.5d+0*density(1)-0.5d+0*density(2))) * e_charge * max(1.5d+0*temperature(1)-0.5d+0*temperature(2),minimum_temperature) / B_field_cb(0) ! we have extrapolated the density linear towards x = 0, i.e. the sheath
         endif
         
      ! the neutral particle diffusion !!!! switch-on in case you want this diagnostic
         ! we do this in the right_hand_side routine itself
         do i = 1, Nx-1
            neutral_flux(i) = - 0.5d+0*(D_neutral(temperature(i),density(i))+D_neutral(temperature(i+1),density(i+1))) * (neutral(i+1)-neutral(i))/delta_x(i)
         enddo
         ! boundary condition at i = 0
         if( L_core_SOL .eq. 0.0d+0 ) then
             ! boundary condition at X-point (zero gradient i.e. at i=0 every equals i=1, i.e. zero flux)
             neutral_flux(0) = 0.0d+0
         elseif( X_core_SOL .eq. 0.0d+0 ) then
             ! apply boundary condition at mid point (i.e. flux = 0)
             neutral_flux(0) = 0.0d+0
         else
             ! sheath boundary condition (= flux of plasma density in case of full recycling)
             neutral_flux(0) = -B_field_cb(0)*Gamma_n(0) * dyn_rec(itime) * (1.0d-0 - dyn_red_frc(itime))
         endif
         ! boundary condition at the sheath (= flux of plasma density in case of full recycling)
         neutral_flux(Nx) = -B_field_cb(Nx)*Gamma_n(Nx) * dyn_rec(itime) * (1.0d-0 - dyn_red_frc(itime))
         ! write(*,*) 'temperature =', temperature
         ! write(*,*) 'q_parallel =', q_parallel

      return
   end subroutine calculate_fluxes


   subroutine initialize_gas_puff(Nx)
   ! this routine needs to be updated if to be used for a case with two targets
      implicit none
      integer,  intent(in)  :: Nx
      real(wp) gas_puff_normalization
      
      ! allocate and initialize the gas puff to zero
      allocate( gas_puff(Nx) )
      gas_puff = 0.0
      if(gas_puff_source .eq. 0.0d+0) return
      
      ! calculate the Gaussian profile of the source and the corresponding normalization factor
      gas_puff = exp( - (x - gas_puff_location)**2 / gas_puff_width**2 )
      gas_puff_normalization = sum(gas_puff * delta_xcb)
      
      ! set the unnormalized gas puff source
      ! gas_puff = gas_puff_source * gas_puff / gas_puff_normalization
      gas_puff = gas_puff / gas_puff_normalization
      return
   end subroutine initialize_gas_puff


   subroutine calculate_sources( Nx, time, density, velocity, temperature, neutral, q_parallel,&
                                 Source_n, Source_v, Source_Q, neutral_source, elm_heat_load, elm_density_change )
   ! this subroutine calculates the source terms of the discretized conservation equations
      implicit none
      integer,  intent(in)  :: Nx
      real(wp), intent(in)  :: time
      integer               :: itime
      real(wp), intent(in)  :: density(Nx), velocity(Nx), temperature(Nx), neutral(Nx),q_parallel(Nx)
      real(wp), intent(out) :: Source_n(Nx), Source_v(Nx), Source_Q(Nx), neutral_source(Nx)
      real(wp) :: rate_cx(Nx), rate_ion(Nx), rate_exc(Nx), rate_rec(Nx), rate_ree(Nx), rate_imp(Nx)
      real(wp) :: radial_sink(Nx)
      real(wp) :: mid_point
      integer  :: ix, iix
      ! input variables for the elm simulation
      real(wp), intent(in)   :: elm_heat_load, elm_density_change
      itime     = time / delta_t  
      Source_n = 0.0d+0
      Source_v = 0.0d+0
      Source_Q = 0.0d+0
      neutral_source = 0.0d+0
      do ix = 1, Nx
         rate_cx(ix)  = density(ix) * neutral(ix) * charge_exchange(temperature(ix))
         rate_ion(ix) = density(ix) * neutral(ix) * ionization(density(ix),temperature(ix))
         rate_exc(ix) = density(ix) * neutral(ix) * excitation(density(ix),temperature(ix))
         rate_rec(ix) = density(ix) * density(ix) * recombination(density(ix),temperature(ix))
         rate_ree(ix) = density(ix) * density(ix) * recombenergy(density(ix),temperature(ix))
         rate_imp(ix) = density(ix) * density(ix) * impurity_radiation(temperature(ix),itime)
      enddo
      ! the particle sources
      neutral_source = rate_rec - rate_ion + gas_puff * dyn_gas(itime)
      Source_n = rate_ion - rate_rec
      ! the momentum sources
      Source_v = - mass * velocity * ( rate_cx + rate_rec )
      ! the energy sources (only internal energy)
      Source_Q = - (1.5d+0 * e_charge * temperature) * (rate_cx + 2.0d+0*rate_rec) ! the factor 2 for the recombination rate accounts for the loss of both an ion and an electron
      ! add energy source associated with particle source from ionization (like a source of internal energy from effective friction with neutrals)
      Source_Q = Source_Q + ( rate_ion * 0.5 * mass * velocity**2 )
      if ( switch_excitation .eq. 0.0d+0 ) then
         Source_Q = Source_Q - rate_ion * e_charge * energy_loss_ion ! note energy loss per ionization is in eV
      else
         Source_Q = Source_Q - switch_excitation * rate_exc * e_charge ! note excitation rate is in eV m^3 / s
      endif
      ! add energy losses associated with radiative and 3-body recombination (the 13.6 eV potential energy per recombination event is added explicitly here)
      Source_Q = Source_Q - ( rate_ree - 13.6 * rate_rec ) * e_charge ! note recombination loss rate is in eV m^3 / s
      ! add impurity radiation losses
      Source_Q = Source_Q - switch_impurity_radiation * rate_imp * e_charge ! note impurity radiation loss rate is in eV m^3 / s
      ! write(*,*) rate_ion, Source_Q
      ! Add the effect of radial losses across the flux tube
      call calculate_radial_losses(Nx,radial_sink,q_parallel, time)
      ! Consecutively, check whether substracting this radial_sink does not yield unphysical results by confirming that the total 
      ! losses over the flux tube are smaller than the incoming flux, so as not to obtain sub-zero fluxes.
      do ix = 1, Nx
         if (abs(Source_Q(ix) - radial_sink(ix)).ge.0.9d+0*abs(q_parallel(ix)/delta_xcb(ix))) then
            do iix = ix, Nx
               radial_sink(iix) = 0
            end do
            exit
         endif
      enddo
      Source_Q = Source_Q - radial_sink ! note impurity radiation loss rate is in eV m^3 / s
      ! remove spikes that cause problems during integration of the ODE
      if( filter_sources ) then
          call spike_filter( Nx, neutral_source )
          call spike_filter( Nx, Source_n )
          call spike_filter( Nx, Source_v )
          call spike_filter( Nx, Source_Q )
      endif
      ! section adding sources along the core-SOL boundary (when present in the grid)
      if( L_core_SOL .gt. 0.0d+0 .and. X_core_SOL .eq. 0.0d+0 ) then
          Source_n(i_Xpoint(1):i_Xpoint(2)) = Source_n(i_Xpoint(1):i_Xpoint(2)) + Gamma_X * (1.0d+0 - (x(i_Xpoint(1):i_Xpoint(2))/L_core_SOL)**2)**alpha_core_profile / normalization_core_profile
          Source_Q(i_Xpoint(1):i_Xpoint(2)) = Source_Q(i_Xpoint(1):i_Xpoint(2)) + (dyn_qparX(itime)+elm_heat_load)  * (1.0d+0 - (x(i_Xpoint(1):i_Xpoint(2))/L_core_SOL)**2)**alpha_core_profile / normalization_core_profile
      elseif( L_core_SOL .gt. 0.0d+0 .and. X_core_SOL .gt. 0.0d+0 ) then
          mid_point = X_core_SOL + 0.5d+0*L_core_SOL
          Source_n(i_Xpoint(1):i_Xpoint(2)) = Source_n(i_Xpoint(1):i_Xpoint(2)) + Gamma_X * (1.0d+0 - 4.0d+0*((x(i_Xpoint(1):i_Xpoint(2))-mid_point)/L_core_SOL)**2)**alpha_core_profile / normalization_core_profile
          Source_Q(i_Xpoint(1):i_Xpoint(2)) = Source_Q(i_Xpoint(1):i_Xpoint(2)) + (dyn_qparX(itime)+elm_heat_load)  * (1.0d+0 - 4.0d+0*((x(i_Xpoint(1):i_Xpoint(2))-mid_point)/L_core_SOL)**2)**alpha_core_profile / normalization_core_profile
      endif
      return
   end subroutine calculate_sources


   subroutine right_hand_side( neq, time, y, ydot )
   ! this subroutine calculates the right hand side ydot of the discretized conservation equations
   ! note that y and ydor are normalized, but the arrays density, velocity, temperature, and neutral are not!
   ! ew 01-03-2021: modified to take account of flux expansion
      implicit none
      integer,  intent(in)  :: neq
      real(wp), intent(in)  :: time, y(neq) 
      integer               :: itime 
      real(wp), intent(out) :: ydot(neq)
      integer               :: Nx, ix
      real(wp)              :: density(neq/4), velocity(neq/4), temperature(neq/4), neutral(neq/4)      ![1/m3] ,[m/s]    ,[eV]   ,[1/m3]
      real(wp)              :: Gamma_n(0:neq/4), Gamma_mom(0:neq/4), q_parallel(0:neq/4), neutral_flux(0:neq/4) ![1/m2s],[kg/ms2] ,[J/m2s],[1/m2s]
      real(wp)              :: Source_n(neq/4), Source_v(neq/4), Source_Q(neq/4), neutral_source(neq/4) ![1/m3s],[kg/m2s2],[J/m3s],[1/m3s] 
      real(wp)              :: Diff_neutral(neq/4)  ! [m2/s]
      real(wp)              :: csound_target(2), q_sheath, v0, Gmom0
      ! input variables for the elm simulation
      real(wp)              :: elm_heat_load, elm_density_change
      Nx        = neq/4
      itime     = time / delta_t
      ! write(*,*) 'RHS called at itime =', itime
      ! write(*,*) 'RHS called at t =', time
      ! write(*,*) 'y =', y
      ! first tranform the solution vector to (density velocity temperature neutral-density)
      call y2nvt( Nx, y, density, velocity, temperature, neutral )
      ! write(*,*) 'density = ', density
      ! write(*,*) 'velocity =', velocity
      ! write(*,*) 'temperature =', temperature
      ! write(*,*) 'neutral density =', neutral

      ! calculate the ELM heat flux and particle flux
      call simulate_elm(elm_heat_load, elm_density_change, time)
      ! calculate the fluxes
      call calculate_fluxes( Nx, time,  density, velocity, temperature, neutral, Gamma_n, Gamma_mom, q_parallel, neutral_flux, elm_heat_load )
      ! calculate the sources
      call calculate_sources( Nx, time, density, velocity, temperature, neutral, q_parallel, Source_n, Source_v, Source_Q, neutral_source, elm_heat_load, elm_density_change )
      ! write(*,*) 'Gamma_n =', Gamma_n
      ! write(*,*) 'Gamma_mom =', Gamma_mom
      ! write(*,*) 'q_parallel =', q_parallel
      ! write(*,*) 'Source_n =', Source_n
      ! write(*,*) 'Source_v =', Source_v
      ! write(*,*) 'Source_Q =', Source_Q
      ! write(*,*) 'neutral_source =', neutral_source

      ! sound velocity at the target(s)
      csound_target(2) = sqrt( 2.0d+0 * e_charge * max(1.5d+0*temperature(Nx)-0.5d+0*temperature(Nx-1),minimum_temperature) / mass )
      if(X_core_SOL .gt. 0.0d+0) csound_target(1) = - sqrt( 2.0d+0 * e_charge * max(1.5d+0*temperature(1)-0.5d+0*temperature(2),minimum_temperature) / mass ) ! note the sign

      ! -------------------------------------------- ydot for the density equation ----------------------------------------------
         ydot(1:Nx) = switch_density_source * Source_n(1:Nx) ! [1/ (m^3 s) ]
         ! add the particle flux term using the flux as calculated in calculate_fluxes
         ! Gamma_n(i) contains the flux at i+1/2
         do ix = 2, Nx
            ydot(ix) = ydot(ix) - B_field(ix) * (Gamma_n(ix)-Gamma_n(ix-1))/delta_xcb(ix)   ! ew 01-03-2021:
         enddo
         if( L_core_SOL .gt. 0.0d+0 ) then
             ! apply boundary condition at x = 0 as calculated in calculate_fluxes
             ydot(1) = ydot(1) - B_field(1) * (Gamma_n(1)-Gamma_n(0))/delta_xcb(1)
         else
             ! apply boundary condition at the X-point, i=1: fixed density with specified ramp rate, elm or perturbation in time	 
             ydot(1) = density_ramp_rate + elm_density_change + dyn_dnu(itime) ! [1/ (m^3 s)]
         endif
      ! write(*,*) 'ydot(density) =', ydot(0*Nx+1:1*Nx) ! ---------------------------------------------------------------------

      ! --------------------------------------------- ydot for the momentum equation ------------------------------------------
         ydot(Nx+1:2*Nx) = switch_momentum_source * Source_v(1:Nx) ![kg/(m^2 s^2]
         ! add the momentum flux term using the flux as calculated in calculate_fluxes
         ! Gamma_mom(i) contains the flux at i+1/2
         do ix = 2, Nx
            ydot(Nx+ix) = ydot(Nx+ix) - B_field(ix) * (Gamma_mom(ix)-Gamma_mom(ix-1))/delta_xcb(ix)   ! ew 01-03-2021:
         enddo
         if( L_core_SOL .gt. 0.0d+0 ) then
             ! apply boundary condition at mid point (i.e. flux = 0) or as calculated for the target at x=0
             ydot(Nx+1) = ydot(Nx+1) - B_field(1) * (Gamma_mom(1)-Gamma_mom(0))/delta_xcb(1)    ! [kg/(m^2 s^2)]
         else
             ! apply boundary condition at the X-point, as following from the constant density n(1)
             ! velocity at i = 0:  v(0) = v(2) - 2 (S_n(1)+density_ramp_rate) delta_xcb(1) / n(1)
             !v0 = velocity(2) - 2.0d+0 * (Source_n(1)+density_ramp_rate) * delta_xcb(1) / density(1) ! GD: no ELM_density_change ?
             v0 = velocity(2) - 2.0d+0 * (Source_n(1)+ydot(1)) * delta_xcb(1) / density(1) ! GD added result from density equation
	         ! momentum flux at i = 0 : Gmom0 = m n (1/4)(v(0)+v(1))**2
             Gmom0 = mass * density(1) * (v0 + velocity(1))**2/4.0d+0       ! [kg/(m   s^2)]
             ydot(Nx+1) = ydot(Nx+1) - (Gamma_mom(1)-Gmom0)/delta_xcb(1)    ! [kg/(m^2 s^2)]
         endif
         ! add the pressure term in the internal region using downwind differencing: NB pressure =2/3 * y(2*Nx+1:3*Nx)
         do ix = 2, Nx-1
            ydot(Nx+ix) = ydot(Nx+ix) - (1.0d+0-central_differencing)*(y(2*Nx+ix+1)-y(2*Nx+ix))*energy_norm/1.5d+0/delta_x(ix)&
                                      - central_differencing*(y(2*Nx+ix+1)-y(2*Nx+ix-1))*energy_norm/1.5d+0/(x(ix+1)-x(ix-1))
            ! add effect of numerical viscosity
            ydot(Nx+ix) = ydot(Nx+ix) + viscosity*(velocity(ix+1) + velocity(ix-1)-2.0d0*velocity(ix))
         enddo
         !add pressure term at the boundaries
         ydot(Nx+1) = ydot(Nx+1) - (y(2*Nx+2)-y(2*Nx+1))*energy_norm/1.5d+0/delta_x(1)
         ydot(2*Nx) = ydot(2*Nx) - (y(3*Nx)  -y(3*Nx-1))*energy_norm/1.5d+0/delta_x(Nx-1)
         ! apply boundary condition at the sheath entrance, i=Nx: (linearly extrapolate velocity beyond the sheath)
         ydot(2*Nx) = ydot(2*Nx) + viscosity*(2.0d0*csound_target(2) + velocity(Nx-1)-3.0d0*velocity(Nx))  ! add numerical viscocity
         if( L_core_SOL .gt. 0.0d+0 .and. X_core_SOL .gt. 0.0d+0 ) then
             ! apply boundary contion at the x=0 sheath for the pressure and viscosity terms
             ydot(Nx+1) = ydot(Nx+1) + viscosity*(2.0d0*csound_target(1) + velocity(2)-3.0d0*velocity(1))  ! add numerical viscocity
         endif
      ! write(*,*) 'ydot(momentum) =', ydot(1*Nx+1:2*Nx) ! -----------------------------------------------------------------------
      
      ! ------------------------------------------------ ydot for the energy equation --------------------------------------------
         ydot(2*Nx+1:3*Nx) = switch_energy_source * Source_Q(1:Nx) ! [J/ (m^3 s)]
         ! add the heat flux term including all boundaries
         ydot(2*Nx+1:3*Nx) = ydot(2*Nx+1:3*Nx) - B_field(1:Nx) * (q_parallel(1:Nx)-q_parallel(0:Nx-1))/delta_xcb(1:Nx)   ! ew 01-03-2021:
         ! add the compression term (we symmetrize this is the internal reagion)
         ydot(2*Nx+1)        = ydot(2*Nx+1)        +        velocity(1)      * (y(2*Nx+2)       -y(2*Nx+1)       )*energy_norm/1.5d+0/delta_x(1)
         ydot(2*Nx+2:3*Nx-1) = ydot(2*Nx+2:3*Nx-1) + 0.5d+0*velocity(2:Nx-1) * (y(2*Nx+3:3*Nx  )-y(2*Nx+2:3*Nx-1))*energy_norm/1.5d+0/delta_x(2:Nx-1)
         ydot(2*Nx+2:3*Nx-1) = ydot(2*Nx+2:3*Nx-1) + 0.5d+0*velocity(2:Nx-1) * (y(2*Nx+2:3*Nx-1)-y(2*Nx+1:3*Nx-2))*energy_norm/1.5d+0/delta_x(1:Nx-2)
         ydot(3*Nx)          = ydot(3*Nx)          +        velocity(Nx)     * (y(3*Nx)         -y(3*Nx-1)       )*energy_norm/1.5d+0/delta_x(Nx-1)
      ! write(*,*) 'ydot(energy) =', ydot(2*Nx+1:3*Nx) ! -------------------------------------------------------------------------

      ! ----------------------------------------------ydot for the neutral density equation --------------------------------------
         ydot(3*Nx+1:4*Nx) = switch_neutral_source * neutral_source(1:Nx) ![1/ (m^3 s)]
         ! add the density diffusion in the internal reagion
         ! the neutral particle diffusion coefficient D == n_n kT / m charge_exchange_rate sin^2theta
         do ix = 1, Nx
            Diff_neutral(ix) = D_neutral( temperature(ix), density(ix) )
         enddo
         ! write(*,*) 'Diff_neutral =', Diff_neutral
         ydot(3*Nx+1:4*Nx) = ydot(3*Nx+1:4*Nx) - (neutral_flux(1:Nx)-neutral_flux(0:Nx-1))/delta_xcb(1:Nx)
         ! add neutral sources and losses from redistribution and finite residence time
         ydot(3*Nx+1:4*Nx) = ydot(3*Nx+1:4*Nx) + B_field_cb(Nx) * Gamma_n(Nx) * dyn_rec(itime) * dyn_red_frc(itime) / L &
                                                                                       - (neutral-dyn_nb(itime)) / neutral_residence_time
      ! write(*,*) 'ydot(neutrals) =', ydot(3*Nx+1:4*Nx) !-------------------------------------------------------------------------

      ! apply evolution switches
         ydot(     1:  Nx) = evolve_density  * ydot(     1:  Nx) / density(1:Nx)
	     ydot(  Nx+1:2*Nx) = evolve_momentum * ydot(  Nx+1:2*Nx) / momentum_norm
      	 ydot(2*Nx+1:3*Nx) = evolve_energy   * ydot(2*Nx+1:3*Nx) / energy_norm
      	 ydot(3*Nx+1:4*Nx) = evolve_neutral  * ydot(3*Nx+1:4*Nx) / neutral(1:Nx)
      ! write(*,*) 'time =', time, 'ydot(Nx) =', ydot(Nx), 'density(Nx) =', density(Nx), 'Source_n(Nx) =', Source_n(Nx), Gamma_n(Nx)
      ! write(*,*) 'ydot =', ydot
      return
   end subroutine right_hand_side

   real(wp) function kappa_parallel( temperature )
   ! function to calculate the parallel heat conductivity
      implicit none
      real(wp) :: temperature
      ! use expression from Stangeby page 187 (Chapter 4.10.1)
      kappa_parallel = 2.0d+3 * max(temperature,minimum_temperature)*max(temperature,minimum_temperature)*sqrt(max(temperature,minimum_temperature))
      return
   end function kappa_parallel

   real(wp) function D_neutral( temperature, density )
   ! function to calculate the neutral particle diffusion coefficient (Ref. Nakazawa et al. 2000 PPCF 42 401 equation 2.14)
      implicit none
      real(wp) :: temperature, density
      ! the neutral particle diffusion coefficient D == n_n kT / m charge_exchange_rate sin^2theta
      !                                               =     kT / m density <sigma v>_cx sin^2theta
      ! D_neutral = e_charge * max(temperature,1.0d+0,minimum_temperature) / (mass * density * charge_exchange(temperature) * sintheta**2)
      D_neutral = e_charge * max(temperature,minimum_temperature) / (mass * density * charge_exchange(temperature) * sintheta**2)
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
   
   subroutine spike_filter(N, array)
   ! simple routine to remove spikes in data arrays
      implicit none
      integer, intent(in)     :: N
      real(wp), intent(inout) :: array(N)
      integer                 :: i
      real(wp)                :: array_input(N)
      array_input = array
      do i = 2, N-1
         array(i) = midvalue( array_input(i-1), array_input(i), array_input(i+1) )
      enddo
      array(1) = midvalue( array_input(1), array_input(2), array_input(3) )
      array(N) = midvalue( array_input(N), array_input(N-1), array_input(N-2) )
      return
   end subroutine spike_filter
   
   
   real(wp) function midvalue( a, b, c )
   ! function returning the median value
      implicit none
      real(wp) :: a, b, c
      if( (a-b)*(b-c) .ge. 0 ) then
         midvalue = b
      elseif( abs(a-b) .le. abs(b-c) ) then
         midvalue = a
      else
         midvalue = c
      endif
      return
   end function midvalue
   

   SUBROUTINE JAC (NEQ, T, Y, ML, MU, PD, NROWPD)
      !dummy subroutine for calculation of Jacobian (dlsode option 21 or 24)
      INTEGER  NEQ, ML, MU, NROWPD
      DOUBLE PRECISION  T, Y(*), PD(NROWPD,*)
      RETURN
   END SUBROUTINE JAC
   
end module physics_routines
