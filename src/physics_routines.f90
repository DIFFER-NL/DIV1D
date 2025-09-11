module physics_routines
! module containing general purpose routines implementing the equations

! License Notice
! This file is part of DIV1D.
! DIV1D is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License 
! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
! DIV1D is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License along with DIV1D. 
! If not, see <https://www.gnu.org/licenses>. 


   use grid_data, only : delta_x, delta_xcb, x, B_field, B_field_cb, &
			i_Xpoint, i_baffle, Area_extern, volumes, gas_puff_profile, &
			I_core_source_profile_Q, I_core_source_profile_n, cc2cb, A_wet, core_source_profile_n, A_int
   use constants, only : e_charge
   use reaction_rates
   use physics_parameters, only : gamma, mass, Gamma_X, q_parX, Gamma_core, Q_core, energy_loss_ion, recycling,	atom_recycle_energy_fraction, &
 				redistributed_fraction, L, neutral_residence_time, molecule_residence_time, &
				sintheta, minimum_density, maximum_density, minimum_temperature, density_ramp_rate, wide_PFR, sigma_nb,&
                                gas_puff_source, gas_puff_location, gas_puff_width, initial_a, &
				mol_rec, D_molecule_puff, D_neutral_puff, pump_rate_n, pump_rate_m, &
				L_core_SOL, X_core_SOL, core_far_sol_ion_loss, core_far_sol_feedthrough, &
				core_ext_neutral_pump, core_ext_molecule_pump, &
                                extern_neutral_ex,  extern_molecule_ex, extern_neutral_volumes, pol_target_angle, &
				neutral_energy, extern_neutral_wall_area, wall_association_probability, &
				core_confinement_time, core_volume , core_sol_neutral_ex, core_sol_molecule_ex , core_ionization_fraction,core_ionization_fraction_mol, &
				D_imp_con, D_neu, D_dneu, D_nb, D_mb, D_gas, D_rec, D_qpar_x, D_red_frc, D_Q_core, &
				D_Gamma_core, D_core_fuelling, internal_istep, D_start_time, D_delta_t, D_nout, &
				D_core_neutral_density, ato_res_asy, mol_res_asy
				! if you want to set other parameters to use as interface:
				! D_extern2sol_flux, D_extern2sol_mol ! for coupling with other reservoir models
           
   use numerics_parameters, only : Nx, evolve_density, evolve_momentum, evolve_energy, evolve_neutral, 		&
				 evolve_neutral_momentum, evolve_molecule, evolve_background, evolve_core, evolve_core_neutral,	&
				 switch_density_source, switch_momentum_source, switch_energy_source, switch_energy_flux, switch_energy_compr,	switch_only_core_source_Q,	&
				  		 					 switch_neutral_source, switch_neutral_momentum_source, switch_molecule_source, 	&
                                 switch_convective_heat, switch_impurity_radiation, 				&
				 viscosity, central_differencing, density_norm, 				&
                                 momentum_norm, energy_norm, neutral_norm, filter_sources, delta_t, lax_switch, &
				nout, D_new, wide_core_profile, detect_nan, max_attempts, switch_neutral_leakage 
 
   ! the only plasma params value that is used here ? if this is an input, then it should be prescribed
   ! use plasma_params, only: core_neutral_density 

   use experiments, only: simulate_elm !, calculate_radial_losses
   !use dvode_f90_m ! for VODE OPTS
   !use div1d_step, only:  


   implicit none
   integer, parameter, private :: wp = KIND(1.0D0)     
   integer :: internal_istep_rhs = 1

   ! add plasma data here

   !type (VODE_OPTS)      :: options 
contains

  subroutine nvt2ys( Nx, density, velocity, temperature, neutral, neutral_velocity, molecule, ys )
   ! subroutine to transform from (density, velocity, temperature) to (density, momentum, energy) in the solution vector y
   !    y(      1, ...,   Nx )     = log( density( 1, ..., Nx ) / density_norm )
   !    y(   Nx+1, ..., 2*Nx )     = Momentum = mass * density * velocity ( 1, ..., Nx ) / momentum_norm
   !    y( 2*Nx+1, ..., 3*Nx )     = internal energy = 3 * density * e_charge * temperature [eV] ( 1, ..., Nx ) / energy_norm
   !    y( 3*Nx+1, ..., 4*Nx )     = log( neutral density( 1, ..., Nx ) / density_norm )
   !	y( 4*Nx+1, ..., 5*Nx )	   = Neutral momentum = mass * density * neutral_velocity / momentum_norm
      implicit none
      integer,  intent(in)  :: Nx
      real(wp), intent(in)  :: density(Nx), velocity(Nx), temperature(Nx), neutral(Nx), neutral_velocity(Nx), molecule(Nx)
      !real(wp), intent(in) :: extern_neutral_density(5), extern_molecule_density(5)
      real(wp), intent(out) :: ys(6*Nx)
      ys(     1:   Nx) = log(density(1:Nx) / density_norm)
      ! momentum
      ys(  Nx+1: 2*Nx) = mass * density(1:Nx) * velocity(1:Nx) / momentum_norm
      ! energy
      ys(2*Nx+1: 3*Nx) = 3.0d+0 * density(1:Nx) * e_charge * temperature(1:Nx) / energy_norm
      ! neutral density
      ys(3*Nx+1: 4*Nx) = log(neutral(1:Nx) / density_norm)
      ! neutral momentum
      ys(4*Nx+1: 5*Nx) = mass * neutral(1:Nx) * neutral_velocity(1:Nx) / momentum_norm
      ! molecule density
      ys(5*Nx+1: 6*Nx) = log(molecule(1:Nx) / density_norm)
	!write(*,*) 'neutral', neutral
	!write(*,*) 'density norm', density_norm
      return
   end subroutine nvt2ys

   subroutine ys2nvt( Nx, ys, density, velocity, temperature, neutral, neutral_velocity, molecule )
   ! subroutine to transform from the solution vector y containing (density, momentum, energy) to (density, velocity, temperature)
   ! note that the solution vector y is normalized
   !    y(      1, ...,   Nx )     = log( density( 1, ..., Nx ) / density_norm )
   !    y(   Nx+1, ..., 2*Nx )     = Momentum = mass * density * velocity ( 1, ..., Nx ) / momentum_norm
   !    y( 2*Nx+1, ..., 3*Nx )     = internal energy = 3 * density * e_charge * temperature [eV] ( 1, ..., Nx ) / energy_norm
   !    y( 3*Nx+1, ..., 4*Nx )     = log( neutral density( 1, ..., Nx ) / density_norm )
   !	y( 4*Nx+1, ..., 5*Nx )	   = Neutral momentum = mass * density * neutral_velocity / momentum_norm
      implicit none
      integer :: i
      real(wp) :: tmp, tmp1, tmpn
      integer,  intent(in)  :: Nx
      real(wp), intent(in)  :: ys(6*Nx)
      real(wp), intent(out) :: density(Nx), velocity(Nx), temperature(Nx), neutral(Nx), neutral_velocity(Nx), molecule(Nx)
      !real(wp), intent(out) :: extern_neutral_density(5), extern_molecule_density(5)
      ! density
      density(1:Nx)     =  exp(ys(1:Nx)) * density_norm
      ! velocity
      velocity(1:Nx)    =  ys(  Nx+1: 2*Nx) * momentum_norm / mass / density(1:Nx)
      ! temperature
      temperature(1:Nx) =  ys(2*Nx+1: 3*Nx) * energy_norm / 3.0d+0 / density(1:Nx) / e_charge
      ! neutral density
      neutral(1:Nx)     =  exp(ys(3*Nx+1: 4*Nx)) * density_norm
      ! neutral velocity
      neutral_velocity(1:Nx)    =  ys(4*Nx+1: 5*Nx) * momentum_norm / mass / neutral(1:Nx)
      ! molecule density
      molecule(1:Nx)     =  exp(ys(5*Nx+1: 6*Nx)) * density_norm
     
	tmp = 0.1d+0
	do i = 1, Nx
	tmp1 = min(temperature(i), 0.1d+0)
	if(tmp1 .lt. tmp) then
	tmp = tmp1
	endif	
	enddo
	
	if ( tmp .lt. 0.1d+0 ) then
	write(*,*)'density(2)', density(2)
	!write(*,*)'velocity', velocity
	!write(*,*)'neutral', neutral
	write(*,*) 'temperature(2)', temperature(2)
	!write(*,*)
	temperature = max(temperature, minimum_temperature) ! clamp to 0.1 eV.
	endif	
 
	! removes NaN as well as constraining solutions
	density = min(max(density,minimum_density), maximum_density )
	neutral = min(max(neutral,minimum_density), maximum_density )
	molecule = min(max(molecule,minimum_density), maximum_density)	

      return
   end subroutine ys2nvt

   subroutine nr2yr(extern_neutral_density, extern_molecule_density, yr)
	implicit none
	real(wp), intent(in) :: extern_neutral_density(5), extern_molecule_density(5)
	real(wp), intent(out) :: yr(10)
        yr(1:5) = log(extern_neutral_density(1:5) / density_norm)
      	yr(6:10) = log(extern_molecule_density(1:5) / density_norm)
	return
   end subroutine nr2yr

  subroutine yr2nr(yr, extern_neutral_density, extern_molecule_density)
	implicit none
	real(wp), intent(in) :: yr(10)
	real(wp), intent(out) :: extern_neutral_density(5), extern_molecule_density(5)
	extern_neutral_density(1:5) = exp(yr(1:5)) * density_norm
        extern_molecule_density(1:5) = exp(yr(6:10)) * density_norm
	! make sure that the guesses for density are not below minima or above maxima ! this could also be specified in solver settings
	extern_neutral_density = min(max(extern_neutral_density,minimum_density), maximum_density )
	extern_molecule_density = min(max(extern_molecule_density,minimum_density), maximum_density)	
	return
   end subroutine yr2nr


   subroutine yc2nc(yc,core_density)
	! y_core to core_density
	implicit none
	real(wp), intent(in) :: yc
	real(wp), intent(out) :: core_density
	core_density     =  exp(yc) * density_norm
	core_density = min(max(core_density, minimum_density),maximum_density)
   end subroutine yc2nc

   subroutine nc2yc(core_density, yc)
	! core density to y_core solution
	implicit none
	real(wp), intent(in) :: core_density
	real(wp), intent(out) :: yc
	yc = log(core_density / density_norm)
   end subroutine nc2yc

   subroutine advection3(Nx, variable, velocity_cb, csound_cb, delta_x, delta_xcb, flux)
   ! this subroutine calculates the advected flux of a variable with a given velocity field
   ! we follow here the discretization as put forward in B. Dudson et al. (2019) PPCF 61 065008  
   ! the temperature is needed to calculate the sound velocity on cell boundaries
   ! the variable and velocity field are given on the (non-equidistant) grid centers
   ! the flux is returned at the grid boundaries i + 1/2 (internal only: i=1,...,Nx-1, BC's must be applied outside this routine)
      implicit none
      integer,  intent(in)  :: Nx
      real(wp), intent(in)  :: variable(Nx), velocity_cb(0:Nx), csound_cb(0:Nx)
      real(wp), intent(in)  :: delta_x(1:Nx-1), delta_xcb(1:Nx)
      real(wp), intent(out) :: flux(0:Nx)
      real(wp)		    :: variable_cb(0:Nx)
      real(wp) 	            :: slopes(0:Nx) ! slopes on cell boundaries
      real(wp)		    :: limslope(1:Nx) ! limited slope for cell centers
      real(wp) 	            :: varL(1:Nx) ! limslope extrapolated to left boundary value of the cell
      real(wp) 		    :: varR(1:Nx) ! limslope extrapolated to right boundary value of cell
      integer               :: i
      flux = 0.0d+0
  
      ! slopes are calculated betweeen cell centers:
      slopes(1:Nx-1) = (variable(2:Nx) - variable(1:Nx-1)) / delta_x(1:Nx-1)
      ! calculate slopes with the extrapolated values for BC
      ! for linear extrapolation this is the same, but for quadratic it is not
      !slopes(0) = (variable(1) - variable_cb(0)) / (delta_xcb(1) / 2.0d+0)
      !slopes(Nx) = (variable_cb(Nx) - variable(Nx)) / (delta_xcb(Nx) / 2.0d+0)
      ! the alternative is to take slope = 0 for the last cells or to simply take the previous value
      slopes(0) = slopes(1) !0.0d+0 !(variable(1) - variable_cb(0)) / (delta_xcb(1) / 2.0d+0)
      slopes(Nx) = slopes(Nx-1) !0.0d+0 !(variable_cb(Nx) - variable(Nx)) / (delta_xcb(Nx) / 2.0d+0)

      ! take the minmod limited slopes per cell value
      do i = 1,Nx 
      limslope(i) = minmod2(slopes(i-1), slopes(i))
      enddo  
      !write(*,*), 'minmod2', minmod2(-1.0,-2.0), minmod2(1.0,2.0), minmod2(-1.0,1.0)

      ! use limited slopes to calculate variable values on Left and Right cell boundary
      varL(1:Nx) = variable(1:Nx) - limslope(1:Nx) * delta_xcb(1:Nx) / 2.0d+0
      varR(1:Nx) = variable(1:Nx) + limslope(1:Nx) * delta_xcb(1:Nx) / 2.0d+0
      !write(*,*), 'varR, cb, varLp', varR(Nx-1), variable_cb(Nx-1), varL(Nx)

      do i = 1, Nx-1      
         if( velocity_cb(i) .gt. csound_cb(i) ) then 
	 ! advection is dominated by the right boundary of the left cell
		flux(i) = varR(i)*velocity_cb(i)

         elseif( velocity_cb(i) .lt. -csound_cb(i) ) then
	 ! advection is dominated by the left boundary of the right cell
		flux(i) = varL(i+1)*velocity_cb(i)
	 else
      	 ! calculate the balance and correct with Lax flux
		flux(i) = (varR(i)+varL(i+1))/2.0d+0 * velocity_cb(i) &  ! physical flux
	        + lax_switch * (varR(i)-varL(i+1))/2.0d+0 * csound_cb(i) ! Lax flux 
         endif
      enddo

   end subroutine advection3
	
   subroutine calculate_sol2extern_ion_flux(sol2extern_ion_flux, sum_sol2extern_ion_flux, density, gamma_core2sol)
	! calculates possible far sol ion losses in the core-sol to have main chamber recycling
	! uses  physics params, only : core_far_sol_ion_loss, core_far_sol_feedthrough 
	! uses  grid data, only : Area extern + core_source_profile_n
	implicit none
	integer i
	real(wp),intent(in) :: density(Nx)
	! todo: make sure to update gamma_core2sol here as well
	real(wp),intent(in) :: gamma_core2sol
	real(wp),intent(out):: sol2extern_ion_flux(Nx)
	real(wp),intent(out) :: sum_sol2extern_ion_flux
	sum_sol2extern_ion_flux = 0.0d+0
	sol2extern_ion_flux = 0.0d+0
	do i = 1,Nx
	if( i .ge. i_baffle(1) .and. i .le. i_baffle(2)) then
		sum_sol2extern_ion_flux = sum_sol2extern_ion_flux + Area_extern(i) * density(i) * core_far_sol_ion_loss 
	endif
	enddo
	! rescale solution to make it smooth following the core-influx distribution
	sol2extern_ion_flux = core_source_profile_n * sum_sol2extern_ion_flux	
	! add the particles that directly come from the core (note that this sum and distribution are asymmetric! 
	sum_sol2extern_ion_flux = sum_sol2extern_ion_flux + gamma_core2sol*core_far_sol_feedthrough 
	! gamma_core2sol = gamma_core2sol*(1.0d0 - core_far_sol_feedthrough) 
	if(core_far_sol_feedthrough .gt. 0.0d0 ) then
	write(*,*) 'warning far sol feedthrough not tested'
	endif
   end subroutine calculate_sol2extern_ion_flux

   subroutine calculate_extern_sol_core_shinethrough_neutral_fluxes(extern2core_st,extern2sol_flux,extern_neutral, sol_neutral, residence_time, core_ion_frac)
	implicit none
	real(wp), intent(out) :: extern2core_st(Nx) !, sum_extern2core_st
 	real(wp), intent(inout):: extern2sol_flux(Nx)
	real(wp), intent(in) :: extern_neutral(5), sol_neutral(Nx), residence_time, core_ion_frac
	! internal vars
	real(wp) :: dummy5(5)
	real(wp) :: dummy3(3), dummy33(3)
	
	call calculate_sol_extern_neutral_fluxes(extern2sol_flux, dummy5, & 
					dummy3, &
					extern_neutral, sol_neutral, residence_time, dummy33 )
	extern2core_st = core_ion_frac*extern2sol_flux
	extern2sol_flux = extern2sol_flux - extern2core_st	
   end subroutine calculate_extern_sol_core_shinethrough_neutral_fluxes

   subroutine calculate_sol_extern_neutral_fluxes(extern2sol_flux, sol2extern_flux, & 
					extern_fluxes, &
					extern_neutral, sol_neutral_ext, residence_time, extern_ex )
		! calculates neutral fluxes between sol and external reservoirs
		!uses grid_data, only: Nx, i_Xpoint, i_baffle, sigma_nb, volumes, delta_xcb, Area_extern		
		implicit none
		integer :: i
		real(wp) ::  	tmp_flux1 = 0.0d+0
         	real(wp) ::  	tmp_flux2 = 0.0d+0
		real(wp) ::  	tmp_flux3 = 0.0d+0
         	real(wp) ::  	tmp_flux4 = 0.0d+0
         	real(wp) ::  	tmp_flux5 = 0.0d+0
	 	real(wp) :: 	red_target_flux(2) = 0.0d+0
		real(wp),intent(in)  :: extern_neutral(5), sol_neutral_ext(Nx)
		real(wp) :: sol_neutral(Nx) 
		real(wp),intent(out) :: extern_fluxes(3)
		real(wp),intent(in)  ::	residence_time, extern_ex(3)
		real(wp),intent(out) :: extern2sol_flux(Nx) !, sol2extern_ion_flux(Nx)
		real(wp),intent(out) :: sol2extern_flux(5)
		sol_neutral = sol_neutral_ext* ato_res_asy  !  apply assymmetric partcle flows?
		!real(wp),intent(out) :: extern2core_flux(5)
		extern2sol_flux = 0.0d+0
		extern_fluxes = 0.0d+0
		!sol2extern_ion_flux = 0.0d+0 ! this can be a source of atoms in the reservoirs or molecules according to the mol_rec fraction
		sol2extern_flux =  (/0.0d+0,0.0d+0,0.0d+0,0.0d+0,0.0d+0/)
	
 	 do i = 1, Nx             
         if( i .lt. i_baffle(1) ) then !(x(i) .lt. X_core_SOL)  then ! first divertor leg bound to PFR and D/CFR
		 tmp_flux1 = Area_extern(i) * ( -sol_neutral(i) + extern_neutral(1)) / residence_time
		 tmp_flux2 = Area_extern(i) * ( -sol_neutral(i) + extern_neutral(2)  &
			- (extern_neutral(2) - extern_neutral(3)) * sigmoid(x(i),x(i_baffle(1)),sigma_nb) ) &
			/ neutral_residence_time ! Added sigmoid function
		 extern2sol_flux(i) = tmp_flux1 + tmp_flux2                  ! flux into sol
		 sol2extern_flux(1) = sol2extern_flux(1) - tmp_flux1 ! MINUS        ! flux out of external volume                        
		 sol2extern_flux(2) = sol2extern_flux(2) - tmp_flux2 ! MINUS
	 elseif( i .ge. i_baffle(1) .and. i .lt. i_Xpoint(1) ) then ! Between baffle and xpoint: PFR and CFR
		 if (wide_PFR == 1) then 
		 	tmp_flux1 = Area_extern(i) * ( -sol_neutral(i) + extern_neutral(1) ) / residence_time
		 else
		 	tmp_flux1 = 0.0d+0
		 endif
		 tmp_flux3 = Area_extern(i) * ( -sol_neutral(i) + extern_neutral(2) & 
			- (extern_neutral(2) - extern_neutral(3)) * sigmoid(x(i),x(i_baffle(1)),sigma_nb) ) &
			/ residence_time
		 extern2sol_flux(i) =  tmp_flux1 + tmp_flux3               ! flux into sol
	
		 sol2extern_flux(1) = sol2extern_flux(1) - tmp_flux1 ! MINUS         ! flux out of external volume                        
		 sol2extern_flux(3) = sol2extern_flux(3) - tmp_flux3 ! MINUS
	         ! wait here we are counting double, there is no flux to the core between baffle and Xpoint	
		 !core2sol_flux(i)   =  -core_ionization_fraction*tmp_flux3 ! leave core ionization fraction the same for now???
		 !sol2core_flux = sol2core_flux + core2sol_flux(i)  
	 elseif( i .ge. i_Xpoint(1) .and. i .le. i_Xpoint(2)) then ! between xpoints: only CFR
		 extern2sol_flux(i) = Area_extern(i) * ( -sol_neutral(i) + extern_neutral(3)) / residence_time
	 	 !core2sol_flux(i) = -core_ionization_fraction*extern2sol_flux(i) ! here you can also use the core_pumping? 
		 sol2extern_flux(3) = sol2extern_flux(3) - extern2sol_flux(i) ! MINUS
		 !sol2core_flux = sol2core_flux + core2sol_flux(i)  
         elseif( i .gt. i_Xpoint(2) .and. i .le. i_baffle(2) ) then  !between xpoint and baffle: PFR and CFR
		 tmp_flux3 = Area_extern(i) * ( -sol_neutral(i) + extern_neutral(3) &
			+ (extern_neutral(4) - extern_neutral(3)) * sigmoid(x(i),x(i_baffle(2)),sigma_nb) )&
			/ residence_time
		 if (wide_PFR == 1) then
		 	tmp_flux5 = Area_extern(i) * ( -sol_neutral(i) + extern_neutral(5) ) / residence_time
		 else
		 	tmp_flux5 = 0.0d+0
		 endif 
		 extern2sol_flux(i) = tmp_flux3 + tmp_flux5                  ! flux into sol
		 sol2extern_flux(3) = sol2extern_flux(3) - tmp_flux3 ! MINUS        ! flux out of external volume                        
		 sol2extern_flux(5) = sol2extern_flux(5) - tmp_flux5 ! MINUS
         else ! second divertor leg bound by PFR and D/CFR 
		 tmp_flux4 = Area_extern(i) * (-sol_neutral(i) + extern_neutral(3) & 
			+ (extern_neutral(4) - extern_neutral(3)) * sigmoid(x(i),x(i_baffle(2)),sigma_nb) ) & 
			/ residence_time
		 tmp_flux5 = Area_extern(i) * (-sol_neutral(i) + extern_neutral(5) ) / residence_time
		 extern2sol_flux(i) = tmp_flux4 + tmp_flux5              ! flux into the SOL       
		 sol2extern_flux(4) = sol2extern_flux(4) - tmp_flux4     ! flux out of external volume                             
		 sol2extern_flux(5) = sol2extern_flux(5) - tmp_flux5     ! should there be a sign different? 
         endif
         enddo
	 ! fluxes between reservoirs
  	 ! volumes go as: (1)left-PFR, (2)left-DFR, (3)CFR, (4)right-DFR, (5)right-PFR 
	 call calculate_extern_fluxes(extern_fluxes,extern_neutral,extern_ex)
     	 !extern_fluxes =          (/(extern_neutral(5)-extern_neutral(1))*extern_ex(1), &   ! from 5->1
          !                        (extern_neutral(2)-extern_neutral(3))*extern_ex(2), &	    ! from 2->3
         !                         (extern_neutral(3)-extern_neutral(4))*extern_ex(3) /)     ! from 3->4
	 return
   end subroutine calculate_sol_extern_neutral_fluxes

   subroutine calculate_extern_fluxes(extern_fluxes,extern_neutral,extern_ex)
	 	implicit none
		real(wp), intent(in) :: extern_neutral(5), extern_ex(3)
		real(wp), intent(out):: extern_fluxes(3) 
 		! fluxes between reservoirs
  		! volumes go as: (1)left-PFR, (2)left-DFR, (3)CFR, (4)right-DFR, (5)right-PFR 
		! fluxes go as : (1) from 5to1, (2) from 2to3, (3) from 3to4
     	 extern_fluxes =          (/(extern_neutral(5)-extern_neutral(1))*extern_ex(1), &   ! from 5->1
                                  (extern_neutral(2)-extern_neutral(3))*extern_ex(2), &	    ! from 2->3
                                  (extern_neutral(3)-extern_neutral(4))*extern_ex(3) /)     ! from 3->4
   end subroutine calculate_extern_fluxes

   subroutine calculate_sol2core_neutral_flux(core2sol_neutral_flux, sol2core_neutral_flux, core_neutral_density , sol_neutral_density, core_sol_neutral_ex_, neutral_core_ionization_flux )
	       	! calculates the neutral fluxes between core and sol (for both atoms and molecules if you will)	
		! uses grid_data, only: Nx, i_Xpoint, i_baffle, sigma_nb, volumes, delta_xcb, Area_extern		
		implicit none
		integer :: i
		real(wp) ::  	tmp_flux1 = 0.0d+0
         	real(wp) ::  	tmp_flux2 = 0.0d+0
		real(wp) ::  	tmp_flux3 = 0.0d+0
         	real(wp) ::  	tmp_flux4 = 0.0d+0
         	real(wp) ::  	tmp_flux5 = 0.0d+0
	 	real(wp) :: 	red_target_flux(2) = 0.0d+0
		real(wp),intent(in)  :: core_neutral_density, sol_neutral_density(Nx)
		real(wp),intent(in)  ::	core_sol_neutral_ex_
		real(wp),intent(in)  :: neutral_core_ionization_flux(Nx)
		real(wp),intent(out) :: core2sol_neutral_flux(Nx)
		real(wp),intent(out) :: sol2core_neutral_flux
	 sol2core_neutral_flux = 0.0d+0
	 core2sol_neutral_flux = 0.0d+0
 	 do i = 1, Nx
	 if( i .ge. i_Xpoint(1) .and. i .le. i_Xpoint(2)) then           
		 core2sol_neutral_flux(i) = Area_extern(i) * ( -sol_neutral_density(i) + core_neutral_density ) * core_sol_neutral_ex_
		 core2sol_neutral_flux(i) = core2sol_neutral_flux(i) - neutral_core_ionization_flux(i)
		 sol2core_neutral_flux = sol2core_neutral_flux - core2sol_neutral_flux(i) ! MINUS
		 ! note if core_neutral_density = 0.0, 
		 ! core_sol_neutral_ex_time = residence time * (neutral_sol-ncore) / (fcore_ion *(neutral sol - neutral reservoir))
		 ! will reproduce the core ionization fraction used in Derks2024
		 !sol2core_flux = sol2core_flux + core2sol_flux(i)  
		 ! setting the core_ionization_fraction will add a neutral flux as function of the extern<->sol fluxes.
	 endif
         enddo
	 return
   end subroutine calculate_sol2core_neutral_flux

   subroutine calculate_core_ionization_neutral_flux(extern2core_shinethr, extern2sol_fl, core_ionization_frac )
	! calculate the flux of neutral (atoms) that passes the SOL into the core
	! core ionization fraction also makes sense in density ramp:
	! if SOL hotter, low den ion populatio -> more ionization, more exter2sol influx -> more into core as well? 
	! if SOL colder, denser ion population -> less extern2sol influx -> less into core? 
	implicit none
	real(wp), intent(in) :: extern2sol_fl(Nx)
        real(wp), intent(in) :: core_ionization_frac ! use from settings? 
	real(wp), intent(out):: extern2core_shinethr(Nx)
	extern2core_shinethr = extern2sol_fl*core_ionization_frac
	return
   end subroutine calculate_core_ionization_neutral_flux

   subroutine calculate_wall_association(atom_sink, molecule_source, wall_area, atom_density, inverse_atom_exchange_speed, association_probability)
	! calculate the number of atoms that hit the wall in the reservoirs.
	! of the atoms that hit the wall, a fraction associates into molecules. 
 	! the sink of atoms that associate is calculated as
	! S_atoms = A_wall*1/4*v_atom*atom_density = A_wall*atom_density*atom_exchange_speed 
	implicit none
	real(wp), intent(in) :: atom_density(5), wall_area(5), inverse_atom_exchange_speed, association_probability
	real(wp), intent(out):: atom_sink(5), molecule_source(5) ! in units D/s and D2/s 
	atom_sink = -atom_density*wall_area*association_probability/inverse_atom_exchange_speed
	molecule_source = -atom_sink/2.0d0 ! factor two for conversion.
	! write(*,*) 'residual =', atom_sink -molecule_source ! should be zero as unit test
   end subroutine calculate_wall_association

   subroutine calculate_extern_pump(extern_pump,extern_neutral,pump_rate)
	implicit none
	real(wp), intent(in) :: extern_neutral(5)
	!real(wp), intent(in) :: volumes(5)
	real(wp), intent(out) :: extern_pump(5)
	real(wp), intent(in) :: pump_rate(5)
	extern_pump = extern_neutral*pump_rate*extern_neutral_volumes
	return
   end subroutine calculate_extern_pump

   subroutine calculate_extern2core_fluxes(extern2core,extern_neut_density,pump_extern2core)
	implicit none
	real(wp), intent(in) :: extern_neut_density(5)
	real(wp), intent(out) :: extern2core(5)
	real(wp), intent(in) :: pump_extern2core(5)
	!write(*,*) 'extern2core_fluxes: pump_core', pump_core
	 extern2core = 0.0d+0
	extern2core = pump_extern2core * extern_neut_density * extern_neutral_volumes
	return
   end subroutine calculate_extern2core_fluxes

   subroutine calculate_sol_recycle_fluxes(Gamma_n0, Gamma_nL, rec_frac, tar2extern_flux) 
	implicit none
	real(wp), intent(in)  :: Gamma_n0, Gamma_nL, rec_frac
	real(wp), intent(out) :: tar2extern_flux(5) ! = 0.0d+0
	real(wp) :: red_target_flux(2) = 0.0d+0
	tar2extern_flux(5)  = 0.0d+0
 	 ! Target fluxes can reflect in different ways
         ! volumes go as: (1)left-PFR, (2)left-DFR, (3)CFR, (4)right-DFR, (5)right-PFR 
	 if( X_core_SOL .eq. 0.0d+0 ) then
	 red_target_flux(1) = 0.0d+0
	 tar2extern_flux(1) = 0.0d+0
	 tar2extern_flux(2) = 0.0d+0
	 else	
	 red_target_flux(1) = - A_wet(1) * B_field_cb(0)  * Gamma_n0  * (1.0d+0 - rec_frac)
	 tar2extern_flux(1) = red_target_flux(1) *           min(max(pol_target_angle(1),0.0d+0),180d+0) / 180d+0
	 tar2extern_flux(2) = red_target_flux(1) * (1.0d+0 - min(max(pol_target_angle(1),0.0d+0),180d+0) / 180d+0 )
	 endif
	 ! always have second target
 	 tar2extern_flux(3) = 0.0d+0! we could have main chamber recycling !
	 red_target_flux(2) = A_wet(2) * B_field_cb(Nx) * Gamma_nL * (1.0d+0 - rec_frac) 
         tar2extern_flux(4) = red_target_flux(2) *           min(max(pol_target_angle(2),0.0d+0),180d+0) / 180d+0
	 tar2extern_flux(5) = red_target_flux(2) * (1.0d+0 - min(max(pol_target_angle(2),0.0d+0),180d+0) / 180d+0 )
	return
   end subroutine calculate_sol_recycle_fluxes
	
   subroutine calculate_core2sol_ion_flux( core_density, Gamma_core2sol)
	! this subroutine calculates the core outflux, only Gamma_core2sol for now
	implicit none
	real(wp), intent(in) :: core_density
	real(wp), intent(out) ::  Gamma_core2sol
	Gamma_core2sol = core_density * core_volume / core_confinement_time 		
	return
   end subroutine calculate_core2sol_ion_flux

   subroutine calculate_core_source( extern2core_flux, extern2core_mol, Gamma_core2sol, sol2core_flux, sol2core_mol, Source_core )
	! this subroutine calculates the core source
	implicit none
	real(wp), intent(in) :: Gamma_core2sol,  extern2core_flux(5), extern2core_mol(5)
	real(wp), intent(in) :: sol2core_flux, sol2core_mol
	real(wp), intent(out) :: Source_core	
	Source_core = (D_core_fuelling(internal_istep_rhs)  &
			+ sum(extern2core_flux) + sol2core_flux  &
	        + 2.0d+0*(sum(extern2core_mol) + sol2core_mol)&
			- Gamma_core2sol ) / core_volume 
	return
   end subroutine calculate_core_source

   subroutine calculate_sol_fluxes( Nx, density, velocity, temperature, neutral, neutral_velocity, molecule, &
				extern_neutral_density,  extern_molecule_density,  &
				Gamma_n, Gamma_mom, q_parallel, Gamma_neutral, Gamma_molecule, Gamma_mom_neutral,& 
				elm_heat_load, & 
                                extern_neutral_flux, extern_molecule_flux, &
				sol2extern_flux, sol2extern_mol, & 
				extern2sol_flux, extern2sol_mol, tar2extern_flux, tar2extern_mol)
   ! this subroutine calculates the 'fluxes' required for the calculation of ydot (the right hand side of the discretized conservation equations)
   ! the fluxes are defined at the halfway point between grid points: i.e. Flux(i) is midway between x(i) and x(i+1) on the cell boundaries xcb
   ! ew 01-03-2021: modified to take account of flux expansion
      implicit none
      integer,  intent(in)  :: Nx
      real(wp), intent(in)  :: density(Nx), velocity(Nx), temperature(Nx), neutral(Nx), neutral_velocity(Nx), molecule(Nx) ! , core_density
      real(wp)		    :: density_cb(0:Nx), velocity_cb(0:Nx), temperature_cb(0:Nx),  csound_cb(0:Nx) !neutral_cb(0:Nx),
      real(wp) 		    :: neutral_cb(0:Nx), neutral_velocity_cb(0:Nx), molecule_cb(0:Nx) 
      real(wp), intent(out) :: Gamma_n(0:Nx), Gamma_mom(0:Nx), q_parallel(0:Nx), Gamma_neutral(0:Nx), Gamma_mom_neutral(0:Nx), Gamma_molecule(0:Nx)
      real(wp), intent(in)  :: extern_neutral_density(5), extern_molecule_density(5)
      real(wp), intent(out) :: extern_neutral_flux(3), sol2extern_flux(5), extern2sol_flux(Nx), tar2extern_flux(5) 
      real(wp), intent(out) :: extern_molecule_flux(3), sol2extern_mol(5), extern2sol_mol(Nx), tar2extern_mol(5) 
      real(wp), intent(in)  :: elm_heat_load
      real(wp)              :: momentum(Nx), enthalpy(Nx), densityB(Nx), neutral_momentum(Nx)
      real(wp)   	    :: average_velocity, vrel_cb
      real(wp) 		    :: red_target_flux(2) !, A_wet(2)
      real(wp) 		    :: leakage_fluxes_atoms(5), leakage_fluxes_mol(5)
      real(wp)              :: tmp_flux1,tmp_flux2,tmp_flux4,tmp_flux5 ! dummy parameters to avoid duplicate calculations
      integer               :: i, ix
	if( detect_nan ) then
      	do ix = 1,Nx
 	if (isnan(density(ix))) stop 'calc_flux: "density" is a NaN'
	if (isnan(velocity(ix))) stop 'calc_flux: "velocity" is a NaN'
	if (isnan(temperature(ix))) stop 'calc_flux: "temperature" is a NaN'
	if (isnan(neutral(ix))) stop 'calc_flux: "neutral" is a NaN'
	if (isnan(neutral_velocity(ix))) stop 'calc_flux: "neutral_velocity" is a NaN'
	if (isnan(molecule(ix))) stop 'calc_flux: "molecule" is a NaN'
      	enddo
	endif
      ! calculate values on cell boundaries
      call cc2cb(Nx,density	,density_cb)   
      call cc2cb(Nx,velocity	,velocity_cb)
      call cc2cb(Nx,temperature	,temperature_cb)
      !call cc2cb(Nx,neutral 	,neutral_cb)
      call cc2cb(Nx,neutral_velocity, neutral_velocity_cb)
      !call cc2cb(Nx,neutral 	,neutral_cb)
      !call cc2cb(Nx,molecule 	,molecule_cb)
      csound_cb = sqrt( 2.0d+0 * e_charge * max(temperature_cb,minimum_temperature) / mass )   

      ! the particle flux = density velocity / B 
	 densityB = density / B_field
         call advection3(Nx, densityB, velocity_cb, csound_cb, delta_x, delta_xcb, Gamma_n)
	 Gamma_n(Nx) = max(minimum_density,density_cb(Nx)) * max(velocity_cb(Nx),csound_cb(Nx)) / B_field_cb(Nx)
         ! boundary condition at i = 0 (not used when X-point density BC is used)
         if( X_core_SOL .eq. 0.0d+0 ) then
             ! case for core-SOL + divertorleg (i=0 is stagnation point)
             Gamma_n(0)  = 0.0d+0
         else
             ! case with two targets so sheath boundary conditions at i=0 as well
	     Gamma_n(0) = max(minimum_density,density_cb(0)) * min(velocity_cb(0),-csound_cb(0)) / B_field_cb(0)
         endif

      ! the momentum flux = momentum * velocity where momentum = density * mass * velocity / B_field
         momentum = density * mass * velocity / B_field
         call advection3(Nx, momentum, velocity_cb, csound_cb, delta_x, delta_xcb, Gamma_mom)
         !call advection(Nx, momentum, velocity, temperature, Gamma_mom)
         ! boundary condition at the sheath
	Gamma_mom(Nx) = max(minimum_density,density_cb(Nx)) * mass * max(velocity_cb(Nx),csound_cb(Nx))**2 / B_field_cb(Nx)
         ! boundary condition at i = 0
         if( X_core_SOL .eq. 0.0d+0 ) then
             ! case for core-SOL + divertorleg (i=0 is stagnation point)
             Gamma_mom(0)  = 0.0d+0
         else
             ! case with two targets so sheath boundary conditions at i=0 as well  
	     Gamma_mom(0) =  max(minimum_density,density_cb(0)) * mass * min(velocity_cb(0),-csound_cb(0))**2 / B_field_cb(0)
         endif
 	
      ! convective heat flux = 5 density k temperature velocity / B_field (i.e. 5/2 pressure)
         enthalpy = 5.0d+0 * density * e_charge * temperature / B_field
         call advection3(Nx, enthalpy, velocity_cb, csound_cb, delta_x, delta_xcb, q_parallel) 
         !write(*,*) 'advection', q_parallel

         ! add the conductive heat flux in the internal region
         do i = 1, Nx-1
	   q_parallel(i) = q_parallel(i) * switch_convective_heat - kappa_parallel(temperature_cb(i)) * (temperature(i+1)-temperature(i))/delta_x(i) / B_field_cb(i)
         enddo

        ! boundary condition at the sheath: sheath heat transmission (note gamma is defined w. csound at the target)
	q_parallel(Nx) = gamma * max(velocity_cb(Nx),csound_cb(Nx)) * max(minimum_density,density_cb(Nx)) * e_charge * max(temperature_cb(Nx),minimum_temperature) / B_field_cb(Nx)   
         ! boundary condition at i = 0
         if( L_core_SOL .eq. 0.0 ) then
             ! case with prescribed heat flux at i = 0 (X-point)
             q_parallel(0)  = D_qpar_x(internal_istep_rhs)+elm_heat_load
         elseif( X_core_SOL .eq. 0.0 ) then
             ! case for core-SOL + divertorleg (i=0 is stagnation point)
             q_parallel(0)  = 0.0d+0
         else
             ! case with two targets so sheath boundary conditions at i=0 as well
  	     q_parallel(0) = gamma * min(velocity_cb(0),-csound_cb(0)) * max(minimum_density,density_cb(0)) * e_charge * max(temperature_cb(0),minimum_temperature) / B_field_cb(0) 
         endif
         !write(*,*) 'temperature =', temperature
         !write(*,*) 'q_parallel =', q_parallel
         
	! Calculate neutral flux
	 call advection3(Nx, neutral, neutral_velocity_cb, csound_cb, delta_x, delta_xcb, Gamma_neutral) 

	 ! neutral particle diffusion
         do i = 1, Nx-1
  	     vrel_cb = velocity_cb(i)-neutral_velocity_cb(i)
	     Gamma_neutral(i) = Gamma_neutral(i) & 
		- D_neutral(temperature_cb(i),density_cb(i),vrel_cb) * (neutral(i+1)-neutral(i))/delta_x(i)
         enddo
	!write(*,*) 'F internal_istep_rhs = ', internal_istep_rhs
	!write(*,*) 'F internal_istep = ', internal_istep
         ! boundary condition at i = 0
         if( L_core_SOL .eq. 0.0d+0 ) then
             ! boundary condition at X-point (zero gradient i.e. at i=0 every equals i=1, i.e. zero flux)
             Gamma_neutral(0) = 0.0d+0
         elseif( X_core_SOL .eq. 0.0d+0 ) then
             ! apply boundary condition at mid point (i.e. flux = 0)
             Gamma_neutral(0) = 0.0d+0
         else
             ! sheath boundary condition (= flux of plasma density in case of full recycling)
	     Gamma_neutral(0) = -B_field_cb(0)*Gamma_n(0) * D_rec(internal_istep_rhs)* (1.0d0 - mol_rec)  
         endif
         ! boundary condition at the sheath (= flux of plasma density in case of full recycling)
         Gamma_neutral(Nx) = -B_field_cb(Nx)*Gamma_n(Nx) * D_rec(internal_istep_rhs)* (1.0d0 - mol_rec) 
	 !write(*,*) 'calc_flux: Gamma_neutral', Gamma_neutral

      ! Neutral momentum flux: 
	neutral_momentum = neutral * mass * neutral_velocity
	call advection3(Nx, neutral_momentum, neutral_velocity_cb, csound_cb, delta_x, delta_xcb, Gamma_mom_neutral)
         ! at target reflection of neutrals leads to Gamma_mom_neutral(Nx) = 2*(m Nn Vn)*Vn)
	!Gamma_mom_neutral(Nx) = 2*max(minimum_density,neutral_cb(Nx)) * mass * max(neutral_velocity_cb(Nx),csound_cb(Nx))**2 
	! see section 3.1 in Stangeby 2000 book: Figure 3.1 atom energy reflection coefficient, note the root because energy is converted to velocity here. E = 0.5*mv^2
	Gamma_mom_neutral(Nx) = Gamma_neutral(Nx) * mass * -1.0d+0 * atom_recycle_energy_fraction**0.5d+0 *csound_cb(Nx) 
	 !write(*,*) 'FLUX: Gamma_mom(Nx-10)', Gamma_mom(Nx-10)
         ! boundary condition at i = 0
         if( X_core_SOL .eq. 0.0d+0 ) then
             ! case for core-SOL + divertorleg (i=0 is stagnation point)
             Gamma_mom_neutral(0)  = 0.0d+0
         else
             ! case with two targets so sheath boundary conditions at i=0 as well   ! why would this be limited to cs?
	     Gamma_mom_neutral(0) = Gamma_neutral(0) * mass * 1.0d+0 * atom_recycle_energy_fraction**0.5d+0 *csound_cb(0) !  2*max(minimum_density,neutral_cb(0)) * mass * min(neutral_velocity_cb(0),-csound_cb(0))**2
         endif

	! Molecular diffusion:
	 do i = 1, Nx-1
	     Gamma_molecule(i) = -D_molecule(temperature_cb(i),density_cb(i), velocity_cb(i)) * (molecule(i+1)-molecule(i))/delta_x(i)
         enddo
         ! boundary condition at i = 0
         if( L_core_SOL .eq. 0.0d+0 ) then
             ! boundary condition at X-point (zero gradient i.e. at i=0 every equals i=1, i.e. zero flux)
             Gamma_molecule(0) = 0.0d+0
         elseif( X_core_SOL .eq. 0.0d+0 ) then
             ! apply boundary condition at mid point (i.e. flux = 0)
             Gamma_molecule(0) = 0.0d+0
         else
             ! sheath boundary condition (= reassociation flux from target could add Eley-Rideal model for reassociation)
             Gamma_molecule(0) = -B_field_cb(0)*Gamma_n(0) * D_rec(internal_istep_rhs) * 0.5d+0 * mol_rec
         endif
         ! sheath boundary condition (= reassociation flux from target could add Eley-Rideal model for reassociation)
         Gamma_molecule(Nx) = -B_field_cb(Nx)*Gamma_n(Nx) * D_rec(internal_istep_rhs) * 0.5d+0 * mol_rec  
	
  	
	 ! Recycling fluxes from sol to external reservoirs
	 call calculate_sol_recycle_fluxes(Gamma_n(0), Gamma_n(Nx), D_rec(internal_istep_rhs), tar2extern_flux) 
	 tar2extern_mol = tar2extern_flux*mol_rec * 0.5d+0 ! 2 ions make 1 molecule
	 tar2extern_flux = tar2extern_flux*(1.0d+0 - mol_rec)
	 !write(*,*) 'calc sol flux: density',extern_neutral_density 
         !write(*,*) 'calc sol flux: exrate',extern_neutral_ex 

	 ! Exchange fluxes between sol and external reservoirs, and the fluxes between the reservoirs themselves
	 call calculate_sol_extern_neutral_fluxes(extern2sol_flux, sol2extern_flux, extern_neutral_flux, & !Out
					extern_neutral_density, neutral, neutral_residence_time, extern_neutral_ex )!, & ! In
					!	core_ionization_fraction ) ! In
	 call calculate_sol_extern_neutral_fluxes(extern2sol_mol, sol2extern_mol, extern_molecule_flux, & ! Out
					extern_molecule_density, molecule, molecule_residence_time, extern_molecule_ex) !, & ! In
					!	core_ionization_fraction ) ! In
	! add leakage from B field as source to sol2extern flows
	call calculate_leakage_fluxes(Gamma_neutral, leakage_fluxes_atoms)
	call calculate_leakage_fluxes(Gamma_molecule, leakage_fluxes_mol)
	sol2extern_flux = sol2extern_flux + leakage_fluxes_atoms * switch_neutral_leakage
	sol2extern_mol = sol2extern_mol + leakage_fluxes_mol * switch_neutral_leakage
         return
   end subroutine calculate_sol_fluxes

   subroutine calculate_leakage_fluxes(Gamma_neut, leakage_fluxes)
	! Calculates the sources to the external volumes due to differences in cell area
	real( wp ), intent(in) :: Gamma_neut(0:Nx)
	real( wp ), intent(out) :: leakage_fluxes(5)
	real( wp ) 		:: G_loss(0:Nx)
	G_loss(0) = 0
	G_loss(Nx) = 0
	G_loss(1:Nx-1) = (A_int(1:Nx-1)-A_int(2:Nx))
	G_loss = G_loss * Gamma_neut
	where (G_loss < 0)
		G_loss = 0
	end where
	! Divertor one (fluxes are divided equally between PFR and CFR
	leakage_fluxes(1) = sum(G_loss(0:i_baffle(1)-1))/2.0d+0 + (A_wet(1)-A_int(1))*Gamma_neut(0) / 2.0d0
	leakage_fluxes(2) = leakage_fluxes(1)
	! Core CFR
	leakage_fluxes(3) = sum(G_loss(i_baffle(1):i_baffle(2)))
	! Divertor 2
	leakage_fluxes(4) = sum(G_loss(i_baffle(2)+1:Nx))/2.0d+0 - (A_wet(2)-A_int(Nx))*Gamma_neut(Nx) / 2.0d0
	leakage_fluxes(5) = leakage_fluxes(4)
	! the target leakage is there because the cb area is larger than the area of the neutral volume at that boundary

	return
   end subroutine calculate_leakage_fluxes
	
  !subroutine calculate_leakage_corrected_Flow_neutral_particles(flux,flow)
	! this subroutines calculates the number of particles that remain in the SOL when a flux establishes on one of the cell boundaries.
	! when other parameters are parsed, momentum, or energy, quanti of those are calculated
	!implicit none
	

 !  end subroutine calculate_leakage_corrected_Gamma_neutral_particles(flux,flow)

   subroutine calculate_sol_leakage_neutral_sources(flux, source)
	real( wp ), intent(in) :: flux(0:Nx)
	real( wp ), intent(out) :: source(Nx)
	integer 		:: ix
	source = 0.0d+0
	! loop over surfaces off cell boundaries
	! the surface of neutral cells that are at the same location xcb are the same.
	! depending on flow direction let particles flow. 
	do ix = 1,Nx-1
	 	if ((flux(ix) .ge. 0) .AND. (A_int(ix) .ge. A_int(ix+1))) then
			! S(n)
		 	source(ix) = source(ix) - A_int(ix) * flux(ix) / volumes(ix)
			! S(n+1)
		 	source(ix+1) = source(ix+1) + A_int(ix+1) * flux(ix) / volumes(ix+1)
		else if ((flux(ix) .lt. 0) .AND. (A_int(ix) .ge. A_int(ix+1))) then
			! S(n)
		 	source(ix) = source(ix) - A_int(ix+1) * flux(ix) / volumes(ix)
			! S(n+1)
		 	source(ix+1) = source(ix+1) + A_int(ix+1) * flux(ix) / volumes(ix+1)
		else if ((flux(ix) .ge. 0) .AND. (A_int(ix) .lt. A_int(ix+1))) then
			! S(n)
		 	source(ix) = source(ix) - A_int(ix) * flux(ix) / volumes(ix)
			! S(n+1)
		 	source(ix+1) = source(ix+1) + A_int(ix) * flux(ix) / volumes(ix+1)
		else 
			! S(n)
		 	source(ix) = source(ix) - A_int(ix) * flux(ix) / volumes(ix)
			! S(n+1)
		 	source(ix+1) = source(ix+1) + A_int(ix+1) * flux(ix) / volumes(ix+1)
		endif	
	 enddo
		source(1) = source(1) + A_int(1)*flux(0) / volumes(1)
	        source(Nx) = source(Nx) - A_int(Nx)*flux(Nx) / volumes(Nx)
   end subroutine calculate_sol_leakage_neutral_sources 

   subroutine calculate_extern_sources(sol2extern_flux, tar2extern_flux, extern_neutral_flux, extern2core_flux, sum_sol2extern_ion_flux, Source_extern, wall_association) ! , leakage_fluxes
	! calculates the source term for external reservoirs	 
	! units all in particles / s 
	implicit none
	! fluxes 2 extern are positive fluxes from extern2other are negative	
	real( wp ), intent(in) :: sol2extern_flux(5), tar2extern_flux(5), extern_neutral_flux(3), extern2core_flux(5)! , leakage_fluxes(5)
	real( wp ), intent(in) :: sum_sol2extern_ion_flux, wall_association(5)
	real( wp ), intent(out) :: Source_extern(5)     
	integer :: i	
	Source_extern = 0.0d+0
	!if( detect_nan) then
	!do i=1,5
	!if(isnan(sol2extern_flux(i))) stop 'sol2extern_flux NaN'
	!if(isnan(tar2extern_flux(i))) stop 'tar2extern_flux NaN'
	!if(isnan(extern2core_flux(i))) stop 'extern2core_flux NaN'
	!enddo
	!endif
        ! Sources of external neutral atoms
        Source_extern = ( sol2extern_flux + tar2extern_flux ) ! + leakage_fluxes * switch_neutral_leakage)
        !Source_extern(1) = (Source_extern(1) + extern_neutral_flux(3))
        Source_extern = Source_extern + &
			(/ extern_neutral_flux(1) , & ! positive direction for neutral flows outside the plasma is anti-clockwise
			  -extern_neutral_flux(2) , &
			   extern_neutral_flux(2) - extern_neutral_flux(3) , &
			   extern_neutral_flux(3) , &
			  -extern_neutral_flux(1) /)	
        Source_extern = Source_extern - extern2core_flux ! in particles per second
	Source_extern(3) = Source_extern(3) + sum_sol2extern_ion_flux
	Source_extern = Source_extern + wall_association
        return
    end subroutine calculate_extern_sources
  
   subroutine calculate_sol_sources( Nx, density, velocity, temperature, neutral, neutral_velocity, molecule, q_parallel, extern2sol_flux, extern2sol_mol, &
                                 Source_n, Source_v, Source_Q, Source_neutral, Source_vn, Source_molecule, elm_heat_load, elm_core_particle_source, &
				 Gamma_core2sol, core2sol_flux, core2sol_mol, sol2extern_ion_flux )
		
   ! this subroutine calculates the source terms of the discretized conservation equations
      implicit none
      integer,  intent(in)  :: Nx
      integer :: istop
      real(wp), intent(in)  :: density(Nx), velocity(Nx), temperature(Nx), neutral(Nx), neutral_velocity(Nx), molecule(Nx), q_parallel(0:Nx),  &
			       extern2sol_flux(Nx),extern2sol_mol(Nx), Gamma_core2sol ! sol2extern_flux(5), tar2extern_flux(5), extern_neutral_flux(3), &
			       !, sol2extern_mol(5), tar2extern_mol(5), extern_molecule_flux(3), & !, core_density, & 
			       !extern2core_flux(5), extern2core_mol(5), 
      real(wp), intent(in)  :: core2sol_flux(Nx), core2sol_mol(Nx), sol2extern_ion_flux(Nx)
      real(wp), intent(out) :: Source_n(Nx), Source_v(Nx), Source_Q(Nx), Source_neutral(Nx), Source_vn(Nx), Source_molecule(Nx) !,&
			       !Source_extern(5), Source_extern_mol(5) !, Source_core

      real(wp) :: rate_cx(Nx), rate_ion(Nx), rate_exc(Nx), rate_rec(Nx), rate_ree(Nx), rate_imp(Nx), rate_mom(Nx), mol_density(Nx), rate_mom_m(Nx), rate_cx_m(Nx), &
		rate_ion_m(Nx), rate_diss_m(Nx), rate_MAI_d(Nx), rate_MAR_Hmin(Nx), rate_MAD_Hmin(Nx), rate_da(Nx), Q_dr(Nx), Q_di(Nx), Q_diss_H2plus(Nx), Q_cx_Hmin(Nx), Q_ion_Hmin(Nx), &
		r_cx_m, r_ion_m, r_diss_m, r_diss_H2plus, r_diss_rec_H2plus, r_diss_ion_H2plus, f_cx, rate_MAR(Nx), rate_MAD(Nx), rate_MAI(Nx), rate_ene_m(Nx), maxrec,&
		r_da, r_ion_Hmin, r_cx_Hmin, relative_velocity, n_Hmin, n_H2plus
      !real(wp) :: Gamma_core2sol
      !real(wp) :: radial_sink(Nx)
      real(wp) :: mid_point, tmp_core_flux
      integer  :: ix, iix, i

      ! input variables for the elm simulation
      real(wp), intent(in)   :: elm_heat_load, elm_core_particle_source
	!write(*,*) "phys rout 460 velocity", velocity
	
	!#IF DEBUG
	if( detect_nan ) then
	do ix = 1,Nx
 	if (isnan(density(ix))) stop 'calc_sour: "density" is a NaN'
	if (isnan(velocity(ix))) stop 'calc_sour: "velocity" is a NaN'
	if (isnan(temperature(ix))) stop 'calc_sour: "temperature" is a NaN'
	if (isnan(neutral(ix))) stop 'calc_sour: "neutral" is a NaN'
	if (isnan(neutral_velocity(ix))) stop 'calc_sour: "neutral_velocity" is a NaN'
	if (isnan(molecule(ix))) stop 'calc_sour: "molecule" is a NaN'
	if (isnan(q_parallel(ix))) stop 'calc_sour: "q_parallel" is a NaN'
	if (isnan(extern2sol_flux(ix))) stop'calc_sour: "extern2sol_flux" is a NaN'
	if (isnan(extern2sol_mol(ix))) stop 'calc_sour: "extern2sol_mol" is a NaN'
	if (isnan(core2sol_flux(ix))) stop 'calc_sour: "core2sol_flux" is a NaN'
	if (isnan(core2sol_mol(ix))) stop 'calc_sour: "core2sol_mol" is a NaN' 
	if (isnan(sol2extern_ion_flux(ix))) stop 'calc_sour: "sol2extern_ion_flux"is a NaN'
	enddo
	endif
      Source_n = 0.0d+0
      Source_v = 0.0d+0
      Source_Q = 0.0d+0
      Source_neutral = 0.0d+0
      Source_vn = 0.0d+0
      Source_molecule = 0.0d+0
      !Source_extern = 0.0d+0
      !Source_extern_mol = 0.0d+0
	!write(*,*)'D_imp_con'
	!write(*,*) D_imp_con
	!write(*,*) 'internal_istep_rhs', internal_istep_rhs
      do ix = 1, Nx
         rate_cx(ix)  = density(ix) * neutral(ix) * charge_exchange(temperature(ix))
         rate_ion(ix) = density(ix) * neutral(ix) * ionization(density(ix),temperature(ix))
         rate_exc(ix) = density(ix) * neutral(ix) * excitation(density(ix),temperature(ix))
         rate_rec(ix) = density(ix) * density(ix) * recombination(density(ix),temperature(ix))
         rate_ree(ix) = density(ix) * density(ix) * recombenergy(density(ix),temperature(ix))
         rate_imp(ix) = density(ix) * density(ix) * impurity_radiation(temperature(ix),internal_istep_rhs,ix)
	 relative_velocity = velocity(ix) - neutral_velocity(ix)
	 rate_mom(ix) = density(ix) * neutral(ix) * momentum_transfer_atoms(temperature(ix),relative_velocity)
	 rate_mom_m(ix) = density(ix) * molecule(ix) * momentum_transfer_molecules(temperature(ix), velocity(ix))
	
	! molecular activated reactions (MAR, MAI and MAD)
	r_cx_m = charge_exchange_m(temperature(ix))
	r_ion_m = ionization_m(density(ix), temperature(ix))
	r_diss_m = dissociation_m(density(ix), temperature(ix))
	
	rate_diss_m(ix) = molecule(ix)*density(ix)*r_diss_m
	rate_cx_m(ix) = molecule(ix)*density(ix)*r_cx_m
	rate_ion_m(ix) = molecule(ix)*density(ix)*r_ion_m

	r_diss_H2plus = dissociation_H2plus(density(ix), temperature(ix))
	r_diss_rec_H2plus = dissociative_recombination_H2plus(density(ix), temperature(ix))
	r_diss_ion_H2plus = dissociative_ionization_H2plus(density(ix), temperature(ix))

	f_cx = r_cx_m/(r_cx_m + r_ion_m)
	if (ISNAN(f_cx)) then
		f_cx = 0.0d0
	endif 
	!if (ix .eq. 400) then
	!	write(*,*) 'f_cx', f_cx
	!endif
	rate_MAR(ix) = density(ix)*molecule(ix)*(r_cx_m+r_ion_m)/(r_diss_rec_H2plus+r_diss_H2plus+r_diss_ion_H2plus)*r_diss_rec_H2plus*f_cx
	rate_MAD(ix) = density(ix)*molecule(ix)*(r_cx_m+r_ion_m)/(r_diss_rec_H2plus+r_diss_H2plus+r_diss_ion_H2plus)*(f_cx*r_diss_H2plus+(1.0d0-f_cx)*r_diss_rec_H2plus)
	rate_MAI(ix) = density(ix)*molecule(ix)*(r_cx_m+r_ion_m)/(r_diss_rec_H2plus+r_diss_H2plus+r_diss_ion_H2plus)*(f_cx*r_diss_ion_H2plus+(1.0d0-f_cx)*(r_diss_H2plus+2.0d0*r_diss_ion_H2plus))	
	rate_ene_m(ix) = density(ix) * molecule(ix) * energy_molecules(density(ix), temperature(ix))

	! correct MAI in neutral particle source (because not all MAI creates a neutral atom):
	rate_MAI_d(ix) = density(ix)*molecule(ix)*(r_cx_m+r_ion_m)/(r_diss_rec_H2plus+r_diss_H2plus+r_diss_ion_H2plus)*(f_cx*r_diss_ion_H2plus+(1.0d0-f_cx)*r_diss_H2plus)
	!rate_MAI_d(ix) = density(ix)*molecule(ix)*(r_cx_m+r_ion_m)/(r_diss_rec_H2plus+r_diss_H2plus+r_diss_ion_H2plus)*(f_cx*r_diss_ion_H2plus+(1.0d0-f_cx)*r_diss_m)
	! MAR and MAD through Hmin
	r_da = dissociative_attachment(temperature(ix))
	r_ion_Hmin = ionization_Hmin(density(ix), temperature(ix))
	r_cx_Hmin=  charge_exchange_Hmin(density(ix), temperature(ix))

	rate_da(ix) = molecule(ix)*density(ix)*r_da
	rate_MAR_Hmin(ix) = molecule(ix)*density(ix)*r_da/(r_ion_Hmin+r_cx_Hmin)*r_cx_Hmin
	rate_MAD_Hmin(ix) = molecule(ix)*density(ix)*r_da/(r_ion_Hmin+r_cx_Hmin)*r_ion_Hmin

	! Energy sources due to reactions with H2plus and Hmin
	! Ionic densities
	n_H2plus = molecule(ix) * (r_cx_m+r_ion_m)/(r_diss_rec_H2plus+r_diss_H2plus+r_diss_ion_H2plus)
	n_Hmin = molecule(ix) * r_da/(r_ion_Hmin+r_cx_Hmin)

	Q_dr(ix) = - n_H2plus*density(ix)*r_diss_rec_H2plus*e_charge * (2.69d0-13.6d0)
	Q_di(ix) = - n_H2plus*density(ix)*r_diss_ion_H2plus*e_charge * (2.69d0+13.6d0)
	Q_diss_H2plus(ix) = - n_H2plus*density(ix)*r_diss_ion_H2plus*e_charge * 2.69d0
	Q_cx_Hmin(ix) = -n_Hmin*density(ix)*r_cx_Hmin*e_charge*0.754d0
	Q_ion_Hmin(ix) = -n_Hmin*density(ix)*r_ion_Hmin*e_charge*(0.754d0+13.6d0)

      enddo
	! source neutral + source_n + source_molecule*2.0 = 0 % da en cx, ion apart. 
	! 
	! below should be the same
	! source:
	!write(*,*) 'ion', rate_ion_m(400)
	!write(*,*) 'cx', rate_cx_m(400)
	! if only MAD active
	!write(*,*) 'mad', rate_MAD(400)
	! if only MAR active  neutrals         ions
	!write(*,*) 'mar', (rate_MAR(400)*3.0d0-rate_MAR(400)) / 2.0d0
	! if only MAI active
	!write(*,*) 'mai', (rate_MAI(400) + rate_MAI_d(400)) / 2.0d0

	! check the Hmin channel
	!write(*,*) 'da', rate_da(400) 
	! the MAR Hmin channel
	!write(*,*)'Hm_MAR', (3.0d0*rate_MAR_Hmin(400)-rate_MAR_Hmin(400))/2.0d0
	! MAD Hmin channel
	!write(*,*)'Hm_MAD', rate_MAD_Hmin(400)

      ! neutral particle sources
      !write(*,*) "core2sol = " !, core2sol_flux / volumes + extern2sol_flux / volumes
      Source_neutral = rate_rec + (3.0d+0*rate_MAR+ 3.0d+0*rate_MAR_Hmin) - rate_ion + 2.0d0*rate_diss_m  + rate_MAI_d + 2.0d+0*(rate_MAD+rate_MAD_Hmin) & 
			+ extern2sol_flux / volumes    + core2sol_flux /  volumes
			! + gas_puff_profile * D_gas(internal_istep_rhs) &
      ! molecule sources (note that cx, ionization and dissociation (and dissociative attachment) are the molecule destroying reactions)
      Source_molecule = - rate_cx_m - rate_ion_m - rate_diss_m - rate_da + gas_puff_profile * D_gas(internal_istep_rhs) & 
			+ extern2sol_mol / volumes  + core2sol_mol / volumes
      ! we are adding the gas as molecules only !	
      
		
      ! plasma particle sources
	!write(*,*) '599 rate ion:', rate_ion    
	!write(*,*) '599 rrate_rec:', rate_rec   
	!write(*,*) '599 rate_MAR:', rate_MAR   
	!write(*,*) '599  rate_MAI:',  rate_MAI    
	!write(*,*) '599  rate_MAR_Hmin:',  rate_MAR_Hmin    
      Source_n = rate_ion - rate_rec - rate_MAR + rate_MAI - rate_MAR_Hmin - sol2extern_ion_flux / volumes
      
      ! the momentum sources
      !write(*,*) "rate_cx = ", maxval(rate_cx), minval(rate_cx)
      !write(*,*) "rate_rec= ", maxval(rate_rec)
      !write(*,*) "rate_ion= ", maxval(rate_ion)
      Source_v = - mass * velocity * ( rate_cx_m + rate_da + 0.6667d0*rate_mom_m + rate_rec ) &
		 - mass * (velocity - neutral_velocity) * ( rate_cx + 0.5d0*rate_mom ) + mass * neutral_velocity * rate_ion &
		 - mass * velocity * sol2extern_ion_flux / volumes

      ! the energy sources (only internal energy)
     if ( neutral_energy .lt. 0.0d+0 ) then ! use energy = neutral_energy * temperature_plasma as neutral energy (i.e. the energy is a multiplier with temperature
      	Source_Q = - ((1.5d+0 + neutral_energy )* e_charge * temperature  ) * ( rate_cx ) ! neutrals have energy and reduce losses from plasma in charge exchange ( energy in e*ev) (1.5 - 1.5 = 0.5)
	! neutral energy also adds to the ionization term        
	Source_Q = Source_Q + ( rate_ion * -1.0d+0 * neutral_energy * e_charge * temperature )
      else	! use fixed energy
      	Source_Q = - (1.5d+0 * e_charge * temperature - e_charge*neutral_energy ) * ( rate_cx ) ! neutrals have energy and reduce losses from plasma in charge exchange ( energy in e*ev)
	! neutral energy also adds to the ionization term  
        Source_Q = Source_Q + ( rate_ion * neutral_energy * e_charge )
      endif
      Source_Q = Source_Q - (1.5d+0 * e_charge * temperature) * ( 2.0d+0 * rate_rec )  ! the factor 2 for the recombination rate is for the loss of both an ion and electron
      ! add energy source associated with particle source from ionization (like a source of internal energy from effective friction with neutrals)
      Source_Q = Source_Q + ( rate_ion * ( 0.5d+0 * mass * velocity**2 ) ) 
      if ( switch_excitation .eq. 0.0d+0 ) then
         Source_Q = Source_Q - rate_ion * e_charge * ( energy_loss_ion ) 
      else
         Source_Q = Source_Q - switch_excitation * rate_exc * e_charge ! note excitation rate is in eV m^3 / s
      endif
      ! add energy losses associated with radiative and 3-body recombination (the 13.6 eV potential energy per recombination event is added explicitly here)
      Source_Q = Source_Q - ( rate_ree - 13.6d+0 * rate_rec ) * e_charge ! note recombination loss rate is in eV m^3 / s
      ! add impurity radiation losses
      Source_Q = Source_Q - switch_impurity_radiation * rate_imp * e_charge ! note impurity radiation loss rate is in eV m^3 / s
      !write(*,*) rate_ion, Source_Q

      ! Add molecular contributions
      Source_Q = Source_Q - (1.5d+0 * e_charge * temperature) *  2.0d+0 * (rate_MAR+rate_MAR_Hmin) & 
			+   ( rate_MAI * ( 0.5d+0 * mass * velocity**2 ) ) &
			-  rate_ene_m*e_charge &
			- 1.786d0*rate_cx_m*e_charge & 
			- 3.724d0*e_charge*rate_da 

      ! Add contributions from H2plus and Hmin
      Source_Q = Source_Q + Q_dr + Q_di + Q_diss_H2plus + Q_cx_Hmin + Q_ion_Hmin

      ! Add the effect of radial particle losses across the core flux tube
      Source_Q = Source_Q - (1.5d+0 * e_charge * temperature) * 2.0d+0 * sol2extern_ion_flux / volumes ! factor 2 accounts for loss of both ion and electron (otherwise not equal charge?)
       
      ! Neutral momentum sources
      Source_vn = mass * velocity * (rate_rec + rate_cx_m + rate_da) &
                  + mass * (velocity - neutral_velocity) * ( rate_cx + rate_mom ) &
	          - mass * neutral_velocity * rate_ion &
                  - mass * neutral_velocity * neutral * Area_extern / volumes / neutral_residence_time * ato_res_asy ! so neutrals actually lose more momentum than the exchanged fluxes 
      ! should the target be included from recycling fluxes? (well for now its all molecules, but still.


      ! remove spikes that cause problems during integration of the ODE
      if( filter_sources ) then
          call spike_filter( Nx, Source_neutral)
          call spike_filter( Nx, Source_n )
          call spike_filter( Nx, Source_v )
          call spike_filter( Nx, Source_Q )
	  call spike_filter( Nx, Source_vn )
	  call spike_filter( Nx, Source_molecule )
      endif

 	! section adding sources along the core-SOL boundary (when present in the grid)
	!write(*,*) 'Source_core', Source_core
	if( evolve_core .eq. 1 ) then
	  tmp_core_flux = Gamma_core2sol
	  !write(*,*) 'Source_core'
	else
	  tmp_core_flux = D_Gamma_core(internal_istep_rhs)
	  !write(*,*) 'DGamma'
	endif
	!write(*,*) 'Gamma_core = ', D_Gamma_core(internal_istep_rhs)
	!write(*,*) 'tmp_core_flux', tmp_core_flux
	!write(*,*) 'I_core_source_profile_n', I_core_source_profile_n
	!write(*,*) 'elm_core_particle_source', elm_core_particle_source
	if( switch_only_core_source_Q .gt. 0.5d+0 ) then
	Source_Q = 0.0d+0
	!write(*,*) 'Source_Q = 0'
	endif
	!write(*,*) '666 source_n:', Source_n
       !write(*,*) 'before Source_Q 10 n-10 = ', Source_Q(10), Source_Q(Nx-10)
       if( L_core_SOL .gt. 0.0d+0 ) then
          if (wide_core_profile == 0) then
          Source_n(i_Xpoint(1):i_Xpoint(2)) = Source_n(i_Xpoint(1):i_Xpoint(2)) & 
		+ tmp_core_flux*(1-core_far_sol_feedthrough) * I_core_source_profile_n / volumes(i_Xpoint(1):i_Xpoint(2))
		!+ (tmp_core_flux + elm_core_particle_source) * I_core_source_profile_n / volumes(i_Xpoint(1):i_Xpoint(2))
		
		
          else
          Source_n(i_baffle(1):i_baffle(2)) = Source_n(i_baffle(1):i_baffle(2)) & 
		+ tmp_core_flux*(1-core_far_sol_feedthrough)  * I_core_source_profile_n / volumes(i_baffle(1):i_baffle(2))
		!+ (tmp_core_flux + elm_core_particle_source) * I_core_source_profile_n / volumes(i_baffle(1):i_baffle(2))
				
          endif    
          Source_Q(i_Xpoint(1):i_Xpoint(2)) = Source_Q(i_Xpoint(1):i_Xpoint(2)) & 
		+ D_Q_core(internal_istep_rhs) * I_core_source_profile_Q / volumes(i_Xpoint(1):i_Xpoint(2)) 
		!+ (D_Q_core(internal_istep_rhs)+elm_heat_load) * I_core_source_profile_Q / volumes(i_Xpoint(1):i_Xpoint(2)) 
		
		
       endif
       !write(*,*) 'Q_core(itime) = ', D_Q_core(internal_istep_rhs)
       !write(*,*) 'after Source_Q 10 n-10 = ', Source_Q(10), Source_Q(Nx-10)
       !write(*,*) 'core source profile Q = ', I_core_source_profile_Q
	!write(*,*) 'core source profile = ', I_core_source_profile_n
	! something goes wrong with Source_n using the wide_core_profile
	!write(*,*) 'size Source_n =', size(Source_n)
	!write(*,*) 'calc_sour: Source_n =', Source_n
	!#IF DEBUG
	!write(*,*) '682 source_n:', Source_n
	if( detect_nan ) then	
	do ix = 1,Nx
 	if (isnan(Source_n(ix))) stop  'calc_sour: "Source_n" is a NaN'
	if (isnan(Source_Q(ix))) stop 'calc_sour: "Source_Q" is a NaN'
	if (isnan(Source_v(ix))) stop 'calc_sour: "Source_v" is a NaN'
	if (isnan(Source_neutral(ix))) stop 'calc_sour: "Source_neutral" is a NaN'
	if (isnan(Source_vn(ix))) stop 'calc_sour: "Source_vn" is a NaN'
	if (isnan(Source_molecule(ix))) stop 'calc_sour: "Source_molecule" is a NaN'
	enddo
	endif	
	!#ENDIF
      return
   end subroutine calculate_sol_sources


   subroutine step_div1d(ys,yr,yc, start_time, end_time, itask, istate, options) 
   ! this subroutine runs dvode using the nested right_hand_side routine
   use dvode_f90_m ! this is only scoped in step_div1d, to call dvode_f90 and pass options 
   implicit none
   integer :: itask, istate, neq, internal_istep_rhs, attempt, ii
   real( wp ), intent(inout) :: ys(6*Nx)
   real( wp ), intent(inout) :: yr(10)
   real( wp ), intent(inout) :: yc
   real( wp ) :: y(6*Nx+10), start_time, end_time, tmp(5)
   type (VODE_OPTS) :: options
   
   ! states of the SOL and surroundings
   real(wp)   :: density(Nx), velocity(Nx), temperature(Nx), neutral(Nx), neutral_velocity(Nx), molecule(Nx) 
   real(wp)   :: extern_neutral_density(5),  extern_molecule_density(5)		  
   real(wp)   :: core_density, tmp_core_neutral
   ! fluxes between sol and extern
   real(wp)   :: extern2sol_tmp, sol2extern_flux(5), sol2extern_mol(5)
   real(wp)   :: dummy5(5), zeroNx(Nx), zero3(3),dummy33(3)
   real(wp)   :: extern2sol_flux(Nx) 
   real(wp)   :: extern2sol_mol(Nx)
   real(wp)   :: extern2core_shinethrough(Nx) 
   real(wp)   :: extern2core_shinethrough_mol(Nx)
   ! fluxes shared as interface between core and sol
   real(wp) :: core2sol_flux(Nx), sol2core_flux ! [1/s]
   real(wp) :: core2sol_mol(Nx), sol2core_mol ! [1/s]
   real(wp) :: Gamma_core2sol, Source_core
   ! fluxes shared as interface between extern and core
   real(wp) :: extern2core_flux(5), extern2core_mol(5)
   !real(wp) :: sol2extern_flux
   ! fluxes shared as interface between sol and extern
   !     -->   backgrounds are solves together with SOL for now, so no shared fluxes there
   internal_istep_rhs = min(int((start_time - D_start_time)/D_delta_t+1),D_nout,nout)
   internal_istep_rhs = max(1,internal_istep_rhs)
   !write(*,*) 'F rout: internal_istep_rhs', internal_istep_rhs
   !write(*,*) 'F rout in T(10)=',temperature(10)
   ! states (log) normalized: THESE STATES ARE USED BY DVODE_F90 and STEP_CORE (although we could also step core_reservoir_rhs)
   extern2sol_flux = 0.0d+0
   sol2extern_flux = 0.0d+0
   extern2core_shinethrough = 0.0d+0
   extern2core_shinethrough_mol = 0.0d+0
   zero3 = 0.0d+0
   zeroNx = 0.0d+0

   ! unpack solution vector
   call ys2nvt( Nx, ys, density, velocity, temperature, neutral, neutral_velocity, molecule)
   call yr2nr(yr, extern_neutral_density, extern_molecule_density)
   call yc2nc(yc, core_density)

   ! use prescribed densities for external reservoirs and core neutrals if they are note evolved.
      do ii = 1, 5
      if( evolve_background(ii) .eq. 0 ) then 
	extern_neutral_density(ii) = D_nb(ii,internal_istep_rhs)
	extern_molecule_density(ii)  = D_mb(ii,internal_istep_rhs)	
      endif
      enddo
      if( evolve_molecule .eq. 0 ) then
	extern_molecule_density = D_mb(1:5,internal_istep_rhs)	 	
      endif
   if( evolve_core_neutral .eq. 0 ) then
   	tmp_core_neutral = D_core_neutral_density(internal_istep_rhs)
	!write(*,*) 'F tmp_core_neutral =', tmp_core_neutral
   else
	stop 'stopped IN step_div1d: core neutral evolution not implemented'
   endif 

   ! Fixate flows between CORE and SOL (these will be fixed when stepping forward in time)
   call calculate_sol_extern_neutral_fluxes(extern2sol_flux, sol2extern_flux, zero3, & !Out
					extern_neutral_density, neutral, neutral_residence_time, extern_neutral_ex )!, & ! In
   ! Fixate molecule flow as well
   call calculate_sol_extern_neutral_fluxes(extern2sol_mol, sol2extern_mol, zero3, & !Out
					extern_molecule_density, molecule, molecule_residence_time, extern_molecule_ex )!, & ! In

	!write(*,*), 'IN div1d step: extern2sol_flux',extern2sol_flux
   ! use core ionization fraction
   !write(*,*) 'extern2sol_mol', extern2sol_mol
   call calculate_core_ionization_neutral_flux(extern2core_shinethrough, extern2sol_flux, core_ionization_fraction )
   call calculate_core_ionization_neutral_flux(extern2core_shinethrough_mol, extern2sol_mol, core_ionization_fraction_mol )
   !write(*,*) 'extern2core_shinethrough_mol', extern2core_shinethrough_mol
   !write(*,*) 'extern2sol_mol', extern2sol_mol


   ! consider the core ionization fraction when calculating the interchange fluxes
   call calculate_sol2core_neutral_flux(core2sol_flux, sol2core_flux, tmp_core_neutral, neutral,  core_sol_neutral_ex, extern2core_shinethrough)
   ! molecules do also shine through but with different factor
   call calculate_sol2core_neutral_flux(core2sol_mol , sol2core_mol , tmp_core_neutral, molecule, core_sol_molecule_ex, extern2core_shinethrough_mol)
   ! fluxes that go directly into the core
   call calculate_extern2core_fluxes(extern2core_flux,extern_neutral_density,core_ext_neutral_pump)
	!write(*,*) 'extern2core_flux', extern2core_flux
	!write(*,*) 'extern_neutral_density', extern_neutral_density
	!write(*,*) 'core_ext_neutral_pump', core_ext_neutral_pump
   call calculate_extern2core_fluxes(extern2core_mol,extern_molecule_density,core_ext_molecule_pump)
	!write(*,*) 'core_ext_neutral_pump',core_ext_neutral_pump
	!write(*,*) 'extern2core_flux',  extern2core_flux
   ! add the core ionization fraction to previous calculations (if the extimes are zero this will only use the ionization fraction and vice versa) 
   !write(*,*) 'F phys rout: c2s-flux', core2sol_flux(100)
   !write(*,*) 'F phys rout: e2c-flst', extern2core_st(100)
   call calculate_core2sol_ion_flux( core_density, Gamma_core2sol)
 
   ! solve the core evolution by means of forward euler
   if( evolve_core .eq. 1 ) then   
     !write(*,*) 'IN step_div1d : core neutral, sol2core, sol2coremol', tmp_core_neutral, sol2core_flux, sol2core_mol
     !write(*,*) 'IN step_div1d : neutral  ', neutral
   
     call calculate_core_source( extern2core_flux, extern2core_mol, Gamma_core2sol, sol2core_flux, sol2core_mol, Source_core )
     !write(*,*) 'IN evolve_core: sol2core_flux', sol2core_flux
     !write(*,*) 'IN evolve_core: sol2core_mol', sol2core_mol
     !write(*,*) 'IN evolve_core: gamma core2sol', Gamma_core2sol
     !write(*,*) 'IN evolve_core: extern2core_flux', sum(extern2core_flux)
     !write(*,*) 'IN evolve_core: extern2core_mol', sum(extern2core_mol)

     !write(*,*) 'IN step_div1d : source_core', Source_core
     call fe_ncore(yc, Source_core, start_time, end_time)
   else
     Gamma_core2sol = D_Gamma_core(internal_istep_rhs) !should be used in rhs. 
     Source_core = 0.0d+0 ! not relevant
   endif
    
   ! solve for the SOL and external volumes with DVODE at once
   write(*,*) 'F rout call DVODE: t0', start_time, 't1', end_time
   neq = 6*Nx + 10
   y(1:6*Nx) = ys ! sol vector
   y(6*Nx +1 : 6*Nx+10) = yr ! reservoir vector
   call dvode_f90( right_hand_side, neq, y, start_time, end_time, itask, istate, options )
        if(istate .ne. 2) then ! restart if unsuccesfull
        attempt = 0
           do while( istate .ne. 2 .and. attempt .le. max_attempts )
                attempt = attempt + 1
                istate = 1
                call dvode_f90( right_hand_side, neq, y, start_time, end_time, itask, istate, options )
           enddo
           if( istate .ne. 2 ) then
                call error_report_timestep( istate )
           else
                write(*,*) 'F success on ', attempt, 'th attempt. continuing'           
           endif     
        endif
   ! write solution back to original vectors
   ys = y(1:6*Nx)
   yr = y(6*Nx +1 : 6*Nx+10)	
   !call nc2yc(core_density, yc)
   ! IF ALL went well, one obtains y and yc as solutions delta_t advanced in time
   contains

   subroutine fe_ncore(yc, Source_core, start_time, end_time) 
	! routines solves for ncore by means of forward euler
   	implicit none
   	real(wp), intent(inout) :: yc
   	real(wp), intent(in) :: Source_core
    	real(wp), intent(in) :: start_time, end_time
   	real(wp) :: core_density
        call	yc2nc(yc,core_density)
   	core_density = core_density + Source_core*( end_time - start_time) ! forward euler
	core_density = max(core_density, minimum_density) ! dont let it go below 0, this will give nan.
   	call 	nc2yc(core_density, yc)
	return ! returns the log normalized yc value for the core (we do not want to return core_density as that should be updated at the end of step_div1d
   end subroutine fe_ncore
   
   subroutine right_hand_side( neq, physical_time, y, ydot )
   ! this subroutine calculates the right hand side ydot of the discretized conservation equations
   ! note that y and ydor are normalized, but the arrays density, velocity, temperature, and neutral are not!
   ! ew 01-03-2021: modified to take account of flux expansion
      use grid_data, only: Nx
      implicit none
      integer :: ii
      integer,  intent(in)  :: neq
      real(wp), intent(in)  :: physical_time, y(neq) 
      real(wp), intent(out) :: ydot(neq)
      real(wp) :: tmpn
!      integer               :: Nx ! Nx is known from numerics_parameters
      real(wp) 		    :: ys(6*Nx), yr(10)
      integer 		    :: ix, i
      ! we here define the internal variables that right_hand_side uses, these do not overwrite the global ones in the DIV1D program
      real(wp)              :: density(Nx), velocity(Nx), temperature(Nx), neutral(Nx)      ![1/m3] ,[m/s]    ,[eV]   ,[1/m3]
      real(wp) 		    :: neutral_velocity(Nx), molecule(Nx)   
      real(wp) 		    :: velocity_cb(0:Nx), temperature_cb(0:Nx), csound_cb(0:Nx), neutral_velocity_cb(0:Nx)
      real(wp)              :: Gamma_n(0:Nx), Gamma_mom(0:Nx), q_parallel(0:Nx) ![1/m2s],[kg/ms2] ,[J/m2s],[1/m2s]
      real(wp) 		    :: Gamma_neutral(0:Nx), Gamma_molecule(0:Nx), Gamma_mom_neutral(0:Nx)
      real(wp)              :: Source_n(Nx), Source_v(Nx), Source_Q(Nx), Source_neutral(Nx) ![1/m3s],[kg/m2s2],[J/m3s],[1/m3s] 
      real(wp) 		    :: Source_vn(Nx), Source_molecule(Nx), Source_extern(5), Source_extern_mol(5)
      real(wp)  	    :: dpdx(Nx), vdpdx(Nx), dpndx(Nx), vndpndx(Nx), flux_source_neutral(Nx), flux_source_molecule(Nx), flux_source_neutral_momentum(Nx)
      real(wp) 		    :: extern_neutral_density(5),  extern_molecule_density(5)
      real(wp)              :: extern_neutral_flux(3), sol2extern_flux(5), tar2extern_flux(5), extern2sol_flux(Nx) ! [1/s]
      real(wp)              :: extern_molecule_flux(3), sol2extern_mol(5), tar2extern_mol(5), extern2sol_mol(Nx) ! [1/s]
      real(wp)              :: sol2extern_ion_flux(Nx), sum_sol2extern_ion_flux
      real(wp)   	    :: neutral_pump(5), molecule_pump(5), nbdot(10)
      !real(wp) 		    :: leakage_fluxes_n(5), leakage_fluxes_m(5)
      real(wp) 		    :: atom_association_sink(5), molecule_association_source(5)
      ! In the above following comes from outside RHS: Gamma_core2sol, extern2core_flux, extern2core_mol
      ! do not redefine here, they will overwrite 
      real(wp)              :: Diff_neutral(Nx),  Diff_mol(Nx)   ! [m2/s]
      real(wp)              :: csound_target(2), q_sheath, v0, Gmom0, v0n, Gmom0n
      real(wp) 	            :: dummy_core2sol_flux(Nx)
      ! input variables for the elm simulation

      real(wp)              :: elm_heat_load, elm_density_change, elm_core_particle_source
      real(wp)		    :: tmp_neutral_energy(Nx)
      internal_istep_rhs = min(int((physical_time - D_start_time)/D_delta_t),D_nout,nout)
      internal_istep_rhs = max(1,internal_istep_rhs)
	!write(*,*) 'RHS internal_istep_rhs = ', internal_istep_rhs	
      
      ! first tranform the solution vector to (density velocity temperature neutral-density)
	ys = y(1:6*Nx)
	yr = y(6*Nx+1:6*Nx+10)
      call ys2nvt( Nx, ys, density, velocity, temperature, neutral, neutral_velocity, molecule)	
      call yr2nr( yr, extern_neutral_density, extern_molecule_density )

      ! use prescribed densities for external reservoirs if they are note evolved.
      !if( evolve_background .eq. 0 ) then 
!	extern_neutral_density	= D_nb(1:5,internal_istep_rhs)
!	extern_molecule_density = D_mb(1:5,internal_istep_rhs)	
 !     endif
      do ii = 1, 5
      if( evolve_background(ii) .eq. 0 ) then 
	extern_neutral_density(ii) = D_nb(ii,internal_istep_rhs)
	extern_molecule_density(ii)  = D_mb(ii,internal_istep_rhs)	
      endif
      enddo
      if( evolve_molecule .eq. 0 ) then
	extern_molecule_density = D_mb(1:5,internal_istep_rhs)	 	
      endif
	! here we also need to unpack external_neutral_density and extern_molecule_density?
	!write(*,*) 'velocity', velocity
      call cc2cb(Nx,velocity,velocity_cb)
      call cc2cb(Nx,temperature,temperature_cb)
      call cc2cb(Nx,neutral_velocity, neutral_velocity_cb)
	 
      !write(*,*) 'density = ', density
      !write(*,*) 'velocity =', velocity
      !write(*,*) 'temperature =', temperature
      !write(*,*) 'neutral density =', neutral
 	
      ! prints below should be constant throughout internal calls 
      !write(*,*) 'RHS Gamma_core2sol = ', Gamma_core2sol
      !write(*,*) 'RHS extern2core_flux = ', extern2core_flux	
      !write(*,*) 'RHS extern2core_mol = ',extern2core_mol
	
      ! calculate the ELM heat flux and particle flux
      !call simulate_elm(elm_heat_load, elm_density_change, elm_core_particle_source, physical_time)
      ! calculate the fluxes
      call calculate_sol_fluxes( Nx, density, velocity, temperature, neutral, neutral_velocity, molecule, &
				extern_neutral_density,  extern_molecule_density,  &
				Gamma_n, Gamma_mom, q_parallel, Gamma_neutral, Gamma_molecule, Gamma_mom_neutral,& 
				elm_heat_load, & 
                                extern_neutral_flux, extern_molecule_flux, &
				sol2extern_flux, sol2extern_mol, & 
				extern2sol_flux, extern2sol_mol, &
				tar2extern_flux, tar2extern_mol)!, &
	!write(*,*) 'rhs: neutral_flux =',extern_neutral_flux
	!write(*,*) 'rhs: molecule =',extern_molecule_flux
      ! 					 Nx			1			Nx
      call calculate_sol2extern_ion_flux(sol2extern_ion_flux, sum_sol2extern_ion_flux, density, Gamma_core2sol) ! note that Gamma_core2sol should not be changed inside RHS
      !						 Nx 		 1			1		Nx		1  
	!write(*,*) 'RHS: Gamma_mom(Nx-10)', Gamma_mom(Nx-10)
      ! calculate the sources ! note here we use the non-dummy fluxes calculated outside rhs and fixed for the entire timestep
      call calculate_sol_sources( Nx, density, velocity, temperature, neutral, neutral_velocity, molecule, q_parallel, extern2sol_flux, extern2sol_mol, & 
                                 Source_n, Source_v, Source_Q, Source_neutral, Source_vn, Source_molecule, elm_heat_load, elm_core_particle_source, &
				 Gamma_core2sol, core2sol_flux, core2sol_mol, sol2extern_ion_flux ) 
      !atoms [#/s]
      !call calculate_leakage_fluxes(Gamma_neutral, leakage_fluxes_n)

      !molecules [#/s]
      !call calculate_leakage_fluxes(Gamma_molecule, leakage_fluxes_m)

      !atoms [D/s] --> molecules [D2/s]
      call calculate_wall_association(atom_association_sink, molecule_association_source,extern_neutral_wall_area, extern_neutral_density,neutral_residence_time,wall_association_probability)

      ! atoms [#/s]
      call calculate_extern_sources(sol2extern_flux, tar2extern_flux, extern_neutral_flux, extern2core_flux, (1.0d+0-mol_rec)*sum_sol2extern_ion_flux, &
				    Source_extern, atom_association_sink)
      ! molecules [#/s]										! factor 0.5 here for 
      call calculate_extern_sources(sol2extern_mol, tar2extern_mol, extern_molecule_flux, extern2core_mol,0.5d+0*mol_rec*sum_sol2extern_ion_flux, &
				    Source_extern_mol, molecule_association_source)
			! In the above following comes from outside RHS: Gamma_core2sol, extern2core_flux, extern2core_mol
			! Far sol recycling fluxes have the same mol_recycle fraction as the target fluxes? 

      ! atoms [#/s]
      call calculate_extern_pump(neutral_pump,extern_neutral_density,pump_rate_n)
		!write(*,*) 'rhs: neutral_pump =', neutral_pump
		!write(*,*) 'rhs: ext_neutral_den =', extern_neutral_density
 
      ! molecules [#/s]
      call calculate_extern_pump(molecule_pump,extern_molecule_density,pump_rate_m)

       !write(*,*) 'Gamma_n =', Gamma_n
       !write(*,*) 'Gamma_mom =', Gamma_mom
       !write(*,*) 'q_parallel =', q_parallel
       !write(*,*) 'RHS: Source_n =', Source_n
       !write(*,*) 'Source_v =', Source_v
        !write(*,*) 'Source_Q =', Source_Q
	!write(*,*) 'Temperature = ', temperature
       !write(*,*) 'density = ', density
       !write(*,*) 'Source_neutral =', Source_neutral

        ! call cc2cb(Nx,density	,density_cb)   
        call cc2cb(Nx,velocity	,velocity_cb)
        call cc2cb(Nx,temperature	,temperature_cb)
        ! call cc2cb(Nx,neutral 	,neutral_cb)
        ! sound velocity at the target(s)
	csound_cb = sqrt( 2.0d+0 * e_charge * max(temperature_cb,minimum_temperature) / mass )   
	!write(*,*) 'temp, temp_cb' 	      
 	!write(*,*) 'csound/velocity', csound_cb(Nx), '/', velocity_cb(Nx)


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
             !ydot(1) = density_ramp_rate + elm_density_change + D_dneu(internal_istep_rhs)  ! [1/ (m^3 s)]
 	     ydot(1) = density_ramp_rate  + D_dneu(internal_istep_rhs)  ! [1/ (m^3 s)]
         endif
      !write(*,*) 'ydot(density) =', ydot(0*Nx+1:1*Nx) ! ---------------------------------------------------------------------

      ! --------------------------------------------- ydot for the momentum equation ------------------------------------------
         ydot(Nx+1:2*Nx) = switch_momentum_source * Source_v(1:Nx) ![kg/(m^2 s^2]
         ! add the momentum flux term using the flux as calculated in calculate_fluxes
         ! Gamma_mom(i) contains the flux at i+1/2
	 !write(*,*) 'RHS2: Gamma_mom(Nx-10)', Gamma_mom(Nx-10)
	 !write(*,*) 'B_f=',B_field(Nx-10), ', delta_xcb=', delta_xcb(Nx-10),', dgamdx=', B_field(Nx-10) * (Gamma_mom(Nx-10)-Gamma_mom(Nx-11))/delta_xcb(Nx-10) 	
         do ix = 2, Nx
            ydot(Nx+ix) = ydot(Nx+ix) - B_field(ix) * (Gamma_mom(ix)-Gamma_mom(ix-1))/delta_xcb(ix) 
         enddo
         
         if( L_core_SOL .gt. 0.0d+0 ) then
             ! apply boundary condition at mid point (i.e. flux = 0) or as calculated for the target at x=0
             ydot(Nx+1) =  ydot(Nx+1) - B_field(1) * (Gamma_mom(1)-Gamma_mom(0))/delta_xcb(1)    ! [kg/(m^2 s^2)]
         else
             ! apply boundary condition at the X-point, as following from the constant density n(1)
             ! velocity at i = 0:  v(0) = v(2) - 2 (S_n(1)+density_ramp_rate) delta_xcb(1) / n(1)
             v0 = velocity(2) - 2.0d+0 * (Source_n(1)+ydot(1)) * delta_xcb(1) / density(1) 
	     ! momentum flux at i = 0 : Gmom0 = m n (1/4)(v(0)+v(1))**2
             Gmom0 = mass * density(1) * (v0 + velocity(1))**2/4.0d+0       ! [kg/(m   s^2)]
             ydot(Nx+1) = ydot(Nx+1) - (Gamma_mom(1)-Gmom0)/delta_xcb(1)    ! [kg/(m^2 s^2)]
         endif
         ! add the pressure term in the internal region using central differencing: NB pressure =2/3 * y(2*Nx+1:3*Nx)
	 dpdx = 0.0d+0
	 dpdx(1) 	= 	  		  (y(2*Nx+2)       -y(2*Nx+1)       )*energy_norm/1.5d+0/delta_x(1)   
	 dpdx(2:Nx-1) 	= 		 0.5d+0 * (y(2*Nx+3:3*Nx  )-y(2*Nx+2:3*Nx-1))*energy_norm/1.5d+0/delta_x(2:Nx-1)
	 dpdx(2:Nx-1)	= dpdx(2:Nx-1) + 0.5d+0 * (y(2*Nx+2:3*Nx-1)-y(2*Nx+1:3*Nx-2))*energy_norm/1.5d+0/delta_x(1:Nx-2)
	 dpdx(Nx) 	= 	  		  (y(3*Nx)         -y(3*Nx-1)       )*energy_norm/1.5d+0/delta_x(Nx-1)
	 ydot(Nx+1:2*Nx) = ydot(Nx+1:2*Nx) - dpdx

         do ix = 2, Nx-1
              ! add effect of numerical viscosity
            ydot(Nx+ix) = ydot(Nx+ix) + viscosity*(velocity(ix+1) + velocity(ix-1)-2.0d+0*velocity(ix))
         enddo

         ! apply boundary condition at the sheath entrance, i=Nx: (linearly extrapolate velocity beyond the sheath)
         ydot(2*Nx) = ydot(2*Nx) + viscosity*(2.0d+0*max(csound_cb(Nx),velocity_cb(Nx)) + velocity(Nx-1)-3.0d+0*velocity(Nx))  ! add numerical viscocity
         if( L_core_SOL .gt. 0.0d+0 .and. X_core_SOL .gt. 0.0d+0 ) then
             ! apply boundary contion at the x=0 sheath for the pressure and viscosity terms
             ydot(Nx+1) = ydot(Nx+1) + viscosity*(2.0d+0*min(-csound_cb(0),velocity_cb(0)) + velocity(2)-3.0d+0*velocity(1))  ! add numerical viscocity
         endif
      !write(*,*) 'ydot(momentum) =', ydot(1*Nx+1:2*Nx) ! -----------------------------------------------------------------------
      
      ! ------------------------------------------------ ydot for the energy equation --------------------------------------------
	!write(*,*) 'rhs, SQ:', Source_Q         
	ydot(2*Nx+1:3*Nx) = switch_energy_source * Source_Q(1:Nx) ! [J/ (m^3 s)]
         !write(*,*) 'rhs, qpar:', q_parallel(Nx-10:Nx)
	! add the heat flux term including all boundaries
         ydot(2*Nx+1:3*Nx) = ydot(2*Nx+1:3*Nx) - switch_energy_flux * B_field(1:Nx) * (q_parallel(1:Nx)-q_parallel(0:Nx-1))/delta_xcb(1:Nx)   
	

         ! add the compression term (we symmetrize this in the internal region)
         ydot(2*Nx+1)        = ydot(2*Nx+1)        +        velocity(1)      * (y(2*Nx+2)       -y(2*Nx+1)       )*energy_norm/1.5d+0/delta_x(1) * switch_energy_compr
         ydot(2*Nx+2:3*Nx-1) = ydot(2*Nx+2:3*Nx-1) + 0.5d+0*velocity(2:Nx-1) * (y(2*Nx+3:3*Nx  )-y(2*Nx+2:3*Nx-1))*energy_norm/1.5d+0/delta_x(2:Nx-1) * switch_energy_compr
         ydot(2*Nx+2:3*Nx-1) = ydot(2*Nx+2:3*Nx-1) + 0.5d+0*velocity(2:Nx-1) * (y(2*Nx+2:3*Nx-1)-y(2*Nx+1:3*Nx-2))*energy_norm/1.5d+0/delta_x(1:Nx-2) * switch_energy_compr
         ydot(3*Nx)          = ydot(3*Nx)          +        velocity(Nx)     * (y(3*Nx)         -y(3*Nx-1)       )*energy_norm/1.5d+0/delta_x(Nx-1) * switch_energy_compr
      !write(*,*) 'ydot(energy) =', ydot(2*Nx+1:3*Nx) ! -------------------------------------------------------------------------

      ! ----------------------------------------------ydot for the neutral density equation --------------------------------------
	 ydot(3*Nx+1:4*Nx) = switch_neutral_source * Source_neutral(1:Nx) ![1/ (m^3 s)]
	 if (switch_neutral_leakage .eq. 1) then
		call calculate_sol_leakage_neutral_sources(Gamma_neutral, flux_source_neutral)
		ydot(3*Nx+1:4*Nx) = ydot(3*Nx+1:4*Nx) + flux_source_neutral
         else
		ydot(3*Nx+1:4*Nx) =  ydot(3*Nx+1:4*Nx) - (Gamma_neutral(1:Nx)-Gamma_neutral(0:Nx-1))/delta_xcb(1:Nx)
	 endif 
      !write(*,*) 'ydot(neutrals) =', ydot(3*Nx+1:4*Nx) !-------------------------------------------------------------------------
	
      ! --------------------------------------------- ydot for the neutral momentum equation ------------------------------------------
	      ! Sources
	      ydot(4*Nx+1:5*Nx) = switch_neutral_momentum_source * Source_vn(1:Nx)
	      ! Flux
		 if (switch_neutral_leakage .eq. 1.0d+0) then

			 call calculate_sol_leakage_neutral_sources(Gamma_mom_neutral, flux_source_neutral_momentum)
			 ydot(4*Nx+2:5*Nx) = ydot(4*Nx+2:5*Nx) + flux_source_neutral_momentum(2:Nx)

		      ! Boundary conditions flux
		      if( L_core_SOL .gt. 0.0d+0 ) then
			     ! apply boundary condition at mid point (i.e. flux = 0) or as calculated for the target at x=0
			     ydot(4*Nx+1) =  ydot(4*Nx+1) + flux_source_neutral_momentum(1)    ! [kg/(m^2 s^2)]
		      else
			     ! apply boundary condition at the X-point, as following from the constant density n(1)
			     ! velocity at i = 0:  v(0) = v(2) - 2 (S_n(1)+density_ramp_rate) delta_xcb(1) / n(1)
			     v0n = neutral_velocity(2) - 2.0d+0 * (Source_neutral(1)+ydot(4*Nx+1)) * delta_xcb(1) / neutral(1)
			     ! momentum flux at i = 0 : Gmom0 = m n (1/4)(v(0)+v(1))**2
			     Gmom0n = mass * neutral(1) * (v0n + neutral_velocity(1))**2/4.0d+0       ! [kg/(m   s^2)]
			     ydot(4*Nx+1) = ydot(4*Nx+1) - (Gamma_mom_neutral(1)-Gmom0n)/delta_xcb(1)    ! [kg/(m^2 s^2)]
		      endif
		 else
			 ydot(4*Nx+2:5*Nx) = ydot(4*Nx+2:5*Nx) - (Gamma_mom_neutral(2:Nx)-Gamma_mom_neutral(1:Nx-1))/delta_xcb(2:Nx)
			! Boundary conditions flux
		      if( L_core_SOL .gt. 0.0d+0 ) then
			     ! apply boundary condition at mid point (i.e. flux = 0) or as calculated for the target at x=0
			     ydot(4*Nx+1) =  ydot(4*Nx+1) - (Gamma_mom_neutral(1)-Gamma_mom_neutral(0))/delta_xcb(1)    ! [kg/(m^2 s^2)]
		      else
			     ! apply boundary condition at the X-point, as following from the constant density n(1)
			     ! velocity at i = 0:  v(0) = v(2) - 2 (S_n(1)+density_ramp_rate) delta_xcb(1) / n(1)
			     v0n = neutral_velocity(2) - 2.0d+0 * (Source_neutral(1)+ydot(4*Nx+1)) * delta_xcb(1) / neutral(1)
			     ! momentum flux at i = 0 : Gmom0 = m n (1/4)(v(0)+v(1))**2
			     Gmom0n = mass * neutral(1) * (v0n + neutral_velocity(1))**2/4.0d+0       ! [kg/(m   s^2)]
			     ydot(4*Nx+1) = ydot(4*Nx+1) - (Gamma_mom_neutral(1)-Gmom0n)/delta_xcb(1)    ! [kg/(m^2 s^2)]
		      endif
		 endif 
	      !do ix = 2, Nx
		!    ydot(4*Nx+ix) = ydot(4*Nx+ix) - (Gamma_mom_neutral(ix)-Gamma_mom_neutral(ix-1))/delta_xcb(ix)   ! ew 01-03-2021:
	      !enddo
		
		
		if(neutral_energy .ge. 0.0d+0) then
			tmp_neutral_energy = neutral_energy
		else
			tmp_neutral_energy = -1.0d+0*neutral_energy*temperature
		endif

		! add the pressure term in the internal region using central differencing: NB pressure = neE_n
		 dpndx = 0.0d+0
		 dpndx(1) 	= 		           (neutral(2) - neutral(1)) * e_charge * neutral_energy / delta_x(1)   
		 dpndx(2:Nx-1) 	= 		  0.5d+0 * (neutral(3:Nx) - neutral(2:Nx-1)) * e_charge * neutral_energy / delta_x(2:Nx-1)
		 dpndx(2:Nx-1)	= dpndx(2:Nx-1) + 0.5d+0 * (neutral(2:Nx-1) - neutral(1:Nx-2)) * e_charge * neutral_energy / delta_x(1:Nx-2)
		 dpndx(Nx) 	= 	  		   (neutral(Nx) - neutral(Nx-1)) * e_charge * neutral_energy / delta_x(Nx-1)
		 ydot(4*Nx+1:5*Nx) = ydot(4*Nx+1:5*Nx) - dpndx

		do ix = 2, Nx-1
		      ! add effect of numerical viscosity
		    ydot(4*Nx+ix) = ydot(4*Nx+ix) + viscosity*(neutral_velocity(ix+1) + neutral_velocity(ix-1)-2.0d+0*neutral_velocity(ix))
		enddo
		! boundary condition v_neutral = 0 at boundary
		!ydot(5*Nx) = 0.0d+0
		ydot(5*Nx) = ydot(5*Nx)  + viscosity*(2.0d+0 * neutral_velocity_cb(0) + neutral_velocity(2)-3.0d+0*neutral_velocity(1)) 
		 if( L_core_SOL .gt. 0.0d+0 .and. X_core_SOL .gt. 0.0d+0 ) then
		     ! apply boundary contion at the x=0 sheath for the pressure and viscosity terms
		     !ydot(4*Nx+1) = ydot(4*Nx+1) + viscosity*(2.0d+0 * neutral_velocity_cb(0) + neutral_velocity(2)-3.0d+0*neutral_velocity(1))  
		     ydot(4*Nx+1) = 0.0d+0	
		 endif

	 ! ----------------------------------------------ydot for the molecule density equation --------------------------------------
		 ydot(5*Nx+1:6*Nx) = switch_molecule_source * Source_molecule(1:Nx) ![1/ (m^3 s)]

	 if (switch_neutral_leakage .eq. 1.0d+0) then
		 call calculate_sol_leakage_neutral_sources(Gamma_molecule, flux_source_molecule)
		 ydot(5*Nx+1:6*Nx) = ydot(5*Nx+1:6*Nx) + flux_source_molecule
         else
		 ydot(5*Nx+1:6*Nx) = ydot(5*Nx+1:6*Nx) - (Gamma_molecule(1:Nx)-Gamma_molecule(0:Nx-1))/delta_xcb(1:Nx)
	 endif 
	       !-------------------------------------------------------------------------
		
 	 ! ---------------------------------------ydot for the neutral background volumes-------------------------
		! Calculate pump rates (in particles per second)
		!write(*,*) 'rhs: D_neutral_puff =', D_neutral_puff(1:5,internal_istep_rhs)
		!write(*,*) 'rhs: D_molecule_puff =', D_molecule_puff(1:5,internal_istep_rhs)
	        nbdot(1:5)  = (D_neutral_puff(1:5,internal_istep_rhs)  + Source_extern 	- neutral_pump  )/extern_neutral_volumes       ! neutrals
		nbdot(6:10) = (D_molecule_puff(1:5,internal_istep_rhs) + Source_extern_mol - molecule_pump ) /extern_neutral_volumes ! molecules
		
	        ! note these are fixed in this routine ! these are already considered in Source_extern
		!nbdot(1:5) = nbdot(1:5) - extern2core_flux / extern_neutral_volumes 
		!nbdot(6:10) = nbdot(6:10) - extern2core_mol / extern_neutral_volumes			
	!------------------------------------------------------------------------------------------------------------------------------------------------------------
	      ! apply evolution switches
		 ydot(     1:  Nx) = evolve_density  * ydot(     1:  Nx) / density(1:Nx)
		 ydot(  Nx+1:2*Nx) = evolve_momentum * ydot(  Nx+1:2*Nx) / momentum_norm
	      	 ydot(2*Nx+1:3*Nx) = evolve_energy   * ydot(2*Nx+1:3*Nx) / energy_norm
	      	 ydot(3*Nx+1:4*Nx) = evolve_neutral  * ydot(3*Nx+1:4*Nx) / neutral(1:Nx)
		 ydot(4*Nx+1:5*Nx) = evolve_neutral_momentum * ydot(4*Nx+1:5*Nx) / momentum_norm
		 ydot(5*Nx+1:6*Nx) = evolve_molecule * ydot(5*Nx+1:6*Nx) / molecule(1:Nx)
		 ydot(6*Nx+1:6*Nx+5)  = evolve_background * evolve_neutral * nbdot(1:5)  / max(extern_neutral_density(1:5),minimum_density)
 		 ydot(6*Nx+6:6*Nx+10) = evolve_background * evolve_molecule * nbdot(6:10)  / max(extern_molecule_density(1:5),minimum_density)
		! for densities we solve in log(density) and the derivative is than solved as ndot/density

	      !write(*,*) 'time =', time, 'ydot(Nx) =', ydot(Nx), 'density(Nx) =', density(Nx), 'Source_n(Nx) =', Source_n(Nx), 
	      !write(*,*) 'ydot molecule =', ydot(5*Nx+1:6*Nx)
	      !write(*,*) 'neutral_velocity = ', neutral_velocity

		!write(*,*) 'molecule = ', molecule
		!write(*,*) 'maxval(ydot_molecule) = ', maxval(ydot(5*Nx+1:6*Nx))
		!write(*,*) 'neutral = ', neutral
		!write(*,*) 'neutral_velocity = ', neutral_velocity
		!write(*,*) 'source_molecule = ', molecule_source
		!write(*,*) 'molecule_flux = ', molecule_flux
		!write(*,*) 'extern2sol_mol = ', extern2sol_mol
		!write(*,*) 'sol2extern_mol = ', sol2extern_mol
		!write(*,*) 'tar2extern_mol = ', tar2extern_mol
		!write(*,*) 'core2sol_mol = ', core2sol_mol
		!write(*,*) 'source_neutral = ', neutral_source
		!write(*,*) 'source_n = ', Source_n
		!write(*,*) 'source_v = ', Source_v
		!write(*,*) 'source_Q = ', Source_Q
		!write(*,*) 'density = ', density
		!write(*,*) 'velocity = ', velocity
		do i = 1,5 ! somehow the neutral density becomes negative (somewhere a minus sign? )
		tmpn = max(extern_neutral_density(i),minimum_density)
		enddo
		!tmp = Dneutral(
		!if(tmpn .lt. minimum_density+1) then
		!write(*,*) 'minimum_density', minimum_density
		
		! CHECK SIGN CONVENTION EXTERN <-> SOL
		!write(*,*) 'nb(1:5) = ', extern_neutral_density(1:5)
		!write(*,*) 'mb(6:10) = ', extern_molecule_density(1:5)
		!write(*,*) 'mean(neutral(1:i_baffle(2))) = ', sum(neutral(1:i_baffle(2)))/i_baffle(2) 
		!write(*,*) 'mean(neutral(i_baffle(2):Nx)) = ', sum(neutral(i_baffle(2):Nx))/(Nx-i_baffle(2)) 
		!write(*,*) 'sol2extern_flux(3:5) = ', sol2extern_flux(3:5)
		!write(*,*) 'sum(extern2sol_flux) = ', sum(extern2sol_flux)
		!tmp(1) = sum(extern2sol_flux) + sum(sol2extern_flux) ! sum should equate to zero
		!write(*,*) 'sum(e2s) - sum(s2e) = ', tmp(1)

		! CHECK SUMMATION EXTERN nbdot
		!write(*,*) '-neutral_pump(3:5) = ', neutral_pump(3:5)
		!write(*,*) '+neutral_puff(3:5) = ', D_neutral_puff(3:5,internal_istep_rhs)
		!write(*,*) '+SourceExt(1:5)    = ', Source_extern(1:5)
	        !write(*,*) 'zero nbdot(3:5)= ', nbdot(3:5)*extern_neutral_volumes(3:5) - ( Source_extern(3:5) +  D_neutral_puff(3:5,internal_istep_rhs) - neutral_pump(3:5)) ! nbdot(3:5)
	
		! Check Source Extern sign conventions
		!write(*,*) 'sol2extern_flux(1:5) = ',sol2extern_flux(1:5)
		!write(*,*) 'tar2extern_mol(1:5) = ',tar2extern_mol(1:5)
		!write(*,*) 'extern2core_flux(1:5) =', extern2core_flux(1:5)
		! Source = sol2extern_flux +  tar2extern_flux -  extern2core_flux + redistribution
		

		!write(*,*) 'tar2extern_flux   = ', tar2extern_flux
		!write(*,*) 'neutral = neutral(Nx-50:Nx)', neutral(50:100)
		!write(*,*) 'sol2extern_flux(3:5) =', sol2extern_flux(3:5) 
		!
		!write(*,*) 'mb(5) = ', extern_molecule_density(5)
		!endif
		
		
		!write(*,*) 'y', y
		!write(*,*) 'ydot', ydot
		!write(*,*) 'ydot(1:100) = ', ydot(1:100)
		!write(*,*) 'ydot(101:200) = ', ydot(101:200)
		!write(*,*) 'ydot(201:300) = ', ydot(201:300)
		!write(*,*) 'ydot(301:400) = ', ydot(301:400)
		!write(*,*) 'ydot(401:500) = ', ydot(401:500)
		!write(*,*) 'ydot(501:600) = ', ydot(501:600)
		!write(*,*) 'ydot(601:610) = ', ydot(601:610)
      !write(*,*) 'time =', time, 'ydot(Nx) =', ydot(Nx), 'density(Nx) =', density(Nx), 'Source_n(Nx) =', Source_n(Nx), Gamma_n(Nx)
      !write(*,*) 'ydot =', ydot
      return
   end subroutine right_hand_side
 
   end subroutine step_div1d
   


   real(wp) function kappa_parallel( temperature )
   ! function to calculate the parallel heat conductivity
      implicit none
      real(wp) :: temperature
      ! use expression from Stangeby page 187 (Chapter 4.10.1)
      kappa_parallel = 2.0d+3 * max(temperature,minimum_temperature)*max(temperature,minimum_temperature)*sqrt(max(temperature,minimum_temperature))
      return
   end function kappa_parallel

    real(wp) function D_neutral( temperature, density, velocity )
   ! function to calculate the neutral particle diffusion coefficient (Ref. Nakazawa et al. 2000 PPCF 42 401 equation 2.14)
      implicit none
      real(wp) :: temperature, density, velocity
      ! the neutral particle diffusion coefficient D == n_n kT / m charge_exchange_rate sin^2theta
      !                                               =     kT / m density <sigma v>_cx sin^2theta
      if (D_new == 1) then
      	D_neutral = e_charge * max((temperature*neutral_energy)**0.5,minimum_temperature) / (mass * density * (charge_exchange(temperature) + momentum_transfer_atoms(temperature, velocity)) * sintheta**2)
      else 
      	D_neutral = e_charge * max(temperature,minimum_temperature) / (mass * density * charge_exchange(temperature) * sintheta**2)
      end if 
      return
   end function D_neutral

   real(wp) function D_molecule( temperature, density, velocity )
   ! function to calculate the neutral particle diffusion coefficient (Ref. Nakazawa et al. 2000 PPCF 42 401 equation 2.14)
      implicit none
      real(wp) :: temperature, density, velocity
      ! the neutral particle diffusion coefficient D == n_n kT / m charge_exchange_rate sin^2theta
      !                                               =     kT / m density <sigma v>_cx sin^2theta
      
      D_molecule = e_charge * max(temperature,minimum_temperature) / (2.0d+0 * mass * density * momentum_transfer_molecules(temperature, velocity) * sintheta**2)
	!write(*,*) 'D_molecule = ', D_molecule
      return
   end function D_molecule

   real(wp) function MC_limit( fm, fc, fp)
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
   
   real(wp) function minmod2( a, b )
   ! function returning 0 when any input has different sign, otherwise minimum of the slopes a, b
      implicit none
      real(wp) :: a, b
	if( a*b .le. 0.0d+0 ) then
	 minmod2 = 0d+0
	elseif(abs(a) .lt. abs(b)) then
	 minmod2 = a
	else
	 minmod2 = b 
        endif
	return
    end function minmod2
      
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

   real(wp) function sigmoid( x, mu, s)
	! Sigmoid function
	implicit none
	real(wp) :: x, mu, s
	sigmoid = 1.0d+0/(exp(-1.0d+0/s*(x - mu)) + 1)
	return
   end function sigmoid

   SUBROUTINE JAC (NEQ, T, Y, ML, MU, PD, NROWPD)
      !dummy subroutine for calculation of Jacobian (dlsode option 21 or 24)
      INTEGER  NEQ, ML, MU, NROWPD
      DOUBLE PRECISION  T, Y(*), PD(NROWPD,*)
      RETURN
   END SUBROUTINE JAC
   

   subroutine error_report_timestep(time_step_error)
        
           implicit none
           integer, parameter :: wp = KIND(1.0D0)
           integer, intent(in) ::  time_step_error
           
           if( time_step_error .ne. 2 ) then
              write(*,*) 'time_step_error = ', time_step_error
              stop 'phys_rout: fatal error in time step'
           endif
        
           return
   end subroutine error_report_timestep

end module physics_routines
