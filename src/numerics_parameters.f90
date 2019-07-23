module numerics_parameters
! module defining the numerics parameters and their default values

   implicit none

   integer, parameter, private :: wp = KIND(1.0D0)
   integer     :: Nx      = 1000     ! number of grid points along the flux tube
   integer     :: ntime   = 1000     ! number of time steps
   integer     :: nout    = 100      ! output is written to file every nout time steps
   integer     :: method  = 0        ! selection of integrator ( 0 = rk4,  >0 = dvode with method_flag = method)
   integer     :: max_step= 100000   ! maximum number of internal steps in dvode
   integer     :: nzswag  = 200000   ! estimated number of nonzero elements in Jacobian (dvode)
   integer     :: evolve_density  = 1   ! evaluate density evolution 1 = yes, 0 = no (multiplier of ydot(1:Nx))
   integer     :: evolve_momentum = 1   ! evaluate momentum evolution 1 = yes, 0 = no (multiplier of ydot(1*Nx+1:2*Nx))
   integer     :: evolve_energy   = 1   ! evaluate energy evolution 1 = yes, 0 = no (multiplier of ydot(2*Nx+1:3*Nx))
   integer     :: evolve_neutral  = 1   ! evaluate neutral density evolution 1 = yes, 0 = no (multiplier of ydot(3*Nx+1:4*Nx))
   real( wp )  :: density_norm            = 0.0d+0   ! normalization of demsities    (when = 0 initial_n is used) only used in normalization of solution vector y
   real( wp )  :: temperature_norm        = 0.0d+0   ! normalization of temperatures (when = 0 1 eV is used) only used to normalize solution vector y
   real( wp )  :: velocity_norm           = 0.0d+0   ! normalization of velocities   (when = 0 sound speed at 1 eV is used) only used in to normalize solution vector y
   real( wp )  :: momentum_norm           = 0.0d+0   ! normalization of momentum (2nd part of y) = mass density_norm velocity_norm
   real( wp )  :: energy_norm             = 0.0d+0   ! normalization of energies (3rd part of y) = density_norm e_charge temperature_norm
   real( wp )  :: switch_density_source   = 1.0d+0   ! multiplier of plasma particle source term
   real( wp )  :: switch_momentum_source  = 1.0d+0   ! multiplier of momentum source term
   real( wp )  :: switch_energy_source    = 1.0d+0   ! multiplier of energy particle source term
   real( wp )  :: switch_neutral_source   = 1.0d+0   ! multiplier of neutral particle source term
   real( wp )  :: switch_charge_exchange  = 1.0d+0   ! multiplier of charge exchange rate
   real( wp )  :: switch_recombination    = 1.0d+0   ! multiplier of recombination rate
   real( wp )  :: switch_ionization       = 1.0d+0   ! multiplier of ionization rate
   real( wp )  :: switch_excitation       = 1.0d+0   ! multiplier of excitation rate
   real( wp )  :: switch_convective_heat   = 1.0d+0   ! multiplier of convective heat transport
   real( wp )  :: switch_impurity_radiation = 1.0d+0   ! multiplier of impurity radiation rate
   real( wp )  :: dxmin   = 1.0d+0   ! grid cell width at the target relative to average cell width: = 1.0 (or larger) for a uniform grid
   real( wp )  :: delta_t = 1.0d-6   ! time step size [s]
   real( wp )  :: abstol  = 1.0d-4   ! required absolute error in integration (used to multiply with initial condition to set abserr_vector)
   real( wp )  :: reltol  = 1.0d-4   ! required relative error in integration
   real( wp )  :: viscosity = 0.0d0  ! numerical viscosity to damp oscillations in velocity
   logical     :: simple_sol = .false.  ! switch for starting from the simple SOL solution when .true.
   logical     :: restart = .false.  ! switch for starting a continuation run when .true.

contains

   subroutine read_numerics_parameters(error)
      implicit none
      integer :: error
      namelist /div1d_numerics/ Nx, dxmin, ntime, nout, delta_t, abstol, reltol, method, evolve_density, evolve_momentum, evolve_energy, evolve_neutral, &
      &                         density_norm, temperature_norm, velocity_norm, &
      &                         switch_density_source, switch_momentum_source, switch_energy_source, switch_neutral_source, &
      &                         switch_charge_exchange, switch_recombination, switch_ionization, switch_excitation, switch_impurity_radiation, &
      &                         switch_convective_heat, viscosity, restart, simple_sol, &
      &                         max_step, nzswag
      error = 0
      read(*, div1d_numerics, IOSTAT = error)
      write(*,*) 'numerics read error =', error
      return
   end subroutine read_numerics_parameters
   
end module numerics_parameters
