module numerics_parameters
! module defining the numerics parameters and their default values

   implicit none

   integer, parameter, private :: wp = KIND(1.0D0)
   integer     :: Nx        = 1000     ! number of grid points along the flux tube
   integer     :: nout_steps = 10      ! number of nout  steps
   integer     :: nout      = 100      ! output is written to file every nout time steps
   integer     :: method  = 0        ! selection of integrator ( 0 = rk4,  >0 = dvode with method_flag = method)
   integer     :: istate_mod = 0     ! number of time steps before restart of dvode with istate = 1
   integer     :: max_step= 100000   ! maximum number of internal steps in dvode
   integer     :: max_attempts = 100 ! maximum number of restarts of dvode after failed integration
   integer     :: nzswag  = 200000   ! estimated number of nonzero elements in Jacobian (dvode)
   integer     :: evolve_density  = 1   ! evaluate density evolution 1 = yes, 0 = no (multiplier of ydot(1:Nx))
   integer     :: evolve_momentum = 1   ! evaluate momentum evolution 1 = yes, 0 = no (multiplier of ydot(1*Nx+1:2*Nx))
   integer     :: evolve_energy   = 1   ! evaluate energy evolution 1 = yes, 0 = no (multiplier of ydot(2*Nx+1:3*Nx))
   integer     :: evolve_neutral  = 1   ! evaluate neutral density evolution 1 = yes, 0 = no (multiplier of ydot(3*Nx+1:4*Nx))
   integer, allocatable :: IAUSER(:)    ! array specifying nonzero elements of jacobian (used by dvode)
   integer     :: NIAUSER = 0           ! dimension of array IAUSER (must be set to Number of odes + 1)
   integer, allocatable :: JAUSER(:)    ! array specifying nonzero elements of the jacobian (used by dvode)
   integer     :: NJAUSER = 0           ! dimension of array JAUSER (must be equal to Number of nonzeros in Jacobian)
   real( wp )  :: density_norm            = 0.0d+0   ! normalization of densities    (when = 0 initial_n is used) only used in normalization of solution vector y
   real( wp )  :: temperature_norm        = 0.0d+0   ! normalization of temperatures (when = 0 1 eV is used) only used to normalize solution vector y
   real( wp )  :: velocity_norm           = 0.0d+0   ! normalization of velocities   (when = 0 sound speed at 1 eV is used) only used in to normalize solution vector y
   real( wp )  :: neutral_norm            = 0.0d+0   ! normalization of neutral densities ( when renormalize = .false. equal to density_norm
   real( wp )  :: momentum_norm           = 0.0d+0   ! normalization of momentum (2nd part of y) = mass density_norm velocity_norm
   real( wp )  :: energy_norm             = 0.0d+0   ! normalization of energies (3rd part of y) = density_norm e_charge temperature_norm
   real( wp )  :: switch_density_source   = 1.0d+0   ! multiplier of plasma particle source term
   real( wp )  :: switch_momentum_source  = 1.0d+0   ! multiplier of momentum source term
   real( wp )  :: switch_energy_source    = 1.0d+0   ! multiplier of energy particle source term
   real( wp )  :: switch_neutral_source   = 1.0d+0   ! multiplier of neutral particle source term
   real( wp )  :: switch_charge_exchange  = 1.0d+0   ! multiplier of charge exchange rate
   real( wp )  :: switch_recombination    = 1.0d+0   ! multiplier of recombination rate
   real( wp )  :: switch_recombenergy     = 1.0d+0   ! multiplier of effective electron colling rate from recombination
   real( wp )  :: switch_ionization       = 1.0d+0   ! multiplier of ionization rate
   real( wp )  :: switch_excitation       = 1.0d+0   ! multiplier of excitation rate
   real( wp )  :: switch_convective_heat   = 1.0d+0   ! multiplier of convective heat transport
   real( wp )  :: switch_impurity_radiation = 1.0d+0   ! multiplier of impurity radiation rate
   real( wp )  :: dxmin   = 1.0d+0   ! grid cell width at the target relative to average cell width: = 1.0 (or larger) for a uniform grid
   real( wp )  :: delta_t = 1.0d-6   ! time step size [s]
   real( wp )  :: abstol  = 1.0d-4   ! required absolute error in integration (used to multiply with initial condition to set abserr_vector)
   real( wp )  :: reltol  = 1.0d-4   ! required relative error in integration
   real( wp )  :: viscosity = 0.0d0  ! numerical viscosity to damp oscillations in velocity
   real( wp )  :: central_differencing = 0.5d0  ! fraction of central differencing contribution in pressure gradient term
   logical     :: filter_sources = .false.  ! switch for spike filtering of source profiles when .true.
   logical     :: simple_sol = .false.  ! switch for starting from the simple SOL solution when .true.
   logical     :: restart = .false.  ! switch for starting a continuation run when .true.
   logical     :: renormalize = .false. ! redefine normalization every nout runs (might speed up ramps through multiple regimes) 

   integer :: ntime = 10 * 100

contains

   subroutine read_numerics_parameters(error)
      implicit none
      integer :: error
      integer :: ntime
      namelist /div1d_numerics/ Nx, dxmin, nout_steps, nout, delta_t, abstol, reltol, method, evolve_density, evolve_momentum, evolve_energy, evolve_neutral, &
      &                         density_norm, temperature_norm, velocity_norm, &
      &                         switch_density_source, switch_momentum_source, switch_energy_source, switch_neutral_source, &
      &                         switch_charge_exchange, switch_recombination, switch_ionization, switch_excitation, switch_impurity_radiation, &
      &                         switch_convective_heat, viscosity, central_differencing, restart, simple_sol, filter_sources, &
      &                         max_step, max_attempts, nzswag, istate_mod
      error = 0

      ntime = nout*nout_steps
      if( istate_mod .eq. 0 ) istate_mod = ntime
      read(*, div1d_numerics, IOSTAT = error)
      write(*,*) 'numerics read error =', error
      return
   end subroutine read_numerics_parameters

   subroutine extern_read_numerics_parameters(floatinnum, intinnum, loginnum)
   implicit none
   real(wp), INTENT(IN) :: floatinnum(24)
   integer, INTENT(IN) :: intinnum(12)
   logical, INTENT(IN) :: loginnum(4)
   
   ! integers
   Nx      = intinnum(1)    ! number of grid points along the flux tube
   ! ntime   = 1000     ! number of time steps is not part of matlab routine
   nout    = intinnum(3)     
   method  = intinnum(4)    !default 226 or 227, look into dvode to see the routine  
   istate_mod = intinnum(5)   ! number of time steps before restart of dvode with istate = 1
   max_step   = intinnum(6)   ! maximum number of internal steps in dvode
   max_attempts = intinnum(7) ! maximum number of restarts of dvode after failed integration
   nzswag  = intinnum(8)      ! estimated number of nonzero elements in Jacobian (dvode)
   evolve_density  = intinnum(9)  
   evolve_momentum = intinnum(10)   
   evolve_energy   = intinnum(11)   
   evolve_neutral  = intinnum(12)  
   !NIAUSER = intinnum(13)          ! dimension of array IAUSER (must be set to Number of odes + 1)
   !NJAUSER = intinnum(14)          ! dimension of array JAUSER (must be equal to Number of nonzeros in Jacobian)

   ! real(wp) 
   density_norm      = floatinnum(1)   ! normalization of densities    (when = 0 initial_n is used) only used in normalization of solution vector y
   temperature_norm  = floatinnum(2)   ! normalization of temperatures (when = 0 1 eV is used) only used to normalize solution vector y
   velocity_norm     = floatinnum(3)   ! normalization of velocities   (when = 0 sound speed at 1 eV is used) only used in to normalize solution vector y
   neutral_norm      = floatinnum(4)   ! normalization of neutral densities ( when renormalize = .false. equal to density_norm
   momentum_norm     = floatinnum(5)   ! normalization of momentum (2nd part of y) = mass density_norm velocity_norm
   energy_norm       = floatinnum(6)   ! normalization of energies (3rd part of y) = density_norm e_charge temperature_norm
   switch_density_source   = floatinnum(8)  ! multiplier of plasma particle source term
   switch_momentum_source  = floatinnum(9)    ! multiplier of momentum source term
   switch_energy_source    = floatinnum(10)   ! multiplier of energy particle source term
   switch_neutral_source   = floatinnum(11)   ! multiplier of neutral particle source term
   switch_charge_exchange  = floatinnum(12)   ! multiplier of charge exchange rate
   switch_recombination    = floatinnum(13)  ! multiplier of recombination rate
   switch_recombenergy     = floatinnum(14) ! multiplier of effective electron colling rate from recombination
   switch_ionization       = floatinnum(15)    ! multiplier of ionization rate
   switch_excitation       = floatinnum(16)   ! multiplier of excitation rate
   switch_convective_heat   = floatinnum(17)   ! multiplier of convective heat transport
   switch_impurity_radiation = floatinnum(18)   ! multiplier of impurity radiation rate
   dxmin   =    floatinnum(19)  ! grid cell width at the target relative to average cell width: = 1.0 (or larger) for a uniform grid
   delta_t =    floatinnum(20)  ! time step size [s]
   abstol  =    floatinnum(21)  ! required absolute error in integration (used to multiply with initial condition to set abserr_vector)
   reltol  =    floatinnum(22)  ! required relative error in integration
   viscosity =  floatinnum(23)  ! numerical viscosity to damp oscillations in velocity
   central_differencing = floatinnum(24)  !fraction of central differencing contribution in pressure gradient term

   ! logical
   filter_sources = loginnum(1)  ! switch for spike filtering of source profiles when .true.
   simple_sol     = loginnum(2)  ! switch for starting from the simple SOL solution when .true.
   restart        = loginnum(3)  ! switch for starting a continuation run when .true.
   renormalize    = loginnum(4)  ! redefine normalization every nout runs (might speed up ramps through multiple regimes) 

   return
   end subroutine extern_read_numerics_parameters


   
   subroutine set_jacobian_sparsity_structure
      implicit none
      integer :: icolumn, jrow, nnonzero

      ! set the dimension of IAUSER to 4*Nx + 1 (the number of equations + 1)
      NIAUSER = 4 * Nx + 1
      ! allocate the array IAUSER
      allocate( IAUSER( NIAUSER ) )

      ! set the total number of nonzero elements in Jacobian = 48*Nx - 32
      NJAUSER = 48*Nx - 32
      ! allocate the array JAUSER
      allocate( JAUSER( NJAUSER ) )
      
      ! set the array IAUSER starting from IAUSER(1) = 1
      IAUSER( 1 ) = 1
      ! for the first 4 columns add each time 8 nonzero elements
      do icolumn = 1, 4
         IAUSER( icolumn + 1 ) = IAUSER( icolumn ) + 8
      enddo
      ! for the next 4*(Nx-2) columns add 12 nonzero elements each
      do icolumn = 5, 4*(Nx-1) 
         IAUSER( icolumn + 1 ) = IAUSER( icolumn ) + 12
      enddo
      ! for the final 4 columns add each time 8 nonzero elements
      do icolumn = 4*Nx-3, 4*Nx
         IAUSER( icolumn + 1 ) = IAUSER( icolumn ) + 8
      enddo

      ! set the array JAUSER specifying the row of the nonzero element
      ! initialize the nonzero element counter to zero
      nnonzero = 0
      ! for the first 4 columns the first 8 rows are nonzero
      do icolumn = 1, 4
         do jrow = 1, 8
            nnonzero = nnonzero + 1
            JAUSER( nnonzero ) = jrow
         enddo
      enddo
      ! for the next 4*(Nx-2) columns locate the 12 rows with nonzero elements as ...
      do icolumn = 5, 4*(Nx-1) 
         do jrow = 4*((icolumn-1)/4) - 3, 4*((icolumn-1)/4) + 8
            nnonzero = nnonzero + 1
            JAUSER( nnonzero ) = jrow
         enddo
      enddo
      ! for the final 4 columns the final 8 rows contain nonzero elements
      do icolumn = 4*Nx-3, 4*Nx
         do jrow = 4*Nx - 7, 4*Nx
            nnonzero = nnonzero + 1
            JAUSER( nnonzero ) = jrow
         enddo
      enddo

      ! check the number of nonzero elements
      if( nnonzero .ne. NJAUSER ) then
         write(*,*) 'error in routine set_jacobian_sparsity_structure', nnonzero, NJAUSER
         stop
      endif

      return
   end subroutine set_jacobian_sparsity_structure
   
   subroutine set_diagonal_jacobian
      implicit none
      integer :: icolumn, jrow, nnonzero

      ! set the dimension of IAUSER to 4*Nx + 1 (the number of equations + 1)
      NIAUSER = 4 * Nx + 1
      ! allocate the array IAUSER
      allocate( IAUSER( NIAUSER ) )

      ! set the total number of nonzero elements in Jacobian = 4*Nx (the diagonal only)
      NJAUSER = 4*Nx
      ! allocate the array JAUSER
      allocate( JAUSER( NJAUSER ) )
      
      ! set the array IAUSER starting from IAUSER(1) = 1
      IAUSER( 1 ) = 1
      ! for the all columns add a single nonzero element each
      do icolumn = 1, 4*Nx
         IAUSER( icolumn + 1 ) = IAUSER( icolumn ) + 1
      enddo

      ! set the array JAUSER specifying the row of the nonzero element
      nnonzero = 4*Nx
      ! for the row number equals the column number for all nonzero elements (i.e. diagonal only)
      do icolumn = 1, 4*Nx
         JAUSER( icolumn ) = icolumn
      enddo

      ! check the number of nonzero elements
      if( nnonzero .ne. NJAUSER ) then
         write(*,*) 'error in routine set_jacobian_sparsity_structure', nnonzero, NJAUSER
         stop
      endif

      return
   end subroutine set_diagonal_jacobian
   
end module numerics_parameters
