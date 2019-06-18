module physics_parameters
! module defining the physics parameters and their default values

   implicit none

   integer, parameter, private :: wp = KIND(1.0D0)
   real( wp ) :: gamma                  = 6.5d+0      ! sheath heat transmission factor [-]
   real( wp ) :: L                      = 5.0d+1      ! lenght along flux tube from X-point to target/sheath [m] (value from Stangeby problem 5.1)
   real( wp ) :: sintheta               = 0.1d+0      ! sinus of angle theta between B-field and divertor target plate [-]
   real( wp ) :: mass                   = 3.3436d-27  ! mass of the dominant ion species (default value representing Deuterium) [kg]
   real( wp ) :: Gamma_X                = 1.0d+23     ! particle flux entering the flux tube at the X-point [/m^2s]
   real( wp ) :: q_parX                 = 1.0d+8      ! parallel heat flux entering the flux tube at the X-point [W/m^2] (value from Stangeby problem 5.1)
   real( wp ) :: initial_n              = 1.0d+20     ! initial plasma particle density (homogeneous) [/m^3]
   real( wp ) :: initial_v              = 0.0d+0      ! initial plasma velocity (homogeneous) [m/s]
   real( wp ) :: initial_T              = 1.0d+2      ! initial plasma temperature (homogeneous) [eV]
   real( wp ) :: initial_a              = 0.0d+4      ! initial neutral density (homogeneous) [/m^3]
   real( wp ) :: energy_loss_ion        = 3.0d+1      ! average loss of plasma energy due to ionization [eV]
   real( wp ) :: recycling              = 0.95d+0     ! fraction of recycled neutrals coming from the target [-]
   real( wp ) :: redistributed_fraction = 0.8d+0      ! fraction of recycled neutrals that is evenly redistributed along the SOL [-]
   real( wp ) :: neutral_residence_time = 1.0d+20     ! time scale on which neutrals are lost from the SOL [s]
   real( wp ) :: minimum_density        = 1.0d+4      ! densities are not allowed to become smaller than this value [/m^3]
   real( wp ) :: minimum_temperature    = 1.0d+0      ! the temperature is not allowed to drop below this value [eV]
   logical    :: case_AMJUEL            = .true.      ! use collision rates from AMJUEL data base

contains

   subroutine read_physics_parameters( error )
      implicit none
      integer :: error
      namelist /div1d_physics/ gamma, L, sintheta, mass, Gamma_X, q_parX, initial_n, initial_v, initial_T, initial_a, energy_loss_ion, neutral_residence_time, redistributed_fraction, recycling, case_AMJUEL
      error = 0
      read(*, div1d_physics, IOSTAT = error)
      write(*,*) 'physics read error =', error
      return
   end subroutine read_physics_parameters
   
end module physics_parameters
