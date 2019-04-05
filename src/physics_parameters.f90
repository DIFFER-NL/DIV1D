module physics_parameters
! module defining the physics parameters and their default values

   implicit none

   integer, parameter, private :: wp = KIND(1.0D0)
   real( wp ) :: gamma   = 6.5d+0      ! sheath heat transmission factor [-]
   real( wp ) :: L       = 5.0d+1      ! lenght along flux tube from X-point to target/sheath [m] (value from Stangeby problem 5.1)
   real( wp ) :: mass    = 3.3444d-27  ! mass of the dominant ion species (default value representing Deuterium) [kg]
   real( wp ) :: Gamma_X = 1.0d+23     ! particle flux entering the flux tube at the X-point [/m^2s]
   real( wp ) :: q_parX  = 1.0d+8      ! parallel heat flux entering the flux tube at the X-point [W/m^2] (value from Stangeby problem 5.1)
   real( wp ) :: initial_n  = 1.0d+20  ! parallel heat flux entering the flux tube at the X-point [/m^3]
   real( wp ) :: initial_v  = 0.0d+0   ! parallel heat flux entering the flux tube at the X-point [m/s]
   real( wp ) :: initial_T  = 1.0d+2   ! parallel heat flux entering the flux tube at the X-point [eV]

contains

   subroutine read_physics_parameters( error )
      implicit none
      integer :: error
      namelist /div1d_physics/ L, Gamma_X, q_parX, initial_n, initial_v, initial_T
      error = 0
      read(*, div1d_physics, IOSTAT = error)
      return
   end subroutine read_physics_parameters
   
end module physics_parameters
