module numerics_parameters
! module defining the numerics parameters and their default values

   implicit none

   integer, parameter, private :: wp = KIND(1.0D0)
   integer     :: Nx      = 1000     ! number of grid points along the flux tube
   integer     :: ntime   = 1000     ! number of time steps
   integer     :: nout    = 100      ! output is written to file every nout time steps
   real( wp )  :: delta_t = 1.0d-6   ! time step size [s]
   real( wp )  :: abstol  = 1.0d-4   ! required absolute error in integration (used to multiply with initial condition to set abserr_vector)
   real( wp )  :: reltol  = 1.0d-4   ! required relative error in integration

contains

   subroutine read_numerics_parameters(error)
      implicit none
      integer :: error
      namelist /div1d_numerics/ Nx, ntime, nout, delta_t, abstol, reltol
      error = 0
      read(*, div1d_numerics, IOSTAT = error)
      return
   end subroutine read_numerics_parameters
   
end module numerics_parameters
