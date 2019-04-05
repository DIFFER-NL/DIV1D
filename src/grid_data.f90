module grid_data
! module defining the grid along the flux tube

   use numerics_parameters, only : Nx
   use physics_parameters, only  : L

   implicit none

   integer, parameter, private :: wp = KIND(1.0D0)
   real( wp ), allocatable :: x(:)        ! distance from the X-point along the flux tube[m]
   real( wp ), allocatable :: delta_x(:)  ! step size in grid defined as delta_x(i) = x(i+1) - x(i) [m]

contains

   subroutine initialize_grid
      implicit none
      real( wp ) :: dx
      integer    :: i
      allocate( x(Nx), delta_x(Nx) )
      ! set-up equidistant grid
      dx = L/Nx
      delta_x = dx
      do i = 1, Nx
         x(i) = dx/2.0d+0 + (i-1)*dx
      end do
      return
   end subroutine initialize_grid

end module grid_data
