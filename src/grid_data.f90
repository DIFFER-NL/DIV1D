module grid_data
! module defining the grid along the flux tube

   use numerics_parameters, only : Nx, dxmin
   use physics_parameters, only  : L

   implicit none

   integer, parameter, private :: wp = KIND(1.0D0)
   real( wp ), allocatable :: x(:)         ! grid cell centre measured as distance from the X-point along the flux tube [m]
   real( wp ), allocatable :: xcb(:)       ! grid cell boundaries [m]: xcb(1) = 0 and xcb(Nx+1) = L
   real( wp ), allocatable :: delta_x(:)   ! step size between grid cell centres delta_x(i) = x(i+1) - x(i) [m]
   real( wp ), allocatable :: delta_xcb(:) ! grid cell size defined as delta_x(i) = xcb(i+1) - xcb(i) [m]
   real( wp ), private, allocatable :: xnorm(:)     ! a normalized array running from 0 to 1 at the cell boundaries, used to calculate the non-uniform grid efficiently

contains

   subroutine initialize_grid
      implicit none
      if( dxmin .ge. 1.0d+0 ) then
         call uniform_grid
      else
         call non_uniform_grid
      endif
      return
   end subroutine initialize_grid

   subroutine uniform_grid
      implicit none
      real( wp ) :: dx
      integer    :: i
      allocate( x(Nx), xcb(Nx+1), delta_x(Nx), delta_xcb(Nx) )
      ! set-up equidistant grid
      dx = L/Nx
      delta_x = dx
      delta_xcb = dx
      do i = 1, Nx
         x(i) = dx/2.0d+0 + (i-1)*dx
         xcb(i) = (i-1)*dx
      end do
      xcb(Nx+1) = L
      return
   end subroutine uniform_grid

   subroutine non_uniform_grid
      implicit none
      integer    :: i
      allocate( x(Nx), xcb(Nx+1), delta_x(Nx), delta_xcb(Nx), xnorm(Nx+1) )
      ! set-up non-equidistant grid as described in SD1D manual
      ! first define a normalized array running from 0 to 1 at the cell boundaries
      do i = 1, Nx+1
         xnorm(i) = dfloat(i-1)/dfloat(Nx)
      end do
      ! calculate the cell boundaries in the grid
      xcb = L * ( (2.0d+0-dxmin)*xnorm - (1.0d+0-dxmin)*xnorm*xnorm)
      ! calculate the grid cell centres
      x = (xcb(1:Nx) + xcb(2:Nx+1)) / 2.0d+0
      ! calcuate the grid cell widths
      delta_xcb = xcb(2:Nx+1) - xcb(1:Nx)
      ! calculate the step size between grid cell centres
      delta_x(1:Nx-1) = x(2:Nx) - x(1:Nx)
      ! extrapolate for step size into the target
      delta_x(Nx) = 4.0d+0*(L-x(Nx))-delta_x(Nx-1)
      return
   end subroutine non_uniform_grid

end module grid_data
