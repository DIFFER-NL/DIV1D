module grid_data
! module defining the grid along the flux tube

   use numerics_parameters, only : Nx, dxmin
   use physics_parameters, only  : L, flux_expansion, L_core_SOL, X_core_SOL, alpha_core_profile, normalization_core_profile

   implicit none

   integer, parameter, private :: wp = KIND(1.0D0)
   real( wp ), allocatable :: x(:)         ! grid cell centre measured as distance from the X-point along the flux tube [m]
   integer                 :: i_Xpoint = 0 ! index of the grid point just above or at the X-point (used in case the core-SOL boundary is included)
   real( wp ), allocatable :: xcb(:)       ! grid cell boundaries [m]: xcb(1) = 0 and xcb(Nx+1) = L
   real( wp ), allocatable :: delta_x(:)   ! step size between grid cell centres delta_x(i) = x(i+1) - x(i) [m]
   real( wp ), allocatable :: delta_xcb(:) ! grid cell size defined as delta_x(i) = xcb(i+1) - xcb(i) [m]
   real( wp ), allocatable :: B_field(:)   ! vector holding the ratio of B/B_target @ x(:)
   real( wp ), allocatable :: B_field_cb(:)! vector holding the ratio of B/B_target @ xcb(:)
   real( wp ), private, allocatable :: xnorm(:)     ! a normalized array running from 0 to 1 at the cell boundaries, used to calculate the non-uniform grid efficiently
   integer, private        :: i            ! array index

contains

   subroutine initialize_grid
      implicit none
      if( dxmin .ge. 1.0d+0 ) then
         call uniform_grid
      else
         if( X_core_SOL .eq. 0.0 ) then
            call non_uniform_grid
         else
            call non_uniform_grid2
         endif
      endif
      ! obtain the index of the grid point just above or at the X-point
      ! and the factor normalizing the heat and particle loss profiles from the core (including the division by L_core_SOL)
      if( L_core_SOL .gt. 0.0 ) then
          normalization_core_profile = 0.0
          do i = 1, Nx
              if( x(i) .le. L_core_SOL ) then
                  i_Xpoint = i
                  normalization_core_profile = normalization_core_profile + (1 - (x(i)/L_core_SOL)**2)**alpha_core_profile * delta_x(i) 
              endif
          enddo
      endif
      call magnetic_field
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
      write(*,*) 'setting-up nonuniform grip option 2'
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
      delta_x(Nx) = 2.0d+0*(L-x(Nx))
      return
   end subroutine non_uniform_grid
   
   subroutine non_uniform_grid2
      implicit none
      integer    :: i
      allocate( x(Nx), xcb(Nx+1), delta_x(Nx), delta_xcb(Nx), xnorm(Nx+1) )
      ! set-up non-equidistant grid with two targets: grid is symmetric
      ! first define a normalized array running from 0 to 1 at the cell boundaries
      ! note that Nx must be even for this to work correctly!
      write(*,*) 'setting-up nonuniform grip option 2'
      do i = 1, Nx/2+1
         xnorm(i) = dfloat(i-1)/dfloat(Nx/2)
      end do
      ! calculate the cell boundaries on the right part of the grid grid
      xcb(Nx/2+1:Nx+1) = (L/2.0d+0) + (L/2.0d+0) * ( (2.0d+0-dxmin)*xnorm(1:Nx/2+1) - (1.0d+0-dxmin)*xnorm(1:Nx/2+1)*xnorm(1:Nx/2+1) )
      ! calculate the cell boundaries on the left part of the grid grid
      xcb(1:Nx/2+1)    = (L/2.0d+0) - (L/2.0d+0) * ( (2.0d+0-dxmin)*xnorm(Nx/2+1:1:-1) - (1.0d+0-dxmin)*xnorm(Nx/2+1:1:-1)*xnorm(Nx/2+1:1:-1) )
      ! calculate the grid cell centres
      x = (xcb(1:Nx) + xcb(2:Nx+1)) / 2.0d+0
      ! calcuate the grid cell widths
      delta_xcb = xcb(2:Nx+1) - xcb(1:Nx)
      ! calculate the step size between grid cell centres
      delta_x(1:Nx-1) = x(2:Nx) - x(1:Nx)
      ! extrapolate for step size into the target
      delta_x(Nx) = 2.0d+0*(L-x(Nx))
      return
   end subroutine non_uniform_grid2

   subroutine magnetic_field
      implicit none
      allocate( B_field(Nx), B_field_cb(Nx+1) )
      ! set the magnetic field values along the grid
      B_field    = 1.0d0
      B_field_cb = 1.0d0
      if( flux_expansion .gt. 1.0d0 ) then
          ! B_field    = 1.0d0 / ( 1.0d0 + (flux_expansion - 1.0d0) * x   / L )
          ! B_field_cb = 1.0d0 / ( 1.0d0 + (flux_expansion - 1.0d0) * xcb / L )
          ! we apply the flux expansion only along the divertor-SOL only
          B_field( i_Xpoint+1 : Nx )    = 1.0d0 / ( 1.0d0 + (flux_expansion - 1.0d0) * (x( i_Xpoint+1 : Nx ) - L_core_SOL)   / (L - L_core_SOL) )
          B_field_cb( i_Xpoint+1 : Nx ) = 1.0d0 / ( 1.0d0 + (flux_expansion - 1.0d0) * (xcb( i_Xpoint+1 : Nx ) - L_core_SOL) / (L - L_core_SOL) )
      endif
   end subroutine magnetic_field

end module grid_data
