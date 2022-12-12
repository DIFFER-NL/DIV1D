module grid_data
! module defining the grid along the flux tube
   use constants, only : pi
   use numerics_parameters, only : Nx, dxmin
   use physics_parameters, only  : L, flux_expansion, L_core_SOL, X_core_SOL, alpha_core_profile,  normalization_core_profile, &
   core_source_location, major_radius, sol_width_omp, i_omp, sintheta

   implicit none

   integer, parameter, private :: wp = KIND(1.0D0)
   real( wp ), allocatable :: x(:)         ! grid cell centre measured as distance from the X-point along the flux tube [m]
   integer                 :: i_Xpoint(2) = 0 ! index of the grid point just above or at the X-point (used in case the core-SOL boundary is included)
   real( wp ), allocatable :: xcb(:)       ! grid cell boundaries [m]: xcb(0) = 0 and xcb(Nx) = L (to be consistent with the flux arrays)
   real( wp ), allocatable :: delta_x(:)   ! step size between grid cell centres delta_x(i) = x(i+1) - x(i) [m]
   real( wp ), allocatable :: delta_xcb(:) ! grid cell size defined as delta_x(i) = xcb(i) - xcb(i-1) [m] i = 1:Nx
   real( wp ), allocatable :: B_field(:)   ! vector holding the ratio of B/B_target @ x(:)
   real( wp ), allocatable :: B_field_cb(:)! vector holding the ratio of B/B_target @ xcb(:) (changed to (0:Nx) to be consistent with flux arrays
   real( wp ), allocatable :: B_trans(:)   ! vector holding the transport based expansion of the main heat channel (Nx)
   real( wp ), allocatable :: B_trans_cb(:) ! vector holding the transport based expansion on the cell boundaries (0:Nx)
   real( wp ), allocatable :: Rcb(:)       ! vector holding the radial position of the cell boundaries (0:Nx)
   real( wp ), allocatable :: Rcc(:)       ! vector holding the radial position of the cell centres (1:Nx)
   real( wp ), allocatable :: A(:)   ! vector holding the surface area between cell and external neutral region (1:Nx)
   real( wp ), allocatable :: sintheta_cc(:) ! vector holding the ratio of B_poloidal/B_toroidal on cell centres (1:Nx)
   real( wp ), allocatable :: sintheta_cb(:) ! vector holding the ratio of B_poloidal/B_toroidal on boundaries (0:Nx)
   real( wp ), allocatable :: w_pol(:)	    ! vector holding the poloidal width of cells in the SOL (1:Nx)
   real( wp ), allocatable :: w_pol_cb(:)   ! vector holding poloidal width of SOL at cell boundaries (0:Nx) 
   real( wp ), allocatable :: volumes(:)    ! vector holding the volume of cells in the SOL (1:Nx)
   real( wp ), private, allocatable :: xnorm(:)  ! a normalized array running from 0 to 1 at the cell boundaries, used to calculate the non-uniform grid efficiently
   integer, private        :: i            ! array index
   real( wp ), private     :: mid_point    ! the mid point of the core_SOL (=0 when X_core_SOL=0; otherwise = X_core_SOL + 0.5*L_core_SOL)

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
          i_Xpoint(1) = Nx
          !mid_point = 0.0d+0
          !if( X_core_SOL .gt. 0.0d+0 ) mid_point = X_core_SOL + 0.5d+0*L_core_SOL
          mid_point = X_core_SOL + min(1.0d+0,max(core_source_location,0.0d+0))*L_core_SOL
          do i = 1, Nx
              if( x(i) .le. X_core_SOL+L_core_SOL ) then
                  if( x(i) .ge. X_core_SOL ) then
                      i_Xpoint(1) = min( i_Xpoint(1), i )
                      i_Xpoint(2) = i
                      if( X_core_SOL .eq. 0.0d+0 ) then
                          normalization_core_profile = normalization_core_profile + (1 - (x(i)/L_core_SOL)**2)**alpha_core_profile * delta_xcb(i)
                      else
                          normalization_core_profile = normalization_core_profile + 0.5*(1 - 4.0d+0*((x(i)-mid_point)/L_core_SOL)**2)**alpha_core_profile * delta_xcb(i) ! factor 0.5 because we need 2 times q_parX, Gamma_X
                      endif
                  endif
              endif
          enddo
      endif
      call magnetic_field
      call width_surface_volume
      !write(*,*) "finished setting up geometry"
      return
   end subroutine initialize_grid

   subroutine uniform_grid
      implicit none
      real( wp ) :: dx
      integer    :: i
      allocate( x(Nx), xcb(0:Nx), delta_x(Nx), delta_xcb(Nx) )
      ! set-up equidistant grid
      dx = L/Nx
      delta_x = dx
      delta_xcb = dx
      do i = 1, Nx
         x(i) = dx/2.0d+0 + dfloat(i-1)*dx
         xcb(i) = dfloat(i)*dx
      end do
      xcb(0) = 0
      return
   end subroutine uniform_grid

   subroutine non_uniform_grid
      implicit none
      integer    :: i
      allocate( x(Nx), xcb(0:Nx), delta_x(Nx), delta_xcb(Nx), xnorm(0:Nx) )
      ! set-up non-equidistant grid as described in SD1D manual
      ! first define a normalized array running from 0 to 1 at the cell boundaries
      write(*,*) 'setting-up nonuniform grid option 2'
      do i = 0, Nx
         xnorm(i) = dfloat(i)/dfloat(Nx)
      end do
      ! calculate the cell boundaries in the grid
      xcb = L * ( (2.0d+0-dxmin)*xnorm - (1.0d+0-dxmin)*xnorm*xnorm)
      ! calculate the grid cell centres
      x = (xcb(0:Nx-1) + xcb(1:Nx)) / 2.0d+0
      ! calcuate the grid cell widths
      delta_xcb = xcb(1:Nx) - xcb(0:Nx-1)
      ! calculate the step size between grid cell centres
      delta_x(1:Nx-1) = x(2:Nx) - x(1:Nx-1)
      ! extrapolate for step size into the target
      delta_x(Nx) = 2.0d+0*(L-x(Nx))
      return
   end subroutine non_uniform_grid
   
   subroutine non_uniform_grid2
      implicit none
      integer    :: i
      allocate( x(Nx), xcb(0:Nx), delta_x(Nx), delta_xcb(Nx), xnorm(0:Nx) )
      ! set-up non-equidistant grid with two targets: grid is symmetric
      ! first define a normalized array running from 0 to 1 at the cell boundaries
      ! note that Nx must be even for this to work correctly!
      write(*,*) 'setting-up nonuniform grip option 2'
      do i = 0, Nx/2
         xnorm(i) = dfloat(i)/dfloat(Nx/2)
      end do
      ! calculate the cell boundaries on the right part of the grid grid
      xcb(Nx/2:Nx) = (L/2.0d+0) + (L/2.0d+0) * ( (2.0d+0-dxmin)*xnorm(0:Nx/2) - (1.0d+0-dxmin)*xnorm(0:Nx/2)*xnorm(0:Nx/2) )
      ! calculate the cell boundaries on the left part of the grid grid
      xcb(0:Nx/2)  = (L/2.0d+0) - (L/2.0d+0) * ( (2.0d+0-dxmin)*xnorm(Nx/2:0:-1) - (1.0d+0-dxmin)*xnorm(Nx/2:0:-1)*xnorm(Nx/2:0:-1) )
      ! calculate the grid cell centres
      x = (xcb(0:Nx-1) + xcb(1:Nx)) / 2.0d+0
      ! calcuate the grid cell widths
      delta_xcb = xcb(1:Nx) - xcb(0:Nx-1)
      ! calculate the step size between grid cell centres
      delta_x(1:Nx-1) = x(2:Nx) - x(1:Nx)
      ! extrapolate for step size into the target
      delta_x(Nx) = 2.0d+0*(L-x(Nx))
      return
   end subroutine non_uniform_grid2

   subroutine magnetic_field
      implicit none
      allocate( B_field(Nx), B_field_cb(0:Nx), B_trans(Nx), B_trans_cb(0:Nx) )
      write(*,*) 'setting-up magnetic fields'
      ! set the magnetic field values along the grid
      ! TODO: seperate  transp_expansion & magn_expansion 
      B_field    = 1.0d0
      B_field_cb = 1.0d0
      B_trans 	 = 1.0d0
      B_trans_cb = 1.0d0
      if( flux_expansion .ne. 1.0d0 ) then ! NB. we staan ook flux contractie toe, i.e. flux_expansion < 1 
          ! B_field    = 1.0d0 / ( 1.0d0 + (flux_expansion - 1.0d0) * x   / L )
          ! B_field_cb = 1.0d0 / ( 1.0d0 + (flux_expansion - 1.0d0) * xcb / L )
          ! we apply the flux expansion only along the divertor-SOL only
          B_field( i_Xpoint(2)+1 : Nx )      = 1.0d0 / ( 1.0d0 + (flux_expansion - 1.0d0) * (x  ( i_Xpoint(2)+1 : Nx ) - X_core_SOL - L_core_SOL) / (L - X_core_SOL - L_core_SOL) )
          B_field_cb( i_Xpoint(2)+1 : Nx )   = 1.0d0 / ( 1.0d0 + (flux_expansion - 1.0d0) * (xcb( i_Xpoint(2)+1 : Nx ) - X_core_SOL - L_core_SOL) / (L - X_core_SOL - L_core_SOL) )
          if( L_core_SOL .gt. 0.0d+0 .and. X_core_SOL .gt. 0.0d+0 ) then
              B_field( 1: i_Xpoint(1)-1 )    = 1.0d0 / ( 1.0d0 + (flux_expansion - 1.0d0) * (X_core_SOL - x  ( 1: i_Xpoint(1)-1 ) ) / X_core_SOL )
              B_field_cb( 0: i_Xpoint(1)-1 ) = 1.0d0 / ( 1.0d0 + (flux_expansion - 1.0d0) * (X_core_SOL - xcb( 0: i_Xpoint(1)-1 ) ) / X_core_SOL )
          endif
      endif
   end subroutine magnetic_field


   subroutine width_surface_volume
           implicit none
           allocate( Rcc(Nx), Rcb(0:Nx), sintheta_cc(Nx), sintheta_cb(0:Nx), A(Nx),w_pol(1:Nx), w_pol_cb(0:Nx), volumes(Nx) )
	   write(*,*) 'setting up width, surface and volume'
           Rcc = major_radius
           Rcb = major_radius
           sintheta_cc = sintheta
	   sintheta_cb = sintheta
	   A = delta_xcb*sintheta_cc*Rcc*2*pi  ! A_extern = dz 2.pi.R = dx*sintheta 2.pi.R
	   ! lambda_poloidal *sintheta * B_transp = constant (as B_field drops out with R dependence) 
	   w_pol_cb = sol_width_omp*sintheta_cb(i_omp)/sintheta_cb * B_trans_cb(i_omp)/B_trans_cb 
	   w_pol = (w_pol_cb(0:Nx-1) + w_pol_cb(1:Nx)) / 2.0d+0         
           volumes = A * w_pol
   end subroutine width_surface_volume

end module grid_data
