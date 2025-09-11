
module grid_data
! module defining the grid along the flux tube
   use constants, only : pi
   use numerics_parameters, only : Nx, dxmin, wide_core_profile
   use physics_parameters, only  : L, flux_expansion, flux_expansion_left, trans_expansion, trans_expansion_left, L_core_SOL, L_baffle, X_core_SOL, & 
				   alpha_core_profile_Q, alpha_core_profile_n,&
  	  		           gas_puff_source, gas_puff_location, gas_puff_width, prf_imp_con, & 
   				   major_radius, sol_width_omp, location_omp, sintheta, case_single_null 

   implicit none

   integer, parameter, private :: wp = KIND(1.0D0)
   real( wp ), allocatable :: x(:)         ! grid cell centre measured as distance from the X-point along the flux tube [m]
   integer                 :: i_Xpoint(2) = 1 ! index of the grid point just above or at the X-point (used in case the core-SOL boundary is included)
 integer 		   :: i_baffle(2) ! index of the grid point just avove or at the Baffle point
   real( wp ) 		   :: A_wet(2) ! wetted area
   real( wp ), allocatable :: xcb(:)       ! grid cell boundaries [m]: xcb(0) = 0 and xcb(Nx) = L (to be consistent with the flux arrays)
   real( wp ), allocatable :: delta_x(:)   ! step size between grid cell centres delta_x(i) = x(i+1) - x(i) [m]
   real( wp ), allocatable :: delta_xcb(:) ! grid cell size defined as delta_x(i) = xcb(i) - xcb(i-1) [m] i = 1:Nx
   real( wp ), allocatable :: B_field(:)   ! vector holding the ratio of B/B_target @ x(:)
   real( wp ), allocatable :: B_field_cb(:)! vector holding the ratio of B/B_target @ xcb(:) (changed to (0:Nx) to be consistent with flux arrays
   real( wp ), allocatable :: B_trans(:)   ! vector holding the transport based expansion of the main heat channel (Nx)
   real( wp ), allocatable :: B_trans_cb(:) ! vector holding the transport based expansion on the cell boundaries (0:Nx)
   real( wp ), allocatable :: R_cb(:)       ! vector holding the radial position of the cell boundaries (0:Nx)
   real( wp ), allocatable :: R_cc(:)       ! vector holding the radial position of the cell centres (1:Nx)
  real( wp ), allocatable :: Area_extern(:)   ! vector holding the surface area between cell and external neutral region (1:Nx)
   real( wp ), allocatable :: sintheta_cc(:) ! vector holding the ratio of B_poloidal/B_total on cell centres (1:Nx)
   real( wp ), allocatable :: sintheta_cb(:) ! vector holding the ratio of B_poloidal/B_total on boundaries (0:Nx)
   real( wp ), allocatable :: sol_width_pol(:)	    ! vector holding the poloidal width of cells in the SOL (1:Nx)
   real( wp ), allocatable :: sol_width_pol_cb(:)   ! vector holding poloidal width of SOL at cell boundaries (0:Nx) 
   real( wp ), allocatable :: volumes(:)    ! vector holding the volume of cells in the SOL (1:Nx)
   real( wp ), allocatable :: gas_puff_profile(:) ! gas puff distribution (1:Nx)
   real( wp ), allocatable :: core_source_profile_Q(:)  ! core source distribution (Xp1:Xp2)
   real( wp ), allocatable :: core_source_profile_n(:)  ! core source distribution (Xp1:Xp2)
   real( wp ), allocatable :: I_core_source_profile_Q(:) ! core source distribution on the grid (1:Nx) 
   real( wp ), allocatable :: I_core_source_profile_n(:) ! core source distribution on the grid (1:Nx)
   real( wp ), private, allocatable :: xnorm(:)  ! a normalized array running from 0 to 1 at the cell boundaries, used to calculate the non-uniform grid efficiently
   integer, private        :: i    ,ix      ! array index
   real( wp )	 	   :: mid_point    ! the mid point of the core_SOL (=0 when X_core_SOL=0; otherwise = X_core_SOL + core_source_location*L_core_SOL)
   integer 		   :: i_omp = 0 ! index of the cell closest to the outer midplane
   real( wp )		   :: X_omp
   real( wp ), allocatable :: A_int(:) ! Areas of neutral cells that are square 
 
   ! vectors not used in computations but that define the geometry
   real( wp ), allocatable :: Z_cc(:) 	    ! vector holding possible vertical position of the flux tube (1:Nx) ( this is used to generate sintheta and can be stored in the output file )
   real( wp ), allocatable :: Z_cb(:) 	    ! vector holding possible vertical position of the flux tube (0:Nx) ( this is used to generate sintheta and can be stored in the output file )
   real( wp ), allocatable :: nr_cc(:)       ! vector direction normal to the flux tube (1:Nx)
   real( wp ), allocatable :: nz_cc(:)       ! vector direction normal to the flux tube (1:Nx)
 !  real( wp ), allocatable :: n_cb(:)       ! vector direction normal to the flux tube (0:Nx)
   real( wp ), allocatable :: vesrz(:,:) ! vector holding vessel r coordinate (max Nx-1 points)
   real( wp ), allocatable :: vesr(:),vesz(:),vesc(:) ! vector holding parts of vesrz.. suboptimal implementation
  ! real( wp ), allocatable :: vesz(:) ! vector holding vessel z coordinate (max Nx-1 points)
  ! real( wp ) :: i_vessel_cutting_points(9) ! vector holding vessel cutting points (pfr cut,target inn, target out, baffle, stag,baffle,target out,target inn, pfr cut)

contains

   subroutine initialize_grid(Nx, x, xcb, B_field, B_field_cb, B_trans, B_trans_cb, &
		            R_cc, R_cb, Area_extern, sintheta_cc, sintheta_cb, sol_width_pol, sol_width_pol_cb, volumes, &
			    gas_puff_profile, E_core_source_profile_Q, E_core_source_profile_n, i_omp, i_Xpoint, i_baffle, A_int)
      implicit none
      integer :: error
      integer, intent(in) :: Nx
      integer :: i_omp, i_Xpoint(2), i_baffle(2), i_tmp(2)
      real( wp ) :: normalization_core_profile = 0.0d+0 ! avoid possible divisions by zero
      !real( wp ) :: edge_rounding_midpoint = 0.0d+0
      real( wp ) :: X_omp = 0.0d+0
      real( wp ) :: L_baffle_real = 0.0d+0 
      real( wp ) :: mid_point = 0.0d+0
      real( wp ) :: L_tmp  = 0.0d0
      real( wp ) :: R_max
	logical :: exists
      ! note the global grid parameters of the program are given as input to this subroutine.
      ! the subroutine is allowed to modify them and then gives them back, however anything 
      ! changed but not passed through the interface is not known to the DIV1D program 
      real(wp) :: x(Nx), xcb(0:Nx), B_field(Nx), B_field_cb(0:Nx), B_trans(Nx), B_trans_cb(0:Nx), &
		            Area_extern(Nx), R_cb(0:Nx), R_cc(Nx), sintheta_cc(Nx), sintheta_cb(0:Nx), &
			    sol_width_pol(Nx), sol_width_pol_cb(0:Nx), volumes(Nx), &
			    gas_puff_profile(Nx), E_core_source_profile_Q(Nx),E_core_source_profile_n(Nx), A_int(1:Nx) 
      !real( wp ), allocatable :: core_source_profile(:)
      write(*,*) 'point 1'
      ! define grid
      call define_grid(delta_x, delta_xcb, x, xcb, xnorm)      
      !Find X_points and OMP
      call locate_X_points(X_omp, i_omp, i_Xpoint, i_baffle, L_baffle, L_baffle_real, x, xcb, X_core_SOL, L_core_SOL)     
            write(*,*) 'point 2'
      !write(*,*) 'Xpoints: ', i_Xpoint
      call magnetic_field(Nx, i_Xpoint, x_core_sol, L_core_SOL, L, flux_expansion, flux_expansion_left, trans_expansion, &
                          trans_expansion_left, B_field, B_field_cb, B_trans, B_trans_cb, x, xcb)
      call width_surface_volume(B_trans_cb, xcb, A_int, R_cc, R_cb, sintheta_cc, sintheta_cb, Area_extern, &
                                delta_xcb, sol_width_pol, sol_width_pol_cb, volumes, A_wet, B_field, B_field_cb, &
                                major_radius, sintheta, sol_width_omp, vesr, vesz, vesc, Z_cc, Z_cb, nr_cc, nz_cc)
      call initialize_gas_puff(gas_puff_profile, gas_puff_location, gas_puff_width, delta_xcb)
      ! this is the only profile that is initialized here as the size is unknown up to this point	
      !allocate(core_source_profile(i_Xpoint(2)-i_Xpoint(1)+1))
     
      !define core_source_profile if there is a core SOL
	!write(*,*) 'l96'
	!write(*,*) 'L_core_SOL, X_core_SOL', L_core_SOL, X_core_SOL
	
	call define_core_source_profiles(wide_core_profile, i_tmp, i_xpoint, L_tmp, i_baffle, L_baffle_real, I_core_source_profile_N, I_core_source_profile_Q)

	call normalize_profile(L_core_sol, X_core_Sol, I_core_source_profile_Q, I_core_source_profile_n, x, i_Xpoint, alpha_core_profile_Q, alpha_core_profile_n, &
	                       i_tmp, L_tmp, mid_point, R_max, R_cc, volumes, normalization_core_profile, core_source_profile_Q, core_source_profile_n, &
	                       E_core_source_profile_Q, E_core_source_profile_n)
	! check if everything is well defined
      do ix = 1,Nx
 	if (isnan(E_core_source_profile_Q(ix))) stop 'grid setup: "E_core_source_profile_Q" is a NaN'
	if (isnan(E_core_source_profile_n(ix))) stop 'grid setup: "E_core_source_profile_n" is a NaN'
	if (isnan(core_source_profile_Q(ix))) stop 'grid setup: "core_source_profile_Q" is a NaN'
	if (isnan(core_source_profile_n(ix))) stop 'grid setup: "core_source_profile_n" is a NaN'
	if (isnan(x(ix))) stop 'grid: "x" is a NaN'
	if (isnan(B_field(ix))) stop 'grid setup: "B_field" is a NaN'
	if (isnan(B_trans(ix))) stop 'grid setup: "B_trans" is a NaN'
	if (isnan(R_cc(ix))) stop 'grid setup: "R_cc" is a NaN'
	if (isnan(Area_extern(ix))) stop 'grid setup: "Area_extern" is a NaN'
	if (isnan(sintheta_cc(ix))) stop 'grid setup: "sintheta_cc" is a NaN'
	if (isnan(sol_width_pol(ix))) stop 'grid setup: "sol_width_pol" is a NaN'
	if (isnan(gas_puff_profile(ix))) stop 'grid setup: "gas_puff_profile" is a NaN'
	if (isnan(volumes(ix))) stop 'grid setup: "volumes" is a NaN'
      enddo	

      write(*,*) "F finished setting up geometry"
       namelist /div1d_grid/ i_omp, i_Xpoint, i_baffle, mid_point, X_omp, x, xcb, B_field, B_field_cb, B_trans, B_trans_cb, &
		            Area_extern, R_cc, R_cb, Z_cc, Z_cb, nr_cc, nz_cc, vesrz, sintheta_cc, sintheta_cb, sol_width_pol, sol_width_pol_cb, volumes, &
			    gas_puff_profile, E_core_source_profile_Q, E_core_source_profile_n, prf_imp_con, A_int, A_wet
	
	inquire(file="div1d_output.txt", exist=exists)
	if(exists) then
      	open( UNIT=10, FILE='div1d_output.txt', status='old', position="append")
	else
	open( UNIT=10, FILE='div1d_output.txt', status='new')
	endif
	write(10,  div1d_grid )
	close(10)
	! namelists read routine in matlab does not consistently read arrays
      return
   end subroutine initialize_grid

   subroutine define_grid(delta_x, delta_xcb, x, xcb, xnorm)
    real(wp) :: x(:), delta_x(:), delta_xcb(:), xcb(0:) 
    real(wp), allocatable :: xnorm(:)
      if( dxmin .ge. 1.0d+0 ) then
         call uniform_grid(delta_x, delta_xcb, x, xcb)
      else
         if( X_core_SOL .eq. 0.0 ) then
            !only calculates X, returns also X
            call non_uniform_grid(delta_x, delta_xcb, x, xcb, xnorm)
         else
            call non_uniform_grid2(delta_x, delta_xcb, x, xcb, xnorm)
         endif
      endif
   end subroutine define_grid

    subroutine locate_X_points(X_omp, i_omp, i_Xpoint, i_baffle, L_baffle, L_baffle_real, x, xcb, X_core_SOL, L_core_SOL)
        integer :: i_omp, i_Xpoint(:), i_baffle(:)
        real(wp) :: x_omp, L_baffle, L_baffle_real,  X_core_sol, L_core_sol
        real(wp) :: x(:), xcb(0:Nx)            
            
          X_omp = X_core_SOL + min(1.0d+0,max(location_omp,0.0d+0))*L_core_SOL
          ! find X-points and OMP
          i_Xpoint = Nx ! make sure it is not out of bounds later
          i_baffle = Nx ! idem dito
          do i = 1, Nx
            if( x(i) .le. X_core_SOL+L_core_SOL ) then
                if( x(i) .ge. X_core_SOL ) then
                    i_Xpoint(1) = min( i_Xpoint(1), i )
                    i_Xpoint(2) = i
	                if( x(i) .le. X_omp ) then
	                    i_omp = i	
	                endif
                endif
	        endif
	        if( x(i) .le. X_core_SOL+L_core_SOL+L_baffle ) then
              if( x(i) .ge. X_core_SOL ) then
                i_baffle(1) = min( i_baffle(1), i )
                i_baffle(2) = i
	          endif
	        endif
          enddo
          X_omp = xcb(i_omp)
          L_baffle_real = xcb(i_baffle(2)) - X_core_SOL - L_core_SOL
    end subroutine locate_X_points

    subroutine define_core_source_profiles(wide_core_profile, i_tmp, i_xpoint, L_tmp, i_baffle, L_baffle_real, I_core_source_profile_N, I_core_source_profile_Q)
    integer :: wide_core_profile 
    integer :: i_tmp(:), i_xpoint(:), i_baffle(:)
    real(wp) :: L_tmp, L_baffle_real 
    real(wp), allocatable :: I_core_source_profile_Q(:), I_core_source_profile_n(:) 
    
    if(wide_core_profile == 0) then
        i_tmp = i_Xpoint
        L_tmp = 0.0d+0
    else
        i_tmp = i_baffle
        L_tmp = L_baffle_real
    endif
    
    write(*,*) 'point 3b'
    
    if( allocated(I_core_source_profile_n) .eqv. .false.) allocate(I_core_source_profile_n(i_tmp(2)-i_tmp(1)+1))
    if( allocated(I_core_source_profile_Q) .eqv. .false.) allocate(I_core_source_profile_Q(i_Xpoint(2)-i_Xpoint(1)+1))
    
    write(*,*) 'point 4'
    
    end subroutine define_core_source_profiles
    
   subroutine cc2cb(Nx, v,v_cb) !,delta_x,delta_xcb )
   ! interpolates cell boundary (cb) values from cell centre (cc)
   ! extrapolates outer cell boundaries
	implicit none
	real( wp ), intent(in)  :: v(Nx)
	!real( wp ), intent(in)  :: delta_x(Nx-1)
	!real( wp ), intent(in)  :: delta_xcb(Nx)
    real( wp ), intent(out) :: v_cb(0:Nx)
    integer, intent(in) :: Nx 
	! interpolate 0.5d+0*
	v_cb(1:Nx-1) = ( v(1:Nx-1)*delta_xcb(2:Nx) + v(2:Nx)*delta_xcb(1:Nx-1) ) / (2.0d+0*delta_x(1:Nx-1))
        ! extrapolate	
	v_cb(0)  = v(1)  -  (v(2)  - v(1)   )/delta_x(1)  * 0.5d+0*delta_xcb(1) 
	v_cb(Nx) = v(Nx) +  (v(Nx) - v(Nx-1))/delta_x(Nx-1) * 0.5d+0*delta_xcb(Nx) 
   end subroutine cc2cb 

    subroutine normalize_profile(L_core_sol, X_core_Sol, I_core_source_profile_Q, I_core_source_profile_n, x, i_Xpoint, alpha_core_profile_Q, alpha_core_profile_n, &
	                       i_tmp, L_tmp, mid_point, R_max, R_cc, volumes, normalization_core_profile, core_source_profile_Q, core_source_profile_n, &
	                       E_core_source_profile_Q, E_core_source_profile_n)
    integer :: i_tmp(:), i_XPoint(:)
    real(wp) :: L_tmp, L_core_sol, X_core_sol, mid_point, R_max, normalization_core_profile, alpha_core_profile_Q, alpha_core_profile_n
    real(wp) :: x(:), R_cc(:), volumes(:), I_core_source_profile_Q(:), I_core_source_profile_n(:), core_source_profile_Q(:), core_source_profile_n(:), & 
                E_core_source_profile_Q(:), E_core_source_profile_n(:)
    
    	!subroutine start 
      I_core_source_profile_Q = 0.0d+0 
      I_core_source_profile_n = 0.0d+0 
	core_source_profile_n = 0.0d+0
	core_source_profile_Q = 0.0d+0
      !write(*,*) 'I_core_source_profile_n', I_core_source_profile_n
      if( L_core_SOL .gt. 0.0d+0 .and. X_core_SOL .eq. 0.0d+0 ) then
	!write(*,*) 'l99'
      	I_core_source_profile_Q = (1.0d+0 - (x(i_Xpoint(1):i_Xpoint(2))/L_core_SOL)**2)**alpha_core_profile_Q 
	I_core_source_profile_n = (1.0d+0 - (x(i_tmp(1):i_tmp(2))/(L_core_SOL+L_tmp))**2)**alpha_core_profile_n 
	!write(*,*) 'I_core_source_profile_n', I_core_source_profile_n
      elseif( L_core_SOL .gt. 0.0d+0 .and. X_core_SOL .gt. 0.0d+0 ) then ! double target!
	!write(*,*) 'l103'
	! if single null, bias to right ( for upper and lower, the left and right should be swapped )
!	if( case_single_null ) then   
		!location_omp
		!do
		!mid_point = X_
		mid_point = X_core_SOL + 0.5d+0*L_core_SOL
		R_max = 0.0d+0
		do i = i_tmp(1),i_tmp(2)
		R_max = max(R_cc(i),R_max)
		enddo
		do i = i_tmp(1),i_tmp(2)
			if (R_cc(i) .eq. R_max)	then
			mid_point = x(i);
			endif
		enddo 
	!	else then
		! if double null, center in the middle:
		!mid_point = X_core_SOL + 0.5d+0*L_core_SOL
     		I_core_source_profile_Q = (1.0d+0 - 4.0d+0*((x(i_Xpoint(1):i_Xpoint(2))-mid_point)/L_core_SOL)**2)**alpha_core_profile_Q
		I_core_source_profile_n = (1.0d+0 - 4.0d+0*((x(i_tmp(1):i_tmp(2))-mid_point)/(L_core_SOL+2.0d+0*L_tmp))**2)**alpha_core_profile_n
		! make sure no negative numbers
		I_core_source_profile_Q = max(I_core_source_profile_Q,0.0d+0)
		I_core_source_profile_n = max(I_core_source_profile_n,0.0d+0)
		! add option for negative alpha_core_power -> load profile from .dat file?
      endif

      ! normalize the profile    
      if( L_core_SOL .gt. 0.0d+0 ) then
      I_core_source_profile_Q = I_core_source_profile_Q * volumes(i_Xpoint(1):i_Xpoint(2))
      I_core_source_profile_n = I_core_source_profile_n * volumes(i_tmp(1):i_tmp(2))
      normalization_core_profile = sum(I_core_source_profile_Q) 
      I_core_source_profile_Q = I_core_source_profile_Q / normalization_core_profile
      !write(*,*) 'normalization core profile', normalization_core_profile	
      normalization_core_profile = sum(I_core_source_profile_n) 
	!write(*,*) 'normalization core profile', normalization_core_profile	
	!write(*,*) 'i_tmp', i_tmp
	!write(*,*) 'i_baffle', i_baffle
	!write(*,*) 'i_Xpoint', i_Xpoint
	!write(*,*) 'I_core_source_profile_n', I_core_source_profile_n	
      I_core_source_profile_n = I_core_source_profile_n / normalization_core_profile
      ! fill the core source profile array that matches the Xgrid with the part that is filled.
      core_source_profile_Q(i_Xpoint(1):i_Xpoint(2)) = I_core_source_profile_Q
      core_source_profile_n(i_tmp(1):i_tmp(2)) = I_core_source_profile_n
      endif
      E_core_source_profile_n = core_source_profile_n
      E_core_source_profile_Q = core_source_profile_Q
    
    end subroutine normalize_profile
    
   subroutine cb2cc(Nx, v,v_cb)
   ! interpolates cell centre values from cell boundaries.
   ! this interpolation is trivial as the cell centre is defined as the center.
	implicit none
	integer, intent(in) :: Nx 
	real( wp ), intent(in)  :: v_cb(0:Nx)
	real( wp ), intent(out) :: v(Nx)
 	v(1:Nx) = 0.5d+0*( v_cb(0:Nx-1) + v_cb(1:Nx) ) 
   end subroutine cb2cc 	
 
   subroutine uniform_grid(delta_x, delta_xcb, x, xcb)
      implicit none
      real( wp ) :: dx
      real(wp) :: x(:)
      !The lower bound MUST be declared. if left implicit, the lower bound will default
      !back to 1, instead of 0. 
      real(wp) :: delta_x(:), delta_xcb(:), xcb(0:) 
      integer    :: i
      !allocate( x(Nx), xcb(0:Nx), delta_x(Nx-1), delta_xcb(Nx) )
      ! set-up equidistant grid
      dx = L/Nx
      delta_x = dx
      delta_xcb = dx
      do i = 1, Nx
         x(i) = dx/2.0d+0 + dfloat(i-1)*dx
         xcb(i) = dfloat(i)*dx
      end do
      xcb(0) = 0.0d+0
      return
   end subroutine uniform_grid

   subroutine non_uniform_grid(delta_x, delta_xcb, x, xcb, xnorm)
      implicit none
      integer :: i
      real(wp) :: x(Nx)
      real(wp) :: delta_x(:), delta_xcb(:), xcb(0:)
      real(wp), allocatable :: xnorm(:)
      !allocate( x(Nx), xcb(0:Nx), delta_x(Nx-1), delta_xcb(Nx), xnorm(0:Nx) )
	! xnorm is the only vector initialized here as it is not used anywhere else
 	if( allocated(xnorm) .eqv. .false.) allocate(xnorm(0:Nx))
      ! set-up non-equidistant grid as described in SD1D manual
      ! first define a normalized array running from 0 to 1 at the cell boundaries
      write(*,*) 'setting-up nonuniform grid option 1'
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
      return
   end subroutine non_uniform_grid
   
   subroutine non_uniform_grid2(delta_x, delta_xcb, x, xcb, xnorm)
      implicit none
      integer    :: i
      real(wp) :: x(Nx)
      real(wp) :: delta_x(:), delta_xcb(:), xcb(0:)
      real(wp), allocatable :: xnorm(:)
      !allocate( x(Nx), xcb(0:Nx), delta_x(Nx-1), delta_xcb(Nx), xnorm(0:Nx) )
	write(*,*) 'checking allocation xnorm'
	if( allocated(xnorm) .eqv. .false.) allocate( xnorm(0:Nx) )
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
      call cb2cc(Nx, x,xcb)
      ! calcuate the grid cell widths
      delta_xcb = xcb(1:Nx) - xcb(0:Nx-1)
      ! calculate the step size between grid cell centres
      delta_x(1:Nx-1) = x(2:Nx) - x(1:Nx-1)
      return
   end subroutine non_uniform_grid2

   subroutine magnetic_field (Nx, i_Xpoint, x_core_sol, L_core_SOL, L, flux_expansion, flux_expansion_left, trans_expansion, &
                              trans_expansion_left, B_field, B_field_cb, B_trans, B_trans_cb, x, xcb)
      integer :: i, Nx, i_Xpoint(:)
      real(wp) :: x_core_sol, L_core_SOL, L, flux_expansion, flux_expansion_left, trans_expansion, trans_expansion_left
      real(wp) :: B_field(:), B_field_cb(0:), B_trans(:), B_trans_cb(0:), x(:), xcb(0:)
      !allocate( B_field(Nx), B_field_cb(0:Nx), B_trans(Nx), B_trans_cb(0:Nx) )
      write(*,*) 'setting-up magnetic fields'
      ! set the magnetic field values along the grid
      B_field    = 1.0d0
      B_field_cb = 1.0d0
      B_trans 	 = 1.0d0
      B_trans_cb = 1.0d0
      if( trans_expansion .ge. 1.0d0 ) then ! NB. we staan geen transport contractie toe, i.e. trans_expansion < 1 
          ! B_field    = 1.0d0 / ( 1.0d0 + (flux_expansion - 1.0d0) * x   / L )
          ! B_field_cb = 1.0d0 / ( 1.0d0 + (flux_expansion - 1.0d0) * xcb / L )
          ! we apply the transport expansion only along the divertor-SOL
          B_trans( i_Xpoint(2)+1 : Nx )      = 1.0d+0 / ( 1.0d+0 + (trans_expansion - 1.0d+0) * (x  ( i_Xpoint(2)+1 : Nx ) - X_core_SOL - L_core_SOL) / (L - X_core_SOL - L_core_SOL) )
          B_trans_cb( i_Xpoint(2)+1 : Nx )   = 1.0d+0 / ( 1.0d+0 + (trans_expansion - 1.0d+0) * (xcb( i_Xpoint(2)+1 : Nx ) - X_core_SOL - L_core_SOL) / (L - X_core_SOL - L_core_SOL) )
          if( L_core_SOL .gt. 0.0d+0 .and. X_core_SOL .gt. 0.0d+0 ) then
	      if( trans_expansion_left .eq. 0.0d+0 ) then
		trans_expansion_left = trans_expansion
	      endif
	      ! use the second value of trans expansion for the inner target
              B_trans( 1: i_Xpoint(1)-1 )    = 1.0d+0 / ( 1.0d+0 + (trans_expansion_left - 1.0d+0) * (X_core_SOL - x  ( 1: i_Xpoint(1)-1 ) ) / X_core_SOL )
              B_trans_cb( 0: i_Xpoint(1)-1 ) = 1.0d+0 / ( 1.0d+0 + (trans_expansion_left - 1.0d+0) * (X_core_SOL - xcb( 0: i_Xpoint(1)-1 ) ) / X_core_SOL )
          endif
      elseif( trans_expansion .lt. 1.0d+0 ) then ! if trans_expansion =<0 the above will fail and DIV1D will try to load B_trans.dat
	  open(1, file = 'B_trans.dat', status = 'old')
        	do i = 1,Nx
          	read(1,*) B_trans(i)
           	B_trans(i) = max(B_trans(i),0.0d+0) ! allow negative values?
         	end do
          close(1)
	  call cc2cb(Nx,B_trans,B_trans_cb) 
      endif
      
      !if( flux_expansion .ne. 1.0d0 ) then ! NB. we staan ook flux contractie toe, i.e. flux_expansion < 1 
	if( flux_expansion .gt. 0.0d0 ) then 
          ! B_field    = 1.0d0 / ( 1.0d0 + (flux_expansion - 1.0d0) * x   / L )
          ! B_field_cb = 1.0d0 / ( 1.0d0 + (flux_expansion - 1.0d0) * xcb / L )
          ! we apply the flux expansion along the divertor-SOL only
          B_field( i_Xpoint(2)+1 : Nx )      = 1.0d0 / ( 1.0d0 + (flux_expansion - 1.0d0) * (x  ( i_Xpoint(2)+1 : Nx ) - X_core_SOL - L_core_SOL) / (L - X_core_SOL - L_core_SOL) )
          B_field_cb( i_Xpoint(2)+1 : Nx )   = 1.0d0 / ( 1.0d0 + (flux_expansion - 1.0d0) * (xcb( i_Xpoint(2)+1 : Nx ) - X_core_SOL - L_core_SOL) / (L - X_core_SOL - L_core_SOL) )
          write(*,*) B_field 
          !stop
          if( L_core_SOL .gt. 0.0d+0 .and. X_core_SOL .gt. 0.0d+0 ) then
	       if( flux_expansion_left .eq. 0.0d+0 ) then
		flux_expansion_left = flux_expansion
	      endif
              B_field( 1: i_Xpoint(1)-1 )    = 1.0d0 / ( 1.0d0 + (flux_expansion_left - 1.0d0) * (X_core_SOL - x  ( 1: i_Xpoint(1)-1 ) ) / X_core_SOL )
              B_field_cb( 0: i_Xpoint(1)-1 ) = 1.0d0 / ( 1.0d0 + (flux_expansion_left - 1.0d0) * (X_core_SOL - xcb( 0: i_Xpoint(1)-1 ) ) / X_core_SOL )
          endif
	else ! if flux_expansion =<0 the above will fail and DIV1D will try to load B_field.dat
	  open(1, file = 'B_field.dat', status = 'old')
        	do i = 1,Nx
          	read(1,*) B_field(i)
           	B_field(i) = max(B_field(i),0.0d+0) ! allow negative values?
         	end do
          close(1)
	call cc2cb(Nx,B_field,B_field_cb) 
	endif
      !endif

      ! magnetic field and transport based field are multiplied and used as one variable for the plasma equations
      B_field_cb = B_field_cb * B_trans_cb
      B_field = B_field * B_trans
   end subroutine magnetic_field

   subroutine width_surface_volume(B_trans_cb, xcb, A_int, R_cc, R_cb, sintheta_cc, sintheta_cb, Area_extern, &
                                   delta_xcb, sol_width_pol, sol_width_pol_cb, volumes, A_wet, B_field, B_field_cb, &
                                   major_radius, sintheta, sol_width_omp, vesr, vesz, vesc, Z_cc, Z_cb, nr_cc, nz_cc)
	   integer :: i,j
	   integer :: n,m
	   logical :: exists
	   real(wp) :: B_trans_cb(0:), xcb(0:), A_int(1:), R_cc(:), R_cb(0:), sintheta_cc(:), sintheta_cb(0:), Area_extern(:), &
	               delta_xcb(:), sol_width_pol(:), sol_width_pol_cb(0:), volumes(:), A_wet(:), B_field(:), B_field_cb(0:), major_radius, sintheta, sol_width_omp
       real(wp), allocatable :: vesr(:), vesz(:), vesc(:), Z_cc(:), Z_cb(:), nr_cc(:), nz_cc(:)

    	   write(*,*) 'setting up width, surface and volume'
	   !write(*,*) 'major radius'
	   !write(*,*) 'major_radius', major_radius
           if( major_radius .gt. 0.0d+0 ) then
           	R_cc = major_radius
           	R_cb = major_radius
		if( allocated(vesc) .eqv. .false.) allocate( vesc(1) )
  			vesc = -1.0d+0
		if( allocated(vesz) .eqv. .false.) allocate( vesz(1) )
  			vesz = -1.0d+0
		if( allocated(vesr) .eqv. .false.) allocate( vesr(1) )
  			vesr = -1.0d+0
		if( allocated(vesrz) .eqv. .false.) allocate( vesrz(1,3) )
  			vesrz = -1.0d+0
	   else
	   	write(*,*) 'point 3'
	   
	     	open(1, file = 'major_radius.dat', status = 'old')
        	do i = 1,Nx
          	read(1,*) R_cc(i)
           	R_cc(i)= max(R_cc(i),0.0d+0) ! allow negative values?
         	end do
             	close(1)
	     	call cc2cb(Nx,R_cc,R_cb) 

	      	inquire(file="major_height.dat", exist=exists)
		if(exists) then
      			open(1, FILE='major_height.dat', status='old')
			do i =1,Nx
			read(1,*) Z_cc(i)
			enddo
			close(1)
			call cc2cb(Nx,Z_cc,Z_cb) 
		else
			Z_cc = 0.0d+0
			Z_cb = 0.0d+0
   		endif

  	      	inquire(file="vessel_r.dat", exist=exists)
		if(exists) then
      			open(1, FILE='vessel_r.dat', status='old')
			read(1,*) (n)
			write(*,*) 'allocated vesr',allocated(vesr)
			if( allocated(vesr) .eqv. .false.) allocate( vesr(n) )
			do i =1,n
			read(1,*) (vesr(i))
			enddo
			close(1) 
		else
			if( allocated(vesr) .eqv. .false.) allocate( vesr(1) )
  			vesr = -1.0d+0
		endif

      		inquire(file="vessel_z.dat", exist=exists)
		if(exists) then
      			open(1, FILE='vessel_z.dat', status='old')
			read(1,*) (n)
	    		write(*,*) 'allocated vesz',allocated(vesz)
			if( allocated(vesz) .eqv. .false.) allocate( vesz(n) )
			do i =1,n
			read(1,*) (vesz(i))
			enddo
			close(1) 
		else
			if( allocated(vesz) .eqv. .false.) allocate( vesz(1) )
  			vesz = -1.0d+0
		endif

      		inquire(file="vessel_c.dat", exist=exists)
		if(exists) then
      			open(1, FILE='vessel_c.dat', status='old')
			read(1,*) (n)
			write(*,*) 'allocated vesc',allocated(vesc)
			if( allocated(vesc) .eqv. .false.)allocate( vesc(n) )
			do i =1,n 
			read(1,*) (vesc(i)) 
			enddo 
			close(1)
 			write(*,*) 'allocated vesrz',allocated(vesrz)
			if( allocated(vesrz) .eqv. .false.) allocate( vesrz(n,3) ) 
			vesrz(1:n,1) = vesr
			vesrz(1:n,2) = vesz
			vesrz(1:n,3) = vesc 
			write(*,*)  'done vesrz'
		else
			if( allocated(vesrz) .eqv. .false.) allocate( vesrz(1,3) )
  			vesrz = -1.0d+0
			if( allocated(vesc) .eqv. .false.) allocate( vesc(1) )
  			vesc = -1.0d+0
	
		endif


            	inquire(file="sol_normal.dat", exist=exists)
		if(exists) then
      			open(1, FILE='sol_normal.dat', status='old')
			do i =1,Nx
			read(1,*) (nr_cc(i), nz_cc(i))
			enddo
			close(1)
			!call cc2cb(Nx,n_cc,n_cb) 
		else
			nr_cc = -1.0d+0
			nz_cc = 0.0d+0
   		endif

	   endif
	  
	   write(*,*) 'sintheta gt 0 and = ', sintheta
	   if( sintheta .gt. 0.0d0 ) then  
           sintheta_cc = sintheta
	   sintheta_cb = sintheta
	   else
	     open(1, file = 'sintheta.dat', status = 'old')
        	do i = 1,Nx
          	read(1,*) sintheta_cc(i)
           	sintheta_cc(i)= max(sintheta_cc(i),0.0d+0) ! allow negative values?
         	end do
             close(1)
	     call cc2cb(Nx,sintheta_cc,sintheta_cb) 	     
	   endif

	
	  ! write(*,*) ' size delta_xcb, sintheta_cc, Rcc', size(delta_xcb), size(sintheta_cc), size(R_cc)
	   Area_extern = delta_xcb*sintheta_cc*R_cc*2*pi  ! A_extern = dz 2.pi.R = dx*sintheta 2.pi.R
	   ! width_pol *sintheta * B_transp = constant (as B_field drops out with R dependence) 
	  !write(*,*) 'i_omp = ', i_omp
	 !write(*,*) 'size_sintheta_cb', size(sintheta_cb)
	   sol_width_pol_cb = sol_width_omp*sintheta_cb(i_omp)/sintheta_cb * B_trans_cb(i_omp)/B_trans_cb 
	   sol_width_pol = (sol_width_pol_cb(0:Nx-1) + sol_width_pol_cb(1:Nx)) / 2.0d+0         
	   !write(*,*) ' sol_wid_pol, Area_extern', size(sol_width_pol), size(Area_extern)        
	   volumes = Area_extern * sol_width_pol

	   A_int = volumes/delta_xcb ! this should always use the volumes from the plasma balance to stay consistent in domains
    	   !write(*,*) 'volumes' , size(volumes)
	   A_wet(1) = A_int(1) * B_field(1) / B_field_cb(0) ! this should always scale with the plasma wetted area to stay consistent.
	   A_wet(2) = A_int(Nx) * B_field(Nx) / B_field_cb(Nx) 
   end subroutine width_surface_volume

   subroutine initialize_gas_puff(gas_puff_profile, gas_puff_location, gas_puff_width, delta_xcb)
   ! this routine needs to be updated if to be used for a case with two targets
      implicit none
      !integer,  intent(in)  :: Nx
      real(wp) :: gas_puff_normalization
      real(wp) :: gas_puff_profile(:), gas_puff_location, gas_puff_width, delta_xcb(:)
      
      ! allocate and initialize the gas puff to zero
      !allocate( gas_puff_profile(Nx) )
      write(*,*) 'F initialize gas_puff_profile'
      !write(*,*) 'gas_puff_location', gas_puff_location
      !write(*,*) 'gas_puff_width', gas_puff_width
      gas_puff_profile = 0.0
      
      !if(gas_puff_source .eq. 0.0d+0) return
      ! calculate the Gaussian profile of the source and the corresponding normalization factor
      gas_puff_profile = exp( - (x - gas_puff_location)**2 / gas_puff_width**2 )
      gas_puff_normalization = sum(gas_puff_profile * delta_xcb)
      ! set the unnormalized gas puff source
      gas_puff_profile = gas_puff_profile / gas_puff_normalization
      !write(*,*) '2 gas_puff_profile', gas_puff_profile
   end subroutine initialize_gas_puff
end module grid_data
