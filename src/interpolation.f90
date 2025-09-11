module interpolation
! module with simple interpolation routine

   implicit none
   integer, parameter, private :: wp = KIND(1.0D0)


contains

   subroutine interpolate(x_in, f_in, N_in, x_out, f_out, N_out)
      ! given x_in and f_in, f_out is calculated with simple linear interpolation and extrapolation
      ! x_in must be in ascending order
      implicit none
      integer,  intent(in)  :: N_in, N_out
      real(wp), intent(in)  :: x_in(N_in), f_in(N_in), x_out(N_out)
      real(wp), intent(out) :: f_out(N_out)
      integer               :: i, j
      f_out = f_in(1)
      do i = 1, N_out
         do j = 1, N_in-1
            if( x_out(i) .ge. x_in(j) .and. x_out(i) .lt. x_in(j+1) ) then
               f_out(i) = f_in(j) + (f_in(j+1)-f_in(j))*(x_out(i)-x_in(j))/(x_in(j+1)-x_in(j))
            endif
         enddo
         if( x_out(i) .ge. x_in(N_in) ) then
            f_out(i) = f_in(N_in)
         endif
      enddo
      return
   end subroutine interpolate
   
end module interpolation
