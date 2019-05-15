subroutine rk4( rhs, neq, y, start_time, end_time )
   ! simple subroutine implmenting 4th order Runge Kutta scheme (after Numerical Recipes)
   implicit none
   integer, parameter :: wp = KIND(1.0D0)
   integer, intent(in) :: neq
   real(wp), intent(in) :: start_time, end_time
   real(wp), intent(inout) :: y(neq)
   real(wp) :: dy1(neq), dym(neq), dyt(neq), yt(neq)
   real(wp) :: delta_t, dt2, dt6

   ! time step sizes
   delta_t = end_time - start_time
   dt2 = delta_t / 2.0d+0
   dt6 = delta_t / 6.0d+0

   ! step one
   call rhs( neq, start_time, y, dy1 )
   ! write(*,*) 'dy1 =', dy1
   yt = y + dt2*dy1

   ! step two
   call rhs( neq, start_time+dt2, yt, dyt )
   yt = y + dt2*dyt
   ! write(*,*) 'dyt =', dyt

   ! step three
   call rhs( neq, start_time+dt2, yt, dym )
   yt = y + delta_t*dym
   dym = dyt + dym
   ! write(*,*) 'dym =', dym

   ! step four
   call rhs( neq, start_time+delta_t, yt, dyt )
   y = y + dt6 * (dy1 + dyt + 2.0d+0*dym)
   ! write(*,*) 'dyt =', dyt

   return
end subroutine rk4
