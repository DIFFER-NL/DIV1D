program outputdata   
implicit none   
   integer :: ntime
   integer, parameter :: wp = KIND(1.0D0)
   real(wp), dimension(10000) :: nu 
   real(wp), dimension(10000) :: tmp  
   real, dimension(10000) :: p
   integer :: i  
   
   ! data  
   ntime = 10000
   do i = 1, 1000 
	! nu(i) = 2e+19
	nu(i) = 0.0D+0
   end do
   
   write(*,*) 2.0d+18

   do i = 1000, ntime   
	tmp(i) = i/1000.0D0 - 1.0D0 !nu(i) = 2e+19 + SIN( tmp )*2e+18 
        nu(i) = cos(tmp(i))*2.0d+18
   end do  
   
   
   ! output data into a file 
   open(1, file = 'nu.dat', status='new')  
   do i = 1,ntime 
      write(1,*) nu(i) !, tmp(i)     
   end do  
   close(1) 

   ! read data from file
   open (2, file = 'nu.dat', status = 'old')
   do i = 1,ntime
	read(2,*) p(i+1)
   end do
   
  print *, 'precision (number of decimal digits of precision) of real(kind(1.0D0)"=', precision(nu) 
  print *, 'range (decimal range of exponant) of real(...)', range(nu)
  print *, 'max exponent of real =', maxexponent(nu)
  print *, 'min exponent of real =', minexponent(nu)


end program outputdata
