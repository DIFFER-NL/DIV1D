program outputdata   
implicit none   
   integer :: ntime		
   real, dimension(1000) :: nu 
   real  :: tmp  
   real, dimension(1000) :: p
   integer :: i  
   
   ! data  
   ntime = 10000
   do i = 1, 1000 
	nu(i) = 2e+19
   end do

   do i = 1000,ntime   
	tmp = i/1000 - 1
      	nu(i) = 2e+19 + SIN( tmp )*2e+18
   end do  
   
   
   ! output data into a file 
   open(1, file = 'nu.dat', status='new')  
   do i = 1,ntime  
      write(1,*) nu(i)     
   end do  
   close(1) 

   ! read data from file
   open (2, file = 'nu.dat', status = 'old')
   do i = 1,ntime
	read(2,*) p(i+1)
   end do
   
end program outputdata
