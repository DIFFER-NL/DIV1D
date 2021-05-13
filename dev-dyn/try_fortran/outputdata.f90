program outputdata   
implicit none   
   integer :: Nx 
   integer :: ntime		
   real, dimension(100) :: t, qpx, nex, gas, rec, red   
   real, dimension(50) :: x, gasdis, qpalos, reddis
   real, dimension(100) :: p, q
   integer :: i  
   
   ! data  
   ntime = 100
   Nx = 50
   do i = 1,ntime  
      t(i) = i  
      qpx(i) = 21000000
      nex(i) = 2e+19
      gas(i) = 0
      rec(i) = 0.8
      red(i) = 0.8 
   end do  
   
   
   do i = 1,Nx
	x(i) = i
	gasdis(i) = 1 / Nx
	qpalos(i) = 0
	reddis(i) = 1 / Nx
   end do

   ! output data into a file 
   open(1, file = 'tdat.dat', status='new')  
      write(1,*) ' tind    ', ' qpx    ', ' nex    ', ' gas    ', ' rec    ', ' red    '
   do i = 1,ntime  
      write(1,*)  t(i), qpx(i), nex(i), gas(i), rec(i), red(i)    
   end do  
   close(1) 


   open(2, file = 'xdat.dat', status='new')
	 write(2,*) ' xind    ', ' gasdis    ', ' qpalos    ', ' reddis    '
   do i = 1,Nx
       write(2,*)  x(i), gasdis(i), qpalos(i), reddis(i) 
   end do
   close(2)

   open (3, file = 'tdat.dat', status = 'old')
   do i = 1,ntime
	read(2,*) p(i+1), q(i+1)
   end do
   
end program outputdata
