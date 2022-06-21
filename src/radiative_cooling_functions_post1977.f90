module radiative_cooling_functions_post1977
        ! module containing radiative cooling rates for low density high temperature plasmas from post et al 1977
        
       ! use physics_parameters, only : impurity_concentration 
        use constants,          only : e_charge



        implicit none
        integer, parameter, private :: wp = KIND(1.0D0)

        real(wp), private, dimension(5) :: tmax = (/2.00E-02, 2.00E-01, 2.00E+00, 2.00E+01, 1.00E+02/)*1.0d+03 ! from keV to eV

    !    real(wp), private, dimension(3,7) :: lithium3_coef = reshape( (/ &
    !            
    !            shape(lithium3_coef), order=(/1,2/) 
               ! tmax, A(0), A(1) ,... , A(6) 
      !  real(wp), private, dimension(3,7) :: berylium4_coef = reshape( (/ &
       !         /)
       !         shape(berylium4_coef), order=(/1,2/) 
               ! tmax, A(0), A(1) ,... , A(6) 
 

      !  real(wp), private, dimension(3,7) :: boron5_coef = reshape( (/ &
      !          /)
      !          shape(boron5_coef), order=(/1,2/) 
               ! tmax, A(0), A(1) ,... , A(6) 

        real(wp), private, dimension(5,6) :: carbon6_coef = reshape( (/ &
                 1.965300E+03, 4.572035E+03, 4.159590E+03, 1.871560E+03, 4.173889E+02, 3.699382E+01, &
                 7.467599E+01, 4.549038E+02, 8.372937E+02, 7.402515E+02, 3.147607E+02, 5.164578E+01, &
                 -2.120151E+01, -3.668933E-01, 7.295099E-01, -1.944827E-01, -1.263576E-01, -1.491027E-01, &
                 -2.121979E+01, -2.346986E-01, 4.093794E-01, 7.874548E-02, -1.841379E-01, 5.590744E-02, &
                 -2.476796E+01, 9.408181E+00, -9.657446E+00, 4.999161E+00, -1.237382E+00, 1.160610E-01 /), &
                shape(carbon6_coef), order=(/2,1/) )
               ! A(0), A(1) ,... , A(6) 

        real(wp), private, dimension(5,6) :: nitrogen7_coef = reshape( (/ &
                -1.967182E+02,-3.615155E+01,-2.093912E+01,-2.093039E+01,-9.452522E+00, &  ! A(0)
                -2.429049E+02,-3.943802E+01,-5.677397E-01,-6.617905E-01,-3.583144E+01, &  ! A(1)
                -7.454123E+01,-5.564129E+00, 7.664689E-01, 1.146777E+00, 4.386446E+01, &  ! A(2)
                 3.126366E+01, 5.140343E+01,-2.610450E-01,-7.390625E-01,-2.639331E+01, &  ! A(3)
                 2.166881E+01, 4.369243E+01, 3.464473E-01, 3.042676E-01, 7.890268E+00, &  ! A(4)
                 3.300054E+00, 1.027448E+01, 6.723385E-01,-6.024562E-02,-9.366682E-01 /), &
                shape(nitrogen7_coef), order=(/1,2/)  )
                ! the order in which it fills the matrix is 1: walk right, 2: walk down

       ! real(wp), private, dimension(3,7) :: oxygen8_coef = reshape( (/ &
       !         /),&
       !         shape(oxygen8_coef), order=(/1,2/) 
               ! tmax, A(0), A(1) ,... , A(6) 

       ! real(wp), private, dimension(3,7) :: neon10_coef = reshape( (/ &
       !         /),&
       !         shape(neon10_coef), order=(/1,2/) 
       !        ! tmax, A(0), A(1) ,... , A(6) 
       ! real(wp), private, dimension(3,7) :: argon18_coef = reshape( (/ &
       !         /),&
       !         shape(argon18_coef), order=(/1,2/) 
               ! tmax, A(0), A(1) ,... , A(6) 
       ! real(wp), private, dimension(3,7) :: krypton36_coef = reshape( (/ &
      !          
       !         shape(krypton36_coef), order=(/1,2/) 
               ! tmax, A(0), A(1) ,... , A(6) 
      !  real(wp), private, dimension(3,7) :: molybdenum42_coef = reshape( (/ &
                
     !           shape(molybdenum42_coef), order=(/1,2/) 
               ! tmax, A(0), A(1) ,... , A(6) 
       ! real(wp), private, dimension(3,7) :: tin50_coef = reshape( (/ &
                
     !           shape(tin50_coef), order=(/1,2/) 
     !          ! tmax, A(0), A(1) ,... , A(6) 
     !   real(wp), private, dimension(3,7) :: xenon54_coef = reshape( (/ &
     !           
     !           shape(xenon54_coef), order=(/1,2/) 
               ! tmax, A(0), A(1) ,... , A(6) 
     !   real(wp), private, dimension(3,7) :: tungsten74_coef = reshape( (/ &
     !           
     !         shape(tungsten74_coef), order=(/1,2/) 
     !          ! tmax, A(0), A(1) ,... , A(6) 
       
contains
      real(wp) function post_radiation( temperature , log_T , impurity_Z )
        ! function to calculate the effective loss rate due to impurity
        ! radiation rate coefficient [eV m^3/s]
        ! temperature is in eV
        ! log_T is in keV

        implicit none
        integer :: Z
        integer :: n, m, impurity_Z
        real(wp) :: temperature, log_T ! post_radiation
        post_radiation = 0.0d+0
        !log_T = log10(max(temperature,minimum_temperature,0.1d+0)/1.0d+3)

        if( temperature .lt. tmax(1) ) then
            m = 1
        elseif( temperature .lt. tmax(2) ) then
            m = 2
        elseif( temperature .lt. tmax(3) ) then
            m = 3
        elseif( temperature .lt. tmax(4) ) then
            m = 4
        elseif( temperature .lt. tmax(5) ) then
            m = 5        
        else              
            post_radiation = 0.0d+1
            return
        endif        
 

        select case (impurity_Z)
        case( 3 )
                ! lithium

        case( 4 )
                ! berylium
        case( 5 )
                ! boron
        case( 6 )
                ! carbon     
                do  n = 6,2,-1       
                post_radiation = carbon6_coef(m,n)*(log_T**(n-1)) + post_radiation

           !     write(*,*) 'm = ', m
           !     write(*,*) 'n = ', n
           !     write(*,*) 'coef  = ', carbon6_coef(m,n) 
              
                enddo
                post_radiation = carbon6_coef(m,1) + post_radiation  
                ! post_radiation = post_radiation*imp_con   ! this multiplication is done in higher function           
        case( 7 )
                ! nitrogen
             do  n = 6,2,-1       
                post_radiation = nitrogen7_coef(m,n)*(log_T**(n-1)) + post_radiation

            !    write(*,*) 'm = ', m
            !    write(*,*) 'n = ', n
            !    write(*,*) 'coef  = ', nitrogen7_coef(m,n)  ! check if the correct coefficient it taken
              
                enddo
                post_radiation = nitrogen7_coef(m,1) + post_radiation  
                ! post_radiation = post_radiation*imp_con   ! this multiplication is done in higher function           
        case( 8 )
                ! oxygen

        case( 10 )
                ! neon

        case( 18 )
                ! argon

        case( 36 )
                ! krypton

        case( 42 )
                ! molybdenum

        case( 50 ) 
                ! tin

        case( 54 )
                ! xenon

        case( 74 )
                ! tungsten             

        case default
           ! no impurity radiation
           post_radiation = 0.0d+1
           return
        end select
        ! translate erg/(cm3.s) to eV/(m3.s) and change log10(Lz) to Lz
        post_radiation = (1.0d-13/e_charge)*1.0d+1**post_radiation
        ! this results in Lz (eV/[m3.s])
        return
      end function post_radiation


end module radiative_cooling_functions_post1977
