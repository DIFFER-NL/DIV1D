module radiative_cooling_functions_post1977
        ! module containing radiative cooling rates for low density high temperature plasmas from post et al 1977
        ! 25-11-2022 updated with similar format fits on rates from T. Putterich et al 2019 Nucl. Fusion 59 056013 for Argon, Xenon, Neon, Krypton
        
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
        ! The following data was fitted on data from T. Putterich et al 2019 Nucl. Fusion 59 056013
        real(wp), private, dimension(5,6) :: neon10_coef = reshape( (/ & 
                -1.846764E+01	, -2.620810E+01	, -1.938579E+01	, -1.981922E+01	, -1.934001E+01	, & 
                1.828627E+00	, -2.346441E+01	, -1.251650E-01	, 4.053980E+00	, -1.190620E-01	, & 
                1.323623E+00	, -1.130072E+01	, 3.085526E+00	, -1.047640E+01	, -1.198263E+00	, & 
                -1.129022E+00	, 2.096207E+01	, -4.705794E+00	, 1.015768E+01	, 1.305294E+00	, & 
                -8.444497E-01	, 2.192461E+01	, -5.342030E+00	, -4.390731E+00	, -4.883286E-01	, & 
                -9.890098E-02	, 5.564161E+00	, 4.331828E+00	, 7.202625E-01	, 6.544618E-02	 /),& 
                shape(neon10_coef), order=(/1,2/) ) 

       ! the Argon rates are from Post et al. but have been extended using ATOMIC data
       ! The ATOMIC data was located via Garland et al. Phys. Plasmas 27 040702 (2020) in oktober 2022 by GL derks
       ! On 01-10-2022 the first column was fitted on data from https://www-amdis.iaea.org/cgi-bin/EFFRATES/ss.pl
       ! The data used was for an electron density of 1e14 [1/cm3] only requesting total radiated power.
       ! real(wp), private, dimension(5,6) :: argon18_coef = reshape( (/ &
       !         -35.605350107279E+00,-2.053043E+01, -1.965204E+01,-1.974883E+01, -2.117935E+01,& !A(0)
       !         -10.502526273607E+00, -2.834287E+00,  -1.172763E-01, 2.964839E+00,  5.191481E+00,& !A(1)
       !          16.127726700431E+00,  1.506902E+01,   7.833220E+00,-8.829391E+00, -7.439717E+00,& !A(2)
       !          17.647892884189E+00,  3.517177E+01,  -6.351577E+00, 9.791004E+00,  4.969023E+00,& !A(3)
       !          6.0667840339606E+00,  2.400122E+01,  -3.058849E+01,-4.960018E+00, -1.553180E+00,& !A(4)
       !          0.7391044788185E+00,  5.072723E+00,  -1.528534E+01, 9.820032E-01,  1.877047E-01 /),&!A(5)
       !         shape(argon18_coef), order=(/1,2/) )
        ! first column was thus FITTED on ATOMIC data.
        ! The following data was fitted on data from T. Putterich et al 2019 Nucl. Fusion 59 056013
        real(wp), private, dimension(5,6) :: argon18_coef = reshape( (/ & 
                -1.846764E+01	, -2.620810E+01	, -1.938579E+01	, -1.981922E+01	, -1.934001E+01	, & 
                1.828627E+00	, -2.346441E+01	, -1.251650E-01	, 4.053980E+00	, -1.190620E-01	, & 
                1.323623E+00	, -1.130072E+01	, 3.085526E+00	, -1.047640E+01	, -1.198263E+00	, & 
                -1.129022E+00	, 2.096207E+01	, -4.705794E+00	, 1.015768E+01	, 1.305294E+00	, & 
                -8.444497E-01	, 2.192461E+01	, -5.342030E+00	, -4.390731E+00	, -4.883286E-01	, & 
                -9.890098E-02	, 5.564161E+00	, 4.331828E+00	, 7.202625E-01	, 6.544618E-02	 /),& 
                shape(argon18_coef), order=(/1,2/) ) 

        ! The following data was fitted on data from T. Putterich et al 2019 Nucl. Fusion 59 056013
        real(wp), private, dimension(5,6) :: krypton36_coef = reshape( (/ & 
                -1.846764E+01	, -2.620810E+01	, -1.938579E+01	, -1.981922E+01	, -1.934001E+01	, & 
                1.828627E+00	, -2.346441E+01	, -1.251650E-01	, 4.053980E+00	, -1.190620E-01	, & 
                1.323623E+00	, -1.130072E+01	, 3.085526E+00	, -1.047640E+01	, -1.198263E+00	, & 
                -1.129022E+00	, 2.096207E+01	, -4.705794E+00	, 1.015768E+01	, 1.305294E+00	, & 
                -8.444497E-01	, 2.192461E+01	, -5.342030E+00	, -4.390731E+00	, -4.883286E-01	, & 
                -9.890098E-02	, 5.564161E+00	, 4.331828E+00	, 7.202625E-01	, 6.544618E-02	 /),& 
                shape(krypton36_coef), order=(/1,2/) ) 

      !  real(wp), private, dimension(3,7) :: molybdenum42_coef = reshape( (/ &
                
     !           shape(molybdenum42_coef), order=(/1,2/) 
               ! tmax, A(0), A(1) ,... , A(6) 
       ! real(wp), private, dimension(3,7) :: tin50_coef = reshape( (/ &
                
     !           shape(tin50_coef), order=(/1,2/) 
     !          ! tmax, A(0), A(1) ,... , A(6) 

        ! The following data was fitted on data from T. Putterich et al 2019 Nucl. Fusion 59 056013
        real(wp), private, dimension(5,6) :: xenon54_coef = reshape( (/ & 
                -1.846764E+01	, -2.620810E+01	, -1.938579E+01	, -1.981922E+01	, -1.934001E+01	, & 
                1.828627E+00	, -2.346441E+01	, -1.251650E-01	, 4.053980E+00	, -1.190620E-01	, & 
                1.323623E+00	, -1.130072E+01	, 3.085526E+00	, -1.047640E+01	, -1.198263E+00	, & 
                -1.129022E+00	, 2.096207E+01	, -4.705794E+00	, 1.015768E+01	, 1.305294E+00	, & 
                -8.444497E-01	, 2.192461E+01	, -5.342030E+00	, -4.390731E+00	, -4.883286E-01	, & 
                -9.890098E-02	, 5.564161E+00	, 4.331828E+00	, 7.202625E-01	, 6.544618E-02	 /),& 
                shape(xenon54_coef), order=(/1,2/) ) 

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
	!write(*,*) 'inside radiative cooling post'
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
                !if( temperature .lt. 3.00E+01 ) then
                !m = 1
                !endif 
                ! the rates from T. Putterich et al 2019 Nucl. Fusion 59 056013 were fitted on the same TLim as post uses for Carbon

                do  n = 6,2,-1       
                post_radiation = argon18_coef(m,n)*(log_T**(n-1)) + post_radiation

            !    write(*,*) 'm = ', m
            !    write(*,*) 'n = ', n
            !    write(*,*) 'coef  = ', argon_18_coef(m,n)  ! check if the correct coefficient it taken
              
                enddo
                post_radiation = argon18_coef(m,1) + post_radiation  
                ! post_radiation = post_radiation*imp_con   ! this multiplication is done in higher function   

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
