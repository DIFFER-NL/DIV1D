module reaction_rates
! module containing routines implementing the reaction rates

   use numerics_parameters, only : Nx, switch_charge_exchange, switch_recombination, switch_ionization, switch_excitation
   use physics_parameters,  only : charge_exchange_model, ionization_model, recombination_model, case_AMJUEL, carbon_concentration, &
                                   minimum_temperature, minimum_density, mass, car_con_prf
   use constants,           only : e_charge

   implicit none
   integer, parameter, private :: wp = KIND(1.0D0)

   ! See page XX of AMJUEL for definition of the fit functions and variables (T in eV; density in 10^14 /m^3; resulting rates <sigma v> in cm^3/s
!   ! the coefficients for the charge exchange rate according to AMJUEL 2.1 reaction 0.1T(Total)
!   real(wp), private, dimension(9)   :: cxa_coef = (/ &
!                -1.833882000000D+01,  2.368705000000D-01, -1.469575000000D-02, -1.139850000000D-02,  6.379644000000D-04,  3.162724000000D-04, -6.681994000000D-05,  3.812123000000D-06,  8.652321000000D-09 /)
   ! the coefficients for the charge exchange rate according to AMJUEL 2.19 reaction 3.1.8 page 43
   real(wp), private, dimension(9)   :: cxa_coef = (/ &
                -1.850280000000D+01,  3.708409000000D-01,  7.949876000000D-03, -6.143769000000D-04, -4.698969000000D-04, -4.096807000000D-04,  1.440382000000D-04, -1.514243000000D-05,  5.122435000000D-07 /)

   ! the coefficients for the total hydrogen recombination rate according to AMJUEL 4.4 reaction 2.1.8 total rate including three-body recombination [m^3 / s]
   real(wp), private, dimension(9,9) :: recomb_coef = reshape( (/  &
                -2.855728479302D+01,  3.488563234375D-02, -2.799644392058D-02,  1.209545317879D-02, -2.436630799820D-03,  2.837893719800D-04, -1.886511169084D-05,  6.752155602894D-07, -1.005893858779D-08, &
                -7.664042607917D-01, -3.583233366133D-03, -7.452514292790D-03,  2.709299760454D-03, -7.745129766167D-04,  1.142444698207D-04, -9.382783518064D-06,  3.902800099653D-07, -6.387411585521D-09, &
                -4.930424003280D-03, -3.620245352252D-03,  6.958711963182D-03, -2.139257298118D-03,  4.603883706734D-04, -5.991636837395D-05,  4.729262545726D-06, -1.993485395689D-07,  3.352589865190D-09, &
                -5.386830982777D-03, -9.532840484460D-04,  4.631753807534D-04, -5.371179699661D-04,  1.543350502150D-04, -2.257565836876D-05,  1.730782954588D-06, -6.618240780594D-08,  1.013364275013D-09, &
                -1.626039237665D-04,  1.888048628708D-04,  1.288577690147D-04, -1.634580516353D-05, -9.601036952725D-06,  3.425262385387D-06, -4.077019941998D-07,  2.042041097083D-08, -3.707977721109D-10, &
                 6.080907650243D-06, -1.014890683861D-05, -1.145028889459D-04,  5.942193980802D-05, -1.211851723717D-05,  1.118965496365D-06, -4.275321573501D-08,  3.708616111085D-10,  7.068450112690D-12, &
                 2.101102051942D-05,  2.245676563601D-05, -2.245624273814D-06, -2.944873763540D-06,  1.002105099354D-06, -1.291320799814D-07,  7.786155463269D-09, -2.441127783437D-10,  3.773208484020D-12, &
                -2.770717597683D-06, -4.695982369246D-06,  3.250878872873D-06, -9.387290785993D-07,  1.392391630459D-07, -1.139093288575D-08,  5.178505597480D-10, -9.452402157390D-12, -4.672724022059D-14, &
                 1.038235939800D-07,  2.523166611507D-07, -2.145390398476D-07,  7.381435237585D-08, -1.299713684966D-08,  1.265189576423D-09, -6.854203970018D-11,  1.836615031798D-12, -1.640492364811D-14 /), &
             shape(recomb_coef), order=(/2,1/) )

   ! the coefficients for the total hydrogen ionization rate according to AMJUEL 4.3 reaction 2.1.5 [m^3 / s]
   real(wp), private, dimension(9,9) :: ionize_coef = reshape( (/  &
                -3.248025330340D+01, -5.440669186583D-02,  9.048888225109D-02, -4.054078993576D-02,  8.976513750477D-03, -1.060334011186D-03,  6.846238436472D-05, -2.242955329604D-06,  2.890437688072D-08, &
                 1.425332391510D+01, -3.594347160760D-02, -2.014729121556D-02,  1.039773615730D-02, -1.771792153042D-03,  1.237467264294D-04, -3.130184159149D-06, -3.051994601527D-08,  1.888148175469D-09, &
                -6.632235026785D+00,  9.255558353174D-02, -5.580210154625D-03, -5.902218748238D-03,  1.295609806553D-03, -1.056721622588D-04,  4.646310029498D-06, -1.479612391848D-07,  2.852251258320D-09, &
                 2.059544135448D+00, -7.562462086943D-02,  1.519595967433D-02,  5.803498098354D-04, -3.527285012725D-04,  3.201533740322D-05, -1.835196889733D-06,  9.474014343303D-08, -2.342505583774D-09, &
                -4.425370331410D-01,  2.882634019199D-02, -7.285771485050D-03,  4.643389885987D-04,  1.145700685235D-06,  8.493662724988D-07, -1.001032516512D-08, -1.476839184318D-08,  6.047700368169D-10, &
                 6.309381861496D-02, -5.788686535780D-03,  1.507382955250D-03, -1.201550548662D-04,  6.574487543511D-06, -9.678782818849D-07,  5.176265845225D-08,  1.291551676860D-09, -9.685157340473D-11, &
                -5.620091829261D-03,  6.329105568040D-04, -1.527777697951D-04,  8.270124691336D-06,  3.224101773605D-08,  4.377402649057D-08, -2.622921686955D-09, -2.259663431436D-10,  1.161438990709D-11, &
                 2.812016578355D-04, -3.564132950345D-05,  7.222726811078D-06,  1.433018694347D-07, -1.097431215601D-07,  7.789031791949D-09, -4.197728680251D-10,  3.032260338723D-11, -8.911076930014D-13, &
                -6.011143453374D-06,  8.089651265488D-07, -1.186212683668D-07, -2.381080756307D-08,  6.271173694534D-09, -5.483010244930D-10,  3.064611702159D-11, -1.355903284487D-12,  2.935080031599D-14 /), &
             shape(ionize_coef), order=(/2,1/) )

   ! the coefficients for the total hydrogen excitation rate according to AMJUEL 10.2 reaction 2.1.5 effective cooling rate by ionization and radiation [eV m^3 / s]
   real(wp), private, dimension(9,9) :: excite_coef = reshape( (/  &
                -2.497580168306D+01,  1.081653961822D-03, -7.358936044605D-04,  4.122398646951D-04, -1.408153300988D-04,  2.469730836220D-05, -2.212823709798D-06,  9.648139704737D-08, -1.611904413846D-09, &
                 1.004448839974D+01, -3.189474633369D-03,  2.510128351932D-03, -7.707040988954D-04,  1.031309578578D-04, -3.716939423005D-06, -4.249704742353D-07,  4.164960852522D-08, -9.893423877739D-10, &
                -4.867952931298D+00, -5.852267850690D-03,  2.867458651322D-03, -8.328668093987D-04,  2.056134355492D-04, -3.301570807523D-05,  2.831739755462D-06, -1.164969298033D-07,  1.785440278790D-09, &
                 1.689422238067D+00,  7.744372210287D-03, -3.087364236497D-03,  4.707676288420D-04, -5.508611815406D-05,  7.305867762241D-06, -6.000115718138D-07,  2.045211951761D-08, -1.790312871690D-10, &
                -4.103532320100D-01, -3.622291213236D-03,  1.327415215304D-03, -1.424078519508D-04,  3.307339563081D-06,  5.256679519499D-09,  7.597020291557D-10,  1.799505288362D-09, -9.280890205774D-11, &
                 6.469718387357D-02,  8.268567898126D-04, -2.830939623802D-04,  2.411848024960D-05,  5.707984861100D-07, -1.016945693300D-07,  3.517154874443D-09, -4.453195673947D-10,  2.002478264932D-11, &
                -6.215861314764D-03, -9.836595524255D-05,  3.017296919092D-05, -1.474253805845D-06, -2.397868837417D-07,  1.518743025531D-08,  4.149084521319D-10, -6.803200444549D-12, -1.151855939531D-12, &
                 3.289809895460D-04,  5.845697922558D-06, -1.479323780613D-06, -4.633029022577D-08,  3.337390374041D-08, -1.770252084837D-09, -5.289806153651D-11,  3.864394776250D-12, -8.694978774411D-15, &
                -7.335808238917D-06, -1.367574486885D-07,  2.423236476442D-08,  5.733871119707D-09, -1.512777532459D-09,  8.733801272834D-11,  7.196798841269D-13, -1.441033650378D-13,  1.734769090475D-15 /), &
             shape(excite_coef), order=(/2,1/) )

   ! the coefficients for the charge exchange rate from Freeman and Jones CLM-R-137 Table 3
   real(wp), private, dimension(9)   :: cxf_coef = (/ &
                -1.841757d+01, 5.282950d-01, -2.200477d-01, 9.750192d-02, -1.749183d-02, 4.954296d-04, 2.174910d-04, -2.530206d-05, 8.230751d-07 /)

   ! the coefficients for the ionization rate from Freeman and Jones CLM-R-137 Table 3
   real(wp), private, dimension(7)   :: ifj_coef = (/ &
                -0.3173850d+2, 0.1143818d+2, -0.3833998d+1, 0.7046692d+0, -0.7431486d-1, 0.4153749d-2, -0.9486967d-4 /)

contains

   real(wp) function charge_exchange( temperature )
   ! function to calculate the total charge_exchange rate coefficient [m^3/s]
      implicit none
      integer  :: j
      real(wp) :: temperature
      real(wp) :: ln_T, xj
      ! note that the temperature must be rescaled with the ratio of proton mass over ion mass used
      select case (charge_exchange_model)
      case ("AMJUEL")
         ! source AMJUEL page 38 2.2 reaction 0.1T
         ln_T = log(max((1.6726d-27/mass)*temperature,minimum_temperature,0.1d+0))
         xj = 1.0d+0
         charge_exchange = 0.0d+0
         do j = 1, 9
            charge_exchange = charge_exchange + cxa_coef(j)*xj
            xj = xj * ln_T
         enddo
         charge_exchange = exp(charge_exchange)*1.0d-6
      case ("Havlickova")
         ! source SD1D manual / Havlickova (2013)
         if( (1.6726d-27/mass)*temperature .le. 1.0 ) then
            charge_exchange = 1.0d-14
         else
            charge_exchange = 1.0d-14 * ((1.6726d-27/mass)*temperature)**(1/3)
         endif
      case ("Freeman")
         ! source Freeman and Jones CLM-R-137 Table 3
         ln_T = log(max((1.6726d-27/mass)*temperature,minimum_temperature,0.1d+0))
         xj = 1.0d+0
         charge_exchange = 0.0d+0
         do j = 1, 9
            charge_exchange = charge_exchange + cxf_coef(j)*xj
            xj = xj * ln_T
         enddo
         charge_exchange = exp(charge_exchange)*1.0d-6
      end select
      charge_exchange = switch_charge_exchange * charge_exchange
      return
   end function charge_exchange

   real(wp) function ionization( density, temperature )
   ! function to calculate the total ionization rate coefficient [m^3/s]
   implicit none
      integer  :: m, n
      real(wp) :: density, temperature
      real(wp) :: ln_n, ln_T, xm, xn
      ionization = 0.0d+0
      select case (ionization_model)
      case ("AMJUEL")
         ! the total hydrogen ionization rate according to AMJUEL 4.3 reaction 2.1.5 [m^3 / s]
         ln_n = log(density*1.0d-14)
         ln_T = log(max(temperature,minimum_temperature,0.1d+0))
         xm = 1.0d+0
         do m = 1, 9
            xn = 1.0d+0
            do n = 1, 9
               ionization = ionization + ionize_coef(n,m) * xm * xn
               xn = xn * ln_T
            enddo
            xm = xm * ln_n
         enddo
         ionization = exp(ionization) * 1.0d-6
      case ("Havlickova")
         ! source SD1D manual / Havlickova (2013)
         ! the discuntinuity at 20 eV has been removed by modifying the exponent of the temperature from -3.054 to -2.987
         if( temperature .le. 1.0d+0 ) then
            ionization = 7.638d-21 
         elseif( temperature .lt. 20.0d+0 ) then
            ionization = 1.0d+1**(-6.0d+0 - 15.72d+0 * exp(-log10(temperature)) + 1.603d+0*exp(-log10(temperature)**2) )* temperature**(-2.987d+0)
         else
            ionization = 5.875d-12 * temperature**(-0.5151d+0) * 10**( -2.563d+0 / log10(temperature) )
         endif
      case ("Freeman")
         ! source Freeman and Jones CLM-R-137 Table 3
         ln_T = log(max(temperature,minimum_temperature,0.1d+0))
         xn = 1.0d+0
         ionization = 0.0d+0
         do n = 1, 7
            ionization = ionization + ifj_coef(n)*xn
            xn = xn * ln_T
         enddo
         ionization = exp(ionization)*1.0d-6
      end select
      ionization = switch_ionization * ionization
      return
   end function ionization

   real(wp) function excitation( density, temperature )
   ! function to calculate the efective excitation rate coefficient [eV m^3/s]
   ! this is a measure of the effective energy loss per ionization event (so the low density high T limit is 13.6 * ionization rate)
   ! source SD1D code / Havlickova (2013)
   implicit none
      integer  :: m, n
      real(wp) :: density, temperature
      real(wp) :: ln_n, ln_T, xm, xn, Y
      excitation = 0.0d+0
      if( case_AMJUEL ) then
         ! the total hydrogen effective cooling rate by ionization and radiation according to AMJUEL 10.2 reaction 2.1.5  [eV m^3 / s]
         ln_n = log(density*1.0d-14)
         ln_T = log(max(temperature,minimum_temperature,0.1d+0))
         xm = 1.0d+0
         do m = 1, 9
            xn = 1.0d+0
            do n = 1, 9
               excitation = excitation + excite_coef(n,m) * xm * xn
               xn = xn * ln_T
            enddo
            xm = xm * ln_n
         enddo
         excitation = exp(excitation) * 1.0d-6
      else
         if( temperature .lt. 1.0d+0 ) then
            Y = 1.02d+1 / 1.0d+0 
         else
            Y = 1.02d+1 / temperature 
         endif
         excitation = 4.90d-13 / (0.28d+0+Y) *exp(-Y)*sqrt(Y*(1.0d+0+Y)) + 13.6d+0 * ionization(density, temperature) ! added to be consistent with SD1D
      endif
      excitation = switch_excitation * excitation
      return
   end function excitation

   real(wp) function recombination( density, temperature )
   ! function to calculate the total recombination rate coefficient [m^3/s]
   ! fit function from AMJUEL 4.4 reaction 2.1.8 total rate including three-body recombination
   ! note that in the AMJUEL fits densities are normalized to 10^+14 m^-3
      implicit none
      integer  :: m, n
      real(wp) :: density, temperature
      real(wp) :: ln_n, ln_T, xm, xn
      recombination = 0.0d+0
      select case (recombination_model)
      case ("AMJUEL")
         ln_n = log(density*1.0d-14)
         ln_T = log(max(temperature,minimum_temperature,0.1d+0))
         xm = 1.0d+0
         do m = 1, 9
            xn = 1.0d+0
            do n = 1, 9
               recombination = recombination + recomb_coef(n,m) * xm * xn
               xn = xn * ln_T
            enddo
            xm = xm * ln_n
         enddo
         recombination = switch_recombination * exp(recombination) * 1.0d-6
      case ("Nakazawa")
         ! use radiative recombination rate from Gordeev et al. 1977 JETP 25 204 (ref 18 of Nakazawa)
         recombination = 1.27d-19 * (13.6d+0/temperature)**1.5d+0 / ((13.6d+0/temperature) + 0.59d+0)
         ! plus the 3 body recombination rate from Hinnov et al. 1962 Phys Rev 125 795 (ref 19 of Nakazawa)
         recombination = recombination + 5.6d-39 * temperature**(-4.5d+0) * density
         recombination = switch_recombination * recombination
      end select
   end function recombination

   real(wp) function impurity_radiation( temperature, ix) 
      ! function to calculate the effective loss rate due to inpurity radiation rate coefficient [eV m^3/s]
      implicit none
      integer :: ix
      real(wp) :: temperature, log_T, carbon_radiation
      impurity_radiation = 0.0d+0
      ! case Carbon
      carbon_radiation = 0.0d+0
      if( case_AMJUEL ) then
         ! fit function Post et al. 1977 extrapolated below its validity range of 3 eV
         log_T = log10(max(temperature,minimum_temperature,0.1d+0)/1.0d+3)
         if( temperature .lt. 2.0d+1 ) then
            carbon_radiation = ( 1.965300d3 +log_T*( 4.572039d3 +log_T*( 4.159590d3 +log_T*( 1.871560d3 +log_T*( 4.173889d2 +log_T* 3.699382d1 )))))
         elseif( temperature .lt. 2.0d+2 ) then
            carbon_radiation = ( 7.467599d1 +log_T*( 4.549038d2 +log_T*( 8.372937d2 +log_T*( 7.402515d2 +log_T*( 3.147607d2 +log_T* 5.164578d1 )))))
         elseif( temperature .lt. 2.0d+3 ) then
            carbon_radiation = (-2.120151d1 +log_T*(-3.668933d-1+log_T*( 7.295099d-1+log_T*(-1.944827d-1+log_T*(-1.263576d-1-log_T* 1.491027d-1)))))
         endif
         ! transform from erg and cm to eV and m
         carbon_radiation = (1.0d-13/e_charge)*1.0d+1**carbon_radiation
      else
         ! use the fit function from SD1D
         carbon_radiation = 2.0d-31/e_charge * (max(temperature,1.0d+0)/1.0d+1)**3 / (1.0d+0 + (max(temperature,1.0d+0)/1.0d+1)**4.5d+0)
      endif

      impurity_radiation = carbon_radiation*car_con_prf(ix)
      return
   end function impurity_radiation

end module reaction_rates
