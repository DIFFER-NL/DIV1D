module reaction_rates
! module containing routines implementing the reaction rates

   use numerics_parameters, only : switch_charge_exchange, switch_recombination, switch_ionization

   implicit none
   integer, parameter, private :: wp = KIND(1.0D0)
   
   ! the coefficients for the total hydrogen recombination rate according to AMJUEL
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

contains

   real(wp) function charge_exchange( temperature )
   ! function to calculate the total charge_exchange rate coefficient
   ! source SD1D manual / Havlickova (2013)
      implicit none
      real(wp) :: temperature
      if( temperature .le. 1.0 ) then
         charge_exchange = 1.0d-14
      else
         charge_exchange = 1.0d-14 * temperature**(1/3)
      endif
      charge_exchange = switch_charge_exchange * charge_exchange
      return
   end function charge_exchange

   real(wp) function ionization( temperature )
   ! function to calculate the total charge_exchange rate coefficient
   ! source SD1D manual / Havlickova (2013)
   ! the discuntinuity at 20 eV has been removed by modifying the exponent of the temperature from -3.054 to -2.987
   implicit none
      real(wp) :: temperature
      if( temperature .le. 1.0d+0 ) then
         ionization = 7.638d-21 
      elseif( temperature .lt. 20.0d+0 ) then
         ionization = 1.0d+1**(-6.0d+0 - 15.72d+0 * exp(-log10(temperature)) + 1.603d+0*exp(-log10(temperature)**2) )* temperature**(-2.987d+0)
      else
         ionization = 5.875d-12 * temperature**(-0.5151d+0) * 10**( -2.563d+0 / log10(temperature) )
      endif
      ionization = switch_ionization * ionization
      return
   end function ionization

   real(wp) function recombination( density, temperature )
   ! function to calculate the total recombination rate coefficient
   ! fit function from AMJUEL
   ! note that in the AMJUEL fits densities are normalized to 10^+14 m^-3
      implicit none
      integer  :: i, j, m, n
      real(wp) :: density
      real(wp) :: temperature
      do m = 1, 9
         do n = 1, 9
            recombination = recombination + recomb_coef(n,m) * log(density*1.0d-14)**(m-1) * log(temperature)**(n-1);
         enddo
      enddo
      recombination = switch_recombination * exp(recombination) * 1.0d-6
   end function recombination

end module reaction_rates
