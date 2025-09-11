program div1d

! DIV1D, a program to solve the time dependant evolution on a 1D fluxtube in a divertor
!
! CopyRight Notice
! author(s) of code: Egbert Westerhof,       DIFFER (c) 2019 - 2024
!                    Jens Peter Frankemolle, DIFFER (c) 2020 - 2021
!                    Gijs Lukas Derks,       DIFFER (c) 2021 - 2024
!                    Stijn Kobussen,         DIFFER (c) 2023 - 2024
! 
! CopyRight Disclaimer
! The Dutch Institute For Fundamental Energy Research (DIFFER) hereby disclaims all copyright interest in the program “DIV1D” 
! (which models the scrape-off layer plasma of magnetically confined fusion reactors) 
! written by Egbert Westerhof, Gijs Derks, Jens-Peter Frankemölle, and Stijn Kobussen being under contract at DIFFER.
!
! signature of Martin van Breukelen May 2024
! Martin van Breukelen, Managing Director of DIFFER
!
! License Notice
! This file is part of DIV1D.
! DIV1D is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License 
! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
! DIV1D is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License along with DIV1D. 
! If not, see <https://www.gnu.org/licenses>. 
!
! Reference Journal Articles
! equations from: Nakazawa et al. Plasma Phys. Control. Fusion 42 (2000) 401
! adjustments: Derks et al. Plasma Physics and Controlled Fusion 64.12 (2022): 125013.
! adjustments: Derks et al. Plasma Physics and Controlled Fusion 66.17 (2024): 055004.

! Start of program


   use constants
   use physics_parameters
   use numerics_parameters
   use grid_data
   use plasma_data
   use physics_routines
   use dvode_f90_m
   use experiments

   use div1d_step
  
   implicit none
   integer :: ibigstep
   integer :: itime 
   ! Time the program exe time excl extern sys calls (ifortran, gfortran) incl extern sys calls (solaris f90)
   real(wp) :: T1,T2
   real(wp), allocatable :: tmp_imp_con(:,:), tmp_neutral_puff(:,:), tmp_molecule_puff(:,:), tmp_nb(:,:), tmp_mb(:,:)
   logical :: exists

   ! Execution phase
   call cpu_time(T1)
   call initialize_div1d_settings( floatinnum, intinnum, loginnum,&                   ! numerics params(INPUTS)
                                   floatinphys, intinphys, &    ! physics params (INPUTS)
                                   call_from_extern)

   allocate( tmp_imp_con(5,nout), tmp_neutral_puff(5,nout), tmp_molecule_puff(5,nout), tmp_nb(5,nout), tmp_mb(5,nout) )
   ! Initialization from FORTRAN will look for .txt and .dat files
   call initialize_div1d_arrays(density, velocity, temperature, neutral, neutral_velocity, molecule,  & ! plasma params (IN/OUT)
                             x, xcb, delta_x, delta_xcb, B_field, B_field_cb, B_trans, B_trans_cb, &
			     R_cc, R_cb, Area_extern, sintheta_cc, sintheta_cb, sol_width_pol, sol_width_pol_cb, volumes, &
			     gas_puff_profile, core_source_profile_Q, core_source_profile_n, & ! grid data     (IN/OUT)  
                             init_grid_fortran, init_prof_fortran, Nx, i_omp, i_Xpoint, i_baffle, mid_point, X_omp )  
   ! Note, here the allocated variables as known by the DIV1D program are passed and will be chagned by the init subroutine

   write(*,*) 'F starting with T(10) =',temperature(10)

   start_time = 0.0
   ! Write the output text file
   inquire(file="div1d_output.txt", exist=exists)
   if(exists) then
   	open( UNIT=10, FILE='div1d_output.txt', status='old', position="append")
	else
	open( UNIT=10, FILE='div1d_output.txt', status='new')
   endif
   !call write_header
   call write_solution( start_time )
     
   write(*,*) 'F starting simulation'   
   itime = 1
   do ibigstep = 1, nout_steps
        istate = 1  
        itime = (ibigstep-1) * nout +1
        start_time = itime*delta_t ! NOTE:if this starts at 0 in external routines ELMS will not work.
   write(*,*) 'F ibigstep=', ibigstep
   write(*,*) 'F itime   =', itime
   !write(*,*) 'nout = ', nout 
   !write(*,*) 'ntime = ', ntime
   !write(*,*) 'delta_t ', delta_t
   !write(*,*) 'size dyn_imp_con', size(dyn_imp_con)
   !write(*,*) 'dyn_Q_core(1:itime+nout-1)',dyn_Q_core(itime:itime+nout-1)
   !write(*,*) 'velocity', velocity	   
   tmp_nb = dyn_nb(1:5,itime:itime+nout-1)
   tmp_mb = dyn_mb(1:5,itime:itime+nout-1)
   tmp_imp_con =  dyn_imp_con(1:5,itime:itime+nout-1)
	!write(*,*) 'neutral puff then'
   tmp_neutral_puff = dyn_neutral_puff(1:5,itime:itime+nout-1)
	!write(*,*) 'molecule puff ?'
   tmp_molecule_puff = dyn_molecule_puff(1:5,itime:itime+nout-1)
   !write(*,*) 'F main in T10)', temperature(10)
   call run_div1d(density, velocity, temperature, neutral, neutral_velocity, molecule, &     ! plasma params  ! input & output
                             Gamma_n, Gamma_mom, q_parallel, Gamma_neutral, Gamma_mom_neutral, Gamma_molecule, &! plasma params  ! Output
			     core2sol_flux, core2sol_mol, sol2extern_ion_flux,	&		     
			     Source_n, Source_v, Source_Q, Source_neutral, Source_vn, Source_molecule, &
			     extern_neutral_density, &
      			     extern_neutral_flux, sol2extern_flux, extern2sol_flux, tar2extern_flux, & ! neutral fluxes related external neutral volumes ! Output
      			     extern2core_flux, sum_sol2extern_ion_flux, & ! fluxes exchanged with the core
			     Source_extern, neutral_pump, &
			     extern_molecule_density, &
			     extern_molecule_flux, sol2extern_mol, extern2sol_mol, tar2extern_mol, & ! neutral fluxes related external neutral volumes ! Output
      			     extern2core_mol, sum_sol2extern_ion_mol,  &! fluxes exchanged with the core                            
			     Source_extern_mol, molecule_pump, &! plasma params  ! Output
			     core_density, Gamma_core2sol,  sol2core_flux, sol2core_mol, Source_core, &
                            start_time, end_time, nout, delta_t, &
	                tmp_imp_con, dyn_nu(itime:itime+nout-1),&
	                dyn_dnu(itime:itime+nout-1), tmp_nb,  tmp_mb,  dyn_gas(itime:itime+nout-1),&
	                dyn_rec(itime:itime+nout-1), dyn_qparX(itime:itime+nout-1), dyn_red_frc(itime:itime+nout-1), &
			dyn_Q_core(itime:itime+nout-1), dyn_Gamma_core(itime:itime+nout-1), dyn_core_neutral_density(itime:itime+nout-1),&
			tmp_neutral_puff, tmp_molecule_puff, dyn_core_fuelling(itime:itime+nout-1) ) 
!	write(*,*) 'dyn_imp_con', dyn_imp_con(1,2)
    !write(*,*) 'F main out T10)', temperature(10)

   ! Write solution in the output file   
   write(*,*) 'F write solution'
   write(*,*) 'end_time = ', end_time
   write(*,*) 'itime=', itime
   write(*,*) 'nout =', nout
   call write_solution( end_time )
   end do
      
   ! Request and write cpu time
   call cpu_time(T2)
   write( 10, * ) '   cpu_time   = ', T2-T1
   close(10) 
   call ys2nvt( Nx, ys, density, velocity, temperature, neutral, neutral_velocity, molecule)	
   call yr2nr( yr, extern_neutral_density, extern_molecule_density )
   call write_restart_file
   stop 'div1d completed'
	! yeah
end program div1d


subroutine write_solution( time )
   use physics_parameters, only : dyn_gas, dyn_nu, dyn_nb, dyn_mb, dyn_rec, dyn_qparX, dyn_red_frc, dyn_imp_con, num_impurities, &
				dyn_Gamma_core, dyn_Q_core, dyn_neutral_puff, dyn_molecule_puff, dyn_core_fuelling, dyn_core_neutral_density, mol_rec
   use numerics_parameters, only : Nx, delta_t
   use grid_data, only : x, xcb, B_field_cb
   use plasma_data, only : density, velocity, temperature, neutral, neutral_velocity, molecule, Gamma_n, Gamma_mom, q_parallel, Gamma_mom_neutral, & 
				extern2sol_flux, extern2sol_mol, Gamma_neutral, Gamma_molecule, &
				Source_n, Source_v, Source_Q, Source_neutral, Source_vn, Source_molecule, Source_core, &
				extern_neutral_density, sol2extern_flux, tar2extern_flux, extern_neutral_flux, sol2core_flux, core2sol_flux, &
				extern_molecule_density, sol2extern_mol, tar2extern_mol, extern_molecule_flux, sol2core_mol, core2sol_mol, &
				core_density, core_neutral_density, Gamma_core2sol, extern2core_flux, extern2core_mol,&
				sol2extern_ion_flux, sum_sol2extern_ion_flux, sum_sol2extern_ion_mol


   implicit none
   integer, parameter :: wp = KIND(1.0D0)
   integer :: i
   real(wp), intent(in) :: time

   real(wp) :: sum_extern2core_flux, sum_extern2core_mol
   integer :: itime 
   itime = time / delta_t  

   ! post process distribution between molecules and atoms
   sum_sol2extern_ion_mol = sum_sol2extern_ion_flux*mol_rec*0.5d+0
   sum_extern2core_flux = sum(extern2core_flux) ! [D/s]
   sum_extern2core_mol = 2.0d+0 * sum(extern2core_mol) ! [D/s]
	!write(*,*) 'itime', itime
	!write(*,*) 'write_solution', dyn_molecule_puff(1:5,itime);
   write( 10, *  ) 'time        = ', time
   write( 10, * ) 'dyn_gas     = ', dyn_gas(itime)
   write( 10, * ) 'dyn_nu      = ', dyn_nu(itime)
   write( 10, *) 'dyn_nb      = ', dyn_nb(1:5,itime)
   write( 10, *  ) 'dyn_mb      = ', dyn_mb(1:5,itime)
   write( 10, * ) 'dyn_rec     = ', dyn_rec(itime)
   write( 10, * ) 'dyn_qparX   = ', dyn_qparX(itime)
   write( 10, * ) 'dyn_imp_con = ', ( dyn_imp_con(i,itime), i = 1,num_impurities )
   write( 10, * ) 'dyn_q_core  = ', dyn_Q_core(itime)
   write( 10, * ) 'dyn_gam_core= ', dyn_Gamma_core(itime)
   write( 10, * ) 'dyn_neutral_puff = ', ( dyn_neutral_puff(i,itime), i = 1,5 )
   write( 10, * ) 'dyn_molecule_puff = ', ( dyn_molecule_puff(i,itime), i = 1,5 )
   write( 10, * ) 'dyn_core_fuelling = ', dyn_core_fuelling(itime)
   write( 10, * ) 'dyn_core_neutral_density = ', dyn_core_neutral_density(itime)
 
   ! SOL quantities, fluxes, sources on Nx cells
   ! write the cell centered properties (variables and sources)
   write( 10, '(A196)' ) 'X [m]  N [/m^3]  V [m/s]  T [eV]  Nn [/m^3]  Vn [m/s] Nm [/m^3] e2s_flux [#/s] e2s_mol[#/s] c2s_flux[#/s] c2s_mol[#/s] s2e_iflu[#/s] Source_n  Source_v  Source_Q [J/m3s]   source_a [1/m3s]  Source_vn source_m '
   write( 10, '(18(1PE16.9))' ) (x(i), density(i), velocity(i), temperature(i), neutral(i), neutral_velocity(i), molecule(i), & ! 7
				extern2sol_flux(i), extern2sol_mol(i), core2sol_flux(i), core2sol_mol(i), sol2extern_ion_flux(i),  & ! 5
				 Source_n(i), Source_v(i), Source_Q(i), Source_neutral(i), Source_vn(i), Source_molecule(i), i=1,Nx ) ! 6 -> 18
   ! write the cell boundary properties (fluxes)
   write( 10, '(A156)' )  '    Xcb [m]     Gamma_n    Gamma_mom [Pa]     q_parallel [J/m2s]   Gamma_neutral [#/s]  Gamma_mom_neutral [Pa]  Gamma_molecule[#/s]'
   write( 10, '(7(1PE16.9))' ) ( xcb(i), Gamma_n(i)*B_field_cb(i),Gamma_mom(i)*B_field_cb(i),q_parallel(i)*B_field_cb(i),Gamma_neutral(i), Gamma_mom_neutral(i), Gamma_molecule(i), i=0,Nx ) ! multiplied by B_field because the code uses normalized values

   ! write external neutral, densities, fluxes and sources (in the 5 chambers)
    write( 10, '(A36)' ) 'extern_neutral_flux [#/s]'
    write( 10,  '(3(1PE16.9))' ) ( extern_neutral_flux(i), i=1,3 )
    write( 10, '(A146)' ) 'sol2extern_flux [#/s]   tar2extern_flux(i) [#/s]  extern2core_flux[D/s]  extern_neutral_density [/m3]'
    write( 10,  '(4(1PE16.9))' ) ( sol2extern_flux(i), tar2extern_flux(i), extern2core_flux(i), extern_neutral_density(i), i=1,5 )
    write( 10, '(A36)' ) 'extern_molecule_flux [#/s]'
    write( 10,  '(3(1PE16.9))' ) ( extern_molecule_flux(i), i=1,3 )
    write( 10, '(A146)' ) 'sol2extern_mol [#/s]   tar2extern_mol(i) [#/s] extern2core_mol[D2/s]   extern_molecule_density [/m3]'
    write( 10,  '(4(1PE16.9))' ) ( sol2extern_mol(i), tar2extern_mol(i), extern2core_mol(i), extern_molecule_density(i), i=1,5 )
   ! the far sol ion fluxes recycle in chamber 3  following distribution mol_rec
    write( 10, '(A76)' ) 'sum sol2extern ion flux [#/s] sum sol2extern ion mol [#/s]'
    write( 10, '(2(1PE16.9))' ),  sum_sol2extern_ion_flux,  sum_sol2extern_ion_mol 
		! sol2extern_ion_flux should be added

    ! Core outcomes
    write( 10, '(A155)' ) 'Core density [/m^3]	Source_core [/s] core_neutral_density[/m3] Gamma_core2sol[#/s] sol2core_flux[D/s] sol2core_mol[D2/s] sum_extern2core_flux [D/s] sum_extern2core_mol [D/s] '
    write( 10, '(8(1PE16.9))' ) , core_density, Source_core, core_neutral_density, Gamma_core2sol,  sol2core_flux, sol2core_mol, sum_extern2core_flux, sum_extern2core_mol
    

   return
end subroutine write_solution


subroutine error_report(input_error, restart_error, time_step_error)

   implicit none
   integer, parameter :: wp = KIND(1.0D0)
   integer, intent(in) :: input_error, restart_error, time_step_error
   
   if( input_error .ne. 0 ) then
      write(*,*) 'input_error = ', input_error
      stop 'fatal error on input'
   elseif( restart_error .ne. 0 ) then
      write(*,*) 'restart_error = ', restart_error
      stop 'fatal error on restart'
   elseif( time_step_error .ne. 2 ) then
      write(*,*) 'time_step_error = ', time_step_error
      stop 'fatal error in time step'
   endif

   return
end subroutine error_report

