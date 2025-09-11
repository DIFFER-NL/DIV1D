module Physics_Routines_Suite
    !Get the test-drive framework
    use testdrive, only : new_unittest, unittest_type, error_type, check, to_string
    use test_helper !, only : testOutputMgr, testToolsMgr, newMgr

    implicit none

    type(testOutputMgr) :: console
    type(testToolsMgr) :: testTools


    logical :: failed

    private test_calcSolRecycleFluxes_2t,test_calcSolRecycleFluxes_1t,test_calcSolRecycleFluxesTol

    public :: Collect_Physics_Tests

    contains

    !Collect all exported unit tests. Also perform initializion actions. 
    subroutine Collect_Physics_Tests(testsuite)
    
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
        new_unittest("calcSolRecycleFluxes_2t", test_calcSolRecycleFluxes_2t), &
        new_unittest("calcSolRecycleFluxes_1t", test_calcSolRecycleFluxes_1t), &
        new_unittest("calcSolRecycleFluxesTol", test_calcSolRecycleFluxesTol) &
        ]

        call initialize()     
    end subroutine Collect_Physics_Tests

    !initializer for the test_helper object. 
    subroutine initialize()
        console = newMgr()
    end subroutine initialize    

    !unit test with GHERKIN method. 
    subroutine test_calcSolRecycleFluxes_2t(error)
        !SCENARIO (optional) describe what you're testing here. Briefly. 
	! test on recycle fluxes which should:
	!  - have correct sign convention w.r.t. pol target angle
	!  - sum up to zero 
	!NOTE:in-SOL recycle flux should be added to this routine, to test particle conservation

        !define the source files here 
        use physics_parameters, only : pol_target_angle, X_core_SOL
        use numerics_parameters, only : Nx
        use grid_data, only : A_wet, B_field_cb
        use physics_routines, only : calculate_sol_recycle_fluxes
        
        !define your variables/types here. The error type is required, and should also be passed as a parameter.
        type(error_type), allocatable, intent(out) :: error
	character(:), allocatable :: message
        integer, parameter :: wp = KIND(1.0D0)
        real(wp) :: actual(5)
        real(wp) :: Gamma_n0, Gamma_nL, rec_frac
        real(wp) :: tar2extern_flux(5)
        real(wp) :: epsilonMargin(5)
	real(wp) :: expected(5)
        Allocate(B_field_cb(0:Nx))
       
       !define your GIVEN values here.
        Nx = 20
	rec_frac = 0.0d+0 ! everything goes to reservoirs
	pol_target_angle = (/ 60 , 60 /)
        A_wet(1) = 0.5_wp
	A_wet(2) = 1.0_wp
        B_field_cb(0) = 1.0_wp
	B_field_cb(Nx) = 1.0_wp
        Gamma_n0 = -18.0_wp ! negative X direction
	Gamma_nL = 90.0_wp
	X_core_SOL = 1.0_wp 

        !define your WHEN here. This should ONLY be one line. 
        call calculate_sol_recycle_fluxes(Gamma_n0, Gamma_nL, rec_frac, tar2extern_flux)
	! question is if this routine now uses the above settings ? 
	actual = tar2extern_flux
	
        !define your THEN here 
	!expected =(/3.0 , 5.0, 0.0, 30.0, 60.0/)
	 expected =(/3.0 , 6.0, 0.0, 30.0, 60.0/) ! the correct one, but for testing purposes we leave the above one.
	epsilonMargin = 10.0**-9
        call check(error, testTools%checkWithMargin(actual, expected, epsilonMargin))
        if(allocated(error)) then
          failed = .true.
        endif
        message = "The actual output from sol_recycle_fluxes_2t differs from the expected output."
        call console%writeOutput(message, actual, expected, failed)
	  !write(*,*) 'Fails on purpose: change expected to 3,5,0,30,60 for the test to pass'
	
        deallocate(B_field_cb)
        failed = .false. 
    end subroutine test_calcSolRecycleFluxes_2t

    subroutine test_calcSolRecycleFluxes_1t(error)
	!SCENARIO (optional) describe what you're testing here. Briefly. 
	! test on recycle fluxes for geometry with single target active.
	!  - have correct sign convention w.r.t. pol target angle
	!  - sum up to zero 
	!NOTE:in-SOL recycle flux should be added to this routine, to test particle conservation

        !define the source files here 
        use physics_parameters, only : pol_target_angle, X_core_SOL
        use numerics_parameters, only : Nx
        use grid_data, only : A_wet, B_field_cb
        use physics_routines, only : calculate_sol_recycle_fluxes
        
        !define your variables/types here. The error type is required, and should also be passed as a parameter.
        type(error_type), allocatable, intent(out) :: error
	character(:), allocatable :: message
        integer, parameter :: wp = KIND(1.0D0)
        real(wp) :: actual(5)
        real(wp) :: Gamma_n0, Gamma_nL, rec_frac
        real(wp) :: tar2extern_flux(5)
        real(wp) :: epsilonMargin(5)
	real(wp) :: expected(5)
        Allocate(B_field_cb(0:Nx))
       
       !define your GIVEN values here.
        Nx = 20
	rec_frac = 0.0d+0 ! everything goes to reservoirs
	pol_target_angle = (/ 60 , 60 /)
        A_wet(1) = 0.5_wp
	A_wet(2) = 1.0_wp
        B_field_cb(0) = 1.0_wp
	B_field_cb(Nx) = 1.0_wp
        Gamma_n0 = 0.0_wp ! negative X direction
	Gamma_nL = 90.0_wp
	! test 2
	X_core_SOL = 0.0_wp
	!WHEN
 	call calculate_sol_recycle_fluxes(Gamma_n0, Gamma_nL, rec_frac, tar2extern_flux)
	actual = tar2extern_flux
	!THEN
	expected = (/0.0 , 0.0, 0.0, 30.0, 60.0/)
	! this is wrong again on purpose, with X_core_SOL = 0, recycling to the first 2 chambers should be zero
	!if(allocated(error)) then deallocate(error) endif
        epsilonMargin = 10.0**-9
	call check(error, testTools%checkWithMargin(actual, expected, epsilonMargin))
        if(allocated(error)) then
            failed = .true.
        end if
          message = "The actual output from sol_recycle_fluxes_1t differs from the expected output."
           call console%writeOutput(message, actual, expected, failed, epsilonMargin)
	   !write(*,*) 'Fails on purpose: change expected to 0,0,0,30,60 for the test to pass'

        deallocate(B_field_cb)
        failed = .false. 
    end subroutine test_calcSolRecycleFluxes_1t

    subroutine test_calcSolRecycleFluxesTol(error)
    	! example test on tolerance of parameters that should add up to zero.
    	!define the source files here 
        use physics_parameters, only : pol_target_angle, X_core_SOL
        use numerics_parameters, only : Nx
        use grid_data, only : A_wet, B_field_cb
        use physics_routines, only : calculate_sol_recycle_fluxes
        
        !define your variables/types here. The error type is required, and should also be passed as a parameter.
        type(error_type), allocatable, intent(out) :: error
	character(:), allocatable :: message
        integer, parameter :: wp = KIND(1.0D0)
        real(wp) :: Gamma_n0, Gamma_nL, rec_frac
        real(wp) :: tar2extern_flux(5)
        real(wp) :: expected_fluxes(5)
        real(wp) :: epsilonMargin(1)
	real(wp) :: actual(1), expected(1)
        Allocate(B_field_cb(0:Nx))
       
       !define your GIVEN values here.
        Nx = 20
	rec_frac = 0.0d+0 ! everything goes to reservoirs
	pol_target_angle = (/ 60 , 60 /)
        A_wet(1) = 0.5_wp
	A_wet(2) = 1.0_wp
        B_field_cb(0) = 1.0_wp
	B_field_cb(Nx) = 1.0_wp
        Gamma_n0 = -18.0_wp ! negative X direction
	Gamma_nL = 90.0_wp
	X_core_SOL = 1.0_wp 

        !define your WHEN here. This should ONLY be one line. 
        call calculate_sol_recycle_fluxes(Gamma_n0, Gamma_nL, rec_frac, tar2extern_flux)
	! question is if this routine now uses the above settings ? 
	
        !define your THEN here 
	expected_fluxes = (/3.0 , 6.0, 0.0, 30.0, 60.0/) 
	expected = 0.0
	actual = sum(abs(expected_fluxes - tar2extern_flux))
	epsilonMargin = 10.0**-9
    	call check(error, testTools%checkWithMargin(actual, expected, epsilonMargin))
	if(allocated(error)) then
            failed = .true.
        end if
          message = "There are particle(s) left over after calculating sol recycle fluxes!"
           call console%writeOutput(message, actual, expected, failed, epsilonMargin)        
        !call check(error, error_actual == 0.0_wp, "There are " // to_string(error_actual) // " particle(s) left over after calculating sol recycle fluxes! Fails on purpose: change expected to 3,5,0,30,60 for the test to pass") 
        deallocate(B_field_cb)
    end subroutine test_calcSolRecycleFluxesTol
end module Physics_Routines_Suite
