module Particle_Conservation_Suite
    !Get the test-drive framework
    use testdrive, only : new_unittest, unittest_type, error_type, check
    !Get the test_helper module. Provides types with tools to make better tests.
    !This inclues the epsilonMargin function, allowing the tester to create a margin in the test
    use test_helper, only : testOutputMgr, testToolsMgr, newMgr
    
    implicit none
    
    !Providing a name to the referenced types
    type(testOutputMgr) :: console
    type(testToolsMgr) :: testTools
    
    !An variable required to check when a test has failed.
    !The error type variable could be used too, but this made the-
    !extra verbose argument difficult to implement.
    logical :: failed

    !privatize all tests, and only allow the "collector" to be accessed from the outside. 
    private :: test_calculate_extern_fluxes
    
    public :: Collect_Conservation_Tests
    real :: epsilonMargin 

    contains

    !Collect all exported unit tests. Also perform initializion actions. 
    subroutine Collect_Conservation_Tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
    
    !Every new test that is created, will be added here.
    !The string/character is the test name that will be displayed in the console.
    ! << TEST SUITE OBJECT BELOW >>
      testsuite = [ &
        new_unittest("calculate_extern_fluxes",test_calculate_extern_fluxes) &
        !new_unittest("cc2cb", test_cc2cb), &
        !new_unittest("cb2cc", test_cb2cc)  &
        ]
        
        call initialize() 
    end subroutine Collect_Conservation_Tests
    
    !initializer for the test_helper object. 
    subroutine initialize()
    !Loads the console variable with the settings provided to the test helper object
        console = newMgr()
    end subroutine initialize    
    

    !Below is a test that can be copied. It provides the same layout as the test above.
    ! 
    !Note. Every test that is created, MUST start with the test_ prefix. 
    !It shows that this function is a test, and prevents duplicates.
    !NO error is thrown if it doesn´t start with this prefix, but this prefix is - 
    !used in the test_helper module. 
    !The test name should, if possible, also mention something about the test
    !Say the cb2cc test is being tested with a margin, it could be: test_cb2ccMargin
    ! 
    subroutine test_calculate_extern_fluxes(error)
        !SCENARIO (optional) describe what you're testing here. Briefly.
	use physics_routines, only: calculate_extern_fluxes
        !Mention your dependencies here. Use the [only] statement to limit the dependencies.
        type(error_type), allocatable, intent(out) :: error
	character(:), allocatable :: message
        integer, parameter :: wp = KIND(1.0D0)
	real(wp) :: extern_ex(3), extern_neutral(5)
        real(wp) :: actual(3), expected(3)
	real(wp) :: epsilonMargin(3)
        !Here the [GIVEN] values are created.
       	extern_neutral = (/1, 2, 3, 4, 5/)
	extern_ex = (/1 ,1, 2/)
        !Here is the [WHEN] segment.   
        call calculate_extern_fluxes(actual,extern_neutral,extern_ex)
        
        !Here is the [THEN] segment. 
	expected = (/ 4 ,-1 ,-2 /)


        !call check(error, actual == expected) 
	epsilonMargin = 1.0E-10

        call check(error,testTools%checkWithMargin(actual,expected,epsilonMargin))
	
        if(allocated(error)) then
            failed = .true.
        endif
        
        message = "The direction of external reservoir flows."
        call console%writeOutput(message, actual, expected, failed, epsilonMargin)
        failed = .false. 
    end subroutine test_calculate_extern_fluxes
end module Particle_Conservation_Suite
