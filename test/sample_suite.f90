module sample_suite
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
    private :: test_cc2cb
    
    public :: Collect_grid_Tests
    real :: epsilonMargin 

    contains

    !Collect all exported unit tests. Also perform initializion actions. 
    subroutine Collect_Grid_Tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
    
    !Every new test that is created, will be added here.
    !The string/character is the test name that will be displayed in the console.
      testsuite = [ &
        new_unittest("cc2cb", test_cc2cb), &
        new_unittest("cb2cc", test_cb2cc)  &
        ]
        
        call initialize() 
    end subroutine Collect_Grid_Tests
    
    !initializer for the test_helper object. 
    subroutine initialize()
    !Loads the console variable with the settings provided to the test helper object
        console = newMgr()
    end subroutine initialize    
    
    subroutine test_cb2cc(error)
        !For developing a test, the gherkin method will be used.
        !This encompasses the [SCENARIO], [GIVEN], [WHEN] and [THEN] statements.
        ![SCENARIO] briefly describes what will occur during this test.
        ![SCENARIO] test the cb2cc function which should convert the cell boundary to cell ??
        
        !It is possible for different scenario's to exist for the same function.
        !This means that different test cases will be covered.
        !Like what would happen if the function receives incorrect values? How will the system respond?
        !That is also a valueable case to cover, but should be in another test. 
    
        !Specify all modules that will be used here.
        !This forces the developer to think about what is required during the test.
        !And prevents the developer from using a sum of global variables.
        use grid_data, only: cb2cc, delta_x, delta_xcb
        
        !Mention all variables used during the test
        type(error_type), allocatable, intent(out) :: error
        character(:), allocatable :: message
        integer, parameter :: wp = KIND(1.0D0)
        integer, parameter :: Nx = 2
        real( wp ) :: epsilonMargin(0:Nx)
        !To make it clear what is being tested, try to use the ACTUAL and EXPECTED names.
        !Here EXPECTED means the value that want to see at the end of the test.
        !ACTUAL being the actual value that was created during the test.
        !If actual and expected are not equal, the test fails. 
        real( wp ) :: Actual(Nx), Expected(Nx)
        real( wp ) :: v_cb(0:Nx)
        
        !Here the [GIVEN] values are created. 
        !These values are basically preconditions to running the test
        !Internal variables that must first have a value, settings that need to be configured, etc.
        !The idea is that through these [GIVEN] values, the test can be run at any time, anywhere-
        !and doesn´t depend on other tests or functions.
        !We also create our expected value beforehand, with pre set values.
        v_cb = [1,2,3]
        Expected = [1.5, 2.5] 
        
        !Here is the [WHEN] segment. 
        !Specifically for a unit test this is quite important. 
        !As here the actual action/function that is being tested, is called.
        !Under no circumstances should a second action exist within a UNIT test.   
        call cb2cc(Nx, Actual, v_cb)
        
        !Here is the [THEN] segment. We compare the ACTUAL and EXPECTED variables. 
        !We already set the expected variable beforehand.
        !The check function is a function provided by test-drive. Multiple generic function calls exist. 
        call check(error, actual(1) == expected(1)) 
        
        !Here we check if an error has occured. If it has, then we know the system has failed. 
        !The reason this is done in a strange manner is because if the system is run in very vorbose mode- 
        !we want to see the output, even if it doesn´t fail. 
        if(allocated(error)) then
            failed = .true.
        endif
        !Here you load the message with what went wrong. This will explain in the console the problem,- 
        !without having to look into the actual output. 
        !These messages will look very similair to eachother, as you typically cover the same test cases for different functions
        message = "The actual output from cc2cb differs from the expected output."
        call console%writeOutput(message, actual, expected, failed)
        
        !At the end we set the failed back to false. If we don´t, the system will think that other tests also failed. 
        failed = .false. 
    end subroutine test_cb2cc
    
    !Below is a test that can be copied. It provides the same layout as the test above.
    ! 
    !Note. Every test that is created, MUST start with the test_ prefix. 
    !It shows that this function is a test, and prevents duplicates.
    !NO error is thrown if it doesn´t start with this prefix, but this prefix is - 
    !used in the test_helper module. 
    !The test name should, if possible, also mention something about the test
    !Say the cb2cc test is being tested with a margin, it could be: test_cb2ccMargin
    ! 
    subroutine test_sample(error)
        !SCENARIO (optional) describe what you're testing here. Briefly.

        !Mention your dependencies here. Use the [only] statement to limit the dependencies.
        integer :: actual, expected
        
        
        !Here the [GIVEN] values are created.
        
        !Here is the [WHEN] segment.   
        call 
        
        !Here is the [THEN] segment. 
        call check(error, actual == expected) 
        
        if(allocated(error)) then
            failed = .true.
        endif
        
        message = ""
        call console%writeOutput(message, actual, expected, failed)
        failed = .false. 
    end subroutine test_cb2cc
end module sample_suite
