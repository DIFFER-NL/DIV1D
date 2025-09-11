module grid_data_suite
    !Get the test-drive framework
    use testdrive, only : new_unittest, unittest_type, error_type, check
    use test_helper, only : testOutputMgr, testToolsMgr, newMgr
    
    implicit none
    
    type(testOutputMgr) :: console
    type(testToolsMgr) :: testTools
    
    logical :: failed

    !privatize all tests, and only allow the "collector" to be accessed from the outside. 
    private :: test_cc2cb

    public :: Collect_grid_Tests
    real :: epsilonMargin 

    contains

    !Collect all exported unit tests. Also perform initializion actions. 
    subroutine Collect_Grid_Tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
        new_unittest("cc2cb", test_cc2cb), &
        new_unittest("cb2cc", test_cb2cc)  &
        ]
        
        call initialize() 
    end subroutine Collect_Grid_Tests
    
    !initializer for the test_helper object. 
    subroutine initialize()
        console = newMgr()
    end subroutine initialize    
    
    subroutine test_cc2cb(error)
    
        use grid_data, only: cc2cb, delta_x, delta_xcb

        type(error_type), allocatable, intent(out) :: error
        character(:), allocatable :: message
        integer, parameter :: wp = KIND(1.0D0)
        integer, parameter :: Nx = 2
        real( wp ) :: epsilonMargin(0:Nx)
        real( wp ) :: v(Nx)
        real( wp ) :: v_cb(0:Nx)
        real( wp ) :: expected(0:Nx), actual(0:Nx)
 
        delta_x = [0.002073600000003, 0.001936799999999]
        delta_xcb = [0.002142000000003, 0.002005199999999, 0.001868399999999]
        v = [17.997129000000001, 17.999065800000000]
        expected = [17.996126400000001, 17.998131600000001, 18.000000000000000]
        
        call cc2cb(Nx, v, actual)

        epsilonMargin = 10.0**-4
        call check(error, testTools%checkWithMargin(Actual, Expected, epsilonMargin))
        if(allocated(error)) then
          failed = .true.
        endif
        message = "The actual output from cc2cb differs from the expected output."
        call console%writeOutput(message, actual, expected, failed)
        failed = .false.
    end subroutine test_cc2cb
    
    subroutine test_cb2cc(error)
    
        use grid_data, only: cb2cc, delta_x, delta_xcb
        
        type(error_type), allocatable, intent(out) :: error
        character(:), allocatable :: message
        integer, parameter :: wp = KIND(1.0D0)
        integer, parameter :: Nx = 2
        real( wp ) :: epsilonMargin(0:Nx)
        real( wp ) :: Actual(Nx), Expected(Nx)
        real( wp ) :: v_cb(0:Nx)
        
        v_cb = [1,2,3]
        Expected = [1.5, 2.5] 
        
        call cb2cc(Nx, Actual, v_cb)
        
        call check(error, actual(1) == expected(1)) 
        if(allocated(error)) then
            failed = .true.
        endif
        message = "The actual output from cc2cb differs from the expected output."
        call console%writeOutput(message, actual, expected, failed)
        failed = .false. 
    end subroutine test_cb2cc
end module grid_data_suite
