program tester
  use, intrinsic :: iso_fortran_env, only : error_unit
  use testdrive, only : run_testsuite, new_testsuite, testsuite_type, collect_interface, junit_output, junit_header
  use Physics_Routines_Suite, only : Collect_Physics_Tests
  use Grid_Data_Suite, only : Collect_Grid_Tests
  use Particle_Conservation_Suite, only: Collect_Conservation_Tests
  use test_helper, only : testOutputMgr, newMgr, testToolsMgr, return_help
  
  implicit none
  integer :: stat, is
  integer :: i, j, n
  type(testsuite_type), allocatable :: testsuites(:)
  integer, allocatable :: positions(:) 
  type(testOutputMgr) :: Mgr
  type(testToolsMgr) :: testTools
  type(junit_output) :: junit
  character(len=*), parameter :: fmt = '("#", *(1x, a))'
  character(len=256) :: command
  character(len=100) :: xmlFile, name
  character(len=10) :: indexstr

    stat = 0
    call junit_header (junit, "global_unit")

    allocate(testsuites(3))
    testsuites(1) = new_testsuite("Grid_Data_Suite", Collect_Grid_Tests)
    testsuites(2) = new_testsuite("Physics_Routines_Suite", Collect_Physics_Tests)
    testsuites(3) = new_testsuite("Particle_Conservation_Suite",Collect_Conservation_Tests)

    xmlfile= 'JUnit'//junit%package//'.xml'

    Mgr = newMgr() 
       if(allocated(Mgr%suites) == .true.) then
            do i = 1, size(testsuites)
                do j = 1, size(Mgr%suites)
                    if(testTools%to_upper(Mgr%suites(j)) == testTools%to_upper(testsuites(i)%name)) then
                        if(.not. allocated(positions)) allocate(positions(size(Mgr%suites)))
                            do n  = 1, size(positions)
                                positions(n) = 0
                                if(positions(n) == 0) positions(n) = i
                            end do
                    end if
                end do
            end do
            if (.not. allocated(positions)) then
                call return_help(4, "suites")
            end if
        else
            allocate(positions(size(testsuites)))
            do i = 1, size(positions)
                positions(i) = i
            end do
        end if
  do is = 1, size(positions)
    write(error_unit, fmt) "Testing:", testsuites(positions(is))%name
    call run_testsuite(testsuites(positions(is))%collect, error_unit, stat, .false., junit)
  end do
  
  do is = 1, size(testsuites)
     write(indexstr, '(I0)') is
     write(command, '(A," ",A," ",A)') 'python3', trim('../XmlManager.py'), trim(xmlfile)//' '//trim(testsuites(is)%name)//' '//indexstr
     call system(trim(command))
  end do 

  if (stat > 0) then
    write(error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
    error stop
  end if
end program tester
