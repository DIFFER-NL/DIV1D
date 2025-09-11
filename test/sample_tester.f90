program tester_integration
  use, intrinsic :: iso_fortran_env, only : error_unit
  use testdrive, only : run_testsuite, new_testsuite, testsuite_type, junit_output, junit_header
  !use SUITE_suite, only : collect_suite_tests
  
  implicit none
  integer :: stat, is
  type(testsuite_type), allocatable :: testsuites(:)
  type(junit_output) :: junit_param
  character(len=*), parameter :: fmt = '("#", *(1x, a))'

  stat = 0
!Here you define the suites, and how many there are. 
!For integration tests, it is advised to only make suite per tester. This is to enforce proper test scoping/
!When creating the collector, always use collect_<SUITE_NAME>_tests
!It works without, but it is proper programming/
allocate(testsuites(1))
testsuites(1) = new_testsuite("SUITENAME", Collect_SUITE_Tests)

!This section is optional. You define the name of the XML file that will be generated.
    call junit_header(junit_param, "XMLFILENAME")

!Here the tester will loop over all suites defined in the testsuites list. 
!The .false. junit_param section is also optional. But the .false. is required if you desire an XML file. 
  do is = 1, size(testsuites)
    write(error_unit, fmt) "Testing:", testsuites(is)%name
    call run_testsuite(testsuites(is)%collect, error_unit, stat, .false., junit_param)
  end do

  if (stat > 0) then
    write(error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
    error stop
  end if

end program tester_intigration
