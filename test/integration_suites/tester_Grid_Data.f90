program tester_intigration
  use, intrinsic :: iso_fortran_env, only : error_unit
  use testdrive, only : run_testsuite, new_testsuite, testsuite_type, junit_output, junit_header
  use grid_data_suite, only : collect_grid_tests
  
  implicit none
  integer :: stat, is
  type(testsuite_type), allocatable :: testsuites(:)
  type(junit_output) :: junit
  character(len=*), parameter :: fmt = '("#", *(1x, a))'
  character(len=256) :: command
  character(len=100) :: xmlFile, name
  character(len=10) :: indexstr

  stat = 0

    allocate(testsuites(1))
    testsuites(1) = new_testsuite("Grid_Data_Suite", Collect_Grid_Tests)

    call junit_header(junit, "Grid_Data")
    xmlfile= 'JUnit'//junit%package//'.xml'
    
  do is = 1, size(testsuites)
    write(error_unit, fmt) "Testing:", testsuites(is)%name
    call run_testsuite(testsuites(is)%collect, error_unit, stat, .false., junit)
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

end program tester_intigration
