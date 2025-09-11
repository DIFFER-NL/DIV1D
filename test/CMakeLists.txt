
cmake_minimum_required(VERSION 3.14) 

project(FortranTests Fortran)

enable_language(Fortran)

include(FetchContent)

FetchContent_Declare(
  test-drive
  GIT_REPOSITORY https://github.com/fortran-lang/test-drive.git
  GIT_TAG        main 
)

FetchContent_MakeAvailable(test-drive)

set(GLOBAL_SOURCES
    test_helper.f90
)

set(UNIT_TEST_SOURCES
    ./unit_suites/tester.f90
    ./unit_suites/Physics_Routines_Suite.f90
    ./unit_suites/Grid_Data_Suite.f90
    ./unit_suites/Particle_Conservation_Suite.f90
)

SET(STEP_TEST_SOURCES
    ./integration_suites/tester_Div1d_Step.f90    
    ./integration_suites/div1d_step_suite.f90    
)

set(GRID_TEST_SOURCES
    ./integration_suites/tester_Grid_Data.f90
    ./integration_suites/grid_data_suite.f90
)

set(MAIN_SOURCES
    ../obj/
    ../obj/physics_routines.o
    ../obj/grid_data.o
    ../obj/constants.o
    ../obj/numerics_parameters.o
    ../obj/physics_parameters.o
    ../obj/reaction_rates.o
    ../obj/radiative_cooling_functions_post1977.o
    ../obj/dvode_f90_m.o
    ../obj/div1d_step.o
    ../obj/experiments.o
    ../obj/interpolation.o
    ../obj/plasma_data.o
)

function (finalize_test_properties test_target_name output_dir)
   set(sources ${ARGN})
   add_executable(${test_target_name} ${sources})
   
   set_target_properties(${test_target_name} PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${output_dir}
   )
   
   target_include_directories(${test_target_name} PRIVATE ../obj)
   target_link_libraries(${test_target_name} test-drive)
endfunction() 

#To add a new suite, create a seperate SET section
#And then add finalize_test_properties(<NAME> <TEST_SOURCE> <GLOBAL_SOURCE> <MAIN_SOURCE>)
#Test source contains the test, global source contains the test helper and other important files 
#Main source contains all files that will be tested by all tests
#NOTE! the name MUST start with run_ and END with _tests. Otherwise the makefile won´t recognize the tests! NOTE!
finalize_test_properties(run_unit_tests "${CMAKE_BINARY_DIR}/unit_tests" ${UNIT_TEST_SOURCES} ${GLOBAL_SOURCES} ${MAIN_SOURCES})

#Integration tests section. Integration tests are split into their own blocks to enforce proper scoping. 
finalize_test_properties(run_step_tests "${CMAKE_BINARY_DIR}/integration_tests" ${STEP_TEST_SOURCES} ${GLOBAL_SOURCES} ${MAIN_SOURCES})
finalize_test_properties(run_grid_tests "${CMAKE_BINARY_DIR}/integration_tests" ${GRID_TEST_SOURCES} ${GLOBAL_SOURCES} ${MAIN_SOURCES})   
