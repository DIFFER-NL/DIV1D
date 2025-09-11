module grid_data_suite 

    use testdrive, only : new_unittest, unittest_type, error_type, check, to_string
    use test_helper, only : testOutputMgr, testToolsMgr, newMgr

   implicit none

    type(testOutputMgr) :: console    
    type(testToolsMgr)  :: testTools

    type, private :: fileVar
        class(*), allocatable :: value
        class(*), allocatable :: array(:)
    end type
    
    logical :: failed
    
    contains 
    
    subroutine Collect_Grid_Tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
        new_unittest("initialize_grid", test_initializeGrid) &
        ]
    
    call initialize()
    end subroutine Collect_Grid_Tests
    
    subroutine initialize()
        console = newMgr()
        call beforeTest()
    end subroutine initialize
    
    !This is a prepare section
    !Prepares everything required for the test to run
    !This specific section initializes the DIV1D settings, through the input file. 
    subroutine beforeTest()
        integer, parameter :: wp = kind(1.0D0) 
        real(wp), dimension(141) :: floatinphys = 1.5
        real(wp), dimension(90) :: floatinnum = 1.5
        integer, dimension(31)  :: intinphys = 1
        integer, dimension(30)  :: intinnum = 1
        logical, dimension(20)  :: loginnum = .true. 
        integer :: call_from_extern = 0
    
        call testTools%fetchInputFile("../inputs/grid_data_inputs/input.txt", "input.txt")
        call initialize_div1d_settings(floatinnum, intinnum, loginnum, floatinphys, intinphys, call_from_extern)
        
    end subroutine beforeTest
    
    subroutine test_initializeGrid(error)
    
        use grid_data, only : initialize_grid, x, xcb, B_field, i_omp, i_Xpoint, i_baffle, mid_point, X_omp, delta_x, delta_xcb, B_field_cb, &
					    R_cc, R_cb, Area_extern, sintheta_cc, sintheta_cb, sol_width_pol, sol_width_pol_cb, volumes, &
					    gas_puff_profile, core_source_profile_Q, core_source_profile_n, Z_cc, Z_cb, nr_cc, nz_cc, vesrz, &
					    prf_imp_con, A_wet
        
        type(error_type), allocatable, intent(out) :: error
        !Variables required to run GRID
        integer, parameter :: wp = selected_real_kind(15)
        integer, parameter :: Nx = 500
        real(wp) :: B_trans(Nx), B_trans_cb(0:Nx), E_core_source_profile_Q(Nx), E_core_source_profile_N(Nx), A_int(1:Nx)
	    
	    !Variables required to go through the different processes
	    integer :: in_unit, out_unit, ios
	    character(len=256) :: lineIn, lineOut
        character(len=30)  :: varName
        character(:), allocatable  :: messageAr
        character(len=100), allocatable :: message
        integer, allocatable  :: intVal
        integer, allocatable  :: intArVal(:)
        real(wp), allocatable :: realVal
        real(wp), allocatable :: realArVal(:)
        logical, allocatable  :: boolVal
        logical :: reachedGrid
        integer :: pos = 1
        integer :: workingPos = 1
        
        !Decleration of actual values. 
        integer, allocatable  :: intValAct
        real(wp), allocatable :: realValAct
        logical, allocatable  :: boolValAct
        integer, allocatable  :: intArValAct(:)
        real(wp), allocatable :: realArValAct(:)
        
        !Extra required variables
        real( wp ) :: epsilonMargin(0:Nx) = 1.0e-12_wp
        real( wp ) :: epsilonMarginInt = 1
        
        call testTools%fetchInputFile("../inputs/grid_data_inputs/B_field.dat", "B_field.dat")
        call testTools%fetchInputFile("../inputs/grid_data_inputs/major_height.dat", "major_height.dat")
        call testTools%fetchInputFile("../inputs/grid_data_inputs/major_radius.dat", "major_radius.dat")
        call testTools%fetchInputFile("../inputs/grid_data_inputs/sintheta.dat", "sintheta.dat")
        call testTools%fetchInputFile("../inputs/grid_data_inputs/sol_normal.dat", "sol_normal.dat")
        call testTools%fetchInputFile("../inputs/grid_data_inputs/vessel_c.dat", "vessel_c.dat")
        call testTools%fetchInputFile("../inputs/grid_data_inputs/vessel_r.dat", "vessel_r.dat")
        call testTools%fetchInputFile("../inputs/grid_data_inputs/vessel_z.dat", "vessel_z.dat")
        
        call initialize_grid(Nx, x, xcb, B_field, B_field_cb, B_trans, B_trans_cb, &
		                R_cc, R_cb, Area_extern, sintheta_cc, sintheta_cb, sol_width_pol, sol_width_pol_cb, volumes, &
			        gas_puff_profile, E_core_source_profile_Q, E_core_source_profile_n, i_omp, i_Xpoint, i_baffle, A_int)
			        
        call testTools%fetchOutFileIndex("../inputs/grid_data_inputs/div1d_output.txt", "div1d_output.txt", in_unit, out_unit)
        
        do      
                if(allocated(intValAct)) then 
                    deallocate(intValAct)
                else if(allocated(realValAct)) then 
                    deallocate(realValAct)
                else if(allocated(boolValAct)) then
                    deallocate(boolValAct)
                else if(allocated(intArValAct)) then
                    deallocate(intArValAct) 
                else if(allocated(realArValAct)) then
                    deallocate(realArValAct)
                end if
                
                read(in_unit, '(A)', iostat=ios) lineIn
                if(ios /= 0) exit
                read(out_unit, '(A)', iostat=ios) lineOut
                if(ios /=0) exit
                !if(index(lineIn, "&DIV1D_GRID") > 0) then
                if(index(lineIn, "&DIV1D_GRID") > 0) then  
                    reachedGrid = .true.
                end if
                if(reachedGrid) then
                        if(.not. pos < workingPos) then
                            call testTools%fetchValFromLine(in_unit, lineIn, pos, workingPos, varName, intVal, realVal, boolVal, intArVal, realArVal)
                        end if               
                end if
                
                message = ""
                messageAr = ""
                
                !Values turned off:
                !X_omp. Nothing is assigned, though the output states it should?
                !B_Trans, B_trans_cb. Improper implementation. Need to tackle the 107*1.0000000. 
                !VESRZ. Different array rank (:,:), which doesn´t work with the current array parser
                !prf_imp_con Different array rank (:,:), which doesn´t work with the current array parser
                select case(varName)
                    case("I_XPOINT")
                        IntArValAct = i_Xpoint
                        messageAr = "The output and actual values assigned to I_Xpoint were different from each other!"

                    case("I_BAFFLE")
                        IntArValAct = i_baffle
                        messageAr = "The output and actual values assigned to I_baffle were different from each other!"

                    case("MID_POINT")
                        realValAct = mid_point
                        message = "The output and actual values assigned to MID_POINT were different from each other!"

                    case("X_OMP")
                        !realValAct = X_omp
                        message = "The output and actual values assigned to X_OMP were different from each other!"

                    case("X")
                        RealArValAct = x
                        messageAr = "The output and actual values assigned to X were different from each other!"

                    case("XCB")
                        RealArValAct = xcb
                        messageAr = "The output and actual values assigned to XCB were different from each other!"

                    case("B_FIELD")
                        RealArValAct = B_field
                        messageAr = "The output and actual values assigned to B_FIELD were different from each other!"

                    case("B_FIELD_CB")
                        RealArValAct = B_field_cb
                        messageAr = "The output and actual values assigned to B_FIELD_CB were different from each other!"

                    case("B_TRANS")
                        !RealArValAct = B_trans
                        messageAr = "The output and actual values assigned to B_TRANS were different from each other!"

                    case("B_TRANS_CB")
                        !RealArValAct = B_trans_cb
                        messageAr = "The output and actual values assigned to B_TRANS_CB were different from each other!"

                    case("AREA_EXTERN")
                        RealArValAct = Area_extern
                        messageAr = "The output and actual values assigned to AREA_EXTERN were different from each other!"

                    case("R_CC")
                        RealArValAct = R_cc
                        messageAr = "The output and actual values assigned to R_CC were different from each other!"

                    case("R_CB")
                        RealArValAct = R_cb
                        messageAr = "The output and actual values assigned to R_CB were different from each other!"

                    case("Z_CC")
                        RealArValAct = Z_cc
                        messageAr = "The output and actual values assigned to Z_CC were different from each other!"

                    case("Z_CB")
                        RealArValAct = Z_cb
                        messageAr = "The output and actual values assigned to Z_CB were different from each other!"

                    case("NR_CC")
                        RealArValAct = nr_cc
                        messageAr = "The output and actual values assigned to NR_CC were different from each other!"

                    case("NZ_CC")
                        RealArValAct = nz_cc
                        messageAr = "The output and actual values assigned to NZ_CC were different from each other!"

                    case("VESRZ")
                        !requires a different array implementation
                        !RealArValAct = vesrz
                        !messageAr = "The output and actual values assigned to VESRZ were different from each other!"

                    case("SINTHETA_CC")
                        RealArValAct = sintheta_cc
                        messageAr = "The output and actual values assigned to SINTHETA_CC were different from each other!"

                    case("SINTHETA_CB")
                        RealArValAct = sintheta_cb
                        messageAr = "The output and actual values assigned to SINTHETA_CB were different from each other!"

                    case("SOL_WIDTH_POL")
                        RealArValAct = sol_width_pol
                        messageAr = "The output and actual values assigned to SOL_WIDTH_POL were different from each other!"

                    case("SOL_WIDTH_POL_CB")
                        RealArValAct = sol_width_pol_cb
                        messageAr = "The output and actual values assigned to SOL_WIDTH_POL_CB were different from each other!"

                    case("VOLUMES")
                        RealArValAct = volumes
                        messageAr = "The output and actual values assigned to VOLUMES were different from each other!"

                    case("GAS_PUFF_PROFILE")
                        RealArValAct = gas_puff_profile
                        messageAr = "The output and actual values assigned to GAS_PUFF_PROFILE were different from each other!"

                    case("E_CORE_SOURCE_PROFILE_Q")
                        RealArValAct = E_core_source_profile_Q
                        messageAr = "The output and actual values assigned to E_CORE_SOURCE_PROFILE_Q were different from each other!"

                    case("E_CORE_SOURCE_PROFILE_N")
                        RealArValAct = E_core_source_profile_n
                        messageAr = "The output and actual values assigned to E_CORE_SOURCE_PROFILE_N were different from each other!"

                    case("PRF_IMP_CON")
                        !This is an physics parameter? Also has a different array rank. 
                        !RealArValAct = prf_imp_con
                        !messageAr = "The output and actual values assigned to PRF_IMP_CON were different from each other!"

                    case("A_INT")
                        RealArValAct = A_int
                        messageAr = "The output and actual values assigned to A_INT were different from each other!"
                        
                    case("A_WET")
                        RealArValAct = A_wet
                        messageAr = "The output and actual values assigned to A_WET were different from each other!"

                    case("I_OMP")
                        IntValAct = i_omp
                        message = "The output and actual values assigned to I_OMP were different from each other!"
                end select 
                
                pos = pos + 1
                
                if(allocated(intValAct) .and. allocated(intVal)) then
                    call check(error, intVal == intValAct)
                    if(allocated(error)) then
                        failed = .true.
                    end if
                    call console%writeOutput(Message, intValAct, intVal, failed)
                else if(allocated(realValAct) .and. allocated(realVal)) then
                    call check(error, realVal == realValAct)
                    if(allocated(error)) then
                        failed = .true.
                    end if
                    call console%writeOutput(Message, realValAct, realVal, failed)
                else if(allocated(boolVal) .and. allocated(boolValAct)) then
                    call check(error, boolVal, boolValAct)
                    if(allocated(error)) then
                        failed = .true.
                    end if
                    call console%writeOutput(Message, boolValAct, boolVal, failed)
                else if(allocated(intArVal) .and. allocated(intArValAct)) then
                    call check(error, all( intArValAct == intArVal))
                    if(allocated(error)) then
                        failed = .true.
                    end if
                    call console%writeOutput(MessageAr, intArVal, intArValAct, failed)
                else if(allocated(realArVal) .and. allocated(realArValAct)) then
                    call check(error, testTools%checkWithMargin(realArValAct, realArVal, epsilonMargin))
                    if(allocated(error)) then
                        failed = .true.
                    end if
                    call console%writeOutput(MessageAr, realArValAct, realArVal, failed)
                end if
                    if(failed) exit
        end do 
        
        close(in_unit)
        close(out_unit)
    end subroutine test_initializeGrid
end module grid_data_suite

