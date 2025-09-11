module div1d_step_suite

    use testdrive, only : new_unittest, unittest_type, error_type, check, to_string
    use test_helper, only : testOutputMgr, testToolsMgr, newMgr

   implicit none

    type(testOutputMgr) :: console    
    type(testToolsMgr)  :: testTools
    
    logical :: failed
    
    contains 
    
    subroutine Collect_Step_Tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
        new_unittest("initialize_div1d_settings", test_initDiv1dSet) &
        ]
    
    call initialize()
    end subroutine Collect_Step_Tests
    
    subroutine initialize()
        console = newMgr()
    end subroutine initialize
    
    subroutine test_initDiv1dSet(error)
        use div1d_step, only : initialize_div1d_settings
        use numerics_parameters, only : Nx, Nout, Ntime, Dxmin, Delta_T, Abstol, Reltol, &
                                        Viscosity, Method, Density_Norm, Temperature_Norm, &
                                        Velocity_Norm, Wide_Core_profile, D_new, Switch_dissociative_attachment, &
                                        Switch_Momentum_Transfer_Atoms, Evolve_Core, Evolve_Background, &
                                        Evolve_Neutral_Momentum, restart, div1d_numerics
        use plasma_data
        use constants  
        use physics_parameters 
        use experiments     
        use grid_data
        
        !Settings for the initialize_div1d_settings function.        
        type(error_type), allocatable, intent(out) :: error
        integer, parameter :: wp = kind(1.0D0) 
        real(wp), dimension(141) :: floatinphys = 1.5
        real(wp), dimension(90) :: floatinnum = 1.5
        integer, dimension(31)  :: intinphys = 1
        integer, dimension(30)  :: intinnum = 1
        logical, dimension(20)  :: loginnum = .true. 
        integer :: call_from_extern = 0
        
        !Values used for fetching input/output files.
        integer :: in_unit, out_unit, ios
        character(len=256) :: lineIn, lineOut
        character(len=100), allocatable :: message
        character(len=30)  :: varName
        integer, allocatable  :: intVal 
        integer, allocatable  :: intArVal(:)
        real(wp), allocatable :: realVal
        real(wp), allocatable  :: realArVal(:)
        logical, allocatable  :: boolVal
        integer :: pos, workingPos
        
        !Values used to automate the testing segment
        real(wp), allocatable :: realExp, realAct
        integer, allocatable  :: intExp, intAct
        logical, allocatable  :: boolExp, boolAct
        
        pos = 1
        
        !GIVEN output file exists. Copies an predetermined output file in the inputs folder.
        call testTools%fetchInputFile("../inputs/div1d_step_inputs/input.txt", "input.txt")
        
        !WHEN the init_div1d_settings is used over the provided input file.
        call initialize_div1d_settings(floatinnum, intinnum, loginnum, floatinphys, intinphys, call_from_extern)
        
        !THIS IS NOT AN WHEN STATEMENT. This fetches the output file that we are going to compare.
        !This output file is also predetermined.
        !Also opens the in and out unit.
        call testTools%fetchOutFileIndex("div1d_output.txt", "../inputs/div1d_step_inputs/div1d_output.txt", in_unit, out_unit) 
        
        !THEN the correct output file is generated.
        !Here the variable name and values are fetched for each line, compared, and an error is thrown if they are not the same.
        !Look to the select case statement to find all values that are being checked
        do 
            if(allocated(intExp) .and. allocated(intAct)) then 
                deallocate(intExp)
                deallocate(intAct)
            else if(allocated(realExp) .and. allocated(realAct)) then 
                deallocate(realExp)
                deallocate(realAct)
            else if(allocated(boolExp) .and. allocated(boolAct)) then
                deallocate(boolExp)
                deallocate(boolAct)
            end if
            read(in_unit, '(A)', iostat=ios) lineIn
            if(ios /= 0) exit
            read(out_unit, '(A)') lineOut
           
            !A check to verify that the lines of the output are similair to eachother. 
            !If this is not the case, the test could not properly check the values. 
            
            call check(error, lineIn == LineOut)
            message = "The following lines were different from eachother. Are you sure you kept the same structure and values?"
            if(allocated(error)) then
                failed = .true.
            end if
            call console%writeOutput(message, lineIn, lineOut, failed)
            if(allocated(error)) exit
            
            call testTools%fetchValFromLine(in_unit, lineIn, pos, workingPos, varName, intVal, realVal, boolVal, intArVal, realArVal)
                                
            !Empty the message, prevents old messages from clogging the system.
            message = ""
            
            !This is the THEN segment. Here you assign the actual, expected and messages.
            
            ![REMOVE ME] Currently, only the Nx and Restart values are verified. 
            !In the future all variables mentioned in the numerics and physics should be verified.
            !Maybe this can be done through a formula? Instead of creating static variables? [REMOVE ME]
            select case (varName)
                case("NX") 
                   intExp = 500
                   intAct = Nx
                   message = "Nx had a different value than was expected"
                case("RESTART")
                   boolExp = .false.
                   boolAct = restart
                   message = "Expected false for restart."
            end select 

            !Automated checker. Must first verify which value to check
            if(allocated(intExp) .and. allocated(intAct)) then
                call check(error, intAct == intExp)
                if(allocated(error)) then
                    failed = .true.
                end if
                call console%writeOutput(Message, intAct, intExp, failed)
            else if(allocated(realExp) .and. allocated(realAct)) then
                call check(error, realExp == realAct)
                if(allocated(error)) then
                    failed = .true.
                end if
                call console%writeOutput(Message, realAct, realExp, failed)
            else if(allocated(boolExp) .and. allocated(boolAct)) then
                call check(error, boolExp, boolAct)
                if(allocated(error)) then
                    failed = .true.
                end if
                call console%writeOutput(Message, boolAct, boolExp, failed)
            end if
            
            !Used to tell the "fetchValFromLine function where in the file it is.
            !Though this isn´t used, as no arrays are being requested.
            pos = pos + 1
        end do 
        
        !Verify that the arrays are scaled to the correct size, which is Nx.
        !Other values do NOT have to be checked. The allocation is an function
        !Built into Fortran, and is already tested. Here we only verify that the
        !Arrays scale to Nx, we do not verify IF they scale/allocate.
        call check(error, size(x) == Nx)
        Message = "array X is not the same size as Nx!"
        if(allocated(error)) then
            failed = .true.
        end if
        call console%writeOutput(Message, size(x), Nx, failed)
        
        !close all units once completed. 
        close(in_unit)
        close(out_unit)
    end subroutine test_initDiv1dSet
end module div1d_step_suite
