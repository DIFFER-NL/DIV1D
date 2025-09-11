module test_helper

    implicit none    
    private 
    public :: newMgr, return_help
    
    type, public :: testOutputMgr
    logical :: verbose
    logical :: veryVerbose
    character(len=50), allocatable, public :: suites(:)
    contains        
        procedure, private :: showOutputScalar
        procedure, private :: showOutputSP
        
        procedure, private :: writeActExpToConsoleScalar
        procedure, private :: writeActExpToConsoleSP
        
        procedure, public :: writeValueScalar
        procedure, public :: writeValueArray
        
        generic, public :: writeOutput => showOutputScalar, showOutputSP
        generic, public :: writeActExp => writeActExpToConsoleScalar, writeActExpToConsoleSP
        generic, public :: writeValue => writeValueScalar, writeValueArray
    end type testOutputMgr
    
    type, public ::  testToolsMgr
    contains
        procedure, public :: to_upper
        procedure, private :: checkWithMargin_Int
        procedure, private :: checkWithMargin_RealSP 
        procedure, private :: checkwithMargin_RealDP
        procedure, private :: checkwithMargin_RealXDP
        procedure, public  :: fetchValFromLine
        procedure, public  :: fetchInputFile
        procedure, public  :: fetchOutFileIndex
        generic, public :: checkWithMargin => checkWithMargin_RealSP, checkWithMargin_RealDP, checkwithMargin_RealXDP, checkWithMargin_Int
    end type testToolsMgr
    
    contains 
    
        function newMgr() result(obj)
            type(testOutputMgr) :: obj
            integer :: i
            integer :: suitesCount
            character(len=50) :: arg
            logical :: capture_suite, capture_test
            
            obj%verbose = .false.
            capture_suite = .false.
            capture_test = .false.
         
            do i = 1, command_argument_count()
                call get_command_argument(i, arg)
                arg = trim(adjustl(arg))
                select case (arg)
                    case("--vv")
                        obj%veryVerbose = .true.
                    case("--v")
                        if(.not. obj%veryVerbose) then
                           obj%verbose = .true.
                        end if
                    case("--s")
                        capture_suite = .true.
                    case("--t")
                        capture_suite = .false.
                        capture_test = .true.
                    case("--help") 
                        call return_help(0)
                        exit
                    case default
                        if( &
                            .not. capture_suite .and. &
                            .not. capture_test .and. &
                            index(arg, "_suite") == 0 .and. &
                            index(arg, "test_") == 0 &
                          ) then
                            call return_help(5, trim(arg)) 
                            exit
                        else if(capture_suite .and. index(arg, "_suite") > 0) then 
                            call add_suite(obj%suites, trim(arg))
                        else if(capture_suite) then 
                            call return_help(1, "suite", trim(arg))
                        else if(capture_test .and. index(arg, "test_") > 0 ) then
                            call return_help(2)
                        else if(capture_test) then 
                            call return_help(1, "test", trim(arg))
                        end if  
                        
                end select
            end do
        end function newMgr
        
        subroutine return_help(id, flag, arg)
            integer, intent(in) :: id
            character(len=*), optional :: flag
            character(len=10) :: message
            character(len=50), optional :: arg
            character(len=100) :: context
            if(present(flag)) flag = trim(flag)
            if(present(arg)) arg = trim(arg)
            select case (id) 
                case(1)
                    context = " Unknown argument encountered for "//flag//": "//arg
                case(2) 
                    context = " not implemented"
                case(3) 
                    context = " Only up to 5 suites can be specified"
                case(4)
                    context = "Could not find suitable "//flag//" with the provided "//flag
                case(5)
                    context = "argument '"//CHAR(27)//"[31m"//flag//CHAR(27)//"[0m"//"' not recognized"
            end select

            if(.not. id == 0) then
                write(*,'(A)')CHAR(27)//"[31m [ERROR] "//CHAR(27)//"[0m"//context
            end if
            write(*,'(A)')"Usage: test_helper [OPTIONS]"
            write(*,'(A)')"Options"
            write(*,'(A)')"--v       Enable verbose output. Shows expected/actual for failed tests."
            write(*,'(A)')"--vv      Enable very verbose output. Shows all expected actual regardless of failure. It is advised to use --s with this."
            write(*,'(A)')"--s       Specify test suite. Files must with ""_suite"" prefix. Only 5 suites can be specified."
            write(*,'(A)')"--t       Specify which test to run from specified suite. --s flag is required! Tests must begin with the ""test_"" prefix."
            
            stop        
        end subroutine 
        
        subroutine add_suite(suites, new_suite)
            character(len=50), dimension(:), allocatable :: suites, container
            character(len=50) :: new_suite
            integer :: new_size
            
            new_size = size(suites) + 1
            if(new_size >= 5) then 
                call return_help(3)
            end if
            allocate(container(size(suites)))
            container = suites
            if(allocated(suites)) deallocate(suites)
            allocate(suites(new_size))
            
            suites(1:size(suites)-1) = container
            suites(new_size) = new_suite
        end subroutine add_suite
                
        subroutine showOutputScalar(this, Message, Actual, Expected, error, epsilonMargin)        
            class(testOutputMgr) :: this
            class(*), intent(in) :: Actual, Expected
            class(*), intent(in), optional :: epsilonMargin
            character(len=*), intent(inout) :: Message
            logical :: error
            
            if(this%veryVerbose == .true.) then       
                if(.not. error) then 
                    message = "No error."
                end if
                call this%writeActExp(Message, Actual, Expected, epsilonMargin)
            else if(error == .true. .and. this%verbose == .true.) then
                call this%writeActExp(Message, Actual, Expected, epsilonMargin)
            end if        
        end subroutine showOutputScalar
        
        subroutine showOutputSP(this, Message, Actual, Expected, error, epsilonMargin)        
            class(testOutputMgr) :: this
            class(*), intent(in), dimension(:) :: Actual, Expected
            class(*), intent(in), dimension(:), optional :: epsilonMargin
            character(len=*), intent(inout) :: Message
            logical :: error
            
            if(this%veryVerbose == .true.) then       
                if(.not. error) then 
                    message = "No error."
                end if
                call this%writeActExp(Message, Actual, Expected, epsilonMargin)
            else if(error == .true. .and. this%verbose == .true.) then
                call this%writeActExp(Message, Actual, Expected, epsilonMargin)
            end if        
        end subroutine showOutputSP        
        
        subroutine writeActExpToConsoleScalar(this, Message, Actual, Expected, epsilonMargin) 
            class(testOutputMgr) :: this
            
            class(*), intent(in) :: Actual, Expected
            class(*), intent(in), optional :: epsilonMargin
            character(len=*), intent(in) :: Message
                     
            call seperator(len(message))
            call this%writeValue(CHAR(27)//"[31m"//message//CHAR(27)//"[0m")
            call seperator(LEN(message))
            
            call this%writeValue("Actual")
            call this%writeValue(Actual)
            call seperator(LEN(message))
            
            call this%writeValue("Expected") 
            call this%writeValue(Expected)
            call seperator(LEN(message))
            
            if(present(epsilonMargin)) then 
                call this%writeValue("epsilonMargin")
                call this%writeValue(epsilonMargin)
                call seperator(LEN(message))
            end if
        end subroutine writeActExpToConsoleScalar
        
        subroutine writeActExpToConsoleSP(this, Message, Actual, Expected, epsilonMargin)
            class(testOutputMgr) :: this
            
            class(*), dimension(:), intent(in) :: Actual, Expected 
            class(*), dimension(:), intent(in), optional :: epsilonMargin
            character(len=*), intent(in) :: Message
                     
            call seperator(len(message))
            call this%writeValue(CHAR(27)//"[31m"//message//CHAR(27)//"[0m")
            call seperator(LEN(message))
            
            call this%writeValue("Actual")
            call this%writeValue(Actual)
            call seperator(LEN(message))
            
            call this%writeValue("Expected") 
            call this%writeValue(Expected)
            call seperator(LEN(message))
            
            if(present(epsilonMargin)) then 
                call this%writeValue("epsilonMargin")
                call this%writeValue(epsilonMargin)
                call seperator(LEN(message))
            end if
        end subroutine writeActExpToConsoleSP
        
        !An function that is accessable from the outside.
        !Multiple variables CAN be passed into this.
        !Example: 
        !write(msg, *) integerVar, .false., realVar
        !call writeValueScalar(msg)
        subroutine writeValueScalar(this, Value, input, formatting) 
            use iso_fortran_env, only: real64
            
            class(*), intent(in) :: Value
            class(testOutputMgr) :: this
            
            integer, intent(in), optional :: input
            character(len=*), intent(in), optional :: formatting
            
            integer :: out_unit
            character(len=100) :: local_fmt
            
            select type(v1 => Value)
                type is (integer)
                    if(present(input) .and. present(formatting)) then
                            write(input, formatting) v1
                        else if(present(input)) then
                            write(input,*) v1
                        else if(present(formatting)) then
                            write(*,formatting) v1
                        else 
                            write(*,*) v1
                    end if    
                type is (logical)
                    if(present(input) .and. present(formatting)) then
                            write(input, formatting) v1
                        else if(present(input)) then
                            write(input,*) v1
                        else if(present(formatting)) then
                            write(*,formatting) v1
                        else 
                            write(*,*) v1
                    end if     
                type is (real(real64))
                    if(present(input) .and. present(formatting)) then
                            write(input, formatting) v1
                        else if(present(input)) then
                            write(input,*) v1
                        else if(present(formatting)) then
                            write(*,formatting) v1
                        else 
                            write(*,*) v1
                    end if  
                type is (character(len=*))
                    if(present(input) .and. present(formatting)) then
                            write(input, formatting) v1
                        else if(present(input)) then
                            write(input,*) v1
                        else if(present(formatting)) then
                            write(*,formatting) v1
                        else 
                            write(*,*) v1
                    end if 
               class default 
                    write(*,*) "Type not recognized" 
            end select
        end subroutine writeValueScalar
        
        !An function that is accessable from the outside.
        !Multiple variables CAN be passed into this.
        !Example: 
        !write(msg, *) integerVarAr, .false., realVarAr
        !call writeValueScalar(msg)
        subroutine writeValueArray(this, Value, input, formatting) 
            use iso_fortran_env, only: real64
        
            class(*), dimension(:), intent(in) :: Value
            class(testOutputMgr) :: this
            
            integer, intent(in), optional :: input
            character(len=*), intent(in), optional :: formatting
            
            integer :: out_unit
            character(len=100) :: local_fmt
            
            select type(v1 => Value)
                type is (integer)
                    if(present(input) .and. present(formatting)) then
                            write(input, formatting) v1
                        else if(present(input)) then
                            write(input,*) v1
                        else if(present(formatting)) then
                            write(*,formatting) v1
                        else 
                            write(*,*) v1
                    end if    
                type is (logical)
                    if(present(input) .and. present(formatting)) then
                            write(input, formatting) v1
                        else if(present(input)) then
                            write(input,*) v1
                        else if(present(formatting)) then
                            write(*,formatting) v1
                        else 
                            write(*,*) v1
                    end if     
                type is (real(real64))
                    if(present(input) .and. present(formatting)) then
                            write(input, formatting) v1
                        else if(present(input)) then
                            write(input,*) v1
                        else if(present(formatting)) then
                            write(*,formatting) v1
                        else 
                            write(*,*) v1
                    end if  
                type is (character(len=*))
                    if(present(input) .and. present(formatting)) then
                            write(input, formatting) v1
                        else if(present(input)) then
                            write(input,*) v1
                        else if(present(formatting)) then
                            write(*,formatting) v1
                        else 
                            write(*,*) v1
                    end if 
               class default 
                    write(*,*) "Type not recognized" 
            end select
        end subroutine writeValueArray
        
        !Turns all characters from the input into capatlized characters
        !Useful for character comparison
        function to_upper(this, input) result(output)
            class(testToolsMgr) :: this
            
            character(len=*), intent(in) :: input
            character(len=len(input)) :: output
            integer :: i, char_code
            
            output = input
            
            do i = 1, len_trim(input)
                char_code = ichar(input(i:i))
                if(char_code >= 97 .and. char_code <= 122) then
                    output(i:i) = char(char_code - 32)
                end if
            end do       
        end function to_upper
        
        function checkWithMargin_Int(this, Actual, Expected, epsilonMargin) result(value)
            class(testToolsMgr) :: this

            logical :: value
            integer, intent(in) :: Actual(:), Expected(:)
            integer, intent(in) :: epsilonMargin(:)   
            
            value = all(abs(actual - expected) < epsilonMargin)
        end function checkWithMargin_Int
        
        function checkWithMargin_RealSP(this, Actual, Expected, epsilonMargin) result(value)
            class(testToolsMgr) :: this

            logical :: value
            integer, parameter :: wp = selected_real_kind(6)
            real( wp ), intent(in) :: Actual(:), Expected(:)
            real( wp ), intent(in) :: epsilonMargin(:)        
            
            value = all(abs(actual - expected) < epsilonMargin)
        end function checkWithMargin_RealSP
        
        function checkWithMargin_RealDP(this, Actual, Expected, epsilonMargin) result(value)
            class(testToolsMgr) :: this
        
            logical :: value
            integer, parameter :: wp = selected_real_kind(15)
            real( wp ), intent(in) :: Actual(:), Expected(:)
            real( wp ), intent(in) :: epsilonMargin(:)        
            
            value = all(abs(actual - expected) < epsilonMargin)
        end function checkWithMargin_RealDP
        
        function checkWithMargin_RealXDP(this, Actual, Expected, epsilonMargin) result(value)
            class(testToolsMgr) :: this

            logical :: value
            integer, parameter :: wp = selected_real_kind(18)
            real( wp ), intent(in) :: Actual(:), Expected(:)
            real( wp ), intent(in) :: epsilonMargin(:)        
            
            value = all(abs(actual - expected) < epsilonMargin)
        end function checkWithMargin_RealXDP

        subroutine seperator(Width)
            integer :: width
            integer, parameter :: max_width = 150

            width = min(width, max_width)
            write(*,'(A)') repeat("-",width)
        end subroutine seperator
        
        subroutine fetchValFromLine(this, fileUnit, line, receivedPos, workingPos, &
                                varName, intVal, realVal, boolVal, intArVal, realArVal)
        
            class(testToolsMgr) :: this
            
            !Variables that will be returned/modified 
            integer, parameter :: wp = kind(1.0D0) 
            character(len=256) :: valueStr
            integer, allocatable  :: intVal
            integer, allocatable, optional  :: intArVal(:)
            real(wp), allocatable :: realVal
            real(wp), allocatable, optional :: realArVal(:) 
            logical, allocatable  :: boolVal
            
            !Varables required to make the subroutine function
            integer  :: eq_pos, ios, i, fileUnit
            integer  :: receivedPos, workingPos
            integer  :: count
            logical  :: isReal
            logical  :: isArray = .false.          
            logical  :: allUnallocated
            character(len=*) :: line, varName
            
            count = 0            
            
            !ensure the variables are NOT allocated.
             if(allocated(intVal)) then
                deallocate(intVal)
             end if
             if(allocated(realVal)) then
                deallocate(realVal) 
             end if
             if(allocated(boolVal)) then
                deallocate(boolVal)
             end if
             if(present(realArVal)) then
                if(allocated(realArVal)) deallocate(realArVal)
             end if 
             if(present(intArVal)) then
                if(allocated(intArVal)) deallocate(intArVal)
             end if
            eq_pos = index(line, "=")
            !check if an equal sign is found.
            if(eq_pos > 0) then
                call splitNameVar(varName, valueStr, line, eq_pos)
                
                call lineValueSize(valueStr, count)
                
                if(count > 1) then          
                    isArray = .true.
                else
                    isArray = .false.
                end if
                isReal = index(valueStr, ".") > 0 .or. &
                         index(valueStr, "e") > 0 .or. &
                         index(valueStr, "E") > 0 .or. &
                         index(valueStr, "d") > 0 .or. &
                         index(valueStr, "D") > 0
                !Check whether an integer or an real value is being read.
                if(.not. isReal) then 
                    if(.not. isArray) then
                        allocate(intVal)
                        read(valueStr, *, iostat=ios) intVal
                    else 
                        call parseArray(fileUnit, isreal, receivedPos, workingPos, line, realArVal, intArVal)
                    end if    
                else if(isReal) then
                    if(.not. isArray) then 
                        allocate(realVal)
                        read(valueStr, *, iostat=ios) realVal
                    else
                        call parseArray(fileUnit, isReal, receivedPos, workingPos, line, realArVal, intArVal)
                    end if
               end if 
                          
               !boolean value has to be handled seperately. 
               allocate(boolVal)
               read(valueStr, *, iostat=ios) boolVal
               if(ios /= 0) then 
                    deallocate(boolVal)
               end if
           
                allUnallocated = .not. allocated(intVal) &
                 .and. .not. allocated(realVal) &
                 .and. .not. allocated(boolVal)
                if (present(realArVal)) then
                    allUnallocated = allUnallocated .and. .not. allocated(realArVal)
                end if

                if (present(intArVal)) then
                    allUnallocated = allUnallocated .and. .not. allocated(intArVal)
                end if

                if (allUnallocated) then
                    write(*,*) "Error in allocating value fetched from file. Stopping..."
                    write(*,*) varName
                    stop
                end if
            end if    
        end subroutine fetchValFromLine
    
        subroutine parseArray(fileUnit, isReal, receivedPos, workingPos, line, realArVal, intArVal)
            integer, parameter :: wp = kind(1.0D0) 
            integer :: fileUnit, count, receivedPos, workingPos, ios, arraySize
            integer :: i, splitIndex, starIndex, currentIndex, repeatCount
            character(len=256) :: line, varName, ValueStr, lineTrimmed, token
            logical :: isReal
            real(wp), allocatable :: realArVal(:) 
            integer, allocatable  :: intArVal(:)
            
                !Force the correct values
                arraySize = 0
                currentIndex = 1
                workingPos = receivedPos
                
                !First, check the size of the array
                do 
                    call splitNameVar(varName, valueStr, line, index(line, "="))
                    call lineValueSize(valueStr, arraySize)
                    lineTrimmed = valueStr
                    read(fileUnit, '(A)', iostat=ios) line
                    if( index(line, "=") > 0 )exit
                    if( index(line, "/") > 0 )exit
                    workingPos = workingPos + 1
                end do
                    !Allocate the correct size, based on it being an integer or real
                    if(isReal) then 
                        allocate(realArVal(ArraySize))
                        realArVal = 0.0_wp
                    end if
                    
                    if(.not. isReal) then 
                        allocate(intArVal(ArraySize))
                        intArVal = 0
                    end if
                    
                !Reset back to the original state. 
                call rewindFile(fileUnit, receivedPos, line)
               !Next go over each line
               do 
                    !Split it into variable, and varname if required
                    call splitNameVar(varName, valueStr, line, index(line, "="))
                    lineTrimmed = valueStr
                    
                    !Loop through the line, add the first number to the array, and trim until the komma
                    !Standard fortran I/O always has the same format
                    do 
                        if(lineTrimmed == "") exit
                        splitIndex = index(lineTrimmed, ",")
                        if (splitIndex > 0) then
                            token = adjustl(lineTrimmed(:splitIndex-1))
                            lineTrimmed = adjustl(lineTrimmed(splitIndex+1:))
                        else
                            token = adjustl(lineTrimmed)
                            lineTrimmed = ""
                        end if

                        !Check for n*value
                        starIndex = index(token, "*")
                        if (starIndex > 0) then
                            read(token(:starIndex-1), *, iostat=ios) repeatCount
                            valueStr = token(starIndex+1:)
                        else
                            repeatCount = 1
                            valueStr = token
                        end if
                        
                        !Read and assign values repeatCount times
                        do i = 1, repeatCount
                            if (allocated(intArVal)) then
                                read(valueStr, *, iostat=ios) intArVal(currentIndex)
                            else if (allocated(realArVal)) then
                                read(valueStr, *, iostat=ios) realArVal(currentIndex)
                            end if
                            currentIndex = currentIndex + 1
                        end do
                    end do
                    read(fileUnit, '(A)', iostat=ios) line
                    if( index(line, "=") > 0 )exit 
                    if( index(line, "/") > 0 )exit
                end do
                
                call rewindFile(fileUnit, receivedPos, line)
        end subroutine parseArray
    
        subroutine rewindFile(fileUnit, receivedPos, line)
            integer :: fileUnit, receivedPos, i
            character(len=*) :: line 
            rewind(fileUnit)
            do i = 1, receivedPos
                read(fileUnit, '(A)') line
            end do
        end subroutine rewindFile
        
        subroutine lineValueSize(valueStr, count)
            character(len=*), intent(in) :: valueStr
            integer, intent(out) :: count
            character(len=256) :: token, rest
            integer :: i, starPos, repeatCount, ios

            rest = trim(valueStr)

            do while (len_trim(rest) > 0)
                i = index(rest, ",")
                if (i > 0) then
                    token = adjustl(rest(:i-1))
                    rest = adjustl(rest(i+1:))
                else
                    token = adjustl(rest)
                    rest = ""
                end if

                starPos = index(token, "*")
                if (starPos > 0) then
                    read(token(:starPos-1), *, iostat=ios) repeatCount
                    if (ios /= 0) repeatCount = 1 
                    count = count + repeatCount
                else
                    count = count + 1
                end if
            end do
        end subroutine lineValueSize
    
        subroutine splitNameVar(varName, valueStr, line, eq_pos)
            character(len=*) :: varName, valueStr, line
            integer :: eq_pos    
            
            varName = adjustl(line(:eq_pos-1))
            varName = trim(varName)
            valueStr = line(eq_pos+1:)
            
        end subroutine splitnameVar
    
        subroutine fetchInputFile(this, srcFilePath, desFilePath)
            class(testToolsMgr) :: this
            integer :: in_unit, out_unit, ios
            logical :: exists = .false.
            character(len=*) :: srcFilePath
            character(len=*) :: desFilePath
            character(len=256) :: line
            
            in_unit = 1
            out_unit = 2
            inquire(file=srcFilePath, exist=exists)
            if(.not. exists) then
                write(*,*) "Error: Could not find source file, Aborting. . .", ios
                stop
            end if
            
            open(unit=in_unit, file=srcFilePath, status="old", action="read", iostat=ios)
            if (ios /= 0) then
                write(*,*) "Error: Problem when trying to use source file.  Aborting. . .", ios
                stop
            end if
            
            
            open(unit=out_unit, file=desFilePath, status="replace", action="write", iostat=ios)
            if (ios /= 0) then
                write(*,*) "Error: problem when trying to use destination file. Aborting. . .", ios
                close(in_unit)
                stop
            end if 
            
            do 
                read(in_unit, '(A)', iostat=ios) line 
                if(ios /= 0) exit
                write(out_unit, '(A)') trim(line) 
            end do 
            
            close(in_unit)
            close(out_unit)
        end subroutine fetchInputFile
        
        subroutine fetchOutFileIndex(this, srcFilePath, desFilePath, in_unit, out_unit)
            class(testToolsMgr) :: this
            integer :: in_unit, out_unit, ios
            logical :: exists = .false.
            character(len=*) :: srcFilePath
            character(len=*) :: desFilePath
            character(len=256) :: lineIn, lineOut
            
            in_unit = 1
            out_unit = 2
            inquire(file=srcFilePath, exist=exists)
            if(.not. exists) then
                write(*,*) "Error: Could not find source file, Aborting. . .", ios
                stop
            end if
            
            open(unit=in_unit, file=srcFilePath, status="old", action="read", iostat=ios)
            if (ios /= 0) then
                write(*,*) "Error: Problem when trying to use source file.  Aborting. . .", ios
                stop
            end if
            
            inquire(file=desFilePath, exist=exists)
            if(.not. exists) then
                write(*,*) "Could not find destination file, stopping. . ."
            end if
            
            open(unit=out_unit, file=desFilePath, status="old", action="read", iostat=ios)
            if (ios /= 0) then
                write(*,*) "Error: problem when trying to use destination file. Aborting. . .", ios
                close(in_unit)
                stop
            end if 
        end subroutine fetchOutFileIndex
end module test_helper
