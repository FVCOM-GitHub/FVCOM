module mod_interp_time
  use lims
  use mod_utils
  use mod_par
  use control
  use mod_prec
  use mod_ncll
  use mod_nctools
  USE MOD_TIME
  USE MOD_FORCE
  implicit none

  CHARACTER(LEN=120) :: FIN, FOUT
  CHARACTER(LEN=120) :: TIMESTEP

  TYPE(NCFILE), POINTER ::NCF_IN
  TYPE(NCFILE), POINTER ::NCF_OUT

  TYPE VARLIST
     TYPE(NCVAR), POINTER :: NEXT
     TYPE(NCVAR), POINTER :: PREV
     TYPE(NCVAR), POINTER :: CURR
  END TYPE VARLIST

  TYPE(VARLIST), ALLOCATABLE :: VARSLIST(:)

contains



  SUBROUTINE GET_COMMANDLINE(CVS_ID,CVS_Date,CVS_Name,CVS_Revision)
    use mod_sng


    character(len=*), INTENT(IN)::CVS_Id  ! [sng] CVS Identification
    character(len=*), INTENT(IN)::CVS_Date ! [sng] Date string
    character(len=*), INTENT(IN)::CVS_Name ! [sng] File name string
    character(len=*), INTENT(IN)::CVS_Revision ! [sng] File revision string

    character(len=*),parameter::nlc=char(0) ! [sng] NUL character = ASCII 0 = char(0)
    ! Command-line parsing
    character(80)::arg_val ! [sng] command-line argument value
    character(200)::cmd_ln ! [sng] command-line
    character(80)::opt_sng ! [sng] Option string
    character(2)::dsh_key ! [sng] command-line dash and switch
    character(200)::prg_ID ! [sng] Program ID

    integer::arg_idx ! [idx] Counting index
    integer::arg_nbr ! [nbr] Number of command-line arguments
    integer::opt_lng ! [nbr] Length of option

    ! Main code
    call ftn_strini(cmd_ln) ! [sng] sng(1:len)=NUL

    call ftn_cmd_ln_sng(cmd_ln) ! [sng] Re-construct command-line into single string
    call ftn_prg_ID_mk(CVS_Id,CVS_Revision,CVS_Date,prg_ID) ! [sng] Program ID

    arg_nbr=command_argument_count() ! [nbr] Number of command-line arguments

    if (arg_nbr .LE. 0 ) then
       if(MSR) WRITE(IPT,*) "You must specify an arugument:"
       if(MSR) Call MYHelpTxt
       call PSHUTDOWN
    end if

    arg_idx=1 ! [idx] Counting index
    do while (arg_idx <= arg_nbr)
       call ftn_getarg_wrp(arg_idx,arg_val) ! [sbr] Call getarg, increment arg_idx
       dsh_key=arg_val(1:2) ! [sng] First two characters of option
       if (dsh_key == "--") then
          opt_lng=ftn_opt_lng_get(arg_val) ! [nbr] Length of option
          if (opt_lng <= 0) then
             if(MSR) write(IPT,*) "Long option has no name"
             call PSHUTDOWN
          end if

          opt_sng=arg_val(3:2+opt_lng) ! [sng] Option string
          if (dbg_lvl >= dbg_io) then
             if(MSR) write (6,"(5a,i3)") prg_nm(1:ftn_strlen(prg_nm)), &
                  ": DEBUG Double hyphen indicates multi-character option: ", &
                  "opt_sng = ",opt_sng(1:ftn_strlen(opt_sng)),", opt_lng = ",opt_lng
          end if
          if (opt_sng == "dbg" .or. opt_sng == "dbg_lvl" ) then
             call ftn_arg_get(arg_idx,arg_val,dbg_lvl) ! [enm] Debugging level

             !          else if (opt_sng == "dbg_par" .or.opt_sng == "Dbg_Par"&
             !               & .or.opt_sng == "DBG_PAR") then

             !             dbg_par = .true.

          else if (opt_sng == "INPUT" .or.opt_sng == "input"&
               & .or.opt_sng == "Input") then

             call ftn_arg_get(arg_idx,arg_val,FIN) ! [sng] Input file
             FIN=FIN(1:ftn_strlen(FIN))
             ! Convert back to a fortran string!

          else if (opt_sng == "OUTPUT" .or.opt_sng == "output"&
               & .or.opt_sng == "Output") then

             call ftn_arg_get(arg_idx,arg_val,FOUT) ! [sng] Input file
             FOUT=FOUT(1:ftn_strlen(FOUT))
             ! Convert back to a fortran string!

          else if (opt_sng == "START" .or.opt_sng == "start"&
               & .or.opt_sng == "Start") then

             call ftn_arg_get(arg_idx,arg_val,Start_Date) ! [sng] Input file
             Start_Date=Start_Date(1:ftn_strlen(Start_Date))
             ! Convert back to a fortran string!

          else if (opt_sng == "END" .or.opt_sng == "end"&
               & .or.opt_sng == "End") then

             call ftn_arg_get(arg_idx,arg_val,End_Date) ! [sng] Input file
             End_Date=End_Date(1:ftn_strlen(End_Date))
             ! Convert back to a fortran string!

          else if (opt_sng == "TIMEZONE" .or.opt_sng == "timezone"&
               & .or.opt_sng == "TimeZone") then

             call ftn_arg_get(arg_idx,arg_val,TIMEZONE) ! [sng] Input file
             TIMEZONE=TIMEZONE(1:ftn_strlen(TIMEZONE))
             ! Convert back to a fortran string

          else if (opt_sng == "TIMESTEP" .or.opt_sng == "timestep"&
               & .or.opt_sng == "TimeStep") then

             call ftn_arg_get(arg_idx,arg_val,TIMESTEP) ! [sng] Input file
             TIMESTEP=TIMESTEP(1:ftn_strlen(TIMESTEP))
             ! Convert back to a fortran string!


          else if (opt_sng == "help" .or.opt_sng == "HELP" .or. opt_sng&
               & == "Help") then

             if(MSR) call MYHelpTxt
             call PSHUTDOWN

          else ! Option not recognized
             arg_idx=arg_idx-1 ! [idx] Counting index
             if(MSR) call ftn_getarg_err(arg_idx,arg_val) ! [sbr] Error handler for getarg()
          endif ! endif option is recognized
          ! Jump to top of while loop
          cycle ! C, F77, and F90 use "continue", "goto", and "cycle"
       endif ! endif long option

       if (dsh_key == "-V" .or.dsh_key == "-v" ) then

          if(MSR) write(IPT,*) prg_id
          call PSHUTDOWN

       else if (dsh_key == "-H" .or.dsh_key == "-h" ) then

          if(MSR) Call MYHelpTxt
          Call PSHUTDOWN

       else ! Option not recognized
          arg_idx=arg_idx-1 ! [idx] Counting index
          if(MSR) call ftn_getarg_err(arg_idx,arg_val) ! [sbr] Error handler for getarg()
       endif ! endif arg_val


    end do ! end while (arg_idx <= arg_nbr)

    CALL dbg_init(IPT_BASE,.false.)

  END SUBROUTINE GET_COMMANDLINE

  SUBROUTINE MYHELPTXT
    IMPLICIT NONE


    write(IPT,*) "! ARGUMENTS FOR Interp_Time:"
    write(IPT,*) "! START  = (a date or time string)"
    write(IPT,*) "! END    = (a date or time string)"
    write(IPT,*) "! TIMEZONE = (the timezone to use)"
    write(IPT,*) "! TIMESTEP = (time step of interpolated data)"
    write(IPT,*) "! INPUT  = The/source/of/data.nc"
    write(IPT,*) "! OUTPUT = The/result/of/interped.nc"
    write(IPT,*) "! "
    WRITE(IPT,*) "! NOTES: Do not run this program in parallel!"


  END SUBROUTINE MYHELPTXT


  SUBROUTINE INTERP_TIME
    IMPLICIT NONE
    LOGICAL FOUND
    integer status

    NCF_IN => NEW_FILE()
    NCF_IN%FNAME=trim(FIN)
    NCF_IN%WRITABLE=.False.
    ! OPEN THE FILE AND LOAD METADATA       

    Call NC_OPEN(NCF_IN)
    CALL NC_LOAD(NCF_IN)

    CALL SET_TIME

    Call Setup_Interp_Memory(NCF_IN)

    NCF_OUT => COPY_FILE(NCF_IN)
    NCF_OUT%FNAME=trim(FOUT)
    NCF_OUT%WRITABLE=.TRUE.
    NCF_OUT%FTIME=>NEW_FTIME()

    IF(DBG_SET(DBG_SBRIO)) THEN
       write(ipt,*) "==============================================="
       write(ipt,*) "============ Debuging Setup ======="
       write(ipt,*) "==============================================="
       CALL PRINT_FILE(NCF_IN)
       CALL PRINT_VAR_LIST(NCF_IN)
       CALL PRINT_DIM_LIST(NCF_IN)
       
       
       write(ipt,*) "==============================================="
       write(ipt,*) "==============================================="
       
       
       call PRINT_FILE(NCF_OUT)
       CALL PRINT_VAR_LIST(NCF_OUT)
       CALL PRINT_DIM_LIST(NCF_OUT)

       write(ipt,*) "==============================================="
       write(ipt,*) "==============================================="
    END IF

    
    CALL STATIC_VARIABLES(NCF_IN,NCF_OUT)

    CALL MAIN_LOOP



  END SUBROUTINE INTERP_TIME

  SUBROUTINE SET_TIME
    USE MOD_SET_TIME
    IMPLICIT NONE
    integer status, NTIMES
    INTEGER(ITIME) :: dummy
    CHARACTER(LEN=4) :: FLAG
    LOGICAL TEST, FOUND
    TYPE(TIME) TIMETEST
    TYPE(NCDIM), POINTER :: DIM

    ! SET DEFAULT TO TRUE FOR REAL TIME MODEL
    use_real_world_time = .TRUE.
    ! TEST FOR IDEALIZED MODEL CASE
    if (timezone == 'none' .or. timezone == "NONE" .or.&
         & timezone == "None") use_real_world_time = .FALSE.
    
    ! CHECK FOR VALID TIME ZONE
    TEST = IS_VALID_TIMEZONE(timezone)
    IF(.not. TEST) call fatal_error("You selected an invalid time zone: "&
         &//trim(timezone),"Time Zones must be CAPITALS",&
         & "see mod_time.F for a list of valid time_zones")
    

    IF (USE_REAL_WORLD_TIME) THEN
       StartTime = read_datetime(START_DATE,DATE_FORMAT,TIMEZONE,status)
       if (status==1) & 
            & Call Fatal_Error("Could not read the date string START_DATE: ", trim(START_DATE))

       ! GET THE END TIME
       EndTime = READ_DateTime(END_DATE,DATE_FORMAT,TIMEZONE,status)
       if (status==1) &
            & Call Fatal_Error("Could not read the date string END_DATE:", trim(END_DATE))

    ELSE

       CALL IDEAL_TIME_STRING2TIME(START_DATE,FLAG,StartTime,dummy)
       IF(FLAG == 'step') CALL FATAL_ERROR&
            &("You must specify a time, not a step, for this restart file", &
            & "The Step will be set by the old restart file...")

       CALL IDEAL_TIME_STRING2TIME(END_DATE,FLAG,EndTime,dummy)
       IF(FLAG == 'step') CALL FATAL_ERROR&
            &("You must specify a time, not a step, for this restart file", &
            & "The Step will be set by the old restart file...")


    END IF


    CALL IDEAL_TIME_STRING2TIME(TIMESTEP,FLAG,IMDTI,dummy)
    IF(FLAG == 'step') CALL FATAL_ERROR&
         &("You must specify a time, not a step, for this restart file", &
         & "The Step will be set by the old restart file...")
    

    CALL PRINT_REAL_TIME(STARTTIME,IPT,"READ: START_DATE",TIMEZONE)
    CALL PRINT_REAL_TIME(ENDTIME,IPT,"READ: END_DATE",TIMEZONE)
    CALL PRINT_REAL_TIME(IMDTI,IPT,"READ: TIMESTEP","none")


    ! SANITY ECK
    if(StartTime .GT. EndTime) &
         & Call Fatal_Error("Runfile Start_Date exceeds or equal to End_Date")

    IF(IMDTI .LE. ZEROTIME) CALL FATAL_ERROR("CAN NOT SET TIMESTEP .LE. ZERO")

    
    ! LOAD TIME AND CHECK TO MAKE SURE THE TIME RANGE IS VALID
    DIM => FIND_UNLIMITED(NCF_IN,FOUND)
    IF (.not. FOUND) CALL FATAL_ERROR("INPUT FILE HAS NO UNLIMITED DIMENSION!")
    NTIMES = DIM%DIM

    TIMETEST = GET_FILE_TIME(NCF_IN,1)
    
    IF(TIMETEST > STARTTIME) THEN
       CALL PRINT_REAL_TIME(TIMETEST,IPT,"FIRST TIME IN FILE",TIMEZONE)

       CALL FATAL_ERROR("The start time preceeds the input file")
    END IF

    TIMETEST = GET_FILE_TIME(NCF_IN,NTIMES)
    
    IF(TIMETEST < ENDTIME) THEN
       CALL PRINT_REAL_TIME(TIMETEST,IPT,"LAST TIME IN FILE",TIMEZONE)
       
       CALL FATAL_ERROR ("The end time exceeds the input file")
    END IF

  END SUBROUTINE SET_TIME

  SUBROUTINE Setup_Interp_Memory(LIST)
    IMPLICIT NONE
    INTEGER :: CNT,STATUS
    LOGICAL FOUND
    TYPE(NCFILE), POINTER :: LIST
    TYPE(NCVARP), POINTER    :: CURRENT, PREVIOUS
    TYPE(NCDIM), POINTER :: DIM
    TYPE(NCVAR), POINTER :: VAR

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START: Setup_Interp_Memory"
    
    PREVIOUS => LIST%VARS
    CURRENT  => PREVIOUS%NEXT

    CNT=COUNT_UNLIMITED_VARS(LIST)

    ALLOCATE(VARSLIST(CNT),STAT=STATUS)
    IF (STATUS/=0) CALL FATAL_ERROR("CAN NOT ALLOCATE VARSLIST")

    CNT=0

    DO 
       IF(.NOT. ASSOCIATED(CURRENT)) RETURN !END OF LIST
       
       IF (HAS_UNLIMITED(CURRENT%VAR)) THEN
          ! VARIABLE IS TIME DEPENDENT
          ! ALOCATE MEMORY
          CALL ALLOC_VAR(CURRENT%VAR)

          CNT=CNT+1

          VAR =>REFERENCE_VAR(CURRENT%VAR)
          CALL ALLOC_VAR(VAR)
          VARSLIST(CNT)%PREV => VAR
          NULLIFY(VAR)

          VAR =>REFERENCE_VAR(CURRENT%VAR)
          CALL ALLOC_VAR(VAR)
          VARSLIST(CNT)%NEXT => VAR
          NULLIFY(VAR)

          VARSLIST(CNT)%CURR => CURRENT%VAR


       ELSE
          ! VARIABLE IS NOT TIME DEPENDENT
       
          CALL ALLOC_VAR(CURRENT%VAR)
       END IF


       PREVIOUS => PREVIOUS%NEXT
       CURRENT  => CURRENT%NEXT
    END DO

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END: Setup_Interp_Memory"

  END SUBROUTINE Setup_Interp_Memory


  SUBROUTINE MAIN_LOOP
    IMPLICIT NONE
    
    INTEGER :: I,NUNLIM
    INTEGER :: CNT,STATUS
    LOGICAL FOUND
    TYPE(NCFTIME), POINTER :: FTM
    TYPE(NCVAR), POINTER :: CURR, NEXT, PREV

    FTM => NCF_IN%FTIME

    CNT = size(varslist)
    IntTime = StartTime
    DO WHILE (IntTime .LE. ENDTIME)
          
       DO I=1,CNT
          CURR=> VARSLIST(I)%CURR 
          NEXT=> VARSLIST(I)%NEXT
          PREV=> VARSLIST(I)%PREV 


          CALL UPDATE_VAR_BRACKET(NCF_IN,PREV,NEXT,IntTime,status)
          IF (STATUS /= 0) THEN
             CALL PRINT_TIME(INTTIME,IPT,"WTF")
             CALL FATAL_ERROR("COULD NOT UPATE THE BRACKET: BOUNDS EXCEEDED?")
          end if
          

          CALL CALCULATE_INTERP(FTM,PREV,CURR,NEXT)

          
       END DO
       
       
       NCF_OUT%FTIME%PREV_STKCNT = NCF_OUT%FTIME%NEXT_STKCNT
       NCF_OUT%FTIME%NEXT_STKCNT = NCF_OUT%FTIME%NEXT_STKCNT + 1
       CALL NC_WRITE_FILE(NCF_OUT)
       
       INTTIME = INTTIME + IMDTI
    END DO


  END SUBROUTINE MAIN_LOOP


  SUBROUTINE STATIC_VARIABLES(IN,OUT)
    IMPLICIT NONE
    INTEGER :: CNT,STATUS
    LOGICAL FOUND
    TYPE(NCFILE), POINTER :: IN,OUT
    TYPE(NCVARP), POINTER    :: CURRENT, PREVIOUS
    TYPE(NCDIM), POINTER :: DIM
    TYPE(NCVAR), POINTER :: VAR

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START: STATIC VARIALBES"
    
    PREVIOUS => IN%VARS
    CURRENT  => PREVIOUS%NEXT

    DO 
       IF(.NOT. ASSOCIATED(CURRENT)) exit !END OF LIST
       
       IF (.NOT. HAS_UNLIMITED(CURRENT%VAR)) THEN
          
          CALL NC_READ_VAR(CURRENT%VAR)
                    
       END IF

       PREVIOUS => PREVIOUS%NEXT
       CURRENT  => CURRENT%NEXT
    END DO

    CALL NC_WRITE_FILE(OUT)

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END: Setup_Interp_Memory"


  END SUBROUTINE STATIC_VARIABLES

  SUBROUTINE CALCULATE_INTERP(FTM,PREV,CURR,NEXT)
    IMPLICIT NONE

    TYPE(NCVAR), POINTER :: PREV,CURR,NEXT
    TYPE(NCFTIME), POINTER :: FTM

    INTEGER :: CNT,STATUS
    LOGICAL FOUND
    TYPE(NCDIM), POINTER :: DIM
    TYPE(NCVAR), POINTER :: VAR
    TYPE(NCVAR), POINTER :: VAR_TMP

    INTEGER :: Ndims
    INTEGER, POINTER :: DIMS(:)

    NULLIFY(DIM,DIMS)

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START: CALCULATE_INTERP"

    IF(.not. ASSOCIATED(PREV))CALL FATAL_ERROR&
            & ("CALCULATE_INTERP: Unassociated PREV!")

    IF(.not. ASSOCIATED(CURR))CALL FATAL_ERROR&
            & ("CALCULATE_INTERP: Unassociated CURR!")

    IF(.not. ASSOCIATED(NEXT))CALL FATAL_ERROR&
            & ("CALCULATE_INTERP: Unassociated NEXT!")
    
    IF(.not. ASSOCIATED(FTM))CALL FATAL_ERROR&
            & ("CALCULATE_INTERP: Unassociated FTM!")


    ! GET RID OF TIME AND SINGLETON DIMENSIONS WHEN ALLOCATING SPACE
    VAR_TMP => COPY_VAR(CURR)

    DIM => FIND_UNLIMITED(VAR_TMP,FOUND)
    IF (FOUND) THEN
       CALL DELETE_DIM_LINK(VAR_TMP,DIM%DIMID,FOUND)
    END IF

    DIMS => MEM_DIMS(VAR_TMP)

    IF(.not. Associated(DIMS)) CALL FATAL_ERROR("ALLOC_VAR: Could not allocate Dims?")
    Ndims=size(DIMS)

    CALL KILL_VAR(VAR_TMP)
    
    IF(DBG_SET(DBG_SBRIO)) WRITE(IPT,*) "RUNNING:"

    select case(CURR%XTYPE)
    case(NF90_CHAR)
       IF(DBG_SET(DBG_SBRIO)) THEN
          WRITE(IPT,*) "XYTPE       :: CHAR"
          WRITE(IPT,*) "DIMS       ::",dims
       END IF

       SELECT CASE(NDIMS)
       CASE(2)
          CURR%VEC_CHR=PREV%VEC_CHR
       CASE(1)
          CURR%SCL_CHR=PREV%SCL_CHR
       CASE(0)
          CALL FATAL_ERROR("Unsupported Character data dimension")
       END SELECT

    case(NF90_BYTE)
       WRITE(IPT,*) "XYTPE       :: BYTE"
       CALL FATAL_ERROR("No Byte Type Available")

    case(NF90_SHORT)
       WRITE(IPT,*) "XYTPE       :: SHORT"
       CALL FATAL_ERROR("No Short Type Available")

    case(NF90_INT)
       IF(DBG_SET(DBG_SBRIO)) THEN
          WRITE(IPT,*) "XYTPE       :: INT"
          WRITE(IPT,*) "dims       ::",dims
       END IF


       SELECT CASE(NDIMS)
       CASE(3)
          CURR%CUB_INT= PREV%CUB_INT
       CASE(2)
          CURR%ARR_INT= PREV%ARR_INT
       CASE(1)
          CURR%VEC_INT= PREV%VEC_INT
       CASE(0)
          CURR%SCL_INT= PREV%SCL_INT
       CASE DEFAULT
          CALL FATAL_ERROR("Unsupported Integer data dimension")
       END SELECT


    case(NF90_FLOAT)
       IF(DBG_SET(DBG_SBRIO)) THEN
          WRITE(IPT,*) "XYTPE       :: FLOAT"
          WRITE(IPT,*) "dims       ::",dims
       END IF

       SELECT CASE(NDIMS)
       CASE(3)
          CURR%CUB_FLT= FTM%PREV_WGHT *PREV%CUB_FLT + FTM%NEXT_WGHT *NEXT%CUB_FLT
       CASE(2)
          CURR%ARR_FLT= FTM%PREV_WGHT *PREV%ARR_FLT + FTM%NEXT_WGHT *NEXT%ARR_FLT
       CASE(1)
          CURR%VEC_FLT= FTM%PREV_WGHT *PREV%VEC_FLT + FTM%NEXT_WGHT *NEXT%VEC_FLT
       CASE(0)
          CURR%SCL_FLT= FTM%PREV_WGHT *PREV%SCL_FLT + FTM%NEXT_WGHT *NEXT%SCL_FLT
       CASE DEFAULT
          CALL FATAL_ERROR("Unsupported Integer data dimension")
       END SELECT
       
    case(NF90_DOUBLE)
       IF(DBG_SET(DBG_SBRIO)) THEN
          WRITE(IPT,*) "XYTPE       :: DOUBLE"
          WRITE(IPT,*) "dims       ::",dims
       END IF

       SELECT CASE(NDIMS)
       CASE(3)
          CURR%CUB_DBL= FTM%PREV_WGHT *PREV%CUB_DBL + FTM%NEXT_WGHT *NEXT%CUB_DBL
       CASE(2)
          CURR%ARR_DBL= FTM%PREV_WGHT *PREV%ARR_DBL + FTM%NEXT_WGHT *NEXT%ARR_DBL
       CASE(1)
          CURR%VEC_DBL= FTM%PREV_WGHT *PREV%VEC_DBL + FTM%NEXT_WGHT *NEXT%VEC_DBL
       CASE(0)
          CURR%SCL_DBL= FTM%PREV_WGHT *PREV%SCL_DBL + FTM%NEXT_WGHT *NEXT%SCL_DBL

       CASE DEFAULT
          CALL FATAL_ERROR("Unsupported Integer data dimension")
       END SELECT
       
    END select
    

    DEALLOCATE(DIMS)

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END: CALCULATE_INTERP"


  END SUBROUTINE CALCULATE_INTERP


end module mod_interp_time
