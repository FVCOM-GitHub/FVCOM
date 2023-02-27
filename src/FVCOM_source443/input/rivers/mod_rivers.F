module mod_rivers
  use lims
  use mod_utils
  use mod_par
  use mod_time
  use control
  use mod_prec
  use mod_nctools
  implicit none
  
  
  REAL(SP), ALLOCATABLE :: RIVER_DISTROBUTION(:,:)

  CHARACTER(len=80) :: RUN_FILE
  CHARACTER(len=80), PARAMETER :: FOUT = 'River_data.nc'

  CHARACTER(len=80) :: river_data_file
  CHARACTER(len=80) :: river_out_file
  CHARACTER(len=80) :: title
  CHARACTER(len=80) :: website
  CHARACTER(len=80) :: history

  integer, parameter :: namelen = 80

  INTEGER :: NRIVERS
  INTEGER :: NTIMES
  CHARACTER(LEN=80) :: INFLOW_TYPE,POINT_ST_TYPE
  CHARACTER(LEN=namelen), ALLOCATABLE :: RIVER_NAMES(:)



  INTEGER,  ALLOCATABLE :: NLOC(:)
  TYPE(TIME), ALLOCATABLE :: TTIME(:)
  REAL(SP), ALLOCATABLE :: QDIS(:,:)
  REAL(SP), ALLOCATABLE :: QTEMP(:,:)
  REAL(SP), ALLOCATABLE :: QSALT(:,:)
  
  INTEGER :: SIGMA_LAYERS

  NAMELIST /NML_RIVER2FVCOM/               &
       & INPUT_DIR,                        &
       & OUTPUT_DIR,                       &
       & RIVER_DATA_FILE,                  &
       & SIGMA_LAYERS,                     &
       & TIMEZONE,                         &
       & DATE_FORMAT,                      &
       & START_DATE,                       &
       & TITLE,                            &
       & WEBSITE,                          &
       & HISTORY

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

       if(MSR) Call HelpTxt
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

          else if (opt_sng == "filename" .or.opt_sng == "FILENAME"&
               & .or.opt_sng == "FileName") then

             call ftn_arg_get(arg_idx,arg_val,RUN_FILE) ! [sng] Input file
             RUN_FILE=RUN_FILE(1:ftn_strlen(RUN_FILE))
             ! Convert back to a fortran string!

!          else if (opt_sng == "Create_NameList" .or.opt_sng == "create_namelist"&
!               & .or.opt_sng == "CREATE_NAMELIST") then
!
!             call ftn_arg_get(arg_idx,arg_val,NAMELIST_NAME)
!             NAMELIST_NAME = NAMELIST_NAME(1:ftn_strlen(NAMELIST_NAME))
!
!             BLANK_NAMELIST = .true.

             
          else if (opt_sng == "help" .or.opt_sng == "HELP" .or. opt_sng&
               & == "Help") then

             if(MSR) call HelpTxt

             call PSHUTDOWN
!!$   THIS DOES NOT SEEM PRACTICAL - MODIFY THE RUN FILE INSTEAD
!!$          else if (opt_sng == "CrashRestart") then
!!$             call ftn_arg_get(arg_idx,arg_val,CrashRestart) ! [lgc] Logical

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
          if(MSR) call HelpTxt

          Call PSHUTDOWN

       else ! Option not recognized
          arg_idx=arg_idx-1 ! [idx] Counting index
          if(MSR) call ftn_getarg_err(arg_idx,arg_val) ! [sbr] Error handler for getarg()
       endif ! endif arg_val


    end do ! end while (arg_idx <= arg_nbr)

    CALL dbg_init(IPT_BASE,.false.)

  END SUBROUTINE GET_COMMANDLINE


  SUBROUTINE HelpTxt
    IMPLICIT NONE
    
    WRITE(IPT,*) "Add better help here!"
    WRITE(IPT,*) "! OPTIONS:"
    WRITE(IPT,*) "! --filename=XXX"
    WRITE(IPT,*) "!   The namelist runfile for the program! "
    WRITE(IPT,*) "!   "
    WRITE(IPT,*) "!   Namelist OPTIONS: "
    WRITE(IPT,*) "!   INPUT_DIR, (Required)"
    WRITE(IPT,*) "!   OUTPUT_DIR, (Required)"
    WRITE(IPT,*) "!   RIVER_DATA_FILE, (Required)"
    WRITE(IPT,*) "!   START_DATE, (Required)"
    WRITE(IPT,*) "!   TIMEZONE, (Required)"
    WRITE(IPT,*) "!   DATE_FORMAT,(Required)"
    WRITE(IPT,*) "!   TITLE, (OPTIONAL)"
    WRITE(IPT,*) "!   WEBSITE, (OPTIONAL)"
    WRITE(IPT,*) "!   HISTORY, (OPTIONAL)"
    WRITE(IPT,*) "!  "
    WRITE(IPT,*) "!  EXAMPLE NAMELSIT:"
    write(UNIT=IPT,NML=NML_RIVER2FVCOM)
    
    WRITE(IPT,*) "! NOTES: Do not run this program in parallel!"

  END SUBROUTINE HelpTxt


  SUBROUTINE READ_NAMELIST
    IMPLICIT NONE
    integer :: ios, i
    Character(Len=120):: FNAME
    if(DBG_SET(dbg_sbr)) &
         & write(IPT,*) "Subroutine Begins: GET_RUNFILE;"

    RIVER_DATA_FILE = "NONE"
    RIVER_OUT_FILE = "NONE"

    TIMEZONE = "NONE"
    DATE_FORMAT = "NONE"

    START_DATE = "NONE"
    
    ! SET DEFAULT VALUES:
    HISTORY = "DEFAULT VALUE : IT IS 1066 AND I AM WILLIAM THE CONQUORER!"

    WEBSITE = "DEFAULT VALUE : USGS.GOV"

    TITLE = "DEFAULT VALUE : LORD OF FVCOM_DOM"

    FNAME = RUN_FILE
    
    if(DBG_SET(dbg_io)) &
         & write(IPT,*) "Read_Name_List: File: ",trim(FNAME)
    
    CALL FOPEN(NMLUNIT,trim(FNAME),'cfr')
    
    !READ NAME LIST FILE
    
    ! Read IO Information
    READ(UNIT=NMLUNIT, NML=NML_RIVER2FVCOM,IOSTAT=ios)
    if(ios .NE. 0 ) CALL FATAL_ERROR&
         &("Can Not Read NameList NML_RIVER2FVCOM from file: "//trim(FNAME))
    REWIND(NMLUNIT)

    if(DBG_SET(dbg_scl)) &
         & write(IPT,*) "Read_Name_List:"

    KBM1 = SIGMA_LAYERS

    if(DBG_SET(dbg_scl)) &
         & write(UNIT=IPT,NML=NML_RIVER2FVCOM)


  END SUBROUTINE READ_NAMELIST



  SUBROUTINE SET_TIME
    USE MOD_SET_TIME
    IMPLICIT NONE
    integer status
    REAL(SP) :: ST
    integer(itime) :: dummy
    CHARACTER(LEN=4) :: flag
    

    IF (USE_REAL_WORLD_TIME) THEN
       StartTime = READ_DATETIME(START_DATE,DATE_FORMAT,TIMEZONE,status)
       if (status /= 1) &
            & Call Fatal_Error("Could not read the date string START_DATE: ",&
            & trim(START_DATE))
       CALL PRINT_REAL_TIME(StartTime,IPT,"FIRST RIVER TIME")
    ELSE

       CALL IDEAL_TIME_STRING2TIME(START_DATE,FLAG,StartTIME,dummy)
       if (FLAG /= 'time') CALL FATAL_ERROR &
            & ('you can not specify a start cycle for the river file:',&
            & "uses 'days = XXXX' or 'seconds = XXXX' ")
       CALL PRINT_TIME(StartTime,IPT,"FIRST RIVER TIME")
          
    END IF

    

  END SUBROUTINE SET_TIME



  SUBROUTINE LOAD_RIVER_DATA
    IMPLICIT NONE
    integer, parameter :: inriv = 201
    INTEGER :: I, J,idx, b,e,IOS
    character(len=160) :: pathnfile,STRING
    real(sp) :: rdummy
    real(DP) :: HOUR, DAY
    REAL(SP), ALLOCATABLE :: MEANQ(:),MEANT(:), MEANS(:)
    REAL(SP), ALLOCATABLE :: MAXQ(:),MAXT(:), MAXS(:)
    REAL(SP), ALLOCATABLE :: MINQ(:),MINT(:), MINS(:)
   
    pathnfile = trim(INPUT_DIR)//trim(RIVER_DATA_FILE)
    CALL FOPEN(INRIV,TRIM(pathnfile),'cfr')

    READ(INRIV,'(A4,2X,A10)') INFLOW_TYPE,POINT_ST_TYPE
    IF(MSR)WRITE(IPT,*) 'River Inflow Information'
    IF(MSR)WRITE(IPT,*) 'INFLOW_TYPE==',TRIM(INFLOW_TYPE)
    IF(MSR)WRITE(IPT,*) 'POINT_ST_TYPE==',TRIM(POINT_ST_TYPE)
    
    IF(INFLOW_TYPE /= 'edge' .AND. INFLOW_TYPE /= 'node') THEN
       CALL FATAL_ERROR("INFLOW TYPE MUST BE 'edge' or 'node'")
    END IF
    
    IF(POINT_ST_TYPE /= 'calculated' .AND. POINT_ST_TYPE /= 'specified') THEN
       CALL FATAL_ERROR("POINT_ST TYPE NOT CORRECT","SHOULD BE 'calculated' or 'specified'")
    END IF
    
    READ(INRIV,*) NRIVERS
    
    WRITE(IPT,*) "NUMBER OF RIVERS: ", NRIVERS

    WRITE(IPT,*) "==================================="
    WRITE(IPT,*) "NOW READING RIVER LOCATION AND RIVER NAME:" 
    WRITE(IPT,*) "==================================="
   
    ALLOCATE(NLOC(NRIVERS))
    ALLOCATE(RIVER_NAMES(NRIVERS))
    ALLOCATE(RIVER_DISTROBUTION(NRIVERS,MAX_LAYERS))
    RIVER_DISTROBUTION =-99.0_SP
    DO I=1,NRIVERS
       READ(INRIV,'(a)') string
       write(ipt,*) "READ STRING: '"//TRIM(string)//"'"

       idx = index(trim(string),":")

       if(idx == 0) call fatal_error&
            & ('You must add the river name to the river file:',& 
            & " 'node : name' ", " '31503 : myrivername' ")

       b = 1
       e = idx -1
       read(string(b:e),*) NLOC(I)

       b = idx + 1
       e = len_trim(string)
       RIVER_NAMES(I) = string(b:e)
       RIVER_NAMES(I) = adjustl( RIVER_NAMES(I))       
       IF (LEN_TRIM(RIVER_NAMES(I)) == 0 ) call fatal_error&
            & ('You must add the river name to the river file:',& 
            & " 'node : name' ", " '31503 : myrivername' ")

       WRITE(IPT,*) "RIVER LOCATION IS: ", NLOC(I)
       WRITE(IPT,*) "RIVER NAME IS: '"//TRIM(RIVER_NAMES(I))//"'"
       WRITE(IPT,*) "==================================="
    END DO
    WRITE(IPT,*) "==================================="
    WRITE(IPT,*) "==================================="
    WRITE(IPT,*) "==================================="
    

    ! READ DISTROBUTION STRING
    DO I = 1,NRIVERS
       READ(INRIV,*,IOSTAT=IOS) J,(RIVER_DISTROBUTION(I,J),J = 1,KBM1)
       IF(IOS /=0) CALL FATAL_ERROR("COULD NOT READ RIVER DISTRIBUTION: KBM1 DOES NOT MATCH FILE?")
    END DO
    
     !----Read in Time Dependent DataSets (DQDIS,DSDIS,DTDIS)------------------------!
     !
    READ(INRIV,*) NTIMES

     WRITE(IPT,*) "NUMBER OF TIMES: ", NTIMES

    ALLOCATE(TTIME(NTIMES))

    ALLOCATE(QDIS(NRIVERS,NTIMES))
    ALLOCATE(QTEMP(NRIVERS,NTIMES))
    ALLOCATE(QSALT(NRIVERS,NTIMES))

    DO I = 1, NTIMES
           READ(INRIV,*) HOUR
           READ(INRIV,*) (QDIS(J,I),J = 1,NRIVERS)
           READ(INRIV,*) (QTEMP(J,I),J = 1,NRIVERS)
           READ(INRIV,*) (QSALT(J,I),J = 1,NRIVERS)

           DAY = REAL(HOUR,DP) / 24.0_DP

           TTIME(I)= DAYS2TIME(DAY) + StartTime

    END DO

    CLOSE(INRIV)


    ALLOCATE(MAXQ(NRIVERS),MAXS(NRIVERS),MAXT(NRIVERS))
    ALLOCATE(MINQ(NRIVERS),MINS(NRIVERS),MINT(NRIVERS))
    ALLOCATE(MEANQ(NRIVERS),MEANS(NRIVERS),MEANT(NRIVERS))


    MAXQ = MAXVAL(QDIS,2)
    MAXS = MAXVAL(QSALT,2)
    MAXT = MAXVAL(QTEMP,2)

    MINQ = MINVAL(QDIS,2)
    MINS = MINVAL(QSALT,2)
    MINT = MINVAL(QTEMP,2)

    MEANQ = SUM(QDIS,2)/REAL(NTIMES,SP)
    MEANS = SUM(QSALT,2)/REAL(NTIMES,SP)
    MEANT = SUM(QTEMP,2)/REAL(NTIMES,SP)

    WRITE(IPT,*) "==================================="
    WRITE(IPT,*) "REPORT RESULTS:"
    WRITE(IPT,*) "==================================="
    
    CALL PRINT_REAL_TIME(TTIME(1),IPT,"TIME OF FIRST DATA")

    CALL PRINT_REAL_TIME(TTIME(NTIMES),IPT,"TIME OF LAST DATA")

    CALL FOPEN(rivernmlunit,"RIVERS_NAMELIST.nml",'ofr')

    DO I = 1,NRIVERS
       
       write(IPT,*) "RIVER NAME: "//TRIM(RIVER_NAMES(I))//";"
       WRITE(IPT,*) "! FLUX MIN/MEAN/MAX :",minq(I),meanq(I),maxq(i)

       write(IPT,*) "! TEMP MIN/MEAN/MAX :",minT(I),meanT(I),maxT(i)

       write(IPT,*) "! SALT MIN/MEAN/MAX :",minS(I),meanS(I),maxS(i)

       WRITE(IPT,*) "==================================="


       !SET VALUES IN NAMELIST OBJECT
       RIVER_NAME = "'"//TRIM(RIVER_NAMES(I))//"'"
       RIVER_FILE = "'"//trim(FOUT)//"'"
       RIVER_GRID_LOCATION = NLOC(I)
# if defined(RIVER_FLOAT)
       RIVER_VERTICAL_DISTRIBUTION = RIVER_DISTROBUTION(I,:)
# else
       RIVER_VERTICAL_DISTRIBUTION = "Set river distribution!"
# endif
       write(rivernmlunit,NML=NML_RIVER)
    END DO
    
    CLOSE(rivernmlunit)
    !
    !--REPORT RESULTS--------------------------------------------------------------!
    !

  END SUBROUTINE LOAD_RIVER_DATA


! NOT IT IS NOT EASY TO CREATE THE NAMES ARRAY USING MOD_NCLL
  !/MOD_NCTOOLS SO I AM DOING IT THE OLD FASHIONED WAY...
  SUBROUTINE DUMP_OUTFILE
    USE MOD_TIME
    IMPLICIT NONE

    TYPE(NCFILE), POINTER :: NCF
    TYPE(NCVAR), POINTER :: VAR,var1,var2
    TYPE(NCATT), POINTER :: ATT

    TYPE(NCDIM), POINTER :: DIM_time
    TYPE(NCDIM), POINTER :: DIM_DateStrLen
    TYPE(NCDIM), POINTER :: DIM_namelen
    TYPE(NCDIM), POINTER :: DIM_rivers

    LOGICAL:: FOUND

    real(SP), ALLOCATABLE :: RF(:),RT(:),RS(:)


    INTEGER :: I,  status

    
    DIM_time       => NC_MAKE_DIM(name='time',len=NF90_UNLIMITED)
    DIM_namelen => NC_MAKE_DIM(name='namelen',len=namelen)
    DIM_rivers => NC_MAKE_DIM(name='rivers',len=nrivers)
    DIM_DateStrLen => NC_MAKE_DIM(name='DateStrLen',len=DateStrLen)
    ALLOCATE(RT(nrivers),RS(Nrivers),RF(nrivers))

    NCF => NEW_FILE()

    
    ALLOCATE(NCF%FTIME)

    NCF%FNAME = TRIM(OUTPUT_DIR)//FOUT

    ATT => NC_MAKE_ATT(name='type',values="FVCOM RIVER FORCING FILE") 
    NCF => ADD(NCF,ATT)

    ATT => NC_MAKE_ATT(name='title',values=TRIM(TITLE)) 
    NCF => ADD(NCF,ATT)

    ATT => NC_MAKE_ATT(name='website',values=TRIM(WEBSITE)) 
    NCF => ADD(NCF,ATT)

    ATT => NC_MAKE_ATT(name='history',values=TRIM(HISTORY)) 
    NCF => ADD(NCF,ATT)

    
    VAR  => NC_MAKE_AVAR(name='river_names',values=river_names,DIM1=DIM_namelen,DIM2=DIM_rivers)
    NCF  => ADD(NCF,VAR)

    ! time
    VAR => FLOAT_TIME_OBJECT &
         &(USE_MJD=use_real_world_time, &
         & DIM=DIM_TIME)

    NCF  => ADD(NCF,VAR)


    ! Itime
    VAR  => ITIME_OBJECT &
         &(Use_MJD=use_real_world_time, &
         & DIM=DIM_TIME)

    NCF  => ADD(NCF,VAR)

    ! Itime2
    VAR => ITIME2_OBJECT &
         &(Use_MJD=use_real_world_time, &
         & DIM=DIM_TIME)

    NCF => ADD(NCF,VAR)

    IF (use_real_world_time) THEN

       VAR => DATETIME_OBJECT &
            &(DIMSTR=DIM_DateStrLen,&
            & DIMTIME=DIM_TIME)
       
       NCF  => ADD(NCF,VAR)
    END IF


    ! RIVER FLUX
    VAR  => NC_MAKE_AVAR(name='river_flux',values=RF, DIM1= DIM_rivers, DIM2= DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='river runoff volume flux') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='m^3s^-1') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! RIVER TEMP
    VAR  => NC_MAKE_AVAR(name='river_temp',values=RT, DIM1= DIM_rivers, DIM2= DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='river runoff temperature') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='Celsius') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)


    ! RIVER SALT
    VAR  => NC_MAKE_AVAR(name='river_salt',values=RS, DIM1= DIM_rivers, DIM2= DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='river runoff salinity') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='PSU') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)
    
    CALL NC_WRITE_FILE(NCF)

    write(ipt,*) "Created file, Dumping data:"

    DO I = 1, nTimes

       RF = qdis(:,I)
       RS = qsalt(:,I)
       RT = qtemp(:,I)

       ! UPDATE TIME STUFF
       VAR1 => FIND_VAR(NCF,"time",FOUND)
       IF(FOUND) CALL UPDATE_FLOAT_TIME(VAR1,TTIME(I))
       
       VAR1 => FIND_VAR(NCF,"Itime",FOUND)
       IF(FOUND) THEN
          VAR2 => FIND_VAR(NCF,"Itime2",FOUND)
          IF (.NOT.FOUND) THEN
             CALL WARNING&
                  & ("FOUND ONLY PART OF INTEGER TIME VARIABLE IN OUT PUT FILE!")
          ELSE
             CALL UPDATE_ITIME(VAR1,VAR2,TTIME(I))
          END IF
       END IF
       VAR1 => FIND_VAR(NCF,"Times",FOUND)
       IF(FOUND) CALL UPDATE_DATETIME(VAR1,TTIME(I))
       
       ! NOW WRITE A NEW FRAME

       NCF%FTIME%NEXT_STKCNT = NCF%FTIME%NEXT_STKCNT + 1
       CALL NC_WRITE_FILE(NCF)

    END DO



  END SUBROUTINE DUMP_OUTFILE

end module mod_rivers
