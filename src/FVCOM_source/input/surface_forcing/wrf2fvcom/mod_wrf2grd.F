MODULE MOD_WRF2GRD
  USE MOD_NCTOOLS
  USE MOD_UTILS
  USE MOD_INPUT
  USE MOD_TIME
  USE CONTROL
  USE LIMS
  USE ALL_VARS
  IMPLICIT NONE


  CHARACTER(len=80) :: RUN_FILE

  CHARACTER(len=80) :: wrf_in_file
  CHARACTER(len=80) :: wrf_out_file
  CHARACTER(len=80) :: OUT_FILE_TYPE
  CHARACTER(len=80) :: OUT_INTERVAL
  CHARACTER(len=80) :: ZERO_TIME

  LOGICAL           :: FEXIST

  TYPE(TIME)  :: INTERVAL

  TYPE(NCFILE), POINTER :: NC_OUT


  NAMELIST /NML_WRF2FVCOM/               &
       & WRF_IN_FILE,                    &
       & WRF_OUT_FILE,                   &
       & OUT_FILE_TYPE,                  &
       & GRID_FILE,                      &
       & GRID_FILE_UNITS,                &
       & PROJECTION_REFERENCE,           &
       & DATE_FORMAT,                    &
       & START_DATE,                     &
       & END_DATE,                       &
       & OUT_INTERVAL,                   &
       & ZERO_TIME

CONTAINS


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
       if(MSR) WRITE(IPT,*) "You must specify: '--filename=<namelist>' "
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

          else if (opt_sng == "dbg_par" .or.opt_sng == "Dbg_Par"&
               & .or.opt_sng == "DBG_PAR") then

             dbg_par = .true.

          else if (opt_sng == "FileName" .or.opt_sng == "FILENAME"&
               & .or.opt_sng == "filename") then

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

             if(MSR) call MYHelpTxt


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

          if(MSR) call MYHelpTxt
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

    WRITE(IPT,*) "Add better help here!"
    WRITE(IPT,*) "! OPTIONS:"
    WRITE(IPT,*) "! --filename=<a namelist file>"
    WRITE(IPT,*) "!   The namelist runfile for the program! "
    WRITE(IPT,*) "!   "
    WRITE(IPT,*) "!   Namelist OPTIONS: "
    WRITE(IPT,*) "!    "
    WRITE(IPT,*) "!   WRF_IN_FILE            # path and file"
    WRITE(IPT,*) "!   WRF_OUT_FILE           # path and file"
    WRITE(IPT,*) "!   OUT_FILE_TYPE          # 'binary','netcdf' or 'both'" 
    WRITE(IPT,*) "!   GRID_FILE              # path and file"
    WRITE(IPT,*) "!   GRID_FILE_UNITS        # 'meters' or 'degrees'"
    WRITE(IPT,*) "!   PROJECTION_REFERENCE   # see proj for help" 
    WRITE(IPT,*) "!   DATE_FORMAT            # 'ymd' or 'dmy' ?"
    WRITE(IPT,*) "!   START_DATE             # Start Here"
    WRITE(IPT,*) "!   END_DATE               # End Here"
    WRITE(IPT,*) "!   OUT_INTERVAL           # 'seconds = n.n' or 'days=n.n'"
    WRITE(IPT,*) "!   ZERO_TIME              # For binary output, when is zero?"

    WRITE(IPT,*) "!    "
    WRITE(IPT,*) "! NOTES: This program is parallel!"
  END SUBROUTINE MYHELPTXT


  SUBROUTINE GET_RUN_FILE
    IMPLICIT NONE
    integer :: ios, i
    Character(Len=120):: FNAME
    if(DBG_SET(dbg_sbr)) &
         & write(IPT,*) "Subroutine Begins: GET_RUNFILE;"


    FNAME = RUN_FILE
    
    if(DBG_SET(dbg_io)) &
         & write(IPT,*) "Read_Name_List: File: ",trim(FNAME)
    
    CALL FOPEN(NMLUNIT,trim(FNAME),'cfr')
    
    !READ NAME LIST FILE
    
    ! Read IO Information
    READ(UNIT=NMLUNIT, NML=NML_WRF2FVCOM,IOSTAT=ios)
    if(ios .NE. 0 ) then
       if(DBG_SET(dbg_log)) &
            & write(IPT,*)"Can Not Read NameList NML_IO from file: "//trim(FNAME)
    end if
    REWIND(NMLUNIT)

    if(DBG_SET(dbg_scl)) &
         & write(IPT,*) "Read_Name_List:"

    if(DBG_SET(dbg_scl)) &
         & write(UNIT=IPT,NML=NML_WRF2FVCOM)


  END SUBROUTINE GET_RUN_FILE




  SUBROUTINE GET_FVCOM_GRID
    USE MOD_SETUP
    IMPLICIT NONE
    CHARACTER(LEN=80) FNAME
    INTEGER STATUS

    ! OPEN AND READ THE FVCOM GRID FILE
    IF (MSR) THEN
       FNAME = GRID_FILE
       WRITE(IPT,*) "OPENING GRIDFILE: "//TRIM(FNAME)
       Call FOPEN(GRIDUNIT,TRIM(FNAME),'cfr')
    END IF
    
    CALL LOAD_COLDSTART_GRID(NVG)
    KB = 1

    CALL SETUP_DOMAIN

    IF(MSR) THEN
       ! ALLOCATE SPACE FOR THE GLOBAL GRID DATA
       ALLOCATE(Y_GBL(0:MGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE Y_GBL")
       ALLOCATE(X_GBL(0:MGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE X_GBL")
    END IF

    ! ALLOCATE SPACE FOR THE LOCAL GRID DATA
    ALLOCATE(Y_LCL(0:MT),stat=status)
    IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE Y_LCL")
    ALLOCATE(X_LCL(0:MT),stat=status)
    IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE X_LCL")


    ALLOCATE(VX(0:MT),VY(0:MT),XM(0:MT),YM(0:MT),LON(0:MT),LAT(0:MT))
    ALLOCATE(XC(0:NT),XMC(0:NT),YC(0:NT),YMC(0:NT))

    CALL LOAD_COLDSTART_COORDS(X_GBL,Y_GBL,X_LCL,Y_LCL)
    CALL COORDINATE_UNITS
    CALL SETUP_COORDS

    DEALLOCATE(X_LCL)
    DEALLOCATE(Y_LCL)

    IF(MSR) THEN
        DEALLOCATE(X_GBL)
        DEALLOCATE(Y_GBL)
     END IF


  END SUBROUTINE GET_FVCOM_GRID


  SUBROUTINE SET_TIME
    USE MOD_SET_TIME
    IMPLICIT NONE
    integer status
    INTEGER(ITIME) :: dummy
    CHARACTER(LEN=4) :: FLAG
    
    CALL READ_TIME_STRING(START_DATE,DATE_FORMAT,StartTime,status)
    if (.not. status) &
         & Call Fatal_Error("Could not read the date string START_DATE: ",&
         & trim(START_DATE))
    
    ! GET THE END TIME
    CALL READ_TIME_STRING(END_DATE,DATE_FORMAT,EndTime,status)
    if (.not. status) &
         & Call Fatal_Error("Could not read the date string END_DATE:",&
         & trim(END_DATE))
    
    ! SANITY CHECK
    if(StartTime .GT. EndTime) &
         & Call Fatal_Error("Runfile Start_Date exceeds or equal to End_Date")
    
    CALL IDEAL_TIME_STRING2TIME(OUT_INTERVAL,FLAG,INTERVAL,dummy)
    
    ! SANITY CHECK
    IF (FLAG /= 'time') CALL FATAL_ERROR&
         & ("TIME INTERVAL SPECIFICATION IS INCORRENT",&
         & "INTERVAL MUST BE 'seconds=' or 'days='")
    
    

  END SUBROUTINE SET_TIME



  SUBROUTINE LOAD_WRF_IN
    USE MOD_FORCE
    IMPLICIT NONE
    TYPE(NCFILE), POINTER :: NCF

    HEATING_ON = .TRUE.
    WIND_ON    = .TRUE.
    PRECIPITATION_ON = .TRUE.

    HEATING_FILE = WRF_IN_FILE
    WIND_FILE    = WRF_IN_FILE
    PRECIPITATION_FILE = WRF_IN_FILE
    

    CALL  NC_INIT(NCF,WRF_IN_FILE)
    
    ! OPEN THE FILE AND LOAD METADATA
    if(.not. NCF%OPEN) then
       Call NC_OPEN(NCF)
       CALL NC_LOAD(NCF)
       FILEHEAD => ADD(FILEHEAD,NCF)
    end if
    
    CALL SETUP_FORCING
    
  END SUBROUTINE LOAD_WRF_IN

  
SUBROUTINE dump_binary_force(NOW)
  USE ALL_VARS
  USE MOD_TIME
  USE MOD_UTILS
  USE MOD_PAR
  IMPLICIT NONE
  TYPE(TIME), INTENT(IN) :: NOW


  LOGICAL, SAVE :: INIT = .FALSE.
  INTEGER, PARAMETER :: wndunit=118 
  INTEGER, PARAMETER :: hfxunit=119
  INTEGER, PARAMETER :: evpunit=120
  REAL(SP) :: hour, mean
  TYPE(TIME), SAVE :: ZEROTIME, MONTH
  INTEGER :: STATUS, I
  REAL(SP), ALLOCATABLE :: WSUGL(:), WSVGL(:), SWRGL(:),NHFGL(:), evapgl(:),precgl(:)


  IF(MSR) THEN
     IF (.NOT. ALLOCATED(WSUGL)) ALLOCATE(WSUGL(0:NGL))
     IF (.NOT. ALLOCATED(WSVGL)) ALLOCATE(WSVGL(0:NGL))

     IF (.NOT. ALLOCATED(evapGL)) ALLOCATE(evapGL(0:MGL))
     IF (.NOT. ALLOCATED(precGL)) ALLOCATE(precGL(0:MGL))

     IF (.NOT. ALLOCATED(SWRGL)) ALLOCATE(SWRGL(0:MGL))
     IF (.NOT. ALLOCATED(NHFGL)) ALLOCATE(NHFGL(0:MGL))
  END IF

!  write(ipt,*)"Wusurf: min/max",minval(wusurf),maxval(wusurf)

  IF (PAR)THEN
# if defined(MULTIPROCESSOR)
     CALL COLLECTA(MYID,MSRID,NPROCS,EMAP,WUSURF,WSUGL)
     CALL COLLECTA(MYID,MSRID,NPROCS,EMAP,WVSURF,WSVGL)
      
      CALL COLLECTA(MYID,MSRID,NPROCS,NMAP,qprec2,precgl)
      CALL COLLECTA(MYID,MSRID,NPROCS,NMAP,qevap2,evapgl)

     CALL COLLECTA(MYID,MSRID,NPROCS,NMAP,SWRAD_WATTS,SWRGL)
     CALL COLLECTA(MYID,MSRID,NPROCS,NMAP,WTSURF_WATTS,NHFGL)
# endif
  ELSE 
     
     WSUGL(1:NT) = WUSURF(1:NT)
     WSVGL(1:NT) = WVSURF(1:NT)
     
      precgl(1:MT) = Qprec2(1:MT)
      evapgl(1:MT) = Qevap2(1:MT)

     SWRGL(1:MT) = SWRAD_WATTS(1:MT)
     NHFGL(1:MT) = WTSURF_WATTS(1:MT)
     
  END IF


  
  IF (.NOT. MSR) RETURN


  IF(.not. INIT) THEN
     
     INIT = .TRUE.
     open(wndunit,file=TRIM(WRF_OUT_FILE)//'_wnd.dat',form='unformatted',status='unknown')
     open(hfxunit,file=TRIM(WRF_OUT_FILE)//'_hfx.dat',form='unformatted',status='unknown')
     open(evpunit,file=TRIM(WRF_OUT_FILE)//'_evp.dat',form='unformatted',status='unknown')

     CALL READ_TIME_STRING(ZERO_TIME,'ymd',ZEROTIME,status)
     IF (.not. STATUS) CALL FATAL_ERROR("COULD NOT READ ZEROTIME STRING")
     

  END IF

  hour = REAL_TIME_DIFF(NOW,ZEROTIME) / 3600.0_SP

  write(ipt,*) "time", hour

  mean = sum(wsugl(1:NGL))/real(NGL,SP)
  if(isnan(mean)) call fatal_error("wsugl is nan")
  write(ipt,*) "WSUGL: min/max/mean",minval(wsugl(1:NGL)),maxval(wsugl(1:NGL)),mean


  mean = sum(wsvgl(1:NGL))/real(NGL,SP)
  if(isnan(mean)) call fatal_error("wsvgl is nan")
  write(ipt,*) "WSVGL: min/max/mean",minval(wsvgl(1:NGL)),maxval(wsvgl(1:NGL)),mean

  mean = sum(swrgl(1:mGL))/real(NGL,SP)
  if(isnan(mean)) call fatal_error("swrgl is nan")
  write(ipt,*) "SWRGL: min/max/mean",minval(swrgl(1:MGL)),maxval(swrgl(1:MGL)),mean

  mean = sum(nhfgl(1:mGL))/real(NGL,SP)
  if(isnan(mean)) call fatal_error("nhfgl is nan")
  write(ipt,*) "NHFGL: min/max/mean",minval(nhfgl(1:MGL)),maxval(nhfgl(1:MGL)),mean

  mean = sum(evapgl(1:mGL))/real(MGL,SP)
  if(isnan(mean)) call fatal_error("EVAPgl is nan")
  write(ipt,*) "EVAPGL: min/max/mean",minval(EVAPgl(1:MGL)),maxval(EVAPgl(1:MGL)),mean

  mean = sum(PRECgl(1:mGL))/real(MGL,SP)
  if(isnan(mean)) call fatal_error("PRECgl is nan")
  write(ipt,*) "PRECGL: min/max/mean",minval(PRECgl(1:MGL)),maxval(PRECgl(1:MGL)),mean


  WRITE(wndunit) hour
  write(wndunit)(WSUGL(i),WSVGL(i),i=1,ngl)


  WRITE(hfxunit) hour
  write(hfxunit)(NHFGL(i),SWRGL(i),i=1,mgl)

  WRITE(evpunit) hour
  write(evpunit)(evapgl(i),precgl(i),i=1,mgl)


END SUBROUTINE dump_binary_force



  SUBROUTINE MAKE_NC_OUTFILE
    USE MOD_NCDIO
    IMPLICIT NONE
    TYPE(NCFILE), POINTER :: NCF
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCATT),  POINTER :: ATT
    TYPE(NCDIM),  POINTER :: DIM1
    TYPE(NCDIM),  POINTER :: DIM2
    TYPE(NCDIM),  POINTER :: DIM3
    TYPE(NCFTIME), POINTER ::FTM
    


  ! ALLOCATE THE NEW FILE OBJECT
  NCF => NEW_FILE()

  NC_OUT => NCF
  
  ALLOCATE(NCF%FTIME)


  NCF%FNAME = TRIM(WRF_OUT_FILE)//'.nc'

  NCF => ADD(NCF,GRID_FILE_OBJECT() )

  NCF => ADD(NCF,SURFACE_HEATING_FILE_OBJECT() )
  NCF => ADD(NCF,WIND_STRESS_FILE_OBJECT() )
  NCF => ADD(NCF,PRECIPITATION_FILE_OBJECT() )

  FTM => NCF%FTIME
  FTM%NEXT_STKCNT = 0
  CALL NC_WRITE_FILE(NC_OUT)
  

END SUBROUTINE MAKE_NC_OUTFILE



END MODULE MOD_WRF2GRD
