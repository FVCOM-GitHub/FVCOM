module mod_ts_clim
  use all_vars
  use mod_utils
  use mod_input
  use mod_prec
  use mod_ncll
  use mod_nctools
  USE MOD_TIME
  USE MOD_NCDIO

  USE SINTER

  implicit none

  INTEGER KSL                                 !!NUMBER OF STANDARD SEA LEVELS 
  REAL(SP), ALLOCATABLE :: DPTHSL(:)          !!DEPTH AT STANDARD SEA LEVEL
  REAL(SP), ALLOCATABLE :: TSL(:,:),SSL(:,:)  !!T/S AT STANDARD SEA LEVEL

  Character(Len=120):: FNAME
  Character(Len=120):: OLD_TS_FILE
  INTEGER, PARAMETER :: TSUNIT = 101
  CHARACTER(len=80) :: OUTPUT_FILE 
  CHARACTER(len=80) :: OUT_INTERVAL
  CHARACTER(len=80) :: ZERO_TIME
  
  LOGICAL :: FEXIST

  TYPE(TIME), SAVE :: INTERVAL

  TYPE(NCFILE), POINTER :: NC_OUT

  ! FOR FVCOM GRID OUTPUT
   REAL(SP), POINTER :: TC(:,:), SC(:,:)

  TYPE(NCFILE), POINTER :: NCF
  TYPE(GRID), SAVE :: MYGRID

  NAMELIST /NML_TS_CLIM/                 &
       & INPUT_DIR,                      &
       & OUTPUT_DIR,                     &
       & OLD_TS_FILE,                    &
       & OUTPUT_FILE,                    &
       & GRID_FILE,                      &
       & GRID_FILE_UNITS,                &
       & SIGMA_LEVELS_FILE,              &
       & DEPTH_FILE,                     &
       & PROJECTION_REFERENCE,           &
       & DATE_FORMAT,                    &
       & START_DATE,                     &
       & END_DATE,                       &
       & OUT_INTERVAL,                   &
       & TIMEZONE,                       &
       & ZERO_TIME

  ! DATA VARIABLES

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

          else if (opt_sng == "Fileame" .or.opt_sng == "filename"&
               & .or.opt_sng == "FILENAME") then

             call ftn_arg_get(arg_idx,arg_val,FName) ! [sng] Input file
             FName=FName(1:ftn_strlen(FName))
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

    WRITE(IPT,*) "Add better help here!"
    WRITE(IPT,*) "! OPTIONS:"
    WRITE(IPT,*) "! --filename=<a namelist file>"
    WRITE(IPT,*) "!   The namelist runfile for the program! "
    WRITE(IPT,*) "!   "
    WRITE(IPT,*) "!   Namelist OPTIONS: "
    WRITE(IPT,*) "!    "
    WRITE(IPT,*) "!   INPUT_DIR              # input path"
    WRITE(IPT,*) "!   OUTPUT_DIR             # output path"
    WRITE(IPT,*) "!   INPUT_FILE             # filename"
    WRITE(IPT,*) "!   OUT_FILE               # filename"
    WRITE(IPT,*) "!   GRID_FILE              # filename"
    WRITE(IPT,*) "!   GRID_FILE_UNITS        # 'meters' or 'degrees'"
    WRITE(IPT,*) "!   SIGMA_LEVELS_FILE      # filename"
    WRITE(IPT,*) "!   DEPTH_FILE             # filename"
    WRITE(IPT,*) "!   PROJECTION_REFERENCE   # see proj for help"
    WRITE(IPT,*) "!   DATE_FORMAT            # 'ymd' or 'dmy' ?"
    WRITE(IPT,*) "!   START_DATE             # Start Here"
    WRITE(IPT,*) "!   END_DATE               # End Here"
    WRITE(IPT,*) "!   OUT_INTERVAL           # 'seconds = n.n' or 'days=n.n'"
    WRITE(IPT,*) "!   TIMEZONE               # The Time zone or none"
    WRITE(IPT,*) "!   ZERO_TIME              # For binary output, when is zero?"

    WRITE(IPT,*) "!    "
    WRITE(IPT,*) "!    Exmaple Namelist"
    WRITE(IPT,*) "    "
    WRITE(IPT,*) "    "
    WRITE(UNIT=IPT,NML=NML_TS_CLIM)
    WRITE(IPT,*) "! NOTES: This program is parallel!"

  END SUBROUTINE MYHELPTXT

  SUBROUTINE READ_NAMELIST
    IMPLICIT NONE
    integer :: ios, i
    if(DBG_SET(dbg_sbr)) &
         & write(IPT,*) "Subroutine Begins: Read_Name_List;"

    if(DBG_SET(dbg_io)) &
         & write(IPT,*) "Read_Name_List: File: ",trim(FNAME)

    CALL FOPEN(NMLUNIT,trim(FNAME),'cfr')

    !READ NAME LIST FILE

    ! Read IO Information
    READ(UNIT=NMLUNIT, NML=NML_TS_CLIM,IOSTAT=ios)
    if(ios .NE. 0 ) then
            CALL FATAL_ERROR("Can Not Read NameList NML_INIT from file: "//trim(FNAME))
    end if
    REWIND(NMLUNIT)

    write(IPT,*) "Read_Name_List:"

    write(UNIT=IPT,NML=NML_TS_CLIM)

    CLOSE(NMLUNIT)


  END SUBROUTINE READ_NAMELIST


  SUBROUTINE OPEN_FILES
    IMPLICIT NONE

    TYPE(NCFILE), POINTER :: NCF
    integer :: ncfileind, datfileind,ios,charnum, i
    logical :: fexist,back,connected
    character(len=100) :: testchar
    character(len=160) :: pathnfile
    character(len=2) :: cios

    back = .true.

    !Check Grid File and open:
    ! TEST FILE NAME
    charnum = index (GRID_FILE,".dat")
    if (charnum /= len_trim(GRID_FILE)-3)&
         & CALL WARNING("GRID FILE does not end in .dat", &
         & trim(GRID_FILE))
    ! OPEN FILE
    pathnfile = trim(INPUT_DIR)//trim(GRID_FILE)
    Call FOPEN(GRIDUNIT,trim(pathnfile),'cfr')

    !Check Sigma File and open:
    ! TEST FILE NAME
    charnum = index (SIGMA_LEVELS_FILE,".dat")
    if (charnum /= len_trim(SIGMA_LEVELS_FILE)-3)&
         & CALL WARNING("SIGMA LEVELS FILE does not end in .dat", &
         & trim(SIGMA_LEVELS_FILE))
    ! OPEN FILE
    pathnfile = trim(INPUT_DIR)//trim(SIGMA_LEVELS_FILE)
    Call FOPEN(SIGMAUNIT,trim(pathnfile),'cfr')

    !Check Depth File and open:
    ! TEST FILE NAME
    charnum = index (DEPTH_FILE,".dat")
    if (charnum /= len_trim(DEPTH_FILE)-3)&
         & CALL WARNING("DEPTH FILE does not end in .dat", &
         & trim(DEPTH_FILE))
    ! OPEN FILE
    pathnfile = trim(INPUT_DIR)//trim(DEPTH_FILE)
    Call FOPEN(DEPTHUNIT,trim(pathnfile),'cfr')

    Pathnfile = trim(INPUT_DIR)//trim(OLD_TS_FILE)
    Call FOPEN(TSUNIT,trim(pathnfile),'cfr')

  END SUBROUTINE OPEN_FILES

  SUBROUTINE ALLOCATE_SPACE
    IMPLICIT NONE

    ALLOCATE(XM(0:MT),YM(0:MT),LON(0:MT),LAT(0:MT))
    ALLOCATE(XC(0:NT),XMC(0:NT),YC(0:NT),YMC(0:NT),LONC(0:NT),LATC(0:NT))


    allocate(H(0:mgl));  !H=0.0_sp
    allocate(el(0:mgl));  el=0.0_sp
    allocate(vx(0:mgl)); !vx=0.0_sp
    allocate(vy(0:mgl)); !vy=0.0_sp

    ALLOCATE(Z(0:MGL,KB)); z=0.0_sp
    ALLOCATE(Z1(0:NGL,KB)); z1=0.0_sp
    ALLOCATE(ZZ(0:MGL,KB)); ZZ=0.0_sp
    ALLOCATE(ZZ1(0:NGL,KB)); ZZ1=0.0_sp
    ALLOCATE(DZ(0:MGL,KB)); DZ=0.0_SP
    ALLOCATE(DZ1(0:NGL,KB)); DZ1=0.0_SP
    ALLOCATE(DZZ(0:MGL,KB)); DZZ=0.0_SP
    ALLOCATE(DZZ1(0:NGL,KB)); DZZ1=0.0_SP

    !ALLOCATE(TC(1:MGL,KBM1)); TC=0.0_SP
    !ALLOCATE(SC(1:mGL,KBM1)); SC=0.0_SP

!    ALLOCATE(TCGL(1:MGL,KBM1)); !TC=0.0_SP
!    ALLOCATE(SCGL(1:mGL,KBM1)); !SC=0.0_SP
!    TC => TCGL
!    SC => SCGL

    ALLOCATE(T1(0:MGL,1:KB))
    ALLOCATE(S1(0:mGL,1:KB))
    ALLOCATE(RHO1(0:mGL,1:KB))
    TC => T1
    SC => S1

    ALLOCATE(T(0:NGL,KB)) 
    ALLOCATE(S(0:NGL,KB))
    ALLOCATE(RHO(0:NGL,1:KB))

    ALLOCATE(GRAV_N(0:MGL)) ;GRAV_N=9.81_SP
    !GRAV_N(I) = 9.78032+0.005172*SIN(VY(I))**2-0.00006*SIN(2.0*VY(I))**2

    ALLOCATE(D(0:MGL))      

  END SUBROUTINE ALLOCATE_SPACE

!-------------------------------------------------------------------
  SUBROUTINE READ_TSC(now,TC,SC)
    IMPLICIT NONE
    REAL(SP), POINTER :: TC(:,:), SC(:,:)
    TYPE(TIME) :: NOW
    REAL(DP) :: hour,Second
    integer :: i,K, SOURCE, ios
    Real(SP),ALLOCATABLE,DIMENSION(:,:)  :: TCI,SCI
    REAL(SP),DIMENSION(KBM1)      :: TI,SI,ZI

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START: READ_TS"

      ALLOCATE(SCI(0:MGL,KB));SCI=0.0_SP
      ALLOCATE(TCI(0:MGL,KB));TCI=0.0_SP

!!------------------------------------------------------------------
!!   read ascii observation and intrepolate-------------------------
       READ(TSUNIT,*) hour
       
       !CALL IOERROR(IOS,"Can't read hour from TSfile")

       DO I=1,MGL 
          READ(TSUNIT,*) (TSL(I,K), K=1,KSL)
          READ(TSUNIT,*) (SSL(I,K), K=1,KSL)
       END DO

    ! NOW INTERPOLATE FROM STANDARD LEVELS TO THE VALUE AT EACH NODE
    DO I=1,MGL
       DO K=1,KBM1
          ZI(K)=ZZ(I,K)*H(I)  !+EL(I)
       END DO

       ! REAL T AND S
       CALL SINTER_EXTRP_UP(DPTHSL,TSL(I,:),ZI,TI,KSL,KBM1)
       CALL SINTER_EXTRP_UP(DPTHSL,SSL(I,:),ZI,SI,KSL,KBM1)

       TCI(I,1:KBM1) = TI(1:KBM1)
       SCI(I,1:KBM1) = SI(1:KBM1)
    END DO

    where(SCI<0.0_SP) SCI=0.0_SP
!!------------------------------------------------------------------
     TC(0:MGL,1:KB) =TCI(0:MGL,1:KB)
     SC(0:MGL,1:KB) =SCI(0:MGL,1:KB)

     TC(1:MGL,1) =TSL(1:MGL,1)
     SC(1:MGL,1) =SSL(1:MGL,1)

    where(SC<0.0_SP) SC=0.0_SP

     Now = seconds2time(hour*3600.) + ZEROTIME
     print*,Now

    DEALLOCATE(TCI)
    DEALLOCATE(SCI)

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END: READ_TS"

  END SUBROUTINE READ_TSC

  SUBROUTINE READ_OLD_TS_INF
    IMPLICIT NONE
    CHARACTER(LEN=80)    :: scan_result
    INTEGER :: ISCAN, I, K, IOS
    REAL(SP), DIMENSION(150) :: TEMP

    IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) "START READING INIT TS FILE"

    ISCAN = SCAN_FILE(TSUNIT,"Standard Levels",FVEC =TEMP  ,NSZE = KSL)
    IF(ISCAN /= 0) then
       write(scan_result,'(I2)') ISCAN
            call fatal_error('Improper formatting of ITS file: ISCAN ERROR# '//trim(scan_result),&
            & 'The header must contain: "Standard Levels ="', &
            & 'Followed by a series of floating point depths.',&
            & 'Upward is positive!')
    END IF


    write(ipt,*) "! # of standard levels:=",KSL

    ALLOCATE(DPTHSL(KSL)); DPTHSL = TEMP(1:KSL)
    ALLOCATE(TSL(MGL,KSL))
    ALLOCATE(SSL(MGL,KSL))

   ! FIND FIRST LINE of )BC ARRAY
    rewind TSUNIT
    read(TSUNIT,*)  !! Skip the first line

  END SUBROUTINE READ_OLD_TS_INF
 
  SUBROUTINE GET_FVCOM_GRID
    USE MOD_SETUP
    IMPLICIT NONE
    CHARACTER(LEN=80) FNAME
    INTEGER STATUS

    ! OPEN AND READ THE FVCOM GRID FILE
    IF (MSR) THEN
       FNAME = TRIM(INPUT_DIR)//TRIM(GRID_FILE)
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
    ALLOCATE(XC(0:NT),XMC(0:NT),YC(0:NT),YMC(0:NT),LONC(0:NT),LATC(0:NT))

    CALL LOAD_COLDSTART_COORDS(X_GBL,Y_GBL,X_LCL,Y_LCL)
    CALL COORDINATE_UNITS(X_LCL,Y_LCL)
    CALL SETUP_CENTER_COORDS

    DEALLOCATE(X_LCL)
    DEALLOCATE(Y_LCL)

    IF(MSR) THEN
       DEALLOCATE(X_GBL)
       DEALLOCATE(Y_GBL)
    END IF

  CALL READ_COLDSTART_SIGMA
  CLOSE(SIGMAUNIT)
  KBM1 = KB - 1
  KBM2 = KB - 2

!  CALL ALLOCATE_SPACE

!  CALL READ_COLDSTART_COORDS(GRIDUNIT,MGL,VX,VY)
!  CLOSE(GRIDUNIT)

  CALL READ_COLDSTART_DEPTH(DEPTHUNIT,MGL,VX,VY,H)
  CLOSE(DEPTHUNIT)

  CALL Setup_Sigma
  CALL SETUP_SIGMA_DERIVATIVES

  END SUBROUTINE GET_FVCOM_GRID

  SUBROUTINE RUN_TGE
    USE MOD_OBCS, only : iobcn
    implicit none
    INTEGER :: NCT

    NCT = NT*3
    IOBCN = 0

    ALLOCATE(NBE(0:NT,3))         ;NBE      = 0 !!INDICES OF ELMNT NEIGHBORS
    ALLOCATE(NTVE(0:MT))          ;NTVE     = 0 
    ALLOCATE(NTSN(MT))            ;NTSN     = 0 
    ALLOCATE(ISONB(0:MT))         ;ISONB    = 0  !!NODE MARKER = 0,1,2
    ALLOCATE(ISBCE(0:NT))         ;ISBCE    = 0 
    ALLOCATE(NIEC(NCT,2))         ;NIEC     = 0
    ALLOCATE(NTRG(NCT))           ;NTRG     = 0    
    ! POSITION OF NODAL CONTROL VOLUME CORNERS 
    ALLOCATE(XIJE(NCT,2))         ;XIJE     = ZERO
    ALLOCATE(YIJE(NCT,2))         ;YIJE     = ZERO 

    ! LENGTH OF NODAL CONTROL VOLUME EDGES
    ALLOCATE(DLTXE(NCT))          ;DLTXE    = ZERO
    ALLOCATE(DLTYE(NCT))          ;DLTYE    = ZERO
    ALLOCATE(DLTXYE(NCT))         ;DLTXYE   = ZERO !! TOTAL LENGTH
    ALLOCATE(SITAE(NCT))          ;SITAE    = ZERO !! ANGLE

    CALL TRIANGLE_GRID_EDGE


  END SUBROUTINE RUN_TGE



  SUBROUTINE SET_TIME
    USE MOD_SET_TIME
    IMPLICIT NONE
    integer status
    INTEGER(ITIME) :: dummy
    CHARACTER(LEN=4) :: FLAG

    IF (USE_REAL_WORLD_TIME) THEN
       StartTime = read_datetime(START_DATE,DATE_FORMAT,TIMEZONE,status)
       if (status /= 1) & 
            & Call Fatal_Error("Could not read the date string START_DATE: ", trim(START_DATE))

       ! GET THE END TIME
       EndTime = READ_DateTime(END_DATE,DATE_FORMAT,TIMEZONE,status)
       if (status /= 1) &
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

    CALL PRINT_REAL_TIME(STARTTIME,IPT,"READ: START_DATE",TIMEZONE)
    CALL PRINT_REAL_TIME(ENDTIME,IPT,"READ: END_DATE",TIMEZONE)
    CALL PRINT_REAL_TIME(ZEROTIME,IPT,"READ: ZERO_TIME",TIMEZONE)


    ! SANITY ECK
    if(StartTime .GT. EndTime) &
         & Call Fatal_Error("Runfile Start_Date exceeds or equal to End_Date")

    CALL IDEAL_TIME_STRING2TIME(OUT_INTERVAL,FLAG,INTERVAL,dummy)
    ! SANITY CHECK
    IF (FLAG /= 'time') CALL FATAL_ERROR&
         & ("TIME INTERVAL SPECIFICATION IS INCORRENT",&
         & "INTERVAL MUST BE 'seconds=' or 'days='")

  END SUBROUTINE SET_TIME

!!-----------------------------------------------------------------------------
!!  Setup the netcdf file for output

  SUBROUTINE MAKE_TS_OUTFILE(FNAME)
    USE MOD_NCDIO
    !USE MOD_FORCE, only : fvcom_cap_grid_source
    IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN) :: FNAME
    TYPE(NCFILE), POINTER :: NCF
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCATT),  POINTER :: ATT
    LOGICAL :: FOUND

    Call SET_FVCOM_GRID(MYGRID)    

    CALL DEFINE_DIMENSIONS(MYGRID)

    ! ALLOCATE THE NEW FILE OBJECT
    NCF => NEW_FILE()

    NC_OUT => NCF

    NCF%FTIME => new_ftime()

    NCF%FNAME = TRIM(FNAME)

    NCF => ADD(NCF,GRID_FILE_OBJECT(MYGRID) )
    
    NCF => ADD(NCF,TIME_FILE_OBJECT() )

    NCF => ADD(NCF,ZETA_FILE_OBJECT() )

    ATT => FIND_ATT(NCF,'source',FOUND)
    IF(.NOT.FOUND) CALL FATAL_ERROR("LOOKING FOR 'source' ATTRIBUTE: NOT FOUND")
    
!    ATT%CHR = fvcom_cap_grid_source

!    NCF => ADD(NCF,FILE_DATE_OBJECT() )

    ! TC
       
    VAR  => NC_MAKE_AVAR(name='temp_clim',&
    !VAR  => NC_MAKE_AVAR(name='temp',&
         & values= T1, DIM1= DIM_node, DIM2= DIM_siglay, DIM3 = DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='temperature')
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='standard_name',values='sea_water_temperature')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='degrees_C')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar)
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data')

    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)
    ! SC
    VAR  => NC_MAKE_AVAR(name='salinity_clim',&
    !VAR  => NC_MAKE_AVAR(name='salinity',&
         & values= S1, DIM1= DIM_node, DIM2= DIM_siglay, DIM3 = DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='salinity')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='standard_name',values='sea_water_salinity')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='1e-3')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar)
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data')
    VAR  => ADD(VAR,ATT)
    NCF  => ADD(NCF,VAR)

    NCF%FTIME%NEXT_STKCNT = 0
    CALL NC_WRITE_FILE(NC_OUT)


  END SUBROUTINE MAKE_TS_OUTFILE


end module mod_ts_clim

