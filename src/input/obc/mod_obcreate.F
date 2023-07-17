module mod_obcreate
  use all_vars
  use mod_utils
  use mod_par
  use mod_input
  use mod_prec
  use mod_ncll
  use mod_nctools
  USE MOD_NCDIO
  USE MOD_TIME
  use bcs
  implicit none

! TIDAL PERIOD IN SECONDS
! REF::H.L. Simmons et al, Deep-Sea Research II 51(2004) 3043-3058;
! HISTORICALLY FVCOM USES TIDAL PERIODS CONSISTANT WITH THE CANONICAL
! FOUR DIGIT HOUR REPRESENTATION OF PERIOD
  REAL(SP), PARAMETER :: S2=43200.0_SP ! 12.00 hours
  REAL(SP), PARAMETER :: M2=44712.0_SP ! 12.42 hours
! REAL(SP), PARAMETER :: M2=44714.165191868036_SP
  REAL(SP), PARAMETER :: N2=45570.0_SP ! 12.66 hours
! REAL(SP), PARAMETER :: N2=45570.053511717721_SP
  REAL(SP), PARAMETER :: K2=43082.0_SP ! 11.97 hours
! REAL(SP), PARAMETER :: K2=43082.050318594716_SP
  REAL(SP), PARAMETER :: K1=86164.0_SP ! 23.93 hours
! REAL(SP), PARAMETER :: K1=86164.077005067069_SP
  REAL(SP), PARAMETER :: P1=86637.0_SP ! 24.07 hours
! REAL(SP), PARAMETER :: P1=86637.199771652769_SP
  REAL(SP), PARAMETER :: O1=92950.0_SP ! 25.82 hours
! REAL(SP), PARAMETER :: O1=92949.635700536543_SP
  REAL(SP), PARAMETER :: Q1=96726.0_SP ! 26.87 hours
! REAL(SP), PARAMETER :: Q1=96726.085702966622_SP


  REAL(SP), PARAMETER :: S2_EQI_AMP=0.112743_sp
  REAL(SP), PARAMETER :: M2_EQI_AMP=0.242334_sp
  REAL(SP), PARAMETER :: N2_EQI_AMP=0.046397_sp
  REAL(SP), PARAMETER :: K2_EQI_AMP=0.030684_sp
  REAL(SP), PARAMETER :: K1_EQI_AMP=0.141565_sp
  REAL(SP), PARAMETER :: P1_EQI_AMP=0.046848_sp
  REAL(SP), PARAMETER :: O1_EQI_AMP=0.100661_sp
  REAL(SP), PARAMETER :: Q1_EQI_AMP=0.019273_sp

  REAL(SP), PARAMETER :: S2_EQI_BETA=0.693_SP
  REAL(SP), PARAMETER :: M2_EQI_BETA=0.693_SP
  REAL(SP), PARAMETER :: N2_EQI_BETA=0.693_SP
  REAL(SP), PARAMETER :: K2_EQI_BETA=0.693_SP
  REAL(SP), PARAMETER :: K1_EQI_BETA=0.736_SP
  REAL(SP), PARAMETER :: P1_EQI_BETA=0.706_SP
  REAL(SP), PARAMETER :: O1_EQI_BETA=0.695_SP
  REAL(SP), PARAMETER :: Q1_EQI_BETA=0.695_SP


! PLEASE ADD MORE CONSTITUENTS HERE AS NEEDED
! NOTE YOU MUST ALSO ADD THEM IN THE CODING BELLOW


  Character(Len=120):: FNAME

  CHARACTER(LEN=80) :: ELEVATION_SOURCE_TYPE
  CHARACTER(LEN=80) :: ELEVATION_SOURCE_FILE
  CHARACTER(LEN=80) :: TS_SOURCE_TYPE
  CHARACTER(LEN=80) :: TS_SOURCE_FILE

  CHARACTER(LEN=80) :: time_origin


  NAMELIST /NML_OBC/                      &
       & INPUT_DIR,                       &
       & OUTPUT_DIR,                      &
       & GRID_FILE,                       &
       & SIGMA_LEVELS_FILE,               &
       & DEPTH_FILE,                      &
       & OBC_NODE_LIST_FILE,              &
       & ELEVATION_SOURCE_TYPE,           &
       & ELEVATION_SOURCE_FILE,           &
       & TIME_ORIGIN,                     &
       & TIMEZONE,                        &
       & TS_SOURCE_TYPE,                  &
       & TS_SOURCE_FILE


  INTEGER, PARAMETER :: ELVUNIT = 101
  INTEGER, PARAMETER :: TSUNIT = 102


  TYPE(NCDIM), POINTER :: DIM_tidal_components
  ! DATA VARIABLES

  TYPE(TIME), SAVE :: NOW
  TYPE(TIME), SAVE :: START
  TYPE(TIME), SAVE :: STEP


  ! SPECTRAL ELEVATION DATA
  INTEGER :: NCOMPS
  REAL(SP), ALLOCATABLE :: Eref(:)
  REAL(SP), ALLOCATABLE :: Eperiod(:)
  REAL(SP), ALLOCATABLE :: Eamp(:,:)
  REAL(SP), ALLOCATABLE :: Ephase(:,:)
  REAL(SP), ALLOCATABLE :: EQI_AMP(:)
  REAL(SP), ALLOCATABLE :: EQI_BETA(:)

  ! USE FROM: BCS
  !  CHARACTER(LEN=80), ALLOCATABLE :: TIDE_TYPE

  CHARACTER(LEN=160) :: COMMENTS, COMPONENTS

  ! TIME SERIES ELEVATION DATA
  REAL(SP), ALLOCATABLE :: ELEVATION(:,:)
  INTEGER :: NTIMES

  ! TS TIME SERIES NUDGIND DATA
  REAL(SP), ALLOCATABLE :: OBC_ZZ(:,:)
  REAL(SP), ALLOCATABLE :: OBC_Z(:,:)
  REAL(SP), ALLOCATABLE :: OBC_H(:)
  REAL(SP), ALLOCATABLE :: OBC_X(:)
  REAL(SP), ALLOCATABLE :: OBC_Y(:)
  REAL(SP), ALLOCATABLE :: OBC_DEPTH(:,:)
  REAL(SP), ALLOCATABLE :: OBC_TEMP(:,:,:)
  REAL(SP), ALLOCATABLE :: OBC_SALT(:,:,:)
  REAL(SP), ALLOCATABLE :: TIMES(:)
  TYPE(TIME), ALLOCATABLE :: MJDS(:)


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
    WRITE(IPT,*) "! --filename=XXX"
    WRITE(IPT,*) "!   The namelist runfile for the program! "
    WRITE(IPT,*) "!   "
    WRITE(IPT,*) "!   Namelist OPTIONS: "
    WRITE(IPT,*) "!    INPUT_DIR, (REQUIRED)"
    WRITE(IPT,*) "!    OUTPUT_DIR, (REQUIRED)"
    WRITE(IPT,*) "!    GRID_FILE, (REQUIRED)"
    WRITE(IPT,*) "!    SIGMA_LEVELS_FILE, (REQUIRED)"
    WRITE(IPT,*) "!    DEPTH_FILE, (REQUIRED)"
    WRITE(IPT,*) "!    OBC_NODE_LIST_FILE, (REQUIRED)"
    WRITE(IPT,*) "!    ELEVATION_SOURCE_FILE (OPTIONAL)"
    WRITE(IPT,*) "!    ELEVATION_SOURCE_TYPE: (OPTIONAL)"
    WRITE(IPT,*) "!    => convert2new_spectral"
    WRITE(IPT,*) "!    => convert2new_julian"
    WRITE(IPT,*) "!    => USER_DEFINED_JULIAN"
    WRITE(IPT,*) "!    => USER_DEFINED_SPECTRAL"
    WRITE(IPT,*) "!    TIME_ORIGIN,(SUGGESTED)"
    WRITE(IPT,*) "!    TIMEZONE,(SUGGESTED)"
    WRITE(IPT,*) "!    TS_SOURCE_FILE: (OPTIONAL)"
    WRITE(IPT,*) "!    TS_SOURCE_TYPE: (OPTIONAL)"
    WRITE(IPT,*) "!    => convert2new"
    WRITE(IPT,*) "!    => USER_DEFINED"
    WRITE(IPT,*) "!    "
    WRITE(IPT,*) "!    EXAMPLE NAMELIST:"
    write(UNIT=IPT,NML=NML_OBC)


    WRITE(IPT,*) "! NOTES: Do not run this program in parallel!"


  END SUBROUTINE MYHELPTXT

  SUBROUTINE READ_NAMELIST
    IMPLICIT NONE
    integer :: ios, i
    if(DBG_SET(dbg_sbr)) &
         & write(IPT,*) "Subroutine Begins: Read_Name_List;"


    if(DBG_SET(dbg_io)) &
         & write(IPT,*) "Read_Name_List: File: ",trim(FNAME)


    ! SET SOME DEFAULTS....
    time_origin="NONE"
    timezone = "NONE"

    CALL FOPEN(NMLUNIT,trim(FNAME),'cfr')

    !READ NAME LIST FILE

    ! Read IO Information
    READ(UNIT=NMLUNIT, NML=NML_OBC,IOSTAT=ios)
    if(ios .NE. 0 )THEN
       write(UNIT=IPT,NML=NML_OBC)
       
       CALL FATAL_ERROR &
         &("Can Not Read NameList NML_OBC from file: "//trim(FNAME))
    END if
    REWIND(NMLUNIT)

    if(DBG_SET(dbg_scl)) &
         & write(IPT,*) "Read_Name_List:"

    if(DBG_SET(dbg_scl)) &
         & write(UNIT=IPT,NML=NML_OBC)

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
    
    charnum = index (OBC_NODE_LIST_FILE,".dat")
    if (charnum /= len_trim(OBC_NODE_LIST_FILE)-3)&
         & CALL WARNING("OBC NODE LIST FILE does not end in .dat", &
         & trim(OBC_NODE_LIST_FILE))
    ! OPEN FILE
    pathnfile = trim(INPUT_DIR)//trim(OBC_NODE_LIST_FILE)
    Call FOPEN(OBCUNIT,trim(pathnfile),'cfr')
    
        
    !Check Sigma File and open:
    ! TEST FILE NAME
    charnum = index (SIGMA_LEVELS_FILE,".dat")
    if (charnum /= len_trim(SIGMA_LEVELS_FILE)-3)&
         & CALL WARNING("SIGMA LEVELS FILE does not end in .dat", &
         & trim(SIGMA_LEVELS_FILE))
    ! OPEN FILE
    pathnfile = trim(INPUT_DIR)//trim(SIGMA_LEVELS_FILE)
    Call FOPEN(SIGMAUNIT,trim(pathnfile),'cfr')
    
    !Check Grid File and open:
    ! TEST FILE NAME
    charnum = index (GRID_FILE,".dat")
    if (charnum /= len_trim(GRID_FILE)-3)&
         & CALL WARNING("GRID FILE does not end in .dat", &
         & trim(GRID_FILE))
    ! OPEN FILE
    pathnfile = trim(INPUT_DIR)//trim(GRID_FILE)
    Call FOPEN(GRIDUNIT,trim(pathnfile),'cfr')
    
    
    !Check Depth File and open:
    ! TEST FILE NAME
    charnum = index (DEPTH_FILE,".dat")
    if (charnum /= len_trim(DEPTH_FILE)-3)&
         & CALL WARNING("DEPTH FILE does not end in .dat", &
         & trim(DEPTH_FILE))
    ! OPEN FILE
    pathnfile = trim(INPUT_DIR)//trim(DEPTH_FILE)
    Call FOPEN(DEPTHUNIT,trim(pathnfile),'cfr')

  END SUBROUTINE OPEN_FILES

  SUBROUTINE OPEN_ELV_SOURCE
    IMPLICIT NONE
    integer :: charnum
    logical :: back
    character(len=160) :: pathnfile

    back = .true.
    !Check ELEVATION SOURCE DATA FILE and open:
    ! TEST FILE NAME
    charnum = index (ELEVATION_SOURCE_FILE,".dat")
    if (charnum /= len_trim(ELEVATION_SOURCE_FILE)-3)&
         & CALL WARNING("ELEVATION SOURCE FILE does not end in .dat", &
         & trim(ELEVATION_SOURCE_FILE))
    ! OPEN FILE
    pathnfile = trim(INPUT_DIR)//trim(ELEVATION_SOURCE_FILE)
    Call FOPEN(ELVUNIT,trim(pathnfile),'cfr')


  END SUBROUTINE OPEN_ELV_SOURCE

  SUBROUTINE OPEN_TS_SOURCE
    IMPLICIT NONE
    integer :: charnum
    logical :: back
    character(len=160) :: pathnfile

    back = .true.
    !Check TS SOURCE DATA FILE and open:
    ! TEST FILE NAME
    charnum = index (TS_SOURCE_FILE,".dat")
    if (charnum /= len_trim(TS_SOURCE_FILE)-3)&
         & CALL WARNING("TS SOURCE FILE does not end in .dat", &
         & trim(TS_SOURCE_FILE))
    ! OPEN FILE
    pathnfile = trim(INPUT_DIR)//trim(TS_SOURCE_FILE)
    Call FOPEN(TSUNIT,trim(pathnfile),'cfr')


  END SUBROUTINE OPEN_TS_SOURCE

  SUBROUTINE ALLOCATE_SPACE
    IMPLICIT NONE

    allocate(H(0:mgl)); H=0.0_sp
    allocate(vx(0:mgl)); vx=0.0_sp
    allocate(vy(0:mgl)); vy=0.0_sp

    ALLOCATE(Z(0:MGL,KB)); z=0.0_sp
    ALLOCATE(Z1(0:NGL,KB)); z1=0.0_sp
    ALLOCATE(ZZ(0:MGL,KB)); ZZ=0.0_sp
    ALLOCATE(ZZ1(0:NGL,KB)); ZZ1=0.0_sp
    ALLOCATE(DZ(0:MGL,KB)); DZ=0.0_SP
    ALLOCATE(DZ1(0:NGL,KB)); DZ1=0.0_SP
    ALLOCATE(DZZ(0:MGL,KB)); DZZ=0.0_SP
    ALLOCATE(DZZ1(0:NGL,KB)); DZZ1=0.0_SP


  END SUBROUTINE ALLOCATE_SPACE

  SUBROUTINE GET_OBC_TYPES
    USE MOD_OBCS
    IMPLICIT NONE

    IOBCN_GL=IOBCN
    IF(IOBCN > 0) THEN
       ALLOCATE(I_OBC_GL(IOBCN_GL))
       I_OBC_GL = I_OBC_N
       ALLOCATE(TYPE_OBC_GL(IOBCN_GL))
       TYPE_OBC_GL=TYPE_OBC

       CALL SETUP_OBCTYPES

    END IF
    

  END SUBROUTINE GET_OBC_TYPES

  SUBROUTINE SET_DIMENSIONS
    USE BCS
    IMPLICIT NONE
    INTEGER :: SBUF, SOURCE
    INTEGER :: RBUF, DEST
    

    DIM_siglev     => NC_MAKE_DIM(name='siglev',len=kb)
    DIM_siglay     => NC_MAKE_DIM(name='siglay',len=kbm1)
    
    DIM_DateStrLen => NC_MAKE_DIM(name='DateStrLen',len=DateStrLen)
    DIM_time       => NC_MAKE_DIM(name='time',len=NF90_UNLIMITED)


    DIM_nobc       => NC_MAKE_DIM(name='nobc',len=IOBCN)

    call print_dim(dim_nobc)

    ! MUST SET tidal Components dimension after reading source file!

  END SUBROUTINE SET_DIMENSIONS
  

  
  SUBROUTINE READ_SPECTRAL
    USE MOD_OBCS
    IMPLICIT NONE
    INTEGER :: I, JN, J, IOS
    
    INTEGER :: nCOLS
    character(LEN=80) :: component_names(100)

    INTEGER :: ITEMP, ISCAN
    CHARACTER(LEN=80) :: cmmnts
    INTEGER,  ALLOCATABLE :: NODE_SBC(:)
    !REAL(SP), ALLOCATABLE :: APT(:,:), PHAI(:,:), EMEAN(:)
    REAL(SP) :: TEMP

    ! READ THE NUMBER OF NODES SPECIFIED IN THE FILE
    ISCAN = SCAN_FILE(ELVUNIT,"Boundary Node Number",ISCAL = ITEMP)
    IF(ISCAN /= 0) then
       write(CMMNTS,'(I2)') ISCAN
       call fatal_error('Improper formatting of Spectral Tide Data File: ISCAN ERROR&
            &# '//trim(CMMNTS),&
            & 'The header must contain: "Boundary Node Number = "', &
            & 'Followed by an integer number of Nodes')
    END IF

    IF(IBCN_GL(1) /= ITEMP)CALL FATAL_ERROR&
         &('THE STATED NUMBER OF OPEN BOUNDARY POINTS WITH SET ELEVATION',&
         &'DOES NOT EQUAL THE NUMBER DETERMINED FROM THE OBC.DAT FILE',&
         &'YOU MAY SPECIFY APLITUDE AND PHASE FOR THE ENTIRE BOUNDARY HERE',&
         &'AND THEN MODIFY YOUR OBC.DAT FILE TO ONLY USE SOME OBC NODES')

!!$    ! CAN NOT READ COMMENTS BECAUSE IT WILL RETURN A VECTOR OF STRINGS...
!!$    ! READ ANY COMMENTS
!!$    ISCAN = SCAN_FILE(ELVUNIT,"COMMENTS",CVAL = cmmnts)
!!$    IF(ISCAN /=0) then
!!$       write(ipt,*) "! No comments found..."
!!$    else
!!$       write(ipt,*) "! Comments found..."
!!$    end IF


    COMMENTS = "Spectral forcing data from:"//TRIM(ELEVATION_SOURCE_FILE)


    ! READ THE COMPONENTS IN USE
    ISCAN = SCAN_FILE(ELVUNIT,"Tidal Constituents",NSZE=NCOLS, CVEC=component_names)
    IF(ISCAN /=0) THEN
       write(CMMNTS,'(I2)') ISCAN
       CALL FATAL_ERROR &
         &('THE NAMES OF THE TIDAL CONSTITUENTS IN EACH COLUMN MUST BE',&
         & 'INCLUDED IN THE HEADER OF THE SPECTRAL DATA INPUT FILE: ISCAN ERROR='//TRIM(CMMNTS),&
         & 'THE HEADER MUST CONTAIN: "Tidal Constituents ="',&
         & 'FOLLOWED BY A KNOWN TWO LETER ABREVIATION; eg S2, M2, K1, or NAN for an empty column')
    END IF

    IF(NCOLS <= 0) CALL FATAL_ERROR &
         &('THE NAMES OF THE TIDAL CONSTITUENTS IN EACH COLUMN MUST BE',&
         & 'INCLUDED IN THE HEADER OF THE SPECTRAL DATA INPUT FILE:',&
         & 'THE HEADER MUST CONTAIN: "Tidal Constituents ="',&
         & 'FOLLOWED BY A KNOWN TWO LETER ABREVIATION; eg S2, M2, K1, or NAN for an empty column')

    WRITE(IPT,*) "COUNT THE NUMBER OF COMPONENTS USED:"

    NCOMPS=0
    DO I = 1, NCOLS
       SELECT CASE(COMPONENT_NAMES(I))
       CASE("S2","s2")
          NCOMPS=NCOMPS+1
       CASE("M2","m2")
          NCOMPS=NCOMPS+1
       CASE("N2","n2")
          NCOMPS=NCOMPS+1
       CASE("K2","k2")
          NCOMPS=NCOMPS+1
       CASE("K1","k1")
          NCOMPS=NCOMPS+1
       CASE("P1","p1")
          NCOMPS=NCOMPS+1
       CASE("O1","o1")
          NCOMPS=NCOMPS+1
       CASE("Q1","q1")
          NCOMPS=NCOMPS+1
       CASE("NAN","nan")
          ! Do Nothing
       CASE DEFAULT
          CALL FATAL_ERROR("ILLEGAL VALUE::"//TRIM(COMPONENT_NAMES(I)),&
               &"IN THE TIDAL CONSTITUENTS HEADER ?")
       END SELECT
    END DO

    WRITE(IPT,*) "USING: ",NCOMPS,"; FROM:",NCOLS,"; COLUMNS OF DATA"

    IF (NCOMPS .LT. 1) CALL FATAL_ERROR &
         & ("SPECIFIED NO VALID TIDAL CONSTITUENT DATA?")


    !
    !----READ IN BOUNDARY POINTS, AMPLITUDES, AND PHASES OF TIDE-------------------|
    !

    ALLOCATE(NODE_SBC(IBCN_GL(1)), EMEAN(IBCN_GL(1)))
    ALLOCATE(APT(IBCN_GL(1),NCOLS), PHAI(IBCN_GL(1),NCOLS))
    APT = 0.0_SP ; PHAI = 0.0_SP ; EMEAN = 0.0_SP

    WRITE(IPT,*) "READING DATA"

    DO WHILE(.TRUE.)
       IF(IBCN_GL(1)==0) EXIT ! BREAK OUT IF THERE ARE NOT SUPPOSED
       ! TO BE ANY DATA POINTS
       READ(ELVUNIT,*,IOSTAT=IOS,END=99) J,temp
       if (IOS == 0) then
          BackSpace ELVUNIT
          exit
       end if
       
       CYCLE
       
99     Call FATAL_ERROR('Improper formatting of SPECTRAL DATA FILE:',&
            &'Reached end of file with out finding first line of data?',&
            &'FORMAT: NODE#,REFERENCE_HEIGHT (INT, FLOAT)')
       
    END DO


    DO I=1,IBCN_GL(1)
       READ(ELVUNIT,*)  NODE_SBC(I),EMEAN(I)
       READ (ELVUNIT,*) (APT(I,J), J=1,NCOLS)
       READ (ELVUNIT,*) (PHAI(I,J), J=1,NCOLS)

       WRITE(IPT,*) NODE_SBC(I),EMEAN(I)
       WRITE(IPT,*) (APT(I,J), J=1,NCOLS)
       WRITE(IPT,*) (PHAI(I,J), J=1,NCOLS)
    END DO


    APT = APT/100.0_SP ! CONVERT FROM CENTIMETERS TO METERS

    WRITE(IPT,*) "ALLOCATE SPACE"

! ALLOCATE SPACE FOR DATA

    ALLOCATE(EREF(IOBCN));          EREF=0.0_SP
    ALLOCATE(EPERIOD(NCOMPS));      EPERIOD=0.0_SP
    ALLOCATE(EAMP(IOBCN,NCOMPS));   EAMP=0.0_SP
    ALLOCATE(EPHASE(IOBCN,NCOMPS)); EPHASE=0.0_SP

! EQUILIBRIUM TIDE DATA
    ALLOCATE(EQI_AMP(NCOMPS))
    ALLOCATE(EQI_BETA(NCOMPS))
    ALLOCATE(TIDE_TYPE(NCOMPS))


    WRITE(IPT,*) "COMPARE NODES SPECIFIED WITH NODES READ AND TRANSFER REFERENCE HEIGHT"

     DO I=1,IBCN(1)
       JN = OBC_LST(1,I)
       IF(NODE_SBC(I) == I_OBC_N(JN))THEN
          EREF(JN) = EMEAN(I)
       ELSE
          CALL FATAL_ERROR &
               &('THE LIST OF SPECIFIED ELEVATION OPEN BOUNDARY NODES',&
               & 'DIFFERS BETWEEN INPUT DATA AND OBC.DAT',&
               & 'SPECIFIED NODES ARE THOSE WITH TYPE 1 OR 2')
       END IF
    END DO

     WRITE(IPT,*) "NOW TRANSFER CONSTITUENT DEPENDENT DATA:"

    COMPONENTS=""
    ITEMP=0
    DO I = 1, NCOLS

       SELECT CASE(COMPONENT_NAMES(I))
       CASE("NAN","nan")
          CYCLE
       CASE("S2","s2")
          ITEMP = ITEMP +1

          EPERIOD(ITEMP)= S2
          EQI_AMP(ITEMP) = S2_EQI_AMP
          EQI_BETA(ITEMP) = S2_EQI_BETA
          DO J=1,IBCN(1)
             JN = OBC_LST(1,J)
             EAMP(JN,ITEMP) = APT(J,I)
             EPHASE(JN,ITEMP) = PHAI(J,I)
          END DO

          TIDE_TYPE(ITEMP) = SEMIDIURNAL
          COMPONENTS= TRIM(COMPONENTS)//TRIM(COMPONENT_NAMES(I))

       CASE("M2","m2")
          ITEMP = ITEMP +1

          EPERIOD(ITEMP)= M2
          EQI_AMP(ITEMP) = M2_EQI_AMP
          EQI_BETA(ITEMP) = M2_EQI_BETA
          DO J=1,IBCN(1)
             JN = OBC_LST(1,J)
             EAMP(JN,ITEMP) = APT(J,I)
             EPHASE(JN,ITEMP) = PHAI(J,I)
          END DO

          TIDE_TYPE(ITEMP) = SEMIDIURNAL
          COMPONENTS= TRIM(COMPONENTS)//TRIM(COMPONENT_NAMES(I))

       CASE("N2","n2")
          ITEMP = ITEMP +1

          EPERIOD(ITEMP)= N2
          EQI_AMP(ITEMP) = N2_EQI_AMP
          EQI_BETA(ITEMP) = N2_EQI_BETA
          DO J=1,IBCN(1)
             JN = OBC_LST(1,J)
             EAMP(JN,ITEMP) = APT(J,I)
             EPHASE(JN,ITEMP) = PHAI(J,I)
          END DO

          TIDE_TYPE(ITEMP) = SEMIDIURNAL
          COMPONENTS= TRIM(COMPONENTS)//TRIM(COMPONENT_NAMES(I))

       CASE("K2","k2")
          ITEMP = ITEMP +1

          EPERIOD(ITEMP)= K2
          EQI_AMP(ITEMP) = K2_EQI_AMP
          EQI_BETA(ITEMP) = K2_EQI_BETA
          DO J=1,IBCN(1)
             JN = OBC_LST(1,J)
             EAMP(JN,ITEMP) = APT(J,I)
             EPHASE(JN,ITEMP) = PHAI(J,I)
          END DO

          TIDE_TYPE(ITEMP) = SEMIDIURNAL
          COMPONENTS= TRIM(COMPONENTS)//TRIM(COMPONENT_NAMES(I))

       CASE("K1","k1")
          ITEMP = ITEMP +1

          EPERIOD(ITEMP)= K1
          EQI_AMP(ITEMP) = K1_EQI_AMP
          EQI_BETA(ITEMP) = K1_EQI_BETA
          DO J=1,IBCN(1)
             JN = OBC_LST(1,J)
             EAMP(JN,ITEMP) = APT(J,I)
             EPHASE(JN,ITEMP) = PHAI(J,I)
          END DO

          TIDE_TYPE(ITEMP) = DIURNAL
          COMPONENTS= TRIM(COMPONENTS)//TRIM(COMPONENT_NAMES(I))

       CASE("P1","p1")
          ITEMP = ITEMP +1

          EPERIOD(ITEMP)= P1
          EQI_AMP(ITEMP) = P1_EQI_AMP
          EQI_BETA(ITEMP) = P1_EQI_BETA
          DO J=1,IBCN(1)
             JN = OBC_LST(1,J)
             EAMP(JN,ITEMP) = APT(J,I)
             EPHASE(JN,ITEMP) = PHAI(J,I)
          END DO

          TIDE_TYPE(ITEMP) = DIURNAL
          COMPONENTS= TRIM(COMPONENTS)//TRIM(COMPONENT_NAMES(I))

       CASE("O1","o1")
          ITEMP = ITEMP +1

          EPERIOD(ITEMP)= O1
          EQI_AMP(ITEMP) = O1_EQI_AMP
          EQI_BETA(ITEMP) = O1_EQI_BETA
          DO J=1,IBCN(1)
             JN = OBC_LST(1,J)
             EAMP(JN,ITEMP) = APT(J,I)
             EPHASE(JN,ITEMP) = PHAI(J,I)
          END DO

          TIDE_TYPE(ITEMP) = DIURNAL
          COMPONENTS= TRIM(COMPONENTS)//TRIM(COMPONENT_NAMES(I))

       CASE("Q1","q1")
          ITEMP = ITEMP +1

          EPERIOD(ITEMP)= Q1
          EQI_AMP(ITEMP) = Q1_EQI_AMP
          EQI_BETA(ITEMP) = Q1_EQI_BETA
          DO J=1,IBCN(1)
             JN = OBC_LST(1,J)
             EAMP(JN,ITEMP) = APT(J,I)
             EPHASE(JN,ITEMP) = PHAI(J,I)
          END DO

          TIDE_TYPE(ITEMP) = DIURNAL
          COMPONENTS= TRIM(COMPONENTS)//TRIM(COMPONENT_NAMES(I))

          ! Do Nothing
       END SELECT

       IF(ITEMP .LT. NCOMPS) COMPONENTS= TRIM(COMPONENTS)//","

    END DO
         
    WRITE(IPT,*) "FINISHED READING DATA FOR SPECTRAL FORCING"
    
  END SUBROUTINE READ_SPECTRAL
  

  SUBROUTINE WRITE_SPECTRAL
    use mod_set_time
    IMPLICIT NONE

    TYPE(NCFILE), POINTER :: NCF
    TYPE(NCVAR),  POINTER :: VAR, var1,var2
    TYPE(NCATT),  POINTER :: ATT
    TYPE(TIME) :: CURRENT
    INTEGER :: status
    CHARACTER(LEN=4) :: bflag
    LOGICAL :: FOUND

    DIM_tidal_components => NC_MAKE_DIM(name='tidal_components',len=size(Eperiod))


    CURRENT = GET_NOW()

    NCF => NEW_FILE()
    
    NCF%FNAME=TRIM(OUTPUT_DIR)//"spectral_obc.nc"


    ! ADD THE FILE ATTRIBUTES
    ATT => NC_MAKE_ATT(name='type',values="FVCOM SPECTRAL ELEVATION FORCING FILE") 
    NCF => ADD(NCF,ATT)

    ATT => NC_MAKE_ATT(name='title',values=TRIM(COMMENTS)) 
    NCF => ADD(NCF,ATT)

    ATT => NC_MAKE_ATT(name='components',values=TRIM(COMPONENTS)) 
    NCF => ADD(NCF,ATT)

    ATT => NC_MAKE_ATT(name='history',values="FILE CREATED: "&
         &//TRIM(WRITE_DATETIME(CURRENT,0,"UTC"))//": UTC") 
    NCF => ADD(NCF,ATT)
    

    IF(IOBCN .gt. 0) THEN
       VAR  => NC_MAKE_AVAR(name='obc_nodes', values=I_OBC_N, DIM1=DIM_nobc)
       
       ATT  => NC_MAKE_ATT(name='long_name',values='Open Boundary Node Number') 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='grid',values='obc_grid') 
       VAR  => ADD(VAR,ATT)
       
       NCF  => ADD(NCF,VAR)
    END IF

    VAR  => NC_MAKE_AVAR(name='tide_period', values=Eperiod, DIM1= DIM_tidal_components)
    
    ATT  => NC_MAKE_ATT(name='long_name',values='tide angular period') 
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='units',values='seconds') 
    VAR  => ADD(VAR,ATT)
    
    NCF  => ADD(NCF,VAR)
 
    IF(IOBCN .gt. 0) THEN
       VAR  => NC_MAKE_AVAR(name='tide_Eref', values=Eref, DIM1= DIM_nobc)
       
       ATT  => NC_MAKE_ATT(name='long_name',values='tidal elevation reference level') 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='units',values='meters') 
       VAR  => ADD(VAR,ATT)
       
       NCF  => ADD(NCF,VAR)
    END IF

    IF(IOBCN .gt. 0) THEN
       VAR  => NC_MAKE_AVAR(name='tide_Ephase', values=EPHASE, DIM1= DIM_nobc,DIM2=DIM_tidal_components)
       
       ATT  => NC_MAKE_ATT(name='long_name',values='tidal elevation phase angle') 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='units',values='degrees, time of maximum &
            &elevation with respect to chosen time origin') 
       VAR  => ADD(VAR,ATT)
       
       NCF  => ADD(NCF,VAR)
    END IF

    IF(IOBCN .gt. 0) THEN
       VAR  => NC_MAKE_AVAR(name='tide_Eamp', values=EAMP, DIM1= DIM_nobc,DIM2=DIM_tidal_components)
       
       ATT  => NC_MAKE_ATT(name='long_name',values='tidal elevation amplitude') 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='units',values='meters') 
       VAR  => ADD(VAR,ATT)
       
       NCF  => ADD(NCF,VAR)
    END IF

    VAR  => NC_MAKE_AVAR(name='equilibrium_tide_Eamp', values=EQI_AMP, DIM1= DIM_tidal_components)
    
    ATT  => NC_MAKE_ATT(name='long_name',values='equilibrium tidal elevation amplitude') 
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='units',values='meters') 
    VAR  => ADD(VAR,ATT)
    
    NCF  => ADD(NCF,VAR)

    VAR  => NC_MAKE_AVAR(name='equilibrium_beta_love', values=EQI_BETA, DIM1= DIM_tidal_components)
    
    ATT  => NC_MAKE_ATT(name='formula',values='beta=1+klove-hlove') 
    VAR  => ADD(VAR,ATT)
    
    NCF  => ADD(NCF,VAR)

    VAR  => NC_MAKE_AVAR(name='equilibrium_tide_type', &
         & values=TIDE_TYPE, DIM1=DIM_DateStrLen, DIM2= DIM_tidal_components)
    
    ATT  => NC_MAKE_ATT(name='long_name',values='formula') 
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='units',values='beta=1+klove-hlove') 
    VAR  => ADD(VAR,ATT)
    
    NCF  => ADD(NCF,VAR)


    if(USE_REAL_WORLD_TIME) then
       
       NOW = READ_DATETIME(time_origin,"YMD",TIMEZONE,status)
       IF(status /= 1) call fatal_error &
            &("could not parse time_origin or time_zone passed for spectral tidal forcing file?")
       
       VAR => DATETIME_OBJECT &
            &(DIMSTR=DIM_DateStrLen,&
            & timezone=timezone)
       ! OVER-RIDE DEFAULT NAME
       VAR%varname="time_origin"

       NCF  => ADD(NCF,VAR)
       
       VAR1 => FIND_VAR(NCF,"time_origin",FOUND)
       CALL UPDATE_DATETIME(VAR1,NOW)
       
    else

       CALL IDEAL_TIME_STRING2TIME(time_origin,BFLAG,NOW,IINT)
       IF(BFLAG == 'step') CALL FATAL_ERROR&
            &("You must specify a time, not a step, for this restart file", &
            & "The Step will be set by the old restart file...")


       ! float time
       VAR => FLOAT_TIME_OBJECT &
            &(Use_MJD=use_real_world_time)
       ! OVER-RIDE DEFAULT NAME
        VAR%varname="time_origin"

       NCF => ADD(NCF,VAR)
       
       VAR1 => FIND_VAR(NCF,"time_origin",FOUND)
       CALL UPDATE_FLOAT_TIME(VAR1,NOW)
       
    END IF
   
    CALL NC_WRITE_FILE(NCF)
    

  END SUBROUTINE WRITE_SPECTRAL

  SUBROUTINE READ_OLD_JULIAN
    USE MOD_SET_TIME
    USE MOD_OBCS
    IMPLICIT NONE
    INTEGER :: NCNT, IOS, ISCAN,j,CNT, i
    INTEGER :: status
    
    REAL(SP), ALLOCATABLE :: line(:)
    REAL(SP) :: tstep, tmp
    CHARACTER(LEN=80) TSTRING
    CHARACTER(LEN=80)      :: FMT, temp
    CHARACTER(LEN=160)      :: pathnfile
    CHARACTER(LEN=4) :: BFLAG

    if(USE_REAL_WORLD_TIME) then
       
       START = READ_DATETIME(time_origin,"YMD",TIMEZONE,status)
       IF(status /= 1) call fatal_error &
            &("could not parse time_origin or time_zone passed for spectral tidal forcing file?")
    else
       CALL IDEAL_TIME_STRING2TIME(time_origin,BFLAG,START,IINT)
       IF(BFLAG == 'step') CALL FATAL_ERROR&
            &("You must specify a time, not a step, for this restart file", &
            & "The Step will be set by the old restart file...")
    END if

    
    !----------------Determine Number of TIMES AND NODES -------------------------!
    
    ISCAN = SCAN_FILE(ELVUNIT,"Time Step",FSCAL = tstep)
    IF(ISCAN /= 0)then
       write(temp,'(I2)') ISCAN
       call fatal_error('Improper formatting of elj.dat file: ISCAN ERROR&
            &#'//trim(temp),&
            & 'The header must contain: "Time Step = "', &
            & 'Followed by a floating point number of seconds.')       
    END IF
    
    STEP = seconds2time(tstep)
    write(IPT,*) 

    ISCAN = SCAN_FILE(ELVUNIT,"Comments",CVAL = comments)
    IF(ISCAN /= 0)then
       Comments = "No comments found... this is mystery data!"
    END IF
    
    COMMENTS = "JULIAN FVCOM TIDAL FORCING DATA CREATED FROM OLD FILE &
         &TYPE: "//COMMENTS

    ! CLOSE AND REOPEN FILE WITH PAD=NO

    WRITE(IPT,*)"CLOSE AND REOPEN FILE TO TEST FOR CORRECT NUMBER OF SPECIFIED NODES"
    WRITE(IPT,*)"Number of nodes per line should match: ", IBCN(1)


    CLOSE(ELVUNIT)
    pathnfile = trim(INPUT_DIR)//trim(ELEVATION_SOURCE_FILE)
    OPEN(UNIT=ELVUNIT,FILE=trim(pathnfile),PAD='no',FORM='formatted',IOSTAT=ios)


    CNT = IBCN_GL(1)
    write(FMT,'(I5)') CNT
    FMT= '('//trim(FMT)//'E14.5)'

    NCNT = 0 
    rewind ELVUNIT
    DO WHILE(.TRUE.)
       NCNT= NCNT+1
       READ(ELVUNIT,FMT,IOSTAT=IOS,END=99) (tmp,J=1,CNT)
       if (IOS == 0) then
          BackSpace ELVUNIT
          ! TRY ADDING ONE MORE
          CNT =CNT +1
          write(FMT,'(I5)') CNT
          FMT= '('//trim(FMT)//'E14.5)'

          READ(ELVUNIT,FMT,IOSTAT=IOS,END=99) (tmp,J=1,(CNT))
          if (IOS /= 0) then
             BACKSPACE ELVUNIT
             exit
          else
             CALL FATAL_ERROR("THE DATA FILE DOES NOT APPEAR TO HAVE T&
                  &HE CORRECT NUMBER OF SPECIFIED NODES!")
          end if

       end if
       
       IF (NCNT > 50) CALL FATAL_ERROR &
            &('The elj.dat file does not have the correct number of nodes specified,',&
            &'or there are more than 50 lines of header?')

       CYCLE
       
99 Call FATAL_ERROR('Improper formatting of GRID FILE:',&
           &'Reached end of file with out finding CONNECTIVITY data?',&
           &'FORMAT: CELL# NODE# NODE# NODE# (ALL INTEGERS)')
          
    END DO

    ! NOW COUNT THE NUMBER OF RECORDS
    CNT = IBCN_GL(1)
    write(FMT,'(I5)') CNT
    FMT= '('//trim(FMT)//'E14.5)'
    NCNT = 0 

    DO WHILE(.TRUE.)
       READ(ELVUNIT,FMT,IOSTAT=IOS) (tmp,J=1,CNT)
       if (IOS == 0) then
          NCNT= NCNT+1
       else
          NTIMES = NCNT
          write(ipt,*) "FOUND",NCNT,"RECORDS!"
          exit
       end if
    END DO

    ALLOCATE(LINE(IBCN_GL(1) ))
    ALLOCATE(ELEVATION(IOBCN,NTIMES))

    ! FIND FIRST RECORD
    rewind(elvunit)
    DO WHILE(.TRUE.)
       READ(ELVUNIT,FMT,IOSTAT=IOS) (tmp,J=1,CNT)
       if (IOS == 0) then
          backspace elvunit
          exit
       end if
    END DO

    !READ RECORDS
    NCNT = 0
    DO WHILE(.TRUE.)
       
       READ(ELVUNIT,FMT,IOSTAT=IOS) line
       if (IOS /= 0) then
          exit
       end if

       ncnt = ncnt + 1

       if(ncnt>NTIMES) CALL FATAL_ERROR("ERROR READING ELEVATION FILE: TO MANY LINES IN FILE???")

       DO I=1,IBCN_GL(1)
          J = OBC_LST(1,I)
          ELEVATION(J,NCNT) = LINE(I)
       END DO

    END DO

    IF (NCNT/=NTIMES) THEN
       CALL FATAL_ERROR("ERROR READING ELEVATION FILE: BAD READ BEFORE END OF FILE")
    END IF



!--REPORT RESULTS--------------------------------------------------------------!
!

    WRITE(IPT,*)'!'
    WRITE(IPT,*)'!  JULIAN TIDE           :    SET'
    WRITE(IPT,*)'!  MAX AMPLITUDE =',maxval(ELEVATION)
    WRITE(IPT,*)'!  MIN AMPLITUDE =',minval(ELEVATION)
    call print_real_time(START,IPT,"Initial Time")
    NOW = start + (Ntimes-1)*STEP


    CALL PRINT_REAL_TIME(NOW, IPT,"END TIME")
    
    
  END SUBROUTINE READ_OLD_JULIAN

  SUBROUTINE WRITE_JULIAN
    USE MOD_NCDIO, only : update_IODATA
    use mod_set_time
    IMPLICIT NONE

    CHARACTER(LEN=8)  D
    CHARACTER(LEN=10) T

    TYPE(NCFILE), POINTER :: NCF
    TYPE(NCVAR),  POINTER :: VAR, VAR1
    TYPE(NCATT),  POINTER :: ATT

    real(SP), ALLOCATABLE :: LINE(:)
    INTEGER I
    INTEGER :: status
    CHARACTER(LEN=4) :: bflag
    LOGICAL :: FOUND

    ALLOCATE(LINE(IOBCN))

    CALL DATE_AND_TIME ( DATE=D,TIME=T )

    NCF => NEW_FILE()
    
    NCF%FNAME=TRIM(OUTPUT_DIR)//"julian_obc.nc"


    ! ADD THE FILE ATTRIBUTES
    ATT => NC_MAKE_ATT(name='type',values="FVCOM TIME SERIES ELEVATION FORCING FILE") 
    NCF => ADD(NCF,ATT)

    ATT => NC_MAKE_ATT(name='title',values=TRIM(COMMENTS)) 
    NCF => ADD(NCF,ATT)

#   if defined (EQUI_TIDE)
    ATT => NC_MAKE_ATT(name='components',values=TRIM(COMPONENTS)) 
    NCF => ADD(NCF,ATT)

    DIM_tidal_components => NC_MAKE_DIM(name='tidal_components',len=size(Eperiod))

    VAR  => NC_MAKE_AVAR(name='tide_period', values=Eperiod, DIM1= DIM_tidal_components)

    ATT  => NC_MAKE_ATT(name='long_name',values='tide angular period')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='seconds')
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    VAR  => NC_MAKE_AVAR(name='equilibrium_tide_Eamp', values=EQI_AMP, DIM1= DIM_tidal_components)

    ATT  => NC_MAKE_ATT(name='long_name',values='equilibrium tidal elevation amplitude')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='meters')
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    VAR  => NC_MAKE_AVAR(name='equilibrium_beta_love', values=EQI_BETA, DIM1= DIM_tidal_components)

    ATT  => NC_MAKE_ATT(name='formula',values='beta=1+klove-hlove')
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    VAR  => NC_MAKE_AVAR(name='equilibrium_tide_type', &
         & values=TIDE_TYPE, DIM1=DIM_DateStrLen, DIM2= DIM_tidal_components)

    ATT  => NC_MAKE_ATT(name='long_name',values='formula')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='beta=1+klove-hlove')
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    if(USE_REAL_WORLD_TIME) then

       NOW = READ_DATETIME(time_origin,"YMD",TIMEZONE,status)
       IF(status /= 1) call fatal_error &
            &("could not parse time_origin or time_zone passed for spectral tidal forcing file?")

       VAR => DATETIME_OBJECT &
            &(DIMSTR=DIM_DateStrLen,&
            & timezone=timezone)
       ! OVER-RIDE DEFAULT NAME
       VAR%varname="time_origin"

       NCF  => ADD(NCF,VAR)

       VAR1 => FIND_VAR(NCF,"time_origin",FOUND)
       CALL UPDATE_DATETIME(VAR1,NOW)

    else

       CALL IDEAL_TIME_STRING2TIME(time_origin,BFLAG,NOW,IINT)
       IF(BFLAG == 'step') CALL FATAL_ERROR&
            &("You must specify a time, not a step, for this restart file", &
            & "The Step will be set by the old restart file...")


       ! float time
       VAR => FLOAT_TIME_OBJECT &
            &(Use_MJD=use_real_world_time)
       ! OVER-RIDE DEFAULT NAME
        VAR%varname="time_origin"

       NCF => ADD(NCF,VAR)

       VAR1 => FIND_VAR(NCF,"time_origin",FOUND)
       CALL UPDATE_FLOAT_TIME(VAR1,NOW)

    END IF


#   endif

    ATT => NC_MAKE_ATT(name='history',values="FILE CREATED: "//D//"T"//T) 
    NCF => ADD(NCF,ATT)
    
    VAR  => NC_MAKE_AVAR(name='obc_nodes', values=I_OBC_N, DIM1= DIM_nobc)
    
    ATT  => NC_MAKE_ATT(name='long_name',values='Open Boundary Node Number') 
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='grid',values='obc_grid') 
    VAR  => ADD(VAR,ATT)
    
    NCF  => ADD(NCF,VAR)
    
    NCF => ADD(NCF, TIME_FILE_OBJECT() )

    ! OPEN BOUNDARY ELEVATION
    VAR  => NC_MAKE_AVAR(name='elevation', values=LINE, DIM1=DIM_nobc, DIM2=DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Open Boundary Elevation') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='meters') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)
    
    ALLOCATE(NCF%FTIME)
    NCF%FTIME%NEXT_STKCNT=0

    CALL NC_WRITE_FILE(NCF)

    NOW = START
    DO I=1,NTIMES
       
       CALL UPDATE_IODATA(NCF,NOW)
       
       LINE(:) = ELEVATION(:,I)
       NCF%FTIME%NEXT_STKCNT = NCF%FTIME%NEXT_STKCNT +1
       
       CALL NC_WRITE_FILE(NCF)
       
       NOW = NOW + STEP
    END DO

  RETURN
        
  END SUBROUTINE WRITE_JULIAN

  SUBROUTINE CREATE_SPECTRAL
    IMPLICIT NONE
    INTEGER :: NCOMPONENTS

    NCOMPONENTS = 2

    COMMENTS = "SPECTRAL TIDAL ELEVATION DATA: USER DEFINED"

    ALLOCATE(EREF(IOBCN));               EREF=0.0_SP
    ALLOCATE(EPERIOD(ncomponents));      EPERIOD=0.0_SP
    ALLOCATE(EAMP(IOBCN,ncomponents));   EAMP=0.0_SP
    ALLOCATE(EPHASE(IOBCN,ncomponents)); EPHASE=0.0_SP



    COMPONENTS = "some components"


  END SUBROUTINE CREATE_SPECTRAL

  SUBROUTINE CREATE_JULIAN
    IMPLICIT NONE
    INTEGER :: NCOMPONENTS, NTMP

    NTIMES = 100
    ! SET START TIME***
    COMMENTS = "JULIAN TIDAL ELEVATION DATA: USER DEFINED"

    ALLOCATE(ELEVATION(IOBCN,NTIMES)); ELEVATION=0.0_SP



  END SUBROUTINE CREATE_JULIAN


  SUBROUTINE READ_OLD_TS
    USE MOD_SET_TIME
    IMPLICIT NONE
    INTEGER :: NTMP, I, J, K,NCNT, ISCAN, STATUS
    INTEGER, ALLOCATABLE :: LYRS(:),NODES(:)
    REAL(SP):: FTEMP1
    CHARACTER(LEN=6) :: STRNG,TEMP
    CHARACTER(LEN=80) ::    TSTRING,fmt
    CHARACTER(LEN=4) :: BFLAG

    COMMENTS = "This data was transformed from an old FVCOM TS nudging file"

    if(USE_REAL_WORLD_TIME) then
       
       START = READ_DATETIME(time_origin,"YMD",TIMEZONE,status)
       IF(status /= 1) call fatal_error &
            &("could not parse time_origin or time_zone passed for spectral tidal forcing file?")
    else
       CALL IDEAL_TIME_STRING2TIME(time_origin,BFLAG,START,IINT)
       IF(BFLAG == 'step') CALL FATAL_ERROR&
            &("You must specify a time, not a step, for this restart file", &
            & "The Step will be set by the old restart file...")
    END if


    NCNT = 0     
    DO WHILE(.TRUE.)
       READ(TSUNIT,*,END=10) FTEMP1
       READ(TSUNIT,*)
       DO J=1,IOBCN
          READ(TSUNIT,*)
          READ(TSUNIT,*)
       ENDDO
       NCNT = NCNT + 1
    END DO
10  CONTINUE
    REWIND(TSUNIT)
    
    IF(NCNT == 0)CALL FATAL_ERROR("NO DATA PROVIDED FOR TEMPERATURE AND SALINITY OBC")
    
    NTIMES = NCNT 
    WRITE(IPT,*) "FOUND: ",NTIMES,"ENTRIES IN THE FILE:"


    !
    !----Read in Data Times and Global Heat Flux/Short Wave Radiation Data---------!
    !
    CALL ALLOCATE_TS
    ALLOCATE(LYRS(KBM1))
    ALLOCATE(NODES(IOBCN))


    ! FOR FOMAT STRING
    write(temp,'(I6)') KBM1
    DO J=1,NTIMES
       
       READ(TSUNIT,*) TIMES(J)               !time(days)
       WRITE(IPT,*) "READ FILE TIME: ",TIMES(J),"(days)"
       MJDS(J) = DAYS2TIME(TIMES(J))
       
       FMT = '(a6,'//trim(temp)//'I7)'
       READ(TSUNIT,FMT) strng, LYRS


       DO K=1,KBM1
          IF (LYRS(K)/=K) CALL FATAL_ERROR("TS OBC DATA LAYERS HEADER IS INCORECT")
       END DO

!       FMT = '(a6,'//trim(temp)//'f7.3)'
       DO I=1,IOBCN
          READ(TSUNIT,*) NODES(I),(OBC_TEMP(I,K,J),K=1,KBM1)   !temp in OB sigma layer
       ENDDO
       
       IF (ANY(NODES /= I_OBC_N)) CALL FATAL_ERROR("THE NODE NUMBER DOES NOT MATCH THE obc.dat FILE?")
       
       DO I=1,IOBCN
          READ(TSUNIT,*) NODES(I),(OBC_SALT(I,K,J),K=1,KBM1)   !sali in OB sigma layer
       ENDDO
       
       IF (ANY(NODES /= I_OBC_N)) CALL FATAL_ERROR("THE NODE NUMBER DOES NOT MATCH THE obc.dat FILE?")
       
    END DO

   


  END SUBROUTINE READ_OLD_TS


  SUBROUTINE ALLOCATE_TS
    IMPLICIT NONE


    ALLOCATE(OBC_TEMP(IOBCN,KBM1,NTIMES))
    ALLOCATE(OBC_SALT(IOBCN,KBM1,NTIMES))
    ALLOCATE(OBC_DEPTH(IOBCN,KBM1))
    ALLOCATE(OBC_H(IOBCN))
    ALLOCATE(OBC_X(IOBCN))
    ALLOCATE(OBC_Y(IOBCN))

    ALLOCATE(OBC_Z(IOBCN,KB))
    ALLOCATE(OBC_ZZ(IOBCN,KBM1))

    ALLOCATE(MJDS(NTIMES))
    ALLOCATE(TIMES(NTIMES))

  END SUBROUTINE ALLOCATE_TS


  SUBROUTINE SET_OBC_DEPTH
    IMPLICIT NONE
    INTEGER I,J

    IF(dbg_set(dbg_sbr)) WRITE(IPT,*) "Begin set_obc_depth"

    DO I = 1,IOBCN
       J = I_OBC_N(I)

       OBC_Z(I,1:KB) = Z(J,1:KB)
       OBC_ZZ(I,1:KBM1) = ZZ(J,1:KBM1)

       OBC_H(I) = H(j)
       OBC_X(I) = VX(j)
       OBC_Y(I) = VY(j)

       OBC_DEPTH(I,1:kbm1) =  ZZ(J,1:KBM1) * H(j) ! SSH = 0



    END DO

    IF(dbg_set(dbg_sbr)) WRITE(IPT,*) "End set_obc_depth"

  END SUBROUTINE SET_OBC_DEPTH


  SUBROUTINE CREATE_TS
    IMPLICIT NONE
    integer :: I, status

    ! SET A DUMMY NUMBER OF TIME POINTS
    NTIMES = 30

    ! SET A START DATE
    START = READ_DATETIME('2008-01-01 00:00:00',"YMD","UTC",status)
    IF(STATUS /= 1) CALL FATAL_ERROR("Bad Initial Time string in create_TS ?")

    CALL ALLOCATE_TS
    CALL SET_OBC_DEPTH

    ! SET A DUMMY TIME
    DO I = 1,NTIMES
       MJDS(I)%MJD=I
       MJDS(I)%MUSOD=0
       
    END DO

    
    ! WRITE SOME FUNCTION TO SET THE VALUE OF OBC_TEMP, OBC_SALT!

!    OBC_SALT = F(depth, x, y, time)
    
  END SUBROUTINE CREATE_TS



  SUBROUTINE WRITE_TSOBC
    USE MOD_NCDIO, only : UPDATE_IODATA
    IMPLICIT NONE

    CHARACTER(LEN=8)  D
    CHARACTER(LEN=10) T

    TYPE(NCFILE), POINTER :: NCF
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCATT),  POINTER :: ATT

    real(SP), ALLOCATABLE :: SALT(:,:)
    real(SP), ALLOCATABLE :: TEMP(:,:)
    INTEGER I

    IF(dbg_set(dbg_sbr)) WRITE(IPT,*) "Begin write_tsobc"

    ALLOCATE(SALT(IOBCN,KBM1))
    ALLOCATE(TEMP(IOBCN,KBM1))

    CALL DATE_AND_TIME ( DATE=D,TIME=T )

    NCF => NEW_FILE()
    
    NCF%FNAME=TRIM(OUTPUT_DIR)//"tsobc.nc"


    ! ADD THE FILE ATTRIBUTES
    ATT => NC_MAKE_ATT(name='type',values="FVCOM TIME SERIES OBC TS FILE") 
    NCF => ADD(NCF,ATT)

    ATT => NC_MAKE_ATT(name='title',values=TRIM(COMMENTS)) 
    NCF => ADD(NCF,ATT)

    ATT => NC_MAKE_ATT(name='history',values="FILE CREATED: "//D//"T"//T) 
    NCF => ADD(NCF,ATT)
    
    ! VARIABLES

    ! NOBC
    VAR  => NC_MAKE_AVAR(name='obc_nodes', values=I_OBC_N, DIM1= DIM_nobc)
    
    ATT  => NC_MAKE_ATT(name='long_name',values='Open Boundary Node Number') 
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='grid',values='obc_grid') 
    VAR  => ADD(VAR,ATT)
    
    NCF  => ADD(NCF,VAR)

    ! obc_h
    VAR  => NC_MAKE_AVAR(name='obc_h', values=obc_h, DIM1= DIM_nobc)
    
    ATT  => NC_MAKE_ATT(name='long_name',values='open boundary depth') 
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='grid',values='obc_grid') 
    VAR  => ADD(VAR,ATT)
    
    NCF  => ADD(NCF,VAR)
    
    ! obc siglay
    VAR  => NC_MAKE_AVAR(name='obc_siglay', values=obc_z, DIM1= DIM_nobc, DIM2=DIM_siglay)
    
    ATT  => NC_MAKE_ATT(name='long_name',values='ocean_sigma/general_coordinate') 
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='grid',values='obc_grid') 
    VAR  => ADD(VAR,ATT)
    
    NCF  => ADD(NCF,VAR)

    ! obc siglev
    VAR  => NC_MAKE_AVAR(name='obc_siglev', values=obc_zz, DIM1= DIM_nobc, DIM2=DIM_siglev)
    
    ATT  => NC_MAKE_ATT(name='long_name',values='ocean_sigma/general_coordinate') 
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='grid',values='obc_grid') 
    VAR  => ADD(VAR,ATT)
    
    NCF  => ADD(NCF,VAR)

    NCF => ADD(NCF, TIME_FILE_OBJECT() )
     
    ! obc temp
    VAR  => NC_MAKE_AVAR(name='obc_temp', values=TEMP, DIM1=DIM_nobc,&
         & DIM2=DIM_siglay,DIM3=DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='sea_water_temperature') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='Celcius') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='obc_grid') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)
 
    ! obc salt
    VAR  => NC_MAKE_AVAR(name='obc_salinity', values=SALT, DIM1=DIM_nobc,&
         & DIM2=DIM_siglay,DIM3=DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='sea_water_salinity') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='PSU') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='obc_grid') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)
    
    ALLOCATE(NCF%FTIME)
    NCF%FTIME%NEXT_STKCNT=0

    CALL NC_WRITE_FILE(NCF)

    NOW = START
    DO I=1,NTIMES

       NOW = START + MJDS(I)

!       CALL PRINT_REAL_TIME(NOW,IPT,"WRITING TIME")
       CALL UPDATE_IODATA(NCF,NOW)

       
       SALT(1:IOBCN,1:KBM1) = OBC_SALT(1:IOBCN,1:KBM1,I)
       TEMP(1:IOBCN,1:KBM1) = OBC_TEMP(1:IOBCN,1:KBM1,I)
       NCF%FTIME%NEXT_STKCNT = NCF%FTIME%NEXT_STKCNT +1
       
       CALL NC_WRITE_FILE(NCF)
       

    END DO

    IF(dbg_set(dbg_sbr)) WRITE(IPT,*) "end write_tsobc"

  END SUBROUTINE WRITE_TSOBC

  SUBROUTINE READ_TIDAL_COMPONENTS

    IMPLICIT NONE
    INTEGER :: I, JN, J, IOS
    
    INTEGER :: nCOLS
    character(LEN=80) :: component_names(100)

    INTEGER :: ITEMP, ISCAN
    CHARACTER(LEN=80) :: cmmnts
    INTEGER,  ALLOCATABLE :: NODE_SBC(:)
    !REAL(SP), ALLOCATABLE :: APT(:,:), PHAI(:,:), EMEAN(:)
    REAL(SP) :: TEMP

    COMMENTS = "Spectral forcing data from:"//TRIM(ELEVATION_SOURCE_FILE)


    ! READ THE COMPONENTS IN USE
    ISCAN = SCAN_FILE(ELVUNIT,"Tidal Constituents",NSZE=NCOLS, CVEC=component_names)
    IF(ISCAN /=0) THEN
       write(CMMNTS,'(I2)') ISCAN
       CALL FATAL_ERROR &
         &('THE NAMES OF THE TIDAL CONSTITUENTS IN EACH COLUMN MUST BE',&
         & 'INCLUDED IN THE HEADER OF THE SPECTRAL DATA INPUT FILE: ISCAN ERROR='//TRIM(CMMNTS),&
         & 'THE HEADER MUST CONTAIN: "Tidal Constituents ="',&
         & 'FOLLOWED BY A KNOWN TWO LETER ABREVIATION; eg S2, M2, K1, or NAN for an empty column')
    END IF

    IF(NCOLS <= 0) CALL FATAL_ERROR &
         &('THE NAMES OF THE TIDAL CONSTITUENTS IN EACH COLUMN MUST BE',&
         & 'INCLUDED IN THE HEADER OF THE SPECTRAL DATA INPUT FILE:',&
         & 'THE HEADER MUST CONTAIN: "Tidal Constituents ="',&
         & 'FOLLOWED BY A KNOWN TWO LETER ABREVIATION; eg S2, M2, K1, or NAN for an empty column')

    WRITE(IPT,*) "COUNT THE NUMBER OF COMPONENTS USED:"

    NCOMPS=0
    DO I = 1, NCOLS
       SELECT CASE(COMPONENT_NAMES(I))
       CASE("S2","s2")
          NCOMPS=NCOMPS+1
       CASE("M2","m2")
          NCOMPS=NCOMPS+1
       CASE("N2","n2")
          NCOMPS=NCOMPS+1
       CASE("K2","k2")
          NCOMPS=NCOMPS+1
       CASE("K1","k1")
          NCOMPS=NCOMPS+1
       CASE("P1","p1")
          NCOMPS=NCOMPS+1
       CASE("O1","o1")
          NCOMPS=NCOMPS+1
       CASE("Q1","q1")
          NCOMPS=NCOMPS+1
       CASE("NAN","nan")
          ! Do Nothing
       CASE DEFAULT
          CALL FATAL_ERROR("ILLEGAL VALUE::"//TRIM(COMPONENT_NAMES(I)),&
               &"IN THE TIDAL CONSTITUENTS HEADER ?")
       END SELECT
    END DO

    WRITE(IPT,*) "USING: ",NCOMPS,"; FROM:",NCOLS,"; COLUMNS OF DATA"

    IF (NCOMPS .LT. 1) CALL FATAL_ERROR &
         & ("SPECIFIED NO VALID TIDAL CONSTITUENT DATA?")


    ALLOCATE(EPERIOD(NCOMPS));      EPERIOD=0.0_SP
! EQUILIBRIUM TIDE DATA
    ALLOCATE(EQI_AMP(NCOMPS))
    ALLOCATE(EQI_BETA(NCOMPS))
    ALLOCATE(TIDE_TYPE(NCOMPS))

    COMPONENTS=""
    ITEMP=0
    DO I = 1, NCOLS

       SELECT CASE(COMPONENT_NAMES(I))
       CASE("NAN","nan")
          CYCLE
       CASE("S2","s2")
          ITEMP = ITEMP +1

          EPERIOD(ITEMP)= S2
          EQI_AMP(ITEMP) = S2_EQI_AMP
          EQI_BETA(ITEMP) = S2_EQI_BETA

          TIDE_TYPE(ITEMP) = SEMIDIURNAL
          COMPONENTS= TRIM(COMPONENTS)//TRIM(COMPONENT_NAMES(I))

       CASE("M2","m2")
          ITEMP = ITEMP +1

          EPERIOD(ITEMP)= M2
          EQI_AMP(ITEMP) = M2_EQI_AMP
          EQI_BETA(ITEMP) = M2_EQI_BETA

          TIDE_TYPE(ITEMP) = SEMIDIURNAL
          COMPONENTS= TRIM(COMPONENTS)//TRIM(COMPONENT_NAMES(I))

       CASE("N2","n2")
          ITEMP = ITEMP +1

          EPERIOD(ITEMP)= N2
          EQI_AMP(ITEMP) = N2_EQI_AMP
          EQI_BETA(ITEMP) = N2_EQI_BETA

          TIDE_TYPE(ITEMP) = SEMIDIURNAL
          COMPONENTS= TRIM(COMPONENTS)//TRIM(COMPONENT_NAMES(I))

       CASE("K2","k2")
          ITEMP = ITEMP +1

          EPERIOD(ITEMP)= K2
          EQI_AMP(ITEMP) = K2_EQI_AMP
          EQI_BETA(ITEMP) = K2_EQI_BETA

          TIDE_TYPE(ITEMP) = SEMIDIURNAL
          COMPONENTS= TRIM(COMPONENTS)//TRIM(COMPONENT_NAMES(I))

       CASE("K1","k1")
          ITEMP = ITEMP +1

          EPERIOD(ITEMP)= K1
          EQI_AMP(ITEMP) = K1_EQI_AMP
          EQI_BETA(ITEMP) = K1_EQI_BETA

          TIDE_TYPE(ITEMP) = DIURNAL
          COMPONENTS= TRIM(COMPONENTS)//TRIM(COMPONENT_NAMES(I))

       CASE("P1","p1")
          ITEMP = ITEMP +1

          EPERIOD(ITEMP)= P1
          EQI_AMP(ITEMP) = P1_EQI_AMP
          EQI_BETA(ITEMP) = P1_EQI_BETA

          TIDE_TYPE(ITEMP) = DIURNAL
          COMPONENTS= TRIM(COMPONENTS)//TRIM(COMPONENT_NAMES(I))

       CASE("O1","o1")
          ITEMP = ITEMP +1

          EPERIOD(ITEMP)= O1
          EQI_AMP(ITEMP) = O1_EQI_AMP
          EQI_BETA(ITEMP) = O1_EQI_BETA

          TIDE_TYPE(ITEMP) = DIURNAL
          COMPONENTS= TRIM(COMPONENTS)//TRIM(COMPONENT_NAMES(I))

       CASE("Q1","q1")
          ITEMP = ITEMP +1

          EPERIOD(ITEMP)= Q1
          EQI_AMP(ITEMP) = Q1_EQI_AMP
          EQI_BETA(ITEMP) = Q1_EQI_BETA

          TIDE_TYPE(ITEMP) = DIURNAL
          COMPONENTS= TRIM(COMPONENTS)//TRIM(COMPONENT_NAMES(I))

          ! Do Nothing
       END SELECT

       IF(ITEMP .LT. NCOMPS) COMPONENTS= TRIM(COMPONENTS)//","

    END DO
         
    WRITE(IPT,*) "FINISHED READING DATA OF TIDAL COMPONENTS"
    
  END SUBROUTINE READ_TIDAL_COMPONENTS
 


end module mod_obcreate
