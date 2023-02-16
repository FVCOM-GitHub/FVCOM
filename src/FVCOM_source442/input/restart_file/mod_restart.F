module mod_restart
  use all_vars
  use mod_utils
  use mod_par
  use mod_input
  use mod_prec
  use mod_ncll
  use mod_nctools
  USE MOD_TIME
  use bcs
  USE MOD_NCDIO
!! ggao
  USE CONTROL  
!--------------------------------------------------------------
# if defined (ICE)
  USE MOD_ICE
  USE MOD_ICE2D
# endif
!--------------------------------------------------------------

# if defined (DYE_RELEASE)
  USE mod_dye
# endif
  
  implicit none



  Character(Len=120):: FNAME
  Character(Len=120):: OLD_RESTART_FILE
  Character(Len=120):: OLD_WD_RESTART_FILE
  Character(Len=120):: OLD_ICE_RESTART_FILE
  INTEGER, PARAMETER :: RESUNIT = 101
  INTEGER, PARAMETER :: RESWDUNIT = 102
  
# if defined(ICE)
  INTEGER, PARAMETER :: ICERESTART = 201
  INTEGER            :: IINT_CHECK
  INTEGER,PARAMETER  :: NARR=5*NCAT ! NUMBER OF STATE VARIABLE ARRAYS
  INTEGER            :: NARRAYS     ! COUNTER FOR NUMBER OF STATE VARIABLE ARRAYS
  
  REAL (SP),ALLOCATABLE :: WORKS(:,:)
  REAL (SP),ALLOCATABLE :: WORKE(:,:)
  REAL (SP),ALLOCATABLE :: WORKF(:,:)
  REAL (SP),ALLOCATABLE :: FISICEC(:)
# endif  

  TYPE(NCFILE), POINTER :: NCF
  TYPE(GRID), SAVE :: MYGRID

  NAMELIST /NML_RST/                      &
       & INPUT_DIR,                       &
       & OUTPUT_DIR,                      &
       & START_DATE,                      &
       & TIMEZONE,                        &
       & DATE_FORMAT,                     &
       & PROJECTION_REFERENCE,            &
       & GRID_FILE_UNITS,                 &
       & GRID_FILE,                       &
       & SIGMA_LEVELS_FILE,               &
       & DEPTH_FILE,                      &
       & OBC_NODE_LIST_FILE,              &
       & CORIOLIS_FILE,                   &
       & SPONGE_FILE,                     &
       & OBC_LONGSHORE_FLOW_FILE,         &
       & OLD_RESTART_FILE,                &
       & OLD_ICE_RESTART_FILE,            &
       & OLD_WD_RESTART_FILE

  ! DATA VARIABLES

  TYPE(TIME), SAVE :: NOW
  LOGICAL CORIOLIS_FILE_ON

  REAL(SP), ALLOCATABLE :: RTP(:) ! THIS VARIABLE NO LONGER EXISTS BUT
  ! IS USED FOR BACKWARD COMPATIBLILITY


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
    WRITE(IPT,*) "!   INPUT_DIR, (Required)"
    WRITE(IPT,*) "!   OUTPUT_DIR, (Required)"
    WRITE(IPT,*) "!   START_DATE, (Required)"
    WRITE(IPT,*) "!   TIMEZONE, (Required)"
    WRITE(IPT,*) "!   GRID_FILE, (Required)"
    WRITE(IPT,*) "!   SIGMA_LEVELS_FILE, (Required)"
    WRITE(IPT,*) "!   DEPTH_FILE, (Required)"
    WRITE(IPT,*) "!   OBC_NODE_LIST_FILE, (Optional)"
    WRITE(IPT,*) "!   CORIOLIS_FILE, (Required if !spherical)"
    WRITE(IPT,*) "!   SPONGE_FILE, (Required)"
    WRITE(IPT,*) "!   OBC_LONGSHORE_FLOW_FILE, (Optional)"
    WRITE(IPT,*) "!   OLD_RESTART_FILE, (Required)"
    WRITE(IPT,*) "!   OLD_WD_RESTART_FILE, (Required if make wd)"
    WRITE(IPT,*) "!    "
    WRITE(IPT,*) "!    Example name list:"
    write(UNIT=IPT,NML=NML_RST)


    WRITE(IPT,*) "! NOTES: Do not run this program in parallel!"


  END SUBROUTINE MYHELPTXT

  SUBROUTINE READ_NAMELIST
    IMPLICIT NONE
    integer :: ios, i
    if(DBG_SET(dbg_sbr)) &
         & write(IPT,*) "Subroutine Begins: Read_Name_List;"


    if(DBG_SET(dbg_io)) &
         & write(IPT,*) "Read_Name_List: File: ",trim(FNAME)

    ! INITIALIZE SOME FIELDS
    INPUT_DIR = "NONE"
    OUTPUT_DIR = "NONE"
    START_DATE = "NONE"
    TIMEZONE = "NONE"
    DATE_FORMAT = "NONE"
    GRID_FILE = "NONE"
    SIGMA_LEVELS_FILE = "NONE"
    DEPTH_FILE = "NONE"
    OBC_NODE_LIST_FILE = "NONE"
    CORIOLIS_FILE = "NONE"
    SPONGE_FILE = "NONE"
    OBC_LONGSHORE_FLOW_FILE = "NONE"
    OLD_RESTART_FILE = "NONE"
    OLD_WD_RESTART_FILE = "NONE"

    CALL FOPEN(NMLUNIT,trim(FNAME),'cfr')

    !READ NAME LIST FILE

    ! Read IO Information
    READ(UNIT=NMLUNIT, NML=NML_RST,IOSTAT=ios)
    if(ios .NE. 0 ) then
       write(UNIT=IPT,NML=NML_RST)
       
       CALL FATAL_ERROR("Can Not Read NameList NML_RST from file: "//trim(FNAME))
    end if
    REWIND(NMLUNIT)
    
    write(IPT,*) "Read_Name_List:"
    
    write(UNIT=IPT,NML=NML_RST)
    
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

    ! LOOK FOR OBC FILE
    pathnfile = trim(INPUT_DIR)//trim(OBC_NODE_LIST_FILE)
    inquire(file=trim(pathnfile),exist=OBC_ON)

    IF(OBC_ON) THEN
       charnum = index (OBC_NODE_LIST_FILE,".dat")
       if (charnum /= len_trim(OBC_NODE_LIST_FILE)-3)&
            & CALL WARNING("OBC NODE LIST FILE does not end in .dat", &
            & trim(OBC_NODE_LIST_FILE))
       ! OPEN FILE
       pathnfile = trim(INPUT_DIR)//trim(OBC_NODE_LIST_FILE)
       Call FOPEN(OBCUNIT,trim(pathnfile),'cfr')


       ! LOOK FOR OBC LONG SHORE FLOW FILE
       pathnfile = trim(INPUT_DIR)//trim(OBC_LONGSHORE_FLOW_FILE)
       inquire(file=trim(pathnfile),exist=OBC_LONGSHORE_FLOW_ON)
       IF(OBC_LONGSHORE_FLOW_ON) THEN
          charnum = index (OBC_LONGSHORE_FLOW_FILE,".dat")
          if (charnum /= len_trim(OBC_LONGSHORE_FLOW_FILE)-3)&
               & CALL WARNING("OBC LONGSHORE FLOW FILE does not end in .dat", &
               & trim(OBC_LONGSHORE_FLOW_FILE))
          ! OPEN FILE
          pathnfile = trim(INPUT_DIR)//trim(OBC_LONGSHORE_FLOW_FILE)
          Call FOPEN(LSFUNIT,trim(pathnfile),'cfr')
       END IF
    END IF

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

    !Check Sponge File and open:
    ! TEST FILE NAME
    charnum = index (SPONGE_FILE,".dat")
    if (charnum /= len_trim(SPONGE_FILE)-3)&
         & CALL WARNING("SPONGE FILE does not end in .dat", &
         & trim(SPONGE_FILE))
    ! OPEN FILE
    pathnfile = trim(INPUT_DIR)//trim(SPONGE_FILE)
    Call FOPEN(SPONGEUNIT,trim(pathnfile),'cfr')

    !Check Coriolis File and open:
    ! TEST FILE NAME
# if !defined(SPHERICAL)
    CORIOLIS_FILE_ON=.TRUE.
    pathnfile = trim(INPUT_DIR)//trim(CORIOLIS_FILE)
    Call FOPEN(CORIOLISUNIT,trim(pathnfile),'cfr')
# endif


    Pathnfile = trim(INPUT_DIR)//trim(OLD_RESTART_FILE)
    Call FOPEN(RESUNIT,trim(pathnfile),'cur')

# if defined(WET_DRY)
    Pathnfile = trim(INPUT_DIR)//trim(OLD_WD_RESTART_FILE)
    Call FOPEN(RESWDUNIT,trim(pathnfile),'cfr')
# endif


# if defined(ICE)
    Pathnfile = trim(INPUT_DIR)//trim(OLD_ICE_RESTART_FILE)
    Call FOPEN(ICERESTART,trim(pathnfile),'cur')
# endif



  END SUBROUTINE OPEN_FILES

  SUBROUTINE ALLOCATE_SPACE
    USE MOD_WD
!--------------------------------------------------------------


    IMPLICIT NONE
    INTEGER :: STATUS

    ! VARIABLES IN RESTART FILE
    ALLOCATE(U(0:NGL,KB)); U= 0.0_SP
    ALLOCATE(V(0:NGL,KB)); V= 0.0_SP
    ALLOCATE(W(0:NGL,KB)); W= 0.0_SP

#    if defined (GOTM)
    ALLOCATE(TKE(0:MGL,KB)); TKE = 0.0_SP
    ALLOCATE(TEPS(0:MGL,KB)); TEPS = 0.0_SP
#    else
    ALLOCATE(Q2(0:MGL,KB)); Q2= 0.0_SP
    ALLOCATE(Q2L(0:MGL,KB)); Q2L=0.0_SP
    ALLOCATE(L(0:MGL,KB)); L = 0.0_SP
#    endif

    ALLOCATE(S(0:NGL,KB)); S=0.0_SP

    ALLOCATE(T(0:NGL,KB)); T=0.0_SP
    ALLOCATE(RHO(0:NGL,KB)); RHO=0.0_SP
    ALLOCATE(TMEAN(0:NGL,KB)); TMEAN=0.0_SP
    ALLOCATE(SMEAN(0:NGL,KB)); SMEAN=0.0_SP
    ALLOCATE(RMEAN(0:NGL,KB)); RMEAN=0.0_SP

    ALLOCATE(S1(0:MGL,KB)); S1=0.0_SP
    ALLOCATE(T1(0:MGL,KB)); T1=0.0_SP
    ALLOCATE(RHO1(0:MGL,KB)); RHO1=0.0_SP
    ALLOCATE(TMEAN1(0:MGL,KB)); TMEAN1=0.0_SP
    ALLOCATE(SMEAN1(0:MGL,KB)); SMEAN1=0.0_SP
    ALLOCATE(RMEAN1(0:MGL,KB)); RMEAN1=0.0_SP
    ALLOCATE(KM(0:MGL,KB)); KM= 0.0_SP
    ALLOCATE(KH(0:MGL,KB)); KH=0.0_SP
    ALLOCATE(KQ(0:MGL,KB)); KQ=0.0_SP

    ALLOCATE(UA(0:NGL)); UA = 0.0_SP
    ALLOCATE(VA(0:NGL)); VA= 0.0_SP
    ALLOCATE(EL1(0:NGL)); EL1 = 0.0_SP
    ALLOCATE(ET1(0:NGL)); ET1 = 0.0_SP
    ALLOCATE(H1(0:NGL)); H1= 0.0_SP
    ALLOCATE(D1(0:NGL)); D1= 0.0_SP
    ALLOCATE(DT1(0:NGL)); DT1 = 0.0_SP
    ALLOCATE(RTP(0:NGL)); RTP = 0.0_SP

    ALLOCATE(EL(0:MGL)); EL = 0.0_SP
    ALLOCATE(ET(0:MGL)); ET = 0.0_SP
    ALLOCATE(DT(0:MGL)); DT = 0.0_SP
    ALLOCATE(D(0:MGL)); D = 0.0_SP


    ! VARIABLE LOADED FROM COLD START FILES
    allocate(H(0:mgl)); H=0.0_sp
    allocate(vx(0:mgl)); vx=0.0_sp
    allocate(vy(0:mgl)); vy=0.0_sp

    allocate(xm(0:mgl)); xm=0.0_sp
    allocate(ym(0:mgl)); ym=0.0_sp

    allocate(lon(0:mgl)); lon=0.0_sp
    allocate(lat(0:mgl)); lat=0.0_sp

    ALLOCATE(KM1(0:NGL,KB)); KM1=0.0_sp

    ALLOCATE(COR(0:NGL)); COR=0.0_SP
    ALLOCATE(CC_SPONGE(0:NGL)); CC_SPONGE=0.0_SP
    ALLOCATE(XC(0:NGL)); XC=0.0_SP
    ALLOCATE(YC(0:NGL)); YC=0.0_SP
    ALLOCATE(LATC(0:NGL)); LATC=0.0_SP
    ALLOCATE(LONC(0:NGL)); LONC=0.0_SP
    ALLOCATE(XMC(0:NGL)); XMC=0.0_SP
    ALLOCATE(YMC(0:NGL)); YMC=0.0_SP


    ALLOCATE(Z(0:MGL,KB)); z=0.0_sp
    ALLOCATE(Z1(0:NGL,KB)); z1=0.0_sp
    ALLOCATE(ZZ(0:MGL,KB)); ZZ=0.0_sp
    ALLOCATE(ZZ1(0:NGL,KB)); ZZ1=0.0_sp
    ALLOCATE(DZ(0:MGL,KB)); DZ=0.0_SP
    ALLOCATE(DZ1(0:NGL,KB)); DZ1=0.0_SP
    ALLOCATE(DZZ(0:MGL,KB)); DZZ=0.0_SP
    ALLOCATE(DZZ1(0:NGL,KB)); DZZ1=0.0_SP

# if defined(WET_DRY)

       allocate(ISWETN(MGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE:ISWETN")
       ISWETN = 0

       allocate(ISWETC(NGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE:ISWETC")
       ISWETC = 0 

       allocate(ISWETNT(MGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE:ISWETNT")
       ISWETN = 0

       allocate(ISWETCT(NGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE:ISWETCT")
       ISWETCT = 0 

       allocate(ISWETCE(NGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE:ISWETCE")
       ISWETCE = 0
    
# endif

# if defined (EQUI_TIDE)
ALLOCATE(EL_EQI(MGL), stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE:EL_EQI")
       EL_EQI = 0.0_SP
# endif

# if defined (ATMO_TIDE)
       ALLOCATE(EL_ATMO(MGL), stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE:EL_ATMO")
       EL_ATMO = 0.0_SP
# endif
!=======================================================================
# if defined (ICE)
       !! ALLOCATE THE ICE HOTSTART VARIABLE
       JMT_LOCAL= 1

      WRITE(IPT,*)'ICE VARIABLES SIZE'
      WRITE(IPT,*)IMT_LOCAL,JMT_LOCAL,NCAT,NARR
      ALLOCATE(AICEN(M,JMT_LOCAL,NCAT));AICEN=0.0_SP  ! CONCENTRATION OF ICE
      ALLOCATE(VICEN(M,JMT_LOCAL,NCAT));VICEN=0.0_SP  ! VOLUME PER UNIT AREA OF ICE          (M)
      ALLOCATE(VSNON(M,JMT_LOCAL,NCAT));VSNON=0.0_SP  ! VOLUME PER UNIT AREA OF SNOW         (M)
      ALLOCATE(TSFCN(M,JMT_LOCAL,NCAT));TSFCN=0.0_SP  ! TEMPERATURE OF ICE/SNOW TOP SURFACE  (C)
      ALLOCATE(ESNON(M,JMT_LOCAL,NCAT));ESNON=0.0_SP  ! ENERGY OF MELT. OF SNOW LAYER    (J/M^2)

      ALLOCATE(AICE0(M,JMT_LOCAL))        ;AICE0=0.0_SP ! CONCENTRATION OF OPEN WATER
      ALLOCATE(EICEN(M,JMT_LOCAL,NTILAY)) ;EICEN=0.0_SP ! ENERGY OF MELTING FOR

      !-----------------------------------------------------------------
      ! OTHER VARIABLES CLOSELY RELATED TO THE STATE OF THE ICE
      !-----------------------------------------------------------------
      ALLOCATE(STRENGTH(M,JMT_LOCAL)); STRENGTH=0.0_SP ! ICE STRENGTH      

      ALLOCATE( FRESH  (1:M,1)); FRESH=0.0_SP    ! FRESH WATER FLUX TO OCEAN (KG/M^2/S)
      ALLOCATE( FSALT  (1:M,1)); FSALT=0.0_SP    ! SALT FLUX TO OCEAN (KG/M^2/S)
      ALLOCATE( FHNET  (1:M,1)); FHNET=0.0_SP    ! NET HEAT FLUX TO OCEAN (W/M^2)
      ALLOCATE( FSWTHRU(1:M,1));FSWTHRU=0.0_SP  ! SHORTWAVE PENETRATING TO OCEAN (W/M^2)
      ALLOCATE(WORKS(0:M,NARR))  ; WORKS=0.0_SP
      ALLOCATE(WORKE(0:M,NTILAY)); WORKE=0.0_SP
      ALLOCATE(WORKF(0:M,1))     ; WORKF=0.0_SP
      ALLOCATE(FISICEC(0:N))     ;FISICEC=0.0_SP
      ALLOCATE(ISICEC(0:N))      ;ISICEC=0
 
      ALLOCATE(UICE2(0:N))      ; UICE2 = 0.0_SP ! ICE VELOCITY
      ALLOCATE(VICE2(0:N))      ; VICE2 = 0.0_SP ! ICE VELOCITY

      ALLOCATE(SIG1(0:M))      ;  SIG1 = 0.0_SP  ! ICE STRESS TENSOR
      ALLOCATE(SIG2(0:M))      ;  SIG2 = 0.0_SP  ! ICE STRESS TENSOR 
      ALLOCATE(SIG12(0:M))     ; SIG12 = 0.0_SP  ! ICE STRESS TENSOR

# endif

!=======================================================================|

  END SUBROUTINE ALLOCATE_SPACE

  SUBROUTINE SET_OBC_TYPES
    USE MOD_OBCS
    IMPLICIT NONE

    IOBCN_GL=IOBCN
    if (IOBCN == 0 ) THEN
       OBC_ON = .False.
    ELSE
       ALLOCATE(I_OBC_GL(IOBCN_GL))
       ALLOCATE(I_OBC_N_OUTPUT(IOBCN_GL))
       I_OBC_GL = I_OBC_N
       I_OBC_N_OUTPUT = I_OBC_N
       
       ALLOCATE(TYPE_OBC_GL(IOBCN_GL))
       TYPE_OBC_GL=TYPE_OBC
       
       CALL SETUP_OBCTYPES
    END if

  END SUBROUTINE SET_OBC_TYPES


 SUBROUTINE READ_OLD_HOTSTART
   USE MOD_WD
   IMPLICIT NONE
   INTEGER :: I,J, K, INT, INTWD
   write(ipt,*) "Begin reading binary restart file:"

   write(ipt,*) "Reading: IINT"
   READ(RESUNIT) INT
   IINT = INT ! CONVERT TO LONG INTEGER

   write(ipt,*) "Reading: U,V,W"
   READ(RESUNIT) ((U(I,K),K=1,KB),I=0,N)
   READ(RESUNIT) ((V(I,K),K=1,KB),I=0,N)
   READ(RESUNIT) ((W(I,K),K=1,KB),I=0,N)


#    if defined (GOTM)
   write(ipt,*) "Reading: TKE,TEPS"
   READ(RESUNIT) ((TKE(I,K),K=1,KB),I=0,M)
   READ(RESUNIT) ((TEPS(I,K),K=1,KB),I=0,M)
#    else
   write(ipt,*) "Reading:Q2,Q2L,L"
   READ(RESUNIT) ((Q2(I,K),K=1,KB),I=0,M)
   READ(RESUNIT) ((Q2L(I,K),K=1,KB),I=0,M)
   READ(RESUNIT) ((L(I,K),K=1,KB),I=0,M)
#    endif


   write(ipt,*) "Reading:S,T,RHO"
   READ(RESUNIT) ((S(I,K),K=1,KB),I=0,N)
   READ(RESUNIT) ((T(I,K),K=1,KB),I=0,N)
   READ(RESUNIT) ((RHO(I,K),K=1,KB),I=0,N)

   write(ipt,*) "Reading: TMEAN,SMEAN,RMEAN"
   READ(RESUNIT) ((TMEAN(I,K),K=1,KB),I=0,N)
   READ(RESUNIT) ((SMEAN(I,K),K=1,KB),I=0,N)
   READ(RESUNIT) ((RMEAN(I,K),K=1,KB),I=0,N)
   
   write(ipt,*) "Reading: S1,T1,RHO1"
   READ(RESUNIT) ((S1(I,K),K=1,KB),I=1,M)
   READ(RESUNIT) ((T1(I,K),K=1,KB),I=1,M)
   READ(RESUNIT) ((RHO1(I,K),K=1,KB),I=1,M)
   write(ipt,*) "Reading: TMEAN1,SMEAN1,RMEAN1"
   READ(RESUNIT) ((TMEAN1(I,K),K=1,KB),I=1,M)
   READ(RESUNIT) ((SMEAN1(I,K),K=1,KB),I=1,M)
   READ(RESUNIT) ((RMEAN1(I,K),K=1,KB),I=1,M)
   write(ipt,*) "Reading:KM,KH,KQ"
   READ(RESUNIT) ((KM(I,K),K=1,KB),I=1,M)
   READ(RESUNIT) ((KH(I,K),K=1,KB),I=1,M)
   READ(RESUNIT) ((KQ(I,K),K=1,KB),I=1,M)
   
   write(ipt,*) "Reading:UA,VA,EL1,ET1"
   READ(RESUNIT) (UA(I), I=0,N)
   READ(RESUNIT) (VA(I), I=0,N)
   READ(RESUNIT) (EL1(I), I=1,N)
   READ(RESUNIT) (ET1(I), I=1,N)

   write(ipt,*) "Reading:H1,D1,DT1,RTP"
   READ(RESUNIT) (H1(I), I=1,N)
   READ(RESUNIT) (D1(I), I=1,N)
   READ(RESUNIT) (DT1(I), I=1,N)
   !    READ(RESUNIT) (DTF1(I), I=1,N)
   READ(RESUNIT) (RTP(I), I=1,N)

   write(ipt,*) "Reading:EL,ET,H,D,DT"   
   READ(RESUNIT) (EL(I), I=1,M)
   READ(RESUNIT) (ET(I), I=1,M)
   READ(RESUNIT) (H(I), I=1,M)
   READ(RESUNIT) (D(I), I=1,M)
   READ(RESUNIT) (DT(I), I=1,M)

   
#    if defined (EQUI_TIDE)
   write(ipt,*) "Reading: EL_EQI"
   READ(RESUNIT) (EL_EQI(I), I=1,M)
#    endif

#    if defined (ATMO_TIDE)
   write(ipt,*) "Reading:EL_ATMO"
   READ(RESUNIT) (EL_ATMO(I), I=1,M)
#    endif
   
#    if defined (WATER_QUALITY)
   DO N1=1,NB
   write(ipt,*) "Reading: WQM"

      READ(RESUNIT) ((WQM(I,K,N1),K=1,KB),I=1,M)
   END DO
#    endif
   
#    if defined (DYE_RELEASE)
!   IF(IINT.GT.DYESTART) THEN
   ! This switch would only read it if the dye has started, but I
   ! don't think the variable DYESTART is set in this program. 

      write(ipt,*) "Reading: DYE,DYEMEAN"

      READ(RESUNIT) ((DYE(I,K),K=1,KB),I=1,M)
!      READ(RESUNIT) ((DYEMEAN(I,K),K=1,KB),I=1,M)
!   ENDIF
#    endif
   
   WRITE(IPT,*) "FINISHED READING BINARY RESTART FILE"

   CLOSE(RESUNIT)

! ============================= READ WETDRY FILE TOO ===============================

# if defined(WET_DRY)

   READ(RESWDUNIT,*) INTWD
   READ(RESWDUNIT,*)
   IF(INTWD /= INT) call Fatal_error &
        & ("The WET DRY RESTART FILE IINT CYCLE DOES NOT MATCH THE MAIN RESTART FILE")

   READ(RESWDUNIT,*) (ISWETC(I), I=1,NGL)
   READ(RESWDUNIT,*) (ISWETN(I), I=1,MGL)

   CLOSE(RESWDUNIT)

   ! USE THE OLD METHOD TO INITIALIZE PREVIOUS STEP DATA...
   ISWETNT = ISWETN
   ISWETCT = ISWETC
   ISWETCE = ISWETC
# endif

# if defined(ICE)
   REWIND(ICERESTART)
   READ(ICERESTART) IINT_CHECK
   WRITE(IPT,*) "READING ICE AT",IINT_CHECK
   IF(IINT_CHECK /= IINT) THEN
      WRITE(IPT,*)'THE TIME OF ICE HOT START DIFFERENT FROM OCEAN MODEL'
      WRITE(IPT,*)'IINT =',IINT,'IINT_ICE =',IINT_CHECK
   END IF
   !-----------------------------------------------------------------
   ! STATE VARIABLES
   !-----------------------------------------------------------------
   READ(ICERESTART) ((WORKS(I,J),J=1,NARR),I=1,M) !NARR
   READ(ICERESTART) ((WORKE(I,J),J=1,NTILAY),I=1,M) !NTILAY
   !-----------------------------------------------------------------
   ! VELOCITY
   !-----------------------------------------------------------------
   
   READ(ICERESTART) (UICE2(I),I=0,N) !
   READ(ICERESTART) (VICE2(I),I=0,N)
   !-----------------------------------------------------------------
   ! FRESH WATER, SALT, AND HEAT FLUX
   !-----------------------------------------------------------------
   READ(ICERESTART) ((WORKF(I,J),J=1,1),I=1,M) !
   FRESH(1:M,1) =WORKF(1:M,1)
   READ(ICERESTART) ((WORKF(I,J),J=1,1),I=1,M) !
   FSALT(1:M,1) =WORKF(1:M,1)
   READ(ICERESTART) ((WORKF(I,J),J=1,1),I=1,M) !
   FHNET(1:M,1) =WORKF(1:M,1)
   !-----------------------------------------------------------------
   ! ICE STRENGTH
   !-----------------------------------------------------------------
   READ(ICERESTART) ((WORKF(I,J),J=1,1),I=1,M) !
   STRENGTH(1:M,1) =WORKF(1:M,1)
   !-----------------------------------------------------------------
   ! ICE MASK FOR DYNAMICS
   !-----------------------------------------------------------------
   READ(ICERESTART) (FISICEC(I),I=1,N)
   
   WRITE(IPT,*) "FINISHED ICE READING"
   CLOSE(ICERESTART)
   
   ISICEC(1:N) =ANINT(FISICEC(1:N))
   
   NARRAYS = 0  !1     ! IF AICE0 IS FIRST ARRAY 
   DO K=1,NCAT
      DO I = 1,M
         DO J = 1,JMT_LOCAL          
            AICEN(I,J,K) = WORKS(I,NARRAYS+1)
            VICEN(I,J,K) = WORKS(I,NARRAYS+2)
            VSNON(I,J,K) = WORKS(I,NARRAYS+3)
            TSFCN(I,J,K) = WORKS(I,NARRAYS+4)
            ESNON(I,J,K) = WORKS(I,NARRAYS+5)
         ENDDO
      ENDDO
      NARRAYS = NARRAYS + 5
   ENDDO
   
   DO K=1,NTILAY
      DO I = 1,M
         DO J = 1,JMT_LOCAL
            EICEN(I,J,K) = WORKE(I,K)
         ENDDO
      ENDDO
   ENDDO
   DEALLOCATE(WORKS,WORKE)
# endif

 END SUBROUTINE READ_OLD_HOTSTART
 
 FUNCTION RST_VERTICAL_VEL_FILE_OBJECT() RESULT(NCF)
   IMPLICIT NONE

    INTEGER :: status
    LOGICAL, SAVE :: IOPROC_ALLOCATED = .FALSE.
    TYPE(NCFILE), POINTER :: NCF
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCATT),  POINTER :: ATT


    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START RST_VERTICAL_VEL_FILE_OBJECT"

  ! ALLOCATE THE NEW FILE OBJECT
    NCF => NEW_FILE()

    VAR  => NC_MAKE_AVAR(name='w',&
         & values=w, DIM1= DIM_nele, DIM2= DIM_siglev, DIM3 = DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Vertical Sigma Coordinate Velocity') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='s-1') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! DIMS SHOULD NOW BE DEALLOCATED BUT I GUESS IT DOESN'T REALLY MATTER.

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END RST_VERTICAL_VEL_FILE_OBJECT"

  END FUNCTION RST_VERTICAL_VEL_FILE_OBJECT

 SUBROUTINE DUMP_RESTART
   IMPLICIT NONE

   IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START DUMP_RESTART"
    WRITE(IPT,*) "START DUMP_RESTART"

    CALL SET_FVCOM_GRID(MYGRID)
    CALL DEFINE_DIMENSIONS(MYGRID)
    


    write(ipt,*)'finished definition'

   NCF => NEW_FILE()

   ALLOCATE(NCF%FTIME)
   
   NCF%FNAME=trim(output_dir)//"restartfile.nc" 
   
   NCF => ADD(NCF,GRID_FILE_OBJECT(MYGRID) )

   NCF => ADD(NCF,TIME_FILE_OBJECT() )
   
   NCF => ADD(NCF,ZETA_FILE_OBJECT() )


   NCF => ADD(NCF,VELOCITY_FILE_OBJECT() )
   
   NCF => ADD(NCF,AVERAGE_VEL_FILE_OBJECT() )
   
   NCF => ADD(NCF,RST_VERTICAL_VEL_FILE_OBJECT() )
   
   NCF => ADD(NCF,TURBULENCE_FILE_OBJECT() )
   
   NCF => ADD(NCF,SALT_TEMP_FILE_OBJECT() )
   
   NCF => ADD(NCF,RESTART_EXTRAS_FILE_OBJECT() )
   
   IF(WETTING_DRYING_ON) THEN
      NCF => ADD(NCF, WET_DRY_FILE_OBJECT() )
   END IF

!!  ggao  0605/2008
# if defined (ICE)
      !-----------------------------------------------------------------
      ! state variables
    NCF => ADD(NCF,ICE_RESTART_STATE_FILE_OBJECT() )

      write(ipt,*)'Finished Ice State variables'
      !-----------------------------------------------------------------
      ! velocity
      !-----------------------------------------------------------------
    NCF => ADD(NCF,ICE_VEL_FILE_OBJECT() )
      write(ipt,*)'Finished Ice velocity'
      !-----------------------------------------------------------------
      ! fresh water, salt, and heat flux
      !-----------------------------------------------------------------
    NCF => ADD(NCF,ICE_EXTRA_FILE_OBJECT() )
      write(ipt,*)'Finished Ice water and heat budget'
# endif
   

!! change end


   ! WRITE THE STATIC VARIABLES
   CALL NC_WRITE_FILE(NCF)
   
   ! WRITE THE CURRENT STATE VARIABLES
   CALL UPDATE_IODATA(NCF,NOW)

   NCF%FTIME%NEXT_STKCNT =1
   CALL NC_WRITE_FILE(NCF)
 

   IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END DUMP_RESTART"
   
 END SUBROUTINE DUMP_RESTART


end module mod_restart

