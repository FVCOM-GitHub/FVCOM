module mod_grid
  use all_vars
  use mod_utils
  use mod_par
  use mod_input
  use mod_prec
  use mod_ncll
  use mod_nctools
  USE MOD_TIME
  USE MOD_NCDIO
  USE CONTROL  
  use mod_setup
  USE MOD_SET_TIME
  USE MOD_OBCS

  
  implicit none

  ! TIME OF DATA
  TYPE(TIME), SAVE :: NOW


  Character(Len=120):: FNAME
  Character(Len=120):: OLD_DATA_FILE
  
  TYPE(NCFILE), POINTER :: NCFIN
  TYPE(NCFILE), POINTER :: NCFOUT


  NAMELIST /NML_MODGRD/                   &
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
       & OLD_DATA_FILE


  ! DATA VARIABLES FOR NEW GRID
  INTEGER :: N_MGL
  INTEGER :: N_NGL
  INTEGER :: N_KB
  INTEGER :: N_KBM1
  
  INTEGER, ALLOCATABLE,TARGET :: N_NV(:,:)

  REAL(SP), ALLOCATABLE,TARGET :: N_H(:)

  REAL(SP), ALLOCATABLE,TARGET :: N_XC(:)
  REAL(SP), ALLOCATABLE,TARGET :: N_YC(:)

  REAL(SP), ALLOCATABLE,TARGET :: N_VX(:)
  REAL(SP), ALLOCATABLE,TARGET :: N_VY(:)

  REAL(SP), ALLOCATABLE,TARGET :: N_Z(:,:)
  REAL(SP), ALLOCATABLE,TARGET :: N_Z1(:,:)

  REAL(SP), ALLOCATABLE,TARGET :: N_ZZ(:,:)
  REAL(SP), ALLOCATABLE,TARGET :: N_ZZ1(:,:)

  INTEGER, ALLOCATABLE ::  N_NBE(:,:)   !!INDICES OF ELMNT NEIGHBORS
  INTEGER, ALLOCATABLE ::  N_NTSN(:)    !! NUMBER OF NODES SURROUNDING EACH NODE
  INTEGER, ALLOCATABLE ::  N_NBSN(:,:)            !! INDICES OF NODES SURROUNDING EACH NODE

  
  LOGICAL, ALLOCATABLE :: N_EFOUND(:) !! ELEMENT OF NEW MESH CONTAINED IN OLD MESH
  LOGICAL, ALLOCATABLE :: N_NFOUND(:) !! NODE OF NEW MESH CONTAINED IN OLD MESH

  INTEGER, ALLOCATABLE :: N2O_EID(:)  !! NUMBER OF NEAREST ELEMENT IN OLD MESH
  INTEGER, ALLOCATABLE :: N2O_NID(:)  !! NUMBER OF NEAREST ELEMENT IN OLD MESH


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
    WRITE(IPT,*) "!   OLD_DATA_FILE, (Required)"
    WRITE(IPT,*) "!    "
    WRITE(IPT,*) "!    Example name list:"
    write(UNIT=IPT,NML=NML_MODGRD)


    WRITE(IPT,*) "! NOTES: Do not run this program in parallel!"


  END SUBROUTINE MYHELPTXT

  SUBROUTINE INITIALIZE_NML
    IMPLICIT NONE

    ! INITIALIZE SOME FIELDS
    INPUT_DIR = "NONE"
    OUTPUT_DIR = "NONE"
    START_DATE = "NONE"
    TIMEZONE = "NONE"
    DATE_FORMAT = "NONE"
    GRID_FILE = "NONE"
    SIGMA_LEVELS_FILE = "NONE"
    DEPTH_FILE = "NONE"
    OLD_DATA_FILE = "NONE"

  END SUBROUTINE INITIALIZE_NML

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
    READ(UNIT=NMLUNIT, NML=NML_MODGRD,IOSTAT=ios)
    if(ios .NE. 0 ) then
       write(UNIT=IPT,NML=NML_MODGRD)
       
       CALL FATAL_ERROR("Can Not Read NameList NML_MODGRD from file: "//trim(FNAME))
    end if
    REWIND(NMLUNIT)
    
    write(IPT,*) "Read_Name_List:"
    
    write(UNIT=IPT,NML=NML_MODGRD)
    
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


    ! GET THE OLD DATA FILE
    pathnfile= trim(INPUT_DIR)//trim(OLD_DATA_FILE)
    NCF => NEW_FILE()
    NCF%FNAME=trim(pathnfile)

    Call NC_OPEN(NCF)
    CALL NC_LOAD(NCF)
       
    NCFIN => NCF

    ! OPEN THE NEW DATA FILE
    pathnfile= trim(INPUT_DIR)//"new_grid.nc"
    NCF => NEW_FILE()
    NCF%FNAME=trim(pathnfile)

    Call NC_CREATE(NCF)
    NCF%FTIME=>NEW_FTIME()
    NCFOUT => NCF
    

  END SUBROUTINE OPEN_FILES

  SUBROUTINE SET_TIME_INDEX
    IMPLICIT NONE
    INTEGER :: STATUS, I
    CHARACTER(LEN=4) :: BFLAG

    IF(USE_REAL_WORLD_TIME)THEN
       
       write(ipt,*) 'DATE_FORMAT' , DATE_FORMAT
       
       NOW = READ_DATETIME(START_DATE,DATE_FORMAT,TIMEZONE,status)
       IF(STATUS /= 1) CALL FATAL_ERROR&
            & ('Bad Start Date format!', &
            & TRIM(START_DATE))
    ELSE
       CALL IDEAL_TIME_STRING2TIME(START_DATE,BFLAG,NOW,IINT)
       IF(BFLAG == 'step') CALL FATAL_ERROR&
            &("You must specify a time, not a step, for this restart file", &
            & "The Step will be set by the old restart file...")
       
    END IF

    call print_time(NOW,ipt,"NOW")

    NC_START => NCFIN
    CALL SET_STARTUP_FILE_STACK(NOW,IINT)
    

    NULLIFY(NC_START)

  END SUBROUTINE SET_TIME_INDEX
  

  subroutine READ_NEW_GRID
    IMPLICIT NONE
    INTEGER :: NCT,I

    ! READ THE DATA FROM THE COLD START FILE USING MOD_MAIN VARAIBLE NAMES
    CALL READ_COLDSTART_SIGMA
    CLOSE(SIGMAUNIT)
    
    CALL READ_COLDSTART_GRID(GRIDUNIT,MGL,NGL,NV)
    m = MGL
    mt = MGL
    n = ngl
    nt = ngl
    NCT = NGL*3
    IOBCN = 0
    
        ! VARIABLE LOADED FROM COLD START FILES
    allocate(H(0:mgl)); H=0.0_sp
    allocate(vx(0:mgl)); vx=0.0_sp
    allocate(vy(0:mgl)); vy=0.0_sp

    allocate(xm(0:mgl)); xm=0.0_sp
    allocate(ym(0:mgl)); ym=0.0_sp

    allocate(lon(0:mgl)); lon=0.0_sp
    allocate(lat(0:mgl)); lat=0.0_sp

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
        

    ! VARS FOR TGE
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


    ALLOCATE(X_LCL(0:MGL),Y_LCL(0:MGL))

    CALL READ_COLDSTART_COORDS(GRIDUNIT,MGL,X_LCL,Y_LCL)
    CLOSE(GRIDUNIT)
    
    CALL COORDINATE_UNITS(X_LCL,Y_LCL)
    CALL SETUP_CENTER_COORDS
    
    CALL READ_COLDSTART_DEPTH(DEPTHUNIT,MGL,X_LCL,Y_LCL,H)
    CLOSE(DEPTHUNIT)
    
    DEALLOCATE(X_LCL,Y_LCL)

    KBM1 = KB - 1
    KBM2 = KB - 2
    CALL Setup_Sigma
    CALL SETUP_SIGMA_DERIVATIVES

    WRITE(IPT,*) "SETUP GRID NEIGHBORS"

    ALLOCATE(NGID(0:MGL))
    DO I =0,MGL
       NGID(I)=I
    END DO
    ALLOCATE(EGID(0:NGL))
    DO I =0,NGL
       EGID(I)=I
    END DO

    CALL TRIANGLE_GRID_EDGE

    DEALLOCATE(EGID)
    DEALLOCATE(NGID)

    WRITE(IPT,*) "TRANSFER TO NEW VARIABLES NAMES"
    
    ! MOVE DATA TO NEW DATA VARIABLE NAMES

    N_MGL  = MGL
    N_NGL  = NGL
    N_KB   = KB
    N_KBM1 = KBM1

    ALLOCATE(N_NV(0:N_NGL,3)); N_NV=0

        ! VARIABLE LOADED FROM COLD START FILES
    allocate(N_H(0:N_mgl));  N_H=0.0_sp
    allocate(N_vx(0:N_mgl)); N_vx=0.0_sp
    allocate(N_vy(0:N_mgl)); N_vy=0.0_sp

    ALLOCATE(N_XC(0:N_NGL)); N_XC=0.0_SP
    ALLOCATE(N_YC(0:N_NGL)); N_YC=0.0_SP

    ALLOCATE(N_Z(0:N_MGL,N_KB));   N_z=0.0_sp
    ALLOCATE(N_Z1(0:N_NGL,N_KB));  N_z1=0.0_sp
    ALLOCATE(N_ZZ(0:N_MGL,N_KB));  N_ZZ=0.0_sp
    ALLOCATE(N_ZZ1(0:N_NGL,N_KB)); N_ZZ1=0.0_sp

    ALLOCATE(N_NBE(0:NT,3))         ;N_NBE    = 0
    ALLOCATE(N_NTSN(MT))            ;N_NTSN   = 0 
    ALLOCATE(N_NBSN(SIZE(NBSN,1),SIZE(NBSN,2)))     ;N_NBSN =0

    N_NV = NV

    N_H  = H
    N_vx = VX + VXMIN
    N_vy = VY + VYMIN
 
    N_XC = XC + VXMIN
    N_YC = YC + VYMIN

    N_Z  = Z
    N_Z1 = Z1

    N_ZZ =ZZ
    N_ZZ1=ZZ1


    N_NBE  = NBE
    N_NTSN = NTSN
    N_NBSN = NBSN

    DEALLOCATE(NV)

    DEALLOCATE(H)

    DEALLOCATE(VX)
    DEALLOCATE(VY)
    DEALLOCATE(XM)
    DEALLOCATE(YM)
    DEALLOCATE(LAT)
    DEALLOCATE(LON)


    DEALLOCATE(XC)
    DEALLOCATE(YC)
    DEALLOCATE(XMC)
    DEALLOCATE(YMC)
    DEALLOCATE(LATC)
    DEALLOCATE(LONC)

    DEALLOCATE(Z)
    DEALLOCATE(Z1)

    DEALLOCATE(ZZ)
    DEALLOCATE(ZZ1)

    DEALLOCATE(DZ)
    DEALLOCATE(DZ1)

    DEALLOCATE(DZZ)
    DEALLOCATE(DZZ1)

    DEALLOCATE(NBE)
    DEALLOCATE(NTVE)
    DEALLOCATE(NTSN)
    DEALLOCATE(ISONB)
    DEALLOCATE(ISBCE)

    DEALLOCATE(NIEC)
    DEALLOCATE(NTRG)

    DEALLOCATE(XIJE)
    DEALLOCATE(YIJE)
    DEALLOCATE(DLTXE)
    DEALLOCATE(DLTYE)
    DEALLOCATE(DLTXYE)
    DEALLOCATE(SITAE)

    ! DELALLCATE MEMORY FROM TGE
    DEALLOCATE(NBSN)
    DEALLOCATE(NBVE)
    DEALLOCATE(NBVT)

    DEALLOCATE(IEC)
    DEALLOCATE(IENODE)
    DEALLOCATE(XIJC)
    DEALLOCATE(YIJC)
    DEALLOCATE(DLTXYC)
    DEALLOCATE(DLTXC)
    DEALLOCATE(DLTYC)
    DEALLOCATE(SITAC)    
    DEALLOCATE(ISBC)


    IF(ALLOCATED(LISBCE_1)) DEALLOCATE(LISBCE_1)
    IF(ALLOCATED(LISBCE_2)) DEALLOCATE(LISBCE_2)
    IF(ALLOCATED(LISBCE_3)) DEALLOCATE(LISBCE_3)

    DEALLOCATE(EPOR)

  end subroutine READ_NEW_GRID
  
  SUBROUTINE READ_OLD_GRID
    IMPLICIT NONE
    INTEGER :: NCT
    

    NC_START => NCFIN

    CALL LOAD_RESTART_GRID(NV)

    m = MGL
    mt = MGL
    n = ngl
    nt = ngl
    NCT = NT*3
    IOBCN = 0
    

    
    ! VARIABLE LOADED FROM THE OLD DATA FILE
    allocate(H(0:mgl)); H=0.0_sp
    allocate(vx(0:mgl)); vx=0.0_sp
    allocate(vy(0:mgl)); vy=0.0_sp

    allocate(xm(0:mgl)); xm=0.0_sp
    allocate(ym(0:mgl)); ym=0.0_sp

    allocate(lon(0:mgl)); lon=0.0_sp
    allocate(lat(0:mgl)); lat=0.0_sp

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



    ALLOCATE(Y_LCL(0:MT)); Y_LCL=0.0_SP
    ALLOCATE(X_LCL(0:MT)); X_LCL=0.0_SP
    ALLOCATE(H_LCL(0:MT)); H_LCL=0.0_SP
    

    ALLOCATE(ART(0:NT))           ;ART  = ZERO   !!AREA OF ELEMENT
    ALLOCATE(ART1(0:MT))          ;ART1 = ZERO   !!AREA OF NODE-BASE CONTROl VOLUME
    ALLOCATE(ART2(MT))            ;ART2 = ZERO   !!AREA OF ELEMENTS AROUND NODE
    
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

    !------------shape coefficient arrays and control volume metrics--------------------!
    
    ALLOCATE(A1U(0:NT,4))         ;A1U   = ZERO
    ALLOCATE(A2U(0:NT,4))         ;A2U   = ZERO 
    ALLOCATE(AWX(0:NT,3))         ;AWX   = ZERO 
    ALLOCATE(AWY(0:NT,3))         ;AWY   = ZERO 
    ALLOCATE(AW0(0:NT,3))         ;AW0   = ZERO 
    ALLOCATE(ALPHA(0:NT))         ;ALPHA = ZERO
    


    CALL LOAD_RESTART_COORDS(X_LCL,Y_LCL)
    CALL COORDINATE_UNITS(X_LCL,Y_LCL)
    CALL SETUP_CENTER_COORDS
    
    DEALLOCATE(X_LCL)
    DEALLOCATE(Y_LCL)
    
    VX = VX + VXMIN
    VY = VY + VYMIN

    XC = XC + VXMIN
    YC = YC + VYMIN


    CALL LOAD_RESTART_DEPTH(H_LCL)
    CALL SETUP_DEPTH
    DEALLOCATE(H_LCL) ! COULD BE LOADED DIRECTLY - MUST SET MAX/MIN
    
    STYPE = STYPE_RESTART
    CALL LOAD_RESTART_SIGMA(Z,Z1) ! LOAD DIRECTLY TO ALL_VARS:Z,Z1
    CALL SETUP_SIGMA_DERIVATIVES


  END SUBROUTINE READ_OLD_GRID

  SUBROUTINE MY_GRID_METRICS
    IMPLICIT NONE
    !============================================
    !Set up fluxes and control Volumes
    !============================================
    CALL TRIANGLE_GRID_EDGE
    
    !============================================
    !Calculate Element and Control Volume Areas
    !============================================
    CALL CELL_AREA            
    
    
    !====================================================
    ! Calculate Shape Coefficients for Flux Construction
    !====================================================
# if defined (GCN)
    CALL SHAPE_COEF_GCN       
# else
    CALL SHAPE_COEF_GCY  
# endif
    
    
    
  END SUBROUTINE MY_GRID_METRICS



  SUBROUTINE INTERP_COEFFICIENTS
    IMPLICIT NONE
    INTEGER :: LOOP_COUNT, I,J, missing,K
    real(SP) :: radlist(MGL)

    ! FIND THE NEW ELEMENTS IN THE OLD MESH

    WRITE(IPT,*) "FINDING LOCATION OF ALL NEW ELEMENTS"

    ALLOCATE(N_EFOUND(0:N_NGL)); N_EFOUND = .false.
    ALLOCATE(N2O_EID(0:N_NGL)); N2O_EID =0 

    DO I = 1,N_NGL
       
       N2O_EID(I) = Find_Element_Containing(N_XC(I),N_YC(I))
       
    END DO


    WHERE(N2O_EID /= 0)
       N_EFOUND = .true.
    END WHERE


    ! FIND ANY CELLS THAT HAVE NEIGHBORS IN THE OLD MESH
    DO I = 1,N_NGL
       
       IF(N_EFOUND(I)) CYCLE
       ! IF a neighbor is in the old mesh just copy the old mesh
       ! cells value
       IF(N_EFOUND(N_NBE(I,1))) THEN
          N2O_EID(I)=N2O_EID(N_NBE(I,1))
          CYCLE
       END IF

       IF(N_EFOUND(N_NBE(I,2))) THEN
          N2O_EID(I)=N2O_EID(N_NBE(I,2))
          CYCLE
       END IF

       IF(N_EFOUND(N_NBE(I,3))) THEN
          N2O_EID(I)=N2O_EID(N_NBE(I,3))
          CYCLE
       END IF

    END DO

    ! KEEP LOOKING UTILL EVERYONE IS FOUND
    LOOP_COUNT = 0
    DO WHILE (MINVAL(N2O_EID(1:N_NGL))==0)
       LOOP_COUNT = LOOP_COUNT+1
       
       missing = 0
       DO I = 1,N_NGL

          IF(N2O_EID(I)/=0) CYCLE

          missing = missing +1
          N2O_EID(I)=max(N2O_EID(N_NBE(I,1)),N2O_EID(I))

          N2O_EID(I)=max(N2O_EID(N_NBE(I,2)),N2O_EID(I))

          N2O_EID(I)=max(N2O_EID(N_NBE(I,3)),N2O_EID(I))
          
       END DO

       WRITE(IPT,*) "LOOP_COUNT: ",LOOP_COUNT, ": missing=",missing

    END DO

    WRITE(IPT,*) "FINDING LOCATION OF ALL NEW NODES"

    ALLOCATE(N_NFOUND(0:N_MGL)); N_NFOUND = .false.
    ALLOCATE(N2O_NID(0:N_MGL)); N2O_NID =0 

    DO I = 1,N_MGL
       
       N2O_NID(I) = Find_Element_Containing(N_VX(I),N_VY(I))
       
    END DO
    
    
    WHERE(N2O_NID /= 0)
       N_NFOUND = .true.
    END WHERE
    
      

    ! FIND ANY NODES THAT HAVE NEIGHBORS IN THE OLD MESH
    outer:DO I = 1,N_MGL
       
       IF(N_NFOUND(I)) CYCLE
       ! IF a neighboring node is in the old mesh just copy the old mesh
       ! nodes number
       
       inner:DO J=1,N_NTSN(I)
          IF(N_NFOUND(N_NBSN(I,J))) THEN
             K=N_NBSN(I,J)
             radlist = abs(vx(1:MGL) - n_vx(k)) + abs(VY(1:MGL) - n_VY(K))
             N2O_NID(I) = minloc(radlist,1)
             CYCLE outer
          END IF
       END DO inner
       
    END DO outer
    


    ! KEEP LOOKING UTILL EVERYONE IS FOUND
    LOOP_COUNT = 0
    DO WHILE (MINVAL(N2O_NID(1:N_MGL))==0)
       LOOP_COUNT = LOOP_COUNT+1

       missing = 0

       outer2:DO I = 1,N_MGL
          
          IF(N2O_NID(I)/=0) CYCLE

          missing = missing +1

          inner2:DO J=1,N_NTSN(I)
!             k=max(N2O_NID(N_NBSN(I,J)),N2O_NID(I))
             k=N2O_NID(N_NBSN(I,J))

             
             IF(K>0) THEN
!                radlist = abs(vx(1:MGL) - n_vx(k)) + abs(VY(1:MGL) - n_VY(K))

!                N2O_NID(I)= minloc(radlist,1)
                N2O_NID(I)=K
                CYCLE outer2
                
             END IF
             
          END DO inner2

          
       END DO outer2

       WRITE(IPT,*) "LOOP_COUNT: ",LOOP_COUNT, ": missing=", missing

    END DO    

  END SUBROUTINE INTERP_COEFFICIENTS


  SUBROUTINE CREATE_NEW_FILE
    USE SINTER
    IMPLICIT NONE
    INTEGER :: I,J,K
    INTEGER :: NVars, STKCNT,NDIMS,MYKB
    LOGICAL :: FOUND

    TYPE(NCVAR), POINTER :: VAR
    TYPE(NCVAR), POINTER :: N_VAR
    TYPE(NCDIM), POINTER :: DIM


!!$    INTEGER, POINTER                :: SCL_INT
!!$    INTEGER, POINTER,DIMENSION(:)   :: LVEC_INT
!!$    INTEGER, POINTER,DIMENSION(:,:) :: LARR_INT
!!$    INTEGER, POINTER,DIMENSION(:,:,:) :: LCUB_INT
!!$    
!!$    REAL(SPA), POINTER                :: SCL_FLT
!!$    REAL(SPA), POINTER,DIMENSION(:)   :: LVEC_FLT
!!$    REAL(SPA), POINTER,DIMENSION(:,:) :: LARR_FLT
!!$    REAL(SPA), POINTER,DIMENSION(:,:,:) :: LCUB_FLT
!!$    
!!$    REAL(DP), POINTER                :: SCL_DBL
!!$    REAL(DP), POINTER,DIMENSION(:)   :: LVEC_DBL
!!$    REAL(DP), POINTER,DIMENSION(:,:) :: LARR_DBL
!!$    REAL(DP), POINTER,DIMENSION(:,:,:) :: LCUB_DBL
!!$    
!!$    CHARACTER(LEN=80), POINTER        :: SCL_CHR
!!$    CHARACTER(LEN=80), POINTER        :: VEC_CHR(:)
    

  
    TYPE(NCVARP), pointer                       :: CURRENT

    REAL(SP), ALLOCATABLE :: XI(:),YI(:),X(:),Y(:)
    
    REAL(SPA) :: XSPA,YSPA,ZSPA
    REAL(DP) :: XDP,YDP,ZDP


    STKCNT = NCFIN%FTIME%PREV_STKCNT


    CURRENT  => NCFIN%VARS%NEXT
    FOUND = .FALSE.

    ! COPY CLOBAL ATTRIBUTES
    NCFOUT%ATTS = NCFIN%ATTS

    DO
       IF(.NOT. ASSOCIATED(CURRENT)) RETURN !END OF LIST
       
       VAR => CURRENT%VAR

       ! FOR DOUBLE PRECISION CASE CHANGE VARIABLE TYPE TO DOUBLE
# if defined(DOUBLE_PRECISION)
       IF(VAR%XTYPE==NF90_FLOAT) VAR%XTYPE=NF90_DOUBLE
# endif

       ! SOME SPECIAL VARIABLES - WE DON'T WANT TO INTERPOLATE!

       SELECT CASE(TRIM(VAR%VARNAME))
       CASE("nv")
          N_VAR => COPY_VAR(VAR)
          DIM => FIND_DIM(N_VAR,"nele",FOUND)
          DIM%DIM = N_NGL
          CALL NC_CONNECT_AVAR(N_VAR,N_NV)

          NCFOUT => ADD(NCFOUT, N_VAR)
          CURRENT  => CURRENT%NEXT          
          CYCLE
          
       CASE("h")
          N_VAR => COPY_VAR(VAR)
          DIM => FIND_DIM(N_VAR,"node",FOUND)
          DIM%DIM = N_MGL
          CALL NC_CONNECT_AVAR(N_VAR,N_H)

          NCFOUT => ADD(NCFOUT, N_VAR)
          CURRENT  => CURRENT%NEXT          
          CYCLE

       CASE("siglay")
          N_VAR => COPY_VAR(VAR)
          DIM => FIND_DIM(N_VAR,"node",FOUND)
          DIM%DIM = N_MGL

          DIM => FIND_DIM(N_VAR,"siglay",FOUND)
          DIM%DIM = N_KBM1
          CALL NC_CONNECT_AVAR(N_VAR,N_ZZ)

          NCFOUT => ADD(NCFOUT, N_VAR)
          CURRENT  => CURRENT%NEXT          
          CYCLE

       CASE("siglev")
          N_VAR => COPY_VAR(VAR)
          DIM => FIND_DIM(N_VAR,"node",FOUND)
          DIM%DIM = N_MGL

          DIM => FIND_DIM(N_VAR,"siglev",FOUND)
          DIM%DIM = N_KB
          CALL NC_CONNECT_AVAR(N_VAR,N_Z)

          NCFOUT => ADD(NCFOUT, N_VAR)
          CURRENT  => CURRENT%NEXT          
          CYCLE

# if defined(SPHERICAL)
       CASE("lon")
          N_VAR => COPY_VAR(VAR)
          DIM => FIND_DIM(N_VAR,"node",FOUND)
          DIM%DIM = N_MGL
          CALL NC_CONNECT_AVAR(N_VAR,N_VX)

          NCFOUT => ADD(NCFOUT, N_VAR)
          CURRENT  => CURRENT%NEXT          
          CYCLE


       CASE("lat")
          N_VAR => COPY_VAR(VAR)
          DIM => FIND_DIM(N_VAR,"node",FOUND)
          DIM%DIM = N_MGL
          CALL NC_CONNECT_AVAR(N_VAR,N_VY)

          NCFOUT => ADD(NCFOUT, N_VAR)
          CURRENT  => CURRENT%NEXT          
          CYCLE

       CASE("lonc")

          N_VAR => COPY_VAR(VAR)
          DIM => FIND_DIM(N_VAR,"nele",FOUND)
          DIM%DIM = N_NGL
          CALL NC_CONNECT_AVAR(N_VAR,N_XC)

          NCFOUT => ADD(NCFOUT, N_VAR)
          CURRENT  => CURRENT%NEXT          
          CYCLE

       CASE("latc")       

          N_VAR => COPY_VAR(VAR)
          DIM => FIND_DIM(N_VAR,"nele",FOUND)
          DIM%DIM = N_NGL
          CALL NC_CONNECT_AVAR(N_VAR,N_YC)

          NCFOUT => ADD(NCFOUT, N_VAR)
          CURRENT  => CURRENT%NEXT          
          CYCLE

       CASE("x")
!!$          N_VAR => COPY_VAR(VAR)
!!$          DIM => FIND_DIM(N_VAR,"node",FOUND)
!!$          DIM%DIM = N_MGL
!!$          CALL NC_CONNECT_AVAR(N_VAR,N_vx)
!!$
!!$          NCFOUT => ADD(NCFOUT, N_VAR)
          CURRENT  => CURRENT%NEXT          
          CYCLE

       CASE("y")
!!$          N_VAR => COPY_VAR(VAR)
!!$          DIM => FIND_DIM(N_VAR,"node",FOUND)
!!$          DIM%DIM = N_MGL
!!$          CALL NC_CONNECT_AVAR(N_VAR,N_vy)
!!$
!!$          NCFOUT => ADD(NCFOUT, N_VAR)
          CURRENT  => CURRENT%NEXT          
          CYCLE

       CASE("xc")
!!$          N_VAR => COPY_VAR(VAR)
!!$          DIM => FIND_DIM(N_VAR,"nele",FOUND)
!!$          DIM%DIM = N_NGL
!!$          CALL NC_CONNECT_AVAR(N_VAR,N_XC)
!!$
!!$          NCFOUT => ADD(NCFOUT, N_VAR)
          CURRENT  => CURRENT%NEXT          
          CYCLE

       CASE("yc")

!!$          N_VAR => COPY_VAR(VAR)
!!$          DIM => FIND_DIM(N_VAR,"nele",FOUND)
!!$          DIM%DIM = N_NGL
!!$          CALL NC_CONNECT_AVAR(N_VAR,N_YC)
!!$
!!$          NCFOUT => ADD(NCFOUT, N_VAR)
          CURRENT  => CURRENT%NEXT          
          CYCLE

# else
       CASE("x")
          N_VAR => COPY_VAR(VAR)
          DIM => FIND_DIM(N_VAR,"node",FOUND)
          DIM%DIM = N_MGL
          CALL NC_CONNECT_AVAR(N_VAR,N_vx)

          NCFOUT => ADD(NCFOUT, N_VAR)
          CURRENT  => CURRENT%NEXT          
          CYCLE

       CASE("y")
          N_VAR => COPY_VAR(VAR)
          DIM => FIND_DIM(N_VAR,"node",FOUND)
          DIM%DIM = N_MGL
          CALL NC_CONNECT_AVAR(N_VAR,N_vy)

          NCFOUT => ADD(NCFOUT, N_VAR)
          CURRENT  => CURRENT%NEXT          
          CYCLE

       CASE("xc")
          N_VAR => COPY_VAR(VAR)
          DIM => FIND_DIM(N_VAR,"nele",FOUND)
          DIM%DIM = N_NGL
          CALL NC_CONNECT_AVAR(N_VAR,N_XC)

          NCFOUT => ADD(NCFOUT, N_VAR)
          CURRENT  => CURRENT%NEXT          
          CYCLE

       CASE("yc")

          N_VAR => COPY_VAR(VAR)
          DIM => FIND_DIM(N_VAR,"nele",FOUND)
          DIM%DIM = N_NGL
          CALL NC_CONNECT_AVAR(N_VAR,N_YC)

          NCFOUT => ADD(NCFOUT, N_VAR)
          CURRENT  => CURRENT%NEXT          
          CYCLE

! allocate these
       CASE("lon")
!!$          N_VAR => COPY_VAR(VAR)
!!$          DIM => FIND_DIM(N_VAR,"node",FOUND)
!!$          DIM%DIM = N_MGL
!!$          CALL NC_CONNECT_AVAR(N_VAR,N_VX)
!!$
!!$          NCFOUT => ADD(NCFOUT, N_VAR)
          CURRENT  => CURRENT%NEXT          
          CYCLE


       CASE("lat")
!!$          N_VAR => COPY_VAR(VAR)
!!$          DIM => FIND_DIM(N_VAR,"node",FOUND)
!!$          DIM%DIM = N_MGL
!!$          CALL NC_CONNECT_AVAR(N_VAR,N_VY)
!!$
!!$          NCFOUT => ADD(NCFOUT, N_VAR)
          CURRENT  => CURRENT%NEXT          
          CYCLE

       CASE("lonc")

!!$          N_VAR => COPY_VAR(VAR)
!!$          DIM => FIND_DIM(N_VAR,"nele",FOUND)
!!$          DIM%DIM = N_NGL
!!$          CALL NC_CONNECT_AVAR(N_VAR,N_XC)
!!$
!!$          NCFOUT => ADD(NCFOUT, N_VAR)
          CURRENT  => CURRENT%NEXT          
          CYCLE

       CASE("latc")       

!!$          N_VAR => COPY_VAR(VAR)
!!$          DIM => FIND_DIM(N_VAR,"nele",FOUND)
!!$          DIM%DIM = N_NGL
!!$          CALL NC_CONNECT_AVAR(N_VAR,N_YC)
!!$
!!$          NCFOUT => ADD(NCFOUT, N_VAR)
          CURRENT  => CURRENT%NEXT          
          CYCLE


# endif

       END SELECT

       ! GET ANY DATA THAT SHOULD ME MOVED TO NEW NODES
       DIM => FIND_DIM(VAR,"node",FOUND)
       IF (FOUND) THEN
          write(ipt,*) "Nodal variable - varname="//TRIM(VAR%VARNAME)

          
          N_VAR => COPY_VAR(VAR)
          DIM => FIND_DIM(N_VAR,"node",FOUND)
          DIM%DIM = N_MGL

          DIM => FIND_DIM(N_VAR,"siglev",FOUND)
          IF(FOUND) THEN
             DIM%DIM = N_KB
             write(ipt,*) "DIM SIGLEV"
          END IF

          DIM => FIND_DIM(N_VAR,"siglay",FOUND)
          IF(FOUND) THEN
             DIM%DIM = N_KBM1
             write(ipt,*) "DIM SIGLAY"
          END IF

          CALL ALLOC_VAR(N_VAR,NDIMS)


          ! READ OLD GRID DATA
          CALL ALLOC_VAR(VAR)
          IF(has_unlimited(VAR)) THEN
             CALL NC_READ_VAR(VAR,STKCNT)
          ELSE
             CALL NC_READ_VAR(VAR)
          END IF
          
          ! INTERP TO NEW GRID
          select case(VAR%XTYPE)
          case(NF90_CHAR)
             CALL FATAL_ERROR("CAN NOT DO CHARACTER DATA!")
          CASE(NF90_BYTE)
             CALL FATAL_ERROR("CAN NOT DO BYTE DATA!")
          CASE(NF90_SHORT)
             CALL FATAL_ERROR("CAN NOT DO SHORT DATA!")
          CASE(NF90_INT)
             ! LEAVE IT EMPTY FOR NOW....
          case(NF90_FLOAT)
             IF(DBG_SET(DBG_SBRIO)) THEN
                WRITE(IPT,*) "XYTPE       :: FLOAT"
                WRITE(IPT,*) "ndims       ::",ndims
             END IF
             SELECT CASE(NDIMS)
             CASE(2)

                MYKB= size(N_VAR%ARR_FLT,2)
                

                IF(MYKB == N_KB )THEN                

                   ALLOCATE(XI(MYKB),YI(MYKB))
                   ALLOCATE(X(KB),Y(KB))


                   DO I = 1,N_MGL
                      k = N2O_NID(I)
                      IF(N_NFOUND(I)) THEN
                         
                         ! K is the element containing the point
                         XSPA = N_VX(I)
                         YSPA = N_VY(I)
                         DO J = 1,MYKB    
                            ZSPA = N_Z(I,J)
                            N_VAR%ARR_FLT(I,J)=INTERP_PNODAL(XSPA,YSPA&
                                 &,ZSPA,KB,k,VAR%ARR_FLT)
                         END DO
                      ELSE
                         
                         ! K is the nearest node
                         X = Z(k,1:KB)
                         Y = VAR%ARR_FLT(k,1:KB)

                         XI = N_Z(I,1:MYKB)

                         CALL SINTER_EXTRP_NONE(X,Y,XI,YI,KB,MYKB)

                         N_VAR%ARR_FLT(I,1:MYKB) = YI
                         

                      END IF
                   END DO

                   DEALLOCATE(XI,YI)
                   DEALLOCATE(X,Y)

                ELSEIF( MYKB == N_KBM1 )THEN 

                   ALLOCATE(XI(MYKB),YI(MYKB))
                   ALLOCATE(X(KBM1),Y(KBM1))

                   DO I = 1,N_MGL
                      k = N2O_NID(I)
                      IF(N_NFOUND(I)) THEN
                         ! K is the element containing the point
                         XSPA = N_VX(I)
                         YSPA = N_VY(I)
                         DO J = 1,MYKB    
                            ZSPA = N_ZZ(I,J)
                            N_VAR%ARR_FLT(I,J)=INTERP_PNODAL(XSPA,YSPA&
                                 &,ZSPA,KBM1,k,VAR%ARR_FLT)
                         END DO

!!$                         IF(I == 60000) THEN
!!$                            write(ipt,*) N_VX(I),N_VY(I)
!!$                            WRITE(IPT,*) VX(NV(K,:))
!!$                            WRITE(IPT,*) VY(NV(K,:))
!!$                            write(ipt,*)  N_VAR%ARR_FLT(I,:)
!!$                            write(ipt,*)  VAR%ARR_FLT(NV(K,3),:)
!!$
!!$                            write(ipt,*) N_ZZ(I,:)
!!$                            write(ipt,*) ZZ(NV(K,3),:)
!!$
!!$                            call pshutdown
!!$                         END IF

                      ELSE

                         ! K is the nearest node
                         X = Z(k,1:KBM1)
                         Y = VAR%ARR_FLT(k,1:KBM1)

                         XI = N_Z(I,1:MYKB)

                         CALL SINTER_EXTRP_NONE(X,Y,XI,YI,KBM1,MYKB)

                         N_VAR%ARR_FLT(I,1:MYKB) = YI

                      END IF
                   END DO

                   DEALLOCATE(XI,YI)
                   DEALLOCATE(X,Y)


                ELSE
                   CALL FATAL_ERROR("INVALID NUMBER OF LAYERS/LEVELS ALLOCATED!")
                END IF

                
             CASE(1)

                DO I = 1,N_MGL
                   
                   IF(N_NFOUND(I))THEN
                      
                      XSPA = N_VX(I)
                      YSPA = N_VY(I)
                      N_VAR%VEC_FLT(I)=INTERP_PNODAL(XSPA,YSPA&
                           &,N2O_NID(I),VAR%VEC_FLT)
                      
                   ELSE
                      
                      N_VAR%VEC_FLT(I)= VAR%VEC_FLT(N2O_NID(I))
                      
                   END IF
                   
                END DO
                
                
             CASE DEFAULT
                CALL FATAL_ERROR("BAD NUMBER OF DIMENSION!")
             END SELECT
             
          case(NF90_DOUBLE)

             IF(DBG_SET(DBG_SBRIO)) THEN
                WRITE(IPT,*) "XYTPE       :: DOUBLE"
                WRITE(IPT,*) "ndims       ::",ndims
             END IF
             SELECT CASE(NDIMS)
             CASE(2)

                MYKB= size(N_VAR%ARR_DBL,2)
                

                IF(MYKB == N_KB )THEN                

                   ALLOCATE(XI(MYKB),YI(MYKB))
                   ALLOCATE(X(KB),Y(KB))


                   DO I = 1,N_MGL
                      k = N2O_NID(I)
                      IF(N_NFOUND(I)) THEN
                         
                         XDP = N_VX(I)
                         YDP = N_VY(I)
                         ! K is the element containing the point
                         DO J = 1,MYKB  
                            ZDP = N_Z(I,J)
                            N_VAR%ARR_DBL(I,J)=INTERP_PNODAL(XDP,YDP&
                                 &,ZDP,KB,k,VAR%ARR_DBL)
                         END DO
                      ELSE
                         
                         ! K is the nearest node
                         X = Z(k,1:KB)
                         Y = VAR%ARR_DBL(k,1:KB)

                         XI = N_Z(I,1:MYKB)

                         CALL SINTER_EXTRP_NONE(X,Y,XI,YI,KB,MYKB)

                         N_VAR%ARR_DBL(I,1:MYKB) = YI
                         

                      END IF
                   END DO

                   DEALLOCATE(XI,YI)
                   DEALLOCATE(X,Y)

                ELSEIF( MYKB == N_KBM1 )THEN 

                   ALLOCATE(XI(MYKB),YI(MYKB))
                   ALLOCATE(X(KBM1),Y(KBM1))

                   DO I = 1,N_MGL
                      k = N2O_NID(I)
                      IF(N_NFOUND(I)) THEN

                         XDP = N_VX(I)
                         YDP = N_VY(I)
                         ! K is the element containing the point
                         DO J = 1,MYKB  
                            ZDP = N_ZZ(I,J)
                            N_VAR%ARR_DBL(I,J)=INTERP_PNODAL(XDP,YDP&
                                 &,ZDP,KBM1,k,VAR%ARR_DBL)
                         END DO

!!$                         IF(I == 60000) THEN
!!$                            write(ipt,*) N_VX(I),N_VY(I)
!!$                            WRITE(IPT,*) VX(NV(K,:))
!!$                            WRITE(IPT,*) VY(NV(K,:))
!!$                            write(ipt,*)  N_VAR%ARR_FLT(I,:)
!!$                            write(ipt,*)  VAR%ARR_FLT(NV(K,3),:)
!!$
!!$                            write(ipt,*) N_ZZ(I,:)
!!$                            write(ipt,*) ZZ(NV(K,3),:)
!!$
!!$                            call pshutdown
!!$                         END IF

                      ELSE

                         ! K is the nearest node
                         X = ZZ(k,1:KBM1)
                         Y = VAR%ARR_DBL(k,1:KBM1)

                         XI = N_ZZ(I,1:MYKB)

                         CALL SINTER_EXTRP_NONE(X,Y,XI,YI,KBM1,MYKB)

                         N_VAR%ARR_DBL(I,1:MYKB) = YI

                      END IF
                   END DO

                   DEALLOCATE(XI,YI)
                   DEALLOCATE(X,Y)


                ELSE
                   CALL FATAL_ERROR("INVALID NUMBER OF LAYERS/LEVELS ALLOCATED!")
                END IF

                
             CASE(1)

                DO I = 1,N_MGL
                   
                   IF(N_NFOUND(I))THEN
                      
                      XDP = N_VX(I)
                      YDP = N_VY(I)
                      N_VAR%VEC_DBL(I)=INTERP_PNODAL(XDP,YDP&
                           &,N2O_NID(I),VAR%VEC_DBL)
                      
                   ELSE
                      
                      N_VAR%VEC_DBL(I)= VAR%VEC_DBL(N2O_NID(I))
                      
                   END IF
                   
                END DO
                
                
             CASE DEFAULT
                CALL FATAL_ERROR("BAD NUMBER OF DIMENSION!")
             END SELECT

          END select

          NCFOUT => ADD(NCFOUT, N_VAR)

          CURRENT  => CURRENT%NEXT          
          CYCLE
       END IF
       

       ! GET ANY DATA THAT SHOULD ME MOVED TO NEW ELEMENTS
       DIM => FIND_DIM(VAR,"nele",FOUND)
       IF (FOUND) THEN
          write(ipt,*) "CELL variable - varname="//TRIM(VAR%VARNAME)

          N_VAR => COPY_VAR(VAR)
          DIM => FIND_DIM(N_VAR,"nele",FOUND)
          DIM%DIM = N_NGL

          DIM => FIND_DIM(N_VAR,"siglev",FOUND)
          IF(FOUND) THEN
             DIM%DIM = N_KB
             write(ipt,*) "DIM SIGLEV"
          END IF
          DIM => FIND_DIM(N_VAR,"siglay",FOUND)
          IF(FOUND) THEN
             DIM%DIM = N_KBM1
             write(ipt,*) "DIM SIGLAY"
          END IF

          CALL ALLOC_VAR(N_VAR,NDIMS)


          ! READ OLD GRID DATA
          CALL ALLOC_VAR(VAR)
          IF(has_unlimited(VAR)) THEN
             CALL NC_READ_VAR(VAR,STKCNT)
          ELSE
             CALL NC_READ_VAR(VAR)
          END IF
          
          ! INTERP TO NEW GRID
          select case(VAR%XTYPE)
          case(NF90_CHAR)
             CALL FATAL_ERROR("CAN NOT DO CHARACTER DATA!")
          CASE(NF90_BYTE)
             CALL FATAL_ERROR("CAN NOT DO BYTE DATA!")
          CASE(NF90_SHORT)
             CALL FATAL_ERROR("CAN NOT DO SHORT DATA!")
          CASE(NF90_INT)
             ! LEAVE IT EMPTY FOR NOW....
          case(NF90_FLOAT)
             IF(DBG_SET(DBG_SBRIO)) THEN
                WRITE(IPT,*) "XYTPE       :: FLOAT"
                WRITE(IPT,*) "ndims       ::",ndims
             END IF
             SELECT CASE(NDIMS)
             CASE(2)

                MYKB= size(N_VAR%ARR_FLT,2)
                

                IF(MYKB == N_KB )THEN                

                   ALLOCATE(XI(MYKB),YI(MYKB))
                   ALLOCATE(X(KB),Y(KB))


                   DO I = 1,N_NGL
                      k = N2O_EID(I)
                      IF(N_EFOUND(I)) THEN
                         
                         XSPA = N_XC(I)
                         YSPA = N_YC(I)
                         
                         DO J = 1,MYKB     
                            ZSPA = N_Z1(I,J)
                            N_VAR%ARR_FLT(I,J)=INTERP_PZONAL(XSPA,YSPA&
                                 &,ZSPA,KB,k,VAR%ARR_FLT)
                         END DO
                      ELSE
                         
                         X = Z1(k,1:KB)
                         Y = VAR%ARR_FLT(k,1:KB)

                         XI = N_Z1(I,1:MYKB)

                         CALL SINTER_EXTRP_NONE(X,Y,XI,YI,KB,MYKB)

                         N_VAR%ARR_FLT(I,1:MYKB) = YI
                         

                      END IF
                   END DO

                   DEALLOCATE(XI,YI)
                   DEALLOCATE(X,Y)

                ELSEIF( MYKB == N_KBM1 )THEN 

                   ALLOCATE(XI(MYKB),YI(MYKB))
                   ALLOCATE(X(KBM1),Y(KBM1))

                   DO I = 1,N_NGL
                      k = N2O_EID(I)
                      IF(N_EFOUND(I)) THEN

                         XSPA = N_XC(I)
                         YSPA = N_YC(I)
                         
                         DO J = 1,MYKB     
                            ZSPA = N_ZZ1(I,J)
                            N_VAR%ARR_FLT(I,J)=INTERP_PZONAL(XSPA,YSPA&
                                 &,ZSPA,KBM1,k,VAR%ARR_FLT)
                         END DO

                      ELSE

                         X = ZZ1(k,1:KBM1)
                         Y = VAR%ARR_FLT(k,1:KBM1)

                         XI = N_ZZ1(I,1:MYKB)

                         CALL SINTER_EXTRP_NONE(X,Y,XI,YI,KBM1,MYKB)

                         N_VAR%ARR_FLT(I,1:MYKB) = YI

                      END IF
                   END DO

                   DEALLOCATE(XI,YI)
                   DEALLOCATE(X,Y)


                ELSE
                   CALL FATAL_ERROR("INVALID NUMBER OF LAYERS/LEVELS ALLOCATED!")
                END IF

                
             CASE(1)

                DO I = 1,N_NGL
                   
                   IF(N_EFOUND(I))THEN
                      XSPA = N_XC(I)
                      YSPA = N_YC(I)
                      N_VAR%VEC_FLT(I)=INTERP_PZONAL(XSPA,YSPA&
                           &,N2O_EID(I),VAR%VEC_FLT)
                      
                   ELSE
                      
                      N_VAR%VEC_FLT(I)= VAR%VEC_FLT(N2O_EID(I))
                      
                   END IF
                   
                END DO
                
                
             CASE DEFAULT
                CALL FATAL_ERROR("BAD NUMBER OF DIMENSION!")
             END SELECT
             
          case(NF90_DOUBLE)

             IF(DBG_SET(DBG_SBRIO)) THEN
                WRITE(IPT,*) "XYTPE       :: DOUBLE"
                WRITE(IPT,*) "ndims       ::",ndims
             END IF
             SELECT CASE(NDIMS)
             CASE(2)

                MYKB= size(N_VAR%ARR_DBL,2)
                

                IF(MYKB == N_KB )THEN                

                   ALLOCATE(XI(MYKB),YI(MYKB))
                   ALLOCATE(X(KB),Y(KB))


                   DO I = 1,N_NGL
                      k = N2O_EID(I)
                      IF(N_EFOUND(I)) THEN
                         
                         XDP = N_XC(I)
                         YDP = N_YC(I)
                         
                         DO J = 1,MYKB 
                            ZDP = N_Z1(I,J)
                            N_VAR%ARR_DBL(I,J)=INTERP_PZONAL(XDP,YDP&
                                 &,ZDP,KB,k,VAR%ARR_DBL)
                         END DO
                      ELSE
                         
                         X = Z1(k,1:KB)
                         Y = VAR%ARR_DBL(k,1:KB)

                         XI = N_Z1(I,1:MYKB)

                         CALL SINTER_EXTRP_NONE(X,Y,XI,YI,KB,MYKB)

                         N_VAR%ARR_DBL(I,1:MYKB) = YI
                         

                      END IF
                   END DO

                   DEALLOCATE(XI,YI)
                   DEALLOCATE(X,Y)

                ELSEIF( MYKB == N_KBM1 )THEN 

                   ALLOCATE(XI(MYKB),YI(MYKB))
                   ALLOCATE(X(KBM1),Y(KBM1))

                   DO I = 1,N_NGL
                      k = N2O_EID(I)
                      IF(N_EFOUND(I)) THEN

                         XDP = N_XC(I)
                         YDP = N_YC(I)
                         
                         DO J = 1,MYKB 
                            ZDP = N_ZZ1(I,J)
                            N_VAR%ARR_DBL(I,J)=INTERP_PZONAL(XDP,YDP&
                                 &,ZDP,KBM1,k,VAR%ARR_DBL)
                         END DO
                      ELSE

                         X = ZZ1(k,1:KBM1)
                         Y = VAR%ARR_DBL(k,1:KBM1)

                         XI = N_ZZ1(I,1:MYKB)

                         CALL SINTER_EXTRP_NONE(X,Y,XI,YI,KBM1,MYKB)

                         N_VAR%ARR_DBL(I,1:MYKB) = YI

                      END IF
                   END DO

                   DEALLOCATE(XI,YI)
                   DEALLOCATE(X,Y)


                ELSE
                   CALL FATAL_ERROR("INVALID NUMBER OF LAYERS/LEVELS ALLOCATED!")
                END IF

                
             CASE(1)

                DO I = 1,N_NGL
                   
                   IF(N_EFOUND(I))THEN
                      
                      XDP = N_XC(I)
                      YDP = N_YC(I)
                      N_VAR%VEC_DBL(I)=INTERP_PZONAL(XDP,YDP&
                           &,N2O_EID(I),VAR%VEC_DBL)
                      
                   ELSE
                      
                      N_VAR%VEC_DBL(I)= VAR%VEC_DBL(N2O_EID(I))
                      
                   END IF
                   
                END DO
                
                
             CASE DEFAULT
                CALL FATAL_ERROR("BAD NUMBER OF DIMENSION!")
             END SELECT

          END select

          NCFOUT => ADD(NCFOUT, N_VAR)

          CURRENT  => CURRENT%NEXT          
          CYCLE
       END IF


       ! JUST COPY ANYTHING ELSE
       call alloc_var(VAR)
       IF(has_unlimited(VAR)) THEN
          CALL NC_READ_VAR(VAR,STKCNT)
       ELSE
          CALL NC_READ_VAR(VAR)
       END IF
       
       NCFOUT => ADD(NCFOUT, COPY_VAR(VAR))
       CURRENT  => CURRENT%NEXT          


    END DO ! END LOOP ON ALL VARIABLES!
    

    CALL NC_CLOSE(NCFIN)


  END SUBROUTINE CREATE_NEW_FILE




  SUBROUTINE DUMP_NEW_GRID_DATA
    IMPLICIT NONE
    TYPE(NCFILE), POINTER :: NCF
    
    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "Dump New Grid_DATA"
    
    
    NCF => NCFOUT
    
    
    ! WRITE THE STATIC VARIABLES
    CALL NC_WRITE_FILE(NCF)
    
    ! WRITE THE CURRENT STATE VARIABLES
    CALL UPDATE_IODATA(NCF,NOW)
    
    NCF%FTIME%NEXT_STKCNT =1
    CALL NC_WRITE_FILE(NCF)
    
    

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "Dump New Grid_DATA"
    
  END SUBROUTINE DUMP_NEW_GRID_DATA


end module mod_grid

