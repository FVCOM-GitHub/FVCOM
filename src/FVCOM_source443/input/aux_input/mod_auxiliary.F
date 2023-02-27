module mod_auxiliary
  use all_vars
  use mod_utils
  use mod_par
  use mod_input
  use mod_prec
  use mod_ncll
  use mod_nctools
  USE MOD_TIME
  USE MOD_NCDIO
  implicit none


  ! SET THE HORIZONTAL MIXING COEFFICIENT
  REAL(SP), PARAMETER :: HMC = 40.00_SP


  Character(Len=120):: FNAME
!  Character(Len=120):: MYINPUT
!  INTEGER, PARAMETER :: MYUNIT = 101


  TYPE(NCFILE), POINTER :: NCF
  TYPE(GRID), SAVE :: MYGRID

  NAMELIST /NML_AUX/                      &
       & INPUT_DIR,                       &
       & OUTPUT_DIR,                      &
       & PROJECTION_REFERENCE,            &
       & GRID_FILE_UNITS,                 &
       & GRID_FILE,                       &
       & SIGMA_LEVELS_FILE,               &
       & DEPTH_FILE
!      & MYINPUT                

  ! DATA VARIABLES

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
    WRITE(IPT,*) "!   GRID_FILE, (Required)"
    WRITE(IPT,*) "!   SIGMA_LEVELS_FILE, (Required)"
    WRITE(IPT,*) "!   DEPTH_FILE, (Required)"
!    WRITE(IPT,*) "!   MYINPUT"
    WRITE(IPT,*) "!   Example namelist:"
    write(UNIT=IPT,NML=NML_AUX)

    WRITE(IPT,*) "! NOTES: Do not run this program in parallel!"


  END SUBROUTINE MYHELPTXT

  SUBROUTINE INITIALIZE_NAMELIST
    IMPLICIT NONE

    TIMEZONE = 'none'

    ! INITIALIZE SOME FIELDS
    INPUT_DIR = "NONE"
    OUTPUT_DIR = "NONE"
    GRID_FILE = "NONE"
    SIGMA_LEVELS_FILE = "NONE"
    DEPTH_FILE = "NONE"

  END SUBROUTINE INITIALIZE_NAMELIST


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
    READ(UNIT=NMLUNIT, NML=NML_AUX,IOSTAT=ios)
    if(ios .NE. 0 ) then
       write(UNIT=IPT,NML=NML_AUX)
       CALL FATAL_ERROR("Can Not Read NameList NML_AUX from file: "//trim(FNAME))
    end if
    REWIND(NMLUNIT)

    write(IPT,*) "Read_Name_List:"

    write(UNIT=IPT,NML=NML_AUX)

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

!    Pathnfile = trim(INPUT_DIR)//trim(MYINPUT)
!    Call FOPEN(MYUNIT,trim(pathnfile),'cur')


  END SUBROUTINE OPEN_FILES

  SUBROUTINE ALLOCATE_SPACE
    IMPLICIT NONE
    INTEGER :: STATUS


    ! VARIABLE LOADED FROM COLD START FILES
    allocate(H(0:mgl)); H=0.0_sp
    allocate(H1(0:ngl)); H1=0.0_sp
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

    ! VARIABLES FOR HVC AND BOTTOM ROUGHNESS
    ALLOCATE(CC_HVC(0:NGL)); CC_HVC=0.0_SP
    ALLOCATE(NN_HVC(0:MGL)); NN_HVC=0.0_SP

    ALLOCATE(CC_Z0B(0:NGL)); CC_Z0B=0.0_SP


  END SUBROUTINE ALLOCATE_SPACE

!====================================================================
! CREATE YOUR SPECIALIZED INPUT DATA HERE AND SAVE TO A FILE!

  ! DATA FOR HORIZONTAL VISCOSITY
  SUBROUTINE CREATE_HVC(NN,CC)
    IMPLICIT NONE
    REAL(SP),ALLOCATABLE :: NN(:),CC(:)
    INTEGER :: I


    WRITE(IPT,*) "! SETTING HORIZONTAL MIXING COEFFICITENT FOR GOM3 MODEL"
    WRITE(IPT,*) "! "

    CC = HMC
    CC(0) = 0.0_SP
    NN = HMC 
    NN(0) = 0.0_SP

    DO I=1,NGL
        IF( H1(I).GE.1499.0 ) CC(I) = 2.0_SP * HMC 
     ENDDO
     DO I=1,MGL
        IF(H(I).GE.1499.0) NN(I) = 5.0_SP * HMC
     ENDDO  

  END SUBROUTINE CREATE_HVC
  ! DATA FOR BOTTOM ROUGHNESS 
  SUBROUTINE CREATE_BOTTOM_ROUGHNESS(Z0)
    IMPLICIT NONE
    REAL(SP),ALLOCATABLE :: Z0(:)
    INTEGER :: I

    WRITE(IPT,*) "! SETTING BOTTOM ROUGHNESS LENGTH SCALE FOR GOM3"
    WRITE(IPT,*) "! "

    ! USE WHERE INSTEAD OF DO,IF
    DO I=1, NGL
       IF(H1(I) <= 40.0_SP) THEN
          Z0(I)=3.E-3_SP
       ELSE IF(H1(I) > 40.0_SP.AND.H1(I) <= 70.0_SP) THEN
          Z0(I)=3.E-3_SP*EXP(-(H1(I)-40.0_SP)/8.8204_SP)
       ELSE IF(H1(I) > 70.0_SP.AND.H1(I) <= 100.0_SP) THEN
          Z0(I)=1.E-4_SP*EXP(-(H1(I)-70.0_SP)/13.0288_SP)
       ELSE
          Z0(I)=1.E-5_SP
       END IF
    END DO

  END SUBROUTINE CREATE_BOTTOM_ROUGHNESS
  
 FUNCTION AUX_FILE_OBJECT() RESULT(NCF)
   USE MOD_CLOCK
   IMPLICIT NONE

    INTEGER :: status
    LOGICAL, SAVE :: IOPROC_ALLOCATED = .FALSE.
    TYPE(NCFILE), POINTER :: NCF
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCATT),  POINTER :: ATT
    character(len=100)    :: timestamp,  temp, netcdf_convention

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START AUX_FILE_OBJECT"

  ! ALLOCATE THE NEW FILE OBJECT
    NCF => NEW_FILE()

! Bottom Roughness Lengthscale
    VAR  => NC_MAKE_AVAR(name='z0b',values=CC_Z0B, DIM1= DIM_nele)

    ATT  => NC_MAKE_ATT(name='long_name',values='Bottom Roughness Lengthscale') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='m') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

! Nodal Horizontal Mixing Coefficient
    VAR  => NC_MAKE_AVAR(name='nn_hvc',values=NN_HVC, DIM1= DIM_node)

    ATT  => NC_MAKE_ATT(name='long_name',values='NN_HVC') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='m+2 s-1') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

! Element Horizontal Mixing Coefficient
    VAR  => NC_MAKE_AVAR(name='cc_hvc',values=CC_HVC, DIM1= DIM_nele)

    ATT  => NC_MAKE_ATT(name='long_name',values='CC_HVC') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='m+2 s-1') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)


    ! DIMS SHOULD NOW BE DEALLOCATED BUT I GUESS IT DOESN'T REALLY MATTER.

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END AUX_FILE_OBJECT"

  END FUNCTION AUX_FILE_OBJECT

 SUBROUTINE DUMP_AUX
   IMPLICIT NONE

   IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START DUMP_RESTART"

   CALL SET_FVCOM_GRID(MYGRID)
   
   CALL DEFINE_DIMENSIONS(MYGRID)
   
   NCF => NEW_FILE()

   NCF%FNAME=trim(output_dir)//"auxiliary_data.nc" 
   
   NCF => ADD(NCF,AUX_FILE_OBJECT() )

!   NCF => ADD(NCF,AUX_FILE_OBJECT() )

   ! WRITE THE STATIC VARIABLES
   CALL NC_WRITE_FILE(NCF)
   
   call kill_dimensions

   IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END DUMP_RESTART"
   
 END SUBROUTINE DUMP_AUX


end module mod_auxiliary

