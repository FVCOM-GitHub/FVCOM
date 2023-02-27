module mod_weather
  use lims
  use mod_utils
  use mod_par
  use control
  use mod_prec
  use mod_ncll
  use mod_nctools
  USE ALL_VARS, only :z,dz,zz,dzz,z1,dz1,zz1,dzz1,vx,vy,xc,yc,NBE,a1u&
       &,a2u,aw0,awx,awy,ntve,nbve,NV
  implicit none


  TYPE(NCDIM), POINTER:: DIMR
  INTEGER, PARAMETER :: nvt = 73   !! times dimension (Kept as Unlimited in file)
  TYPE(NCDIM), POINTER:: DIMNVT
  INTEGER, PARAMETER :: ndv= 10   !! data_variables
  TYPE(NCDIM), POINTER:: DIMDV
  INTEGER, PARAMETER :: nCPL = 11 !! charsPerLevel
  TYPE(NCDIM), POINTER:: DIMCPL
  INTEGER, PARAMETER :: nml = 132!! namelen
  TYPE(NCDIM), POINTER:: DIMNML
  INTEGER, PARAMETER :: nx = 650  !! X dimension
!  INTEGER, PARAMETER :: nx = 65  !! X dimension
  TYPE(NCDIM), POINTER:: DIMX
  INTEGER, PARAMETER :: ny = 550  !! Y dimension
!  INTEGER, PARAMETER :: ny = 55  !! Y dimension
  TYPE(NCDIM), POINTER:: DIMY
  INTEGER, PARAMETER :: nlvs = 1  !! levels_1
  TYPE(NCDIM), POINTER:: DIMLVS
  
  CHARACTER(LEN=120) :: FIN, FOUT,FSET

  TYPE(NCFILE), POINTER ::NCF_IN
  TYPE(NCFILE), POINTER ::NCF_OUT

  real(sp), parameter :: fillval=-99999.0_SP
  real(sp), parameter :: badval=3.402823e+38_SP

  ! SOME DATA POINTERS

  ! REGULAR GRID STUFF
  REAL(SP), PARAMETER :: LLLON=-72.4
  REAL(SP), PARAMETER :: LLLAT=40.2

  REAL(SP), PARAMETER :: URLON=-66.8
  REAL(SP), PARAMETER :: URLAT=45.0

  REAL(SP), PARAMETER :: DLON = (URLON-LLLON)/(real(nx-1))
  REAL(SP), PARAMETER :: DLAT = (URLAT-LLLAT)/(real(ny-1))

  REAL(SP), POINTER :: GRIDLON(:,:), GRIDLAT(:,:)
  REAL(SP), POINTER :: GRIDXCE(:,:), GRIDYCE(:,:)
  REAL(SP), POINTER :: GRIDXTMERC(:,:), GRIDYTMERC(:,:)

  INTEGER, POINTER :: GRIDCELLID(:,:)

!  CHARACTER(LEN=80),PARAMETER :: tmerc="proj=tmerc +datum=NAD83 +lon_0=-70d10 lat_0=42d50 k=.9999666666666667 x_0=900000 y_0=0 no_defs"

  CHARACTER(LEN=80),PARAMETER :: tmerc="init=nad83:1802"
  CHARACTER(LEN=80),PARAMETER :: eqc="proj=eqc"

  ! FVCOM Data
  INTEGER :: MAXELEM, ntimes

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


          else if (opt_sng == "gridlist" .or.opt_sng == "GRIDLIST"&
               & .or.opt_sng == "GridList") then

             call ftn_arg_get(arg_idx,arg_val,FSET) ! [sng] Input file
             FSET=FSET(1:ftn_strlen(FSET))


!!$          else if (opt_sng == "START" .or.opt_sng == "start"&
!!$               & .or.opt_sng == "Start") then
!!$
!!$             call ftn_arg_get(arg_idx,arg_val,Start_Date) ! [sng] Input file
!!$             Start_Date=Start_Date(1:ftn_strlen(Start_Date))
!!$             ! Convert back to a fortran string!
!!$
!!$          else if (opt_sng == "END" .or.opt_sng == "end"&
!!$               & .or.opt_sng == "End") then
!!$
!!$             call ftn_arg_get(arg_idx,arg_val,End_Date) ! [sng] Input file
!!$             End_Date=End_Date(1:ftn_strlen(End_Date))
!!$             ! Convert back to a fortran string!


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

# if defined(DOUBLE_PRECISION)

    CALL FATAL_ERROR("THIS PROGRAM CAN NOT RUN IN DOUBLE PRECISION MODE")

# endif

  END SUBROUTINE GET_COMMANDLINE

  SUBROUTINE MYHELPTXT
    IMPLICIT NONE


    write(IPT,*) "! ARGUMENTS FOR WEATHER_DATA:"
!!$    write(IPT,*) "! START  = (a date or time string)"
!!$    write(IPT,*) "! END    = (a date or time string)"
!!$    write(IPT,*) "! TIMEZONE = (the timezone to use)"
    write(IPT,*) "! INPUT  = The/source/of/data.nc"
    write(IPT,*) "! OUTPUT = The/result/of/extraction.nc"
    write(IPT,*) "! "
    WRITE(IPT,*) "! NOTES: Do not run this program in parallel!"


  END SUBROUTINE MYHELPTXT
!================================================================================
!================================================================================
  SUBROUTINE BUILD_GRID
    IMPLICIT NONE
    INTEGER :: I,J, status
    CHARACTER(LEN=80) :: proj_ref

    ALLOCATE(GRIDLON(nx,ny),stat=status)
    if(status/=0) call fatal_error("Could not allocate: gridlon")

    ALLOCATE(GRIDLAT(nx,ny),stat=status)
    if(status/=0) call fatal_error("Could not allocate: gridlat")

    ALLOCATE(GRIDXCE(nx,ny),stat=status)
    if(status/=0) call fatal_error("Could not allocate: gridXCE")

    ALLOCATE(GRIDYCE(nx,ny),stat=status)
    if(status/=0) call fatal_error("Could not allocate: gridYCE")

    ALLOCATE(GRIDXTMERC(nx,ny),stat=status)
    if(status/=0) call fatal_error("Could not allocate: gridCTMERC")

    ALLOCATE(GRIDYTMERC(nx,ny),stat=status)
    if(status/=0) call fatal_error("Could not allocate: gridYTMERC")

    
    DO I=1,nx
       GRIDLON(I,:)=LLLON + (I-1)*DLON
    END DO

    DO I=1,ny
       GRIDLAT(:,I)=LLLAT + (I-1)*DLAT
    END DO

    CALL DEGREES2METERS(GRIDLON,GRIDLAT,tmerc,GRIDXTMERC,GRIDYTMERC,nx,ny)
    CALL DEGREES2METERS(GRIDLON,GRIDLAT,eqc,GRIDXCE,GRIDYCE,nx,ny)

    IF (DBG_SET(DBG_IO)) THEN
       write(ipt,*) "D{LAT/LON}:",DLAT,DLON
       write(ipt,*) "LLC{LAT/LON}=",LLLAT,LLLON
       write(ipt,*) "URC{LAT/LON}=",URLAT,URLON
       
       write(ipt,*) "GLON({1,nx},1)=",GRIDLON(1,1),GRIDLON(nx,1)
       write(ipt,*) "GLAT(1,{1,ny})=",GRIDLAT(1,1),GRIDLAT(1,ny)

       write(ipt,*) "GXtmerc({1,nx},1)=",GRIDXTMERC(1,1),GRIDXTMERC(nx,1)
       write(ipt,*) "GYtmerc(1,{1,ny})=",GRIDYTMERC(1,1),GRIDYTMERC(1,ny)
       
       write(ipt,*) "GXCE({1,nx},1)=",GRIDXCE(1,1),GRIDXCE(nx,1)
       write(ipt,*) "GYCE(1,{1,ny})=",GRIDYCE(1,1),GRIDYCE(1,ny)
       
       write(ipt,*) "DXCE=",(GRIDXCE(1,1)-GRIDXCE(nx,1))/(nx-1)
       write(ipt,*) "DYCE=",(GRIDYCE(1,1)-GRIDYCE(1,ny))/(ny-1)
       
       write(ipt,*) "DXCE=",GRIDXCE(1,1)-GRIDXCE(2,1)
       write(ipt,*) "DXCE=",GRIDXCE(1,ny)-GRIDXCE(2,ny)
       write(ipt,*) "DYCE=",GRIDYCE(1,1)-GRIDYCE(1,2)
    END IF


  END SUBROUTINE BUILD_GRID
!================================================================================
!================================================================================
  SUBROUTINE BUILD_OUTPUT
    IMPLICIT NONE
    
    TYPE(NCDIM), POINTER:: DIM
    TYPE(NCVAR), POINTER:: VAR
    TYPE(NCATT), POINTER:: ATT
    TYPE(NCFILE),POINTER:: NCF

    CHARACTER(LEN=80),POINTER :: CHR_VEC(:)
    CHARACTER(LEN=80),POINTER :: CHR_SCL
    REAL(SP), POINTER :: FLT_CUB(:,:,:)
    REAL(SP), POINTER :: FLT_ARR(:,:)
    INTEGER, POINTER  :: INT_VEC(:)
    REAL(DP), POINTER :: DBL_SCL
    TYPE(TIME) :: REF_TIME
    TYPE(TIME) :: MODEL_TIME

    INTEGER :: STATUS, I

    REAL(SP), ALLOCATABLE :: att_range(:)

    NCF => NEW_FILE()

    ! DEFINE DIMENSIONS FOR THE FILE
    DIMR   => NC_MAKE_DIM(name='record',len=NF90_UNLIMITED)
    DIMNVT => NC_MAKE_DIM(name='n_valtimes',len=nvt)
    DIMDV  => NC_MAKE_DIM(name='data_variables',len=ndv)
    DIMCPL => NC_MAKE_DIM(name='charsPerLevel',len=nCPL)
    DIMNML => NC_MAKE_DIM(name='namelen',len=nml)
    DIMX   => NC_MAKE_DIM(name='x',len=nx)
    DIMY   => NC_MAKE_DIM(name='y',len=ny)
    DIMLVS => NC_MAKE_DIM(name='levels_1',len=nlvs)    

    ! DEFINE GlOBAL ATTRIBUTES FOR THE FILE
    allocate(att_range(2))

    ATT => NC_MAKE_ATT(name='cdlDate',values="20070914") 
    NCF => ADD(NCF,ATT)

    ATT => NC_MAKE_ATT(name='depictorName',values="fvcomUMD@20048527") 
    NCF => ADD(NCF,ATT)

    ATT => NC_MAKE_ATT(name='projIndex',values=8) 
    NCF => ADD(NCF,ATT)

    ATT => NC_MAKE_ATT(name='projName',values="CYLINDRICAL_EQUIDISTANT") 
    NCF => ADD(NCF,ATT)

    ATT => NC_MAKE_ATT(name='centralLat',values=0.0) 
    NCF => ADD(NCF,ATT)

    ATT => NC_MAKE_ATT(name='centralLon',values=0.0) 
    NCF => ADD(NCF,ATT)

    ATT => NC_MAKE_ATT(name='rotation',values=0.0) 
    NCF => ADD(NCF,ATT)

    ATT => NC_MAKE_ATT(name='xMin',values=-72.4) 
    NCF => ADD(NCF,ATT)

    ATT => NC_MAKE_ATT(name='xMax',values=-66.8) 
    NCF => ADD(NCF,ATT)

    ATT => NC_MAKE_ATT(name='yMax',values=45.0) 
    NCF => ADD(NCF,ATT)

    ATT => NC_MAKE_ATT(name='yMin',values=40.2) 
    NCF => ADD(NCF,ATT)

    ATT => NC_MAKE_ATT(name='lat00',values=40.2) 
    NCF => ADD(NCF,ATT)

    ATT => NC_MAKE_ATT(name='lon00',values=-72.4) 
    NCF => ADD(NCF,ATT)

    ATT => NC_MAKE_ATT(name='latNxNy',values=45.) 
    NCF => ADD(NCF,ATT)

    ATT => NC_MAKE_ATT(name='lonNxNy',values=-66.8) 
    NCF => ADD(NCF,ATT)

    ATT => NC_MAKE_ATT(name='dxKm',values=0.9586438) 
    NCF => ADD(NCF,ATT)

    ATT => NC_MAKE_ATT(name='dyKm',values=0.9713659) 
    NCF => ADD(NCF,ATT)

    ATT => NC_MAKE_ATT(name='latDxDy',values=42.6) 
    NCF => ADD(NCF,ATT)

    ATT => NC_MAKE_ATT(name='lonDxDy',values=-69.60001) 
    NCF => ADD(NCF,ATT)


    ! DEFINE VARIABLES FOR THE FILE

    !! VAR UW
    ALLOCATE(FLT_CUB(nx,ny,nlvs),stat=status)
    if(status/=0) call fatal_error("Could not allocate: UW")
    VAR  => NC_MAKE_PVAR(name='uw',&
         &DIM1=DIMX,&
         &DIM2=DIMY,&
         &DIM3=DIMLVS,&
         &DIM4=DIMR,&
         &values=FLT_CUB)
    NULLIFY(FLT_CUB)

    ATT  => NC_MAKE_ATT(name='long_name',values='u-component wind') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='m/s') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='udunits',values='meter/sec') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='uiname',values='uWind') 
    VAR  => ADD(VAR,ATT)

    att_range=(/-150., 150./)
    ATT  => NC_MAKE_ATT(name='valid_range',values=att_range) 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='_FillValue',values=fillval) 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='_n3D',values=1) 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='levels',values="SFC") 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ALLOCATE(CHR_VEC(nLVS))
    CHR_VEC(1) = "SFC"
    VAR  => NC_MAKE_PVAR(name='uwLevels',&
         &DIM1=DIMCPL,&
         &DIM2=DIMLVS,&
         &values=CHR_VEC)
    NULLIFY(CHR_VEC)
    NCF  => ADD(NCF,VAR)

    ALLOCATE(CHR_VEC(NVT))
    CHR_VEC(:)="1"
    VAR  => NC_MAKE_PVAR(name='uwInventory',&
         &DIM1=DIMLVS,&
         &DIM2=DIMNVT,&
         &values=CHR_VEC)
    NULLIFY(CHR_VEC)
    NCF  => ADD(NCF,VAR)

    !! VAR VW
    ALLOCATE(FLT_CUB(nx,ny,nlvs),stat=status)
    if(status/=0) call fatal_error("Could not allocate: VW")
    VAR  => NC_MAKE_PVAR(name='vw',&
         &DIM1=DIMX,&
         &DIM2=DIMY,&
         &DIM3=DIMLVS,&
         &DIM4=DIMR,&
         &values=FLT_CUB)
    NULLIFY(FLT_CUB)

    ATT  => NC_MAKE_ATT(name='long_name',values='v-component wind') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='m/s') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='udunits',values='meter/sec') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='uiname',values='vWind') 
    VAR  => ADD(VAR,ATT)

    att_range=(/-150., 150./)
    ATT  => NC_MAKE_ATT(name='valid_range',values=att_range) 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='_FillValue',values=fillval) 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='_n3D',values=1) 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='levels',values="SFC") 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ALLOCATE(CHR_VEC(nLVS))
    CHR_VEC(1) = "SFC"
    VAR  => NC_MAKE_PVAR(name='vwLevels',&
         &DIM1=DIMCPL,&
         &DIM2=DIMLVS,&
         &values=CHR_VEC)
    NULLIFY(CHR_VEC)
    NCF  => ADD(NCF,VAR)

    ALLOCATE(CHR_VEC(NVT))
    CHR_VEC(:)="1"
    VAR  => NC_MAKE_PVAR(name='vwInventory',&
         &DIM1=DIMLVS,&
         &DIM2=DIMNVT,&
         &values=CHR_VEC)
    NULLIFY(CHR_VEC)
    NCF  => ADD(NCF,VAR)

    !! VAR UC
    ALLOCATE(FLT_CUB(nx,ny,nlvs),stat=status)
    if(status/=0) call fatal_error("Could not allocate: UC")
    VAR  => NC_MAKE_PVAR(name='uc',&
         &DIM1=DIMX,&
         &DIM2=DIMY,&
         &DIM3=DIMLVS,&
         &DIM4=DIMR,&
         &values=FLT_CUB)
    NULLIFY(FLT_CUB)

    ATT  => NC_MAKE_ATT(name='long_name',values='u-component current') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='m/s') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='udunits',values='meter/sec') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='uiname',values='uCurrent') 
    VAR  => ADD(VAR,ATT)

    att_range=(/-150., 150./)
    ATT  => NC_MAKE_ATT(name='valid_range',values=att_range) 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='_FillValue',values=fillval) 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='_n3D',values=1) 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='levels',values="MSL") 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ALLOCATE(CHR_VEC(nLVS))
    CHR_VEC(1) = "MSL"
    VAR  => NC_MAKE_PVAR(name='ucLevels',&
         &DIM1=DIMCPL,&
         &DIM2=DIMLVS,&
         &values=CHR_VEC)
    NULLIFY(CHR_VEC)
    NCF  => ADD(NCF,VAR)

    ALLOCATE(CHR_VEC(NVT))
    CHR_VEC(:)="1"
    VAR  => NC_MAKE_PVAR(name='ucInventory',&
         &DIM1=DIMLVS,&
         &DIM2=DIMNVT,&
         &values=CHR_VEC)
    NULLIFY(CHR_VEC)
    NCF  => ADD(NCF,VAR)

    !! VAR VC
    ALLOCATE(FLT_CUB(nx,ny,nlvs),stat=status)
    if(status/=0) call fatal_error("Could not allocate: VC")
    VAR  => NC_MAKE_PVAR(name='vc',&
         &DIM1=DIMX,&
         &DIM2=DIMY,&
         &DIM3=DIMLVS,&
         &DIM4=DIMR,&
         &values=FLT_CUB)
    NULLIFY(FLT_CUB)

    ATT  => NC_MAKE_ATT(name='long_name',values='v-component current') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='m/s') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='udunits',values='meter/sec') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='uiname',values='vCurrent') 
    VAR  => ADD(VAR,ATT)

    att_range=(/-150., 150./)
    ATT  => NC_MAKE_ATT(name='valid_range',values=att_range) 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='_FillValue',values=fillval) 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='_n3D',values=1) 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='levels',values="MSL") 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ALLOCATE(CHR_VEC(nLVS))
    CHR_VEC(1) = "MSL"
    VAR  => NC_MAKE_PVAR(name='vcLevels',&
         &DIM1=DIMCPL,&
         &DIM2=DIMLVS,&
         &values=CHR_VEC)
    NULLIFY(CHR_VEC)
    NCF  => ADD(NCF,VAR)

    ALLOCATE(CHR_VEC(NVT))
    CHR_VEC(:)="1"
    VAR  => NC_MAKE_PVAR(name='vcInventory',&
         &DIM1=DIMLVS,&
         &DIM2=DIMNVT,&
         &values=CHR_VEC)
    NULLIFY(CHR_VEC)
    NCF  => ADD(NCF,VAR)

    !! VAR SST
    ALLOCATE(FLT_CUB(nx,ny,nlvs),stat=status)
    if(status/=0) call fatal_error("Could not allocate: SST")
    VAR  => NC_MAKE_PVAR(name='sst',&
         &DIM1=DIMX,&
         &DIM2=DIMY,&
         &DIM3=DIMLVS,&
         &DIM4=DIMR,&
         &values=FLT_CUB)
    NULLIFY(FLT_CUB)

    ATT  => NC_MAKE_ATT(name='long_name',values='sea-surface temperature') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='K') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='udunits',values='degree_Kelvin') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='uiname',values='SST') 
    VAR  => ADD(VAR,ATT)

    att_range=(/-10., 40./)
    ATT  => NC_MAKE_ATT(name='valid_range',values=att_range) 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='_FillValue',values=fillval) 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='_n3D',values=0) 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='levels',values="MSL") 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ALLOCATE(CHR_VEC(nLVS))
    CHR_VEC(1) = "MSL"
    VAR  => NC_MAKE_PVAR(name='sstLevels',&
         &DIM1=DIMCPL,&
         &DIM2=DIMLVS,&
         &values=CHR_VEC)
    NULLIFY(CHR_VEC)
    NCF  => ADD(NCF,VAR)

    ALLOCATE(CHR_VEC(NVT))
    CHR_VEC(:)="1"
    VAR  => NC_MAKE_PVAR(name='sstInventory',&
         &DIM1=DIMLVS,&
         &DIM2=DIMNVT,&
         &values=CHR_VEC)
    NULLIFY(CHR_VEC)
    NCF  => ADD(NCF,VAR)

    !! VAR TEMPERATURE
    ALLOCATE(FLT_CUB(nx,ny,nlvs),stat=status)
    if(status/=0) call fatal_error("Could not allocate: T")
    VAR  => NC_MAKE_PVAR(name='t',&
         &DIM1=DIMX,&
         &DIM2=DIMY,&
         &DIM3=DIMLVS,&
         &DIM4=DIMR,&
         &values=FLT_CUB)
    NULLIFY(FLT_CUB)

    ATT  => NC_MAKE_ATT(name='long_name',values='air temperature') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='K') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='udunits',values='degree_Kelvin') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='uiname',values='T') 
    VAR  => ADD(VAR,ATT)

    att_range=(/-40., 60./)
    ATT  => NC_MAKE_ATT(name='valid_range',values=att_range) 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='_FillValue',values=fillval) 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='_n3D',values=0) 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='levels',values="SFC") 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ALLOCATE(CHR_VEC(nLVS))
    CHR_VEC(1) = "SFC"
    VAR  => NC_MAKE_PVAR(name='tLevels',&
         &DIM1=DIMCPL,&
         &DIM2=DIMLVS,&
         &values=CHR_VEC)
    NULLIFY(CHR_VEC)
    NCF  => ADD(NCF,VAR)

    ALLOCATE(CHR_VEC(NVT))
    CHR_VEC(:)="1"
    VAR  => NC_MAKE_PVAR(name='tInventory',&
         &DIM1=DIMLVS,&
         &DIM2=DIMNVT,&
         &values=CHR_VEC)
    NULLIFY(CHR_VEC)
    NCF  => ADD(NCF,VAR)

    !! VAR icef0
    ALLOCATE(FLT_CUB(nx,ny,nlvs),stat=status)
    if(status/=0) call fatal_error("Could not allocate: ice0")
    VAR  => NC_MAKE_PVAR(name='icef0',&
         &DIM1=DIMX,&
         &DIM2=DIMY,&
         &DIM3=DIMLVS,&
         &DIM4=DIMR,&
         &values=FLT_CUB)
    NULLIFY(FLT_CUB)

    ATT  => NC_MAKE_ATT(name='long_name',values='Icing forecast') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='uiname',values='ICEF') 
    VAR  => ADD(VAR,ATT)

    att_range=(/0., 100./)
    ATT  => NC_MAKE_ATT(name='valid_range',values=att_range) 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='_FillValue',values=fillval) 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='_n3D',values=0) 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='levels',values="MSL") 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ALLOCATE(CHR_VEC(nLVS))
    CHR_VEC(1) = "MSL"
    VAR  => NC_MAKE_PVAR(name='icef0Levels',&
         &DIM1=DIMCPL,&
         &DIM2=DIMLVS,&
         &values=CHR_VEC)
    NULLIFY(CHR_VEC)
    NCF  => ADD(NCF,VAR)

    ALLOCATE(CHR_VEC(NVT))
    CHR_VEC(:)="1"
    VAR  => NC_MAKE_PVAR(name='icef0Inventory',&
         &DIM1=DIMLVS,&
         &DIM2=DIMNVT,&
         &values=CHR_VEC)
    NULLIFY(CHR_VEC)
    NCF  => ADD(NCF,VAR)

    !! VAR icef10
    ALLOCATE(FLT_CUB(nx,ny,nlvs),stat=status)
    if(status/=0) call fatal_error("Could not allocate: ice10")
    VAR  => NC_MAKE_PVAR(name='icef10',&
         &DIM1=DIMX,&
         &DIM2=DIMY,&
         &DIM3=DIMLVS,&
         &DIM4=DIMR,&
         &values=FLT_CUB)
    NULLIFY(FLT_CUB)

    ATT  => NC_MAKE_ATT(name='long_name',values='Icing forecast') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='uiname',values='ICEF') 
    VAR  => ADD(VAR,ATT)

    att_range=(/0., 100./)
    ATT  => NC_MAKE_ATT(name='valid_range',values=att_range) 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='_FillValue',values=fillval) 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='_n3D',values=0) 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='levels',values="MSL") 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ALLOCATE(CHR_VEC(nLVS))
    CHR_VEC(1) = "MSL"
    VAR  => NC_MAKE_PVAR(name='icef10Levels',&
         &DIM1=DIMCPL,&
         &DIM2=DIMLVS,&
         &values=CHR_VEC)
    NULLIFY(CHR_VEC)
    NCF  => ADD(NCF,VAR)

    ALLOCATE(CHR_VEC(NVT))
    CHR_VEC(:)="1"
    VAR  => NC_MAKE_PVAR(name='icef10Inventory',&
         &DIM1=DIMLVS,&
         &DIM2=DIMNVT,&
         &values=CHR_VEC)
    NULLIFY(CHR_VEC)
    NCF  => ADD(NCF,VAR)

    !! VAR SSH
    ALLOCATE(FLT_CUB(nx,ny,nlvs),stat=status)
    if(status/=0) call fatal_error("Could not allocate: ssh")
    VAR  => NC_MAKE_PVAR(name='ssh',&
         &DIM1=DIMX,&
         &DIM2=DIMY,&
         &DIM3=DIMLVS,&
         &DIM4=DIMR,&
         &values=FLT_CUB)
    NULLIFY(FLT_CUB)

    ATT  => NC_MAKE_ATT(name='long_name',values='Sea Surface Height') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='meters') 
    VAR  => ADD(VAR,ATT)

    att_range=(/-20., 20./)
    ATT  => NC_MAKE_ATT(name='valid_range',values=att_range) 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='_FillValue',values=fillval) 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='_n3D',values=0) 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='levels',values="MSL") 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ALLOCATE(CHR_VEC(nLVS))
    CHR_VEC(1) = "MSL"
    VAR  => NC_MAKE_PVAR(name='sshLevels',&
         &DIM1=DIMCPL,&
         &DIM2=DIMLVS,&
         &values=CHR_VEC)
    NULLIFY(CHR_VEC)
    NCF  => ADD(NCF,VAR)

    ALLOCATE(CHR_VEC(NVT))
    CHR_VEC(:)="1"
    VAR  => NC_MAKE_PVAR(name='sshInventory',&
         &DIM1=DIMLVS,&
         &DIM2=DIMNVT,&
         &values=CHR_VEC)
    NULLIFY(CHR_VEC)
    NCF  => ADD(NCF,VAR)

    ! VAR time
    ALLOCATE(INT_VEC(NVT))
    DO I=1,NVT
       INT_VEC(I)= 3600*(I-1)
    END DO
    VAR  => NC_MAKE_PVAR(name='valtimeMINUSreftime',&
         &DIM1=DIMNVT,&
         &values=INT_VEC)
    NULLIFY(INT_VEC)
    
    ATT  => NC_MAKE_ATT(name='units',values="seconds") 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! VAR reftime
    ALLOCATE(DBL_SCL)
    MODEL_TIME = GET_FILE_TIME(NCF_IN,ntimes-72)
    ref_time = READ_DATETIME("1970-1-1 00:00:00.0","YMD","GMT",status)
    DBL_SCL = SECONDS(MODEL_TIME - REF_TIME)

    VAR  => NC_MAKE_PVAR(name='reftime',&
         &values=DBL_SCL)
    NULLIFY(DBL_SCL)
    
    ATT  => NC_MAKE_ATT(name='long_name',values="reference time") 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values="seconds since (1970-1-1 00:00:00.0)") 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)
    
    ! Var origin
    ALLOCATE(CHR_SCL)
    CHR_SCL = "SMAST/ UMASSD"
    VAR  => NC_MAKE_PVAR(name='origin',&
         &DIM1=DIMNML,&
         &values=CHR_SCL)
    NULLIFY(CHR_SCL)
    NCF  => ADD(NCF,VAR)


    ! Var model
    ALLOCATE(CHR_SCL)
    CHR_SCL = "FVCOM"
    VAR  => NC_MAKE_PVAR(name='model',&
         &DIM1=DIMNML,&
         &values=CHR_SCL)
    NULLIFY(CHR_SCL)
    NCF  => ADD(NCF,VAR)
    
    ! VAR staticTopo
    ALLOCATE(FLT_ARR(nx,ny))
    VAR  => NC_MAKE_PVAR(name='staticTopo',&
         &DIM1=DIMX,&
         &DIM2=DIMY,&
         &values=FLT_ARR)
    NULLIFY(FLT_ARR)

    ATT  => NC_MAKE_ATT(name='long_name',values='Topography') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values="meters") 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='_FillValue',values=fillval) 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! VAR staticCoriolis
    ALLOCATE(FLT_ARR(nx,ny))
    VAR  => NC_MAKE_PVAR(name='staticCoriolis',&
         &DIM1=DIMX,&
         &DIM2=DIMY,&
         &values=FLT_ARR)
    NULLIFY(FLT_ARR)

    ATT  => NC_MAKE_ATT(name='long_name',values='Coriolis parameter') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values="/seconds") 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='_FillValue',values=fillval) 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! VAR staticSpacing
    ALLOCATE(FLT_ARR(nx,ny))
    VAR  => NC_MAKE_PVAR(name='staticSpacing',&
         &DIM1=DIMX,&
         &DIM2=DIMY,&
         &values=FLT_ARR)
    NULLIFY(FLT_ARR)

    ATT  => NC_MAKE_ATT(name='long_name',values='Grid spacing') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values="meters") 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='_FillValue',values=fillval) 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    NCF%FTIME => NEW_FTIME()
    

    CALL PRINT_FILE(NCF)

    NCF%FNAME=fout

    NCF_OUT => NCF

  END SUBROUTINE BUILD_OUTPUT
!================================================================================
!================================================================================

  SUBROUTINE READ_INPUT
    USE MOD_INPUT
    USE MOD_SETUP
    IMPLICIT NONE
    TYPE(NCFILE), POINTER :: NCF
    LOGICAL :: FOUND
    TYPE(NCVAR), POINTER :: VAR
    TYPE(NCDIM), POINTER :: DIM


    NCF => NEW_FILE()
    NCF%FNAME=trim(FIN)

    Call NC_OPEN(NCF)
    CALL NC_LOAD(NCF)

    NCF_IN => NCF

    CALL PRINT_FILE(NCF)

    DIM => find_unlimited(NCF_IN,FOUND)
    IF(.not.FOUND) CALL FATAL_ERROR("FILE MUST HAVE UNLIMITED")
    ntimes=dim%dim


    NC_START => NCF
    CALL LOAD_RESTART_GRID(NV)
    m = MGL
    mt = MGL
    n = ngl
    nt = ngl
    
    DIM => FIND_DIM(NCF,"maxelem",FOUND)
    IF(.not. Found) CALL FATAL_ERROR &
         & ("THIS FILE DOES NOT HAVE ALL THE REQUIRED DATA TO DO INTERPOLATION")
    MAXELEM = DIM%DIM

    CALL ALLOCATE_SPACE

    VAR => FIND_VAR(NCF,'x',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'x' IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, vx)
    CALL NC_READ_VAR(VAR)
    CALL NC_DISCONNECT(VAR)

    VAR => FIND_VAR(NCF,'y',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'y' IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, vy)
    CALL NC_READ_VAR(VAR)
    CALL NC_DISCONNECT(VAR)

    VAR => FIND_VAR(NCF,'xc',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'xc' IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, xc)
    CALL NC_READ_VAR(VAR)
    CALL NC_DISCONNECT(VAR)

    VAR => FIND_VAR(NCF,'yc',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'yc' IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, yc)
    CALL NC_READ_VAR(VAR)
    CALL NC_DISCONNECT(VAR)

    VAR => FIND_VAR(NCF,'siglev',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'siglev' IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, z)
    CALL NC_READ_VAR(VAR)
    CALL NC_DISCONNECT(VAR)
    
    CALL n2e3d(z,z1)

    CALL SETUP_SIGMA_DERIVATIVES

    VAR => FIND_VAR(NCF,'a1u',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'a1u' IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, a1u)
    CALL NC_READ_VAR(VAR)
    CALL NC_DISCONNECT(VAR)

    VAR => FIND_VAR(NCF,'a2u',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'a2u' IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, a2u)
    CALL NC_READ_VAR(VAR)
    CALL NC_DISCONNECT(VAR)

    VAR => FIND_VAR(NCF,'nbe',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'nbe' IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, nbe)
    CALL NC_READ_VAR(VAR)
    CALL NC_DISCONNECT(VAR)

    VAR => FIND_VAR(NCF,'ntve',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'ntve' IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, ntve)
    CALL NC_READ_VAR(VAR)
    CALL NC_DISCONNECT(VAR)

    VAR => FIND_VAR(NCF,'nbve',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'nbve' IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, nbve)
    CALL NC_READ_VAR(VAR)
    CALL NC_DISCONNECT(VAR)

    VAR => FIND_VAR(NCF,'awx',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'awx' IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, awx)
    CALL NC_READ_VAR(VAR)
    CALL NC_DISCONNECT(VAR)

    VAR => FIND_VAR(NCF,'awy',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'awy' IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, awy)
    CALL NC_READ_VAR(VAR)
    CALL NC_DISCONNECT(VAR)

    VAR => FIND_VAR(NCF,'aw0',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'aw0' IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, aw0)
    CALL NC_READ_VAR(VAR)
    CALL NC_DISCONNECT(VAR)

    IF (DBG_SET(DBG_IO))THEN
       WRITE(IPT,*) "AW0(1,:)=",AW0(1,:)
       WRITE(IPT,*) "AWx(1,:)=",AWx(1,:)
       WRITE(IPT,*) "AWy(1,:)=",AWy(1,:)

       WRITE(IPT,*) "A1u(10,:)=",A1u(1,:)
       WRITE(IPT,*) "A2u(10,:)=",A2u(1,:)

       
       WRITE(IPT,*) "NBE(1,:)=",NBE(1,:)
       WRITE(IPT,*) "NTVE(1)=",NTVE(1)
       WRITE(IPT,*) "NBVE(1,:)=",NBVE(1,:)

       WRITE(IPT,*) "Z(1,:)=",Z(1,:)
       WRITE(IPT,*) "ZZ(1,:)=",ZZ(1,:)
       WRITE(IPT,*) "Z1(1,:)=",Z1(1,:)
       WRITE(IPT,*) "ZZ1(1,:)=",ZZ1(1,:)

    END IF

  END SUBROUTINE READ_INPUT

  SUBROUTINE ALLOCATE_SPACE
    IMPLICIT NONE
    

    allocate(vx(0:mgl)); vx=0.0_sp
    allocate(vy(0:mgl)); vy=0.0_sp

    allocate(xc(0:ngl)); xc=0.0_sp
    allocate(yc(0:ngl)); yc=0.0_sp
    
    ALLOCATE(Z(0:MGL,KB)); z=0.0_sp
    ALLOCATE(Z1(0:NGL,KB)); z1=0.0_sp
    ALLOCATE(ZZ(0:MGL,KB)); ZZ=0.0_sp
    ALLOCATE(ZZ1(0:NGL,KB)); ZZ1=0.0_sp
    ALLOCATE(DZ(0:MGL,KB)); DZ=0.0_SP
    ALLOCATE(DZ1(0:NGL,KB)); DZ1=0.0_SP
    ALLOCATE(DZZ(0:MGL,KB)); DZZ=0.0_SP
    ALLOCATE(DZZ1(0:NGL,KB)); DZZ1=0.0_SP

    ALLOCATE(NBE(0:ngl,3));
    ALLOCATE(NTVE(0:MGL));
    ALLOCATE(NBVE(0:MGL,MAXELEM));

    ALLOCATE(A1U(0:NGL,4)); A1U   = 0.0_SP
    ALLOCATE(A2U(0:NGL,4)); A2U   = 0.0_SP
    ALLOCATE(AWX(0:NGL,3)); AWX   = 0.0_SP
    ALLOCATE(AWY(0:NGL,3)); AWY   = 0.0_SP
    ALLOCATE(AW0(0:NGL,3)); AW0   = 0.0_SP
    
  END SUBROUTINE ALLOCATE_SPACE



  SUBROUTINE SETUP_CONVERT
    IMPLICIT NONE
    INTEGER :: I,J, status
    integer :: nlast,nthis,nrow
    integer,parameter :: myunit=900
    integer, allocatable :: test(:,:)
    LOGICAL :: FEXIST

    IF(DBG_SET(DBG_LOG)) WRITE(ipt,*) "START SETUP_CONVERT"

    ALLOCATE(GRIDCELLID(nx,ny),stat=status)
    IF(status/=0) CALL FATAL_ERROR("COULD NOT ALLOCATE: GRIDCELLID")

    IF(LEN_TRIM(FSET)==0) THEN
       FEXIST = .FALSE.
    ELSE
       INQUIRE(FILE=FSET,EXIST=FEXIST)
    END IF

    ! GET TRIANGLE DATA!!!
    IF (FEXIST) THEN
       call FOPEN(MYUNIT,FSET,'cfr')
       READ(myunit,*) GRIDCELLID
       close(myunit)
    ELSE

       nrow = 0
       DO I =1,nx
          nlast = nrow
          
          WRITE(IPT,*) "PERCENT FOUND TRIANGLE =",real(i-1)/real(nx)*100.0
          nrow =0
          DO J = 1,ny
             nthis = Find_Element_Containing(GRIDXTMERC(i,j),GRIDYTMERC(i,j),nlast)
             
             GRIDCELLID(i,j) = nthis
             ! SAVE THIS CELL AND TO GUESS ON THE NEXT CELL
             IF(Nthis /=0) nlast = nthis
             ! SAVE FIRST VALUE FROM THIS ROW FOR THE NEXT ROW
             IF(nthis/=0 .and. nrow==0) nrow = nthis
             
          END DO
          
       END DO

       IF(LEN_TRIM(FSET)/=0) THEN
          call FOPEN(MYUNIT,FSET,'ofr')
          
          write(myunit,*) GRIDCELLID
          close(myunit)
       END IF

    END IF


    IF(DBG_SET(DBG_LOG)) WRITE(ipt,*) "END SETUP_CONVERT"

  END SUBROUTINE SETUP_CONVERT

  SUBROUTINE EXTRACT
    IMPLICIT NONE

    TYPE(NCVAR), POINTER :: FICE0,GICE0
    TYPE(NCVAR), POINTER :: FICE10,GICE10
    TYPE(NCVAR), POINTER :: FSAT,GSAT
    TYPE(NCVAR), POINTER :: FSST,GSST
    TYPE(NCVAR), POINTER :: FSSH,GSSH
    TYPE(NCVAR), POINTER :: FU,GU
    TYPE(NCVAR), POINTER :: FV,GV

    TYPE(NCVAR), POINTER :: FWX,GWX
    TYPE(NCVAR), POINTER :: FWY,GWY

    INTEGER, ALLOCATABLE :: STRT(:),CNT(:),STRD(:)
    REAL(SP), PARAMETER :: K2C = 273.15

    TYPE(NCDIM),POINTER :: DIM
    INTEGER :: i,j,k,stack
    LOGICAL :: FOUND

    ! CONNECT TO AND ALLOCATE FVCOM DATA MEMORY
    FICE0=>find_var(NCF_IN,'icing_0kts',FOUND)
    IF(.not.FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIALBE icing_0kts")
    ALLOCATE(FICE0%vec_flt(0:mgl))

    FICE10=>find_var(NCF_IN,'icing_10kts',FOUND)
    IF(.not.FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIALBE icing_10kts")
    ALLOCATE(FICE10%vec_flt(0:mgl))   

    FSAT=>find_var(NCF_IN,'icing_satmp',FOUND)
    IF(.not.FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIALBE icing_satmp")
    ALLOCATE(FSAT%vec_flt(0:mgl))   

    FSST=>find_var(NCF_IN,'temp',FOUND)
    IF(.not.FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIALBE temp")
    ALLOCATE(FSST%vec_flt(0:mgl))  

    FSSH=>find_var(NCF_IN,'zeta',FOUND)
    IF(.not.FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIALBE zeta")
    ALLOCATE(FSSH%vec_flt(0:mgl))  

    FU=>find_var(NCF_IN,'u',FOUND)
    IF(.not.FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIALBE u")
    ALLOCATE(FU%vec_flt(0:ngl))  

    FV=>find_var(NCF_IN,'v',FOUND)
    IF(.not.FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIALBE v")
    ALLOCATE(FV%vec_flt(0:ngl))  

    FWX=>find_var(NCF_IN,'icing_wndx',FOUND)
    IF(.not.FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIALBE icing_wndx")
    ALLOCATE(FWX%vec_flt(0:mgl))  

    FWY=>find_var(NCF_IN,'icing_wndy',FOUND)
    IF(.not.FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIALBE icing_wndy")
    ALLOCATE(FWY%vec_flt(0:mgl))  


    ! CONNECT TO EXISTING MEMORY FOR GRID VARIABLES
    GICE0=>find_var(NCF_OUT,'icef0',FOUND)
    IF(.not.FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIALBE icef0")
    
    GICE10=>find_var(NCF_OUT,'icef10',FOUND)
    IF(.not.FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIALBE icef10")

    GSAT=>find_var(NCF_OUT,'t',FOUND)
    IF(.not.FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIALBE t")

    GSST=>find_var(NCF_OUT,'sst',FOUND)
    IF(.not.FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIALBE sst")

    GSSH=>find_var(NCF_OUT,'ssh',FOUND)
    IF(.not.FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIALBE ssh")

    GU=>find_var(NCF_OUT,'uc',FOUND)
    IF(.not.FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIALBE uc")

    GV=>find_var(NCF_OUT,'vc',FOUND)
    IF(.not.FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIALBE vc")

    GWX=>find_var(NCF_OUT,'uw',FOUND)
    IF(.not.FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIALBE uw")

    GWY=>find_var(NCF_OUT,'vw',FOUND)
    IF(.not.FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIALBE vw")
    

    CALL NC_WRITE_FILE(NCF_OUT)
    
    ALLOCATE(STRT(3),CNT(3),STRD(3))

# if !defined(DOUBLE_PRECISION)
    DO i=1,73

       stack = ntimes-73+i

       call print_real_time(get_file_time(ncf_in,stack),ipt,"stack time is:")

       write(ipt,*) "READING ICE0",i
       CALL NC_READ_VAR(FICE0,stack)

       write(ipt,*) "Interp ICE0",i
       DO j=1,nx
          DO k=1,ny
             IF(gridcellid(j,k) /=0) THEN
                GICE0%cub_flt(j,k,1)=interp_pnodal(gridxtmerc(j,k)&
                     &,gridytmerc(j,k),gridcellid(j,k),FICE0%vec_flt)
             ELSE
                GICE0%cub_flt(j,k,1)=badval
             END IF
          END DO
       END DO


       write(ipt,*) "READING ICE10",i
       CALL NC_READ_VAR(FICE10,stack)

       write(ipt,*) "Interp ICE10",i
       DO j=1,nx
          DO k=1,ny
             IF(gridcellid(j,k) /=0) THEN
                GICE10%cub_flt(j,k,1)=interp_pnodal(gridxtmerc(j,k)&
                     &,gridytmerc(j,k),gridcellid(j,k),FICE10%vec_flt)
             ELSE
                GICE10%cub_flt(j,k,1)=badval
             END IF
          END DO
       END DO

       write(ipt,*) "READING SAT",i
       CALL NC_READ_VAR(FSAT,stack)

       write(ipt,*) "Interp SAT",i
       DO j=1,nx
          DO k=1,ny
             IF(gridcellid(j,k) /=0) THEN
                GSAT%cub_flt(j,k,1)=interp_pnodal(gridxtmerc(j,k)&
                     &,gridytmerc(j,k),gridcellid(j,k),FSAT%vec_flt)+K2C
             ELSE
                GSAT%cub_flt(j,k,1)=badval
             END IF
          END DO
       END DO

       write(ipt,*) "READING SSH",i
       CALL NC_READ_VAR(FSSH,stack)

       write(ipt,*) "Interp SSH",i
       DO j=1,nx
          DO k=1,ny
             IF(gridcellid(j,k) /=0) THEN
                GSSH%cub_flt(j,k,1)=interp_pnodal(gridxtmerc(j,k)&
                     &,gridytmerc(j,k),gridcellid(j,k),FSSH%vec_flt)
             ELSE
                GSSH%cub_flt(j,k,1)=badval
             END IF
          END DO
       END DO

       write(ipt,*) "READING SST",i
       STRT=(/1,1,i/)
       CNT = (/MGL,1,1/)
       STRD=1
       CALL NC_READ_VAR(FSST,IOSTART=STRT,IOCOUNT=CNT,IOSTRIDE=STRD)

       write(ipt,*) "Interp SST",i
       DO j=1,nx
          DO k=1,ny
             IF(gridcellid(j,k) /=0) THEN
                GSST%cub_flt(j,k,1)=interp_pnodal(gridxtmerc(j,k)&
                     &,gridytmerc(j,k),gridcellid(j,k),FSST%vec_flt)+K2C
             ELSE
                GSST%cub_flt(j,k,1)=badval
             END IF
          END DO
       END DO

       write(ipt,*) "READING UC",i
       STRT=(/1,1,i/)
       CNT = (/NGL,1,1/)
       STRD=1
       CALL NC_READ_VAR(FU,IOSTART=STRT,IOCOUNT=CNT,IOSTRIDE=STRD)

       write(ipt,*) "Interp UC",i
       DO j=1,nx
          DO k=1,ny
             IF(gridcellid(j,k) /=0) THEN
                GU%cub_flt(j,k,1)=interp_pzonal(gridxtmerc(j,k)&
                     &,gridytmerc(j,k),gridcellid(j,k),FU%vec_flt)
             ELSE
                GU%cub_flt(j,k,1)=badval
             END IF
          END DO
       END DO

       write(ipt,*) "READING VC",i
       STRT=(/1,1,i/)
       CNT = (/NGL,1,1/)
       STRD=1
       CALL NC_READ_VAR(FV,IOSTART=STRT,IOCOUNT=CNT,IOSTRIDE=STRD)

       write(ipt,*) "Interp VC",i
       DO j=1,nx
          DO k=1,ny
             IF(gridcellid(j,k) /=0) THEN
                GV%cub_flt(j,k,1)=interp_pzonal(gridxtmerc(j,k)&
                     &,gridytmerc(j,k),gridcellid(j,k),FV%vec_flt)
             ELSE
                GV%cub_flt(j,k,1)=badval
             END IF
          END DO
       END DO

       write(ipt,*) "READING wndx",i
       CALL NC_READ_VAR(FWX,stack)

       write(ipt,*) "Interp wndx",i
       DO j=1,nx
          DO k=1,ny
             IF(gridcellid(j,k) /=0) THEN
                GWX%cub_flt(j,k,1)=interp_pnodal(gridxtmerc(j,k)&
                     &,gridytmerc(j,k),gridcellid(j,k),FWX%vec_flt)
             ELSE
                GWX%cub_flt(j,k,1)=badval
             END IF
          END DO
       END DO

       write(ipt,*) "READING wndy",i
       CALL NC_READ_VAR(FWY,stack)

       write(ipt,*) "Interp wndy",i
       DO j=1,nx
          DO k=1,ny
             IF(gridcellid(j,k) /=0) THEN
                GWY%cub_flt(j,k,1)=interp_pnodal(gridxtmerc(j,k)&
                     &,gridytmerc(j,k),gridcellid(j,k),FWY%vec_flt)
             ELSE
                GWY%cub_flt(j,k,1)=badval
             END IF
          END DO
       END DO


       write(ipt,*) "WRITE STACK",i
       NCF_OUT%FTIME%NEXT_STKCNT = i
       CALL NC_WRITE_FILE(NCF_OUT)
    END DO

# endif


  END SUBROUTINE EXTRACT

  
end module mod_weather
