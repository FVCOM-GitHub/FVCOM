module mod_init
  use all_vars
  use mod_utils
  use mod_input
  use mod_prec
  use mod_ncll
  use mod_nctools
  USE MOD_TIME
  USE MOD_NCDIO
  implicit none

  INTEGER KSL                                 !!NUMBER OF STANDARD SEA LEVELS 
  REAL(SP), ALLOCATABLE :: DPTHSL(:)          !!DEPTH AT STANDARD SEA LEVEL
  REAL(SP), ALLOCATABLE :: TSL(:,:),SSL(:,:)  !!T/S AT STANDARD SEA LEVEL

  Character(Len=120):: FNAME
  Character(Len=120):: OLD_INIT_FILE
  INTEGER, PARAMETER :: INITUNIT = 101


  TYPE(NCFILE), POINTER :: NCF
  TYPE(GRID), SAVE :: MYGRID

  NAMELIST /NML_INIT/                      &
       & INPUT_DIR,                       &
       & OUTPUT_DIR,                      &
       & START_DATE,                      &
       & TIMEZONE,                        &
       & DATE_FORMAT,                     &
       & PROJECTION_REFERENCE,            &
       & GRID_FILE,                       &
       & GRID_FILE_UNITS,                 &
       & OLD_INIT_FILE

  ! DATA VARIABLES

  TYPE(TIME), SAVE :: NOW


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
    WRITE(IPT,*) "!   DATE_FORMAT, (Required)"
    WRITE(IPT,*) "!   GRID_FILE, (Required)"
    WRITE(IPT,*) "!   GRID_FILE_UNITS, (Required)"
    WRITE(IPT,*) "!   OLD_INIT_FILE, (Required)"
    WRITE(IPT,*) "!    "
    WRITE(IPT,*) "!    Example Namelist:"
    write(UNIT=IPT,NML=NML_INIT)


    WRITE(IPT,*) "! NOTES: Do not run this program in parallel!"


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
    READ(UNIT=NMLUNIT, NML=NML_INIT,IOSTAT=ios)
    if(ios .NE. 0 ) then
            CALL FATAL_ERROR("Can Not Read NameList NML_INIT from file: "//trim(FNAME))
    end if
    REWIND(NMLUNIT)

    write(IPT,*) "Read_Name_List:"

    write(UNIT=IPT,NML=NML_INIT)

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



    Pathnfile = trim(INPUT_DIR)//trim(OLD_INIT_FILE)
    Call FOPEN(INITUNIT,trim(pathnfile),'cfr')


  END SUBROUTINE OPEN_FILES


  SUBROUTINE READ_OLD_INIT_TS
    IMPLICIT NONE
    CHARACTER(LEN=80)    :: scan_result
    INTEGER :: ISCAN, I, K, IOS
    REAL(SP), DIMENSION(150) :: TEMP

    IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) "START READING INIT TS FILE"

    ISCAN = SCAN_FILE(INITUNIT,"Standard Levels",FVEC =TEMP  ,NSZE = KSL)
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
    rewind INITUNIT
    DO WHILE(.TRUE.)
       READ(INITUNIT,*,IOSTAT=IOS,END=99) (TSL(1,K), K=1,KSL)
       if (IOS == 0) then
          BackSpace INITUNIT
          exit
       end if

       CYCLE
       
99 Call FATAL_ERROR('Improper formatting of INIT file:',&
           &'Reached end of file with out finding t&s data?')
    END DO

    

    DO I=1,MGL 
       READ(INITUNIT,*) (TSL(I,K), K=1,KSL)
       READ(INITUNIT,*) (SSL(I,K), K=1,KSL)
    END DO

    write(ipt,*) TSL(1,1:5)

    IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) "FINISHED READING INIT TS FILE"


  END SUBROUTINE READ_OLD_INIT_TS
 
  
  SUBROUTINE DUMP_INIT
    IMPLICIT NONE
    TYPE(NCDIM), POINTER :: DIM_KSL
    TYPE(NCVAR), POINTER :: VAR
    TYPE(NCATT), POINTER :: ATT

   IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START DUMP_INIT"

   CALL SET_FVCOM_GRID(MYGRID)

   CALL DEFINE_DIMENSIONS(MYGRID)

   DIM_ksl => NC_MAKE_DIM(name="ksl",len=ksl)

   NCF => NEW_FILE()

   ALLOCATE(NCF%FTIME)
   
   NCF%FNAME=trim(output_dir)//"initfile.nc" 
   
   NCF => ADD(NCF,GRID_FILE_OBJECT(mygrid) )

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
           & DIMTIME=DIM_TIME,&
           TIMEZONE=TIMEZONE)

      NCF  => ADD(NCF,VAR)
   END IF


   VAR => NC_MAKE_AVAR(name='zsl',values=dpthsl,DIM1=DIM_ksl)
   ATT  => NC_MAKE_ATT(name='long_name',values='Standard Depths') 
   VAR  => ADD(VAR,ATT)
   NCF  => ADD(NCF,VAR)

   VAR => NC_MAKE_AVAR(name='tsl',values=tsl, &
        & DIM1=DIM_node,DIM2=DIM_ksl,DIM3=DIM_time)
   ATT  => NC_MAKE_ATT(name='long_name',values='Temperature') 
   VAR  => ADD(VAR,ATT)
   NCF  => ADD(NCF,VAR)

   VAR => NC_MAKE_AVAR(name='ssl',values=ssl,&
        & DIM1=DIM_node,DIM2=DIM_ksl,DIM3=DIM_time)
   ATT  => NC_MAKE_ATT(name='long_name',values='Salinity') 
   VAR  => ADD(VAR,ATT)
   NCF  => ADD(NCF,VAR)
   
   ! WRITE THE STATIC VARIABLES
   CALL NC_WRITE_FILE(NCF)
   
   ! WRITE THE CURRENT STATE VARIABLES
   CALL UPDATE_IODATA(NCF,NOW)

   NCF%FTIME%NEXT_STKCNT =1
   CALL NC_WRITE_FILE(NCF)

   IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END DUMP_RESTART"
   
 END SUBROUTINE DUMP_INIT


end module mod_init

