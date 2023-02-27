MODULE MOD_SST2GRD
  USE MOD_NCTOOLS
  USE MOD_NCDIO
  USE MOD_UTILS
  USE MOD_INPUT
  USE MOD_TIME
  USE CONTROL
  IMPLICIT NONE



  REAL(SP) :: CVAL
  LOGICAL :: CONSTANT

  INTEGER, PARAMETER       :: NX = 600
  INTEGER, PARAMETER       :: NY = 320
  INTEGER, PARAMETER       :: FUNIT = 46

  ! MICRO_SECONDS PER DAY
  TYPE(TIME), SAVE :: SST_INTERVAL 


  REAL(SP), ALLOCATABLE, DIMENSION(:,:) :: SSTin,Xs,Ys, lons, lats
  REAL(SP), ALLOCATABLE, DIMENSION(:) :: SSTout
  INTEGER,  ALLOCATABLE, DIMENSION(:,:) :: MASK

  INTEGER :: IDAY


  CHARACTER CH*3,month(12)*3,yn*1,cyear*4


  CHARACTER(len=80) :: Lat_data_name
  CHARACTER(len=80) :: Lon_data_name
  CHARACTER(len=80) :: Data_prefix

  NAMELIST /NML_SST/         &
       & INPUT_DIR,            &
       & LAT_DATA_NAME,        &
       & LON_DATA_NAME,        &
       & DATA_PREFIX,          &
       & OUTPUT_DIR,           & 
       & GRID_FILE,            &
       & GRID_FILE_UNITS,      &
       & PROJECTION_REFERENCE, &
       & START_DATE,           &
       & TIMEZONE           

  CHARACTER(len=120) :: FNAME
  LOGICAL           :: FEXIST


  TYPE(INTERP_WEIGHTS) :: WEIGHTS

  TYPE(NCFILE), POINTER :: NC_OUT
  TYPE(NCFTIME), POINTER :: FTM

  TYPE(GRID), SAVE :: MYGRID

  TYPE(TIME), SAVE :: NOW


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

    CONSTANT = .false.
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

          else if (opt_sng == "dbg_par" .or.opt_sng == "Dbg_Par"&
               & .or.opt_sng == "DBG_PAR") then

             dbg_par = .true.


          else if (opt_sng == "constant" .or. opt_sng == "Constant" ) then
             call ftn_arg_get(arg_idx,arg_val,cval) ! [enm] Debugging level

             write(ipt,*) "Seeting a constant value:", cval
             CONSTANT = .TRUE.

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
    WRITE(IPT,*) "! --constant=X.X (Set constant value)"
    WRITE(IPT,*) "!   "
    WRITE(IPT,*) "!   "
    WRITE(IPT,*) "!   The namelist runfile for the program! "
    WRITE(IPT,*) "!   "
    WRITE(IPT,*) "!   Namelist OPTIONS: "
    WRITE(IPT,*) "!    INPUT_DIR"
    WRITE(IPT,*) "!    LAT_DATA_NAME"
    WRITE(IPT,*) "!    LON_DATA_NAME"
    WRITE(IPT,*) "!    DATA_PREFIX"
    WRITE(IPT,*) "!    OUTPUT_DIR"
    WRITE(IPT,*) "!    GRID_DIR"
    WRITE(IPT,*) "!    GRID_FILE"
    WRITE(IPT,*) "!    PROJECTION_REFERENCE"
    WRITE(IPT,*) "!    START_DATE"
    WRITE(IPT,*) "!    TIMEZONE"
    WRITE(IPT,*) "!    "
    WRITE(IPT,*) "!    "
    WRITE(IPT,*) "!    "
    WRITE(IPT,*) "!    Exmaple Namelist"
    write(UNIT=IPT,NML=NML_SST)


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
    READ(UNIT=NMLUNIT, NML=NML_SST,IOSTAT=ios)
   if(ios .NE. 0 ) THEN
       if(DBG_SET(dbg_log)) write(UNIT=IPT,NML=NML_SST)
       CALL FATAL_ERROR("Can Not Read NameList NML_INPUT from file: "//trim(FNAME))
    end if

    REWIND(NMLUNIT)

    if(DBG_SET(dbg_scl)) &
         & write(IPT,*) "Read_Name_List:"

    if(DBG_SET(dbg_scl)) &
         & write(UNIT=IPT,NML=NML_SST)

    CLOSE(NMLUNIT)


  END SUBROUTINE READ_NAMELIST


  SUBROUTINE SET_TIME
    USE mod_set_time

    IMPLICIT NONE
    character(len=4) :: bflag
    INTEGER :: STATUS

    ! assume sst data interval is one day!
    SST_INTERVAL = days2time(1.0_DP)
    
    if(USE_REAL_WORLD_TIME) then
       
       NOW = READ_DATETIME(Start_Date,"YMD",TIMEZONE,status)
       IF(status /= 1) call fatal_error &
            &("could not parse time_origin or time_zone passed for spectral tidal forcing file?")
       
       CALL PRINT_REAL_TIME(NOW,IPT,'START DATE')
       
    else
       
       CALL IDEAL_TIME_STRING2TIME(Start_date,BFLAG,NOW,IINT)
       IF(BFLAG == 'step') CALL FATAL_ERROR&
            &("You must Secify a time, not a step, for this restart file", &
            & "The Step will be set by the old restart file...")
       
       CALL PRINT_TIME(NOW,IPT,'START DATE')
       
    END if
    
  END SUBROUTINE SET_TIME

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


  END SUBROUTINE GET_FVCOM_GRID

  SUBROUTINE RUN_TGE
    USE MOD_OBCS, only : iobcn
    USE ALL_VARS
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

  SUBROUTINE CREATE_INTERP
    USE ALL_VARS, only : XM,YM
    IMPLICIT NONE

    INTEGER :: j,i,source,IERR


    IF(MSR) THEN
       ! open the grid file for the sst data longitude
       CALL FOPEN(FUNIT,TRIM(INPUT_DIR)//TRIM(LON_DATA_NAME),'cfr')
       DO J=1,NY
          READ(FUNIT,*) (LONS(I,J), I=1,NX)
       ENDDO
       CLOSE(FUNIT)
       
       ! open the grid file for the sst data latitude
       CALL FOPEN(FUNIT,TRIM(INPUT_DIR)//TRIM(LAT_DATA_NAME),'cfr')
       
       DO J=1,NY
          READ(FUNIT,*) (LATS(I,J), I=1,NX)
       ENDDO
       CLOSE(FUNIT)
       
       
       write(ipt,*) "READ SST DATA lon and lat: "
       WRITE(ipt,*) "min/max(lon)",minval(lons),maxval(lons)
       WRITE(ipt,*) "min/max(lat)",minval(lats),maxval(lats)
       
       CALL DEGREES2METERS(LONS,LATS,PROJECTION_REFERENCE,XS,YS,nx,ny)
       
       write(ipt,*) "CONVERTED TO METERS: "  
       WRITE(ipt,*) "min/max(Xs)",minval(Xs),maxval(Xs)
       WRITE(ipt,*) "min/max(Ys)",minval(Ys),maxval(Ys)
       
       
       !READ FIRST SST MAP TO GET MASK OF INPUT DATA  
       write(CH,'(i3.3)') 1  ! 
       
       FNAME = TRIM(INPUT_DIR)//TRIM(data_prefix)//CH
       CALL FOPEN(FUNIT,TRIM(FNAME),'cfr')


       DO J=1,NY
          READ(FUNIT,*) (SSTin(I,J),I=1,NX)   
       ENDDO
       CLOSE(FUNIT)

       
       ! BUILD LIST OF VALID DATA LOCATIONS
       mask = 1
       WHERE(SSTin .GT. -90) mask = 0

    END IF

    IF(PAR) THEN

# if defined(MULTIPROCESSOR)
       
       IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "SENDING COORDS DATA"
       
       SOURCE = MSRID -1

       CALL MPI_BCAST(xs,nx*ny,MPI_F,SOURCE,MPI_FVCOM_GROUP,ierr)

       CALL MPI_BCAST(ys,nx*ny,MPI_F,SOURCE,MPI_FVCOM_GROUP,ierr)


       CALL MPI_BCAST(lons,nx*ny,MPI_F,SOURCE,MPI_FVCOM_GROUP,ierr)
       CALL MPI_BCAST(lats,nx*ny,MPI_F,SOURCE,MPI_FVCOM_GROUP,ierr)

       CALL MPI_BCAST(mask,nx*ny,MPI_INTEGER,SOURCE,MPI_FVCOM_GROUP,ierr)
     
# endif


    END IF

    CALL RUN_TGE

    CALL SETUP_INTERP_BILINEAR_A(XS,YS,XM,YM,WEIGHTS,MASK)


  END SUBROUTINE CREATE_INTERP


  SUBROUTINE ALLOCATE_SPACE
    IMPLICIT NONE

    
    ALLOCATE(SSTin(NX,NY))
    ALLOCATE(XS(NX,NY))
    ALLOCATE(YS(NX,NY))
    ALLOCATE(LATS(NX,NY))
    ALLOCATE(LONS(NX,NY))
    ALLOCATE(MASK(NX,NY))



  END SUBROUTINE ALLOCATE_SPACE


  SUBROUTINE MY_OUTFILE
    USE ALL_VARS
    IMPLICIT NONE
    TYPE(NCFILE), POINTER :: NCF
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCATT),  POINTER :: ATT

    TYPE(NCDIM),  POINTER :: DIM_node
    TYPE(NCDIM),  POINTER :: DIM_nele
    TYPE(NCDIM),  POINTER :: DIM_three
    TYPE(NCDIM),  POINTER :: DIM_DateStrLen
    TYPE(NCDIM),  POINTER :: DIM_time


    CALL SET_FVCOM_GRID(MYGRID)
    CALL DEFINE_DIMENSIONS(MYGRID)

    ! ALLOCATE THE NEW FILE OBJECT
    NCF => NEW_FILE()

    NC_OUT => NCF

    ALLOCATE(NCF%FTIME)


    NCF%FNAME = TRIM(OUTPUT_DIR)//'sst.nc'

    NCF => ADD(NCF, GRID_FILE_OBJECT(MYGRID) ) 

    NCF => ADD(NCF, TIME_FILE_OBJECT() )

    ALLOCATE(SSTOUT(0:MT))


    ! SST
    VAR  => NC_MAKE_AVAR(name='sst',&
         & values=SSTout, DIM1= DIM_node, DIM2= DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Sea Surface Temperature') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='celcius') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)


  END SUBROUTINE MY_OUTFILE


  SUBROUTINE UPDATE_SST
    IMPLICIT NONE
    INTEGER :: I,J, VALUE

     
    
    ! input data
    write(CH,'(i3.3)') iday  !
    FNAME = TRIM(INPUT_DIR)//TRIM(DATA_PREFIX)//CH
    INQUIRE(FILE=FNAME,EXIST=FEXIST)
    IF(.NOT. FEXIST) THEN
       CALL WARNING("COULD NOT FIND FILE: "//TRIM(FNAME),&
            & "Incriment Year!:")
       
       iday=1
       write(CH,'(i3.3)') iday
       read(DATA_PREFIX(4:7),*) value
       value = value + 1
       write(ipt,*) "Year=",value
       write(DATA_PREFIX(4:7),'(I4.4)') value
       
       
       FNAME = TRIM(INPUT_DIR)//TRIM(DATA_PREFIX)//CH
       INQUIRE(FILE=FNAME,EXIST=FEXIST)
       IF(.NOT. FEXIST)  THEN
          WRITE(IPT,*) "COULD NOT FIND ANY MORE DATA!"
          Write(IPT,*) "Last File Name: ",TRIM(FNAME)
          CALL PSHUTDOWN
       END IF
       
    END IF
    
    CALL FOPEN(FUNIT,TRIM(FNAME),'cfr')
    
    
    DO J=1,NY
       READ(FUNIT,*) (SSTin(I,J),I=1,NX)   
    ENDDO
    CLOSE(FUNIT)
    
    
    CALL INTERP_BILINEAR_A(SSTin,WEIGHTS,SSTout)
    


  END SUBROUTINE UPDATE_SST


END MODULE MOD_SST2GRD
