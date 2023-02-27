PROGRAM Init_create
  use mod_input
  use mod_init
  use mod_setup
  USE MOD_SET_TIME
  IMPLICIT NONE
  INTEGER :: STATUS, I
  CHARACTER(LEN=4) :: BFLAG
  character(len=*),parameter::CVS_Id="$Id$" ! [sng] CVS Identification
  character(len=*),parameter::CVS_Date="$Date$" ! [sng] Date string
  character(len=*),parameter::CVS_Name="$Name$" ! [sng] File name string
  character(len=*),parameter::CVS_Revision="$Revision$" ! [sng] File revision string



  ! INTIALIZE MEMORY FROM libfvcom.a
  CALL INITIALIZE_CONTROL("Init TS Create")

# if defined (MULTIPROCESSOR)
  CALL INIT_MPI_ENV(MYID,NPROCS,SERIAL,PAR,MSR,MSRID)
  MPI_FVCOM_GROUP = MPI_COMM_WORLD ! FOR NOW MAKE THEM EQUAL
# endif

  IF (PAR) CALL FATAL_ERROR &
       & ("THIS PROGRAM IS WAY TO SHORT TO BOTHER WITH WRITING AS A PARALLE CODE.",&
       & "STOP AND RERUN WITH JUST ONE PROCESSOR!")


  CALL GET_COMMANDLINE(CVS_ID,CVS_Date,CVS_Name,CVS_Revision)

  IF (DBG_SET(DBG_LOG)) THEN
     WRITE(IPT,*) "! ========================="
     WRITE(IPT,*) "!   BEGIN INIT CREATE"
     WRITE(IPT,*) "! ========================="
  END IF

  CALL NAME_LIST_INITIALIZE

  CALL READ_NAMELIST

  CALL CNTRL_PRMTRS

  IF(USE_REAL_WORLD_TIME)THEN

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


  CALL PRINT_REAL_TIME(NOW,IPT,"Init TS Time Is:")

  CALL CHECK_IO_DIRS

  CALL OPEN_FILES

  CALL READ_COLDSTART_GRID(GRIDUNIT,MGL,NGL,NV)
  m = MGL
  mt = MGL
  n = ngl
  nt = ngl

  ALLOCATE(X_LCL(0:MGL),Y_LCL(0:MGL))
  CALL READ_COLDSTART_COORDS(GRIDUNIT,MGL,X_LCL,Y_LCL)
  CLOSE(GRIDUNIT)

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

!  CALL COORDINATE_UNITS(X_LCL,Y_LCL)
!  CALL SETUP_CENTER_COORDS

  DEALLOCATE(X_LCL,Y_LCL)
  
  CALL READ_OLD_INIT_TS

  CALL DUMP_INIT

  IF (DBG_SET(DBG_LOG)) WRITE(IPT,*) "! TADA!!!"
  CALL PSHUTDOWN


END PROGRAM Init_create
