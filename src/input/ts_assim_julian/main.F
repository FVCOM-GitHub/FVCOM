!!!  the is program can convert the grided standard level TS data to model grid for data assimilation
!!!  file  ggao 032009

PROGRAM TS_clim
  use mod_input
  use mod_ts_clim
  use mod_setup
  USE MOD_NCDIO
  USE EQS_OF_STATE

  IMPLICIT NONE
  INTEGER :: STATUS, I
  CHARACTER(LEN=4) :: BFLAG
  character(len=*),parameter::CVS_Id="$Id$" ! [sng] CVS Identification
  character(len=*),parameter::CVS_Date="$Date$" ! [sng] Date string
  character(len=*),parameter::CVS_Name="$Name$" ! [sng] File name string
  character(len=*),parameter::CVS_Revision="$Revision$" ! [sng] File revision string

  CHARACTER(len=120) :: TSC_OUT_NAME
  TYPE(TIME) :: NOW
  TYPE(NCFTIME), POINTER ::FTM


  ! INTIALIZE MEMORY FROM libfvcom.a
  CALL INITIALIZE_CONTROL("Climatology TS Create")

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
     WRITE(IPT,*) "!  CLIMATOLOGY TS CREATE"
     WRITE(IPT,*) "! ========================="
  END IF

  CALL NAME_LIST_INITIALIZE

  CALL READ_NAMELIST

  CALL CNTRL_PRMTRS

! CALL PRINT_REAL_TIME(NOW,IPT,"Init TS Time Is:")
!!  ggao 
  CALL CHECK_IO_DIRS

  CALL OPEN_FILES

!  CALL GET_FVCOM_GRID
!   if(.false.) then
 CALL READ_COLDSTART_GRID(GRIDUNIT,MGL,NGL,NV)
  m = MGL
  mt = MGL
  n = ngl
  nt = ngl

  CALL READ_COLDSTART_SIGMA
  CLOSE(SIGMAUNIT)
  KBM1 = KB - 1
  KBM2 = KB - 2

  CALL ALLOCATE_SPACE

  ALLOCATE(X_LCL(0:MGL),Y_LCL(0:MGL))
  CALL READ_COLDSTART_COORDS(GRIDUNIT,MGL,X_LCL,Y_LCL)
  CLOSE(GRIDUNIT)

  CALL COORDINATE_UNITS(X_LCL,Y_LCL)
  CALL SETUP_CENTER_COORDS

  CALL READ_COLDSTART_DEPTH(DEPTHUNIT,MGL,X_LCL,Y_LCL,H)
  D = H
  CLOSE(DEPTHUNIT)


  CALL Setup_Sigma
  CALL SETUP_SIGMA_DERIVATIVES
  DEALLOCATE(X_LCL,Y_LCL)
!   endif

!  CALL SET_PARAMETERS

  CALL SET_TIME
!!   read the TS-CLIM data information
   CALL READ_OLD_TS_INF

   !CALL LOAD_TSC(OLD_TS_FILE)
!   CALL LOAD_TSC

!!   set up the output file
   TSC_OUT_NAME = TRIM(output_dir)//trim(output_file)
   CALL MAKE_TS_OUTFILE(TSC_OUT_NAME)
   FTM => NC_OUT%FTIME
   

  now = StartTime
  IINT = 1

!  CALL READ_TSC(TC0=TC,SC0=SC)
!  Now = TC%dtm

  DO WHILE(IINT <= 4)

      CALL READ_TSC(now,TC,SC)

        print*,t1(205,1:2)
        print*,tsl(205,1:2)
!!

!  IF(CONVECTIVE_OVERTURNING)THEN
     CALL DENS2
     CALL CONV_OVER
     CALL DENS2
     CALL CONV_OVER

!     IF(.NOT.BAROTROPIC)THEN
!        SELECT CASE(SEA_WATER_DENSITY_FUNCTION)
!        CASE(SW_DENS1)
!           CALL DENS1
!        CASE(SW_DENS2)
!           CALL DENS2
!        CASE(SW_DENS3)
!           CALL DENS3
!        CASE DEFAULT
!           CALL FATAL_ERROR("INVALID DENSITY FUNCTION SELECTED:",&
!                & "   "//TRIM(SEA_WATER_DENSITY_FUNCTION) )
!        END SELECT
!
!     END IF
!  END IF
!!
        FTM%NEXT_STKCNT = FTM%NEXT_STKCNT + 1
        CALL UPDATE_IODATA(NC_OUT,NOW)
        CALL NC_WRITE_FILE(NC_OUT)

!     NOW = NOW + INTERVAL
     IINT = IINT + 1

  END DO


  IF (DBG_SET(DBG_LOG)) WRITE(IPT,*) "! TADA!!!"
  CALL PSHUTDOWN

END PROGRAM TS_clim
