PROGRAM obc_create
  use mod_input
  use mod_obcreate
  use mod_setup
  IMPLICIT NONE

  character(len=*),parameter::CVS_Id="$Id$" ! [sng] CVS Identification
  character(len=*),parameter::CVS_Date="$Date$" ! [sng] Date string
  character(len=*),parameter::CVS_Name="$Name$" ! [sng] File name string
  character(len=*),parameter::CVS_Revision="$Revision$" ! [sng] File revision string



  ! INTIALIZE MEMORY FROM libfvcom.a
  CALL INITIALIZE_CONTROL("obcreate")

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
     WRITE(IPT,*) "!   BEGIN OBC CREATE"
     WRITE(IPT,*) "! ========================="
  END IF

  CALL NAME_LIST_INITIALIZE
  CALL READ_NAMELIST


  !Set some variables we don't really need but are check in
  ! CNTRL_PRMTRS
  VERTICAL_MIXING_COEFFICIENT=1.0_SP
  VERTICAL_PRANDTL_NUMBER=1.0_SP
  HORIZONTAL_MIXING_COEFFICIENT=1.0_SP
  HORIZONTAL_PRANDTL_NUMBER=1.0_SP
  BOTTOM_ROUGHNESS_MINIMUM=1.0_SP
  BOTTOM_ROUGHNESS_LENGTHSCALE=1.0_SP
  CALL CNTRL_PRMTRS

  CALL CHECK_IO_DIRS

  CALL OPEN_FILES

  CALL READ_COLDSTART_SIGMA
  CLOSE(SIGMAUNIT)

  CALL READ_COLDSTART_GRID(GRIDUNIT,MGL,NGL,NV)
  m = MGL
  mt = MGL
  n = ngl
  nt = ngl


  CALL READ_COLDSTART_OBC_GRID(OBCUNIT,MGL,IOBCN,I_OBC_N, TYPE_OBC)
  CLOSE(OBCUNIT)
  
  CALL ALLOCATE_SPACE

  CALL READ_COLDSTART_COORDS(GRIDUNIT,MGL,VX,VY)
  CLOSE(GRIDUNIT)

  CALL READ_COLDSTART_DEPTH(DEPTHUNIT,MGL,VX,VY,H)
  CLOSE(DEPTHUNIT)

  KBM1 = KB - 1
  KBM2 = KB - 2
  CALL Setup_Sigma
  CALL SETUP_SIGMA_DERIVATIVES

  CALL SET_DIMENSIONS

!========================================================
!========================================================
!========================================================

  SELECT CASE(ELEVATION_SOURCE_TYPE)
  CASE("convert2new_spectral")
     write(ipt,*) "ELEVATION DATA: CONVERT 2 SPECTRAL NEW FVCOM FORMAT"

     CALL GET_OBC_TYPES
     write(ipt,*)'pass 1'

     CALL OPEN_ELV_SOURCE
     write(ipt,*)'pass 2'

     CALL READ_SPECTRAL
     write(ipt,*)'pass 3'

     CALL WRITE_SPECTRAL
     CLOSE(ELVUNIT)

  CASE("convert2new_julian")
     write(ipt,*) "ELEVATION DATA: CONVERT JULIAN 2 NEW FVCOM FORMAT"

     CALL GET_OBC_TYPES

     CALL OPEN_ELV_SOURCE

#    if defined (EQUI_TIDE)
     CALL READ_TIDAL_COMPONENTS
#    endif

     CALL READ_OLD_JULIAN

     CALL WRITE_JULIAN
     CLOSE(ELVUNIT)

  CASE("USER_DEFINED_JULIAN")
     write(ipt,*) "ELEVATION DATA: USER DEFINED JULIAN METHOD"

     CALL CREATE_JULIAN

     CALL WRITE_JULIAN

  CASE("USER_DEFINED_SPECTRAL")
     write(ipt,*) "ELEVATION DATA: USER DEFINED SPECTRAL METHOD"

     CALL CREATE_SPECTRAL

     CALL WRITE_SPECTRAL
 
  CASE DEFAULT
     CALL WARNING("NO ELEVATION TYPE SELECTED: MOVING ON")
  END SELECT


!========================================================
!========================================================
!========================================================
  SELECT CASE(TS_SOURCE_TYPE)
  CASE("convert2new")
     write(ipt,*) "TS DATA: CONVERT 2 NEW FVCOM FORMAT"

     CALL OPEN_TS_SOURCE
     
     CALL READ_OLD_TS

     CALL SET_OBC_DEPTH
     
     CALL WRITE_TSOBC

  CASE("USER_DEFINED")
     write(ipt,*) "TS DATA: USER DEFINED METHOD"

     CALL CREATE_TS
     
     CALL WRITE_TSOBC


  CASE DEFAULT
     CALL WARNING("NO TS TYPE SELECTED: MOVING ON")
  END SELECT

!========================================================
!========================================================
!========================================================

  
  IF (DBG_SET(DBG_LOG)) WRITE(IPT,*) "! TADA!!!"
  CALL PSHUTDOWN


END PROGRAM obc_create
