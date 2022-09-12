PROGRAM Modify_grid
  use mod_grid
  IMPLICIT NONE
  character(len=*),parameter::CVS_Id="$Id$" ! [sng] CVS Identification
  character(len=*),parameter::CVS_Date="$Date$" ! [sng] Date string
  character(len=*),parameter::CVS_Name="$Name$" ! [sng] File name string
  character(len=*),parameter::CVS_Revision="$Revision$" ! [sng] File revision string



  ! INTIALIZE MEMORY FROM libfvcom.a
  CALL INITIALIZE_CONTROL("Modify Grid")

# if defined (MULTIPROCESSOR)
  CALL INIT_MPI_ENV(MYID,NPROCS,SERIAL,PAR,MSR,MSRID)
  MPI_FVCOM_GROUP = MPI_COMM_WORLD ! FOR NOW MAKE THEM EQUAL
# endif

  IF (PAR) CALL FATAL_ERROR &
       & ("THIS PROGRAM IS WAY TO SHORT TO BOTHER WITH WRITING AS A PARALLE CODE.",&
       & "STOP AND RERUN WITH JUST ONE PROCESSOR!")


  ! INITIALIZE BOTH THE FVCOM MOD_MAIN VARIABLES AND THE LOCAL VARS
  CALL NAME_LIST_INITIALIZE
  CALL INITIALIZE_NML

  CALL GET_COMMANDLINE(CVS_ID,CVS_Date,CVS_Name,CVS_Revision)

  IF (DBG_SET(DBG_LOG)) THEN
     WRITE(IPT,*) "! ========================="
     WRITE(IPT,*) "!   BEGIN Modify Grid"
     WRITE(IPT,*) "! ========================="
  END IF


  ! READ THE NAME LIST
  CALL READ_NAMELIST

  CALL CNTRL_PRMTRS

  CALL CHECK_IO_DIRS

  CALL OPEN_FILES

  CALL SET_TIME_INDEX

  WRITE(IPT,*) "! ========================="
  WRITE(IPT,*) "! LOADING NEW GRID DATA"
  WRITE(IPT,*) "! ========================="
  
  CALL READ_NEW_GRID

  WRITE(IPT,*) "! ========================="
  WRITE(IPT,*) "! FINISHED LOADING NEW GRID DATA"
  WRITE(IPT,*) "! LOAD OLD GRID DATA!"
  WRITE(IPT,*) "! ========================="
  
  CALL READ_OLD_GRID

  WRITE(IPT,*) "! ========================="
  WRITE(IPT,*) "! FINISHED LOADING OLD GRID DATA!"
  WRITE(IPT,*) "! CREATE GRID SETUP!"
  WRITE(IPT,*) "! ========================="
  
  CALL MY_GRID_METRICS

  WRITE(IPT,*) "! ========================="
  WRITE(IPT,*) "! FINISHED GRID SETUP!"
  WRITE(IPT,*) "! CREATE COEFFICENTS FOR INTERPOLATION!"
  WRITE(IPT,*) "! ========================="

  CALL INTERP_COEFFICIENTS


!  write(ipt,*) n2o_nid

  CALL CREATE_NEW_FILE

  
  CALL DUMP_NEW_GRID_DATA


  IF (DBG_SET(DBG_LOG)) WRITE(IPT,*) "! TADA!!!"
  CALL PSHUTDOWN


END PROGRAM Modify_grid
