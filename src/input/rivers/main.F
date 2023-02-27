PROGRAM river_files
  use mod_input
  use mod_rivers
  IMPLICIT NONE

  character(len=*),parameter::CVS_Id="$Id$" ! [sng] CVS Identification
  character(len=*),parameter::CVS_Date="$Date$" ! [sng] Date string
  character(len=*),parameter::CVS_Name="$Name$" ! [sng] File name string
  character(len=*),parameter::CVS_Revision="$Revision$" ! [sng] File revision string



  ! INTIALIZE MEMORY FROM libfvcom.a
  CALL INITIALIZE_CONTROL("RIVER2FVCOM")

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
     WRITE(IPT,*) "!   BEGIN RIVER 2 FVCOM"
     WRITE(IPT,*) "! ========================="
  END IF
  
  CALL NAME_LIST_INITIALIZE

  CALL READ_NAMELIST

  CALL CNTRL_PRMTRS

  CALL SET_TIME

  CALL CHECK_IO_DIRS

  CALL LOAD_RIVER_DATA

  CALL DUMP_OUTFILE


  IF (DBG_SET(DBG_LOG)) WRITE(IPT,*) "! TADA!!!"
  CALL PSHUTDOWN


END PROGRAM river_files
