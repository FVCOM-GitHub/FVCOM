PROGRAM sms2grid
  use mod_input
  use mod_sms2grid
  use mod_setup
  IMPLICIT NONE

  character(len=*),parameter::CVS_Id="$Id$" ! [sng] CVS Identification
  character(len=*),parameter::CVS_Date="$Date$" ! [sng] Date string
  character(len=*),parameter::CVS_Name="$Name$" ! [sng] File name string
  character(len=*),parameter::CVS_Revision="$Revision$" ! [sng] File revision string



  ! INTIALIZE MEMORY FROM libfvcom.a
  CALL INITIALIZE_CONTROL("sms2grid")

# if defined (MULTIPROCESSOR)
  CALL INIT_MPI_ENV(MYID,NPROCS,SERIAL,PAR,MSR,MSRID)
  MPI_FVCOM_GROUP = MPI_COMM_WORLD ! FOR NOW MAKE THEM EQUAL
# endif

  IF (PAR) CALL FATAL_ERROR &
       & ("THIS PROGRAM IS WAY TO SHORT TO BOTHER WITH WRITING AS A PARALLE CODE.",&
       & "STOP AND RERUN WITH JUST ONE PROCESSOR!")

  CALL INITIALIZE_NML

  CALL GET_COMMANDLINE(CVS_ID,CVS_Date,CVS_Name,CVS_Revision)

  IF (DBG_SET(DBG_LOG)) THEN
     WRITE(IPT,*) "! ========================="
     WRITE(IPT,*) "!   BEGIN SMS2GRID"
     WRITE(IPT,*) "! ========================="
  END IF

  CALL READ_NAMELIST

  CALL CHECK_IO_DIRS

  CALL OPEN_FILES

  CALL CONVERT_SMS2DAT

  
  IF (DBG_SET(DBG_LOG)) WRITE(IPT,*) "! TADA!!!"
  CALL PSHUTDOWN


END PROGRAM sms2grid
