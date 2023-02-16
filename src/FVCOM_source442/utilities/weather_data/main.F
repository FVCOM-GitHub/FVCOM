PROGRAM Weather_Data
  use mod_weather
  IMPLICIT NONE

  character(len=*),parameter::CVS_Id="$Id$" ! [sng] CVS Identification
  character(len=*),parameter::CVS_Date="$Date$" ! [sng] Date string
  character(len=*),parameter::CVS_Name="$Name$" ! [sng] File name string
  character(len=*),parameter::CVS_Revision="$Revision$" ! [sng] File revision string



  ! INTIALIZE MEMORY FROM libfvcom.a
  CALL INITIALIZE_CONTROL("Weather Data")

# if defined (MULTIPROCESSOR)
  CALL INIT_MPI_ENV(MYID,NPROCS,SERIAL,PAR,MSR,MSRID)
# endif

  IF (PAR) CALL FATAL_ERROR &
       & ("THIS PROGRAM IS WAY TO SHORT TO BOTHER WITH WRITING AS A PARALLE CODE.",&
       & "STOP AND RERUN WITH JUST ONE PROCESSOR!")


  FIN=""
  FOUT=""
  FSET=""

  CALL GET_COMMANDLINE(CVS_ID,CVS_Date,CVS_Name,CVS_Revision)

  IF (DBG_SET(DBG_LOG)) THEN
     WRITE(IPT,*) "! ========================="
     WRITE(IPT,*) "!   BEGIN Weather Data"
     WRITE(IPT,*) "! ========================="
  END IF

  Call READ_INPUT

  Call BUILD_GRID

  Call BUILD_OUTPUT

  CALL SETUP_CONVERT

  CALL EXTRACT

  
  IF (DBG_SET(DBG_LOG)) WRITE(IPT,*) "! TADA!!!"
  CALL PSHUTDOWN


END PROGRAM Weather_Data
