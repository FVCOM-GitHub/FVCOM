PROGRAM Interp_Data
  use mod_interp_time
  USE MOD_INPUT
  IMPLICIT NONE

  character(len=*),parameter::CVS_Id="$Id$" ! [sng] CVS Identification
  character(len=*),parameter::CVS_Date="$Date$" ! [sng] Date string
  character(len=*),parameter::CVS_Name="$Name$" ! [sng] File name string
  character(len=*),parameter::CVS_Revision="$Revision$" ! [sng] File revision string



  ! INTIALIZE MEMORY FROM libfvcom.a
  CALL INITIALIZE_CONTROL("Interp Time")

# if defined (MULTIPROCESSOR)
  CALL INIT_MPI_ENV(MYID,NPROCS,SERIAL,PAR,MSR,MSRID)
# endif

  IF (PAR) CALL FATAL_ERROR &
       & ("THIS PROGRAM IS WAY TO SHORT TO BOTHER WITH WRITING AS A PARALLE CODE.",&
       & "STOP AND RERUN WITH JUST ONE PROCESSOR!")

  CALL GET_COMMANDLINE(CVS_ID,CVS_Date,CVS_Name,CVS_Revision)

  IF (DBG_SET(DBG_LOG)) THEN
     WRITE(IPT,*) "! ===================================="
     WRITE(IPT,*) "!   BEGIN Interpolation of Time data"
     WRITE(IPT,*) "! ===================================="
  END IF

  write(IPT,*) "! READ COMMAND LINE:"
  write(IPT,*) "! INPUT    :"//TRIM(FIN)
  write(IPT,*) "! OUTPUT   :"//TRIM(FOUT)
  write(IPT,*) "! START    :"//TRIM(START_DATE)
  write(IPT,*) "! END      :"//TRIM(END_DATE)
  write(IPT,*) "! TIMEZONE :"//TRIM(TIMEZONE)
  write(IPT,*) "! TIMESTEP :"//TRIM(TIMESTEP)

  CALL INTERP_TIME

  
  IF (DBG_SET(DBG_LOG)) WRITE(IPT,*) "! TADA!!!"
  CALL PSHUTDOWN


END PROGRAM Interp_Data
