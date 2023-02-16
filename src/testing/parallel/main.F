PROGRAM Parallel
  USE MOD_parallel
  ! SEE PARAMETER SETTINGS IN MOD_SST2GRD.F90
  IMPLICIT NONE
 
  character(len=*),parameter::CVS_Id="$Id$" ! [sng] CVS Identification
  character(len=*),parameter::CVS_Date="$Date$" ! [sng] Date string
  character(len=*),parameter::CVS_Name="$Name$" ! [sng] File name string
  character(len=*),parameter::CVS_Revision="$Revision$" ! [sng] File revision string



  INTEGER :: I

  ! INTIALIZE MEMORY FROM libfvcom.a
  CALL INITIALIZE_CONTROL("Parallel")

#if defined(MULTIPROCESSOR)

  CALL INIT_MPI_ENV(MYID,NPROCS,SERIAL,PAR,MSR,MSRID)
# endif
  
  CALL GET_COMMANDLINE(CVS_ID,CVS_Date,CVS_Name,CVS_Revision)


  IF (.NOT. PAR) CALL FATAL_ERROR("THIS IS A TEST PROGRAM DESIGNED ONLY TO RUN IN PARRAL!")


  IF (DBG_SET(DBG_LOG)) THEN
     WRITE(IPT,*) "! ===================================="
     WRITE(IPT,*) "!   FVCOM PARALLEL COMMUNICATION TEST"
     WRITE(IPT,*) "! ===================================="
  END IF
  
  CALL GET_FVCOM_GRID

  IF (DBG_SET(DBG_LOG)) WRITE(IPT,*) "! LOADED GRID AND SET UP PARALLEL"
  
  IF (DBG_SET(DBG_LOG)) WRITE(IPT,*) "! RUNNING Par_Test:"

  DO I =0,4
     KB = I *10

     WRITE(IPT,*) "! ===================================="
     WRITE(IPT,*) "!   FVCOM PARALLEL TEST: ",I
     WRITE(IPT,*) "! ===================================="

     Call Par_test


  END DO

  IF (DBG_SET(DBG_LOG)) WRITE(IPT,*) "! TADA: All done!"
  CALL PSHUTDOWN



END PROGRAM Parallel

