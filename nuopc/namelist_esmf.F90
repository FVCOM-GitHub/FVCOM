!
!! @FVCOM and WAVE Coupling NUOPC Cap
!! @author Jianhua Qi (jqi@umassd.edu)
!! @date 12/20/2019 Original documentation
!------------------------------------------------------
!
SUBROUTINE NAMELIST_ESMF(NAME_DOMAIN)
  USE MOD_DRIVER
  USE MOD_UTILS
  USE CONTROL
  USE MOD_INPUT
  USE MOD_NESTING
  USE MOD_STATION_TIMESERIES  

  IMPLICIT NONE

  CHARACTER(LEN=20) NAME_DOMAIN
  INTEGER :: S_YEAR_I , E_YEAR_I
  INTEGER :: S_MONTH_I, E_MONTH_I
  INTEGER :: S_DAY_I  , E_DAY_I
  INTEGER :: S_HOUR_I , E_HOUR_I
  INTEGER :: S_MIN_I  , E_MIN_I

  !==============================================================================!
  ! SET DEFAULT VALUES IN NAME LIST                                                   
  !==============================================================================!
  CALL NAME_LIST_INITIALIZE

  CALL NAME_LIST_INITIALIZE_NEST
  
  CALL STATION_NAME_LIST_INITIALIZE 

  ! IF FVCOM IS ONLY PRINTING A BLANK NAME LIST FOR A NEW CASE:
  if (BLANK_NAMELIST) then
     CALL NAME_LIST_PRINT
     CALL STATION_NAME_LIST_PRINT 
     CALL PSHUTDOWN
  end if

  !==============================================================================!
  !   SETUP MODEL RUN PARAMETERS                                                 !
  !==============================================================================!

  !READ DATA IN THE NAME LIST FILE
  CASENAME = NAME_DOMAIN
  CALL NAME_LIST_READ ! ALL PROCS READ THIS

  READ(START_DATE(1:4),*)   S_YEAR_I
  READ(START_DATE(6:7),*)   S_MONTH_I
  READ(START_DATE(9:10),*)  S_DAY_I
  READ(START_DATE(12:13),*) S_HOUR_I
  READ(START_DATE(15:16),*) S_MIN_I
  
  READ(END_DATE(1:4),*)   E_YEAR_I
  READ(END_DATE(6:7),*)   E_MONTH_I
  READ(END_DATE(9:10),*)  E_DAY_I
  READ(END_DATE(12:13),*) E_HOUR_I
  READ(END_DATE(15:16),*) E_MIN_I
 
  IF(SYEAR  /= S_YEAR_I)  CALL FATAL_ERROR("The start year in domain "//trim(casename)//" does not match with the NUOPC/ESMF start year.")
  IF(SMONTH /= S_MONTH_I) CALL FATAL_ERROR("The start month in domain "//trim(casename)//" does not match with the NUOPC/ESMF start month.")
  IF(SDAY   /= S_DAY_I)   CALL FATAL_ERROR("The start day in domain "//trim(casename)//" does not match with the NUOPC/ESMF start day.")
  IF(SHOUR  /= S_HOUR_I)  CALL FATAL_ERROR("The start hour in domain "//trim(casename)//" does not match with the NUOPC/ESMF start hour.")
  IF(SMINUTE   /= S_MIN_I)   CALL FATAL_ERROR("The start minute in domain "//trim(casename)//" does not match with the NUOPC/ESMF start minute.")

!  IF(EYEAR  /= E_YEAR_I)  CALL FATAL_ERROR("The end year in domain "//trim(casename)//" does not match with the NUOPC/ESMF end year.")
!  IF(EMONTH /= E_MONTH_I) CALL FATAL_ERROR("The end month in domain "//trim(casename)//" does not match with the NUOPC/ESMF end month.")
!  IF(EDAY   /= E_DAY_I)   CALL FATAL_ERROR("The end day in domain "//trim(casename)//" does not match with the NUOPC/ESMF end day.")
!  IF(EHOUR  /= E_HOUR_I)  CALL FATAL_ERROR("The end hour in domain "//trim(casename)//" does not match with the NUOPC/ESMF end hour.")
!  IF(EMINUTE   /= E_MIN_I)   CALL FATAL_ERROR("The end minute in domain "//trim(casename)//" does not match with the NUOPC/ESMF end minute.")

  CALL NAME_LIST_READ_NEST
  IF(NESTING_ON .AND. SERIAL)THEN
    IF(MSR) WRITE(*,*) 'PLEASE USE MORE THAN ONE PROCESSOR TO RUN NESTING. STOP RUNNING...'
    CALL PSTOP
  END IF
  IF(NCNEST_ON .AND. SERIAL)THEN
    IF(MSR) WRITE(*,*) 'PLEASE USE MORE THAN ONE PROCESSOR TO RUN NCNEST. STOP RUNNING...'
    CALL PSTOP
  END IF

  CALL STATION_NAME_LIST_READ
  
  !PRINT THE NAME LIST DATA TO THE SCREEN FOR THE LOG
  IF(DBG_SET(DBG_LOG)) CALL NAME_LIST_PRINT 

END SUBROUTINE NAMELIST_ESMF
