!/===========================================================================/
! Copyright (c) 2007, The University of Massachusetts Dartmouth 
! Produced at the School of Marine Science & Technology 
! Marine Ecosystem Dynamics Modeling group
! All rights reserved.
!
! FVCOM has been developed by the joint UMASSD-WHOI research team. For 
! details of authorship and attribution of credit please see the FVCOM
! technical manual or contact the MEDM group.
!
! 
! This file is part of FVCOM. For details, see http://fvcom.smast.umassd.edu 
! The full copyright notice is contained in the file COPYRIGHT located in the 
! root directory of the FVCOM code. This original header must be maintained
! in all distributed versions.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
! AND ANY EXPRESS OR  IMPLIED WARRANTIES, INCLUDING,  BUT NOT  LIMITED TO,
! THE IMPLIED WARRANTIES OF MERCHANTABILITY AND  FITNESS FOR A PARTICULAR
! PURPOSE ARE DISCLAIMED.  
!
!/---------------------------------------------------------------------------/
! CVS VERSION INFORMATION
! $Id$
! $Name$
! $Revision$
!/===========================================================================/

SUBROUTINE NAMELIST
  USE MOD_UTILS
  USE CONTROL
  USE MOD_INPUT
  USE MOD_NESTING
  USE MOD_STATION_TIMESERIES  
  USE MOD_SPARSE_TIMESERIES
# if defined (SEMI_IMPLICIT)
  USE MOD_SEMI_IMPLICIT
# endif
# if defined (BALANCE_2D)
  USE MOD_BALANCE_2D
# endif  
# if defined (ONE_D_MODEL)
  USE MOD_ONEDTIDE
# endif
# if defined (WATER_QUALITY)
  USE MOD_WQM
# endif    
# if defined (DYE_RELEASE)
  USE MOD_DYE
# endif
# if defined (BioGen)
!  USE MOD_BIO_3D
# endif    
# if defined (NH)
  USE NON_HYDRO, ONLY: NAME_LIST_INITIALIZE_NH, NAME_LIST_PRINT_NH, NAME_LIST_READ_NH
# endif
# if defined (PWP)
  USE MOD_PWP, ONLY: NAME_LIST_INITIALIZE_PWP, NAME_LIST_PRINT_PWP, NAME_LIST_READ_PWP
# endif
# if defined (RRKF)
  USE MOD_RRK
# endif 
# if defined (ENKF)
  USE MOD_ENKF
# endif 
# if defined (DATA_ASSIM)
  USE MOD_ASSIM
# endif 
# if defined (HEATING_CALCULATED)
  USE MOD_HEATFLUX
# endif
# if defined (HEATING_SOLAR)
  USE MOD_SOLAR
# endif

  IMPLICIT NONE


  !==============================================================================!
  ! SET DEFAULT VALUES IN NAME LIST                                                   
  !==============================================================================!
  CALL NAME_LIST_INITIALIZE

  CALL NAME_LIST_INITIALIZE_NEST

# if defined (SEMI_IMPLICIT)    
  CALL NAME_LIST_INITIALIZE_SEMI
# endif
# if defined (DYE_RELEASE)    
  CALL NAME_LIST_INITIALIZE_DYE
# endif
# if defined (BALANCE_2D)    
  CALL NAME_LIST_INITIALIZE_BALANCE
# endif
# if defined (ONE_D_MODEL)    
  CALL NAME_LIST_INITIALIZE_ONED
# endif
# if defined (WATER_QUALITY)    
  CALL NAME_LIST_INITIALIZE_WQM
# endif
# if defined (BioGen)    
!  CALL NAME_LIST_INITIALIZE_BIO
# endif
# if defined (NH)
  CALL NAME_LIST_INITIALIZE_NH
# endif
# if defined (PWP)
  CALL NAME_LIST_INITIALIZE_PWP
# endif
# if defined (RRKF)
  CALL NAME_LIST_INITIALIZE_RRKF
# endif
# if defined (ENKF)
  CALL NAME_LIST_INITIALIZE_ENKF
# endif  
# if defined (DATA_ASSIM)
  CALL NAME_LIST_INITIALIZE_ASSIM
# endif 
# if defined (HEATING_CALCULATED)
  CALL HEATING_CALCULATE_NAMELIST_INITIALIZE
# endif
# if defined (HEATING_SOLAR)
  CALL HEATING_SOLAR_NAMELIST_INITIALIZE
# endif 

  CALL STATION_NAME_LIST_INITIALIZE 
# if defined (WAVE_CURRENT_INTERACTION)
  CALL SPARSE_NAME_LIST_INITIALIZE
# endif 

  ! IF FVCOM IS ONLY PRINTING A BLANK NAME LIST FOR A NEW CASE:
  if (BLANK_NAMELIST) then
     CALL NAME_LIST_PRINT

# if defined(MULTIPROCESSOR)
     ! NESTING ONLY WORKS IN PARALLEL
     CALL NAME_LIST_PRINT_NEST
# endif

#    if defined (SEMI_IMPLICIT)    
     CALL NAME_LIST_PRINT_SEMI
#    endif
#    if defined (DYE_RELEASE)    
     CALL NAME_LIST_PRINT_DYE
#    endif
#    if defined (BALANCE_2D)    
     CALL NAME_LIST_PRINT_BALANCE
#    endif
#    if defined (ONE_D_MODEL)    
     CALL NAME_LIST_PRINT_ONED
#    endif
#    if defined (WATER_QUALITY)    
     CALL NAME_LIST_PRINT_WQM
#    endif
#    if defined (BioGen)    
!     CALL NAME_LIST_PRINT_BIO
#    endif
#    if defined (NH)
     CALL NAME_LIST_PRINT_NH
#    endif
#    if defined (PWP)
     CALL NAME_LIST_PRINT_PWP
#    endif
#    if defined (RRKF)
     CALL NAME_LIST_PRINT_RRKF
#    endif
#    if defined (ENKF)
     CALL NAME_LIST_PRINT_ENKF
#    endif     
#    if defined (DATA_ASSIM)
     CALL NAME_LIST_PRINT_ASSIM
#    endif     
#    if defined (HEATING_CALCULATED)
     CALL HEATING_CALCULATE_NAMELIST_PRINT
#    endif     
#    if defined (HEATING_SOLAR)
     CALL HEATING_SOLAR_NAMELIST_PRINT
#    endif     

     CALL STATION_NAME_LIST_PRINT 
#    if defined (WAVE_CURRENT_INTERACTION)
     CALL SPARSE_NAME_LIST_PRINT
#    endif 

     CALL PSHUTDOWN
  end if

  !==============================================================================!
  !   SETUP MODEL RUN PARAMETERS                                                 !
  !==============================================================================!

  !READ DATA IN THE NAME LIST FILE
  CALL NAME_LIST_READ ! ALL PROCS READ THIS

  CALL NAME_LIST_READ_NEST
  IF(NESTING_ON .AND. SERIAL)THEN
    IF(MSR) WRITE(*,*) 'PLEASE USE MORE THAN ONE PROCESSOR TO RUN NESTING. STOP RUNNING...'
    CALL PSTOP
  END IF
  IF(NCNEST_ON .AND. SERIAL)THEN
    IF(MSR) WRITE(*,*) 'PLEASE USE MORE THAN ONE PROCESSOR TO RUN NCNEST. STOP RUNNING...'
    CALL PSTOP
  END IF

# if defined (WAVE_CURRENT_INTERACTION)
  IF(NESTING_ON_WAVE .AND. SERIAL)THEN
    IF(MSR) WRITE(*,*) 'PLEASE USE MORE THAN ONE PROCESSOR TO RUN WAVE NESTING. STOP RUNNING...'
    CALL PSTOP
  END IF
  IF(NCNEST_ON_WAVE .AND. SERIAL)THEN
    IF(MSR) WRITE(*,*) 'PLEASE USE MORE THAN ONE PROCESSOR TO RUN WAVE NCNEST. STOP RUNNING...'
    CALL PSTOP
  END IF
# endif  

# if defined (SEMI_IMPLICIT)    
  CALL NAME_LIST_READ_SEMI
# endif
# if defined (DYE_RELEASE)    
  CALL NAME_LIST_READ_DYE
# endif
# if defined (BALANCE_2D)    
  CALL NAME_LIST_READ_BALANCE
# endif
# if defined (ONE_D_MODEL)    
  CALL GET_UMAX
# endif
# if defined (WATER_QUALITY)    
  CALL NAME_LIST_READ_WQM
# endif
# if defined (NH)
  CALL NAME_LIST_READ_NH
# endif
# if defined (RRKF)
  CALL NAME_LIST_READ_RRKF
# endif
# if defined (ENKF)
  CALL NAME_LIST_READ_ENKF
# endif
# if defined (HEATING_CALCULATED)
  CALL HEATING_CALCULATE_NAMELIST_READ
# endif
# if defined (HEATING_SOLAR)
  CALL HEATING_SOLAR_NAMELIST_READ
# endif

  CALL STATION_NAME_LIST_READ
# if defined (WAVE_CURRENT_INTERACTION)
  CALL SPARSE_NAME_LIST_READ  
# endif

  !PRINT THE NAME LIST DATA TO THE SCREEN FOR THE LOG
  IF(DBG_SET(DBG_LOG)) CALL NAME_LIST_PRINT 


END SUBROUTINE NAMELIST
