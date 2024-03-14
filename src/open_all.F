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

SUBROUTINE OPEN_ALL
!===============================================================================!
! OPEN FILES
! UNLESS OTHERWISE SPECIFED THE ROUTINES CALLED HERE ARE IN mod_input.F
!===============================================================================!
  USE CONTROL
  USE MOD_INPUT
  USE MOD_NESTING
  IMPLICIT NONE
  CHARACTER(LEN=160) :: FNAME
  LOGICAL :: FEXIST
  

  CALL NULLIFY_FILE_POINTERS

  SELECT CASE(STARTUP_TYPE)
     !=================================================
     ! HOTSTART 
  CASE(STARTUP_TYPE_HOTSTART)
     !=================================================
     if(DBG_SET(dbg_log)) then 
        WRITE(IPT,*)'!               OPEN INIT FILES FOR HOTSTART                       !'
        WRITE(IPT,*)'!                                                                !'
     end if
     
     CALL CHECK_IO_DIRS
     
     CALL OPEN_STARTUP_FILE
     
     CALL OPEN_FORCING

     CALL OPEN_NEW_OUTPUT
          
     IF(NESTING_ON) CALL OPEN_NESTING_FILE
#    if defined (WAVE_CURRENT_INTERACTION)
     IF(NESTING_ON_WAVE) CALL OPEN_NESTING_FILE_WAVE
#    endif     
#    if defined (VEGETATION)
     CALL OPEN_FORCING_VEGETATION
#    endif     
     !=================================================
     ! CRASHSTART 
  CASE(STARTUP_TYPE_CRASHRESTART)
     !=================================================
     if(DBG_SET(dbg_log)) then 
        WRITE(IPT,*)'!              OPENING FILES FOR CRASHRESTART                    !'
        WRITE(IPT,*)'!                                                                !'
     end if
 
     CALL CHECK_IO_DIRS

     CALL OPEN_CRASHSTART
     
     CALL OPEN_FORCING
     
     IF(NESTING_ON) CALL OPEN_NESTING_FILE
#    if defined (WAVE_CURRENT_INTERACTION)
     IF(NESTING_ON_WAVE) CALL OPEN_NESTING_FILE_WAVE
#    endif     
#    if defined (VEGETATION)
     CALL OPEN_FORCING_VEGETATION
#    endif     

     !=================================================
     ! COLDSTART 
  CASE(STARTUP_TYPE_COLDSTART)
     !=================================================
     if(DBG_SET(dbg_log)) then 
        WRITE(IPT,*)'!              OPENING FILES FOR COLDSTART                       !'
        WRITE(IPT,*)'!                                                                !'
     end if
          

     CALL CHECK_IO_DIRS

        ! MAKE SURE THE RUN FILE DOES NOT REQUEST A START FILE
        

     IF (STARTUP_TS_TYPE .eq. STARTUP_TYPE_OBSERVED) THEN
        CALL OPEN_STARTUP_FILE

     ELSE IF (STARTUP_TS_TYPE .eq. STARTUP_TYPE_SETVALUES) THEN
        CALL OPEN_STARTUP_FILE
        
     ELSE IF (STARTUP_UV_TYPE .eq. STARTUP_TYPE_SETVALUES) THEN
        CALL OPEN_STARTUP_FILE
        
     ELSE IF (STARTUP_TURB_TYPE .eq. STARTUP_TYPE_SETVALUES) THEN
        CALL OPEN_STARTUP_FILE
     ELSE 
        if(dbg_set(dbg_log)) write(ipt,*) "! No Startup file needed fo&
             &r this cold start"

     END IF
    
     ! OPEN THE OTHER COLD START FILES (GRID,DEPTH SPONGE, ETC)
     IF (MSR) CALL OPEN_COLDSTART ! ONLY MASTER READS THESE FILES
     
     CALL OPEN_FORCING

     CALL OPEN_NEW_OUTPUT 

     IF(NESTING_ON) CALL OPEN_NESTING_FILE
#    if defined (WAVE_CURRENT_INTERACTION)
     IF(NESTING_ON_WAVE) CALL OPEN_NESTING_FILE_WAVE
#    endif     
#    if defined (VEGETATION)
     CALL OPEN_FORCING_VEGETATION
#    endif     

  END SELECT



END SUBROUTINE OPEN_ALL
