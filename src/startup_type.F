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

SUBROUTINE SET_STARTUP_TYPE
  USE CONTROL
  USE MOD_UTILS
# if defined (WATER_QUALITY)
  USE MOD_WQM
# endif    
# if defined (BioGen)
  USE MOD_BIO_3D
# endif    
  IMPLICIT NONE
  INTEGER :: II

  IF(STARTUP_TYPE .EQ. STARTUP_TYPE_FORECAST) THEN
     FORECAST_MODE = .TRUE.
     
     IF (dbg_set(dbg_log)) then
        WRITE(IPT,*) "! ======================================================="
        WRITE(IPT,*) "!!!!!!      STARTING FVCOM FORECAST MODE            !!!!!"
        WRITE(IPT,*) "! ======================================================="
     END IF
     STARTUP_TYPE = STARTUP_TYPE_CRASHRESTART


  ELSE
     FORECAST_MODE = .FALSE.
  END IF



  IF(CMDLN_RESTART)THEN

     IF (FORECAST_MODE) THEN
        CALL WARNING("CRASH RESTART DOES NOT WORK IN FORECAST MODE!")

     ELSE
        IF (dbg_set(dbg_log)) then
           WRITE(IPT,*) "! ======================================================="
           WRITE(IPT,*) "!!!!!!OVER RIDING NAMELIST WITH COMMAND LINE RESTART!!!!!"
           WRITE(IPT,*) "! ======================================================="
        END IF
        STARTUP_TYPE = STARTUP_TYPE_CRASHRESTART
     END IF
  END IF



  SELECT CASE(STARTUP_TYPE)
     !=================================================
     ! HOTSTART 
  CASE(STARTUP_TYPE_HOTSTART)
     !=================================================
     if(DBG_SET(dbg_log)) then 
        WRITE(IPT,*)'!               RUNNING HOTSTART                                 !'
        WRITE(IPT,*)'!                                                                !'
     end if
     
          
     !=================================================
     ! CRASHSTART 
  CASE(STARTUP_TYPE_CRASHRESTART)
     !=================================================

     IF (.not. FORECAST_MODE) THEN
        if(DBG_SET(dbg_log)) then 
           WRITE(IPT,*)'!              RUNNING CRASHRESTART                              !'
           WRITE(IPT,*)'!                                                                !'
        end if
     END IF

     STARTUP_TS_TYPE = STARTUP_TYPE_SETVALUES
     STARTUP_UV_TYPE = STARTUP_TYPE_SETVALUES
     STARTUP_TURB_TYPE = STARTUP_TYPE_SETVALUES

     !=================================================
     ! COLDSTART 
  CASE(STARTUP_TYPE_COLDSTART)
     !=================================================
     if(DBG_SET(dbg_log)) then 
        WRITE(IPT,*)'!                 RUNNING COLDSTART                              !'
        WRITE(IPT,*)'!                                                                !'
     end if
          

  END SELECT
  
  ! CHECK FOR UNKNOW STARTUP TYPES


  SELECT CASE (STARTUP_UV_TYPE)
  CASE (STARTUP_TYPE_OBSERVED)
     CALL FATAL_ERROR("I DON'T KNOW HOW TO DO THAT KIND OF STARTUP_TYPE_OBSERVED1")
  CASE(STARTUP_TYPE_LINEAR) 
     CALL FATAL_ERROR("I DON'T KNOW HOW TO DO THAT KIND OF STARTUP_TYPE_LINEAR2")
  CASE(STARTUP_TYPE_CONSTANT)
     !OKAY 
     !CALL FATAL_ERROR("I DON'T KNOW HOW TO DO THAT KIND OF STARTUP_TYPE_CONSTANT3")
  CASE(STARTUP_TYPE_DEFAULT)
     ! OKAY
  CASE(STARTUP_TYPE_SETVALUES)
     ! OKAY
  CASE DEFAULT
     CALL FATAL_ERROR("UNKNOWN STARTUP_UV_TYPE")
  END SELECT



  SELECT CASE(STARTUP_TURB_TYPE)
  CASE(STARTUP_TYPE_OBSERVED) 
     CALL FATAL_ERROR("I DON'T KNOW HOW TO DO THAT KIND OF STARTUP_TYPE_OBSERVED4")
  CASE(STARTUP_TYPE_LINEAR)
     CALL FATAL_ERROR("I DON'T KNOW HOW TO DO THAT KIND OF STARTUP_TYPE_LINEAR5")
  CASE(STARTUP_TYPE_CONSTANT)
     CALL FATAL_ERROR("I DON'T KNOW HOW TO DO THAT KIND OF STARTUP_TYPE_CONSTANT6")
  CASE(STARTUP_TYPE_DEFAULT)
     ! OKAY
  CASE(STARTUP_TYPE_SETVALUES)
     ! OKAY
  CASE DEFAULT
     CALL FATAL_ERROR("UNKNOWN STARTUP_TURB_TYPE")
  END SELECT
   

  SELECT CASE(STARTUP_TS_TYPE)
  CASE(STARTUP_TYPE_OBSERVED) 
     IF(BAROTROPIC) CALL FATAL_ERROR("CAN'T SET OBSERVERD TS VALUES IN BAROTROPIC MODE!")
     ! OKAY
  CASE(STARTUP_TYPE_LINEAR)
     IF(BAROTROPIC) CALL FATAL_ERROR("CAN'T SET LINEAR TS VALUES IN BAROTROPIC MODE!")
     ! OKAY
     IF (STARTUP_T_VALS(1) == -99.0_SP) CALL FATAL_ERROR("STARTUP_T_VAL not set in run file",&
          & "two values required for linear startup (default is -99.0)")
     IF (STARTUP_T_VALS(2) == -99.0_SP)CALL FATAL_ERROR("STARTUP_T_VAL not set in run file",&
          & "two values required for linear startup (default is -99.0)")

     IF (STARTUP_S_VALS(1) == -99.0_SP)CALL FATAL_ERROR("STARTUP_S_VAL not set in run file",&
          & "two values required for linear startup (default is -99.0)")
     IF (STARTUP_S_VALS(2) == -99.0_SP)CALL FATAL_ERROR("STARTUP_S_VAL not set in run file",&
          & "two values required for linear startup (default is -99.0)")

     IF (STARTUP_DMAX == -99.0_SP)CALL FATAL_ERROR("STARTUP_DMAX not set in run file",&
          & "two values required for linear startup (default is -99.0)")


  CASE(STARTUP_TYPE_CONSTANT)
     ! OKAY
     IF (STARTUP_T_VALS(1) == -99.0_SP)CALL FATAL_ERROR("STARTUP_T_VAL not set in run file",&
          & "one values required for constant startup (default is -99.0)") 
     IF (STARTUP_T_VALS(2) /= -99.0_SP)CALL FATAL_ERROR("STARTUP_T_VAL is incorrect run file:",&
          & "only allowed for constant startup!")

     IF (STARTUP_S_VALS(1) == -99.0_SP) CALL FATAL_ERROR("STARTUP_S_VAL not set in run file",&
          & "one values required for constant startup (default is -99.0)")
     IF (STARTUP_S_VALS(2) /= -99.0_SP)CALL FATAL_ERROR("STARTUP_S_VAL is incorrect run file:",&
          & "only allowed for constant startup!")

  CASE(STARTUP_TYPE_DEFAULT)
     CALL FATAL_ERROR("I DON'T KNOW HOW TO DO THAT KIND OF STARTUP_TYPE_DEFAULT7")
  CASE(STARTUP_TYPE_SETVALUES)
     IF(BAROTROPIC) CALL WARNING("YOU ARE RESTARTING A BAROTROPIC CASE&
          &: T&S BETTER BE CONSTANT IN THE RESTART FILE!")
     ! OKAY
  CASE DEFAULT
     CALL FATAL_ERROR("UNKNOWN STARTUP_TS_TYPE")
  END SELECT
   
# if defined (WATER_QUALITY)
  SELECT CASE(STARTUP_WQM_TYPE)
  CASE(STARTUP_TYPE_OBSERVED) 
     IF(BAROTROPIC) CALL FATAL_ERROR("CAN'T SET OBSERVERD WQM VALUES IN BAROTROPIC MODE!")
     ! OKAY
  CASE(STARTUP_TYPE_LINEAR)
     IF(BAROTROPIC) CALL FATAL_ERROR("CAN'T SET LINEAR WQM VALUES IN BAROTROPIC MODE!")
     ! OKAY
     DO II = 1,NB
     IF (STARTUP_WQM_VALS(1,II) == -99.0_SP) CALL FATAL_ERROR("STARTUP_WQM_VAL not set in run file",&
          & "two values required for linear startup (default is -99.0)")
     IF (STARTUP_WQM_VALS(2,II) == -99.0_SP) CALL FATAL_ERROR("STARTUP_WQM_VAL not set in run file",&
          & "two values required for linear startup (default is -99.0)")
     END DO
  CASE(STARTUP_TYPE_CONSTANT)
     ! OKAY
     DO II = 1,NB
     IF (STARTUP_WQM_VALS(1,II) == -99.0_SP)CALL FATAL_ERROR("STARTUP_WQM_VAL not set in run file",&
          & "one values required for constant startup (default is -99.0)") 
     IF (STARTUP_WQM_VALS(2,II) /= -99.0_SP)CALL FATAL_ERROR("STARTUP_WQM_VAL is incorrect run file:",&
          & "only allowed for constant startup!")
     END DO	  

  CASE(STARTUP_TYPE_DEFAULT)
     CALL FATAL_ERROR("I DON'T KNOW HOW TO DO THAT KIND OF STARTUP_TYPE_DEFAULT8")
  CASE(STARTUP_TYPE_SETVALUES)
     IF(BAROTROPIC) CALL WARNING("YOU ARE RESTARTING A BAROTROPIC CASE&
          &: WQM VARIABLES BETTER BE CONSTANT IN THE RESTART FILE!")
     ! OKAY
  CASE DEFAULT
     CALL FATAL_ERROR("UNKNOWN STARTUP_WQM_TYPE")
  END SELECT
# endif   

# if defined (BioGen)
  IF(BIOLOGICAL_MODEL)THEN

  SELECT CASE(STARTUP_BIO_TYPE)
  CASE(STARTUP_TYPE_OBSERVED) 
     IF(BAROTROPIC) CALL FATAL_ERROR("CAN'T SET OBSERVERD BIO VALUES IN BAROTROPIC MODE!")
     ! OKAY
  CASE(STARTUP_TYPE_LINEAR)
     IF(BAROTROPIC) CALL FATAL_ERROR("CAN'T SET LINEAR BIO VALUES IN BAROTROPIC MODE!")
     ! OKAY
  CASE(STARTUP_TYPE_CONSTANT)
     ! OKAY
  CASE(STARTUP_TYPE_DEFAULT)
     CALL FATAL_ERROR("I DON'T KNOW HOW TO DO THAT KIND OF STARTUP_TYPE_DEFAULT9")
  CASE(STARTUP_TYPE_SETVALUES)
     IF(BAROTROPIC) CALL WARNING("YOU ARE RESTARTING A BAROTROPIC CASE&
          &: BIO VARIABLES BETTER BE CONSTANT IN THE RESTART FILE!")
     ! OKAY
  CASE DEFAULT
     CALL FATAL_ERROR("UNKNOWN STARTUP_BIO_TYPE")
  END SELECT

  END IF
# endif   



END SUBROUTINE SET_STARTUP_TYPE
