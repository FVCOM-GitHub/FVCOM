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

SUBROUTINE CNTRL_PRMTRS
  USE LIMS
  USE CONTROL
  USE MOD_UTILS
  USE MOD_ASSIM
  USE MOD_TIME
  USE MOD_NCDIO
  USE MOD_NESTING
# if defined (HEATING_CALCULATED)
  USE MOD_HEATFLUX
# endif 
  IMPLICIT NONE
  LOGICAL TEST
  integer :: stat,ierr


  ! SIMPLE SETTINGS FROM RUN FILE ETC SHOULD GO HERE!

  IF (MSR) USE_PROJ = HAVE_PROJ(PROJECTION_REFERENCE)
# if defined (MULTIPROCESSOR)
  IF(PAR) CALL MPI_BCAST(USE_PROJ,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
# endif
  IF (USE_PROJ) THEN
     if(dbg_set(dbg_log)) write(IPT,*) "! PROJ 4 CARTOGRAPHIC PROJECTION ENABLED!"
  ELSE
     if(dbg_set(dbg_log)) write(IPT,*) "! PROJ 4 CARTOGRAPHIC PROJECTION DISABLED!"
  END IF


# if !defined(WET_DRY)
  If(WETTING_DRYING_ON) CALL FATAL_ERROR &
       & ('You must compile with wet dry to use it!',&
       & 'Recompile after selecting it in the make file,',&
       & 'or turn it off in the run file')
# endif

# if !defined(ICE)
  If(ICE_MODEL) CALL FATAL_ERROR &
       & ('You must compile with ICE to use it!',&
       & 'Recompile after selecting it in the make file,',&
       & 'or turn it off in the run file')
# endif



!  write(ipt,*)  "PRG_NAME='"//trim(prg_name)//"'"


  IF(TRIM(PRG_NAME) == "FVCOM") THEN

     ! TRANSLATE TO OLD FVCOM PARAMTER NAMES
     
     ! SET THE FVCOM NAMES FOR MIXING TERMS
     UMOL   = VERTICAL_MIXING_COEFFICIENT
     if (UMOL .LT. 0) CALL FATAL_ERROR&
          &("Name List Error: VERTICAL_MIXING_COEFFICIENT outside valid range")
     VPRNU  = VERTICAL_PRANDTL_NUMBER
     if (VPRNU .LT. 0) CALL FATAL_ERROR&
          &("Name List Error: VERTICAL_PRANDTL_NUMBER outside valid range")
     
     if (HORIZONTAL_MIXING_KIND .eq. STTC) THEN
        if (HORIZONTAL_MIXING_COEFFICIENT .NE. -1.0_SP)THEN
           CALL WARNING("The name list variable: HORIZONTAL_MIXING_COEFFICIENT",&
                "should have the default value when the horizontal mixing is set using 'static'",&
                "(ie set from an input file!)")
           
           HORIZONTAL_MIXING_COEFFICIENT = -1.0_SP
        end if
        
     else if(HORIZONTAL_MIXING_KIND .eq. CNSTNT) THEN
        
        if (HORIZONTAL_MIXING_COEFFICIENT .LT. 0.0_SP) CALL FATAL_ERROR&
             &("The name list variable: HORIZONTAL_MIXING_COEFFICIENT",&
             & "Can not have the default value -1.0 when the horizontal mixing is set using 'constant'",&
             & "(ie set the viscosity in the model!)")
     else
        CALL FATAL_ERROR&
             &("The name list variable: HORIZONTAL_MIXING_KIND",&
             & "Is set to an invalid option: (choices: constant, static)")     
     end if
     
     HPRNU  = HORIZONTAL_PRANDTL_NUMBER
     if (HPRNU .LT. 0) CALL FATAL_ERROR&
          &("Name List Error: HORIZONTAL_PRANDTL_NUMBER outside valid range")
     
     
     CBCMIN = BOTTOM_ROUGHNESS_MINIMUM
     if (CBCMIN .LT. 0) CALL FATAL_ERROR&
          &("Name List Error: BOTTOM_ROUGHNESS_MINIMUM outside valid range")
     

     IF (BOTTOM_ROUGHNESS_KIND .eq. STTC) THEN
        if (BOTTOM_ROUGHNESS_LENGTHSCALE .NE. -1.0_SP)THEN
           CALL WARNING("The name list variable: BOTTOM_ROUGHNESS_LENGTHSCALE",&
                "should have the default value when the bottom roughness is set using 'static'",&
                "(ie set from an input file!)")
           
           BOTTOM_ROUGHNESS_LENGTHSCALE = BOTTOM_ROUGHNESS_LENGTHSCALE
        end if
        
     else if(BOTTOM_ROUGHNESS_KIND .eq. CNSTNT) THEN
        
        if (BOTTOM_ROUGHNESS_LENGTHSCALE .LE. 0.0_SP) CALL FATAL_ERROR&
             &("The name list variable: BOTTOM_ROUGHNESS_LENGTHSCALE",&
             & "Can not have the default value -1.0 when the bottom roughness is set using 'constant'")
     else
        CALL FATAL_ERROR&
             &("The name list variable: BOTTOM_ROUGHNESS_KIND",&
             & "Is set to an invalid option: (choices: constant, static)")     
     end if

     
     
     ! FIX SCALAR_POSITIVITY_CONTROL and HEATING_TYPE SETTINGS
#    if !defined (HEATING_CALCULATED)
     IF (HEATING_ON) THEN

        ! SET THE FVCOM NAMES FOR SURFACE HEATING
        RHEAT = HEATING_LONGWAVE_PERCTAGE
        if (RHEAT .LT. 0) CALL FATAL_ERROR&
             &("Name List Error: HEATING_LONGWAVE_PERCTAGE outside valid range")
        
        ZETA1 = HEATING_LONGWAVE_LENGTHSCALE
        if (ZETA1 .LT. 0) CALL FATAL_ERROR&
             &("Name List Error: HEATING_LONGWAVE_LENGTHSCALE outside valid range")
        
        ZETA2 = HEATING_SHORTWAVE_LENGTHSCALE
        if (ZETA2 .LT. 0) CALL FATAL_ERROR&
             &("Name List Error: HEATING_SHORTWAVE_LENGTHSCALE outside valid range")
        
        
        IF(SCALAR_POSITIVITY_CONTROL .AND. (HEATING_TYPE == 'body')) THEN
           CALL FATAL_ERROR &
                &("YOU CAN NOT USE SCALAR POSITIVITY CONTROL WITH BODY HEATING")
        END IF
     ELSE 
        HEATING_TYPE = 'none'
        RHEAT=0.0_SP
        ZETA1=0.0_SP
        ZETA2=0.0_SP
     END IF
#    endif
     
#    if defined (HEATING_CALCULATED)
     IF (HEATING_CALCULATE_ON) THEN
        
        ! SET THE FVCOM NAMES FOR SURFACE HEATING
        RHEAT = HEATING_LONGWAVE_PERCTAGE_IN_HEATFLUX
        if (RHEAT .LT. 0) CALL FATAL_ERROR&
             &("Name List Error: HEATING_LONGWAVE_PERCTAGE_IN_HEATFLUX outside valid range")
        
        ZETA1 = HEATING_LONGWAVE_LENGTHSCALE_IN_HEATFLUX
        if (ZETA1 .LT. 0) CALL FATAL_ERROR&
             &("Name List Error: HEATING_LONGWAVE_LENGTHSCALE_IN_HEATFLUX outside valid range")
        
        ZETA2 = HEATING_SHORTWAVE_LENGTHSCALE_IN_HEATFLUX
        if (ZETA2 .LT. 0) CALL FATAL_ERROR&
             &("Name List Error: HEATING_SHORTWAVE_LENGTHSCALE_IN_HEATFLUX outside valid range")
        
        
        IF(SCALAR_POSITIVITY_CONTROL .AND. (HEATING_CALCULATE_TYPE == 'body')) THEN
           CALL FATAL_ERROR &
                &("YOU CAN NOT USE SCALAR POSITIVITY CONTROL WITH BODY HEATING")
        END IF
     ELSE 
        HEATING_CALCULATE_TYPE = 'none'
        HEATING_TYPE = 'none'
        RHEAT=0.0_SP
        ZETA1=0.0_SP
        ZETA2=0.0_SP
     END IF
#    endif     
     
     IF( .not. OBC_ON .or. .not. OBC_ELEVATION_FORCING_ON) THEN
        IF (OBC_LONGSHORE_FLOW_ON) THEN
           CALL FATAL_ERROR &
                &("YOU CAN NOT USE THE LONGSHORE FLOW BOUNDRY ADJUSTMENT &
                &WITHOUT OPEN BOUNDARY ELEVATION FORCING!")
        END IF
     END IF
     
     ! SET DEFAULT RUN MODE - DATA ASSIMILATION IS OFF
     FVCOM_RUN_MODE = FVCOM_PURE_SIM
     !==============================================================================!
     !  READ DATA ASSIMILATION NAMELIST AND SET PARAMETERS                          !
     !==============================================================================!
     IF (DATA_ASSIMILATION) THEN
        if(DBG_SET(dbg_log)) WRITE(IPT,*) "! STARTING DATA ASSIMILATION MODE "
#       if defined (DATA_ASSIM)
        CALL SET_ASSIM_PARAM
#       endif
     ELSE
        if(DBG_SET(dbg_log)) WRITE(IPT,*) "! DATA ASSIMILATION MODE IS OFF "
     END IF
     
  END IF


  ! SET DEFAULT TO TRUE FOR REAL TIME MODEL
  use_real_world_time = .TRUE.
  ! TEST FOR IDEALIZED MODEL CASE
  if (timezone == 'none' .or. timezone == "NONE" .or.&
       & timezone == "None") use_real_world_time = .FALSE.


  ! CHECK FOR VALID TIME ZONE
  TEST = IS_VALID_TIMEZONE(timezone)
  IF(.not. TEST) call fatal_error("You selected an invalid time zone: "&
             &//trim(timezone),"Time Zones must be CAPITALS",&
             & "see mod_time.F for a list of valid time_zones")

  call REGISTER_FUNC(DUMP_NC_DAT,NC_CODE,stat)
  IF (stat/=0) CALL FATAL_ERROR("REGISTER_FUNC: FAILED TO REGISTER:: DUMP_NC_DAT")

  call REGISTER_FUNC(DUMP_NC_SF,NCSF_CODE,stat)
  IF (stat/=0) CALL FATAL_ERROR("REGISTER_FUNC: FAILED TO REGISTER:: DUMP_NC_SF")

  call REGISTER_FUNC(DUMP_NC_RST,RESTART_CODE,stat)
  IF (stat/=0) CALL FATAL_ERROR("REGISTER_FUNC: FAILED TO REGISTER:: DUMP_NC_RST")

  call REGISTER_FUNC(DUMP_NC_AVG,NCAV_CODE,stat)
  IF (stat/=0) CALL FATAL_ERROR("REGISTER_FUNC: FAILED TO REGISTER:: DUMP_NC_AVG")

  call REGISTER_FUNC(INIT_NCDIO,INIT_CODE,stat)
  IF (stat/=0) CALL FATAL_ERROR("REGISTER_FUNC: FAILED TO REGISTER:: INIT_NCDIO")


# if defined (MULTIPROCESSOR)
  call REGISTER_FUNC(SETUP_NEST_DOMAIN,INITNEST_CODE,stat)
  IF (stat/=0) CALL FATAL_ERROR("REGISTER_FUNC: FAILED TO REGISTER:: SETUP_NEST_DOMAIN")

  call REGISTER_FUNC(DUMP_NCNEST_FILE,NESTING_CODE,stat)
  IF (stat/=0) CALL FATAL_ERROR("REGISTER_FUNC: FAILED TO REGISTER:: DUMP_NC_NEST")

# if defined (WAVE_CURRENT_INTERACTION)
  call REGISTER_FUNC(DUMP_NCNEST_FILE_WAVE,NESTING_CODE_WAVE,stat)
  IF (stat/=0) CALL FATAL_ERROR("REGISTER_FUNC: FAILED TO REGISTER:: DUMP_NC_NEST_WAVE")
# endif
# endif

  !//////////////////////////////////////////////////////////////////
  !//////////////////////////////////////////////////////////////////
  !==============================================================================!
  !  LOAD THE VISIT LIBRARIES AND DUMP SIM FILE                                  !
  !==============================================================================!
  Call INIT_VISIT


  !//////////////////////////////////////////////////////////////////
  !//////////////////////////////////////////////////////////////////


END SUBROUTINE CNTRL_PRMTRS



