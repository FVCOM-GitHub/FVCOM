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

SUBROUTINE COORDS_N_CONST

  USE MOD_INPUT
  USE MOD_SETUP
  IMPLICIT NONE
  INTEGER STATUS, I

!===============================================================================!
! OPEN FILES
! UNLESS OTHERWISE SPECIFED THE ROUTINES CALLED HERE ARE IN mod_input.F
!===============================================================================!
    
  SELECT CASE(STARTUP_TYPE)
     !=================================================
     ! HOTSTART 
  CASE("hotstart")
     !=================================================
     if(DBG_SET(dbg_log)) then 
        WRITE(IPT,*)'!               READING FILES FOR HOTSTART                       !'
        WRITE(IPT,*)'!                                                                !'
     end if
               
     ! ALLOCATE SPACE FOR THE LOCAL GRID DATA
     ALLOCATE(Y_LCL(0:MT),stat=status)
     IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE Y_LCL")
     Y_LCL=0.0_SP
     ALLOCATE(X_LCL(0:MT),stat=status)
     IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE X_LCL")
     X_LCL=0.0_SP
     ALLOCATE(H_LCL(0:MT),stat=status)
     IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE H_LCL")
     H_LCL=0.0_SP

     CALL LOAD_RESTART_COORDS(X_LCL,Y_LCL)
     CALL COORDINATE_UNITS(X_LCL,Y_LCL)
     CALL SETUP_CENTER_COORDS

#    if defined (SPHERICAL)
     IF(EQUATOR_BETA_PLANE)THEN
      CALL N2E2D(Y_LCL,COR)
      F_ALFA = 1.0_SP-0.8_SP*EXP(-(COR/2.2_SP)**2)
     END IF 
#    endif

     DEALLOCATE(X_LCL)
     DEALLOCATE(Y_LCL)


     CALL LOAD_RESTART_DEPTH(H_LCL)
     CALL SETUP_DEPTH
     DEALLOCATE(H_LCL) ! COULD BE LOADED DIRECTLY - MUST SET MAX/MIN

     CALL LOAD_RESTART_CORIOLIS(COR)! LOAD DIRECTLY

    !!  ggao for equatoral min (4deg)
!     IF(.NOT. EQUATOR_BETA_PLANE)THEN
!      WHERE(COR <  1.e-5_SP .AND. COR > 0.0_SP) COR =  1.e-5_SP
!      WHERE(COR > -1.e-5_SP .AND. COR < 0.0_SP) COR = -1.e-5_SP
!     END IF 

     CALL SETUP_GRAVITY
     
     CALL LOAD_RESTART_SPONGE(CC_SPONGE) ! LOAD DIRECTLY
     
     STYPE = STYPE_RESTART
     CALL LOAD_RESTART_SIGMA(Z,Z1) ! LOAD DIRECTLY TO ALL_VARS:Z,Z1
     CALL SETUP_SIGMA_DERIVATIVES

     CALL SETUP_OBCTYPES

     CALL SETUP_HORIZONTAL_MIXING_COEFFICIENT

     CALL SETUP_BOTTOM_ROUGHNESS

     !=================================================
     ! CRASHSTART 
  CASE("crashrestart")
     !=================================================
     if(DBG_SET(dbg_log)) then 
        WRITE(IPT,*)'!              READING FILES FOR CRASHRESTART                    !'
        WRITE(IPT,*)'!                                                                !'
     end if
     ! ALLOCATE SPACE FOR THE LOCAL GRID DATA
     ALLOCATE(Y_LCL(0:MT),stat=status)
     IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE Y_LCL")
     Y_LCL=0.0_SP
     ALLOCATE(X_LCL(0:MT),stat=status)
     IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE X_LCL")
     X_LCL=0.0_SP
     ALLOCATE(H_LCL(0:MT),stat=status)
     IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE H_LCL")
     H_LCL=0.0_SP

          
     CALL LOAD_RESTART_COORDS(X_LCL,Y_LCL)
     CALL COORDINATE_UNITS(X_LCL,Y_LCL)
     CALL SETUP_CENTER_COORDS

#    if defined (SPHERICAL)
     IF(EQUATOR_BETA_PLANE)THEN
      CALL N2E2D(Y_LCL,COR)
      F_ALFA = 1.0_SP-0.8_SP*EXP(-(COR/2.2_SP)**2)
     END IF 
#    endif

     DEALLOCATE(X_LCL)
     DEALLOCATE(Y_LCL)

     
     CALL LOAD_RESTART_DEPTH(H_LCL)
     CALL SETUP_DEPTH
     DEALLOCATE(H_LCL) ! COULD BE LOADED DIRECTLY - MUST SET MAX/MIN

     CALL LOAD_RESTART_CORIOLIS(COR)

     CALL SETUP_GRAVITY
     
     CALL LOAD_RESTART_SPONGE(CC_SPONGE) ! LOAD DIRECTLY TO ALL_VARS:CC_SPONGE


     STYPE = STYPE_RESTART
     CALL LOAD_RESTART_SIGMA(Z,Z1) ! LOAD DIRECTLY TO ALL_VARS:Z,Z1
     CALL SETUP_SIGMA_DERIVATIVES

     CALL SETUP_OBCTYPES

     CALL SETUP_HORIZONTAL_MIXING_COEFFICIENT

     CALL SETUP_BOTTOM_ROUGHNESS

     !=================================================
     ! COLDSTART 
  CASE("coldstart")
     !=================================================
     if(DBG_SET(dbg_log)) then 
        WRITE(IPT,*)'!              READING FILES FOR COLDSTART                       !'
        WRITE(IPT,*)'!                                                                !'
     end if
     
     IF(MSR) THEN
        ! ALLOCATE SPACE FOR THE GLOBAL GRID DATA
        ALLOCATE(Y_GBL(0:MGL),stat=status)
        IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE Y_GBL")
        Y_GBL=0.0_SP
        ALLOCATE(X_GBL(0:MGL),stat=status)
        IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE X_GBL")
        X_GBL=0.0_SP
     END IF

    ! ALLOCATE SPACE FOR THE LOCAL GRID DATA
    ALLOCATE(Y_LCL(0:MT),stat=status)
    IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE Y_LCL")
    Y_LCL=0.0_SP
    ALLOCATE(X_LCL(0:MT),stat=status)
    IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE X_LCL")
    X_LCL=0.0_SP

    ALLOCATE(H_LCL(0:MT),stat=status)
    IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE H_LCL")
    H_LCL=0.0_SP

    ALLOCATE(C_LCL(0:MT),stat=status)
    IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE C_LCL")
    C_LCL=0.0_SP

     
     CALL LOAD_COLDSTART_COORDS(X_GBL,Y_GBL,X_LCL,Y_LCL)
     CALL COORDINATE_UNITS(X_LCL,Y_LCL)
     CALL SETUP_CENTER_COORDS


     DEALLOCATE(X_LCL)
     DEALLOCATE(Y_LCL)

     CALL LOAD_COLDSTART_DEPTH(X_GBL,Y_GBL,H_LCL)
     CALL SETUP_DEPTH

     DEALLOCATE(H_LCL)

     CALL LOAD_COLDSTART_CORIOLIS(X_GBL,Y_GBL,C_LCL)
     CALL SETUP_CORIOLIS

     DEALLOCATE(C_LCL)

     CALL SETUP_GRAVITY

     ! THESE VARIABLES ARE NOT ALLOCATED YET
     CALL LOAD_COLDSTART_SPONGE(X_GBL,Y_GBL,NSPONGE,N_SPG,R_SPG,C_SPG,X_SPG,Y_SPG)
     CALL SETUP_SPONGE

 
!     CALL LOAD_COLDSTART_SIGMA
     CALL SETUP_SIGMA
     CALL SETUP_SIGMA_DERIVATIVES

     IF(MSR) THEN
        DEALLOCATE(X_GBL)
        DEALLOCATE(Y_GBL)
     END IF

     CALL SETUP_OBCTYPES

     CALL SETUP_HORIZONTAL_MIXING_COEFFICIENT

     CALL SETUP_BOTTOM_ROUGHNESS

  END SELECT
  
  
END SUBROUTINE COORDS_N_CONST
