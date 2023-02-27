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

!==============================================================================|
!  FLUX CONTROL FOR TEMPERATURE                                                     |
!==============================================================================|

SUBROUTINE FCT_T
  !#  if defined (WET_DRY)

  !==============================================================================|
  USE ALL_VARS
  USE MOD_UTILS
  USE BCS
  USE MOD_OBCS
  IMPLICIT NONE
  REAL(SP):: TMAX,TMIN
  INTEGER :: I,J,K
  !==============================================================================|


  IF(HEATING_TYPE == 'body') RETURN

  IF(DBG_SET(DBG_SBR)) WRITE(IPT,*)"Start: fct_t"

  nodes: DO I=1,M
     ! SKIP OPEN BOUNDARY NODES
     IF(IOBCN > 0)THEN
        DO J=1,IOBCN
           IF(I == I_OBC_N(J)) CYCLE nodes
        END DO
     END IF

     ! SKIP RIVER INFLOW POINTS
     IF(NUMQBC > 0)THEN
        DO J=1,NUMQBC
           IF(RIVER_INFLOW_LOCATION == 'node')THEN
              IF(I == INODEQ(J)) CYCLE nodes
           END IF
           IF(RIVER_INFLOW_LOCATION == 'edge')THEN
              IF(I == N_ICELLQ(J,1) .OR. I == N_ICELLQ(J,2)) CYCLE nodes
           END IF
        END DO
     END IF

     ! SKIP GROUND WATER INFLOW POINTS
     IF(BFWDIS(I) .GT. 0.0_SP .and. GROUNDWATER_TEMP_ON) CYCLE nodes

     DO K=1,KBM1
        TMAX = MAXVAL(T1(NBSN(I,1:NTSN(I)),K))
        TMIN = MINVAL(T1(NBSN(I,1:NTSN(I)),K))

        IF(K == 1)THEN
           TMAX = MAX(TMAX,(T1(I,K)*DZ(I,K+1)+T1(I,K+1)*DZ(I,K))/  &
                (DZ(I,K)+DZ(I,K+1)))
           TMIN = MIN(TMIN,(T1(I,K)*DZ(I,K+1)+T1(I,K+1)*DZ(I,K))/  &
                (DZ(I,K)+DZ(I,K+1)))
        ELSE IF(K == KBM1)THEN
           TMAX = MAX(TMAX,(T1(I,K)*DZ(I,K-1)+T1(I,K-1)*DZ(I,K))/  &
                (DZ(I,K)+DZ(I,K-1)))
           TMIN = MIN(TMIN,(T1(I,K)*DZ(I,K-1)+T1(I,K-1)*DZ(I,K))/  &
                (DZ(I,K)+DZ(I,K-1)))
        ELSE
           TMAX = MAX(TMAX,(T1(I,K)*DZ(I,K-1)+T1(I,K-1)*DZ(I,K))/  &
                (DZ(I,K)+DZ(I,K-1)),                             &
                (T1(I,K)*DZ(I,K+1)+T1(I,K+1)*DZ(I,K))/           &
                (DZ(I,K)+DZ(I,K+1)))
           TMIN = MIN(TMIN,(T1(I,K)*DZ(I,K-1)+T1(I,K-1)*DZ(I,K))/  &
                (DZ(I,K)+DZ(I,K-1)),                             &
                (T1(I,K)*DZ(I,K+1)+T1(I,K+1)*DZ(I,K))/           &
                (DZ(I,K)+DZ(I,K+1)))
        END IF

        IF(TMIN-TF1(I,K) > 0.0_SP)TF1(I,K) = TMIN
        IF(TF1(I,K)-TMAX > 0.0_SP)TF1(I,K) = TMAX

     END DO

  END DO nodes

  IF(DBG_SET(DBG_SBR)) WRITE(IPT,*)"End: fct_t"
  !#  endif
END SUBROUTINE FCT_T
!==============================================================================|


