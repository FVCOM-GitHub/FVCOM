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
!     ADJUST TEMPERATURE NEAR RIVER MOUTHS USING ADJACENT NODES                |
!     ADJUST SALINITY AT RIVER MOUTHS
!==============================================================================|

   SUBROUTINE ADJUST_TS

!==============================================================================|
   USE ALL_VARS
   USE MOD_UTILS
   use mod_par
   USE BCS
   IMPLICIT NONE
   REAL(SP) :: TAVE,TAVE1,TAVE2
   INTEGER :: I,K,JJ,I1,J,J1,J2,NUM_TAVE,NUM_TAVE1,NUM_TAVE2
!==============================================================================|

   IF (DBG_SET(DBG_SBR)) WRITE(IPT,*) "Start: adjust_ts"
   
   IF(NUMQBC > 0)THEN   

!     IF(RIVER_INFLOW_LOCATION == 'node')THEN
!       DO K=1,KBM1
!         DO I=1,NUMQBC
!           JJ=INODEQ(I)
!           TAVE = 0.0_SP
!           NUM_TAVE = 0
!           DO J=2,NTSN(JJ)-1
!             I1=NBSN(JJ,J)
!	     IF(NUMQBC == 1)THEN
!	       NUM_TAVE = NUM_TAVE + 1
!	       TAVE = TAVE + T1(I1,K)
!	     ELSE
!	       IF(I == 1)THEN
!	         IF(I1 /= INODEQ(I+1))THEN 
!	         NUM_TAVE = NUM_TAVE + 1
!	         TAVE = TAVE + T1(I1,K)
!		 END IF
!	       ELSE IF(I == NUMQBC)THEN
!	         IF(I1 /= INODEQ(I-1))THEN	 
!	         NUM_TAVE = NUM_TAVE + 1
!	         TAVE = TAVE + T1(I1,K)
!                 END IF
!	       ELSE IF(I1 /= INODEQ(I-1) .AND. I1 /= INODEQ(I+1))THEN
!	         NUM_TAVE = NUM_TAVE + 1
!                 TAVE = TAVE + T1(I1,K)
!	       END IF	 
!	     END IF
!           END DO
!           T1(JJ,K) = TAVE/FLOAT(NUM_TAVE)
!         END DO
!       END DO
!     ELSE IF(RIVER_INFLOW_LOCATION == 'edge')THEN
!       DO K=1,KBM1
!         DO I=1,NUMQBC
!           J1=N_ICELLQ(I,1)
!           J2=N_ICELLQ(I,2)
!	   TAVE1 = 0.0_SP
!	   TAVE2 = 0.0_SP
!	   NUM_TAVE1 = 0
!	   NUM_TAVE2 = 0

!           DO J=2,NTSN(J1)-1
!             I1=NBSN(J1,J)
!	     IF(NUMQBC == 1)THEN
!	       IF(I1 /= J2)THEN
!	         NUM_TAVE1 = NUM_TAVE1 + 1
!                 TAVE1 = TAVE1 + T1(I1,K)
!	       END IF
!	     ELSE IF(I == 1)THEN
!	       IF(I1 /= J2 .AND. I1 /= N_ICELLQ(I+1,1) .AND. &
!	         I1 /= N_ICELLQ(I+1,2))THEN
!	         NUM_TAVE1 = NUM_TAVE1 + 1
!                 TAVE1 = TAVE1 + T1(I1,K)
!	       END IF
!	     ELSE IF(I == NUMQBC)THEN
!	       IF(I1 /= J2 .AND. I1 /= N_ICELLQ(I-1,1) .AND. &
!	         I1 /= N_ICELLQ(I-1,2))THEN
!	         NUM_TAVE1 = NUM_TAVE1 + 1
!                 TAVE1 = TAVE1 + T1(I1,K)
!	       END IF
!             ELSE IF(I1 /= J2 .AND. &
!	        I1 /= N_ICELLQ(I-1,1) .AND. I1 /= N_ICELLQ(I-1,2) .AND.  &
!		I1 /= N_ICELLQ(I+1,1) .AND. I1 /= N_ICELLQ(I+1,2))THEN
!	       NUM_TAVE1 = NUM_TAVE1 + 1
!               TAVE1 = TAVE1 + T1(I1,K)
!	     END IF
!           END DO
!           T1(J1,K) = TAVE1/FLOAT(NUM_TAVE1)

!           DO J=2,NTSN(J2)-1
!             I1=NBSN(J2,J)
!	     IF(NUMQBC == 1)THEN
!	       IF(I1 /= J1)THEN
!	         NUM_TAVE2 = NUM_TAVE2 + 1
!                 TAVE2 = TAVE2 + T1(I1,K)
!	       END IF
!	     ELSE IF(I == 1)THEN
!	       IF(I1 /= J1 .AND. I1 /= N_ICELLQ(I+1,1) .AND. &
!	         I1 /= N_ICELLQ(I+1,2))THEN
!	         NUM_TAVE2 = NUM_TAVE2 + 1
!                 TAVE2 = TAVE2 + T1(I1,K)
!	       END IF
!	     ELSE IF(I == NUMQBC)THEN  
!	       IF(I1 /= J1 .AND. I1 /= N_ICELLQ(I-1,1) .AND. &
!	         I1 /= N_ICELLQ(I-1,2))THEN
!	         NUM_TAVE2 = NUM_TAVE2 + 1
!                 TAVE2 = TAVE2 + T1(I1,K)
!	       END IF
!	     ELSE IF(I1 /= J1 .AND. &
!	        I1 /= N_ICELLQ(I-1,1) .AND. I1 /= N_ICELLQ(I-1,2) .AND.  &
!		I1 /= N_ICELLQ(I+1,1) .AND. I1 /= N_ICELLQ(I+1,2))THEN
!	       NUM_TAVE2 = NUM_TAVE2 + 1
!               TAVE2 = TAVE2 + T1(I1,K)
!	     END IF
!           END DO
!           T1(J2,K) = TAVE2/FLOAT(NUM_TAVE2)
!	   
!         END DO
!       END DO	 
!     END IF
   
     DO I=1,NUMQBC
       IF(RIVER_INFLOW_LOCATION == 'node')THEN
         J = INODEQ(I)
         DO K=1,KBM1
           T1(J,K) = MAX(T1(J,K),TDIS(I))
           S1(J,K) = MAX(S1(J,K),SDIS(I))
        END DO
       ELSE IF(RIVER_INFLOW_LOCATION == 'edge')THEN
         J1 = N_ICELLQ(I,1)
         J2 = N_ICELLQ(I,2)
         DO K=1,KBM1
           T1(J1,K) = MAX(T1(J1,K),TDIS(I))
           T1(J2,K) = MAX(T1(J2,K),TDIS(I))
           S1(J1,K) = MAX(S1(J1,K),SDIS(I))
           S1(J2,K) = MAX(S1(J2,K),SDIS(I))
         END DO
       END IF
     END DO

     CALL N2E3D(T1,T)
     CALL N2E3D(S1,S)

   END IF

  IF (DBG_SET(DBG_SBR)) WRITE(IPT,*) "End: adjust_ts"
 

   RETURN
   END SUBROUTINE ADJUST_TS
!==============================================================================|

