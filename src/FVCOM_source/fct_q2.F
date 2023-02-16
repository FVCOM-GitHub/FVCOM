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
!  FLUX CONTROL FOR SALINITY                                                        |
!==============================================================================|

   SUBROUTINE FCT_Q2
!#  if defined (WET_DRY)

!==============================================================================|
   USE ALL_VARS
   USE MOD_UTILS
   USE BCS
   USE MOD_OBCS
   IMPLICIT NONE
   REAL(SP):: Q2MAX,Q2MIN
   INTEGER :: I,J,K
!==============================================================================|

   IF(HEATING_TYPE == 'body') return
   
   IF(DBG_SET(DBG_SBR)) write(ipt,*) "Start: fct_q2"


nodes:DO I=1,M

     IF(IOBCN > 0)THEN
       DO J=1,IOBCN
         IF(I == I_OBC_N(J)) CYCLE nodes
       END DO
     END IF  	
 
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
     DO K=2,KBM1
       Q2MAX = MAXVAL(Q2(NBSN(I,1:NTSN(I)),K))
       Q2MIN = MINVAL(Q2(NBSN(I,1:NTSN(I)),K))

       IF(K == 2)THEN
         Q2MAX = MAX(Q2MAX,(Q2(I,K)*DZZ(I,K)+Q2(I,K+1)*DZZ(I,K-1))/  &
                 (DZZ(I,K)+DZZ(I,K-1)))
         Q2MIN = MIN(Q2MIN,(Q2(I,K)*DZZ(I,K)+Q2(I,K+1)*DZZ(I,K-1))/  &
                 (DZZ(I,K)+DZZ(I,K-1)))
       ELSE IF(K == KBM1)THEN
         Q2MAX = MAX(Q2MAX,(Q2(I,K)*DZZ(I,K-2)+Q2(I,K-1)*DZZ(I,K-1))/  &
                 (DZZ(I,K-1)+DZZ(I,K-2)))
         Q2MIN = MIN(Q2MIN,(Q2(I,K)*DZZ(I,K-2)+Q2(I,K-1)*DZZ(I,K-1))/  &
                 (DZZ(I,K-1)+DZZ(I,K-2)))
       ELSE
         Q2MAX = MAX(Q2MAX,(Q2(I,K)*DZZ(I,K-2)+Q2(I,K-1)*DZZ(I,K-1))/  &
                 (DZZ(I,K-1)+DZZ(I,K-2)),                              &
                 (Q2(I,K)*DZZ(I,K)+Q2(I,K+1)*DZZ(I,K-1))/              &
                 (DZZ(I,K)+DZZ(I,K-1)))  
         Q2MIN = MIN(Q2MIN,(Q2(I,K)*DZZ(I,K-2)+Q2(I,K-1)*DZZ(I,K-1))/  &
                 (DZZ(I,K-1)+DZZ(I,K-2)),                              &
                 (Q2(I,K)*DZZ(I,K)+Q2(I,K+1)*DZZ(I,K-1))/              &
                 (DZZ(I,K)+DZZ(I,K-1)))
       END IF

       IF(Q2MIN-Q2F(I,K) > 0.0_SP)Q2F(I,K) = Q2MIN
       IF(Q2F(I,K)-Q2MAX > 0.0_SP)Q2F(I,K) = Q2MAX

     END DO

   END DO nodes

   IF(DBG_SET(DBG_SBR)) write(ipt,*) "End: fct_q2"

   END SUBROUTINE FCT_Q2
!==============================================================================|


