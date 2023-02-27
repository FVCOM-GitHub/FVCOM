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
!  ADJUST THE VERTICAL WATER COLUMN WHEN UNSTABLE                              |
!==============================================================================|

   SUBROUTINE CONV_OVER           

!==============================================================================|
   USE ALL_VARS
   USE MOD_UTILS
   IMPLICIT NONE
   REAL(SP):: AVE_T,AVE_S,AVE_R 
   INTEGER :: I,K,KK,J1,J2,J3
!==============================================================================|

   IF (DBG_SET(DBG_SBR)) WRITE(IPT,*) "Start: conv_over"
 

!--APPROXIMATE CONVECTIVE OVERTURNING------------------------------------------!

   DO I=1,M
     DO K=KBM1,2,-1
       DO KK=K-1,1,-1
         IF(RHO1(I,K) < RHO1(I,KK)) THEN
           AVE_T = SUM(  T1(I,KK:K))/FLOAT(K-KK+1) 
           AVE_S = SUM(  S1(I,KK:K))/FLOAT(K-KK+1) 
           AVE_R = SUM(RHO1(I,KK:K))/FLOAT(K-KK+1) 
           T1(I,KK:K)   = AVE_T
           S1(I,KK:K)   = AVE_S
           RHO1(I,KK:K) = AVE_R
         END IF
       END DO
     END DO
   END DO
       
!-----RECALCULATE ELEMENT-BASED VALUES OF SALINITY/TEMP/DENSITY----------------!

   DO I=1,N
     J1=NV(I,1) ; J2 = NV(I,2) ; J3 = NV(I,3)
     DO K=1,KBM1
       T(I,K)  = ONE_THIRD*(  T1(J1,K)+  T1(J2,K)+  T1(J3,K))
       S(I,K)  = ONE_THIRD*(  S1(J1,K)+  S1(J2,K)+  S1(J3,K))
     END DO
   END DO   
          
   IF (DBG_SET(DBG_SBR)) WRITE(IPT,*) "End: conv_over"
       
   RETURN
   END SUBROUTINE CONV_OVER       
!==============================================================================|

