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
!   Calculate Advection and Horizontal Diffusion Terms for Temperature         |
!==============================================================================|

   SUBROUTINE VISCOF_H               

!------------------------------------------------------------------------------|
   USE MOD_UTILS
   USE ALL_VARS

   IMPLICIT NONE
   REAL(SP) :: PUPX,PUPY,PVPX,PVPY
   REAL(SP) :: tmp1,tmp2
   INTEGER  :: I,I1,K,J


   IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "Start: viscofh"

   SELECT CASE(HORIZONTAL_MIXING_TYPE)
   CASE ('closure')
      ! Run Subroutine
   CASE('constant')
      IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "End: viscofh (constant)"
   CASE DEFAULT
      CALL FATAL_ERROR("UNKNOW HORIZONTAL MIXING TYPE:",&
           & TRIM(HORIZONTAL_MIXING_TYPE) )
   END SELECT


   DO K=1,KBM1
     DO I=1,M
       PUPX=0.0_SP
       PUPY=0.0_SP
       PVPX=0.0_SP
       PVPY=0.0_SP

       J=1
       I1=NBVE(I,J)


       PUPX=PUPX+U(I1,K)*DLTYECEC(I,J)
       PUPY=PUPY+U(I1,K)*DLTXECEC(I,J)
       PVPX=PVPX+V(I1,K)*DLTYECEC(I,J)
       PVPY=PVPY+V(I1,K)*DLTXECEC(I,J)
       

       IF(ISONB(I) /= 0) THEN

         PUPX=PUPX+U(I1,K)*DLTYNEC(I,J)
         PUPY=PUPY+U(I1,K)*DLTXNEC(I,J)
         PVPX=PVPX+V(I1,K)*DLTYNEC(I,J)
         PVPY=PVPY+V(I1,K)*DLTXNEC(I,J)

       END IF

       DO J=2,NTVE(I)-1
         I1=NBVE(I,J)
         
         PUPX=PUPX+U(I1,K)*DLTYECEC(I,J)
         PUPY=PUPY+U(I1,K)*DLTXECEC(I,J)
         PVPX=PVPX+V(I1,K)*DLTYECEC(I,J)
         PVPY=PVPY+V(I1,K)*DLTXECEC(I,J)

       END DO

       J=NTVE(I)
       I1=NBVE(I,J)


       PUPX=PUPX+U(I1,K)*DLTYECEC(I,J)
       PUPY=PUPY+U(I1,K)*DLTXECEC(I,J)
       PVPX=PVPX+V(I1,K)*DLTYECEC(I,J)
       PVPY=PVPY+V(I1,K)*DLTXECEC(I,J)
       

       IF(ISONB(I) /= 0) THEN

          PUPX=PUPX+U(I1,K)*(-DLTYNEC(I,J))
          PUPY=PUPY+U(I1,K)*(-DLTXNEC(I,J))
          PVPX=PVPX+V(I1,K)*(-DLTYNEC(I,J))
          PVPY=PVPY+V(I1,K)*(-DLTXNEC(I,J))

       END IF

       PUPX=PUPX/ART1(I)
       PUPY=PUPY/ART1(I)
       PVPX=PVPX/ART1(I)
       PVPY=PVPY/ART1(I)
       TMP1=PUPX**2+PVPY**2
       TMP2=0.5_SP*(PUPY+PVPX)**2
       VISCOFH(I,K)=SQRT(TMP1+TMP2)*ART1(I)
       
     END DO
   END DO  
    

   IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "End: viscofh"

   END SUBROUTINE VISCOF_H
!==============================================================================|
