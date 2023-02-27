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

   SUBROUTINE ADJUST2D3D(ADJUST_TYPE)
# if !defined (SEMI_IMPLICIT)
!==============================================================================|
!    ADJUST 3D VELOCITY USING DEFECT BETWEEN UPDATED AND CURRENT VERTICALLY    !
!    AVERAGED VELOCITIES						       !
! 									       !
!    FORMULA IS:							       !
!									       !
!      U_adjusted = U_orig + eps*(U_avg_new - U_avg_current)		       !
!      eps = 0 : no adjustment						       !
!      eps = 1 : full adjustment					       !
!==============================================================================|
#  if defined (WET_DRY)
   USE MOD_WD
#  endif
   USE ALL_VARS
   USE MOD_UTILS
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: ADJUST_TYPE
   INTEGER :: I,K
   REAL(SP), PARAMETER :: EPS = 1.0_SP
   REAL(SP) :: UAC,VAC,UTMP,VTMP
!==============================================================================!

   if(dbg_set(dbg_sbr)) write(ipt,*) "Start: adjust2d3d"

   SELECT CASE(ADJUST_TYPE)

   CASE(1)
   DO I=1,NT
     UAC    = SUM(U(I,1:KBM1)*DZ1(I,1:KBM1))
     VAC    = SUM(V(I,1:KBM1)*DZ1(I,1:KBM1))
     U(I,1:KBM1) = U(I,1:KBM1) + EPS*(UA(I) - UAC) 
     V(I,1:KBM1) = V(I,1:KBM1) + EPS*(VA(I) - VAC) 
   END DO

   CASE(2)
   UARD = UARD/FLOAT(ISPLIT)
   VARD = VARD/FLOAT(ISPLIT)
!!#  if defined (WET_DRY)
!!   UARDS = UARDS/FLOAT(ISPLIT)
!!   VARDS = VARDS/FLOAT(ISPLIT)
!!#  endif

                                                                                                                         
   DO I=1,NT
#    if defined (WET_DRY)
     IF(ISWETCT(I)*ISWETC(I) == 1)THEN
#    endif
       UTMP = 0.0_SP ; VTMP = 0.0_SP
       DO K=1,KBM1
         UTMP = UTMP + U(I,K)*DZ1(I,K)
         VTMP = VTMP + V(I,K)*DZ1(I,K)
       END DO
       UTMP = UTMP*DT1(I)
       VTMP = VTMP*DT1(I)
       DO K=1,KBM1
         U(I,K) = U(I,K) - (UTMP-UARD(I))/DT1(I)
         V(I,K) = V(I,K) - (VTMP-VARD(I))/DT1(I)
       END DO
#    if defined (WET_DRY)
     END IF
#    endif
   END DO

!!#  if defined (WET_DRY)
!!   DO I=1,NT
!!     UTMP = 0.0_SP ; VTMP = 0.0_SP
!!     DO K=1,KBM1
!!       UTMP = UTMP + U(I,K)*DZ1(I,K)
!!       VTMP = VTMP + V(I,K)*DZ1(I,K)
!!     END DO
!!     UTMP = UTMP*DT1(I)
!!     VTMP = VTMP*DT1(I)
!!     DO K=1,KBM1
!!       US(I,K) = U(I,K) - (UTMP-UARDS(I))/DT1(I)
!!       VS(I,K) = V(I,K) - (VTMP-VARDS(I))/DT1(I)
!!     END DO
!!   END DO
!!#  endif

   END SELECT

   if(dbg_set(dbg_sbr)) write(ipt,*) "End: adjust2d3d"

#  endif
   END SUBROUTINE ADJUST2D3D
!==============================================================================|
