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
! this subroutine uses an implicit method to calculate the vertical            !
! diffusion term in the x and y momentum equations                             !
!==============================================================================|

   SUBROUTINE VDIF_UV             

!------------------------------------------------------------------------------|

   USE ALL_VARS
   USE MOD_UTILS
#  if defined (WET_DRY)
   USE MOD_WD
#  endif
#  if defined (SEMI_IMPLICIT)
   USE MOD_SEMI_IMPLICIT, only: STAGE, KSTAGE_UV, UC, VC
#  endif
#  if defined (THIN_DAM)
   USE MOD_DAM
#  endif
   IMPLICIT NONE
   INTEGER I,K,ITMP1,ITMP2,ITMP3,KI
   REAL(SP), DIMENSION(0:NT,KB) :: C,A,VHU,VHPU,VHV,VHPV
   REAL(SP) :: TMP

   IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "Start: vdif_uv"

   C = KM1

   DO K = 2, KBM1
     DO I = 1, N
#  if !defined (WET_DRY)
       IF (D1(I) > 0.0_SP) THEN
#  else
       IF(ISWETC(I) == 1)THEN
#  endif

         A(I,K-1) = -DTI *(C(I,K)+UMOL) / (DZ1(I,K-1)*DZZ1(I,K-1)* D1(I)*D1(I))
         C(I,K) = -DTI *(C(I,K)+UMOL) / (DZ1(I,K)*DZZ1(I,K-1)* D1(I)*D1(I))
       END IF
     END DO
   END DO

   DO I = 1, N
#  if !defined (WET_DRY)
     IF (D1(I) > 0.0_SP) THEN
#  else
     IF(ISWETC(I) == 1)THEN
#  endif
       ITMP1=ISONB(NV(I,1))
       ITMP2=ISONB(NV(I,2))
       ITMP3=ISONB(NV(I,3))
       IF(ITMP1 == 2 .OR. ITMP2 == 2 .OR. ITMP3 == 2) THEN
         TMP=0.0_SP
       ELSE
         TMP=1.0_SP
       END IF
       VHU(I,1) = A(I,1) / (A(I,1)-1.)
       VHV(I,1) = A(I,1) / (A(I,1)-1.)
       VHPU(I,1) = (-DTI*WUSURF(I)*TMP/(-DZ1(I,1)*D1(I))-UF(I,1)) / (A(I,1)-1.)
       VHPV(I,1) = (-DTI*WVSURF(I)*TMP/(-DZ1(I,1)*D1(I))-VF(I,1)) / (A(I,1)-1.)
     END IF
   END DO

   DO  K = 2, KBM2
     DO I = 1, N
#  if !defined (WET_DRY)
       IF (D1(I) > 0.0_SP) THEN
#  else
       IF(ISWETC(I) == 1)THEN
#  endif
         VHPU(I,K) = 1.0_SP / (A(I,K)+C(I,K)*(1.-VHU(I,K-1))-1.)
         VHPV(I,K) = 1.0_SP / (A(I,K)+C(I,K)*(1.-VHV(I,K-1))-1.)
         VHU(I,K)  = A(I,K) * VHPU(I,K)
         VHV(I,K)  = A(I,K) * VHPV(I,K)
         VHPU(I,K) = (C(I,K)*VHPU(I,K-1)-UF(I,K))*VHPU(I,K)
         VHPV(I,K) = (C(I,K)*VHPV(I,K-1)-VF(I,K))*VHPV(I,K)
       END IF
     END DO
   END DO

   DO  I = 1, N
#  if !defined (WET_DRY)
     IF (D1(I) > 0.0_SP) THEN
#  else
     IF(ISWETC(I) == 1)THEN
#  endif
#      if !defined (SEMI_IMPLICIT)
       TPS(I) = CBC(I) * SQRT(U(I,KBM1)**2+V(I,KBM1)**2)
#      else
       TPS(I) = CBC(I) * SQRT(UC(I,KBM1)**2+VC(I,KBM1)**2)
#      endif
       UF(I,KBM1) = (C(I,KBM1)*VHPU(I,KBM2)-UF(I,KBM1))/ &
       (TPS(I)*DTI/(-DZ1(I,KBM1)*D1(I))-1.-(VHU(I,KBM2)-1.)*C(I,KBM1))
       VF(I,KBM1) = (C(I,KBM1)*VHPV(I,KBM2)-VF(I,KBM1))/ &
       (TPS(I)*DTI/(-DZ1(I,KBM1)*D1(I))-1.-(VHV(I,KBM2)-1.)*C(I,KBM1))
     END IF
   END DO

   DO  K = 2, KBM1
     KI = KB - K
     DO  I = 1, N
#  if !defined (WET_DRY)
       IF (D1(I) > 0.0_SP) THEN
#  else
       IF(ISWETC(I) == 1)THEN
#  endif
         UF(I,KI) = (VHU(I,KI)*UF(I,KI+1)+VHPU(I,KI))
         VF(I,KI) = (VHV(I,KI)*VF(I,KI+1)+VHPV(I,KI))
#  if defined (WET_DRY)
       ELSE
         UF(I,KI) = 0.0_SP
         VF(I,KI) = 0.0_SP
#  endif
       END IF
     END DO
   END DO

!
!--Damp Velocity In Sponge Region----------------------------------------------!
!
   DO K = 1, KBM1
     DO I = 1, N
!#  if defined (THIN_DAM)
!       UF(I,K) = UF(I,K)-DAM_SPONGE(I)*UF(I,K)
!       VF(I,K) = VF(I,K)-DAM_SPONGE(I)*VF(I,K)
!#  endif
!old:       UF(I,K) = UF(I,K)-CC_SPONGE(I)*UF(I,K)
!old:       VF(I,K) = VF(I,K)-CC_SPONGE(I)*VF(I,K)
! ---- new: Karsten Lettmann: 2012.06.25 -------
        UF(I,K) = UF(I,K)/(1.0_SP+CC_SPONGE(I)*UF(I,K)**2.0_SP)
        VF(I,K) = VF(I,K)/(1.0_SP+CC_SPONGE(I)*VF(I,K)**2.0_SP)
! ------- end new -------------------------------
     END DO
   END DO

#  if defined (SEMI_IMPLICIT)
   IF(STAGE==KSTAGE_UV) THEN
#  endif
   DO  I = 1, N
#  if !defined (WET_DRY)
     IF (D1(I) > 0.0_SP) THEN
#  else
     IF(ISWETC(I) == 1)THEN
#  endif
       WUBOT(I) = -TPS(I) * U(I,KBM1)
       WVBOT(I) = -TPS(I) * V(I,KBM1)
#  if defined (WET_DRY)
     ELSE
       WUBOT(I) = 0.0_SP
       WVBOT(I) = 0.0_SP
#  endif
     END IF
   END DO
#  if defined (SEMI_IMPLICIT)
   ENDIF
#  endif

   IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "Start: vdif_uv"

   END SUBROUTINE VDIF_UV
!==============================================================================|
