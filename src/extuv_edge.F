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
!     ACCUMLATE FLUXES FOR EXTERNAL MODE                                       |
!==============================================================================|

   SUBROUTINE EXTUV_EDGE(K)       
# if !defined (SEMI_IMPLICIT)

!==============================================================================|
   USE ALL_VARS
   USE MOD_UTILS
   USE MOD_WD

   USE MOD_NORTHPOLE

#  if defined (BALANCE_2D)
   USE MOD_BALANCE_2D
#  endif

#  if defined (NH)
   USE NON_HYDRO, ONLY: NHQ2DX, NHQ2DY
#  endif

#  if defined (WAVE_CURRENT_INTERACTION)
   USE MOD_WAVE_CURRENT_INTERACTION
#  endif

#  if defined (VEGETATION)
   USE MOD_VEGETATION
#  endif

   IMPLICIT NONE
   INTEGER, INTENT(IN) :: K
   REAL(SP), DIMENSION(0:NT) :: RESX,RESY,TMP
   REAL(SP) :: UAFT,VAFT
   INTEGER  :: I

!==============================================================================|

   if(dbg_set(dbg_sbr)) write(ipt,*) "Start: extuv_edge.F"

!
!--ACCUMULATE RESIDUALS FOR EXTERNAL MODE EQUATIONS----------------------------|
!
   UAFT = UAF(0)
   VAFT = VAF(0)

   ! THIS APPEARS TO BE TO PREVENT DIVISION BY ZERO, BUT IT IS A
   ! STRANGE WAY TO DO IT!
   H1(0)= H1(1)

!!#  if defined (WET_DRY)
!!   IF(K == 3)THEN

!!#    if !defined (NH)
!!     RESX = ADX2D + ADVUA + DRX2D + PSTX - COR*VA*D1*ART  &
!!            -(WUSURF2 + WUBOT)*ART
!!     RESY = ADY2D + ADVVA + DRY2D + PSTY + COR*UA*D1*ART  &
!!            -(WVSURF2 + WVBOT)*ART
!!#    else
!!     RESX = ADX2D + ADVUA + DRX2D + PSTX - COR*VA*D1*ART  &
!!            -(WUSURF2 + WUBOT)*ART + NHQ2DX
!!     RESY = ADY2D + ADVVA + DRY2D + PSTY + COR*UA*D1*ART  &
!!            -(WVSURF2 + WVBOT)*ART + NHQ2DY
!!#    endif

!!#  if defined (SPHERICAL)
!!     RESX = RESX -UA*VA/REARTH*TAN(DEG2RAD*YC)*D1*ART
!!     RESY = RESY +UA*UA/REARTH*TAN(DEG2RAD*YC)*D1*ART
!!#  endif

!!!
!!!--UPDATE----------------------------------------------------------------------|
!!!

!!     UAF = (UARK*(H1+ELRK1)-ALPHA_RK(K)*DTE*RESX/ART)/(H1+ELF1)
!!     VAF = (VARK*(H1+ELRK1)-ALPHA_RK(K)*DTE*RESY/ART)/(H1+ELF1)
!!     UAS = UAF
!!     VAS = VAF
!!   END IF
!!#  endif

   DO I=1,NT
#  if defined (WET_DRY)
     IF(ISWETCE(I)*ISWETC(I) == 1)THEN
#  endif

#      if !defined (NH)
       RESX(I) = ADX2D(I)+ADVUA(I)+DRX2D(I)+PSTX(I)-COR(I)*VA(I)*D1(I)*ART(I)  &
                 -(WUSURF2(I)+WUBOT(I))*ART(I)
       RESY(I) = ADY2D(I)+ADVVA(I)+DRY2D(I)+PSTY(I)+COR(I)*UA(I)*D1(I)*ART(I)  &
                 -(WVSURF2(I)+WVBOT(I))*ART(I)
# if defined (VEGETATION)
       IF(VEGETATION_MODEL.and.VEG_DRAG)THEN
         RESX(I) = RESX(I) + step2d_uveg(i)*ART(I)
         RESY(I) = RESY(I) + step2d_vveg(i)*ART(I)
       END IF
# endif
#      else
       RESX(I) = ADX2D(I)+ADVUA(I)+DRX2D(I)+PSTX(I)-COR(I)*VA(I)*D1(I)*ART(I)  &
                 -(WUSURF2(I)+WUBOT(I))*ART(I)+NHQ2DX(I)
       RESY(I) = ADY2D(I)+ADVVA(I)+DRY2D(I)+PSTY(I)+COR(I)*UA(I)*D1(I)*ART(I)  &
                 -(WVSURF2(I)+WVBOT(I))*ART(I)+NHQ2DY(I)
# if defined (VEGETATION)
       if(VEGETATION_MODEL.and.VEG_DRAG)then
         RESX(I) = RESX(I) + step2d_uveg(i)*ART(I)
         RESY(I) = RESY(I) + step2d_vveg(i)*ART(I)
       endif
# endif
#      endif

#  if defined (SPHERICAL)
       RESX(I) = RESX(I)     &
                 -UA(I)*VA(I)/REARTH*TAN(DEG2RAD*YC(I))*D1(I)*ART(I)
       RESY(I) = RESY(I)     &
                 +UA(I)*UA(I)/REARTH*TAN(DEG2RAD*YC(I))*D1(I)*ART(I)
#  endif


#  if defined (WAVE_CURRENT_INTERACTION) || (NUOPC)
       RESX(I) = RESX(I) + WAVESTRX_2D(I)
       RESY(I) = RESY(I) + WAVESTRY_2D(I)
#  endif  
!
!--UPDATE----------------------------------------------------------------------|
!

       UAF(I)  =  (UARK(I)*(H1(I)+ELRK1(I))-ALPHA_RK(K)*DTE*RESX(I)/ART(I))/(H1(I)+ELF1(I))
       VAF(I)  =  (VARK(I)*(H1(I)+ELRK1(I))-ALPHA_RK(K)*DTE*RESY(I)/ART(I))/(H1(I)+ELF1(I))
#  if defined (WET_DRY)
     ELSE
       UAF(I) = 0.0_SP
       VAF(I) = 0.0_SP
     END IF
#  endif
   END DO

#  if defined (SPHERICAL)
   CALL EXTUV_EDGE_XY(K)
#  endif
   
   VAF(0) = VAFT
   UAF(0) = UAFT

!
!--ADJUST EXTERNAL VELOCITY IN SPONGE REGION-----------------------------------|
!
!old:   UAF = UAF-CC_SPONGE*UAF
!old:   VAF = VAF-CC_SPONGE*VAF
! ---- new: Karsten Lettmann: 2012.06.25 -------
   UAF = UAF/(1.0_SP+CC_SPONGE*UAF**2.0_SP)
   VAF = VAF/(1.0_SP+CC_SPONGE*VAF**2.0_SP)
! ------- end new -------------------------------


!
!--STORE VARIABLES FOR MOMENTUM BALANCE CHECK----------------------------------|
!
#  if defined (BALANCE_2D)
   IF(K == 4) THEN
     TMP=ART*(H1+ELF1)
     ADVUA2 = ADVUA2 + (ADVUA-ADFXA)/TMP/FLOAT(ISPLIT)    !X- HORIZONTAL ADVECTION (m/s^2)
     ADVVA2 = ADVVA2 + (ADVVA-ADFYA)/TMP/FLOAT(ISPLIT)    !X- HORIZONTAL ADVECTION (m/s^2)
     ADFX2  = ADFX2 + ADFXA/TMP/FLOAT(ISPLIT)
     ADFY2  = ADFY2 + ADFYA/TMP/FLOAT(ISPLIT)
     DRX2D2 = DRX2D2 + DRX2D/TMP/FLOAT(ISPLIT)            !X- BAROCLINIC PRESURE GRADIENT FORCE
     DRY2D2 = DRY2D2 + DRY2D/TMP/FLOAT(ISPLIT)            !Y- BAROCLINIC PRESURE GRADIENT FORCE
     CORX2  = CORX2  - COR*VA*D1/(H1+ELF1)/FLOAT(ISPLIT)  !X- CORIOLIS FORCE
     CORY2  = CORY2  + COR*UA*D1/(H1+ELF1)/FLOAT(ISPLIT)  !Y- CORIOLIS FORCE
     PSTX2  = PSTX2  + PSTX/TMP/FLOAT(ISPLIT)             !X- BAROTROPIC PRESURE GRADIENT FORCE
     PSTY2  = PSTY2  + PSTY/TMP/FLOAT(ISPLIT)             !Y- BAROTROPIC PRESURE GRADIENT FORCE
     ADX2D2 = ADX2D2 + ADX2D/TMP/FLOAT(ISPLIT)            !GX (m/s^2)
     ADY2D2 = ADY2D2 + ADY2D/TMP/FLOAT(ISPLIT)            !GY (m/s^2)
     WUSURBF2=WUSURBF2-(WUSURF2+WUBOT)/(H1+ELF1)/FLOAT(ISPLIT) !X-SURFACE & BOTTOM FRICTION
     WVSURBF2=WVSURBF2-(WVSURF2+WVBOT)/(H1+ELF1)/FLOAT(ISPLIT) !Y-SURFACE & BOTTOM FRICTION
     DUDT2  = DUDT2 + (UAF-UARK*(H1+ELRK1)/(H1+ELF1))/DTE/FLOAT(ISPLIT)
     DVDT2  = DVDT2 + (VAF-VARK*(H1+ELRK1)/(H1+ELF1))/DTE/FLOAT(ISPLIT)
     
   END IF     
#  endif

   if(dbg_set(dbg_sbr)) write(ipt,*) "End: extuv_edge.F"

#  endif
   END SUBROUTINE EXTUV_EDGE
!==============================================================================|
