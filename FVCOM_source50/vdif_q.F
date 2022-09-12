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
!   this subroutine is used to calculate the q2 and q2l by solving             !
!   the vertical diffusion equation implicitly.                                !
!==============================================================================|

   SUBROUTINE VDIF_Q              

!------------------------------------------------------------------------------|
   USE ALL_VARS
   USE MOD_UTILS
   USE MOD_WD
   USE MOD_PAR
#  if defined (SEMI_IMPLICIT)
   USE MOD_SEMI_IMPLICIT, ONLY : STAGE, KSTAGE_TE
#  endif

#  if defined (WAVE_CURRENT_INTERACTION)
   USE MOD_WAVE_CURRENT_INTERACTION
#  endif  

#  if defined (THIN_DAM)
   USE MOD_DAM
#  endif

# if defined (VEGETATION)
   USE MOD_VEGETATION
# endif

   IMPLICIT NONE
   INTEGER :: I,J,K,KI
   REAL(SP)  :: CONST1,COEF1,COEF2,COEF3,COEF4,COEF5,LMAX
   REAL(SP), DIMENSION(0:MT,KB) :: GH,SM,SH,PROD,DTEF,KN,BOYGR,A,C,VHP,VH,STF
   REAL(SP), DIMENSION(0:MT) :: L0
   REAL(SP)  :: UTAU2
   REAL(SP)  :: WUSURF_NODE,WVSURF_NODE,WUBOT_NODE,WVBOT_NODE
   REAL(SP)  :: UU_NODE_K,VV_NODE_K,UU_NODE_KM1,VV_NODE_KM1

!--Stability Function Coefficients-------------
   REAL(SP), PARAMETER :: A1    = 0.92_SP
   REAL(SP), PARAMETER :: B1    = 16.6_SP
   REAL(SP), PARAMETER :: A2    = 0.74_SP
   REAL(SP), PARAMETER :: B2    = 10.1_SP
   REAL(SP), PARAMETER :: C1    = 0.08_SP

!--Source/Sink Term Coefficients---------------
   REAL(SP), PARAMETER :: E1    = 1.80_SP
   REAL(SP), PARAMETER :: E2    = 1.33_SP
   REAL(SP), PARAMETER :: SEF   = 1.00_SP

!--Von Karmans Constant------------------------
   REAL(SP), PARAMETER :: KAPPA = 0.40_SP

!--Gravitational Constant----------------------
   REAL(SP), PARAMETER :: GEE   = 9.806_SP

!--Limiting Values of Q2/Q2L-------------------
   REAL(SP), PARAMETER :: SMALL = 1.E-8_SP

!--Coefficient of buoyancy length scale--------
   REAL(SP), PARAMETER :: CB    = 1.E-8_SP

!--Denominator small number threshold----------
   REAL(SP), PARAMETER :: DEN_S = 1.E-10_SP

!--Upper bound of GH for unstable flow---------
   REAL(SP), PARAMETER :: GH_MAX = .0233_SP

!--Lower bound for GH for stable flow----------
   REAL(SP), PARAMETER :: GH_MIN = -.281_SP

#  if defined (WAVE_CURRENT_INTERACTION)
   REAL(SP) :: WIND10,CDD,USTAR,CPP,WAVE_AGE,ALFA_CB
   REAL(SP) :: WUWIND_NODE,WVWIND_NODE 
#  endif  

#  if defined (WAVE_CURRENT_INTERACTION)
!--??????????????????????????????????----------
!   REAL(SP), PARAMETER :: HLC = 2.4_SP 
#  endif  

!----------------------------------------------
   REAL(SP), PARAMETER :: CBCNST = 100.0_SP
   REAL(SP), PARAMETER :: SURFL = 2.E+5_SP
   REAL(SP), PARAMETER :: SHIW = 0.0_SP
   REAL(SP), PARAMETER :: GHC = -6.0_SP

#  if defined (WET_DRY)
   REAL(SP)  :: KQ_TMP,KM_TMP,KH_TMP
#  endif

   IF(dbg_set(dbg_sbr)) write(ipt,*) "Start: vdif_q"
!------------------------------------------------------------------------------|

   ! INITIALIZE MEMORY
   BOYGR = 0.0_SP
   GH    = 0.0_SP
   SM    = 0.0_SP
   SH    = 0.0_SP
   PROD  = 0.0_SP
   DTEF  = 0.0_SP
   KN    = 0.0_SP
   A     = 0.0_SP
   C     = 0.0_SP
   VHP   = 0.0_SP
   VH    = 0.0_SP
   STF   = 0.0_SP
   L0    = 0.0_SP

!------------------------------------------------------------------------------|
!  place a threshold on the turbulence quantities following advection          |
!------------------------------------------------------------------------------|

   DO K = 2, KBM1
     DO I = 1, M
       IF (Q2F(I,K) <= small .OR. Q2LF(I,K) <= small) THEN
!         Q2F(I,K)  = 0.0_SP
!         Q2LF(I,K) = 0.0_SP
         Q2F(I,K)  = small
         Q2LF(I,K) = small
       END IF
     END DO
   END DO

!------------------------------------------------------------------------------!
!  set up coefficients for implicit calculation of vertical diffusion          !
!------------------------------------------------------------------------------!

   DO I=1,M
#  if defined (WET_DRY)
       IF(ISWETN(I) == 1)THEN
#  endif
     DO K=2,KBM1
       A(I,K) = -0.5_SP* DTI * (KQ(I,K+1)+KQ(I,K)+2.0_SP*UMOL)/ &
             (DZZ(I,K-1)*DZ(I,K)*D(I)*D(I))
       C(I,K) = -0.5_SP* DTI*(KQ(I,K-1)+KQ(I,K)+2.0_SP*UMOL)/ &
             (DZZ(I,K-1)*DZ(I,K-1)*D(I)*D(I))
     END DO
#  if defined (WET_DRY)
     END IF
#  endif
   END DO

!------------------------------------------------------------------------------!
!    the following section solves the equation                                 !
!    dti*(kq*q2')' - q2*(2.*dti*dtef+1.) = -q2b                                !
!------------------------------------------------------------------------------!

   CONST1 = 16.6_SP ** .6666667_SP * SEF
   DO I = 1, M
!       VHP(I,1) = SQRT(WUSURF(I)**2+WVSURF(I)**2) * CONST1
!      VH(I,1)  = 0.0_SP

     WUSURF_NODE = 0.0_SP
     WVSURF_NODE = 0.0_SP
     DO J=1,NTVE(I)
       WUSURF_NODE = WUSURF_NODE + WUSURF(NBVE(I,J))
       WVSURF_NODE = WVSURF_NODE + WVSURF(NBVE(I,J))
     END DO
     WUSURF_NODE = WUSURF_NODE/FLOAT(NTVE(I))
     WVSURF_NODE = WVSURF_NODE/FLOAT(NTVE(I))

!#    if defined (WAVE_CURRENT_INTERACTION) && !defined (WAVE_OFFLINE)
#    if defined (WAVE_CURRENT_INTERACTION)
!!$     UTAU2 = SQRT(WUSURF_NODE**2+WVSURF_NODE**2)
!!$     VHP(I,1) = (15.8_SP*CBCNST)**0.6666667_SP*UTAU2
!!$     VH(I,1)  = 0.0_SP
!!$!     L0(I)    = SURFL*UTAU2/GEE
!!$     L0(I)    = SURFL*UTAU2/GRAV_N(I)
     UTAU2 = SQRT(WUSURF_NODE**2+WVSURF_NODE**2)

     WUWIND_NODE = 0.0_SP
     WVWIND_NODE = 0.0_SP
     DO J=1,NTVE(I)
#    if defined (WAVE_OFFLINE)
       WUWIND_NODE = WUWIND_NODE + UUWIND(NBVE(I,J))
       WVWIND_NODE = WVWIND_NODE + VVWIND(NBVE(I,J))
#    else
       WUWIND_NODE = WUWIND_NODE + UWWIND(NBVE(I,J))
       WVWIND_NODE = WVWIND_NODE + VWWIND(NBVE(I,J))
#    endif
     END DO
     WUWIND_NODE = WUWIND_NODE/FLOAT(NTVE(I))
     WVWIND_NODE = WVWIND_NODE/FLOAT(NTVE(I))

     WIND10=WUWIND_NODE*WUWIND_NODE + WVWIND_NODE*WVWIND_NODE
     WIND10=SQRT(WIND10)
     IF(WIND10 < 7.5)THEN
       CDD = 1.2875E-3_SP
     ELSE
       CDD = (0.8+0.065*WIND10)*1.0E-3_SP
     END IF
     USTAR=CDD*WIND10*WIND10
     USTAR=SQRT(USTAR)
     
     CPP = WLEN(I)/(TPEAK(I) + small)
     
     WAVE_AGE = CPP/(USTAR + small)
     ALFA_CB = 0.04_SP*WAVE_AGE
     ALFA_CB = ALFA_CB**4.0_SP
     ALFA_CB = EXP(-ALFA_CB)
     ALFA_CB = 15.0_SP*WAVE_AGE*ALFA_CB
     
     VHP(I,1) = ALFA_CB*UTAU2**1.5_SP    !1.5=3.0/2.0
     VH(I,1) = 0.0_SP
     L0(I) = 0.85*HSC1(I)
     
#    else
     UTAU2 = SQRT(WUSURF_NODE**2+WVSURF_NODE**2)
     VHP(I,1) = (15.8_SP*CBCNST)**0.6666667_SP*UTAU2
     VH(I,1)  = 0.0_SP
!     L0(I)    = SURFL*UTAU2/GEE
     L0(I)    = SURFL*UTAU2/GRAV_N(I)
!     HLC   = 2.4_SP          
!      L0(I) = HLC*HSC1(I)            
#    endif


     WUBOT_NODE = 0.0_SP
     WVBOT_NODE = 0.0_SP
     DO J=1,NTVE(I)
       WUBOT_NODE = WUBOT_NODE + WUBOT(NBVE(I,J))
       WVBOT_NODE = WVBOT_NODE + WVBOT(NBVE(I,J))
     END DO
     WUBOT_NODE = WUBOT_NODE/FLOAT(NTVE(I))
     WVBOT_NODE = WVBOT_NODE/FLOAT(NTVE(I))

     Q2F(I,KB) = SQRT(WUBOT_NODE**2+WVBOT_NODE**2) * CONST1
   END DO

   Q2  = ABS(Q2+small)
   Q2L = ABS(Q2L)
 
!------------------------------------------------------------------------------!
!  calculate boygr = -Brunt Vaisala Frequency^2                                !
!  calculate internal wave shear energy contribution: prod = 200*(b.v.freq)**3 !
!------------------------------------------------------------------------------!

   DO K=2,KBM1
     DO I=1,M
#  if defined (WET_DRY)
      IF(ISWETN(I) == 1)THEN
#  endif
!       BOYGR(I,K) = GEE * (RHO1(I,K-1)-RHO1(I,K))/(DZZ(I,K-1)*D(I))
       BOYGR(I,K) = GRAV_N(I) * (RHO1(I,K-1)-RHO1(I,K))/(DZZ(I,K-1)*D(I))
       PROD(I,K)  = KM(I,K) * 0.0_SP * (.5_SP*(-BOYGR(I,K)+ABS(BOYGR(I,K)))) ** 1.5_SP
#  if defined (WET_DRY)
      END IF
#  endif
     END DO
   END DO

!------------------------------------------------------------------------------!
!  calculate shear production source term = prod                               !
!  calculate buoyancy production source term = kh*boygr                        !
!------------------------------------------------------------------------------!
   DO  K = 2, KBM1
     DO  I = 1, M
#  if !defined (WET_DRY)
       IF (D(I) > 0.0_SP) THEN
#  else
       IF(ISWETN(I) == 1)THEN
#  endif

         UU_NODE_K   = 0.0_SP
         VV_NODE_K   = 0.0_SP
         UU_NODE_KM1 = 0.0_SP
         VV_NODE_KM1 = 0.0_SP
         DO J=1,NTVE(I)
           UU_NODE_K   = UU_NODE_K + U(NBVE(I,J),K)
           VV_NODE_K   = VV_NODE_K + V(NBVE(I,J),K)
           UU_NODE_KM1 = UU_NODE_KM1 + U(NBVE(I,J),K-1)
           VV_NODE_KM1 = VV_NODE_KM1 + V(NBVE(I,J),K-1)
         END DO
         UU_NODE_K   = UU_NODE_K/FLOAT(NTVE(I))
         VV_NODE_K   = VV_NODE_K/FLOAT(NTVE(I))
         UU_NODE_KM1 = UU_NODE_KM1/FLOAT(NTVE(I))
         VV_NODE_KM1 = VV_NODE_KM1/FLOAT(NTVE(I))

         PROD(I,K) = PROD(I,K) + KM(I,K) * SEF * &
                     ((UU_NODE_K-UU_NODE_KM1)**2+(VV_NODE_K-VV_NODE_KM1)**2)/ &
                     (DZZ(I,K-1)*D(I))**2
         PROD(I,K) = PROD(I,K) + KH(I,K) * BOYGR(I,K)
       END IF
     END DO
   END DO
!------------------------------------------------------------------------------!
!  solve for turbulent length scale l = q2l/q2                                 !
!------------------------------------------------------------------------------!

!   L = Q2L/Q2
    DO K=2,KBM1
      DO I=1,M
        L(I,K) = Q2L(I,K)/Q2(I,K)
      END DO
    END DO

!------------------------------------------------------------------------------!
!  in stably stratified regions, length scale is limited by buoyancy length    !
!  scale, l_b = c_b*(k^.5/N).  length scale limitation also puts an effective  !
!  lower bound of -.281 on gh in stably stratified regions                     !
!  length scale near surface is modified with wind mixed length scale          !
!------------------------------------------------------------------------------!

   DO K=2,KBM1
     DO I=1,M
       IF(Z(I,K) > -0.5_SP)L(I,K) = MAX(L(I,K),KAPPA*L0(I)) 
       IF(BOYGR(I,K) < 0.0_SP)THEN
         LMAX      = SQRT(ABS(GH_MIN) * Q2(I,K) / (-BOYGR(I,K) + SMALL) )
         L(I,K)    = MIN(L(I,K),LMAX)
         Q2L(I,K) = Q2(I,K)*L(I,K)
       END IF
     END DO
   END DO

!------------------------------------------------------------------------------!
!  solve for gh = (-l^2*N^2/q^2)                                               !
!------------------------------------------------------------------------------!

!!   GH = (L**2/Q2)*BOYGR
!   GH(1:MT,1:KBM1) = (L(1:MT,1:KBM1)**2/Q2(1:MT,1:KBM1))*BOYGR(1:MT,1:KBM1)
   GH(1:M,1:KBM1) = (L(1:M,1:KBM1)**2/Q2(1:M,1:KBM1))*BOYGR(1:M,1:KBM1)

!------------------------------------------------------------------------------!
!  limit maximum of GH in unstable regions to reasonable physical value        !
!  limiting GH also constrains stability functions SH/SM to finite values      !
!------------------------------------------------------------------------------!

   GH = MIN(GH,GH_MAX)

   L(:,1)   = KAPPA*L0(:)         !0.0_SP
   L(:,KB)  = 0.0_SP
   GH(:,1)  = 0.0_SP
   GH(:,KB) = 0.0_SP


!------------------------------------------------------------------------------!
!  calculate eddy kinetic energy dissipation rate: dtef                        !
!------------------------------------------------------------------------------!

   STF = 1.0_SP

   IF(.not.SURFACE_WAVE_MIXING)THEN
   DO K = 1,KB
     DO I = 1,M
       IF(GH(I,K) < 0.0_SP) STF(I,K)=1.0-0.9*(GH(I,K)/GHC)**1.5
       IF(GH(I,K) < GHC) STF(I,K) = 0.1
     END DO
   END DO
   END IF

!   DTEF = Q2*SQRT(Q2)/(B1*Q2L + small)*STF
   DTEF = SQRT(Q2)/(B1*L + small)*STF

   DO K = 2, KBM1
     DO I = 1, M
       VHP(I,K) = 1. / (A(I,K)+C(I,K)*(1.-VH(I,K-1))-(2.*DTI*DTEF(I,K)+1.))
       VH(I,K) = A(I,K) * VHP(I,K)
       VHP(I,K) = (-2.*DTI*PROD(I,K)+C(I,K)*VHP(I,K-1)-Q2F(I,K))*VHP(I,K)
     END DO
   END DO

   DO K=1,KBM1
     KI = KB-K
     DO I = 1, M
#  if defined (WET_DRY)
     IF(ISWETN(I) == 1)THEN
#  endif
       Q2F(I,KI) = VH(I,KI) * Q2F(I,KI+1) + VHP(I,KI)
#  if defined (WET_DRY)
     END IF
#  endif
     END DO
   END DO

!------------------------------------------------------------------------------!
!      the following section solves the equation                               !
!      dti*(kq*q2l')' - q2l*(dti*dtef+1.) = -q2lb                              !
!------------------------------------------------------------------------------!

   VH(:,1)  = 0.0_SP
   VHP(:,1) = 0.0_SP

   DO  K = 2, KBM1
     DO  I = 1, M
#  if !defined (WET_DRY)
       IF (D(I) > 0.0_SP) THEN
#  else
       IF(ISWETN(I) == 1)THEN
#  endif
         DTEF(I,K) = DTEF(I,K) * (1.+E2*((1./ABS(Z(I,K)-Z(I,1))+1./ &
                     ABS(Z(I,K)-Z(I,KB)))*L(I,K)/(D(I)*KAPPA))**2)
         VHP(I,K)  = 1. / (A(I,K)+C(I,K)*(1.-VH(I,K-1))- &
                     (DTI*DTEF(I,K)+1.))
         VH(I,K)   = A(I,K) * VHP(I,K)
         VHP(I,K)  = (DTI*(-PROD(I,K)*L(I,K)*E1)+C(I,K)* &
                     VHP(I,K-1)-Q2LF(I,K)) * VHP(I,K)
       END IF
     END DO
   END DO


   DO K=1,KBM1
     KI = KB - K
     DO I = 1, M
#  if defined (WET_DRY)
     IF(ISWETN(I) == 1)THEN
#  endif
       Q2LF(I,KI) = VH(I,KI) * Q2LF(I,KI+1) + VHP(I,KI)
#  if defined (WET_DRY)
     END IF
#  endif
     END DO
   END DO

!------------------------------------------------------------------------------|
!  place a threshold on the turbulence quantities following vertical diffusion |
!------------------------------------------------------------------------------|

#  if defined (SEMI_IMPLICIT)
   IF(STAGE==KSTAGE_TE) THEN
#  endif

   DO K = 2, KBM1
     DO I = 1, M
       IF (Q2F(I,K) <= small .OR. Q2LF(I,K) <= small) THEN
         Q2F(I,K)  = small
         Q2LF(I,K) = small
       END IF
     END DO
   END DO


!===========THRESHOLD ON LENGTH SCALE ALT STYLE (Active in Barotropic Case)====
!   DO K=2,KBM1
!     DO I=1,N
!       LMAX      = SQRT(ABS(GH_MIN) * Q2F(I,K) / MAX(0.0_SP,-BOYGR(I,K)+SMALL))
!       L(I,K)    = MIN(L(I,K),LMAX)
!       Q2LF(I,K) = Q2F(I,K)*L(I,K)
!     END DO
!   END DO
!------------------------------------------------------------------------------!

!   L(:,1)   = 0.0_SP
!   L(:,KB)  = 0.0_SP
!   GH(:,1)  = 0.0_SP
!   GH(:,KB) = 0.0_SP


!------------------------------------------------------------------------------!
!  calculate stability functions sh + sm using Galperin 1988 formulation       !
!------------------------------------------------------------------------------!
   COEF4 = 18.0_SP * A1 * A1 + 9.0_SP * A1 * A2
   COEF5 = 9.0_SP * A1 * A2

   DO K = 1,KB
     DO I = 1,M
       COEF1 = A2 * (1.0_SP-6.0_SP*A1/B1*STF(i,k))
       COEF2 = 3.0_SP * A2 * B2/STF(i,k) + 18.0_SP * A1 * A2
       COEF3 = A1 * (1.0_SP-3.0_SP*C1-6.0_SP*A1/B1*STF(i,k))

       SH(i,k) = COEF1/(1.0_SP-COEF2*GH(i,k))
       SM(i,k) = (COEF3+SH(i,k)*COEF4*GH(i,k))/(1.-COEF5*GH(i,k))
     END DO
   END DO

!------------------------------------------------------------------------------!
!  calculate turbulent eddy viscosities/diffusivities kq,km,kh                 !
!------------------------------------------------------------------------------!

   KN = L*SQRT(Q2)
   KQ = 0.5_SP*(KN*KAPPA*SM+KQ)
   KM = 0.5_SP*(KN*SM+KM)
   KH = 0.5_SP*(KN*SH*VPRNU+KH)

#  if defined (WET_DRY)
    DO I=1,M
     IF(H(I) <= STATIC_SSH_ADJ)THEN
      KQ_TMP = 0.0_SP
      KM_TMP = 0.0_SP
      KH_TMP = 0.0_SP
      DO K=2,KBM1
        KQ_TMP = KQ_TMP + KQ(I,K)
        KM_TMP = KM_TMP + KM(I,K)
        KH_TMP = KH_TMP + KH(I,K)
      END DO
        KQ_TMP = KQ_TMP/(KBM1-1)
        KM_TMP = KM_TMP/(KBM1-1)
        KH_TMP = KH_TMP/(KBM1-1)
      DO K=2,KBM1
        KQ(I,K) = KQ_TMP
        KM(I,K) = KM_TMP
        KH(I,K) = KH_TMP
      END DO
     END IF
    END DO
#  endif

# if defined (VEGETATION)
    !
    !-----------------------------------------------------------------------
    !  Add the effect of vegetation on horizontal viscosity.
    !-----------------------------------------------------------------------
    !
    if(VEGETATION_MODEL.and.VEG_HMIXING)then
        KM = KM  + VEG % visc3d_r_veg
    end if
#  endif

#  if defined (THIN_DAM)
   DO K=1,KBM1
     DO I=1,NODE_DAM1_N
       IF(K<=KDAM(I_NODE_DAM1_N(I,1)).AND.K<=KDAM(I_NODE_DAM1_N(I,2)))THEN
         KH(I_NODE_DAM1_N(I,1),K) = SUM(KH(I_NODE_DAM1_N(I,1:2),K))/2.0
         KH(I_NODE_DAM1_N(I,2),K) = KH(I_NODE_DAM1_N(I,1),K)
       END IF
     END DO
     DO I=1,NODE_DAM2_N
       IF(K<=KDAM(I_NODE_DAM2_N(I,1)).AND.K<=KDAM(I_NODE_DAM2_N(I&
            &,2)).AND.K<=KDAM(I_NODE_DAM2_N(I,3)))THEN
         KH(I_NODE_DAM2_N(I,1),K) = SUM(KH(I_NODE_DAM2_N(I,1:3),K))/3.0
         KH(I_NODE_DAM2_N(I,2),K) = KH(I_NODE_DAM2_N(I,1),K)
         KH(I_NODE_DAM2_N(I,3),K) = KH(I_NODE_DAM2_N(I,1),K)
       END IF
     END DO
     DO I=1,NODE_DAM3_N
       IF(K<=KDAM(I_NODE_DAM3_N(I,1)).AND.K<=KDAM(I_NODE_DAM3_N(I&
      &,2)).AND.K<=KDAM(I_NODE_DAM3_N(I,3)).AND.K<=KDAM(I_NODE_DAM3_N(I,4)))THEN
         KH(I_NODE_DAM3_N(I,1),K) = SUM(KH(I_NODE_DAM3_N(I,1:4),K))/4.0
         KH(I_NODE_DAM3_N(I,2),K) = KH(I_NODE_DAM3_N(I,1),K)
         KH(I_NODE_DAM3_N(I,3),K) = KH(I_NODE_DAM3_N(I,1),K)
         KH(I_NODE_DAM3_N(I,4),K) = KH(I_NODE_DAM3_N(I,1),K)
       END IF
     END DO
   END DO
#  endif


#  if defined (SEMI_IMPLICIT)
   ENDIF
#  endif

   IF(dbg_set(dbg_sbr)) write(ipt,*) "End: vdif_q"

   END SUBROUTINE VDIF_Q
!==============================================================================|
