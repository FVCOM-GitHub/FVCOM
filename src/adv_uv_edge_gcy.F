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

!==============================================================================!

   SUBROUTINE ADV_UV_EDGE_GCY

!==============================================================================!
! this subroutine calculate advective, coriolis, pressure gradient, etc in     !
! x and y momentum equations except vertical diffusion terms for internal mode ! 
!==============================================================================!

   USE ALL_VARS
   USE MOD_UTILS
   USE MOD_SPHERICAL
   USE MOD_NORTHPOLE
   USE BCS
   USE MOD_WD

#  if defined (ONE_D_MODEL)
   USE MOD_ONEDTIDE
#  endif   

#  if defined (SEMI_IMPLICIT)
   USE MOD_SEMI_IMPLICIT
#  endif

#  if defined (NH)
   USE NON_HYDRO, ONLY: NHQDRX, NHQDRY
#  endif

#  if defined (WAVE_CURRENT_INTERACTION)
   USE MOD_WAVE_CURRENT_INTERACTION
#  endif

   IMPLICIT NONE
   REAL(SP) :: XFLUX(0:NT,KB),YFLUX(0:NT,KB)
   REAL(SP) :: PSTX_TM(0:NT,KB),PSTY_TM(0:NT,KB)
   REAL(SP) :: COFA1,COFA2,COFA3,COFA4,COFA5,COFA6,COFA7,COFA8
   REAL(SP) :: XADV,YADV,TXXIJ,TYYIJ,TXYIJ,UN
   REAL(SP) :: VISCOF,VISCOF1,VISCOF2,TEMP,TPA,TPB
   REAL(SP) :: XIJA,YIJA,XIJB,YIJB,UIJ,VIJ
   REAL(SP) :: DIJ,ELIJ,TMPA,TMPB,TMP,XFLUXV,YFLUXV
   REAL(SP) :: FACT,FM1,EXFLUX,ISWETTMP
   INTEGER  :: I,IA,IB,J1,J2,K1,K2,K3,K4,K5,K6,K,II,J,I1,I2

# if defined (SPHERICAL)
   REAL(DP) :: XTMP,XTMP1
   REAL(SP) :: U_TMP,V_TMP,UF_TMP,VF_TMP
# endif

#  if defined (LIMITED_NO)
   REAL(SP) :: UIJ1,VIJ1,UIJ2,VIJ2,FXX,FYY
   REAL(SP) :: UK1(KB),UK2(KB),UK3(KB),UK4(KB),UK5(KB),UK6(KB), &
               VK1(KB),VK2(KB),VK3(KB),VK4(KB),VK5(KB),VK6(KB)
#  else
   REAL(SP),ALLOCATABLE,DIMENSION(:) :: UIJ1,VIJ1,UIJ2,VIJ2,FXX,FYY
   REAL(SP),ALLOCATABLE,DIMENSION(:) :: UALFA,VALFA
   REAL(SP) :: UALFA_TMP,VALFA_TMP
   INTEGER :: ERROR
   REAL(SP) :: EPS

   REAL(SP) :: UK1,UK2,UK3,UK4,UK5,UK6,VK1,VK2,VK3,VK4,VK5,VK6
#  endif

#  if defined (NH) || defined (LIMITER_VER_ADV)
   REAL(SP)  CONV_U(1:KB), CONV_V(1:KB), DISS_U(1:KB), DISS_V(1:KB)
   REAL(SP)  SL_H(0:KB), U_TEMP(0:KB), V_TEMP(0:KB)
   REAL(SP)  SL_U, SL_F
#  endif

   IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "Start: adv_uv_edge_gcy"

!------------------------------------------------------------------------------!
   SELECT CASE(HORIZONTAL_MIXING_TYPE)
   CASE ('closure')
      FACT = 1.0_SP
      FM1  = 0.0_SP
   CASE('constant')
      FACT = 0.0_SP
      FM1  = 1.0_SP
   CASE DEFAULT
      CALL FATAL_ERROR("UNKNOW HORIZONTAL MIXING TYPE:",&
           & TRIM(HORIZONTAL_MIXING_TYPE) )
   END SELECT

!
!-----Initialize Flux Variables------------------------------------------------!
!
   VISCOFM = 0.0_SP
   
   XFLUX  = 0.0_SP
   YFLUX  = 0.0_SP
   PSTX_TM = 0.0_SP
   PSTY_TM = 0.0_SP
# if defined (SEMI_IMPLICIT)
   XFLUX3= 0.0_SP
   YFLUX3= 0.0_SP
# endif

!
!-----Loop Over Edges and Accumulate Flux--------------------------------------!
!
#  if !defined (LIMITED_NO)
   ALLOCATE(UIJ1(NE),VIJ1(NE),UIJ2(NE),VIJ2(NE),STAT=ERROR)
   IF(ERROR /= 0) &
   & CALL FATAL_ERROR("The arrays UIJ1,VIJ1,UIJ2 and VIJ2 can not be allocated.")
   
   ALLOCATE(UALFA(0:NT),VALFA(0:NT),STAT=ERROR)
   IF(ERROR /= 0) &
   & CALL FATAL_ERROR("The arrays UALFA,VALFA can not be allocated.")
   
   ALLOCATE(FXX(NE),FYY(NE),STAT=ERROR)
   IF(ERROR /= 0) &
   & CALL FATAL_ERROR("The arrays FXX,FYY can not be allocated.")

   
   DO K=1, KBM1
     UIJ1=0.0_SP;VIJ1=0.0_SP;UIJ2=0.0_SP;VIJ2=0.0_SP
     UALFA=1.0_SP;VALFA=1.0_SP
     FXX=0.0_SP;FYY=0.0_SP

     DO I=1,NE
       IA=IEC(I,1)
       IB=IEC(I,2)
       J1=IENODE(I,1)
       J2=IENODE(I,2)

       K1=NBE(IA,1)
       K2=NBE(IA,2)
       K3=NBE(IA,3)
       K4=NBE(IB,1)
       K5=NBE(IB,2)
       K6=NBE(IB,3)
#      if defined (SPHERICAL)
       XIJA=DLTXNE(I,1)
       YIJA=DLTYNE(I,1)
       XIJB=DLTXNE(I,2)
       YIJB=DLTYNE(I,2)
#      else
       XIJA=XIJC(I)-XC(IA)
       YIJA=YIJC(I)-YC(IA)
       XIJB=XIJC(I)-XC(IB)
       YIJB=YIJC(I)-YC(IB)
#      endif

       UK1 = U(K1,K)
       UK2 = U(K2,K)
       UK3 = U(K3,K)
       UK4 = U(K4,K)
       UK5 = U(K5,K)
       UK6 = U(K6,K)
       VK1 = V(K1,K)
       VK2 = V(K2,K)
       VK3 = V(K3,K)
       VK4 = V(K4,K)
       VK5 = V(K5,K)
       VK6 = V(K6,K)
       IF(K1 == 0) CALL GHOSTUV3(IA,1,UK1,VK1,K)
       IF(K2 == 0) CALL GHOSTUV3(IA,2,UK2,VK2,K)
       IF(K3 == 0) CALL GHOSTUV3(IA,3,UK3,VK3,K)
       IF(K4 == 0) CALL GHOSTUV3(IB,1,UK4,VK4,K)
       IF(K5 == 0) CALL GHOSTUV3(IB,2,UK5,VK5,K)
       IF(K6 == 0) CALL GHOSTUV3(IB,3,UK6,VK6,K)

       DIJ=0.5_SP*(DT(J1)*DZ(J1,K)+DT(J2)*DZ(J2,K))
#      if defined (WET_DRY)
#      if !defined (SEMI_IMPLICIT)
       IF(ISWETCT(IA)*ISWETC(IA) == 1 .OR. ISWETCT(IB)*ISWETC(IB) == 1)THEN
#      else
       IF(ISWETCT(IA) == 1 .OR. ISWETCT(IB) == 1)THEN
#      endif
#      endif

       COFA1=A1U(IA,1)*U(IA,K)+A1U(IA,2)*UK1+A1U(IA,3)*UK2+A1U(IA,4)*UK3
       COFA2=A2U(IA,1)*U(IA,K)+A2U(IA,2)*UK1+A2U(IA,3)*UK2+A2U(IA,4)*UK3
       COFA5=A1U(IA,1)*V(IA,K)+A1U(IA,2)*VK1+A1U(IA,3)*VK2+A1U(IA,4)*VK3
       COFA6=A2U(IA,1)*V(IA,K)+A2U(IA,2)*VK1+A2U(IA,3)*VK2+A2U(IA,4)*VK3
!       COFA1=A1U(IA,1)*U(IA,K)+A1U(IA,2)*U(K1,K)+A1U(IA,3)*U(K2,K)+A1U(IA,4)*U(K3,K)
!       COFA2=A2U(IA,1)*U(IA,K)+A2U(IA,2)*U(K1,K)+A2U(IA,3)*U(K2,K)+A2U(IA,4)*U(K3,K)
!       COFA5=A1U(IA,1)*V(IA,K)+A1U(IA,2)*V(K1,K)+A1U(IA,3)*V(K2,K)+A1U(IA,4)*V(K3,K)
!       COFA6=A2U(IA,1)*V(IA,K)+A2U(IA,2)*V(K1,K)+A2U(IA,3)*V(K2,K)+A2U(IA,4)*V(K3,K)

       UIJ1(I)=COFA1*XIJA+COFA2*YIJA
       VIJ1(I)=COFA5*XIJA+COFA6*YIJA
       UALFA_TMP=ABS(U(IA,K)-U(IB,K))/ABS(UIJ1(I)+EPSILON(EPS))
       VALFA_TMP=ABS(V(IA,K)-V(IB,K))/ABS(VIJ1(I)+EPSILON(EPS))
       IF(UALFA_TMP > 1)UALFA_TMP = 1.0_SP
       IF(VALFA_TMP > 1)VALFA_TMP = 1.0_SP
       UALFA(IA)=MIN(UALFA(IA),UALFA_TMP)
       VALFA(IA)=MIN(VALFA(IA),VALFA_TMP)

       COFA3=A1U(IB,1)*U(IB,K)+A1U(IB,2)*UK4+A1U(IB,3)*UK5+A1U(IB,4)*UK6
       COFA4=A2U(IB,1)*U(IB,K)+A2U(IB,2)*UK4+A2U(IB,3)*UK5+A2U(IB,4)*UK6
       COFA7=A1U(IB,1)*V(IB,K)+A1U(IB,2)*VK4+A1U(IB,3)*VK5+A1U(IB,4)*VK6
       COFA8=A2U(IB,1)*V(IB,K)+A2U(IB,2)*VK4+A2U(IB,3)*VK5+A2U(IB,4)*VK6
!       COFA3=A1U(IB,1)*U(IB,K)+A1U(IB,2)*U(K4,K)+A1U(IB,3)*U(K5,K)+A1U(IB,4)*U(K6,K)
!       COFA4=A2U(IB,1)*U(IB,K)+A2U(IB,2)*U(K4,K)+A2U(IB,3)*U(K5,K)+A2U(IB,4)*U(K6,K)
!       COFA7=A1U(IB,1)*V(IB,K)+A1U(IB,2)*V(K4,K)+A1U(IB,3)*V(K5,K)+A1U(IB,4)*V(K6,K)
!       COFA8=A2U(IB,1)*V(IB,K)+A2U(IB,2)*V(K4,K)+A2U(IB,3)*V(K5,K)+A2U(IB,4)*V(K6,K)

       UIJ2(I)=COFA3*XIJB+COFA4*YIJB
       VIJ2(I)=COFA7*XIJB+COFA8*YIJB
       UALFA_TMP=ABS(U(IA,K)-U(IB,K))/ABS(UIJ2(I)+EPSILON(EPS))
       VALFA_TMP=ABS(V(IA,K)-V(IB,K))/ABS(VIJ2(I)+EPSILON(EPS))
       IF(UALFA_TMP > 1)UALFA_TMP = 1.0_SP
       IF(VALFA_TMP > 1)VALFA_TMP = 1.0_SP
       UALFA(IB)=MIN(UALFA(IB),UALFA_TMP)
       VALFA(IB)=MIN(VALFA(IB),VALFA_TMP)

!
!-------ADD THE VISCOUS TERM & ADVECTION TERM---------------------------------!
!

       VISCOF1=ART(IA)*SQRT(COFA1**2+COFA6**2+0.5_SP*(COFA2+COFA5)**2)
       VISCOF2=ART(IB)*SQRT(COFA3**2+COFA8**2+0.5_SP*(COFA4+COFA7)**2)

!       VISCOF=FACT*0.5_SP*HORCON*(VISCOF1+VISCOF2)/HPRNU + FM1*HORCON
!       VISCOF=FACT*0.5_SP*HORCON*(VISCOF1+VISCOF2)/HPRNU + FM1*HORCON/HPRNU
       ! David moved HPRNU and added HVC
       VISCOF=(FACT*0.5_SP*(VISCOF1*CC_HVC(IA)+VISCOF2*CC_HVC(IB)) + FM1*0.5_SP*(CC_HVC(IA)+CC_HVC(IB)))/HPRNU
       VISCOFM(IA,K) = VISCOFM(IA,K) + VISCOF
       VISCOFM(IB,K) = VISCOFM(IB,K) + VISCOF

       TXXIJ=(COFA1+COFA3)*VISCOF
       TYYIJ=(COFA6+COFA8)*VISCOF
       TXYIJ=0.5_SP*(COFA2+COFA4+COFA5+COFA7)*VISCOF
       FXX(I)=DIJ*(TXXIJ*DLTYC(I)-TXYIJ*DLTXC(I))
       FYY(I)=DIJ*(TXYIJ*DLTYC(I)-TYYIJ*DLTXC(I))

# if defined (PLBC)
       TXYIJ=0.5_SP*(COFA5+COFA7)*VISCOF
       FXX(I)=DIJ*(TXXIJ*DLTYC(I)-0.0_SP*DLTXC(I))
       FYY(I)=DIJ*(TXYIJ*DLTYC(I)-0.0_SP*DLTXC(I))
# endif

#      if defined (WET_DRY)
       END IF
#      endif
     END DO

     DO I=1,NE
       IA=IEC(I,1)
       IB=IEC(I,2)
       J1=IENODE(I,1)
       J2=IENODE(I,2)

# if !defined (SEMI_IMPLICIT)

       ELIJ=0.5_SP*(EGF(J1)+EGF(J2))

#      if defined (AIR_PRESSURE)
       ELIJ=ELIJ-0.5_SP*(EGF_AIR(J1)+EGF_AIR(J2))
#      endif

#      if defined (EQUI_TIDE)
       ELIJ=ELIJ-0.5_SP*(EGF_EQI(J1)+EGF_EQI(J2))
#      endif
#      if defined (ATMO_TIDE)
       ELIJ=ELIJ-0.5_SP*(EGF_ATMO(J1)+EGF_ATMO(J2))
#      endif       

# else

       IF(STAGE<KSTAGE_UV) THEN
         ELIJ=0.5_SP*(ET(J1)+ET(J2))
       ELSE
         ELIJ=(1.0_SP-IFCETA)*0.5_SP*(ET(J1)+ET(J2))
       ENDIF
#      if defined (AIR_PRESSURE)
       ELIJ=ELIJ-( (1.0_SP-IFCETA)*0.5_SP*(EL_AIR(J1)+EL_AIR(J2))+IFCETA*0.5_SP*(ELF_AIR(J1)+ELF_AIR(J2)) ) 
#      endif
#      if defined (EQUI_TIDE)
       ELIJ=ELIJ-( (1.0_SP-IFCETA)*0.5_SP*(EL_EQI(J1)+EL_EQI(J2))+IFCETA*0.5_SP*(ELF_EQI(J1)+ELF_EQI(J2)) )
#      endif
#      if defined (ATMO_TIDE)
       ELIJ=ELIJ-( (1.0_SP-IFCETA)*0.5_SP*(EL_ATMO(J1)+EL_ATMO(J2))+IFCETA*0.5_SP*(ELF_ATMO(J1)+ELF_ATMO(J2)) )
#      endif

# endif

       DIJ=0.5_SP*(DT(J1)*DZ(J1,K)+DT(J2)*DZ(J2,K))
#      if defined (WET_DRY)
#      if !defined (SEMI_IMPLICIT)
       IF(ISWETCT(IA)*ISWETC(IA) == 1 .OR. ISWETCT(IB)*ISWETC(IB) == 1)THEN
#      else
       IF(ISWETCT(IA) == 1 .OR. ISWETCT(IB) == 1)THEN
#      endif
#      endif
       UIJ1(I)=U(IA,K)+UALFA(IA)*UIJ1(I)
       VIJ1(I)=V(IA,K)+VALFA(IA)*VIJ1(I)
       UIJ2(I)=U(IB,K)+UALFA(IB)*UIJ2(I)
       VIJ2(I)=V(IB,K)+VALFA(IB)*VIJ2(I)

#      if defined (LIMITED_1)
       IF(UIJ1(I) > MAX(U(IA,K),U(IB,K)) .OR. UIJ1(I) < MIN(U(IA,K),U(IB,K)) .OR. &
          UIJ2(I) > MAX(U(IA,K),U(IB,K)) .OR. UIJ2(I) < MIN(U(IA,K),U(IB,K)))THEN
         UIJ1(I)=U(IA,K)
         UIJ2(I)=U(IB,K)
       END IF

       IF(VIJ1(I) > MAX(V(IA,K),V(IB,K)) .OR. VIJ1(I) < MIN(V(IA,K),V(IB,K)) .OR. &
          VIJ2(I) > MAX(V(IA,K),V(IB,K)) .OR. VIJ2(I) < MIN(V(IA,K),V(IB,K)))THEN
         VIJ1(I)=V(IA,K)
         VIJ2(I)=V(IB,K)
       END IF
#      endif

!      NORMAL VELOCITY
       UIJ=0.5_SP*(UIJ1(I)+UIJ2(I))
       VIJ=0.5_SP*(VIJ1(I)+VIJ2(I))
       EXFLUX = DIJ*(-UIJ*DLTYC(I) + VIJ*DLTXC(I))

#     if defined (PLBC)
       EXFLUX = DIJ*(-UIJ*DLTYC(I) + 0.0_SP*DLTXC(I))
#     endif

       XADV=EXFLUX*((1.0_SP-SIGN(1.0_SP,EXFLUX))*UIJ2(I)+(1.0_SP+SIGN(1.0_SP,EXFLUX))*UIJ1(I))*0.5_SP
       YADV=EXFLUX*((1.0_SP-SIGN(1.0_SP,EXFLUX))*VIJ2(I)+(1.0_SP+SIGN(1.0_SP,EXFLUX))*VIJ1(I))*0.5_SP

       !!CALCULATE BOUNDARY FLUX AUGMENTERS
       TPA = FLOAT(1-ISBC(I))*EPOR(IA)
       TPB = FLOAT(1-ISBC(I))*EPOR(IB)

       !!ACCUMULATE ADVECTIVE + DIFFUSIVE + BAROTROPIC PRESSURE GRADIENT TERMS
       XFLUX(IA,K)=XFLUX(IA,K)+XADV*TPA+FXX(I)*TPA
       YFLUX(IA,K)=YFLUX(IA,K)+YADV*TPA+FYY(I)*TPA
       XFLUX(IB,K)=XFLUX(IB,K)-XADV*TPB-FXX(I)*TPB
       YFLUX(IB,K)=YFLUX(IB,K)-YADV*TPB-FYY(I)*TPB

#  if defined (WET_DRY)
    END IF
#  endif
! for spherical coordinator and domain across 360^o latitude
#    if defined (SPHERICAL)
        XTMP  = VX(J2)*TPI-VX(J1)*TPI
        XTMP1 = VX(J2)-VX(J1)
        IF(XTMP1 >  180.0_SP)THEN
          XTMP = -360.0_SP*TPI+XTMP
        ELSE IF(XTMP1 < -180.0_SP)THEN
          XTMP =  360.0_SP*TPI+XTMP
        END IF
  
!        PSTX_TM(IA,K)=PSTX_TM(IA,K)-GRAV*DT1(IA)*ELIJ*DLTYC(I)
!        PSTY_TM(IA,K)=PSTY_TM(IA,K)+GRAV*DT1(IA)*ELIJ*XTMP*COS(DEG2RAD*YC(IA)) 
!        PSTX_TM(IB,K)=PSTX_TM(IB,K)+GRAV*DT1(IB)*ELIJ*DLTYC(I)
!        PSTY_TM(IB,K)=PSTY_TM(IB,K)-GRAV*DT1(IB)*ELIJ*XTMP*COS(DEG2RAD*YC(IB)) 
!%        PSTX_TM(IA,K)=PSTX_TM(IA,K)-GRAV_E(IA)*DT1(IA)*DZ1(IA,K)*ELIJ*DLTYC(I)
!%        PSTY_TM(IA,K)=PSTY_TM(IA,K)+GRAV_E(IA)*DT1(IA)*DZ1(IA,K)*ELIJ*XTMP*COS(DEG2RAD*YC(IA)) 
!%        PSTX_TM(IB,K)=PSTX_TM(IB,K)+GRAV_E(IB)*DT1(IB)*DZ1(IB,K)*ELIJ*DLTYC(I)
!%        PSTY_TM(IB,K)=PSTY_TM(IB,K)-GRAV_E(IB)*DT1(IB)*DZ1(IB,K)*ELIJ*XTMP*COS(DEG2RAD*YC(IB)) 

#        if !defined (SEMI_IMPLICIT)
         PSTX_TM(IA,K)=PSTX_TM(IA,K)-F_ALFA(IA)*GRAV_E(IA)*DT1(IA)*DZ1(IA,K)*ELIJ*DLTYC(I)
         PSTY_TM(IA,K)=PSTY_TM(IA,K)+F_ALFA(IA)*GRAV_E(IA)*DT1(IA)*DZ1(IA,K)*ELIJ*XTMP*COS(DEG2RAD*YC(IA))
         PSTX_TM(IB,K)=PSTX_TM(IB,K)+F_ALFA(IB)*GRAV_E(IB)*DT1(IB)*DZ1(IB,K)*ELIJ*DLTYC(I)
         PSTY_TM(IB,K)=PSTY_TM(IB,K)-F_ALFA(IB)*GRAV_E(IB)*DT1(IB)*DZ1(IB,K)*ELIJ*XTMP*COS(DEG2RAD*YC(IB))
#        else

         IF(STAGE<KSTAGE_UV) THEN
           PSTX_TM(IA,K)=PSTX_TM(IA,K)-F_ALFA(IA)*GRAV_E(IA)*DT1(IA)*DZ1(IA,K)*ELIJ*DLTYC(I)
           PSTY_TM(IA,K)=PSTY_TM(IA,K)+F_ALFA(IA)*GRAV_E(IA)*DT1(IA)*DZ1(IA,K)*ELIJ*XTMP*COS(DEG2RAD*YC(IA))
           PSTX_TM(IB,K)=PSTX_TM(IB,K)+F_ALFA(IB)*GRAV_E(IB)*DT1(IB)*DZ1(IB,K)*ELIJ*DLTYC(I)
           PSTY_TM(IB,K)=PSTY_TM(IB,K)-F_ALFA(IB)*GRAV_E(IB)*DT1(IB)*DZ1(IB,K)*ELIJ*XTMP*COS(DEG2RAD*YC(IB))
         ELSE
           PSTX_TM(IA,K)=PSTX_TM(IA,K)-GRAV_E(IA)*DT1(IA)*DZ1(IA,K)*ELIJ*DLTYC(I)
           PSTY_TM(IA,K)=PSTY_TM(IA,K)+GRAV_E(IA)*DT1(IA)*DZ1(IA,K)*ELIJ*XTMP*COS(DEG2RAD*YC(IA))
           PSTX_TM(IB,K)=PSTX_TM(IB,K)+GRAV_E(IB)*DT1(IB)*DZ1(IB,K)*ELIJ*DLTYC(I)
           PSTY_TM(IB,K)=PSTY_TM(IB,K)-GRAV_E(IB)*DT1(IB)*DZ1(IB,K)*ELIJ*XTMP*COS(DEG2RAD*YC(IB))
         ENDIF

#        endif

#    else
!        PSTX_TM(IA,K)=PSTX_TM(IA,K)-GRAV*DT1(IA)*ELIJ*DLTYC(I)
!        PSTY_TM(IA,K)=PSTY_TM(IA,K)+GRAV*DT1(IA)*ELIJ*DLTXC(I)
!        PSTX_TM(IB,K)=PSTX_TM(IB,K)+GRAV*DT1(IB)*ELIJ*DLTYC(I)
!        PSTY_TM(IB,K)=PSTY_TM(IB,K)-GRAV*DT1(IB)*ELIJ*DLTXC(I)
        PSTX_TM(IA,K)=PSTX_TM(IA,K)-GRAV_E(IA)*DT1(IA)*DZ1(IA,K)*ELIJ*DLTYC(I)
        PSTY_TM(IA,K)=PSTY_TM(IA,K)+GRAV_E(IA)*DT1(IA)*DZ1(IA,K)*ELIJ*DLTXC(I)
        PSTX_TM(IB,K)=PSTX_TM(IB,K)+GRAV_E(IB)*DT1(IB)*DZ1(IB,K)*ELIJ*DLTYC(I)
        PSTY_TM(IB,K)=PSTY_TM(IB,K)-GRAV_E(IB)*DT1(IB)*DZ1(IB,K)*ELIJ*DLTXC(I)
#    endif     

     END DO
   END DO

   DEALLOCATE(UIJ1,VIJ1,UIJ2,VIJ2,STAT=ERROR)
   IF(ERROR /= 0) &
   & CALL FATAL_ERROR("Unexpected deallocation error for UIJ1,VIJ1,UIJ2 and VIJ2.")
   DEALLOCATE(UALFA,VALFA,STAT=ERROR)
   IF(ERROR /= 0) &
   & CALL FATAL_ERROR("Unexpected deallocation error for UALFA,VALFA.")
   DEALLOCATE(FXX,FYY,STAT=ERROR)
   IF(ERROR /= 0) &
   & CALL FATAL_ERROR("Unexpected deallocation error for FXX,FYY.")

#  else

   DO I=1,NE
     IA=IEC(I,1)
     IB=IEC(I,2)
     J1=IENODE(I,1)
     J2=IENODE(I,2)

# if !defined (SEMI_IMPLICIT)

     ELIJ=0.5_SP*(EGF(J1)+EGF(J2))

# if defined (AIR_PRESSURE)
     ELIJ=ELIJ-0.5_SP*(EGF_AIR(J1)+EGF_AIR(J2))
# endif

#    if defined (EQUI_TIDE)
     ELIJ=ELIJ-0.5_SP*(EGF_EQI(J1)+EGF_EQI(J2))
#    endif
#    if defined (ATMO_TIDE)
     ELIJ=ELIJ-0.5_SP*(EGF_ATMO(J1)+EGF_ATMO(J2))
#    endif       

# else

     IF(STAGE<KSTAGE_UV) THEN
       ELIJ=0.5_SP*(ET(J1)+ET(J2))
     ELSE
       ELIJ=(1.0_SP-IFCETA)*0.5_SP*(ET(J1)+ET(J2))
     ENDIF
#    if defined (AIR_PRESSURE)
     ELIJ=ELIJ-( (1.0_SP-IFCETA)*0.5_SP*(EL_AIR(J1)+EL_AIR(J2))+IFCETA*0.5_SP*(ELF_AIR(J1)+ELF_AIR(J2)) ) 
#    endif
#    if defined (EQUI_TIDE)
     ELIJ=ELIJ-( (1.0_SP-IFCETA)*0.5_SP*(EL_EQI(J1)+EL_EQI(J2))+IFCETA*0.5_SP*(ELF_EQI(J1)+ELF_EQI(J2)) )
#    endif
#    if defined (ATMO_TIDE)
     ELIJ=ELIJ-( (1.0_SP-IFCETA)*0.5_SP*(EL_ATMO(J1)+EL_ATMO(J2))+IFCETA*0.5_SP*(ELF_ATMO(J1)+ELF_ATMO(J2)) )
#    endif

# endif

     K1=NBE(IA,1)
     K2=NBE(IA,2)
     K3=NBE(IA,3)
     K4=NBE(IB,1)
     K5=NBE(IB,2)
     K6=NBE(IB,3)
#    if defined (SPHERICAL)
     XIJA=DLTXNE(I,1)
     YIJA=DLTYNE(I,1)
     XIJB=DLTXNE(I,2)
     YIJB=DLTYNE(I,2)
#    else
     XIJA=XIJC(I)-XC(IA)
     YIJA=YIJC(I)-YC(IA)
     XIJB=XIJC(I)-XC(IB)
     YIJB=YIJC(I)-YC(IB)
#    endif

     UK1 = U(K1,:)
     UK2 = U(K2,:)
     UK3 = U(K3,:)
     UK4 = U(K4,:)
     UK5 = U(K5,:)
     UK6 = U(K6,:)
     VK1 = V(K1,:)
     VK2 = V(K2,:)
     VK3 = V(K3,:)
     VK4 = V(K4,:)
     VK5 = V(K5,:)
     VK6 = V(K6,:)
     IF(K1 == 0) CALL GHOSTUV3(IA,1,UK1,VK1)
     IF(K2 == 0) CALL GHOSTUV3(IA,2,UK2,VK2)
     IF(K3 == 0) CALL GHOSTUV3(IA,3,UK3,VK3)
     IF(K4 == 0) CALL GHOSTUV3(IB,1,UK4,VK4)
     IF(K5 == 0) CALL GHOSTUV3(IB,2,UK5,VK5)
     IF(K6 == 0) CALL GHOSTUV3(IB,3,UK6,VK6)

     DO K=1,KBM1
       DIJ=0.5_SP*(DT(J1)*DZ(J1,K)+DT(J2)*DZ(J2,K))
#      if defined (WET_DRY)
#      if !defined (SEMI_IMPLICIT)
       IF(ISWETCT(IA)*ISWETC(IA) == 1 .OR. ISWETCT(IB)*ISWETC(IB) == 1)THEN
#      else
       IF(ISWETCT(IA) == 1 .OR. ISWETCT(IB) == 1)THEN
#      endif
#      endif

       COFA1=A1U(IA,1)*U(IA,K)+A1U(IA,2)*UK1(K)+A1U(IA,3)*UK2(K)+A1U(IA,4)*UK3(K)
       COFA2=A2U(IA,1)*U(IA,K)+A2U(IA,2)*UK1(K)+A2U(IA,3)*UK2(K)+A2U(IA,4)*UK3(K)
       COFA5=A1U(IA,1)*V(IA,K)+A1U(IA,2)*VK1(K)+A1U(IA,3)*VK2(K)+A1U(IA,4)*VK3(K)
       COFA6=A2U(IA,1)*V(IA,K)+A2U(IA,2)*VK1(K)+A2U(IA,3)*VK2(K)+A2U(IA,4)*VK3(K)
!       COFA1=A1U(IA,1)*U(IA,K)+A1U(IA,2)*U(K1,K)+A1U(IA,3)*U(K2,K)+A1U(IA,4)*U(K3,K)
!       COFA2=A2U(IA,1)*U(IA,K)+A2U(IA,2)*U(K1,K)+A2U(IA,3)*U(K2,K)+A2U(IA,4)*U(K3,K)
!       COFA5=A1U(IA,1)*V(IA,K)+A1U(IA,2)*V(K1,K)+A1U(IA,3)*V(K2,K)+A1U(IA,4)*V(K3,K)
!       COFA6=A2U(IA,1)*V(IA,K)+A2U(IA,2)*V(K1,K)+A2U(IA,3)*V(K2,K)+A2U(IA,4)*V(K3,K)

       UIJ1=U(IA,K)+COFA1*XIJA+COFA2*YIJA
       VIJ1=V(IA,K)+COFA5*XIJA+COFA6*YIJA

       COFA3=A1U(IB,1)*U(IB,K)+A1U(IB,2)*UK4(K)+A1U(IB,3)*UK5(K)+A1U(IB,4)*UK6(K)
       COFA4=A2U(IB,1)*U(IB,K)+A2U(IB,2)*UK4(K)+A2U(IB,3)*UK5(K)+A2U(IB,4)*UK6(K)
       COFA7=A1U(IB,1)*V(IB,K)+A1U(IB,2)*VK4(K)+A1U(IB,3)*VK5(K)+A1U(IB,4)*VK6(K)
       COFA8=A2U(IB,1)*V(IB,K)+A2U(IB,2)*VK4(K)+A2U(IB,3)*VK5(K)+A2U(IB,4)*VK6(K)
!       COFA3=A1U(IB,1)*U(IB,K)+A1U(IB,2)*U(K4,K)+A1U(IB,3)*U(K5,K)+A1U(IB,4)*U(K6,K)
!       COFA4=A2U(IB,1)*U(IB,K)+A2U(IB,2)*U(K4,K)+A2U(IB,3)*U(K5,K)+A2U(IB,4)*U(K6,K)
!       COFA7=A1U(IB,1)*V(IB,K)+A1U(IB,2)*V(K4,K)+A1U(IB,3)*V(K5,K)+A1U(IB,4)*V(K6,K)
!       COFA8=A2U(IB,1)*V(IB,K)+A2U(IB,2)*V(K4,K)+A2U(IB,3)*V(K5,K)+A2U(IB,4)*V(K6,K)

       UIJ2=U(IB,K)+COFA3*XIJB+COFA4*YIJB
       VIJ2=V(IB,K)+COFA7*XIJB+COFA8*YIJB

       UIJ=0.5_SP*(UIJ1+UIJ2)
       VIJ=0.5_SP*(VIJ1+VIJ2)
       EXFLUX = DIJ*(-UIJ*DLTYC(I) + VIJ*DLTXC(I))

#   if defined (PLBC)
       EXFLUX = DIJ*(-UIJ*DLTYC(I) + 0.0_SP*DLTXC(I))
#   endif

!
!-------ADD THE VISCOUS TERM & ADVECTION TERM---------------------------------!
!

       VISCOF1=ART(IA)*SQRT(COFA1**2+COFA6**2+0.5_SP*(COFA2+COFA5)**2)
       VISCOF2=ART(IB)*SQRT(COFA3**2+COFA8**2+0.5_SP*(COFA4+COFA7)**2)

!       VISCOF=FACT*0.5_SP*HORCON*(VISCOF1+VISCOF2)/HPRNU + FM1*HORCON
!       VISCOF=FACT*0.5_SP*HORCON*(VISCOF1+VISCOF2) + FM1*HORCON

       ! David moved HPRNU and added HVC
       VISCOF=(FACT*0.5_SP*(VISCOF1*CC_HVC(IA)+VISCOF2*CC_HVC(IB)) + FM1*0.5_SP*(CC_HVC(IA)+CC_HVC(IB)))/HPRNU
       VISCOFM(IA,K) = VISCOFM(IA,K) + VISCOF
       VISCOFM(IB,K) = VISCOFM(IB,K) + VISCOF

       TXXIJ=(COFA1+COFA3)*VISCOF
       TYYIJ=(COFA6+COFA8)*VISCOF
       TXYIJ=0.5_SP*(COFA2+COFA4+COFA5+COFA7)*VISCOF
       FXX=DIJ*(TXXIJ*DLTYC(I)-TXYIJ*DLTXC(I))
       FYY=DIJ*(TXYIJ*DLTYC(I)-TYYIJ*DLTXC(I))

#    if defined (PLBC)
       TXYIJ=0.5_SP*(COFA5+COFA7)*VISCOF
       FXX=DIJ*(TXXIJ*DLTYC(I)-0.0_SP*DLTXC(I))
       FYY=DIJ*(TXYIJ*DLTYC(I)-0.0_SP*DLTXC(I))
#    endif

       XADV=EXFLUX*((1.0_SP-SIGN(1.0_SP,EXFLUX))*UIJ2+(1.0_SP+SIGN(1.0_SP,EXFLUX))*UIJ1)*0.5_SP
       YADV=EXFLUX*((1.0_SP-SIGN(1.0_SP,EXFLUX))*VIJ2+(1.0_SP+SIGN(1.0_SP,EXFLUX))*VIJ1)*0.5_SP

       !!CALCULATE BOUNDARY FLUX AUGMENTERS
       TPA = FLOAT(1-ISBC(I))*EPOR(IA)
       TPB = FLOAT(1-ISBC(I))*EPOR(IB)


       !!ACCUMULATE ADVECTIVE + DIFFUSIVE + BAROTROPIC PRESSURE GRADIENT TERMS
       XFLUX(IA,K)=XFLUX(IA,K)+XADV*TPA+FXX*TPA
       YFLUX(IA,K)=YFLUX(IA,K)+YADV*TPA+FYY*TPA
       XFLUX(IB,K)=XFLUX(IB,K)-XADV*TPB-FXX*TPB
       YFLUX(IB,K)=YFLUX(IB,K)-YADV*TPB-FYY*TPB

#  if defined (WET_DRY)
    END IF
#  endif
! for spherical coordinator and domain across 360^o latitude
#    if defined (SPHERICAL)
        XTMP  = VX(J2)*TPI-VX(J1)*TPI
        XTMP1 = VX(J2)-VX(J1)
        IF(XTMP1 >  180.0_SP)THEN
	  XTMP = -360.0_SP*TPI+XTMP
        ELSE IF(XTMP1 < -180.0_SP)THEN
	  XTMP =  360.0_SP*TPI+XTMP
	END IF
	  
!        PSTX_TM(IA,K)=PSTX_TM(IA,K)-GRAV*DT1(IA)*ELIJ*DLTYC(I)
!        PSTY_TM(IA,K)=PSTY_TM(IA,K)+GRAV*DT1(IA)*ELIJ*XTMP*COS(DEG2RAD*YC(IA)) 
!        PSTX_TM(IB,K)=PSTX_TM(IB,K)+GRAV*DT1(IB)*ELIJ*DLTYC(I)
!        PSTY_TM(IB,K)=PSTY_TM(IB,K)-GRAV*DT1(IB)*ELIJ*XTMP*COS(DEG2RAD*YC(IB)) 
!%        PSTX_TM(IA,K)=PSTX_TM(IA,K)-GRAV_E(IA)*DT1(IA)*DZ1(IA,K)*ELIJ*DLTYC(I)
!%        PSTY_TM(IA,K)=PSTY_TM(IA,K)+GRAV_E(IA)*DT1(IA)*DZ1(IA,K)*ELIJ*XTMP*COS(DEG2RAD*YC(IA)) 
!%        PSTX_TM(IB,K)=PSTX_TM(IB,K)+GRAV_E(IB)*DT1(IB)*DZ1(IB,K)*ELIJ*DLTYC(I)
!%        PSTY_TM(IB,K)=PSTY_TM(IB,K)-GRAV_E(IB)*DT1(IB)*DZ1(IB,K)*ELIJ*XTMP*COS(DEG2RAD*YC(IB)) 

#        if !defined (SEMI_IMPLICIT)
         PSTX_TM(IA,K)=PSTX_TM(IA,K)-F_ALFA(IA)*GRAV_E(IA)*DT1(IA)*DZ1(IA,K)*ELIJ*DLTYC(I)
         PSTY_TM(IA,K)=PSTY_TM(IA,K)+F_ALFA(IA)*GRAV_E(IA)*DT1(IA)*DZ1(IA,K)*ELIJ*XTMP*COS(DEG2RAD*YC(IA))
         PSTX_TM(IB,K)=PSTX_TM(IB,K)+F_ALFA(IB)*GRAV_E(IB)*DT1(IB)*DZ1(IB,K)*ELIJ*DLTYC(I)
         PSTY_TM(IB,K)=PSTY_TM(IB,K)-F_ALFA(IB)*GRAV_E(IB)*DT1(IB)*DZ1(IB,K)*ELIJ*XTMP*COS(DEG2RAD*YC(IB))
#        else

         IF(STAGE<KSTAGE_UV) THEN
           PSTX_TM(IA,K)=PSTX_TM(IA,K)-F_ALFA(IA)*GRAV_E(IA)*DT1(IA)*DZ1(IA,K)*ELIJ*DLTYC(I)
           PSTY_TM(IA,K)=PSTY_TM(IA,K)+F_ALFA(IA)*GRAV_E(IA)*DT1(IA)*DZ1(IA,K)*ELIJ*XTMP*COS(DEG2RAD*YC(IA))
           PSTX_TM(IB,K)=PSTX_TM(IB,K)+F_ALFA(IB)*GRAV_E(IB)*DT1(IB)*DZ1(IB,K)*ELIJ*DLTYC(I)
           PSTY_TM(IB,K)=PSTY_TM(IB,K)-F_ALFA(IB)*GRAV_E(IB)*DT1(IB)*DZ1(IB,K)*ELIJ*XTMP*COS(DEG2RAD*YC(IB))
         ELSE
           PSTX_TM(IA,K)=PSTX_TM(IA,K)-GRAV_E(IA)*DT1(IA)*DZ1(IA,K)*ELIJ*DLTYC(I)
           PSTY_TM(IA,K)=PSTY_TM(IA,K)+GRAV_E(IA)*DT1(IA)*DZ1(IA,K)*ELIJ*XTMP*COS(DEG2RAD*YC(IA))
           PSTX_TM(IB,K)=PSTX_TM(IB,K)+GRAV_E(IB)*DT1(IB)*DZ1(IB,K)*ELIJ*DLTYC(I)
           PSTY_TM(IB,K)=PSTY_TM(IB,K)-GRAV_E(IB)*DT1(IB)*DZ1(IB,K)*ELIJ*XTMP*COS(DEG2RAD*YC(IB))
         ENDIF

#        endif

#    else
!        PSTX_TM(IA,K)=PSTX_TM(IA,K)-GRAV*DT1(IA)*ELIJ*DLTYC(I)
!        PSTY_TM(IA,K)=PSTY_TM(IA,K)+GRAV*DT1(IA)*ELIJ*DLTXC(I)
!        PSTX_TM(IB,K)=PSTX_TM(IB,K)+GRAV*DT1(IB)*ELIJ*DLTYC(I)
!        PSTY_TM(IB,K)=PSTY_TM(IB,K)-GRAV*DT1(IB)*ELIJ*DLTXC(I)
        PSTX_TM(IA,K)=PSTX_TM(IA,K)-GRAV_E(IA)*DT1(IA)*DZ1(IA,K)*ELIJ*DLTYC(I)
        PSTY_TM(IA,K)=PSTY_TM(IA,K)+GRAV_E(IA)*DT1(IA)*DZ1(IA,K)*ELIJ*DLTXC(I)
        PSTX_TM(IB,K)=PSTX_TM(IB,K)+GRAV_E(IB)*DT1(IB)*DZ1(IB,K)*ELIJ*DLTYC(I)
        PSTY_TM(IB,K)=PSTY_TM(IB,K)-GRAV_E(IB)*DT1(IB)*DZ1(IB,K)*ELIJ*DLTXC(I)
#    endif     

     END DO
   END DO
#  endif

# if defined (PLBC)
    PSTY_TM=0.0_SP
# endif

#     if !defined (ONE_D_MODEL) 
      DO I=1,N
#     if defined (WET_DRY)
#       if !defined (SEMI_IMPLICIT)
        ISWETTMP = ISWETCT(I)*ISWETC(I)
#       else
        ISWETTMP = ISWETCT(I)
#       endif
        DO K=1,KBM1
         XFLUX(I,K)  = XFLUX(I,K)*ISWETTMP
         YFLUX(I,K)  = YFLUX(I,K)*ISWETTMP
         PSTX_TM(I,K)= PSTX_TM(I,K)*ISWETTMP
         PSTY_TM(I,K)= PSTY_TM(I,K)*ISWETTMP
        END DO
#     endif
       DO K=1,KBM1
        XFLUX(I,K)=XFLUX(I,K)+PSTX_TM(I,K)
        YFLUX(I,K)=YFLUX(I,K)+PSTY_TM(I,K)
       END DO
      END DO
#     else
        CALL ONEDTIDE(PSTX_TM,PSTY_TM)
       DO I=1,N
         DO K=1,KBM1
            XFLUX(I,K)=PSTX_TM(I,K)*ART(I)
            YFLUX(I,K)=PSTY_TM(I,K)*ART(I)
         END DO
       END DO
#     endif

!
!-------ADD VERTICAL CONVECTIVE FLUX, CORIOLIS TERM AND BAROCLINIC PG TERM----!
!
   DO I=1,N
#    if defined (WET_DRY)
#    if !defined (SEMI_IMPLICIT)
     IF(ISWETCT(I)*ISWETC(I) == 1)THEN
#    else
     IF(ISWETCT(I) == 1)THEN
#    endif
#    endif

#    if defined (NH) || defined (LIMITER_VER_ADV)
     U_TEMP(0)  = -U(I,1)
     V_TEMP(0)  = -V(I,1)
     U_TEMP(KB) = -U(I,KBM1)
     V_TEMP(KB) = -V(I,KBM1)
     SL_H(0)    = DZ1(I,1)
     SL_H(KB)   = DZ1(I,KBM1)
     DO K=1, KBM1
       U_TEMP(K) = U(I,K)
       V_TEMP(K) = V(I,K)
       SL_H(K)   = DZ1(I,K)
     ENDDO

     DO K=2, KBM1
#    if defined (WAVE_CURRENT_INTERACTION) && (VORTEX_FORCE)
       CONV_U(K) = (W(I,K)+W_STOKES_3D(I,K))*(U_TEMP(K)+U_TEMP(K-1))*0.5_SP    !Xia et al. (2020)
#    else
       CONV_U(K) = W(I,K)*(U_TEMP(K)+U_TEMP(K-1))*0.5_SP
#    endif
       SL_U = 2.0_SP*(U_TEMP(K)-U_TEMP(K+1))/(SL_H(K)+SL_H(K+1))
       SL_F = 2.0_SP*(U_TEMP(K-2)-U_TEMP(K-1))/(SL_H(K-2)+SL_H(K-1))
       DISS_U(K) = 0.5_SP*ABS(W(I,K))*(U_TEMP(K-1)-U_TEMP(K)-0.5_SP*LIMLED2(SL_U,SL_F,2.0_SP)*(SL_H(K-1)+SL_H(K)))

#    if defined (WAVE_CURRENT_INTERACTION) && (VORTEX_FORCE)
       CONV_V(K) = (W(I,K)+W_STOKES_3D(I,K))*(V_TEMP(K)+V_TEMP(K-1))*0.5_SP
#    else
       CONV_V(K) = W(I,K)*(V_TEMP(K)+V_TEMP(K-1))*0.5_SP
#    endif
       SL_U = 2.0_SP*(V_TEMP(K)-V_TEMP(K+1))/(SL_H(K)+SL_H(K+1))
       SL_F = 2.0_SP*(V_TEMP(K-2)-V_TEMP(K-1))/(SL_H(K-2)+SL_H(K-1))
       DISS_V(K) = 0.5_SP*ABS(W(I,K))*(V_TEMP(K-1)-V_TEMP(K)-0.5_SP*LIMLED2(SL_U,SL_F,2.0_SP)*(SL_H(K-1)+SL_H(K)))
     ENDDO
     CONV_U(1)  = 0.0_SP
     DISS_U(1)  = 0.0_SP
     CONV_U(KB) = 0.0_SP
     DISS_U(KB) = 0.0_SP

     CONV_V(1)  = 0.0_SP
     DISS_V(1)  = 0.0_SP
     CONV_V(KB) = 0.0_SP
     DISS_V(KB) = 0.0_SP
#    endif

     DO K=1,KBM1

#      if defined (NH) || defined (LIMITER_VER_ADV)
       XFLUXV = CONV_U(K)-CONV_U(K+1)+DISS_U(K+1)-DISS_U(K)
       YFLUXV = CONV_V(K)-CONV_V(K+1)+DISS_V(K+1)-DISS_V(K)
#      endif
#      if !defined (NH) && !defined (LIMITER_VER_ADV)
#      if !defined (WAVE_CURRENT_INTERACTION)
       IF(K == 1) THEN
         XFLUXV=-W(I,K+1)*(U(I,K)*DZ1(I,K+1)+U(I,K+1)*DZ1(I,K))/&
                 (DZ1(I,K)+DZ1(I,K+1))
         YFLUXV=-W(I,K+1)*(V(I,K)*DZ1(I,K+1)+V(I,K+1)*DZ1(I,K))/&
                 (DZ1(I,K)+DZ1(I,K+1))
       ELSE IF(K == KBM1) THEN
         XFLUXV= W(I,K)*(U(I,K)*DZ1(I,K-1)+U(I,K-1)*DZ1(I,K))/&
                 (DZ1(I,K)+DZ1(I,K-1))
         YFLUXV= W(I,K)*(V(I,K)*DZ1(I,K-1)+V(I,K-1)*DZ1(I,K))/&
                 (DZ1(I,K)+DZ1(I,K-1))
      ELSE
         XFLUXV= W(I,K)*(U(I,K)*DZ1(I,K-1)+U(I,K-1)*DZ1(I,K))/&
                 (DZ1(I,K)+DZ1(I,K-1))-&
                 W(I,K+1)*(U(I,K)*DZ1(I,K+1)+U(I,K+1)*DZ1(I,K))/&
                 (DZ1(I,K)+DZ1(I,K+1))
         YFLUXV= W(I,K)*(V(I,K)*DZ1(I,K-1)+V(I,K-1)*DZ1(I,K))/&
                 (DZ1(I,K)+DZ1(I,K-1))-&
                 W(I,K+1)*(V(I,K)*DZ1(I,K+1)+V(I,K+1)*DZ1(I,K))/&
                 (DZ1(I,K)+DZ1(I,K+1))
       END IF
#      else
!START CHANGE Xia et al. (2020)
#      if  !defined (VORTEX_FORCE)
       IF(K == 1) THEN
         XFLUXV=-W(I,K+1)*(U(I,K)*DZ1(I,K+1)+U(I,K+1)*DZ1(I,K))/&
                 (DZ1(I,K)+DZ1(I,K+1))
         YFLUXV=-W(I,K+1)*(V(I,K)*DZ1(I,K+1)+V(I,K+1)*DZ1(I,K))/&
                 (DZ1(I,K)+DZ1(I,K+1))
       ELSE IF(K == KBM1) THEN
         XFLUXV= W(I,K)*(U(I,K)*DZ1(I,K-1)+U(I,K-1)*DZ1(I,K))/&
                 (DZ1(I,K)+DZ1(I,K-1))
         YFLUXV= W(I,K)*(V(I,K)*DZ1(I,K-1)+V(I,K-1)*DZ1(I,K))/&
                 (DZ1(I,K)+DZ1(I,K-1))
      ELSE
         XFLUXV= W(I,K)*(U(I,K)*DZ1(I,K-1)+U(I,K-1)*DZ1(I,K))/&
                 (DZ1(I,K)+DZ1(I,K-1))-&
                 W(I,K+1)*(U(I,K)*DZ1(I,K+1)+U(I,K+1)*DZ1(I,K))/&
                 (DZ1(I,K)+DZ1(I,K+1))
         YFLUXV= W(I,K)*(V(I,K)*DZ1(I,K-1)+V(I,K-1)*DZ1(I,K))/&
                 (DZ1(I,K)+DZ1(I,K-1))-&
                 W(I,K+1)*(V(I,K)*DZ1(I,K+1)+V(I,K+1)*DZ1(I,K))/&
                 (DZ1(I,K)+DZ1(I,K+1))
       END IF
#      else
       IF(K == 1) THEN
         XFLUXV=-(W(I,K+1)+W_STOKES_3D(I,K+1))*(U(I,K)*DZ1(I,K+1)+U(I,K+1)*DZ1(I,K))/&
                 (DZ1(I,K)+DZ1(I,K+1))
         YFLUXV=-(W(I,K+1)+W_STOKES_3D(I,K+1))*(V(I,K)*DZ1(I,K+1)+V(I,K+1)*DZ1(I,K))/&
                 (DZ1(I,K)+DZ1(I,K+1))
       ELSE IF(K == KBM1) THEN
         XFLUXV= (W(I,K)+W_STOKES_3D(I,K))*(U(I,K)*DZ1(I,K-1)+U(I,K-1)*DZ1(I,K))/&
                 (DZ1(I,K)+DZ1(I,K-1))
         YFLUXV= (W(I,K)+W_STOKES_3D(I,K))*(V(I,K)*DZ1(I,K-1)+V(I,K-1)*DZ1(I,K))/&
                 (DZ1(I,K)+DZ1(I,K-1))
       ELSE
         XFLUXV= (W(I,K)+W_STOKES_3D(I,K))*(U(I,K)*DZ1(I,K-1)+U(I,K-1)*DZ1(I,K))/&
                 (DZ1(I,K)+DZ1(I,K-1))-&
                 (W(I,K+1)+W_STOKES_3D(I,K+1))*(U(I,K)*DZ1(I,K+1)+U(I,K+1)*DZ1(I,K))/&
                 (DZ1(I,K)+DZ1(I,K+1))
         YFLUXV= (W(I,K)+W_STOKES_3D(I,K))*(V(I,K)*DZ1(I,K-1)+V(I,K-1)*DZ1(I,K))/&
                 (DZ1(I,K)+DZ1(I,K-1))-&
                 (W(I,K+1)+W_STOKES_3D(I,K+1))*(V(I,K)*DZ1(I,K+1)+V(I,K+1)*DZ1(I,K))/&
                 (DZ1(I,K)+DZ1(I,K+1))
       END IF
!END CHANGE Xia et al. (2020)
#      endif
#      endif
#      endif
#      if defined (SPHERICAL)
!       XFLUX(I,K)=XFLUX(I,K)+XFLUXV/DZ(K)*ART(I)&
!                 +DRHOX(I,K)-COR(I)*V(I,K)*DT1(I)*ART(I)&
!                 -U(I,K)*V(I,K)/REARTH*TAN(YC(I)*DEG2RAD)*DT1(I)*ART(I)&
!                 +0.5_SP*U(I,K)*(W(I,K+1)+W(I,K))/REARTH*DT1(I)*ART(I)
!       YFLUX(I,K)=YFLUX(I,K)+YFLUXV/DZ(K)*ART(I)&
!                 +DRHOY(I,K)+COR(I)*U(I,K)*DT1(I)*ART(I)&
!                 +U(I,K)*U(I,K)/REARTH*TAN(YC(I)*DEG2RAD)*DT1(I)*ART(I)&
!                 +0.5_SP*V(I,K)*(W(I,K+1)+W(I,K))/REARTH*DT1(I)*ART(I)
#      if !defined (WAVE_CURRENT_INTERACTION) && !defined (VORTEX_FORCE)
       XFLUX(I,K)=XFLUX(I,K)+XFLUXV*ART(I)&
                 +DRHOX(I,K)-COR(I)*V(I,K)*DT1(I)*DZ1(I,K)*ART(I)&
                 -U(I,K)*V(I,K)/REARTH*TAN(YC(I)*DEG2RAD)*DT1(I)*DZ1(I,K)*ART(I)&
                 +0.5_SP*U(I,K)*(W(I,K+1)+W(I,K))/REARTH*DT1(I)*DZ1(I,K)*ART(I)
       YFLUX(I,K)=YFLUX(I,K)+YFLUXV*ART(I)&
                 +DRHOY(I,K)+COR(I)*U(I,K)*DT1(I)*DZ1(I,K)*ART(I)&
                 +U(I,K)*U(I,K)/REARTH*TAN(YC(I)*DEG2RAD)*DT1(I)*DZ1(I,K)*ART(I)&
                 +0.5_SP*V(I,K)*(W(I,K+1)+W(I,K))/REARTH*DT1(I)*DZ1(I,K)*ART(I)
#      else
       XFLUX(I,K)=XFLUX(I,K)+XFLUXV*ART(I)&
                 +DRHOX(I,K)-COR(I)*V(I,K)*DT1(I)*DZ1(I,K)*ART(I)&
                 -U(I,K)*V(I,K)/REARTH*TAN(YC(I)*DEG2RAD)*DT1(I)*DZ1(I,K)*ART(I)&
                 +0.5_SP*U(I,K)*(W(I,K+1)+W(I,K)+W_STOKES_3D(I,K+1)+W_STOKES_3D(I,K))/REARTH*DT1(I)*DZ1(I,K)*ART(I)
       YFLUX(I,K)=YFLUX(I,K)+YFLUXV*ART(I)&
                 +DRHOY(I,K)+COR(I)*U(I,K)*DT1(I)*DZ1(I,K)*ART(I)&
                 +U(I,K)*U(I,K)/REARTH*TAN(YC(I)*DEG2RAD)*DT1(I)*DZ1(I,K)*ART(I)&
                 +0.5_SP*V(I,K)*(W(I,K+1)+W(I,K)+W_STOKES_3D(I,K+1)+W_STOKES_3D(I,K))/REARTH*DT1(I)*DZ1(I,K)*ART(I)
#      endif
#      else
#      if !defined (NH)   
!       XFLUX(I,K)=XFLUX(I,K)+XFLUXV/DZ(K)*ART(I)&
!                 +DRHOX(I,K)-COR(I)*V(I,K)*DT1(I)*ART(I)
!       YFLUX(I,K)=YFLUX(I,K)+YFLUXV/DZ(K)*ART(I)&
!                 +DRHOY(I,K)+COR(I)*U(I,K)*DT1(I)*ART(I)
       XFLUX(I,K)=XFLUX(I,K)+XFLUXV*ART(I)&
                 +DRHOX(I,K)-COR(I)*V(I,K)*DT1(I)*DZ1(I,K)*ART(I)
       YFLUX(I,K)=YFLUX(I,K)+YFLUXV*ART(I)&
                 +DRHOY(I,K)+COR(I)*U(I,K)*DT1(I)*DZ1(I,K)*ART(I)
#      else
       XFLUX(I,K)=XFLUX(I,K)+XFLUXV*ART(I)&
                 +DRHOX(I,K)-COR(I)*V(I,K)*DT1(I)*DZ1(I,K)*ART(I)+NHQDRX(I,K)
       YFLUX(I,K)=YFLUX(I,K)+YFLUXV*ART(I)&
                 +DRHOY(I,K)+COR(I)*U(I,K)*DT1(I)*DZ1(I,K)*ART(I)+NHQDRY(I,K)
#      endif
#      endif

#  if defined (WAVE_CURRENT_INTERACTION) || (NUOPC)
       XFLUX(I,K) = XFLUX(I,K) + WAVESTRX_3D(I,K)
       YFLUX(I,K) = YFLUX(I,K) + WAVESTRY_3D(I,K)
#  endif   

     END DO
#  if defined (WET_DRY)
    END IF
#  endif
   END DO

#  if defined (SPHERICAL)
#  if !defined (SEMI_IMPLICIT)
     CALL ADV_UV_EDGE_XY(XFLUX,YFLUX,0.0_SP,0,0)
#  else
     CALL ADV_UV_EDGE_XY(XFLUX,YFLUX,IFCETA,STAGE,KSTAGE_UV)
#  endif
#  endif
   
   DO I=1,N
     IF(ISBCE(I) == 2) THEN
       DO K=1,KBM1
#        if !defined (SEMI_IMPLICIT)
         XFLUX(I,K)=0.0_SP
         YFLUX(I,K)=0.0_SP
#        else

#        if !defined (NH)
         XFLUX(I,K)=PSTX_TM(I,K)
         YFLUX(I,K)=PSTY_TM(I,K)
#        else
         XFLUX(I,K)=PSTX_TM(I,K)+NHQDRX(I,K)
         YFLUX(I,K)=PSTY_TM(I,K)+NHQDRY(I,K)
#        endif

#        endif
       END DO
     END IF
   END DO


   !ADJUST FLUX AT RIVER INFLOWS
   IF(NUMQBC >= 1) THEN
     IF(RIVER_INFLOW_LOCATION == 'node') THEN
       DO II=1,NUMQBC
         J=INODEQ(II)
         I1=NBVE(J,1)
         I2=NBVE(J,NTVE(J))
         DO K=1,KBM1
           VLCTYQ(II)=QDIS(II)/QAREA(II)
!           TEMP=0.5_SP*QDIS(II)*VQDIST(II,K)*VLCTYQ(II)
           TEMP=0.5_SP*QDIS(II)*VQDIST(II,K)*VQDIST(II,K)*VLCTYQ(II)/DZ(J,K)
!           XFLUX(I1,K)=XFLUX(I1,K)-TEMP/DZ(J,K)*COS(ANGLEQ(II))
!           XFLUX(I2,K)=XFLUX(I2,K)-TEMP/DZ(J,K)*COS(ANGLEQ(II))
!           YFLUX(I1,K)=YFLUX(I1,K)-TEMP/DZ(J,K)*SIN(ANGLEQ(II))
!           YFLUX(I2,K)=YFLUX(I2,K)-TEMP/DZ(J,K)*SIN(ANGLEQ(II))
           XFLUX(I1,K)=XFLUX(I1,K)-TEMP*COS(ANGLEQ(II))
           XFLUX(I2,K)=XFLUX(I2,K)-TEMP*COS(ANGLEQ(II))
           YFLUX(I1,K)=YFLUX(I1,K)-TEMP*SIN(ANGLEQ(II))
           YFLUX(I2,K)=YFLUX(I2,K)-TEMP*SIN(ANGLEQ(II))
         END DO
       END DO
     ELSE IF(RIVER_INFLOW_LOCATION == 'edge') THEN
       DO II=1,NUMQBC
         I1=ICELLQ(II)
         DO K=1,KBM1
           VLCTYQ(II)=QDIS(II)/QAREA(II)
!           TEMP=QDIS(II)*VQDIST(II,K)*VLCTYQ(II)
           TEMP=QDIS(II)*VQDIST(II,K)*VQDIST(II,K)*VLCTYQ(II)/DZ1(I1,K)
!           XFLUX(I1,K)=XFLUX(I1,K)-TEMP/DZ1(I1,K)*COS(ANGLEQ(II))
!           YFLUX(I1,K)=YFLUX(I1,K)-TEMP/DZ1(I1,K)*SIN(ANGLEQ(II))
           XFLUX(I1,K)=XFLUX(I1,K)-TEMP*COS(ANGLEQ(II))
           YFLUX(I1,K)=YFLUX(I1,K)-TEMP*SIN(ANGLEQ(II))
         END DO
       END DO
     ELSE
       PRINT*,'RIVER_INFLOW_LOCATION NOT CORRECT'
       CALL PSTOP
     END IF
   END IF


   DO I=1,N
      DO K=1,KBM1
#      if !defined (SEMI_IMPLICIT)
         !       UF(I,K)=U(I,K)*DT1(I,K)/D1(I,K)-DTI*XFLUX(I,K)/ART(I)/D1(I,K)
         !       VF(I,K)=V(I,K)*DT1(I,K)/D1(I,K)-DTI*YFLUX(I,K)/ART(I)/D1(I,K)
         UF(I,K)=U(I,K)*DT1(I)/D1(I)-DTI*XFLUX(I,K)/ART(I)/(D1(I)*DZ1(I,K))
         VF(I,K)=V(I,K)*DT1(I)/D1(I)-DTI*YFLUX(I,K)/ART(I)/(D1(I)*DZ1(I,K))
         IF(ADCOR_ON) THEN
         ! FOR ADCOR
           UBETA(I,K)=XFLUX(I,K) +COR(I)*V(I,K)*DT1(I)*DZ1(I,K)*ART(I)*EPOR(I)
           VBETA(I,K)=YFLUX(I,K) -COR(I)*U(I,K)*DT1(I)*DZ1(I,K)*ART(I)*EPOR(I)
         ENDIF
#      else
       IF(STAGE<KSTAGE_UV) THEN
         UF(I,K)=UC(I,K)-RK_UV(STAGE)*DTI*XFLUX(I,K)/ART(I)/(DT1(I)*DZ1(I,K))-BEDF*U_N(I,K)
         VF(I,K)=VC(I,K)-RK_UV(STAGE)*DTI*YFLUX(I,K)/ART(I)/(DT1(I)*DZ1(I,K))-BEDF*V_N(I,K)
       ELSE
         XFLUX3(I,K)=XFLUX(I,K)
         YFLUX3(I,K)=YFLUX(I,K)
         IF(ADCOR_ON) THEN
         ! FOR ADCOR
           UBETA(I,K)=XFLUX(I,K) +COR(I)*V(I,K)*DT1(I)*DZ1(I,K)*ART(I)*EPOR(I)
           VBETA(I,K)=YFLUX(I,K) -COR(I)*U(I,K)*DT1(I)*DZ1(I,K)*ART(I)*EPOR(I)
         ENDIF
       END IF 
#      endif
      END DO
   END DO

#  if defined (SPHERICAL)
  ! IF YOUR PROCESSOR HAS THE NORTH POLE
   DO J =1,NP
      I = NP_LST(J)
      DO K=1,KBM1
#      if !defined (SEMI_IMPLICIT)
         U_TMP = -V(I,K)*COS(XC(I)*DEG2RAD)-U(I,K)*SIN(XC(I)*DEG2RAD)
         V_TMP = -V(I,K)*SIN(XC(I)*DEG2RAD)+U(I,K)*COS(XC(I)*DEG2RAD)
         !       UF_TMP=U_TMP*DT1(I,K)/D1(I,K)-DTI*XFLUX(I,K)/ART(I)/D1(I,K)
         !       VF_TMP=V_TMP*DT1(I,K)/D1(I,K)-DTI*YFLUX(I,K)/ART(I)/D1(I,K)
         UF_TMP=U_TMP*DT1(I)/D1(I)-DTI*XFLUX(I,K)/ART(I)/(D1(I)*DZ1(I,K))
         VF_TMP=V_TMP*DT1(I)/D1(I)-DTI*YFLUX(I,K)/ART(I)/(D1(I)*DZ1(I,K))
         
       IF(ADCOR_ON)THEN
         UF(I,K)  = UF_TMP
         VF(I,K)  = VF_TMP
         UBETA(I,K)=XFLUX(I,K) +COR(I)*V_TMP*DT1(I)*DZ1(I,K)*ART(I)
         VBETA(I,K)=YFLUX(I,K) -COR(I)*U_TMP*DT1(I)*DZ1(I,K)*ART(I)
       ELSE
         UF(I,K)  = VF_TMP*COS(XC(I)*DEG2RAD)-UF_TMP*SIN(XC(I)*DEG2RAD)
         VF(I,K)  = -UF_TMP*COS(XC(I)*DEG2RAD)-VF_TMP*SIN(XC(I)*DEG2RAD)
       END IF 
#      else
       IF(STAGE<KSTAGE_UV) THEN
         U_TMP = -VC(I,K)*COS(XC(I)*DEG2RAD)-UC(I,K)*SIN(XC(I)*DEG2RAD)
         V_TMP = -VC(I,K)*SIN(XC(I)*DEG2RAD)+UC(I,K)*COS(XC(I)*DEG2RAD)

         UF_TMP=U_TMP-RK_UV(STAGE)*DTI*XFLUX3_NP(I,K)/ART(I)/(DT1(I)*DZ1(I,K))
         VF_TMP=V_TMP-RK_UV(STAGE)*DTI*YFLUX3_NP(I,K)/ART(I)/(DT1(I)*DZ1(I,K))
         UF(I,K)  = VF_TMP*COS(XC(I)*DEG2RAD)-UF_TMP*SIN(XC(I)*DEG2RAD)
         VF(I,K)  = UF_TMP*COS(XC(I)*DEG2RAD)+VF_TMP*SIN(XC(I)*DEG2RAD)
         VF(I,K)  = -VF(I,K)
       ELSE
         XFLUX3(I,K)=XFLUX(I,K)
         YFLUX3(I,K)=YFLUX(I,K)
         IF (ADCOR_ON) THEN
           U_TMP = -V(I,K)*COS(XC(I)*DEG2RAD)-U(I,K)*SIN(XC(I)*DEG2RAD)
           V_TMP = -V(I,K)*SIN(XC(I)*DEG2RAD)+U(I,K)*COS(XC(I)*DEG2RAD)

           UBETA(I,K)=XFLUX(I,K) +COR(I)*V(I,K)*DT1(I)*DZ1(I,K)*ART(I)
           VBETA(I,K)=YFLUX(I,K) -COR(I)*U(I,K)*DT1(I)*DZ1(I,K)*ART(I)
           UBETA_NP(I,K)=XFLUX3_NP(I,K) +COR(I)*V_TMP*DT1(I)*DZ1(I,K)*ART(I)
           VBETA_NP(I,K)=YFLUX3_NP(I,K) -COR(I)*U_TMP*DT1(I)*DZ1(I,K)*ART(I)
         ENDIF
       ENDIF
#      endif
         
      END DO
   END DO
#  endif   
   
#  if defined (WET_DRY)
   DO I =1,N
#    if !defined (SEMI_IMPLICIT)
     IF(ISWETCT(I)*ISWETC(I) .NE. 1)THEN
#    else
     IF(ISWETCT(I) .NE. 1)THEN
#    endif
       DO K=1,KBM1
         UF(I,K)=0.0_SP
         VF(I,K)=0.0_SP
         UBETA(I,K)=0.0_SP
         VBETA(I,K)=0.0_SP
       END DO
     END IF
   END DO
#  endif
   
   DO K=1,KB
    VISCOFM(:,K) = VISCOFM(:,K)/ART(:)
   END DO 

   IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "End: adv_uv_edge_gcy"

   END SUBROUTINE ADV_UV_EDGE_GCY
!==============================================================================!
