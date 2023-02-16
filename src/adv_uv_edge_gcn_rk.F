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

   SUBROUTINE ADV_UV_EDGE_GCN_RK(UB,VB)

!==============================================================================!
! this subroutine calculate advective, coriolis, pressure gradient, etc in     !
! x and y momentum equations except vertical diffusion terms for internal mode ! 
!==============================================================================!

   USE ALL_VARS
   USE BCS
   USE MOD_UTILS
   USE MOD_SPHERICAL
   USE MOD_NORTHPOLE
   USE MOD_WD

#  if defined (WAVE_CURRENT_INTERACTION)
   USE MOD_WAVE_CURRENT_INTERACTION
#  endif

#  if defined (THIN_DAM)
   USE MOD_DAM   !Jadon
#  endif

   IMPLICIT NONE
   REAL(SP) :: XFLUX(0:NT,KB),YFLUX(0:NT,KB)
   REAL(SP) :: PSTX_TM(0:NT,KB),PSTY_TM(0:NT,KB)
   REAL(SP) :: COFA1,COFA2,COFA3,COFA4,COFA5,COFA6,COFA7,COFA8
   REAL(SP) :: XADV,YADV,TXXIJ,TYYIJ,TXYIJ
   REAL(SP) :: VISCOF,VISCOF1,VISCOF2,TEMP,TPA,TPB
   REAL(SP) :: XIJA,YIJA,XIJB,YIJB,UIJ,VIJ
   REAL(SP) :: DIJ,ELIJ,TMPA,TMPB,TMP,XFLUXV,YFLUXV
   REAL(SP) :: FACT,FM1,EXFLUX,ISWETTMP
   INTEGER  :: I,IA,IB,J1,J2,K1,K2,K3,K4,K5,K6,K,II,J,I1,I2

#  if defined (SPHERICAL)
   REAL(DP) :: XTMP,XTMP1
   REAL(SP) :: U_TMP,V_TMP,UF_TMP,VF_TMP
#  endif

#  if defined (LIMITED_NO)
   REAL(SP) :: UIJ1,VIJ1,UIJ2,VIJ2,FXX,FYY
#  else
   REAL(SP),ALLOCATABLE,DIMENSION(:) :: UIJ1,VIJ1,UIJ2,VIJ2,FXX,FYY
   REAL(SP),ALLOCATABLE,DIMENSION(:) :: UALFA,VALFA
   REAL(SP) :: UALFA_TMP,VALFA_TMP
   INTEGER :: ERROR
   REAL(SP) :: EPS
#  endif

#  if defined (LIMITER_VER_ADV)
   REAL(SP)  CONV_U(1:KB), CONV_V(1:KB), DISS_U(1:KB), DISS_V(1:KB)
   REAL(SP)  SL_H(0:KB), U_TEMP(0:KB), V_TEMP(0:KB)
   REAL(SP)  SL_U, SL_F
#  endif

!#  if defined (THIN_DAM)
   REAL(SP) :: A1UIA1,A1UIA2,A1UIA3,A1UIA4,A2UIA1,A2UIA2,A2UIA3,A2UIA4
   REAL(SP) :: A1UIB1,A1UIB2,A1UIB3,A1UIB4,A2UIB1,A2UIB2,A2UIB3,A2UIB4   
   INTEGER  :: J11,J12,J21,J22,E1,E2,ISBCE1,ISBC_TMP,IB_TMP
   LOGICAL  :: ISMATCH
!#  endif

   REAL(SP) :: UB(0:NT,KB),VB(0:NT,KB)

   IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "Start: adv_uv_edge_gcn_rk.F"

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
   
   DO K=1,KBM1
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
#      if defined (THIN_DAM)
       IF(IB==0.AND.E_DAM_MATCH(IA)/=0.AND.K<=KDAM1(IA))XIJB=DLTXNE_DAM_MATCH(I)
       IF(IB==0.AND.E_DAM_MATCH(IA)/=0.AND.K<=KDAM1(IA))YIJB=DLTYNE_DAM_MATCH(I)
#      endif
#      else
       XIJA=XIJC(I)-XC(IA)
       YIJA=YIJC(I)-YC(IA)
       XIJB=XIJC(I)-XC(IB)
       YIJB=YIJC(I)-YC(IB)
#      if defined (THIN_DAM)
       IF(IB==0.AND.E_DAM_MATCH(IA)/=0.AND.K<=KDAM1(IA))XIJB=XIJC(I)-XC(E_DAM_MATCH(IA))
       IF(IB==0.AND.E_DAM_MATCH(IA)/=0.AND.K<=KDAM1(IA))YIJB=YIJC(I)-YC(E_DAM_MATCH(IA))
#      endif
#      endif

       DIJ=0.5_SP*(DT(J1)*DZ(J1,K)+DT(J2)*DZ(J2,K))
#      if defined (WET_DRY)
       IF(ISWETCT(IA)*ISWETC(IA) == 1 .OR. ISWETCT(IB)*ISWETC(IB) == 1)THEN
#      endif

       IB_TMP = IB
!----------------------Used for Dam Model By Jadon--------------------
#      if defined (THIN_DAM)
       A1UIA1 = A1U(IA,1)
       A1UIA2 = A1U(IA,2)
       A1UIA3 = A1U(IA,3)
       A1UIA4 = A1U(IA,4)
       A2UIA1 = A2U(IA,1)
       A2UIA2 = A2U(IA,2)
       A2UIA3 = A2U(IA,3)
       A2UIA4 = A2U(IA,4)
       
       A1UIB1 = A1U(IB,1)
       A1UIB2 = A1U(IB,2)
       A1UIB3 = A1U(IB,3)
       A1UIB4 = A1U(IB,4)
       A2UIB1 = A2U(IB,1)
       A2UIB2 = A2U(IB,2)
       A2UIB3 = A2U(IB,3)
       A2UIB4 = A2U(IB,4)
       
       IF(ISBCE(IA) == 1 .AND. K <= KDAM1(IA))THEN
         A1UIA1 = A1U_DAM(IA,1)
         A1UIA2 = A1U_DAM(IA,2)
         A1UIA3 = A1U_DAM(IA,3)
         A1UIA4 = A1U_DAM(IA,4)
         A2UIA1 = A2U_DAM(IA,1)
         A2UIA2 = A2U_DAM(IA,2)
         A2UIA3 = A2U_DAM(IA,3)
         A2UIA4 = A2U_DAM(IA,4)
	 IF(K1 == 0)K1 = NBE_DAM(IA)
	 IF(K2 == 0)K2 = NBE_DAM(IA)
	 IF(K3 == 0)K3 = NBE_DAM(IA)
       END IF

       IF(IB==0.AND.E_DAM_MATCH(IA)/=0.AND.K<=KDAM1(IA))IB_TMP=E_DAM_MATCH(IA)
       IF(ISBCE(IB_TMP) == 1 .AND. K <= KDAM1(IB_TMP))THEN
         K4=NBE(IB_TMP,1)
         K5=NBE(IB_TMP,2)
         K6=NBE(IB_TMP,3)
         A1UIB1 = A1U_DAM(IB_TMP,1)
         A1UIB2 = A1U_DAM(IB_TMP,2)
         A1UIB3 = A1U_DAM(IB_TMP,3)
         A1UIB4 = A1U_DAM(IB_TMP,4)
         A2UIB1 = A2U_DAM(IB_TMP,1)
         A2UIB2 = A2U_DAM(IB_TMP,2)
         A2UIB3 = A2U_DAM(IB_TMP,3)
         A2UIB4 = A2U_DAM(IB_TMP,4)
	 IF(K4 == 0)K4 = NBE_DAM(IB_TMP)
	 IF(K5 == 0)K5 = NBE_DAM(IB_TMP)
	 IF(K6 == 0)K6 = NBE_DAM(IB_TMP)
       END IF
#      else
       A1UIA1 = A1U(IA,1)
       A1UIA2 = A1U(IA,2)
       A1UIA3 = A1U(IA,3)
       A1UIA4 = A1U(IA,4)
       A2UIA1 = A2U(IA,1)
       A2UIA2 = A2U(IA,2)
       A2UIA3 = A2U(IA,3)
       A2UIA4 = A2U(IA,4)
       
       A1UIB1 = A1U(IB_TMP,1)
       A1UIB2 = A1U(IB_TMP,2)
       A1UIB3 = A1U(IB_TMP,3)
       A1UIB4 = A1U(IB_TMP,4)
       A2UIB1 = A2U(IB_TMP,1)
       A2UIB2 = A2U(IB_TMP,2)
       A2UIB3 = A2U(IB_TMP,3)
       A2UIB4 = A2U(IB_TMP,4)
#      endif
!---------------------------------------------------------------------
       COFA1=A1UIA1*U(IA,K)+A1UIA2*U(K1,K)+A1UIA3*U(K2,K)+A1UIA4*U(K3,K)
       COFA2=A2UIA1*U(IA,K)+A2UIA2*U(K1,K)+A2UIA3*U(K2,K)+A2UIA4*U(K3,K)
       COFA5=A1UIA1*V(IA,K)+A1UIA2*V(K1,K)+A1UIA3*V(K2,K)+A1UIA4*V(K3,K)
       COFA6=A2UIA1*V(IA,K)+A2UIA2*V(K1,K)+A2UIA3*V(K2,K)+A2UIA4*V(K3,K)

       UIJ1(I)=COFA1*XIJA+COFA2*YIJA
       VIJ1(I)=COFA5*XIJA+COFA6*YIJA
       UALFA_TMP=ABS(U(IA,K)-U(IB_TMP,K))/ABS(UIJ1(I)+EPSILON(EPS))
       VALFA_TMP=ABS(V(IA,K)-V(IB_TMP,K))/ABS(VIJ1(I)+EPSILON(EPS))
       IF(UALFA_TMP > 1)UALFA_TMP = 1.0_SP
       IF(VALFA_TMP > 1)VALFA_TMP = 1.0_SP
       UALFA(IA)=MIN(UALFA(IA),UALFA_TMP)
       VALFA(IA)=MIN(VALFA(IA),VALFA_TMP)

       COFA3=A1UIB1*U(IB_TMP,K)+A1UIB2*U(K4,K)+A1UIB3*U(K5,K)+A1UIB4*U(K6,K)
       COFA4=A2UIB1*U(IB_TMP,K)+A2UIB2*U(K4,K)+A2UIB3*U(K5,K)+A2UIB4*U(K6,K)
       COFA7=A1UIB1*V(IB_TMP,K)+A1UIB2*V(K4,K)+A1UIB3*V(K5,K)+A1UIB4*V(K6,K)
       COFA8=A2UIB1*V(IB_TMP,K)+A2UIB2*V(K4,K)+A2UIB3*V(K5,K)+A2UIB4*V(K6,K)

       UIJ2(I)=COFA3*XIJB+COFA4*YIJB
       VIJ2(I)=COFA7*XIJB+COFA8*YIJB
       UALFA_TMP=ABS(U(IA,K)-U(IB_TMP,K))/ABS(UIJ2(I)+EPSILON(EPS))
       VALFA_TMP=ABS(V(IA,K)-V(IB_TMP,K))/ABS(VIJ2(I)+EPSILON(EPS))
       IF(UALFA_TMP > 1)UALFA_TMP = 1.0_SP
       IF(VALFA_TMP > 1)VALFA_TMP = 1.0_SP
       UALFA(IB_TMP)=MIN(UALFA(IB_TMP),UALFA_TMP)
       VALFA(IB_TMP)=MIN(VALFA(IB_TMP),VALFA_TMP)
       
!
!-------ADD THE VISCOUS TERM & ADVECTION TERM---------------------------------!
!

       VISCOF1=ART(IA)*SQRT(COFA1**2+COFA6**2+0.5_SP*(COFA2+COFA5)**2)
       VISCOF2=ART(IB_TMP)*SQRT(COFA3**2+COFA8**2+0.5_SP*(COFA4+COFA7)**2)

       ! David moved HPRNU and added HVC
       VISCOF=(FACT*0.5_SP*(VISCOF1*CC_HVC(IA)+VISCOF2*CC_HVC(IB_TMP)) + FM1*0.5_SP*(CC_HVC(IA)+CC_HVC(IB_TMP)))/HPRNU
       VISCOFM(IA,K)     = VISCOFM(IA,K)     + VISCOF
       VISCOFM(IB_TMP,K) = VISCOFM(IB_TMP,K) + VISCOF

       TXXIJ=(COFA1+COFA3)*VISCOF
       TYYIJ=(COFA6+COFA8)*VISCOF
       TXYIJ=0.5_SP*(COFA2+COFA4+COFA5+COFA7)*VISCOF
       FXX(I)=DIJ*(TXXIJ*DLTYC(I)-TXYIJ*DLTXC(I))
       FYY(I)=DIJ*(TXYIJ*DLTYC(I)-TYYIJ*DLTXC(I))
#      if defined (WET_DRY)
       END IF
#      endif
     END DO

     DO I=1,NE
       IA=IEC(I,1)
       IB=IEC(I,2)

       J1=IENODE(I,1)
       J2=IENODE(I,2)
       IB_TMP = IB

#      if defined (THIN_DAM)
       IF(IB==0.AND.E_DAM_MATCH(IA)/=0.AND.K<=KDAM1(IA))IB_TMP=E_DAM_MATCH(IA) 
#      endif

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

       DIJ=0.5_SP*(DT(J1)*DZ(J1,K)+DT(J2)*DZ(J2,K))
#      if defined (WET_DRY)
       IF(ISWETCT(IA)*ISWETC(IA) == 1 .OR. ISWETCT(IB)*ISWETC(IB) == 1)THEN
#      endif       
       UIJ1(I)=U(IA,K)+UALFA(IA)*UIJ1(I)
       VIJ1(I)=V(IA,K)+VALFA(IA)*VIJ1(I)
       UIJ2(I)=U(IB_TMP,K)+UALFA(IB_TMP)*UIJ2(I)
       VIJ2(I)=V(IB_TMP,K)+VALFA(IB_TMP)*VIJ2(I)

#      if defined (LIMITED_1)
       IF(UIJ1(I) > MAX(U(IA,K),U(IB_TMP,K)) .OR. UIJ1(I) < MIN(U(IA,K),U(IB_TMP,K)) .OR. &
          UIJ2(I) > MAX(U(IA,K),U(IB_TMP,K)) .OR. UIJ2(I) < MIN(U(IA,K),U(IB_TMP,K)))THEN
         UIJ1(I)=U(IA,K)
         UIJ2(I)=U(IB_TMP,K)
       END IF

       IF(VIJ1(I) > MAX(V(IA,K),V(IB_TMP,K)) .OR. VIJ1(I) < MIN(V(IA,K),V(IB_TMP,K)) .OR. &
          VIJ2(I) > MAX(V(IA,K),V(IB_TMP,K)) .OR. VIJ2(I) < MIN(V(IA,K),V(IB_TMP,K)))THEN
         VIJ1(I)=V(IA,K)
         VIJ2(I)=V(IB_TMP,K)
       END IF
#      endif

!      NORMAL VELOCITY              
       UIJ=0.5_SP*(UIJ1(I)+UIJ2(I))
       VIJ=0.5_SP*(VIJ1(I)+VIJ2(I))
       EXFLUX = DIJ*(-UIJ*DLTYC(I) + VIJ*DLTXC(I))

       XADV=EXFLUX*((1.0_SP-SIGN(1.0_SP,EXFLUX))*UIJ2(I)+(1.0_SP+SIGN(1.0_SP,EXFLUX))*UIJ1(I))*0.5_SP
       YADV=EXFLUX*((1.0_SP-SIGN(1.0_SP,EXFLUX))*VIJ2(I)+(1.0_SP+SIGN(1.0_SP,EXFLUX))*VIJ1(I))*0.5_SP

       !!CALCULATE BOUNDARY FLUX AUGMENTERS
#  if defined (THIN_DAM)
       IF(IB==0.AND.E_DAM_MATCH(IA)/=0.AND.K <= KDAM1(IA))THEN
          ISBC_TMP = 0
       ELSE      
          ISBC_TMP = ISBC(I)
       ENDIF
#  else
       ISBC_TMP = ISBC(I)
#  endif
       TPA = FLOAT(1-ISBC_TMP)*EPOR(IA)
       TPB = FLOAT(1-ISBC_TMP)*EPOR(IB_TMP)

       !!ACCUMULATE ADVECTIVE + DIFFUSIVE + BAROTROPIC PRESSURE GRADIENT TERMS
!       XFLUX(IA,K)=XFLUX(IA,K)+XADV*TPA+FXX*TPA
!       YFLUX(IA,K)=YFLUX(IA,K)+YADV*TPA+FYY*TPA
!       XFLUX(IB,K)=XFLUX(IB,K)-XADV*TPB-FXX*TPB
!       YFLUX(IB,K)=YFLUX(IB,K)-YADV*TPB-FYY*TPB
       XFLUX(IA,K)=XFLUX(IA,K)+XADV*TPA+(FXX(I)+3.0_SP*FXX(I)*FLOAT(ISBC_TMP))*EPOR(IA)
       YFLUX(IA,K)=YFLUX(IA,K)+YADV*TPA+(FYY(I)+3.0_SP*FYY(I)*FLOAT(ISBC_TMP))*EPOR(IA)
       XFLUX(IB,K)=XFLUX(IB,K)-XADV*TPB-(FXX(I)+3.0_SP*FXX(I)*FLOAT(ISBC_TMP))*EPOR(IB)
       YFLUX(IB,K)=YFLUX(IB,K)-YADV*TPB-(FYY(I)+3.0_SP*FYY(I)*FLOAT(ISBC_TMP))*EPOR(IB)

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
  
         PSTX_TM(IA,K)=PSTX_TM(IA,K)-F_ALFA(IA)*GRAV_E(IA)*DT1(IA)*DZ1(IA,K)*ELIJ*DLTYC(I)
         PSTY_TM(IA,K)=PSTY_TM(IA,K)+F_ALFA(IA)*GRAV_E(IA)*DT1(IA)*DZ1(IA,K)*ELIJ*XTMP*COS(DEG2RAD*YC(IA))
         PSTX_TM(IB,K)=PSTX_TM(IB,K)+F_ALFA(IB)*GRAV_E(IB)*DT1(IB)*DZ1(IB,K)*ELIJ*DLTYC(I)
         PSTY_TM(IB,K)=PSTY_TM(IB,K)-F_ALFA(IB)*GRAV_E(IB)*DT1(IB)*DZ1(IB,K)*ELIJ*XTMP*COS(DEG2RAD*YC(IB))

#    else
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

     K1=NBE(IA,1)
     K2=NBE(IA,2)
     K3=NBE(IA,3)
     K4=NBE(IB,1)
     K5=NBE(IB,2)
     K6=NBE(IB,3)

     DO K=1,KBM1
       DIJ=0.5_SP*(DT(J1)*DZ(J1,K)+DT(J2)*DZ(J2,K))
#    if defined (WET_DRY)
      IF(ISWETCT(IA)*ISWETC(IA) == 1 .OR. ISWETCT(IB)*ISWETC(IB) == 1)THEN
#    endif

#     if defined (SPHERICAL)
       XIJA=DLTXNE(I,1)
       YIJA=DLTYNE(I,1)
       XIJB=DLTXNE(I,2)
       YIJB=DLTYNE(I,2)
#     if defined (THIN_DAM)
       IF(IB==0.AND.E_DAM_MATCH(IA)/=0.AND.K<=KDAM1(IA))XIJB=DLTXNE_DAM_MATCH(I)
       IF(IB==0.AND.E_DAM_MATCH(IA)/=0.AND.K<=KDAM1(IA))YIJB=DLTYNE_DAM_MATCH(I)
#     endif
#     else
       XIJA=XIJC(I)-XC(IA)
       YIJA=YIJC(I)-YC(IA)
       XIJB=XIJC(I)-XC(IB)
       YIJB=YIJC(I)-YC(IB)
#     if defined (THIN_DAM)
       IF(IB==0.AND.E_DAM_MATCH(IA)/=0.AND.K<=KDAM1(IA))XIJB=XIJC(I)-XC(E_DAM_MATCH(IA))
       IF(IB==0.AND.E_DAM_MATCH(IA)/=0.AND.K<=KDAM1(IA))YIJB=YIJC(I)-YC(E_DAM_MATCH(IA))
#     endif
#     endif

       IB_TMP = IB
!----------------------Used for Dam Model By Jadon--------------------
#      if defined (THIN_DAM)
       A1UIA1 = A1U(IA,1)
       A1UIA2 = A1U(IA,2)
       A1UIA3 = A1U(IA,3)
       A1UIA4 = A1U(IA,4)
       A2UIA1 = A2U(IA,1)
       A2UIA2 = A2U(IA,2)
       A2UIA3 = A2U(IA,3)
       A2UIA4 = A2U(IA,4)
        
       A1UIB1 = A1U(IB,1)
       A1UIB2 = A1U(IB,2)
       A1UIB3 = A1U(IB,3)
       A1UIB4 = A1U(IB,4)
       A2UIB1 = A2U(IB,1)
       A2UIB2 = A2U(IB,2)
       A2UIB3 = A2U(IB,3)
       A2UIB4 = A2U(IB,4)
       
       IF(ISBCE(IA) == 1 .AND. K <= KDAM1(IA))THEN
         A1UIA1 = A1U_DAM(IA,1)
         A1UIA2 = A1U_DAM(IA,2)
         A1UIA3 = A1U_DAM(IA,3)
         A1UIA4 = A1U_DAM(IA,4)
         A2UIA1 = A2U_DAM(IA,1)
         A2UIA2 = A2U_DAM(IA,2)
         A2UIA3 = A2U_DAM(IA,3)
         A2UIA4 = A2U_DAM(IA,4)
         IF(K1 == 0)K1 = NBE_DAM(IA)
	 IF(K2 == 0)K2 = NBE_DAM(IA)
	 IF(K3 == 0)K3 = NBE_DAM(IA)
       END IF
       
       IF(IB==0.AND.E_DAM_MATCH(IA)/=0.AND.K<=KDAM1(IA))IB_TMP=E_DAM_MATCH(IA) 
       IF(ISBCE(IB_TMP) == 1 .AND. K <= KDAM1(IB_TMP))THEN
         K4=NBE(IB_TMP,1)
         K5=NBE(IB_TMP,2)
         K6=NBE(IB_TMP,3)
         A1UIB1 = A1U_DAM(IB_TMP,1)
         A1UIB2 = A1U_DAM(IB_TMP,2)
         A1UIB3 = A1U_DAM(IB_TMP,3)
         A1UIB4 = A1U_DAM(IB_TMP,4)
         A2UIB1 = A2U_DAM(IB_TMP,1)
         A2UIB2 = A2U_DAM(IB_TMP,2)
         A2UIB3 = A2U_DAM(IB_TMP,3)
         A2UIB4 = A2U_DAM(IB_TMP,4)
 	 IF(K4 == 0)K4 = NBE_DAM(IB_TMP)
	 IF(K5 == 0)K5 = NBE_DAM(IB_TMP)
	 IF(K6 == 0)K6 = NBE_DAM(IB_TMP)
       END IF

#      else
       A1UIA1 = A1U(IA,1)
       A1UIA2 = A1U(IA,2)
       A1UIA3 = A1U(IA,3)
       A1UIA4 = A1U(IA,4)
       A2UIA1 = A2U(IA,1)
       A2UIA2 = A2U(IA,2)
       A2UIA3 = A2U(IA,3)
       A2UIA4 = A2U(IA,4)
        
       A1UIB1 = A1U(IB_TMP,1)
       A1UIB2 = A1U(IB_TMP,2)
       A1UIB3 = A1U(IB_TMP,3)
       A1UIB4 = A1U(IB_TMP,4)
       A2UIB1 = A2U(IB_TMP,1)
       A2UIB2 = A2U(IB_TMP,2)
       A2UIB3 = A2U(IB_TMP,3)
       A2UIB4 = A2U(IB_TMP,4)
#      endif
!---------------------------------------------------------------------

       COFA1=A1UIA1*U(IA,K)+A1UIA2*U(K1,K)+A1UIA3*U(K2,K)+A1UIA4*U(K3,K)
       COFA2=A2UIA1*U(IA,K)+A2UIA2*U(K1,K)+A2UIA3*U(K2,K)+A2UIA4*U(K3,K)
       COFA5=A1UIA1*V(IA,K)+A1UIA2*V(K1,K)+A1UIA3*V(K2,K)+A1UIA4*V(K3,K)
       COFA6=A2UIA1*V(IA,K)+A2UIA2*V(K1,K)+A2UIA3*V(K2,K)+A2UIA4*V(K3,K)

       UIJ1=U(IA,K)+COFA1*XIJA+COFA2*YIJA
       VIJ1=V(IA,K)+COFA5*XIJA+COFA6*YIJA

       COFA3=A1UIB1*U(IB_TMP,K)+A1UIB2*U(K4,K)+A1UIB3*U(K5,K)+A1UIB4*U(K6,K)
       COFA4=A2UIB1*U(IB_TMP,K)+A2UIB2*U(K4,K)+A2UIB3*U(K5,K)+A2UIB4*U(K6,K)
       COFA7=A1UIB1*V(IB_TMP,K)+A1UIB2*V(K4,K)+A1UIB3*V(K5,K)+A1UIB4*V(K6,K)
       COFA8=A2UIB1*V(IB_TMP,K)+A2UIB2*V(K4,K)+A2UIB3*V(K5,K)+A2UIB4*V(K6,K)

       UIJ2=U(IB_TMP,K)+COFA3*XIJB+COFA4*YIJB
       VIJ2=V(IB_TMP,K)+COFA7*XIJB+COFA8*YIJB

       UIJ=0.5_SP*(UIJ1+UIJ2)
       VIJ=0.5_SP*(VIJ1+VIJ2)
       EXFLUX = DIJ*(-UIJ*DLTYC(I) + VIJ*DLTXC(I))

!
!-------ADD THE VISCOUS TERM & ADVECTION TERM---------------------------------!
!

       VISCOF1=ART(IA)*SQRT(COFA1**2+COFA6**2+0.5_SP*(COFA2+COFA5)**2)
       VISCOF2=ART(IB_TMP)*SQRT(COFA3**2+COFA8**2+0.5_SP*(COFA4+COFA7)**2)

       ! David moved HPRNU and added HVC
       VISCOF=(FACT*0.5_SP*(VISCOF1*CC_HVC(IA)+VISCOF2*CC_HVC(IB_TMP)) + FM1*0.5_SP*(CC_HVC(IA)+CC_HVC(IB_TMP)))/HPRNU
       VISCOFM(IA,K) = VISCOFM(IA,K) + VISCOF
       VISCOFM(IB_TMP,K) = VISCOFM(IB_TMP,K) + VISCOF

       TXXIJ=(COFA1+COFA3)*VISCOF
       TYYIJ=(COFA6+COFA8)*VISCOF
       TXYIJ=0.5_SP*(COFA2+COFA4+COFA5+COFA7)*VISCOF
       FXX=DIJ*(TXXIJ*DLTYC(I)-TXYIJ*DLTXC(I))
       FYY=DIJ*(TXYIJ*DLTYC(I)-TYYIJ*DLTXC(I))

       XADV=EXFLUX*((1.0_SP-SIGN(1.0_SP,EXFLUX))*UIJ2+(1.0_SP+SIGN(1.0_SP,EXFLUX))*UIJ1)*0.5_SP
       YADV=EXFLUX*((1.0_SP-SIGN(1.0_SP,EXFLUX))*VIJ2+(1.0_SP+SIGN(1.0_SP,EXFLUX))*VIJ1)*0.5_SP

       !!CALCULATE BOUNDARY FLUX AUGMENTERS
#  if defined (THIN_DAM)
       IF(IB==0.AND.E_DAM_MATCH(IA)/=0.AND.K<=KDAM1(IA))THEN
          ISBC_TMP = 0
       ELSE      
          ISBC_TMP = ISBC(I)
       ENDIF
#  else
       ISBC_TMP = ISBC(I)
#  endif
       TPA = FLOAT(1-ISBC_TMP)*EPOR(IA)
       TPB = FLOAT(1-ISBC_TMP)*EPOR(IB_TMP)

       !!ACCUMULATE ADVECTIVE + DIFFUSIVE + BAROTROPIC PRESSURE GRADIENT TERMS

       XFLUX(IA,K)=XFLUX(IA,K)+XADV*TPA+(FXX+3.0_SP*FXX*FLOAT(ISBC_TMP))*EPOR(IA)
       YFLUX(IA,K)=YFLUX(IA,K)+YADV*TPA+(FYY+3.0_SP*FYY*FLOAT(ISBC_TMP))*EPOR(IA)
       XFLUX(IB,K)=XFLUX(IB,K)-XADV*TPB-(FXX+3.0_SP*FXX*FLOAT(ISBC_TMP))*EPOR(IB)
       YFLUX(IB,K)=YFLUX(IB,K)-YADV*TPB-(FYY+3.0_SP*FYY*FLOAT(ISBC_TMP))*EPOR(IB)

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
	  
         PSTX_TM(IA,K)=PSTX_TM(IA,K)-F_ALFA(IA)*GRAV_E(IA)*DT1(IA)*DZ1(IA,K)*ELIJ*DLTYC(I)
         PSTY_TM(IA,K)=PSTY_TM(IA,K)+F_ALFA(IA)*GRAV_E(IA)*DT1(IA)*DZ1(IA,K)*ELIJ*XTMP*COS(DEG2RAD*YC(IA))
         PSTX_TM(IB,K)=PSTX_TM(IB,K)+F_ALFA(IB)*GRAV_E(IB)*DT1(IB)*DZ1(IB,K)*ELIJ*DLTYC(I)
         PSTY_TM(IB,K)=PSTY_TM(IB,K)-F_ALFA(IB)*GRAV_E(IB)*DT1(IB)*DZ1(IB,K)*ELIJ*XTMP*COS(DEG2RAD*YC(IB))
#    else
        PSTX_TM(IA,K)=PSTX_TM(IA,K)-GRAV_E(IA)*DT1(IA)*DZ1(IA,K)*ELIJ*DLTYC(I)
        PSTY_TM(IA,K)=PSTY_TM(IA,K)+GRAV_E(IA)*DT1(IA)*DZ1(IA,K)*ELIJ*DLTXC(I)
        PSTX_TM(IB,K)=PSTX_TM(IB,K)+GRAV_E(IB)*DT1(IB)*DZ1(IB,K)*ELIJ*DLTYC(I)
        PSTY_TM(IB,K)=PSTY_TM(IB,K)-GRAV_E(IB)*DT1(IB)*DZ1(IB,K)*ELIJ*DLTXC(I)
#    endif     

     END DO
   END DO

#  endif

      DO I=1,N
#     if defined (WET_DRY)
       ISWETTMP = ISWETCT(I)*ISWETC(I)
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

!
!-------ADD VERTICAL CONVECTIVE FLUX, CORIOLIS TERM AND BAROCLINIC PG TERM----!
!
   DO I=1,N
#    if defined (WET_DRY)
     IF(ISWETCT(I)*ISWETC(I) == 1)THEN
#    endif

#    if defined (LIMITER_VER_ADV)
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
       CONV_U(K) = W(I,K)*(U_TEMP(K)+U_TEMP(K-1))*0.5_SP
       SL_U = 2.0_SP*(U_TEMP(K)-U_TEMP(K+1))/(SL_H(K)+SL_H(K+1))
       SL_F = 2.0_SP*(U_TEMP(K-2)-U_TEMP(K-1))/(SL_H(K-2)+SL_H(K-1))
       DISS_U(K) = 0.5_SP*ABS(W(I,K))*(U_TEMP(K-1)-U_TEMP(K)-0.5_SP*LIMLED2(SL_U,SL_F,2.0_SP)*(SL_H(K-1)+SL_H(K)))

       CONV_V(K) = W(I,K)*(V_TEMP(K)+V_TEMP(K-1))*0.5_SP
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

#      if defined (LIMITER_VER_ADV)
       XFLUXV = CONV_U(K)-CONV_U(K+1)+DISS_U(K+1)-DISS_U(K)
       YFLUXV = CONV_V(K)-CONV_V(K+1)+DISS_V(K+1)-DISS_V(K)
#      endif

#      if !defined (LIMITER_VER_ADV)
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
#      endif
#      if defined (SPHERICAL)
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
                 +DRHOX(I,K)-COR(I)*V(I,K)*DT1(I)*DZ1(I,K)*ART(I)
       YFLUX(I,K)=YFLUX(I,K)+YFLUXV*ART(I)&
                 +DRHOY(I,K)+COR(I)*U(I,K)*DT1(I)*DZ1(I,K)*ART(I)
#      endif


#  if defined (WAVE_CURRENT_INTERACTION)
       XFLUX(I,K) = XFLUX(I,K) + WAVESTRX_3D(I,K)
       YFLUX(I,K) = YFLUX(I,K) + WAVESTRY_3D(I,K)
#  endif   

     END DO
#  if defined (WET_DRY)
    END IF
#  endif
   END DO

#  if defined (SPHERICAL)
     CALL ADV_UV_EDGE_XY(XFLUX,YFLUX,0.0_SP,0,0)
#  endif

      DO I=1,N
         IF(ISBCE(I) == 2) THEN
            DO K=1,KBM1
               XFLUX(I,K)=0.0_SP
               YFLUX(I,K)=0.0_SP
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
           TEMP=0.5_SP*QDIS(II)*VQDIST(II,K)*VQDIST(II,K)*VLCTYQ(II)/DZ(J,K)
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
           TEMP=QDIS(II)*VQDIST(II,K)*VQDIST(II,K)*VLCTYQ(II)/DZ1(I1,K)
           XFLUX(I1,K)=XFLUX(I1,K)-TEMP*COS(ANGLEQ(II))
           YFLUX(I1,K)=YFLUX(I1,K)-TEMP*SIN(ANGLEQ(II))
         END DO
       END DO
     ELSE
       PRINT*,'RIVER_INFLOW_LOCATION NOT CORRECT'
       CALL PSTOP
     END IF
   END IF

   DO I =1,N
     DO K=1,KBM1
       UF(I,K)=UB(I,K)*DT1(I)/D1(I)-DTI*XFLUX(I,K)/ART(I)/(D1(I)*DZ1(I,K))
       VF(I,K)=VB(I,K)*DT1(I)/D1(I)-DTI*YFLUX(I,K)/ART(I)/(D1(I)*DZ1(I,K))
       IF(ADCOR_ON) THEN
         UBETA(I,K)=XFLUX(I,K) +COR(I)*VB(I,K)*DT1(I)*DZ1(I,K)*ART(I)*EPOR(I)
         VBETA(I,K)=YFLUX(I,K) -COR(I)*UB(I,K)*DT1(I)*DZ1(I,K)*ART(I)*EPOR(I)
       ENDIF

     END DO
  END DO

#  if defined (SPHERICAL)
  ! IF YOUR PROCESSOR HAS THE NORTH POLE
  DO J =1,NP
     I = NP_LST(J)
     DO K=1,KBM1
        U_TMP = -VB(I,K)*COS(XC(I)*DEG2RAD)-UB(I,K)*SIN(XC(I)*DEG2RAD)
        V_TMP = -VB(I,K)*SIN(XC(I)*DEG2RAD)+UB(I,K)*COS(XC(I)*DEG2RAD)
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

     END DO
  END DO
#  endif   


#  if defined (WET_DRY)
   DO I =1,N
    IF(ISWETCT(I)*ISWETC(I) .NE. 1)THEN
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

   IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "End: adv_uv_edge_gcn_rk"

   END SUBROUTINE ADV_UV_EDGE_GCN_RK
!==============================================================================!
