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
!   CALCULATE CONVECTION AND DIFFUSION FLUXES FOR EXTERNAL MODE                !
!   Ghost cell boundary conditions are used in here
!==============================================================================|
   SUBROUTINE ADVAVE_EDGE_GCY(XFLUX,YFLUX)
!==============================================================================|

   USE ALL_VARS
   USE MOD_UTILS
   USE MOD_SPHERICAL
   USE MOD_NORTHPOLE
   USE BCS
   USE MOD_OBCS
   USE MOD_WD

#  if defined (SEMI_IMPLICIT)
   USE MOD_SEMI_IMPLICIT
#  endif

   IMPLICIT NONE
   INTEGER  :: I,J,K,IA,IB,J1,J2,K1,K2,K3,I1,I2
   REAL(SP) :: DIJ,ELIJ,XIJ,YIJ,UIJ,VIJ
   REAL(SP) :: COFA1,COFA2,COFA3,COFA4,COFA5,COFA6,COFA7,COFA8
   REAL(SP) :: XADV,YADV,TXXIJ,TYYIJ,TXYIJ,UN
   REAL(SP) :: VISCOF,VISCOF1,VISCOF2,TEMP
   REAL(SP) :: XFLUX(0:NT),YFLUX(0:NT)
   REAL(SP) :: FACT,FM1,ISWETTMP
   REAL(SP) :: UAK1,UAK2,UAK3,VAK1,VAK2,VAK3

#  if defined (SPHERICAL)
   REAL(DP) :: XTMP,XTMP1
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

   REAL(SP) :: BTPS
   REAL(SP) :: U_TMP,V_TMP,UAC_TMP,VAC_TMP,WUSURF_TMP,WVSURF_TMP,WUBOT_TMP,WVBOT_TMP,UAF_TMP,VAF_TMP

  if(dbg_set(dbg_sbr)) write(ipt,*) "Start: advave_gcy.F"

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
!-------------------------INITIALIZE FLUXES------------------------------------!
!
   XFLUX = 0.0_SP
   YFLUX = 0.0_SP
   PSTX  = 0.0_SP
   PSTY  = 0.0_SP

!
!-------------------------ACCUMULATE FLUX OVER ELEMENT EDGES-------------------!
!
#  if !defined (LIMITED_NO)
   ALLOCATE(UIJ1(NE),VIJ1(NE),UIJ2(NE),VIJ2(NE),STAT=ERROR)
   IF(ERROR /= 0) &
   & CALL FATAL_ERROR("The arrays UIJ1,VIJ1,UIJ2 and VIJ2 can not be allocated.")
   UIJ1=0.0_SP;VIJ1=0.0_SP;UIJ2=0.0_SP;VIJ2=0.0_SP
   
   ALLOCATE(UALFA(0:NT),VALFA(0:NT),STAT=ERROR)
   IF(ERROR /= 0) &
   & CALL FATAL_ERROR("The arrays UALFA,VALFA can not be allocated.")
   UALFA=1.0_SP;VALFA=1.0_SP
   
   ALLOCATE(FXX(NE),FYY(NE),STAT=ERROR)
   IF(ERROR /= 0) &
   & CALL FATAL_ERROR("The arrays FXX,FYY can not be allocated.")
   FXX=0.0_SP;FYY=0.0_SP

   DO I=1,NE
     IA=IEC(I,1)
     IB=IEC(I,2)
     J1=IENODE(I,1)
     J2=IENODE(I,2)

#    if !defined (SEMI_IMPLICIT)
     DIJ=0.5_SP*(D(J1)+D(J2))
#    else
     DIJ=0.5_SP*(DT(J1)+DT(J2))
#    endif

#    if defined (WET_DRY)
#    if !defined (SEMI_IMPLICIT)
     IF(ISWETCE(IA)*ISWETC(IA) == 1 .OR. ISWETCE(IB)*ISWETC(IB) == 1)THEN
#    else
     IF(ISWETCT(IA) == 1 .OR. ISWETCT(IB) == 1)THEN
#    endif
#    endif
!    FLUX FROM LEFT
     K1=NBE(IA,1)
     K2=NBE(IA,2)
     K3=NBE(IA,3)

     UAK1 = UA(K1)
     UAK2 = UA(K2)
     UAK3 = UA(K3)
     VAK1 = VA(K1)
     VAK2 = VA(K2)
     VAK3 = VA(K3)

     IF(K1 == 0) CALL GHOSTUV2(IA,1,UAK1,VAK1)
     IF(K2 == 0) CALL GHOSTUV2(IA,2,UAK2,VAK2)
     IF(K3 == 0) CALL GHOSTUV2(IA,3,UAK3,VAK3)

     COFA1=A1U(IA,1)*UA(IA)+A1U(IA,2)*UAK1+A1U(IA,3)*UAK2+A1U(IA,4)*UAK3
     COFA2=A2U(IA,1)*UA(IA)+A2U(IA,2)*UAK1+A2U(IA,3)*UAK2+A2U(IA,4)*UAK3
     COFA5=A1U(IA,1)*VA(IA)+A1U(IA,2)*VAK1+A1U(IA,3)*VAK2+A1U(IA,4)*VAK3
     COFA6=A2U(IA,1)*VA(IA)+A2U(IA,2)*VAK1+A2U(IA,3)*VAK2+A2U(IA,4)*VAK3

#    if defined (SPHERICAL)
     UIJ1(I)=COFA1*DLTXNE(I,1)+COFA2*DLTYNE(I,1)
     VIJ1(I)=COFA5*DLTXNE(I,1)+COFA6*DLTYNE(I,1)
#    else
     XIJ=XIJC(I)-XC(IA)
     YIJ=YIJC(I)-YC(IA)
     UIJ1(I)=COFA1*XIJ+COFA2*YIJ
     VIJ1(I)=COFA5*XIJ+COFA6*YIJ
#    endif
     UALFA_TMP=ABS(UA(IA)-UA(IB))/ABS(UIJ1(I)+EPSILON(EPS))
     VALFA_TMP=ABS(VA(IA)-VA(IB))/ABS(VIJ1(I)+EPSILON(EPS))
     IF(UALFA_TMP > 1)UALFA_TMP = 1.0_SP
     IF(VALFA_TMP > 1)VALFA_TMP = 1.0_SP
     UALFA(IA)=MIN(UALFA(IA),UALFA_TMP)
     VALFA(IA)=MIN(VALFA(IA),VALFA_TMP)

!    FLUX FROM RIGHT
     K1=NBE(IB,1)
     K2=NBE(IB,2)
     K3=NBE(IB,3)

     UAK1 = UA(K1)
     UAK2 = UA(K2)
     UAK3 = UA(K3)
     VAK1 = VA(K1)
     VAK2 = VA(K2)
     VAK3 = VA(K3)

     IF(K1 == 0) CALL GHOSTUV2(IB,1,UAK1,VAK1)
     IF(K2 == 0) CALL GHOSTUV2(IB,2,UAK2,VAK2)
     IF(K3 == 0) CALL GHOSTUV2(IB,3,UAK3,VAK3)

     COFA3=A1U(IB,1)*UA(IB)+A1U(IB,2)*UAK1+A1U(IB,3)*UAK2+A1U(IB,4)*UAK3
     COFA4=A2U(IB,1)*UA(IB)+A2U(IB,2)*UAK1+A2U(IB,3)*UAK2+A2U(IB,4)*UAK3
     COFA7=A1U(IB,1)*VA(IB)+A1U(IB,2)*VAK1+A1U(IB,3)*VAK2+A1U(IB,4)*VAK3
     COFA8=A2U(IB,1)*VA(IB)+A2U(IB,2)*VAK1+A2U(IB,3)*VAK2+A2U(IB,4)*VAK3

#    if defined (SPHERICAL)
     UIJ2(I)=COFA3*DLTXNE(I,2)+COFA4*DLTYNE(I,2)
     VIJ2(I)=COFA7*DLTXNE(I,2)+COFA8*DLTYNE(I,2)
#    else
     XIJ=XIJC(I)-XC(IB)
     YIJ=YIJC(I)-YC(IB)
     UIJ2(I)=COFA3*XIJ+COFA4*YIJ
     VIJ2(I)=COFA7*XIJ+COFA8*YIJ
#    endif
     UALFA_TMP=ABS(UA(IA)-UA(IB))/ABS(UIJ2(I)+EPSILON(EPS))
     VALFA_TMP=ABS(VA(IA)-VA(IB))/ABS(VIJ2(I)+EPSILON(EPS))
     IF(UALFA_TMP > 1)UALFA_TMP = 1.0_SP
     IF(VALFA_TMP > 1)VALFA_TMP = 1.0_SP
     UALFA(IB)=MIN(UALFA(IB),UALFA_TMP)
     VALFA(IB)=MIN(VALFA(IB),VALFA_TMP)

!    VISCOSITY COEFFICIENT
     VISCOF1=ART(IA)*SQRT(COFA1**2+COFA6**2+0.5_SP*(COFA2+COFA5)**2)
     VISCOF2=ART(IB)*SQRT(COFA3**2+COFA8**2+0.5_SP*(COFA4+COFA7)**2)
!     VISCOF=HORCON*(FACT*0.5_SP*(VISCOF1+VISCOF2)/HPRNU + FM1)
!     VISCOF=HORCON*(FACT*0.5_SP*(VISCOF1+VISCOF2) + FM1)/HPRNU
     ! David moved HPRNU and added HVC
     VISCOF=(FACT*0.5_SP*(VISCOF1*CC_HVC(IA)+VISCOF2*CC_HVC(IB)) + FM1*0.5_SP*(CC_HVC(IA)+CC_HVC(IB)))/HPRNU

!    SHEAR STRESSES
     TXXIJ=(COFA1+COFA3)*VISCOF
     TYYIJ=(COFA6+COFA8)*VISCOF
     TXYIJ=0.5_SP*(COFA2+COFA4+COFA5+COFA7)*VISCOF
     FXX(I)=DIJ*(TXXIJ*DLTYC(I)-TXYIJ*DLTXC(I))
     FYY(I)=DIJ*(TXYIJ*DLTYC(I)-TYYIJ*DLTXC(I))

# if defined(PLBC)
     TXYIJ=0.5_SP*(COFA5+COFA7)*VISCOF
     FXX(I)=DIJ*(TXXIJ*DLTYC(I)-0.0_SP*DLTXC(I))
     FYY(I)=DIJ*(TXYIJ*DLTYC(I)-0.0_SP*DLTXC(I))
# endif

#    if defined (WET_DRY)
     ENDIF
#    endif
   END DO
   
   DO I=1,NE
     IA=IEC(I,1)
     IB=IEC(I,2)
     J1=IENODE(I,1)
     J2=IENODE(I,2)

#    if !defined (SEMI_IMPLICIT)
     DIJ=0.5_SP*(D(J1)+D(J2))
     ELIJ=0.5_SP*(EL(J1)+EL(J2))

#    if defined (AIR_PRESSURE)
     ELIJ =ELIJ-0.5_SP*(EL_AIR(J1)+EL_AIR(J2))  !*RAMP
#    endif

#    if defined (EQUI_TIDE)
     ELIJ=ELIJ-0.5_SP*(EL_EQI(J1)+EL_EQI(J2))
#    endif
#    if defined (ATMO_TIDE)
     ELIJ=ELIJ-0.5_SP*(EL_ATMO(J1)+EL_ATMO(J2))
#    endif       

#    else

     DIJ=0.5_SP*(DT(J1)+DT(J2))
     IF(STAGE < KSTAGE_UV) THEN
       ELIJ=0.5_SP*(ET(J1)+ET(J2))
     ELSE
       ELIJ=(1.0_SP-IFCETA)*0.5_SP*(ET(J1)+ET(J2))
     ENDIF
#    if defined (AIR_PRESSURE)
     ELIJ =ELIJ-( (1.0_SP-IFCETA)*0.5_SP*(EL_AIR(J1)+EL_AIR(J2))+IFCETA*0.5_SP*(ELF_AIR(J1)+ELF_AIR(J2)) )  !*RAMP
#    endif
#    if defined (EQUI_TIDE)
     ELIJ=ELIJ-( (1.0_SP-IFCETA)*0.5_SP*(EL_EQI(J1)+EL_EQI(J2))+IFCETA*0.5_SP*(ELF_EQI(J1)+ELF_EQI(J2)) )
#    endif
#    if defined (ATMO_TIDE)
     ELIJ=ELIJ-( (1.0_SP-IFCETA)*0.5_SP*(EL_ATMO(J1)+EL_ATMO(J2))+IFCETA*0.5_SP*(ELF_ATMO(J1)+ELF_ATMO(J2)) )
#    endif

#    endif

     UIJ1(I)=UA(IA)+UALFA(IA)*UIJ1(I)
     VIJ1(I)=VA(IA)+VALFA(IA)*VIJ1(I)
     UIJ2(I)=UA(IB)+UALFA(IB)*UIJ2(I)
     VIJ2(I)=VA(IB)+VALFA(IB)*VIJ2(I)

#    if defined (LIMITED_1)
     IF(UIJ1(I) > MAX(UA(IA),UA(IB)) .OR. UIJ1(I) < MIN(UA(IA),UA(IB)) .OR. &
        UIJ2(I) > MAX(UA(IA),UA(IB)) .OR. UIJ2(I) < MIN(UA(IA),UA(IB)))THEN
       UIJ1(I)=UA(IA)
       UIJ2(I)=UA(IB)
     END IF

     IF(VIJ1(I) > MAX(VA(IA),VA(IB)) .OR. VIJ1(I) < MIN(VA(IA),VA(IB)) .OR. &
        VIJ2(I) > MAX(VA(IA),VA(IB)) .OR. VIJ2(I) < MIN(VA(IA),VA(IB)))THEN
       VIJ1(I)=VA(IA)
       VIJ2(I)=VA(IB)
     END IF
#    endif

!    NORMAL VELOCITY
     UIJ=0.5_SP*(UIJ1(I)+UIJ2(I))
     VIJ=0.5_SP*(VIJ1(I)+VIJ2(I))
     UN=-UIJ*DLTYC(I) + VIJ*DLTXC(I)
     
#    if defined (PLBC)
     UN=-UIJ*DLTYC(I) + 0.0_SP*DLTXC(I)
#    endif


#    if defined (WET_DRY)
#    if !defined (SEMI_IMPLICIT)
     IF(ISWETCE(IA)*ISWETC(IA) == 1 .OR. ISWETCE(IB)*ISWETC(IB) == 1)THEN
#    else
     IF(ISWETCT(IA) == 1 .OR. ISWETCT(IB) == 1)THEN
#    endif
#    endif
!    ADD CONVECTIVE AND VISCOUS FLUXES
     XADV=DIJ*UN*&
          ((1.0_SP-SIGN(1.0_SP,UN))*UIJ2(I)+(1.0_SP+SIGN(1.0_SP,UN))*UIJ1(I))*0.5_SP
     YADV=DIJ*UN* &
          ((1.0_SP-SIGN(1.0_SP,UN))*VIJ2(I)+(1.0_SP+SIGN(1.0_SP,UN))*VIJ1(I))*0.5_SP

!    ACCUMULATE FLUX
!     XFLUX(IA)=XFLUX(IA)+XADV
!     YFLUX(IA)=YFLUX(IA)+YADV
!     XFLUX(IB)=XFLUX(IB)-XADV
!     YFLUX(IB)=YFLUX(IB)-YADV

     XFLUX(IA)=XFLUX(IA)+(XADV+FXX(I)*EPOR(IA))*(1.0_SP-ISBC(I))*IUCP(IA)
     YFLUX(IA)=YFLUX(IA)+(YADV+FYY(I)*EPOR(IA))*(1.0_SP-ISBC(I))*IUCP(IA)
     XFLUX(IB)=XFLUX(IB)-(XADV+FXX(I)*EPOR(IB))*(1.0_SP-ISBC(I))*IUCP(IB)
     YFLUX(IB)=YFLUX(IB)-(YADV+FYY(I)*EPOR(IB))*(1.0_SP-ISBC(I))*IUCP(IB)

#    if defined (WET_DRY)
     END IF
#    endif

!    ACCUMULATE BAROTROPIC FLUX
! for spherical coordinator and domain across 360^o latitude         
#    if defined (SPHERICAL)
     XTMP  = VX(J2)*TPI-VX(J1)*TPI
     XTMP1 = VX(J2)-VX(J1)
     IF(XTMP1 >  180.0_SP)THEN
       XTMP = -360.0_SP*TPI+XTMP
     ELSE IF(XTMP1 < -180.0_SP)THEN
       XTMP =  360.0_SP*TPI+XTMP
     END IF  

#    if !defined (SEMI_IMPLICIT)
!     PSTX(IA)=PSTX(IA)-GRAV_E(IA)*D1(IA)*ELIJ*DLTYC(I)
!     PSTY(IA)=PSTY(IA)+GRAV_E(IA)*D1(IA)*ELIJ*XTMP*COS(DEG2RAD*YC(IA))  
!     PSTX(IB)=PSTX(IB)+GRAV_E(IB)*D1(IB)*ELIJ*DLTYC(I)
!     PSTY(IB)=PSTY(IB)-GRAV_E(IB)*D1(IB)*ELIJ*XTMP*COS(DEG2RAD*YC(IB))  

     PSTX(IA)=PSTX(IA)-F_ALFA(IA)*GRAV_E(IA)*D1(IA)*ELIJ*DLTYC(I)
     PSTY(IA)=PSTY(IA)+F_ALFA(IA)*GRAV_E(IA)*D1(IA)*ELIJ*XTMP*COS(DEG2RAD*YC(IA))
     PSTX(IB)=PSTX(IB)+F_ALFA(IB)*GRAV_E(IB)*D1(IB)*ELIJ*DLTYC(I)
     PSTY(IB)=PSTY(IB)-F_ALFA(IB)*GRAV_E(IB)*D1(IB)*ELIJ*XTMP*COS(DEG2RAD*YC(IB))
#    else

     IF(STAGE<KSTAGE_UV) THEN
       PSTX(IA)=PSTX(IA)-F_ALFA(IA)*GRAV_E(IA)*DT1(IA)*ELIJ*DLTYC(I)
       PSTY(IA)=PSTY(IA)+F_ALFA(IA)*GRAV_E(IA)*DT1(IA)*ELIJ*XTMP*COS(DEG2RAD*YC(IA))
       PSTX(IB)=PSTX(IB)+F_ALFA(IB)*GRAV_E(IB)*DT1(IB)*ELIJ*DLTYC(I)
       PSTY(IB)=PSTY(IB)-F_ALFA(IB)*GRAV_E(IB)*DT1(IB)*ELIJ*XTMP*COS(DEG2RAD*YC(IB))
     ELSE
       PSTX(IA)=PSTX(IA)-GRAV_E(IA)*DT1(IA)*ELIJ*DLTYC(I)
       PSTY(IA)=PSTY(IA)+GRAV_E(IA)*DT1(IA)*ELIJ*XTMP*COS(DEG2RAD*YC(IA))
       PSTX(IB)=PSTX(IB)+GRAV_E(IB)*DT1(IB)*ELIJ*DLTYC(I)
       PSTY(IB)=PSTY(IB)-GRAV_E(IB)*DT1(IB)*ELIJ*XTMP*COS(DEG2RAD*YC(IB))
     ENDIF

#    endif
#    else
#    if !defined (SEMI_IMPLICIT)
     PSTX(IA)=PSTX(IA)-GRAV_E(IA)*D1(IA)*ELIJ*DLTYC(I)
     PSTY(IA)=PSTY(IA)+GRAV_E(IA)*D1(IA)*ELIJ*DLTXC(I)
     PSTX(IB)=PSTX(IB)+GRAV_E(IB)*D1(IB)*ELIJ*DLTYC(I)
     PSTY(IB)=PSTY(IB)-GRAV_E(IB)*D1(IB)*ELIJ*DLTXC(I)
#    else
     PSTX(IA)=PSTX(IA)-GRAV_E(IA)*DT1(IA)*ELIJ*DLTYC(I)
     PSTY(IA)=PSTY(IA)+GRAV_E(IA)*DT1(IA)*ELIJ*DLTXC(I)
     PSTX(IB)=PSTX(IB)+GRAV_E(IB)*DT1(IB)*ELIJ*DLTYC(I)
     PSTY(IB)=PSTY(IB)-GRAV_E(IB)*DT1(IB)*ELIJ*DLTXC(I)
#    endif
#    endif     

   END DO

!lwu for no y gradient
#    if defined (PLBC)
     PSTY = 0.0_SP
#    endif

#  else
   
   DO I=1,NE
     IA=IEC(I,1)
     IB=IEC(I,2)
     J1=IENODE(I,1)
     J2=IENODE(I,2)

#    if !defined (SEMI_IMPLICIT)
     DIJ=0.5_SP*(D(J1)+D(J2))
     ELIJ=0.5_SP*(EL(J1)+EL(J2))

#    if defined (AIR_PRESSURE)
     ELIJ =ELIJ-0.5_SP*(EL_AIR(J1)+EL_AIR(J2))  !*RAMP
#    endif

#    if defined (EQUI_TIDE)
     ELIJ=ELIJ-0.5_SP*(EL_EQI(J1)+EL_EQI(J2))
#    endif
#    if defined (ATMO_TIDE)
     ELIJ=ELIJ-0.5_SP*(EL_ATMO(J1)+EL_ATMO(J2))
#    endif       

#    else

     DIJ=0.5_SP*(DT(J1)+DT(J2))
     IF(STAGE < KSTAGE_UV) THEN
       ELIJ=0.5_SP*(ET(J1)+ET(J2))
     ELSE
       ELIJ=(1.0_SP-IFCETA)*0.5_SP*(ET(J1)+ET(J2))
     ENDIF
#    if defined (AIR_PRESSURE)
     ELIJ =ELIJ-( (1.0_SP-IFCETA)*0.5_SP*(EL_AIR(J1)+EL_AIR(J2))+IFCETA*0.5_SP*(ELF_AIR(J1)+ELF_AIR(J2)) )  !*RAMP
#    endif
#    if defined (EQUI_TIDE)
     ELIJ=ELIJ-( (1.0_SP-IFCETA)*0.5_SP*(EL_EQI(J1)+EL_EQI(J2))+IFCETA*0.5_SP*(ELF_EQI(J1)+ELF_EQI(J2)) )
#    endif
#    if defined (ATMO_TIDE)
     ELIJ=ELIJ-( (1.0_SP-IFCETA)*0.5_SP*(EL_ATMO(J1)+EL_ATMO(J2))+IFCETA*0.5_SP*(ELF_ATMO(J1)+ELF_ATMO(J2)) )
#    endif

#    endif

#    if defined (WET_DRY)
#    if !defined (SEMI_IMPLICIT)
     IF(ISWETCE(IA)*ISWETC(IA) == 1 .OR. ISWETCE(IB)*ISWETC(IB) == 1)THEN
#    else
     IF(ISWETCT(IA) == 1 .OR. ISWETCT(IB) == 1)THEN
#    endif
#    endif
!    FLUX FROM LEFT
     K1=NBE(IA,1)
     K2=NBE(IA,2)
     K3=NBE(IA,3)

     UAK1 = UA(K1)
     UAK2 = UA(K2)
     UAK3 = UA(K3)
     VAK1 = VA(K1)
     VAK2 = VA(K2)
     VAK3 = VA(K3)

     IF(K1 == 0) CALL GHOSTUV2(IA,1,UAK1,VAK1)
     IF(K2 == 0) CALL GHOSTUV2(IA,2,UAK2,VAK2)
     IF(K3 == 0) CALL GHOSTUV2(IA,3,UAK3,VAK3)

     COFA1=A1U(IA,1)*UA(IA)+A1U(IA,2)*UAK1+A1U(IA,3)*UAK2+A1U(IA,4)*UAK3
     COFA2=A2U(IA,1)*UA(IA)+A2U(IA,2)*UAK1+A2U(IA,3)*UAK2+A2U(IA,4)*UAK3
     COFA5=A1U(IA,1)*VA(IA)+A1U(IA,2)*VAK1+A1U(IA,3)*VAK2+A1U(IA,4)*VAK3
     COFA6=A2U(IA,1)*VA(IA)+A2U(IA,2)*VAK1+A2U(IA,3)*VAK2+A2U(IA,4)*VAK3

#    if defined (SPHERICAL)
     UIJ1=UA(IA)+COFA1*DLTXNE(I,1)+COFA2*DLTYNE(I,1)
     VIJ1=VA(IA)+COFA5*DLTXNE(I,1)+COFA6*DLTYNE(I,1)
#    else
     XIJ=XIJC(I)-XC(IA)
     YIJ=YIJC(I)-YC(IA)
     UIJ1=UA(IA)+COFA1*XIJ+COFA2*YIJ
     VIJ1=VA(IA)+COFA5*XIJ+COFA6*YIJ
#    endif

!    FLUX FROM RIGHT
     K1=NBE(IB,1)
     K2=NBE(IB,2)
     K3=NBE(IB,3)

     UAK1 = UA(K1)
     UAK2 = UA(K2)
     UAK3 = UA(K3)
     VAK1 = VA(K1)
     VAK2 = VA(K2)
     VAK3 = VA(K3)

     IF(K1 == 0) CALL GHOSTUV2(IB,1,UAK1,VAK1)
     IF(K2 == 0) CALL GHOSTUV2(IB,2,UAK2,VAK2)
     IF(K3 == 0) CALL GHOSTUV2(IB,3,UAK3,VAK3)

     COFA3=A1U(IB,1)*UA(IB)+A1U(IB,2)*UAK1+A1U(IB,3)*UAK2+A1U(IB,4)*UAK3
     COFA4=A2U(IB,1)*UA(IB)+A2U(IB,2)*UAK1+A2U(IB,3)*UAK2+A2U(IB,4)*UAK3
     COFA7=A1U(IB,1)*VA(IB)+A1U(IB,2)*VAK1+A1U(IB,3)*VAK2+A1U(IB,4)*VAK3
     COFA8=A2U(IB,1)*VA(IB)+A2U(IB,2)*VAK1+A2U(IB,3)*VAK2+A2U(IB,4)*VAK3

#    if defined (SPHERICAL)
     UIJ2=UA(IB)+COFA3*DLTXNE(I,2)+COFA4*DLTYNE(I,2)
     VIJ2=VA(IB)+COFA7*DLTXNE(I,2)+COFA8*DLTYNE(I,2)
#    else
     XIJ=XIJC(I)-XC(IB)
     YIJ=YIJC(I)-YC(IB)
     UIJ2=UA(IB)+COFA3*XIJ+COFA4*YIJ
     VIJ2=VA(IB)+COFA7*XIJ+COFA8*YIJ
#    endif

!    NORMAL VELOCITY
     UIJ=0.5_SP*(UIJ1+UIJ2)
     VIJ=0.5_SP*(VIJ1+VIJ2)
     UN=-UIJ*DLTYC(I) + VIJ*DLTXC(I)

!lwu
#    if defined (PLBC)
     UN=-UIJ*DLTYC(I) + 0.0_SP*DLTXC(I)
#    endif

!    VISCOSITY COEFFICIENT
     VISCOF1=ART(IA)*SQRT(COFA1**2+COFA6**2+0.5_SP*(COFA2+COFA5)**2)
     VISCOF2=ART(IB)*SQRT(COFA3**2+COFA8**2+0.5_SP*(COFA4+COFA7)**2)
!     VISCOF=HORCON*(FACT*0.5_SP*(VISCOF1+VISCOF2)/HPRNU + FM1)
!     VISCOF=HORCON*(FACT*0.5_SP*(VISCOF1+VISCOF2) + FM1)
     ! David moved HPRNU and added HVC
     VISCOF=(FACT*0.5_SP*(VISCOF1*CC_HVC(IA)+VISCOF2*CC_HVC(IB)) + FM1*0.5_SP*(CC_HVC(IA)+CC_HVC(IB)))/HPRNU

!    SHEAR STRESSES
     TXXIJ=(COFA1+COFA3)*VISCOF
     TYYIJ=(COFA6+COFA8)*VISCOF
     TXYIJ=0.5_SP*(COFA2+COFA4+COFA5+COFA7)*VISCOF
     FXX=DIJ*(TXXIJ*DLTYC(I)-TXYIJ*DLTXC(I))
     FYY=DIJ*(TXYIJ*DLTYC(I)-TYYIJ*DLTXC(I))

!lwu
#    if defined (PLBC)
     FXX=DIJ*(TXXIJ*DLTYC(I)-0.0_SP*DLTXC(I))
     FYY=DIJ*(TXYIJ*DLTYC(I)-0.0_SP*DLTXC(I))
#    endif

!    ADD CONVECTIVE AND VISCOUS FLUXES
     XADV=DIJ*UN*&
          ((1.0_SP-SIGN(1.0_SP,UN))*UIJ2+(1.0_SP+SIGN(1.0_SP,UN))*UIJ1)*0.5_SP
     YADV=DIJ*UN* &
          ((1.0_SP-SIGN(1.0_SP,UN))*VIJ2+(1.0_SP+SIGN(1.0_SP,UN))*VIJ1)*0.5_SP

!    ACCUMULATE FLUX
!     XFLUX(IA)=XFLUX(IA)+XADV
!     YFLUX(IA)=YFLUX(IA)+YADV
!     XFLUX(IB)=XFLUX(IB)-XADV
!     YFLUX(IB)=YFLUX(IB)-YADV

     XFLUX(IA)=XFLUX(IA)+(XADV+FXX*EPOR(IA))*(1.0_SP-ISBC(I))*IUCP(IA)
     YFLUX(IA)=YFLUX(IA)+(YADV+FYY*EPOR(IA))*(1.0_SP-ISBC(I))*IUCP(IA)
     XFLUX(IB)=XFLUX(IB)-(XADV+FXX*EPOR(IB))*(1.0_SP-ISBC(I))*IUCP(IB)
     YFLUX(IB)=YFLUX(IB)-(YADV+FYY*EPOR(IB))*(1.0_SP-ISBC(I))*IUCP(IB)

#    if defined (WET_DRY)
     END IF
#    endif


!    ACCUMULATE BAROTROPIC FLUX
! for spherical coordinator and domain across 360^o latitude         
#    if defined (SPHERICAL)
     XTMP  = VX(J2)*TPI-VX(J1)*TPI
     XTMP1 = VX(J2)-VX(J1)
     IF(XTMP1 >  180.0_SP)THEN
       XTMP = -360.0_SP*TPI+XTMP
     ELSE IF(XTMP1 < -180.0_SP)THEN
       XTMP =  360.0_SP*TPI+XTMP
     END IF  

#    if !defined (SEMI_IMPLICIT)
!     PSTX(IA)=PSTX(IA)-GRAV_E(IA)*D1(IA)*ELIJ*DLTYC(I)
!     PSTY(IA)=PSTY(IA)+GRAV_E(IA)*D1(IA)*ELIJ*XTMP*COS(DEG2RAD*YC(IA))  
!     PSTX(IB)=PSTX(IB)+GRAV_E(IB)*D1(IB)*ELIJ*DLTYC(I)
!     PSTY(IB)=PSTY(IB)-GRAV_E(IB)*D1(IB)*ELIJ*XTMP*COS(DEG2RAD*YC(IB))  

     PSTX(IA)=PSTX(IA)-F_ALFA(IA)*GRAV_E(IA)*D1(IA)*ELIJ*DLTYC(I)
     PSTY(IA)=PSTY(IA)+F_ALFA(IA)*GRAV_E(IA)*D1(IA)*ELIJ*XTMP*COS(DEG2RAD*YC(IA))
     PSTX(IB)=PSTX(IB)+F_ALFA(IB)*GRAV_E(IB)*D1(IB)*ELIJ*DLTYC(I)
     PSTY(IB)=PSTY(IB)-F_ALFA(IB)*GRAV_E(IB)*D1(IB)*ELIJ*XTMP*COS(DEG2RAD*YC(IB))
#    else

     IF(STAGE<KSTAGE_UV) THEN
       PSTX(IA)=PSTX(IA)-F_ALFA(IA)*GRAV_E(IA)*DT1(IA)*ELIJ*DLTYC(I)
       PSTY(IA)=PSTY(IA)+F_ALFA(IA)*GRAV_E(IA)*DT1(IA)*ELIJ*XTMP*COS(DEG2RAD*YC(IA))
       PSTX(IB)=PSTX(IB)+F_ALFA(IB)*GRAV_E(IB)*DT1(IB)*ELIJ*DLTYC(I)
       PSTY(IB)=PSTY(IB)-F_ALFA(IB)*GRAV_E(IB)*DT1(IB)*ELIJ*XTMP*COS(DEG2RAD*YC(IB))
     ELSE
       PSTX(IA)=PSTX(IA)-GRAV_E(IA)*DT1(IA)*ELIJ*DLTYC(I)
       PSTY(IA)=PSTY(IA)+GRAV_E(IA)*DT1(IA)*ELIJ*XTMP*COS(DEG2RAD*YC(IA))
       PSTX(IB)=PSTX(IB)+GRAV_E(IB)*DT1(IB)*ELIJ*DLTYC(I)
       PSTY(IB)=PSTY(IB)-GRAV_E(IB)*DT1(IB)*ELIJ*XTMP*COS(DEG2RAD*YC(IB))
     ENDIF

#    endif
#    else
#    if !defined (SEMI_IMPLICIT)
     PSTX(IA)=PSTX(IA)-GRAV_E(IA)*D1(IA)*ELIJ*DLTYC(I)
     PSTY(IA)=PSTY(IA)+GRAV_E(IA)*D1(IA)*ELIJ*DLTXC(I)
     PSTX(IB)=PSTX(IB)+GRAV_E(IB)*D1(IB)*ELIJ*DLTYC(I)
     PSTY(IB)=PSTY(IB)-GRAV_E(IB)*D1(IB)*ELIJ*DLTXC(I)
#    else
     PSTX(IA)=PSTX(IA)-GRAV_E(IA)*DT1(IA)*ELIJ*DLTYC(I)
     PSTY(IA)=PSTY(IA)+GRAV_E(IA)*DT1(IA)*ELIJ*DLTXC(I)
     PSTX(IB)=PSTX(IB)+GRAV_E(IB)*DT1(IB)*ELIJ*DLTYC(I)
     PSTY(IB)=PSTY(IB)-GRAV_E(IB)*DT1(IB)*ELIJ*DLTXC(I)
#    endif
#    endif     

   END DO
#  endif

!#  if !defined (SEMI_IMPLICIT)
#  if defined (SPHERICAL)
#  if !defined (SEMI_IMPLICIT)
   CALL ADVAVE_EDGE_XY(XFLUX,YFLUX,0.0_SP)
#  else
   CALL ADVAVE_EDGE_XY(XFLUX,YFLUX,IFCETA)
#  endif  
#  endif

#  if defined (WET_DRY)
   DO I = 1,N
#    if !defined (SEMI_IMPLICIT)
     ISWETTMP = ISWETCE(I)*ISWETC(I)
#    else
     ISWETTMP = ISWETCT(I)
#    endif
     XFLUX(I) = XFLUX(I)*ISWETTMP
     YFLUX(I) = YFLUX(I)*ISWETTMP
   END DO
#  endif   


#    if defined (PLBC)
     PSTY = 0.0_SP
#    endif


!
!-------------------------SET BOUNDARY VALUES----------------------------------!
!

!  MODIFY BOUNDARY FLUX
   DO I=1,N
     IF(ISBCE(I) == 2) THEN
#      if !defined (SEMI_IMPLICIT)
       XFLUX(I)=XFLUX(I)+Fluxobn(I)*UA(I)*IUCP(I)
       YFLUX(I)=YFLUX(I)+Fluxobn(I)*VA(I)*IUCP(I)
#      else
       XFLUX(I)=0.0_SP
       YFLUX(I)=0.0_SP
#      endif
     ENDIF
   END DO

!  ADJUST FLUX FOR RIVER INFLOW
   IF(NUMQBC > 0) THEN
     IF(RIVER_INFLOW_LOCATION == 'node')THEN
       DO K=1,NUMQBC
         J=INODEQ(K)
         I1=NBVE(J,1)
         I2=NBVE(J,NTVE(J))
         VLCTYQ(K)=QDIS(K)/QAREA(K)
         XFLUX(I1)=XFLUX(I1)-0.5_SP*QDIS(K)*VLCTYQ(K)*COS(ANGLEQ(K))
         YFLUX(I1)=YFLUX(I1)-0.5_SP*QDIS(K)*VLCTYQ(K)*SIN(ANGLEQ(K))
         XFLUX(I2)=XFLUX(I2)-0.5_SP*QDIS(K)*VLCTYQ(K)*COS(ANGLEQ(K))
         YFLUX(I2)=YFLUX(I2)-0.5_SP*QDIS(K)*VLCTYQ(K)*SIN(ANGLEQ(K))
       END DO
     ELSE IF(RIVER_INFLOW_LOCATION == 'edge') THEN
       DO K=1,NUMQBC
         I1=ICELLQ(K)
         VLCTYQ(K)=QDIS(K)/QAREA(K)
         TEMP=QDIS(K)*VLCTYQ(K)
         XFLUX(I1)=XFLUX(I1)-TEMP*COS(ANGLEQ(K))
         YFLUX(I1)=YFLUX(I1)-TEMP*SIN(ANGLEQ(K))
       END DO
     END IF
   END IF

#  if defined (SEMI_IMPLICIT) 
   IF(STAGE==0) THEN

   ELSE

     DO I=1, N

#      if defined (WET_DRY)
       IF(ISWETCT(I) == 1) THEN
#      endif

#        if defined (TWO_D_MODEL)

         XFLUX(I) = ADX2D(I) + XFLUX(I) + DRX2D(I) + PSTX(I) - COR(I)*VA(I)*DT1(I)*ART(I)*EPOR(I)
         YFLUX(I) = ADY2D(I) + YFLUX(I) + DRY2D(I) + PSTY(I) + COR(I)*UA(I)*DT1(I)*ART(I)*EPOR(I)
#        if defined (SPHERICAL)
           XFLUX(I) = XFLUX(I) - UA(I)*VA(I)/REARTH*TAN(DEG2RAD*YC(I))*DT1(I)*ART(I)*EPOR(I)
           YFLUX(I) = YFLUX(I) + UA(I)*UA(I)/REARTH*TAN(DEG2RAD*YC(I))*DT1(I)*ART(I)*EPOR(I)
#        endif

         IF(STAGE<KSTAGE_UV) THEN
           UAF(I) = UAC(I) - RK_UV(STAGE)*DTI*( XFLUX(I)/ART(I)-(-WUSURF(I) + WUBOT(I)) )/DT1(I) - BEDF*UA_N(I)
           VAF(I) = VAC(I) - RK_UV(STAGE)*DTI*( YFLUX(I)/ART(I)-(-WVSURF(I) + WVBOT(I)) )/DT1(I) - BEDF*VA_N(I)

!old:           UAF(I) = UAF(I)-CC_SPONGE(I)*UAF(I)
!old:           VAF(I) = VAF(I)-CC_SPONGE(I)*VAF(I)
! ---- new: Karsten Lettmann: 2012.06.25 -------
            UAF(I) = UAF(I)/(1.0_SP+CC_SPONGE(I)*UAF(I)**2.0_SP)
            VAF(I) = VAF(I)/(1.0_SP+CC_SPONGE(I)*VAF(I)**2.0_SP)
! ------- end new -------------------------------

           IF(DT1(I) > 0.0_SP) THEN
             BTPS= CBC(I)*SQRT(UAF(I)**2+VAF(I)**2)
             WUBOT(I) = -BTPS * UAF(I)
             WVBOT(I) = -BTPS * VAF(I)
           ELSE
             WUBOT(I) = 0.0_SP
             WVBOT(I) = 0.0_SP
           END IF
         ELSE
           IF(ADCOR_ON) THEN
             UBETA2D(I) = XFLUX(I) + COR(I)*VA(I)*DT1(I)*ART(I)*EPOR(I)
             VBETA2D(I) = YFLUX(I) - COR(I)*UA(I)*DT1(I)*ART(I)*EPOR(I)
           ENDIF
         ENDIF

#        else

#        if defined (SPHERICAL)

#        if defined (NORTHPOLE)
         IF(CELL_NORTHAREA(I)==1) THEN

           U_TMP = -VA(I)*COS(XC(I)*DEG2RAD)-UA(I)*SIN(XC(I)*DEG2RAD)
           V_TMP = -VA(I)*SIN(XC(I)*DEG2RAD)+UA(I)*COS(XC(I)*DEG2RAD)

           XFLUX(I) = ADX2D(I) + XFLUX(I) + DRX2D(I) + PSTX(I) - COR(I)*V_TMP*DT1(I)*ART(I)*EPOR(I)
           YFLUX(I) = ADY2D(I) + YFLUX(I) + DRY2D(I) + PSTY(I) + COR(I)*U_TMP*DT1(I)*ART(I)*EPOR(I)

           UAC_TMP = -VAC(I)*COS(XC(I)*DEG2RAD)-UAC(I)*SIN(XC(I)*DEG2RAD)
           VAC_TMP = -VAC(I)*SIN(XC(I)*DEG2RAD)+UAC(I)*COS(XC(I)*DEG2RAD)
           WUSURF_TMP = -WVSURF(I)*COS(XC(I)*DEG2RAD)-WUSURF(I)*SIN(XC(I)*DEG2RAD)
           WVSURF_TMP = -WVSURF(I)*SIN(XC(I)*DEG2RAD)+WUSURF(I)*COS(XC(I)*DEG2RAD)
           WUBOT_TMP = -WVBOT(I)*COS(XC(I)*DEG2RAD)-WUBOT(I)*SIN(XC(I)*DEG2RAD)
           WVBOT_TMP = -WVBOT(I)*SIN(XC(I)*DEG2RAD)+WUBOT(I)*COS(XC(I)*DEG2RAD)

           UAF_TMP = UAC_TMP - RK_UV(STAGE)*DTI*( XFLUX(I)/ART(I)-(-WUSURF_TMP + WUBOT_TMP) )/DT1(I)
           VAF_TMP = VAC_TMP - RK_UV(STAGE)*DTI*( YFLUX(I)/ART(I)-(-WVSURF_TMP + WVBOT_TMP) )/DT1(I)
           UAF(I)  = VAF_TMP*COS(XC(I)*DEG2RAD)-UAF_TMP*SIN(XC(I)*DEG2RAD)
           VAF(I)  = UAF_TMP*COS(XC(I)*DEG2RAD)+VAF_TMP*SIN(XC(I)*DEG2RAD)
           VAF(I)  = -VAF(I)

         ELSE
#        endif
           XFLUX(I) = ADX2D(I) + XFLUX(I) + DRX2D(I) + PSTX(I) - COR(I)*VA(I)*DT1(I)*ART(I)*EPOR(I)
           YFLUX(I) = ADY2D(I) + YFLUX(I) + DRY2D(I) + PSTY(I) + COR(I)*UA(I)*DT1(I)*ART(I)*EPOR(I)
           XFLUX(I) = XFLUX(I) - UA(I)*VA(I)/REARTH*TAN(DEG2RAD*YC(I))*DT1(I)*ART(I)*EPOR(I)
           YFLUX(I) = YFLUX(I) + UA(I)*UA(I)/REARTH*TAN(DEG2RAD*YC(I))*DT1(I)*ART(I)*EPOR(I)

           UAF(I)   = UAC(I) - RK_UV(STAGE)*DTI*( XFLUX(I)/ART(I)-(-WUSURF(I) + WUBOT(I)) )/DT1(I) - BEDF*UA_N(I)
           VAF(I)   = VAC(I) - RK_UV(STAGE)*DTI*( YFLUX(I)/ART(I)-(-WVSURF(I) + WVBOT(I)) )/DT1(I) - BEDF*VA_N(I)

!old:           UAF(I) = UAF(I)-CC_SPONGE(I)*UAF(I)
!old:           VAF(I) = VAF(I)-CC_SPONGE(I)*VAF(I)
! ---- new: Karsten Lettmann: 2012.06.25 -------
            UAF(I) = UAF(I)/(1.0_SP+CC_SPONGE(I)*UAF(I)**2.0_SP)
            VAF(I) = VAF(I)/(1.0_SP+CC_SPONGE(I)*VAF(I)**2.0_SP)
! ------- end new -------------------------------

#        if defined (NORTHPOLE)
         ENDIF
#        endif     

#        else
           XFLUX(I) = ADX2D(I) + XFLUX(I) + DRX2D(I) + PSTX(I) - COR(I)*VA(I)*DT1(I)*ART(I)*EPOR(I)
           YFLUX(I) = ADY2D(I) + YFLUX(I) + DRY2D(I) + PSTY(I) + COR(I)*UA(I)*DT1(I)*ART(I)*EPOR(I)

           UAF(I)   = UAC(I) - RK_UV(STAGE)*DTI*( XFLUX(I)/ART(I)-(-WUSURF(I) + WUBOT(I)) )/DT1(I) - BEDF*UA_N(I)
           VAF(I)   = VAC(I) - RK_UV(STAGE)*DTI*( YFLUX(I)/ART(I)-(-WVSURF(I) + WVBOT(I)) )/DT1(I) - BEDF*VA_N(I)

!old:           UAF(I) = UAF(I)-CC_SPONGE(I)*UAF(I)
!old:           VAF(I) = VAF(I)-CC_SPONGE(I)*VAF(I)
! ---- new: Karsten Lettmann: 2012.06.25 -------
            UAF(I) = UAF(I)/(1.0_SP+CC_SPONGE(I)*UAF(I)**2.0_SP)
            VAF(I) = VAF(I)/(1.0_SP+CC_SPONGE(I)*VAF(I)**2.0_SP)
! ------- end new -------------------------------

#        endif

#        endif
! defined (TWO_D_MODEL)

#      if defined (WET_DRY)
       ELSE
         IF(STAGE<KSTAGE_UV) THEN
           UAF(I) = 0.0_SP
           VAF(I) = 0.0_SP
         ELSE
           XFLUX(I) = 0.0_SP
           YFLUX(I) = 0.0_SP
         ENDIF
       ENDIF
#      endif

     ENDDO

   ENDIF
#  endif

#  if !defined (LIMITED_NO) 
   DEALLOCATE(UIJ1,VIJ1,UIJ2,VIJ2,STAT=ERROR)
   IF(ERROR /= 0) &
   & CALL FATAL_ERROR("Unexpected deallocation error for UIJ1,VIJ1,UIJ2 and VIJ2.")
   DEALLOCATE(UALFA,VALFA,STAT=ERROR)
   IF(ERROR /= 0) &
   & CALL FATAL_ERROR("Unexpected deallocation error for UALFA,VALFA.")
   DEALLOCATE(FXX,FYY,STAT=ERROR)
   IF(ERROR /= 0) &
   & CALL FATAL_ERROR("Unexpected deallocation error for FXX,FYY.")
#  endif
   
   RETURN
   END SUBROUTINE ADVAVE_EDGE_GCY
!==============================================================================|
