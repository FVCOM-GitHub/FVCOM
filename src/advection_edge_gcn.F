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
   SUBROUTINE ADVECTION_EDGE_GCN(XFLUX,YFLUX)
!==============================================================================|
!   Calculate the Advection and Diffusion Terms of 3D Velocity Field           |
!   These Terms will be vertically integrated to form the Mean Terms in        |
!   the Gx and Gy Terms of the External Mode Equation                          |
!==============================================================================|

   USE ALL_VARS
   USE MOD_UTILS
   USE BCS
   USE MOD_SPHERICAL
   USE MOD_NORTHPOLE
   USE MOD_WD

#  if defined (MEAN_FLOW)
   USE MOD_MEANFLOW
   USE MOD_OBCS2
   USE MOD_OBCS3
#  endif

#  if defined (THIN_DAM)
   USE MOD_DAM   !Jadon
#  endif

   IMPLICIT NONE
   REAL(SP), INTENT(OUT), DIMENSION(0:NT,KB) :: XFLUX,YFLUX
   REAL(SP) :: DIJ
   REAL(SP) :: COFA1,COFA2,COFA3,COFA4,COFA5,COFA6,COFA7,COFA8
   REAL(SP) :: XADV,YADV,TXXIJ,TYYIJ,TXYIJ,UN_TMP
   REAL(SP) :: VISCOF,VISCOF1,VISCOF2,TEMP,TPA,TPB
   REAL(SP) :: XIJA,YIJA,XIJB,YIJB,UIJ,VIJ
   REAL(SP) :: FACT,FM1
   INTEGER  :: I,IA,IB,J1,J2,K1,K2,K3,K4,K5,K6,K,II,J,I1,I2
   REAL(SP) :: ISWETTMP

!#  if defined (THIN_DAM)
   REAL(SP) :: A1UIA1,A1UIA2,A1UIA3,A1UIA4,A2UIA1,A2UIA2,A2UIA3,A2UIA4
   REAL(SP) :: A1UIB1,A1UIB2,A1UIB3,A1UIB4,A2UIB1,A2UIB2,A2UIB3,A2UIB4
   INTEGER  :: J11,J12,J21,J22,E1,E2,ISBCE1,ISBC_TMP,IB_TMP
   LOGICAL  :: ISMATCH
!#  endif

#  if defined (LIMITED_NO)
   REAL(SP) :: UIJ1,VIJ1,UIJ2,VIJ2,FXX,FYY
#  else
   REAL(SP),ALLOCATABLE,DIMENSION(:) :: UIJ1,VIJ1,UIJ2,VIJ2,FXX,FYY
   REAL(SP),ALLOCATABLE,DIMENSION(:) :: UALFA,VALFA
   REAL(SP) :: UALFA_TMP,VALFA_TMP
   INTEGER :: ERROR
   REAL(SP) :: EPS
#  endif
!------------------------------------------------------------------------------|
   
   if(dbg_set(dbg_sbr)) write(ipt,*) "Start: advection_edge_gcn.F"

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
!--Initialize Variables--------------------------------------------------------|
!
   XFLUX = 0.0_SP
   YFLUX = 0.0_SP

!
!--Loop Over Edges and Accumulate Fluxes-For Each Element----------------------|
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
     UIJ1=0.0_SP;  VIJ1=0.0_SP;  UIJ2=0.0_SP;  VIJ2=0.0_SP
     UALFA=1.0_SP; VALFA=1.0_SP
     FXX=0.0_SP;   FYY=0.0_SP

     DO I=1,NE
       IA=IEC(I,1)
       IB=IEC(I,2)

#      if defined (WET_DRY)
#      if !defined (SEMI_IMPLICIT)
       IF(ISWETCT(IA)*ISWETC(IA) == 1 .OR. ISWETCT(IB)*ISWETC(IB) == 1)THEN
#      else
       IF(ISWETCT(IA) == 1 .OR. ISWETCT(IB) == 1)THEN
#      endif
#      endif
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

       DIJ= 0.5_SP*(DT(J1)*DZ(J1,K)+DT(J2)*DZ(J2,K))
       IB_TMP = IB

!--------------Used for Dam Model by Jadon----------------------
# if defined (THIN_DAM)
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
# else
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
# endif
!---------------------------------------------------------------
       !!FORM THE LEFT FLUX
       COFA1=A1UIA1*U(IA,K)+A1UIA2*U(K1,K)+A1UIA3*U(K2,K)+A1UIA4*U(K3,K)
       COFA2=A2UIA1*U(IA,K)+A2UIA2*U(K1,K)+A2UIA3*U(K2,K)+A2UIA4*U(K3,K)
       COFA5=A1UIA1*V(IA,K)+A1UIA2*V(K1,K)+A1UIA3*V(K2,K)+A1UIA4*V(K3,K)
       COFA6=A2UIA1*V(IA,K)+A2UIA2*V(K1,K)+A2UIA3*V(K2,K)+A2UIA4*V(K3,K)
!       COFA1=A1U(IA,1)*U(IA,K)+A1U(IA,2)*U(K1,K)+A1U(IA,3)*U(K2,K)+A1U(IA,4)*U(K3,K)
!       COFA2=A2U(IA,1)*U(IA,K)+A2U(IA,2)*U(K1,K)+A2U(IA,3)*U(K2,K)+A2U(IA,4)*U(K3,K)
!       COFA5=A1U(IA,1)*V(IA,K)+A1U(IA,2)*V(K1,K)+A1U(IA,3)*V(K2,K)+A1U(IA,4)*V(K3,K)
!       COFA6=A2U(IA,1)*V(IA,K)+A2U(IA,2)*V(K1,K)+A2U(IA,3)*V(K2,K)+A2U(IA,4)*V(K3,K)

       UIJ1(I)=COFA1*XIJA+COFA2*YIJA
       VIJ1(I)=COFA5*XIJA+COFA6*YIJA
       UALFA_TMP=ABS(U(IA,K)-U(IB_TMP,K))/ABS(UIJ1(I)+EPSILON(EPS))
       VALFA_TMP=ABS(V(IA,K)-V(IB_TMP,K))/ABS(VIJ1(I)+EPSILON(EPS))
       IF(UALFA_TMP > 1)UALFA_TMP = 1.0_SP
       IF(VALFA_TMP > 1)VALFA_TMP = 1.0_SP
       UALFA(IA)=MIN(UALFA(IA),UALFA_TMP)
       VALFA(IA)=MIN(VALFA(IA),VALFA_TMP)
     
       !!FORM THE RIGHT FLUX
       COFA3=A1UIB1*U(IB_TMP,K)+A1UIB2*U(K4,K)+A1UIB3*U(K5,K)+A1UIB4*U(K6,K)
       COFA4=A2UIB1*U(IB_TMP,K)+A2UIB2*U(K4,K)+A2UIB3*U(K5,K)+A2UIB4*U(K6,K)
       COFA7=A1UIB1*V(IB_TMP,K)+A1UIB2*V(K4,K)+A1UIB3*V(K5,K)+A1UIB4*V(K6,K)
       COFA8=A2UIB1*V(IB_TMP,K)+A2UIB2*V(K4,K)+A2UIB3*V(K5,K)+A2UIB4*V(K6,K)
!       COFA3=A1U(IB,1)*U(IB,K)+A1U(IB,2)*U(K4,K)+A1U(IB,3)*U(K5,K)+A1U(IB,4)*U(K6,K)
!       COFA4=A2U(IB,1)*U(IB,K)+A2U(IB,2)*U(K4,K)+A2U(IB,3)*U(K5,K)+A2U(IB,4)*U(K6,K)
!       COFA7=A1U(IB,1)*V(IB,K)+A1U(IB,2)*V(K4,K)+A1U(IB,3)*V(K5,K)+A1U(IB,4)*V(K6,K)
!       COFA8=A2U(IB,1)*V(IB,K)+A2U(IB,2)*V(K4,K)+A2U(IB,3)*V(K5,K)+A2U(IB,4)*V(K6,K)

       UIJ2(I)=COFA3*XIJB+COFA4*YIJB
       VIJ2(I)=COFA7*XIJB+COFA8*YIJB
       UALFA_TMP=ABS(U(IA,K)-U(IB_TMP,K))/ABS(UIJ2(I)+EPSILON(EPS))
       VALFA_TMP=ABS(V(IA,K)-V(IB_TMP,K))/ABS(VIJ2(I)+EPSILON(EPS))
       IF(UALFA_TMP > 1)UALFA_TMP = 1.0_SP
       IF(VALFA_TMP > 1)VALFA_TMP = 1.0_SP
       UALFA(IB_TMP)=MIN(UALFA(IB_TMP),UALFA_TMP)
       VALFA(IB_TMP)=MIN(VALFA(IB_TMP),VALFA_TMP)
     
       VISCOF1=ART(IA)*SQRT(COFA1**2+COFA6**2+0.5_SP*(COFA2+COFA5)**2)
       VISCOF2=ART(IB_TMP)*SQRT(COFA3**2+COFA8**2+0.5_SP*(COFA4+COFA7)**2)
        
!       VISCOF = HORCON*(FACT*0.5_SP*(VISCOF1+VISCOF2)/HPRNU + FM1)
!       VISCOF = HORCON*(FACT*0.5_SP*(VISCOF1+VISCOF2) + FM1)/HPRNU
       ! David moved HPRNU and added HVC
       VISCOF=(FACT*0.5_SP*(VISCOF1*CC_HVC(IA)+VISCOF2*CC_HVC(IB_TMP)) + FM1*0.5_SP*(CC_HVC(IA)+CC_HVC(IB_TMP)))/HPRNU

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

#      if defined (WET_DRY)
#      if !defined (SEMI_IMPLICIT)
       IF(ISWETCT(IA)*ISWETC(IA) == 1 .OR. ISWETCT(IB)*ISWETC(IB) == 1)THEN
#      else
       IF(ISWETCT(IA) == 1 .OR. ISWETCT(IB) == 1)THEN
#      endif
#      endif
       J1=IENODE(I,1)
       J2=IENODE(I,2)

       IB_TMP = IB

#      if defined (THIN_DAM)
       IF(IB==0.AND.E_DAM_MATCH(IA)/=0.AND.K<=KDAM1(IA))IB_TMP=E_DAM_MATCH(IA) 
#      endif

       DIJ= 0.5_SP*(DT(J1)*DZ(J1,K)+DT(J2)*DZ(J2,K))
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

       !!COMPUTE THE NORMAL VELOCITY ACROSS THE EDGE
       UIJ=0.5_SP*(UIJ1(I)+UIJ2(I))
       VIJ=0.5_SP*(VIJ1(I)+VIJ2(I))
       UN_TMP=VIJ*DLTXC(I) - UIJ*DLTYC(I)

       !!UPWIND THE ADVECTIVE FLUX
       XADV=DIJ*UN_TMP*((1.0_SP-SIGN(1.0_SP,UN_TMP))*UIJ2(I)+(1.0_SP+SIGN(1.0_SP,UN_TMP))*UIJ1(I))*0.5_SP
       YADV=DIJ*UN_TMP*((1.0_SP-SIGN(1.0_SP,UN_TMP))*VIJ2(I)+(1.0_SP+SIGN(1.0_SP,UN_TMP))*VIJ1(I))*0.5_SP
        

       !!COMPUTE BOUNDARY FLUX AUGMENTERS   
#      if !defined (MEAN_FLOW)

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
       !!ACCUMULATE THE FLUX
!       XFLUX(IA,K)=XFLUX(IA,K)+XADV*TPA+FXX*TPA
!       YFLUX(IA,K)=YFLUX(IA,K)+YADV*TPA+FYY*TPA
!       XFLUX(IB,K)=XFLUX(IB,K)-XADV*TPB-FXX*TPB
!       YFLUX(IB,K)=YFLUX(IB,K)-YADV*TPB-FYY*TPB

        XFLUX(IA,K)=XFLUX(IA,K)+XADV*TPA+(FXX(I)+3.0_SP*FXX(I)*FLOAT(ISBC_TMP))*EPOR(IA)
        YFLUX(IA,K)=YFLUX(IA,K)+YADV*TPA+(FYY(I)+3.0_SP*FYY(I)*FLOAT(ISBC_TMP))*EPOR(IA)
        XFLUX(IB,K)=XFLUX(IB,K)-XADV*TPB-(FXX(I)+3.0_SP*FXX(I)*FLOAT(ISBC_TMP))*EPOR(IB)
        YFLUX(IB,K)=YFLUX(IB,K)-YADV*TPB-(FYY(I)+3.0_SP*FYY(I)*FLOAT(ISBC_TMP))*EPOR(IB)
#       else
        TPA = FLOAT(1-ISBC(I))
        TPB = FLOAT(1-ISBC(I))
        XFLUX(IA,K)=XFLUX(IA,K)+(XADV*TPA+(FXX(I)+3.0_SP*FXX(I)*FLOAT(ISBC(I))))*IUCP(IA)
        YFLUX(IA,K)=YFLUX(IA,K)+(YADV*TPA+(FYY(I)+3.0_SP*FYY(I)*FLOAT(ISBC(I))))*IUCP(IA)
        XFLUX(IB,K)=XFLUX(IB,K)-(XADV*TPB+(FXX(I)+3.0_SP*FXX(I)*FLOAT(ISBC(I))))*IUCP(IB)
        YFLUX(IB,K)=YFLUX(IB,K)-(YADV*TPB+(FYY(I)+3.0_SP*FYY(I)*FLOAT(ISBC(I))))*IUCP(IB)
#       endif

#       if defined (WET_DRY)
        END IF
#       endif
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

#  if defined (WET_DRY)
#   if !defined (SEMI_IMPLICIT)
    IF(ISWETCT(IA)*ISWETC(IA) == 1 .OR. ISWETCT(IB)*ISWETC(IB) == 1)THEN
#   else
    IF(ISWETCT(IA) == 1 .OR. ISWETCT(IB) == 1)THEN
#   endif
#  endif
     J1=IENODE(I,1)
     J2=IENODE(I,2)

     K1=NBE(IA,1)
     K2=NBE(IA,2)
     K3=NBE(IA,3)
     K4=NBE(IB,1)
     K5=NBE(IB,2)
     K6=NBE(IB,3)
!#    if defined (SPHERICAL)
!     XIJA=DLTXNE(I,1)
!     YIJA=DLTYNE(I,1)
!     XIJB=DLTXNE(I,2)
!     YIJB=DLTYNE(I,2)
!#    if defined (THIN_DAM)
!     IF(IB==0.AND.E_DAM_MATCH(IA)/=0)XIJB=DLTXNE_DAM_MATCH(I)
!     IF(IB==0.AND.E_DAM_MATCH(IA)/=0)YIJB=DLTYNE_DAM_MATCH(I)
!#    endif
!#    else
!     XIJA=XIJC(I)-XC(IA)
!     YIJA=YIJC(I)-YC(IA)
!     XIJB=XIJC(I)-XC(IB)
!     YIJB=YIJC(I)-YC(IB)
!#    if defined (THIN_DAM)
!     IF(IB==0.AND.E_DAM_MATCH(IA)/=0)XIJB=XIJC(I)-XC(E_DAM_MATCH(IA))
!     IF(IB==0.AND.E_DAM_MATCH(IA)/=0)YIJB=YIJC(I)-YC(E_DAM_MATCH(IA))
!#    endif
!#    endif

     DO K=1,KBM1

#    if defined (SPHERICAL)
     XIJA=DLTXNE(I,1)
     YIJA=DLTYNE(I,1)
     XIJB=DLTXNE(I,2)
     YIJB=DLTYNE(I,2)
#    if defined (THIN_DAM)
     IF(IB==0.AND.E_DAM_MATCH(IA)/=0.AND.K<=KDAM1(IA))XIJB=DLTXNE_DAM_MATCH(I)
     IF(IB==0.AND.E_DAM_MATCH(IA)/=0.AND.K<=KDAM1(IA))YIJB=DLTYNE_DAM_MATCH(I)
#    endif
#    else
     XIJA=XIJC(I)-XC(IA)
     YIJA=YIJC(I)-YC(IA)
     XIJB=XIJC(I)-XC(IB)
     YIJB=YIJC(I)-YC(IB)
#    if defined (THIN_DAM)
     IF(IB==0.AND.E_DAM_MATCH(IA)/=0.AND.K<=KDAM1(IA))XIJB=XIJC(I)-XC(E_DAM_MATCH(IA))
     IF(IB==0.AND.E_DAM_MATCH(IA)/=0.AND.K<=KDAM1(IA))YIJB=YIJC(I)-YC(E_DAM_MATCH(IA))
#    endif
#    endif

       IB_TMP = IB
       DIJ= 0.5_SP*(DT(J1)*DZ(J1,K)+DT(J2)*DZ(J2,K))
!--------------Used for Dam Model by Jadon----------------------
# if defined (THIN_DAM)
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
       
       IF(IB==0.AND.E_DAM_MATCH(IA)/=0)IB_TMP=E_DAM_MATCH(IA) 
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
# else
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
# endif
!---------------------------------------------------------------
       !!FORM THE LEFT FLUX
       COFA1=A1UIA1*U(IA,K)+A1UIA2*U(K1,K)+A1UIA3*U(K2,K)+A1UIA4*U(K3,K)
       COFA2=A2UIA1*U(IA,K)+A2UIA2*U(K1,K)+A2UIA3*U(K2,K)+A2UIA4*U(K3,K)
       COFA5=A1UIA1*V(IA,K)+A1UIA2*V(K1,K)+A1UIA3*V(K2,K)+A1UIA4*V(K3,K)
       COFA6=A2UIA1*V(IA,K)+A2UIA2*V(K1,K)+A2UIA3*V(K2,K)+A2UIA4*V(K3,K)
!       COFA1=A1U(IA,1)*U(IA,K)+A1U(IA,2)*U(K1,K)+A1U(IA,3)*U(K2,K)+A1U(IA,4)*U(K3,K)
!       COFA2=A2U(IA,1)*U(IA,K)+A2U(IA,2)*U(K1,K)+A2U(IA,3)*U(K2,K)+A2U(IA,4)*U(K3,K)
!       COFA5=A1U(IA,1)*V(IA,K)+A1U(IA,2)*V(K1,K)+A1U(IA,3)*V(K2,K)+A1U(IA,4)*V(K3,K)
!       COFA6=A2U(IA,1)*V(IA,K)+A2U(IA,2)*V(K1,K)+A2U(IA,3)*V(K2,K)+A2U(IA,4)*V(K3,K)
       UIJ1=U(IA,K)+COFA1*XIJA+COFA2*YIJA
       VIJ1=V(IA,K)+COFA5*XIJA+COFA6*YIJA

       !!FORM THE RIGHT FLUX
       COFA3=A1UIB1*U(IB_TMP,K)+A1UIB2*U(K4,K)+A1UIB3*U(K5,K)+A1UIB4*U(K6,K)
       COFA4=A2UIB1*U(IB_TMP,K)+A2UIB2*U(K4,K)+A2UIB3*U(K5,K)+A2UIB4*U(K6,K)
       COFA7=A1UIB1*V(IB_TMP,K)+A1UIB2*V(K4,K)+A1UIB3*V(K5,K)+A1UIB4*V(K6,K)
       COFA8=A2UIB1*V(IB_TMP,K)+A2UIB2*V(K4,K)+A2UIB3*V(K5,K)+A2UIB4*V(K6,K)
!       COFA3=A1U(IB,1)*U(IB,K)+A1U(IB,2)*U(K4,K)+A1U(IB,3)*U(K5,K)+A1U(IB,4)*U(K6,K)
!       COFA4=A2U(IB,1)*U(IB,K)+A2U(IB,2)*U(K4,K)+A2U(IB,3)*U(K5,K)+A2U(IB,4)*U(K6,K)
!       COFA7=A1U(IB,1)*V(IB,K)+A1U(IB,2)*V(K4,K)+A1U(IB,3)*V(K5,K)+A1U(IB,4)*V(K6,K)
!       COFA8=A2U(IB,1)*V(IB,K)+A2U(IB,2)*V(K4,K)+A2U(IB,3)*V(K5,K)+A2U(IB,4)*V(K6,K)
       UIJ2=U(IB_TMP,K)+COFA3*XIJB+COFA4*YIJB
       VIJ2=V(IB_TMP,K)+COFA7*XIJB+COFA8*YIJB

       !!COMPUTE THE NORMAL VELOCITY ACROSS THE EDGE
       UIJ=0.5_SP*(UIJ1+UIJ2)
       VIJ=0.5_SP*(VIJ1+VIJ2)
       UN_TMP=VIJ*DLTXC(I) - UIJ*DLTYC(I)

       VISCOF1=ART(IA)*SQRT(COFA1**2+COFA6**2+0.5_SP*(COFA2+COFA5)**2)
       VISCOF2=ART(IB_TMP)*SQRT(COFA3**2+COFA8**2+0.5_SP*(COFA4+COFA7)**2)
        
!       VISCOF = HORCON*(FACT*0.5_SP*(VISCOF1+VISCOF2)/HPRNU + FM1)
!       VISCOF = HORCON*(FACT*0.5_SP*(VISCOF1+VISCOF2) + FM1)
       ! David moved HPRNU and added HVC
       VISCOF=(FACT*0.5_SP*(VISCOF1*CC_HVC(IA)+VISCOF2*CC_HVC(IB_TMP)) + FM1*0.5_SP*(CC_HVC(IA)+CC_HVC(IB_TMP)))/HPRNU

       TXXIJ=(COFA1+COFA3)*VISCOF
       TYYIJ=(COFA6+COFA8)*VISCOF
       TXYIJ=0.5_SP*(COFA2+COFA4+COFA5+COFA7)*VISCOF
       FXX=DIJ*(TXXIJ*DLTYC(I)-TXYIJ*DLTXC(I))
       FYY=DIJ*(TXYIJ*DLTYC(I)-TYYIJ*DLTXC(I))

       !!UPWIND THE ADVECTIVE FLUX
       XADV=DIJ*UN_TMP*((1.0_SP-SIGN(1.0_SP,UN_TMP))*UIJ2+(1.0_SP+SIGN(1.0_SP,UN_TMP))*UIJ1)*0.5_SP
       YADV=DIJ*UN_TMP*((1.0_SP-SIGN(1.0_SP,UN_TMP))*VIJ2+(1.0_SP+SIGN(1.0_SP,UN_TMP))*VIJ1)*0.5_SP
        

       !!COMPUTE BOUNDARY FLUX AUGMENTERS   
#  if !defined (MEAN_FLOW)

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
       !!ACCUMULATE THE FLUX
!       XFLUX(IA,K)=XFLUX(IA,K)+XADV*TPA+FXX*TPA
!       YFLUX(IA,K)=YFLUX(IA,K)+YADV*TPA+FYY*TPA
!       XFLUX(IB,K)=XFLUX(IB,K)-XADV*TPB-FXX*TPB
!       YFLUX(IB,K)=YFLUX(IB,K)-YADV*TPB-FYY*TPB

        XFLUX(IA,K)=XFLUX(IA,K)+XADV*TPA+(FXX+3.0_SP*FXX*FLOAT(ISBC_TMP))*EPOR(IA)
        YFLUX(IA,K)=YFLUX(IA,K)+YADV*TPA+(FYY+3.0_SP*FYY*FLOAT(ISBC_TMP))*EPOR(IA)
        XFLUX(IB,K)=XFLUX(IB,K)-XADV*TPB-(FXX+3.0_SP*FXX*FLOAT(ISBC_TMP))*EPOR(IB)
        YFLUX(IB,K)=YFLUX(IB,K)-YADV*TPB-(FYY+3.0_SP*FYY*FLOAT(ISBC_TMP))*EPOR(IB)

!        XFLUX(IA,K)=XFLUX(IA,K)+XADV*TPA+(FXX+3.0_SP*FXX*FLOAT(ISBC(I)))*EPOR(IA)
!        YFLUX(IA,K)=YFLUX(IA,K)+YADV*TPA+(FYY+3.0_SP*FYY*FLOAT(ISBC(I)))*EPOR(IA)
!        XFLUX(IB,K)=XFLUX(IB,K)-XADV*TPB-(FXX+3.0_SP*FXX*FLOAT(ISBC(I)))*EPOR(IB)
!        YFLUX(IB,K)=YFLUX(IB,K)-YADV*TPB-(FYY+3.0_SP*FYY*FLOAT(ISBC(I)))*EPOR(IB)
#  else
        TPA = FLOAT(1-ISBC(I))
        TPB = FLOAT(1-ISBC(I))
        XFLUX(IA,K)=XFLUX(IA,K)+(XADV*TPA+(FXX+3.0_SP*FXX*FLOAT(ISBC(I))))*IUCP(IA)
        YFLUX(IA,K)=YFLUX(IA,K)+(YADV*TPA+(FYY+3.0_SP*FYY*FLOAT(ISBC(I))))*IUCP(IA)
        XFLUX(IB,K)=XFLUX(IB,K)-(XADV*TPB+(FXX+3.0_SP*FXX*FLOAT(ISBC(I))))*IUCP(IB)
        YFLUX(IB,K)=YFLUX(IB,K)-(YADV*TPB+(FYY+3.0_SP*FYY*FLOAT(ISBC(I))))*IUCP(IB)
#  endif

     END DO

#  if defined (WET_DRY)
    END IF
#  endif
   END DO

#  endif

#  if defined (SPHERICAL)
   CALL ADVECTION_EDGE_XY(XFLUX,YFLUX)
#  endif  

#  if defined (WET_DRY)
   DO I=1,N
#    if !defined (SEMI_IMPLICIT)
     ISWETTMP = ISWETCT(I)*ISWETC(I)
#    else
     ISWETTMP = ISWETCT(I)
#    endif
     DO K=1,KBM1
       XFLUX(I,K) = XFLUX(I,K)*ISWETTMP
       YFLUX(I,K) = YFLUX(I,K)*ISWETTMP
     END DO
   END DO
#  endif       	 


!
!--Boundary Conditions on Flux-------------------------------------------------|
!
#  if !defined (MEAN_FLOW)
      DO I=1,N
        IF(ISBCE(I) == 2)THEN
           DO K=1,KBM1
             XFLUX(I,K)=0.0_SP
             YFLUX(I,K)=0.0_SP
           END DO
        END IF
      END DO
#  else
   IF (nmfcell_i > 0) THEN
     DO II=1,nmfcell_i
        I1=I_MFCELL_N(II)
        DO K=1,KBM1
           XFLUX(I1,K) = XFLUX(I1,K) + FLUXOBC3D_X_2(II,K)*IUCP(I1)
           YFLUX(I1,K) = YFLUX(I1,K) + FLUXOBC3D_Y_2(II,K)*IUCP(I1)
        END DO
     END DO
   END IF
#  endif

!
!--Adjust Flux for Fresh Water Inflow------------------------------------------|
!

   IF(NUMQBC > 0) THEN
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
     END IF
   END IF

!
!--Adjust Flux for Open Boundary Inflow/Outflow------------------------------------------|
!
#  if defined (MEAN_FLOW)
   IF(nmfcell_i > 0) THEN
     DO II=1,nmfcell_i
       I1=I_MFCELL_N(II)
       DO K=1,KBM1
         VLCTYMF(II)=MFQDIS(II)/MFAREA(II)
!         TEMP=MFQDIS(II)*MFDIST(II,K)*VLCTYMF(II)
         TEMP=MFQDIS(II)*MFDIST(II,K)*MFDIST(II,K)*VLCTYMF(II)/DZ1(I1,K)
!         XFLUX(I1,K)=XFLUX(I1,K)-TEMP/DZ1(I1,K)*COS(ANGLEMF(II))
!         YFLUX(I1,K)=YFLUX(I1,K)-TEMP/DZ1(I1,K)*SIN(ANGLEMF(II))
         XFLUX(I1,K)=XFLUX(I1,K)-TEMP*COS(ANGLEMF(II))
         YFLUX(I1,K)=YFLUX(I1,K)-TEMP*SIN(ANGLEMF(II))
       END DO
     END DO
   END IF
#  endif

   if(dbg_set(dbg_sbr)) write(ipt,*) "End: advection_edge_gcn.F"

   RETURN
   END SUBROUTINE ADVECTION_EDGE_GCN
!==============================================================================|
