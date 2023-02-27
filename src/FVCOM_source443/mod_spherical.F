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
!  MODULE CONTAINING SUBROUTINES USED TO SET UP SPHERICAL COORDINATE SYSTEM    |
!  FLUXES                                                                      |
!==============================================================================|

MODULE MOD_SPHERICAL
  USE CONTROL
  USE MOD_PREC
  IMPLICIT NONE
  SAVE
  REAL(SP), ALLOCATABLE :: DLTXNE(:,:)
  REAL(SP), ALLOCATABLE :: DLTYNE(:,:)

  REAL(SP), ALLOCATABLE :: DELTUX(:,:)
  REAL(SP), ALLOCATABLE :: DELTUY(:,:)
  REAL(SP), ALLOCATABLE :: SITAU(:,:)


  INTERFACE ARC
     MODULE PROCEDURE ARC_FLT
     MODULE PROCEDURE ARC_DBL
  END INTERFACE

  INTERFACE ARCC
     MODULE PROCEDURE ARCC_FLT
     MODULE PROCEDURE ARCC_DBL
  END INTERFACE

  INTERFACE AREA
     MODULE PROCEDURE AREA_FLT
     MODULE PROCEDURE AREA_DBL
  END INTERFACE

  INTERFACE ARCX
     MODULE PROCEDURE ARCX_FLT
     MODULE PROCEDURE ARCX_DBL
  END INTERFACE

  INTERFACE ARCY
     MODULE PROCEDURE ARCY_FLT
     MODULE PROCEDURE ARCY_DBL
  END INTERFACE

  INTERFACE ARCX_BACK
     MODULE PROCEDURE ARCX_BACK_FLT
     MODULE PROCEDURE ARCX_BACK_DBL
  END INTERFACE


  !===================================================================================|
CONTAINS   !!INCLUDED SUBROUTINES FOLLOW
  !===================================================================================|

  SUBROUTINE ARC_DBL(XX1,YY1,XX2,YY2,ARCL)
    !----------------------------------------------------------------------------
    !      function:
    !           calculate the arc lenth for given two point on the spherical plane
    !      input:
    !           xx1,yy1,xx2,yy2 :are longitude and latitude of two points
    !      output:
    !           arcl :  arc lenth of two points in spherical plane
    !-----------------------------------------------------------------------------       

    !  solve the arc length through the earth center
    IMPLICIT NONE
    REAL(DP) :: X1,Y1,X2,Y2,XA,YA,ZA,XB,YB,ZB,AB,AOB
    REAL(DP),INTENT(OUT) :: ARCL
    REAL(DP),INTENT(IN) :: XX1,YY1,XX2,YY2

    X1=XX1*DEG2RAD
    Y1=YY1*DEG2RAD

    X2=XX2*DEG2RAD
    Y2=YY2*DEG2RAD

    ! USE DOUBLE PRECISION COS AND SIN
    XA=DCOS(Y1)*DCOS(X1)
    YA=DCOS(Y1)*DSIN(X1)
    ZA=DSIN(Y1)

    XB=DCOS(Y2)*DCOS(X2)
    YB=DCOS(Y2)*DSIN(X2)
    ZB=DSIN(Y2)

    AB=DSQRT((XB-XA)**2+(YB-YA)**2+(ZB-ZA)**2)
    AOB=(2.0_DP -AB*AB)/2.0_DP
    AOB=DACOS(AOB)
    ARCL=REARTH*AOB

    RETURN
  END SUBROUTINE ARC_DBL

  SUBROUTINE ARC_FLT(XX1,YY1,XX2,YY2,ARCL)
    IMPLICIT NONE
    REAL(SPA), INTENT(IN) :: XX1,YY1,XX2,YY2
    REAL(SPA), INTENT(OUT) :: ARCL
    REAL(DP) ARCL_DP
    CALL ARC_DBL(DBLE(xx1),DBLE(YY1),DBLE(XX2),DBLE(YY2),ARCL_DP)
    ARCL = ARCL_DP
  END SUBROUTINE ARC_FLT


  SUBROUTINE AREA_DBL(SIDEA,SIDEB,SIDEC,AREA1)
    !--------------------------------------------------------------------
    !      function:
    !           calculate the area of a triangle on a spherical plane
    !      input:
    !           side1,side2 and side3: are 3 arc lenth for one triangle
    !      output:
    !           areal: is area of a triangle on a spherical plane
    !--------------------------------------------------------------------
    IMPLICIT NONE
    REAL(DP), INTENT(IN) :: SIDEA,SIDEB,SIDEC
    REAL(DP), INTENT(OUT) :: AREA1
    REAL(DP) :: SIDE1,SIDE2,SIDE3
    REAL(DP) :: PSUM,PM,QMJC

    SIDE1=SIDEA/REARTH
    SIDE2=SIDEB/REARTH
    SIDE3=SIDEC/REARTH

    ! SLOWER TO CHECK THEN TO CALCULATE
    !   IF(SIDE1 == 0.0_DP .OR. SIDE2 == 0.0_DP .OR. SIDE3 == 0.0_DP)THEN
    !     AREA1=0.0_DP
    !   ELSE

    PSUM=0.5_DP*(SIDE1+SIDE2+SIDE3)
    PM=DSIN(PSUM)*DSIN(PSUM-SIDE1)*DSIN(PSUM-SIDE2)*DSIN(PSUM-SIDE3)
    PM=DSQRT(PM)/(2.0_DP*DCOS(SIDE1*0.5_DP)*DCOS(SIDE2*0.5_DP)*DCOS(SIDE3*0.5_DP))
    QMJC = 2.0_DP*DASIN(PM)

    AREA1=REARTH*REARTH*QMJC

    !   END IF

    RETURN
  END SUBROUTINE AREA_DBL

  SUBROUTINE AREA_FLT(SIDE1,SIDE2,SIDE3,AREA1)
    IMPLICIT NONE
    REAL(SPA), INTENT(IN) :: SIDE1,SIDE2,SIDE3
    REAL(SPA), INTENT(OUT) :: AREA1
    REAL(DP) :: AREA_DP

    CALL AREA_DBL(DBLE(SIDE1),DBLE(SIDE2),DBLE(SIDE3),AREA_DP)
    AREA1=AREA_DP

  END SUBROUTINE AREA_FLT


  SUBROUTINE ARCC_DBL(XX1,YY1,XX2,YY2,XXC,YYC)
    IMPLICIT NONE
    REAL(DP), INTENT(OUT) :: XXC,YYC
    REAL(DP), INTENT(IN) :: XX1,YY1,XX2,YY2
    REAL(DP) :: X1,Y1,X2,Y2

    X1=XX1*DEG2RAD
    Y1=YY1*DEG2RAD

    X2=XX2*DEG2RAD
    Y2=YY2*DEG2RAD

    XXC=DCOS(Y1)*DSIN(X1)+DCOS(Y2)*DSIN(X2)
    !   XXC=XXC/(COS(Y1)*COS(X1)+COS(Y2)*COS(X2))
    !   XXC=ATAN(XXC)
    XXC=DATAN2(XXC,(DCOS(Y1)*DCOS(X1)+DCOS(Y2)*DCOS(X2)))
    XXC=XXC/DEG2RAD

    !   IF(XXC .LT. 0.0) XXC=180.0+XXC
    IF(XXC < 0.0_DP) XXC=360.0_DP+XXC

    YYC=DCOS(Y1)*DCOS(Y1)+DCOS(Y2)*DCOS(Y2)+2.0_DP*DCOS(Y1)*DCOS(Y2)*DCOS(X1-X2)
    !   YYC=SQRT(YYC)/(SIN(Y1)+SIN(Y2))
    YYC=DATAN2(DSQRT(YYC),(DSIN(Y1)+DSIN(Y2)))
    !   YYC=ATAN(YYC)
    YYC=90.0_DP-YYC/DEG2RAD

    RETURN
  END SUBROUTINE ARCC_DBL

  SUBROUTINE ARCC_FLT(XX1,YY1,XX2,YY2,XXC,YYC)
    IMPLICIT NONE
    REAL(SPA) :: XX1,YY1,XX2,YY2
    REAL(SPA) :: XXC,YYC
    REAL(DP) :: XXC_DP,YYC_DP

    CALL ARCC_DBL(DBLE(XX1),DBLE(YY1),DBLE(XX2),DBLE(YY2),XXC_DP,YYC_DP)
    XXC = XXC_DP
    YYC = YYC_DP

  END SUBROUTINE ARCC_FLT


  SUBROUTINE ARCX_DBL(XX1,YY1,XX2,YY2,ARCX1)
    IMPLICIT NONE
    REAL(DP), INTENT(IN) :: XX1,YY1,XX2,YY2
    REAL(DP), INTENT(OUT)::ARCX1

    REAL(DP) :: X1,Y1,X2,Y2,TY
    REAL(DP) :: XTMP	      

    IF(XX1 == XX2)THEN
       ARCX1=0.0_DP
    ELSE
       X1=XX1*DEG2RAD
       Y1=YY1*DEG2RAD

       X2=XX2*DEG2RAD
       Y2=YY2*DEG2RAD

       XTMP  = X2-X1
       IF(XTMP >  PI)THEN
          XTMP = REAL(-2*PI,DP)+XTMP
       ELSE IF(XTMP < -PI)THEN
          XTMP =  REAL(2*PI,DP)+XTMP
       END IF

       TY=0.5_DP*(Y2+Y1)
       ARCX1=REARTH*DCOS(TY)*XTMP
    END IF

    RETURN
  END SUBROUTINE ARCX_DBL

  SUBROUTINE ARCY_DBL(XX1,YY1,XX2,YY2,ARCY1)
    IMPLICIT NONE
    REAL(DP), INTENT(IN) :: XX1,YY1,XX2,YY2
    REAL(DP), INTENT(OUT)::ARCY1

    REAL(DP) :: X1,Y1,X2,Y2,TY
    REAL(DP) :: YTMP	      

    IF(YY1 == YY2)THEN
       ARCY1=0.0_DP
    ELSE
       X1=XX1*DEG2RAD
       Y1=YY1*DEG2RAD

       X2=XX2*DEG2RAD
       Y2=YY2*DEG2RAD

       YTMP  = Y2-Y1
       IF(YTMP >  PI)THEN
          YTMP = REAL(-2*PI,DP)+YTMP
       ELSE IF(YTMP < -PI)THEN
          YTMP =  REAL(2*PI,DP)+YTMP
       END IF

       ARCY1=REARTH*YTMP
    END IF

    RETURN
  END SUBROUTINE ARCY_DBL

  SUBROUTINE ARCY_FLT(XX1,YY1,XX2,YY2,ARCY1)
    IMPLICIT NONE
    REAL(SPA), INTENT(IN) :: XX1,YY1,XX2,YY2
    REAL(SPA), INTENT(OUT)::ARCY1
    
    REAL(DP) ::ARCY_DP

    CALL ARCY_DBL(DBLE(XX1),DBLE(YY1),DBLE(XX2),DBLE(YY2),ARCY_DP)
    ARCY1 = ARCY_DP

  END SUBROUTINE ARCY_FLT

  SUBROUTINE ARCX_FLT(XX1,YY1,XX2,YY2,ARCX1)
    IMPLICIT NONE
    REAL(SPA), INTENT(IN) :: XX1,YY1,XX2,YY2
    REAL(SPA), INTENT(OUT)::ARCX1
    
    REAL(DP) ::ARCX_DP

    CALL ARCX_DBL(DBLE(XX1),DBLE(YY1),DBLE(XX2),DBLE(YY2),ARCX_DP)
    ARCX1 = ARCX_DP

  END SUBROUTINE ARCX_FLT

  SUBROUTINE ARCX_BACK_DBL(XX1,YY1,XX2,YY2,ARCX1)
    IMPLICIT NONE
    REAL(DP), INTENT(IN) :: XX1,YY1,XX2,YY2
    REAL(DP), INTENT(OUT) :: ARCX1

    INTEGER I
    INTEGER,PARAMETER ::NX=500
    REAL(DP) :: X1,Y1,X2,Y2,TY,A1,A2,B1,B2,C1,C2,A,B,C,X(NX+1),Y(NX+1)
    REAL(DP) :: XTMP	      

    IF(XX1 == XX2)THEN
       ARCX1=0.
    ELSE
       X1=XX1*DEG2RAD
       Y1=YY1*DEG2RAD

       X2=XX2*DEG2RAD
       Y2=YY2*DEG2RAD

       X(1)=X1
       Y(1)=Y1
       X(NX+1)=X2
       Y(NX+1)=Y2

       XTMP=X(NX+1)-X(1)
       IF(XTMP >  PI)THEN
          XTMP = REAL(-2*PI,DP)+XTMP
       ELSE IF(XTMP < -PI)THEN
          XTMP =  REAL(2*PI,DP)+XTMP
       END IF

       DO I=2,NX
          X(I)=X(I-1)+XTMP/FLOAT(NX)
          !       x(i)=x(i-1)+(x(nx+1)-x(1))/float(nx)
       END DO

       A1=DCOS(Y(1))*DCOS(X(1))
       A2=DCOS(Y(NX+1))*DCOS(X(NX+1))

       B1=DCOS(Y(1))*DSIN(X(1))
       B2=DCOS(Y(NX+1))*DSIN(X(NX+1))

       C1=DSIN(Y(1))
       C2=DSIN(Y(NX+1))

       A=A1*B2-A2*B1
       B=B1*C2-B2*C1
       C=A2*C1-A1*C2

       DO I=2,NX
          Y(I)=-B*DCOS(X(I))-C*DSIN(X(I))
          Y(I)=Y(I)/A
          Y(I)=DATAN(Y(I))
       END DO

       ARCX1=0.
       DO I=1,NX
          TY=0.5*(Y(I)+Y(I+1))
          XTMP=X(I+1)-X(I)
          IF(XTMP >  PI)THEN
             XTMP = real(-2*PI,DP)+XTMP
          ELSE IF(XTMP < -PI)THEN
             XTMP =  real(2*PI,DP)+XTMP
          END IF
          ARCX1=ARCX1+REARTH*DCOS(TY)*XTMP
          !       arcx1=arcx1+rearth*cos(ty)*(x(i+1)-x(i))
       END DO
    END IF

    RETURN
  END SUBROUTINE ARCX_BACK_DBL

    SUBROUTINE ARCX_BACK_FLT(XX1,YY1,XX2,YY2,ARCX1)
    IMPLICIT NONE
    REAL(SPA), INTENT(IN) :: XX1,YY1,XX2,YY2
    REAL(SPA), INTENT(OUT)::ARCX1
    
    REAL(DP) ::ARCX_DP

    CALL ARCX_BACK_DBL(DBLE(XX1),DBLE(YY1),DBLE(XX2),DBLE(YY2),ARCX_DP)
    ARCX1 = ARCX_DP

  END SUBROUTINE ARCX_BACK_FLT

  SUBROUTINE ALLOC_SPHERE_VARS
    USE LIMS
    INTEGER NCT
    INTEGER NDB
# if !defined (DOUBLE_PRECISION)
    NDB = 1       !!GWC BASE THIS ON KIND
# else
    NDB = 2
# endif

    NCT = NT*3
    ALLOCATE(DLTXNE(NCT,2))       ;DLTXNE    = ZERO
    ALLOCATE(DLTYNE(NCT,2))       ;DLTYNE    = ZERO
    ALLOCATE(DELTUX(NT,3))        ;DELTUX    = ZERO
    ALLOCATE(DELTUY(NT,3))        ;DELTUY    = ZERO
    ALLOCATE(SITAU(NT,3))         ;SITAU     = ZERO

    memcnt = memcnt + NCT*4*NDB + NT*9*NDB
    RETURN
  END SUBROUTINE ALLOC_SPHERE_VARS


END MODULE  MOD_SPHERICAL
