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

!===========================================================================!  
! n-dimentional optimal interpolation, Fortran module                       !
! Modified from the original program optimal_interpolation.F which was      !
! written by Alexander Barth <abarth@marine.usf.edu>                        !
! Dependencies: LAPACK (dsyev)                                              !
!===========================================================================!

   MODULE MOD_OPTIMAL_INTERPOLATION
#  if defined (DATA_ASSIM)
   USE MOD_PREC
   USE CONTROL
   
   CONTAINS

!--select the m observations closest to x

   SUBROUTINE SELECT_NEAREST(X,OX,PARAM,MM,INDEX,DISTANCE)
   IMPLICIT NONE

   REAL(SP),INTENT(IN) :: X(:),OX(:,:),PARAM(:,:)
   INTEGER, INTENT(IN) :: MM
   REAL(SP),INTENT(OUT) :: DISTANCE(MM)
   INTEGER, INTENT(OUT) :: INDEX(MM)
   
   REAL(SP) :: D(SIZE(OX,2))
   INTEGER  :: I

   DO I=1,SIZE(OX,2)
     D(I) = SUM(((X - OX(:,I)) * PARAM(:,I))**2)
   END DO

   CALL SORT(D,MM,INDEX)

   DISTANCE = D(INDEX)

   RETURN
   END SUBROUTINE SELECT_NEAREST


!--returns the indexes of the m smallest elements in d

   SUBROUTINE SORT(D,MM,PANNIER)
   IMPLICIT NONE
   REAL(SP),INTENT(IN)  :: D(:)
   INTEGER, INTENT(IN)  :: MM
   INTEGER, INTENT(OUT) :: PANNIER(MM)

   INTEGER :: I,MAX_PANNIER(MM)
   
   DO I=1,MM
     PANNIER(I) = I
   END DO

   MAX_PANNIER = MAXLOC(D(PANNIER))

   DO I=MM+1,SIZE(D)
     IF(D(I) < D(PANNIER(MAX_PANNIER(1))))THEN
       PANNIER(MAX_PANNIER(1)) = I
       MAX_PANNIER = MAXLOC(D(PANNIER))
     END IF
   END DO

   RETURN
   END SUBROUTINE SORT

   SUBROUTINE OBSERVATION_COVARIANCE(OVAR,INDEX,R)
   IMPLICIT NONE
   REAL(SP),INTENT(IN)  :: OVAR(:)
   INTEGER, INTENT(IN)  :: INDEX(:)
   REAL(SP),INTENT(OUT) :: R(SIZE(INDEX),SIZE(INDEX))

   INTEGER :: I
   
   R = 0.0_SP

   DO I=1,SIZE(INDEX)
     R(I,I) = OVAR(INDEX(I))
   END DO

   RETURN
   END SUBROUTINE OBSERVATION_COVARIANCE


   FUNCTION BACKGROUND_COVARIANCE(X1,X2,PARAM) RESULT(C)
   IMPLICIT NONE
   REAL(SP),INTENT(IN) :: X1(:),X2(:),PARAM(:)
   REAL(SP) :: C

   REAL(SP) :: DIS(SIZE(X1))

   DIS = (X1 - X2)*PARAM

   C = EXP(-SUM(DIS**2))

   END FUNCTION BACKGROUND_COVARIANCE 

!--compute pseudo-inverse

   SUBROUTINE PINV(A,TOLERANCE,WORK,C)
   IMPLICIT NONE
   REAL(SP),INTENT(IN) :: A(:,:),TOLERANCE,WORK(:)
   REAL(SP) :: C(SIZE(A,1),SIZE(A,1))

   INTEGER :: I,J,K,INFO,NP
   REAL(SP) :: WP(SIZE(A,1)), UP(SIZE(A,1),SIZE(A,1))

   NP = SIZE(A,1)
   UP = A

#  if defined (DOUBLE_PRECISION)   
   CALL DSYEV('V','U', NP, UP,NP, WP, WORK, SIZE(WORK), INFO)
#  else   
   CALL SSYEV('V','U', NP, UP,NP, WP, WORK, SIZE(WORK), INFO)
#  endif   

   DO I=1,NP
     IF(WP(I) > TOLERANCE)THEN
       WP(I) = 1.0_SP/WP(I)
     ELSE
       WP(I) = 0.0_SP
     END IF
   END DO

   DO K=1,NP
     DO J=1,NP
       C(J,K) = 0.0_SP
       DO I=1,NP
         C(J,K) = C(J,K) + UP(J,I) * WP(I) * UP(K,I)
       END DO
     END DO
   END DO

   RETURN
   END SUBROUTINE PINV

!--query the necessary workspace

   FUNCTION PINV_WORKSPACE(N) RESULT(LWORK)
   IMPLICIT NONE
   INTEGER,INTENT(IN) :: N
   INTEGER :: LWORK

   INTEGER :: INFO
   REAL(SP):: DUMMY,RWORK

#  if defined (DOUBLE_PRECISION)   
   CALL DSYEV('V','U', N, DUMMY, N, DUMMY, RWORK, -1, INFO)
#  else   
   CALL SSYEV('V','U', N, DUMMY, N, DUMMY, RWORK, -1, INFO)
#  endif   
   LWORK = CEILING(RWORK)

   END FUNCTION


!--main optimal interpolation routine

   SUBROUTINE OPTIMINTERP(OX,OF,OVAR,PARAM,MM,GX,GF,GVAR)
   IMPLICIT NONE
   REAL(SP),INTENT(IN)  :: GX(:,:),OX(:,:),OF(:),OVAR(:),PARAM(:,:)
   INTEGER, INTENT(IN)  :: MM
   REAL(SP),INTENT(OUT) :: GF(:),GVAR(:)

   REAL(SP) :: HPH(MM,MM), R(MM,MM), PH(MM), IA(MM,MM), PHIA(MM)
   INTEGER  :: GN,INDEX(MM)
   REAL(SP) :: DISTANCE(MM)

   INTEGER  :: I,J1,J2,LWORK
   REAL(SP) :: TOLERANCE = 1E-5

   REAL(SP),ALLOCATABLE :: WORK(:)

   GN = SIZE(GX,2)
   
!--query and allocate workspace for pseudo-inverse
   LWORK = PINV_WORKSPACE(MM)

!$omp parallel private(work,i,iA,PHiA,index,distance,HPH,j1,j2)

   ALLOCATE(WORK(LWORK))
      
!$omp do 
   DO I=1,GN

!--get the indexes of the nearest observations

     CALL SELECT_NEAREST(GX(:,I),OX,PARAM,MM,INDEX,DISTANCE)
     
!--form compute the error covariance matrix of the observation 

     CALL OBSERVATION_COVARIANCE(OVAR,INDEX,R)

!--form the error covariance matrix background field

     DO J1=1,MM
       DO J2=1,MM
         HPH(J1,J2) =     &
	 BACKGROUND_COVARIANCE(OX(:,INDEX(J1)),OX(:,INDEX(J2)),PARAM(:,INDEX(J2)))
       END DO

       PH(J1) = BACKGROUND_COVARIANCE(GX(:,I),OX(:,INDEX(J1)),PARAM(:,INDEX(J1)))
     END DO

!--covariance matrix of the innovation

     IA = HPH + R

!--pseudo inverse of the covariance matrix of the innovation

     CALL PINV(IA,TOLERANCE,WORK,IA)

     PHIA = MATMUL(PH,IA)

!--compute the analysis

     GF(I) = DOT_PRODUCT(PHIA,OF(INDEX))

!--compute the error variance of the analysis

     GVAR(I) = 1. - DOT_PRODUCT(PHIA,PH)

   END DO
!$omp end do

   DEALLOCATE(WORK)
 
!$omp end parallel 

   RETURN
   END SUBROUTINE OPTIMINTERP
#  endif
   END MODULE MOD_OPTIMAL_INTERPOLATION
