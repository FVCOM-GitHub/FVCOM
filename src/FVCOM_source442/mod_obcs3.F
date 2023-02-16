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

MODULE MOD_OBCS3

# if !defined (SEMI_IMPLICIT)
   USE ALL_VARS
   USE MOD_PREC
   USE MOD_OBCS

   USE MOD_MEANFLOW
   USE MOD_OBCS2

   IMPLICIT NONE
   SAVE

!--Nonlinear Velocity Open Boundary Condition Arrays
   REAL(SP), ALLOCATABLE :: FLUXOBN2(:),CXOBC(:),CYOBC(:)
   REAL(SP), ALLOCATABLE :: FLUXOBC2D_X(:),     FLUXOBC2D_Y(:)
   REAL(SP), ALLOCATABLE :: OBC2D_X_TIDE(:),    OBC2D_Y_TIDE(:)
   REAL(SP), ALLOCATABLE :: FLUXOBC3D_X(:,:),   FLUXOBC3D_Y(:,:)
   REAL(SP), ALLOCATABLE :: FLUXOBC3D_X_2(:,:), FLUXOBC3D_Y_2(:,:)
   CONTAINS


!=========================================================================|
   SUBROUTINE ALLOC_OBC3_DATA

   IMPLICIT NONE

   ALLOCATE(FLUXOBN2(0:IOBCN))                 ;FLUXOBN2      = ZERO
   ALLOCATE(CXOBC(0:NT))                       ;CXOBC         = ZERO
   ALLOCATE(CYOBC(0:NT))                       ;CYOBC         = ZERO
   ALLOCATE(FLUXOBC2D_X(0:nmfcell_i))          ;FLUXOBC2D_X   = ZERO
   ALLOCATE(FLUXOBC2D_Y(0:nmfcell_i))          ;FLUXOBC2D_Y   = ZERO
   ALLOCATE(FLUXOBC3D_X(0:nmfcell_i,1:KBM1))   ;FLUXOBC3D_X   = ZERO
   ALLOCATE(FLUXOBC3D_Y(0:nmfcell_i,1:KBM1))   ;FLUXOBC3D_Y   = ZERO
   ALLOCATE(FLUXOBC3D_X_2(0:nmfcell_i,1:KBM1)) ;FLUXOBC3D_X_2 = ZERO
   ALLOCATE(FLUXOBC3D_Y_2(0:nmfcell_i,1:KBM1)) ;FLUXOBC3D_Y_2 = ZERO

   ALLOCATE(OBC2D_X_TIDE(0:nmfcell),OBC2D_Y_TIDE(0:nmfcell)) 

   RETURN
   END SUBROUTINE ALLOC_OBC3_DATA
!==========================================================================|

!==========================================================================|
   SUBROUTINE ZERO_OBC3

   IMPLICIT NONE

   integer :: I

   IF (nmfcell > 0) THEN
      DO I = 0, nmfcell
         OBC2D_X_TIDE(I) = 0.0_SP
         OBC2D_Y_TIDE(I) = 0.0_SP
      END DO
   END IF

   RETURN
   END SUBROUTINE ZERO_OBC3
!==========================================================================|


!==============================================================================|
   SUBROUTINE SETUP_OBC3

#  if defined (SPHERICAL)
   USE MOD_SPHERICAL
#  endif
# if defined (MULTIPROCESSOR)
   USE MOD_PAR  
# endif
   IMPLICIT NONE

   INTEGER  :: I, I1, I2, J, IC, N1, N2, N3
   REAL(SP) :: DXN,DYN,DXC,DYC,CROSS
#  if defined (SPHERICAL)
   REAL(DP) :: x1_dp,y1_dp,x2_dp,y2_dp,side
#  endif

 IF (IOBCN > 0) THEN
   DO I=1,IOBCN
     I1 = I_OBC_N(I)
     I2 = ADJN_OBC(I,1)
     DO J = 1, NTVE(I1)
       IC = NBVE(I1,J)
       N1 = NV(IC,1) ; N2 = NV(IC,2) ; N3 = NV(IC,3)
       IF( N1-I2 == 0 .OR. N2-I2 == 0 .OR. N3-I2 == 0)THEN
# if defined (SPHERICAL)
         x1_dp=VX(I1); y1_dp=VY(I1)
         x2_dp=VX(I2); y2_dp=VY(I2)
         CALL ARCX(x1_dp,y1_dp,x2_dp,y2_dp,side)
         DXN = side; DYN = (VY(I2)-VY(I1))*TPI

         x2_dp=XC(IC); y2_dp=YC(IC)
         CALL ARCX(x1_dp,y1_dp,x2_dp,y2_dp,side)
         DXC = side; DYC = (YC(IC)-VY(I1))*TPI
# else
         DXN = VX(I2)-VX(I1) ; DYN = VY(I2)-VY(I1)
         DXC = XC(IC)-VX(I1) ; DYC = YC(IC)-VY(I1)
# endif
         CROSS     =  SIGN(1.0_SP,DXC*DYN - DYC*DXN)
         CXOBC(IC) =  CROSS*DYN/SQRT(DXN**2 +DYN**2)
         CYOBC(IC) = -CROSS*DXN/SQRT(DXN**2 +DYN**2)
       END IF
     END DO

     IF(NADJN_OBC(I) > 1)THEN
       I2 = ADJN_OBC(I,2)
       DO J = 1, NTVE(I1)
         IC = NBVE(I1,J)
         N1 = NV(IC,1) ; N2 = NV(IC,2) ; N3 = NV(IC,3)
         IF( N1-I2 == 0 .OR. N2-I2 == 0 .OR. N3-I2 == 0)THEN
# if defined (SPHERICAL)
         x1_dp=VX(I1); y1_dp=VY(I1)
         x2_dp=VX(I2); y2_dp=VY(I2)
         CALL ARCX(x1_dp,y1_dp,x2_dp,y2_dp,side)
         DXN = side; DYN = (VY(I2)-VY(I1))*TPI

         x2_dp=XC(IC); y2_dp=YC(IC)
         CALL ARCX(x1_dp,y1_dp,x2_dp,y2_dp,side)
         DXC = side; DYC = (YC(IC)-VY(I1))*TPI
# else
         DXN = VX(I2)-VX(I1) ; DYN = VY(I2)-VY(I1)
         DXC = XC(IC)-VX(I1) ; DYC = YC(IC)-VY(I1)
# endif
           CROSS     = SIGN(1.0_SP,DXC*DYN - DYC*DXN)
           CXOBC(IC) = CROSS*DYN/SQRT(DXN**2 +DYN**2)
           CYOBC(IC) =-CROSS*DXN/SQRT(DXN**2 +DYN**2)
         END IF
       END DO
     END IF
   END DO
 END IF

   RETURN
   END SUBROUTINE SETUP_OBC3

!==============================================================================|
   SUBROUTINE FLUX_OBN2D(KKT)

   IMPLICIT NONE

   INTEGER, INTENT(IN)  :: KKT
   INTEGER  :: I,I1,N1,N2,N3,C1,C2
   REAL(SP) :: E1_X,E1_Y,E2_X,E2_Y
   REAL(SP) :: FLUXF_OBC_1,FLUXF_OBC_2
# if defined (SPHERICAL)
   REAL(DP) X1_DP,Y1_DP,X2_DP,Y2_DP,SIDE
# endif


!   IF (IOBCN > 0) THEN
!   DO I = 1, IOBCN
!     IF(NADJN_OBC(I) == 1)THEN
!       N1 = I_OBC_N(I)
!       N2 = ADJN_OBC(I,1)
!       I1 = I_OBC_NODE(N1)
!
!       E1_Y = VY(N1)-VY(N2)
!#      if defined (SPHERICAL)
!       X1_DP = VX(N2)
!       Y1_DP = VY(N2)
!       X2_DP = VX(N1)
!       Y2_DP = VY(N1)
!       CALL ARCX(X1_DP,Y1_DP,X2_DP,Y2_DP,SIDE)
!       E1_X = SIDE
!       E1_Y = TPI*E1_Y
!#      else
!       E1_X = VX(N1)-VX(N2)
!#      endif
!       
!       C1 = I_OBC_CELL2(ADJC_OBC(I,1)) 
!       FLUXF_OBC_1 = 0.5_SP*SQRT(E1_X**2+E1_Y**2)*  &
!                   (UANP(C1)*CXOBC(ADJC_OBC(I,1))+VANP(C1)*CYOBC(ADJC_OBC(I,1)))
!
!       FLUXOBN2(I) = -FLUXF_OBC_1*(H(N1)+ELT(I1)+ELP(I1)) 
!
!     ELSE
!       N1   = I_OBC_N(I)
!       N2   = ADJN_OBC(I,1)
!       N3   = ADJN_OBC(I,2)
!       I1   = I_OBC_NODE(N1)
!
!       E1_Y = VY(N1)-VY(N2)
!#      if defined (SPHERICAL)
!       X1_DP = VX(N2)
!       Y1_DP = VY(N2)
!       X2_DP = VX(N1)
!       Y2_DP = VY(N1)
!       CALL ARCX(X1_DP,Y1_DP,X2_DP,Y2_DP,SIDE)
!       E1_X = SIDE
!       E1_Y = TPI*E1_Y
!#      else
!       E1_X = VX(N1)-VX(N2)
!#      endif
!
!       E2_Y = VY(N1)-VY(N3)
!#      if defined (SPHERICAL)
!       X1_DP = VX(N3)
!       Y1_DP = VY(N3)
!       X2_DP = VX(N1)
!       Y2_DP = VY(N1)
!       CALL ARCX(X1_DP,Y1_DP,X2_DP,Y2_DP,SIDE)
!       E2_X = SIDE
!       E2_Y = TPI*E2_Y
!#      else
!       E2_X = VX(N1)-VX(N3)
!#      endif
!
!       C1 = I_OBC_CELL2(ADJC_OBC(I,1)) 
!       C2 = I_OBC_CELL2(ADJC_OBC(I,2)) 
!       FLUXF_OBC_1 = 0.5_SP*SQRT(E1_X**2+E1_Y**2)*  &
!                   (UANP(C1)*CXOBC(ADJC_OBC(I,1))+VANP(C1)*CYOBC(ADJC_OBC(I,1)))
!       FLUXF_OBC_2 = 0.5_SP*SQRT(E2_X**2+E2_Y**2)*  &
!                   (UANP(C2)*CXOBC(ADJC_OBC(I,2))+VANP(C2)*CYOBC(ADJC_OBC(I,2)))
!
!       FLUXOBN2(I) = -(FLUXF_OBC_1+FLUXF_OBC_2)*(H(N1)+ELT(I1)+ELP(I1)) 
!     END IF
!   END DO
!   END IF   

   IF (nmfcell > 0) THEN
     DO I = 1, nmfcell
        C1= I_OBC_NODE2(NODE_MFCELL(I,1))
        C2= I_OBC_NODE2(NODE_MFCELL(I,2))
        FLUXOBN2(C1) = FLUXOBN2(C1) - MFQDIS(I)*RDISMF(I,1)
        FLUXOBN2(C2) = FLUXOBN2(C2) - MFQDIS(I)*RDISMF(I,2)
     END DO
   END IF

   RETURN
   END SUBROUTINE FLUX_OBN2D


!==============================================================================|
   SUBROUTINE FLUX_OBC2D

#  if defined (SPHERICAL)
   USE MOD_SPHERICAL
#  endif
   IMPLICIT NONE

   INTEGER  :: I,II,I1,I2,J,J1,J2,ITMP,JTMP
   REAL(DP)  DX12,DY12,TMP1,DTMP
   REAL(SP) :: FLUXF_OBC_1,FLUXF_OBC_2
# if defined (SPHERICAL)
   REAL(DP) X1_DP,Y1_DP,X2_DP,Y2_DP,SIDE
# endif

   IF (nmfcell_i > 0) THEN
     DO I = 1, nmfcell_i
       II= I_MFCELL_N(I)
       ITMP=0
       DO J=1,3
         IF(NBE(II,J) == 0) THEN
           JTMP=J
           ITMP=ITMP+1
         END IF
       END DO
       IF (ITMP /= 1) THEN
         PRINT*,'something is wrong here 2'
         CALL PSTOP
       END IF
    
       J1=JTMP+1-INT((JTMP+1)/4)*3
       J2=JTMP+2-INT((JTMP+2)/4)*3
       I1=NV(II,J1)
       I2=NV(II,J2)
       DY12=VY(I1)-VY(I2)
#      if defined (SPHERICAL)
       X1_DP = VX(I2)
       Y1_DP = VY(I2)
       X2_DP = VX(I1)
       Y2_DP = VY(I1)
       CALL ARCX(X1_DP,Y1_DP,X2_DP,Y2_DP,SIDE)
       DX12 = SIDE
       DY12 = TPI*DY12
#      else
       DX12=VX(I1)-VX(I2)
#      endif

       DTMP = 0.5_SP*(H(I1)+H(I2)+ELT(I_OBC_NODE(I1))+ELT(I_OBC_NODE(I2)))    ! May be a problem, should be replaced by D
!       TMP1 = -(UANT(I)*cos(ANGLEMF(I))+VANT(I)*sin(ANGLEMF(I))*(SQRT(DX12**2+DY12**2))
       TMP1 = UANT(I) * DY12 - VANT(I) * DX12
       FLUXOBC2D_X(I) = DTMP * TMP1 * UANT(I)
       FLUXOBC2D_Y(I) = DTMP * TMP1 * VANT(I)

     END DO
   END IF

   END SUBROUTINE FLUX_OBC2D


!==============================================================================|
   SUBROUTINE FLUX_OBC3D

#  if defined (SPHERICAL)
   USE MOD_SPHERICAL
#  endif
   IMPLICIT NONE

   INTEGER  :: I,II,I1,I2,J,J1,J2,K,ITMP,JTMP
   REAL(DP)  DX12,DY12,TMP1,DTMP
   REAL(SP) :: UTMP,VTMP
# if defined (SPHERICAL)
   REAL(DP) X1_DP,Y1_DP,X2_DP,Y2_DP,SIDE
# endif

   IF (nmfcell > 0) THEN

! 2-D and 3-D adjustment
      OBC2D_X_TIDE = OBC2D_X_TIDE/FLOAT(ISPLIT)
      OBC2D_Y_TIDE = OBC2D_Y_TIDE/FLOAT(ISPLIT)
                                                                                                                         
      DO I = 1, nmfcell
         II= I_MFCELL_N(I)
         ITMP=0
         DO J=1,3
           IF(NBE(II,J) == 0 .and. ISONB(nv(II,J)) /= 2) THEN
             JTMP=J
             ITMP=ITMP+1
           END IF
         END DO
         IF (ITMP /= 1) THEN
           PRINT*,'something is wrong here 3'
           CALL PSTOP
         END IF

         J1=JTMP+1-INT((JTMP+1)/4)*3
         J2=JTMP+2-INT((JTMP+2)/4)*3
         I1=NV(II,J1)
         I2=NV(II,J2)
         DTMP =	0.5_SP*(H(I1)+H(I2)+ELTDT(I_OBC_NODE(I1))+ELTDT(I_OBC_NODE(I2)))      ! May be a problem, should be replaced by D

         UTMP = 0.0_SP ; VTMP = 0.0_SP
         DO K=1,KBM1
            UTMP = UTMP + UNT(I,K)*DZ1(II,K)
            VTMP = VTMP + VNT(I,K)*DZ1(II,K)
         END DO
         UTMP = UTMP * DTMP
         VTMP = VTMP * DTMP
         DO K=1,KBM1
           UNT(I,K) = UNT(I,K) - (UTMP-OBC2D_X_TIDE(I))/DTMP
           VNT(I,K) = VNT(I,K) - (VTMP-OBC2D_Y_TIDE(I))/DTMP
         END DO
      END DO
   END IF

   IF (nmfcell_i > 0) THEN
     DO I = 1, nmfcell_i
       II= I_MFCELL_N(I)
       ITMP=0
       DO J=1,3
         IF(NBE(II,J) == 0) THEN
           JTMP=J
           ITMP=ITMP+1
         END IF
       END DO
       J1=JTMP+1-INT((JTMP+1)/4)*3
       J2=JTMP+2-INT((JTMP+2)/4)*3
       I1=NV(II,J1)
       I2=NV(II,J2)
       DY12=VY(I1)-VY(I2)
#      if defined (SPHERICAL)
       X1_DP = VX(I2)
       Y1_DP = VY(I2)
       X2_DP = VX(I1)
       Y2_DP = VY(I1)
       CALL ARCX(X1_DP,Y1_DP,X2_DP,Y2_DP,SIDE)
       DX12 = SIDE
       DY12 = TPI*DY12
#      else
       DX12=VX(I1)-VX(I2)
#      endif

       DTMP = 0.5_SP*(H(I1)+H(I2)+ELTDT(I_OBC_NODE(I1))+ELTDT(I_OBC_NODE(I2)))       ! May be a problem, should be replaced by D
       DO K =1, KBM1
!          TMP1 = -(UNT(I,K)*cos(ANGLEMF(I))+VNT(I,K)*sin(ANGLEMF(I))*(SQRT(DX12**2+DY12**2))
          TMP1 = UNT(I,K) * DY12 - VNT(I,K) * DX12
          FLUXOBC3D_X(I,K) = DTMP * TMP1 * UNT(I,K)
          FLUXOBC3D_Y(I,K) = DTMP * TMP1 * VNT(I,K)
       END DO
     END DO

   END IF

   END SUBROUTINE FLUX_OBC3D


!==============================================================================|
   SUBROUTINE FLUX_OBC3D_2

   USE ALL_VARS
   USE MOD_OBCS
   USE MOD_OBCS2

#  if defined (SPHERICAL)
   USE MOD_SPHERICAL
#  endif
   IMPLICIT NONE

   INTEGER  :: I,II,I1,I2,J,J1,J2,K,ITMP,JTMP
   REAL(DP)  DX12,DY12,TMP1,DTMP
   REAL(SP) :: UTMP,VTMP
# if defined (SPHERICAL)
   REAL(DP) X1_DP,Y1_DP,X2_DP,Y2_DP,SIDE
# endif

   IF (nmfcell > 0) THEN

! 2-D and 3-D adjustment
      DO I = 1, nmfcell
       II= I_MFCELL_N(I)
         UTMP   = SUM(UNT(I,1:KBM1)*DZ1(II,1:KBM1))
         VTMP   = SUM(VNT(I,1:KBM1)*DZ1(II,1:KBM1))
         UNT(I,1:KBM1) = UNT(I,1:KBM1) + (UANT(I) - UTMP) 
         VNT(I,1:KBM1) = VNT(I,1:KBM1) + (VANT(I) - VTMP) 
      END DO
   END IF

   IF (nmfcell_i > 0) THEN
     DO I = 1, nmfcell_i
       II= I_MFCELL_N(I)
       ITMP=0
       DO J=1,3
         IF(NBE(II,J) == 0) THEN
           JTMP=J
           ITMP=ITMP+1
         END IF
       END DO
       J1=JTMP+1-INT((JTMP+1)/4)*3
       J2=JTMP+2-INT((JTMP+2)/4)*3
       I1=NV(II,J1)
       I2=NV(II,J2)
       DY12=VY(I1)-VY(I2)
#      if defined (SPHERICAL)
       X1_DP = VX(I2)
       Y1_DP = VY(I2)
       X2_DP = VX(I1)
       Y2_DP = VY(I1)
       CALL ARCX(X1_DP,Y1_DP,X2_DP,Y2_DP,SIDE)
       DX12 = SIDE
       DY12 = TPI*DY12
#      else
       DX12=VX(I1)-VX(I2)
#      endif

       DTMP = 0.5_SP*(H(I1)+H(I2)+ELTDT(I_OBC_NODE(I1))+ELTDT(I_OBC_NODE(I2)))          ! May be a problem, should be replaced by D
       DO K =1, KBM1
!          TMP1 = -(UNT(I,K)*cos(ANGLEMF(I))+VNT(I,K)*sin(ANGLEMF(I))*(SQRT(DX12**2+DY12**2))
          TMP1 = UNT(I,K) * DY12 - VNT(I,K) * DX12
          FLUXOBC3D_X_2(I,K) = DTMP * TMP1 * UNT(I,K)
          FLUXOBC3D_Y_2(I,K) = DTMP * TMP1 * VNT(I,K)
       END DO
     END DO

   END IF

   END SUBROUTINE FLUX_OBC3D_2

# endif
!========================================================================
END MODULE MOD_OBCS3
