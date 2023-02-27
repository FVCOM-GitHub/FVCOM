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

!--------------------------------------------------------------------|
! Nonhydrostatic Module                                              |
!--------------------------------------------------------------------|
MODULE NON_HYDRO
# if defined (NH)

  USE MOD_PREC, ONLY : SP, DP
# if defined (WET_DRY)
  USE MOD_WD, ONLY: ISWETC, ISWETCT, ISWETN, ISWETNT, WET_JUDGE, WD_UPDATE
# endif

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: ALLOC_VARS_NH
  PUBLIC :: ADV_W
  PUBLIC :: VDIF_W
  PUBLIC :: GEN_POISSON
  PUBLIC :: UPDATE_QUVW
  PUBLIC :: NAME_LIST_INITIALIZE_NH
  PUBLIC :: NAME_LIST_PRINT_NH
  PUBLIC :: NAME_LIST_READ_NH

! must include all necessary headers in petsc fortran interface
#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscda.h"
#include "include/finclude/petscmat.h"
#include "include/finclude/petscksp.h"
#include "include/finclude/petscpc.h"
#include "include/finclude/petscis.h"
#include "include/finclude/petscis.h90"
#include "include/finclude/petscao.h"
#include "include/finclude/petscvec.h90"
#include "include/finclude/petscviewer.h"

  INTEGER, PUBLIC                ::  PROJ_SWITCH
  REAL(SP), PUBLIC, ALLOCATABLE  ::  W4ZT(:,:), RHS(:,:), QP(:,:)
  REAL(SP), PUBLIC, ALLOCATABLE  ::  NHQDRX(:,:), NHQDRY(:,:), NHQDRZ(:,:), NHQ2DX(:), NHQ2DY(:)

  REAL(SP), ALLOCATABLE  ::  AVG_U(:,:), AVG_V(:,:)
  REAL(SP), ALLOCATABLE  ::  NHA1(:,:) , NHA2(:,:)
  REAL(SP), ALLOCATABLE  ::  NHA1D(:)  , NHA2D(:)
  REAL(SP), ALLOCATABLE  ::  A1AVG(:,:), A2AVG(:,:)
  REAL(SP), ALLOCATABLE  ::  AVELX(:)  , AVELY(:)
  REAL(SP), ALLOCATABLE  ::  A1AVGD(:) , A2AVGD(:)
  REAL(SP), ALLOCATABLE  ::  W4Z(:,:)
  REAL(SP), ALLOCATABLE  ::  NHQ(:,:)  , QP1(:,:), NHQ1(:,:)

  NAMELIST /NML_NH/          &
          & PROJ_SWITCH 

  PetscViewer :: viewer

  !----------------------------------------------------------------------------
  CONTAINS
  ! ALLOC_VARS_NH            : 
  ! A1A2                     :
  ! ADV_W, VDIF_W            : w momentum calculation
  ! POISSON_SOLVE            : solve poisson problem for NH pressure q
  !   --> PETSc_SetICS       : set initial condition for petsc solution vec x
  !   --> PETSc_Solve        : solve Ax=b for x
  ! UPDATE_QUVW              : correct divergent free velocity
  ! UPDATE_ZETA              : correct free surface 
  !   --> PETSc_SetICS_EL    :
  !   --> PETSc_Solve_EL     :  
  ! PETSc_GetQ               :
  ! TRIDAG                   :
  !----------------------------------------------------------------------------

  SUBROUTINE NAME_LIST_INITIALIZE_NH
  USE CONTROL, ONLY : IPT

  IMPLICIT NONE

  !--Parameters in NameList NML_NH
  PROJ_SWITCH    = -1 

  END SUBROUTINE NAME_LIST_INITIALIZE_NH

  SUBROUTINE NAME_LIST_PRINT_NH
  USE CONTROL, ONLY : IPT

  IMPLICIT NONE

  WRITE (UNIT=IPT,NML=NML_NH)

  RETURN
  END SUBROUTINE NAME_LIST_PRINT_NH

  SUBROUTINE NAME_LIST_READ_NH

  USE CONTROL
  USE MOD_UTILS
  IMPLICIT NONE
  INTEGER  :: ISCAN
  CHARACTER(LEN=120) :: FNAME
  INTEGER :: ios

  ios = 0

  FNAME = "./"//trim(casename)//"_run.nml"

  if(DBG_SET(dbg_io)) &
         & write(IPT,*) "Read_Name_List: File: ",trim(FNAME)

  CALL FOPEN(NMLUNIT,trim(FNAME),'cfr')

  ! Read NH Settings
  READ(UNIT=NMLUNIT, NML=NML_NH,IOSTAT=ios)
  if(ios .NE. 0 ) then
    if(DBG_SET(dbg_log)) write(UNIT=IPT,NML=NML_NH)
    Call Fatal_Error("Can Not Read NameList NML_NH from file: "//trim(FNAME))
  endif

  REWIND(NMLUNIT)
  CLOSE(NMLUNIT)

  RETURN
  END SUBROUTINE NAME_LIST_READ_NH

  SUBROUTINE ALLOC_VARS_NH

  USE MOD_PREC, ONLY : SP, DP
  USE ALL_VARS, ONLY : GRAV, M, N, MT, NT, MGL, KBM1, KB, NTSN, NBSN, ISONB
  IMPLICIT NONE

  ALLOCATE(NHA1(0:NT,KBM1))           ; NHA1     = 0.0_SP
  ALLOCATE(NHA2(0:NT,KBM1))           ; NHA2     = 0.0_SP
  ALLOCATE(NHA1D(0:NT))               ; NHA1D    = 0.0_SP
  ALLOCATE(NHA2D(0:NT))               ; NHA2D    = 0.0_SP
  ALLOCATE(A1AVG(0:MT,KBM1))          ; A1AVG    = 0.0_SP
  ALLOCATE(A2AVG(0:MT,KBM1))          ; A2AVG    = 0.0_SP
  ALLOCATE(AVELX(0:MT))               ; AVELX    = 0.0_SP
  ALLOCATE(AVELY(0:MT))               ; AVELY    = 0.0_SP
  ALLOCATE(A1AVGD(0:MT))              ; A1AVGD   = 0.0_SP
  ALLOCATE(A2AVGD(0:MT))              ; A2AVGD   = 0.0_SP

  ALLOCATE(AVG_U(0:MT,KBM1))          ; AVG_U    = 0.0_SP
  ALLOCATE(AVG_V(0:MT,KBM1))          ; AVG_V    = 0.0_SP
  ALLOCATE(W4ZT(0:MT,KB))             ; W4ZT     = 0.0_SP
  ALLOCATE(W4Z(0:MT,KB))              ; W4Z      = 0.0_SP

  ALLOCATE(RHS(0:MT,KBM1))            ; RHS      = 0.0_SP

  ALLOCATE(NHQDRX(0:NT,KBM1))         ; NHQDRX   = 0.0_SP
  ALLOCATE(NHQDRY(0:NT,KBM1))         ; NHQDRY   = 0.0_SP
  ALLOCATE(NHQDRZ(0:MT,KBM1))         ; NHQDRZ   = 0.0_SP
  ALLOCATE(NHQ2DX(0:NT))              ; NHQ2DX   = 0.0_SP
  ALLOCATE(NHQ2DY(0:NT))              ; NHQ2DY   = 0.0_SP
  ALLOCATE(QP(0:MT,KBM1))             ; QP       = 0.0_SP
  ALLOCATE(QP1(0:NT,KBM1))            ; QP1      = 0.0_SP
  ALLOCATE(NHQ(0:MT,KBM1))            ; NHQ      = 0.0_SP
  ALLOCATE(NHQ1(0:NT,KBM1))           ; NHQ1     = 0.0_SP

  RETURN
  END SUBROUTINE ALLOC_VARS_NH

  SUBROUTINE A1A2(EL_TMP, D_TMP)

  USE LIMS
  USE CONTROL
  USE ALL_VARS
  IMPLICIT NONE

  INTEGER N1, N2, N3, I, J, K, I1, I2
  REAL(SP) ::  EL_TMP(0:MT), D_TMP(0:MT), F1, F2, F3
  REAL(SP), ALLOCATABLE  ::  DELX(:), DELY(:), DDX(:), DDY(:)

  ALLOCATE(DELX(0:NT))                ; DELX     = 0.0_SP
  ALLOCATE(DELY(0:NT))                ; DELY     = 0.0_SP
  ALLOCATE(DDX(0:NT))                 ; DDX      = 0.0_SP
  ALLOCATE(DDY(0:NT))                 ; DDY      = 0.0_SP

! calculate A1, A2
  DO K=1, KBM1
    DO I=1, NT
      N1 = NV(I,1)
      N2 = NV(I,2)
      N3 = NV(I,3)
    
      IF(K==1) THEN
        DELX(I) = AWX(I,1)*EL_TMP(N1)+AWX(I,2)*EL_TMP(N2)+AWX(I,3)*EL_TMP(N3)
        DELY(I) = AWY(I,1)*EL_TMP(N1)+AWY(I,2)*EL_TMP(N2)+AWY(I,3)*EL_TMP(N3)
      ENDIF
      DDX(I)  = AWX(I,1)*D_TMP(N1)*ZZ(N1,K)+AWX(I,2)*D_TMP(N2)*ZZ(N2,K)+AWX(I,3)*D_TMP(N3)*ZZ(N3,K)
      DDY(I)  = AWY(I,1)*D_TMP(N1)*ZZ(N1,K)+AWY(I,2)*D_TMP(N2)*ZZ(N2,K)+AWY(I,3)*D_TMP(N3)*ZZ(N3,K)
      NHA1(I,K) = -(DELX(I)+DDX(I))
      NHA2(I,K) = -(DELY(I)+DDY(I))
      IF(K==KBM1) THEN
        NHA1D(I)  = AWX(I,1)*H(N1)+AWX(I,2)*H(N2)+AWX(I,3)*H(N3)
        NHA2D(I)  = AWY(I,1)*H(N1)+AWY(I,2)*H(N2)+AWY(I,3)*H(N3)
      ENDIF
    ENDDO
  ENDDO

  AVELX  = 0.0_SP
  AVELY  = 0.0_SP
  A1AVG  = 0.0_SP
  A2AVG  = 0.0_SP
  A1AVGD = 0.0_SP
  A2AVGD = 0.0_SP
  DO K=1, KBM1
    DO I=1, M
      DO J=1, NTVE(I)
        A1AVG(I,K) = A1AVG(I,K) + NHA1(NBVE(I,J),K)*ART(NBVE(I,J))
        A2AVG(I,K) = A2AVG(I,K) + NHA2(NBVE(I,J),K)*ART(NBVE(I,J))
        IF(K==KBM1) THEN
          AVELX(I)  = AVELX(I)  + DELX(NBVE(I,J))*ART(NBVE(I,J))
          AVELY(I)  = AVELY(I)  + DELY(NBVE(I,J))*ART(NBVE(I,J))    

          A1AVGD(I) = A1AVGD(I) + NHA1D(NBVE(I,J))*ART(NBVE(I,J))
          A2AVGD(I) = A2AVGD(I) + NHA2D(NBVE(I,J))*ART(NBVE(I,J))
        ENDIF
      ENDDO
      A1AVG(I,K) = A1AVG(I,K)/ART2(I)
      A2AVG(I,K) = A2AVG(I,K)/ART2(I)
      IF(K==KBM1) THEN
        AVELX(I)  = AVELX(I)/ART2(I)
        AVELY(I)  = AVELY(I)/ART2(I)

        A1AVGD(I) = A1AVGD(I)/ART2(I)
        A2AVGD(I) = A2AVGD(I)/ART2(I)
      ENDIF
    ENDDO
  ENDDO

  DEALLOCATE(DELX,DELY,DDX,DDY)
  RETURN
  END SUBROUTINE A1A2

  SUBROUTINE ADV_W                
                            
  USE LIMS
  USE CONTROL
  USE ALL_VARS
  USE MOD_UTILS
# if defined(MULTIPROCESSOR)
  USE MOD_PAR, ONLY: NBN, BN_MLT, BN_LOC, BNC, NODE_MATCH 
# endif
  IMPLICIT NONE

  INTEGER  ::  I, J, K, I1, I2, IA, IB, IERR
  REAL(SP) ::  F1, XI, YI, DXA, DYA, DXB, DYB, FIJ1, FIJ2
  REAL(SP) ::  TXX, TYY, FXX, FYY, VISCOF, UN, EXFLUX
  REAL(SP) ::  W1MIN, W1MAX, W2MIN, W2MAX
  REAL(SP) ::  FACT, FM1
  REAL(SP), ALLOCATABLE  :: U_TMP(:,:), V_TMP(:,:)
  REAL(SP), ALLOCATABLE  :: PWPX(:), PWPY(:)
  REAL(SP), ALLOCATABLE  :: UVN(:,:)
  REAL(SP), ALLOCATABLE  :: XFLUX(:,:)
  REAL(SP), ALLOCATABLE  :: VISCOFF(:)

  REAL(SP)  CONV_W(1:KBM1), DISS_W(1:KBM1)
  REAL(SP)  SL_H(-1:KB), W_TMP(0:KB+1)
  REAL(SP)  SL_U, SL_F, TEMP1

  ALLOCATE(U_TMP(0:NT,KB))
  ALLOCATE(V_TMP(0:NT,KB))

  ALLOCATE(PWPX(0:MT))
  ALLOCATE(PWPY(0:MT))

  ALLOCATE(UVN(3*NT,KB))
  ALLOCATE(VISCOFF(0:MT))
  ALLOCATE(XFLUX(0:MT,KB))       ;  XFLUX  = 0.0_SP

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

! step1: calculate the advection and diffusion term
  DO I=1, NT
    DO K=2, KBM1
      U_TMP(I,K) = (U(I,K-1)*DZ1(I,K)+U(I,K)*DZ1(I,K-1))/(DZ1(I,K-1)+DZ1(I,K))
      V_TMP(I,K) = (V(I,K-1)*DZ1(I,K)+V(I,K)*DZ1(I,K-1))/(DZ1(I,K-1)+DZ1(I,K))
    ENDDO
  ENDDO

  DO I=1, NCV
    I1 = NTRG(I)
    DO K=2, KBM1
      UVN(I,K) = V_TMP(I1,K)*DLTXE(I) - U_TMP(I1,K)*DLTYE(I)
    ENDDO
  ENDDO

  DO K=2, KBM1        ! main loop
    PWPX = 0.0_SP
    PWPY = 0.0_SP
    DO I=1,M
      DO J=1,NTSN(I)-1
        I1 = NBSN(I,J)
        I2 = NBSN(I,J+1)
        F1 = 0.5_SP*(W4ZT(I1,K)+W4ZT(I2,K))

        PWPX(I) = PWPX(I) + F1*(VY(I1)-VY(I2))
        PWPY(I) = PWPY(I) + F1*(VX(I2)-VX(I1))
      END DO
      PWPX(I) = PWPX(I)/ART2(I)
      PWPY(I) = PWPY(I)/ART2(I)

      VISCOFF(I) = (VISCOFH(I,K-1)*DZ(I,K)+VISCOFH(I,K)*DZ(I,K-1))/(DZ(I,K-1)+DZ(I,K))
    END DO

    DO I=1,NCV_I
      I1=NTRG(I)
      IA=NIEC(I,1)
      IB=NIEC(I,2)
      XI=0.5_SP*(XIJE(I,1)+XIJE(I,2))
      YI=0.5_SP*(YIJE(I,1)+YIJE(I,2))

      DXA=XI-VX(IA)
      DYA=YI-VY(IA)
      DXB=XI-VX(IB)
      DYB=YI-VY(IB)

      FIJ1=W4ZT(IA,K)+DXA*PWPX(IA)+DYA*PWPY(IA)
      FIJ2=W4ZT(IB,K)+DXB*PWPX(IB)+DYB*PWPY(IB)

      W1MIN=MINVAL(W4ZT(NBSN(IA,1:NTSN(IA)-1),K))
      W1MIN=MIN(W1MIN, W4ZT(IA,K))
      W1MAX=MAXVAL(W4ZT(NBSN(IA,1:NTSN(IA)-1),K))
      W1MAX=MAX(W1MAX, W4ZT(IA,K))
      W2MIN=MINVAL(W4ZT(NBSN(IB,1:NTSN(IB)-1),K))
      W2MIN=MIN(W2MIN, W4ZT(IB,K))
      W2MAX=MAXVAL(W4ZT(NBSN(IB,1:NTSN(IB)-1),K))
      W2MAX=MAX(W2MAX, W4ZT(IB,K))
      IF(FIJ1 < W1MIN) FIJ1=W1MIN
      IF(FIJ1 > W1MAX) FIJ1=W1MAX
      IF(FIJ2 < W2MIN) FIJ2=W2MIN
      IF(FIJ2 > W2MAX) FIJ2=W2MAX

      UN=UVN(I,K)

      VISCOF=(FACT*0.5_SP*(VISCOFF(IA)*NN_HVC(IA)+VISCOFF(IB)*NN_HVC(IB)) + FM1*0.5_SP*(NN_HVC(IA)+NN_HVC(IB)))/HPRNU

      TXX=0.5_SP*(PWPX(IA)+PWPX(IB))*VISCOF
      TYY=0.5_SP*(PWPY(IA)+PWPY(IB))*VISCOF

      FXX=-DT1(I1)*DZZ1(I1,K-1)*TXX*DLTYE(I)
      FYY= DT1(I1)*DZZ1(I1,K-1)*TYY*DLTXE(I)

      EXFLUX=-UN*DT1(I1)*DZZ1(I1,K-1)* &
         ((1.0_SP+SIGN(1.0_SP,UN))*FIJ2+(1.0_SP-SIGN(1.0_SP,UN))*FIJ1)*0.5_SP+FXX+FYY

      XFLUX(IA,K)=XFLUX(IA,K)+EXFLUX
      XFLUX(IB,K)=XFLUX(IB,K)-EXFLUX
    END DO

  ENDDO
# if defined (MULTIPROCESSOR)
  IF(PAR) CALL NODE_MATCH(0,NBN,BN_MLT,BN_LOC,BNC,MT,KB,MYID,NPROCS,XFLUX)
# endif

! step2: calculate vertical flux term
  DO I=1, M
#   if defined (WET_DRY)
    IF(ISWETN(I)*ISWETNT(I) ==1 ) THEN
#   endif 
      SL_H(-1) = DZ(I,1)
      SL_H(0)  = DZ(I,1)
      W_TMP(0) = -W4ZT(I,1)
      DO K=1, KBM1-1
        SL_H(K) = DZZ(I,K)
        W_TMP(K) = W4ZT(I,K)
      ENDDO
      SL_H(KBM1) = DZ(I,KBM1)
      SL_H(KB) = DZ(I,KBM1)
      W_TMP(KBM1) = W4ZT(I,KBM1)
      W_TMP(KB) = W4ZT(I,KB)
      W_TMP(KB+1) = -W4ZT(I,KB)

      DO K=1, KBM1
        CONV_W(K) = 0.25_SP*(WTS(I,K)+WTS(I,K+1))*(W4ZT(I,K)+W4ZT(I,K+1))
        SL_U = 2.0_SP*(W_TMP(K+1)-W_TMP(K+2))/(SL_H(K)+SL_H(K+1))
        SL_F = 2.0_SP*(W_TMP(K-1)-W_TMP(K))/(SL_H(K-2)+SL_H(K-1))
        DISS_W(K) = 0.25_SP*ABS(WTS(I,K)+WTS(I,K+1))*(W_TMP(K)-W_TMP(K+1)-0.5_SP*LIMLED2(SL_U,SL_F,2.0_SP)*(SL_H(K-1)+SL_H(K)))
      ENDDO
        
      DO K=2, KBM1
        TEMP1 = CONV_W(K-1)-CONV_W(K)+DISS_W(K)-DISS_W(K-1)
        XFLUX(I,K) = XFLUX(I,K) + ( TEMP1 + NHQDRZ(I,K))*ART1(I)
      ENDDO
#   if defined (WET_DRY)
    ENDIF
#   endif
  ENDDO

!  step3: update w
  DO I=1, M
#   if defined (WET_DRY)
    IF(ISWETN(I)*ISWETNT(I) == 1) THEN
#   endif
      DO K=2, KBM1
        W4ZT(I,K) = (W4ZT(I,K)-XFLUX(I,K)/ART1(I)/DZZ(I,K-1)*(DTI/DT(I)))*(DT(I)/D(I))
      ENDDO
#   if defined (WET_DRY)
    ENDIF
#   endif
  ENDDO

  DEALLOCATE(U_TMP,V_TMP,PWPX,PWPY,UVN,VISCOFF,XFLUX)
  RETURN
  END SUBROUTINE ADV_W

  SUBROUTINE VDIF_W

  USE LIMS
  USE CONTROL
  USE ALL_VARS
# if defined(MULTIPROCESSOR)
  USE MOD_PAR
# endif
  IMPLICIT NONE

  INTEGER  I, I1, I2, J, K, KI, ITMP
  REAL(SP) A(KBM1-1), B(KBM1-1), C(KBM1-1), R(KBM1-1), RU(KBM1-1)
  REAL(SP) F1, AVG_UTMP, AVG_VTMP
  REAL(SP) A1AVG_TMP, A2AVG_TMP

  CALL A1A2(EL,D)

# if defined(MULTIPROCESSOR)
  IF(PAR) CALL AEXCHANGE(EC,MYID,NPROCS,UF,VF)
# endif
  AVG_U  = 0.0_SP
  AVG_V  = 0.0_SP
  DO K=1, KBM1
    DO I=1, M
      DO J=1, NTVE(I)
        AVG_U(I,K) = AVG_U(I,K) + UF(NBVE(I,J),K)*ART(NBVE(I,J))
        AVG_V(I,K) = AVG_V(I,K) + VF(NBVE(I,J),K)*ART(NBVE(I,J))
      ENDDO
      AVG_U(I,K) = AVG_U(I,K)/ART2(I)
      AVG_V(I,K) = AVG_V(I,K)/ART2(I)
    ENDDO
  ENDDO

  DO I=1, M
!   W4ZT(I,1) = (EL(I)-ET(I))/DTI + AVG_U(I,1)*AVELX(I) + AVG_V(I,1)*AVELY(I)
    W4ZT(I,KB) = -AVG_U(I,KBM1)*A1AVGD(I) - AVG_V(I,KBM1)*A2AVGD(I)
  ENDDO

  DO K=2, KBM1
    DO I=1, M
      AVG_UTMP  = (AVG_U(I,K-1)*DZ(I,K)+AVG_U(I,K)*DZ(I,K-1))/(DZ(I,K-1)+DZ(I,K))
      AVG_VTMP  = (AVG_V(I,K-1)*DZ(I,K)+AVG_V(I,K)*DZ(I,K-1))/(DZ(I,K-1)+DZ(I,K))
      A1AVG_TMP = (A1AVG(I,K-1)*DZ(I,K)+A1AVG(I,K)*DZ(I,K-1))/(DZ(I,K-1)+DZ(I,K))
      A2AVG_TMP = (A2AVG(I,K-1)*DZ(I,K)+A2AVG(I,K)*DZ(I,K-1))/(DZ(I,K-1)+DZ(I,K))
      W4Z(I,K) = W4ZT(I,K) + AVG_UTMP*A1AVG_TMP + AVG_VTMP*A2AVG_TMP - (1.0_SP+Z(I,K))*(EL(I)-ET(I))/DTI
    ENDDO
  ENDDO

  DO I=1, M
#   if defined (WET_DRY)
    IF(ISWETN(I)*ISWETNT(I) ==1) THEN
#   endif
    A(1) = 0.0_SP
    C(1) = -DTI * 0.5_SP*(KM(I,2)+KM(I,3)+2.0_SP*UMOL)/(DZZ(I,1)*DZ(I,2)*D(I)*D(I))
    B(1) = 1.0_SP - C(1) + DTI*0.5_SP*(KM(I,1)+KM(I,2)+2.0_SP*UMOL)/(DZZ(I,1)*DZ(I,1)*D(I)*D(I))
    R(1) = W4Z(I,2) !+ DTI*0.5_SP*(KM(I,1)+KM(I,2)+2.0_SP*UMOL)/(DZZ(I,1)*DZ(I,1)*D(I)*D(I))*W4ZT(I,1)

    DO K=2, KBM1-2
      A(K) = -DTI * 0.5_SP*(KM(I,K)+KM(I,K+1)+2.0_SP*UMOL)/(DZZ(I,K)*DZ(I,K)*D(I)*D(I))
      C(K) = -DTI * 0.5_SP*(KM(I,K+1)+KM(I,K+2)+2.0_SP*UMOL)/(DZZ(I,K)*DZ(I,K+1)*D(I)*D(I))
      B(K) = 1.0_SP - A(K) - C(K)
      R(K) = W4Z(I,K+1)
    ENDDO

    A(KBM1-1) = -DTI * 0.5_SP*(KM(I,KBM1-1)+KM(I,KBM1)+2.0_SP*UMOL)/(DZZ(I,KBM1-1)*DZ(I,KBM1-1)*D(I)*D(I))
    C(KBM1-1) = 0.0_SP
    B(KBM1-1) = 1.0_SP - A(KBM1-1) + DTI*0.5_SP*(KM(I,KBM1)+KM(I,KB)+2.0_SP*UMOL)/(DZZ(I,KBM1-1)*DZ(I,KBM1)*D(I)*D(I))
    R(KBM1-1) = W4Z(I,KBM1) !+ DTI*0.5_SP*(KM(I,KBM1)+KM(I,KB)+2.0_SP*UMOL)/(DZZ(I,KBM1-1)*DZ(I,KBM1)*D(I)*D(I))*W4ZT(I,KB)

    CALL TRIDAG(A,B,C,R,RU,KBM1-1)

    DO K=2, KBM1
      W4Z(I,K) = RU(K-1)
    ENDDO
#   if defined (WET_DRY)
    ENDIF
#   endif
  ENDDO

  DO I=1, M
#   if defined (WET_DRY)
    IF(ISWETN(I)*ISWETNT(I)==1) THEN
#   endif
    DO K=2, KBM1
      AVG_UTMP  = (AVG_U(I,K-1)*DZ(I,K)+AVG_U(I,K)*DZ(I,K-1))/(DZ(I,K-1)+DZ(I,K))
      AVG_VTMP  = (AVG_V(I,K-1)*DZ(I,K)+AVG_V(I,K)*DZ(I,K-1))/(DZ(I,K-1)+DZ(I,K))
      A1AVG_TMP = (A1AVG(I,K-1)*DZ(I,K)+A1AVG(I,K)*DZ(I,K-1))/(DZ(I,K-1)+DZ(I,K))
      A2AVG_TMP = (A2AVG(I,K-1)*DZ(I,K)+A2AVG(I,K)*DZ(I,K-1))/(DZ(I,K-1)+DZ(I,K))
      W4ZT(I,K) = W4Z(I,K) - AVG_UTMP*A1AVG_TMP - AVG_VTMP*A2AVG_TMP + (1.0_SP+Z(I,K))*(EL(I)-ET(I))/DTI
    ENDDO
#   if defined (WET_DRY)
    ELSE
      W4ZT(I,1:KB) = 0.0_SP
    ENDIF
#   endif
  ENDDO

  CALL SET_W_OBC

# if defined(MULTIPROCESSOR)
  IF(PAR) CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,KB,MYID,NPROCS,W4ZT)
# endif

  RETURN
  END SUBROUTINE VDIF_W

  SUBROUTINE GEN_POISSON
  !============================================================================
  !   setup poisson discretization for non-hydrostatic pressure
  !   
  !   basic function:  Ax=b
  !============================================================================

  USE LIMS
  USE ALL_VARS
  USE MOD_OBCS
# if defined(MULTIPROCESSOR)
  USE MOD_PAR , ONLY :EMAP, NMAP, NODE_MATCH, NBN, BN_MLT, BN_LOC, BNC
# endif
  USE MOD_PETSc, ONLY: PDEBUG, PUNIT, BLOCAL, A, B, ALO_2_PLO, PSIZE, L2G, USE_LAST, PETSc_SOLVER 
  IMPLICIT NONE

  INTEGER  I, I1, I2, J, J1, J2, JJ, K, LN, IA, IB, PETSc_POS, CHK_POS, NODE, NREC
  INTEGER  TMP2, TMP3, PROW1, PROW2, PROW, PCOL1, PCOL2, PCOL3
  REAL(SP), ALLOCATABLE  ::  XFL_EL(:),YFL_EL(:),XFL_D(:,:),YFL_D(:,:)
  REAL(SP), ALLOCATABLE  ::  XFLUX(:,:), XFLUX_EXT(:,:)
  REAL(SP), ALLOCATABLE  ::  CA(:,:), CAD(:)
  REAL(SP)  :: AVG_U_TMP, AVG_V_TMP, U1, U2, V1, V2
  REAL(SP)  :: DIJ, UIJ, VIJ, EXFLUX, DX, DY, V_TMP1, DD
  REAL(SP)  :: AA1, AB1, AA2, AB2, AA3, AB3, BA, DL, FAI, TMP4, F1
  REAL(SP)  :: DZ_1K, DZ_K, DZ_K1, DZ_1KK, DZ_KK1

  PetscReal :: STERM, VCOL1, VCOL2
  PetscReal :: VAL1_(24), VAL2_(24)
  PetscInt  :: COL1_(24), COL2_(24)
  PetscInt  :: IERR
  CHARACTER(LEN=20) :: SUBNAME

  !============================================================================
  ! basic function:  set b
  ! set local b vector with source term of poisson problem to insert we 
  ! transform from (i,k) in application local order to petsc local order.
  ! Note, in transformation, we consider k the faster running index as we
  ! want nearby diagonals to come from d^2/dz^2 discretization
  ! index argument is also in C counting style (from 0)
  !============================================================================
  SUBNAME = 'PETSc_setRHS'
  IF(PDEBUG) WRITE(PUNIT,*) 'STARTING: ', TRIM(SUBNAME)

  ALLOCATE(XFL_EL(0:MT))        ; XFL_EL   = 0.0_SP
  ALLOCATE(YFL_EL(0:MT))        ; YFL_EL   = 0.0_SP
  ALLOCATE(XFL_D(0:MT,KB))      ; XFL_D    = 0.0_SP
  ALLOCATE(YFL_D(0:MT,KB))      ; YFL_D    = 0.0_SP
  DO I=1,NCV
    I1=NTRG(I)
    IA=NIEC(I,1)
    IB=NIEC(I,2)
    DIJ=(D(NV(I1,1))+D(NV(I1,2))+D(NV(I1,3)))/3.0_SP
    DO K=1, KB
      XFL_D(IA,K) = XFL_D(IA,K) + DIJ*Z1(I1,K)*DLTYE(I)
      XFL_D(IB,K) = XFL_D(IB,K) - DIJ*Z1(I1,K)*DLTYE(I)
      YFL_D(IA,K) = YFL_D(IA,K) - DIJ*Z1(I1,K)*DLTXE(I)
      YFL_D(IB,K) = YFL_D(IB,K) + DIJ*Z1(I1,K)*DLTXE(I)

      IF(ISONB(IA)>=1 .AND. ISONB(IB)>=1) THEN
      DO J1=1, 3
        IF(NV(I1,J1)/=IA .AND. NV(I1,J1)/=IB) THEN
          IF(NBE(I1,J1)==0) THEN

            DX = VX(IB) - VX(IA)
            DY = VY(IB) - VY(IA)
            XFL_D(IA,K) = XFL_D(IA,K) + 0.25_SP*(-DY)*(D(IA)*Z(IA,K)+D(IB)*Z(IB,K))
            XFL_D(IB,K) = XFL_D(IB,K) + 0.25_SP*(-DY)*(D(IA)*Z(IA,K)+D(IB)*Z(IB,K))
            YFL_D(IA,K) = XFL_D(IA,K) + 0.25_SP*DX*(D(IA)*Z(IA,K)+D(IB)*Z(IB,K))
            YFL_D(IB,K) = XFL_D(IB,K) + 0.25_SP*DX*(D(IA)*Z(IA,K)+D(IB)*Z(IB,K))

          ENDIF
        ENDIF
      ENDDO
      ENDIF
    ENDDO

    DIJ=(EL(NV(I1,1))+EL(NV(I1,2))+EL(NV(I1,3)))/3.0_SP
    XFL_EL(IA) = XFL_EL(IA) + DIJ*DLTYE(I)
    XFL_EL(IB) = XFL_EL(IB) - DIJ*DLTYE(I)
    YFL_EL(IA) = YFL_EL(IA) - DIJ*DLTXE(I)
    YFL_EL(IB) = YFL_EL(IB) + DIJ*DLTXE(I)

    IF(ISONB(IA)>=1 .AND. ISONB(IB)>=1) THEN
    DO J1=1, 3
      IF(NV(I1,J1)/=IA .AND. NV(I1,J1)/=IB) THEN
         IF(NBE(I1,J1)==0) THEN

            DX = VX(IB) - VX(IA)
            DY = VY(IB) - VY(IA)
            XFL_EL(IA) = XFL_EL(IA) + 0.25_SP*(-DY)*(EL(IA)+EL(IB))
            XFL_EL(IB) = XFL_EL(IB) + 0.25_SP*(-DY)*(EL(IA)+EL(IB))
            YFL_EL(IA) = XFL_EL(IA) + 0.25_SP*DX*(EL(IA)+EL(IB))
            YFL_EL(IB) = XFL_EL(IB) + 0.25_SP*DX*(EL(IA)+EL(IB))

         ENDIF
      ENDIF
    ENDDO
    ENDIF

  ENDDO

  ALLOCATE(XFLUX(0:MT,KBM1)) ;  XFLUX = 0.0_SP
  DO I=1,NCV
    I1=NTRG(I)
    IA=NIEC(I,1)
    IB=NIEC(I,2)
    DIJ=(D(NV(I1,1))+D(NV(I1,2))+D(NV(I1,3)))/3.0_SP
    DO K=1,KBM1
      UIJ=UF(I1,K)
      VIJ=VF(I1,K)
      EXFLUX=DIJ*DZ1(I1,K)*(-UIJ*DLTYE(I)+VIJ*DLTXE(I))
      XFLUX(IA,K)=XFLUX(IA,K)-EXFLUX
      XFLUX(IB,K)=XFLUX(IB,K)+EXFLUX
    END DO
  ENDDO

! external source terms
  ALLOCATE(XFLUX_EXT(0:MT,KBM1)) ;  XFLUX_EXT = 0.0_SP
  IF(NUMQBC >= 1) THEN
    IF(RIVER_INFLOW_LOCATION == 'node') THEN
      DO J=1,NUMQBC
        JJ=INODEQ(J)
        DO K=1,KBM1
          XFLUX_EXT(JJ,K)=XFLUX_EXT(JJ,K)-QDIS(J)*VQDIST(J,K)
        END DO
      END DO
    ELSE IF(RIVER_INFLOW_LOCATION == 'edge') THEN
      DO J=1,NUMQBC
        J1=N_ICELLQ(J,1)
        J2=N_ICELLQ(J,2)
        DO K=1,KBM1
          XFLUX_EXT(J1,K)=XFLUX_EXT(J1,K)-QDIS(J)*RDISQ(J,1)*VQDIST(J,K)
          XFLUX_EXT(J2,K)=XFLUX_EXT(J2,K)-QDIS(J)*RDISQ(J,2)*VQDIST(J,K)
        END DO
      END DO
    END IF
  END IF
! end add source terms

  CALL VecSet(BLOCAL,ZERO,IERR);CHKERRQ(IERR)
  CALL VecSet(B,ZERO,IERR);CHKERRQ(IERR)
  RHS = 0.0_SP
 
  DO I=1, M

    PETSc_POS = ALO_2_PLO(I,1)
    IF( PETSc_POS > PSIZE ) CYCLE
 
    DO K=1, KBM1
      
      IF(K==1) THEN
        AVG_U_TMP = (-ZZ(I,2))*(-ZZ(I,3))/(ZZ(I,1)-ZZ(I,2))/(ZZ(I,1)-ZZ(I,3))*AVG_U(I,1)+(-ZZ(I,1))*(-ZZ(I,3))/(ZZ(I,2)-ZZ(I,1))/(ZZ(I,2)-ZZ(I,3))*AVG_U(I,2)+(-ZZ(I,1))*(-ZZ(I,2))/(ZZ(I,3)-ZZ(I,1))/(ZZ(I,3)-ZZ(I,2))*AVG_U(I,3)
        AVG_V_TMP = (-ZZ(I,2))*(-ZZ(I,3))/(ZZ(I,1)-ZZ(I,2))/(ZZ(I,1)-ZZ(I,3))*AVG_V(I,1)+(-ZZ(I,1))*(-ZZ(I,3))/(ZZ(I,2)-ZZ(I,1))/(ZZ(I,2)-ZZ(I,3))*AVG_V(I,2)+(-ZZ(I,1))*(-ZZ(I,2))/(ZZ(I,3)-ZZ(I,1))/(ZZ(I,3)-ZZ(I,2))*AVG_V(I,3)

        U2 = (AVG_U(I,1)*DZ(I,2)+AVG_U(I,2)*DZ(I,1))/(DZ(I,1)+DZ(I,2))
        V2 = (AVG_V(I,1)*DZ(I,2)+AVG_V(I,2)*DZ(I,1))/(DZ(I,1)+DZ(I,2))

        V_TMP1 = ((U2-AVG_U_TMP)*XFL_EL(I)+U2*XFL_D(I,2)-AVG_U_TMP*XFL_D(I,1)) + ((V2-AVG_V_TMP)*YFL_EL(I)+V2*YFL_D(I,2)-AVG_V_TMP*YFL_D(I,1))

      ELSE IF(K==KBM1) THEN

        AVG_U_TMP = (-1.0D0-ZZ(I,KBM1-1))*(-1.0D0-ZZ(I,KBM1-2))/(ZZ(I,KBM1)-ZZ(I,KBM1-1))/(ZZ(I,KBM1)-ZZ(I,KBM1-2))*AVG_U(I,KBM1)+(-1.0D0-ZZ(I,KBM1))*(-1.0D0-ZZ(I,KBM1-2))/(ZZ(I,KBM1-1)-ZZ(I,KBM1))/(ZZ(I,KBM1-1)-ZZ(I,KBM1-2))*AVG_U(I,KBM1-1)+(-1.0D0-ZZ(I,KBM1))*(-1.0D0-ZZ(I,KBM1-1))/(ZZ(I,KBM1-2)-ZZ(I,KBM1))/(ZZ(I,KBM1-2)-ZZ(I,KBM1-1))*AVG_U(I,KBM1-2)
        AVG_V_TMP = (-1.0D0-ZZ(I,KBM1-1))*(-1.0D0-ZZ(I,KBM1-2))/(ZZ(I,KBM1)-ZZ(I,KBM1-1))/(ZZ(I,KBM1)-ZZ(I,KBM1-2))*AVG_V(I,KBM1)+(-1.0D0-ZZ(I,KBM1))*(-1.0D0-ZZ(I,KBM1-2))/(ZZ(I,KBM1-1)-ZZ(I,KBM1))/(ZZ(I,KBM1-1)-ZZ(I,KBM1-2))*AVG_V(I,KBM1-1)+(-1.0D0-ZZ(I,KBM1))*(-1.0D0-ZZ(I,KBM1-1))/(ZZ(I,KBM1-2)-ZZ(I,KBM1))/(ZZ(I,KBM1-2)-ZZ(I,KBM1-1))*AVG_V(I,KBM1-2)

        U1 = (AVG_U(I,KBM1-1)*DZ(I,KBM1)+AVG_U(I,KBM1)*DZ(I,KBM1-1))/(DZ(I,KBM1-1)+DZ(I,KBM1))
        V1 = (AVG_V(I,KBM1-1)*DZ(I,KBM1)+AVG_V(I,KBM1)*DZ(I,KBM1-1))/(DZ(I,KBM1-1)+DZ(I,KBM1))

        V_TMP1 = ((AVG_U_TMP-U1)*XFL_EL(I)+AVG_U_TMP*XFL_D(I,KB)-U1*XFL_D(I,KBM1)) + ((AVG_V_TMP-V1)*YFL_EL(I)+AVG_V_TMP*YFL_D(I,KB)-V1*YFL_D(I,KBM1))

      ELSE

        U1 = (AVG_U(I,K-1)*DZ(I,K)+AVG_U(I,K)*DZ(I,K-1))/(DZ(I,K-1)+DZ(I,K))
        V1 = (AVG_V(I,K-1)*DZ(I,K)+AVG_V(I,K)*DZ(I,K-1))/(DZ(I,K-1)+DZ(I,K))
        U2 = (AVG_U(I,K)*DZ(I,K+1)+AVG_U(I,K+1)*DZ(I,K))/(DZ(I,K)+DZ(I,K+1))
        V2 = (AVG_V(I,K)*DZ(I,K+1)+AVG_V(I,K+1)*DZ(I,K))/(DZ(I,K)+DZ(I,K+1))

        V_TMP1 = ((U2-U1)*XFL_EL(I)+U2*XFL_D(I,K+1)-U1*XFL_D(I,K)) + ((V2-V1)*YFL_EL(I)+V2*YFL_D(I,K+1)-V1*YFL_D(I,K))
      ENDIF

      V_TMP1   = V_TMP1 + (W4ZT(I,K) - W4ZT(I,K+1))*ART1(I)
      RHS(I,K) = 1.0E3_SP*(XFLUX(I,K) + V_TMP1)/DTI
      
    ENDDO

    RHS(I,1) =  (ZZ(I,1)-ZZ(I,3))*(ZZ(I,1)-ZZ(I,4))/(ZZ(I,2)-ZZ(I,3))/(ZZ(I,2)-ZZ(I,4))*RHS(I,2)+(ZZ(I,1)-ZZ(I,2))*(ZZ(I,1)-ZZ(I,4))/(ZZ(I,3)-ZZ(I,2))/(ZZ(I,3)-ZZ(I,4))*RHS(I,3)+(ZZ(I,1)-ZZ(I,2))*(ZZ(I,1)-ZZ(I,3))/(ZZ(I,4)-ZZ(I,2))/(ZZ(I,4)-ZZ(I,3))*RHS(I,4)

    DO K=1, KBM1
      RHS(I,K) = RHS(I,K) + 1.0E3_SP*XFLUX_EXT(I,K)/DTI
    ENDDO

  ENDDO

# if defined(MULTIPROCESSOR)
   IF(PAR) CALL NODE_MATCH(0,NBN,BN_MLT,BN_LOC,BNC,MT,KBM1,MYID,NPROCS,RHS)
# endif

  DO I=1, M
    DO K=1, KBM1

      PETSc_POS = ALO_2_PLO(I,K)
      IF( PETSc_POS > PSIZE ) CYCLE
#     if !defined (WET_DRY)
      IF(ISONB(I)/=2 ) THEN
#     else
      IF(ISONB(I)/=2 .AND. ISWETN(I)==1) THEN
#     endif
        STERM = RHS(I,K)
      ELSE
        STERM = 0.0D0
      ENDIF 

      CALL VecSetValues(BLOCAL,1,PETSc_POS-1,STERM,INSERT_VALUES,IERR);CHKERRQ(IERR)

    ENDDO
  ENDDO

  !-----------------------------------------------------------------------------
  !load global vectors using a scatter with the l2g index set  ==> b
  !-----------------------------------------------------------------------------
# if defined (OLD_PETSC)
  CALL VecScatterBegin(BLOCAL,B,INSERT_VALUES,SCATTER_FORWARD,L2G,IERR);CHKERRQ(IERR)
  !could put code here
  call VecScatterEnd(BLOCAL,B,INSERT_VALUES,SCATTER_FORWARD,L2G,IERR);CHKERRQ(IERR)
# else
  CALL VecScatterBegin(L2G,BLOCAL,B,INSERT_VALUES,SCATTER_FORWARD,IERR);CHKERRQ(IERR)
  !could put code here
  call VecScatterEnd(L2G,BLOCAL,B,INSERT_VALUES,SCATTER_FORWARD,IERR);CHKERRQ(IERR)
# endif
  IF(PDEBUG) WRITE(PUNIT,*) 'FINISHING: ',TRIM(SUBNAME)

  !============================================================================
  ! petsc_setA
  ! discretize Poisson equation to calculate factors load into petsc matrix A
  ! let petsc assemble the matrix across partitions
  ! assemble matrix, use local vars (requires that MatSetLocalToGlobalMapping is setup)
  ! note: row and column are C order (starts from 0) see page 51 of petsc manual, 
  ! bottom paragraph we will set entries row by row
  !============================================================================
  SUBNAME = 'PETSc_setA'
  IF(PDEBUG) WRITE(PUNIT,*) 'STARTING: ',TRIM(SUBNAME)  
 
  CALL MatZeroEntries(A,IERR);CHKERRQ(IERR)

  DO I=1, NCV_I
    I1 = NTRG(I)
    IA = NIEC(I,1)
    IB = NIEC(I,2)

!    IF((ISONB(IA)+ISONB(IB))<4) THEN

#   if defined (WET_DRY)
    IF(ISWETC(I1) == 1) THEN
#   endif

    NREC  = 0
    COL1_ = 0
    COL2_ = 0
    VAL1_ = 0.0_SP
    VAL2_ = 0.0_SP

    PROW1 = ALO_2_PLO(IA,1)
    PROW2 = ALO_2_PLO(IB,1)

    DO J=1, 3
      IF(NV(I1,J)==IA) THEN
        TMP2 = NV(I1,J+1-INT((J+1)/4)*3)
        TMP3 = NV(I1,J+2-INT((J+2)/4)*3)
        AA1  = -0.5_SP*D(IA)*(VX(TMP3)-VX(TMP2))/ART(I1)
        AB1  =  0.5_SP*D(IA)*(VY(TMP3)-VY(TMP2))/ART(I1)
        AA2  = -0.5_SP*D(TMP2)*(VX(IA)-VX(TMP3))/ART(I1)
        AB2  =  0.5_SP*D(TMP2)*(VY(IA)-VY(TMP3))/ART(I1)
        AA3  = -0.5_SP*D(TMP3)*(VX(TMP2)-VX(IA))/ART(I1)
        AB3  =  0.5_SP*D(TMP3)*(VY(TMP2)-VY(IA))/ART(I1)
      ENDIF
    ENDDO         

    PCOL1 = ALO_2_PLO(IA,1)
    PCOL2 = ALO_2_PLO(IB,1)
    PCOL3 = ALO_2_PLO(TMP3,1)
! for the first flux term
    NREC  = NREC + 3
    VCOL1 = (-AA1*DLTXE(I)+AB1*DLTYE(I))*DZ(IA,1)
    COL1_(NREC-2) = PCOL1-1
    VAL1_(NREC-2) = VCOL1
    COL2_(NREC-2) = PCOL1-1
    VAL2_(NREC-2) = -VCOL1
 
    VCOL1 = (-AA2*DLTXE(I)+AB2*DLTYE(I))*DZ(IB,1)
    COL1_(NREC-1) = PCOL2-1
    VAL1_(NREC-1) = VCOL1
    COL2_(NREC-1) = PCOL2-1
    VAL2_(NREC-1) = -VCOL1
 
    VCOL1 = (-AA3*DLTXE(I)+AB3*DLTYE(I))*DZ(TMP3,1)
    COL1_(NREC)   = PCOL3-1
    VAL1_(NREC)   = VCOL1
    COL2_(NREC)   = PCOL3-1
    VAL2_(NREC)   = -VCOL1

! for the second flux term
    NREC  = NREC + 3
    PCOL1 = ALO_2_PLO(IA,1)
    PCOL2 = ALO_2_PLO(IB,1)
    PCOL3 = ALO_2_PLO(TMP3,1)
    VCOL1 = -2.0_SP*(-NHA2(I1,1)*DLTXE(I)+NHA1(I1,1)*DLTYE(I))*DZ1(I1,2)/(DZ1(I1,1)+DZ1(I1,2))/3.0_SP
    COL1_(NREC-2) = PCOL1-1
    VAL1_(NREC-2) = VCOL1
    COL2_(NREC-2) = PCOL1-1
    VAL2_(NREC-2) = -VCOL1

    COL1_(NREC-1) = PCOL2-1
    VAL1_(NREC-1) = VCOL1
    COL2_(NREC-1) = PCOL2-1
    VAL2_(NREC-1) = -VCOL1

    COL1_(NREC)   = PCOL3-1
    VAL1_(NREC)   = VCOL1
    COL2_(NREC)   = PCOL3-1
    VAL2_(NREC)   = -VCOL1

    NREC  = NREC + 3
    PCOL1 = ALO_2_PLO(IA,2)
    PCOL2 = ALO_2_PLO(IB,2)
    PCOL3 = ALO_2_PLO(TMP3,2)
    VCOL1 = -2.0_SP*(-NHA2(I1,2)*DLTXE(I)+NHA1(I1,2)*DLTYE(I))*DZ1(I1,1)/(DZ1(I1,1)+DZ1(I1,2))/3.0_SP
    COL1_(NREC-2) = PCOL1-1
    VAL1_(NREC-2) = VCOL1
    COL2_(NREC-2) = PCOL1-1
    VAL2_(NREC-2) = -VCOL1

    COL1_(NREC-1) = PCOL2-1
    VAL1_(NREC-1) = VCOL1
    COL2_(NREC-1) = PCOL2-1
    VAL2_(NREC-1) = -VCOL1
    
    COL1_(NREC)   = PCOL3-1
    VAL1_(NREC)   = VCOL1
    COL2_(NREC)   = PCOL3-1
    VAL2_(NREC)   = -VCOL1

    IF(ISONB(IA)>=1 .AND. ISONB(IB)>=1) THEN
      DO J1=1, 3
        IF(NV(I1,J1)/=IA .AND. NV(I1,J1)/=IB) THEN
          IF(NBE(I1,J1)==0) THEN

! accumulate the first flux term at solid edge
            DX = VX(IB) - VX(IA)
            DY = VY(IB) - VY(IA)
            DL = SQRT(DX**2+DY**2)
            FAI = ATAN2(DY,DX)

            DZ_K  = 0.5_SP*(DZ(IA,1)+DZ(IB,1))
            DZ_K1 = 0.5_SP*(DZ(IA,2)+DZ(IB,2))
            DZ_KK1= 0.5_SP*(DZ(IA,1)+DZ(IB,1)+DZ(IA,2)+DZ(IB,2))

            NREC  = NREC + 4
            PCOL1 = ALO_2_PLO(IA,1)
            PCOL2 = ALO_2_PLO(IB,1)
            VCOL1 = 0.25_SP*(-SIN(FAI)*NHA1(I1,1)+COS(FAI)*NHA2(I1,1))*DL*DZ_K1/DZ_KK1
            COL1_(NREC-3) = PCOL1-1
            VAL1_(NREC-3) = VCOL1
            COL2_(NREC-3) = PCOL1-1
            VAL2_(NREC-3) = VCOL1
            COL1_(NREC-2) = PCOL2-1
            VAL1_(NREC-2) = VCOL1
            COL2_(NREC-2) = PCOL2-1
            VAL2_(NREC-2) = VCOL1

            PCOL1 = ALO_2_PLO(IA,2)
            PCOL2 = ALO_2_PLO(IB,2)
            VCOL1 = 0.25_SP*(-SIN(FAI)*NHA1(I1,2)+COS(FAI)*NHA2(I1,2))*DL*DZ_K/DZ_KK1
            COL1_(NREC-1) = PCOL1-1
            VAL1_(NREC-1) = VCOL1
            COL2_(NREC-1) = PCOL1-1
            VAL2_(NREC-1) = VCOL1
            COL1_(NREC)   = PCOL2-1
            VAL1_(NREC)   = VCOL1
            COL2_(NREC)   = PCOL2-1
            VAL2_(NREC)   = VCOL1

! accumulate the second flux term at solid edge
!           DX = VX(IA) - VX(IB)
!           DY = VY(IA) - VY(IB)

            NREC  = NREC + 4
            PCOL1 = ALO_2_PLO(IA,1)
            PCOL2 = ALO_2_PLO(IB,1)
            VCOL1 = 0.5_SP*(-DX*NHA2(I1,1)+DY*NHA1(I1,1))*DZ_K1/DZ_KK1
            COL1_(NREC-3) = PCOL1-1
            VAL1_(NREC-3) = VCOL1
            COL2_(NREC-3) = PCOL1-1
            VAL2_(NREC-3) = VCOL1
            COL1_(NREC-2) = PCOL2-1
            VAL1_(NREC-2) = VCOL1
            COL2_(NREC-2) = PCOL2-1
            VAL2_(NREC-2) = VCOL1

            PCOL1 = ALO_2_PLO(IA,2)
            PCOL2 = ALO_2_PLO(IB,2)
            VCOL1 = 0.5_SP*(-DX*NHA2(I1,2)+DY*NHA1(I1,2))*DZ_K/DZ_KK1
            COL1_(NREC-1) = PCOL1-1
            VAL1_(NREC-1) = VCOL1
            COL2_(NREC-1) = PCOL1-1
            VAL2_(NREC-1) = VCOL1
            COL1_(NREC)   = PCOL2-1
            VAL1_(NREC)   = VCOL1
            COL2_(NREC)   = PCOL2-1
            VAL2_(NREC)   = VCOL1

          ENDIF
        ENDIF
      ENDDO
    ENDIF

    CALL MatSetValuesLocal(A,1,PROW1-1,NREC,COL1_,VAL1_,ADD_VALUES,IERR);CHKERRQ(IERR)
    CALL MatSetValuesLocal(A,1,PROW2-1,NREC,COL2_,VAL2_,ADD_VALUES,IERR);CHKERRQ(IERR)

    DO K=2, KBM1-1

      NREC  = 0
      COL1_ = 0
      COL2_ = 0
      VAL1_ = 0.0_SP
      VAL2_ = 0.0_SP

      PROW1 = ALO_2_PLO(IA,K)
      PROW2 = ALO_2_PLO(IB,K)

      PCOL1 = ALO_2_PLO(IA,K)
      PCOL2 = ALO_2_PLO(IB,K)
      PCOL3 = ALO_2_PLO(TMP3,K)
! for the first flux term
      NREC  = NREC + 3
      VCOL1 = (-AA1*DLTXE(I)+AB1*DLTYE(I))*DZ(IA,K)
      COL1_(NREC-2) = PCOL1-1
      VAL1_(NREC-2) = VCOL1
      COL2_(NREC-2) = PCOL1-1
      VAL2_(NREC-2) = -VCOL1

      VCOL1 = (-AA2*DLTXE(I)+AB2*DLTYE(I))*DZ(IB,K)
      COL1_(NREC-1) = PCOL2-1
      VAL1_(NREC-1) = VCOL1
      COL2_(NREC-1) = PCOL2-1
      VAL2_(NREC-1) = -VCOL1

      VCOL1 = (-AA3*DLTXE(I)+AB3*DLTYE(I))*DZ(TMP3,K)
      COL1_(NREC)   = PCOL3-1
      VAL1_(NREC)   = VCOL1
      COL2_(NREC)   = PCOL3-1
      VAL2_(NREC)   = -VCOL1

! for the second flux term
      NREC  = NREC + 3
      PCOL1 = ALO_2_PLO(IA,K-1)
      PCOL2 = ALO_2_PLO(IB,K-1)
      PCOL3 = ALO_2_PLO(TMP3,K-1)
      VCOL1 = 2.0_SP*(-NHA2(I1,K-1)*DLTXE(I)+NHA1(I1,K-1)*DLTYE(I))*DZ1(I1,K)/(DZ1(I1,K-1)+DZ1(I1,K))/3.0_SP
      COL1_(NREC-2) = PCOL1-1
      VAL1_(NREC-2) = VCOL1
      COL2_(NREC-2) = PCOL1-1
      VAL2_(NREC-2) = -VCOL1

      COL1_(NREC-1) = PCOL2-1
      VAL1_(NREC-1) = VCOL1
      COL2_(NREC-1) = PCOL2-1
      VAL2_(NREC-1) = -VCOL1

      COL1_(NREC)   = PCOL3-1
      VAL1_(NREC)   = VCOL1
      COL2_(NREC)   = PCOL3-1
      VAL2_(NREC)   = -VCOL1

      NREC  = NREC + 3
      PCOL1 = ALO_2_PLO(IA,K)
      PCOL2 = ALO_2_PLO(IB,K)
      PCOL3 = ALO_2_PLO(TMP3,K)
      VCOL1 = 2.0_SP*(-NHA2(I1,K)*DLTXE(I)+NHA1(I1,K)*DLTYE(I))*DZ1(I1,K-1)/(DZ1(I1,K-1)+DZ1(I1,K))/3.0_SP - 2.0_SP*(-NHA2(I1,K)*DLTXE(I)+NHA1(I1,K)*DLTYE(I))*DZ1(I1,K+1)/(DZ1(I1,K)+DZ1(I1,K+1))/3.0_SP
      COL1_(NREC-2) = PCOL1-1
      VAL1_(NREC-2) = VCOL1
      COL2_(NREC-2) = PCOL1-1
      VAL2_(NREC-2) = -VCOL1

      COL1_(NREC-1) = PCOL2-1
      VAL1_(NREC-1) = VCOL1
      COL2_(NREC-1) = PCOL2-1
      VAL2_(NREC-1) = -VCOL1

      COL1_(NREC)   = PCOL3-1
      VAL1_(NREC)   = VCOL1
      COL2_(NREC)   = PCOL3-1
      VAL2_(NREC)   = -VCOL1

      NREC  = NREC + 3
      PCOL1 = ALO_2_PLO(IA,K+1)
      PCOL2 = ALO_2_PLO(IB,K+1)
      PCOL3 = ALO_2_PLO(TMP3,K+1)
      VCOL1 = -2.0_SP*(-NHA2(I1,K+1)*DLTXE(I)+NHA1(I1,K+1)*DLTYE(I))*DZ1(I1,K)/(DZ1(I1,K)+DZ1(I1,K+1))/3.0_SP
      COL1_(NREC-2) = PCOL1-1
      VAL1_(NREC-2) = VCOL1
      COL2_(NREC-2) = PCOL1-1
      VAL2_(NREC-2) = -VCOL1

      COL1_(NREC-1) = PCOL2-1
      VAL1_(NREC-1) = VCOL1
      COL2_(NREC-1) = PCOL2-1
      VAL2_(NREC-1) = -VCOL1

      COL1_(NREC)   = PCOL3-1
      VAL1_(NREC)   = VCOL1
      COL2_(NREC)   = PCOL3-1
      VAL2_(NREC)   = -VCOL1

      IF(ISONB(IA)>=1 .AND. ISONB(IB)>=1) THEN
        DO J1=1, 3
          IF(NV(I1,J1)/=IA .AND. NV(I1,J1)/=IB) THEN
            IF(NBE(I1,J1)==0) THEN
! accumulate the first flux term at solid edge
              DX = VX(IB) - VX(IA)
              DY = VY(IB) - VY(IA)
              DL = SQRT(DX**2+DY**2)
              FAI = ATAN2(DY,DX)

              DZ_1K = 0.5_SP*(DZ(IA,K-1)+DZ(IB,K-1))
              DZ_K  = 0.5_SP*(DZ(IA,K)+DZ(IB,K))
              DZ_K1 = 0.5_SP*(DZ(IA,K+1)+DZ(IB,K+1))
              DZ_1KK= 0.5_SP*(DZ(IA,K-1)+DZ(IB,K-1)+DZ(IA,K)+DZ(IB,K))
              DZ_KK1= 0.5_SP*(DZ(IA,K)+DZ(IB,K)+DZ(IA,K+1)+DZ(IB,K+1))

              NREC  = NREC + 6
              PCOL1 = ALO_2_PLO(IA,K-1)
              PCOL2 = ALO_2_PLO(IB,K-1)
              PCOL3 = ALO_2_PLO(TMP3,K-1)
              VCOL1 = 0.25_SP*(SIN(FAI)*NHA1(I1,K-1)-COS(FAI)*NHA2(I1,K-1))*DL*DZ_K/DZ_1KK
              COL1_(NREC-5) = PCOL1-1
              VAL1_(NREC-5) = VCOL1
              COL2_(NREC-5) = PCOL1-1
              VAL2_(NREC-5) = VCOL1
              COL1_(NREC-4) = PCOL2-1
              VAL1_(NREC-4) = VCOL1
              COL2_(NREC-4) = PCOL2-1
              VAL2_(NREC-4) = VCOL1

              PCOL1 = ALO_2_PLO(IA,K)
              PCOL2 = ALO_2_PLO(IB,K)
              PCOL3 = ALO_2_PLO(TMP3,K)
              VCOL1 = 0.25_SP*( (SIN(FAI)*NHA1(I1,K)-COS(FAI)*NHA2(I1,K))*DL*DZ_1K/DZ_1KK+(-SIN(FAI)*NHA1(I1,K)+COS(FAI)*NHA2(I1,K))*DL*DZ_K1/DZ_KK1 )
              COL1_(NREC-3) = PCOL1-1
              VAL1_(NREC-3) = VCOL1
              COL2_(NREC-3) = PCOL1-1
              VAL2_(NREC-3) = VCOL1
              COL1_(NREC-2) = PCOL2-1
              VAL1_(NREC-2) = VCOL1
              COL2_(NREC-2) = PCOL2-1
              VAL2_(NREC-2) = VCOL1

              PCOL1 = ALO_2_PLO(IA,K+1)
              PCOL2 = ALO_2_PLO(IB,K+1)
              PCOL3 = ALO_2_PLO(TMP3,K+1)
              VCOL1 = 0.25_SP*(-SIN(FAI)*NHA1(I1,K+1)+COS(FAI)*NHA2(I1,K+1))*DL*DZ_K/DZ_KK1
              COL1_(NREC-1) = PCOL1-1
              VAL1_(NREC-1) = VCOL1
              COL2_(NREC-1) = PCOL1-1
              VAL2_(NREC-1) = VCOL1
              COL1_(NREC)   = PCOL2-1
              VAL1_(NREC)   = VCOL1
              COL2_(NREC)   = PCOL2-1
              VAL2_(NREC)   = VCOL1 

! accumulate the second flux term at solid edge
!             DX = VX(IA) - VX(IB)
!             DY = VY(IA) - VY(IB)

              NREC  = NREC + 6
              PCOL1 = ALO_2_PLO(IA,K-1)
              PCOL2 = ALO_2_PLO(IB,K-1)
              VCOL1 = 0.5_SP*(DX*NHA2(I1,K-1)-DY*NHA1(I1,K-1))*DZ_K/DZ_1KK
              COL1_(NREC-5) = PCOL1-1
              VAL1_(NREC-5) = VCOL1
              COL2_(NREC-5) = PCOL1-1
              VAL2_(NREC-5) = VCOL1
              COL1_(NREC-4) = PCOL2-1
              VAL1_(NREC-4) = VCOL1
              COL2_(NREC-4) = PCOL2-1
              VAL2_(NREC-4) = VCOL1
 
              PCOL1 = ALO_2_PLO(IA,K)
              PCOL2 = ALO_2_PLO(IB,K)
              VCOL1 = 0.5_SP*( (DX*NHA2(I1,K)-DY*NHA1(I1,K))*DZ_1K/DZ_1KK+(-DX*NHA2(I1,K)+DY*NHA1(I1,K))*DZ_K1/DZ_KK1 )
              COL1_(NREC-3) = PCOL1-1
              VAL1_(NREC-3) = VCOL1
              COL2_(NREC-3) = PCOL1-1
              VAL2_(NREC-3) = VCOL1
              COL1_(NREC-2) = PCOL2-1
              VAL1_(NREC-2) = VCOL1
              COL2_(NREC-2) = PCOL2-1
              VAL2_(NREC-2) = VCOL1

              PCOL1 = ALO_2_PLO(IA,K+1)
              PCOL2 = ALO_2_PLO(IB,K+1)
              VCOL1 = 0.5_SP*(-DX*NHA2(I1,K+1)+DY*NHA1(I1,K+1))*DZ_K/DZ_KK1
              COL1_(NREC-1) = PCOL1-1
              VAL1_(NREC-1) = VCOL1
              COL2_(NREC-1) = PCOL1-1
              VAL2_(NREC-1) = VCOL1
              COL1_(NREC)   = PCOL2-1
              VAL1_(NREC)   = VCOL1
              COL2_(NREC)   = PCOL2-1
              VAL2_(NREC)   = VCOL1
            ENDIF
          ENDIF
        ENDDO
      ENDIF

      CALL MatSetValuesLocal(A,1,PROW1-1,NREC,COL1_,VAL1_,ADD_VALUES,IERR);CHKERRQ(IERR)
      CALL MatSetValuesLocal(A,1,PROW2-1,NREC,COL2_,VAL2_,ADD_VALUES,IERR);CHKERRQ(IERR)

    ENDDO

    NREC  = 0
    COL1_ = 0
    COL2_ = 0
    VAL1_ = 0.0_SP
    VAL2_ = 0.0_SP

    PROW1 = ALO_2_PLO(IA,KBM1)
    PROW2 = ALO_2_PLO(IB,KBM1)

    PCOL1 = ALO_2_PLO(IA,KBM1)
    PCOL2 = ALO_2_PLO(IB,KBM1)
    PCOL3 = ALO_2_PLO(TMP3,KBM1)
! for the first flux term
    NREC  = NREC + 3
    VCOL1 = (-AA1*DLTXE(I)+AB1*DLTYE(I))*DZ(IA,KBM1)
    COL1_(NREC-2) = PCOL1-1
    VAL1_(NREC-2) = VCOL1
    COL2_(NREC-2) = PCOL1-1
    VAL2_(NREC-2) = -VCOL1

    VCOL1 = (-AA2*DLTXE(I)+AB2*DLTYE(I))*DZ(IB,KBM1)
    COL1_(NREC-1) = PCOL2-1
    VAL1_(NREC-1) = VCOL1
    COL2_(NREC-1) = PCOL2-1
    VAL2_(NREC-1) = -VCOL1
 
    VCOL1 = (-AA3*DLTXE(I)+AB3*DLTYE(I))*DZ(TMP3,KBM1)
    COL1_(NREC)   = PCOL3-1
    VAL1_(NREC)   = VCOL1
    COL2_(NREC)   = PCOL3-1
    VAL2_(NREC)   = -VCOL1

! for the second flux term
    NREC  = NREC + 3
    PCOL1 = ALO_2_PLO(IA,KBM1-1)
    PCOL2 = ALO_2_PLO(IB,KBM1-1)
    PCOL3 = ALO_2_PLO(TMP3,KBM1-1)
    VCOL1 = 2.0_SP*(-NHA2(I1,KBM1-1)*DLTXE(I)+NHA1(I1,KBM1-1)*DLTYE(I))*DZ1(I1,KBM1)/(DZ1(I1,KBM1-1)+DZ1(I1,KBM1))/3.0_SP
    COL1_(NREC-2) = PCOL1-1
    VAL1_(NREC-2) = VCOL1
    COL2_(NREC-2) = PCOL1-1
    VAL2_(NREC-2) = -VCOL1

    COL1_(NREC-1) = PCOL2-1
    VAL1_(NREC-1) = VCOL1
    COL2_(NREC-1) = PCOL2-1
    VAL2_(NREC-1) = -VCOL1

    COL1_(NREC)   = PCOL3-1
    VAL1_(NREC)   = VCOL1
    COL2_(NREC)   = PCOL3-1
    VAL2_(NREC)   = -VCOL1

    NREC  = NREC + 3
    PCOL1 = ALO_2_PLO(IA,KBM1)
    PCOL2 = ALO_2_PLO(IB,KBM1)
    PCOL3 = ALO_2_PLO(TMP3,KBM1)
    VCOL1 = 2.0_SP*(-NHA2(I1,KBM1)*DLTXE(I)+NHA1(I1,KBM1)*DLTYE(I))*DZ1(I1,KBM1-1)/(DZ1(I1,KBM1-1)+DZ1(I1,KBM1))/3.0_SP - 2.0_SP*(-NHA2D(I1)*DLTXE(I)+NHA1D(I1)*DLTYE(I))/3.0_SP
    COL1_(NREC-2) = PCOL1-1
    VAL1_(NREC-2) = VCOL1
    COL2_(NREC-2) = PCOL1-1
    VAL2_(NREC-2) = -VCOL1

    COL1_(NREC-1) = PCOL2-1
    VAL1_(NREC-1) = VCOL1
    COL2_(NREC-1) = PCOL2-1
    VAL2_(NREC-1) = -VCOL1

    COL1_(NREC)   = PCOL3-1
    VAL1_(NREC)   = VCOL1
    COL2_(NREC)   = PCOL3-1
    VAL2_(NREC)   = -VCOL1

    IF(ISONB(IA)>=1 .AND. ISONB(IB)>=1) THEN
      DO J1=1, 3
        IF(NV(I1,J1)/=IA .AND. NV(I1,J1)/=IB) THEN
          IF(NBE(I1,J1)==0) THEN

! accumulate the first flux term at solid edge
            DX = VX(IB) - VX(IA)
            DY = VY(IB) - VY(IA)
            DL = SQRT(DX**2+DY**2)
            FAI = ATAN2(DY,DX)

            DZ_1K = 0.5_SP*(DZ(IA,KBM1-1)+DZ(IB,KBM1-1))
            DZ_K  = 0.5_SP*(DZ(IA,KBM1)+DZ(IB,KBM1))
            DZ_1KK= 0.5_SP*(DZ(IA,KBM1-1)+DZ(IB,KBM1-1)+DZ(IA,KBM1)+DZ(IB,KBM1))

            NREC  = NREC + 4
            PCOL1 = ALO_2_PLO(IA,KBM1-1)
            PCOL2 = ALO_2_PLO(IB,KBM1-1)
            VCOL1 = 0.25_SP*(SIN(FAI)*NHA1(I1,KBM1-1)-COS(FAI)*NHA2(I1,KBM1-1))*DL*DZ_K/DZ_1KK
            COL1_(NREC-3) = PCOL1-1
            VAL1_(NREC-3) = VCOL1
            COL2_(NREC-3) = PCOL1-1
            VAL2_(NREC-3) = VCOL1
            COL1_(NREC-2) = PCOL2-1
            VAL1_(NREC-2) = VCOL1
            COL2_(NREC-2) = PCOL2-1
            VAL2_(NREC-2) = VCOL1

            PCOL1 = ALO_2_PLO(IA,KBM1)
            PCOL2 = ALO_2_PLO(IB,KBM1)
            VCOL1 = 0.25_SP*( (SIN(FAI)*NHA1(I1,KBM1)-COS(FAI)*NHA2(I1,KBM1))*DL*DZ_1K/DZ_1KK - SIN(FAI)*NHA1D(I1)*DL + COS(FAI)*NHA2D(I1)*DL )
            COL1_(NREC-1) = PCOL1-1
            VAL1_(NREC-1) = VCOL1
            COL2_(NREC-1) = PCOL1-1
            VAL2_(NREC-1) = VCOL1
            COL1_(NREC)   = PCOL2-1
            VAL1_(NREC)   = VCOL1
            COL2_(NREC)   = PCOL2-1
            VAL2_(NREC)   = VCOL1

! accumulate the second flux term at solid edge
!           DX = VX(IA) - VX(IB)
!           DY = VY(IA) - VY(IB)

            NREC  = NREC + 4
            PCOL1 = ALO_2_PLO(IA,KBM1-1)
            PCOL2 = ALO_2_PLO(IB,KBM1-1)
            VCOL1 = 0.5_SP*(DX*NHA2(I1,KBM1-1)-DY*NHA1(I1,KBM1-1))*DZ_K/DZ_1KK
            COL1_(NREC-3) = PCOL1-1
            VAL1_(NREC-3) = VCOL1
            COL2_(NREC-3) = PCOL1-1
            VAL2_(NREC-3) = VCOL1
            COL1_(NREC-2) = PCOL2-1
            VAL1_(NREC-2) = VCOL1
            COL2_(NREC-2) = PCOL2-1
            VAL2_(NREC-2) = VCOL1

            PCOL1 = ALO_2_PLO(IA,KBM1)
            PCOL2 = ALO_2_PLO(IB,KBM1)
            VCOL1 = 0.5_SP*( (DX*NHA2(I1,KBM1)-DY*NHA1(I1,KBM1))*DZ_1K/DZ_1KK - DX*NHA2D(I1) + DY*NHA1D(I1) )
            COL1_(NREC-1) = PCOL1-1
            VAL1_(NREC-1) = VCOL1
            COL2_(NREC-1) = PCOL1-1
            VAL2_(NREC-1) = VCOL1
            COL1_(NREC)   = PCOL2-1
            VAL1_(NREC)   = VCOL1
            COL2_(NREC)   = PCOL2-1
            VAL2_(NREC)   = VCOL1

          ENDIF
        ENDIF
      ENDDO
    ENDIF

    CALL MatSetValuesLocal(A,1,PROW1-1,NREC,COL1_,VAL1_,ADD_VALUES,IERR);CHKERRQ(IERR)
    CALL MatSetValuesLocal(A,1,PROW2-1,NREC,COL2_,VAL2_,ADD_VALUES,IERR);CHKERRQ(IERR)

#   if defined (WET_DRY)
    ENDIF
#   endif
!    ENDIF
  ENDDO

!  step3:
  ALLOCATE(CA(0:MT,KBM1))       ; CA       = 0.0_SP
  ALLOCATE(CAD(0:MT))           ; CAD      = 0.0_SP
  DO K=1, KBM1
    DO I=1, NCV_I
      I1 = NTRG(I)
      IA = NIEC(I,1)
      IB = NIEC(I,2)
      DX = VX(IA) - VX(IB)
      DY = VY(IA) - VY(IB)
      DD = (D(NV(I1,1))+D(NV(I1,2))+D(NV(I1,3)))/3.0_SP
 
      CA(IA,K) = CA(IA,K) + (-NHA1(I1,K)*DLTYE(I)+NHA2(I1,K)*DLTXE(I))/(DD*DZ1(I1,K))
      CA(IB,K) = CA(IB,K) - (-NHA1(I1,K)*DLTYE(I)+NHA2(I1,K)*DLTXE(I))/(DD*DZ1(I1,K))

      IF(ISONB(IA)>=1 .AND. ISONB(IB)>=1) THEN
        DO J1=1, 3
          IF(NV(I1,J1)/=IA .AND. NV(I1,J1)/=IB) THEN
            IF(NBE(I1,J1)==0) THEN
               CA(IA,K) = CA(IA,K) + (-NHA1(I1,K)*DY+NHA2(I1,K)*DX)/(D(IA)*DZ(IA,K)+D(IB)*DZ(IB,K))
               CA(IB,K) = CA(IB,K) + (-NHA1(I1,K)*DY+NHA2(I1,K)*DX)/(D(IA)*DZ(IA,K)+D(IB)*DZ(IB,K))
            ENDIF
          ENDIF
        ENDDO
      ENDIF

      IF(K==KBM1) THEN
        CAD(IA) = CAD(IA) + (-NHA1D(I1)*DLTYE(I)+NHA2D(I1)*DLTXE(I))/DD
        CAD(IB) = CAD(IB) - (-NHA1D(I1)*DLTYE(I)+NHA2D(I1)*DLTXE(I))/DD
        IF(ISONB(IA)>=1 .AND. ISONB(IB)>=1) THEN
          DO J1=1, 3
            IF(NV(I1,J1)/=IA .AND. NV(I1,J1)/=IB) THEN
              IF(NBE(I1,J1)==0) THEN
                CAD(IA) = CAD(IA) + (-NHA1D(I1)*DY+NHA2D(I1)*DX)/(D(IA)+D(IB))
                CAD(IB) = CAD(IB) + (-NHA1D(I1)*DY+NHA2D(I1)*DX)/(D(IA)+D(IB))
              ENDIF
            ENDIF
          ENDDO
        ENDIF
      ENDIF

    ENDDO
  ENDDO
# if defined(MULTIPROCESSOR)
  IF(PAR) CALL NODE_MATCH(0,NBN,BN_MLT,BN_LOC,BNC,MT,KBM1,MYID,NPROCS,CA)
  IF(PAR) CALL NODE_MATCH(0,NBN,BN_MLT,BN_LOC,BNC,MT,1,MYID,NPROCS,CAD)
# endif

  DO I=1, M

#   if defined (WET_DRY)
    IF(ISWETN(I) == 1) THEN
#   endif

    PROW  = ALO_2_PLO(I,1)
    IF(PROW > PSIZE ) CYCLE

    COL1_ = 0
    VAL1_ = 0.0_SP

    PROW  = ALO_2_PLO(I,1)
    PCOL1 = ALO_2_PLO(I,1)
    VCOL1 = -CA(I,1)*D(I)*DZ(I,1)*DZ(I,2)/(DZ(I,1)+DZ(I,2))
    VCOL1 = VCOL1 + ( 2.0_SP*(AVELX(I)*A1AVG(I,1)+AVELY(I)*A2AVG(I,1))*ART1(I)/D(I)/DZ(I,1) - ( (A1AVG(I,1)*DZ(I,2)+A1AVG(I,2)*DZ(I,1))/(DZ(I,1)+DZ(I,2))*A1AVG(I,1)+(A2AVG(I,1)*DZ(I,2)+A2AVG(I,2)*DZ(I,1))/(DZ(I,1)+DZ(I,2))*A2AVG(I,1) )/D(I)/DZZ(I,1)*ART1(I) )
    VAL1_(1) = VCOL1 + (-2.0_SP*ART1(I)/D(I)/DZ(I,1) - ART1(I)/D(I)/DZZ(I,1) )
    COL1_(1) = PCOL1-1

    PCOL1 = ALO_2_PLO(I,2)
    VCOL1 = -CA(I,2)*D(I)*DZ(I,2)*DZ(I,1)/(DZ(I,1)+DZ(I,2))
    VCOL1 = VCOL1 + ( (A1AVG(I,1)*DZ(I,2)+A1AVG(I,2)*DZ(I,1))/(DZ(I,1)+DZ(I,2))*A1AVG(I,2)+(A2AVG(I,1)*DZ(I,2)+A2AVG(I,2)*DZ(I,1))/(DZ(I,1)+DZ(I,2))*A2AVG(I,2) )/D(I)/DZZ(I,1)*ART1(I)
    VAL1_(2) = VCOL1 + ART1(I)/D(I)/DZZ(I,1)
    COL1_(2) = PCOL1-1
    CALL MatSetValuesLocal(A,1,PROW-1,2,COL1_,VAL1_,ADD_VALUES,IERR);CHKERRQ(IERR)

    DO K=2, KBM1-1

      COL1_ = 0
      VAL1_ = 0.0_SP

      PROW  = ALO_2_PLO(I,K)
      PCOl1 = ALO_2_PLO(I,K-1)
      VCOL1 = CA(I,K-1)*D(I)*DZ(I,K-1)*DZ(I,K)/(DZ(I,K-1)+DZ(I,K))
      VCOL1 = VCOL1 + ( (A1AVG(I,K-1)*DZ(I,K)+A1AVG(I,K)*DZ(I,K-1))/(DZ(I,K-1)+DZ(I,K))*A1AVG(I,K-1)+(A2AVG(I,K-1)*DZ(I,K)+A2AVG(I,K)*DZ(I,K-1))/(DZ(I,K-1)+DZ(I,K))*A2AVG(I,K-1) )/D(I)/DZZ(I,K-1)*ART1(I)
      VAL1_(1) = VCOL1 + ART1(I)/D(I)/DZZ(I,K-1)
      COL1_(1) = PCOL1-1

      PCOL1 = ALO_2_PLO(I,K)
      VCOL1 = CA(I,K)*D(I)*DZ(I,K)*DZ(I,K-1)/(DZ(I,K-1)+DZ(I,K)) - CA(I,K)*D(I)*DZ(I,K)*DZ(I,K+1)/(DZ(I,K)+DZ(I,K+1))
      VCOL1 = VCOL1 + ( -( (A1AVG(I,K-1)*DZ(I,K)+A1AVG(I,K)*DZ(I,K-1))/(DZ(I,K-1)+DZ(I,K))*A1AVG(I,K)+(A2AVG(I,K-1)*DZ(I,K)+A2AVG(I,K)*DZ(I,K-1))/(DZ(I,K-1)+DZ(I,K))*A2AVG(I,K) )/D(I)/DZZ(I,K-1)*ART1(I) - ( (A1AVG(I,K)*DZ(I,K+1)+A1AVG(I,K+1)*DZ(I,K))/(DZ(I,K)+DZ(I,K+1))*A1AVG(I,K)+(A2AVG(I,K)*DZ(I,K+1)+A2AVG(I,K+1)*DZ(I,K))/(DZ(I,K)+DZ(I,K+1))*A2AVG(I,K) )/D(I)/DZZ(I,K)*ART1(I) )
      VAL1_(2) = VCOL1 + ( -ART1(I)/D(I)/DZZ(I,K-1) - ART1(I)/D(I)/DZZ(I,K) )
      COL1_(2) = PCOL1-1

      PCOL1 = ALO_2_PLO(I,K+1)
      VCOL1 = -CA(I,K+1)*D(I)*DZ(I,K+1)*DZ(I,K)/(DZ(I,K)+DZ(I,K+1))
      VCOL1 = VCOL1 + ( (A1AVG(I,K)*DZ(I,K+1)+A1AVG(I,K+1)*DZ(I,K))/(DZ(I,K)+DZ(I,K+1))*A1AVG(I,K+1)+(A2AVG(I,K)*DZ(I,K+1)+A2AVG(I,K+1)*DZ(I,K))/(DZ(I,K)+DZ(I,K+1))*A2AVG(I,K+1) )/D(I)/DZZ(I,K)*ART1(I)
      VAL1_(3) = VCOL1 + ART1(I)/D(I)/DZZ(I,K)
      COL1_(3) = PCOL1-1
      CALL MatSetValuesLocal(A,1,PROW-1,3,COL1_,VAL1_,ADD_VALUES,IERR);CHKERRQ(IERR)
    ENDDO

    COL1_ = 0
    VAL1_ = 0.0_SP
 
    PROW  = ALO_2_PLO(I,KBM1)  
    PCOL1 = ALO_2_PLO(I,KBM1-1)
    VCOL1 = CA(I,KBM1-1)*D(I)*DZ(I,KBM1-1)*DZ(I,KBM1)/(DZ(I,KBM1-1)+DZ(I,KBM1))
    VCOL1 = VCOL1 + ( (A1AVG(I,KBM1-1)*DZ(I,KBM1)+A1AVG(I,KBM1)*DZ(I,KBM1-1))/(DZ(I,KBM1-1)+DZ(I,KBM1))*A1AVG(I,KBM1-1)+(A2AVG(I,KBM1-1)*DZ(I,KBM1)+A2AVG(I,KBM1)*DZ(I,KBM1-1))/(DZ(I,KBM1-1)+DZ(I,KBM1))*A2AVG(I,KBM1-1) )/D(I)/DZZ(I,KBM1-1)*ART1(I)
    VAL1_(1) = VCOL1 + ART1(I)/D(I)/DZZ(I,KBM1-1)
    COL1_(1) = PCOL1-1
 
    PCOL1 = ALO_2_PLO(I,KBM1)
    VCOL1 = CA(I,KBM1)*D(I)*DZ(I,KBM1)*DZ(I,KBM1-1)/(DZ(I,KBM1-1)+DZ(I,KBM1)) - CAD(I)*D(I)
    VCOL1 = VCOL1 + ( -( (A1AVG(I,KBM1-1)*DZ(I,KBM1)+A1AVG(I,KBM1)*DZ(I,KBM1-1))/(DZ(I,KBM1-1)+DZ(I,KBM1))*A1AVG(I,KBM1)+(A2AVG(I,KBM1-1)*DZ(I,KBM1)+A2AVG(I,KBM1)*DZ(I,KBM1-1))/(DZ(I,KBM1-1)+DZ(I,KBM1))*A2AVG(I,KBM1) )/D(I)/DZZ(I,KBM1-1)*ART1(I) + A1AVGD(I)*ART1(I)*A1AVGD(I)/H(I) + A2AVGD(I)*ART1(I)*A2AVGD(I)/H(I) )
    VAL1_(2) = VCOL1 + ( - ART1(I)/D(I)/DZZ(I,KBM1-1) )
    COL1_(2) = PCOL1-1
    CALL MatSetValuesLocal(A,1,PROW-1,2,COL1_,VAL1_,ADD_VALUES,IERR);CHKERRQ(IERR)

    DO LN =1,NTSN(I)-1
      I1 = NBSN(I,LN)
      I2 = NBSN(I,LN+1)

      PCOL1 = ALO_2_PLO(I1,KBM1)
      PCOL2 = ALO_2_PLO(I2,KBM1)
 
      VCOL1 = -(0.5_SP*(VY(I1)-VY(I2))*COS(SITA_GD(I)) + 0.5_SP*(VX(I2)-VX(I1))*SIN(SITA_GD(I)))*ART1(I)*PHPN(I)/ART2(I)
      CALL MatSetValuesLocal(A,1,PROW-1,1,PCOL1-1,VCOL1,ADD_VALUES,IERR);CHKERRQ(IERR)
      CALL MatSetValuesLocal(A,1,PROW-1,1,PCOL2-1,VCOL1,ADD_VALUES,IERR);CHKERRQ(IERR)
    ENDDO

#   if defined (WET_DRY)
    ENDIF
#   endif
  ENDDO
! end step3

  !distribute components where they are needed among partitions
  CALL MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,IERR);CHKERRQ(IERR)
  !--> could put coding here
  CALL MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,IERR);CHKERRQ(IERR)

  CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

# if defined (WET_DRY)
  DO I=1, M
    IF(ISWETN(I) /= 1) THEN

      PROW    = ALO_2_PLO(I,1)
      IF( PROW > PSIZE ) CYCLE

      VCOL1 = 0.0D0 
      DO J=1, NTSN(I)-1
        PROW  = ALO_2_PLO(I,1)

        PCOL1 = ALO_2_PLO(NBSN(I,J),1)
        CALL MatSetValuesLocal(A,1,PROW-1,1,PCOL1-1,VCOL1,INSERT_VALUES,IERR);CHKERRQ(IERR)
        PCOL1 = ALO_2_PLO(NBSN(I,J),2)
        CALL MatSetValuesLocal(A,1,PROW-1,1,PCOL1-1,VCOL1,INSERT_VALUES,IERR);CHKERRQ(IERR)

        DO K=2, KBM1-1
          PROW  = ALO_2_PLO(I,K)

          PCOL1 = ALO_2_PLO(NBSN(I,J),K-1)
          CALL MatSetValuesLocal(A,1,PROW-1,1,PCOL1-1,VCOL1,INSERT_VALUES,IERR);CHKERRQ(IERR)
          PCOL1 = ALO_2_PLO(NBSN(I,J),K)
          CALL MatSetValuesLocal(A,1,PROW-1,1,PCOL1-1,VCOL1,INSERT_VALUES,IERR);CHKERRQ(IERR)
          PCOL1 = ALO_2_PLO(NBSN(I,J),K+1)
          CALL MatSetValuesLocal(A,1,PROW-1,1,PCOL1-1,VCOL1,INSERT_VALUES,IERR);CHKERRQ(IERR)
        ENDDO

        PROW  = ALO_2_PLO(I,KBM1)

        PCOL1 = ALO_2_PLO(NBSN(I,J),KBM1-1)
        CALL MatSetValuesLocal(A,1,PROW-1,1,PCOL1-1,VCOL1,INSERT_VALUES,IERR);CHKERRQ(IERR)
        PCOL1 = ALO_2_PLO(NBSN(I,J),KBM1)
        CALL MatSetValuesLocal(A,1,PROW-1,1,PCOL1-1,VCOL1,INSERT_VALUES,IERR);CHKERRQ(IERR)
      ENDDO

      VCOL1 = 1.0D0
      DO K=1, KBM1
        PROW = ALO_2_PLO(I,K)
        CALL MatSetValuesLocal(A,1,PROW-1,1,PROW-1,VCOL1,INSERT_VALUES,IERR);CHKERRQ(IERR)
      ENDDO

    ENDIF
  ENDDO
# endif


!  DO I=1, IOBCN

!    PROW = ALO_2_PLO(I_OBC_N(I),1)
!    IF(PROW>PSIZE) CYCLE

!    VCOL1 = 0.0D0
!    NODE  = I_OBC_N(I)
!    DO J=1, NTSN(NODE)-1
!      PROW  = ALO_2_PLO(NODE,1)

!      PCOL1 = ALO_2_PLO(NBSN(NODE,J),1)
!      CALL MatSetValuesLocal(A,1,PROW-1,1,PCOL1-1,VCOL1,INSERT_VALUES,IERR);CHKERRQ(IERR)
!      PCOL1 = ALO_2_PLO(NBSN(NODE,J),2)
!      CALL MatSetValuesLocal(A,1,PROW-1,1,PCOL1-1,VCOL1,INSERT_VALUES,IERR);CHKERRQ(IERR)

!      DO K=2, KBM1-1
!        PROW  = ALO_2_PLO(NODE,K)

!        PCOL1 = ALO_2_PLO(NBSN(NODE,J),K-1)
!        CALL MatSetValuesLocal(A,1,PROW-1,1,PCOL1-1,VCOL1,INSERT_VALUES,IERR);CHKERRQ(IERR)
!        PCOL1 = ALO_2_PLO(NBSN(NODE,J),K)
!        CALL MatSetValuesLocal(A,1,PROW-1,1,PCOL1-1,VCOL1,INSERT_VALUES,IERR);CHKERRQ(IERR)
!        PCOL1 = ALO_2_PLO(NBSN(NODE,J),K+1)
!        CALL MatSetValuesLocal(A,1,PROW-1,1,PCOL1-1,VCOL1,INSERT_VALUES,IERR);CHKERRQ(IERR)
!      ENDDO

!      PROW  = ALO_2_PLO(NODE,KBM1)

!      PCOL1 = ALO_2_PLO(NBSN(NODE,J),KBM1-1)
!      CALL MatSetValuesLocal(A,1,PROW-1,1,PCOL1-1,VCOL1,INSERT_VALUES,IERR);CHKERRQ(IERR)
!      PCOL1 = ALO_2_PLO(NBSN(NODE,J),KBM1)
!      CALL MatSetValuesLocal(A,1,PROW-1,1,PCOL1-1,VCOL1,INSERT_VALUES,IERR);CHKERRQ(IERR)
!    ENDDO

!    VCOL1 = 1.0D0
!    DO K=1, KBM1
!      PROW = ALO_2_PLO(NODE,K)
!      CALL MatSetValuesLocal(A,1,PROW-1,1,PROW-1,VCOL1,INSERT_VALUES,IERR);CHKERRQ(IERR)
!    ENDDO

!  ENDDO

  CALL MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,IERR);CHKERRQ(IERR)
  CALL MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,IERR);CHKERRQ(IERR)

  IF(PDEBUG) WRITE(PUNIT,*) 'SET VALUES IN SETA COMPLETE: ', TRIM(SUBNAME)

  !at this point we can view the matrix.
  !options are:
  ! PETSC_VIEWER_STDOUT_WORLD:  view as text the column number and value by row
  !           this option is useful for debugging but not good for big mats
  ! PETSC_VIEWER_DRAW_WORLD:  draw the matrix nonzero pattern to an X11 screen
  !           this option is good for large matrices but you will need to
  !           put a pause or sleep call afterword to see anything if the prog.
  !           exits quickly.
  !
!  CALL PetscViewerASCIIOpen(PETSC_COMM_WORLD,"mat.output",viewer,IERR);CHKERRQ(IERR)
!  CALL MatView(A,viewer,IERR);CHKERRQ(IERR)
!  CALL MatView(A,PETSC_VIEWER_STDOUT_WORLD,IERR)
!  CALL MatView(A,PETSC_VIEWER_DRAW_WORLD,IERR)
!  CALL SYSTEM('sleep 20')

  IF(PDEBUG) WRITE(PUNIT,*) 'FINISHING: ',TRIM(SUBNAME)

  IF(USE_LAST) CALL PETSc_SETICS
  CALL PETSc_SOLVER

  DEALLOCATE(XFL_EL,YFL_EL,XFL_D,YFL_D)
  DEALLOCATE(XFLUX,XFLUX_EXT,CA,CAD)
  RETURN
  END SUBROUTINE GEN_POISSON

  SUBROUTINE UPDATE_QUVW

  USE LIMS
  USE CONTROL
  USE MOD_OBCS
  USE ALL_VARS
# if defined(MULTIPROCESSOR)
  USE MOD_PAR
# endif
  IMPLICIT NONE

  INTEGER  I, J, JJ, K, I1, I2, IA, IB, J1, J2
  REAL(SP) TMP2, F1, DIJ, QIJ, UIJ, VIJ, EXFLUX, DX, DY
  REAL(SP), ALLOCATABLE :: PQPPX(:), PQPPY(:)
  REAL(SP), ALLOCATABLE :: XFLUX(:,:)
  REAL(SP), ALLOCATABLE :: XFL_EL(:),YFL_EL(:),XFL_D(:,:),YFL_D(:,:)
  REAL(SP) AVG_U_TMP, AVG_V_TMP, V_TMP1
  REAL(SP) UTMP, VTMP, U1, U2, V1, V2
  
!==============================================================================|
  ALLOCATE(PQPPX(0:NT))        ; PQPPX   = 0.0_SP
  ALLOCATE(PQPPY(0:NT))        ; PQPPY   = 0.0_SP
  ALLOCATE(XFLUX(MT,KBM1))     ; XFLUX   = 0.0_SP

  CALL PETSc_GETQ
  IF(IINT<=PROJ_SWITCH) THEN
    NHQ = QP
  ELSE
    NHQ = NHQ + QP
  ENDIF

! update next time step pressure gradient
  DO I=1,N
#   if defined (WET_DRY)
    IF(ISWETC(I)==1) THEN
#   endif
    DO K = 1, KBM1
      QP1(I,K)  = ( QP(NV(I,1),K)+ QP(NV(I,2),K)+ QP(NV(I,3),K) )/3.0_SP
    END DO
#   if defined (WET_DRY)
    ELSE
      QP1(I,1:KBM1) = 0.0_SP
    ENDIF
#   endif
  END DO
  IF(IINT<=PROJ_SWITCH) THEN
    NHQ1 = QP1
  ELSE
    NHQ1 = NHQ1 + QP1
  ENDIF

  DO K=1, KBM1
    PQPPX = 0.0_SP
    PQPPY = 0.0_SP
    DO I=1, NE
      IA=IEC(I,1)
      IB=IEC(I,2)
      J1=IENODE(I,1)
      J2=IENODE(I,2)
      DIJ=0.5_SP*(D(J1)*DZ(J1,K)+D(J2)*DZ(J2,K))
      QIJ=0.5_SP*(QP(J1,K)+QP(J2,K))

      PQPPX(IA)=PQPPX(IA)-QIJ*DIJ*DLTYC(I)
      PQPPY(IA)=PQPPY(IA)+QIJ*DIJ*DLTXC(I)
      PQPPX(IB)=PQPPX(IB)+QIJ*DIJ*DLTYC(I)
      PQPPY(IB)=PQPPY(IB)-QIJ*DIJ*DLTXC(I)
    ENDDO

    DO I=1, N

      IF(K==1) THEN
        TMP2 = -(QP1(I,1)*NHA1(I,1)*DZ1(I,2)+QP1(I,2)*NHA1(I,2)*DZ1(I,1))/(DZ1(I,1)+DZ1(I,2))
      ELSE IF(K==KBM1) THEN
        TMP2 = (QP1(I,KBM2)*NHA1(I,KBM2)*DZ1(I,KBM1)+QP1(I,KBM1)*NHA1(I,KBM1)*DZ1(I,KBM2))/(DZ1(I,KBM2)+DZ1(I,KBM1)) - NHA1D(I)*QP1(I,KBM1)
      ELSE
        TMP2 = (QP1(I,K-1)*NHA1(I,K-1)*DZ1(I,K)+QP1(I,K)*NHA1(I,K)*DZ1(I,K-1))/(DZ1(I,K-1)+DZ1(I,K)) - (QP1(I,K)*NHA1(I,K)*DZ1(I,K+1)+QP1(I,K+1)*NHA1(I,K+1)*DZ1(I,K))/(DZ1(I,K)+DZ1(I,K+1))
      ENDIF
      TMP2 = 1.0E-3_SP*(PQPPX(I)+TMP2*ART(I))
      UF(I,K)   = UF(I,K) - (DTI/D1(I)/DZ1(I,K))*TMP2/ART(I)

      IF(K==1) THEN
        TMP2 = -(QP1(I,1)*NHA2(I,1)*DZ1(I,2)+QP1(I,2)*NHA2(I,2)*DZ1(I,1))/(DZ1(I,1)+DZ1(I,2))
      ELSE IF(K==KBM1) THEN
        TMP2 = (QP1(I,KBM2)*NHA2(I,KBM2)*DZ1(I,KBM1)+QP1(I,KBM1)*NHA2(I,KBM1)*DZ1(I,KBM2))/(DZ1(I,KBM2)+DZ1(I,KBM1)) - NHA2D(I)*QP1(I,KBM1)
      ELSE
        TMP2 = (QP1(I,K-1)*NHA2(I,K-1)*DZ1(I,K)+QP1(I,K)*NHA2(I,K)*DZ1(I,K-1))/(DZ1(I,K-1)+DZ1(I,K)) - (QP1(I,K)*NHA2(I,K)*DZ1(I,K+1)+QP1(I,K+1)*NHA2(I,K+1)*DZ1(I,K))/(DZ1(I,K)+DZ1(I,K+1))
      ENDIF
      TMP2 = 1.0E-3_SP*(PQPPY(I)+TMP2*ART(I))
      VF(I,K)   = VF(I,K) - (DTI/D1(I)/DZ1(I,K))*TMP2/ART(I)
    ENDDO

  ENDDO

# if defined (WET_DRY)
  DO I=1, N
    IF(ISWETCT(I)*ISWETC(I) /= 1) THEN
      UF(I,1:KB) = 0.0_SP
      VF(I,1:KB) = 0.0_SP
    ENDIF
  ENDDO
# endif

  CALL SET_UV_OBC

# if !defined (SEMI_IMPLICIT)
# if defined (GCN)
    CALL BCOND_GCN(3,0)    ! Boundary Condition on U/V At River Input           !
# else
    CALL BCOND_GCY(3,0)    ! Boundary Condition on U/V At River Input           !
# endif
# endif

# if defined (WET_DRY)
  DO I=1, NT
    IF(H1(I) <= STATIC_SSH_ADJ ) THEN
      UTMP = SUM(UF(I,1:KBM1)*DZ1(I,1:KBM1))
      VTMP = SUM(VF(I,1:KBM1)*DZ1(I,1:KBM1))
      DO K=1, KBM1
        UF(I,K)=UTMP
        VF(I,K)=VTMP
      ENDDO
    ENDIF
  ENDDO
# endif

# if defined(MULTIPROCESSOR)
  IF(PAR) CALL AEXCHANGE(EC,MYID,NPROCS,UF,VF)
# endif

  CALL UPDATE_ZETA

  NHQ2DX = 0.0_SP
  NHQ2DY = 0.0_SP
  CALL A1A2(EL,D)

  DO K=1, KBM1
    PQPPX = 0.0_SP
    PQPPY = 0.0_SP
    DO I=1, NE
      IA=IEC(I,1)
      IB=IEC(I,2)
      J1=IENODE(I,1)
      J2=IENODE(I,2)
      DIJ=0.5_SP*(D(J1)*DZ(J1,K)+D(J2)*DZ(J2,K))
      QIJ=0.5_SP*(NHQ(J1,K)+NHQ(J2,K))

      PQPPX(IA)=PQPPX(IA)-QIJ*DIJ*DLTYC(I)
      PQPPY(IA)=PQPPY(IA)+QIJ*DIJ*DLTXC(I)
      PQPPX(IB)=PQPPX(IB)+QIJ*DIJ*DLTYC(I)
      PQPPY(IB)=PQPPY(IB)-QIJ*DIJ*DLTXC(I)
    ENDDO

    DO I=1, N

      IF(K==1) THEN
        TMP2 = -(NHQ1(I,1)*NHA1(I,1)*DZ1(I,2)+NHQ1(I,2)*NHA1(I,2)*DZ1(I,1))/(DZ1(I,1)+DZ1(I,2))
      ELSE IF(K==KBM1) THEN
        TMP2 = (NHQ1(I,KBM2)*NHA1(I,KBM2)*DZ1(I,KBM1)+NHQ1(I,KBM1)*NHA1(I,KBM1)*DZ1(I,KBM2))/(DZ1(I,KBM2)+DZ1(I,KBM1)) - NHA1D(I)*NHQ1(I,KBM1)
      ELSE
        TMP2 = (NHQ1(I,K-1)*NHA1(I,K-1)*DZ1(I,K)+NHQ1(I,K)*NHA1(I,K)*DZ1(I,K-1))/(DZ1(I,K-1)+DZ1(I,K)) - (NHQ1(I,K)*NHA1(I,K)*DZ1(I,K+1)+NHQ1(I,K+1)*NHA1(I,K+1)*DZ1(I,K))/(DZ1(I,K)+DZ1(I,K+1))
      ENDIF
      NHQ2DX(I) = NHQ2DX(I) + 1.0E-3_SP*PQPPX(I)
      NHQDRX(I,K) = 1.0E-3_SP*(PQPPX(I)+TMP2*ART(I))

      IF(K==1) THEN
        TMP2 = -(NHQ1(I,1)*NHA2(I,1)*DZ1(I,2)+NHQ1(I,2)*NHA2(I,2)*DZ1(I,1))/(DZ1(I,1)+DZ1(I,2))
      ELSE IF(K==KBM1) THEN
        TMP2 = (NHQ1(I,KBM2)*NHA2(I,KBM2)*DZ1(I,KBM1)+NHQ1(I,KBM1)*NHA2(I,KBM1)*DZ1(I,KBM2))/(DZ1(I,KBM2)+DZ1(I,KBM1)) - NHA2D(I)*NHQ1(I,KBM1)
      ELSE
        TMP2 = (NHQ1(I,K-1)*NHA2(I,K-1)*DZ1(I,K)+NHQ1(I,K)*NHA2(I,K)*DZ1(I,K-1))/(DZ1(I,K-1)+DZ1(I,K)) - (NHQ1(I,K)*NHA2(I,K)*DZ1(I,K+1)+NHQ1(I,K+1)*NHA2(I,K+1)*DZ1(I,K))/(DZ1(I,K)+DZ1(I,K+1))
      ENDIF
      NHQ2DY(I) = NHQ2DY(I) + 1.0E-3_SP*PQPPY(I)
      NHQDRY(I,K) = 1.0E-3_SP*(PQPPY(I)+TMP2*ART(I))
    ENDDO

  ENDDO

# if defined (WET_DRY)
  DO I=1, N
    IF(ISWETC(I) /= 1) THEN
      NHQ2DX(I) = 0.0_SP
      NHQ2DY(I) = 0.0_SP
      NHQDRX(I,1:KBM1) = 0.0_SP
      NHQDRY(I,1:KBM1) = 0.0_SP
    ENDIF
  ENDDO
# endif

  DO I=1, M
#   if defined (WET_DRY)
    IF(ISWETN(I) ==1 ) THEN
#   endif
      DO K=2, KBM1
        NHQDRZ(I,K) = 1.0E-3_SP*( NHQ(I,K-1) - NHQ(I,K) )
      ENDDO
#   if defined (WET_DRY)
    ELSE
      NHQDRZ(I,2:KBM1) = 0.0_SP
    ENDIF
#   endif 
  ENDDO

  IF(IINT<=PROJ_SWITCH) THEN
    NHQDRX = 0.0_SP
    NHQDRY = 0.0_SP
    NHQDRZ = 0.0_SP
  ENDIF

  DO I=1, N
#   if defined (WET_DRY)
    IF(ISWETC(I)==1) THEN
#   endif    
      NHQ2DX(I) = NHQ2DX(I) - 1.0E-3_SP*ART(I)*NHQ1(I,KBM1)*NHA1D(I)
      NHQ2DY(I) = NHQ2DY(I) - 1.0E-3_SP*ART(I)*NHQ1(I,KBM1)*NHA2D(I)
#   if defined (WET_DRY)
    ENDIF
#   endif
  ENDDO

  CALL SET_NH_OBC

  ALLOCATE(XFL_EL(0:MT))        ; XFL_EL   = 0.0_SP
  ALLOCATE(YFL_EL(0:MT))        ; YFL_EL   = 0.0_SP
  ALLOCATE(XFL_D(0:MT,KB))      ; XFL_D    = 0.0_SP
  ALLOCATE(YFL_D(0:MT,KB))      ; YFL_D    = 0.0_SP
  DO I=1,NCV
    I1=NTRG(I)
    IA=NIEC(I,1)
    IB=NIEC(I,2)
    DIJ=(D(NV(I1,1))+D(NV(I1,2))+D(NV(I1,3)))/3.0_SP
    DO K=1, KB
      XFL_D(IA,K) = XFL_D(IA,K) + DIJ*Z1(I1,K)*DLTYE(I)
      XFL_D(IB,K) = XFL_D(IB,K) - DIJ*Z1(I1,K)*DLTYE(I)
      YFL_D(IA,K) = YFL_D(IA,K) - DIJ*Z1(I1,K)*DLTXE(I)
      YFL_D(IB,K) = YFL_D(IB,K) + DIJ*Z1(I1,K)*DLTXE(I)

      IF(ISONB(IA)>=1 .AND. ISONB(IB)>=1) THEN
      DO J1=1, 3
        IF(NV(I1,J1)/=IA .AND. NV(I1,J1)/=IB) THEN
          IF(NBE(I1,J1)==0) THEN

            DX = VX(IB) - VX(IA)
            DY = VY(IB) - VY(IA)
            XFL_D(IA,K) = XFL_D(IA,K) + 0.25_SP*(-DY)*(D(IA)*Z(IA,K)+D(IB)*Z(IB,K))
            XFL_D(IB,K) = XFL_D(IB,K) + 0.25_SP*(-DY)*(D(IA)*Z(IA,K)+D(IB)*Z(IB,K))
            YFL_D(IA,K) = XFL_D(IA,K) + 0.25_SP*DX*(D(IA)*Z(IA,K)+D(IB)*Z(IB,K))
            YFL_D(IB,K) = XFL_D(IB,K) + 0.25_SP*DX*(D(IA)*Z(IA,K)+D(IB)*Z(IB,K))

          ENDIF
        ENDIF
      ENDDO
      ENDIF
    ENDDO

    DIJ=(EL(NV(I1,1))+EL(NV(I1,2))+EL(NV(I1,3)))/3.0_SP
    XFL_EL(IA) = XFL_EL(IA) + DIJ*DLTYE(I)
    XFL_EL(IB) = XFL_EL(IB) - DIJ*DLTYE(I)
    YFL_EL(IA) = YFL_EL(IA) - DIJ*DLTXE(I)
    YFL_EL(IB) = YFL_EL(IB) + DIJ*DLTXE(I)

    IF(ISONB(IA)>=1 .AND. ISONB(IB)>=1) THEN
    DO J1=1, 3
      IF(NV(I1,J1)/=IA .AND. NV(I1,J1)/=IB) THEN
         IF(NBE(I1,J1)==0) THEN

            DX = VX(IB) - VX(IA)
            DY = VY(IB) - VY(IA)
            XFL_EL(IA) = XFL_EL(IA) + 0.25_SP*(-DY)*(EL(IA)+EL(IB))
            XFL_EL(IB) = XFL_EL(IB) + 0.25_SP*(-DY)*(EL(IA)+EL(IB))
            YFL_EL(IA) = XFL_EL(IA) + 0.25_SP*DX*(EL(IA)+EL(IB))
            YFL_EL(IB) = XFL_EL(IB) + 0.25_SP*DX*(EL(IA)+EL(IB))

         ENDIF
      ENDIF
    ENDDO
    ENDIF

  ENDDO

  AVG_U  = 0.0_SP
  AVG_V  = 0.0_SP
  DO K=1, KBM1
    DO I=1, M
      DO J=1, NTVE(I)
        AVG_U(I,K) = AVG_U(I,K) + UF(NBVE(I,J),K)*ART(NBVE(I,J))
        AVG_V(I,K) = AVG_V(I,K) + VF(NBVE(I,J),K)*ART(NBVE(I,J))
      ENDDO
      AVG_U(I,K) = AVG_U(I,K)/ART2(I)
      AVG_V(I,K) = AVG_V(I,K)/ART2(I)
    ENDDO
  ENDDO

  XFLUX = 0.0_SP
  DO I=1,NCV
    I1=NTRG(I)
    IA=NIEC(I,1)
    IB=NIEC(I,2)
    DIJ=(D(NV(I1,1))+D(NV(I1,2))+D(NV(I1,3)))/3.0_SP
    DO K=1,KBM1
      UIJ=UF(I1,K)
      VIJ=VF(I1,K)
      EXFLUX=DIJ*DZ1(I1,K)*(-UIJ*DLTYE(I)+VIJ*DLTXE(I))
      XFLUX(IA,K)=XFLUX(IA,K)-EXFLUX
      XFLUX(IB,K)=XFLUX(IB,K)+EXFLUX
    END DO
  ENDDO
! add external source 
  IF(NUMQBC >= 1) THEN
    IF(RIVER_INFLOW_LOCATION == 'node') THEN
      DO J=1,NUMQBC
        JJ=INODEQ(J)
        DO K=1,KBM1
          XFLUX(JJ,K)=XFLUX(JJ,K)-QDIS(J)*VQDIST(J,K)    !/DZ(JJ,K)
        END DO
      END DO
    ELSE IF(RIVER_INFLOW_LOCATION == 'edge') THEN
      DO J=1,NUMQBC
        J1=N_ICELLQ(J,1)
        J2=N_ICELLQ(J,2)
        DO K=1,KBM1
          XFLUX(J1,K)=XFLUX(J1,K)-QDIS(J)*RDISQ(J,1)*VQDIST(J,K)    !/DZ1(J1,K)
          XFLUX(J2,K)=XFLUX(J2,K)-QDIS(J)*RDISQ(J,2)*VQDIST(J,K)    !/DZ1(J2,K)
        END DO
      END DO
    END IF
  END IF
! add external source

  DO I=1, M
    W4ZT(I,KB) = -AVG_U(I,KBM1)*A1AVGD(I) - AVG_V(I,KBM1)*A2AVGD(I)
  ENDDO

  DO I=1, M

#   if defined (WET_DRY)
    IF(ISWETN(I)*ISWETNT(I) ==1) THEN
#   endif
    DO K= KBM1, 1, -1

      IF(K==1) THEN
        AVG_U_TMP = (-ZZ(I,2))*(-ZZ(I,3))/(ZZ(I,1)-ZZ(I,2))/(ZZ(I,1)-ZZ(I,3))*AVG_U(I,1)+(-ZZ(I,1))*(-ZZ(I,3))/(ZZ(I,2)-ZZ(I,1))/(ZZ(I,2)-ZZ(I,3))*AVG_U(I,2)+(-ZZ(I,1))*(-ZZ(I,2))/(ZZ(I,3)-ZZ(I,1))/(ZZ(I,3)-ZZ(I,2))*AVG_U(I,3)
        AVG_V_TMP = (-ZZ(I,2))*(-ZZ(I,3))/(ZZ(I,1)-ZZ(I,2))/(ZZ(I,1)-ZZ(I,3))*AVG_V(I,1)+(-ZZ(I,1))*(-ZZ(I,3))/(ZZ(I,2)-ZZ(I,1))/(ZZ(I,2)-ZZ(I,3))*AVG_V(I,2)+(-ZZ(I,1))*(-ZZ(I,2))/(ZZ(I,3)-ZZ(I,1))/(ZZ(I,3)-ZZ(I,2))*AVG_V(I,3)

        U2 = (AVG_U(I,1)*DZ(I,2)+AVG_U(I,2)*DZ(I,1))/(DZ(I,1)+DZ(I,2))
        V2 = (AVG_V(I,1)*DZ(I,2)+AVG_V(I,2)*DZ(I,1))/(DZ(I,1)+DZ(I,2))
        V_TMP1 = ((U2-AVG_U_TMP)*XFL_EL(I)+U2*XFL_D(I,2)-AVG_U_TMP*XFL_D(I,1)) + ((V2-AVG_V_TMP)*YFL_EL(I)+V2*YFL_D(I,2)-AVG_V_TMP*YFL_D(I,1))
      ELSE IF(K==KBM1) THEN
        AVG_U_TMP = (-1.0d0-ZZ(I,KBM1-1))*(-1.0d0-ZZ(I,KBM1-2))/(ZZ(I,KBM1)-ZZ(I,KBM1-1))/(ZZ(I,KBM1)-ZZ(I,KBM1-2))*AVG_U(I,KBM1)+(-1.0d0-ZZ(I,KBM1))*(-1.0d0-ZZ(I,KBM1-2))/(ZZ(I,KBM1-1)-ZZ(I,KBM1))/(ZZ(I,KBM1-1)-ZZ(I,KBM1-2))*AVG_U(I,KBM1-1)+(-1.0d0-ZZ(I,KBM1))*(-1.0d0-ZZ(I,KBM1-1))/(ZZ(I,KBM1-2)-ZZ(I,KBM1))/(ZZ(I,KBM1-2)-ZZ(I,KBM1-1))*AVG_U(I,KBM1-2)
        AVG_V_TMP = (-1.0d0-ZZ(I,KBM1-1))*(-1.0d0-ZZ(I,KBM1-2))/(ZZ(I,KBM1)-ZZ(I,KBM1-1))/(ZZ(I,KBM1)-ZZ(I,KBM1-2))*AVG_V(I,KBM1)+(-1.0d0-ZZ(I,KBM1))*(-1.0d0-ZZ(I,KBM1-2))/(ZZ(I,KBM1-1)-ZZ(I,KBM1))/(ZZ(I,KBM1-1)-ZZ(I,KBM1-2))*AVG_V(I,KBM1-1)+(-1.0d0-ZZ(I,KBM1))*(-1.0d0-ZZ(I,KBM1-1))/(ZZ(I,KBM1-2)-ZZ(I,KBM1))/(ZZ(I,KBM1-2)-ZZ(I,KBM1-1))*AVG_V(I,KBM1-2)

        U1 = (AVG_U(I,KBM1-1)*DZ(I,KBM1)+AVG_U(I,KBM1)*DZ(I,KBM1-1))/(DZ(I,KBM1-1)+DZ(I,KBM1))
        V1 = (AVG_V(I,KBM1-1)*DZ(I,KBM1)+AVG_V(I,KBM1)*DZ(I,KBM1-1))/(DZ(I,KBM1-1)+DZ(I,KBM1))
        V_TMP1 = ((AVG_U_TMP-U1)*XFL_EL(I)+AVG_U_TMP*XFL_D(I,KB)-U1*XFL_D(I,KBM1)) + ((AVG_V_TMP-V1)*YFL_EL(I)+AVG_V_TMP*YFL_D(I,KB)-V1*YFL_D(I,KBM1))
      ELSE
        U1 = (AVG_U(I,K-1)*DZ(I,K)+AVG_U(I,K)*DZ(I,K-1))/(DZ(I,K-1)+DZ(I,K))
        V1 = (AVG_V(I,K-1)*DZ(I,K)+AVG_V(I,K)*DZ(I,K-1))/(DZ(I,K-1)+DZ(I,K))
        U2 = (AVG_U(I,K)*DZ(I,K+1)+AVG_U(I,K+1)*DZ(I,K))/(DZ(I,K)+DZ(I,K+1))
        V2 = (AVG_V(I,K)*DZ(I,K+1)+AVG_V(I,K+1)*DZ(I,K))/(DZ(I,K)+DZ(I,K+1))

        V_TMP1 = ((U2-U1)*XFL_EL(I)+U2*XFL_D(I,K+1)-U1*XFL_D(I,K)) + ((V2-V1)*YFL_EL(I)+V2*YFL_D(I,K+1)-V1*YFL_D(I,K))
      ENDIF

      W4ZT(I,K) = W4ZT(I,K+1) - (XFLUX(I,K) + V_TMP1)/ART1(I)

    ENDDO
#   if defined (WET_DRY)
    ELSE
      W4ZT(I,1:KB) = 0.0_SP
    ENDIF
#   endif
  ENDDO

  CALL SET_W_OBC

# if defined(MULTIPROCESSOR)
  IF(PAR) CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,KB,MYID,NPROCS,W4ZT)
  IF(PAR) CALL AEXCHANGE(NC,MYID,NPROCS,W4ZT)
# endif

  DO K=1, KBM1
    DO I=1, NT
      WW(I,K) = (W4ZT(NV(I,1),K)+W4ZT(NV(I,2),K)+W4ZT(NV(I,3),K)+W4ZT(NV(I,1),K+1)+W4ZT(NV(I,2),K+1)+W4ZT(NV(I,3),K+1))/6.0_SP
    ENDDO
  ENDDO

  DEALLOCATE(XFL_EL,YFL_EL,XFL_D,YFL_D)
  DEALLOCATE(PQPPX,PQPPY,XFLUX)

  RETURN
  END SUBROUTINE UPDATE_QUVW

  SUBROUTINE UPDATE_ZETA

  USE LIMS
  USE CONTROL
  USE MOD_OBCS
  USE ALL_VARS
# if defined(MULTIPROCESSOR)
  USE MOD_PAR
# endif
  USE MOD_PETSc, ONLY: BL_EL, A_EL, B_EL, X_EL, XL_EL, ALO_2_PLO_NODE, N_VERTS, L2G_EL, G2L_EL, XVALS_EL, PLO_2_ALO_NODE, PETSc_SOLVER_EL
# if defined (SEMI_IMPLICIT)
  USE MOD_SEMI_IMPLICIT, ONLY: IFCETA
# endif
  IMPLICIT NONE

  INTEGER  I, J, J1, J2, JJ, IK, K, I1, I2, JN, IA, IB, TMP2, TMP3
  INTEGER  NNZ2, LN, LOC
  INTEGER  PETSc_POS, NODE, PROW1, PROW2, PCOL1, PCOL2, PCOL3
  REAL(SP) DIJ, DIJ1, UIJ, VIJ, EXFLUX, CC
  REAL(SP), ALLOCATABLE ::  UF_AVG(:), VF_AVG(:)
  REAL(SP), ALLOCATABLE ::  XFLUX(:)

  PetscReal :: STERM, VCOL1, VCOL2
  PetscInt  :: IERR

  ALLOCATE(UF_AVG(0:NT))   ; UF_AVG = 0.0_SP
  ALLOCATE(VF_AVG(0:NT))   ; VF_AVG = 0.0_SP
  ALLOCATE(XFLUX(0:MT))

  DO I=1,NT
    UF_AVG(I) = SUM(UF(I,1:KBM1)*DZ1(I,1:KBM1))
    VF_AVG(I) = SUM(VF(I,1:KBM1)*DZ1(I,1:KBM1))
  END DO

  CALL VecSet(BL_EL,ZERO,IERR);CHKERRQ(IERR)
  CALL VecSet(B_EL,ZERO,IERR);CHKERRQ(IERR)

  XFLUX   = 0.0_SP
  DO I=1,NCV
    I1=NTRG(I)
    IA=NIEC(I,1)
    IB=NIEC(I,2)
    DIJ=(H(NV(I1,1))+H(NV(I1,2))+H(NV(I1,3)))/3.0_SP
#   if defined (SEMI_IMPLICIT)
    DIJ1 = (DT(NV(I1,1))+DT(NV(I1,2))+DT(NV(I1,3)))/3.0_SP
#   endif

    UIJ=UF_AVG(I1)
    VIJ=VF_AVG(I1)
#   if !defined (SEMI_IMPLICIT)
    EXFLUX=DIJ*(-UIJ*DLTYE(I)+VIJ*DLTXE(I))
#   else
    EXFLUX=(1.0_SP-IFCETA)*DIJ1*(-UA(I1)*DLTYE(I)+VA(I1)*DLTXE(I)) + IFCETA*DIJ*(-UIJ*DLTYE(I)+VIJ*DLTXE(I))
#   endif
    XFLUX(IA)=XFLUX(IA)-EXFLUX
    XFLUX(IB)=XFLUX(IB)+EXFLUX
  ENDDO
! add external source
!  XFLUX = XFLUX+(QEVAP-QPREC)*ROFVROS*ART1
  XFLUX = XFLUX-(QEVAP+QPREC)*ROFVROS*ART1

  IF(GROUNDWATER_ON)THEN
    DO I=1,M
      XFLUX(I)=XFLUX(I)+BFWDIS(I)*ROFVROS*ART1(I)
    END DO
  END IF

  IF(NUMQBC >= 1) THEN
    IF(RIVER_INFLOW_LOCATION == 'node') THEN
      DO J=1,NUMQBC
        JJ=INODEQ(J)
        XFLUX(JJ)=XFLUX(JJ)-QDIS(J)
      END DO
    ELSE IF(RIVER_INFLOW_LOCATION == 'edge') THEN
      DO J=1,NUMQBC
        J1=N_ICELLQ(J,1)
        J2=N_ICELLQ(J,2)
        XFLUX(J1)=XFLUX(J1)-QDIS(J)*RDISQ(J,1)
        XFLUX(J2)=XFLUX(J2)-QDIS(J)*RDISQ(J,2)
      END DO
    END IF
  END IF
! add external source

  DO I=1, M
    PETSc_POS = ALO_2_PLO_NODE(I)
    IF (PETSc_POS > N_VERTS) CYCLE
    IF( ISONB(I)/=2 ) THEN
      STERM = ET(I) - DTI/ART1(I)*XFLUX(I)
    ELSE
      STERM = EL(I)
    ENDIF
    CALL VecSetValues(BL_EL,1,PETSc_POS-1,STERM,INSERT_VALUES,IERR);CHKERRQ(IERR)
  ENDDO

! radiation obc
  DO J = 1,IBCN(4)
    JN = OBC_LST(4,J)
    I1 = I_OBC_N(JN)
    I2 = NEXT_OBC(JN)
    PETSc_POS = ALO_2_PLO_NODE(I1)
    IF (PETSc_POS > N_VERTS) CYCLE

    CC = SQRT(GRAV_N(I1)*H(I1))*DTI/DLTN_OBC(JN)
    STERM = ET(I1)*(1.0_SP-DTI/10800.0_SP)/(1.0_SP+CC)

    CALL VecSetValues(BL_EL,1,PETSc_POS-1,STERM,INSERT_VALUES,IERR);CHKERRQ(IERR)
  END DO
!

# if defined (OLD_PETSC)
  CALL VecScatterBegin(BL_EL,B_EL,INSERT_VALUES,SCATTER_FORWARD,L2G_EL,IERR);CHKERRQ(IERR)
  CALL VecScatterEnd(BL_EL,B_EL,INSERT_VALUES,SCATTER_FORWARD,L2G_EL,IERR);CHKERRQ(IERR)
# else
  CALL VecScatterBegin(L2G_EL,BL_EL,B_EL,INSERT_VALUES,SCATTER_FORWARD,IERR);CHKERRQ(IERR)
  CALL VecScatterEnd(L2G_EL,BL_EL,B_EL,INSERT_VALUES,SCATTER_FORWARD,IERR);CHKERRQ(IERR)
# endif

  CALL MatZeroEntries(A_EL,IERR);CHKERRQ(IERR)   

  DO I=1, NCV_I
    I1 = NTRG(I)
    IA = NIEC(I,1)
    IB = NIEC(I,2)

    IF((ISONB(IA)+ISONB(IB))<4) THEN

#     if defined (WET_DRY)
      IF(ISWETCT(I1) == 1) THEN
#     endif

      DO J=1, 3
        IF(NV(I1,J)==IA) THEN
          TMP2 = NV(I1,J+1-INT((J+1)/4)*3)
          TMP3 = NV(I1,J+2-INT((J+2)/4)*3)
        ENDIF
      ENDDO

      PROW1 = ALO_2_PLO_NODE(IA)
      PROW2 = ALO_2_PLO_NODE(IB)

      PCOL1 = ALO_2_PLO_NODE(IA)
      PCOL2 = ALO_2_PLO_NODE(IB)
      PCOL3 = ALO_2_PLO_NODE(TMP3)

#     if !defined (SEMI_IMPLICIT)
      VCOL1 = DTI*(-VF_AVG(I1)*DLTXE(I)+UF_AVG(I1)*DLTYE(I))/ART1(IA)/3.0_SP
      VCOL2 = -DTI*(-VF_AVG(I1)*DLTXE(I)+UF_AVG(I1)*DLTYE(I))/ART1(IB)/3.0_SP
#     else
      VCOL1 = IFCETA*DTI*(-VF_AVG(I1)*DLTXE(I)+UF_AVG(I1)*DLTYE(I))/ART1(IA)/3.0_SP
      VCOL2 = -IFCETA*DTI*(-VF_AVG(I1)*DLTXE(I)+UF_AVG(I1)*DLTYE(I))/ART1(IB)/3.0_SP
#     endif
      CALL MatSetValuesLocal(A_EL,1,PROW1-1,1,PCOL1-1,VCOL1,ADD_VALUES,IERR);CHKERRQ(IERR)
      CALL MatSetValuesLocal(A_EL,1,PROW2-1,1,PCOL1-1,VCOL2,ADD_VALUES,IERR);CHKERRQ(IERR)

      CALL MatSetValuesLocal(A_EL,1,PROW1-1,1,PCOL2-1,VCOL1,ADD_VALUES,IERR);CHKERRQ(IERR)
      CALL MatSetValuesLocal(A_EL,1,PROW2-1,1,PCOL2-1,VCOL2,ADD_VALUES,IERR);CHKERRQ(IERR) 

      CALL MatSetValuesLocal(A_EL,1,PROW1-1,1,PCOL3-1,VCOL1,ADD_VALUES,IERR);CHKERRQ(IERR)
      CALL MatSetValuesLocal(A_EL,1,PROW2-1,1,PCOL3-1,VCOL2,ADD_VALUES,IERR);CHKERRQ(IERR)

#     if defined (WET_DRY)
      ENDIF
#     endif
    ENDIF
  ENDDO

  DO I=1, M

    PROW1  = ALO_2_PLO_NODE(I)
    IF(PROW1 > N_VERTS) CYCLE
    PCOL1  = ALO_2_PLO_NODE(I)
    
    VCOL1 = 1.0D0
    CALL MatSetValuesLocal(A_EL,1,PROW1-1,1,PCOL1-1,VCOL1,ADD_VALUES,IERR);CHKERRQ(IERR)    

  ENDDO

  CALL MatAssemblyBegin(A_EL,MAT_FINAL_ASSEMBLY,IERR);CHKERRQ(IERR)
  CALL MatAssemblyEnd(A_EL,MAT_FINAL_ASSEMBLY,IERR);CHKERRQ(IERR)
  CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

  DO I=1, IOBCN

    PROW1 = ALO_2_PLO_NODE(I_OBC_N(I))
    IF(PROW1>N_VERTS) CYCLE

    VCOL1 = 0.0D0
    NODE  = I_OBC_N(I)
    DO J=1, NTSN(NODE)-1

      PCOL1 = ALO_2_PLO_NODE(NBSN(NODE,J))
      CALL MatSetValuesLocal(A_EL,1,PROW1-1,1,PCOL1-1,VCOL1,INSERT_VALUES,IERR);CHKERRQ(IERR)

    ENDDO

    VCOL1 = 1.0D0
    PROW1 = ALO_2_PLO_NODE(NODE)
    CALL MatSetValuesLocal(A_EL,1,PROW1-1,1,PROW1-1,VCOL1,INSERT_VALUES,IERR);CHKERRQ(IERR)

  ENDDO

! radiation obc
  DO J = 1, IBCN(4)
    JN = OBC_LST(4,J)
    I1 = I_OBC_N(JN)
    I2 = NEXT_OBC(JN)
    PROW1 = ALO_2_PLO_NODE(I1)
    IF (PROW1 > N_VERTS) CYCLE

      PCOL1 = ALO_2_PLO_NODE(I2)
      CC = SQRT(GRAV_N(I1)*H(I1))*DTI/DLTN_OBC(JN)
      VCOL1 = -CC/(1.0_SP+CC)
      CALL MatSetValuesLocal(A_EL,1,PROW1-1,1,PCOL1-1,VCOL1,INSERT_VALUES,IERR);CHKERRQ(IERR)
  END DO
!

  CALL MatAssemblyBegin(A_EL,MAT_FINAL_ASSEMBLY,IERR);CHKERRQ(IERR)
  CALL MatAssemblyEnd(A_EL,MAT_FINAL_ASSEMBLY,IERR);CHKERRQ(IERR)

  CALL PETSc_SETICS_EL
  CALL PETSc_SOLVER_EL

# if defined (OLD_PETSC)
  CALL VecScatterBegin(X_EL,XL_EL,INSERT_VALUES,SCATTER_FORWARD,G2L_EL,IERR);CHKERRQ(IERR)
  CALL VecScatterEnd(X_EL,XL_EL,INSERT_VALUES,SCATTER_FORWARD,G2L_EL,IERR);CHKERRQ(IERR)
# else
  CALL VecScatterBegin(G2L_EL,X_EL,XL_EL,INSERT_VALUES,SCATTER_FORWARD,IERR);CHKERRQ(IERR)
  CALL VecScatterEnd(G2L_EL,X_EL,XL_EL,INSERT_VALUES,SCATTER_FORWARD,IERR);CHKERRQ(IERR)
# endif
  CALL VecGetArrayF90(XL_EL,XVALS_EL,IERR);CHKERRQ(IERR)

  EL = 0.0_SP
  DO I=1,N_VERTS
    IK = PLO_2_ALO_NODE(I)
    EL(IK) = XVALS_EL(I)
  ENDDO

# if defined (MULTIPROCESSOR)
    IF(PAR) CALL NODE_MATCH(0,NBN,BN_MLT,BN_LOC,BNC,MT,1,MYID,NPROCS,EL)
    IF(PAR) CALL AEXCHANGE(NC,MYID,NPROCS,EL)
# endif

# if defined (WET_DRY)
  IF(WETTING_DRYING_ON) THEN
    ELF = EL
    CALL WET_JUDGE
    IF(PAR) CALL AEXCHANGE(NC,MYID,NPROCS,ELF)
    EL = ELF
  ENDIF
# endif

  CALL N2E2D(EL,EL1)
  D   = H + EL
  D1  = H1 + EL1
# if !defined (SEMI_IMPLICIT)
  UA = UF_AVG
  VA = VF_AVG
  DTFA = D
# endif

  CALL VERTVL_EDGE
# if defined (WET_DRY)
  IF(WETTING_DRYING_ON) CALL WD_UPDATE(2)
# endif

  DEALLOCATE(UF_AVG,VF_AVG,XFLUX)

  RETURN
  END SUBROUTINE UPDATE_ZETA

  SUBROUTINE PETSc_SETICS
!============================================================================
! set petsc x vector with initial conditions from last solution note, petsc
! by default uses x=0 unless you specify KSPInitialGuessNonzero, which we
! do optionally using the use_last=.true. In this routine, if we have set
! use_last=.false., we just set x to zero and exit.  Setting it up would
! be a waste of time if Petsc KSP solver doesnt use it.
!============================================================================
  USE ALL_VARS, ONLY : MYID
  USE MOD_PETSc, ONLY : USE_LAST, X, PDEBUG, PUNIT, PSIZE, PLO_2_ALO, XLOCAL, L2G
  IMPLICIT NONE
  INTEGER :: I,IK(2),ICNT
  PetscReal :: QTERM
  PetscScalar :: ZERO    =  0.0D0
  PetscInt    :: IERR
  CHARACTER(LEN=20) :: SUBNAME = 'PETSc_SETICS'

  IF(.NOT.USE_LAST) THEN
    CALL VecSet(X,ZERO,IERR);CHKERRQ(IERR)
    RETURN
  ENDIF
   
  IF(PDEBUG) WRITE(PUNIT,*) 'STARTING: ',TRIM(SUBNAME)
  !-----------------------------------------------------------------------------
  ! set local solution vector using initial conditions from last Q 
  ! vertical coordinate is inner loop!  we want nearby diagonals to come
  ! from d^2/dz^2 discretization.
  ! index argument (icnt) is in C counting style (from 0)
  !-----------------------------------------------------------------------------
  ICNT = 0
  DO I=1,PSIZE
    !get local i,k ALO indices from PLO index 
    IK = PLO_2_ALO(I)

    !set value (should promote to double if FVCOM is running single precision)
    QTERM = QP(IK(1),IK(2))
    
    !insert the value into the local PLO-sized vector xlocal
    CALL VecSetValues(XLOCAL,1,I-1,QTERM,INSERT_VALUES,IERR);CHKERRQ(IERR)
  ENDDO

  !-----------------------------------------------------------------------------
  !load global vector using a scatter with the l2g index set  ==> x
  !-----------------------------------------------------------------------------
# if defined (OLD_PETSC)
  CALL VecScatterBegin(XLOCAL,X,INSERT_VALUES,SCATTER_FORWARD,L2G,IERR);CHKERRQ(IERR)
    !could put code here
  CALL VecScatterEnd(XLOCAL,X,INSERT_VALUES,SCATTER_FORWARD,L2G,IERR);CHKERRQ(IERR)
# else
  CALL VecScatterBegin(L2G,XLOCAL,X,INSERT_VALUES,SCATTER_FORWARD,IERR);CHKERRQ(IERR)
  !could put code here
  CALL VecScatterEnd(L2G,XLOCAL,X,INSERT_VALUES,SCATTER_FORWARD,IERR);CHKERRQ(IERR)
# endif

  IF(PDEBUG) WRITE(PUNIT,*) 'FINISHING: ',TRIM(SUBNAME)

  RETURN
  END SUBROUTINE PETSc_SETICS

  SUBROUTINE PETSc_SETICS_EL
!=========================================================================
!
!=========================================================================
  USE ALL_VARS,  ONLY: MYID, EL
  USE MOD_PETSc, ONLY: USE_LAST, XL_EL, N_VERTS, PLO_2_ALO_NODE, X_EL, L2G_EL

  IMPLICIT NONE
  INTEGER :: I, IK
  PetscReal :: QTERM
  PetscScalar :: ZERO    =  0.0D0
  PetscInt    :: IERR
  CHARACTER(LEN=20) :: SUBNAME = 'PETSc_SETICS'

  IF(.NOT.USE_LAST) THEN
    CALL VecSet(X_EL,ZERO,IERR);CHKERRQ(IERR)
    RETURN
  ENDIF

  DO I=1,N_VERTS
    IK = PLO_2_ALO_NODE(I)
    QTERM = EL(IK)
    CALL VecSetValues(XL_EL,1,I-1,QTERM,INSERT_VALUES,IERR);CHKERRQ(IERR)
  ENDDO

# if defined (OLD_PETSC)
  CALL VecScatterBegin(XL_EL,X_EL,INSERT_VALUES,SCATTER_FORWARD,L2G_EL,IERR);CHKERRQ(IERR)
  CALL VecScatterEnd(XL_EL,X_EL,INSERT_VALUES,SCATTER_FORWARD,L2G_EL,IERR);CHKERRQ(IERR)
# else
  CALL VecScatterBegin(L2G_EL,XL_EL,X_EL,INSERT_VALUES,SCATTER_FORWARD,IERR);CHKERRQ(IERR)
  CALL VecScatterEnd(L2G_EL,XL_EL,X_EL,INSERT_VALUES,SCATTER_FORWARD,IERR);CHKERRQ(IERR)
# endif

  RETURN
  END SUBROUTINE PETSc_SETICS_EL

  SUBROUTINE PETSc_GETQ
  !============================================================================
  !                         PETSc_GETQ
  ! scatter solution (x) back to petsc local (xlocal)
  ! map to the non-hydrostatic pressure (q)
  ! ensure boundary nodes have correct values
  !
  ! basic function:  set FVCOM field Q
  !
  ! dependencies:  x,xlocal,g2l,ierr,xvals,plo_2_alo,Psize from this module
  !                par,q,mt,kb,myid,nprocs from all_vars
  !                node_set,exchange,+various vars from mod_par
  !============================================================================
  USE ALL_VARS, ONLY : PAR,MT,KBM1,MYID,NPROCS,NT
# if defined(MULTIPROCESSOR)
  USE MOD_PAR,  ONLY : NBN,BN_MLT,BN_LOC,BNC,NC,AEXCHANGE,NODE_MATCH
# endif
  USE MOD_PETSc, ONLY: PDEBUG, PUNIT, X, XLOCAL, G2L, XVALS, PSIZE, PLO_2_ALO

  IMPLICIT NONE
  INTEGER :: I,IK(2)
  CHARACTER(LEN=20) :: SUBNAME = 'PETSc_GETQ'
  PetscInt  :: IERR

  IF(PDEBUG) WRITE(PUNIT,*) 'STARTING: ',TRIM(SUBNAME)
  !--------------------------------------------------------------------------
  ! scatter the solution vector back to local array and set back in q solution
  ! in theory, xvals is just a pointer and will not require additional memory.
  !--------------------------------------------------------------------------
# if defined (OLD_PETSC)
  CALL VecScatterBegin(X,XLOCAL,INSERT_VALUES,SCATTER_FORWARD,G2L,IERR);CHKERRQ(IERR)
  CALL VecScatterEnd(X,XLOCAL,INSERT_VALUES,SCATTER_FORWARD,G2L,IERR);CHKERRQ(IERR)
# else
  CALL VecScatterBegin(G2L,X,XLOCAL,INSERT_VALUES,SCATTER_FORWARD,IERR);CHKERRQ(IERR)
  CALL VecScatterEnd(G2L,X,XLOCAL,INSERT_VALUES,SCATTER_FORWARD,IERR);CHKERRQ(IERR)
# endif
  CALL VecGetArrayF90(XLOCAL,XVALS,IERR);CHKERRQ(IERR)

  !--------------------------------------------------------------------------
  !zero out q so that we can set properly on boundary nodes
  !see below for explanation
  !--------------------------------------------------------------------------
  QP = 0.0_SP

  !--------------------------------------------------------------------------
  !put the solution from a given petsc partition into local 3D q array
  !p2f_loc gives the ALO from the PLO
  !we need this transform because partitions are not the same as Petsc
  !does not allow for overlapping partitions
  !should we index xvals from 0 or starting from 1?
  !we will first code , assuming it starts from 1
  !--------------------------------------------------------------------------
  DO I=1,PSIZE
    IK = PLO_2_ALO(I)
    QP(IK(1),IK(2)) = XVALS(I)
  ENDDO

  !--------------------------------------------------------------------------
  ! set on boundary nodes issue is that petsc cannot have overlapping domains so
  ! interprocessor so an interprocessor boundary nodes from FVCOM partition will
  ! lie in only one Petsc partition and thus solution on that node from Petsc 
  ! solve is only valid in that one partition we must set that value in other
  ! partitions containing the boundary node we do this by summing the value at
  ! each boundary node.  In procs not containing the BN in the petsc partition,
  ! the initial value will be zero thus the final sum will be equal to the single
  ! correct value in all FVCOM partitions containing a given boundary node.
  ! finally, exchange on q to make sure halo nodes setup properly
  ! note: in single processor case, we do not have to do this, petsc
  ! local partition is same as application local partition
  !--------------------------------------------------------------------------
# if defined (MULTIPROCESSOR)
    IF(PAR) CALL NODE_MATCH(0,NBN,BN_MLT,BN_LOC,BNC,MT,KBM1,MYID,NPROCS,QP)
    IF(PAR) CALL AEXCHANGE(NC,MYID,NPROCS,QP)
# endif

# if defined (WET_DRY)
  DO I=1, NT
   IF(ISWETN(I)/=1) THEN
     QP(I,1:KBM1) = 0.0_SP
   ENDIF
  ENDDO
# endif

  IF(PDEBUG) WRITE(PUNIT,*) 'FINISHING: ',TRIM(SUBNAME)

  RETURN
  END SUBROUTINE PETSc_GETQ

  SUBROUTINE TRIDAG(a,b,c,r,u,n)
    INTEGER n,NMAX
    REAL(sp) a(n),b(n),c(n),r(n),u(n)
    PARAMETER (NMAX=500)
    INTEGER j
    REAL(sp) bet,gam(NMAX)
    if(b(1).eq.0.)pause 'tridag: rewrite equations'
    bet=b(1)
    u(1)=r(1)/bet
    do 11 j=2,n
      gam(j)=c(j-1)/bet
      bet=b(j)-a(j)*gam(j)
      if(bet.eq.0.)pause 'tridag failed'
      u(j)=(r(j)-a(j)*u(j-1))/bet
11  continue
    do 12 j=n-1,1,-1
      u(j)=u(j)-gam(j+1)*u(j+1)
12  continue
  return
  END SUBROUTINE TRIDAG

# endif
END MODULE NON_HYDRO

