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

MODULE MOD_SEMI_IMPLICIT
# if defined (SEMI_IMPLICIT)

  USE MOD_PREC, ONLY : SP, DP
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: NAME_LIST_READ_SEMI
  PUBLIC :: ALLOC_VARS_SEMI
  PUBLIC :: SEMI_IMPLICIT_EL
  PUBLIC :: UPDATES_SEMI
  PUBLIC :: UV2D_SBD
  PUBLIC :: UV3D_SBD
  PUBLIC :: NAME_LIST_INITIALIZE_SEMI
  PUBLIC :: NAME_LIST_PRINT_SEMI

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

  INTEGER , PUBLIC  :: STAGE, KSTAGE_UV, KSTAGE_TE, KSTAGE_TS
  REAL(SP), PUBLIC  :: IFCETA, BEDF
  CHARACTER(LEN=80), PUBLIC :: MSTG
  REAL(SP), PUBLIC, ALLOCATABLE :: RK_UV(:), RK_TE(:), RK_TS(:)
  REAL(SP), PUBLIC, ALLOCATABLE :: XFLUX3(:,:), YFLUX3(:,:)
  REAL(SP), PUBLIC, ALLOCATABLE :: UC(:,:), VC(:,:), UAC(:), VAC(:)
  REAL(SP), PUBLIC, ALLOCATABLE :: Q2C(:,:), Q2LC(:,:), TE_TMP(:,:)
  REAL(SP), PUBLIC, ALLOCATABLE :: TSC(:,:)
  REAL(SP), PUBLIC, ALLOCATABLE :: U_N(:,:), V_N(:,:), UA_N(:), VA_N(:) 
  
  NAMELIST /NML_SEMI/        &
          & IFCETA,          &
          & BEDF,            &
          & KSTAGE_UV,       &
          & KSTAGE_TE,       &
          & KSTAGE_TS,       &
          & MSTG

  PetscViewer :: viewer

  CONTAINS

!!  SUBROUTINE SET_IMPLICIT_PARAM
!!
!!  USE CONTROL
!!!  USE MOD_INP
!!  IMPLICIT NONE
!!  INTEGER  :: ISCAN
!!  CHARACTER(LEN=120) :: FNAME

!!  FNAME = "./"//trim(casename)//"_run.dat"

!!!----------------------------------------------------------------------------|
!!!   implicit factor ceta 
!!!----------------------------------------------------------------------------|
!!  ISCAN = SCAN_FILE(TRIM(FNAME),"IFCETA",FSCAL = IFCETA)
!!  IF(ISCAN /= 0)THEN
!!    WRITE(IPT,*)'ERROR READING IFCETA: ',ISCAN
!!    CALL PSTOP
!!  END IF

!!!----------------------------------------------------------------------------|
!!!   implicit factor ceta
!!!----------------------------------------------------------------------------|
!!  ISCAN = SCAN_FILE(TRIM(FNAME),"BEDF",FSCAL = BEDF)
!!  IF(ISCAN /= 0)THEN
!!    WRITE(IPT,*)'ERROR READING BEDF: ',ISCAN
!!    CALL PSTOP
!!  END IF

!!!==============================================================================|
!!!            SCREEN REPORT                                                     !
!!!==============================================================================|
!!  IF(MSR) THEN
!!    WRITE(IPT,*) '!                                                    !'
!!    WRITE(IPT,*) '!-------------- SPECIFY IMPLICIT FACTOR -------------!'
!!    WRITE(IPT,*) '!                                                    !'
!!    WRITE(IPT,*) '!  # IMPLICIT FACTOR SELECTED                        :',IFCETA
!!    WRITE(IPT,*) '!  # SOLID BOUNDARY ENERGY DAMPING FACTOR SELECTED   :',BEDF
!!    WRITE(IPT,*) '!                                                    !'
!!  ENDIF

!!  RETURN
!!  END SUBROUTINE SET_IMPLICIT_PARAM

  SUBROUTINE NAME_LIST_READ_SEMI

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

  ! Read SEMI Settings
  READ(UNIT=NMLUNIT, NML=NML_SEMI,IOSTAT=ios)
  if(ios .NE. 0 ) then
    if(DBG_SET(dbg_log)) write(UNIT=IPT,NML=NML_SEMI)
    Call Fatal_Error("Can Not Read NameList NML_SEMI from file: "//trim(FNAME))
  end if

  REWIND(NMLUNIT)

  if(DBG_SET(dbg_scl)) &
       & write(IPT,*) "Read_Name_List: NML_SEMI"

  if(DBG_SET(dbg_scl)) &
       & write(UNIT=IPT,NML=NML_SEMI)

  CLOSE(NMLUNIT)

  RETURN
  END SUBROUTINE NAME_LIST_READ_SEMI

  SUBROUTINE ALLOC_VARS_SEMI

  USE MOD_PREC, ONLY : SP, DP
  USE ALL_VARS, ONLY : NT,MT,KBM1,KB
# if defined (SPHERICAL)  
  USE MOD_NORTHPOLE, ONLY : XFLUX3_NP, YFLUX3_NP, UBETA_NP, VBETA_NP
# endif
  IMPLICIT NONE
  INTEGER K

  ALLOCATE(XFLUX3(0:NT,KBM1))         ; XFLUX3    = 0.0_SP
  ALLOCATE(YFLUX3(0:NT,KBM1))         ; YFLUX3    = 0.0_SP
# if defined (SPHERICAL) 
  ALLOCATE(XFLUX3_NP(0:NT,KBM1))      ; XFLUX3_NP = 0.0_SP
  ALLOCATE(YFLUX3_NP(0:NT,KBM1))      ; YFLUX3_NP = 0.0_SP
  ALLOCATE(UBETA_NP(0:NT,KB))         ; UBETA_NP  = 0.0_SP   
  ALLOCATE(VBETA_NP(0:NT,KB))         ; VBETA_NP  = 0.0_SP 
# endif
  ALLOCATE(UA_N(0:NT))                ; UA_N      = 0.0_SP
  ALLOCATE(VA_N(0:NT))                ; VA_N      = 0.0_SP
  ALLOCATE(U_N(0:NT,KBM1))            ; U_N       = 0.0_SP
  ALLOCATE(V_N(0:NT,KBM1))            ; V_N       = 0.0_SP

  ALLOCATE(UC(0:NT,KB))               ; UC        = 0.0_SP
  ALLOCATE(VC(0:NT,KB))               ; VC        = 0.0_SP
  ALLOCATE(UAC(0:NT))                 ; UAC       = 0.0_SP
  ALLOCATE(VAC(0:NT))                 ; VAC       = 0.0_SP

  ALLOCATE(Q2C(0:MT,KB))              ; Q2C       = 0.0_SP
  ALLOCATE(Q2LC(0:MT,KB))             ; Q2LC      = 0.0_SP
  ALLOCATE(TE_TMP(0:MT,KB))           ; TE_TMP    = 0.0_SP

  ALLOCATE(TSC(0:MT,KB))              ; TSC       = 0.0_SP

  ALLOCATE(RK_UV(KSTAGE_UV))          ; RK_UV     = 1.0_SP
  DO K=1, KSTAGE_UV
    RK_UV(K) = 1.0_SP/REAL(KSTAGE_UV+1-K)
  ENDDO

  ALLOCATE(RK_TE(KSTAGE_TE))          ; RK_TE     = 1.0_SP
  DO K=1, KSTAGE_TE
    RK_TE(K) = 1.0_SP/REAL(KSTAGE_TE+1-K)
  ENDDO 

  ALLOCATE(RK_TS(KSTAGE_TS))          ; RK_TS     = 1.0_SP
  DO K=1, KSTAGE_TS
    RK_TS(K) = 1.0_SP/REAL(KSTAGE_TS+1-K)
  ENDDO

  RETURN
  END SUBROUTINE ALLOC_VARS_SEMI

  SUBROUTINE SEMI_IMPLICIT_EL

  USE MOD_PREC, ONLY : SP, DP
  USE LIMS
  USE CONTROL
  USE MOD_OBCS
  USE ALL_VARS
  USE BCS
  USE MOD_PETSC, ONLY: A_EL,X_EL,B_EL,BL_EL,XL_EL,N_VERTS,L2G_EL,G2L_EL,XVALS_EL,ALO_2_PLO_NODE,PLO_2_ALO_NODE,PETSc_SOLVER_EL
# if defined(MULTIPROCESSOR)
  USE MOD_PAR
# endif
# if defined (SPHERICAL)   
  USE MOD_SPHERICAL
  USE MOD_NORTHPOLE,   ONLY :  NPE, NPEDGE_LST, NPCV,  NCEDGE_LST, NODE_NORTHPOLE, NP, NP_LST, CELL_NORTHAREA, XFLUX3_NP, YFLUX3_NP
# endif   
# if defined (WET_DRY)
  USE MOD_WD, ONLY : ISWETC, ISWETCT, WET_JUDGE
# endif

  IMPLICIT NONE 
  INTEGER II, I, J, JJ, J1, J2, IK, K, I1, I2, JN, IA, IB, TMP2, TMP3
  INTEGER  PETSc_POS, NODE, PROW1, PROW2, PCOL1, PCOL2, PCOL3
  REAL(SP) DIJ, UIJ, VIJ, EXFLUX, COEF, CC
  REAL(SP) Y3Y2,Y1Y3,Y2Y1,X3X2,X1X3,X2X1
  REAL(SP) :: XFLUX(0:MT)
  REAL(SP), ALLOCATABLE, DIMENSION(:) :: XFLUX_TMP, YFLUX_TMP, XFLUX_TMP_NP, YFLUX_TMP_NP

  REAL(SP) :: UA_TMP, VA_TMP, WUSURF_TMP, WVSURF_TMP, WUBOT_TMP, WVBOT_TMP
  REAL(SP) :: U_TMP, V_TMP, UF_TMP, VF_TMP
  REAL(SP) :: UIJ_TMP,VIJ_TMP,EXFLUX_TMP,FLUX_TMP, ELIJ, ELIJ2, ELIJ3
  REAL(SP) :: DLTXE_TMP,DLTYE_TMP,XTMP, XTMP1, DLTXC_TMP, DLTYC_TMP
  REAL(SP) :: VX1_TMP,VX2_TMP,VY1_TMP,VY2_TMP,VX3_TMP,VY3_TMP
  REAL(SP) :: XXX1, XXX2, XXX3, YYY1, YYY2, YYY3

  PetscReal :: STERM, VCOL1,VCOL2,VCOL3
  PetscInt  :: IERR
  PetscScalar :: ZERO1   =  0.0D0  

# if defined (TWO_D_MODEL)
  ALLOCATE(XFLUX_TMP(0:NT));   XFLUX_TMP = 0.0_SP
  ALLOCATE(YFLUX_TMP(0:NT));   YFLUX_TMP = 0.0_SP
  DO I=1, N
#   if defined (WET_DRY)
    IF(ISWETCT(I) == 1) THEN
#   endif
    XFLUX_TMP(I) = UA(I)*DT1(I) - DTI*IFCETA*ADVUA(I)/ART(I) + DTI*IFCETA*( -WUSURF(I) + WUBOT(I) ) - IFCETA*BEDF*UA_N(I)*DT1(I)
    YFLUX_TMP(I) = VA(I)*DT1(I) - DTI*IFCETA*ADVVA(I)/ART(I) + DTI*IFCETA*( -WVSURF(I) + WVBOT(I) ) - IFCETA*BEDF*VA_N(I)*DT1(I)
#   if defined (WET_DRY)
    ENDIF
#   endif
  ENDDO
# if defined (MULTIPROCESSOR)
  IF(PAR) CALL AEXCHANGE(EC,MYID,NPROCS,XFLUX_TMP,YFLUX_TMP)
# endif

# else
  ALLOCATE(XFLUX_TMP(0:NT));   XFLUX_TMP = 0.0_SP
  ALLOCATE(YFLUX_TMP(0:NT));   YFLUX_TMP = 0.0_SP
  DO I=1, N
#   if defined (WET_DRY)
    IF(ISWETCT(I) == 1) THEN
#   endif
    DO K=1, KBM1
      XFLUX_TMP(I) = XFLUX_TMP(I) - DTI*IFCETA*XFLUX3(I,K)/ART(I) - IFCETA*BEDF*U_N(I,K)*DT1(I)*DZ1(I,K)
      YFLUX_TMP(I) = YFLUX_TMP(I) - DTI*IFCETA*YFLUX3(I,K)/ART(I) - IFCETA*BEDF*V_N(I,K)*DT1(I)*DZ1(I,K)
    ENDDO
    XFLUX_TMP(I) = XFLUX_TMP(I) + UA(I)*DT1(I) + DTI*IFCETA*( -WUSURF(I) + WUBOT(I) )
    YFLUX_TMP(I) = YFLUX_TMP(I) + VA(I)*DT1(I) + DTI*IFCETA*( -WVSURF(I) + WVBOT(I) )
#   if defined (WET_DRY)
    ENDIF
#   endif
  ENDDO
# if defined (MULTIPROCESSOR)
  IF(PAR) CALL AEXCHANGE(EC,MYID,NPROCS,XFLUX_TMP,YFLUX_TMP)
# endif

!# if defined (SPHERICAL) && (NORTHPOLE)
# if defined (SPHERICAL)
  ALLOCATE(XFLUX_TMP_NP(0:NT));   XFLUX_TMP_NP = 0.0_SP
  ALLOCATE(YFLUX_TMP_NP(0:NT));   YFLUX_TMP_NP = 0.0_SP
  DO II=1, NP
    I = NP_LST(II)
    IF(CELL_NORTHAREA(I) ==1 ) THEN
      DO K=1, KBM1
        XFLUX_TMP_NP(I) = XFLUX_TMP_NP(I) - DTI*IFCETA*XFLUX3_NP(I,K)/ART(I)
        YFLUX_TMP_NP(I) = YFLUX_TMP_NP(I) - DTI*IFCETA*YFLUX3_NP(I,K)/ART(I)
      ENDDO
      UA_TMP = -VA(I)*COS(XC(I)*DEG2RAD)-UA(I)*SIN(XC(I)*DEG2RAD)
      VA_TMP = -VA(I)*SIN(XC(I)*DEG2RAD)+UA(I)*COS(XC(I)*DEG2RAD)
      WUSURF_TMP = -WVSURF(I)*COS(XC(I)*DEG2RAD)-WUSURF(I)*SIN(XC(I)*DEG2RAD)
      WVSURF_TMP = -WVSURF(I)*SIN(XC(I)*DEG2RAD)+WUSURF(I)*COS(XC(I)*DEG2RAD)
      WUBOT_TMP  = -WVBOT(I)*COS(XC(I)*DEG2RAD)-WUBOT(I)*SIN(XC(I)*DEG2RAD)
      WVBOT_TMP  = -WVBOT(I)*SIN(XC(I)*DEG2RAD)+WUBOT(I)*COS(XC(I)*DEG2RAD)
      XFLUX_TMP_NP(I) = XFLUX_TMP_NP(I) + UA_TMP*DT1(I) + DTI*IFCETA*( -WUSURF_TMP + WUBOT_TMP )
      YFLUX_TMP_NP(I) = YFLUX_TMP_NP(I) + VA_TMP*DT1(I) + DTI*IFCETA*( -WVSURF_TMP + WVBOT_TMP )
    ENDIF
  ENDDO
# if defined (MULTIPROCESSOR)
  IF(PAR) CALL AEXCHANGE(EC,MYID,NPROCS,XFLUX_TMP_NP,YFLUX_TMP_NP)
# endif
# endif
!JQI# if defined (MULTIPROCESSOR)
!JQI  IF(PAR) CALL AEXCHANGE(EC,MYID,NPROCS,XFLUX_TMP_NP,YFLUX_TMP_NP)
!JQI# endif

# endif

  CALL VecSet(BL_EL,ZERO1,IERR);CHKERRQ(IERR)
  CALL VecSet(B_EL,ZERO1,IERR);CHKERRQ(IERR)

  XFLUX   = 0.0_SP
  DO I=1,NCV
    I1=NTRG(I)
    IA=NIEC(I,1)
    IB=NIEC(I,2)

    UIJ=XFLUX_TMP(I1)
    VIJ=YFLUX_TMP(I1)
    EXFLUX= -UIJ*DLTYE(I)+VIJ*DLTXE(I)
    XFLUX(IA)=XFLUX(IA)-EXFLUX
    XFLUX(IB)=XFLUX(IB)+EXFLUX
  ENDDO

!# if defined (SPHERICAL) && (NORTHPOLE)
# if defined (SPHERICAL) 
  IF(NODE_NORTHPOLE/=0) XFLUX(NODE_NORTHPOLE) = 0.0_SP   

  DO II=1,NPCV
    I = NCEDGE_LST(II)
    I1  = NTRG(I)
    IA  = NIEC(I,1)
    IB  = NIEC(I,2)

    UIJ = XFLUX_TMP_NP(I1)
    VIJ = YFLUX_TMP_NP(I1)

    IF(IA == NODE_NORTHPOLE .OR. IB == NODE_NORTHPOLE)THEN
      VX1_TMP = REARTH * COS(YIJE(I,1)*DEG2RAD) * COS(XIJE(I,1)*DEG2RAD)!&
!                * 2._SP /(1._SP+sin(YIJE(I,1)*DEG2RAD))
      VY1_TMP = REARTH * COS(YIJE(I,1)*DEG2RAD) * SIN(XIJE(I,1)*DEG2RAD)!&
!                * 2._SP /(1._SP+sin(YIJE(I,1)*DEG2RAD))

      VX2_TMP = REARTH * COS(YIJE(I,2)*DEG2RAD) * COS(XIJE(I,2)*DEG2RAD)!&
!                * 2._SP /(1._SP+sin(YIJE(I,2)*DEG2RAD))
      VY2_TMP = REARTH * COS(YIJE(I,2)*DEG2RAD) * SIN(XIJE(I,2)*DEG2RAD)!&
!                * 2._SP /(1._SP+sin(YIJE(I,2)*DEG2RAD))

      DLTXE_TMP = VX2_TMP-VX1_TMP
      DLTYE_TMP = VY2_TMP-VY1_TMP

      EXFLUX_TMP = -UIJ*DLTYE_TMP+VIJ*DLTXE_TMP
    END IF

    IF(IA == NODE_NORTHPOLE) THEN
      XFLUX(IA) = XFLUX(IA)-EXFLUX_TMP
    ELSE IF(IB == NODE_NORTHPOLE)THEN
      XFLUX(IB) = XFLUX(IB)+EXFLUX_TMP
    END IF

  END DO 
  
# if !defined (TWO_D_MODEL)
  DEALLOCATE(XFLUX_TMP_NP,YFLUX_TMP_NP)
# endif
# endif


IF (PRECIPITATION_ON) THEN       !!yding
!  XFLUX = XFLUX+(QEVAP-QPREC)*ROFVROS*ART1
  XFLUX = XFLUX-(QEVAP+QPREC)*ROFVROS*ART1
END IF

#   if defined (ICE_EMBEDDING)
    DO I=1,M
      XFLUX(I)=XFLUX(I) + QFLICE(I)*(917.0_SP/((RHO1(I,1)+1.0_SP)*1.0E3_SP))*ART1(I)  !!yding
    END DO
#   endif






  IF(GROUNDWATER_ON)THEN
    DO I=1,M
!      XFLUX(I)=XFLUX(I)+BFWDIS(I)*ROFVROS*ART1(I)
      XFLUX(I)=XFLUX(I)-BFWDIS(I)
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

  DO I=1, M
    PETSc_POS = ALO_2_PLO_NODE(I)
    IF (PETSc_POS > N_VERTS) CYCLE
!!$    IF( ISONB(I)/=2 ) THEN
!!$    IF( ISONB(I) < 2 ) THEN
    IF( ISONB(I) /= 2 .AND. ISONB_3(I) /= 3) THEN
      STERM = ET(I) - DTI/ART1(I)*XFLUX(I)
    ELSE
      STERM = ELF(I)
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

  COEF  = DTI**2*IFCETA**2
  DO I=1, NCV_I
    I1 = NTRG(I)
    IA = NIEC(I,1)
    IB = NIEC(I,2)

!!$    IF((ISONB(IA)+ISONB(IB))<4) THEN

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

#     if defined (SPHERICAL)

!#     if defined (NORTHPOLE)
      IF (CELL_NORTHAREA(I1)==1) THEN

        VX1_TMP = REARTH * COS(VY(IA)*DEG2RAD) * COS(VX(IA)*DEG2RAD)&
                * 2._SP /(1._SP+sin(VY(IA)*DEG2RAD))
        VY1_TMP = REARTH * COS(VY(IA)*DEG2RAD) * SIN(VX(IA)*DEG2RAD)&
                * 2._SP /(1._SP+sin(VY(IA)*DEG2RAD))
    
        VX2_TMP = REARTH * COS(VY(IB)*DEG2RAD) * COS(VX(IB)*DEG2RAD)&
                * 2._SP /(1._SP+sin(VY(IB)*DEG2RAD))
        VY2_TMP = REARTH * COS(VY(IB)*DEG2RAD) * SIN(VX(IB)*DEG2RAD)&
                * 2._SP /(1._SP+sin(VY(IB)*DEG2RAD))
    
        VX3_TMP = REARTH * COS(VY(TMP3)*DEG2RAD) * COS(VX(TMP3)*DEG2RAD)&
                * 2._SP /(1._SP+sin(VY(TMP3)*DEG2RAD))
        VY3_TMP = REARTH * COS(VY(TMP3)*DEG2RAD) * SIN(VX(TMP3)*DEG2RAD)&
                * 2._SP /(1._SP+sin(VY(TMP3)*DEG2RAD))

        X3X2  = VX3_TMP-VX2_TMP
        X1X3  = VX1_TMP-VX3_TMP
        X2X1  = VX2_TMP-VX1_TMP
        Y3Y2  = VY3_TMP-VY2_TMP
        Y1Y3  = VY1_TMP-VY3_TMP
        Y2Y1  = VY2_TMP-VY1_TMP
      
        XXX1  =( X3X2*COS(XC(I1)*DEG2RAD)+Y3Y2*SIN(XC(I1)*DEG2RAD) )/(2._SP/(1._SP+sin(YC(I1)*DEG2RAD)))
        XXX2  =( X1X3*COS(XC(I1)*DEG2RAD)+Y1Y3*SIN(XC(I1)*DEG2RAD) )/(2._SP /(1._SP+sin(YC(I1)*DEG2RAD)))
        XXX3  =( X2X1*COS(XC(I1)*DEG2RAD)+Y2Y1*SIN(XC(I1)*DEG2RAD) )/(2._SP /(1._SP+sin(YC(I1)*DEG2RAD)))
        YYY1  =(-X3X2*SIN(XC(I1)*DEG2RAD)+Y3Y2*COS(XC(I1)*DEG2RAD) )/(2._SP /(1._SP+sin(YC(I1)*DEG2RAD)))
        YYY2  =(-X1X3*SIN(XC(I1)*DEG2RAD)+Y1Y3*COS(XC(I1)*DEG2RAD) )/(2._SP /(1._SP+sin(YC(I1)*DEG2RAD)))
        YYY3  =(-X2X1*SIN(XC(I1)*DEG2RAD)+Y2Y1*COS(XC(I1)*DEG2RAD) )/(2._SP /(1._SP+sin(YC(I1)*DEG2RAD)))

        VCOL1 = 0.5_SP*GRAV_E(I1)*COEF*(XXX1*DLTYE(I)-YYY1*DLTXE(I))*DT1(I1)/ART(I1)/ART1(IA)
        VCOL2 = 0.5_SP*GRAV_E(I1)*COEF*(XXX2*DLTYE(I)-YYY2*DLTXE(I))*DT1(I1)/ART(I1)/ART1(IA)
        VCOL3 = 0.5_SP*GRAV_E(I1)*COEF*(XXX3*DLTYE(I)-YYY3*DLTXE(I))*DT1(I1)/ART(I1)/ART1(IA)

      ELSE
!#     endif
      XTMP  = VX(TMP3)-VX(TMP2)
      IF(XTMP >  180.0_SP)THEN
          X3X2 = -360.0_SP*TPI+XTMP*TPI
        ELSE IF(XTMP < -180.0_SP)THEN
          X3X2 =  360.0_SP*TPI+XTMP*TPI
        ELSE
          X3X2 = XTMP*TPI
      END IF
      XTMP  = VX(IA)  -VX(TMP3)
      IF(XTMP >  180.0_SP)THEN
          X1X3 = -360.0_SP*TPI+XTMP*TPI
        ELSE IF(XTMP < -180.0_SP)THEN
          X1X3 =  360.0_SP*TPI+XTMP*TPI
        ELSE
          X1X3 = XTMP*TPI
      END IF
      XTMP  = VX(TMP2)-VX(IA)
      IF(XTMP >  180.0_SP)THEN
          X2X1 = -360.0_SP*TPI+XTMP*TPI
        ELSE IF(XTMP < -180.0_SP)THEN
          X2X1 =  360.0_SP*TPI+XTMP*TPI
        ELSE
          X2X1 = XTMP*TPI
      END IF
      Y3Y2  = (VY(TMP3)-VY(TMP2))*TPI
      Y1Y3  = (VY(IA)  -VY(TMP3))*TPI
      Y2Y1  = (VY(TMP2)-VY(IA))*TPI

        VCOL1 = 0.5_SP*GRAV_E(I1)*COEF*(-X3X2*COS(DEG2RAD*YC(I1))*DLTXE(I)-Y3Y2*DLTYE(I))*DT1(I1)/ART(I1)/ART1(IA)
        VCOL2 = 0.5_SP*GRAV_E(I1)*COEF*(-X1X3*COS(DEG2RAD*YC(I1))*DLTXE(I)-Y1Y3*DLTYE(I))*DT1(I1)/ART(I1)/ART1(IA)
        VCOL3 = 0.5_SP*GRAV_E(I1)*COEF*(-X2X1*COS(DEG2RAD*YC(I1))*DLTXE(I)-Y2Y1*DLTYE(I))*DT1(I1)/ART(I1)/ART1(IA)
!#     if defined (NORTHPOLE)
      ENDIF
!#     endif      

#     else
      X3X2  = VX(TMP3)-VX(TMP2)
      X1X3  = VX(IA)  -VX(TMP3)
      X2X1  = VX(TMP2)-VX(IA)
      Y3Y2  = VY(TMP3)-VY(TMP2)
      Y1Y3  = VY(IA)  -VY(TMP3)
      Y2Y1  = VY(TMP2)-VY(IA)
      VCOL1 = 0.5_SP*GRAV_E(I1)*COEF*(-X3X2*DLTXE(I)-Y3Y2*DLTYE(I))*DT1(I1)/ART(I1)/ART1(IA)
      VCOL2 = 0.5_SP*GRAV_E(I1)*COEF*(-X1X3*DLTXE(I)-Y1Y3*DLTYE(I))*DT1(I1)/ART(I1)/ART1(IA)
      VCOL3 = 0.5_SP*GRAV_E(I1)*COEF*(-X2X1*DLTXE(I)-Y2Y1*DLTYE(I))*DT1(I1)/ART(I1)/ART1(IA)
#     endif

      CALL MatSetValuesLocal(A_EL,1,PROW1-1,1,PCOL1-1,VCOL1,ADD_VALUES,IERR);CHKERRQ(IERR)
      CALL MatSetValuesLocal(A_EL,1,PROW1-1,1,PCOL2-1,VCOL2,ADD_VALUES,IERR);CHKERRQ(IERR)
      CALL MatSetValuesLocal(A_EL,1,PROW1-1,1,PCOL3-1,VCOL3,ADD_VALUES,IERR);CHKERRQ(IERR)

#     if defined (SPHERICAL)

!#     if defined (NORTHPOLE)
      IF (CELL_NORTHAREA(I1)==1) THEN
        VCOL1 = -0.5_SP*GRAV_E(I1)*COEF*(XXX1*DLTYE(I)-YYY1*DLTXE(I))*DT1(I1)/ART(I1)/ART1(IB)
        VCOL2 = -0.5_SP*GRAV_E(I1)*COEF*(XXX2*DLTYE(I)-YYY2*DLTXE(I))*DT1(I1)/ART(I1)/ART1(IB)
        VCOL3 = -0.5_SP*GRAV_E(I1)*COEF*(XXX3*DLTYE(I)-YYY3*DLTXE(I))*DT1(I1)/ART(I1)/ART1(IB)
      ELSE
!#     endif
        VCOL1 = -0.5_SP*GRAV_E(I1)*COEF*(-X3X2*COS(DEG2RAD*YC(I1))*DLTXE(I)-Y3Y2*DLTYE(I))*DT1(I1)/ART(I1)/ART1(IB)
        VCOL2 = -0.5_SP*GRAV_E(I1)*COEF*(-X1X3*COS(DEG2RAD*YC(I1))*DLTXE(I)-Y1Y3*DLTYE(I))*DT1(I1)/ART(I1)/ART1(IB)
        VCOL3 = -0.5_SP*GRAV_E(I1)*COEF*(-X2X1*COS(DEG2RAD*YC(I1))*DLTXE(I)-Y2Y1*DLTYE(I))*DT1(I1)/ART(I1)/ART1(IB)
!#     if defined (NORTHPOLE)
      ENDIF
!#     endif

#     else
      VCOL1 = -0.5_SP*GRAV_E(I1)*COEF*(-X3X2*DLTXE(I)-Y3Y2*DLTYE(I))*DT1(I1)/ART(I1)/ART1(IB)
      VCOL2 = -0.5_SP*GRAV_E(I1)*COEF*(-X1X3*DLTXE(I)-Y1Y3*DLTYE(I))*DT1(I1)/ART(I1)/ART1(IB)
      VCOL3 = -0.5_SP*GRAV_E(I1)*COEF*(-X2X1*DLTXE(I)-Y2Y1*DLTYE(I))*DT1(I1)/ART(I1)/ART1(IB)
#     endif

      CALL MatSetValuesLocal(A_EL,1,PROW2-1,1,PCOL1-1,VCOL1,ADD_VALUES,IERR);CHKERRQ(IERR)
      CALL MatSetValuesLocal(A_EL,1,PROW2-1,1,PCOL2-1,VCOL2,ADD_VALUES,IERR);CHKERRQ(IERR)
      CALL MatSetValuesLocal(A_EL,1,PROW2-1,1,PCOL3-1,VCOL3,ADD_VALUES,IERR);CHKERRQ(IERR)

#     if defined (WET_DRY)
      ENDIF
#     endif
!!$    ENDIF
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
# if defined (MULTIPROCESSOR)
  IF(PAR) CALL MPI_BARRIER(MPI_FVCOM_GROUP,IERR)
# endif

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

  DO I=1, M
    IF(ISONB_3(I)==3) THEN
      PROW1  = ALO_2_PLO_NODE(I)
      IF(PROW1 > N_VERTS) CYCLE
      VCOL1 = 0.0D0
      DO J=1, NTSN(I)-1
        PCOL1 = ALO_2_PLO_NODE(NBSN(I,J))
        CALL MatSetValuesLocal(A_EL,1,PROW1-1,1,PCOL1-1,VCOL1,INSERT_VALUES,IERR);CHKERRQ(IERR)
      ENDDO
      VCOL1 = 1.0D0
      PROW1 = ALO_2_PLO_NODE(I)
      CALL MatSetValuesLocal(A_EL,1,PROW1-1,1,PROW1-1,VCOL1,INSERT_VALUES,IERR);CHKERRQ(IERR)
    ENDIF
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

!# if defined(SPHERICAL) && (NORTHPOLE)
# if defined(SPHERICAL) 
  IF(NODE_NORTHPOLE/=0) THEN
    PROW1 = ALO_2_PLO_NODE(NODE_NORTHPOLE)
    IF(PROW1<=N_VERTS) THEN
      VCOL1 = 0.0D0
      NODE  = NODE_NORTHPOLE
      DO J=1, NTSN(NODE)-1
        PCOL1 = ALO_2_PLO_NODE(NBSN(NODE,J))
        CALL MatSetValuesLocal(A_EL,1,PROW1-1,1,PCOL1-1,VCOL1,INSERT_VALUES,IERR);CHKERRQ(IERR)
      ENDDO
      PROW1 = ALO_2_PLO_NODE(NODE)
      CALL MatSetValuesLocal(A_EL,1,PROW1-1,1,PROW1-1,VCOL1,INSERT_VALUES,IERR);CHKERRQ(IERR)
    ENDIF
  ENDIF
# endif

  CALL MatAssemblyBegin(A_EL,MAT_FINAL_ASSEMBLY,IERR);CHKERRQ(IERR)
  CALL MatAssemblyEnd(A_EL,MAT_FINAL_ASSEMBLY,IERR);CHKERRQ(IERR)

!# if defined(SPHERICAL) && (NORTHPOLE)
# if defined(SPHERICAL) 
  COEF  = DTI**2*IFCETA**2
  DO II=1,NPCV
    I = NCEDGE_LST(II)
    I1  = NTRG(I)
    IA  = NIEC(I,1)
    IB  = NIEC(I,2)

    IF(IA==NODE_NORTHPOLE .OR. IB==NODE_NORTHPOLE) THEN
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

      VX1_TMP = REARTH * COS(VY(IA)*DEG2RAD) * COS(VX(IA)*DEG2RAD)&
              * 2._SP /(1._SP+sin(VY(IA)*DEG2RAD)) 
      VY1_TMP = REARTH * COS(VY(IA)*DEG2RAD) * SIN(VX(IA)*DEG2RAD)&
              * 2._SP /(1._SP+sin(VY(IA)*DEG2RAD))

      VX2_TMP = REARTH * COS(VY(IB)*DEG2RAD) * COS(VX(IB)*DEG2RAD)&
              * 2._SP /(1._SP+sin(VY(IB)*DEG2RAD)) 
      VY2_TMP = REARTH * COS(VY(IB)*DEG2RAD) * SIN(VX(IB)*DEG2RAD)&
              * 2._SP /(1._SP+sin(VY(IB)*DEG2RAD))

      VX3_TMP = REARTH * COS(VY(TMP3)*DEG2RAD) * COS(VX(TMP3)*DEG2RAD)&
              * 2._SP /(1._SP+sin(VY(TMP3)*DEG2RAD))
      VY3_TMP = REARTH * COS(VY(TMP3)*DEG2RAD) * SIN(VX(TMP3)*DEG2RAD)&
              * 2._SP /(1._SP+sin(VY(TMP3)*DEG2RAD))

      X3X2  = VX3_TMP-VX2_TMP
      X1X3  = VX1_TMP-VX3_TMP
      X2X1  = VX2_TMP-VX1_TMP
      Y3Y2  = VY3_TMP-VY2_TMP
      Y1Y3  = VY1_TMP-VY3_TMP
      Y2Y1  = VY2_TMP-VY1_TMP

      VX1_TMP = REARTH * COS(YIJE(I,1)*DEG2RAD) * COS(XIJE(I,1)*DEG2RAD)!&
!              * 2._SP /(1._SP+sin(YIJE(I,1)*DEG2RAD))
      VY1_TMP = REARTH * COS(YIJE(I,1)*DEG2RAD) * SIN(XIJE(I,1)*DEG2RAD)!&
!              * 2._SP /(1._SP+sin(YIJE(I,1)*DEG2RAD))

      VX2_TMP = REARTH * COS(YIJE(I,2)*DEG2RAD) * COS(XIJE(I,2)*DEG2RAD)!&
!              * 2._SP /(1._SP+sin(YIJE(I,2)*DEG2RAD))
      VY2_TMP = REARTH * COS(YIJE(I,2)*DEG2RAD) * SIN(XIJE(I,2)*DEG2RAD)!&
!              * 2._SP /(1._SP+sin(YIJE(I,2)*DEG2RAD))

      DLTXE_TMP = VX2_TMP-VX1_TMP
      DLTYE_TMP = VY2_TMP-VY1_TMP

      IF(IA == NODE_NORTHPOLE) THEN
        VCOL1 = 0.5_SP*GRAV_E(I1)*COEF*(-X3X2*DLTXE_TMP-Y3Y2*DLTYE_TMP)*DT1(I1)/ART(I1)/ART1(IA)/(2._SP /(1._SP+sin(YC(I1)*DEG2RAD)))
        VCOL2 = 0.5_SP*GRAV_E(I1)*COEF*(-X1X3*DLTXE_TMP-Y1Y3*DLTYE_TMP)*DT1(I1)/ART(I1)/ART1(IA)/(2._SP /(1._SP+sin(YC(I1)*DEG2RAD)))
        VCOL3 = 0.5_SP*GRAV_E(I1)*COEF*(-X2X1*DLTXE_TMP-Y2Y1*DLTYE_TMP)*DT1(I1)/ART(I1)/ART1(IA)/(2._SP /(1._SP+sin(YC(I1)*DEG2RAD)))
        CALL MatSetValuesLocal(A_EL,1,PROW1-1,1,PCOL1-1,VCOL1,ADD_VALUES,IERR);CHKERRQ(IERR)
        CALL MatSetValuesLocal(A_EL,1,PROW1-1,1,PCOL2-1,VCOL2,ADD_VALUES,IERR);CHKERRQ(IERR)
        CALL MatSetValuesLocal(A_EL,1,PROW1-1,1,PCOL3-1,VCOL3,ADD_VALUES,IERR);CHKERRQ(IERR)
      ELSE IF (IB == NODE_NORTHPOLE) THEN
        VCOL1 = -0.5_SP*GRAV_E(I1)*COEF*(-X3X2*DLTXE_TMP-Y3Y2*DLTYE_TMP)*DT1(I1)/ART(I1)/ART1(IB)/(2._SP /(1._SP+sin(YC(I1)*DEG2RAD)))
        VCOL2 = -0.5_SP*GRAV_E(I1)*COEF*(-X1X3*DLTXE_TMP-Y1Y3*DLTYE_TMP)*DT1(I1)/ART(I1)/ART1(IB)/(2._SP /(1._SP+sin(YC(I1)*DEG2RAD)))
        VCOL3 = -0.5_SP*GRAV_E(I1)*COEF*(-X2X1*DLTXE_TMP-Y2Y1*DLTYE_TMP)*DT1(I1)/ART(I1)/ART1(IB)/(2._SP /(1._SP+sin(YC(I1)*DEG2RAD)))
        CALL MatSetValuesLocal(A_EL,1,PROW2-1,1,PCOL1-1,VCOL1,ADD_VALUES,IERR);CHKERRQ(IERR)
        CALL MatSetValuesLocal(A_EL,1,PROW2-1,1,PCOL2-1,VCOL2,ADD_VALUES,IERR);CHKERRQ(IERR)
        CALL MatSetValuesLocal(A_EL,1,PROW2-1,1,PCOL3-1,VCOL3,ADD_VALUES,IERR);CHKERRQ(IERR)
      ENDIF
    ENDIF

  ENDDO

  IF(NODE_NORTHPOLE/=0) THEN
    PROW1 = ALO_2_PLO_NODE(NODE_NORTHPOLE)
    IF(PROW1<=N_VERTS) THEN
      VCOL1 = 1.0D0
      CALL MatSetValuesLocal(A_EL,1,PROW1-1,1,PROW1-1,VCOL1,ADD_VALUES,IERR);CHKERRQ(IERR)
    ENDIF
  ENDIF
 
  CALL MatAssemblyBegin(A_EL,MAT_FINAL_ASSEMBLY,IERR);CHKERRQ(IERR)
  CALL MatAssemblyEnd(A_EL,MAT_FINAL_ASSEMBLY,IERR);CHKERRQ(IERR)
# endif

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

  DO I = 1, NT
    EL1(I) = (EL(NV(I,1))+EL(NV(I,2))+EL(NV(I,3)))/3.0_SP
  END DO

#   if defined (ICE_EMBEDDING)
  DO I = 1, NT
    QTHICE1(I) = (QTHICE(NV(I,1))+QTHICE(NV(I,2))+QTHICE(NV(I,3)))/3.0_SP   !!yding
  END DO

  D   = H + EL - QTHICE     
  D1  = H1 + EL1 - QTHICE1

!!---------Adjust total depth if the node/cell is dry------------!!yding
     DO I = 1, MT
      IF (D(I)-MIN_DEPTH<1.0E-5_SP) THEN
       D(I) = MIN_DEPTH
      END IF
     END DO

    DO I = 1, NT
      IF (D1(I)-MIN_DEPTH<1.0E-5_SP) THEN
       D1(I) = MIN_DEPTH
      END IF
     END DO

!!=====================yding change end

#   else

 D   = H + EL
 D1  = H1 + EL1

#   endif



  DO I=1, NE
    IA=IEC(I,1)
    IB=IEC(I,2)
    J1=IENODE(I,1)
    J2=IENODE(I,2)
    ELIJ=IFCETA*0.5_SP*(EL(J1)+EL(J2))

    ELIJ2=(1.0_SP-IFCETA)*0.5_SP*(ET(J1)+ET(J2))
#   if defined (AIR_PRESSURE)
    ELIJ2=ELIJ2-( (1.0_SP-IFCETA)*0.5_SP*(EL_AIR(J1)+EL_AIR(J2))+IFCETA*0.5_SP*(ELF_AIR(J1)+ELF_AIR(J2)) )
#   endif
#   if defined (EQUI_TIDE)
    ELIJ2=ELIJ2-( (1.0_SP-IFCETA)*0.5_SP*(EL_EQI(J1)+EL_EQI(J2))+IFCETA*0.5_SP*(ELF_EQI(J1)+ELF_EQI(J2)) )
#   endif
#   if defined (ATMO_TIDE)
    ELIJ2=ELIJ2-( (1.0_SP-IFCETA)*0.5_SP*(EL_ATMO(J1)+EL_ATMO(J2))+IFCETA*0.5_SP*(ELF_ATMO(J1)+ELF_ATMO(J2)) )
#   endif    

    ELIJ3= (1.0_SP-IFCETA)*0.5_SP*(ET(J1)+ET(J2)) + IFCETA*0.5_SP*(EL(J1)+EL(J2))
#   if defined (AIR_PRESSURE)
    ELIJ3=ELIJ3-( (1.0_SP-IFCETA)*0.5_SP*(EL_AIR(J1)+EL_AIR(J2))+IFCETA*0.5_SP*(ELF_AIR(J1)+ELF_AIR(J2)) )
#   endif
#   if defined (EQUI_TIDE)
    ELIJ3=ELIJ3-( (1.0_SP-IFCETA)*0.5_SP*(EL_EQI(J1)+EL_EQI(J2))+IFCETA*0.5_SP*(ELF_EQI(J1)+ELF_EQI(J2)) )
#   endif
#   if defined (ATMO_TIDE)
    ELIJ3=ELIJ3-( (1.0_SP-IFCETA)*0.5_SP*(EL_ATMO(J1)+EL_ATMO(J2))+IFCETA*0.5_SP*(ELF_ATMO(J1)+ELF_ATMO(J2)) )
#   endif

#   if defined (TWO_D_MODEL)

#   if defined (SPHERICAL)
    XTMP  = VX(J2)*TPI-VX(J1)*TPI
    XTMP1 = VX(J2)-VX(J1)
    IF(XTMP1 >  180.0_SP)THEN
      XTMP = -360.0_SP*TPI+XTMP
    ELSE IF(XTMP1 < -180.0_SP)THEN
      XTMP =  360.0_SP*TPI+XTMP
    END IF

!    ADVUA(IA)=ADVUA(IA)-GRAV_E(IA)*DT1(IA)*ELIJ*DLTYC(I)
!    ADVVA(IA)=ADVVA(IA)+GRAV_E(IA)*DT1(IA)*ELIJ*XTMP*COS(DEG2RAD*YC(IA)) 
!    ADVUA(IB)=ADVUA(IB)+GRAV_E(IB)*DT1(IB)*ELIJ*DLTYC(I)
!    ADVVA(IB)=ADVVA(IB)-GRAV_E(IB)*DT1(IB)*ELIJ*XTMP*COS(DEG2RAD*YC(IB)) 

    ADVUA(IA)=ADVUA(IA)+GRAV_E(IA)*DT1(IA)*ELIJ2*DLTYC(I)  &
             -F_ALFA(IA)*GRAV_E(IA)*DT1(IA)*ELIJ3*DLTYC(I)
    ADVVA(IA)=ADVVA(IA)-GRAV_E(IA)*DT1(IA)*ELIJ2*XTMP*COS(DEG2RAD*YC(IA)) &
             +F_ALFA(IA)*GRAV_E(IA)*DT1(IA)*ELIJ3*XTMP*COS(DEG2RAD*YC(IA))

    ADVUA(IB)=ADVUA(IB)-GRAV_E(IB)*DT1(IB)*ELIJ2*DLTYC(I) &
             +F_ALFA(IB)*GRAV_E(IB)*DT1(IB)*ELIJ3*DLTYC(I)
    ADVVA(IB)=ADVVA(IB)+GRAV_E(IB)*DT1(IB)*ELIJ2*XTMP*COS(DEG2RAD*YC(IB)) &
             -F_ALFA(IB)*GRAV_E(IB)*DT1(IB)*ELIJ3*XTMP*COS(DEG2RAD*YC(IB))
#   else
    ADVUA(IA)=ADVUA(IA)-GRAV_E(IA)*DT1(IA)*ELIJ*DLTYC(I)
    ADVVA(IA)=ADVVA(IA)+GRAV_E(IA)*DT1(IA)*ELIJ*DLTXC(I)
    ADVUA(IB)=ADVUA(IB)+GRAV_E(IB)*DT1(IB)*ELIJ*DLTYC(I)
    ADVVA(IB)=ADVVA(IB)-GRAV_E(IB)*DT1(IB)*ELIJ*DLTXC(I)
#   endif

#   else
    
    DO K=1, KBM1
#     if defined (SPHERICAL)
      XTMP  = VX(J2)*TPI-VX(J1)*TPI
      XTMP1 = VX(J2)-VX(J1)
      IF(XTMP1 >  180.0_SP)THEN
        XTMP = -360.0_SP*TPI+XTMP
      ELSE IF(XTMP1 < -180.0_SP)THEN
        XTMP =  360.0_SP*TPI+XTMP
      END IF

!      XFLUX3(IA,K)=XFLUX3(IA,K)-GRAV_E(IA)*DT1(IA)*DZ1(IA,K)*ELIJ*DLTYC(I)
!      YFLUX3(IA,K)=YFLUX3(IA,K)+GRAV_E(IA)*DT1(IA)*DZ1(IA,K)*ELIJ*XTMP*COS(DEG2RAD*YC(IA)) 
!      XFLUX3(IB,K)=XFLUX3(IB,K)+GRAV_E(IB)*DT1(IB)*DZ1(IB,K)*ELIJ*DLTYC(I)
!      YFLUX3(IB,K)=YFLUX3(IB,K)-GRAV_E(IB)*DT1(IB)*DZ1(IB,K)*ELIJ*XTMP*COS(DEG2RAD*YC(IB)) 

      XFLUX3(IA,K)=XFLUX3(IA,K)+GRAV_E(IA)*DT1(IA)*DZ1(IA,K)*ELIJ2*DLTYC(I) &
                  -F_ALFA(IA)*GRAV_E(IA)*DT1(IA)*DZ1(IA,K)*ELIJ3*DLTYC(I)
      YFLUX3(IA,K)=YFLUX3(IA,K)-GRAV_E(IA)*DT1(IA)*DZ1(IA,K)*ELIJ2*XTMP*COS(DEG2RAD*YC(IA)) &
                  +F_ALFA(IA)*GRAV_E(IA)*DT1(IA)*DZ1(IA,K)*ELIJ3*XTMP*COS(DEG2RAD*YC(IA))

      XFLUX3(IB,K)=XFLUX3(IB,K)-GRAV_E(IB)*DT1(IB)*DZ1(IB,K)*ELIJ2*DLTYC(I) &
                  +F_ALFA(IB)*GRAV_E(IB)*DT1(IB)*DZ1(IB,K)*ELIJ3*DLTYC(I)
      YFLUX3(IB,K)=YFLUX3(IB,K)+GRAV_E(IB)*DT1(IB)*DZ1(IB,K)*ELIJ2*XTMP*COS(DEG2RAD*YC(IB)) &
                  -F_ALFA(IB)*GRAV_E(IB)*DT1(IB)*DZ1(IB,K)*ELIJ3*XTMP*COS(DEG2RAD*YC(IB))

#     else
      XFLUX3(IA,K)=XFLUX3(IA,K)-GRAV_E(IA)*DT1(IA)*DZ1(IA,K)*ELIJ*DLTYC(I)
      YFLUX3(IA,K)=YFLUX3(IA,K)+GRAV_E(IA)*DT1(IA)*DZ1(IA,K)*ELIJ*DLTXC(I)
      XFLUX3(IB,K)=XFLUX3(IB,K)+GRAV_E(IB)*DT1(IB)*DZ1(IB,K)*ELIJ*DLTYC(I)
      YFLUX3(IB,K)=YFLUX3(IB,K)-GRAV_E(IB)*DT1(IB)*DZ1(IB,K)*ELIJ*DLTXC(I)
#     endif
    ENDDO

#   endif
  ENDDO

!# if defined (SPHERICAL) && (NORTHPOLE)
# if defined (SPHERICAL) 
  DO II=1,NPE
    I=NPEDGE_LST(II)
    IA=IEC(I,1)
    IB=IEC(I,2)
    J1=IENODE(I,1)
    J2=IENODE(I,2)
    ELIJ=IFCETA*0.5_SP*(EL(J1)+EL(J2))

    DO K=1, KBM1
      VX1_TMP = REARTH * COS(VY(IENODE(I,1))*DEG2RAD) * COS(VX(IENODE(I,1))*DEG2RAD) &
                * 2._SP /(1._SP+sin(VY(IENODE(I,1))*DEG2RAD))
      VY1_TMP = REARTH * COS(VY(IENODE(I,1))*DEG2RAD) * SIN(VX(IENODE(I,1))*DEG2RAD) &
                * 2._SP /(1._SP+sin(VY(IENODE(I,1))*DEG2RAD))

      VX2_TMP = REARTH * COS(VY(IENODE(I,2))*DEG2RAD) * COS(VX(IENODE(I,2))*DEG2RAD) &
                * 2._SP /(1._SP+sin(VY(IENODE(I,2))*DEG2RAD))
      VY2_TMP = REARTH * COS(VY(IENODE(I,2))*DEG2RAD) * SIN(VX(IENODE(I,2))*DEG2RAD) &
                * 2._SP /(1._SP+sin(VY(IENODE(I,2))*DEG2RAD))

      DLTXC_TMP = VX2_TMP-VX1_TMP
      DLTYC_TMP = VY2_TMP-VY1_TMP

      IF(CELL_NORTHAREA(IA) == 1 .AND. CELL_NORTHAREA(IB) == 1)THEN
        XFLUX3_NP(IA,K)=XFLUX3_NP(IA,K)-GRAV_E(IA)*DT1(IA)*DZ1(IA,K)*ELIJ*DLTYC_TMP/  &
                       (2._SP /(1._SP+sin(YC(IA)*DEG2RAD)))
        YFLUX3_NP(IA,K)=YFLUX3_NP(IA,K)+GRAV_E(IA)*DT1(IA)*DZ1(IA,K)*ELIJ*DLTXC_TMP/  &
                       (2._SP /(1._SP+sin(YC(IA)*DEG2RAD)))
        XFLUX3_NP(IB,K)=XFLUX3_NP(IB,K)+GRAV_E(IB)*DT1(IB)*DZ1(IB,K)*ELIJ*DLTYC_TMP/  &
                       (2._SP /(1._SP+sin(YC(IB)*DEG2RAD)))
        YFLUX3_NP(IB,K)=YFLUX3_NP(IB,K)-GRAV_E(IB)*DT1(IB)*DZ1(IB,K)*ELIJ*DLTXC_TMP/  &
                       (2._SP /(1._SP+sin(YC(IB)*DEG2RAD)))
      ELSE IF(CELL_NORTHAREA(IA) == 1 .AND. CELL_NORTHAREA(IB) /= 1)THEN
        XFLUX3_NP(IA,K)=XFLUX3_NP(IA,K)-GRAV_E(IA)*DT1(IA)*DZ1(IA,K)*ELIJ*DLTYC_TMP/   &
                      (2._SP /(1._SP+sin(YC(IA)*DEG2RAD)))
        YFLUX3_NP(IA,K)=YFLUX3_NP(IA,K)+GRAV_E(IA)*DT1(IA)*DZ1(IA,K)*ELIJ*DLTXC_TMP/   &
                      (2._SP /(1._SP+sin(YC(IA)*DEG2RAD)))
      ELSE IF(CELL_NORTHAREA(IB) == 1 .AND. CELL_NORTHAREA(IA) /= 1)THEN
        XFLUX3_NP(IB,K)=XFLUX3_NP(IB,K)+GRAV_E(IB)*DT1(IB)*DZ1(IB,K)*ELIJ*DLTYC_TMP/  &
                      (2._SP /(1._SP+sin(YC(IB)*DEG2RAD)))
        YFLUX3_NP(IB,K)=YFLUX3_NP(IB,K)-GRAV_E(IB)*DT1(IB)*DZ1(IB,K)*ELIJ*DLTXC_TMP/  &
                      (2._SP /(1._SP+sin(YC(IB)*DEG2RAD)))
      END IF
    ENDDO

  ENDDO
# endif

  DO I=1, N

# if defined (WET_DRY)
  IF(ISWETCT(I)*ISWETC(I) == 1) THEN
# endif
    
# if defined (TWO_D_MODEL)

#   if defined (SPHERICAL)
    UAF(I) = UA(I)*DT1(I)/D1(I) - DTI*ADVUA(I)/ART(I)/D1(I) + DTI*( -WUSURF(I) + WUBOT(I) )/D1(I) - BEDF*UA_N(I)*DT1(I)/D1(I) 
    VAF(I) = VA(I)*DT1(I)/D1(I) - DTI*ADVVA(I)/ART(I)/D1(I) + DTI*( -WVSURF(I) + WVBOT(I) )/D1(I) - BEDF*VA_N(I)*DT1(I)/D1(I)
#   else
    UAF(I) = UA(I)*DT1(I)/D1(I) - DTI*ADVUA(I)/ART(I)/D1(I) + DTI*( -WUSURF(I) + WUBOT(I) )/D1(I) - BEDF*UA_N(I)*DT1(I)/D1(I)
    VAF(I) = VA(I)*DT1(I)/D1(I) - DTI*ADVVA(I)/ART(I)/D1(I) + DTI*( -WVSURF(I) + WVBOT(I) )/D1(I) - BEDF*VA_N(I)*DT1(I)/D1(I)
#   endif

# else

    DO K=1, KBM1
#     if defined (SPHERICAL)
!#     if defined (NORTHPOLE)
      IF(CELL_NORTHAREA(I)==1) THEN      
        U_TMP = -V(I,K)*COS(XC(I)*DEG2RAD)-U(I,K)*SIN(XC(I)*DEG2RAD)
        V_TMP = -V(I,K)*SIN(XC(I)*DEG2RAD)+U(I,K)*COS(XC(I)*DEG2RAD)  
        UF_TMP = U_TMP*DT1(I)/D1(I) - DTI*XFLUX3_NP(I,K)/ART(I)/D1(I)/DZ1(I,K)
        VF_TMP = V_TMP*DT1(I)/D1(I) - DTI*YFLUX3_NP(I,K)/ART(I)/D1(I)/DZ1(I,K)

        UF(I,K) = UF_TMP
        VF(I,K) = VF_TMP  
      ELSE
!#     endif
       UF(I,K) = U(I,K)*DT1(I)/D1(I) - DTI*XFLUX3(I,K)/ART(I)/D1(I)/DZ1(I,K) - BEDF*U_N(I,K)*DT1(I)/D1(I)
       VF(I,K) = V(I,K)*DT1(I)/D1(I) - DTI*YFLUX3(I,K)/ART(I)/D1(I)/DZ1(I,K) - BEDF*V_N(I,K)*DT1(I)/D1(I)
!#     if defined (NORTHPOLE)
      ENDIF
!#     endif

#     else
      UF(I,K) = U(I,K)*DT1(I)/D1(I) - DTI*XFLUX3(I,K)/ART(I)/D1(I)/DZ1(I,K) - BEDF*U_N(I,K)*DT1(I)/D1(I)
      VF(I,K) = V(I,K)*DT1(I)/D1(I) - DTI*YFLUX3(I,K)/ART(I)/D1(I)/DZ1(I,K) - BEDF*V_N(I,K)*DT1(I)/D1(I)
#     endif
    ENDDO
    
# endif

# if defined (WET_DRY)
  ELSE
#   if defined (TWO_D_MODEL)
    UAF(I) = 0.0_SP
    VAF(I) = 0.0_SP
#   else
    DO K=1, KBM1
      UF(I,K) = 0.0_SP
      VF(I,K) = 0.0_SP
    ENDDO  
#   endif
  ENDIF
# endif

  ENDDO

# if defined (TWO_D_MODEL)
!old:  UAF = UAF-CC_SPONGE*UAF
!old:  VAF = VAF-CC_SPONGE*VAF
! ---- new: Karsten Lettmann: 2012.06.25 -------
    UAF = UAF/(1.0_SP+CC_SPONGE*UAF**2.0_SP)
    VAF = VAF/(1.0_SP+CC_SPONGE*VAF**2.0_SP)
! ------- end new -------------------------------
# endif

  RETURN
  END SUBROUTINE SEMI_IMPLICIT_EL

  SUBROUTINE UPDATES_SEMI

  USE MOD_PREC, ONLY : SP, DP
  USE LIMS
  USE CONTROL 
  USE MOD_OBCS
  USE ALL_VARS
  USE BCS
# if defined(MULTIPROCESSOR)
  USE MOD_PAR
# endif
# if defined (SPHERICAL)   
  USE MOD_SPHERICAL
  USE MOD_NORTHPOLE,   ONLY :  NPCV,  NCEDGE_LST, NODE_NORTHPOLE
# endif   

  IMPLICIT NONE
  INTEGER  II, I, J, JJ, J1, J2, IK, K, I1, IA, IB, TMP2, TMP3
  REAL(SP) DIJ, DIJ1, UIJ, VIJ, EXFLUX, COEF
  REAL(SP) :: XFLUX(0:MT), UF_AVG(0:NT), VF_AVG(0:NT)

  REAL(SP) :: UIJ_TMP,VIJ_TMP,EXFLUX_TMP,FLUX_TMP
  REAL(SP) :: DLTXE_TMP,DLTYE_TMP
  REAL(SP) :: VX1_TMP,VX2_TMP,VY1_TMP,VY2_TMP

  UF_AVG = 0.0_SP
  VF_AVG = 0.0_SP
# if defined(MULTIPROCESSOR)
  IF(PAR) CALL AEXCHANGE(EC,MYID,NPROCS,UF,VF)
# endif
  DO I=1,NT
    UF_AVG(I) = SUM(UF(I,1:KBM1)*DZ1(I,1:KBM1))
    VF_AVG(I) = SUM(VF(I,1:KBM1)*DZ1(I,1:KBM1))
  END DO

  XFLUX   = 0.0_SP
  DO I=1,NCV
    I1=NTRG(I)
    IA=NIEC(I,1)
    IB=NIEC(I,2)
    DIJ=(DT(NV(I1,1))+DT(NV(I1,2))+DT(NV(I1,3)))/3.0_SP
    DIJ1=(D(NV(I1,1))+D(NV(I1,2))+D(NV(I1,3)))/3.0_SP

    UIJ=UF_AVG(I1)
    VIJ=VF_AVG(I1)
    EXFLUX=(1.0_SP-IFCETA)*DIJ*(-UA(I1)*DLTYE(I)+VA(I1)*DLTXE(I))+IFCETA*DIJ1*(-UIJ*DLTYE(I)+VIJ*DLTXE(I))
    XFLUX(IA)=XFLUX(IA)-EXFLUX
    XFLUX(IB)=XFLUX(IB)+EXFLUX
  ENDDO
!# if defined (SPHERICAL) && (NORTHPOLE)
# if defined (SPHERICAL) 
  IF(NODE_NORTHPOLE/=0) XFLUX(NODE_NORTHPOLE) = 0.0_SP

  DO II=1,NPCV
    I = NCEDGE_LST(II)
    I1  = NTRG(I)
    IA  = NIEC(I,1)
    IB  = NIEC(I,2)

    DIJ = (DT(NV(I1,1))+DT(NV(I1,2))+DT(NV(I1,3)))/3.0_SP
    DIJ1= (D(NV(I1,1))+D(NV(I1,2))+D(NV(I1,3)))/3.0_SP
    UIJ = -VA(I1)*COS(XC(I1)*DEG2RAD)-UA(I1)*SIN(XC(I1)*DEG2RAD)
    VIJ = -VA(I1)*SIN(XC(I1)*DEG2RAD)+UA(I1)*COS(XC(I1)*DEG2RAD)
    UIJ_TMP = -VF_AVG(I1)*COS(XC(I1)*DEG2RAD)-UF_AVG(I1)*SIN(XC(I1)*DEG2RAD)
    VIJ_TMP = -VF_AVG(I1)*SIN(XC(I1)*DEG2RAD)+UF_AVG(I1)*COS(XC(I1)*DEG2RAD)

    IF(IA == NODE_NORTHPOLE .OR. IB == NODE_NORTHPOLE)THEN
      VX1_TMP = REARTH * COS(YIJE(I,1)*DEG2RAD) * COS(XIJE(I,1)*DEG2RAD)!&
!                * 2._SP /(1._SP+sin(YIJE(I,1)*DEG2RAD))
      VY1_TMP = REARTH * COS(YIJE(I,1)*DEG2RAD) * SIN(XIJE(I,1)*DEG2RAD)!&
!                * 2._SP /(1._SP+sin(YIJE(I,1)*DEG2RAD))

      VX2_TMP = REARTH * COS(YIJE(I,2)*DEG2RAD) * COS(XIJE(I,2)*DEG2RAD)!&
!                * 2._SP /(1._SP+sin(YIJE(I,2)*DEG2RAD))
      VY2_TMP = REARTH * COS(YIJE(I,2)*DEG2RAD) * SIN(XIJE(I,2)*DEG2RAD)!&
!                * 2._SP /(1._SP+sin(YIJE(I,2)*DEG2RAD))

      DLTXE_TMP = VX2_TMP-VX1_TMP
      DLTYE_TMP = VY2_TMP-VY1_TMP

      EXFLUX_TMP =(1.0_SP-IFCETA)*DIJ*(-UIJ*DLTYE_TMP+VIJ*DLTXE_TMP)+IFCETA*DIJ1*(-UIJ_TMP*DLTYE_TMP+VIJ_TMP*DLTXE_TMP)
    END IF

    IF(IA == NODE_NORTHPOLE) THEN
      XFLUX(IA) = XFLUX(IA)-EXFLUX_TMP
    ELSE IF(IB == NODE_NORTHPOLE)THEN
      XFLUX(IB) = XFLUX(IB)+EXFLUX_TMP
    END IF

  END DO
# endif


  IF (PRECIPITATION_ON) THEN       !!yding
! this part may need be changed to semi-implicit form 
!  XFLUX = XFLUX+(QEVAP-QPREC)*ROFVROS*ART1
  XFLUX = XFLUX-(QEVAP+QPREC)*ROFVROS*ART1
  END IF

!!======================yding======================
#  if defined (ICE)

   XFLUX = XFLUX-QPICE*ROFVROS*ART1

#  endif



#   if defined (ICE_EMBEDDING)
    DO I=1,M
      XFLUX(I)=XFLUX(I) + QFLICE(I)*(917.0_SP/((RHO1(I,1)+1.0_SP)*1.0E3_SP))*ART1(I)  !!yding
    END DO

#   endif

  IF(GROUNDWATER_ON)THEN
    DO I=1,M
!      XFLUX(I)=XFLUX(I)+BFWDIS(I)*ROFVROS*ART1(I)
      XFLUX(I)=XFLUX(I)-BFWDIS(I)
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

  IF(IOBCN > 0) THEN
    DO I=1,IOBCN
      XFLUX_OBCN(I)=XFLUX(I_OBC_N(I))
      XFLUX(I_OBC_N(I)) = 0.0_SP
    END DO
  END IF

  UA = UF_AVG
  VA = VF_AVG
  DTFA = D

  RETURN
  END SUBROUTINE UPDATES_SEMI

  SUBROUTINE PETSc_SETICS_EL
  USE ALL_VARS, ONLY: MYID, ET
  USE MOD_PETSc, ONLY: USE_LAST, X_EL, XL_EL, L2G_EL, N_VERTS, PLO_2_ALO_NODE
  IMPLICIT NONE
  INTEGER :: I, IK
  PetscReal :: QTERM
  PetscInt  :: IERR
  PetscScalar :: ZERO   =  0.0D0
  CHARACTER(LEN=20) :: SUBNAME = 'PETSc_SETICS'
  
  IF(.NOT.USE_LAST) THEN
    CALL VecSet(X_EL,ZERO,IERR);CHKERRQ(IERR)
    RETURN
  ENDIF
  
  DO I=1,N_VERTS
    IK = PLO_2_ALO_NODE(I)
    QTERM = ET(IK)
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

  SUBROUTINE UV2D_SBD
  USE MOD_PREC, ONLY : SP, DP
  USE ALL_VARS
  USE BCS
  USE MOD_OBCS
  IMPLICIT NONE

  INTEGER  :: I,J,I1,J1,J2
  REAL(SP) :: ALPHA1,UNTMP,VNTMP,UI,VI

  UA_N = 0.0_SP
  VA_N = 0.0_SP

  DO I=1, N

    IF(ISBCE(I) == 1) THEN
      ALPHA1=ALPHA(I)
      IF(NUMQBC > 0) THEN
        IF(RIVER_INFLOW_LOCATION == 'node') THEN
          DO J=1,NUMQBC
            I1=INODEQ(J)
            J1=NBVE(I1,1)
            J2=NBVE(I1,NTVE(I1))
            IF((I == J1).OR.(I == J2)) THEN
              UNTMP=UA(I)*COS(ANGLEQ(J))+VA(I)*SIN(ANGLEQ(J))
              IF(UNTMP<0.0_SP) THEN
                UA_N(I) = UNTMP*COS(ANGLEQ(J))
                VA_N(I) = UNTMP*SIN(ANGLEQ(J))
              ENDIF
              GOTO 21
            END IF
          END DO
        ELSE IF(RIVER_INFLOW_LOCATION == 'edge') THEN
          DO J=1,NUMQBC
            J1=ICELLQ(J)
            IF(I == J1) THEN
              UNTMP=UA(I)*COS(ANGLEQ(J))+VA(I)*SIN(ANGLEQ(J))
              IF(UNTMP<0.0_SP) THEN
                UA_N(I) = UNTMP*COS(ANGLEQ(J))
                VA_N(I) = UNTMP*SIN(ANGLEQ(J))
              ENDIF
              GOTO 21
            END IF
          END DO
        END IF
      END IF

      UNTMP = UA(I)*COS(ALPHA1)+VA(I)*SIN(ALPHA1)
      IF(UNTMP>0.0_SP) THEN
        UA_N(I)= UNTMP*COS(ALPHA1)
        VA_N(I)= UNTMP*SIN(ALPHA1)
      ENDIF

21  CONTINUE
    END IF
  END DO

  RETURN
  END SUBROUTINE UV2D_SBD

  SUBROUTINE UV3D_SBD
  USE MOD_PREC, ONLY : SP, DP
  USE ALL_VARS
  USE BCS
  USE MOD_OBCS
# if defined (THIN_DAM)
  USE MOD_DAM
# endif
  IMPLICIT NONE

  INTEGER  :: I,J,K,I1,J1,J2
  REAL(SP) :: ALPHA1,UNTMP,VNTMP,UI,VI

  U_N = 0.0_SP
  V_N = 0.0_SP

  DO I=1, N
    DO K=1, KBM1

      IF(ISBCE(I).EQ.1) THEN
        ALPHA1=ALPHA(I)
        IF(NUMQBC.GE.1) THEN
          IF(RIVER_INFLOW_LOCATION .EQ. 'node') THEN
            DO J=1,NUMQBC
              I1=INODEQ(J)
              J1=NBVE(I1,1)
              J2=NBVE(I1,NTVE(I1))
              IF((I.EQ.J1).OR.(I.EQ.J2)) THEN
                UNTMP=U(I,K)*COS(ANGLEQ(J))+V(I,K)*SIN(ANGLEQ(J))
                IF(UNTMP<0.0_SP) THEN
                  U_N(I,K) = UNTMP*COS(ANGLEQ(J))
                  V_N(I,K) = UNTMP*SIN(ANGLEQ(J))
                ENDIF
                GOTO 31
              ENDIF
            ENDDO
          ELSE IF(RIVER_INFLOW_LOCATION .EQ. 'edge') THEN
            DO J=1,NUMQBC
              J1=ICELLQ(J)
              IF(I.EQ.J1) THEN
                UNTMP=U(I,K)*COS(ANGLEQ(J))+V(I,K)*SIN(ANGLEQ(J))
                IF(UNTMP<0.0_SP) THEN
                  U_N(I,K) = UNTMP*COS(ANGLEQ(J))
                  V_N(I,K) = UNTMP*SIN(ANGLEQ(J))
                ENDIF
                GOTO 31
              ENDIF
            ENDDO
          ELSE
             PRINT*, 'inflow_type not correct'
             CALL PSTOP
          ENDIF
        ENDIF

#         if defined (THIN_DAM)
          IF(NBE_DAM(I) == 0.OR.K > KDAM1(I))THEN
!          IF(NBE_DAM(I) == 0)THEN
#         endif

        UNTMP=U(I,K)*COS(ALPHA1)+V(I,K)*SIN(ALPHA1)
        IF(UNTMP>0.0_SP) THEN
          U_N(I,K) = UNTMP*COS(ALPHA1)
          V_N(I,K) = UNTMP*SIN(ALPHA1)
        ENDIF

#         if defined (THIN_DAM)
          end if
#         endif

31    CONTINUE
      ENDIF
    ENDDO
  ENDDO

  RETURN
  END SUBROUTINE UV3D_SBD

  SUBROUTINE NAME_LIST_INITIALIZE_SEMI
  USE CONTROL, ONLY : IPT

  IMPLICIT NONE

  !--Parameters in NameList NML_SEMI
  IFCETA    = 0.55
  BEDF      = 1.0
  KSTAGE_UV = 1
  KSTAGE_TE = 1
  KSTAGE_TS = 1
  MSTG      = 'slow'
    
  END SUBROUTINE NAME_LIST_INITIALIZE_SEMI

  SUBROUTINE NAME_LIST_PRINT_SEMI
  USE CONTROL, ONLY : IPT
  
  IMPLICIT NONE

  write(UNIT=IPT,NML=NML_SEMI)

  RETURN
  END SUBROUTINE NAME_LIST_PRINT_SEMI

# endif
END MODULE MOD_SEMI_IMPLICIT
