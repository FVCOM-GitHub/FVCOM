!--------------------------------------------------------------------|
!--------------------------------------------------------------------|
  MODULE MOD_ACTION_IM
# if defined (WAVE_CURRENT_INTERACTION) && !defined (EXPLICIT)
# if defined (MULTIPROCESSOR)
  IMPLICIT NONE

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

  PRIVATE
  PUBLIC :: ALGORITHM_CRANK_NICOLSON
  PUBLIC :: ALGORITHM_FCT
  PUBLIC :: ACTION_ALLO
  PUBLIC :: ACTION_DEALLO
  PUBLIC :: ADV_N
  REAL, PUBLIC, ALLOCATABLE  ::  RHS(:)

  REAL, ALLOCATABLE,PUBLIC :: N31(:,:)
  REAL, ALLOCATABLE,PUBLIC :: N32(:,:,:)            
   REAL :: DS  !,DZETA    

!============ END PETSc BLOCK ========================================================

  CONTAINS
!==================================================================================
!
!==========================================================================|
!
   SUBROUTINE ACTION_ALLO
   
   USE ALL_VARS, ONLY : MT
   USE SWCOMM3, ONLY : MDC,MSC
!   USE MOD_USGRID, ONLY : MDC,MSC
   IMPLICIT NONE
   
   ALLOCATE(N31(MDC,MSC))   ;  N31 = 0.0
   ALLOCATE(N32(MDC,MSC,0:MT));  N32 = 0.0
   
   RETURN
   END SUBROUTINE ACTION_ALLO
!
!==========================================================================|
!
!==========================================================================|
!
   SUBROUTINE ACTION_DEALLO
   
   IMPLICIT NONE
   
   DEALLOCATE(N31)
   DEALLOCATE(N32)
   
   RETURN
   END SUBROUTINE ACTION_DEALLO
!==============================================================================
!
!==============================================================================
   SUBROUTINE ALGORITHM_FCT(CAS,IG,DTW,IDCMIN,IDCMAX)

   USE SWCOMM3, ONLY : MDC,MSC
   USE M_GENARR, ONLY : SPCSIG
!   USE MOD_USGRID, ONLY : SPCSIG,MDC,MSC
   USE ALL_VARS, ONLY : MT
   USE VARS_WAVE, ONLY :AC2
# if defined(MULTIPROCESSOR)
  USE MOD_PAR
# endif
   IMPLICIT NONE

   INTEGER :: ISS,ID,IG
   REAL :: CASR,CASL,FLUX1,FLUX2,FLUXLP,FLUXLM,FLUXHP,FLUXHM
   REAL :: MIN11,MIN22,MIN33,ADLP,ADLM
   REAL, DIMENSION(MDC,MSC) :: ADP,ADM,NL
   REAL :: CAS(MDC,MSC,10)
   REAL :: DTW
   
   INTEGER, DIMENSION(MSC) :: IDCMIN,IDCMAX
   
   N31 = 0.0

   DO ISS = 1,MSC
     IF(ISS == 1)THEN
       DS   = SPCSIG(ISS+1) - SPCSIG(ISS)                              
       DO ID = 1,MDC
         CASR = 0.5*(CAS(ID,ISS,1)+CAS(ID,ISS+1,1))
         FLUX1 = 0.5*(CASR+ABS(CASR))*AC2(ID,ISS,IG)
         FLUX2 = 0.5*(CASR-ABS(CASR))*AC2(ID,ISS+1,IG)
         FLUXLP = FLUX1+FLUX2

         FLUXLM = 0.0
         
         FLUXHP = CASR*0.5*(AC2(ID,ISS,IG)+AC2(ID,ISS+1,IG))
         FLUXHM = 0.0
   
         NL(ID,ISS) = AC2(ID,ISS,IG)-(FLUXLP-FLUXLM)*DTW/DS

         ADP(ID,ISS) = FLUXHP-FLUXLP
         ADM(ID,ISS) = FLUXHM-FLUXLM
       END DO	 
     ELSE IF(ISS == MSC)THEN
       DS   = SPCSIG(ISS) - SPCSIG(ISS-1)                              
       DO ID = 1,MDC
         CASR = CAS(ID,ISS,1)
         FLUX1 = 0.5*(CASR+ABS(CASR))*AC2(ID,ISS,IG)
         FLUX2 = 0.0
         FLUXLP = FLUX1+FLUX2

         CASL = CAS(ID,ISS-1,1)
	 FLUX1 = 0.5*(CASL+ABS(CASL))*AC2(ID,ISS-1,IG)
         FLUX2 = 0.5*(CASL-ABS(CASL))*AC2(ID,ISS,IG)
         FLUXLM = FLUX1+FLUX2
         
         FLUXHP = CASR*AC2(ID,ISS,IG)
         FLUXHM = CASL*AC2(ID,ISS-1,IG)
   
         NL(ID,ISS) = AC2(ID,ISS,IG)-(FLUXLP-FLUXLM)*DTW/DS

         ADP(ID,ISS) = FLUXHP-FLUXLP
         ADM(ID,ISS) = FLUXHM-FLUXLM
       END DO	 
     ELSE
       DS   = SPCSIG(ISS) - SPCSIG(ISS-1)                              
       DO ID = 1,MDC
         CASR = 0.5*(CAS(ID,ISS,1)+CAS(ID,ISS+1,1))
         FLUX1 = 0.5*(CASR+ABS(CASR))*AC2(ID,ISS,IG)
         FLUX2 = 0.5*(CASR-ABS(CASR))*AC2(ID,ISS+1,IG)
         FLUXLP = FLUX1+FLUX2

         CASL = 0.5*(CAS(ID,ISS,1)+CAS(ID,ISS-1,1))
	 FLUX1 = 0.5*(CASL+ABS(CASL))*AC2(ID,ISS-1,IG)
         FLUX2 = 0.5*(CASL-ABS(CASL))*AC2(ID,ISS,IG)
         FLUXLM = FLUX1+FLUX2
         
         FLUXHP = CASR*0.5*(AC2(ID,ISS,IG)+AC2(ID,ISS+1,IG))
         FLUXHM = CASL*0.5*(AC2(ID,ISS,IG)+AC2(ID,ISS-1,IG))
   
         NL(ID,ISS) = AC2(ID,ISS,IG)-(FLUXLP-FLUXLM)*DTW/DS

         ADP(ID,ISS) = FLUXHP-FLUXLP
         ADM(ID,ISS) = FLUXHM-FLUXLM
       END DO	 
     END IF
   END DO
   
   DO ISS = 1,MSC
     IF(ISS == 1)THEN
       DS   = SPCSIG(ISS+1) - SPCSIG(ISS)                              
       DO ID = 1,MDC
         MIN11 = ABS(ADP(ID,ISS))
         MIN22 = SIGN(1.,ADP(ID,ISS))*(NL(ID,ISS+2)-NL(ID,ISS+1))*DS/DTW
         ADLP = MIN(MIN11,MIN22)
         ADLP = MAX(0.,ADLP)
         ADLP = SIGN(1.,ADP(ID,ISS))*ADLP

         MIN11 = ABS(ADM(ID,ISS))
         MIN22 = SIGN(1.,ADM(ID,ISS))*(NL(ID,ISS+1)-NL(ID,ISS))*DS/DTW
         ADLM = MIN(MIN11,MIN22)
         ADLM = MAX(0.,ADLM)
         ADLM = SIGN(1.,ADM(ID,ISS))*ADLM

         N31(ID,ISS) = NL(ID,ISS)-(ADLP-ADLM)*DTW/DS
       END DO
     ELSE IF(ISS == 2)THEN   
       DS   = SPCSIG(ISS) - SPCSIG(ISS-1)                              
       DO ID = 1,MDC
         MIN11 = ABS(ADP(ID,ISS))
         MIN22 = SIGN(1.,ADP(ID,ISS))*(NL(ID,ISS+2)-NL(ID,ISS+1))*DS/DTW
         MIN33 = SIGN(1.,ADP(ID,ISS))*(NL(ID,ISS)-NL(ID,ISS-1))*DS/DTW
         ADLP = MIN(MIN11,MIN22,MIN33)
         ADLP = MAX(0.,ADLP)
         ADLP = SIGN(1.,ADP(ID,ISS))*ADLP

         MIN11 = ABS(ADM(ID,ISS))
         MIN22 = SIGN(1.,ADM(ID,ISS))*(NL(ID,ISS+1)-NL(ID,ISS))*DS/DTW
!         MIN33 = SIGN(1.,ADM(ID,ISS))*(NL(ID,ISS-1)-NL(ID,ISS-2))*DS/DTW
!         ADLM = MIN(MIN11,MIN22,MIN33)
         ADLM = MIN(MIN11,MIN22)
         ADLM = MAX(0.,ADLM)
         ADLM = SIGN(1.,ADM(ID,ISS))*ADLM
    
         N31(ID,ISS) = NL(ID,ISS)-(ADLP-ADLM)*DTW/DS
       END DO
     ELSE IF(ISS == MSC-1)THEN   
       DS   = SPCSIG(ISS) - SPCSIG(ISS-1)                              
       DO ID = 1,MDC
         MIN11 = ABS(ADP(ID,ISS))
!         MIN22 = SIGN(1.,ADP(ID,ISS))*(NL(ID,ISS+2)-NL(ID,ISS+1))*DS/DTW
         MIN33 = SIGN(1.,ADP(ID,ISS))*(NL(ID,ISS)-NL(ID,ISS-1))*DS/DTW
!         ADLP = MIN(MIN11,MIN22,MIN33)
         ADLP = MIN(MIN11,MIN33)
         ADLP = MAX(0.,ADLP)
         ADLP = SIGN(1.,ADP(ID,ISS))*ADLP

         MIN11 = ABS(ADM(ID,ISS))
         MIN22 = SIGN(1.,ADM(ID,ISS))*(NL(ID,ISS+1)-NL(ID,ISS))*DS/DTW
         MIN33 = SIGN(1.,ADM(ID,ISS))*(NL(ID,ISS-1)-NL(ID,ISS-2))*DS/DTW
         ADLM = MIN(MIN11,MIN22,MIN33)
         ADLM = MAX(0.,ADLM)
         ADLM = SIGN(1.,ADM(ID,ISS))*ADLM
    
         N31(ID,ISS) = NL(ID,ISS)-(ADLP-ADLM)*DTW/DS
       END DO
     ELSE IF(ISS == MSC)THEN
       DS   = SPCSIG(ISS) - SPCSIG(ISS-1)                              
       DO ID = 1,MDC
         MIN11 = ABS(ADP(ID,ISS))
         MIN33 = SIGN(1.,ADP(ID,ISS))*(NL(ID,ISS)-NL(ID,ISS-1))*DS/DTW
         ADLP = MIN(MIN11,MIN33)
         ADLP = MAX(0.,ADLP)
         ADLP = SIGN(1.,ADP(ID,ISS))*ADLP

         MIN11 = ABS(ADM(ID,ISS))
         MIN33 = SIGN(1.,ADM(ID,ISS))*(NL(ID,ISS-1)-NL(ID,ISS-2))*DS/DTW
         ADLM = MIN(MIN11,MIN33)
         ADLM = MAX(0.,ADLM)
         ADLM = SIGN(1.,ADM(ID,ISS))*ADLM
    
         N31(ID,ISS) = NL(ID,ISS)-(ADLP-ADLM)*DTW/DS
       END DO
     ELSE    
       DS   = SPCSIG(ISS) - SPCSIG(ISS-1)                              
       DO ID = 1,MDC
         MIN11 = ABS(ADP(ID,ISS))
         MIN22 = SIGN(1.,ADP(ID,ISS))*(NL(ID,ISS+2)-NL(ID,ISS+1))*DS/DTW
         MIN33 = SIGN(1.,ADP(ID,ISS))*(NL(ID,ISS)-NL(ID,ISS-1))*DS/DTW
         ADLP = MIN(MIN11,MIN22,MIN33)
         ADLP = MAX(0.,ADLP)
         ADLP = SIGN(1.,ADP(ID,ISS))*ADLP

         MIN11 = ABS(ADM(ID,ISS))
         MIN22 = SIGN(1.,ADM(ID,ISS))*(NL(ID,ISS+1)-NL(ID,ISS))*DS/DTW
         MIN33 = SIGN(1.,ADM(ID,ISS))*(NL(ID,ISS-1)-NL(ID,ISS-2))*DS/DTW
         ADLM = MIN(MIN11,MIN22,MIN33)
         ADLM = MAX(0.,ADLM)
         ADLM = SIGN(1.,ADM(ID,ISS))*ADLM
    
         N31(ID,ISS) = NL(ID,ISS)-(ADLP-ADLM)*DTW/DS
       END DO
     END IF  
   END DO    

   RETURN
   END SUBROUTINE ALGORITHM_FCT
!============================================================================|
!============================================================================|
   SUBROUTINE ALGORITHM_CRANK_NICOLSON(CAD,IG,DTW,IDCMIN,IDCMAX,DD)

   USE SWCOMM3, ONLY : MDC,MSC
!   USE MOD_USGRID, ONLY : MDC,MSC
   IMPLICIT NONE
   INTEGER :: ISS,ID,IDM1,IDM2,MDCM,IG,II,IDP1
   INTEGER :: IDDUM
   INTEGER, DIMENSION(MSC) :: IDCMIN,IDCMAX
   REAL, PARAMETER :: ZETA = 0.5
   REAL :: CAD(:,:,:)
   REAL :: DTW,DD
   REAL :: N32M,N32P
   
   REAL,DIMENSION(MDC) :: A,B,C,R,U
   
!   IF(ALLOCATED(N32)) DEALLOCATE(N32)
!   ALLOCATE(N32(MDC,MSC,MCGRD)); N32 = 0.0

   DO ISS = 1,MSC
     DO IDDUM = IDCMIN(ISS),IDCMAX(ISS)
       ID = MOD(IDDUM-1+MDC,MDC)+1
       IDP1 = MOD(IDDUM+MDC,MDC)+1
       IDM1 = MOD(IDDUM-2+MDC,MDC)+1
 
       B(ID) = 1.0
       IF(ID == 1)THEN
         A(ID) = 0.0
       ELSE
         A(ID) = -0.5*ZETA*DTW*CAD(IDM1,ISS,1)/DD
       END IF
       IF(ID == MDC)THEN
         C(ID) = 0.0
       ELSE
         C(ID) = 0.5*ZETA*DTW*CAD(IDP1,ISS,1)/DD
       END IF
 
       IF(ID == 1)THEN
         R(ID) = CAD(IDP1,ISS,1)*N31(IDP1,ISS) 
         R(ID) = (1.0-ZETA)*0.5*DTW*R(ID)/DD
         R(ID) = N31(ID,ISS)-R(ID)
       ELSE IF(ID == MDC)THEN
         R(ID) = -CAD(IDM1,ISS,1)*N31(IDM1,ISS)
         R(ID) = (1.0-ZETA)*0.5*DTW*R(ID)/DD
         R(ID) = N31(ID,ISS)-R(ID)       
       ELSE
         R(ID) = CAD(IDP1,ISS,1)*N31(IDP1,ISS)-CAD(IDM1,ISS,1)*N31(IDM1,ISS) 
         R(ID) = (1.0-ZETA)*0.5*DTW*R(ID)/DD
         R(ID) = N31(ID,ISS)-R(ID)     
       END IF
     END DO
       
     CALL TRIDAG(A,B,C,R,U,MDC)
       
     DO IDDUM = IDCMIN(ISS),IDCMAX(ISS)
       ID = MOD(IDDUM-1+MDC,MDC)+1
       N32(ID,ISS,IG) = U(ID)
     END DO	 
   END DO
   
   RETURN
   END SUBROUTINE ALGORITHM_CRANK_NICOLSON
!==========================================================================|
!
!==================================================================================!
    SUBROUTINE TRIDAG(A,B,C,R,U,N)
    IMPLICIT NONE
    INTEGER  :: N,J
    REAL,DIMENSION(N) :: A,B,C,R,U
    INTEGER, PARAMETER :: NMAX = 500
    REAL BET,GAM(NMAX)
    
    IF(B(1) == 0.)PAUSE 'TRIDAG: REWRITE EQUATIONS'
    BET = B(1)
    U(1) = R(1)/BET
    DO J=2,N
      GAM(J) = C(J-1)/BET
      BET = B(J)-A(J)*GAM(J)
      IF(BET == 0.)PAUSE 'TRIDAG FAILED'
      U(J) = (R(J)-A(J)*U(J-1))/BET
    END DO
    DO J=N-1,1,-1
      U(J) = U(J)-GAM(J+1)*U(J+1)
    END DO
    
    RETURN
    END SUBROUTINE TRIDAG  
!==========================================================================|
!
!==================================================================================!
  SUBROUTINE ADV_N(DTW)

  USE VARS_WAVE
# if defined(MULTIPROCESSOR)
  USE MOD_PAR
# endif
   USE SWCOMM3
   USE M_GENARR
   USE MOD_OBCS
!  USE MOD_USGRID
#  if defined (SPHERICAL)   
   USE MOD_SPHERICAL
#  endif   

# if defined(PLBC)
  USE MOD_PERIODIC_LBC
# endif

  USE MOD_PETSC, ONLY : ALO_2_PLO_NODE,BL_WAVE,B_WAVE,N_VERTS,L2G_WAVE,A_WAVE,  &
                        PLO_2_PGO_NODE,X_WAVE,XL_WAVE,G2L_WAVE,XVALS_WAVE,      &
                        PLO_2_ALO_NODE,PETSc_SOLVER_WAVE
  
  IMPLICIT NONE

  INTEGER  I, J, IK, K, I1, IA, IB, TMP2, TMP3, ISS, ID
  INTEGER  NNZ2, LN, LOC
  INTEGER  PETSc_POS, NODE, PROW1, PROW2, PCOL1, PCOL2, PCOL3
  REAL DIJ, UIJ, VIJ, EXFLUX, DTW
  REAL, ALLOCATABLE ::  UF_AVG(:), VF_AVG(:)
  REAL, ALLOCATABLE ::  XFLUX(:)
!  REAL :: DEP2(MT),AC2LOC(0:MT)
!  REAL, ALLOCATABLE :: DEP2(:),AC2LOC(:)   
  REAL, ALLOCATABLE :: DEP2(:)   
  REAL(SP), ALLOCATABLE :: AC2LOC(:)   
  REAL :: CANX,CANY,UN
  REAL    :: DEPLOC,KWAVELOC,CGLOC,NN,ND,SPCSIGL
  REAL :: UTMP,VTMP,DLTXE_TMP,DLTYE_TMP

   REAL, DIMENSION(0:M) :: PSPX,PSPY
   REAL :: FF1,XI,YI
   REAL :: DXA,DYA,DXB,DYB,FIJ1,FIJ2
   INTEGER :: I2
#  if defined (SPHERICAL)
!!$   REAL(8) :: TY,TXPI,TYPI
   REAL(8) :: XTMP1,XTMP
!!$   REAL(8) :: X1_DP,Y1_DP,X2_DP,Y2_DP,XII,YII
!!$   REAL(8) :: X11_TMP,Y11_TMP,X33_TMP,Y33_TMP
#  endif
   REAL(SP) :: UA_NODE,VA_NODE
   INTEGER  :: CNT,JJ
 
  PetscReal :: STERM
  PetscReal :: VCOL1,VCOL2
  PetscInt:: IERR
  
  ALLOCATE( DEP2(MT),AC2LOC(0:MT)  )
  DEP2(1:MT) = COMPDA(1:MT,JDP2)

# if defined(MULTIPROCESSOR)
!  IF(PAR) CALL EXCHANGE(NC,MT,1,MYID,NPROCS,DEP2)
  IF(PAR) CALL AEXCHANGE(NC,MYID,NPROCS,DEP2)     
# endif  

  DO ISS = 1,MSC
    DO ID = 1,MDC

     PSPX  = 0.0 
     PSPY  = 0.0 
!periodic 
# if defined(PLBC)
       CALL replace_N32(N32,ID,ISS)
# endif

     DO I=1,M
       DO J=1,NTSN(I)-1
         I1=NBSN(I,J)
         I2=NBSN(I,J+1)
!!$         IF(DEP2(I1) <= DEPMIN .AND. DEP2(I2) > DEPMIN)THEN
!!$          FF1=0.5*(N32(ID,ISS,I)+N32(ID,ISS,I2))
!!$         ELSE IF(DEP2(I1) > DEPMIN .AND. DEP2(I2) <= DEPMIN)THEN
!!$          FF1=0.5*(N32(ID,ISS,I1)+N32(ID,ISS,I))
!!$         ELSE IF(DEP2(I1) <= DEPMIN .AND. DEP2(I2) <= DEPMIN)THEN
!!$          FF1=N32(ID,ISS,I)
!!$         ELSE
!!$          FF1=0.5*(N32(ID,ISS,I1)+N32(ID,ISS,I2))
!!$         END IF
         FF1=0.5*(N32(ID,ISS,I1)+N32(ID,ISS,I2))
!!$#        if defined (SPHERICAL)
!!$         XTMP  = VX(I2)*TPI-VX(I1)*TPI
!!$	 XTMP1 = VX(I2)-VX(I1)
!!$	 IF(XTMP1 >  180.0)THEN
!!$	   XTMP = -360.0*TPI+XTMP
!!$	 ELSE IF(XTMP1 < -180.0)THEN
!!$	   XTMP =  360.0*TPI+XTMP
!!$	 END IF  
!!$         TXPI=XTMP*COS(DEG2RAD*VY(I))
!!$         TYPI=(VY(I1)-VY(I2))*TPI
!!$         PSPX(I)=PSPX(I)+FF1*TYPI
!!$         PSPY(I)=PSPY(I)+FF1*TXPI
!!$!         PSPXD(I)=PSPXD(I)+FFD*TYPI
!!$!         PSPYD(I)=PSPYD(I)+FFD*TXPI
!!$#        else
!!$         PSPX(I)=PSPX(I)+FF1*(VY(I1)-VY(I2))
!!$         PSPY(I)=PSPY(I)+FF1*(VX(I2)-VX(I1))
!!$#        endif
          PSPX(I)=PSPX(I)+FF1*DLTYTRIE(I,J)
          PSPY(I)=PSPY(I)+FF1*DLTXTRIE(I,J)
       END DO
       PSPX(I)=PSPX(I)/ART2(I)
       PSPY(I)=PSPY(I)/ART2(I)
     END DO

    CANX = 0.0
    CANY = 0.0

    CALL VecSet(BL_WAVE,ZERO,IERR);CHKERRQ(IERR)
    CALL VecSet(B_WAVE,ZERO,IERR);CHKERRQ(IERR)

  DO I=1, M
    PETSc_POS = ALO_2_PLO_NODE(I)
    IF (PETSc_POS > N_VERTS) CYCLE

    IF(ISONB_W(I) /= 2) THEN
      STERM = N32(ID,ISS,I)       
    ELSE
      STERM = AC2(ID,ISS,I)  
    ENDIF
    CALL VecSetValues(BL_WAVE,1,PETSc_POS-1,STERM,INSERT_VALUES,IERR);CHKERRQ(IERR)
  ENDDO

! PETSC changed the VecScatterBegin/VecScatterEnd interfaces after 2.3.2
# if defined (OLD_PETSC)
  CALL VecScatterBegin(BL_WAVE,B_WAVE,INSERT_VALUES,SCATTER_FORWARD,L2G_WAVE,IERR);CHKERRQ(IERR)
  CALL VecScatterEnd(BL_WAVE,B_WAVE,INSERT_VALUES,SCATTER_FORWARD,L2G_WAVE,IERR);CHKERRQ(IERR)
# else
  CALL VecScatterBegin(L2G_WAVE,BL_WAVE,B_WAVE,INSERT_VALUES,SCATTER_FORWARD,IERR);CHKERRQ(IERR) 
  CALL VecScatterEnd(L2G_WAVE,BL_WAVE,B_WAVE,INSERT_VALUES,SCATTER_FORWARD,IERR);CHKERRQ(IERR)   
# endif

  CALL MatZeroEntries(A_WAVE,IERR);CHKERRQ(IERR)

  DO I=1, NCV_I
    I1 = NTRG(I)
    IA = NIEC(I,1)
    IB = NIEC(I,2)

!!$    XI=0.5*(XIJE(I,1)+XIJE(I,2))
!!$    YI=0.5*(YIJE(I,1)+YIJE(I,2))
!!$#   if defined (SPHERICAL)
!!$    X1_DP=XIJE(I,1)
!!$    Y1_DP=YIJE(I,1)
!!$    X2_DP=XIJE(I,2)
!!$    Y2_DP=YIJE(I,2)
!!$!    CALL ARCC(X2_DP,Y2_DP,X1_DP,Y1_DP,XII,YII)
!!$    XII = XCG2(I)
!!$    YII = YCG2(I)
!!$    XI=XII		
!!$    XTMP  = XI*TPI-VX(IA)*TPI
!!$    XTMP1 = XI-VX(IA)
!!$    IF(XTMP1 >  180.0)THEN
!!$      XTMP = -360.0*TPI+XTMP
!!$    ELSE IF(XTMP1 < -180.0)THEN
!!$      XTMP =  360.0*TPI+XTMP
!!$    END IF	 

!!$!   DXA=XTMP*COS(DEG2RAD*VY(IA))    
!!$    DXA=XTMP*VAL_COS_VY(IA)    
!!$    DYA=(YI-VY(IA))*TPI
!!$    XTMP  = XI*TPI-VX(IB)*TPI
!!$    XTMP1 = XI-VX(IB)
!!$    IF(XTMP1 >  180.0)THEN
!!$      XTMP = -360.0*TPI+XTMP
!!$    ELSE IF(XTMP1 < -180.0)THEN
!!$      XTMP =  360.0*TPI+XTMP
!!$    END IF	 

!!$!    DXB=XTMP*COS(DEG2RAD*VY(IB)) 
!!$    DXB=XTMP*VAL_COS_VY(IB) 
!!$    DYB=(YI-VY(IB))*TPI
!!$#   else
!!$    DXA=XI-VX(IA)
!!$    DYA=YI-VY(IA)
!!$    DXB=XI-VX(IB)
!!$    DYB=YI-VY(IB)
!!$#   endif

!!$    FIJ1=DXA*PSPX(IA)+DYA*PSPY(IA)
!!$    FIJ2=DXB*PSPX(IB)+DYB*PSPY(IB)
    FIJ1=DLTXNCVE(I,1)*PSPX(IA)+DLTYNCVE(I,1)*PSPY(IA)
    FIJ2=DLTXNCVE(I,2)*PSPX(IB)+DLTYNCVE(I,2)*PSPY(IB)

    CALL SWAPAR1(I1,ISS,ID,DEP2(1),KWAVELOC,CGLOC)
	 
#   if !defined (TWO_D_MODEL) 
    UTMP = (COMPDA(NV(I1,1),JVX2)+COMPDA(NV(I1,2),JVX2)+        &
	    COMPDA(NV(I1,3),JVX2))/3.0
    VTMP = (COMPDA(NV(I1,1),JVY2)+COMPDA(NV(I1,2),JVY2)+        &
	    COMPDA(NV(I1,3),JVY2))/3.0
#   else
    UTMP = UA(I1)
    VTMP = VA(I1)
#   endif    
	  
    CALL SPROXY(I1     ,ISS         ,ID          ,CANX   ,CANY   , &
	        CGLOC  ,SPCDIR(ID,2),SPCDIR(ID,3),UTMP   ,VTMP   )

!    DIJ = DEP2(NV(I1,1))+DEP2(NV(I1,2))+DEP2(NV(I1,3))
!    DIJ = DIJ/3.0
    
    IF((ISONB_W(IA)+ISONB_W(IB)) < 4) THEN
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
      
      VCOL1 = DTW*(-CANY*DLTXE(I)+CANX*DLTYE(I))/ART1(IA)
      VCOL2 = -DTW*(-CANY*DLTXE(I)+CANX*DLTYE(I))/ART1(IB)

      IF(VCOL1 > 0.0)THEN
        CALL MatSetValuesLocal(A_WAVE,1,PROW1-1,1,PCOL1-1,VCOL1,ADD_VALUES,IERR);CHKERRQ(IERR)
        CALL MatSetValuesLocal(A_WAVE,1,PROW2-1,1,PCOL1-1,VCOL2,ADD_VALUES,IERR);CHKERRQ(IERR)
        
        CALL VecSetValues(B_WAVE,1,PLO_2_PGO_NODE(PROW1)-1,-VCOL1*FIJ1,ADD_VALUES,IERR);CHKERRQ(IERR)
        CALL VecSetValues(B_WAVE,1,PLO_2_PGO_NODE(PROW2)-1,-VCOL2*FIJ1,ADD_VALUES,IERR);CHKERRQ(IERR)
     ELSE
        CALL MatSetValuesLocal(A_WAVE,1,PROW1-1,1,PCOL2-1,VCOL1,ADD_VALUES,IERR);CHKERRQ(IERR)
        CALL MatSetValuesLocal(A_WAVE,1,PROW2-1,1,PCOL2-1,VCOL2,ADD_VALUES,IERR);CHKERRQ(IERR)

        CALL VecSetValues(B_WAVE,1,PLO_2_PGO_NODE(PROW1)-1,-VCOL1*FIJ2,ADD_VALUES,IERR);CHKERRQ(IERR)
        CALL VecSetValues(B_WAVE,1,PLO_2_PGO_NODE(PROW2)-1,-VCOL2*FIJ2,ADD_VALUES,IERR);CHKERRQ(IERR)
     END IF
    ENDIF
  ENDDO

  CALL VecAssemblyBegin(B_WAVE,IERR);CHKERRQ(IERR)
  CALL VecAssemblyEnd(B_WAVE,IERR);CHKERRQ(IERR)

  DO I=1, M
    PETSc_POS = ALO_2_PLO_NODE(I)
    IF (PETSc_POS > N_VERTS) CYCLE

    IF(ISONB_W(I) == 1)THEN
      DEPLOC = DEP2(I)
      IF(DEPLOC <= DEPMIN)THEN
!      IF(DEPLOC < 0.05)THEN
!       *** depth is negative ***
        KWAVELOC = -1.                                           
        CGLOC   = 0.                                            
      ELSE
!       *** call KSCIP1 to compute KWAVE and CGO ***
        SPCSIGL = SPCSIG(ISS)
        CALL KSCIP1(1,SPCSIGL,DEPLOC,KWAVELOC,CGLOC,NN,ND)                                 
      ENDIF
      CANX = CGLOC * SPCDIR(ID,2)
      CANY = CGLOC * SPCDIR(ID,3)
!
!     --- adapt the velocities in case of diffraction
!
      IF(IDIFFR == 1 .AND. PDIFFR(3) /= 0.)THEN 
!JQI       CANX = CNAX*DIFPARAM(I)      
!JQI       CANY = CNAY*DIFPARAM(I)      
      END IF
!
!     --- ambient currents added
!
      IF(ICUR == 1)THEN 
#     if !defined (TWO_D_MODEL)
        CANX = CANX + COMPDA(I,JVX2)
        CANY = CANY + COMPDA(I,JVY2)
#     else
        UA_NODE = 0.0_SP
	VA_NODE = 0.0_SP
	CNT = 0
        DO JJ=1,NTVE(I)
         CNT =CNT + 1
         UA_NODE = UA_NODE + UA(NBVE(I,JJ))
         VA_NODE = VA_NODE + VA(NBVE(I,JJ))
        ENDDO
        UA_NODE = UA_NODE/CNT
        VA_NODE = VA_NODE/CNT
        CANX = CANX + UA_NODE
        CANY = CANY + VA_NODE
#     endif	    
      END IF
	    
      PROW1 = ALO_2_PLO_NODE(I)

      PCOL1 = ALO_2_PLO_NODE(I)
      
!      IF(NBSN(I,NTSN(I)-1) > M)THEN
!#      if defined (SPHERICAL)
!        XTMP1 = VX(I)-VX(NBSN(I,2))
!        XTMP = XTMP1*TPI
!	IF(XTMP1 > 180.0)THEN
!	  XTMP = -360.0*TPI+XTMP
!	ELSE IF(XTMP1 < -180.0)THEN
!	  XTMP =  360.0*TPI+XTMP
!	END IF
!	DLTXE_TMP = XTMP*COS(DEG2RAD*VY(I))
!        DLTYE_TMP = (VY(I)-VY(NBSN(I,2)))*TPI
!#      else
!        DLTXE_TMP = VX(I)-VX(NBSN(I,2))
!        DLTYE_TMP = VY(I)-VY(NBSN(I,2))
!#      endif
!      ELSE IF(NBSN(I,2) > M)THEN
!#      if defined (SPHERICAL)
!        XTMP1 = VX(NBSN(I,NTSN(I)-1))-VX(I)
!        XTMP = XTMP1*TPI
!	IF(XTMP1 > 180.0)THEN
!	  XTMP = -360.0*TPI+XTMP
!	ELSE IF(XTMP1 < -180.0)THEN
!	  XTMP =  360.0*TPI+XTMP
!	END IF
!	DLTXE_TMP = XTMP*COS(DEG2RAD*VY(I))
!        DLTYE_TMP = (VY(NBSN(I,NTSN(I)-1))-VY(I))*TPI
!#      else
!        DLTXE_TMP = VX(NBSN(I,NTSN(I)-1))-VX(I)
!        DLTYE_TMP = VY(NBSN(I,NTSN(I)-1))-VY(I)
!#      endif
!      ELSE 
#      if defined (SPHERICAL)
        XTMP1 = VX(NBSN(I,NTSN(I)-1))-VX(NBSN(I,2))
        XTMP = XTMP1*TPI
	IF(XTMP1 > 180.0)THEN
	  XTMP = -360.0*TPI+XTMP
	ELSE IF(XTMP1 < -180.0)THEN
	  XTMP =  360.0*TPI+XTMP
	END IF
	DLTXE_TMP = XTMP*COS(DEG2RAD*VY(I))
        DLTYE_TMP = (VY(NBSN(I,NTSN(I)-1))-VY(NBSN(I,2)))*TPI
#      else
        DLTXE_TMP = VX(NBSN(I,NTSN(I)-1))-VX(NBSN(I,2))
        DLTYE_TMP = VY(NBSN(I,NTSN(I)-1))-VY(NBSN(I,2))
#      endif
!      END IF
      
!      VCOL1 = DTW*(-CANY*DLTXE_TMP+CANX*DLTYE_TMP)/ART1(I)
      VCOL1 = -DTW*(CANY*DLTXE_TMP-CANX*DLTYE_TMP)/ART1(I)
!      VCOL1 = MAX(0.0,-VCOL1)
      VCOL1 = MAX(0.0,0.5*VCOL1)
      
      CALL MatSetValuesLocal(A_WAVE,1,PROW1-1,1,PCOL1-1,VCOL1,ADD_VALUES,IERR);CHKERRQ(IERR)

      VCOL1 = 0.0D0
      DO J=1, NTSN(I)-1
        PCOL1 = ALO_2_PLO_NODE(NBSN(I,J))
        CALL MatSetValuesLocal(A_WAVE,1,PROW1-1,1,PCOL1-1,VCOL1,ADD_VALUES,IERR);CHKERRQ(IERR)
      END DO
    END IF
  END DO

  DO I=1, M
   IF(ISONB_W(I) /= 2)THEN
    PROW1  = ALO_2_PLO_NODE(I)
    IF(PROW1 > N_VERTS) CYCLE
    PCOL1  = ALO_2_PLO_NODE(I)
    
    VCOL1 = 1.0D0
    CALL MatSetValuesLocal(A_WAVE,1,PROW1-1,1,PCOL1-1,VCOL1,ADD_VALUES,IERR);CHKERRQ(IERR)
   END IF
  ENDDO

  CALL MatAssemblyBegin(A_WAVE,MAT_FINAL_ASSEMBLY,IERR);CHKERRQ(IERR)
  CALL MatAssemblyEnd(A_WAVE,MAT_FINAL_ASSEMBLY,IERR);CHKERRQ(IERR)
  CALL MPI_BARRIER(MPI_FVCOM_GROUP,IERR)

  DO I=1, IOBCN_W

    PROW1 = ALO_2_PLO_NODE(I_OBC_N_W(I))
    IF(PROW1>N_VERTS) CYCLE

    VCOL1 = 0.0D0
    NODE  = I_OBC_N_W(I)
    DO J=1, NTSN(NODE)-1
      PCOL1 = ALO_2_PLO_NODE(NBSN(NODE,J))
      CALL MatSetValuesLocal(A_WAVE,1,PROW1-1,1,PCOL1-1,VCOL1,INSERT_VALUES,IERR);CHKERRQ(IERR)
    ENDDO

    VCOL1 = 1.0D0
    PROW1 = ALO_2_PLO_NODE(NODE)
    CALL MatSetValuesLocal(A_WAVE,1,PROW1-1,1,PROW1-1,VCOL1,INSERT_VALUES,IERR);CHKERRQ(IERR)

  ENDDO

  CALL MatAssemblyBegin(A_WAVE,MAT_FINAL_ASSEMBLY,IERR);CHKERRQ(IERR)
  CALL MatAssemblyEnd(A_WAVE,MAT_FINAL_ASSEMBLY,IERR);CHKERRQ(IERR)

  AC2LOC = 0.0
  DO I=1,M
   AC2LOC(I) = AC2(ID,ISS,I)
  END DO
   
  CALL PETSc_SETICS_WAVE(AC2LOC)
  CALL PETSc_SOLVER_WAVE

# if defined (OLD_PETSC)
  CALL VecScatterBegin(X_WAVE,XL_WAVE,INSERT_VALUES,SCATTER_FORWARD,G2L_WAVE,IERR);CHKERRQ(IERR)
  CALL VecScatterEnd(X_WAVE,XL_WAVE,INSERT_VALUES,SCATTER_FORWARD,G2L_WAVE,IERR);CHKERRQ(IERR)
# else
  CALL VecScatterBegin(G2L_WAVE,X_WAVE,XL_WAVE,INSERT_VALUES,SCATTER_FORWARD,IERR);CHKERRQ(IERR)
  CALL VecScatterEnd(G2L_WAVE,X_WAVE,XL_WAVE,INSERT_VALUES,SCATTER_FORWARD,IERR);CHKERRQ(IERR)
# endif

  CALL VecGetArrayF90(XL_WAVE,XVALS_WAVE,IERR);CHKERRQ(IERR)

  AC2(ID,ISS,:) = 0.0
  DO I=1,N_VERTS
    IK = PLO_2_ALO_NODE(I)
    AC2(ID,ISS,IK) = MAX(0.0, XVALS_WAVE(I))
!    AC2(ID,ISS,IK) = XVALS_WAVE(I)
!    IF(AC2(ID,ISS,IK) < 0.0)AC2(ID,ISS,IK) = 0.0
  ENDDO

  DO I=1,M
   IF(DEP2(I) <= DEPMIN)THEN
     AC2(ID,ISS,I) = 0.0_SP
   END IF
  END DO   

# if defined (MULTIPROCESSOR)
  IF(PAR)THEN
    AC2LOC = 0.0
    DO I=1,M          !MT
      AC2LOC(I) = AC2(ID,ISS,I)
    END DO
   
    CALL NODE_MATCH(0,NBN,BN_MLT,BN_LOC,BNC,MT,1,MYID,NPROCS,AC2LOC)
!    CALL EXCHANGE(NC,MT,1,MYID,NPROCS,AC2LOC)
    CALL AEXCHANGE(NC,MYID,NPROCS,AC2LOC)   
    AC2(ID,ISS,:) = 0.0
    DO I=1,MT
     AC2(ID,ISS,I) = AC2LOC(I)
    END DO
  END IF
   
# endif
# if defined(PLBC)
  CALL replace_ac2(ID,ISS)
# endif

  END DO
  END DO

  RETURN
  END SUBROUTINE ADV_N
!
!==========================================================================|
!
!==========================================================================|
!
  SUBROUTINE PETSc_SETICS_WAVE(ELL)
  USE MOD_PREC
  USE ALL_VARS, ONLY : MYID,MT
  USE MOD_PETSC, ONLY : USE_LAST,X_WAVE,XL_WAVE,L2G_WAVE,N_VERTS,PLO_2_ALO_NODE
  IMPLICIT NONE
  INTEGER :: I, IK
  PetscReal :: QTERM
  PetscInt  :: IERR
  PetscScalar :: ZERO   =  0.0D0
  CHARACTER(LEN=20) :: SUBNAME = 'PETSc_SETICS'
  REAL(SP) :: ELL(0:MT)

  IF(.NOT.USE_LAST)THEN
    CALL VecSet(X_WAVE,ZERO,IERR);CHKERRQ(IERR)
    RETURN
  END IF

  DO I=1,N_VERTS
    IK = PLO_2_ALO_NODE(I)
    QTERM = ELL(IK)
    CALL VecSetValues(XL_WAVE,1,I-1,QTERM,INSERT_VALUES,IERR);CHKERRQ(IERR)
  END DO

# if defined (OLD_PETSC)
  CALL VecScatterBegin(XL_WAVE,X_WAVE,INSERT_VALUES,SCATTER_FORWARD,L2G_WAVE,IERR);CHKERRQ(IERR)
  CALL VecScatterEnd(XL_WAVE,X_WAVE,INSERT_VALUES,SCATTER_FORWARD,L2G_WAVE,IERR);CHKERRQ(IERR)
# else
  CALL VecScatterBegin(L2G_WAVE,XL_WAVE,X_WAVE,INSERT_VALUES,SCATTER_FORWARD,IERR);CHKERRQ(IERR)
  CALL VecScatterEnd(L2G_WAVE,XL_WAVE,X_WAVE,INSERT_VALUES,SCATTER_FORWARD,IERR);CHKERRQ(IERR)
# endif


  RETURN
  END SUBROUTINE PETSc_SETICS_WAVE
!=========================================================================================
!
!=========================================================================================
# endif

# endif

END MODULE MOD_ACTION_IM
