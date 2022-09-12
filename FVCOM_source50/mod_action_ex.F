MODULE MOD_ACTION_EX
      
  IMPLICIT NONE
  SAVE
   
  REAL, ALLOCATABLE :: N31(:,:)
  REAL, ALLOCATABLE :: N32(:,:,:)            

  CONTAINS
!
!==============================================================================|
!
!==========================================================================|
!
  SUBROUTINE ACTION_ALLO
   
  USE ALL_VARS, ONLY : MT
  USE SWCOMM3, ONLY : MDC,MSC
!  USE MOD_USGRID, ONLY : MDC,MSC
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
!
!==========================================================================|
!
!==========================================================================|
!
  SUBROUTINE ALGORITHM_FCT(CAS,IG,DTW,IDCMIN,IDCMAX)

  USE SWCOMM3, ONLY : MDC,MSC
  USE M_GENARR, ONLY : SPCSIG
!  USE MOD_USGRID, ONLY : SPCSIG,MDC,MSC
  USE VARS_WAVE, ONLY : MT,AC2
!  USE ALL_VARS, ONLY : MT,AC2
# if defined(MULTIPROCESSOR)
  USE MOD_PAR
# endif
  IMPLICIT NONE

  INTEGER :: ISS,ID,IG
  REAL :: CASR,CASL,FLUX1,FLUX2,FLUXLP,FLUXLM,FLUXHP,FLUXHM
  REAL :: MIN11,MIN22,MIN33,ADLP,ADLM
  REAL, DIMENSION(MDC,MSC) :: ADP,ADM,NL
  REAL :: CAS(MDC,MSC,10)
  REAL :: DTW,DS
   
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
      DS = SPCSIG(ISS+1) - SPCSIG(ISS)                              
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
        ADLM = MIN(MIN11,MIN22)
        ADLM = MAX(0.,ADLM)
        ADLM = SIGN(1.,ADM(ID,ISS))*ADLM
    
        N31(ID,ISS) = NL(ID,ISS)-(ADLP-ADLM)*DTW/DS
      END DO
    ELSE IF(ISS == MSC-1)THEN   
      DS   = SPCSIG(ISS) - SPCSIG(ISS-1)                              
      DO ID = 1,MDC
        MIN11 = ABS(ADP(ID,ISS))
        MIN33 = SIGN(1.,ADP(ID,ISS))*(NL(ID,ISS)-NL(ID,ISS-1))*DS/DTW
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
!  USE MOD_USGRID, ONLY : MDC,MSC
  IMPLICIT NONE
  INTEGER :: ISS,ID,IDM1,IDM2,IDP1,MDCM,IG,II
  INTEGER :: IDDUM
  INTEGER, DIMENSION(MSC) :: IDCMIN,IDCMAX
  REAL, PARAMETER :: ZETA = 0.5
  REAL :: CAD(:,:,:)
  REAL :: DTW,DD
  REAL :: N32M,N32P
   
  REAL,DIMENSION(MDC) :: A,B,C,R,U
   
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
!==========================================================================|

  SUBROUTINE ADV_N(DTW)                   

!------------------------------------------------------------------------------|

  USE VARS_WAVE
!  USE ALL_VARS
# if defined(MULTIPROCESSOR)
  USE MOD_PAR
# endif
# if defined (SPHERICAL)   
  USE MOD_SPHERICAL
  USE MOD_NORTHPOLE    
# endif   
  USE SWCOMM3
  USE M_GENARR
# if defined(PLBC)
  USE MOD_PERIODIC_LBC
# endif
# if defined(ICE)
  USE ICE_STATE        
# endif

  IMPLICIT NONE
  
  REAL(SP),ALLOCATABLE :: EXFLUX(:,:),UN(:,:),PSPX(:,:,:),PSPY(:,:,:),FF1(:,:),FIJ1(:,:,:)
  REAL(SP),ALLOCATABLE :: XFLUX(:,:,:)
  REAL(SP),ALLOCATABLE :: XFLUX_TMP(:)
  REAL(SP) :: UTMP,VTMP,UN1,X11,Y11,X22,Y22,X33,Y33,TMP1,TMP2,XI,YI
  REAL(SP) :: DXA,DYA,DXB,DYB,FIJ1_TMP
  REAL(SP) :: TXX,TYY,FXX,FYY,VISCOF,TEMP,STPOINT
  REAL(SP) :: FACT,FM1
  INTEGER  :: I,I1,I2,IA,IB,J,J1,J2,K,JTMP,JJ,II
  REAL(SP) :: AC1MIN, AC1MAX, AC2MIN, AC2MAX
  INTEGER  :: ID,ISS,IG,IG2 
  REAL(SP) :: XIN,YIN,XIC,YIC,CANX,CANY
  REAL(SP) :: CAX(MDC,MSC),CAY(MDC,MSC)
!  REAL(SP) :: DTW,RF(MDC,MSC),DF(MDC,MSC)
  REAL(SP) :: RF(MDC,MSC),DF(MDC,MSC)
  REAL     :: DTW
  REAL(SP) :: SPCDIR2(MDC),SPCDIR3(MDC),CANX_2D(MDC,MSC),CANY_2D(MDC,MSC),CGO_1D(MSC),CGO_2D(MDC,MSC)
  INTEGER  :: NTSN_T
  REAL(SP) :: N32_T,XFLUX_T
  REAL(SP),ALLOCATABLE :: DEP2(:),AC2LOC(:)
  REAL(SP) :: DEPLOC,KWAVELOC,CGLOC,NN,ND,SPCSIGL
# if defined (ICE)
  REAL(SP),ALLOCATABLE :: Sice_TMP(:),SiceLOC(:)
# endif  
# if defined (SPHERICAL)
  REAL(8) :: XTMP1,XTMP
  REAL(8) :: DLTXE_TMP,DLTYE_TMP
# endif
  REAL(SP) :: UA_NODE,VA_NODE
  INTEGER  :: CNT

  REAL,ALLOCATABLE ::EXFLUXA(:,:),EXFLUXB(:,:),UNA(:,:),UNBB(:,:),FIJ2(:,:,:),&
                     FIJ1A(:,:,:),FIJ1B(:,:,:),FIJ2A(:,:,:),FIJ2B(:,:,:),     &
		     N32_TMPP2(:,:),N32_TMPP3(:,:)  
  REAL :: DLTXEA,DLTXEB,DLTXETMP,DLTYETMP,C_ICE
  REAL     :: CANXAB(MDC,MSC),CANYA(MDC,MSC),CANYB(MDC,MSC)
  REAL(SP)   :: XFLUX_1D(0:MT),PSPX_1D(M),PSPY_1D(M)
  REAL   :: RATIO(MSC)
  REAL   :: EXFLUX_TMPA(MDC,MSC),EXFLUX_TMPB(MDC,MSC)
  REAL   :: EXFLUX_TMP(MDC,MSC) 

!------------------------------------------------------------------------------!
  ALLOCATE(DEP2(MT),AC2LOC(0:MT),XFLUX_TMP(0:MT))
# if defined (ICE)
  ALLOCATE(Sice_TMP(0:MT),SiceLOC(0:MT))
  ALLOCATE(Sice(MDC,MSC,0:MT));Sice=0.0
# endif  
  ALLOCATE(EXFLUX(MDC,MSC),FF1(MDC,MSC),FIJ1(MDC,MSC,2),UN(MDC,MSC))
  ALLOCATE(XFLUX(MDC,MSC,0:MT));XFLUX=0.0
  ALLOCATE(PSPX(MDC,MSC,M));PSPX=0.0
  ALLOCATE(PSPY(MDC,MSC,M));PSPY=0.0
  IF(HIGH_LATITUDE_WAVE)THEN 
    ALLOCATE(EXFLUXA(MDC,MSC),EXFLUXB(MDC,MSC),FIJ2(MDC,MSC,0:M),&
             FIJ1A(MDC,MSC,0:M),FIJ1B(MDC,MSC,0:M),FIJ2A(MDC,MSC,0:M),FIJ2B(MDC,MSC,0:M),UNA(MDC,MSC),UNBB(MDC,MSC))
    ALLOCATE(N32_TMPP2(MDC,MSC),N32_TMPP3(MDC,MSC))
  END IF
  DEP2(1:MT) = COMPDA(1:MT,JDP2)
# if defined(MULTIPROCESSOR)
!  IF(PAR) CALL EXCHANGE(NC,MT,1,MYID,NPROCS,DEP2)
  IF(PAR)CALL aexchange(nc,myid,nprocs,dep2)
# endif  

  CAX = 0.0
  CAY = 0.0
    
!  DO ISS = 1,MSC
!    DO ID = 1,MDC
!!
!!--Initialize Fluxes-----------------------------------------------------------!
!!
!      XFLUX = 0.0
!      PSPX  = 0.0 
!      PSPY  = 0.0 
!# if defined(PLBC)
!   call replace_N32(N32,ID,ISS)
!# endif
!
  IF(.NOT. HIGH_LATITUDE_WAVE)THEN
  
  DO I = 1,M
    DO J=1,NTSN(I)-1
      I1=NBSN(I,J)
      I2=NBSN(I,J+1)
!!$          IF(DEP2(I1) <= DEPMIN .AND. DEP2(I2) > DEPMIN)THEN
!!$           FF1=0.5*(N32(ID,ISS,I)+N32(ID,ISS,I2))
!!$          ELSE IF(DEP2(I1) > DEPMIN .AND. DEP2(I2) <= DEPMIN)THEN
!!$           FF1=0.5*(N32(ID,ISS,I1)+N32(ID,ISS,I))
!!$          ELSE IF(DEP2(I1) <= DEPMIN .AND. DEP2(I2) <= DEPMIN)THEN
!!$           FF1=N32(ID,ISS,I)
!!$          ELSE
!!$           FF1=0.5*(N32(ID,ISS,I1)+N32(ID,ISS,I2))
!!$          END IF
      FF1(:,:)=0.5*(N32(:,:,I1)+N32(:,:,I2))
!!$#         if defined (SPHERICAL)
!!$          XTMP  = VX(I2)*TPI-VX(I1)*TPI
!!$       XTMP1 = VX(I2)-VX(I1)
!!$       IF(XTMP1 >  180.0)THEN
!!$         XTMP = -360.0*TPI+XTMP
!!$       ELSE IF(XTMP1 < -180.0)THEN
!!$         XTMP =  360.0*TPI+XTMP
!!$       END IF  
!!$          TXPI=XTMP*COS(DEG2RAD*VY(I))
!!$          TYPI=(VY(I1)-VY(I2))*TPI
!!$          PSPX(I)=PSPX(I)+FF1*TYPI
!!$          PSPY(I)=PSPY(I)+FF1*TXPI
!!$#         else
!!$          PSPX(I)=PSPX(I)+FF1*(VY(I1)-VY(I2))
!!$          PSPY(I)=PSPY(I)+FF1*(VX(I2)-VX(I1))
!!$#         endif
      PSPX(:,:,I)=PSPX(:,:,I)+FF1(:,:)*DLTYTRIE(I,J)
      PSPY(:,:,I)=PSPY(:,:,I)+FF1(:,:)*DLTXTRIE(I,J)
    END DO
    PSPX(:,:,I)=PSPX(:,:,I)/ART2(I)
    PSPY(:,:,I)=PSPY(:,:,I)/ART2(I)
  END DO
!# if defined (PLBC)
!  PSPY=0.0_SP
!# endif
  DO I=1,NCV_I
    I1 = NTRG(I)
    IA = NIEC(I,1)
    IB = NIEC(I,2)
!!$     XI = 0.5*(XIJE(I,1)+XIJE(I,2))
!!$     YI = 0.5*(YIJE(I,1)+YIJE(I,2))
!!$      
!!$#       if defined (SPHERICAL)
!!$        X1_DP=XIJE(I,1)
!!$        Y1_DP=YIJE(I,1)
!!$        X2_DP=XIJE(I,2)
!!$        Y2_DP=YIJE(I,2)
!!$        XII = XCG2(I)
!!$        YII = YCG2(I)
!!$        XI=XII               
!!$        XTMP  = XI*TPI-VX(IA)*TPI
!!$        XTMP1 = XI-VX(IA)
!!$        IF(XTMP1 >  180.0)THEN
!!$          XTMP = -360.0*TPI+XTMP
!!$        ELSE IF(XTMP1 < -180.0)THEN
!!$          XTMP =  360.0*TPI+XTMP
!!$        END IF        
!!$
!!$        DXA=XTMP*VAL_COS_VY(IA)    
!!$        DYA=(YI-VY(IA))*TPI
!!$        XTMP  = XI*TPI-VX(IB)*TPI
!!$        XTMP1 = XI-VX(IB)
!!$        IF(XTMP1 >  180.0)THEN
!!$          XTMP = -360.0*TPI+XTMP
!!$        ELSE IF(XTMP1 < -180.0)THEN
!!$          XTMP =  360.0*TPI+XTMP
!!$        END IF        
!!$
!!$        DXB=XTMP*VAL_COS_VY(IB) 
!!$        DYB=(YI-VY(IB))*TPI
!!$#       else
!!$        DXA = XI - VX(IA)
!!$        DYA = YI - VY(IA)
!!$     DXB = XI - VX(IB)
!!$     DYB = YI - VY(IB)
!!$#       endif
!!$      
!!$        FIJ1=N32(ID,ISS,IA)+DXA*PSPX(IA)+DYA*PSPY(IA)
!!$        FIJ2=N32(ID,ISS,IB)+DXB*PSPX(IB)+DYB*PSPY(IB)

! DEVELOPMENT TESTING - FIRST ORDER WAVE ACTION ADVECTION IN GEOG. SPACE
# if defined (DEVELOP1)
    FIJ1(:,:,1)=N32(:,:,IA)
    FIJ1(:,:,2)=N32(:,:,IB)
# else
    FIJ1(:,:,1)=N32(:,:,IA)+DLTXNCVE(I,1)*PSPX(:,:,IA)+DLTYNCVE(I,1)*PSPY(:,:,IA)
    FIJ1(:,:,2)=N32(:,:,IB)+DLTXNCVE(I,2)*PSPX(:,:,IB)+DLTYNCVE(I,2)*PSPY(:,:,IB)
# endif
	 
!DO ISS=1,MSC
!DO ID=1,MDC
!        AC1MIN=MINVAL(N32(ID,ISS,NBSN(IA,1:NTSN(IA)-1)))
!        AC1MIN=MIN(AC1MIN, N32(ID,ISS,IA))
!        AC1MAX=MAXVAL(N32(ID,ISS,NBSN(IA,1:NTSN(IA)-1)))
!        AC1MAX=MAX(AC1MAX, N32(ID,ISS,IA))
!        AC2MIN=MINVAL(N32(ID,ISS,NBSN(IB,1:NTSN(IB)-1)))
!        AC2MIN=MIN(AC2MIN, N32(ID,ISS,IB))
!        AC2MAX=MAXVAL(N32(ID,ISS,NBSN(IB,1:NTSN(IB)-1)))
!        AC2MAX=MAX(AC2MAX, N32(ID,ISS,IB))
!        IF(FIJ1(ID,ISS,1) < AC1MIN) FIJ1(ID,ISS,1) = AC1MIN
!        IF(FIJ1(ID,ISS,1) > AC1MAX) FIJ1(ID,ISS,1) = AC1MAX
!        IF(FIJ1(ID,ISS,2) < AC2MIN) FIJ1(ID,ISS,2) = AC2MIN
!        IF(FIJ1(ID,ISS,2) > AC2MAX) FIJ1(ID,ISS,2) = AC2MAX
!END DO
!END DO

    DEPLOC = (DEP2(NV(I1,1))+DEP2(NV(I1,2))+DEP2(NV(I1,3)))/3.0
    IF(DEPLOC <= DEPMIN)THEN
!   *** depth is negative ***
!     KWAVEL = -1.
     CGO_1D   = 0.
    ELSE
!   *** call KSCIP1 to compute KWAVE and CGO ***
     CALL KSCIP2(SPCSIG,DEPLOC,CGO_1D)
    END IF

!        CALL SWAPAR1(I1,ISS,ID,DEP2(1),KWAVELOC,CGLOC)
	  
#   if !defined (TWO_D_MODEL)
    UTMP = (COMPDA(NV(I1,1),JVX2)+COMPDA(NV(I1,2),JVX2)+        &
            COMPDA(NV(I1,3),JVX2))/3.0
    VTMP = (COMPDA(NV(I1,1),JVY2)+COMPDA(NV(I1,2),JVY2)+        &
            COMPDA(NV(I1,3),JVY2))/3.0
#   else
    UTMP=UA(I1)
    VTMP=VA(I1)
#   endif	
	  
    DO ID=1,MDC
      CGO_2D(ID,:)=CGO_1D
    END DO  
    SPCDIR2=SPCDIR(:,2)
    SPCDIR3=SPCDIR(:,3)
    CALL SPROXY2(CANX_2D ,CANY_2D ,   &
                 CGO_2D  ,SPCDIR2,SPCDIR3,UTMP ,VTMP )

    UN = CANY_2D*DLTXE(I) - CANX_2D*DLTYE(I)  
!# if defined (PLBC)
!    UN =  - CANX*DLTYE(I)
!# endif
    EXFLUX=-UN*((1.0+SIGN(REAL(1.0,SP),UN))*FIJ1(:,:,2)+(1.0-SIGN(REAL(1.0,SP),UN))*FIJ1(:,:,1))*0.5

#   if defined (WAVE_ONLY)
    IF(ICEIN_ON)THEN
      C_ICE=(CICE_TEST(IA)+CICE_TEST(IB))/2.0_SP
      if(C_ICE.gt.1e-8)then
        CALL ACTT2(IA,IB,RATIO)
        DO ID=1,MDC
          DO ISS=1,MSC
            EXFLUX_TMP(ID,ISS)=EXFLUX(ID,ISS)*RATIO(ISS)
            if (EXFLUX(ID,ISS).ge.0)then   
              XFLUX(ID,ISS,IA)=XFLUX(ID,ISS,IA)+EXFLUX(ID,ISS)
              XFLUX(ID,ISS,IB)=XFLUX(ID,ISS,IB)-EXFLUX_TMP(ID,ISS)
            else
              XFLUX(ID,ISS,IA)=XFLUX(ID,ISS,IA)+EXFLUX_TMP(ID,ISS)
              XFLUX(ID,ISS,IB)=XFLUX(ID,ISS,IB)-EXFLUX(ID,ISS)
            end if
          END DO
        END DO
      else
        XFLUX(:,:,IA)=XFLUX(:,:,IA)+EXFLUX
        XFLUX(:,:,IB)=XFLUX(:,:,IB)-EXFLUX
      endif
    ELSE
      XFLUX(:,:,IA)=XFLUX(:,:,IA)+EXFLUX
      XFLUX(:,:,IB)=XFLUX(:,:,IB)-EXFLUX     
    END IF  
#   elif defined (ICE)                                 
    C_ICE=(aice(IA,1)+aice(IB,1))/2.0_SP
!    C_ICE=0.0
    if(C_ICE.gt.1e-8)then
      CALL ACTT2(IA,IB,RATIO)
      DO ID=1,MDC
        DO ISS=1,MSC
          EXFLUX_TMP(ID,ISS)=EXFLUX(ID,ISS)*RATIO(ISS)
          if (EXFLUX(ID,ISS).ge.0)then   
            XFLUX(ID,ISS,IA)=XFLUX(ID,ISS,IA)+EXFLUX(ID,ISS)
            XFLUX(ID,ISS,IB)=XFLUX(ID,ISS,IB)-EXFLUX_TMP(ID,ISS)
          else
            XFLUX(ID,ISS,IA)=XFLUX(ID,ISS,IA)+EXFLUX_TMP(ID,ISS)
            XFLUX(ID,ISS,IB)=XFLUX(ID,ISS,IB)-EXFLUX(ID,ISS)
          end if
        END DO
      END DO
    else
      XFLUX(:,:,IA)=XFLUX(:,:,IA)+EXFLUX
      XFLUX(:,:,IB)=XFLUX(:,:,IB)-EXFLUX
    endif
#   else
    XFLUX(:,:,IA)=XFLUX(:,:,IA)+EXFLUX
    XFLUX(:,:,IB)=XFLUX(:,:,IB)-EXFLUX
#   endif
 
  END DO

  ELSE
   
# if defined(SPHERICAL)
  DO I = 1,M
    DO J=1,NTSN(I)-1
      I1=NBSN(I,J)
      I2=NBSN(I,J+1)
      FF1(:,:)=0.5*(N32(:,:,I1)+N32(:,:,I2))
!--wave energy advection at high latitudes yzhang-------------------------------!
      IF(VY(I1)>80.AND.VY(I2)>80) THEN
        CALL ADV_HL(FF1,I,I1,I2)
      END IF
!--end--------------------------------------------------------------------------!
      PSPX(:,:,I)=PSPX(:,:,I)+FF1(:,:)*DLTYTRIE(I,J)
      PSPY(:,:,I)=PSPY(:,:,I)+FF1(:,:)*DLTXTRIE(I,J)
    END DO
    PSPX(:,:,I)=PSPX(:,:,I)/ART2(I)
    PSPY(:,:,I)=PSPY(:,:,I)/ART2(I)
  END DO  
  DO I=1,NCV_I
    I1 = NTRG(I)
    IA = NIEC(I,1)
    IB = NIEC(I,2)

!!--wave energy at high latitudes yzhang--------------------------------------
    IF((VY(IA) > 80.0_SP) .AND. (VY(IB) > 80.0))THEN
      CALL ADV_HL_FIJ(N32_TMPP3,N32_TMPP2,I,IA,IB)
      FIJ1A(:,:,IA)=N32(:,:,IA)+DLTXNCVE(I,1)*PSPX(:,:,IA)+DLTYNCVE(I,1)*PSPY(:,:,IA)
      FIJ2A(:,:,IB)=N32_TMPP3(:,:)+DLTXNCVE(I,2)*PSPX(:,:,IB)+DLTYNCVE(I,2)*PSPY(:,:,IB)
      FIJ1B(:,:,IA)=N32_TMPP2(:,:)+DLTXNCVE(I,1)*PSPX(:,:,IA)+DLTYNCVE(I,1)*PSPY(:,:,IA)
      FIJ2B(:,:,IB)=N32(:,:,IB)+DLTXNCVE(I,2)*PSPX(:,:,IB)+DLTYNCVE(I,2)*PSPY(:,:,IB)
    ELSE
      FIJ1A(:,:,IA)=N32(:,:,IA)+DLTXNCVE(I,1)*PSPX(:,:,IA)+DLTYNCVE(I,1)*PSPY(:,:,IA)
      FIJ2A(:,:,IB)=N32(:,:,IB)+DLTXNCVE(I,2)*PSPX(:,:,IB)+DLTYNCVE(I,2)*PSPY(:,:,IB)
      FIJ1B(:,:,IA)=N32(:,:,IA)+DLTXNCVE(I,1)*PSPX(:,:,IA)+DLTYNCVE(I,1)*PSPY(:,:,IA)
      FIJ2B(:,:,IB)=N32(:,:,IB)+DLTXNCVE(I,2)*PSPX(:,:,IB)+DLTYNCVE(I,2)*PSPY(:,:,IB)
    END IF
    DEPLOC = (DEP2(NV(I1,1))+DEP2(NV(I1,2))+DEP2(NV(I1,3)))/3.0
    IF(DEPLOC <= DEPMIN)THEN
!   *** depth is negative ***
!        KWAVEL = -1.
      CGO_1D   = 0.
    ELSE
!   *** call KSCIP1 to compute KWAVE and CGO ***
      CALL KSCIP2(SPCSIG,DEPLOC,CGO_1D)
    END IF
        
#   if !defined (TWO_D_MODEL)
    UTMP = (COMPDA(NV(I1,1),JVX2)+COMPDA(NV(I1,2),JVX2)+        &
            COMPDA(NV(I1,3),JVX2))/3.0
    VTMP = (COMPDA(NV(I1,1),JVY2)+COMPDA(NV(I1,2),JVY2)+        &
            COMPDA(NV(I1,3),JVY2))/3.0
#   else
    UTMP=UA(I1)
    VTMP=VA(I1)
#   endif   
    DLTXETMP=DLTXE(I)                                 !*COS(DEG2RAD*YC(I1))       
    DLTXEA=DLTXE(I)/COS(DEG2RAD*YC(I1))*COS(DEG2RAD*VY(IA))
    DLTXEB=DLTXE(I)/COS(DEG2RAD*YC(I1))*COS(DEG2RAD*VY(IB))
    DLTYETMP=DLTYE(I)
    DO ID=1,MDC
      CGO_2D(ID,:)=CGO_1D
    END DO
    SPCDIR2=SPCDIR(:,2)
    SPCDIR3=SPCDIR(:,3)
    CALL SPROXY3(CANXAB ,CANYA, CANYB ,   &           !yzhang_w3
            CGO_2D  ,SPCDIR2,SPCDIR3,UTMP ,VTMP,DLTYETMP,DLTXETMP,DLTXEA,DLTXEB)
    UNA = CANYA - CANXAB
    UNBB = CANYB - CANXAB
    EXFLUXA=-UNA*((1.0+SIGN(REAL(1.0,SP),UNA))*FIJ2A(:,:,IB)+(1.0-SIGN(REAL(1.0,SP),UNA))*FIJ1A(:,:,IA))*0.5
    EXFLUXB=-UNBB*((1.0+SIGN(REAL(1.0,SP),UNBB))*FIJ2B(:,:,IB)+(1.0-SIGN(REAL(1.0,SP),UNBB))*FIJ1B(:,:,IA))*0.5

#   if defined (ICE)                                 
    C_ICE=(aice(IA,1)+aice(IB,1))/2
        if(C_ICE.gt.1e-8)then
        CALL ACTT2(IA,IB,RATIO)
        DO ID=1,MDC
         DO ISS=1,MSC
        EXFLUX_TMPA(ID,ISS)=EXFLUXA(ID,ISS)*RATIO(ISS)
        EXFLUX_TMPB(ID,ISS)=EXFLUXB(ID,ISS)*RATIO(ISS)
        if (EXFLUXA(ID,ISS).ge.0)then   !  .ge.
        Sice(ID,ISS,IB)=Sice(ID,ISS,IB)+EXFLUXB(ID,ISS)*(1.0-RATIO(ISS))
        else
        Sice(ID,ISS,IA)=Sice(ID,ISS,IA)-EXFLUXA(ID,ISS)*(1.0-RATIO(ISS))
        end if
        END DO
        END DO
!    IF(C_ICE > 1.0e-8)THEN
!      CALL ACTT2(IA,IB,RATIO)
!      DO ID=1,MDC
!        DO ISS=1,MSC
!          EXFLUX_TMPA(ID,ISS)=EXFLUXA(ID,ISS)*RATIO(ISS)
!          EXFLUX_TMPB(ID,ISS)=EXFLUXB(ID,ISS)*RATIO(ISS)
!          IF(EXFLUXA(ID,ISS) > 0.0_SP)THEN
!            XFLUX(ID,ISS,IA)=XFLUX(ID,ISS,IA)+EXFLUXA(ID,ISS)
!            XFLUX(ID,ISS,IB)=XFLUX(ID,ISS,IB)-EXFLUX_TMPB(ID,ISS)
!          ELSE
!            XFLUX(ID,ISS,IA)=XFLUX(ID,ISS,IA)+EXFLUX_TMPA(ID,ISS)
!            XFLUX(ID,ISS,IB)=XFLUX(ID,ISS,IB)-EXFLUXB(ID,ISS)
!          END IF
!        END DO
!      END DO
      XFLUX(:,:,IA)=XFLUX(:,:,IA)+EXFLUXA
      XFLUX(:,:,IB)=XFLUX(:,:,IB)-EXFLUXB
    ELSE
      XFLUX(:,:,IA)=XFLUX(:,:,IA)+EXFLUXA
      XFLUX(:,:,IB)=XFLUX(:,:,IB)-EXFLUXB
    END IF
#   else
    XFLUX(:,:,IA)=XFLUX(:,:,IA)+EXFLUXA
    XFLUX(:,:,IB)=XFLUX(:,:,IB)-EXFLUXB
#   endif
!--end wave energy at high latitudes yzhang--------------------------------------
  END DO

!--wave energy advection at the north pole yzhang--------------------------------!
  DO ISS=1,MSC
    DO ID=1,MDC
      XFLUX_1D=XFLUX(ID,ISS,:)
      PSPX_1D=PSPX(ID,ISS,:)
      PSPY_1D=PSPY(ID,ISS,:)
      CALL ADV_N_XY(XFLUX_1D,PSPX_1D,PSPY_1D,ISS,ID,DEP2,SPCDIR,N32)  !yzhang
    END DO
  END DO
!--end--------------------------------------------------------------------------!
# endif                     
  
  END IF  
	 
  DO ISS=1,MSC
    DO ID=1,MDC
      DO I = 1,M
        IF(ISONB_W(I) /= 0)THEN
!#          if !defined(PLBC)
           DEPLOC = DEP2(I)  
!#          else        
!           CALL replace_node_wave(I,IG2)
!           !print*,'I,IG2===',I,IG2
!           DEPLOC = DEP2(IG2)
!#          endif           
          IF(DEPLOC <= DEPMIN)THEN
!         *** depth is negative ***
            KWAVELOC = -1.                                           
            CGLOC   = 0.                                            
          ELSE
!         *** call KSCIP1 to compute KWAVE and CGO ***
            SPCSIGL = SPCSIG(ISS)
            CALL KSCIP1(1,SPCSIGL,DEPLOC,KWAVELOC,CGLOC,NN,ND)                                 
          ENDIF
          CAX(ID,ISS) = CGLOC * SPCDIR(ID,2)
          CAY(ID,ISS) = CGLOC * SPCDIR(ID,3)
!
!         --- adapt the velocities in case of diffraction
!
          IF(IDIFFR == 1 .AND. PDIFFR(3) /= 0.)THEN 
!JQI            CAX(ID,ISS) = CAX(ID,ISS)*DIFPARAM(I)      
!JQI            CAY(ID,ISS) = CAY(ID,ISS)*DIFPARAM(I)      
          END IF
!
!         --- ambient currents added
!
          IF(ICUR == 1)THEN 
#         if !defined (TWO_D_MODEL)
            CAX(ID,ISS) = CAX(ID,ISS) + COMPDA(I,JVX2)
            CAY(ID,ISS) = CAY(ID,ISS) + COMPDA(I,JVY2)
#         else
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
            CAX(ID,ISS) = CAX(ID,ISS) + UA_NODE
            CAY(ID,ISS) = CAY(ID,ISS) + VA_NODE
#         endif	    
          END IF
	    
          CANX = CAX(ID,ISS)
	  CANY = CAY(ID,ISS)
!          write(100,*)I,ID,ISS,CANX,CANY
	 
	  FIJ1_TMP = N32(ID,ISS,I)
	 
          IF(NBSN(I,NTSN(I)-1) > M)THEN
#           if defined (SPHERICAL)
            XTMP1 = VX(I)-VX(NBSN(I,2))
            XTMP = XTMP1*TPI
	    IF(XTMP1 > 180.0)THEN
	      XTMP = -360.0*TPI+XTMP
            ELSE IF(XTMP1 < -180.0)THEN
	      XTMP =  360.0*TPI+XTMP
	    END IF
	    DLTXE_TMP = XTMP*COS(DEG2RAD*VY(I))
            DLTYE_TMP = (VY(I)-VY(NBSN(I,2)))*TPI
	    UN1 = CANY*DLTXE_TMP-CANX*DLTYE_TMP
#           else

	    UN1 = CANY*(VX(I)-VX(NBSN(I,2)))      &
	        -CANX*(VY(I)-VY(NBSN(I,2)))
!# if defined (PLBC)
!            UN1 =-CANX*(VY(I)-VY(NBSN(I,2)))
!# endif
#           endif
          ELSE IF(NBSN(I,2) > M)THEN
#           if defined (SPHERICAL)
            XTMP1 = VX(NBSN(I,NTSN(I)-1))-VX(I)
            XTMP = XTMP1*TPI
	    IF(XTMP1 > 180.0)THEN
	      XTMP = -360.0*TPI+XTMP
	    ELSE IF(XTMP1 < -180.0)THEN
	      XTMP =  360.0*TPI+XTMP
	    END IF
	    DLTXE_TMP = XTMP*COS(DEG2RAD*VY(I))
            DLTYE_TMP = (VY(NBSN(I,NTSN(I)-1))-VY(I))*TPI
	    UN1 = CANY*DLTXE_TMP-CANX*DLTYE_TMP
#           else
	    UN1 = CANY*(VX(NBSN(I,NTSN(I)-1))-VX(I))      &
	        -CANX*(VY(NBSN(I,NTSN(I)-1))-VY(I))
!# if defined (PLBC)
!            UN1 = -CANX*(VY(NBSN(I,NTSN(I)-1))-VY(I))
!# endif
#           endif
          ELSE
#           if defined (SPHERICAL)
            XTMP1 = VX(NBSN(I,NTSN(I)-1))-VX(NBSN(I,2))
            XTMP = XTMP1*TPI
	    IF(XTMP1 > 180.0)THEN
	      XTMP = -360.0*TPI+XTMP
	    ELSE IF(XTMP1 < -180.0)THEN
	      XTMP =  360.0*TPI+XTMP
	    END IF
	    DLTXE_TMP = XTMP*COS(DEG2RAD*VY(I))
            DLTYE_TMP = (VY(NBSN(I,NTSN(I)-1))-VY(NBSN(I,2)))*TPI
	    UN1 = CANY*DLTXE_TMP-CANX*DLTYE_TMP
#           else
	    UN1 = CANY*(VX(NBSN(I,NTSN(I)-1))-VX(NBSN(I,2)))      &
	        -CANX*(VY(NBSN(I,NTSN(I)-1))-VY(NBSN(I,2)))
!# if defined (PLBC)
!            UN1 =-CANX*(VY(NBSN(I,NTSN(I)-1))-VY(NBSN(I,2)))
!# endif
#           endif
          END IF
	     
	  UN1 = 0.5*UN1
	     
!          EXFLUX   = MAX(0.0,-UN1*FIJ1_TMP)
          XFLUX(ID,ISS,I) = XFLUX(ID,ISS,I)+MAX(0.0,-UN1*FIJ1_TMP)
        END IF	 
      END DO  
    END DO
  END DO
!
!-Accumulate Fluxes at Boundary Nodes
!
  DO ISS=1,MSC
    DO ID=1,MDC
#     if defined (MULTIPROCESSOR)
      DO I=1,MT
        XFLUX_TMP(I)=XFLUX(ID,ISS,I)
      END DO
      IF(PAR)CALL NODE_MATCH(0,NBN,BN_MLT,BN_LOC,BNC,MT,1,MYID,NPROCS,XFLUX_TMP)
      DO I=1,MT
        XFLUX(ID,ISS,I)=XFLUX_TMP(I)
      END DO
#     endif
    END DO
  END DO

  DO I = 1,M
!--Update Action Density ------------------------------------------------------!
!
    IF(DEP2(I) > DEPMIN)THEN
      AC2(:,:,I) = N32(:,:,I)-XFLUX(:,:,I)/ART1(I)*DTW
      AC2(:,:,I) = MAX(0.0,AC2(:,:,I))
    ELSE
      AC2(:,:,I) = 0.0_SP
    END IF	
  END DO

!# if defined(PLBC)
!   CALL replace_ac2(ID,ISS)
!# endif

  DO ISS=1,MSC
    DO ID=1,MDC
#     if defined (MULTIPROCESSOR)
      IF(PAR)THEN
        AC2LOC = 0.0 
        DO I = 1,MT
          AC2LOC(I) = AC2(ID,ISS,I)
        END DO
   
        CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,1,MYID,NPROCS,AC2LOC)
!        CALL EXCHANGE(NC,MT,1,MYID,NPROCS,AC2LOC)
        CALL AEXCHANGE(NC,MYID,NPROCS,AC2LOC)

        AC2(ID,ISS,:) = 0.0
        DO I = 1,MT
          AC2(ID,ISS,I) = AC2LOC(I)
        END DO
      END IF
#     endif
    END DO
  END DO

!----yzhang----
# if defined (ICE)
# if defined (MULTIPROCESSOR)
  DO ISS=1,MSC
    DO ID=1,MDC
      DO I=1,MT
        Sice_TMP(I)=Sice(ID,ISS,I)
      END DO
      IF(PAR)CALL NODE_MATCH(0,NBN,BN_MLT,BN_LOC,BNC,MT,1,MYID,NPROCS,Sice_TMP)
      DO I=1,MT
        Sice(ID,ISS,I)=Sice_TMP(I)
      END DO
    END DO
  END DO
# endif

!--Update Sice ------------------------------------------------------!
!
  DO I = 1,M
    IF(DEP2(I) > DEPMIN)THEN
     Sice(:,:,I)=Sice(:,:,I)/ART1(I)*DTW
    END IF
  END DO

# if defined (MULTIPROCESSOR)
  IF(PAR)THEN
    DO ISS=1,MSC
      DO ID=1,MDC
        SiceLOC= 0.0
        DO I = 1,MT
          SiceLOC(I)=Sice(ID,ISS,I)
        END DO
        CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,1,MYID,NPROCS,SiceLOC)
!        CALL EXCHANGE(NC,MT,1,MYID,NPROCS,AC2LOC)
        CALL AEXCHANGE(NC,MYID,NPROCS,SiceLOC)
        Sice(ID,ISS,:)=0.0
        DO I = 1,MT
          Sice(ID,ISS,I) = SiceLOC(I)
        END DO
      END DO
    END DO
  END IF
# endif
  DEALLOCATE(SiceLOC,Sice_TMP)
# endif
!----yzhang----
  DEALLOCATE(DEP2,AC2LOC,XFLUX_TMP,EXFLUX,FF1,FIJ1,UN)
  DEALLOCATE(XFLUX,PSPX,PSPY)
  IF(HIGH_LATITUDE_WAVE)THEN
    DEALLOCATE(EXFLUXA,EXFLUXB,FIJ2)
    DEALLOCATE(FIJ1A,FIJ1B,FIJ2A,FIJ2B,UNA,UNBB)
    DEALLOCATE(N32_TMPP2,N32_TMPP3)
  END IF
   
  RETURN
  END SUBROUTINE ADV_N

!==============================================================================|
  SUBROUTINE ADV_HL(FF1,I,I1,I2)
! This subroutine is for wave advection at high latitude, 
! in order to solve the invalid scalar assumption problem.
! yzhang
!==========YZHANG_NORTHPOLE_TESTING================================
  USE VARS_WAVE
# if defined(MULTIPROCESSOR)
  USE MOD_PAR
# endif
# if defined (SPHERICAL)   
  USE MOD_SPHERICAL
# endif   
  USE SWCOMM3
  USE M_GENARR
  USE MOD_NORTHPOLE
  
  IMPLICIT NONE
  REAL(SP) :: FF1(MDC,MSC),N32_TMPP(MDC,MSC),N32_TMPP4(MDC,MSC)
  INTEGER  :: ISS, ID,I1,I2,IDT,IDD,III,I,IDD1,IDD2
  REAL     :: UL_DEGREE,CENTER_DEGREE,DL_DEGREE,DIFLAT1,DIFLAT2,RATELAT1,RATELAT2


  DO ISS = 1,MSC  
    DO ID = 1,MDC
      IF(I1==NODE_NORTHPOLE)THEN
        UL_DEGREE=22.5
        CENTER_DEGREE=0
        DL_DEGREE=337.5

        IF (VX(I)<UL_DEGREE .OR. VX(I)>DL_DEGREE)THEN
          IDT=(MDC*2/8)
          IDD=ID+IDT
!          IDD=ID
          IF(IDD>MDC)THEN
            IDD=IDD-MDC
          END IF
          N32_TMPP(ID,ISS)=N32(IDD,ISS,I1)
        END IF
        DO III=1,7
          CENTER_DEGREE=III*45.0
          UL_DEGREE=CENTER_DEGREE+22.5
          DL_DEGREE=CENTER_DEGREE-22.5
          IF(VX(I)<UL_DEGREE .AND. VX(I)>DL_DEGREE) THEN
            IDT=(MDC*(III+2)/8)
!            IDT=(MDC*III/8)
            IF(IDT>=MDC)THEN
              IDT=IDT-MDC
            END IF
            IDD=ID+IDT
            IF(IDD>MDC)THEN
              IDD=IDD-MDC
            END IF
            N32_TMPP(ID,ISS)=N32(IDD,ISS,I1)
          END IF
        END DO
        FF1(ID,ISS)=0.5*(N32_TMPP(ID,ISS)+N32(ID,ISS,I2))
      ELSE IF(I2==NODE_NORTHPOLE)THEN
        UL_DEGREE=22.5
        CENTER_DEGREE=0
        DL_DEGREE=337.5
        IF (VX(I)<UL_DEGREE .OR. VX(I)>DL_DEGREE)THEN
          IDT=(MDC*2/8)
          IDD=ID+IDT
!          IDD=ID
          IF(IDD>MDC)THEN
            IDD=IDD-MDC
          END IF
          N32_TMPP(ID,ISS)=N32(IDD,ISS,I2)
        END IF
        DO III=1,7
          CENTER_DEGREE=III*45.0
          UL_DEGREE=CENTER_DEGREE+22.5
          DL_DEGREE=CENTER_DEGREE-22.5
          IF(VX(I)<UL_DEGREE .AND. VX(I)>DL_DEGREE) THEN
            IDT=(MDC*(III+2)/8)
!            IDT=(MDC*III/8)
            IF(IDT>=MDC)THEN
              IDT=IDT-MDC
            END IF
            IDD=ID+IDT
            IF(IDD>MDC)THEN
              IDD=IDD-MDC
            END IF
            N32_TMPP(ID,ISS)=N32(IDD,ISS,I2)
          END IF
        END DO
        FF1(ID,ISS)=0.5*(N32_TMPP(ID,ISS)+N32(ID,ISS,I1))
      ELSE IF(VY(I)>80.AND.I1/=NODE_NORTHPOLE.AND.I2/=NODE_NORTHPOLE)THEN
        DIFLAT1=VX(I)-VX(I1)
        DIFLAT2=VX(I)-VX(I2)
        IF(DIFLAT1 >  180.0)THEN
          DIFLAT1 = -360.0+DIFLAT1
        ELSE IF(DIFLAT1 < -180.0)THEN
          DIFLAT1 =  360.0+DIFLAT1
        END IF
        IF(DIFLAT2 >  180.0)THEN
          DIFLAT2 = -360.0+DIFLAT2
        ELSE IF(DIFLAT2 < -180.0)THEN
          DIFLAT2 =  360.0+DIFLAT2
        END IF
        RATELAT1=DIFLAT1/360.
        RATELAT2=DIFLAT2/360.
        IDD1=ANINT(RATELAT1*MDC)+ID
        IDD2=ANINT(RATELAT2*MDC)+ID
        IF(IDD1<=0)THEN
          IDD1=IDD1+MDC
        END IF

        IF(IDD1>MDC)THEN
          IDD1=IDD1-MDC
        END IF

        IF(IDD2<=0)THEN
          IDD2=IDD2+MDC
        END IF

        IF(IDD2>MDC)THEN
          IDD2=IDD2-MDC
        END IF
        N32_TMPP(ID,ISS)=N32(IDD1,ISS,I1)
        N32_TMPP4(ID,ISS)=N32(IDD2,ISS,I2)
        FF1(ID,ISS)=0.5*(N32_TMPP(ID,ISS)+N32_TMPP4(ID,ISS))
!      ELSE
!        FF1(ID,ISS,I)=0.5*(N32(ID,ISS,I1)+N32(ID,ISS,I2))
      END IF
    END DO 
  END DO

  RETURN
  END SUBROUTINE ADV_HL

!=========================================================================
  SUBROUTINE ADV_HL_FIJ(N32_TMPP3,N32_TMPP2,I,IA,IB)
! This subroutine is for wave energy trasport at high latitudes
! yzhang
!=========================================================================
  USE VARS_WAVE
# if defined(MULTIPROCESSOR)
  USE MOD_PAR
# endif
# if defined (SPHERICAL)   
  USE MOD_SPHERICAL
# endif   
  USE SWCOMM3
  USE M_GENARR
  USE MOD_NORTHPOLE

  IMPLICIT NONE 
  REAL     ::FF1(MDC,MSC),N32_TMPP2(MDC,MSC),N32_TMPP3(MDC,MSC)
  INTEGER  :: ISS, ID,IA,IB,IDT,IDD,III,I,IDD1,IDD2
  REAL     :: UL_DEGREE,CENTER_DEGREE,DL_DEGREE,DIFLAT1,DIFLAT2,RATELAT1,RATELAT2

 
  N32_TMPP2=N32(:,:,IA)
  N32_TMPP3=N32(:,:,IB)

  DO ISS=1,MSC
    DO ID=1,MDC
      IF(IA==NODE_NORTHPOLE)THEN
        UL_DEGREE=22.5
        CENTER_DEGREE=0
        DL_DEGREE=337.5
        IF (VX(IB)<UL_DEGREE .OR. VX(IB)>DL_DEGREE)THEN
          IDT=(MDC*2/8)
          IDD=ID+IDT
!          IDD=ID
          IF(IDD>MDC)THEN
            IDD=IDD-MDC
          END IF
          N32_TMPP2(ID,ISS)=N32(IDD,ISS,IA)
        END IF
        DO III=1,7
          CENTER_DEGREE=III*45.0
          UL_DEGREE=CENTER_DEGREE+22.5
          DL_DEGREE=CENTER_DEGREE-22.5
          IF(VX(IB)<UL_DEGREE .AND. VX(IB)>DL_DEGREE) THEN
            IDT=(MDC*(III+2)/8)
!            IDT=(MDC*III/8)
            IF(IDT>=MDC)THEN
              IDT=IDT-MDC
            END IF
            IDD=ID+IDT
            IF(IDD>MDC)THEN
              IDD=IDD-MDC
            END IF
            N32_TMPP2(ID,ISS)=N32(IDD,ISS,IA)
          END IF
        END DO
      ELSE IF(IB==NODE_NORTHPOLE)THEN
        UL_DEGREE=22.5
        CENTER_DEGREE=0
        DL_DEGREE=337.5
        IF (VX(IA)<UL_DEGREE .OR. VX(IA)>DL_DEGREE)THEN
          IDT=(MDC*2/8)
          IDD=ID+IDT
!          IDD=ID
          IF(IDD>MDC)THEN
            IDD=IDD-MDC
          END IF
          N32_TMPP2(ID,ISS)=N32(IDD,ISS,IB)
        END IF
        DO III=1,7
          CENTER_DEGREE=III*45.0
          UL_DEGREE=CENTER_DEGREE+22.5
          DL_DEGREE=CENTER_DEGREE-22.5
          IF(VX(IA)<UL_DEGREE .AND. VX(IA)>DL_DEGREE) THEN
            IDT=(MDC*(III+2)/8)
!            IDT=(MDC*III/8)
            IF(IDT>=MDC)THEN
              IDT=IDT-MDC
            END IF
            IDD=ID-IDT
            IF(IDD<1)THEN
              IDD=IDD+MDC
            END IF
            N32_TMPP2(ID,ISS)=N32(IDD,ISS,IA)
          END IF
        END DO
      ELSE IF((VY(IA)>80.OR.VY(IB)>80).AND.IB/=NODE_NORTHPOLE.AND.IA/=NODE_NORTHPOLE)THEN
        DIFLAT1=VX(IB)-VX(IA)
        IF(DIFLAT1 >  180.0)THEN
          DIFLAT1 = -360.0+DIFLAT1
        ELSE IF(DIFLAT1 < -180.0)THEN
          DIFLAT1 =  360.0+DIFLAT1
        END IF
        RATELAT1=DIFLAT1/360.
        IDD1=ANINT(RATELAT1*MDC)+ID
        IF(IDD1<=0)THEN
          IDD1=IDD1+MDC
        END IF

        IF(IDD1>MDC)THEN
          IDD1=IDD1-MDC
        END IF

        N32_TMPP2(ID,ISS)=N32(IDD1,ISS,IA)
!       ELSE
!        N32_TMPP2(ID,ISS,IA)=N32(ID,ISS,IA)
      END IF
      IF(IB==NODE_NORTHPOLE)THEN
        UL_DEGREE=22.5
        CENTER_DEGREE=0
        DL_DEGREE=337.5
        IF (VX(IA)<UL_DEGREE .OR. VX(IA)>DL_DEGREE)THEN
          IDT=(MDC*2/8)
          IDD=ID+IDT
!          IDD=ID
          IF(IDD>MDC)THEN
            IDD=IDD-MDC
          END IF
          N32_TMPP3(ID,ISS)=N32(IDD,ISS,IB)
        END IF
        DO III=1,7
          CENTER_DEGREE=III*45.0
          UL_DEGREE=CENTER_DEGREE+22.5
          DL_DEGREE=CENTER_DEGREE-22.5
          IF(VX(IA)<UL_DEGREE .AND. VX(IA)>DL_DEGREE) THEN
            IDT=(MDC*(III+2)/8)
!            IDT=(MDC*III/8)
            IF(IDT>=MDC)THEN
              IDT=IDT-MDC
            END IF
            IDD=ID+IDT
            IF(IDD>MDC)THEN
              IDD=IDD-MDC
            END IF
            N32_TMPP3(ID,ISS)=N32(IDD,ISS,IB)
          END IF
        END DO
      ELSE IF(IA==NODE_NORTHPOLE)THEN
        UL_DEGREE=22.5
        CENTER_DEGREE=0
        DL_DEGREE=337.5
        IF (VX(IB)<UL_DEGREE .OR. VX(IB)>DL_DEGREE)THEN
          IDT=(MDC*2/8)
          IDD=ID+IDT
!          IDD=ID
          IF(IDD>MDC)THEN
            IDD=IDD-MDC
          END IF
          N32_TMPP3(ID,ISS)=N32(IDD,ISS,IA)
        END IF
        DO III=1,7
          CENTER_DEGREE=III*45.0
          UL_DEGREE=CENTER_DEGREE+22.5
          DL_DEGREE=CENTER_DEGREE-22.5
          IF(VX(IB)<UL_DEGREE .AND. VX(IB)>DL_DEGREE) THEN
            IDT=(MDC*(III+2)/8)
!            IDT=(MDC*III/8)
            IF(IDT>=MDC)THEN
              IDT=IDT-MDC
            END IF
            IDD=ID-IDT
            IF(IDD<1)THEN
              IDD=IDD+MDC
            END IF
            N32_TMPP3(ID,ISS)=N32(IDD,ISS,IB)
          END IF
        END DO
      ELSE IF((VY(IA)>80.OR.VY(IB)>80).AND.IB/=NODE_NORTHPOLE.AND.IA/=NODE_NORTHPOLE)THEN
        DIFLAT1=VX(IA)-VX(IB)
        IF(DIFLAT1 >  180.0)THEN
          DIFLAT1 = -360.0+DIFLAT1
        ELSE IF(DIFLAT1 < -180.0)THEN
          DIFLAT1 =  360.0+DIFLAT1
        END IF
        RATELAT1=DIFLAT1/360.
        IDD1=ANINT(RATELAT1*MDC)+ID
        IF(IDD1<=0)THEN
          IDD1=IDD1+MDC
        END IF

        IF(IDD1>MDC)THEN
          IDD1=IDD1-MDC
        END IF

        N32_TMPP3(ID,ISS)=N32(IDD1,ISS,IB)
!       ELSE
!       N32_TMPP3(ID,ISS,IB)=N32(ID,ISS,IB)
      END IF
!      FIJ1A(ID,ISS,IA)=N32(ID,ISS,IA)+DLTXNCVE(I,1)*PSPX(ID,ISS,IA)+DLTYNCVE(I,1)*PSPY(ID,ISS,IA)
!      FIJ2A(ID,ISS,IB)=N32_TMPP2(ID,ISS,IB)+DLTXNCVE(I,2)*PSPX(ID,ISS,IB)+DLTYNCVE(I,2)*PSPY(ID,ISS,IB)
    END DO
  END DO

  RETURN
  END SUBROUTINE ADV_HL_FIJ
 
!==========================================================================
!
  SUBROUTINE ACTT2(IA,IB,RATIO)
! This subroutine is for Ice induced wave attenuation 
! The parameter 'RATIO' is the dimentional wave attenuation coefficient
! yzhang
!
!****************************************************************
!
  USE M_GENARR
  USE SWCOMM2
  USE SWCOMM3
  USE SWCOMM4
  USE OCPCOMM4
  USE M_PARALL
  USE ALL_VARS
  USE VARS_WAVE
# if defined (SPHERICAL)
  USE MOD_SPHERICAL
# endif
  USE ICE_STATE

!  USE ice_state  , only :aice , vice!, floe_diam
!  USE ICE_THERM_VERTICAL ,only:hicen_old
!  USE ice_model_size ,only: ncat
# if defined (MULTIPROCESSOR)
  USE MOD_PAR
# endif

  IMPLICIT NONE 
  REAL    :: C_ICE,COES(7,9),COES2(7,9),ICETHI,SPCSIGL,DEPLOC,TWAVE
  REAL    :: ICEATT,ICEATTS(7),ICEATT1,ICEATT2,RATE,RATIO(MSC)
  REAL    :: mean_diam,mean_diam_u,mean_diam_d,floe_diam
  REAL ,parameter   :: ee=2,ff=0.9,dmin=20
  REAL(SP):: DEP3(MT),DIS,DTMP_DP,X1_DP,X2_DP,Y1_DP,Y2_DP
  INTEGER :: I,ISS,ID,ni,mm,IA,IB
 
  DEP3(1:MT) = COMPDA(1:MT,JDP2)
  COES2(1,:)=(/0.0,0.0,0.0,0.0,-0.003509,0.1473,-2.121,11.24,-22.42/)
  COES2(2,:)=(/0.0,0.0,0.0,0.0,-0.003723,0.173,-2.833,18.33,-43.89/)
  COES2(3,:)=(/0.0,0.0,0.0,0.0,-0.001264,0.07156,-1.364,9.616,-25.34/)
  COES2(4,:)=(/0.0,0.0,0.0,0.0,0.001401,-0.05012,0.5915,-3.329,5.372/)
  COES2(5,:)=(/0.0,0.0,0.0,0.0,0.0009987,-0.0403,0.5441,-3.471,6.658/)
  COES2(6,:)=(/0.0,0.0,0.0,0.0,0.0002636,-0.01267,0.1853,-1.403,2.459/)
  COES2(7,:)=(/0.0,0.0,0.0,0.0,1.551E-5,-0.002103,0.03268,-0.4304,0.3751/)

# if defined (WAVE_ONLY)
  C_ICE=(CICE_TEST(IA)+CICE_TEST(IB))/2.0_SP
# else
  C_ICE=(aice(IA,1)+aice(IB,1))/2.0_SP
# endif
  floe_diam=200
# if defined (SPHERICAL)
  Y1_DP=vy(IA)
  Y2_DP=vy(IB)
  X1_DP=vx(IA)
  X2_DP=vx(IB)
  CALL ARC(X1_DP,Y1_DP,X2_DP,Y2_DP,DTMP_DP)
  DIS=DTMP_DP
# else
  DIS=((vy(IA)-vy(IB))**2+(vx(IA)-vx(IB))**2)**0.5*2.0
# endif
  if(C_ICE>0.00001)then
    mean_diam_u=0.0
    mean_diam_d=0.0
    mean_diam=0.0
    DO mm=1,3
      mean_diam_u=mean_diam_u+(ee**2*ff)**mm*ee**(-mm)*(floe_diam+floe_diam)/2
      mean_diam_d=mean_diam_d+(ee**2*ff)**mm
    END DO
    mean_diam=mean_diam_u/mean_diam_d
    if(C_ICE.eq.1)then
      ICETHI=3.0
    else
      ICETHI=0.3*(1.0+4*C_ICE)
    endif
    DO ISS=1,MSC
      SPCSIGL=SPCSIG(ISS)
      DEPLOC=DEP3(I)
!      CALL KSCIP1(1,SPCSIGL,DEPLOC,KWAVEL,CGOL,NN,ND)
      TWAVE=1/(SPCSIGL/PI2_W)
      IF(TWAVE<=16.0)THEN
        IF(TWAVE<=6.0)TWAVE=6.0
          ICEATTS(1)=COES2(1,1)*TWAVE**8+COES2(1,2)*TWAVE**7+COES2(1,3)*TWAVE**6+COES2(1,4)*TWAVE**5+COES2(1,5)*TWAVE**4+&
          COES2(1,6)*TWAVE**3+COES2(1,7)*TWAVE**2+COES2(1,8)*TWAVE+COES2(1,9)

          ICEATTS(2)=COES2(2,1)*TWAVE**8+COES2(2,2)*TWAVE**7+COES2(2,3)*TWAVE**6+COES2(2,4)*TWAVE**5+COES2(2,5)*TWAVE**4+&
          COES2(2,6)*TWAVE**3+COES2(2,7)*TWAVE**2+COES2(2,8)*TWAVE+COES2(2,9)

          ICEATTS(3)=COES2(3,1)*TWAVE**8+COES2(3,2)*TWAVE**7+COES2(3,3)*TWAVE**6+COES2(3,4)*TWAVE**5+COES2(3,5)*TWAVE**4+&
          COES2(3,6)*TWAVE**3+COES2(3,7)*TWAVE**2+COES2(3,8)*TWAVE+COES2(3,9)

          ICEATTS(4)=COES2(4,1)*TWAVE**8+COES2(4,2)*TWAVE**7+COES2(4,3)*TWAVE**6+COES2(4,4)*TWAVE**5+COES2(4,5)*TWAVE**4+&
          COES2(4,6)*TWAVE**3+COES2(4,7)*TWAVE**2+COES2(4,8)*TWAVE+COES2(4,9)

          ICEATTS(5)=COES2(5,1)*TWAVE**8+COES2(5,2)*TWAVE**7+COES2(5,3)*TWAVE**6+COES2(5,4)*TWAVE**5+COES2(5,5)*TWAVE**4+&
          COES2(5,6)*TWAVE**3+COES2(5,7)*TWAVE**2+COES2(5,8)*TWAVE+COES2(5,9)

          ICEATTS(6)=COES2(6,1)*TWAVE**8+COES2(6,2)*TWAVE**7+COES2(6,3)*TWAVE**6+COES2(6,4)*TWAVE**5+COES2(6,5)*TWAVE**4+&
          COES2(6,6)*TWAVE**3+COES2(6,7)*TWAVE**2+COES2(6,8)*TWAVE+COES2(6,9)

          ICEATTS(7)=COES2(7,1)*TWAVE**8+COES2(7,2)*TWAVE**7+COES2(7,3)*TWAVE**6+COES2(7,4)*TWAVE**5+COES2(7,5)*TWAVE**4+&
          COES2(7,6)*TWAVE**3+COES2(7,7)*TWAVE**2+COES2(7,8)*TWAVE+COES2(7,9)
          IF (ICETHI<=0.4)THEN
            ICEATT=ICEATTS(1)
          ELSE IF(ICETHI==0.6)THEN
            ICEATT=ICEATTS(2)
          ELSE IF(ICETHI==0.8)THEN
            ICEATT=ICEATTS(3)
          ELSE IF(ICETHI==1.2)THEN
            ICEATT=ICEATTS(4)
          ELSE IF(ICETHI==1.6)THEN
            ICEATT=ICEATTS(5)
          ELSE IF(ICETHI==2.4)THEN
            ICEATT=ICEATTS(6)
          ELSE IF(ICETHI>=3.2)THEN
            ICEATT=ICEATTS(7)
          ELSE IF(ICETHI>0.4.AND.ICETHI<0.6)THEN
            ICEATT1=ICEATTS(1)
            ICEATT2=ICEATTS(2)
            RATE=(ICETHI-0.4)/0.2
            ICEATT=ICEATT1*(1-RATE)+ICEATT2*RATE
          ELSE IF(ICETHI>0.6.AND.ICETHI<0.8)THEN
            ICEATT1=ICEATTS(2)
            ICEATT2=ICEATTS(3)
            RATE=(ICETHI-0.6)/0.2
            ICEATT=ICEATT1*(1-RATE)+ICEATT2*RATE
          ELSEIF(ICETHI>0.8.AND.ICETHI<1.2)THEN
            ICEATT1=ICEATTS(3)
            ICEATT2=ICEATTS(4)
            RATE=(ICETHI-0.8)/0.4
            ICEATT=ICEATT1*(1-RATE)+ICEATT2*RATE
          ELSE IF(ICETHI>1.2.AND.ICETHI<1.6)THEN
            ICEATT1=ICEATTS(4)
            ICEATT2=ICEATTS(5)
            RATE=(ICETHI-1.2)/0.4
            ICEATT=ICEATT1*(1-RATE)+ICEATT2*RATE
          ELSE IF(ICETHI>1.6.AND.ICETHI<2.4)THEN
            ICEATT1=ICEATTS(5)
            ICEATT2=ICEATTS(6)
            RATE=(ICETHI-1.6)/0.8
            ICEATT=ICEATT1*(1-RATE)+ICEATT2*RATE
          ELSE IF(ICETHI>2.4.AND.ICETHI<3.2)THEN
            ICEATT1=ICEATTS(6)
            ICEATT2=ICEATTS(7)
            RATE=(ICETHI-2.4)/0.8
            ICEATT=ICEATT1*(1-RATE)+ICEATT2*RATE
          ENDIF
          ICEATT=exp(ICEATT)
          ICEATT=ICEATT*C_ICE/mean_diam
          RATIO(ISS)=exp(-ICEATT*DIS)
        ELSE
          ICEATT=0.02*exp(-0.386*TWAVE)*C_ICE
          RATIO(ISS)=exp(-ICEATT*DIS)
        ENDIF
      END DO
    else
      DO ISS=1,MSC
      RATIO(ISS)=1.0
    end do
  endif

  RETURN
  END SUBROUTINE ACTT2
!==============================================================================|

END MODULE MOD_ACTION_EX
