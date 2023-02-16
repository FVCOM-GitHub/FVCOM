#  if defined (WAVE_CURRENT_INTERACTION)
MODULE MOD_WAVE_CURRENT_INTERACTION
   USE MOD_PREC
   USE VARS_WAVE
#  if defined(MULTIPROCESSOR)
   USE MOD_PAR
#  endif
#  if defined (SPHERICAL)
   USE MOD_SPHERICAL
#  endif
#  if defined (WET_DRY)
   USE MOD_WD
#  endif
   USE MOD_STATION_TIMESERIES
   USE MOD_SPARSE_TIMESERIES
   
   IMPLICIT NONE
   
   REAL(SP), ALLOCATABLE :: WAVESTRX_2D(:),WAVESTRY_2D(:)
   REAL(SP), ALLOCATABLE :: WAVESTRX_3D(:,:),WAVESTRY_3D(:,:)
   REAL(SP), ALLOCATABLE :: U_STOKES_2D(:),V_STOKES_2D(:)
   REAL(SP), ALLOCATABLE :: U_STOKES_3D(:,:),V_STOKES_3D(:,:)
#  if defined (VORTEX_FORCE)
   REAL(SP), ALLOCATABLE :: W_STOKES_3D(:,:) !Xia et al. (2020)
   REAL(SP), ALLOCATABLE :: WW_STOKES_3D(:,:) !Xia et al. (2020)
   REAL(SP), ALLOCATABLE :: W_STOKES_3D_TMP(:,:) !Xia et al. (2020)
   REAL(SP), ALLOCATABLE :: WW_STOKES_3D_TMP(:,:) !Xia et al. (2020)
#  endif
   REAL(SP), ALLOCATABLE :: U_STOKES_2D_TMP(:),V_STOKES_2D_TMP(:)
   REAL(SP), ALLOCATABLE :: U_STOKES_3D_TMP(:,:),V_STOKES_3D_TMP(:,:)
   REAL(SP), ALLOCATABLE :: TPZDIST(:,:)
   REAL(SP), ALLOCATABLE :: UW10(:),VW10(:)
   REAL(SP), ALLOCATABLE :: USTW(:),USTP(:)
   REAL(SP), ALLOCATABLE :: TPX0(:),TPY0(:)
   REAL(SP), ALLOCATABLE :: TTX0(:),TTY0(:)
   REAL(SP), ALLOCATABLE :: TPX(:,:),TPY(:,:)
   REAL(SP), ALLOCATABLE :: UDOP(:),VDOP(:)
   REAL(SP), ALLOCATABLE :: UNODE(:,:),VNODE(:,:)
! for roller
   REAL(SP), ALLOCATABLE :: GAMW(:),OROLLER(:),ROLLA(:)
#  if defined (VORTEX_FORCE)
   REAL(SP), ALLOCATABLE :: ROLLA_VF(:) !Xia et al. (2020)
#  endif
! end
   REAL(SP), PARAMETER :: KDMAX = 3.0_SP ! Based on MELLOR(2015). For old code: KDMAX = 5.0_SP
   REAL(SP), PARAMETER :: WAVE_LENGTH_MIN = 0.01_SP
   REAL(SP), PARAMETER :: eps1 = 1E-14_SP

#  if defined (VORTEX_FORCE)
!  Start add
   REAL(SP), ALLOCATABLE :: DISSIP_BREAK(:),DISSIP_WCAP(:),DISSIP_ROLLER(:)
   REAL(SP), ALLOCATABLE :: XSTFLUX(:,:) !Xia et al. (2020)
   REAL(SP), ALLOCATABLE :: DISSIP_FRIC(:)
!  End add
#  endif

! for bbl
!   The Options MB_Z0BL and MB_Z0RIP should be activated concurrently.      **

 LOGICAL :: MB_BBL_USE        !  use if Meinte Blaas BBL closure                       **
 LOGICAL :: MB_CALC_ZNOT      !  use if computing bottom roughness internally          **
 LOGICAL :: MB_CALC_UB        !  use if computing bottom orbital velocity internally   **
 LOGICAL :: MB_Z0BIO          !  use if biogenic bedform roughness for ripples         **
 LOGICAL :: MB_Z0BL           !  use if bedload roughness for ripples                  **
 LOGICAL :: MB_Z0RIP          !  use if bedform roughness for ripples                  **
 
! OPTIONS for Styles and Glenn (2000) bottom boundary layer closure:        **

 LOGICAL :: SG_BBL_USE        !  use if Styles and Glenn (2000) BBL closure            **
 LOGICAL :: SG_CALC_ZNOT      !  use if computing bottom roughness internally          **
 LOGICAL :: SG_CALC_UB        !  use if computing bottom orbital velocity internally   **
 LOGICAL :: SG_LOGINT         !  use if logarithmic interpolation of (Ur,Vr)           **

! OPTIONS for the Sherwood/Signell/Warner bottom boundary layer closure:    **

 LOGICAL :: SSW_BBL_USE       !  use if Sherwood et al. BBL closure                    **
 LOGICAL :: SSW_CALC_ZNOT     !  use if computing bottom roughness internally          **
 LOGICAL :: SSW_LOGINT        !  use if logarithmic interpolation of (Ur,Vr)           **
 LOGICAL :: SSW_CALC_UB       !  use if computing bottom orbital velocity internally   **
 LOGICAL :: SSW_FORM_DRAG_COR !  use to activate form drag coefficient                 **
 LOGICAL :: SSW_ZOBIO         !  use if biogenic bedform roughness from ripples        **
 LOGICAL :: SSW_ZOBL          !  use if bedload roughness for ripples                  **
 LOGICAL :: SSW_ZORIP         !  use if bedform roughness from ripples                 **

 LOGICAL :: SGWC
 LOGICAL :: M94WC

 LOGICAL :: GM82_RIPRUF 
 LOGICAL :: N92_RIPRUF 
 LOGICAL :: R88_RIPRUF 

!=======================================================================
!                                                                      !
!  Ubot         Wind-induced, bed wave orbital U-velocity (m/s) at     !
!                 RHO-points.                                          !
!  Ur           Bottom U-momentum above bed (m/s) at RHO-points.       !
!  Vbot         Wind-induced, bed wave orbital V-velocity (m/s) at     !
!                 RHO-points.                                          !
!  Vr           Bottom V-momentum above bed (m/s) at RHO-points.       !
!  bustrc       Kinematic bottom stress (m2/s2) due currents in the    !
!                 XI-direction at RHO-points.                          !
!  bustrw       Kinematic bottom stress (m2/s2) due to wind-induced    !
!                 waves the XI-direction at horizontal RHO-points.     !
!  bustrcwmax   Kinematic bottom stress (m2/s2) due to maximum wind    !
!                 and currents in the XI-direction at RHO-points.      !
!  bvstrc       Kinematic bottom stress (m2/s2) due currents in the    !
!                 ETA-direction at RHO-points.                         !
!  bvstrw       Kinematic bottom stress (m2/s2) due to wind-induced    !
!                 waves the ETA-direction at horizontal RHO-points.    !
!  bvstrcwmax   Kinematic bottom stress (m2/s2) due to maximum wind    !
!                 and currents in the ETA-direction RHO-points.        !
!                                                                      !
!=======================================================================
  real(SP), allocatable :: Ubot(:)
  real(SP), allocatable :: Vbot(:)
  real(SP), allocatable :: Ur(:)
  real(SP), allocatable :: Vr(:)
  real(SP), allocatable :: bustrc(:)
  real(SP), allocatable :: bvstrc(:)
  real(SP), allocatable :: bustrw(:)
  real(SP), allocatable :: bvstrw(:)
  real(SP), allocatable :: bustrcwmax(:)
  real(SP), allocatable :: bvstrcwmax(:)
  real(SP), allocatable :: taucwmax(:)
  real(SP), allocatable :: bustr(:)
  real(SP), allocatable :: bvstr(:)


   CONTAINS
!====================================================================================|
   SUBROUTINE WAVE_CURRENT_SETUP
   USE MOD_SPHERICAL
   USE VARS_WAVE
   USE SWCOMM3
   IMPLICIT NONE
   
   ALLOCATE(WAVESTRX_2D(0:NT));      WAVESTRX_2D = 0.0_SP
   ALLOCATE(WAVESTRY_2D(0:NT));      WAVESTRY_2D = 0.0_SP
   ALLOCATE(WAVESTRX_3D(0:NT,KB));   WAVESTRX_3D = 0.0_SP
   ALLOCATE(WAVESTRY_3D(0:NT,KB));   WAVESTRY_3D = 0.0_SP
   ALLOCATE(U_STOKES_2D(0:NT));      U_STOKES_2D = 0.0_SP
   ALLOCATE(V_STOKES_2D(0:NT));      V_STOKES_2D = 0.0_SP
   ALLOCATE(U_STOKES_3D(0:NT,KB));   U_STOKES_3D = 0.0_SP
   ALLOCATE(V_STOKES_3D(0:NT,KB));   V_STOKES_3D = 0.0_SP
#  if defined (VORTEX_FORCE)
   ALLOCATE(W_STOKES_3D(0:NT,KB));   W_STOKES_3D = 0.0_SP !Xia et al. (2020)
   ALLOCATE(WW_STOKES_3D(0:NT,KB));  WW_STOKES_3D = 0.0_SP !Xia et al. (2020)
   ALLOCATE(XSTFLUX(0:MT,KBM1)); XSTFLUX = 0.0_SP !Xia et al. (2020)
   ALLOCATE(DISSIP_BREAK(0:MT));     DISSIP_BREAK = 0.0_SP !Xia et al. (2020)
   ALLOCATE(DISSIP_WCAP(0:MT));      DISSIP_WCAP = 0.0_SP  !Xia et al. (2020)
   ALLOCATE(DISSIP_ROLLER(0:MT));    DISSIP_ROLLER = 0.0_SP !Xia et al. (2020)
   ALLOCATE(DISSIP_FRIC(0:MT));      DISSIP_FRIC = 0.0_SP !Xia et al. (2020)
#  endif
   ALLOCATE(U_STOKES_2D_TMP(0:MT));  U_STOKES_2D_TMP = 0.0_SP
   ALLOCATE(V_STOKES_2D_TMP(0:MT));  V_STOKES_2D_TMP = 0.0_SP
   ALLOCATE(U_STOKES_3D_TMP(0:MT,KB));   U_STOKES_3D_TMP = 0.0_SP
   ALLOCATE(V_STOKES_3D_TMP(0:MT,KB));   V_STOKES_3D_TMP = 0.0_SP
#  if defined (VORTEX_FORCE)
   ALLOCATE(W_STOKES_3D_TMP(0:MT,KB)); W_STOKES_3D_TMP = 0.0_SP !Xia et al. (2020)
   ALLOCATE(WW_STOKES_3D_TMP(0:MT,KB)); WW_STOKES_3D_TMP = 0.0_SP !Xia et al. (2020)
#  endif
   ALLOCATE(GAMW(0:MT));             GAMW        = 0.0_SP
   ALLOCATE(OROLLER(0:MT));          OROLLER     = 0.0_SP
   ALLOCATE(ROLLA(0:MT));            ROLLA      = 0.0_SP
#  if defined (VORTEX_FORCE)
   ALLOCATE(ROLLA_VF(0:MT));         ROLLA_VF   = 0.0_SP !Xia et al. (2020)
#  endif
   ALLOCATE(TPZDIST(0:NT,KB));       TPZDIST     = 0.0_SP
   ALLOCATE(UW10(0:NT));             UW10        = 0.0_SP
   ALLOCATE(VW10(0:NT));             VW10        = 0.0_SP
   ALLOCATE(USTW(0:NT));             USTW        = 0.0_SP
   ALLOCATE(USTP(0:NT));             USTP        = 0.0_SP
   ALLOCATE(TPX0(0:NT));             TPX0        = 0.0_SP
   ALLOCATE(TPY0(0:NT));             TPY0        = 0.0_SP
   ALLOCATE(TTX0(0:NT));             TTX0        = 0.0_SP
   ALLOCATE(TTY0(0:NT));             TTY0        = 0.0_SP
   ALLOCATE(TPX(0:NT,KB));           TPX         = 0.0_SP
   ALLOCATE(TPY(0:NT,KB));           TPY         = 0.0_SP
   ALLOCATE(UDOP(0:MT));             UDOP        = 0.0_SP
   ALLOCATE(VDOP(0:MT));             VDOP        = 0.0_SP
   ALLOCATE(UNODE(0:MT,KB));         UNODE       = 0.0_SP
   ALLOCATE(VNODE(0:MT,KB));         VNODE       = 0.0_SP

   ALLOCATE(HSC1(0:MT));             HSC1        = 0.0_SP       
   ALLOCATE(DIRDEG1(0:MT));          DIRDEG1     = 0.0_SP
   ALLOCATE(TPEAK(0:MT))  ;          TPEAK       = 0.0_SP
   ALLOCATE(WLEN(0:MT))   ;          WLEN        = 0.0_SP
   ALLOCATE(QB1(0:MT))    ;          QB1         = 0.0_SP
   ALLOCATE(Pwave_bot(0:MT)) ;       Pwave_bot   = 0.0_SP
   ALLOCATE(Ub_swan(0:MT))  ;        Ub_swan     = 0.0_SP
   ALLOCATE(Dwave(0:MT));            Dwave       = 0.0_SP
   ALLOCATE(DIRBOT(0:MT));           DIRBOT      = 0.0_SP
   ALLOCATE(SPEC_DENSITY(0:MT,MSC)); SPEC_DENSITY= 0.0_SP

#  if !defined (WAVE_OFFLINE)
   IF(OUT_WAVE_PARTITION .OR. OUT_WAVE_PARTITION_SPARSE)THEN
    ALLOCATE(HS_WIND(0:MT));                HS_WIND             = -999       
    ALLOCATE(DIRDEG_WIND(0:MT));            DIRDEG_WIND         = -999
    ALLOCATE(TPEAK_WIND(0:MT));             TPEAK_WIND          = -999
    ALLOCATE(TPEAK_WIND_POS(0:MT));         TPEAK_WIND_POS      = -999
    ALLOCATE(HS_SWELL_ALL(0:MT,50));        HS_SWELL_ALL        = -999
    ALLOCATE(DIRDEG_SWELL_ALL(0:MT,50));    DIRDEG_SWELL_ALL    = -999
    ALLOCATE(TPEAK_SWELL_ALL(0:MT,50));     TPEAK_SWELL_ALL     = -999
    ALLOCATE(TPEAK_SWELL_POS_ALL(0:MT,50)); TPEAK_SWELL_POS_ALL = -999              
   END IF 
#  endif
!---------------Coordinates of Center Pionts around the Nodes-----------------------!
!!$# if defined (SPHERICAL)
!!$   ALLOCATE(XCA(M))             ;XCA       = ZERO
!!$   ALLOCATE(YCA(M))             ;YCA       = ZERO
!!$   ALLOCATE(XCB(M))             ;XCB       = ZERO
!!$   ALLOCATE(YCB(M))             ;YCB       = ZERO
!!$   ALLOCATE(XCC(M,20))          ;XCC       = ZERO 
!!$   ALLOCATE(YCC(M,20))          ;YCC       = ZERO 
!!$   ALLOCATE(XCD(M,20))          ;XCD       = ZERO 
!!$   ALLOCATE(YCD(M,20))          ;YCD       = ZERO 
!!$   ALLOCATE(XCE(M))             ;XCE       = ZERO
!!$   ALLOCATE(YCE(M))             ;YCE       = ZERO
!!$   ALLOCATE(XCF(M))             ;XCF       = ZERO
!!$   ALLOCATE(YCF(M))             ;YCF       = ZERO
!!$   ALLOCATE(VAL_COS_VY(M))      ;VAL_COS_VY= ZERO
!!$   CALL CAL_CENTER  
!!$# endif

   MCGRD = MGL
   RETURN
   END SUBROUTINE WAVE_CURRENT_SETUP 

!====================================================================================|  
   SUBROUTINE RADIATION_STRESS_3D

   IMPLICIT NONE

   REAL(SP), ALLOCATABLE        :: SXX(:,:),SXY(:,:),SYY(:,:)   !Jianzhong

!   REAL(SP), DIMENSION(0:MT   ) :: SXXA,SXYA,SYYA
   REAL(SP), DIMENSION(0:NT,KB) :: PSXXPX,PSXYPX,PSXYPY,PSYYPY
   REAL(SP), DIMENSION(0:NT   ) :: PSXXPXA,PSXYPXA,PSXYPYA,PSYYPYA
   REAL(SP), DIMENSION(0:NT,KB) :: PSPXPZ,PSPYPZ

   REAL(SP), DIMENSION(0:MT) :: WAVE_NUMBER,WAVE_NUMBER_X,WAVE_NUMBER_Y,SIN_DIR,COS_DIR
   REAL(SP), DIMENSION(0:MT) :: WAVE_ENERGY,KD,WAVE_C
   REAL(SP), DIMENSION(0:MT) :: O_WAVE_NUMBER
   REAL(SP), DIMENSION(0:MT) :: O_COSH,O_SINH,O_2SINH
   REAL(SP) :: EXFLUX
   INTEGER  :: I,K,IA,IB,J1,J2
   REAL(SP) :: FSS,FCS,FSC,FCC
   REAL(SP) :: CFF1,CFF2,CFF3,CFF4,CFF5,CFF6,FAC2,sum3dsxx
   REAL(SP) :: SXXIJ,SXYIJ,SYYIJ,DIJ
   
   REAL(SP) :: XTMP,XTMP1
!-------------------------------------------------------------------------------------|
  
!---------------Jianzhong----------------------------
   IF(.NOT.ALLOCATED(SXX)) ALLOCATE(SXX(0:MT,KB))
   IF(.NOT.ALLOCATED(SXY)) ALLOCATE(SXY(0:MT,KB))
   IF(.NOT.ALLOCATED(SYY)) ALLOCATE(SYY(0:MT,KB))
!----------------------------------------------------
 
   WAVE_NUMBER   = 0.0_SP   ;WAVE_NUMBER_X = 0.0_SP   ;WAVE_NUMBER_Y = 0.0_SP
   O_COSH        = 0.0_SP   ;O_SINH        = 0.0_SP   ;O_2SINH       = 0.0_SP
   O_WAVE_NUMBER = 0.0_SP

!
!  Compute wave numbers and wave energy.
!
   DO I=1,MT
    WAVE_NUMBER(I) = 2.0_SP*PI/MAX(WLEN(I),WAVE_LENGTH_MIN)
   END DO 
   O_WAVE_NUMBER = 1.0_SP/WAVE_NUMBER
   SIN_DIR       = SIN(DIRDEG1*DEG2RAD)
   COS_DIR       = COS(DIRDEG1*DEG2RAD)
   WAVE_NUMBER_X = WAVE_NUMBER*COS_DIR
   WAVE_NUMBER_Y = WAVE_NUMBER*SIN_DIR
   WAVE_ENERGY   = 0.0625_SP*GRAV_N*HSC1*HSC1
!
!  Compute wave celerity and phase velocity.
!
   DO I=1,MT
!     KD(I) = MIN(WAVE_NUMBER(I)*D(I)+eps1,KDMAX)
     KD(I) = WAVE_NUMBER(I)*D(I)+eps1
   END DO 
   
   WHERE(KD <= KDMAX) 
    WAVE_C = SQRT(GRAV_N*O_WAVE_NUMBER*TANH(KD))

    O_COSH  = 1.0_SP/COSH(KD)
    O_SINH  = 1.0_SP/SINH(KD)
    O_2SINH = 1.0_SP/SINH(2.0_SP*KD)
   ELSEWHERE
    WAVE_C = SQRT(GRAV_N*O_WAVE_NUMBER*TANH(KD))
   END WHERE
    
#if defined(WAVE_ROLLER)
   OROLLER = 0.0_SP;GAMW = 0.0_SP; ROLLA = 0.0_SP
   DO I=1,MT
     GAMW(I) = MIN(D(I)/(HSC1(I)+eps1),5.0_SP) 
     DO K=1,KBM1
        OROLLER(I)=OROLLER(I)+D(I)*DZ(I,K)*(1.0_SP-TANH((2.0_SP*ZZ(I,K)*GAMW(I))**4))
     END DO
     OROLLER(I)=1.0_SP/(OROLLER(I)+eps1)
     ROLLA(I)=0.0424*HSC1(I)*QB1(I)*WLEN(I)
   END DO
#endif
   
!----------INITIALIZE STRESS ARRAY ----------------------------------------------!

   SXX    = 0.0_SP   ;SXY    = 0.0_SP   ;SYY    = 0.0_SP
!   SXXA   = 0.0_SP   ;SXYA   = 0.0_SP   ;SYYA   = 0.0_SP
   PSXXPX = 0.0_SP   ;PSXYPX = 0.0_SP   ;PSXYPY = 0.0_SP   ;PSYYPY = 0.0_SP
   PSPXPZ = 0.0_SP   ;PSPYPZ = 0.0_SP

   DO I=1,M
     sum3dsxx=0
    IF(KD(I) <= KDMAX)THEN 
     DO K=1,KBM1
       FAC2 = 1.0_SP+ZZ(I,K)
       FCC  = COSH(KD(I)*FAC2)*O_COSH(I)
       FCS  = COSH(KD(I)*FAC2)*O_SINH(I)
       FSC  = SINH(KD(I)*FAC2)*O_COSH(I)
       FSS  = SINH(KD(I)*FAC2)*O_SINH(I)

       CFF1 = WAVE_NUMBER(I)*WAVE_ENERGY(I)
!       CFF4 = CFF1*FCS*FSS
       CFF4 = CFF1*FSC*FSS
       CFF6 = CFF1*FCS*FCC
       CFF5 = CFF1*FCS*FCC*O_WAVE_NUMBER(I)*O_WAVE_NUMBER(I)
#if defined(WAVE_ROLLER)
       CFF3 = 1.0_SP-TANH((2.0_SP*ZZ(I,K)*GAMW(I))**4)
       CFF3 = CFF3*OROLLER(I)*ROLLA(I)/(WLEN(I)+eps1)*WAVE_C(I)**2
       SXX(I,K) = CFF5*WAVE_NUMBER_X(I)*WAVE_NUMBER_X(I)-CFF4 + &
                  + CFF6 + CFF3*COS_DIR(I)*COS_DIR(I)
       SYY(I,K) = CFF5*WAVE_NUMBER_Y(I)*WAVE_NUMBER_Y(I)-CFF4 + &
                  + CFF6 + CFF3*SIN_DIR(I)*SIN_DIR(I)
       SXY(I,K) = CFF5*WAVE_NUMBER_X(I)*WAVE_NUMBER_Y(I)      + &
                  CFF3*SIN_DIR(I)*COS_DIR(I)
#else
       SXX(I,K) = CFF5*WAVE_NUMBER_X(I)*WAVE_NUMBER_X(I)+CFF6-CFF4
       SYY(I,K) = CFF5*WAVE_NUMBER_Y(I)*WAVE_NUMBER_Y(I)+CFF6-CFF4
       SXY(I,K) = CFF5*WAVE_NUMBER_X(I)*WAVE_NUMBER_Y(I)
#endif       
       TPZDIST(I,K) = FCC*FSS
     END DO  
    ELSE
     DO K=1,KBM1
       FAC2 = ZZ(I,K)
       FCC  = EXP(KD(I)*FAC2)
       FCS  = FCC
       FSC  = FCC
       FSS  = FCC

       CFF1 = WAVE_NUMBER(I)*WAVE_ENERGY(I)
!       CFF4 = CFF1*FCS*FSS
       CFF4 = CFF1*FSC*FSS
       CFF6 = CFF1*FCS*FCC
       CFF5 = CFF1*FCS*FCC*O_WAVE_NUMBER(I)*O_WAVE_NUMBER(I)
#if defined(WAVE_ROLLER)
       CFF3 = 1.0_SP-TANH((2.0_SP*ZZ(I,K)*GAMW(I))**4)
       CFF3 = CFF3*OROLLER(I)*ROLLA(I)/(WLEN(I)+eps1)*WAVE_C(I)**2
       SXX(I,K) = CFF5*WAVE_NUMBER_X(I)*WAVE_NUMBER_X(I)-CFF4 + &
                  + CFF6 + CFF3*COS_DIR(I)*COS_DIR(I)
       SYY(I,K) = CFF5*WAVE_NUMBER_Y(I)*WAVE_NUMBER_Y(I)-CFF4 + &
                  + CFF6 + CFF3*SIN_DIR(I)*SIN_DIR(I)
       SXY(I,K) = CFF5*WAVE_NUMBER_X(I)*WAVE_NUMBER_Y(I)      + &
                  CFF3*SIN_DIR(I)*COS_DIR(I)
#else
       SXX(I,K) = CFF5*WAVE_NUMBER_X(I)*WAVE_NUMBER_X(I)+CFF6-CFF4
       SYY(I,K) = CFF5*WAVE_NUMBER_Y(I)*WAVE_NUMBER_Y(I)+CFF6-CFF4
       SXY(I,K) = CFF5*WAVE_NUMBER_X(I)*WAVE_NUMBER_Y(I)
#endif       
       TPZDIST(I,K) = FCC*FSS
     END DO  
    END IF
   END DO  

# if defined(MULTIPROCESSOR)
   IF(PAR)CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,KB,MYID,NPROCS,SXX,SXY,SYY)
   IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,SXX,SXY,SYY)   !Jianzhong
# endif

   DO I=1,NE
     IA=IEC(I,1)
     IB=IEC(I,2)
     J1=IENODE(I,1)
     J2=IENODE(I,2)
#    if defined (WET_DRY)
     IF(ISWETC(IA) == 1 .OR. ISWETC(IB) == 1)THEN
#    endif
       DO K=1,KBM1
         SXXIJ = 0.5_SP*(SXX(J1,K)+SXX(J2,K))
         SXYIJ = 0.5_SP*(SXY(J1,K)+SXY(J2,K))
         SYYIJ = 0.5_SP*(SYY(J1,K)+SYY(J2,K))
         DIJ = 0.5_SP*(D(J1)*DZ(J1,K)+D(J2)*DZ(J2,K))
#        if defined (SPHERICAL)
         !for spherical coordinator and domain across 360^o          
         XTMP  = VX(J2)*TPI-VX(J1)*TPI
         XTMP1 = VX(J2)-VX(J1)
         IF(XTMP1 >  180.0_SP)THEN
           XTMP = -360.0_SP*TPI+XTMP
         ELSE IF(XTMP1 < -180.0_SP)THEN
           XTMP =  360.0_SP*TPI+XTMP
         END IF  

         EXFLUX       = DIJ*SXXIJ*DLTYC(I)
         PSXXPX(IA,K) = PSXXPX(IA,K) - EXFLUX
         PSXXPX(IB,K) = PSXXPX(IB,K) + EXFLUX

         EXFLUX       = DIJ*SXYIJ*XTMP*COS(DEG2RAD*YC(IA))
         PSXYPY(IA,K) = PSXYPY(IA,K) + EXFLUX
         EXFLUX       = DIJ*SXYIJ*XTMP*COS(DEG2RAD*YC(IB))
         PSXYPY(IB,K) = PSXYPY(IB,K) - EXFLUX

         EXFLUX       = DIJ*SXYIJ*DLTYC(I)
         PSXYPX(IA,K) = PSXYPX(IA,K) - EXFLUX
         PSXYPX(IB,K) = PSXYPX(IB,K) + EXFLUX

         EXFLUX       = DIJ*SYYIJ*XTMP*COS(DEG2RAD*YC(IA))
         PSYYPY(IA,K) = PSYYPY(IA,K) + EXFLUX
         EXFLUX       = DIJ*SYYIJ*XTMP*COS(DEG2RAD*YC(IB))
         PSYYPY(IB,K) = PSYYPY(IB,K) - EXFLUX

#        else
         EXFLUX       = DIJ*SXXIJ*DLTYC(I)
         PSXXPX(IA,K) = PSXXPX(IA,K) - EXFLUX
         PSXXPX(IB,K) = PSXXPX(IB,K) + EXFLUX

         EXFLUX       = DIJ*SXYIJ*DLTXC(I)
         PSXYPY(IA,K) = PSXYPY(IA,K) + EXFLUX
         PSXYPY(IB,K) = PSXYPY(IB,K) - EXFLUX

         EXFLUX       = DIJ*SXYIJ*DLTYC(I)
         PSXYPX(IA,K) = PSXYPX(IA,K) - EXFLUX
         PSXYPX(IB,K) = PSXYPX(IB,K) + EXFLUX

         EXFLUX       = DIJ*SYYIJ*DLTXC(I)
         PSYYPY(IA,K) = PSYYPY(IA,K) + EXFLUX
         PSYYPY(IB,K) = PSYYPY(IB,K) - EXFLUX

#        endif     
       END DO
#    if defined (WET_DRY)
     END IF
#    endif
   END DO
   
   WAVESTRX_3D = 0.0_SP
   WAVESTRY_3D = 0.0_SP
   
   CALL RADIATION_STRESS_Z(WAVE_ENERGY,KD,KDMAX,PSPXPZ,PSPYPZ)

   WAVESTRX_3D = PSXXPX + PSXYPY - PSPXPZ
   WAVESTRY_3D = PSXYPX + PSYYPY - PSPYPZ
#  if defined (WET_DRY)
   DO I = 1,NT
     WAVESTRX_3D(I,:) = WAVESTRX_3D(I,:)*ISWETC(I)
     WAVESTRY_3D(I,:) = WAVESTRY_3D(I,:)*ISWETC(I)
     IF(ISBCE(I) == 2)THEN
       WAVESTRX_3D(I,:) = 0.0_SP
       WAVESTRY_3D(I,:) = 0.0_SP
     END IF
   END DO  
#  endif
!#  if !defined (TWO_D_MODEL)
   WAVESTRX_2D(:) = 0.0_SP; WAVESTRY_2D(:) = 0.0_SP
   DO I = 1,NT
     DO K=1,KBM1
        WAVESTRX_2D(I) = WAVESTRX_2D(I)+WAVESTRX_3D(I,K)
        WAVESTRY_2D(I) = WAVESTRY_2D(I)+WAVESTRY_3D(I,K)
     END DO
   END DO
!#  else
!   CALL RADIATION_STRESS_2D
!#  endif
   WAVESTRX_2D = WAVESTRX_2D*RAMP
   WAVESTRY_2D = WAVESTRY_2D*RAMP
   WAVESTRX_3D = WAVESTRX_3D*RAMP
   WAVESTRY_3D = WAVESTRY_3D*RAMP

#  if defined(MULTIPROCESSOR)
   IF(PAR)CALL AEXCHANGE(EC,MYID,NPROCS,WAVESTRX_3D,WAVESTRY_3D) 
#  endif
!Calculate stokes velocity
   U_STOKES_3D_TMP = 0.0_SP; U_STOKES_3D = 0.0_SP; U_STOKES_2D = 0.0_SP
   V_STOKES_3D_TMP = 0.0_SP; V_STOKES_3D = 0.0_SP; V_STOKES_2D = 0.0_SP
   DO I=1,M
    IF(KD(I) <= KDMAX)THEN
     DO K=1,KBM1
        FAC2 = 1.0_SP+ZZ(I,K)
# if defined (WAVE_ROLLER)
        CFF2=2/WAVE_C(I)*COSH(2*KD(I)*FAC2)/SINH(2*KD(I))*(WAVE_ENERGY(I)+D(I)*GRAV_N(I)*ROLLA(I)/(WLEN(I)+eps1))
# else
        CFF2=2/WAVE_C(I)*COSH(2*KD(I)*FAC2)/SINH(2*KD(I))*WAVE_ENERGY(I)
# endif
        U_STOKES_3D_TMP(I,K)=CFF2*WAVE_NUMBER_X(I)
        V_STOKES_3D_TMP(I,K)=CFF2*WAVE_NUMBER_Y(I)
     ENDDO
    ELSE
     DO K=1,KBM1
        FAC2 = ZZ(I,K)
# if defined (WAVE_ROLLER)
        CFF2=2/WAVE_C(I)*EXP(KD(I)*FAC2)*(WAVE_ENERGY(I)+D(I)*GRAV_N(I)*ROLLA(I)/(WLEN(I)+eps1))
# else
        CFF2=2/WAVE_C(I)*EXP(KD(I)*FAC2)*WAVE_ENERGY(I)
# endif
        U_STOKES_3D_TMP(I,K)=CFF2*WAVE_NUMBER_X(I)
        V_STOKES_3D_TMP(I,K)=CFF2*WAVE_NUMBER_Y(I)
     ENDDO
    END IF 
   ENDDO
   DO I=1,NT
     DO K=1,KBM1
       U_STOKES_3D(I,K)=(U_STOKES_3D_TMP(NV(I,1),K)+U_STOKES_3D_TMP(NV(I,2),K)+U_STOKES_3D_TMP(NV(I,3),K))/3.0_SP
       V_STOKES_3D(I,K)=(V_STOKES_3D_TMP(NV(I,1),K)+V_STOKES_3D_TMP(NV(I,2),K)+V_STOKES_3D_TMP(NV(I,3),K))/3.0_SP
       U_STOKES_2D(I)=U_STOKES_2D(I)+U_STOKES_3D(I,K)*DZ1(I,K)
       V_STOKES_2D(I)=V_STOKES_2D(I)+V_STOKES_3D(I,K)*DZ1(I,K)
     ENDDO
   ENDDO
#  if defined(MULTIPROCESSOR)
   IF(PAR)CALL AEXCHANGE(EC,MYID,NPROCS,U_STOKES_3D,V_STOKES_3D) 
   IF(PAR)CALL AEXCHANGE(EC,MYID,NPROCS,U_STOKES_2D,V_STOKES_2D) 
#  endif

  
   RETURN
   END SUBROUTINE RADIATION_STRESS_3D
!==============================================================================|
!
!==============================================================================|
   SUBROUTINE RADIATION_STRESS_2D       
!------------------------------------------------------------------------------|

   IMPLICIT NONE

   REAL(SP), ALLOCATABLE     :: SXX(:),SXY(:),SYY(:)
   REAL(SP), DIMENSION(0:NT) :: PSXXPX,PSXYPX,PSXYPY,PSYYPY

   REAL(SP), DIMENSION(0:MT) :: WAVE_NUMBER,WAVE_NUMBER_X,WAVE_NUMBER_Y,SIN_DIR,COS_DIR
   REAL(SP), DIMENSION(0:MT) :: WAVE_ENERGY,KD,WAVE_C,WAVE_CGDC
   REAL(SP), DIMENSION(0:MT) :: O_WAVE_NUMBER
   REAL(SP) :: EXFLUX

   INTEGER  :: I,K,IA,IB,J1,J2
   
   REAL(SP) :: CFF1,CFF2,CFF3,SXXIJ,SXYIJ,SYYIJ
   
   REAL(SP) :: XTMP,XTMP1
!==============================================================================|

!---------------Jianzhong----------------------------
   IF(.NOT.ALLOCATED(SXX)) ALLOCATE(SXX(0:MT))
   IF(.NOT.ALLOCATED(SXY)) ALLOCATE(SXY(0:MT))
   IF(.NOT.ALLOCATED(SYY)) ALLOCATE(SYY(0:MT))
!----------------------------------------------------
 
   WAVE_NUMBER = 0.0_SP   ;WAVE_NUMBER_X = 0.0_SP   ;WAVE_NUMBER_Y = 0.0_SP
   WAVE_ENERGY = 0.0_SP   ;KD            = 0.0_SP   ;WAVE_C        = 0.0_SP
   WAVE_CGDC   = 0.0_SP   ;O_WAVE_NUMBER = 0.0_SP

!
!  Compute wave numbers and wave energy.
!
   DO I=1,MT
    WAVE_NUMBER(I) = 2.0_SP*PI/MAX(WLEN(I),WAVE_LENGTH_MIN)
   END DO 
   O_WAVE_NUMBER = 1.0_SP/WAVE_NUMBER
   SIN_DIR       = SIN(DIRDEG1*DEG2RAD)
   COS_DIR       = COS(DIRDEG1*DEG2RAD)
   WAVE_NUMBER_X = WAVE_NUMBER*COS_DIR
   WAVE_NUMBER_Y = WAVE_NUMBER*SIN_DIR
   WAVE_ENERGY   = 0.0625_SP*GRAV_N*HSC1*HSC1
!
!  Compute wave celerity and phase velocity.
!
   DO I=1,MT
     KD(I) = MIN(WAVE_NUMBER(I)*D(I),KDMAX)
!     KD(I) = WAVE_NUMBER(I)*D(I)
   END DO
   
    WAVE_C    = SQRT(GRAV_N*O_WAVE_NUMBER*TANH(KD))
    WAVE_CGDC = 0.5_SP+KD/SINH(2.0_SP*KD)  !Cg/C
    
# if defined(WAVE_ROLLER)
   DO I=1,MT
     ROLLA(I)=0.0424_SP*HSC1(I)*QB1(I)*WLEN(I)
   ENDDO
# endif
!----------INITIALIZE STRESS ARRAY ----------------------------------------------!

   SXX    = 0.0_SP   ;SXY    = 0.0_SP   ;SYY    = 0.0_SP
   PSXXPX = 0.0_SP   ;PSXYPX = 0.0_SP   ;PSXYPY = 0.0_SP   ;PSYYPY = 0.0_SP

   DO I=1,M
     CFF1   = O_WAVE_NUMBER(I)*O_WAVE_NUMBER(I)
     
# if defined(WAVE_ROLLER)
     !ROLLA(I)=0.0424_SP*HSC1(I)*QB1(I)*WLEN(I)
     CFF3   = ROLLA(I)/(WLEN(I)+eps1)*WAVE_C(I)**2
     SXX(I) = WAVE_ENERGY(I)*(WAVE_CGDC(I)*                                       &
              (WAVE_NUMBER_X(I)*WAVE_NUMBER_X(I)*CFF1+1.0_SP)-0.5_SP)         +   &
              CFF3*COS_DIR(I)*COS_DIR(I)
     
     SXY(I) = WAVE_ENERGY(I)*WAVE_CGDC(I)*WAVE_NUMBER_X(I)*WAVE_NUMBER_Y(I)*CFF1+ &
              CFF3*SIN_DIR(I)*COS_DIR(I)

     SYY(I) = WAVE_ENERGY(I)*(WAVE_CGDC(I)*                                       &
              (WAVE_NUMBER_Y(I)*WAVE_NUMBER_Y(I)*CFF1+1.0_SP)-0.5_SP)         +  &
              CFF3*SIN_DIR(I)*SIN_DIR(I)
# else
     SXX(I) = WAVE_ENERGY(I)*(WAVE_CGDC(I)*                                       &
      &         (WAVE_NUMBER_X(I)*WAVE_NUMBER_X(I)*CFF1+1.0_SP)-0.5_SP)
     SXY(I) = WAVE_ENERGY(I)*WAVE_CGDC(I)*WAVE_NUMBER_X(I)*WAVE_NUMBER_Y(I)*CFF1
     SYY(I) = WAVE_ENERGY(I)*(WAVE_CGDC(I)*                                       &
      &         (WAVE_NUMBER_Y(I)*WAVE_NUMBER_Y(I)*CFF1+1.0_SP)-0.5_SP)		
# endif
   !if(I==35) write(789,'(f20.8)')SXX(I)
   END DO  

# if defined(MULTIPROCESSOR)
   IF(PAR)CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,1,MYID,NPROCS,SXX,SXY,SYY)
   IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,SXX,SXY,SYY) !Jianzhong
# endif

   DO I=1,NE
     IA=IEC(I,1)
     IB=IEC(I,2)
     J1=IENODE(I,1)
     J2=IENODE(I,2)
     
#    if defined (WET_DRY)
     IF(ISWETC(IA) == 1 .OR. ISWETC(IB) == 1)THEN
#    endif
 
       SXXIJ=0.5_SP*(SXX(J1)+SXX(J2))
       SXYIJ=0.5_SP*(SXY(J1)+SXY(J2))
       SYYIJ=0.5_SP*(SYY(J1)+SYY(J2))

#      if defined (SPHERICAL)
       !for spherical coordinator and domain across 360^o          
       XTMP  = VX(J2)*TPI-VX(J1)*TPI
       XTMP1 = VX(J2)-VX(J1)
       IF(XTMP1 >  180.0_SP)THEN
         XTMP = -360.0_SP*TPI+XTMP
       ELSE IF(XTMP1 < -180.0_SP)THEN
         XTMP =  360.0_SP*TPI+XTMP
       END IF  

       EXFLUX     = SXXIJ*DLTYC(I)
       PSXXPX(IA) = PSXXPX(IA) - EXFLUX
       PSXXPX(IB) = PSXXPX(IB) + EXFLUX

       EXFLUX     = SXYIJ*XTMP*COS(DEG2RAD*YC(IA))
       PSXYPY(IA) = PSXYPY(IA) + EXFLUX
       EXFLUX     = SXYIJ*XTMP*COS(DEG2RAD*YC(IB))
       PSXYPY(IB) = PSXYPY(IB) - EXFLUX

       EXFLUX     = SXYIJ*DLTYC(I)
       PSXYPX(IA) = PSXYPX(IA) - EXFLUX
       PSXYPX(IB) = PSXYPX(IB) + EXFLUX

       EXFLUX     = SYYIJ*XTMP*COS(DEG2RAD*YC(IA))
       PSYYPY(IA) = PSYYPY(IA) + EXFLUX
       EXFLUX     = SYYIJ*XTMP*COS(DEG2RAD*YC(IB))
       PSYYPY(IB) = PSYYPY(IB) - EXFLUX

#      else
       EXFLUX     = SXXIJ*DLTYC(I)
       PSXXPX(IA) = PSXXPX(IA) - EXFLUX
       PSXXPX(IB) = PSXXPX(IB) + EXFLUX

       EXFLUX     = SXYIJ*DLTXC(I)
       PSXYPY(IA) = PSXYPY(IA) + EXFLUX
       PSXYPY(IB) = PSXYPY(IB) - EXFLUX

       EXFLUX     = SXYIJ*DLTYC(I)
       PSXYPX(IA) = PSXYPX(IA) - EXFLUX
       PSXYPX(IB) = PSXYPX(IB) + EXFLUX

       EXFLUX     = SYYIJ*DLTXC(I)
       PSYYPY(IA) = PSYYPY(IA) + EXFLUX
       PSYYPY(IB) = PSYYPY(IB) - EXFLUX


#      endif     
 
#    if defined (WET_DRY)
     END IF
#    endif
   END DO
   
# if defined (PLBC)
   PSXYPY = 0.0_SP
   PSYYPY = 0.0_SP
# endif


   WAVESTRX_2D = PSXXPX + PSXYPY
   WAVESTRY_2D = PSXYPX + PSYYPY
#  if defined (WET_DRY)
   WAVESTRX_2D = WAVESTRX_2D*ISWETC
   WAVESTRY_2D = WAVESTRY_2D*ISWETC
#  endif
   WAVESTRX_2D = WAVESTRX_2D*RAMP
   WAVESTRY_2D = WAVESTRY_2D*RAMP
   WHERE(ISBCE == 2)
     WAVESTRX_2D = 0.0_SP
     WAVESTRY_2D = 0.0_SP
   END WHERE
#  if defined(MULTIPROCESSOR)
   IF(PAR)CALL AEXCHANGE(EC,MYID,NPROCS,WAVESTRX_2D,WAVESTRY_2D) !Jianzhong
#  endif

!Calculate stokes velocity
   U_STOKES_2D_TMP = 0.0_SP; U_STOKES_2D = 0.0_SP;
   V_STOKES_2D_TMP = 0.0_SP; V_STOKES_2D = 0.0_SP;
   DO I=1,M
        CFF2=1/(WAVE_C(I)*KD(I))*WAVE_ENERGY(I)
# if defined(WAVE_ROLLER)
        CFF3=GRAV_N(I)*ROLLA(I)/(WAVE_C(I)*WLEN(I)+eps1)
        U_STOKES_2D_TMP(I)=CFF2*WAVE_NUMBER_X(I)+CFF3*COS(DIRDEG1(I)*DEG2RAD)
        V_STOKES_2D_TMP(I)=CFF2*WAVE_NUMBER_Y(I)+CFF3*SIN(DIRDEG1(I)*DEG2RAD)
# else
        U_STOKES_2D_TMP(I)=CFF2*WAVE_NUMBER_X(I)
        V_STOKES_2D_TMP(I)=CFF2*WAVE_NUMBER_Y(I)
#endif
   ENDDO
   DO I=1,NT
       U_STOKES_2D(I)=(U_STOKES_2D_TMP(NV(I,1))+U_STOKES_2D_TMP(NV(I,2))+U_STOKES_2D_TMP(NV(I,3)))/3.0_SP
       V_STOKES_2D(I)=(V_STOKES_2D_TMP(NV(I,1))+V_STOKES_2D_TMP(NV(I,2))+V_STOKES_2D_TMP(NV(I,3)))/3.0_SP
   ENDDO
#  if defined(MULTIPROCESSOR)
   IF(PAR)CALL AEXCHANGE(EC,MYID,NPROCS,U_STOKES_2D,V_STOKES_2D)
#  endif

   RETURN
   END SUBROUTINE RADIATION_STRESS_2D	
!==============================================================================|
!
!==============================================================================|

   SUBROUTINE RADIATION_STRESS_Z(WAVE_ENERGY,KD,KDMAX,PSPXPZ,PSPYPZ) 

!==============================================================================|

   IMPLICIT NONE
   REAL(SP), INTENT(IN)  :: WAVE_ENERGY(0:MT),KD(0:MT),KDMAX
   REAL(SP), INTENT(OUT) :: PSPXPZ(0:NT,KB),PSPYPZ(0:NT,KB)
   REAL(SP)              :: SPX(KB),SPY(KB)

   REAL(SP) :: WAVE_ENERGY1(0:NT), KD1(0:NT)
   REAL(SP), DIMENSION(0:MT) :: O_COSH,O_SINH
   INTEGER  :: I,K,J,J1,J2,I1,I2,I3
   REAL(SP) :: FSS1,FCS1,FSC1,FCC1
   REAL(SP) :: CFF1,CFF2,FAC1,FAC2,FAC3
   REAL(SP) :: WEIJ,KDIJ,DIJ,SIJ
#  if defined (SPHERICAL)
   REAL(SP) :: XTMP,XTMP1
#  endif
!==============================================================================|

!----------INITIALIZE ARRAYS---------------------------------------------------!
    CALL N2E2D(WAVE_ENERGY,WAVE_ENERGY1)
    CALL N2E2D(KD,KD1)
    PSPXPZ  = 0.0_SP   ;PSPYPZ  = 0.0_SP  
    O_COSH  = 1.0_SP/COSH(KD)
    O_SINH  = 1.0_SP/SINH(KD)

   DO I = 1, N
     SPX = 0.0_SP; SPY = 0.0_SP
#    if defined (WET_DRY)
     IF(ISWETCT(I)*ISWETC(I) == 1)THEN
#    endif
       I1=NV(I,1);I2=NV(I,2);I3=NV(I,3)
       IF(KD1(I) <= KDMAX)THEN
       
       DO K=1,KBM1
! Calculate some coefficients
         FAC1 = 1.0_SP+Z(I1,K) 
         FAC2 = 1.0_SP+Z(I2,K)
         FAC3 = 1.0_SP+Z(I3,K)
         FCC1 = (COSH(KD(I1)*FAC1)*O_COSH(I1)+COSH(KD(I2)*FAC2)*O_COSH(I2)+COSH(KD(I3)*FAC3)*O_COSH(I3))/3
         FCS1 = (COSH(KD(I1)*FAC1)*O_SINH(I1)+COSH(KD(I2)*FAC2)*O_SINH(I2)+COSH(KD(I3)*FAC3)*O_SINH(I3))/3
         FSC1 = (SINH(KD(I1)*FAC1)*O_COSH(I1)+SINH(KD(I2)*FAC2)*O_COSH(I2)+SINH(KD(I3)*FAC3)*O_COSH(I3))/3
         FSS1 = (SINH(KD(I1)*FAC1)*O_SINH(I1)+SINH(KD(I2)*FAC2)*O_SINH(I2)+SINH(KD(I3)*FAC3)*O_SINH(I3))/3
         CFF1=(FCC1-FSS1)*(FSS1*0.5_SP)
         CFF2=(FCC1-FSS1)*(FCS1*(1+Z1(I,K))*WAVE_ENERGY1(I)-WAVE_ENERGY1(I)*FSS1/TANH(KD1(I)))
         DO J = 1, 3
           J1=J+1-INT((J+1)/4)*3
           J2=J+2-INT((J+2)/4)*3
           WEIJ=0.5_SP*(WAVE_ENERGY(NV(I,J1))+WAVE_ENERGY(NV(I,J2)))*CFF1
           KDIJ=0.5_SP*(KD(NV(I,J1))+KD(NV(I,J2)))*CFF2
           SIJ=WEIJ+KDIJ
#          if defined (SPHERICAL)
           SPX(K)=SPX(K)-DELTUY(I,J)*SIJ
#          else
           SPX(K)=SPX(K)-(VY(NV(I,J2))-VY(NV(I,J1)))*SIJ
#          endif

#          if defined (SPHERICAL)
           XTMP  = VX(NV(I,J2))*TPI-VX(NV(I,J1))*TPI
           XTMP1 = VX(NV(I,J2))-VX(NV(I,J1))
           IF(XTMP1 >  180.0_SP)THEN
             XTMP = -360.0_SP*TPI+XTMP
           ELSE IF(XTMP1 < -180.0_SP)THEN
             XTMP =  360.0_SP*TPI+XTMP
           END IF  

           SPY(K)=SPY(K)+XTMP*COS(DEG2RAD*YC(I))*SIJ
#          else
           SPY(K)=SPY(K)+(VX(NV(I,J2))-VX(NV(I,J1)))*SIJ
#          endif
         END DO
       END DO
       
       ELSE
       
       DO K=1,KBM1
! Calculate some coefficients
         FAC1 = Z(I1,K) 
         FAC2 = Z(I2,K)
         FAC3 = Z(I3,K)
         FCC1 = (EXP(KD(I1)*FAC1)+EXP(KD(I2)*FAC2)+EXP(KD(I3)*FAC3))/3.0_SP
	 FCS1 = FCC1
         FSC1 = FCC1
         FSS1 = FCC1
	 
         CFF1=(FCC1-FSS1)*(FSS1*0.5_SP)
         CFF2=(FCC1-FSS1)*(FCS1*(1+Z1(I,K))*WAVE_ENERGY1(I)-WAVE_ENERGY1(I)*FSS1)
         DO J = 1, 3
           J1=J+1-INT((J+1)/4)*3
           J2=J+2-INT((J+2)/4)*3
           WEIJ=0.5_SP*(WAVE_ENERGY(NV(I,J1))+WAVE_ENERGY(NV(I,J2)))*CFF1
           KDIJ=0.5_SP*(KD(NV(I,J1))+KD(NV(I,J2)))*CFF2
           SIJ=WEIJ+KDIJ
#          if defined (SPHERICAL)
           SPX(K)=SPX(K)-DELTUY(I,J)*SIJ
#          else
           SPX(K)=SPX(K)-(VY(NV(I,J2))-VY(NV(I,J1)))*SIJ
#          endif

#          if defined (SPHERICAL)
           XTMP  = VX(NV(I,J2))*TPI-VX(NV(I,J1))*TPI
           XTMP1 = VX(NV(I,J2))-VX(NV(I,J1))
           IF(XTMP1 >  180.0_SP)THEN
             XTMP = -360.0_SP*TPI+XTMP
           ELSE IF(XTMP1 < -180.0_SP)THEN
             XTMP =  360.0_SP*TPI+XTMP
           END IF  

           SPY(K)=SPY(K)+XTMP*COS(DEG2RAD*YC(I))*SIJ
#          else
           SPY(K)=SPY(K)+(VX(NV(I,J2))-VX(NV(I,J1)))*SIJ
#          endif
         END DO
       END DO
       
       END IF
       
       DO K = 1,KBM1
         PSPXPZ(I,K) = SPX(K)-SPX(K+1) 
	 PSPYPZ(I,K) = SPY(K)-SPY(K+1)
       END DO
#    if defined (WET_DRY)
     END IF
#    endif
   END DO

   RETURN
   END SUBROUTINE RADIATION_STRESS_Z
!==============================================================================|
!
!==============================================================================|
 SUBROUTINE CAL_S(SXX,SXY,SYY)
   IMPLICIT NONE
   REAL(SP), ALLOCATABLE         :: SXX(:,:),SXY(:,:),SYY(:,:)  !Jianzhong

   REAL(SP), DIMENSION(0:MT) :: WAVE_NUMBER,WAVE_NUMBER_X,WAVE_NUMBER_Y ,   &
                                SIN_DIR,COS_DIR
   REAL(SP), DIMENSION(0:MT) :: WAVE_ENERGY,KD,WAVE_C
   REAL(SP), DIMENSION(0:MT) :: O_WAVE_NUMBER
   REAL(SP), DIMENSION(0:MT) :: O_COSH,O_SINH,O_2SINH
   INTEGER  :: I,K,IA,IB,J1,J2
   REAL(SP) :: FSS,FCS,FSC,FCC
   REAL(SP) :: CFF1,CFF2,CFF3,CFF4,CFF5,FAC2
   REAL(SP) :: SXXIJ,SXYIJ,SYYIJ,DIJ,DZD
   
   REAL(SP) :: XTMP,XTMP1
!-------------------------------------------------------------------------------------|

!--------------------Jianzhong---------------------!
   IF(.NOT.ALLOCATED(SXX)) ALLOCATE(SXX(0:MT,KB))
   IF(.NOT.ALLOCATED(SXY)) ALLOCATE(SXY(0:MT,KB))
   IF(.NOT.ALLOCATED(SYY)) ALLOCATE(SYY(0:MT,KB))
!--------------------------------------------------!

 
   WAVE_NUMBER   = 0.0_SP   ;WAVE_NUMBER_X = 0.0_SP   ;WAVE_NUMBER_Y = 0.0_SP
   WAVE_ENERGY   = 0.0_SP   ;KD            = 0.0_SP   ;WAVE_C        = 0.0_SP
   O_COSH        = 0.0_SP   ;O_SINH        = 0.0_SP   ;O_2SINH       = 0.0_SP
   O_WAVE_NUMBER = 0.0_SP
!
!  Compute wave numbers and wave energy.
!
   DO I=1,MT
    WAVE_NUMBER(I) = 2.0_SP*PI/MAX(WLEN(I),WAVE_LENGTH_MIN)
   END DO 
   O_WAVE_NUMBER = 1.0_SP/WAVE_NUMBER
   SIN_DIR       = SIN(DIRDEG1*DEG2RAD)
   COS_DIR       = COS(DIRDEG1*DEG2RAD)
   WAVE_NUMBER_X = WAVE_NUMBER*COS_DIR
   WAVE_NUMBER_Y = WAVE_NUMBER*SIN_DIR
   WAVE_ENERGY   = 0.0625_SP*GRAV_N*HSC1*HSC1
!
!  Compute wave celerity and phase velocity.
!
   DO I=1,MT
      KD(I) = MIN(WAVE_NUMBER(I)*D(I),KDMAX)
   END DO  
   WAVE_C = SQRT(GRAV_N*O_WAVE_NUMBER*TANH(KD))
   O_COSH  = 1.0_SP/COSH(KD)
   O_SINH  = 1.0_SP/SINH(KD)
   O_2SINH = 1.0_SP/SINH(2.0_SP*KD)

   
!----------INITIALIZE STRESS ARRAY ----------------------------------------------!
   SXX    = 0.0_SP   ;SXY    = 0.0_SP   ;SYY    = 0.0_SP
   DO I=1,M
     DO K=1,KB
       FAC2 = 1.0_SP+Z(I,K)
       FCC  = COSH(KD(I)*FAC2)*O_COSH(I)
       FCS  = COSH(KD(I)*FAC2)*O_SINH(I)
       FSC  = SINH(KD(I)*FAC2)*O_COSH(I)
       FSS  = SINH(KD(I)*FAC2)*O_SINH(I)
       CFF1 = WAVE_NUMBER(I)*WAVE_ENERGY(I)
       CFF4 = CFF1*FSC*FSS
       CFF5 = CFF1*FCS*FCC*O_WAVE_NUMBER(I)*O_WAVE_NUMBER(I)
# if defined (WAVE_ROLLER)
       CFF3 = 1.0_SP-TANH((2.0_SP*Z(I,K)*GAMW(I))**4)
       CFF3 = CFF3*OROLLER(I)*ROLLA(I)/(WLEN(I)+eps1)*WAVE_C(I)**2
       SXX(I,K) = CFF5*WAVE_NUMBER_X(I)*WAVE_NUMBER_X(I)-CFF4 + &
                  CFF3*COS_DIR(I)*COS_DIR(I)

       SYY(I,K) = CFF5*WAVE_NUMBER_Y(I)*WAVE_NUMBER_Y(I)-CFF4 + &
                  CFF3*SIN_DIR(I)*SIN_DIR(I)

       SXY(I,K) = CFF5*WAVE_NUMBER_X(I)*WAVE_NUMBER_Y(I)      + &
                  CFF3*SIN_DIR(I)*COS_DIR(I)
# else
       SXX(I,K) = CFF5*WAVE_NUMBER_X(I)*WAVE_NUMBER_X(I)-CFF4
       SYY(I,K) = CFF5*WAVE_NUMBER_Y(I)*WAVE_NUMBER_Y(I)-CFF4
       SXY(I,K) = CFF5*WAVE_NUMBER_X(I)*WAVE_NUMBER_Y(I)
# endif
     END DO  
     DZD = DZ(I,1)*D(I)
!     SXX(I,1) = SXX(I,1)+0.5_SP*WAVE_ENERGY(I)/DZD
!     SYY(I,1) = SYY(I,1)+0.5_SP*WAVE_ENERGY(I)/DZD
     SXX(I,1) = SXX(I,1)+0.5_SP*WAVE_ENERGY(I)/DZD*2.0_SP
     SYY(I,1) = SYY(I,1)+0.5_SP*WAVE_ENERGY(I)/DZD*2.0_SP

   END DO  
# if defined(MULTIPROCESSOR)
   IF(PAR)CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,KB,MYID,NPROCS,SXX,SXY,SYY)
   IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,SXX,SXY,SYY) !Jianzhong
# endif
END SUBROUTINE CAL_S

   SUBROUTINE WAVE_STRESS

   IMPLICIT NONE
   
   INTEGER  :: I,K
   REAL(SP) :: U1,U2,Z0W,Z0T,CDT,CDW 
   REAL(SP) :: WAVE_HEIGHT_E(NT), SIGP_E(NT),SIGP(0:MT)
!      real ageinv(im,jm),usdif(im,jm),vsdif(im,jm)  
!      real cd(im,jm),tpx0(im,jm),tpy0(im,jm),ustw(im,jm),ustp(im,jm) 
!      real ttx0(im,jm),tty0(im,jm)
!      data z0w/1.e-6/,z0t/1.e-6/
      
!   REAL(SP), PARAMETER :: R = 0.0011630_SP, KAPPA = 0.41_SP
    REAL(SP), PARAMETER :: R = 0.0_SP, KAPPA = 0.41_SP  !R=0 for zero wind
   REAL(SP), PARAMETER :: NU = 1.8E-6_SP, CDMAX = 0.005_SP
   
   REAL(SP) :: UVW10,USDIF,VSDIF,AGEINV

!  Calculate inverse of wave age
!   DO I = 1,N
!     if (fsm(i,j).gt.0.) then
!            if(u10(i,j).gt.4.0) then
!     AGEINV(I) = U10(I)*SIGP(I)/GRAV_E(I)       
!            endif
!          endif
!   END DO

   SIGP = 1.0_SP/TPEAK
   DO I = 1,N
     WAVE_HEIGHT_E(I) = (HSC1(NV(I,1))+HSC1(NV(I,2))+HSC1(NV(I,3)))/3.0_SP
     SIGP_E(I) = (SIGP(NV(I,1))+SIGP(NV(I,2))+SIGP(NV(I,3)))/3.0_SP

     UVW10 = SQRT(UW10(I)**2+VW10(I)**2)
     USDIF = UW10(I)-U(I,1)
     VSDIF = VW10(I)-V(I,1)
     U2  = USDIF**2+VSDIF**2
     U1  = SQRT(U2)
     AGEINV = UVW10*SIGP_E(I)/GRAV_E(I)  
          
     Z0W = 1.38E-4_SP*WAVE_HEIGHT_E(I)*(AGEINV)**2.66_SP+1.E-5_SP   !Donelan
     Z0T = 0.18_SP*NU/(USTW(I) +.000001_SP)
     CDW = (KAPPA/LOG(10._SP/Z0W))**2
     CDW = MIN(CDW,CDMAX)
     CDT = (KAPPA/LOG(10._SP/Z0T))**2
! The transition from smooth surface turbulent flow and friction
! drag to a wave surface and form drag is abrupt; this probably
! needs improvement but nevertheless gives a continuos Cd vs U10
! result that agrees with data.
     IF(CDW > CDT)THEN
       CDT = 0.0_SP
     ELSE
       CDW = 0.0_SP
     END IF
!     CD(I) = CDW + CDT
     TPX0(I) = R*CDW*U1*USDIF
     TPY0(I) = R*CDW*U1*VSDIF
     TTX0(I) = R*CDT*U1*USDIF
     TTY0(I) = R*CDT*U1*VSDIF   
     USTW(I) = SQRT(R*(CDW+CDT)*U2)  
     USTP(I) = SQRT(SQRT(TPX0(I)**2+TPY0(I)**2))

!     WUSURF(I) = -TTX0(I)
!     WVSURF(I) = -TTY0(I)
     DO K = 1,KB      !WAVE  Wave pressure momentum transfer
       TPX(I,K) = TPX0(I)*TPZDIST(I,K)
       TPY(I,K) = TPY0(I)*TPZDIST(I,K)
     END DO
   END DO
     
   RETURN
   END SUBROUTINE WAVE_STRESS

#  if defined (VORTEX_FORCE)
!==============================================================================|
! Cite Xia, M., Mao, M., Niu, Q., 2020. Implementation and comparison of the recent
! three-dimensional radiation stress theory and vortex-force formalism in an
! unstructured-grid coastal circulation model. Estuar. Coast. Shelf Sci. 240, 106771.
! https://doi.org/10.1016/j.ecss.2020.106771.
!====================================================================================|  
   SUBROUTINE VORTEX_FORMALISM_3D
   USE MOD_NORTHPOLE
    IMPLICIT NONE
    REAL(SP), DIMENSION(0:MT) :: WAVE_NUMBER,WAVE_NUMBER_X,WAVE_NUMBER_Y,SIN_DIR,COS_DIR
    REAL(SP), DIMENSION(0:MT) :: WAVE_ENERGY,KD,WAVE_C
    REAL(SP), DIMENSION(0:MT) :: DSTP,ODSTP,WAVEAA,XIGMA,OXIGMA,ODISS,GAMR
    REAL(SP), DIMENSION(0:MT) :: KWD,owd
    REAL(SP), DIMENSION(0:MT) :: rukvf2d,rvkvf2d,rubrk2d,rvbrk2d,rurol2d,rvrol2d
    REAL(SP), DIMENSION(0:MT,KB) :: BBXbrk,BBYbrk,BBXrol,BBYrol
    REAL(SP), DIMENSION(0:NT,KB) :: BBXbrk_E,BBYbrk_E,BBXrol_E,BBYrol_E
    REAL(SP), DIMENSION(0:MT,KB) :: BBXwbs,BBYwbs
    REAL(SP), DIMENSION(0:NT,KB) :: BBXwbs_E,BBYwbs_E
    REAL(SP), DIMENSION(0:MT,KB) :: BBXwss,BBYwss
    REAL(SP), DIMENSION(0:NT,KB) :: BBXwss_E,BBYwss_E
    REAL(SP), DIMENSION(0:NT,KB) :: DZ_E
    REAL(SP), DIMENSION(0:NT) :: D_E
    REAL(SP), DIMENSION(0:MT,KB) :: ku,kv,shear1,shear2
    REAL(SP), DIMENSION(0:NT,KB) :: HASt_X,HASt_Y,HVF_X,HVF_Y,WVF_X,WVF_Y,BH_X,BH_Y,shear3,shear4
    REAL(SP), DIMENSION(0:NT,KB) :: DIJ_X,DIJ_Y
    REAL(SP), ALLOCATABLE :: bh(:),qsp(:)
    REAL(SP) :: fac1,fac2,vi1,vi2,wec_alpha
    REAL(SP) :: DDDX,DDDY,DEDX,DEDY,ETF1AA,WW1_STOKES_3D,WW2_STOKES_3D
    REAL(SP), DIMENSION(0:MT) :: O_WAVE_NUMBER
    REAL(SP), DIMENSION(0:MT) :: O_COSH,O_SINH,O_2SINH
    REAL(SP), PARAMETER :: eps1 = 1E-14_SP
    INTEGER :: I,K,IA,IB,J1,J2
    INTEGER :: I1,J3 
    INTEGER :: CNT,CC,JJ
    REAL(SP) :: SIGMAT,EXFLUX
    REAL(SP) :: CFF11,CFF12,CFF13,CFF14,CFF15,CFF16,CFF17,CFF18,CFF19, &
    CFF20,CFF21,CFF22,CFF23,CFF24,CFF25,CFF26,CFF27,CFF28,CFF29,CFF30,CFF31,CFF32
    REAL(SP) :: CFF2,UNODEIJ,VNODEIJ,DZIJ
    ! for roller
    REAL(SP) :: CFF1,CFF3,CFF4
    REAL(SP) :: CFF40,CFF41,CFF42,CFF43
    REAL(SP), PARAMETER :: GAMMAW = 0.31_SP
    REAL(SP) :: USTIJ,VSTIJ
    REAL(SP), PARAMETER :: sinb = 0.1_SP
    REAL(SP), PARAMETER :: KWDmax = 200.0_SP
    REAL(SP), PARAMETER :: awd = 1.0_SP
    REAL(SP), PARAMETER :: ks = 0.045_SP
    REAL(SP), ALLOCATABLE :: GAMW(:),OROLLER(:),ROLLA(:)
    INTEGER :: FLAG_DISSIP
    REAL(SP) :: XTMP,XTMP1,EXFLUXX,EXFLUXY,EXFLUXXA,EXFLUXXB,EXFLUXYA,EXFLUXYB,DIJ,U_STOKES_3DIJ,V_STOKES_3DIJ,W_STOKES_3DIJ,WW_STOKES_3DIJ,bhIJ,qspIJ

   WAVE_NUMBER = 0.0_SP; WAVE_NUMBER_X = 0.0_SP; WAVE_NUMBER_Y = 0.0_SP;
   O_WAVE_NUMBER = 0.0_SP;
   O_COSH = 0.0_SP; O_SINH = 0.0_SP; O_2SINH = 0.0_SP;
   HASt_X =0.0_SP; HASt_Y = 0.0_SP; HVF_X = 0.0_SP; HVF_Y = 0.0_SP
   WVF_X = 0.0_SP; WVF_Y = 0.0_SP; DIJ_X =0.0_SP; DIJ_Y = 0.0_SP
   BH_X = 0.0_SP; BH_Y = 0.0_SP
! Compute total depth
    DSTP = D
  DO I=0,MT
    ODSTP(I) = 1.0_SP/(DSTP(I)+eps1)
  END DO

! Compte wave amplitude (0.5*Hrsm), wave number, instrinsic frequency
   DO I=1,MT
   WAVE_NUMBER(I) = 2.0_SP*PI/MAX(WLEN(I),WAVE_LENGTH_MIN)
   END DO

   O_WAVE_NUMBER = 1.0_SP/WAVE_NUMBER
   SIN_DIR = SIN(DIRDEG1*DEG2RAD)
   COS_DIR = COS(DIRDEG1*DEG2RAD)
   WAVE_NUMBER_X = WAVE_NUMBER*COS_DIR
   WAVE_NUMBER_Y = WAVE_NUMBER*SIN_DIR
   WAVE_ENERGY   = 0.0625_SP*GRAV_N*HSC1*HSC1

   DO I=1,MT
   CFF11 = 0.25_SP*SQRT(2.0_SP)*HSC1(I)
   WAVEAA(I) = CFF11*CFF11
   END DO

# if defined (WET_DRY)
  DO I = 1,MT
      WAVEAA(I) = WAVEAA(I)*ISWETN(I)
  END DO
# endif

  DO I=1,MT
  XIGMA(I) = SQRT(GRAV_N(I)*WAVE_NUMBER(I)*TANH(WAVE_NUMBER(I)*DSTP(I)))
  OXIGMA(I) = 1.0_SP/XIGMA(I)
  END DO

! Computer wave celerity and nonlinear water depth
   DO I=1,MT
    KD(I) = MIN(WAVE_NUMBER(I)*DSTP(I),KDMAX)
  END DO

  WHERE(KD <= KDMAX)
  WAVE_C = SQRT(GRAV_N*O_WAVE_NUMBER*TANH(KD))
    O_COSH  = 1.0_SP/COSH(KD)
    O_SINH  = 1.0_SP/SINH(KD)
    O_2SINH = 1.0_SP/SINH(2.0_SP*KD)
   ELSEWHERE
    WAVE_C = SQRT(GRAV_N*O_WAVE_NUMBER*TANH(KD))
   END WHERE

   ALLOCATE(GAMW(0:MT));             GAMW        = 0.0_SP
   ALLOCATE(OROLLER(0:MT));          OROLLER     = 0.0_SP
   ALLOCATE(ROLLA(0:MT));            ROLLA      = 0.0_SP

#if defined(WAVE_ROLLER)
   OROLLER = 0.0_SP;GAMW = 0.0_SP; ROLLA = 0.0_SP
   DO I=1,MT
     GAMW(I) = MIN(D(I)/(HSC1(I)+eps1),5.0_SP)
     DO K=1,KBM1
        OROLLER(I)=OROLLER(I)+D(I)*DZ(I,K)*(1.0_SP-TANH((2.0_SP*ZZ(I,K)*GAMW(I))**4))
     ENDDO
   OROLLER(I)=1.0_SP/(OROLLER(I)+eps1)
   ROLLA(I)=0.0424*HSC1(I)*QB1(I)*WLEN(I)
   ENDDO
#endif

! Compute metrics for vertical dissipation distribution
  DO I=1,MT
   ODISS(I) = 0.0_SP
    GAMR(I) = MIN(0.707_SP*DSTP(I)/(1.25_SP*HSC1(I)+eps1),1.0_SP)
  DO K=KBM1,1,-1
     CFF12 = (ZZ(I,K)-Z(I,KB))*D(I)*ODSTP(I)*GAMR(I)
     ODISS(I)=ODISS(I)+DZ(I,K)*D(I)*COSH(2.0_SP*PI*CFF12)
  END DO
  ODISS(I)=1.0_SP/(ODISS(I)+eps1)
  END DO

! ------------------------------------------------------------------
! Compute Bernoulli's Head (BH), represented as the
! symbol cursive K in Eqn. 5
! -------------------------------------------------------------------
  DO I=1,M
    DO K=1,KBM1
     UNODE(I,K) = 0.0_SP
     VNODE(I,K) = 0.0_SP
     CNT = 0
     DO JJ=1,NTVE(I)
        CNT =CNT + 1
        CC = NBVE(I,JJ)
        UNODE(I,K) = UNODE(I,K) + U(CC,K)
        VNODE(I,K) = VNODE(I,K) + V(CC,K)
     ENDDO
     UNODE(I,K) = UNODE(I,K) / CNT
     VNODE(I,K) = VNODE(I,K) / CNT
    ENDDO
   ENDDO
!------------------------------------------------------------------------
   ALLOCATE(bh(0:MT)); bh = 0.0_SP
   ALLOCATE(qsp(0:MT)); qsp = 0.0_SP
! Compute k*u and k*v at rho pts.
 DO I=1,MT
  DO K=KBM1,1,-1
 ku(I,K) = UNODE(I,K)*WAVE_NUMBER_X(I)
 kv(I,K) = VNODE(I,K)*WAVE_NUMBER_Y(I)
 END DO
  END DO

! Compute gradients for shear d(k_x*u)/dz at cell faces
 DO I=1,MT
 DO K=KBM1,2,-1
 CFF13 = 1.0_SP/(D(I)*(ZZ(I,K-1)-ZZ(I,K)))
  shear1(I,K) = CFF13*((ku(I,K-1)-ku(I,K))+(kv(I,K-1)-kv(I,K)))
 END DO
   END DO

! Applying boundary conditions
 DO I=1,MT
  shear1(I,KB)= shear1(I,KB-1)
  shear1(I,1) = shear1(I,2)
 END DO

! Compute second derivative d(k_x*u)^2/dz^2
 DO I=1,MT
 DO K=KBM1,1,-1
! CFF14 = 1.0_SP/(D(I)*(ZZ(I,K)-ZZ(I,K+1)))
  CFF14 = 1.0_SP/(D(I)*(Z(I,K)-Z(I,K+1)))
 shear2(I,K) = CFF14*(shear1(I,K)-shear1(I,K+1))
 END DO
 END DO

! Perform integration on the shear2 term calculated above
 DO I=1,MT
 CFF15 = D(I)*DZ(I,KBM1)
 fac2 = D(I)*(ZZ(I,1)-ZZ(I,KBM1))*ODSTP(I)
 vi1 = CFF15*shear2(I,KBM1)*SINH(2.0_SP*KD(I)*fac2)  

 DO K=KBM1-1,1,-1
 CFF16 = D(I)*DZ(I,K)
 fac2 = D(I)*(ZZ(I,1)-ZZ(I,K))*ODSTP(I)
 vi1 = vi1+CFF16*shear2(I,K)*SINH(2.0_SP*KD(I)*fac2)
 END DO
! Calculate Bernoulli's Head at the surface
 bh(I)=0.25_SP*vi1*XIGMA(I)*WAVEAA(I)*O_SINH(I)*O_SINH(I)/(WAVE_NUMBER(I))
 END DO
!---------------------------------------------------------------------------
! Compute quasi-static components of pressure terms (i.e. Cursive P)
! in Eqn. 9
!---------------------------------------------------------------------------

  DO I=1,MT
 fac2 = D(I)*(ZZ(I,KBM1)-Z(I,KBM1+1))*ODSTP(I)
 CFF17 = D(I)*DZ(I,KBM1)
 vi2 = CFF17*shear2(I,KBM1)*COSH(2.0_SP*KD(I)*fac2)
  DO K=KBM1-1,1,-1
  fac2 = D(I)*(ZZ(I,K)-Z(I,KBM1+1))*ODSTP(I)
  CFF18 = D(I)*DZ(I,K)
  vi2 = vi2+CFF18*shear2(I,K)*COSH(2.0_SP*KD(I)*fac2)

 END DO
 CFF19 = 0.5_SP*WAVEAA(I)*OXIGMA(I)
 CFF20 = TANH(KD(I))*O_2SINH(I)

 qsp(I)=CFF19*CFF20*(-shear1(I,1)+COSH(2*KD(I))*shear1(I,KBM1+1)+ &
         vi2)-2.0_SP*WAVE_NUMBER(I)*TANH(KD(I))*(ku(I,1)+kv(I,1))
 qsp(I)=0.0_SP 
 END DO

! ---------------------------------------------------------------------
! Compute Vertical Vortex Force terms, denoted as K in Eqn. 5
! This next section needs to be re-writeen for piped J. 
! ---------------------------------------------------------------------
! Calculate Stokes velocity first
   U_STOKES_3D_TMP = 0.0_SP; U_STOKES_3D = 0.0_SP; U_STOKES_2D = 0.0_SP
   V_STOKES_3D_TMP = 0.0_SP; V_STOKES_3D = 0.0_SP; V_STOKES_2D = 0.0_SP
   W_STOKES_3D = 0.0_SP; WW_STOKES_3D = 0.0_SP

DO I=1,M
     DO K=1,KBM1
        FAC2 = 1.0_SP+ZZ(I,K)
        IF (KD(I)<=1.0_SP) THEN ! Xia et al. (2020)
# if defined (WAVE_ROLLER)
        CFF2=2/WAVE_C(I)*COSH(2*KD(I)*FAC2)/SINH(2*KD(I))*(WAVE_ENERGY(I)+D(I)*GRAV_N(I)*ROLLA(I)/(WLEN(I)+eps1))
# else
        CFF2=2/WAVE_C(I)*COSH(2*KD(I)*FAC2)/SINH(2*KD(I))*WAVE_ENERGY(I)
# endif
        ELSE IF (KD(I)>=10.0_SP) THEN
# if defined (WAVE_ROLLER)
        CFF2=2/WAVE_C(I)*EXP(2.0_SP*ZZ(I,K)*KD(I))*(WAVE_ENERGY(I)+D(I)*GRAV_N(I)*ROLLA(I)/(WLEN(I)+eps1))
# else
        CFF2=2/WAVE_C(I)*EXP(2.0_SP*ZZ(I,K)*KD(I))*WAVE_ENERGY(I)
# endif
       ELSE
# if defined (WAVE_ROLLER)
        CFF2=((10.0_SP-KD(I))/9.0_SP)*(2/WAVE_C(I)*COSH(2*KD(I)*FAC2)/SINH(2*KD(I))*(WAVE_ENERGY(I)+D(I)*GRAV_N(I)*ROLLA(I)/(WLEN(I)+eps1))) + &
             ((KD(I)-1.0_SP)/9.0_SP)*(2/WAVE_C(I)*EXP(2.0_SP*ZZ(I,K)*KD(I))*(WAVE_ENERGY(I)+D(I)*GRAV_N(I)*ROLLA(I)/(WLEN(I)+eps1)))
# else
        CFF2=((10.0_SP-KD(I))/9.0_SP)*(2/WAVE_C(I)*COSH(2*KD(I)*FAC2)/SINH(2*KD(I))*WAVE_ENERGY(I)) + &
             ((KD(I)-1.0_SP)/9.0_SP)*(2/WAVE_C(I)*EXP(2.0_SP*ZZ(I,K)*KD(I))*WAVE_ENERGY(I))
# endif
       END IF
        U_STOKES_3D_TMP(I,K)=CFF2*WAVE_NUMBER_X(I)
        V_STOKES_3D_TMP(I,K)=CFF2*WAVE_NUMBER_Y(I)
     END DO
   END DO

  DO I=1,NT
     DO K=KBM1,1,-1
       U_STOKES_3D(I,K)=(U_STOKES_3D_TMP(NV(I,1),K)+U_STOKES_3D_TMP(NV(I,2),K)+U_STOKES_3D_TMP(NV(I,3),K))/3.0_SP
       V_STOKES_3D(I,K)=(V_STOKES_3D_TMP(NV(I,1),K)+V_STOKES_3D_TMP(NV(I,2),K)+V_STOKES_3D_TMP(NV(I,3),K))/3.0_SP

#  if defined (WET_DRY)
     U_STOKES_3D(I,K) = U_STOKES_3D(I,K)*ISWETC(I)
     V_STOKES_3D(I,K) = V_STOKES_3D(I,K)*ISWETC(I)
     IF(ISBCE(I) == 2)THEN
       U_STOKES_3D(I,K) = 0.0_SP
       V_STOKES_3D(I,K) = 0.0_SP
     END IF
#  endif

       U_STOKES_2D(I)=U_STOKES_2D(I)+U_STOKES_3D(I,K)*DZ1(I,K)
       V_STOKES_2D(I)=V_STOKES_2D(I)+V_STOKES_3D(I,K)*DZ1(I,K)
     END DO
   END DO

# if defined(MULTIPROCESSOR)
   IF(PAR)CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,KB,MYID,NPROCS,U_STOKES_3D_TMP,V_STOKES_3D_TMP) 
   IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,U_STOKES_3D_TMP,V_STOKES_3D_TMP)
# endif

  EXFLUXXA = 0.0_SP;
  EXFLUXXB = 0.0_SP;
  EXFLUXYA = 0.0_SP;
  EXFLUXYB = 0.0_SP;
  
   DO I=1,NE
     IA=IEC(I,1)
     IB=IEC(I,2)
     J1=IENODE(I,1)
     J2=IENODE(I,2)

#    if defined (WET_DRY)
     IF(ISWETC(IA) == 1 .OR. ISWETC(IB) == 1)THEN
#    endif 

    DO K=KBM1,1,-1
     U_STOKES_3DIJ=0.5_SP*(U_STOKES_3D_TMP(J1,K)+U_STOKES_3D_TMP(J2,K))
     V_STOKES_3DIJ=0.5_SP*(V_STOKES_3D_TMP(J1,K)+V_STOKES_3D_TMP(J2,K))
     UNODEIJ=0.5_SP*(UNODE(J1,K)+UNODE(J2,K))
     VNODEIJ=0.5_SP*(VNODE(J1,K)+VNODE(J2,K))
     DIJ=0.5_SP*(D(J1)*DZ(J1,K)+D(J2)*DZ(J2,K))
     DZIJ=0.5_SP*(DZ(J1,K)+DZ(J2,K))
#      if defined (SPHERICAL)
       !for spherical coordinator and domain across 360^o          
       XTMP  = VX(J2)*TPI-VX(J1)*TPI
       XTMP1 = VX(J2)-VX(J1)
       IF(XTMP1 >  180.0_SP)THEN
         XTMP = -360.0_SP*TPI+XTMP
       ELSE IF(XTMP1 < -180.0_SP)THEN
         XTMP =  360.0_SP*TPI+XTMP
       END IF
      EXFLUXXA=DIJ*U_STOKES_3DIJ*DLTYC(I)-DIJ*V_STOKES_3DIJ*XTMP*COS(DEG2RAD*YC(IA))
      EXFLUXXB=DIJ*U_STOKES_3DIJ*DLTYC(I)-DIJ*V_STOKES_3DIJ*XTMP*COS(DEG2RAD*YC(IB))
      EXFLUXYA=DIJ*U_STOKES_3DIJ*DLTYC(I)-DIJ*V_STOKES_3DIJ*XTMP*COS(DEG2RAD*YC(IA))
      EXFLUXYB=DIJ*U_STOKES_3DIJ*DLTYC(I)-DIJ*V_STOKES_3DIJ*XTMP*COS(DEG2RAD*YC(IB))

#      else
      EXFLUXXA=DIJ*U_STOKES_3DIJ*DLTYC(I)-DIJ*V_STOKES_3DIJ*DLTXC(I)
      EXFLUXXB=DIJ*U_STOKES_3DIJ*DLTYC(I)-DIJ*V_STOKES_3DIJ*DLTXC(I)
      EXFLUXYA=DIJ*U_STOKES_3DIJ*DLTYC(I)-DIJ*V_STOKES_3DIJ*DLTXC(I)
      EXFLUXYB=DIJ*U_STOKES_3DIJ*DLTYC(I)-DIJ*V_STOKES_3DIJ*DLTXC(I)
# endif
      HASt_X(IA,K)=HASt_X(IA,K)-EXFLUXXA
      HASt_X(IB,K)=HASt_X(IB,K)+EXFLUXXB
      HASt_Y(IA,K)=HASt_Y(IA,K)-EXFLUXYA
      HASt_Y(IB,K)=HASt_Y(IB,K)+EXFLUXYB
 END DO

#    if defined (WET_DRY)
     END IF
#    endif 
   END DO

      WAVESTRX_3D = 0.0_SP
      WAVESTRY_3D = 0.0_SP

  DO I=1,NT
   DO K=KBM1,1,-1
     WAVESTRX_3D(I,K) = WAVESTRX_3D(I,K) + HASt_X(I,K)*U(I,K)
     WAVESTRY_3D(I,K) = WAVESTRY_3D(I,K) + HASt_Y(I,K)*V(I,K)
  ENDDO
 ENDDO

#  if defined (WET_DRY)
   DO I = 1,NT
     WAVESTRX_3D(I,:) = WAVESTRX_3D(I,:)*ISWETC(I)
     WAVESTRY_3D(I,:) = WAVESTRY_3D(I,:)*ISWETC(I)
     IF(ISBCE(I) == 2)THEN
       WAVESTRX_3D(I,:) = 0.0_SP
       WAVESTRY_3D(I,:) = 0.0_SP
     END IF
   END DO
#  endif
!3D1
!-----------------Calculate the omega stokes velocity-------------
!-----------------------------------------INITIALIZE FLUX-----------------------------!
 XSTFLUX = 0.0_SP
!----------------------------------ACCUMULATE FLUX-------------------------------------!
  DO I=1,NCV
 I1=NTRG(I)
 IA=NIEC(I,1)
 IB=NIEC(I,2)

 DO K=KBM1,1,-1

 DIJ = DT1(I1)*DZ1(I1,K)
 USTIJ = U_STOKES_3D(I1,K)
 VSTIJ = V_STOKES_3D(I1,K)
 EXFLUX = DIJ*(-USTIJ*DLTYE(I)+VSTIJ*DLTXE(I))

# if defined (WET_DRY)
 EXFLUX =  EXFLUX*ISWETCT(I1)*ISWETC(I1)
# endif

 XSTFLUX(IA,K)=XSTFLUX(IA,K)-EXFLUX
 XSTFLUX(IB,K)=XSTFLUX(IB,K)+EXFLUX

 ENDDO

 ENDDO

# if defined (SPHERICAL)
 CALL VERTVL_EDGE_XY(XSTFLUX,0.0_SP)
#  endif
 !--------------------------NULIFY BOUNDARY FLUX-----------------------!
# if !defined (MEAN_FLOW)
 DO I=1,M
  DO K=KBM1,1,-1
   IF(ISONB(I) == 2) XSTFLUX(I,K)=0.0_SP
  ENDDO
 ENDDO
# endif

# if defined (ONE_D_MODEL)
 XSTFLUX = 0.0_SP
# endif

!CLEAR OLD VALUES
 W_STOKES_3D_TMP = 0.0_SP

!---------------------CALCULATE STOKES OMEGA --------------------------!
 DO I = 1,M
# if defined (WET_DRY)
 IF(ISWETNT(I)*ISWETN(I) == 1) THEN
# endif

 DO K = KBM1,1,-1
 W_STOKES_3D_TMP(I,K+1) = W_STOKES_3D_TMP(I,K)+XSTFLUX(I,K)/ART1(I)
 ENDDO

# if defined (WET_DRY)
 ELSE
 DO K = KBM1,1,-1
 W_STOKES_3D_TMP(I,K+1)=0.0_SP
 ENDDO

 ENDIF
# endif
 ENDDO

 DO I=1,M
 W_STOKES_3D_TMP(I,KBM1) = 0.0_SP
 W_STOKES_3D_TMP(I,1) = 0.0_SP
 END DO

# if defined (MULTIPROCESSOR)
 IF(PAR) CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,KB,MYID,NPROCS,W_STOKES_3D_TMP)
 IF(PAR) CALL AEXCHANGE(NC,MYID,NPROCS,W_STOKES_3D_TMP)
# endif
!---------------------TRANSFER STOKES OMEGA TO FACE CENTER----------------------!
 DO I = 1,N
 DO K = KB,1,-1
 W_STOKES_3D(I,K) = ONE_THIRD*(W_STOKES_3D_TMP(NV(I,1),K)+W_STOKES_3D_TMP(NV(I,2),K)+W_STOKES_3D_TMP(NV(I,3),K))
 ENDDO
 ENDDO

# if defined (WET_DRY)
 DO I = 1, N
 DO K = KB,1,-1
 W_STOKES_3D(I,K) = FLOAT(ISWETC(I))*W_STOKES_3D(I,K)
 ENDDO
 ENDDO
# endif

  DO I=1,NT 
#  if defined (WET_DRY)
    IF(ISWETC(I) == 1)THEN
#  endif
     J1=NV(I,1)
     J2=NV(I,2)
     J3=NV(I,3)
     DDDX=AWX(I,1) * D(J1)+AWX(I,2) * D(J2)+AWX(I,3)*D(J3)
     DDDY=AWY(I,1) * D(J1)+AWY(I,2) * D(J2)+AWY(I,3)*D(J3)
     DEDX=AWX(I,1)*ELF(J1)+AWX(I,2)*ELF(J2)+AWX(I,3)*ELF(J3)
     DEDY=AWY(I,1)*ELF(J1)+AWY(I,2)*ELF(J2)+AWY(I,3)*ELF(J3)
     ETF1AA=ONE_THIRD*(EL(NV(I,1))+EL(NV(I,2))+EL(NV(I,3)))

  DO K = KBM1,1,-1
 WW1_STOKES_3D = 0.5_SP*(W_STOKES_3D(I,K)+W_STOKES_3D(I,K+1))+ &
  U_STOKES_3D(I,K)*(ZZ1(I,K)*DDDX+DEDX)+ &
  V_STOKES_3D(I,K)*(ZZ1(I,K)*DDDY+DEDY)
 WW2_STOKES_3D = 0.0_SP
 WW_STOKES_3D(I,K) = WW1_STOKES_3D+WW2_STOKES_3D
 ENDDO

#  if defined (WET_DRY)
    ELSE
     DO K=1,KBM1
      WW_STOKES_3D(I,K)=0.0_SP !Xia et al. (2020)
     END DO
    END IF
#  endif
   END DO

# if defined(MULTIPROCESSOR)
   IF(PAR)CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,KB,MYID,NPROCS,UNODE,VNODE)
   IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,UNODE,VNODE)
# endif

    DO I=1,NE
     IA=IEC(I,1)
     IB=IEC(I,2)
     J1=IENODE(I,1)
     J2=IENODE(I,2)
#    if defined (WET_DRY)
     IF(ISWETC(IA) == 1 .OR. ISWETC(IB) == 1)THEN
#    endif 

    DO K=KBM1,1,-1

    UNODEIJ=0.5_SP*(UNODE(J1,K)+UNODE(J2,K))
    VNODEIJ=0.5_SP*(VNODE(J1,K)+VNODE(J2,K))
    DIJ=0.5_SP*(D(J1)*DZ(J1,K)+D(J2)*DZ(J2,K))
    DZIJ=0.5_SP*(DZ(J1,K)+DZ(J2,K))

#    if defined (SPHERICAL)
       !for spherical coordinator and domain across 360^o          
       XTMP  = VX(J2)*TPI-VX(J1)*TPI
       XTMP1 = VX(J2)-VX(J1)
       IF(XTMP1 >  180.0_SP)THEN
         XTMP = -360.0_SP*TPI+XTMP
       ELSE IF(XTMP1 < -180.0_SP)THEN
         XTMP =  360.0_SP*TPI+XTMP
       END IF

      EXFLUXXA=-(VNODEIJ*DLTYC(I)+UNODEIJ*XTMP*(COS(DEG2RAD*YC(IA))))*DZIJ
      EXFLUXXB=-(VNODEIJ*DLTYC(I)+UNODEIJ*XTMP*(COS(DEG2RAD*YC(IB))))*DZIJ
      EXFLUXYA= (VNODEIJ*DLTYC(I)+UNODEIJ*XTMP*(COS(DEG2RAD*YC(IA))))*DZIJ
      EXFLUXYB= (VNODEIJ*DLTYC(I)+UNODEIJ*XTMP*(COS(DEG2RAD*YC(IB))))*DZIJ
# else
      EXFLUXXA=-(VNODEIJ*DLTYC(I)+UNODEIJ*DLTXC(I))*DZIJ
      EXFLUXXB=-(VNODEIJ*DLTYC(I)+UNODEIJ*DLTXC(I))*DZIJ
      EXFLUXYA= (VNODEIJ*DLTYC(I)+UNODEIJ*DLTXC(I))*DZIJ
      EXFLUXYB= (VNODEIJ*DLTYC(I)+UNODEIJ*DLTXC(I))*DZIJ
# endif       
      HVF_X(IA,K)=HVF_X(IA,K)-EXFLUXXA
      HVF_X(IB,K)=HVF_X(IB,K)+EXFLUXXB
      HVF_Y(IA,K)=HVF_Y(IA,K)-EXFLUXYA
      HVF_Y(IB,K)=HVF_Y(IB,K)+EXFLUXYB
  END DO

#    if defined (WET_DRY)
     ENDIF
#    endif 
  END DO

  CALL N2E3D(DZ,DZ_E)
  CALL N2E2D(D,D_E)
  DO I = 1, NT
     DO K = KBM1,1,-1
     WAVESTRX_3D(I,K) = WAVESTRX_3D(I,K) + HVF_X(I,K)*V_STOKES_3D(I,K)*D_E(I)
     WAVESTRY_3D(I,K) = WAVESTRY_3D(I,K) + HVF_Y(I,K)*U_STOKES_3D(I,K)*D_E(I)
  ENDDO
 ENDDO
!3D2
!---------------------------------------------------
!Compute non conservative wave acceleration terms here
!----------------------------------------------------
 DO I=1,MT
  rubrk2d(I)=0.0_SP
  rvbrk2d(I)=0.0_SP
#if defined(WAVE_ROLLER)
  rurol2d(I)=0.0_SP
  rvrol2d(I)=0.0_SP
#endif
  END DO
!NOTE:CALCULATE THE DISWCP AND DISBRK TERMS BASED ON ROMS wec_dissip.F
!read in the flag of dissipation equations (FLAG_DISSIP)
!FLAG_DISSIP=1:Thorton and Guza (1986)
!FLAG_DISSIP=2:Church and Thorton (1993)
!Default:FLAG_DISSIP=1:Thorton and Guza (1986)
!NOTE:NEED TO CHECK, AS DISSIP_wWCAP = 0.0 IN THIS METHOD 

! ININTH('FLAG_DISSIP',FLAG_DISSIP,'STA',1)

 FLAG_DISSIP = 1
  IF(FLAG_DISSIP<1 .OR. FLAG_DISSIP>2) CALL FATAL_ERROR ("FLAG_DISSIP should be defined as integer 1 or 2")
  IF(FLAG_DISSIP==1) THEN
! Calculate based on Thorton and Guza (1986)
  fac1 = 3.0_SP*9.81_SP**SQRT(PI)/16.0_SP
  DO I = 1, MT
  CFF1=0.707_SP*HSC1(I)
  SIGMAT = MIN(1.0_SP/(TPEAK(I)+eps1),1.0_SP)
  CFF2 = 1.0_SP/((GAMMAW**4)*(DSTP(I)**5))
  DISSIP_BREAK(I)=0.2621_SP*fac1*SIGMAT*(CFF1**7.0_SP)*CFF2
  ENDDO
 ELSEIF(FLAG_DISSIP==2) THEN
! Calculate based on Church and Thorton (1993)
 fac1 = 3.0_SP*SQRT(PI)/16.0_SP
 DO I = 1, MT
 CFF1 = 0.707_SP*HSC1(I) 
 CFF2 = 1.0_SP/(GAMMAW*DSTP(I))
 CFF3 = 1.0_SP+TANH(8.0_SP*((CFF1*CFF2)-1.0_SP))
 CFF4 = 1.0_SP-(1.0_SP+(CFF1*CFF2)**2.0_SP)**(-2.5_SP)
 SIGMAT = MIN(1.0_SP/(TPEAK(I)+eps1),1.0_SP)

 DISSIP_BREAK(I) = (0.2621_SP/DSTP(I))*fac1*SIGMAT*(CFF1**3.0_SP)*CFF3*CFF4
 DISSIP_WCAP(I) = 0.0_SP

 ENDDO
 ENDIF
!NOTE:CALCULATE THE DISSIP_ROLLER BASED ON ROMS wec_roller.F
# if defined(WAVE_ROLLER)

 DO I = 1, MT
  CFF1 = 0.0424_SP*HSC1(I)*QB1(I)
  ROLLA_VF(I) = CFF1*WAVE_C(I)*WAVE_C(I)*OXIGMA(I)
  DISSIP_ROLLER(I) = GRAV_N(I)*sinb*ROLLA_VF(I)*XIGMA(I)/WAVE_C(I)
 ENDDO

# endif

 CALL INREAL('BRE_WEC_ALPHA',wec_alpha,'STA',0.0)
 
 BBXbrk = 0.0_SP; BBYbrk = 0.0_SP
 BBXrol = 0.0_SP; BBYrol = 0.0_SP 
 DO I=1,MT 
  DO K=KBM1,1,-1
  fac2=D(I)*(ZZ(I,K)-Z(I,KBM1+1))*ODSTP(I)
  CFF28 = fac2*GAMR(I)
  CFF29 = COSH(2.0_SP*PI*CFF28)*ODISS(I)
  CFF30 = DISSIP_WCAP(I)+DISSIP_BREAK(I)
  CFF31 = ((1.0_SP-wec_alpha)*CFF30+wec_alpha*DISSIP_WCAP(I))*OXIGMA(I)
  BBXbrk(I,K) = CFF29*CFF31*WAVE_NUMBER_X(I)
  BBYbrk(I,K) = CFF29*CFF31*WAVE_NUMBER_Y(I)

#  if defined(WAVE_ROLLER)
 CFF32 = DISSIP_ROLLER(I)*OXIGMA(I)
 BBXrol(I,K) = CFF29*CFF32*WAVE_NUMBER_X(I)
 BBYrol(I,K) = CFF29*CFF32*WAVE_NUMBER_Y(I)
#  endif
   END DO 
  END DO

 BBXbrk = BBXbrk + BBXrol
 BBYbrk = BBYbrk + BBYrol

 CALL N2E3D(BBXbrk,BBXbrk_E)
 CALL N2E3D(BBYbrk,BBYbrk_E)
 CALL N2E3D(BBXrol,BBXrol_E)
 CALL N2E3D(BBYrol,BBYrol_E)
! Compute contribution to U-momentum  and V-momentum
 DO I=1,NT
 DO K=KBM1,1,-1
  WAVESTRX_3D(I,K)=WAVESTRX_3D(I,K)-BBXbrk_E(I,K)*D_E(I)*DZ1(I,K)*ART(I)
  WAVESTRY_3D(I,K)=WAVESTRY_3D(I,K)-BBYbrk_E(I,K)*D_E(I)*DZ1(I,K)*ART(I)
 END DO 
 END DO

! Wave-induced bottom streaming based on Xu and Bowen (1994) and Uchiyama et al. (2010) using Type I f function
! CFF40: wave orbital velocity CFF41: decay length
! CFF43: wave friction factor
 DO I=1,MT
   CFF40=ABS(0.25_SP*SQRT(2.0_SP)*XIGMA(I)*HSC1(I)/SINH(KD(I)+eps1))
   CFF41=awd*0.09_SP*ks*(CFF40/(ks*XIGMA(I)))**0.82_SP
   KWD(I)=MAX(MIN(DSTP(I)/(CFF41+eps1),KWDmax),eps1)
 
  owd(I) = 0.0_SP

 DO K=KBM1,1,-1
 CFF42=D(I)*(ZZ(I,K)-Z(I,KBM1+1))*ODSTP(I)
 owd(I)=owd(I)+DZ(I,K)*D(I)*(1-TANH((KWD(I)*CFF42)**4.0_SP))
 END DO

 owd(I)=1.0_SP/(owd(I)+eps1)

   CFF43=MIN(1.39_SP*(XIGMA(I)*(ks/30.0_SP)/CFF40)**0.52_SP,0.2_SP)
   DISSIP_FRIC(I)=(0.5_SP/SQRT(PI))*CFF43*(CFF40**3.0_SP)
 END DO

 BBXwbs = 0.0_SP; BBYwbs = 0.0_SP

 DO I=1,MT
  DO K=KBM1,1,-1 
  fac2=D(I)*(ZZ(I,1)-ZZ(I,K))*ODSTP(I)
  CFF42=1-TANH((KWD(I)*fac2)**4.0_SP)

  CFF43=DISSIP_FRIC(I)*OXIGMA(I)
  BBXwbs(I,K)=CFF42*CFF43*WAVE_NUMBER_X(I)*owd(I)
  BBYwbs(I,K)=CFF42*CFF43*WAVE_NUMBER_Y(I)*owd(I)
  END DO
 END DO 

 CALL N2E3D(BBXwbs,BBXwbs_E)
 CALL N2E3D(BBYwbs,BBYwbs_E)
! Compute contribution to U-momentum  and V-momentum
 DO I=1,NT
 DO K=KBM1,1,-1
  WAVESTRX_3D(I,K)=WAVESTRX_3D(I,K)-BBXwbs_E(I,K)*D_E(I)*DZ1(I,K)*ART(I)
  WAVESTRY_3D(I,K)=WAVESTRY_3D(I,K)-BBYwbs_E(I,K)*D_E(I)*DZ1(I,K)*ART(I)
 END DO
 END DO

! Wave-induced surface streaming effects based on Xu and Bowen (1994)

 BBXwss = 0.0_SP; BBYwss = 0.0_SP

 DO I=1,MT
 BBXwss(I,1)=(0.25_SP/TANH(KD(I)+eps1))*KM(I,1)*(HSC1(I)**2.0_SP)*XIGMA(I)*WAVE_NUMBER(I)*WAVE_NUMBER_X(I)/(2.0_SP*DZ(I,1)*D(I))
 BBYwss(I,1)=(0.25_SP/TANH(KD(I)+eps1))*KM(I,1)*(HSC1(I)**2.0_SP)*XIGMA(I)*WAVE_NUMBER(I)*WAVE_NUMBER_Y(I)/(2.0_SP*DZ(I,1)*D(I))
 END DO

 CALL N2E3D(BBXwss,BBXwss_E)
 CALL N2E3D(BBYwss,BBYwss_E)

 DO I=1,NT
  WAVESTRX_3D(I,1)=WAVESTRX_3D(I,1)-BBXwbs_E(I,1)*D_E(I)*DZ1(I,1)*ART(I)
  WAVESTRY_3D(I,1)=WAVESTRY_3D(I,1)-BBYwbs_E(I,1)*D_E(I)*DZ1(I,1)*ART(I)
 END DO

!3D3

  DO I=1,NT
  DO K=KBM1,1,-1
    WAVESTRX_3D(I,K)=WAVESTRX_3D(I,K)-COR(I)*V_STOKES_3D(I,K)*D_E(I)*DZ1(I,K)*ART(I)
    WAVESTRY_3D(I,K)=WAVESTRY_3D(I,K)+COR(I)*U_STOKES_3D(I,K)*D_E(I)*DZ1(I,K)*ART(I)
  END DO
 END DO

!3D4

 ! Compute  shear3 and shear4
 DO I=1,NT
 DO K=KBM1,2,-1
  shear3(I,K) = (U(I,K-1)*W_STOKES_3D(I,K-1)-U(I,K)*W_STOKES_3D(I,K))/DZZ1(I,K-1)
  shear4(I,K) = (V(I,K-1)*W_STOKES_3D(I,K-1)-V(I,K)*W_STOKES_3D(I,K))/DZZ1(I,K-1)
 END DO
   END DO

! Applying boundary conditions
 DO I=1,NT
  shear3(I,KB)= shear3(I,KB-1)
  shear3(I,1) = shear3(I,2)
  shear4(I,KB)= shear4(I,KB-1)
  shear4(I,1) = shear4(I,2)
 END DO

 DO I=1,NT
 DO K=KBM1,1,-1
 WAVESTRX_3D(I,K) = WAVESTRX_3D(I,K) + shear3(I,K)*DZ1(I,K)*ART(I)
 WAVESTRY_3D(I,K) = WAVESTRY_3D(I,K) + shear4(I,K)*DZ1(I,K)*ART(I)
 END DO
 END DO

! 3D5

# if defined(MULTIPROCESSOR)
   IF(PAR)CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,1,MYID,NPROCS,bh,qsp)
   IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,bh,qsp)   
# endif

 qsp = 0.0_SP

    DO I=1,NE
     IA=IEC(I,1)
     IB=IEC(I,2)
     J1=IENODE(I,1)
     J2=IENODE(I,2)
#    if defined (WET_DRY)
     IF(ISWETC(IA) == 1 .OR. ISWETC(IB) == 1)THEN
#    endif 

    DO K=KBM1,1,-1

     bhIJ=0.5_SP*(bh(J1)+bh(J2))
     qspIJ=0.5_SP*(qsp(J1)+qsp(J2))
     DIJ=0.5_SP*(D(J1)*DZ(J1,K)+D(J2)*DZ(J2,K))
     DZIJ=0.5_SP*(DZ(J1,K)+DZ(J2,K)) 
#      if defined (SPHERICAL)
       !for spherical coordinator and domain across 360^o          
       XTMP  = VX(J2)*TPI-VX(J1)*TPI
       XTMP1 = VX(J2)-VX(J1)
       IF(XTMP1 >  180.0_SP)THEN
         XTMP = -360.0_SP*TPI+XTMP
       ELSE IF(XTMP1 < -180.0_SP)THEN
         XTMP =  360.0_SP*TPI+XTMP
       END IF

      EXFLUXXA=DZIJ*(bhIJ+qspIJ)*DLTYC(I)
      EXFLUXXB=DZIJ*(bhIJ+qspIJ)*DLTYC(I)
      EXFLUXYA=-DZIJ*(bhIJ+qspIJ)*XTMP*(COS(DEG2RAD*YC(IA)))
      EXFLUXYB=-DZIJ*(bhIJ+qspIJ)*XTMP*(COS(DEG2RAD*YC(IB)))
# else
      EXFLUXXA=DZIJ*(bhIJ+qspIJ)*DLTYC(I)
      EXFLUXXB=DZIJ*(bhIJ+qspIJ)*DLTYC(I)
      EXFLUXYA=-DZIJ*(bhIJ+qspIJ)*DLTXC(I)
      EXFLUXYB=-DZIJ*(bhIJ+qspIJ)*DLTXC(I)
# endif
      BH_X(IA,K)=BH_X(IA,K)-EXFLUXXA
      BH_X(IB,K)=BH_X(IB,K)+EXFLUXXB                       
      BH_Y(IA,K)=BH_Y(IA,K)-EXFLUXYA
      BH_Y(IB,K)=BH_Y(IB,K)+EXFLUXYB
  END DO
#    if defined (WET_DRY)
     ENDIF
#    endif 
  END DO

 DO I=1,NT
 DO K=KBM1,1,-1
 WAVESTRX_3D(I,K)=WAVESTRX_3D(I,K)+BH_X(I,K)*D_E(I)
 WAVESTRY_3D(I,K)=WAVESTRY_3D(I,K)+BH_Y(I,K)*D_E(I)
 ENDDO
 ENDDO
 
! Applying boundary conditions
 DO I=1,NT
  WAVESTRX_3D(I,KB)= WAVESTRX_3D(I,KB-1)
  WAVESTRX_3D(I,1) = WAVESTRX_3D(I,2)
  WAVESTRY_3D(I,KB)= WAVESTRY_3D(I,KB-1)
  WAVESTRY_3D(I,1) = WAVESTRY_3D(I,2) 
END DO

#  if defined (WET_DRY)
   DO I = 1,NT
     WAVESTRX_3D(I,:) = WAVESTRX_3D(I,:)*ISWETC(I)
     WAVESTRY_3D(I,:) = WAVESTRY_3D(I,:)*ISWETC(I)
     IF(ISBCE(I) == 2)THEN
       WAVESTRX_3D(I,:) = 0.0_SP
       WAVESTRY_3D(I,:) = 0.0_SP
     END IF
   END DO
#  endif

!#  if !defined (TWO_D_MODEL)
   WAVESTRX_2D(:) = 0.0_SP; WAVESTRY_2D(:) = 0.0_SP
   DO I = 1,NT
     DO K=1,KBM1
        WAVESTRX_2D(I) = WAVESTRX_2D(I)+WAVESTRX_3D(I,K)
        WAVESTRY_2D(I) = WAVESTRY_2D(I)+WAVESTRY_3D(I,K)
     END DO
   END DO
!#  else
!   CALL RADIATION_STRESS_2D
!#  endif
   WAVESTRX_2D = WAVESTRX_2D*RAMP
   WAVESTRY_2D = WAVESTRY_2D*RAMP  
   WAVESTRX_3D = WAVESTRX_3D*RAMP
   WAVESTRY_3D = WAVESTRY_3D*RAMP

! 3D6
#  if defined(MULTIPROCESSOR)
   IF(PAR)CALL AEXCHANGE(EC,MYID,NPROCS,WAVESTRX_3D,WAVESTRY_3D)
#  endif

 RETURN
 END SUBROUTINE VORTEX_FORMALISM_3D
# endif
! ------------------------------------------------------------------
!
!====================================================================================|  
   SUBROUTINE CURRENT2WAVE       

   IMPLICIT NONE

   REAL(SP), DIMENSION(0:MT) :: WAVE_NUMBER,KD
   REAL(SP), DIMENSION(0:MT) :: O_COSH,O_SINH
   INTEGER  :: I,K,CNT,CC,JJ
   
   REAL(SP) :: FSS,FCS,FSC,FCC
   REAL(SP) :: CFF1,CFF2,CFF3,CFF4,CFF5,FAC2
!-------------------------------------------------------------------------------------|
#  if !defined (TWO_D_MODEL)
   DO I=1,M
    DO K=1,KBM1
     UNODE(I,K) = 0.0_SP
     VNODE(I,K) = 0.0_SP
     CNT = 0
     DO JJ=1,NTVE(I)
        CNT =CNT + 1
        CC = NBVE(I,JJ)
        UNODE(I,K) = UNODE(I,K) + U(CC,K)
        VNODE(I,K) = VNODE(I,K) + V(CC,K)
     END DO
     UNODE(I,K) = UNODE(I,K) / CNT
     VNODE(I,K) = VNODE(I,K) / CNT
    END DO
   END DO
!-------------------------------------------------------------------------------------|
 
   WAVE_NUMBER = 0.0_SP
   KD          = 0.0_SP
   O_COSH      = 0.0_SP
   O_SINH      = 0.0_SP
   
!
!  Compute wave numbers and wave energy.
!
   DO I=1,MT
    WAVE_NUMBER(I) = 2.0_SP*PI/MAX(WLEN(I),WAVE_LENGTH_MIN)
   END DO 
!
!  Compute wave celerity and phase velocity.
!
   DO I=1,MT
!     KD(I) = MIN(WAVE_NUMBER(I)*D(I)+eps1,KDMAX)
     KD(I) = WAVE_NUMBER(I)*D(I)+eps1
   END DO  

   WHERE(KD <= KDMAX)
    O_COSH  = 1.0_SP/COSH(KD)
    O_SINH  = 1.0_SP/SINH(KD)
   END WHERE 
   
!----------INITIALIZE STRESS ARRAY ----------------------------------------------!

   UDOP = 0.0_SP;  VDOP = 0.0_SP

   DO I=1,M
    IF(KD(I) <= KDMAX)THEN
     DO K=1,KBM1
       FAC2 = KD(I)*(1.0_SP+ZZ(I,K))
       FCC  = COSH(FAC2)*O_COSH(I)
       FCS  = COSH(FAC2)*O_SINH(I)
       FSC  = SINH(FAC2)*O_COSH(I)
       FSS  = SINH(FAC2)*O_SINH(I)
       
       CFF1 = 0.5*(FCS*FCC+FSS*FSC)
       CFF4 = FCS*FSS
       
       UDOP(I) = UDOP(I) + UNODE(I,K)*(CFF1+CFF4)*DZ(I,K)
       VDOP(I) = VDOP(I) + VNODE(I,K)*(CFF1+CFF4)*DZ(I,K)
     END DO  
     UDOP(I) = UDOP(I)*KD(I)
     VDOP(I) = VDOP(I)*KD(I)
    ELSE
     DO K=1,KBM1
       FAC2 = KD(I)*ZZ(I,K)
       FCC  = EXP(FAC2)
       FCS  = FCC
       FSC  = FCC
       FSS  = FCC
       
       CFF1 = 0.5*(FCS*FCC+FSS*FSC)
       CFF4 = FCS*FSS
       
       UDOP(I) = UDOP(I) + UNODE(I,K)*(CFF1+CFF4)*DZ(I,K)
       VDOP(I) = VDOP(I) + VNODE(I,K)*(CFF1+CFF4)*DZ(I,K)
     END DO  
     UDOP(I) = UDOP(I)*KDMAX
     VDOP(I) = VDOP(I)*KDMAX
    END IF 
   END DO  
#  else
   CALL E2N2D(UA,UDOP)
   CALL E2N2D(VA,VDOP)
#  endif   
# if defined(MULTIPROCESSOR)
   IF(PAR)CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,1,MYID,NPROCS,UDOP,VDOP)
   IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,UDOP,VDOP)
# endif

   RETURN
   END SUBROUTINE CURRENT2WAVE
!==============================================================================|
!==============================================================================|

END MODULE MOD_WAVE_CURRENT_INTERACTION
#  endif   
 
