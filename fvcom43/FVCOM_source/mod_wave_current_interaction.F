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
! end
   REAL(SP), PARAMETER :: KDMAX = 5.0_SP
   REAL(SP), PARAMETER :: WAVE_LENGTH_MIN = 0.01_SP
   REAL(SP), PARAMETER :: eps1 = 1E-14_SP

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
   ALLOCATE(U_STOKES_2D_TMP(0:MT));  U_STOKES_2D_TMP = 0.0_SP
   ALLOCATE(V_STOKES_2D_TMP(0:MT));  V_STOKES_2D_TMP = 0.0_SP
   ALLOCATE(U_STOKES_3D_TMP(0:MT,KB));   U_STOKES_3D_TMP = 0.0_SP
   ALLOCATE(V_STOKES_3D_TMP(0:MT,KB));   V_STOKES_3D_TMP = 0.0_SP
   ALLOCATE(GAMW(0:MT));             GAMW        = 0.0_SP
   ALLOCATE(OROLLER(0:MT));          OROLLER     = 0.0_SP
   ALLOCATE(ROLLA(0:MT));            ROLLA      = 0.0_SP
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
     ENDDO
     OROLLER(I)=1.0_SP/(OROLLER(I)+eps1)
     ROLLA(I)=0.0424*HSC1(I)*QB1(I)*WLEN(I)
   ENDDO
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
     ENDDO
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
   ! need add MPI part
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
   ! need add MPI part
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
!==============================================================================|
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
     ENDDO
     UNODE(I,K) = UNODE(I,K) / CNT
     VNODE(I,K) = VNODE(I,K) / CNT
    ENDDO
   ENDDO

   
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
 
