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

SUBROUTINE INTERNAL_STEP
  USE MOD_NESTING
  USE ALL_VARS
  USE MOD_UTILS
  USE MOD_OBCS
  USE MOD_TIME
  USE EQS_OF_STATE
  USE MOD_WD
  USE MOD_ASSIM
  USE MOD_PAR
  USE MOD_ICE
  USE MOD_ICE2D
  USE MOD_NORTHPOLE
! lwang added for SST mld assimilation
  USE MOD_MLD_RHO
! lwang
# if defined (SEMI_IMPLICIT)
  USE MOD_SEMI_IMPLICIT
# endif
# if defined (BALANCE_2D)
  USE MOD_BALANCE_2D
# endif
# if defined (WATER_QUALITY)
  USE MOD_WQM
# endif
# if defined (GOTM)
  USE MOD_GOTM
# endif  
# if defined (BioGen)
  USE MOD_BIO_3D
# endif      
# if defined (NH)
  USE NON_HYDRO
# endif
# if defined (SEDIMENT)
# if defined (ORIG_SED)
  USE MOD_SED
# elif defined (CSTMS_SED)
  USE MOD_SED_CSTMS
# endif
# endif
# if defined (VEGETATION)
  USE MOD_VEGETATION
# endif
# if defined (DYE_RELEASE)
  USE MOD_DYE
# endif
# if defined (MEAN_FLOW)
  USE MOD_MEANFLOW
  USE MOD_OBCS2
  USE MOD_OBCS3
# endif  
# if defined (WAVE_CURRENT_INTERACTION)
  USE TIMECOMM
  USE SWCOMM3
  USE VARS_WAVE
  USE MOD_WAVE_CURRENT_INTERACTION
# endif 

# if defined (PLBC)
  USE MOD_PERIODIC_LBC
# endif

# if defined (THIN_DAM)
  USE MOD_DAM
# endif
# if defined (ENKF)
  use MOD_ENKFA , only : RECALC_RHO_MEAN_ENKF
  USE ENKFVAL,  ONLY : ENKF_NENS
  USE MOD_ENKF,   ONLY : ENKF_memberCONTR
# endif

# if defined (PWP)
  USE MOD_PWP
# endif
  IMPLICIT NONE
  INTEGER :: I,J,K
  REAL(SP) :: UTMP,VTMP
# if defined (ICE)
    real(SP) :: DUVI,ice_ocnDX,ice_ocnDY
# endif
# if defined (SEDIMENT)
  REAL(DP) :: DPDAYS
  REAL(SP) :: SPDAYS
# endif
  integer :: i1,i2

  REAL(SP), ALLOCATABLE :: TB1(:,:),SB1(:,:),Q2B(:,:),Q2LB(:,:),UB(:,:),VB(:,:)
  REAL(SP), ALLOCATABLE :: DYEB(:,:)
  REAL(SP) :: RCOFT
  INTEGER :: ITER
!------------------------------------------------------------------------------|

  if(dbg_set(dbg_sbr)) write(ipt,*)&
       & "Start: internal_step"

!----SET RAMP FACTOR TO EASE SPINUP--------------------------------------------!
  RAMP = 0.0_SP
  IF(IRAMP /= 0) THEN
     RAMP=TANH(real(IINT,sp)/real(IRAMP,sp))
  ELSE
     RAMP = 1.0_SP
  END IF

!----OFFLINE SEDIMENT UPDATE
# if defined (OFFLINE_SEDIMENT) || (OFFLINE_BIOLOGY)
  !----SPECIFY THE SURFACE FORCING OF INTERNAL MODES-----------------------------!
# if defined (GCN)
  CALL BCOND_GCN(8,0)
#   else
  CALL BCOND_GCY(8,0)
# endif

# if defined (WAVE_CURRENT_INTERACTION)
# if defined (WAVE_OFFLINE)
  HSC1      = WHS
  DIRDEG1   = WDIR
  TPEAK     = WPER
  WLEN      = WLENGTH
  Pwave_bot = WPER_BOT
  Ub_swan   = WUB_BOT
  Dwave = DIRDEG1*pi/180.0
# endif
# endif

  U   = OFFLINE_U
  V   = OFFLINE_V 
  WTS = OFFLINE_W
  S1  = OFFLINE_S1
  T1  = OFFLINE_T1
  KH  = OFFLINE_KH
  EL  = OFFLINE_EL
# if defined (GOTM)
  TEPS = OFFLINE_TEPS
# else
  Q2  = OFFLINE_Q2
  Q2L = OFFLINE_Q2L
# endif
  
# if defined (WET_DRY)
  iswetn = INT(OFFLINE_WN)
  iswetnt = iswetn
  iswetc = INT(OFFLINE_WC)
  iswetct = iswetc
# endif
  D = H + EL
  DT = D
  DTFA = D
  CALL N2E2D(D,D1)
  CALL N2E2D(DT,DT1)

# if defined(MULTIPROCESSOR)
  IF(PAR)THEN
    CALL AEXCHANGE(EC,MYID,NPROCS,U,V)
    CALL AEXCHANGE(NC,MYID,NPROCS,S1,T1)
    CALL AEXCHANGE(NC,MYID,NPROCS,KH,WTS)
    CALL AEXCHANGE(NC,MYID,NPROCS,EL,D,DT)
    CALL AEXCHANGE(NC,MYID,NPROCS,DTFA)
    CALL AEXCHANGE(EC,MYID,NPROCS,D1,DT1)
# if defined (WET_DRY)
    CALL AEXCHANGE(NC,MYID,NPROCS,iswetn,iswetnt)
    CALL AEXCHANGE(EC,MYID,NPROCS,iswetc,iswetct)
# endif
  END IF
# endif

!# if defined (WET_DRY)
!  US=U
!  VS=V
!# endif

# if defined (THIN_DAM)
  CALL GET_KDAM
  CALL GET_KDAM_INTERNAL
# endif



!==============================================================================!
!    SEDIMENT MODEL SECTION  
!      Advance Sed Model (Erode/Deposit/Advect/Diffuse)      
!      Change bathymetry (if MORPHO_MODEL=T)  
!        Note:  morph array returned from sed:  Deposition:  morph > 0
!        0 < morpho_factor < 1  
!      morpho_model and morpho_factor are stored and set in mod_sed.F              
!==============================================================================!
# if defined (SEDIMENT)

  CALL BOTTOM_ROUGHNESS

  DPDAYS = days(inttime)
  SPDAYS = DPDAYS
  IF(SEDIMENT_MODEL)CALL ADVANCE_SED(DTI,spdays,taubm_n(1:m))

  IF(MORPHO_MODEL)THEN
    if(mod(sed_its,morpho_incr)==0 .and. sed_its > morpho_strt)then
       CALL UPDATE_DEPTH
    endif
  ENDIF

# endif

# if defined (BioGen)
  IF(BIOLOGICAL_MODEL)CALL BIO_3D1D                           
# endif


  RETURN

# endif
!----loop end for offline sediment



!----SET UP WATER QUALITY MODEL COEFFICIENTS-----------------------------------!
# if defined (WATER_QUALITY)
  IF(WATER_QUALITY_MODEL)THEN
     TIME_R=MOD(IINT*DTI/3600.0_SP-14.743_SP, 24.0_SP)+6.0_SP
     CALL WQMCONST
  END IF
# endif


  
# if !defined (TWO_D_MODEL)

# if !defined(SEMI_IMPLICIT) && !defined(NH)
  !----ADJUST CONSISTENCY BETWEEN 3-D VELOCITY AND VERT-AVERAGED VELOCITIES------!
#   if !defined (ONE_D_MODEL)
    CALL ADJUST2D3D(1)
#   endif

  !----SPECIFY THE SOLID BOUNDARY CONDITIONS OF U&V INTERNAL MODES---------------!
#   if defined (GCN)
    CALL BCOND_GCN(5,0)
#   else
    CALL BCOND_GCY(5,0)
#   endif

#   if defined(MULTIPROCESSOR)
    IF(PAR)CALL AEXCHANGE(EC,MYID,NPROCS,U,V)
#   endif
# endif
! end defined semi-implicit && nh

  !----SPECIFY THE SURFACE FORCING OF INTERNAL MODES-----------------------------!
# if defined (GCN)
  CALL BCOND_GCN(8,0)
#   else
  CALL BCOND_GCY(8,0)
# endif

  ! New Open Boundary Condition ----3
# if defined (MEAN_FLOW)
  CALL BCOND_MEANFLOW
  CALL BCOND_TIDE_3D
  !     CALL BCOND_NG_3D       !change BKI
  !     CALL BCOND_NG_2D       !change BKI
  CALL BCOND_BKI_3D(1)  
  !     CALL BCOND_BKI_2D(4)

  CALL FLUX_OBC3D_2
  CALL FLUX_OBC2D
# endif

# endif    
  !end !defined (TWO_D_MODEL)
  
  !----SPECIFY THE BOTTOM ROUGHNESS AND CALCULATE THE BOTTOM STRESSES------------!
  CALL BOTTOM_ROUGHNESS

# if defined (ICE)

  !----ITERATE THE ICE DYNAMIC/THERMODYNAMIC MODEL ------------!
  
  IF(ICE_MODEL) THEN
     
     IF(dbg_set(dbg_log))   WRITE(IPT,*)IINT, 'Calculating ICE'
     !!  the thermodynamics part are modified from CICE
     !!  the dynamics and advection part are use FVCOM-ice
     
     
     !==================================================================
     ! UPDATE ICE FORCING
     !==================================================================
     IF(DBG_SET(DBG_IO)) WRITE(IPT,*) "START UPDATE ICE FORCING"
     CALL UPDATE_ICE(now=IntTime,SAT=T_AIR,SWV=DSW_AIR,SPQ=QA_AIR,CLD=CLOUD)
     IF(DBG_SET(DBG_IO)) WRITE(IPT,*) "FINISHED UPDATE ICE FORCING"
     
     CALL ICE_1D_2D
     
     IF(DBG_SET(DBG_IO)) WRITE(IPT,*) "FINISHED ICE MODEL"
     
  END IF

     !----------------------------------------------------------------------------------!
     !----------------------------------------------------------------------------------!
     !----------------------------------------------------------------------------------!
#endif

!==============================================================================!
!  CALCULATE DISPERSION (GX/GY) AND BAROCLINIC PRESSURE GRADIENT TERMS         !
!==============================================================================!
# if !defined (TWO_D_MODEL)

# if !defined (SEMI_IMPLICIT)
#   if defined (GCN)
    CALL ADVECTION_EDGE_GCN(ADVX,ADVY)          !Calculate 3-D Adv/Diff       !
#   if defined (THIN_DAM)
    IF(PAR)CALL AEXCHANGE(EC,MYID,NPROCS,ADVX,ADVY) 
#   endif

#   else
    CALL ADVECTION_EDGE_GCY(ADVX,ADVY)          !Calculate 3-D Adv/Diff       !
#   endif
# endif
# if defined (ENKF)
# if !defined (DATA_ASSIM)
  IF(RECALCULATE_RHO_MEAN) THEN
     IF(RECALC_RHO_MEAN .LT. IntTime)THEN
         IF (RECALC_RHO_MEAN_ENKF(ENKF_memberCONTR)==0) THEN
            CALL RHO_PMEAN 
            RECALC_RHO_MEAN_ENKF(ENKF_memberCONTR)=1
            IF (RECALC_RHO_MEAN_ENKF(ENKF_NENS)==1) THEN
                RECALC_RHO_MEAN_ENKF=0
                RECALC_RHO_MEAN = RECALC_RHO_MEAN + DELT_RHO_MEAN
            END IF
         END IF
     END IF
  END IF
# endif
# else
  IF(RECALCULATE_RHO_MEAN) THEN
     IF(RECALC_RHO_MEAN .LT. IntTime)THEN
        RECALC_RHO_MEAN = RECALC_RHO_MEAN + DELT_RHO_MEAN
        CALL RHO_PMEAN    
     END IF
  END IF

# endif

  IF(.NOT. BAROTROPIC)THEN                    !Barotropic Flow ?            !
     SELECT CASE(BAROCLINIC_PRESSURE_GRADIENT)
     CASE ("sigma levels")
        CALL BAROPG      !Sigma Level Pressure Gradient!
     CASE('z coordinates')
        CALL PHY_BAROPG  !Z Level Pressure Gradient    !
     CASE DEFAULT
        CALL FATAL_ERROR("UNKNOW BAROCLINIC PRESURE GRADIENT TYPE",&
             & TRIM(BAROCLINIC_PRESSURE_GRADIENT))
     END SELECT

  END IF                                      !                             !

# if !defined (SEMI_IMPLICIT)

  ADX2D = 0.0_SP ; ADY2D = 0.0_SP             !Initialize GX/GY Terms       !
  DRX2D = 0.0_SP ; DRY2D = 0.0_SP             !Initialize BCPG for Ext Mode !

  DO K=1,KBM1
     DO I=1, N
        ADX2D(I)=ADX2D(I)+ADVX(I,K)   !*DZ1(I,K)
        ADY2D(I)=ADY2D(I)+ADVY(I,K)   !*DZ1(I,K)
        DRX2D(I)=DRX2D(I)+DRHOX(I,K)  !*DZ1(I,K)
        DRY2D(I)=DRY2D(I)+DRHOY(I,K)  !*DZ1(I,K)
     END DO
  END DO

# if defined (GCN)
  CALL ADVAVE_EDGE_GCN(ADVUA,ADVVA)           !Compute Ext Mode Adv/Diff
# else
  CALL ADVAVE_EDGE_GCY(ADVUA,ADVVA)           !Compute Ext Mode Adv/Diff
# endif
  ADX2D = ADX2D - ADVUA                       !Subtract to Form GX
  ADY2D = ADY2D - ADVVA                       !Subtract to Form GY

  !----INITIALIZE ARRAYS USED TO CALCULATE AVERAGE UA/E  OVER EXTERNAL STEPS-----!
  UARD = 0.0_SP
  VARD = 0.0_SP
  EGF  = 0.0_SP

# if defined (EQUI_TIDE)
  EGF_EQI = 0.0_SP
# endif
# if defined (ATMO_TIDE)
  EGF_ATMO = 0.0_SP
# endif       
# if defined (AIR_PRESSURE)
  EGF_AIR = 0.0_SP
# endif

!!# if defined (WET_DRY)
!!  UARDS = 0.0_SP
!!  VARDS = 0.0_SP
!!# endif

  IF(IOBCN > 0) THEN
     UARD_OBCN(1:IOBCN)=0.0_SP
  END IF
# endif
! end defined semi-implicit
# endif    
! end defined (TWO_D_MODEL)

# if defined (BALANCE_2D)
  ADVUA2    = 0.0_SP
  ADVVA2    = 0.0_SP
  ADFX2     = 0.0_SP
  ADFY2     = 0.0_SP
  DRX2D2    = 0.0_SP
  DRY2D2    = 0.0_SP
  CORX2     = 0.0_SP
  CORY2     = 0.0_SP
  PSTX2     = 0.0_SP
  PSTY2     = 0.0_SP
  ADX2D2    = 0.0_SP
  ADY2D2    = 0.0_SP
  WUSURBF2  = 0.0_SP
  WVSURBF2  = 0.0_SP 
  DUDT2     = 0.0_SP 
  DVDT2     = 0.0_SP

  DIVX2D2   = ZERO 
  DIVY2D2   = ZERO 
  DEDT2     = ZERO  
# endif

  ! New Open Boundary Condition ----4
# if defined (MEAN_FLOW)
  CALL ZERO_OBC3
# endif   

# if defined (WAVE_CURRENT_INTERACTION)
# if defined (DATA_ASSIM)
  IF((ASSIM_FLAG == 0 .AND. (.NOT. SST_ASSIM .AND. .NOT. SSTGRD_ASSIM) .AND. .NOT. SSHGRD_ASSIM .AND. .NOT. TSGRD_ASSIM) .OR. &
     (ASSIM_FLAG == 1 .AND. (SST_ASSIM .OR. SSTGRD_ASSIM .OR. SSHGRD_ASSIM .OR. TSGRD_ASSIM)))THEN
# endif     
# if defined (WAVE_OFFLINE)
  HSC1      = WHS
  DIRDEG1   = WDIR
  TPEAK     = WPER
  WLEN      = WLENGTH
  Pwave_bot = WPER_BOT
  Ub_swan   = WUB_BOT
# else
  IF(IINT == 1 .OR. MOD(IINT,NINT(DTW/DTI)) == 0.0_SP)THEN
    ITW = ITW + 1
    CALL SWMAIN_LOOP
  END IF
# endif
# if defined (SEDIMENT) && (WAVE_CURRENT_INTERACTION)
  Dwave = DIRDEG1*pi/180.0
# endif

# if !defined (WAVE_ONLY)
# if !defined (TWO_D_MODEL)
# if !defined (VORTEX_FORCE)
  CALL RADIATION_STRESS_3D
# else
!Xia et al. (2020)
  CALL VORTEX_FORMALISM_3D
# endif
# else
  CALL RADIATION_STRESS_2D
# endif
# endif
# if defined (DATA_ASSIM)
  END IF
# endif     
  
# endif

# if defined (WAVE_ONLY)
  goto 5501
# endif

# if !defined(ONE_D_MODEL) && !defined (SEMI_IMPLICIT)
  !==============================================================================!
  !  LOOP OVER EXTERNAL TIME STEPS                                               !
  !==============================================================================!
  DO IEXT=1,ISPLIT

     IF (DBG_SET(DBG_SBRIO)) WRITE(IPT,*) "/// EXT SETP: ",IEXT

     ExtTime = ExtTime + IMDTE

     CALL EXTERNAL_STEP
     
     
  END DO
# endif

# if defined (THIN_DAM)
  CALL GET_KDAM   
  CALL GET_KDAM_INTERNAL
# endif

!  new update  !  jqi ggao 0730/2007  for E-P

#  if defined (ONE_D_MODEL)
!   EL = EL - DTI*(QEVAP-QPREC)*ROFVROS
   EL = EL + DTI*(QEVAP+QPREC)*ROFVROS
   !!INTERPOLATE DEPTH FROM NODE-BASED TO ELEMENT-BASED VALUES
   CALL N2E2D(EL,EL1)

   !UPDATE DEPTH AND VERTICALLY AVERAGED VELOCITY FIELD
   D   = H + EL
   D1  = H1 + EL1
   DTFA = D

#  endif
!  new update  !  jqi ggao 0730/2007

# if defined (SEMI_IMPLICIT)

!! David for VISIT
!!=======================================================
#  if defined (VISIT)
   Call VisitCheck
#  endif
!!=======================================================

#  if defined (EQUI_TIDE)
   CALL ELEVATION_EQUI
#  endif
#  if defined (ATMO_TIDE)
   CALL ELEVATION_ATMO
#  endif
#  if defined (AIR_PRESSURE)
   CALL BCOND_PA_AIR   
#  endif

#  if defined (GCN)
   CALL BCOND_GCN(1,0)
#  else
   CALL BCOND_GCY(1,0)
#  endif

   IF(NESTING_ON )THEN
      CALL SET_VAR(IntTime,EL=ELF)
   END IF

# endif
! end defined semi-implicit

!==============================================================================!
!==============================================================================!
!                     BEGIN THREE D ADJUSTMENTS
# if !defined (TWO_D_MODEL)
!==============================================================================!
!==============================================================================!

  !==============================================================================!
  !    ADJUST INTERNAL VELOCITY FIELD TO CORRESPOND TO EXTERNAL                  !
  !==============================================================================!
# if defined (NH)
   EGF = ET
# endif

# if !defined (SEMI_IMPLICIT) && !defined (NH)
#   if !defined (ONE_D_MODEL)
    CALL ADJUST2D3D(2)
#   endif
# endif

  ! New Open Boundary Condition ----9
# if defined (MEAN_FLOW)
  !     CALL BCOND_NG_3D       ! change BWI
  CALL BCOND_BKI_3D(2)      
  CALL FLUX_OBC3D
# endif
  
  !==============================================================================!
  !     CALCULATE INTERNAL VELOCITY FLUXES                                       |
  !==============================================================================!
  !                                                    !
!!
# if defined (ICE)
       !  ggao 0104/2008 for ice ocean couple
       DO I=1,NT
         IF(ISICEC(I)==1) THEN

         !! (magnitude of relative ocean current)*rhow*drag*aice
!         DUVI = DRAGW*SQRT((U(I,1)-UICE2(I))**2+(V(I,1)-VICE2(I))**2)  ! m/s
         DUVI = DRAGW*min(SQRT((U(I,1)-UICE2(I))**2+(V(I,1)-VICE2(I))**2),1.0_SP)  ! m/s
         !! ice/ocean stress
!         ice_ocnDX = DUVI*((U(I,1)*COSW-V(I,1)*SINW)-(uice2(i)*cosw-vice2(i)*sinw))
!         ice_ocnDY = DUVI*((V(I,1)*COSW+U(I,1)*SINW)-(vice2(i)*cosw+uice2(i)*sinw))


         ice_ocnDX = DUVI*max(min(((U(I,1)*COSW-V(I,1)*SINW)-(uice2(i)*cosw-vice2(i)*sinw)),1.0_SP),-1.0_SP)
         ice_ocnDY = DUVI*max(min(((V(I,1)*COSW+U(I,1)*SINW)-(vice2(i)*cosw+uice2(i)*sinw)),1.0_SP),-1.0_SP)

         wusurf(I) =wusurf(I) *(1.0_SP-AIU(I)) +ice_ocnDX*aiu(I)*1.0E-3
         wvsurf(I) =wvsurf(I) *(1.0_SP-AIU(I)) +ice_ocnDY*aiu(I)*1.0E-3

         ENDIF
       END DO
#   if defined (MULTIPROCESSOR)
         IF(PAR)CALL AEXCHANGE(EC,MYID,NPROCS,WUSURF,WVSURF)
#   endif
# endif

# endif
!end !defined (TWO_D_MODEL)

!==============================================================================!
!     CALCULATE INTERNAL VELOCITY FLUXES                                       |
!==============================================================================!
# if !defined (SEMI_IMPLICIT)

#  if !defined (TWO_D_MODEL)

    CALL VERTVL_EDGE     ! Calculate/Update Sigma Vertical Velocity (Omega)   !

#   if defined (WET_DRY)
    IF(WETTING_DRYING_ON) CALL WD_UPDATE(2)
#   endif

    CALL VISCOF_H        ! Calculate horizontal diffusion coefficient scalars !
  
  IF(.NOT. RK_3D_ON)THEN
#   if defined (GCN)
    CALL ADV_UV_EDGE_GCN ! Horizontal Advect/Diff + Vertical Advection        !
#   else
    CALL ADV_UV_EDGE_GCY ! Horizontal Advect/Diff + Vertical Advection        !
#   endif

#    if defined (WAVE_CURRENT_INTERACTION)
     CALL VDIF_UV         ! Implicit Integration of Vertical Diffusion of U/V  !
     !    CALL VDIF_UV_WC     
#    else
     CALL VDIF_UV      ! Implicit Integration of Vertical Diffusion of U/V  !
#    endif

    IF(ADCOR_ON) THEN
      CALL ADCOR
#    if defined (WAVE_CURRENT_INTERACTION)
     CALL VDIF_UV         ! Implicit Integration of Vertical Diffusion of U/V  !
     !    CALL VDIF_UV_WC     
#    else
     CALL VDIF_UV      ! Implicit Integration of Vertical Diffusion of U/V  !
#    endif
!      CALL VDIF_UV   ! Implicit Integration of Vertical Diffusion of U/V  !

    ENDIF

#   if defined (WET_DRY)
    DO I=1,N
      IF(H1(I) <= STATIC_SSH_ADJ ) THEN
        DO K=1,KBM1
           UF(I,K)=UA(I)
           VF(I,K)=VA(I)
        END DO
      END IF
    END DO
#   endif

#   if !defined (NH)  
#   if defined (GCN)
    CALL BCOND_GCN(3,0)    ! Boundary Condition on U/V At River Input           !
#   else
    CALL BCOND_GCY(3,0)    ! Boundary Condition on U/V At River Input           !
#   endif

# if defined (PLBC)
     CALL replace_vel_3D
# endif 

#   endif

#   if defined (NH)
    CALL ADV_W
    CALL VDIF_W

    CALL GEN_POISSON
    CALL UPDATE_QUVW
#   endif

  ELSE
   ALLOCATE(UB(0:NT,KB))      ; UB = 0.0_SP
   ALLOCATE(VB(0:NT,KB))      ; VB = 0.0_SP
   UB = U
   VB = V
   RCOFT = 0.5_SP
   DO ITER = 1,2
#   if defined (GCN)
    CALL ADV_UV_EDGE_GCN_RK(UB,VB) ! Horizontal Advect/Diff + Vertical Advection        !
#   else
    CALL ADV_UV_EDGE_GCY ! Horizontal Advect/Diff + Vertical Advection        !
#   endif

     CALL VDIF_UV         ! Implicit Integration of Vertical Diffusion of U/V  !

    IF(ADCOR_ON) THEN
      CALL ADCOR
      CALL VDIF_UV         ! Implicit Integration of Vertical Diffusion of U/V  !
    ENDIF

#   if defined (WET_DRY)
    DO I=1,N
      IF(H1(I) <= STATIC_SSH_ADJ ) THEN
        DO K=1,KBM1
           UF(I,K)=UA(I)
           VF(I,K)=VA(I)
        END DO
      END IF
    END DO
#   endif

#   if defined (GCN)
    CALL BCOND_GCN(3,0)    ! Boundary Condition on U/V At River Input
#   else
    CALL BCOND_GCY(3,0)    ! Boundary Condition on U/V At River Input
#   endif

    U = UF
    V = VF
    IF(ITER /= 2)THEN
      U = RCOFT*U+(1.0_SP-RCOFT)*UB
      V = RCOFT*V+(1.0_SP-RCOFT)*VB
    END IF
   END DO
   DEALLOCATE(UB,VB)
  END IF

#   if !defined (NH)
    IF(NESTING_ON )THEN
      CALL SET_VAR(intTime,U=UF)
      CALL SET_VAR(intTime,V=VF)
    END IF  
#   endif

# endif

# else

#  if defined (TWO_D_MODEL) 

#   if defined (GCN)
    CALL BCOND_GCN(8,0)
#   else
    CALL BCOND_GCY(8,0)
#   endif

    UAC = UA
    VAC = VA
    CALL UV2D_SBD
    DO STAGE=1, KSTAGE_UV
#     if defined (GCN)
      CALL ADVAVE_EDGE_GCN(ADVUA,ADVVA)
#     else
      CALL ADVAVE_EDGE_GCY(ADVUA,ADVVA)
#     endif     
      IF(STAGE<KSTAGE_UV) THEN
        UA = UAF
        VA = VAF
#       if defined (MULTIPROCESSOR)
        IF(PAR)CALL AEXCHANGE(EC,MYID,NPROCS,UA,VA)
        IF(PAR)CALL AEXCHANGE(EC,MYID,NPROCS,WUBOT,WVBOT)
#       endif        
      ELSE
        UA = UAC
        VA = VAC
        CALL SEMI_IMPLICIT_EL
        IF(ADCOR_ON) THEN
          CALL ADCOR
          CALL SEMI_IMPLICIT_EL
        ENDIF
      ENDIF
    ENDDO

    UA  = UAF
    VA  = VAF
#   if defined (MULTIPROCESSOR)
    IF(PAR)CALL AEXCHANGE(EC,MYID,NPROCS,UA,VA)
#   endif

#  else

    UC = U
    VC = V
    CALL UV3D_SBD
    IF(MSTG == 'fast' .AND. KSTAGE_UV>1) THEN
      STAGE = 0
      UAC = UA
      VAC = VA
      CALL UV2D_SBD
#     if defined (GCN)
      CALL ADVECTION_EDGE_GCN(ADVX,ADVY)
#     else
      CALL ADVECTION_EDGE_GCY(ADVX,ADVY)
#     endif

      ADX2D = 0.0_SP ; ADY2D = 0.0_SP
      DRX2D = 0.0_SP ; DRY2D = 0.0_SP
      DO K=1,KBM1
        DO I=1, N
          ADX2D(I)=ADX2D(I)+ADVX(I,K)
          ADY2D(I)=ADY2D(I)+ADVY(I,K)
          DRX2D(I)=DRX2D(I)+DRHOX(I,K)
          DRY2D(I)=DRY2D(I)+DRHOY(I,K)
        ENDDO
      ENDDO

#     if defined (GCN)
      CALL ADVAVE_EDGE_GCN(ADVUA,ADVVA)           !Compute Ext Mode Adv/Diff
#     else
      CALL ADVAVE_EDGE_GCY(ADVUA,ADVVA)           !Compute Ext Mode Adv/Diff
#     endif     
      ADX2D = ADX2D - ADVUA
      ADY2D = ADY2D - ADVVA
    ENDIF

    DO STAGE=1, KSTAGE_UV-1
      IF(MSTG == 'fast') THEN
#       if defined (GCN)
        CALL ADVAVE_EDGE_GCN(ADVUA,ADVVA)
#       else
        CALL ADVAVE_EDGE_GCY(ADVUA,ADVVA)
#       endif
        UA = UAF
        VA = VAF
#       if defined (MULTIPROCESSOR)
        IF(PAR)CALL AEXCHANGE(EC,MYID,NPROCS,UA,VA)
#       endif
      ELSE
#       if defined (GCN)
        CALL ADV_UV_EDGE_GCN
#       else
        CALL ADV_UV_EDGE_GCY
#       endif
        CALL VDIF_UV
#       if defined (WET_DRY)
        DO I=1, N
          IF(H1(I) <= STATIC_SSH_ADJ ) THEN
            UTMP = SUM(UF(I,1:KBM1)*DZ1(I,1:KBM1))
            VTMP = SUM(VF(I,1:KBM1)*DZ1(I,1:KBM1))
            DO K=1, KBM1
              UF(I,K)=UTMP
              VF(I,K)=VTMP
            ENDDO
          ENDIF
        ENDDO
#       endif
#       if defined(MULTIPROCESSOR)
        IF(PAR) CALL AEXCHANGE(EC,MYID,NPROCS,UF,VF)
#       endif
        U = UF
        V = VF
      ENDIF
    ENDDO

    STAGE = KSTAGE_UV
    IF(MSTG == 'fast' .AND. KSTAGE_UV>1) THEN
      DO I=1,NT
        UTMP = SUM(U(I,1:KBM1)*DZ1(I,1:KBM1))
        VTMP = SUM(V(I,1:KBM1)*DZ1(I,1:KBM1))
        U(I,1:KBM1) = U(I,1:KBM1) + (UA(I) - UTMP)
        V(I,1:KBM1) = V(I,1:KBM1) + (VA(I) - VTMP)
      ENDDO
    ENDIF
#   if defined (GCN)
    CALL ADV_UV_EDGE_GCN 
#   else
    CALL ADV_UV_EDGE_GCY
#   endif
    U = UC
    V = VC 
    IF(MSTG == 'fast' .AND. KSTAGE_UV>1) THEN
      UA = UAC
      VA = VAC
    ENDIF
    CALL SEMI_IMPLICIT_EL
# if defined (THIN_DAM)
    IF(NODE_DAM1_N > 0)CALL ELE_DAM1
    IF(NODE_DAM2_N > 0)CALL ELE_DAM2
    IF(NODE_DAM3_N > 0)CALL ELE_DAM3
    IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,ELF) 
    EL=ELF
    CALL N2E2D(ELF,ELF1)
# endif
    CALL VDIF_UV
    IF(ADCOR_ON) THEN
#     if defined (WET_DRY)
      DO I=1, N
        IF(H1(I) <= STATIC_SSH_ADJ ) THEN
          UTMP = SUM(UF(I,1:KBM1)*DZ1(I,1:KBM1))
          VTMP = SUM(VF(I,1:KBM1)*DZ1(I,1:KBM1))
          DO K=1, KBM1
            UF(I,K)=UTMP
            VF(I,K)=VTMP
          ENDDO
        ENDIF
      ENDDO
#     endif
      CALL ADCOR
      CALL SEMI_IMPLICIT_EL
# if defined (THIN_DAM)
      IF(NODE_DAM1_N > 0)CALL ELE_DAM1
      IF(NODE_DAM2_N > 0)CALL ELE_DAM2
      IF(NODE_DAM3_N > 0)CALL ELE_DAM3
      IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,ELF)
      EL=ELF
      CALL N2E2D(ELF,ELF1)
# endif
      CALL VDIF_UV
    ENDIF

!#  if defined (SPHERICAL) !&& (NORTHPOLE)
#   if defined (SPHERICAL) 
    DO J=1, NP
      I = NP_LST(J)
      IF(CELL_NORTHAREA(I)==1) THEN
        DO K=1, KBM1
          UTMP=UF(I,K)
          VTMP=VF(I,K)
          UF(I,K) = VTMP*COS(XC(I)*PI/180.0_SP)-UTMP*SIN(XC(I)*PI/180.0_SP)
          VF(I,K) = -( UTMP*COS(XC(I)*PI/180.0_SP)+VTMP*SIN(XC(I)*PI/180.0_SP) )
        ENDDO
      ENDIF
    ENDDO
#   endif

#   if !defined (NH)
    IF(NESTING_ON )THEN
      CALL SET_VAR(intTime,U=UF)
      CALL SET_VAR(intTime,V=VF)
    END IF
#   endif

#   if defined (WET_DRY)
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
!!    US = U
!!    VS = V
#   endif

    CALL VERTVL_EDGE
#   if defined (WET_DRY)
    IF(WETTING_DRYING_ON) CALL WD_UPDATE(2)
#   endif

    CALL VISCOF_H        ! Calculate horizontal diffusion coefficient scalars !

#   if defined (NH)
    CALL ADV_W
    CALL VDIF_W

    CALL GEN_POISSON
    CALL UPDATE_QUVW
#   endif

    CALL UPDATES_SEMI
    DO I=1,IOBCN            ! do linear T/S case, must redo this part, carefule!
      J=I_OBC_N(I)
      UARD_OBCN(I)=-(EL(J)-ET(J))*ART1(J)/DTI-XFLUX_OBCN(I)
    END DO

#  endif
!  if defined (TWO_D_MODEL) 

#   if defined (EQUI_TIDE)
    EL_EQI = ELF_EQI
#   endif
#   if defined (ATMO_TIDE)
    EL_ATMO = ELF_ATMO
#   endif
#   if defined (AIR_PRESSURE)
    EL_AIR = ELF_AIR
#   endif

#   if !defined (WET_DRY)
    CALL DEPTH_CHECK
#   endif

# endif
!if !defined (SEMI_IMPLICIT)

# if !defined (TWO_D_MODEL)
           
# if !defined (NH)
  CALL WREAL           ! Calculate True Vertical Velocity (W)               !
# endif

!  CALL REPORT("before VISCOF_H")

  !==============================================================================!
  !    TURBULENCE MODEL SECTION                                                  |
  !==============================================================================!
#    if defined (MULTIPROCESSOR)
!     IF(PAR)CALL EXCHANGE(EC,NT,1,MYID,NPROCS,WUSURF,WVSURF)
!     IF(PAR)CALL EXCHANGE(EC,NT,1,MYID,NPROCS,WUBOT,WVBOT)
     IF(PAR)CALL AEXCHANGE(EC,MYID,NPROCS,WUSURF,WVSURF)
     IF(PAR)CALL AEXCHANGE(EC,MYID,NPROCS,WUBOT,WVBOT)
#    endif

  
  SELECT CASE(VERTICAL_MIXING_TYPE)
  CASE('closure')
     !=================General Ocean Turbulence Model==========================!
#   if defined (GOTM)
     
     ! There is a bug in the advection of turbulent kinetic energy
     ! when using the GOTM MODULE
     
     !     CALL ADV_Q(TKE,TKEF)     !!Advection of Tubulent Kinetic Energy
     !     CALL ADV_Q(TEPS,TEPSF)   !!Advection of Turbulent Dissipation Rate 
     !#    if defined(MULTIPROCESSOR)
     !     IF(PAR)CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,KB,MYID,NPROCS,TKEF)
     !     IF(PAR)CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,KB,MYID,NPROCS,TEPSF)
     !#    endif
     !     TKE  = TKEF
     !     TEPS = TEPSF
     CALL ADVANCE_GOTM            !!Solve TKE/TEPS eqns for KH/KM/KQ
     
#   else     

#  if defined (VEGETATION)
      if(VEGETATION_MODEL.and.VEG_TURB)then
      !
      !  Add the effect of vegetation on Q2 and Q2L
      !      
        Q2  = Q2  +  VEG % tke_veg
        Q2L = Q2L +  VEG % tke_veg*L
      end if
#  endif

     !===================Original FVCOM MY-2.5/Galperin 1988 Model=============!
   IF(.NOT. RK_3D_ON)THEN  
#    if defined (SEMI_IMPLICIT)     
     Q2C  = Q2
     Q2LC = Q2L 
     DO STAGE = 1, KSTAGE_TE
#    endif

#     if !defined (ONE_D_MODEL)
#    if defined (SEMI_IMPLICIT)
     TE_TMP = Q2C
#    endif
     CALL ADV_Q(Q2,Q2F)       !!Advection of Q2 

#    if defined (SEMI_IMPLICIT)
     TE_TMP = Q2LC
#    endif
     CALL ADV_Q(Q2L,Q2LF) 

#       if defined(MULTIPROCESSOR)
     IF(PAR)CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,KB,MYID,NPROCS,Q2F)
     IF(PAR)CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,KB,MYID,NPROCS,Q2LF)
#       endif


#    if defined (SEMI_IMPLICIT)
     Q2  = Q2C
     Q2L = Q2LC
#    endif
     IF(SCALAR_POSITIVITY_CONTROL) CALL FCT_Q2             !Conservation Correction   !
     IF(SCALAR_POSITIVITY_CONTROL) CALL FCT_Q2L            !Conservation Correction   !
#     endif
!end !defined (ONE_D_MODEL)


# if defined (PLBC)
  CALL replace_q2
  CALL replace_q2l
# endif

     CALL VDIF_Q                  !! Solve Q2,Q2*L eqns for KH/KM/KQ 
#     if defined (MULTIPROCESSOR)
!     IF(PAR)CALL EXCHANGE(NC,MT,KB,MYID,NPROCS,Q2F,Q2LF,L) !Interprocessor Exchange   !
     IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,Q2F,Q2LF,L) !Interprocessor Exchange   !
#     endif
      Q2  = Q2F
      Q2L = Q2LF

#    if defined (SEMI_IMPLICIT)
     ENDDO
#    endif
  ELSE
    ALLOCATE(Q2B(0:MT,KB))      ;  Q2B = 0.0_SP
    ALLOCATE(Q2LB(0:MT,KB))     ;  Q2LB = 0.0_SP
    Q2B  = Q2
    Q2LB = Q2L
    RCOFT = 0.5_SP
    DO ITER = 1,2
     CALL ADV_Q_RK(Q2,Q2B,Q2F)       !!Advection of Q2 
     CALL ADV_Q_RK(Q2L,Q2LB,Q2LF)
#    if defined(MULTIPROCESSOR)
     IF(PAR)CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,KB,MYID,NPROCS,Q2F)
     IF(PAR)CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,KB,MYID,NPROCS,Q2LF)
#    endif

     IF(SCALAR_POSITIVITY_CONTROL) CALL FCT_Q2             !Conservation Correction   !
     IF(SCALAR_POSITIVITY_CONTROL) CALL FCT_Q2L            !Conservation Correction   !

     CALL VDIF_Q                  !! Solve Q2,Q2*L eqns for KH/KM/KQ 
#    if defined (MULTIPROCESSOR)
     IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,Q2F,Q2LF,L) !Interprocessor Exchange !
#    endif
      Q2  = Q2F
      Q2L = Q2LF

      IF(ITER /= 2)THEN
        Q2 = RCOFT*Q2+(1.0_SP-RCOFT)*Q2B
        Q2L = RCOFT*Q2L+(1.0_SP-RCOFT)*Q2LB
      END IF
    END DO
    DEALLOCATE(Q2B,Q2LB) 
  END IF
# endif
! end if defined(GOTM)
  CASE('constant')
     KM = UMOL
     KH = UMOL*VPRNU
  END SELECT

  
#   if defined (MULTIPROCESSOR)
!     IF(PAR)CALL EXCHANGE(NC,MT,KB,MYID,NPROCS,KM,KQ,KH)
  IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,KM,KQ,KH)
#   endif
  CALL N2E3D(KM,KM1)

# if defined (WAVE_CURRENT_INTERACTION) && !defined (VORTEX_FORCE)
!     U=U-U_STOKES_3D
!     V=V-V_STOKES_3D
# endif



!==============================================================================!
!    SEDIMENT MODEL SECTION  
!      Advance Sed Model (Erode/Deposit/Advect/Diffuse)      
!      Change bathymetry (if MORPHO_MODEL=T)  
!        Note:  morph array returned from sed:  Deposition:  morph > 0
!        0 < morpho_factor < 1  
!      morpho_model and morpho_factor are stored and set in mod_sed.F              
!==============================================================================!
# if defined (SEDIMENT)
# if defined (DATA_ASSIM)
  IF(ASSIM_FLAG.AND.(SST_ASSIM.OR.SSTGRD_ASSIM.OR.SSHGRD_ASSIM.OR.TSGRD_ASSIM))THEN
# endif

   DPDAYS = days(inttime)
   SPDAYS = DPDAYS
# if defined (WAVE_CURRENT_INTERACTION)
   IF(SEDIMENT_MODEL)CALL ADVANCE_SED(DTI,spdays,taucwmax(1:m))
# else
   IF(SEDIMENT_MODEL)CALL ADVANCE_SED(DTI,spdays,taubm_n(1:m))
# endif
  IF(MORPHO_MODEL)THEN
    if(mod(sed_its,morpho_incr)==0 .and. sed_its > morpho_strt)then
       !H(1:m) = H(1:m) - bottom(1:m,morph)
       !CALL N2E2D(H,H1)
       !bottom(1:m,morph) = 0.0
       CALL UPDATE_DEPTH
    endif

!# if defined (MULTIPROCESSOR)
!    IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,H)
!# endif

  ENDIF

# if defined (DATA_ASSIM)
  ENDIF
# endif

# endif
  
#   if defined (BioGen)
  IF(BIOLOGICAL_MODEL)CALL BIO_3D1D
                           
#   endif

# if defined (VEGETATION)
  if(VEGETATION_MODEL)then
     CALL UPDATE_VEGETATION
  end if
# endif

  !==============================================================================!
  !    UPDATE TEMPERATURE IN NON-BAROTROPIC CASE                                 !
  !==============================================================================!
  IF(TEMPERATURE_ACTIVE)THEN
# if !defined (SEMI_IMPLICIT)
  IF(.NOT. RK_3D_ON)THEN
# endif 
     
#    if defined (SEMI_IMPLICIT)
     TSC  = T1
     DO STAGE = 1, KSTAGE_TS
#    endif
#   if !defined (ONE_D_MODEL)                                   
     CALL ADV_T                                     !Advection                 !

#     if defined(MULTIPROCESSOR)
     IF(PAR)CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,KB,MYID,NPROCS,TF1)
#     endif


#    if defined (SEMI_IMPLICIT)
     T1 = TSC
#    endif
     !#                                                   if !defined (DOUBLE_PRECISION)
     IF(SCALAR_POSITIVITY_CONTROL) CALL FCT_T            !Conservation Correction   !
     !#                                                   endif
#   endif
        !end !defined (ONE_D_MODEL)                                    !                          !

     IF(CASENAME(1:3) == 'gom')THEN
        CALL VDIF_TS_GOM(1,TF1)
     ELSE  
#   if !defined (ONE_D_MODEL)
        CALL VDIF_TS(1,TF1)                            !Vertical Diffusion        !
#   else
        CALL VDIF_TS(1,T1)                            !Vertical Diffusion        !
#   endif
     END IF

    
#   if defined (MULTIPROCESSOR)
     IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,TF1) !Interprocessor Exchange   !
#   endif
     CALL BCOND_TS(1)                               !Boundary Conditions       !

# if !defined (SEMI_IMPLICIT)
  ELSE
    ALLOCATE(TB1(0:MT,KB))        ; TB1 = 0.0_SP
    TB1 = T1
    RCOFT = 0.5_SP
    DO ITER = 1,2
     CALL ADV_T_RK(TB1)                                     !Advection                 !
#    if defined(MULTIPROCESSOR)
     IF(PAR)CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,KB,MYID,NPROCS,TF1)
#    endif

     !#                                                   if !defined !(DOUBLE_PRECISION)
     IF(SCALAR_POSITIVITY_CONTROL) CALL FCT_T            !Conservation Correction   !
     !#                                                   endif

     IF(CASENAME(1:3) == 'gom')THEN
        CALL VDIF_TS_GOM(1,TF1)
     ELSE
        CALL VDIF_TS(1,TF1)                            !Vertical Diffusion
     END IF

#   if defined (MULTIPROCESSOR)
     IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,TF1) !Interprocessor Exchange   !
#   endif
     CALL BCOND_TS(1)                               !Boundary Conditions       !
     
     T1 = TF1
     IF(ITER /= 2)THEN
       T1 = RCOFT*T1+(1.0_SP-RCOFT)*TB1
     END IF
    END DO
    DEALLOCATE(TB1)
  ENDIF
# endif

     IF(NESTING_ON )THEN
       CALL SET_VAR(intTime,T1=TF1)
     END IF

!!$!QXU{
!!$#  if defined (SST_GRID_ASSIM)
!!$      IF(ASSIM_FLAG==0 .AND. .NOT. SST_ASSIM_GRD)CALL TEMP_NUDGING
!!$      IF(ASSIM_FLAG==1 .AND. SST_ASSIM_GRD)CALL TEMP_NUDGING
!!$#  else
!!$      IF(ASSIM_FLAG==0 .AND. .NOT. SST_ASSIM)CALL TEMP_NUDGING
!!$      IF(ASSIM_FLAG==1 .AND. SST_ASSIM)CALL TEMP_NUDGING
!!$#  endif
!!$!QXU}
     

#    if defined (DATA_ASSIM)
!!    CALL TEMP_ASSIMILATION                              !Temperature Assimilation        !
!     IF(TS_NGASSIM .OR. TS_OIASSIM)THEN
!      IF(ASSIM_FLAG==0 .AND. .NOT. SST_ASSIM)CALL TEMP_NUDGING
!      IF(ASSIM_FLAG==1 .AND. SST_ASSIM)CALL TEMP_NUDGING
!     END IF
#    endif

#    if defined (DATA_ASSIM_OI)
!     IF(TS_OIASSIM.AND.MOD(IINT,TS_NSTEP_OI).EQ.0)THEN
!      IF(ASSIM_FLAG==0 .AND. .NOT. SST_OIASSIM)CALL TEMP_OI
!      IF(ASSIM_FLAG==1 .AND. SST_OIASSIM)CALL TEMP_OI
!     END IF
#    endif

#    if defined (DATA_ASSIM)
!---> Siqi, SST_improved
     IF(SST_ASSIM .OR. SSTGRD_ASSIM)THEN
!    IF(ASSIM_FLAG==0) CALL SST_SAVE
    IF(ASSIM_FLAG==1) CALL SSTGRD_NUDGE
  END IF
!<--- Siqi

     IF(TS_NGASSIM .OR. (TS_OIASSIM .AND. MOD(IINT,TS_NSTEP_OI) == 0))THEN
      IF(ASSIM_FLAG==0 .AND. (.NOT.SST_ASSIM .AND. .NOT.SSTGRD_ASSIM .AND. &
         .NOT.SSHGRD_ASSIM .AND. .NOT.TSGRD_ASSIM))CALL TEMP_ASSIMILATION
      IF(ASSIM_FLAG==1 .AND. (SST_ASSIM .OR. SSTGRD_ASSIM .OR. SSHGRD_ASSIM &
         .OR. TSGRD_ASSIM))CALL TEMP_ASSIMILATION
     END IF
#    endif

!     IF(NESTING_ON )THEN
!       CALL SET_VAR(intTime,T1=TF1)
!     END IF

#   if !defined (ONE_D_MODEL)
!J. Ge for tracer advection
     IF(BACKWARD_ADVECTION .EQV. .TRUE.)THEN
       IF(BACKWARD_STEP==1)THEN
        T0 = T1
       ELSEIF(BACKWARD_STEP==2)THEN
         T2 = T0
         T0 = T1
       ENDIF
     ENDIF
!J. Ge for tracer advection
     T1 = TF1                                       !Update to new time level  !
#   endif

     CALL N2E3D(T1,T)                               !Shift to Elements         !

#    if defined (SEMI_IMPLICIT)
     ENDDO
#    endif

  END IF                                         !                          !

  !==============================================================================!
  !    UPDATE SALINITY IN NON-BAROTROPIC CASE                                    !
  !==============================================================================!

  IF(SALINITY_ACTIVE)THEN
# if !defined (SEMI_IMPLICIT)
  IF(.NOT. RK_3D_ON)THEN
# endif                   

#    if defined (SEMI_IMPLICIT)
     TSC  = S1
     DO STAGE = 1, KSTAGE_TS
#    endif

#   if !defined (ONE_D_MODEL)
     CALL ADV_S                                     !Advection                 !

#     if defined(MULTIPROCESSOR)
     IF(PAR)CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,KB,MYID,NPROCS,SF1)
#     endif
     

#    if defined (SEMI_IMPLICIT)
     S1=TSC
#    endif
     !#                                                   if !defined (DOUBLE_PRECISION)
     IF(SCALAR_POSITIVITY_CONTROL) CALL FCT_S       !Conservation Correction   !
     !#                                                   endif
#   endif
     !end !defined (ONE_D_MODEL)                                    !                          !

     IF(CASENAME(1:3) == 'gom')THEN
        CALL VDIF_TS_GOM(2,SF1)
     ELSE  
#   if !defined (ONE_D_MODEL)
        CALL VDIF_TS(2,SF1)                            !Vertical Diffusion        !
#   else
        CALL VDIF_TS(2,S1)                            !Vertical Diffusion        !
#   endif
     END IF

#   if defined (MULTIPROCESSOR)
     IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,SF1) !Interprocessor Exchange   !
#   endif

     CALL BCOND_TS(2)                               !Boundary Conditions       !

# if !defined (SEMI_IMPLICIT)
  ELSE
    ALLOCATE(SB1(0:MT,KB))         ; SB1 = 0.0_SP
    SB1 = S1
    RCOFT = 0.5_SP
    DO ITER = 1,2
     CALL ADV_S_RK(SB1)                                     !Advection                 !
#    if defined(MULTIPROCESSOR)
     IF(PAR)CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,KB,MYID,NPROCS,SF1)
#    endif

     !#                                                   if !defined     !(DOUBLE_PRECISION)
     IF(SCALAR_POSITIVITY_CONTROL) CALL FCT_S       !Conservation Correction   !
     !#                                                   endif

     IF(CASENAME(1:3) == 'gom')THEN
        CALL VDIF_TS_GOM(2,SF1)
     ELSE
        CALL VDIF_TS(2,SF1)                            !Vertical Diffusion
     END IF

#   if defined (MULTIPROCESSOR)
     IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,SF1) !Interprocessor Exchange   !
#   endif

     CALL BCOND_TS(2)                               !Boundary Conditions       !

     S1 = SF1
     IF(ITER /= 2)THEN
       S1 = RCOFT*S1+(1.0_SP-RCOFT)*SB1
     END IF
    END DO
    DEALLOCATE(SB1)
  END IF
# endif

     IF(NESTING_ON )THEN
       CALL SET_VAR(intTime,S1=SF1)
     END IF

#   if defined (DATA_ASSIM)
!!    CALL SALT_NUDGING                              !Nudge Salinity            !
!     IF(TS_ASSIM)THEN
!      IF(ASSIM_FLAG==0 .AND. .NOT. SST_ASSIM)CALL SALT_NUDGING
!      IF(ASSIM_FLAG==1 .AND. SST_ASSIM)CALL SALT_NUDGING
!     END IF
#   endif

#   if defined (DATA_ASSIM_OI)
!     IF(TS_OIASSIM.AND.MOD(IINT,TS_NSTEP_OI).EQ.0)THEN
!       IF(ASSIM_FLAG==0 .AND. .NOT. SST_OIASSIM)CALL SALT_OI
!       IF(ASSIM_FLAG==1 .AND. SST_OIASSIM)CALL SALT_OI
!     END IF
#   endif

#   if defined (DATA_ASSIM)
     IF(TS_NGASSIM .OR. (TS_OIASSIM .AND. MOD(IINT,TS_NSTEP_OI) == 0))THEN
       IF(ASSIM_FLAG==0 .AND. (.NOT.SST_ASSIM .AND. .NOT.SSTGRD_ASSIM .AND. &
          .NOT.SSHGRD_ASSIM .AND. .NOT.TSGRD_ASSIM))CALL SALT_ASSIMILATION
       IF(ASSIM_FLAG==1 .AND. (SST_ASSIM .OR. SSTGRD_ASSIM .OR. SSHGRD_ASSIM &
          .OR. TSGRD_ASSIM))CALL SALT_ASSIMILATION
     END IF
#   endif

!     IF(NESTING_ON )THEN
!       CALL SET_VAR(intTime,S1=SF1)
!     END IF

#   if !defined (ONE_D_MODEL)
!J. Ge for tracer advection
     IF(BACKWARD_ADVECTION .EQV. .TRUE.)THEN
       IF(BACKWARD_STEP==1)THEN
        S0 = S1
       ELSEIF(BACKWARD_STEP==2)THEN
         S2 = S0
         S0 = S1
       ENDIF
     ENDIF
!J. Ge for tracer advection
     S1 = SF1                                     !Update to new time level  !
#   endif

     CALL N2E3D(S1,S)                               !Shift to Elements         !

#    if defined (SEMI_IMPLICIT)
     ENDDO
#    endif

  END IF                      

  !==============================================================================!
#   if defined (DYE_RELEASE)
  
  !==============================================================================!
  !    UPDATE DYE IN NON-BAROTROPIC CASE                                         !
  !==============================================================================!
  !     IF(DYE_ON.AND.IINT.GE.IINT_SPE_DYE_B) THEN     !                          !                          !
  IF(DYE_ON) THEN     !
   IF(.NOT. RK_3D_ON)THEN                          !                          !
#    if defined (SEMI_IMPLICIT)
     TSC  = DYE
     DO STAGE = 1, KSTAGE_TS
#    endif
     CALL ADV_DYE                                   !Advection                 !
     
#     if defined(MULTIPROCESSOR)
     IF(PAR)CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,KB,MYID,NPROCS,DYEF)
#     endif
     
#    if defined (SEMI_IMPLICIT)
     DYE = TSC
#    endif
     
#     if !defined (DOUBLE_PRECISION)
     !     IF(TS_FCT) CALL AVER_DYE                      !Conservation Correction   !
#     endif
     
     CALL VDIF_DYE(DYEF)                            !Vertical Diffusion        !
     
     !check
     !!     DYE = DYEF                                     !Update to new time level  !
     !!     IF(IINT.GE.IINT_SPE_DYE_B) CALL ARCHIVE
     !check    
#     if defined (MULTIPROCESSOR)
     IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,DYEF) !Interprocessor Exchange   !
#     endif
     CALL BCOND_DYE                                 !Boundary Conditions       !
     DYE = DYEF                                     !Update to new time level  !
     
#    if defined (SEMI_IMPLICIT)
     ENDDO
#    endif

      !!     IF(MSR) WRITE(IPT,*) 'CALL Dye_on--iint=',iint,IINT_SPE_DYE_B
   ELSE
     ALLOCATE(DYEB(0:MT,KB))     ; DYEB=0.0_SP
     DYEB = DYE
     RCOFT = 0.5_SP
     DO ITER = 1,2
     CALL ADV_DYE_RK(DYEB)                                   !Advection                 !
#    if defined(MULTIPROCESSOR)
     IF(PAR)CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,KB,MYID,NPROCS,DYEF)
#    endif

     CALL VDIF_DYE(DYEF)                            !Vertical Diffusion        !

#    if defined (MULTIPROCESSOR)
     IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,DYEF) !Interprocessor Exchange   !
#    endif
     CALL BCOND_DYE                                 !Boundary Conditions       !
     DYE = DYEF                                     !Update to new time level  !

     IF(ITER /= 2)THEN
       DYE = RCOFT*DYE+(1.0_SP-RCOFT)*DYEB
     ENDIF
    END DO
    DEALLOCATE(DYEB)
   END IF
  END IF                                         !                          !
  !==================================================================================!
#   endif  
  !endif defined (DYE)

!==================================================================================!
!    ADJUST TEMPERATURE AND SALINITY AT RIVER MOUTHS
!==================================================================================!
#   if !defined (MPDATA)
  IF( RIVER_TS_SETTING == 'calculated')THEN
     CALL ADJUST_TS
  END IF
#   endif   
  

#   if defined (WATER_QUALITY)
  !==============================================================================!
  !    CALCULATE WATER QUALITY VARIABLES CONCENTRATIONS                          |
  !==============================================================================!
  !                          !
  IF(WATER_QUALITY_MODEL)THEN                                 !Water Quality Active?     !
     CALL ADV_WQM                                   !Advection                 !
     CALL VDIF_WQM(WQM_F)                           !Vertical Diffusion        !
     CALL EXCHANGE_WQM                              !Interprocessor Exchange   !
     CALL BCOND_WQM                                 !Boundary Conditions       !
!!$     WQM(1:M,1:KBM2,1:NB) = WQM_F(1:M,1:KBM2,1:NB)  !Update                    !
     WQM(1:M,1:KBM1,1:NB) = WQM_F(1:M,1:KBM1,1:NB)  !Update                    !
  END IF                                         !                          !
  !==============================================================================!
#   endif

# if defined (DATA_ASSIM)

# if defined (PWP)
  IF(ASSIM_FLAG==1) CALL PWPGO
# endif

  IF(TSGRD_ASSIM)THEN  ! keep the order
    CALL TSGRD_SAVE
    IF(ASSIM_FLAG==1) CALL TSGRD_NUDGE
  END IF

  IF(SST_ASSIM .OR. SSTGRD_ASSIM)THEN
    IF(ASSIM_FLAG==0) CALL SST_SAVE
! lwang added for SST mld assimilation
    IF(SSTGRD_ASSIM .AND. SSTGRD_MLD .AND. IntTime%MuSOD==0)THEN
      IF(MSR) write(ipt,*) 'RUN MLD_RHO CALCULATION', IntTime%MJD,IntTime%MuSOD
!---> Siqi Li, @20210809
!      CALL MLD_RHO_CAL(RHO_AVG,MLD_RHO)
      CALL MLD_RHO_CAL(M, KBM1, RHO_AVG,MLD_RHO)
!<--- Siqi Li, @20210809
    END IF
! lwang
!    IF(ASSIM_FLAG==1) CALL SSTGRD_NUDGE   ! Siqi, SST_improved
  END IF

  IF(SSHGRD_ASSIM) THEN
    IF(ASSIM_FLAG==0) CALL SSH_SAVE
    IF(ASSIM_FLAG==1) CALL SSHGRD_NUDGE
  ENDIF
# endif

# if defined (DATA_ASSIM)
!==============================================================================!
!     ASSIMILATE SEA SURFACE TEMPERATURE DATA                                  |
!==============================================================================!
!  IF(SST_ASSIM .AND. ASSIM_FLAG == 1) CALL SST_NUDGING !NUDGING SST         !
!  IF(SST_ASSIM .AND. ASSIM_FLAG == 0) CALL SST_INT     !STORE HOURLY SST    !
!  IF(SST_ASSIM .AND. ASSIM_FLAG == 1) CALL N2E3D(T1,T) !RECALCULATE T       !
  
!
!----Recalculate Element Based Temperatures Based on OI Temp Field T1------!
!
!  IF(SST_ASSIM .AND. ASSIM_FLAG == 1)CALL N2E3D(T1,T)
  
!==============================================================================!
# endif

# if defined (DATA_ASSIM_OI)
!==============================================================================!
!     ASSIMILATE SEA SURFACE TEMPERATURE DATA (OI)                                  |
!==============================================================================!
!  IF(SST_OIASSIM .AND. ASSIM_FLAG == 1) CALL SST_OI      !OI  SST           !
!  IF(SST_OIASSIM .AND. ASSIM_FLAG == 0) CALL SST_INT     !STORE HOURLY SST  !
!  IF(SST_OIASSIM .AND. ASSIM_FLAG == 1) CALL N2E3D(T1,T) !RECALCULATE T     !
!
!----Recalculate Element Based Temperatures Based on OI Temp Field T1------!
!
!  IF(SST_OIASSIM .AND. ASSIM_FLAG == 1)CALL N2E3D(T1,T)
  
!==============================================================================!
# endif
  
  
  !==============================================================================!
  !     UPDATE THE DENSITY IN NON-BAROTROPIC CASE                                |
  !==============================================================================!
  IF(.NOT.BAROTROPIC)THEN
   SELECT CASE(SEA_WATER_DENSITY_FUNCTION)
   CASE(SW_DENS1)
      CALL DENS1
   CASE(SW_DENS2)
      CALL DENS2
   CASE(SW_DENS3)
      CALL DENS3
   CASE DEFAULT
      CALL FATAL_ERROR("INVALID DENSITY FUNCTION SELECTED:",&
           & "   "//TRIM(SEA_WATER_DENSITY_FUNCTION) )
   END SELECT
  END IF
  !==============================================================================!
  !     MIMIC CONVECTIVE OVERTURNING TO STABILIZE VERTICAL DENSITY PROFILE       |
  !==============================================================================!
  
  IF(CONVECTIVE_OVERTURNING)THEN
     CALL CONV_OVER
     IF(.NOT.BAROTROPIC)THEN
        SELECT CASE(SEA_WATER_DENSITY_FUNCTION)
        CASE(SW_DENS1)
           CALL DENS1
        CASE(SW_DENS2)
           CALL DENS2
        CASE(SW_DENS3)
           CALL DENS3
        CASE DEFAULT
           CALL FATAL_ERROR("INVALID DENSITY FUNCTION SELECTED:",&
                & "   "//TRIM(SEA_WATER_DENSITY_FUNCTION) )
        END SELECT
        
     END IF
  END IF
!==============================================================================!
!==============================================================================!
# endif     
  ! end if !defined (TWO_D_MODEL)
!==============================================================================!
!==============================================================================!
  
# if defined (DATA_ASSIM)
!==============================================================================!
!     DATA ASSIMILATION FOR CURRENT FIELD                                      |
!==============================================================================!

  IF(CUR_NGASSIM .OR. (CUR_OIASSIM .AND. MOD(IINT,CUR_NSTEP_OI) == 0))THEN
    IF(ASSIM_FLAG==0 .AND. (.NOT. SST_ASSIM .AND. .NOT. SSTGRD_ASSIM .AND. .NOT. SSHGRD_ASSIM .AND. .NOT. TSGRD_ASSIM) )CALL CURRENT_ASSIMILATION
    IF(ASSIM_FLAG==1 .AND. (SST_ASSIM .OR. SSTGRD_ASSIM .OR. SSHGRD_ASSIM .OR. TSGRD_ASSIM)) CALL CURRENT_ASSIMILATION
#   if defined (GCN)
    CALL BCOND_GCN(5,0)
#   else
    CALL BCOND_GCY(5,0)
#   endif
  END IF
!==============================================================================!
# endif
  
# if defined (DATA_ASSIM_OI)
!==============================================================================!
!     DATA ASSIMILATION FOR CURRENT FIELD (OI or NUDGING)                      |
!==============================================================================!
!  IF(CUR_OIASSIM)THEN
!    IF(ASSIM_FLAG==0 .AND. .NOT. SST_OIASSIM)CALL CURRENT_OI
!    IF(ASSIM_FLAG==1 .AND. SST_OIASSIM)CALL CURRENT_OI
!#   if defined (GCN)
!    CALL BCOND_GCN(5)
!#   else
!    CALL BCOND_GCY(5)
!#   endif
!  END IF
!  IF(CUR_ASSIM2)THEN
!    IF(ASSIM_FLAG==0 .AND. .NOT. SST_OIASSIM)CALL CURRENT_NUDGING
!    IF(ASSIM_FLAG==1 .AND. SST_OIASSIM)CALL CURRENT_NUDGING
!#   if defined (GCN)
!    CALL BCOND_GCN(5)
!#   else
!    CALL BCOND_GCY(5)
!#   endif
!  END IF
!==============================================================================!
# endif

# if !defined (TWO_D_MODEL)
  !==============================================================================!
  !     UPDATE VELOCITY FIELD (NEEDED TO WAIT FOR SALINITY/TEMP/TURB/TRACER)     |
  !==============================================================================!
  U = UF
  V = VF
  !==============================================================================!
  !    PERFORM DATA EXCHANGE FOR ELEMENT BASED INFORMATION AT PROC BNDRIES       |
  !==============================================================================!
  
#   if defined (MULTIPROCESSOR)
  IF(PAR)THEN
     CALL AEXCHANGE(EC,MYID,NPROCS,U,V)
#     if defined (GOTM)
     CALL AEXCHANGE(NC,MYID,NPROCS,TKE,TEPS)
#     else 
     CALL AEXCHANGE(NC,MYID,NPROCS,Q2,Q2L)
#     endif
     CALL AEXCHANGE(EC,MYID,NPROCS,RHO,T,S)
     CALL AEXCHANGE(NC,MYID,NPROCS,S1,T1,RHO1)
     
#     if defined (DYE_RELEASE)
     CALL AEXCHANGE(NC,MYID,NPROCS,DYE)
#     endif       
     
     CALL AEXCHANGE(EC,MYID,NPROCS,VISCOFM)
     CALL AEXCHANGE(NC,MYID,NPROCS,VISCOFH)
  END IF
#   endif

!---> Siqi Li, @20210809  
  IF (NC_MLD) THEN
    IF (abs(NC_DAT%FTIME%NEXT_IO -IntTime)<0.1_SP*IMDTI) THEN
      CALL MLD_RHO_CAL(M, KBM1, RHO1(1:M,:), MLD_OUT(1:M)) 
#   if defined (MULTIPROCESSOR)
      IF (PAR) CALL AEXCHANGE(NC,MYID,NPROCS,MLD_OUT) 
#   endif
    END IF
  END IF
!<--- 

  !==============================================================================!
  !     PERFORM DATA EXCHANGE FOR WATER QUALITY VARIABLES                        |
  !==============================================================================!
#   if defined (WATER_QUALITY)
  CALL EXCHANGE_WQM
#   endif
# endif      
  ! end if !defined (TWO_D_MODEL)
  


  !
  !----SHIFT SEA SURFACE ELEVATION AND DEPTH TO CURRENT TIME LEVEL---------------!
  !
  ET  = EL  
  DT  = D 
  ET1 = EL1
  DT1 = D1
  
  IF(WETTING_DRYING_ON) CALL WD_UPDATE(3)
  
  ! New Open Boundary Condition ----10
# if defined (MEAN_FLOW)
  ELTDT = ELT
# endif

# if defined (WAVE_ONLY)
   5501 continue
# endif

!!$  CALL DUMP_PROBE_DATA 


  if(dbg_set(dbg_sbr)) write(ipt,*)&
       & "End: internal_step"

END SUBROUTINE INTERNAL_STEP
