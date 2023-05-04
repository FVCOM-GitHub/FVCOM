!
!! @FVCOM NUOPC Cap
!! @author Jianhua Qi (jqi@umassd.edu)
!! @date 12/20/2019 Original documentation
!------------------------------------------------------
!
#define MULTIPROCESSOR
#define WET_DRY
#define SEMI_IMPLICIT_disabled
#define TVD
#define MPDATA

MODULE NUOPC_FVCOM

  !==============================================================================!
  !  INCLUDE MODULES                                                             !
  !==============================================================================!

  USE MOD_DRIVER
  USE MOD_UTILS
  USE CONTROL
  USE MOD_PAR  
  USE MOD_STARTUP
  USE MOD_TIME
  USE MOD_CLOCK
  USE MOD_INPUT
  USE MOD_NCDIO
  USE MOD_NCLL
  USE MOD_SETUP
  USE MOD_SET_TIME
  USE MOD_FORCE
  USE MOD_OBCS
  USE MOD_REPORT
  USE PROBES
  USE MOD_BOUNDSCHK 
  USE MOD_STATION_TIMESERIES  
  USE MOD_WD
  USE EQS_OF_STATE
  USE MOD_NESTING
! Added by researchers at Akvaplan-niva 2018, idealized tests give promising results.
# if defined (TVD)
  USE MOD_TVD
# endif


  !------------------------------------------------------------------------------|
  IMPLICIT NONE
  PUBLIC 

  CONTAINS
  
  SUBROUTINE NUOPC_FVCOM_init(fvcom_comm,name_domain)

  IMPLICIT NONE

  INTEGER :: IERR,iptt
  INTEGER :: fvcom_comm
  CHARACTER(LEN=20) :: name_domain

  type(watch) Timer
  
  type(TIME) ::GET_BEGIN
  integer status

  !==============================================================================!
  ! INITIALIZE ALL CONTROL VARIABLES
  !==============================================================================!
  CALL INITIALIZE_CONTROL_NUOPC("FVCOM")

  ! INTIALIZE MPI CONTROL VARIABLES
  MPI_FVCOM_GROUP = FVCOM_COMM ! FOR NOW MAKE THEM EQUAL
  MPI_COMM_FVCOM  = FVCOM_COMM ! FOR NOW MAKE THEM EQUAL

  CALL INIT_MPI_ENV_ESMF(MYID,NPROCS,SERIAL,PAR,MSR,MSRID)

!  MPI_FVCOM_GROUP = FVCOM_COMM ! FOR NOW MAKE THEM EQUAL
!  MPI_COMM_FVCOM  = FVCOM_COMM ! FOR NOW MAKE THEM EQUAL

  !==============================================================================!
  !   INITIALIZE A STOP WATCH TIMER FOR TESTING SUBROUTINE EFFICENCY             !
  !==============================================================================!
  CALL WATCH_INIT(TIMER)

  if(DBG_SET(dbg_log)) Call WRITE_BANNER(PAR,NPROCS,MYID)

  !==============================================================================!
  ! SET DEFAULT VALUES AND READ NAME LISTS                                            
  !==============================================================================!

!  print*,"NAMELIST_ESMF"
  CALL NAMELIST_ESMF(NAME_DOMAIN)

  !==============================================================================!
  !   SET MODEL CONTROL PARAMTERS BASED ON NAME LIST HERE                        !
  !==============================================================================!
!  print*,"CNTRL_PRMTRS"
  CALL CNTRL_PRMTRS

  !==============================================================================!
  !   SET THE STARTUP TYPE TO BE USED!                                           !
  !==============================================================================!
!  print*,"SET_STARTUP_TYPE"
  CALL SET_STARTUP_TYPE ! see: startup_type.F

  !==============================================================================!
  !   OPEN ALL FILES NEEDED BASED ON THE RUN PARAMETERS                          !
  !==============================================================================!
!  print*,"OPEN_ALL"
  CALL OPEN_ALL

  !==============================================================================!
  !   SET MODEL TIME BASED ON THE NAMELIST TIME STRINGS OR RESTART FILE          !
  !==============================================================================!
!  print*,"SETUP_TIME"
  CALL SETUP_TIME

  !==============================================================================!
  !   LOAD GRID CONNECTIVITY AND OBC LIST FOR METIS DECOMPOSITION                !
  !==============================================================================!
!  print*,"LOAD_GRID"
  CALL LOAD_GRID

  !==============================================================================!
  !   SETUP THE DOMAIN FOR PARALLEL OR SERIAL RUNNING                            !
  !==============================================================================!
!  PRINT*,"SETUP_DOMAIN"
  CALL SETUP_DOMAIN

  !==============================================================================!
  !   ALLOCATE ALL DOMAIN SIZE VARIABLES HERE                                    !
  !==============================================================================!
!  print*,"ALLOCATE_ALL"
  CALL ALLOCATE_ALL

  !==============================================================================!
  !   LOAD/SETUP PHYSICAL QUANTITIES (CORIOLIS, GRAVITY, SPONGE LAYER, XY/LATLON)!
  !==============================================================================!
!  print*,"COORDS_N_CONST"
  CALL COORDS_N_CONST

  !==============================================================================!
  ! CALCULATE GRID METRICS - NEIGHBORS, GRADIENTS, CELL AREA, INTERP COEFF'S     !
  !==============================================================================!
!  print*,"GRID_METRICS"
  CALL GRID_METRICS
  
   !SETUP TVD ADVECTION
# if defined (TVD)
   CALL SETUP_TVD
# endif

  !==============================================================================!
  !  SETUP THE MODEL FORCING                                                     !
  !==============================================================================!
!  print*,"SETUP_FORCING"
  CALL SETUP_FORCING

# if !defined (SEMI_IMPLICIT)
  IntStep_Seconds = ExtStep_seconds * Isplit
# endif
  print*, cpl_int, IntStep_Seconds
  
  !IF(cpl_int /= IntStep_Seconds) &
  !   CALL FATAL_ERROR("The time step INTSTEP_SECONDS for "//trim(name_domain)//  &
  !                    " is not consistent with NUOPC/ESMF system time step T_STEP_NUOPC")

  !==============================================================================!
  !  SET THE INITIAL CONDITIONS FOR THE MODEL                                    !
  !==============================================================================!
  CALL STARTUP

  !==============================================================================!
  !  CALL ARCHIVE TO SETUP THE OUTPUT AND DUMP CONSTANT VALUES                   !
  !==============================================================================!
  
  CALL ARCHIVE
  ! ORDER MATTERS - ARCHIVE_NEST MUST GO AFTER ARCHIVE DURING SETUP!
  CALL ARCHIVE_NEST

  CALL SET_PROBES(PROBES_ON,PROBES_NUMBER,PROBES_FILE)
  IF(OUT_STATION_TIMESERIES_ON)CALL READ_STATION_FILE

  ! Setup Bounds checking (shutdown if variables exceed threshold)
  CALL SETUP_BOUNDSCHK !bounds checking

  IF(OUT_STATION_TIMESERIES_ON)CALL GET_OUTPUT_FILE_INTERVAL(TRIM(OUT_INTERVAL),INTERVAL_TIME_SERIES)
  IF(OUT_STATION_TIMESERIES_ON)CALL OUT_STATION_TIMESERIES
  IF(OUT_STATION_TIMESERIES_ON) TIME_SERIES = STARTTIME + INTERVAL_TIME_SERIES

     !==============================================================================!
     !  PREPARE TO START FVCOM'S MAIN LOOP                                          !
     !==============================================================================!
     if(DBG_SET(dbg_log)) THEN
        write(IPT,*) "!===================================================="
        write(IPT,*) "!===================================================="
        write(IPT,*) "!============== STARTING MAIN LOOP AT:==============="
        if(DBG_SET(dbg_log)) &
             & Call REPORT_TIME(IINT,ISTART,IEND,IntTime)
        write(IPT,*) "!===================================================="
     end if

     CALL REPORT('INITIAL CONDITIONS')

     if(DBG_SET(dbg_log)) THEN
        write(IPT,*) "!===================================================="
        write(IPT,*) "!===================================================="
        write(IPT,*) "!===================================================="
     end if

  RETURN
  
  END SUBROUTINE NUOPC_FVCOM_init

!===============================================================================
!
!===============================================================================

  SUBROUTINE NUOPC_FVCOM_run(nCplFVCOM)      
  
  USE MOD_NORTHPOLE, only : NP,NP_LST,CELL_NORTHAREA
  USE MOD_NESTING, ONLY : NESTING_ON, SET_VAR,NESTING_DATA
  USE ALL_VARS, ONLY : U, V, T1, S1
  IMPLICIT NONE
!  TYPE(GRID), POINTER :: G1, G4
!  integer :: IOB,IOB_NODE,IOB_CELL,split,iplt1 
  integer :: split 
  INTEGER :: nCplFVCOM
  INTEGER :: I,J,K,ii
  REAL(SP) :: UTMP,VTMP
  
  if(dbg_set(dbg_sbr)) write(ipt,*)&
       & "Start: internal_step"

  DO II = 1, nCplFVCOM
!----SET RAMP FACTOR TO EASE SPINUP--------------------------------------------!
  RAMP = 0.0_SP
  IF(IRAMP /= 0) THEN
     RAMP=TANH(real(IINT,sp)/real(IRAMP,sp))
  ELSE
     RAMP = 1.0_SP
  END IF

  IntTime=IntTime + IMDTI

  !----ADJUST CONSISTENCY BETWEEN 3-D VELOCITY AND VERT-AVERAGED VELOCITIES------!
  CALL ADJUST2D3D(1)

  !----SPECIFY THE SOLID BOUNDARY CONDITIONS OF U&V INTERNAL MODES---------------!
  CALL BCOND_GCN(5,0)

!  print*,'EXCHANGE ',myid,nprocs,par,MPI_FVCOM_GROUP
# if defined(MULTIPROCESSOR)
  IF(PAR)CALL AEXCHANGE(EC,MYID,NPROCS,U,V)
# endif
!  print*,'AFTER EXCHANGE ',myid,nprocs,par,MPI_FVCOM_GROUP

  !----SPECIFY THE SURFACE FORCING OF INTERNAL MODES-----------------------------!
  CALL BCOND_GCN(8,0)

  !----SPECIFY THE BOTTOM ROUGHNESS AND CALCULATE THE BOTTOM STRESSES------------!
  CALL BOTTOM_ROUGHNESS

!==============================================================================!
!  CALCULATE DISPERSION (GX/GY) AND BAROCLINIC PRESSURE GRADIENT TERMS         !
!==============================================================================!
  CALL ADVECTION_EDGE_GCN(ADVX,ADVY)          !Calculate 3-D Adv/Diff       !

  IF(RECALCULATE_RHO_MEAN) THEN
     IF(RECALC_RHO_MEAN .LT. IntTime)THEN
        RECALC_RHO_MEAN = RECALC_RHO_MEAN + DELT_RHO_MEAN
        CALL RHO_PMEAN    
     END IF
  END IF

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

  CALL ADVAVE_EDGE_GCN(ADVUA,ADVVA)           !Compute Ext Mode Adv/Diff

  ADX2D = ADX2D - ADVUA                       !Subtract to Form GX
  ADY2D = ADY2D - ADVVA                       !Subtract to Form GY

  !----INITIALIZE ARRAYS USED TO CALCULATE AVERAGE UA/E  OVER EXTERNAL STEPS-----!
  UARD = 0.0_SP
  VARD = 0.0_SP
  EGF  = 0.0_SP

  IF(IOBCN > 0) THEN
     UARD_OBCN(1:IOBCN)=0.0_SP
  END IF

  !==============================================================================!
  !  LOOP OVER EXTERNAL TIME STEPS                                               !
  !==============================================================================!
  DO IEXT=1,ISPLIT

     IF (DBG_SET(DBG_SBRIO)) WRITE(IPT,*) "/// EXT SETP: ",IEXT

     ExtTime = ExtTime + IMDTE

     CALL EXTERNAL_STEP_ESMF(split)
     
     
  END DO

!==============================================================================!
!==============================================================================!
!                     BEGIN THREE D ADJUSTMENTS
!==============================================================================!
!==============================================================================!

  !==============================================================================!
  !    ADJUST INTERNAL VELOCITY FIELD TO CORRESPOND TO EXTERNAL                  !
  !==============================================================================!
    CALL ADJUST2D3D(2)

!==============================================================================!
!     CALCULATE INTERNAL VELOCITY FLUXES                                       |
!==============================================================================!
    CALL VERTVL_EDGE     ! Calculate/Update Sigma Vertical Velocity (Omega)   !

#   if defined (WET_DRY)
    IF(WETTING_DRYING_ON) CALL WD_UPDATE(2)
#   endif

    CALL VISCOF_H        ! Calculate horizontal diffusion coefficient scalars !
  
    CALL ADV_UV_EDGE_GCN ! Horizontal Advect/Diff + Vertical Advection        !

     CALL VDIF_UV      ! Implicit Integration of Vertical Diffusion of U/V  !

    IF(ADCOR_ON) THEN
      CALL ADCOR
      CALL VDIF_UV      ! Implicit Integration of Vertical Diffusion of U/V  !
    ENDIF

!DEC2021
    IF(NESTING_ON )THEN
      CALL SET_VAR(intTime,U=UF)
      CALL SET_VAR(intTime,V=VF)
    END IF
!DEC2021

#   if defined (WET_DRY)
    IF(WETTING_DRYING_ON)THEN
    DO I=1,N
      IF(H1(I) <= STATIC_SSH_ADJ ) THEN
        DO K=1,KBM1
           UF(I,K)=UA(I)
           VF(I,K)=VA(I)
        END DO
      END IF
    END DO
    END IF
#   endif

    CALL BCOND_GCN(3,0)    ! Boundary Condition on U/V At River Input           !

  CALL WREAL           ! Calculate True Vertical Velocity (W)               !

  !==============================================================================!
  !    TURBULENCE MODEL SECTION                                                  |
  !==============================================================================!
# if defined (MULTIPROCESSOR)
  IF(PAR)CALL AEXCHANGE(EC,MYID,NPROCS,WUSURF,WVSURF)
  IF(PAR)CALL AEXCHANGE(EC,MYID,NPROCS,WUBOT,WVBOT)
# endif
  
  SELECT CASE(VERTICAL_MIXING_TYPE)
  CASE('closure')
     !===================Original FVCOM MY-2.5/Galperin 1988 Model=============!

     CALL ADV_Q(Q2,Q2F)       !!Advection of Q2 

     CALL ADV_Q(Q2L,Q2LF) 

# if defined (MULTIPROCESSOR)
     IF(PAR)CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,KB,MYID,NPROCS,Q2F)
     IF(PAR)CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,KB,MYID,NPROCS,Q2LF)
# endif

     IF(SCALAR_POSITIVITY_CONTROL) CALL FCT_Q2             !Conservation Correction   !
     IF(SCALAR_POSITIVITY_CONTROL) CALL FCT_Q2L            !Conservation Correction   !

     CALL VDIF_Q                  !! Solve Q2,Q2*L eqns for KH/KM/KQ 
# if defined (MULTIPROCESSOR)
     IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,Q2F,Q2LF,L) !Interprocessor Exchange   !
# endif
      Q2  = Q2F
      Q2L = Q2LF

  CASE('constant')
     KM = UMOL
     KH = UMOL*VPRNU
  END SELECT

# if defined (MULTIPROCESSOR)
  IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,KM,KQ,KH)
# endif
  CALL N2E3D(KM,KM1)

  !==============================================================================!
  !    UPDATE TEMPERATURE IN NON-BAROTROPIC CASE                                 !
  !==============================================================================!
  IF(TEMPERATURE_ACTIVE)THEN 
     
     CALL ADV_T                                     !Advection                 !

# if defined(MULTIPROCESSOR)
     IF(PAR)CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,KB,MYID,NPROCS,TF1)
# endif

     !#                                                   if !defined (DOUBLE_PRECISION)
     IF(SCALAR_POSITIVITY_CONTROL) CALL FCT_T            !Conservation Correction   !
     !#                                                   endif

        CALL VDIF_TS(1,TF1)                            !Vertical Diffusion        !

# if defined (MULTIPROCESSOR)
     IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,TF1) !Interprocessor Exchange   !
# endif

     CALL BCOND_TS(1)                               !Boundary Conditions       !

     IF(NESTING_ON )THEN
       CALL SET_VAR(intTime,T1=TF1)
     END IF

!qxu{
!    WHERE (TF1 < -2.0) TF1=-2.0_SP
!qxu}    

!     TT1 = T1
     T1 = TF1                                       !Update to new time level  !

     CALL N2E3D(T1,T)                               !Shift to Elements         !

  END IF                                         !                          !

  !==============================================================================!
  !    UPDATE SALINITY IN NON-BAROTROPIC CASE                                    !
  !==============================================================================!

  IF(SALINITY_ACTIVE)THEN                            !                          !   

     CALL ADV_S                                     !Advection                 !

# if defined(MULTIPROCESSOR)
     IF(PAR)CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,KB,MYID,NPROCS,SF1)
# endif

     !#                                                   if !defined (DOUBLE_PRECISION)
     IF(SCALAR_POSITIVITY_CONTROL) CALL FCT_S       !Conservation Correction   !
     !#                                                   endif

        CALL VDIF_TS(2,SF1)                            !Vertical Diffusion        !

# if defined(MULTIPROCESSOR)
     IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,SF1) !Interprocessor Exchange   !
# endif

     CALL BCOND_TS(2)                               !Boundary Conditions       !

     IF(NESTING_ON )THEN
       CALL SET_VAR(intTime,S1=SF1)
     END IF

!qxu{
!    WHERE (SF1 < 0.0) SF1=0.0_SP
!qxu}

!     ST1 = S1
     S1 = SF1                                       !Update to new time level  !

     CALL N2E3D(S1,S)                               !Shift to Elements         !

  END IF                      

!==================================================================================!
!    ADJUST TEMPERATURE AND SALINITY AT RIVER MOUTHS
!==================================================================================!
# if !defined (MPDATA)
  IF( RIVER_TS_SETTING == 'calculated')THEN
     CALL ADJUST_TS
  END IF
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
  
  !==============================================================================!
  !     UPDATE VELOCITY FIELD (NEEDED TO WAIT FOR SALINITY/TEMP/TURB/TRACER)     |
  !==============================================================================!

  U = UF
  V = VF
  !==============================================================================!
  !    PERFORM DATA EXCHANGE FOR ELEMENT BASED INFORMATION AT PROC BNDRIES       |
  !==============================================================================!
  
# if defined (MULTIPROCESSOR)
  IF(PAR)THEN
     CALL AEXCHANGE(EC,MYID,NPROCS,U,V)
     CALL AEXCHANGE(NC,MYID,NPROCS,Q2,Q2L)
     CALL AEXCHANGE(EC,MYID,NPROCS,RHO,T,S)
     CALL AEXCHANGE(NC,MYID,NPROCS,S1,T1,RHO1)
     
     CALL AEXCHANGE(EC,MYID,NPROCS,VISCOFM)
     CALL AEXCHANGE(NC,MYID,NPROCS,VISCOFH)
  END IF
# endif

  !
  !----SHIFT SEA SURFACE ELEVATION AND DEPTH TO CURRENT TIME LEVEL---------------!
  !
  ET  = EL  
  DT  = D 
  ET1 = EL1
  DT1 = D1
  
  IF(WETTING_DRYING_ON) CALL WD_UPDATE(3)
  
  if(dbg_set(dbg_sbr)) write(ipt,*)&
       & "End: internal_step"

        !==============================================================================!
        !    OUTPUT SCREEN REPORT/TIME SERIES DATA/OUTPUT FILES                        |
        !==============================================================================!
        if(DBG_SET(dbg_log)) &
             & Call REPORT_TIME(IINT,ISTART,IEND,IntTime)
        
        IF(REPORT_NOW(IINT,IREPORT)) CALL REPORT('FLOW FIELD STATS')

        !==============================================================================!
        !  CALL ARCHIVE TO WRITE THE OUTPUT (SELECTED BASED ON INTTIME)                !
        !==============================================================================!

        CALL ARCHIVE

        CALL DUMP_PROBE_DATA
        IF(OUT_STATION_TIMESERIES_ON)CALL OUT_STATION_TIMESERIES
        !==============================================================================!
        !  CALL SHUTDOWN CHECK TO LOOK FOR BAD VALUES                                  !
        !==============================================================================!
        CALL SHUTDOWN_CHECK(D1)

        !==============================================================================!
        !  CALL BOUNDS CHECK TO SEE IF VARIABLES EXCEED USER-DEFINED THRESHOLDS 
        !==============================================================================!
        CALL BOUNDSCHK  !bounds checking

        !==============================================================================!
        !    NESTING OUTPUT                                                            !
        !==============================================================================!
        IF(NCNEST_ON)      CALL ARCHIVE_NEST


        IINT = IINT + 1
	
	if(msr) print*,"IINT = ",IINT
  
     END DO

  RETURN
  
  END SUBROUTINE NUOPC_FVCOM_run

!===============================================================================
!
!===============================================================================

  SUBROUTINE NUOPC_FVCOM_final
  
  IMPLICIT NONE

  if(DBG_SET(dbg_log)) write(IPT,*)"TADA!"
  CALL PSHUTDOWN

  RETURN
  END SUBROUTINE NUOPC_FVCOM_final

  SUBROUTINE INIT_MPI_ENV_ESMF(MYID,NPROCS,SERIAL,PAR,MSR,MSRID) 
!===================================================================================|
!  INITIALIZE MPI ENVIRONMENT                                                       |
!===================================================================================|
     USE CONTROL, ONLY : MPI_COMM_FVCOM
     IMPLICIT NONE
     INTEGER, INTENT(OUT) :: MYID,NPROCS,MSRID
     LOGICAL, INTENT(OUT) :: SERIAL,PAR,MSR
     INTEGER IERR

!     INTEGER :: MPI_COMM_FVCOM
     
     if(DBG_SET(dbg_sbr)) &
          & write(IPT,*) "STARTING INIT_MPI_ENV_ESMF"
     
     IERR=0
     
     CALL MPI_COMM_RANK(MPI_COMM_FVCOM,MYID,IERR)
     IF(IERR/=0) WRITE(*,*) "BAD MPI_COMM_RANK"
     CALL MPI_COMM_SIZE(MPI_COMM_FVCOM,NPROCS,IERR)
     IF(IERR/=0) WRITE(*,*) "BAD MPI_COMM_SIZE"

     print*,"NPROCS = ",NPROCS,MYID,MPI_COMM_FVCOM
          
     MYID = MYID + 1
     MSRID = 1
     IF(NPROCS > 1) SERIAL=.FALSE.
     IF(NPROCS > 1) PAR   =.TRUE.
     IF(MYID /=  1) MSR   =.FALSE.
     
     ! INITIALIZE THE LIST OF MAPS
     nullify(halo_maps%next)
     nullify(halo_maps%map)
     
     nullify(internal_maps%next)
     nullify(internal_maps%map)
     
     ! USE MPI TYPES TO EXCHANGE FVCOM TIME TYPE
     CALL CREATE_MPI_TIME

     if(DBG_SET(dbg_sbr)) &
          & write(IPT,*) "END INIT_MPI_ENV_ESMF"
     
     RETURN
   END SUBROUTINE INIT_MPI_ENV_ESMF
!==============================================================================|
  SUBROUTINE EXTERNAL_STEP_ESMF(split)

  USE ALL_VARS
  USE MOD_NESTING, ONLY : NESTING_ON,SET_VAR,NESTING_DATA

  IMPLICIT NONE
  REAL(SP) :: TMP
  INTEGER :: K, I, J, JN, J1,i1,i2
  TYPE(GRID), POINTER :: G1, G4
  integer :: IOB,IOB_NODE,IOB_CELL,split,iplt1,iplt2 

!------------------------------------------------------------------------------|

  if(dbg_set(dbg_sbr)) write(ipt,*) "Start: external_step_esmf"
 
!----SET RAMP FACTOR TO EASE SPINUP--------------------------------------------!
  IF(IRAMP /= 0) THEN
     TMP = real(IINT-1,sp)+real(IEXT,sp)/real(ISPLIT,sp)
     RAMP=TANH(TMP/real(IRAMP,sp))
  ELSE
     RAMP = 1.0_SP
  END IF

!
!------SURFACE BOUNDARY CONDITIONS FOR EXTERNAL MODEL--------------------------!
!
  CALL BCOND_GCN(9,0)

!
!------SAVE VALUES FROM CURRENT TIME STEP--------------------------------------!
!
  ELRK1 = EL1
  ELRK  = EL
  UARK  = UA
  VARK  = VA

!
!------BEGIN MAIN LOOP OVER EXTERNAL MODEL 4 STAGE RUNGE-KUTTA INTEGRATION-----!
!

  DO K=1,4
     
     RKTIME = ExtTime + IMDTE * (ALPHA_RK(K) - 1.0_DP)
!     print*, 'EXTTIME = ', EXTTIME, RKTIME

!FREE SURFACE AMPLITUDE UPDATE  --> ELF
     CALL EXTEL_EDGE(K)
     IF(PAR) CALL AEXCHANGE(NC,MYID,NPROCS,ELF)
     
     
     ! VALUES FOR THE OPEN BOUNDARY ARE ONLY UPDATED IN THE LOCAL DOMAIN
     ! THE HALO IS NOT SET HERE
     CALL BCOND_GCN(1,K)
 
     IF(NESTING_ON)THEN
!     print*,'BEFORE NESTING ELF'
        CALL SET_VAR(ExtTime,EL=ELF)
!     print*,'AFTER NESTING ELF'
     END IF
     
     DO I=1,IBCN(1)
        JN = OBC_LST(1,I)
        J=I_OBC_N(JN)
        ELF(J)=ELRK(J)+ALPHA_RK(K)*(ELF(J)-ELRK(J))
     END DO

     
     ! DAVID ADDED THIS EXCHANGE CALL:
     ! IT SEEMS LIKELY THAT THE HALO VALUES OF ELF WILL BE USED
     ! BEFORE THEY ARE SET CORRECTLY OTHERWISE
     IF(PAR) CALL AEXCHANGE(NC,MYID,NPROCS,ELF)

     CALL N2E2D(ELF,ELF1)
          
     IF(WETTING_DRYING_ON)CALL WET_JUDGE

     CALL FLUX_OBN(K)

     !CALCULATE ADVECTIVE, DIFFUSIVE, AND BAROCLINIC MODES --> UAF ,VAF
     CALL ADVAVE_EDGE_GCN(ADVUA,ADVVA)           !Compute Ext Mode Adv/Diff

     CALL EXTUV_EDGE(K)

     CALL BCOND_GCN(2,K)

     IF(NESTING_ON )THEN
       CALL SET_VAR(ExtTime,UA=UAF)
       CALL SET_VAR(ExtTime,VA=VAF)
     END IF

     IF(PAR)CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,1,MYID,NPROCS,ELF)


     !UPDATE WATER SURFACE ELEVATION
     CALL ASSIGN_ELM1_TO_ELM2

     EL  = ELF
     EL1 = ELF1

     !!INTERPOLATE DEPTH FROM NODE-BASED TO ELEMENT-BASED VALUES
     CALL N2E2D(EL,EL1)
     
     !UPDATE DEPTH AND VERTICALLY AVERAGED VELOCITY FIELD
     D   = H + EL
     D1  = H1 + EL1
     UA  = UAF
     VA  = VAF
     DTFA = D

     !EXCHANGE ELEMENT-BASED VALUES ACROSS THE INTERFACE
     IF(PAR)CALL AEXCHANGE(EC,MYID,NPROCS,UA,VA,D1)

     !SAVE VALUES FOR 3D MOMENTUM CORRECTION AND UPDATE
     IF(K == 3)THEN
        UARD = UARD + UA*D1
        VARD = VARD + VA*D1
        EGF  = EGF  + EL/ISPLIT
     END IF
     
     !CALCULATE VALUES USED FOR SALINITY/TEMP BOUNDARY CONDITIONS
     IF(K == 4.AND.IOBCN > 0) THEN
        DO I=1,IOBCN
           J=I_OBC_N(I)
           TMP=-(ELF(J)-ELRK(J))*ART1(J)/DTE-XFLUX_OBCN(I)
           UARD_OBCN(I)=UARD_OBCN(I)+TMP/FLOAT(ISPLIT)
        END DO
     END IF

     
     !UPDATE WET/DRY FACTORS
     IF(WETTING_DRYING_ON)CALL WD_UPDATE(1)

  END DO     !! END RUNGE-KUTTA LOOP
  

  if(dbg_set(dbg_sbr)) write(ipt,*) "End: external_step_esmf"

END SUBROUTINE EXTERNAL_STEP_ESMF

  SUBROUTINE INITIALIZE_CONTROL_NUOPC(NAME)
    USE LIMS
    USE CONTROL
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: NAME
    
!==============================================================================!
!  FVCOM VERSION                                                               !
!==============================================================================!

   FVCOM_VERSION     = 'FVCOM_5.0'
   FVCOM_WEBSITE     = 'http://fvcom.smast.umassd.edu, http://codfish.smast.umassd.edu'
   INSTITUTION       = 'School for Marine Science and Technology'

   ! Set the IO UNIT value to screen output for now:
   IPT = 6


!==============================================================================!
!   SETUP PARALLEL ENVIRONMENT                                                 !
!==============================================================================!
   ! DEFAULT SETTINGS
                 SERIAL = .TRUE. 
                    PAR = .FALSE. 
                    MSR = .TRUE.
         IN_MPI_IO_LOOP = .FALSE.
        USE_MPI_IO_MODE = .FALSE.
                   MYID = 1
                  MSRID = 1
                 NPROCS = 1
               PRG_NAME = NAME
            MX_NBR_ELEM = 0
           NPROCS_TOTAL => NPROCS
               IOPROCID =-1

         ZEROTIME%MUSOD = 0
           ZEROTIME%MJD = 0
!# if defined(MULTIPROCESSOR)
        MPI_FVCOM_GROUP = MPI_COMM_FVCOM ! FOR NOW MAKE THEM EQUAL
        MPI_IO_GROUP    = MPI_COMM_FVCOM ! FOR NOW MAKE THEM EQUAL
!# else
!        MPI_FVCOM_GROUP = -1
!        MPI_IO_GROUP    = -1
!# endif
  END SUBROUTINE INITIALIZE_CONTROL_NUOPC

 
END MODULE NUOPC_FVCOM
