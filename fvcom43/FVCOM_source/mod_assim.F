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

MODULE MOD_ASSIM
# if defined (DATA_ASSIM)

  USE MOD_PREC
  USE CONTROL
  USE ALL_VARS
  USE MOD_NCTOOLS
  USE MOD_UTILS
  USE MOD_WD
# if defined (DYE_RELEASE)
  USE MOD_DYE
# endif 
  USE MOD_STARTUP
  USE MOD_INPUT
#  if defined (ENKF)
      USE ENKFVAL,  ONLY : ENKF_NENS
      USE MOD_ENKF,   ONLY : ENKF_memberCONTR
#  endif
#  if defined (SEDIMENT)
      USE MOD_SED
#  endif

  IMPLICIT NONE
  SAVE

   TYPE(NCFILE), POINTER ::NC_RST_ASSIM
   TYPE(NCFILE), POINTER ::NC_TSAVG_ASSIM
   Integer :: Pyear,Pmonth,Pmdays
   !Pmonth : present month
   !Pmdays : total days in present month

  !
  !--Data Assimilation Parameters for SST Assimilation
  !
  LOGICAL  :: SST_ASSIM                    !!TRUE IF SST ASSIMILATION ACTIVE
  CHARACTER(LEN=80) SST_ASSIM_FILE         !!FILE NAME FOR SST DATA
  REAL(SP) :: SST_RADIUS                   !!SEARCH RADIUS FOR INTERPOLATION POINTS
  REAL(SP) :: SST_WEIGHT_MAX                    
  REAL(SP) :: SST_TIMESCALE
  REAL(SP) :: SST_TIME_WINDOW              !!TIME WINDOW FOR OBSERVATION ASSIMILATION (in hours)
  INTEGER  :: SST_N_PER_INTERVAL           !! ASSUMING DAILY DATA:
  !! AVERAGE OVER

  NAMELIST /NML_SST_ASSIMILATION/      &
       & SST_ASSIM,                    &
       & SST_ASSIM_FILE,               &
       & SST_RADIUS,                   &
       & SST_WEIGHT_MAX,               &
       & SST_TIMESCALE,                &
       & SST_TIME_WINDOW,              &
       & SST_N_PER_INTERVAL

  !
  !--Data Assimilation Parameters for SST GRID Assimilation
  !
  LOGICAL  :: SSTGRD_ASSIM                    !!TRUE IF SST ASSIMILATION ACTIVE
  CHARACTER(LEN=80) SSTGRD_ASSIM_FILE         !!FILE NAME FOR SST DATA
  REAL(SP) :: SSTGRD_WEIGHT_MAX                    
  REAL(SP) :: SSTGRD_TIMESCALE
  REAL(SP) :: SSTGRD_TIME_WINDOW              !!TIME WINDOW FOR OBSERVATION ASSIMILATION (in hours)
  INTEGER  :: SSTGRD_N_PER_INTERVAL           !! ASSUMING DAILY DATA:
  !! AVERAGE OVER

  NAMELIST /NML_SSTGRD_ASSIMILATION/      &
       & SSTGRD_ASSIM,                    &
       & SSTGRD_ASSIM_FILE,               &
       & SSTGRD_WEIGHT_MAX,               &
       & SSTGRD_TIMESCALE,                &
       & SSTGRD_TIME_WINDOW,              &
       & SSTGRD_N_PER_INTERVAL


  ! SST ASSIM IS EITHER GRID OR NOT - CAN NOT DO BOTH AT ONCE!

  !--Data Assimilation Parameters for SSH GRID Assimilation
  !
  LOGICAL  :: SSHGRD_ASSIM                    !!TRUE IF SSH ASSIMILATION ACTIVE
  CHARACTER(LEN=80) SSHGRD_ASSIM_FILE         !!FILE NAME FOR SSH DATA
  REAL(SP) :: SSHGRD_WEIGHT_MAX                    
  REAL(SP) :: SSHGRD_TIMESCALE
  REAL(SP) :: SSHGRD_TIME_WINDOW            !!TIME WINDOW FOR OBSERVATION ASSIMILATION (in hours)
  INTEGER  :: SSHGRD_N_PER_INTERVAL           !! ASSUMING DAILY DATA:
  !! AVERAGE OVER

  NAMELIST /NML_SSHGRD_ASSIMILATION/      &
       & SSHGRD_ASSIM,                    &
       & SSHGRD_ASSIM_FILE,               &
       & SSHGRD_WEIGHT_MAX,               &
       & SSHGRD_TIMESCALE,                &
       & SSHGRD_TIME_WINDOW,              &
       & SSHGRD_N_PER_INTERVAL

  !
  !--Data Assimilation Parameters for Climatological TS GRID Assimilation
  !
  LOGICAL  :: TSGRD_ASSIM                    !!TRUE IF SST ASSIMILATION ACTIVE
  CHARACTER(LEN=80) TSGRD_ASSIM_FILE         !!FILE NAME FOR TS DATA
  REAL(SP) :: TSGRD_WEIGHT_MAX                    
  REAL(SP) :: TSGRD_TIMESCALE
  REAL(SP) :: TSGRD_TIME_WINDOW              !!TIME WINDOW FOR OBSERVATION ASSIMILATION (in hours)
  INTEGER  :: TSGRD_N_PER_INTERVAL           !! ASSUMING DAILY DATA:
  !! AVERAGE OVER
  CHARACTER(LEN=80) TSGRD_ASSIM_SAVE_FILE    !!FILE NAME FOR TS ASSIMILATION RESULT

  NAMELIST /NML_TSGRD_ASSIMILATION/      &
       & TSGRD_ASSIM,                    &
       & TSGRD_ASSIM_FILE,               &
       & TSGRD_WEIGHT_MAX,               &
       & TSGRD_TIMESCALE,                &
       & TSGRD_TIME_WINDOW,              &
       & TSGRD_N_PER_INTERVAL

  REAL(SP),ALLOCATABLE :: AW0G(:,:),AWXG(:,:),AWYG(:,:)

!--Data Assimilation Parameters for Current Nudging  Assimilation
!
  LOGICAL  :: CUR_NGASSIM            !!TRUE IF CURRENT ASSIMILATION ACTIVE(nudging)
  CHARACTER(LEN=80) CUR_NGASSIM_FILE   !!FILE NAME FOR CURRENT DATA
  REAL(SP) :: CUR_NG_RADIUS          !!SEARCH RADIUS FOR INTERPOLATION POINTS 
  REAL(SP) :: CUR_GAMA
  REAL(SP) :: CUR_GALPHA
  REAL(SP) :: CUR_NG_ASTIME_WINDOW       !(in hours)

  INTEGER  :: CUR_MAX_LAYER          !!MAXIMUM NUMBER OF VERTICAL DATA FROM ANY OBS POINT

  NAMELIST /NML_CUR_NGASSIMILATION/     &
       & CUR_NGASSIM,                   &
       & CUR_NGASSIM_FILE,              &
       & CUR_NG_RADIUS,                 &
       & CUR_GAMA,                      &
       & CUR_GALPHA,                    &
       & CUR_NG_ASTIME_WINDOW

!
!--Data Assimilation Parameters for Current OI Assimilation
!
  LOGICAL  :: CUR_OIASSIM            !!TRUE IF CURRENT OI ASSIMILATION ACTIVE
  CHARACTER(LEN=80) CUR_OIASSIM_FILE !!FILE NAME FOR CURRENT DATA
  REAL(SP) :: CUR_OI_RADIUS             !!SEARCH RADIUS FOR INTERPOLATION POINTS 
  REAL(SP) :: CUR_OIGALPHA
  REAL(SP) :: CUR_OI_ASTIME_WINDOW    !(in hours)
  INTEGER  :: CUR_N_INFLU            !!NUMBER OF INFLUENTIAL OBSERVATIONS
  REAL(SP),ALLOCATABLE,DIMENSION(:,:)  :: CUR_PARAM
  INTEGER(ITIME)  :: CUR_NSTEP_OI

  NAMELIST /NML_CUR_OIASSIMILATION/    &
       & CUR_OIASSIM,                  &
       & CUR_OIASSIM_FILE,             &
       & CUR_OI_RADIUS,                &
       & CUR_OIGALPHA,                    &
       & CUR_OI_ASTIME_WINDOW,         &
!       & CUR_MAX_LAYER,                &
       & CUR_N_INFLU,                  &
       & CUR_NSTEP_OI

  !
  !--Data Assimilation Parameters for Temp/Salinity Data Nudging Assimilation
  !
  LOGICAL  :: TS_NGASSIM                  !!TRUE IF TEMP/SAL ASSIMILATION ACTIVE
  CHARACTER(LEN=80) TS_NGASSIM_FILE         !!FILE NAME FOR TEMP/SAL DATA
  REAL(SP) :: TS_NG_RADIUS                !!SEARCH RADIUS FOR INTERPOLATION POINTS
  REAL(SP) :: TS_GAMA
  REAL(SP) :: TS_GALPHA
!  REAL(SP) :: TS_WEIGHT_MAX  
!  REAL(SP) :: TS_TIMESCALE  
  REAL(SP) :: TS_NG_ASTIME_WINDOW   !(in hours)

  NAMELIST /NML_TS_NGASSIMILATION/     &
       & TS_NGASSIM,                   &
       & TS_NGASSIM_FILE,              &
       & TS_NG_RADIUS,                 &
!       & TS_WEIGHT_MAX,                &
       & TS_GAMA,                      &
       & TS_GALPHA,                    &
!       & TS_TIMESCALE,                 &
       & TS_NG_ASTIME_WINDOW

!
!--Data Assimilation Parameters for Temp/Salinity Data OI Assimilation
!
  LOGICAL  :: TS_OIASSIM              !!TRUE IF TEMP/SAL OI ASSIMILATION ACTIVE
  CHARACTER(LEN=80):: TS_OIASSIM_FILE !!FILE NAME FOR TEMP/SAL DATA
  REAL(SP) :: TS_OI_RADIUS            !!SEARCH RADIUS FOR INTERPOLATION POINTS
  REAL(SP) :: TS_OIGALPHA  
  REAL(SP) :: TS_OI_ASTIME_WINDOW    !(in hours)
  INTEGER  :: TS_MAX_LAYER         !!MAXIMUM NUMBER OF VERTICAL DATA FROM ANY OBS POINT
  INTEGER  :: TS_N_INFLU           !!NUMBER OF INFLUENTIAL OBSERVATIONS
  REAL(SP),ALLOCATABLE :: TS_PARAM(:,:)

  INTEGER(ITIME)  :: TS_NSTEP_OI
  NAMELIST /NML_TS_OIASSIMILATION/    &
       & TS_OIASSIM,                &
       & TS_OIASSIM_FILE,           &
       & TS_OI_RADIUS,              &
       & TS_OIGALPHA,               &
       & TS_OI_ASTIME_WINDOW,          &
       & TS_MAX_LAYER,              &
       & TS_N_INFLU,                &
       & TS_NSTEP_OI

!qxu{
!
!--Current Assimilation Object Type
!
  TYPE ASSIM_OBJ_CUR
   INTEGER  :: N_TIMES                     !!NUMBER OF DATA TIMES
   INTEGER  :: N_INTPTS                    !!NUMBER OF INTERPOLATION POINTS 
   INTEGER  :: N_T_WEIGHT                  !!DATA TIME FOR CURRENT OBSERVATION WEIGHTING
   INTEGER  :: N_LAYERS                    !!NUMBER OF OBSERVATIONS IN THE VERTICAL
   INTEGER  :: N_CELL                      !!CELL NUMBER OF OBSERVATION
   REAL(SP) :: X,Y                         !!X AND Y COORDINATES OF OBSERVATION
   REAL(SP) :: T_WEIGHT                    !!TIME WEIGHT
   REAL(SP) :: DEPTH                       !!DEPTH AT OBSERVATION STATION (MOORING)
   REAL(SP) :: SITA                        !!LOCAL ISOBATH ANGLE AT OBSERVATION
   REAL(SP),ALLOCATABLE,DIMENSION(:)  :: ODEPTH      !!OBSERVATION DEPTHS
!   REAL(SP),ALLOCATABLE,DIMENSION(:)  :: TIMES       !!DATA TIMES
   TYPE(TIME),ALLOCATABLE,DIMENSION(:)  :: TIMES       !!DATA TIMES
   REAL(SP),ALLOCATABLE,DIMENSION(:,:):: UO,VO       !!OBSERVATION DATA FOR X,Y VELOCITY COMPONENTS
   INTEGER,ALLOCATABLE,DIMENSION(:)   :: INTPTS      !!POINTS USED TO INTERPOLATE TO OBSERVATION LOC 
   INTEGER,ALLOCATABLE,DIMENSION(:,:) :: S_INT       !!SIGMA INTERVALS SURROUNDING CURRENT MEASUREMENT
   REAL(SP),ALLOCATABLE,DIMENSION(:,:):: S_WEIGHT    !!SIGMA WEIGHTING                               
   REAL(SP),ALLOCATABLE,DIMENSION(:)  :: X_WEIGHT    !!SPATIAL WEIGHTING FOR INTERPOLATION POINTS
  END TYPE ASSIM_OBJ_CUR

!
!--Temperature/Salinity Assimilation Object Type
!
   TYPE ASSIM_OBJ_TS  
     INTEGER  :: N_TIMES                                  !!NUMBER OF DATA TIMES
     INTEGER  :: N_INTPTS                                 !!NUMBER OF INTERPOLATION POINTS
     INTEGER  :: N_T_WEIGHT                               !!DATA TIME FOR CURRENT OBSERVATION WEIGHTING
     INTEGER  :: N_LAYERS                                 !!NUMBER OF OBSERVATIONS IN THE VERTICAL
     INTEGER  :: N_CELL                                   !!CELL NUMBER OF OBSERVATION
     REAL(SP) :: X,Y                                      !!X AND Y COORDINATES OF OBSERVATION
     REAL(SP) :: T_WEIGHT                                 !!TIME WEIGHT
     REAL(SP) :: DEPTH                                    !!DEPTH AT OBSERVATION STATION (MOORING)
     REAL(SP) :: SITA                                     !!LOCAL ISOBATH ANGLE AT OBSERVATION
     REAL(SP), ALLOCATABLE, DIMENSION(:)   :: ODEPTH      !!OBSERVATION DEPTHS
!     REAL(SP), ALLOCATABLE, DIMENSION(:)   :: TIMES       !!DATA TIMES
     TYPE(TIME),ALLOCATABLE,DIMENSION(:)  :: TIMES        !!DATA TIMES
     REAL(SP), ALLOCATABLE, DIMENSION(:,:) :: TEMP        !!OBSERVATION DATA FOR TEMPERATURE 
     REAL(SP), ALLOCATABLE, DIMENSION(:,:) :: SAL         !!OBSERVATION DATA FOR SALINITY 
     INTEGER,  ALLOCATABLE, DIMENSION(:)   :: INTPTS      !!POINTS USED TO INTERPOLATE TO OBSERVATION LOC
     INTEGER,  ALLOCATABLE, DIMENSION(:,:) :: S_INT       !!SIGMA INTERVALS SURROUNDING CURRENT MEASUREMENT
     REAL(SP), ALLOCATABLE, DIMENSION(:,:) :: S_WEIGHT    !!SIGMA WEIGHTING
     REAL(SP), ALLOCATABLE, DIMENSION(:)   :: X_WEIGHT    !!SPATIAL WEIGHTING FOR INTERPOLATION POINTS
   END TYPE ASSIM_OBJ_TS    

!
!--Current Data Assimilation Variables
!
   INTEGER                          :: N_ASSIM_CUR !!NUMBER OF CURRENT OBSERVATIONS 
   TYPE(ASSIM_OBJ_CUR), ALLOCATABLE :: CUR_OBS(:)  !!CURRENT ASSIMILATION DATA OBJECTS
   INTEGER, ALLOCATABLE             :: DA_CUR(:)   !!FLAG IF ELEMENT IS USED FOR CURRENT DA INTERP 

!
!--Salinity/Temperature Data Assimilation Variables
!
   INTEGER                                N_ASSIM_TS   !!NUMBER OF TEMPERATURE OBSERVATIONS
   TYPE(ASSIM_OBJ_TS), ALLOCATABLE     :: TS_OBS(:)    !!TEMP ASSIMILATION DATA OBJECTS
   INTEGER, ALLOCATABLE                :: DA_TS(:)     !!FLAG IF NODE IS USED FOR CURRENT DA INTERP

!qxu}
   INTEGER                                N_ASSIM_TS_W   !!NUMBER OF TEMPERATURE O BSERVATIONS
   INTEGER, ALLOCATABLE, DIMENSION(:)   :: ID_STATION


  TYPE ASSIM_OBS_LOC
     INTEGER  :: N_TIMES
     TYPE(TIME), ALLOCATABLE, DIMENSION(:) :: OBSTIME(:)
     INTEGER  :: N_LOCS
     REAL(SP), ALLOCATABLE, DIMENSION(:,:) :: X ! (N_LOCS,N_TIMES)
     REAL(SP), ALLOCATABLE, DIMENSION(:,:) :: Y ! (N_LOCS,N_TIMES)
     REAL(SP), ALLOCATABLE, DIMENSION(:,:) :: Z ! (N_LOCS,N_TIMES)
     REAL(SP), ALLOCATABLE, DIMENSION(:)   :: H
     REAL(SP), ALLOCATABLE, DIMENSION(:)   :: SITA

  END TYPE ASSIM_OBS_LOC

  TYPE ASSIM_OBS_DATA
     REAL(SP), ALLOCATABLE, DIMENSION(:,:) :: DATA
     INTEGER, ALLOCATABLE, DIMENSION(:,:) :: MASK
  END TYPE ASSIM_OBS_DATA


  TYPE ASSIM_OBJ
     CHARACTER(LEN=80) :: VARNAME
     REAL(SP), POINTER, DIMENSION(:,:) :: DATA
     REAL(SP), POINTER, DIMENSION(:) :: X,Y,EL,H
     REAL(SP), POINTER, DIMENSION(:,:) :: SIGMA
     REAL(SP):: RADIUS
     REAL(SP):: WEIGHT_MAX
     REAL(SP):: TIMESCALE  
     REAL(SP):: TIME_WINDOW
     TYPE(NCFILE), POINTER :: FILE
     INTEGER :: N_OBS
     TYPE(ASSIM_OBS_LOC), POINTER :: OBS_LOC
     TYPE(ASSIM_OBS_DATA), POINTER :: OBS_DAT
  END TYPE ASSIM_OBJ

  !
  !--Current Data Assimilation Variables
  !
  TYPE(ASSIM_OBJ), POINTER :: ASSIM_U  !!CURRENT ASSIMILATION DATA OBJECTS
  TYPE(ASSIM_OBJ), POINTER :: ASSIM_V  !!CURRENT ASSIMILATION DATA OBJECTS

  !
  !--Salinity/Temperature Data Assimilation Variables
  !
  TYPE(ASSIM_OBJ), POINTER     :: ASSIM_T    !!TEMP ASSIMILATION DATA OBJECTS
  TYPE(ASSIM_OBJ), POINTER     :: ASSIM_S    !!TEMP ASSIMILATION DATA OBJECTS

  !
  !--SST Data Assimilation Variables
  !
  TYPE(ASSIM_OBJ), POINTER     :: ASSIM_SST   !!SST ASSIMILATION DATA OBJECTS

  !
  !-- GRID SST Data Assimilation Variables
  !

  TYPE(NCFILE), POINTER :: SST_FILE

  INTEGER                             N_TIMES_SST !!NUMBER OF SST OBSERVATIONS TIMES

  INTEGER    :: SST_SAVE_INDEX
  TYPE(TIME) :: SST_SAVE_TIME
  TYPE(TIME)  :: SST_SAVE_INTERVAL
  TYPE(NCVAR), POINTER :: VAR_SST
  TYPE(TIME)  :: ASSIM_RESTART_TIME 

  INTEGER                             ASSIM_FLAG  !!TRUE IF ON ASSIMILATION SWEEP   

  REAL(SP), ALLOCATABLE, TARGET, DIMENSION(:)  :: SST_OBS   !!SST OBS ON GRID
  REAL(SP), ALLOCATABLE, DIMENSION(:,:)  :: SST_SAVED  !!SST ON EACH HOUR DURING SIM 
  REAL(SP), ALLOCATABLE, DIMENSION(:)    :: SST_AVG   !!SIM SST AVERAGED OVER OBSERVATION PERIOD 
!# if defined (ENKF)
  REAL(SP), ALLOCATABLE, DIMENSION(:,:)    :: SST_AVG_ENKF   !!SIM SST AVERAGED OVER OBSERVATION PERIOD FOR EACH ENSEMBLE MEMBER
  REAL(SP), ALLOCATABLE, DIMENSION(:,:,:)  :: SST_SAVED_ENKF  !!SST ON EACH HOUR DURING SIM FOR EACH ENSEMBLE MEMBER
  REAL(SP), ALLOCATABLE, DIMENSION(:,:,:)  :: SSH_SAVED_ENKF  !!SSH ON EACH HOUR DURING SIM FOR EACH MEMEBER
  REAL(SP), ALLOCATABLE, DIMENSION(:,:)    :: SSH_AVG_ENKF   !!SIM SSH AVERAGED OVER OBSERVATION PERIOD FOR EACH MEMBER
!# endif
!qxu{
  TYPE(NCFILE), POINTER :: SSH_FILE

  INTEGER                             N_TIMES_SSH !!NUMBER OF SST OBSERVATIONS TIMES

  INTEGER    :: SSH_SAVE_INDEX
  TYPE(TIME) :: SSH_SAVE_TIME
  TYPE(TIME)  :: SSH_SAVE_INTERVAL
  TYPE(NCVAR), POINTER :: VAR_SSH

  REAL(SP), ALLOCATABLE, TARGET, DIMENSION(:)  :: SSH_OBS   !!SST OBS ON GRID
  REAL(SP), ALLOCATABLE, DIMENSION(:,:)  :: SSH_SAVED  !!SST ON EACH HOUR DURING SIM 
  REAL(SP), ALLOCATABLE, DIMENSION(:)    :: SSH_AVG   !!SIM SST AVERAGED OVER OBSERVATION PERIOD 

  TYPE(NCFILE), POINTER :: TS_FILE

  INTEGER  ::                          N_TIMES_TSC !!NUMBER OF TSC OBSERVATIONS TIMES

  INTEGER    :: TSC_SAVE_INDEX, TSC_SAVE_INDEX2
  TYPE(TIME) :: TSC_SAVE_TIME, TSC_SAVE_TIME2
  TYPE(TIME)  :: TSC_SAVE_INTERVAL
!  TYPE(NCVAR), POINTER :: VAR_TC
!  TYPE(NCVAR), POINTER :: VAR_SC
  !TYPE(TIME)  :: ASSIM_RESTART_TIME 

  !INTEGER                             ASSIM_FLAG  !!TRUE IF ON ASSIMILATION SWEEP   

  REAL(SP), ALLOCATABLE, TARGET, DIMENSION(:,:)  :: TC_OBS   !!Temp OBS ON GRID
  REAL(SP), ALLOCATABLE, TARGET, DIMENSION(:,:)  :: SC_OBS   !!Salinity OBS ON GRID

  TYPE(NCVAR), POINTER :: TSC_T1_N, TSC_T1_P
  TYPE(NCVAR), POINTER :: TSC_S1_N, TSC_S1_P

  REAL(SP), ALLOCATABLE, DIMENSION(:,:,:)  :: TC_SAVED  !!SST ON EACH HOUR DURING SIM 
  REAL(SP), ALLOCATABLE, DIMENSION(:,:,:)  :: SC_SAVED  !!SST ON EACH HOUR DURING SIM 
  REAL(SP), ALLOCATABLE, DIMENSION(:,:)    :: TC_AVG    !!SIM TEMP AVERAGED OVER OBSERVATION PERIOD 
  REAL(SP), ALLOCATABLE, DIMENSION(:,:)    :: SC_AVG    !!SIM SAlinity AVERAGED OVER OBSERVATION PERIOD 

  REAL(SP), ALLOCATABLE, DIMENSION(:,:)    :: TC0_AVG   !!SIM TEMP AVERAGED OVER MONTHLY PERIOD
  REAL(SP), ALLOCATABLE, DIMENSION(:,:)    :: SC0_AVG   !!SIM SALINITY AVERAGED OVER MONTHLY PERIOD

  ! Some IINT VARIABLES TO KEEP TRACK OF PROGRESS
  INTEGER(ITIME) :: INT_COUNT, INT_START, INT_END


  !--save all_var to the *_BUF 

  TYPE(TIME)                              ::IntTime_BUF 
  TYPE(TIME)                              ::ExtTime_BUF 
  INTEGER(ITIME)                          ::IINT_BUF
  REAL(SP), ALLOCATABLE, DIMENSION(:,:)   ::U_BUF
  REAL(SP), ALLOCATABLE, DIMENSION(:,:)   ::V_BUF
  REAL(SP), ALLOCATABLE, DIMENSION(:,:)   ::WTS_BUF
  REAL(SP), ALLOCATABLE, DIMENSION(:,:)   ::S1_BUF
  REAL(SP), ALLOCATABLE, DIMENSION(:,:)   ::T1_BUF
!  REAL(SP), ALLOCATABLE, DIMENSION(:,:)   ::RHO_BUF

  REAL(SP), ALLOCATABLE, DIMENSION(:,:)   ::KM_BUF
  REAL(SP), ALLOCATABLE, DIMENSION(:,:)   ::KH_BUF
  REAL(SP), ALLOCATABLE, DIMENSION(:,:)   ::KQ_BUF

# if defined (EQUI_TIDE)
  REAL(SP), ALLOCATABLE, DIMENSION(:)   :: EL_EQI_BUF
# endif

# if defined (ATMO_TIDE)
  REAL(SP), ALLOCATABLE, DIMENSION(:)   :: EL_ATMO_BUF
# endif

# if defined (GOTM)
  REAL(SP), ALLOCATABLE, DIMENSION(:,:)   :: TKE_BUF
  REAL(SP), ALLOCATABLE, DIMENSION(:,:)   :: TEPS_BUF
# else
  REAL(SP), ALLOCATABLE, DIMENSION(:,:)   ::Q2_BUF
  REAL(SP), ALLOCATABLE, DIMENSION(:,:)   ::Q2L_BUF  
  REAL(SP), ALLOCATABLE, DIMENSION(:,:)   ::L_BUF
# endif

  REAL(SP), ALLOCATABLE, DIMENSION(:)     ::UA_BUF
  REAL(SP), ALLOCATABLE, DIMENSION(:)     ::VA_BUF

  REAL(SP), ALLOCATABLE, DIMENSION(:)     ::EL_BUF
  REAL(SP), ALLOCATABLE, DIMENSION(:)     ::ET_BUF
  REAL(SP), ALLOCATABLE, DIMENSION(:)     ::H_BUF

# if defined (SEDIMENT)
  REAL(SP), ALLOCATABLE, DIMENSION(:,:,:)   ::CONC_BUF
# endif
 

# if defined (DYE_RELEASE)  
  REAL(SP), ALLOCATABLE, DIMENSION(:,:)   ::DYE_BUF
  REAL(SP), ALLOCATABLE, DIMENSION(:,:)   ::DYEF_BUF
!  REAL(SP), ALLOCATABLE, DIMENSION(:,:)   ::DYEMEAN_BUF
# endif 


#   if defined(ICE)
    REAL(SP), ALLOCATABLE, DIMENSION(:,:,:) ::aicen_buf 
    REAL(SP), ALLOCATABLE, DIMENSION(:,:,:) ::vicen_buf
    REAL(SP), ALLOCATABLE, DIMENSION(:,:,:) ::vsnon_buf
    REAL(SP), ALLOCATABLE, DIMENSION(:,:,:) ::tsfcn_buf
    REAL(SP), ALLOCATABLE, DIMENSION(:,:,:) ::esnon_buf
    REAL(SP), ALLOCATABLE, DIMENSION(:,:,:) ::eicen_buf
    REAL(SP), ALLOCATABLE, DIMENSION(:)     ::uice2_buf
    REAL(SP), ALLOCATABLE, DIMENSION(:)     ::vice2_buf
    REAL(SP), ALLOCATABLE, DIMENSION(:,:)   ::fresh_buf
    REAL(SP), ALLOCATABLE, DIMENSION(:,:)   ::fsalt_buf
    REAL(SP), ALLOCATABLE, DIMENSION(:,:)   ::fhnet_buf
    REAL(SP), ALLOCATABLE, DIMENSION(:,:)   ::strength_buf
    INTEGER , ALLOCATABLE, DIMENSION(:)     ::isicec_buf
    REAL(SP), ALLOCATABLE, DIMENSION(:)   ::sig1_buf
    REAL(SP), ALLOCATABLE, DIMENSION(:)   ::sig2_buf
    REAL(SP), ALLOCATABLE, DIMENSION(:)   ::sig12_buf
# endif




CONTAINS !------------------------------------------------------------------!
  ! NAME_LIST_INITIALIZE_ASSIM : INITIALIZE NAMELIST VARIABLES VALUE !
  ! SET_ASSIM_PARAM     :   READ ASSIMILATION PARAMETERS FROM INPUT  !
  ! SETUP_DATA_ASSIMILATION                                          !
  ! OPEN_SST_ASSIM      :                                            !
  ! OPEN_SSTGRD_ASSIM   :                                            !
  ! OPEN_SSHGRD_ASSIM   :                                            !
  ! OPEN_TSGRD_ASSIM    :                                            !
  ! LOAD_SST_ASSIM      :                                            !
  ! SET_SSTGRD_ASSIM    :                                            !
  ! SET_SSHGRD_ASSIM    :                                            !
  ! SET_TSGRD_ASSIM     :                                            !
  ! SET_CUR_ASSIM_DATA  :   READ AND SET CURRENT ASSIMILATION DATA   ! 
  ! SET_TS_ASSIM_DATA   :   READ AND SET TEMP/SAL ASSIMILATION DATA  ! 
  ! CURRENT_ASSIMILATION:   NUDGE CURRENT USING ASSIMILATION         !
  ! TEMP ASSIMILATION   :   NUDGE TEMP USING ASSIMILATION            !
  ! SALT_ASSIMILATION   :   NUDGE SALT USING ASSIMILATION            !
  ! CUR_OPTIMINTERP     :                                            !
  ! TS_OPTIMINTERP      :                                            !
  ! SSTGRD_OBSERVATION_UPDATE                                        !
  ! SSHGRD_OBSERVATION_UPDATE                                        !
  ! SSTGRD_NUDGE        :                                            !
  ! SSHGRD_NUDGE        :                                            !
  ! TSGRD_NUDGE         :                                            !
  ! TSGRD_OBSERVATION_UPDATE                                         !
  ! SST_SAVE            :                                            !
  ! SSH_SAVE            :                                            !
  ! TSGRD_SAVE          :                                            !
  ! RESTORE_STATE       :                                            !
  ! EXCHANGE_ALL        :                                            !
  ! SAVE_STATE          :                                            !
  ! ALLOC_BUFFER        :                                            !
  ! SETUP_RSTFILE_ASSIM :                                            !
  ! DUMP_DATA_ASSIM     :                                            !
  ! SET_MONTHLY_TS_MEAN :                                            ! 
  ! DUMP_TSAVG_ASSIM    :                                            !
  ! -----------------------------------------------------------------!


  !==============================================================================|
  !==============================================================================|

  SUBROUTINE NAME_LIST_INITIALIZE_ASSIM 
    IMPLICIT NONE

    !    NML_SST_ASSIMILIATION
    SST_ASSIM          = .FALSE.
    SST_ASSIM_FILE     = trim(CASENAME)//"_sst.nc"
    SST_RADIUS         = 0.0_SP
    SST_WEIGHT_MAX     = 0.0_SP
    SST_TIMESCALE      = 0.0_SP
    SST_TIME_WINDOW    = 0.0_SP
    SST_N_PER_INTERVAL = 0

    !    NML_SSTGRD_ASSIMILIATION
    SSTGRD_ASSIM          = .FALSE.
    SSTGRD_ASSIM_FILE     = trim(CASENAME)//"_sstgrd.nc"
    SSTGRD_WEIGHT_MAX     = 0.0_SP
    SSTGRD_TIMESCALE      = 0.0_SP
    SSTGRD_TIME_WINDOW    = 0.0_SP
    SSTGRD_N_PER_INTERVAL = 0

    !    NML_SSHGRD_ASSIMILIATION
    SSHGRD_ASSIM          = .FALSE.
    SSHGRD_ASSIM_FILE     = trim(CASENAME)//"_sshgrd.nc"
    SSHGRD_WEIGHT_MAX     = 0.0_SP
    SSHGRD_TIMESCALE      = 0.0_SP
    SSHGRD_TIME_WINDOW    = 0.0_SP
    SSHGRD_N_PER_INTERVAL = 0

    !    NML_TSGRD_ASSIMILIATION
    TSGRD_ASSIM          = .FALSE.
    TSGRD_ASSIM_FILE     = trim(CASENAME)//"_tsgrd.nc"
    TSGRD_WEIGHT_MAX     = 0.0_SP
    TSGRD_TIMESCALE      = 0.0_SP
    TSGRD_TIME_WINDOW    = 0.0_SP
    TSGRD_N_PER_INTERVAL = 0

!    !    NML_CUR_ASSIMILIATION
!    CUR_ASSIM       = .FALSE.
!    CUR_ASSIM_FILE  = trim(CASENAME)//"_cur.nc"
!    CUR_RADIUS      = 0.0_SP
!    CUR_WEIGHT_MAX  = 0.0_SP
!    CUR_TIMESCALE   = 0.0_SP
!    CUR_TIME_WINDOW = 0.0_SP

    !    NML_CUR_NGASSIMILIATION
    CUR_NGASSIM            = .FALSE.
    CUR_NGASSIM_FILE       = trim(CASENAME)//"_cur"
    CUR_GAMA               = 0.0_SP
    CUR_NG_RADIUS          = 0.0_SP
    CUR_GALPHA             = 0.0_SP
    CUR_NG_ASTIME_WINDOW   = 0.0_SP

    !    NML_CUR_OIASSIMILIATION
    CUR_OIASSIM            = .FALSE.
    CUR_OIASSIM_FILE       = trim(CASENAME)//"_cur"
    CUR_OI_RADIUS          = 0.0_SP 
    CUR_OI_ASTIME_WINDOW   = 0.0_SP
    CUR_MAX_LAYER          = 0.0_SP
    CUR_N_INFLU            = 0.0_SP

    !    NML_TS_NGASSIMILIATION
    TS_NGASSIM             = .FALSE.
    TS_NGASSIM_FILE        = trim(CASENAME)//"_ts"
    TS_NG_RADIUS           = 0.0_SP
    TS_GAMA                = 0.0_SP
    TS_GALPHA              = 0.0_SP
    TS_NG_ASTIME_WINDOW    = 0.0_SP

    !    NML_TS_OIASSIMILIATION
    TS_OIASSIM             = .FALSE.
    TS_OIASSIM_FILE        = trim(CASENAME)//"_ts"
    TS_OI_RADIUS           = 0.0_SP
    TS_OIGALPHA            = 0.0_SP
    TS_OI_ASTIME_WINDOW    = 0.0_SP
    TS_MAX_LAYER           = 0.0_SP
    TS_N_INFLU             = 0.0_SP
    TS_NSTEP_OI            = 0

  END SUBROUTINE NAME_LIST_INITIALIZE_ASSIM 

  SUBROUTINE SET_ASSIM_PARAM 

    !------------------------------------------------------------------------------|
    !  READ IN PARAMETERS CONTROLLING ASSIMILATION                                 |
    !------------------------------------------------------------------------------|

    IMPLICIT NONE
    CHARACTER(LEN=120) :: FNAME, pathnfile
    INTEGER :: IOS, charnum
    LOGICAL FEXIST


    IF (DBG_SET(DBG_SBR)) WRITE(IPT,*) "START: SET_ASSIM_PARAM"

    charnum = index (DATA_ASSIMILATION_FILE,".nml")
    if (charnum /= len_trim(DATA_ASSIMILATION_FILE)-3)&
         & CALL WARNING("DATA ASSIMILATION FILE does not end in .nml", &
         & trim(DATA_ASSIMILATION_FILE))
    ! OPEN FILE - try both with appending input dir and without!
    pathnfile = trim(INPUT_DIR)//trim(DATA_ASSIMILATION_FILE)
    INQUIRE(FILE=PATHNFILE,EXIST=FEXIST)
    IF(FEXIST) THEN
       Call FOPEN(ASSIMUNIT,trim(pathnfile),'cfr')
    ELSE
       pathnfile = trim(DATA_ASSIMILATION_FILE)
       Call FOPEN(ASSIMUNIT,trim(pathnfile),'cfr')
    END IF


    CALL NAME_LIST_INITIALIZE_ASSIM

    ! SAVE FILE NAME USED TO FNAME FOR ERROR CHECKING
    FNAME =  pathnfile

    ! Read SST ASSIMILATION Settings
    READ(UNIT=ASSIMUNIT, NML=NML_SST_ASSIMILATION,IOSTAT=ios)
    if(ios .NE. 0 ) then
       if(DBG_SET(dbg_log)) write(UNIT=IPT,NML=NML_SST_ASSIMILATION)
       CALL FATAL_ERROR("Can Not Read NameList NML_SST_ASSIMILATION from file: "//trim(FNAME))
    end if
    SST_TIME_WINDOW = SST_TIME_WINDOW*3600.0   ! convert hours -> seconds
    REWIND(ASSIMUNIT)

    if(DBG_SET(dbg_scl)) &
         & write(IPT,*) "Read_Name_List:"

    if(DBG_SET(dbg_scl)) &
         & write(UNIT=IPT,NML=NML_SST_ASSIMILATION)


    ! Read SST GRID ASSIMILATION Settings
    READ(UNIT=ASSIMUNIT, NML=NML_SSTGRD_ASSIMILATION,IOSTAT=ios)
    if(ios .NE. 0 ) then
       if(DBG_SET(dbg_log)) write(UNIT=IPT,NML=NML_SSTGRD_ASSIMILATION)
       CALL FATAL_ERROR("Can Not Read NameList NML_SSTGRD_ASSIMILATION from file: "//trim(FNAME))
    end if
   SSTGRD_TIME_WINDOW = SSTGRD_TIME_WINDOW*3600.0   ! convert hours -> seconds
    REWIND(ASSIMUNIT)


    if(DBG_SET(dbg_scl)) &
         & write(IPT,*) "Read_Name_List:"

    if(DBG_SET(dbg_scl)) &
         & write(UNIT=IPT,NML=NML_SSTGRD_ASSIMILATION)

    ! Read SSH GRID ASSIMILATION Settings
    READ(UNIT=ASSIMUNIT, NML=NML_SSHGRD_ASSIMILATION,IOSTAT=ios)
    if(ios .NE. 0 ) then
       if(DBG_SET(dbg_log)) write(UNIT=IPT,NML=NML_SSHGRD_ASSIMILATION)
       CALL FATAL_ERROR("Can Not Read NameList NML_SSHGRD_ASSIMILATION from file: "//trim(FNAME))
    end if
   SSHGRD_TIME_WINDOW = SSHGRD_TIME_WINDOW*3600.0   ! convert hours -> seconds
    REWIND(ASSIMUNIT)


    if(DBG_SET(dbg_scl)) &
         & write(IPT,*) "Read_Name_List:"

    if(DBG_SET(dbg_scl)) &
         & write(UNIT=IPT,NML=NML_SSHGRD_ASSIMILATION)

    ! Read TS GRID ASSIMILATION Settings
    READ(UNIT=ASSIMUNIT, NML=NML_TSGRD_ASSIMILATION,IOSTAT=ios)
    if(ios .NE. 0 ) then
       if(DBG_SET(dbg_log)) write(UNIT=IPT,NML=NML_TSGRD_ASSIMILATION)
       CALL FATAL_ERROR("Can Not Read NameList NML_TSGRD_ASSIMILATION from file: "//trim(FNAME))
    end if
    TSGRD_TIME_WINDOW = TSGRD_TIME_WINDOW*3600.0   ! convert hours -> seconds
    REWIND(ASSIMUNIT)
    
    if(DBG_SET(dbg_scl)) &
         & write(IPT,*) "Read_Name_List:"

    if(DBG_SET(dbg_scl)) &
         & write(UNIT=IPT,NML=NML_TSGRD_ASSIMILATION)

    ! Read CUR NUDGING ASSIMILATION Settings
    rewind(ASSIMUNIT)
    READ(UNIT=ASSIMUNIT, NML=NML_CUR_NGASSIMILATION,IOSTAT=ios)
    if(ios .NE. 0 ) then
       if(DBG_SET(dbg_log)) write(UNIT=IPT,NML=NML_CUR_NGASSIMILATION)
       CALL FATAL_ERROR("Can Not Read NameList NML_CUR_NGASSIMILATION from file: "//trim(FNAME))
    end if
    CUR_NG_ASTIME_WINDOW = CUR_NG_ASTIME_WINDOW*3600.0      !convert hours -> seconds
    REWIND(ASSIMUNIT)

    if(DBG_SET(dbg_scl)) &
         & write(IPT,*) "Read_Name_List:"

    if(DBG_SET(dbg_scl)) &
         & write(UNIT=IPT,NML=NML_CUR_NGASSIMILATION)

    ! Read CUR OI ASSIMILATION Settings
    READ(UNIT=ASSIMUNIT, NML=NML_CUR_OIASSIMILATION,IOSTAT=ios)
    if(ios .NE. 0 ) then
       if(DBG_SET(dbg_log)) write(UNIT=IPT,NML=NML_CUR_OIASSIMILATION)
       CALL FATAL_ERROR("Can Not Read NameList NML_CUR_OIASSIMILATION from file: "//trim(FNAME))
    end if
    CUR_OI_ASTIME_WINDOW = CUR_OI_ASTIME_WINDOW*3600.0      !convert hours -> seconds
    REWIND(ASSIMUNIT)

    if(DBG_SET(dbg_scl)) &
         & write(IPT,*) "Read_Name_List:"

    if(DBG_SET(dbg_scl)) &
         & write(UNIT=IPT,NML=NML_CUR_OIASSIMILATION)

    ! Read TS NUDGING ASSIMILATION Settings
    READ(UNIT=ASSIMUNIT, NML=NML_TS_NGASSIMILATION,IOSTAT=ios)
    if(ios .NE. 0 ) then
       if(DBG_SET(dbg_log)) write(UNIT=IPT,NML=NML_TS_NGASSIMILATION)
       CALL FATAL_ERROR("Can Not Read NameList NML_TS_NGASSIMILATION from file: "//trim(FNAME))
    end if
    TS_NG_ASTIME_WINDOW = TS_NG_ASTIME_WINDOW*3600.0      !convert hours -> seconds
    REWIND(ASSIMUNIT)

    if(DBG_SET(dbg_scl)) &
         & write(IPT,*) "Read_Name_List:"

    if(DBG_SET(dbg_scl)) &
         & write(UNIT=IPT,NML=NML_TS_NGASSIMILATION)

    ! Read TS OI ASSIMILATION Settings
    READ(UNIT=ASSIMUNIT, NML=NML_TS_OIASSIMILATION,IOSTAT=ios)
    if(ios /= 0 ) then
       if(DBG_SET(dbg_log)) write(UNIT=IPT,NML=NML_TS_OIASSIMILATION)
       CALL FATAL_ERROR("Can Not Read NameList NML_TS_OIASSIMILATION from file: "//trim(FNAME))
    end if
    TS_OI_ASTIME_WINDOW = TS_OI_ASTIME_WINDOW*3600.0      !convert hours -> seconds
    REWIND(ASSIMUNIT)

    if(DBG_SET(dbg_scl)) &
         & write(IPT,*) "Read_Name_List:"

    if(DBG_SET(dbg_scl)) &
         & write(UNIT=IPT,NML=NML_TS_OIASSIMILATION)

!------------------------------------------------------------------------------|
!   MODIFY VARIABLES TO CORRESPOND TO CORRECT UNITS              
!------------------------------------------------------------------------------|
!   IF(CUR_NGASSIM)CUR_NG_ASTIME_WINDOW = CUR_NG_ASTIME_WINDOW*3600. !!CONVERT HOURS --> SECONDS
!   IF(CUR_OIASSIM)CUR_OI_ASTIME_WINDOW = CUR_OI_ASTIME_WINDOW*3600. !!CONVERT HOURS --> SECONDS

!   IF(TS_NGASSIM)TS_NG_ASTIME_WINDOW = TS_NG_ASTIME_WINDOW*3600. !!CONVERT HOURS --> SECONDS
!   IF(TS_OIASSIM)TS_OI_ASTIME_WINDOW = TS_OI_ASTIME_WINDOW*3600. !!CONVERT HOURS --> SECONDS

    IF (DBG_SET(DBG_LOG)) THEN
       WRITE(IPT,*)''
       WRITE(IPT,*)'!        DATA ASSIMILATION PARAMETERS       '

!       IF(CUR_ASSIM)THEN
!          WRITE(IPT,*)'!  # CUR_ASSIM       :  ACTIVE'
!          WRITE(IPT,*)'!  # CUR_ASSIM_FILE      : '//TRIM(CUR_ASSIM_FILE)
!          WRITE(IPT,*)'!  # CUR_RADIUS          :',CUR_RADIUS
!          WRITE(IPT,*)'!  # CUR_WEIGHT_MAX      :',CUR_WEIGHT_MAX
!          WRITE(IPT,*)'!  # CUR_TIMESCALE       :',CUR_TIMESCALE   
!          WRITE(IPT,*)'!  # CUR_TIME_WINDOW     :',CUR_TIME_WINDOW
!          WRITE(IPT,*)''
!       ELSE
!          WRITE(IPT,*)'!  # CUR_ASSIM           :  NOT ACTIVE'
!          WRITE(IPT,*)''
!       END IF

!qxu{
       IF(CUR_NGASSIM)THEN
          WRITE(IPT,*)'!  # CUR_NGASSIM         :  ACTIVE'
          WRITE(IPT,*)'!  # CUR_NGASSIM_FILE    : '//TRIM(CUR_NGASSIM_FILE)
          WRITE(IPT,*)'!  # CUR_NG_RADIUS       :',CUR_NG_RADIUS
          WRITE(IPT,*)'!  # CUR_GALPHA          :',CUR_GALPHA
          WRITE(IPT,*)'!  # CUR_GAMA            :',CUR_GAMA   
          WRITE(IPT,*)'!  # CUR_NG_ASTIME_WINDOW (seconds):',CUR_NG_ASTIME_WINDOW
       ELSE
          WRITE(IPT,*)'!  # CUR_NGASSIM         :  NOT ACTIVE'
       END IF
!qxu}       

       IF(CUR_OIASSIM)THEN
          WRITE(IPT,*)'!  # CUR_OIASSIM         :  ACTIVE'
          WRITE(IPT,*)'!  # CUR_OIASSIM_FILE    : '//TRIM(CUR_OIASSIM_FILE)
          WRITE(IPT,*)'!  # CUR_OI_RADIUS       :',CUR_OI_RADIUS
          WRITE(IPT,*)'!  # CUR_OI_ASTIME_WINDOW (seconds):',CUR_OI_ASTIME_WINDOW
!          WRITE(IPT,*)'!  # CUR_MAX_LAYER       :',CUR_MAX_LAYER
          WRITE(IPT,*)'!  # CUR_N_INFLU         :',CUR_N_INFLU
       ELSE
          WRITE(IPT,*)'!  # CUR_OIASSIM          :  NOT ACTIVE'
       END IF
       
       IF(SST_ASSIM)THEN
          WRITE(IPT,*)'!  # SST_ASSIM           :  ACTIVE'
          WRITE(IPT,*)'!  # SST_RADIUS          :',SST_RADIUS
          WRITE(IPT,*)'!  # SST_ASSIM_FILE      : '//TRIM(SST_ASSIM_FILE)
          WRITE(IPT,*)'!  # SST_TIMESCALE       :',SST_TIMESCALE
          WRITE(IPT,*)'!  # SST_WEIGHT_MAX      :',SST_WEIGHT_MAX   
          WRITE(IPT,*)'!  # SST_TIME_WINDOW(seconds)     :',SST_TIME_WINDOW
          WRITE(IPT,*)'!  # SST_N_PER_INTERVAL  :',SST_N_PER_INTERVAL
          WRITE(IPT,*)''
       ELSE
          WRITE(IPT,*)'!  # SST_ASSIM           :  NOT ACTIVE'
          WRITE(IPT,*)''
       END IF

       IF(SSTGRD_ASSIM)THEN
          WRITE(IPT,*)'!  # SSTGRD_ASSIM           :  ACTIVE'
          WRITE(IPT,*)'!  # SSTGRD_ASSIM_FILE      : '//TRIM(SSTGRD_ASSIM_FILE)
          WRITE(IPT,*)'!  # SSTGRD_TIMESCALE       :',SSTGRD_TIMESCALE
          WRITE(IPT,*)'!  # SSTGRD_WEIGHT_MAX      :',SSTGRD_WEIGHT_MAX   
          WRITE(IPT,*)'!  # SSTGRD_TIME_WINDOW (seconds)     :',SSTGRD_TIME_WINDOW
          WRITE(IPT,*)'!  # SSTGRD_N_PER_INTERVAL  :',SSTGRD_N_PER_INTERVAL
          WRITE(IPT,*)''
       ELSE
          WRITE(IPT,*)'!  # SSTGRD_ASSIM           :  NOT ACTIVE'
          WRITE(IPT,*)''
       END IF

       IF(TSGRD_ASSIM)THEN
          WRITE(IPT,*)'!  # TSGRD_ASSIM           :  ACTIVE'
          WRITE(IPT,*)'!  # TSGRD_ASSIM_FILE      : '//TRIM(TSGRD_ASSIM_FILE)
          WRITE(IPT,*)'!  # TSGRD_TIMESCALE       :',TSGRD_TIMESCALE
          WRITE(IPT,*)'!  # TSGRD_WEIGHT_MAX      :',TSGRD_WEIGHT_MAX  
          WRITE(IPT,*)'!  # TSGRD_TIME_WINDOW (seconds)     :',TSGRD_TIME_WINDOW
          WRITE(IPT,*)'!  # TSGRD_N_PER_INTERVAL  :',TSGRD_N_PER_INTERVAL
          WRITE(IPT,*)''
       ELSE
          WRITE(IPT,*)'!  # TSGRD_ASSIM           :  NOT ACTIVE'
          WRITE(IPT,*)''
       END IF

       IF(TS_NGASSIM)THEN
          WRITE(IPT,*)'!  # TS_NGASSIM       :  ACTIVE'
          WRITE(IPT,*)'!  # TS_NGASSIM_FILE     : '//TRIM(TS_NGASSIM_FILE)
          WRITE(IPT,*)'!  # TS_NG_RADIUS           :',TS_NG_RADIUS
          WRITE(IPT,*)'!  # TS_GAMA             :',TS_GAMA
          WRITE(IPT,*)'!  # TS_GALPHA           :',TS_GALPHA   
!          WRITE(IPT,*)'!  # TS_WEIGHT_MAX      :',TS_WEIGHT_MAX
!          WRITE(IPT,*)'!  # TS_TIMESCALE       :',TS_TIMESCALE   
          WRITE(IPT,*)'!  # TS_NG_ASTIME_WINDOW(seconds)     :',TS_NG_ASTIME_WINDOW
          WRITE(IPT,*)''
       ELSE
          WRITE(IPT,*)'!  # TS_ASSIM           :  NOT ACTIVE'
          WRITE(IPT,*)''
       END IF

       IF(TS_OIASSIM)THEN
          WRITE(IPT,*)'!  # TS_OIASSIM          :  ACTIVE'
          WRITE(IPT,*)'!  # TS_OIASSIM_FILE     : '//TRIM(TS_OIASSIM_FILE)
          WRITE(IPT,*)'!  # TS_OI_RADIUS        :',TS_OI_RADIUS
          WRITE(IPT,*)'!  # TS_OIGALPHA         :',TS_OIGALPHA
          WRITE(IPT,*)'!  # TS_OI_ASTIME_WINDOW(seconds) :',TS_OI_ASTIME_WINDOW
!         WRITE(IPT,*)'!  # TS_MAX_LAYER        :',TS_MAX_LAYER
          WRITE(IPT,*)'!  # TS_N_INFLU          :',TS_N_INFLU
          WRITE(IPT,*)'!  # TS_NSTEP_OI         :',TS_NSTEP_OI
          WRITE(IPT,*)''
       ELSE
          WRITE(IPT,*)'!  # TS_OIASSIM         :  NOT ACTIVE'
          WRITE(IPT,*)''
       END IF
    END IF
    
    IF(SST_ASSIM .AND. SSTGRD_ASSIM) CALL FATAL_ERROR &
         & ('Using two kinds of sst assimilation does not make any sense!')

    IF(TS_NGASSIM .AND. TSGRD_ASSIM) CALL FATAL_ERROR &
         & ('Using two kinds of ts assimilation does not make any sense!')

    IF(TS_OIASSIM .AND. TSGRD_ASSIM) CALL FATAL_ERROR &
         & ('Using two kinds of ts assimilation does not make any sense!')

    IF(CUR_OIASSIM.AND.CUR_NGASSIM) CALL FATAL_ERROR &
         & ('Using two kinds of current assimilation does not make any sense!')

    IF(TS_OIASSIM .AND. TS_NGASSIM) CALL FATAL_ERROR &
         & ('Using two kinds of TS assimilation does not make any sense!')

    ! SET THE FVCOM RUN MODE BASED ON WHAT TYPE OF ASSIMILATION WE ARE DOING
!    IF(SST_ASSIM .OR. SSTGRD_ASSIM) THEN
!       FVCOM_RUN_MODE = FVCOM_NUDGE_AVG_SST
!    END IF
      
    ! SET THE FVCOM RUN MODE BASED ON WHAT TYPE OF ASSIMILATION WE ARE DOING
!    IF(TS_ASSIM .OR. TSGRD_ASSIM) THEN
!       FVCOM_RUN_MODE = FVCOM_NUDGE_AVG_TSGRD
!    END IF
 
!    IF(SST_ASSIM .OR. SSTGRD_ASSIM) THEN
!      FVCOM_RUN_MODE = FVCOM_NUDGE_ASSIM
!    END IF  

    FVCOM_RUN_MODE = FVCOM_NUDGE_OI_ASSIM

!    IF(CUR_ASSIM) THEN
!       IF(CUR_TIMESCALE * DTI >= 1.0_sp) CALL FATAL_ERROR &
!           ('The e-folding time scale for the cur nudging is to short for the model time step!') 
!    END IF

    IF(CUR_NGASSIM) THEN
       IF(CUR_GALPHA * DTI >= 1.0_sp) CALL FATAL_ERROR &
           ('The e-folding time scale for the cur nudging is to short for the model time step!') 
    END IF
           
    IF(SST_ASSIM) THEN

       IF (SST_TIMESCALE * DTI >= 1.0_sp) CALL FATAL_ERROR &
            ('The efolding time scale for the sst nudging is to short for the model time step!') 

    ELSE IF(SSTGRD_ASSIM) THEN
       IF (SSTGRD_TIMESCALE * DTI >= 1.0_sp) CALL FATAL_ERROR &
            ('The efolding time scale for the sst nudging is to short for the model time step!') 

    END IF

    IF(TS_NGASSIM) THEN
       IF(TS_GALPHA * DTI >= 1.0_sp) CALL FATAL_ERROR &
           ('The e-folding time scale for the TS nudging is to short for the model time step!') 
    END IF
           
!    IF(TS_NGASSIM) THEN
!
!       IF (TS_TIMESCALE * DTI >= 1.0_sp) CALL FATAL_ERROR &
!            ('The efolding time scale for the TS nudging is to short for the model time step!') 
!
!    ELSE IF(TSGRD_ASSIM) THEN
     IF(TSGRD_ASSIM) THEN
       IF (TSGRD_TIMESCALE * DTI >= 1.0_sp) CALL FATAL_ERROR &
            ('The efolding time scale for the TS nudging is to short for the model time step!') 

    END IF


    IF (DBG_SET(DBG_SBR)) WRITE(IPT,*) "END: SET_ASSIM_PARAM"

    RETURN
  END SUBROUTINE SET_ASSIM_PARAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|

!==============================================================================!
  SUBROUTINE SETUP_DATA_ASSIMILATION
!==============================================================================!
    IMPLICIT NONE

    IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) "START SETUP_DATA_ASSIMILATION"

    IF(SST_ASSIM) THEN

       ALLOCATE(ASSIM_SST)
       ASSIM_SST%RADIUS     = SST_RADIUS
       ASSIM_SST%WEIGHT_MAX = SST_WEIGHT_MAX
       ASSIM_SST%TIMESCALE  = SST_TIMESCALE
       ASSIM_SST%TIME_WINDOW= SST_TIME_WINDOW
       
       CALL OPEN_SST_ASSIM(ASSIM_SST%FILE)

       ALLOCATE(ASSIM_SST%OBS_LOC)
       ALLOCATE(ASSIM_SST%OBS_DAT)

       CALL LOAD_SST_ASSIM


    END IF

    IF(SSTGRD_ASSIM) THEN
       CALL OPEN_SSTGRD_ASSIM

       CALL SET_SSTGRD_ASSIM
    END IF

    IF(SSHGRD_ASSIM) THEN
       CALL OPEN_SSHGRD_ASSIM
       CALL SET_SSHGRD_ASSIM
    END IF

    IF(TSGRD_ASSIM) THEN
       CALL OPEN_TSGRD_ASSIM
       CALL SET_TSGRD_ASSIM
       CALL SET_MONTHLY_TS_MEAN(NC_TSAVG_ASSIM)

       IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) "SETUP ASSIMULTION RSTFILE"
       CALL SETUP_RSTFILE_ASSIM
    END IF

!    IF(TS_ASSIM) THEN
!       CALL OPEN_TS_ASSIM(ASSIM_SST%FILE)
!
!       CALL SET_TS_ASSIM
!    END IF

!    IF(CUR_ASSIM) THEN
!       CALL OPEN_CUR_ASSIM(ASSIM_SST%FILE)

!       CALL SET_CUR_ASSIM
!    END IF

!qxu{ 
    IF(TS_NGASSIM .OR. TS_OIASSIM) CALL SET_TS_ASSIM_DATA
    IF(CUR_NGASSIM .OR. CUR_OIASSIM) CALL SET_CUR_ASSIM_DATA
!qxu}

    IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) "END SETUP_DATA_ASSIMILATION"
  END SUBROUTINE SETUP_DATA_ASSIMILATION
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|

!==============================================================================!
  SUBROUTINE OPEN_SST_ASSIM(NCF)
!==============================================================================!
    IMPLICIT NONE
    TYPE(NCFILE), POINTER :: NCF,NCT
    character(len=160) :: pathnfile
    
    pathnfile= trim(INPUT_DIR)//trim(SST_ASSIM_FILE)
    
    
    ! INITIALIZE TYPE TO HOLD FILE METADATA
    
    CALL NC_INIT(NCF,pathnfile)
    
!!$    ! OPEN THE FILE AND LOAD METADATA       
!!$    If(.not. NCF%OPEN) then
!!$       Call NC_OPEN(NCF)
!!$       CALL NC_LOAD(NCF)
!!$       NCT => NCF
!!$       FILEHEAD => ADD(FILEHEAD,NCT)
!!$    end if
    
  END SUBROUTINE OPEN_SST_ASSIM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|

!==============================================================================!
  SUBROUTINE OPEN_SSTGRD_ASSIM
!==============================================================================!
    IMPLICIT NONE
    TYPE(NCFILE), POINTER :: NCF
    integer :: charnum
    logical :: fexist,back
    character(len=160) :: pathnfile
    
    ! INITIALIZE TYPE TO HOLD FILE METADATA
    pathnfile= trim(INPUT_DIR)//trim(SSTGRD_ASSIM_FILE)

    CALL NC_INIT(NCF,pathnfile)
    
    ! OPEN THE FILE AND LOAD METADATA       
    If(.not. NCF%OPEN) then
       Call NC_OPEN(NCF)
       CALL NC_LOAD(NCF)
       FILEHEAD => ADD(FILEHEAD,NCF)
    end if

  END SUBROUTINE OPEN_SSTGRD_ASSIM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|

!==============================================================================!
  SUBROUTINE OPEN_SSHGRD_ASSIM
!==============================================================================!
    IMPLICIT NONE
    TYPE(NCFILE), POINTER :: NCF
    integer :: charnum
    logical :: fexist,back
    character(len=160) :: pathnfile
    
    ! INITIALIZE TYPE TO HOLD FILE METADATA
    pathnfile= trim(INPUT_DIR)//trim(SSHGRD_ASSIM_FILE)

    CALL NC_INIT(NCF,pathnfile)
    
    ! OPEN THE FILE AND LOAD METADATA       
    If(.not. NCF%OPEN) then
       Call NC_OPEN(NCF)
       CALL NC_LOAD(NCF)
       FILEHEAD => ADD(FILEHEAD,NCF)
    end if

  END SUBROUTINE OPEN_SSHGRD_ASSIM

!==============================================================================!
  SUBROUTINE OPEN_TSGRD_ASSIM
!==============================================================================!
    IMPLICIT NONE
    TYPE(NCFILE), POINTER :: NCF
    integer :: charnum
    logical :: fexist,back
    character(len=160) :: pathnfile
    
    ! INITIALIZE TYPE TO HOLD FILE METADATA
    pathnfile= trim(INPUT_DIR)//trim(TSGRD_ASSIM_FILE)

    CALL NC_INIT(NCF,pathnfile)
   
    ! OPEN THE FILE AND LOAD METADATA       
    If(.not. NCF%OPEN) then
       Call NC_OPEN(NCF)
       CALL NC_LOAD(NCF)
       FILEHEAD => ADD(FILEHEAD,NCF)
    end if

  END SUBROUTINE OPEN_TSGRD_ASSIM


!==============================================================================!
!  SUBROUTINE OPEN_CUR_ASSIM(NCF)
!==============================================================================!
!    IMPLICIT NONE
!    TYPE(NCFILE), POINTER :: NCF,NCT
!    character(len=160) :: pathnfile
    
!    pathnfile= trim(INPUT_DIR)//trim(CUR_NGASSIM_FILE)
    
    
!    ! INITIALIZE TYPE TO HOLD FILE METADATA
    
!    CALL NC_INIT(NCF,pathnfile)
    
!!$    ! OPEN THE FILE AND LOAD METADATA       
!!$    If(.not. NCF%OPEN) then
!!$       Call NC_OPEN(NCF)
!!$       CALL NC_LOAD(NCF)
!!$       NCT => NCF
!!$       FILEHEAD => ADD(FILEHEAD,NCT)
!!$    end if


!  END SUBROUTINE OPEN_CUR_ASSIM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|

!==============================================================================!
!  SUBROUTINE OPEN_TS_ASSIM(NCF)
!==============================================================================!
!    IMPLICIT NONE
!    TYPE(NCFILE), POINTER :: NCF,NCT
!    character(len=160) :: pathnfile
    
!    pathnfile= trim(INPUT_DIR)//trim(TS_NGASSIM_FILE)
    
    
!    ! INITIALIZE TYPE TO HOLD FILE METADATA
    
!    CALL NC_INIT(NCF,pathnfile)
    
!!$    ! OPEN THE FILE AND LOAD METADATA       
!!$    If(.not. NCF%OPEN) then
!!$       Call NC_OPEN(NCF)
!!$       CALL NC_LOAD(NCF)
!!$       NCT => NCF
!!$       FILEHEAD => ADD(FILEHEAD,NCT)
!!$    end if

!  END SUBROUTINE OPEN_TS_ASSIM
!==============================================================================|

   SUBROUTINE LOAD_SST_ASSIM
     IMPLICIT NONE

     INTEGER I,TCNT,NCNT
     CHARACTER(LEN=160) :: FNAME
     LOGICAL :: FEXIST
     REAL(SP):: TEMPF
     INTEGER,  ALLOCATABLE, DIMENSION(:) :: ITEMP
     REAL(SP), ALLOCATABLE, DIMENSION(:) :: FTEMP
     REAL(SP):: X0,Y0,T0

     INTEGER J,K,ECNT,ITMP,IOS
     REAL(SP)::DX,DY,RD
       
!------------------------------------------------------------------------------!
       
     FNAME = ASSIM_SST%FILE%FNAME
!
!--Make Sure SST Assimilation Data File Exists---------------------------------!
!
   INQUIRE(FILE=TRIM(FNAME),EXIST=FEXIST)
   IF(.NOT.FEXIST)THEN
      CALL FATAL_ERROR('THE SST OBSERVATION FILE: ',&
           & TRIM(FNAME),&
           &' DOES NOT EXIST')
   END IF
     
!
!--Read Number SST Measurements and Data Set Size------------------------------!
!


   OPEN(1,FILE=TRIM(FNAME),STATUS='OLD') ; REWIND(1)
   READ(1,*) TCNT
   READ(1,*) TEMPF,NCNT


   ASSIM_SST%OBS_LOC%N_TIMES = TCNT
   ASSIM_SST%OBS_LOC%N_LOCS  = NCNT

   ALLOCATE(ASSIM_SST%OBS_LOC%X(NCNT,TCNT))
   ALLOCATE(ASSIM_SST%OBS_LOC%Y(NCNT,TCNT))
   ALLOCATE(ASSIM_SST%OBS_LOC%OBSTIME(TCNT))

   ALLOCATE(ASSIM_SST%OBS_DAT%DATA(NCNT,TCNT))
   ALLOCATE(ASSIM_SST%OBS_DAT%MASK(NCNT,TCNT))

   REWIND(1)  ; READ(1,*)
   
!
!--Read SST Data---------------------------------------------------------------!
!
   DO I=1,TCNT
     READ(1,*) TEMPF,NCNT 
     DO J=1,NCNT
       READ(1,*) X0,Y0,T0
       ASSIM_SST%OBS_LOC%X(J,I)=X0
       ASSIM_SST%OBS_LOC%Y(J,I)=Y0

       ASSIM_SST%OBS_DAT%DATA(J,I)=T0

       ASSIM_SST%OBS_LOC%OBSTIME(I) = SECONDS2TIME(TEMPF)

     END DO
   END DO

!
!--Close SST Observation Data File---------------------------------------------!
!
   CLOSE(1)

!
!--Shift Coordinates-of SST Observation Locations------------------------------!
!
   ASSIM_SST%OBS_LOC%X = ASSIM_SST%OBS_LOC%X - VXMIN
   ASSIM_SST%OBS_LOC%Y = ASSIM_SST%OBS_LOC%Y - VYMIN


   ASSIM_SST%OBS_LOC%OBSTIME = ASSIM_SST%OBS_LOC%OBSTIME + RUNFILE_StartTime


   END SUBROUTINE LOAD_SST_ASSIM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   SUBROUTINE SET_SSTGRD_ASSIM
!------------------------------------------------------------------------------!
!  SET UP ASSIMILATION DATA FOR SST OBSERVATIONS                               |
!------------------------------------------------------------------------------!

   IMPLICIT NONE
   INTEGER I
   REAL(DP) :: TEMP      
   ! SOME NC POINTERS

   TYPE(NCATT), POINTER :: ATT
   TYPE(NCDIM), POINTER :: DIM
   TYPE(NCVAR), POINTER :: VAR
   
   ! SOME HANDY VARIABLES TO PLAY WITH
   LOGICAL FOUND  
   TYPE(TIME) :: FSTART, FEND,TCNT,ONEDAY
   !------------------------------------------------------------------------------!
       
   ! GET THE FILE POINTER
   SST_FILE => FIND_FILE(FILEHEAD,trim(SSTGRD_ASSIM_FILE),FOUND)
   IF(.not. FOUND) CALL FATAL_ERROR &
        & ("COULD NOT FIND SSTGRD_ASSIM_FILE FILE OBJECT",&
        & "FILE NAME: "//TRIM(SSTGRD_ASSIM_FILE))
   
   ! CHECK THE ATTRIBUTES
   ATT => FIND_ATT(SST_FILE,"source",FOUND)
   IF(.not. FOUND) CALL FATAL_ERROR &
        & ("IN SST_ASSIMGRD_FILE FILE OBJECT",&
        & "FILE NAME: "//TRIM(SSTGRD_ASSIM_FILE),&
        &"COULD NOT FIND GLOBAL ATTRIBURE: 'source'")
   
   IF(ATT%CHR(1)(1:5) /= "FVCOM") CALL FATAL_ERROR &
        &('IN SSTGRD_ASSIM_FILE FILE OBJECT',&
        & "FILE NAME: "//TRIM(SSTGRD_ASSIM_FILE),&
        & "THE GLOBAL ATTRIBURE 'source' SHOULD BE 'FVCOM' ")
   
   ! LOOK FOR THE DIMENSIONS 
   DIM => FIND_DIM(SST_FILE,'node',FOUND) ! NODE
   IF(.not. FOUND) CALL FATAL_ERROR &
        & ("IN SSTGRD_ASSIM_FILE OBJECT",&
        & "FILE NAME: "//TRIM(SSTGRD_ASSIM_FILE),&
        &"COULD NOT FIND DIMENSION 'node'")
   
   IF (DIM%DIM /= MGL) CALL FATAL_ERROR&
        & ("THE NUMBER OF NODES IN THE SST GRID ASSIM FILE DOES NOT MATCH THE GRID FILE!")
   
   DIM => FIND_DIM(SST_FILE,'time',FOUND) ! TIME
   IF(.not. FOUND) CALL FATAL_ERROR &
        & ("IN SSTGRD_ASSIM_FILE OBJECT",&
        & "FILE NAME: "//TRIM(SSTGRD_ASSIM_FILE),&
        &"COULD NOT FIND DIMENSION 'time'")
   
   N_TIMES_SST = DIM%DIM ! SET THE VARIABLE VALUE
   
   ONEDAY = days2time(1.0_DP)
   
   ! CHECK THE BEGIN TIME
   FSTART = GET_FILE_TIME(SST_FILE,1)
   
   IF (FSTART > STARTTIME + oneday/2) THEN
      IF(DBG_SET(DBG_LOG)) THEN
         CALL PRINT_REAL_TIME(FSTART,IPT,"FIRST TIME IN SST GRID ASSIM FILE")
         CALL PRINT_REAL_TIME(STARTTIME,IPT,"MODEL START TIME")
      END IF
      CALL FATAL_ERROR("THE SST GRID ASSIM START TIME IS INCORRECT",&
           & "THE MODEL STARTS MORE THAN 12 HOURS BEFORE THE FIRST SST DATA")
   END IF
   
   ! CHECK THE END TIME
   FEND = GET_FILE_TIME(SST_FILE,N_TIMES_SST)
   
   IF (FEND < ENDTIME - oneday/2) THEN
      IF(DBG_SET(DBG_LOG)) THEN
         CALL PRINT_REAL_TIME(FEND,IPT,"LAST TIME IN SST GRID ASSIM FILE")
         CALL PRINT_REAL_TIME(ENDTIME,IPT,"MODEL END TIME")
      END IF
      CALL FATAL_ERROR("THE SST GRID ASSIM END TIME IS INCORRECT",&
           & "THE MODEL ENDS MORE THAN 12 HOURS AFTER THE LAST SST DATA")
   END IF
   
   ! CHECK THE FILE INTERVAL - EXPECTING DAILY DATA
   TCNT = FSTART + ONEDAY *(N_TIMES_SST-1)
   
   IF (.NOT. FEND == TCNT) THEN
      IF(DBG_SET(DBG_LOG)) THEN
         CALL PRINT_REAL_TIME(FEND,IPT,"LAST TIME IN SST GRID ASSIM FILE")
         CALL PRINT_REAL_TIME(TCNT,IPT,"EXPECTED END TIME FOR DAILY DATA")
      END IF
      CALL FATAL_ERROR("THE SST GRID ASSIM FILE TIME STEP IS NOT DAILY",&
           & "THE START TIME PLUS NUMBER OF DAYS IN THE FILE DOES N&
           &OT EQUAL THE END TIME?")
   END IF
   
   ! GET INTERVAL FOR SAVING SST STATE 
   SST_SAVE_INTERVAL = ONEDAY / SSTGRD_N_PER_INTERVAL

   ! ALLOCATE MEMORY
   ALLOCATE(SST_OBS(M))   
   ALLOCATE(SST_AVG(M))
   ALLOCATE(SST_SAVED(M,1:SSTGRD_N_PER_INTERVAL))
# if defined (ENKF)
   ALLOCATE(SST_AVG_ENKF(M,ENKF_NENS))
   ALLOCATE(SST_SAVED_ENKF(M,1:SSTGRD_N_PER_INTERVAL,ENKF_NENS))
# endif

   VAR => FIND_VAR(SST_FILE,'sst',FOUND)
   IF (.NOT.FOUND) CALL FATAL_ERROR &
        ("THE VARIABLE 'sst' IS MISSING FROM THE SSTGRD_ASSIM FILE?")

   CALL NC_CONNECT_AVAR(VAR,SST_OBS)
   VAR_SST =>VAR
   

   ! SET THE FILE STACK COUNT
   ! LOOK FOR TIME WITHIN ONE DTI OF MID-DAY

   FSTART = StartTime + ONEDAY/2

   DO I =1 ,N_TIMES_SST

      TCNT = GET_FILE_TIME(SST_FILE,I)

      ! IF TIME IS BEFORE MID POINT THE FIRST DAY
      IF (TCNT < FSTART - IMDTI) CYCLE
      
      ! IF IT IS ALSO LESS THEN THIS
      IF (TCNT < FSTART + IMDTI) THEN
         ! WE FOUND IT
         SST_FILE%FTIME%NEXT_STKCNT = I
         EXIT

      ELSE ! THRE IS NO VALID TIME IN THE FILE

         CALL FATAL_ERROR &
              & ("THERE IS NO VALID FIRST TIME IN THE SST OBSERVATION FILE")
      END IF
   END DO


   RETURN
   END SUBROUTINE SET_SSTGRD_ASSIM
!==============================================================================|

   SUBROUTINE SET_SSHGRD_ASSIM
!------------------------------------------------------------------------------!
!  SET UP ASSIMILATION DATA FOR SST OBSERVATIONS                               |
!------------------------------------------------------------------------------!

   IMPLICIT NONE
   INTEGER I
   REAL(DP) :: TEMP      
   ! SOME NC POINTERS

   TYPE(NCATT), POINTER :: ATT
   TYPE(NCDIM), POINTER :: DIM
   TYPE(NCVAR), POINTER :: VAR
   
   ! SOME HANDY VARIABLES TO PLAY WITH
   LOGICAL FOUND  
   TYPE(TIME) :: FSTART, FEND,TCNT,ONEDAY
   !------------------------------------------------------------------------------!
       
   ! GET THE FILE POINTER
   SSH_FILE => FIND_FILE(FILEHEAD,trim(SSHGRD_ASSIM_FILE),FOUND)
   IF(.not. FOUND) CALL FATAL_ERROR &
        & ("COULD NOT FIND SSHGRD_ASSIM_FILE FILE OBJECT",&
        & "FILE NAME: "//TRIM(SSHGRD_ASSIM_FILE))
   
   ! CHECK THE ATTRIBUTES
   ATT => FIND_ATT(SSH_FILE,"source",FOUND)
   IF(.not. FOUND) CALL FATAL_ERROR &
        & ("IN SSH_ASSIMGRD_FILE FILE OBJECT",&
        & "FILE NAME: "//TRIM(SSHGRD_ASSIM_FILE),&
        &"COULD NOT FIND GLOBAL ATTRIBURE: 'source'")
   
   IF(ATT%CHR(1)(1:5) /= "FVCOM") CALL FATAL_ERROR &
        &('IN SSHGRD_ASSIM_FILE FILE OBJECT',&
        & "FILE NAME: "//TRIM(SSHGRD_ASSIM_FILE),&
        & "THE GLOBAL ATTRIBURE 'source' SHOULD BE 'FVCOM' ")
   
   ! LOOK FOR THE DIMENSIONS 
   DIM => FIND_DIM(SSH_FILE,'node',FOUND) ! NODE
   IF(.not. FOUND) CALL FATAL_ERROR &
        & ("IN SSHGRD_ASSIM_FILE OBJECT",&
        & "FILE NAME: "//TRIM(SSHGRD_ASSIM_FILE),&
        &"COULD NOT FIND DIMENSION 'node'")
   
   IF (DIM%DIM /= MGL) CALL FATAL_ERROR&
        & ("THE NUMBER OF NODES IN THE SSH GRID ASSIM FILE DOES NOT MATCH THE GRID FILE!")
   
   DIM => FIND_DIM(SSH_FILE,'time',FOUND) ! TIME
   IF(.not. FOUND) CALL FATAL_ERROR &
        & ("IN SSHGRD_ASSIM_FILE OBJECT",&
        & "FILE NAME: "//TRIM(SSHGRD_ASSIM_FILE),&
        &"COULD NOT FIND DIMENSION 'time'")
   
   N_TIMES_SSH = DIM%DIM ! SET THE VARIABLE VALUE
   
   ONEDAY = days2time(1.0_DP)
   
   ! CHECK THE BEGIN TIME
   FSTART = GET_FILE_TIME(SSH_FILE,1)
   
   IF (FSTART > STARTTIME + oneday/2) THEN
      IF(DBG_SET(DBG_LOG)) THEN
         CALL PRINT_REAL_TIME(FSTART,IPT,"FIRST TIME IN SSH GRID ASSIM FILE")
         CALL PRINT_REAL_TIME(STARTTIME,IPT,"MODEL START TIME")
      END IF
      CALL FATAL_ERROR("THE SSH GRID ASSIM START TIME IS INCORRECT",&
           & "THE MODEL STARTS MORE THAN 12 HOURS BEFORE THE FIRST SSH DATA")
   END IF
   
   ! CHECK THE END TIME
   FEND = GET_FILE_TIME(SSH_FILE,N_TIMES_SSH)
   
   IF (FEND < ENDTIME - oneday/2) THEN
      IF(DBG_SET(DBG_LOG)) THEN
         CALL PRINT_REAL_TIME(FEND,IPT,"LAST TIME IN SSH GRID ASSIM FILE")
         CALL PRINT_REAL_TIME(ENDTIME,IPT,"MODEL END TIME")
      END IF
      CALL FATAL_ERROR("THE SSH GRID ASSIM END TIME IS INCORRECT",&
           & "THE MODEL ENDS MORE THAN 12 HOURS AFTER THE LAST SSH DATA")
   END IF
   
   ! CHECK THE FILE INTERVAL - EXPECTING DAILY DATA
   TCNT = FSTART + ONEDAY *(N_TIMES_SSH-1)

   IF (.NOT. FEND == TCNT) THEN
      IF(DBG_SET(DBG_LOG)) THEN
         CALL PRINT_REAL_TIME(FEND,IPT,"LAST TIME IN SSH GRID ASSIM FILE")
         CALL PRINT_REAL_TIME(TCNT,IPT,"EXPECTED END TIME FOR DAILY DATA")
      END IF
      CALL FATAL_ERROR("THE SSH GRID ASSIM FILE TIME STEP IS NOT DAILY",&
           & "THE START TIME PLUS NUMBER OF DAYS IN THE FILE DOES N&
           &OT EQUAL THE END TIME?")
   END IF
   
   ! GET INTERVAL FOR SAVING SST STATE 
   SSH_SAVE_INTERVAL = ONEDAY / SSHGRD_N_PER_INTERVAL

   ! ALLOCATE MEMORY
   ALLOCATE(SSH_OBS(M))   
   ALLOCATE(SSH_AVG(M))
   ALLOCATE(SSH_SAVED(M,0:SSHGRD_N_PER_INTERVAL))
#  if defined (ENKF)
   ALLOCATE(SSH_AVG_ENKF(M,ENKF_NENS))
   ALLOCATE(SSH_SAVED_ENKF(M,1:SSTGRD_N_PER_INTERVAL,ENKF_NENS))
#  endif   

   VAR => FIND_VAR(SSH_FILE,'ssh',FOUND)
   IF (.NOT.FOUND) CALL FATAL_ERROR &
        ("THE VARIABLE 'ssh' IS MISSING FROM THE SSTGRD_ASSIM FILE?")

   CALL NC_CONNECT_AVAR(VAR,SSH_OBS)
   VAR_SSH =>VAR
   

   ! SET THE FILE STACK COUNT
   ! LOOK FOR TIME WITHIN ONE DTI OF MID-DAY

   FSTART = StartTime + ONEDAY/2

   DO I =1 ,N_TIMES_SSH

      TCNT = GET_FILE_TIME(SSH_FILE,I)

      ! IF TIME IS BEFORE MID POINT THE FIRST DAY
      IF (TCNT < FSTART - IMDTI) CYCLE
      
      ! IF IT IS ALSO LESS THEN THIS
      IF (TCNT < FSTART + IMDTI) THEN
         ! WE FOUND IT
         SSH_FILE%FTIME%NEXT_STKCNT = I
         EXIT

      ELSE ! THRE IS NO VALID TIME IN THE FILE

         CALL FATAL_ERROR &
              & ("THERE IS NO VALID FIRST TIME IN THE SSH OBSERVATION FILE")
      END IF
   END DO


   RETURN
 END SUBROUTINE SET_SSHGRD_ASSIM

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   SUBROUTINE SET_TSGRD_ASSIM
!------------------------------------------------------------------------------!
!  SET UP ASSIMILATION DATA FOR TS OBSERVATIONS                               |
!------------------------------------------------------------------------------!

   IMPLICIT NONE
   INTEGER I
   REAL(DP) :: TEMP      
   ! SOME NC POINTERS

   TYPE(NCATT), POINTER :: ATT
   TYPE(NCDIM), POINTER :: DIM
   TYPE(NCVAR), POINTER :: VAR
   REAL(SP), POINTER :: STORAGE_ARR(:,:)
   
   ! SOME HANDY VARIABLES TO PLAY WITH
   LOGICAL FOUND  
   TYPE(TIME) :: FSTART, FEND,TCNT,ONEDAY
   INTEGER :: ERROR, STATUS
   !------------------------------------------------------------------------------!
       
   ! GET THE FILE POINTER
   TS_FILE => FIND_FILE(FILEHEAD,trim(TSGRD_ASSIM_FILE),FOUND)
   IF(.not. FOUND) CALL FATAL_ERROR &
        & ("COULD NOT FIND TSGRD_ASSIM_FILE FILE OBJECT",&
        & "FILE NAME: "//TRIM(TSGRD_ASSIM_FILE))
   
   ! CHECK THE ATTRIBUTES
   ATT => FIND_ATT(TS_FILE,"source",FOUND)
   IF(.not. FOUND) CALL FATAL_ERROR &
        & ("IN TS_ASSIMGRD_FILE FILE OBJECT",&
        & "FILE NAME: "//TRIM(TSGRD_ASSIM_FILE),&
        &"COULD NOT FIND GLOBAL ATTRIBURE: 'source'")
   
   IF(ATT%CHR(1)(1:5) /= "FVCOM") CALL FATAL_ERROR &
        &('IN TSGRD_ASSIM_FILE FILE OBJECT',&
        & "FILE NAME: "//TRIM(TSGRD_ASSIM_FILE),&
        & "THE GLOBAL ATTRIBURE 'source' SHOULD BE 'FVCOM' ")
   
   ! LOOK FOR THE DIMENSIONS 
   DIM => FIND_DIM(TS_FILE,'node',FOUND) ! NODE
   IF(.not. FOUND) CALL FATAL_ERROR &
        & ("IN TSGRD_ASSIM_FILE OBJECT",&
        & "FILE NAME: "//TRIM(SSTGRD_ASSIM_FILE),&
        &"COULD NOT FIND DIMENSION 'node'")
   
   IF (DIM%DIM /= MGL) CALL FATAL_ERROR&
        & ("THE NUMBER OF NODES IN THE SST GRID ASSIM FILE DOES NOT MATCH THE GRID FILE!")
   
   
   DIM => FIND_DIM(TS_FILE,'time',FOUND) ! TIME
   IF(.not. FOUND) CALL FATAL_ERROR &
        & ("IN TSGRD_ASSIM_FILE OBJECT",&
        & "FILE NAME: "//TRIM(SSTGRD_ASSIM_FILE),&
        &"COULD NOT FIND DIMENSION 'time'")
   
   N_TIMES_TSC = DIM%DIM ! SET THE VARIABLE VALUE
   
   ONEDAY = days2time(1.0_DP)
   
!   ! CHECK THE BEGIN TIME
!   FSTART = GET_FILE_TIME(TS_FILE,1)
!   
!   IF (FSTART > STARTTIME + oneday/2) THEN
!      IF(DBG_SET(DBG_LOG)) THEN
!         CALL PRINT_REAL_TIME(FSTART,IPT,"FIRST TIME IN TSC GRID ASSIM FILE")
!         CALL PRINT_REAL_TIME(STARTTIME,IPT,"MODEL START TIME")
!      END IF
!      CALL FATAL_ERROR("THE TS GRID ASSIM START TIME IS INCORRECT",&
!           & "THE MODEL STARTS MORE THAN 12 HOURS BEFORE THE FIRST TS DATA")
!   END IF
   
   ! CHECK THE END TIME
!   FEND = GET_FILE_TIME(TS_FILE,N_TIMES_TSC)
   
!   IF (FEND < ENDTIME - oneday/2) THEN
!      IF(DBG_SET(DBG_LOG)) THEN
!         CALL PRINT_REAL_TIME(FEND,IPT,"LAST TIME IN TS GRID ASSIM FILE")
!         CALL PRINT_REAL_TIME(ENDTIME,IPT,"MODEL END TIME")
!      END IF
!      CALL FATAL_ERROR("THE TSC  GRID ASSIM END TIME IS INCORRECT",&
!           & "THE MODEL ENDS MORE THAN 12 HOURS AFTER THE LAST TS DATA")
!   END IF
   
   ! CHECK THE FILE INTERVAL - EXPECTING DAILY DATA
!   TCNT = FSTART + ONEDAY *(N_TIMES_TSC-1)
!   
!   IF (.NOT. FEND == TCNT) THEN
!      IF(DBG_SET(DBG_LOG)) THEN
!         CALL PRINT_REAL_TIME(FEND,IPT,"LAST TIME IN TS GRID ASSIM FILE")
!         CALL PRINT_REAL_TIME(TCNT,IPT,"EXPECTED END TIME FOR DAILY DATA")
!      END IF
!      CALL FATAL_ERROR("THE TS GRID ASSIM FILE TIME STEP IS NOT DAILY",&
!           & "THE START TIME PLUS NUMBER OF DAYS IN THE FILE DOES N&
!           &OT EQUAL THE END TIME?")
!   END IF
    
   !TSGRD_N_PER_INTERVAL 
 
   ! GET INTERVAL FOR SAVING TS STATE 
   TSC_SAVE_INTERVAL = ONEDAY / TSGRD_N_PER_INTERVAL

!   IF(DBG_SET(DBG_LOG)) WRITE(IPT,*)"TSGRD_N_PER_INTERVAL:",TSGRD_N_PER_INTERVAL
!   IF(DBG_SET(DBG_LOG)) WRITE(IPT,*)"Memory:",   MEMCNT ! MT*KB*4*NDB 
   ! ALLOCATE MEMORY

!   TSGRD_N_PER_INTERVAL=TSGRD_N_PER_INTERVAL*Pmdays

!!  IF not allocate !
     !allocate here
   IF (ALLOCATED(TC_OBS)) DEALLOCATE(TC_OBS,STAT=ERROR)
   IF (ALLOCATED(TC_AVG)) DEALLOCATE(TC_AVG,STAT=ERROR)
   IF (ALLOCATED(TC0_AVG)) DEALLOCATE(TC0_AVG,STAT=ERROR)

   ALLOCATE(TC_OBS(M,KBM1))   
   ALLOCATE(TC_AVG(M,KBM1))
   ALLOCATE(TC0_AVG(M,KBM1))   ; TC0_AVG = 0.0_SP

   ALLOCATE(TC_SAVED(M,KBM1,1:TSGRD_N_PER_INTERVAL))

   MEMCNT = M*KB*3 + M*KB*TSGRD_N_PER_INTERVAL + MEMCNT
   IF(DBG_SET(DBG_LOG)) WRITE(IPT,*)"Memory:",   MEMCNT ! MT*KB*4*NDB 
   IF(DBG_SET(DBG_LOG)) WRITE(IPT,*)"TSGRD_N_PER_INTERVAL:",TSGRD_N_PER_INTERVAL
   
   IF (ALLOCATED(SC_OBS)) DEALLOCATE(SC_OBS,STAT=ERROR)
   IF (ALLOCATED(SC_AVG)) DEALLOCATE(SC_AVG,STAT=ERROR)
   IF (ALLOCATED(SC0_AVG)) DEALLOCATE(SC0_AVG,STAT=ERROR)

   ALLOCATE(SC_OBS(M,KBM1))   
   ALLOCATE(SC_AVG(M,KBM1))
   ALLOCATE(SC0_AVG(M,KBM1))  ; SC0_AVG = 0.0_SP

   ALLOCATE(SC_SAVED(M,KBM1,1:TSGRD_N_PER_INTERVAL))

   MEMCNT = M*KB*3 + M*KB*TSGRD_N_PER_INTERVAL + MEMCNT
   IF(DBG_SET(DBG_LOG)) WRITE(IPT,*)"Memory:",   MEMCNT ! MT*KB*4*NDB

   NULLIFY(TSC_T1_N, TSC_T1_P, TSC_S1_N, TSC_S1_P)

   VAR => FIND_VAR(TS_FILE,'temp_clim',FOUND)
   IF (.NOT.FOUND) CALL FATAL_ERROR &
        ("THE VARIABLE 'temp_clim' IS MISSING FROM THE TSGRD_ASSIM FILE?")

   ! MAKE SPACE FOR THE DATA FROM THE FILE
   TSC_T1_N => reference_var(var)
   ALLOCATE(STORAGE_ARR(1:M,KBM1), stat = status)
   IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN TSGRD_ASSIM FORCING")
   CALL NC_CONNECT_PVAR(TSC_T1_N,STORAGE_ARR)
   NULLIFY(STORAGE_ARR)

   ! MAKE SPACE FOR THE DATA FROM THE FILE
   TSC_T1_P => reference_var(var)
   ALLOCATE(STORAGE_ARR(1:M,KBM1), stat = status)
   IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN TSGRD_ASSIM FORCING")
   CALL NC_CONNECT_PVAR(TSC_T1_P,STORAGE_ARR)
   NULLIFY(STORAGE_ARR)

!   CALL NC_CONNECT_AVAR(VAR,TC_OBS)
!   VAR_TC =>VAR

   VAR => FIND_VAR(TS_FILE,'salinity_clim',FOUND)
   IF (.NOT.FOUND) CALL FATAL_ERROR &
        ("THE VARIABLE 'salinity_clim' IS MISSING FROM THE SSTGRD_ASSIM FILE?")

   ! MAKE SPACE FOR THE DATA FROM THE FILE
   TSC_S1_N => reference_var(var)
   ALLOCATE(STORAGE_ARR(1:M,KBM1), stat = status)
   IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN OFFLINE FORCING")
   CALL NC_CONNECT_PVAR(TSC_S1_N,STORAGE_ARR)
   NULLIFY(STORAGE_ARR)

   ! MAKE SPACE FOR THE DATA FROM THE FILE
   TSC_S1_P => reference_var(var)
   ALLOCATE(STORAGE_ARR(1:M,KBM1), stat = status)
   IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN OFFLINE FORCING")
   CALL NC_CONNECT_PVAR(TSC_S1_P,STORAGE_ARR)
   NULLIFY(STORAGE_ARR)

!   CALL NC_CONNECT_AVAR(VAR,SC_OBS)
!   VAR_SC =>VAR

   RETURN
   END SUBROUTINE SET_TSGRD_ASSIM
!==============================================================================|
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!   SUBROUTINE SET_TSGRD_ASSIM_julian
!------------------------------------------------------------------------------!
!  SET UP ASSIMILATION DATA FOR TS OBSERVATIONS                               |
!------------------------------------------------------------------------------!
!   USE ALL_VARS
!   USE MOD_PAR 

!   IMPLICIT NONE
!   INTEGER I
!   REAL(DP) :: TEMP      
!   ! SOME NC POINTERS

!   TYPE(NCATT), POINTER :: ATT
!   TYPE(NCDIM), POINTER :: DIM
!   TYPE(NCVAR), POINTER :: VAR
   
!   ! SOME HANDY VARIABLES TO PLAY WITH
!   LOGICAL FOUND  
!   TYPE(TIME) :: FSTART, FEND,TCNT,ONEDAY
!   INTEGER :: ERROR
   !------------------------------------------------------------------------------!
   
!!   ONEDAY = days2time(1.0_DP)
   
!   ! GET INTERVAL FOR SAVING TS STATE 
!   !TSC_SAVE_INTERVAL = ONEDAY / TSGRD_N_PER_INTERVAL

!   IF(DBG_SET(DBG_LOG)) WRITE(IPT,*)"TSGRD_N_PER_INTERVAL:",TSGRD_N_PER_INTERVAL
!   IF(DBG_SET(DBG_LOG)) WRITE(IPT,*)"Memory:",   MEMCNT ! MT*KB*4*NDB 
!   ! ALLOCATE MEMORY
   
  
!   TSC_SAVE_TIME =ASSIM_RESTART_TIME-days2time(real(Pmdays,kind=DP))/2.0
!   !TSGRD_N_PER_INTERVAL=TSGRD_N_PER_INTERVAL*Pmdays

!   TC_AVG=0.0_SP
!   SC_AVG=0.0_SP

!   RETURN
!   END SUBROUTINE SET_TSGRD_ASSIM_julian
!==============================================================================|

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!   SUBROUTINE SET_TS_ASSIM
!------------------------------------------------------------------------------!
!!     IMPLICIT NONE

  
!!   END SUBROUTINE SET_TS_ASSIM
!==============================================================================|

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!   SUBROUTINE SET_CUR_ASSIM
!------------------------------------------------------------------------------!
!!     IMPLICIT NONE

  
!!   END SUBROUTINE SET_CUR_ASSIM
!==============================================================================|


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
!==============================================================================!
   SUBROUTINE SET_CUR_ASSIM_DATA 
!==============================================================================!

!------------------------------------------------------------------------------!
!  SET UP ASSIMILATION DATA FOR CURRENT OBSERVATIONS                           |
!------------------------------------------------------------------------------!
   IMPLICIT NONE
   INTEGER I,J,K,NN,ECNT,ITMP,NCNT,IOS,NLAY
   CHARACTER(LEN=120) :: FNAME,ONAME
   CHARACTER(LEN= 2 ) :: NAC   
   INTEGER,  ALLOCATABLE, DIMENSION(:) :: ITEMP
   REAL(SP), ALLOCATABLE, DIMENSION(:) :: FTEMP
   REAL(SP):: X0,Y0,DX,DY,RD,SIGMA_C,ISOBATH_ANGLE,D_ANGLE,ANG_OBS_SIM,DIR_WEIGHT
   REAL(SP), PARAMETER :: ALF = 0.05_SP
   LOGICAL :: FEXIST
   INTEGER :: MAXEL,NBD_CNT
   REAL(SP) :: CUR_RADIUS 
   INTEGER :: JJ

   REAL(SP), DIMENSION(3) :: XTRI,YTRI 
   REAL(SP) :: RDLAST 
   INTEGER :: LOCIJ(2),MIN_LOC,IERR,Nsite_tmp 
   INTEGER :: ND1,ND2,ND3 
   REAL(SP) :: DELTA,COFA,COFB,COFC 
   REAL(SP) ::S11,S22,S33,RTMP,RRTMP 
   REAL(SP), DIMENSION(1:NGL,1) :: RDLIST 
   REAL(SP), DIMENSION(KB) :: ZZ_OB 
   REAL(SP), ALLOCATABLE :: ZZ_G(:,:)   
   REAL(SP) :: X0C,Y0C,HX,HY,H0,TEMP_H
   REAL(DP) :: ODAYS

   REAL(SP) :: X11_TMP, X22_TMP, X33_TMP

   IF (DBG_SET(DBG_SBR)) WRITE(IPT,*) "START: SET_CUR_ASSIM_DATA"

!------------------------------------------------------------------------------!
!  Read Number of Current Observations and Coordinates of Each                 !
!------------------------------------------------------------------------------!
    
!   FNAME = "./"//TRIM(INPUT_DIR)//"/"//trim(casename)//"_cur.xy"
   IF(CUR_NGASSIM .AND. .NOT. CUR_OIASSIM)THEN
     FNAME = "./"//TRIM(INPUT_DIR)//"/"//trim(CUR_NGASSIM_FILE)//".xy"
   ELSE IF(CUR_OIASSIM .AND. .NOT. CUR_NGASSIM)THEN
     FNAME = "./"//TRIM(INPUT_DIR)//"/"//trim(CUR_OIASSIM_FILE)//".xy"
   ELSE
     CALL FATAL_ERROR("CUR_NGASSIM and CUR_OIASSIM cannot be both true or false")
   END IF  
!
!--Make Sure Current Assimilation Data File Exists-----------------------------!
!
   INQUIRE(FILE=TRIM(FNAME),EXIST=FEXIST)
   IF(MSR .AND. .NOT.FEXIST)THEN
     WRITE(IPT,*)'CURRENT OBSERVATION FILE: ',FNAME,' DOES NOT EXIST'
     WRITE(IPT,*)'HALTING.....'
     CALL PSTOP
   END IF
     
!
!--Read Number of Current Measurement Stations---------------------------------!
!
   OPEN(1,FILE=TRIM(FNAME),STATUS='OLD')
   READ(1,*) N_ASSIM_CUR

   ALLOCATE(CUR_OBS(N_ASSIM_CUR))

!
!--Read X,Y Coordinate of Measurement Stations---------------------------------!
!

   DO I=1,N_ASSIM_CUR
     READ(1,*)ITMP,CUR_OBS(I)%X,CUR_OBS(I)%Y,CUR_OBS(I)%DEPTH,NLAY,CUR_OBS(I)%SITA
     CUR_OBS(I)%N_LAYERS = NLAY

     ALLOCATE(CUR_OBS(I)%ODEPTH(NLAY))
     DO J=1,NLAY
       READ(1,*)CUR_OBS(I)%ODEPTH(J)
       IF(CUR_OBS(I)%ODEPTH(J) > CUR_OBS(I)%DEPTH)THEN
         IF(MSR)WRITE(IPT,*)'CURRENT OBSERVATION LAYER',J,'OF CURRENT MOORING',I
         IF(MSR)WRITE(IPT,*)'EXCEEDS BATHYMETRIC DEPTH=',CUR_OBS(I)%ODEPTH(J),CUR_OBS(I)%DEPTH
         IF(MSR)WRITE(IPT,*)'HALTING...........'
         CALL PSTOP
       END IF
     END DO
   END DO
   CUR_MAX_LAYER = MAXVAL(CUR_OBS(1:N_ASSIM_CUR)%N_LAYERS)

!
!--Shift Coordinates-----------------------------------------------------------!
!
# if !defined (SPHERICAL) 
   CUR_OBS(:)%X = CUR_OBS(:)%X - VXMIN 
   CUR_OBS(:)%Y = CUR_OBS(:)%Y - VYMIN 
# endif

!# if defined (SPHERICAL) 
   IF(SERIAL)THEN
     XCG = XC
     YCG = YC
     XG = VX
     YG = VY
   END IF
#  if defined(MULTIPROCESSOR)
   IF(PAR)CALL ACOLLECT(MYID,MSRID,NPROCS,EMAP,XC,XCG)
   IF(PAR)CALL ACOLLECT(MYID,MSRID,NPROCS,EMAP,YC,YCG)
   IF(PAR)CALL ACOLLECT(MYID,MSRID,NPROCS,NMAP,VX,XG)
   IF(PAR)CALL ACOLLECT(MYID,MSRID,NPROCS,NMAP,VY,YG)
   IF(PAR)CALL MPI_BCAST(XCG,NGL+1,MPI_F,0,MPI_COMM_WORLD,IERR)
   IF(PAR)CALL MPI_BCAST(YCG,NGL+1,MPI_F,0,MPI_COMM_WORLD,IERR)
   IF(PAR)CALL MPI_BCAST(XG,MGL+1,MPI_F,0,MPI_COMM_WORLD,IERR)
   IF(PAR)CALL MPI_BCAST(YG,MGL+1,MPI_F,0,MPI_COMM_WORLD,IERR)
#  endif
!# endif
!  


!
!--find cell number (CUR_OBS(J)%N_CELL) of Obs.
!
   RRTMP = 100000.0     !100km
!   IF(CUR_NGASSIM)THEN
!         RRTMP = CUR_NG_RADIUS
!   ELSE IF(CUR_OIASSIM)THEN
!        RRTMP = CUR_OI_RADIUS
!   ELSE
!         WRITE(IPT,*) 'IN SUBROUTINE SET_CUR_ASSIM_DATA, ONE OF CUR_NGASSIM AND &
!                       CUR_OIASSIM MUST BE DEFINED AS TRUE.'
!         WRITE(IPT,*) 'HALTING.....'            
!         CALL PSTOP
!   END IF

   DO J= 1,N_ASSIM_CUR 
     X0 = CUR_OBS(J)%X
     Y0 = CUR_OBS(J)%Y
     DO I=1,NGL
#    if !defined (SPHERICAL)
       Rtmp = SQRT((XCG(I)-X0)*(XCG(I)-X0)+(YCG(I)-Y0)*(YCG(I)-Y0))
#    else
       CALL ARC(X0,Y0,XCG(I),YCG(I),Rtmp)
#    endif
       IF(Rtmp < RRTMP) THEN
         X11_TMP = XG(NVG(I,1))-X0
         X22_TMP = XG(NVG(I,2))-X0
         X33_TMP = XG(NVG(I,3))-X0
#        if defined (SPHERICAL)
         IF(X11_TMP >  180.0_SP)THEN
           X11_TMP = -360.0_SP+X11_TMP
         ELSEIF(X11_TMP < -180.0_SP)THEN
           X11_TMP =  360.0_SP+X11_TMP
         END IF
         IF(X22_TMP >  180.0_SP)THEN
           X22_TMP = -360.0_SP+X22_TMP
         ELSEIF(X22_TMP < -180.0_SP)THEN
           X22_TMP =  360.0_SP+X22_TMP
         END IF
         IF(X33_TMP >  180.0_SP)THEN
           X33_TMP = -360.0_SP+X33_TMP
         ELSEIF(X33_TMP < -180.0_SP)THEN
           X33_TMP =  360.0_SP+X33_TMP
         END IF
#        endif
         S11 = X22_TMP*(YG(NVG(I,3))-Y0)-X33_TMP*(YG(NVG(I,2))-Y0)
         S22 = X33_TMP*(YG(NVG(I,1))-Y0)-X11_TMP*(YG(NVG(I,3))-Y0)
         S33 = X11_TMP*(YG(NVG(I,2))-Y0)-X22_TMP*(YG(NVG(I,1))-Y0)
         IF(S11 <= 0._SP .AND. S22 <= 0._SP .AND. S33<= 0._SP) THEN
           CUR_OBS(J)%N_CELL = I
           GOTO 300
         ELSE
           CUR_OBS(J)%N_CELL = 0
         ENDIF
       ELSE
         CUR_OBS(J)%N_CELL = -1
       ENDIF
     ENDDO
300  CONTINUE
!     IF(MSR) WRITE(200,*) j,X0,Y0,CUR_OBS(J)%N_CELL
     IF(CUR_OBS(J)%N_CELL <= 0) THEN
       IF(MSR) WRITE(IPT,*)'ERROR--CURRENT OBS SITE:',J,' OUT OF DOMAN',&
               CUR_OBS(J)%N_CELL 
       CALL PSTOP
     ENDIF 
   ENDDO
             
!--Gather AW0G,AWXG & AWYG use for interp grid to Obs station
   IF(.NOT. ALLOCATED(AW0G))   ALLOCATE(AW0G(0:NGL,3)) 
   IF(.NOT. ALLOCATED(AWXG))   ALLOCATE(AWXG(0:NGL,3))
   IF(.NOT. ALLOCATED(AWYG))   ALLOCATE(AWYG(0:NGL,3))

   IF(SERIAL)THEN
     AW0G = AW0
     AWXG = AWX
     AWYG = AWY
   END IF
# if defined(MULTIPROCESSOR)
   IF(PAR)CALL ACOLLECT(MYID,MSRID,NPROCS,EMAP,AW0,AW0G)
   IF(PAR)CALL ACOLLECT(MYID,MSRID,NPROCS,EMAP,AWX,AWXG)
   IF(PAR)CALL ACOLLECT(MYID,MSRID,NPROCS,EMAP,AWY,AWYG)
   IF(PAR)CALL MPI_BCAST(AW0G,(NGL+1)*3,MPI_F,0,MPI_COMM_WORLD,IERR)
   IF(PAR)CALL MPI_BCAST(AWXG,(NGL+1)*3,MPI_F,0,MPI_COMM_WORLD,IERR)
   IF(PAR)CALL MPI_BCAST(AWYG,(NGL+1)*3,MPI_F,0,MPI_COMM_WORLD,IERR)
# endif
!
!--Close Current Observation Global File---------------------------------------!
!
   CLOSE(1)   

   IF(CUR_OIASSIM)THEN
!------------------------------------------------------------------------------!
!  Read Correlation Length of Current Observations                             !
!------------------------------------------------------------------------------!
       
!JQI   FNAME = "./"//TRIM(INPDIR)//"/"//trim(casename)//"_radius_cur.dat"

!JQI   INQUIRE(FILE=TRIM(FNAME),EXIST=FEXIST)
!JQI   IF(MSR .AND. .NOT.FEXIST)THEN
!JQI     WRITE(IPT,*)'CURRENT OBSERVATION FILE: ',FNAME,' DOES NOT EXIST'
!JQI     WRITE(IPT,*)'HALTING.....'
!JQI     CALL PSTOP
!JQI   END IF
     
!JQI   OPEN(1,FILE=TRIM(FNAME),STATUS='OLD')

   ALLOCATE(CUR_PARAM(2,N_ASSIM_CUR))

!JQI   DO I=1,N_ASSIM_CUR
!JQI     READ(1,*)PARAM_CUR(1,I),PARAM_CUR(2,I)
!JQI   END DO

!JQI   CLOSE(1)
     CUR_PARAM = CUR_OI_RADIUS       !30000.0_SP
   END IF
       
!
!------------------------------------------------------------------------------!
!  Open Current Observation Files for Each Observation Point and Read Data     !
!------------------------------------------------------------------------------!
   IF(CUR_NGASSIM .AND. .NOT. CUR_OIASSIM)THEN
     ONAME = "./"//TRIM(INPUT_DIR)//"/"//trim(CUR_NGASSIM_FILE)//".dat"
   ELSE IF(CUR_OIASSIM .AND. .NOT. CUR_NGASSIM)THEN
     ONAME = "./"//TRIM(INPUT_DIR)//"/"//trim(CUR_OIASSIM_FILE)//".dat"
   ELSE
     CALL FATAL_ERROR("CUR_NGASSIM and CUR_OIASSIM cannot be both true or false")
   END IF  
   INQUIRE(FILE=TRIM(ONAME),EXIST=FEXIST)
   IF(MSR .AND. .NOT.FEXIST)THEN
     WRITE(IPT,*)'CURRENT OBSERVATION FILE: ',ONAME,' DOES NOT EXIST'
     WRITE(IPT,*)'HALTING.....'
     CALL PSTOP
   END IF
   OPEN(1,FILE=ONAME,STATUS='old')  ; REWIND(1)

   DO I=1,N_ASSIM_CUR
!----Count Number of Data Entries for Observation I----------------------------!
!     DO WHILE(.TRUE.)
!       READ(1,*,IOSTAT=IOS)
!       IF(IOS < 0)EXIT
!       NCNT = NCNT + 1
!     END DO

     READ(1,*,IOSTAT=IOS) Nsite_tmp,NCNT     
     IF(IOS<0) then
       WRITE(IPT,*) 'ERROR in read ',trim(casename),'_cur.dat at site number:',I
       CALL PSTOP
     END IF
     CUR_OBS(I)%N_TIMES = NCNT
     
     IF(NCNT == 0)THEN
       IF(MSR)WRITE(IPT,*)'NO DATA FOR CURRENT OBSERVATION',I
       CALL PSTOP
     END IF

!----Allocate Arrays to Hold Current (UA,VA) and Time (TIME)-------------------!
     ALLOCATE(CUR_OBS(I)%TIMES(CUR_OBS(I)%N_TIMES))
     ALLOCATE(CUR_OBS(I)%UO( CUR_OBS(I)%N_TIMES , CUR_OBS(I)%N_LAYERS ))
     ALLOCATE(CUR_OBS(I)%VO( CUR_OBS(I)%N_TIMES , CUR_OBS(I)%N_LAYERS ))

!----Read in Current Data for Observation I------------------------------------!
     NLAY = CUR_OBS(I)%N_LAYERS
     DO J=1,CUR_OBS(I)%N_TIMES
       READ(1,*,IOSTAT=IOS) ODAYS,(CUR_OBS(I)%UO(J,K),CUR_OBS(I)%VO(J,K),K=1,NLAY)
       IF(IOS < 0)THEN         
         WRITE(IPT,*)'ERROR in read ',trim(casename),'_cur.dat at site:',I,'Time No:',J
         CALL PSTOP
       END IF  
       CUR_OBS(I)%TIMES(J) = DAYS2TIME(ODAYS)
     END DO

!----Convert Time to Seconds---------------------------------------------------!
!----Shift Jan 1 Based Time Data to Dec 1 Based Time Data-----CASESPECIFIC-----!
!     IF(trim(CASENAME) == 'gom')THEN
!       CUR_OBS(I)%TIMES = (CUR_OBS(I)%TIMES-1.0_SP)*24.0_SP*3600.0_SP !after 1996 
!     ELSE   
!       CUR_OBS(I)%TIMES = CUR_OBS(I)%TIMES*3600.0_SP*24.0_SP
!     END IF
 
!----Convert Current Data from cm/s to m/s-------------------------------------!
     CUR_OBS(I)%UO = CUR_OBS(I)%UO * .01_SP
     CUR_OBS(I)%VO = CUR_OBS(I)%VO * .01_SP
  END DO
  CLOSE(1)

!------------------------------------------------------------------------------!
!  Count Number of Points with Bad Data (UO = 0. + V0 = 0.)         
!------------------------------------------------------------------------------!
  NBD_CNT = 0
  DO I=1,N_ASSIM_CUR
    DO J=1,CUR_OBS(I)%N_TIMES
      DO K=1,CUR_OBS(I)%N_LAYERS
        IF(ABS(CUR_OBS(I)%UO(J,K)) + ABS(CUR_OBS(I)%VO(J,K)) < .0001) THEN
          NBD_CNT = NBD_CNT + 1
        END IF
      END DO
    END DO
  END DO
 
!------------------------------------------------------------------------------!
!  Compute Spatial Interpolation Weights for each Mooring Location 
!------------------------------------------------------------------------------!
!   DO I=1,N_ASSIM_CUR
!   LMIN = 100000000.
!     X0 = CUR_OBS(I)%X
!     Y0 = CUR_OBS(I)%Y
!     DO J=1,MGL
!       DX = ABS(XG(J)-X0)
!       DY = ABS(YG(J)-Y0)
!       IF(SQRT(DX**2 + DY**2) < LMIN)THEN
!         LMIN = SQRT(DX**2 + DY**2)--dbg=7
!         JMIN = J
!       END IF
!      END DO
!      CUR_OBS(I)%SITA = SITA_GD(JMIN) + 3.14159_SP/2.0_SP
!    END DO

   ALLOCATE(ITEMP(NGL),FTEMP(NGL),DA_CUR(NGL))     ; DA_CUR = 0
   DO I=1,N_ASSIM_CUR
     X0 = CUR_OBS(I)%X
     Y0 = CUR_OBS(I)%Y
     ISOBATH_ANGLE = CUR_OBS(I)%SITA*DEG2RAD       !/180.0_SP*3.1415926_SP
     ECNT = 0

     DO J=1,NGL
#      if defined (SPHERICAL)
       CALL ARCX(X0,Y0,XCG(J),YCG(J),DX)
       CALL ARCY(X0,Y0,XCG(J),YCG(J),DY)
       CALL ARC(X0,Y0,XCG(J),YCG(J),RD)
#      else
       DX = ABS(XCG(J)-X0)
       DY = ABS(YCG(J)-Y0)
       RD = SQRT(DX**2 + DY**2)
#      endif

       IF(CUR_NGASSIM)THEN
         CUR_RADIUS = CUR_NG_RADIUS
       ELSE IF(CUR_OIASSIM)THEN
         CUR_RADIUS = CUR_OI_RADIUS
       ELSE
         WRITE(IPT,*) 'IN SUBROUTINE SET_CUR_ASSIM_DATA, ONE OF CUR_NGASSIM AND ', &
                       'CUR_OIASSIM MUST BE DEFINED AS TRUE.'
         WRITE(IPT,*) 'HALTING.....'            
         CALL PSTOP
       END IF

       IF(RD <= CUR_RADIUS)THEN         
         DA_CUR(J)   = 1
         ECNT        = ECNT + 1      
         ITEMP(ECNT) =  J
         FTEMP(ECNT) = (CUR_RADIUS**2 - RD**2) / (CUR_RADIUS**2 + RD**2)

         IF(ISOBATH_ANGLE==0.0_SP)THEN
           DIR_WEIGHT = 1.0_SP
         ELSE
           ANG_OBS_SIM = ATAN2(DY,DX)
           D_ANGLE     = ANG_OBS_SIM - ISOBATH_ANGLE 
           D_ANGLE     = D_ANGLE - INT(D_ANGLE/PI)*PI      !3.1415926_SP)*3.1415926_SP
           D_ANGLE     = ABS(D_ANGLE)
!           DIR_WEIGHT  = (ABS(D_ANGLE-0.5*3.1415926_SP)+ALF*3.1415926_SP)/ &
!                        ((0.5_SP+ALF)*3.1415926_SP)
           DIR_WEIGHT  = (ABS(D_ANGLE-0.5*PI)+ALF*PI)/((0.5_SP+ALF)*PI)
           IF(J == CUR_OBS(J)%N_CELL) DIR_WEIGHT = 1.0_SP
         END IF
! qxu       FTEMP(ECNT) = FTEMP(ECNT)*DIR_WEIGHT
       END IF
     END DO
     IF(ECNT == 0)THEN
       WRITE(IPT,*)'ERROR SETTING UP CURRENT DATA ASSIMILATION'
       WRITE(IPT,*)'NO ELEMENTS LIE WITHIN RADIUS',CUR_RADIUS
       WRITE(IPT,*)'OF OBSERVATION POINT',I
       CALL PSTOP   
     ELSE
       CUR_OBS(I)%N_INTPTS = ECNT
       ALLOCATE(CUR_OBS(I)%INTPTS(ECNT))
       ALLOCATE(CUR_OBS(I)%X_WEIGHT(ECNT))
       CUR_OBS(I)%INTPTS(1:ECNT)  = ITEMP(1:ECNT)
       CUR_OBS(I)%X_WEIGHT(1:ECNT) = FTEMP(1:ECNT)
     END IF
   END DO
   DEALLOCATE(FTEMP,ITEMP)


     
!------------------------------------------------------------------------------!
!  Compute Sigma Layer Weights for Vertical Interpolation                      
!------------------------------------------------------------------------------!
   ALLOCATE(ZZ_G(0:MGL,KB))
   IF(.NOT. ALLOCATED(HG)) ALLOCATE(HG(0:MGL))
   IF(SERIAL)THEN
     ZZ_G = ZZ
     HG = H
   END IF  
#  if defined (MULTIPROCESSOR)   
   IF(PAR)CALL ACOLLECT(MYID,MSRID,NPROCS,NMAP,ZZ,ZZ_G) 
   IF(PAR)CALL ACOLLECT(MYID,MSRID,NPROCS,NMAP,H,HG) 

   IF(PAR)CALL MPI_BCAST(ZZ_G,(MGL+1)*KB,MPI_F,0,MPI_COMM_WORLD,IERR)
   IF(PAR)CALL MPI_BCAST(HG,MGL+1,MPI_F,0,MPI_COMM_WORLD,IERR)
#  endif
   DO I=1,N_ASSIM_CUR
     NLAY = CUR_OBS(I)%N_LAYERS
     ALLOCATE(CUR_OBS(I)%S_INT(NLAY,2))
     ALLOCATE(CUR_OBS(I)%S_WEIGHT(NLAY,2))
     X0 = CUR_OBS(I)%X
     Y0 = CUR_OBS(I)%Y
# if defined (SPHERICAL)
     DO NN=1,NGL
       CALL ARC(X0,Y0,XCG(NN),YCG(NN),RDLIST(NN,1))
     END DO
# else
     RDLIST(1:NGL,1) = SQRT((XCG(1:NGL)-X0)**2 + (YCG(1:NGL)-Y0)**2)
# endif
     RDLAST = -1.0_SP
in:  DO WHILE(.TRUE.)
       LOCIJ = MINLOC(RDLIST,RDLIST>RDLAST)
       MIN_LOC = LOCIJ(1)
       IF(MIN_LOC == 0)THEN
         EXIT in
       END IF
       XTRI = XG(NVG(MIN_LOC,1:3))
       YTRI = YG(NVG(MIN_LOC,1:3))
       RDLAST = RDLIST(MIN_LOC,1)
       IF(ISINTRIANGLE1(XTRI,YTRI,X0,Y0))THEN
         JJ = MIN_LOC
         EXIT IN
       END IF
       RDLAST = RDLIST(MIN_LOC,1)
     END DO IN

     ND1=NVG(JJ,1)
     ND2=NVG(JJ,2)
     ND3=NVG(JJ,3)

# if defined (SPHERICAL)
     CALL ARCX(X0,Y0,XCG(JJ),YCG(JJ),X0C)
     CALL ARCX(X0,Y0,XCG(JJ),YCG(JJ),Y0C)
# else
     X0C = X0 - XCG(JJ)
     Y0C = Y0 - YCG(JJ)
# endif
     !----Linear Interpolation of Bathymetry------------------------------------------!
     H0 = AW0G(JJ,1)*HG(ND1)+AW0G(JJ,2)*HG(ND2)+AW0G(JJ,3)*HG(ND3)
     HX = AWXG(JJ,1)*HG(ND1)+AWXG(JJ,2)*HG(ND2)+AWXG(JJ,3)*HG(ND3)
     HY = AWYG(JJ,1)*HG(ND1)+AWYG(JJ,2)*HG(ND2)+AWYG(JJ,3)*HG(ND3)
     TEMP_H = H0 + HX*X0C + HY*Y0C !here is the interpolated observation depth

     DO K=1,KBM1
      !----Linear Interpolation of Bathymetry------------------------------------------!
      H0 = AW0G(JJ,1)*ZZ_G(ND1,K)*HG(ND1)+AW0G(JJ,2)*ZZ_G(ND2,K)*HG(ND2)+AW0G(JJ,3)*ZZ_G(ND3,K)*HG(ND3)
      HX = AWXG(JJ,1)*ZZ_G(ND1,K)*HG(ND1)+AWXG(JJ,2)*ZZ_G(ND2,K)*HG(ND2)+AWXG(JJ,3)*ZZ_G(ND3,K)*HG(ND3)
      HY = AWYG(JJ,1)*ZZ_G(ND1,K)*HG(ND1)+AWYG(JJ,2)*ZZ_G(ND2,K)*HG(ND2)+AWYG(JJ,3)*ZZ_G(ND3,K)*HG(ND3)
      ZZ_OB(K) = (H0 + HX*X0C + HY*Y0C)/TEMP_H
     END DO  

     DO J=1,NLAY
       SIGMA_C = -CUR_OBS(I)%ODEPTH(J)/TEMP_H
       DO K=2,KBM1
         IF(ZZ_OB(K) <= SIGMA_C .AND. ZZ_OB(K-1) > SIGMA_C)THEN 
           CUR_OBS(I)%S_INT(J,1) = K-1
           CUR_OBS(I)%S_INT(J,2) = K
           CUR_OBS(I)%S_WEIGHT(J,1) = (SIGMA_C-ZZ_OB(K))/(ZZ_OB(K-1)-ZZ_OB(K))
           CUR_OBS(I)%S_WEIGHT(J,2) = 1.0_SP - CUR_OBS(I)%S_WEIGHT(J,1) 
         END IF  
       END DO
       IF(ZZ_OB(1) <= SIGMA_C)THEN  !!OBSERVATION ABOVE CENTROID OF FIRST SIGMA LAYER
         CUR_OBS(I)%S_INT(J,1) = 1
         CUR_OBS(I)%S_INT(J,2) = 1
         CUR_OBS(I)%S_WEIGHT(J,1) = 1.0_SP
         CUR_OBS(I)%S_WEIGHT(J,2) = 0.0_SP
       END IF
       IF(ZZ_OB(KBM1) >= SIGMA_C)THEN !!OBSERVATION BELOW CENTROID OF BOTTOM SIGMA LAYER
         CUR_OBS(I)%S_INT(J,1) = KBM1
         CUR_OBS(I)%S_INT(J,2) = KBM1
         CUR_OBS(I)%S_WEIGHT(J,1) = 1.0_SP
         CUR_OBS(I)%S_WEIGHT(J,2) = 0.0_SP
       END IF
       if (msr) print *, "obs",i,"lay", j,SIGMA_C,CUR_OBS(I)%S_INT(J,1),CUR_OBS(I)%S_INT(J,2), CUR_OBS(I)%ODEPTH(J),CUR_OBS(I)%DEPTH,'scoor',zz_ob(1:kbm1)
     END DO
   END DO
!#  if defined (MULTIPROCESSOR)   
   DEALLOCATE(ZZ_G)  
!#  endif
!------------------------------------------------------------------------------!
!  Report Number of Interpolation Points, Location and Number of Data 
!------------------------------------------------------------------------------!
   IF(.NOT. MSR)RETURN

   WRITE(IPT,*)
   WRITE(IPT,*)'!            CURRENT OBSERVATION DATA           '
# if defined (SPHERICAL)
   WRITE(IPT,*)" MOORING#    LON       LAT    #INTERP PTS  #DATA TIMES  NEAR_CELL  SITA"
# else
   WRITE(IPT,*)" MOORING#   X(KM)      Y(KM)  #INTERP PTS  #DATA TIMES  NEAR_CELL  SITA"
# endif
   DO I=1,N_ASSIM_CUR
     MAXEL = MAXLOC(CUR_OBS(I)%X_WEIGHT,DIM=1)
     WRITE(IPT,'(2X,I5,3X,F8.1,3X,F8.1,3X,I6,5X,I6,5X,I6,5X,F8.1)') &
# if defined (SPHERICAL)
     I,CUR_OBS(I)%X,CUR_OBS(I)%Y, &
# else
     I,(CUR_OBS(I)%X+VXMIN)/1000.,(CUR_OBS(I)%Y+VYMIN)/1000., &
# endif
       CUR_OBS(I)%N_INTPTS,CUR_OBS(I)%N_TIMES,CUR_OBS(I)%INTPTS(MAXEL),&
       CUR_OBS(I)%SITA
   END DO
   WRITE(IPT,*)
   WRITE(IPT,*)'NUMBER OF BAD CURRENT DATA POINTS: ',NBD_CNT
   WRITE(IPT,*)" MOORING #   BEGIN TIME  END TIME"
   DO I=1,N_ASSIM_CUR
!   WRITE(IPT,*)I,CUR_OBS(I)%TIMES(1)/(24.*3600.),&
!       CUR_OBS(I)%TIMES(CUR_OBS(I)%N_TIMES)/(24.*3600.)
   WRITE(IPT,*)I,CUR_OBS(I)%TIMES(1)%MJD,&
       CUR_OBS(I)%TIMES(CUR_OBS(I)%N_TIMES)%MJD
   END DO
   IF (DBG_SET(DBG_SBR)) WRITE(IPT,*) "END: SET_CUR_ASSIM_DATA"

 
   RETURN
   END SUBROUTINE SET_CUR_ASSIM_DATA

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|

   SUBROUTINE CURRENT_ASSIMILATION
!==============================================================================|
!  USE CURRENT OBSERVATION DATA TO ADJUST VELOCITY COMPONENTS                  |
!==============================================================================|
   IMPLICIT NONE
   REAL(SP), ALLOCATABLE, DIMENSION(:,:) :: UINT,VINT,UCORR,VCORR,UG,VG,TWGHT
   REAL(SP), ALLOCATABLE, DIMENSION(:,:) :: UCORR1,VCORR1
   REAL(SP), ALLOCATABLE, DIMENSION(:)   :: FTEMP
   REAL(SP) :: WEIGHT,DEFECT,CORRECTION,DT_MIN,SIMTIME,T_THRESH,WGHT,TOT_WGHT
   REAL(SP) :: U1,U2,V1,V2,W1,W2,WEIGHT1,WEIGHT2
   TYPE(TIME) ::DIFTIME
   INTEGER I,J,K,J1,K1,K2,NLAY,ITIMEA,NTIME,IERR
   INTEGER COUNT
   INTRINSIC MINLOC

   INTEGER IP,JN1,JN2,JN3
   REAL(SP) XSC,YSC,COFT0,COFTX,COFTY
!==============================================================================|

   IF (DBG_SET(DBG_SBR)) WRITE(IPT,*) "START: CURRENT_ASSIMILATION" 

!------------------------------------------------------------------------------!
!  Check the Temporal and Spatial OBS Availability                             ! 
!------------------------------------------------------------------------------!
   CUR_OBS%T_WEIGHT = 0. 
   IF(CUR_NGASSIM)THEN
     T_THRESH = CUR_NG_ASTIME_WINDOW   
   ELSE IF(CUR_OIASSIM)THEN  
     T_THRESH = CUR_OI_ASTIME_WINDOW   
   ELSE
     WRITE(IPT,*) "EITHER CUR_NGASSIM OR CUR_OIASSIM MUST BE SET TRUE"
     CALL PSTOP
   END IF       
!   SIMTIME = TIME*86400
   SIMTIME = DTI*FLOAT(IINT)

   COUNT = 0
   DO I=1,N_ASSIM_CUR       
     NTIME = CUR_OBS(I)%N_TIMES
     ALLOCATE(FTEMP(NTIME)) 
!    FTEMP(1:NTIME) = ABS(SIMTIME - CUR_OBS(I)%TIMES(1:NTIME))
     DO j=1,NTIME
        FTEMP(J) = SECONDS(IntTime - CUR_OBS(I)%TIMES(J))
        FTEMP(J) = ABS(FTEMP(J))
     ENDDO
     DT_MIN = MINVAL(FTEMP(1:NTIME))
     CUR_OBS(I)%N_T_WEIGHT = MINLOC(FTEMP,DIM=1)

     IF(DT_MIN < T_THRESH)THEN     
       IF(DT_MIN < .5_SP*T_THRESH) THEN
         CUR_OBS(I)%T_WEIGHT = 1.0_SP
       ELSE
         CUR_OBS(I)%T_WEIGHT = (T_THRESH-DT_MIN)/T_THRESH*2.0_SP
       END IF
       COUNT = COUNT + 1
     END IF

     DEALLOCATE(FTEMP)
   END DO
   IF(COUNT==0)RETURN  ! No OBS data 

!------------------------------------------------------------------------------!
!  Gather U and V Fields to Master Processor                                   ! 
!------------------------------------------------------------------------------!
   ALLOCATE(UG(0:NGL,KB))    ;UG = 0.0_SP
   ALLOCATE(VG(0:NGL,KB))    ;VG = 0.0_SP
#  if defined (MULTIPROCESSOR)
   IF(PAR)THEN
!     CALL GATHER(LBOUND(UF,1),  UBOUND(UF,1),  N,NGL,KB,MYID,NPROCS,EMAP,UF,UG)
!     CALL GATHER(LBOUND(VF,1),  UBOUND(VF,1),  N,NGL,KB,MYID,NPROCS,EMAP,VF,VG)
     CALL ACOLLECT(MYID,MSRID,NPROCS,EMAP,UF,UG)
     CALL ACOLLECT(MYID,MSRID,NPROCS,EMAP,VF,VG)
   END IF
#  endif
   IF(SERIAL)THEN
     UG(1:NGL,1:KBM1) = UF(1:NGL,1:KBM1)
     VG(1:NGL,1:KBM1) = VF(1:NGL,1:KBM1)
   END IF
!------------------------------------------------------------------------------!
!  Calculate Temporal Weight of Measurement (I) at Time(TIME)                  ! 
!------------------------------------------------------------------------------!

   IF(MSR)THEN
!!JQI   CUR_OBS%T_WEIGHT = 0. 
!!JQI   T_THRESH         = CUR_NG_ASTIME_WINDOW   
!   SIMTIME          = TIME*86400
!!JQI   SIMTIME          = DTI*FLOAT(IINT)

!!JQI   DO I=1,N_ASSIM_CUR       
!!JQI     NTIME = CUR_OBS(I)%N_TIMES
!!JQI     ALLOCATE(FTEMP(NTIME)) 
!     FTEMP(1:NTIME) = ABS(SIMTIME - CUR_OBS(I)%TIMES(1:NTIME))
!!JQI     DO J=1,NTIME
!!JQI        FTEMP(J) = SECONDS(IntTime - CUR_OBS(I)%TIMES(J))
!!JQI        FTEMP(J) = ABS(FTEMP(J))
!!JQI     ENDDO
!!JQI     DT_MIN = MINVAL(FTEMP(1:NTIME))
!!JQI     CUR_OBS(I)%N_T_WEIGHT = MINLOC(FTEMP,DIM=1)

!!JQI     IF(DT_MIN < T_THRESH)THEN     
!!JQI       IF(DT_MIN < .5_SP*T_THRESH) THEN
!!JQI         CUR_OBS(I)%T_WEIGHT = 1.0_SP
!!JQI       ELSE
!!JQI         CUR_OBS(I)%T_WEIGHT = (T_THRESH-DT_MIN)/T_THRESH*2.0_SP
!!JQI       END IF
!!JQI     END IF

!!JQI     DEALLOCATE(FTEMP)
!!JQI   END DO
   
       
!------------------------------------------------------------------------------!
!  Interpolate Simulation Data to Local Observation Point                      ! 
!------------------------------------------------------------------------------!
       
   ALLOCATE(UINT(N_ASSIM_CUR,CUR_MAX_LAYER)) ; UINT = 0. 
   ALLOCATE(VINT(N_ASSIM_CUR,CUR_MAX_LAYER)) ; VINT = 0.

   DO I=1,N_ASSIM_CUR  
     DO J=1,CUR_OBS(I)%N_INTPTS
       J1        = CUR_OBS(I)%INTPTS(J)
       WGHT      = CUR_OBS(I)%X_WEIGHT(J)
       NLAY      = CUR_OBS(I)%N_LAYERS
       DO K=1,NLAY
         U1 = UG(J1,CUR_OBS(I)%S_INT(K,1))
         U2 = UG(J1,CUR_OBS(I)%S_INT(K,2))
         V1 = VG(J1,CUR_OBS(I)%S_INT(K,1))
         V2 = VG(J1,CUR_OBS(I)%S_INT(K,2))
         W1 = CUR_OBS(I)%S_WEIGHT(K,1)
         W2 = CUR_OBS(I)%S_WEIGHT(K,2)
         UINT(I,K) = UINT(I,K) + (U1*W1 + U2*W2)*WGHT 
         VINT(I,K) = VINT(I,K) + (V1*W1 + V2*W2)*WGHT 
       END DO
     END DO
     TOT_WGHT = SUM(CUR_OBS(I)%X_WEIGHT(1:CUR_OBS(I)%N_INTPTS))
     UINT(I,1:NLAY) = UINT(I,1:NLAY)/TOT_WGHT
     VINT(I,1:NLAY) = VINT(I,1:NLAY)/TOT_WGHT
   END DO
     
!------------------------------------------------------------------------------!
!  Compute Local Correction by Interpolating Observed/Computed Defect          ! 
!------------------------------------------------------------------------------!

   ALLOCATE(TWGHT(NGL,KBM1))   ; TWGHT = 0.
   ALLOCATE(UCORR(NGL,KBM1))   ; UCORR   = 0.
   ALLOCATE(VCORR(NGL,KBM1))   ; VCORR   = 0.

   IF(CUR_NGASSIM)THEN
    DO I=1,N_ASSIM_CUR
     DO J=1,CUR_OBS(I)%N_INTPTS
       J1     = CUR_OBS(I)%INTPTS(J)
       ITIMEA  = CUR_OBS(I)%N_T_WEIGHT
       NLAY   = CUR_OBS(I)%N_LAYERS
       DO K=1,NLAY
         K1           = CUR_OBS(I)%S_INT(K,1)
         K2           = CUR_OBS(I)%S_INT(K,2)
         W1           = CUR_OBS(I)%S_WEIGHT(K,1)
         W2           = CUR_OBS(I)%S_WEIGHT(K,2)
         WEIGHT1      = CUR_OBS(I)%T_WEIGHT*CUR_OBS(I)%X_WEIGHT(J)*W1
         WEIGHT2      = CUR_OBS(I)%T_WEIGHT*CUR_OBS(I)%X_WEIGHT(J)*W2
         TWGHT(J1,K1) = TWGHT(J1,K1) + WEIGHT1   
         TWGHT(J1,K2) = TWGHT(J1,K2) + WEIGHT2   
         DEFECT       = CUR_OBS(I)%UO(ITIMEA,K) - UINT(I,K)
         CORRECTION   = CUR_GAMA*DEFECT

!         UCORR(J1,K1) = UCORR(J1,K1) + CORRECTION*WEIGHT1
!         UCORR(J1,K2) = UCORR(J1,K2) + CORRECTION*WEIGHT2
         UCORR(J1,K1) = UCORR(J1,K1) + CORRECTION*WEIGHT1*WEIGHT1
         UCORR(J1,K2) = UCORR(J1,K2) + CORRECTION*WEIGHT2*WEIGHT2

         DEFECT       = CUR_OBS(I)%VO(ITIMEA,K) - VINT(I,K)
         CORRECTION   = CUR_GAMA*DEFECT
!         VCORR(J1,K1) = VCORR(J1,K1) + CORRECTION*WEIGHT1
!         VCORR(J1,K2) = VCORR(J1,K2) + CORRECTION*WEIGHT2
         VCORR(J1,K1) = VCORR(J1,K1) + CORRECTION*WEIGHT1*WEIGHT1
         VCORR(J1,K2) = VCORR(J1,K2) + CORRECTION*WEIGHT2*WEIGHT2

!        GEOFF NEW
         IF(ABS(CUR_OBS(I)%UO(ITIMEA,K)) + ABS(CUR_OBS(I)%VO(ITIMEA,K)) < .0001)THEN
           TWGHT(J1,K1) = 0.
           TWGHT(J1,K2) = 0.
         END IF
       END DO
     END DO
    END DO
   
!------------------------------------------------------------------------------!
!  Nudge Simulation Data Using Local Corrections                               ! 
!------------------------------------------------------------------------------!

    DO I=1,NGL
     DO K=1,KBM1
       IF(DA_CUR(I) == 1 .AND. TWGHT(I,K) > 1.0E-08)THEN
         UG(I,K) = UG(I,K) + DTI*CUR_GALPHA*UCORR(I,K)/TWGHT(I,K)
         VG(I,K) = VG(I,K) + DTI*CUR_GALPHA*VCORR(I,K)/TWGHT(I,K)
       END IF
     END DO
    END DO

    DEALLOCATE(TWGHT,UCORR,VCORR,UINT,VINT)
   ELSE IF(CUR_OIASSIM)THEN
     DO I=1,N_ASSIM_CUR
       ITIMEA  = CUR_OBS(I)%N_T_WEIGHT
       NLAY   = CUR_OBS(I)%N_LAYERS
       DO K=1,NLAY
         K1           = CUR_OBS(I)%S_INT(K,1)
         K2           = CUR_OBS(I)%S_INT(K,2)
         W1           = CUR_OBS(I)%S_WEIGHT(K,1)
         W2           = CUR_OBS(I)%S_WEIGHT(K,2)
         WEIGHT1      = CUR_OBS(I)%T_WEIGHT*W1
         WEIGHT2      = CUR_OBS(I)%T_WEIGHT*W2
         TWGHT(I,K1) = TWGHT(I,K1) + WEIGHT1   
         TWGHT(I,K2) = TWGHT(I,K2) + WEIGHT2   
         DEFECT       = CUR_OBS(I)%UO(ITIMEA,K) - UINT(I,K)
         CORRECTION   = DEFECT
         UCORR(I,K1) = UCORR(I,K1) + CORRECTION*WEIGHT1
         UCORR(I,K2) = UCORR(I,K2) + CORRECTION*WEIGHT2
         DEFECT       = CUR_OBS(I)%VO(ITIMEA,K) - VINT(I,K)
         CORRECTION   = DEFECT
         VCORR(I,K1) = VCORR(I,K1) + CORRECTION*WEIGHT1
         VCORR(I,K2) = VCORR(I,K2) + CORRECTION*WEIGHT2
         IF(ABS(CUR_OBS(I)%UO(ITIMEA,K)) + ABS(CUR_OBS(I)%VO(ITIMEA,K)) < .0001)THEN
           TWGHT(I,K1) = 0.
           TWGHT(I,K2) = 0.
         END IF
       END DO
     END DO

     DO I=1,N_ASSIM_CUR
       DO K=1,KBM1
         IF(TWGHT(I,K) > 1.0E-8)THEN
           UCORR(I,K)=UCORR(I,K)/TWGHT(I,K)
	   VCORR(I,K)=VCORR(I,K)/TWGHT(I,K)
         END IF
       END DO
     END DO
       	 
!------------------------------------------------------------------------------!
!  Simulation Data Using Local Corrections                               ! 
!------------------------------------------------------------------------------!

     ALLOCATE(UCORR1(NGL,KBM1));UCORR1=0.0_SP
     ALLOCATE(VCORR1(NGL,KBM1));VCORR1=0.0_SP
   
     DO K=1,KBM1
       CALL CUR_OPTIMINTERP(UCORR(:,K),VCORR(:,K),UCORR1(:,K),VCORR1(:,K))
     END DO

     DO I=1,NGL
       DO K=1,KBM1
         IF(DA_CUR(I) == 1)THEN
!          UG(I,K) = UG(I,K)+UCORR1(I,K)
!          VG(I,K) = VG(I,K)+VCORR1(I,K)
           UG(I,K) = UG(I,K)+ CUR_OIGALPHA*UCORR1(I,K)
           VG(I,K) = VG(I,K)+ CUR_OIGALPHA*VCORR1(I,K)
         END IF
       END DO
     END DO

     DEALLOCATE(TWGHT,UCORR,VCORR,UINT,VINT)
     DEALLOCATE(UCORR1,VCORR1)
   ELSE
     WRITE(IPT,*)'EITHER CUR_NGASSIM OR CUR_OIASSIM SHOULD BE SET AS TRUE'
     CALL PSTOP
   END IF
   END IF  !!MASTER
!------------------------------------------------------------------------------!
!  Disperse New Data Fields to Slave Processors                                ! 
!------------------------------------------------------------------------------!
   IF(SERIAL)THEN
     UF(1:N,1:KBM1) = UG(1:N,1:KBM1)
     VF(1:N,1:KBM1) = VG(1:N,1:KBM1)
   END IF
#  if defined (MULTIPROCESSOR) 
!QXU{
   CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!QXU}
   CALL MPI_BCAST(UG,(NGL+1)*KB,MPI_F,0,MPI_COMM_WORLD,IERR)
   CALL MPI_BCAST(VG,(NGL+1)*KB,MPI_F,0,MPI_COMM_WORLD,IERR)
!QXU{
   CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!QXU}

   IF(PAR)THEN
     DO I=1,N
       UF(I,1:KBM1) = UG(EGID(I),1:KBM1)
       VF(I,1:KBM1) = VG(EGID(I),1:KBM1)
     END DO
   END IF
#  endif

   DEALLOCATE(UG,VG)
          
   IF (DBG_SET(DBG_SBR)) WRITE(IPT,*) "END: CURRENT_ASSIMILATION" 
     
   RETURN
   END SUBROUTINE CURRENT_ASSIMILATION

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|

!==============================================================================|
!==============================================================================|
   SUBROUTINE CUR_OPTIMINTERP(F1,F2,FI1,FI2)
   USE MOD_OPTIMAL_INTERPOLATION
   IMPLICIT NONE

!------------------------------------------------------------------------------|
!  xi(1,:) and xi(2,:) represent the x and y coordindate of the grid of the    |
!  interpolated field                                                          |
!  fi and vari are the interpolated field and its error variance resp.         |
!------------------------------------------------------------------------------|
   REAL(SP) :: XI(2,NGL),FI1(NGL),FI2(NGL),VARI(NGL)

!------------------------------------------------------------------------------|
!  x(1,:) and x(2,:) represent the x and y coordindate of the observations     |
!  f and var are observations and their error variance resp.                   |
!------------------------------------------------------------------------------|
   REAL(SP) :: X(2,N_ASSIM_CUR),VAR(N_ASSIM_CUR),F1(N_ASSIM_CUR),F2(N_ASSIM_CUR)

!------------------------------------------------------------------------------|
!  param: inverse of the correlation length                                    |
!------------------------------------------------------------------------------|
   REAL(SP) :: PARAM(2,N_ASSIM_CUR)

   INTEGER  :: I,J,MM

   IF (DBG_SET(DBG_SBR)) WRITE(IPT,*) "START: CUR_OPTIMINTERP" 
!------------------------------------------------------------------------------|
!  create a regular 2D grid                                                    |
!------------------------------------------------------------------------------|
   DO I=1,NGL
     XI(1,I) = XCG(I)
     XI(2,I) = YCG(I)
   END DO	
   
!------------------------------------------------------------------------------|   
!  param is the inverse of the correlation length                              |
!------------------------------------------------------------------------------|
   PARAM = 1.0_SP/CUR_PARAM 
   MM = CUR_N_INFLU
 
!------------------------------------------------------------------------------| 
!  the error variance of the observations                                      |
!------------------------------------------------------------------------------|
   VAR = 0.0_SP   

   DO I=1,N_ASSIM_CUR
     X(1,I) = CUR_OBS(I)%X
     X(2,I) = CUR_OBS(I)%Y
   END DO  

!------------------------------------------------------------------------------|
!  fi is the interpolated function and vari its error variance                 |
!------------------------------------------------------------------------------|
   CALL OPTIMINTERP(X,F1,VAR,PARAM,MM,XI,FI1,VARI)
   CALL OPTIMINTERP(X,F2,VAR,PARAM,MM,XI,FI2,VARI)

   IF (DBG_SET(DBG_SBR)) WRITE(IPT,*) "END: CUR_OPTIMINTERP" 
   RETURN
   END SUBROUTINE CUR_OPTIMINTERP
!==============================================================================|
!==============================================================================|

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|


   SUBROUTINE SET_TS_ASSIM_DATA 

!------------------------------------------------------------------------------!
!  SET UP ASSIMILATION DATA FOR TEMP/SAL OBSERVATIONS                          |
!------------------------------------------------------------------------------!
   IMPLICIT NONE
   INTEGER I,J,K,ECNT,ITMP,NCNT,IOS,NLAY
   CHARACTER(LEN=120) :: FNAME,ONAME
   CHARACTER(LEN= 2 ) :: NAC   
   INTEGER,  ALLOCATABLE, DIMENSION(:) :: ITEMP
   REAL(SP), ALLOCATABLE, DIMENSION(:) :: FTEMP
   REAL(SP):: X0,Y0,DX,DY,RD,SIGMA_C,ISOBATH_ANGLE,D_ANGLE,ANG_OBS_SIM,DIR_WEIGHT
   REAL(SP), PARAMETER :: ALF = 0.05_SP
   LOGICAL :: FEXIST
   INTEGER :: MAXEL,NBD_CNT
   REAL(SP) :: TS_RADIUS 
   INTEGER :: JMIN 
   REAL(SP):: LMIN

   REAL(SP), DIMENSION(1:NGL,1) :: RDLIST
   REAL(SP), DIMENSION(3) :: XTRI,YTRI
   REAL(SP) :: RDLAST
   INTEGER :: LOCIJ(2),MIN_LOC,JJ,IERR,Nsite_tmp
   INTEGER :: ND1,ND2,ND3,NN
   REAL(SP) :: DELTA,COFA,COFB,COFC
   REAL(SP) ::S11,S22,S33,RTMP,RRTMP
   REAL(SP), DIMENSION(KB) :: ZZ_OB
!#  if defined (MULTIPROCESSOR)   
   REAL(SP), ALLOCATABLE :: ZZ_G(:,:)  
!# endif
   REAL(SP) :: X0C,Y0C,HX,HY,H0,TEMP_H
   REAL(DP) :: ODAYS

   REAL(SP) :: X11_TMP, X22_TMP, X33_TMP
      
   IF(DBG_SET(DBG_SBR)) write(ipt,*) "START: SET_TS_ASSIM_DATA" 

!------------------------------------------------------------------------------!
!  Read Number of Scalar Observations and Coordinates of Each                  !
!------------------------------------------------------------------------------!
       
!   FNAME = "./"//TRIM(INPUT_DIR)//"/"//trim(casename)//"_ts.xy"
   IF(TS_NGASSIM .AND. .NOT. TS_OIASSIM)THEN
     FNAME = "./"//TRIM(INPUT_DIR)//"/"//trim(TS_NGASSIM_FILE)//".xy"
   ELSE IF(TS_OIASSIM .AND. .NOT. TS_NGASSIM)THEN
     FNAME = "./"//TRIM(INPUT_DIR)//"/"//trim(TS_OIASSIM_FILE)//".xy"
   ELSE
     CALL FATAL_ERROR("TS_NGASSIM and TS_OIASSIM cannot be both true or false")
   END IF  
!
!--Make Sure Temperature and Salinity Assimilation Data File Exists-----------------------------!
!
   INQUIRE(FILE=TRIM(FNAME),EXIST=FEXIST)
   IF(MSR .AND. .NOT.FEXIST)THEN
     WRITE(IPT,*)'TEMP/SALINITY OBSERVATION FILE: ',FNAME,' DOES NOT EXIST'
     WRITE(IPT,*)'HALTING.....'
     CALL PSTOP
   END IF
     
!
!--Read Number of T/S Measurement Stations---------------------------------!
!
   OPEN(1,FILE=TRIM(FNAME),STATUS='OLD')
   READ(1,*) N_ASSIM_TS                  !nomber of TS Obs station
   ALLOCATE(TS_OBS(N_ASSIM_TS))          !Type for TS_OBS

!
!--Read X,Y Coordinate of Measurement Stations---------------------------------!
!
   DO I=1,N_ASSIM_TS  
     READ(1,*)ITMP,TS_OBS(I)%X,TS_OBS(I)%Y,TS_OBS(I)%DEPTH,NLAY,TS_OBS(I)%SITA
     TS_OBS(I)%N_LAYERS = NLAY
     ALLOCATE(TS_OBS(I)%ODEPTH(NLAY))
     DO J=1,NLAY
       READ(1,*)TS_OBS(I)%ODEPTH(J)
       IF(TS_OBS(I)%ODEPTH(J) > TS_OBS(I)%DEPTH)THEN
         IF(MSR)WRITE(IPT,*)'OBSERVATION DEPTH',J,'OF TEMP/SALINITY OBS',I
         IF(MSR)WRITE(IPT,*)'EXCEEDS BATHYMETRIC DEPTH=',TS_OBS(I)%ODEPTH(J),TS_OBS(I)%DEPTH
         IF(MSR)WRITE(IPT,*)'HALTING...........'
         CALL PSTOP
       END IF
     END DO
   END DO
   TS_MAX_LAYER = MAXVAL(TS_OBS(1:N_ASSIM_TS)%N_LAYERS)

!
!--Shift Coordinates-----------------------------------------------------------!
!
# if !defined (SPHERICAL) 
   TS_OBS(:)%X = TS_OBS(:)%X - VXMIN 
   TS_OBS(:)%Y = TS_OBS(:)%Y - VYMIN 
# endif
   
!# if defined (SPHERICAL) 
   IF(SERIAL)THEN
     XCG = XC
     YCG = YC
     XG = VX
     YG = VY
   END IF
#  if defined(MULTIPROCESSOR)
   IF(PAR)CALL ACOLLECT(MYID,MSRID,NPROCS,EMAP,XC,XCG)
   IF(PAR)CALL ACOLLECT(MYID,MSRID,NPROCS,EMAP,YC,YCG)
   IF(PAR)CALL ACOLLECT(MYID,MSRID,NPROCS,NMAP,VX,XG)
   IF(PAR)CALL ACOLLECT(MYID,MSRID,NPROCS,NMAP,VY,YG)
   IF(PAR)CALL MPI_BCAST(XCG,NGL+1,MPI_F,0,MPI_COMM_WORLD,IERR)
   IF(PAR)CALL MPI_BCAST(YCG,NGL+1,MPI_F,0,MPI_COMM_WORLD,IERR)
   IF(PAR)CALL MPI_BCAST(XG,MGL+1,MPI_F,0,MPI_COMM_WORLD,IERR)
   IF(PAR)CALL MPI_BCAST(YG,MGL+1,MPI_F,0,MPI_COMM_WORLD,IERR)
#  endif
!# endif
!  
!
!--find the cell number (TS_OBS(:)%N_CELL) of Obs station---------------------
!
   RRTMP = 100000.0     !100km
!   IF(TS_NGASSIM)THEN
!       RRTMP = TS_NG_RADIUS
!   ELSE IF(TS_OIASSIM)THEN
!       RRTMP = TS_OI_RADIUS
!   ELSE
!       WRITE(IPT,*) 'IN SUBROUTINE SET_TS_ASSIM_DATA, ONE OF TS_NGASSIM AND &
!                         TS_OIASSIM MUST BE DEFINED AS TRUE.'
!      WRITE(IPT,*) 'HALTING.....'            
!       CALL PSTOP
!   END IF

   DO J= 1,N_ASSIM_TS 
     X0 = TS_OBS(J)%X
     Y0 = TS_OBS(J)%Y
     DO I=1,NGL
#    if !defined (SPHERICAL)
       Rtmp = SQRT((XCG(I)-X0)*(XCG(I)-X0)+(YCG(I)-Y0)*(YCG(I)-Y0))
#    else
       CALL ARC(X0,Y0,XCG(I),YCG(I),Rtmp)
#    endif
       IF(Rtmp < RRTMP) THEN
         X11_TMP = XG(NVG(I,1))-X0
         X22_TMP = XG(NVG(I,2))-X0
         X33_TMP = XG(NVG(I,3))-X0 
#        if defined (SPHERICAL)
         IF(X11_TMP >  180.0_SP)THEN
           X11_TMP = -360.0_SP+X11_TMP
         ELSEIF(X11_TMP < -180.0_SP)THEN
           X11_TMP =  360.0_SP+X11_TMP
         END IF
         IF(X22_TMP >  180.0_SP)THEN
           X22_TMP = -360.0_SP+X22_TMP
         ELSEIF(X22_TMP < -180.0_SP)THEN
           X22_TMP =  360.0_SP+X22_TMP
         END IF
         IF(X33_TMP >  180.0_SP)THEN
           X33_TMP = -360.0_SP+X33_TMP
         ELSEIF(X33_TMP < -180.0_SP)THEN
           X33_TMP =  360.0_SP+X33_TMP
         END IF
#        endif
         S11 = X22_TMP*(YG(NVG(I,3))-Y0)-X33_TMP*(YG(NVG(I,2))-Y0)
         S22 = X33_TMP*(YG(NVG(I,1))-Y0)-X11_TMP*(YG(NVG(I,3))-Y0)
         S33 = X11_TMP*(YG(NVG(I,2))-Y0)-X22_TMP*(YG(NVG(I,1))-Y0)
         IF(S11 <= 0._SP .AND. S22 <= 0._SP .AND. S33 <= 0._SP) THEN
           TS_OBS(J)%N_CELL = I
           GOTO 300
         ELSE
           TS_OBS(J)%N_CELL = 0
         ENDIF
       ELSE
         TS_OBS(J)%N_CELL = -1
       ENDIF
     ENDDO
300  CONTINUE
     IF(TS_OBS(J)%N_CELL <= 0) THEN
       IF(MSR) WRITE(IPT,*)'ERROR--T/S OBS SITE:',J,' OUT OF DOMAN',&
                            TS_OBS(J)%N_CELL 
       CALL PSTOP
     ENDIF

   ENDDO             
!--Gather AW0G,AWXG & AWYG use for interp grid to Obs station
   IF(.NOT. ALLOCATED(AW0G)) ALLOCATE(AW0G(0:NGL,3))
   IF(.NOT. ALLOCATED(AWXG)) ALLOCATE(AWXG(0:NGL,3))
   IF(.NOT. ALLOCATED(AWYG)) ALLOCATE(AWYG(0:NGL,3))
   
   IF(SERIAL)THEN
     AW0G = AW0
     AWXG = AWX
     AWYG = AWY
   END IF
#  if defined (MULTIPROCESSOR)
   IF(PAR)THEN     
     CALL ACOLLECT(MYID,MSRID,NPROCS,EMAP,AW0,AW0G)
     CALL ACOLLECT(MYID,MSRID,NPROCS,EMAP,AWX,AWXG)
     CALL ACOLLECT(MYID,MSRID,NPROCS,EMAP,AWY,AWYG)
     CALL MPI_BCAST(AW0G,(NGL+1)*3,MPI_F,0,MPI_COMM_WORLD,IERR)
     CALL MPI_BCAST(AWXG,(NGL+1)*3,MPI_F,0,MPI_COMM_WORLD,IERR)
     CALL MPI_BCAST(AWYG,(NGL+1)*3,MPI_F,0,MPI_COMM_WORLD,IERR)
   END IF
#  endif     
!
!--Close Current Observation Global File---------------------------------------!
!
   CLOSE(1)
       
   IF(TS_OIASSIM)THEN
!------------------------------------------------------------------------------!
!  Read Correlation Length of Scalar Observations                             !
!------------------------------------------------------------------------------!
       
!JQI   FNAME = "./"//TRIM(INPDIR)//"/"//trim(casename)//"_radius_ts.dat"

!JQI   INQUIRE(FILE=TRIM(FNAME),EXIST=FEXIST)
!JQI   IF(MSR .AND. .NOT.FEXIST)THEN
!JQI     WRITE(IPT,*)'CURRENT OBSERVATION FILE: ',FNAME,' DOES NOT EXIST'
!JQI     WRITE(IPT,*)'HALTING.....'
!JQI     CALL PSTOP
!JQI   END IF
     
!JQI   OPEN(1,FILE=TRIM(FNAME),STATUS='OLD')

   ALLOCATE(TS_PARAM(2,N_ASSIM_TS))

!JQI   DO I=1,N_ASSIM_TS
!JQI     READ(1,*)TS_PARAM(1,I),TS_PARAM(2,I)
!JQI   END DO

!JQI   CLOSE(1)
     TS_PARAM = TS_OI_RADIUS        !30000.0_SP
   END IF
   
!------------------------------------------------------------------------------!
!  Open Temp/Sal Observation Files for Each Observation Point and Read Data    !
!------------------------------------------------------------------------------!

!----Make Sure Temperature/Salinity Observation File Exists--------------------!
   IF(TS_NGASSIM .AND. .NOT. TS_OIASSIM)THEN
     ONAME = "./"//TRIM(INPUT_DIR)//"/"//trim(TS_NGASSIM_FILE)//".dat"
   ELSE IF(TS_OIASSIM .AND. .NOT. TS_NGASSIM)THEN
     ONAME = "./"//TRIM(INPUT_DIR)//"/"//trim(TS_OIASSIM_FILE)//".dat"
   ELSE
     CALL FATAL_ERROR("TS_NGASSIM and TS_OIASSIM cannot be both true or false")
   END IF  
   INQUIRE(FILE=TRIM(ONAME),EXIST=FEXIST)
   IF(MSR .AND. .NOT.FEXIST)THEN
     WRITE(IPT,*)'TEMP/SALINITY OBSERVATION FILE: ',ONAME,' DOES NOT EXIST'
     WRITE(IPT,*)'HALTING.....'
     CALL PSTOP
   END IF
!----Open Temp/Salinity Observation File for Read------------------------------!
   OPEN(1,FILE=ONAME,STATUS='old')  ; REWIND(1)

   DO I=1,N_ASSIM_TS
     READ(1,*,IOSTAT=IOS) nsite_tmp,NCNT
     IF(IOS<0) then
       WRITE(IPT,*) 'ERROR in read ',trim(casename),'_ts.dat at site number:',I
       CALL PSTOP
     ENDIF
     TS_OBS(I)%N_TIMES = NCNT
      
     IF(NCNT == 0)THEN
       IF(MSR)WRITE(IPT,*)'NO DATA FOR TS OBSERVATION',I
       CALL PSTOP
     END IF

!----Allocate Arrays to Hold Temp/Salinity (TEMP/SAL) and Time (TIME)----------!
     ALLOCATE(TS_OBS(I)%TIMES(TS_OBS(I)%N_TIMES))
     ALLOCATE(TS_OBS(I)%TEMP( TS_OBS(I)%N_TIMES , TS_OBS(I)%N_LAYERS ))
     ALLOCATE(TS_OBS(I)%SAL(  TS_OBS(I)%N_TIMES , TS_OBS(I)%N_LAYERS ))

!----Read in Current Data for Observation I------------------------------------!
     NLAY = TS_OBS(I)%N_LAYERS
     DO J=1,TS_OBS(I)%N_TIMES
       READ(1,*,IOSTAT=IOS)ODAYS,(TS_OBS(I)%TEMP(J,K),TS_OBS(I)%SAL(J,K),K=1,NLAY)
       IF(IOS < 0) THEN         ! ios=0 if all goes ok.
         WRITE(IPT,*)'ERROR in read ',trim(casename),'_ts.dat at site:',I,  &
                     'Time No:',J
         CALL PSTOP
       ENDIF  
       TS_OBS(I)%TIMES(J) = DAYS2TIME(ODAYS)
     END DO
     

!----Convert Time to Seconds---------------------------------------------------!
!----Shift Jan 1 Based Time Data to Dec 1 Based Time Data-----CASESPECIFIC-----!
!     IF(trim(CASENAME) == 'gom')THEN
!       TS_OBS(I)%TIMES = ((TS_OBS(I)%TIMES-1.0_SP)*24.0_SP+744.0_SP)*3600.0_SP
!     ELSE   
!       TS_OBS(I)%TIMES = TS_OBS(I)%TIMES*3600.0_SP*24.0_SP
!     END IF
 
!----Convert Temperature and Salinity to PSU/Celsius-(If Necessary)------------!
     TS_OBS(I)%TEMP = TS_OBS(I)%TEMP 
     TS_OBS(I)%SAL  = TS_OBS(I)%SAL  
  END DO
  CLOSE(1)

!------------------------------------------------------------------------------!
!  Count Number of Points with Bad Data (TEMP = 0. .OR. SAL = 0.)         
!------------------------------------------------------------------------------!
  NBD_CNT = 0
  DO I=1,N_ASSIM_TS
    DO J=1,TS_OBS(I)%N_TIMES
      DO K=1,TS_OBS(I)%N_LAYERS
        IF(TS_OBS(I)%TEMP(J,K) < -90. .OR. TS_OBS(I)%SAL(J,K) < -90.) THEN
          NBD_CNT = NBD_CNT + 1
        END IF
      END DO
    END DO
  END DO
 
!------------------------------------------------------------------------------!
!  Compute Spatial Interpolation Weights for each Mooring Location 
!------------------------------------------------------------------------------!
!   LMIN = 100000000.
!   DO I=1,N_ASSIM_TS
!     X0 = TS_OBS(I)%X
!     Y0 = TS_OBS(I)%Y
!     DO J=1,MGL
!       DX = ABS(XG(J)-X0)
!       DY = ABS(YG(J)-Y0)
!       IF(SQRT(DX**2 + DY**2) < LMIN)THEN
!         LMIN = SQRT(DX**2 + DY**2)
!         JMIN = J
!       END IF
!      END DO
!      TS_OBS(I)%SITA = SITA_GD(JMIN) + 3.14159_SP/2.0_SP
!    end do


   ALLOCATE(ITEMP(MGL),FTEMP(MGL),DA_TS(MGL))     ; DA_TS = 0

   DO I=1,N_ASSIM_TS
     X0 = TS_OBS(I)%X
     Y0 = TS_OBS(I)%Y
     ISOBATH_ANGLE = TS_OBS(I)%SITA*DEG2RAD        !/180.0_SP*3.1415926_SP
     ECNT = 0
     DO J=1,MGL
#      if defined (SPHERICAL)
       CALL ARCX(X0,Y0,XG(J),YG(J),DX)
       CALL ARCY(X0,Y0,XG(J),YG(J),DY)
       CALL ARC(X0,Y0,XG(J),YG(J),RD)
#      else
       DX = ABS(XG(J)-X0)
       DY = ABS(YG(J)-Y0)
       RD = SQRT(DX**2 + DY**2)
#      endif

       IF(TS_NGASSIM)THEN
         TS_RADIUS = TS_NG_RADIUS
       ELSE IF(TS_OIASSIM)THEN
         TS_RADIUS = TS_OI_RADIUS
       ELSE
!JQI         WRITE(IPT,*) 'IN SUBROUTINE SET_TS_ASSIM_DATA, ONE OF TS_NGASSIM AND ',&
!JQI	              'TS_OIASSIM MUST BE DEFINED AS TRUE.'
         WRITE(IPT,*) 'IN SUBROUTINE SET_TS_ASSIM_DATA, ONE OF TS_NGASSIM AND TS_OIASSIM MUST BE DEFINED AS TRUE.'
         WRITE(IPT,*) 'HALTING.....'       
         CALL PSTOP
       END IF 
         
       IF(RD <= TS_RADIUS)THEN
         DA_TS(J)   = 1
         ECNT        = ECNT + 1      
         ITEMP(ECNT) =  J
         FTEMP(ECNT) = (TS_RADIUS**2 - RD**2) / (TS_RADIUS**2 + RD**2)
         IF(ISOBATH_ANGLE==0)THEN
           DIR_WEIGHT = 1.0_SP
         ELSE
           ANG_OBS_SIM = ATAN2(DY,DX)
           D_ANGLE     = ANG_OBS_SIM - ISOBATH_ANGLE
           D_ANGLE     = D_ANGLE - INT(D_ANGLE/PI)*PI      !3.1415926_SP)*3.1415926_SP
           D_ANGLE     = ABS(D_ANGLE)
!           DIR_WEIGHT  = (ABS(D_ANGLE-0.5*3.1415926_SP)+ALF*3.1415926_SP)/ &
!                         ((0.5_SP+ALF)*3.1415926_SP)
           DIR_WEIGHT  = (ABS(D_ANGLE-0.5*PI)+ALF*PI)/((0.5_SP+ALF)*PI)
!           IF(J==TS_OBS(J)%N_CELL)DIR_WEIGHT = 1.0_SP
         END IF
         FTEMP(ECNT) = FTEMP(ECNT)*DIR_WEIGHT
       END IF
     END DO
     IF(ECNT == 0)THEN
       WRITE(IPT,*)'ERROR SETTING UP TEMP/SAL DATA ASSIMILATION'
       WRITE(IPT,*)'NO NODES LIE WITHIN RADIUS',TS_RADIUS
       WRITE(IPT,*)'OF OBSERVATION POINT',I
       CALL PSTOP   
     ELSE
       TS_OBS(I)%N_INTPTS = ECNT
       ALLOCATE(TS_OBS(I)%INTPTS(ECNT))
       ALLOCATE(TS_OBS(I)%X_WEIGHT(ECNT))
       TS_OBS(I)%INTPTS(1:ECNT)  = ITEMP(1:ECNT)
       TS_OBS(I)%X_WEIGHT(1:ECNT) = FTEMP(1:ECNT)
     END IF
   END DO
   DEALLOCATE(FTEMP,ITEMP)

!------------------------------------------------------------------------------!
!  Compute Sigma Layer Weights for Vertical Interpolation            
!------------------------------------------------------------------------------!
   ALLOCATE(ZZ_G(0:MGL,KB))  
   IF(.NOT. ALLOCATED(HG)) ALLOCATE(HG(0:MGL))
   IF(SERIAL)THEN
     ZZ_G = ZZ
     HG = H
   END IF      
#  if defined (MULTIPROCESSOR)   
   IF(PAR)CALL ACOLLECT(MYID,MSRID,NPROCS,NMAP,ZZ,ZZ_G)
   IF(PAR)CALL ACOLLECT(MYID,MSRID,NPROCS,NMAP,H,HG) 
   IF(PAR)CALL MPI_BCAST(ZZ_G,(MGL+1)*KB,MPI_F,0,MPI_COMM_WORLD,IERR)
   IF(PAR)CALL MPI_BCAST(HG,MGL+1,MPI_F,0,MPI_COMM_WORLD,IERR)
#  endif

   DO I=1,N_ASSIM_TS
     NLAY = TS_OBS(I)%N_LAYERS
     ALLOCATE(TS_OBS(I)%S_INT(NLAY,2))
     ALLOCATE(TS_OBS(I)%S_WEIGHT(NLAY,2))
     
     X0 = TS_OBS(I)%X
     Y0 = TS_OBS(I)%Y
# if defined (SPHERICAL)
     DO NN=1,NGL
       CALL ARC(X0,Y0,XCG(NN),YCG(NN),RDLIST(NN,1))
     END DO
# else
     RDLIST(1:NGL,1) = SQRT((XCG(1:NGL)-X0)**2 + (YCG(1:NGL)-Y0)**2)
# endif
     RDLAST = -1.0_SP
in:  DO WHILE(.TRUE.)
       LOCIJ = MINLOC(RDLIST,RDLIST>RDLAST)
       MIN_LOC = LOCIJ(1)
       IF(MIN_LOC == 0)THEN
         EXIT in
       END IF
       XTRI = XG(NVG(MIN_LOC,1:3))
       YTRI = YG(NVG(MIN_LOC,1:3))
       RDLAST = RDLIST(MIN_LOC,1)
       IF(ISINTRIANGLE1(XTRI,YTRI,X0,Y0))THEN
         JJ = MIN_LOC
         EXIT IN
       END IF
       RDLAST = RDLIST(MIN_LOC,1)
     END DO IN      

     ND1=NVG(JJ,1)
     ND2=NVG(JJ,2)
     ND3=NVG(JJ,3)

# if defined (SPHERICAL)
     CALL ARCX(X0,Y0,XCG(JJ),YCG(JJ),X0C)
     CALL ARCX(X0,Y0,XCG(JJ),YCG(JJ),Y0C)
# else
     X0C = X0 - XCG(JJ)
     Y0C = Y0 - YCG(JJ)
# endif
     !----Linear Interpolation of Bathymetry------------------------------------------!
     H0 = AW0G(JJ,1)*HG(ND1)+AW0G(JJ,2)*HG(ND2)+AW0G(JJ,3)*HG(ND3)
     HX = AWXG(JJ,1)*HG(ND1)+AWXG(JJ,2)*HG(ND2)+AWXG(JJ,3)*HG(ND3)
     HY = AWYG(JJ,1)*HG(ND1)+AWYG(JJ,2)*HG(ND2)+AWYG(JJ,3)*HG(ND3)
     TEMP_H = H0 + HX*X0C + HY*Y0C !here is the interpolated observation depth

     DO K=1,KBM1
      !----Linear Interpolation of Bathymetry------------------------------------------!
      H0 = AW0G(JJ,1)*ZZ_G(ND1,K)*HG(ND1)+AW0G(JJ,2)*ZZ_G(ND2,K)*HG(ND2)+AW0G(JJ,3)*ZZ_G(ND3,K)*HG(ND3)
      HX = AWXG(JJ,1)*ZZ_G(ND1,K)*HG(ND1)+AWXG(JJ,2)*ZZ_G(ND2,K)*HG(ND2)+AWXG(JJ,3)*ZZ_G(ND3,K)*HG(ND3)
      HY = AWYG(JJ,1)*ZZ_G(ND1,K)*HG(ND1)+AWYG(JJ,2)*ZZ_G(ND2,K)*HG(ND2)+AWYG(JJ,3)*ZZ_G(ND3,K)*HG(ND3)
      ZZ_OB(K) = (H0 + HX*X0C + HY*Y0C)/TEMP_H
     END DO  

!     DELTA=(XG(ND2)-XG(ND1))*(YG(ND3)-YG(ND1))-     &
!           (XG(ND3)-XG(ND1))*(YG(ND2)-YG(ND1))

!     IF(SERIAL)THEN
!     DO K=1,KBM1
!       COFA=(YG(ND3)-YG(ND1))*(ZZ(ND2,K)-ZZ(ND1,K))-   &
!            (YG(ND2)-YG(ND1))*(ZZ(ND3,K)-ZZ(ND1,K))
!       COFB=(XG(ND2)-XG(ND1))*(ZZ(ND3,K)-ZZ(ND1,K))-   &
!            (XG(ND3)-XG(ND1))*(ZZ(ND2,K)-ZZ(ND1,K))
!       COFA=COFA/DELTA
!       COFB=COFB/DELTA
!       COFC=ZZ(ND1,K)-COFA*XG(ND1)-COFB*YG(ND1)
!       ZZ_OB(K)=COFA*X0+COFB*Y0+COFC
!     END DO  
!     END IF

!#  if defined (MULTIPROCESSOR)   
!     IF(PAR)THEN
!     DO K=1,KBM1
!       COFA=(YG(ND3)-YG(ND1))*(ZZ_G(ND2,K)-ZZ_G(ND1,K))-   &
!            (YG(ND2)-YG(ND1))*(ZZ_G(ND3,K)-ZZ_G(ND1,K))
!       COFB=(XG(ND2)-XG(ND1))*(ZZ_G(ND3,K)-ZZ_G(ND1,K))-   &
!            (XG(ND3)-XG(ND1))*(ZZ_G(ND2,K)-ZZ_G(ND1,K))
!       COFA=COFA/DELTA
!       COFB=COFB/DELTA
!       COFC=ZZ_G(ND1,K)-COFA*XG(ND1)-COFB*YG(ND1)
!       ZZ_OB(K)=COFA*X0+COFB*Y0+COFC
!     END DO  
!     END IF
!#  endif
     
     DO J=1,NLAY
       SIGMA_C = -TS_OBS(I)%ODEPTH(J)/TEMP_H         !TS_OBS(I)%DEPTH
       DO K=2,KBM1
         IF(ZZ_OB(K) <= SIGMA_C .AND. ZZ_OB(K-1) > SIGMA_C)THEN 
           TS_OBS(I)%S_INT(J,1) = K-1
           TS_OBS(I)%S_INT(J,2) = K
           TS_OBS(I)%S_WEIGHT(J,1) = (SIGMA_C-ZZ_OB(K))/(ZZ_OB(K-1)-ZZ_OB(K))
           TS_OBS(I)%S_WEIGHT(J,2) = 1.0_SP - TS_OBS(I)%S_WEIGHT(J,1) 
         END IF  
       END DO
       IF(ZZ_OB(1) <= SIGMA_C)THEN  !!OBSERVATION ABOVE CENTROID OF FIRST SIGMA LAYER
         TS_OBS(I)%S_INT(J,1) = 1
         TS_OBS(I)%S_INT(J,2) = 1
         TS_OBS(I)%S_WEIGHT(J,1) = 1.0_SP
         TS_OBS(I)%S_WEIGHT(J,2) = 0.0_SP
       END IF
       IF(ZZ_OB(KBM1) > SIGMA_C)THEN !!OBSERVATION BELOW CENTROID OF BOTTOM SIGMA LAYER
         TS_OBS(I)%S_INT(J,1) = KBM1
         TS_OBS(I)%S_INT(J,2) = KBM1
         TS_OBS(I)%S_WEIGHT(J,1) = 1.0_SP
         TS_OBS(I)%S_WEIGHT(J,2) = 0.0_SP
       END IF

     END DO
   END DO

   DEALLOCATE(ZZ_G)  

!------------------------------------------------------------------------------!
!  Report Number of Interpolation Points, Location and Number of Data 
!------------------------------------------------------------------------------!
   IF(.NOT. MSR)RETURN

   WRITE(IPT,*)
   WRITE(IPT,*)'!            TEMP/SALINITY OBSERVATION DATA           '
#  if defined (SPHERICAL)
   WRITE(IPT,*)" MOORING#    LON       LAT    #INTERP PTS  #DATA TIMES  NEAR_NODE  SITA"
#  else
   WRITE(IPT,*)" MOORING#   X(KM)      Y(KM)  #INTERP PTS  #DATA TIMES  NEAR_NODE  SITA"
#  endif
   DO I=1,N_ASSIM_TS
     MAXEL = MAXLOC(TS_OBS(I)%X_WEIGHT,DIM=1)
     WRITE(IPT,'(2X,I5,3X,F8.1,3X,F8.1,3X,I6,5X,I6,5X,I6,5X,F8.1)') &
#  if defined (SPHERICAL)
     I,TS_OBS(I)%X,TS_OBS(I)%Y, &
#  else
     I,(TS_OBS(I)%X+VXMIN)/1000.,(TS_OBS(I)%Y+VYMIN)/1000., &
#  endif
       TS_OBS(I)%N_INTPTS,TS_OBS(I)%N_TIMES,TS_OBS(I)%INTPTS(MAXEL),&
       TS_OBS(I)%SITA
   END DO
   WRITE(IPT,*)
   WRITE(IPT,*)'NUMBER OF BAD TS DATA POINTS: ',NBD_CNT
   WRITE(IPT,*)" MOORING #   BEGIN TIME  END TIME"
   DO I=1,N_ASSIM_TS
!   WRITE(IPT,*)I,TS_OBS(I)%TIMES(1)/(24.*3600.),&
!       TS_OBS(I)%TIMES(TS_OBS(I)%N_TIMES)/(24.*3600.)
    WRITE(IPT,*)I,TS_OBS(I)%TIMES(1)%MJD,&
        TS_OBS(I)%TIMES(TS_OBS(I)%N_TIMES)%MJD
   END DO
   IF (DBG_SET(DBG_SBR)) WRITE(IPT,*) "END: SET_TS_ASSIM_DATA"

   RETURN
   END SUBROUTINE SET_TS_ASSIM_DATA

   SUBROUTINE SET_TS_ASSIM_DATA_WINDOW

!------------------------------------------------------------------------------!
!  SET UP ASSIMILATION DATA FOR TEMP/SAL OBSERVATIONS in ASSIMILATION WINDOW   |
!------------------------------------------------------------------------------!
   IMPLICIT NONE
   INTEGER I,J,IW,NLAY,NTIME
   REAL(SP), ALLOCATABLE, DIMENSION(:) :: FTEMP
   REAL(SP) :: T_THRESH,DT_MIN

!   INTEGER, ALLOCATABLE, DIMENSION(:)   :: ID_STATION

   IF(DBG_SET(DBG_SBR)) write(ipt,*) "START: SET_TS_ASSIM_DATA_WINDOW"

   IF(ALLOCATED(ID_STATION)) DEALLOCATE(ID_STATION)
   ALLOCATE(ID_STATION(N_ASSIM_TS))
!   T_THRESH = TS_OI_ASTIME_WINDOW+12.0_SP*3600.0_SP   
   T_THRESH = TS_OI_ASTIME_WINDOW
   IF(MSR) write(555,*) '-------------test---1--',T_THRESH
   IW = 0
   N_ASSIM_TS_W = 0
   DO I=1,N_ASSIM_TS
     NLAY  = TS_OBS(I)%N_LAYERS
     NTIME = TS_OBS(I)%N_TIMES
     ALLOCATE(FTEMP(NTIME))
     DO J=1,NTIME
       FTEMP(J) = SECONDS(IntTime - TS_OBS(I)%TIMES(J))
       FTEMP(J) = ABS(FTEMP(J))
     END DO
     DT_MIN = MINVAL(FTEMP(1:NTIME))
     IF(DT_MIN < T_THRESH) THEN
       IW = IW+1
       ID_STATION(IW)=I
     ENDIF
     DEALLOCATE(FTEMP)
   ENDDO
   IF(MSR) write(555,*) '-------------test---2-IW=',IW
   N_ASSIM_TS_W = IW
   RETURN
   END SUBROUTINE SET_TS_ASSIM_DATA_WINDOW

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|

   SUBROUTINE TEMP_ASSIMILATION
!==============================================================================|
!  USE TEMP OBSERVATION DATA TO ADJUST TEMP FIELD                              |
!==============================================================================|
#  if defined (PWP)
   USE MOD_PWP, only : MLD
#  endif
   IMPLICIT NONE
   REAL(SP), ALLOCATABLE, DIMENSION(:,:) :: TINT,TCORR,TG,TCORR1
   REAL(SP), ALLOCATABLE, DIMENSION(:,:) :: TWGHT_T
   REAL(SP), ALLOCATABLE, DIMENSION(:)   :: FTEMP
   REAL(SP) :: WEIGHT,DEFECT,CORRECTION,DT_MIN,SIMTIME,T_THRESH,WGHT,TOT_WGHT
   REAL(SP) :: U1,U2,V1,V2,W1,W2,WEIGHT1,WEIGHT2
   TYPE(TIME) ::DIFTIME
   INTEGER I,J,K,J1,K1,K2,NLAY,ITIMEA,NTIME,IERR
   INTRINSIC MINLOC
   INTEGER ID
   INTEGER IP,JN1,JN2,JN3
   REAL(SP) XSC,YSC,COFT0,COFTX,COFTY

#  if defined (PWP)
   INTEGER, ALLOCATABLE, DIMENSION(:) :: MLD_GLB
#  endif
!==============================================================================|

   IF (DBG_SET(DBG_SBR)) WRITE(IPT,*) "START: TEMP_ASSIMILATION" 

   CALL  SET_TS_ASSIM_DATA_WINDOW
   IF(MSR) then
     DO I=1, N_ASSIM_TS_W
       ID = ID_STATION(I)
     ENDDO
   ENDIF
   IF(N_ASSIM_TS_W.le.0) return
       
!------------------------------------------------------------------------------!
!  Gather T Fields to Master Processor                                   ! 
!------------------------------------------------------------------------------!
   ALLOCATE(TG(0:MGL,KB))        ;TG = 0.0_SP
#  if defined (MULTIPROCESSOR)
   IF(PAR)THEN
!     CALL GATHER(LBOUND(TF1,1),UBOUND(TF1,1),M,MGL,KB,MYID,NPROCS,NMAP,TF1,TG)
     CALL ACOLLECT(MYID,MSRID,NPROCS,NMAP,TF1,TG)
   END IF
#  endif
   IF(SERIAL)THEN
     TG(1:MGL,1:KBM1) = TF1(1:MGL,1:KBM1)
   END IF

#  if defined (PWP)
   ALLOCATE(MLD_GLB(0:MGL))
#  if defined (MULTIPROCESSOR)
   IF(PAR) CALL ACOLLECT(MYID,MSRID,NPROCS,NMAP,MLD,MLD_GLB)
#  endif
   IF(SERIAL) MLD_GLB(1:MGL) = MLD(1:MGL)
#  endif
!------------------------------------------------------------------------------!
!  Calculate Temporal Weight of Measurement (I) at Time(TIME)                  ! 
!------------------------------------------------------------------------------!

   IF(MSR)THEN
   TS_OBS%T_WEIGHT = 0. 
   IF(TS_NGASSIM)THEN
     T_THRESH = TS_NG_ASTIME_WINDOW   
   ELSE IF(TS_OIASSIM)THEN  
     T_THRESH = TS_OI_ASTIME_WINDOW   
   ELSE
     WRITE(IPT,*) "EITHER TS_NGASSIM OR TS_OIASSIM MUST BE SET TRUE"
     CALL PSTOP
   END IF       
!   SIMTIME          = TIME*86400
   SIMTIME = DTI*FLOAT(IINT)

   DO I=1,N_ASSIM_TS_w
     ID = ID_STATION(I)
     NTIME = TS_OBS(Id)%N_TIMES
     ALLOCATE(FTEMP(NTIME)) 
!     FTEMP(1:NTIME) = ABS(SIMTIME - TS_OBS(I)%TIMES(1:NTIME))
     DO j=1,NTIME
        FTEMP(J) = SECONDS(IntTime - TS_OBS(Id)%TIMES(J))
        FTEMP(J) = ABS(FTEMP(J))
     ENDDO
     DT_MIN = MINVAL(FTEMP(1:NTIME))
     TS_OBS(Id)%N_T_WEIGHT = MINLOC(FTEMP,DIM=1)

     IF(DT_MIN < T_THRESH)THEN     
       IF(DT_MIN < .5_SP*T_THRESH) THEN
         TS_OBS(Id)%T_WEIGHT = 1.0_SP
       ELSE
         TS_OBS(Id)%T_WEIGHT = (T_THRESH-DT_MIN)/T_THRESH*2.0_SP
       END IF
     END IF

     DEALLOCATE(FTEMP)
   END DO
   
       
!------------------------------------------------------------------------------!
!  Interpolate Simulation Data to Local Observation Point                      ! 
!------------------------------------------------------------------------------!
       
   ALLOCATE(TINT(N_ASSIM_TS,TS_MAX_LAYER)) ; TINT = 0. 

!   IF(TS_NGASSIM)THEN
   DO I=1,N_ASSIM_TS_w
       ID = ID_STATION(I) 
     DO J=1,TS_OBS(Id)%N_INTPTS
       J1        = TS_OBS(Id)%INTPTS(J)
       WGHT      = TS_OBS(Id)%X_WEIGHT(J)
       NLAY      = TS_OBS(Id)%N_LAYERS
       DO K=1,NLAY
         U1 = TG(J1,TS_OBS(Id)%S_INT(K,1))
         U2 = TG(J1,TS_OBS(Id)%S_INT(K,2))
         W1 = TS_OBS(Id)%S_WEIGHT(K,1)
         W2 = TS_OBS(Id)%S_WEIGHT(K,2)
         TINT(I,K) = TINT(I,K) + (U1*W1 + U2*W2)*WGHT 
       END DO
     END DO
     TOT_WGHT = SUM(TS_OBS(Id)%X_WEIGHT(1:TS_OBS(Id)%N_INTPTS))
     TINT(I,1:NLAY) = TINT(I,1:NLAY)/TOT_WGHT
   END DO

!   ELSE IF(TS_OIASSIM)THEN
!   DO I=1,N_ASSIM_TS
!     IP = TS_OBS(I)%N_CELL
!     NLAY      = TS_OBS(I)%N_LAYERS
!       XSC = TS_OBS(I)%X-XCG(IP)
!       YSC = TS_OBS(I)%Y-YCG(IP)
!       JN1 = NVG(IP,1)
!       JN2 = NVG(IP,2)
!       JN3 = NVG(IP,3)
      
!       DO K=1,NLAY
!         K1 = TS_OBS(I)%S_INT(K,1)
!         K2 = TS_OBS(I)%S_INT(K,2)
!         COFT0 = AW0G(IP,1)*TG(JN1,K1)+AW0G(IP,2)*TG(JN2,K1)+AW0G(IP,3)*TG(JN3,K1)
!         COFTX = AWXG(IP,1)*TG(JN1,K1)+AWXG(IP,2)*TG(JN2,K1)+AWXG(IP,3)*TG(JN3,K1)
!         COFTY = AWYG(IP,1)*TG(JN1,K1)+AWYG(IP,2)*TG(JN2,K1)+AWYG(IP,3)*TG(JN3,K1)
!         U1 = COFT0 + COFTX*XSC + COFTY*YSC
!         COFT0 = AW0G(IP,1)*TG(JN1,K2)+AW0G(IP,2)*TG(JN2,K2)+AW0G(IP,3)*TG(JN3,K2)
!         COFTX = AWXG(IP,1)*TG(JN1,K2)+AWXG(IP,2)*TG(JN2,K2)+AWXG(IP,3)*TG(JN3,K2)
!         COFTY = AWYG(IP,1)*TG(JN1,K2)+AWYG(IP,2)*TG(JN2,K2)+AWYG(IP,3)*TG(JN3,K2)
!         U2 = COFT0 + COFTX*XSC + COFTY*YSC
!         W1 = TS_OBS(I)%S_WEIGHT(K,1)
!         W2 = TS_OBS(I)%S_WEIGHT(K,2)
!	 TINT(I,K) = (U1*W1 + U2*W2)
!       ENDDO
!     ENDDO   	 
!   ELSE
!     WRITE(IPT,*)'EITHER TS_NGASSIM OR TS_OIASSIM SHOULD BE TRUE'
!     CALL PSTOP
!  END IF  
!------------------------------------------------------------------------------!
!  Compute Local Correction by Interpolating Observed/Computed Defect          ! 
!------------------------------------------------------------------------------!

   ALLOCATE(TWGHT_T(MGL,KBM1))   ; TWGHT_T = 0.
   ALLOCATE(TCORR(MGL,KBM1))     ; TCORR   = 0.

   IF(TS_NGASSIM)THEN
     DO I=1,N_ASSIM_TS_w
     ID = ID_STATION(I) 
       DO J=1,TS_OBS(Id)%N_INTPTS
         J1     = TS_OBS(Id)%INTPTS(J)
         ITIMEA  = TS_OBS(Id)%N_T_WEIGHT
         NLAY   = TS_OBS(Id)%N_LAYERS
         DO K=1,NLAY
           K1             = TS_OBS(Id)%S_INT(K,1)
           K2             = TS_OBS(Id)%S_INT(K,2)
           W1             = TS_OBS(Id)%S_WEIGHT(K,1)
           W2             = TS_OBS(Id)%S_WEIGHT(K,2)

           DEFECT         = TS_OBS(Id)%TEMP(ITIMEA,K) - TINT(I,K)
           IF(ABS(DEFECT) < 20.0_SP)THEN     !quality control
             WEIGHT1        = TS_OBS(Id)%T_WEIGHT*TS_OBS(Id)%X_WEIGHT(J)*W1
             WEIGHT2        = TS_OBS(Id)%T_WEIGHT*TS_OBS(Id)%X_WEIGHT(J)*W2
             TWGHT_T(J1,K1) = TWGHT_T(J1,K1) + WEIGHT1   
             TWGHT_T(J1,K2) = TWGHT_T(J1,K2) + WEIGHT2   

             CORRECTION     = TS_GAMA*DEFECT
!           TCORR(J1,K1)   = TCORR(J1,K1) + CORRECTION*WEIGHT1
!           TCORR(J1,K2)   = TCORR(J1,K2) + CORRECTION*WEIGHT2
             TCORR(J1,K1)   = TCORR(J1,K1) + CORRECTION*WEIGHT1*WEIGHT1
             TCORR(J1,K2)   = TCORR(J1,K2) + CORRECTION*WEIGHT2*WEIGHT2
           ENDIF
         END DO
       END DO
     END DO

!------------------------------------------------------------------------------!
!  Nudge Simulation Data Using Local Corrections                               ! 
!------------------------------------------------------------------------------!

     DO I=1,MGL

       IF(SSTGRD_ASSIM) THEN
         K1 = 2
       ELSE
         K1 = 1
       ENDIF
#      if defined (PWP)
       IF(MLD_GLB(I)<0) THEN
         K1 = 2
       ELSE
         K1 = MLD_GLB(I)+1 
       ENDIF
#      endif

       DO K=K1,KBM1
         IF(DA_TS(I) == 1 .AND. TWGHT_T(I,K) > 1.0E-08)THEN
           TG(I,K) = TG(I,K) + DTI*TS_GALPHA*TCORR(I,K)/TWGHT_T(I,K)
         END IF
       END DO
     END DO

     DEALLOCATE(TWGHT_T,TCORR,TINT)
   ELSE IF(TS_OIASSIM)THEN
     DO I=1,N_ASSIM_TS_w 
     ID = ID_STATION(I)

       ITIMEA  = TS_OBS(Id)%N_T_WEIGHT
       NLAY   = TS_OBS(Id)%N_LAYERS
       DO K=1,NLAY
         K1             = TS_OBS(Id)%S_INT(K,1)
         K2             = TS_OBS(Id)%S_INT(K,2)
         W1             = TS_OBS(Id)%S_WEIGHT(K,1)
         W2             = TS_OBS(Id)%S_WEIGHT(K,2)

         DEFECT         = TS_OBS(Id)%TEMP(ITIMEA,K) - TINT(I,K)
         IF(ABS(DEFECT) < 20.0)THEN     !quality control
           WEIGHT1        = TS_OBS(Id)%T_WEIGHT*W1
           WEIGHT2        = TS_OBS(Id)%T_WEIGHT*W2
           TWGHT_T(I,K1) = TWGHT_T(I,K1) + WEIGHT1   
           TWGHT_T(I,K2) = TWGHT_T(I,K2) + WEIGHT2   
           CORRECTION     = DEFECT
           TCORR(I,K1)   = TCORR(I,K1) + CORRECTION*WEIGHT1
           TCORR(I,K2)   = TCORR(I,K2) + CORRECTION*WEIGHT2
         ENDIF
       END DO
     END DO

     DO I=1,N_ASSIM_TS_w
       DO K=1,KBM1
         IF(TWGHT_T(I,K) > 1.0E-8)THEN
           TCORR(I,K)=TCORR(I,K)/TWGHT_T(I,K)
         END IF
       END DO
     END DO
        
!------------------------------------------------------------------------------!
!  'OI' Simulation Data Using Local Corrections                               ! 
!------------------------------------------------------------------------------!

     ALLOCATE(TCORR1(MGL,KBM1));TCORR1=0.0_SP

     DO K=1,KBM1
       CALL TS_OPTIMINTERP(TCORR(:,K),TCORR1(:,K))
     END DO

     DO I=1,MGL

       IF(SSTGRD_ASSIM) THEN
         K1 = 2
       ELSE 
         K1 = 1
       ENDIF
#      if defined (PWP) 
       IF(MLD_GLB(I)<0) THEN
         K1 = 2
       ELSE 
         K1 = MLD_GLB(I)+1
       ENDIF
#      endif

       DO K=K1,KBM1
         IF(DA_TS(I) == 1)THEN
!QXU         TG(I,K) = TG(I,K) + TCORR1(I,K)
!           TG(I,K) = TG(I,K) + TS_GALPHA*TCORR1(I,K)
           TG(I,K) = TG(I,K) + TS_OIGALPHA*TCORR1(I,K)
         END IF
       END DO
     END DO

     DEALLOCATE(TWGHT_T,TCORR,TINT,TCORR1)
   ELSE
     WRITE(IPT,*)'EITHER TS_NGASSIM OR TS_OIASSIM SHOULD BE TRUE'
     CALL PSTOP
   END IF

   END IF  !!MASTER


!------------------------------------------------------------------------------!
!  Disperse New Data Fields to Slave Processors                                ! 
!------------------------------------------------------------------------------!
   IF(SERIAL)THEN
     TF1(1:M,1:KBM1) = TG(1:M,1:KBM1)
   END IF
#  if defined (MULTIPROCESSOR) 
   CALL MPI_BCAST(TG,(MGL+1)*KB,MPI_F,0,MPI_COMM_WORLD,IERR)
   IF(PAR)THEN
     DO I=1,M
       TF1(I,1:KBM1) = TG(NGID(I),1:KBM1)
     END DO
   END IF
#  endif

   DEALLOCATE(TG)
#  if defined (PWP)
   DEALLOCATE(MLD_GLB)
#  endif
          
   IF (DBG_SET(DBG_SBR)) WRITE(IPT,*) "END: TEMP_ASSIMILATION" 

   RETURN
   END SUBROUTINE TEMP_ASSIMILATION
!==============================================================================|
!==============================================================================|

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|

   SUBROUTINE SALT_ASSIMILATION
!==============================================================================|
!  USE SALT OBSERVATION DATA TO ADJUST SALINITY FIELD                          |
!==============================================================================|
#  if defined (PWP)
   USE MOD_PWP, only : MLD
#  endif
   IMPLICIT NONE
   REAL(SP), ALLOCATABLE, DIMENSION(:,:) :: SINT,SCORR,SG,SCORR1
   REAL(SP), ALLOCATABLE, DIMENSION(:,:) :: TWGHT_S
   REAL(SP), ALLOCATABLE, DIMENSION(:)   :: FTEMP
   REAL(SP) :: WEIGHT,DEFECT,CORRECTION,DT_MIN,SIMTIME,T_THRESH,WGHT,TOT_WGHT
   REAL(SP) :: U1,U2,V1,V2,W1,W2,WEIGHT1,WEIGHT2
   TYPE(TIME) ::DIFTIME
   INTEGER I,J,K,J1,K1,K2,NLAY,ITIMEA,NTIME,IERR
   INTRINSIC MINLOC

   INTEGER IP,JN1,JN2,JN3
   REAL(SP) XSC,YSC,COFT0,COFTX,COFTY
   INTEGER ID
#  if defined (PWP)
   INTEGER, ALLOCATABLE, DIMENSION(:) :: MLD_GLB
#  endif
!==============================================================================|

   IF (DBG_SET(DBG_SBR)) WRITE(IPT,*) "START: SALT_ASSIMILATION" 
       IF(N_ASSIM_TS_W.le.0) return     
!------------------------------------------------------------------------------!
!  Gather S Fields to Master Processor                                         ! 
!------------------------------------------------------------------------------!
   ALLOCATE(SG(0:MGL,KB))       ;SG = 0.0_SP
#  if defined (MULTIPROCESSOR)
   IF(PAR)THEN
!     CALL GATHER(LBOUND(SF1,1),UBOUND(SF1,1),M,MGL,KB,MYID,NPROCS,NMAP,SF1,SG)
     CALL ACOLLECT(MYID,MSRID,NPROCS,NMAP,SF1,SG)
   END IF
#  endif
   IF(SERIAL)THEN
     SG(1:MGL,1:KBM1) = SF1(1:MGL,1:KBM1)
   END IF

#  if defined (PWP)
   ALLOCATE(MLD_GLB(0:MGL))
   IF(PAR) CALL ACOLLECT(MYID,MSRID,NPROCS,NMAP,MLD,MLD_GLB)
   IF(SERIAL) MLD_GLB(1:MGL) = MLD(1:MGL)
#  endif
!------------------------------------------------------------------------------!
!  Calculate Temporal Weight of Measurement (I) at Time(TIME)                  ! 
!------------------------------------------------------------------------------!

   IF(MSR)THEN
   TS_OBS%T_WEIGHT = 0. 
   IF(TS_NGASSIM)THEN
     T_THRESH = TS_NG_ASTIME_WINDOW   
   ELSE IF(TS_OIASSIM)THEN  
     T_THRESH = TS_OI_ASTIME_WINDOW   
   ELSE
     WRITE(IPT,*) "EITHER TS_NGASSIM OR TS_OIASSIM MUST BE SET TRUE"
     CALL PSTOP
   END IF       
!   SIMTIME          = TIME*86400
   SIMTIME = DTI*FLOAT(IINT)

   DO I=1,N_ASSIM_TS_w
     ID = ID_STATION(I)
     NTIME = TS_OBS(Id)%N_TIMES
     ALLOCATE(FTEMP(NTIME)) 
!     FTEMP(1:NTIME) = ABS(SIMTIME - TS_OBS(I)%TIMES(1:NTIME))
     DO j=1,NTIME
        FTEMP(J) = SECONDS(IntTime - TS_OBS(Id)%TIMES(J))
        FTEMP(J) = ABS(FTEMP(J))
     ENDDO
     DT_MIN = MINVAL(FTEMP(1:NTIME))
     TS_OBS(Id)%N_T_WEIGHT = MINLOC(FTEMP,DIM=1)

     IF(DT_MIN < T_THRESH)THEN     
       IF(DT_MIN < .5_SP*T_THRESH) THEN
         TS_OBS(Id)%T_WEIGHT = 1.0_SP
       ELSE
         TS_OBS(Id)%T_WEIGHT = (T_THRESH-DT_MIN)/T_THRESH*2.0_SP
       END IF
     END IF

     DEALLOCATE(FTEMP)
   END DO
   
       
!------------------------------------------------------------------------------!
!  Interpolate Simulation Data to Local Observation Point                      ! 
!------------------------------------------------------------------------------!
       
   ALLOCATE(SINT(N_ASSIM_TS,TS_MAX_LAYER)) ; SINT = 0.

!   IF(TS_NGASSIM)THEN
     DO I=1,N_ASSIM_TS_w
     ID = ID_STATION(I)   
       DO J=1,TS_OBS(Id)%N_INTPTS
         J1        = TS_OBS(Id)%INTPTS(J)
         WGHT      = TS_OBS(Id)%X_WEIGHT(J)
         NLAY      = TS_OBS(Id)%N_LAYERS
         DO K=1,NLAY
           V1 = SG(J1,TS_OBS(Id)%S_INT(K,1))
           V2 = SG(J1,TS_OBS(Id)%S_INT(K,2))
           W1 = TS_OBS(Id)%S_WEIGHT(K,1)
           W2 = TS_OBS(Id)%S_WEIGHT(K,2)
           SINT(I,K) = SINT(I,K) + (V1*W1 + V2*W2)*WGHT 
         END DO
       END DO
       TOT_WGHT = SUM(TS_OBS(Id)%X_WEIGHT(1:TS_OBS(Id)%N_INTPTS))
       SINT(I,1:NLAY) = SINT(I,1:NLAY)/TOT_WGHT
     END DO
!   ELSE IF(TS_OIASSIM)THEN
!     DO I=1,N_ASSIM_TS
!       IP = TS_OBS(I)%N_CELL
!       NLAY      = TS_OBS(I)%N_LAYERS
!       XSC = TS_OBS(I)%X-XCG(IP)
!       YSC = TS_OBS(I)%Y-YCG(IP)
!       JN1 = NVG(IP,1)
!       JN2 = NVG(IP,2)
!       JN3 = NVG(IP,3)
      
!       DO K=1,NLAY
!         K1 = TS_OBS(I)%S_INT(K,1)
!         K2 = TS_OBS(I)%S_INT(K,2)
!         COFT0 = AW0G(IP,1)*SG(JN1,K1)+AW0G(IP,2)*SG(JN2,K1)+AW0G(IP,3)*SG(JN3,K1)
!         COFTX = AWXG(IP,1)*SG(JN1,K1)+AWXG(IP,2)*SG(JN2,K1)+AWXG(IP,3)*SG(JN3,K1)
!         COFTY = AWYG(IP,1)*SG(JN1,K1)+AWYG(IP,2)*SG(JN2,K1)+AWYG(IP,3)*SG(JN3,K1)
!         U1 = COFT0 + COFTX*XSC + COFTY*YSC
!         COFT0 = AW0G(IP,1)*SG(JN1,K2)+AW0G(IP,2)*SG(JN2,K2)+AW0G(IP,3)*SG(JN3,K2)
!         COFTX = AWXG(IP,1)*SG(JN1,K2)+AWXG(IP,2)*SG(JN2,K2)+AWXG(IP,3)*SG(JN3,K2)
!         COFTY = AWYG(IP,1)*SG(JN1,K2)+AWYG(IP,2)*SG(JN2,K2)+AWYG(IP,3)*SG(JN3,K2)
!         U2 = COFT0 + COFTX*XSC + COFTY*YSC
!         W1 = TS_OBS(I)%S_WEIGHT(K,1)
!         W2 = TS_OBS(I)%S_WEIGHT(K,2)
!	 SINT(I,K) = (U1*W1 + U2*W2)
!       ENDDO
!     ENDDO   	 
!   ELSE
!     WRITE(IPT,*)'EITHER TS_NGASSIM OR TS_OIASSIM SHOULD BE TRUE'
!     CALL PSTOP
!   END IF  
!------------------------------------------------------------------------------!
!  Compute Local Correction by Interpolating Observed/Computed Defect          ! 
!------------------------------------------------------------------------------!

   ALLOCATE(TWGHT_S(MGL,KBM1))   ; TWGHT_S   = 0.
   ALLOCATE(SCORR(MGL,KBM1))   ; SCORR   = 0.

   IF(TS_NGASSIM)THEN
     DO I=1,N_ASSIM_TS_w
     ID = ID_STATION(I) 
       DO J=1,TS_OBS(Id)%N_INTPTS
         J1     = TS_OBS(Id)%INTPTS(J)
         ITIMEA  = TS_OBS(Id)%N_T_WEIGHT
         NLAY   = TS_OBS(Id)%N_LAYERS
         DO K=1,NLAY
           K1             = TS_OBS(Id)%S_INT(K,1)
           K2             = TS_OBS(Id)%S_INT(K,2)
           W1             = TS_OBS(Id)%S_WEIGHT(K,1)
           W2             = TS_OBS(Id)%S_WEIGHT(K,2)

           DEFECT         = TS_OBS(Id)%SAL(ITIMEA,K) - SINT(I,K)
           IF(ABS(DEFECT) < 20.0_SP) THEN             !quality control
             WEIGHT1        = TS_OBS(Id)%T_WEIGHT*TS_OBS(Id)%X_WEIGHT(J)*W1
             WEIGHT2        = TS_OBS(Id)%T_WEIGHT*TS_OBS(Id)%X_WEIGHT(J)*W2
             TWGHT_S(J1,K1) = TWGHT_S(J1,K1) + WEIGHT1   
             TWGHT_S(J1,K2) = TWGHT_S(J1,K2) + WEIGHT2   

             CORRECTION     = TS_GAMA*DEFECT

!          SCORR(J1,K1)   = SCORR(J1,K1) + CORRECTION*WEIGHT1
!          SCORR(J1,K2)   = SCORR(J1,K2) + CORRECTION*WEIGHT2
             SCORR(J1,K1)   = SCORR(J1,K1) + CORRECTION*WEIGHT1*WEIGHT1
             SCORR(J1,K2)   = SCORR(J1,K2) + CORRECTION*WEIGHT2*WEIGHT2

           ENDIF
         END DO
       END DO
     END DO

!------------------------------------------------------------------------------!
!  Nudge Simulation Data Using Local Corrections                               ! 
!------------------------------------------------------------------------------!

     DO I=1,MGL

       IF(SSTGRD_ASSIM) THEN
         K1 = 2
       ELSE 
         K1 = 1
       ENDIF
#      if defined (PWP) 
       IF(MLD_GLB(I)<0) THEN
         K1 = 2
       ELSE 
         K1 = MLD_GLB(I)+1
       ENDIF
#      endif

       DO K=K1,KBM1
         IF(DA_TS(I) == 1 .AND. TWGHT_S(I,K) > 1.0E-08)THEN
           SG(I,K) = SG(I,K) + DTI*TS_GALPHA*SCORR(I,K)/TWGHT_S(I,K)
         END IF
       END DO
     END DO

     DEALLOCATE(TWGHT_S,SCORR,SINT)
   ELSE IF(TS_OIASSIM)THEN
     DO I=1,N_ASSIM_TS_w 
     ID = ID_STATION(I)
       ITIMEA  = TS_OBS(Id)%N_T_WEIGHT
       NLAY   = TS_OBS(Id)%N_LAYERS
       DO K=1,NLAY
         K1             = TS_OBS(Id)%S_INT(K,1)
         K2             = TS_OBS(Id)%S_INT(K,2)
         W1             = TS_OBS(Id)%S_WEIGHT(K,1)
         W2             = TS_OBS(Id)%S_WEIGHT(K,2)

         DEFECT         = TS_OBS(Id)%SAL(ITIMEA,K) - SINT(I,K)
         IF(ABS(DEFECT) < 20.0_SP) THEN             !quality control
           WEIGHT1        = TS_OBS(Id)%T_WEIGHT*W1
           WEIGHT2        = TS_OBS(Id)%T_WEIGHT*W2
           TWGHT_S(I,K1) = TWGHT_S(I,K1) + WEIGHT1   
           TWGHT_S(I,K2) = TWGHT_S(I,K2) + WEIGHT2   
           CORRECTION     = DEFECT
           SCORR(I,K1)   = SCORR(I,K1) + CORRECTION*WEIGHT1
           SCORR(I,K2)   = SCORR(I,K2) + CORRECTION*WEIGHT2
         ENDIF
       END DO
     END DO

     DO I=1,N_ASSIM_TS_w
       DO K=1,KBM1
         IF(TWGHT_S(I,K) > 1.0E-8)THEN
           SCORR(I,K)=SCORR(I,K)/TWGHT_S(I,K)
         END IF
       END DO
     END DO

!------------------------------------------------------------------------------!
!  'OI' Simulation Data Using Local Corrections                               ! 
!------------------------------------------------------------------------------!

     ALLOCATE(SCORR1(MGL,KBM1));SCORR1=0.0_SP

     DO K=1,KBM1
       CALL TS_OPTIMINTERP(SCORR(:,K),SCORR1(:,K))
     END DO
   
     DO I=1,MGL

       IF(SSTGRD_ASSIM) THEN
         K1 = 2
       ELSE 
         K1 = 1
       ENDIF
#      if defined (PWP) 
       IF(MLD_GLB(I)<0) THEN
         K1 = 2
       ELSE 
         K1 = MLD_GLB(I)+1
       ENDIF
#      endif

       DO K=K1,KBM1
         IF(DA_TS(I) == 1)THEN
!QXU       SG(I,K)= SG(I,K) + SCORR1(I,K)
           SG(I,K) = SG(I,K) + TS_OIGALPHA*SCORR1(I,K)
         END IF
       END DO
     END DO

     DEALLOCATE(TWGHT_S,SCORR,SINT,SCORR1)
   ELSE
     WRITE(IPT,*)'EITHER TS_NGASSIM OR TS_OIASSIM SHOULD BE TRUE'
     CALL PSTOP
   END IF

   END IF  !!MASTER
   
!------------------------------------------------------------------------------!
!  Disperse New Data Fields to Slave Processors                                ! 
!------------------------------------------------------------------------------!
   IF(SERIAL)THEN
     SF1(1:M,1:KBM1) = SG(1:M,1:KBM1)
   END IF
#  if defined (MULTIPROCESSOR) 
   CALL MPI_BCAST(SG,(MGL+1)*KB,MPI_F,0,MPI_COMM_WORLD,IERR)
   IF(PAR)THEN
     DO I=1,M
       SF1(I,1:KBM1) = SG(NGID(I),1:KBM1)
     END DO
   END IF
#  endif

   DEALLOCATE(SG)
#  if defined (PWP)
   DEALLOCATE(MLD_GLB)
#  endif
          
   IF (DBG_SET(DBG_SBR)) WRITE(IPT,*) "END: SALT_ASSIMILATION" 

   RETURN
   END SUBROUTINE SALT_ASSIMILATION
!==============================================================================|
!==============================================================================|

!==============================================================================|
!==============================================================================|
   SUBROUTINE TS_OPTIMINTERP(F,FI)
   USE MOD_OPTIMAL_INTERPOLATION
   IMPLICIT NONE

!------------------------------------------------------------------------------|
! xi(1,:) and xi(2,:) represent the x and y coordindate of the grid of the     |
! interpolated field                                                           |
! fi and vari are the interpolated field and its error variance resp.          |
!------------------------------------------------------------------------------|
   REAL(SP) :: XI(2,MGL),FI(MGL),VARI(MGL)

!------------------------------------------------------------------------------|
! x(1,:) and x(2,:) represent the x and y coordindate of the observations      |
! f and var are observations and their error variance resp.                    |
!------------------------------------------------------------------------------|
!   REAL(SP) :: X(2,N_ASSIM_TS),VAR(N_ASSIM_TS),F(N_ASSIM_TS)
  REAL(SP) :: X(2,N_ASSIM_TS_W),VAR(N_ASSIM_TS_W),F(N_ASSIM_TS_W)

!------------------------------------------------------------------------------|
! param: inverse of the correlation length                                     |
!------------------------------------------------------------------------------|
   REAL(SP) :: PARAM(2,N_ASSIM_TS)

   INTEGER  :: I,J,MM   
   INTEGER ID
!------------------------------------------------------------------------------|
! create a regular 2D grid                                                     |
!------------------------------------------------------------------------------|
   DO I=1,MGL
     XI(1,I) = XG(I)
     XI(2,I) = YG(I)
   END DO	

!------------------------------------------------------------------------------|   
! param is the inverse of the correlation length                               |
!------------------------------------------------------------------------------|
   PARAM = 1.0_SP/TS_PARAM             
 
   MM = TS_N_INFLU

!------------------------------------------------------------------------------|   
! the error variance of the observations                                       |
!------------------------------------------------------------------------------|
   VAR = 0.0_SP   

!------------------------------------------------------------------------------|
! location of observations                                                     |
!------------------------------------------------------------------------------|
   DO I=1,N_ASSIM_TS_w
     ID = ID_STATION(I)
     X(1,I) = TS_OBS(Id)%X
     X(2,I) = TS_OBS(Id)%Y
   END DO  

!------------------------------------------------------------------------------|
! fi is the interpolated function and vari its error variance                  |
!------------------------------------------------------------------------------|
   CALL OPTIMINTERP(X,F,VAR,PARAM,MM,XI,FI,VARI)

   RETURN
   END SUBROUTINE TS_OPTIMINTERP
!==============================================================================|
!==============================================================================|



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   SUBROUTINE SSTGRD_OBSERVATION_UPDATE
!------------------------------------------------------------------------------!
   IMPLICIT NONE

   CALL NC_READ_VAR(VAR_SST,SST_FILE%FTIME%NEXT_STKCNT)

   SST_FILE%FTIME%NEXT_STKCNT = SST_FILE%FTIME%NEXT_STKCNT + 1


   ! RESET THE INDEX AND THE TIME FOR USE IN THE ASSIMILATION CYCLE
   SST_SAVE_INDEX = 1
   SST_SAVE_TIME = IntTime_BUF + SST_SAVE_INTERVAL
   

   END SUBROUTINE SSTGRD_OBSERVATION_UPDATE
!==============================================================================|

   SUBROUTINE SSHGRD_OBSERVATION_UPDATE
!------------------------------------------------------------------------------!
   IMPLICIT NONE

   !! specified to certain month
   CALL NC_READ_VAR(VAR_SSH,SSH_FILE%FTIME%NEXT_STKCNT)

   SSH_FILE%FTIME%NEXT_STKCNT = SSH_FILE%FTIME%NEXT_STKCNT + 1


   ! RESET THE INDEX AND THE TIME FOR USE IN THE ASSIMILATION CYCLE
   SSH_SAVE_INDEX = 1
   SSH_SAVE_TIME = IntTime_BUF + SSH_SAVE_INTERVAL
   

   END SUBROUTINE SSHGRD_OBSERVATION_UPDATE

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
   SUBROUTINE SSTGRD_NUDGE
  
!==============================================================================|
!  USE SST OBSERVATION DATA TO NUDGE SURFACE TEMPERATURE                       |
!==============================================================================|
#  if defined (PWP)
   USE MOD_PWP, only : MLD
#  endif
   IMPLICIT NONE
   REAL(SP) :: TDIFF,T_WEIGHT,TRUTH,ADJUSTMENT, S_AVG
   INTEGER I,K

!==============================================================================|


!------------------------------------------------------------------------------!
!  Calculate Temporal Weight of Measurement (I) at Time(TIME)                  ! 
!------------------------------------------------------------------------------!
# if defined (ENKF)
     SST_SAVED(1:M,:)=SST_SAVED_ENKF(1:M,:,ENKF_memberCONTR)
     SST_AVG(:)=SST_AVG_ENKF(:,ENKF_memberCONTR)
# endif

   IF(IntTime > SST_SAVE_TIME + SST_SAVE_INTERVAL/2) THEN
   
      SST_SAVE_INDEX = SST_SAVE_INDEX + 1

      SST_SAVE_TIME = SST_SAVE_TIME + SST_SAVE_INTERVAL
      
      IF(DBG_SET(DBG_LOG)) THEN
         WRITE(IPT,*) "Setting sst state # ",SST_SAVE_INDEX
!         CALL PRINT_REAL_TIME(SST_SAVE_TIME,IPT,"Set at")
      END IF

   END IF

   ! GET TIME DIFFERENCE IN SECONDS
   TDIFF = SECONDS(IntTime - SST_SAVE_TIME)
   TDIFF = ABS(TDIFF)

   IF(TDIFF < 0.5*SSTGRD_TIME_WINDOW) THEN
     T_WEIGHT = 1.0
   ELSE IF(TDIFF < SSTGRD_TIME_WINDOW) THEN
     T_WEIGHT = (SSTGRD_TIME_WINDOW-TDIFF)/SSTGRD_TIME_WINDOW*2.0
   ELSE
     T_WEIGHT = 0.0
   END IF
  
   IF(T_WEIGHT == 0.0) RETURN

   DO I=1,M
# if defined (WET_DRY)
     IF(WETTING_DRYING_ON.and.ISWETN(i)==0)cycle
# endif
      !TRUTH = MODEL_SAVED@INTERVAL + OBSERVED_AVG  - MODEL_AVG
      TRUTH = SST_SAVED(I,SST_SAVE_INDEX) + SST_OBS(I) - SST_AVG(I)
      ! ADJUSMENT = MAX_WEIGHT * (TRUTH-CURRENT) * (WEIGHT_BY_TIME)
      ADJUSTMENT = SSTGRD_WEIGHT_MAX*(TRUTH-T1(i,1))*T_WEIGHT
      ! RESULT = CURRENT + DTI*TIMESCALE * ADJUSTMENT 
      T1(I,1)=T1(I,1)+DTI*SSTGRD_TIMESCALE*ADJUSTMENT

#     if defined (PWP)
      IF(MLD(I)>0) THEN
        S_AVG = 0.0_SP
        DO K=1, MLD(I)
          S_AVG = S_AVG + S1(I,K)/FLOAT(MLD(I))
        ENDDO
        S1(I,1) = S_AVG
        DO K=2,MLD(I)
          T1(I,K)=T1(I,1)
          S1(I,K)=S_AVG
        END DO
      ENDIF
#     endif

   ENDDO    

#    if defined (MULTIPROCESSOR)
     IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,T1) !Interprocessor Exchange   !
     IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,S1)
#    endif

   
   CALL N2E3D(T1,T)
   CALL N2E3D(S1,S)

   RETURN
   END SUBROUTINE SSTGRD_NUDGE
!==============================================================================|

   SUBROUTINE SSHGRD_NUDGE
  
!==============================================================================|
!  USE SST OBSERVATION DATA TO NUDGE SURFACE TEMPERATURE                       |
!==============================================================================|
   IMPLICIT NONE
   REAL(SP) :: TDIFF,T_WEIGHT,TRUTH,ADJUSTMENT
   INTEGER I,K
   REAL(SP) :: F_SPATIAL, SSH_SAVED_R
   REAL(SP) :: R0, R1, R2, VX_TMP

!==============================================================================|

# if defined (ENKF)

     SSH_SAVED(1:M,:)=SSH_SAVED_ENKF(1:M,:,ENKF_memberCONTR)
     SSH_AVG(:)=SSH_AVG_ENKF(:,ENKF_memberCONTR)
# endif

!------------------------------------------------------------------------------!
!  Calculate Temporal Weight of Measurement (I) at Time(TIME)                  ! 
!------------------------------------------------------------------------------!
   IF(IntTime > SSH_SAVE_TIME + SSH_SAVE_INTERVAL/2) THEN
   
      SSH_SAVE_INDEX = SSH_SAVE_INDEX + 1

      SSH_SAVE_TIME = SSH_SAVE_TIME + SSH_SAVE_INTERVAL
      
      IF(DBG_SET(DBG_LOG)) THEN
         WRITE(IPT,*) "Setting ssh state # ",SSH_SAVE_INDEX
!         CALL PRINT_REAL_TIME(SSH_SAVE_TIME,IPT,"Set at")
      END IF

   END IF

   ! GET TIME DIFFERENCE IN SECONDS
   TDIFF = SECONDS(IntTime - SSH_SAVE_TIME)
   TDIFF = ABS(TDIFF)

   IF(TDIFF < 0.5*SSHGRD_TIME_WINDOW) THEN
     T_WEIGHT = 1.0
   ELSE IF(TDIFF < SSHGRD_TIME_WINDOW) THEN
     T_WEIGHT = (SSHGRD_TIME_WINDOW-TDIFF)/SSHGRD_TIME_WINDOW*2.0
   ELSE
     T_WEIGHT = 0.0
   END IF
  
   IF(T_WEIGHT == 0.0) RETURN

   R0 = SECONDS(IntTime)
   R1 = SECONDS(SSH_SAVE_TIME)
   R2 = SECONDS(SSH_SAVE_INTERVAL)
   DO I=1,M

#     if defined (SPHERICAL)
      F_SPATIAL = 1.0_SP
#     endif

      IF(IntTime>SSH_SAVE_TIME .AND. IntTime<= SSH_SAVE_TIME + SSH_SAVE_INTERVAL/2 ) THEN
        SSH_SAVED_R = ( SSH_SAVED(I,SSH_SAVE_INDEX)*(R1+R2-R0)+SSH_SAVED(I,SSH_SAVE_INDEX+1)*(R0-R1) )/R2
      ELSE
        SSH_SAVED_R = ( SSH_SAVED(I,SSH_SAVE_INDEX-1)*(R1-R0)+SSH_SAVED(I,SSH_SAVE_INDEX)*(R0-R1+R2) )/R2 
      ENDIF

      !TRUTH = MODEL_SAVED@INTERVAL + OBSERVED_AVG  - MODEL_AVG
!      TRUTH = SSH_SAVED(I,SSH_SAVE_INDEX) + SSH_OBS(I) - SSH_AVG(I)
      TRUTH = SSH_SAVED_R + SSH_OBS(I) - SSH_AVG(I)
      ! ADJUSMENT = MAX_WEIGHT * (TRUTH-CURRENT) * (WEIGHT_BY_TIME)
      ADJUSTMENT = SSHGRD_WEIGHT_MAX*(TRUTH-EL(i))*T_WEIGHT*F_SPATIAL
      ! RESULT = CURRENT + DTI*TIMESCALE * ADJUSTMENT 
      EL(I)=EL(I)+DTI*SSHGRD_TIMESCALE*ADJUSTMENT
   ENDDO    

   
#  if defined (MULTIPROCESSOR)
   IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,EL) !Interprocessor Exchange   !
#  endif

   CALL N2E2D(EL,EL1)

#   if defined (ICE_EMBEDDING)
  DO I = 1, NT                                !!yding
    QTHICE1(I) = (QTHICE(NV(I,1))+QTHICE(NV(I,2))+QTHICE(NV(I,3)))/3.0_SP
  END DO

  D   = H + EL - QTHICE     !!yding
  D1  = H1 + EL1 - QTHICE1  !!yding

#   else

   D   = H + EL
   D1  = H1 + EL1

#   endif

   RETURN
   END SUBROUTINE SSHGRD_NUDGE

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   SUBROUTINE TSGRD_OBSERVATION_UPDATE(NOW)
!------------------------------------------------------------------------------!
   IMPLICIT NONE

   TYPE(TIME), INTENT(IN) :: NOW
   TYPE(TIME)             :: WTIME
   TYPE(NCFTIME), POINTER :: FTM
   REAL(SP), POINTER      :: VNP_ARR(:,:), VPP_ARR(:,:)
   INTEGER :: STATUS

!   TS_FILE%FTIME%NEXT_STKCNT = Pmonth  !!!  the month

!   CALL NC_READ_VAR(VAR_TC,TS_FILE%FTIME%NEXT_STKCNT)
!   CALL NC_READ_VAR(VAR_SC,TS_FILE%FTIME%NEXT_STKCNT)

!!   TS_FILE%FTIME%NEXT_STKCNT = TS_FILE%FTIME%NEXT_STKCNT + 1

   WTIME = NOW

   FTM => TS_FILE%FTIME

   ! T OBSERVATION UPDATE
   CALL UPDATE_VAR_BRACKET(TS_FILE,TSC_T1_P,TSC_T1_N,WTIME,STATUS)
   IF (STATUS /= 0) THEN
     CALL FATAL_ERROR("COULD NOT UPATE TS BRACKET: BOUNDS EXCEEDED?")
   end if
    
   CALL NC_POINT_VAR(TSC_T1_N,VNP_ARR)
   CALL NC_POINT_VAR(TSC_T1_P,VPP_ARR)   
   TC_OBS = FTM%NEXT_WGHT * VNP_ARR + FTM%PREV_WGHT * VPP_ARR

   ! S OBSERVATION UPDATE
   CALL UPDATE_VAR_BRACKET(TS_FILE,TSC_S1_P,TSC_S1_N,WTIME,STATUS)
   IF (STATUS /= 0) THEN
     CALL FATAL_ERROR("COULD NOT UPATE TS BRACKET: BOUNDS EXCEEDED?")
   end if
    
   CALL NC_POINT_VAR(TSC_S1_N,VNP_ARR)
   CALL NC_POINT_VAR(TSC_S1_P,VPP_ARR)   
   SC_OBS = FTM%NEXT_WGHT * VNP_ARR + FTM%PREV_WGHT * VPP_ARR

   ! RESET THE INDEX AND THE TIME FOR USE IN THE ASSIMILATION CYCLE
   TSC_SAVE_INDEX = 1
   TSC_SAVE_TIME = IntTime_BUF + TSC_SAVE_INTERVAL

   END SUBROUTINE TSGRD_OBSERVATION_UPDATE

!==============================================================================|

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
   SUBROUTINE TSGRD_NUDGE
  
!==============================================================================|
!  USE TS OBSERVATION DATA TO NUDGE SURFACE TEMPERATURE                       |
!==============================================================================|
#  if defined (PWP)
   USE MOD_PWP, only : MLD
#  endif
   IMPLICIT NONE
   REAL(SP) :: TDIFF,T_WEIGHT,TRUTH,ADJUSTMENT
   INTEGER I,K,K1

!==============================================================================|

!     IF(.false.) THEN

!------------------------------------------------------------------------------!
!  Calculate Temporal Weight of Measurement (I) at Time(TIME)                  ! 
!------------------------------------------------------------------------------!
   IF(IntTime > TSC_SAVE_TIME + TSC_SAVE_INTERVAL/2) THEN
   
      TSC_SAVE_INDEX = TSC_SAVE_INDEX + 1

      TSC_SAVE_TIME = TSC_SAVE_TIME + TSC_SAVE_INTERVAL
      
      IF(DBG_SET(DBG_LOG)) THEN
         WRITE(IPT,*) "Setting T/S state # ",TSC_SAVE_INDEX
!         CALL PRINT_REAL_TIME(TSC_SAVE_TIME,IPT,"Set at")
      END IF

   END IF

   ! GET TIME DIFFERENCE IN SECONDS
   TDIFF = SECONDS(IntTime - TSC_SAVE_TIME)
   TDIFF = ABS(TDIFF)

   IF(TDIFF < 0.5*TSGRD_TIME_WINDOW) THEN
     T_WEIGHT = 1.0
   ELSE IF(TDIFF < TSGRD_TIME_WINDOW) THEN
     T_WEIGHT = (TSGRD_TIME_WINDOW-TDIFF)/TSGRD_TIME_WINDOW*2.0
   ELSE
     T_WEIGHT = 0.0
   END IF
  
   IF(T_WEIGHT == 0.0) RETURN

   DO I=1,M

      IF(SSTGRD_ASSIM) THEN
        K1 = 2
      ELSE
        K1 = 1
      ENDIF
#     if defined (PWP)
      IF(MLD(I)<0) THEN
        K1 = 2
      ELSE
        K1 = MLD(I)+1
      ENDIF
#     endif

      DO K=K1,KBM1
      !TRUTH = MODEL_SAVED@INTERVAL + OBSERVED_AVG  - MODEL_AVG
      TRUTH = TC_SAVED(I,K,TSC_SAVE_INDEX) + TC_OBS(I,K) - TC_AVG(I,K)
      ! ADJUSMENT = MAX_WEIGHT * (TRUTH-CURRENT) * (WEIGHT_BY_TIME)
      ADJUSTMENT = TSGRD_WEIGHT_MAX*(TRUTH-T1(I,K))*T_WEIGHT
      ! RESULT = CURRENT + DTI*TIMESCALE * ADJUSTMENT 
      T1(I,K)=T1(I,K)+DTI*TSGRD_TIMESCALE*ADJUSTMENT
      ENDDO
   ENDDO

   DO I=1,M
      DO K=1, KBM1
      !TRUTH = MODEL_SAVED@INTERVAL + OBSERVED_AVG  - MODEL_AVG
      TRUTH = SC_SAVED(I,K,TSC_SAVE_INDEX) + SC_OBS(I,K) - SC_AVG(I,K)
      ! ADJUSMENT = MAX_WEIGHT * (TRUTH-CURRENT) * (WEIGHT_BY_TIME)
      ADJUSTMENT = TSGRD_WEIGHT_MAX*(TRUTH-S1(I,K))*T_WEIGHT
      ! RESULT = CURRENT + DTI*TIMESCALE * ADJUSTMENT 
      S1(I,K)=S1(I,K)+DTI*TSGRD_TIMESCALE*ADJUSTMENT
      END DO
   ENDDO    

!   DO I=1,M
!      DO K=1,KBM1
!      !TRUTH = MODEL_SAVED@INTERVAL + OBSERVED_AVG  - MODEL_AVG
!      TRUTH = SC_SAVED(I,K,TSC_SAVE_INDEX) + SC_OBS(I,K) - SC_AVG(I,K)
!      ! ADJUSMENT = MAX_WEIGHT * (TRUTH-CURRENT) * (WEIGHT_BY_TIME)
!      ADJUSTMENT = TSGRD_WEIGHT_MAX*(TRUTH-S1(i,k))*T_WEIGHT
!      ! RESULT = CURRENT + DTI*TIMESCALE * ADJUSTMENT 
!      S1(I,k)=S1(I,k)+DTI*TSGRD_TIMESCALE*ADJUSTMENT
!      END DO
!   ENDDO      
   
!    endif !! change to the new methon

!!  Average the simulation TS over the assimulation perion
  
!------------------------------------------------------------------------------!
!  Calculate Temporal Weight of Measurement (I) at Time(TIME)                  ! 
!------------------------------------------------------------------------------!


   ! GET TIME DIFFERENCE IN SECONDS
!   TDIFF = SECONDS(IntTime - TSC_SAVE_TIME)
!   TDIFF = ABS(TDIFF)
!   !TDIFF = ABS(TDIFF)/(10.0_SP*144.*dti)  !!  5 days windows 5.0*
!   TDIFF = ABS(TDIFF)/(10.0_SP*86400.)  !! 20 days windows 5.0*

!   !T_WEIGHT =DTI * ISQPI2*EXP(-(TDIFF*TDIFF))
!   T_WEIGHT = ISQPI2*EXP(-(TDIFF*TDIFF)/2.0_SP)
!   !T_WEIGHT = T_WEIGHT/(10.0_SP*144.)
!   T_WEIGHT = T_WEIGHT/(10.0_SP*86400./dti)
   
!   DO I=1,M
!      DO K=1,KBM1
!      TRUTH = TC_OBS(I,K) - TC0_AVG(I,K)
!      ADJUSTMENT = TRUTH*T_WEIGHT
!      T1(I,k)=T1(I,k)+ADJUSTMENT
      
!      TRUTH = SC_OBS(I,K) - SC0_AVG(I,K)
!      ADJUSTMENT = TRUTH*T_WEIGHT
!      S1(I,k)=S1(I,k)+ADJUSTMENT
!      END DO
!   ENDDO    
!!  ggao 2009
!     WHERE(S1 < 0.0_SP)S1=0.0_SP

   
#  if defined (MULTIPROCESSOR)
   IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,T1) !Interprocessor Exchange   !
   IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,S1) !Interprocessor Exchange   !
#  endif

   
   CALL N2E3D(T1,T)
   CALL N2E3D(S1,S)
   

   RETURN
   END SUBROUTINE TSGRD_NUDGE
!==============================================================================|
  SUBROUTINE SET_SST_SAVE_TIME_INDEX(FLAG)
   IMPLICIT NONE
   INTEGER I,FLAG
    TYPE(TIME) :: TEMP_TIME
  IF (FLAG==0) THEN
    TEMP_TIME%MJD=Inttime%mjd
    TEMP_TIME%musod=0.0
   DO I=1,SSTGRD_N_PER_INTERVAL
      IF  (inttime<SST_SAVE_INTERVAL*I+TEMP_TIME) THEN
           SST_SAVE_TIME=SST_SAVE_INTERVAL*I+TEMP_TIME
           SST_SAVE_INDEX=I
           GOTO 386
      END IF
   END DO 
   PRINT *, "FAILED TO FIND THE RIGHT SST SAVE TIME INDEX , STOP"
   PRINT *, 'CURRERNT INTTIME IS: ', INTTIME
   CALL PSTOP
386 RETURN
  ELSEIF (FLAG==1) THEN
    TEMP_TIME%MJD=Inttime%mjd
    TEMP_TIME%musod=0.0
   DO I=1,SSTGRD_N_PER_INTERVAL
      IF  (inttime<=SST_SAVE_INTERVAL*I+TEMP_TIME+SST_SAVE_INTERVAL/2.0) THEN
           SST_SAVE_TIME=SST_SAVE_INTERVAL*I+TEMP_TIME
           SST_SAVE_INDEX=I
           GOTO 387
      END IF
   END DO
   PRINT *, "FAILED TO FIND THE RIGHT SST SAVE TIME INDEX , STOP"
   PRINT *, 'CURRERNT INTTIME IS: ', INTTIME
   CALL PSTOP
387 RETURN
  ELSE
     PRINT *, 'FLAG IN SET_SST_SAVE_TIME_INDEX CANNOT BE :',FLAG
     CALL PSTOP
  END IF
  
  END SUBROUTINE SET_SST_SAVE_TIME_INDEX


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
   SUBROUTINE SST_SAVE     
!==============================================================================|
!  SAVE INTERVAL-AVERAGAED SST DATA ON GRID TO USE AS OBSERVATION              | 
!==============================================================================|
   IMPLICIT NONE
   INTEGER I
   
!==============================================================================|

   IF(IntTime >= SST_SAVE_TIME) THEN
   
      SST_SAVE_INDEX = SST_SAVE_INDEX + 1
      
      IF(DBG_SET(DBG_LOG)) THEN
         WRITE(IPT,*) "Saving sst state # ",SST_SAVE_INDEX
      END IF
!------------------------------------------------------------------------------!
!  Save Hourly SST on grid in Simulation running                               ! 
!------------------------------------------------------------------------------!
      SST_SAVED(1:M,SST_SAVE_INDEX) = T1(1:M,1)
# if defined (ENKF)
     SST_SAVED_ENKF(1:M,SST_SAVE_INDEX,ENKF_memberCONTR)=SST_SAVED(1:M,SST_SAVE_INDEX)
# endif
      ! INCRIMENT THE SAVE_TIME
      SST_SAVE_TIME = SST_SAVE_TIME + SST_SAVE_INTERVAL

!------------------------------------------------------------------------------!
!  Average Hourly Data Over Integration Period                                 ! 
!------------------------------------------------------------------------------!

      IF(SST_SAVE_INDEX == SSTGRD_N_PER_INTERVAL )THEN
      
         DO I=1,M
            SST_AVG(I) = SUM(SST_SAVED(I,1:SSTGRD_N_PER_INTERVAL))&
                 &/real(SSTGRD_N_PER_INTERVAL,SP)
# if defined (ENKF)          
            SST_AVG_ENKF(I,ENKF_memberCONTR)=SST_AVG(I)
# endif

         END DO
         
      END IF

   END IF
   
   RETURN
   END SUBROUTINE SST_SAVE

!==============================================================================|
  SUBROUTINE SET_SSH_SAVE_TIME_INDEX(FLAG)
   IMPLICIT NONE
   INTEGER I,FLAG
    TYPE(TIME) :: TEMP_TIME
  IF (FLAG==0) THEN
    TEMP_TIME%MJD=Inttime%mjd
    TEMP_TIME%musod=0.0
   DO I=1,SSHGRD_N_PER_INTERVAL
      IF  (inttime<SSH_SAVE_INTERVAL*I+TEMP_TIME) THEN
           SSH_SAVE_TIME=SSH_SAVE_INTERVAL*I+TEMP_TIME
           SSH_SAVE_INDEX=I
           GOTO 386
      END IF
   END DO
   PRINT *, "FAILED TO FIND THE RIGHT SSH SAVE TIME INDEX , STOP"
   PRINT *, 'CURRERNT INTTIME IS: ', INTTIME
   CALL PSTOP
386 RETURN
  ELSEIF (FLAG==1) THEN
        TEMP_TIME%MJD=Inttime%mjd
    TEMP_TIME%musod=0.0
   DO I=1,SSHGRD_N_PER_INTERVAL
      IF  (inttime<=SSH_SAVE_INTERVAL*I+TEMP_TIME+SSH_SAVE_INTERVAL/2.0) THEN
           SSH_SAVE_TIME=SSH_SAVE_INTERVAL*I+TEMP_TIME
           SSH_SAVE_INDEX=I
           GOTO 387
      END IF
   END DO
   PRINT *, "FAILED TO FIND THE RIGHT SSH SAVE TIME INDEX , STOP"
   PRINT *, 'CURRERNT INTTIME IS: ', INTTIME
   CALL PSTOP
387 RETURN
  ELSE
    PRINT *, 'FLAG IN SET_SSH_SAVE_TIME_INDEX CANNOT BE: ',FLAG
    CALL PSTOP
  END IF
  END SUBROUTINE SET_SSH_SAVE_TIME_INDEX

   SUBROUTINE SSH_SAVE     
!==============================================================================|
!  SAVE INTERVAL-AVERAGAED SSH DATA ON GRID TO USE AS OBSERVATION              | 
!==============================================================================|
   IMPLICIT NONE
   INTEGER I
   
!==============================================================================|

   IF(IntTime >= SSH_SAVE_TIME) THEN
   
      SSH_SAVE_INDEX = SSH_SAVE_INDEX + 1
      
      IF(DBG_SET(DBG_LOG)) THEN
         WRITE(IPT,*) "Saving ssh state # ",SSH_SAVE_INDEX
      END IF
!------------------------------------------------------------------------------!
!  Save Hourly SSH on grid in Simulation running                               ! 
!------------------------------------------------------------------------------!
      SSH_SAVED(1:M,SSH_SAVE_INDEX) = EL(1:M)
# if defined (ENKF)
     SSH_SAVED_ENKF(1:M,SSH_SAVE_INDEX,ENKF_memberCONTR)=SSH_SAVED(1:M,SSH_SAVE_INDEX)
# endif
      
      ! INCRIMENT THE SAVE_TIME
      SSH_SAVE_TIME = SSH_SAVE_TIME + SSH_SAVE_INTERVAL

!------------------------------------------------------------------------------!
!  Average Hourly Data Over Integration Period                                 ! 
!------------------------------------------------------------------------------!

      IF(SSH_SAVE_INDEX == SSHGRD_N_PER_INTERVAL )THEN
      
         DO I=1,M
            SSH_AVG(I) = SUM(SSH_SAVED(I,1:SSHGRD_N_PER_INTERVAL))&
                 &/real(SSHGRD_N_PER_INTERVAL,SP)
# if defined (ENKF)
            SSH_AVG_ENKF(I,ENKF_memberCONTR)=SSH_AVG(I)
# endif

         END DO
         
      END IF

   END IF
   
   RETURN
 END SUBROUTINE SSH_SAVE

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
   SUBROUTINE TSGRD_SAVE     
!==============================================================================|
!  SAVE INTERVAL-AVERAGAED TS DATA ON GRID TO USE AS OBSERVATION              | 
!==============================================================================|
!   USE MOD_PREC
!   USE ALL_VARS
!   USE MOD_PAR
!   IMPLICIT NONE
!   INTEGER I,K
   
!==============================================================================|

!   IF(IntTime >= TSC_SAVE_TIME) THEN
   
!      TSC_SAVE_INDEX = TSC_SAVE_INDEX + 1
      
!      IF(DBG_SET(DBG_LOG)) THEN
!         WRITE(IPT,*) "Saving ts state # ",TSC_SAVE_INDEX
!      END IF
!------------------------------------------------------------------------------!
!  Save Hourly TC on grid in Simulation running                               ! 
!------------------------------------------------------------------------------!
!      TC_SAVED(1:M,1:KBM1,TSC_SAVE_INDEX) = T1(1:M,1:KBM1)
!      SC_SAVED(1:M,1:KBM1,TSC_SAVE_INDEX) = S1(1:M,1:KBM1)

      ! INCRIMENT THE SAVE_TIME
!!      TSC_SAVE_TIME = TSC_SAVE_TIME + TSC_SAVE_INTERVAL

!------------------------------------------------------------------------------!
!  Average Hourly Data Over Integration Period                                 ! 
!------------------------------------------------------------------------------!

!      IF(TSC_SAVE_INDEX == TSGRD_N_PER_INTERVAL )THEN
      
!         DO I=1,M
!            DO K=1,KBM1
!            TC_AVG(I,K) = SUM(TC_SAVED(I,K,1:TSGRD_N_PER_INTERVAL))&
!                 &/real(TSGRD_N_PER_INTERVAL,SP)
!            SC_AVG(I,K) = SUM(SC_SAVED(I,K,1:TSGRD_N_PER_INTERVAL))&
!                 &/real(TSGRD_N_PER_INTERVAL,SP)
!            END DO
!         END DO
         
!      END IF

!   END IF

   IMPLICIT NONE
   INTEGER I,K
   
!==============================================================================|

   IF(ASSIM_FLAG==0) THEN

     IF(IntTime >= TSC_SAVE_TIME) THEN

       TSC_SAVE_INDEX = TSC_SAVE_INDEX + 1
      
       IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) "Saving T/S hourly state # ",TSC_SAVE_INDEX
!------------------------------------------------------------------------------!
!  Save Hourly TC on grid in Simulation running                               ! 
!------------------------------------------------------------------------------!
       TC_SAVED(1:M,1:KBM1,TSC_SAVE_INDEX) = T1(1:M,1:KBM1)
       SC_SAVED(1:M,1:KBM1,TSC_SAVE_INDEX) = S1(1:M,1:KBM1)

       ! INCRIMENT THE SAVE_TIME
       TSC_SAVE_TIME = TSC_SAVE_TIME + TSC_SAVE_INTERVAL

!------------------------------------------------------------------------------!
!  Average Hourly Data Over Integration Period                                 ! 
!------------------------------------------------------------------------------!
   
       IF(TSC_SAVE_INDEX == TSGRD_N_PER_INTERVAL )THEN
      
         DO I=1,M
           DO K=1,KBM1
             TC_AVG(I,K) = SUM(TC_SAVED(I,K,1:TSGRD_N_PER_INTERVAL))&
                         &/real(TSGRD_N_PER_INTERVAL,SP)
             SC_AVG(I,K) = SUM(SC_SAVED(I,K,1:TSGRD_N_PER_INTERVAL))&
                         &/real(TSGRD_N_PER_INTERVAL,SP)
           END DO
         END DO
         
       END IF

     END IF

   ENDIF

!------------------------------------------------------------------------------!
!  Save Monthly averaged T/S on grid in Assimilation run                       !
!------------------------------------------------------------------------------!
   IF(ASSIM_FLAG==1) THEN

     DO I=1,M
       DO K=1,KBM1
         TC0_AVG(I,K) = TC0_AVG(I,K)+T1(I,K)
         SC0_AVG(I,K) = SC0_AVG(I,K)+S1(I,K)
       END DO
     END DO

!------------------------------------------------------------------------------!
!  Average Hourly Data Over Integration Period                                 !
!------------------------------------------------------------------------------!
     IF(IntTime >= TSC_SAVE_TIME2 )THEN
       IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) "Save averaged T/S state "
       DO I=1,M
         DO K=1,KBM1
           TC0_AVG(I,K) = TC0_AVG(I,K)/(real(Pmdays,kind=DP)*24.*3600./DTI)
           SC0_AVG(I,K) = SC0_AVG(I,K)/(real(Pmdays,kind=DP)*24.*3600./DTI)
         END DO
       END DO
     END IF

   ENDIF
   
   RETURN
   END SUBROUTINE TSGRD_SAVE


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
!==============================================================================|
!   READ IN DATA FROM SIMULATION STAGE AND FOR ASSIMILATION STAGE STARTUP      |
!==============================================================================|

   SUBROUTINE RESTORE_STATE

    !------------------------------------------------------------------------------|
    USE EQS_OF_STATE
# if defined (ICE)
    USE mod_ICE
    USE mod_ICE2D, only: sig1, sig2, sig12
# endif

    IMPLICIT NONE

    INTEGER :: I

    if(dbg_set(dbg_sbr)) WRITE(IPT,*) "start RESTORE_STATE"

    IINT   = IINT_BUF
    ExtTime=ExtTime_BUF
    IntTime=IntTime_BUF

    IF(RECALCULATE_RHO_MEAN) THEN
       RECALC_RHO_MEAN =IntTime_BUF
    END IF

    U      = U_BUF
    V      = V_BUF
    WTS    = WTS_BUF
    S1     = S1_BUF
    T1     = T1_BUF

    KM     = KM_BUF
    KH     = KH_BUF
    KQ     = KQ_BUF

    UA     = UA_BUF
    VA     = VA_BUF

    EL     = EL_BUF
    ET     = ET_BUF
    H      = H_BUF

#  if defined (SEDIMENT)
   DO I=1,NSED
     sed(i)%conc = conc_buf(i,:,:)
   END DO
#  endif

#  if defined (EQUI_TIDE)
    EL_EQI = EL_EQI_BUF
#  endif

#  if defined (ATMO_TIDE)
    EL_ATMO= EL_ATMO_BUF
#  endif

# if defined(GOTM)
    TEPS   = TEPS_BUF
    TKE    = TKE_BUF
# else
    Q2     = Q2_BUF
    Q2L    = Q2L_BUF
    L      = L_BUF
# endif

#  if defined (DYE_RELEASE)
    DYE    = DYE_BUF
    DYEF   = DYEF_BUF
!    DYEMEAN = DYEMEAN_BUF
#  endif

    ! TO SAVE MEMEORY WE SAVED ONLY THE INTERNAL VALUES: CALL EXCHANGE
    ! TO SET THE HALO ELEMENTS
    IF(PAR) CALL EXCHANGE_ALL

    ! SET DEPTH AN INTERPOLATE VALUES TO CELL CENTERS
#   if defined (ICE_EMBEDDING)
    D   = H + EL - QTHICE     !!yding
    D   = H + ET - QTHICERK   !!yding
#   else
    D  = H + EL
    DT = H + ET
#   endif


    CALL N2E2D(H,H1)
    CALL N2E2D(EL,EL1)
    CALL N2E2D(D,D1)
    CALL N2E2D(DT,DT1)

    ! SET THE DENSITY
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

    CALL N2E3D(T1,T)
    CALL N2E3D(S1,S)
    !DENSITY IS ALREADY INTERPOLATED TO ELEMENT IN DENSX

    ! SET THE VERTICAL SIGMA VELOCITY
    CALL N2E3D(WTS,W)

    ! SET THE TURBULENT QUANTITIES
#    if defined (GOTM)
    L = .001
    L(1:MT,2:KBM1) = (.5544**3)*TKE(1:MT,2:KBM1)**1.5/TEPS(1:MT,2:KBM1)
# endif

    CALL N2E3D(KM,KM1)


#   if defined(ICE)
    aicen=aicen_buf
    vicen=vicen_buf
    vsnon=vsnon_buf
    tsfcn=tsfcn_buf
    esnon=esnon_buf
    eicen=eicen_buf
    uice2=uice2_buf
    vice2=vice2_buf
    fresh=fresh_buf
    fsalt=fsalt_buf
    fhnet=fhnet_buf
    strength=strength_buf
    isicec=isicec_buf
    sig1=sig1_buf
    sig2=sig2_buf
    sig12=sig12_buf
    CALL  AGGREGATE
# endif

   if(dbg_set(dbg_sbr)) WRITE(IPT,*) "END RESTORE_STATE"

   RETURN
   END SUBROUTINE RESTORE_STATE
!==============================================================================|

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!==============================================================================|
!   EXCHANGE HALO'S BECAUSE WE ONLY SAVED INTERNAL VALUES TO SAVE MEMORY       |
!==============================================================================|
   SUBROUTINE EXCHANGE_ALL

!------------------------------------------------------------------------------|
!!$# if defined (WATER_QUALITY)
!!$     USE MOD_WQM
!!$# endif

#  if defined (DYE_RELEASE)
     USE MOD_DYE
#  endif

     IMPLICIT NONE

!==============================================================================|

   if(dbg_set(dbg_sbr)) WRITE(IPT,*) "001 EXCHANGE_ALL"
#  if defined (MULTIPROCESSOR)

#  if defined (GOTM)
   CALL AEXCHANGE(NC,MYID,NPROCS,TKE,TEPS)
#  else
   CALL AEXCHANGE(NC,MYID,NPROCS,Q2,Q2L,L)
#  endif

   if(dbg_set(dbg_sbr)) WRITE(IPT,*) "002 EXCHANGE_ALL"
   !CALL AEXCHANGE(NC,MYID,NPROCS,KM,KQ,KH)

   CALL AEXCHANGE(NC,MYID,NPROCS,KM)
   if(dbg_set(dbg_sbr)) WRITE(IPT,*) "0021 EXCHANGE_ALL"
   CALL AEXCHANGE(NC,MYID,NPROCS,KQ)
   if(dbg_set(dbg_sbr)) WRITE(IPT,*) "0022 EXCHANGE_ALL"
   CALL AEXCHANGE(NC,MYID,NPROCS,KH)
   if(dbg_set(dbg_sbr)) WRITE(IPT,*) "002a EXCHANGE_ALL"
   CALL AEXCHANGE(EC,MYID,NPROCS,U,V)
   if(dbg_set(dbg_sbr)) WRITE(IPT,*) "002b EXCHANGE_ALL"
   CALL AEXCHANGE(EC,MYID,NPROCS,UA,VA)
   if(dbg_set(dbg_sbr)) WRITE(IPT,*) "002c EXCHANGE_ALL"

   CALL AEXCHANGE(NC,MYID,NPROCS,S1,T1,WTS)
   if(dbg_set(dbg_sbr)) WRITE(IPT,*) "002d EXCHANGE_ALL"
   CALL AEXCHANGE(NC,MYID,NPROCS,EL,H,ET)
   if(dbg_set(dbg_sbr)) WRITE(IPT,*) "002e EXCHANGE_ALL"

   if(dbg_set(dbg_sbr)) WRITE(IPT,*) "003 EXCHANGE_ALL"
#  if defined (EQUI_TIDE)
   CALL AEXCHANGE(NC,MYID,NPROCS,EL_EQI)
#  endif

#  if defined (ATMO_TIDE)
   CALL AEXCHANGE(NC,MYID,NPROCS,EL_ATMO)
#  endif

#  if defined (WATER_QUALITY)
   CALL FATAL_ERROR("SST ASSIMILATION IS NOT SET UP FOR WATER QUALITY MODULE")
#  endif

#  if defined (DYE_RELEASE)
!   CALL AEXCHANGE(NC,MYID,NPROCS,DYEMEAN)
   CALL AEXCHANGE(NC,MYID,NPROCS,DYE)
#  endif

   if(dbg_set(dbg_sbr)) WRITE(IPT,*) "555 EXCHANGE_ALL"
#  endif

   RETURN
   END SUBROUTINE EXCHANGE_ALL

!==============================================================================|

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!==============================================================================|
!   DUMP DATA FILE FOR ASSIMILATION RESTART                                    |
!==============================================================================|

   SUBROUTINE SAVE_STATE

!------------------------------------------------------------------------------|

# if defined (ICE)
   USE mod_ICE
   USE mod_ICE2D, only: sig1, sig2, sig12
# endif

   IMPLICIT NONE

   INTEGER :: I

   ExtTime_BUF = ExtTime
   IntTime_BUF = IntTime
   IINT_BUF   = IINT

   U_BUF      = U
   V_BUF      = V
   WTS_BUF    = WTS
   S1_BUF     = S1
   T1_BUF     = T1

   KM_BUF     = KM
   KH_BUF     = KH
   KQ_BUF     = KQ

# if defined (GOTM)
   TKE_BUF    = TKE
   TEPS_BUF   = TEPS
# else
   Q2_BUF     = Q2
   Q2L_BUF    = Q2L
   L_BUF      = L
# endif

   UA_BUF     = UA
   VA_BUF     = VA

   EL_BUF     = EL
   ET_BUF     = ET
   H_BUF      = H

#  if defined (SEDIMENT)
   do i=1,nsed
      CONC_BUF(i,:,:) = sed(i)%conc
   end do 
#  endif

#  if defined (DYE_RELEASE)
   DYE_BUF    = DYE
   DYEF_BUF   = DYEF
!   DYEMEAN_BUF = DYEMEAN
#  endif

#  if defined (equi_tide)
   EL_EQI_BUF = EL_EQI
#  endif

#  if defined (atmo_tide)
   EL_ATMO_BUF = EL_ATMO
#  endif

#   if defined(ICE)
    aicen_buf=aicen
    vicen_buf=vicen
    vsnon_buf=vsnon
    tsfcn_buf=tsfcn
    esnon_buf=esnon
    eicen_buf=eicen
    uice2_buf=uice2
    vice2_buf=vice2
    fresh_buf=fresh
    fsalt_buf=fsalt
    fhnet_buf=fhnet
    strength_buf=strength
    isicec_buf=isicec
    sig1_buf=sig1
    sig2_buf=sig2
    sig12_buf=sig12
# endif

   IF (TSGRD_ASSIM) THEN
   IF(TSC_SAVE_INDEX2 == 1) THEN

     IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) "DUMPING MONTHLY RSTFILE"
     CALL DUMP_DATA_ASSIM(NC_RST_ASSIM)

     IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) "DUMPING MONTHLY T/S MEAN"
     CALL DUMP_TSAVG_ASSIM(NC_TSAVG_ASSIM)

     TC0_AVG = 0.0_SP
     SC0_AVG = 0.0_SP

   ENDIF
   ENDIF

   IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END SAVE_STATE"

   RETURN
   END SUBROUTINE SAVE_STATE
!==============================================================================|


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   SUBROUTINE ALLOC_BUFFER
!------------------------------------------------------------------------------|
#  if defined (ICE)
   USE ice_model_size
   USE mod_ICE2D, only :sig1,sig2,sig12
#  endif

   IMPLICIT NONE

   IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START SAVE_STATE"

   ALLOCATE(U_BUF(0:NT,KB));   U_BUF = 0.0_SP
   ALLOCATE(V_BUF(0:NT,KB));   V_BUF = 0.0_SP
   ALLOCATE(WTS_BUF(0:MT,KB)); WTS_BUF = 0.0_SP

   ALLOCATE(S1_BUF(0:MT,KB));  S1_BUF = 0.0_SP
   ALLOCATE(T1_BUF(0:MT,KB));  T1_BUF = 0.0_SP

   ALLOCATE(KM_BUF(0:MT,KB));  KM_BUF = 0.0_SP
   ALLOCATE(KH_BUF(0:MT,KB));  KH_BUF = 0.0_SP
   ALLOCATE(KQ_BUF(0:MT,KB));  KQ_BUF = 0.0_SP

   ALLOCATE(UA_BUF(0:NT));     UA_BUF = 0.0_SP
   ALLOCATE(VA_BUF(0:NT));     VA_BUF = 0.0_SP

   ALLOCATE(EL_BUF(0:MT));     EL_BUF = 0.0_SP
   ALLOCATE(ET_BUF(0:MT));     ET_BUF = 0.0_SP
   ALLOCATE(H_BUF(0:MT));      H_BUF = 0.0_SP

#  if defined (SEDIMENT)
   ALLOCATE(CONC_BUF(NSED,0:MT,KB));  CONC_BUF = 0.0_SP
#  endif

#  if defined (EQUI_TIDE)
   ALLOCATE(EL_EQI_BUF(0:MT)); EL_EQI_BUF = 0.0_SP
#  endif

#  if defined (ATMO_TIDE)
   ALLOCATE(EL_ATMO_BUF(0:MT)); EL_ATMO_BUF = 0.0_SP
#  endif

#  if defined (GOTM)
   ALLOCATE(TKE_BUF(0:MT,KB)); TKE_BUF = 0.0_SP
   ALLOCATE(TEPS_BUF(0:MT,KB));TEPS_BUF = 0.0_SP
#  else
   ALLOCATE(Q2_BUF(0:MT,KB));  Q2_BUF = 0.0_SP
   ALLOCATE(Q2L_BUF(0:MT,KB)); Q2L_BUF = 0.0_SP
   ALLOCATE(L_BUF(0:MT,KB));   L_BUF = 0.0_SP
#  endif

#  if defined (DYE_RELEASE)
   ALLOCATE(DYE_BUF(0:MT,KB)); DYE_BUF = 0.0_SP
   ALLOCATE(DYEF_BUF(0:MT,KB)); DYEF_BUF = 0.0_SP
!   ALLOCATE(DYEMEAN_BUF(0:M,KB)); DYEMEAN_BUF = 0.0_SP
#  endif   

#  if defined(ICE)
   ALLOCATE(aicen_buf(m,1,ncat));aicen_buf=0.0_SP
   ALLOCATE(vicen_buf(m,1,ncat));vicen_buf=0.0_SP
   ALLOCATE(vsnon_buf(m,1,ncat));vsnon_buf=0.0_SP
   ALLOCATE(tsfcn_buf(m,1,ncat));tsfcn_buf=0.0_SP
   ALLOCATE(esnon_buf(m,1,ncat));esnon_buf=0.0_SP
   ALLOCATE(eicen_buf(m,1,ntilay));eicen_buf=0.0_SP
   ALLOCATE(uice2_buf(0:NT))    ;uice2_buf=0.0_SP
   ALLOCATE(vice2_buf(0:NT))    ;vice2_buf=0.0_SP
   ALLOCATE(fresh_buf(m,1))     ;fresh_buf=0.0_SP
   ALLOCATE(fsalt_buf(m,1))     ;fsalt_buf=0.0_SP
   ALLOCATE(fhnet_buf(m,1))     ;fhnet_buf=0.0_SP
   ALLOCATE(strength_buf(m,1))  ;strength_buf=0.0_SP
   ALLOCATE(isicec_buf(0:NT))   ;isicec_buf=0.0_SP
   ALLOCATE(sig1_buf(0:mt))     ;sig1_buf=0.0_SP
   ALLOCATE(sig2_buf(0:mt))     ;sig2_buf=0.0_SP
   ALLOCATE(sig12_buf(0:mt))    ;sig12_buf=0.0_SP
#  endif


   RETURN
   END SUBROUTINE ALLOC_BUFFER

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!==============================================================================|
!   DUMP  DATA FILE FOR ASSIMILATION RESTART                                    |
!==============================================================================|

!   SUBROUTINE TS_SAVE_STATE
!------------------------------------------------------------------------------|
!=============================================================
!# if defined (ICE)
!   use mod_ice2d, only :sig1,sig2,sig12
!# endif


!    IMPLICIT NONE
!    INTEGER :: ERROR

!   ExtTime_BUF = ExtTime
!   IntTime_BUF = IntTime
!   IINT_BUF   = IINT

!!# if defined (ICE)
!!   DEALLOCATE(SIG1_BUF,STAT=ERROR)
!!   DEALLOCATE(SIG2_BUF,STAT=ERROR)
!!   DEALLOCATE(SIG12_BUF,STAT=ERROR)
!!   ALLOCATE(SIG1_BUF(0:MT))
!!   ALLOCATE(SIG2_BUF(0:MT))
!!   ALLOCATE(SIG12_BUF(0:MT))
!!   SIG1_BUF=SIG1
!!   SIG2_BUF=SIG2
!!   SIG12_BUF=SIG12
!!# endif

!    !NULLIFY(NC_RST_ASSIM)

!    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START TS_SAVE_STATE"

!!    IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) "SETUP ASSIMULTION RSTFILE"
!!      CALL SETUP_RSTFILE_ASSIM

!    IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) "DUMPING ASSIMULTION RSTFILE"

!      CALL DUMP_DATA_ASSIM(NC_RST_ASSIM)

!    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END TS_SAVE_STATE"

!   END SUBROUTINE TS_SAVE_STATE
!------------------------------------------------------------------------------|

!    SUBROUTINE TS_RESTORE_STATE
!    use MOD_SET_TIME
!    use MOD_STARTUP
!#  if defined (ICE)
!    use mod_ice2d,only : sig1,sig2,sig12
!#  endif
!    USE CONTROL
!    IMPLICIT NONE
!    TYPE(NCFILE), POINTER :: NCF
!    integer :: ncfileind, datfileind,ios,charnum, i
!    logical :: fexist,back,connected
!    character(len=100) :: testchar
!    character(len=160) :: pathnfile
!    character(len=2) :: cios
!    ! CHECK FOR INPUT AND OUTPUT DIRECTORIES
!!  OPEN THE RESTARTE  FILE
!    back = .true.

!    IINT   = IINT_BUF
!    ExtTime=ExtTime_BUF
!    IntTime=IntTime_BUF

!!# if defined (ICE)
!!    SIG1=SIG1_BUF
!!    SIG2=SIG2_BUF
!!    SIG12=SIG12_BUF
!!    deallocate(sig1_buf,sig2_buf,sig12_buf)
!!# endif

!    ! INITIALIZE TYPE TO HOLD FILE METADATA
!    !pathnfile= trim(OUTPUT_DIR)//'TS_assim_state.nc'
!    !pathnfile= trim(INPUT_DIR)//'TS_assim_state.nc'

!    !   IF(.NOT. ASSOCIATED(NC_START)) CALL FATAL_ERROR&
!    !        & ('STARUP FILE IS NOT ASSOCIATED IN SETUP_TIME!')

!    Nullify(NCF)
!    NULLIFY(NC_START)

!    NCF  => NEW_FILE()
!    NCF%FNAME=NC_RST_ASSIM%FNAME !trim(pathnfile)

!    Call NC_OPEN(NCF)
!    CALL NC_LOAD(NCF)
!    NC_START => NCF

!    NCF%FTIME%PREV_STKCNT=1
!    CALL  SET_LAST_FILE_TIME(NC_START,IntTime, IINT)

!!=============================================================
!!  restore the date for the temp restart file
!! =================================================
!#    if defined (DATA_ASSIM)
!     CALL STARTUP_ASSIM
!#    endif

!    RETURN

!    END SUBROUTINE TS_RESTORE_STATE

!------------------------------------------------------------------------------|
!=============================================================
    SUBROUTINE SETUP_RSTFILE_ASSIM
    USE MOD_NCDIO
    IMPLICIT NONE
    TYPE(NCFILE), POINTER :: NCF1, NCF2

    character(len=4) :: cyear

    TYPE(GRID), SAVE :: MYGRID

    NULLIFY(NCF1)

    IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) "START SETUP_ASSIMULTION RSTFILE"

    CALL SET_FVCOM_GRID(MYGRID)

    CALL DEFINE_DIMENSIONS(MYGRID)

    ! ALLOCATE THE NEW FILE OBJECT
    NC_RST_ASSIM => NEW_FILE()

    NC_RST_ASSIM%FTIME => new_ftime()

    !NC_RST_ASSIM%FNAME=trim(OUTPUT_DIR)//'TS_assim_state.nc'
    write(cyear,'(I4.4)')PYear
    NC_RST_ASSIM%FNAME=trim(OUTPUT_DIR)//'TS_assim_state_'//cyear//'.nc'

    IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) "START_ASSIMULTION RSTFILE",trim(OUTPUT_DIR)//'TS_assim_state_'//cyear//'.nc'

    NCF2 => GRID_FILE_OBJECT(MYGRID)
    NC_RST_ASSIM => ADD(NC_RST_ASSIM,NCF2)
!!$    NC_RST_ASSIM => ADD(NC_RST_ASSIM,GRID_FILE_OBJECT(MYGRID) )

    NCF2 => TIME_FILE_OBJECT()
    NC_RST_ASSIM => ADD(NC_RST_ASSIM,NCF2)
!!$    NC_RST_ASSIM => ADD(NC_RST_ASSIM,TIME_FILE_OBJECT() )

    NCF2 => ZETA_FILE_OBJECT()
    NC_RST_ASSIM => ADD(NC_RST_ASSIM,NCF2)
!!$    NC_RST_ASSIM => ADD(NC_RST_ASSIM,ZETA_FILE_OBJECT() )

    NCF2 => FILE_DATE_OBJECT()
    NC_RST_ASSIM => ADD(NC_RST_ASSIM,NCF2)
!!$    NC_RST_ASSIM => ADD(NC_RST_ASSIM,FILE_DATE_OBJECT() )

    NCF2 => VELOCITY_FILE_OBJECT()
    NC_RST_ASSIM => ADD(NC_RST_ASSIM,NCF2)
!!$    NC_RST_ASSIM => ADD(NC_RST_ASSIM,VELOCITY_FILE_OBJECT() )

    NCF2 => AVERAGE_VEL_FILE_OBJECT()
    NC_RST_ASSIM => ADD(NC_RST_ASSIM,NCF2)
!!$    NC_RST_ASSIM => ADD(NC_RST_ASSIM,AVERAGE_VEL_FILE_OBJECT() )

    NCF2 => VERTICAL_VEL_FILE_OBJECT()
    NC_RST_ASSIM => ADD(NC_RST_ASSIM,NCF2)
!!$    NC_RST_ASSIM => ADD(NC_RST_ASSIM,VERTICAL_VEL_FILE_OBJECT() )

    NCF2 => TURBULENCE_FILE_OBJECT()
    NC_RST_ASSIM => ADD(NC_RST_ASSIM,NCF2)
!!$    NC_RST_ASSIM => ADD(NC_RST_ASSIM,TURBULENCE_FILE_OBJECT() )

    NCF2 => SALT_TEMP_FILE_OBJECT()
    NC_RST_ASSIM => ADD(NC_RST_ASSIM,NCF2)
!!$    NC_RST_ASSIM => ADD(NC_RST_ASSIM,SALT_TEMP_FILE_OBJECT() )

    NCF2 => RESTART_EXTRAS_FILE_OBJECT()
    NC_RST_ASSIM => ADD(NC_RST_ASSIM,NCF2)
!!$    NC_RST_ASSIM => ADD(NC_RST_ASSIM,RESTART_EXTRAS_FILE_OBJECT() )

    NCF2 => DENSITY_FILE_OBJECT()
    NC_RST_ASSIM => ADD(NC_RST_ASSIM,NCF2)
!!$    NC_RST_ASSIM => ADD(NC_RST_ASSIM,DENSITY_FILE_OBJECT() )

    IF(WETTING_DRYING_ON) THEN
       NCF2 => WET_DRY_FILE_OBJECT()
       NC_RST_ASSIM => ADD(NC_RST_ASSIM, NCF2)
!!$       NC_RST_ASSIM => ADD(NC_RST_ASSIM, WET_DRY_FILE_OBJECT() )
    END IF

#   if defined (BioGen)
    NCF2 => BIO_FILE_OBJECT()
    NC_RST_ASSIM => ADD(NC_RST_ASSIM,NCF2)
!!$    NC_RST_ASSIM => ADD(NC_RST_ASSIM,BIO_FILE_OBJECT() )
#   endif

#   if defined (WATER_QUALITY)
    NCF2 => WQM_FILE_OBJECT()
    NC_RST_ASSIM => ADD(NC_RST_ASSIM,NCF2)
!!$    NC_RST_ASSIM => ADD(NC_RST_ASSIM,WQM_FILE_OBJECT() )
#   endif
#   if defined (SEDIMENT)
    IF(SEDIMENT_MODEL)THEN
      NCF2 => SED_FILE_OBJECT()
      NC_RST_ASSIM => ADD(NC_RST_ASSIM,NCF2)
!!$      NC_RST_ASSIM => ADD(NC_RST_ASSIM,SED_FILE_OBJECT() )
    ENDIF
#   endif

#   if defined (NH)
    NCF2 => NH_RST_FILE_OBJECT()
    NC_RST_ASSIM => ADD(NC_RST_ASSIM,NCF2)
!!$    NC_RST_ASSIM => ADD(NC_RST_ASSIM,NH_RST_FILE_OBJECT() )
#   endif

# if defined (ICE)
      !-----------------------------------------------------------------
      ! state variables
    NCF2 => ICE_RESTART_STATE_FILE_OBJECT()
    NC_RST_ASSIM => ADD(NC_RST_ASSIM,NCF2)
!!$    NC_RST_ASSIM => ADD(NC_RST_ASSIM,ICE_RESTART_STATE_FILE_OBJECT() )
      !-----------------------------------------------------------------
      ! velocity
      !-----------------------------------------------------------------
    NCF2 => ICE_VEL_FILE_OBJECT()
    NC_RST_ASSIM => ADD(NC_RST_ASSIM,NCF2)
!!$    NC_RST_ASSIM => ADD(NC_RST_ASSIM,ICE_VEL_FILE_OBJECT() )
      !-----------------------------------------------------------------
      ! fresh water, salt, and heat flux
      !-----------------------------------------------------------------
    NCF2 => ICE_EXTRA_FILE_OBJECT()
    NC_RST_ASSIM => ADD(NC_RST_ASSIM,NCF2)
!!$    NC_RST_ASSIM => ADD(NC_RST_ASSIM,ICE_EXTRA_FILE_OBJECT() )
# endif

       CALL NC_WRITE_FILE(NC_RST_ASSIM)
       NC_RST_ASSIM%FTIME%NEXT_STKCNT = 1 !NC_RST_ASSIM%FTIME%NEXT_STKCNT+1

     NC_RST_ASSIM%CONNECTED = .TRUE.

    CALL KILL_DIMENSIONS

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END ASSIMULATION RSTFILE SETUP"


    END SUBROUTINE SETUP_RSTFILE_ASSIM

!------------------------------------------------------------------------------|
!=============================================================
    SUBROUTINE DUMP_DATA_ASSIM(NCF)
    USE MOD_NCDIO
    IMPLICIT NONE
    TYPE(NCFILE), POINTER ::NCF
    TYPE(NCFTIME), POINTER :: FTM

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START DUMP_DATA_ASSIM"


    FTM => NCF%FTIME


    ! IF UPDATE IODATA BECOMES SPECIFIC TO DIFFERENT DATA SETS IT
    ! WILL HAVE TO BE MOVED INSIDE OF THE PARTICULAR OUTPUT STATEMENTS

    IF (FTM%MAX_STKCNT .NE. 0 .AND. &
         & Pmonth==1) THEN   !!! new year first month
         !& FTM%NEXT_STKCNT > FTM%MAX_STKCNT) THEN

       FTM%NEXT_STKCNT=0
       CALL INCRIMENT_FNAME(NCF%FNAME)
       NCF%CONNECTED=.FALSE.


       ! WRITE NEW FILE'S CONSTANT DATA (GRID ETC)
       CALL NC_WRITE_FILE(NCF)

       ! INCRIMENT THE STACK COUNT
       FTM%NEXT_STKCNT = 1

    END IF


    ! IF UPDATE IODATA BECOMES SPECIFIC TO DIFFERENT DATA SETS IT
    ! WILL HAVE TO BE MOVED INSIDE OF THE PARTICULAR OUTPUT STATEMENTS
    CALL UPDATE_IODATA(NCF,IntTime)

    CALL NC_WRITE_FILE(NCF)

    ! ONCE THE FILE IS WRITEN INCRIMENT THE FILE OBJECT TIME
    FTM%PREV_IO = IntTime
    FTM%NEXT_IO = FTM%NEXT_IO + FTM%INTERVAL

    ! INCRIMENT THE STACK COUNT
    FTM%PREV_STKCNT = FTM%NEXT_STKCNT
    FTM%NEXT_STKCNT = FTM%NEXT_STKCNT + 1

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END DUMP_DATA"


    END SUBROUTINE DUMP_DATA_ASSIM
!============================================================= 

!!=============================================================================
!!-----------------------------------------------------------------------------
!!  Setup the netcdf file for output

!    SUBROUTINE MAKE_TSAVG_SAVE
!    USE MOD_NCDIO
!    !USE MOD_FORCE, only : fvcom_cap_grid_source
!    IMPLICIT NONE
!    TYPE(GRID), SAVE :: MYGRID

!    TYPE(NCFILE), POINTER :: NCF
!    TYPE(NCVAR),  POINTER :: VAR
!    TYPE(NCATT),  POINTER :: ATT
!    LOGICAL :: FOUND

!    Call SET_FVCOM_GRID(MYGRID)

!    CALL DEFINE_DIMENSIONS(MYGRID)

!    ! ALLOCATE THE NEW FILE OBJECT
!    NCF => NEW_FILE()


!    NCF%FTIME => new_ftime()

!    !NCF%FNAME = trim(INPUT_DIR)//'TS_save_state.nc'
!    NCF%FNAME = trim(OUTPUT_DIR)//'TS_save_state.nc'

!    NCF => ADD(NCF,GRID_FILE_OBJECT(MYGRID) )

!    NCF => ADD(NCF,TIME_FILE_OBJECT() )

!    NCF => ADD(NCF,ZETA_FILE_OBJECT() )

!    ATT => FIND_ATT(NCF,'source',FOUND)
!    IF(.NOT.FOUND) CALL FATAL_ERROR("LOOKING FOR 'source' ATTRIBUTE: NOT FOUND")

!!    ATT%CHR = fvcom_cap_grid_source

!!    NCF => ADD(NCF,FILE_DATE_OBJECT() )

    ! T_save

!    VAR  => NC_MAKE_AVAR(name='temp_save',&
!    !VAR  => NC_MAKE_AVAR(name='temp',&
!         & values= TC_AVG, DIM1= DIM_node, DIM2= DIM_siglay, DIM3 = DIM_time)

!    ATT  => NC_MAKE_ATT(name='long_name',values='temperature')
!    VAR  => ADD(VAR,ATT)

!    ATT  => NC_MAKE_ATT(name='standard_name',values='sea_water_temperature')
!    VAR  => ADD(VAR,ATT)

!    ATT  => NC_MAKE_ATT(name='units',values='degrees_C')
!    VAR  => ADD(VAR,ATT)

!    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid')
!    VAR  => ADD(VAR,ATT)

!    ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar)
!    VAR  => ADD(VAR,ATT)

!    ATT  => NC_MAKE_ATT(name='type',values='data')

!    VAR  => ADD(VAR,ATT)

!    NCF  => ADD(NCF,VAR)
    ! S_save
!    VAR  => NC_MAKE_AVAR(name='salinity_save',&
!    !VAR  => NC_MAKE_AVAR(name='salinity',&
!         & values= SC_AVG, DIM1= DIM_node, DIM2= DIM_siglay, DIM3 = DIM_time)

!    ATT  => NC_MAKE_ATT(name='long_name',values='salinity')
!    VAR  => ADD(VAR,ATT)

!    ATT  => NC_MAKE_ATT(name='standard_name',values='sea_water_salinity')
!    VAR  => ADD(VAR,ATT)

!    ATT  => NC_MAKE_ATT(name='units',values='1e-3')
!    VAR  => ADD(VAR,ATT)

!    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid')
!    VAR  => ADD(VAR,ATT)

!    ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar)
!    VAR  => ADD(VAR,ATT)

!    ATT  => NC_MAKE_ATT(name='type',values='data')
!    VAR  => ADD(VAR,ATT)
!    NCF  => ADD(NCF,VAR)

!    CALL NC_WRITE_FILE(NCF)
!    NCF%FTIME%NEXT_STKCNT = 1

!    NCF%FTIME%NEXT_STKCNT = 1

!    CALL UPDATE_IODATA(NCF,IntTime)

!    CALL NC_WRITE_FILE(NCF)

!    END SUBROUTINE MAKE_TSAVG_SAVE

!!=============================================================
!!=============================================================


    SUBROUTINE SET_MONTHLY_TS_MEAN(NCF) 
    USE MOD_NCDIO
    !USE MOD_FORCE, only : fvcom_cap_grid_source
    IMPLICIT NONE
    TYPE(GRID), SAVE :: MYGRID

    TYPE(NCFILE), POINTER :: NCF
    TYPE(NCFILE), POINTER :: NCF1
    TYPE(NCFILE), POINTER :: NCF2
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCATT),  POINTER :: ATT
    LOGICAL :: FOUND

    NULLIFY(NCF)
    NULLIFY(NCF1)

    Call SET_FVCOM_GRID(MYGRID)    

    CALL DEFINE_DIMENSIONS(MYGRID)

    ! ALLOCATE THE NEW FILE OBJECT
    NCF => NEW_FILE()


    NCF%FTIME => new_ftime()

    NCF%FNAME = trim(OUTPUT_DIR)//'TSAVG_assim_save.nc'

    NCF2 => GRID_FILE_OBJECT(MYGRID)
    NCF => ADD(NCF,NCF2)
!!$    NCF => ADD(NCF,GRID_FILE_OBJECT(MYGRID) )

    NCF2 => TIME_FILE_OBJECT()
    NCF => ADD(NCF,NCF2)
!!$    NCF => ADD(NCF,TIME_FILE_OBJECT() )
    
    NCF2 => ZETA_FILE_OBJECT()
    NCF => ADD(NCF,NCF2)
!!$    NCF => ADD(NCF,ZETA_FILE_OBJECT() )

    ATT => FIND_ATT(NCF,'source',FOUND)
    IF(.NOT.FOUND) CALL FATAL_ERROR("LOOKING FOR 'source' ATTRIBUTE: NOT FOUND")
    
!    ATT%CHR = fvcom_cap_grid_source

!    NCF => ADD(NCF,FILE_DATE_OBJECT() )

    ! T_save
       
    VAR  => NC_MAKE_AVAR(name='temp_save',&
    !VAR  => NC_MAKE_AVAR(name='temp',&
         & values= TC0_AVG, DIM1= DIM_node, DIM2= DIM_siglay, DIM3 = DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='temperature')
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='standard_name',values='sea_water_temperature')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='degrees_C')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar)
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data')

    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)
    ! S_save
    VAR  => NC_MAKE_AVAR(name='salinity_save',&
    !VAR  => NC_MAKE_AVAR(name='salinity',&
         & values= SC0_AVG, DIM1= DIM_node, DIM2= DIM_siglay, DIM3 = DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='salinity')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='standard_name',values='sea_water_salinity')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='1e-3')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar)
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data')
    VAR  => ADD(VAR,ATT)
    NCF  => ADD(NCF,VAR)

    NCF%FTIME%NEXT_STKCNT = 0
    CALL NC_WRITE_FILE(NCF)
 
    END SUBROUTINE SET_MONTHLY_TS_MEAN

!=============================================================  
    SUBROUTINE DUMP_TSAVG_ASSIM(NCF)
    USE MOD_NCDIO
    IMPLICIT NONE  
    TYPE(NCFILE), POINTER ::NCF
    TYPE(NCFTIME), POINTER :: FTM

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START DUMP_TSAVG_ASSIM"

    FTM => NCF%FTIME
    FTM%NEXT_STKCNT = FTM%NEXT_STKCNT + 1

    ! IF UPDATE IODATA BECOMES SPECIFIC TO DIFFERENT DATA SETS IT
    ! WILL HAVE TO BE MOVED INSIDE OF THE PARTICULAR OUTPUT STATEMENTS
    CALL UPDATE_IODATA(NCF,IntTime)
    CALL NC_WRITE_FILE(NCF)

    ! ONCE THE FILE IS WRITEN INCRIMENT THE FILE OBJECT TIME
    
    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END DUMP_TSAVG_ASSIM"
    
    END SUBROUTINE DUMP_TSAVG_ASSIM
!=============================================================
!==============================================================================|
   LOGICAL FUNCTION ISINTRIANGLE1(XT,YT,X0,Y0) 
!==============================================================================|
!  determine if point (x0,y0) is in triangle defined by nodes (xt(3),yt(3))    |
!  using algorithm used for scene rendering in computer graphics               |
!  algorithm works well unless particle happens to lie in a line parallel      |
!  to the edge of a triangle.                                                  |
!  This can cause problems if you use a regular grid, say for idealized        |
!  modelling and you happen to see particles right on edges or parallel to     |
!  edges.                                                                      |
!==============================================================================|

   USE MOD_PREC
   IMPLICIT NONE
   REAL(SP), INTENT(IN) :: X0,Y0
   REAL(SP), INTENT(IN) :: XT(3),YT(3)
   REAL(SP) :: F1,F2,F3
   REAL(SP) :: X1(2)
   REAL(SP) :: X2(2)
   REAL(SP) :: X3(2)
   REAL(SP) :: P(2)

!------------------------------------------------------------------------------|
   IF (DBG_SET(DBG_SBR)) WRITE(IPT,*) "START: ISINTRIANGLE1" 

   ISINTRIANGLE1 = .FALSE.  

   IF(Y0 < MINVAL(YT) .OR. Y0 > MAXVAL(YT))THEN
     ISINTRIANGLE1 = .FALSE.
     RETURN
   END IF
   IF(X0 < MINVAL(XT) .OR. X0 > MAXVAL(XT))THEN
     ISINTRIANGLE1 = .FALSE.
     RETURN
   END IF

   F1 = (Y0-YT(1))*(XT(2)-XT(1)) - (X0-XT(1))*(YT(2)-YT(1))
   F2 = (Y0-YT(3))*(XT(1)-XT(3)) - (X0-XT(3))*(YT(1)-YT(3))
   F3 = (Y0-YT(2))*(XT(3)-XT(2)) - (X0-XT(2))*(YT(3)-YT(2))
   IF(F1*F3 >= 0.0_SP .AND. F3*F2 >= 0.0_SP) ISINTRIANGLE1 = .TRUE.

   IF (DBG_SET(DBG_SBR)) WRITE(IPT,*) "END: ISINTRIANGLE1" 
   RETURN
   END FUNCTION ISINTRIANGLE1
!==============================================================================|
!
!==============================================================================|
  SUBROUTINE NAME_LIST_PRINT_ASSIM
    USE CONTROL

    IMPLICIT NONE

    write(UNIT=IPT,NML=NML_SST_ASSIMILATION)
    write(UNIT=IPT,NML=NML_SSTGRD_ASSIMILATION)
    write(UNIT=IPT,NML=NML_SSHGRD_ASSIMILATION)
    write(UNIT=IPT,NML=NML_TSGRD_ASSIMILATION)
    write(UNIT=IPT,NML=NML_CUR_NGASSIMILATION)
    write(UNIT=IPT,NML=NML_CUR_OIASSIMILATION)
    write(UNIT=IPT,NML=NML_TS_NGASSIMILATION)
    write(UNIT=IPT,NML=NML_TS_OIASSIMILATION)

  RETURN
  END SUBROUTINE NAME_LIST_PRINT_ASSIM

# endif
END MODULE MOD_ASSIM
