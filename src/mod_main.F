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

!==============================================================================|
!   GLOBAL LIMITS AND ARRAY SIZING PARAMETERS                                  !
!==============================================================================|

MODULE LIMS
   USE MOD_PREC  
   IMPLICIT NONE
   SAVE
 
   INTEGER NGL                !!GLOBAL NUMBER OF ELEMENTS
   INTEGER MGL                !!GLOBAL NUMBER OF NODES
   INTEGER NUMQBC_GL          !!GLOBAL NUMBER OF FRESHWATER INFLOW NODES (RIVERS)
   INTEGER NOBCLSF_GL         !!GLOBAL NUMBER OF LONGSHORE FLOW ADJUSTED OPEN BOUNDARY NODES
   INTEGER NDRFT_GL           !!GLOBAL NUMBER OF LAGRANGIAN TRACKING PARTICLES 

   INTEGER N                  !!LOCAL NUMBER OF ELEMENTS 
   INTEGER M                  !!LOCAL NUMBER OF NODES
   INTEGER NUMQBC             !!LOCAL NUMBER OF FRESHWATER INFLOW NODES
   INTEGER NOBCLSF            !!LOCAL NUMBER OF LONGSHORE FLOW ADJUSTED OPEN BOUNDARY NODES
   INTEGER NDRFT              !!LOCAL NUMBER OF LAGRANGIAN TRACKING PARTICLES
   INTEGER NISBCE_1           !!LOCAL NUMBER OF ELEMENTS WITH ISBCE = 1
   INTEGER NISBCE_2           !!LOCAL NUMBER OF ELEMENTS WITH ISBCE = 2
   INTEGER NISBCE_3           !!LOCAL NUMBER OF ELEMENTS WITH ISBCE = 3

   INTEGER KB                 !!NUMBER OF SIGMA LEVELS
   INTEGER KBM1               !!NUMBER OF SIGMA LEVELS-1
   INTEGER KBM2               !!NUMBER OF SIGMA LEVELS-2
   INTEGER MYID               !!UNIQUE PROCESSOR ID (1 => NPROCS)
   INTEGER MSRID
   INTEGER IOPROCID           !!PROCESSOR ID OF THE PRINTER NODE IF PRESENT
!   INTEGER KSL               !!NUMBER OF STANDARD SEA LEVELS 
   INTEGER,POINTER:: NPROCS_TOTAL !!TOTAL NUMBER OF PROCESSORS
   INTEGER,TARGET :: NPROCS       !!NUMBER OF PROCESSORS IN FVCOM
   INTEGER NE                 !!NUMBER OF UNIQUE EDGES (LOCAL DOMAIN ONLY)
   INTEGER NCV                !!NUMBER OF INTERNAL CONTROL VOLUMES (EXTENDED LOCAL ONLY)
   
   INTEGER NCV_I              !!NUMBER OF INTERNAL CONTROL VOLUMES (LOCAL ONLY)
   INTEGER NT                 !!TOTAL OF LOCAL INTERNAL + HALO ELEMENTS
   INTEGER MT                 !!TOTAL OF LOCAL INTERNAL + HALO NODES
   INTEGER MX_NBR_ELEM        !!MAX NUMBER OF ELEMENTS SURROUNDING A NODE

   REAL(SP) :: MEMTOT,MEMCNT


END MODULE LIMS



!==============================================================================|
!   CONTROL VARIABLES                                                          |
!==============================================================================|

MODULE CONTROL
   USE MOD_TIME
   USE MOD_PREC   
   USE NETCDF
   IMPLICIT NONE
   SAVE

   ! Stuff set inside FVCOM at compile time:
   LOGICAL :: SERIAL         !!TRUE IF SINGLE PROCESSOR
   LOGICAL ::    MSR         !!TRUE IF MASTER PROCESSOR (MYID==1)
   LOGICAL ::    PAR         !!TRUE IF MULTIPROCESSOR RUN
   LOGICAL ::   IOPROC         !!TRUE IF PROCESSOR IS THE IOPROC
   LOGICAL :: IN_MPI_IO_LOOP !!TRUE IF OUTNODE IS ALREADY IN THE LOOP
   CHARACTER(LEN=80) PRG_NAME

   INTEGER MPI_FVCOM_GROUP         !! PROCESSORS WORKING ON FVCOM
!   INTEGER MPI_OTHER_GROUP        !! ADD OTHER GROUPS HERE 

   ! FVCOM info strings
   CHARACTER(LEN=80) FVCOM_VERSION !!STRING DESCRIBING VERSION
   CHARACTER(LEN=80) INSTITUTION   !!STRING DESCRIBING FVCOM MOTHER SHIP
   CHARACTER(LEN=80) FVCOM_WEBSITE !!STRING DESCRIBING WEBSITE FOR FVCOM INFO 

   ! Set at command line
   CHARACTER(LEN=80) CASENAME      !!LETTER ACRONYM SPECIFYING CASE
   CHARACTER(LEN=80) INFOFILE      !!INFO  FILE            

   ! FOR LOADING RUNFILE
   LOGICAL BLANK_NAMELIST
   CHARACTER(LEN=80) NAMELIST_NAME

    !--Parameters in NameList NML_CASE
    CHARACTER(LEN=80) CASE_TITLE     !!CASE TITLE                                 
    CHARACTER(Len=80) DATE_FORMAT    !!DATE FORMAT: ie 2007-11-19 => YMD 
    CHARACTER(Len=80) TIMEZONE       !!Name of TIME ZONE OR NONE
    CHARACTER(Len=80) START_DATE     !!DATE TO START MODEL
    CHARACTER(Len=80) END_DATE       !!DATE TO END MODEL 
    CHARACTER(Len=80) DATE_REFERENCE !!USER DEFINED REFERENCE DATE OR DEFAULT
    
    LOGICAL USE_REAL_WORLD_TIME
    NAMELIST /NML_CASE/        &
         & CASE_TITLE,         &
         & TIMEZONE,           &
         & DATE_FORMAT,        &
	 & DATE_REFERENCE,     &
         & START_DATE,         &
         & END_DATE

    !--Paramters in NameList NML_STARTUP
    CHARACTER(LEN=80) STARTUP_TYPE  !!'hotstart' or 'coldstart' 
    CHARACTER(LEN=80) STARTUP_FILE  !!NAME OF START FILE
    CHARACTER(LEN=80) STARTUP_TS_TYPE  !!TYPE OF TS START
    CHARACTER(LEN=80) STARTUP_UV_TYPE  !!TYPE OF TS START
    CHARACTER(LEN=80) STARTUP_TURB_TYPE  !!TYPE OF TS START
    REAL(SP) :: STARTUP_T_VALS(2)
    REAL(SP) :: STARTUP_S_VALS(2)
    REAL(SP) :: STARTUP_U_VALS
    REAL(SP) :: STARTUP_V_VALS
    REAL(SP) :: STARTUP_DMAX
        
    LOGICAL CMDLN_RESTART

    CHARACTER(LEN=80), parameter :: STARTUP_TYPE_COLDSTART = 'coldstart'
    CHARACTER(LEN=80), parameter :: STARTUP_TYPE_HOTSTART = 'hotstart'
    CHARACTER(LEN=80), parameter :: STARTUP_TYPE_CRASHRESTART = 'crashrestart'
    CHARACTER(LEN=80), parameter :: STARTUP_TYPE_FORECAST = 'forecast'

    LOGICAL FORECAST_MODE

    CHARACTER(LEN=80), parameter :: STARTUP_TYPE_DEFAULT = 'default'
    CHARACTER(LEN=80), parameter :: STARTUP_TYPE_CONSTANT = 'constant'
    CHARACTER(LEN=80), parameter :: STARTUP_TYPE_LINEAR = 'linear'
    CHARACTER(LEN=80), parameter :: STARTUP_TYPE_OBSERVED = 'observed'
    CHARACTER(LEN=80), parameter :: STARTUP_TYPE_SETVALUES = 'set values'


    NAMELIST /NML_STARTUP/     &
         & STARTUP_TYPE,       & ! coldstart, hotstart, crashrestart
         & STARTUP_FILE,       &
         & STARTUP_UV_TYPE,    &  ! default, set values
         & STARTUP_TURB_TYPE,  &  ! default, set values
         & STARTUP_TS_TYPE,    &  ! constant, linear, observed, set values
         & STARTUP_T_VALS,     &
         & STARTUP_S_VALS,     &
         & STARTUP_U_VALS,     &
         & STARTUP_V_VALS,     &
         & STARTUP_DMAX



    !--Parameters in NameList NML_IO
    CHARACTER(LEN=80) INPUT_DIR         !!MAIN   INPUT DIRECTORY
    CHARACTER(LEN=80) OUTPUT_DIR        !!PARENT OUTPUT DIRECTORY
    INTEGER IREPORT                !!INTERVAL (IINT) FOR REPORTING OF FLOWFIELD STATISTICS
    LOGICAL VISIT_ALL_VARS         !!SET THE LEVEL OF COMPLEXITY IN VISIT
    LOGICAL WAIT_FOR_VISIT         !!WAIT FOR VISIT TO CONNECT BEFORE BEGINING INTEGRATION
    LOGICAL USE_MPI_IO_MODE        !!TURN ON THE PRINTER NODE FOR MPI JOBS

    NAMELIST /NML_IO/         &
         & INPUT_DIR,         &
         & OUTPUT_DIR,        &
         & IREPORT,           &
         & VISIT_ALL_VARS,    &
         & WAIT_FOR_VISIT,    &
         & USE_MPI_IO_MODE
        
    !--Parameters in NameList NML_INTEGRATION

! FOR EXPLICIT TIME STEP MODEL
    REAL(DP)  EXTSTEP_SECONDS!!EXTERNAL TIME STEP IN SECONDS;
                             !! THIS VALUE WILL BE TRUNCATED TO MICROSECONDS!!
    INTEGER   ISPLIT         !!NUMBER OF ITERATIONS OF EXTERNAL MODE/INTERNAL STEP

! FOR SEMI-IMPLICIT TIME STEP MODEL
    REAL(DP)  INTSTEP_SECONDS!!INTERNAL TIME STEP IN SECONDS FOR SEMI IMPLICIT;

    INTEGER   IRAMP          !!RAMP FACTOR USED TO EASE STARTUP = f(IINT)
    REAL(SP)  STATIC_SSH_ADJ !!WATER LEVEL CONSTANT ADJUSTMENT
    REAL(SP)  MIN_DEPTH      !!MINIMUM ALLOWABLE DEPTH         

    NAMELIST /NML_INTEGRATION/  &
# if defined(SEMI_IMPLICIT)
         & INTSTEP_SECONDS,     &
# else
         & EXTSTEP_SECONDS,     &
         & ISPLIT,              &
# endif
         & IRAMP,               &
         & MIN_DEPTH,           &
         & STATIC_SSH_ADJ
    


    !--Parameters in NameList NML_RESTART
    LOGICAL RST_ON                  !!TRUE IF OUTPUT RESART FILES
    CHARACTER(LEN=80) RST_FIRST_OUT !!DATE OF FIRST RESTART FILE OUTPUT
    CHARACTER(LEN=80) RST_OUT_INTERVAL       !!INTERVAL IN DECIMAL DAYS FOR CREATING A RESTART FILE
    INTEGER  RST_OUTPUT_STACK       !!NUMBER OF TIMESTEPS PER FILE

    CHARACTER(LEN=80) RESTART_FILE_NAME !!NAME OF RESTART FILE
    ! ONLY VALID TILL THE NAME IS INCRIMENTED IF THE STKCNT IS EXCEEDED

    NAMELIST /NML_RESTART/   &
         & RST_ON,           &
         & RST_FIRST_OUT,    &
         & RST_OUT_INTERVAL, &
         & RST_OUTPUT_STACK

        
    !--Parameters in NameList NML_NETCDF
    LOGICAL NC_ON
    CHARACTER(LEN=80) NC_FIRST_OUT        !! DATE TO START NETCDF OUTPUT
    CHARACTER(LEN=80) NC_OUT_INTERVAL     !! OUTPUT INTERVAL IN DECIMAL DAYS
    INTEGER  NC_OUTPUT_STACK              !! NUMBER OF TIMESTEPS PER FILE
    CHARACTER(LEN=120) NC_SUBDOMAIN_FILES !! DOMAIN OUTPUT SPECS
    ! GRID STUFF
    LOGICAL NC_GRID_METRICS
    LOGICAL NC_FILE_DATE
    ! MODEL DATA
    LOGICAL NC_VELOCITY
    LOGICAL NC_SALT_TEMP
    LOGICAL NC_TURBULENCE
    LOGICAL NC_VERTICAL_VEL
    LOGICAL NC_AVERAGE_VEL
    LOGICAL NC_VORTICITY
    LOGICAL NC_NH_QP
    LOGICAL NC_NH_RHS
    LOGICAL NC_ICE
    
    ! FORCING DATA
    LOGICAL NC_WIND_VEL
    LOGICAL NC_WIND_STRESS
#   if defined (AIR_PRESSURE) || (HEATING_CALCULATED)
    LOGICAL NC_ATM_PRESS  
#   endif    
    LOGICAL NC_WAVE_PARA     
    LOGICAL NC_WAVE_STRESS   
    LOGICAL NC_EVAP_PRECIP
    LOGICAL NC_SURFACE_HEAT
    LOGICAL NC_GROUNDWATER
    LOGICAL NC_BIO
    LOGICAL NC_WQM

    CHARACTER(LEN=80) NC_FILE_NAME !!NAME OF NC FILE
    ! ONLY VALID TILL THE NAME IS INCRIMENTED IF THE STKCNT IS EXCEEDED

    NAMELIST /NML_NETCDF/    &
         & NC_ON,            &
         & NC_FIRST_OUT,     &
         & NC_OUT_INTERVAL,  &
         & NC_OUTPUT_STACK,  &
         & NC_SUBDOMAIN_FILES,&
         & NC_GRID_METRICS,  &
         & NC_FILE_DATE,     &
         & NC_VELOCITY,      &
         & NC_SALT_TEMP,     &
         & NC_TURBULENCE,    &
         & NC_AVERAGE_VEL,   &
         & NC_VERTICAL_VEL,  &
# if defined(NH)
         & NC_NH_QP,         &
         & NC_NH_RHS,        &
# endif
         & NC_WIND_VEL,      &
         & NC_WIND_STRESS,   &
# if defined (AIR_PRESSURE) || (HEATING_CALCULATED)
         & NC_ATM_PRESS,     &
# endif
# if defined(WAVE_CURRENT_INTERACTION)
         & NC_WAVE_PARA,     &    
         & NC_WAVE_STRESS,   &
# endif
         & NC_EVAP_PRECIP,   &
         & NC_SURFACE_HEAT,  &
         & NC_GROUNDWATER,   &
# if defined(ICE)
         & NC_ICE,           &
# endif
         & NC_BIO,           &
         & NC_WQM,           &
	 & NC_VORTICITY

    !--Parameters in NameList NML_NETCDF_SURFACE
    LOGICAL NCSF_ON
    CHARACTER(LEN=80) NCSF_FIRST_OUT        !! DATE TO START NETCDF OUTPUT
    CHARACTER(LEN=80) NCSF_OUT_INTERVAL     !! OUTPUT INTERVAL IN DECIMAL DAYS
    INTEGER  NCSF_OUTPUT_STACK              !! NUMBER OF TIMESTEPS PER FILE
    CHARACTER(LEN=120) NCSF_SUBDOMAIN_FILES !! DOMAIN OUTPUT SPECS
    ! GRID STUFF
    LOGICAL NCSF_GRID_METRICS
    LOGICAL NCSF_FILE_DATE
    ! MODEL DATA
    LOGICAL NCSF_VELOCITY
    LOGICAL NCSF_SALT_TEMP
    LOGICAL NCSF_TURBULENCE
    LOGICAL NCSF_ICE
    LOGICAL NCSF_WAVE_PARA     
    
    ! FORCING DATA
    LOGICAL NCSF_WIND_VEL
    LOGICAL NCSF_WIND_STRESS
#   if defined (AIR_PRESSURE) || (HEATING_CALCULATED)
    LOGICAL NCSF_ATM_PRESS  
#   endif    
    LOGICAL NCSF_EVAP_PRECIP
    LOGICAL NCSF_SURFACE_HEAT

    CHARACTER(LEN=80) NCSF_FILE_NAME !!NAME OF NC FILE
    ! ONLY VALID TILL THE NAME IS INCRIMENTED IF THE STKCNT IS EXCEEDED

    NAMELIST /NML_NETCDF_SURFACE/    &
         & NCSF_ON,            &
         & NCSF_FIRST_OUT,     &
         & NCSF_OUT_INTERVAL,  &
         & NCSF_OUTPUT_STACK,  &
         & NCSF_SUBDOMAIN_FILES,&
         & NCSF_GRID_METRICS,  &
         & NCSF_FILE_DATE,     &
         & NCSF_VELOCITY,      &
         & NCSF_SALT_TEMP,     &
         & NCSF_TURBULENCE,    &
         & NCSF_WIND_VEL,      &
         & NCSF_WIND_STRESS,   &
# if defined (AIR_PRESSURE) || (HEATING_CALCULATED)
         & NCSF_ATM_PRESS,     &
# endif
# if defined(WAVE_CURRENT_INTERACTION)
         & NCSF_WAVE_PARA,     &    
# endif
# if defined(ICE)
         & NCSF_ICE,           &
# endif
         & NCSF_EVAP_PRECIP,   &
         & NCSF_SURFACE_HEAT

    !--Parameters in NameList NML_NETCDF_AV
    LOGICAL NCAV_ON                         !! TURN ON NETCDF AVERAGING
    CHARACTER(LEN=80) NCAV_FIRST_OUT        !! DATE TO START NETCDF AVERAGE OUTPUT
    CHARACTER(LEN=80) NCAV_OUT_INTERVAL     !! OUTPUT INTERVAL IN DECIMAL DAYS
    INTEGER  NCAV_OUTPUT_STACK              !! NUMBER OF TIMESTEPS PER FILE
    CHARACTER(LEN=120) NCAV_SUBDOMAIN_FILES !! DOMAIN OUTPUT SPECS
    ! GRID STUFF
    LOGICAL NCAV_GRID_METRICS
    LOGICAL NCAV_FILE_DATE
    ! MODEL DATA
    LOGICAL NCAV_VELOCITY
    LOGICAL NCAV_SALT_TEMP
    LOGICAL NCAV_TURBULENCE
    LOGICAL NCAV_VERTICAL_VEL
    LOGICAL NCAV_AVERAGE_VEL
    LOGICAL NCAV_VORTICITY
    LOGICAL NCAV_NH_QP
    LOGICAL NCAV_NH_RHS
    LOGICAL NCAV_ICE

    ! FORCING DATA
    LOGICAL NCAV_WIND_VEL
    LOGICAL NCAV_WIND_STRESS
# if defined (AIR_PRESSURE) || (HEATING_CALCULATED)
    LOGICAL NCAV_ATM_PRESS  
# endif
    LOGICAL NCAV_WAVE_PARA    
    LOGICAL NCAV_WAVE_STRESS  
    LOGICAL NCAV_EVAP_PRECIP
    LOGICAL NCAV_SURFACE_HEAT
    LOGICAL NCAV_GROUNDWATER
    LOGICAL NCAV_BIO
    LOGICAL NCAV_WQM

    CHARACTER(LEN=80) NCAV_FILE_NAME !!NAME OF NCAV FILE
    ! ONLY VALID TILL THE NAME IS INCRIMENTED IF THE STKCNT IS EXCEEDED

    NAMELIST /NML_NETCDF_AV/  &
         & NCAV_ON,           &
         & NCAV_FIRST_OUT,    &
         & NCAV_OUT_INTERVAL, &
         & NCAV_OUTPUT_STACK, &
         & NCAV_SUBDOMAIN_FILES,&
         & NCAV_GRID_METRICS,   &
         & NCAV_FILE_DATE,      &
         & NCAV_VELOCITY,       &
         & NCAV_SALT_TEMP,      &
         & NCAV_TURBULENCE,     &
         & NCAV_AVERAGE_VEL,    &
         & NCAV_VERTICAL_VEL,   &
# if defined(NH)
         & NCAV_NH_QP,          &
         & NCAV_NH_RHS,         &
# endif
         
# if defined(ICE)
         & NCAV_ICE,          &
# endif
         & NCAV_WIND_VEL,       &
         & NCAV_WIND_STRESS,    &
# if defined (AIR_PRESSURE) || (HEATING_CALCULATED)
         & NCAV_ATM_PRESS,      &
# endif
# if defined(WAVE_CURRENT_INTERACTION)
         & NCAV_WAVE_PARA,      &    
         & NCAV_WAVE_STRESS,    &
# endif
         & NCAV_EVAP_PRECIP,    &
         & NCAV_SURFACE_HEAT,   &
         & NCAV_GROUNDWATER,    &
         & NCAV_BIO,            &
         & NCAV_WQM,            &
	 & NCAV_VORTICITY

    !--Parameters in NameList NML_PHYSICS
    CHARACTER(LEN=80) HORIZONTAL_MIXING_TYPE
    CHARACTER(LEN=80) HORIZONTAL_MIXING_FILE
    CHARACTER(LEN=80) HORIZONTAL_MIXING_KIND
    REAL(SP)  HORIZONTAL_MIXING_COEFFICIENT
    REAL(SP)  HORIZONTAL_PRANDTL_NUMBER

! REMOVED HORCON, THIS IS NOW STORED IN THE ARRAYS, NN_HVC AND CC_HVC
!    REAL(SP) HORCON ! FVCOM NAME
    REAL(SP) HPRNU  ! FVCOM NAME

    CHARACTER(LEN=80) VERTICAL_MIXING_TYPE
    REAL(SP) VERTICAL_MIXING_COEFFICIENT
    REAL(SP) VERTICAL_PRANDTL_NUMBER

    REAL(SP) UMOL   ! FVCOM NAME
    REAL(SP) VPRNU  ! FVCOM NAME

    CHARACTER(LEN=80) BOTTOM_ROUGHNESS_KIND 
    CHARACTER(LEN=80) BOTTOM_ROUGHNESS_TYPE 
    CHARACTER(LEN=80) BOTTOM_ROUGHNESS_FILE 
    REAL(SP)  BOTTOM_ROUGHNESS_MINIMUM          !!MINIMUM BOTTOM DRAG COEFFICIENT
    REAL(SP)  BOTTOM_ROUGHNESS_LENGTHSCALE      !!BOTTOM FRICTION DEPTH LENGTH SCALE

    REAL(SP) CBCMIN ! FVCOM NAME

    CHARACTER(LEN=80),parameter :: BR_ORIG   = 'orig'
    CHARACTER(LEN=80),parameter :: BR_GOTM   = 'gotm'

    LOGICAL CONVECTIVE_OVERTURNING
    LOGICAL SCALAR_POSITIVITY_CONTROL
    LOGICAL BAROTROPIC
    CHARACTER(LEN=80) SEA_WATER_DENSITY_FUNCTION
    CHARACTER(LEN=80) BAROCLINIC_PRESSURE_GRADIENT
    LOGICAL TEMPERATURE_ACTIVE
    LOGICAL SALINITY_ACTIVE
!J. Ge
    ! for tracer advection
    LOGICAL BACKWARD_ADVECTION
    INTEGER BACKWARD_STEP
!J. Ge
    LOGICAL SURFACE_WAVE_MIXING

    ! AT PRESENT YOU CAN NOT COMPILE WITH WET DRY WITHOUT USING IT!
    LOGICAL WETTING_DRYING_ON

    LOGICAL  RECALCULATE_RHO_MEAN
    CHARACTER(LEN=80) INTERVAL_RHO_MEAN

    CHARACTER(LEN=80),parameter :: SW_DENS1 = 'dens1'
    CHARACTER(LEN=80),parameter :: SW_DENS2 = 'dens2'
    CHARACTER(LEN=80),parameter :: SW_DENS3 = 'dens3'

    LOGICAL ADCOR_ON
    LOGICAL EQUATOR_BETA_PLANE
    LOGICAL NOFLUX_BOT_CONDITION

    NAMELIST /NML_PHYSICS/               &
         & HORIZONTAL_MIXING_TYPE,       &
         & HORIZONTAL_MIXING_FILE,       &
         & HORIZONTAL_MIXING_KIND,       &
         & HORIZONTAL_MIXING_COEFFICIENT,&
         & HORIZONTAL_PRANDTL_NUMBER,    & 
         & VERTICAL_MIXING_TYPE,         &
         & VERTICAL_MIXING_COEFFICIENT,  &
         & VERTICAL_PRANDTL_NUMBER,      &
         & BOTTOM_ROUGHNESS_TYPE,        &
         & BOTTOM_ROUGHNESS_KIND,        &
         & BOTTOM_ROUGHNESS_FILE,        &
         & BOTTOM_ROUGHNESS_LENGTHSCALE, &
         & BOTTOM_ROUGHNESS_MINIMUM,     &
         & CONVECTIVE_OVERTURNING,       &
         & SCALAR_POSITIVITY_CONTROL,    &
         & BAROTROPIC,                   &
         & BAROCLINIC_PRESSURE_GRADIENT, &
         & SEA_WATER_DENSITY_FUNCTION,   &
         & RECALCULATE_RHO_MEAN,         &
         & INTERVAL_RHO_MEAN,            &
         & TEMPERATURE_ACTIVE,           &
         & SALINITY_ACTIVE,              &
         & SURFACE_WAVE_MIXING,          &
!J. Ge
         & BACKWARD_ADVECTION,           &
         & BACKWARD_STEP,                &
!J. Ge
         & WETTING_DRYING_ON,            &
         & ADCOR_ON,                     &
	 & EQUATOR_BETA_PLANE,           &
         & NOFLUX_BOT_CONDITION


    !--Parameters in NameList NML_SURFACE_FORCING
  
    LOGICAL  WIND_ON
    CHARACTER(LEN=80) WIND_TYPE
    CHARACTER(LEN=80) WIND_FILE
    CHARACTER(LEN=80) WIND_KIND
    REAL(SP) WIND_X
    REAL(SP) WIND_Y
    
    LOGICAL HEATING_ON
    CHARACTER(LEN=80) HEATING_TYPE
    CHARACTER(LEN=80) HEATING_FILE
    CHARACTER(LEN=80) HEATING_KIND
    REAL(SP) HEATING_LONGWAVE_PERCTAGE
    REAL(SP) HEATING_LONGWAVE_LENGTHSCALE
    REAL(SP) HEATING_SHORTWAVE_LENGTHSCALE
    REAL(SP) HEATING_RADIATION
    REAL(SP) HEATING_NETFLUX

    !FVCOM NAMES
    REAL(SP) RHEAT ! LONG WAVE PERCENTAGE
    REAL(SP) ZETA1 ! LONG WAVE LENGTH
    REAL(SP) ZETA2 ! SHORT WAVE LENGTH
    

    LOGICAL PRECIPITATION_ON
    CHARACTER(LEN=80) PRECIPITATION_FILE    
    CHARACTER(LEN=80) PRECIPITATION_KIND    
    REAL(SP) :: PRECIPITATION_EVP
    REAL(SP) :: PRECIPITATION_PRC

    LOGICAL AIRPRESSURE_ON
    CHARACTER(LEN=80) AIRPRESSURE_FILE    
    CHARACTER(LEN=80) AIRPRESSURE_KIND
    REAL(SP) :: AIRPRESSURE_VALUE

    LOGICAL WAVE_ON               
    CHARACTER(LEN=80) WAVE_FILE    
    CHARACTER(LEN=80) WAVE_KIND
    REAL(SP) :: WAVE_HEIGHT
    REAL(SP) :: WAVE_LENGTH
    REAL(SP) :: WAVE_DIRECTION
    REAL(SP) :: WAVE_PERIOD
    REAL(SP) :: WAVE_PER_BOT
    REAL(SP) :: WAVE_UB_BOT

    ! FORCING KINDS
    CHARACTER(LEN=80),parameter:: CNSTNT   = "constant"
    CHARACTER(LEN=80),parameter:: STTC     = "static"
    CHARACTER(LEN=80),parameter:: TMDPNDNT = "time dependant"
    CHARACTER(LEN=80),parameter:: PRDC     = "periodic"
    CHARACTER(LEN=80),parameter:: VRBL     = "variable"

    CHARACTER(LEN=80), parameter:: SPEED  = "speed"
    CHARACTER(LEN=80), parameter:: STRESS = "stress"

!J. Ge
    CHARACTER(LEN=80),parameter:: UNIFORM = "uniform"
    CHARACTER(LEN=80),parameter:: NON_UNIFORM = "non-uniform"
!J. Ge

    NAMELIST /NML_SURFACE_FORCING/          &
         & WIND_ON,                         &
         & WIND_TYPE,                       &
         & WIND_FILE,                       &
         & WIND_KIND,                       &
         & WIND_X,                          &
         & WIND_Y,                          &
         & HEATING_ON,                      &
         & HEATING_TYPE,                    &
         & HEATING_KIND,                    &
         & HEATING_FILE,                    &
         & HEATING_LONGWAVE_LENGTHSCALE,    &
         & HEATING_LONGWAVE_PERCTAGE,       &
         & HEATING_SHORTWAVE_LENGTHSCALE,   & 
         & HEATING_RADIATION,               &
         & HEATING_NETFLUX,                 & 
         & PRECIPITATION_ON,                &
         & PRECIPITATION_KIND,              &
         & PRECIPITATION_FILE,              &
         & PRECIPITATION_PRC,               &
         & PRECIPITATION_EVP,               &
	 & AIRPRESSURE_ON,                  &
	 & AIRPRESSURE_KIND,                &
	 & AIRPRESSURE_FILE,                &
	 & AIRPRESSURE_VALUE,               &
         & WAVE_ON,                         &   
         & WAVE_FILE,                       &
         & WAVE_KIND,                       &
         & WAVE_HEIGHT,                     &
         & WAVE_LENGTH,                     &
         & WAVE_DIRECTION,                  &
         & WAVE_PERIOD,                     &
         & WAVE_PER_BOT,                    &
         & WAVE_UB_BOT  

    !--Parameters in NameList NML_RIVER_TYPE
    CHARACTER(LEN=80) RIVER_TS_SETTING ! METHOD TO CALCULATE T&S AT
    ! THE RIVER MOUTH
    CHARACTER(LEN=80) RIVER_INFLOW_LOCATION ! IS THE FLUX
    CHARACTER(LEN=80) RIVER_KIND ! PERIODIC OR VARIABLE
    CHARACTER(LEN=80) RIVER_INFO_FILE ! RIVERS INFORMATION FILE
    ! SPECIFIED FOR NODES OR EDGES
    INTEGER RIVER_NUMBER ! THE NUMBER OF RIVERS

    NAMELIST /NML_RIVER_TYPE/               &
         & RIVER_NUMBER,                    &
         & RIVER_KIND,                      &
         & RIVER_TS_SETTING,                &
         & RIVER_INFO_FILE,                 &
         & RIVER_INFLOW_LOCATION
    

    !--Parameters in NameList NML_RIVER
    INTEGER, PARAMETER :: RIVER_CHAR_LEN=60

    INTEGER, PARAMETER :: MAX_LAYERS = 100

    CHARACTER(LEN=80) RIVER_NAME
    CHARACTER(LEN=80) RIVER_FILE
    INTEGER RIVER_GRID_LOCATION

# if defined(RIVER_FLOAT)
    ! ASSUMING MAXIMUM NUMBER OF LAYERS!
    REAL(SP) :: RIVER_VERTICAL_DISTRIBUTION(MAX_LAYERS)
# else
    CHARACTER(LEN=RIVER_CHAR_LEN) RIVER_VERTICAL_DISTRIBUTION
# endif

    NAMELIST /NML_RIVER/                    &
         & RIVER_NAME,                      &
         & RIVER_FILE,                      &
         & RIVER_GRID_LOCATION,             &        
         & RIVER_VERTICAL_DISTRIBUTION

    
    ! THIS TYPE PROVIDES INTERMEDIATE STORAGE FOR RIVER INFO UNTILL
    ! THE RIVER FORCING IS SET UP IN MOD_FORCE
    type RIVER
       CHARACTER(LEN=80) NAME
       CHARACTER(LEN=80) FILE
       INTEGER LOCATION
# if defined(RIVER_FLOAT)
       REAL(SP) :: DISTRIBUTION(MAX_LAYERS)
# else
       CHARACTER(LEN=RIVER_CHAR_LEN) DISTRIBUTION
# endif
    end type RIVER

    type(RIVER), Allocatable, DIMENSION(:) :: RIVERS

    !--Parameters in NameList NML_OPEN_BOUNDARY
    LOGICAL OBC_ON
    CHARACTER(LEN=80) OBC_NODE_LIST_FILE
    LOGICAL OBC_ELEVATION_FORCING_ON
    CHARACTER(LEN=80) OBC_ELEVATION_FILE
    INTEGER OBC_TS_TYPE
    LOGICAL OBC_TEMP_NUDGING
    CHARACTER(LEN=80) OBC_TEMP_FILE
    REAL(SP) OBC_TEMP_NUDGING_TIMESCALE
    LOGICAL OBC_SALT_NUDGING
    CHARACTER(LEN=80) OBC_SALT_FILE
    REAL(SP) OBC_SALT_NUDGING_TIMESCALE
#   if defined (WATER_QUALITY)
    LOGICAL OBC_WQM_NUDGING
    CHARACTER(LEN=80) OBC_WQM_FILE
    REAL(SP) OBC_WQM_NUDGING_TIMESCALE
#   endif    
#   if defined (BioGen)
    LOGICAL OBC_BIO_NUDGING
    CHARACTER(LEN=80) OBC_BIO_FILE
    REAL(SP) OBC_BIO_NUDGING_TIMESCALE
#   endif    
    LOGICAL OBC_MEANFLOW
    CHARACTER(LEN=80) OBC_MEANFLOW_FILE
    LOGICAL OBC_LONGSHORE_FLOW_ON
    CHARACTER(LEN=80) OBC_LONGSHORE_FLOW_FILE
    INTEGER OBC_TIDEOUT_INITIAL,OBC_TIDEOUT_INTERVAL
    LOGICAL OBC_DEPTH_CONTROL_ON

    NAMELIST /NML_OPEN_BOUNDARY_CONTROL/    &
         & OBC_ON,                          &
         & OBC_NODE_LIST_FILE,              &
         & OBC_ELEVATION_FORCING_ON,        &
         & OBC_ELEVATION_FILE,              &
         & OBC_TS_TYPE,                     &
         & OBC_TEMP_NUDGING,                &
         & OBC_TEMP_FILE,                   &
         & OBC_TEMP_NUDGING_TIMESCALE,      &
         & OBC_SALT_NUDGING,                &
         & OBC_SALT_FILE,                   &
         & OBC_SALT_NUDGING_TIMESCALE,      &
#   if defined (WATER_QUALITY)
         & OBC_WQM_NUDGING,                 &
         & OBC_WQM_FILE,                    &
         & OBC_WQM_NUDGING_TIMESCALE,       &
#   endif    
#   if defined (BioGen)
         & OBC_BIO_NUDGING,                 &
         & OBC_BIO_FILE,                    &
         & OBC_BIO_NUDGING_TIMESCALE,       &
#   endif    
         & OBC_MEANFLOW,                    &
         & OBC_MEANFLOW_FILE,               &
         & OBC_TIDEOUT_INITIAL,             &
         & OBC_TIDEOUT_INTERVAL,            &
         & OBC_LONGSHORE_FLOW_ON,           &
         & OBC_LONGSHORE_FLOW_FILE,         &
         & OBC_DEPTH_CONTROL_ON


    !--Parameters in NameList GRID_COORDINATES
    CHARACTER(LEN=80) GRID_FILE
    CHARACTER(LEN=80) GRID_FILE_UNITS
    CHARACTER(LEN=200) PROJECTION_REFERENCE
    CHARACTER(LEN=80) SIGMA_LEVELS_FILE
    CHARACTER(LEN=80) DEPTH_FILE
    CHARACTER(LEN=80) CORIOLIS_FILE
    CHARACTER(LEN=80) SPONGE_FILE

    LOGICAL USE_PROJ
! THESE ARE THE FVCOM NAMES FOR 
     !SIGMA_LEVELS = KB
     !SIGMA_LAYERS = KBM1 = KB - 1
     !SIGMA_LEVELS - 2 = KBM2

    NAMELIST /NML_GRID_COORDINATES/        &
         & GRID_FILE,                      &
         & GRID_FILE_UNITS,                &
         & PROJECTION_REFERENCE,           &
         & SIGMA_LEVELS_FILE,              &
         & DEPTH_FILE,                     &
         & CORIOLIS_FILE,                  &
         & SPONGE_FILE

    !--Parameters in NameList NML_GROUNDWATER
    LOGICAL GROUNDWATER_ON
    CHARACTER(LEN=80) GROUNDWATER_KIND
    CHARACTER(LEN=80) GROUNDWATER_FILE
    REAL(SP) GROUNDWATER_FLOW
    REAL(SP) GROUNDWATER_TEMP
    LOGICAL GROUNDWATER_TEMP_ON
    REAL(SP) GROUNDWATER_SALT
    LOGICAL GROUNDWATER_SALT_ON



    NAMELIST /NML_GROUNDWATER/            &
         & GROUNDWATER_ON,                &
         & GROUNDWATER_TEMP_ON,           &
         & GROUNDWATER_SALT_ON,           &
         & GROUNDWATER_KIND,              &
         & GROUNDWATER_FILE,              &
         & GROUNDWATER_FLOW,              &
         & GROUNDWATER_TEMP,              &
         & GROUNDWATER_SALT

    !--Parameters in NameList NML_LAG_PART
    LOGICAL LAG_PARTICLES_ON
    CHARACTER(LEN=80) LAG_START_FILE
    CHARACTER(LEN=80) LAG_OUT_FILE
    CHARACTER(LEN=80) LAG_FIRST_OUT
    CHARACTER(LEN=80) LAG_RESTART_FILE
    CHARACTER(LEN=80) LAG_OUT_INTERVAL 
    CHARACTER(LEN=80) LAG_SCAL_CHOICE

!!$    LOGICAL LAG_TEMPERATURE
!!$    LOGICAL LAG_SALINITY
!!$    LOGICAL LAG_DENSITY
!!$    LOGICAL LAG_EDDY_VISCOSITY
!!$    LOGICAL LAG_DIFFUSIVITY
  
    NAMELIST /NML_LAG/                    &
         & LAG_PARTICLES_ON,              &
         & LAG_START_FILE,                &
         & LAG_OUT_FILE,                  &
         & LAG_FIRST_OUT,                 &
         & LAG_RESTART_FILE,              &
         & LAG_OUT_INTERVAL,              &
         & LAG_SCAL_CHOICE


!!$         & LAG_TEMPERATURE,               &
!!$         & LAG_SALINITY,                  &
!!$         & LAG_DENSITY,                   &
!!$         & LAG_EDDY_VISCOSITY,            &
!!$         & LAG_DIFFUSIVITY

    !--Parameters in NameList NML_ADDITIONAL_MODELS
!    LOGICAL WATER_QUALITY_MODEL
!    CHARACTER(LEN=80) WATER_QUALITY_MODEL_FILE
    LOGICAL DATA_ASSIMILATION
    CHARACTER(LEN=80) DATA_ASSIMILATION_FILE
    LOGICAL BIOLOGICAL_MODEL
! J. Ge for online biology
    CHARACTER(LEN=80) BIOLOGICAL_MODEL_FILE
!--------------------------
# if defined (OFFLINE_BIOLOGY)
    CHARACTER(LEN=80) OFFLINE_BIOLOGY_FILE
# endif

    CHARACTER(LEN=80) STARTUP_BIO_TYPE  !!TYPE OF BIO START
    LOGICAL SEDIMENT_MODEL
    CHARACTER(LEN=80) SEDIMENT_MODEL_FILE
# if defined (OFFLINE_SEDIMENT)
    CHARACTER(LEN=80) OFFLINE_SEDIMENT_FILE
# endif
    CHARACTER(LEN=80) SEDIMENT_PARAMETER_TYPE
    CHARACTER(LEN=80) SEDIMENT_PARAMETER_FILE
    CHARACTER(LEN=80) BEDFLAG_TYPE 
    CHARACTER(LEN=80) BEDFLAG_FILE 

    LOGICAL ICING_MODEL
    CHARACTER(LEN=80) ICING_FORCING_FILE
    CHARACTER(LEN=80) ICING_FORCING_KIND
    REAL(SP) ::ICING_AIR_TEMP
    REAL(SP) ::ICING_WSPD
    
    LOGICAL ICE_MODEL
    CHARACTER(LEN=80) ICE_FORCING_FILE
    CHARACTER(LEN=80) ICE_FORCING_KIND
    REAL(SP) ::ICE_SEA_LEVEL_PRESSURE
    REAL(SP) ::ICE_AIR_TEMP
    REAL(SP) ::ICE_SPEC_HUMIDITY
    REAL(SP) ::ICE_CLOUD_COVER
    REAL(SP) ::ICE_SHORTWAVE
    CHARACTER(LEN=80) ICE_LONGWAVE_TYPE
    
    LOGICAL HIGH_LATITUDE_WAVE

    !ADDITIONAL MODEL DATA
    INTEGER            :: N_SED
    INTEGER, PARAMETER :: N_SED_MAX = 10
    CHARACTER(LEN=20)  :: SED_NAMES(N_SED_MAX)
    REAL(SP), ALLOCATABLE :: SEDDIS(:,:)

#   if defined (WATER_QUALITY)
    INTEGER, PARAMETER :: NB     = 8
    CHARACTER(LEN=80) :: WQM_NAME(NB,4)
#   endif

    REAL(SP), ALLOCATABLE :: BIODIS(:,:)


    ! FVCOM RUN MODE PARAMETERS
    CHARACTER(LEN=80) FVCOM_RUN_MODE
    CHARACTER(LEN=80),parameter :: FVCOM_PURE_SIM         = 'pure sim'
!    CHARACTER(LEN=80),parameter :: FVCOM_NUDGE_AVG_SST = 'nudge avg sst'
!    CHARACTER(LEN=80),parameter :: FVCOM_NUDGE_AVG_TSGRD = 'nudge avg ts'
    CHARACTER(LEN=80),parameter :: FVCOM_NUDGE_OI_ASSIM   = 'nudge or OI assim'   
!    CHARACTER(LEN=80),parameter :: FVCOM_OI_ASSIM      = 'OI ASSIM'
!    CHARACTER(LEN=80),parameter :: FVCOM_KALMAN_RRKF   = 'Kalman RRKF'
    CHARACTER(LEN=80),parameter :: FVCOM_RRKF_WITHOUT_SSA = 'RRKF WITHOUT SSH/SST'
    CHARACTER(LEN=80),parameter :: FVCOM_RRKF_WITH_SSA    = 'RRKF WITH SSH/SST'
    CHARACTER(LEN=80),parameter :: FVCOM_ENKF_WITHOUT_SSA = 'ENKF WITHOUT SSH/SST'
    CHARACTER(LEN=80),parameter :: FVCOM_ENKF_WITH_SSA    = 'ENKF WITH SSH/SST'
    CHARACTER(LEN=80),parameter :: FVCOM_KALMAN_4         = 'Kalman 4'

    NAMELIST /NML_ADDITIONAL_MODELS/      &
!         & WATER_QUALITY_MODEL,           &
!         & WATER_QUALITY_MODEL_FILE,      &
         & DATA_ASSIMILATION,             &
         & DATA_ASSIMILATION_FILE,        &
         & BIOLOGICAL_MODEL,              &
         & STARTUP_BIO_TYPE,              &
!--------- J. Ge for biology --------------
         & BIOLOGICAL_MODEL_FILE,         &
!------------------------------------------
# if defined (OFFLINE_BIOLOGY)
         & OFFLINE_BIOLOGY_FILE,         &
# endif
         & SEDIMENT_MODEL,                &
         & SEDIMENT_MODEL_FILE,           &
# if defined (OFFLINE_SEDIMENT)
         & OFFLINE_SEDIMENT_FILE,         &
# endif
         & SEDIMENT_PARAMETER_TYPE,       &
         & SEDIMENT_PARAMETER_FILE,       &
         & BEDFLAG_TYPE,                  &
         & BEDFLAG_FILE,                  &
         & ICING_MODEL,                   &
         & ICING_FORCING_FILE,            &
         & ICING_FORCING_KIND,            &
         & ICING_AIR_TEMP,                &
         & ICING_WSPD,                    &
         & ICE_MODEL,                     &
         & ICE_FORCING_FILE,              &
         & ICE_FORCING_KIND,              &
         & ICE_SEA_LEVEL_PRESSURE,        &
         & ICE_AIR_TEMP,                  &
         & ICE_SPEC_HUMIDITY,             &
         & ICE_SHORTWAVE,                 &
         & ICE_LONGWAVE_TYPE,                 &
         & ICE_CLOUD_COVER,               &
         & HIGH_LATITUDE_WAVE

    !--Parameters in NameList NML_PROBE
    LOGICAL PROBES_ON
    INTEGER PROBES_NUMBER
    CHARACTER(len=80) PROBES_FILE
    
    NAMELIST /NML_PROBES/      &
         & PROBES_ON,          &
         & PROBES_NUMBER,      &
         & PROBES_FILE

    !--Parameters in NameList NML_BOUNDSCHK
    !=> bounds checking
    LOGICAL  :: FORCE_ARCHIVE = .false.
    LOGICAL  :: BOUNDSCHK_ON = .false.
    INTEGER  :: CHK_INTERVAL
    REAL(SP) :: VELOC_MAG_MAX
    REAL(SP) :: ZETA_MAG_MAX
    REAL(SP) :: TEMP_MAX
    REAL(SP) :: TEMP_MIN
    REAL(SP) :: SALT_MAX
    REAL(SP) :: SALT_MIN
#   if defined (WAVE_CURRENT_INTERACTION)
    REAL(SP) :: HS_MAX = 1e9
#   endif
  
    NAMELIST /NML_BOUNDSCHK/      &
           & BOUNDSCHK_ON,        &
           & CHK_INTERVAL,        &
           & VELOC_MAG_MAX,       &
           & ZETA_MAG_MAX,        &
#  if defined(WAVE_CURRENT_INTERACTION)
           & HS_MAX,              &
#  endif
           & TEMP_MAX,            &
           & TEMP_MIN,            &
           & SALT_MAX,            &
           & SALT_MIN
    !<= bounds checking

!--Time Variables for FVCOM-----------------------------------------!
    TYPE(TIME) :: IntTime
    TYPE(TIME) :: ExtTime
    Type(TIME) :: RKTime
    TYPE(TIME) :: ZEROTIME

    TYPE(TIME) :: EndTime
    Type(TIME) :: StartTime
    Type(TIME) :: RUNFILE_StartTime
    Type(TIME) :: ReferenceDate

    TYPE(TIME) :: DELT_RHO_MEAN
    TYPE(TIME) :: RECALC_RHO_MEAN

    ! ZERO PHASE TIME FOR SPECTRAL (NON JULIAN) TIDE
    TYPE(TIME) :: SPECTIME
    
    REAL(SP)   :: DTE        !!EXTERNAL TIME STEP (Seconds)
    REAL(SP)   :: DTI        !!INTERNAL TIME STEP (Seconds)
    REAL(SP)   :: RAMP

    TYPE(TIME) :: IMDTE      !!EXTERNAL TIME STEP 
    TYPE(TIME) :: IMDTI      !!INTERNAL TIME STEP 
    
    INTEGER(itime) :: IINT     !!INTERNAL TIME STEP ITERATION NUMBER (ISTART => IEND)
    INTEGER(itime) :: IEXT     !!EXTERNAL TIME STEP ITERATION NUMBER (1 => ISPLIT)
    INTEGER(itime) :: ISTART   !!STARTING INTERNAL TIME STEP ITERATION NUMBER
    INTEGER(itime) :: IEND     !!ENDING INTERNAL TIME STEP ITERATION NUMBER
    INTEGER(itime) :: NSTEPS   !!Number OF INTERAL TIME STEPS IN SIMULATION
    
    ! Time variables for File IO
    
    INTEGER, PARAMETER:: TIMEPREC = 6
    !! THIS IS THE LENGHT THAT LOOKS NICE FOR GIVEN TIME PREC
    INTEGER, PARAMETER:: DATESTRLEN = 20+TIMEPREC 

!    CHARACTER(LEN=80) :: IO_FILE_DATE
!    CHARACTER(LEN=80) :: IO_timestr
!    real(SP)          :: IO_days
!    integer           :: IO_mjd
!    integer           :: IO_msec
!    integer           :: IO_IINT



!--Constants-------------------------------------------------------------------!

    ! SELECT WHETHER YOU WANT TO MAKE TIME MORE EXACT (SP CAUSES
    ! ROUND OF ERROR BUT IS THE TRADITIONAL PRECISION)
   REAL(DP), PARAMETER, DIMENSION(4) :: ALPHA_RK = (/0.2500_DP,1.0_DP/3.0_DP,0.5000_DP,1.0_DP/)
!   REAL(SP), PARAMETER, DIMENSION(4) :: ALPHA_RK = (/0.2500_DP,1.0_DP/3.0_DP,0.5000_DP,1.0_DP/)


   REAL(DP), PARAMETER :: GRAV      = 9.81_SP
   REAL(DP), PARAMETER :: PI        = 3.141592653589793238_DP
   REAL(DP), PARAMETER :: PI2       = 2.0_DP * 3.141592653589793238_DP
   REAL(DP), PARAMETER :: ZERO      = 0.0_DP 
   REAL(DP), PARAMETER :: ONE_THIRD = 1.0_DP/3.0_DP 
   REAL(DP), PARAMETER :: REARTH    = 6371.0E03_DP   !!Earth Radius in Meters
   REAL(DP), PARAMETER :: DEG2RAD   = PI2/360.0_DP   !!Radians/Degree
   REAL(DP), PARAMETER :: TPI       = DEG2RAD*REARTH !TPI=pi*rearth/180.=3.14159265/180.0*6371.*1000.
   REAL(DP), PARAMETER :: ROFVROS   = 0.9775171065_DP!!RATIO OF THE DENSITY OF FRESH AND SEA WATER 1000./1023.   
   real(DP), parameter :: SLP0      = 101325.0_sp    !! mean sea surface pressure (Pa)


!--Parameter Controlling Vertical Coordinate Distribution------------
   !----------!
   CHARACTER(LEN=80) :: STYPE
   CHARACTER(LEN=80), PARAMETER :: STYPE_UNIFORM= "UNIFORM"
   CHARACTER(LEN=80), PARAMETER :: STYPE_GEOMETRIC= "GEOMETRIC"
   CHARACTER(LEN=80), PARAMETER :: STYPE_TANH= "TANH"
   CHARACTER(LEN=80), PARAMETER :: STYPE_GENERALIZED= "GENERALIZED"
   CHARACTER(LEN=80), PARAMETER :: STYPE_RESTART= "RESTART"
!--Sigma Level Parameters for case GEOMETRIC OF UNIFORM------------------------!
   REAL(SP) :: P_SIGMA      !!PARAMETER CONTROLLING SIGMA LEVEL DISTRIBUTION

!--General Vertical Level Parameters for case TANH ----------------------------!
   REAL(SP) :: DU2          !!PARAMETER CONTROLLING LEVEL DISTRIBUTION OF SURFACE
   REAL(SP) :: DL2          !!PARAMETER CONTROLLING LEVEL DISTRIBUTION OF BOTTOM

!--General Vertical Level Parameters for case GENERALIZED ---------------------!
   REAL(SP) :: DUU          !!THE UPPER BOUNDARY OF PARALLEL COORDINATE
   REAL(SP) :: DLL          !!THE LOWER BOUNDARY OF PARALLEL COORDINATE
   REAL(SP) :: HMIN1        !!THE MIN DEPTH AT WHICH THE LAYERS ARE CONSTANT
   
   INTEGER  :: KU           !!THE NUMBERS OF LAYERS ABOVE UPPER BOUNDARY
   INTEGER  :: KL           !!THE NUMBVER OF LAYERS BELOW LOWER BOUNDARY
   
   REAL(SP), ALLOCATABLE :: ZKU(:)       !!THE DEPTHS OF PARALLEL LAYERS ABOVE UPPER BOUNDARY
   REAL(SP), ALLOCATABLE :: ZKL(:)       !!THE DEPTHS OF PARALLEL LAYERS BELOW LOWER BOUNDARY

   REAL(SP) :: HMAX         !!GLOBAL MAXIMUM DEPTH IN DEPTH FILE
   REAL(SP) :: HMIN         !!GLOBAL MINIMUM DEPTH IN DEPTH FILE


   ! All file identifiers go here
   INTEGER :: IPT                     !! IUNIT FOR LOG FILE OUTPUT
   INTEGER, PUBLIC, PARAMETER :: IPT_BASE= 7000 ! FOR PAR LOG FILES

   INTEGER, PARAMETER :: TESTUNIT = 200 !! TEST TO SEE IF OUTPUT DIR EXISTS/WRITABLE
   INTEGER, PARAMETER :: NMLUNIT = 10 !! NAMELIST RUN FILE
   INTEGER, PARAMETER :: ITSUNIT=11
   INTEGER, PARAMETER :: OBCUNIT=12
   INTEGER, PARAMETER :: GRIDUNIT=13
   INTEGER, PARAMETER :: SIGMAUNIT=14
   INTEGER, PARAMETER :: DEPTHUNIT=15
   INTEGER, PARAMETER :: CORIOLISUNIT=16
   INTEGER, PARAMETER :: SPONGEUNIT=17
   INTEGER, PARAMETER :: LSFUNIT=18
   INTEGER, PARAMETER :: ASSIMUNIT=19
   INTEGER, PARAMETER :: OIASSIMUNIT=23
   INTEGER, PARAMETER :: PROBEUNIT=20
   INTEGER, PARAMETER :: JULOBCUNIT=21
   INTEGER, PARAMETER :: KFUNIT=22
   INTEGER, PARAMETER :: NESTUNIT=30
   INTEGER, PARAMETER :: SUBDUNIT=31

   INTEGER            :: RIVERNMLUNIT

# if defined (PLBC) 
   INTEGER, PARAMETER :: INTCELL1 =46
   INTEGER, PARAMETER :: INTNODE1 =47
# endif

END MODULE CONTROL

!==============================================================================|
!==============================================================================|
!==============================================================================|
!     ALL VARS                                                                 |
!     CONATINS:
!      FVCOM VARIABLES
!      N2E3D: simple average from nodes to elements for 3D variables
!      N2E2D: simple average from nodes to elements for 2D variables
!      E2N3D: simple average from elements to nodes for 3D variables
!      E2N2D: simple average from elements to nodes for 2D variables
!==============================================================================|
MODULE ALL_VARS
   USE MOD_PREC
   USE LIMS
   USE CONTROL 
   IMPLICIT NONE
   SAVE


!--------------------------Temporary Array------------------------------------------!

  INTEGER, ALLOCATABLE :: NVG(:,:)

!--------------------------Global Grid Variables------------------------------------!

!!  REAL(SP), POINTER     :: XG(:)               !!GLOBAL X-COORD AT NODE 
!!  REAL(SP), POINTER     :: YG(:)               !!GLOBAL X-COORD AT NODE 
  REAL(SP), ALLOCATABLE     :: XG(:)               !!GLOBAL X-COORD AT NODE 
  REAL(SP), ALLOCATABLE     :: YG(:)               !!GLOBAL X-COORD AT NODE 
  REAL(SP), ALLOCATABLE :: HG(:)               !!GLOBAL DEPTH AT NODE
!  REAL(SP), ALLOCATABLE :: CORG(:)             !!GLOBAL COORIOLIS AT NODE 
  REAL(SP), ALLOCATABLE :: XCG(:)              !!GLOBAL X-COORD AT FACE CENTER 
  REAL(SP), ALLOCATABLE :: YCG(:)              !!GLOBAL X-COORD AT FACE CENTER 

!--------------------------Grid Metrics---------------------------------------------!

  ! JUST DON'T ALLOCATE THINGS THAT YOU DON'T NEED IN A PARTICULAR
  ! MODULE IT IS A PAIN IN THE NECK TO NO HAVE THE VARIABLE DECLARED
  

!# if !defined(SPHERICAL)
   REAL(SP)              :: VXMIN,VYMIN,VXMAX,VYMAX
!# endif
   REAL(SP), ALLOCATABLE,TARGET :: XM(:)               !!X-COORD AT NODE IN METERS
   REAL(SP), ALLOCATABLE,TARGET :: YM(:)               !!Y-COORD AT NODE IN METERS
   REAL(SP), ALLOCATABLE,TARGET :: XMC(:)              !!X-COORD AT CELL CENTER IN METERS
   REAL(SP), ALLOCATABLE,TARGET :: YMC(:)              !!Y-COORD AT CELL CENTER IN METERS
   REAL(SP), ALLOCATABLE,TARGET :: LON(:)              !!LONGITUDE AT THE NODE
   REAL(SP), ALLOCATABLE,TARGET :: LAT(:)              !!LATITUDE AT THE NODE
   REAL(SP), ALLOCATABLE,TARGET :: LONC(:)              !!LONGITUDE AT THE NODE
   REAL(SP), ALLOCATABLE,TARGET :: LATC(:)              !!LATITUDE AT THE NODE
   ! VX,VY and XC,YC are used for either meters or spherical depending on
   ! make file option 'SPHERICAL'
   REAL(SP), ALLOCATABLE,TARGET :: VX(:)               !!X-COORD AT GRID POINT
   REAL(SP), ALLOCATABLE,TARGET :: VY(:)               !!Y-COORD AT GRID POINT
   REAL(SP), ALLOCATABLE,TARGET :: XC(:)               !!X-COORD AT FACE CENTER 
   REAL(SP), ALLOCATABLE,TARGET :: YC(:)               !!Y-COORD AT FACE CENTER

   ! NOTES: SHOULD MAKE AN ARRAY TO STORE 1/ART, 1/ART2 and 1/ART2
   ! IT is faster and safer

   REAL(SP), ALLOCATABLE,TARGET :: ART(:)              !!AREA OF ELEMENT
   REAL(SP), ALLOCATABLE,TARGET :: ART1(:)             !!AREA OF NODE-BASE CONTROl VOLUME
   REAL(SP), ALLOCATABLE,TARGET :: ART2(:)             !!AREA OF ELEMENTS AROUND NODE
!--Gravity (vary with latitute) -----------------------------------------------!
   REAL(SP),ALLOCATABLE,TARGET :: GRAV_N(:),GRAV_E(:) ! CALCULATED AS A
   ! FUNCTION OF LATITUDE IN SPHERICAL COORDINATES MODEL
   
!----------------Node, Boundary Condition, and Control Volume-----------------------!

   INTEGER, ALLOCATABLE,TARGET :: NV(:,:)             !!NODE NUMBERING FOR ELEMENTS
!   INTEGER, ALLOCATABLE,TARGET :: NVGL(:,:)             !!NODE GLOBAL NUMBERING OF LOCAL ELEMENTS
   INTEGER, ALLOCATABLE,TARGET :: NBE(:,:)            !!INDICES OF ELMNT NEIGHBORS
!   INTEGER, POINTER            :: NBEGL(:,:)          !!GLOBAL INDICES OF LOCAL ELEMENT NEIGHBORS
   INTEGER, ALLOCATABLE,TARGET :: NTVE(:)             !! NUMBER OF ELEMENTS SURROUNDING EACH NODE
   INTEGER, ALLOCATABLE,TARGET :: NTSN(:)             !! NUMBER OF NODES SURROUNDING EACH NODE
   INTEGER, ALLOCATABLE,TARGET :: ISONB(:)            !!NODE MARKER = 0,1,2
   INTEGER, ALLOCATABLE,TARGET :: ISONB_W(:)          !!NODE MARKER = 0,1,2
   INTEGER, ALLOCATABLE,TARGET :: ISBC(:)     
   INTEGER, ALLOCATABLE,TARGET :: ISBCE(:)     
   INTEGER, ALLOCATABLE,TARGET :: IEC(:,:)
   INTEGER, ALLOCATABLE,TARGET :: IENODE(:,:)
   INTEGER, ALLOCATABLE,TARGET :: NBSN(:,:)            !! INDICES OF NODES SURROUNDING EACH NODE
!   INTEGER, POINTER            :: NBSNGL(:,:)          !! GLOBAL INDICIES OF NODES SURROUNDING EACH LOCAL NODE
   INTEGER, ALLOCATABLE,TARGET :: NIEC(:,:)
   INTEGER, ALLOCATABLE,TARGET :: NTRG(:)
   INTEGER, ALLOCATABLE,TARGET :: NBVE(:,:)            !! INDICIES OF ELEMENTS SURROUNDING EACH NODE
!   INTEGER, POINTER            :: NBVEGL(:,:)          !! GLOBAL INDICIES OF ELEMENTS SURROUNDING EACH LOCAL NODE
   INTEGER, ALLOCATABLE,TARGET :: NBVT(:,:)
   INTEGER, ALLOCATABLE,TARGET :: LISBCE_1(:)          !!LIST OF ELEMENTS WITH ISBCE=1
   INTEGER, ALLOCATABLE,TARGET :: LISBCE_2(:)          !!LIST OF ELEMENTS WITH ISBCE=2
   INTEGER, ALLOCATABLE,TARGET :: LISBCE_3(:)          !!LIST OF ELEMENTS WITH ISBCE=3
   REAL(SP),ALLOCATABLE,TARGET :: DLTXC(:)
   REAL(SP),ALLOCATABLE,TARGET :: DLTYC(:)
   REAL(SP),ALLOCATABLE,TARGET :: DLTXYC(:)
   REAL(SP),ALLOCATABLE,TARGET :: SITAE(:) 
   REAL(SP),ALLOCATABLE,TARGET :: XIJC(:) 
   REAL(SP),ALLOCATABLE,TARGET :: YIJC(:)
   ! POSITION OF NODAL CONTROL VOLUME CORNERS 
   REAL(SP),ALLOCATABLE,TARGET :: XIJE(:,:) 
   REAL(SP),ALLOCATABLE,TARGET :: YIJE(:,:) 
   ! LENGTH OF NODAL CONTROL VOLUME EDGES
   REAL(SP),ALLOCATABLE,TARGET :: DLTXE(:)
   REAL(SP),ALLOCATABLE,TARGET :: DLTYE(:)
   REAL(SP),ALLOCATABLE,TARGET :: DLTXYE(:)
   REAL(SP),ALLOCATABLE,TARGET :: SITAC(:) 


   REAL(SP),ALLOCATABLE,TARGET :: EPOR(:)            !!ELEMENT FLUX POROSITY (=0. IF ISBCE = 2)

   ! LENGTH BETWEEN NODE AND CONTROL VOLUMEN EDGE CENTER
   REAL(SP),ALLOCATABLE,TARGET :: DLTXNCVE(:,:)!! DeLTa X Node to Control Volume Edge 
   REAL(SP),ALLOCATABLE,TARGET :: DLTYNCVE(:,:)!! DeLTa Y Node to Control Volume Edge 


#  if defined (SPHERICAL)
   REAL(DP),ALLOCATABLE,TARGET :: DLTYTRIE(:,:)!! DeLTa Y TRIangle Edge
   REAL(DP),ALLOCATABLE,TARGET :: DLTXTRIE(:,:)!! DeLTa X TRIangle Edge

   REAL(DP),ALLOCATABLE,TARGET :: DLTXECEC(:,:)!! DeLTa X Edge Center to Edge Center
   REAL(DP),ALLOCATABLE,TARGET :: DLTYECEC(:,:)!! DeLTa Y Edge Center to Edge Center

   REAL(DP),ALLOCATABLE,TARGET :: DLTXNEC(:,:)!! DeLTa X Node to Edge Center
   REAL(DP),ALLOCATABLE,TARGET :: DLTYNEC(:,:)!! DeLTa Y Node to Edge Center
# else
   REAL(SP),ALLOCATABLE,TARGET :: DLTYTRIE(:,:)!! DeLTa Y TRIangle Edge
   REAL(SP),ALLOCATABLE,TARGET :: DLTXTRIE(:,:)!! DeLTa X TRIangle Edge
   
   REAL(SP),ALLOCATABLE,TARGET :: DLTXECEC(:,:)!! DeLTa X Edge Center to Edge Center
   REAL(SP),ALLOCATABLE,TARGET :: DLTYECEC(:,:)!! DeLTa Y Edge Center to Edge Center

   REAL(SP),ALLOCATABLE,TARGET :: DLTXNEC(:,:)!! DeLTa X Node to Edge Center
   REAL(SP),ALLOCATABLE,TARGET :: DLTYNEC(:,:)!! DeLTa Y Node to Edge Center
# endif   


   ! LONG SHORE FLOW VARIABLES
   INTEGER, ALLOCATABLE,TARGET :: IBCLSF_GL(:)     !!GLOBAL NODE NUMBER OF LSF BOUNDARY   
   REAL(SP),ALLOCATABLE,TARGET :: RBC_GEO_GL(:)     !!GLOBAL GEOSTROPHIC FRICTION CORRECTION NODES
   REAL(SP),ALLOCATABLE,TARGET :: RBC_WDF_GL(:)     !!GLOBAL WIND DRIVEN FLOW CORRECTION NODES

   REAL(SP),ALLOCATABLE,TARGET :: WDF_ANG(:)       !!ANGLE ALLONG THE OPEN BOUNDARY
   REAL(SP),ALLOCATABLE,TARGET :: WDF_DIST(:)      !!DISTANCE ALLONG THE OPEN BOUNDARY 
   INTEGER, ALLOCATABLE,TARGET :: IBCLSF(:)        !!LOCAL NODE NUMBER OF LSF BOUNDARY   
   INTEGER, ALLOCATABLE,TARGET :: IBCLSF_OUTPUT(:) !! LIST OF LOCAL LSF NODES GLOBAL NUMBER FOR OUTPUT   
   INTEGER, ALLOCATABLE,TARGET :: NBCLSF(:)        !!LOCAL NODE NUMBER OF THE NEXT LSF BOUNDARY   
   REAL(SP),ALLOCATABLE,TARGET :: RBC_GEO(:)        !!LOCAL GEOSTROPHIC FRICTION CORRECTION NODES
   REAL(SP),ALLOCATABLE,TARGET :: RBC_WDF(:)        !!LOCAL WIND DRIVEN FLOW CORRECTION NODES

!   INTEGER, ALLOCATABLE,TARGET :: N_ICELLQ(:,:)    !!FLUX ANGLE 

!----------------2-d arrays for the general vertical coordinate -------------------------------!

   REAL(SP), ALLOCATABLE,TARGET :: Z(:,:)                    !!SIGMA COORDINATE VALUE 
   REAL(SP), ALLOCATABLE,TARGET :: ZZ(:,:)                   !!INTRA LEVEL SIGMA VALUE
   REAL(SP), ALLOCATABLE,TARGET :: DZ(:,:)                   !!DELTA-SIGMA VALUE
   REAL(SP), ALLOCATABLE,TARGET :: DZZ(:,:)                  !!DELTA OF INTRA LEVEL SIGMA 
   REAL(SP), ALLOCATABLE,TARGET :: Z1(:,:)                   !!SIGMA COORDINATE VALUE 
   REAL(SP), ALLOCATABLE,TARGET :: ZZ1(:,:)                  !!INTRA LEVEL SIGMA VALUE
   REAL(SP), ALLOCATABLE,TARGET :: DZ1(:,:)                  !!DELTA-SIGMA VALUE
   REAL(SP), ALLOCATABLE,TARGET :: DZZ1(:,:)                 !!DELTA OF INTRA LEVEL SIGMA 
!   REAL(SP), ALLOCATABLE,TARGET :: DPTHSL(:)               !!Z-DEPTHS FOR SALINITY/TEMP ICs


!---------------2-d flow variable arrays at elements-------------------------------!

   REAL(SP), ALLOCATABLE,TARGET :: UA(:)            !!VERTICALLY AVERAGED X-VELOC
   REAL(SP), ALLOCATABLE,TARGET :: VA(:)            !!VERTICALLY AVERAGED Y-VELOC
   REAL(SP), ALLOCATABLE,TARGET :: UAF(:)           !!UA FROM PREVIOUS RK STAGE 
   REAL(SP), ALLOCATABLE,TARGET :: VAF(:)           !!VA FROM PREVIOUS RK STAGE 
!# if !defined (SEMI_IMPLICIT)
   REAL(SP), ALLOCATABLE,TARGET :: UARK(:)          !!UA FROM PREVIOUS TIMESTEP 
   REAL(SP), ALLOCATABLE,TARGET :: VARK(:)          !!VA FROM PREVIOUS TIMESTEP 
   REAL(SP), ALLOCATABLE,TARGET :: UARD(:)          !!UA AVERAGED OVER EXTERNAL INT
   REAL(SP), ALLOCATABLE,TARGET :: VARD(:)          !!VA AVERAGED OVER EXTERNAL INT
!# endif
   REAL(SP), ALLOCATABLE,TARGET :: COR(:)           !!CORIOLIS PARAMETER
   REAL(SP), ALLOCATABLE,TARGET :: F_ALFA(:)
   REAL(SP), ALLOCATABLE,TARGET :: H1(:)            !!BATHYMETRIC DEPTH   
   REAL(SP), ALLOCATABLE,TARGET :: D1(:)            !!CURRENT DEPTH
   REAL(SP), ALLOCATABLE,TARGET :: DT1(:)           !!DEPTH AT PREVIOUS TIME STEP
   REAL(SP), ALLOCATABLE,TARGET :: EL1(:)           !!CURRENT SURFACE ELEVATION
   REAL(SP), ALLOCATABLE,TARGET :: ET1(:)           !!SURFACE ELEVATION AT PREVIOUS TIME STEP
!# if !defined (SEMI_IMPLICIT)
   REAL(SP), ALLOCATABLE,TARGET :: ELRK1(:)         !!SURFACE ELEVATION AT BEGINNING OF RK INT 
!# endif
   REAL(SP), ALLOCATABLE,TARGET :: ELF1(:)          !!SURFACE ELEVATION STORAGE FOR RK INT
   REAL(SP), ALLOCATABLE,TARGET :: DTFA(:)          !!ADJUSTED DEPTH FOR MASS CONSERVATION


   REAL(SP), ALLOCATABLE,TARGET :: CC_SPONGE(:)     !!SPONGE DAMPING COEFFICIENT FOR MOMENTUM
   
!---------------2-d flow variable arrays at nodes----------------------------------!

   REAL(SP), ALLOCATABLE,TARGET :: H(:)     !!BATHYMETRIC DEPTH   
   REAL(SP), ALLOCATABLE,TARGET :: D(:)             !!CURRENT DEPTH   
   REAL(SP), ALLOCATABLE,TARGET :: DT(:)            !!DEPTH AT PREVIOUS TIME STEP
   REAL(SP), ALLOCATABLE,TARGET :: EL(:)            !!CURRENT SURFACE ELEVATION
   REAL(SP), ALLOCATABLE,TARGET :: ET(:)            !!SURFACE ELEVATION AT PREVIOUS TIME STEP
   REAL(SP), ALLOCATABLE,TARGET :: EGF(:)           !!AVERAGE SURFACE ELEVATION OVER EXTERNAL INT
!# if !defined (SEMI_IMPLICIT)
   REAL(SP), ALLOCATABLE,TARGET :: ELRK(:)          !!SURFACE ELEVATION AT BEGINNING OF RK INT
!# endif
   REAL(SP), ALLOCATABLE,TARGET :: ELF(:)           !!SURFACE ELEVATION STORAGE FOR RK INT

   ! DEFINED HERE, BUT ONLY USED IF EQI&ATMO ARE DEFINED
   REAL(SP), ALLOCATABLE,TARGET :: ELF_EQI(:)
!# if !defined (SEMI_IMPLICIT)
   REAL(SP), ALLOCATABLE,TARGET :: ELRK_EQI(:)
!# endif
   REAL(SP), ALLOCATABLE,TARGET :: EL_EQI(:)
   REAL(SP), ALLOCATABLE,TARGET :: EGF_EQI(:)

   REAL(SP), ALLOCATABLE,TARGET :: ELF_ATMO(:)
!# if !defined (SEMI_IMPLICIT)
   REAL(SP), ALLOCATABLE,TARGET :: ELRK_ATMO(:)
!# endif
   REAL(SP), ALLOCATABLE,TARGET :: EL_ATMO(:)
   REAL(SP), ALLOCATABLE,TARGET :: EGF_ATMO(:)

   ! DEFINED HERE, BUT ONLY USED IF AIR PRESSURE ARE DEFINED
   REAL(SP), ALLOCATABLE,TARGET :: ELF_AIR(:)
!# if !defined (SEMI_IMPLICIT)
   REAL(SP), ALLOCATABLE,TARGET :: ELRK_AIR(:)
!# endif
   REAL(SP), ALLOCATABLE,TARGET :: EL_AIR(:)
   REAL(SP), ALLOCATABLE,TARGET :: EGF_AIR(:)
   
   REAL(SP), ALLOCATABLE,TARGET :: VORT(:)


!---------------surface/bottom boundary conditions---------------------------------!

   REAL(SP), ALLOCATABLE,TARGET :: CBC(:)           !!BOTTOM FRICTION     
   REAL(SP), ALLOCATABLE,TARGET :: CC_Z0B(:)        !!BOTTOM ROUGHNESS VARIABLE

   REAL(SP), ALLOCATABLE,TARGET :: SWRAD_WATTS(:)  !!SURFACE INCIDENT RADIATION
   REAL(SP), ALLOCATABLE,TARGET :: WTSURF_WATTS(:) !!NET HEAT FLUX AT SURFACE
#  if defined (HEATING_CALCULATED) || (HEATING_SOLAR) || (ICE)
   !EJA edit to add Sensible, Latent, Longwave heat fluxes for output 05/03/2016
   REAL(SP), ALLOCATABLE,TARGET :: HSENS_WATTS(:)  !!SENSIBLE HEAT FLUX AT SURFACE
   REAL(SP), ALLOCATABLE,TARGET :: HLAT_WATTS(:)   !!LATENT HEAT FLUX AT SURFACE
   REAL(SP), ALLOCATABLE,TARGET :: LWRAD_WATTS(:)  !!LONGWAVE RADIATION AT SURFACE
#  endif

   ! THESE TWO ARE DIVIDED BY THE SPECFIC HEAT AND AVERAGE DENSITY!
   REAL(SP), ALLOCATABLE,TARGET :: SWRAD(:)         !!SURFACE INCIDENT RADIATION
   REAL(SP), ALLOCATABLE,TARGET :: WTSURF(:)        !!NET HEAT FLUX AT SURFACE

   ! SUSPECTED UNITS FOR SURFACE STRESS: N/M * 1/RHO_BAR ... see bcond_gcn.F
!# if !defined (SEMI_IMPLICIT)
   REAL(SP), ALLOCATABLE,TARGET :: WUSURF2(:)       !!SURFACE FRICTION FOR EXT
   REAL(SP), ALLOCATABLE,TARGET :: WVSURF2(:)       !!SURFACE FRICTION FOR EXT
!# endif
   REAL(SP), ALLOCATABLE,TARGET :: WUBOT(:)         !!BOTTOM FRICTION
   REAL(SP), ALLOCATABLE,TARGET :: WVBOT(:)         !!BOTTOM FRICTION
   REAL(SP), ALLOCATABLE,TARGET :: TAUBM(:)         !!BOTTOM FRICTION MAGNITUDE(Caution is Tau' [no Rho])
   REAL(SP), ALLOCATABLE,TARGET :: WUBOT_N(:)       !!BOTTOM FRICTION ON NODES (Caution is Tau' [no Rho])
   REAL(SP), ALLOCATABLE,TARGET :: WVBOT_N(:)       !!BOTTOM FRICTION ON NODES (Caution is Tau' [no Rho])
   REAL(SP), ALLOCATABLE,TARGET :: TAUBM_N(:)       !!BOTTOM FRICTION MAGNITUDE ON NODES (Caution is Tau' [no Rho])
   REAL(SP), ALLOCATABLE,TARGET :: WUSURF(:)        !!SURFACE FRICTION FOR INT
   REAL(SP), ALLOCATABLE,TARGET :: WVSURF(:)        !!SURFACE FRICTION FOR INT
   REAL(SP), ALLOCATABLE,TARGET :: WUSURF_save(:)   !!SURFACE FRICTION FOR INT
   REAL(SP), ALLOCATABLE,TARGET :: WVSURF_save(:)   !!SURFACE FRICTION FOR INT
   ! BFWDIS - UNITS: m3s-1 Cubic meters per second
   REAL(SP), ALLOCATABLE,TARGET :: BFWDIS(:)        !!GROUNDWATER FLUX AT CURRENT TIME
!# if !defined (SEMI_IMPLICIT)
   REAL(SP), ALLOCATABLE,TARGET :: BFWDIS2(:)       !!GROUNDWATER FLUX FOR EXT
!# endif
   REAL(SP), ALLOCATABLE,TARGET :: BFWTMP(:)        !!GROUNDWATER TEMP AT CURRENT TIME
   REAL(SP), ALLOCATABLE,TARGET :: BFWSLT(:)        !!GROUNDWATER SALT AT CURRENT TIME


   ! ICING MODEL DATA
   REAL(SP), ALLOCATABLE,TARGET :: ICING_WNDX(:)
   REAL(SP), ALLOCATABLE,TARGET :: ICING_WNDY(:)
   REAL(SP), ALLOCATABLE,TARGET :: ICING_SATMP(:)
   REAL(SP), ALLOCATABLE,TARGET :: ICING_0kts(:)
   REAL(SP), ALLOCATABLE,TARGET :: ICING_10kts(:)

   ! NOT SURE WHETHER RIVER STUFF BELONGS HERE OR IN ALL_VARS???

   ! RIVER STUFF
   INTEGER,  ALLOCATABLE,TARGET :: INODEQ(:)        !!LOCAL FRESH WATER INFLOW NODES
   INTEGER,  ALLOCATABLE,TARGET :: ICELLQ(:)        !!LOCAL FRESH WATER INFLOW ELEMENT 
   INTEGER,  ALLOCATABLE,TARGET :: N_ICELLQ(:,:)    !!TWO NODES BOUNDING THE INFLOW ELEMENT
   REAL(SP), ALLOCATABLE,TARGET :: VQDIST(:,:)      !!DISCHARGE VERTICAL DISTRIBUTION
   INTEGER,  ALLOCATABLE,TARGET :: RIV_GL2LOC(:)

   REAL(SP), ALLOCATABLE,TARGET :: QDIS(:)          !!RIVER FLUX AT CURRENT TIME
!# if !defined (SEMI_IMPLICIT)
   REAL(SP), ALLOCATABLE,TARGET :: QDIS2(:)         !!RIVER FLUX (EXT MODE, NOT USED) 
!# endif
   REAL(SP), ALLOCATABLE,TARGET :: TDIS(:)          !!RIVER WATER TEMP AT CURRENT TIME
   REAL(SP), ALLOCATABLE,TARGET :: SDIS(:)          !!RIVER WATER SLNT AT CURRENT TIME
   REAL(SP), ALLOCATABLE,TARGET :: QAREA(:)         !!AREA OF RIVER DISCHARGE 
   REAL(SP), ALLOCATABLE,TARGET :: RDISQ(:,:)       !!AREA OF FLUX                    
   REAL(SP), ALLOCATABLE,TARGET :: ANGLEQ(:)        !!RIVER DISCHARGE ANGLE           
   REAL(SP), ALLOCATABLE,TARGET :: VLCTYQ(:)        !!RIVER DISCHARGE VELOCITY
   
   ! SURFACE MET STUFF
   REAL(SP), ALLOCATABLE,TARGET :: UUWIND(:)        !!SURFACE X-WIND 
   REAL(SP), ALLOCATABLE,TARGET :: VVWIND(:)        !!SURFACE Y-WIND
   ! PRECIP/EVAP are in units of meters/second
!# if !defined (SEMI_IMPLICIT)
   REAL(SP), ALLOCATABLE,TARGET :: QPREC2(:)        !!SURFACE PRECIPITATION FOR EXT
   REAL(SP), ALLOCATABLE,TARGET :: QEVAP2(:)        !!SURFACE EVAPORATION FOR EXT
!# endif
   REAL(SP), ALLOCATABLE,TARGET :: QPREC(:)        !!SURFACE PRECIPITATION FOR INT
   REAL(SP), ALLOCATABLE,TARGET :: QEVAP(:)        !!SURFACE EVAPORATION FOR INT

   REAL(SP), ALLOCATABLE,TARGET :: WHS(:)          !!SURFACE WAVE HEIGHT
   REAL(SP), ALLOCATABLE,TARGET :: WDIR(:)         !!SURFACE WAVE DIRECTION
   REAL(SP), ALLOCATABLE,TARGET :: WPER(:)         !!SURFACE WAVE PERIOD
   REAL(SP), ALLOCATABLE,TARGET :: WLENGTH(:)      !!SURFACE WAVE LENGTH
   REAL(SP), ALLOCATABLE,TARGET :: WPER_BOT(:)     !!BOTTOM WAVE PERIOD
   REAL(SP), ALLOCATABLE,TARGET :: WUB_BOT(:)      !!BOTTOM ORBITAL VELOCITY

# if defined (OFFLINE_SEDIMENT) || (OFFLINE_BIOLOGY)
   REAL(SP), ALLOCATABLE,TARGET :: OFFLINE_U(:,:)  !! OFFLINE U VELOCITY
   REAL(SP), ALLOCATABLE,TARGET :: OFFLINE_V(:,:)  !! OFFLINE V VELOCITY
   REAL(SP), ALLOCATABLE,TARGET :: OFFLINE_W(:,:)  !! OFFLINE W VELOCITY
   REAL(SP), ALLOCATABLE,TARGET :: OFFLINE_S1(:,:) !! OFFLINE SALINITY
   REAL(SP), ALLOCATABLE,TARGET :: OFFLINE_T1(:,:) !! OFFLINE TEMPERATURE
   REAL(SP), ALLOCATABLE,TARGET :: OFFLINE_EL(:)   !! OFFLINE ELEVATION
   REAL(SP), ALLOCATABLE,TARGET :: OFFLINE_KH(:,:) !! OFFLINE KH
# if defined (GOTM)
   REAL(SP), ALLOCATABLE,TARGET :: OFFLINE_TEPS(:,:) !! OFFLINE TEPS
# else
   REAL(SP), ALLOCATABLE,TARGET :: OFFLINE_Q2(:,:) !! OFFLINE Q2
   REAL(SP), ALLOCATABLE,TARGET :: OFFLINE_Q2L(:,:)!! OFFLINE Q2L
# endif
   INTEGER , ALLOCATABLE,TARGET :: OFFLINE_WN(:)   !! OFFLINE WET NODES
   INTEGER , ALLOCATABLE,TARGET :: OFFLINE_WC(:)   !! OFFLINE WET CELLS
# endif

!----------------boundary conditions: meteo conditions-----------------!
# if defined (ICE) || (HEATING_CALCULATED) || (HEATING_SOLAR)
   REAL(SP), ALLOCATABLE,TARGET  :: T_AIR(:)         !!BULK AIR TEMPERATURE( deg C ) AT HEIGHT ?(m)
# endif
# if defined (ICE) || (HEATING_CALCULATED)
   REAL(SP), ALLOCATABLE,TARGET  :: RH_AIR(:)        !!RELATIVE HUMIDITY
   REAL(SP), ALLOCATABLE,TARGET  :: DLW_AIR(:)       !!DOWNWARD LONGWAVE RADIATION
   REAL(SP), ALLOCATABLE,TARGET  :: DSW_AIR(:)       !!DOWNWARD SHORTWAVE RADIATION
# endif
# if defined (ICE)
   REAL(SP), ALLOCATABLE,TARGET  :: QA_AIR(:)        !!SPECIFIC HUMIDITY  (KG/KG)
!   REAL(SP), ALLOCATABLE,TARGET  :: CLOUD(:)         !!Total Cloud Cover
   REAL(SP), ALLOCATABLE,TARGET  :: WSSURF(:)        !!  sea surface salinity    flux
   REAL(SP), ALLOCATABLE,TARGET  :: WHSURF(:)        !!  sea surface heat        flux
   REAL(SP), ALLOCATABLE,TARGET  :: FWtSURF(:)       !!  sea surface fresh water flux

!!============yding=============================================================
   REAL(SP), ALLOCATABLE,TARGET :: QPICE(:)         !!SURFACE ice volume FOR INT
!!============yding end===========================================================

#   if defined (ICE_EMBEDDING)
   REAL(SP), ALLOCATABLE,TARGET :: QFLICE(:)         !!water volume FOR INT

   REAL(SP), ALLOCATABLE,TARGET :: QTHICE(:)         !!ice thickness in the water
   REAL(SP), ALLOCATABLE,TARGET :: mlt_lat(:)         !!ice thickness in the water
   REAL(SP), ALLOCATABLE,TARGET :: QTHICERK(:)       !!ice thickness in the water
   REAL(SP), ALLOCATABLE,TARGET :: QTHICE_incr(:)    !!ice thickness increament
   REAL(SP), ALLOCATABLE,TARGET :: QTHICE1(:)         !!ice thickness in the water

#   endif

# endif
# if defined (AIR_PRESSURE) || (HEATING_CALCULATED)
   REAL(SP), ALLOCATABLE,TARGET  :: PA_AIR(:)        !!SURFACE PRESSURE
# endif
# if defined (ICE) || (HEATING_SOLAR)
   REAL(SP), ALLOCATABLE,TARGET  :: CLOUD(:)         !!Total Cloud Cover
# endif
# if defined (HEATING_SOLAR)
   REAL(SP), ALLOCATABLE,TARGET  :: TDEW(:)         !!Dew Point
# endif

!-----------------------2-d flow fluxes--------------------------------------------!

   REAL(SP), ALLOCATABLE,TARGET :: PSTX(:)           !!EXT MODE BAROTROPIC TERMS
   REAL(SP), ALLOCATABLE,TARGET :: PSTY(:)           !!EXT MODE BAROTROPIC TERMS
   REAL(SP), ALLOCATABLE,TARGET :: ADVUA(:) 
   REAL(SP), ALLOCATABLE,TARGET :: ADVVA(:) 
   REAL(SP), ALLOCATABLE,TARGET :: ADX2D(:) 
   REAL(SP), ALLOCATABLE,TARGET :: ADY2D(:) 
   REAL(SP), ALLOCATABLE,TARGET :: DRX2D(:) 
   REAL(SP), ALLOCATABLE,TARGET :: DRY2D(:) 
   REAL(SP), ALLOCATABLE,TARGET :: TPS(:)            !!WORKING ARRAY             
   REAL(SP), ALLOCATABLE,TARGET :: ADVX(:,:)      
   REAL(SP), ALLOCATABLE,TARGET :: ADVY(:,:)     

!---------------- internal mode   arrays-(element based)----------------------------!

   REAL(SP), ALLOCATABLE,TARGET :: U(:,:)         !X-VELOCITY
   REAL(SP), ALLOCATABLE,TARGET :: V(:,:)         !Y-VELOCITY

   REAL(SP), ALLOCATABLE,TARGET :: UBETA(:,:)     !X-VELOCITY temp time step
   REAL(SP), ALLOCATABLE,TARGET :: VBETA(:,:)     !Y-VELOCITY temp time step

   REAL(SP), ALLOCATABLE :: UBETA2D(:)
   REAL(SP), ALLOCATABLE :: VBETA2D(:)

   REAL(SP), ALLOCATABLE, TARGET :: partition(:) !gwc

   REAL(SP), ALLOCATABLE,TARGET :: W(:,:)         !VERTICAL VELOCITY IN SIGMA SYSTEM
   REAL(SP), ALLOCATABLE,TARGET :: WW(:,:)        !Z-VELOCITY
   REAL(SP), ALLOCATABLE,TARGET :: UF(:,:)        !X-VELOCITY FROM PREVIOUS TIMESTEP
   REAL(SP), ALLOCATABLE,TARGET :: VF(:,:)        !Y-VELOCITY FROM PREVIOUS TIMESTEP
   REAL(SP), ALLOCATABLE,TARGET :: WT(:,:)        !Z-VELOCITY FROM PREVIOUS TIMESTEP
   REAL(SP), ALLOCATABLE,TARGET :: RHO(:,:)       !DENSITY AT ELEMENTS
   REAL(SP), ALLOCATABLE,TARGET :: RMEAN(:,:)     !INITIAL DENSITY AT ELEMENTS
   REAL(SP), ALLOCATABLE,TARGET :: T(:,:)         !TEMPERATURE AT ELEMENTS
   REAL(SP), ALLOCATABLE,TARGET :: TMEAN(:,:)     !INITIAL TEMPERATURE AT ELEMENTS
   REAL(SP), ALLOCATABLE,TARGET :: S(:,:)         !SALINITY AT ELEMENTS
   REAL(SP), ALLOCATABLE,TARGET :: SMEAN(:,:)     !INITIAL SALINITY AT ELEMENTS
   REAL(SP), ALLOCATABLE,TARGET :: Q2(:,:)        !2 X TURBULENT KINETIC ENERGY AT NODES
   REAL(SP), ALLOCATABLE,TARGET :: L(:,:)         !TURBULENT LENGTH MACROSCALE 
   REAL(SP), ALLOCATABLE,TARGET :: Q2L(:,:)       !2 X TURBULENT KE X LENGTH AT NODES
#  if defined (GOTM)
  REAL(SP), ALLOCATABLE,TARGET :: TKE(:,:)       !TURBULENT KINETIC ENERGY AT NODES
  REAL(SP), ALLOCATABLE,TARGET :: TKEF(:,:)      !TURBULENT KINETIC ENERGY AT PREVIOUS TIME
  REAL(SP), ALLOCATABLE,TARGET :: TEPS(:,:)      !TURBULENT DISSIPATION AT NODES
  REAL(SP), ALLOCATABLE,TARGET :: TEPSF(:,:)     !TURBULENT DISSIPATION AT PREVIOUS TIME
#  endif
  REAL(SP), ALLOCATABLE,TARGET :: KM(:,:)        !TURBULENT EDDY VISCOSITY FOR MOMENTUM 
  REAL(SP), ALLOCATABLE,TARGET :: KH(:,:)        !TURBULENT DIFFUSIVITY FOR SALINITY/TEMP 
  REAL(SP), ALLOCATABLE,TARGET :: KQ(:,:)        !TURBULENT DIFFUSIVITY FOR Q2/Q2L 
  REAL(SP), ALLOCATABLE,TARGET :: AAM(:,:)       !STORAGE FOR OUTPUT OF HORIZONTAL VISCOSITY
  REAL(SP), ALLOCATABLE,TARGET :: Q2F(:,:)       !WORKING ARRAY FOR UPDATING Q2
  REAL(SP), ALLOCATABLE,TARGET :: Q2LF(:,:)      !WORKING ARRAY FOR UPDATING Q2F 
  REAL(SP), ALLOCATABLE,TARGET :: KM1(:,:)       !TURBULENT EDDY VISCOSITY FOR MOMENTUM 

  ! VARIABLE HORIZONTAL VISCOSITY COEFFICENTS  
  REAL(SP), ALLOCATABLE :: CC_HVC(:)
  REAL(SP), ALLOCATABLE :: NN_HVC(:)

  !-----------------------3d variable arrays-(node based)-----------------------------!

  REAL(SP), ALLOCATABLE,TARGET :: T1(:,:)         !!TEMPERATURE AT NODES               
  REAL(SP), ALLOCATABLE,TARGET :: S1(:,:)         !!SALINITY AT NODES               
  REAL(SP), ALLOCATABLE,TARGET :: RHO1(:,:)       !!DENSITY AT NODES               
  REAL(SP), ALLOCATABLE,TARGET :: TF1(:,:)        !!TEMPERATURE FROM PREVIOUS TIME
  REAL(SP), ALLOCATABLE,TARGET :: SF1(:,:)        !!SALINITY FROM PREVIOUS TIME 
!J. Ge for tracer advection
  REAL(SP), ALLOCATABLE,TARGET :: T0(:,:)         !!TEMPERATURE AT NODES AT PREVIOUS TIME STEP
  REAL(SP), ALLOCATABLE,TARGET :: T2(:,:)         !!TEMPERATURE AT NODES AT PREVIOUS TWO STEP
  REAL(SP), ALLOCATABLE,TARGET :: S0(:,:)         !!SALINITY AT NODES AT PREVIOUS TIME STEP
  REAL(SP), ALLOCATABLE,TARGET :: S2(:,:)         !!SALINITY AT NODES AT PREVIOUS TWO STEP
!J. Ge for tracer advection
  REAL(SP), ALLOCATABLE,TARGET :: TMEAN1(:,:)     !!MEAN INITIAL TEMP
  REAL(SP), ALLOCATABLE,TARGET :: SMEAN1(:,:)     !!MEAN INITIAL SALINITY 
  REAL(SP), ALLOCATABLE,TARGET :: RMEAN1(:,:)     !!MEAN INITIAL DENSITY 
  REAL(SP), ALLOCATABLE,TARGET :: WTS(:,:)        !!VERTICAL VELOCITY IN SIGMA SYSTEM
  REAL(SP), ALLOCATABLE,TARGET :: WTTS(:,:)       !!WTS FROM PREVIOUS TIMESTEP        

  !---------------------------internal mode fluxes------------------------------------!

  REAL(SP), ALLOCATABLE,TARGET :: DRHOX(:,:)      !!BAROCLINIC PG IN X DIRECTION
  REAL(SP), ALLOCATABLE,TARGET :: DRHOY(:,:)      !!BAROCLINIC PG IN Y DIRECTION

  !------------shape coefficient arrays and control volume metrics--------------------!

  REAL(SP), ALLOCATABLE,TARGET :: A1U(:,:)      
  REAL(SP), ALLOCATABLE,TARGET :: A2U(:,:)     
  REAL(SP), ALLOCATABLE,TARGET :: AWX(:,:)   
  REAL(SP), ALLOCATABLE,TARGET :: AWY(:,:)  
  REAL(SP), ALLOCATABLE,TARGET :: AW0(:,:) 
  REAL(SP), ALLOCATABLE,TARGET :: ALPHA(:)


  !-----salinity and temperature bottom diffusion condition/bottom depth gradients----!

  REAL(SP), ALLOCATABLE,TARGET :: PHPN(:)   
  REAL(SP), ALLOCATABLE,TARGET :: PFPXB(:)   
  REAL(SP), ALLOCATABLE,TARGET :: PFPYB(:)   
  REAL(SP), ALLOCATABLE,TARGET :: SITA_GD(:)  
  REAL(SP), ALLOCATABLE,TARGET :: AH_BOTTOM(:)  

  !-----arrays used for averaging of flow quantities for output-----------------------!

!  REAL(SP), ALLOCATABLE,TARGET :: U_AVE(:,:)   !U AVERAGED OVER INT_AVGE ITERATIONS
!  REAL(SP), ALLOCATABLE,TARGET :: V_AVE(:,:)   !V AVERAGED OVER INT_AVGE ITERATIONS
!  REAL(SP), ALLOCATABLE,TARGET :: W_AVE(:,:)   !WW AVERAGED OVER INT_AVGE ITERATIONS
!  REAL(SP), ALLOCATABLE,TARGET :: KM_AVE(:,:)  !KM AVERAGED OVER INT_AVGE ITERATIONS
!  REAL(SP), ALLOCATABLE,TARGET :: KH_AVE(:,:)  !KH AVERAGED OVER INT_AVGE ITERATIONS
!  REAL(SP), ALLOCATABLE,TARGET :: T_AVE(:,:)   !T1 AVERAGED OVER INT_AVGE ITERATIONS
!  REAL(SP), ALLOCATABLE,TARGET :: S_AVE(:,:)   !S1 AVERAGED OVER INT_AVGE ITERATIONS
!  REAL(SP), ALLOCATABLE,TARGET :: R_AVE(:,:)   !RHO1 AVERAGED OVER INT_AVGE ITERATIONS
!  REAL(SP), ALLOCATABLE,TARGET :: EL_AVE(:)    !EL AVERAGED OVER INT_AVGE ITERATIONS

  !-----arrays used for surface of flow quantities for output-----------------------!

  REAL(SP), ALLOCATABLE,TARGET :: U_SURFACE(:)    !!SURFACE X-VELOC
  REAL(SP), ALLOCATABLE,TARGET :: V_SURFACE(:)    !!SURFACE Y-VELOC
  REAL(SP), ALLOCATABLE,TARGET :: T1_SURFACE(:)   !!SURFACE TEMPERATURE AT NODES               
  REAL(SP), ALLOCATABLE,TARGET :: S1_SURFACE(:)   !!SURFACE SALINITY AT NODES               
  REAL(SP), ALLOCATABLE,TARGET :: VISCOFH_SURFACE(:)
  REAL(SP), ALLOCATABLE,TARGET :: VISCOFM_SURFACE(:)
 
  !---------------------------------------------------------------------------------!
  REAL(SP), ALLOCATABLE,TARGET :: VISCOFH(:,:)
  REAL(SP), ALLOCATABLE,TARGET :: VISCOFM(:,:)

  REAL(SP), ALLOCATABLE,TARGET :: HYW(:,:)

# if defined (EQUI_TIDE)
  INTEGER,  ALLOCATABLE :: EQUI_ID(:)
  REAL*8,   ALLOCATABLE :: EQUI_VN(:)
  REAL*8,   ALLOCATABLE :: EQUI_U(:,:)
  REAL*8,   ALLOCATABLE :: EQUI_F(:,:)

  INTEGER :: month_equi = 0
  TYPE(TIME):: TIME_16DAY

  TYPE :: ASTR
    INTEGER  :: II   = 0
    INTEGER  :: JJ   = 0
    INTEGER  :: KK   = 0
    INTEGER  :: LL   = 0
    INTEGER  :: MM   = 0
    INTEGER  :: NN   = 0
    REAL*8   :: SEMI = 0.0D0
    INTEGER  :: NJ   = 0
    INTEGER, DIMENSION(10) ::   LDEL = 0
    INTEGER, DIMENSION(10) ::   MDEL = 0
    INTEGER, DIMENSION(10) ::   NDEL = 0
    REAL*8,  DIMENSION(10) ::   PH   = 0.0D0
    REAL*8,  DIMENSION(10) ::   EE   = 0.0D0
    INTEGER, DIMENSION(10) ::   IR   = 0
  END TYPE ASTR

  TYPE(ASTR) EQUI_ASTR(8)
# endif

CONTAINS


!==============================================================================|
   SUBROUTINE N2E3D(NVAR,EVAR)
!==============================================================================|
   IMPLICIT NONE 
   REAL(SP), DIMENSION(0:MT,1:KB), INTENT(IN)  :: NVAR  
   REAL(SP), DIMENSION(0:NT,1:KB), INTENT(INOUT) :: EVAR  
   INTEGER I,K
!------------------------------------------------------------------------------|
   DO K=1,KB
     DO I = 1, NT
       EVAR(I,K) = ONE_THIRD*(NVAR(NV(I,1),K)+NVAR(NV(I,2),K)+NVAR(NV(I,3),K))
     END DO
   END DO
   RETURN
   END SUBROUTINE N2E3D 
!==============================================================================|

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!==============================================================================|
   SUBROUTINE N2E2D(NVAR,EVAR)
!==============================================================================|
   IMPLICIT NONE
   REAL(SP), DIMENSION(0:MT), INTENT(IN)  :: NVAR 
   REAL(SP), DIMENSION(0:NT), INTENT(INOUT) :: EVAR  
   INTEGER I,K
!------------------------------------------------------------------------------|
   DO I = 1, NT
     EVAR(I) = ONE_THIRD*(NVAR(NV(I,1))+NVAR(NV(I,2))+NVAR(NV(I,3)))
   END DO
   RETURN
   END SUBROUTINE N2E2D
!==============================================================================|

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!==============================================================================|
   SUBROUTINE E2N2D(EVAR,NVAR)
!==============================================================================|
   IMPLICIT NONE
   REAL(SP), DIMENSION(0:NT), INTENT(IN ) :: EVAR  
   REAL(SP), DIMENSION(0:MT), INTENT(INOUT) :: NVAR 

   INTEGER I,K
!------------------------------------------------------------------------------|
   DO I=1,M
     NVAR(I) = SUM(EVAR(NBVE(I,1:NTVE(I))))/FLOAT(NTVE(I))
   END DO
   RETURN
   END SUBROUTINE E2N2D
!==============================================================================|

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!==============================================================================|
   SUBROUTINE E2N3D(EVAR,NVAR)
!==============================================================================|
   IMPLICIT NONE
   REAL(SP), DIMENSION(0:NT,1:KB), INTENT(IN ) :: EVAR  
   REAL(SP), DIMENSION(0:MT,1:KB), INTENT(INOUT) :: NVAR 
   INTEGER I,K
!------------------------------------------------------------------------------|
   DO K=1,KB
     DO I=1,M
       NVAR(I,K) = SUM(EVAR(NBVE(I,1:NTVE(I)),K))/FLOAT(NTVE(I))
     END DO
   END DO
   RETURN
   END SUBROUTINE E2N3D
!==============================================================================|
   
!==============================================================================|
!    Allocate and Initialize Most Arrays                                       !
!==============================================================================|
   
   SUBROUTINE ALLOC_VARS(dbg_set)
     
!==============================================================================!
   IMPLICIT NONE
   logical,intent(in) :: dbg_set
   INTEGER NCT,NDB
!==============================================================================!
# if !defined (DOUBLE_PRECISION)
   NDB = 1       !!GWC BASE THIS ON KIND
# else
   NDB = 2
# endif

   NCT = NT*3 ! A DIMENSION USED FOR ALLOCATION!
   
!==============================================================================!
!  ALLOCATE:                                                                   !
!==============================================================================!
   
!--------------------------Grid Metrics---------------------------------------------!
   
   ALLOCATE(LON(0:MT))           ;LON  = ZERO   !!LONGITUDE AT THE NODE
   ALLOCATE(LAT(0:MT))           ;LAT  = ZERO   !!LATITUDE AT THE NODE
   ALLOCATE(LONC(0:NT))          ;LONC = ZERO   !!LONGITUDE AT THE NODE
   ALLOCATE(LATC(0:NT))          ;LATC = ZERO   !!LATITUDE AT THE NODE
   ALLOCATE(XM(0:MT))            ;XM   = ZERO   !!X-COORD AT NODE IN METERS
   ALLOCATE(YM(0:MT))            ;YM   = ZERO   !!Y-COORD AT NODE IN METERS
   ALLOCATE(XMC(0:NT))           ;XMC  = ZERO   !!X-COORD AT FACE CENTER IN METERS 
   ALLOCATE(YMC(0:NT))           ;YMC  = ZERO   !!Y-COORD AT FACE CENTER IN METERS
   ALLOCATE(XC(0:NT))            ;XC   = ZERO   !!X-COORD AT FACE CENTER 
   ALLOCATE(YC(0:NT))            ;YC   = ZERO   !!Y-COORD AT FACE CENTER
   ALLOCATE(VX(0:MT))            ;VX   = ZERO   !!X-COORD AT GRID POINT
   ALLOCATE(VY(0:MT))            ;VY   = ZERO   !!Y-COORD AT GRID POINT
   ALLOCATE(ART(0:NT))           ;ART  = ZERO   !!AREA OF ELEMENT
   ALLOCATE(ART1(0:MT))          ;ART1 = ZERO   !!AREA OF NODE-BASE CONTROl VOLUME
   ALLOCATE(ART2(0:MT))          ;ART2 = ZERO   !!AREA OF ELEMENTS AROUND NODE
   ALLOCATE(GRAV_N(0:MT))        ;GRAV_N = zero !! VARIABLE GRAVITY
   ALLOCATE(GRAV_E(0:NT))        ;GRAV_E = zero !! VARIABLE GRAVITY

   MEMCNT = MT*9*NDB + NT*8*NDB + MEMCNT

!----------------Node, Boundary Condition, and Control Volume-----------------------!

! ALLOCATED IN setup_domain.F But count the memory here
!   ALLOCATE(NV(0:NT,4))          ;NV       = 0  !!NODE NUMBERING FOR ELEMENTS
!   ALLOCATE(NVGL(0:NT,3))        ;NVGL     = 0  !!GLOBAL NODE NUMBERING FOR LOCAL ELEMENTS
   ALLOCATE(NBE(0:NT,3))         ;NBE      = 0 !!INDICES OF ELMNT NEIGHBORS
   ALLOCATE(NTVE(0:MT))          ;NTVE     = 0 
   ALLOCATE(NTSN(0:MT))          ;NTSN     = 0 
   ALLOCATE(ISONB(0:MT))         ;ISONB    = 0  !!NODE MARKER = 0,1,2
   ALLOCATE(ISONB_W(0:MT))       ;ISONB_W  = 0  !!NODE MARKER = 0,1,2
   ALLOCATE(ISBCE(0:NT))         ;ISBCE    = 0 
   ALLOCATE(NIEC(NCT,2))         ;NIEC     = 0
   ALLOCATE(NTRG(NCT))           ;NTRG     = 0

   ! POSITION OF NODAL CONTROL VOLUME CORNERS 
   ALLOCATE(XIJE(NCT,2))         ;XIJE     = ZERO
   ALLOCATE(YIJE(NCT,2))         ;YIJE     = ZERO 
   ! LENGTH OF NODAL CONTROL VOLUME EDGES
   ALLOCATE(DLTXE(NCT))          ;DLTXE    = ZERO
   ALLOCATE(DLTYE(NCT))          ;DLTYE    = ZERO
   ALLOCATE(DLTXYE(NCT))         ;DLTXYE   = ZERO !! TOTAL LENGTH
   ALLOCATE(SITAE(NCT))          ;SITAE    = ZERO !! ANGLE

   ! LENGTH BETWEEN NODE AND CONTROL VOLUMEN EDGE CENTER
   ALLOCATE(DLTXNCVE(NCT,2))     ; DLTXNCVE   = ZERO  !! DeLTa X Node to Control Volume Edge 
   ALLOCATE(DLTYNCVE(NCT,2))     ; DLTYNCVE   = ZERO  !! DeLTa Y Node to Control Volume Edge 


   ! THE FOLLOWING ARRAYS COULD BE REPLACED WITH (N,3) ARRAYS
   ! BUT THE INDEXING WOULD BE COMPLEX AND THERE ARE SIGN ISSUES!

   ! TRIANGLE EDGE LENGTH FOR EDGES SURROUNDING EACH NODE
   ! NTSN Can not be greater than 13!
   ALLOCATE(DLTXTRIE(M,12))         ;DLTXTRIE   = ZERO !! DeLTa X TRIangle Edge
   ALLOCATE(DLTYTRIE(M,12))         ;DLTYTRIE   = ZERO !! DeLTa Y TRIangle Edge

   ALLOCATE(DLTXNEC(M,12))         ;DLTXNEC   = ZERO !! DeLTa X Node to Edge Center
   ALLOCATE(DLTYNEC(M,12))         ;DLTYNEC   = ZERO !! DeLTa Y Node to Edge Center

   ! DISTANCE BETWEEN TRIANGLE EDGE CENTERS FOR EACH TRIANGLE AROUND A NODE
   ! NTVE Can not be greater than 13!  
   ALLOCATE(DLTXECEC(M,12))         ;DLTXECEC   = ZERO  !! DeLTa X Edge Center to Edge Center
   ALLOCATE(DLTYECEC(M,12))         ;DLTYECEC   = ZERO  !! DeLTa Y Edge Center to Edge Center


   
   MEMCNT = NT*4 + MT*3 +M*6*12*NDB  + NCT*3 + NCT*12*NDB  + MEMCNT

!----------------2-d arrays for the general vertical coordinate -------------------------------!

   ALLOCATE(Z(0:MT,KB))               ; Z      = ZERO    !!SIGMA COORDINATE VALUE 
   ALLOCATE(ZZ(0:MT,KB))              ; ZZ     = ZERO    !!INTRA LEVEL SIGMA VALUE
   ALLOCATE(DZ(0:MT,KB))              ; DZ     = ZERO    !!DELTA-SIGMA VALUE
   ALLOCATE(DZZ(0:MT,KB))             ; DZZ    = ZERO    !!DELTA OF INTRA LEVEL SIGMA 
   ALLOCATE(Z1(0:NT,KB))              ; Z1     = ZERO    !!SIGMA COORDINATE VALUE 
   ALLOCATE(ZZ1(0:NT,KB))             ; ZZ1    = ZERO    !!INTRA LEVEL SIGMA VALUE
   ALLOCATE(DZ1(0:NT,KB))             ; DZ1    = ZERO    !!DELTA-SIGMA VALUE
   ALLOCATE(DZZ1(0:NT,KB))            ; DZZ1   = ZERO    !!DELTA OF INTRA LEVEL SIGMA 
   MEMCNT = MT*KB*4*NDB + NT*KB*4*NDB +MEMCNT

!---------------2-d flow variable arrays at elements-------------------------------!

   ALLOCATE(UA(0:NT))            ;UA        = ZERO  !!VERTICALLY AVERAGED X-VELOC
   ALLOCATE(VA(0:NT))            ;VA        = ZERO  !!VERTICALLY AVERAGED Y-VELOC
   ALLOCATE(UAF(0:NT))           ;UAF       = ZERO  !!UA FROM PREVIOUS RK STAGE 
   ALLOCATE(VAF(0:NT))           ;VAF       = ZERO  !!VA FROM PREVIOUS RK STAGE 
# if !defined (SEMI_IMPLICIT)
   ALLOCATE(UARK(0:NT))          ;UARK      = ZERO  !!UA FROM PREVIOUS TIMESTEP 
   ALLOCATE(VARK(0:NT))          ;VARK      = ZERO  !!VA FROM PREVIOUS TIMESTEP 
   ALLOCATE(UARD(0:NT))          ;UARD      = ZERO  !!UA AVERAGED OVER EXTERNAL INT
   ALLOCATE(VARD(0:NT))          ;VARD      = ZERO  !!VA AVERAGED OVER EXTERNAL INT
# endif
   ALLOCATE(COR(0:NT))           ;COR       = ZERO  !!CORIOLIS PARAMETER
   ALLOCATE(F_ALFA(0:NT))        ;F_ALFA    = 1.0_SP  !!EQUATORIAL BETA PLANE PARAMETER
   ALLOCATE(H1(0:NT))            ;H1        = ZERO  !!BATHYMETRIC DEPTH   
   ALLOCATE(D1(0:NT))            ;D1        = ZERO  !!DEPTH
   ALLOCATE(DT1(0:NT))           ;DT1       = ZERO  !!DEPTH  
   ALLOCATE(EL1(0:NT))           ;EL1       = ZERO  !!SURFACE ELEVATION
   ALLOCATE(ELF1(0:NT))          ;ELF1      = ZERO  !!SURFACE ELEVATION
   ALLOCATE(DTFA(0:MT))          ;DTFA      = ZERO  !!ADJUSTED DEPTH FOR MASS CONSERVATION
   ALLOCATE(ET1(0:NT))           ;ET1       = ZERO  !!SURFACE ELEVATION
# if !defined (SEMI_IMPLICIT)
   ALLOCATE(ELRK1(0:NT))         ;ELRK1     = ZERO  !!SURFACE ELEVATION
# endif
   ALLOCATE(CC_SPONGE(0:NT))     ;CC_SPONGE = ZERO  !!SPONGE DAMPING COEFFICIENT FOR MOMENTUM
                 MEMCNT = NT*17*NDB + MT*NDB + MEMCNT

!---------------2-d flow variable arrays at nodes----------------------------------!

   ALLOCATE(H(0:MT))             ;H    = ZERO       !!BATHYMETRIC DEPTH   
   ALLOCATE(D(0:MT))             ;D    = ZERO       !!DEPTH   
   ALLOCATE(DT(0:MT))            ;DT   = ZERO       !!DEPTH   
   ALLOCATE(EL(0:MT))            ;EL   = ZERO       !!SURFACE ELEVATION
   ALLOCATE(ELF(0:MT))           ;ELF  = ZERO       !!SURFACE ELEVATION
   ALLOCATE(ET(0:MT))            ;ET   = ZERO       !!SURFACE ELEVATION
   ALLOCATE(EGF(0:MT))           ;EGF  = ZERO       !!SURFACE ELEVATION
# if !defined (SEMI_IMPLICIT)
   ALLOCATE(ELRK(0:MT))          ;ELRK = ZERO       !!SURFACE ELEVATION
# endif
               MEMCNT = MT*8*NDB + MEMCNT

# if defined (EQUI_TIDE)
   ALLOCATE(ELF_EQI(0:MT))       ; ELF_EQI  = ZERO
# if !defined (SEMI_IMPLICIT)
   ALLOCATE(ELRK_EQI(0:MT))      ; ELRK_EQI = ZERO
# endif
   ALLOCATE(EL_EQI(0:MT))        ; EL_EQI   = ZERO
   ALLOCATE(EGF_EQI(0:MT))       ; EGF_EQI  = ZERO
               MEMCNT = MT*4*NDB + MEMCNT
# endif              

# if defined (ATMO_TIDE)
   ALLOCATE(ELF_ATMO(0:MT))       ; ELF_ATMO  = ZERO
# if !defined (SEMI_IMPLICIT)
   ALLOCATE(ELRK_ATMO(0:MT))      ; ELRK_ATMO = ZERO
# endif
   ALLOCATE(EL_ATMO(0:MT))        ; EL_ATMO   = ZERO
   ALLOCATE(EGF_ATMO(0:MT))       ; EGF_ATMO  = ZERO
               MEMCNT = MT*4*NDB + MEMCNT
# endif              
   
# if defined (AIR_PRESSURE)
   ALLOCATE(ELF_AIR(0:MT))       ; ELF_AIR  = ZERO
# if !defined (SEMI_IMPLICIT)
   ALLOCATE(ELRK_AIR(0:MT))      ; ELRK_AIR = ZERO
# endif
   ALLOCATE(EL_AIR(0:MT))        ; EL_AIR   = ZERO
   ALLOCATE(EGF_AIR(0:MT))       ; EGF_AIR  = ZERO
               MEMCNT = MT*4*NDB + MEMCNT
# endif  

   ALLOCATE(VORT(0:MT))          ; VORT     = ZERO            
               MEMCNT = MT*NDB + MEMCNT

!---------------surface/bottom/edge boundary conditions-----------------------------!

   ALLOCATE(CBC(0:NT))           ;CBC    = ZERO     !!BOTTOM FRICTION     
   ALLOCATE(CC_Z0B(0:NT))        ;CC_Z0B = ZERO     !!BOTTOM ROUGHNESS VARIABLE     
!# if !defined (SEMI_IMPLICIT)
   ALLOCATE(WUSURF2(0:NT))       ;WUSURF2= ZERO     !!SURFACE FRICTION FOR EXT
   ALLOCATE(WVSURF2(0:NT))       ;WVSURF2= ZERO     !!SURFACE FRICTION FOR EXT
!# endif
   ALLOCATE(WUBOT(0:NT))         ;WUBOT  = ZERO     !!BOTTOM FRICTION
   ALLOCATE(WVBOT(0:NT))         ;WVBOT  = ZERO     !!BOTTOM FRICTION
   ALLOCATE(TAUBM(0:NT))         ;TAUBM  = ZERO     !!BOTTOM FRICTION
   ALLOCATE(WUBOT_N(0:MT))       ;WUBOT_N  = ZERO   !!U-Component bottom shear stress on nodes  
   ALLOCATE(WVBOT_N(0:MT))       ;WVBOT_N  = ZERO   !!V-Component bottom shear stress on nodes
   ALLOCATE(TAUBM_N(0:MT))       ;TAUBM_N  = ZERO   !!Magnitude bottom shear stress on nodes
   ALLOCATE(WUSURF(0:NT))        ;WUSURF = ZERO     !!SURFACE FRICTION FOR INT
   ALLOCATE(WVSURF(0:NT))        ;WVSURF = ZERO     !!SURFACE FRICTION FOR INT
   ALLOCATE(WUSURF_save(0:NT))   ;WUSURF_save = ZERO!!SURFACE FRICTION FOR INT
   ALLOCATE(WVSURF_save(0:NT))   ;WVSURF_save = ZERO!!SURFACE FRICTION FOR INT
   ALLOCATE(UUWIND(0:NT))        ;UUWIND = ZERO     !!SURFACE X-WIND 
   ALLOCATE(VVWIND(0:NT))        ;VVWIND = ZERO     !!SURFACE Y-WIND 
   ALLOCATE(SWRAD(0:MT))         ;SWRAD  = ZERO     !!SURFACE INCIDENT RADIATION
   ALLOCATE(WTSURF(0:MT))        ;WTSURF = ZERO 
   ALLOCATE(SWRAD_WATTS(0:MT))  ;SWRAD_WATTS  = ZERO !!SURFACE INCIDENT RADIATION
   ALLOCATE(WTSURF_WATTS(0:MT)) ;WTSURF_WATTS = ZERO 
# if defined (HEATING_CALCULATED) || (HEATING_SOLAR) || (ICE)
   !EJA edit to add SENSIBLE, LATENT, LONGWAVE HEAT FLUX
   ALLOCATE(HSENS_WATTS(0:MT))  ;HSENS_WATTS  = ZERO !!SENSIBLE HEAT FLUX
   ALLOCATE(HLAT_WATTS(0:MT))   ;HLAT_WATTS   = ZERO !!LATENT HEAT FLUX
   ALLOCATE(LWRAD_WATTS(0:MT))  ;LWRAD_WATTS  = ZERO !!LONGWAVE RADIATION
#  endif

# if !defined (SEMI_IMPLICIT)
   ALLOCATE(QPREC2(0:MT))        ;QPREC2 = ZERO     !!SURFACE PRECIPITATION FOR EXT
   ALLOCATE(QEVAP2(0:MT))        ;QEVAP2 = ZERO     !!SURFACE EVAPORATION FOR EXT
# endif
   ALLOCATE(QPREC(0:MT))         ;QPREC = ZERO     !!SURFACE PRECIPITATION FOR INT
   ALLOCATE(QEVAP(0:MT))         ;QEVAP = ZERO     !!SURFACE EVAPORATION FOR INT

# if defined (WAVE_CURRENT_INTERACTION) && (WAVE_OFFLINE)
   ALLOCATE(WHS(0:MT))           ;WHS      = ZERO
   ALLOCATE(WDIR(0:MT))          ;WDIR     = ZERO
   ALLOCATE(WPER(0:MT))          ;WPER     = ZERO
   ALLOCATE(WLENGTH(0:MT))       ;WLENGTH  = ZERO
   ALLOCATE(WPER_BOT(0:MT))      ;WPER_BOT = ZERO
   ALLOCATE(WUB_BOT(0:MT))       ;WUB_BOT  = ZERO
# endif

# if defined (OFFLINE_SEDIMENT) || (OFFLINE_BIOLOGY)
   ALLOCATE(OFFLINE_U(0:NT,KB) )    ;OFFLINE_U   = ZERO
   ALLOCATE(OFFLINE_V(0:NT,KB) )    ;OFFLINE_V   = ZERO
   ALLOCATE(OFFLINE_W(0:MT,KB) )    ;OFFLINE_W   = ZERO
   ALLOCATE(OFFLINE_S1(0:MT,KB))    ;OFFLINE_S1  = ZERO
   ALLOCATE(OFFLINE_T1(0:MT,KB))    ;OFFLINE_T1  = ZERO
   ALLOCATE(OFFLINE_EL(0:MT)   )    ;OFFLINE_EL  = ZERO
   ALLOCATE(OFFLINE_KH(0:MT,KB))    ;OFFLINE_KH  = ZERO
   ALLOCATE(OFFLINE_WN(0:MT)   )    ;OFFLINE_WN  = ZERO
   ALLOCATE(OFFLINE_WC(0:NT)   )    ;OFFLINE_WC  = ZERO

# if defined (GOTM)
   ALLOCATE(OFFLINE_TEPS(0:MT,KB))  ;OFFLINE_TEPS= ZERO
# else
   ALLOCATE(OFFLINE_Q2(0:MT,KB))    ;OFFLINE_Q2  = ZERO
   ALLOCATE(OFFLINE_Q2L(0:MT,KB))   ;OFFLINE_Q2L = ZERO
# endif
# endif

   MEMCNT = NT*10*NDB + MT*8*NDB + MEMCNT

   IF (ICING_MODEL) THEN
      ALLOCATE(ICING_WNDX(0:MT)) ;ICING_WNDX = ZERO
      ALLOCATE(ICING_WNDY(0:MT)) ;ICING_WNDY = ZERO
      ALLOCATE(ICING_SATMP(0:MT));ICING_SATMP  = ZERO
      ALLOCATE(ICING_0kts(0:MT)) ;ICING_0kts = ZERO
      ALLOCATE(ICING_10kts(0:MT));ICING_10kts= ZERO

      MEMCNT = MEMCNT + MT*4*NDB
   END IF
      
# if defined (ICE) || (HEATING_CALCULATED) || (HEATING_SOLAR)
   ALLOCATE(T_AIR(0:MT))        !!BULK AIR TEMPERATURE
   MEMCNT =  MT*NDB + MEMCNT
# endif
# if defined (ICE) || (HEATING_CALCULATED)
   ALLOCATE(RH_AIR(0:MT))       !!RELATIVE HUMIDITY
   ALLOCATE(DLW_AIR(0:MT))       !!DOWNWARD LONGWAVE RADIATION
   ALLOCATE(DSW_AIR(0:MT))       !!DOWNWARD SHORTWAVE RADIATION
   MEMCNT =  MT*3*NDB + MEMCNT
# endif
# if defined (ICE)
   ALLOCATE(QA_AIR(0:MT))        !!SPECIFIC HUMIDITY
!   ALLOCATE(CLOUD(0:MT))         !!Total Cloud Cover
   ALLOCATE(WSSURF(0:MT))        ;WSSURF = Zero      !  sea surface salinity    flux
   ALLOCATE(WHSURF(0:MT))        ;WHSURF = ZERO      !  sea surface heat        flux
   ALLOCATE(FWtSURF(0:MT))       ;FWtSURF= ZERO      !  sea surface fresh water flux
   MEMCNT =  MT*4*NDB + MEMCNT

!!yding
   ALLOCATE(QPICE(0:MT))         ;QPICE = ZERO      !!SURFACE ice volume FOR INT
!!yding
#   if defined (ICE_EMBEDDING)
   ALLOCATE(QFLICE(0:MT))         ;QFLICE = ZERO      !!SURFACE ice volume FOR INT
   ALLOCATE(QTHICE(0:MT))         ;QTHICE = ZERO      !!ice thickness for current step
   ALLOCATE(mlt_lat(0:MT))         ;mlt_lat = ZERO      !!ice thickness for current step
   ALLOCATE(QTHICERK(0:MT))       ;QTHICERK = ZERO    !!ice thickness for previous step
   ALLOCATE(QTHICE_incr(0:MT))    ;QTHICE_incr = ZERO    !!ice thickness increament
   ALLOCATE(QTHICE1(0:NT))         ;QTHICE1 = ZERO       !!ice thickness for current step
#  endif

# endif
# if defined (AIR_PRESSURE) || (HEATING_CALCULATED)
   ALLOCATE(PA_AIR(0:MT))        !!SURFACE PRESSURE
   MEMCNT =  MT*NDB + MEMCNT
# endif
# if defined (ICE) || (HEATING_SOLAR)
   ALLOCATE(CLOUD(0:MT))         !!Total Cloud Cover
   MEMCNT =  MT*NDB + MEMCNT
# endif
# if defined (HEATING_SOLAR)
   ALLOCATE(TDEW(0:MT))         !!Dew Point
   MEMCNT =  MT*NDB + MEMCNT
# endif

!--------------------------------------------------------------
!--------------------------------------------------------------

   ALLOCATE(BFWDIS(0:MT))        ;BFWDIS = ZERO     !!GROUNDWATER FLUX FOR INT
# if !defined (SEMI_IMPLICIT)
   ALLOCATE(BFWDIS2(0:MT))       ;BFWDIS2= ZERO     !!GROUNDWATER FLUX FOR EXT
# endif
   ALLOCATE(BFWSLT(0:MT))        ;BFWSLT = ZERO     !!GROUNDWATER SALT AT CURRENT TIME
   ALLOCATE(BFWTMP(0:MT))        ;BFWTMP = ZERO     !!GROUNDWATER TEMP AT CURRENT TIME

   MEMCNT =  MT*4*NDB + MEMCNT

!-----------------------2-d flow fluxes---------------------------------------------!

   ALLOCATE(PSTX(0:NT))          ;PSTX  = ZERO       !!EXT MODE BAROTROPIC TERMS
   ALLOCATE(PSTY(0:NT))          ;PSTY  = ZERO       !!EXT MODE BAROTROPIC TERMS
   ALLOCATE(ADVUA(0:NT))         ;ADVUA = ZERO 
   ALLOCATE(ADVVA(0:NT))         ;ADVVA = ZERO
   ALLOCATE(ADX2D(0:NT))         ;ADX2D = ZERO
   ALLOCATE(ADY2D(0:NT))         ;ADY2D = ZERO
   ALLOCATE(DRX2D(0:NT))         ;DRX2D = ZERO
   ALLOCATE(DRY2D(0:NT))         ;DRY2D = ZERO
   ALLOCATE(ADVX(0:NT,KB))       ;ADVX  = ZERO 
   ALLOCATE(ADVY(0:NT,KB))       ;ADVY  = ZERO 
   ALLOCATE(TPS(0:NT))           ;TPS   = ZERO      !!WORKING ARRAY             
   MEMCNT = NT*9*NDB + NT*KB*2*NDB + MEMCNT


!---------------- internal mode   arrays-(element based)----------------------------!

   ALLOCATE(U(0:NT,KB))          ;U     = ZERO   !!X-VELOCITY
   ALLOCATE(V(0:NT,KB))          ;V     = ZERO   !!Y-VELOCITY
   ALLOCATE(UBETA(0:NT,KB))      ;UBETA = ZERO   !!X-VELOCITY temp time step
   ALLOCATE(VBETA(0:NT,KB))      ;VBETA = ZERO   !!X-VELOCITY temp time step
   ALLOCATE(UBETA2D(0:NT))       ;UBETA2D = ZERO
   ALLOCATE(VBETA2D(0:NT))       ;VBETA2D = ZERO

   ALLOCATE(W(0:NT,KB))          ;W     = ZERO   !!VERTICAL VELOCITY IN SIGMA SYSTEM
   ALLOCATE(WW(0:NT,KB))         ;WW    = ZERO   !!Z-VELOCITY
   ALLOCATE(UF(0:NT,KB))         ;UF    = ZERO   !!X-VELOCITY FROM PREVIOUS TIMESTEP
   ALLOCATE(VF(0:NT,KB))         ;VF    = ZERO   !!Y-VELOCITY FROM PREVIOUS TIMESTEP
   ALLOCATE(WT(0:NT,KB))         ;WT    = ZERO   !!Z-VELOCITY FROM PREVIOUS TIMESTEP
   ALLOCATE(RHO(0:NT,KB))        ;RHO   = ZERO   !!DENSITY AT ELEMENTS
   ALLOCATE(RMEAN(0:NT,KB))      ;RMEAN = ZERO   !!MEAN INITIAL DENSITY AT ELEMENTS
   ALLOCATE(T(0:NT,KB))          ;T     = ZERO   !!TEMPERATURE AT ELEMENTS
   ALLOCATE(TMEAN(0:NT,KB))      ;TMEAN = ZERO   !!MEAN INITIAL TEMPERATURE AT ELEMENTS
   ALLOCATE(S(0:NT,KB))          ;S     = ZERO   !!SALINITY AT ELEMENTS
   ALLOCATE(SMEAN(0:NT,KB))      ;SMEAN = ZERO   !!MEAN INITIAL SALINITY AT ELEMENTS
               MEMCNT = NT*KB*13*NDB + MEMCNT

!-----------------------3d variable arrays-(node based)-----------------------------!

   ALLOCATE(T1(0:MT,KB))         ;T1     = ZERO  !!TEMPERATURE AT NODES               
   ALLOCATE(S1(0:MT,KB))         ;S1     = ZERO  !!SALINITY AT NODES               
!J. Ge for tracer advection
   ALLOCATE(T0(0:MT,KB))         ;T0    = ZERO   !!TEMPERATURE FROM PREVIOUS TIME STEP 
   ALLOCATE(T2(0:MT,KB))         ;T2    = ZERO   !!TEMPERATURE FROM PREVIOUS TWO STEP
   ALLOCATE(S0(0:MT,KB))         ;S0    = ZERO   !!SALINITY FROM PREVIOUS TIME STEP 
   ALLOCATE(S2(0:MT,KB))         ;S2    = ZERO   !!SALINITY FROM PREVIOUS TWO STEP 
!J. Ge for tracer advection
   ALLOCATE(RHO1(0:MT,KB))       ;RHO1   = ZERO  !!DENSITY AT NODES               
   ALLOCATE(TF1(0:MT,KB))        ;TF1    = ZERO  !!TEMPERATURE FROM PREVIOUS TIME
   ALLOCATE(SF1(0:MT,KB))        ;SF1    = ZERO  !!SALINITY FROM PREVIOUS TIME 
   ALLOCATE(TMEAN1(0:MT,KB))     ;TMEAN1 = ZERO  !!MEAN INITIAL TEMP
   ALLOCATE(SMEAN1(0:MT,KB))     ;SMEAN1 = ZERO  !!MEAN INITIAL SALINITY 
   ALLOCATE(RMEAN1(0:MT,KB))     ;RMEAN1 = ZERO  !!MEAN INITIAL DENSITY 
   ALLOCATE(WTS(0:MT,KB))        ;WTS    = ZERO  !!VERTICAL VELOCITY IN SIGMA SYSTEM
   ALLOCATE(WTTS(0:MT,KB))       ;WTTS   = ZERO  !!WTS FROM PREVIOUS TIMESTEP        
   ALLOCATE(Q2(0:MT,KB))         ;Q2    = ZERO   !!TURBULENT KINETIC ENERGY AT NODES
   ALLOCATE(Q2L(0:MT,KB))        ;Q2L   = ZERO   !!TURBULENT KE*LENGTH AT NODES
   ALLOCATE(L(0:MT,KB))          ;L     = ZERO   !!TURBULENT LENGTH SCALE AT ELEMENTS
   ALLOCATE(KM(0:MT,KB))         ;KM    = ZERO   !!TURBULENT QUANTITY
   ALLOCATE(KH(0:MT,KB))         ;KH    = ZERO   !!TURBULENT QUANTITY
   ALLOCATE(KQ(0:MT,KB))         ;KQ    = ZERO   !!TURBULENT QUANTITY
   ALLOCATE(AAM(0:MT,KB))        ;AAM   = ZERO   !!??
   ALLOCATE(KM1(0:NT,KB))        ;KM1   = ZERO   !!TURBULENT QUANTITY AT ELEMENTS


   ALLOCATE(CC_HVC(0:NT))        ;CC_HVC = ZERO  !!VISCOSITY COEFFICIENT AT ELEMENTS
   ALLOCATE(NN_HVC(0:MT))        ;NN_HVC = ZERO  !!VISCOSITY COEFFICIENT AT NODES

   MEMCNT = MT*KB*18*NDB + NT*KB*NDB + MEMCNT + (MT+NT)*NDB

#  if defined (GOTM)
   ALLOCATE(TKE(0:MT,KB))        ;TKE   = ZERO   !!TURBULENT KINETIC ENERGY AT ELEMENTS
   ALLOCATE(TKEF(0:MT,KB))       ;TKEF  = ZERO   !!TURBULENT KINETIC ENERGY AT ELEMENTS
   ALLOCATE(TEPS(0:MT,KB))       ;TEPS  = ZERO   !!TURBULENT KINETIC ENERGY AT ELEMENTS
   ALLOCATE(TEPSF(0:MT,KB))      ;TEPSF = ZERO   !!TURBULENT KE*LENGTH AT ELEMENTS
   MEMCNT = MT*KB*4*NDB + MEMCNT
#  endif
  
!---------------------------internal mode fluxes------------------------------------!

   ALLOCATE(DRHOX(0:NT,KB))      ;DRHOX  = ZERO 
   ALLOCATE(DRHOY(0:NT,KB))      ;DRHOY  = ZERO 
   ALLOCATE(Q2F(0:MT,KB))        ;Q2F    = ZERO 
   ALLOCATE(Q2LF(0:MT,KB))       ;Q2LF   = ZERO
   MEMCNT = NT*KB*2*NDB + MT*KB*2*NDB + MEMCNT

!------------shape coefficient arrays and control volume metrics--------------------!

   ALLOCATE(A1U(0:NT,4))         ;A1U   = ZERO
   ALLOCATE(A2U(0:NT,4))         ;A2U   = ZERO 
   ALLOCATE(AWX(0:NT,3))         ;AWX   = ZERO 
   ALLOCATE(AWY(0:NT,3))         ;AWY   = ZERO 
   ALLOCATE(AW0(0:NT,3))         ;AW0   = ZERO 
   ALLOCATE(ALPHA(0:NT))         ;ALPHA = ZERO
   MEMCNT = NT*4*2*NDB + NT*3*3*NDB + NT*NDB + MEMCNT

!-----salinity and temperature bottom diffusion condition/bottom depth gradients----!

   ALLOCATE(PHPN(0:MT))          ;PHPN      = ZERO 
   ALLOCATE(PFPXB(MT))           ;PFPXB     = ZERO
   ALLOCATE(PFPYB(MT))           ;PFPYB     = ZERO
   ALLOCATE(SITA_GD(0:MT))       ;SITA_GD   = ZERO 
   ALLOCATE(AH_BOTTOM(MT))       ;AH_BOTTOM = ZERO 
   MEMCNT = MT*5*NDB + MEMCNT

   ALLOCATE(VISCOFH(0:MT,KB))    ;VISCOFH = ZERO
   ALLOCATE(VISCOFM(0:NT,KB))    ;VISCOFM = ZERO
   MEMCNT = MT*KB*NDB + NT*KB*NDB + MEMCNT

  !-----arrays used for surface of flow quantities for output-----------------------!

   ALLOCATE(U_SURFACE(0:NT))          ;U_SURFACE     = ZERO   !!X-VELOCITY
   ALLOCATE(V_SURFACE(0:NT))          ;V_SURFACE     = ZERO   !!Y-VELOCITY
   ALLOCATE(T1_SURFACE(0:MT))         ;T1_SURFACE    = ZERO   !!SURFACE TEMPERATURE AT NODES               
   ALLOCATE(S1_SURFACE(0:MT))         ;S1_SURFACE    = ZERO   !!SURFACE SALINITY AT NODES               
   ALLOCATE(VISCOFH_SURFACE(0:MT))    ;VISCOFH_SURFACE = ZERO
   ALLOCATE(VISCOFM_SURFACE(0:NT))    ;VISCOFM_SURFACE = ZERO
   MEMCNT = MT*3*NDB + NT*3*NDB + MEMCNT


   ALLOCATE(HYW(0:MT,KB))        ;HYW     = ZERO
!-----special initialization which probably do nothing------------------------------!

   DT1(0)   = 100.0_SP

!---------------report approximate memory usage-------------------------------------!


   RETURN
   END SUBROUTINE ALLOC_VARS
!==============================================================================|


END MODULE ALL_VARS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|

MODULE BCS
   USE MOD_TYPES
   USE MOD_PREC
   IMPLICIT NONE
   SAVE

!----------------boundary conditions: Julian tidal forcing--------------------------!

    TYPE(BC)                 :: ELO_TM        !!TIME MAP FOR SURFACE ELEVATION DATA


!----------------GLobal Boundary Condition Information ----------------------------! 


   INTEGER,  ALLOCATABLE :: I_OBC_GL(:)      !!GLOBAL ID OF OPEN BOUNDARY NODES
   INTEGER,  ALLOCATABLE :: I_OBC_GL_W(:)    !!GLOBAL ID OF OPEN BOUNDARY NODES FOR WAVE
   INTEGER               :: IOBCN_GL         !!LOCAL NUMBER OF OPEN BOUNDARY NODES FOR FVCOM
   INTEGER               :: IOBCN_GL_W       !!LOCAL NUMBER OF OPEN BOUNDARY NODES FOR WAVE
   INTEGER               :: IOBCN            !!LOCAL NUMBER OF OPEN BOUNDARY NODES FOR FVCOM
   INTEGER               :: IOBCN_W          !!LOCAL NUMBER OF OPEN BOUNDARY NODES FOR SWAVE
   INTEGER,  ALLOCATABLE :: I_OBC_N(:)       !!OPEN BOUNDARY NODE LIST FOR FVCOM
   INTEGER,  ALLOCATABLE :: I_OBC_N_W(:)     !!OPEN BOUNDARY NODE LIST FOR SWAVE
   INTEGER,  ALLOCATABLE :: I_OBC_N_OUTPUT(:)!!LIST OF LOCAL OBC GLOBAL NODES FOR OUTPUT
   INTEGER,  ALLOCATABLE :: TYPE_OBC_GL(:)   !!OUTER BOUNDARY NODE TYPE (FOR SURFACE ELEVATION)
   INTEGER,  ALLOCATABLE :: TYPE_OBC(:)      !!OUTER BOUNDARY NODE TYPE (FOR SURFACE ELEVATION)
   INTEGER           ::  OBC_NTIME

!----------------boundary conditions: ground water----------------------------------!

  INTEGER, ALLOCATABLE :: NODE_BFW(:)      !!LOCAL GROUNDWATER NODES
  INTEGER, ALLOCATABLE :: BFW_GL2LOC(:)    !!GLOBAL TO LOCAL MAPPING OF GWATER NODES
  REAL(SP),  ALLOCATABLE :: BFWQDIS(:,:)   !!GROUNDWATER FRESH WATER FLUX DATA
!  TYPE(BC)      :: BFW_TM                  !!TIME MAP FOR GROUNDWATER DATA

!----------------boundary conditions: spectral tidal forcing----------------------!
   INTEGER :: nTideComps
   REAL(SP),  ALLOCATABLE :: PERIOD(:)     !!TIDE PERIOD
   REAL(SP),  ALLOCATABLE :: APT(:,:)      !!TIDE AMPLITUDE
   REAL(SP),  ALLOCATABLE :: PHAI(:,:)     !!TIDE PHASE 
   REAL(SP),  ALLOCATABLE :: EMEAN(:)      !!MEAN SURFACE ELEVATION

   ! ONLY USED IF COMPILED WITH EQUI_TIDE
   REAL(SP),  ALLOCATABLE :: APT_EQI(:)    !! EQUILIBIRUIM TIDE AMPLITUDE
    REAL(SP),  ALLOCATABLE :: BETA_EQI(:)   !! EQUILIBIRUIM TIDE LOVE NUMBER
   CHARACTER(LEN=80), ALLOCATABLE :: TIDE_TYPE(:)    !! EQUILIBIRUIM TIDE AMPLITUDE

   CHARACTER(LEN=12), PARAMETER :: DIURNAL="DIURNAL"
   CHARACTER(LEN=12), PARAMETER :: SEMIDIURNAL="SEMIDIURNAL"


!-- Old Tidal Periods before they became part of the forcing file ---------:    
!              s2         m2           n2          k1          p1         o1 
!PERIOD = (/43200.0_SP, 44712.0_SP, 45570.0_SP, 86164.0_SP, 86637.0_SP, 92950.0_SP/)


!----------------boundary conditions: Julian tidal forcing--------------------------!

   REAL(SP), ALLOCATABLE    :: ELSBC(:,:)    !!INPUT SURFACE ELEVATION

END MODULE  BCS
