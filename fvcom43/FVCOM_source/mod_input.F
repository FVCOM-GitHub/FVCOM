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

MODULE MOD_INPUT
  USE MOD_NCTOOLS
  USE MOD_UTILS
  IMPLICIT NONE

  PUBLIC
  SAVE

  TYPE(NCFILE), POINTER ::NC_DAT
  TYPE(NCFILE), POINTER ::NC_AVG
  TYPE(NCFILE), POINTER ::NC_RST
  TYPE(NCFILE), POINTER ::NC_START
  TYPE(NCFILE), POINTER ::NC_SF

CONTAINS

  !==================================================================
  ! PUBLIC SUBROUTINES
  !
  !    INIT_FVCOM
  !    COMMANDLINEIO
  !==================================================================

  !=========================================================================


  SUBROUTINE COMMANDLINEIO(CVS_ID,CVS_Date,CVS_Name,CVS_Revision)
    use mod_sng ! [mdl] String manipulatio
    USE control


    implicit none
    ! Parameters
    character(len=*), INTENT(IN)::CVS_Id  ! [sng] CVS Identification
    character(len=*), INTENT(IN)::CVS_Date ! [sng] Date string
    character(len=*), INTENT(IN)::CVS_Name ! [sng] File name string
    character(len=*), INTENT(IN)::CVS_Revision ! [sng] File revision string
    character(len=*),parameter::nlc=char(0) ! [sng] NUL character = ASCII 0 = char(0)

    ! Command-line parsing
    character(80)::arg_val ! [sng] command-line argument value
    character(200)::cmd_ln ! [sng] command-line
    character(80)::opt_sng ! [sng] Option string
    character(2)::dsh_key ! [sng] command-line dash and switch
    character(200)::prg_ID ! [sng] Program ID

    integer::arg_idx ! [idx] Counting index
    integer::arg_nbr ! [nbr] Number of command-line arguments
    integer::opt_lng ! [nbr] Length of option

    logical :: outtofile = .false.

    ! Initialize values from in variables 
    CaseName="empty"//nlc ! declared in mod_main - all_vars
    dbg_lvl=0 ! declared in mod_dbg
    INFOFILE="screen"//nlc ! Get output file parameter
    NameList_Name = "empty"//nlc ! Name of blank name list file-
    ! declared in mod_main

    ! Main code
    call ftn_strini(cmd_ln) ! [sng] sng(1:len)=NUL

    call ftn_cmd_ln_sng(cmd_ln) ! [sng] Re-construct command-line into single string
    call ftn_prg_ID_mk(CVS_Id,CVS_Revision,CVS_Date,prg_ID) ! [sng] Program ID

    arg_nbr=command_argument_count() ! [nbr] Number of command-line arguments

    if (arg_nbr .LE. 0 ) then
       if(MSR) WRITE(IPT,*) "You must specify a case name: "
       if(MSR) Call HelpTxt(IPT)
       call PSHUTDOWN
    end if

    arg_idx=1 ! [idx] Counting index
    do while (arg_idx <= arg_nbr)
       call ftn_getarg_wrp(arg_idx,arg_val) ! [sbr] Call getarg, increment arg_idx
       dsh_key=arg_val(1:2) ! [sng] First two characters of option
       if (dsh_key == "--") then
          opt_lng=ftn_opt_lng_get(arg_val) ! [nbr] Length of option
          if (opt_lng <= 0) then
             if(MSR) write(IPT,*) "Long option has no name"
             call PSHUTDOWN
          end if

          opt_sng=arg_val(3:2+opt_lng) ! [sng] Option string
          if (dbg_lvl >= dbg_io) then
             if(MSR) write (6,"(5a,i3)") prg_nm(1:ftn_strlen(prg_nm)), &
                  ": DEBUG Double hyphen indicates multi-character option: ", &
                  "opt_sng = ",opt_sng(1:ftn_strlen(opt_sng)),", opt_lng = ",opt_lng
          end if
          if (opt_sng == "dbg" .or. opt_sng == "dbg_lvl" ) then
             call ftn_arg_get(arg_idx,arg_val,dbg_lvl) ! [enm] Debugging level

          else if (opt_sng == "dbg_par" .or.opt_sng == "Dbg_Par"&
               & .or.opt_sng == "DBG_PAR") then

             dbg_par = .true.

          else if (opt_sng == "crashrestart" .or.opt_sng == "CrashRestart"&
               & .or.opt_sng == "CRASHRESTART") then

             CMDLN_RESTART = .true.
             !           call ftn_arg_get(arg_idx,arg_val,dbg_par) ! [sng] Input file

          else if (opt_sng == "CaseName" .or.opt_sng == "casename"&
               & .or.opt_sng == "CASENAME") then

             call ftn_arg_get(arg_idx,arg_val,CaseName) ! [sng] Input file
             CaseName=CaseName(1:ftn_strlen(CaseName))
             ! Convert back to a fortran string!

          else if (opt_sng == "Create_NameList" .or.opt_sng == "create_namelist"&
               & .or.opt_sng == "CREATE_NAMELIST") then

             call ftn_arg_get(arg_idx,arg_val,NAMELIST_NAME)
             NAMELIST_NAME = NAMELIST_NAME(1:ftn_strlen(NAMELIST_NAME))

             BLANK_NAMELIST = .true.

          else if (opt_sng == "LogFile" .or.opt_sng == "logfile"&
               & .or.opt_sng == "LOGFILE") then

             call ftn_arg_get(arg_idx,arg_val,INFOFILE) 
             INFOFILE=INFOFILE(1:ftn_strlen(INFOFILE))
             
             OUTTOFILE=.true.
             
          else if (opt_sng == "help" .or.opt_sng == "HELP" .or. opt_sng&
               & == "Help") then

             if(MSR) call HelpTxt(IPT)
             call PSHUTDOWN
!!$   THIS DOES NOT SEEM PRACTICAL - MODIFY THE RUN FILE INSTEAD
!!$          else if (opt_sng == "CrashRestart") then
!!$             call ftn_arg_get(arg_idx,arg_val,CrashRestart) ! [lgc] Logical

          else ! Option not recognized
             arg_idx=arg_idx-1 ! [idx] Counting index
             if(MSR) call ftn_getarg_err(arg_idx,arg_val) ! [sbr] Error handler for getarg()
          endif ! endif option is recognized
          ! Jump to top of while loop
          cycle ! C, F77, and F90 use "continue", "goto", and "cycle"
       endif ! endif long option

       if (dsh_key == "-V" .or.dsh_key == "-v" ) then

          if(MSR) write(IPT,*) prg_id
          call PSHUTDOWN

       else if (dsh_key == "-H" .or.dsh_key == "-h" ) then

          if(MSR) Call helptxt(IPT)
          Call PSHUTDOWN

       else ! Option not recognized
          arg_idx=arg_idx-1 ! [idx] Counting index
          if(MSR) call ftn_getarg_err(arg_idx,arg_val) ! [sbr] Error handler for getarg()
       endif ! endif arg_val


    end do ! end while (arg_idx <= arg_nbr)

    if (blank_namelist) then
       INFOFILE =trim(NameList_Name)//"_run.nml"
       CaseName= NAMELIST_NAME
       outtofile = .true.
    end if

    CALL dbg_init(IPT_BASE,outtofile)


  END SUBROUTINE COMMANDLINEIO




  SUBROUTINE NAME_LIST_INITIALIZE
    USE CONTROL

    IMPLICIT NONE

    !--Parameters in NameList NML_CASE
    CASE_TITLE = "'AN FVCOM CASE DESCRIPTION' - note string must be in 'quotes'"
    TIMEZONE = "Select Time Zone or for idealized case select 'none' (start time=0.0)" 
    USE_REAL_WORLD_TIME=.FALSE.
    DATE_FORMAT="A three letter string specify date format: 'YMD' or 'DMY'"
    DATE_REFERENCE= "Date (specified as a string -- example '2007-11-05 00:00:00') or 'default'"
    START_DATE= "Date and Time are specified as a string (example '2007-11-05 00:00:00')"
    END_DATE= "For an idealized case specify 'seconds=(flt)','days=(flt)', or 'cycles=(int)'"

    
    STARTUP_TYPE = "'hotstart', 'coldstart', 'forecast' or 'crashrestart'"
    STARTUP_FILE= trim(casename)//"_restart.nc"
    STARTUP_UV_TYPE = "'default' or 'set values'"
    STARTUP_TURB_TYPE = "'default' or 'set values'"
    STARTUP_TS_TYPE = "'constant' 'linear' 'observed' or 'set values'"
    STARTUP_T_VALS = -99.0_sp
    STARTUP_S_VALS = -99.0_sp
    STARTUP_U_VALS = -99.0_sp
    STARTUP_V_VALS = -99.0_sp
    STARTUP_DMAX = -99.0_sp
    
    !--Parameters in NameList NML_IO
    INPUT_DIR = "/Your/relative/path/to/input/files"
    OUTPUT_DIR = "/Your/relative/path/to/output/files :Must already exist!"
    IREPORT = 0
    VISIT_ALL_VARS = .False.
    WAIT_FOR_VISIT = .False.
    USE_MPI_IO_MODE = .FALSE.


    !--Parameters in NameList NML_INTEGRATION
    ExtStep_Seconds = 0.0
    ISplit = 0
    IRamp = 0 
    Static_SSH_Adj = 0.0
    Min_Depth = 0.0

    !--Parameters in NameList NML_RESTART
    RST_ON            = .False.
    RST_FIRST_OUT     = 'Date to start RESTART OUTPUT: Format the same as START_DATE'
    RST_OUT_INTERVAL  = "A length of time: 'seconds= ','days= ', or 'cycles= '"
    RST_OUTPUT_STACK  = 0

    !--Parameters in NameList NML_NETCDF
    NC_ON = .False.
    NC_FIRST_OUT = 'Date to start NETCDF OUTPUT: Format the same as START_DATE'
    NC_OUT_INTERVAL = "A length of time: 'seconds= ','days= ', or 'cycles= '"
    NC_OUTPUT_STACK = 0
    NC_SUBDOMAIN_FILES= "FVCOM"
    NC_GRID_METRICS = .False.
    NC_FILE_DATE    = .False.   
    NC_VELOCITY     = .False.
    NC_SALT_TEMP    = .False.
    NC_TURBULENCE   = .False.
    NC_AVERAGE_VEL  = .False.
    NC_VERTICAL_VEL = .False.
    NC_NH_QP        = .False.
    NC_NH_RHS       = .False.
    NC_WIND_VEL     = .False.
# if defined (AIR_PRESSURE) || (HEATING_CALCULATED)
    NC_ATM_PRESS    = .False.
# endif    
    NC_WIND_STRESS  = .False.
# if defined (WAVE_CURRENT_INTERACTION)
    NC_WAVE_PARA    = .False. 
    NC_WAVE_STRESS  = .False. 
# endif
    NC_EVAP_PRECIP  = .False.
    NC_SURFACE_HEAT = .False.
    NC_GROUNDWATER  = .False.
# if defined (ICE)
    NC_ICE          = .False.
# endif
    NC_VORTICITY    = .False.
    NC_WQM          = .False.
    NC_BIO          = .False.
    
    ! OUTPUT VARIABLES DEFAULT TO OFF

    !--Parameters in NameList NML_NETCDF_SURFACE
    NCSF_ON = .False.
    NCSF_FIRST_OUT = 'Date to start NETCDF OUTPUT: Format the same as START_DATE'
    NCSF_OUT_INTERVAL = "A length of time: 'seconds= ','days= ', or 'cycles= '"
    NCSF_OUTPUT_STACK = 0
    NCSF_SUBDOMAIN_FILES= "FVCOM"
    NCSF_GRID_METRICS = .False.
    NCSF_FILE_DATE    = .False.   
    NCSF_VELOCITY     = .False.
    NCSF_SALT_TEMP    = .False.
    NCSF_TURBULENCE   = .False.
    NCSF_WIND_VEL     = .False.
# if defined (AIR_PRESSURE) || (HEATING_CALCULATED)
    NCSF_ATM_PRESS    = .False.
# endif    
    NCSF_WIND_STRESS  = .False.
# if defined (WAVE_CURRENT_INTERACTION)
    NCSF_WAVE_PARA    = .False.
# endif
# if defined (ICE)
    NCSF_ICE          = .False.
# endif
    NCSF_EVAP_PRECIP  = .False.
    NCSF_SURFACE_HEAT = .False.
    
    ! OUTPUT VARIABLES DEFAULT TO OFF

    !--Parameters in NameList NML_NETCDF_AV
    NCAV_ON = .False.
    NCAV_FIRST_OUT = "Date to start NETCDF interval averaged output: Format the same as START_DATE"
    NCAV_OUT_INTERVAL = "A length of time: 'seconds= ','days= ', or 'cycles= '"
    NCAV_OUTPUT_STACK = 0
    NCAV_SUBDOMAIN_FILES= "FVCOM"
    NCAV_GRID_METRICS = .False.
    NCAV_FILE_DATE    = .False.   
    NCAV_VELOCITY     = .False.
    NCAV_SALT_TEMP    = .False.
    NCAV_TURBULENCE   = .False.
    NCAV_AVERAGE_VEL  = .False.
    NCAV_VERTICAL_VEL = .False.
    NCAV_NH_QP        = .False.
    NCAV_NH_RHS       = .False.
    NCAV_WIND_VEL     = .False.
# if defined (AIR_PRESSURE) || (HEATING_CALCULATED)
    NCAV_ATM_PRESS    = .False.  
# endif    
    NCAV_WIND_STRESS  = .False.
    NCAV_WAVE_PARA    = .False.  !Jadon
    NCAV_WAVE_STRESS  = .False.  !Jadon
    NCAV_EVAP_PRECIP  = .False.
    NCAV_SURFACE_HEAT = .False.
    NCAV_GROUNDWATER  = .False.
    NCAV_VORTICITY    = .False.
    NCAV_WQM          = .False.
    NCAV_BIO          = .False.


    ! OUTPUT VARIABLES DEFAULT TO OFF

    !--Parameters in NameList NML_PHYSICS
    HORIZONTAL_MIXING_TYPE = "'closure' or 'constant'"
    HORIZONTAL_MIXING_FILE = trim(casename)//"_hvc.nc"
    HORIZONTAL_MIXING_KIND = "Options:"//TRIM(CNSTNT)//","//TRIM(STTC)
    HORIZONTAL_MIXING_COEFFICIENT =-1.0_SP
    HORIZONTAL_PRANDTL_NUMBER = -1.0_SP

    VERTICAL_MIXING_TYPE = "'closure' or 'constant'"
    VERTICAL_MIXING_COEFFICIENT = -1.0_SP
    VERTICAL_PRANDTL_NUMBER = -1.0_SP


    BOTTOM_ROUGHNESS_MINIMUM = -1.0_SP
    BOTTOM_ROUGHNESS_LENGTHSCALE = -1.0_SP
    BOTTOM_ROUGHNESS_TYPE   = "'"//TRIM(BR_ORIG)//"', or '"&
         &//TRIM(BR_GOTM)//"'; Select your bottom roughness equation (brough.F)"
    BOTTOM_ROUGHNESS_KIND = "Options:"//TRIM(CNSTNT)//","//TRIM(STTC)
    BOTTOM_ROUGHNESS_FILE =trim(casename)//"_brf.nc"
    
    CONVECTIVE_OVERTURNING     = .FALSE.  
    SCALAR_POSITIVITY_CONTROL  = .FALSE.
    BAROTROPIC                 = .FALSE. 

    BAROCLINIC_PRESSURE_GRADIENT = "'sigma levels' or 'z coordinates'; select method of calculation"

    SEA_WATER_DENSITY_FUNCTION = "'"//TRIM(SW_DENS1)//"', '"&
         &//TRIM(SW_DENS2)//"', or '"//TRIM(SW_DENS3)//"; Select your equation of state (eqs_of_state.F)"

    RECALCULATE_RHO_MEAN = .FALSE.
    INTERVAL_RHO_MEAN = "A length of time or number of cycles in standard format"

    TEMPERATURE_ACTIVE  = .FALSE.          
    SALINITY_ACTIVE     = .FALSE.          
    SURFACE_WAVE_MIXING = .FALSE.          
    WETTING_DRYING_ON   = .FALSE.
!J. Ge
    ! for tracer advection
    BACKWARD_ADVECTION = .FALSE.
    BACKWARD_STEP = -1
!J. Ge

    ADCOR_ON = .TRUE.
    EQUATOR_BETA_PLANE = .FALSE.
    NOFLUX_BOT_CONDITION = .TRUE.

    !--Parameters in NameList NML_SURFACE_FORCING
    WIND_ON = .FALSE.
    WIND_TYPE = "Options::"//TRIM(SPEED)//","//TRIM(STRESS)
    WIND_FILE = trim(casename)//"_wnd.nc"
    WIND_KIND = "Options:"//TRIM(CNSTNT)//","//TRIM(STTC)//","//TRIM(TMDPNDNT)//","//TRIM(PRDC)//","//TRIM(VRBL)
    WIND_X = 0.0_SP
    WIND_Y = 0.0_SP
    HEATING_ON = .FALSE.
    HEATING_TYPE = "'body' or 'flux'"
    HEATING_FILE = trim(casename)//"_hfx.nc"
    HEATING_KIND = "Options:"//TRIM(CNSTNT)//","//TRIM(STTC)//","//TRIM(TMDPNDNT)//","//TRIM(PRDC)//","//TRIM(VRBL)
    HEATING_RADIATION = 0.0_SP
    HEATING_NETFLUX = 0.0_SP
    HEATING_LONGWAVE_PERCTAGE = 0.78_SP
    HEATING_LONGWAVE_LENGTHSCALE = 1.4_SP
    HEATING_SHORTWAVE_LENGTHSCALE= 6.3_SP
    PRECIPITATION_ON = .FALSE.
    PRECIPITATION_FILE = trim(casename)//"_emp.nc"
    PRECIPITATION_KIND = "Options:"//TRIM(CNSTNT)//","//TRIM(STTC)//","//TRIM(TMDPNDNT)//","//TRIM(PRDC)//","//TRIM(VRBL)
    PRECIPITATION_PRC = 0.0_SP
    PRECIPITATION_EVP = 0.0_SP
    AIRPRESSURE_ON    = .FALSE.
    AIRPRESSURE_KIND  = "Options:"//TRIM(CNSTNT)//","//TRIM(STTC)//","//TRIM(TMDPNDNT)//","//TRIM(PRDC)//","//TRIM(VRBL)
    AIRPRESSURE_FILE  = trim(casename)//"_aip.nc"
    AIRPRESSURE_VALUE = 0.0_SP
!    Jadon
    WAVE_ON = .FALSE.
    WAVE_FILE = trim(casename)//"_wav.nc"
    WAVE_KIND = "Options:"//TRIM(CNSTNT)//","//TRIM(STTC)//","//TRIM(TMDPNDNT)//","//TRIM(PRDC)//","//TRIM(VRBL)
    WAVE_HEIGHT     = 0.0_SP
    WAVE_LENGTH     = 0.0_SP
    WAVE_DIRECTION  = 0.0_SP 
    WAVE_PERIOD     = 0.0_SP
    WAVE_PER_BOT    = 0.0_SP
    WAVE_UB_BOT     = 0.0_SP

    !--Parameters in NameList NML_RIVER_TYPE
    RIVER_NUMBER = -1
    RIVER_KIND = "Options:"//TRIM(PRDC)//" or "//TRIM(VRBL)
    RIVER_TS_SETTING = "'calculated' or 'specified'"
    RIVER_INFLOW_LOCATION = "'node' or 'edge'"
    RIVER_INFO_FILE = "'default' or 'filename'"


    !--Parameters in NameList NML_RIVERS
    RIVER_NAME = "River Name in netcdf data file; use mulitple namelists for multiple rivers!"
    RIVER_FILE = trim(casename)//"_riv.nc"
    RIVER_GRID_LOCATION = -1
# if defined(RIVER_FLOAT)
    RIVER_VERTICAL_DISTRIBUTION = -99
# else
    RIVER_VERTICAL_DISTRIBUTION = "FUNCTIONAL VERTICAL RIVER DISTROBUTION: SEE FVCOM MANUAL FOR HELP"
# endif

    !--Parameters in NameList NML_OPEN_BOUNDARY
    OBC_ON = .FALSE.
    OBC_NODE_LIST_FILE = trim(casename)//"_obc.dat"
    OBC_ELEVATION_FORCING_ON = .False.
    OBC_ELEVATION_FILE = trim(casename)//"_obc.nc "
    OBC_TS_TYPE = -1
    OBC_TEMP_NUDGING = .False.
    OBC_TEMP_FILE = trim(casename)//"_obc.nc "
    OBC_TEMP_NUDGING_TIMESCALE = 0.0
    OBC_SALT_NUDGING = .False.
    OBC_SALT_FILE = trim(casename)//"_obc.nc "
    OBC_SALT_NUDGING_TIMESCALE = 0.0
#   if defined (WATER_QUALITY)
    OBC_WQM_NUDGING = .False.
    OBC_WQM_FILE = trim(casename)//"_obc.nc "
    OBC_WQM_NUDGING_TIMESCALE = 0.0
#   endif    
#   if defined (BioGen)
    OBC_BIO_NUDGING = .False.
    OBC_BIO_FILE = trim(casename)//"_obc.nc "
    OBC_BIO_NUDGING_TIMESCALE = 0.0
#   endif    
    OBC_MEANFLOW = .FALSE.
    OBC_MEANFLOW_FILE = trim(casename)//"_obc.nc"
    OBC_LONGSHORE_FLOW_ON = .FALSE.
    OBC_LONGSHORE_FLOW_FILE = trim(casename)//"_lsf.dat"
    OBC_TIDEOUT_INITIAL  = 0       !TIME STEPS
    OBC_TIDEOUT_INTERVAL = 0       !TIME STEPS
    OBC_DEPTH_CONTROL_ON = .TRUE.

    !--Parameters in NameList GRID_COORDINATES
    GRID_FILE = trim(casename)//"_grd.dat"
    GRID_FILE_UNITS = "Can be 'degrees' or 'meters'; certain make options required"
    PROJECTION_REFERENCE = "none: A recognized reference coordinate for proj&
         &tion for PROJ4"
    SIGMA_LEVELS_FILE = trim(casename)//"_sigma.dat"
    CORIOLIS_FILE = trim(casename)//"_cor.dat"
    DEPTH_FILE = trim(casename)//"_dep.dat"
    SPONGE_FILE = trim(casename)//"_spg.dat"


    !--Parameters in NameList NML_GROUNDWATER
    GROUNDWATER_ON = .False.
    GROUNDWATER_SALT_ON = .False.
    GROUNDWATER_TEMP_ON = .False.
    GROUNDWATER_KIND ="Options:"//TRIM(CNSTNT)//","//TRIM(STTC)&
         &//","//TRIM(TMDPNDNT)//","//TRIM(PRDC)//","//TRIM(VRBL)
    GROUNDWATER_FILE = trim(casename)//"_grndwtr.nc"
    GROUNDWATER_FLOW = 0.0
    GROUNDWATER_TEMP = 0.0
    GROUNDWATER_SALT = 0.0

    !--Parameters in NameList NML_LAG_PART
    LAG_PARTICLES_ON = .False.
    LAG_START_FILE   = "init_lag.nc"
    LAG_OUT_FILE     = "lag_out.nc"
    LAG_FIRST_OUT    = "A Date or time"
    LAG_RESTART_FILE = "lag_restart.nc"
    LAG_OUT_INTERVAL = "A length of time: 'seconds= ','days= ', or 'cycles= '"
    LAG_SCAL_CHOICE  = "none"


    !--Parameters in NameList NML_ADDITIONAL_MODELS
!    WATER_QUALITY_MODEL = .FALSE.
!    WATER_QUALITY_MODEL_FILE = "DO NOT ADD UNTILL FVCOM IS RUNNING BY ITS SELF FIRST"
    DATA_ASSIMILATION  = .FALSE.
    DATA_ASSIMILATION_FILE = "./"//trim(casename)//"_run.nml"
    BIOLOGICAL_MODEL= .FALSE.
!--------- J. Ge for biology --------------
    BIOLOGICAL_MODEL_FILE = "DO NOT ADD UNTILL FVCOM IS RUNNING BY ITS SELF FIRST"
!------------------------------------------
    STARTUP_BIO_TYPE = "'observed' use this option only now"   ! constant, linear, observed, set values
    SEDIMENT_MODEL= .FALSE.
    SEDIMENT_MODEL_FILE = "DO NOT ADD UNTILL FVCOM IS RUNNING BY ITS SELF FIRST"
# if defined (OFFLINE_SEDIMENT)
    OFFLINE_SEDIMENT_FILE = "DO NOT ADD UNTILL FVCOM IS RUNNING BY ITS SELF FIRST"
# endif
# if defined (OFFLINE_BIOLOGY)
    OFFLINE_BIOLOGY_FILE = "DO NOT ADD UNTILL FVCOM IS RUNNING BY ITS SELF FIRST"
# endif
    SEDIMENT_PARAMETER_TYPE = "DO NOT ADD UNTILL FVCOM IS RUNNING BY ITS SELF FIRST"
    SEDIMENT_PARAMETER_FILE = "DO NOT ADD UNTILL FVCOM IS RUNNING BY ITS SELF FIRST"
    BEDFLAG_TYPE = "DO NOT ADD UNTIL FVCOM IS RUNNING BY ITS SELF FIRST"
    BEDFLAG_FILE = "DO NOT ADD UNTILL FVCOM IS RUNNING BY ITS SELF FIRST"
    ICE_MODEL= .FALSE.
    ICE_FORCING_FILE = "DO NOT ADD UNTILL FVCOM IS RUNNING BY ITS SELF FIRST"
    ICE_FORCING_KIND ="Options:"//TRIM(CNSTNT)//","//TRIM(STTC)&
         &//","//TRIM(TMDPNDNT)//","//TRIM(PRDC)//","//TRIM(VRBL)
    ICE_SEA_LEVEL_PRESSURE = 0.0
    ICE_AIR_TEMP           = 0.0
    ICE_SPEC_HUMIDITY      = 0.0
    ICE_CLOUD_COVER        = 0.0
    ICE_SHORTWAVE          = 0.0
    ICE_LONGWAVE_TYPE      = "'PW' or 'RM'" !PW -- longwave as in Parkinson and Washington (1979)
                                            !RM -- longwave, Rosati and Miyakoda, JPO 18, p. 1607 (1988)

    ICING_MODEL= .FALSE.
    ICING_FORCING_FILE = "DO NOT ADD UNTILL FVCOM IS RUNNING BY ITS SELF FIRST"
    ICING_FORCING_KIND ="Options:"//TRIM(CNSTNT)//","//TRIM(STTC)&
         &//","//TRIM(TMDPNDNT)//","//TRIM(PRDC)//","//TRIM(VRBL)
    ICING_AIR_TEMP = 0.0
    ICING_WSPD    = 0.0

    PROBES_ON = .FALSE.
    PROBES_NUMBER =0
    PROBES_FILE = "Probe namelist file name"

    HIGH_LATITUDE_WAVE = .FALSE.

  END SUBROUTINE NAME_LIST_INITIALIZE


  SUBROUTINE NAME_LIST_PRINT
    USE CONTROL

    IMPLICIT NONE


    ! MODIFY THIS ROUTINE TO PRINT EACH NAME LIST TO A CHARCTER STRING:
    !  PARSE AND FORMAT THE STRING TO MAKE IT LOOK 'PRETTY'!

    write(UNIT=IPT,NML=NML_CASE)
    write(UNIT=IPT,NML=NML_STARTUP)
    write(UNIT=IPT,NML=NML_IO)
    write(UNIT=IPT,NML=NML_INTEGRATION)
    write(UNIT=IPT,NML=NML_RESTART)
    write(UNIT=IPT,NML=NML_NETCDF)
    write(UNIT=IPT,NML=NML_NETCDF_SURFACE)
    write(UNIT=IPT,NML=NML_NETCDF_AV)
    write(UNIT=IPT,NML=NML_SURFACE_FORCING)
    write(UNIT=IPT,NML=NML_PHYSICS)
    write(UNIT=IPT,NML=NML_RIVER_TYPE)
    write(UNIT=IPT,NML=NML_RIVER)
    write(UNIT=IPT,NML=NML_OPEN_BOUNDARY_CONTROL)
    write(UNIT=IPT,NML=NML_GRID_COORDINATES)
    write(UNIT=IPT,NML=NML_GROUNDWATER)
    write(UNIT=IPT,NML=NML_LAG)
    write(UNIT=IPT,NML=NML_ADDITIONAL_MODELS)
    write(UNIT=IPT,NML=NML_PROBES)
    write(UNIT=IPT,NML=NML_BOUNDSCHK)  !bounds checking

  RETURN
  END SUBROUTINE NAME_LIST_PRINT


  SUBROUTINE NAME_LIST_READ
     USE CONTROL
    IMPLICIT NONE
    integer :: ios, i
    Character(Len=120):: FNAME
    character(len=160) :: pathnfile
    if(DBG_SET(dbg_sbr)) &
         & write(IPT,*) "Subroutine Begins: Read_Name_List;"

    ios = 0

    FNAME = "./"//trim(casename)//"_run.nml"

    if(DBG_SET(dbg_io)) &
         & write(IPT,*) "Read_Name_List: File: ",trim(FNAME)

    CALL FOPEN(NMLUNIT,trim(FNAME),'cfr')

    !READ NAME LIST FILE

    ! Read IO Information
    READ(UNIT=NMLUNIT, NML=NML_IO,IOSTAT=ios)
    if(ios .NE. 0 ) Then
       if(DBG_SET(dbg_log)) write(UNIT=IPT,NML=NML_IO)
       Call Fatal_error("Can Not Read NameList NML_IO from file: "//trim(FNAME))
    end if
    CALL CHECK_IO_DIRS

    REWIND(NMLUNIT)

    if(DBG_SET(dbg_scl)) &
         & write(IPT,*) "Read_Name_List:"

    if(DBG_SET(dbg_scl)) &
         & write(UNIT=IPT,NML=NML_IO)

    ! Read Case Information
    READ(UNIT=NMLUNIT, NML=NML_CASE,IOSTAT=ios)
    if(ios .NE. 0 ) THEN
       if(DBG_SET(dbg_log)) write(UNIT=IPT,NML=NML_CASE)
       CALL FATAL_ERROR("Can Not Read NameList NML_CASE from file: "//trim(FNAME))
    end if

    REWIND(NMLUNIT)

    if(DBG_SET(dbg_scl)) &
         & write(IPT,*) "Read_Name_List:"

    if(DBG_SET(dbg_scl)) &
         & write(UNIT=IPT,NML=NML_CASE)

    ! Read STARTUP TYPE INFORMATION
    READ(UNIT=NMLUNIT, NML=NML_STARTUP,IOSTAT=ios)
    if(ios .NE. 0 ) Then
       if(DBG_SET(dbg_log)) write(UNIT=IPT,NML=NML_STARTUP)
       Call Fatal_Error("Can Not Read NameList NML_STARTUP from file: "//trim(FNAME))
    end if

    REWIND(NMLUNIT)

    if(DBG_SET(dbg_scl)) &
         & write(IPT,*) "Read_Name_List:"

    if(DBG_SET(dbg_scl)) &
         & write(UNIT=IPT,NML=NML_STARTUP)

    ! Read Integration Settings
    READ(UNIT=NMLUNIT, NML=NML_INTEGRATION,IOSTAT=ios)
    if(ios .NE. 0 ) Then
       if(DBG_SET(dbg_log)) write(UNIT=IPT,NML=NML_INTEGRATION)
       Call Fatal_Error("Can Not Read NameList NML_INTEGRATION from file: "//trim(FNAME))
    end if

    REWIND(NMLUNIT)

    if(DBG_SET(dbg_scl)) &
         & write(IPT,*) "Read_Name_List:"

    if(DBG_SET(dbg_scl)) &
         & write(UNIT=IPT,NML=NML_INTEGRATION)

    ! Read Netcdf Output Settings
    READ(UNIT=NMLUNIT, NML=NML_RESTART,IOSTAT=ios)
    if(ios .NE. 0 ) then
       if(DBG_SET(dbg_log)) write(UNIT=IPT,NML=NML_RESTART)
       Call Fatal_Error("Can Not Read NameList NML_RESTART from file: "//trim(FNAME))
    end if

    REWIND(NMLUNIT)

    if(DBG_SET(dbg_scl)) &
         & write(IPT,*) "Read_Name_List:"

    if(DBG_SET(dbg_scl)) &
         & write(UNIT=IPT,NML=NML_RESTART)

    ! Read Netcdf Output Settings
    READ(UNIT=NMLUNIT, NML=NML_NETCDF,IOSTAT=ios)
    if(ios .NE. 0 ) then
       if(DBG_SET(dbg_log)) write(UNIT=IPT,NML=NML_NETCDF)
       Call Fatal_Error("Can Not Read NameList NML_NETCDF from file: "//trim(FNAME))
    end if

    REWIND(NMLUNIT)

    if(DBG_SET(dbg_scl)) &
         & write(IPT,*) "Read_Name_List:"

    if(DBG_SET(dbg_scl)) &
         & write(UNIT=IPT,NML=NML_NETCDF)

    ! Read Netcdf Output Settings
    READ(UNIT=NMLUNIT, NML=NML_NETCDF_SURFACE,IOSTAT=ios)
    if(ios .NE. 0 ) then
       if(DBG_SET(dbg_log)) write(UNIT=IPT,NML=NML_NETCDF_SURFACE)
       Call Fatal_Error("Can Not Read NameList NML_NETCDF_SURFACE from file: "//trim(FNAME))
    end if

    REWIND(NMLUNIT)

    if(DBG_SET(dbg_scl)) &
         & write(IPT,*) "Read_Name_List:"

    if(DBG_SET(dbg_scl)) &
         & write(UNIT=IPT,NML=NML_NETCDF_SURFACE)

    ! Read Netcdf Average Output Settings
    READ(UNIT=NMLUNIT, NML=NML_NETCDF_AV,IOSTAT=ios)
    if(ios .NE. 0 ) then
       if(DBG_SET(dbg_log)) write(UNIT=IPT,NML=NML_NETCDF_AV)
       Call Fatal_Error("Can Not Read NameList NML_NETCDF_AV from file: "//trim(FNAME))
    end if

    REWIND(NMLUNIT)

    if(DBG_SET(dbg_scl)) &
         & write(IPT,*) "Read_Name_List:"

    if(DBG_SET(dbg_scl)) &
         & write(UNIT=IPT,NML=NML_NETCDF_AV)

    ! Read Model Physics Settings
    READ(UNIT=NMLUNIT, NML=NML_PHYSICS,IOSTAT=ios)
    if(ios .NE. 0 ) then
       if(DBG_SET(dbg_log)) write(UNIT=IPT,NML=NML_PHYSICS)
       Call Fatal_Error("Can Not Read NameList NML_PHYSICS from file: "//trim(FNAME))
    end if

    REWIND(NMLUNIT)

    if(DBG_SET(dbg_scl)) &
         & write(IPT,*) "Read_Name_List:"

    if(DBG_SET(dbg_scl)) &
         & write(UNIT=IPT,NML=NML_PHYSICS)


    ! Read Surface Forcing Settings
    READ(UNIT=NMLUNIT, NML=NML_SURFACE_FORCING,IOSTAT=ios)
    if(ios .NE. 0 ) then
       if(DBG_SET(dbg_log)) write(UNIT=IPT,NML=NML_SURFACE_FORCING)
       Call Fatal_Error("Can Not Read NameList NML_SURFACE_FORCING from file: "//trim(FNAME))
    end if

    REWIND(NMLUNIT)

    if(DBG_SET(dbg_scl)) &
         & write(IPT,*) "Read_Name_List:"

    if(DBG_SET(dbg_scl)) &
         & write(UNIT=IPT,NML=NML_SURFACE_FORCING)


    ! Read River settings
    READ(UNIT=NMLUNIT, NML=NML_RIVER_TYPE,IOSTAT=ios)
    if(ios .NE. 0 ) then
       if(DBG_SET(dbg_log)) write(UNIT=IPT,NML=NML_RIVER_TYPE)
       Call Fatal_Error("Can Not Read NameList NML_RIVER_TYPE from file: "//trim(FNAME))
    end if

    if(DBG_SET(dbg_scl)) &
         & write(IPT,*) "Read_Name_List:"

    if(DBG_SET(dbg_scl)) &
         & write(UNIT=IPT,NML=NML_RIVER_TYPE)

    IF(RIVER_NUMBER > 0) THEN
       SELECT CASE (RIVER_INFO_FILE(2:8))
       CASE('default')
          REWIND(NMLUNIT)
          RIVERNMLUNIT=NMLUNIT
         !!  'default' use present runfile namelist
       CASE DEFAULT
         !!  use the specified river file information
          RIVERNMLUNIT=800
          pathnfile = trim(INPUT_DIR)//trim(RIVER_INFO_FILE)
          Call FOPEN(RIVERNMLUNIT,trim(pathnfile),'cfr')
       END SELECT

       ALLOCATE(RIVERS(RIVER_NUMBER))

       ! Read River Namelists...
       ios = 0
       i = 0
       DO 
          READ(UNIT=RIVERNMLUNIT, NML=NML_RIVER,IOSTAT=ios)

          if (IOS /= 0 ) exit
          I = I +1
          if  (I > RIVER_NUMBER) exit ! To prevent sigsev...
    
          RIVERS(i)%NAME=RIVER_NAME
          RIVERS(i)%FILE=RIVER_FILE
          RIVERS(i)%LOCATION=RIVER_GRID_LOCATION
          RIVERS(i)%DISTRIBUTION=RIVER_VERTICAL_DISTRIBUTION
       END DO

       IF(I .NE. RIVER_NUMBER) THEN
          if(DBG_SET(dbg_log))  then
             write(ipt,*)"Bad River data in the Name List!"
             write(ipt,*)"Specified number of rivers:",RIVER_NUMBER
             write(ipt,*)"But Found",I, "; Valid river name list objects.(Printing Last)"
             write(UNIT=IPT,NML=NML_RIVER)
          end if
          
          CALL FATAL_ERROR('PLEASE REPAIR THE NAME LIST SO IT IS CONSISTANT... see above')
       END IF
       
    ELSE IF (RIVER_NUMBER .eq. 0) THEN
       
       READ(UNIT=NMLUNIT, NML=NML_RIVER,IOSTAT=ios)
       !THERE SHOULD BE NO RIVER NAME LISTS
       if (IOS == 0 ) CALL FATAL_ERROR &
            & ('THERE ARE ONE OR MORE RIVER NAME LISTS, BUT RIVER TYPE SPECIFIED ZERO RIVERS?')
       
    ELSE
       if(DBG_SET(dbg_log)) write(UNIT=IPT,NML=NML_RIVER_TYPE)
       CALL FATAL_ERROR("YOU CAN'T HAVE A NEGATIVE NUMBER OF RIVERS!")
    END IF


    REWIND(NMLUNIT)


    ! Read Open Boundary Control Settings
    READ(UNIT=NMLUNIT, NML=NML_OPEN_BOUNDARY_CONTROL,IOSTAT=ios)
    if(ios .NE. 0 ) then
       if(DBG_SET(dbg_log)) write(UNIT=IPT,NML=NML_OPEN_BOUNDARY_CONTROL)
       Call Fatal_Error("Can Not Read NameList NML_OPEN_BOUNDARY_CONTROL from file: "//trim(FNAME))
    end if

    REWIND(NMLUNIT)

    if(DBG_SET(dbg_scl)) &
         & write(IPT,*) "Read_Name_List:"

    if(DBG_SET(dbg_scl)) &
         & write(UNIT=IPT,NML=NML_OPEN_BOUNDARY_CONTROL)


    ! Read Grid Coordinates Settings
    READ(UNIT=NMLUNIT, NML=NML_GRID_COORDINATES,IOSTAT=ios)
    if(ios .NE. 0 ) then
       if(DBG_SET(dbg_log)) write(UNIT=IPT,NML=NML_GRID_COORDINATES)
       Call Fatal_Error("Can Not Read NameList NML_GRID_COORDINATES from file: "//trim(FNAME))
    end if

    REWIND(NMLUNIT)

    if(DBG_SET(dbg_scl)) &
         & write(IPT,*) "Read_Name_List:"

    if(DBG_SET(dbg_scl)) &
         & write(UNIT=IPT,NML=NML_GRID_COORDINATES)


    ! Read Groundwater Settings
    READ(UNIT=NMLUNIT, NML=NML_GROUNDWATER,IOSTAT=ios)
    if(ios .NE. 0 ) then
       if(DBG_SET(dbg_log)) write(UNIT=IPT,NML=NML_GROUNDWATER)
       Call Fatal_Error("Can Not Read NameList NML_GROUNDWATER from file: "//trim(FNAME))
    end if

    REWIND(NMLUNIT)

    if(DBG_SET(dbg_scl)) &
         & write(IPT,*) "Read_Name_List:"

    if(DBG_SET(dbg_scl)) &
         & write(UNIT=IPT,NML=NML_GROUNDWATER)

    ! Read LAG
    READ(UNIT=NMLUNIT, NML=NML_LAG,IOSTAT=ios)
    if(ios .NE. 0 ) then
       if(DBG_SET(dbg_log)) write(UNIT=IPT,NML=NML_LAG)
       Call Fatal_Error("Can Not Read NameList NML_LAG from file: "//trim(FNAME))
    end if
       
    REWIND(NMLUNIT)

    if(DBG_SET(dbg_scl)) &
         & write(IPT,*) "Read_Name_List:"

    if(DBG_SET(dbg_scl)) &
         & write(UNIT=IPT,NML=NML_LAG)


    ! Read Additional Models Settings
    READ(UNIT=NMLUNIT, NML=NML_ADDITIONAL_MODELS,IOSTAT=ios)
    if(ios .NE. 0 ) then
       if(DBG_SET(dbg_log)) write(UNIT=IPT,NML=NML_ADDITIONAL_MODELS)
       Call Fatal_Error("Can Not Read NameList NML_ADDITIONAL_MODELS from file: "//trim(FNAME))
    end if

    REWIND(NMLUNIT)

    if(DBG_SET(dbg_scl)) &
         & write(IPT,*) "Read_Name_List:"

    if(DBG_SET(dbg_scl)) &
         & write(UNIT=IPT,NML=NML_ADDITIONAL_MODELS)

    ! Read PROBE Settings
    READ(UNIT=NMLUNIT, NML=NML_PROBES,IOSTAT=ios)
    if(ios .NE. 0 ) then
       if(DBG_SET(dbg_log)) write(UNIT=IPT,NML=NML_PROBES)
       Call Fatal_Error("Can Not Read NameList NML_PROBES from file: "//trim(FNAME))
    end if

    ! Read BOUNDS CHECK (THRESHOLD SHUTDOWN)  Settings
    !=> bounds checking
    READ(UNIT=NMLUNIT, NML=NML_BOUNDSCHK,IOSTAT=ios)
    if(ios .NE. 0 ) then
       if(DBG_SET(dbg_log)) write(UNIT=IPT,NML=NML_BOUNDSCHK)
       Call Fatal_Error("Can Not Read NameList NML_BOUNDSCHK from file: "//trim(FNAME))
    end if
    !<= bounds checking

    REWIND(NMLUNIT)

    if(DBG_SET(dbg_scl)) &
         & write(IPT,*) "Read_Name_List:"

    if(DBG_SET(dbg_scl)) &
         & write(UNIT=IPT,NML=NML_PROBES)
    

    CLOSE(NMLUNIT)

! ----- END
    if(DBG_SET(dbg_sbr)) &
         & write(IPT,*) "Subroutine Ends: Read_Name_List;"    
  END SUBROUTINE NAME_LIST_READ
  

  SUBROUTINE CHECK_IO_DIRS
    USE CONTROL
    IMPLICIT NONE
    integer :: ncfileind, datfileind,ios,charnum, i,ios2
    logical :: fexist,back,connected
    character(len=100) :: testchar
    character(len=160) :: pathnfile
    character(len=2) :: cios
    ! CHECK FOR INPUT AND OUTPUT DIRECTORIES


    ! Check for '/' at the end of directory strings

    back = .true.
    charnum=index(trim(input_dir),"/",back)
    if (charnum /= len_trim(input_dir)) then
       input_dir=trim(input_dir)//"/"
       if(DBG_SET(dbg_io)) &
            write(ipt,*) "Added '/' to input_dir: ",trim(input_dir)
    end if

    charnum= index(trim(output_dir),"/",back)
    if ( charnum /= len_trim(output_dir)) then
       output_dir=trim(output_dir)//"/"
       if(DBG_SET(dbg_io)) &
            write(ipt,*) "Added '/' to output_dir: ",trim(output_dir)
    end if

    ! CHECK for Existance of output_dir and input_dir
    if (MSR) then
       ! OUTPUT_DIR TEST FILE - must exist + Read/write permissions
       testchar = trim(output_dir)//".fvcomtestfile"

       OPEN(UNIT=TESTUNIT,FILE=trim(testchar),&
            & FORM="formatted",STATUS="unknown",IOSTAT=ios)
       write(cios,'(i2.2)') ios
       if (ios == 9) then
          CALL FATAL_ERROR("Unable to OPEN the test file:",&
               & trim(testchar), "IOSTAT ERROR#"//cios//"; suggests ba&
               &d permissions in:", trim(output_dir))

       elseif (ios ==29) then
          CALL FATAL_ERROR("Unable to OPEN the test file:",&
               & trim(testchar), "IOSTAT ERROR#"//cios//"; suggests ba&
               &d directory path:", trim(output_dir))
   
       else if (ios /= 0) then
          CALL FATAL_ERROR("Unable to OPEN the test file:",&
               & trim(testchar), "Unknown IOSTAT error# "//cios)
       end if
          
       write(TESTUNIT,*)"This is a test file created by FVCOM. You may delete it."
       write(TESTUNIT,*)"Have a nice day."
       
       CLOSE(UNIT=TESTUNIT,IOSTAT=ios2)
       if (ios2 /= 0) then
          write(cios,'(i2.2)') ios2
          CALL FATAL_ERROR("Unable to CLOSE the test file:",&
               & trim(testchar), "Unknown IOSTAT error# "//cios)
       end if


       ! INPUT_DIR TEST FILE
       testchar = trim(input_dir)//".fvcomtestfile"
       OPEN(UNIT=TESTUNIT,FILE=trim(testchar),&
            & FORM="formatted",STATUS="unknown",IOSTAT=ios)
       write(cios,'(i2.2)') ios
       if (ios == 9) then
          CALL WARNING("Unable to OPEN the test file:",&
               & trim(testchar), "IOSTAT ERROR#"//cios//"; suggests bad permissions&
               & in:", trim(input_dir))
       elseif (ios ==29) then
          CALL WARNING("Unable to OPEN the test file:",&
               & trim(testchar), "IOSTAT ERROR#"//cios//"; suggests ba&
               &d directory path:",trim(input_dir))
       else if (ios /= 0) then
          CALL FATAL_ERROR("Unable to OPEN the test file:",&
               & trim(testchar), "Unknown IOSTAT error# "//cios)

       else if (ios == 0) then
          write(TESTUNIT,*)"This is a test file created by FVCOM. You may delete it."
          write(TESTUNIT,*)"Have a nice day."
          CLOSE(UNIT=TESTUNIT,IOSTAT=ios2)
          if (ios2 /= 0) then
             write(cios,'(i2.2)') ios2
             CALL FATAL_ERROR("Unable to CLOSE the test file:",&
                  & trim(testchar), "Unknown IOSTAT error# "//cios)
          end if
       end if
          

    end if

  END SUBROUTINE CHECK_IO_DIRS  
!=============================================================  
  SUBROUTINE NULLIFY_FILE_POINTERS
    IMPLICIT NONE

    NULLIFY(NC_DAT)
    NULLIFY(NC_AVG)
    NULLIFY(NC_RST)
    NULLIFY(NC_START)
    NULLIFY(NC_SF)
  END SUBROUTINE NULLIFY_FILE_POINTERS
!=============================================================  
  SUBROUTINE INCRIMENT_FNAME(FNAME)
    IMPLICIT NONE
    CHARACTER(LEN=160) FNAME
    CHARACTER(LEN=4) NUMSTR
    integer LENGTH,start,end,value
    if(DBG_SET(dbg_sbr)) &
         & write(IPT,*) "STARTING INCRIMENT_FNAME"

    if(DBG_SET(dbg_SBRIO)) write(IPT,*) 'INCRIMENTING OLD FILE NAME: '//TRIM(FNAME)

    ! GET POSITION OF FILE NUMBER IN FNAME
    LENGTH = LEN_TRIM(FNAME)
    start = LENGTH - 6
    end = LENGTH - 3

    ! READ FILE NUMBER AND INCRIMENT BY ONE
    NUMSTR=FNAME(start:end)
    read(NUMSTR,*) value
    value = value + 1
    write(NUMSTR,'(I4.4)') value

    ! INSERT BACK INTO FNAME
    FNAME(START:END) = NUMSTR

    if(DBG_SET(dbg_SBRIO)) write(IPT,*) 'NEW FILE NAME: '//TRIM(FNAME)

   if(DBG_SET(dbg_sbr)) &
         & write(IPT,*) "END INCRIMENT_FNAME"

  END SUBROUTINE INCRIMENT_FNAME
!=============================================================  
  SUBROUTINE OPEN_STARTUP_FILE
    USE CONTROL
    IMPLICIT NONE
    TYPE(NCFILE), POINTER :: NCF
    integer :: ncfileind, datfileind,ios,charnum, i
    logical :: fexist,back,connected
    character(len=100) :: testchar
    character(len=160) :: pathnfile
    character(len=2) :: cios
    ! CHECK FOR INPUT AND OUTPUT DIRECTORIES

    back = .true.
    
    ! TEST FILE NAME
    charnum = index (STARTUP_FILE,".nc",back)
    if (charnum /= len_trim(STARTUP_FILE)-2)&
         & CALL WARNING("STARUP FILE NAME does not end in .nc", &
         & trim(STARTUP_FILE))
    
    ! INITIALIZE TYPE TO HOLD FILE METADATA
    pathnfile= trim(INPUT_DIR)//trim(STARTUP_FILE)

    NCF => NEW_FILE()
    NCF%FNAME=trim(pathnfile)

    Call NC_OPEN(NCF)
    CALL NC_LOAD(NCF)
       
    NC_START => NCF
    
  END SUBROUTINE OPEN_STARTUP_FILE
!=============================================================  
  SUBROUTINE OPEN_CRASHSTART
    USE CONTROL
    IMPLICIT NONE
    TYPE(NCFILE), POINTER :: NCF
    integer :: ncfileind, datfileind,ios,charnum, i
    logical :: fexist,back,connected
    character(len=100) :: testchar
    character(len=160) :: pathnfile
    character(len=160) :: nextpathnfile
    character(len=2) :: cios

    Nullify(NCF)

    IF(RST_ON) then
       !RESTART FILE
       RESTART_FILE_NAME = trim(casename)//"_restart_0001.nc"
       pathnfile= trim(OUTPUT_DIR)//TRIM(RESTART_FILE_NAME)
       
       CALL SEARCH_FOR_LAST_MATCHING_NAME(PATHNFILE)
       
       ! OPEN THE FILE AND LOAD FOR STARTUP
       NCF => NEW_FILE()
       NCF%FNAME=trim(pathnfile)
       Call NC_OPEN(NCF)
       CALL NC_LOAD(NCF)
       NC_START => NCF
       
       Nullify(NCF)
       

       ! NOW CREATE ANOTHER FILE OBJECT FOR OUTPUT
       ! INITIALIZE TYPE TO HOLD FILE METADATA
       NCF => NEW_FILE()
       NCF%FNAME=trim(pathnfile)
       NCF%writable = .true.
       NC_RST => NCF
       FILEHEAD => ADD(FILEHEAD,NCF)
       
    END IF

    IF(NC_ON) then
       ! DATA FILE
       NC_FILE_NAME = trim(casename)//"_0001.nc"
       pathnfile= trim(OUTPUT_DIR)//TRIM(NC_FILE_NAME)
       
       CALL SEARCH_FOR_LAST_MATCHING_NAME(PATHNFILE)
       NCF => NEW_FILE()
       NCF%FNAME=trim(pathnfile)
       NCF%writable = .true.
       NC_DAT => NCF
       FILEHEAD => ADD(FILEHEAD,NCF)

    END IF
       

     IF(NCSF_ON) then
       ! DATA FILE
       NCSF_FILE_NAME = trim(casename)//"_surface_0001.nc"
       pathnfile= trim(OUTPUT_DIR)//TRIM(NCSF_FILE_NAME)
       
       CALL SEARCH_FOR_LAST_MATCHING_NAME(PATHNFILE)
       NCF => NEW_FILE()
       NCF%FNAME=trim(pathnfile)
       NCF%writable = .true.
       NC_SF => NCF
       FILEHEAD => ADD(FILEHEAD,NCF)

    END IF
       

   IF(NCAV_ON) then
       ! DATA FILE
       NCAV_FILE_NAME = trim(casename)//"_avg_0001.nc"
       pathnfile= trim(OUTPUT_DIR)//TRIM(NCAV_FILE_NAME)
       
       CALL SEARCH_FOR_LAST_MATCHING_NAME(PATHNFILE)
       NCF => NEW_FILE()
       NCF%FNAME=trim(pathnfile)
       NCF%writable = .true.
       NC_AVG => NCF
       FILEHEAD => ADD(FILEHEAD,NCF)

    END IF



  END SUBROUTINE OPEN_CRASHSTART
!=============================================================  

  ! THIS ROUTINE SEARCHES FOR THE LAST NAME SUFFIX IN THE DIRECTORY
  SUBROUTINE SEARCH_FOR_LAST_MATCHING_NAME(FNAME)
    IMPLICIT NONE
    CHARACTER(LEN=160), INTENT(INOUT) :: FNAME 
    CHARACTER(LEN=160) :: FNAME_NEXT
    logical :: fexist

    inquire(file=trim(fname),exist=fexist)
    IF(.not. Fexist) CALL FATAL_ERROR &
         & ("Base name can not be found while searching for crashrestart file:",&
         & TRIM(FNAME), "If there is no output yet a crashrestart does&
         & not make much sense...", "there is something wrong with your model")

    DO 
       FNAME_NEXT = FNAME
       CALL INCRIMENT_FNAME(FNAME_NEXT)
       ! ADD A MORE DEFINIATIVE CHECK THAN EXISTANCE OF FILE!
       inquire(file=trim(fname_next),exist=fexist)

       IF(.not. FEXIST) THEN
          if(dbg_set(dbg_log))&
               & write(ipt,*) "FOUND LAST FILE: "//TRIM(FNAME)
       
          RETURN
       ELSE
          FNAME = FNAME_NEXT
       END IF
    END DO

  END SUBROUTINE SEARCH_FOR_LAST_MATCHING_NAME
!=============================================================  
  SUBROUTINE OPEN_COLDSTART
    USE CONTROL
    IMPLICIT NONE
    TYPE(NCFILE), POINTER :: NCF
    integer :: ncfileind, datfileind,ios,charnum, i
    logical :: fexist,back,connected
    character(len=100) :: testchar
    character(len=160) :: pathnfile
    character(len=2) :: cios

    back = .true.
    
    ! TEST FILE NAME
    IF(OBC_ON) THEN
       charnum = index (OBC_NODE_LIST_FILE,".dat")
       if (charnum /= len_trim(OBC_NODE_LIST_FILE)-3)&
            & CALL WARNING("OBC NODE LIST FILE does not end in .dat", &
            & trim(OBC_NODE_LIST_FILE))
       ! OPEN FILE
       pathnfile = trim(INPUT_DIR)//trim(OBC_NODE_LIST_FILE)
       Call FOPEN(OBCUNIT,trim(pathnfile),'cfr')

       IF(OBC_LONGSHORE_FLOW_ON) THEN
          charnum = index (OBC_LONGSHORE_FLOW_FILE,".dat")
          if (charnum /= len_trim(OBC_LONGSHORE_FLOW_FILE)-3)&
               & CALL WARNING("OBC LONGSHORE FLOW FILE does not end in .dat", &
               & trim(OBC_LONGSHORE_FLOW_FILE))
          ! OPEN FILE
          pathnfile = trim(INPUT_DIR)//trim(OBC_LONGSHORE_FLOW_FILE)
          Call FOPEN(LSFUNIT,trim(pathnfile),'cfr')
       END IF


    END IF

        
    !Check Sigma File and open:
    ! TEST FILE NAME
    charnum = index (SIGMA_LEVELS_FILE,".dat")
    if (charnum /= len_trim(SIGMA_LEVELS_FILE)-3)&
         & CALL WARNING("SIGMA LEVELS FILE does not end in .dat", &
         & trim(SIGMA_LEVELS_FILE))
    ! OPEN FILE
    pathnfile = trim(INPUT_DIR)//trim(SIGMA_LEVELS_FILE)
    Call FOPEN(SIGMAUNIT,trim(pathnfile),'cfr')
    
    !Check Grid File and open:
    ! TEST FILE NAME
    charnum = index (GRID_FILE,".dat")
    if (charnum /= len_trim(GRID_FILE)-3)&
         & CALL WARNING("GRID FILE does not end in .dat", &
         & trim(GRID_FILE))
    ! OPEN FILE
    pathnfile = trim(INPUT_DIR)//trim(GRID_FILE)
    Call FOPEN(GRIDUNIT,trim(pathnfile),'cfr')
    
# if defined (PLBC) 
    !Check plbc_cell File and open:
    ! OPEN FILE
    pathnfile = trim(INPUT_DIR)//trim('plb_plbc_cell.dat')
    Call FOPEN(INTCELL1,trim(pathnfile),'cfr')

    !Check plbc_node File and open:
    ! OPEN FILE
    pathnfile = trim(INPUT_DIR)//trim('plb_plbc_node.dat')
    Call FOPEN(INTNODE1,trim(pathnfile),'cfr')
# endif
    
    !Check Depth File and open:
    ! TEST FILE NAME
    charnum = index (DEPTH_FILE,".dat")
    if (charnum /= len_trim(DEPTH_FILE)-3)&
         & CALL WARNING("DEPTH FILE does not end in .dat", &
         & trim(DEPTH_FILE))
    ! OPEN FILE
    pathnfile = trim(INPUT_DIR)//trim(DEPTH_FILE)
    Call FOPEN(DEPTHUNIT,trim(pathnfile),'cfr')
    
    !Check Sponge File and open:
    ! TEST FILE NAME
    charnum = index (SPONGE_FILE,".dat")
    if (charnum /= len_trim(SPONGE_FILE)-3)&
         & CALL WARNING("SPONGE FILE does not end in .dat", &
         & trim(SPONGE_FILE))
    ! OPEN FILE
    pathnfile = trim(INPUT_DIR)//trim(SPONGE_FILE)
    Call FOPEN(SPONGEUNIT,trim(pathnfile),'cfr')
    
    IF (GRID_FILE_UNITS == 'meters') THEN
       !Check Coriolis File and open:
       ! TEST FILE NAME
       charnum = index (CORIOLIS_FILE,".dat")
       if (charnum /= len_trim(CORIOLIS_FILE)-3)&
            & CALL WARNING("CORIOLIS FILE does not end in .dat", &
            & trim(CORIOLIS_FILE))
       ! OPEN FILE
       pathnfile = trim(INPUT_DIR)//trim(CORIOLIS_FILE)
       Call FOPEN(CORIOLISUNIT,trim(pathnfile),'cfr')
    END IF
    
  END SUBROUTINE OPEN_COLDSTART

 SUBROUTINE OPEN_NEW_OUTPUT
    USE CONTROL
    IMPLICIT NONE
    TYPE(NCFILE), POINTER :: NCF
    integer :: ncfileind, datfileind,ios,charnum, i
    logical :: fexist,back,connected
    character(len=100) :: testchar
    character(len=160) :: pathnfile
    character(len=2) :: cios


    back = .true.
    ! SETUP AND CREATE DATA OUTPUT FILES!
    IF(NC_ON) THEN
       NC_FILE_NAME = trim(casename)//"_0001.nc"
       pathnfile = trim(OUTPUT_DIR)//trim(NC_FILE_NAME)
       CALL NC_INIT(NCF,pathnfile)
       if(msr) then
          CALL NC_CREATE(NCF)
          CALL NC_CLOSE(NCF)
       else
          NCF%writable = .true.
       end if

       NC_DAT => NCF

       FILEHEAD => ADD(FILEHEAD,NCF)

    END IF
    
     IF(NCSF_ON) THEN
       NCSF_FILE_NAME = trim(casename)//"_surface_0001.nc"
       pathnfile = trim(OUTPUT_DIR)//trim(NCSF_FILE_NAME)
       CALL NC_INIT(NCF,pathnfile)
       if(msr) then
          CALL NC_CREATE(NCF)
          CALL NC_CLOSE(NCF)
       else
          NCF%writable = .true.
       end if

       NC_SF => NCF

       FILEHEAD => ADD(FILEHEAD,NCF)

    END IF
    
   IF(NCAV_ON) THEN
       NCAV_FILE_NAME = trim(casename)//"_avg_0001.nc"
       pathnfile = trim(OUTPUT_DIR)//trim(NCAV_FILE_NAME)
       CALL NC_INIT(NCF,pathnfile)
       if(msr) then
          CALL NC_CREATE(NCF)
          CALL NC_CLOSE(NCF)
       else
          NCF%writable = .true.
       end if

       NC_AVG => NCF

       FILEHEAD => ADD(FILEHEAD,NCF)
    END IF
    
    
    if(RST_ON) then
       RESTART_FILE_NAME = trim(casename)//"_restart_0001.nc"
       pathnfile = trim(OUTPUT_DIR)//trim(RESTART_FILE_NAME)
       CALL NC_INIT(NCF,pathnfile)
       if(msr) then
          CALL NC_CREATE(NCF)
          CALL NC_CLOSE(NCF)
       else
          NCF%writable = .true.
       end if

       NC_RST => NCF

       FILEHEAD => ADD(FILEHEAD,NCF)
    end if

  END SUBROUTINE OPEN_NEW_OUTPUT

  SUBROUTINE OPEN_FORCING
    USE CONTROL
#   if defined (HEATING_CALCULATED)
    USE MOD_HEATFLUX, ONLY : HEATING_CALCULATE_ON, HEATING_CALCULATE_KIND, HEATING_CALCULATE_FILE
#   endif
#   if defined (HEATING_SOLAR)
    USE MOD_SOLAR, ONLY : HEATING_SOLAR_ON, HEATING_SOLAR_KIND, HEATING_SOLAR_FILE
#   endif

    IMPLICIT NONE
    TYPE(NCFILE), POINTER :: NCF
    integer :: ncfileind, datfileind,ios,charnum, i
    logical :: fexist,back,connected
    character(len=100) :: testchar
    character(len=160) :: pathnfile
    character(len=2) :: cios

    character(len=3) :: ftype
    integer :: fid, status

    back = .true.

    ! Check air pressure file and open
    if (AIRPRESSURE_ON .and. AIRPRESSURE_KIND/= CNSTNT) then

       ! TEST FILE NAME
       charnum = index (AIRPRESSURE_FILE,".nc",back)
       if (charnum /= len_trim(AIRPRESSURE_FILE)-2)&
            & CALL WARNING("AIRPRESSURE FILE does not end in .nc", &
            & trim(AIRPRESSURE_FILE))

       ! INITIALIZE TYPE TO HOLD FILE METAData
       pathnfile = trim(INPUT_DIR)//trim(AIRPRESSURE_FILE)
       CALL NC_INIT(NCF,pathnfile)

       ! OPEN THE FILE AND LOAD METADATA
       If(.not. NCF%OPEN) then
          Call NC_OPEN(NCF)
          CALL NC_LOAD(NCF)
          FILEHEAD => ADD(FILEHEAD,NCF)
       end if

    end if

    ! Check wind stress file and open
    if (WIND_ON .and. Wind_KIND /= CNSTNT) then
       
       ! TEST FILE NAME
       charnum = index (WIND_FILE,".nc",back)
       if (charnum /= len_trim(WIND_FILE)-2)&
            & CALL WARNING("WIND FILE does not end in .nc", &
            & trim(WIND_FILE))
       
       ! INITIALIZE TYPE TO HOLD FILE METAData
       pathnfile = trim(INPUT_DIR)//trim(WIND_FILE)
       CALL NC_INIT(NCF,pathnfile)

       ! OPEN THE FILE AND LOAD METADATA       
       If(.not. NCF%OPEN) then
          Call NC_OPEN(NCF)
          CALL NC_LOAD(NCF)
          FILEHEAD => ADD(FILEHEAD,NCF)
       end if

    end if

    ! Check Heat file and open
#   if !defined (HEATING_CALCULATED) && !defined (HEATING_SOLAR)
    if (HEATING_ON .and. HEATING_KIND /= CNSTNT) then

       ! TEST FILE NAME
       charnum = index (HEATING_FILE,".nc",back)
       if (charnum /= len_trim(HEATING_FILE)-2)&
            & CALL WARNING("HEATING FILE does not end in .nc", &
            & trim(HEATING_FILE))

       ! INITIALIZE TYPE TO HOLD FILE METADATA
       pathnfile= trim(INPUT_DIR)//trim(HEATING_FILE)
       CALL  NC_INIT(NCF,pathnfile)

       ! OPEN THE FILE AND LOAD METADATA
       if(.not. NCF%OPEN) then
          Call NC_OPEN(NCF)
          CALL NC_LOAD(NCF)

          FILEHEAD => ADD(FILEHEAD,NCF)
       end if

    end if
#   endif    

#   if defined (HEATING_CALCULATED)
    if (HEATING_CALCULATE_ON .and. HEATING_CALCULATE_KIND /= CNSTNT) then

       ! TEST FILE NAME
       charnum = index (HEATING_CALCULATE_FILE,".nc",back)
       if (charnum /= len_trim(HEATING_CALCULATE_FILE)-2)&
            & CALL WARNING("HEATING FILE does not end in .nc", &
            & trim(HEATING_CALCULATE_FILE))

       ! INITIALIZE TYPE TO HOLD FILE METADATA
       pathnfile= trim(INPUT_DIR)//trim(HEATING_CALCULATE_FILE)
       CALL  NC_INIT(NCF,pathnfile)

       ! OPEN THE FILE AND LOAD METADATA
       if(.not. NCF%OPEN) then
          Call NC_OPEN(NCF)
          CALL NC_LOAD(NCF)

          FILEHEAD => ADD(FILEHEAD,NCF)
       end if

    end if
#   endif    

#   if defined (HEATING_SOLAR)
    if (HEATING_SOLAR_ON .and. HEATING_SOLAR_KIND /= CNSTNT) then

       ! TEST FILE NAME
       charnum = index (HEATING_SOLAR_FILE,".nc",back)
       if (charnum /= len_trim(HEATING_SOLAR_FILE)-2)&
            & CALL WARNING("HEATING FILE does not end in .nc", &
            & trim(HEATING_SOLAR_FILE))

       ! INITIALIZE TYPE TO HOLD FILE METADATA
       pathnfile= trim(INPUT_DIR)//trim(HEATING_SOLAR_FILE)
       CALL  NC_INIT(NCF,pathnfile)

       ! OPEN THE FILE AND LOAD METADATA
       if(.not. NCF%OPEN) then
          Call NC_OPEN(NCF)
          CALL NC_LOAD(NCF)

          FILEHEAD => ADD(FILEHEAD,NCF)
       end if

    end if
#   endif    

    ! Check Precip file and open
    if (PRECIPITATION_ON .and. PRECIPITATION_KIND /= CNSTNT) then

       ! TEST FILE NAME
       charnum = index (PRECIPITATION_FILE,".nc",back)
       if (charnum /= len_trim(PRECIPITATION_FILE)-2)&
            & CALL WARNING("PRECIPITATION FILE does not end in .nc", &
            & trim(PRECIPITATION_FILE))
       
       ! INITIALIZE TYPE TO HOLD FILE METADATA
       pathnfile= trim(INPUT_DIR)//trim(PRECIPITATION_FILE)
       CALL  NC_INIT(NCF,pathnfile)
       
       ! OPEN THE FILE AND LOAD METADATA
       if(.not. NCF%OPEN) then
          Call NC_OPEN(NCF)
          CALL NC_LOAD(NCF)
           FILEHEAD => ADD(FILEHEAD,NCF)
        end if
    end if

    ! Check Wave file and open
    if (WAVE_ON .and. WAVE_KIND /= CNSTNT) then

       ! TEST FILE NAME
       charnum = index (WAVE_FILE,".nc",back)
       if (charnum /= len_trim(WAVE_FILE)-2)&
            & CALL WARNING("WAVE FILE does not end in .nc", &
            & trim(WAVE_FILE))
       
       ! INITIALIZE TYPE TO HOLD FILE METADATA
       pathnfile= trim(INPUT_DIR)//trim(WAVE_FILE)
       CALL  NC_INIT(NCF,pathnfile)
       
       ! OPEN THE FILE AND LOAD METADATA
       if(.not. NCF%OPEN) then
          Call NC_OPEN(NCF)
          CALL NC_LOAD(NCF)
           FILEHEAD => ADD(FILEHEAD,NCF)
        end if
    end if

# if defined (OFFLINE_SEDIMENT)
    ! Check Wave file and open
 
       ! TEST FILE NAME
       charnum = index (OFFLINE_SEDIMENT_FILE,".nc",back)
       if (charnum /= len_trim(OFFLINE_SEDIMENT_FILE)-2)&
            & CALL WARNING("OFFLINE SEDIMENT FILE does not end in .nc", &
            & trim(OFFLINE_SEDIMENT_FILE))
       
       ! INITIALIZE TYPE TO HOLD FILE METADATA
       pathnfile= trim(INPUT_DIR)//trim(OFFLINE_SEDIMENT_FILE)
       CALL  NC_INIT(NCF,pathnfile)
       
       ! OPEN THE FILE AND LOAD METADATA
       if(.not. NCF%OPEN) then
          Call NC_OPEN(NCF)
          CALL NC_LOAD(NCF)
           FILEHEAD => ADD(FILEHEAD,NCF)
        end if

# endif


# if defined (OFFLINE_BIOLOGY)
    ! Check Wave file and open
 
       ! TEST FILE NAME
       charnum = index (OFFLINE_BIOLOGY_FILE,".nc",back)
       if (charnum /= len_trim(OFFLINE_BIOLOGY_FILE)-2)&
            & CALL WARNING("OFFLINE BIOLOGY FILE does not end in .nc", &
            & trim(OFFLINE_BIOLOGY_FILE))
       
       ! INITIALIZE TYPE TO HOLD FILE METADATA
       pathnfile= trim(INPUT_DIR)//trim(OFFLINE_BIOLOGY_FILE)
       CALL  NC_INIT(NCF,pathnfile)
       
       ! OPEN THE FILE AND LOAD METADATA
       if(.not. NCF%OPEN) then
          Call NC_OPEN(NCF)
          CALL NC_LOAD(NCF)
           FILEHEAD => ADD(FILEHEAD,NCF)
        end if

# endif

    !Check RIVER files and open
    
    do i =1, river_number

       ! TEST FILE NAME
       charnum = index (RIVERS(i)%FILE,".nc",back)
       if (charnum /= len_trim(RIVERS(i)%FILE)-2)&
            & CALL WARNING("RIVER FILE does not end in .nc", &
            & trim(RIVERS(i)%FILE))
       
       ! INITIALIZE TYPE TO HOLD FILE METADATA
       pathnfile= trim(INPUT_DIR)//trim(RIVERS(i)%FILE)
       CALL  NC_INIT(NCF,pathnfile)
       
       ! OPEN THE FILE AND LOAD METADATA IF NOT ALREADY DONE
       if(.not. NCF%OPEN) then
          Call NC_OPEN(NCF)
          CALL NC_LOAD(NCF)
          FILEHEAD => ADD(FILEHEAD,NCF)
       end if

    end do

    
    !Check OBC files and open:

# if !defined (EQUI_TIDE)
    if (OBC_ON) then

       IF(OBC_ELEVATION_FORCING_ON) THEN
# endif

          ! Determine file type:
          pathnfile= trim(INPUT_DIR)//trim(OBC_ELEVATION_FILE)
          

          ! TRY opening file as a netcdf file
          IF(DBG_SET(DBG_LOG))write(ipt,*) "! Trying to open Boundary Forcing file: "//TRIM(OBC_ELEVATION_FILE)
          
          status = nf90_open(trim(pathnfile), NF90_NOWRITE, fid)
          if(status == nf90_noerr) then
             status = nf90_close(fid)

             IF(DBG_SET(DBG_LOG))write(ipt,*) "! Open Boundary Forcing file is a NETCDF FILE"

             ! INITIALIZE TYPE TO HOLD FILE METADATA
             CALL  NC_INIT(NCF,pathnfile)
             
             ! OPEN THE FILE AND LOAD METADATA
             if(.not. NCF%OPEN) then
                Call NC_OPEN(NCF)
                CALL NC_LOAD(NCF)
                FILEHEAD => ADD(FILEHEAD,NCF)
             end if
          else

             IF(DBG_SET(DBG_LOG)) write(ipt,*) "! Open Boundary Forcing file is not a NETCDF file"
             IF(MSR) THEN
                
                Call FOPEN(JULOBCUNIT,trim(pathnfile),'cfr')

                write(ipt,*) "! Open Boundary Forcing file is an ASCII file"
             END IF

             ! INITIALIZE DUMMY FILE TYPE FOR ASCII FILE
             NCF => NEW_FILE(pathnfile)
             NCF => ADD(NCF,NC_MAKE_ATT("type","ASCII FILE DUMMY ATTRIBUTE"))
             FILEHEAD => ADD(FILEHEAD,NCF)
             
          end if
          
# if !defined (EQUI_TIDE)
       END IF
# else
   if (OBC_ON) then
# endif      
      
      
      if (OBC_TEMP_NUDGING) then
         
         ! TEST FILE NAME
         charnum = index (OBC_TEMP_FILE,".nc",back)
         if (charnum /= len_trim(OBC_TEMP_FILE)-2)&
              & CALL WARNING("OBC TEMP FILE does not end in .nc", &
              & trim(OBC_TEMP_FILE))
         
         ! INITIALIZE TYPE TO HOLD FILE METADATA
         pathnfile= trim(INPUT_DIR)//trim(OBC_TEMP_FILE)
         CALL  NC_INIT(NCF,pathnfile)
         
         ! OPEN THE FILE AND LOAD METADATA
         if(.not. NCF%OPEN)then 
            Call NC_OPEN(NCF)
            CALL NC_LOAD(NCF)
            FILEHEAD => ADD(FILEHEAD,NCF)
         end if
      end if
      
      if (OBC_SALT_NUDGING) then
         
         ! TEST FILE NAME
         charnum = index (OBC_SALT_FILE,".nc",back)
         if (charnum /= len_trim(OBC_SALT_FILE)-2)&
              & CALL WARNING("OBC SALT FILE does not end in .nc", &
              & trim(OBC_SALT_FILE))
         
         ! INITIALIZE TYPE TO HOLD FILE METADATA
         pathnfile= trim(INPUT_DIR)//trim(OBC_SALT_FILE)
         CALL  NC_INIT(NCF,pathnfile)
         
         ! OPEN THE FILE AND LOAD METADATA
         if(.not. NCF%OPEN) then
            Call NC_OPEN(NCF)
            CALL NC_LOAD(NCF)
            FILEHEAD => ADD(FILEHEAD,NCF)
         end if
      end if
      
      if (OBC_MEANFLOW) then
         
         ! TEST FILE NAME
         charnum = index (OBC_MEANFLOW_FILE,".nc",back)
         if (charnum /= len_trim(OBC_MEANFLOW_FILE)-2)&
              & CALL WARNING("OBC MEANFLOW FILE does not end in .nc", &
              & trim(OBC_MEANFLOW_FILE))
         
         ! INITIALIZE TYPE TO HOLD FILE METADATA
         pathnfile= trim(INPUT_DIR)//trim(OBC_MEANFLOW_FILE)
         CALL  NC_INIT(NCF,pathnfile)
         
         ! OPEN THE FILE AND LOAD METADATA
         if(.not. NCF%OPEN) then
            Call NC_OPEN(NCF)
            CALL NC_LOAD(NCF)
            FILEHEAD => ADD(FILEHEAD,NCF)
         end if
      end if
      
      
   end if
   

    !Check Ground Water File and open:
    if (GROUNDWATER_ON .and. GROUNDWATER_KIND /= CNSTNT) then

       ! TEST FILE NAME
       charnum = index (GROUNDWATER_FILE,".nc",back)
       if (charnum /= len_trim(GROUNDWATER_FILE)-2)&
            & CALL WARNING("GROUNDWATER FILE does not end in .nc", &
            & trim(GROUNDWATER_FILE))
       
       ! INITIALIZE TYPE TO HOLD FILE METADATA
       pathnfile= trim(INPUT_DIR)//trim(GROUNDWATER_FILE)
       CALL  NC_INIT(NCF,pathnfile)
       
       ! OPEN THE FILE AND LOAD METADATA
       if(.not. NCF%OPEN) then
          Call NC_OPEN(NCF)
          CALL NC_LOAD(NCF)
          FILEHEAD => ADD(FILEHEAD,NCF)
       end if
    end if


    IF (ICING_MODEL .and. ICING_FORCING_KIND /= CNSTNT) THEN

       ! TEST FILE NAME
       charnum = index (ICING_FORCING_FILE,".nc",back)
       if (charnum /= len_trim(ICING_FORCING_FILE)-2)&
            & CALL WARNING("ICING MODEL FILE does not end in .nc", &
            & trim(ICING_FORCING_FILE))
       
       ! INITIALIZE TYPE TO HOLD FILE METADATA
       pathnfile= trim(INPUT_DIR)//trim(ICING_FORCING_FILE)
       CALL  NC_INIT(NCF,pathnfile)
       
       ! OPEN THE FILE AND LOAD METADATA
       if(.not. NCF%OPEN) then
          Call NC_OPEN(NCF)
          CALL NC_LOAD(NCF)
          FILEHEAD => ADD(FILEHEAD,NCF)
       end if

    END IF

    IF (ICE_MODEL .and. ICE_FORCING_KIND /= CNSTNT) THEN

       ! TEST FILE NAME
       charnum = index (ICE_FORCING_FILE,".nc",back)
       if (charnum /= len_trim(ICE_FORCING_FILE)-2)&
            & CALL WARNING("ICE MODEL FILE does not end in .nc", &
            & trim(ICE_FORCING_FILE))
       
       ! INITIALIZE TYPE TO HOLD FILE METADATA
       pathnfile= trim(INPUT_DIR)//trim(ICE_FORCING_FILE)
       CALL  NC_INIT(NCF,pathnfile)
       
       ! OPEN THE FILE AND LOAD METADATA
       if(.not. NCF%OPEN) then
          Call NC_OPEN(NCF)
          CALL NC_LOAD(NCF)
          FILEHEAD => ADD(FILEHEAD,NCF)
       end if
       
    END IF

    ! LOAD HORIZONTAL MIXING FILE
    IF (HORIZONTAL_MIXING_KIND /= CNSTNT) THEN

       ! TEST FILE NAME
       charnum = index (HORIZONTAL_MIXING_FILE,".nc",back)
       if (charnum /= len_trim(HORIZONTAL_MIXING_FILE)-2)&
            & CALL WARNING("Horizontal Mixing File does not end in .nc", &
            & trim(HORIZONTAL_MIXING_FILE))
       
       ! INITIALIZE TYPE TO HOLD FILE METADATA
       pathnfile= trim(INPUT_DIR)//trim(HORIZONTAL_MIXING_FILE)
       CALL  NC_INIT(NCF,pathnfile)
       
       ! OPEN THE FILE AND LOAD METADATA
       if(.not. NCF%OPEN) then
          Call NC_OPEN(NCF)
          CALL NC_LOAD(NCF)
          FILEHEAD => ADD(FILEHEAD,NCF)
       end if
       
    END IF


    ! LOAD BOTTOM ROUGHNESS LENGTH SCALE
    IF (BOTTOM_ROUGHNESS_KIND /= CNSTNT) THEN

       ! TEST FILE NAME
       charnum = index (BOTTOM_ROUGHNESS_FILE,".nc",back)
       if (charnum /= len_trim(BOTTOM_ROUGHNESS_FILE)-2)&
            & CALL WARNING("Bottom Roughness File does not end in .nc", &
            & trim(BOTTOM_ROUGHNESS_FILE))
       
       ! INITIALIZE TYPE TO HOLD FILE METADATA
       pathnfile= trim(INPUT_DIR)//trim(BOTTOM_ROUGHNESS_FILE)
       CALL  NC_INIT(NCF,pathnfile)
       
       ! OPEN THE FILE AND LOAD METADATA
       if(.not. NCF%OPEN) then
          Call NC_OPEN(NCF)
          CALL NC_LOAD(NCF)
          FILEHEAD => ADD(FILEHEAD,NCF)
       end if
       
    END IF

# if defined (SEDIMENT)
    !LOAD BEDFLAG FLAG
    IF (BEDFLAG_TYPE /= CNSTNT) THEN

       ! TEST FILE NAME
       charnum = index (BEDFLAG_FILE,".nc",back)
       if (charnum /= len_trim(BEDFLAG_FILE)-2)&
            & CALL WARNING("BedFlag file does not end in .nc", &
            & trim(BEDFLAG_FILE)) 

       ! INITIALIZE TYPE TO HOLD FILE METADATA
       pathnfile= trim(INPUT_DIR)//trim(BEDFLAG_FILE) 
       CALL  NC_INIT(NCF,pathnfile)

       ! OPEN THE FILE AND LOAD METADATA
       if(.not. NCF%OPEN) then
          Call NC_OPEN(NCF)
          CALL NC_LOAD(NCF)
          FILEHEAD => ADD(FILEHEAD,NCF)
       end if

    END IF

    ! LOAD SEDIMENT PARAMETER
    IF (SEDIMENT_PARAMETER_TYPE /= UNIFORM) THEN

       ! TEST FILE NAME
       charnum = index (SEDIMENT_PARAMETER_FILE,".nc",back)
       if (charnum /= len_trim(SEDIMENT_PARAMETER_FILE)-2)&
            & CALL WARNING("Sediment parameter file does not end in .nc", &
            & trim(SEDIMENT_PARAMETER_FILE))
       
       ! INITIALIZE TYPE TO HOLD FILE METADATA
       pathnfile= trim(INPUT_DIR)//trim(SEDIMENT_PARAMETER_FILE)
       CALL  NC_INIT(NCF,pathnfile)
       
       ! OPEN THE FILE AND LOAD METADATA
       if(.not. NCF%OPEN) then
          Call NC_OPEN(NCF)
          CALL NC_LOAD(NCF)
          FILEHEAD => ADD(FILEHEAD,NCF)
       end if
       
    END IF
# endif

  END SUBROUTINE OPEN_FORCING

  SUBROUTINE LOAD_HORIZONTAL_MIXING_COEFFICIENT(NN,CC)
    USE CONTROL
    IMPLICIT NONE
    REAL(SP),ALLOCATABLE :: NN(:),CC(:)
    TYPE(NCFILE), POINTER :: NCF
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCDIM),  POINTER :: DIM1
    TYPE(NCDIM),  POINTER :: DIM2
    integer status,I,IERR
    
    LOGICAL FOUND

    ! FIND THE HVC FILE OBJECT
    NCF => FIND_FILE(FILEHEAD,trim(HORIZONTAL_MIXING_FILE),FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR &
         & ("COULD NOT FIND HORIZONTAL_MIXING_FILE FILE OBJECT",&
         & "FILE NAME: "//TRIM(HORIZONTAL_MIXING_FILE))

    DIM1 => FIND_DIM(NCF,'nele',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR &
         & ("COULD NOT FIND HORIZONTAL_MIXING_FILE DIMENSION 'nele' in:",&
         & "FILE NAME: "//TRIM(HORIZONTAL_MIXING_FILE))
    IF (DIM1%DIM /= NGL)CALL FATAL_ERROR &
         & ("Dimension 'nele' in the HORIZONTAL_MIXING_FILE does not match NGL for this model?",&
         & "FILE NAME: "//TRIM(HORIZONTAL_MIXING_FILE)) 

    DIM2 => FIND_DIM(NCF,'node',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR &
         & ("COULD NOT FIND HORIZONTAL_MIXING_FILE DIMENSION 'node' in:",&
         & "FILE NAME: "//TRIM(HORIZONTAL_MIXING_FILE))
    IF (DIM2%DIM /= MGL)CALL FATAL_ERROR &
         & ("Dimension 'node' in the HORIZONTAL_MIXING_FILE does not match MGL for this model?",&
         & "FILE NAME: "//TRIM(HORIZONTAL_MIXING_FILE)) 

    ! FIND THE 'nn_hvc' variable
    VAR => FIND_VAR(NCF,'nn_hvc',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR &
         & ("COULD NOT FIND HORIZONTAL_MIXING_FILE VARIABLE 'nn_hvc' in:",&
         & "FILE NAME: "//TRIM(HORIZONTAL_MIXING_FILE))

    CALL NC_CONNECT_AVAR(VAR,NN)
    CALL NC_READ_VAR(VAR)
    CALL NC_DISCONNECT(VAR)

    ! FIND THE 'cc_hvc' variable
    VAR => FIND_VAR(NCF,'cc_hvc',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR &
         & ("COULD NOT FIND HORIZONTAL_MIXING_FILE VARIABLE 'cc_hvc' in:",&
         & "FILE NAME: "//TRIM(HORIZONTAL_MIXING_FILE))

    CALL NC_CONNECT_AVAR(VAR,CC)
    CALL NC_READ_VAR(VAR)
    CALL NC_DISCONNECT(VAR)


  END SUBROUTINE LOAD_HORIZONTAL_MIXING_COEFFICIENT

  SUBROUTINE LOAD_BOTTOM_ROUGHNESS(Z0)
    USE CONTROL
    IMPLICIT NONE
    REAL(SP),ALLOCATABLE :: Z0(:)
    TYPE(NCFILE), POINTER :: NCF
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCDIM),  POINTER :: DIM1
    TYPE(NCDIM),  POINTER :: DIM2
    integer status,I,IERR
    
    LOGICAL FOUND

    ! FIND THE Bottom Roughness FILE OBJECT
    NCF => FIND_FILE(FILEHEAD,trim(BOTTOM_ROUGHNESS_FILE),FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR &
         & ("COULD NOT FIND BOTTOM_ROUGHNESS_FILE FILE OBJECT",&
         & "FILE NAME: "//TRIM(BOTTOM_ROUGHNESS_FILE))

    DIM1 => FIND_DIM(NCF,'nele',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR &
         & ("COULD NOT FIND BOTTOM_ROUGHNESS_FILE DIMENSION 'nele' in:",&
         & "FILE NAME: "//TRIM(BOTTOM_ROUGHNESS_FILE))
    IF (DIM1%DIM /= NGL)CALL FATAL_ERROR &
         & ("Dimension 'nele' in the BOTTOM_ROUGHNESS_FILE does not match NGL for this model?",&
         & "FILE NAME: "//TRIM(BOTTOM_ROUGHNESS_FILE)) 


    ! FIND THE 'z0B' variable
    VAR => FIND_VAR(NCF,'z0b',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR &
         & ("COULD NOT FIND BOTTOM_ROUGHNESS_FILE VARIABLE 'z0b' in:",&
         & "FILE NAME: "//TRIM(BOTTOM_ROUGHNESS_FILE))

    CALL NC_CONNECT_AVAR(VAR,Z0)
    CALL NC_READ_VAR(VAR)
    CALL NC_DISCONNECT(VAR)

    
  END SUBROUTINE LOAD_BOTTOM_ROUGHNESS

# if defined (SEDIMENT)

  SUBROUTINE LOAD_BEDFLAG(BEDFLAG)
    USE CONTROL
    IMPLICIT NONE
    REAL(SP),ALLOCATABLE :: BEDFLAG(:)
    TYPE(NCFILE), POINTER :: NCF
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCDIM),  POINTER :: DIM1
    TYPE(NCDIM),  POINTER :: DIM2
    integer status,I,IERR

    LOGICAL FOUND

    ! FIND THE BedFlag FILE OBJECT
    NCF => FIND_FILE(FILEHEAD,trim(BEDFLAG_FILE),FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR &
         & ("COULD NOT FIND BEDFLAG FILE OBJECT",&
         & "FILE NAME: "//TRIM(BEDFLAG_FILE)) 

    DIM1 => FIND_DIM(NCF,'node',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR &
         & ("COULD NOT FIND BEDFLAG_FILE DIMENSION 'node' in:",&
         & "FILE NAME: "//TRIM(BEDFLAG_FILE))
    IF (DIM1%DIM /= MGL)CALL FATAL_ERROR &
         & ("Dimension 'node' in the BEDFLAG_FILE does not match MGL for this model?",&
         & "FILE NAME: "//TRIM(BEDFLAG_FILE))


    ! FIND THE 'beflag' variable
    VAR => FIND_VAR(NCF,'bedflag',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR &
         & ("COULD NOT FIND BEDFLAG_FILE VARIABLE 'bedflag' in:",&
         & "FILE NAME: "//TRIM(BEDFLAG_FILE))

    CALL NC_CONNECT_AVAR(VAR,BEDFLAG)
    CALL NC_READ_VAR(VAR)
    CALL NC_DISCONNECT(VAR)


  END SUBROUTINE LOAD_BEDFLAG


  SUBROUTINE LOAD_SEDIMENT_PARAMETER(t_cd,t_ce,erate)
    USE CONTROL
    IMPLICIT NONE
    REAL(SP),ALLOCATABLE :: t_cd(:),t_ce(:),erate(:)
    TYPE(NCFILE), POINTER :: NCF
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCDIM),  POINTER :: DIM1
    TYPE(NCDIM),  POINTER :: DIM2
    integer status,I,IERR
    
    LOGICAL FOUND

    ! FIND THE Sediment Parameter FILE OBJECT
    NCF => FIND_FILE(FILEHEAD,trim(SEDIMENT_PARAMETER_FILE),FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR &
         & ("COULD NOT FIND SEDIMENT_PARAMETER_FILE OBJECT",&
         & "FILE NAME: "//TRIM(SEDIMENT_PARAMETER_FILE))

    DIM1 => FIND_DIM(NCF,'node',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR &
         & ("COULD NOT FIND SEDIMENT_PARAMETER_FILE DIMENSION 'node' in:",&
         & "FILE NAME: "//TRIM(SEDIMENT_PARAMETER_FILE))
    IF (DIM1%DIM /= MGL)CALL FATAL_ERROR &
         & ("Dimension 'node' in the SEDIMENT_PARAMETER_FILE does not match MGL for this model?",&
         & "FILE NAME: "//TRIM(SEDIMENT_PARAMETER_FILE)) 


    ! FIND THE 'tau_cd' variable
    VAR => FIND_VAR(NCF,'tau_cd',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR &
         & ("COULD NOT FIND SEDIMENT_PARAMETER_FILE VARIABLE 'tau_cd' in:",&
         & "FILE NAME: "//TRIM(SEDIMENT_PARAMETER_FILE))

    CALL NC_CONNECT_AVAR(VAR,t_cd)
    CALL NC_READ_VAR(VAR)
    CALL NC_DISCONNECT(VAR)

    ! FIND THE 'tau_ce' variable
    VAR => FIND_VAR(NCF,'tau_ce',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR &
         & ("COULD NOT FIND SEDIMENT_PARAMETER_FILE VARIABLE 'tau_ce' in:",&
         & "FILE NAME: "//TRIM(SEDIMENT_PARAMETER_FILE))

    CALL NC_CONNECT_AVAR(VAR,t_ce)
    CALL NC_READ_VAR(VAR)
    CALL NC_DISCONNECT(VAR)

    ! FIND THE 'tau_cd' variable
    VAR => FIND_VAR(NCF,'erate',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR &
         & ("COULD NOT FIND SEDIMENT_PARAMETER_FILE VARIABLE 'erate' in:",&
         & "FILE NAME: "//TRIM(SEDIMENT_PARAMETER_FILE))

    CALL NC_CONNECT_AVAR(VAR,erate)
    CALL NC_READ_VAR(VAR)
    CALL NC_DISCONNECT(VAR)

    
  END SUBROUTINE LOAD_SEDIMENT_PARAMETER

# endif


  SUBROUTINE LOAD_GRID_TYPE(NCF,G)
    IMPLICIT NONE

    TYPE(NCFILE), POINTER :: NCF
    TYPE(GRID), POINTER :: G
    
    
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCDIM),  POINTER :: DIM1
    TYPE(NCDIM),  POINTER :: DIM2
    integer status,I,IERR
    
    LOGICAL FOUND
    
    !==============================================================================|
    
    
    ! GET GLOBAL DIMENSION!
    DIM1 => FIND_DIM(NCF,"node",FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR&
         &("COULD NOT FIND DIMENSION 'node' IN THE FILE OBJECT:"//TRIM(NCF%FNAME))
    
    IF(G%MGL==0) THEN
       G%MGL = DIM1%DIM
    ELSEIF(G%MGL/= DIM1%DIM)THEN
       CALL FATAL_ERROR &
            &("THE GRID TYPE DIMENSION MGL DOES NOT MATCH THE FILE DIMENSION:"//TRIM(NCF%FNAME))
    END IF

    DIM1 => FIND_DIM(NCF,"nele",FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR&
         &("COULD NOT FIND DIMENSION 'nele' IN THE FILE OBJECT:"//TRIM(NCF%FNAME))

    IF(G%NGL==0) THEN
       G%NGL = DIM1%DIM
    ELSEIF(G%NGL/= DIM1%DIM)THEN
       CALL FATAL_ERROR &
            &("THE GRID TYPE DIMENSION NGL DOES NOT MATCH THE FILE DIMENSION:"//TRIM(NCF%FNAME))
    END IF

    IF(G%MT == 0) G%MT = G%MGL
    IF(G%NT == 0) G%NT = G%NGL

    
    ! CHECK TO MAKE SURE DIMENSION THREE IS THERE
    DIM1 => FIND_DIM(NCF,"three",FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR&
         &("COULD NOT FIND DIMENSION 'three' IN THE FILE OBJECT:"//TRIM(NCF%FNAME))

    IF(DIM1%DIM /=3) CALL FATAL_ERROR&
         &("DIMENSION 'three' IS NOT 3 IN THE FILE OBJECT:"//TRIM(NCF%FNAME))
    

    ! SIGLEV/Z
    DIM1 => FIND_DIM(NCF,"siglev",FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR&
         &("COULD NOT FIND DIMENSION 'siglev' IN THE FILE OBJECT:"//TRIM(NCF%FNAME))
    
    IF(G%KB==0) THEN
       G%KB = DIM1%DIM
    ELSEIF(G%KB/= DIM1%DIM)THEN
       CALL FATAL_ERROR &
            &("THE GRID TYPE DIMENSION KB DOES NOT MATCH THE FILE DIMENSION:"//TRIM(NCF%FNAME))
    END IF
    
    
    !SIGLAY/ZZ
    DIM1 => FIND_DIM(NCF,"siglay",FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR&
         &("COULD NOT FIND DIMENSION 'siglay' IN THE FILE OBJECT:"//TRIM(NCF%FNAME))
    
    G%KBM1 = DIM1%DIM
    IF(G%KBM1 /= G%KB -1) CALL FATAL_ERROR&
         &("KB and KBM1 DO NOT MATCH IN THE FILE:"//TRIM(NCF%FNAME))
    
    G%KBM2 = G%KB - 2
    

    IF(.NOT. IOPROC)THEN
       
       ! READ THE GRID DATA
       ALLOCATE(G%NV(0:G%NGL,4),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE G%NV WHILE READING:"//TRIM(NCF%FNAME))
       G%nv=0
       
       VAR => FIND_VAR(NCF,'nv',FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR&
            &("COULD NOT FIND VARIABLE 'nv' IN THE FILE OBJECT:"//TRIM(NCF%FNAME))
       
       VAR%ARR_INT => G%nv(1:G%NGL,1:3)
       CALL NC_READ_VAR(VAR)
       CALL NC_DISCONNECT(VAR)
       
       G%NV(:,4)=G%NV(:,1)
       
       
       ! READ THE COORDINATES
       ! X - meters
       ALLOCATE(G%XM(0:G%MT),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE G%XM WHILE READING:"//TRIM(NCF%FNAME))
       G%XM=0.0_SP
       
       VAR => FIND_VAR(NCF,'x',FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR&
            &("COULD NOT FIND VARIABLE 'x' IN THE FILE OBJECT:"//TRIM(NCF%FNAME))
       
       CALL NC_CONNECT_PVAR(VAR,G%XM)
       CALL NC_READ_VAR(VAR)
       CALL NC_DISCONNECT(VAR)
       
       ! Y - meters
       ALLOCATE(G%YM(0:G%MT),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE G%YM WHILE READING:"//TRIM(NCF%FNAME))
       G%YM=0.0_SP
       
       VAR => FIND_VAR(NCF,'y',FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR&
            &("COULD NOT FIND VARIABLE 'y' IN THE FILE OBJECT:"//TRIM(NCF%FNAME))
       
       CALL NC_CONNECT_PVAR(VAR,G%YM)
       CALL NC_READ_VAR(VAR)
       CALL NC_DISCONNECT(VAR)
       
       ! X element - meters
       ALLOCATE(G%XMC(0:G%NT),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE G%XMC WHILE READING:"//TRIM(NCF%FNAME))
       G%XMC=0.0_SP
       
       VAR => FIND_VAR(NCF,'xc',FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR&
            &("COULD NOT FIND VARIABLE 'xc' IN THE FILE OBJECT:"//TRIM(NCF%FNAME))
       
       CALL NC_CONNECT_PVAR(VAR,G%XMC)
       CALL NC_READ_VAR(VAR)
       CALL NC_DISCONNECT(VAR)
       
       ! Y element - meters
       ALLOCATE(G%YMC(0:G%NT),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE G%YMC WHILE READING:"//TRIM(NCF%FNAME))
       G%YMC=0.0_SP

       VAR => FIND_VAR(NCF,'yc',FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR&
            &("COULD NOT FIND VARIABLE 'yc' IN THE FILE OBJECT:"//TRIM(NCF%FNAME))
       
       CALL NC_CONNECT_PVAR(VAR,G%YMC)
       CALL NC_READ_VAR(VAR)
       CALL NC_DISCONNECT(VAR)
       
       
       ! Longitude
       ALLOCATE(G%LON(0:G%MT),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE G%LON WHILE READING:"//TRIM(NCF%FNAME))
       G%LON=0.0_SP
       
       VAR => FIND_VAR(NCF,'lon',FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR&
            &("COULD NOT FIND VARIABLE 'lon' IN THE FILE OBJECT:"//TRIM(NCF%FNAME))
       
       CALL NC_CONNECT_PVAR(VAR,G%LON)
       CALL NC_READ_VAR(VAR)
       CALL NC_DISCONNECT(VAR)
       
       ! Latitude
       ALLOCATE(G%LAT(0:G%MT),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE G%LAT WHILE READING:"//TRIM(NCF%FNAME))
       G%LAT=0.0_SP
       
       VAR => FIND_VAR(NCF,'lat',FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR&
            &("COULD NOT FIND VARIABLE 'lat' IN THE FILE OBJECT:"//TRIM(NCF%FNAME))
       
       CALL NC_CONNECT_PVAR(VAR,G%LAT)
       CALL NC_READ_VAR(VAR)
       CALL NC_DISCONNECT(VAR)
       
       ! Longitude - element
       ALLOCATE(G%LONC(0:G%NT),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE G%LONC WHILE READING:"//TRIM(NCF%FNAME))
       G%LONC=0.0_SP
       
       VAR => FIND_VAR(NCF,'lonc',FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR&
            &("COULD NOT FIND VARIABLE 'lonc' IN THE FILE OBJECT:"//TRIM(NCF%FNAME))
       
       CALL NC_CONNECT_PVAR(VAR,G%LONC)
       CALL NC_READ_VAR(VAR)
       CALL NC_DISCONNECT(VAR)
       
       ! Latitude - element
       ALLOCATE(G%LATC(0:G%NT),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE G%LATC WHILE READING:"//TRIM(NCF%FNAME))
       G%LATC=0.0_SP
       
       VAR => FIND_VAR(NCF,'latc',FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR&
            &("COULD NOT FIND VARIABLE 'latc' IN THE FILE OBJECT:"//TRIM(NCF%FNAME))
       
       CALL NC_CONNECT_PVAR(VAR,G%LATC)
       CALL NC_READ_VAR(VAR)
       CALL NC_DISCONNECT(VAR)
       
       
       ! h
       ALLOCATE(G%H(0:G%MT),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE G%H WHILE READING:"//TRIM(NCF%FNAME))
       G%H=0.0_SP
       
       VAR => FIND_VAR(NCF,'h',FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR&
            &("COULD NOT FIND VARIABLE 'h' IN THE FILE OBJECT:"//TRIM(NCF%FNAME))
       
       CALL NC_CONNECT_PVAR(VAR,G%H)
       CALL NC_READ_VAR(VAR)
       CALL NC_DISCONNECT(VAR)
       
       ! h1
       ALLOCATE(G%H1(0:G%NT),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE G%H1 WHILE READING:"//TRIM(NCF%FNAME))
       G%H1=0.0_SP
       
       VAR => FIND_VAR(NCF,'h_center',FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR&
            &("COULD NOT FIND VARIABLE 'h_center' IN THE FILE OBJECT:"//TRIM(NCF%FNAME))
       
       CALL NC_CONNECT_PVAR(VAR,G%H1)
       CALL NC_READ_VAR(VAR)
       CALL NC_DISCONNECT(VAR)
       
       ! zz
       ALLOCATE(G%ZZ(0:G%MT,G%KBM1),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE G%ZZ WHILE READING:"//TRIM(NCF%FNAME))
       G%ZZ=0.0_SP
       
       VAR => FIND_VAR(NCF,'siglay',FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR&
            &("COULD NOT FIND VARIABLE 'siglay' IN THE FILE OBJECT:"//TRIM(NCF%FNAME))
       
       CALL NC_CONNECT_PVAR(VAR,G%ZZ)
       CALL NC_READ_VAR(VAR)
       CALL NC_DISCONNECT(VAR)
       
       ! zz1
       ALLOCATE(G%ZZ1(0:G%NT,G%KBM1),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE G%ZZ1 WHILE READING:"//TRIM(NCF%FNAME))
       G%ZZ1=0.0_SP
       
       VAR => FIND_VAR(NCF,'siglay_center',FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR&
            &("COULD NOT FIND VARIABLE 'siglay_center' IN THE FILE OBJECT:"//TRIM(NCF%FNAME))
       
       CALL NC_CONNECT_PVAR(VAR,G%ZZ1)
       CALL NC_READ_VAR(VAR)
       CALL NC_DISCONNECT(VAR)
       
    END IF
       
       
    
  END SUBROUTINE LOAD_GRID_TYPE



  SUBROUTINE LOAD_RESTART_GRID(NVG)
    USE CONTROL
    IMPLICIT NONE
    INTEGER, ALLOCATABLE, TARGET, INTENT(OUT) :: NVG(:,:)
    
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCDIM),  POINTER :: DIM1
    TYPE(NCDIM),  POINTER :: DIM2
    integer status,I,IERR
    
    LOGICAL FOUND
    
    !==============================================================================|
    
    
    ! FIND THE RESTART FILE OBJECT
    
    
    ! GET GLOBAL DIMENSION!
    DIM1 => FIND_DIM(NC_START,"node",FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND DIMENSION NODE&
         & IN THE HOTSTART FILE OBJECT")
    
    MGL = DIM1%DIM
    
    DIM1 => FIND_DIM(NC_START,"nele",FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND DIMENSION NODE&
         & IN THE HOTSTART FILE OBJECT")
    
    NGL= DIM1%DIM
    
    ! CHECK TO MAKE SURE DIMENSION THREE IS THERE
    DIM1 => FIND_DIM(NC_START,"three",FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND DIMENSION NODE&
         & IN THE HOTSTART FILE OBJECT")
    
    DIM1 => FIND_DIM(NC_START,"siglev",FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND DIMENSION SIGLEV&
         & IN THE HOTSTART FILE OBJECT")
    
    KB = DIM1%DIM

    DIM1 => FIND_DIM(NC_START,"siglay",FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND DIMENSION SIGLEV&
         & IN THE HOTSTART FILE OBJECT")
    
    KBM1 = DIM1%DIM

    KBM2 = KB - 2
    

    ALLOCATE(NVG(0:NGL,4),stat=status)
    IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE YG_GRD")
    nvg=0

    VAR => FIND_VAR(NC_START,'nv',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'nv'&
         & IN THE HOTSTART FILE OBJECT")
    ! THIS IS DANGEROUS BUT IT WORKS... IN GENERAL IT IS NOT A GOOD
    ! IDEA TO PASS AROUND NON-CONTIGUOUS MEMORY POINTERS!
    VAR%ARR_INT => nvg(1:NGL,1:3)

    CALL NC_READ_VAR(VAR)
    CALL NC_DISCONNECT(VAR)
    
    NVG(:,4)=NVG(:,1)
    
    if(DBG_SET(dbg_log)) then
       WRITE(IPT,*)'!  Finished Reading Grid from HOTSTART'
       WRITE(IPT,*)'!  # OF NODES            :',MGL
       WRITE(IPT,*)'!  # OF CELLS            :',NGL
       WRITE(IPT,*)'!  # OF LEVELS           :',KB
       WRITE(IPT,*)'!'
    end if
    
  END SUBROUTINE LOAD_RESTART_GRID


 SUBROUTINE LOAD_RESTART_OBC_GRID(IOBCN_GL,I_OBC_GL,TYPE_OBC_GL)
   USE CONTROL
   IMPLICIT NONE
   INTEGER, INTENT(OUT)  :: IOBCN_GL
   INTEGER, INTENT(OUT), Allocatable, TARGET :: I_OBC_GL(:), TYPE_OBC_GL(:)
   TYPE(NCVAR),  POINTER :: VAR
   TYPE(NCDIM),  POINTER :: DIM
   LOGICAL :: FOUND


   IF(.NOT. OBC_ON) THEN
      if(DBG_SET(dbg_log)) then
         WRITE(IPT,*)'!  OBC IS OFF  '
         WRITE(IPT,*)'!'
      end if
      RETURN
   END IF
   
   
   DIM => FIND_DIM(NC_START,"nobc",FOUND)
   IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND DIMENSION 'nobc'&
         & IN THE HOTSTART FILE OBJECT")

   IOBCN_GL=DIM%DIM

   if(IOBCN_GL==0) return

   ALLOCATE(I_OBC_GL(IOBCN_GL))
   ALLOCATE(TYPE_OBC_GL(IOBCN_GL))

   VAR => FIND_VAR(NC_START,'obc_nodes',FOUND)
   IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'obc_nodes'&
        & IN THE HOTSTART FILE OBJECT")
   CALL NC_CONNECT_AVAR(VAR, I_OBC_GL)
   CALL NC_READ_VAR(VAR)
    CALL NC_DISCONNECT(VAR)

   VAR => FIND_VAR(NC_START,'obc_type',FOUND)
   IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'obc_type'&
        & IN THE HOTSTART FILE OBJECT")
   CALL NC_CONNECT_AVAR(VAR, type_obc_gl)
   CALL NC_READ_VAR(VAR)
    CALL NC_DISCONNECT(VAR)

   if(DBG_SET(dbg_log)) then
      WRITE(IPT,*)'!  FINISHED READING OBC GRID FROM HOTSTART:'
      WRITE(IPT,*)'!  OBC NODES  =           :',IOBCN_GL
      WRITE(IPT,*)'!'
   end if


 END SUBROUTINE LOAD_RESTART_OBC_GRID

 SUBROUTINE LOAD_RESTART_LSF_GRID(N_GL,I_GL,GEO_GL,WDF_GL)
   USE CONTROL
   IMPLICIT NONE
   INTEGER, INTENT(OUT)  :: N_GL
   INTEGER, INTENT(OUT), Allocatable, TARGET :: I_GL(:)
   REAL(SP), INTENT(OUT), Allocatable, TARGET :: GEO_GL(:),WDF_GL(:)

   TYPE(NCVAR),  POINTER :: VAR
   TYPE(NCDIM),  POINTER :: DIM
   LOGICAL :: FOUND


   IF(.NOT. OBC_LONGSHORE_FLOW_ON) THEN
      if(DBG_SET(dbg_log)) then
         WRITE(IPT,*)'!  OPEN BOUNDARY LONGSHORE FLOW IS OFF  '
         WRITE(IPT,*)'!'
      end if
      RETURN
   END IF

   DIM => FIND_DIM(NC_START,"nlsf",FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND DIMENSION 'nlsf'&
         & IN THE HOTSTART FILE OBJECT")

   N_GL=DIM%DIM

   ALLOCATE(I_GL(N_GL))
   ALLOCATE(GEO_GL(N_GL))
   ALLOCATE(WDF_GL(N_GL))


   VAR => FIND_VAR(NC_START,'lsf_nodes',FOUND)
   IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'lsf_nodes'&
        & IN THE HOTSTART FILE OBJECT")
   CALL NC_CONNECT_AVAR(VAR, I_GL)
   CALL NC_READ_VAR(VAR)
    CALL NC_DISCONNECT(VAR)

   VAR => FIND_VAR(NC_START,'wdf',FOUND)
   IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'wdf'&
        & IN THE HOTSTART FILE OBJECT")
   CALL NC_CONNECT_AVAR(VAR, wdf_gl)
   CALL NC_READ_VAR(VAR)
   CALL NC_DISCONNECT(VAR)

   VAR => FIND_VAR(NC_START,'geo',FOUND)
   IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'geo'&
        & IN THE HOTSTART FILE OBJECT")
   CALL NC_CONNECT_AVAR(VAR, geo_gl)
   CALL NC_READ_VAR(VAR)
   CALL NC_DISCONNECT(VAR)


   if(DBG_SET(dbg_log)) then
      WRITE(IPT,*)'!  FINISHED READING LSF GRID FROM HOTSTART:'
      WRITE(IPT,*)'!  LSF NODES  =           :',N_GL
      WRITE(IPT,*)'!'
   end if


 END SUBROUTINE LOAD_RESTART_LSF_GRID

  
  SUBROUTINE LOAD_RESTART_COORDS(X_LCL,Y_LCL)
    USE CONTROL
    IMPLICIT NONE
    REAL(SP), ALLOCATABLE,TARGET :: X_LCL(:),Y_LCL(:)

    TYPE(NCVAR),  POINTER :: VAR
    INTEGER STATUS
    LOGICAL FOUND

    ! NEED LOGIC TO DECIDE WHICH TO LOAD LAT/LON or X/Y
    SELECT CASE(GRID_FILE_UNITS)
    
    CASE('meters')
       ! CONNECT THE VARIABLE OBJECTS AND LOAD THE DATA
       VAR => FIND_VAR(NC_START,'x',FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'x'&
            & IN THE HOTSTART FILE OBJECT")
       CALL NC_CONNECT_AVAR(VAR, X_LCL)
       CALL NC_READ_VAR(VAR)
       CALL NC_DISCONNECT(VAR)
    
       
       VAR => FIND_VAR(NC_START,'y',FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'y'&
            & IN THE HOTSTART FILE OBJECT")
       CALL NC_CONNECT_AVAR(VAR, Y_LCL)
       CALL NC_READ_VAR(VAR)
       CALL NC_DISCONNECT(VAR)

    CASE('degrees')
       ! CONNECT THE VARIABLE OBJECTS AND LOAD THE DATA
       VAR => FIND_VAR(NC_START,'lon',FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'lon'&
            & IN THE HOTSTART FILE OBJECT")
       CALL NC_CONNECT_AVAR(VAR, X_LCL)
       CALL NC_READ_VAR(VAR)
       CALL NC_DISCONNECT(VAR)
       
       
       VAR => FIND_VAR(NC_START,'lat',FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'lat'&
            & IN THE HOTSTART FILE OBJECT")
       
       CALL NC_CONNECT_AVAR(VAR, Y_LCL)
       CALL NC_READ_VAR(VAR)
       CALL NC_DISCONNECT(VAR)

    CASE DEFAULT
       
       CALL FATAL_ERROR("UNKNOWN GRID_FILE_UNITS: "//TRIM(GRID_FILE_UNITS))
       
    END SELECT
    
  END SUBROUTINE LOAD_RESTART_COORDS


  SUBROUTINE LOAD_RESTART_DEPTH(H_LCL)
    IMPLICIT NONE
    REAL(SP), ALLOCATABLE, TARGET:: H_LCL(:)
    TYPE(NCVAR),  POINTER :: VAR
    INTEGER STATUS
    LOGICAL FOUND

    VAR => FIND_VAR(NC_START,'h',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'h'&
         & IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, H_LCL)
    CALL NC_READ_VAR(VAR)
    CALL NC_DISCONNECT(VAR)

  END SUBROUTINE LOAD_RESTART_DEPTH

  SUBROUTINE LOAD_RESTART_CORIOLIS(C_LCL)
    IMPLICIT NONE
    REAL(SP), ALLOCATABLE, TARGET:: C_LCL(:)
    TYPE(NCVAR),  POINTER :: VAR
    INTEGER STATUS
    LOGICAL FOUND

    VAR => FIND_VAR(NC_START,'cor',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'cor'&
         & IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, C_LCL)
    CALL NC_READ_VAR(VAR)
    CALL NC_DISCONNECT(VAR)

  END SUBROUTINE LOAD_RESTART_CORIOLIS

  SUBROUTINE LOAD_RESTART_SPONGE(SPG)
    IMPLICIT NONE
    REAL(SP), ALLOCATABLE, TARGET:: SPG(:)
    TYPE(NCVAR),  POINTER :: VAR
    INTEGER STATUS
    LOGICAL FOUND


    VAR => FIND_VAR(NC_START,'cc_sponge',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'cc_sponge'&
         & IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, SPG)
    CALL NC_READ_VAR(VAR)
    CALL NC_DISCONNECT(VAR)

  END SUBROUTINE LOAD_RESTART_SPONGE

  SUBROUTINE LOAD_RESTART_SIGMA(Z,Z1)
    USE ALL_VARS, ONLY : N2E3D
    IMPLICIT NONE
    REAL(SP), ALLOCATABLE, TARGET:: Z(:,:),Z1(:,:)
    TYPE(NCVAR),  POINTER :: VAR
    INTEGER STATUS
    LOGICAL FOUND


    VAR => FIND_VAR(NC_START,'siglev',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'siglev'&
         & IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, Z)
    CALL NC_READ_VAR(VAR)
    CALL NC_DISCONNECT(VAR)

    CALL N2E3D(Z,Z1)

  END SUBROUTINE LOAD_RESTART_SIGMA



  SUBROUTINE LOAD_COLDSTART_GRID(NVG)
    USE CONTROL
    IMPLICIT NONE
    INTEGER, ALLOCATABLE, INTENT(OUT) :: NVG(:,:)
    integer :: status,I,IERR, SENDER, nvals
    
    IF(MSR) THEN
       CALL READ_COLDSTART_GRID(GRIDUNIT,MGL,NGL,NVG)
       ! DO NOT CLOSE THE GRID FILE HERE
       ! READ THE COORDS FIRST
    END IF
    
    ! BROADCAST TO OTHER PROCS
    IF(PAR) THEN
#  if defined (MULTIPROCESSOR) 
       
       SENDER = 0 ! SEND FROM MASTER
       CALL MPI_BCAST(MGL,1,MPI_INTEGER,SENDER,MPI_COMM_WORLD,IERR)
       CALL MPI_BCAST(NGL,1,MPI_INTEGER,SENDER,MPI_COMM_WORLD,IERR)
       
       if(.not. MSR) then
          ALLOCATE(NVG(0:NGL,4),STAT=STATUS)
          IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE CONNECTIVITY")
          NVG = 0
       end if
       
       nvals= (NGL+1)*4
       CALL MPI_BCAST(NVG,nvals,MPI_INTEGER,SENDER,MPI_COMM_WORLD,IERR)
# endif
    END IF
    
  END SUBROUTINE LOAD_COLDSTART_GRID
  

  SUBROUTINE LOAD_COLDSTART_OBC_GRID(IOBCN_GL,I_OBC_GL, TYPE_OBC_GL)
    USE CONTROL
    IMPLICIT NONE
    INTEGER, INTENT(OUT)  :: IOBCN_GL
    INTEGER, INTENT(OUT), Allocatable :: I_OBC_GL(:), TYPE_OBC_GL(:)
    INTEGER :: IERR, SENDER

    IF(.NOT. OBC_ON) THEN
       if(DBG_SET(dbg_log)) then
          WRITE(IPT,*)'!  OBC IS OFF  '
          WRITE(IPT,*)'!'
       end if
       RETURN
    END IF
    
    IF(MSR) THEN
       CALL READ_COLDSTART_OBC_GRID(OBCUNIT,MGL,IOBCN_GL,I_OBC_GL, TYPE_OBC_GL)
       CLOSE(OBCUNIT)
    END IF

    IF(PAR) THEN
# if defined (MULTIPROCESSOR)
       SENDER = MSRID -1 ! SEND FROM MASTER
       CALL MPI_BCAST(IOBCN_GL,1,MPI_INTEGER,SENDER,MPI_COMM_WORLD,IERR)
       
       IF(IOBCN_GL == 0) RETURN
       
       IF(.NOT. MSR)THEN
          ALLOCATE(I_OBC_GL(IOBCN_GL))
          ALLOCATE(TYPE_OBC_GL(IOBCN_GL))
       END IF
       
       CALL MPI_BCAST(I_OBC_GL,IOBCN_GL,MPI_INTEGER,SENDER,MPI_COMM_WORLD,IERR)
       CALL MPI_BCAST(TYPE_OBC_GL,IOBCN_GL,MPI_INTEGER,SENDER,MPI_COMM_WORLD,IERR)    
# endif
    END IF
    
  END SUBROUTINE LOAD_COLDSTART_OBC_GRID
  
    SUBROUTINE LOAD_COLDSTART_LSF(N_GL,I_GL, GEO_GL,WDF_GL)
    USE CONTROL
    IMPLICIT NONE
    INTEGER, INTENT(OUT)  :: N_GL
    INTEGER, INTENT(OUT), Allocatable :: I_GL(:)
    REAL(SP), INTENT(OUT), Allocatable :: GEO_GL(:),WDF_GL(:)
    INTEGER :: IERR, SENDER

    IF(.NOT. OBC_LONGSHORE_FLOW_ON) THEN
       if(DBG_SET(dbg_log)) then
          WRITE(IPT,*)'!  OPEN BOUNDARY LONGSHORE FLOW IS OFF  '
          WRITE(IPT,*)'!'
       end if
       RETURN
    END IF
    
    IF(MSR) THEN
       CALL READ_COLDSTART_LSF(LSFUNIT,N_GL,I_GL, GEO_GL,WDF_GL)
       CLOSE(OBCUNIT)
    END IF

    IF(PAR) THEN
# if defined (MULTIPROCESSOR)
       SENDER = MSRID -1 ! SEND FROM MASTER
       CALL MPI_BCAST(N_GL,1,MPI_INTEGER,SENDER,MPI_COMM_WORLD,IERR)
       
       IF(N_GL == 0) RETURN
       
       IF(.NOT. MSR)THEN
          ALLOCATE(I_GL(N_GL))
          ALLOCATE(WDF_GL(N_GL))
          ALLOCATE(GEO_GL(N_GL))
       END IF
       
       CALL MPI_BCAST(I_GL,N_GL,MPI_INTEGER,SENDER,MPI_COMM_WORLD,IERR)
       CALL MPI_BCAST(GEO_GL,N_GL,MPI_F,SENDER,MPI_COMM_WORLD,IERR)    
       CALL MPI_BCAST(WDF_GL,N_GL,MPI_F,SENDER,MPI_COMM_WORLD,IERR)    
# endif
    END IF
    
  END SUBROUTINE LOAD_COLDSTART_LSF


  SUBROUTINE LOAD_COLDSTART_COORDS(X_GBL,Y_GBL,X_LCL,Y_LCL)
    USE CONTROL
    IMPLICIT NONE
    REAL(SP), ALLOCATABLE :: X_GBL(:),Y_GBL(:),X_LCL(:),Y_LCL(:)
    INTEGER :: SENDID

    integer status,I,IERR

    IF (MSR) THEN
       CALL READ_COLDSTART_COORDS(GRIDUNIT,MGL,X_GBL,Y_GBL)
       CLOSE(GRIDUNIT)
    END IF

    IF (SERIAL) THEN
       X_LCL = X_GBL
       Y_LCL = Y_GBL
       
    ! BROADCAST TO OTHER PROCS
    ELSE
#  if defined (MULTIPROCESSOR) 
       
       ! DEAL GLOBAL COORDINATES TO LOCAL PROCESSORS 
       SENDID=MSRID
       CALL ADEAL(MYID,SENDID,NPROCS,NXMAP,X_GBL,X_LCL)
       CALL ADEAL(MYID,SENDID,NPROCS,NXMAP,Y_GBL,Y_LCL)

# endif
    END IF

    ! DO NOT DEALLOCATE THE GLOBAL ARRAYS!

  END SUBROUTINE LOAD_COLDSTART_COORDS

  SUBROUTINE LOAD_COLDSTART_DEPTH(X_GBL,Y_GBL,H_LCL)
    USE CONTROL
    USE ALL_VARS
    IMPLICIT NONE
    integer status,I,IERR
    REAL(SP), ALLOCATABLE             :: H_LCL(:)
    REAL(SP), ALLOCATABLE, INTENT(IN) :: X_GBL(:),Y_GBL(:)
!    REAL(SP), ALLOCATABLE :: HG(:)
    INTEGER :: SENDID,SENDER

    IF (MSR) THEN
       ALLOCATE(HG(0:MGL)); HG=0.0_SP
       CALL READ_COLDSTART_DEPTH(DEPTHUNIT,MGL,X_GBL,Y_GBL,HG)
       CLOSE(DEPTHUNIT)
    END IF

# if defined (MULTIPROCESSOR) 
    IF(PAR)THEN 
      IF(.NOT.ALLOCATED(HG))  ALLOCATE(HG(0:MGL))
      SENDER = MSRID - 1
      CALL MPI_BCAST(HG,MGL,MPI_F,SENDER,MPI_COMM_WORLD,IERR)
      HG(0)=0.0_SP
    END IF
# endif

   IF (SERIAL) THEN
      H_LCL(0:MGL) = HG(0:MGL)
      
   ELSE
#  if defined (MULTIPROCESSOR)
      SENDID=MSRID
      CALL ADEAL(MYID,SENDID,NPROCS,NXMAP,HG,H_LCL)
#  endif
    END IF
    
!    IF(MSR) DEALLOCATE(HG)

  END SUBROUTINE LOAD_COLDSTART_DEPTH

  SUBROUTINE LOAD_COLDSTART_CORIOLIS(X_GBL,Y_GBL,C_LCL)
    USE CONTROL
    IMPLICIT NONE
    integer :: SENDID
    REAL(SP), ALLOCATABLE             :: C_LCL(:)
    REAL(SP), ALLOCATABLE, INTENT(IN) :: X_GBL(:),Y_GBL(:)
    REAL(SP), ALLOCATABLE :: C_GBL(:)
    
    ! WE ONLY USE A CORIOLIS FILE IF THE GRID FILE UNITS ARE METERS
    IF (MSR) THEN
       ALLOCATE(C_GBL(0:MGL)); C_GBL= 0.0_SP
       
       IF (GRID_FILE_UNITS == 'degrees') THEN
          C_GBL = Y_GBL
          
       ELSE IF (GRID_FILE_UNITS == 'meters') THEN
          
          CALL READ_COLDSTART_CORIOLIS(CORIOLISUNIT,MGL,X_GBL,Y_GBL,C_GBL)
          CLOSE (CORIOLISUNIT)
       END IF
    END IF

    IF (SERIAL) THEN
       C_LCL(0:MT) = C_GBL(0:MT)
       
    ELSE
#  if defined (MULTIPROCESSOR)
       SENDID=MSRID
       CALL ADEAL(MYID,SENDID,NPROCS,NXMAP,C_GBL,C_LCL)
#  endif
    END IF
    
    IF(MSR) DEALLOCATE(C_GBL)
    
  END SUBROUTINE LOAD_COLDSTART_CORIOLIS

  SUBROUTINE LOAD_COLDSTART_SPONGE(X_GBL,Y_GBL,NSPONGE,N_SPG,R_SPG,C_SPG,X_SPG,Y_SPG)
    USE CONTROL
    IMPLICIT NONE
    REAL(SP), ALLOCATABLE, INTENT(IN) :: X_GBL(:),Y_GBL(:)
    INTEGER, INTENT(OUT)              :: NSPONGE
    INTEGER, ALLOCATABLE, INTENT(OUT) :: N_SPG(:)
    REAL(SP), ALLOCATABLE, INTENT(OUT):: R_SPG(:),C_SPG(:),X_SPG(:),Y_SPG(:)
    INTEGER :: SENDER,IERR

    IF(MSR) THEN
       CALL READ_COLDSTART_SPONGE(SPONGEUNIT,MGL,NSPONGE,N_SPG,R_SPG,C_SPG)
       CLOSE(SPONGEUNIT)
    END IF


# if defined (MULTIPROCESSOR)
    IF(PAR)THEN
       SENDER = MSRID-1 ! SEND FROM MASTER
       CALL MPI_BCAST(NSPONGE,1,MPI_INTEGER,SENDER,MPI_FVCOM_GROUP,IERR)
    END IF
# endif
    
    IF(NSPONGE == 0) RETURN
    
    ALLOCATE(X_SPG(NSPONGE)); X_SPG = 0.0 
    ALLOCATE(Y_SPG(NSPONGE)); Y_SPG = 0.0 

    IF(MSR) THEN
       X_SPG= X_GBL(N_SPG) 
       Y_SPG= Y_GBL(N_SPG)
    END IF


# if defined (MULTIPROCESSOR)
    IF(PAR)THEN

       IF(.NOT. MSR)THEN
          ALLOCATE(R_SPG(NSPONGE)); R_SPG = 0.0 
          ALLOCATE(C_SPG(NSPONGE)); C_SPG = 0.0
          ALLOCATE(N_SPG(NSPONGE)); N_SPG = 0 
       END IF
       
       SENDER = 0 ! SEND FROM MASTER
       CALL MPI_BCAST(X_SPG,NSPONGE,MPI_F,SENDER,MPI_FVCOM_GROUP,IERR)
       CALL MPI_BCAST(Y_SPG,NSPONGE,MPI_F,SENDER,MPI_FVCOM_GROUP,IERR)

       CALL MPI_BCAST(R_SPG,NSPONGE,MPI_F,SENDER,MPI_FVCOM_GROUP,IERR)
       CALL MPI_BCAST(C_SPG,NSPONGE,MPI_F,SENDER,MPI_FVCOM_GROUP,IERR)

       CALL MPI_BCAST(N_SPG,NSPONGE,MPI_INTEGER,SENDER,MPI_FVCOM_GROUP,IERR)

    END IF
# endif
    
  END SUBROUTINE LOAD_COLDSTART_SPONGE

  SUBROUTINE LOAD_COLDSTART_SIGMA
    USE CONTROL
    IMPLICIT NONE
    INTEGER :: STYPE_LEN, IERR, SENDER, STATUS

    stype=""

    IF (MSR) THEN
       CALL READ_COLDSTART_SIGMA
       CLOSE(SIGMAUNIT)
    END IF


# if defined (MULTIPROCESSOR)
    IF(PAR)THEN
       SENDER = MSRID-1 ! SEND FROM MASTER
       CALL MPI_BCAST(KB,1,MPI_INTEGER,SENDER,MPI_FVCOM_GROUP,IERR)


       STYPE_LEN = LEN_TRIM(STYPE)
       CALL MPI_BCAST(STYPE_LEN,1,MPI_INTEGER,SENDER,MPI_FVCOM_GROUP,IERR)

       CALL MPI_BCAST(STYPE(1:STYPE_LEN),STYPE_LEN,MPI_CHARACTER,SENDER,MPI_FVCOM_GROUP,IERR)

       ! SELECT CASE BASED ON SIGMA COORDINATE TYPE
       select case(TRIM(STYPE))
          
       case(STYPE_UNIFORM)

          CALL MPI_BCAST(P_SIGMA,1,MPI_F,SENDER,MPI_FVCOM_GROUP,IERR)

       case(STYPE_GEOMETRIC)

          CALL MPI_BCAST(P_SIGMA,1,MPI_F,SENDER,MPI_FVCOM_GROUP,IERR)

       case(STYPE_TANH)

          CALL MPI_BCAST(DL2,1,MPI_F,SENDER,MPI_FVCOM_GROUP,IERR)
          CALL MPI_BCAST(DU2,1,MPI_F,SENDER,MPI_FVCOM_GROUP,IERR)

       case(STYPE_GENERALIZED)

          CALL MPI_BCAST(DLL,1,MPI_F,SENDER,MPI_FVCOM_GROUP,IERR)
          CALL MPI_BCAST(DUU,1,MPI_F,SENDER,MPI_FVCOM_GROUP,IERR)

          CALL MPI_BCAST(HMIN1,1,MPI_F,SENDER,MPI_FVCOM_GROUP,IERR)

          CALL MPI_BCAST(KU,1,MPI_INTEGER,SENDER,MPI_FVCOM_GROUP,IERR)
          CALL MPI_BCAST(KL,1,MPI_INTEGER,SENDER,MPI_FVCOM_GROUP,IERR)
          
          IF (.NOT. MSR) THEN
             IF(KU>0) THEN
                ALLOCATE(ZKU(KU),stat=status)
                if(status/=0) CALL FATAL_ERROR("LOAD_SIGMA COULD NOT ALLOCATE")
                ZKU=0.0_SP
             END IF
             IF(KL>0) THEN
                ALLOCATE(ZKL(KL),stat=status)
                if(status/=0) CALL FATAL_ERROR("LOAD_SIGMA COULD NOT ALLOCATE")
                ZKL=0.0_SP
             END IF
          END IF
             
          IF(KU>0) CALL MPI_BCAST(ZKU,KU,MPI_F,SENDER,MPI_FVCOM_GROUP,IERR)
          IF(KL>0) CALL MPI_BCAST(ZKL,KL,MPI_F,SENDER,MPI_FVCOM_GROUP,IERR)

       CASE DEFAULT

           call fatal_error('STYPE CHARACTER STRING PASSED FROM MASTER&
                & PROCESSOR DOES NOT MATCH', 'THIS IS A PROGRAM ERROR &
                &SINCE IT PASSED THE READ FOR THE MASTER')

       END select

    END IF
# endif

    KBM1 = KB -1 
    KBM2 = KB -2


  END SUBROUTINE LOAD_COLDSTART_SIGMA


  SUBROUTINE LOAD_JULIAN_OBC(NTC,NAMES,PRD,EQ_AMP,EQ_BETA,EQ_TYPE&
       &,MPTD,PHS,RF,TORG)
    USE CONTROL
    USE BCS
    IMPLICIT NONE
    INTEGER,INTENT(OUT) :: NTC  ! Number of Tidal Components
    REAL(SP),INTENT(OUT), ALLOCATABLE :: PRD(:) ! Tidal Period
    REAL(SP),INTENT(OUT), ALLOCATABLE :: EQ_AMP(:) ! Tidal Period
    REAL(SP),INTENT(OUT), ALLOCATABLE :: EQ_BETA(:) ! Tidal Period
    CHARACTER(LEN=*),ALLOCATABLE :: EQ_TYPE(:)   ! Name of Components

    REAL(SP),INTENT(OUT), ALLOCATABLE :: MPTD(:,:) ! Amplitude
    REAL(SP),INTENT(OUT), ALLOCATABLE :: PHS(:,:)! Phase
    REAL(SP),INTENT(OUT), ALLOCATABLE :: RF(:)   ! Reference Height
    CHARACTER(LEN=*), ALLOCATABLE  :: Names(:)   ! Name of Components
    TYPE(TIME),INTENT(OUT) :: TORG


    REAL(SP), ALLOCATABLE :: MPTD_GL(:,:) ! Amplitude
    REAL(SP), ALLOCATABLE :: PHS_GL(:,:)! Phase
    REAL(SP), ALLOCATABLE :: RF_GL(:)   ! Reference Height

    INTEGER :: CHAR_LEN,IERR, SENDER, STATUS, I

    IF(SERIAL) THEN
       CALL READ_JULIAN_OBC(JULOBCUNIT,NTC,NAMES,PRD,EQ_AMP,EQ_BETA&
            &,EQ_TYPE,MPTD,PHS,RF,TORG)
       CLOSE(JULOBCUNIT)
       
       IF(OBC_ELEVATION_FORCING_ON .and. SIZE(RF) /= IOBCN_GL) CALL FATAL_ERROR&
            &("LOAD_JULIAN_OBC: THE NUMBER OF  OBC NODES DOES NOT MATCH",&
            & "THE NON JULIAN ASCII FORCING FILE")
       
    END IF


# if defined (MULTIPROCESSOR)
    IF(PAR) THEN

       IF(MSR) THEN
          CALL READ_JULIAN_OBC(JULOBCUNIT,NTC,NAMES,PRD,EQ_AMP,EQ_BETA&
               &,EQ_TYPE,MPTD_GL,PHS_GL,RF_GL,TORG)
          CLOSE(JULOBCUNIT)
          
          IF(OBC_ELEVATION_FORCING_ON .and. SIZE(RF_GL) /= IOBCN_GL) CALL FATAL_ERROR&
               &("LOAD_JULIAN_OBC: THE NUMBER OF  OBC NODES DOES NOT MATCH",&
               & "THE NON JULIAN ASCII FORCING FILE")
          
       END IF

       ! SEND THE COMPONENT DATA

       SENDER = MSRID-1 ! SEND FROM MASTER
       CALL MPI_BCAST(NTC,1,MPI_INTEGER,SENDER,MPI_FVCOM_GROUP,IERR)
       
       IF(.NOT. MSR) THEN
          ALLOCATE(PRD(NTC))
          ALLOCATE(NAMES(NTC))
          
       END IF

       ! PERIOD
       CALL MPI_BCAST(PRD,NTC,MPI_F,SENDER,MPI_FVCOM_GROUP,IERR)
       
       ! NAMES
       DO I = 1, NTC
          CHAR_LEN = LEN_TRIM(NAMES(I))
          CALL MPI_BCAST(CHAR_LEN,1,MPI_INTEGER,SENDER,MPI_FVCOM_GROUP,IERR)
          
          CALL MPI_BCAST(NAMES(I)(1:CHAR_LEN),CHAR_LEN,MPI_CHARACTER,SENDER,MPI_FVCOM_GROUP,IERR)

       END DO
       
       ! TIME ORIGIN
       CALL MPI_BCAST(TORG,1,MPI_TIME,SENDER,MPI_FVCOM_GROUP,IERR)

       
# if !defined (EQUI_TIDE)
       
       IF(.NOT. MSR)THEN
          ALLOCATE(EQ_AMP(0))
          ALLOCATE(EQ_BETA(0))
          ALLOCATE(EQ_TYPE(0))
       END IF
# else
       ! SEMD THE EQUILIBRIUM TIDE DATA  
       IF(.NOT. MSR)THEN
          ALLOCATE(EQ_AMP(NTC)); EQ_AMP= 0.0_SP
          ALLOCATE(EQ_BETA(NTC)); EQ_BETA= 0.0_SP
          ALLOCATE(EQ_TYPE(NTC)); EQ_TYPE= ""
       END IF
       
       ! EQ_AMP
       CALL MPI_BCAST(EQ_AMP,NTC,MPI_F,SENDER,MPI_FVCOM_GROUP,IERR)
       
       ! EQ_BETA
       CALL MPI_BCAST(EQ_BETA,NTC,MPI_F,SENDER,MPI_FVCOM_GROUP,IERR)
       
       ! NAMES
       DO I = 1, NTC
          CHAR_LEN = LEN_TRIM(EQ_TYPE(I))
          CALL MPI_BCAST(CHAR_LEN,1,MPI_INTEGER,SENDER,MPI_FVCOM_GROUP,IERR)
          
          CALL MPI_BCAST(EQ_TYPE(I)(1:CHAR_LEN),CHAR_LEN,MPI_CHARACTER,SENDER,MPI_FVCOM_GROUP,IERR)
       END DO


# endif

       IF(.not. OBC_ELEVATION_FORCING_ON) THEN
          RETURN
       END IF
       

       ! ALLOCATE THE LOCAL ARRAY
    
       ALLOCATE(MPTD(IOBCN,NTC)); MPTD= 0.0_SP
       ALLOCATE(PHS(IOBCN,NTC)); PHS= 0.0_SP
       ALLOCATE(RF(IOBCN)); RF= 0.0_SP
    
       ! SEND THE OBC DATA      
   
       ! AMPLITUDE
       CALL ADEAL(MYID,MSRID,NPROCS,BCMAP,MPTD_GL,MPTD)

       ! PHASE
       CALL ADEAL(MYID,MSRID,NPROCS,BCMAP,PHS_GL,PHS)

       ! REFERENCE HEIGHT
       CALL ADEAL(MYID,MSRID,NPROCS,BCMAP,RF_GL,RF)


    END IF
# endif

    
  END SUBROUTINE LOAD_JULIAN_OBC



  SUBROUTINE READ_JULIAN_OBC(JULOBCUNIT_TEMP,NTC,NAMES,PRD,EQ_AMP,EQ_BETA,EQ_TYPE,MPTD,PHS,RF,TORG)
    USE CONTROL
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: JULOBCUNIT_TEMP
    INTEGER,INTENT(OUT) :: NTC  ! Number of Tidal Components
    REAL(SP),INTENT(OUT), ALLOCATABLE :: PRD(:) ! Tidal Period
    REAL(SP),INTENT(OUT), ALLOCATABLE :: EQ_AMP(:) ! Equilibrium Amplitude
    REAL(SP),INTENT(OUT), ALLOCATABLE :: EQ_BETA(:) ! Equilibrium Beta
    CHARACTER(LEN=*),INTENT(OUT),ALLOCATABLE :: EQ_TYPE(:)   ! Equilibrium Type
    REAL(SP),INTENT(OUT), ALLOCATABLE :: MPTD(:,:) ! Amplitude
    REAL(SP),INTENT(OUT), ALLOCATABLE :: PHS(:,:)! Phase
    REAL(SP),INTENT(OUT), ALLOCATABLE :: RF(:)   ! Reference Height
    CHARACTER(LEN=*),INTENT(OUT), ALLOCATABLE  :: Names(:)   ! Name of Components
    TYPE(TIME),INTENT(OUT) :: TORG

    REAL(DP) :: TORGDP

    REAL(SP), ALLOCATABLE :: RFTMP(:,:)   ! Dummy REF
    INTEGER :: ISCAN, I, MYOBC, IOS, CNT, J
    CHARACTER(LEN=80)      :: Iserr,COMPN

    CHARACTER(LEN=80)     :: line
    CHARACTER(LEN=500)     :: long_line
    CHARACTER(LEN=20), allocatable :: item(:)
    CHARACTER(LEN=80) :: TEST

    CHARACTER(LEN=80), Parameter :: line_amp = "Amplitude"
    CHARACTER(LEN=80), Parameter :: line_pha = "Phase"
    CHARACTER(LEN=80), Parameter :: line_ref = "Eref"

     LOGICAL :: ISFLOAT

    if(DBG_SET(DBG_IO)) WRITE(IPT,*) "! READING NON-JULIAN TIDAL FORCING FILE"

    ISCAN = SCAN_FILE(JULOBCUNIT_TEMP,"Tidal Component Number",ISCAL = NTC)
    IF(ISCAN /= 0) then
       write(iserr,'(I2)') ISCAN
       call fatal_error('Improper formatting of Non-Julian Tidal Forcing File: ISCAN ERROR&
            &# '//trim(iserr),&
            & 'The header must contain: "Tidal Component Number = "', &
            & 'Followed by an integer number')
    END IF
    

    if(DBG_SET(DBG_IO)) write(IPT,*) "Tidal Component Number", NTC

 
# if !defined (EQUI_TIDE)
    ALLOCATE(PRD(NTC))
    ALLOCATE(NAMES(NTC))
    ALLOCATE(EQ_AMP(0))
    ALLOCATE(EQ_BETA(0))
    ALLOCATE(EQ_TYPE(0))
# else
    ALLOCATE(PRD(NTC))
    ALLOCATE(NAMES(NTC))
    ALLOCATE(EQ_AMP(NTC))
    ALLOCATE(EQ_BETA(NTC))
    ALLOCATE(EQ_TYPE(NTC))
# endif

    DO I = 1,NTC

       !write(COMPN,'(I)') I
       write(COMPN,*) I
       compn= adjustl(compn)
       ISCAN = SCAN_FILE(JULOBCUNIT_TEMP,TRIM(COMPN),CVAL = line)
       IF(ISCAN /= 0) then
          write(iserr,'(I2)') ISCAN
          call fatal_error('Improper formatting of Non-Julian Tidal Forcing File: ISCAN ERROR&
               &# '//trim(iserr),&
               & 'The header must contain: "'//TRIM(COMPN)//' = "', &
               & 'Followed by: Name,Period,Eq. Amp.*, Eq. Beta*, Type*')
       END IF
       
       CALL SPLIT_STRING(line," ",item)

# if !defined (EQUI_TIDE)
       IF(SIZE(ITEM) >= 2) THEN
          
          ! GET THE NAME OF THE COMPONENT
          Names(I) = TRIM(item(1))
          
          ! GET THE PERIOD OF THE COMPONENT
          PRD(I) = READ_FLOAT(item(2), IOS)
          ! -------- new: Karsten Lettmann, 2016, march -------
          ! initialize Test:
          TEST = item(2)
          ! -------------- end new ----------------------------
          IF(IOS /=0) CALL FATAL_ERROR&
               &("INVALID DATA IN ASCII NON_JULIAN FORCING FILE",&
               & "Non Floating Point Value: '"//TRIM(TEST)//"' ; in line: "//TRIM(LINE))

       ELSE
          CALL FATAL_ERROR&
               & ("Improper Line in Non Julian Tidal forcing file:",&
               &  "Line :"//TRIM(LINE),&
               &  "The Tidal Component section must conatin:",&
               &  "Component# = Name, Period")
       END IF
# else
       IF(SIZE(ITEM)>=5) THEN

          ! GET THE NAME OF THE COMPONENT
          Names(I) = TRIM(item(1))
          
          ! GET THE PERIOD OF THE COMPONENT
          PRD(I) = READ_FLOAT(item(2), IOS)
          ! -------- new: Karsten Lettmann, 2016, march -------
          ! initialize Test:
          TEST = item(2)
          ! -------------- end new ----------------------------
          IF(IOS /=0) CALL FATAL_ERROR&
               &("INVALID DATA IN ASCII NON_JULIAN FORCING FILE",&
               & "Non Floating Point Value: '"//TRIM(TEST)//"' ; in line: "//TRIM(LINE))

          ! GET THE EQILIBRIUM AMPLITUDE OF THE COMPONENT
          EQ_AMP(I) = READ_FLOAT(item(3), IOS)
          ! -------- new: Karsten Lettmann, 2016, march -------
          ! initialize Test:
          TEST = item(3)
          ! -------------- end new ----------------------------
          IF(IOS /=0) CALL FATAL_ERROR&
               &("INVALID DATA IN ASCII NON_JULIAN FORCING FILE",&
               & "Non Floating Point Value: '"//TRIM(TEST)//"' ; in line: "//TRIM(LINE))


          ! GET THE EQILIBRIUM BETA OF THE COMPONENT
          ! -------- new: Karsten Lettmann, 2016, march -------
          !EQ_BETA(I) = READ_FLOAT(item(3), IOS)  ! original line
          EQ_BETA(I) = READ_FLOAT(item(4), IOS)
          ! initialize Test:
          TEST = item(4)
          ! -------------- end new ----------------------------
          IF(IOS /=0) CALL FATAL_ERROR&
               &("INVALID DATA IN ASCII NON_JULIAN FORCING FILE",&
               & "Non Floating Point Value: '"//TRIM(TEST)//"' ; in line: "//TRIM(LINE))

          ! GET THE EQILIBRIUM TYPE OF THE COMPONENT
          EQ_TYPE(I) = TRIM(item(5))

       ELSE

          CALL FATAL_ERROR&
               & ("Improper Line in Non Julian Tidal forcing file:",&
               &  "Line :"//TRIM(LINE),&
               &  "The Tidal Component section must conatin:",&
               &  "Component# = Name, Period, Equilibrium Amp, Eq. Beta, Eq. Type")
          
       END IF
# endif

       DEALLOCATE(ITEM)

    END DO

    ! GET THE TIME ORIGIN
    ISCAN = SCAN_FILE(JULOBCUNIT_TEMP,"Time Origin",CVAL = LINE)
    IF(ISCAN /= 0) then
       write(iserr,'(I2)') ISCAN
       call fatal_error('Improper formatting of Non-Julian Tidal Forcing File: ISCAN ERROR&
            &# '//trim(iserr),&
            & 'The header must contain: "Time Origin = "', &
            & 'Followed by a data or time')
    END IF
    
    IF (USE_REAL_WORLD_TIME) THEN

       TORG = READ_DATETIME(LINE,DATE_FORMAT,TIMEZONE,IOS)
       IF(IOS == 0) Call Fatal_Error&
            &("NON JULIAN TIDAL FORCING FILE - Time Origin error"&
            &"Could not read the date string Time Origin",&
            & "The model is running using real dates")

       if(DBG_SET(DBG_IO)) CALL PRINT_REAL_TIME(TORG,IPT,"NON JULIAN T0")

    ELSE

       CALL SPLIT_STRING(line," ",item)

       IF (size(item)>1 ) THEN
          TORGDP = READ_FLOAT(ITEM(1),IOS)
          IF(IOS /=0) CALL FATAL_ERROR&
               & ("NON JULIAN TIDAL FORCING FILE - Time Origin error",&
               &  "Could not read the floating point number",&
               &  "Line :"//TRIM(LINE),&
               &  "The model is running using ideal time (starting from 0.0)")
       ELSE
          CALL FATAL_ERROR&
            &("NON JULIAN TIDAL FORCING FILE - Time Origin error",&
            & "Could not read the string Time Origin",&
            & "Line :"//TRIM(LINE),&
            & "The model is running using ideal time (starting from 0.0)")
       END IF

       IF(ITEM(2)=="seconds") THEN
            TORG = seconds2time(TORGDP)
       ELSEIF(ITEM(2)=="days" .and. size(item)>1 ) THEN
            TORG = days2time(TORGDP)
       ELSE
          CALL FATAL_ERROR&
            &("NON JULIAN TIDAL FORCING FILE - Time Origin error",&
            & "Could not read the string Time Origin",&
            & "The model is running using ideal time (starting from 0.0)")
       END IF

       if(DBG_SET(DBG_IO)) CALL PRINT_TIME(TORG,IPT,"NON JULIAN T0")


    END IF

    ISCAN = SCAN_FILE(JULOBCUNIT_TEMP,"OBC Node Number",ISCAL = MYOBC)
    IF(ISCAN /= 0) then
       write(iserr,'(I2)') ISCAN
       call fatal_error&
            &('Improper formatting of Non-Julian Tidal Forcing File: ISCAN ERROR# '//trim(iserr),&
            & 'The header must contain: "OBC Node Number = "', &
            & 'Followed by and integer number of boundary nodes')
    END IF
    
    if(DBG_SET(DBG_IO)) write(IPT,*) "OBC NODE NUMBER =",MYOBC
    
    ALLOCATE(MPTD(MYOBC,NTC))
    ALLOCATE(PHS(MYOBC,NTC))
    ALLOCATE(RFTMP(MYOBC,1))
    ALLOCATE(RF(MYOBC))

    IF(MYOBC == 0) RETURN



    ! READ THE AMPLITUDE

    if(DBG_SET(DBG_IO)) write(IPT,*) "READING AMPLITUDE DATA"
    rewind JULOBCUNIT_TEMP

    DO WHILE(.TRUE.)
       READ(JULOBCUNIT_TEMP,*,IOSTAT=IOS) line
       if (IOS /= 0) CALL FATAL_ERROR&
            &("Could not read Non Julian Tidal Forcing file. no keyword: 'Amplitude'")
       
       IF(Line == Line_amp) Exit
       
    END DO

    CNT = 0
    DO
       READ(JULOBCUNIT_TEMP,'(a)',IOSTAT=IOS) Long_LINE
       IF(IOS /=0) CALL FATAL_ERROR&
            &("While Reading Non Julian Tidal forcing Amplitude:",&
            & "Invalid line or end of file reached with out end of section!")


       IF(Long_Line == Line_amp) THEN

          ! IF FINISHED READING SECTION => EXIT
          IF(CNT == MYOBC) EXIT

          ! OTHERWISE CALL AN ERROR
          CALL FATAL_ERROR&
            &("Unexpected end of section Amplitude in Non Julian Tidal forcing file",&
            & "Check the number of nodes in the list")
       END IF

       Call Parse_tide(Long_line,CNT,NTC,MPTD,IOS)
       ! IF THIS IS A COMMENT LINE OR BLANK
       IF(IOS /=0) CYCLE

       ! ELSE CONTINUE TO NEXT
    END DO


    ! READ THE PHASE
    if(DBG_SET(DBG_IO)) write(IPT,*) "READING PHASE DATA"
    rewind JULOBCUNIT_TEMP

    DO WHILE(.TRUE.)
       READ(JULOBCUNIT_TEMP,*,IOSTAT=IOS) line
       if (IOS /= 0) CALL FATAL_ERROR&
            &("Could not read Non Julian Tidal Forcing file. no keyword: 'Phase'")
       
       IF(Line == Line_PHA) Exit
       
    END DO

    CNT = 0
    DO WHILE(CNT <= MYOBC)
       READ(JULOBCUNIT_TEMP,'(a)',IOSTAT=IOS) Long_LINE
       IF(IOS /=0) CALL FATAL_ERROR&
            &("While Reading Non Julian Tidal forcing Phase:",&
            & "Invalid line or end of file reached with out end of section!")


       IF(Long_Line == Line_pha) THEN

          ! IF FINISHED READING SECTION => EXIT
          IF(CNT == MYOBC) EXIT

          ! OTHERWISE CALL AN ERROR
          CALL FATAL_ERROR&
            &("Unexpected end of section Phase in Non Julian Tidal forcing file",&
            & "Check the number of nodes in the list")
       END IF


       Call Parse_tide(Long_line,CNT,NTC,PHS,IOS)
       ! IF THIS IS A COMMENT LINE OR BLANK
       IF(IOS /=0) CYCLE

    END DO

    ! READ THE REFERENCE HEIGHT

    if(DBG_SET(DBG_IO)) write(IPT,*) "READING REFERENCE HEIGHT DATA"
    rewind JULOBCUNIT_TEMP

    DO WHILE(.TRUE.)
       READ(JULOBCUNIT_TEMP,*,IOSTAT=IOS) line
       if (IOS /= 0) CALL FATAL_ERROR&
            &("Could not read Non Julian Tidal Forcing file. no keyword: 'Eref'")
       
       IF(Line == Line_ref) Exit
       
    END DO

    CNT = 0
    DO WHILE(CNT <= MYOBC)
       READ(JULOBCUNIT_TEMP,'(a)',IOSTAT=IOS) Long_LINE
       IF(IOS /=0) CALL FATAL_ERROR&
            &("While Reading Non Julian Tidal forcing Eref:",&
            & "Invalid line or end of file reached with out end of section!")

       IF(Long_Line == Line_ref) THEN

          ! IF FINISHED READING SECTION => EXIT
          IF(CNT == MYOBC) EXIT

          ! OTHERWISE CALL AN ERROR
          CALL FATAL_ERROR&
            &("Unexpected end of section Eref in Non Julian Tidal forcing file",&
            & "Check the number of nodes in the list")
       END IF



       Call Parse_tide(Long_line,CNT,1,RFTMP,IOS)
       ! IF THIS IS A COMMENT LINE OR BLANK
       IF(IOS /=0) CYCLE


    END DO

    RF = RFTMP(:,1)
    DEALLOCATE(RFTMP)


    if(DBG_SET(DBG_IO)) WRITE(IPT,*) "! FINISHED READING NON-JULIAN TIDAL FORCING FILE"

  END SUBROUTINE READ_JULIAN_OBC


  SUBROUTINE Parse_tide(line,cnt,ntc,data,ierr)
    implicit none
    
    CHARACTER(LEN=*) :: line
    Integer :: cnt, ntc, ierr
    Real(sp), allocatable :: data(:,:) 
    
    CHARACTER(LEN=20), allocatable :: item(:) !Automatically deallocate on exit!
    CHARACTER(LEN=80) :: TEST
    INTEGER :: I, VAL, J
    LOGICAL :: ISFLOAT

    ! If the string is empty
    ierr = -1
    IF(LEN_TRIM(line)<=1) return

    CALL SPLIT_STRING(line," ",item)

    VAL = READ_INT(item(1),IERR)
    IF(IERR /= 0) RETURN
       

    CNT = CNT + 1

    IF(CNT >SIZE(DATA,1)) CALL FATAL_ERROR&
         &("THERE IS A MISTAKE IN THE NON JULIAN TIDAL FORCING INPUT FILE",&
         & "THERE ARE MORE BOUNDARY POINTS LISTED THAN THE STATED NUMBER?",&
         & "Line :"//TRIM(LINE))

    IF(VAL /= CNT) CALL FATAL_ERROR&
         &("THERE IS A MISTAKE IN THE NON JULIAN TIDAL FORCING INPUT FILE",&
         & "THE LIST OF BOUNDARY POINTS IS OUT OF ORDER OR CAN NOT BE READ?",&
         & "Line :"//TRIM(LINE))


    ! NOW READ THE DATA
    DO J = 1, NTC

       IF(J+1 > SIZE(ITEM)) CALL FATAL_ERROR&
            &("INVALID LINE IN NON JULIAN TIDAL FORCING FILE",&
            & "Incorrect number of tidal compontents",&
            & "Line: "//TRIM(line))

       DATA(CNT,J) = READ_FLOAT(item(J+1),ierr)
       if(IERR/=0)CALL FATAL_ERROR&
            &("INVALID DATA IN ASCII NON_JULIAN FORCING FILE",&
            & "Non Floating Point Value: '"//item(J+1)//"' ; in li&
            &ne: "//TRIM(LINE))

    END DO
    
    
    ! SUCCESS
    IERR = 0
    
  END SUBROUTINE Parse_tide


  SUBROUTINE READ_COLDSTART_GRID(GRIDUNIT,MGL,NGL,NVG)
    IMPLICIT NONE
    INTEGER, INTENT(OUT) :: NGL, MGL
    INTEGER, INTENT(IN)    :: GRIDUNIT
    INTEGER, ALLOCATABLE   :: NVG(:,:)
    CHARACTER(LEN=80)      :: temp
    INTEGER                :: I,LM1,J,IOS,LMAX,CellCount,NodeCount
    INTEGER                :: N1, N2, N3, ISCAN,K
    INTEGER :: SENDER,nvals, IERR, STATUS
    real(SP)               :: X1, Y1
    !----------------Determine Number of Nodes -------------------------------!
    ISCAN = SCAN_FILE(GRIDUNIT,"Node Number",ISCAL = MGL)
    IF(ISCAN /= 0) then
       write(temp,'(I2)') ISCAN
       call fatal_error('Improper formatting of GRID FILE: ISCAN ERROR&
            &# '//trim(temp),&
            & 'The header must contain: "Node Number = "', &
            & 'Followed by an integer number of Nodes')
    END IF
    !----------------Determine Number of Elements -------------------------------!
    
    ISCAN = SCAN_FILE(GRIDUNIT,"Cell Number",ISCAL = NGL)
    IF(ISCAN /= 0)then
       write(temp,'(I2)') ISCAN
       call fatal_error('Improper formatting of GRID FILE: ISCAN ERROR&
            &#'//trim(temp),&
            & 'The header must contain: "Cell Number = "', &
            & 'Followed by an integer number of Cells.')
       
    END IF
!----------------Read Connectivity Array----------------------------------!
!

    ! FIND FIRST LINE of CONNECTIVITY ARRAY
    rewind GRIDUNIT
    DO WHILE(.TRUE.)
       READ(GRIDUNIT,*,IOSTAT=IOS,END=99)J,N1,N2,N3
       if (IOS == 0) then
          BackSpace GRIDUNIT
          exit
       end if

       CYCLE
       
99 Call FATAL_ERROR('Improper formatting of GRID FILE:',&
           &'Reached end of file with out finding CONNECTIVITY data?',&
           &'FORMAT: CELL# NODE# NODE# NODE# (ALL INTEGERS)')
          
    END DO

    ! READ IN CONNECTIVITY

    ALLOCATE(NVG(0:NGL,4)); NVG=0
    J = 0
    I = 1
    LM1=1
    DO WHILE(.TRUE.)

       READ(GRIDUNIT,*,IOSTAT=IOS)J,N1,N2,N3
       IF(IOS < 0) CALL FATAL_ERROR('ERROR READING GRID FILE CONNECTIV&
            &ITY LIST')
       

       IF(J == 1 .AND. LM1 /= 1)THEN
          CellCount = LM1
          BACKSPACE GRIDUNIT

          READ(GRIDUNIT,*) J
          IF(J .GT. 0) THEN
             DO K=1,J
                BACKSPACE GRIDUNIT
             END DO
             EXIT
          ELSE
             CALL WARNING('Trouble reading grid file!')
             
             EXIT
          END IF
       END IF

       
       IF(I > NGL) CALL FATAL_ERROR &
            &('Number of rows of data in the grid file CONNECTIVITY data exceeds the stated number of Cells ?')

       ! LIST IS REORDERD!
       NVG(I,1)=N1
       NVG(I,2)=N3
       NVG(I,3)=N2
       NVG(I,4)=N1
       

       I = I + 1
       LM1 = J
    END DO

    if ( CELLCOUNT .NE. NGL) CALL FATAL_ERROR&
         ('The number of rows of data in the grid file CONNECTIVITY does  not equal the stated number?')

    NODECOUNT = MAX(MAXVAL(NVG(:,1)), MAXVAL(NVG(:,2)), MAXVAL(NVG(:,3)))
    if ( NODECOUNT .NE. MGL) &
         & CALL FATAL_ERROR('The number of nodes in the grid file CONNECTIVITY does not equal the stated number ?')

    
   if(DBG_SET(dbg_log)) then
      WRITE(IPT,*)'!  Finished Reading Grid File'
      WRITE(IPT,*)'!  # OF NODES            :',MGL
      WRITE(IPT,*)'!  # OF CELLS            :',NGL
      WRITE(IPT,*)'!'
   end if

  END SUBROUTINE READ_COLDSTART_GRID


 SUBROUTINE READ_COLDSTART_OBC_GRID(OBCUNIT,MGL,IOBCN_GL,I_OBC_GL, TYPE_OBC_GL)
    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: OBCUNIT, MGL
    INTEGER, INTENT(OUT)  :: IOBCN_GL
    INTEGER, INTENT(OUT), Allocatable :: I_OBC_GL(:), TYPE_OBC_GL(:)
    CHARACTER(LEN=80)    :: temp,temp2
    INTEGER              :: I,J,IOS,NodeCount,ISCAN
    INTEGER              :: N1, N2, N3
    !==============================================================================|

    !
    !----------------Determine Number of Nodes -------------------------------!

    ISCAN = SCAN_FILE(OBCUNIT,"OBC Node Number",ISCAL = IOBCN_GL)
    IF(ISCAN /= 0) then
       write(temp,'(I2)') ISCAN
       call fatal_error('Improper formatting of OBC FILE: ISCAN ERROR&
            &# '//trim(temp),&
            & 'The header must contain: "OBC Node Number ="', &
            & 'Followed by an integer number of boundary nodes')
    END IF


    if(IOBCN_GL==0) then

       if(DBG_SET(dbg_log)) then
          WRITE(IPT,*)'!  Finished Reading OBC File: NO OPEN BOUNDARY'
          WRITE(IPT,*)'!  OBC NODES  =           :',IOBCN_GL
          WRITE(IPT,*)'!'
       end if
       return
    end if

!----------------Read OBC Array----------------------------------!
!
   ! FIND FIRST LINE of )BC ARRAY
    rewind OBCUNIT
    DO WHILE(.TRUE.)
       READ(OBCUNIT,*,IOSTAT=IOS,END=99)N1,N2,N3
       if (IOS == 0) then
          BackSpace OBCUNIT
          exit
       end if

       CYCLE
       
99 Call FATAL_ERROR('Improper formatting of OBC FILE:',&
           &'Reached end of file with out finding OBC data?',&
           &'FORMAT: OBCNODE# GLNODE# TYPE# (ALL INTEGERS)')
          
    END DO
    ALLOCATE(I_OBC_GL(IOBCN_GL))
    ALLOCATE(TYPE_OBC_GL(IOBCN_GL))

    
    I = 0
    DO WHILE(.TRUE.)

       READ(OBCUNIT,*,IOSTAT=IOS) N1,N2,N3
       IF(IOS < 0) exit
       
       I = I + 1
       IF(I > IOBCN_GL) CALL FATAL_ERROR('Number of rows of data in the OBC file &
            &exceeds the stated number of boundary nodes in the header ?')

       IF( 1 > N2 .or. N2 > MGL) then
          write(temp,'(I8)') I
          write(temp2,'(I8)') MGL
          CALL FATAL_ERROR('OPEN BOUNDARY NODE NUMBER'//trim(temp)//&
               & 'IS NOT IN THE GLOBAL DOMAIN',&
               & 'CHECK INPUT FILE AND ENSURE OPEN BOUNDARY NODES <= '//trim(temp2))
       END IF

       IF( 1 > N3 .or. N3 > 10) then
          write(temp,'(I8)') I
          CALL FATAL_ERROR('OPEN BOUNDARY NODE NUMBER'//trim(temp)//&
               & ' IS NOT IN THE VALID RANGE',&
               & 'THE OPEN BOUNDARY NODE TYPE MUST BE GREATER THAN 0',&
               & 'AND LESS THAN 11. SEE MOD_OBC.F FOR DESCRIPTION')
       END IF


       I_OBC_GL(I)    = N2
       TYPE_OBC_GL(I) = N3

    END DO

    if ( I .NE. IOBCN_GL) &
         & CALL FATAL_ERROR('The number of rows of data in the OBC file does&
         & not equal the number of nodes in the header?')
   

   if(DBG_SET(dbg_log)) then
      WRITE(IPT,*)'!  Finished Reading OBC File'
      WRITE(IPT,*)'!  OBC NODES  =           :',IOBCN_GL
      WRITE(IPT,*)'!'
   end if
   
    
 END SUBROUTINE READ_COLDSTART_OBC_GRID

 SUBROUTINE READ_COLDSTART_LSF(LSFUNIT,N_GL,I_GL,GEO_GL,WDF_GL)
    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: LSFUNIT
    INTEGER, INTENT(OUT)  :: N_GL
    INTEGER, INTENT(OUT), Allocatable :: I_GL(:)
    REAL(SP), INTENT(OUT), Allocatable :: GEO_GL(:),WDF_GL(:)
    CHARACTER(LEN=80)    :: temp,temp2
    INTEGER              :: I,J,IOS,NodeCount,ISCAN
    INTEGER              :: N1, N2, N3, N4
    REAL(SP)             :: R1, R2
    !==============================================================================|

    !
    !----------------Determine Number of Nodes -------------------------------!

    ISCAN = SCAN_FILE(LSFUNIT,"Longshore Flow Node Number",ISCAL = N_GL)
    IF(ISCAN /= 0) then
       write(temp,'(I2)') ISCAN
       call fatal_error('Improper formatting of LONGSHORE FLOW FILE: ISCAN ERROR&
            &# '//trim(temp),&
            & 'The header must contain: "Longshre Flow Node Number ="', &
            & 'Followed by an integer number of nodes')
    END IF


    if(N_GL==0) then

       if(DBG_SET(dbg_log)) then
          WRITE(IPT,*)'!  Finished Reading LSF file: No Long Shore Flow Nodes!'
          WRITE(IPT,*)'!  LSF NODES  =           :',N_GL
          WRITE(IPT,*)'!'
       end if
       return
    end if

!----------------Read OBC Array----------------------------------!
!
   ! FIND FIRST LINE of )BC ARRAY
    rewind LSFUNIT
    DO WHILE(.TRUE.)
       READ(LSFUNIT,*,IOSTAT=IOS,END=99)N1,N2,R1,R2
       if (IOS == 0) then
          BackSpace LSFUNIT
          exit
       end if

       CYCLE
       
99 Call FATAL_ERROR('Improper formatting of LongShore Flow FILE:',&
           &'Reached end of file with out finding LSF data?',&
           &'FORMAT: LSFNODE# GLNODE# GEO WND (last two are real 0<=X<=1)')
          
    END DO
    ALLOCATE(I_GL(N_GL))
    ALLOCATE(GEO_GL(N_GL))
    ALLOCATE(WDF_GL(N_GL))

    
    I = 0
    DO WHILE(.TRUE.)

       READ(LSFUNIT,*,IOSTAT=IOS) N1,N2,R1,R2
       IF(IOS < 0) exit
       
       I = I + 1
       IF(I > N_GL) CALL FATAL_ERROR('Number of rows of data in the LongShore Flow file &
            &exceeds the stated number of boundary nodes in the header ?')


       I_GL(I)    = N2
       GEO_GL(I) = R1
       WDF_GL(I) = R2

    END DO

    if ( I .NE. N_GL) &
         & CALL FATAL_ERROR('The number of rows of data in the LONGSHORE FLOW file does&
         & not equal the number of nodes in the header?')
   

   if(DBG_SET(dbg_log)) then
      WRITE(IPT,*)'!  Finished Reading LSF File'
      WRITE(IPT,*)'!  LSF NODES  =           :',N_GL
      WRITE(IPT,*)'!'
   end if
   
    
 END SUBROUTINE READ_COLDSTART_LSF


  SUBROUTINE READ_COLDSTART_COORDS(GRIDUNIT,MGL,XG2,YG2)
    USE ALL_VARS, ONLY : NVG
    IMPLICIT NONE
    INTEGER, INTENT(IN)    :: MGL
    INTEGER, INTENT(IN)    :: GRIDUNIT
    REAL(SP), ALLOCATABLE  :: XG2(:),YG2(:)
    CHARACTER(LEN=80)      :: temp
    INTEGER                :: I,LM1,J,IOS,LMAX,CellCount,NodeCount
    INTEGER                :: N1, N2, N3, ISCAN
    real(SP)               :: X1, Y1
#  if defined (SPHERICAL)
   INTEGER K,ITMP
   REAL(DP) VX1,VY1,VX2,VY2,VX3,VY3,EVX12,EVX13,EVX23,&
            EVY12,EVY13,EVY23,EVXY,VX12,VY12,VX23,VY23,VX31,VY31,&
	    X1_DP,Y1_DP,X2_DP,Y2_DP,DTMP_DP
#  endif
    REAL(SP)               :: X_CW(3),Y_CW(3)
    INTEGER                :: IS_CW
    !==============================================================================|
    
    ! DO NOT REWIND - SAVED LOCATION FROM READING CONNECTIVITY
    I =0
    DO WHILE (.TRUE.)
       
       READ(GRIDUNIT,*,IOSTAT=IOS)J,X1,Y1
       IF(IOS<0) exit
       
       I = I + 1
       IF(I > MGL) THEN
          write(ipt,*) "Read ", I, "; lines of coordiante data with out reaching EOF?"
          CALL FATAL_ERROR('Number of rows of data in the grid file coordinates exceeds the stated number of nodes ?')
       END IF

#  if defined (SPHERICAL)
       IF(X1 < 0.0) X1 = X1 + 360.0_sp
#  endif
       
       XG2(I) = X1
       YG2(I) = Y1
    END DO
    
    if( I .NE. MGL) THEN
       write(ipt,*) "Read, ", I, "rows of data but mgl= ",mgl 
       CALL FATAL_ERROR('Number of rows of data in the grid file coordinates does not equal the stated number of nodes ?')
    END if
    
    IF(MSR)THEN
      X_CW(1:3)=XG2(NVG(1,1:3))
      Y_CW(1:3)=YG2(NVG(1,1:3))
      CALL IS_TRI_CW(X_CW,Y_CW,IS_CW)
      IF(IS_CW /= 1) CALL FATAL_ERROR('Three nodes on each triangular cell in file *_grd.dat should be identified', &
                                      & 'counter-clockwise...')
    END IF

      if(DBG_SET(dbg_log)) then
      WRITE(IPT,*)'!  Finished Reading coordinates from Grid File'
      WRITE(IPT,*)'!  Max/Min(X)  =          :',maxval(XG2(1:MGL)),minval(XG2(1:MGL))
      WRITE(IPT,*)'!  Max/Min(Y)  =          :',maxval(YG2(1:MGL)),minval(YG2(1:MGL))
      WRITE(IPT,*)'!'
   end if
    
  END SUBROUTINE READ_COLDSTART_COORDS
  
  
  SUBROUTINE READ_COLDSTART_DEPTH(DEPTHUNIT,MGL,XG2,YG2,HG2)
!    USE ALL_VARS !Jadon Ge
    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: DEPTHUNIT, MGL
    REAL(SP),ALLOCATABLE,INTENT(IN) :: XG2(:),YG2(:)
    REAL(SP),ALLOCATABLE ::HG2(:)
    CHARACTER(LEN=80)    :: temp,XCHR,YCHR
    INTEGER              :: I,J,IOS,NodeCount, ISCAN
    Real(sp)             :: X1, Y1, HDEP,DiffxMax,DiffyMax
    logical              :: back, test
    !==============================================================================|
    test=.false.
    DIFFXMAX=0.0_SP
    DIFFYMAX=0.0_SP
    !
    !----------------Determine Number of Nodes -------------------------------!
    
    ISCAN = SCAN_FILE(DEPTHUNIT,"Node Number",ISCAL = NODECOUNT)
    IF(ISCAN /= 0) then
       write(temp,'(I2)') ISCAN
       call fatal_error('Improper formatting of DEPTH FILE: ISCAN ERROR&
            &# '//trim(temp),&
            & 'The header must contain: "Node Number ="', &
            & 'Followed by an integer number of nodes')
    END IF

    if ( NODECOUNT .NE. MGL) &
         & CALL FATAL_ERROR('The stated number of nodes in the depth file',&
         ' does not match the number in the grid file')



!----------------Read Depth Array----------------------------------!
!

   ! FIND FIRST LINE of CONNECTIVITY ARRAY
    rewind DEPTHUNIT
    DO WHILE(.TRUE.)
       READ(DEPTHUNIT,*,IOSTAT=IOS,END=99)X1,Y1,HDEP
       if (IOS == 0) then
          BackSpace DEPTHUNIT
          exit
       end if

       CYCLE
       
99 Call FATAL_ERROR('Improper formatting of DEPTH FILE:',&
           &'Reached end of file with out finding DEPTH data?',&
           &'FORMAT: X Y H (ALL REALS)')
          
    END DO

    I = 0
    DO WHILE(.TRUE.)

       READ(DEPTHUNIT,*,IOSTAT=IOS) X1,Y1,HDEP
       IF(IOS < 0) exit
       
       I = I + 1
       IF(I > MGL) CALL FATAL_ERROR('Number of rows of data in the depth file &
            &exceeds the number of nodes ?')

#  if defined (SPHERICAL)
       IF(X1 < 0.0) X1 = X1 + 360.0_SP
#  endif

       ! THIS SHOULD SCREEN OUT MOST ROUNDOFF ERRORS
       IF (XG2(I) .NE. X1 .or. YG2(I) .NE. Y1) then
          test = .true.
          
          DIFFXMAX=max(diffxmax,abs(XG2(I)-X1))
          DIFFYMAX=max(diffymax,abs(YG2(I)-Y1))
          
!!$  UNCOMMNET FOR STRICT MATCHING
!!$             write(temp,'(I8)') I
!!$             Call FATAL_ERROR('Grid Coordinates do not match between the&
!!$                  & grid file and the depth file','The bad value occurs&
!!$                  & at Node Number: '//trim(temp))

       END IF

       HG2(I)=HDEP
    END DO

    WRITE(XCHR,*)DIFFXMAX
    WRITE(YCHR,*)DIFFYMAX

    IF(TEST) CALL WARNING("THE GRID FILE AND DEPTH FILE COORDINATES DO NOT MATCH EXACTLY",&
         & "LARGEST DIFFERENCE IN X-COORDINATE:"//TRIM(XCHR),&
         & "LARGEST DIFFERENCE IN Y-COORDINATE:"//TRIM(YCHR),&
         & "See mod_input.F::READ_COLDSTART_DEPTH for details")

    if ( I .NE. MGL) &
         & CALL FATAL_ERROR('The number of rows of data in the depth file does&
         & not equal the number of nodes?')

   if(DBG_SET(dbg_log)) then
      WRITE(IPT,*)'!  Finished Reading DEPTH File'
      WRITE(IPT,*)'!  Max DEPTH =           :',maxval(HG2(1:MGL))
      WRITE(IPT,*)'!  Min DEPTH =           :',minval(HG2(1:MGL))
      WRITE(IPT,*)'!'
   end if

!%# if defined (DATA_ASSIM_OI)
! Jadon Ge added for the HG support in mod_assim_io.F
!   ALLOCATE(HG(0:MGL))  ; HG = 0.0_SP
!%   HG(1:MGL) = HG2(1:MGL)
!%# endif

 END SUBROUTINE READ_COLDSTART_DEPTH
!==============================================================================|
  SUBROUTINE READ_COLDSTART_CORIOLIS(CORIOLISUNIT,MGL,XG,YG,CORG)
    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: CORIOLISUNIT, MGL
    REAL(SP),ALLOCATABLE,INTENT(IN) :: XG(:),YG(:)
    REAL(SP),ALLOCATABLE ::CORG(:)
    CHARACTER(LEN=80)    :: INPLINE,temp,xchr,ychr
    INTEGER              :: I,J,IOS,NodeCount,ISCAN
    Real(sp)             :: X1, Y1, C1,DiffxMax,DiffyMax
    logical              :: test
    !==============================================================================|
    test=.false.  
    DIFFXMAX=0.0_SP
    DIFFYMAX=0.0_SP
    !
    !----------------Determine Number of Nodes -------------------------------!
    
    ISCAN = SCAN_FILE(CORIOLISUNIT,"Node Number",ISCAL = NODECOUNT)
    IF(ISCAN /= 0) then
       write(temp,'(I2)') ISCAN
       call fatal_error('Improper formatting of CORIOLIS FILE: ISCAN ERROR&
            &# '//trim(temp),&
            & 'The header must contain: "Node Number = "', &
            & 'Followed by an integer number of nodes')
    END IF

    if ( NODECOUNT .NE. MGL) &
         & CALL FATAL_ERROR('The stated number of nodes in the coriolis file',&
         ' does not match the number in the gird file')

!----------------Read Depth Array----------------------------------!
!

    ! FIND FIRST LINE of CORIOLIS DATA
    rewind CORIOLISUNIT
    DO WHILE(.TRUE.)
       READ(CORIOLISUNIT,*,IOSTAT=IOS,END=99)X1,Y1,C1
       if (IOS == 0) then
          BackSpace CORIOLISUNIT
          exit
       end if

       CYCLE
99 Call FATAL_ERROR('Improper formatting of CORIOLIS FILE:',&
           &'Reached end of file with out finding CORIOLIS data?',&
           &'FORMAT: X Y COR (ALL REALS)')
    END DO

    ! READ IN CORIOLIS DATA

    I = 0
    DO WHILE(.TRUE.)

       READ(CORIOLISUNIT,*,IOSTAT=IOS) X1,Y1,C1
       IF(IOS < 0) exit
       
       I = I + 1
       IF(I > MGL) CALL FATAL_ERROR('Number of rows of data in the Coriolis file &
            &exceeds the number of nodes ?')

#  if defined (SPHERICAL)
       IF(X1 < 0.0) X1 = X1 + 360.0_sp
#  endif

       ! THIS SHOULD SCREEN OUT MOST ROUNDOFF ERRORS
       IF (XG(I) .NE. X1 .or. YG(I) .NE. Y1) then

          test=.true.
          DIFFXMAX=max(diffxmax,abs(XG(I)-X1))
          DIFFYMAX=max(diffymax,abs(YG(I)-Y1))

!!$  UNCOMMNET FOR STRICT MATCHING 
!!$             write(temp,'(I8)') I
!!$             Call FATAL_ERROR('Grid Coordinates do not match between the&
!!$                  & grid file and the coriolis file','The bad value occurs&
!!$                  & at Node Number: '//trim(temp))


       END IF

       CORG(I)=C1
    END DO

    WRITE(XCHR,*)DIFFXMAX
    WRITE(YCHR,*)DIFFYMAX
    IF(TEST) CALL WARNING("THE GRID FILE AND CORIOLIS FILE COORDINATES DO NOT MATCH EXACTLY",&
         & "LARGEST DIFFERENCE IN X-COORDINATE:"//TRIM(XCHR),&
         & "LARGEST DIFFERENCE IN Y-COORDINATE:"//TRIM(YCHR),&
         & "See mod_input.F::READ_COLDSTART_CORIOLIS for details")

    if ( I .NE. MGL) &
         & CALL FATAL_ERROR('The number of rows of data in the coriolis file does&
         & not equal the number of nodes?')

   if(DBG_SET(dbg_log)) then
      WRITE(IPT,*)'!  Finished Reading Coriolis File'
      WRITE(IPT,*)'!  Max Coriolis =           :',maxval(CORG(1:MGL))
      WRITE(IPT,*)'!  Min Coriolis =           :',minval(CORG(1:MGL))
      WRITE(IPT,*)'!'
   end if

    
 END SUBROUTINE READ_COLDSTART_CORIOLIS
 


 SUBROUTINE READ_COLDSTART_SPONGE(SPONGEUNIT,MGL,NSPONGE, N_SPG,R_SPG,C_SPG)
    IMPLICIT NONE
    REAL(SP), INTENT(OUT), ALLOCATABLE :: R_SPG(:),C_SPG(:) 
    INTEGER, INTENT(OUT), ALLOCATABLE  :: N_SPG(:)
    INTEGER, INTENT(IN)   :: SPONGEUNIT,MGL
    INTEGER, INTENT(OUT)  :: NSPONGE
    CHARACTER(LEN=80)     :: temp,temp2
    INTEGER               :: I,J,IOS,NodeCount,ISCAN
    INTEGER               :: N1
    REAL(SP)              :: R1,R2
    !==============================================================================|

    !
    !----------------Determine Number of Nodes -------------------------------

    ISCAN = SCAN_FILE(SPONGEUNIT,"Sponge Node Number",ISCAL = NSPONGE)
    IF(ISCAN /= 0) then
       write(temp,'(I2)') ISCAN
       call fatal_error('Improper formatting of SPONGE FILE: ISCAN ERROR&
            &# '//trim(temp),&
            & 'The header must contain: "Sponge Node Number ="', &
            & 'Followed by an integer number of nodes where the sponge&
            & value is set')
    END IF

!----------------Read SPONGE Array----------------------------------!
!

    if(NSPONGE==0) then
       if(DBG_SET(dbg_log)) then
          WRITE(IPT,*)'!  Finished Reading SPONGE File: NO SPONGE NODES'
          WRITE(IPT,*)'!  SPONGE NODES  =',NSPONGE
          WRITE(IPT,*)'!'
       end if
       return
    end if



   ! FIND FIRST LINE of )BC ARRAY
    rewind SPONGEUNIT
    DO WHILE(.TRUE.)
       READ(SPONGEUNIT,*,IOSTAT=IOS,END=99)N1,R1,R2
       if (IOS == 0) then
          BackSpace SPONGEUNIT
          exit
       end if

       CYCLE
       
99 Call FATAL_ERROR('Improper formatting of SPONGE FILE:',&
           &'Reached end of file with out finding sponge data?',&
           &'FORMAT: GLBNODE#  RADIUS SPGVAL (INT, REAL, REAL)')
          
    END DO

    ALLOCATE(N_SPG(NSPONGE)); N_SPG = 0
    ALLOCATE(R_SPG(NSPONGE)); R_SPG = 0.0
    ALLOCATE(C_SPG(NSPONGE)); C_SPG = 0.0


    
    I = 0
    DO WHILE(.TRUE.)

       READ(SPONGEUNIT,*,IOSTAT=IOS) N1,R1,R2
       IF(IOS < 0) exit
       
       I = I + 1
       IF(I > NSPONGE) CALL FATAL_ERROR('Number of rows of data in the SPONGE file &
            &exceeds the stated number in the header ?')

       IF( 1 > N1 .or. N1 > MGL) then
          write(temp,'(I8)') I
          write(temp2,'(I8)') MGL
          CALL FATAL_ERROR('SPONGE NODE NUMBER'//trim(temp)//&
               & 'IS NOT IN THE GLOBAL DOMAIN',&
               & 'CHECK INPUT FILE AND ENSURE SPONGE NODE# is <= '//trim(temp2))
       END IF
       N_SPG(I) = N1
       R_SPG(I) = R1
       C_SPG(I) = R2

    END DO

    if ( I .NE. NSPONGE) &
         & CALL FATAL_ERROR('The number of rows of data in the sponge file does&
         & not equal the number of nodes in the header?')

   if(DBG_SET(dbg_log)) then
      WRITE(IPT,*)'!  Finished Reading Sponge File'
      WRITE(IPT,*)'!  SPONGE NODES  =        :',NSPONGE
      WRITE(IPT,*)'!'
   end if


    
 END SUBROUTINE READ_COLDSTART_SPONGE


 SUBROUTINE READ_COLDSTART_SIGMA
    USE CONTROL ! Several variables for sigma coords!
    IMPLICIT NONE
    CHARACTER(LEN=80)    :: INPLINE,temp,temp2
    INTEGER              :: I,J,IOS,NodeCount, ISCAN
    INTEGER              :: N1, N2, N3

    !==============================================================================|
    
    ! Get data from Sigma file for the following variables - depends
    ! on coordinate type, some may not exist!

    ! GET SIGMA COORDINATE TYPE
    ISCAN = SCAN_FILE(SIGMAUNIT,"NUMBER OF SIGMA LEVELS",ISCAL = KB)
    IF(ISCAN /= 0) then
       write(temp,'(I2)') ISCAN
       call fatal_error('Improper formatting of SIGMA FILE: ISCAN ERROR&
            &# '//trim(temp),&
            & 'The header must contain: "NUMBER OF SIGMA LEVELS"', &
            & 'Followed by an integer number of levels.')
    END IF


    ! GET SIGMA COORDINATE TYPE
    ISCAN = SCAN_FILE(SIGMAUNIT,"SIGMA COORDINATE TYPE",CVAL = STYPE)
    IF(ISCAN /= 0) then
       write(temp,'(I2)') ISCAN
       call fatal_error('Improper formatting of SIGMA FILE: ISCAN ERROR&
            &# '//trim(temp),&
            & 'The header must contain: "SIGMA COORDINATE TYPE"', &
            & 'Followed by one of four defined types:',&
            & '"UNIFORM" or "GEOMETRIC" or "TANH" or "GENERALIZED" ')
    END IF

    ! SELECT CASE BASED ON SIGMA COORDINATE TYPE
    select case(TRIM(STYPE))
! DEGENERATE CASE OF GEOMETRIC SIGMA COORDINATES - UNIFROM DISTRIBUTION
    case(TRIM(STYPE_UNIFORM))  
       P_SIGMA = 1.0
       if(DBG_SET(dbg_log)) then
          WRITE(IPT,*)'!  Finished Reading SIGMA File'
          WRITE(IPT,*)'!  SIGMA COORDINATE TYPE = : '//trim(STYPE)
          WRITE(IPT,*)'!  P_SIGMA (UNIFORM)     = : ',P_SIGMA
          WRITE(IPT,*)'!  # OF SIGMA LEVELS(KB) = : ',KB
          WRITE(IPT,*)'!'
       end if

! GEOMETRIC SIGMA COORDINATES- P_SIGMA = 2 => Quadratic distribution
    case(TRIM(STYPE_GEOMETRIC))

       ISCAN = SCAN_FILE(SIGMAUNIT,"SIGMA POWER",FSCAL = P_SIGMA)
       IF(ISCAN /= 0) then
          write(temp,'(I2)') ISCAN
          call fatal_error('Improper formatting of SIGMA FILE: ISCAN ERROR&
               &# '//trim(temp),&
               & 'For GEOMETRIC SIGMA COORDINATE TYPE', &
               & 'The header must conatain "SIGMA POWER"',&
               & 'Followed by a real value (1.0 is uniform sigma coordinates)')
       END IF
       
       
       if(DBG_SET(dbg_log)) then
          WRITE(IPT,*)'!  Finished Reading SIGMA File'
          WRITE(IPT,*)'!  SIGMA COORDINATE TYPE = : '//trim(STYPE)
          WRITE(IPT,*)'!  P_SIGMA               = : ',P_SIGMA
          WRITE(IPT,*)'!  # OF SIGMA LEVELS(KB) = : ',KB
          WRITE(IPT,*)'!'
       end if
              
! HYPERBOLIC TANGENT DISTRIBUTION OF SURFACE AND BOTTOM INTESIFIED LAYERS
    case(TRIM(STYPE_TANH))
       
       ISCAN = SCAN_FILE(SIGMAUNIT,"DU",FSCAL = DU2)
       IF(ISCAN /= 0) then
          write(temp,'(I2)') ISCAN
          call fatal_error('Improper formatting of SIGMA FILE: ISCAN ERROR&
               &# '//trim(temp),&
               & 'For TANH SIGMA COORDINATE TYPE', &
               & 'The header must conatain "DU"',&
               & 'Followed by a real value (See set_sigma.F')
       END IF

       ISCAN = SCAN_FILE(SIGMAUNIT,"DL",FSCAL = DL2)
       IF(ISCAN /= 0) then
          write(temp,'(I2)') ISCAN
          call fatal_error('Improper formatting of SIGMA FILE: ISCAN ERROR&
               &# '//trim(temp),&
               & 'For TANH SIGMA COORDINATE TYPE', &
               & 'The header must conatain "DL"',&
               & 'Followed by a real value (See set_sigma.F')
       END IF
       
       if(DBG_SET(dbg_log)) then
          WRITE(IPT,*)'!  Finished Reading SIGMA File'
          WRITE(IPT,*)'!  SIGMA COORDINATE TYPE = : '//trim(STYPE)
          WRITE(IPT,*)'!  # OF SIGMA LEVELS(KB) = : ',KB
          WRITE(IPT,*)'!  DU                    = : ',DU2
          WRITE(IPT,*)'!  DL                    = : ',DL2
          WRITE(IPT,*)'!'
       end if

! A SPATIALLY DEPENDENT DISTRIBUTION OF LAYER THICKNESS BASED ON DEPTH
    case(TRIM(STYPE_GENERALIZED))

       ISCAN = SCAN_FILE(SIGMAUNIT,"DU",FSCAL = DUU)
       IF(ISCAN /= 0) then
          write(temp,'(I2)') ISCAN
          call fatal_error('Improper formatting of SIGMA FILE: ISCAN ERROR&
               &# '//trim(temp),&
               & 'For GENERALIZED SIGMA COORDINATE TYPE', &
               & 'The header must conatain "DU"',&
               & 'Followed by a real value (See set_sigma.F')
       END IF

       ISCAN = SCAN_FILE(SIGMAUNIT,"DL",FSCAL = DLL)
       IF(ISCAN /= 0) then
          write(temp,'(I2)') ISCAN
          call fatal_error('Improper formatting of SIGMA FILE: ISCAN ERROR&
               &# '//trim(temp),&
               & 'For GENERALIZED SIGMA COORDINATE TYPE', &
               & 'The header must conatain "DL"',&
               & 'Followed by a real value (See set_sigma.F')
       END IF

       ISCAN = SCAN_FILE(SIGMAUNIT,"MIN CONSTANT DEPTH",FSCAL = HMIN1)
       IF(ISCAN /= 0) then
          write(temp,'(I2)') ISCAN
          call fatal_error('Improper formatting of SIGMA FILE: ISCAN ERROR&
               &# '//trim(temp),&
               & 'For GENERALIZED SIGMA COORDINATE TYPE', &
               & 'The header must conatain "MIN CONSTANT DEPTH"',&
               & 'Followed by a real value (See set_sigma.F')
       END IF

       ISCAN = SCAN_FILE(SIGMAUNIT,"KU",ISCAL = KU)
       IF(ISCAN /= 0) then
          write(temp,'(I2)') ISCAN
          call fatal_error('Improper formatting of SIGMA FILE: ISCAN ERROR&
               &# '//trim(temp),&
               & 'For GENERALIZED SIGMA COORDINATE TYPE', &
               & 'The header must conatain "KU"',&
               & 'Followed by a real value (See set_sigma.F')
       END IF

       ISCAN = SCAN_FILE(SIGMAUNIT,"KL",ISCAL = KL)
       IF(ISCAN /= 0) then
          write(temp,'(I2)') ISCAN
          call fatal_error('Improper formatting of SIGMA FILE: ISCAN ERROR&
               &# '//trim(temp),&
               & 'For GENERALIZED SIGMA COORDINATE TYPE', &
               & 'The header must conatain "KL"',&
               & 'Followed by a real value (See set_sigma.F')
       END IF

!------------------------------------------------------------------------------|
!     "ZKU"   !! 
!------------------------------------------------------------------------------|
   IF(KU .ge. 1 .and. KU .LE. 150)THEN

     ALLOCATE(ZKU(KU)); ZKU=0.0_SP
     ISCAN = SCAN_FILE(SIGMAUNIT,"ZKU",FVEC = ZKU ,NSZE = N1)
     IF(ISCAN /= 0) then
        write(temp,'(I2)') ISCAN
        call fatal_error('Improper formatting of SIGMA FILE: ISCAN ERROR&
             &# '//trim(temp),&
             & 'For GENERALIZED SIGMA COORDINATE TYPE', &
             & 'The header must conatain "ZKU"',&
             & 'Followed by a real values (See set_sigma.F')
     END IF
     
     
     IF(N1 /= KU)THEN
        call fatal_error('Improper formatting of SIGMA FILE:',&
             & 'For GENERALIZED SIGMA COORDINATE TYPE', &
             & 'THE NUMBER OF SPECIFIED DEPTHS IN ZKU IS NOT EQUAL TO KU') 
     END IF
     


  ELSE IF( KU .NE. 0) THEN
     call fatal_error('Improper formatting of SIGMA FILE:',&
          & 'For GENERALIZED SIGMA COORDINATE TYPE', &
          & 'Requirement: 1<= KL <= 150;') 
     
  END IF
!------------------------------------------------------------------------------|
!     "ZKL"   !! 
!--------------------------------------------------------------------
   IF(KL .ge. 1 .and. KL .LE. 150)THEN

     ALLOCATE(ZKL(KL)); ZKL=0.0_SP
     ISCAN = SCAN_FILE(SIGMAUNIT,"ZKL",FVEC = ZKL ,NSZE = N1)
     IF(ISCAN /= 0) then
        write(temp,'(I2)') ISCAN
        call fatal_error('Improper formatting of SIGMA FILE: ISCAN ERROR&
             &# '//trim(temp),&
             & 'For GENERALIZED SIGMA COORDINATE TYPE', &
             & 'The header must conatain "ZKL"',&
             & 'Followed by a real values (See set_sigma.F')
     END IF
     
     
     IF(N1 /= KL)THEN
        call fatal_error('Improper formatting of SIGMA FILE:',&
             & 'For GENERALIZED SIGMA COORDINATE TYPE', &
             & 'THE NUMBER OF SPECIFIED DEPTHS IN ZKL IS NOT EQUAL TO KL') 
     END IF
     
  ELSE IF (KL .NE. 0) THEN
     call fatal_error('Improper formatting of SIGMA FILE:',&
          & 'For GENERALIZED SIGMA COORDINATE TYPE', &
          & 'Requirement: 1<= KL <= 150;') 
     
  END IF

   ! END GET VARIABLES
  
  if(DBG_SET(dbg_log)) then
     WRITE(IPT,*)'!  Finished Reading SIGMA File'
     WRITE(IPT,*)'!  SIGMA COORDINATE TYPE = : '//trim(STYPE)
     WRITE(IPT,*)'!  # OF SIGMA LEVELS(KB) = : ',KB
     WRITE(IPT,*)'!  DU                    = : ',DUU
     WRITE(IPT,*)'!  DL                    = : ',DLL
     WRITE(IPT,*)'!  MIN CONSTANT DEPTH    = : ',HMIN1
     WRITE(IPT,*)'!  KU                    = : ',KU
     WRITE(IPT,*)'!  KL                    = : ',KL
     WRITE(IPT,*)'!'
  end if
  
  
  
case default
   call fatal_error('Improper formatting of SIGMA FILE',&
        & 'Allowed SIGMA COORDINATE TYPEs are:',&
         & '"UNIFORM" or "GEOMETRIC" or "TANH" or "GENERALIZED"',&
         & 'See Set_Sigma.F for a description')
    end select
    
    
 END SUBROUTINE READ_COLDSTART_SIGMA

  SUBROUTINE HelpTxt(IPT)
    implicit none
    INTEGER, INTENT(IN) :: IPT
    write(IPT,*) "Need to put something here!"
    write(IPT,*) "This is not a very helpful help message!"
    write(IPT,*) "LONG INPUT OPTIONS"
    write(IPT,*) "--HELP => PRINT THIS MESSAGE"
    write(IPT,*) "--CASENAME=<YOUR_CASE> (REQUIRED)"
    write(IPT,*) "--CREATE_NAMELIST => PRINT BLANK NAMELIST AND RETURN"
    write(IPT,*) "--LOGFILE=<FILENAME> => TO OUTPUT TO A LOG FILE"
    write(IPT,*) "--CRASHRESTART  => RUN FROM CURRENT TIME IN RESTART FILE"
    write(IPT,*) "SHORT INPUT OPTIONS"
    write(IPT,*) "-V => PRINT FVCOM VERSION INFO AND RETURN"
    write(IPT,*) "-H => PRINT THIS MESSAGE AND RETURN"
    write(IPT,*) ""
    write(IPT,*) "DEBUG LEVELS"
    write(IPT,*) "--dbg=0 => DBG LOG (DEFAULT"
    write(IPT,*) "--dbg=1 => DBG IO FILENAMES"
    write(IPT,*) "--dbg=2 => DBG SCALARS"
    write(IPT,*) "--dbg=4 => DBG SUBROUTINE NAMES"
    write(IPT,*) "--dbg=5 => DBG SUBROUTINE IO"
    write(IPT,*) "--dbg=6 => DBG VECTORS"
    write(IPT,*) "--dbg=7 => DBG EVERYTHING"
    write(IPT,*) "--dbg_par => WRITE LOG FOR EACH PROCESSOR"
    write(IPT,*) ""

  END SUBROUTINE HelpTxt

  SUBROUTINE IS_TRI_CW(X, Y, IS_CW)

  !===================================================
  ! Subroutine to check if a triangle is clockwise.
  ! Input  : x(3), the x-coordinate 
  !          y(3), the y-coordinate
  ! Output : is_cw, integer result. 1, clockwise
  !                                -1, anti-clockwise
  !                                 0, line
  ! Siqi Li, SMAST.
  ! 
  !===================================================
    implicit none
    !
    real(sp), intent(in) :: x(3), y(3)
    integer, intent(out) :: is_cw
    real(sp)             :: cross_prod
    !
    cross_prod=(x(2)-x(1))*(y(3)-y(1))-(y(2)-y(1))*(x(3)-x(1))
    !
    if (abs(cross_prod)<1e-6) then
      is_cw=0
!      print*, 'This is a line, rather than a triangle.'
    else
      if (cross_prod>0) then
        is_cw=-1
!        print*, 'This triangle is anti-clockwise.'
      else
        is_cw=1
!        print*, 'This triangle is clockwise.'
      end if
    end if
    !
  END SUBROUTINE IS_TRI_CW
 
END MODULE MOD_INPUT
