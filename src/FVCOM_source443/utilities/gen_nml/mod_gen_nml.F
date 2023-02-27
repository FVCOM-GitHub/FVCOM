module mod_gen_nml
  use lims
  use mod_utils
  use mod_par
  use control
  use mod_prec
  use mod_ncll
  use mod_nctools
  USE MOD_TIME
  implicit none

  

  CHARACTER(len=80), PARAMETER :: mode_forecast = "forecast"
  CHARACTER(len=80), PARAMETER :: mode_hindcast = "hindcast"
  CHARACTER(len=80) :: mode
  CHARACTER(LEN=120) :: filename

  ! DUMMY VARIABLE - REMOVE ME!
  CHARACTER(LEN=120) :: BOTTOM_ROUGHNESS_FUNCTION


contains



  SUBROUTINE GET_COMMANDLINE(CVS_ID,CVS_Date,CVS_Name,CVS_Revision)
    use mod_sng


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

    ! Main code
    call ftn_strini(cmd_ln) ! [sng] sng(1:len)=NUL

    call ftn_cmd_ln_sng(cmd_ln) ! [sng] Re-construct command-line into single string
    call ftn_prg_ID_mk(CVS_Id,CVS_Revision,CVS_Date,prg_ID) ! [sng] Program ID

    arg_nbr=command_argument_count() ! [nbr] Number of command-line arguments

    if (arg_nbr .LE. 0 ) then
       if(MSR) WRITE(IPT,*) "You must specify an arugument:"
       if(MSR) Call MYHelpTxt
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

             !          else if (opt_sng == "dbg_par" .or.opt_sng == "Dbg_Par"&
             !               & .or.opt_sng == "DBG_PAR") then

             !             dbg_par = .true.

          else if (opt_sng == "NMLNAME" .or.opt_sng == "NmlName"&
               & .or.opt_sng == "nmlname") then

             call ftn_arg_get(arg_idx,arg_val,filename) ! [sng] Input file
             filename=filename(1:ftn_strlen(filename))
             ! Convert back to a fortran string!

          else if (opt_sng == "STARTDATE" .or.opt_sng == "startdate"&
               & .or.opt_sng == "StartDate") then

             call ftn_arg_get(arg_idx,arg_val,START_DATE) ! [sng] Input file
             START_DATE=START_DATE(1:ftn_strlen(START_DATE))
             ! Convert back to a fortran string!

          else if (opt_sng == "ENDDATE" .or.opt_sng == "enddate"&
               & .or.opt_sng == "EndDate") then

             call ftn_arg_get(arg_idx,arg_val,END_DATE) ! [sng] Input file
             END_DATE=END_DATE(1:ftn_strlen(END_DATE))
             ! Convert back to a fortran string!


          else if (opt_sng == "MODE" .or.opt_sng == "mode"&
               & .or.opt_sng == "Mode") then

             call ftn_arg_get(arg_idx,arg_val,MODE) ! [sng] Input file
             MODE=MODE(1:ftn_strlen(MODE))
             ! Convert back to a fortran string!

          else if (opt_sng == "help" .or.opt_sng == "HELP" .or. opt_sng&
               & == "Help") then

             if(MSR) call MYHelpTxt
             call PSHUTDOWN

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

          if(MSR) Call MYHelpTxt
          Call PSHUTDOWN

       else ! Option not recognized
          arg_idx=arg_idx-1 ! [idx] Counting index
          if(MSR) call ftn_getarg_err(arg_idx,arg_val) ! [sbr] Error handler for getarg()
       endif ! endif arg_val


    end do ! end while (arg_idx <= arg_nbr)

    CALL dbg_init(IPT_BASE,.false.)

  END SUBROUTINE GET_COMMANDLINE

  SUBROUTINE MYHELPTXT
    IMPLICIT NONE


    write(IPT,*) "! ARGUMENTS FOR gen_nml:"
    write(IPT,*) "! "
    write(IPT,*) "! NMLNMAE : the name of the namelist file to create"
    write(IPT,*) "! "
    write(IPT,*) "! MODE : 'forecast' or 'hindcast'"
    write(IPT,*) "! "
    write(IPT,*) "! StartDate : the starting date of the fvcom run"
    write(IPT,*) "! "
    write(IPT,*) "! EndDate : the ending date of the fvcom run"
    write(IPT,*) "! "

    WRITE(IPT,*) "! NOTES: Do not run this program in parallel!"


  END SUBROUTINE MYHELPTXT


  SUBROUTINE CHECK_INPUT
    IMPLICIT NONE
    INTEGER :: STATUS

    DATE_FORMAT="YMD"
    TIMEZONE="UTC"
    StartTime= READ_DATETIME(START_DATE,DATE_FORMAT,TIMEZONE,status)
    if (.not. status) &
         & Call Fatal_Error("Could not read the date string START_DATE: "//trim(START_DATE))
        
    ! GET THE END TIME
    EndTime = READ_DATETIME(END_DATE,DATE_FORMAT,TIMEZONE,status)
    if (.not. status) &
         & Call Fatal_Error("Could not read the date string END_DATE: "//trim(END_DATE))
    
    IF (ENDTIME <= STARTTIME) CALL FATAL_ERROR("INVALID START AND END TIME!")

    
  END SUBROUTINE CHECK_INPUT


  SUBROUTINE NML_HINDCAST
    USE MOD_ASSIM
    IMPLICIT NONE
    CHARACTER(len=8), parameter :: none = "'none'"

    CALL FOPEN(NMLUNIT,trim(FILENAME),'ofr')

    !--Parameters in NameList NML_CASE
    CASE_TITLE = "'FVCOM GOM NECOFS HINDCAST UPDATE'"
    TIMEZONE = "'UTC'"
    DATE_FORMAT="'YMD'"
    START_DATE = "'"//TRIM(START_DATE)//"'"
    END_DATE = "'"//TRIM(END_DATE)//"'"
    ! START_DATE SET AT COMMAND LINE
    ! END_DATE SET AT COMMAND LINE

    write(UNIT=NMLUNIT,NML=NML_CASE)

    STARTUP_TYPE = "'forecast'"
    STARTUP_FILE = none
    STARTUP_UV_TYPE= none
    STARTUP_TURB_TYPE = none
    STARTUP_TS_TYPE = none
    write(UNIT=NMLUNIT,NML=NML_STARTUP)

    !--Parameters in NameList NML_IO
    INPUT_DIR = "'/data01/necofs/FVCOM/INPUT'"
    OUTPUT_DIR = "'/data01/necofs/FVCOM/RESULTS'"
    IREPORT = 60

    write(UNIT=NMLUNIT,NML=NML_IO)

    !--Parameters in NameList NML_INTEGRATION
    ExtStep_Seconds = 4.0
    ISplit = 10
    IRamp = 0 
    Static_SSH_Adj = 0.0
    Min_Depth = 0.05

    write(UNIT=NMLUNIT,NML=NML_INTEGRATION)

    !--Parameters in NameList NML_RESTART
    RST_ON            = .TRUE.
    RST_FIRST_OUT     = START_DATE
    RST_OUT_INTERVAL  = "'days= 1.0'"
    RST_OUTPUT_STACK  = 0
    
    write(UNIT=NMLUNIT,NML=NML_RESTART)

    !--Parameters in NameList NML_NETCDF
    NC_ON = .True.
    NC_FIRST_OUT = START_DATE
    NC_OUT_INTERVAL = "'seconds= 3600.0'"
    NC_OUTPUT_STACK = 0
    NC_VELOCITY     = .TRUE.
    NC_GRID_METRICS = .TRUE.
    NC_FILE_DATE    = .TRUE.
    NC_SALT_TEMP    = .TRUE.
    NC_TURBULENCE   = .TRUE.
    NC_VERTICAL_VEL = .TRUE.
    NC_AVERAGE_VEL  = .TRUE.
    NC_WIND_STRESS  = .TRUE.
    NC_EVAP_PRECIP  = .FALSE.
    NC_SURFACE_HEAT = .TRUE.

    write(UNIT=NMLUNIT,NML=NML_NETCDF)

    ! OUTPUT VARIABLES DEFAULT TO OFF

    !--Parameters in NameList NML_NETCDF_AV
    NCAV_ON = .False.
    NCAV_FIRST_OUT = none
    NCAV_OUT_INTERVAL = none
    write(UNIT=NMLUNIT,NML=NML_NETCDF_AV)

    !--Parameters in NameList NML_PHYSICS
    HORIZONTAL_MIXING_TYPE = "'closure'"
    HORIZONTAL_MIXING_COEFFICIENT = 4.0E-1_SP
    HORIZONTAL_PRANDTL_NUMBER = 1.0_SP

    VERTICAL_MIXING_TYPE = "'closure'"
    VERTICAL_MIXING_COEFFICIENT = 1.0E-5
    VERTICAL_PRANDTL_NUMBER = 1.0_SP


    BOTTOM_ROUGHNESS_MINIMUM = 0.0025_SP
    BOTTOM_ROUGHNESS_LENGTHSCALE = 0.03_SP
    BOTTOM_ROUGHNESS_FUNCTION  = "'user_defined'"
    SEA_WATER_DENSITY_FUNCTION = "'dens2'"
    BAROCLINIC_PRESSURE_GRADIENT = "'sigma levels'" 
  
    CONVECTIVE_OVERTURNING          = .TRUE.
    SCALAR_POSITIVITY_CONTROL       = .TRUE.
    BAROTROPIC                      = .FALSE.
    TEMPERATURE_ACTIVE              = .TRUE.
    SALINITY_ACTIVE                 = .TRUE.
    
    write(UNIT=NMLUNIT,NML=NML_PHYSICS)

    !--Parameters in NameList NML_SURFACE_FORCING
    WIND_ON = .TRUE.
    WIND_TYPE = "'stress'"
    WIND_FILE = "'wrf_hnd.nc'"
    WIND_KIND = "'variable'"
    WIND_X = 0.0
    WIND_Y = 0.0
    HEATING_ON = .TRUE.
    HEATING_TYPE =  "'flux'"
    HEATING_KIND = "'variable'"
    HEATING_FILE = "'wrf_hnd.nc'"
    HEATING_LONGWAVE_LENGTHSCALE = 6.3_SP
    HEATING_LONGWAVE_PERCTAGE = 0.78_SP
    HEATING_SHORTWAVE_LENGTHSCALE= 1.4_SP
    PRECIPITATION_ON = .FALSE.
    PRECIPITATION_KIND = none
    PRECIPITATION_FILE = none
    write(UNIT=NMLUNIT,NML=NML_SURFACE_FORCING)

    !--Parameters in NameList NML_RIVER_TYPE
    RIVER_NUMBER = 9
    RIVER_TS_SETTING = "'calculated'"
    RIVER_INFLOW_LOCATION = "'node'"

    write(UNIT=NMLUNIT,NML=NML_RIVER_TYPE)


    !--Parameters in NameList NML_RIVERS
    RIVER_NAME                  = "'St_John_basin_ME'"
    RIVER_FILE                  = "'RIVER/St_John_basin_ME.nc'"
    RIVER_GRID_LOCATION         = 30315
    RIVER_VERTICAL_DISTRIBUTION = "'uniform'"
    write(UNIT=NMLUNIT,NML=NML_RIVER)
    
    RIVER_NAME                  = "'St_John_dickey_ME'"
    RIVER_FILE                  = "'RIVER/St_John_dickey_ME.nc'"
    RIVER_GRID_LOCATION         = 30315
    RIVER_VERTICAL_DISTRIBUTION = "'uniform'"
    write(UNIT=NMLUNIT,NML=NML_RIVER)

    RIVER_NAME                  = "'St_John_fish_ME'"
    RIVER_FILE                  = "'RIVER/St_John_fish_ME.nc'"
    RIVER_GRID_LOCATION         = 30315
    RIVER_VERTICAL_DISTRIBUTION = "'uniform'"
    write(UNIT=NMLUNIT,NML=NML_RIVER)

    RIVER_NAME                  = "'St_Croix_ME'"
    RIVER_FILE                  = "'RIVER/St_Croix_ME.nc'"
    RIVER_GRID_LOCATION         = 30139
    RIVER_VERTICAL_DISTRIBUTION = "'uniform'"
    write(UNIT=NMLUNIT,NML=NML_RIVER)

    RIVER_NAME                  = "'Penobscot_ME'"
    RIVER_FILE                  = "'RIVER/Penobscot_ME.nc'"
    RIVER_GRID_LOCATION         = 29435
    RIVER_VERTICAL_DISTRIBUTION = "'uniform'"
    write(UNIT=NMLUNIT,NML=NML_RIVER)

    RIVER_NAME                  = "'Kennebec_MA'"
    RIVER_FILE                  = "'RIVER/Kennebec_MA.nc'"
    RIVER_GRID_LOCATION         = 32623
    RIVER_VERTICAL_DISTRIBUTION = "'uniform'"
    write(UNIT=NMLUNIT,NML=NML_RIVER)

    RIVER_NAME                  = "'Androscoggin_ME'"
    RIVER_FILE                  = "'RIVER/Androscoggin_ME.nc'"
    RIVER_GRID_LOCATION         = 32448
    RIVER_VERTICAL_DISTRIBUTION = "'uniform'"
    write(UNIT=NMLUNIT,NML=NML_RIVER)

    RIVER_NAME                  = "'Saco_ME'"
    RIVER_FILE                  = "'RIVER/Saco_ME.nc'"
    RIVER_GRID_LOCATION         = 24594
    RIVER_VERTICAL_DISTRIBUTION = "'uniform'"
    write(UNIT=NMLUNIT,NML=NML_RIVER)

    RIVER_NAME                  = "'Merrimac_MA'"
    RIVER_FILE                  = "'RIVER/Merrimac_MA.nc'"
    RIVER_GRID_LOCATION         = 26525
    RIVER_VERTICAL_DISTRIBUTION = "'uniform'"
    write(UNIT=NMLUNIT,NML=NML_RIVER)




    !--Parameters in NameList NML_OPEN_BOUNDARY
    OBC_ON = .TRUE.
    OBC_NODE_LIST_FILE = none
    OBC_ELEVATION_FORCING_ON = .TRUE.
    OBC_ELEVATION_FILE = "'julian_obc.nc'"
    OBC_TS_TYPE      = 3
    OBC_TEMP_NUDGING = .TRUE.
    OBC_TEMP_FILE = "'tsobc.nc'"
    OBC_TEMP_NUDGING_TIMESCALE = 0.1
    OBC_SALT_NUDGING = .TRUE.
    OBC_SALT_FILE = "'tsobc.nc'"
    OBC_SALT_NUDGING_TIMESCALE = 0.1
    OBC_MEANFLOW = .FALSE.
    OBC_MEANFLOW_FILE = none
    OBC_LONGSHORE_FLOW_ON = .TRUE.
    OBC_LONGSHORE_FLOW_FILE = none
    write(UNIT=NMLUNIT,NML=NML_OPEN_BOUNDARY_CONTROL)

    !--Parameters in NameList GRID_COORDINATES
    GRID_FILE = none
    SIGMA_LEVELS_FILE = none
    DEPTH_FILE = none
    CORIOLIS_FILE = none
    SPONGE_FILE = none
    GRID_FILE_UNITS = "'meters'"
    PROJECTION_REFERENCE = "'init=nad83:1802'"
    write(UNIT=NMLUNIT,NML=NML_GRID_COORDINATES)

    !--Parameters in NameList NML_GROUNDWATER
    GROUNDWATER_ON = .False.
    GROUNDWATER_FILE = none
    
    write(UNIT=NMLUNIT,NML=NML_GROUNDWATER)

    !--Parameters in NameList NML_LAG_PART
    LAG_PARTICLES_ON = .False.
    LAG_START_FILE = none
    LAG_OUT_FILE = none
    LAG_RESTART_FILE = none
    LAG_OUT_INTERVAL = none
    LAG_SCAL_CHOICE = none
    write(UNIT=NMLUNIT,NML=NML_LAG)

    !--Parameters in NameList NML_ADDITIONAL_MODELS
    DATA_ASSIMILATION  = .TRUE.
    DATA_ASSIMILATION_FILE = "'"//TRIM(filename)//"'"
!    WATER_QUALITY_MODEL_FILE = none
!    BIOLOGICAL_MODEL_FILE = none
    SEDIMENT_MODEL_FILE = none

    ICE_MODEL = .FALSE.
    ICE_FORCING_FILE = none
    ICE_FORCING_KIND = none
    ICE_SEA_LEVEL_PRESSURE = 0.0_SP
    ICE_AIR_TEMP           = 0.0_SP
    ICE_SPEC_HUMIDITY      = 0.0_SP
    ICE_CLOUD_COVER        = 0.0_SP

    ICING_MODEL = .TRUE.
    ICING_FORCING_FILE = "'wrf_hnd.nc'"
    ICING_FORCING_KIND = "'variable'"
    ICING_AIR_TEMP = 0.0_SP
    ICING_WSPD     = 0.0_SP

    write(UNIT=NMLUNIT,NML=NML_ADDITIONAL_MODELS)

    !--Parameters in NameList NML_PROBES
    PROBES_ON=.False.
    PROBES_NUMBER=-1
    PROBES_FILE=none
    write(UNIT=NMLUNIT,NML=NML_PROBES)


    !
    !--Data Assimilation Parameters for SST Assimilation
    !
    SST_ASSIM = .FALSE.
    SST_ASSIM_FILE = none
    write(UNIT=NMLUNIT,NML=NML_SST_ASSIMILATION)

    !
    !--Data Assimilation Parameters for SST GRID Assimilation
    !
    SSTGRD_ASSIM           = .TRUE.
    SSTGRD_ASSIM_FILE      = "'sst.nc'"
    SSTGRD_WEIGHT_MAX      = 0.5
    SSTGRD_TIMESCALE       = 3.e-3
    SSTGRD_TIME_WINDOW     = 1800.0
    SSTGRD_N_PER_INTERVAL  = 24
    write(UNIT=NMLUNIT,NML=NML_SSTGRD_ASSIMILATION)

    !
    !--Data Assimilation Parameters for Current Assimilation
    !
    CUR_ASSIM = .FALSE.
    CUR_ASSIM_FILE = none
    write(UNIT=NMLUNIT,NML=NML_CUR_ASSIMILATION)
    !
    !--Data Assimilation Parameters for Temp/Salinity Data Assimilation
    !
    TS_ASSIM = .FALSE.
    TS_ASSIM_FILE = none
    write(UNIT=NMLUNIT,NML=NML_TS_ASSIMILATION)

  END SUBROUTINE NML_HINDCAST

  SUBROUTINE NML_FORECAST
    IMPLICIT NONE
    CHARACTER(len=8), parameter :: none = "'none'"

    CALL FOPEN(NMLUNIT,trim(FILENAME),'ofr')

    !--Parameters in NameList NML_CASE
    CASE_TITLE = "'FVCOM GOM NECOFS FORECAST UPDATE'"
    TIMEZONE = "'UTC'"
    DATE_FORMAT="'YMD'"
    START_DATE = "'"//TRIM(START_DATE)//"'"
    END_DATE = "'"//TRIM(END_DATE)//"'"

    write(UNIT=NMLUNIT,NML=NML_CASE)

    STARTUP_TYPE = "'forecast'"
    STARTUP_FILE = none
    STARTUP_UV_TYPE= none
    STARTUP_TURB_TYPE = none
    STARTUP_TS_TYPE = none
    write(UNIT=NMLUNIT,NML=NML_STARTUP)

    !--Parameters in NameList NML_IO
    INPUT_DIR = "'/data01/necofs/FVCOM/INPUT'"
    OUTPUT_DIR = "'/data01/necofs/FVCOM/RESULTS'"
    IREPORT = 60

    write(UNIT=NMLUNIT,NML=NML_IO)

    !--Parameters in NameList NML_INTEGRATION
    ExtStep_seconds = 4.0
    ISplit = 10
    IRamp = 0 
    Static_SSH_Adj = 0.0
    Min_Depth = 0.05

    write(UNIT=NMLUNIT,NML=NML_INTEGRATION)

    !--Parameters in NameList NML_RESTART
    RST_ON            = .TRUE.
    RST_FIRST_OUT     = START_DATE
    RST_OUT_INTERVAL  = "'days= 1.0'"
    RST_OUTPUT_STACK  = 0
    
    write(UNIT=NMLUNIT,NML=NML_RESTART)

    !--Parameters in NameList NML_NETCDF
    NC_ON = .True.
    NC_FIRST_OUT = START_DATE
    NC_OUT_INTERVAL = "'seconds= 3600.0'"
    NC_OUTPUT_STACK = 0
    NC_VELOCITY     = .TRUE.
    NC_GRID_METRICS = .TRUE.
    NC_FILE_DATE    = .TRUE.
    NC_SALT_TEMP    = .TRUE.
    NC_TURBULENCE   = .TRUE.
    NC_VERTICAL_VEL = .TRUE.
    NC_AVERAGE_VEL  = .TRUE.
    NC_WIND_STRESS  = .TRUE.
    NC_EVAP_PRECIP  = .FALSE.
    NC_SURFACE_HEAT = .TRUE.

    write(UNIT=NMLUNIT,NML=NML_NETCDF)

    ! OUTPUT VARIABLES DEFAULT TO OFF

    !--Parameters in NameList NML_NETCDF_AV
    NCAV_ON = .False.
    NCAV_FIRST_OUT = none
    NCAV_OUT_INTERVAL = none
    write(UNIT=NMLUNIT,NML=NML_NETCDF_AV)

    !--Parameters in NameList NML_PHYSICS
    HORIZONTAL_MIXING_TYPE = "'closure'"
    HORIZONTAL_MIXING_COEFFICIENT = 4.0E-1_SP
    HORIZONTAL_PRANDTL_NUMBER = 1.0_SP

    VERTICAL_MIXING_TYPE = "'closure'"
    VERTICAL_MIXING_COEFFICIENT = 1.0E-5
    VERTICAL_PRANDTL_NUMBER = 1.0_SP


    BOTTOM_ROUGHNESS_MINIMUM = 0.0025_SP
    BOTTOM_ROUGHNESS_LENGTHSCALE = 0.03_SP
    BOTTOM_ROUGHNESS_FUNCTION  = "'user_defined'"
    SEA_WATER_DENSITY_FUNCTION = "'dens2'"
    BAROCLINIC_PRESSURE_GRADIENT = "'sigma levels'" 
  
    CONVECTIVE_OVERTURNING          = .TRUE.
    SCALAR_POSITIVITY_CONTROL       = .TRUE.
    BAROTROPIC                      = .FALSE.
    TEMPERATURE_ACTIVE              = .TRUE.
    SALINITY_ACTIVE                 = .TRUE.
    
    write(UNIT=NMLUNIT,NML=NML_PHYSICS)

    !--Parameters in NameList NML_SURFACE_FORCING
    WIND_ON = .TRUE.
    WIND_TYPE = "'stress'"
    WIND_FILE = "'wrf_for.nc'"
    WIND_KIND = "'variable'"
    WIND_X = 0.0
    WIND_Y = 0.0
    HEATING_ON = .TRUE.
    HEATING_TYPE =  "'flux'"
    HEATING_KIND = "'variable'"
    HEATING_FILE = "'wrf_for.nc'"
    HEATING_LONGWAVE_LENGTHSCALE = 6.3_SP
    HEATING_LONGWAVE_PERCTAGE = 0.78_SP
    HEATING_SHORTWAVE_LENGTHSCALE= 1.4_SP
    PRECIPITATION_ON = .FALSE.
    PRECIPITATION_KIND = none
    PRECIPITATION_FILE = none
    write(UNIT=NMLUNIT,NML=NML_SURFACE_FORCING)

    !--Parameters in NameList NML_RIVER_TYPE
    RIVER_NUMBER = 9
    RIVER_TS_SETTING = "'calculated'"
    RIVER_INFLOW_LOCATION = "'node'"

    write(UNIT=NMLUNIT,NML=NML_RIVER_TYPE)


    !--Parameters in NameList NML_RIVERS
    RIVER_NAME                  = "'St_John_basin_ME'"
    RIVER_FILE                  = "'RIVER/St_John_basin_ME.nc'"
    RIVER_GRID_LOCATION         = 30315
    RIVER_VERTICAL_DISTRIBUTION = "'uniform'"
    write(UNIT=NMLUNIT,NML=NML_RIVER)
    
    RIVER_NAME                  = "'St_John_dickey_ME'"
    RIVER_FILE                  = "'RIVER/St_John_dickey_ME.nc'"
    RIVER_GRID_LOCATION         = 30315
    RIVER_VERTICAL_DISTRIBUTION = "'uniform'"
    write(UNIT=NMLUNIT,NML=NML_RIVER)

    RIVER_NAME                  = "'St_John_fish_ME'"
    RIVER_FILE                  = "'RIVER/St_John_fish_ME.nc'"
    RIVER_GRID_LOCATION         = 30315
    RIVER_VERTICAL_DISTRIBUTION = "'uniform'"
    write(UNIT=NMLUNIT,NML=NML_RIVER)

    RIVER_NAME                  = "'St_Croix_ME'"
    RIVER_FILE                  = "'RIVER/St_Croix_ME.nc'"
    RIVER_GRID_LOCATION         = 30139
    RIVER_VERTICAL_DISTRIBUTION = "'uniform'"
    write(UNIT=NMLUNIT,NML=NML_RIVER)

    RIVER_NAME                  = "'Penobscot_ME'"
    RIVER_FILE                  = "'RIVER/Penobscot_ME.nc'"
    RIVER_GRID_LOCATION         = 29435
    RIVER_VERTICAL_DISTRIBUTION = "'uniform'"
    write(UNIT=NMLUNIT,NML=NML_RIVER)

    RIVER_NAME                  = "'Kennebec_MA'"
    RIVER_FILE                  = "'RIVER/Kennebec_MA.nc'"
    RIVER_GRID_LOCATION         = 32623
    RIVER_VERTICAL_DISTRIBUTION = "'uniform'"
    write(UNIT=NMLUNIT,NML=NML_RIVER)

    RIVER_NAME                  = "'Androscoggin_ME'"
    RIVER_FILE                  = "'RIVER/Androscoggin_ME.nc'"
    RIVER_GRID_LOCATION         = 32448
    RIVER_VERTICAL_DISTRIBUTION = "'uniform'"
    write(UNIT=NMLUNIT,NML=NML_RIVER)

    RIVER_NAME                  = "'Saco_ME'"
    RIVER_FILE                  = "'RIVER/Saco_ME.nc'"
    RIVER_GRID_LOCATION         = 24594
    RIVER_VERTICAL_DISTRIBUTION = "'uniform'"
    write(UNIT=NMLUNIT,NML=NML_RIVER)

    RIVER_NAME                  = "'Merrimac_MA'"
    RIVER_FILE                  = "'RIVER/Merrimac_MA.nc'"
    RIVER_GRID_LOCATION         = 26525
    RIVER_VERTICAL_DISTRIBUTION = "'uniform'"
    write(UNIT=NMLUNIT,NML=NML_RIVER)




    !--Parameters in NameList NML_OPEN_BOUNDARY
    OBC_ON = .TRUE.
    OBC_NODE_LIST_FILE = none
    OBC_ELEVATION_FORCING_ON = .TRUE.
    OBC_ELEVATION_FILE = "'julian_obc.nc'"
    OBC_TS_TYPE      = 3
    OBC_TEMP_NUDGING = .TRUE.
    OBC_TEMP_FILE = "'tsobc.nc'"
    OBC_TEMP_NUDGING_TIMESCALE = 0.1
    OBC_SALT_NUDGING = .TRUE.
    OBC_SALT_FILE = "'tsobc.nc'"
    OBC_SALT_NUDGING_TIMESCALE = 0.1
    OBC_MEANFLOW = .FALSE.
    OBC_MEANFLOW_FILE = none
    OBC_LONGSHORE_FLOW_ON = .TRUE.
    OBC_LONGSHORE_FLOW_FILE = none
    write(UNIT=NMLUNIT,NML=NML_OPEN_BOUNDARY_CONTROL)

    !--Parameters in NameList GRID_COORDINATES
    GRID_FILE = none
    SIGMA_LEVELS_FILE = none
    DEPTH_FILE = none
    CORIOLIS_FILE = none
    SPONGE_FILE = none
    GRID_FILE_UNITS = "'meters'"
    PROJECTION_REFERENCE = "'init=nad83:1802'"
    write(UNIT=NMLUNIT,NML=NML_GRID_COORDINATES)

    !--Parameters in NameList NML_GROUNDWATER
    GROUNDWATER_ON = .False.
    GROUNDWATER_FILE = none
    
    write(UNIT=NMLUNIT,NML=NML_GROUNDWATER)

    !--Parameters in NameList NML_LAG_PART
    LAG_PARTICLES_ON = .False.
    LAG_START_FILE = none
    LAG_OUT_FILE = none
    LAG_RESTART_FILE = none
    LAG_OUT_INTERVAL = none
    LAG_SCAL_CHOICE = none
    write(UNIT=NMLUNIT,NML=NML_LAG)

    !--Parameters in NameList NML_PROBES
    PROBES_ON=.False.
    PROBES_NUMBER=-1
    PROBES_FILE=none
    write(UNIT=NMLUNIT,NML=NML_PROBES)

    !--Parameters in NameList NML_ADDITIONAL_MODELS
    DATA_ASSIMILATION  = .FALSE.
    DATA_ASSIMILATION_FILE  = none
!    WATER_QUALITY_MODEL_FILE = none
!    BIOLOGICAL_MODEL_FILE = none
    SEDIMENT_MODEL_FILE = none

    ICE_MODEL = .FALSE.
    ICE_FORCING_FILE = none
    ICE_FORCING_KIND = none
    ICE_SEA_LEVEL_PRESSURE = 0.0_SP
    ICE_AIR_TEMP           = 0.0_SP
    ICE_SPEC_HUMIDITY      = 0.0_SP
    ICE_CLOUD_COVER        = 0.0_SP

    ICING_MODEL = .TRUE.
    ICING_FORCING_FILE = "'wrf_for.nc'"
    ICING_FORCING_KIND = "'variable'"
    ICING_AIR_TEMP = 0.0_SP
    ICING_WSPD     = 0.0_SP

    write(UNIT=NMLUNIT,NML=NML_ADDITIONAL_MODELS)



  END SUBROUTINE NML_FORECAST


end module mod_gen_nml
