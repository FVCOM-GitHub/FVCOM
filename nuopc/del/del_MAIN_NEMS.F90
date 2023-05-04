#include "./ESMFVersionDefine.h"

#if ESMF_VERSION_MAJOR >= 8
! Do things for ESMF 8 and later
#else
#ifndef ESMF_VERSION_MAJOR
#error ESMF_VERSION_MAJOR is unset or empty; I do not know what version of ESMF you are using.
#else
! Do things for ESMF versions before 8
#endif
#endif

!-----------------------------------------------------------------------
!
      PROGRAM MAIN_NEMS
!
!-----------------------------------------------------------------------
!***  Main Program for NEMS.
!***  Define ESMF data types and procedures.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
! PROGRAM HISTORY LOG:
!   2007-       Black   - Modified from Wei-yu's version
!   2007-09     Black   - Create the Clock here.
!   2009-08     Colon   - Unified NEM-NMM & NEMS-GFS
!   2009-06-29  Black   - Modified for addition of NMM nesting;
!                         added new ATM Driver Component.
!   2009-09     Lu      - Add w3tage calls for resource statistics
!   2009-08     W. Yang - Ensemble GEFS Concurrency Code.
!   2010-03     Jovic/Black - Revised to create NEMS gridded component
!                             for new structure.
!   2010-04     Yang    - Add GEFS and GFS for the revised NEMS.
!   2010-11     Yang    - Add the "Generic Core" to NEMS
!   2011-02     Yang    - Updated to use both the ESMF 4.0.0rp2 library,
!                         ESMF 5 series library and the the 
!                         ESMF 3.1.0rp2 library.
!   2011-05     Theurich & Yang - Modified for using the ESMF 5.2.0r_beta_snapshot_07.
!   2011-10     Yang    - Modified for using the ESMF 5.2.0r library.
!   2013-07     Theurich - Macro based ESMF error handling
!   2016-11     Trahan  - Resource usage reporting
!
!-----------------------------------------------------------------------
!
      USE ESMF
!
!-----------------------------------------------------------------------
!***  USE the NEMS gridded component module.  Although it
!***  contains the calls to Register and the top level Initialize,
!***  Run, and Finalize, only the Register routine is public.
!-----------------------------------------------------------------------
!
       USE module_NEMS_GRID_COMP, ONLY: NEMS_REGISTER
!
!-----------------------------------------------------------------------
!***  The following module contains error-checking, and other utilities
!-----------------------------------------------------------------------
!
       USE module_NEMS_UTILS, ONLY: check_esmf_pet, err_msg, message_check
!
!-----------------------------------------------------------------------
!***  This module calculates resource usage across all ranks.
!-----------------------------------------------------------------------
!
       USE module_NEMS_Rusage,ONLY: NEMS_Rusage

#ifdef FRONT_FVCOM
       USE mod_driver,  ONLY: syear, smonth,sday, shour, sminute, ssecond
#endif 
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
      INCLUDE 'mpif.h'
!
!-----------------------------------------------------------------------
!***  Local Variables.
!-----------------------------------------------------------------------
!
      INTEGER :: MYPE                                                   &  !<-- The MPI task ID
                ,NHOURS_FCST                                            &  !<-- Length of forecast in hours
                ,NSECONDS_FCST                                          &  !<-- Length of forecast in seconds
                ,TIMESTEP_SEC_WHOLE                                     &  !<-- Integer part of timestep
                ,TIMESTEP_SEC_NUMERATOR                                 &  !<-- Numerator of fractional part
                ,TIMESTEP_SEC_DENOMINATOR                               &  !<-- Denominator of fractional part
                ,YY,MM,DD                                               &  !<-- Time variables for date
                ,HH,MNS,SEC                                                !<-- Time variables for time of day
!

      TYPE(NEMS_Rusage) :: rusage                                          !<-- Resource usage tracking object

      TYPE(ESMF_TimeInterval) :: RUNDURATION                            &  !<-- The ESMF time. The total forecast hours.
                                ,TIMESTEP                                  !<-- The ESMF timestep length (we only need a dummy here)
!
      TYPE(ESMF_Time) :: CURRTIME                                       &  !<-- The ESMF current time.
                        ,STARTTIME                                         !<-- The ESMF start time.
!
      TYPE(ESMF_VM) :: VM                                                  !<-- The ESMF virtual machine,
                                                                           !    which contains and manages
                                                                           !    the computer CPU resource
                                                                           !    for the ESMF grid components.
!
      TYPE(ESMF_GridComp) :: NEMS_GRID_COMP                                !<-- The NEMS gridded component.
!
      TYPE(ESMF_State) :: NEMS_EXP_STATE                                &  !<-- The NEMS export state
                         ,NEMS_IMP_STATE                                   !<-- The NEMS import state
!      
      TYPE(ESMF_Clock) :: CLOCK_MAIN                                       !<-- The ESMF time management clock
!
      TYPE(ESMF_Config) :: CF_MAIN                                         !<-- The Configure object
!
      LOGICAL :: RUN_CONTINUE                                           &  !<-- Flag for more than one NEMS run.
                ,PRINT_ESMF                                                !<-- Flag for ESMF PET files
!
      INTEGER :: HH_START                                               &
                ,HH_FINAL
!
      INTEGER :: RC, RC_USER                                               !<-- The running error signal
      INTEGER :: RUSAGE_RC                                                 !<-- Resource usage collection flag
!
      CHARACTER(LEN=MPI_MAX_PROCESSOR_NAME) :: PROCNAME                    !<-- The processor(host) name
      INTEGER :: PROCNAME_LEN                                              !<-- Actual PROCRNAME string length
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Check if we want ESMF PET files or not
!-----------------------------------------------------------------------
!
      CALL CHECK_ESMF_PET(PRINT_ESMF)
!
!-----------------------------------------------------------------------
!***  Initialize the ESMF framework. 
!-----------------------------------------------------------------------
! 
      IF(PRINT_ESMF) THEN
        CALL ESMF_Initialize(VM             =VM                         & !<-- The ESMF Virtual Machine
                            ,defaultCalKind =ESMF_CALKIND_GREGORIAN     & !<-- Set up the default calendar.
                            ,logkindflag    =ESMF_LOGKIND_MULTI         & !<-- Define multiple log error output files;
                            ,rc             =RC)
        ESMF_ERR_ABORT(RC)
      ELSE
        CALL ESMF_Initialize(VM             =VM                         & !<-- The ESMF Virtual Machine
                            ,defaultCalKind =ESMF_CALKIND_GREGORIAN     & !<-- Set up the default calendar.
                            ,logkindflag    =ESMF_LOGKIND_NONE          & !<-- Define no log error output files;
                            ,rc             =RC)
        ESMF_ERR_ABORT(RC)
      ENDIF
!
!-----------------------------------------------------------------------
!***  Get the processor(host) name
!-----------------------------------------------------------------------
!
      CALL MPI_Get_processor_name(PROCNAME,PROCNAME_LEN,RC)
!     write(0,*)'processor_name =',trim(PROCNAME)
!
!-----------------------------------------------------------------------
!***  Extract the MPI task ID.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="MAIN: Obtain the local task ID"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_VMGet(vm      =VM                                       &  !<-- The ESMF Virtual Machine
                     ,localpet=MYPE                                     &  !<-- The local MPI task ID
                     ,rc      =RC)
      ESMF_ERR_ABORT(RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Print subversion version and other status information.
!-----------------------------------------------------------------------
!
!#if defined (SVN_INFO) && defined (CMP_YEAR) && defined (CMP_JD)
!      if (mype==0) call w3tagb('NEMS '//SVN_INFO,                       &
!                               CMP_YEAR, CMP_JD, 0000, 'NEMS')
!#else
!      if (mype==0) call w3tagb('nems     ',0000,0000,0000,'np23   ')
!#endif
!
!-----------------------------------------------------------------------
!***  Start gathering resource usage information.
!-----------------------------------------------------------------------
!
      ! Note intentional use of MPI_COMM_WORLD since we are getting
      ! resource usage for ALL ranks on ALL machines.
      call rusage%start(MPI_COMM_WORLD,PROCNAME,PROCNAME_LEN,RC)
      ! It is safe to ignore RC since rusage%is_valid will tell us if
      ! the start succeeded.
!
!-----------------------------------------------------------------------
!***  Set up the default log.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Set Up ESMF Log"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      IF(PRINT_ESMF) THEN
        CALL ESMF_LogSet(flush      =.true.                             &
                        ,trace      =.false.                            &
! --> do not abort inside of ESMF, or else no ESMF backtrace will be in Log!!!!
!                        ,logmsgAbort=(/ ESMF_LOGMSG_ERROR /)            &
                        ,rc         =RC)
        ESMF_ERR_ABORT(RC)
      ENDIF

! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Create and load the Configure object which will hold the contents
!***  of the Main configure file.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create/Load the Main Configure Object"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CF_MAIN=ESMF_ConfigCreate(rc=RC)
!
      CALL ESMF_ConfigLoadFile(config  =CF_MAIN                         &  !<-- The Configure object
                              ,filename='model_configure'               &  !<-- The name of the configure file 
                              ,rc      =RC)
      ESMF_ERR_ABORT(RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Create the NEMS gridded component which will create and
!***  control the ATM (atmoshpere), OCN (ocean), ICE (sea ice), etc.
!***  gridded components.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create the NEMS Gridded Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      NEMS_GRID_COMP=ESMF_GridCompCreate(name        ='NEMS Grid Comp'  &  !<-- NEMS component name
                                        ,configFile  ='model_configure' &  !<-- Link the user-created configure file.
                                        ,rc          =RC)
      ESMF_ERR_ABORT(RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Register the NEMS gridded component's Initialize, Run and
!***  Finalize routines.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Register NEMS Gridded Component Init, Run, Finalize"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompSetServices(NEMS_GRID_COMP                      &  !<-- The NEMS component
                                   ,NEMS_REGISTER                       &  !<-- User's subroutineName
                                   ,rc=RC)
      ESMF_ERR_ABORT(RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Create the main ESMF Clock.
!***  The Clock is needed for all calls to Init, Run, and Finalize.
!
!***  A timestep is needed to create the Clock but actual timesteps
!***  are handled by the individual subcomponents therefore a dummy
!***  value is used here.
!-----------------------------------------------------------------------
!
      TIMESTEP_SEC_WHOLE      =1                                           !<-- Dummy timestep values
      TIMESTEP_SEC_NUMERATOR  =0                                           !
      TIMESTEP_SEC_DENOMINATOR=1                                           !<--
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Set up Time Step Interval in Main Clock"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_TimeIntervalSet(timeinterval=TIMESTEP                   &  !<-- Main Clock's timestep
                               ,s           =TIMESTEP_SEC_WHOLE         &  !<-- Whole part of timestep
                               ,sn          =TIMESTEP_SEC_NUMERATOR     &  !<-- Numerator of fractional part
                               ,sd          =TIMESTEP_SEC_DENOMINATOR   &  !<-- Denominator of fractional part
                               ,rc          =RC)
      ESMF_ERR_ABORT(RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Set the start time in the Main Clock.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="MAIN: Extract Starting Year from Config File"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF_MAIN                       &
                                  ,value =YY                            &
                                  ,label ='start_year:'                 &
                                  ,rc    =RC)
      ESMF_ERR_ABORT(RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="MAIN: Extract Starting Month from Config File"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF_MAIN                       &
                                  ,value =MM                            &
                                  ,label ='start_month:'                &
                                  ,rc    =RC)
      ESMF_ERR_ABORT(RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="MAIN: Extract Starting Day from Config File"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF_MAIN                       &
                                  ,value =DD                            &
                                  ,label ='start_day:'                  &
                                  ,rc    =RC)
      ESMF_ERR_ABORT(RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="MAIN: Extract Starting Hour from Config File"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF_MAIN                       &
                                  ,value =HH                            &
                                  ,label ='start_hour:'                 &
                                  ,rc    =RC)
      ESMF_ERR_ABORT(RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="MAIN: Extract Starting Minute from Config File"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF_MAIN                       &
                                  ,value =MNS                           &
                                  ,label ='start_minute:'               &
                                  ,rc    =RC)
      ESMF_ERR_ABORT(RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="MAIN: Extract Starting Second from Config File"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF_MAIN                       &
                                  ,value =SEC                           &
                                  ,label ='start_second:'               &
                                  ,rc    =RC)
      ESMF_ERR_ABORT(RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="MAIN: Set the Forecast Start Time"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_TimeSet(time=STARTTIME                                  &  !<-- The start time of the forecast (ESMF)
                       ,yy  =YY                                         &  !<-- Year from config file
                       ,mm  =MM                                         &  !<-- Month from config file
                       ,dd  =DD                                         &  !<-- Day from config file
                       ,h   =HH                                         &  !<-- Hour from config file
                       ,m   =MNS                                        &  !<-- Minute from config file
                       ,s   =SEC                                        &  !<-- Second from config file
                       ,rc  =RC)
      ESMF_ERR_ABORT(RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Set the run duration in the Main Clock.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="MAIN: Extract Forecast Length from Config File"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF_MAIN                       &
                                  ,value =NHOURS_FCST                   &
                                  ,label ='nhours_fcst:'                &
                                  ,rc    =RC)
      ESMF_ERR_ABORT(RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      NSECONDS_FCST=NHOURS_FCST*3600                                       !<-- The forecast length (sec) (REAL)
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="MAIN: Set the Forecast Length"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_TimeIntervalSet(timeinterval=RUNDURATION                &  !<-- The forecast length (s) (ESMF)
                               ,s           =NSECONDS_FCST              &
                               ,rc          =RC)
      ESMF_ERR_ABORT(RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
#ifdef FRONT_FVCOM
      SYEAR   = YY
      SMONTH  = MM
      SDAY    = DD
      SHOUR   = HH
      SMINUTE = MNS
      SSECOND = SEC
#endif
!-----------------------------------------------------------------------
!***  Now the Main Clock can be created.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create the Main Clock"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!

! NUOPC requires a correct timeStep to function. Here TIMESTEP was just created
! pro forma, but now is replaced with RUNDURATION to become meaningful.

      TIMESTEP = RUNDURATION

      CLOCK_MAIN=ESMF_ClockCreate(name       ='CLOCK_MAIN'              &  !<-- The top-level ESMF Clock
                                 ,timeStep   =TIMESTEP                  &  !<-- Timestep needed by the Clock (ESMF)
                                 ,startTime  =STARTTIME                 &  !<-- The integration start time (ESMF)
                                 ,runDuration=RUNDURATION               &  !<-- The integration duration (ESMF)
                                 ,rc         =RC)
      ESMF_ERR_ABORT(RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Create the NEMS component's import/export states.
!***  Currently they are not required to perform an actual function.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create the NEMS Import/Export States"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      NEMS_IMP_STATE=ESMF_StateCreate(name='NEMS Import State'     &
                                     ,rc       =RC)
      ESMF_ERR_ABORT(RC)
!
      NEMS_EXP_STATE=ESMF_StateCreate(name='NEMS Export State'     &
                                     ,rc       =RC)
      ESMF_ERR_ABORT(RC)
!      
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Execute the INITIALIZE step for the NEMS component.
!***  The Initialize routine that is called here as well as the
!***  Run and Finalize routines invoked below are those specified
!***  in the Register routine called in ESMF_GridCompSetServices above.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Execute the NEMS Component Initialize Step"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompInitialize(gridcomp   =NEMS_GRID_COMP           &  !<-- The NEMS component
                                  ,importState=NEMS_IMP_STATE           &  !<-- The NEMS import state
                                  ,exportState=NEMS_EXP_STATE           &  !<-- The NEMS export state
                                  ,clock      =CLOCK_MAIN               &  !<-- The ESMF clock
                                  ,phase      =1                        &
                                  ,userRc     =RC_USER                  &
                                  ,rc         =RC)
      ESMF_ERR_ABORT(RC)
      ESMF_ERR_ABORT(RC_USER)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Execute the RUN step for the NEMS component.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Execute the NEMS Component Run Step"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompRun(gridcomp   =NEMS_GRID_COMP                  &  !<-- The NEMS component
                           ,importState=NEMS_IMP_STATE                  &  !<-- The NEMS import state
                           ,exportState=NEMS_EXP_STATE                  &  !<-- The NEMS export state
                           ,clock      =CLOCK_MAIN                      &  !<-- The ESMF clock
                           ,phase      =1                               &
                           ,userRc     =RC_USER                         &
                           ,rc         =RC)
      ESMF_ERR_ABORT(RC)
      ESMF_ERR_ABORT(RC_USER)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  The RUN_CONTINUE flag tells us if the RUN step of the
!***  NEMS component must be called multiple times for ensembles.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract the RUN_CONTINUE flag from Config File"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config = CF_MAIN                     &
                                  ,value  = RUN_CONTINUE                &
                                  ,label  = 'RUN_CONTINUE:'             &
                                  ,default= .false.                     &
                                  ,rc     = RC)
      ESMF_ERR_ABORT(RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Update the Main clock.  This is for calling the RUN step
!***  of the NEMS component multiple times.
!-----------------------------------------------------------------------
!
      IF(RUN_CONTINUE) THEN
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract the Ensemble Clock Parameters from Config File"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ConfigGetAttribute(config = CF_MAIN                   &
                                    ,value  = hh_start                  &
                                    ,label  = 'HH_START:'               &
                                    ,rc     = RC)
        ESMF_ERR_ABORT(RC)
!
        CALL ESMF_ConfigGetAttribute(config = CF_MAIN                   &
                                    ,value  = hh_final                  &
                                    ,label  = 'HH_FINAL:'               &
                                    ,rc     = RC)
        ESMF_ERR_ABORT(RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        NHOURS_FCST = HH_FINAL - HH_START
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="MAIN: Re-set the clock after the ensemble run cycles."
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_TimeIntervalSet(timeInterval = RUNDURATION            &
                                 ,h            = NHOURS_FCST            &
                                 ,rc           = RC)
        ESMF_ERR_ABORT(RC)
!
        CALL ESMF_ClockGet(clock    = CLOCK_MAIN                        &
                          ,currTime = CURRTIME                          &
                          ,rc       = rc)
        ESMF_ERR_ABORT(RC)
!
        CURRTIME = CURRTIME + RUNDURATION
!
        CALL ESMF_ClockSet(clock       = CLOCK_MAIN                     &
                          ,RunDuration = RUNDURATION                    &
                          ,currTime    = CURRTIME                       &
                          ,rc          = RC)
        ESMF_ERR_ABORT(RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      END IF
!
!-----------------------------------------------------------------------
!***  Execute the FINALIZE step for the NEMS component.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Execute the NEMS Component Finalize Step"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompFinalize(gridcomp   =NEMS_GRID_COMP             &  !<-- The NEMS component
                                ,importState=NEMS_IMP_STATE             &  !<-- The NEMS component import state
                                ,exportState=NEMS_EXP_STATE             &  !<-- The NEMS component export state
                                ,clock      =CLOCK_MAIN                 &  !<-- The Main ESMF clock
                                ,phase      =1                          &
                                ,userRc     =RC_USER                    &
                                ,rc         =RC)
      ESMF_ERR_ABORT(RC)
      ESMF_ERR_ABORT(RC_USER)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Destroy the Main Clock.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Destroy the Main Clock"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ClockDestroy(clock=CLOCK_MAIN                           &
                            ,rc   =RC)
      ESMF_ERR_ABORT(RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Destroy the Main Configure object.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Destroy the Main Configure Object"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigDestroy(config=CF_MAIN                            &
                             ,rc    =RC)
      ESMF_ERR_ABORT(RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK = "Destroy the ESMF states"
!      CALL ESMF_LogWrite(MESSAGE_CHECK, ESMF_LOGMSG_INFO, rc = RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

      CALL ESMF_StateDestroy(state=NEMS_IMP_STATE                       &
                            ,rc   =RC)
      ESMF_ERR_ABORT(RC)
!
      CALL ESMF_StateDestroy(state=NEMS_EXP_STATE                       &
                            ,rc   =RC)
      ESMF_ERR_ABORT(RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK = "Destroy ESMF Grid Comp and Cpl Comp"
!      CALL ESMF_LogWrite(MESSAGE_CHECK, ESMF_LOGMSG_INFO, rc = RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompDestroy(gridcomp=NEMS_GRID_COMP                 &
                               ,rc      =RC)
      ESMF_ERR_ABORT(RC)

! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
! Resource usage reporting.  We only do this if the rusage%start was
! successful.
      if(rusage%is_valid()) then
         call rusage%stop(rc)
         if(rc==0) call rusage%report(rc)
      endif

!-----------------------------------------------------------------------
!***  Shut down the ESMF system.
!-----------------------------------------------------------------------
!
      CALL ESMF_Finalize()
!
!-----------------------------------------------------------------------
!
!      if(mype==0)call w3tage('nems     ')
!
!-----------------------------------------------------------------------
!
      END PROGRAM MAIN_NEMS
!
!-----------------------------------------------------------------------


#ifndef IBM
        REAL(8) FUNCTION RTC()
          RTC = 0.d0
        END FUNCTION

        REAL(8) FUNCTION TIMEF()
          TIMEF = 0.d0
        END FUNCTION
#endif

