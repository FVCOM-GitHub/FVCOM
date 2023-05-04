#include "./ESMFVersionDefine.h"

!-----------------------------------------------------------------------
!
      MODULE module_EARTH_GRID_COMP
!
!-----------------------------------------------------------------------
!***  This module contains codes directly related to the EARTH component.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!  2010-03-24  Black - Created Earth component module.
!  2010-04     Yang  - Added Ensemble capability.
!  2011-05-11  Theurich & Yang - Modified for using the ESMF 5.2.0r_beta_snapshot_07.
!  2011-10-04  Yang - Modified for using the ESMF 5.2.0r library.
!  2012-02     Tripp - Added ESMF superstructure to support an OCN model
!  2013-06     Theurich - Reworked OCN dependency to be NUOPC based
!  2013-07     Theurich - Macro based ESMF error handling
!-----------------------------------------------------------------------
!
!***  The EARTH component lies in the hierarchy seen here:
!
!          Main program
!               |
!               |
!          NEMS component
!               |     |________________________.
!               |                              |
!          EARTH component        Ensemble Coupler component
!              /|\
!             / | \
!          ATM/OCN/ICE/WAV/LND/IPM/HYD .. components
!          |    |   |
!          |    |   (CICE, etc.)
!          |    |
!          |    (MOM5, MOM6, HYCOM, POM, etc.)
!          |
!          CORE component (GSM, NMM, FV3, etc.)
!
!-----------------------------------------------------------------------
!
      USE ESMF

      use NUOPC
      use NUOPC_Driver, &
        Driver_routine_SS             => SetServices, &
        Driver_label_SetModelServices => label_SetModelServices, &
        Driver_label_SetRunSequence   => label_SetRunSequence, &
        Driver_label_SetRunClock      => label_SetRunClock, &
        Driver_label_Finalize         => label_Finalize
      use NUOPC_Connector, only: conSS => SetServices
  ! - Handle build time ATM options:
#ifdef FRONT_SATM
      use FRONT_SATM,       only: SATM_SS   => SetServices
#endif
#ifdef FRONT_XATM
      use FRONT_XATM,       only: XATM_SS   => SetServices
#endif
#ifdef FRONT_DATAWAM
      use FRONT_DATAWAM,    only: DATAWAM_SS=> SetServices
#endif
#ifdef FRONT_GSM
      use FRONT_GSM,        only: GSM_SS   => SetServices
#endif
#ifdef FRONT_NMMB
      use FRONT_NMMB,       only: NMMB_SS   => SetServices
#endif
#ifdef FRONT_FV3
      use FRONT_FV3,        only: FV3_SS   => SetServices
#endif
#ifdef FRONT_DATM
      use FRONT_DATM,       only: DATM_SS  => SetServices
#endif
  ! - Handle build time OCN options:
#ifdef FRONT_SOCN
      use FRONT_SOCN,       only: SOCN_SS   => SetServices
#endif
#ifdef FRONT_XOCN
      use FRONT_XOCN,       only: XOCN_SS   => SetServices
#endif
#ifdef FRONT_HYCOM
      use FRONT_HYCOM,      only: HYCOM_SS  => SetServices
#endif
#ifdef FRONT_MOM5
      use FRONT_MOM5,       only: MOM5_SS   => SetServices
#endif
#ifdef FRONT_MOM6
      use FRONT_MOM6,       only: MOM6_SS   => SetServices
#endif
#ifdef FRONT_POM
      use FRONT_POM,        only: POM_SS    => SetServices
#endif
#ifdef FRONT_SCHISM
      use FRONT_SCHISM,     only: SCHISM_SS  => SetServices
#endif
#ifdef FRONT_ADCIRC
      use FRONT_ADCIRC,     only: ADCIRC_SS  => SetServices
#endif
#ifdef FRONT_HWRFDATA
      use FRONT_HWRFDATA,   only: HWRFDATA_SS  => SetServices
#endif
#ifdef FRONT_WW3DATA
      use FRONT_WW3DATA,    only: WW3DATA_SS  => SetServices
#endif
#ifdef FRONT_BARDATA
      use FRONT_BARDATA,     only: BARDATA_SS  => SetServices
#endif
#ifdef FRONT_ATMESH
      use FRONT_ATMESH,     only: ATMESH_SS  => SetServices
#endif
#ifdef FRONT_PAHM
      use FRONT_PAHM,   only: PAHM_SS  => SetServices
#endif
  ! - Handle build time ICE options:
#ifdef FRONT_SICE
      use FRONT_SICE,       only: SICE_SS  => SetServices
#endif
#ifdef FRONT_XICE
      use FRONT_XICE,       only: XICE_SS  => SetServices
#endif
#ifdef FRONT_CICE
      use FRONT_CICE,       only: CICE_SS  => SetServices
#endif
#ifdef FRONT_CICE6
      use FRONT_CICE6,      only: CICE6_SS => SetServices
#endif
  ! - Handle build time WAV options:
#ifdef FRONT_SWAV
      use FRONT_SWAV,       only: SWAV_SS  => SetServices
#endif
#ifdef FRONT_XWAV
      use FRONT_XWAV,       only: XWAV_SS  => SetServices
#endif
#ifdef FRONT_WW3
      use FRONT_WW3,        only: WW3_SS  => SetServices
#endif
  ! - Handle build time LND options:
#ifdef FRONT_SLND
      use FRONT_SLND,       only: SLND_SS  => SetServices
#endif
#ifdef FRONT_XLND
      use FRONT_XLND,       only: XLND_SS  => SetServices
#endif
#ifdef FRONT_NOAH
      use FRONT_NOAH,       only: NOAH_SS  => SetServices
#endif
#ifdef FRONT_LIS
      use FRONT_LIS,        only: LIS_SS   => SetServices
#endif
  ! - Handle build time IPM options:
#ifdef FRONT_SIPM
      use FRONT_SIPM,       only: SIPM_SS  => SetServices
#endif
#ifdef FRONT_XIPM
      use FRONT_XIPM,       only: XIPM_SS  => SetServices
#endif
#ifdef FRONT_IPE
      use FRONT_IPE,        only: IPE_SS   => SetServices
#endif
#ifdef FRONT_DATAIPE
      use FRONT_DATAIPE,    only: DATAIPE_SS=> SetServices
#endif
  ! - Handle build time HYD options:
#ifdef FRONT_SHYD
      use FRONT_SHYD,       only: SHYD_SS  => SetServices
#endif
#ifdef FRONT_XHYD
      use FRONT_XHYD,       only: XHYD_SS  => SetServices
#endif
#ifdef FRONT_WRFHYDRO
      use FRONT_WRFHYDRO,   only: WRFHYDRO_SS  => SetServices
#endif
#ifdef FRONT_GSDCHEM
      use FRONT_GSDCHEM,    only: GSDCHEM_SS  => SetServices
#endif
#ifdef FRONT_NWM
      use FRONT_NWM,        only: NWM_SS  => SetServices
#endif
#ifdef FRONT_FVCOM
      use FVCOM_CAP,        only: FVCOMSS => SetServices
      use FVCOM_CAP,        only: fvcom_name
      use mod_driver,       only: fvcom_pet_num
#endif
  ! - Mediator
#ifdef CMEPS
      use MED, only : MED_SS     => SetServices 
#else
      use module_MEDIATOR,        only: MED_SS     => SetServices
#endif
      use module_MEDSpaceWeather, only: MEDSW_SS   => SetServices

      USE module_EARTH_INTERNAL_STATE,ONLY: EARTH_INTERNAL_STATE        &
                                           ,WRAP_EARTH_INTERNAL_STATE
!
!      USE module_ATM_GRID_COMP
!
      USE module_NEMS_UTILS,ONLY: ERR_MSG,MESSAGE_CHECK
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
      PRIVATE
!
      PUBLIC :: EARTH_REGISTER
      PUBLIC :: VERBOSE_DIAGNOSTICS
!
!-----------------------------------------------------------------------
!

      LOGICAL, PRIVATE :: flag_verbose_diagnostics = .false.


      CONTAINS

      logical function verbose_diagnostics(set)
        !! Mutator for the verbose diagnostics flag; returns true if
        !! verbose diagnostics should be used, and false otherwise.
        !! If the "set" argument is present, then the flag is set to
        !! the given value.
        logical, optional :: set
        if(present(set)) then
           flag_verbose_diagnostics=set
        endif
        verbose_diagnostics=flag_verbose_diagnostics
      end function verbose_diagnostics

!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE EARTH_REGISTER(EARTH_GRID_COMP,RC_REG)
!
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      TYPE(ESMF_GridComp) :: EARTH_GRID_COMP                               !<-- The EARTH component
!
      INTEGER,INTENT(OUT) :: RC_REG                                        !<-- Error return code
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER :: RC
      type(ESMF_Config)             :: config
      
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC_REG = ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      ! Derive from NUOPC_Driver
      call NUOPC_CompDerive(EARTH_GRID_COMP, Driver_routine_SS, rc=RC)
      ESMF_ERR_RETURN(RC,RC_REG)

      ! specializations:

      call NUOPC_CompSpecialize(EARTH_GRID_COMP, &
        specLabel=Driver_label_SetModelServices, specRoutine=SetModelServices, &
        rc=RC)
      ESMF_ERR_RETURN(RC,RC_REG)
      
      call NUOPC_CompSpecialize(EARTH_GRID_COMP, &
        specLabel=Driver_label_SetRunSequence, specRoutine=SetRunSequence, &
        rc=RC)
      ESMF_ERR_RETURN(RC,RC_REG)

      ! The NEMS Earth component is currently the top-level driver and
      ! does not need to coordinate Clocks with its parent.
      call ESMF_MethodRemove(EARTH_GRID_COMP, Driver_label_SetRunClock, rc=RC_REG)
      ESMF_ERR_RETURN(RC,RC_REG)
      call NUOPC_CompSpecialize(EARTH_GRID_COMP, &
        specLabel=Driver_label_SetRunClock, specRoutine=NUOPC_NoOp, rc=RC_REG)
      ESMF_ERR_RETURN(RC,RC_REG)
      
      call NUOPC_CompSpecialize(EARTH_GRID_COMP, &
        specLabel=Driver_label_Finalize, specRoutine=Finalize, &
        rc=RC)
      ESMF_ERR_RETURN(RC,RC_REG)
      
      ! register an internal initialization method
      call NUOPC_CompSetInternalEntryPoint(EARTH_GRID_COMP, ESMF_METHOD_INITIALIZE, &
        phaseLabelList=(/"IPDv04p2"/), userRoutine=ModifyCplLists, rc=rc)
      ESMF_ERR_RETURN(RC,RC_REG)

      ! create, open, and set the config
      config = ESMF_ConfigCreate(rc=RC)
      ESMF_ERR_RETURN(RC,RC_REG)
      call ESMF_ConfigLoadFile(config, "nems.configure", rc=RC)
      ESMF_ERR_RETURN(RC,RC_REG)
      call ESMF_GridCompSet(EARTH_GRID_COMP, config=config, rc=RC)
      ESMF_ERR_RETURN(RC,RC_REG)
      
      ! Added the following Field Dictionary block to the EARTH component level
      ! in order to prevent different dictionary definitions in the lower
      ! components. Doing this here isn't without problems because it
      ! potentially makes the components (ATM & OCN) depend on this environment,
      ! which lowers their transferability to other coupled systems. However,
      ! extending the Field Dictionary is a temporary solution anyway (see the
      ! TODO: below), so this isn't going to stay for ever this way.
      
      ! Extend the NUOPC Field Dictionary to hold required entries.
      !TODO: In the long run this section will not be needed when we have
      !TODO: absorbed the needed standard names into the default dictionary.
      ! -> 20 fields identified as exports by the GSM component
#ifdef CMEPS
      call NUOPC_FieldDictionarySetup("fd_nems.yaml", rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
#else
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "air_density_height_lowest")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="air_density_height_lowest", &
          canonicalUnits="kg m-3", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "mean_zonal_moment_flx")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_zonal_moment_flx", &
          canonicalUnits="N m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "mean_zonal_moment_flx_atm")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_zonal_moment_flx_atm", &
          canonicalUnits="N m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "mean_merid_moment_flx")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_merid_moment_flx", &
          canonicalUnits="N m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "mean_merid_moment_flx_atm")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_merid_moment_flx_atm", &
          canonicalUnits="N m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "mean_sensi_heat_flx")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_sensi_heat_flx", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "mean_sensi_heat_flx_atm")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_sensi_heat_flx_atm", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "mean_sensi_heat_flx_atm_into_ice")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_sensi_heat_flx_atm_into_ice", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "mean_sensi_heat_flx_atm_into_ocn")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_sensi_heat_flx_atm_into_ocn", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "mean_laten_heat_flx")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_laten_heat_flx", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "mean_laten_heat_flx_atm")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_laten_heat_flx_atm", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "mean_laten_heat_flx_atm_into_ice")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_laten_heat_flx_atm_into_ice", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "mean_laten_heat_flx_atm_into_ocn")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_laten_heat_flx_atm_into_ocn", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "mean_down_lw_flx")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_down_lw_flx", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "mean_down_sw_flx")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_down_sw_flx", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "mean_fprec_rate")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_fprec_rate", &
          canonicalUnits="kg s m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "mean_prec_rate")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_prec_rate", &
          canonicalUnits="kg s m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "mean_evap_rate")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_evap_rate", &
          canonicalUnits="kg s m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "mean_evap_rate_atm_into_ice")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_evap_rate_atm_into_ice", &
          canonicalUnits="kg s m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "mean_evap_rate_atm_into_ocn")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_evap_rate_atm_into_ocn", &
          canonicalUnits="kg s m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "inst_zonal_moment_flx")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_zonal_moment_flx", &
          canonicalUnits="N m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "inst_merid_moment_flx")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_merid_moment_flx", &
          canonicalUnits="N m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "inst_sensi_heat_flx")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_sensi_heat_flx", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "inst_laten_heat_flx")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_laten_heat_flx", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "inst_down_lw_flx")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_down_lw_flx", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "inst_down_sw_flx")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_down_sw_flx", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "inst_temp_height2m")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_temp_height2m", &
          canonicalUnits="K", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "inst_spec_humid_height2m")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_spec_humid_height2m", &
          canonicalUnits="kg kg-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "inst_u_wind_height10m")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_u_wind_height10m", &
          canonicalUnits="m s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "inst_v_wind_height10m")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_v_wind_height10m", &
          canonicalUnits="m s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "inst_zonal_wind_height10m")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_zonal_wind_height10m", &
          canonicalUnits="m s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "inst_merid_wind_height10m")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_merid_wind_height10m", &
          canonicalUnits="m s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif

      !For MOM6 and WW3 variables to match:
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "eastward_partitioned_stokes_drift_1")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="eastward_partitioned_stokes_drift_1", &
          canonicalUnits="m s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif 

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "northward_partitioned_stokes_drift_1")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="northward_partitioned_stokes_drift_1", &
          canonicalUnits="m s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "eastward_partitioned_stokes_drift_2")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="eastward_partitioned_stokes_drift_2", &
          canonicalUnits="m s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "northward_partitioned_stokes_drift_2")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="northward_partitioned_stokes_drift_2", &
          canonicalUnits="m s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif



     if (.not.NUOPC_FieldDictionaryHasEntry( &
        "eastward_partitioned_stokes_drift_3")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="eastward_partitioned_stokes_drift_3", &
          canonicalUnits="m s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "northward_partitioned_stokes_drift_3")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="northward_partitioned_stokes_drift_3", &
          canonicalUnits="m s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      ! end of MOM6 and WW3 variables to match

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "inst_temp_height_surface")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_temp_height_surface", &
          canonicalUnits="K", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "inst_pres_height_surface")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_pres_height_surface", &
          canonicalUnits="Pa", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "inst_surface_height")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_surface_height", &
          canonicalUnits="m", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      ! -> Additional fields identified as needed by MOM5 and others...
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "mean_down_sw_vis_dir_flx")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_down_sw_vis_dir_flx", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "mean_down_sw_vis_dif_flx")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_down_sw_vis_dif_flx", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "mean_down_sw_ir_dir_flx")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_down_sw_ir_dir_flx", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "mean_down_sw_ir_dif_flx")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_down_sw_ir_dif_flx", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "inst_down_sw_vis_dir_flx")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_down_sw_vis_dir_flx", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "inst_down_sw_vis_dif_flx")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_down_sw_vis_dif_flx", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "inst_down_sw_ir_dir_flx")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_down_sw_ir_dir_flx", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "inst_down_sw_ir_dif_flx")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_down_sw_ir_dif_flx", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "mean_net_sw_vis_dir_flx")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_net_sw_vis_dir_flx", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "mean_net_sw_vis_dif_flx")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_net_sw_vis_dif_flx", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "mean_net_sw_ir_dir_flx")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_net_sw_ir_dir_flx", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "mean_net_sw_ir_dif_flx")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_net_sw_ir_dif_flx", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "inst_net_sw_vis_dir_flx")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_net_sw_vis_dir_flx", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "inst_net_sw_vis_dif_flx")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_net_sw_vis_dif_flx", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "inst_net_sw_ir_dir_flx")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_net_sw_ir_dir_flx", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "inst_net_sw_ir_dif_flx")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_net_sw_ir_dif_flx", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "mean_salt_rate")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_salt_rate", &
          canonicalUnits="kg psu m-2 s", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "mean_runoff_rate")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_runoff_rate", &
          canonicalUnits="kg m-2 s", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "mean_calving_rate")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_calving_rate", &
          canonicalUnits="kg m-2 s", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "mean_runoff_heat_flx")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_runoff_heat_flx", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry(  &
        "mean_calving_heat_flx")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_calving_heat_flx", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "ice_fraction")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="ice_fraction", &
          canonicalUnits="1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "mean_sw_pen_to_ocn")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_sw_pen_to_ocn", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "mean_up_lw_flx")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_up_lw_flx", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "mass_of_overlying_sea_ice")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mass_of_overlying_sea_ice", &
          canonicalUnits="kg", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "s_surf")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="s_surf", &
          canonicalUnits="psu", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "freezing_melting_potential")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="freezing_melting_potential", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      ! following two added for export from MOM6
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "accum_heat_frazil")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="accum_heat_frazil", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "inst_melt_potential")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_melt_potential", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "u_surf")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="u_surf", &
          canonicalUnits="m s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "v_surf")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="v_surf", &
          canonicalUnits="m s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "sea_lev")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="sea_lev", &
          canonicalUnits="m", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif 
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "wind_stress_zonal")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="wind_stress_zonal", &
          canonicalUnits="N m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif 
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "wind_stress_merid")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="wind_stress_merid", &
          canonicalUnits="N m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif 
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "ocn_current_zonal")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="ocn_current_zonal", &
          canonicalUnits="m s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif 
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "ocn_current_merid")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="ocn_current_merid", &
          canonicalUnits="m s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif 
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "ocn_current_idir")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="ocn_current_idir", &
          canonicalUnits="m s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif 
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "ocn_current_jdir")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="ocn_current_jdir", &
          canonicalUnits="m s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif 
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "sea_surface_slope_zonal")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="sea_surface_slope_zonal", &
          canonicalUnits="m m-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif 
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "sea_surface_slope_merid")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="sea_surface_slope_merid", &
          canonicalUnits="m m-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif 
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "sea_surface_slope_zonal")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="sea_surface_slope_zonal", &
          canonicalUnits="m m-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "sea_surface_slope_merid")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="sea_surface_slope_merid", &
          canonicalUnits="m m-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "stress_on_air_ice_zonal")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="stress_on_air_ice_zonal", &
          canonicalUnits="N m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif 
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "stress_on_air_ice_merid")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="stress_on_air_ice_merid", &
          canonicalUnits="N m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif 
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "stress_on_air_ocn_zonal")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="stress_on_air_ocn_zonal", &
          canonicalUnits="N m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif 
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "stress_on_air_ocn_merid")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="stress_on_air_ocn_merid", &
          canonicalUnits="N m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif 
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "stress_on_ocn_ice_zonal")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="stress_on_ocn_ice_zonal", &
          canonicalUnits="N m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif 
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "stress_on_ocn_ice_merid")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="stress_on_ocn_ice_merid", &
          canonicalUnits="N m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif 
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "stress_on_ocn_ice_idir")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="stress_on_ocn_ice_idir", &
          canonicalUnits="N m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif 
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "stress_on_ocn_ice_jdir")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="stress_on_ocn_ice_jdir", &
          canonicalUnits="N m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif 
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "mixed_layer_depth")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mixed_layer_depth", &
          canonicalUnits="m", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif 
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "mean_net_lw_flx")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_net_lw_flx", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif 
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "mean_net_lw_flx_atm")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_net_lw_flx_atm", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif 
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "mean_net_sw_flx")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_net_sw_flx", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif 
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "mean_up_lw_flx_ice")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_up_lw_flx_ice", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif 
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "mean_up_lw_flx_ocn")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_up_lw_flx_ocn", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif 
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "inst_net_lw_flx")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_net_lw_flx", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif 
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "inst_net_sw_flx")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_net_sw_flx", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "mean_sw_pen_to_ocn_vis_dir_flx")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_sw_pen_to_ocn_vis_dir_flx", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "mean_sw_pen_to_ocn_vis_dif_flx")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_sw_pen_to_ocn_vis_dif_flx", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "mean_sw_pen_to_ocn_ir_dir_flx")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_sw_pen_to_ocn_ir_dir_flx", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "mean_sw_pen_to_ocn_ir_dif_flx")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_sw_pen_to_ocn_ir_dif_flx", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "inst_ir_dir_albedo")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_ir_dir_albedo", &
          canonicalUnits="1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif 
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "inst_ir_dif_albedo")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_ir_dif_albedo", &
          canonicalUnits="1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif 
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "inst_vis_dir_albedo")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_vis_dir_albedo", &
          canonicalUnits="1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif 
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "inst_vis_dif_albedo")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_vis_dif_albedo", &
          canonicalUnits="1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif 
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "inst_ocn_ir_dir_albedo")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_ocn_ir_dir_albedo", &
          canonicalUnits="1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif 
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "inst_ocn_ir_dif_albedo")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_ocn_ir_dif_albedo", &
          canonicalUnits="1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif 
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "inst_ocn_vis_dir_albedo")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_ocn_vis_dir_albedo", &
          canonicalUnits="1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif 
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "inst_ocn_vis_dif_albedo")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_ocn_vis_dif_albedo", &
          canonicalUnits="1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "inst_ice_ir_dir_albedo")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_ice_ir_dir_albedo", &
          canonicalUnits="1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif 
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "inst_ice_ir_dif_albedo")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_ice_ir_dif_albedo", &
          canonicalUnits="1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif 
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "inst_ice_vis_dir_albedo")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_ice_vis_dir_albedo", &
          canonicalUnits="1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif 
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "inst_ice_vis_dif_albedo")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_ice_vis_dif_albedo", &
          canonicalUnits="1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "inst_land_sea_mask")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_land_sea_mask", &
          canonicalUnits="1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "inst_temp_height_lowest")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_temp_height_lowest", &
          canonicalUnits="K", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "inst_spec_humid_height_lowest")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_spec_humid_height_lowest", &
          canonicalUnits="kg kg-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "humidity_2m")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="humidity_2m", &
          canonicalUnits="kg kg-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "inst_zonal_wind_height_lowest")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_zonal_wind_height_lowest", &
          canonicalUnits="m s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "inst_merid_wind_height_lowest")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_merid_wind_height_lowest", &
          canonicalUnits="m s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "inst_pres_height_lowest")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_pres_height_lowest", &
          canonicalUnits="Pa", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "inst_height_lowest")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_height_lowest", &
          canonicalUnits="m", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "ocean_mask")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="ocean_mask", &
          canonicalUnits="1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "ice_mask")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="ice_mask", &
          canonicalUnits="1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "land_mask")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="land_mask", &
          canonicalUnits="1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      ! special HYCOM exports
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "surface_downward_eastward_stress")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="surface_downward_eastward_stress", &
          canonicalUnits="Pa", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "surface_downward_northward_stress")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="surface_downward_northward_stress", &
          canonicalUnits="Pa", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "wind_speed_height10m")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="wind_speed_height10m", &
          canonicalUnits="m s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "wind_speed_squared_10m")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="wind_speed_squared_10m", &
          canonicalUnits="m2 s-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "friction_speed")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="friction_speed", &
          canonicalUnits="m s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "mean_lat_flx")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_lat_flx", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "mean_sens_flx")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_sens_flx", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "water_flux_into_sea_water")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="water_flux_into_sea_water", &
          canonicalUnits="kg m-2 s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif 
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "frozen_water_flux_into_sea_water")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="frozen_water_flux_into_sea_water", &
          canonicalUnits="kg m-2 s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif 
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "surface_temperature")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="surface_temperature", &
          canonicalUnits="K", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "air_surface_temperature")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="air_surface_temperature", &
          canonicalUnits="K", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "temperature_2m")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="temperature_2m", &
          canonicalUnits="K", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "upward_sea_ice_basal_available_heat_flux")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="upward_sea_ice_basal_available_heat_flux", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      ! special HYCOM imports
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "sea_ice_area_fraction")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="sea_ice_area_fraction", &
          canonicalUnits="1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "downward_x_stress_at_sea_ice_base")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="downward_x_stress_at_sea_ice_base", &
          canonicalUnits="Pa", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "downward_y_stress_at_sea_ice_base")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="downward_y_stress_at_sea_ice_base", &
          canonicalUnits="Pa", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "downward_sea_ice_basal_solar_heat_flux")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="downward_sea_ice_basal_solar_heat_flux", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "upward_sea_ice_basal_heat_flux")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="upward_sea_ice_basal_heat_flux", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "downward_sea_ice_basal_salt_flux")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="downward_sea_ice_basal_salt_flux", &
          canonicalUnits="kg m-2 s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "downward_sea_ice_basal_water_flux")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="downward_sea_ice_basal_water_flux", &
          canonicalUnits="kg m-2 s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "sea_ice_surface_temperature")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="sea_ice_surface_temperature", &
          canonicalUnits="K", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "sea_ice_temperature")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="sea_ice_temperature", &
          canonicalUnits="K", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "sea_ice_thickness")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="sea_ice_thickness", &
          canonicalUnits="m", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "sea_ice_x_velocity")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="sea_ice_x_velocity", &
          canonicalUnits="m s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "sea_ice_y_velocity")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="sea_ice_y_velocity", &
          canonicalUnits="m s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif

      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "net_heat_flx_to_ocn")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="net_heat_flx_to_ocn", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "mean_fresh_water_to_ocean_rate")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_fresh_water_to_ocean_rate", &
          canonicalUnits="kg m-2 s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "mean_ice_volume")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_ice_volume", &
          canonicalUnits="m", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "mean_snow_volume")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_snow_volume", &
          canonicalUnits="m", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
     
      !Mass flux of liquid runoff
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "Foxx_rofl")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="Foxx_rofl", &
          canonicalUnits="kg m-2 s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      !Mass flux of frozen runoff
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "Foxx_rofi")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="Foxx_rofi", &
          canonicalUnits="kg m-2 s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      !Ocean surface boundary layer depth
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "So_bldepth")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="So_bldepth", &
          canonicalUnits="m", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
 
      ! Synonyms for HYCOM fields
      call NUOPC_FieldDictionarySetSyno( &
        standardNames = (/"surface_downward_eastward_stress",&
                          "mean_zonal_moment_flx           "/), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      call NUOPC_FieldDictionarySetSyno( &
        standardNames = (/"surface_downward_northward_stress",&
                          "mean_merid_moment_flx            "/), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      call NUOPC_FieldDictionarySetSyno( &
        standardNames = (/"mean_lat_flx       ",&
                          "mean_laten_heat_flx"/), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      call NUOPC_FieldDictionarySetSyno( &
        standardNames = (/"mean_sens_flx      ",&
                          "mean_sensi_heat_flx"/), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

      ! DCR - Fields added for Regional Application
      ! ATM-OCN-ICE-LND-HYD
      ! List of exisitng fields
      ! ice_mask, inst_down_lw_flx, inst_down_sw_flx, inst_height_lowest,
      ! inst_merid_wind_height_lowest, inst_pres_height_lowest,
      ! inst_pres_height_surface, inst_spec_humid_height_lowest,
      ! inst_temp_height_lowest, inst_temp_height_surface,
      ! inst_zonal_wind_height_lowest, mean_down_lw_flx, mean_down_sw_flx,
      ! mean_fprec_rate, mean_laten_heat_flx, mean_net_lw_flx, mean_net_sw_flx,
      ! mean_prec_rate, mean_sensi_heat_flx

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "aerodynamic_roughness_length")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="aerodynamic_roughness_length", &
          canonicalUnits="m", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "canopy_moisture_storage")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="canopy_moisture_storage", &
          canonicalUnits="kg m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "carbon_dioxide")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="carbon_dioxide", &
          canonicalUnits="ppmv", & ! Units must be clarified
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "cosine_zenith_angle")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="cosine_zenith_angle", &
          canonicalUnits="degree", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "exchange_coefficient_heat")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="exchange_coefficient_heat", &
          canonicalUnits="W m-2 K-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "exchange_coefficient_heat_height2m")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="exchange_coefficient_heat_height2m", &
          canonicalUnits="W m-2 K-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "exchange_coefficient_moisture_height2m")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="exchange_coefficient_moisture_height2m", &
          canonicalUnits="kg m-2 s-1 Pa-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "inst_wind_speed_height_lowest")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_wind_speed_height_lowest", &
          canonicalUnits="m s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "mean_cprec_rate")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_cprec_rate", &
          canonicalUnits="kg m-2 s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "mean_grnd_sensi_heat_flx")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_grnd_sensi_heat_flx", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "mean_laten_heat_flx_kinematic")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_laten_heat_flx_kinematic", &
          canonicalUnits="Kg m-2 s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "mean_surface_albedo")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_surface_albedo", &
          canonicalUnits="1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "mean_surface_skin_temp")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_surface_skin_temp", &
          canonicalUnits="K", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "mixing_ratio_surface")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mixing_ratio_surface", &
          canonicalUnits="kg kg-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "root_moisture")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="root_moisture", &
          canonicalUnits="kg m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "saturated_mixing_ratio")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="saturated_mixing_ratio", &
          canonicalUnits="kg kg-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "surface_snow_area_fraction")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="surface_snow_area_fraction", &
          canonicalUnits="1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "surface_snow_thickness")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="surface_snow_thickness", &
          canonicalUnits="m", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "surface_snow_melt_flux")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="surface_snow_melt_flux", &
          canonicalUnits="kg m-2 s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "liquid_water_content_of_surface_snow")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="liquid_water_content_of_surface_snow", &
          canonicalUnits="kg m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "soil_depth")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="soil_depth", &
          canonicalUnits="m", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "soil_hydraulic_conductivity_at_saturation")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="soil_hydraulic_conductivity_at_saturation", &
          canonicalUnits="m s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "moisture_content_of_soil_layer")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="moisture_content_of_soil_layer", &
          canonicalUnits="kg m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "moisture_content_of_soil_layer_1")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="moisture_content_of_soil_layer_1", &
          canonicalUnits="kg m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "moisture_content_of_soil_layer_2")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="moisture_content_of_soil_layer_2", &
          canonicalUnits="kg m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "moisture_content_of_soil_layer_3")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="moisture_content_of_soil_layer_3", &
          canonicalUnits="kg m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "moisture_content_of_soil_layer_4")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="moisture_content_of_soil_layer_4", &
          canonicalUnits="kg m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "soil_porosity")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="soil_porosity", &
          canonicalUnits="1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "temperature_of_soil_layer")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="temperature_of_soil_layer", &
          canonicalUnits="K", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "temperature_of_soil_layer_1")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="temperature_of_soil_layer_1", &
          canonicalUnits="K", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "temperature_of_soil_layer_2")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="temperature_of_soil_layer_2", &
          canonicalUnits="K", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "temperature_of_soil_layer_3")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="temperature_of_soil_layer_3", &
          canonicalUnits="K", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "temperature_of_soil_layer_4")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="temperature_of_soil_layer_4", &
          canonicalUnits="K", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "soil_temperature_bottom")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="soil_temperature_bottom", &
          canonicalUnits="K", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "soil_type")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="soil_type", &
          canonicalUnits="1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "soil_moisture_content")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="soil_moisture_content", &
          canonicalUnits="kg m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "subsurface_basin_mask")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="subsurface_basin_mask", &
          canonicalUnits="1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "subsurface_runoff_flux")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="subsurface_runoff_flux", &
          canonicalUnits="kg m-2 s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "surface_microwave_emissivity")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="surface_microwave_emissivity", &
          canonicalUnits="1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "surface_runoff_flux")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="surface_runoff_flux", &
          canonicalUnits="kg m-2 s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "vegetation_type")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="vegetation_type", &
          canonicalUnits="1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "volume_fraction_of_frozen_water_in_soil")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="volume_fraction_of_frozen_water_in_soil", &
          canonicalUnits="m3 m-3", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "liquid_water_content_of_soil_layer")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="liquid_water_content_of_soil_layer", &
          canonicalUnits="kg m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "liquid_water_content_of_soil_layer_1")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="liquid_water_content_of_soil_layer_1", &
          canonicalUnits="kg m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "liquid_water_content_of_soil_layer_2")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="liquid_water_content_of_soil_layer_2", &
          canonicalUnits="kg m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "liquid_water_content_of_soil_layer_3")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="liquid_water_content_of_soil_layer_3", &
          canonicalUnits="kg m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "liquid_water_content_of_soil_layer_4")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="liquid_water_content_of_soil_layer_4", &
          canonicalUnits="kg m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "volume_fraction_of_total_water_in_soil")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="volume_fraction_of_total_water_in_soil", &
          canonicalUnits="m3 m-3", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "volume_fraction_of_total_water_in_soil_at_critical_point")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="volume_fraction_of_total_water_in_soil_at_critical_point", &
          canonicalUnits="m3 m-3", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "volume_fraction_of_total_water_in_soil_at_field_capacity")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="volume_fraction_of_total_water_in_soil_at_field_capacity", &
          canonicalUnits="m3 m-3", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "volume_fraction_of_total_water_in_soil_at_wilting_point")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="volume_fraction_of_total_water_in_soil_at_wilting_point", &
          canonicalUnits="m3 m-3", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "water_surface_height_above_reference_datum")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="water_surface_height_above_reference_datum", &
          canonicalUnits="m", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "mean_sensi_heat_flx_atm_into_lnd")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_sensi_heat_flx_atm_into_lnd", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "mean_laten_heat_flx_atm_into_lnd")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="mean_laten_heat_flx_atm_into_lnd", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif

      ! Fields from and to WW3

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "eastward_wind_at_10m_height")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="eastward_wind_at_10m_height", &
          canonicalUnits="m s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      call NUOPC_FieldDictionarySetSyno( &
        standardNames = (/"eastward_wind_at_10m_height",&
                          "inst_zonal_wind_height10m  "/), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "northward_wind_at_10m_height")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="northward_wind_at_10m_height", &
          canonicalUnits="m s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      call NUOPC_FieldDictionarySetSyno( &
        standardNames = (/"northward_wind_at_10m_height",&
                          "inst_merid_wind_height10m   "/), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "sea_ice_concentration")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="sea_ice_concentration", &
          canonicalUnits="1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      call NUOPC_FieldDictionarySetSyno( &
        standardNames = (/"ice_fraction         ",&
                          "sea_ice_concentration"/), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

      !For MOM6 and WW3 variables to match: 
      call NUOPC_FieldDictionarySetSyno( &
        standardNames = (/"surface_eastward_sea_water_velocity",&
                          "ocn_current_zonal                  "/), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      call NUOPC_FieldDictionarySetSyno( &
        standardNames = (/"surface_northward_sea_water_velocity",&
                          "ocn_current_merid                   "/), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "eastward_stokes_drift_current")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="eastward_stokes_drift_current", &
          canonicalUnits="m s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "northward_stokes_drift_current")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="northward_stokes_drift_current", &
          canonicalUnits="m s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "eastward_wave_bottom_current")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="eastward_wave_bottom_current", &
          canonicalUnits="m s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "northward_wave_bottom_current")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="northward_wave_bottom_current", &
          canonicalUnits="m s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "wave_bottom_current_radian_frequency")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="wave_bottom_current_radian_frequency", &
          canonicalUnits="rad s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "eastward_wave_radiation_stress_gradient")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="eastward_wave_radiation_stress_gradient", &
          canonicalUnits="Pa", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "northward_wave_radiation_stress_gradient")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="northward_wave_radiation_stress_gradient", &
          canonicalUnits="Pa", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "eastward_wave_radiation_stress")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="eastward_wave_radiation_stress", &
          canonicalUnits="N m-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "eastward_northward_wave_radiation_stress")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="eastward_northward_wave_radiation_stress", &
          canonicalUnits="N m-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "northward_wave_radiation_stress")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="northward_wave_radiation_stress", &
          canonicalUnits="N m-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "wave_induced_charnock_parameter")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="wave_induced_charnock_parameter", &
          canonicalUnits="1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "wave_z0_roughness_length")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="wave_z0_roughness_length", &
          canonicalUnits="m", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "wave_bottom_current_period")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="wave_bottom_current_period", &
          canonicalUnits="s", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif


      ! Fields from WAM to IPE

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "northward_wind_neutral")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="northward_wind_neutral", &
          canonicalUnits="m s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "eastward_wind_neutral")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="eastward_wind_neutral", &
          canonicalUnits="m s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "upward_wind_neutral")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="upward_wind_neutral", &
          canonicalUnits="m s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "temp_neutral")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="temp_neutral", &
          canonicalUnits="K", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "O_Density")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="O_Density", &
          canonicalUnits="m-3", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "O2_Density")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="O2_Density", &
          canonicalUnits="m-3", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "N2_Density")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="N2_Density", &
          canonicalUnits="m-3", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "height")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="height", &
          canonicalUnits="km", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif

      ! Chemistry fields

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "inst_pres_interface")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_pres_interface", &
          canonicalUnits="Pa", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "inst_pres_levels")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_pres_levels", &
          canonicalUnits="Pa", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "inst_geop_interface")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_geop_interface", &
          canonicalUnits="m2 s-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "inst_geop_levels")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_geop_levels", &
          canonicalUnits="m2 s-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "inst_temp_levels")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_temp_levels", &
          canonicalUnits="K", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "inst_zonal_wind_levels")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_zonal_wind_levels", &
          canonicalUnits="m s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "inst_merid_wind_levels")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_merid_wind_levels", &
          canonicalUnits="m s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "inst_omega_levels")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_omega_levels", &
          canonicalUnits="Pa s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "inst_tracer_mass_frac")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_tracer_mass_frac", &
          canonicalUnits="kg kg-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "inst_tracer_up_surface_flx")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_tracer_up_surface_flx", &
          canonicalUnits="kg m-2 s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "inst_tracer_down_surface_flx")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_tracer_down_surface_flx", &
          canonicalUnits="kg m-2 s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "inst_tracer_clmn_mass_dens")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_tracer_clmn_mass_dens", &
          canonicalUnits="g m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "inst_tracer_anth_biom_flx")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_tracer_anth_biom_flx", &
          canonicalUnits="ug m-2 s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "inst_pbl_height")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_pbl_height", &
          canonicalUnits="m", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "surface_cell_area")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="surface_cell_area", &
          canonicalUnits="m2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "inst_convective_rainfall_amount")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_convective_rainfall_amount", &
          canonicalUnits="kg m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "inst_exchange_coefficient_heat_levels")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_exchange_coefficient_heat_levels", &
          canonicalUnits="W m-2 K-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "inst_spec_humid_conv_tendency_levels")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_spec_humid_conv_tendency_levels", &
          canonicalUnits="kg kg-1 s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "inst_friction_velocity")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_friction_velocity", &
          canonicalUnits="m s-1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "inst_rainfall_amount")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_rainfall_amount", &
          canonicalUnits="kg m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "inst_soil_moisture_content")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_soil_moisture_content", &
          canonicalUnits="kg m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "inst_up_sensi_heat_flx")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_up_sensi_heat_flx", &
          canonicalUnits="W m-2", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "inst_lwe_snow_thickness")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_lwe_snow_thickness", &
          canonicalUnits="m", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "inst_vegetation_area_frac")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_vegetation_area_frac", &
          canonicalUnits="1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "inst_surface_roughness")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="inst_surface_roughness", &
          canonicalUnits="1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif

!JQI20211004
      ! Fields from WW3 to FVCOM

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "significant_wave_height")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="significant_wave_height", &
          canonicalUnits="m", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "average_wave_length")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="average_wave_length", &
          canonicalUnits="m", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif

      if (.not.NUOPC_FieldDictionaryHasEntry( &
        "average_wave_direction")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="average_wave_direction", &
          canonicalUnits="degree", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif

!JQI20211004

      ! Dummy fields

      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "dummyfield")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="dummyfield", &
          canonicalUnits="1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "dummyfield1")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="dummyfield1", &
          canonicalUnits="1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      if (.not. NUOPC_FieldDictionaryHasEntry( &
        "dummyfield2")) then
        call NUOPC_FieldDictionaryAddEntry( &
          standardName="dummyfield2", &
          canonicalUnits="1", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
#endif

!-----------------------------------------------------------------------
!
      END SUBROUTINE EARTH_REGISTER
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!

      subroutine SetModelServices(driver, rc)
#ifdef CMEPS
        use med_internalstate_mod , only : med_id
#endif
        type(ESMF_GridComp)  :: driver
        integer, intent(out) :: rc

        ! local variables
        integer                         :: localrc, stat, i, j, petCount
        character(ESMF_MAXSTR)          :: name
        type(WRAP_EARTH_INTERNAL_STATE) :: is
        type(ESMF_GridComp)             :: comp
        type(ESMF_Config)               :: config
        character(len=32), allocatable  :: compLabels(:)
        integer, allocatable            :: petList(:)
        character(len=10)               :: value
        character(len=20)               :: model, prefix
        character(len=160)              :: msg
        integer                         :: petListBounds(2)
        integer                         :: componentCount
        type(NUOPC_FreeFormat)          :: attrFF, fdFF
#ifdef CMEPS
        logical                         :: read_restart
        character(ESMF_MAXSTR)          :: cvalue
        character(len=5)                :: inst_suffix
        logical                         :: isPresent
#endif
        rc = ESMF_SUCCESS

        ! query the Component for info
        call ESMF_GridCompGet(driver, name=name, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out

        ! allocate memory for the internal state and store in Component
        allocate(is%EARTH_INT_STATE, stat=stat)
        if (ESMF_LogFoundAllocError(statusToCheck=stat, &
          msg="Allocation of internal state memory failed.", &
          line=__LINE__, file=trim(name)//":"//__FILE__, rcToReturn=rc)) &
          return  ! bail out
        call ESMF_GridCompSetInternalState(driver, is, rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
        
        ! get petCount and config
        call ESMF_GridCompGet(driver, petCount=petCount, config=config, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
        
        ! read and ingest free format driver attributes
        attrFF = NUOPC_FreeFormatCreate(config, label="EARTH_attributes::", &
          relaxedflag=.true., rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
        call NUOPC_CompAttributeIngest(driver, attrFF, addFlag=.true., rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
        call NUOPC_FreeFormatDestroy(attrFF, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
        
        ! dump the current field dictionary into the Log file
        call ESMF_AttributeGet(driver, name="DumpFieldDictionary", &
          value=value, defaultValue="false", &
          convention="NUOPC", purpose="Instance", rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
        if (trim(value)=="true") then
          call ESMF_LogWrite( &
            "===>===>===>===> Begin Dumping Field Dictionary <===<===<===<===",&
            ESMF_LOGMSG_INFO, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
          call NUOPC_FieldDictionaryEgest(fdFF, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
          call NUOPC_FreeFormatLog(fdFF, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
          call ESMF_LogWrite( &
            "===>===>===>===> Done Dumping Field Dictionary <===<===<===<===", &
            ESMF_LOGMSG_INFO, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
        endif
        
        ! determine the generic component labels
        componentCount = ESMF_ConfigGetLen(config, &
          label="EARTH_component_list:", rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
        allocate(compLabels(componentCount), stat=stat)
        if (ESMF_LogFoundAllocError(statusToCheck=stat, &
          msg="Allocation of compLabels failed.", &
          line=__LINE__, file=trim(name)//":"//__FILE__, rcToReturn=rc)) &
          return  ! bail out
        call ESMF_ConfigGetAttribute(config, valueList=compLabels, &
          label="EARTH_component_list:", count=componentCount, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out

#ifdef CMEPS
        inst_suffix = ""

        ! obtain driver attributes (for CMEPS)
        call ReadAttributes(driver, config, "DRIVER_attributes::", formatprint=.true., rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out

        call ReadAttributes(driver, config, "ALLCOMP_attributes::", formatprint=.true., rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
#endif

        ! determine information for each component and add to the driver
        do i=1, componentCount
          ! construct component prefix
          prefix=trim(compLabels(i))
          ! read in petList bounds
          call ESMF_ConfigGetAttribute(config, petListBounds, &
            label=trim(prefix)//"_petlist_bounds:", default=-1, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
          ! handle the default situation
          if (petListBounds(1)==-1 .or. petListBounds(2)==-1) then
            petListBounds(1) = 0
            petListBounds(2) = petCount - 1
          endif
          ! read in model instance name
          call ESMF_ConfigGetAttribute(config, model, &
            label=trim(prefix)//"_model:", default="none", rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
          ! check that there was a model instance specified
          if (trim(model) == "none") then
            ! Error condition: no model was specified
            write (msg, *) "No model was specified for component: ",trim(prefix)
            call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msg, line=__LINE__, &
              file=__FILE__, rcToReturn=rc)
            return  ! bail out
          endif
          ! set petList for this component
          allocate(petList(petListBounds(2)-petListBounds(1)+1))
          do j=petListBounds(1), petListBounds(2)
            petList(j-petListBounds(1)+1) = j ! PETs are 0 based
          enddo

          if (trim(model) == "satm") then
#ifdef FRONT_SATM
            call NUOPC_DriverAddComp(driver, trim(prefix), SATM_SS, &
              petList=petList, comp=comp, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
#else
            write (msg, *) "Model '", trim(model), "' was requested, "// &
              "but is not available in the executable!"
            call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msg, line=__LINE__, &
              file=__FILE__, rcToReturn=rc)
            return  ! bail out
#endif
          elseif (trim(model) == "xatm") then
#ifdef FRONT_XATM
            call NUOPC_DriverAddComp(driver, trim(prefix), XATM_SS, &
              petList=petList, comp=comp, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
#else
            write (msg, *) "Model '", trim(model), "' was requested, "// &
              "but is not available in the executable!"
            call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msg, line=__LINE__, &
              file=__FILE__, rcToReturn=rc)
            return  ! bail out
#endif
          elseif (trim(model) == "datawam") then
#ifdef FRONT_DATAWAM
            call NUOPC_DriverAddComp(driver, trim(prefix), DATAWAM_SS, &
              petList=petList, comp=comp, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
#else
            write (msg, *) "Model '", trim(model), "' was requested, "// &
              "but is not available in the executable!"
            call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msg, line=__LINE__, &
              file=__FILE__, rcToReturn=rc)
            return  ! bail out
#endif
          elseif (trim(model) == "gsm") then
#ifdef FRONT_GSM
            call NUOPC_DriverAddComp(driver, trim(prefix), GSM_SS, &
              petList=petList, comp=comp, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
#else
            write (msg, *) "Model '", trim(model), "' was requested, "// &
              "but is not available in the executable!"
            call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msg, line=__LINE__, &
              file=__FILE__, rcToReturn=rc)
            return  ! bail out
#endif
          elseif (trim(model) == "nmmb") then
#ifdef FRONT_NMMB
            call NUOPC_DriverAddComp(driver, trim(prefix), NMMB_SS, &
              petList=petList, comp=comp, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
#else
            write (msg, *) "Model '", trim(model), "' was requested, "// &
              "but is not available in the executable!"
            call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msg, line=__LINE__, &
              file=__FILE__, rcToReturn=rc)
            return  ! bail out
#endif
          elseif (trim(model) == "fv3") then
#ifdef FRONT_FV3
            call NUOPC_DriverAddComp(driver, trim(prefix), FV3_SS, &
              petList=petList, comp=comp, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
#else
            write (msg, *) "Model '", trim(model), "' was requested, "// &
              "but is not available in the executable!"
            call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msg, line=__LINE__, &
              file=__FILE__, rcToReturn=rc)
            return  ! bail out
#endif
          elseif (trim(model) == "datm") then
#ifdef FRONT_DATM
            call NUOPC_DriverAddComp(driver, trim(prefix), DATM_SS, &
              petList=petList, comp=comp, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
#else
            write (msg, *) "Model '", trim(model), "' was requested, "// &
              "but is not available in the executable!"
            call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msg, line=__LINE__, &
              file=__FILE__, rcToReturn=rc)
            return  ! bail out
#endif
          elseif (trim(model) == "socn") then
#ifdef FRONT_SOCN
            call NUOPC_DriverAddComp(driver, trim(prefix), SOCN_SS, &
              petList=petList, comp=comp, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
#else
            write (msg, *) "Model '", trim(model), "' was requested, "// &
              "but is not available in the executable!"
            call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msg, line=__LINE__, &
              file=__FILE__, rcToReturn=rc)
            return  ! bail out
#endif
          elseif (trim(model) == "xocn") then
#ifdef FRONT_XOCN
            call NUOPC_DriverAddComp(driver, trim(prefix), XOCN_SS, &
              petList=petList, comp=comp, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
#else
            write (msg, *) "Model '", trim(model), "' was requested, "// &
              "but is not available in the executable!"
            call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msg, line=__LINE__, &
              file=__FILE__, rcToReturn=rc)
            return  ! bail out
#endif
          elseif (trim(model) == "hycom") then
#ifdef FRONT_HYCOM
            call NUOPC_DriverAddComp(driver, trim(prefix), HYCOM_SS, &
              petList=petList, comp=comp, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
#else
            write (msg, *) "Model '", trim(model), "' was requested, "// &
              "but is not available in the executable!"
            call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msg, line=__LINE__, &
              file=__FILE__, rcToReturn=rc)
            return  ! bail out
#endif
          elseif (trim(model) == "mom5") then
#ifdef FRONT_MOM5
            call NUOPC_DriverAddComp(driver, trim(prefix), MOM5_SS, &
              petList=petList, comp=comp, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
#else
            write (msg, *) "Model '", trim(model), "' was requested, "// &
              "but is not available in the executable!"
            call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msg, line=__LINE__, &
              file=__FILE__, rcToReturn=rc)
            return  ! bail out
#endif
          elseif (trim(model) == "mom6") then
#ifdef FRONT_MOM6
            call NUOPC_DriverAddComp(driver, trim(prefix), MOM6_SS, &
              petList=petList, comp=comp, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
#else
            write (msg, *) "Model '", trim(model), "' was requested, "// &
              "but is not available in the executable!"
            call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msg, line=__LINE__, &
              file=__FILE__, rcToReturn=rc)
            return  ! bail out
#endif
          elseif (trim(model) == "pom") then
#ifdef FRONT_POM
            call NUOPC_DriverAddComp(driver, trim(prefix), POM_SS, &
              petList=petList, comp=comp, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
#else
            write (msg, *) "Model '", trim(model), "' was requested, "// &
              "but is not available in the executable!"
            call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msg, line=__LINE__, &
              file=__FILE__, rcToReturn=rc)
            return  ! bail out
#endif
          elseif (trim(model) == "sice") then
#ifdef FRONT_SICE
            call NUOPC_DriverAddComp(driver, trim(prefix), SICE_SS, &
              petList=petList, comp=comp, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
#else
            write (msg, *) "Model '", trim(model), "' was requested, "// &
              "but is not available in the executable!"
            call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msg, line=__LINE__, &
              file=__FILE__, rcToReturn=rc)
            return  ! bail out
#endif
          elseif (trim(model) == "xice") then
#ifdef FRONT_XICE
            call NUOPC_DriverAddComp(driver, trim(prefix), XICE_SS, &
              petList=petList, comp=comp, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
#else
            write (msg, *) "Model '", trim(model), "' was requested, "// &
              "but is not available in the executable!"
            call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msg, line=__LINE__, &
              file=__FILE__, rcToReturn=rc)
            return  ! bail out
#endif
          elseif (trim(model) == "cice") then
#ifdef FRONT_CICE
            call NUOPC_DriverAddComp(driver, trim(prefix), CICE_SS, &
              petList=petList, comp=comp, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
#else
            write (msg, *) "Model '", trim(model), "' was requested, "// &
              "but is not available in the executable!"
            call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msg, line=__LINE__, &
              file=__FILE__, rcToReturn=rc)
            return  ! bail out
#endif
          elseif (trim(model) == "cice6") then
#ifdef FRONT_CICE6
            call NUOPC_DriverAddComp(driver, trim(prefix), CICE6_SS, &
              petList=petList, comp=comp, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
#else
            write (msg, *) "Model '", trim(model), "' was requested, "// &
              "but is not available in the executable!"
            call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msg, line=__LINE__, &
              file=__FILE__, rcToReturn=rc)
            return  ! bail out
#endif
          elseif (trim(model) == "swav") then
#ifdef FRONT_SWAV
            call NUOPC_DriverAddComp(driver, trim(prefix), SWAV_SS, &
              petList=petList, comp=comp, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
#else
            write (msg, *) "Model '", trim(model), "' was requested, "// &
              "but is not available in the executable!"
            call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msg, line=__LINE__, &
              file=__FILE__, rcToReturn=rc)
            return  ! bail out
#endif
          elseif (trim(model) == "xwav") then
#ifdef FRONT_XWAV
            call NUOPC_DriverAddComp(driver, trim(prefix), XWAV_SS, &
              petList=petList, comp=comp, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
#else
            write (msg, *) "Model '", trim(model), "' was requested, "// &
              "but is not available in the executable!"
            call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msg, line=__LINE__, &
              file=__FILE__, rcToReturn=rc)
            return  ! bail out
#endif
          elseif (trim(model) == "ww3") then
#ifdef FRONT_WW3
            call NUOPC_DriverAddComp(driver, trim(prefix), WW3_SS, &
              petList=petList, comp=comp, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
#else
            write (msg, *) "Model '", trim(model), "' was requested, "// &
              "but is not available in the executable!"
            call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msg, line=__LINE__, &
              file=__FILE__, rcToReturn=rc)
            return  ! bail out
#endif
          elseif (trim(model) == "slnd") then
#ifdef FRONT_SLND
            call NUOPC_DriverAddComp(driver, trim(prefix), SLND_SS, &
              petList=petList, comp=comp, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
#else
            write (msg, *) "Model '", trim(model), "' was requested, "// &
              "but is not available in the executable!"
            call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msg, line=__LINE__, &
              file=__FILE__, rcToReturn=rc)
            return  ! bail out
#endif
          elseif (trim(model) == "xlnd") then
#ifdef FRONT_XLND
            call NUOPC_DriverAddComp(driver, trim(prefix), XLND_SS, &
              petList=petList, comp=comp, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
#else
            write (msg, *) "Model '", trim(model), "' was requested, "// &
              "but is not available in the executable!"
            call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msg, line=__LINE__, &
              file=__FILE__, rcToReturn=rc)
            return  ! bail out
#endif
          elseif (trim(model) == "noah") then
#ifdef FRONT_NOAH
            call NUOPC_DriverAddComp(driver, trim(prefix), NOAH_SS, &
              petList=petList, comp=comp, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
#else
            write (msg, *) "Model '", trim(model), "' was requested, "// &
              "but is not available in the executable!"
            call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msg, line=__LINE__, &
              file=__FILE__, rcToReturn=rc)
            return  ! bail out
#endif
          elseif (trim(model) == "lis") then
#ifdef FRONT_LIS
            call NUOPC_DriverAddComp(driver, trim(prefix), LIS_SS, &
              petList=petList, comp=comp, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
#else
            write (msg, *) "Model '", trim(model), "' was requested, "// &
              "but is not available in the executable!"
            call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msg, line=__LINE__, &
              file=__FILE__, rcToReturn=rc)
            return  ! bail out
#endif
          elseif (trim(model) == "sipm") then
#ifdef FRONT_SIPM
            call NUOPC_DriverAddComp(driver, trim(prefix), SIPM_SS, &
              petList=petList, comp=comp, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
#else
            write (msg, *) "Model '", trim(model), "' was requested, "// &
              "but is not available in the executable!"
            call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msg, line=__LINE__, &
              file=__FILE__, rcToReturn=rc)
            return  ! bail out
#endif
          elseif (trim(model) == "xipm") then
#ifdef FRONT_XIPM
            call NUOPC_DriverAddComp(driver, trim(prefix), XIPM_SS, &
              petList=petList, comp=comp, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
#else
            write (msg, *) "Model '", trim(model), "' was requested, "// &
              "but is not available in the executable!"
            call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msg, line=__LINE__, &
              file=__FILE__, rcToReturn=rc)
            return  ! bail out
#endif
          elseif (trim(model) == "ipe") then
#ifdef FRONT_IPE
            call NUOPC_DriverAddComp(driver, trim(prefix), IPE_SS, &
              petList=petList, comp=comp, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
#else
            write (msg, *) "Model '", trim(model), "' was requested, "// &
              "but is not available in the executable!"
            call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msg, line=__LINE__, &
              file=__FILE__, rcToReturn=rc)
            return  ! bail out
#endif
          elseif (trim(model) == "dataipe") then
#ifdef FRONT_DATAIPE
            call NUOPC_DriverAddComp(driver, trim(prefix), DATAIPE_SS, &
              petList=petList, comp=comp, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
#else
            write (msg, *) "Model '", trim(model), "' was requested, "// &
              "but is not available in the executable!"
            call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msg, line=__LINE__, &
              file=__FILE__, rcToReturn=rc)
            return  ! bail out
#endif
          elseif (trim(model) == "shyd") then
#ifdef FRONT_SHYD
            call NUOPC_DriverAddComp(driver, trim(prefix), SHYD_SS, &
              petList=petList, comp=comp, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
#else
            write (msg, *) "Model '", trim(model), "' was requested, "// &
              "but is not available in the executable!"
            call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msg, line=__LINE__, &
              file=__FILE__, rcToReturn=rc)
            return  ! bail out
#endif
          elseif (trim(model) == "xhyd") then
#ifdef FRONT_XHYD
            call NUOPC_DriverAddComp(driver, trim(prefix), XHYD_SS, &
              petList=petList, comp=comp, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
#else
            write (msg, *) "Model '", trim(model), "' was requested, "// &
              "but is not available in the executable!"
            call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msg, line=__LINE__, &
              file=__FILE__, rcToReturn=rc) 
            return  ! bail out
#endif
          elseif (trim(model) == "wrfhydro") then
#ifdef FRONT_WRFHYDRO
            call NUOPC_DriverAddComp(driver, trim(prefix), WRFHYDRO_SS, &
              petList=petList, comp=comp, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
#else
            write (msg, *) "Model '", trim(model), "' was requested, "// &
              "but is not available in the executable!"
            call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msg, line=__LINE__, &
              file=__FILE__, rcToReturn=rc)
            return  ! bail out
#endif
          elseif (trim(model) == "gsdchem") then
#ifdef FRONT_GSDCHEM
            call NUOPC_DriverAddComp(driver, trim(prefix), GSDCHEM_SS, &
              petList=petList, comp=comp, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, file=trim(name)//":"//__FILE__)) return  !  bail out
#else
            write (msg, *) "Model '", trim(model), "' was requested, "// &
              "but is not available in the executable!"
            call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msg, line=__LINE__, &
              file=__FILE__, rcToReturn=rc)
            return  ! bail out
#endif
          elseif (trim(model) == "schism") then
#ifdef FRONT_SCHISM
            call NUOPC_DriverAddComp(driver, trim(prefix), SCHISM_SS, &
              petList=petList, comp=comp, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, file=trim(name)//":"//__FILE__)) return  !  bail out
#else
            write (msg, *) "Model '", trim(model), "' was requested, "// &
              "but is not available in the executable!"
            call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msg, line=__LINE__, &
              file=__FILE__, rcToReturn=rc)
            return  ! bail out
#endif
          elseif (trim(model) == "adcirc") then
#ifdef FRONT_ADCIRC
            call NUOPC_DriverAddComp(driver, trim(prefix), ADCIRC_SS, &
              petList=petList, comp=comp, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, file=trim(name)//":"//__FILE__)) return  !  bail out
#else
            write (msg, *) "Model '", trim(model), "' was requested, "// &
              "but is not available in the executable!"
            call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msg, line=__LINE__, &
              file=__FILE__, rcToReturn=rc)
            return  ! bail out
#endif
          elseif (trim(model) == "hwrfdata") then
#ifdef FRONT_HWRFDATA
            call NUOPC_DriverAddComp(driver, trim(prefix), HWRFDATA_SS, &
              petList=petList, comp=comp, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, file=trim(name)//":"//__FILE__)) return  !  bail out
#else
            write (msg, *) "Model '", trim(model), "' was requested, "// &
              "but is not available in the executable!"
            call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msg, line=__LINE__, &
              file=__FILE__, rcToReturn=rc)
            return  ! bail out
#endif
          elseif (trim(model) == "ww3data") then
#ifdef FRONT_WW3DATA
            call NUOPC_DriverAddComp(driver, trim(prefix), WW3DATA_SS, &
              petList=petList, comp=comp, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, file=trim(name)//":"//__FILE__)) return  !  bail out
#else
            write (msg, *) "Model '", trim(model), "' was requested, "// &
              "but is not available in the executable!"
            call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msg, line=__LINE__, &
              file=__FILE__, rcToReturn=rc)
            return  ! bail out
#endif
          elseif (trim(model) == "bardata") then
#ifdef FRONT_BARDATA
            call NUOPC_DriverAddComp(driver, trim(prefix), BARDATA_SS, &
              petList=petList, comp=comp, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, file=trim(name)//":"//__FILE__)) return  !  bail out
#else
            write (msg, *) "Model '", trim(model), "' was requested, "// &
              "but is not available in the executable!"
            call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msg, line=__LINE__, &
              file=__FILE__, rcToReturn=rc)
            return  ! bail out
#endif
          elseif (trim(model) == "atmesh") then
#ifdef FRONT_ATMESH
            call NUOPC_DriverAddComp(driver, trim(prefix), ATMESH_SS, &
              petList=petList, comp=comp, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, file=trim(name)//":"//__FILE__)) return  !  bail out
#else
            write (msg, *) "Model '", trim(model), "' was requested, "// &
              "but is not available in the executable!"
            call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msg, line=__LINE__, &
              file=__FILE__, rcToReturn=rc)
            return  ! bail out
#endif
          elseif (trim(model) == "pahm") then
#ifdef FRONT_PAHM
            call NUOPC_DriverAddComp(driver, trim(prefix), PAHM_SS, &
              petList=petList, comp=comp, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, file=trim(name)//":"//__FILE__)) return  !  bail out
#else
            write (msg, *) "Model '", trim(model), "' was requested, "// &
              "but is not available in the executable!"
            call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msg, line=__LINE__, &
              file=__FILE__, rcToReturn=rc)
            return  ! bail out
#endif
          elseif (trim(model) == "nwm") then
#ifdef FRONT_NWM
            call NUOPC_DriverAddComp(driver, trim(prefix), NWM_SS, &
              petList=petList, comp=comp, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, file=trim(name)//":"//__FILE__)) return  !  bail out
#else
            write (msg, *) "Model '", trim(model), "' was requested, "// &
              "but is not available in the executable!"
            call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msg, line=__LINE__, &
              file=__FILE__, rcToReturn=rc)
            return  ! bail out
#endif
          elseif (trim(model) == "fvcom_cap") then
#ifdef FRONT_FVCOM
            call NUOPC_DriverAddComp(driver, trim(prefix), FVCOMSS, &
              petList=petList, comp=comp, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, file=trim(name)//":"//__FILE__)) return  !bail out

            call ESMF_ConfigGetAttribute(config, fvcom_name, &
                label="fvcom_name:",default=' ', rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, file=trim(name)//":"//__FILE__)) return !bail out 

            fvcom_pet_num = size(petList)
            print*,'FVCOM pet = ', fvcom_pet_num
#else
            write (msg, *) "Model '", trim(model), "' was requested, "// &
              "but is not available in the executable!"
            call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msg,line=__LINE__, &
              file=__FILE__, rcToReturn=rc)
            return  ! bail out
#endif

! - Two mediator choices currently built into NEMS from internal
          elseif (trim(model) == "nems") then
#ifdef CMEPS
            med_id = i+1
#endif
            call NUOPC_DriverAddComp(driver, trim(prefix), MED_SS, &
              petList=petList, comp=comp, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
          elseif (trim(model) == "spaceweather") then
            call NUOPC_DriverAddComp(driver, trim(prefix), MEDSW_SS, &
              petList=petList, comp=comp, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
             line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
          else
            ! Error condition: unknown model requested
            write (msg, *) "The requested model '", trim(model), &
              "' is an invalid choice!"
            call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msg, line=__LINE__, &
              file=__FILE__, rcToReturn=rc)
            return  ! bail out
          endif
          
          ! read and ingest free format component attributes
          attrFF = NUOPC_FreeFormatCreate(config, &
            label=trim(prefix)//"_attributes::", relaxedflag=.true., rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
          call NUOPC_CompAttributeIngest(comp, attrFF, addFlag=.true., rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
          call NUOPC_FreeFormatDestroy(attrFF, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
          
          ! clean-up
          deallocate(petList)

#ifdef CMEPS
        ! Perform restarts if appropriate
        call InitRestart(driver, rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out

        call AddAttributes(comp, driver, config, i+1, trim(prefix), inst_suffix, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
#endif          
        enddo

#if ESMF_VERSION_MAJOR < 8
!TODOgjt: REMOVE THIS BLOCK ONCE SHOWN TO WORK WITHOUT
        ! SetServices for Connectors
        call SetFromConfig(driver, mode="setServicesConnectors", rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
#endif

        ! clean-up
        deallocate(compLabels)
        
      end subroutine

  !-----------------------------------------------------------------------------
  
  subroutine SetRunSequence(driver, rc)
    type(ESMF_GridComp)  :: driver
    integer, intent(out) :: rc
    
    ! local variables
    character(ESMF_MAXSTR)          :: name
#if ESMF_VERSION_MAJOR >= 8
    type(ESMF_Config)               :: config
    type(NUOPC_FreeFormat)          :: runSeqFF
#endif

    rc = ESMF_SUCCESS

    ! query the Component for info
    call ESMF_GridCompGet(driver, name=name, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out

#if ESMF_VERSION_MAJOR >= 8
    ! read free format run sequence from config
    call ESMF_GridCompGet(driver, config=config, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
    runSeqFF = NUOPC_FreeFormatCreate(config, label="runSeq::", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out

    ! ingest FreeFormat run sequence
    call NUOPC_DriverIngestRunSequence(driver, runSeqFF, &
      autoAddConnectors=.true., rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
#else
    ! access runSeq in the config
    call SetFromConfig(driver, mode="setRunSequence", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
#endif
    
    ! Diagnostic output
    if(verbose_diagnostics()) then
       call NUOPC_DriverPrint(driver, orderflag=.true., rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
    endif
    
  end subroutine
    
  !-----------------------------------------------------------------------------

#if ESMF_VERSION_MAJOR < 8
!TODOgjt: REMOVE THIS BLOCK ONCE SHOWN TO WORK WITHOUT
  subroutine SetFromConfig(driver, mode, rc)
    type(ESMF_GridComp)   :: driver
    character(len=*)      :: mode
    integer, intent(out)  :: rc
    
    ! local variables
    character(ESMF_MAXSTR)          :: name
    type(ESMF_Config)               :: config
    integer                         :: lineCount, columnCount, i, slotCount
    integer, allocatable            :: count(:)
    character(len=32), allocatable  :: line(:)
    character(len=32)               :: tempString
    logical                         :: phaseFlag
    integer                         :: level, slot, slotHWM
    real(ESMF_KIND_R8)              :: seconds
    integer, allocatable            :: slotStack(:)
    type(ESMF_TimeInterval)         :: timeStep
    type(ESMF_Clock)                :: internalClock, subClock
    character(len=60), allocatable  :: connectorInstance(:)
    integer                         :: connectorCount, j
    type(ESMF_CplComp)              :: conn

    character(len=ESMF_MAXSTR) :: msgString
    character(len=10)          :: value

    !can set to 'max' to recover intro/extro CurrGarbInfo for 
    !all connectors 
    character(len=10)          :: defaultVerbosity = "0"
    !character(len=10)          :: defaultVerbosity = "max"

    rc = ESMF_SUCCESS
    
    ! query the Component for info
    call ESMF_GridCompGet(driver, name=name, config=config, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out

    ! reset config to beginning of runSeq:: block
    call ESMF_ConfigFindLabel(config, label="runSeq::", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
    call ESMF_ConfigGetDim(config, lineCount, columnCount, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
    
    allocate(count(lineCount))
    
    if (trim(mode)=="setServicesConnectors") then
      allocate(connectorInstance(lineCount))  ! max number of connectors
      connectorCount = 0 ! reset
    write(msgString,'(a,i6)')'max number of connectors ',lineCount
    call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
    endif
    
    ! reset config to beginning of runSeq:: block
    call ESMF_ConfigFindLabel(config, label="runSeq::", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out

    ! determine number of entries on each line
    do i=1, lineCount
      call ESMF_ConfigNextLine(config)
      count(i) = ESMF_ConfigGetLen(config) ! entries on line i
    enddo
    
    ! reset config to beginning of runSeq:: block
    call ESMF_ConfigFindLabel(config, label="runSeq::", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out

    ! read each line and determine slotCount
    slotCount = 0
    do i=1, lineCount
      call ESMF_ConfigNextLine(config)
      allocate(line(count(i)))
      call ESMF_ConfigGetAttribute(config, line, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
      
      ! process the configuration line
      if (size(line) == 1) then
        if (index(trim(line(1)),"@") == 1) then
          slotCount = slotCount + 1
        endif
      elseif ((size(line) == 3) .or. (size(line) == 4)) then
        if (trim(mode)=="setServicesConnectors") then
          ! a connector if the second element is "->"
          if (trim(line(2)) /= "->") then
            call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_BAD, &
              msg="Configuration line incorrectly formatted.", &
              line=__LINE__, &
              file=__FILE__)
            return  ! bail out
          else
            ! found a connector entry, see if it is the first instance
            do j=1, connectorCount
              if (trim(connectorInstance(j)) == &
                trim(line(1))//trim(line(2))//trim(line(3))) exit
            enddo
            if (j>connectorCount) then
              ! this is a new Connector instance
              connectorCount = j
              connectorInstance(j) = trim(line(1))//trim(line(2))//trim(line(3))
              write(msgString,'(a,i4,a,i4,4a)')'Connector j = ',j,&
                                                ' line number ', i,&
                                                '  ',trim(connectorInstance(j)),&
                                                ' Verbosity = ',trim(defaultVerbosity)
              call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
              ! SetServices for new Connector instance
              call NUOPC_DriverAddComp(driver, &
                srcCompLabel=trim(line(1)), dstCompLabel=trim(line(3)), &
                compSetServicesRoutine=conSS, comp=conn, rc=rc)
              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
                line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail
              call NUOPC_CompAttributeSet(conn, name="Verbosity", value=defaultVerbosity, &
                rc=rc)
              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
                line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail
              if (size(line) == 4) then
                ! there are additional connection options specified
                ! -> set as Attribute for now on the connector object
                call ESMF_AttributeSet(conn, name="ConnectionOptions", &
                  value=trim(line(4)), rc=rc)
                if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
                  line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail
              endif
            endif
          endif
        endif
      endif
      ! clean-up
      deallocate(line)
    enddo
    slotCount = (slotCount+1) / 2
    slotCount = max(slotCount, 1) ! at least one slot
    
    if (trim(mode)=="setRunSequence") then
    
      allocate(slotStack(slotCount))

      ! Replace the default RunSequence with a customized one
      call NUOPC_DriverNewRunSequence(driver, slotCount=slotCount, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

      ! Get driver intenalClock
      call ESMF_GridCompGet(driver, clock=internalClock, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

      ! reset config to beginning of runSeq:: block
      call ESMF_ConfigFindLabel(config, label="runSeq::", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out

      level = 0
      slot = 0
      slotHWM = 0
      do i=1, lineCount
        call ESMF_ConfigNextLine(config)
        allocate(line(count(i)))
        call ESMF_ConfigGetAttribute(config, line, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
        
        ! process the configuration line
        if ((size(line) < 1) .or. (size(line) > 4)) then
          call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_BAD, &
            msg="Configuration line incorrectly formatted.", &
            line=__LINE__, &
            file=__FILE__)
          return  ! bail out
        elseif (size(line) == 1) then
          ! either a model or a time step indicator
          if (index(trim(line(1)),"@") == 1) then
            ! time step indicator
            tempString=trim(line(1))
            if (len(trim(tempString)) > 1) then
              ! entering new time loop level
              level = level + 1
              slotStack(level)=slot
              slot = slotHWM + 1
              slotHWM = slotHWM + 1
              read(tempString(2:len(tempString)), *) seconds
              !print *, "found time step indicator: ", seconds
              call ESMF_TimeIntervalSet(timeStep, s_r8=seconds, rc=rc)
              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
                line=__LINE__, &
                file=__FILE__)) &
                return  ! bail out
              if (slot==1) then
                ! Set the timeStep of the internalClock
                call ESMF_ClockSet(internalClock, timeStep=timeStep, rc=rc)
                if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
                  line=__LINE__, &
                  file=__FILE__)) &
                  return  ! bail out
              else
                ! Insert the link to a new slot, and set the timeStep
                call NUOPC_DriverAddRunElement(driver, slot=slotStack(level), &
                  linkSlot=slot, rc=rc)
                if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
                  line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
                subClock = ESMF_ClockCreate(internalClock, rc=rc)  ! make a copy first
                if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
                  line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
                call ESMF_ClockSet(subClock, timeStep=timeStep, rc=rc)
                if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
                  line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
                call NUOPC_DriverSetRunSequence(driver, slot=slot, &
                  clock=subClock, rc=rc)
                if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
                  line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
              endif
            else
              ! exiting time loop level
              slot = slotStack(level)
              level = level - 1
            endif
          else
            ! model
            slot = max(slot, 1) ! model outside of a time loop
            call NUOPC_DriverAddRunElement(driver, slot=slot, &
              compLabel=trim(line(1)), rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, &
              file=__FILE__)) &
              return  ! bail out
          endif
        elseif (size(line) == 2) then
          ! a model with a specific phase label
          call NUOPC_DriverAddRunElement(driver, slot=slot, &
            compLabel=trim(line(1)), phaseLabel=trim(line(2)), rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
        elseif ((size(line) == 3) .or. (size(line) == 4)) then
          ! a connector if the second element is "->", with options if 4th part
          if (trim(line(2)) /= "->") then
            call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_BAD, &
              msg="Configuration line incorrectly formatted.", &
              line=__LINE__, &
              file=__FILE__)
            return  ! bail out
          endif
          call NUOPC_DriverAddRunElement(driver, slot=slot, &
            srcCompLabel=trim(line(1)), dstCompLabel=trim(line(3)), rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
        endif    
        
        ! clean-up
        deallocate(line)
      enddo
      ! clean-up
      deallocate(slotStack)
    endif

    ! clean-up
    deallocate(count)
    if (trim(mode)=="setServicesConnectors") then
      deallocate(connectorInstance)
    endif

  end subroutine
#endif
  !-----------------------------------------------------------------------------

  subroutine Finalize(driver, rc)
    type(ESMF_GridComp)  :: driver
    integer, intent(out) :: rc
    
    ! local variables
    integer                         :: localrc, stat
    type(WRAP_EARTH_INTERNAL_STATE) :: is
    logical                         :: existflag
    character(ESMF_MAXSTR)          :: name

    rc = ESMF_SUCCESS

    ! query the Component for info
    call ESMF_GridCompGet(driver, name=name, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
    
    ! query Component for this internal State
    nullify(is%EARTH_INT_STATE)
    call ESMF_GridCompGetInternalState(driver, is, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
      
    ! deallocate internal state memory
    deallocate(is%EARTH_INT_STATE, stat=stat)
    if (ESMF_LogFoundDeallocError(statusToCheck=stat, &
      msg="Deallocation of internal state memory failed.", &
      line=__LINE__, file=trim(name)//":"//__FILE__, rcToReturn=rc)) &
      return  ! bail out
      
  end subroutine
      
  !-----------------------------------------------------------------------------
  
  recursive subroutine ModifyCplLists(driver, importState, exportState, clock, &
    rc)
    type(ESMF_GridComp)  :: driver
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    character(len=160)              :: name, msg
    type(ESMF_CplComp), pointer     :: connectorList(:)
    integer                         :: i, j, cplListSize
    character(len=160), allocatable :: cplList(:)
    character(len=160)              :: value
    type(WRAP_EARTH_INTERNAL_STATE) :: is

    rc = ESMF_SUCCESS
    
    ! query Component for this internal State
    nullify(is%EARTH_INT_STATE)
    call ESMF_GridCompGetInternalState(driver, is, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_LogWrite("Driver is in ModifyCplLists()", ESMF_LOGMSG_INFO, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    nullify(connectorList)
    call NUOPC_DriverGetComp(driver, compList=connectorList, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    write (msg,*) "Found ", size(connectorList), " Connectors."// &
      " Modifying CplList Attribute...."
    call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
      
    do i=1, size(connectorList)
      ! query Connector i for its name
      call ESMF_CplCompGet(connectorList(i), name=name, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      ! access CplList for Connector i
      call NUOPC_CompAttributeGet(connectorList(i), name="CplList", &
        itemCount=cplListSize, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      if (cplListSize>0) then
        allocate(cplList(cplListSize))
        call NUOPC_CompAttributeGet(connectorList(i), name="CplList", &
          valueList=cplList, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        ! go through all of the entries in the cplList and add options
        do j=1, cplListSize
          cplList(j) = trim(cplList(j))//":DumpWeights=true"
          cplList(j) = trim(cplList(j))//":SrcTermProcessing=1:TermOrder=SrcSeq"
          ! add connection options read in from configuration file
          call ESMF_AttributeGet(connectorList(i), name="ConnectionOptions", &
            value=value, defaultValue="", rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
          cplList(j) = trim(cplList(j))//trim(value)
        enddo
        ! store the modified cplList in CplList attribute of connector i
        call NUOPC_CompAttributeSet(connectorList(i), &
          name="CplList", valueList=cplList, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        deallocate(cplList)
      endif
    enddo
      
    deallocate(connectorList)
    
  end subroutine

  !-----------------------------------------------------------------------------

#ifdef CMEPS
  subroutine ReadAttributes(gcomp, config, label, relaxedflag, formatprint, rc)

    use ESMF  , only : ESMF_GridComp, ESMF_Config, ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use NUOPC , only : NUOPC_FreeFormatCreate, NUOPC_CompAttributeIngest
    use NUOPC , only : NUOPC_FreeFormatDestroy, NUOPC_FreeFormat

    ! input/output arguments
    type(ESMF_GridComp) , intent(inout)        :: gcomp
    type(ESMF_Config)   , intent(in)           :: config
    character(len=*)    , intent(in)           :: label
    logical             , intent(in), optional :: relaxedflag
    logical             , intent(in), optional :: formatprint
    integer             , intent(inout)        :: rc

    ! local variables
    type(NUOPC_FreeFormat)  :: attrFF
    character(len=*), parameter :: subname = "(module_EARTH_GRID_COMP.F90:ReadAttributes)"
    !-------------------------------------------

    rc = ESMF_SUCCESS

    if (present(relaxedflag)) then
       attrFF = NUOPC_FreeFormatCreate(config, label=trim(label), relaxedflag=.true., rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__)) &
         return  ! bail out
    else
       attrFF = NUOPC_FreeFormatCreate(config, label=trim(label), rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__)) &
         return  ! bail out
    end if

    call NUOPC_CompAttributeIngest(gcomp, attrFF, addFlag=.true., rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call NUOPC_FreeFormatDestroy(attrFF, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

  end subroutine ReadAttributes

  subroutine InitRestart(driver, rc)

    !-----------------------------------------------------
    ! Determine if will restart and read pointer file if appropriate
    !-----------------------------------------------------

    use ESMF         , only : ESMF_GridComp, ESMF_VM, ESMF_GridCompGet, ESMF_VMGet, ESMF_SUCCESS
    use ESMF         , only : ESMF_LogSetError, ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_RC_NOT_VALID
    use NUOPC        , only : NUOPC_CompAttributeGet, NUOPC_CompAttributeSet, NUOPC_CompAttributeAdd

    ! input/output variables
    type(ESMF_GridComp)    , intent(inout) :: driver
    integer                , intent(out)   :: rc

    ! local variables
    logical           :: read_restart   ! read the restart file, based on start_type
    character(len=ESMF_MAXSTR) :: cvalue         ! temporary
    character(len=ESMF_MAXSTR) :: rest_case_name ! Short case identification
    character(len=*) , parameter :: subname = "(module_EARTH_GRID_COMP.F90:InitRestart)"
    !-------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=rc)

    !-----------------------------------------------------
    ! Carry out restart if appropriate
    !-----------------------------------------------------

    read_restart = IsRestart(driver, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! Add rest_case_name and read_restart to driver attributes
    call NUOPC_CompAttributeAdd(driver, attrList=(/'rest_case_name','read_restart  '/), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    rest_case_name = ' '
    call NUOPC_CompAttributeSet(driver, name='rest_case_name', value=rest_case_name, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    write(cvalue,*) read_restart
    call NUOPC_CompAttributeSet(driver, name='read_restart', value=trim(cvalue), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

  end subroutine InitRestart

  function IsRestart(gcomp, rc)

    use ESMF         , only : ESMF_GridComp, ESMF_SUCCESS
    use ESMF         , only : ESMF_LogSetError, ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_RC_NOT_VALID
    use NUOPC        , only : NUOPC_CompAttributeGet

    ! input/output variables
    logical                                :: IsRestart
    type(ESMF_GridComp)    , intent(inout) :: gcomp
    integer                , intent(out)   :: rc

    ! locals
    character(len=ESMF_MAXSTR)            :: start_type     ! Type of startup
    character(len=ESMF_MAXSTR)            :: msgstr
    character(len=*) , parameter :: start_type_start = "startup"
    character(len=*) , parameter :: start_type_cont  = "continue"
    character(len=*) , parameter :: start_type_brnch = "branch"
    character(len=*) , parameter  :: subname = "(module_EARTH_GRID_COMP.F90:IsRestart)"
    !---------------------------------------

    rc = ESMF_SUCCESS

    ! First Determine if restart is read
    call NUOPC_CompAttributeGet(gcomp, name='start_type', value=start_type, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    if ((trim(start_type) /= start_type_start) .and.  &
        (trim(start_type) /= start_type_cont ) .and.  &
        (trim(start_type) /= start_type_brnch)) then
       write (msgstr, *) subname//': start_type invalid = '//trim(start_type)
       call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msgstr, line=__LINE__, file=__FILE__, rcToReturn=rc)
       return
    end if

    !TODO: this is hard-wired to CIME start/continue types in terms of gcomp
    IsRestart = .false.
    if (trim(start_type) == trim(start_type_cont) .or. trim(start_type) == trim(start_type_brnch)) then
       IsRestart = .true.
    end if

  end function IsRestart

  subroutine AddAttributes(gcomp, driver, config, compid, compname, inst_suffix, rc)

    ! Add specific set of attributes to components from driver attributes

    use ESMF  , only : ESMF_GridComp, ESMF_Config, ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use ESMF  , only : ESMF_LogFoundAllocError, ESMF_ConfigGetLen, ESMF_ConfigGetAttribute
    use NUOPC , only : NUOPC_CompAttributeAdd, NUOPC_CompAttributeGet, NUOPC_CompAttributeSet

    ! input/output parameters
    type(ESMF_GridComp) , intent(inout) :: gcomp
    type(ESMF_GridComp) , intent(in)    :: driver
    type(ESMF_Config)   , intent(inout) :: config
    integer             , intent(in)    :: compid
    character(len=*)    , intent(in)    :: compname
    character(len=*)    , intent(in)    :: inst_suffix
    integer             , intent(inout) :: rc

    ! local variables
    integer                        :: n
    integer                        :: stat
    integer                        :: inst_index
    logical                        :: is_present
    character(len=ESMF_MAXSTR)     :: cvalue
    character(len=32), allocatable :: compLabels(:)
    character(len=32), allocatable :: attrList(:)
    integer                        :: componentCount
    character(len=*), parameter    :: subname = "(module_EARTH_GRID_COMP.F90:AddAttributes)"
    logical                        :: lvalue = .false.
    !-------------------------------------------

    rc = ESMF_Success
    call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)

    !------
    ! Add compid to gcomp attributes
    !------
    !write(cvalue,*) compid
    !call NUOPC_CompAttributeAdd(gcomp, attrList=(/'MCTID'/), rc=rc)
    !if (chkerr(rc,__LINE__,u_FILE_u)) return
    !call NUOPC_CompAttributeSet(gcomp, name='MCTID', value=trim(cvalue), rc=rc)
    !if (chkerr(rc,__LINE__,u_FILE_u)) return

    !------
    ! Add all the other attributes in AttrList (which have already been added to driver attributes)
    !------
    !allocate(attrList(5))
    !attrList =  (/"read_restart", "orb_eccen   ", "orb_obliqr  ", "orb_lambm0  ", "orb_mvelpp  "/)
    ! TODO: orb_obliqr and orb_lambm0 not exist
    allocate(attrList(1))
    attrList =  (/"read_restart"/)

    call NUOPC_CompAttributeAdd(gcomp, attrList=attrList, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    do n = 1,size(attrList)
       if (trim(attrList(n)) == "read_restart") then
          call NUOPC_CompAttributeGet(driver, name="mediator_read_restart", value=cvalue, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out

          read(cvalue,*) lvalue

          if (.not. lvalue) then
            call NUOPC_CompAttributeGet(driver, name=trim(attrList(n)), value=cvalue, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, &
              file=__FILE__)) &
              return  ! bail out
          end if

          call NUOPC_CompAttributeSet(gcomp, name=trim(attrList(n)), value=trim(cvalue), rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
       else
          print*, trim(attrList(n))
          call NUOPC_CompAttributeGet(driver, name=trim(attrList(n)), value=cvalue, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
          call NUOPC_CompAttributeSet(gcomp, name=trim(attrList(n)), value=trim(cvalue), rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
       end if
    enddo
    deallocate(attrList)

    !------
    ! Add component specific attributes
    !------
    call ReadAttributes(gcomp, config, trim(compname)//"_attributes::", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ReadAttributes(gcomp, config, "ALLCOMP_attributes::", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    !------
    ! Add multi-instance specific attributes
    !------
    call NUOPC_CompAttributeAdd(gcomp, attrList=(/'inst_index'/), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! add inst_index attribute (inst_index is not required for cime internal components)
    ! for now hard-wire inst_index to 1
    inst_index = 1
    write(cvalue,*) inst_index
    call NUOPC_CompAttributeSet(gcomp, name='inst_index', value=trim(cvalue), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! add inst_suffix attribute
    if (len_trim(inst_suffix) > 0) then
       call NUOPC_CompAttributeAdd(gcomp, attrList=(/'inst_suffix'/), rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__)) &
         return  ! bail out
       call NUOPC_CompAttributeSet(gcomp, name='inst_suffix', value=inst_suffix, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__)) &
         return  ! bail out
    end if

  end subroutine AddAttributes
#endif
!
!-----------------------------------------------------------------------
!
      END MODULE module_EARTH_GRID_COMP
!
!-----------------------------------------------------------------------
