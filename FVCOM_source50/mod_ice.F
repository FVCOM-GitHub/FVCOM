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


! VERSIONS.  Finite volume 
!
! Comments:     Sea Ice Module 
! connect the thermodynamic module(CICE) and Finite-volume dynamics module(EVP)
!
! Current FVCOM dependency
!   archive.F:   - calls archive_ice
!   mod_ncdio.F: - netcdf output includes some variables
!   us_fvcom.F:  - main calls ice setup 
!=======================================================================


MODULE MOD_ICE 
#  if defined (ICE)
   USE ALL_VARS
   USE MOD_PAR

   USE MOD_UTILS
   use mod_ice2d
! afm 20180723 -- begin
#  if defined (HEATING_CALCULATED)
   USE MOD_HEATFLUX, ONLY : HEATING_FRESHWATER
#  endif
! --end

!====================================================================
!--------------------------------------------------------------
!   for cice code
! !USES:
!
      use ice_model_size
      use ice_albedo
      use ice_calendar
!      use ice_coupling
!      use ice_diagnostics
      use ice_domain
!      use ice_dyn_evp
      use ice_fileunits
      use ice_flux_in
      use ice_grid
!     use ice_history
      use ice_init
      use ice_itd
      use ice_itd_linear
      use ice_kinds_mod
      use ice_mechred
# if !defined (ICE_FRESHWATER)
! afm 20160513 & EJA 20160921 not used for freshwater
!      use ice_ocean
# endif
      use ice_atmo

      use ice_scaling
      use ice_therm_vertical
      use ice_therm_itd
      use ice_work
!====================================================================

    IMPLICIT NONE

    CONTAINS

!===============================================================================|
    SUBROUTINE ice_init_0
!
!==============================================================

  IMPLICIT NONE
  integer :: I,J,K
!--------------------------------------------------------------
!   for cice code
! !USES:
!=======================================================================
!=======================================================================|

!      ice_domain_def      
      IF(DBG_SET(DBG_SBR)) write(ipt,*)'START ice_init_0'
       
       ilo = 1  ! beg index of actual physical subdomain
       ihi = M !T  ! end index
       imt_local =  M !T !, &

       jlo = 1   ! beg index of actual physical subdomain
       jhi =1 
       jmt_local= 1

       allocate(index_global(ilo:ihi,jlo:jhi), rndex_global(ilo:ihi,jlo:jhi))

!       ice_work     
        allocate( work_l1(imt_local,jmt_local))
        allocate( work_l2(imt_local,jmt_local))
        allocate( worka(ilo:ihi,jlo:jhi))
        allocate( workb(ilo:ihi,jlo:jhi))

!      ice_state_def
 
      allocate(aice(imt_local,jmt_local))  ;aice=ZERO  ! concentration of ice 
      allocate(vice(imt_local,jmt_local))  ;vice=ZERO  ! volume per unit area of ice          (m)
      allocate(vsno(imt_local,jmt_local))  ;vsno=ZERO  ! volume per unit area of snow         (m)
      allocate(Tsfc(imt_local,jmt_local))  ;Tsfc=ZERO  ! temperature of ice/snow top surface  (C)
      allocate(eice(imt_local,jmt_local))  ;eice=ZERO  ! energy of melt. of ice           (J/m^2)
      allocate(esno(imt_local,jmt_local))  ;esno=ZERO  ! energy of melt. of snow layer    (J/m^2)
      allocate(tmass(imt_local,jmt_local)) ;tmass=ZERO ! total mass of ice and snow
      allocate(aice_init(imt_local,jmt_local));aice_init=ZERO ! concentration of ice at beginning of dt (for diagnostics)
! - begin afm 20171113 for ice open boundary -------------
      allocate(xflux_ice(0:MT));xflux_ice=ZERO
! - end afm 20171113 for ice open boundary -------------

      !-----------------------------------------------------------------
      ! state of the ice for each category
      !-----------------------------------------------------------------

      allocate(aicen(imt_local,jmt_local,ncat)) ;aicen=ZERO   ! concentration of ice 
      allocate(vicen(imt_local,jmt_local,ncat)) ;vicen=ZERO   ! volume per unit area of ice          (m)
      allocate(vsnon(imt_local,jmt_local,ncat)) ;vsnon=ZERO   ! volume per unit area of snow         (m)
      allocate(Tsfcn(imt_local,jmt_local,ncat)) ;Tsfcn=ZERO   ! temperature of ice/snow top surface  (C)
      allocate(esnon(imt_local,jmt_local,ncat)) ;esnon=ZERO   ! energy of melt. of snow layer    (J/m^2)

      allocate(aice0(imt_local,jmt_local))       ;aice0=C1i  ! concentration of open water
      allocate(eicen(imt_local,jmt_local,ntilay));eicen=ZERO  ! energy of melting for

      !-----------------------------------------------------------------
      ! other variables closely related to the state of the ice
      !-----------------------------------------------------------------
      !allocate(uvel(imt_local,jmt_local))      ;uvel=ZERO! x-component of velocity (m/s)
      !allocate(vvel(imt_local,jmt_local))      ;vvel=ZERO! y-component of velocity (m/s)
      allocate(strength(imt_local,jmt_local))  ;strength=ZERO   ! ice strength (N/m)

!=======================================================================
!      subroutine ice_albedo_def

      allocate(alvdrn (ilo:ihi,jlo:jhi,ncat)) ;alvdrn=ZERO ! visible, direct   (fraction)
      allocate(alidrn (ilo:ihi,jlo:jhi,ncat)) ;alidrn=ZERO ! near-ir, direct   (fraction)
      allocate(alvdfn (ilo:ihi,jlo:jhi,ncat)) ;alvdfn=ZERO ! visible, diffuse  (fraction)
      allocate(alidfn (ilo:ihi,jlo:jhi,ncat)) ;alidfn=ZERO ! near-ir, diffuse  (fraction)

      ! albedos aggregated over categories
      allocate(alvdr (ilo:ihi,jlo:jhi)) ;alvdr=ZERO  ! visible, direct   (fraction)
      allocate(alidr (ilo:ihi,jlo:jhi)) ;alidr=ZERO ! near-ir, direct   (fraction)
      allocate(alvdf (ilo:ihi,jlo:jhi)) ;alvdf=ZERO ! visible, diffuse  (fraction)
      allocate(alidf (ilo:ihi,jlo:jhi)) ;alidf=ZERO ! near-ir, diffuse  (fraction)

!=======================================================================
!      subroutine ice_flux_def
     !-----------------------------------------------------------------
      ! Dynamics component
      !-----------------------------------------------------------------
!   Dynamic component use FVCOM code

!      real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi)  :: &

       ! in from ocean
!       allocate(ss_tltx(ilo:ihi,jlo:jhi))    ;ss_tltx=ZERO   ! sea surface slope, x-direction (m/m)
!       allocate(ss_tlty(ilo:ihi,jlo:jhi))    ;ss_tlty=ZERO   ! sea surface slope, y-direction
!       allocate(uocn   (ilo:ihi,jlo:jhi))    ;uocn   =ZERO   ! ocean current, x-direction (m/s)
!       allocate(vocn   (ilo:ihi,jlo:jhi))    ;vocn   =ZERO   ! ocean current, y-direction (m/s)

       ! out to atmosphere

       allocate(strairxT(ilo:ihi,jlo:jhi))  ;strairxT=ZERO ! stress on ice by air, x-direction
       allocate(strairyT(ilo:ihi,jlo:jhi))  ;strairyT=ZERO ! stress on ice by air, y-direction
                     
       ! out to ocean          T-cell (kg/m s^2)
       allocate(strocnxT(ilo:ihi,jlo:jhi))  ;strocnxT=ZERO ! ice-ocean stress, x-direction
       allocate(strocnyT(ilo:ihi,jlo:jhi))  ;strocnyT=ZERO ! ice-ocean stress, y-direction

       ! diagnostic
       allocate(daidtd(ilo:ihi,jlo:jhi))  ;daidtd=ZERO ! ice area tendency due to transport   (s^-1)
       allocate(dvidtd(ilo:ihi,jlo:jhi))  ;dvidtd=ZERO ! ice volume tendency due to transport (m/s)

       ! internal              U-cell (kg/m s^2)
       allocate(strocnx(ilo:ihi,jlo:jhi))  ;strocnx=ZERO ! ice-ocean stress, x-direction
       allocate(strocny(ilo:ihi,jlo:jhi))  ;strocny=ZERO ! ice-ocean stress, y-direction
       allocate(strairx(ilo:ihi,jlo:jhi))  ;strairx=ZERO ! stress on ice by air, x-direction
       allocate(strairy(ilo:ihi,jlo:jhi))  ;strairy=ZERO ! stress on ice by air, y-direction
!       allocate(strtltx(ilo:ihi,jlo:jhi))  ;strtltx=ZERO ! stress due to sea surface slope, x-direction
!       allocate(strtlty(ilo:ihi,jlo:jhi))  ;strtlty=ZERO ! stress due to sea surface slope, y-direction
                                         
      !-----------------------------------------------------------------
      ! Thermodynamics component
      !-----------------------------------------------------------------

       ! in from atmosphere

!      real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi)  :: &
      allocate(zlvl (ilo:ihi,jlo:jhi))     ;zlvl =ZERO  ! atm level height (m)
      allocate(uatm (ilo:ihi,jlo:jhi))     ;uatm =ZERO  ! wind speed (m/s)
      allocate(vatm (ilo:ihi,jlo:jhi))     ;vatm =ZERO  !
      allocate(potT (ilo:ihi,jlo:jhi))     ;potT =ZERO  ! air potential temperature  (K)
      allocate(Tair (ilo:ihi,jlo:jhi))     ;Tair =ZERO  ! air temperature  (K)
      allocate(Qa   (ilo:ihi,jlo:jhi))     ;Qa   =ZERO  ! specific humidity (kg/kg)
      allocate(rhoa (ilo:ihi,jlo:jhi))     ;rhoa =1.3_SP! air density (kg/m^3)
      allocate(swvdr(ilo:ihi,jlo:jhi))     ;swvdr=ZERO  ! sw down, visible, direct  (W/m^2)
      allocate(swvdf(ilo:ihi,jlo:jhi))     ;swvdf=ZERO  ! sw down, visible, diffuse (W/m^2)
      allocate(swidr(ilo:ihi,jlo:jhi))     ;swidr=ZERO  ! sw down, near IR, direct  (W/m^2)
      allocate(swidf(ilo:ihi,jlo:jhi))     ;swidf=ZERO  ! sw down, near IR, diffuse (W/m^2)
      allocate(flw  (ilo:ihi,jlo:jhi))     ;flw  =ZERO  ! incoming longwave radiation (W/m^2)
      allocate(frain(ilo:ihi,jlo:jhi))     ;frain=ZERO  ! rainfall rate (kg/m^2 s)
      allocate(fsnow(ilo:ihi,jlo:jhi))     ;fsnow=ZERO  ! snowfall rate (kg/m^2 s)
                                                
       ! in from ocean                          
                                                
      allocate(frzmlt(ilo:ihi,jlo:jhi))   ;frzmlt  =ZERO  ! freezing/melting potential (W/m^2)
      allocate(sss   (ilo:ihi,jlo:jhi))   ;sss     =ZERO   ! sea surface salinity (ppt)
      allocate(sst   (ilo:ihi,jlo:jhi))   ;sst     =ZERO    ! sea surface temperature (C)
      allocate( Tf   (ilo:ihi,jlo:jhi))   ; Tf     =ZERO   ! freezing temperature (C)
      allocate( hmix (ilo:ihi,jlo:jhi))   ; hmix   =ZERO  ! mixed layer depth (m)
      allocate( qdp (ilo:ihi,jlo:jhi))    ; qdp    =ZERO  ! deep ocean heat flux (W/m^2)

       ! out to atmosphere 
       ! note albedos are in ice_albedo.F, Tsfc in ice_state.F

!      real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi)  :: &
       allocate( fsens (ilo:ihi,jlo:jhi))   ;fsens    =ZERO   ! sensible heat flux (W/m^2)
       allocate( flat  (ilo:ihi,jlo:jhi))   ;flat     =ZERO   ! latent heat flux   (W/m^2)
       allocate( fswabs(ilo:ihi,jlo:jhi))   ;fswabs   =ZERO  ! shortwave flux absorbed in ice and ocean (W/m^2)
       allocate( flwout(ilo:ihi,jlo:jhi))   ;flwout   =ZERO  ! outgoing longwave radiation (W/m^2)
       allocate( evap  (ilo:ihi,jlo:jhi))   ;evap     =ZERO ! evaporative water flux (kg/m^2/s)
       allocate( Tref  (ilo:ihi,jlo:jhi))   ;Tref     =ZERO ! 2m atm reference temperature (K)
       allocate( Qref  (ilo:ihi,jlo:jhi))   ;Qref     =ZERO  ! 2m atm reference sp humidity (kg/kg)
                                                     
       ! out to ocean

       allocate( fresh  (ilo:ihi,jlo:jhi))   ;fresh   =ZERO  ! fresh water flux to ocean (kg/m^2/s)
       allocate( fsalt  (ilo:ihi,jlo:jhi))   ;fsalt   =ZERO  ! salt flux to ocean (kg/m^2/s)
       allocate( fhnet  (ilo:ihi,jlo:jhi))   ;fhnet   =ZERO  ! net heat flux to ocean (W/m^2)
       allocate( fswthru  (ilo:ihi,jlo:jhi)) ;fswthru =ZERO   ! shortwave penetrating to ocean (W/m^2)
       ! afm 20150513 for outputs of heat fluxes.
!       allocate( fsh(ilo:ihi,jlo:jhi))   ;fsh =ZERO    ! sensible heat flux at water surface (W/m^2)
!       allocate( flh(ilo:ihi,jlo:jhi))   ;flh =ZERO    ! latent heat flux at water surface (W/m^2)
!       allocate( flwup(ilo:ihi,jlo:jhi)) ;flwup =ZERO  ! outgoing longwave radiation at water surface (W/m^2)
!       allocate( flwtot(ilo:ihi,jlo:jhi)) ;flwtot =ZERO  ! net longwave radiation at water surface (W/m^2)
! afm 20190807 - bounds fix
       allocate( fsh(0:MT,jlo:jhi))   ;fsh =ZERO    ! sensible heat flux at water surface (W/m^2)
       allocate( flh(0:MT,jlo:jhi))   ;flh =ZERO    ! latent heat flux at water surface (W/m^2)
       allocate( flwup(0:MT,jlo:jhi)) ;flwup =ZERO  ! outgoing longwave radiation at water surface (W/m^2)
       allocate( flwtot(0:MT,jlo:jhi)) ;flwtot =ZERO  ! net longwave radiation at water surface (W/m^2)

! afm 20170125 - hflx
!       allocate( flwitot(ilo:ihi,jlo:jhi)) ;flwitot =ZERO  ! net longwave radiation at ice surface (W/m^2)
! afm 20190807 - bounds fix
       allocate( flwitot(0:MT,jlo:jhi)) ;flwitot =ZERO  ! net longwave radiation at ice surface (W/m^2)

       ! diagnostic 

!      real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi)  :: &       
       allocate(congel(ilo:ihi,jlo:jhi))     ;congel =ZERO! basal ice growth         (m/step-->cm/day)
       allocate(frazil(ilo:ihi,jlo:jhi))     ;frazil =ZERO! frazil ice growth        (m/step-->cm/day)
       allocate(snoice(ilo:ihi,jlo:jhi))     ;snoice =ZERO! snow-ice formation       (m/step-->cm/day)
       allocate(meltt (ilo:ihi,jlo:jhi))     ;meltt  =ZERO! top ice melt             (m/step-->cm/day)
       allocate(meltb (ilo:ihi,jlo:jhi))     ;meltb  =ZERO! basal ice melt           (m/step-->cm/day)
       allocate(meltl (ilo:ihi,jlo:jhi))     ;meltl  =ZERO! lateral ice melt         (m/step-->cm/day)
       allocate(daidtt(ilo:ihi,jlo:jhi))     ;daidtt =ZERO! ice area tendency thermo.   (s^-1)
       allocate(dvidtt(ilo:ihi,jlo:jhi))     ;dvidtt =ZERO! ice volume tendency thermo. (m/s)
       allocate(mlt_onset(ilo:ihi,jlo:jhi))  ;mlt_onset =ZERO! day of year that sfc melting begins
       allocate(frz_onset(ilo:ihi,jlo:jhi))  ;frz_onset =ZERO! day of year that freezing begins (congel or frazil)

       ! NOTE: The following ocean diagnostic fluxes measure
       ! the same thing as their coupler counterparts but over
       ! different time intervals.  The coupler variables are 
       ! computed from one to_coupler call to the next, whereas
       ! the diagnostic variables are computed over a time step.
       allocate(fresh_hist  (ilo:ihi,jlo:jhi))  ;fresh_hist  =ZERO! fresh water flux to ocean (kg/m^2/s)
       allocate(fsalt_hist  (ilo:ihi,jlo:jhi))  ;fsalt_hist  =ZERO! salt flux to ocean (kg/m^2/s)
       allocate(fhnet_hist  (ilo:ihi,jlo:jhi))  ;fhnet_hist  =ZERO! net heat flux to ocean (W/m^2)
       allocate(fswthru_hist(ilo:ihi,jlo:jhi))  ;fswthru_hist=ZERO! shortwave penetrating to ocean (W/m^2)

       ! internal

!      real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi)  :: &

       allocate(fsw    (ilo:ihi,jlo:jhi))   ;fsw   =ZERO!    ! incoming shortwave radiation (W/m^2)
       allocate(wind   (ilo:ihi,jlo:jhi))   ;wind  =ZERO    ! wind speed (m/s)
       allocate(shcoef (ilo:ihi,jlo:jhi))   ;shcoef=ZERO  ! transfer coefficient for sensible heat
       allocate(lhcoef (ilo:ihi,jlo:jhi))   ;lhcoef=ZERO   ! transfer coefficient for latent heat

!c=======================================================================
!        subroutine ice_flux_in_def

      allocate(cldf(ilo:ihi,jlo:jhi)) 
!      real (kind=dbl_kind), dimension(ilo:ihi,jlo:jhi,2)  :: &
      !allocate(fsw_data  (ilo:ihi,jlo:jhi,2))  ;fsw_data =ZERO! field values at 2 temporal data points
      !allocate(cldf_data (ilo:ihi,jlo:jhi,2))  ;cldf_data=ZERO    
      !allocate(fsnow_data(ilo:ihi,jlo:jhi,2))  ;fsnow_data=ZERO     
      !allocate(Tair_data (ilo:ihi,jlo:jhi,2))  ;Tair_data=ZERO     
      !allocate(uatm_data (ilo:ihi,jlo:jhi,2))  ;uatm_data=ZERO    
      !allocate(vatm_data (ilo:ihi,jlo:jhi,2))  ;vatm_data=ZERO    
      !allocate( Qa_data  (ilo:ihi,jlo:jhi,2))  ; Qa_data =ZERO   
      !allocate(rhoa_data (ilo:ihi,jlo:jhi,2))  ;rhoa_data=ZERO    
      !allocate(potT_data (ilo:ihi,jlo:jhi,2))  ;potT_data=ZERO    
      !allocate(zlvl_data (ilo:ihi,jlo:jhi,2))  ;zlvl_data=ZERO    
      !allocate( flw_data (ilo:ihi,jlo:jhi,2))  ; flw_data=ZERO    
      !allocate(sst_data  (ilo:ihi,jlo:jhi,2))  ;sst_data =ZERO   
      !allocate(sss_data  (ilo:ihi,jlo:jhi,2))  ;sss_data =ZERO
                                               
!=======================================================================
!       subroutine ice_mechred_def
      
!      real (kind=dbl_kind), dimension(ilo:ihi,jlo:jhi) :: &

      allocate(dardg1dt(ilo:ihi,jlo:jhi)) ;dardg1dt=ZERO ! rate of fractional area loss by ridging ice (1/s)
      allocate(dardg2dt(ilo:ihi,jlo:jhi)) ;dardg2dt=ZERO ! rate of fractional area gain by new ridges (1/s)
      allocate(dvirdgdt(ilo:ihi,jlo:jhi)) ;dvirdgdt=ZERO ! rate of ice volume ridged (m/s)
      allocate(opening (ilo:ihi,jlo:jhi)) ;opening =ZERO ! rate of opening due to divergence/shear (1/s)

!      real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi)  :: &    
      allocate(asum (ilo:ihi,jlo:jhi)) ;asum =ZERO ! sum of total ice and open water area
      allocate(aksum(ilo:ihi,jlo:jhi)) ;aksum=ZERO ! ratio of area removed to area ridged

!      real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi,0:ncat) :: &
      allocate(apartic(ilo:ihi,jlo:jhi,0:ncat)) ; apartic =ZERO ! participation function; fraction of ridging/

!      real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi,ncat)  :: &    
      allocate(hrmin (ilo:ihi,jlo:jhi,ncat))  ;hrmin=ZERO ! minimum ridge thickness
      allocate(hrmax (ilo:ihi,jlo:jhi,ncat))  ;hrmax=ZERO ! maximum ridge thickness
      allocate(krdg  (ilo:ihi,jlo:jhi,ncat))  ;krdg =ZERO ! mean ridge thickness/thickness of ridging ice 
      allocate(hrexp (ilo:ihi,jlo:jhi,ncat))  ;hrexp=ZERO ! ridge e-folding thickness (krdg_redist = 1)

!=======================================================================
!      subroutine ice_therm_vertical_def
 
!      real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi,ncat)  ::  &
       allocate(hicen_old(ilo:ihi,jlo:jhi,ncat))    ;hicen_old =ZERO    ! old ice thickness (m), used by ice_therm_itd
!      real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi)  ::  &
       allocate(rside(ilo:ihi,jlo:jhi))  ;rside =ZERO ! fraction of ice that melts laterally
!=======================================================================
      !subroutine ice_therm_itd

!     real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi,ncat) :: &
      allocate(hicen(imt_local,jmt_local,ncat)) ;hicen=ZERO   !  ice thickness (m)

!=======================================================================
!      subroutine  ice_grid_def
      !allocate(TLAT_G(imt_global,jmt_global))     ;TLAT_G=ZERO! latitude of cell center T grid
      !allocate(TLON_G(imt_global,jmt_global))     ;TLON_G=ZERO! longitude of cell center T grid
!      real (kind=dbl_kind), dimension (imt_local,jmt_local)  ::  &
      !allocate(dxt(imt_local,jmt_local))    ;dxt=ZERO! width of T-cell through the middle (m)
      !allocate(dyt(imt_local,jmt_local))    ;dyt=ZERO! height of T-cell through the middle (m)
      !allocate(dxu(imt_local,jmt_local))    ;dxu=ZERO! width of U-cell through the middle (m)
      !allocate(dyu(imt_local,jmt_local))    ;dyu=ZERO! height of U-cell through the middle (m)
      !allocate(HTE(imt_local,jmt_local))    ;HTE=ZERO! length of eastern edge of T-cell (m)
      !allocate(HTN(imt_local,jmt_local))    ;HTN=ZERO! length of northern edge of T-cell (m)
      !allocate(HTS(imt_local,jmt_local))    ;HTS=ZERO! length of southern edge of T-cell
      !allocate(HTW(imt_local,jmt_local))    ;HTW=ZERO! length of western edge of T-cell
      !allocate(tarea(imt_local,jmt_local))  ;tarea=ZERO! area of T-cell (m^2)
      !allocate(uarea(imt_local,jmt_local))  ;uarea=ZERO! area of U-cell (m^2)
      allocate(ULON (imt_local,jmt_local))  ;ULON =ZERO! longitude of velocity pts (radians)
      allocate(ULAT (imt_local,jmt_local))  ;ULAT =ZERO! latitude of velocity pts (radians)
      allocate(TLON (imt_local,jmt_local))  ;TLON =ZERO!  longitude of temp pts (radians)
      allocate(TLAT (imt_local,jmt_local))  ;TLAT =ZERO! latitude of temp pts (radians)
      !allocate(dxhy (imt_local,jmt_local))  ;dxhy =ZERO! 0.5*(HTE - HTE)
      !allocate(dyhx (imt_local,jmt_local))  ;dyhx =ZERO! 0.5*(HTN - HTN)
                                             
!      real (kind=dbl_kind), dimension (ilo:;ihi,jlo:jhi)  ::  &
      !allocate(cyp (ilo:ihi,jlo:jhi))        ;cyp =ZERO! 1.5*HTE - 0.5*HTE
      !allocate(cxp (ilo:ihi,jlo:jhi))        ;cxp =ZERO! 1.5*HTN - 0.5*HTN
      !allocate(cym (ilo:ihi,jlo:jhi))        ;cym =ZERO! 0.5*HTE - 1.5*HTE
      !allocate(cxm (ilo:ihi,jlo:jhi))        ;cxm =ZERO! 0.5*HTN - 1.5*HTN
      !allocate(dxt2(ilo:ihi,jlo:jhi))        ;dxt2=ZERO! 0.5*dxt
      !allocate(dyt2(ilo:ihi,jlo:jhi))        ;dyt2=ZERO! 0.5*dyt
      !allocate(dxt4(ilo:ihi,jlo:jhi))        ;dxt4=ZERO! 0.25*dxt
      !allocate(dyt4(ilo:ihi,jlo:jhi))        ;dyt4=ZERO! 0.25*dyt
      !allocate(tarear  (ilo:ihi,jlo:jhi))    ;tarear  =ZERO! 1/tarea
      !allocate(uarear  (ilo:ihi,jlo:jhi))    ;uarear  =ZERO! 1/uarea
      !allocate(tinyarea(ilo:ihi,jlo:jhi))    ;tinyarea=ZERO! puny*tarea
      !allocate(ANGLE   (ilo:ihi,jlo:jhi))    ;ANGLE   =ZERO! for conversions between POP grid and lat/lon
      !allocate(ANGLET  (ilo:ihi,jlo:jhi))    ;ANGLET  =ZERO! ANGLE converted to T-cells
      !allocate(tarean  (ilo:ihi,jlo:jhi))    ;tarean  =ZERO! area of NH T-cells
      !allocate(tareas  (ilo:ihi,jlo:jhi))    ;tareas  =ZERO! area of SH T-cells

      ! Masks
!     real (kind=dbl_kind), dimension (imt_local,jmt_local)  ::  &
      allocate(hm    (imt_local,jmt_local))    ; hm     =zero! land/boundary mask, thickness (T-cell)
      !allocate(uvm   (imt_local,jmt_local))    ; uvm    =ZERO! land/boundary mask, velocity (U-cell)
      !allocate(mask_n(imt_local,jmt_local))    ; mask_n =ZERO! northern hemisphere
      !allocate(mask_s(imt_local,jmt_local))    ; mask_s =ZERO! southern hemisphere

!      logical (kind=log_kind)  ::  &
      allocate(tmask(imt_local,jmt_local))   ;tmask =.true.  ! land/boundary mask, thickness (T-cell)
      allocate(umask(imt_local,jmt_local))   ;umask =.true.  ! land/boundary mask, velocity (U-cell)
!      allocate(tmask(imt_local,jmt_local))  ;tmask =c1i  ! land/boundary mask, thickness (T-cell)
!      allocate(umask(imt_local,jmt_local))  ;umask =c1i  ! land/boundary mask, velocity (U-cell)
      !allocate(icetmask(ilo:ihi,jlo:jhi))   ;icetmask =ZERO! ice extent mask (T-cell)
      !allocate(iceumask(ilo:ihi,jlo:jhi))   ;iceumask =ZERO! ice extent mask (U-cell)

     !====================================================================

!       my_task,      &          ! task id for local process
!       master_task            ! task id for master process
        my_task=MYID    ! processor id
        master_task = 1  ! master processer id
!--------------------------------------------------------------
!====================================================================

      ! timing devices
!      call ice_timer_clear(-1) ! initialize all timers
!      call ice_timer_start(0)  ! begin timing entire run

      ! parameters, grid, parameterizations, etc
      call init_constants      ! pi, etc
      call input_data          ! namelist variables
      call init_grid           ! grid variables

!      if (advection == 'remap') 
!     &     call init_remap     ! more grid variables
!      call init_calendar       ! initialize some calendar stuff
!      call init_hist           ! initialize output history file
!      call init_evp            ! define evp dynamics parameters, variables
      call init_flux           ! initialize coupler fluxes
      call init_thermo_vertical ! initialize vertical thermodynamics
      call init_mechred        ! initialize ridging variables
      call init_itd            ! initialize ice thickness distribution

!     set the initial ice field base on the latitude
!      if (restart) call restartfile       ! start from restart data
!      call calendar(time)      ! determine the initial date

!#ifdef coupled
!      call init_cpl            ! initialize message passing 
!#else
      call init_getflux        ! initialize forcing data info
!#endif

       Do I=1,M !T
         Tair(I,1)=273.15+T_air(I)
       END DO

       call init_state          ! initialize the ice state

      dtice = dti
      IF(DBG_SET(DBG_LOG)) write(ipt,*)' ICE time step set to DTI dtice =' , dtice
      call ALLOC_UVICE

      call albedos             ! albedo based on initial ice distribution

      IF(DBG_SET(DBG_SBR)) write(ipt,*)'END ice_init_0'
      return
      end subroutine ice_init_0

!=======================================================================
! !IROUTINE: from_coupler - input from coupler to sea ice model
!
! !INTERFACE:
!
      subroutine from_coupler
!
!!DESCRIPTION:
!
! Reads input data from coupler to sea ice model
!

  IMPLICIT NONE
      integer :: i,j,k,icn,n2   ! local loop indices
      integer :: II,i1,I2,j1,j2   ! local loop indices
      real  ::   gsum, workx, worky
      REAL(SP), allocatable, dimension (:,:) :: rbuf

!     first step initialized the ice variables

      zlvl(:,:)    = C10I
      uatm(:,:)    = 0.0_SP
      vatm(:,:)    = 0.0_SP
      potT(:,:)    = 0.0_SP
      Tair(:,:)    = 0.0_SP
      Qa(:,:)      = 0.0_SP
      rhoa(:,:)    = 1.3_SP
      swvdr(:,:)   = 0.0_SP
      swvdf(:,:)   = 0.0_SP
      swidr(:,:)   = 0.0_SP
      swidf(:,:)   = 0.0_SP
      flw(:,:)     = 0.0_SP
      frain(:,:)   = 0.0_SP
      fsnow(:,:)   = 0.0_SP
      sst(:,:)     = 0.0_SP
      sss(:,:)     = 0.0_SP
      !uocn(:,:)    = 0.0_SP
      !vocn(:,:)    = 0.0_SP
      !ss_tltx(:,:) = 0.0_SP
      !ss_tlty(:,:) = 0.0_SP
      frzmlt(:,:)  = 0.0_SP


      do j=jlo,jhi
      do i=ilo,ihi
         !--- ocn states--
         sst  (i,j) = T1(I,1) !+Tffresh
         sss  (i,j) = S1(I,1)
         Tf   (i,j) = -depressT*sss(i,j) ! freezing temp (C)
         Tf   (i,j) = max(Tf(i,j),-1.85) ! freezing temp (C)

         !--- atm states-
         zlvl (i,j) = c10i
         Tair (i,j) = T_air(I) +Tffresh  !  ^0C--->"K"
         potT (i,j) = Tair(i,j)          ! K 
         Qa   (i,j) = QA_AIR(I)
         rhoa (i,j) = 1.3_SP !rhoair

         !--- ocn states--
         ! compute potential to freeze or melt ice

         frzmlt(i,j) = (Tf(i,j)-sst(i,j))*cp_ocn*rhow*(D(I))*DZ(I,1)/dtice
# if !defined (ICE_FRESHWATER)
         frzmlt(i,j) = min(max(frzmlt(i,j),-c1000),c1000)
# endif

         if (sst(i,j) <= Tf(i,j)) then
            sst(i,j) = Tf(i,j)

# if defined (ICE_FRESHWATER)
            ! afm 20151112  & EJA 20160921
            ! When supercooling (freezing, i.e. frzmlt>0), 
            ! T1 (temperature variable in hydrodynamic modules) is 
            ! set to Tf in subroutine to_coupler.
            ! The original method only counts surface layer k=1
            do k=2,kbm1
             IF(T1(I,K) .LT. Tf(I,j)) THEN
              frzmlt(i,j) = frzmlt(i,j) + (Tf(i,j)-T1(i,k))*cp_ocn*rhow*D(I)*DZ(I,k)/dtice
             END IF 
            end do
# else
            T1(I,1)=Tf(i,j)
# endif

          endif
       
       !--- atm fluxes--
       !!======================================================
# if defined (ICE_FRESHWATER)
! afm 20151001 & EJA 20160921 - For solar. For non-solar, shortwave from forcing data is used. 
          fsw(i,j)   = swrad_watts(i) 
# else
          fsw (i,j)  = DSW_AIR(I)
# endif
          cldf(i,j)  = CLOUD(I)
       !!======================================================

          fsnow(i,j)=0.0_SP
          fsnow(i,j) = (QPREC(I)+QEVAP(I))*RHOFRESH   !m/s---> kg/m^2/s
          frain(i,j) = c0i
        if (Tair(i,j) >= Tffresh) then
          frain(i,j) = fsnow(i,j)
          fsnow(i,j) = c0i
        endif

      end do
      end do
!!======================================================
         !!=============================================
         !! Special treatment the river mouth          !
         !! No new ice add at river mouth
         !! Caused by the salinity changing and the
         !! freezing potential increasing
         !!=============================================
     IF(RIVER_INFLOW_LOCATION == 'node') THEN
       IF(NUMQBC > 0) THEN
         DO J=1,NUMQBC
           I1=INODEQ(J)
           IF(frzmlt(I1,1)>0) frzmlt(I1,1)=0.0
         END DO
        END IF
     ELSE IF(RIVER_INFLOW_LOCATION == 'edge') THEN
       IF(NUMQBC > 0) THEN
         DO J=1,NUMQBC
           J1=N_ICELLQ(J,1)
           J2=N_ICELLQ(J,2)
           IF(frzmlt(J1,1)>0) frzmlt(J1,1)=0.0
           IF(frzmlt(J2,1)>0) frzmlt(J2,1)=0.0
         II=ICELLQ(J)
          DO I1=1,3
            IF(isonb(NV(II,I1))==0 ) then
              I2=NV(II,I1)
            END IF
          END DO
           IF(ISICEN(I2)==1 ) then
           !IF(frzmlt(J1,1)>0) 
            frzmlt(J1,1)=frzmlt(I2,1)
           !IF(frzmlt(J2,1)>0) 
            frzmlt(J2,1)=frzmlt(I2,1)
           END IF
         END DO
       END IF
     END IF

!!======================================================

      uatm =0.0_SP
      vatm =0.0_SP
      wind =0.0_SP
   DO I=1,M
      DO I1=1, NTVE(I)
         Uatm(I,1) =Uatm(I,1)+UUWIND(NBVE(I,I1))
         Vatm(I,1) =Vatm(I,1)+VVWIND(NBVE(I,I1))
         wind (i,1) =wind(i,1)+ sqrt(UUWIND(NBVE(I,I1))**2 + VVWIND(NBVE(I,I1))**2) ! wind speed, m/s
      END DO
      Uatm(I,1) = Uatm(I,1) /float(NTVE(I))
      Vatm(I,1) = Vatm(I,1) /float(NTVE(I))
      Wind(I,1) = sqrt(Uatm(I,1) * Uatm(I,1)+Vatm(I,1)*Vatm(I,1))
   END DO

#   if defined (MULTIPROCESSOR)
     allocate(rbuf(0:MT,3))
     rbuf =0.0
     rbuf(1:M,1)=Uatm(1:M,1)
     rbuf(1:M,2)=Vatm(1:M,1)
     rbuf(1:M,3)=Wind(1:M,1)

     if (PAR) CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,3,MYID,NPROCS,rbuf)
     IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,rbuf)

     Uatm(1:M,1)=rbuf(1:M,1)
     Vatm(1:M,1)=rbuf(1:M,2)
     Wind(1:M,1)=rbuf(1:M,3)
     deallocate(rbuf)
#   endif

!-----+----------------------------------------------------------------+
!    some parameter correction
!    Qa, Fsw , SW-->SWvdr SWVDF SWidr SWidf , FLW, frain, frain

       call prepare_forcing  !scale input forcing data

      end subroutine from_coupler

!=======================================================================
!BOP
!
! !IROUTINE: to_coupler - send data from sea ice model to coupler
!
! !INTERFACE:
!
      subroutine to_coupler
!
! !DESCRIPTION:
!
! Sea ice model to coupler
!
      use ice_flux
      integer(kind=int_kind) :: i,j,k,icn,n2     ! local loop indices
      real (kind=dbl_kind) ::    &
         gsum, workx, worky      &     ! tmps for converting grid
      ,  Tsrf (ilo:ihi,jlo:jhi)  &     ! surface temperature
      ,  tauxa(ilo:ihi,jlo:jhi)  &     ! atmo/ice stress
      ,  tauya(ilo:ihi,jlo:jhi)  &             
      ,  tauxo(ilo:ihi,jlo:jhi)  &     ! ice/ocean stress
      ,  tauyo(ilo:ihi,jlo:jhi)  &             
      ,  ailohi(ilo:ihi,jlo:jhi)      ! fractional ice area
      logical :: flag

!!=================================================================
      REAL(SP) :: SPRO,SPCP,ROSEA
! afm 20171115 for heat flux breakdown consistency-- begin
      REAL(SP) :: SPRO1(0:MT,jlo:jhi),ROSEA1(0:MT,jlo:jhi)
! afm 20171115 for heat flux breakdown consistency-- end
      REAL(SP), allocatable, dimension (:,:) :: rbuf
      REAL (SP) :: Naice   ! 

!!=================================================================
!--------ice_ocean---------------------------------------------------
      real (kind=dbl_kind), parameter ::  &
        cphm = cp_ocn*rhow*20 !hmix

      REAL (SP) :: DUVI,ice_ocnDX, ice_ocnDY  ! for ice-ocean drag 
!--------atmo-ocean---------------------------------------------------
      real (kind=dbl_kind), dimension(ilo:ihi,jlo:jhi) ::&
         delt     &! potential temperature difference   (K)
      ,  delq     &! specific humidity difference   (kg/kg)
      ,  dummy1, dummy2, dummy3, dummy4  ! dummy arrays

      real (kind=dbl_kind) :: &
         TsfK     &! surface temperature (K)
! afm for heat flux outputs
!      ,  fsh      &! sensible heat flux  (W/m^2)
!      ,  flh      &! latent heat flux    (W/m^2)
      ,  swabs    !&! surface absorbed shortwave heat flux (W/m^2)
!      ,  flwup     ! long-wave upward heat flux  (W/m^2)
!!====================================================================
     flag  =.false.
     SPCP  = 4.2174E3_SP
     ROSEA = 1.023E3_SP
     SPRO  = SPCP*ROSEA

! afm 20171115 for heat flux breakdown consistency-- begin
# if defined (HEATING_CALCULATED)
!     SPCP = 4.21734E3_SP
!     ROSEA = 1.023E3_SP
!     SPRO  = SPCP*ROSEA
!# endif
!# if defined (HEATING_CALCULATED_GL)
!     SPCP = 4.186E3_SP
!     ROSEA = 1.000E3_SP
!     SPRO  = SPCP*ROSEA
! afm 20180723 -- HEATING_CALCULATED_GL option is merged to HEATING_CALCULATED
         IF(.NOT. HEATING_FRESHWATER)THEN
           SPCP=4.2174E3_SP
           ROSEA = 1.023E3_SP
           SPRO = SPCP*ROSEA
         ELSE
           SPCP= 4.186E3_SP
           ROSEA = 1.000E3_SP
           SPRO = SPCP*ROSEA
         END IF
# endif
# if defined (HEATING_SOLAR)
     SPCP = 4.186E3_SP
     ROSEA1(:,1) = (-.0057_SP)*T1(:,1)*T1(:,1) + 0.0249_SP*T1(:,1) + 1000.0_SP
     SPRO1  = SPCP*ROSEA1
#endif
! afm 20171115 for heat flux breakdown consistency-- end

!  calculate the heat budget at the surface of ocean

!!====================================================================
!!====================================================================
!      if(.false.) then
! afm 20160622 & EJA 20160921 -- begin
# if defined (HEATING_CALCULATED) || (HEATING_SOLAR) 
      ! Use net heat flux (wtsurf) from COARE or SOLAR.
      ! Import latent, sensible heat fluxes and net longwave radiation
      ! for output purpose.
! afm 20190807 - bounds fix
      flh(1:M,jlo) = HLAT_WATTS(1:M)
      fsh(1:M,jlo) = HSENS_WATTS(1:M)
      flwtot(1:M,jlo) = LWRAD_WATTS(1:M)
# else
      ! Calculate net heat flux using CICE's algorithm. 
      ! Overwrite wtsurf.
      call atmo_boundary_layer (1, 'ocn', sst, &
           dummy1, dummy2, dummy3, dummy4, delt,  delq)
      do j = jlo,jhi
        do i = ilo,ihi
         ! ocean surface temperature in Kelvin
         TsfK = sst(i,j) + Tffresh

         ! shortwave radiative flux
# if defined (ICE_FRESHWATER)
! afm 20151112 & EJA 20160921 - albedo included in the solar subroutine or forcing file
         swabs =  fsw(i,j)
# else
         swabs = (c1i - albocn) * fsw(i,j)
# endif

         ! longwave radiative flux
         !flwup  = -stefan_boltzmann * TsfK**4
         ! afm for heat flux outputs
         flwup(i,j)  = -stefan_boltzmann * TsfK**4


         ! downward latent and sensible heat fluxes
         !flh = lhcoef(i,j) * delq(i,j)
         !fsh = shcoef(i,j) * delt(i,j)
! afm for heat flux outputs
         flh(i,j) = lhcoef(i,j) * delq(i,j)
         fsh(i,j) = shcoef(i,j) * delt(i,j)

# if defined (ICE_FRESHWATER)
! afm 20151112 & EJA 20160921 - Bug fix to include emissivity, assume similar to ice surface.
!         WTSURF(I) = -(fsh + flh + flwup + flw(i,j)+swabs ) /SPRO 
    !WTSURF(I) = -( fsh + flh + emissivity*( flwup + flw(i,j) ) + swabs ) / SPRO 
         ! for heat flux outputs 
         WTSURF(I) = -( fsh(i,j) + flh(i,j) + emissivity*( flwup(i,j) + flw(i,j)) + swabs ) / SPRO
! afm 20170125 - hflx
! emissivity included in the output
         FLWTOT(I,J) = emissivity * ( flwup(i,j) + flw(i,j) )
# else
         WTSURF(I) = -(fsh(i,j) + flh(i,j) + flwup(i,j) + flw(i,j)+swabs ) /SPRO 
# endif

# if !defined (ICE_FRESHWATER)
! afm 20151112 albedo included in the solar subroutine
         SWRAD(I)  =  (C1i-albocn)*SWRAD(I)
# endif
       end do
     end do


# endif
! afm 20160622 -- end 
!!====================================================================

      do j = jlo,jhi
      do i = ilo,ihi
       if (aice(i,j) > puny) then


! afm 20171115 for heat flux breakdown consistency-- begin
# if defined (HEATING_SOLAR)
        SPRO = SPRO1(i,1)
# endif
! afm 20171115 for heat flux breakdown consistency-- end

         !!==========================================================
         ! When the ocean is coupled to a sea-ice model, the temperature
         ! surface boundary condiction switches from a Neuman to a Dirichlet
         ! boundary condiction: the surface heat flux is no more specified,
         ! but diagnosed from the setting of the freezing temperature at the
         ! first ocean level
         ! the net heat including the second part, tha

           WTSURF(I) =  WTSURF(I)*(c1i-aice(I,1)) &
                       -fhnet(I,J) /SPRO & 
                       -fswthru(i,j)/SPRO 
! include wssurf correction due to 'fresh' here for salt water application.
        !  the short wave two parts
        !  first:  the open region
        !  second: penetrating through the ice to ocean
           SWRAD(I)  = SWRAD(I)*(c1i-aice(I,1))  -fswthru(i,j)/SPRO 

# if defined (HEATING_CALCULATED) || defined (HEATING_SOLAR)
! afm 20170125 - hflx
! ------ begin
flwitot(I,1) = flw(I,1)*emissivity + flwout(I,1)
HSENS_WATTS(I) = fsens(I,1)*aice(I,1) + fsh(I,1)*(c1i-aice(i,jlo))
HLAT_WATTS(I)  = flat(I,1)*aice(I,1) + flh(I,1)*(c1i-aice(i,jlo))
LWRAD_WATTS(I) = flwitot(I,1)*aice(I,1) + flwtot(I,1)*(c1i-aice(i,jlo))
! ------ end
# endif

        end if
        enddo                   ! i
      enddo

# if defined (HEATING_CALCULATED) || defined (HEATING_SOLAR)
! afm 20170125 - hflx
! ------ begin                      ! j
# if defined (MULTIPROCESSOR)
         IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,HSENS_WATTS,HLAT_WATTS,LWRAD_WATTS)
         IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,FLWITOT)
# endif
! ------ end
# endif



      DO I=1,M 
         IF(T1(I,1) .LT. Tf(I,1)) THEN
           T1(I,1) = Tf(I,1)
         END IF

         DO K=2,KBM1
         IF(T1(I,K) .LT. Tf(I,1)) THEN
           T1(I,K) = Tf(I,1)
         END IF
         END DO
      END DO
# if defined (MULTIPROCESSOR)
       IF(PAR) CALL AEXCHANGE(NC,MYID,NPROCS,T1)
# endif

      CALL N2E3D(T1,T)                               !Shift to Elements         !

!!======================================================
!!======================================================

# if defined (MULTIPROCESSOR)
       IF(PAR) CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,1,MYID,NPROCS,WTSURF,SWRAD)
       IF(PAR) CALL AEXCHANGE(NC,MYID,NPROCS,WTSURF,SWRAD)
# endif

! afm 20171115 for heat flux breakdown consistency-- begin
# if defined (HEATING_SOLAR)
      wtsurf_watts = -1.0_SP*wtsurf*SPRO1(:,1)
      swrad_watts = -1.0_SP*swrad *SPRO1(:,1)
# else
      wtsurf_watts = -1.0_SP*wtsurf*SPRO
      swrad_watts = -1.0_SP*swrad *SPRO
# endif

# if defined (MULTIPROCESSOR)
       IF(PAR) CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,1,MYID,NPROCS,WTSURF_WATTS,SWRAD_WATTS)
       IF(PAR) CALL AEXCHANGE(NC,MYID,NPROCS,WTSURF_WATTS,SWRAD_WATTS)
# endif
!! afm 20171115 -- end

       DO I=1,M
        IF(ISICEN(i)==1) THEN
! afm 20180801 moved QEVAP ice fractioning to the non ICE_FRESHWATER flag branch
!         QEVAP(I)=(1-aice(i,1))*QEVAP(I) !+EVAP(I,1)/RHOFRESH
# if !defined (ICE_FRESHWATER)
         QEVAP(I)=(1-aice(i,1))*QEVAP(I) !+EVAP(I,1)/RHOFRESH
# if !defined (SEMI_IMPLICIT)
         QPREC(I)=(1-aice(i,1))*QPREC(I)+(fresh(I,1)-evap(I,1))/rhofresh
# else
         QPREC(I)=(1-aice(i,1))*QPREC(I)-evap(I,1)/rhofresh   !!yding moved fresh out
         QPICE(I)=fresh(I,1)/rhofresh                         !!yding treat fresh seperatelly
# endif
# if defined (ICE_EMBEDDING)
         QFLICE(I) = fresh(I,1)/rhofresh     !!((RHO1(I,1)+1.0_SP)*1.0E3_SP)
         QTHICE(I)= vice(I,1)*(RHOI/((RHO1(I,1)+1.0_SP)*1.0E3_SP))   !!yding ICE THICKNESS IN THE WATER

         mlt_lat(i) = D(i)/QTHICE(i)           !!yding for lateral melting adjustment
# endif

# else
!         QPREC(I)=(1-aice(i,1))*QPREC(I)+(fresh(I,1)-evap(I,1))/rhofresh
! afm 20151112 & EJA 20160921 - Change in QPREC calculation.
! no mass flux due to ice melt/freezing 
! this causes large omega (WTS) at surface, resulting in unrealistic water
! temperature. Temperature flux due to ice formation/melting still are taken
! account of in wtsurf.
! Note: fresh & evap are already weighted by aicen, and fresh is a summation
!       of freshwater mass flux and evap, hence (fresh-evap) in the eq. above.

  ! EJA test comment out line below to prevent ice-fractioning of precip
! 20180801 afm incorporated the change in 3.2.5
!         QPREC(I)=(1-aice(i,1))*QPREC(I)
!         QEVAP(I)=(1-aice(i,1))*QEVAP(I) !+EVAP(I,1)/RHOFRESH

# endif

          ! evap  :: flux of vapor, atmos to ice   (kg m-2 s-1)
          !! this is the part lost from the ocean and ice system
          !          < 0 => sublimation
          !          > 0 => condensation
          ! Fresh :: flux of water, ice to ocean   (kg/m^2/s)
          !          negetive  ---> ice freezing
          !          positive  ---> ice melting

#        if !defined (SEMI_IMPLICIT)
         QPREC2(I) = QPREC(I)
         QEVAP2(I) = QEVAP(I)
!         QPICE2(I) = QPICE(I)    !!yding
!#   if defined (ICE_EMBEDDING)
!         QFLICE2(I)= QFLICE(I)   !!yding
!#   endif
#        endif

        END IF
       END DO

#   if defined (ICE_EMBEDDING)
         QTHICE_incr = QTHICE - QTHICERK    !!yding

         QTHICERK = QTHICE                  !!SAVE ICE THICKNESS FOR CURRENT STEP

#   endif


# if defined (MULTIPROCESSOR)
       if(PAR) CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,1,MYID,NPROCS,QPREC,QEVAP)
       IF(PAR) CALL AEXCHANGE(NC,MYID,NPROCS,QPREC,QEVAP)

# if !defined (ICE_FRESHWATER)
# if defined (SEMI_IMPLICIT)
!!yding
       if(PAR) CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,1,MYID,NPROCS,QPICE)
       IF(PAR) CALL AEXCHANGE(NC,MYID,NPROCS,QPICE)
!!yding
# endif

#   if defined (ICE_EMBEDDING)
       if(PAR) CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,1,MYID,NPROCS,QFLICE)
       IF(PAR) CALL AEXCHANGE(NC,MYID,NPROCS,QFLICE)

       if(PAR) CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,1,MYID,NPROCS,QTHICE)
       IF(PAR) CALL AEXCHANGE(NC,MYID,NPROCS,QTHICE)

       if(PAR) CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,1,MYID,NPROCS,QTHICERK)
       IF(PAR) CALL AEXCHANGE(NC,MYID,NPROCS,QTHICERK)

       if(PAR) CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,1,MYID,NPROCS,QTHICE_incr)
       IF(PAR) CALL AEXCHANGE(NC,MYID,NPROCS,QTHICE_incr)

       if(PAR) CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,1,MYID,NPROCS,mlt_lat)
       IF(PAR) CALL AEXCHANGE(NC,MYID,NPROCS,mlt_lat)

#   endif
#   endif

#      if !defined (SEMI_IMPLICIT)
       if(PAR) CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,1,MYID,NPROCS,QPREC2,QEVAP2)
       IF(PAR) CALL AEXCHANGE(NC,MYID,NPROCS,QPREC2,QEVAP2)

!!yding
!       if(PAR) CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,1,MYID,NPROCS,QPICE2)
!       IF(PAR) CALL AEXCHANGE(NC,MYID,NPROCS,QPICE2)
!!yding

#   if defined (ICE_EMBEDDING)
!       if(PAR) CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,1,MYID,NPROCS,QFLICE2)
!       IF(PAR) CALL AEXCHANGE(NC,MYID,NPROCS,QFLICE2)
#   endif

#      endif
# endif

      end subroutine to_coupler

!==============================================================
!==============================================================
        ! Main driver routine for FVCOM and UG-CICE and FVCOM
        !
        !
    subroutine ice_1D_2D

# if defined (MULTIPROCESSOR)
   USE MOD_PAR
# endif
   USE BCS
   USE MOD_OBCS
  
     IMPLICIT NONE
     integer :: i,ii,j,ij,K,Ni, JJ
     REAL(SP), allocatable, dimension (:,:) :: rbuf
     REAL(SP), allocatable, dimension (:) :: Irbuf
!--------------------------------------------------------------
!        for advection of ice variables
!--------------------------------------------------------------
      integer (kind=int_kind), parameter ::   &
        narr = 1 + 5*ncat     ! number of state variable arrays

      integer (kind=int_kind) :: &
       narrays          ! counter for number of state variable arrays

      real (kind=dbl_kind) :: works(imt_local,narr)  !works(imt_local,narr)
      real (kind=dbl_kind) :: worke(imt_local,ntilay)
      real (kind=dbl_kind) :: TTmin,TTmax
      integer :: iout
      integer  :: kice
!--------------------------------------------------------------
      !-----------------------------------------------------------------
      ! initializations
      !-----------------------------------------------------------------

      IF(dbg_set(dbg_SBR)) WRITE(IPT,*) "START ICE_1d_2d"

    ALLOCATE(AICETMP(0:MT));   AICETMP = 0.0_SP
    ALLOCATE(IMASS1(0:NT)) ;   IMASS1  = 0.0_SP
    ALLOCATE(IMASS(0:MT))  ;   IMASS   = 0.0_SP

    ISICEC_OLD = ISICEC  !  ice in last step
    ISICEC = 0
    ISICEN = 0
    IMASS  =0.0_SP
    IMASS1 =0.0_SP

     DO I=1,M !T
       IF (AICE(I,1) >p001) then
       IMASS(I) = (RHOI*VICE(I,1) + RHOS*VSNO(I,1)) ! kg/m^2
       IF (IMASS(I) >p01) ISICEN(I) =1
       END IF
     END DO

# if defined (MULTIPROCESSOR)
       IF(PAR) CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,1,MYID,NPROCS,IMASS)
       IF(PAR) CALL AEXCHANGE(NC,MYID,NPROCS,IMASS)
# endif

   CALL  N2E2D(IMASS,IMASS1)

   DO I=1,M 
     AICETMP(I) = AICE(I,1)
   END DO
   AIU =0.0_SP

# if defined (MULTIPROCESSOR)
       if(PAR) CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,1,MYID,NPROCS,AICETMP)
       IF(PAR) CALL AEXCHANGE(NC,MYID,NPROCS,AICETMP)
# endif

   CALL N2E2D(AICETMP,AIU)

   DO I=1,N 
      IF(sum(ISICEN(NV(I,1:3))) > 0 .and.IMASS1(I)>p01  &
         .and. AIU(I)>P001) &
      ISICEC(i)=1
   END DO

   Deallocate(AICETMP)
   DEALLOCATE(IMASS1,IMASS)

      !-----------------------------------------------------------------
      !-----------------------------------------------------------------

      call from_coupler        ! get updated info from flux couple

      Call thermo_vertical     ! thermodynamic growth rates and fluxes

      call to_coupler          ! collect/send data to flux coupler

      call thermo_itd          ! thermodynamics and associated itd changes

!!===========================================================================!!

#  if defined (ONE_D_MODEL)
      call zap_small_areas     ! remove categories with very small areas
      if (ncat > 1) call rebin ! make sure thicknesses are in bounds
      call aggregate           ! aggregate state variables over categories
      call albedos             ! compute ice albedos
      return                   ! NO advection for spinup run
# endif

      !!  DYNAMICS FOR SEA ICE
       DO KICE=1,NDYN_DT  !!  sub loop for ice model
       CALL ICE_RUNGE_KUTTA
!!===================================================================
      !-----------------------------------------------------------------
      ! Make sure state variables are not too close to zero.
      !-----------------------------------------------------------------

# if defined (ICE_FRESHWATER)
! afm 20160516 not needed.
! commented out as mass loss seems significant
!      call check_state
# else
      call check_state
# endif      

      !-----------------------------------------------------------------
      ! fill work arrays with fields to be advected
      !-----------------------------------------------------------------
      ! two arrays are used for performance (balance memory/cache vs
      ! number of bound calls);  one array or more than two may perform
      ! better depending on the machine used, number of processors, etc.
      ! --tested on SGI R2000, using 4 pes for the ice model under MPI
      !-----------------------------------------------------------------
      !!================================================================
 
      works =0.0_SP
      worke =0.0_SP

      do i = 1,imt_local
        do j = 1,jmt_local
         works(I,1) = min(max(aice0(i,j),0.0),1.0)
        enddo
      enddo

      narrays = 1
      do ni=1,ncat
       do i = 1,imt_local
         do j = 1,jmt_local
          ! this is the hice adjust the convergence to hice change
             works(I,narrays+1) = aicen(i,j,ni)  
             works(I,narrays+2) = vicen(i,j,ni) 
             works(I,narrays+3) = vsnon(i,j,ni) 
             works(I,narrays+4) = -aicen(i,j,ni)*Tsfcn(i,j,ni) 
             works(I,narrays+5) = -(esnon(i,j,ni)   &
                               +rhos*Lfresh*vsnon(i,j,ni))   ! for Tsnow<0
           enddo
         enddo
!         ! if there were 3 arrays in loop, use 3 instead and change narr
         narrays = narrays + 5
         enddo
!-----------------------------------------------------------------------------
  
     do k=1,ntilay
        do i = 1,imt_local
           do j = 1,jmt_local
           worke(I,k) = -eicen(i,j,k)
           enddo
       enddo
     enddo
      if (narr /= narrays)  &
     & write(ipt,*) "Wrong number of arrays in transport bound call"

!!--------------------------------------------------------------------------
     allocate(Irbuf(0:MT))
!!--------------------------------------------------------------------------

     Irbuf=0.0_SP
     DO K=1,narr
       Irbuf(1:M)=works(1:M,K)
# if defined (MULTIPROCESSOR)
       if (PAR) CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,1,MYID,NPROCS,Irbuf)
       IF (PAR) CALL AEXCHANGE(NC,MYID,NPROCS,Irbuf) 
#endif

       call ICE_ADV_fvcom(Irbuf,UICE2,VICE2)
! - begin afm 20171113 for ice open boundary -------------
! if vice, save XFLUX to XFLUX_ICE for the open boundary treatment.
! xflux_ice includes contributions from all the categories.
        IF ( MOD( K, 5 ) == 3 ) then
         DO I=1,IOBCN
          XFLUX_ICE_OBC(I) = XFLUX_ICE_OBC(I) + XFLUX_ICE(I_OBC_N(I))
         END DO
        END IF
! - end afm 20171113 for ice open boundary -------------

# if defined (MULTIPROCESSOR)
       if (PAR) CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,1,MYID,NPROCS,Irbuf)
       IF (PAR) CALL AEXCHANGE(NC,MYID,NPROCS,Irbuf) 
#endif
       works(1:M,K)=Irbuf(1:M)
     END DO

!!--------------------------------------------------------------------------

     Irbuf=0.0_SP
     do k=1,ntilay
     Irbuf(1:M)=worke(1:M,K)
# if defined (MULTIPROCESSOR)
       IF (PAR) CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,1,MYID,NPROCS,Irbuf)
       IF (PAR) CALL AEXCHANGE(NC,MYID,NPROCS,Irbuf) 
#endif

       call ICE_ADV_fvcom(Irbuf,UICE2,VICE2)
# if defined (MULTIPROCESSOR)
       IF (PAR) CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,1,MYID,NPROCS,Irbuf)
       IF (PAR) CALL AEXCHANGE(NC,MYID,NPROCS,Irbuf) 
#endif
      worke(1:M,K)=Irbuf(1:M)
      end do
      deallocate(Irbuf)

!!==========================================================================
      !-----------------------------------------------------------------
      ! retrieve advected fields from work array
      !-----------------------------------------------------------------
      do i = 1,imt_local
          do j = 1,jmt_local
          aice0(i,j) = min(max(works(I,1),C0i),c1i)
          enddo
      enddo

      narrays = 1               ! aice0 is first array
      do ni=1,ncat
        do i = 1,imt_local
           do j = 1,jmt_local
              aicen(i,j,ni) = min(max(works(I,narrays+1),C0i),c1i)
              vicen(i,j,ni) = works(I,narrays+2) 
             if ((aicen(i,j,ni)<puny .or.vicen(i,j,ni) < puny))then
                 aicen(i,j,ni)=0.0
                 vicen(i,j,ni)=0.0
             endif
              vsnon(i,j,ni) = works(I,narrays+3) 
             if (vsnon(i,j,ni) < puny) vsnon(i,j,ni) =0.0_SP

             if (aicen(i,j,ni) > puny) then
               Tsfcn(i,j,ni) = -works(i,narrays+4)/aicen(i,j,ni)
             else
                Tsfcn(i,j,ni) = Tf(i,j)
             endif
             if (aicen(i,j,ni) > puny) &
             esnon(i,j,ni) = -works(i,narrays+5)  &
     &                       -rhos*Lfresh*vsnon(i,j,ni) ! for Tsnow<0

           enddo
         enddo
         narrays = narrays + 5
      enddo

      J=1
      do i = 1,imt_local
          do ni=1,ncat
# if defined (ICE_FRESHWATER)
                Tsfcn(i,j,ni) = max( min( Tf(i,j), Tsfcn(i,j,ni)), T_air(i) - 10. )
# endif
           !IF(Tsfcn(i,j,ni) .le.T_air(I)) Tsfcn(i,j,ni) =T_air(i)
          enddo
      enddo

       !=============================================================================
      do k=1,ntilay
        do i=1,imt_local
         do j=1,jmt_local
             eicen(i,j,k) = -worke(i,k)
          enddo
         enddo
      enddo

150     continue
      call check_state

     call aggregate           ! aggregate state variables over categories

      J=1
      do i = 1,imt_local
         IF (aice(I,j) > 1.0)then
         do ni=1,ncat
             aicen(i,j,ni) = aicen(i,j,ni)/aice(i,j)
             vicen(i,j,ni) = vicen(i,j,ni)/aice(i,j)
             vsnon(i,j,ni) = vsnon(i,j,ni)/aice(i,j)
             esnon(i,j,ni) = esnon(i,j,ni)/aice(i,j)
         do k = 1, nilyr
             eicen(i,j,ilyr1(ni)+k-1) =eicen(i,j,ilyr1(ni)+k-1)/aice(i,j)
         end do

         enddo
         aice0 (i,j)=0.0
         aice(i,j)=1.0
         END IF
      enddo
!
     call ridge_ice(Delta,divu) ! ridging

      ENDDO  !! KTICE  sub loop for ice dynamics

      call zap_small_areas     ! remove categories with very small areas
      if (ncat > 1) call rebin ! make sure thicknesses are in bounds
      call aggregate           ! aggregate state variables over categories
      call albedos             ! compute ice albedos

! hhu ice bc 2/8/2017
      call bcond_gcn(4)
! hhu

      IF(dbg_set(dbg_SBR)) WRITE(IPT,*) "END ICE_1d_2d"

       return

       end subroutine ice_1D_2D

!=======================================================================!
#  endif

END MODULE MOD_ICE
