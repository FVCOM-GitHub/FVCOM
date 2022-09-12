!/===========================================================================/
! CVS VERSION INFORMATION
! $Id$
! $Name$
! $Revision$
!/===========================================================================/

!=========================================================================
!BOP
!
! !MODULE: ice_therm_vertical - thermo calculations before call to coupler
!
! !DESCRIPTION:
!
! Update ice and snow internal temperatures and compute
! thermodynamic growth rates and atmospheric fluxes.
!     
! See Bitz, C.M., and W.H. Lipscomb, 1999: 
! An energy-conserving thermodynamic model of sea ice,
! J. Geophys. Res., 104, 15,669-15,677.
!
! NOTE: The thermodynamic calculation is split in two for load balancing. 
!       First ice\_therm\_vertical computes vertical growth rates and coupler
!       fluxes.  Then ice\_therm\_itd does thermodynamic calculations not 
!       needed for coupling.
!
! !REVISION HISTORY:
!
! authors: William H. Lipscomb (LANL) and C.M. Bitz (UW)
!
! Vectorized by Clifford Chen (Fujitsu) and William H. Lipscomb (LANL)
!
! !INTERFACE:
!
      module ice_therm_vertical
!
! !USES:
! 
      use ice_model_size
      use ice_kinds_mod
      use ice_domain
      use ice_fileunits
      use ice_constants
      use ice_calendar
      use ice_grid
      use ice_state
      use ice_flux
      use ice_itd
      use mod_utils

#   if defined (ICE_EMBEDDING)
      use ALL_VARS, only : mlt_lat, QTHICE           !!yding
#   endif


!      use ice_diagnostics, only: print_state
!      use ice_diagnostics  ! print_state

!
!EOP
!
      implicit none
      save

      real (kind=dbl_kind), parameter ::   &
         ferrmax = 1.0e-3_dbl_kind & ! max  allowed energy flux error (W m-2)
                                   ! rec&ommend ferrmax < 0.01 W m-2
      ,  hsnomin = 1.0e-6_dbl_kind ! min  thickness for which Tsno computed (m)
                                         
      real (kind=dbl_kind) ::           &
         salin(nilyr+1)  &! salinity (ppt)
      ,  Tmlt(nilyr+1)   &! melting temp, -mu * salinity
      ,  ustar_scale     ! scaling for i ce-ocean heat flux
                                         
!      real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi,ncat) ::&
      real (kind=dbl_kind), dimension (:,:,:),allocatable,save  ::&
         hicen_old       ! old ice thick ness (m), used by ice_therm_itd
                                         
!      real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi) ::&
      real (kind=dbl_kind), dimension (:,:) ,allocatable,save ::&
         rside           ! fraction of i ce that melts laterally
                                         
      logical (kind=log_kind) ::        &
         l_brine         ! if true, treat brine pocket effects

      character (char_len), private :: stoplabel

      real (kind=dbl_kind) ::  i0vis,      floediam 
!=======================================================================

      contains

!=======================================================================
!BOP
!
! !ROUTINE: init_thermo_vertical - initialize salinity and melting temp
!
! !DESCRIPTION:
!
! Initialize the vertical profile of ice salinity and melting temperature.
!
! !REVISION HISTORY:
!
! authors: C. M. Bitz, UW
!          William H. Lipscomb, LANL
!
! !INTERFACE:
!
      subroutine init_thermo_vertical
!
! !USES:
! 
      use ice_itd
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      real (kind=dbl_kind), parameter ::&
         nsal      = 0.407_dbl_kind&
      ,  msal      = 0.573_dbl_kind&
# if defined (ICE_FRESHWATER)
! afm 20151112 & EJA 20160921
      ,  saltmax   = 0.0_dbl_kind  &! max salinity at ice base (ppt)
# else
      ,  saltmax   = 3.2_dbl_kind  &! max salinity at ice base (ppt)
# endif

      ,  min_salin = 0.1_dbl_kind  ! threshold for brine pocket treatment 

      integer (kind=int_kind) :: k        ! ice layer index
      real (kind=dbl_kind)    :: zn       ! normalized ice thickness

      !-----------------------------------------------------------------
      ! Determine l_brine based on saltmax.
      ! Thermodynamic solver will not converge if l_brine is true and
      !  saltmax is close to zero.
      !-----------------------------------------------------------------

      if (saltmax > min_salin) then
         l_brine = .true.
      else
         l_brine = .false.
      endif

      ! salinity and melting temperature profile
      if (l_brine) then
         do k = 1, nilyr
            zn = (real(k,kind=dbl_kind)-p5) / real(nilyr,kind=dbl_kind)
            salin(k)=(saltmax/c2i)*(c1i-cos(pi*zn**(nsal/(msal+zn))))
!            salin(k)=saltmax ! for isosaline ice
         enddo
         salin(nilyr+1) = saltmax
         do k = 1, nilyr+1
            Tmlt(k) = -salin(k)*depressT
         enddo
      else
         do k = 1, nilyr+1
            salin(k) = c0i
            Tmlt(k) = c0i
         enddo
      endif
      ustar_scale = c1i           ! for nonzero currents

      end subroutine init_thermo_vertical

!=======================================================================
!BOP
!
! !ROUTINE: thermo_vertical - driver for pre-coupler thermodynamics
!
! !DESCRIPTION:
!
! Driver for updating ice and snow internal temperatures and 
! computing thermodynamic growth rates and atmospheric fluxes.
!
! NOTE: The wind stress is computed here for later use.
!
! !REVISION HISTORY:
!
! authors: William H. Lipscomb, LANL
!          C. M. Bitz, UW
!          
! !INTERFACE:
!
      subroutine thermo_vertical
!
! !USES:
! 
      use ice_atmo
!     use ice_timers
!      use ice_ocean
!      ggao 
      use ice_work, only: worka, workb

!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) ::   & 
         i, j            &! horizontal indices
      ,  ij              &! horizontal index, combines i and j loops
      ,  ni               &! thickness category index
      ,  k               &! ice layer index
      ,  icells           ! number of cells with aicen > puny

      integer (kind=int_kind), dimension(1:(ihi-ilo+1)*(jhi-jlo+1)) ::   & 
        indxi, indxj    ! compressed indices for cells with aicen > puny

      real (kind=dbl_kind) ::   &
         dhi            & ! change in ice thickness
      ,  dhs             ! change in snow thickness

! 2D coupler variables (computed for each category, then aggregated) 

      real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi) ::   &
         strxn          & ! air/ice zonal  strss,              (N/m^2)
      ,  stryn          & ! air/ice merdnl strss,              (N/m^2)
      ,  fsensn         & ! surface downward sensible heat     (W/m^2)
      ,  flatn          & ! surface downward latent heat       (W/m^2)
      ,  fswabsn        & ! shortwave absorbed by ice          (W/m^2)
      ,  flwoutn        & ! upward LW at surface               (W/m^2)
      ,  evapn          & ! flux of vapor, atmos to ice   (kg m-2 s-1)
      ,  Trefn          & ! air tmp reference level                (K)
      ,  Qrefn          & ! air sp hum reference level         (kg/kg)
      ,  freshn         & ! flux of water, ice to ocean     (kg/m^2/s)
      ,  fsaltn         & ! flux of salt, ice to ocean      (kg/m^2/s)
      ,  fhnetn         & ! fbot corrected for leftover energy (W/m^2)
      ,  fswthrun         ! SW through ice to ocean            (W/m^2)

! 2D state variables (thickness, temperature, enthalpy)

      real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi) ::   &
         hin            & ! ice thickness (m)
      ,  hsn            & ! snow thickness (m)
      ,  hlyr           & ! ice layer thickness
      ,  hin_init       & ! initial value of hin
      ,  hsn_init       & ! initial value of hsn
      ,  hsn_new        & ! thickness of new snow (m) 
      ,  qsn            & ! snow enthalpy, qsn < 0 (J m-3) 
      ,  Tsn            & ! internal snow temperature (deg C)
      ,  Tsf              ! ice/snow top surface temp, same as Tsfcn (deg C)

      real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi,nilyr) ::   &
         hnew           & ! new ice layer thicknesses (m)
      ,  qin            & ! ice layer enthalpy, qin < 0 (J m-3) 
      ,  Tin             ! internal ice layer temperatures

! other 2D flux and energy variables

      real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi) ::   &
         Tbot           & ! ice bottom surface temperature (deg C)
      ,  fbot           & ! ice-ocean heat flux at bottom surface (W/m^2)
      ,  fsurf          & ! net flux to top surface, not including fcondtop
      ,  fcondtop       & ! downward cond flux at top surface (W m-2)
      ,  fcondbot       & ! downward cond flux at bottom surface (W m-2)
      ,  fswint         & ! SW absorbed in ice interior, below surface (W m-2)
      ,  einit          & ! initial energy of melting (J m-2)
      ,  efinal         & ! final energy of melting (J m-2)
      ,  mvap            ! ice/snow mass sublimated/condensed (kg m-2)

!      call ice_timer_start(4)   ! column model
!      call ice_timer_start(5)   ! thermodynamics

      !-----------------------------------------------------------------
      ! Initialize coupler variables sent to the atmosphere. 
      !
      ! We cannot initialize the ocean coupler fluxes (fresh, fsalt,
      ! fhnet, and fswthru) because these may have changed during the
      ! previous time step after the call to_coupler.  The ocean 
      ! coupler fluxes are initialized by init_flux_ocn after
      ! the call to_coupler.
      !-----------------------------------------------------------------
      call init_flux_atm

      !-----------------------------------------------------------------
      ! Initialize diagnostic variables for history files.
      !-----------------------------------------------------------------
      call init_diagnostics

      !-----------------------------------------------------------------
      ! Adjust frzmlt to account for ice-ocean heat fluxes since last 
      !  call to coupler.
      ! Compute lateral and bottom heat fluxes.
      !-----------------------------------------------------------------
      call frzmlt_bottom_lateral (Tbot, fbot, rside)

      !-----------------------------------------------------------------
      ! Vertical thermodynamics: growth rates, updated ice state, and
      ! coupler variables
      !-----------------------------------------------------------------

      do ni = 1, ncat

      !-----------------------------------------------------------------
      ! Save initial ice thickness (used later for ITD remapping)
      !-----------------------------------------------------------------

         do j = jlo, jhi
         do i = ilo, ihi
            if (aicen(i,j,ni) > puny) then
               hicen_old(i,j,ni) = vicen(i,j,ni) / aicen(i,j,ni)
            else
               hicen_old(i,j,ni) = c0i
            endif
         enddo                  ! i
         enddo                  ! j


      !-----------------------------------------------------------------
      ! Initialize the single-category versions of coupler fluxes,
      ! along with other thermodynamic arrays.
      !-----------------------------------------------------------------

         call init_thermo_vars (strxn,    stryn,    Trefn,    Qrefn,      &
                                fsensn,   flatn,    fswabsn,  flwoutn,    &
                                evapn,    freshn,   fsaltn,   fhnetn,     &
                                fswthrun, fsurf,    fcondtop, fcondbot,   &
                                fswint,   einit,    efinal,   mvap)

      !-----------------------------------------------------------------
      ! prep for air-to-ice heat, momentum, radiative and water fluxes 
      !-----------------------------------------------------------------

         call atmo_boundary_layer (ni, 'ice', Tsfcn(ilo:ihi,jlo:jhi,ni),&
                          strxn, stryn, Trefn, Qrefn, worka, workb)      
                                                                         
      !----------------------------------------------------------------- 
      ! Identify cells with nonzero ice area                             
      !----------------------------------------------------------------- 
                                                                         
         icells = 0                                                      
         do j = jlo, jhi                                                 
         do i = ilo, ihi                                                 
            if (aicen(i,j,ni) > puny) then                               
               icells = icells + 1                                       
               indxi(icells) = i                                         
               indxj(icells) = j                                         
            endif                                                        
         enddo                  ! i                                      
         enddo                  ! j                                      
                                                                         
      !----------------------------------------------------------------- 
      ! Compute variables needed for vertical thermo calculation         
      !----------------------------------------------------------------- 
                                                                         
         call init_vertical_profile                                     &
             (ni,    icells, indxi, indxj,                              &
              hin,  hlyr,   hsn,   hin_init, hsn_init,                  &
              qin,  Tin,    qsn,   Tsn,      Tsf,   einit)

      !-----------------------------------------------------------------
      ! Compute new surface temperature and internal ice and snow 
      !  temperatures  
      !-----------------------------------------------------------------

         call temperature_changes                                    &
             (ni,       icells,   indxi,  indxj,                     &
              hlyr,    hsn,      qin,    Tin,    qsn,      Tsn,      &
              Tsf,     Tbot,     fbot,   fsensn, flatn,    fswabsn,  &
              flwoutn, fswthrun, fhnetn, fsurf,  fcondtop, fcondbot, &
              fswint,  einit)

      !-----------------------------------------------------------------
      ! Compute growth and/or melting at the top and bottom surfaces.
      ! Repartition ice into equal-thickness layers, conserving energy.
      !-----------------------------------------------------------------

         call thickness_changes                           &
             (ni,     icells,   indxi,    indxj,          &
              hin,   hlyr,     hsn,      qin,    qsn,     &
              mvap,  Tbot,     fbot,     flatn,  fhnetn,  &
              fsurf, fcondtop, fcondbot, efinal) 

      !-----------------------------------------------------------------
      ! Check for energy conservation by comparing the change in energy 
      ! to the net energy input
      !-----------------------------------------------------------------

         call conservation_check_vthermo                              &
             (ni,     icells, indxi,  indxj,                           &
              fsurf, flatn,  fhnetn, fswint, einit, efinal)            
                                                                       
      !-----------------------------------------------------------------
      ! Let it snow                                                    
      !-----------------------------------------------------------------
                                                                       
         call add_new_snow                                            &
             (ni,    icells,  indxi, indxj,                            &
              hsn,  hsn_new, qsn,   Tsf)                               
                                                                       
      !-----------------------------------------------------------------
      ! Compute fluxes of water and salt from ice to ocean             
      !-----------------------------------------------------------------
                                                                       
         do ij = 1, icells                                             
            i = indxi(ij)                                              
            j = indxj(ij)                                              
                                                                       
            dhi = hin(i,j) - hin_init(i,j)                             
            dhs = hsn(i,j) - hsn_init(i,j)                             
                                                                       
!            evapn(i,j)  = mvap(i,j)/dt  ! mvap < 0 => sublimation,     
            evapn(i,j)  = mvap(i,j)/dtice  ! mvap < 0 => sublimation,     
                                        ! mvap > 0 => condensation     
            freshn(i,j) = evapn(i,j) -                                &
                 (rhoi*dhi + rhos*(dhs-hsn_new(i,j))) / dtice 
            fsaltn(i,j) = -rhoi*dhi*ice_ref_salinity*p001/dtice
!                 (rhoi*dhi + rhos*(dhs-hsn_new(i,j))) / dt
!            fsaltn(i,j) = -rhoi*dhi*ice_ref_salinity*p001/dt

         enddo                  ! ij

      !-----------------------------------------------------------------
      ! Increment area-weighted fluxes.
      ! Note: Must not change aicen before calling this subroutine.
      !-----------------------------------------------------------------

         call merge_fluxes (ni,                                          &
           strxn,   stryn,  fsensn, flatn, fswabsn,                     &
           flwoutn, evapn,  Trefn,  Qrefn,                              &
           freshn,  fsaltn, fhnetn, fswthrun)                            

      !----------------------------------------------------------------- 
      !  Given the vertical thermo state variables (hin, hsn, qin,       
      !   qsn, compute the new ice state variables (vicen, vsnon,        
      !   eicen, esnon).                                                 
      !----------------------------------------------------------------- 
                                                                         
         call update_state_vthermo                                      &
             (ni,    icells, indxi, indxj,                               &
              hin,  hsn,    qin,   qsn,  Tsf)

      enddo                     ! ncat

      !-----------------------------------------------------------------
      ! Update mixed layer with heat from ice 
      !-----------------------------------------------------------------

!      if (oceanmixed_ice) then
!         do j=jlo,jhi
!         do i=ilo,ihi
!          if (hmix(i,j) > c0i) sst(i,j) = sst(i,j)                    &
!               + (fhnet(i,j)+fswthru(i,j))*dt/(cp_ocn*rhow*hmix(i,j))
!         enddo
!         enddo
!      endif

!      call ice_timer_stop(5)    ! thermodynamics
!      call ice_timer_stop(4)    ! column model
!      ggao

      end subroutine thermo_vertical

!=======================================================================
!BOP
!
! !ROUTINE: frzmlt_bottom_lateral - bottom and lateral heat fluxes
!
! !DESCRIPTION:
!
! Adjust frzmlt to account for changes in fhnet since from\_coupler.
! Compute heat flux to bottom surface.
! Compute fraction of ice that melts laterally.
!
! !REVISION HISTORY:
!
! authors C. M. Bitz, UW
!         William H. Lipscomb, LANL
!         Elizabeth C. Hunke, LANL
!
! !INTERFACE:
!
      subroutine frzmlt_bottom_lateral (Tbot, fbot, rside)
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind), dimension(ilo:ihi,jlo:jhi), intent(out) ::   &
         Tbot       & ! ice bottom surface temperature (deg C)
      ,  fbot       & ! heat flux to ice bottom  (W/m^2)
      ,  rside       ! fraction of ice that melts laterally
!
!EOP
!
      integer (kind=int_kind) ::   &
         i, j       & ! horizontal indices
      ,  ni          & ! thickness category index
      ,  k          & ! layer index
      ,  ij         & ! horizontal index, combines i and j loops
      ,  icells      ! number of cells with ice melting

      integer (kind=int_kind), dimension(1:(ihi-ilo+1)*(jhi-jlo+1)) ::   & 
         indxi, indxj    ! compressed indices for cells with ice melting

      real (kind=dbl_kind), dimension(ilo:ihi,jlo:jhi) ::   &
         etot       & ! total energy in column
      ,  fside       ! lateral heat flux (W/m^2)     

      real (kind=dbl_kind) ::   &
         deltaT      &! SST - Tbot >= 0
      ,  ustar       &! skin friction velocity for fbot (m/s)
      ,  wlat        &! lateral melt rate (m/s)
      ,  xtmp        ! temporarty variable

      ! Parameters for bottom melting
 
      ! 0.006 = unitless param for basal heat flx ala McPhee and Maykut

      real (kind=dbl_kind), parameter ::   & 
         cpchr = -cp_ocn*rhow*0.006_dbl_kind  &
      ,  ustar_min = 5.e-3_dbl_kind

      ! Parameters for lateral melting 

      real (kind=dbl_kind), parameter ::   &
!         floediam = 300.0_dbl_kind   & ! effective floe diameter (m)
        alpha    = 0.66_dbl_kind    & ! constant from Steele (unitless)
      ,  m1 = 1.6e-6_dbl_kind        & ! constant from Maykut & Perovich
                                      ! (m/s/deg^(-m2))
      ,  m2 = 1.36_dbl_kind           ! constant from Maykut & Perovich
                                      ! (unitless)

      !-----------------------------------------------------------------
      ! Set bottom ice surface temperature and initialize fluxes.
      !-----------------------------------------------------------------


      do j = jlo, jhi
      do i = ilo, ihi
         Tbot(i,j)  = Tf(i,j)
         fbot(i,j)  = c0i
         rside(i,j) = c0i
         fside(i,j) = c0i
         etot(i,j)  = c0i
      enddo                     ! i
      enddo                     ! j

      !-----------------------------------------------------------------
      ! Identify grid cells where ice can melt.
      !-----------------------------------------------------------------

      icells = 0
      do j = jlo, jhi
      do i = ilo, ihi
         if (aice(i,j) > puny .and. frzmlt(i,j) < c0i) then ! ice can melt
            icells = icells + 1
            indxi(icells) = i
            indxj(icells) = j
         endif
      enddo                     ! i
      enddo                     ! j

      do ij = 1, icells  ! cells where ice can melt
         i = indxi(ij)
         j = indxj(ij)

      !-----------------------------------------------------------------
      ! Use boundary layer theory for fbot.
      ! See Maykut and McPhee (1995): JGR, 100, 24,691-24,703.
      !-----------------------------------------------------------------

         deltaT = max((sst(i,j)-Tbot(i,j)),c0i) 

         ! strocnx has units N/m^2 so strocnx/rho has units m^2/s^2 
         ustar = sqrt (sqrt(strocnxT(i,j)**2+strocnyT(i,j)**2)/rhow)
         ustar = max (ustar,ustar_min*ustar_scale)

         fbot(i,j) = cpchr * deltaT * ustar ! < 0
#   if defined (ICE_EMBEDDING)
         fbot(i,j) = fbot(i,j) + cpchr * deltaT * ustar * QTHICE(i)     
#   endif

         fbot(i,j) = max (fbot(i,j), frzmlt(i,j)) ! frzmlt < fbot < 0

!!! uncomment to use all frzmlt for standalone runs
!         fbot(i,j) = min (c0i, frzmlt(i,j))  

      !-----------------------------------------------------------------
      ! Compute rside.  See these references:
      !    Maykut and Perovich (1987): JGR, 92, 7032-7044
      !    Steele (1992): JGR, 97, 17,729-17,738  
      !-----------------------------------------------------------------
            
         wlat = m1 * deltaT**m2 ! Maykut & Perovich
!         rside(i,j) = wlat*dt*pi/(alpha*floediam) ! Steele 
         rside(i,j) = wlat*dtice*pi/(alpha*floediam) ! Steele

#   if defined (ICE_EMBEDDING)
         if (rside(i,j)<0.5) then                
          rside(i,j) = rside(i,j)*mlt_lat(i)     
         else if (rside(i,j)>=0.5) then
          rside(i,j) = rside(i,j)*QTHICE(i)    
         end if
#   endif
 
         rside(i,j) = max(c0i,min(rside(i,j),c1i))

      enddo                     ! ij

      !-----------------------------------------------------------------
      ! Compute heat flux associated with this value of rside.
      !-----------------------------------------------------------------

      do ni = 1, ncat

         ! melting energy/unit area in each column, etot < 0
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)
            etot(i,j) = esnon(i,j,ni)
         enddo                  ! ij

         do k = 1, nilyr         
            do ij = 1, icells
               i = indxi(ij)
               j = indxj(ij)
               etot(i,j) = etot(i,j) + eicen(i,j,ilyr1(ni)+k-1) 
            enddo               ! ij
         enddo                  ! k         

         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)
            ! lateral heat flux
!            fside(i,j) = fside(i,j) + rside(i,j)*etot(i,j)/dt ! fside < 0
            fside(i,j) = fside(i,j) + rside(i,j)*etot(i,j)/dtice ! fside < 0
         enddo                  ! ij

      enddo                     ! n

      !-----------------------------------------------------------------
      ! Limit bottom and lateral heat fluxes if necessary.
      !-----------------------------------------------------------------


      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)

         xtmp = frzmlt(i,j)/(fbot(i,j) + fside(i,j) + puny) 
         xtmp = min(xtmp, c1i)
         fbot(i,j)  = fbot(i,j)  * xtmp
         rside(i,j) = rside(i,j) * xtmp
      enddo                     ! ij

      end subroutine frzmlt_bottom_lateral

!=======================================================================
!BOP
!
! !ROUTINE: init_thermo_vars - initialize thermodynamic variables
!
! !DESCRIPTION:
!
! For current thickness category, initialize the thermodynamic
! variables that are aggregated and sent to the coupler, along
! with other fluxes passed among subroutines.
!
! !REVISION HISTORY:
!
! author William H. Lipscomb, LANL
!        C. M. Bitz, UW
!
! !INTERFACE:
!
      subroutine init_thermo_vars(strxn,    stryn,   Trefn,    Qrefn,    &   
     &                            fsensn,   flatn,   fswabsn,  flwoutn,  &
     &                            evapn,    freshn,  fsaltn,   fhnetn,   &
     &                            fswthrun, fsurf,   fcondtop, fcondbot, &
     &                            fswint,   einit,   efinal,   mvap)
!
! !USES:
! 
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind), dimension(ilo:ihi,jlo:jhi), intent(out) ::   &

         ! fields aggregated for coupler
         strxn          & ! air/ice zonal  strss                (N/m^2)
      ,  stryn          & ! air/ice merdnl strss                (N/m^2)
      ,  Trefn          & ! air tmp rfrnc level                     (K)
      ,  Qrefn          & ! air sp hum rfrnc level              (kg/kg)
      ,  fsensn         & ! surface downward sensible heat      (W m-2)
      ,  flatn          & ! surface downward latent heat        (W m-2)
      ,  fswabsn        & ! SW absorbed by ice                  (W m-2)
      ,  flwoutn        & ! upward LW at surface                (W m-2)
      ,  evapn          & ! flux of vapor, atmos to ice    (kg m-2 s-1)
      ,  freshn         & ! flux of water, ice to ocean    (kg m-2 s-1)
      ,  fsaltn         & ! flux of salt, ice to ocean     (kg m-2 s-1)
      ,  fhnetn         & ! fbot, corrected for surplus energy  (W m-2)
      ,  fswthrun       & ! SW through ice to ocean             (W m-2)
                         
         ! other fields  passed among subroutines
      ,  fsurf          & ! net flux to top surface, not including fcondtop
      ,  fcondtop       & ! downward cond flux at top surface    (W m-2)
      ,  fcondbot       & ! downward cond flux at bottom surface (W m-2)
      ,  fswint         & ! SW absorbed in ice interior, below surface (W m-2)
      ,  einit          & ! initial energy of melting            (J m-2)
      ,  efinal         & ! final energy of melting              (J m-2)
      ,  mvap            ! ice/snow mass sublimated/condensed   (kg m-2)
!
!EOP
!
      integer (kind=int_kind) ::   &
        i, j            ! horizontal indices

      do j = jlo, jhi
      do i = ilo, ihi

         strxn(i,j)      = c0i
         stryn(i,j)      = c0i
         Trefn(i,j)      = c0i
         Qrefn(i,j)      = c0i
         fsensn(i,j)     = c0i
         flatn(i,j)      = c0i
         fswabsn(i,j)    = c0i
         flwoutn(i,j)    = c0i
         evapn(i,j)      = c0i
         freshn(i,j)     = c0i
         fsaltn(i,j)     = c0i
         fhnetn(i,j)     = c0i
         fswthrun(i,j)   = c0i

         fsurf(i,j)      = c0i
         fcondtop(i,j)   = c0i
         fcondbot(i,j)   = c0i
         fswint(i,j)     = c0i
         einit(i,j)      = c0i
         efinal(i,j)     = c0i
         mvap(i,j)       = c0i

      enddo                     ! i
      enddo                     ! j

      end subroutine init_thermo_vars

!=======================================================================
!BOP
!
! !ROUTINE: init_vertical_profile - initial thickness, enthalpy, temperature
!
! !DESCRIPTION:
!
! Given the state variables (vicen, vsnon, eicen, esnon, Tsfcn),
! compute variables needed for the vertical thermodynamics
! (hin, hsn, qin, qsn, Tin, Tsn, Tsf).
!
! !REVISION HISTORY:
!
! authors William H. Lipscomb, LANL
!         C. M. Bitz, UW
!
! !INTERFACE:
!
      subroutine init_vertical_profile                   &
          (ni,    icells, indxi, indxj,                   &
           hin,  hlyr,   hsn,   hin_init, hsn_init,      &
           qin,  Tin,    qsn,   Tsn,      Tsf,   einit)
!
! !USES:
! 
!      use ice_exit
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) ::   &
         ni           &    ! thickness category index
      ,  icells          ! number of cells with aicen > puny

      integer (kind=int_kind), dimension(1:(ihi-ilo+1)*(jhi-jlo+1)),   &
         intent(in) ::   & 
         indxi, indxj    ! compressed indices for cells with aicen > puny

      real (kind=dbl_kind), dimension(ilo:ihi,jlo:jhi), intent(out) ::   &
         hin             &! ice thickness (m)
      ,  hsn             &! snow thickness (m)
      ,  hlyr            &! ice layer thickness
      ,  hin_init        &! initial value of hin
      ,  hsn_init        &! initial value of hsn
      ,  qsn             &! snow enthalpy
      ,  Tsn             &! snow temperature
      ,  Tsf             ! ice/snow surface temperature, Tsfcn

      real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi,nilyr),         &
         intent(out) ::   &
         qin             &! ice layer enthalpy (J m-3)
      ,  Tin             ! internal ice layer temperatures

      real (kind=dbl_kind), dimension(ilo:ihi,jlo:jhi), intent(inout) ::   &
         einit           ! initial energy of melting (J m-2)
!
!EOP
!
      real (kind=dbl_kind), parameter ::   &
         Tmin = -100._dbl_kind ! min allowed internal temperature (deg C)

      integer (kind=int_kind) ::   & 
         i, j           & ! horizontal indices
      ,  ij             & ! horizontal index, combines i and j loops
      ,  k                ! ice layer index

      real (kind=dbl_kind), dimension(ilo:ihi,jlo:jhi) ::   &
         Tmax            ! maximum allowed snow temperature (deg C)

      real (kind=dbl_kind) ::   &
         aa1, bb1, cc1   ! terms in quadratic formula

      logical (kind=log_kind) ::   &   ! for vector-friendly error checks
         tsno_high      & ! flag for Tsn > Tmax 
      ,  tice_high      & ! flag for Tin > Tmlt
      ,  tsno_low       & ! flag for Tsn < Tmin
      ,  tice_low        ! flag for Tin < Tmin

      !-----------------------------------------------------------------
      ! Initialize arrays
      !-----------------------------------------------------------------

      do j = jlo, jhi
      do i = ilo, ihi
         hin(i,j)      = c0i
         hsn(i,j)      = c0i
         hlyr(i,j)     = c0i
         hin_init(i,j) = c0i
         hsn_init(i,j) = c0i
         qsn(i,j)      = c0i
         Tsn(i,j)      = c0i
         Tsf(i,j)      = c0i
         Tmax(i,j)     = c0i
      enddo
      enddo

      do k = 1, nilyr
         do j = jlo, jhi
         do i = ilo, ihi
            qin(i,j,k) = c0i
            Tin(i,j,k) = c0i
         enddo
         enddo
      enddo

      !-----------------------------------------------------------------
      ! Initialize flags
      !-----------------------------------------------------------------

      tsno_high = .false.
      tice_high = .false.
      tsno_low  = .false.
      tice_low  = .false.

      !-----------------------------------------------------------------
      ! Load arrays for vertical thermo calculation.
      !-----------------------------------------------------------------
      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)

      !-----------------------------------------------------------------
      ! Surface temperature, ice and snow thickness, snow enthalpy
      !
      ! Tmax based on the idea that dT ~ dq / (rhos*cp_ice)
      !                             dq ~ q dv / v
      !                             dv ~ eps11
      ! where 'd' denotes an error due to roundoff.
      !-----------------------------------------------------------------

         Tsf(i,j)  = Tsfcn(i,j,ni)
         hin(i,j) = vicen(i,j,ni) / aicen(i,j,ni)
         hlyr(i,j) = hin(i,j) / real(nilyr,kind=dbl_kind)
         hsn(i,j) = vsnon(i,j,ni) / aicen(i,j,ni)
         hin_init(i,j) = hin(i,j)
         hsn_init(i,j) = hsn(i,j)
    
         if (hsn(i,j) > hsnomin) then
            qsn(i,j) = esnon(i,j,ni) / vsnon(i,j,ni)  ! qsn, esnon < 0
            Tmax(i,j) = -qsn(i,j)*eps11 / (rhos*cp_ice*vsnon(i,j,ni))
         else
            qsn(i,j) = -rhos * Lfresh
            Tmax(i,j) = puny
         endif

      !-----------------------------------------------------------------
      ! Compute snow temperatures from enthalpies.  
      !-----------------------------------------------------------------
         Tsn(i,j) = (Lfresh + qsn(i,j)/rhos)/cp_ice  ! qsn <= -rhos*Lfresh

      !-----------------------------------------------------------------
      ! Check for Tsn > Tmax (allowing for roundoff error)
      !  and Tsn < Tmin.
      !-----------------------------------------------------------------
         if (Tsn(i,j) > Tmax(i,j)) then
            tsno_high = .true.
         elseif (Tsn(i,j) < Tmin) then
            tsno_low  = .true.
         endif

      enddo

      !-----------------------------------------------------------------
      ! If Tsn is out of bounds, print diagnostics and exit.
      !-----------------------------------------------------------------
      if (tsno_high) then
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)

            if (Tsn(i,j) > Tmax(i,j)) then ! allowing for roundoff error
!               write(nu_diag,*) ' '
!               write(nu_diag,*) 'Starting thermo, Tsn > Tmax, cat',ni
!               write(nu_diag,*) 'Tsn=',Tsn(i,j)
!               write(nu_diag,*) 'Tmax=',Tmax(i,j)
!               write(nu_diag,*) 'istep1, my_task, i, j:',  &
!                                 istep1, my_task, i, j
!               write(nu_diag,*) 'qsn',qsn(i,j)
!               stoplabel = 'ice state at stop'
!               call print_state(stoplabel,i,j)
!               call abort_ice('vertical thermo: Tsn must be <= 0')
            endif
            
         enddo                  ! ij
      endif                     ! tsno_high

      if (tsno_low) then
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)

!            if (Tsn(i,j) < Tmin) then ! allowing for roundoff error
!               write(nu_diag,*) ' '
!               write(nu_diag,*) 'Starting thermo, Tsn < Tmin, cat',ni
!               write(nu_diag,*) 'Tsn=', Tsn(i,j)
!               write(nu_diag,*) 'Tmin=', Tmin
!               write(nu_diag,*) 'istep1, my_task, i, j:',   &
!                                 istep1, my_task, i, j
!               write(nu_diag,*) 'qsn', qsn(i,j)
!               stoplabel = 'ice state at stop'
!               call print_state(stoplabel,i,j)
!               call abort_ice('vertical thermo: Tsn must be <= 0')
!            endif 
            
         enddo                  ! ij
      endif                     ! tsno_low

      !-----------------------------------------------------------------
      ! initial energy per unit area of ice/snow, relative to 0 C
      ! incremented later for ice      
      !-----------------------------------------------------------------

      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)

         if (Tsn(i,j) > c0i) then   ! correct roundoff error
            Tsn(i,j) = c0i
            qsn(i,j) = -rhos*Lfresh
         endif

         einit(i,j) = hsn(i,j)*qsn(i,j)

      enddo                     ! ij

      do k = 1, nilyr
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)

      !-----------------------------------------------------------------
      ! Compute ice enthalpy
      !-----------------------------------------------------------------

            qin(i,j,k) = eicen(i,j,ilyr1(ni)+k-1)   & ! qin < 0
                        * real(nilyr,kind=dbl_kind)/vicen(i,j,ni)

      !-----------------------------------------------------------------
      ! Compute ice temperatures from enthalpies using quadratic formula
      !-----------------------------------------------------------------

            if (l_brine) then
               aa1 = cp_ice
               bb1 = (cp_ocn-cp_ice)*Tmlt(k) - qin(i,j,k)/rhoi - Lfresh 
               cc1 = Lfresh * Tmlt(k)
               Tin(i,j,k) =  (-bb1 - sqrt(bb1*bb1 - c4i*aa1*cc1)) /   &
                             (c2i*aa1)
               Tmax(i,j) = Tmlt(k)

            else                ! fresh ice
               Tin(i,j,k) = (Lfresh + qin(i,j,k)/rhoi) / cp_ice
               Tmax(i,j) = -qin(i,j,k)*eps11/(rhos*cp_ice*vicen(i,j,ni))
                         ! as above for snow
            endif

           IF (Tin(i,j,k) < Tmin)THEN
              !(Tin=Tmlt(k))
               eicen(i,j,ilyr1(ni)+k-1) =              &
                    (rhoi *cp_ocn*Tmlt(k)) &
                    * vicen(i,j,ni)/real(nilyr,kind=dbl_kind)
              qin(i,j,k) = eicen(i,j,ilyr1(ni)+k-1)   & ! qin < 0
                        * real(nilyr,kind=dbl_kind)/vicen(i,j,ni)
              aa1 = cp_ice
              bb1 = (cp_ocn-cp_ice)*Tmlt(k) - qin(i,j,k)/rhoi - Lfresh
              cc1 = Lfresh * Tmlt(k)
              Tin(i,j,k) =  (-bb1 - sqrt(bb1*bb1 - c4i*aa1*cc1)) /   &
                             (c2i*aa1)
           end if


            
      !-----------------------------------------------------------------
      ! Check for Tin > Tmax and Tin < Tmin
      !-----------------------------------------------------------------

            if (Tin(i,j,k) > Tmax(i,j)) then
               tice_high = .true.
            elseif (Tin(i,j,k) < Tmin) then
               tice_low  = .true.
            endif

         enddo                  ! ij

      !-----------------------------------------------------------------
      ! If Tin is out of bounds, print diagnostics and exit.
      !-----------------------------------------------------------------
     
         if (tice_high) then
            do ij = 1, icells
               i = indxi(ij)
               j = indxj(ij)

               if (Tin(i,j,k) > Tmax(i,j)) then
!                  write(nu_diag,*) ' '
!                  write(nu_diag,*) 'Starting thermo, T > Tmax, cat',&
!                       ni,', layer',k                                 
!                  write(nu_diag,*) 'Tin=',Tin(i,j,k),', Tmax=',Tmax(i,j)
!                  write(nu_diag,*) 'istep1, my_task, i, j:',        &
!                       istep1, my_task, i, j
!                  write(nu_diag,*) 'qin',qin(i,j,k)
!                  stoplabel = 'ice state at stop'
!                  call print_state(stoplabel,i,j)
!                  call abort_ice('vertical thermo: Tin must be <= Tmax')
               endif
               
            enddo               ! ij
         endif                  ! tice_high   

         if (tice_low) then
            do ij = 1, icells
               i = indxi(ij)
               j = indxj(ij)
               
!               if (Tin(i,j,k) < Tmin) then
!                  write(nu_diag,*) ' '
!                  write(nu_diag,*) 'Starting thermo T < Tmin, cat',  &
!                       ni,', layer',k                                  
!                  write(nu_diag,*) 'Tin =', Tin(i,j,k)                
!                  write(nu_diag,*) 'Tmin =', Tmin                     
!!                  write(nu_diag,*) 'istep1, my_task, i, j:',         &
!                       istep1, my_task, i, j
!                  stoplabel = 'ice state at stop'
!                  call print_state(stoplabel,i,j)
!                  call abort_ice('vertical_thermo: Tin must be > Tmin') 
!               endif
            enddo               ! ij
         endif                  ! tice_low

      !-----------------------------------------------------------------
      ! initial energy per unit area of ice/snow, relative to 0 C
      !-----------------------------------------------------------------

         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)

            if (Tin(i,j,k) > c0i) then ! correct roundoff error
               Tin(i,j,k) = c0i
               qin(i,j,k) = -rhoi*Lfresh
            endif
            
            einit(i,j) = einit(i,j) + hlyr(i,j)*qin(i,j,k) 

         enddo                  ! ij
      enddo                     ! k
      
      end subroutine init_vertical_profile

!=======================================================================
!BOP
!
! !ROUTINE: temperature_changes  - new vertical temperature profile
!
! !DESCRIPTION:
!
! Compute new surface temperature and internal ice and snow 
! temperatures.  Include effects of salinity on sea ice heat 
! capacity in a way that conserves energy (Bitz and Lipscomb, 1999).
!
! New temperatures are computed iteratively by solving a tridiagonal
! system of equations; heat capacity is updated with each iteration.
! Finite differencing is backward implicit.
!
! !REVISION HISTORY:
!
! authors William H. Lipscomb, LANL
!         C. M. Bitz, UW
!
! !INTERFACE:
!
      subroutine temperature_changes                                    &
             (ni,       icells,   indxi,  indxj,                         &
              hlyr,    hsn,      qin,    Tin,    qsn,      Tsn,         &
              Tsf,     Tbot,     fbot,   fsensn, flatn,    fswabsn,     &
              flwoutn, fswthrun, fhnetn, fsurf,  fcondtop, fcondbot,    &
              fswint,  einit)
!
! !USES:
!
!      use ice_exit
! 
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) ::   &
         ni             &  ! thickness category index
      ,  icells          ! number of cells with aicen > puny

      integer (kind=int_kind), dimension(1:(ihi-ilo+1)*(jhi-jlo+1)),   &
        intent(in) ::   & 
        indxi, indxj    ! compressed indices for cells with aicen > puny

      real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi), intent(inout)::   &
         hlyr          &  ! ice layer thickness
      ,  hsn           &  ! snow thickness (m)
      ,  Tbot          &  ! ice bottom surface temperature (deg C)
      ,  fbot          &  ! ice-ocean heat flux at bottom surface (W/m^2)
      ,  qsn           &  ! snow enthalpy
      ,  Tsn           &  ! internal snow temperature
      ,  Tsf           &  ! ice/snow surface temperature, Tsfcn
      ,  fsensn        &  ! surface downward sensible heat (W m-2)
      ,  flatn         &  ! surface downward latent heat (W m-2)
      ,  fswabsn       &  ! shortwave absorbed by ice (W m-2)
      ,  flwoutn       &  ! upward LW at surface (W m-2)
      ,  fswthrun      &  ! SW through ice to ocean (W m-2)
      ,  fhnetn        &  ! fbot, corrected for any surplus energy
      ,  fsurf         &  ! net flux to top surface, not including fcondtop
      ,  fcondtop      &  ! downward cond flux at top surface (W m-2)
      ,  fcondbot      &  ! downward cond flux at bottom surface (W m-2)
      ,  fswint          ! SW absorbed in ice interior, below surface (W m-2)

      real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi,nilyr), &
         intent(inout) ::   &
         qin           &  ! ice layer enthalpy (J m-3)
      ,  Tin             ! internal ice layer temperatures

      real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi), intent(in)::   &
         einit           ! initial energy of melting (J m-2)
!
!EOP
!
      integer (kind=int_kind), parameter ::   &
# if defined (ICE_FRESHWATER)
! afm 20151113 & EJA 20160921
! Solution converges in the first 10 iternation. 
! Output did not change signigicantly with nitermax=50-->10,
! but saved some computational time.
         nitermax = 10   ! max number of iterations in temperature solver
# else
         nitermax = 50   ! max number of iterations in temperature solver
# endif

      real (kind=dbl_kind), parameter ::   &
         alph =  c3i     & ! constant used to get 2nd order accurate fluxes
      ,  bet  = -p333   & ! constant used to get 2nd order accurate fluxes
      ,  Tsf_errmax = 5.e-4 ! max allowed error in Tsf
                            ! recommend Tsf_errmax < 0.01 K

      integer (kind=int_kind) ::   & 
         i, j           & ! horizontal indices
      ,  ij             & ! horizontal index, combines i and j loops
      ,  k              & ! ice layer index
      ,  nmat           & ! matrix dimension
      ,  niter          & ! iteration counter in temperature solver
      ,  isolve          ! number of cells with temps not converged

      integer (kind=int_kind), dimension(1:(ihi-ilo+1)*(jhi-jlo+1)) ::   &
         indxii, indxjj  ! compressed indices for cells not converged

      real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi) ::   &
         Tsn_init      &  ! Tsn at beginning of time step
      ,  Tsn_start     &  ! Tsn at start of iteration
      ,  Tsf_start     &  ! Tsf at start of iteration
      ,  dTsf          &  ! Tsf - Tsf_start
      ,  dTsf_prev     &  ! dTsf from previous iteration
      ,  dfsurf_dT     &  ! derivative of fsurf wrt Tsf
      ,  dfsens_dT     &  ! deriv of fsens wrt Tsf (W m-2 deg-1)
      ,  dflat_dT      &  ! deriv of flat wrt Tsf (W m-2 deg-1)
      ,  dflwout_dT    &  ! deriv of flwout wrt Tsf (W m-2 deg-1)
      ,  fswsfc          ! SW absorbed at ice/snow surface (W m-2)

      real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi) ::   &
         khs            & ! ksno / hsn
      ,  etas           & ! dt / (rhos * cp_ice * hsn)
      ,  dt_rhoi_hlyr   & ! dt/(rhoi*hlyr)
      ,  avg_Tin        & ! = 1. if Tin averaged w/Tin_start, else = 0.
      ,  enew            ! new energy of melting after temp change (J m-2)

      real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi,nilyr) ::   &
         Tin_init      &  ! Tin at beginning of time step
      ,  Tin_start     &  ! Tin at start of iteration
      ,  Iabs          &  ! SW absorbed in particular layer (W m-2) 
      ,  etai            ! dt / (rho * cp * h) for ice layers

      real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi,nilyr+1) ::   &
         khi             ! ki / hlyr

      real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi,nilyr+2) ::   &
         diag          &  ! diagonal matrix elements
      ,  sbdiag        &  ! sub-diagonal matrix elements
      ,  spdiag        &  ! super-diagonal matrix elements
      ,  rhs           &  ! rhs of tri-diagonal matrix equation
      ,  Tmat            ! matrix output temperatures

      real (kind=dbl_kind) ::   &
         ci            & ! specific heat of sea ice (J kg-1 deg-1)
      ,  spdiag2       & ! term to right of superdiag on top row
      ,  avg_Tsf       & ! = 1. if Tsf averaged w/Tsf_start, else = 0.
      ,  avg_Tsn       & ! = 1. if Tsn averaged w/Tsn_start, else = 0.
      ,  ferr          & ! energy conservation error (W m-2)
      ,  w1,w2,w3,w4,w5,w6,w7  ! temporary variables

      logical (kind=log_kind), dimension (ilo:ihi,jlo:jhi) ::   &
         converged      ! = true when local solution has converged

      logical (kind=log_kind) ::   &
         all_converged  ! = true when all cells have converged

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

         all_converged   = .false.

      do ij = 1, icells

         i = indxi(ij)
         j = indxj(ij)

         converged(i,j)  = .false.
         khs      (i,j)  = c0i
         dTsf     (i,j)  = c0i
         dTsf_prev(i,j)  = c0i
         Tsf_start(i,j)  = c0i
         enew     (i,j)  = c0i

         dfsurf_dT (i,j) = c0i
         dfsens_dT (i,j) = c0i
         dflat_dT  (i,j) = c0i
         dflwout_dT(i,j) = c0i

         fhnetn(i,j) = fbot(i,j)  ! ocean energy used by the ice, <= 0

         fswsfc(i,j)     = c0i
!         dt_rhoi_hlyr(i,j) = dt / (rhoi*hlyr(i,j)) ! used by matrix solver
         dt_rhoi_hlyr(i,j) = dtice / (rhoi*hlyr(i,j)) ! used by matrix solver

         if (hsn(i,j) > hsnomin) then
!            etas(i,j) = dt / (rhos*cp_ice*hsn(i,j))  ! used by matrix solver
            etas(i,j) = dtice / (rhos*cp_ice*hsn(i,j))  ! used by matrix solver
         else
            etas(i,j) = c0i
         endif

         Tsn_init (i,j) = Tsn(i,j)   ! beginning of time step
         Tsn_start(i,j) = Tsn(i,j)   ! beginning of iteration

      enddo                     ! ij

      do k = 1, nilyr
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)

            Iabs     (i,j,k) = c0i
            Tin_init (i,j,k) = Tin(i,j,k)  ! beginning of time step
            Tin_start(i,j,k) = Tin(i,j,k)  ! beginning of iteration
         enddo                  ! ij
      enddo                     ! k

      !-----------------------------------------------------------------
      ! Compute thermal conductivity at interfaces (held fixed during 
      !  the subsequent iteration).
      !-----------------------------------------------------------------
      call conductivity                                   &
          (icells, indxi, indxj,                          &
           hlyr,   hsn,   Tin,  Tbot, khi, khs)            
                                                           
      !-----------------------------------------------------------------
      ! Compute solar radiation absorbed in ice and penetrating to ocean
      !-----------------------------------------------------------------
      call absorbed_solar                                 &
          (ni,    icells, indxi,  indxj,                   &
           hlyr, hsn,    fswsfc, fswint, fswthrun, Iabs)

      !-----------------------------------------------------------------
      ! Solve for new temperatures.  
      ! Iterate until temperatures converge with minimal energy error.
      !-----------------------------------------------------------------

      do niter = 1, nitermax      

         if (all_converged) then  ! thermo calculation is done
            exit
         else                     ! identify cells not yet converged
            isolve = 0
            do ij = 1, icells
               i = indxi(ij)
               j = indxj(ij)
               if (.not.converged(i,j)) then
                  isolve = isolve + 1
                  indxii(isolve) = i
                  indxjj(isolve) = j
               endif
            enddo               ! ij
         endif

      !-----------------------------------------------------------------
      ! Update radiative and turbulent fluxes and their derivatives
      ! with respect to Tsf.  
      !-----------------------------------------------------------------

         call surface_fluxes                                &
              (isolve,     indxii,    indxjj,               &
               Tsf,        fswsfc,                          &
               flwoutn,    fsensn,    flatn,    fsurf,      &
               dflwout_dT, dfsens_dT, dflat_dT, dfsurf_dT) 

         do ij = 1, isolve
            i = indxii(ij)   ! NOTE: not indxi and indxj
            j = indxjj(ij)

      !-----------------------------------------------------------------
      ! Compute conductive flux at top surface, fcondtop.
      ! If fsurf < fcondtop and Tsf = 0, then reset Tsf to slightly less 
      !  than zero (but not less than -puny).
      !-----------------------------------------------------------------

            if (hsn(i,j) > hsnomin) then
               fcondtop(i,j) = c2i * khs(i,j) * (Tsf(i,j) - Tsn(i,j))
            else
               fcondtop(i,j) = khi(i,j,1)                  &
                    * (alph*(Tsf(i,j) - Tin(i,j,1))        &
                      + bet*(Tsf(i,j) - Tin(i,j,2)))        
            endif                                           
            if (fsurf(i,j) < fcondtop(i,j))                &
                 Tsf(i,j) = min (Tsf(i,j), -puny)
             
      !-----------------------------------------------------------------
      ! Save surface temperature at start of iteration
      !-----------------------------------------------------------------
            Tsf_start(i,j) = Tsf(i,j)

         enddo                  ! ij

      !-----------------------------------------------------------------
      ! Compute specific heat of each ice layer
      !-----------------------------------------------------------------

         do k = 1, nilyr
            do ij = 1, isolve
               i = indxii(ij)
               j = indxjj(ij)

               if (l_brine) then
                  ci = cp_ice - Lfresh*Tmlt(k) /             &
                                (Tin_init(i,j,k)*Tin(i,j,k))
               else
                  ci = cp_ice
               endif
               etai(i,j,k) = dt_rhoi_hlyr(i,j) / ci

            enddo               ! ij
         enddo                  ! k

         do ij = 1, isolve
            i = indxii(ij)
            j = indxjj(ij)

      !-----------------------------------------------------------------
      ! Compute matrix elements
      !
      ! Four possible cases to solve:
      !   (1) Cold surface (Tsf < 0), snow present
      !   (2) Melting surface (Tsf = 0), snow present
      !   (3) Cold surface (Tsf < 0), no snow
      !   (4) Melting surface (Tsf = 0), no snow
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! first row
      ! Tsf for case 1; dummy equation for cases 2, 3 and 4
      !-----------------------------------------------------------------

            sbdiag(i,j,1) = c0i
            spdiag2 = c0i

            if (hsn(i,j) > hsnomin) then
               if (Tsf(i,j) <= -puny) then
                  spdiag(i,j,1) = c2i*khs(i,j)  
                  diag(i,j,1)   = dfsurf_dT(i,j) - c2i*khs(i,j)
                  rhs(i,j,1)    = dfsurf_dT(i,j)*Tsf(i,j) - fsurf(i,j)
               else
                  spdiag(i,j,1) = c0i
                  diag(i,j,1)   = c1i
                  rhs(i,j,1)    = c0i
               endif
            else                ! hsn <= hsnomin
               if (Tsf(i,j) <= -puny) then
                  spdiag(i,j,1) = c0i
                  diag(i,j,1)   = c1i
                  rhs(i,j,1)    = c0i
               else
                  spdiag(i,j,1) = c0i
                  diag(i,j,1)   = c1i
                  rhs(i,j,1)    = c0i
               endif
            endif

      !-----------------------------------------------------------------
      ! second row
      ! Tsn for cases 1 and 2; Tsf for case 3; dummy for case 4
      !-----------------------------------------------------------------

            if (hsn(i,j) > hsnomin) then
               if (Tsf(i,j) <= -puny) then
                  sbdiag(i,j,2) = -etas(i,j) * c2i * khs(i,j)
                  spdiag(i,j,2) = -etas(i,j) * khi(i,j,1)
                  spdiag2       = c0i
                  diag(i,j,2)   = c1i                                      &
                                + etas(i,j) * (c2i*khs(i,j) + khi(i,j,1))   
                  rhs(i,j,2)    = Tsn_init(i,j)                            
               else                                                        
                  sbdiag(i,j,2) = c0i                                       
                  spdiag(i,j,2) = -etas(i,j) * khi(i,j,1)                  
                  spdiag2       = c0i                                       
                  diag(i,j,2)   = c1i                                      &
                                + etas(i,j) * (c2i*khs(i,j) + khi(i,j,1))   
                  rhs(i,j,2)    = Tsn_init(i,j)                           &
                                + etas(i,j)*c2i*khs(i,j)*Tsf(i,j)           
               endif                                                       
            else                ! hsn <= hsnomin                           
               if (Tsf(i,j) <= -puny) then                                 
                  sbdiag(i,j,2) = c0i                                       
                  spdiag(i,j,2) = alph * khi(i,j,1)                        
                  spdiag2       =  bet * khi(i,j,1)                        
                  diag(i,j,2)   = dfsurf_dT(i,j) - alph*khi(i,j,1)        &
                                                 -  bet*khi(i,j,1)
                  rhs(i,j,2)    = dfsurf_dT(i,j)*Tsf(i,j) - fsurf(i,j)
               else
                  sbdiag(i,j,2) = c0i
                  spdiag(i,j,2) = c0i
                  spdiag2       = c0i
                  diag(i,j,2)   = c1i
                  rhs(i,j,2)    = c0i
               endif
            endif

      !-----------------------------------------------------------------
      ! third row, Tin(i,j,1)
      !
      ! For each internal ice layer compute the specific heat ci, 
      ! a function of the starting temperature and of the latest
      ! guess for the final temperature.
      !-----------------------------------------------------------------

            if (hsn(i,j) > hsnomin) then
               if (Tsf(i,j) <= -puny) then
                  sbdiag(i,j,3) = -etai(i,j,1) * khi(i,j,1)
                  spdiag(i,j,3) = -etai(i,j,1) * khi(i,j,2)
                  diag(i,j,3)   = c1i                                    &
                                + etai(i,j,1)*(khi(i,j,1) + khi(i,j,2))  
                  rhs(i,j,3)    = Tin_init(i,j,1)                       &
                                + etai(i,j,1)*Iabs(i,j,1)                
               else                                                      
                  sbdiag(i,j,3) = -etai(i,j,1) * khi(i,j,1)              
                  spdiag(i,j,3) = -etai(i,j,1) * khi(i,j,2)              
                  diag(i,j,3)   = c1i                                    &
                                + etai(i,j,1)*(khi(i,j,1) + khi(i,j,2))  
                  rhs(i,j,3)    = Tin_init(i,j,1)                       &
                                + etai(i,j,1)*Iabs(i,j,1)                
               endif                                                     
            else                ! hsn <= hsnomin                         
               if (Tsf(i,j) <= -puny) then                               
                  sbdiag(i,j,3) = -(alph+bet) * etai(i,j,1) * khi(i,j,1) 
                  spdiag(i,j,3) = -etai(i,j,1)                          &
                                * (khi(i,j,2) - bet*khi(i,j,1))          
                  diag(i,j,3)   = c1i + etai(i,j,1)                      &
                                     * (khi(i,j,2) + alph*khi(i,j,1))    
                  rhs(i,j,3)    = Tin_init(i,j,1)                       &
                                + etai(i,j,1)*Iabs(i,j,1)                
               else                                                      
                  sbdiag(i,j,3) = c0i                                     
                  spdiag(i,j,3) = -etai(i,j,1)                          &
                                * (khi(i,j,2) - bet*khi(i,j,1))          
                  diag(i,j,3)   = c1i + etai(i,j,1)                      &
                                * (khi(i,j,2) + alph*khi(i,j,1))         
                  rhs(i,j,3)    = Tin_init(i,j,1)                       &
                                + etai(i,j,1)*Iabs(i,j,1)
!!!     &                 + (alph+bet)*etai(i,j,k)*khi(i,j,1)*Tsf(i,j) 
                       ! commented line needed if Tsf were in Kelvin
               endif
            endif

            if (hsn(i,j) <= hsnomin .and. Tsf(i,j) <= -puny) then ! case 3
               ! remove spdiag2 in row 2 by adding rows 2 and 3 
               diag(i,j,2)   = spdiag(i,j,3) * diag(i,j,2)              &
                                   - spdiag2 * sbdiag(i,j,3)             
               spdiag(i,j,2) = spdiag(i,j,3) * spdiag(i,j,2)            &
                                   - spdiag2 * diag(i,j,3)               
               rhs(i,j,2)    = spdiag(i,j,3) * rhs(i,j,2)               &
                                   - spdiag2 * rhs(i,j,3)                
            endif                                                        
                                                                         
      !----------------------------------------------------------------- 
      ! Bottom row, Tin(i,j,nilyr)                                       
      !----------------------------------------------------------------- 
                                                                         
            k = nilyr                                                    
            sbdiag(i,j,k+2) = -etai(i,j,k)                              &
                            * (khi(i,j,k) - bet*khi(i,j,k+1))            
            spdiag(i,j,k+2) =  c0i                                        
            diag(i,j,k+2)   = c1i + etai(i,j,k)                          &
                            * (khi(i,j,k) + alph*khi(i,j,k+1))           
            rhs(i,j,k+2)    = Tin_init(i,j,k) + etai(i,j,k)*Iabs(i,j,k) &
                            + (alph+bet)*etai(i,j,k)                    &
                             *khi(i,j,k+1)*Tbot(i,j)

         enddo                  ! ij

      !-----------------------------------------------------------------
      ! Ice interior
      !-----------------------------------------------------------------
         do k = 2, nilyr-1
            do ij = 1, isolve
               i = indxii(ij)
               j = indxjj(ij)

               sbdiag(i,j,k+2) =  -etai(i,j,k) * khi(i,j,k)
               spdiag(i,j,k+2) =  -etai(i,j,k) * khi(i,j,k+1)
               diag(i,j,k+2)   =  c1i - sbdiag(i,j,k+2) - spdiag(i,j,k+2)
               rhs(i,j,k+2)    =  Tin_init(i,j,k)               &
                                + etai(i,j,k) * Iabs(i,j,k)

            enddo               ! ij
         enddo                  ! k

      !-----------------------------------------------------------------
      ! Solve tridiagonal matrix to obtain the new temperatures.
      !-----------------------------------------------------------------

         nmat = nilyr + 2

         call tridiag_solver                &
             (isolve, indxii, indxjj,       &
              nmat,   sbdiag, diag,  spdiag, rhs, Tmat)
# if defined (ICE_FRESHWATER)
! afm 20151114 & EJA 20160921 - isnan check
        do k = 1, nilyr+2
         do ij = 1, isolve
          i= indxii(ij); j = indxjj(ij)
          if (isnan(Tmat(i,j,k)) /= 0 ) then
           print*, 'isnan tisfc :', i,j,k,tsf(i,j), tmat(i,j,k)
!           tmat(i,j,k) = potT(i,j) - Tffresh
           tmat(i,j,k) = c0i
           end if
          end do
         end do
# endif

      !-----------------------------------------------------------------
      ! Reload temperatures from matrix solution.
      !-----------------------------------------------------------------

         do ij = 1, isolve
            i = indxii(ij)
            j = indxjj(ij)

            if (hsn(i,j) > hsnomin) then
               if (Tsf(i,j) <= -puny) then
                  Tsf(i,j) = Tmat(i,j,1)
                  Tsn(i,j) = Tmat(i,j,2)
               else
                  Tsf(i,j) = c0i
                  Tsn(i,j) = Tmat(i,j,2)
               endif
            else                ! hsn <= hsnomin
               if (Tsf(i,j) <= -puny) then
                  Tsf(i,j) = Tmat(i,j,2)
               else
                  Tsf(i,j) = c0i
               endif
            endif

         enddo

         do k = 1, nilyr
            do ij = 1, isolve
               i = indxii(ij)
               j = indxjj(ij)

               Tin(i,j,k) = Tmat(i,j,k+2)
            enddo
         enddo

      !-----------------------------------------------------------------
      ! Determine whether the computation has converged to an acceptable
      ! solution.  Five conditions must be satisfied:
      !
      !    (1) Tsf <= 0 C.
      !    (2) Tsf is not oscillating; i.e., if both dTsf(niter) and
      !        dTsf(niter-1) have magnitudes greater than puny, then
      !        dTsf(niter)/dTsf(niter-1) cannot be a negative number
      !        with magnitude greater than 0.5.  
      !    (3) abs(dTsf) < Tsf_errmax
      !    (4) If Tsf = 0 C, then the downward turbulent/radiative 
      !        flux, fsurf, must be greater than or equal to the downward
      !        conductive flux, fcondtop.
      !    (5) The net energy added to the ice per unit time must equal 
      !        the net change in internal ice energy per unit time,
      !        within the prescribed error ferrmax.
      !
      ! For briny ice (the standard case), Tsn and Tin are limited
      !  to prevent them from exceeding their melting temperatures.
      !  (Note that the specific heat formula for briny ice assumes
      !  that T < Tmlt.)  
      ! For fresh ice there is no limiting, since there are cases
      !  when the only convergent solution has Tsn > 0 and/or Tin > 0.
      !  Above-zero temperatures are then reset to zero (with melting 
      !  to conserve energy) in the thickness_changes subroutine.
      !-----------------------------------------------------------------

         ! initialize global convergence flag
         all_converged = .true.

         do ij = 1, isolve
            i = indxii(ij)
            j = indxjj(ij)

            if (l_brine) Tsn(i,j) = min(Tsn(i,j), c0i)

      !-----------------------------------------------------------------
      ! Initialize convergence flag (true until proven false), dTsf,
      !  and temperature-averaging coefficients.
      ! Average only if test 1 or 2 fails.
      !-----------------------------------------------------------------

            converged(i,j) = .true.
            dTsf(i,j) = Tsf(i,j) - Tsf_start(i,j)
            avg_Tsf      = c0i
            avg_Tsn      = c0i
            avg_Tin(i,j) = c0i

      !-----------------------------------------------------------------
      ! Condition 1: check for Tsf > 0
      ! If Tsf > 0, set Tsf = 0, then average Tsn and Tin to force 
      ! internal temps below their melting temps.
      !-----------------------------------------------------------------

            if (Tsf(i,j) > puny) then
               Tsf(i,j) = c0i
               dTsf(i,j) = -Tsf_start(i,j)
               if (l_brine) then ! average with starting temp
                  avg_Tsn = c1i  
                  avg_Tin(i,j) = c1i
               endif
               converged(i,j) = .false.
               all_converged = .false.

      !-----------------------------------------------------------------
      ! Condition 2: check for oscillating Tsf
      ! If oscillating, average all temps to increase rate of convergence.
      !-----------------------------------------------------------------

            elseif (niter > 1                  &! condition (2)
              .and. Tsf_start(i,j) <= -puny    &
              .and. abs(dTsf(i,j)) > puny      &
              .and. abs(dTsf_prev(i,j)) > puny &
              .and. -dTsf(i,j)/(dTsf_prev(i,j)+puny*puny) > p5) then       

               if (l_brine) then ! average with starting temp
                  avg_Tsf = c1i 
                  avg_Tsn = c1i   
                  avg_Tin(i,j) = c1i
               endif
               dTsf(i,j) = p5 * dTsf(i,j)
               converged(i,j) = .false.
               all_converged = .false.
            endif

      !-----------------------------------------------------------------
      ! If condition 1 or 2 failed, average new surface/snow 
      !  temperatures with their starting values.
      ! (No change if both tests passed)
      !-----------------------------------------------------------------
            Tsf(i,j)  = Tsf(i,j)                                      &
                      + avg_Tsf * p5 * (Tsf_start(i,j) - Tsf(i,j))     
            Tsn(i,j) = Tsn(i,j)                                       &
                      + avg_Tsn * p5 * (Tsn_start(i,j) - Tsn(i,j))
           
      !-----------------------------------------------------------------
      ! Compute qsn and increment new energy.
      !-----------------------------------------------------------------
            qsn(i,j) = -rhos * (Lfresh - cp_ice*Tsn(i,j))
            enew(i,j) = hsn(i,j) * qsn(i,j)

            Tsn_start(i,j) = Tsn(i,j) ! for next iteration

         enddo  ! ij

         do k = 1, nilyr
            do ij = 1, isolve
               i = indxii(ij)
               j = indxjj(ij)

               if (l_brine) Tin(i,j,k) = min (Tin(i,j,k), Tmlt(k))

      !-----------------------------------------------------------------
      ! If condition 1 or 2 failed, average new ice layer 
      !  temperatures with their starting values.
      !-----------------------------------------------------------------
               Tin(i,j,k) = Tin(i,j,k)                                  &
                    + avg_Tin(i,j)*p5*(Tin_start(i,j,k)-Tin(i,j,k))      
                                                                         
      !----------------------------------------------------------------- 
      ! Compute qin and increment new energy.                            
      !----------------------------------------------------------------- 
# if defined (ICE_FRESHWATER)
! 20151112 afm avoid zero division
               if ( Tin(i,j,k) /= c0i ) then 
               qin(i,j,k) = -rhoi * (cp_ice*(Tmlt(k)-Tin(i,j,k))        &
                                   + Lfresh*(c1i-Tmlt(k)/Tin(i,j,k))     &
                                   - cp_ocn*Tmlt(k)) 
               else
               qin(i,j,k) = -rhoi * Lfresh     
               end if
!------
# else
               qin(i,j,k) = -rhoi * (cp_ice*(Tmlt(k)-Tin(i,j,k))        &
                                   + Lfresh*(c1i-Tmlt(k)/Tin(i,j,k))     &
                                   - cp_ocn*Tmlt(k)) 
# endif

               enew(i,j) = enew(i,j) + hlyr(i,j)*qin(i,j,k)

               Tin_start(i,j,k) = Tin(i,j,k) ! for next iteration

            enddo               ! ij
         enddo                  ! k               

         do ij = 1, isolve
            i = indxii(ij)
            j = indxjj(ij)

      !-----------------------------------------------------------------
      ! Condition 3: check for large change in Tsf
      !-----------------------------------------------------------------

            if (abs(dTsf(i,j)) > Tsf_errmax) then
               converged(i,j) = .false.
               all_converged = .false.
            endif

      !-----------------------------------------------------------------
      ! Condition 4: check for fsurf < fcondtop with Tsf > 0
      !-----------------------------------------------------------------

            fsurf(i,j) = fsurf(i,j) + dTsf(i,j)*dfsurf_dT(i,j)
            if (hsn(i,j) > hsnomin) then
               fcondtop(i,j) = c2i * khs(i,j) * (Tsf(i,j)-Tsn(i,j))
            else
               fcondtop(i,j) = khi(i,j,1)            &
                    * (alph*(Tsf(i,j)-Tin(i,j,1))    &
                      + bet*(Tsf(i,j)-Tin(i,j,2)))
            endif

            if (Tsf(i,j) > -puny .and. fsurf(i,j) < fcondtop(i,j)) then
               converged(i,j) = .false.                 
               all_converged = .false.
            endif

      !-----------------------------------------------------------------
      ! Condition 5: check for energy conservation error
      ! Change in internal ice energy should equal net energy input.
      !-----------------------------------------------------------------

            fcondbot(i,j) = khi(i,j,nilyr+1) *                   &
                       (alph*(Tin(i,j,nilyr)   - Tbot(i,j))      &
                       + bet*(Tin(i,j,nilyr-1) - Tbot(i,j)))      
                                                                  
!            ferr = abs( (enew(i,j)-einit(i,j))/dt                &
            ferr = abs( (enew(i,j)-einit(i,j))/dtice              &
                      - (fcondtop(i,j) - fcondbot(i,j) + fswint(i,j)) )

            ! factor of 0.9 allows for roundoff errors later
            if (ferr > 0.9_dbl_kind*ferrmax) then           ! condition (5)
               converged(i,j) = .false.
               all_converged = .false.
            endif
         
            dTsf_prev(i,j) = dTsf(i,j)

         enddo                  ! ij

      enddo                     ! temperature iteration

      if (.not.all_converged .and.DBG_SET(DBG_SBRIO) ) then
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)

      !-----------------------------------------------------------------
      ! Check for convergence failures.
      !-----------------------------------------------------------------
            if (.not.converged(i,j)) then
!               write(nu_diag,*) 'Thermo iteration does not converge,',     &
!                                'istep1, my_task, i, j, n:',               &
!                                 istep1, my_task, i, j, ni                   
!               write(nu_diag,*) 'Ice thickness:', hlyr(i,j)*nilyr           
!               write(nu_diag,*) 'Snow thickness:', hsn(i,j)                 
!               write(nu_diag,*) 'dTsf, Tsf_errmax:',dTsf(i,j),Tsf_errmax    
!               write(nu_diag,*) 'Tsf:', Tsf(i,j)                            
!               write(nu_diag,*) 'fsurf:', fsurf(i,j)                        
!               write(nu_diag,*) 'fcondtop, fcondbot, fswint',              &
!                               fcondtop(i,j), fcondbot(i,j), fswint(i,j)    
!               write(nu_diag,*) 'Flux conservation error =', ferr           
!               write(nu_diag,*) 'Initial snow and ice temperatures:'        
!               write(nu_diag,*) Tsn_init(i,j),                             &
!                               (Tin_init(i,j,k),k=1,nilyr)
!               write(nu_diag,*) 'Final temperatures:'
!               write(nu_diag,*) Tsn(i,j), (Tin(i,j,k),k=1,nilyr)
!               stoplabel = 'ice state at stop'
!               call print_state(stoplabel,i,j)
!               call abort_ice('vertical thermo: convergence error')
            endif
         enddo                  ! ij
      endif                     ! all_converged

      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)

         ! update fluxes that depend on Tsf
         flwoutn(i,j) = flwoutn(i,j) + dTsf(i,j) * dflwout_dT(i,j)
         fsensn(i,j)  = fsensn(i,j)  + dTsf(i,j) * dfsens_dT(i,j)
         flatn(i,j)   = flatn(i,j)   + dTsf(i,j) * dflat_dT(i,j)

         ! absorbed shortwave flux for coupler
         fswabsn(i,j) = fswsfc(i,j) + fswint(i,j) + fswthrun(i,j)

      enddo                     ! ij

      end subroutine temperature_changes
      
!=======================================================================
!BOP
!
! !ROUTINE: conductivity - compute ice thermal conductivity
!
! !DESCRIPTION:
!
! Compute thermal conductivity at interfaces (held fixed during 
!  the subsequent iteration).
!
! NOTE: Ice conductivity must be >= kimin
!
! !REVISION HISTORY:
!
! authors William H. Lipscomb, LANL
!         C. M. Bitz, UW
!
! !INTERFACE:
!
      subroutine conductivity          & 
         (icells, indxi, indxj,        &
          hlyr,   hsn,   Tin,  Tbot, khi, khs)
!
! !USES:
! 
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) ::   & 
         icells          ! number of cells with aicen > puny

      integer (kind=int_kind), dimension(1:(ihi-ilo+1)*(jhi-jlo+1)),  &
         intent(in) ::   & 
         indxi, indxj    ! compressed indices for cells with aicen > puny

      real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi), intent(in) :: &
         hlyr          &  ! ice layer thickness
      ,  hsn           &  ! snow layer thickness
      ,  Tbot            ! ice bottom surface temperature (deg C)

      real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi,nilyr),       &
         intent(in) ::   &
         Tin             ! internal ice layer temperatures

      real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi,nilyr+1),   &
         intent(out) ::   &
         khi             ! ki / hlyr

      real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi), intent(out) ::   &
         khs             ! ksno / hsn
!
!EOP
!
      real (kind=dbl_kind), parameter ::   &
         betak   = 0.13_dbl_kind &! constant in formula for k (W m-1 ppt-1)
      ,  kimin   = 0.10_dbl_kind ! min conductivity of saline ice (W m-1 deg-1)

      integer (kind=int_kind) ::   & 
         i, j           & ! horizontal indices
      ,  ij             & ! horizontal index, combines i and j loops
      ,  k               ! ice layer index

      real (kind=dbl_kind) ::   &
         ki              ! thermal cond of sea ice (W m-1 deg-1)

      khi(:,:,:) = c0i
      khs(:,:) = c0i

      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)

         ! top surface of top ice layer
         ki = kice + betak*salin(1) / min(-puny,Tin(i,j,1))
         ki = max (ki, kimin)
         khi(i,j,1) = ki / hlyr(i,j)

         ! bottom surface of bottom ice layer
         ki = kice + betak*salin(nilyr+1) / min(-puny,Tbot(i,j)) 
         ki = max (ki, kimin)
         khi(i,j,nilyr+1) = ki / hlyr(i,j)

         ! if snow is present: top snow surface, snow-ice interface
         if (hsn(i,j) > hsnomin) then
            khs(i,j) = ksno / hsn(i,j)
            khi(i,j,1) = c2i*khi(i,j,1)*khs(i,j) /  &
                        (khi(i,j,1) + khs(i,j))
         endif

      enddo                     ! ij

      ! interior ice interfaces
      do k = 2, nilyr
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)

            ki = kice + betak*p5*(salin(k-1)+salin(k)) /  &
                 min (-puny, p5*(Tin(i,j,k-1)+Tin(i,j,k)))
            ki = max (ki, kimin)
            khi(i,j,k) = ki / hlyr(i,j)

         enddo                  ! ij
      enddo                     ! k

      end subroutine conductivity

!=======================================================================
!BOP
!
! !ROUTINE: absorbed_solar - shortwave radiation absorbed by ice, ocean
!
! !DESCRIPTION:
!
! Compute solar radiation absorbed in ice and penetrating to ocean
!
! !REVISION HISTORY:
!
! authors William H. Lipscomb, LANL
!         C. M. Bitz, UW
!
! !INTERFACE:
!
      subroutine absorbed_solar                 &
          (ni,    icells,  indxi,  indxj,       &
           hlyr, hsn,     fswsfc, fswint, fswthrun, Iabs)
!
! !USES:
! 
      use ice_albedo
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) ::   &
         ni           &    ! thickness category index
      ,  icells          ! number of cells with aicen > puny

      integer (kind=int_kind), dimension(1:(ihi-ilo+1)*(jhi-jlo+1)), &
         intent(in) ::   & 
         indxi, indxj    ! compressed indices for cells with aicen > puny

      real (kind=dbl_kind), dimension(ilo:ihi,jlo:jhi), intent(in) ::   &
         hlyr         &   ! ice layer thickness
      ,  hsn             ! snow thickness (m)

      real (kind=dbl_kind), dimension(ilo:ihi,jlo:jhi), intent(inout) ::   &
         fswsfc         & ! SW absorbed at ice/snow surface (W m-2)
      ,  fswint         & ! SW absorbed in ice interior, below surface (W m-2)
      ,  fswthrun        ! SW through ice to ocean (W m-2)

      real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi,nilyr),       &
           intent(inout) ::   &
         Iabs            ! SW absorbed in particular layer (W m-2) 
!
!EOP
!

# if defined (ICE_FRESHWATER)
! EJA 201609212016 
! i0vis=0.7-->0.068 (Parkinson & Washington 1979, JGR) 17%*0.4=0.068
      real (kind=dbl_kind) ::   &
         i0vis                      ! fraction of penetrating solar rad (visible) - set in namelist
         !i0vis   = 0.068_dbl_kind  ! fraction of penetrating solar rad (visible) - set in namelist
# else
      real (kind=dbl_kind) ::   &
         i0vis                     ! fraction of penetrating solar rad (visible) - set in namelist
         !i0vis   = 0.70_dbl_kind  ! fraction of penetrating solar rad (visible) - set in namelist
# endif

      integer (kind=int_kind) ::   & 
         i, j           & ! horizontal indices
      ,  ij             & ! horizontal index, combines i and j loops
      ,  k               ! ice layer index

      real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi) ::   &
         fswpen        &  ! SW penetrating beneath surface (W m-2)
      ,  trantop       &  ! transmitted frac of penetrating SW at layer top
      ,  tranbot         ! transmitted frac of penetrating SW at layer bot

      real (kind=dbl_kind) ::   &
         swabs         &  ! net SW down at surface (W m-2)
      ,  swabsv        &  ! swabs in vis (wvlngth < 700nm)  (W/m^2)
      ,  swabsi        &  ! swabs in nir (wvlngth > 700nm)  (W/m^2)
      ,  frsnow          ! fractional snow coverage


      fswpen (:,:) = c0i
      trantop(:,:) = c0i
      tranbot(:,:) = c0i

      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)

      !-----------------------------------------------------------------
      ! Fractional snow cover
      !-----------------------------------------------------------------

         if (hsn(i,j) > puny) then
            frsnow = hsn(i,j) / (hsn(i,j) + snowpatch)
         else
            frsnow = c0i
         endif

      !-----------------------------------------------------------------
      ! Shortwave flux absorbed at surface, absorbed internally,
      !  and penetrating to mixed layer.  
      ! This parameterization assumes that all IR is absorbed at the
      !  surface; only visible is absorbed in the ice interior or
      !  transmitted to the ocean.
      !-----------------------------------------------------------------

         swabsv  = swvdr(i,j)*(c1i-alvdrn(i,j,ni))        &
                 + swvdf(i,j)*(c1i-alvdfn(i,j,ni))         
         swabsi  = swidr(i,j)*(c1i-alidrn(i,j,ni))        &
                 + swidf(i,j)*(c1i-alidfn(i,j,ni))
         swabs   = swabsv + swabsi

         fswpen(i,j) = swabsv * (c1i-frsnow) * i0vis
!c    &               + swabsi * (c1i-frsnow) * i0nir  ! i0nir = 0
         fswsfc(i,j) = swabs - fswpen(i,j)

         trantop(i,j) = c1i  ! transmittance at top of ice

      enddo                     ! ij

      !-----------------------------------------------------------------
      ! penetrating SW absorbed in each ice layer 
      !-----------------------------------------------------------------

      do k = 1, nilyr
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)

            tranbot(i,j) = exp(-kappav*hlyr(i,j)*real(k,kind=dbl_kind))
            Iabs(i,j,k) = fswpen(i,j) * (trantop(i,j)-tranbot(i,j))

            ! bottom of layer k = top of layer k+1
            trantop(i,j) = tranbot(i,j)

         enddo                  ! ij
      enddo                     ! k

      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)

         ! SW penetrating thru ice into ocean
         fswthrun(i,j) = fswpen(i,j) * tranbot(i,j)

         ! SW absorbed in ice interior
         fswint(i,j)  = fswpen(i,j) - fswthrun(i,j)

      enddo                     ! ij

      end subroutine absorbed_solar

!=======================================================================
!BOP
!
! !ROUTINE: surface_fluxes - surface radiative and turbulent fluxes
!
! !DESCRIPTION:
!
! Compute radiative and turbulent fluxes and their derivatives
! with respect to Tsf.  
!
! !REVISION HISTORY:
!
! authors William H. Lipscomb, LANL
!         C. M. Bitz, UW
!
! !INTERFACE:
!
      subroutine surface_fluxes                      &
          (isolve,     indxii,    indxjj,            &
           Tsf,        fswsfc,                       &
           flwoutn,    fsensn,    flatn,    fsurf,   &
           dflwout_dT, dfsens_dT, dflat_dT, dfsurf_dT) 
!
! !USES:
! 
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) ::   & 
         isolve          ! number of cells with temps not converged

      integer (kind=int_kind), dimension(1:(ihi-ilo+1)*(jhi-jlo+1)),   &
         intent(in) ::   & 
         indxii, indxjj  ! compressed indices for cells not converged

      real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi), intent(in) :: &
         Tsf           &  ! ice/snow surface temperature, Tsfcn
      ,  fswsfc          ! SW absorbed at ice/snow surface (W m-2)

      real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi), intent(inout):: &
         fsensn         & ! surface downward sensible heat (W m-2)
      ,  flatn          & ! surface downward latent heat (W m-2)
      ,  flwoutn        & ! upward LW at surface (W m-2)
      ,  fsurf          & ! net flux to top surface, not including fcondtop
      ,  dfsens_dT      & ! deriv of fsens wrt Tsf (W m-2 deg-1)
      ,  dflat_dT       & ! deriv of flat wrt Tsf (W m-2 deg-1)
      ,  dflwout_dT     & ! deriv of flwout wrt Tsf (W m-2 deg-1)
      ,  dfsurf_dT       ! derivative of fsurf wrt Tsf
!
!EOP
!
      integer (kind=int_kind) ::   & 
         i, j           & ! horizontal indices
      ,  ij             & ! horizontal index, combines i and j loops
      ,  k                ! ice layer index

      real (kind=dbl_kind) ::   &
         TsfK          &  ! ice/snow surface temperature (K)
      ,  Qsfc          &  ! saturated surface specific humidity (kg/kg)
      ,  dQsfcdT       &  ! derivative of Qsfc wrt surface temperature
      ,  qsat          &  ! the saturation humidity of air (kg/m^3)
      ,  flwdabs       &  ! downward longwave absorbed heat flx (W/m^2)
      ,  tmpvar           ! 1/TsfK

      do ij = 1, isolve
         i = indxii(ij)         ! NOTE: not indxi and indxj
         j = indxjj(ij)


         ! ice surface temperature in Kelvin
         TsfK = Tsf(i,j) + Tffresh
         tmpvar = c1i/TsfK

         ! saturation humidity 
         qsat    = qqqice * exp(-TTTice*tmpvar)
         Qsfc    = qsat / rhoa(i,j)
         dQsfcdT = TTTice * tmpvar*tmpvar * Qsfc

         ! longwave radiative flux
         flwdabs =  emissivity * flw(i,j)
         flwoutn(i,j) = -emissivity * stefan_boltzmann * TsfK**4

         ! downward latent and sensible heat fluxes
         fsensn(i,j) = shcoef(i,j) * (potT(i,j) - TsfK)
         flatn(i,j)  = lhcoef(i,j) * (Qa(i,j) - Qsfc)

         ! derivatives wrt surface temp
         dflwout_dT(i,j) = - emissivity*stefan_boltzmann * c4i*TsfK**3 
         dfsens_dT(i,j)  = - shcoef(i,j)
         dflat_dT(i,j)   = - lhcoef(i,j) * dQsfcdT
         
         fsurf(i,j) = fswsfc(i,j) + flwdabs + flwoutn(i,j)      &
                  + fsensn(i,j) + flatn(i,j)                     
         dfsurf_dT(i,j) = dflwout_dT(i,j)                       &
                         + dfsens_dT(i,j) + dflat_dT(i,j)
      enddo                     

      end subroutine surface_fluxes

!=======================================================================
!BOP
!
! !ROUTINE: tridiag_solver - tridiagonal matrix solver
!
! !DESCRIPTION:
!
! Tridiagonal matrix solver--used to solve the implicit vertical heat 
! equation in ice and snow
!
! !REVISION HISTORY:
!
! authors William H. Lipscomb, LANL
!         C. M. Bitz, UW
!
! !INTERFACE:
!
      subroutine tridiag_solver         &
          (isolve, indxii, indxjj,      &
           nmat,   sbdiag, diag,  spdiag, rhs, xout)
!
! !USES:
! 
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) ::   & 
         isolve          ! number of cells with temps not converged

      integer (kind=int_kind), dimension(1:(ihi-ilo+1)*(jhi-jlo+1)),&
         intent(in) ::   & 
         indxii, indxjj  ! compressed indices for cells not converged

      integer (kind=int_kind), intent(in) ::   &
         nmat            ! matrix dimension

      real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi,nmat), &
           intent(in) ::   &
         diag           & ! diagonal matrix elements
      ,  sbdiag         & ! sub-diagonal matrix elements
      ,  spdiag         & ! super-diagonal matrix elements
      ,  rhs              ! rhs of tri-diagonal matrix eqn.

      real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi,nmat),        &
           intent(out) ::   &
         xout            ! solution vector
!
!EOP
!
      integer (kind=int_kind) ::   & 
         i, j          & ! horizontal indices
      ,  ij            & ! horizontal index, combines i and j loops
      ,  k               ! row counter

      real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi) ::   &
         wbeta           ! temporary matrix variable

      real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi,nmat) ::   &
         wgamma          ! temporary matrix variable


      do ij = 1, isolve
         i = indxii(ij)
         j = indxjj(ij)

         wbeta(i,j) = diag(i,j,1)
         xout(i,j,1) = rhs(i,j,1) / wbeta(i,j)
         
      enddo                     ! ij

      do k = 2, nmat
         do ij = 1, isolve
            i = indxii(ij)
            j = indxjj(ij)

            wgamma(i,j,k) = spdiag(i,j,k-1) / wbeta(i,j)
            wbeta(i,j) = diag(i,j,k) - sbdiag(i,j,k)*wgamma(i,j,k)
            xout(i,j,k) = (rhs(i,j,k) - sbdiag(i,j,k)*xout(i,j,k-1))  &
                         / wbeta(i,j)

         enddo                  ! ij
      enddo                     ! k

      do k = nmat-1, 1, -1
         do ij = 1, isolve
            i = indxii(ij)
            j = indxjj(ij)

            xout(i,j,k) = xout(i,j,k) - wgamma(i,j,k+1)*xout(i,j,k+1)
            
         enddo                  ! ij
      enddo                     ! k

      end subroutine tridiag_solver

!=======================================================================
!BOP
!
! !ROUTINE: thickness changes - top and bottom growth/melting
!
! !DESCRIPTION:
!
! Compute growth and/or melting at the top and bottom surfaces.
!
! !REVISION HISTORY:
!
! authors William H. Lipscomb, LANL
!         C. M. Bitz, UW
!
! !INTERFACE:
!
      subroutine thickness_changes                      &
          (ni,     icells,   indxi,    indxj,            &
           hin,   hlyr,     hsn,      qin,    qsn,      &
           mvap,  Tbot,     fbot,     flatn,  fhnetn,   &
           fsurf, fcondtop, fcondbot, efinal) 
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) ::   &
         ni            &  ! thickness category index
      ,  icells          ! number of cells with aicen > puny

      integer (kind=int_kind), dimension(1:(ihi-ilo+1)*(jhi-jlo+1)),       &
         intent(in) ::   & 
         indxi, indxj    ! compressed indices for cells with aicen > puny

      real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi), intent(in) ::   &
         fbot           & ! ice-ocean heat flux at bottom surface (W/m^2)
      ,  Tbot           & ! ice bottom surface temperature (deg C)
      ,  flatn          & ! surface downward latent heat (W m-2)
      ,  fsurf          & ! net flux to top surface, not including fcondtop
      ,  fcondtop       & ! downward cond flux at top surface (W m-2)
      ,  fcondbot        ! downward cond flux at bottom surface (W m-2)

      real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi,nilyr),       &
         intent(inout) ::   &
         qin             ! ice layer enthalpy (J m-3)

      real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi), intent(inout)::   &
         fhnetn        & ! fbot, corrected for any surplus energy
      ,  hlyr          & ! ice layer thickness
      ,  hsn           & ! snow thickness (m)
      ,  qsn           & ! snow enthalpy
      ,  efinal        & ! final energy of melting (J m-2)
      ,  mvap            ! ice/snow mass sublimated/condensed (kg m-2) 

      real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi), intent(out) ::   &
         hin             ! total ice thickness
!
!EOP
!
      real (kind=dbl_kind), parameter ::   &
         qbotmax = -p5*rhoi*Lfresh  ! max enthalpy of ice growing at bottom

      integer (kind=int_kind) ::   & 
         i, j          &  ! horizontal indices
      ,  ij            &  ! horizontal index, combines i and j loops
      ,  k, kold         ! ice layer indices

      real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi) ::   &
         esub            &! energy for sublimation, > 0    (J m-2)
      ,  etop_mlt        &! energy for top melting, > 0    (J m-2)
      ,  ebot_mlt        &! energy for bottom melting, > 0 (J m-2)
      ,  rhlyr            ! reciprocal ice layer thickness

      real (kind=dbl_kind) ::   &
         dhi            & ! change in ice thickness
      ,  dhs            & ! change in snow thickness
      ,  Ti             & ! ice temperature
      ,  Ts             & ! snow temperature
      ,  econ           & ! energy for condensation, < 0   (J m-2)
      ,  ebot_gro       & ! energy for bottom growth, < 0  (J m-2)
      ,  qbot           & ! enthalpy of ice growing at bottom surface (J m-3)
      ,  qsub           & ! energy/unit volume to sublimate ice/snow (J m-3) 
      ,  hqtot          & ! sum of h*q for two layers
      ,  w1             & ! temporary variable
      ,  hovlp           ! overlap between old and new layers (m)

      real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi,nilyr) ::   &
         hnew          &  ! new layer thickness
      ,  hq              ! h * q for a layer

      real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi,0:nilyr) ::   &
         zold          &  ! depth of layer boundaries (m)
      ,  znew            ! adjusted depths, with equal hlyr (m)

      !-----------------------------------------------------------------
      ! Initialize ice layer thicknesses 
      !-----------------------------------------------------------------
     
      do k = 1, nilyr
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)
            hnew(i,j,k) = hlyr(i,j)
         enddo  ! if
      enddo

      !-----------------------------------------------------------------
      ! For l_brine = false (fresh ice), check for temperatures > 0.
      !  Melt ice or snow as needed to bring temperatures back to 0.
      ! For l_brine = true, this should not be necessary.
      !-----------------------------------------------------------------

      if (.not. l_brine) then 
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)

            Ts = (Lfresh + qsn(i,j)/rhos) / cp_ice
            if (Ts > c0i) then
               dhs = cp_ice*Ts*hsn(i,j) / Lfresh
               hsn(i,j) = hsn(i,j) - dhs
               qsn(i,j) = -rhos*Lfresh
            endif
         enddo

         do k = 1, nilyr
            do ij = 1, icells
               i = indxi(ij)
               j = indxj(ij)

               Ti = (Lfresh + qin(i,j,k)/rhoi) / cp_ice
               if (Ti > c0i) then
                  dhi = cp_ice*Ti*hnew(i,j,k) / Lfresh
                  hnew(i,j,k) = hnew(i,j,k) - dhi
                  qin(i,j,k) = -rhoi*Lfresh
               endif
            enddo               ! ij
         enddo                  ! k

      endif                     ! .not. l_brine

      !-----------------------------------------------------------------
      ! Sublimation/condensation of ice/snow
      !-----------------------------------------------------------------

      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)

!         w1 = -flatn(i,j) * dt          
         w1 = -flatn(i,j) * dtice          
         esub(i,j) = max(w1, c0i)    ! energy for sublimation, > 0
         econ      = min(w1, c0i)    ! energy for condensation, < 0

         !--------------------------------------------------------------
         ! Condensation (mvap > 0)
         !--------------------------------------------------------------

         if (hsn(i,j) > puny) then    ! add snow with enthalpy qsn
            dhs = econ / (qsn(i,j) - rhos*Lvap)      ! econ < 0, dhs > 0
            hsn(i,j) = hsn(i,j) + dhs
            mvap(i,j) = mvap(i,j) + dhs*rhos
         else                   ! add ice with enthalpy qin(i,j,1)
            dhi = econ / (qin(i,j,1) - rhoi*Lvap)    ! econ < 0, dhi > 0  
            hnew(i,j,1) = hnew(i,j,1) + dhi
            mvap(i,j) = mvap(i,j) + dhi*rhoi
         endif

         !--------------------------------------------------------------
         ! Sublimation of snow (mvap < 0)
         !--------------------------------------------------------------
         qsub = qsn(i,j) - rhos*Lvap                 ! qsub < 0
         dhs  = max (-hsn(i,j), esub(i,j)/qsub)      ! esub > 0, dhs < 0
         hsn(i,j) = hsn(i,j) + dhs
         esub(i,j) = esub(i,j) - dhs*qsub
         esub(i,j) = max(esub(i,j), c0i)   ! in case of roundoff error
         mvap(i,j) = mvap(i,j) + dhs*rhos               
      enddo                     ! ij

         !--------------------------------------------------------------
         ! Sublimation of ice (mvap < 0)
         !--------------------------------------------------------------

      do k = 1, nilyr
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)
            qsub = qin(i,j,k) - rhoi*Lvap            ! qsub < 0
            dhi = max (-hnew(i,j,k), esub(i,j)/qsub) ! esub < 0, dhi < 0
            hnew(i,j,k) = hnew(i,j,k) + dhi   
            esub(i,j) = esub(i,j) - dhi*qsub  
            esub(i,j) = max(esub(i,j), c0i)
            mvap(i,j) = mvap(i,j) + dhi*rhoi
         enddo                  ! ij
      enddo                     ! k

      !-----------------------------------------------------------------
      ! Top melt 
      ! (There is no top growth because there is no liquid water to
      !  freeze at the top.)
      !-----------------------------------------------------------------

         !--------------------------------------------------------------
         ! Melt snow
         !--------------------------------------------------------------
      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)

!         w1 = (fsurf(i,j) - fcondtop(i,j)) * dt
         w1 = (fsurf(i,j) - fcondtop(i,j)) * dtice
         etop_mlt(i,j) = max(w1, c0i)                  ! etop_mlt > 0
         dhs = max(-hsn(i,j), etop_mlt(i,j)/qsn(i,j)) ! qsn < 0, dhs < 0 
         hsn(i,j) = hsn(i,j) + dhs
         etop_mlt(i,j) = etop_mlt(i,j) - dhs*qsn(i,j)
         etop_mlt(i,j) = max(etop_mlt(i,j), c0i) ! in case of roundoff error

         ! history diagnostics
         if (dhs < -puny .and. mlt_onset(i,j) < puny)      &
              mlt_onset(i,j) = yday

      enddo                     ! ij

         !--------------------------------------------------------------
         ! Melt ice
         !--------------------------------------------------------------
         
      do k = 1, nilyr
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)

            dhi = max(-hnew(i,j,k), etop_mlt(i,j)/qin(i,j,k)) 
            hnew(i,j,k) = hnew(i,j,k) + dhi           ! qin < 0, dhi < 0
            etop_mlt(i,j) = etop_mlt(i,j) - dhi*qin(i,j,k)
            etop_mlt(i,j) = max(etop_mlt(i,j), c0i)

            ! history diagnostics
            if (dhi < -puny .and. mlt_onset(i,j) < puny) &
                 mlt_onset(i,j) = yday
            meltt(i,j) = meltt(i,j) - dhi*aicen(i,j,ni)  

         enddo                  ! ij
      enddo                     ! k

       !----------------------------------------------------------------
       ! Bottom growth and melt
       !----------------------------------------------------------------

      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)

!         w1 = (fcondbot(i,j) - fbot(i,j)) * dt
         w1 = (fcondbot(i,j) - fbot(i,j)) * dtice
         ebot_mlt(i,j) = max(w1, c0i)           ! ebot_mlt > 0
         ebot_gro = min(w1, c0i)                ! ebot_gro < 0

         !--------------------------------------------------------------
         ! Grow ice
         !--------------------------------------------------------------
   
         ! enthalpy of new ice growing at bottom surface
# if defined (ICE_FRESHWATER)
! afm 20151112 & EJA 20160921
! put if statement to avoid zero division for fresh water ice
         if (Tbot(i,j) /= c0i ) then
                qbot = -rhoi * (cp_ice * (Tmlt(nilyr+1)-Tbot(i,j))           &
                       + Lfresh * (c1i-Tmlt(nilyr+1)/Tbot(i,j))        &
                       - cp_ocn * Tmlt(nilyr+1)) 
         else
                qbot = -rhoi * Lfresh
         endif 
# else
         qbot = -rhoi * (cp_ice * (Tmlt(nilyr+1)-Tbot(i,j))           &
                       + Lfresh * (c1i-Tmlt(nilyr+1)/Tbot(i,j))        &
                       - cp_ocn * Tmlt(nilyr+1))                       
# endif
                
         qbot = min (qbot, qbotmax) ! in case Tbot is close to Tmlt    
         dhi  = ebot_gro / qbot     ! dhi > 0                          
                                                                       
         hqtot = hnew(i,j,nilyr)*qin(i,j,nilyr) + dhi*qbot             
         hnew(i,j,nilyr) = hnew(i,j,nilyr) + dhi                       
                                                                       
         if (hnew(i,j,nilyr) > puny)                                  &
              qin(i,j,nilyr) = hqtot / hnew(i,j,nilyr)                 
                                                                       
         ! history diagnostics                                         
         congel(i,j) = congel(i,j) + dhi*aicen(i,j,ni)                 
         if (dhi > puny .and. frz_onset(i,j) < puny)                  &
                 frz_onset(i,j) = yday

      enddo                     ! ij


         !--------------------------------------------------------------
         ! Melt ice
         !--------------------------------------------------------------
              
      do k = nilyr, 1, -1
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)

            dhi = max(-hnew(i,j,k), ebot_mlt(i,j)/qin(i,j,k))
            hnew(i,j,k) = hnew(i,j,k) + dhi           ! qin < 0, dhi < 0
            ebot_mlt(i,j) = ebot_mlt(i,j) - dhi*qin(i,j,k)
            ebot_mlt(i,j) = max(ebot_mlt(i,j), c0i)

            ! history diagnostics
            meltb(i,j) = meltb(i,j) - dhi*aicen(i,j,ni)   

         enddo                  ! ij
      enddo                     ! k

         !--------------------------------------------------------------
         ! Melt snow (only if all the ice has melted)
         !--------------------------------------------------------------

      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)

         dhs = max(-hsn(i,j), ebot_mlt(i,j)/qsn(i,j))
         hsn(i,j) = hsn(i,j) + dhs                    ! qsn < 0, dhs < 0
         ebot_mlt(i,j) = ebot_mlt(i,j) - dhs*qsn(i,j)
         ebot_mlt(i,j) = max(ebot_mlt(i,j), c0i)

         ! snow-to-ice conversion (very rare)
         if (hsn(i,j) > puny .and. ebot_mlt(i,j) > puny*Lfresh) then 
            hnew(i,j,1) = hsn(i,j) * rhos/rhoi ! reduce thickness
            qin(i,j,1)  = qsn(i,j) * rhoi/rhos ! increase q to conserve energy
            hsn(i,j) = c0i
         endif

      enddo                     ! ij

      !-----------------------------------------------------------------
      ! Give the ocean any energy left over after melting
      !-----------------------------------------------------------------

      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)

         fhnetn(i,j) = fhnetn(i,j)    &
                     + (esub(i,j) + etop_mlt(i,j) + ebot_mlt(i,j))/dtice
!                     + (esub(i,j) + etop_mlt(i,j) + ebot_mlt(i,j))/dt

      enddo                     ! ij

!---!-------------------------------------------------------------------
!---! Repartition the ice into equal-thickness layers, conserving energy.
!---!-------------------------------------------------------------------

      !-----------------------------------------------------------------
      ! Initialize the new ice thickness.
      ! Initialize the final energy: snow energy plus energy of 
      !  sublimated/condensed ice.
      !-----------------------------------------------------------------

      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)
         hin(i,j) = c0i
         efinal(i,j) = hsn(i,j)*qsn(i,j) - mvap(i,j)*Lvap
      enddo

      !-----------------------------------------------------------------
      ! Compute the new ice thickness.
      ! Initialize h*q for new layers.
      !-----------------------------------------------------------------

      do k = 1, nilyr
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)
# if defined (ICE_FRESHWATER)
! afm 20151114 & EJA 20160921 isnan check
if ( isnan( hin(i,j) ) /= 0 .or. isnan( hnew(i,j,k) ) /= 0 ) then
print*, 'hin isnan:', i, k,hlyr(i,j), hnew(i,j,k), dhi, dhs, qbot, &
qsn(i,j), hsn(i,j), esub(i,j),etop_mlt(i,j), ebot_mlt(i,j), ebot_gro 
end if
# endif

            hin(i,j) = hin(i,j) + hnew(i,j,k)
            hq(i,j,k) = c0i
         enddo                  ! ij
      enddo                     ! k

      !-----------------------------------------------------------------
      ! Make sure hin > puny.
      ! Compute depths zold of old layers (unequal thicknesses).
      ! Compute depths znew of new layers (all the same thickness). 
      !-----------------------------------------------------------------

      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)

         if (hin(i,j) > puny) then
            hlyr(i,j) = hin(i,j) / real(nilyr,kind=dbl_kind)
            rhlyr(i,j) = c1i / hlyr(i,j)
         else
            hin(i,j) = c0i
            hlyr(i,j) = c0i
            rhlyr(i,j) = c0i
         endif

         zold(i,j,0) = c0i
         znew(i,j,0) = c0i

         zold(i,j,nilyr) = hin(i,j)
         znew(i,j,nilyr) = hin(i,j)

      enddo                     ! ij

      do k = 1, nilyr-1
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)

            zold(i,j,k) = zold(i,j,k-1) + hnew(i,j,k) ! old unequal layers
            znew(i,j,k) = znew(i,j,k-1) + hlyr(i,j)   ! new equal layers

         end do                 ! ij
      enddo                     ! k
              
      !-----------------------------------------------------------------
      ! Compute h*q for new layer (k) given overlap with old layers (kold)
      !-----------------------------------------------------------------

      do k = 1, nilyr
         do kold = 1, nilyr
            do ij = 1, icells
               i = indxi(ij)
               j = indxj(ij)

               hovlp = min (zold(i,j,kold),   znew(i,j,k))    &
                     - max (zold(i,j,kold-1), znew(i,j,k-1))
               hovlp = max (hovlp, c0i)

               hq(i,j,k) = hq(i,j,k) + hovlp*qin(i,j,kold)

            enddo               ! ij
         enddo                  ! kold
      enddo                     ! k

      !-----------------------------------------------------------------
      ! Compute new values of qin and increment ice-snow energy.
      ! Note: Have to finish previous loop before updating qin
      !       in this loop.
      !-----------------------------------------------------------------

      do k = 1, nilyr
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)

            qin(i,j,k) = hq(i,j,k) * rhlyr(i,j)
            efinal(i,j) = efinal(i,j) + hlyr(i,j)*qin(i,j,k)

         enddo                  ! ij
      enddo                     ! 

      end subroutine thickness_changes

!=======================================================================
!BOP
!
! !ROUTINE: conservation_check_vthermo - energy conservation check
!
! !DESCRIPTION:
!
! Check for energy conservation by comparing the change in energy 
! to the net energy input.
!
! !REVISION HISTORY:
!
! authors William H. Lipscomb, LANL
!         C. M. Bitz, UW
!
! !INTERFACE:
!
      subroutine conservation_check_vthermo        & 
          (ni,     icells, indxi,  indxj,          & 
           fsurf, flatn,  fhnetn, fswint, einit, efinal)
!
! !USES:
! 
!      use ice_exit
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) ::   &
         ni             &  ! thickness category index (diagnostic only)
      ,  icells          ! number of cells with aicen > puny

      integer (kind=int_kind), dimension(1:(ihi-ilo+1)*(jhi-jlo+1)),      &
         intent(in) ::   & 
         indxi, indxj    ! compressed indices for cells with aicen > puny

      real (kind=dbl_kind), dimension (ilo:ihi, jlo:jhi), intent(in) ::   &
         fsurf          & ! net flux to top surface, not including fcondtop
      ,  flatn          & ! surface downward latent heat (W m-2)
      ,  fhnetn         & ! fbot, corrected for any surplus energy
      ,  fswint         & ! SW absorbed in ice interior, below surface (W m-2)
      ,  einit          & ! initial energy of melting (J m-2)
      ,  efinal          ! final energy of melting (J m-2)
!
!EOP
!
      integer (kind=int_kind) ::   & 
         i, j         &   ! horizontal indices
      ,  ij              ! horizontal index, combines i and j loops

      real (kind=dbl_kind) ::   &
         einp          & ! energy input during timestep (J m-2)
      ,  ferr            ! energy conservation error (W m-2)

      logical (kind=log_kind) ::   &   ! for vector-friendly error checks
         ferr_flag       ! flag for energy error, ferr > ferrmax

      ferr_flag = .false.     ! ferr <= ferrmax, initialization

      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)

         einp = (fsurf(i,j) - flatn(i,j) + fswint(i,j) - fhnetn(i,j)) &
              * dtice
!              * dt
!         ferr = abs(efinal(i,j)-einit(i,j)-einp) / dt
         ferr = abs(efinal(i,j)-einit(i,j)-einp) / dtice
         if (ferr > ferrmax) then
           ferr_flag = .true.
         endif
      enddo
      !----------------------------------------------------------------
      ! If energy is not conserved, print diagnostics and exit.
      !----------------------------------------------------------------
      if (ferr_flag .and. DBG_SET(DBG_SBRIO)) then
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)

      !-----------------------------------------------------------------
      ! Note that fsurf - flat = fsw + flw + fsens; i.e., the latent
      ! heat is not included in the energy input, since de is the energy 
      ! change in the system ice + vapor, and the latent heat lost by 
      ! the ice is equal to that gained by the vapor.
      !-----------------------------------------------------------------

            einp = (fsurf(i,j) - flatn(i,j) + fswint(i,j) -         &
                    fhnetn(i,j)) * dtice                           
!                    fhnetn(i,j)) * dt                                
!            ferr = abs(efinal(i,j)-einit(i,j)-einp) / dt             
            ferr = abs(efinal(i,j)-einit(i,j)-einp) / dtice             
            if (ferr > ferrmax) then                                 
!               write(nu_diag,*) 'Energy error, cat',ni,              &
!                    'my_task,i,j',my_task,i,j                        
!               write(nu_diag,*)'wind(I,j)',wind(i,j),shcoef(i,j), potT(i,j)
!               write(nu_diag,*) 'Flux error (W/m^2) =', ferr         
!               write(nu_diag,*) 'Energy error (J) =', ferr*dt        
!               write(nu_diag,*) 'Energy error (J) =', ferr*dtice        
!               write(nu_diag,*) 'Initial energy =', einit(i,j)       
!               write(nu_diag,*) 'Final energy =', efinal(i,j)        
!               write(nu_diag,*) 'efinal - einit =',                 &
!                                 efinal(i,j)-einit(i,j)              
!               write(nu_diag,*) 'Input energy =', einp               
!               stoplabel = 'ice state at stop'                       
!               call print_state(stoplabel,i,j)                       
!               call abort_ice                                       &
!                    ('vertical thermo: energy conservation error')
            endif
         enddo
      endif

      end subroutine conservation_check_vthermo

!=======================================================================
!BOP
!
! !ROUTINE: add_new_snow - add new snow
!
! !DESCRIPTION:
!
! Add new snow at top surface
!
! !REVISION HISTORY:
!
! authors William H. Lipscomb, LANL
!         C. M. Bitz, UW
!
! !INTERFACE:
!
      subroutine add_new_snow             &
          (ni,    icells,  indxi, indxj,  &
           hsn,  hsn_new, qsn,   Tsf)
!
! !USES:
! 
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) ::   &
         ni           &    ! thickness category index
      ,  icells          ! number of cells with aicen > puny

      integer (kind=int_kind), dimension(1:(ihi-ilo+1)*(jhi-jlo+1)), &
         intent(in) ::   & 
         indxi, indxj    ! compressed indices for cells with aicen > puny

      real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi), intent(in)::   &
         Tsf             ! ice/snow surface temperature, Tsfcn

      real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi), intent(inout)::   &
         hsn          &   ! snow thickness (m)
      ,  qsn             ! snow enthalpy

      real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi), intent(out) ::   &
         hsn_new        ! thickness of new snow (m) 
!
!EOP
!
      real (kind=dbl_kind) ::   &
         qsnew       &   ! enthalpy of new snow (J kg-1) 
      ,  hstot           ! snow thickness including new snow (J kg-1)

      integer (kind=int_kind) ::   & 
         i, j         &  ! horizontal indices
      ,  ij              ! horizontal index, combines i and j loops


      hsn_new(:,:) = c0i

      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)

      !----------------------------------------------------------------
      ! NOTE: If heat flux diagnostics are to work, new snow should
      !       have T = 0 (i.e. q = -rhos*Lfresh) and should not be
      !       converted to rain.
      !----------------------------------------------------------------

         if (fsnow(i,j) > c0i) then

!            hsn_new(i,j) = fsnow(i,j)/rhos * dt
            hsn_new(i,j) = fsnow(i,j)/rhos * dtice
            qsnew = -rhos*Lfresh
            hstot = hsn(i,j) + hsn_new(i,j)
            if (hstot > puny) then
               qsn(i,j) = (hsn(i,j)    *qsn(i,j)   &
                         + hsn_new(i,j)*qsnew   ) / hstot
               ! avoid roundoff errors if hsn=0
               qsn(i,j) = min(qsn(i,j), -rhos*Lfresh) 
               hsn(i,j) = hstot
            endif
         endif

      enddo                     ! ij

      end subroutine add_new_snow

!=======================================================================
!BOP
!
! !ROUTINE: update_state_vthermo - new state variables
!
! !DESCRIPTION:
!
! Given the vertical thermo state variables (hin, hsn, qin,
!  qsn, Tsf), compute the new ice state variables (vicen, vsnon, 
!  eicen, esnon, Tsfcn).
! Zero out state variables if ice has melted entirely.
!
! !REVISION HISTORY:
!
! authors William H. Lipscomb, LANL
!         C. M. Bitz, UW
!
! !INTERFACE:
!
      subroutine update_state_vthermo  & 
         (ni,    icells, indxi, indxj, &
          hin,  hsn,    qin,   qsn,   Tsf)
!
! !USES:
! 
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) ::   &
         ni             &  ! thickness category index
      ,  icells          ! number of cells with aicen > puny

      integer (kind=int_kind), dimension(1:(ihi-ilo+1)*(jhi-jlo+1)),         &
         intent(in) ::   & 
         indxi, indxj    ! compressed indices for cells with aicen > puny


      real (kind=dbl_kind), dimension(ilo:ihi,jlo:jhi), intent(in) ::   & 
         hin           &  ! ice thickness (m)
      ,  hsn           &  ! snow thickness (m)
      ,  qsn           &  ! snow enthalpy
      ,  Tsf             ! ice/snow surface temperature, Tsfcn

      real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi,nilyr),    &
         intent(in) ::   &
         qin             ! ice layer enthalpy (J m-3)
!
!EOP
!
      integer (kind=int_kind) ::   & 
         i, j         &   ! horizontal indices
      ,  ij           &   ! horizontal index, combines i and j loops
      ,  k               ! ice layer index

      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)

         if (hin(i,j) > c0i) then
            ! aicen is already up to date
            vicen(i,j,ni) = aicen(i,j,ni) * hin(i,j)
            vsnon(i,j,ni) = aicen(i,j,ni) * hsn(i,j)
            esnon(i,j,ni) = qsn(i,j) * vsnon(i,j,ni)
            Tsfcn(i,j,ni) = Tsf(i,j)
         else  ! (hin(i,j) == c0i)
            aicen(i,j,ni) = c0i
            vicen(i,j,ni) = c0i
            vsnon(i,j,ni) = c0i
            esnon(i,j,ni) = c0i
            Tsfcn(i,j,ni) = Tf(i,j)
         endif

      enddo                     ! ij

      do k = 1, nilyr
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)

            if (hin(i,j) > c0i) then
               eicen(i,j,ilyr1(ni)+k-1) =                    &
                    qin(i,j,k) * vicen(i,j,ni)/real(nilyr,kind=dbl_kind)  
            else
               eicen(i,j,ilyr1(ni)+k-1) = c0i
            endif

         enddo                  ! ij
      enddo                     ! k

      end subroutine update_state_vthermo

!=======================================================================

      end module ice_therm_vertical

!=======================================================================
