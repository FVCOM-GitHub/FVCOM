!/===========================================================================/
! CVS VERSION INFORMATION
! $Id$
! $Name$
! $Revision$
!/===========================================================================/

!=======================================================================
!BOP
!
! !MODULE: ice_therm_itd - thermo calculations after call to coupler
!
! !DESCRIPTION:
!
! Thermo calculations after call to coupler, mostly related to ITD:
! ice thickness redistribution, lateral growth and melting, and
! freeboard adjustment.
! 
! NOTE: The thermodynamic calculation is split in two for load balancing. 
!       First ice\_therm\_vertical computes vertical growth rates and coupler
!       fluxes.  Then ice\_therm\_itd does thermodynamic calculations not 
!       needed for coupling.
!
! !REVISION HISTORY:
!
! authors C. M. Bitz, UW
!         Elizabeth C. Hunke, LANL
!         William H. Lipscomb, LANL
!
! Vectorized by Clifford Chen (Fujitsu) and William H. Lipscomb (LANL)
!
! !INTERFACE:
!
      module ice_therm_itd
!
! !USES:
!
      use ice_kinds_mod
      use ice_model_size
      use ice_constants
      use ice_domain
      use ice_state
      use ice_flux
!      use ice_diagnostics
      use ice_calendar
      use ice_grid
      use ice_itd
!
!EOP
!
      implicit none
      save

!      real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi,ncat) ::
      real (kind=dbl_kind), dimension (:,:,:),allocatable,save :: &
     &   hicen           ! ice thickness (m)

!=======================================================================

      contains

!=======================================================================
!BOP
!
! !ROUTINE: thermo_itd - driver for post-coupler thermodynamics
!
! !DESCRIPTION:
!
!-----------------------------------------------------------------------
! Driver for thermodynamic changes not needed for coupling: 
! transport in thickness space, lateral growth and melting, and 
! freeboard adjustment.
!
! NOTE: Ocean fluxes are initialized here.
!
! !REVISION HISTORY:
!
! authors:     C. M. Bitz, UW
!              Elizabeth C. Hunke, LANL
!              William H. Lipscomb, LANL
!
! !INTERFACE:
!
      subroutine thermo_itd
!
! !USES:
! 
!      use ice_timers
      use ice_itd_linear
      use ice_therm_vertical, only: hicen_old, rside
!
!EOP
!
      integer (kind=int_kind) ::  &
     &   i, j          &  ! horizontal indices
     &,  ni             &  ! thickness category index
     &,  k               ! ice layer index

!      call ice_timer_start(4)  ! column model

      !-----------------------------------------------------------------
      ! Save the ice area passed to the coupler.
      ! This is needed to make the history fields consistent with
      !  the coupler fields.
      !-----------------------------------------------------------------
      aice_init = aice
   
      !-----------------------------------------------------------------
      ! Initialize ocean fluxes sent to the coupler.
      !-----------------------------------------------------------------
      call init_flux_ocn

      !-----------------------------------------------------------------
      ! Let rain drain through to the ocean. 
      !-----------------------------------------------------------------

      do j=jlo,jhi
      do i=ilo,ihi
         fresh(i,j)      = fresh(i,j)      + frain(i,j)*aice(i,j)
         fresh_hist(i,j) = fresh_hist(i,j) + frain(i,j)*aice(i,j)
      enddo
      enddo

      !-----------------------------------------------------------------
      ! Update ice thickness.
      !-----------------------------------------------------------------

!      call ice_timer_start(5)   ! thermodynamics
      do ni = 1, ncat
         do j=jlo,jhi
         do i=ilo,ihi
            if (aicen(i,j,ni) > puny) then
               hicen(i,j,ni) = vicen(i,j,ni) / aicen(i,j,ni)
            else
               hicen(i,j,ni) = c0i
               hicen_old(i,j,ni) = c0i
            endif
         enddo                  ! i
         enddo                  ! j
      enddo                     ! n
!      call ice_timer_stop(5)    ! thermodynamics

      !-----------------------------------------------------------------
      ! Given thermodynamic growth rates, transport ice between 
      ! thickness categories.
      !-----------------------------------------------------------------

!      call ice_timer_start(7)   ! category conversions (transport in h)
      if (kitd == 1) call linear_itd (hicen_old, hicen)
!      call ice_timer_stop(7)    ! category conversions 

      !-----------------------------------------------------------------
      ! Add frazil ice growing in leads. 
      !-----------------------------------------------------------------

!      call ice_timer_start(5)   ! thermodynamics
      call add_new_ice

      !-----------------------------------------------------------------
      ! Melt ice laterally.
      !-----------------------------------------------------------------
      call lateral_melt (rside)

      !-----------------------------------------------------------------
      ! Convert snow below freeboard to ice.
      !-----------------------------------------------------------------
      call freeboard
!      call ice_timer_stop(5)    ! thermodynamics

      !-----------------------------------------------------------------
      ! Make sure ice in each category is within its thickness bounds.
      ! NOTE: The rebin subroutine is needed only in the rare cases 
      !       when the linear_itd subroutine cannot transfer ice 
      !       correctly (e.g., very fast ice growth).
      !-----------------------------------------------------------------

!      call ice_timer_start(7)   ! category conversions
      if (ncat==1) then
         call reduce_area (hicen_old(:,:,1), hicen(:,:,1))
      else
         call rebin
      endif   ! ncat = 1
!      call ice_timer_stop(7)    ! category conversions 

      !-----------------------------------------------------------------
      ! Zero out ice categories with very small areas.
      !-----------------------------------------------------------------
      call zap_small_areas

      !-----------------------------------------------------------------
      ! Aggregate cell values over thickness categories. 
      !-----------------------------------------------------------------
      call aggregate

      !-----------------------------------------------------------------
      ! Compute thermodynamic area and volume tendencies.
      !-----------------------------------------------------------------

      do j=jlo,jhi
      do i=ilo,ihi
!         daidtt(i,j) = (aice(i,j) - daidtt(i,j)) / dt
!         dvidtt(i,j) = (vice(i,j) - dvidtt(i,j)) / dt
         daidtt(i,j) = (aice(i,j) - daidtt(i,j)) / dtice
         dvidtt(i,j) = (vice(i,j) - dvidtt(i,j)) / dtice

         daidtd(i,j) = aice(i,j) ! temporarily used for initial area
         dvidtd(i,j) = vice(i,j) ! temporarily used for initial volume
      enddo                     ! i
      enddo                     ! j

!      call ice_timer_stop(4)     ! column model

      end subroutine thermo_itd

!=======================================================================
!BOP
!
! !ROUTINE: add_new_ice - add frazil ice to ice thickness distribution 
!
! !DESCRIPTION:
!
! Given the volume of new ice grown in open water, compute its area
! and thickness and add it to the appropriate category or categories.
!
! NOTE: Usually all the new ice is added to category 1.  An exception is
!       made if there is no open water or if the new ice is too thick
!       for category 1, in which case ice is distributed evenly over the 
!       entire cell.  Subroutine rebin should be called in case the ice 
!       thickness lies outside category bounds after new ice formation.  
!
! !REVISION HISTORY:
!
! authors William H. Lipscomb, LANL
!         Elizabeth C. Hunke, LANL
!
! !INTERFACE:
!
      subroutine add_new_ice
!
! !USES:
! 
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) :: &
     &   i, j          &   ! horizontal indices
     &,  ni            &    ! ice category index
     &,  k                ! ice layer index

      real (kind=dbl_kind), dimension(ilo:ihi,jlo:jhi) :: &
     &   ai0new         &  ! area of new ice added to cat 1
     &,  vi0new         &  ! volume of new ice added to cat 1
     &,  hsurp          &  ! thickness of new ice added to each cat
     &,  vlyr             ! ice layer volume

      real (kind=dbl_kind), dimension (imt_local,jmt_local) :: &    
     &   vice_init, vice_final  ! ice volume summed over categories

      real (kind=dbl_kind) :: &
     &   fnew            & ! heat flx to open water for new ice (W/m^2)
     &,  hi0new          & ! thickness of new ice
     &,  hi0max          & ! max allowed thickness of new ice
     &,  qi0(nilyr)      & ! frazil ice enthalpy
     &,  qi0av           & ! mean value of qi0 for new ice (J kg-1)
     &,  vsurp           & ! volume of new ice added to each cat
     &,  area1           & ! starting fractional area of existing ice
     &,  rnilyr          & ! real(nilyr)
     &,  dfresh          & ! change in fresh
     &,  dfsalt          & ! change in fsalt
     &,  vi_frzmlt       & ! ice vol formed by frzmlt acting alone
     &,  vi_diff          ! vi0new - vi_frzmlt 

      real (kind=dbl_kind), parameter :: &
# if defined (ICE_FRESHWATER)
! afm 20151112 & EJA 20160921 hfrazilmin 5cm-->1cm.
     &   hfrazilmin = 0.01_dbl_kind  ! min thickness of new frazil ice (m)
# else
     &   hfrazilmin = 0.05_dbl_kind  ! min thickness of new frazil ice (m)
# endif

      integer (kind=int_kind) :: &
     &   icells, jcells, kcells &! grid cell counters
     &,  ij                     ! combined i/j horizontal index

      integer (kind=int_kind),                         &
     &        dimension (1:(ihi-ilo+1)*(jhi-jlo+1)) ::  &
     &   indxi,  indxj         & ! compressed i/j indices
     &,  indxi2, indxj2        &
     &,  indxi3, indxj3

      character (len=char_len) :: &
     &   fieldid           ! field identifier

      if (ncat > 1) then
         hi0max = hin_max(1)*0.9_dbl_kind  ! not too close to boundary
      else
         hi0max = 1.e8_dbl_kind            ! big number
      endif

      ! initial ice volume in each grid cell
      call column_sum (ncat, vicen, vice_init)

      !-----------------------------------------------------------------
      ! Compute average enthalpy of new ice.
      !
      ! POP assumes new ice is fresh.  Otherwise, it would be better
      ! to do something like this:
      !  qi0(i,j,k) = -rhoi * (cp_ice*(Tmlt(k)-Tf(i,j))
      !             + Lfresh*(1.-Tmlt(k)/Tf(i,j)) - cp_ocn*Tmlt(k))
      !-----------------------------------------------------------------

      rnilyr = real(nilyr,kind=dbl_kind)
      qi0av = c0i
      do k = 1, nilyr 
         qi0(k) = -rhoi*Lfresh  ! note sign convention, qi < 0
         qi0av  = qi0av + qi0(k)
      enddo
      qi0av = qi0av/rnilyr

      !-----------------------------------------------------------------
      ! Identify ice/ocean grid points.
      !-----------------------------------------------------------------
      icells = 0
      do j = jlo, jhi
      do i = ilo, ihi
         if (tmask(i,j)) then
            icells = icells + 1
            indxi(icells) = i
            indxj(icells) = j
         endif
      enddo       ! i
      enddo       ! j

      !-----------------------------------------------------------------
      ! Compute the volume, area, and thickness of new ice.
      !-----------------------------------------------------------------

      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)

         fnew = max (frzmlt(i,j), c0i)   ! fnew > 0 iff frzmlt > 0
         vi0new(i,j) = -fnew*dtice / qi0av ! note sign convention, qi < 0
!         vi0new(i,j) = -fnew*dt / qi0av ! note sign convention, qi < 0

         ! increment ice volume
         vice_init(i,j) = vice_init(i,j) + vi0new(i,j)

         ! history diagnostics
         frazil(i,j) = vi0new(i,j)      
         if (frazil(i,j) > puny .and. frz_onset(i,j) < puny)  &
     &           frz_onset(i,j) = yday

      !-----------------------------------------------------------------
      ! Update fresh water and salt fluxes.
      !
      ! NOTE: POP assumes fresh water and salt flux due to frzmlt > 0 
      !       is NOT included in fluxes fresh and fsalt.
      !-----------------------------------------------------------------

!!!         dfresh = -rhoi*vi0new(i,j)/dt  ! if POP had not already adjusted
                                           ! itself based on frzmlt
!!!         dfsalt = ice_ref_salinity*p001*dfresh

!!!         fresh(i,j)      = fresh(i,j)      + dfresh
!!!         fresh_hist(i,j) = fresh_hist(i,j) + dfresh
!!!         fsalt(i,j)      = fsalt(i,j)      + dfsalt
!!!         fsalt_hist(i,j) = fsalt_hist(i,j) + dfsalt

!!  There is no adjust in FVCOM when it is coupled with CICE model 
!!  be careful to adjust the fresh water balance
         !dfresh = -rhoi*vi0new(i,j)/dt  ! if POP had not already adjusted
         dfresh = -rhoi*vi0new(i,j)/dtice  ! if POP had not already adjusted
                                           ! itself based on frzmlt
         dfsalt = ice_ref_salinity*p001*dfresh

         fresh(i,j)      = fresh(i,j)      + dfresh
         fresh_hist(i,j) = fresh_hist(i,j) + dfresh
         fsalt(i,j)      = fsalt(i,j)      + dfsalt
         fsalt_hist(i,j) = fsalt_hist(i,j) + dfsalt





      !-----------------------------------------------------------------
      ! Decide how to distribute the new ice.
      !-----------------------------------------------------------------

         hsurp(i,j)  = c0i
         ai0new(i,j) = c0i

         if (vi0new(i,j) > c0i) then

            ! new ice area and thickness
            ! hin_max(0) < new ice thickness < hin_max(1)
            if (aice0(i,j) > puny) then
               hi0new = max(vi0new(i,j)/aice0(i,j), hfrazilmin)
               if (hi0new > hi0max .and. aice0(i,j)+puny < c1i) then
                  ! distribute excess volume over all categories (below)
                  hi0new = hi0max
                  ai0new(i,j) = aice0(i,j)
                  vsurp       = vi0new(i,j) - ai0new(i,j)*hi0new
                  hsurp(i,j)  = vsurp / aice(i,j)
                  vi0new(i,j) = ai0new(i,j)*hi0new
               else
                  ! put ice in a single category, with hsurp = 0
                  ai0new(i,j) = vi0new(i,j)/hi0new
               endif
            else                ! aice0 < puny
               hsurp(i,j) = vi0new(i,j)/aice(i,j) ! new thickness in each cat
               vi0new(i,j) = c0i
            endif               ! aice0 > puny
         endif                  ! vi0new > puny

      enddo                     ! ij

      !-----------------------------------------------------------------
      ! Identify grid cells receiving new ice.
      !-----------------------------------------------------------------
      jcells = 0
      kcells = 0

      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)

         if (vi0new(i,j) > c0i) then  ! add ice to category 1
            jcells = jcells + 1
            indxi2(jcells) = i
            indxj2(jcells) = j
         endif

         if (hsurp(i,j) > c0i) then   ! add ice to all categories 
            kcells = kcells + 1
            indxi3(kcells) = i
            indxj3(kcells) = j
         endif

      enddo

      !-----------------------------------------------------------------
      ! Distribute excess ice volume among ice categories by increasing
      ! ice thickness, leaving ice area unchanged.
      !-----------------------------------------------------------------

      do ni = 1, ncat

         do ij = 1, kcells
            i = indxi3(ij)
            j = indxj3(ij)

            vicen(i,j,ni) = vicen(i,j,ni) + aicen(i,j,ni)*hsurp(i,j)
            vlyr(i,j) = hsurp(i,j)/rnilyr * aicen(i,j,ni)
         enddo                  ! ij

         do k=1,nilyr
            do ij = 1, kcells
               i = indxi3(ij)
               j = indxj3(ij)

               eicen(i,j,ilyr1(ni)+k-1) =                      &
     &              eicen(i,j,ilyr1(ni)+k-1) + qi0(k)*vlyr(i,j) 
            enddo               ! ij
         enddo                  ! k

      enddo                     ! n
           
      !-----------------------------------------------------------------
      ! Combine new ice grown in open water with category 1 ice.
      ! NOTE: vsnon and esnon are unchanged.
      !-----------------------------------------------------------------
      do ij = 1, jcells
         i = indxi2(ij)
         j = indxj2(ij)

         area1 = aicen(i,j,1)   ! save
         aicen(i,j,1) = aicen(i,j,1) + ai0new(i,j)
         aice0(i,j)   = aice0(i,j)   - ai0new(i,j)
         vicen(i,j,1) = vicen(i,j,1) + vi0new(i,j)
         Tsfcn(i,j,1) = (Tf(i,j)*ai0new(i,j) + Tsfcn(i,j,1)*area1) &
     &                / aicen(i,j,1)
         Tsfcn(i,j,1) = min (Tsfcn(i,j,1), c0i)
         vlyr(i,j)    = vi0new(i,j)/rnilyr
      enddo                     ! ij

      do k = 1, nilyr
         do ij = 1, jcells
            i = indxi2(ij)
            j = indxj2(ij)
            eicen(i,j,k) = eicen(i,j,k) + qi0(k)*vlyr(i,j)
         enddo
      enddo
         
      call column_sum (ncat, vicen, vice_final)
      fieldid = 'vice, add_new_ice'
      call column_conservation_check(vice_init, vice_final,  &
     &                               puny, fieldid)

      end subroutine add_new_ice

!=======================================================================
!BOP
!
! !ROUTINE: lateral_melt - melt ice laterally
!
! !DESCRIPTION:
!
! Given the fraction of ice melting laterally in each grid cell 
!  (computed in subroutine frzmlt\_bottom\_lateral), melt ice.
!
! !REVISION HISTORY:
! 
! author:      C. M. Bitz, UW
! modified by: Elizabeth C. Hunke, LANL
!              William H. Lipscomb, LANL
!
! !INTERFACE:
!
      subroutine lateral_melt (rside)
!
! !USES:
! 
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind), dimension(ilo:ihi,jlo:jhi),   &
     &   intent(in) ::     &
     &   rside       ! fraction of ice that melts laterally
!
!EOP
!
      integer (kind=int_kind) ::   &
     &   i, j        & ! horizontal indices
     &,  ni          & ! thickness category index
     &,  k           &! layer index
     &,  ij          &! horizontal index, combines i and j loops
     &,  icells      ! number of cells with aice > puny

      integer (kind=int_kind), dimension(1:(ihi-ilo+1)*(jhi-jlo+1)) :: &
     &   indxi, indxj    ! compressed indices for cells with aice > puny

      real (kind=dbl_kind) ::  &
     &   dfhnet     & ! change in fhnet
     &,  dfresh     & ! change in fresh
     &,  dfsalt      ! change in fsalt

      do ni = 1, ncat

      !-----------------------------------------------------------------
      ! Identify grid cells with lateral melting. 
      !-----------------------------------------------------------------

         icells = 0
         do j = jlo, jhi
         do i = ilo, ihi
            if (rside(i,j) > c0i) then
               icells = icells + 1
               indxi(icells) = i
               indxj(icells) = j
            endif
         enddo                  ! i
         enddo                  ! j

      !-----------------------------------------------------------------
      ! Melt the ice and increment fluxes.
      !-----------------------------------------------------------------

         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)
               
            ! fluxes to coupler (except heat flux for ice melt)
            ! dfhnet < 0, dfresh > 0, dfsalt > 0
               
            dfhnet = esnon(i,j,ni)*rside(i,j) / dtice
!            dfhnet = esnon(i,j,n)*rside(i,j) / dt
            dfresh = (rhos*vsnon(i,j,ni) + rhoi*vicen(i,j,ni))  &
!     &             * rside(i,j) / dt
     &             * rside(i,j) / dtice
            dfsalt = rhoi*vicen(i,j,ni)*ice_ref_salinity*p001   &
     &             * rside(i,j) / dtice
!     &             * rside(i,j) / dt

            fhnet(i,j)      = fhnet(i,j)      + dfhnet
            fhnet_hist(i,j) = fhnet_hist(i,j) + dfhnet 

            fresh(i,j)      = fresh(i,j)      + dfresh
            fresh_hist(i,j) = fresh_hist(i,j) + dfresh

            fsalt(i,j)      = fsalt(i,j)      + dfsalt
            fsalt_hist(i,j) = fsalt_hist(i,j) + dfsalt

            ! history diagnostics
            meltl(i,j) = meltl(i,j) + vicen(i,j,ni)*rside(i,j)

            ! state variables (except ice energy)
            aicen(i,j,ni) = aicen(i,j,ni) * (c1i - rside(i,j))
            vicen(i,j,ni) = vicen(i,j,ni) * (c1i - rside(i,j))
            vsnon(i,j,ni) = vsnon(i,j,ni) * (c1i - rside(i,j))
            esnon(i,j,ni) = esnon(i,j,ni) * (c1i - rside(i,j))

         enddo                  ! ij

         do k = 1, nilyr
            do ij = 1, icells
               i = indxi(ij)
               j = indxj(ij)
               
               ! heat flux to coupler for ice melt (dfhnet < 0)

!               dfhnet = eicen(i,j,ilyr1(ni)+k-1)*rside(i,j) / dt
               dfhnet = eicen(i,j,ilyr1(ni)+k-1)*rside(i,j) / dtice
               fhnet(i,j)      = fhnet(i,j)      + dfhnet 
               fhnet_hist(i,j) = fhnet_hist(i,j) + dfhnet 

               ! ice energy
               eicen(i,j,ilyr1(ni)+k-1) = eicen(i,j,ilyr1(ni)+k-1)  &
     &                                 * (c1i - rside(i,j))
            enddo               ! ij
         enddo                  ! k

      enddo  ! n

      end subroutine lateral_melt
!=======================================================================
!BOP
!
! !ROUTINE: freeboard - snow-ice conversion
!
! !DESCRIPTION:
!
! If there is enough snow to lower the ice/snow interface below 
! sea level, convert enough snow to ice to bring the interface back 
! to sea level.
!
! NOTE: Subroutine rebin should be called after freeboard to make sure
!       ice thicknesses are within category bounds.
!
! !REVISION HISTORY:
!
! authors William H. Lipscomb, LANL
!         Elizabeth C. Hunke, LANL
!
! !INTERFACE:
!
      subroutine freeboard
!
! !USES:
! 
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) :: i, j, ni, k

      real (kind=dbl_kind) ::   &
     &   hi       & ! ice thickness (m)
     &,  hs       & ! snow thickness (m)
     &,  dhi      & ! change in ice thickness (m)
     &,  dhs      & ! change in snow thickness (m)
     &,  dz       & ! distance freeboard below SL (m)
     &,  fs         ! salt flux due to snow-ice conversion

      real (kind=dbl_kind), dimension(ilo:ihi,jlo:jhi) ::  &
     &   de        ! energy transferred to each ice layer(J/m^2)

      integer (kind=int_kind) ::      &
     &   icells   & ! number of cells with ice 
     &,  jcells   & ! number of cells with freeboard adjustment
     &,  ij        ! combined i/j horizontal index

      integer (kind=int_kind),      &
     &         dimension (1:(ihi-ilo+1)*(jhi-jlo+1)) :: &
     &  indxi     & ! compressed i/j indices
     &, indxj     &
     &, indxi2    & 
     &, indxj2

      do ni = 1, ncat

      !-----------------------------------------------------------------
      ! Identify grid cells with ice.
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

         jcells = 0             ! freeboard adjustment counter

         do ij = 1, icells      ! aicen > puny
            i = indxi(ij)
            j = indxj(ij)

      !-----------------------------------------------------------------
      ! Determine whether snow lies below freeboard.
      !-----------------------------------------------------------------

            hi = vicen(i,j,ni) / aicen(i,j,ni)
            hs = vsnon(i,j,ni) / aicen(i,j,ni)

            dz = hs - hi*(rhow-rhoi)/rhos

            if (dz > puny .and. hs > puny) then ! snow below freeboard
               jcells = jcells + 1
               indxi2(jcells) = i
               indxj2(jcells) = j

               dhs = min(dz*rhoi/rhow, hs) ! snow to remove
               dhi = dhs*rhos/rhoi         ! ice to add

      !-----------------------------------------------------------------
      ! Compute energy transferred from snow to ice.
      ! NOTE: It would be more realistic to transfer energy only to
      !       the top ice layer, but it is simpler to transfer equal 
      !       energy to all layers.)
      !-----------------------------------------------------------------

               de(i,j) = esnon(i,j,ni)*dhs/hs
               esnon(i,j,ni) = esnon(i,j,ni) - de(i,j)
               de(i,j) = de(i,j)/real(nilyr,kind=dbl_kind) ! energy to each ice layer

      !-----------------------------------------------------------------
      ! Adjust snow and ice volume.
      !-----------------------------------------------------------------

               hi = hi + dhi
               hs = hs - dhs
               vicen(i,j,ni) = hi * aicen(i,j,ni)
               vsnon(i,j,ni) = hs * aicen(i,j,ni)

      !-----------------------------------------------------------------
      ! Update history and coupler variables.
      !-----------------------------------------------------------------

            ! history diagnostic
               snoice(i,j) = snoice(i,j) + dhi*aicen(i,j,ni) 
                  
            ! Remove salt from the ocean.
            ! This is not physically realistic but is needed to 
            ! conserve salt, because salt will be returned to the ocean 
            ! when the ice melts. 

!               fs = -ice_ref_salinity*p001*aicen(i,j,n)*dhi*rhoi/dt 
               fs = -ice_ref_salinity*p001*aicen(i,j,ni)*dhi*rhoi/dtice
               fsalt(i,j)      = fsalt(i,j)      + fs                 
               fsalt_hist(i,j) = fsalt_hist(i,j) + fs

            endif               ! dz > puny and hs > puny
         enddo                  ! ij

      !-----------------------------------------------------------------
      ! Adjust ice energy.
      !-----------------------------------------------------------------
         do k = 1, nilyr
            do ij = 1, jcells   ! just cells with freeboard adjustment
               i = indxi2(ij)
               j = indxj2(ij)
               eicen(i,j,ilyr1(ni)+k-1) =                &
     &              eicen(i,j,ilyr1(ni)+k-1) + de(i,j)
            enddo               ! ij
         enddo                  ! k

      enddo                     ! n

      end subroutine freeboard

!=======================================================================

      end module ice_therm_itd

!=======================================================================
