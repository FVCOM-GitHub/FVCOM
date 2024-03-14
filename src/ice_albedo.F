!/===========================================================================/
! CVS VERSION INFORMATION
! $Id$
! $Name$
! $Revision$
!/===========================================================================/

!=======================================================================
!BOP
!
! !MODULE: ice_albedo - snow and ice albedo parameterization
!
! !DESCRIPTION:
!
! The albedo parameterization
!
! !REVISION HISTORY:
!
! authors:  Bruce P. Briegleb, NCAR 
!           Elizabeth C. Hunke, LANL
!
! Vectorized by Clifford Chen (Fujitsu) and William H. Lipscomb (LANL)
!
! !INTERFACE:
!
      module ice_albedo
!
! !USES:
!
      use ice_kinds_mod
      use ice_domain
!
!EOP
!
      implicit none

      real (kind=dbl_kind), parameter  ::  &
         albocn = 0.06_dbl_kind  ! ocean albedo

      ! weights for albedos match those for isccp shortwave forcing
      real (kind=dbl_kind), parameter  ::  &             ! currently used only
         awtvdr = 0.29_dbl_kind  &! visible, direct   ! for history and
      ,  awtidr = 0.31_dbl_kind  &! near IR, direct   ! diagnostics
      ,  awtvdf = 0.24_dbl_kind  &! visible, diffuse
      ,  awtidf = 0.16_dbl_kind   ! near IR, diffuse

      ! parameter for fractional snow area 
      real (kind=dbl_kind), parameter  ::  &
         snowpatch = 0.02_dbl_kind

!      ! albedos for ice in each category
!      real (kind=dbl_kind)  ::  &
!         alvdrn (ilo:ihi,jlo:jhi,ncat) &! visible, direct   (fraction)
!      ,  alidrn (ilo:ihi,jlo:jhi,ncat) &! near-ir, direct   (fraction)
!      ,  alvdfn (ilo:ihi,jlo:jhi,ncat) &! visible, diffuse  (fraction)
!      ,  alidfn (ilo:ihi,jlo:jhi,ncat)  ! near-ir, diffuse  (fraction)

      ! albedos aggregated over categories
!      real (kind=dbl_kind)  ::  &
!         alvdr (ilo:ihi,jlo:jhi)  &! visible, direct   (fraction)
!      ,  alidr (ilo:ihi,jlo:jhi)  &! near-ir, direct   (fraction)
!      ,  alvdf (ilo:ihi,jlo:jhi)  &! visible, diffuse  (fraction)
!      ,  alidf (ilo:ihi,jlo:jhi)   ! near-ir, diffuse  (fraction)

     real (kind=dbl_kind),dimension(:,:,:),allocatable,save  ::  &
         alvdrn &! visible, direct   (fraction)
      ,  alidrn &! near-ir, direct   (fraction)
      ,  alvdfn &! visible, diffuse  (fraction)
      ,  alidfn   ! near-ir, diffuse  (fraction)

      ! albedos aggregated over categories
      real (kind=dbl_kind),dimension(:,:),allocatable,save  ::  &
         alvdr   &! visible, direct   (fraction)
      ,  alidr   &! near-ir, direct   (fraction)
      ,  alvdf   &! visible, diffuse  (fraction)
      ,  alidf    ! near-ir, diffuse  (fraction)


      ! baseline albedos for thick cases
      real (kind=dbl_kind)  ::  &
         albicev   &! visible ice albedo for h > ahmax
      ,  albicei   &! near-ir ice albedo for h > ahmax
      ,  albsnowv  &! cold snow albedo, visible
      ,  albsnowi   ! cold snow albedo, near IR

!=======================================================================
 
      contains

!=======================================================================
!BOP
!
! !IROUTINE: albedos - compute snow/ice albedos and aggregate
!
! !INTERFACE:
!
      subroutine albedos
!
! !DESCRIPTION:
!
! Compute albedos and aggregate them \\
! note: ice albedo is zero if no ice present
!
! !REVISION HISTORY:
!
! authors:  Bruce P. Briegleb, NCAR 
!           Elizabeth C. Hunke, LANL
!
! !USES:
!
      use ice_constants
      use ice_grid
      use ice_state
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      real (kind=dbl_kind), parameter  ::  & 
         ahmax     = 0.5_dbl_kind    &! thickness above which ice albedo constant (m)
      ,  dT_mlt    = 1._dbl_kind     &! change in temp to give dalb_mlt albedo change
      ,  dalb_mlt  = -0.075_dbl_kind &! albedo change per dT_mlt change 
                                      ! in temp for ice
      ,  dalb_mltv = -0.100_dbl_kind &! albedo vis change per dT_mlt change 
                                      ! in temp for snow
      ,  dalb_mlti = -0.150_dbl_kind  ! albedo nir change per dT_mlt change 
                                     ! in temp for snow

      integer (kind=int_kind) :: i, j, ni

      real (kind=dbl_kind)  ::  & 
         hi      &! ice thickness  (m)
      ,  hs      &! snow thickness (m)
      ,  albo    &! effective ocean albedo, function of ice thickness
      ,  asnow   &! snow-covered area fraction
      ,  asnwv   &! snow albedo, visible 
      ,  asnwi   &! snow albedo, near IR
      ,  fh      &! piecewise linear function of thickness 
      ,  fT      &! piecewise linear function of surface temperature
      ,  dTs     &! difference of Tsfc and Timelt
      ,  fhtan    ! factor used in albedo dependence on ice thickness

      integer (kind=int_kind)  ::  &
         icells  &! number of ice/ocean grid cells
      ,  ij       ! horizontal index, combines i and j loops

      integer (kind=int_kind), dimension (1:(ihi-ilo+1)*(jhi-jlo+1))  ::  &
         indxi   &! compressed indices for ice/ocean cells
      ,  indxj

      fhtan = atan(ahmax*c4i)

      icells = 0
      do j = jlo, jhi
      do i = ilo, ihi
         if (tmask(i,j)) then
            icells = icells + 1
            indxi(icells) = i
            indxj(icells) = j
         endif                  ! tmask
      enddo
      enddo

      !-----------------------------------------------------------------
      ! albedo for each thickness category
      !-----------------------------------------------------------------

      do ni = 1, ncat
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)

            if (aicen(i,j,ni) > puny) then
               hi = vicen(i,j,ni) / aicen(i,j,ni)
               hs = vsnon(i,j,ni) / aicen(i,j,ni)

               ! bare ice, thickness dependence
               fh = min(atan(hi*c4i)/fhtan,c1i)
               albo = albocn*(c1i-fh)
               alvdfn(i,j,ni) = albicev*fh + albo
               alidfn(i,j,ni) = albicei*fh + albo

               ! bare ice, temperature dependence
               dTs = Timelt - Tsfcn(i,j,ni)
               fT = min(dTs/dT_mlt-c1i,c0i)
               alvdfn(i,j,ni) = alvdfn(i,j,ni) - dalb_mlt*fT
               alidfn(i,j,ni) = alidfn(i,j,ni) - dalb_mlt*fT

               ! avoid negative albedos for thin, bare, melting ice
               alvdfn(i,j,ni) = max (alvdfn(i,j,ni), albocn)
               alidfn(i,j,ni) = max (alidfn(i,j,ni), albocn)

               if( hs > puny ) then

                  ! fractional area of snow on ice (thickness dependent)
                  asnow = hs / ( hs + snowpatch ) 
                  asnwv = albsnowv
                  asnwi = albsnowi

                  ! snow on ice, temperature dependence
                  asnwv = asnwv - dalb_mltv*fT
                  asnwi = asnwi - dalb_mlti*fT

                  ! combine ice and snow albedos
                  alvdfn(i,j,ni) = alvdfn(i,j,ni)*(c1i-asnow) + &
                                  asnwv*asnow
                  alidfn(i,j,ni) = alidfn(i,j,ni)*(c1i-asnow) + &
                                  asnwi*asnow
               endif            ! hs > puny

               alvdrn(i,j,ni) = alvdfn(i,j,ni)
               alidrn(i,j,ni) = alidfn(i,j,ni)

            else                ! no ice
               alvdfn(i,j,ni) = albocn
               alidfn(i,j,ni) = albocn
               alvdrn(i,j,ni) = albocn
               alidrn(i,j,ni) = albocn
            endif               ! aicen > puny
         enddo                  ! ij
      enddo                     ! ncat

      !-----------------------------------------------------------------
      ! aggregate
      !-----------------------------------------------------------------

      do j = jlo, jhi
      do i = ilo, ihi
         alvdf(i,j) = c0i
         alidf(i,j) = c0i
         alvdr(i,j) = c0i
         alidr(i,j) = c0i
      enddo
      enddo

      do ni = 1, ncat
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)
            if (aice(i,j) > puny) then
               alvdf(i,j) = alvdf(i,j) + alvdfn(i,j,ni)*aicen(i,j,ni)
               alidf(i,j) = alidf(i,j) + alidfn(i,j,ni)*aicen(i,j,ni)
               alvdr(i,j) = alvdr(i,j) + alvdrn(i,j,ni)*aicen(i,j,ni)
               alidr(i,j) = alidr(i,j) + alidrn(i,j,ni)*aicen(i,j,ni)
            endif
         enddo                  ! ij
      enddo                     ! ncat

      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)
         if (aice(i,j) > puny) then
            alvdf(i,j) = alvdf(i,j) / aice(i,j)
            alidf(i,j) = alidf(i,j) / aice(i,j)
            alvdr(i,j) = alvdr(i,j) / aice(i,j)
            alidr(i,j) = alidr(i,j) / aice(i,j)
         endif                  ! aicen > puny
      enddo                     ! ij

      end subroutine albedos

!=======================================================================

      end module ice_albedo

!=======================================================================
