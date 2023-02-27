!/===========================================================================/
! CVS VERSION INFORMATION
! $Id$
! $Name$
! $Revision$
!/===========================================================================/

!=======================================================================
!BOP
!
! !MODULE: ice_ocean - ocean mixed layer internal to sea ice model
!
! !DESCRIPTION:
!
! Ocean mixed layer calculation (internal to sea ice model).
! Allows heat storage in ocean for uncoupled runs.
!
! !REVISION HISTORY:
!
! authors:   John Weatherly, CRREL
!            C.M. Bitz, UW
!            Elizabeth C. Hunke, LANL
!            Bruce P. Briegleb, NCAR
!            William H. Lipscomb, LANL
!
! !INTERFACE:
!
      module ice_ocean
!
! !USES:
!
      use ice_kinds_mod
      use ice_constants
!      use ice_calendar, only: dt
      use ice_calendar, only: dtice
!
!EOP
!
      implicit none
      save

      logical (kind=log_kind) :: &
        oceanmixed_ice           ! if true, use ocean mixed layer

!=======================================================================

      contains

!=======================================================================
!BOP
!
! !ROUTINE: mixed_layer - compute SST and freeze/melt potential
!
! !DESCRIPTION:
!
! Compute the mixed layer heat balance and update the SST.
! Compute the energy available to freeze or melt ice.
! NOTE: SST changes due to fluxes through the ice are computed in 
!       ice_therm_vertical.
!
! !REVISION HISTORY:
!
! authors:   John Weatherly, CRREL
!            C.M. Bitz, UW
!            Elizabeth C. Hunke, LANL
!            Bruce P. Briegleb, NCAR
!            William H. Lipscomb, LANL
!
! !INTERFACE:
!
      subroutine mixed_layer
!
! !USES:
! 
      use ice_flux
      use ice_grid, only: tmask
      use ice_atmo
      use ice_state
      use ice_albedo
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) :: &
         i, j            ! horizontal indices

      real (kind=dbl_kind), dimension(ilo:ihi,jlo:jhi) ::&
         delt    & ! potential temperature difference   (K)
      ,  delq    & ! specific humidity difference   (kg/kg)
      ,  dummy1, dummy2, dummy3, dummy4  ! dummy arrays

      real (kind=dbl_kind) :: &
         TsfK     &! surface temperature (K)
      ,  fsh_t      &! sensible heat flux  (W/m^2)
      ,  flh_t      &! latent heat flux    (W/m^2)
      ,  swabs    &! surface absorbed shortwave heat flux (W/m^2)
      ,  flwup_t    &! long-wave upward heat flux  (W/m^2)
      ,  ft       &! fraction reduction of positive qdp
      ,  dtcprho  ! dt/cp_ocn*rhow

      dtcprho = dtice/(cp_ocn*rhow)

      call atmo_boundary_layer (1, 'ocn', sst, &
           dummy1, dummy2, dummy3, dummy4, delt,  delq)

      do j = jlo,jhi
      do i = ilo,ihi
       if (tmask(i,j)) then
!       if (hmix(i,j) > puny) then

         ! ocean surface temperature in Kelvin
         TsfK = sst(i,j) + Tffresh

         ! shortwave radiative flux
         swabs = (c1i - albocn) * fsw(i,j)

         ! longwave radiative flux
         flwup_t  = -stefan_boltzmann * TsfK**4

         ! downward latent and sensible heat fluxes
         flh_t = lhcoef(i,j) * delq(i,j)
         fsh_t = shcoef(i,j) * delt(i,j)

         ! compute sst change due to exchange with atm/ice above
         ! Note: fhnet, fswthru are added in ice_therm_vertical.F
           sst(i,j) = sst(i,j) + (c1i-aice(i,j))*dtcprho/hmix(i,j) &
                         *(fsh_t + flh_t + flwup_t + flw(i,j) + swabs)

         ! adjust qdp if cooling of mixed layer would occur when sst le Tf
         if( sst(i,j) <= Tf(i,j) .and. qdp(i,j) > c0i ) qdp(i,j) = c0i

         ! computed T change due to exchange with deep layers:
           sst(i,j) = sst(i,j) - qdp(i,j)*dtcprho/hmix(i,j)

         ! compute potential to freeze or melt ice
         frzmlt(i,j) = (Tf(i,j)-sst(i,j))/dtcprho*hmix(i,j)
         frzmlt(i,j) = min(max(frzmlt(i,j),-c1000),c1000)

         ! if sst is below freezing, reset sst to Tf
         if (sst(i,j) <= Tf(i,j)) sst(i,j) = Tf(i,j)

       else

         frzmlt(i,j) = c0i
         sst(i,j) = Tf(i,j)

!       endif                    ! hmix > puny
       endif                    ! tmask
      enddo                     ! i
      enddo                     ! j

      end subroutine mixed_layer

!=======================================================================

      end module ice_ocean

!=======================================================================
