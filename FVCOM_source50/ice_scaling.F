!/===========================================================================/
! CVS VERSION INFORMATION
! $Id$
! $Name$
! $Revision$
!/===========================================================================/

!=======================================================================
!BOP
!
! !MODULE: ice_scaling
!
! !DESCRIPTION:
!
!  Scale ice fluxes by ice area
!
! !REVISION HISTORY:
!
! author: C.M.Bitz ?
! Modified: 16 Oct 03 For CICE-CSIM merge and vectorization
!
! !INTERFACE:
!
      module ice_scaling
!
! !USES:
!
      use ice_domain
      use ice_kinds_mod
      use ice_constants
      use ice_state
      use ice_flux
      use ice_grid, only: tmask
!
!EOP
!
      implicit none

!=======================================================================

      contains

!=======================================================================
!
!BOP
!
! !IROUTINE: scale_fluxes
!
! !DESCRIPTION:
!
!  Divide ice fluxes by ice area before sending them to the
!  coupler, since the coupler multiplies by ice area. This
!  is the ice area at the beginning of the timestep, i.e.
!  the value sent to the coupler.
!
! !INTERFACE:
!
      subroutine scale_fluxes
!
! !REVISION HISTORY:
!
! author: C.M.Bitz ?
!
! !USES:
!
      use ice_albedo
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      real (kind=dbl_kind) :: ar   ! 1/aice

      integer (kind=int_kind) :: i, j

      do j=jlo,jhi
       do i=ilo,ihi
          if (tmask(i,j) .and. aice(i,j) > c0i) then
             ar = c1i/aice(i,j)
!             strairxT(i,j) = strairxT(i,j) * ar
!             strairyT(i,j) = strairyT(i,j) * ar
!     the dynamic is not considered

             fsens   (i,j) = fsens   (i,j) * ar
             flat    (i,j) = flat    (i,j) * ar
             fswabs  (i,j) = fswabs  (i,j) * ar
             flwout  (i,j) = flwout  (i,j) * ar
             evap    (i,j) = evap    (i,j) * ar
             Tref    (i,j) = Tref    (i,j) * ar  
             Qref    (i,j) = Qref    (i,j) * ar  
             fresh   (i,j) = fresh   (i,j) * ar
             fsalt   (i,j) = fsalt   (i,j) * ar
             fhnet   (i,j) = fhnet   (i,j) * ar
             fswthru (i,j) = fswthru (i,j) * ar
          else  ! zero out fluxes 
             strairxT(i,j) = c0i
             strairyT(i,j) = c0i
!    dynamic part using fvcom code
             fsens   (i,j) = c0i
             flat    (i,j) = c0i
             fswabs  (i,j) = c0i
             flwout  (i,j) = -stefan_boltzmann*(Tffresh - Tf(i,j))**4 
                                ! to make upward longwave over ocean
                                ! reasonable for history file
             evap    (i,j) = c0i
             Tref    (i,j) = Tair(i,j)
             Qref    (i,j) = Qa(i,j)
             fresh   (i,j) = c0i
             fsalt   (i,j) = c0i
             fhnet   (i,j) = c0i
             fswthru (i,j) = c0i    
             alvdf   (i,j) = c0i ! zero out albedo where ice is absent
             alidf   (i,j) = c0i
             alvdr   (i,j) = c0i
             alidr   (i,j) = c0i
          endif 
       enddo  
      enddo

      end subroutine scale_fluxes

!=======================================================================
!
!BOP
!
! !IROUTINE: scale_hist_fluxes - scale history fluxes
!
! !INTERFACE:
!
      subroutine scale_hist_fluxes
!
! !DESCRIPTION:
!
!  Divide ice fluxes by ice area used by the coupler before writing out
!  diagnostics. aice\_init is the ice area saved from coupling.  This
!  makes the fluxes written to the history file consistent with those
!  sent to the coupler.
!
! !REVISION HISTORY:
!
! author: C.M.Bitz ?
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      real (kind=dbl_kind) :: ar  ! 1/aice

      integer (kind=int_kind) :: i, j

      do j=jlo,jhi
       do i=ilo,ihi
          if (tmask(i,j) .and. aice_init(i,j) > c0i) then
             ar = c1i/aice_init(i,j)
             fresh_hist  (i,j)  = fresh_hist  (i,j) * ar
             fsalt_hist  (i,j)  = fsalt_hist  (i,j) * ar
             fhnet_hist  (i,j)  = fhnet_hist  (i,j) * ar
             fswthru_hist(i,j)  = fswthru_hist(i,j) * ar
          else
             fresh_hist  (i,j)  = c0i
             fsalt_hist  (i,j)  = c0i
             fhnet_hist  (i,j)  = c0i
             fswthru_hist(i,j)  = c0i
          endif 
       enddo  
      enddo

      end subroutine scale_hist_fluxes

!=======================================================================

      end module ice_scaling

!=======================================================================
