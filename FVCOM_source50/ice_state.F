!/===========================================================================/
! CVS VERSION INFORMATION
! $Id$
! $Name$
! $Revision$
!/===========================================================================/

!=======================================================================
!BOP
!
! !MODULE: ice_state - primary state variables
!
! !DESCRIPTION:
!
! Primary state variables in various configurations
! Note: other state variables are at the end of this...
! The primary state variable names are: 
!-------------------------------------------------------------------
! for each category   aggregated over     units
!                       categories      
!-------------------------------------------------------------------
! aicen(i,j,n)         aice(i,j)           ---
! vicen(i,j,n)         vice(i,j)           m
! vsnon(i,j,n)         vsno(i,j)           m
! eicen(i,j,k)         eice(i,j)           J/m$^2$
! esnon(i,j,n)         esno(i,j)           J/$m^2$
! Tsfcn(i,j,n)         Tsfc(i,j)           deg
!
! Area is dimensionless because aice is the fractional area 
! (normalized so that the sum over all categories, including open 
! water, is 1.0).  That is why vice/vsno have units of m instead of 
! m$^3$, and eice/esno have units of J/m$^2$ instead of J.
!
! Variable names follow these rules:
!
! (1) For 3D variables (indices i,j,n), write 'ice' or 'sno' or 
!     'sfc' and put an 'n' at the end.
! (2) For 2D variables (indices i,j) aggregated over all categories, 
!     write 'ice' or 'sno' or 'sfc' without the 'n'.
! (3) For 2D variables (indices i,j) associated with an individual 
!     category, write 'i' or 's' instead of 'ice' or 'sno' and put an 'n'
!     at the end: e.g. hin, hsn.  These are not declared here
!     but in individual modules (e.g., ice\_therm\_vertical). 
!
! !REVISION HISTORY:
!
! authors C. M. Bitz, UW
!         Elizabeth C. Hunke, LANL
!         William H. Lipscomb, LANL
!
! !INTERFACE:
!
      module ice_state
!
! !USES:
!
      use ice_kinds_mod
      use ice_model_size
      use ice_domain
!
!EOP
!
      implicit none

      !-----------------------------------------------------------------
      ! state of the ice aggregated over all categories
      !-----------------------------------------------------------------

!      real (kind=dbl_kind), dimension (imt_local,jmt_local) :: &
!      real (kind=dbl_kind),dimension(:,:),allocatable,save :: &
      real (kind=dbl_kind),dimension(:,:),allocatable,save,TARGET :: &
         aice     &! concentration of ice 
      ,  vice     &! volume per unit area of ice          (m)
      ,  vsno     &! volume per unit area of snow         (m)
      ,  Tsfc     &! temperature of ice/snow top surface  (C)
      ,  eice     &! energy of melt. of ice           (J/m^2)
      ,  esno     &! energy of melt. of snow layer    (J/m^2)
      ,  tmass    &! total mass of ice and snow
      ,  aice_init ! concentration of ice at beginning of dt (for diagnostics)
! - begin afm 20171113 for ice open boundary -------------
 real (kind=dbl_kind),dimension(:),allocatable,save,TARGET ::  xflux_ice
! - end afm 20171113 for ice open boundary -------------

      !-----------------------------------------------------------------
      ! state of the ice for each category
      !-----------------------------------------------------------------

!      real (kind=dbl_kind), dimension (imt_local,jmt_local,ncat) :: &
      real (kind=dbl_kind),dimension(:,:,:),allocatable,save,TARGET :: &
         aicen    &! concentration of ice 
      ,  vicen    &! volume per unit area of ice          (m)
      ,  vsnon    &! volume per unit area of snow         (m)
      ,  Tsfcn    &! temperature of ice/snow top surface  (C)
      ,  esnon     ! energy of melt. of snow layer    (J/m^2)

!      real (kind=dbl_kind) ::                                      &
      real (kind=dbl_kind),dimension(:,:),allocatable,save,TARGET :: &
         aice0       ! concentration of open water

      real (kind=dbl_kind),dimension(:,:,:),allocatable,save,TARGET :: &
        eicen        ! energy of melting for
      
!         aice0(imt_local,jmt_local)       &! concentration of open water
!      ,  eicen(imt_local,jmt_local,ntilay) ! energy of melting for
                                           ! each ice layer (J/m^2)

      !-----------------------------------------------------------------
      ! other variables closely related to the state of the ice
      !-----------------------------------------------------------------

!      real (kind=dbl_kind), dimension (imt_local,jmt_local) :: &
      real (kind=dbl_kind),dimension(:,:),allocatable,save :: &
         uvel        &! x-component of velocity (m/s)
      ,  vvel        &! y-component of velocity (m/s)
      ,  strength     ! ice strength (N/m)



      end module ice_state

!=======================================================================
