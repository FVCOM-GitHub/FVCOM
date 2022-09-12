!/===========================================================================/
! CVS VERSION INFORMATION
! $Id$
! $Name$
! $Revision$
!/===========================================================================/

!=======================================================================
!BOP
!
! !MODULE: ice_constants - sets physical constants
!
! !DESCRIPTION:
!
! This module defines a variety of physical and numerical constants
! used throughout the ice model \\
!
! Code originally based on constants.F in POP
!
! !REVISION HISTORY:
!
! author Elizabeth C. Hunke, LANL
!
! !INTERFACE:

      module ice_constants
      use control, only : pi,pi2     

!
! !USES:
!
!#ifdef CCSMcoupled
!      use shr_const_mod
!#endif
      use ice_kinds_mod
      use ice_domain
!
!EOP
!
      implicit none
      save

      !-----------------------------------------------------------------
      ! physical constants
      !-----------------------------------------------------------------

!#ifdef CCSMcoupled
!      real (kind=dbl_kind), parameter :: &
!     &   pi        = SHR_CONST_PI     ! pi
!     &,  gravit    = SHR_CONST_G      ! gravitational acceleration (m/s^2)
!     &,  secday    = SHR_CONST_CDAY   ! seconds in calendar day
!     &,  omega     = SHR_CONST_OMEGA  ! angular velocity of earth (rad/sec)
!     &,  radius    = SHR_CONST_REARTH ! earth radius (m)
!     &,  rhoi      = SHR_CONST_RHOICE ! density of ice (kg/m^3)
!     &,  rhow      = SHR_CONST_RHOSW  ! density of seawater (kg/m^3)
!     &,  rhofresh  = SHR_CONST_RHOFW  ! density of fresh water (kg/m^3)
!     &,  zvir      = SHR_CONST_ZVIR   ! rh2o/rair - 1.0
!     &,  vonkar    = SHR_CONST_KARMAN ! von Karman constant
!     &,  cp_air    = SHR_CONST_CPDAIR ! specific heat of air (J/kg/K)
!     &,  cp_wv     = SHR_CONST_CPWV   ! specific heat of water vapor (J/kg/K)
!     &,  stefan_boltzmann = SHR_CONST_STEBOL !  W/m^2/K^4
!     &,  Tffresh   = SHR_CONST_TKFRZ  ! freezing temp of fresh ice (K)
!     &,  cp_ice    = SHR_CONST_CPICE  ! specific heat of fresh ice (J/kg/K)
!     &,  cp_ocn    = SHR_CONST_CPSW   ! specific heat of ocn    (J/kg/K)
!     &,  Lsub      = SHR_CONST_LATSUB ! latent heat, sublimation freshwater (J/kg)
!     &,  Lvap      = SHR_CONST_LATVAP ! latent heat, vaporization freshwater (J/kg)
!     &,  Lfresh    = SHR_CONST_LATICE ! latent heat of melting of fresh ice (J/kg)
!     &,  Timelt    = SHR_CONST_TKFRZ-SHR_CONST_TKFRZ ! melting temp. ice top surface  (C)
!     &,  Tsmelt    = SHR_CONST_TKFRZ-SHR_CONST_TKFRZ ! melting temp. snow top surface (C)
!     &,  ice_ref_salinity = SHR_CONST_ICE_REF_SAL  ! (psu)
!     &,  spval     = SHR_CONST_SPVAL  ! special value
!c    &,  ocn_ref_salinity = SHR_CONST_OCN_REF_SAL  ! (psu)
!c    &,  rho_air   = SHR_CONST_RHODAIR  ! ambient air density (kg/m^3)
!#else
      real (kind=dbl_kind), parameter :: &
         gravit    = 9.80616_dbl_kind  &! gravitational acceleration (m/s^2)
      ,  secday    = 86400.0_dbl_kind  &! seconds in calendar day
      ,  omega     = 7.292e-5_dbl_kind &! angular velocity of earth (rad/sec)
      ,  radius    = 6.37e6_dbl_kind   &! earth radius (m)
      ,  rhoi      = 917.0_dbl_kind    &! density of ice (kg/m^3)
      ,  rhow      = 1026.0_dbl_kind   &! density of seawater (kg/m^3)
      ,  rhofresh  = 1000.0_dbl_kind   &! density of fresh water (kg/m^3)
      ,  zvir      = 0.606_dbl_kind    &! rh2o/rair - 1.0
      ,  vonkar    = 0.4_dbl_kind      &! von Karman constant
      ,  cp_air    = 1005.0_dbl_kind   &! specific heat of air (J/kg/K)
      ,  cp_wv     = 1.81e3_dbl_kind   &! specific heat of water vapor (J/kg/K)
      ,  stefan_boltzmann = 567.0e-10_dbl_kind &!  W/m^2/K^4
      ,  Tffresh   = 273.15_dbl_kind   &! freezing temp of fresh ice (K)
!#ifdef CCSM
!     &,  cp_ice = 2054._dbl_kind       ! specific heat of fresh ice (J/kg/K)
!#else
      ,  cp_ice = 2106._dbl_kind      & ! specific heat of fresh ice (J/kg/K)
!#endif
      ,  cp_ocn    = 4218._dbl_kind    &! specific heat of ocn    (J/kg/K)
      ,  Lsub      = 2.835e6_dbl_kind  &! latent heat, sublimation freshwater (J/kg)
      ,  Lvap      = 2.501e6_dbl_kind  &! latent heat, vaporization freshwater (J/kg)
      ,  Lfresh    = Lsub-Lvap         &! latent heat of melting of fresh ice (J/kg)
      ,  Timelt    = 0.0_dbl_kind      &! melting temperature, ice top surface  (C)
      ,  Tsmelt    = 0.0_dbl_kind      &! melting temperature, snow top surface (C)
#  if defined (ICE_FRESHWATER)
! afm 20151112 & EJA 20160920
      ,  ice_ref_salinity = 0._dbl_kind   ! (psu)
#  else
      ,  ice_ref_salinity = 4._dbl_kind   ! (psu)
#  endif
!c    &,  ocn_ref_salinity   = 34.7_dbl_kind  ! (psu)
!c    &,  rho_air   = 1.2_dbl_kind      ! ambient air density (kg/m^3)

      real (kind=real_kind), parameter :: &
     &   spval     = 1.0e30_real_kind  ! special value
!#endif

      real (kind=dbl_kind), parameter :: &
         rhos     = 330.0_dbl_kind     &! density of snow (kg/m^3)
      ,  depressT = 0.054_dbl_kind     &! Tf:brine salinity ratio (C/ppt)
      ,  iceruf   = 0.0005_dbl_kind    & ! ice surface roughness (m)

         ! (Ebert, Schramm and Curry JGR 100 15965-15975 Aug 1995)
      ,  kappav = 1.4_dbl_kind   &! vis extnctn coef in ice, wvlngth<700nm (1/m)
      ,  kappan = 17.6_dbl_kind  &! vis extnctn coef in ice, wvlngth<700nm (1/m)

         ! (Briegleb JGR 97 11475-11485  July 1992)
      ,  emissivity = 0.95_dbl_kind  & ! emissivity of snow and ice

      ,  kice   = 2.03_dbl_kind   &! thermal conductivity of fresh ice(W/m/deg)
      ,  ksno   = 0.30_dbl_kind   &! thermal conductivity of snow  (W/m/deg)
      ,  zref   = 10._dbl_kind     ! reference height for stability (m)

      real (kind=dbl_kind), parameter :: &
        qqqice       = 11637800._dbl_kind       &! for qsat over ice
      , TTTice       = 5897.8_dbl_kind          &! for qsat over ice
      , qqqocn       = 627572.4_dbl_kind        &! for qsat over ocn
      , TTTocn       = 5107.4_dbl_kind          ! for qsat over ocn

      !-----------------------------------------------------------------
      ! numbers
      !-----------------------------------------------------------------

!      real (kind=dbl_kind), dimension(ilo:ihi,jlo:jhi) :: &
!        one                 ! array of ones which is often useful
!     defined in other place ??



      real (kind=dbl_kind), parameter :: &
        c0i   = 0.0_dbl_kind,             &
        c1i   = 1.0_dbl_kind,             &
        c2i   = 2.0_dbl_kind,             &
        c3i   = 3.0_dbl_kind,             &
        c4i   = 4.0_dbl_kind,             &
        c5i   = 5.0_dbl_kind,             &
        c6i   = 6.0_dbl_kind,             &
        c7i   = 7.0_dbl_kind,             &
        c8i   = 8.0_dbl_kind,             &
        c9i   = 9.0_dbl_kind,             &
        c10i  = 10.0_dbl_kind,            &
        c12  = 12.0_dbl_kind,             &
        c15  = 15.0_dbl_kind,             &
        c16  = 16.0_dbl_kind,             &
        c20  = 20.0_dbl_kind,             &
        c25  = 25.0_dbl_kind,             &
        c90  = 90.0_dbl_kind,             &
        c100 = 100.0_dbl_kind,            &
        c360 = 360.0_dbl_kind,            &
        c365 = 365.0_dbl_kind,            &
        c1000= 1000.0_dbl_kind,           &
        c1p5 = 1.5_dbl_kind,              &
        p001 = 0.001_dbl_kind,            &
        p01  = 0.01_dbl_kind,             &
        p1   = 0.1_dbl_kind,              &
        p2   = 0.2_dbl_kind,              &
        p4   = 0.4_dbl_kind,              &
        p5   = 0.5_dbl_kind,              &
        p6   = 0.6_dbl_kind,              &
        p15  = 0.15_dbl_kind,             &
        p25  = 0.25_dbl_kind,             &
        p166 = c1i/c6i,                     &
        p333 = c1i/c3i,                     &
        p666 = c2i/c3i,                     &
        p111 = c1i/c9i,                     &
        p055 = p111*p5,                   &
        p027 = p055*p5,                   &
        p222 = c2i/c9i,                     &
        eps04  = 1.0e-4_dbl_kind,         &
        eps11  = 1.0e-11_dbl_kind,        &
        eps12  = 1.0e-12_dbl_kind,        &
        eps13  = 1.0e-13_dbl_kind,        &
        eps15  = 1.0e-15_dbl_kind,        &
        eps16  = 1.0e-16_dbl_kind,        &
!c     &  puny = eps13
!        puny = eps11


#  if !defined (DOUBLE_PRECISION)
        puny = 1.0e-5_dbl_kind
#  else
        puny = 1.0e-8_dbl_kind
# endif


!#ifndef CCSMcoupled
!      real (kind=dbl_kind) :: pi
!#endif

      real (kind=dbl_kind) :: & 
        pih !, pi2         ! pi/2 and 2pi

      !-----------------------------------------------------------------
      ! conversion factors
      !-----------------------------------------------------------------

      real (kind=dbl_kind), parameter :: & 
        cm_to_m       = 0.01_dbl_kind      &  ! cm to meters
      , m_to_cm       = 100._dbl_kind      &  ! meters to cm
      , m2_to_km2     = 1.e-6_dbl_kind     &  ! m^2 to km^2
      , kg_to_g       = 1000._dbl_kind     &  ! kilograms to grams
      , mps_to_cmpdy  = 8.64e6_dbl_kind &! m per s to cm per day
      , mps_to_cmpyr  = mps_to_cmpdy*365._dbl_kind ! m per s to cm per yr

      real (kind=dbl_kind) :: & 
        rad_to_deg               ! degree-radian conversion

!=======================================================================

      contains

!=======================================================================
!BOP
!
! !IROUTINE: init_constants - initialize constants defined at run time
!
! !INTERFACE:
!
      subroutine init_constants
!
! !DESCRIPTION:
!
! Initializes constants that are best defined at run time (e.g. pi)
!
! !REVISION HISTORY:
!
! author Elizabeth C. Hunke, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) :: i, j

      do j=jlo,jhi
       do i=ilo,ihi
!        one(i,j) = c1i
       enddo
      enddo

!#ifndef CCSMcoupled
!      pi  = c4*atan(c1)
!#endif
!      pi2 = c2*pi
      pih = p5*pi

      rad_to_deg = 180.0_dbl_kind/pi

      end subroutine init_constants

!=======================================================================

      end module ice_constants

!=======================================================================
