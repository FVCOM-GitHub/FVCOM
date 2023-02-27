!/===========================================================================/
! CVS VERSION INFORMATION
! $Id$
! $Name$
! $Revision$
!/===========================================================================/

!=======================================================================
!BOP
!
! !MODULE: ice_atmo - atm-ice interface: stability based flux calculations
!
! !DESCRIPTION:
!
! Atmospheric boundary interface (stability based flux calculations)
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
!
! Vectorized by Clifford Chen (Fujitsu) and William H. Lipscomb (LANL)
!
! !INTERFACE:
!
      module ice_atmo
!
! !USES:
!
      use ice_domain
      use ice_constants
      use ice_flux
      use ice_state
!
!EOP
!
      implicit none

!=======================================================================
 
      contains

!=======================================================================
!BOP
!
! !IROUTINE: atmo_boundary_layer - compute coefficients for atm-ice fluxes, stress and Tref
!
! !INTERFACE:
!
      subroutine atmo_boundary_layer (ni, sfctype, Tsf,  &
                                     strx, stry, Trf, Qrf, delt, delq)
!
! !DESCRIPTION:
!
! Compute coefficients for atm/ice fluxes, stress, and reference 
! temperature. NOTE: \\
! (1) all fluxes are positive downward, \\
! (2) here, tstar = (WT)/U*, and qstar = (WQ)/U*, \\
! (3) wind speeds should all be above a minimum speed (eg. 1.0 m/s). \\
!
! ASSUME: \\
!  The saturation humidity of air at T(K): qsat(T)  (kg/m**3) \\
!
! Code originally based on CSM1 \\
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
!
! !USES:
!
      use ice_grid       ! for tmask
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
        ni        ! thickness category index

      character (len=3), intent(in) ::&
        sfctype  ! ice or ocean

      real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi), intent(in) :: &
        Tsf      ! surface temperature of ice or ocean

      real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi), intent(out) ::&
         strx     &! x surface stress (N)
      ,  stry     &! y surface stress (N)
      ,  Trf      &! reference height temperature  (K)
      ,  Qrf      &! reference height specific humidity (kg/kg)
      ,  delt     &! potential T difference   (K)
      ,  delq     ! humidity difference      (kg/kg)
!
!EOP
!
       integer (kind=int_kind) ::&
         k       & ! iteration index
      ,  i, j     ! horizontal indices

      real (kind=dbl_kind) :: &
         TsfK     &! surface temperature in Kelvin (K)
      ,  xqq      &! temporary variable
      ,  psimh    &! stability function at zlvl   (momentum)
      ,  tau      &! stress at zlvl
      ,  fac      &! interpolation factor
      ,  al2      &! ln(z10   /zTrf)
      ,  psix2    &! stability function at zTrf   (heat and water)
      ,  psimhs   &! stable profile 
      ,  ssq      &! sat surface humidity     (kg/kg)
      ,  qqq      &! for qsat, dqsatdt
      ,  TTT      &! for qsat, dqsatdt
      ,  qsat     &! the saturation humidity of air (kg/m^3)
      ,  Lheat    ! Lvap or Lsub, depending on surface type

      real (kind=dbl_kind), dimension (1:(ihi-ilo+1)*(jhi-jlo+1)) ::&
         ustar    &! ustar (m/s)
      ,  tstar    &! tstar
      ,  qstar    &! qstar
      ,  rdn      &! sqrt of neutral exchange coefficient (momentum)
      ,  rhn      &! sqrt of neutral exchange coefficient (heat)
      ,  ren      &! sqrt of neutral exchange coefficient (water)
      ,  rd       &! sqrt of exchange coefficient (momentum)
      ,  re       &! sqrt of exchange coefficient (water)
      ,  rh       &! sqrt of exchange coefficient (heat)
      ,  vmag     &! surface wind magnitude   (m/s)
      ,  alz      &! ln(zlvl  /z10)
      ,  thva     &! virtual temperature      (K)
      ,  cp       &! specific heat of moist air
      ,  hol      &! H (at zlvl  ) over L
      ,  stable   &! stability factor
      ,  psixh    ! stability function at zlvl   (heat and water)

      integer (kind=int_kind), dimension (1:(ihi-ilo+1)*(jhi-jlo+1)) ::&
         indxi   & ! compressed index in i-direction
      ,  indxj    ! compressed index in j-direction

      integer (kind=int_kind) ::&
         icells  & ! number of cells with aicen > puny or tmask true
      ,  ij       ! combined ij index

      real (kind=dbl_kind), parameter ::   &
         cpvir = cp_wv/cp_air - c1i &! defined as cp_wv/cp_air - 1.
      ,  zTrf  = c2i                &! reference height for air temp (m)
      ,  umin  = c1i                 ! minimum wind speed (m/s)

      ! local functions
      real (kind=dbl_kind) :: &
         xd      & ! dummy argument  
      ,  psimhu  & ! unstable part of psimh
      ,  psixhu   ! unstable part of psimx

      !------------------------------------------------------------
      ! Define functions
      !------------------------------------------------------------

      psimhu(xd)  = log((c1i+xd*(c2i+xd))*(c1i+xd*xd)/c8i)  &
                  - c2i*atan(xd) + pi*p5
!ech     &            - c2i*atan(xd) + 1.571_dbl_kind

      psixhu(xd)  =  c2i * log((c1i + xd*xd)/c2i)

      al2 = log(zref/zTrf)

      !------------------------------------------------------------
      ! Initialize
      !------------------------------------------------------------

      do j = jlo,jhi
      do i = ilo,ihi
         lhcoef(i,j) = c0i
         shcoef(i,j) = c0i
         strx(i,j)   = c0i
         stry(i,j)   = c0i
         Trf(i,j)    = c0i
         Qrf(i,j)    = c0i
         delt(i,j)   = c0i
         delq(i,j)   = c0i
      enddo
      enddo

      !------------------------------------------------------------
      ! Compute turbulent flux coefficients, wind stress, and 
      ! reference temperature and humidity.
      !------------------------------------------------------------

      if ( sfctype(1:3)=='ice' ) then

      !------------------------------------------------------------
      ! identify grid cells with ice
      !------------------------------------------------------------
         icells = 0
         do j = jlo,jhi
         do i = ilo,ihi
            if ( aicen(i,j,ni) > puny) then
               icells = icells + 1
               indxi(icells) = i
               indxj(icells) = j
            endif
         enddo
         enddo

      !------------------------------------------------------------
      ! define some needed variables
      !------------------------------------------------------------
         qqq  = qqqice          ! for qsat
         TTT  = TTTice          ! for qsat
         Lheat = Lsub           ! ice to vapor
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)
            vmag(ij) = max(umin, wind(i,j))
            rdn(ij)  = vonkar/log(zref/iceruf) ! neutral coefficient
         enddo   ! ij

      elseif ( sfctype(1:3)=='ocn' ) then

      !------------------------------------------------------------
      ! identify ocean cells
      !------------------------------------------------------------
         icells = 0
         do j = jlo,jhi
         do i = ilo,ihi
            if ( tmask(i,j) ) then
               icells = icells + 1
               indxi(icells) = i
               indxj(icells) = j
            endif
         enddo
         enddo

       !------------------------------------------------------------
       ! define some needed variables
       !------------------------------------------------------------
         qqq  = qqqocn
         TTT  = TTTocn
         Lheat = Lvap       ! liquid to vapor
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)
            vmag(ij) = max(umin, wind(i,j))
            rdn(ij)  = sqrt(0.0027_dbl_kind/vmag(ij)    &
                    + .000142_dbl_kind + .0000764_dbl_kind*vmag(ij))
         enddo   ! ij

      endif   ! sfctype

      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)

      !------------------------------------------------------------
      ! define some more needed variables
      !------------------------------------------------------------

         TsfK       = Tsf(i,j) + Tffresh     ! surface temp (K)
         qsat       = qqq * exp(-TTT/TsfK)   ! saturation humidity (kg/m^3)
         ssq        = qsat / rhoa(i,j)       ! sat surf hum (kg/kg)

         thva(ij)   = potT(i,j) * (c1i + zvir * Qa(i,j)) ! virtual pot temp (K)
         delt(i,j)  = potT(i,j) - TsfK       ! pot temp diff (K)
         delq(i,j)  = Qa(i,j) - ssq          ! spec hum dif (kg/kg)
         alz(ij)    = log(zlvl(i,j)/zref) 
         cp(ij)     = cp_air*(c1i + cpvir*ssq)

      !------------------------------------------------------------
      ! first estimate of Z/L and ustar, tstar and qstar
      !------------------------------------------------------------

         ! neutral coefficients, z/L = 0.0 
         rhn(ij) = rdn(ij)
         ren(ij) = rdn(ij)

         ! ustar,tstar,qstar
         ustar(ij) = rdn(ij) * vmag(ij)
         tstar(ij) = rhn(ij) * delt(i,j)
         qstar(ij) = ren(ij) * delq(i,j)

      enddo                     ! ij

      !------------------------------------------------------------
      ! iterate to converge on Z/L, ustar, tstar and qstar
      !------------------------------------------------------------

      do k=1,5

         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)

        ! compute stability & evaluate all stability functions 
            hol(ij) = vonkar * gravit * zlvl(i,j)                     &
                   * (tstar(ij)/thva(ij)                              &
                    + qstar(ij)/(c1i/zvir+Qa(i,j)))                    &
                   / ustar(ij)**2                                      
            hol(ij)    = sign( min(abs(hol(ij)),c10i), hol(ij) )        
            stable(ij) = p5 + sign(p5 , hol(ij))                       
            xqq    = max(sqrt(abs(c1i - c16*hol(ij))) , c1i)             
            xqq    = sqrt(xqq)                                         
                                                                       
            ! Jordan et al 1999                                        
            psimhs = -(0.7_dbl_kind*hol(ij)                           &
                   + 0.75_dbl_kind*(hol(ij)-14.3_dbl_kind)            &
                   * exp(-0.35_dbl_kind*hol(ij)) + 10.7_dbl_kind)      
            psimh  = psimhs*stable(ij)                                &
                    + (c1i - stable(ij))*psimhu(xqq)                    
            psixh(ij)  = psimhs*stable(ij)                            &
                    + (c1i - stable(ij))*psixhu(xqq)

        ! shift all coeffs to measurement height and stability
            rd(ij) = rdn(ij) / (c1i+rdn(ij)/vonkar*(alz(ij)-psimh))
            rh(ij) = rhn(ij) / (c1i+rhn(ij)/vonkar*(alz(ij)-psixh(ij)))
            re(ij) = ren(ij) / (c1i+ren(ij)/vonkar*(alz(ij)-psixh(ij)))

        ! update ustar, tstar, qstar using updated, shifted coeffs 
            ustar(ij) = rd(ij) * vmag(ij)
            tstar(ij) = rh(ij) * delt(i,j)
            qstar(ij) = re(ij) * delq(i,j)

         enddo                  ! ij
      enddo                     ! end iteration


      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)

      !------------------------------------------------------------
      ! coefficients for turbulent flux calculation
      !------------------------------------------------------------

!ccsm      shcoef(i,j) = rhoa(i,j) * ustar(ij) * cp(ij) * rh(ij)
      ! add windless coefficient for sensible heat flux
      ! as in Jordan et al (JGR, 1999)
         shcoef(i,j) = rhoa(i,j) * ustar(ij) * cp(ij) * rh(ij) + c1i
         lhcoef(i,j) = rhoa(i,j) * ustar(ij) * Lheat  * re(ij)
      !------------------------------------------------------------
      ! momentum flux
      !------------------------------------------------------------
      ! tau = rhoa(i,j) * ustar * ustar 
      ! strx = tau * uatm(i,j) / vmag 
      ! stry = tau * vatm(i,j) / vmag 
      !------------------------------------------------------------

         tau = rhoa(i,j) * ustar(ij) * rd(ij) ! not the stress at zlvl(i,j)
         strx(i,j) = tau * uatm(i,j)
         stry(i,j) = tau * vatm(i,j)

      !------------------------------------------------------------
      ! Compute diagnostics: 2m ref T & Q
      !------------------------------------------------------------
         hol(ij) = hol(ij)*zTrf/zlvl(i,j)
         xqq     = max( c1i, sqrt(abs(c1i-c16*hol(ij))) )
         xqq     = sqrt(xqq)
         psix2   = -c5i*hol(ij)*stable(ij) + (c1i-stable(ij))*psixhu(xqq)  
         fac     = (rh(ij)/vonkar) * (alz(ij) + al2 - psixh(ij) + psix2)
         Trf(i,j)= potT(i,j) - delt(i,j)*fac
         Trf(i,j)= Trf(i,j) - p01*zTrf ! pot temp to temp correction
         fac     = (re(ij)/vonkar) * (alz(ij) + al2 - psixh(ij) + psix2)
         Qrf(i,j)= Qa(i,j) - delq(i,j)*fac

      enddo                     ! ij

      end subroutine atmo_boundary_layer

!=======================================================================

      end module ice_atmo

!=======================================================================
