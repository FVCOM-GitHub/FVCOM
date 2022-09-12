!/===========================================================================/
! CVS VERSION INFORMATION
! $Id$
! $Name$
! $Revision$
!/===========================================================================/

!=======================================================================
!BOP
!
! !MODULE: ice_mechred - mechanical redestribution and ice strength
!
! !DESCRIPTION:
!
! Ice mechanical redistribution (ridging) and strength computations 
!
! See these references:
!
! Flato, G. M., and W. D. Hibler III, 1995: Ridging and strength
!  in modeling the thickness distribution of Arctic sea ice,
!  J. Geophys. Res., 100, 18,611-18,626.
!
! Hibler, W. D. III, 1980: Modeling a variable thickness sea ice
!  cover, Mon. Wea. Rev., 108, 1943-1973, 1980.
!
! Lipscomb, W. H., E. C. Hunke, W. Maslowski, and J. Jakacki, 2006:
!  Ridging, strength, and stability in sea ice models, submitted
!  to J. Geophys. Res.
!
! Rothrock, D. A., 1975: The energetics of the plastic deformation of
!  pack ice by ridging, J. Geophys. Res., 80, 4514-4519.
!
! Thorndike, A. S., D. A. Rothrock, G. A. Maykut, and R. Colony, 
!  1975: The thickness distribution of sea ice, J. Geophys. Res., 
!  80, 4501-4513. 
!
! !REVISION HISTORY:
!
! authors William H. Lipscomb, LANL
!         Elizabeth C. Hunke (LANL)
!
! 2004: Vectorized by Clifford Chen (Fujitsu) and William Lipscomb
! 2006: New options for participation and redistribution (WHL)
!
! !INTERFACE:
!
      module ice_mechred
!
! !USES:
! 
      use ice_model_size
      use ice_constants
      use ice_state
      use ice_itd
      use ice_grid
      use ice_fileunits
      use ice_domain
      use ice_calendar, only: istep1, dyn_dt
      use ice_work, only:  worka
!
!EOP
!
      implicit none
      save

!-----------------------------------------------------------------------
! Ridging parameters
!-----------------------------------------------------------------------

      integer (kind=int_kind) :: &  ! defined in namelist
         kstrength       & ! 0 for simple Hibler (1979) formulation 
                           ! 1 for Rothrock (1975) pressure formulation 
      ,  krdg_partic     & ! 0 for Thorndike et al. (1975) formulation 
                           ! 1 for exponential participation function 
      ,  krdg_redist       ! 0 for Hibler (1980) formulation
                          ! 1 for exponential redistribution function
      
      real (kind=dbl_kind), parameter :: &
         Cf = 17._dbl_kind  & ! ratio of ridging work to PE change in ridging
      ,  Cs = p25           & ! fraction of shear energy contrbtng to ridging
      ,  Cp = p5*gravit*(rhow-rhoi)*rhoi/rhow  &! proport const for PE
      ,  fsnowrdg = p5           &! snow fraction that survives in ridging
      ,  Gstar = p15             &! max value of G(h) that participates
                                  ! (krdg_partic = 0)
      ,  astar = 0.05_dbl_kind   &! e-folding scale for G(h) participation
                                  ! (krdg_partic = 1)
      ,  maxraft = c1i            &! max value of hrmin - hi = max thickness
                                  ! of ice that rafts (m)
      ,  Hstar = c25             &! determines mean thickness of ridged ice (m)
                                  ! (krdg_redist = 0)
                                  ! Flato & Hibler (1995) have Hstar = 100
# if defined (ICE_FRESHWATER)
! afm 20160613 & EJA 20160921 -  best value suggested from the sensitivity study 
! from mu_rdg=0.5-4
      ,  mu_rdg = c2i             &! gives e-folding scale of ridged ice (m^.5)
                                  ! (krdg_redist = 1)
# else
      ,  mu_rdg = c4i             &! gives e-folding scale of ridged ice (m^.5)
                                  ! (krdg_redist = 1)
# endif
      ,  Pstar = 2.75e4_dbl_kind &! constant in Hibler strength formula
                                  ! (kstrength = 0)
      ,  Cstar = c20              ! constant in Hibler strength formula
                                 ! (kstrength = 0)

!-----------------------------------------------------------------------
!     Ridging diagnostic arrays for history files
!-----------------------------------------------------------------------

!      real (kind=dbl_kind), dimension(ilo:ihi,jlo:jhi) :: &
!     &   dardg1dt         ! rate of fractional area loss by ridging ice (1/s)
!     &,  dardg2dt         ! rate of fractional area gain by new ridges (1/s)
!     &,  dvirdgdt         ! rate of ice volume ridged (m/s)
!     &,  opening          ! rate of opening due to divergence/shear (1/s)

      real (kind=dbl_kind), dimension(:,:),allocatable,save :: &
         dardg1dt        & ! rate of fractional area loss by ridging ice (1/s)
      ,  dardg2dt        & ! rate of fractional area gain by new ridges (1/s)
      ,  dvirdgdt        & ! rate of ice volume ridged (m/s)
      ,  opening           ! rate of opening due to divergence/shear (1/s)



!-----------------------------------------------------------------------
! Variables shared among ridging subroutines
!-----------------------------------------------------------------------

!      real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi) :: &    
      real (kind=dbl_kind), dimension (:,:),allocatable,save :: &    
         asum           &  ! sum of total ice and open water area
      ,  aksum             ! ratio of area removed to area ridged

!      real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi,0:ncat) :: &    
      real (kind=dbl_kind), dimension (:,:,:) ,allocatable,save:: &    
         apartic          ! participation function; fraction of ridging/
                          !  closing associated w/ category n

!      real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi,ncat) :: &    
      real (kind=dbl_kind), dimension (:,:,:),allocatable,save :: &    
         hrmin           & ! minimum ridge thickness
      ,  hrmax           & ! maximum ridge thickness
      ,  hrexp           & ! ridge e-folding thickness (krdg_redist = 1) 
      ,  krdg              ! mean ridge thickness/thickness of ridging ice 

      logical (kind=log_kind), parameter :: &
         l_conservation_check = .true.  ! if true, check conservation 
                                        ! (useful for debugging)
                                       
!=======================================================================

      contains

!=======================================================================
!BOP
!
! !ROUTINE: init_mechred
!
! !DESCRIPTION:
!
! Initialize some variables written to the history file
! 
! !REVISION HISTORY:
!
! author Elizabeth C. Hunke, LANL
!
! !INTERFACE:
!
      subroutine init_mechred
!
! !USES:
! 
! !INPUT/OUTPUT PARAMETERS:
!
!
!EOP
!
      integer (kind=int_kind) :: i,j

      !-----------------------------------------------------------------
      ! Initialize history fields.
      !-----------------------------------------------------------------

      dardg1dt(:,:) = c0i
      dardg2dt(:,:) = c0i
      dvirdgdt(:,:) = c0i
      opening (:,:) = c0i

      end subroutine init_mechred

!=======================================================================
!BOP
!
! !ROUTINE: ridge_ice - driver for mechanical redistribution
!
! !DESCRIPTION:
!
! Compute changes in the ice thickness distribution due to divergence
! and shear.
!
! !REVISION HISTORY:
!
! authors William H. Lipscomb, LANL
!
! !INTERFACE:
!
      subroutine ridge_ice (Delta, divu)
!
! !USES:
! 
!      use ice_timers

      use ice_flux
!      use ice_exit
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind),                                &
         dimension (imt_local,jmt_local), intent(in) ::    &
         Delta   & ! term in the denominator of zeta, eta          (1/s)
                  ! = sqrt (epsI^2 + (1/e^2)*epsII^2)
      ,  divu     ! strain rate I component, velocity divergence  (1/s)
!
!EOP
!
      integer (kind=int_kind), parameter :: &
        nitermax = 20     ! max number of ridging iterations

      integer (kind=int_kind) :: &
         i,j               &! horizontal indices
      ,  ni                 &! thickness category index
      ,  niter              ! iteration counter
                            
      real (kind=dbl_kind),  dimension(ilo:ihi,jlo:jhi) :: &
         closing_net       &! net rate at which area is removed    (1/s)
                            ! (ridging ice area - area of new ridges) / dyn_dt
      ,  divu_adv          &! divu as implied by transport scheme  (1/s)
      ,  opning            &! rate of opening due to divergence/shear
      ,  closing_gross     &! rate at which area removed, not counting
                            ! area of new ridges
      ,  msnow_mlt         &! mass of snow added to ocean (kg m-2)
      ,  esnow_mlt          ! energy needed to melt snow in ocean (J m-2)

      real (kind=dbl_kind) :: &     
         w1              &  ! temporary variable
      ,  tmpfac          &  ! factor by which opening/closing rates are cut
!      ,  dti                ! 1 / dyn_dt
      ,  dti_ice                ! 1 / dyn_dt  ! conflict  with main dti
!  ggao 

      logical (kind=log_kind) :: &
         iterate_ridging   ! if true, repeat the ridging

      real (kind=dbl_kind), parameter :: &    
         big = 1.0e+8_dbl_kind

      logical (kind=log_kind) ::      &
         asum_error        ! flag for asum .ne. 1

!      call ice_timer_start(6)  ! ridging 

      !-----------------------------------------------------------------
      ! Set hin_max(ncat) to a big value to ensure that all ridged ice 
      ! is thinner than hin_max(ncat).
      !-----------------------------------------------------------------
      hin_max(ncat) = big

      !-----------------------------------------------------------------
      ! Compute the ice strength, the thickness distribution of ridging 
      ! ice, and various quantities associated with the new ridged ice.
      !-----------------------------------------------------------------
      call ridge_prep

      !-----------------------------------------------------------------
      ! Compute total area of ice plus open water.
      ! This may not be equal to one because of divergence during
      !  transport.
      !-----------------------------------------------------------------
      call asum_ridging

      !-----------------------------------------------------------------
      ! Initialize arrays.
      !-----------------------------------------------------------------

      msnow_mlt(:,:) = c0i
      esnow_mlt(:,:) = c0i
      dardg1dt (:,:) = c0i
      dardg2dt (:,:) = c0i
      dvirdgdt (:,:) = c0i
      opening  (:,:) = c0i

      do j = jlo, jhi
      do i = ilo, ihi

      !-----------------------------------------------------------------
      ! Compute the net rate of closing due to convergence and shear,
      ! based on Flato and Hibler (1995).
      ! 
      ! The energy dissipation rate is equal to the net closing rate
      ! times the ice strength.
      !
      ! NOTE: The NET closing rate is equal to the rate that open water 
      !  area is removed, plus the rate at which ice area is removed by 
      !  ridging, minus the rate at which area is added in new ridges.
      !  The GROSS closing rate is equal to the first two terms (open
      !  water closing and thin ice ridging) without the third term
      !  (thick, newly ridged ice).
      !-----------------------------------------------------------------

         closing_net(i,j) =            &
              Cs*p5*(Delta(i,j)-abs(divu(i,j))) - min(divu(i,j),c0i)

      !-----------------------------------------------------------------
      ! Compute divu_adv, the divergence rate given by the transport/
      ! advection scheme, which may not be equal to divu as computed 
      ! from the velocity field.
      !
      ! If divu_adv < 0, make sure the closing rate is large enough
      ! to give asum = 1.0 after ridging.
      !-----------------------------------------------------------------

         divu_adv(i,j) = (c1i-asum(i,j)) / dyn_dt

         if (divu_adv(i,j) < c0i)     &
              closing_net(i,j) = max(closing_net(i,j), -divu_adv(i,j))

      !-----------------------------------------------------------------
      ! Compute the (non-negative) opening rate that will give 
      ! asum = 1.0 after ridging.
      !-----------------------------------------------------------------
         opning(i,j) = closing_net(i,j) + divu_adv(i,j)
      enddo
      enddo


      do niter = 1, nitermax    ! iteration counter

      do j = jlo, jhi
      do i = ilo, ihi

      !-----------------------------------------------------------------
      ! Based on the ITD of ridging and ridged ice, convert the net
      !  closing rate to a gross closing rate.  
      ! NOTE: 0 < aksum <= 1
      !-----------------------------------------------------------------


         closing_gross(i,j) = closing_net(i,j) / aksum(i,j)

      !-----------------------------------------------------------------
      ! Reduce the closing rate if more than 100% of the open water 
      ! would be removed.  Reduce the opening rate proportionately.
      !-----------------------------------------------------------------

         if (apartic(i,j,0) > c0i) then
            w1 = apartic(i,j,0) * closing_gross(i,j) * dyn_dt
            if (w1 > aice0(i,j)) then
               tmpfac = aice0(i,j) / w1
               closing_gross(i,j) = closing_gross(i,j) * tmpfac
               opning(i,j) = opning(i,j) * tmpfac
            endif
         endif
      enddo                     ! i
      enddo                     ! j

      !-----------------------------------------------------------------
      ! Reduce the closing rate if more than 100% of any ice category 
      ! would be removed.  Reduce the opening rate proportionately.
      !-----------------------------------------------------------------
      do ni = 1, ncat
         do j = jlo, jhi
         do i = ilo, ihi

            if (aicen(i,j,ni) > puny .and. apartic(i,j,ni) > c0i) then
               w1 = apartic(i,j,ni) * closing_gross(i,j) * dyn_dt
               if (w1 > aicen(i,j,ni)) then
                  tmpfac = aicen(i,j,ni) / w1
                  closing_gross(i,j) = closing_gross(i,j) * tmpfac
                  opning(i,j) = opning(i,j) * tmpfac
               endif
            endif

         enddo                  ! i
         enddo                  ! j
      enddo                     ! n

      !-----------------------------------------------------------------
      ! Redistribute area, volume, and energy.
      !-----------------------------------------------------------------
      call ridge_shift (opning,    closing_gross,     &
                        msnow_mlt, esnow_mlt)

      !-----------------------------------------------------------------
      ! Compute total area of ice plus open water after ridging.
      !-----------------------------------------------------------------
      call asum_ridging

      !-----------------------------------------------------------------
      ! Check whether asum = 1.  If not (because the closing and opening
      ! rates were reduced above), ridge again with new rates.
      !-----------------------------------------------------------------

      iterate_ridging = .false.

      do j = jlo, jhi
      do i = ilo, ihi
         if (abs(asum(i,j) - c1i) < puny) then
            closing_net(i,j) = c0i   ! no ridging the next time through
            opning(i,j) = c0i
         else
            iterate_ridging = .true.
            divu_adv(i,j) = (c1i - asum(i,j)) / dyn_dt
            closing_net(i,j) = max(c0i, -divu_adv(i,j))
            opning(i,j) = max(c0i, divu_adv(i,j))
         endif
      enddo
      enddo

      !-----------------------------------------------------------------
      ! Repeat if necessary.
      !-----------------------------------------------------------------

      if (iterate_ridging) then
         if (niter > nitermax) then
            write(nu_diag,*) 'istep1, my_task, nitermax =',  &
                              istep1, my_task, nitermax
!            call abort_ice('Exceeded max number of ridging iterations')
         endif
!         write(nu_diag,*) 'REPEAT RIDGING, istep1, my_task, niter =', &
!                                           istep1, my_task, niter
         call ridge_prep
      else
         exit
      endif

      enddo                     ! niter

      !-----------------------------------------------------------------
      ! Convert ridging rate diagnostics to correct units.
      !
      ! Update fresh water and heat fluxes due to snow melt.
      !-----------------------------------------------------------------

      dti_ice = c1i/dyn_dt

      asum_error = .false. 

      do j = jlo, jhi
      do i = ilo, ihi

         if (abs(asum(i,j) - c1i) > puny) asum_error = .true.

         dardg1dt(i,j) = dardg1dt(i,j) * dti_ice
         dardg2dt(i,j) = dardg2dt(i,j) * dti_ice
         dvirdgdt(i,j) = dvirdgdt(i,j) * dti_ice
         opening (i,j) = opening (i,j) * dti_ice

         ! fresh water source for ocean
         fresh(i,j)      = fresh(i,j)      + msnow_mlt(i,j)*dti_ice
         fresh_hist(i,j) = fresh_hist(i,j) + msnow_mlt(i,j)*dti_ice
      
         ! heat sink for ocean
         fhnet(i,j)      = fhnet(i,j)      + esnow_mlt(i,j)*dti_ice
         fhnet_hist(i,j) = fhnet_hist(i,j) + esnow_mlt(i,j)*dti_ice

      enddo
      enddo

      !-----------------------------------------------------------------
      ! Abort if area does not add up to one.
      !-----------------------------------------------------------------

      if (asum_error) then
         do j = jlo, jhi
         do i = ilo, ihi
            if (abs(asum(i,j) - c1i) > puny) then ! there is a bug
               write(nu_diag,*) ' '
               write(nu_diag,*) 'Ridging error: total area =', asum(i,j)
               write(nu_diag,*) 'istep1, my_task, i, j:',  &
                                istep1, my_task, i, j
              !                  istep1, my_task, i, ngid(j),aice(i,j)
               write(nu_diag,*) 'n, aicen, apartic:'
               write(nu_diag,*)  0, aice0(i,j), apartic(i,j,0)
               do ni = 1, ncat
                  write(nu_diag,*) ni, aicen(i,j,ni), apartic(i,j,ni)
               enddo
!               call abort_ice('ridging: total area must be <= 1')
            endif
         enddo
         enddo
      endif

!      call ice_timer_stop(6)  ! ridging 

      end subroutine ridge_ice

!=======================================================================
!BOP
!
! !ROUTINE: ice_strength - compute ice strength
!
! !DESCRIPTION:
!
! Compute the strength of the ice pack, defined as the energy (J m-2) 
! dissipated per unit area removed from the ice pack under compression,
! and assumed proportional to the change in potential energy caused
! by ridging.
!
! See Rothrock (1975) and Hibler (1980).
!
! For simpler strength parameterization, see this reference:
! Hibler, W. D. III, 1979: A dynamic-thermodynamic sea ice model,
!  J. Phys. Oceanog., 9, 817-846.
!
! !REVISION HISTORY:
!
! authors William H. Lipscomb, LANL
!         Elizabeth C. Hunke, LANL
!
! !INTERFACE:
!
      subroutine ice_strength (kstrngth)
!
! !USES:
! 
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
         kstrngth    ! = 1 for Rothrock formulation, 0 for Hibler (1979)
!
!EOP
!
      integer (kind=int_kind) :: &
         i,j           &   ! horizontal indices
      ,  ni                ! thickness category index

      real (kind=dbl_kind) :: &    
         hi              & ! ice thickness (m)
      ,  h2rdg           & ! mean value of h^2 in new ridge
      ,  dh2rdg            ! change in mean value of h^2 per unit area
                          ! consumed by ridging 


      ! initialize
      strength(:,:) = c0i   

      !-----------------------------------------------------------------
      ! Compute thickness distribution of ridging and ridged ice.
      !-----------------------------------------------------------------
      call ridge_prep 

      if (kstrngth == 1) then

      !-----------------------------------------------------------------
      ! Compute ice strength based on change in potential energy,
      ! as in Rothrock (1975)
      !-----------------------------------------------------------------

         if (krdg_redist==0) then ! Hibler redistribution function

            do ni = 1, ncat
               do j = jlo, jhi
               do i = ilo, ihi

                  if(aicen(i,j,ni) > puny .and. apartic(i,j,ni) > c0i) then
                     hi = vicen(i,j,ni) / aicen(i,j,ni)
                     h2rdg = p333 * (hrmax(i,j,ni)**3 - hrmin(i,j,ni)**3) &
                                  / (hrmax(i,j,ni) - hrmin(i,j,ni))        
                     dh2rdg = -hi*hi + h2rdg/krdg(i,j,ni)                 
                     strength(i,j) = strength(i,j)                      &
                                   + apartic(i,j,ni) * dh2rdg
                  endif
               enddo            ! i
               enddo            ! j
            enddo               ! n

         elseif (krdg_redist==1) then  ! exponential function

            do ni = 1, ncat
               do j = jlo, jhi
               do i = ilo, ihi

                  if(aicen(i,j,ni) > puny .and. apartic(i,j,ni) > c0i) then
                     hi = vicen(i,j,ni) / aicen(i,j,ni)
                     h2rdg =    hrmin(i,j,ni)*hrmin(i,j,ni)  &
                           + c2i*hrmin(i,j,ni)*hrexp(i,j,ni)  &
                           + c2i*hrexp(i,j,ni)*hrexp(i,j,ni)   
                     dh2rdg = -hi*hi + h2rdg/krdg(i,j,ni)    
                     strength(i,j) = strength(i,j)         &
                                   + apartic(i,j,ni) * dh2rdg
                  endif
               enddo            ! i
               enddo            ! j
            enddo               ! n

         endif                  ! krdg_redist

         do j = jlo, jhi
         do i = ilo, ihi

            strength(i,j) = Cf * Cp * strength(i,j) / aksum(i,j) 
                          ! Cp = (g/2)*(rhow-rhoi)*(rhoi/rhow)
                          ! Cf accounts for frictional dissipation

         enddo                  ! j
         enddo                  ! i

      else                      ! kstrngth ne 1:  Hibler (1979) form

      !-----------------------------------------------------------------
      ! Compute ice strength as in Hibler (1979)
      !-----------------------------------------------------------------
         do j = jlo, jhi
         do i = ilo, ihi
            strength(i,j) = Pstar*vice(i,j)*exp(-Cstar*(c1i-aice(i,j)))
         enddo                  ! j
         enddo                  ! i

      endif                     ! kstrngth

      end subroutine ice_strength

!=======================================================================
!BOP
!
! !ROUTINE: ridge_prep - preparation for ridging and strength calculations
!
! !DESCRIPTION:
!
! Compute the thickness distribution of the ice and open water 
! participating in ridging and of the resulting ridges.
!
! This version includes new options for ridging participation and
!  redistribution.
! The new participation scheme (krdg_partic = 1) improves model
!  stability by increasing the time scale for large changes in ice strength.
! The new exponential redistribution function (krdg_redist = 1) improves 
!  agreement between ITDs of modeled and observed ridges.   
!
! !REVISION HISTORY:
!
! authors William H. Lipscomb, LANL
! 
! 2006: Added new options for ridging participation and redistribution.  
!
! !INTERFACE:
!
      subroutine ridge_prep 
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) :: &
         i,j              &! horizontal indices
      ,  ni                ! thickness category index
                           
      real (kind=dbl_kind) , parameter :: &    
         Gstari   = c1i/Gstar     &
      ,  astari   = c1i/astar
                           
      real (kind=dbl_kind) , dimension(ilo:ihi,jlo:jhi,-1:ncat) :: &    
         Gsum              ! Gsum(n) = sum of areas in categories 0 to n
                           
      real (kind=dbl_kind)  :: &    
         hi               &! ice thickness for each cat (m)
      ,  hieff            &! effective ice thickness (m) (krdg_redist = 2)
      ,  hrmean           &! mean ridge thickness (m)
      ,  xtmp              ! temporary variable

      !-----------------------------------------------------------------
      ! Initialize 
      !----------------------------------------------------------------- 
 
      aksum(:,:) = c0i 
 
      do ni = 0, ncat 
         apartic(:,:,ni) = c0i 
      enddo 
 
      do ni = 1, ncat
         hrmin (:,:,ni) = c0i 
         hrmax (:,:,ni) = c0i
         hrexp (:,:,ni) = c0i
         krdg  (:,:,ni) = c1i
      enddo 

      !-----------------------------------------------------------------
      ! Compute the thickness distribution of ice participating in ridging.
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! First compute the cumulative thickness distribution function Gsum,
      !  where Gsum(n) is the fractional area in categories 0 to n.
      ! Ignore categories with very small areas.
      !-----------------------------------------------------------------

      Gsum(:,:,-1) = c0i

      do j = jlo, jhi
      do i = ilo, ihi
         if (aice0(i,j) > puny) then
            Gsum(i,j,0) = aice0(i,j)
         else
            Gsum(i,j,0) = Gsum(i,j,-1)
         endif
      enddo  
      enddo

      do ni = 1, ncat
         do j = jlo, jhi
         do i = ilo, ihi
            if (aicen(i,j,ni) > puny) then
               Gsum(i,j,ni) = Gsum(i,j,ni-1) + aicen(i,j,ni)
            else
               Gsum(i,j,ni) = Gsum(i,j,ni-1)
            endif
         enddo
         enddo
      enddo

      ! normalize

      worka(:,:) = c1i / Gsum(:,:,ncat)

      do ni = 0, ncat
         do j = jlo, jhi
         do i = ilo, ihi
            Gsum(i,j,ni) = Gsum(i,j,ni) * worka(i,j)
         enddo
         enddo
      enddo

      !-----------------------------------------------------------------
      ! Compute the participation function apartic; this is analogous to
      ! a(h) = b(h)g(h) as defined in Thorndike et al. (1975).
      !
      !                area lost from category n due to ridging/closing
      !  apartic(n) = ---------------------------------------------------
      !                    total area lost due to ridging/closing
      !
      !-----------------------------------------------------------------

      if (krdg_partic==0) then  ! Thorndike et al. 1975

      !-----------------------------------------------------------------
      ! b(h) = (2/Gstar) * (1 - G(h)/Gstar). 
      ! The expressions for apartic are found by integrating b(h)g(h)
      ! between the category boundaries.
      !-----------------------------------------------------------------

         do ni = 0, ncat
            do j = jlo, jhi
            do i = ilo, ihi

               if (Gsum(i,j,ni) < Gstar) then
                  apartic(i,j,ni) = Gstari * (Gsum(i,j,ni)-Gsum(i,j,ni-1))    &
                            * (c2i - (Gsum(i,j,ni-1)+Gsum(i,j,ni))*Gstari)     
               elseif (Gsum(i,j,ni-1) < Gstar) then                          
                  apartic(i,j,ni) = Gstari * (Gstar-Gsum(i,j,ni-1))          &
                            * (c2i - (Gsum(i,j,ni-1)+Gstar)*Gstari)
               endif

            enddo               ! i
            enddo               ! j
         enddo                  ! ni

      elseif (krdg_partic==1) then ! exponential dependence on G(h)

      !-----------------------------------------------------------------
      ! b(h) = exp(-G(h)/astar)
      ! apartic(n) = [exp(-G(ni-1)/astar - exp(-G(n)/astar] / [1-exp(-1/astar)]. 
      ! The expression for apartic is found by integrating b(h)g(h)
      ! between the category boundaries.
      !-----------------------------------------------------------------

         ! precompute exponential terms using Gsum as work array
 
         xtmp = c1i / (c1i-exp(-astari))
 
         do ni = -1, ncat
            do j = jlo, jhi
            do i = ilo, ihi
               Gsum(i,j,ni) = exp(-Gsum(i,j,ni)*astari) * xtmp
            enddo               ! i
            enddo               ! j
         enddo                  ! n
 
         ! compute apartic

         do ni = 0, ncat
            do j = jlo, jhi
            do i = ilo, ihi
               apartic(i,j,ni) = Gsum(i,j,ni-1) - Gsum(i,j,ni)
            enddo               ! i
            enddo               ! j
         enddo                  ! n

      endif   ! krdg_partic

      !----------------------------------------------------------------- 
      ! Compute variables related to ITD of ridged ice: 
      ! 
      ! krdg = mean ridge thickness / thickness of ridging ice
      ! hrmin = min ridge thickness 
      ! hrmax = max ridge thickness (krdg_redist = 0) 
      ! hrexp = ridge e-folding scale (krdg_redist = 1) 
      !---------------------------------------------------------------- 

      if (krdg_redist == 0) then  ! Hibler 1980 formulation 
 
      !----------------------------------------------------------------- 
      ! Assume ridged ice is uniformly distributed between hrmin and hrmax. 
      ! 
      ! This parameterization is a modified version of Hibler (1980). 
      ! In the original paper the min ridging thickness is hrmin = 2*hi,
      !  and the max thickness is hrmax = 2*sqrt(hi*Hstar).
      !
      ! Here the min thickness is hrmin = min(2*hi, hi+maxraft),
      !  so thick ridging ice is not required to raft.
      ! 
      !----------------------------------------------------------------- 
 
         do ni = 1, ncat
            do j = jlo, jhi
            do i = ilo, ihi

               if (aicen(i,j,ni) > puny) then
                  hi = vicen(i,j,ni) / aicen(i,j,ni)
                  hrmin(i,j,ni) = min(c2i*hi, hi + maxraft)
                  hrmax(i,j,ni) = c2i*sqrt(Hstar*hi)
                  hrmax(i,j,ni) = max(hrmax(i,j,ni), hrmin(i,j,ni)+puny)
                  hrmean = p5 * (hrmin(i,j,ni) + hrmax(i,j,ni))
                  krdg(i,j,ni) = hrmean / hi
               endif

            enddo               ! i
            enddo               ! j
         enddo                  ! n

      else               ! krdg_redist = 1; exponential redistribution
 
      !----------------------------------------------------------------- 
      ! The ridge ITD is a negative exponential: 
      ! 
      !  g(h) ~ exp[-(h-hrmin)/hrexp], h >= hrmin 
      ! 
      ! where hrmin is the minimum thickness of ridging ice and 
      ! hrexp is the e-folding thickness.
      ! 
      ! Here, assume as above that hrmin = min(2*hi, hi+maxraft).
      ! That is, the minimum ridge thickness results from rafting,
      !  unless the ice is thicker than maxraft.
      !
      ! Also, assume that hrexp = mu_rdg*sqrt(hi).
      ! The parameter mu_rdg is tuned to give e-folding scales mostly
      !  in the range 2-4 m as observed by upward-looking sonar.
      !
      ! Values of mu_rdg in the right column give ice strengths
      !  roughly equal to values of Hstar in the left column
      !  (within ~10 kN/m for typical ITDs):
      !
      !   Hstar      mu_rdg
      !
      !     25        3.0
      !     50        4.0
      !     75        5.0
      !    100        6.0
      !----------------------------------------------------------------- 
 
         do ni = 1, ncat 
            do j = jlo, jhi
            do i = ilo, ihi 
               if (aicen(i,j,ni) > puny) then 
                  hi = vicen(i,j,ni) / aicen(i,j,ni)
                  hrmin(i,j,ni) = min(c2i*hi, hi + maxraft)
                  hrexp(i,j,ni) = mu_rdg * sqrt(hi)
                  krdg(i,j,ni) = (hrmin(i,j,ni) + hrexp(i,j,ni)) / hi 
               endif 
            enddo 
            enddo 
         enddo 
 
      endif                     ! krdg_redist 

      !----------------------------------------------------------------
      ! Compute aksum = net ice area removed / total area participating.
      ! For instance, if a unit area of ice with h = 1 participates in
      !  ridging to form a ridge with a = 1/3 and h = 3, then
      !  aksum = 1 - 1/3 = 2/3.
      !---------------------------------------------------------------- 

      do j = jlo, jhi
      do i = ilo, ihi
         aksum(i,j) = apartic(i,j,0)   ! area participating = area removed
      enddo
      enddo

      do ni = 1, ncat
         do j = jlo, jhi
         do i = ilo, ihi
            ! area participating > area removed
            aksum(i,j) = aksum(i,j)      &
                       + apartic(i,j,ni) * (c1i - c1i/krdg(i,j,ni)) 
         enddo
         enddo
      enddo

      end subroutine ridge_prep

!=======================================================================
!BOP
!
! !ROUTINE: ridge_shift - shift ridging ice among thickness categories
!
! !DESCRIPTION:
!
! Remove area, volume, and energy from each ridging category
! and add to thicker ice categories.
!
! !REVISION HISTORY:
!
! author William H. Lipscomb, LANL
!
! !INTERFACE:
!
      subroutine ridge_shift (opning,    closing_gross, &
                              msnow_mlt, esnow_mlt)
!
! !USES:
!
!      use ice_exit
!      ggao

! 
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind), dimension(ilo:ihi,jlo:jhi), intent(in) :: &
         opning          & ! rate of opening due to divergence/shear
      ,  closing_gross     ! rate at which area removed, not counting
                           ! area of new ridges

      real (kind=dbl_kind), dimension(ilo:ihi,jlo:jhi), intent(inout) :: &
         msnow_mlt      &  ! mass of snow added to ocean (kg m-2)
      ,  esnow_mlt        ! energy needed to melt snow in ocean (J m-2)
!
!EOP
!
      integer (kind=int_kind) :: &
         i,j              & ! horizontal indices
      ,  ni, n1, n2        & ! thickness category indices
      ,  k                & ! ice layer index
      ,  ij               & ! horizontal index, combines i and j loops
      ,  icells             ! number of cells with aicen > puny

      integer (kind=int_kind), dimension(1:(ihi-ilo+1)*(jhi-jlo+1)) :: &
         indxi, indxj      ! compressed indices

      real (kind=dbl_kind), dimension (imt_local,jmt_local) :: &    
         vice_init, vice_final & ! ice volume summed over categories
      ,  eice_init, eice_final  ! ice energy summed over layers

      real (kind=dbl_kind), dimension (imt_local,jmt_local,ncat) :: &    
         aicen_init       & ! ice area before ridging
      ,  vicen_init       & ! ice volume before ridging
      ,  vsnon_init       & ! snow volume before ridging
      ,  esnon_init         ! snow energy before ridging

      real (kind=dbl_kind), dimension (imt_local,jmt_local,ntilay) :: &    
        eicen_init        ! ice energy before ridging

      real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi) :: &    
         afrac            & ! fraction of category area ridged
      ,  ardg1            & ! area of ice ridged
      ,  ardg2            & ! area of new ridges
      ,  virdg            & ! ice volume of ridging ice
      ,  vsrdg            & ! snow volume of ridging ice
      ,  esrdg            & ! snow energy of ridging ice
      ,  dhr              & ! hrmax - hrmin
      ,  dhr2             & ! hrmax^2 - hrmin^2
      ,  farea            & ! fraction of new ridge area going to n2
      ,  fvol               ! fraction of new ridge volume going to n2

      real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi,nilyr) :: &    
         eirdg             ! ice energy of ridging ice

      real (kind=dbl_kind) :: &    
         hi1             &  ! thickness of ridging ice
      ,  hexp            &  ! ridge e-folding thickness
      ,  hL, hR          &  ! left and right limits of integration
      ,  expL, expR        ! exponentials involving hL, hR

      character (len=char_len) :: &
         fieldid           ! field identifier

      logical (kind=log_kind) :: &
         neg_aice0       &  ! flag for aice0(i,j) < -puny
      ,  large_ardg        ! flag for ardg > aicen_init

      !-----------------------------------------------------------------
      ! Compute quantities that ridging should conserve
      ! (not done for snow because snow may be dumped in ocean)
      !-----------------------------------------------------------------

      if (l_conservation_check) then
         call column_sum (ncat,   vicen, vice_init)
         call column_sum (ntilay, eicen, eice_init)
      endif

      !-----------------------------------------------------------------
      ! Compute change in open water area due to closing and opening.
      !-----------------------------------------------------------------

      neg_aice0 = .false.

      do j = jlo, jhi
      do i = ilo, ihi
         aice0(i,j) = aice0(i,j)                               &
                    - apartic(i,j,0)*closing_gross(i,j)*dyn_dt &
                    + opning(i,j)*dyn_dt                       
 
         if (aice0(i,j) < -puny) then                           
            neg_aice0 = .true.                                  
            aice0(i,j)= c0i
         elseif (aice0(i,j) < c0i) then    ! roundoff error      
            aice0(i,j) = c0i                                     
         endif                                                  
      enddo                                                     
      enddo                                                     

      IF(.FALSE.)THEN
      if (neg_aice0) then       ! there is a bug                
         do j = jlo, jhi                                        
         do i = ilo, ihi                                        
            if (aice0(i,j) < -puny) then                        
               write (nu_diag,*) ' '                            
               write (nu_diag,*) 'Ridging error: aice0 < 0'     
               write (nu_diag,*) 'istep1, my_task, i, j:',     &
                                  istep1, my_task, i, j
               !                 istep1, my_task, i, ngid(j),aice(i,j)
               write (nu_diag,*) 'aice0:', aice0(i,j)
!               call abort_ice('ridging: aice0 must be >= 0')
            endif               ! aice0 < -puny
         enddo                  ! i
         enddo                  ! j
      endif                     ! neg_aice0
      ENDIF

      !-----------------------------------------------------------------
      ! Save initial state variables
      !-----------------------------------------------------------------

      aicen_init(:,:,:) = aicen(:,:,:)
      vicen_init(:,:,:) = vicen(:,:,:)
      vsnon_init(:,:,:) = vsnon(:,:,:)
      aicen_init(:,:,:) = aicen(:,:,:)
      esnon_init(:,:,:) = esnon(:,:,:)
      eicen_init(:,:,:) = eicen(:,:,:)
            
      !-----------------------------------------------------------------
      ! Compute the area, volume, and energy of ice ridging in each
      !  category, along with the area of the resulting ridge.
      !-----------------------------------------------------------------

      do n1 = 1, ncat

      !-----------------------------------------------------------------
      ! Identify grid cells with nonzero ridging
      !-----------------------------------------------------------------

         icells = 0
         do j = jlo, jhi
         do i = ilo, ihi
            if (aicen_init(i,j,n1) > puny .and. apartic(i,j,n1) > c0i  &
                 .and. closing_gross(i,j) > c0i) then
               icells = icells + 1
               indxi(icells) = i
               indxj(icells) = j
            endif
         enddo                  ! i
         enddo                  ! j

         large_ardg = .false.

         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)

      !-----------------------------------------------------------------
      ! Compute area of ridging ice (ardg1) and of new ridge (ardg2).
      ! Make sure ardg1 <= aiceninit.
      !-----------------------------------------------------------------

            ardg1(i,j) = apartic(i,j,n1)*closing_gross(i,j)*dyn_dt

            if (ardg1(i,j) > aicen_init(i,j,n1) + puny) then
               large_ardg = .true.
            else           ! correct for roundoff error
               ardg1(i,j) = min (ardg1(i,j),aicen_init(i,j,n1))
            endif

            ardg2(i,j) = ardg1(i,j) / krdg(i,j,n1)
            afrac(i,j) = ardg1(i,j) / aicen_init(i,j,n1)

      !-----------------------------------------------------------------
      ! Subtract area, volume, and energy from ridging category n1.
      ! (Ice energy in separate loop for vector friendliness)
      !-----------------------------------------------------------------

            virdg(i,j) = vicen_init(i,j,n1) * afrac(i,j)
            vsrdg(i,j) = vsnon_init(i,j,n1) * afrac(i,j)
            esrdg(i,j) = esnon_init(i,j,n1) * afrac(i,j)

            aicen(i,j,n1) = aicen(i,j,n1) - ardg1(i,j)
            vicen(i,j,n1) = vicen(i,j,n1) - virdg(i,j)
            vsnon(i,j,n1) = vsnon(i,j,n1) - vsrdg(i,j)
            esnon(i,j,n1) = esnon(i,j,n1) - esrdg(i,j)

      !-----------------------------------------------------------------
      ! Increment ridging diagnostics
      !-----------------------------------------------------------------

            dardg1dt(i,j) = dardg1dt(i,j) + ardg1(i,j)
            dardg2dt(i,j) = dardg2dt(i,j) + ardg2(i,j)
            dvirdgdt(i,j) = dvirdgdt(i,j) + virdg(i,j)
            opening(i,j)  = opening (i,j) + opning(i,j)*dyn_dt

      !-----------------------------------------------------------------
      !  Place part of the snow lost by ridging into the ocean. 
      !  Note that esnow_mlt < 0; the ocean must cool to melt snow.
      !  If the ocean temp = Tf already, new ice must grow.
      !-----------------------------------------------------------------
               
            msnow_mlt(i,j) = msnow_mlt(i,j)                        &
                           + rhos*vsrdg(i,j)*(c1i-fsnowrdg)          
            esnow_mlt(i,j) = esnow_mlt(i,j)                        &
                           + esrdg(i,j)*(c1i-fsnowrdg)               
                                                                    
      !-----------------------------------------------------------------
      ! Compute quantities used to apportion ice among categories  
      ! in the n2 loop below                                       
      !-----------------------------------------------------------------
                                                                    
            dhr(i,j)  = hrmax(i,j,n1) - hrmin(i,j,n1)               
            dhr2(i,j) = hrmax(i,j,n1) * hrmax(i,j,n1)              &
                      - hrmin(i,j,n1) * hrmin(i,j,n1)

         enddo                  ! ij

      !-----------------------------------------------------------------
      ! Subtract ice energy from ridging category n1. 
      !-----------------------------------------------------------------

         do k = 1, nilyr
            do ij = 1, icells
               i = indxi(ij)
               j = indxj(ij)
               
               eirdg(i,j,k) = eicen_init(i,j,ilyr1(n1)+k-1) * afrac(i,j)
               eicen(i,j,ilyr1(n1)+k-1) = eicen(i,j,ilyr1(n1)+k-1)  &   
                                        - eirdg(i,j,k)                  
            enddo                                                       
         enddo                                                          
                                                                        
         if (large_ardg) then  ! there is a bug                         
            do ij = 1, icells                                           
               i = indxi(ij)                                            
               j = indxj(ij)                                            
               if (ardg1(i,j) > aicen_init(i,j,n1) + puny) then         
                  write (nu_diag,*) ''                                  
                  write (nu_diag,*) 'ardg > aicen'                     
                  write (nu_diag,*) 'istep1, my_task, i, j, n:',      & 
                                     istep1, my_task, i, j, n1          
                  write (nu_diag,*) 'ardg, aicen_init:',             &
                                     ardg1(i,j), aicen_init(i,j,n1)
!                  call abort_ice ('ridging: ardg must be <= aicen')
               endif            ! ardg1 > aicen_init
            enddo               ! ij
         endif                  ! large_ardg

      !-----------------------------------------------------------------
      ! Add area, volume, and energy of new ridge to each category n2.
      !-----------------------------------------------------------------

         do n2 = 1, ncat

            if (krdg_redist == 0) then ! Hibler 1980 formulation

               do ij = 1, icells
                  i = indxi(ij)
                  j = indxj(ij)

      !-----------------------------------------------------------------
      ! Compute the fraction of ridged ice area and volume going to 
      !  thickness category n2.
      !-----------------------------------------------------------------
               
                  if (hrmin(i,j,n1) >= hin_max(n2) .or.  & 
                      hrmax(i,j,n1) <= hin_max(n2-1)) then
                     hL = c0i
                     hR = c0i
                  else
                     hL = max (hrmin(i,j,n1), hin_max(n2-1))
                     hR = min (hrmax(i,j,n1), hin_max(n2))
                  endif
                  
                  ! fraction of ridged ice area and volume going to n2
                  farea(i,j) = (hR-hL) / dhr(i,j) 
                  fvol (i,j) = (hR*hR - hL*hL) / dhr2(i,j)
                  
               enddo            ! ij

            else         ! krdg_redist = 1; exponential formulation

      !-----------------------------------------------------------------
      ! Compute the fraction of ridged ice area and volume going to
      !  thickness category n2.
      !-----------------------------------------------------------------

               if (n2 < ncat) then

                  do ij = 1, icells
                     i = indxi(ij)
                     j = indxj(ij)

                     hi1  = hrmin(i,j,n1)
                     hexp = hrexp(i,j,n1)

                     if (hi1 >= hin_max(n2)) then
                        farea(i,j) = c0i
                        fvol (i,j) = c0i
                     else
                        hL = max (hi1, hin_max(n2-1))
                        hR = hin_max(n2)
                        expL = exp(-(hL-hi1)/hexp)
                        expR = exp(-(hR-hi1)/hexp)
                        farea(i,j) = expL - expR
                        fvol (i,j) = ((hL + hexp)*expL         &
                                    - (hR + hexp)*expR) / (hi1 + hexp)
                     endif
                  enddo         ! ij

               else             ! n2 = ncat

                  do ij = 1, icells
                     i = indxi(ij)
                     j = indxj(ij)

                     hi1  = hrmin(i,j,n1)
                     hexp = hrexp(i,j,n1)

                     hL = max (hi1, hin_max(n2-1))
                     expL = exp(-(hL-hi1)/hexp)
                     farea(i,j) = expL
                     fvol (i,j) = (hL + hexp)*expL / (hi1 + hexp)

                  enddo

               endif            ! n2 < ncat

            endif               ! krdg_redist

      !-----------------------------------------------------------------
      ! Transfer ice area, ice and snow volume, and snow energy to 
      ! category n2.
      !-----------------------------------------------------------------

            do ij = 1, icells
               i = indxi(ij)
               j = indxj(ij)
               aicen(i,j,n2) = aicen(i,j,n2) + farea(i,j)*ardg2(i,j)
               vicen(i,j,n2) = vicen(i,j,n2) + fvol(i,j) *virdg(i,j)
               vsnon(i,j,n2) = vsnon(i,j,n2)                    &
                             + fvol(i,j)*vsrdg(i,j)*fsnowrdg     
               esnon(i,j,n2) = esnon(i,j,n2)                    &
                             + fvol(i,j)*esrdg(i,j)*fsnowrdg

            enddo
            
      !-----------------------------------------------------------------
      ! Transfer ice energy to category n2
      !-----------------------------------------------------------------
            do k = 1, nilyr
               do ij = 1, icells
                  i = indxi(ij)
                  j = indxj(ij)
                  eicen(i,j,ilyr1(n2)+k-1) = eicen(i,j,ilyr1(n2)+k-1)     &
                                           + fvol(i,j)*eirdg(i,j,k)        
               enddo            ! ij                                       
            enddo               ! k                                        
                                                                           
         enddo                  ! n2 (new ridges)                          
      enddo                     ! n1 (ridging categories)                  
                                                                           
      !-----------------------------------------------------------------   
      ! Check volume and energy conservation                               
      !-----------------------------------------------------------------   
                                                                           
      if (l_conservation_check) then                                       
                                                                           
         call column_sum (ncat,   vicen, vice_final)                       
         fieldid = 'vice, ridging'                                         
         call column_conservation_check (vice_init, vice_final,           &
                                         puny,      fieldid)               
                                                                           
         call column_sum (ntilay, eicen, eice_final)                       
         fieldid = 'eice, ridging'                                         
         call column_conservation_check (eice_init, eice_final,           &
                                         puny*Lfresh*rhoi, fieldid)

      endif

      end subroutine ridge_shift

!=======================================================================
!BOP
!
! !ROUTINE: asum_ridging - find total fractional area
!
! !DESCRIPTION:
!
! Find the total area of ice plus open water in each grid cell.
!
! This is similar to the aggregate_area subroutine except that the
! total area can be greater than 1, so the open water area is 
! included in the sum instead of being computed as a residual. 
!
! !REVISION HISTORY:
!
! author William H. Lipscomb, LANL
!
! !INTERFACE:
!
      subroutine asum_ridging
!
! !USES:
! 
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) :: i, j, ni

      !-----------------------------------------------------------------
      ! open water
      !-----------------------------------------------------------------

      do j = jlo, jhi
      do i = ilo, ihi
         asum(i,j) = aice0(i,j)
      enddo
      enddo

      !-----------------------------------------------------------------
      ! ice categories
      !-----------------------------------------------------------------

      do ni = 1, ncat
         do j = jlo, jhi
         do i = ilo, ihi
            asum(i,j) = asum(i,j) + aicen(i,j,ni)
         enddo                  ! i
         enddo                  ! j
      enddo                     ! ni 

      end subroutine asum_ridging
      
!=======================================================================

      end module ice_mechred

!=======================================================================
