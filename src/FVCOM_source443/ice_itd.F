!/===========================================================================/
! CVS VERSION INFORMATION
! $Id$
! $Name$
! $Revision$
!/===========================================================================/

!=======================================================================
!BOP
!
! !MODULE: ice_itd - initialize and redistribute ice in the ITD
!
! !DESCRIPTION:
!
! Routines to initialize the ice thickness distribution and 
! utilities to redistribute ice among categories. These routines 
! are not specific to a particular numerical implementation. \\
!
! See Bitz, C.M., and W.H. Lipscomb, 1999: 
! An energy-conserving thermodynamic model of sea ice,
! J. Geophys. Res., 104, 15,669--15,677. \\
!     
! See Bitz, C.M., M.M. Holland, A.J. Weaver, M. Eby, 2001: 
! Simulating the ice-thickness distribution in a climate model,
! J. Geophys. Res., 106, 2441--2464. \\
!
! !REVISION HISTORY:
!
! author: C. M. Bitz, UW
!         Elizabeth C. Hunke, LANL
!         William H. Lipscomb, LANL
!
! Summer 2003: Vectorized by Clifford Chen (Fujitsu) and William Lipscomb (LANL)
!
! !INTERFACE:
!
      module ice_itd
!
! !USES:
!
      use ice_kinds_mod
      use ice_model_size
      use ice_constants
      use ice_state
      use ice_fileunits
!      use ice_exit

!
!EOP
!
      implicit none
      save

      integer (kind=int_kind) :: &
         kitd              & ! type of itd conversions 
                             !   0 = delta function
                             !   1 = linear remap
      ,  kcatbound         & !   0 = old category boundary formula
                             !   1 = new formula giving round numbers
      ,  ilyr1 (ncat)      & ! position of the top layer in each cat
      ,  ilyrn (ncat)        ! position of the bottom layer in each cat

      real (kind=dbl_kind), parameter :: &
         hi_min = p01       ! minimum ice thickness allowed (m)

      real (kind=dbl_kind) :: &
         hin_max(0:ncat)    ! category limits                (m)

!-------------------------------------------------------------------
! a note regarding hi_min and hin_max(0):
! both represent a minimum ice thickness.  hin_max(0) is
! intended to be used for particular numerical implementations
! of category conversions in the ice thickness distribution.
! hi_min is a more general purpose parameter, but is specifically 
! for maintaining stability in the thermodynamics.  Currently,
! hi_min = 0.1 m
! hin_max(0) = 0.1 m for the delta function itd
! hin_max(0) = 0.0 m for linear remapping
!
! similarly, there are two values of minimum snow thickness
! (the other is defined in ice_vthermo.H since it is used only
! for thermo.) 
!
! Also note that the upper limit on the thickest category
! is only used for the linear remapping scheme
! and it is not a true upper limit on the thickness
!-------------------------------------------------------------------

!=======================================================================

      contains

!=======================================================================
!BOP
!
! !IROUTINE: init_itd - initalize area fraction and thickness boundaries for ITD
!
! !INTERFACE:
!
      subroutine init_itd
!
! !DESCRIPTION:
!
! Initialize area fraction and thickness boundaries for the itd model
!
! !REVISION HISTORY:
!
! authors: William H. Lipscomb, LANL
!          Elizabeth C. Hunke LANL
!          C. M. Bitz UW
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) :: &
           ni    ! thickness category index

      real (kind=dbl_kind) :: &
           cc1, cc2, cc3       &        ! parameters for kcatbound = 0
      ,    x1

      real (kind=dbl_kind) :: &
           rn        &! real(n)
      ,    rncat     &! real(ncat)
      ,    d1        &                  ! parameters for kcatbound = 1 (m)
      ,    d2

      if (ncat == 1 .and. kitd == 1) then
         write (nu_diag,*) 'Remapping the ITD is not allowed for ncat=1'
         write (nu_diag,*) 'Use the delta function ITD option instead'
!         call abort_ice ('(init_itd) Linear remapping not allowed')
      endif

      rncat = real(ncat, kind=dbl_kind)
#  if defined (ICE_FRESHWATER)
! ayumi 20151112 & EJA 20160920 - thickness categories for the Great Lakes
! (5 categories)
!       d1 = 1.0_dbl_kind / rncat
!       d2 = 0.5_dbl_kind / rncat
      d1 = 0.25_dbl_kind / rncat
      d2 = 0.75_dbl_kind / rncat
#  else
      d1 = 3.0_dbl_kind / rncat
      d2 = 0.5_dbl_kind / rncat
#  endif

      !-----------------------------------------------------------------
      ! Choose category boundaries based on one of two formulas.
      !
      ! The first formula (kcatbound = 0) was used in Lipscomb (2001) 
      !  and in CICE versions 3.0 and 3.1.
      !
      ! The second formula is more user-friendly in the sense that it
      !  is easy to obtain round numbers for category boundaries:
      !
      !    H(n) = n * [d1 + d2*(n-1)] 
      ! 
      ! Default values are d1 = 300/ncat, d2 = 50/ncat.
      ! For ncat = 5, boundaries in cm are 60, 140, 240, 360, which are 
      !  close to the standard values given by the first formula.
      ! For ncat = 10, boundaries in cm are 30, 70, 120, 180, 250, 330,
      !  420, 520, 630.    
      !-----------------------------------------------------------------

      if (kcatbound == 0) then   ! original scheme

         if (kitd == 1) then
            ! linear remapping itd category limits
            cc1 = c3i/rncat
            cc2 = c15*cc1
            cc3 = c3i

            hin_max(0) = c0i     ! minimum ice thickness, m
         else
            ! delta function itd category limits
            cc1 = max(1.1_dbl_kind/rncat,c1i*hi_min)
            cc2 = c25*cc1
            cc3 = 2.25_dbl_kind

            ! hin_max(0) should not be zero
            ! use some caution in making it less than 0.10
            hin_max(0) = hi_min ! minimum ice thickness, m
         endif                  ! kitd

         do ni = 1, ncat
            x1 = real(ni-1,kind=dbl_kind) / rncat
            hin_max(ni) = hin_max(ni-1)   &
                       + cc1 + cc2*(c1i + tanh(cc3*(x1-c1i)))
         enddo

      elseif (kcatbound == 1) then  ! new scheme

         hin_max(0) = c0i
         do ni = 1, ncat
            rn = real(ni, kind=dbl_kind)
            hin_max(ni) = rn * (d1 + (rn-c1i)*d2)
         enddo

      endif

      if (my_task == master_task) then
         write (nu_diag,*) ' '
         write (nu_diag,*) 'The thickness categories are:'
         write (nu_diag,*) ' '
         write (nu_diag,*) 'hin_max(n-1) < Cat n < hin_max(n)'
         do ni = 1, ncat
         write (nu_diag,*) hin_max(ni-1),' < Cat ',ni, ' < ',hin_max(ni)
         enddo
         write (nu_diag,*) ' '
      endif

      !-----------------------------------------------------------------
      ! vectors identifying first and last layer in each category
      !-----------------------------------------------------------------
      ilyr1(1) = 1                       ! if nilyr  = 4
      ilyrn(1) = nilyr                   !   ilyr1 = { 1,5,9 }
      do ni = 2,ncat                      !   ilyrn = { 4,8,12} etc
         ilyr1(ni) = ilyrn(ni-1) + 1
         ilyrn(ni) = ilyrn(ni-1) + nilyr
      enddo

      end subroutine init_itd

!=======================================================================
!BOP
!
! !IROUTINE: aggregate - aggregate ice state variables
!
! !INTERFACE:
!
      subroutine aggregate
!
! !DESCRIPTION:
!
! Aggregate ice state variables over thickness categories.
!
! !REVISION HISTORY:
!
! authors: C. M. Bitz, UW
!          W. H. Lipscomb, LANL
!
! !USES:
!
      use ice_domain
      use ice_flux, only : Tf 
      use ice_grid
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) :: i, j, k, ni 

      integer (kind=int_kind), dimension (1:(ihi-ilo+1)*(jhi-jlo+1)) :: &
        indxi              &    ! compressed indices in i/j directions
      , indxj

      integer (kind=int_kind) :: &
        icells                 & ! number of ocean/ice cells
      , ij                    ! combined i/j horizontal index

!!  ggao 60162008
      real (kind=dbl_kind) ::EPS
!! change end

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------
      do j=jlo,jhi
      do i=ilo,ihi
         aice0(i,j) = c1i
         aice(i,j) = c0i
         vice(i,j) = c0i
         vsno(i,j) = c0i
         eice(i,j) = c0i
         esno(i,j) = c0i
         Tsfc(i,j) = c0i
      enddo
      enddo

      !-----------------------------------------------------------------
      ! Aggregate
      !-----------------------------------------------------------------

      icells = 0
      do j = jlo, jhi
      do i = ilo, ihi
        if (tmask(i,j)) then
          icells = icells + 1
          indxi(icells) = i
          indxj(icells) = j
        endif   ! tmask
      enddo
      enddo

      do ni = 1, ncat

         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)
            aice(i,j) = aice(i,j) + aicen(i,j,ni)
            vice(i,j) = vice(i,j) + vicen(i,j,ni)
            vsno(i,j) = vsno(i,j) + vsnon(i,j,ni)
            esno(i,j) = esno(i,j) + esnon(i,j,ni)
            Tsfc(i,j) = Tsfc(i,j) + Tsfcn(i,j,ni)*aicen(i,j,ni)
         enddo                  ! ij

         do k = 1, nilyr
            do ij = 1, icells
               i = indxi(ij)
               j = indxj(ij)
               eice(i,j) = eice(i,j) + eicen(i,j,(ni-1)*nilyr+k)
            enddo
         enddo                  ! k

      enddo                     ! n

      !-----------------------------------------------------------------
      ! Temperature, open water fraction
      !-----------------------------------------------------------------
      do j=jlo,jhi
      do i=ilo,ihi
         if (aice(i,j) > c0i)  then
            aice0(i,j) = max (c1i - aice(i,j), c0i)
!            Tsfc(i,j)  = Tsfc(i,j) / aice(i,j)
!!  ggao change 0616-2008
            Tsfc(i,j)  = Tsfc(i,j) /(aice(i,j)+EPSILON(EPS))
         else
            Tsfc(i,j)  = Tf(i,j)
         endif
      enddo                     ! i
      enddo                     ! j

      end subroutine aggregate

!=======================================================================
!BOP
!
! !IROUTINE: aggregate_area - aggregate ice area
!
! !INTERFACE:
!
      subroutine aggregate_area
!
! !DESCRIPTION:
!
! Aggregate ice area (but not other state variables) over thickness categories
!
! !REVISION HISTORY:
!
! authors: William H. Lipscomb, LANL
!          modified Jan 2004 by Clifford Chen, Fujitsu
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!

      use ice_grid, only : tlat,tlon
!     ggao
      integer (kind=int_kind) :: i, j, ni

      logical (kind=log_kind) :: &
         outbound          ! = .true. if aggregate ice area out of bounds


      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------
      do j=jlo,jhi
      do i=ilo,ihi
         aice(i,j) = c0i
      enddo
      enddo

      !-----------------------------------------------------------------
      ! Aggregate
      !-----------------------------------------------------------------
      do ni = 1, ncat
         do j=jlo,jhi
         do i=ilo,ihi
            aice(i,j) = aice(i,j) + aicen(i,j,ni)
         enddo                  ! i
         enddo                  ! j
      enddo                     ! n

      outbound = .false.

      do j = jlo, jhi
      do i = ilo, ihi

      !-----------------------------------------------------------------
      ! Bug check
      !-----------------------------------------------------------------
         if (aice(i,j) > c1i+puny .or.  &
             aice(i,j) < -puny) then
            outbound = .true.
         endif

      !-----------------------------------------------------------------
      ! open water fraction
      !-----------------------------------------------------------------
         aice0(i,j) = max (c1i - aice(i,j), c0i)

      enddo                     ! i
      enddo                     ! j

      if ( outbound ) then      ! area out of bounds
        do j = jlo, jhi
        do i = ilo, ihi

         if (aice(i,j) > c1i+puny .or.  &
             aice(i,j) < -puny) then
!            write(nu_diag,*) ' '
!            write(nu_diag,*) 'aggregate ice area out of bounds'
!            write(nu_diag,*) 'my_task, i, j, aice:',     &
!                              my_task, i, j, aice(i,j)

!!            write(nu_diag,*) 'aicen =', aicen(i,j,:) 
!            call abort_ice('ice_itd: aggregate_area')
         endif

        enddo                   ! i
        enddo                   ! j
      endif                     ! outbound

      end subroutine aggregate_area

!=======================================================================
!BOP
!
! !IROUTINE: bound_aggregate - bound calls for aggregate ice state
!
! !INTERFACE:
!
      subroutine bound_aggregate
!
! !DESCRIPTION:
!
! Get ghost cell values for aggregate ice state variables
! NOTE: This subroutine is called only at initialization.  It could be
!       eliminated if the aggregate variables were defined only in
!       physical grid cells.
!
! !REVISION HISTORY:
!
! authors: William H. Lipscomb, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      use ice_grid

!      call bound (aice0)
!      call bound (aice)
!      call bound (vice)
!      call bound (vsno)
!      call bound (eice)
!      call bound (esno)
!      call bound (Tsfc)

!    need inter processors exchange 
!   ggao


      end subroutine bound_aggregate

!=======================================================================
!BOP
!
! !IROUTINE: bound_state - bound calls for ice state variables
!
! !INTERFACE:
!
      subroutine bound_state
!
! !DESCRIPTION:
!
! Get ghost cell values for ice state variables in each thickness category
!
! !REVISION HISTORY:
!
! authors: William H. Lipscomb, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      use ice_grid

!      call bound_narr (ncat,aicen)
!      call bound_narr (ncat,vicen)
!      call bound_narr (ncat,vsnon)
!      call bound_narr (ntilay,eicen)
!      call bound_narr (ncat,esnon)
!      call bound_narr (ncat,Tsfcn)
!    need inter processors exchange 
!   ggao


      end subroutine bound_state

!=======================================================================
!BOP
!
! !IROUTINE: check_state - require certain fields to be monotone
!
! !INTERFACE:
!
      subroutine check_state
!
! !DESCRIPTION:
!
!  Insist that certain fields are monotone.
!  Should not be necessary if all is well, 
!  but best to keep going. Model will not conserve
!  energy and water if fields are zeroed here.
!
! !REVISION HISTORY:
!
! author: C. M. Bitz, UW
!
! !USES:
!
      use ice_flux
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) :: &
         i, j         &   ! horizontal indices
      ,  ni           &    ! thickness category index
      ,  k               ! ice layer index

      integer (kind=int_kind), dimension (ilo:ihi,jlo:jhi) :: &
         zerout    ! 0=false, 1=true

      do ni = 1, ncat

         do j = jlo, jhi
         do i = ilo, ihi
            if (aicen(i,j,ni) < puny .or. vicen(i,j,ni) < puny) then
               zerout(i,j) = 1
            else
               zerout(i,j) = 0
            endif
         enddo
         enddo

         do k = 1, nilyr
            do j = jlo, jhi
            do i = ilo, ihi
               if (eicen(i,j,ilyr1(ni)+k-1) > -puny) zerout(i,j) = 1
            enddo               ! i
            enddo               ! j
         enddo                  ! k

         do j = jlo, jhi
         do i = ilo, ihi

            if (zerout(i,j)==1) then
               aice0(i,j) = aice0(i,j) + aicen(i,j,ni) 
               aicen(i,j,ni) = c0i
               vicen(i,j,ni) = c0i
               vsnon(i,j,ni) = c0i
               esnon(i,j,ni) = c0i
               Tsfcn(i,j,ni) = Tf(i,j)
            elseif (vsnon(i,j,ni) <= puny) then
               vsnon(i,j,ni) = c0i
               esnon(i,j,ni) = c0i
            endif

            if (vsnon(i,j,ni) > puny) then
               if (-esnon(i,j,ni)/vsnon(i,j,ni)-Lfresh*rhos < eps04)  &
                    esnon(i,j,ni) = -vsnon(i,j,ni)*(Lfresh*rhos + eps04)
            endif

         enddo                  ! i
         enddo                  ! j

         do k = 1, nilyr
            do j = jlo, jhi
            do i = ilo, ihi
               if (zerout(i,j)==1) eicen(i,j,ilyr1(ni)+k-1) = c0i
            enddo
            enddo
         enddo                  ! k

      enddo                     ! n

      do j=jlo,jhi
      do i=ilo,ihi
         if (aice0(i,j) < puny) aice0(i,j) = c0i
      enddo
      enddo

      end subroutine check_state

!=======================================================================
!BOP
!
! !IROUTINE: rebin - rebins thicknesses into defined categories
!
! !INTERFACE:
!
      subroutine rebin
!
! !DESCRIPTION:
!
! Rebins thicknesses into defined categories
!
! !REVISION HISTORY:
!
! authors: William H. Lipscomb, LANL
!          Elizabeth C. Hunke LANL 
!
! !USES:
!
      use ice_grid
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) :: &
         i,j            &    ! horizontal indices
      ,  ni                  ! category index

      logical (kind=log_kind) :: &
         shiftflag          ! = .true. if ice must be shifted

      real (kind=dbl_kind), dimension(ilo:ihi,jlo:jhi,ncat) :: &
         hicen              ! ice thickness for each cat        (m)

      integer (kind=int_kind), dimension(ilo:ihi,jlo:jhi,ncat) :: &
         donor              ! donor category index

      real (kind=dbl_kind), dimension(ilo:ihi,jlo:jhi,ncat) :: &
         daice            & ! ice area transferred
      ,  dvice              ! ice volume transferred

      !-----------------------------------------------------------------
      ! Compute ice thickness.
      !-----------------------------------------------------------------
      do ni = 1, ncat
         do j = jlo, jhi
         do i = ilo, ihi
            if (aicen(i,j,ni) > puny) then
               hicen(i,j,ni) = vicen(i,j,ni) / aicen(i,j,ni)
            else
               hicen(i,j,ni) = c0i
            endif
         enddo                  ! i
         enddo                  ! j
      enddo                     ! n

      !-----------------------------------------------------------------
      ! make sure thickness of cat 1 is at least hin_max(0)
      !-----------------------------------------------------------------
      do j = jlo, jhi
      do i = ilo, ihi
         if (aicen(i,j,1) > puny) then
            if (hicen(i,j,1) <= hin_max(0) .and. hin_max(0) > c0i ) then
               aicen(i,j,1) = vicen(i,j,1) / hin_max(0)
               hicen(i,j,1) = hin_max(0)
            endif
         endif
      enddo                     ! i
      enddo                     ! j

      !-----------------------------------------------------------------
      ! If a category thickness is not in bounds, shift the
      ! entire area, volume, and energy to the neighboring category
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! Initialize shift arrays
      !-----------------------------------------------------------------

      donor(:,:,:) = 0
      daice(:,:,:) = c0i
      dvice(:,:,:) = c0i

      !-----------------------------------------------------------------
      ! Move thin categories up
      !-----------------------------------------------------------------

      do ni = 1, ncat-1          ! loop over category boundaries

      !-----------------------------------------------------------------
      ! identify thicknesses that are too big
      !-----------------------------------------------------------------
         shiftflag = .false.
         do j = jlo, jhi
         do i = ilo, ihi
            if (aicen(i,j,ni) > puny .and.  & 
                hicen(i,j,ni) > hin_max(ni)) then
               shiftflag = .true.
               donor(i,j,ni) = ni
               daice(i,j,ni) = aicen(i,j,ni)
               dvice(i,j,ni) = vicen(i,j,ni)
            endif
         enddo                  ! i
         enddo                  ! j

         if (shiftflag) then

      !-----------------------------------------------------------------
      ! shift ice between categories
      !-----------------------------------------------------------------
            call shift_ice (donor, daice, dvice, hicen)
             
      !-----------------------------------------------------------------
      ! reset shift parameters
      !-----------------------------------------------------------------
            do j = jlo, jhi
            do i = ilo, ihi
               donor(i,j,ni) = 0
               daice(i,j,ni) = c0i
               dvice(i,j,ni) = c0i
            enddo
            enddo

         endif                  ! shiftflag

      enddo                     ! n

      !-----------------------------------------------------------------
      ! Move thick categories down
      !-----------------------------------------------------------------

      do ni = ncat-1, 1, -1      ! loop over category boundaries

      !-----------------------------------------------------------------
      ! identify thicknesses that are too small
      !-----------------------------------------------------------------
         shiftflag = .false.
         do j = jlo, jhi
         do i = ilo, ihi
            if (aicen(i,j,ni+1) > puny .and.  &
                hicen(i,j,ni+1) <= hin_max(ni)) then
               shiftflag = .true.
               donor(i,j,ni) = ni+1
               daice(i,j,ni) = aicen(i,j,ni+1)
               dvice(i,j,ni) = vicen(i,j,ni+1)
            endif
         enddo                  ! i
         enddo                  ! j

         if (shiftflag) then

      !-----------------------------------------------------------------
      ! shift ice between categories
      !-----------------------------------------------------------------
            call shift_ice (donor, daice, dvice, hicen)

      !-----------------------------------------------------------------
      ! reset shift parameters
      !-----------------------------------------------------------------
            do j = jlo, jhi
            do i = ilo, ihi
               donor(i,j,ni) = 0
               daice(i,j,ni) = c0i
               dvice(i,j,ni) = c0i
            enddo
            enddo

         endif                  ! shiftflag

      enddo                     ! n


      end subroutine rebin

!=======================================================================
!BOP
!
! !IROUTINE: reduce_area - reduce area when ice melts for special case ncat=1
!
! !INTERFACE:
!
      subroutine reduce_area(hice1_old, hice1)
!
! !DESCRIPTION:
!
! Reduce area when ice melts for special case of ncat=1
!
! Use CSM 1.0-like method of reducing ice area
! when melting occurs: assume only half the ice volume
! change goes to thickness decrease, the other half
! to reduction in ice fraction
!
! !REVISION HISTORY:
!
! authors: C. M. Bitz, UW 
! modified by: Elizabeth C. Hunke, LANL
!
! !USES:
!
      use ice_grid
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind), dimension(ilo:ihi,jlo:jhi), &
           intent(in) :: &
         hice1_old   ! old ice thickness for category 1 (m)

      real (kind=dbl_kind), dimension(ilo:ihi,jlo:jhi), &
           intent(inout) :: &
         hice1       ! new ice thickness for category 1 (m)
!
!EOP
!
      integer (kind=int_kind) :: &
         i, j        ! horizontal indices

      real (kind=dbl_kind) :: &
         hi0       & ! current hi for ice fraction adjustment
      ,  dai0      & ! change in aice for ice fraction adjustment
      ,  dhi         ! hice1 - hice1_old

      do j=jlo,jhi
      do i=ilo,ihi
         if (tmask(i,j)) then

      !-----------------------------------------------------------------
      ! make sure thickness of cat 1 is at least hin_max(0)
      !-----------------------------------------------------------------

            if (hice1(i,j) <= hin_max(0) .and. hin_max(0) > c0i ) then
               aicen(i,j,1) = vicen(i,j,1) / hin_max(0)
               hice1(i,j) = hin_max(0)
            endif

            if (aicen(i,j,1) > c0i) then
               dhi = hice1(i,j) - hice1_old(i,j)
               if (dhi < c0i) then  
                  hi0  = vicen(i,j,1) / aicen(i,j,1)
                  dai0 = vicen(i,j,1) / (hi0-p5*dhi) &
                       - aicen(i,j,1)
                  aicen(i,j,1) = aicen(i,j,1) + dai0 
               endif
            endif

         endif                  ! tmask
      enddo                     ! i
      enddo                     ! j

      end subroutine reduce_area

!=======================================================================
!BOP
!
! !IROUTINE: shift_ice - shift ice across category boundaries 
!
! !INTERFACE:
!
      subroutine shift_ice (donor, daice, dvice, hicen)
!
! !DESCRIPTION:
!
! Shift ice across category boundaries, conserving area, volume, and
! energy.
!
! !REVISION HISTORY:
!
! authors: William H. Lipscomb, LANL
!          Elizabeth C. Hunke, LANL
!
! !USES:
!
      use ice_flux
      use ice_work, only: worka
!
! !INPUT/OUTPUT PARAMETERS:
! 
      ! NOTE: Third index of donor, daice, dvice should be ncat-1,
      !       except that compilers would have trouble when ncat = 1 
      integer (kind=int_kind), dimension(ilo:ihi,jlo:jhi,ncat), &
         intent(in) :: &
         donor             ! donor category index

      real (kind=dbl_kind), dimension(ilo:ihi,jlo:jhi,ncat), &
           intent(inout) :: &
         daice          &  ! ice area transferred across boundary
      ,  dvice             ! ice volume transferred across boundary

      real (kind=dbl_kind), dimension(ilo:ihi,jlo:jhi,ncat), &
           intent(inout) :: &
         hicen             ! ice thickness for each cat        (m)
!
!EOP
!
      integer (kind=int_kind) :: &
         i, j            & ! horizontal indices
      ,  ni              & ! thickness category index
      ,  n2              & ! receiver category
      ,  n1              & ! donor category
      ,  k                 ! ice layer index

      real (kind=dbl_kind), dimension(ilo:ihi,jlo:jhi,ncat) :: &
         aTsfn

      real (kind=dbl_kind) :: &
         dvsnow          & ! snow volume transferred
      ,  desnow          & ! snow energy transferred
      ,  deice           & ! ice energy transferred
      ,  daTsf             ! aicen*Tsfcn transferred

      integer (kind=int_kind), dimension (1:(ihi-ilo+1)*(jhi-jlo+1)) :: &
        indxi         & ! compressed indices for i/j directions
      , indxj

      integer (kind=int_kind) :: &
        icells        & ! number of cells with ice to transfer
      , ij              ! combined i/j horizontal index

      logical (kind=log_kind) :: &
        daice_negative      & ! true if daice < -puny
      , dvice_negative      & ! true if dvice < -puny
      , daice_greater_aicen & ! true if daice > aicen
      , dvice_greater_vicen  ! true if dvice > vicen

!!  ggao 60162008
      real (kind=dbl_kind) ::EPS
!! change end


      !-----------------------------------------------------------------
      ! Define a variable equal to aicen*Tsfcn
      !-----------------------------------------------------------------
      do ni = 1, ncat
         do j = jlo,jhi
         do i = ilo,ihi
            aTsfn(i,j,ni) = aicen(i,j,ni)*Tsfcn(i,j,ni)
         enddo                  ! i
         enddo                  ! j
      enddo                     ! n

      !-----------------------------------------------------------------
      ! Check for daice or dvice out of range, allowing for roundoff error
      !-----------------------------------------------------------------

      do ni = 1, ncat-1

         daice_negative = .false.
         dvice_negative = .false.
         daice_greater_aicen = .false.
         dvice_greater_vicen = .false.

         do j = jlo,jhi
         do i = ilo,ihi

            if (donor(i,j,ni) > 0) then 
               n1 = donor(i,j,ni)

               if (daice(i,j,ni) < c0i) then
                  if (daice(i,j,ni) > -puny*aicen(i,j,n1)) then   
                     daice(i,j,ni) = c0i ! shift no ice
                     dvice(i,j,ni) = c0i
                  else
                     daice_negative = .true.
                  endif
               endif
         
               if (dvice(i,j,ni) < c0i) then
                  if (dvice(i,j,ni) > -puny*vicen(i,j,n1)) then   
                     daice(i,j,ni) = c0i ! shift no ice
                     dvice(i,j,ni) = c0i
                  else
                     dvice_negative = .true.
                  endif
               endif

               if (daice(i,j,ni) > aicen(i,j,n1)*(c1i-puny)) then
                  if (daice(i,j,ni) < aicen(i,j,n1)*(c1i+puny)) then
                     daice(i,j,ni) = aicen(i,j,n1)
                     dvice(i,j,ni) = vicen(i,j,n1)
                  else
                     daice_greater_aicen = .true.
                  endif
               endif    

               if (dvice(i,j,ni) > vicen(i,j,n1)*(c1i-puny)) then
                  if (dvice(i,j,ni) < vicen(i,j,n1)*(c1i+puny)) then
                     daice(i,j,ni) = aicen(i,j,n1)
                     dvice(i,j,ni) = vicen(i,j,n1)
                  else
                     dvice_greater_vicen = .true.
                  endif
               endif
               
            endif               ! donor > 0 
         enddo                  ! i
         enddo                  ! j

      !-----------------------------------------------------------------
      ! error messages
      !-----------------------------------------------------------------

         if (daice_negative) then
            do j = jlo,jhi
            do i = ilo,ihi
               if (donor(i,j,ni) > 0 .and.                            &
                   daice(i,j,ni) <= -puny*aicen(i,j,n1)) then          
                  write(nu_diag,*) my_task,':',i,j,                   &
                       'ITD Neg daice =',daice(i,j,ni),' boundary',ni   
!                  call abort_ice ('ice: ITD Neg daice')                
               endif                                                   
            enddo                                                      
            enddo                                                      
         endif                                                         
                                                                       
         if (dvice_negative) then                                      
            do j = jlo,jhi                                             
            do i = ilo,ihi                                             
               if (donor(i,j,ni) > 0 .and.                            &
                   dvice(i,j,ni) <= -puny*vicen(i,j,n1)) then          
                  write(nu_diag,*) my_task,':',i,j,                   &
                       'ITD Neg dvice =',dvice(i,j,ni),' boundary',ni    
!                  call abort_ice ('ice: ITD Neg dvice')                
               endif                                                   
            enddo                                                      
            enddo                                                      
         endif                                                         
                                                                       
         if (daice_greater_aicen) then                                 
            do j = jlo,jhi                                             
            do i = ilo,ihi                                             
               if (donor(i,j,ni) > 0) then                              
                  n1 = donor(i,j,ni)                                    
                  if (daice(i,j,ni) >= aicen(i,j,n1)*(c1i+puny)) then    
                     write(nu_diag,*) my_task,':',i,j,                &
                          'ITD daice > aicen, cat',n1                  
                     write(nu_diag,*) my_task,':',i,j,                &
                          'daice =', daice(i,j,ni),                    &
                          'aicen =', aicen(i,j,n1)
!                     call abort_ice ('ice: ITD daice > aicen')
                  endif
               endif
            enddo
            enddo
         endif

         if (dvice_greater_vicen) then
            do j = jlo,jhi
            do i = ilo,ihi
               if (donor(i,j,ni) > 0) then
                  n1 = donor(i,j,ni)
                  if (dvice(i,j,ni) >= vicen(i,j,n1)*(c1i+puny)) then
                     write(nu_diag,*) my_task,':',i,j,          &
                          'ITD dvice > vicen, cat',n1            
                     write(nu_diag,*) my_task,':',i,j,          &
                          'dvice =', dvice(i,j,ni),              &
                          'vicen =', vicen(i,j,n1)
!                     call abort_ice ('ice: ITD dvice > vicen')
                  endif
               endif
            enddo
            enddo
         endif

      !-----------------------------------------------------------------
      ! transfer volume and energy between categories
      !-----------------------------------------------------------------

         icells = 0
         do j = jlo, jhi
         do i = ilo, ihi
           if (daice(i,j,ni) > c0i) then ! daice(n) can be < puny
             icells = icells + 1
             indxi(icells) = i
             indxj(icells) = j
           endif   ! tmask
         enddo
         enddo

         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)

            n1 = donor(i,j,ni)
           !worka(i,j) = dvice(i,j,ni) / vicen(i,j,n1)
!! ggao change
           worka(i,j) = dvice(i,j,ni) / (vicen(i,j,n1)+EPSILON(EPS))
!! ggao change end 0616-2008

            if (n1  ==  ni) then
               n2 = n1+1
            else                ! n1 = n+1
               n2 = ni
            endif
            
            aicen(i,j,n1) = aicen(i,j,n1) - daice(i,j,ni)
            aicen(i,j,n2) = aicen(i,j,n2) + daice(i,j,ni)
            vicen(i,j,n1) = vicen(i,j,n1) - dvice(i,j,ni)
            vicen(i,j,n2) = vicen(i,j,n2) + dvice(i,j,ni)
            
            dvsnow = vsnon(i,j,n1) * worka(i,j)
            vsnon(i,j,n1) = vsnon(i,j,n1) - dvsnow
            vsnon(i,j,n2) = vsnon(i,j,n2) + dvsnow
            
            daTsf = daice(i,j,ni)*Tsfcn(i,j,n1)
            aTsfn(i,j,n1) = aTsfn(i,j,n1) - daTsf
            aTsfn(i,j,n2) = aTsfn(i,j,n2) + daTsf 

            desnow = esnon(i,j,n1) * worka(i,j)
            esnon(i,j,n1) = esnon(i,j,n1) - desnow
            esnon(i,j,n2) = esnon(i,j,n2) + desnow

         enddo                  ! ij

         do k = 1, nilyr
            do ij = 1, icells
               i = indxi(ij)
               j = indxj(ij)

               n1 = donor(i,j,ni)
               if (n1  ==  ni) then
                  n2 = n1+1
               else             ! n1 = n+1
                  n2 = ni
               endif

               deice = eicen(i,j,ilyr1(n1)+k-1) * worka(i,j)
               eicen(i,j,ilyr1(n1)+k-1) =           &
                    eicen(i,j,ilyr1(n1)+k-1) - deice
               eicen(i,j,ilyr1(n2)+k-1) =           &
                    eicen(i,j,ilyr1(n2)+k-1) + deice
            enddo               ! ij
         enddo                  ! k

      enddo                     ! boundaries, 1 to ncat-1

      !-----------------------------------------------------------------
      ! Update ice thickness and temperature
      !-----------------------------------------------------------------

      do ni = 1, ncat
         do j = jlo,jhi
         do i = ilo,ihi
            if (aicen(i,j,ni) > puny) then
!               hicen(i,j,ni) = vicen(i,j,ni) / aicen(i,j,ni)
!               Tsfcn(i,j,ni) = aTsfn(i,j,ni) / aicen(i,j,ni)
!!  ggao change
               hicen(i,j,ni) = vicen(i,j,ni) /( aicen(i,j,ni)+EPSILON(EPS))
               Tsfcn(i,j,ni) = aTsfn(i,j,ni) /( aicen(i,j,ni)+EPSILON(EPS))
!! ggao change end 06162008

            else
               hicen(i,j,ni) = c0i
               Tsfcn(i,j,ni) = Tf(i,j)
            endif
         enddo                  ! i
         enddo                  ! j
      enddo                     ! n

      end subroutine shift_ice

!=======================================================================
!BOP
!
! !IROUTINE: column_sum - sum field over all ice categories
!
! !INTERFACE:
!
      subroutine column_sum (nsum, xin, xout)
!
! !DESCRIPTION:
!
! For each grid cell, sum field over all ice categories.
!
! !REVISION HISTORY:
!
! author: William H. Lipscomb, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
           nsum                              ! number of categories/layers

      real (kind=dbl_kind), intent(in) :: &
           xin  (imt_local,jmt_local,nsum)   ! input field


      real (kind=dbl_kind), intent(out) :: &
           xout (imt_local,jmt_local)        ! output field
!
!EOP
!
      integer (kind=int_kind) :: &
           i, j               &  ! horizontal indices
      ,    ni                    ! category/layer index

      do j = 1, jmt_local
      do i = 1, imt_local
         xout(i,j) = c0i
      enddo
      enddo

      do ni = 1, nsum
         do j = 1, jmt_local
         do i = 1, imt_local
            xout(i,j) = xout(i,j) + xin(i,j,ni)
         enddo                  ! i
         enddo                  ! j
      enddo                     ! n

      end subroutine column_sum

!=======================================================================
!BOP
!
! !IROUTINE: column_conservation_check
!
! !INTERFACE:
!
      subroutine column_conservation_check (x1, x2, max_err, fieldid)
!
! !DESCRIPTION:
!
! For each physical grid cell, check that initial and final values
! of a conserved field are equal to within a small value.
!
! !REVISION HISTORY:
!
! author: William H. Lipscomb, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind), intent(in) :: &
           x1 (imt_local,jmt_local)    ! initial field

      real (kind=dbl_kind), intent(in) :: &
           x2 (imt_local,jmt_local)    ! final field

      real (kind=dbl_kind), intent(in) :: & 
          max_err                     ! max allowed error

      character (len=char_len), intent(in) :: &
           fieldid                     ! field identifier
!
!EOP
!
      integer (kind=int_kind) :: &
           i, j                        ! horizontal indices      

      logical (kind=log_kind) :: &
         conserv_err          ! = .true. if conservation check failed

      conserv_err = .false.

      do j = jlo, jhi
      do i = ilo, ihi
         if (abs(x2(i,j) - x1(i,j)) > max_err) then
            conserv_err = .true.
         endif
      enddo
      enddo

      if ( conserv_err ) then
        do j = jlo, jhi
        do i = ilo, ihi
         if (abs(x2(i,j) - x1(i,j)) > max_err) then
            write (nu_diag,*) ' '
            write (nu_diag,*) 'Conservation error: ', fieldid
            write (nu_diag,*)  my_task, ':', i, j
            write (nu_diag,*) 'Initial value =', x1(i,j)
            write (nu_diag,*) 'Final value =',   x2(i,j)
            write (nu_diag,*) 'Difference =', x2(i,j) - x1(i,j)
!            call abort_ice ('ice: Conservation error')
         endif
        enddo
        enddo
      endif

      end subroutine column_conservation_check

!=======================================================================
!BOP
!
! !IROUTINE: zap_small_areas - eliminate very small ice areas
!
! !INTERFACE:
!
      subroutine zap_small_areas
!
! !DESCRIPTION:
!
! For each ice category in each grid cell, remove ice if the fractional
! area is less than puny.
!
! !REVISION HISTORY:
!
! author: William H. Lipscomb, LANL
! Nov 2003:  Modified by Julie Schramm to conserve volume and energy
! Sept 2004: Modified by William Lipscomb; replaced normalize_state with
!            additions to local freshwater, salt, and heat fluxes
!            
!
! !USES:
!
      use ice_flux
!      use ice_calendar, only: dt
      use ice_calendar, only: dtice

!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) :: &
           i,j      & ! horizontal indices
      ,    ni       &  ! ice category index
      ,    k        & ! ice layer index
      ,    icells   & ! number of cells with ice to zap
      ,    ij        ! combined i/j horizontal index

      integer (kind=int_kind), dimension (1:(ihi-ilo+1)*(jhi-jlo+1)) :: &
           indxi    & ! compressed indices for i/j directions
      ,    indxj

      real (kind=dbl_kind) :: &
           xtmp      ! temporary variable

      do ni = 1, ncat

      !-----------------------------------------------------------------
      ! Count categories to be zapped.
      ! Abort model in case of negative area.
      !-----------------------------------------------------------------

         icells = 0
         do j = jlo, jhi
         do i = ilo, ihi
            if (aicen(i,j,ni) < -puny) then
               write (nu_diag,*) 'Negative ice area: i, j, n:', i, j, ni
               write (nu_diag,*) 'aicen =', aicen(i,j,ni)
!               call abort_ice ('zap: negative ice area')
            elseif ((aicen(i,j,ni) >= -puny .and. aicen(i,j,ni) < c0i) &
                                         .or.                         &
                    (aicen(i,j,ni) > c0i .and. aicen(i,j,ni) <= puny)) then
               icells = icells + 1
               indxi(icells) = i
               indxj(icells) = j
            endif
         enddo
         enddo

      !-----------------------------------------------------------------
      ! Zap ice energy and use ocean heat to melt ice
      !-----------------------------------------------------------------

         do k = 1, nilyr
            do ij = 1, icells
               i = indxi(ij)
               j = indxj(ij)

!               xtmp = eicen(i,j,ilyr1(ni)+k-1) / dt ! < 0
               xtmp = eicen(i,j,ilyr1(ni)+k-1) / dtice ! < 0
               fhnet(i,j)      = fhnet(i,j)      + xtmp
               fhnet_hist(i,j) = fhnet_hist(i,j) + xtmp
               eicen(i,j,ilyr1(ni)+k-1) = c0i

            enddo               ! ij
         enddo                  ! k

         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)

      !-----------------------------------------------------------------
      ! Zap snow energy and use ocean heat to melt snow
      !-----------------------------------------------------------------

!            xtmp = esnon(i,j,ni) / dt ! < 0
            xtmp = esnon(i,j,ni) / dtice ! < 0
            fhnet(i,j)      = fhnet(i,j)      + xtmp
            fhnet_hist(i,j) = fhnet_hist(i,j) + xtmp
            esnon(i,j,ni) = c0i

      !-----------------------------------------------------------------
      ! zap ice and snow volume, add water and salt to ocean
      !-----------------------------------------------------------------

!            xtmp = (rhoi*vicen(i,j,ni) + rhos*vsnon(i,j,ni)) / dt
            xtmp = (rhoi*vicen(i,j,ni) + rhos*vsnon(i,j,ni)) / dtice
            fresh(i,j)      = fresh(i,j)      + xtmp
            fresh_hist(i,j) = fresh_hist(i,j) + xtmp

!            xtmp = rhoi*vicen(i,j,n)*ice_ref_salinity*p001/ dt
            xtmp = rhoi*vicen(i,j,ni)*ice_ref_salinity*p001/ dtice
            fsalt(i,j)      = fsalt(i,j)      + xtmp
            fsalt_hist(i,j) = fsalt_hist(i,j) + xtmp

            aice0(i,j) = aice0(i,j) + aicen(i,j,ni)
            aicen(i,j,ni) = c0i
            vicen(i,j,ni) = c0i
            vsnon(i,j,ni) = c0i
            Tsfcn(i,j,ni) = Tf(i,j)

         enddo                  ! ij
      enddo                     ! n

      end subroutine zap_small_areas

!=======================================================================

      end module ice_itd

!=======================================================================
