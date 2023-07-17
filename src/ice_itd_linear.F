!/===========================================================================/
! CVS VERSION INFORMATION
! $Id$
! $Name$
! $Revision$
!/===========================================================================/

!=======================================================================
!BOP
!
! !MODULE: ice_itd_linear - linear remapping scheme for ITD
!
! !DESCRIPTION:
!
! Linear remapping scheme for the ice thickness distribution
!
! See Lipscomb, W. H.  Remapping the thickness distribution of sea  \\
!     ice. 2001, J. Geophys. Res., Vol 106, 13989--14000.
!
! !REVISION HISTORY:
!
! authors: William H. Lipscomb, LANL
!          Elizabeth C. Hunke, LANL
!
! Vectorized by Clifford Chen (Fujitsu) and William Lipscomb (LANL)
!
! !INTERFACE:
!
      module ice_itd_linear
!
! !USES:
!
      use ice_model_size
      use ice_kinds_mod
      use ice_domain
      use ice_constants
      use ice_state
      use ice_itd
      use ice_calendar
      use ice_fileunits
!
!EOP
!

      implicit none

!=======================================================================

      contains

!=======================================================================
!BOP
!
! !IROUTINE: linear_itd - ITD scheme that shifts ice among categories
!
! !INTERFACE:
!
      subroutine linear_itd (hicen_old, hicen)
!
! !DESCRIPTION:
!
! Ice thickness distribution scheme that shifts ice among categories. \\
!
! The default scheme is linear remapping, which works as follows.  See
! Lipscomb (2001) for more details. \\
!
!   Using the thermodynamic "velocities", interpolate to find the 
!   velocities in thickness space at the category boundaries, and 
!   compute the new locations of the boundaries.  Then for each 
!   category, compute the thickness distribution function,  g(h), 
!   between hL and hR, the left and right boundaries of the category.
!   Assume g(h) is a linear polynomial that satisfies two conditions: \\
!
!   (1) The ice area implied by g(h) equals aicen(n).
!   (2) The ice volume implied by g(h) equals aicen(n)*hicen(n).
!
! Given g(h), at each boundary compute the ice area and volume lying 
! between the original and new boundary locations.  Transfer area 
! and volume across each boundary in the appropriate direction, thus
! restoring the original boundaries. See Lipscomb (2001) for details. 
!
! !REVISION HISTORY:
!
! authors: William H. Lipscomb, LANL
!          Elizabeth C. Hunke, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind), dimension(ilo:ihi,jlo:jhi,ncat), &
           intent(in) :: &
         hicen_old        ! starting value of hicen (m)

      real (kind=dbl_kind), dimension(ilo:ihi,jlo:jhi,ncat), &
           intent(inout) :: &
         hicen            ! ice thickness for each cat        (m)
!
!EOP
!
      integer (kind=int_kind) :: & 
         i, j           &  ! horizontal indices
      ,  ni, nd         &   ! category indices
      ,  k                ! ice layer index

      real (kind=dbl_kind) :: &
         slope            &! rate of change of dhice with hice
      ,  dh0              &! change in ice thickness at h = 0
      ,  da0              &! area melting from category 1
      ,  damax            &! max allowed reduction in category 1 area
      ,  etamin, etamax   &! left and right limits of integration 
      ,  x1               &! etamax - etamin
      ,  x2               &! (etamax^2 - etamin^2) / 2
      ,  x3               &! (etamax^3 - etamin^3) / 3
      ,  wk1, wk2          ! temporary variables

      real (kind=dbl_kind), dimension(ilo:ihi,jlo:jhi,0:ncat) :: &
         Hbnew            ! new boundary locations

      real (kind=dbl_kind), dimension(ilo:ihi,jlo:jhi) :: &
         Hb0            &  ! hin_max(0)
      ,  Hb1              ! hin_max(1)

      real (kind=dbl_kind), dimension(ilo:ihi,jlo:jhi,ncat) :: &
         dhicen          & ! thickness change for remapping (m)
      ,  g0              & ! constant coefficient in g(h)
      ,  g1              & ! linear coefficient in g(h) 
      ,  hL              & ! left end of range over which g(h) > 0
      ,  hR               ! right end of range over which g(h) > 0

      real (kind=dbl_kind), dimension(imt_local, jmt_local) :: &
         vice_init, vice_final  &! ice volume summed over categories
      ,  vsno_init, vsno_final  &! snow volume summed over categories
      ,  eice_init, eice_final  &! ice energy summed over categories
      ,  esno_init, esno_final  ! snow energy summed over categories

      ! NOTE: Third index of donor, daice, dvice should be ncat-1,
      !       except that compilers would have trouble when ncat = 1 
      integer (kind=int_kind), dimension(ilo:ihi,jlo:jhi,ncat) :: & 
         donor            ! donor category index

      real (kind=dbl_kind), dimension(ilo:ihi,jlo:jhi,ncat) :: &
         daice          &  ! ice area transferred across boundary
      ,  dvice            ! ice volume transferred across boundary

      logical (kind=log_kind), dimension(ilo:ihi,jlo:jhi) :: &
         remap_flag       ! remap ITD if remap_flag(i,j) is true

      character (len=char_len) :: &
         fieldid           ! field identifier

      logical (kind=log_kind), parameter :: &
         l_conservation_check = .true.   ! if true, check conservation 
                                         ! (useful for debugging)

       integer (kind=int_kind) :: &
         icells         &   ! number of grid cells with ice
      ,  ij                ! combined horizontal index

       integer (kind=int_kind), dimension (1:(ihi-ilo+1)*(jhi-jlo+1)) :: &
        indxi            &  ! compressed i/j indices
      , indxj

      logical (kind=log_kind) :: &
         flag_changed

!!  ggao 60162008
      real (kind=dbl_kind) ::EPS
!! change end


      !-----------------------------------------------------------------
      ! Compute volume and energy sums that linear remapping should 
      !  conserve.
      !-----------------------------------------------------------------

      if (l_conservation_check) then
         call column_sum (ncat,   vicen, vice_init)
         call column_sum (ncat,   vsnon, vsno_init)
         call column_sum (ntilay, eicen, eice_init)
         call column_sum (ncat,   esnon, esno_init)
      endif

      !-----------------------------------------------------------------
      ! Compute thickness change in each category. 
      !-----------------------------------------------------------------

      dhicen = c0i
      do ni = 1, ncat
         do j = jlo,jhi
         do i = ilo,ihi
            if (aicen(i,j,ni) > puny) then
               dhicen(i,j,ni) = hicen(i,j,ni) - hicen_old(i,j,ni)
            endif               ! aicen > puny
         enddo                  ! i
         enddo                  ! j
      enddo                     ! n

      !-----------------------------------------------------------------
      ! Compute fractional ice area in each grid cell.
      !-----------------------------------------------------------------
      call aggregate_area

      !-----------------------------------------------------------------
      ! Identify grid cells with ice and initialize remapping flag.
      ! Remapping is done wherever remap_flag = .true.
      ! In rare cases the category boundaries may shift too far for the 
      !  remapping algorithm to work, and remap_flag is set to .false.  
      ! In these cases the simpler 'rebin' subroutine will shift ice 
      !  between categories if needed.
      !-----------------------------------------------------------------

      icells = 0
      do j = jlo,jhi
      do i = ilo,ihi
         if (aice(i,j) > puny) then
            remap_flag(i,j) = .true.
            icells = icells + 1
            indxi(icells) = i
            indxj(icells) = j
         else
            remap_flag(i,j) = .false.
         endif
      enddo
      enddo

      !-----------------------------------------------------------------
      ! Compute new category boundaries, Hbnew, based on changes in
      ! ice thickness from vertical thermodynamics.
      !-----------------------------------------------------------------

      hin_max(ncat) = 999.9_dbl_kind ! arbitrary big number
      Hbnew = c0i

      do ni = 1, ncat-1 

         do ij = 1, icells       ! aice(i,j) > puny
            i = indxi(ij)
            j = indxj(ij)

            if (hicen_old(i,j,ni)   > puny .and.         &
                hicen_old(i,j,ni+1) > puny) then
                 ! interpolate between adjacent category growth rates
               slope = (dhicen(i,j,ni+1)-dhicen(i,j,ni)) /  & 
!                       (hicen_old(i,j,ni+1)-hicen_old(i,j,ni))
!!  ggao change 0616-2008
                       (hicen_old(i,j,ni+1)-hicen_old(i,j,ni)+EPSILON(EPS))
               Hbnew(i,j,ni) = hin_max(ni) + dhicen(i,j,ni)  &
                            + slope * (hin_max(ni) - hicen_old(i,j,ni))
            elseif (hicen_old(i,j,ni) > puny) then ! hicen_old(ni+1)=0
               Hbnew(i,j,ni) = hin_max(ni) + dhicen(i,j,ni)
            elseif (hicen_old(i,j,ni+1) > puny) then ! hicen_old(ni)=0
               Hbnew(i,j,ni) = hin_max(ni) + dhicen(i,j,ni+1)
            else
               Hbnew(i,j,ni) = hin_max(ni)
            endif
         enddo                  ! ij

      !-----------------------------------------------------------------
      ! Check that each boundary lies between adjacent values of hicen.
      ! If not, set remap_flag = .false.
      !-----------------------------------------------------------------

         flag_changed = .false.
         do ij = 1, icells       ! aice(i,j) > puny
            i = indxi(ij)
            j = indxj(ij)

            if (aicen(i,j,ni) > puny .and.            &
                hicen(i,j,ni) >= Hbnew(i,j,ni)) then
               remap_flag(i,j) = .false.
               flag_changed = .true.
            elseif (aicen(i,j,ni+1) > puny .and.      &
                    hicen(i,j,ni+1) <= Hbnew(i,j,ni)) then
               remap_flag(i,j) = .false.
               flag_changed = .true.
            endif

      !-----------------------------------------------------------------
      ! Check that Hbnew(n) lies between hin_max(n-1) and hin_max(n+1).
      ! If not, set remap_flag = .false.
      ! (In principle we could allow this, but it would make the code 
      ! more complicated.)
      !-----------------------------------------------------------------

            if (Hbnew(i,j,ni) > hin_max(ni+1)) then
               remap_flag(i,j) = .false.
               flag_changed = .true.
            endif

            if (Hbnew(i,j,ni) < hin_max(ni-1)) then
               remap_flag(i,j) = .false.
               flag_changed = .true.
            endif

         enddo                  ! ij

      !-----------------------------------------------------------------
      ! Write diagnosis outputs if remap_flag was changed to false
      !-----------------------------------------------------------------

         if (flag_changed) then
            do ij = 1, icells   ! aice(i,j) > puny

               i = indxi(ij)
               j = indxj(ij)

               if (aicen(i,j,ni) > puny .and.                  &
                   hicen(i,j,ni) >= Hbnew(i,j,ni)) then         
                  write(nu_diag,*) 'istep1 = ',istep1           
                  write(nu_diag,*) my_task,':',i,j,            &
                       'ITD: hicen(ni) > Hbnew(ni)'             
                  write(nu_diag,*) 'cat ',ni                    
                  write(nu_diag,*) my_task,':',i,j,            &
                       'hicen(ni) =', hicen(i,j,ni)             
                  write(nu_diag,*) my_task,':',i,j,            &
                       'Hbnew(ni) =', Hbnew(i,j,ni)             
               elseif (aicen(i,j,ni+1) > puny .and.            &
                       hicen(i,j,ni+1) <= Hbnew(i,j,ni)) then   
                  write(nu_diag,*) 'istep1 = ',istep1           
                  write(nu_diag,*) my_task,':',i,j,            &
                       'ITD: hicen(ni+1) < Hbnew(ni)'           
                  write(nu_diag,*) 'cat ',ni                     
                  write(nu_diag,*) my_task,':',i,j,            &
                       'hicen(ni+1) =', hicen(i,j,ni+1)         
                  write(nu_diag,*) my_task,':',i,j,            &
                       'Hbnew(ni) =', Hbnew(i,j,ni)             
               endif                                            
                                                                
               if (Hbnew(i,j,ni) > hin_max(ni+1)) then            
                  write(nu_diag,*) 'istep1 = ',istep1           
                  write(nu_diag,*) my_task,':',i,j,            &
                       'ITD Hbnew(ni) > hin_max(ni+1)'
                  write(nu_diag,*) 'cat ',ni
                  write(nu_diag,*) my_task,':',i,j,           &
                       'Hbnew(ni) =', Hbnew(i,j,ni)              
                  write(nu_diag,*) my_task,':',i,j,           &
                       'hin_max(ni+1) =', hin_max(ni+1)          
               endif                                           
                                                               
               if (Hbnew(i,j,ni) < hin_max(ni-1)) then           
                  write(nu_diag,*) 'istep1 = ',istep1          
                  write(nu_diag,*) my_task,':',i,j,           &
                       'ITD: Hbnew(ni) < hin_max(ni-1)'          
                  write(nu_diag,*) 'cat ',ni                    
                  write(nu_diag,*) my_task,':',i,j,           &
                       'Hbnew(ni) =', Hbnew(i,j,ni)              
                  write(nu_diag,*) my_task,':',i,j,           &
                       'hin_max(ni-1) =', hin_max(ni-1)
               endif

            enddo               ! ij
         endif                  ! flag_changed

      enddo                     ! boundaries, 1 to ncat-1

      !-----------------------------------------------------------------
      ! Identify cells where the ITD is to be remapped
      !-----------------------------------------------------------------

      icells = 0
      do j = jlo,jhi
      do i = ilo,ihi
         if (remap_flag(i,j)) then
            icells = icells + 1
            indxi(icells) = i
            indxj(icells) = j
         endif
      enddo
      enddo

      !-----------------------------------------------------------------
      ! Fill arrays with initial boundaries of category 1
      ! Prescribe Hbnew(0) and Hbnew(ncat)
      !-----------------------------------------------------------------

      do j = jlo,jhi
      do i = ilo,ihi
         Hb0(i,j) = hin_max(0)
         Hb1(i,j) = hin_max(1)

         Hbnew(i,j,0) = c0i

         if (aicen(i,j,ncat) > puny) then
           Hbnew(i,j,ncat) = c3i*hicen(i,j,ncat) - c2i*Hbnew(i,j,ncat-1)
         else
            Hbnew(i,j,ncat) = hin_max(ncat)
         endif

         if (Hbnew(i,j,ncat) < hin_max(ncat-1))   &
              Hbnew(i,j,ncat) = hin_max(ncat-1)
      enddo
      enddo

      !-----------------------------------------------------------------
      ! Compute g(h) for category 1 at start of time step 
      ! (hicen = hicen_old)
      !-----------------------------------------------------------------

      call fit_line(1,         Hb0,       Hb1,       hicen_old(:,:,1), &
                    g0(:,:,1), g1(:,:,1), hL(:,:,1), hR(:,:,1),        &
                    remap_flag)
 
      !-----------------------------------------------------------------
      ! Find area lost due to melting of thin (category 1) ice
      !-----------------------------------------------------------------

      do ij = 1, icells    ! remap_flag = .true.
         i = indxi(ij)
         j = indxj(ij)

         if (aicen(i,j,1) > puny) then

            dh0 = dhicen(i,j,1)         

            if (dh0 < c0i) then   ! remove area from category 1   
               dh0 = min(-dh0,hin_max(1))   ! dh0 --> |dh0| 
	
      !-----------------------------------------------------------------
      ! Integrate g(1) from 0 to dh0 to estimate area melted
      !-----------------------------------------------------------------

               ! right integration limit (left limit = 0)
               etamax = min(dh0,hR(i,j,1)) - hL(i,j,1) 
                                              
               if (etamax > c0i) then
                  x1 = etamax
                  x2 = p5 * etamax*etamax
                  da0 = g1(i,j,1)*x2 + g0(i,j,1)*x1 ! ice area removed

               ! constrain new thickness <= hicen_old
                  damax = aicen(i,j,1)                       &
                        * (c1i-hicen(i,j,1)/hicen_old(i,j,1)) ! damax > 0
                  da0 = min (da0, damax)
      
               ! remove area, conserving volume
                  hicen(i,j,1) = hicen(i,j,1)                    &
!                               * aicen(i,j,1) / (aicen(i,j,1)-da0)
!!  gao change 0616-2008
                               * aicen(i,j,1) / (aicen(i,j,1)-da0+EPSILON(EPS))
                  aicen(i,j,1) = aicen(i,j,1) - da0
               endif            ! etamax > 0

            else                ! dh0 >= 0
               Hbnew(i,j,0) = min(dh0,hin_max(1))  ! shift Hbnew(0) to right
            endif

         endif                  ! aicen(i,j,1) > puny
      enddo                     ! ij

      !-----------------------------------------------------------------
      ! Compute g(h) for each ice thickness category.
      !-----------------------------------------------------------------

      do ni = 1, ncat
        call fit_line(ni, Hbnew(:,:,ni-1), Hbnew(:,:,ni), hicen(:,:,ni),& 
                      g0(:,:,ni), g1(:,:,ni), hL(:,:,ni),    hR(:,:,ni),& 
                      remap_flag)
      enddo

      !-----------------------------------------------------------------
      ! Compute area and volume to be shifted across each boundary.
      !-----------------------------------------------------------------

      donor(:,:,:) = 0
      daice(:,:,:) = c0i
      dvice(:,:,:) = c0i

      do ni = 1, ncat-1
         do ij = 1, icells   ! remap_flag = .true.
            i = indxi(ij)
            j = indxj(ij)

            if (Hbnew(i,j,ni) > hin_max(ni)) then ! transfer from n to n+1

               ! left and right integration limits in eta space
               etamin = max(hin_max(ni), hL(i,j,ni)) - hL(i,j,ni)
               etamax = min(Hbnew(i,j,ni), hR(i,j,ni)) - hL(i,j,ni) 
               donor(i,j,ni) = ni

            else             ! Hbnew(n) <= hin_max(n); transfer from n+1 to n

               ! left and right integration limits in eta space
               etamin = c0i                                 
               etamax = min(hin_max(ni), hR(i,j,ni+1)) - hL(i,j,ni+1)
               donor(i,j,ni) = ni+1

            endif            ! Hbnew(n) > hin_max(n)

            if (etamax > etamin) then
               x1  = etamax - etamin
               wk1 = etamin*etamin
               wk2 = etamax*etamax
               x2  = p5 * (wk2 - wk1)
               wk1 = wk1*etamin
               wk2 = wk2*etamax
               x3  = p333 * (wk2 - wk1)
               nd  = donor(i,j,ni)
               daice(i,j,ni) = g1(i,j,nd)*x2 + g0(i,j,nd)*x1        
               if (daice(i,j,ni) > c0i) then
                  dvice(i,j,ni) = g1(i,j,nd)*x3 + g0(i,j,nd)*x2 &
                               + daice(i,j,ni)*hL(i,j,nd)
               else
                  daice(i,j,ni) = c0i
                  donor(i,j,ni) = 0
               endif
            endif


            ! If daice or dvice is very small, shift no ice.

            nd = donor(i,j,ni)

!!!===============================
            IF(ND>0) THEN
!!!===============================

            if (daice(i,j,ni) < aicen(i,j,nd)*puny) then
               daice(i,j,ni) = c0i
               dvice(i,j,ni) = c0i
               donor(i,j,ni) = 0
            endif 

            if (dvice(i,j,ni) < vicen(i,j,nd)*puny) then
               daice(i,j,ni) = c0i
               dvice(i,j,ni) = c0i
               donor(i,j,ni) = 0
            endif

            ! If daice is close to aicen or dvice is close to vicen,
            ! shift entire category

            if (daice(i,j,ni) > aicen(i,j,nd)*(c1i-puny)) then
               daice(i,j,ni) = aicen(i,j,nd)
               dvice(i,j,ni) = vicen(i,j,nd)
            endif

            if (dvice(i,j,ni) > vicen(i,j,nd)*(c1i-puny)) then
               daice(i,j,ni) = aicen(i,j,nd)
               dvice(i,j,ni) = vicen(i,j,nd)
            endif
!!!   ggao fix the bug 02102008
            ELSE

               daice(i,j,ni) = c0i
               dvice(i,j,ni) = c0i
               donor(i,j,ni) = 0
            ENDIF
!!! change end

         enddo                  ! ij
      enddo                     ! boundaries, 1 to ncat-1

      !-----------------------------------------------------------------
      ! Shift ice between categories as necessary  
      !-----------------------------------------------------------------

      call shift_ice (donor, daice, dvice, hicen)

      !-----------------------------------------------------------------
      ! Make sure hice(i,j,1) >= minimum ice thickness hi_min.
      !-----------------------------------------------------------------

      do ij = 1, icells          ! remap_flag = .true.
         i = indxi(ij)
         j = indxj(ij)
         if (hi_min > c0i .and.       &
              aicen(i,j,1) > puny .and. hicen(i,j,1) < hi_min) then
            aicen(i,j,1) = aicen(i,j,1) * hicen(i,j,1)/hi_min
            hicen(i,j,1) = hi_min
         endif
      enddo                     ! ij

      !-----------------------------------------------------------------
      ! Update ice and open water area.
      !-----------------------------------------------------------------
      call aggregate_area

      !-----------------------------------------------------------------
      ! Check volume and energy conservation.
      !-----------------------------------------------------------------

      if (l_conservation_check) then

         call column_sum (ncat,   vicen, vice_final)
         fieldid = 'vice, ITD remap'
         call column_conservation_check (vice_init, vice_final,      &
                                         puny,      fieldid)          
                                                                      
         call column_sum (ncat,   vsnon, vsno_final)                  
         fieldid = 'vsno, ITD remap'                                  
         call column_conservation_check (vsno_init, vsno_final,      &
                                         puny,      fieldid)          
                                                                      
         call column_sum (ntilay, eicen, eice_final)                  
         fieldid = 'eice, ITD remap'                                  
         call column_conservation_check (eice_init, eice_final,      &
                                         puny*Lfresh*rhoi, fieldid)   
                                                                      
         call column_sum (ncat,   esnon, esno_final)                  
         fieldid = 'esno, ITD remap'                                  
         call column_conservation_check (esno_init, esno_final,      &
                                         puny*Lfresh*rhos, fieldid)

      endif                     ! conservation check

      end subroutine linear_itd

!=======================================================================
!BOP
!
! !IROUTINE: fit_line - fit g(h) with a line using area, volume constraints
!
! !INTERFACE:
!
      subroutine fit_line (ni,  HbL, HbR, hice,  &
                          g0, g1,  hL,  hR, remap_flag)
!
! !DESCRIPTION:
!
! Fit g(h) with a line, satisfying area and volume constraints.
! To reduce roundoff errors caused by large values of g0 and g1,
! we actually compute g(eta), where eta = h - hL, and hL is the
! left boundary.
!
! !REVISION HISTORY:
!
! authors: William H. Lipscomb, LANL
!          Elizabeth C. Hunke, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: ni      ! category index

      real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi), &
           intent(in) :: & 
         HbL, HbR        ! left and right category boundaries

      real (kind=dbl_kind), dimension(ilo:ihi,jlo:jhi), intent(in) :: & 
         hice            ! ice thickness

      real (kind=dbl_kind), dimension(ilo:ihi,jlo:jhi), intent(out) :: &
         g0, g1       &   ! coefficients in linear equation for g(eta)
      ,  hL           &   ! min value of range over which g(h) > 0
      ,  hR              ! max value of range over which g(h) > 0

      logical (kind=log_kind), dimension (ilo:ihi,jlo:jhi),      &
         intent(in) :: &
         remap_flag
!
!EOP
!
      integer (kind=int_kind) :: &
         i,j             ! horizontal indices

      real  (kind=dbl_kind) :: &
         h13           &  ! HbL + 1/3 * (HbR - HbL)
      ,  h23           &  ! HbL + 2/3 * (HbR - HbL)
      ,  dhr           &  ! 1 / (hR - hL)
      ,  wk1, wk2        ! temporary variables

!!  ggao 60162008
      real (kind=dbl_kind) ::EPS
!! change end


      do j = jlo,jhi
      do i = ilo,ihi

         if (remap_flag(i,j) .and. aicen(i,j,ni) > puny    &
                      .and. HbR(i,j) - HbL(i,j) > puny) then

         ! Initialize hL and hR

            hL(i,j) = HbL(i,j)
            hR(i,j) = HbR(i,j)

         ! Change hL or hR if hicen(n) falls outside central third of range

            h13 = p333 * (c2i*hL(i,j) + hR(i,j))
            h23 = p333 * (hL(i,j) + c2i*hR(i,j))
            if (hice(i,j) < h13) then
               hR(i,j) = c3i*hice(i,j) - c2i*hL(i,j)
            elseif (hice(i,j) > h23) then
               hL(i,j) = c3i*hice(i,j) - c2i*hR(i,j)
            endif

         ! Compute coefficients of g(eta) = g0 + g1*eta

!            dhr = c1i / (hR(i,j) - hL(i,j))
!!  ggao change 0616-2008
            dhr = c1i / (hR(i,j) - hL(i,j)+EPSILON(EPS))
!! change end
            wk1 = c6i * aicen(i,j,ni) * dhr
            wk2 = (hice(i,j) - hL(i,j)) * dhr
            g0(i,j) = wk1 * (p666 - wk2)
            g1(i,j) = c2i*dhr * wk1 * (wk2 - p5)
               
         else                   ! remap_flag = .false. or aicen < puny
                                ! or hbR <= hbL
            hL(i,j) = c0i
            hR(i,j) = c0i
            g0(i,j) = c0i
            g1(i,j) = c0i
            
         endif                  ! aicen > puny

      enddo                     ! i
      enddo                     ! j

      end subroutine fit_line

!=======================================================================

      end module ice_itd_linear

!=======================================================================
