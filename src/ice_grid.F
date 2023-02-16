!/===========================================================================/
! CVS VERSION INFORMATION
! $Id$
! $Name$
! $Revision$
!/===========================================================================/

!=======================================================================
!BOP
!
! !MODULE: ice_grid - spatial grids, masks and boundary conditions
!
! !DESCRIPTION:
!
! Spatial grids, masks, and boundary conditions
!
! !REVISION HISTORY:
!
! authors: Elizabeth C. Hunke, LANL
!          Tony Craig, NCAR
!
! !INTERFACE:
!
      module ice_grid
!
! !USES:
!
      use ice_kinds_mod
      use ice_constants
      use ice_domain
      use ice_fileunits
!      use ice_mpi_internal
      use ice_work, only: work_g1, work_g2, work_l1, worka
!====================================================================
!     link to FVCOM
!      use LIMS
      USE MOD_UTILS
      use ALL_VARS , only: VX,VY,ART1,vxmin,vymin
      ! need the area of t-cell
!       ggao
!
!EOP
!
      implicit none

      character (len=char_len)  ::  & 
!         grid_file & !  input file for POP grid info
!      ,  kmt_file  & !  input file for POP grid info
        grid_type   !  rectangular (default) or displaced_pole 

!      real (kind=dbl_kind), dimension (imt_global,jmt_global)  ::  &
!         TLAT_G     &! latitude of cell center T grid
!      ,  TLON_G      ! longitude of cell center T grid

!      real (kind=dbl_kind), dimension (imt_local,jmt_local)  ::  &
!         dxt        &! width of T-cell through the middle (m)
!      ,  dyt        &! height of T-cell through the middle (m)
!      ,  dxu        &! width of U-cell through the middle (m)
!      ,  dyu        &! height of U-cell through the middle (m)
!      ,  HTE        &! length of eastern edge of T-cell (m)
!      ,  HTN        &! length of northern edge of T-cell (m)
!      ,  HTS        &! length of southern edge of T-cell
!      ,  HTW        &! length of western edge of T-cell
!      ,  tarea      &! area of T-cell (m^2)
!      ,  uarea      &! area of U-cell (m^2)
!      ,  ULON       &! longitude of velocity pts (radians)
!      ,  ULAT       &! latitude of velocity pts (radians)
!      ,  TLON       &! longitude of temp pts (radians)
!      ,  TLAT       &! latitude of temp pts (radians)
!      ,  dxhy       &! 0.5*(HTE - HTE)
!      ,  dyhx        ! 0.5*(HTN - HTN)

!      real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi)  ::  &
!         cyp        &! 1.5*HTE - 0.5*HTE
!      ,  cxp        &! 1.5*HTN - 0.5*HTN
!      ,  cym        &! 0.5*HTE - 1.5*HTE
!      ,  cxm        &! 0.5*HTN - 1.5*HTN
!      ,  dxt2       &! 0.5*dxt
!      ,  dyt2       &! 0.5*dyt
!      ,  dxt4       &! 0.25*dxt
!      ,  dyt4       &! 0.25*dyt
!      ,  tarear     &! 1/tarea
!      ,  uarear     &! 1/uarea
!      ,  tinyarea   &! puny*tarea
!      ,  ANGLE      &! for conversions between POP grid and lat/lon
!      ,  ANGLET     &! ANGLE converted to T-cells
!      ,  tarean     &! area of NH T-cells
!      ,  tareas      ! area of SH T-cells

      ! Masks
!      real (kind=dbl_kind), dimension (imt_local,jmt_local)  ::  &
!         hm         &! land/boundary mask, thickness (T-cell)
!      ,  uvm        &! land/boundary mask, velocity (U-cell)
!      ,  mask_n     &! northern hemisphere
!      ,  mask_s      ! southern hemisphere

!      logical (kind=log_kind)  ::  &
!         tmask(imt_local,jmt_local) &! land/boundary mask, thickness (T-cell)
!      ,  umask(imt_local,jmt_local) &! land/boundary mask, velocity (U-cell)
!      ,  icetmask(ilo:ihi,jlo:jhi)  &! ice extent mask (T-cell)
!      ,  iceumask(ilo:ihi,jlo:jhi)   ! ice extent mask (U-cell)

!      real (kind=dbl_kind)  ::  &
!         shlat = -40.0_dbl_kind     &! artificial masking edge
!      ,  nhlat =  35.0_dbl_kind      ! artificial masking edge

      real (kind=dbl_kind), dimension(:,:),allocatable   ::  &
         TLAT_G     &! latitude of cell center T grid
      ,  TLON_G      ! longitude of cell center T grid

      real (kind=dbl_kind), dimension (:,:),allocatable  ::  &
         dxt        &! width of T-cell through the middle (m)
      ,  dyt        &! height of T-cell through the middle (m)
      ,  dxu        &! width of U-cell through the middle (m)
      ,  dyu        &! height of U-cell through the middle (m)
      ,  HTE        &! length of eastern edge of T-cell (m)
      ,  HTN        &! length of northern edge of T-cell (m)
      ,  HTS        &! length of southern edge of T-cell
      ,  HTW        &! length of western edge of T-cell
      ,  tarea      &! area of T-cell (m^2)
      ,  uarea      &! area of U-cell (m^2)
      ,  ULON       &! longitude of velocity pts (radians)
      ,  ULAT       &! latitude of velocity pts (radians)
      ,  TLON       &! longitude of temp pts (radians)
      ,  TLAT       &! latitude of temp pts (radians)
      ,  dxhy       &! 0.5*(HTE - HTE)
      ,  dyhx        ! 0.5*(HTN - HTN)

      real (kind=dbl_kind), dimension  (:,:),allocatable  ::  &
         cyp        &! 1.5*HTE - 0.5*HTE
      ,  cxp        &! 1.5*HTN - 0.5*HTN
      ,  cym        &! 0.5*HTE - 1.5*HTE
      ,  cxm        &! 0.5*HTN - 1.5*HTN
      ,  dxt2       &! 0.5*dxt
      ,  dyt2       &! 0.5*dyt
      ,  dxt4       &! 0.25*dxt
      ,  dyt4       &! 0.25*dyt
      ,  tarear     &! 1/tarea
      ,  uarear     &! 1/uarea
      ,  tinyarea   &! puny*tarea
      ,  ANGLE      &! for conversions between POP grid and lat/lon
      ,  ANGLET     &! ANGLE converted to T-cells
      ,  tarean     &! area of NH T-cells
      ,  tareas      ! area of SH T-cells

      ! Masks
      real (kind=dbl_kind), dimension  (:,:),allocatable ::  &
         hm         &! land/boundary mask, thickness (T-cell)
      ,  uvm        &! land/boundary mask, velocity (U-cell)
      ,  mask_n     &! northern hemisphere
      ,  mask_s      ! southern hemisphere

      logical (kind=log_kind), dimension (:,:),allocatable  ::  &
         tmask     &! land/boundary mask, thickness (T-cell)
      ,  umask     &! land/boundary mask, velocity (U-cell)
      ,  icetmask  &! ice extent mask (T-cell)
      ,  iceumask   ! ice extent mask (U-cell)


      real (kind=dbl_kind)  ::  &
         shlat = -40.0_dbl_kind     &! artificial masking edge
      ,  nhlat =  35.0_dbl_kind      ! artificial masking edge

!=======================================================================

      contains


!=======================================================================



!BOP
!
! !IROUTINE: init_grid - horizontal grid initialization
!
! !INTERFACE:
!
      subroutine init_grid
!
! !DESCRIPTION:
!
! Horizontal grid initialization:
!
!     HT{N,E} = cell widths on {N,E} sides of T cell;
!     U{LAT,LONG} = true {latitude,longitude} of U points;
!     D{X,Y}{T,U} = {x,y} spacing centered at {T,U} points.
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
!
! !USES:
!
!      use ice_exit
!      ggao
!
! !INPUT/OUTPUT PARAMETERS:
!
!
!EOP
!
      integer (kind=int_kind) :: i, j

      real (kind=dbl_kind)  ::  &
           angle_0, angle_w, angle_s, angle_sw

      logical (kind=log_kind), dimension(ilo:ihi,jlo:jhi):: out_of_range

!      if (grid_type == 'displaced_pole') then
!        call popgrid              ! read POP grid lengths directly
!      elseif (grid_type == 'column') then
!        call columngrid           ! column model grid
!      else
!        call rectgrid             ! regular rectangular grid
!      endif

!!--------------------------------------------------------------------------
      !! All FVCOM grid are in water  !  ggao 
      !!--------------------------------------------------------------------------
      !!--------------------------------------------------------------------------
      do j=jlo,jhi
       do i=ilo,ihi
         HM(I,J)= 1.0  

#  if !defined (SPHERICAL) 
      !! 
      write(nu_diag,*)'Please check the domain'
      write(nu_diag,*)' and specify initial condiction for ICE MODEL'
      call PSTOP
!      ULAT(i,J)=90.0_SP/rad_to_deg !! Specify the latitude  !!yding

     
# else
         ULAT(i,J)=VY(I)/rad_to_deg  !vy(I)  !!  spherical coordinate
         ULON(i,J)=VX(I)/rad_to_deg  ! vx(I)
# endif

         TLON(I,J)=ULON(I,J) ! VX(I)/rad_to_deg ! longitude of temp pts (radians)
         TLAT(I,J)=ULAT(I,j) ! VY(I)/rad_to_deg ! latitude of temp pts (radians)
       end do
      end do

!      call bound(HTN)
!      call bound(HTE)
!      call bound(ULAT)
!      call bound(ULON)

      !-----------------------------------------------------------------
      ! construct T-grid cell and U-grid cell widths
      !-----------------------------------------------------------------

      do j=jlo,jhi
       do i=ilo,ihi

!        dxt(i,j) = p5*(HTN(i,j) + HTN(i,j-1))
!        dyt(i,j) = p5*(HTE(i,j) + HTE(i-1,j))

!        tarea(i,j) = dxt(i,j)*dyt(i,j)
!     ggao  AREA OF NODE-BASE CONTROl VOLUME
!        tarea(i,j) = ART1(j)

!        dxu(i,j) = p5*(HTN(i,j) + HTN(i+1,j))
!        dyu(i,j) = p5*(HTE(i,j) + HTE(i,j+1))

!        HTS(i,j) = HTN(i,j-1)
!        HTW(i,j) = HTE(i-1,j)

       enddo
      enddo

      do j=1,jmt_local
       do i=ilo,ihi
!     ggao  AREA OF NODE-BASE CONTROl VOLUME
!        tarea(i,j) = ART1(j)
!        uarea(i,j) = ART1(j)
       enddo
      enddo


!      call bound(dxt)
!      call bound(dyt)
!      call bound(dxu)
!      call bound(dyu)
!      call bound(tarea)
!      call bound(HTS)
!      call bound(HTW)

      do j=jlo,jhi
       do i=ilo,ihi
!        uarea(i,j) = p25*(tarea(i,j) + tarea(i+1,j)  &
!               + tarea(i,j+1) + tarea(i+1,j+1))        ! m^2
       enddo
      enddo
!      call bound(uarea)

      ! grid length combinations
      do j=jlo,jhi
       do i=ilo,ihi
!        dxt2(i,j) = 0.5*dxt(i,j)
!        dyt2(i,j) = 0.5*dyt(i,j)
!        dxt4(i,j) = 0.25*dxt(i,j)
!        dyt4(i,j) = 0.25*dyt(i,j)
!        tarear(i,j) = 1./tarea(i,j)
!        uarear(i,j) = 1./uarea(i,j)
!        tinyarea(i,j) = puny*tarea(i,j)

!        cyp(i,j) = (1.5*HTE(i,j) - 0.5*HTE(i-1,j))
!        cxp(i,j) = (1.5*HTN(i,j) - 0.5*HTN(i,j-1))
!        cym(i,j) = (0.5*HTE(i,j) - 1.5*HTE(i-1,j))
!        cxm(i,j) = (0.5*HTN(i,j) - 1.5*HTN(i,j-1))

!        dxhy(i,j) = 0.5*(HTE(i,j) - HTE(i-1,j))
!        dyhx(i,j) = 0.5*(HTN(i,j) - HTN(i,j-1))
       enddo
      enddo

!      call bound(dxhy)
!      call bound(dyhx)

!-----------------------------------------------------------------
! Calculate ANGLET to be compatible with POP ocean model
! First, ensure that -pi <= ANGLE <= pi
!-----------------------------------------------------------------

!      out_of_range = .false.
!      where (ANGLE < -pi .or. ANGLE > pi) out_of_range = .true.
!      if (count(out_of_range) > 0) then
!         call abort_ice ('init_grid: ANGLE out of expected range')
!      endif
!     ggao


!-----------------------------------------------------------------
! Pad ANGLE so ghost cells can be used in averaging to get ANGLET
!-----------------------------------------------------------------

!    ggao 0105-2007
!      work_l1(ilo:ihi,jlo:jhi) = ANGLE(ilo:ihi,jlo:jhi)
!      call bound (work_l1)
!    change end

!-----------------------------------------------------------------
! Compute ANGLE on T-grid
!-----------------------------------------------------------------
!      do j=jlo,jhi
!        do i=ilo,ihi
!          angle_0  = work_l1(i  ,j  )     !   w----0
!          angle_w  = work_l1(i-1,j  )     !   |    |
!          angle_s  = work_l1(i,  j-1)     !   |    |
!          angle_sw = work_l1(i-1,j-1)     !   sw---s

!          if ( angle_0 < c0i ) then              
!            if ( abs(angle_w - angle_0) > pi)   &
!              angle_w  = angle_w  - pi2          
!            if ( abs(angle_s - angle_0) > pi)   &
!              angle_s  = angle_s  - pi2          
!            if ( abs(angle_sw - angle_0) > pi)  &
!              angle_sw = angle_sw - pi2
!          endif

!          ANGLET(i,j) = angle_0 * p25 + angle_w * p25  &
!                      + angle_s * p25 + angle_sw* p25
!         ANGLET(i,j) =0.0
!        enddo
!      enddo
!   ggao

!      call Tlatlon           ! get lat, lon on the T grid
!      call makemask          ! velocity mask, hemisphere masks

!      do j=1,jmt_global
!      do i=1,imt_global
!        work_g1(i,j) = float((j-1)*imt_global + i)
!      enddo
!      enddo
!      call global_scatter(work_g1,rndex_global)
!      index_global = nint(rndex_global)
!   ggao

      end subroutine init_grid

!=======================================================================
!BOP
!
! !IROUTINE: popgrid - reads and sets pop displaced pole grid and land mask
!
! !INTERFACE:
!
      subroutine popgrid
!
! !DESCRIPTION:
!
! POP displaced pole grid and land mask.
! Grid record number, field and units are:
! (1) ULAT  (radians)
! (2) ULON  (radians)
! (3) HTN   (cm)  
! (4) HTE   (cm)
! (5) HUS   (cm)
! (6) HUW   (cm)
! (7) ANGLE (radians)
!
! Land mask record number and field is (1) KMT.
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
!
! !USES:
!
!      use ice_read_write
!     ggao
!
! !INPUT/OUTPUT PARAMETERS:
!
!
!EOP
!
      integer (kind=int_kind) :: i, j
      logical (kind=log_kind) :: scatter, diag

!      change 
!      call ice_open(nu_grid,grid_file,64)
!      call ice_open(nu_kmt,kmt_file,32)

      scatter = .true.    ! scatter data to all processors
      diag = .true.       ! write diagnostic info

!      call ice_read(nu_grid,1,worka,'rda8',scatter,diag)
!      ULAT(ilo:ihi,jlo:jhi)=worka(ilo:ihi,jlo:jhi)
!      call ice_read(nu_grid,2,worka,'rda8',scatter,diag)
!      ULON(ilo:ihi,jlo:jhi)=worka(ilo:ihi,jlo:jhi)
!      call ice_read(nu_grid,3,worka,'rda8',scatter,diag)
!      HTN(ilo:ihi,jlo:jhi)=worka(ilo:ihi,jlo:jhi)*cm_to_m
!      call ice_read(nu_grid,4,worka,'rda8',scatter,diag)
!      HTE(ilo:ihi,jlo:jhi)=worka(ilo:ihi,jlo:jhi)*cm_to_m
!      call ice_read(nu_grid,7,worka,'rda8',scatter,diag)
!      ANGLE(ilo:ihi,jlo:jhi)=worka(ilo:ihi,jlo:jhi)
!      call ice_read(nu_kmt,1,worka,'ida4',scatter,diag)

      if (my_task == master_task) then
!         close (nu_grid) 
!         close (nu_kmt) 
      endif

!       ggao change end

      do j=jlo,jhi
       do i=ilo,ihi
!         hm(i,j) = worka(i,j)
!         if (hm(i,j) >= c1i) hm(i,j) = c1i 

         ! uncomment to mask out tropics
         ! Do this only if running uncoupled
!!!          if (ULAT(i,j) > shlat/rad_to_deg .and.
!!!     &        ULAT(i,j) < nhlat/rad_to_deg) hm(i,j) = c0i
       enddo
      enddo

      end subroutine popgrid

!=======================================================================
!BOP
!
! !IROUTINE: columngrid - column grid and mask
!
! !INTERFACE:
!
      subroutine columngrid
!
! !DESCRIPTION:
!
! Column grid and mask
!
! !REVISION HISTORY:
!
! author: C. M. Bitz UW, (based on rectgrid by Hunke)
!
! modified Nov. 2003 by William H. Lipscomb, LANL
!
! !USES:
!
      use ice_model_size
!      use ice_exit
!     ggao

!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) :: i, j

      !-----------------------------------------------------------------
      ! Calculate various geometric 2d arrays
      !-----------------------------------------------------------------

!      do j=jlo,jhi
!      do i=ilo,ihi
!         HTN  (i,j) = 1.6e4_dbl_kind   ! constant longitude spacing = 
!                                       ! meaningless
!         HTE  (i,j) = 1.6e4_dbl_kind   ! constant latitude  spacing = 
!                                       ! meaningless
!         ULAT (i,j) = 75.5/rad_to_deg  ! used to determine hemisphere and 
!         ULON (i,j) = 170.0/rad_to_deg ! init_state, need not be exact
!         ANGLE(i,j) = c0i               ! "square with the world"
!      enddo
!      enddo

      !-----------------------------------------------------------------
      ! Verify that imt_global and jmt_global are 1
      !-----------------------------------------------------------------

!      if ((imt_global /= 1).or. (jmt_global /= 1)) then
!        write(nu_diag,*)  & 
!             'Because you have selected the column model flag'
!         write(nu_diag,*) 'Please set imt_global=jmt_global=1 in file'
!         write(nu_diag,*) 'ice_model_size.F and recompile'
!         call abort_ice ('columngrid: check imt_global and jmt_global')
!      endif

      !-----------------------------------------------------------------
      ! Construct T-cell land mask
      !-----------------------------------------------------------------

!      do j=1,jmt_global
!      do i=1,imt_global
!         work_g1(i,j) = c1i
!      enddo
!      enddo

!      call global_scatter(work_g1,worka)

      do j=jlo,jhi
      do i=ilo,ihi
!         hm(i,j) = worka(i,j)
!   for fvcom 1D
         hm(i,j) = c1i  !worka(i,j)
      enddo
      enddo

      end subroutine columngrid

!=======================================================================
!BOP
!
! !IROUTINE: rectgrid - regular rectangular grid and mask
!
! !INTERFACE:
!
      subroutine rectgrid
!
! !DESCRIPTION:
!
! Regular rectangular grid and mask
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
!
! !USES:
!
      use ice_model_size
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) :: i, j

      !-----------------------------------------------------------------
      ! Calculate various geometric 2d arrays
      !-----------------------------------------------------------------

      do j=jlo,jhi
      do i=ilo,ihi
!!!         HTN  (i,j) = 3.1e4_dbl_kind  ! constant longitude spacing =  
                                         ! POP <4/3> min, m
!!!         HTE  (i,j) = 3.1e4_dbl_kind  ! constant latitude  spacing = 
                                         ! POP <4/3> min, m
!         HTN  (i,j) = 1.6e4_dbl_kind   ! constant longitude spacing = 
                                       ! POP <2/3> min, m
!         HTE  (i,j) = 1.6e4_dbl_kind   ! constant latitude  spacing = 
                                       ! POP <2/3> min, m
!         ULAT (i,j) = c0i               ! remember to set Coriolis !
!         ULON (i,j) = c0i
!         ANGLE(i,j) = c0i               ! "square with the world"
      enddo
      enddo

      !-----------------------------------------------------------------
      ! Construct T-cell land mask
      !-----------------------------------------------------------------

      do j=1,jmt_global         ! initialize hm as land
      do i=1,imt_global
!         work_g1(i,j) = c0i
      enddo
      enddo

!!!      do j=1,jmt_global        ! open
!!!      do i=1,imt_global        ! open
      do j=3,jmt_global-2       ! closed: NOTE jmt_global > 5
      do i=3,imt_global-2       ! closed: NOTE imt_global > 5
!         work_g1(i,j) = c1i
      enddo
      enddo

!      call global_scatter(work_g1,worka)

      do j=jlo,jhi
      do i=ilo,ihi
!         hm(i,j) = worka(i,j)
      enddo
      enddo

      end subroutine rectgrid

!=======================================================================
!BOP
!
! !IROUTINE: makemask - makes logical land masks (T,U) and hemispheric masks
!
! !INTERFACE:
!
      subroutine makemask
!
! !DESCRIPTION:
!
! Sets the boundary values for the T cell land mask (hm) and
! makes the logical land masks for T and U cells (tmask, umask).
! Also creates hemisphere masks (mask-n northern, mask-s southern)
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) :: i, j

!      call bound(hm)  !!! use real arrays to get boundary conditions

      !-----------------------------------------------------------------
      ! construct T-cell and U-cell masks
      !-----------------------------------------------------------------

      do j=jlo,jhi
       do i=ilo,ihi
         !! T and UV at the node 
         !! All FVCOM grid are in water  !  ggao 
!         UVM(i,J) =HM(i,J) 
!        uvm(i,j) = min(hm(i,j),hm(i+1,j),hm(i,j+1),hm(i+1,j+1))
!        
       enddo
      enddo
!      call bound(uvm)  !!! use real arrays to get boundary conditions



      do j=1,jmt_local
       do i=1,imt_local
!        tmask(i,j) = .false.
!        umask(i,j) = .false.
!        if ( hm(i,j) > p5) tmask(i,j) = .true. 
!        if (uvm(i,j) > p5) umask(i,j) = .true. 
!!--------------------------------------------------------------------------
      !! All FVCOM grid are in water  !  ggao 
      !!--------------------------------------------------------------------------
      !!--------------------------------------------------------------------------
!         tmask(i,j) = .true.
!         umask(i,j) = .true.

       enddo
      enddo

      !-----------------------------------------------------------------
      ! create hemisphere masks
      !-----------------------------------------------------------------

      do j=1,jmt_local
       do i=1,imt_local
!        mask_n(i,j) = c0i
!        mask_s(i,j) = c0i
       enddo
      enddo
      do j=jlo,jhi
       do i=ilo,ihi
!        if (ULAT(i,j) >= -puny) mask_n(i,j) = c1i  ! northern hemisphere
!        if (ULAT(i,j) <  -puny) mask_s(i,j) = c1i  ! southern hemisphere

!        tarean(i,j) = tarea(i,j)*mask_n(i,j)  ! N hemisphere area mask (m^2)
!        tareas(i,j) = tarea(i,j)*mask_s(i,j)  ! S hemisphere area mask (m^2)
       enddo
      enddo

      end subroutine makemask

!=======================================================================
!BOP
!
! !IROUTINE: Tlatlon - initializes latitude and longitudes on T grid
!
! !INTERFACE:
!
      subroutine Tlatlon
!
! !DESCRIPTION:
!
! Initializes latitude and longitude on T grid
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL; code originally based on POP grid 
! generation routine
!
! !USES:
!
      use ice_model_size
!      use ice_read_write ! if reading ULAT, ULON directly from file
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) ::  &
        i, j                           ! horizontal indices

      integer (kind=int_kind) :: im1
      real (kind=dbl_kind) ::  &
        z1,x1,y1,z2,x2,y2,z3,x3,y3,z4,x4,y4,tx,ty,tz,da

!      allocate (work_g2(imt_global,jmt_global))

!      call global_gather(work_g1,ULON(ilo:ihi,jlo:jhi))
!      call global_gather(work_g2,ULAT(ilo:ihi,jlo:jhi))

      if (my_task == master_task) then

      do j=2,jmt_global
        do i=1,imt_global

            if (i==1) then
               im1=imt_global
            else
               im1=i-1
            endif

!            z1 = cos(work_g2(im1,j-1))
!            x1 = cos(work_g1(im1,j-1))*z1
!            y1 = sin(work_g1(im1,j-1))*z1
!            z1 = sin(work_g2(im1,j-1))

!            z2 = cos(work_g2(i,j-1))
!            x2 = cos(work_g1(i,j-1))*z2
!            y2 = sin(work_g1(i,j-1))*z2
!            z2 = sin(work_g2(i,j-1))

!            z3 = cos(work_g2(im1,j))
!            x3 = cos(work_g1(im1,j))*z3
!            y3 = sin(work_g1(im1,j))*z3
!            z3 = sin(work_g2(im1,j))

!            z4 = cos(work_g2(i,j))
!            x4 = cos(work_g1(i,j))*z4
!            y4 = sin(work_g1(i,j))*z4
!            z4 = sin(work_g2(i,j))

!            tx = (x1+x2+x3+x4)/c4i
!            ty = (y1+y2+y3+y4)/c4i
!            tz = (z1+z2+z3+z4)/c4i
!            da = sqrt(tx**2+ty**2+tz**2)

!            tz = tz/da

            ! TLON_G in radians East
!            TLON_G(i,j) = c0i
!            if (tx /= c0i .or. ty /= c0i) TLON_G(i,j) = atan2(ty,tx)

            ! TLAT_G in radians North
!            TLAT_G(i,j) = asin(tz)

        end do
      end do

      ! j=1: linear approximation
      do i=1,imt_global
!         TLON_G(i,1) = TLON_G(i,2)
!         TLAT_G(i,1) = c2i*TLAT_G(i,2) - TLAT_G(i,3)
      end do

!      write(nu_diag,*) ''
!      write(nu_diag,*) 'min/max ULON_G:',minval(work_g1),maxval(work_g1)
!      write(nu_diag,*) 'min/max TLON_G:',minval(TLON_G),maxval(TLON_G)
!      write(nu_diag,*) 'min/max ULAT_G:',minval(work_g2),maxval(work_g2)
!      write(nu_diag,*) 'min/max TLAT_G:',minval(TLAT_G),maxval(TLAT_G)

      endif ! master_task

!      deallocate (work_g2)

!      call global_scatter(TLON_G,worka)
      do j=jlo,jhi
      do i=ilo,ihi
!        TLON(i,j) = worka(i,j)
      enddo
      enddo

!      call global_scatter(TLAT_G,worka)
      do j=jlo,jhi
      do i=ilo,ihi
!        TLAT(i,j) = worka(i,j)
      enddo
      enddo

!      call bound(TLON)
!      call bound(TLAT)

      end subroutine Tlatlon

!=======================================================================
!BOP
!
! !IROUTINE: t2ugrid - transfer from T-cell centers to U-cell centers
!
! !INTERFACE:
!
      subroutine t2ugrid(work)
!
! !DESCRIPTION:
!
! Transfer from T-cell centers to U-cell centers. Writes work into another 
! array that has ghost cells
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind) :: work (ilo:ihi,jlo:jhi)
!
!EOP
!
      integer (kind=int_kind) :: i, j

      do j=jlo,jhi
       do i=ilo,ihi
!        work_l1(i,j) = work(i,j)
       enddo
      enddo
!      call bound(work_l1)
!      call to_ugrid(work_l1,work)

      end subroutine t2ugrid

!=======================================================================
!BOP
!
! !IROUTINE: to_ugrid - shift from T-cell to U-cell midpoints
!
! !INTERFACE:
!
      subroutine to_ugrid(work1,work2)
!
! !DESCRIPTION:
!
! Shifts quantities from the T-cell midpoint (work1) to the U-cell 
! midpoint (work2)
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind) :: work1(imt_local,jmt_local)  &
      ,                       work2(ilo:ihi,jlo:jhi)
!
!EOP
!
      integer (kind=int_kind) :: i, j

      do j=jlo,jhi
       do i=ilo,ihi
!       work2(i,j) = p25*(work1(i,j)*tarea(i,j)                      &
!                         + work1(i+1,j)*tarea(i+1,j)                &
!                         + work1(i,j+1)*tarea(i,j+1)                &
!                         + work1(i+1,j+1)*tarea(i+1,j+1))/uarea(i,j)
       enddo
      enddo

      end subroutine to_ugrid

!=======================================================================
!BOP
!
! !IROUTINE: u2tgrid - transfer from U-cell centers to T-cell centers
!
! !INTERFACE:
!
      subroutine u2tgrid(work)
!
! !DESCRIPTION:
!
! Transfer from U-cell centers to T-cell centers. Writes work into 
! another array that has ghost cells
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind) :: work (ilo:ihi,jlo:jhi)
!
!EOP
!
      integer (kind=int_kind) :: i, j

      do j=jlo,jhi
       do i=ilo,ihi
!        work_l1(i,j) = work(i,j)
       enddo
      enddo
!      call bound(work_l1)
!      call to_tgrid(work_l1,work)

      end subroutine u2tgrid

!=======================================================================
!BOP
!
! !IROUTINE: to_tgrid - shifts array from U-cell to T-cell midpoints
!
! !INTERFACE:
!
      subroutine to_tgrid(work1,work2)
!
! !DESCRIPTION:
!
! Shifts quantities from the U-cell midpoint (work1) to the T-cell 
! midpoint (work2)
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind) :: work1(imt_local,jmt_local) &
      ,                       work2(ilo:ihi,jlo:jhi)
!
!EOP
!
      integer (kind=int_kind) :: i, j

      do j=jlo,jhi
       do i=ilo,ihi
!       work2(i,j) = p25*(work1(i,j)    * uarea(i,j)     &
!                       + work1(i-1,j)  * uarea(i-1,j)   &
!                       + work1(i,j-1)  * uarea(i,j-1)   &
!                       + work1(i-1,j-1)* uarea(i-1,j-1))/tarea(i,j)
       enddo
      enddo

      end subroutine to_tgrid

!=======================================================================
!BOP
!
! !IROUTINE: bound - fills ghost cells with boundary information
!
! !INTERFACE:
!
      subroutine bound(work1)
!
! !DESCRIPTION:
!
! Fills ghost cells with boundary information
!
! !REVISION HISTORY:
!
! author: Tony Craig, NCAR
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind) :: work1(1)
!
!EOP
!
!      call bound_ijn(1,work1,.true.,.true.,.true.,.true.)

      end subroutine bound

!=======================================================================
!BOP
!
! !IROUTINE: bound_sw - fills south and west ghost cells
!
! !INTERFACE:
!
      subroutine bound_sw(work1)
!
! !DESCRIPTION:
!
! Fills south and west ghost cells with boundary information
!
! !REVISION HISTORY:
!
! author: Tony Craig, NCAR
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind) :: work1(1)
!
!EOP
!
      call bound_ijn(1,work1,.false.,.true.,.false.,.true.)

      end subroutine bound_sw

!=======================================================================
!BOP
!
! !IROUTINE: bound_narr - fills neighboring ghost cells with boundary info
!
! !INTERFACE:
!
      subroutine bound_narr(narrays,work1)
!
! !DESCRIPTION:
!
! Fills neighboring ghost cells with boundary information; 
! several arrays at once (for performance)
!
! !REVISION HISTORY:
!
! authors: Tony Craig, NCAR
!          Elizabeth C. Hunke, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind) :: narrays
      real (kind=dbl_kind) :: work1(1)
!
!EOP
!
!      call bound_ijn(narrays,work1,.true.,.true.,.true.,.true.)

      end subroutine bound_narr

!=======================================================================
!BOP
!
! !IROUTINE: bound_narr_ne - fills north and east ghost cells
!
! !INTERFACE:
!
      subroutine bound_narr_ne(narrays,work1)
!
! !DESCRIPTION:
!
! Fills north and east ghost cells with boundary information;
! several arrays at once (for performance)
!
! !REVISION HISTORY:
!
! authors: Tony Craig, NCAR
!          Elizabeth C. Hunke, LANL
!
! modified Nov. 2003 by William H. Lipscomb, LANL 
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind) :: narrays
      real (kind=dbl_kind) :: work1(1)
!
!EOP
!
!      call bound_ijn(narrays,work1,.true.,.false.,.true.,.false.)

      end subroutine bound_narr_ne

!=======================================================================
!BOP
!
! !IROUTINE: bound_ijn - Periodic/Neumann boundary conditions
!
! !INTERFACE:
!
      subroutine bound_ijn(nd,work1,north,south,east,west)
!
! !DESCRIPTION:
!
! Periodic/Neumann conditions for global domain boundaries. \\
! Assumptions:  a *single* row of ghost cells (num-ghost-cells=1); \\
! work1 array has form (i-index,j-index,number-arrays)
!
! !REVISION HISTORY:
!
! authors: Tony Craig, NCAR
!          Elizabeth Hunke, LANL
!
! !USES:
!
!      use ice_timers

!      ggao
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind) :: nd
      real (kind=dbl_kind) :: work1(imt_local,jmt_local,nd)
      logical (kind=log_kind) :: north,south,east,west
!
!EOP
!
      integer (kind=int_kind) :: i, j, ni

!#ifdef _MPI
!      integer (kind=int_kind) :: icnt,jcnt  &
!      ,   status(MPI_STATUS_SIZE),request(4)
!      real (kind=dbl_kind) :: workw(jlo:jhi,nd),worke(jlo:jhi,nd) &
!      ,  workn(ilo-1:ihi+1,nd),works(ilo-1:ihi+1,nd)
!#endif

!      call ice_timer_start(10) ! bound

!#ifdef _MPI
!      jcnt=(jhi-jlo+1)*nd
!      icnt=(ihi-ilo+1+2*num_ghost_cells)*nd
!
      !-----------------------------------------------------------------
      ! west data to east data, west shift
      !-----------------------------------------------------------------
      if (east) then

      do ni=1,nd
      do j=jlo,jhi
!      workw(j,ni)=work1(ilo,j,ni)
      enddo
      enddo

!      call MPI_SENDRECV(workw,jcnt,MPI_REAL8,nbr_west,my_task,  &
!                        worke,jcnt,MPI_REAL8,nbr_east,nbr_east, &
!                        MPI_COMM_ICE,status,ierr)

      do ni=1,nd
      do j=jlo,jhi
!      work1(ihi+1,j,ni)=worke(j,ni)
      enddo
      enddo

      endif

      !-----------------------------------------------------------------
      ! east data to west data, east shift
      !-----------------------------------------------------------------
!      if (west) then

      do ni=1,nd
      do j=jlo,jhi
!      worke(j,ni)=work1(ihi,j,ni)
      enddo
      enddo

!      call MPI_SENDRECV(worke,jcnt,MPI_REAL8,nbr_east,my_task, &
!                        workw,jcnt,MPI_REAL8,nbr_west,nbr_west,&
!                        MPI_COMM_ICE,status,ierr)

      do ni=1,nd
      do j=jlo,jhi
!      work1(ilo-1,j,n)=workw(j,ni)
      enddo
      enddo

!      endif

      !-----------------------------------------------------------------
      ! north data to south data, north shift
      !-----------------------------------------------------------------
!      if (south) then
!      if (nbr_south /= -1) then
!        call MPI_IRECV(works,                              &
!                       icnt,MPI_REAL8,nbr_south,nbr_south, &
!                       MPI_COMM_ICE,request(1),ierr)
!      else

!        do ni=1,nd
!        do i=ilo-1,ihi+1
!        work1(i,jlo-1,ni)=work1(i,jlo,ni)
!        enddo
!        enddo

!      endif
!      if (nbr_north /= -1) then

!        do ni=1,nd
!        do i=ilo-1,ihi+1
!        workn(i,n)=work1(i,jhi,ni)
!        enddo
!        enddo

!        call MPI_ISEND (workn,                          &
!                       icnt,MPI_REAL8,nbr_north,my_task,&
!                       MPI_COMM_ICE,request(2),ierr)
!      endif
!      if (nbr_north /= -1) then
!        call MPI_WAIT(request(2), status, ierr)
!      endif
!      if (nbr_south /= -1) then
!        call MPI_WAIT(request(1), status, ierr)

!        do ni=1,nd
!        do i=ilo-1,ihi+1
!        work1(i,jlo-1,ni)=works(i,ni)
!        enddo
!        enddo

!      endif
!      endif

      !-----------------------------------------------------------------
      ! south data to north data, south shift
      !-----------------------------------------------------------------
!      if (north) then
!      if (nbr_north /= -1) then
!        call MPI_IRECV(workn,                              &
!                       icnt,MPI_REAL8,nbr_north,nbr_north, &
!                       MPI_COMM_ICE,request(3),ierr)
!      else

        do ni=1,nd
        do i=ilo-1,ihi+1
!        work1(i,jhi+1,ni)=work1(i,jhi,ni)
        enddo
        enddo

!      endif
!      if (nbr_south /= -1) then

        do ni=1,nd
        do i=ilo-1,ihi+1
!        works(i,ni)=work1(i,jlo,ni)
        enddo
        enddo

!        call MPI_ISEND (works,                           &
!                       icnt,MPI_REAL8,nbr_south,my_task, &
!                       MPI_COMM_ICE,request(4),ierr)
!      endif
!      if (nbr_south /= -1) then
!        call MPI_WAIT(request(4), status, ierr)
!      endif
!      if (nbr_north /= -1) then
!        call MPI_WAIT(request(3), status, ierr)

        do ni=1,nd
        do i=ilo-1,ihi+1
!        work1(i,jhi+1,ni)=workn(i,ni)
        enddo
        enddo

!      endif
!      endif

!#else
      !-----------------------------------------------------------------
      ! single domain
      !-----------------------------------------------------------------
 
      do ni=1,nd

      ! Periodic conditions
      do j=jlo,jhi
!       work1(ilo-1,j,ni) = work1(ihi,j,ni)
!       work1(ihi+1,j,ni) = work1(ilo,j,ni)
      enddo

      ! Neumann conditions (POP grid land points)
      do i=ilo-1,ihi+1
!        work1(i,jlo-1,ni) = work1(i,jlo,ni)
!        work1(i,jhi+1,ni) = work1(i,jhi,ni)
      enddo

      enddo  ! ni
!#endif
!      call ice_timer_stop(10)  ! bound

      end subroutine bound_ijn

!=======================================================================

      end module ice_grid

!=======================================================================
