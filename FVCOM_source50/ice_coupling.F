!/===========================================================================/
! CVS VERSION INFORMATION
! $Id$
! $Name$
! $Revision$
!/===========================================================================/

!!=====================================================================
!!=====================================================================
!!   the modular is modified to couple with FVCOM
!!
!!=====================================================================
!!=====================================================================

!=======================================================================
!
!BOP
!
! !MODULE: ice_coupling - message passing to and from the coupler
!
! !DESCRIPTION:
!
! Message passing to and from the coupler
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
!         Tony Craig, NCAR, Dec-30-2002, modified for cpl6
!
! !INTERFACE:
!
      module ice_coupling
!
! !USES:
!
      use ice_kinds_mod
      use ice_model_size
      use ice_constants
      use ice_calendar
      use ice_grid
      use ice_state
      use ice_flux
      use ice_albedo
!      use ice_mpi_internal
!      use ice_timers
      use ice_fileunits
      use ice_work, only: worka, work_l1
!#ifdef coupled
!      use shr_sys_mod, only : shr_sys_flush
!      use ice_history, only : runtype
!      use cpl_contract_mod
!      use cpl_interface_mod
!      use cpl_fields_mod
!#endif
!
!EOP
!
         use ice_init
         use ice_flux_in

!----------------------------------------------------------------------
!----------------------------------------------------------------------
     use all_vars
#if defined (MULTIPROCESSOR)
     USE MOD_PAR
     USE lims, only: m,mt,n,nt
#  endif
!     USE MOD_ICE2D, ONLY: UICE2,VICE2

!====================================================================

!#ifdef coupled  
! #endif at end of module

      implicit none

!      integer (kind=int_kind), dimension (cpl_fields_ibuf_total) ::&
!        isbuf                                                     &
!      ,  irbuf                                                      
                                                                    
      real (kind=dbl_kind), allocatable ::                         &
         sbuf(:,:)                                                  
                                                                    
!      real (kind=dbl_kind) ::                                      &
!         buffs((jhi-jlo+1)*(ihi-ilo+1),cpl_fields_i2c_total)       &
!      ,  buffr((jhi-jlo+1)*(ihi-ilo+1),cpl_fields_c2i_total)        
                                                                    
!      type(cpl_contract) ::                                        &
!         contractS                                                 &
!      ,  contractR                                                  
                                                                    
      integer(kind=int_kind), save ::                              &
         nadv_i                                                    &
      ,  info_dbug

!=======================================================================

      contains

!=======================================================================
!BOP
!
! !IROUTINE: ice_coupling_setup - sets mpi communicators and task ids
!
! !INTERFACE:
!
!      subroutine ice_coupling_setup(in_model_name,model_comm)
!
! !DESCRIPTION:
!
! This routine uses get the model communicator from ccsm share code
!
! !REVISION HISTORY:
!
! author: T. Craig, NCAR, Dec 30, 2002: for cpl6
!
!
! !INPUT/OUTPUT PARAMETERS:
!
!      character (3), intent(in) :: in_model_name   

!      integer, intent(out) ::  model_comm     ! communicator for model
!
!EOP
!
!      write(nu_diag,*) 'calling cpl_interface_init for model: ',
!     &     in_model_name,' ', trim(cpl_fields_icename)
!
!      call cpl_interface_init(cpl_fields_icename,model_comm)

!      call shr_sys_flush(nu_diag)

!      end subroutine ice_coupling_setup

!=======================================================================
!BOP
!
! !IROUTINE: init_cpl - initializes message passing between ice and coupler
!
! !INTERFACE:
!
      subroutine init_cpl
!
! !DESCRIPTION:
!
! Initializes message passing between ice and coupler
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer(kind=int_kind) :: i,j,icn     ! local loop indices

!      write(nu_diag,*) 
!     &     '(ice_coupling,init_cpl) send initial msg. set contract'
!      call shr_sys_flush(nu_diag)

!      nadv_i = nint(secday/dt)

!      isbuf                          = 0         ! default info-buffer value
!      isbuf(cpl_fields_ibuf_cdate  ) = idate     ! initial date (coded: yyyymmdd)

!      isbuf(cpl_fields_ibuf_sec    ) = sec       ! elapsed seconds into date
!      isbuf(cpl_fields_ibuf_stopnow) = stop_now  ! stop now flag
!      isbuf(cpl_fields_ibuf_userest) = 0         ! use model restart data initally
!      isbuf(cpl_fields_ibuf_ncpl   ) = nadv_i    ! number of comms per day
!      isbuf(cpl_fields_ibuf_lsize  ) = (ihi-ilo+1)*(jhi-jlo+1) ! size of local grid
!      isbuf(cpl_fields_ibuf_lisize ) = (ihi-ilo+1) ! local size wrt i-index
!      isbuf(cpl_fields_ibuf_ljsize ) = (jhi-jlo+1) ! local size wrt i-index
!      isbuf(cpl_fields_ibuf_gsize  ) = imt_global*jmt_global ! size of global grid
!      isbuf(cpl_fields_ibuf_gisize ) = imt_global  ! global size wrt i-index
!      isbuf(cpl_fields_ibuf_gjsize ) = jmt_global  ! global size wrt j-index
!      isbuf(cpl_fields_ibuf_nfields) = cpl_fields_grid_total
!      isbuf(cpl_fields_ibuf_dead   ) = 0           ! not a dead model

!      allocate(sbuf((ihi-ilo+1)*(jhi-jlo+1),cpl_fields_grid_total))
      sbuf = -888.0
      icn=0
      do j=jlo,jhi
      do i=ilo,ihi
!         in=in+1
!         sbuf(in,cpl_fields_grid_lon  ) = TLON(i,j)*rad_to_deg
!         sbuf(in,cpl_fields_grid_lat  ) = TLAT(i,j)*rad_to_deg
!         sbuf(in,cpl_fields_grid_area ) = tarea(i,j)/(radius*radius)
!         sbuf(in,cpl_fields_grid_mask ) = float(nint(hm(i,j)))
!         sbuf(in,cpl_fields_grid_index) = rndex_global(i,j)
      enddo
      enddo

!      call cpl_interface_contractInit
!     &     (contractS, cpl_fields_icename, cpl_fields_cplname,
!     &      cpl_fields_i2c_fields, isbuf, sbuf)

!      call cpl_interface_contractInit
!     &     (contractR, cpl_fields_icename, cpl_fields_cplname,
!     &      cpl_fields_c2i_fields, isbuf, sbuf)

      write(nu_diag,*) '(init_cpl) Initialized contracts with coupler'
!      call shr_sys_flush(nu_diag)

      !-----------------------------------------------------------------
      ! Receive initial message from coupler.
      !-----------------------------------------------------------------

!      call cpl_interface_ibufRecv(cpl_fields_cplname,irbuf)

!      if (my_task==master_task) then
!         write(nu_diag,*)
!     &        '(init_cpl) Received control buffer from coupler'
!         call shr_sys_flush(nu_diag)

!         if (trim(runtype)=='startup' .or.
!     &       trim(runtype)== 'hybrid') then
!            idate = irbuf(cpl_fields_ibuf_cdate)
!            write(nu_diag,*) '(init_cpl) idate from coupler = ',idate
!            nyr   = (idate/10000)               ! integer year of basedate
!            month = (idate-nyr*10000)/100       ! integer month of basedate
!            mday  = idate-nyr*10000-month*100-1 ! day of year of basedate
!            time  = ((nyr-1)*daycal(13)+daycal(month)+mday)*secday
!            call calendar(time)                 ! recompute calendar info
!            time_forc = time
!            call shr_sys_flush(nu_diag)
!         endif

!      endif                     ! my_task==master_task

!      call ice_bcast_iscalar(idate)
!      call ice_bcast_rscalar(time)
!      call ice_bcast_rscalar(time_forc)

!      deallocate(sbuf)

!      write(nu_diag,*) '(ice_coupling,init_cpl) done setting contract'

      !-----------------------------------------------------------------
      ! Send initial state info to coupler.
      !-----------------------------------------------------------------

!      call to_coupler

      end subroutine init_cpl

!=======================================================================
!BOP
!
! !IROUTINE: from_coupler - input from coupler to sea ice model
!
! !INTERFACE:
!
      subroutine from_coupler
!
! !DESCRIPTION:
!
! Reads input data from coupler to sea ice model
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) :: i,j,icn,n2   ! local loop indices

      real (kind=dbl_kind) ::  &
        gsum, workx, worky
!--------ice_ocean---------------------------------------------------
      real (kind=dbl_kind), parameter ::  &
        cphm = cp_ocn*rhow*20 !hmix
!         hmix = 20._dbl_kind         &    ! ocean mixed layer depth (m)
!--------ice_ocean---------------------------------------------------

       REAL, dimension (nt) :: utmpc,vtmpc
       REAL, dimension(mt) :: utmp,vtmp
       REAL, allocatable, dimension (:,:) :: rbuf

!      call ice_timer_start(8)  ! time spent coupling

      !-----------------------------------------------------------------
      ! Zero stuff while waiting, only filling in active cells.
      !-----------------------------------------------------------------

      zlvl(:,:)    = c0i
      uatm(:,:)    = c0i
      vatm(:,:)    = c0i
      potT(:,:)    = c0i
      Tair(:,:)    = c0i
      Qa(:,:)      = c0i
      rhoa(:,:)    = c0i
      swvdr(:,:)   = c0i
      swvdf(:,:)   = c0i
      swidr(:,:)   = c0i
      swidf(:,:)   = c0i
      flw(:,:)     = c0i
      frain(:,:)   = c0i
      fsnow(:,:)   = c0i
      sst(:,:)     = c0i
      sss(:,:)     = c0i
!      uocn(:,:)    = c0i
!      vocn(:,:)    = c0i
!      ss_tltx(:,:) = c0i
!      ss_tlty(:,:) = c0i
      frzmlt(:,:)  = c0i

!      first step initialized the ice variables

      !-----------------------------------------------------------------
      ! recv input field msg
      !-----------------------------------------------------------------
     
!      call ice_timer_start(16)  ! time spent receiving

!      call cpl_interface_contractRecv
!     &     (cpl_fields_cplname, contractR, irbuf, buffr)

!      call ice_timer_stop(16)
!      call ice_timer_start(17)  ! time spent cr-unpacking

      !--- unpack message
      icn=0
      do j=jlo,jhi
      do i=ilo,ihi
         icn=icn+1

         !--- ocn states--
!         sst  (i,j) = buffr(icn,cpl_fields_c2i_ot)
!         sss  (i,j) = buffr(icn,cpl_fields_c2i_os)
!         uocn (i,j) = buffr(icn,cpl_fields_c2i_ou)
!         vocn (i,j) = buffr(icn,cpl_fields_c2i_ov)

         !--- atm states-
!         zlvl (i,j) = buffr(icn,cpl_fields_c2i_z)
!         uatm (i,j) = buffr(icn,cpl_fields_c2i_u)
!         vatm (i,j) = buffr(icn,cpl_fields_c2i_v)
!         potT (i,j) = buffr(icn,cpl_fields_c2i_ptem)
!         Tair (i,j) = buffr(icn,cpl_fields_c2i_tbot)
!         Qa   (i,j) = buffr(icn,cpl_fields_c2i_shum)
!         rhoa (i,j) = buffr(icn,cpl_fields_c2i_dens)

         !--- ocn states--
!         ss_tltx(i,j) = buffr(icn,cpl_fields_c2i_dhdx)
!         ss_tlty(i,j) = buffr(icn,cpl_fields_c2i_dhdy)
!         frzmlt (i,j) = buffr(icn,cpl_fields_c2i_q)

         !--- atm fluxes--
!         swvdr(i,j) = buffr(icn,cpl_fields_c2i_swvdr)
!         swidr(i,j) = buffr(icn,cpl_fields_c2i_swndr)
!         swvdf(i,j) = buffr(icn,cpl_fields_c2i_swvdf)
!         swidf(i,j) = buffr(icn,cpl_fields_c2i_swndf)
!         flw  (i,j) = buffr(icn,cpl_fields_c2i_lwdn)
!         frain(i,j) = buffr(icn,cpl_fields_c2i_rain)
!         fsnow(i,j) = buffr(icn,cpl_fields_c2i_snow)

!!======================================================
!!         get flux from FVCOM
!!======================================================
         !--- ocn states--
         sst  (i,j) = T1(J,1) !+Tffresh
         sss  (i,j) = S1(J,1)
         Tf   (i,j) = -depressT*sss(i,j) ! freezing temp (C)

!         for test case no flux to fvcom   ggao


         !--- atm states-
         zlvl (i,j) = c10i
         Tair (i,j) = T_air(J) +Tffresh  !  ^0C--->"K"
         potT (i,j) = Tair(i,j)          ! K 
         Qa   (i,j) = QA_AIR(J)
         rhoa (i,j) = 1.3_SP !rhoair

         !--- ocn states--
         ss_tltx(i,j) = utmp(j) ! sea surface slope, x-direction (m/m) 
         ss_tlty(i,j) = vtmp(j) ! sea surface slope, y-direction (m/m)

         ! compute potential to freeze or melt ice
         frzmlt(i,j) = (Tf(i,j)-sst(i,j))*cphm/dtice
         frzmlt(i,j) = (Tf(i,j)-sst(i,j))*cp_ocn*rhow*HM(I,J)/dtice
         frzmlt(i,j) = min(max(frzmlt(i,j),-c1000),c1000)

         !--- atm fluxes--
!         swvdr(i,j) = buffr(icn,cpl_fields_c2i_swvdr)
!         swidr(i,j) = buffr(icn,cpl_fields_c2i_swndr)
!         swvdf(i,j) = buffr(icn,cpl_fields_c2i_swvdf)
!         swidf(i,j) = buffr(icn,cpl_fields_c2i_swndf)

          fsw (i,j)  = DSW_AIR(J)
          flw(i,j) = stefan_boltzmann*Tair(i,j)**4 !!! downward longwave !!!
!          downward longwave calculate in prepare_forcing
          frain(i,j) = QPREC(J)
          fsnow(i,j) = QPREC(J)   
          cldf(i,j)  = CLOUD(J)
!!======================================================
!!======================================================

      end do
      end do

!-----+----------------------------------------------------------------+
!    some parameter correction
!    Qa, Fsw , SW-->SWvdr SWVDF SWidr SWidf , FLW, frain, frain

       call prepare_forcing  !scale input forcing data
!!======================================================
!!     calculate the downward longwave radiation
!!======================================================


!      call ice_timer_stop(17)  ! time spent cr-unpacking

      !-----------------------------------------------------------------
      ! broadcast dbug diagnostic level
      !-----------------------------------------------------------------
!      if (irbuf(cpl_fields_ibuf_infobug) >= 2 ) then
!         if (my_task == master_task) write (nu_diag,*)
!     &        '(from_coupler) dbug level >= 2'
!         info_dbug = 1
!      endif

      !-----------------------------------------------------------------
      ! broadcast write_restart flag
      !-----------------------------------------------------------------
!      if (irbuf(cpl_fields_ibuf_resteod) == 1 .AND. new_day) then
!         if (my_task == master_task) write (nu_diag,*)
!     &        '(from_coupler) received write restart signal'
!         write_restart = 1
!      endif

      !-----------------------------------------------------------------
      ! broadcast cpl_write_history flag
      !-----------------------------------------------------------------
!      if (irbuf(cpl_fields_ibuf_histeod) == 1 .AND. new_day) then
!         if (my_task == master_task) write (nu_diag,*)
!     &        '(from_coupler) received write history signal'
!         cpl_write_history = 1
!      endif

      !-----------------------------------------------------------------
      ! broadcast stop_now flag
      !-----------------------------------------------------------------
!      if (irbuf(cpl_fields_ibuf_stopnow) == 1) then
!         if (my_task==master_task) write (nu_diag,*)
!     &        '(from_coupler) received terminate signal'
!         stop_now = 1
!      endif

!      if (info_dbug == 1 .AND. stop_now /= 1) then

!        do j=jlo,jhi
!        do i=ilo,ihi
!           worka(i,j) = tarea(i,j)
!        enddo
!        enddo

!        do n=1,cpl_fields_c2i_total
!           work_l1 = c0i
!           n2   = 0
!           do j=jlo,jhi
!           do i=ilo,ihi
!              n2 = n2 + 1 
!              if (hm(i,j) > p5) work_l1(i,j) = buffr(n2,n)
!           enddo
!           enddo
!           call bound(work_l1)
!           call get_sum(0, worka, one, work_l1, gsum)
!           if (my_task == master_task) then
!              write (nu_diag,100) 'ice', 'recv', n, gsum
!           endif
!        enddo                   ! cpl_fields_c2i_total
!      endif
 100  format ('comm_diag',1x,a3,1x,a4,1x,i3,es26.19)

      !-----------------------------------------------------------------
      ! rotate zonal/meridional vectors to local coordinates
      ! compute data derived quantities
      !-----------------------------------------------------------------

      ! Vector fields come in on T grid, but are oriented geographically
      ! need to rotate to pop-grid FIRST using ANGLET
      ! then interpolate to the U-cell centers  (otherwise we
      ! interpolate across the pole)
      ! use ANGLET which is on the T grid !

      do j=jlo,jhi
      do i=ilo,ihi
         ! ocean
!         workx      = uocn  (i,j) ! currents, m/s 
!         worky      = vocn  (i,j)
!         uocn(i,j) = workx*cos(ANGLET(i,j))  &  ! convert to POP grid 
!                   + worky*sin(ANGLET(i,j))
!         vocn(i,j) = worky*cos(ANGLET(i,j))  &
!                   - workx*sin(ANGLET(i,j))
!
!         workx      = ss_tltx  (i,j)           ! sea sfc tilt, m/m
!         worky      = ss_tlty  (i,j)
!         ss_tltx(i,j) = workx*cos(ANGLET(i,j)) & ! convert to POP grid 
!                      + worky*sin(ANGLET(i,j))
!         ss_tlty(i,j) = worky*cos(ANGLET(i,j)) &
!                      - workx*sin(ANGLET(i,j))

!         sst(i,j) = sst(i,j) !- Tffresh         ! sea sfc temp (C)
!         Tf (i,j) = -1.8_dbl_kind              ! hardwired for NCOM
!c        Tf (i,j) = -depressT*sss(i,j)         ! freezing temp (C)
!c        Tf (i,j) = -depressT*max(sss(i,j),ice_ref_salinity)

      enddo
      enddo

      ! Interpolate ocean dynamics variables from T-cell centers to 
      ! U-cell centers.

!      call t2ugrid(uocn)
!      call t2ugrid(vocn)
!      call t2ugrid(ss_tltx)
!      call t2ugrid(ss_tlty)

      ! Atmosphere variables are needed in T cell centers in
      ! subroutine stability and are interpolated to the U grid
      ! later as necessary.

      do j=jlo,jhi
      do i=ilo,ihi
         ! atmosphere
!         workx      = uatm(i,j) ! wind velocity, m/s
!         worky      = vatm(i,j) 
!         uatm (i,j) = workx*cos(ANGLET(i,j)) ! convert to POP grid
!     &              + worky*sin(ANGLET(i,j)) ! note uatm, vatm, wind
!         vatm (i,j) = worky*cos(ANGLET(i,j)) !  are on the T-grid here
!     &              - workx*sin(ANGLET(i,j))

!         wind (i,j) = sqrt(uatm(i,j)**2 + vatm(i,j)**2) ! wind speed, m/s
!         fsw  (i,j) = swvdr(i,j) + swvdf(i,j)
!     &              + swidr(i,j) + swidf(i,j)

!!-------------------------------------------------------------------------
!!-------------------------------------------------------------------------
         ! atmosphere      fvcom forcing

        ! workx      = uatm(i,j) ! wind velocity, m/s
        ! worky      = vatm(i,j)
        ! uatm (i,j) = workx*cos(ANGLET(i,j))  & ! convert to POP grid
        !           + worky*sin(ANGLET(i,j))     ! note uatm, vatm, wind
        ! vatm (i,j) = worky*cos(ANGLET(i,j))  & !  are on the T-grid here
        !           - workx*sin(ANGLET(i,j))

        ! wind (i,j) = sqrt(uatm(i,j)**2 + vatm(i,j)**2) ! wind speed, m/s
        ! fsw  (i,j) = swvdr(i,j) + swvdf(i,j)  &
        !           + swidr(i,j) + swidf(i,j)
!       if(msr) write(107,'(5f10.2)')fsw  (i,j),vatm (i,j),uatm (i,j),Tair(i,J),flw(i,j)

!!-------------------------------------------------------------------------
!!-------------------------------------------------------------------------
      enddo
      enddo

!      time_forc=time

!      call ice_timer_stop(8)   ! time spent coupling

      end subroutine from_coupler

!=======================================================================
!BOP
!
! !IROUTINE: to_coupler - send data from sea ice model to coupler
!
! !INTERFACE:
!
      subroutine to_coupler
!
! !DESCRIPTION:
!
! Sea ice model to coupler
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
      use ice_flux
!      use ice_calendar, only: dt
!      use ice_grid, only: tmask

      use ice_atmo
      use ice_state
      use ice_albedo

      integer(kind=int_kind) :: i,j,icn,n2     ! local loop indices

      real (kind=dbl_kind) ::    &
         gsum, workx, worky      &     ! tmps for converting grid
      ,  Tsrf (ilo:ihi,jlo:jhi)  &     ! surface temperature
      ,  tauxa(ilo:ihi,jlo:jhi)  &     ! atmo/ice stress
      ,  tauya(ilo:ihi,jlo:jhi)  &             
      ,  tauxo(ilo:ihi,jlo:jhi)  &     ! ice/ocean stress
      ,  tauyo(ilo:ihi,jlo:jhi)  &             
      ,  ailohi(ilo:ihi,jlo:jhi)      ! fractional ice area
!YY+
      logical :: flag

!!=================================================================
!!=================================================================
      !!  just as ice_ocean--->couple to ocean mixing 
      real (kind=dbl_kind), dimension(ilo:ihi,jlo:jhi) ::&
         delt    & ! potential temperature difference   (K)
      ,  delq    & ! specific humidity difference   (kg/kg)
      ,  dummy1, dummy2, dummy3, dummy4  ! dummy arrays

      real (kind=dbl_kind) ::  &
         TsfK    & ! surface temperature (K)
      ,  fsh     & ! sensible heat flux  (W/m^2)
      ,  flh     & ! latent heat flux    (W/m^2)
      ,  swabs   & ! surface absorbed shortwave heat flux (W/m^2)
      ,  flwup   & ! long-wave upward heat flux  (W/m^2)
      ,  qdp     & ! deep ocean heat flux
      ,  ft        ! fraction reduction of positive qdp

!!=================================================================
      REAL(SP) :: SPRO,SPCP,ROSEA
      REAL(SP), allocatable, dimension (:,:) :: rbuf

!!=================================================================
!--------ice_ocean---------------------------------------------------
      real (kind=dbl_kind), parameter ::  &
        cphm = cp_ocn*rhow*20 !hmix
!         hmix = 20._dbl_kind         &    ! ocean mixed layer depth (m)
!--------ice_ocean---------------------------------------------------


      flag=.false.
!YY-

!      call ice_timer_start(8)  ! time spent coupling

      do j=jlo,jhi
      do i=ilo,ihi

        ! ice fraction
!        ailohi(i,j) = aice(i,j)

        ! surface temperature
!        Tsrf(i,j)  = Tffresh + Tsfc(i,j)                    !K

        ! wind stress  (on POP T-grid:  convert to lat-lon)
!        workx = strairxT(i,j)                               ! N/m^2
!        worky = strairyT(i,j)                               ! N/m^2
!        tauxa(i,j) = workx*cos(ANGLET(i,j)) - worky*sin(ANGLET(i,j))
!        tauya(i,j) = worky*cos(ANGLET(i,j)) + workx*sin(ANGLET(i,j))

        ! ice/ocean stress (on POP T-grid:  convert to lat-lon)
!        workx = -strocnxT(i,j)                              ! N/m^2
!        worky = -strocnyT(i,j)                              ! N/m^2
!        tauxo(i,j) = workx*cos(ANGLET(i,j)) - worky*sin(ANGLET(i,j))
!        tauyo(i,j) = worky*cos(ANGLET(i,j)) + workx*sin(ANGLET(i,j))
!!====================================================================
!!====================================================================
        ! ice fraction
!        ailohi(i,j) = aice(i,j)

        ! surface temperature
!        Tsrf(i,j)  = Tffresh + Tsfc(i,j)                    !K

        ! wind stress  (on POP T-grid:  convert to lat-lon)
!        workx = strairxT(i,j)                               ! N/m^2
!        worky = strairyT(i,j)                               ! N/m^2
!        tauxa(i,j) = workx*cos(ANGLET(i,j)) - worky*sin(ANGLET(i,j))
!        tauya(i,j) = worky*cos(ANGLET(i,j)) + workx*sin(ANGLET(i,j))

        ! ice/ocean stress (on POP T-grid:  convert to lat-lon)
!        workx = -strocnxT(i,j)                              ! N/m^2
!        worky = -strocnyT(i,j)                              ! N/m^2
!        tauxo(i,j) = workx*cos(ANGLET(i,j)) - worky*sin(ANGLET(i,j))
!        tauyo(i,j) = worky*cos(ANGLET(i,j)) + workx*sin(ANGLET(i,j))


!!====================================================================
!!====================================================================


      enddo
      enddo

      !--- set info buffer flags ---
!      isbuf                          = 0       ! unused
!      isbuf(cpl_fields_ibuf_stopnow) = 0       ! stop flag: 0 <=> able to continue
!      isbuf(cpl_fields_ibuf_cdate)   = idate   ! model date, coded: yyyymmdd
!      isbuf(cpl_fields_ibuf_sec)     = sec     ! elapsed seconds on model date
!      isbuf(cpl_fields_ibuf_lisize)  = (ihi-ilo+1) 
!      isbuf(cpl_fields_ibuf_ljsize)  = (jhi-jlo+1) 
!      isbuf(cpl_fields_ibuf_ncpl)    = nadv_i  ! number of msg-pairs per day
!      isbuf(cpl_fields_ibuf_lsize)   = (jhi-jlo+1)*(ihi-ilo+1)
!      isbuf(cpl_fields_ibuf_dead)    = 0       ! not a dead model

!      call ice_timer_start(18)      ! Time spent packing

      !--- pack & send msg buffer ---
!YY+
!      do j=jlo,jhi
!      do i=ilo,ihi
!         if (tmask(i,j) .and. ailohi(i,j) < c0i ) then
!            flag = .true.
!         endif
!      end do
!      end do
!      if (flag) then
!        do j=jlo,jhi
!        do i=ilo,ihi
!          if (tmask(i,j) .and. ailohi(i,j) < c0i ) then
!            write(nu_diag,*)            &
!                 ' (ice) send: ERROR ailohi < 0.0 ',i,j,ailohi(i,j)
!            call shr_sys_flush(nu_diag)
!          endif
!        end do
!        end do
!      endif
!YY-
!      buffs(:,:)=spval
      icn=0
      do j=jlo,jhi
      do i=ilo,ihi
         icn=icn+1
!YY         if (tmask(i,j) .and. ailohi(i,j) < c0i ) then
!YY            write(nu_diag,*)
!YY     &           ' (ice) send: ERROR ailohi < 0.0 ',i,j,ailohi(i,j)
!YY            call shr_sys_flush(nu_diag)
!YY         endif

         !--- ice states
!         buffs(in,cpl_fields_i2c_ifrac) = ailohi(i,j)    ! frac 

!!====================================================================
!!====================================================================
!YY         if (tmask(i,j) .and. ailohi(i,j) < c0i ) then
!YY            write(nu_diag,*)
!YY     &           ' (ice) send: ERROR ailohi < 0.0 ',i,j,ailohi(i,j)
!YY            call shr_sys_flush(nu_diag)
!YY         endif

         !--- ice states
!         buffs(in,cpl_fields_i2c_ifrac) = ailohi(i,j)    ! frac

!!====================================================================
!!====================================================================


!         endif  ! tmask and ailohi > c0i
      end do
      end do


 100  format('comm_diag',1x,a3,1x,a4,1x,i3,es26.19)

      call atmo_boundary_layer (1, 'ocn', sst,   &
          dummy1, dummy2, dummy3, dummy4, delt,  delq)

      do j = jlo,jhi
      do i = ilo,ihi
!       if (tmask(i,j)) then

!       if (aice(i,j) > puny) then

         ! specify deep ocean heat flux as constant for now
!         qdp = -c2i              ! negative upward
!
!         ! ocean surface temperature in Kelvin
!         TsfK = sst(i,j) + Tffresh
!
!         ! shortwave radiative flux
!         swabs = (c1i - albocn)    &
!             * (swvdr(i,j) + swvdf(i,j) + swidr(i,j) + swidf(i,j))
!
!         ! longwave radiative flux
!         flwup  = -emissivity*stefan_boltzmann * TsfK**4
!
!         ! downward latent and sensible heat fluxes
!         flh = lhcoef(i,j) * delq(i,j)
!         fsh = shcoef(i,j) * delt(i,j)
!
!         ! first, compute sst change due to exchange with atm/ice above
!         sst(i,j) = sst(i,j) +                         &
!              (fsh + flh + flwup + flw(i,j) + swabs)  &
!              * (c1i-aice(i,j)) * dtice / cphm
!
!         ! computed T change due to exchange with deep layers:
!         sst(i,j) = sst(i,j) - qdp*dtice/cphm
!
!         ! compute potential to freeze or melt ice
!         frzmlt(i,j) = (Tf(i,j)-sst(i,j))*cphm/dtice
!         frzmlt(i,j) = min(max(frzmlt(i,j),-c1000),c1000)
!
!         ! if sst is below freezing, reset sst to Tf
!         if (sst(i,j) <= Tf(i,j)) then
!            sst(i,j) = Tf(i,j)
!            T1(j,1) =sst(i,j)
!         end if
!
!       endif                    ! tmask
!       endif                    ! tmask

       !!!==============================================================
       !!!==============================================================
           !!  heat for the ocean  model
!           WTSURF(J) = (fsh + flh + flwup + flw(i,j))*(c1i-aice(i,j))
!           SWRAD(J)  = swabs*(c1i-aice(i,j))


       !!!==============================================================
       !!!==============================================================
!           WTSURF(J) = -WTSURF(J)/SPRO  !*RAMP
!           SWRAD(J)  = -SWRAD(J)/SPRO   !*RAMP

           WTSURF(J) = WTSURF(J)*(c1i-aice(i,j))  !*RAMP
           SWRAD(J)  = SWRAD(J)*(c1i-aice(i,j))   !*RAMP

           if(sst(i,j)<=Tf(i,j)) then         
           sst(i,j)=Tf(i,j)
           T1(j,1)= sst(i,j)
           end if

      enddo                     ! i
      enddo                     ! j

       !!!==============================================================
       !!!==============================================================
       !!! modify the fresh water and salinity exchange between the ice and ocean
       !        reculate(reset) the heat flux and fresh water for the ocean
       !===============================================================!!
       !     &,   freshn      ! fresh water flux to ocean       (kg/m2/s)
       !     &,   fsaltn      ! salt flux to ocean              (kg/m2/s)
       !         fresh     (i,j) = fresh     (i,j) + freshn(i,j)  * aicen(i,j,n)
       !         fsalt     (i,j) = fsalt     (i,j) + fsaltn(i,j)  * aicen(i,j,n)

       DO J=1,MT
!         QPREC3(J)=frain(1,J)+fresh(1,J)  !  precipation and fresh from ice melt
!         QEVAP3(J)=QEVAP3(J) -fsalt(1,j)/1000./sss(1,J)
         !   convert the reject salt  ---as evap
!       if(MSR)write(201,*)frain(1,J),fresh(1,J),fsnow(1,j),QEVAP3(J)
       END DO

       !!!==============================================================
       !!!==============================================================
          !! momentum for the ocean (wind stress and ice drag)
          !! modify the wind stress over the ocean according to the ice concentration
          !!------------------------------------------------------------
          !! Update to new time

          !! momentum for the ocean (wind stress and ice drag)
          !! modify the wind stress over the ocean according to the ice concentration
!          allocate(rbuf(nne,2))
!          rbuf(1:N,2) =  cbcice(1:N) * &
!               sqrt((U(1:N,1)-uice2(1:N))**2+(V(1:N,1)-vice2(1:N))**2)
!          rbuf(1:N,1) = rbuf(1:N,2)*(U(1:N,1)-fuicec(1:N))
!          rbuf(1:N,2) = rbuf(1:N,2)*(V(1:N,1)-fvicec(1:N))
!!
!          wusurf(1:N)=wusurf(1:N)*(1.0_SP-faicec)+rbuf(1:N,1)
!          wvsurf(1:N)=wvsurf(1:N)*(1.0_SP-faicec)+rbuf(1:N,2)
!          wusurf2(1:N)=wusurf2(1:N)*(1.0_SP-faicec)+rbuf(1:N,1)
!          wvsurf2(1:N)=wvsurf2(1:N)*(1.0_SP-faicec)+rbuf(1:N,2)
          DO I=1,N
!          write(205,'(4f20.8)')wusurf(I),wusurf(I),wusurf2(I),wusurf2(I)
          END DO
!          IF(PAR)CALL EXCHANGE(EC,NT,1,MYID,NPROCS,WUSURF,WVSURF)
!          IF(PAR)CALL EXCHANGE(EC,NT,1,MYID,NPROCS,WUSURF2,WVSURF2)

          deallocate(rbuf)

       !!!==============================================================

      end subroutine to_coupler

!=======================================================================
!BOP
!
! !IROUTINE: exit_coupler - exit from coupled/mpi environment
!
! !INTERFACE:
!
      subroutine exit_coupler
!
! !DESCRIPTION:
!
! Exit from coupled/MPI environment
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
!      include "mpif.h"         ! MPI library definitions

      integer (kind=int_kind) ::  &
        ierr  ! error flag

!      if (my_task == master_task) then
!         if (irbuf(cpl_fields_ibuf_stopnow) == 1) then
!            write (nu_diag,*) '(ice) received final coupler msg'
!         else
!            write (nu_diag,*) '(ice) terminating before coupler'
!            call MPI_ABORT(MPI_FVCOM_GROUP,-1,ierr)
!         endif
!      endif

!      call cpl_interface_finalize (cpl_fields_icename)

!      if (my_task == master_task) then
!         write(nu_diag,*) '(ice) exit_coupler finished',my_task
!      endif

      end subroutine exit_coupler

!=======================================================================

!#endif                          ! coupled

      end module ice_coupling

!=======================================================================
