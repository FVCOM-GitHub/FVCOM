!/===========================================================================/
! CVS VERSION INFORMATION
! $Id$
! $Name$
! $Revision$
!/===========================================================================/

!=======================================================================
!BOP
!
! !MODULE:   ice_init - parameter and variable initializations
!
! !DESCRIPTION:
! 
! parameter and variable initializations
!
! !REVISION HISTORY:
! 
! authors Elizabeth C. Hunke, LANL
!         C. M. Bitz, UW
!
! !INTERFACE:
!
      module ice_init
!
! !USES:
!
      use ice_domain
!
!EOP
!
      implicit none
      save

      character (len=char_len) :: & 
         advection   ! type of advection algorithm used
                     ! 'upwind'  => 1st order mpdata scheme (donor cell)
                     ! 'mpdata'  => 2nd order mpdata scheme
                     ! 'remap' (or anything else) => remapping scheme      

      character(len=char_len) :: & 
         ice_ic      ! method of ice cover initialization
                     ! 'default'  => latitude and sst dependent
                     ! 'none'     => no ice 
                     ! note:  restart = .true. overwrites

!=======================================================================

      contains

!=======================================================================
!BOP
!
! !IROUTINE: input_data - namelist variables
!
! !INTERFACE:
!
      subroutine input_data
!
! !DESCRIPTION:
!
! Namelist variables, set to default values; may be altered
! at run time
!
! !REVISION HISTORY:
! 
! author Elizabeth C. Hunke, LANL
!
! !USES:
!
      use ice_albedo
      use ice_mechred, only: kstrength, krdg_partic, krdg_redist
!      use ice_diagnostics
!      use ice_history
      use ice_calendar
!      use ice_dyn_evp
      use ice_itd, only: kitd, kcatbound
!      use ice_ocean, only: oceanmixed_ice
      use ice_flux_in, only:                            &
           dbug, ycycle, fyear_init, precip_units,      &
           atm_data_type, atm_data_dir,                 &
           sss_data_type, sst_data_type, ocn_data_dir,  &
           oceanmixed_file, restore_sst, trestore

       use ice_grid

       use ice_fileunits
       use control
       use ice_therm_vertical, only :  i0vis, floediam 
       use all_vars, only : MSR
!      use ice_exit
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) :: & 
         ni       &  ! index for processor count
      ,  nml_error ! namelist i/o error flag

      character (len=6) :: chartmp

!#ifdef CCSM
      !-----------------------------------------------------------------
      ! Declare namelist variables used for CCSM coupled runs and not
      ! declared elsewhere in CICE
      !-----------------------------------------------------------------

      integer (kind=int_kind) :: kcolumn  ! 1 for column model

      character(len=char_len_long) :: & 
         runid                 ! identifier for coupled run
!#endif
 
      integer(kind=int_kind) :: ios, i
      character(Len=120):: FNAME
 
      !-----------------------------------------------------------------
      ! Namelist variables.
      ! NOTE: Not all of these are used by both models.
      !-----------------------------------------------------------------

      namelist/nml_ice/ kstrength,  krdg_partic,      krdg_redist      & 
      , albicev,        albicei,         albsnowv,      albsnowi       &
      , kitd,           kcatbound,          i0vis,      floediam     

!!     &  year_init,      istep0,          dt,            npt
!        year_init,      istep0,                     npt                 &
!!      , diagfreq,       print_points,    print_global,  diag_type         &
!!      , diag_file                                                         &
!!      , histfreq,       hist_avg,        history_dir,   history_file      &
!!      , dumpfreq,       dumpfreq_n,      dump_file                        &
!!      , restart,        restart_dir,     pointer_file,  ice_ic            &
!!      , grid_type,      grid_file,       kmt_file                         &
!      , grid_type                                 &
!      , kitd,           kcatbound                                      &
!!      , kdyn,           ndyn_dt,         ndte                          &
!      , kstrength,      krdg_partic,     krdg_redist,   advection      &
!      , albicev,        albicei,         albsnowv,      albsnowi       &
!      , fyear_init,     ycycle                                         &
!!      , atm_data_type,  atm_data_dir,    precip_units                     &
!!      , oceanmixed_ice, sss_data_type,   sst_data_type                    &
!!      , ocn_data_dir,   oceanmixed_file, restore_sst,   trestore          &
!!      , dbug                                                               
!!#ifdef CCSM
!! These variables are used in CCSM, but not CICE
!!     &, runid,           runtype,             kcolumn
!!     , runid,           runtype,      
!      , kcolumn
!!#endif



      !-----------------------------------------------------------------
      ! default values
      !-----------------------------------------------------------------

      year_init = 1      ! initial year
      istep0 = 0         ! number of steps taken in previous integrations,
                         ! real (dumped) or imagined (use to set calendar)
!      dt = 3600.0_dbl_kind        ! time step, s 
      dtice = 3600.0_dbl_kind        ! time step, s 
      npt = 99999        ! total number of time steps (dt) 
!      diagfreq = 24      ! how often diag output is written
!      print_points = .false.      ! if true, print point data
!      print_global = .true.       ! if true, print global diagnostic data
!      diag_type = 'stdout'        ! 'file' writes to diag_file
!      diag_file = 'ice_diag.d'    ! if diag_type /= 'file'
!      histfreq='m'       ! output frequency
!      hist_avg = .true.  ! if true, write time-averages rather than snapshots
!      history_dir  = ' '          ! default = executable directory
!      history_file = 'iceh'       ! history file name prefix
!      dumpfreq='y'       ! restart frequency option
!      dumpfreq_n= 1      ! restart frequency
!      dump_file    = 'iced'       ! restart file name prefix
!      restart = .false.  ! if true, read restart files for initialization
!      restart_dir  = ' ' ! default = executable directory
!      pointer_file = 'ice.restart_file'
!      ice_ic       = 'default'     ! latitude and sst-dependent
!      grid_type    = 'rectangular' ! define rectangular grid internally
!      grid_file    = 'unknown_grid_file'
!      kmt_file     = 'unknown_kmt_file'
       kitd = 1           ! type of itd conversions (0 = delta, 1 = linear)
       kcatbound     = 1  ! category boundary formula (0 = old, 1 = new)
!      kdyn = 1           ! type of dynamics (1 = evp, 0 = off)
      ndyn_dt = 2        ! dynamics subcycles per thermodynamics timestep

!      ndte = 120         ! subcycles per dynamics timestep:  ndte=dyn_dt/dte
!      evp_damping = .false.       ! if true, use damping procedure in evp dynamics
# if defined (ICE_FRESHWATER)
! afm 20160717 & EJA 20160921 
!classical Hibler's P formula produces better thickness for the Great Lakes
      kstrength = 0      ! 1 = Rothrock 75 strength, 0 = Hibler 79 strength
# else
      kstrength = 1      ! 1 = Rothrock 75 strength, 0 = Hibler 79 strength
# endif
      krdg_partic   = 1  ! 1 = LH06 participation, 0 = Thorndike et al 75
      krdg_redist   = 1  ! 1 = LH06 redistribution, 0 = Hibler 80
!      advection  = 'remap'        ! incremental remapping transport scheme
      albicev   = 0.78_dbl_kind   ! visible ice albedo for h > ahmax
      albicei   = 0.36_dbl_kind   ! near-ir ice albedo for h > ahmax
      albsnowv  = 0.98_dbl_kind   ! cold snow albedo, visible
      albsnowi  = 0.70_dbl_kind   ! cold snow albedo, near IR
      fyear_init = 1900           ! first year of forcing cycle
      ycycle        = 1           ! number of years in forcing cycle
!      atm_data_type = 'default'   ! see ice_flux_in.F for other options
!      atm_data_dir  = ' '
!      precip_units  = 'mks'       ! 'mm_per_month' or
                                  ! 'mm_per_sec' = 'mks' = kg/m^2 s
!      oceanmixed_ice = .false.    ! if true, use internal ocean mixed layer
!      sss_data_type = 'default'
!      sst_data_type = 'default'
!      ocn_data_dir  = ' '
!      oceanmixed_file = 'unknown_oceanmixed_file' ! ocean forcing data
!      restore_sst   = .false.     ! restore sst if true
!      trestore      = 90          ! restoring timescale, days (0 instantaneous)
!      dbug          = .false.     ! true writes diagnostics for input forcing

!#ifdef CCSM
!      ! The following are in CCSM but not CICE.
!      ! Note: CCSM acts as if print_global = .true., kcatbound = 0
!      runid   = 'unknown'   ! default run ID
!      runtype = 'unknown'   ! default runtype
      kcolumn = 1           ! 1 = column model
!#endif

      i0vis   = 0.70_dbl_kind  ! fraction of penetrating solar rad (visible)
      floediam = 300.0_dbl_kind   ! effective floe diameter (m)

      !-----------------------------------------------------------------
      ! read from input file
      !-----------------------------------------------------------------


      ! some of the above are overwritten 
      ! after reading the namelist file 
      ios = 0
      FNAME = "./"//trim(casename)//"_run.nml"
      CALL FOPEN(NMLUNIT,trim(FNAME),'cfr')
      READ(UNIT=NMLUNIT, NML=nml_ice,IOSTAT=ios)
      if(MSR) write(UNIT=IPT,NML=NML_ICE)
      if(ios .NE. 0 ) then
       Call Fatal_Error("Can Not Read NameList NML_ICE from file: "//trim(FNAME))
      end if
      REWIND(NMLUNIT)
      CLOSE(NMLUNIT)


!      if (my_task == master_task) then
!        open (nu_nml, file='ice_in', status='old')
!   10   continue  !*** keep reading until right namelist is found
!        read(nu_nml, nml=ice_nml, iostat=nml_error)
!        if (nml_error > 0) goto 10 ! An error occurred
!        if (nml_error < 0) goto 20 ! End of file condition
!        close(nu_nml)
!   20   continue
!      endif
!      call ice_bcast_iscalar(nml_error)

!      if (nml_error /= 0) then
!         call abort_ice ('ice: Namelist read error in ice_init.F')
!      endif

!      if (trim(diag_type) == 'file') then
!         nu_diag = 48
!      else
!         nu_diag = 6
!      endif

!      if (histfreq == '1') hist_avg = .false. ! potential conflict
!      chartmp = advection(1:6)
!      if (chartmp /= 'upwind' .and. chartmp /= 'mpdata')
!     &     advection = 'remap'

!#ifdef CCSM
      if (kcolumn == 1) grid_type    = 'column'
!#endif

      !-----------------------------------------------------------------
      ! broadcast to all processors
      !-----------------------------------------------------------------

!      call ice_bcast_iscalar(year_init)
!      call ice_bcast_iscalar(istep0)
!      call ice_bcast_rscalar(dt)
!      call ice_bcast_iscalar(npt)
!      call ice_bcast_iscalar(diagfreq)
!      call ice_bcast_logical(print_points)
!      call ice_bcast_logical(print_global)
!      call ice_bcast_char   (diag_type)
!      call ice_bcast_char   (diag_file)
!!      call ice_bcast_iscalar(nu_diag) ! only master_task writes to file
!      call ice_bcast_char   (histfreq)
!      call ice_bcast_logical(hist_avg)
!      call ice_bcast_char   (history_dir)
!      call ice_bcast_char   (history_file)
!      call ice_bcast_char   (dumpfreq)
!      call ice_bcast_iscalar(dumpfreq_n)
!      call ice_bcast_char   (dump_file)
!      call ice_bcast_logical(restart)
!      call ice_bcast_char   (restart_dir)
!      call ice_bcast_char   (pointer_file)
!      call ice_bcast_char   (ice_ic)
!      call ice_bcast_char   (grid_type)
!      call ice_bcast_char   (grid_file)
!      call ice_bcast_char   (kmt_file)
!      call ice_bcast_iscalar(kitd)
!      call ice_bcast_iscalar(kcatbound)
!      call ice_bcast_iscalar(kdyn)
!      call ice_bcast_iscalar(ndyn_dt)
!      call ice_bcast_iscalar(ndte)
!      call ice_bcast_logical(evp_damping)
!      call ice_bcast_iscalar(kstrength)
!      call ice_bcast_iscalar(krdg_partic)
!      call ice_bcast_iscalar(krdg_redist)
!      call ice_bcast_char   (advection)
!      call ice_bcast_rscalar(albicev)
!      call ice_bcast_rscalar(albicei)
!      call ice_bcast_rscalar(albsnowv)
!      call ice_bcast_rscalar(albsnowi)
!      call ice_bcast_iscalar(fyear_init)
!      call ice_bcast_iscalar(ycycle)
!      call ice_bcast_char   (atm_data_type)
!      call ice_bcast_char   (atm_data_dir)
!      call ice_bcast_char   (precip_units)
!      call ice_bcast_logical(oceanmixed_ice)
!      call ice_bcast_char   (sss_data_type)
!      call ice_bcast_char   (sst_data_type)
!      call ice_bcast_char   (ocn_data_dir)
!      call ice_bcast_char   (oceanmixed_file)
!      call ice_bcast_logical(restore_sst)
!      call ice_bcast_iscalar(trestore)
!      call ice_bcast_logical(dbug)

!#ifdef CCSM
!      ! The following are in CCSM but not CICE.
!      call ice_bcast_char   (runid)
!      call ice_bcast_char   (runtype)
!      call ice_bcast_iscalar(kcolumn)
!#endif
      !-----------------------------------------------------------------
      ! write namelist variables to diagnostic file
      !-----------------------------------------------------------------

!      if (my_task == master_task) then
!         if (trim(diag_type) == 'file') then
!            write(6,*) 'Diagnostic output will be in file ', diag_file
!            open (nu_diag, file=diag_file, status='unknown')
!         endif

!         write(nu_diag,*) ' '
!         write(nu_diag,*) '--------------------------------'
!         write(nu_diag,*) '  CICE model diagnostic output  '
!         write(nu_diag,*) '--------------------------------'
!         write(nu_diag,*) ' '
!#ifdef CCSM
!         if (trim(runid) /= 'unknown')
!     &    write(nu_diag,*)    ' runid                     = ',
!     &                         trim(runid)
!         if (trim(runtype) /= 'unknown')
!     &    write(nu_diag,*) ' runtype                   = ',
!     &                         trim(runtype)
!#endif
!         write(nu_diag,*) ' year_init                 = ', year_init
!         write(nu_diag,*) ' istep0                    = ', istep0
!         write(nu_diag,*) ' dt                        = ', dt
!         write(nu_diag,*) ' npt                       = ', npt
!         write(nu_diag,*) ' diagfreq                  = ', diagfreq
!         write(nu_diag,*) ' print_global               = ',             &
!                               print_global                              
!         write(nu_diag,*) ' print_points              = ',              &
!                               print_points                              
!         write(nu_diag,*) ' histfreq                  = ',              &
!                               trim(histfreq)                            
!         write(nu_diag,*) ' hist_avg                  = ', hist_avg      
!         if (hist_avg) then                                              
!           write (nu_diag,*) ' History data will be averaged over 1 ',  &
!                               histfreq                                  
!         else                                                            
!           write (nu_diag,*) ' history data will be snapshots'           
!         endif                                                           
!         write(nu_diag,*)    ' history_dir               = ',           &
!                               trim(history_dir)                         
!         write(nu_diag,*)    ' history_file              = ',           &
!                               trim(history_file)                        
!         write(nu_diag,*) ' dumpfreq                  = ',              &
!                               trim(dumpfreq)                            
!         write(nu_diag,*) ' dumpfreq_n                = ', dumpfreq_n    
!         write(nu_diag,*)    ' dump_file                 = ',           &
!                               trim(dump_file)                           
!         write(nu_diag,*) ' restart                   = ', restart       
!         write(nu_diag,*)    ' restart_dir               = ',           &
!                               trim(restart_dir)                         
!         write(nu_diag,*)    ' pointer_file              = ',           &
!                               trim(pointer_file)                        
!         write(nu_diag,*)    ' ice_ic                    = ', ice_ic     
!         write(nu_diag,*) ' grid_type                 = ',              &
!                               trim(grid_type)                           
!         if (trim(grid_type) /= 'rectangular' .or.                      &
!             trim(grid_type) /= 'column') then                           
!         write(nu_diag,*) ' grid_file                 = ',              &
!                                  trim(grid_file)                        
!         write(nu_diag,*) ' kmt_file                  = ',              &
!                                  trim(kmt_file)                         
!         endif                                                           
!         write(nu_diag,*) ' kitd                      = ', kitd          
!         write(nu_diag,*) ' kcatbound                 = ',              &
!                               kcatbound                                 
!         write(nu_diag,*) ' kdyn                      = ', kdyn          
!         write(nu_diag,*) ' ndyn_dt                   = ', ndyn_dt       
!         write(nu_diag,*) ' ndte                      = ', ndte          
!         write(nu_diag,*) ' evp_damping               = ',              &
!                               evp_damping                               
!         write(nu_diag,*) ' kstrength                 = ', kstrength     
!         write(nu_diag,*) ' krdg_partic               = ',              &
!                               krdg_partic                               
!         write(nu_diag,*) ' krdg_redist               = ',              &
!                               krdg_redist                               
!         write(nu_diag,*) ' advection                 = ',              &
!                               trim(advection)                           
!         write(nu_diag,*) ' albicev                   = ', albicev       
!         write(nu_diag,*) ' albicei                   = ', albicei       
!         write(nu_diag,*) ' albsnowv                  = ', albsnowv      
!         write(nu_diag,*) ' albsnowi                  = ', albsnowi      
!         write(nu_diag,*) ' fyear_init                = ',              &
!                               fyear_init                                
!         write(nu_diag,*) ' ycycle                    = ', ycycle        
!         write(nu_diag,*)    ' atm_data_type             = ',           &
!                               trim(atm_data_type)                       
!         if (trim(atm_data_type) /= 'default') then                      
!         write(nu_diag,*)    ' atm_data_dir              = ',           &
!                               trim(atm_data_dir)                        
!         write(nu_diag,*)    ' precip_units              = ',           &
!                               trim(precip_units)                        
!         endif                                                           
!         write(nu_diag,*) ' oceanmixed_ice            = ',              &
!                               oceanmixed_ice                            
!         write(nu_diag,*)    ' sss_data_type             = ',           &
!                               trim(sss_data_type)                       
!         write(nu_diag,*)    ' sst_data_type             = ',           &
!                               trim(sst_data_type)                       
!         if (trim(sss_data_type) /= 'default' .or.                      &
!             trim(sst_data_type) /= 'default') then                      
!         write(nu_diag,*)    ' ocn_data_dir              = ',           &
!                               trim(ocn_data_dir)                        
!         endif                                                           
!         if (trim(sss_data_type) == 'ncar' .or.                         &
!             trim(sst_data_type) == 'ncar') then                         
!         write(nu_diag,*)    ' oceanmixed_file           = ',           &
!                                  trim(oceanmixed_file)
!         endif
#ifdef coupled
!         if( oceanmixed_ice ) then
!            write (nu_diag,*) 'WARNING WARNING WARNING WARNING '
!            write (nu_diag,*) '*Coupled and oceanmixed flags are  *'
!            write (nu_diag,*) '*BOTH ON.  Ocean data received from*'
!            write (nu_diag,*) '*coupler will be altered by mixed  *'
!            write (nu_diag,*) '*layer routine!                    *'
!            write (nu_diag,*) ' '
!         endif
#endif

#ifdef CCSM
!            write(nu_diag,*) ' kcolumn                   = ', kcolumn
#endif

!        if (grid_type  /=  'displaced_pole' .and.
!     &      grid_type  /=  'column' .and.
!     &      grid_type  /=  'rectangular') then
!          call abort_ice('ice_init: unknown grid_type')
!        endif

      !-----------------------------------------------------------------
      ! Document grid and subdomain sizes
      !-----------------------------------------------------------------

!        write(nu_diag,*) ' '
!        write(nu_diag,*) ' Grid and subdomain sizes:'
!        write(nu_diag,*) ' ------------------------ '
!        write(nu_diag,*) ' '
!        write(nu_diag,1050) imt_global, jmt_global
!        write(nu_diag,1060) nproc_s, nproc_x, nproc_y
!        write(nu_diag,1070) imt_local,jmt_local
!        write(nu_diag,1080) ihi-ilo+1,jhi-jlo+1
!        write(nu_diag,1090) ilo,jlo
!        write(nu_diag,1095) ihi,jhi
!        write(nu_diag,*) ' Global i start for each processor: ', &
!                           (local_start(1,n),n=1,nproc_s)
!        write(nu_diag,*) ' Global j start for each processor: ', &
!                           (local_start(2,n),n=1,nproc_s)
!        write(nu_diag,*) ' '

 1050 format(' Global problem size:',2x,i6,' x ',i6)
 1060 format(' Using ',i6,' processors in a ',i6,' x ',i6,      &
                 ' Cartesian decomposition')
 1070 format(' Local array size is:',2x,i6,' x ',i6)
 1080 format(' Physical domain is (approximately):',2x,i6,' x ',i6)
 1090 format(' Local i,j start for each processor:',2x,i6,2x,i6)
 1095 format(' Local i,j end   for each processor:',2x,i6,2x,i6)

!      endif  ! my_task = master_task

      end subroutine input_data

!=======================================================================
!BOP
!
! !IROUTINE: init_state - initialize state for itd 
!
! !INTERFACE:
!
      subroutine init_state
!
! !DESCRIPTION:
!
! Initialize state for the itd model
!
! !REVISION HISTORY:
! 
! author C. M. Bitz
!
! !USES:
!
      use ice_model_size
      use ice_constants
      use ice_flux
      use ice_therm_vertical
      use ice_grid
      use ice_state
      use ice_itd
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) :: &
         i, j           & ! horizontal indices
      ,  ij             & ! horizontal index, combines i and j loops
      ,  k              & ! ice layer index
      ,  ni              & ! thickness category index
      ,  icells           ! number of cells initialized with ice

      integer (kind=int_kind), dimension(1:(ihi-ilo+1)*(jhi-jlo+1)) :: & 
          indxi, indxj    ! compressed indices for cells with aicen > puny

      real (kind=dbl_kind) :: &
         slope, Ti, sum, hbar  &
      ,  ainit(ncat)           &
      ,  hinit(ncat)           

      real (kind=dbl_kind), parameter :: &
         hsno_init = 0.20_dbl_kind  ! initial snow thickness (m)

      if (trim(ice_ic) == 'none') then
      ! Initialize grid with no ice.
      ! If restarting, these values are overwritten.

      do ni = 1,ncat
         do j = 1,jmt_local
         do i = 1,imt_local
            aicen(i,j,ni) = c0i
            vicen(i,j,ni) = c0i
            vsnon(i,j,ni) = c0i
            esnon(i,j,ni) = c0i
         enddo
         enddo
         do j = jlo,jhi
         do i = ilo,ihi
            Tsfcn(i,j,ni) = Tf(i,j)  ! Tf not defined for ghost cells
         enddo
         enddo
!         call bound(Tsfcn(:,:,n))
      enddo

      do k = 1,ntilay
         do j = 1,jmt_local
         do i = 1,imt_local
            eicen(i,j,k) = c0i
         enddo
         enddo
      enddo

      else ! ice_ic = 'default'

      ! initial category areas in cells with ice
      hbar = c3i  ! initial ice thickness with greatest area
                 ! Note: the resulting average ice thickness 
                 ! tends to be less than hbar due to the
                 ! nonlinear distribution of ice thicknesses 
      sum = c0i
      do ni = 1, ncat
         if (ni < ncat) then
            hinit(ni) = p5*(hin_max(ni-1) + hin_max(ni)) ! m
         else                ! n=ncat
            hinit(ni) = (hin_max(ni-1) + c1i) ! m
         endif
         ! parabola, max at h=hbar, zero at h=0, 2*hbar
         ainit(ni) = max(c0i, (c2i*hbar*hinit(ni) - hinit(ni)**2))
         sum = sum + ainit(ni)
      enddo
      do ni = 1, ncat
         ainit(ni) = ainit(ni) / (sum + puny/ncat) ! normalize
      enddo

      ! place ice at high latitudes where ocean sfc is cold
      icells = 0
      do j = jlo,jhi
      do i = ilo,ihi
        if (tmask(i,j)) then
          if ((sst (i,j) <= Tf(i,j)+p2) .and.             &
              (ULAT(i,j) < -64.0_dbl_kind/rad_to_deg .or. &
!               ULAT(i,j) >  70.0_dbl_kind/rad_to_deg)) then
                ULAT(i,j) >  65.0_dbl_kind/rad_to_deg)) then
             icells = icells + 1
             indxi(icells) = i
             indxj(icells) = j
          endif          ! cold surface
        endif            ! tmask
      enddo              ! i
      enddo              ! j
      do ni = 1,ncat

         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)

            aicen(i,j,ni) = ainit(ni)
            vicen(i,j,ni) = hinit(ni) * ainit(ni) ! m
            vsnon(i,j,ni) = min(aicen(i,j,ni)*hsno_init, p2*vicen(i,j,ni))
            Tsfcn(i,j,ni) = min(Tsmelt, Tair(i,j) - Tffresh) ! deg C

            ! snow
            Ti = min(c0i, Tsfcn(i,j,ni))
            esnon(i,j,ni) = -rhos*(Lfresh - cp_ice*Ti)*vsnon(i,j,ni)
         enddo                  ! ij

         do k = 1, nilyr
            do ij = 1, icells
               i = indxi(ij)
               j = indxj(ij)

               ! assume linear temp profile and compute enthalpy
               slope = Tf(i,j) - Tsfcn(i,j,ni)
               Ti = Tsfcn(i,j,ni)                       & 
                  + slope*(real(k,kind=dbl_kind)-p5)   &
                  /real(nilyr,kind=dbl_kind)
      
               eicen(i,j,ilyr1(ni)+k-1) =              &
                    -(rhoi * (cp_ice*(Tmlt(k)-Ti)     &
                    + Lfresh*(c1i-Tmlt(k)/Ti) - cp_ocn*Tmlt(k))) &
                    * vicen(i,j,ni)/real(nilyr,kind=dbl_kind)

            enddo               ! ij
         enddo                  ! k
      enddo                     ! n

      endif ! ice_ic

      ! compute aggregate ice state and open water area
      call aggregate
      call bound_aggregate

      do j = 1, jmt_local
        do i = 1, imt_local
          aice_init(i,j) = aice(i,j)
        enddo
      enddo

      end subroutine init_state

!=======================================================================
!BOP
!
! !IROUTINE: init_flux - initialize fluxes exchanged with coupler
!
! !INTERFACE:
!
      subroutine init_flux
!
! !DESCRIPTION:
!
! Initialize all fluxes exchanged with flux coupler
! and some data derived fields
!
! !REVISION HISTORY:
! 
! author Elizabeth C. Hunke, LANL
!
! !USES:
!
      use ice_constants
      use ice_flux
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer i,j

      do j=jlo,jhi
      do i=ilo,ihi
      !-----------------------------------------------------------------
      ! fluxes received
      !-----------------------------------------------------------------
        zlvl  (i,j) = c10i             ! atm level height (m)
        uatm  (i,j) = c0i              ! wind velocity    (m/s)
        vatm  (i,j) = c0i
        potT  (i,j) = 273._dbl_kind   ! air potential temperature  (K)
        Tair  (i,j) = 273._dbl_kind   ! air temperature  (K)
        Qa    (i,j) = 0.014_dbl_kind  ! specific humidity (kg/kg)
        rhoa  (i,j) = 1.3_dbl_kind    ! air density (kg/m^3)
        fsnow (i,j) = 3.3e-6_dbl_kind ! snowfall rate (kg/m2/s)
        frain (i,j) = c0i              ! rainfall rate (kg/m2/s)
        fsw   (i,j) = c0i              ! shortwave radiation (W/m^2)
        swvdr (i,j) = c0i              ! shortwave radiation (W/m^2)
        swvdf (i,j) = c0i              ! shortwave radiation (W/m^2)
        swidr (i,j) = c0i              ! shortwave radiation (W/m^2)
        swidf (i,j) = c0i              ! shortwave radiation (W/m^2)
        flw   (i,j) = 280.0_dbl_kind  ! incoming longwave radiation (W/m^2)
        sss   (i,j) = 34.0_dbl_kind   ! sea surface salinity (o/oo)
!        uocn  (i,j) = c0i              ! surface ocean currents (m/s)
!        vocn  (i,j) = c0i
        frzmlt(i,j) = c0i              ! freezing/melting potential (W/m^2)
!        qdp   (i,j) = c0i              ! deep ocean heat flux

      !-----------------------------------------------------------------
      ! derived or computed fields
      !-----------------------------------------------------------------

        Tf      (i,j) = -depressT*sss(i,j) ! freezing temp (C)
        sst     (i,j) = Tf(i,j)            ! sea surface temp (C)

        wind    (i,j) = sqrt(uatm(i,j)**2 + vatm(i,j)**2) ! wind speed, (m/s)

        strocnx (i,j) = c0i        ! ice-ocean stress, x-direction (U-cell)
        strocny (i,j) = c0i        ! ice-ocean stress, y-direction (U-cell)
        strocnxT(i,j) = c0i        ! ice-ocean stress, x-direction (T-cell)
        strocnyT(i,j) = c0i        ! ice-ocean stress, y-direction (T-cell)

      enddo
      enddo

      call init_flux_atm
      call init_flux_ocn

      end subroutine init_flux

!=======================================================================
!BOP
!
! !IROUTINE: setup_mpi - initialize mpi
!
! !INTERFACE:
!
!      subroutine setup_mpi
!
! !DESCRIPTION:
!
! This routine initializes mpi for either internal parallel
! processing or for message passing with the coupler
!
! !REVISION HISTORY:
! 
! author Elizabeth C. Hunke, LANL
! code originally based on POP routine
!
! !USES:
!
!      use ice_mpi_internal
!      use ice_coupling
!      use ice_exit
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
!      integer (kind=int_kind) :: &
!        coords1, coords2, n, ilen, jlen &
!      , interior_i, interior_j  ! dummies for interior blocks

!      master_task = 0

!#ifdef coupled
!      !  if running in coupled mode
!#ifdef fcd_coupled
!      ! direct ice-ocean coupling
!      call MPI_COMM_DUP(MPI_FVCOM_GROUP, MPI_COMM_ICE, ierr)
!#else
!      ! CCSM coupling
!      call ice_coupling_setup('ice',MPI_COMM_ICE)
!#endif
!#else
!      !  if running in stand-alone MPI mode
!#ifdef _MPI
!#if fcd_coupled
!#else
!      call MPI_INIT(ierr)
!#endif
!      call MPI_COMM_DUP(MPI_FVCOM_GROUP, MPI_COMM_ICE, ierr)
!#endif
!#endif
!
!#ifdef _MPI
!      call MPI_COMM_SIZE (MPI_COMM_ICE, nb_tasks, ierr)
!      call MPI_COMM_RANK (MPI_COMM_ICE, my_task, ierr)

!      if (nb_tasks /= nproc_s) then
!         write (6,*) 'nb_tasks, nproc_s =', nb_tasks, nproc_s
!         call abort_ice ('nb_tasks must equal nproc_s')
!      endif

!      if (real(imt_global,kind=dbl_kind)/real(nproc_x,kind=dbl_kind) 
!     &       /= int(imt_global/nproc_x)) then
!        write (6,*) 'nproc_x, imt_global =', nproc_x, imt_global
!        call abort_ice
!     &       ('number of pes in x must evenly divide imt_global')
!      endif

!      if (real(jmt_global,kind=dbl_kind)/real(nproc_y,kind=dbl_kind) 
!     &       /= int(jmt_global/nproc_y)) then
!        write (6,*) 'nproc_y, jmt_global =', nproc_y, jmt_global
!        call abort_ice
!     &       ('number of pes in y must evenly divide jmt_global')
!      endif

!      if ( ierr /= MPI_SUCCESS ) then 
!         call abort_ice('(setup_mpi) ERROR after MPI_COMM_xxx')
!      endif

!      coords1 = mod(my_task,nproc_x)
!      coords2 = my_task/nproc_x
!      nbr_east = coords2*nproc_x+mod(my_task+1,nproc_x)
!      nbr_west = coords2*nproc_x+mod(my_task-1+nproc_x,nproc_x)
!      nbr_north = my_task+nproc_x
!      nbr_south = my_task-nproc_x
!      if (nbr_south < 0) nbr_south = -1
!      if (nbr_north > nproc_s-1) nbr_north=-1
!
!      ilen = ihi-ilo+1
!      jlen = jhi-jlo+1
!
!      do n=1,nproc_s
!
!      local_start(1,n)=((imt_global-1)/nproc_x+1)*mod((n-1),nproc_x)+1
!      local_start(2,n)=((jmt_global-1)/nproc_y+1)*((n-1)/nproc_x)+1
!
!      call MPI_TYPE_VECTOR(jlen, ilen, ilen, 
!     &     mpi_integer, mpi_interior_int(n), ierr)
!      call MPI_TYPE_COMMIT(mpi_interior_int(n), ierr)
!
!      call MPI_TYPE_VECTOR(jlen, ilen, ilen, 
!     &     mpi_real8, mpi_interior_real(n), ierr)
!      call MPI_TYPE_COMMIT(mpi_interior_real(n), ierr)
!
!      call MPI_TYPE_VECTOR(jlen, ilen, imt_global, 
!     &     mpi_integer, mpi_interior_int_global(n), ierr)
!      call MPI_TYPE_COMMIT(mpi_interior_int_global(n), ierr)
!
!      call MPI_TYPE_VECTOR(jlen, ilen, imt_global, 
!     &     mpi_real8, mpi_interior_real_global(n), ierr)
!      call MPI_TYPE_COMMIT(mpi_interior_real_global(n), ierr)
!
!      enddo

!      do n=1,nproc_s
!         if (my_task == n-1) then
!            write (6,*) ' my_task,e,w,n,s ', my_task,
!     &           nbr_east, nbr_west, nbr_north, nbr_south
!         endif
!      enddo
!      write (6,*) ' '
!
!#else
!      ! not MPI
!      local_start(1,1)= 1
!      local_start(2,1)= 1
!      my_task = master_task
!      nb_tasks = 1 
!#endif

!      end subroutine setup_mpi

!=======================================================================

      end module ice_init

!=======================================================================
