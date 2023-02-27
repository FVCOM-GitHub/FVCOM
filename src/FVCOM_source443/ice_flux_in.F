!/===========================================================================/
! CVS VERSION INFORMATION
! $Id$
! $Name$
! $Revision$
!/===========================================================================/

!c=======================================================================
!
!BOP
!
! !MODULE: ice_flux_in - reads and interpolates input forcing data
!
! !DESCRIPTION:
!
! Reads and interpolates forcing data for atmospheric and oceanic quantities.
!
! !REVISION HISTORY:
!
! authors Elizabeth C. Hunke, LANL
!         William H. Lipscomb, LANL
!
! !INTERFACE:
!
      module ice_flux_in
!
! !USES:
!
      use ice_kinds_mod
      use ice_domain
      use ice_constants
      use ice_flux
      use ice_calendar
!      use ice_read_write
      use ice_fileunits
      USE CONTROL, ONLY : ICE_LONGWAVE_TYPE
!
!EOP
!
      implicit none
      save

      integer (kind=int_kind) ::  &
         ycycle             & ! number of years in forcing cycle
      ,  fyear_init         & ! first year of data in forcing cycle
      ,  fyear_final        & ! last year in cycle
      ,  fyear                ! current year in forcing cycle

!      real (kind=dbl_kind), dimension(ilo:ihi,jlo:jhi) :: &
      real (kind=dbl_kind), dimension(:,:),allocatable,save:: &
          cldf                ! cloud fraction

      character (char_len_long) :: &        ! input data file names
         height_file      &
      ,   uwind_file      &
      ,   vwind_file      &
      ,    potT_file      &
      ,    tair_file      &
      ,   humid_file      &
      ,    rhoa_file      &
      ,     fsw_file      &
      ,     flw_file      &
      ,    rain_file      &
      ,     sst_file      &
      ,     sss_file       

      real (kind=dbl_kind) :: & 
           c1intp, c2intp    & ! interpolation coefficients  
      ,    ftime              ! forcing time (for restart)

      integer (kind=int_kind) :: & 
           oldrecnum = 0      ! old record number (save between steps)

!      real (kind=dbl_kind), dimension(ilo:ihi,jlo:jhi,2)  :: &
      real (kind=dbl_kind), dimension(:,:,:),allocatable,save  :: &
            fsw_data     & ! field values at 2 temporal data points
      ,    cldf_data     &
      ,   fsnow_data     &
      ,    Tair_data     &
      ,    uatm_data     &
      ,    vatm_data     &
      ,      Qa_data     &
      ,    rhoa_data     &
      ,    potT_data     &
      ,    zlvl_data     &
      ,     flw_data     &
      ,     sst_data     &
      ,     sss_data      

      character (char_len) :: &
         atm_data_type   & ! 'default', 'ncar' or 'LYq'
      ,  sss_data_type   & ! 'default', 'clim' or 'ncar'
      ,  sst_data_type   & ! 'default', 'clim' or 'ncar'
      ,  precip_units     ! 'mm_per_month', 'mm_per_sec', 'mks'

      character(char_len_long) :: & 
         atm_data_dir     &! top directory for atmospheric data
      ,  ocn_data_dir     &! top directory for ocean data
      ,  oceanmixed_file   ! netCDF file name for ocean forcing data

      integer (kind=int_kind), parameter :: & 
         nfld = 8    ! number of fields to search for in netCDF file

!      real (kind=dbl_kind) :: & 
!         ocn_frc_m(ilo:ihi,jlo:jhi,nfld,12) ! ocn data for 12 months

      logical (kind=log_kind) :: &
         restore_sst                 ! restore sst if true

      integer (kind=int_kind) :: &
         trestore                    ! restoring time scale (days)

      real (kind=dbl_kind) :: & 
         trest                       ! restoring time scale (sec)

      logical (kind=log_kind) :: &
         dbug             ! prints debugging output if true

!c=======================================================================

      contains

!c=======================================================================
!
!BOP
!
! !IROUTINE: init_getflux - initialize input forcing files
!
! !INTERFACE:
!
      subroutine init_getflux
!
! !DESCRIPTION:
!
! Determines the current and final years of the forcing cycle based on
! namelist input, initializes the forcing data filenames, and initializes
! surface temperature and salinity from data.
!
! !REVISION HISTORY:
!
! authors Elizabeth C. Hunke, LANL
!
! !USES:
!
!      use ice_history, only: restart
!      use ice_therm_vertical, only: ustar_scale
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
!      fyear_final = fyear_init + ycycle - 1  ! last year in forcing cycle
!      fyear = fyear_init + mod(nyr-1,ycycle) ! current year

!      if (trim(atm_data_type) /= 'default' .and.         &
!                          my_task == master_task) then
!      write (nu_diag,*) ' Initial forcing data year = ',fyear_init
!      write (nu_diag,*) ' Final   forcing data year = ',fyear_final
!      endif

!      if (restore_sst) then
!        if (trestore == 0) then
!!          trest = dt                ! use data instantaneously
!          trest = dtice                ! use data instantaneously
!        else
!          trest = real(trestore,kind=dbl_kind) * secday ! seconds
!        endif
!      endif

      ! default forcing values from init_flux_atm
!      if (trim(atm_data_type) == 'ncar') then
!         call NCAR_files(fyear)     ! data for individual years
!      elseif (trim(atm_data_type) == 'LYq') then
!         call LY_files(fyear)       ! data for individual years
!      endif

      ! default forcing values from init_flux_ocn
!      if (trim(sss_data_type) == 'clim') then
!         call sss_clim              ! climatology (12-month avg)
!      endif
!      if (trim(sst_data_type) == 'clim' .and. .not. (restart)) then
!         call sst_ic                ! not interpolated but depends on sss
!      endif
!      if (trim(sst_data_type) == 'ncar' .or.  &
!          trim(sss_data_type) == 'ncar') then
!         call getflux_ocn_ncar_init
!      endif

! default ustar_scale = c1i (for nonzero currents) set in init_thermo_vertical
!      if (trim(sst_data_type) /= 'ncar' .or.     &
!         trim(sss_data_type) /= 'ncar') then
!         ustar_scale = c10i          ! for zero currents
!      endif

      end subroutine init_getflux

!c=======================================================================
!
!BOP
!
! !IROUTINE: getflux - Get forcing data and interpolate as necessary
!
! !INTERFACE:
!
      subroutine getflux
!
! !DESCRIPTION:
!
! Get forcing data and interpolate as necessary! 
!
! !REVISION HISTORY:
!
! authors: Elizabeth C. Hunke, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP

!
!      fyear = fyear_init + mod(nyr-1,ycycle)  ! current year
!      if (trim(atm_data_type) /= 'default' .and. istep <= 1 &
!                   .and. my_task == master_task) then
!      write (nu_diag,*) ' Current forcing data year = ',fyear
!      endif

!      ftime = time         ! forcing time
!      time_forc = ftime    ! for restarting

    !-------------------------------------------------------------------
    ! Read and interpolate annual climatologies of SSS and SST.
    ! Restore model SST to data.
    ! Interpolate ocean fields to U grid.
    !-------------------------------------------------------------------

!      if (trim(sst_data_type) == 'clim' .or.    &
!          trim(sss_data_type) == 'clim') then    
!         call getflux_ocn_clim                   
!      elseif (trim(sst_data_type) == 'ncar' .or.& 
!              trim(sss_data_type) == 'ncar') then
!         call getflux_ocn_ncar      
!      endif

    !-------------------------------------------------------------------
    ! Read and interpolate atmospheric data
    !-------------------------------------------------------------------

!      if (trim(atm_data_type) == 'ncar') then
!         call NCAR_bulk_data
!      elseif (trim(atm_data_type) == 'LYq') then
!         call LY_bulk_data
!     else ! default values set in init_flux*
!      endif

!      call prepare_forcing

      end subroutine getflux

!c=======================================================================
!
!BOP
!
! !IROUTINE: read_data - Read data needed for interpolation
!
! !INTERFACE:
!
!      subroutine read_data (flag, recd, yr, imx, ixx, ipx, &
!                             maxrec, data_file, field_data)
!
! !DESCRIPTION:
!
! If data is at the beginning of a one-year record, get data from
!  the previous year.
! If data is at the end of a one-year record, get data from the 
!  following year.
! If no earlier data exists (beginning of fyear\_init), then \\
!  (1) For monthly data, get data from the end of fyear\_final. \\
!  (2) For more frequent data, let the imx value equal the 
!      first value of the year. \\
! If no later data exists (end of fyear\_final), then \\
!  (1) For monthly data, get data from the beginning of fyear\_init. \\
!  (2) For more frequent data, let the ipx 
!      value equal the last value of the year. \\
! In other words, we assume persistence when daily or 6-hourly
!   data is missing, and we assume periodicity when monthly data
!   is missing. \\
!
! !REVISION HISTORY:
!
! authors: same as module
!
! !USES:
!
!      use ice_diagnostics
!
! !INPUT/OUTPUT PARAMETERS:
!
!      logical (kind=log_kind), intent(in) :: flag
!
!      integer (kind=int_kind), intent(in) :: & 
!        recd                    &! baseline record number
!      , yr                      &! year of forcing data
!      , imx, ixx, ipx           &! record numbers of 3 data values
!                                 ! relative to recd
!      , maxrec                   ! maximum record value
!
!      real (kind=dbl_kind), dimension(ilo:ihi,jlo:jhi,2), intent(out) :: &
!        field_data              ! 2 values needed for interpolation
!!
!!EOP
!!
!      character (char_len_long) :: &
!        data_file               ! data file to be read 
!
!      integer (kind=int_kind) :: & 
!        nbits              &  ! = 32 for single precision, 64 for double
!      , nrec               &  ! record number to read
!      , n2, n4             &  ! like imx and ipx, but 
!                              ! adjusted at beginning and end of data
!      , arg                 ! value of time argument in field_data

!      nbits = 64                ! double precision data

!      if (istep1 > check_step) dbug = .true.   !! debugging
!      if (my_task==master_task .and. (dbug))  &
!         write(nu_diag,*) '  ',data_file

!      if (flag) then
      !-----------------------------------------------------------------
      ! Initialize record counters
      ! (n2, n4 will change only at the very beginning or end of
      !  a forcing cycle.)
      !-----------------------------------------------------------------
!      n2 = imx
!      n4 = ipx
!      arg = 0
      
      !-----------------------------------------------------------------
      ! read data
      !-----------------------------------------------------------------

!      if (imx  /=  99) then
!      ! currently in first half of data interval
!        if (ixx<=1) then
!          if (yr > fyear_init) then  ! get data from previous year
!            call file_year (data_file, yr-1)
!          else                   ! yr = fyear_init, no prior data exists
!            if (maxrec > 12) then  ! extrapolate from first record
!              if (ixx==1) n2 = ixx
!            else                 ! go to end of fyear_final
!              call file_year (data_file, fyear_final)
!            endif
!          endif                  ! yr > fyear_init
!        endif                    ! ixx <= 1
!        call ice_open (nu_forcing, data_file, nbits)

!        arg = 1
!        nrec = recd + n2
!        call ice_read (nu_forcing, nrec, field_data(:,:,arg),       &
!                      'rda8', dbug)
! 
!        if (ixx==1 .and. my_task==master_task) close(nu_forcing)
!      endif  ! imx ne 99

!      ! always read ixx data from data file for current year
!      call file_year (data_file, yr)
!      call ice_open (nu_forcing, data_file, nbits)

!      arg = arg + 1
!      nrec = recd + ixx
!      call ice_read (nu_forcing, nrec, field_data(:,:,arg), &
!                      'rda8', dbug)

!      if (ipx  /=  99) then 
!      ! currently in latter half of data interval
!        if (ixx==maxrec) then
!          if (yr < fyear_final) then ! get data from following year
!            if (my_task == master_task) close(nu_forcing)
!            call file_year (data_file, yr+1)
!            call ice_open (nu_forcing, data_file, nbits)
!          else                     ! yr = fyear_final, no more data exists
!            if (maxrec > 12) then  ! extrapolate from ixx
!              n4 = ixx
!            else                   ! go to beginning of fyear_init
!              if (my_task == master_task) close(nu_forcing)
!              call file_year (data_file, fyear_init)
!              call ice_open (nu_forcing, data_file, nbits)
!            endif
!          endif                    ! yr < fyear_final
!        endif                      ! ixx = maxrec

!        arg = arg + 1
!        nrec = recd + n4
!        call ice_read (nu_forcing, nrec, field_data(:,:,arg),&
!                     'rda8', dbug)
!      endif  ! ipx ne 99

!      if (my_task == master_task) close(nu_forcing)
!      endif  ! flag

!      end subroutine read_data

!c=======================================================================
!
!BOP
!
! !IROUTINE: read_clim_data - read annual climatological data
!
! !INTERFACE:
!
!      subroutine read_clim_data (readflag, recd, imx, ixx, ipx,&
!                                data_file, field_data)
!
! !DESCRIPTION:
!
! Read data needed for interpolation, as in read\_data.
! Assume a one-year cycle of climatological data, so that there is
!  no need to get data from other years or to extrapolate data beyond
!  the forcing time period.
!
! !REVISION HISTORY:
!
! authors: same as module
!
! !USES:
!
!      use ice_diagnostics  !! debugging
!
! !INPUT/OUTPUT PARAMETERS:
!!
!      logical (kind=log_kind),intent(in) :: readflag
!
!      integer (kind=int_kind), intent(in) :: & 
!       recd              &  ! baseline record number
!      , imx,ixx,ipx         ! record numbers of 3 data values
!                            ! relative to recd
!
!      character (char_len_long), intent(in) ::  data_file
!
!      real (kind=dbl_kind), dimension(ilo:ihi,jlo:jhi,2), intent(out) :: &
!        field_data         ! 2 values needed for interpolation
!!
!!EOP
!!
!      integer (kind=int_kind) :: & 
!        nbits             & ! = 32 for single precision, 64 for double
!      , nrec              & ! record number to read
!      , arg                ! value of time argument in field_data
!
!      nbits = 64                ! double precision data

!      if (istep1 > check_step) dbug = .true. !! debugging

!      if (my_task==master_task .and. (dbug))       &
!        write(nu_diag,*) '  ', trim(data_file)

!      if (readflag) then
      !-----------------------------------------------------------------
      ! read data
      !-----------------------------------------------------------------
      
!      call ice_open (nu_forcing, data_file, nbits)

!      arg = 0
!      if (imx  /=  99) then
!        arg = 1
!        nrec = recd + imx
!        call ice_read (nu_forcing, nrec, field_data(:,:,arg), &
!                       'rda8', dbug)                           
!      endif                                                    
                                                               
!      arg = arg + 1                                            
!      nrec = recd + ixx                                        
!      call ice_read (nu_forcing, nrec, field_data(:,:,arg),   &
!                       'rda8', dbug)                           
                                                               
!      if (ipx  /=  99) then                                    
!        arg = arg + 1                                          
!        nrec = recd + ipx                                      
!        call ice_read (nu_forcing, nrec, field_data(:,:,arg), &
!                       'rda8', dbug)
!      endif

!      if (my_task == master_task) close (nu_forcing)
!      endif  ! readflag

!      end subroutine read_clim_data

!c=======================================================================
!
!BOP
!
! !IROUTINE: interp_coeff_monthly - Compute monthly data interpolation coefficients
!
! !INTERFACE:
!
      subroutine interp_coeff_monthly (recslot)
!
! !DESCRIPTION:
!
! Compute coefficients for interpolating monthly data to current time step.
!
! !REVISION HISTORY:
!
! authors: same as module
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
           recslot         ! slot (1 or 2) for current record
!
!EOP
!
      real (kind=dbl_kind) :: &
          tt             &   ! seconds elapsed in current year
      ,   t1, t2           ! seconds elapsed at month midpoint

      real (kind=dbl_kind) :: & 
          daymid(0:13)     ! month mid-points

      daymid(1:13) = 14._dbl_kind   ! time frame ends 0 sec into day 15
      daymid(0)    = -17._dbl_kind  ! Dec 15, 0 sec

      ! make time cyclic
      tt = mod(ftime/secday,c365)
 
      ! Find neighboring times
      
      if (recslot==2) then      ! first half of month
        t2 = daycal(month) + daymid(month)   ! midpoint, current month
        if (month == 1) then
          t1 = daymid(0)                 ! Dec 15 (0 sec)
        else
          t1 = daycal(month-1) + daymid(month-1) ! midpoint, previous month
        endif
      else                      ! second half of month
        t1 = daycal(month) + daymid(month)    ! midpoint, current month
        t2 = daycal(month+1) + daymid(month+1)! day 15 of next month (0 sec)
      endif
 
      ! Compute coefficients
      c1intp = (t2 - tt) / (t2 - t1)
      c2intp =  c1i - c1intp

      end subroutine interp_coeff_monthly

!c=======================================================================
!
!BOP
!
! !IROUTINE: interp_coeff
!
! !INTERFACE:
!
      subroutine interp_coeff (recnum, recslot, secint)
!
! !DESCRIPTION:
!
! Compute coefficients for interpolating data to current time step.
! Works for any data interval that divides evenly into a 365-day 
!  year (daily, 6-hourly, etc.)
! Use interp\_coef\_monthly for monthly data.
!
! !REVISION HISTORY:
!
! authors: same as module
!
! !USES:
!
      integer (kind=int_kind), intent(in) :: & 
          recnum         & ! record number for current data value
      ,   recslot         ! spline slot for current record 

      real (kind=dbl_kind), intent(in) :: &
           secint                    ! seconds in data interval
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      real (kind=dbl_kind), parameter :: &
         secyr = c365 * secday     ! seconds in a 365-day year 

      real (kind=dbl_kind) :: &
          tt           &    ! seconds elapsed in current year
      ,   t1, t2           ! seconds elapsed at data points

!      tt = mod(ftime,secyr)
 
      ! Find neighboring times
!      if (recslot==2) then          ! current record goes in slot 2 (NCEP)
!         t2 = real(recnum,kind=dbl_kind)*secint
!         t1 = t2 - secint           !  - 1 interval
!      else                          ! recslot = 1
!         t1 = real(recnum,kind=dbl_kind)*secint
!         t2 = t1 + secint           !  + 1 interval
!      endif
 
      ! Compute coefficients
!      c1intp =  abs((t2 - tt) / (t2 - t1))
!      c2intp =  c1i - c1intp

      end subroutine interp_coeff

!c=======================================================================
!
!BOP
!
! !IROUTINE: file_year - 
!
! !INTERFACE:
!
      subroutine file_year (data_file, yr)
!
! !DESCRIPTION:
!
! Construct the correct name of the atmospheric data file
! to be read, given the year and assuming the naming convention
! that filenames end with 'yyyy.dat'.
!
! !REVISION HISTORY:
!
! authors Elizabeth C. Hunke, LANL
!
! !USES:
!
!
! !INPUT/OUTPUT PARAMETERS:
!
      character (char_len_long), intent(inout) ::  data_file
!
!EOP
!
      integer (kind=int_kind), intent(in) :: yr

      character (char_len_long) :: tmpname

      integer (kind=int_kind) :: i

      i = index(data_file,'.dat') - 5
      tmpname = data_file
      write(data_file,'(a,i4.4,a)') tmpname(1:i), yr, '.dat'

      end subroutine file_year

!c=======================================================================
!
!BOP
!
! !IROUTINE: NCAR_files - construct filenames for NCAR_bulk atmospheric data
!
! !INTERFACE:
!
      subroutine NCAR_files (yr)
!
! !DESCRIPTION:
!
! This subroutine is based on the LANL file naming conventions.
! Edit for other directory structures or filenames.
! Note: The year number in these filenames does not matter, because
!       subroutine file\_year will insert the correct year.
! 
!
! !REVISION HISTORY:
!
! authors Elizabeth C. Hunke, LANL
!
! !USES:
!
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) ::  & 
           yr                   ! current forcing year
!
!EOP
!
      fsw_file =                                                      &    
           trim(atm_data_dir)//'ISCCPM/MONTHLY/RADFLX/swdn.1996.dat'   
      call file_year(fsw_file,yr)                                      
                                                                       
      flw_file =                                                      &
           trim(atm_data_dir)//'ISCCPM/MONTHLY/RADFLX/cldf.1996.dat'   
      call file_year(flw_file,yr)                                      
                                                                       
      rain_file =                                                     &
           trim(atm_data_dir)//'MXA/MONTHLY/PRECIP/prec.1996.dat'      
      call file_year(rain_file,yr)                                     
                                                                       
      uwind_file =                                                    &
           trim(atm_data_dir)//'NCEP/4XDAILY/STATES/u_10.1996.dat'     
      call file_year(uwind_file,yr)                                    
                                                                       
      vwind_file =                                                    &
           trim(atm_data_dir)//'NCEP/4XDAILY/STATES/v_10.1996.dat'     
      call file_year(vwind_file,yr)                                    
                                                                       
      tair_file =                                                     &
           trim(atm_data_dir)//'NCEP/4XDAILY/STATES/t_10.1996.dat'     
      call file_year(tair_file,yr)                                     
                                                                       
      humid_file =                                                    &
           trim(atm_data_dir)//'NCEP/4XDAILY/STATES/q_10.1996.dat'     
      call file_year(humid_file,yr)                                    
                                                                       
      rhoa_file =                                                     &
           trim(atm_data_dir)//'NCEP/4XDAILY/STATES/dn10.1996.dat'
      call file_year(rhoa_file,yr)

      if (my_task == master_task) then
         write (nu_diag,*) ''
         write (nu_diag,*) 'Initial atmospheric data files:'
         write (nu_diag,*) trim(fsw_file)
         write (nu_diag,*) trim(flw_file)
         write (nu_diag,*) trim(rain_file)
         write (nu_diag,*) trim(uwind_file)
         write (nu_diag,*) trim(vwind_file)
         write (nu_diag,*) trim(tair_file)
         write (nu_diag,*) trim(humid_file)
         write (nu_diag,*) trim(rhoa_file)
      endif                     ! master_task

      end subroutine NCAR_files

!c=======================================================================
!
!BOP
!
! !IROUTINE: LY_files - construct filenames for LY atmospheric data
!
! !INTERFACE:
!
      subroutine LY_files (yr)
!
! !DESCRIPTION:
!
! This subroutine is based on the LANL file naming conventions.
! Edit for other directory structures or filenames.
! Note: The year number in these filenames does not matter, because
!       subroutine file\_year will insert the correct year.
! 
!
! !REVISION HISTORY:
!
! authors Elizabeth C. Hunke, LANL
!
! !USES:
!
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) ::             &
           yr                   ! current forcing year    
!                                                         
!EOP                                                      
!                                                         
      flw_file =                                         &
           trim(atm_data_dir)//'MONTHLY/cldf.omip.dat'    
                                                          
      rain_file =                                        &
           trim(atm_data_dir)//'MONTHLY/prec.nmyr.dat'    
                                                          
      uwind_file =                                       &
           trim(atm_data_dir)//'4XDAILY/u_10.1996.dat'    
      call file_year(uwind_file,yr)                       
                                                          
      vwind_file =                                       &
           trim(atm_data_dir)//'4XDAILY/v_10.1996.dat'    
      call file_year(vwind_file,yr)                       
                                                          
      tair_file =                                        &
           trim(atm_data_dir)//'4XDAILY/t_10.1996.dat'    
      call file_year(tair_file,yr)                        
                                                          
      humid_file =                                       &
           trim(atm_data_dir)//'4XDAILY/q_10.1996.dat'
      call file_year(humid_file,yr)

      if (my_task == master_task) then
         write (nu_diag,*) ''
         write (nu_diag,*) 'Forcing data year = ', fyear
         write (nu_diag,*) 'Atmospheric data files:'
         write (nu_diag,*) trim(flw_file)
         write (nu_diag,*) trim(rain_file)
         write (nu_diag,*) trim(uwind_file)
         write (nu_diag,*) trim(vwind_file)
         write (nu_diag,*) trim(tair_file)
         write (nu_diag,*) trim(humid_file)
      endif                     ! master_task

      end subroutine LY_files

!c=======================================================================
!
!BOP
!
! !IROUTINE: NCAR_bulk_data - read NCAR_bulk atmospheric data
!
! !INTERFACE:
!
      subroutine NCAR_bulk_data
!
! !DESCRIPTION:
!
! !REVISION HISTORY:
!
! authors: same as module
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) ::& 
          i, j                  &
      ,   imx,ixx,ipx           &! record numbers for neighboring months
      ,   recnum                &! record number
      ,   maxrec                &! maximum record number
      ,   recslot               &! spline slot for current record
      ,   midmonth               ! middle day of month

      real (kind=dbl_kind) ::  &
          sec6hr              ! number of seconds in 6 hours

      logical (kind=log_kind) :: readm, read6

    !-------------------------------------------------------------------
    ! monthly data 
    !
    ! Assume that monthly data values are located in the middle of the 
    ! month.
    !-------------------------------------------------------------------
      
      midmonth = 15  ! data is given on 15th of every month
!      midmonth = fix(p5 * real(daymo(month),kind=dbl_kind))  ! exact middle

      ! Compute record numbers for surrounding months
      maxrec = 12
      imx  = mod(month+maxrec-2,maxrec) + 1
      ipx  = mod(month,         maxrec) + 1
      if (mday >= midmonth) imx = 99  ! other two points will be used
      if (mday <  midmonth) ipx = 99

      ! Determine whether interpolation will use values 1:2 or 2:3
      ! recslot = 2 means we use values 1:2, with the current value (2)
      !  in the second slot
      ! recslot = 1 means we use values 2:3, with the current value (2)
      !  in the first slot
      recslot = 1                             ! latter half of month
      if (mday < midmonth) recslot = 2        ! first half of month

      ! Find interpolation coefficients
      call interp_coeff_monthly (recslot)

      ! Read 2 monthly values 
      readm = .false.
      if (istep==1 .or. (mday==midmonth .and. sec==0)) readm = .true.

!      call read_data (readm, 0, fyear, imx, month, ipx, &
!            maxrec, fsw_file, fsw_data)                  
!      call read_data (readm, 0, fyear, imx, month, ipx, &
!            maxrec, flw_file, cldf_data)                 
!      call read_data (readm, 0, fyear, imx, month, ipx, &
!            maxrec, rain_file, fsnow_data)

!      call interpolate_data (fsw_data, fsw)
!      call interpolate_data (cldf_data, cldf)
!      call interpolate_data (fsnow_data, fsnow)

    !-------------------------------------------------------------------
    ! 6-hourly data
    ! 
    ! Assume that the 6-hourly value is located at the end of the
    !  6-hour period.  This is the convention for NCEP reanalysis data.
    !  E.g. record 1 gives conditions at 6 am GMT on 1 January.
    !-------------------------------------------------------------------

      sec6hr = secday/c4i        ! seconds in 6 hours
      maxrec = 1460             ! 365*4

      ! current record number
      recnum = 4*int(yday) - 3 + int(real(sec,kind=dbl_kind)/sec6hr)

      ! Compute record numbers for surrounding data (2 on each side)

      imx = mod(recnum+maxrec-2,maxrec) + 1
      ixx = mod(recnum-1,       maxrec) + 1
!      ipx = mod(recnum,         maxrec) + 1

      ! Compute interpolation coefficients
      ! If data is located at the end of the time interval, then the
      !  data value for the current record goes in slot 2

      recslot = 2
      ipx = 99
      call interp_coeff (recnum, recslot, sec6hr)

      ! Read
      read6 = .false.
      if (istep==1 .or. oldrecnum  /=  recnum) read6 = .true.

!      call read_data (read6, 0, fyear, imx, ixx, ipx, maxrec,   &
!                      tair_file, Tair_data)                      
!      call read_data (read6, 0, fyear, imx, ixx, ipx, maxrec,   &
!                      uwind_file, uatm_data)                     
!      call read_data (read6, 0, fyear, imx, ixx, ipx, maxrec,   &
!                      vwind_file, vatm_data)                     
!      call read_data (read6, 0, fyear, imx, ixx, ipx, maxrec,   &
!                      rhoa_file, rhoa_data)                      
!      call read_data (read6, 0, fyear, imx, ixx, ipx, maxrec,   &
!                      humid_file, Qa_data)

      ! Interpolate
!      call interpolate_data (Tair_data, Tair)
!      call interpolate_data (uatm_data, uatm)
!      call interpolate_data (vatm_data, vatm)
!      call interpolate_data (rhoa_data, rhoa)
!      call interpolate_data (Qa_data, Qa)

      ! Save record number
      oldrecnum = recnum

      end subroutine NCAR_bulk_data

!c=======================================================================
!
!BOP
!
! !IROUTINE: LY_bulk_data - read LY atmospheric data
!
! !INTERFACE:
!
      subroutine LY_bulk_data
!
! !DESCRIPTION:
!
! !REVISION HISTORY:
!
! authors: Elizabeth C. Hunke, LANL
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
          i, j              &    
      ,   imx,ixx,ipx       &  ! record numbers for neighboring months
      ,   recnum            &  ! record number
      ,   maxrec            &  ! maximum record number
      ,   recslot           &  ! spline slot for current record
      ,   midmonth          &  ! middle day of month
      ,   ng                 

      real (kind=dbl_kind) ::  &
          sec6hr              ! number of seconds in 6 hours

      logical (kind=log_kind) :: readm, read6

    !-------------------------------------------------------------------
    ! monthly data 
    !
    ! Assume that monthly data values are located in the middle of the 
    ! month.
    !-------------------------------------------------------------------

      midmonth = 15  ! data is given on 15th of every month
!      midmonth = fix(p5 * real(daymo(month),kind=dbl_kind))  ! exact middle

      ! Compute record numbers for surrounding months
      maxrec = 12
      imx  = mod(month+maxrec-2,maxrec) + 1
      ipx  = mod(month,         maxrec) + 1
      if (mday >= midmonth) imx = 99  ! other two points will be used
      if (mday <  midmonth) ipx = 99

      ! Determine whether interpolation will use values 1:2 or 2:3
      ! recslot = 2 means we use values 1:2, with the current value (2)
      !  in the second slot
      ! recslot = 1 means we use values 2:3, with the current value (2)
      !  in the first slot
      recslot = 1                             ! latter half of month
      if (mday < midmonth) recslot = 2        ! first half of month

      ! Find interpolation coefficients
      call interp_coeff_monthly (recslot)

      ! Read 2 monthly values 
      readm = .false.
      if (istep==1 .or. (mday==midmonth .and. sec==0)) readm = .true.

!      call read_clim_data (readm, 0, imx, month, ipx, &
!             flw_file, cldf_data)                      
!      call read_clim_data (readm, 0, imx, month, ipx, &
!             rain_file, fsnow_data)
!
!      call interpolate_data (cldf_data, cldf)
!      call interpolate_data (fsnow_data, fsnow)  ! units mm/s = kg/m^2/s

    !-------------------------------------------------------------------
    ! 6-hourly data
    ! 
    ! Assume that the 6-hourly value is located at the end of the
    !  6-hour period.  This is the convention for NCEP reanalysis data.
    !  E.g. record 1 gives conditions at 6 am GMT on 1 January.
    !-------------------------------------------------------------------

      sec6hr = secday/c4i        ! seconds in 6 hours
      maxrec = 1460             ! 365*4

      ! current record number
      recnum = 4*int(yday) - 3 + int(real(sec,kind=dbl_kind)/sec6hr)

      ! Compute record numbers for surrounding data (2 on each side)

      imx = mod(recnum+maxrec-2,maxrec) + 1
      ixx = mod(recnum-1,       maxrec) + 1
!      ipx = mod(recnum,         maxrec) + 1

      ! Compute interpolation coefficients
      ! If data is located at the end of the time interval, then the
      !  data value for the current record goes in slot 2

      recslot = 2
      ipx = 99
      call interp_coeff (recnum, recslot, sec6hr)

      ! Read
      read6 = .false.
      if (istep==1 .or. oldrecnum  /=  recnum) read6 = .true.

!      call read_data (read6, 0, fyear, imx, ixx, ipx, maxrec, &
!                      tair_file, Tair_data)                    
!      call read_data (read6, 0, fyear, imx, ixx, ipx, maxrec, &
!                      uwind_file, uatm_data)                   
!      call read_data (read6, 0, fyear, imx, ixx, ipx, maxrec, &
!                      vwind_file, vatm_data)                   
!      call read_data (read6, 0, fyear, imx, ixx, ipx, maxrec, &
!                     humid_file, Qa_data)

      ! Interpolate
!      call interpolate_data (Tair_data, Tair)
!      call interpolate_data (uatm_data, uatm)
!      call interpolate_data (vatm_data, vatm)
!      call interpolate_data (Qa_data, Qa)

!      call Qa_fixLY

      do j=jlo,jhi
       do i=ilo,ihi
        Qa(i,j)   = Qa(i,j)   * hm(i,j)
        Tair(i,j) = Tair(i,j) * hm(i,j)
        uatm(i,j) = uatm(i,j) * hm(i,j)
        vatm(i,j) = vatm(i,j) * hm(i,j)
       enddo
      enddo

      call compute_shortwave ! AOMIP

      ! Save record number
      oldrecnum = recnum

         if (dbug) then
           ! cice
           if (my_task == master_task) write (nu_diag,*) 'LY_bulk_dat'
           ng = (ihi-ilo+1)*(jhi-jlo+1)
!           call ice_global_real_minmax(ng,Fsw,'A Fsw   ')
!           call ice_global_real_minmax(ng,cldf,'A cld   ')
!           call ice_global_real_minmax(ng,Fsnow,'A Fsnow ')
!           call ice_global_real_minmax(ng,Tair,'A Tair  ')
!           call ice_global_real_minmax(ng,uatm,'A uatm  ')
!           call ice_global_real_minmax(ng,vatm,'A vatm  ')
!           call ice_global_real_minmax(ng,Qa,'A Qa    ')
!           call ice_global_real_minmax(ng,rhoa,'A rhoa  ')
         endif

      end subroutine LY_bulk_data

!c=======================================================================
!
!BOP
!
! !IROUTINE: compute_shortwave - AOMIP shortwave forcing
!
! !INTERFACE:
!
      subroutine compute_shortwave
!
! !DESCRIPTION:
!
! Computes downwelling shortwave radiation based on sun declination
! and cloud fraction, as in AOMIP, following Zillman (1972) and
! Parkinson and Washington (1979).
!
! !REVISION HISTORY:
!
! authors: Elizabeth C. Hunke, LANL
!
! !USES:
!
      use ice_grid
      use ice_albedo
      use ice_constants
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      real (kind=dbl_kind) :: &
          hour_angle          &
      ,   solar_time          &
      ,   declin              &
      ,   cosZ                &
      ,   e_i, d_i                &
      ,   sw0                 &
      ,   deg2rad_i

      integer (kind=int_kind) :: i,j

      do j=jlo,jhi
       do i=ilo,ihi
        deg2rad_i = pi/180._dbl_kind
        solar_time = mod(real(sec,kind=dbl_kind),secday)/3600._dbl_kind &
                   + c12*sin(p5*TLON(i,j))                               
        hour_angle = (c12 - solar_time)*pi/c12                           
        declin = 23.44_dbl_kind*cos((172._dbl_kind-yday)                &
                 * c2i*pi/c365)*deg2rad_i                                
        cosZ = sin(TLAT(i,j))*sin(declin)                               &
             + cos(TLAT(i,j))*cos(declin)*cos(hour_angle)                
        cosZ = max(cosZ,c0i)                                              
        e_i = 1.e5_dbl_kind*Qa(i,j)                                       &
            /(0.622_dbl_kind + 0.378_dbl_kind*Qa(i,j))
        d_i = (cosZ+2.7_dbl_kind)*e_i*1.e-5_dbl_kind+1.085_dbl_kind*cosZ+p1
        sw0 = 1353._dbl_kind*cosZ**2/d_i
        sw0 = max(sw0,c0i)

        ! total downward shortwave for cice
        Fsw(i,j) = sw0*(c1i-p6*cldf(i,j)**3) 
        Fsw(i,j) = Fsw(i,j)*hm(i,j)
       enddo
      enddo

      end subroutine compute_shortwave

!c=======================================================================
!
!BOP
!
! !IROUTINE: Qa_fixLY
!
! !INTERFACE:
!
      subroutine Qa_fixLY
!
! !DESCRIPTION:
!
! Limits relative humidity <= 100%
!
! !REVISION HISTORY:
!
! authors: Elizabeth C. Hunke, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      real (kind=dbl_kind),            &
        dimension(ilo:ihi,jlo:jhi) ::  &
          work, slope

      work = Tair - Tffresh
      work = c2i + (0.7859_dbl_kind + 0.03477_dbl_kind*work)                    &
                     /(c1i + 0.00412_dbl_kind*work) &! 2+ converts ea mb -> Pa
                + 0.00422_dbl_kind*work            ! for ice
      ! vapor pressure
      work = (c10i**work)      ! saturated 
      work = max(work,puny)   ! puny over land to prevent division by zero
      ! specific humidity
      work = 0.622_dbl_kind*work/(1.e5_dbl_kind - 0.378_dbl_kind*work)

      Qa = min(Qa, work)

      end subroutine Qa_fixLY

!c=======================================================================
!
!BOP
!
! !IROUTINE: prepare_forcing - finish manipulating forcing
!
! !INTERFACE:
!
      subroutine prepare_forcing
!
! !DESCRIPTION:
!
! !REVISION HISTORY:
!
! authors Elizabeth C. Hunke, LANL
!
! !USES:
!
      use ice_grid, only: ANGLET, t2ugrid, hm
      use ice_state
      use mod_utils, only: Fatal_Error
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) :: i, j
      real (kind=dbl_kind) ::   &
          workx, worky          &
      ,   fcc,sstk,rtea,qlwm,ptem ! terms needed for lwdn computation

      do j=jlo,jhi
       do i=ilo,ihi

      !-----------------------------------------------------------------
      ! Fix interpolation errors
      !-----------------------------------------------------------------
        fsw (i,j)  = max(fsw(i,j),c0i)
        cldf(i,j)  = max(min(cldf(i,j),c1i),c0i)
        fsnow(i,j) = max(fsnow(i,j),c0i)   
        rhoa(i,j)  = max(rhoa(i,j),c0i)
        Qa  (i,j)  = max(Qa(i,j),c0i)

       enddo
      enddo

      !-----------------------------------------------------------------
      ! calculations specific to data sets
      !-----------------------------------------------------------------    

      if (trim(atm_data_type) == 'ncar') then

         ! correct known biases in NCAR data (as in CCSM latm)
         do j=jlo,jhi
         do i=ilo,ihi
            Qa(i,j)  = Qa(i,j)  * 0.94_dbl_kind
            fsw(i,j) = fsw(i,j) * 0.92_dbl_kind
         enddo
         enddo

      endif                     ! atm_data_type

      !-----------------------------------------------------------------
      ! Compute other fields needed by model
      !-----------------------------------------------------------------

      if(trim(ICE_LONGWAVE_TYPE) == 'PW')then

      do j=jlo,jhi
       do i = ilo, ihi

        zlvl(i,j) = c10i
        !wind (i,j) = sqrt(uatm(i,j)**2 + vatm(i,j)**2) ! wind speed, m/s
        potT(i,j) = Tair(i,j)

        ! divide shortwave into spectral bands
        swvdr(i,j) = fsw(i,j)*(.28_dbl_kind)         ! visible direct
        swvdf(i,j) = fsw(i,j)*(.24_dbl_kind)         ! visible diffuse
        swidr(i,j) = fsw(i,j)*(.31_dbl_kind)         ! near IR direct
        swidf(i,j) = fsw(i,j)*(.17_dbl_kind)         ! near IR diffuse
                                            ! as in the dummy atm (latm)
        ! longwave as in Parkinson and Washington (1979)
        flw(i,j) = stefan_boltzmann*Tair(i,j)**4 &! downward longwave
     &          * (c1i - 0.261_dbl_kind*         &
     &                  exp(-7.77e-4_dbl_kind*(Tffresh - Tair(i,j))**2)) &
     &          * (c1i + 0.275_dbl_kind*cldf(i,j))


        ! longwave, Rosati and Miyakoda, JPO 18, p. 1607 (1988) - sort of -
!        fcc = c1i - 0.8_dbl_kind * cldf(i,j)
!        sstk = (Tsfc(i,j) * aice(i,j)                         &
!               + sst(i,j) * (c1i - aice(i,j))) + Tffresh        
!        rtea = sqrt(c1000*Qa(i,j) / (0.622_dbl_kind           &
!             + 0.378_dbl_kind*Qa(i,j)))                        
!        ptem = Tair(i,j)  ! get this from stability?           
!        qlwm = ptem * ptem * ptem *                           &
!             ( ptem*(0.39_dbl_kind-0.05_dbl_kind*rtea)*fcc    &
!                                  + c4i*(sstk-ptem) )
!        flw(i,j) = emissivity * stefan_boltzmann * ( sstk**4 - qlwm )

!        flw(i,j) = flw(i,j) * hm(i,j) ! land mask

        ! determine whether precip is rain or snow
!        if (trim(precip_units) == 'mm_per_month') then
!          fsnow(i,j) = fsnow(i,j)/2.592e+06_dbl_kind  ! mm/month -> kg/m^2 s
!       elseif (trim(precip_units) == 'mm_per_sec' .or.
!               trim(precip_units) == 'mks') then
!               ! no change:  mm/sec = kg/m^2 s
!        endif
!        frain(i,j) = c0i                     
!        if (Tair(i,j) >= Tffresh) then
!            frain(i,j) = fsnow(i,j)
!            fsnow(i,j) = c0i
!        endif

      !-----------------------------------------------------------------
      ! rotate zonal/meridional vectors to local coordinates
      ! Vector fields come in on T grid, but are oriented geographically
      ! need to rotate to pop-grid FIRST using ANGLET
      ! then interpolate to the U-cell centers  (otherwise we
      ! interpolate across the pole)
      ! use ANGLET which is on the T grid !
      ! atmo variables are needed in T cell centers in subroutine stability,
      ! and are interpolated to the U grid later as necessary
      !-----------------------------------------------------------------
        !workx      = uatm(i,j)                ! wind velocity, m/s
        !worky      = vatm(i,j) 
        !uatm (i,j) = workx*cos(ANGLET(i,j))  & ! convert to POP grid
        !           + worky*sin(ANGLET(i,j))    ! note uatm, vatm, wind
        !vatm (i,j) = worky*cos(ANGLET(i,j))  & !  are on the T-grid here
        !           - workx*sin(ANGLET(i,j))

       enddo
      enddo
      
      else if(trim(ICE_LONGWAVE_TYPE) == 'RM')then

      do j=jlo,jhi
       do i = ilo, ihi

        zlvl(i,j) = c10i
        !wind (i,j) = sqrt(uatm(i,j)**2 + vatm(i,j)**2) ! wind speed, m/s
        potT(i,j) = Tair(i,j)

        ! divide shortwave into spectral bands
        swvdr(i,j) = fsw(i,j)*(.28_dbl_kind)         ! visible direct
        swvdf(i,j) = fsw(i,j)*(.24_dbl_kind)         ! visible diffuse
        swidr(i,j) = fsw(i,j)*(.31_dbl_kind)         ! near IR direct
        swidf(i,j) = fsw(i,j)*(.17_dbl_kind)         ! near IR diffuse
                                            ! as in the dummy atm (latm)

        ! longwave, Rosati and Miyakoda, JPO 18, p. 1607 (1988) - sort of -
        fcc = c1i - 0.8_dbl_kind * cldf(i,j)
        sstk = (Tsfc(i,j) * aice(i,j)                         &
               + sst(i,j) * (c1i - aice(i,j))) + Tffresh
        rtea = sqrt(c1000*Qa(i,j) / (0.622_dbl_kind           &
             + 0.378_dbl_kind*Qa(i,j)))
        ptem = Tair(i,j)  ! get this from stability?           
        qlwm = ptem * ptem * ptem *                           &
             ( ptem*(0.39_dbl_kind-0.05_dbl_kind*rtea)*fcc    &
                                  + c4i*(sstk-ptem) )
        flw(i,j) = emissivity * stefan_boltzmann * ( sstk**4 - qlwm )

        flw(i,j) = flw(i,j) * hm(i,j) ! land mask

        ! determine whether precip is rain or snow
!        if (trim(precip_units) == 'mm_per_month') then
!          fsnow(i,j) = fsnow(i,j)/2.592e+06_dbl_kind  ! mm/month -> kg/m^2 s
!       elseif (trim(precip_units) == 'mm_per_sec' .or.
!               trim(precip_units) == 'mks') then
!               ! no change:  mm/sec = kg/m^2 s
!        endif
!        frain(i,j) = c0i                     
!        if (Tair(i,j) >= Tffresh) then
!            frain(i,j) = fsnow(i,j)
!            fsnow(i,j) = c0i
!        endif

      !-----------------------------------------------------------------
      ! rotate zonal/meridional vectors to local coordinates
      ! Vector fields come in on T grid, but are oriented geographically
      ! need to rotate to pop-grid FIRST using ANGLET
      ! then interpolate to the U-cell centers  (otherwise we
      ! interpolate across the pole)
      ! use ANGLET which is on the T grid !
      ! atmo variables are needed in T cell centers in subroutine stability,
      ! and are interpolated to the U grid later as necessary
      !-----------------------------------------------------------------
        !workx      = uatm(i,j)                ! wind velocity, m/s
        !worky      = vatm(i,j) 
        !uatm (i,j) = workx*cos(ANGLET(i,j))  & ! convert to POP grid
        !           + worky*sin(ANGLET(i,j))    ! note uatm, vatm, wind
        !vatm (i,j) = worky*cos(ANGLET(i,j))  & !  are on the T-grid here
        !           - workx*sin(ANGLET(i,j))

       enddo
      enddo
      
      else
        CALL FATAL_ERROR("ICE_LONGWAVE_TYPE should be PW or RM")
      end if

      end subroutine prepare_forcing

!c=======================================================================
!
!BOP
!
! !IROUTINE: interpolate_data 
!
! !INTERFACE:
!
      subroutine interpolate_data (field_data, field)
!
! !DESCRIPTION:
!
! Linear interpolation
!
! !REVISION HISTORY:
!
! authors Elizabeth C. Hunke, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind), dimension(ilo:ihi,jlo:jhi,2), intent(in) ::&
        field_data    ! 2 values used for interpolation

      real (kind=dbl_kind), dimension(ilo:ihi,jlo:jhi), intent(out) ::&
        field         ! interpolated field
!
!EOP
!
      integer (kind=int_kind) :: i,j

      do j=jlo,jhi
       do i=ilo,ihi
          field(i,j) = c1intp * field_data(i,j,1)      &
                     + c2intp * field_data(i,j,2) 
       enddo
      enddo

      end subroutine interpolate_data

!c=======================================================================
!
!BOP
!
! !IROUTINE: sss_clim - annual mean climatology for Levitus sss
!
! !INTERFACE:
!
!      subroutine sss_clim
!
! !DESCRIPTION:
!
! Creates annual mean climatology for Levitus sss from 12-month climatology.
!
! !REVISION HISTORY:
!
! authors Elizabeth C. Hunke, LANL
!
! !USES:
!
!      use ice_work, only:  worka
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
!      integer (kind=int_kind) :: i, j, nbits, k
!
!      nbits = 64                ! double precision data
!
!      sss_file = trim(ocn_data_dir)//'sss_Lev.mm'
!
!      if (my_task == master_task) then
!         write (nu_diag,*) ''
!         write (nu_diag,*) 'SSS climatology computed from:'
!         write (nu_diag,*) trim(sss_file)
!      endif

!      if (my_task == master_task) &
!          call ice_open (nu_forcing, sss_file, nbits)
!
!    !-------------------------------------------------------------------
!    ! create surface salinity climatology from monthly data
!    !-------------------------------------------------------------------
!
!      do j = jlo,jhi
!       do i = ilo,ihi
!        worka(i,j) = c0i
!       enddo
!      enddo

!      do k = 1,12   ! loop over 12 months
!
!         call ice_read (nu_forcing, k, worka, 'rda8', dbug)
!         do j = jlo,jhi
!          do i = ilo,ihi
!             sss(i,j) = worka(i,j) + sss(i,j)
!          enddo
!        enddo

!      enddo  ! k

!      do j = jlo,jhi
!       do i = ilo,ihi
!        sss(i,j) = sss(i,j) / c12       ! annual average salinity
!        sss(i,j) = max(sss(i,j), c0i)
!        Tf (i,j) = -depressT*sss(i,j) ! deg C 
!       enddo
!      enddo

!      ! close file
!      if (my_task == master_task) close(nu_forcing)

!      end subroutine sss_clim

!c=======================================================================
!
!BOP
!
! !IROUTINE: sst_ic - sst initial condition
!
! !INTERFACE:
!
      subroutine sst_ic
!
! !DESCRIPTION:
!
! Reads sst data for current month, and adjusts sst based on freezing 
! temperature.  Does not interpolate.
!
! !REVISION HISTORY:
!
! authors Elizabeth C. Hunke, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) :: i, j, nbits, k

      nbits = 64                ! double precision data

!      if (imt_global == 320) then   ! gx1
!         sst_file = trim(ocn_data_dir)//'sst_clim_hurrell.dat'
!      else                          ! gx3
!         sst_file = trim(ocn_data_dir)//'sst_Lev.mm'
!      endif

!      if (my_task == master_task) then
!         write (nu_diag,*) ''
!         write (nu_diag,*) 'SST initial condition:'
!         write (nu_diag,*) trim(sst_file)
!      endif

!      if (my_task == master_task)  &
!           call ice_open (nu_forcing, sst_file, nbits)
!
!      call ice_read (nu_forcing, month, sst, 'rda8', dbug)
!
!      if (my_task == master_task) close(nu_forcing)
!
!      do j=jlo,jhi
!       do i=ilo,ihi
!         ! Make sure sst is not less than Tf
!         sst(i,j) = max(sst(i,j),Tf(i,j))
!       enddo
!      enddo

      end subroutine sst_ic

!c=======================================================================
!
!BOP
!
! !IROUTINE: getflux_ocn_clim - interpolates sss, sst; restores sst
!
! !INTERFACE:
!
      subroutine getflux_ocn_clim
!
! !DESCRIPTION:
!
! Interpolate monthly sss, sst data to timestep.
! Note: Restoring is done only if sss_data_type and/or 
!       sst_data_type are set (not default) in namelist.
! 
! !REVISION HISTORY:
!
! authors: same as module
!
! !USES:
!
!      use ice_ocean
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) :: &
          i, j               &
      ,   imx,ipx            & ! record numbers for neighboring months
      ,   maxrec             & ! maximum record number
      ,   recslot            & ! spline slot for current record
      ,   midmonth            ! middle day of month

      real (kind=dbl_kind), dimension(ilo:ihi,jlo:jhi) ::& 
         sstdat              ! data value toward which SST is restored

      logical (kind=log_kind) :: readm

!      if (trim(sss_data_type) == 'clim') then
!         sss_file = trim(ocn_data_dir)//'sss_Lev.mm'
!      endif 
!      if (trim(sst_data_type) == 'clim') then
!         sst_file = trim(ocn_data_dir)//'sst_clim_hurrell.dat'
!      endif

!      if (my_task == master_task .and. istep == 1) then
!         write (nu_diag,*) ' '
!         if (trim(sss_data_type) == 'clim') then 
!            write (nu_diag,*) 'SSS data interpolated to timestep:' 
!            write (nu_diag,*) trim(sss_file) 
!         endif 
!         if (trim(sst_data_type) == 'clim') then 
!            write (nu_diag,*) 'SST data interpolated to timestep:' 
!            write (nu_diag,*) trim(sst_file) 
!            if (restore_sst) write (nu_diag,*) &
!            'SST restoring timescale = ',trestore,' days' 
!         endif
!      endif                     ! my_task, istep 

    !-------------------------------------------------------------------
    ! monthly data 
    !
    ! Assume that monthly data values are located in the middle of the 
    ! month.
    !-------------------------------------------------------------------
      
      midmonth = 15  ! data is given on 15th of every month
!      midmonth = fix(p5 * real(daymo(month),kind=dbl_kind))  ! exact middle

      ! Compute record numbers for surrounding months
      maxrec = 12
      imx  = mod(month+maxrec-2,maxrec) + 1
      ipx  = mod(month,         maxrec) + 1
      if (mday >= midmonth) imx = 99  ! other two points will be used
      if (mday <  midmonth) ipx = 99

      ! Determine whether interpolation will use values 1:2 or 2:3
      ! recslot = 2 means we use values 1:2, with the current value (2)
      !  in the second slot
      ! recslot = 1 means we use values 2:3, with the current value (2)
      !  in the first slot
!      recslot = 1                             ! latter half of month
!      if (mday < midmonth) recslot = 2        ! first half of month

      ! Find interpolation coefficients
!      call interp_coeff_monthly (recslot)

!      readm = .false.
!      if (istep==1 .or. (mday==midmonth .and. sec==0)) readm = .true.

    !------------------------------------------------------------------- 
    ! Read two monthly SSS values and interpolate. 
    ! Note: SSS is restored instantaneously to data. 
    !------------------------------------------------------------------- 
 
!      if (trim(sss_data_type) == 'clim') then 
!         call read_clim_data (readm, 0, imx, month, ipx, &
!                             sss_file, sss_data) 
!         call interpolate_data (sss_data, sss) 
!         do j = jlo,jhi 
!           do i = ilo,ihi 
!             sss(i,j) = max(sss(i,j), c0i) 
!             Tf (i,j) = -depressT*sss(i,j) ! deg C 
!           enddo 
!         enddo 
!      endif 

    !------------------------------------------------------------------- 
    ! Read two monthly SST values and interpolate. 
    ! Restore toward interpolated value. 
    ! Make sure SST is not below Tf. 
    !------------------------------------------------------------------- 
 
!      if (trim(sst_data_type) == 'clim') then 
!         call read_clim_data (readm, 0, imx, month, ipx,   &
!                             sst_file, sst_data) 
!         call interpolate_data (sst_data, sstdat) 
 
    !------------------------------------------------------------------- 
    ! Restore sst to data 
    !------------------------------------------------------------------- 
!        if (restore_sst) then
!         do j = jlo,jhi 
!          do i = ilo,ihi 
!            sst(i,j) = sst(i,j) + (sstdat(i,j)-sst(i,j))*dt/trest 
!          enddo 
!         enddo 
!        endif
!      endif 

      end subroutine getflux_ocn_clim

!c=======================================================================
!
!BOP
!
! !IROUTINE: getflux_ocn_ncar_init - reads data set
!
! !INTERFACE:
!
      subroutine getflux_ocn_ncar_init
!
! !DESCRIPTION:
!
! Reads NCAR pop ocean forcing data set 'pop_frc_gx1v3_010815.nc'
! 
! List of ocean forcing fields: Note that order is important!
! (order is determined by field list in vname).
! 
! For ocean mixed layer-----------------------------units 
! 
! 1  sst------temperature---------------------------(C)
! 2  sss------salinity------------------------------(ppt)
! 3  hbl------depth---------------------------------(m) 
! 4  u--------surface u current---------------------(m/s)
! 5  v--------surface v current---------------------(m/s)
! 6  dhdx-----surface tilt x direction--------------(m/m)
! 7  dhdy-----surface tilt y direction--------------(m/m)
! 8  qdp------ocean sub-mixed layer heat flux-------(W/m2) 
!
! Fields 4, 5, 6, 7 are on the U-grid; 1, 2, 3, and 8 are
! on the T-grid.
!
! !REVISION HISTORY:
!
! authors: Bruce Briegleb, NCAR
!          Elizabeth Hunke, LANL
!
! !USES:
!
      use ice_grid, only: t2ugrid
!      use ice_ocean
!      use ice_mpi_internal
!      use ice_exit
!      use ice_history, only: restart
      use ice_work, only: work_g1
!#ifdef CCSM
!      use shr_sys_mod, only : shr_sys_flush
!#endif
!      include "netcdf.inc"
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) :: &
       i, j     &
      , ni        &! field     index
      , mi         ! month     index

      character(len=16) ::&
       vname(nfld) ! variable names to search for on netCDF file
      data vname /                                         &
          'T',      'S',      'hblt',  'U',     'V',       &
          'dhdx',   'dhdy',   'qdp' /

      integer (kind=int_kind) ::&
        fid         &! file id for netCDF routines
      , dimid       &! dimension id for netCDF file
      , varid(nfld) &! variable id for field in netCDF file
      , ntim         ! number of times of data for netCDF file

      integer (kind=int_kind) ::  &
        status   &! status variable from netCDF routine 
      , nlat     &! number of longitudes of data for netCDF file
      , nlon     &! number of latitudes  of data for netCDF file
      , start(3) &! start location for netCDF data reads
      , count(3)  ! number of data items to read in

!      if (my_task == master_task) then

!         write (nu_diag,*) 'WARNING: evp_prep calculates surface tilt'
!         write (nu_diag,*) 'WARNING: stress from geostrophic currents,'
!         write (nu_diag,*) 'WARNING: not data from ocean forcing file.'
!         write (nu_diag,*) 'WARNING: Alter ice_dyn_evp.F if desired.'

!         if (restore_sst) write (nu_diag,*)                      &
!             'SST restoring timescale = ',trestore,' days' 
      
!         sst_file = trim(ocn_data_dir)//oceanmixed_file ! not just sst
      
        !---------------------------------------------------------------
        ! Read in ocean forcing data from an existing local netCDF file
        !---------------------------------------------------------------
!        write (nu_diag,*) 'ocean mixed layer forcing data file = ',   &
!                           sst_file

!        status = nf_open(sst_file, NF_NOWRITE, fid)
!        if (status /= NF_NOERR) then
!          call abort_ice ('ice: no netCDF file with ocn forcing data')
!        endif
!        write(nu_diag,*) 'Successful open of ocean forcing file'

!        status = nf_inq_dimid(fid,'nlon',dimid)
!        status = nf_inq_dimid(fid,'ni',dimid)
!        status = nf_inq_dimlen(fid,dimid,nlon)

!        status = nf_inq_dimid(fid,'nlat',dimid)
!        status = nf_inq_dimid(fid,'nj',dimid)
!        status = nf_inq_dimlen(fid,dimid,nlat)
!
!        if( nlon .ne. imt_global ) then
!          call abort_ice ('ice: ocn frc file nlon ne imt_global')
!        endif
!        if( nlat .ne. jmt_global ) then
!          call abort_ice ('ice: ocn frc file nlat ne jmt_global')
!        endif
!
!        do ni=1,nfld
!          status = nf_inq_varid(fid,vname(n),varid(n))
!          if (status /= NF_NOERR ) then
!            write(nu_diag,*) 'ice error- cannot find field ',vname(n)
!            call abort_ice ('ice: cannot find ocn frc field')
!          endif
!          write (nu_diag,*) 'Ocean forcing field found = ',vname(n)
!        enddo

!#ifdef CCSM
!        call shr_sys_flush(nu_diag)
!#endif

!        allocate (work_g1(imt_global,jmt_global))

!      endif ! master_task

      ! Read in ocean forcing data for all 12 months
!      do ni=1,nfld
!        do mi=1,12
                
!          if (my_task == master_task) then
!          start(1) = 1
!          start(2) = 1
!          start(3) = mi
!          count(1) = imt_global
!          count(2) = jmt_global
!          count(3) = 1

          ! Note: netCDF does single to double conversion if necessary
!          status = nf_get_vara_double(fid,varid(n),start,count,work_g1)
!          if (status /= NF_NOERR ) then
!            write(nu_diag,*) 'could not read in ocn forcing field ',&
!                                vname(ni)
!            call abort_ice ('ice read failed')
!          endif

!          endif ! master_task

!          call global_scatter(work_g1,ocn_frc_m(:,:,ni,mi)) 

!        enddo               ! month loop
!      enddo               ! field loop

!      if (my_task == master_task) then
!        status = nf_close(fid)
!        deallocate (work_g1)
!      endif

      end subroutine getflux_ocn_ncar_init

!c=======================================================================
!
!BOP
!
! !IROUTINE: getflux_ocn_ncar - interpolates data to timestep
!
! !INTERFACE:
!
      subroutine getflux_ocn_ncar
!
! !DESCRIPTION:
!
! Interpolate monthly ocean data to timestep.
! Restore sst if desired. sst is updated with surface fluxes in ice_ocean.F.
! 
! !REVISION HISTORY:
!
! authors: same as module
!
! !USES:
!
!      use ice_ocean
!      use ice_history, only: restart
      use ice_grid, only: hm
      use ice_work, only: worka
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) :: &
          i, j, ni             &
      ,   imx,ipx             &! record numbers for neighboring months
      ,   maxrec              &! maximum record number
      ,   recslot             &! spline slot for current record
      ,   midmonth            &! middle day of month
      ,   ng                   

      real (kind=dbl_kind), dimension(ilo:ihi,jlo:jhi) :: &
         sstdat              ! data value toward which SST is restored

    !-------------------------------------------------------------------
    ! monthly data 
    !
    ! Assume that monthly data values are located in the middle of the 
    ! month.
    !-------------------------------------------------------------------
      
      midmonth = 15  ! data is given on 15th of every month
!      midmonth = fix(p5 * real(daymo(month),kind=dbl_kind))  ! exact middle

      ! Compute record numbers for surrounding months
      maxrec = 12
      imx  = mod(month+maxrec-2,maxrec) + 1
      ipx  = mod(month,         maxrec) + 1
      if (mday >= midmonth) imx = 99  ! other two points will be used
      if (mday <  midmonth) ipx = 99

      ! Determine whether interpolation will use values 1:2 or 2:3
      ! recslot = 2 means we use values 1:2, with the current value (2)
      !  in the second slot
      ! recslot = 1 means we use values 2:3, with the current value (2)
      !  in the first slot
!      recslot = 1                             ! latter half of month
!      if (mday < midmonth) recslot = 2        ! first half of month

      ! Find interpolation coefficients
!      call interp_coeff_monthly (recslot)

!      do ni = nfld, 1, -1
        ! use sst_data arrays as temporary work space until n=1
!        if (imx /= 99) then  ! first half of month
!          sst_data(:,:,1) = ocn_frc_m(:,:,n,imx)
!          sst_data(:,:,2) = ocn_frc_m(:,:,n,month)
!        else                 ! second half of month
!          sst_data(:,:,1) = ocn_frc_m(:,:,n,month)
!          sst_data(:,:,2) = ocn_frc_m(:,:,n,ipx)
!        endif
!        call interpolate_data (sst_data,worka)
        ! masking by hm is necessary due to NaNs in the data file
!        do j = jlo,jhi 
!          do i = ilo,ihi 
!            if (ni == 2) sss    (i,j) = c0i
!            if (ni == 3) hmix   (i,j) = c0i
!            if (ni == 4) uocn   (i,j) = c0i
!            if (ni == 5) vocn   (i,j) = c0i
!            if (ni == 6) ss_tltx(i,j) = c0i
!            if (ni == 7) ss_tlty(i,j) = c0i
!            if (ni == 8) qdp    (i,j) = c0i
!            if (hm(i,j) == c1i) then
!              if (ni == 2) sss    (i,j) = worka(i,j)
!              if (ni == 3) hmix   (i,j) = worka(i,j)
!              if (ni == 4) uocn   (i,j) = worka(i,j)
!              if (ni == 5) vocn   (i,j) = worka(i,j)
!              if (ni == 6) ss_tltx(i,j) = worka(i,j)
!              if (ni == 7) ss_tlty(i,j) = worka(i,j)
!              if (ni == 8) qdp    (i,j) = worka(i,j)
!
              ! debugging
!              if (ni == 3 .and. hmix(i,j) == c0i) then
!                print*,i,j,hm(i,j),hmix(i,j),sss(i,j),uocn(i,j),qdp(i,j)
!                stop
!              endif
!            endif
!          enddo
!        enddo
!      enddo

      do j = jlo,jhi 
        do i = ilo,ihi 
!          sss (i,j) = max (sss(i,j), c0i) 
!          Tf  (i,j) = -depressT*sss(i,j) ! deg C 
!          hmix(i,j) = max(hmix(i,j), c0i) 
        enddo 
      enddo 

!      if (restore_sst) then
!        do j = jlo,jhi 
!         do i = ilo,ihi 
!           sst(i,j) = sst(i,j) + (sstdat(i,j)-worka(i,j))*dt/trest 
!         enddo 
!        enddo 
!     else sst is only updated in ice_ocean.F
!      endif

      ! initialize sst properly on first step
!      if (istep1 <= 1 .and. .not. (restart)) then
!        call interpolate_data (sst_data,sst)
!        do j = jlo,jhi 
!          do i = ilo,ihi 
!            if (hm(i,j) == c1i) then
!              sst(i,j) =  max (sst(i,j), Tf(i,j)) 
!            else
!              sst(i,j) = c0i
!            endif
!          enddo 
!        enddo 
!      endif

!      if (dbug) then
!         if (my_task == master_task)  & 
!              write (nu_diag,*) 'getflux_ocn_ncar'
!         ng = (ihi-ilo+1)*(jhi-jlo+1)
!         call ice_global_real_minmax(ng,Tf,     'Tf      ')
!         call ice_global_real_minmax(ng,sst,    'sst     ')
!         call ice_global_real_minmax(ng,sss,    'sss     ')
!         call ice_global_real_minmax(ng,hmix,   'hmix    ')
!         call ice_global_real_minmax(ng,uocn,   'uocn    ')
!         call ice_global_real_minmax(ng,vocn,   'vocn    ')
!         call ice_global_real_minmax(ng,ss_tltx,'tltx    ')
!         call ice_global_real_minmax(ng,ss_tlty,'tlty    ')
!         call ice_global_real_minmax(ng,qdp,    'qdp     ')
!      endif

      end subroutine getflux_ocn_ncar

!c=======================================================================

      end module ice_flux_in

!c=======================================================================
