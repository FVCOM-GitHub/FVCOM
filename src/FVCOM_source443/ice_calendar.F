!/===========================================================================/
! CVS VERSION INFORMATION
! $Id$
! $Name$
! $Revision$
!/===========================================================================/

!=======================================================================
!BOP
!
! !MODULE: ice_calendar - calendar routines for managing time 
!
! !DESCRIPTION:
!
! Calendar routines for managing time
!
! !REVISION HISTORY:
!
! authors: Elizabeth C. Hunke, LANL
!          Tony Craig, NCAR
!
! !INTERFACE:
!
      module ice_calendar
!
! !USES:
!
      use ice_constants


!
!EOP
!
      implicit none
      save

      integer (kind=int_kind) ::&
         daymo(12)              &  ! number of days in each month
      ,  daycal(13)               ! day number at end of month

      data daymo /   31,28,31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
      data daycal/ 0,31,59,90,120,151,181,212,243,273,304,334,365/

      integer (kind=int_kind) ::  &
         istep     &! local step counter for time loop
      ,  istep0    &! step counter, number of steps taken in previous run
      ,  istep1    &! step counter, number of steps at current timestep
      ,  mday      &! day of the month
      ,  week      &! week of the year
      ,  month     &! 1 \le month \le 12
      ,  monthp    &! last month
      ,  year_init &! initial year
      ,  nyr       &! year number
      ,  idate     &! date (yyyymmdd)
      ,  sec       &! elapsed seconds into date
      ,  npt       &! total number of time steps (dt) 
      ,  ndyn_dt   &! reduced timestep for dynamics: ndyn_dt=dt/dyn_dt
      ,  stop_now           &! if 1, end program execution
      ,  write_restart      &! if 1, write restart now
      ,  cpl_write_history  &! if 1, write history on command from cpl
      ,  diagfreq           &! diagnostic output frequency (10 = once per 10 dt)
      ,  dumpfreq_n          ! restart output frequency (10 = once per 10 d,m,y)

      real (kind=dbl_kind) ::& 
!         dt                & ! thermodynamics timestep (s)
         dtice,                & ! thermodynamics timestep (s)
!       defined in old ice model  dtice
        dyn_dt            & ! dynamics/transport/ridging timestep (s)
      ,  dtei              & ! 1/dte, where dte is the EVP timestep (1/s)
!      ,  time              & ! total elapsed time (s)
!      ggao conflict
      ,  time_forc         & ! time of last forcing update (s)
      ,  yday               ! day of the year

      logical (kind=log_kind) ::  &
         new_year           &! new year = .true.
      ,  new_month          &! new month = .true.
      ,  new_week           &! new week = .true.
      ,  new_day            &! new day = .true.
      ,  write_ic           &! write initial condition now
      ,  write_history       ! write history now
                             
      character (len=1) ::  & 
         histfreq           &! history output frequency, 'y','m','d','1'
      ,  dumpfreq            ! restart frequency, 'y','m','d'

!=======================================================================

      contains

!=======================================================================
!BOP
!
! !IROUTINE: init_calendar - initialize calendar variables
!
! !INTERFACE:
!
      subroutine init_calendar
!
! !DESCRIPTION:
!
! Initialize calendar variables
!
! !REVISION HISTORY:
!
! authors: Elizabeth C. Hunke, LANL
!          Tony Craig, NCAR
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      istep1 = istep0   ! number of steps at current timestep
                        ! real (dumped) or imagined (use to set calendar)
      istep = 0         ! local timestep number
!      time=istep0*dt    ! s
!      time=istep0*dtice    ! s
!       ggao be careful 

      yday=c0i           ! absolute day number
      mday=0            ! day of the month
      month=0           ! month
      nyr=0             ! year
      idate=00000101    ! date
      sec=0             ! seconds into date
      stop_now = 0      ! end program execution if stop_now=1
!      dyn_dt = dt/real(ndyn_dt,kind=dbl_kind) ! dynamics et al timestep
      dyn_dt = dtice/real(ndyn_dt,kind=dbl_kind) ! dynamics et al timestep
!      ggao 
      end subroutine init_calendar

!=======================================================================
!BOP
!
! !IROUTINE: calendar - computes date at the end of the time step
!
! !INTERFACE:
!
      subroutine calendar(ttime)
!
! !DESCRIPTION:
!
! Determine the date at the end of the time step
!
! !REVISION HISTORY:
!
! authors: Elizabeth C. Hunke, LANL
!          Tony Craig, NCAR
!
! !USES:
      use ice_fileunits
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind), intent(in) :: &
         ttime                          ! time variable
!
!EOP
!
      integer (kind=int_kind) ::     &
         k                           &
      ,  nyrp,mdayp,weekp            &  ! previous year, day, week
      ,  elapsed_days                &  ! since beginning this run
      ,  elapsed_months                 ! since beginning this run
                                      
      real (kind=dbl_kind) ::        &
         tday                        &  ! absolute day number
      ,  dayyr                          ! number of days per year

      dayyr = 365.0_dbl_kind

      nyrp=nyr
      monthp=month
      weekp=week
      mdayp=mday
      new_year=.false.
      new_month=.false.
      new_week=.false.
      new_day=.false.
      write_history=.false.
      write_restart=0

      sec = mod(ttime,secday)           ! elapsed seconds into date at 
                                        ! end of dt
      tday = (ttime-sec)/secday + c1i    ! absolute day number 
      yday = mod(tday-c1i,dayyr) + c1i   ! day of the year 
      week = int(yday/c7i) + c1i         ! week of the year
      do k = 1, 12
        if (yday > real(daycal(k),kind=dbl_kind)) month = k ! month
      enddo
      mday = int(yday) - daycal(month)  ! day of the month
      nyr = int((tday-c1i)/dayyr) + 1    ! year number
      elapsed_months = (nyr - 1)*12 + month - 1
      elapsed_days = int(tday) - 1 

      idate = (nyr+year_init-1)*10000 + month*100 + mday ! date (yyyymmdd) 

      if (istep >= npt+1) stop_now = 1
      if (nyr /= nyrp) new_year = .true.
      if (month /= monthp) new_month = .true.
      if (week /= weekp) new_week = .true.
      if (mday /= mdayp) new_day = .true.

      if (histfreq == '1') write_history=.true.
      if (istep > 1) then
        if (((histfreq == 'y'.or.histfreq == 'Y').and.new_year).or.  & 
            ((histfreq == 'm'.or.histfreq == 'M').and.new_month).or. &
            ((histfreq == 'w'.or.histfreq == 'W').and.new_week).or.  &
            ((histfreq == 'd'.or.histfreq == 'D').and.new_day))      &
            write_history=.true.
        select case (dumpfreq)
        case ("y", "Y")
          if (new_year  .and. mod(nyr, dumpfreq_n)==0)           &
                write_restart = 1                                 
        case ("m", "M")                                           
          if (new_month .and. mod(elapsed_months,dumpfreq_n)==0) &
                write_restart=1                                   
        case ("d", "D")                                           
          if (new_day   .and. mod(elapsed_days, dumpfreq_n)==0)  &
                write_restart = 1                                 
        end select                                                
      endif                                                       
                                                                  
      if (my_task == master_task.and.mod(istep,diagfreq) == 0    &
        .and.stop_now /= 1) then
        write(nu_diag,*) ' '
        write(nu_diag,'(a7,i10,4x,a6,i10,4x,a4,i10)')            &
            'istep1:', istep1, 'idate:', idate, 'sec:', sec
      endif

      end subroutine calendar

!=======================================================================

      end module ice_calendar

!=======================================================================
