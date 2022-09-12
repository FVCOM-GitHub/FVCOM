MODULE MOD_binary
  USE MOD_BULK
  USE MOD_SURFACEFORCE
  USE MOD_SET_TIME
  implicit none
  
  TYPE bin_data
     type(TIME) :: dtm
     REAL(SP), POINTER :: data(:)
  END TYPE bin_data
  
  TYPE(bin_data), POINTER :: SWR_NEXT, SWR_PREV
  TYPE(bin_data), POINTER :: NHF_NEXT, NHF_PREV
  TYPE(bin_data), POINTER :: LWR_NEXT, LWR_PREV
  TYPE(bin_data), POINTER :: WNDX_NEXT, WNDX_PREV
  TYPE(bin_data), POINTER :: WNDY_NEXT, WNDY_PREV
  TYPE(bin_data), POINTER :: EVP_NEXT, EVP_PREV
  TYPE(bin_data), POINTER :: PRC_NEXT, PRC_PREV

  TYPE(bin_data), POINTER :: AIP_NEXT, AIP_PREV
  TYPE(bin_data), POINTER :: SAT_NEXT, SAT_PREV
  TYPE(bin_data), POINTER :: SPQ_NEXT, SPQ_PREV
  TYPE(bin_data), POINTER :: RH_NEXT, RH_PREV
  TYPE(bin_data), POINTER :: CLD_NEXT, CLD_PREV
! EJA <
  TYPE(bin_data), POINTER :: DPT_NEXT, DPT_PREV
! EJA >

  !  STORAGE FOR OUTPUT OF BINARY FILES
  REAL(SP), POINTER :: WNDYGL(:), WNDXGL(:)
  REAL(SP), POINTER :: SWRGL(:),  NHFGL(:),  LWRGL(:)
  REAL(SP), POINTER :: PRCGL(:),  EVPGL(:)

  REAL(SP), POINTER :: AIPGL(:), SATGL(:), SPQGL(:), CLDGL(:), RHGL(:), DPTGL(:)

  LOGICAL ::SAT_ON, SPQ_ON, CLD_ON, DPT_ON

  INTEGER IERR
  
CONTAINS

  SUBROUTINE UPDATE_BINARY(NOW)
    IMPLICIT NONE
    TYPE(TIME) :: NOW
    
    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START: UPDATE_BINARY"

    IF(PRECIPITATION_ON) THEN
       CALL UPDATE_EVP(NOW,EVP,PRC)
    END IF

    IF(WIND_ON) THEN

       IF(WIND_TYPE == 'stress') THEN
          CALL UPDATE_WND(NOW,STRX,STRY)
       ELSE IF(WIND_TYPE == 'speed') THEN
          CALL UPDATE_WND(NOW,SPDX,SPDY)
          CALL PSIMPLE_DRAG(SPDX,SPDY,STRX,STRY)
       END IF
       
!    WRITE(IPT,*) "MIN/MAX(SPDX)::",minval(spdx),maxval(spdx)
!    WRITE(IPT,*) "MIN/MAX(SPDy)::",minval(spdy),maxval(spdy)

!    WRITE(IPT,*) "MIN/MAX(STRX)::",minval(strx),maxval(strx)
!    WRITE(IPT,*) "MIN/MAX(STRy)::",minval(stry),maxval(stry)

    END IF

    IF(HEATING_ON) THEN
       CALL UPDATE_HFX(NOW,SWR,NHF,LWR)
    END IF

    IF(AIRPRESSURE_ON) THEN
       CALL UPDATE_AIP(NOW,AIP)
    END IF

    IF(ICE_MODEL) THEN
       CALL UPDATE_ICE(NOW,SAT,SPQ,CLD,RH)
    END IF
! EJA <
    IF(HEATING_SOLAR) THEN
       CALL UPDATE_SOLAR(NOW,SAT,DPT,CLD)
    END IF
! EJA >
    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END: UPDATE_BINARY"
    
  END SUBROUTINE UPDATE_BINARY


  FUNCTION NEW_DATA(dims)
    IMPLICIT NONE
    TYPE(bin_data), POINTER :: NEW_DATA
    integer, intent(IN) :: DIMS
    INTEGER :: STATUS

    ALLOCATE(NEW_DATA,stat=status)
    IF(status /=0) CALL FATAL_ERROR("NEW_DATA: COULD NOT ALLOCATE TYPE POINTER?")
    
    ALLOCATE(NEW_DATA%DATA(0:DIMS), STAT=STATUS)
    IF(status /=0) CALL FATAL_ERROR("NEW_DATA: COULD NOT ALLOCATE DATA POINTER?")
    
  END FUNCTION NEW_DATA


  SUBROUTINE LOAD_BINARY(WND,HFX,EVP,AIPF,SATF,SPQF,CLDF,DPTF)
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: WND,HFX,EVP,AIPF,SATF,SPQF,CLDF,DPTF
    INTEGER :: STATUS
    INTEGER(ITIME) :: dummy
    CHARACTER(LEN=4) :: FLAG

    
    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START:LOAD_BINARY"
    
    IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) "Looking for Wind Stress file:"&
         &//TRIM(WND)
    inquire(file=trim(WND),exist=WIND_ON)
    IF(WIND_ON) THEN
       
       IF(WIND_TYPE /= 'stress' .and. WIND_TYPE/='speed')THEN
          CALL FATAL_ERROR("To convert a binary wind file,",&
               &"you must specify 'binary speed' or 'binary stress'")
       END IF
       
       IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) "FOUND WIND FILE: OPEN AND READ"
       
       IF(MSR) CALL FOPEN(WNDUNIT,WND,'cfr')
       
       WNDX_NEXT => NEW_DATA(N)
       WNDX_PREV => NEW_DATA(N)
       WNDY_NEXT => NEW_DATA(N)
       WNDY_PREV => NEW_DATA(N)
       
       CALL READ_WND(WNDX=WNDX_PREV,WNDY=WNDY_PREV)
       CALL READ_WND(WNDX=WNDX_NEXT,WNDY=WNDY_NEXT)
       
       IF(DBG_SET(DBG_LOG)) THEN
          CALL PRINT_REAL_TIME(WNDX_PREV%dtm,IPT,"FIRST TIME POINT",timezone)
          CALL PRINT_REAL_TIME(WNDX_NEXT%dtm,IPT,"SECOND TIME POINT",timezone)
       END IF
       
       IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) "FOUND WIND FILE: READ FIRST DATA POINTS"
       
    ELSE
       
       IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) "! NO WIND FILE FOUND"
       
    END IF
    
    IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) "Looking for HeatFlux file:"&
         &//TRIM(HFX)
    inquire(file=trim(HFX),exist=HEATING_ON)
    IF(HEATING_ON) THEN
       
       IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) "FOUND HEATING FILE: OPEN AND READ"
       
       IF(MSR) CALL FOPEN(HFXUNIT,HFX,'cfr')
       
       SWR_NEXT => NEW_DATA(M)
       SWR_PREV => NEW_DATA(M)
       NHF_NEXT => NEW_DATA(M)
       NHF_PREV => NEW_DATA(M)
       LWR_NEXT => NEW_DATA(M)
       LWR_PREV => NEW_DATA(M)

       CALL READ_HFX(SWR=SWR_PREV,NHF=NHF_PREV,LWR=LWR_PREV)
       CALL READ_HFX(SWR=SWR_NEXT,NHF=NHF_NEXT,LWR=LWR_NEXT)

       IF(DBG_SET(DBG_LOG)) THEN
          CALL PRINT_REAL_TIME(SWR_PREV%dtm,IPT,"FIRST TIME POINT",timezone)
          CALL PRINT_REAL_TIME(SWR_NEXT%dtm,IPT,"SECOND TIME POINT",timezone)
       END IF
       
       IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) "FOUND HEATING FILE: READ FIRST DATA POINTS"
    ELSE
       
       IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) "! NO HEATING FILE FOUND"
       
       
    END IF
    
    IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) "Looking for Evaporation file:"&
         &//TRIM(EVP)
    inquire(file=trim(EVP),exist=PRECIPITATION_ON)
    IF(PRECIPITATION_ON) THEN
 
       IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) "FOUND PRECIPITATION FILE: OPEN AND READ"
       
       IF(MSR) CALL FOPEN(EVPUNIT,EVP,'cur')

       EVP_NEXT => NEW_DATA(M)
       EVP_PREV => NEW_DATA(M)
       PRC_NEXT => NEW_DATA(M)
       PRC_PREV => NEW_DATA(M)

       CALL READ_EVP(EVP=EVP_PREV,PRC=PRC_PREV)
       CALL READ_EVP(EVP=EVP_NEXT,PRC=PRC_NEXT)
       
       IF(DBG_SET(DBG_LOG)) THEN
          CALL PRINT_REAL_TIME(SWR_PREV%dtm,IPT,"FIRST TIME POINT",timezone)
          CALL PRINT_REAL_TIME(SWR_NEXT%dtm,IPT,"SECOND TIME POINT",timezone)
       END IF
       
       IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) "FOUND PRECIPITATION FILE: READ FIRST DATA POINTS"
       
    ELSE
       
       IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) "! NO PRECIPITATION FILE FOUND"
       
    END IF


# if defined (AIR_PRESSURE)
    IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) "Looking for Air Pressure file:"&
         &//TRIM(AIPF)
    inquire(file=trim(AIPF),exist=AIRPRESSURE_ON)
    IF(AIRPRESSURE_ON) THEN
 
       IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) "FOUND AIR PRESSURE FILE: OPEN AND READ"
       
       IF(MSR) CALL FOPEN(AIPUNIT,AIPF,'cfr')

       AIP_NEXT => NEW_DATA(M)
       AIP_PREV => NEW_DATA(M)

       CALL READ_AIP(AIP=AIP_PREV)
       CALL READ_AIP(AIP=AIP_NEXT)
       
       IF(DBG_SET(DBG_LOG)) THEN
          CALL PRINT_REAL_TIME(AIP_PREV%dtm,IPT,"FIRST TIME POINT",timezone)
          CALL PRINT_REAL_TIME(AIP_NEXT%dtm,IPT,"SECOND TIME POINT",timezone)
       END IF
       
       IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) "FOUND AIR PRESSURE FILE: READ FIRST DATA POINTS"
       
    ELSE
       
       IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) "! NO AIR PRESSURE FILE FOUND"
       
    END IF
# endif


# if defined (ICE)
    IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) "Looking for Air Temperature file:"&
         &//TRIM(SATF)
    inquire(file=trim(SATF),exist=SAT_ON)
    IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) "Looking for Specific Humidity file:"&
         &//TRIM(SPQF)
    inquire(file=trim(SPQF),exist=SPQ_ON)
    IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) "Looking for Cloud Cover file:"&
         &//TRIM(CLDF)
    inquire(file=trim(CLDF),exist=CLD_ON)
    
    IF(SAT_ON .AND. SPQ_ON .AND. CLD_ON) ICE_MODEL = .TRUE.
    
    IF(ICE_MODEL) THEN
 
       IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) "FOUND ICE MODEL FILES: OPEN AND READ"
       
       IF(MSR) CALL FOPEN(SATUNIT,SATF,'cfr')
       IF(MSR) CALL FOPEN(SPQUNIT,SPQF,'cfr')
       IF(MSR) CALL FOPEN(CLDUNIT,CLDF,'cfr')

       SAT_NEXT => NEW_DATA(M)
       SAT_PREV => NEW_DATA(M)
       SPQ_NEXT => NEW_DATA(M)
       SPQ_PREV => NEW_DATA(M)
       CLD_NEXT => NEW_DATA(M)
       CLD_PREV => NEW_DATA(M)
       RH_NEXT => NEW_DATA(M)
       RH_PREV => NEW_DATA(M)

       CALL READ_ICE(SAT=SAT_PREV,SPQ=SPQ_PREV,CLD=CLD_PREV,RH=RH_PREV)
       CALL READ_ICE(SAT=SAT_NEXT,SPQ=SPQ_NEXT,CLD=CLD_NEXT,RH=RH_NEXT)
       
       IF(DBG_SET(DBG_LOG)) THEN
          CALL PRINT_REAL_TIME(SAT_PREV%dtm,IPT,"FIRST TIME POINT",timezone)
          CALL PRINT_REAL_TIME(SAT_NEXT%dtm,IPT,"SECOND TIME POINT",timezone)
       END IF
       
       IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) "FOUND ICE MODEL FILES: READ FIRST DATA POINTS"
       
    ELSE
       
       IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) "! NO ICE MODEL FILES FOUND OR SOME OF FILES MISSING"
       
    END IF

# endif

! EJA <
# if defined (HEATING_SOLAR)
    IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) "Looking for Air Temperature file:"&
         &//TRIM(SATF)
    inquire(file=trim(SATF),exist=SAT_ON)
    IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) "Looking for Dew Point file:"&
         &//TRIM(DPTF)
    inquire(file=trim(DPTF),exist=DPT_ON)
    IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) "Looking for Cloud Cover file:"&
         &//TRIM(CLDF)
    inquire(file=trim(CLDF),exist=CLD_ON)
    
    IF(SAT_ON .AND. DPT_ON .AND. CLD_ON) HEATING_SOLAR_ON = .TRUE.
    
    IF(HEATING_SOLAR_ON) THEN
 
       IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) "FOUND HEATING SOLAR MODEL FILES: OPEN AND READ"
       
       IF(MSR) CALL FOPEN(SATUNIT,SATF,'cfr')
       IF(MSR) CALL FOPEN(DPTUNIT,DPTF,'cfr')
       IF(MSR) CALL FOPEN(CLDUNIT,CLDF,'cfr')

       SAT_NEXT => NEW_DATA(M)
       SAT_PREV => NEW_DATA(M)
       DPT_NEXT => NEW_DATA(M)
       DPT_PREV => NEW_DATA(M)
       CLD_NEXT => NEW_DATA(M)
       CLD_PREV => NEW_DATA(M)

       CALL READ_SOLAR(SAT=SAT_PREV,DPT=DPT_PREV,CLD=CLD_PREV)
       CALL READ_SOLAR(SAT=SAT_NEXT,DPT=DPT_NEXT,CLD=CLD_NEXT)
       
       IF(DBG_SET(DBG_LOG)) THEN
          CALL PRINT_REAL_TIME(SAT_PREV%dtm,IPT,"FIRST TIME POINT",timezone)
          CALL PRINT_REAL_TIME(SAT_NEXT%dtm,IPT,"SECOND TIME POINT",timezone)
       END IF
       
       IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) "FOUND HEATING SOLAR MODEL FILES: READ FIRST DATA POINTS"
       
    ELSE
       
       IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) "! NO HEATING SOLAR MODEL FILES FOUND OR SOME OF FILES MISSING"
       
    END IF

# endif
! EJA >

    IF (.not. PRECIPITATION_ON .and. .not. HEATING_ON .and. &
      & .not. WIND_ON .and. .not. AIRPRESSURE_ON .and. .not. ICE_MODEL) &
      &  CALL FATAL_ERROR("FOUND NO BINARY FORCING INPUT FILES?")
    
    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END:LOAD_BINARY"
    
  END SUBROUTINE LOAD_BINARY
!------------------------------------------------------------------
  SUBROUTINE IOERROR(IOS,MSG)
    IMPLICIT NONE
    INTEGER IOS
    CHARACTER(LEN=*) MSG
    CHARACTER(LEN=4) IOSC
    
    IF(IOS ==0) RETURN

    WRITE(IOSC,'(I4)') IOS
    
    CALL FATAL_ERROR("ERROR DURING FILE IO:"//TRIM(IOSC),&
         TRIM(MSG))
  
  END SUBROUTINE IOERROR

!-------------------------------------------------------------------
  SUBROUTINE READ_WND(WNDX,WNDY)
    IMPLICIT NONE
    TYPE(bin_data) :: WNDX,WNDY
    REAL(SP) :: hour,yr,day,iter
    integer :: i, SOURCE, ios
    Real(SP), POINTER :: WNDYGL(:),WNDXGL(:)
   
    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START: READ_WND"


    IF(MSR) THEN
 
       IF(PAR) THEN
          ALLOCATE(WNDYGL(NGL))
          ALLOCATE(WNDXGL(NGL))
       ELSE
          WNDXGL => WNDX%DATA(1:NGL)
          WNDYGL => WNDY%DATA(1:NGL)
       END IF

!	*** EJA Edit 1/11/2011 - alter for formatted (ascii) input file
       !READ(UNIT=wndunit,IOSTAT=ios) hour
       READ(wndunit,*,IOSTAT=ios) hour,yr,day,iter
       CALL IOERROR(IOS,"Can't read hour from wind file")
!	*** EJA Edit 1/11/2011 - alter for formatted (ascii) input file
       !READ(UNIT= wndunit,IOSTAT=ios)(WNDXGL(i),WNDYGL(i),i=1,ngl)
       DO I=1,NGL
          READ(wndunit,*,IOSTAT=ios) WNDXGL(i),WNDYGL(i)
       END DO
       CALL IOERROR(IOS,"Can't read data from wind file")
       

       IF(DBG_SET(DBG_IO)) THEN
          WRITE(IPT,*) "MIN/MAX/MEAN(WNDX)::",MINVAL(WNDXgl),maxval(WNDXgl),sum(WNDXgl)/real(Ngl)
          WRITE(IPT,*) "MIN/MAX/MEAN(WNDY)::",MINVAL(WNDYgl),maxval(WNDYgl),sum(WNDYgl)/real(Ngl)
       END IF

    END IF


    ! IF NOT PAR, VARIABLES ARE ALREADY POINTED CORRECTLY
    IF (PAR)THEN
# if defined(MULTIPROCESSOR)
       
       SOURCE = MSRID -1

       CALL MPI_BCAST(hour,1,MPI_F,SOURCE,MPI_FVCOM_GROUP,IERR)
       CALL PDEAL(MYID,MSRID,NPROCS,EMAP,WNDXGL,WNDX%DATA)
       CALL PDEAL(MYID,MSRID,NPROCS,EMAP,WNDYGL,WNDY%DATA)

       IF(MSR) THEN
          DEALLOCATE(WNDXGL)
          DEALLOCATE(WNDYGL)
       END IF
# endif
    END IF

    NULLIFY(WNDXGL,WNDYGL)

    WNDX%dtm = seconds2time(hour*3600.0_SP) + ZEROTIME

    WNDY%dtm = WNDX%dtm

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END: READ_WND"
  END SUBROUTINE READ_WND

!-------------------------------------------------------------------
  SUBROUTINE READ_HFX(SWR,NHF,LWR)
    IMPLICIT NONE
    TYPE(bin_data) :: SWR,NHF,LWR
    REAL(SP) :: hour,yr,day,iter
    integer :: i, SOURCE, ios
    Real(SP), POINTER :: NHFGL(:),SWRGL(:),LWRGL(:)
    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START: READ_HFX"
    

    IF(MSR) THEN
       IF(PAR) THEN
          ALLOCATE(NHFGL(MGL))
          ALLOCATE(SWRGL(MGL))
          ALLOCATE(LWRGL(MGL))
       ELSE
          SWRGL => SWR%DATA(1:MGL)
          NHFGL => NHF%DATA(1:MGL)
          LWRGL => LWR%DATA(1:MGL)
       END IF



       READ(hfxunit,*,IOSTAT=ios) hour,yr,day,iter
       CALL IOERROR(IOS,"Can't read hour from heatflux file")
! *** EJA Edit 1/11/2011 - adjust for Schwab input files
!       READ(UNIT=hfxunit,IOSTAT=ios)(NHFGL(i),SWRGL(i),i=1,mgl)
       DO I=1,MGL
          READ(hfxunit,*,IOSTAT=ios) NHFGL(i),SWRGL(i),LWRGL(i)
       END DO
       CALL IOERROR(IOS,"Can't read data from heatflux file")

       IF(DBG_SET(DBG_IO)) THEN
          WRITE(IPT,*) "MIN/MAX/MEAN(SWR)::",MINVAL(swrgl),maxval(swrgl),sum(swrgl)/real(Mgl)
          WRITE(IPT,*) "MIN/MAX/MEAN(NHF)::",MINVAL(NHFgl),maxval(NHFgl),sum(NHFgl)/real(Mgl)
          WRITE(IPT,*) "MIN/MAX/MEAN(LWR)::",MINVAL(LWRgl),maxval(LWRgl),sum(LWRgl)/real(Mgl)
       END IF


    END IF


    ! IF NOT PAR, VARIABLES ARE ALREADY POINTED CORRECTLY
    IF (PAR)THEN
# if defined(MULTIPROCESSOR)
       
       SOURCE = MSRID -1

       CALL MPI_BCAST(hour,1,MPI_F,SOURCE,MPI_FVCOM_GROUP,IERR)
       CALL PDEAL(MYID,MSRID,NPROCS,NMAP,SWRGL,SWR%DATA)
       CALL PDEAL(MYID,MSRID,NPROCS,NMAP,NHFGL,NHF%DATA)
       CALL PDEAL(MYID,MSRID,NPROCS,NMAP,LWRGL,LWR%DATA)

       IF(MSR) THEN
          DEALLOCATE(SWRGL)
          DEALLOCATE(NHFGL)
          DEALLOCATE(LWRGL)
       END IF
# endif
    END IF

    NULLIFY(SWRGL,NHFGL,LWRGL)

    SWR%dtm = seconds2time(hour*3600.0_SP) + ZEROTIME

    NHF%dtm = SWR%dtm
    
    LWR%dtm = SWR%dtm
    
    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END: READ_HFX"

  END SUBROUTINE READ_HFX

!-------------------------------------------------------------------
  SUBROUTINE READ_EVP(EVP,PRC)
    IMPLICIT NONE
    TYPE(bin_data) :: EVP,PRC
    REAL(SP) :: hour
    integer :: i, SOURCE,ios
    Real(SP), POINTER :: EVPGL(:),PRCGL(:)
  
    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START: READ_EVP"

    IF(MSR) THEN

         IF(PAR) THEN
          ALLOCATE(EVPGL(MGL))
          ALLOCATE(PRCGL(MGL))
       ELSE
          EVPGL => EVP%DATA(1:MGL)
          PRCGL => PRC%DATA(1:MGL)
       END IF

       READ(UNIT=evpunit,IOSTAT=ios) hour
       CALL IOERROR(IOS,"Can't read hour from evap/prec file")
       READ(UNIT=evpunit,IOSTAT=ios)(evpgl(i),prcgl(i),i=1,mgl)
       CALL IOERROR(IOS,"Can't read data from evap/prec file")

       IF(DBG_SET(DBG_IO)) THEN
          WRITE(IPT,*) "MIN/MAX/MEAN(EVAP)::",MINVAL(evpgl),maxval(evpgl),sum(evpgl)/real(Mgl)
          WRITE(IPT,*) "MIN/MAX/MEAN(PREC)::",MINVAL(prcgl),maxval(prcgl),sum(prcgl)/real(Mgl)
       END IF

    END IF


    ! IF NOT PAR, VARIABLES ARE ALREADY POINTED CORRECTLY
    IF (PAR)THEN
# if defined(MULTIPROCESSOR)
       
       SOURCE = MSRID -1
       CALL MPI_BCAST(hour,1,MPI_F,SOURCE,MPI_FVCOM_GROUP,IERR)
       CALL PDEAL(MYID,MSRID,NPROCS,NMAP,prcgl,PRC%DATA)
       CALL PDEAL(MYID,MSRID,NPROCS,NMAP,evpgl,EVP%DATA)
       
       IF(MSR) THEN
          DEALLOCATE(PRCGL)
          DEALLOCATE(EVPGL)
       END IF

# endif
    END IF

    NULLIFY(EVPGL,PRCGL)
    
    PRC%dtm = seconds2time(hour*3600.0_SP) + ZEROTIME

    EVP%dtm = PRC%dtm

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END: READ_EVP"

  END SUBROUTINE READ_EVP

  SUBROUTINE READ_AIP(AIP)
    IMPLICIT NONE
    TYPE(bin_data) :: AIP
    REAL(SP) :: hour
    integer :: i, SOURCE, ios
    Real(SP), POINTER :: AIPGL(:)
    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START: READ_AIP"
    

    IF(MSR) THEN
       IF(PAR) THEN
          ALLOCATE(AIPGL(MGL))
       ELSE
          AIPGL => AIP%DATA(1:MGL)
       END IF
!*** Edit by EJA Anderson - 2/24/2011 for formatted pressure input ***
!       READ(UNIT=aipunit,IOSTAT=ios) hour
       READ(aipunit,*,IOSTAT=ios) hour
       CALL IOERROR(IOS,"Can't read hour from air pressure file")
!       READ(UNIT=aipunit,IOSTAT=ios)(AIPGL(i),i=1,mgl)
       DO I=1,MGL
          READ(aipunit,*,IOSTAT=ios) AIPGL(i)
       END DO
       CALL IOERROR(IOS,"Can't read data from air pressure file")

       IF(DBG_SET(DBG_IO)) THEN
          WRITE(IPT,*) "MIN/MAX/MEAN(AIP)::",MINVAL(aipgl),maxval(aipgl),sum(aipgl)/real(Mgl)
       END IF
    END IF

    ! IF NOT PAR, VARIABLES ARE ALREADY POINTED CORRECTLY
    IF (PAR)THEN
# if defined(MULTIPROCESSOR)
       
       SOURCE = MSRID -1

       CALL MPI_BCAST(hour,1,MPI_F,SOURCE,MPI_FVCOM_GROUP,IERR)
       CALL PDEAL(MYID,MSRID,NPROCS,NMAP,AIPGL,AIP%DATA)

       IF(MSR) THEN
          DEALLOCATE(AIPGL)
       END IF
# endif
    END IF

    NULLIFY(AIPGL)

    AIP%dtm = seconds2time(hour*3600.0_SP) + ZEROTIME

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END: READ_AIP"

  END SUBROUTINE READ_AIP

  SUBROUTINE READ_ICE(SAT,SPQ,CLD,RH)
    IMPLICIT NONE
    TYPE(bin_data) :: SAT,SPQ,CLD,RH
    REAL(SP) :: hour,yr,day,iter
    integer :: i, SOURCE, ios
    Real(SP), POINTER :: SATGL(:),SPQGL(:),CLDGL(:),RHGL(:)
    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START: READ_ICE"
    

    IF(MSR) THEN
       IF(PAR) THEN
          ALLOCATE(SATGL(MGL))
          ALLOCATE(SPQGL(MGL))
          ALLOCATE(CLDGL(MGL))
          ALLOCATE(RHGL(MGL))
       ELSE
          SATGL => SAT%DATA(1:MGL)
          SPQGL => SPQ%DATA(1:MGL)
          CLDGL => CLD%DATA(1:MGL)
          RHGL => RH%DATA(1:MGL)
       END IF

! *** EJA Edit 4/29/2011 - adjust for formatted input files
       READ(satunit,*,IOSTAT=ios) hour,yr,day,iter
       CALL IOERROR(IOS,"Can't read hour from air temperature file")
       DO I=1,MGL
          READ(satunit,*,IOSTAT=ios) SATGL(i)
       END DO
!       READ(UNIT=satunit,IOSTAT=ios)(SATGL(i),i=1,mgl)
       CALL IOERROR(IOS,"Can't read data from air temperature file")
       
       READ(spqunit,*,IOSTAT=ios) hour,yr,day,iter
       CALL IOERROR(IOS,"Can't read hour from specific humidity file")
       DO I=1,MGL
          READ(spqunit,*,IOSTAT=ios) SPQGL(i),RHGL(i)
       END DO
!       READ(UNIT=spqunit,IOSTAT=ios)(SPQGL(i),i=1,mgl)
       CALL IOERROR(IOS,"Can't read data from specific humidity file")

       READ(cldunit,*,IOSTAT=ios) hour,yr,day,iter
       CALL IOERROR(IOS,"Can't read hour from cloud cover file")
       DO I=1,MGL
          READ(cldunit,*,IOSTAT=ios) CLDGL(i)
       END DO
!       READ(UNIT=cldunit,IOSTAT=ios)(CLDGL(i),i=1,mgl)
       CALL IOERROR(IOS,"Can't read data from cloud cover file")

       IF(DBG_SET(DBG_IO)) THEN
          WRITE(IPT,*) "MIN/MAX/MEAN(SAT)::",MINVAL(satgl),maxval(satgl),sum(satgl)/real(Mgl)
          WRITE(IPT,*) "MIN/MAX/MEAN(SPQ)::",MINVAL(spqgl),maxval(spqgl),sum(spqgl)/real(Mgl)
          WRITE(IPT,*) "MIN/MAX/MEAN(CLD)::",MINVAL(cldgl),maxval(cldgl),sum(cldgl)/real(Mgl)
          WRITE(IPT,*) "MIN/MAX/MEAN(RH)::",MINVAL(rhgl),maxval(rhgl),sum(rhgl)/real(Mgl)
       END IF


    END IF


    ! IF NOT PAR, VARIABLES ARE ALREADY POINTED CORRECTLY
    IF (PAR)THEN
# if defined(MULTIPROCESSOR)
       
       SOURCE = MSRID -1

       CALL MPI_BCAST(hour,1,MPI_F,SOURCE,MPI_FVCOM_GROUP,IERR)
       CALL PDEAL(MYID,MSRID,NPROCS,NMAP,SATGL,SAT%DATA)
       CALL PDEAL(MYID,MSRID,NPROCS,NMAP,SPQGL,SPQ%DATA)
       CALL PDEAL(MYID,MSRID,NPROCS,NMAP,CLDGL,CLD%DATA)
       CALL PDEAL(MYID,MSRID,NPROCS,NMAP,RHGL,RH%DATA)

       IF(MSR) THEN
          DEALLOCATE(SATGL)
          DEALLOCATE(SPQGL)
          DEALLOCATE(CLDGL)
          DEALLOCATE(RHGL)
       END IF
# endif
    END IF

    NULLIFY(SATGL,SPQGL,CLDGL,RHGL)

    SAT%dtm = seconds2time(hour*3600.0_SP) + ZEROTIME

    SPQ%dtm = SAT%dtm
    CLD%dtm = SAT%dtm
    RH%dtm = SAT%dtm
    
    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END: READ_ICE"

  END SUBROUTINE READ_ICE

! EJA <
  SUBROUTINE READ_SOLAR(SAT,DPT,CLD)
    IMPLICIT NONE
    TYPE(bin_data) :: SAT,DPT,CLD
    REAL(SP) :: hour,yr,day,iter
    integer :: i, SOURCE, ios
    Real(SP), POINTER :: SATGL(:),DPTGL(:),CLDGL(:)
    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START: READ_SOLAR"
    

    IF(MSR) THEN
       IF(PAR) THEN
          ALLOCATE(SATGL(MGL))
          ALLOCATE(DPTGL(MGL))
          ALLOCATE(CLDGL(MGL))
       ELSE
          SATGL => SAT%DATA(1:MGL)
          DPTGL => DPT%DATA(1:MGL)
          CLDGL => CLD%DATA(1:MGL)
       END IF

! *** EJA Edit 4/29/2011 - adjust for formatted input files
       READ(satunit,*,IOSTAT=ios) hour,yr,day,iter
       CALL IOERROR(IOS,"Can't read hour from air temperature file")
       DO I=1,MGL
          READ(satunit,*,IOSTAT=ios) SATGL(i)
       END DO
!       READ(UNIT=satunit,IOSTAT=ios)(SATGL(i),i=1,mgl)
       CALL IOERROR(IOS,"Can't read data from air temperature file")
       
       READ(dptunit,*,IOSTAT=ios) hour,yr,day,iter
       CALL IOERROR(IOS,"Can't read hour from dew point file")
       DO I=1,MGL
          READ(dptunit,*,IOSTAT=ios) DPTGL(i)
       END DO
!       READ(UNIT=spqunit,IOSTAT=ios)(SPQGL(i),i=1,mgl)
       CALL IOERROR(IOS,"Can't read data from dew point file")

       READ(cldunit,*,IOSTAT=ios) hour,yr,day,iter
       CALL IOERROR(IOS,"Can't read hour from cloud cover file")
       DO I=1,MGL
          READ(cldunit,*,IOSTAT=ios) CLDGL(i)
       END DO
!       READ(UNIT=cldunit,IOSTAT=ios)(CLDGL(i),i=1,mgl)
       CALL IOERROR(IOS,"Can't read data from cloud cover file")

       IF(DBG_SET(DBG_IO)) THEN
          WRITE(IPT,*) "MIN/MAX/MEAN(SAT)::",MINVAL(satgl),maxval(satgl),sum(satgl)/real(Mgl)
          WRITE(IPT,*) "MIN/MAX/MEAN(DPT)::",MINVAL(dptgl),maxval(dptgl),sum(dptgl)/real(Mgl)
          WRITE(IPT,*) "MIN/MAX/MEAN(CLD)::",MINVAL(cldgl),maxval(cldgl),sum(cldgl)/real(Mgl)
       END IF


    END IF


    ! IF NOT PAR, VARIABLES ARE ALREADY POINTED CORRECTLY
    IF (PAR)THEN
# if defined(MULTIPROCESSOR)
       
       SOURCE = MSRID -1

       CALL MPI_BCAST(hour,1,MPI_F,SOURCE,MPI_FVCOM_GROUP,IERR)
       CALL PDEAL(MYID,MSRID,NPROCS,NMAP,SATGL,SAT%DATA)
       CALL PDEAL(MYID,MSRID,NPROCS,NMAP,DPTGL,DPT%DATA)
       CALL PDEAL(MYID,MSRID,NPROCS,NMAP,CLDGL,CLD%DATA)

       IF(MSR) THEN
          DEALLOCATE(SATGL)
          DEALLOCATE(DPTGL)
          DEALLOCATE(CLDGL)
       END IF
# endif
    END IF

    NULLIFY(SATGL,DPTGL,CLDGL)

    SAT%dtm = seconds2time(hour*3600.0_SP) + ZEROTIME

    DPT%dtm = SAT%dtm
    CLD%dtm = SAT%dtm
    
    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END: READ_SOLAR"

  END SUBROUTINE READ_SOLAR
! EJA >  
  

  SUBROUTINE UPDATE_WND(NOW,WNDX,WNDY)
    IMPLICIT NONE
    TYPE(TIME) :: NOW
    REAL(SP), POINTER :: WNDX(:), WNDY(:)
    TYPE(BIN_DATA), POINTER :: A, B
    REAL(DP)     :: denom, numer
    REAL(SP)     :: fw, bw

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START:UPDATE_WND "

    DO       
       IF(NOW .LT. WNDX_PREV%dtm) THEN

          CALL PRINT_REAL_TIME(NOW,IPT,"OUTPUT TIME",timezone)
          CALL PRINT_REAL_TIME(WNDX_PREV%dtm,IPT,"DATA TIME",timezone)

          CALL FATAL_ERROR("CAN NOT REWIND BINARY FILES",&
               & "SOMETHING IS WRONG WITH TIME IN THE WIND FILE")
          
       ELSE IF(NOW .gt. WNDX_NEXT%dtm)THEN
          
          A=> WNDX_PREV
          WNDX_PREV => WNDX_NEXT
          
          B => WNDY_PREV
          WNDY_PREV => WNDY_NEXT
          
          CALL READ_WND(WNDX=A,WNDY=B)

          WNDX_NEXT => A
          WNDY_NEXT => B
          
       ELSE
          EXIT
          
       END IF
       
    END DO

    NUMER = SECONDS(NOW - WNDX_PREV%dtm)

    DENOM = SECONDS(WNDY_NEXT%dtm - WNDX_PREV%dtm)
    
    fw = NUMER/DENOM
    bw = 1.0_DP - NUMER/DENOM

    WNDX= WNDX_NEXT%data *fw + WNDX_PREV%data *bw 

    WNDY= WNDY_NEXT%data *fw + WNDY_PREV%data *bw 

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END:UPDATE_WND "

  END SUBROUTINE UPDATE_WND

  SUBROUTINE UPDATE_HFX(NOW,SWR,NHF,LWR)
    IMPLICIT NONE
    TYPE(TIME) :: NOW
    REAL(SP), POINTER :: SWR(:), NHF(:), LWR(:)
    TYPE(BIN_DATA), POINTER :: A, B, C
    REAL(DP)     :: denom, numer
    REAL(SP)     :: fw, bw

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START:UPDATE_HFX "

    DO       
       IF(NOW .LT. SWR_PREV%dtm) THEN

          CALL PRINT_REAL_TIME(NOW,IPT,"OUTPUT TIME",timezone)
          CALL PRINT_REAL_TIME(SWR_PREV%dtm,IPT,"DATA TIME",timezone)

          CALL FATAL_ERROR("CAN NOT REWIND BINARY FILES",&
               & "SOMETHING IS WRONG WITH TIME IN THE HEATING FILE")
          
       ELSE IF(NOW .gt. SWR_NEXT%dtm)THEN
          
          A=> SWR_PREV
          SWR_PREV => SWR_NEXT
          
          B => NHF_PREV
          NHF_PREV => NHF_NEXT
          
          C => LWR_PREV
          LWR_PREV => LWR_NEXT
          
          CALL READ_HFX(SWR=A,NHF=B,LWR=C)
          SWR_NEXT => A
          NHF_NEXT => B
          LWR_NEXT => C
!!zhuxm 2013-05-15  NHF_NEXT => C
          
       ELSE
          
          EXIT
          
       END IF
       
    END DO
       
    NUMER = SECONDS(NOW - SWR_PREV%dtm)

    DENOM = SECONDS(SWR_NEXT%dtm - SWR_PREV%dtm)
    
    fw = NUMER/DENOM
    bw = 1.0_DP - NUMER/DENOM

    SWR= SWR_NEXT%data *fw + SWR_PREV%data *bw 

    NHF= NHF_NEXT%data *fw + NHF_PREV%data *bw 

    LWR= LWR_NEXT%data *fw + LWR_PREV%data *bw 


    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END:UPDATE_HFX "

  END SUBROUTINE UPDATE_HFX

  SUBROUTINE UPDATE_EVP(NOW,EVP,PRC)
    IMPLICIT NONE
    TYPE(TIME) :: NOW
    REAL(SP), POINTER :: EVP(:), PRC(:)
    TYPE(BIN_DATA), POINTER :: A, B
    REAL(DP)     :: denom, numer
    REAL(SP)     :: fw, bw
    
    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START:UPDATE_EVP "

    DO       
       IF(NOW .LT. EVP_PREV%dtm) THEN

          CALL PRINT_REAL_TIME(NOW,IPT,"OUTPUT TIME",timezone)
          CALL PRINT_REAL_TIME(EVP_PREV%dtm,IPT,"DATA TIME",timezone)


          CALL FATAL_ERROR("CAN NOT REWIND BINARY FILES",&
               & "SOMETHING IS WRONG WITH TIME IN THE EVAP/PREC FILE")
          
       ELSE IF(NOW .gt. EVP_NEXT%dtm)THEN
          
          A=> EVP_PREV
          EVP_PREV => EVP_NEXT
          
          B => PRC_PREV
          PRC_PREV => PRC_NEXT
          
          CALL READ_EVP(EVP=A,PRC=B)
          EVP_NEXT => A
          PRC_NEXT => B
          
       ELSE
          
          EXIT
          
       END IF
       
    END DO
       
    NUMER = SECONDS(NOW - EVP_PREV%dtm)

    DENOM = SECONDS(PRC_NEXT%dtm - EVP_PREV%dtm)
    
    fw = NUMER/DENOM
    bw = 1.0_DP - NUMER/DENOM

    EVP= EVP_NEXT%data *fw + EVP_PREV%data *bw 

    PRC= PRC_NEXT%data *fw + PRC_PREV%data *bw 

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END:UPDATE_EVP "

  END SUBROUTINE UPDATE_EVP


  SUBROUTINE UPDATE_AIP(NOW,AIP)
    IMPLICIT NONE
    TYPE(TIME) :: NOW
    REAL(SP), POINTER :: AIP(:)
    TYPE(BIN_DATA), POINTER :: A
    REAL(DP)     :: denom, numer
    REAL(SP)     :: fw, bw
    
    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START:UPDATE_AIP "

    DO       
       IF(NOW .LT. AIP_PREV%dtm) THEN

          CALL PRINT_REAL_TIME(NOW,IPT,"OUTPUT TIME",timezone)
          CALL PRINT_REAL_TIME(AIP_PREV%dtm,IPT,"DATA TIME",timezone)


          CALL FATAL_ERROR("CAN NOT REWIND BINARY FILES",&
               & "SOMETHING IS WRONG WITH TIME IN THE AIR PRESSURE FILE")
          
       ELSE IF(NOW .gt. AIP_NEXT%dtm)THEN
          
          A=> AIP_PREV
          AIP_PREV => AIP_NEXT
          
          CALL READ_AIP(AIP=A)
          AIP_NEXT => A
          
       ELSE
          
          EXIT
          
       END IF
       
    END DO
       
    NUMER = SECONDS(NOW - AIP_PREV%dtm)

    DENOM = SECONDS(AIP_NEXT%dtm - AIP_PREV%dtm)
    
    fw = NUMER/DENOM
    bw = 1.0_DP - NUMER/DENOM

    AIP= AIP_NEXT%data *fw + AIP_PREV%data *bw 

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END:UPDATE_AIP "

  END SUBROUTINE UPDATE_AIP


  SUBROUTINE UPDATE_ICE(NOW,SAT,SPQ,CLD,RH)
    IMPLICIT NONE
    TYPE(TIME) :: NOW
    REAL(SP), POINTER :: SAT(:),SPQ(:),CLD(:),RH(:)
    TYPE(BIN_DATA), POINTER :: A,B,C,D
    REAL(DP)     :: denom, numer
    REAL(SP)     :: fw, bw
    
    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START:UPDATE_ICE "

    DO       
       IF(NOW .LT. SAT_PREV%dtm) THEN

          CALL PRINT_REAL_TIME(NOW,IPT,"OUTPUT TIME",timezone)
          CALL PRINT_REAL_TIME(SAT_PREV%dtm,IPT,"DATA TIME",timezone)


          CALL FATAL_ERROR("CAN NOT REWIND BINARY FILES",&
               & "SOMETHING IS WRONG WITH TIME IN THE ICE MODEL FILES")
          
       ELSE IF(NOW .gt. SAT_NEXT%dtm)THEN
          
          A=> SAT_PREV
          SAT_PREV => SAT_NEXT
          B=> SPQ_PREV
          SPQ_PREV => SPQ_NEXT
          C=> CLD_PREV
          CLD_PREV => CLD_NEXT
          D=> RH_PREV
          RH_PREV => RH_NEXT
          
          CALL READ_ICE(SAT=A,SPQ=B,CLD=C,RH=D)
          SAT_NEXT => A
          SPQ_NEXT => B
          CLD_NEXT => C
          RH_NEXT => D
          
       ELSE
          
          EXIT
          
       END IF
       
    END DO
       
    NUMER = SECONDS(NOW - SAT_PREV%dtm)

    DENOM = SECONDS(SAT_NEXT%dtm - SAT_PREV%dtm)
    
    fw = NUMER/DENOM
    bw = 1.0_DP - NUMER/DENOM

    SAT= SAT_NEXT%data *fw + SAT_PREV%data *bw 
    SPQ= SPQ_NEXT%data *fw + SPQ_PREV%data *bw 
    CLD= CLD_NEXT%data *fw + CLD_PREV%data *bw 
    RH= RH_NEXT%data *fw + RH_PREV%data *bw 

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END:UPDATE_ICE "

  END SUBROUTINE UPDATE_ICE

! EJA <
  SUBROUTINE UPDATE_SOLAR(NOW,SAT,DPT,CLD)
    IMPLICIT NONE
    TYPE(TIME) :: NOW
    REAL(SP), POINTER :: SAT(:),DPT(:),CLD(:)
    TYPE(BIN_DATA), POINTER :: A,B,C,D
    REAL(DP)     :: denom, numer
    REAL(SP)     :: fw, bw
    
    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START:UPDATE_SOLAR "

    DO       
       IF(NOW .LT. SAT_PREV%dtm) THEN

          CALL PRINT_REAL_TIME(NOW,IPT,"OUTPUT TIME",timezone)
          CALL PRINT_REAL_TIME(SAT_PREV%dtm,IPT,"DATA TIME",timezone)


          CALL FATAL_ERROR("CAN NOT REWIND BINARY FILES",&
               & "SOMETHING IS WRONG WITH TIME IN THE ICE MODEL FILES")
          
       ELSE IF(NOW .gt. SAT_NEXT%dtm)THEN
          
          A=> SAT_PREV
          SAT_PREV => SAT_NEXT
          B=> DPT_PREV
          DPT_PREV => DPT_NEXT
          C=> CLD_PREV
          CLD_PREV => CLD_NEXT
          
          CALL READ_SOLAR(SAT=A,DPT=B,CLD=C)
          SAT_NEXT => A
          DPT_NEXT => B
          CLD_NEXT => C
          
       ELSE
          
          EXIT
          
       END IF
       
    END DO
       
    NUMER = SECONDS(NOW - SAT_PREV%dtm)

    DENOM = SECONDS(SAT_NEXT%dtm - SAT_PREV%dtm)
    
    fw = NUMER/DENOM
    bw = 1.0_DP - NUMER/DENOM

    SAT= SAT_NEXT%data *fw + SAT_PREV%data *bw 
    DPT= DPT_NEXT%data *fw + DPT_PREV%data *bw 
    CLD= CLD_NEXT%data *fw + CLD_PREV%data *bw 

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END:UPDATE_SOLAR "

  END SUBROUTINE UPDATE_SOLAR
! EJA >






END MODULE MOD_BINARY
