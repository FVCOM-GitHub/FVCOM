MODULE MOD_MM5 
  USE MOD_BULK
  USE MOD_SURFACEFORCE
 implicit none
  

 ! HARD CODED PATH TO MEDM MM5 RESULTS GRID
 CHARACTER(LEN=100),parameter :: mm5ll = "/hosts/medm-vault/data01/shu/data_source/raw/util/mm5_domain2.ll" 

 INTEGER, parameter :: mx=129,ny=102

 TYPE(INTERP_WEIGHTS) :: INTP_N
 TYPE(INTERP_WEIGHTS) :: INTP_C

 REAL(SP), ALLOCATABLE, DIMENSION(:,:) :: mm5lon,mm5lat,mm5x,mm5y



 TYPE mm5_data
    type(TIME) :: dtm
    REAL(SP)   :: data(mx,ny)
 END TYPE mm5_data
  
  TYPE(mm5_data), POINTER :: SWR_NEXT, SWR_PREV ! SHORT WAVE RADIATION
  TYPE(mm5_data), POINTER :: NHF_NEXT, NHF_PREV ! NET HEAT FLUX
  TYPE(mm5_data), POINTER :: LNG_NEXT, LNG_PREV ! LONG WAVE RADIATION
  TYPE(mm5_data), POINTER :: SEN_NEXT, SEN_PREV ! SENSIBLE HEAT FLUX
  TYPE(mm5_data), POINTER :: LAT_NEXT, LAT_PREV ! LATENT HEAT FLUX


  TYPE(mm5_data), POINTER :: WNDX_NEXT, WNDX_PREV ! WIND SPEED, EASTWARD
  TYPE(mm5_data), POINTER :: WNDY_NEXT, WNDY_PREV ! WIND SPEED, NORTHWARD
  TYPE(mm5_data), POINTER :: STRX_NEXT, STRX_PREV ! WIND STRESS, EASTWARD
  TYPE(mm5_data), POINTER :: STRY_NEXT, STRY_PREV ! WIND STRESS, NORTHWARD


  TYPE(mm5_data), POINTER :: EVP_NEXT, EVP_PREV ! EVAPORATION
  TYPE(mm5_data), POINTER :: PRC_NEXT, PRC_PREV ! PRECIPITATION
  TYPE(mm5_data), POINTER :: PRS_NEXT, PRS_PREV ! PRESSURE
  TYPE(mm5_data), POINTER :: SST_NEXT, SST_PREV ! SEA SURFACE TEMPERATURE
  TYPE(mm5_data), POINTER :: SAT_NEXT, SAT_PREV ! AIR SURFACE TEMPERATURE


CONTAINS


    SUBROUTINE UPDATE_MM5(NOW)
    IMPLICIT NONE
    TYPE(TIME) :: NOW
    
    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START: UPDATE_BINARY"

    IF(PRECIPITATION_ON) THEN
       IF(OUTPUT_TYPE == struct) THEN
          CALL UPDATE_EVP_RG(NOW,EVP_RG,PRC_RG,PRS_RG,SST_RG,SAT_RG)
       ELSE
          CALL UPDATE_EVP_US(NOW,EVP,PRC,PRS,SST,SAT)
       END IF

    END IF

    IF(WIND_ON) THEN

       IF(OUTPUT_TYPE == struct) THEN
          CALL UPDATE_WND_RG(NOW,SPDX_RG,SPDY_RG,STRX_RG,STRY_RG)
       ELSE
          CALL UPDATE_WND_US(NOW,SPDX,SPDY,STRX,STRY)
       END IF
    END IF


    IF(HEATING_ON) THEN
       IF(OUTPUT_TYPE == struct) THEN
          CALL UPDATE_HFX_RG(NOW,SWR_RG,NHF_RG,LNGWV_RG,SENSE_RG,LATENT_RG)
       ELSE
          CALL UPDATE_HFX_US(NOW,SWR,NHF,LNGWV,SENSE,LATENT)
       END IF
    END IF


    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END: UPDATE_BINARY"
    
  END SUBROUTINE UPDATE_MM5
  



  SUBROUTINE LOAD_MM5(WND,HFX,EVP)
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: WND
    CHARACTER(LEN=*), INTENT(IN) :: HFX
    CHARACTER(LEN=*), INTENT(IN) :: EVP
    INTEGER :: STATUS, I,J, SOURCE
    INTEGER(ITIME) :: dummy
    CHARACTER(LEN=4) :: FLAG
    LOGICAL :: EXST

    real(SP), allocatable :: XTMP(:),YTMP(:),lontmp(:),lattmp(:)


      IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START:LOAD_MM5"
    
      
      allocate(mm5x(mx,ny))
      allocate(mm5y(mx,ny))


      allocate(mm5lon(mx,ny))
      allocate(mm5lat(mx,ny))

      IF(MSR) THEN
         
         inquire(file=trim(mm5ll),exist=EXST )
         IF(.not.EXST) CALL FATAL_ERROR&
              &("PLEASE SPECIFY THE LOCATION OF THE MM5 LAT/LON data file in mod_mm5.F and recompile!",&
              & "The file:"//TRIM(mm5ll)//"; does not exist!")
         
         
         allocate(XTMP(mx*ny))
         allocate(YTMP(mx*ny))

         allocate(lonTMP(mx*ny))
         allocate(latTMP(mx*ny))
         

         CALL FOPEN(mm5unit,mm5ll,'cfr')
         do i=1,mx*ny
            read(mm5unit,*)lontmp(i),lattmp(i)
         enddo

         CALL DEGREES2METERS(lontmp,lattmp,PROJECTION_REFERENCE,xtmp,ytmp,mx*ny)        

         DO J = 1,ny
            DO I = 1,mx

               mm5x(I,J) = xtmp(I + (J-1)*mx)
               mm5y(I,J) = ytmp(I + (J-1)*mx)

               mm5lon(I,J) = lontmp(I + (J-1)*mx)
               mm5lat(I,J) = lattmp(I + (J-1)*mx)
            END DO
         END DO

         DEALLOCATE(xtmp,ytmp,lontmp,lattmp)

         ! TO DUMP THE XY OR LAT LON LOCATION DATA
!         CALL FOPEN(84,"XXXX",'ofr')
!         DO J = 1,ny
!            WRITE(UNIT= 84,FMT='(102F12.0)') mm5lon(:,j)
!         END DO
!         CLOSE(84)
!         CALL FOPEN(85,"YYYY",'ofr')
!         DO J = 1,ny
!            WRITE(UNIT= 85,FMT='(102F12.0)') mm5lat(:,j)
!         END DO
!         CLOSE(85)


     END IF


     IF (PAR)THEN
# if defined(MULTIPROCESSOR)

      IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "SENDING COORDS DATA"

        SOURCE = MSRID -1
        CALL MPI_BCAST(mm5x,mx*ny,MPI_F,SOURCE,MPI_FVCOM_GROUP)
        CALL MPI_BCAST(mm5y,mx*ny,MPI_F,SOURCE,MPI_FVCOM_GROUP)

        CALL MPI_BCAST(mm5lon,mx*ny,MPI_F,SOURCE,MPI_FVCOM_GROUP)
        CALL MPI_BCAST(mm5lat,mx*ny,MPI_F,SOURCE,MPI_FVCOM_GROUP)
     
# endif
     END IF

     IF(OUTPUT_TYPE /= struct)THEN
        IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "CREATING INTERP MATRIX"
        
        CALL SETUP_INTERP_BILINEAR_A(MM5X,MM5Y,XM,YM,INTP_N)
        
        CALL SETUP_INTERP_BILINEAR_A(MM5X,MM5Y,XMC,YMC,INTP_C)
     END IF
     

     IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "LOADING FILES..."

      IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) "Looking for Wind Stress file:"&
           &//TRIM(WND)
      inquire(file=trim(WND),exist=WIND_ON)
      IF(WIND_ON) THEN
       
       IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) "FOUND WIND FILE: OPEN AND READ"
       
       IF(MSR) CALL FOPEN(WNDUNIT,WND,'cur')
       
       WNDX_NEXT => NEW_DATA()
       WNDX_PREV => NEW_DATA()
       WNDY_NEXT => NEW_DATA()
       WNDY_PREV => NEW_DATA()

       STRX_NEXT => NEW_DATA()
       STRX_PREV => NEW_DATA()
       STRY_NEXT => NEW_DATA()
       STRY_PREV => NEW_DATA()
       
       CALL READ_WND(WNDX=WNDX_PREV,WNDY=WNDY_PREV,STRX=STRX_PREV,STRY=STRY_PREV)
       CALL READ_WND(WNDX=WNDX_NEXT,WNDY=WNDY_NEXT,STRX=STRX_NEXT,STRY=STRY_NEXT)
       
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
       
       IF(MSR) CALL FOPEN(HFXUNIT,HFX,'cur')
       
       SWR_NEXT => NEW_DATA()
       SWR_PREV => NEW_DATA()

       NHF_NEXT => NEW_DATA()
       NHF_PREV => NEW_DATA()

       LNG_NEXT => NEW_DATA()
       LNG_PREV => NEW_DATA()

       SEN_NEXT => NEW_DATA()
       SEN_PREV => NEW_DATA()

       LAT_NEXT => NEW_DATA()
       LAT_PREV => NEW_DATA()

       CALL READ_HFX(SWR=SWR_PREV,NHF=NHF_PREV,LNG=LNG_PREV,SEN=SEN_PREV,LAT=LAT_PREV)

       CALL READ_HFX(SWR=SWR_NEXT,NHF=NHF_NEXT,LNG=LNG_NEXT,SEN=SEN_NEXT,LAT=LAT_NEXT )

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

       EVP_NEXT => NEW_DATA()
       EVP_PREV => NEW_DATA()

       PRC_NEXT => NEW_DATA()
       PRC_PREV => NEW_DATA()

       PRS_NEXT => NEW_DATA()
       PRS_PREV => NEW_DATA()
       
       SST_NEXT => NEW_DATA()
       SST_PREV => NEW_DATA()

       SAT_NEXT => NEW_DATA()
       SAT_PREV => NEW_DATA()




       CALL READ_EVP(EVP=EVP_PREV,PRC=PRC_PREV,PRS=PRS_PREV,SST=SST_PREV,SAT=SAT_PREV)
       CALL READ_EVP(EVP=EVP_NEXT,PRC=PRC_NEXT,PRS=PRS_NEXT,SST=SST_NEXT,SAT=SAT_NEXT)
       
       IF(DBG_SET(DBG_LOG)) THEN
          CALL PRINT_REAL_TIME(SWR_PREV%dtm,IPT,"FIRST TIME POINT",timezone)
          CALL PRINT_REAL_TIME(SWR_NEXT%dtm,IPT,"SECOND TIME POINT",timezone)
       END IF
       
       IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) "FOUND PRECIPITATION FILE: READ FIRST DATA POINTS"
       
    ELSE
       
       IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) "! NO PRECIPITATION FILE FOUND"
       
    END IF

    IF (.not. PRECIPITATION_ON .and. .not. HEATING_ON .and. &
      & .not. WIND_ON ) &
      &  CALL FATAL_ERROR("FOUND NO MM5 FORCING INPUT FILES?")
    
    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END:LOAD_MM5"


  END SUBROUTINE LOAD_MM5
  
  FUNCTION NEW_DATA()
    IMPLICIT NONE
    TYPE(mm5_data), POINTER :: NEW_DATA
!    integer, intent(IN) :: DIMS
    INTEGER :: STATUS

    ALLOCATE(NEW_DATA,stat=status)
    IF(status /=0) CALL FATAL_ERROR("NEW_DATA: COULD NOT ALLOCATE TYPE POINTER?")
    
    NEW_DATA%DTM = ZEROTIME

 !   ALLOCATE(NEW_DATA%DATA(0:DIMS), STAT=STATUS)
 !   IF(status /=0) CALL FATAL_ERROR("NEW_DATA: COULD NOT ALLOCATE DATA POINTER?")
    
  END FUNCTION NEW_DATA

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
  SUBROUTINE READ_WND(WNDX,WNDY,STRX,STRY)
    IMPLICIT NONE
    TYPE(mm5_data) :: WNDX,WNDY,STRX,STRY
    REAL(SP) :: hour
    integer :: i,j, SOURCE, ios
   
    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START: READ_WND"


    IF(MSR) THEN
 

       READ(UNIT=wndunit,IOSTAT=ios) hour
       CALL IOERROR(IOS,"Can't read hour from wind file")
       READ(UNIT= wndunit,IOSTAT=ios)&
            &((WNDX%DATA(i,j),WNDY%DATA(i,j),STRX%DATA(i,j),STRY%DATA(i,j),i=1,mx),j=1,ny)
       CALL IOERROR(IOS,"Can't read data from wind file")

       IF(DBG_SET(DBG_IO)) THEN
          WRITE(IPT,*) "MIN/MAX/MEAN(WNDX)::",MINVAL(WNDX%DATA),maxval(WNDX%DATA),sum(WNDX%DATA)/real(mx*ny)
          WRITE(IPT,*) "MIN/MAX/MEAN(STRX)::",MINVAL(STRX%DATA),maxval(STRX%DATA),sum(STRX%DATA)/real(mx*ny)

          WRITE(IPT,*) "MIN/MAX/MEAN(WNDY)::",MINVAL(WNDY%DATA),maxval(WNDY%DATA),sum(WNDY%DATA)/real(mx*ny)
          WRITE(IPT,*) "MIN/MAX/MEAN(STRY)::",MINVAL(STRY%DATA),maxval(STRY%DATA),sum(STRY%DATA)/real(mx*ny)

       END IF

    END IF


    ! IF NOT PAR, VARIABLES ARE ALREADY POINTED CORRECTLY
    IF (PAR)THEN
# if defined(MULTIPROCESSOR)
       
       SOURCE = MSRID -1

       CALL MPI_BCAST(hour,1,MPI_F,SOURCE,MPI_FVCOM_GROUP)
       CALL MPI_BCAST(WNDX%DATA,mx*ny,MPI_F,SOURCE,MPI_FVCOM_GROUP)
       CALL MPI_BCAST(WNDY%DATA,mx*ny,MPI_F,SOURCE,MPI_FVCOM_GROUP)
       CALL MPI_BCAST(STRX%DATA,mx*ny,MPI_F,SOURCE,MPI_FVCOM_GROUP)
       CALL MPI_BCAST(STRY%DATA,mx*ny,MPI_F,SOURCE,MPI_FVCOM_GROUP)

# endif
    END IF


    WNDX%dtm = seconds2time(hour*3600.0_SP) + ZEROTIME

    WNDY%dtm = WNDX%dtm

    STRY%dtm = WNDX%dtm

    STRX%dtm = WNDX%dtm

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END: READ_WND"
  END SUBROUTINE READ_WND

  !-------------------------------------------------------------------
  SUBROUTINE READ_HFX(SWR,NHF,LNG,SEN,LAT)
    IMPLICIT NONE
    TYPE(mm5_data) :: SWR,NHF,LNG,SEN,LAT
    REAL(SP) :: hour
    integer :: i,j, SOURCE, ios
   
    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START: READ_HFX"

    IF(MSR) THEN
 

       READ(UNIT=HFXunit,IOSTAT=ios) hour
       CALL IOERROR(IOS,"Can't read hour from heat file")
       READ(UNIT= HFXunit,IOSTAT=ios)&
            &((NHF%DATA(i,j),SWR%DATA(i,j),LNG%DATA(i,j),SEN%DATA(i,j),LAT%DATA(i,j),i=1,mx),j=1,ny)
       CALL IOERROR(IOS,"Can't read data from wind file")

       IF(DBG_SET(DBG_IO)) THEN
          WRITE(IPT,*) "MIN/MAX/MEAN(NHF)::",MINVAL(NHF%DATA),maxval(NHF%DATA),sum(NHF%DATA)/real(mx*ny)
          WRITE(IPT,*) "MIN/MAX/MEAN(SWR)::",MINVAL(SWR%DATA),maxval(SWR%DATA),sum(SWR%DATA)/real(mx*ny)

          WRITE(IPT,*) "MIN/MAX/MEAN(LNG)::",MINVAL(LNG%DATA),maxval(LNG%DATA),sum(LNG%DATA)/real(mx*ny)
          WRITE(IPT,*) "MIN/MAX/MEAN(SEN)::",MINVAL(SEN%DATA),maxval(SEN%DATA),sum(SEN%DATA)/real(mx*ny)
          WRITE(IPT,*) "MIN/MAX/MEAN(LAT)::",MINVAL(LAT%DATA),maxval(LAT%DATA),sum(LAT%DATA)/real(mx*ny)

       END IF

    END IF


    ! IF NOT PAR, VARIABLES ARE ALREADY POINTED CORRECTLY
    IF (PAR)THEN
# if defined(MULTIPROCESSOR)
       
       SOURCE = MSRID -1

       CALL MPI_BCAST(hour,1,MPI_F,SOURCE,MPI_FVCOM_GROUP)
       CALL MPI_BCAST(NHF%DATA,mx*ny,MPI_F,SOURCE,MPI_FVCOM_GROUP)
       CALL MPI_BCAST(SWR%DATA,mx*ny,MPI_F,SOURCE,MPI_FVCOM_GROUP)
       CALL MPI_BCAST(LNG%DATA,mx*ny,MPI_F,SOURCE,MPI_FVCOM_GROUP)
       CALL MPI_BCAST(SEN%DATA,mx*ny,MPI_F,SOURCE,MPI_FVCOM_GROUP)
       CALL MPI_BCAST(LAT%DATA,mx*ny,MPI_F,SOURCE,MPI_FVCOM_GROUP)

# endif
    END IF


    NHF%dtm = seconds2time(hour*3600.0_SP) + ZEROTIME

    SWR%dtm = NHF%dtm

    LNG%dtm = NHF%dtm

    SEN%dtm = NHF%dtm

    LAT%dtm = NHF%dtm

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END: READ_HFX"
  END SUBROUTINE READ_HFX

  !-------------------------------------------------------------------
  SUBROUTINE READ_EVP(EVP,PRC,PRS,SST,SAT)
    IMPLICIT NONE
    TYPE(mm5_data) :: EVP,PRC,PRS,SST,SAT
    REAL(SP) :: hour
    integer :: i,j, SOURCE, ios
   
    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START: READ_EVP"

    IF(MSR) THEN
 

       READ(UNIT=EVPunit,IOSTAT=ios) hour
       CALL IOERROR(IOS,"Can't read hour from heat file")
       READ(UNIT= EVPunit,IOSTAT=ios)&
            &((PRS%DATA(i,j),SST%DATA(i,j),SAT%DATA(i,j),PRC%DATA(i,j),EVP%DATA(i,j),i=1,mx),j=1,ny)
       CALL IOERROR(IOS,"Can't read data from wind file")

       ! CONVERT UNITS
       PRC%DATA = PRC%DATA/ real(100 * 3600)
       EVP%DATA = EVP%DATA/ real(100 * 3600)

       IF(DBG_SET(DBG_IO)) THEN
          WRITE(IPT,*) "MIN/MAX/MEAN(PRC)::",MINVAL(PRC%DATA),maxval(PRC%DATA),sum(PRC%DATA)/real(mx*ny)
          WRITE(IPT,*) "MIN/MAX/MEAN(EVP)::",MINVAL(EVP%DATA),maxval(EVP%DATA),sum(EVP%DATA)/real(mx*ny)

          WRITE(IPT,*) "MIN/MAX/MEAN(PRS)::",MINVAL(PRS%DATA),maxval(PRS%DATA),sum(PRS%DATA)/real(mx*ny)
          WRITE(IPT,*) "MIN/MAX/MEAN(SST)::",MINVAL(SST%DATA),maxval(SST%DATA),sum(SST%DATA)/real(mx*ny)
          WRITE(IPT,*) "MIN/MAX/MEAN(SAT)::",MINVAL(SAT%DATA),maxval(SAT%DATA),sum(SAT%DATA)/real(mx*ny)

       END IF

    END IF


    ! IF NOT PAR, VARIABLES ARE ALREADY POINTED CORRECTLY
    IF (PAR)THEN
# if defined(MULTIPROCESSOR)
       
       SOURCE = MSRID -1

       CALL MPI_BCAST(hour,1,MPI_F,SOURCE,MPI_FVCOM_GROUP)
       CALL MPI_BCAST(PRC%DATA,mx*ny,MPI_F,SOURCE,MPI_FVCOM_GROUP)
       CALL MPI_BCAST(EVP%DATA,mx*ny,MPI_F,SOURCE,MPI_FVCOM_GROUP)
       CALL MPI_BCAST(PRS%DATA,mx*ny,MPI_F,SOURCE,MPI_FVCOM_GROUP)
       CALL MPI_BCAST(SST%DATA,mx*ny,MPI_F,SOURCE,MPI_FVCOM_GROUP)
       CALL MPI_BCAST(SAT%DATA,mx*ny,MPI_F,SOURCE,MPI_FVCOM_GROUP)

# endif
    END IF


    write(ipt,*) "READ HOUR:", HOUR

    PRC%dtm = seconds2time(hour*3600.0_SP) + ZEROTIME

    EVP%dtm = PRC%dtm

    PRS%dtm = PRC%dtm

    SST%dtm = PRC%dtm

    SAT%dtm = PRC%dtm

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END: READ_EVP"
  END SUBROUTINE READ_EVP


  SUBROUTINE UPDATE_WND_US(NOW,SPDX,SPDY,STRX,STRY)
    IMPLICIT NONE
    TYPE(TIME) :: NOW
    REAL(SP), POINTER :: STRX(:), STRY(:),SPDX(:), SPDY(:)
    TYPE(MM5_DATA), POINTER :: A, B,C, D
    REAL(DP)     :: denom, numer
    REAL(SP)     :: fw, bw
   REAL(SP), POINTER :: TEMP(:,:)

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START:UPDATE_WND "

    DO       
       IF(NOW .LT. STRX_PREV%dtm) THEN

          CALL PRINT_REAL_TIME(NOW,IPT,"OUTPUT TIME",timezone)
          CALL PRINT_REAL_TIME(STRX_PREV%dtm,IPT,"DATA TIME",timezone)

          CALL FATAL_ERROR("CAN NOT REWIND MM5 FILES",&
               & "SOMETHING IS WRONG WITH TIME IN THE WIND FILE")
          
       ELSE IF(NOW .gt. STRX_NEXT%dtm)THEN
          
          A=> WNDX_PREV
          WNDX_PREV => WNDX_NEXT
          
          B => WNDY_PREV
          WNDY_PREV => WNDY_NEXT

          C=> STRX_PREV
          STRX_PREV => STRX_NEXT
          
          D => STRY_PREV
          STRY_PREV => STRY_NEXT
          
          CALL READ_WND(WNDX=A,WNDY=B,STRX=C,STRY=D)

          WNDX_NEXT => A
          WNDY_NEXT => B

          STRX_NEXT => C
          STRY_NEXT => D
          
       ELSE
          
          EXIT
          
       END IF
       
    END DO

    NUMER = SECONDS(NOW - WNDX_PREV%dtm)

    DENOM = SECONDS(WNDX_NEXT%dtm - WNDX_PREV%dtm)
    
    fw = NUMER/DENOM
    bw = 1.0_DP - NUMER/DENOM

    ALLOCATE(TEMP(mx,ny))

    TEMP= WNDX_NEXT%data *fw + WNDX_PREV%data *bw 
    CALL INTERP_BILINEAR_P(TEMP,INTP_C,SPDX)

    TEMP= WNDY_NEXT%data *fw + WNDY_PREV%data *bw 
    CALL INTERP_BILINEAR_P(TEMP,INTP_C,SPDY)

    TEMP= STRX_NEXT%data *fw + STRX_PREV%data *bw 
    CALL INTERP_BILINEAR_P(TEMP,INTP_C,STRX)

    TEMP= STRY_NEXT%data *fw + STRY_PREV%data *bw 
    CALL INTERP_BILINEAR_P(TEMP,INTP_C,STRY)

    DEALLOCATE(TEMP)
 
    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END:UPDATE_WND "

  END SUBROUTINE UPDATE_WND_US

  SUBROUTINE UPDATE_WND_RG(NOW,SPDX,SPDY,STRX,STRY)
    IMPLICIT NONE
    TYPE(TIME) :: NOW
    REAL(SP), POINTER :: STRX(:,:), STRY(:,:),SPDX(:,:), SPDY(:,:)
    TYPE(MM5_DATA), POINTER :: A, B,C, D
    REAL(DP)     :: denom, numer
    REAL(SP)     :: fw, bw

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START:UPDATE_WND "

    DO       
       IF(NOW .LT. STRX_PREV%dtm) THEN

          CALL PRINT_REAL_TIME(NOW,IPT,"OUTPUT TIME",timezone)
          CALL PRINT_REAL_TIME(STRX_PREV%dtm,IPT,"DATA TIME",timezone)

          CALL FATAL_ERROR("CAN NOT REWIND MM5 FILES",&
               & "SOMETHING IS WRONG WITH TIME IN THE WIND FILE")
          
       ELSE IF(NOW .gt. STRX_NEXT%dtm)THEN
          
          A=> WNDX_PREV
          WNDX_PREV => WNDX_NEXT
          
          B => WNDY_PREV
          WNDY_PREV => WNDY_NEXT

          C=> STRX_PREV
          STRX_PREV => STRX_NEXT
          
          D => STRY_PREV
          STRY_PREV => STRY_NEXT
          
          CALL READ_WND(WNDX=A,WNDY=B,STRX=C,STRY=D)

          WNDX_NEXT => A
          WNDY_NEXT => B

          STRX_NEXT => C
          STRY_NEXT => D
          
       ELSE
          
          EXIT
          
       END IF
       
    END DO

    NUMER = SECONDS(NOW - WNDX_PREV%dtm)

    DENOM = SECONDS(WNDX_NEXT%dtm - WNDX_PREV%dtm)
    
    fw = NUMER/DENOM
    bw = 1.0_DP - NUMER/DENOM

    SPDX= WNDX_NEXT%data *fw + WNDX_PREV%data *bw 

    SPDY= WNDY_NEXT%data *fw + WNDY_PREV%data *bw 

    STRX= STRX_NEXT%data *fw + STRX_PREV%data *bw 

    STRY= STRY_NEXT%data *fw + STRY_PREV%data *bw 

 
    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END:UPDATE_WND "

  END SUBROUTINE UPDATE_WND_RG


  SUBROUTINE UPDATE_EVP_US(NOW,EVP,PRC,PRS,SST,SAT)
    IMPLICIT NONE
    TYPE(TIME) :: NOW
    REAL(SP), POINTER :: EVP(:), PRC(:),PRS(:), SST(:), SAT(:)
    TYPE(MM5_DATA), POINTER :: A, B,C, D,E
    REAL(DP)     :: denom, numer
    REAL(SP)     :: fw, bw
    REAL(SP), POINTER :: TEMP(:,:)

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START:UPDATE_EVP "

    DO       
       IF(NOW .LT. EVP_PREV%dtm) THEN

          CALL PRINT_REAL_TIME(NOW,IPT,"OUTPUT TIME",timezone)
          CALL PRINT_REAL_TIME(EVP_PREV%dtm,IPT,"DATA TIME",timezone)

          CALL FATAL_ERROR("CAN NOT REWIND MM5 FILES",&
               & "SOMETHING IS WRONG WITH TIME IN THE PREC/EVAP FILE")
          
       ELSE IF(NOW .gt. EVP_NEXT%dtm)THEN
          
          A=> PRS_PREV
          PRS_PREV => PRS_NEXT
          
          B => SST_PREV
          SST_PREV => SST_NEXT

          C=> EVP_PREV
          EVP_PREV => EVP_NEXT
          
          D => PRC_PREV
          PRC_PREV => PRC_NEXT

          E => SAT_PREV
          SAT_PREV => SAT_NEXT
          
          CALL READ_EVP(EVP=C,PRC=D,PRS=A,SST=B,SAT=E)

          PRS_NEXT => A
          SST_NEXT => B

          EVP_NEXT => C
          PRC_NEXT => D
          SAT_NEXT => E
          
       ELSE
          
          EXIT
          
       END IF
       
    END DO

    NUMER = SECONDS(NOW - PRS_PREV%dtm)

    DENOM = SECONDS(SST_NEXT%dtm - PRS_PREV%dtm)
    
    fw = NUMER/DENOM
    bw = 1.0_DP - NUMER/DENOM

    ALLOCATE(TEMP(mx,ny))

    TEMP= PRS_NEXT%data *fw + PRS_PREV%data *bw 
    CALL INTERP_BILINEAR_P(TEMP,INTP_N,PRS)

    TEMP= SST_NEXT%data *fw + SST_PREV%data *bw 
    CALL INTERP_BILINEAR_P(TEMP,INTP_N,SST)

    TEMP= SAT_NEXT%data *fw + SAT_PREV%data *bw 
    CALL INTERP_BILINEAR_P(TEMP,INTP_N,SAT)

    TEMP= EVP_NEXT%data *fw + EVP_PREV%data *bw 
    CALL INTERP_BILINEAR_P(TEMP,INTP_N,EVP)

    TEMP= PRC_NEXT%data *fw + PRC_PREV%data *bw 
    CALL INTERP_BILINEAR_P(TEMP,INTP_N,PRC)

    DEALLOCATE(TEMP)

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END:UPDATE_EVP "

  END SUBROUTINE UPDATE_EVP_US

  SUBROUTINE UPDATE_EVP_RG(NOW,EVP,PRC,PRS,SST,SAT)
    IMPLICIT NONE
    TYPE(TIME) :: NOW
    REAL(SP), POINTER :: EVP(:,:), PRC(:,:),PRS(:,:), SST(:,:), SAT(:,:)
    TYPE(MM5_DATA), POINTER :: A, B,C, D,E
    REAL(DP)     :: denom, numer
    REAL(SP)     :: fw, bw

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START:UPDATE_EVP "

    DO       
       IF(NOW .LT. EVP_PREV%dtm) THEN

          CALL PRINT_REAL_TIME(NOW,IPT,"OUTPUT TIME",timezone)
          CALL PRINT_REAL_TIME(EVP_PREV%dtm,IPT,"DATA TIME",timezone)

          CALL FATAL_ERROR("CAN NOT REWIND MM5 FILES",&
               & "SOMETHING IS WRONG WITH TIME IN THE PREC/EVAP FILE")
          
       ELSE IF(NOW .gt. EVP_NEXT%dtm)THEN
          
          A=> PRS_PREV
          PRS_PREV => PRS_NEXT
          
          B => SST_PREV
          SST_PREV => SST_NEXT

          C=> EVP_PREV
          EVP_PREV => EVP_NEXT
          
          D => PRC_PREV
          PRC_PREV => PRC_NEXT

          E => SAT_PREV
          SAT_PREV => SAT_NEXT
          
          CALL READ_EVP(EVP=C,PRC=D,PRS=A,SST=B,SAT=E)

          PRS_NEXT => A
          SST_NEXT => B

          EVP_NEXT => C
          PRC_NEXT => D
          PRC_NEXT => E
          
       ELSE
          
          EXIT
          
       END IF
       
    END DO

    NUMER = SECONDS(NOW - PRS_PREV%dtm)

    DENOM = SECONDS(SST_NEXT%dtm - PRS_PREV%dtm)
    
    fw = NUMER/DENOM
    bw = 1.0_DP - NUMER/DENOM


    PRS= PRS_NEXT%data *fw + PRS_PREV%data *bw 

    SST= SST_NEXT%data *fw + SST_PREV%data *bw 

    SAT= SAT_NEXT%data *fw + SAT_PREV%data *bw 
    
    EVP= EVP_NEXT%data *fw + EVP_PREV%data *bw 

    PRC= PRC_NEXT%data *fw + PRC_PREV%data *bw 

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END:UPDATE_EVP "

  END SUBROUTINE UPDATE_EVP_RG


  SUBROUTINE UPDATE_HFX_US(NOW,SWR,NHF,LNG,SEN,LAT)
    IMPLICIT NONE
    TYPE(TIME) :: NOW
    REAL(SP), POINTER :: SWR(:), NHF(:),LNG(:), SEN(:), LAT(:)
    TYPE(MM5_DATA), POINTER :: A, B,C, D,E
    REAL(DP)     :: denom, numer
    REAL(SP)     :: fw, bw

    REAL(SP), POINTER :: TEMP(:,:)

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START:UPDATE_HFX "

    DO       
       IF(NOW .LT. SWR_PREV%dtm) THEN

          CALL PRINT_REAL_TIME(NOW,IPT,"OUTPUT TIME",timezone)
          CALL PRINT_REAL_TIME(SWR_PREV%dtm,IPT,"DATA TIME",timezone)

          CALL FATAL_ERROR("CAN NOT REWIND MM5 FILES",&
               & "SOMETHING IS WRONG WITH TIME IN THE PREC/EVAP FILE")
          
       ELSE IF(NOW .gt. SWR_NEXT%dtm)THEN
          
          A=> LNG_PREV
          LNG_PREV => LNG_NEXT
          
          B => SEN_PREV
          SEN_PREV => SEN_NEXT

          C=> SWR_PREV
          SWR_PREV => SWR_NEXT
          
          D => NHF_PREV
          NHF_PREV => NHF_NEXT

          E => LAT_PREV
          LAT_PREV => LAT_NEXT
          
          CALL READ_HFX(SWR=C,NHF=D,LNG=A,SEN=B,LAT=E)

          LNG_NEXT => A
          SEN_NEXT => B

          SWR_NEXT => C
          NHF_NEXT => D
          LAT_NEXT => E
          
       ELSE
          
          EXIT
          
       END IF
       
    END DO

    NUMER = SECONDS(NOW - LNG_PREV%dtm)

    DENOM = SECONDS(SEN_NEXT%dtm - LNG_PREV%dtm)
    
    fw = NUMER/DENOM
    bw = 1.0_DP - NUMER/DENOM

    ALLOCATE(TEMP(mx,ny))

    TEMP= LNG_NEXT%data *fw + LNG_PREV%data *bw 
    CALL INTERP_BILINEAR_P(TEMP,INTP_N,LNG)

    TEMP= SEN_NEXT%data *fw + SEN_PREV%data *bw 
    CALL INTERP_BILINEAR_P(TEMP,INTP_N,SEN)

    TEMP= LAT_NEXT%data *fw + LAT_PREV%data *bw 
    CALL INTERP_BILINEAR_P(TEMP,INTP_N,LAT)

    TEMP= SWR_NEXT%data *fw + SWR_PREV%data *bw 
    CALL INTERP_BILINEAR_P(TEMP,INTP_N,SWR)    

    TEMP= NHF_NEXT%data *fw + NHF_PREV%data *bw 
    CALL INTERP_BILINEAR_P(TEMP,INTP_N,NHF)

    DEALLOCATE(TEMP)

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END:UPDATE_HFX "

  END SUBROUTINE UPDATE_HFX_US

  SUBROUTINE UPDATE_HFX_RG(NOW,SWR,NHF,LNG,SEN,LAT)
    IMPLICIT NONE
    TYPE(TIME) :: NOW
    REAL(SP), POINTER :: SWR(:,:), NHF(:,:),LNG(:,:), SEN(:,:), LAT(:,:)
    TYPE(MM5_DATA), POINTER :: A, B,C, D,E
    REAL(DP)     :: denom, numer
    REAL(SP)     :: fw, bw

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START:UPDATE_HFX "

    DO       
       IF(NOW .LT. SWR_PREV%dtm) THEN

          CALL PRINT_REAL_TIME(NOW,IPT,"OUTPUT TIME",timezone)
          CALL PRINT_REAL_TIME(SWR_PREV%dtm,IPT,"DATA TIME",timezone)

          CALL FATAL_ERROR("CAN NOT REWIND MM5 FILES",&
               & "SOMETHING IS WRONG WITH TIME IN THE PREC/EVAP FILE")
          
       ELSE IF(NOW .gt. SWR_NEXT%dtm)THEN
          
          A=> LNG_PREV
          LNG_PREV => LNG_NEXT
          
          B => SEN_PREV
          SEN_PREV => SEN_NEXT

          C=> SWR_PREV
          SWR_PREV => SWR_NEXT
          
          D => NHF_PREV
          NHF_PREV => NHF_NEXT

          E => LAT_PREV
          LAT_PREV => LAT_NEXT
          
          CALL READ_HFX(SWR=C,NHF=D,LNG=A,SEN=B,LAT=E)

          LNG_NEXT => A
          SEN_NEXT => B

          SWR_NEXT => C
          NHF_NEXT => D
          LAT_NEXT => E
          
       ELSE
          
          EXIT
          
       END IF
       
    END DO

    NUMER = SECONDS(NOW - LNG_PREV%dtm)

    DENOM = SECONDS(SEN_NEXT%dtm - LNG_PREV%dtm)
    
    fw = NUMER/DENOM
    bw = 1.0_DP - NUMER/DENOM

    LNG= LNG_NEXT%data *fw + LNG_PREV%data *bw 

    SEN= SEN_NEXT%data *fw + SEN_PREV%data *bw 
    
    LAT= LAT_NEXT%data *fw + LAT_PREV%data *bw 
    
    SWR= SWR_NEXT%data *fw + SWR_PREV%data *bw 
    
    NHF= NHF_NEXT%data *fw + NHF_PREV%data *bw 
    

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END:UPDATE_HFX "

  END SUBROUTINE UPDATE_HFX_RG





END MODULE MOD_MM5
