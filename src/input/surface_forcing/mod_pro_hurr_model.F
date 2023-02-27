MODULE MOD_PRO_HURR_MODEL
  USE MOD_BULK
  USE MOD_SURFACEFORCE
  USE MOD_SET_TIME

  IMPLICIT NONE
  
  TYPE bin_data
     type(TIME) :: dtm
     REAL(SP) :: X_E, Y_E, P_E, V_E, R_E
  END TYPE bin_data
  
  TYPE(bin_data), POINTER :: HURR_NEXT, HURR_PREV

  CONTAINS

  SUBROUTINE UPDATE_PRO_HURR_MODEL(NOW)
    IMPLICIT NONE
    TYPE(TIME) :: NOW
    
    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START: UPDATE_PRO_HURR_MODEL"

    IF(WIND_ON) THEN

       IF(WIND_TYPE == 'stress') THEN
         CALL FATAL_ERROR("TO USE PROTOTYPICAL HURRICANE MODEL,",&
	      & "WIND_TYPE SHOULD BE speed.")
       ELSE IF(WIND_TYPE == 'speed') THEN
          CALL UPDATE_HURR(NOW,SPDX,SPDY,AIP)
       END IF
       
    END IF

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END: UPDATE_PRO_HURR_MODEL"
    
  END SUBROUTINE UPDATE_PRO_HURR_MODEL


  FUNCTION NEW_DATA(dims)
    IMPLICIT NONE
    TYPE(bin_data), POINTER :: NEW_DATA
    integer, intent(IN) :: DIMS
    INTEGER :: STATUS

    ALLOCATE(NEW_DATA,stat=status)
    IF(status /=0) CALL FATAL_ERROR("NEW_DATA: COULD NOT ALLOCATE TYPE POINTER?")
    
  END FUNCTION NEW_DATA


  SUBROUTINE LOAD_PRO_HURR_MODEL(WND)
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: WND
    INTEGER :: STATUS
    INTEGER(ITIME) :: dummy
    CHARACTER(LEN=4) :: FLAG
    integer :: ios

    
    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START:LOAD_PRO_HURR_MODEL"
    
    IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) "Looking for hurricane path file:"&
         &//TRIM(WND)
    inquire(file=trim(WND),exist=WIND_ON)
    IF(WIND_ON) THEN
       
       IF(WIND_TYPE/='speed')THEN
          CALL FATAL_ERROR("To get hurricane wind and pressure file,",&
               &"you must specify 'binary speed'")
       END IF
       
       IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) "FOUND HURRICANE PATH FILE: OPEN AND READ"
       
       IF(MSR) CALL FOPEN(WNDUNIT,WND,'cfs')
       
       HURR_NEXT => NEW_DATA(N)
       HURR_PREV => NEW_DATA(N)
       
       READ(wndunit,'(a)',IOSTAT=ios)
       CALL IOERROR(IOS,"error on reading the first line of hurricane path file")
       READ(wndunit,'(a)',IOSTAT=ios)
       CALL IOERROR(IOS,"error on reading the second line of hurricane path file")

       CALL READ_HURR(HURR=HURR_PREV)
       CALL READ_HURR(HURR=HURR_NEXT)
       
       IF(DBG_SET(DBG_LOG)) THEN
          CALL PRINT_REAL_TIME(HURR_PREV%dtm,IPT,"FIRST TIME POINT",timezone)
          CALL PRINT_REAL_TIME(HURR_NEXT%dtm,IPT,"SECOND TIME POINT",timezone)
       END IF
       
       IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) "FOUND HURRICANE PATH FILE: READ FIRST DATA POINTS"
       
    ELSE
       
       IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) "! NO HURRICANE PATH FILE FOUND"
       
    END IF
    
    IF (.not. PRECIPITATION_ON .and. .not. HEATING_ON .and. &
      & .not. WIND_ON .and. .not. AIRPRESSURE_ON ) &
      &  CALL FATAL_ERROR("FOUND NO ASCII HURRICANE PSTH INPUT FILES?")
    
    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END:LOAD_PRO_HURR_MODEL"
    
  END SUBROUTINE LOAD_PRO_HURR_MODEL
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
  SUBROUTINE READ_HURR(HURR)
    IMPLICIT NONE
    TYPE(BIN_DATA) :: HURR
    REAL(SP) :: hour,X_E,Y_E,P_E,V_E,R_E
    integer :: i, SOURCE, ios
    INTEGER :: IYEAR,IMONTH,IDAY,IHOUR
    CHARACTER(LEN=4) :: CYEAR,CHOUR
    CHARACTER(LEN=2) :: CMONTH,CDAY
    CHARACTER(LEN=25) :: TS
    TYPE(TIME) :: GET_TIME
    integer :: status
   
    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START: READ_HURR"


    IF(MSR) THEN
 
       READ(wndunit,*,IOSTAT=ios) IYEAR,IMONTH,IDAY,IHOUR,Y_E,X_E,P_E,V_E,R_E
       CALL IOERROR(IOS,"Can't read data from hurricane path file")
       
       WRITE(CYEAR,  '(I4.4)') IYEAR
       WRITE(CMONTH, '(I2.2)') IMONTH
       WRITE(CDAY,   '(I2.2)') IDAY
       WRITE(CHOUR,  '(I4.4)') IHOUR
       
       TS = CYEAR//"/"//CMONTH//"/"//CDAY//"/"//" "//CHOUR(1:2)//":"//CHOUR(3:4)//":"//"00"
       
       GET_TIME = READ_DATETIME(TRIM(TS),'ymd','UTC',status)
       CALL PRINT_TIME(GET_TIME,IPT,TS)
       
       R_E = R_E*1.852_SP * 1000.0_SP      ! mile -> m
       V_E = V_E*0.51444444_SP             ! kt -> m/s

    END IF

    HURR%dtm = GET_TIME      
    HURR%X_E = X_E
    HURR%Y_E = Y_E
    HURR%P_E = P_E
    HURR%V_E = V_E
    HURR%R_E = R_E

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END: READ_WND"
  END SUBROUTINE READ_HURR

!-------------------------------------------------------------------
  SUBROUTINE UPDATE_HURR(NOW,WNDX,WNDY,AIP)
    IMPLICIT NONE
    TYPE(TIME) :: NOW
    REAL(SP), POINTER :: WNDX(:), WNDY(:), AIP(:)
    TYPE(BIN_DATA), POINTER :: A, B
    REAL(DP)     :: denom, numer
    REAL(SP)     :: fw, bw
    REAL(SP) :: X_CENTER, Y_CENTER, P_C, V_MAX, R_MAX
    REAL(SP) :: X_TEMP, Y_TEMP, RRR,BBB,COR
    REAL(SP), ALLOCATABLE, TARGET :: WIND_S(:), WIND_D(:),AIP_T(:)
    REAL(SP), ALLOCATABLE :: LONC_H(:), LON_H(:)
    INTEGER  :: I
    
    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START:UPDATE_HURR "

    DO       
      IF(NOW < HURR_PREV%dtm) THEN

        CALL PRINT_REAL_TIME(NOW,IPT,"OUTPUT TIME",timezone)
        CALL PRINT_REAL_TIME(HURR_PREV%dtm,IPT,"DATA TIME",timezone)

        CALL FATAL_ERROR("CAN NOT REWIND INPUT FILES",&
             & "SOMETHING IS WRONG WITH TIME IN THE INPUT FILE")
          
      ELSE IF(NOW > HURR_NEXT%dtm)THEN
          
        A=> HURR_PREV
        HURR_PREV => HURR_NEXT
          
        CALL READ_HURR(HURR=A)

        HURR_NEXT => A
          
      ELSE
          
        EXIT
          
      END IF
       
    END DO

    NUMER = SECONDS(NOW - HURR_PREV%dtm)

    DENOM = SECONDS(HURR_NEXT%dtm - HURR_PREV%dtm)
    
    fw = NUMER/DENOM
    bw = 1.0_DP - NUMER/DENOM

    ALLOCATE(LONC_H(0:NT))
    ALLOCATE(LON_H(0:MT))
    LONC_H = LONC
    WHERE(LONC_H < 0.0_SP) LONC_H = LONC_H + 360.0_SP
    LON_H  = LON
    WHERE(LON_H < 0.0_SP) LON_H = LON_H + 360.0_SP

    X_CENTER = HURR_NEXT%X_E*fw + HURR_PREV%X_E*bw 
    Y_CENTER = HURR_NEXT%Y_E*fw + HURR_PREV%Y_E*bw 
    P_C      = HURR_NEXT%P_E*fw + HURR_PREV%P_E*bw 
    V_MAX    = HURR_NEXT%V_E*fw + HURR_PREV%V_E*bw 
    R_MAX    = HURR_NEXT%R_E*fw + HURR_PREV%R_E*bw 

    IF(X_CENTER < 0.0_SP) X_CENTER = X_CENTER + 360.0_SP
    BBB=1.15_SP*EXP(1.0_SP)*V_MAX**2/(SLP0-P_C*100.0_SP)     !RHO_A = 1.15_SP
    IF(BBB < 1.0_SP) BBB=1.0_SP
    IF(BBB >= 2.5_SP) BBB=2.5_SP

    ALLOCATE(WIND_S(0:NT))         ;WIND_S = 0.0_SP
    ALLOCATE(WIND_D(0:NT))         ;WIND_D = 0.0_SP
#   if defined (AIR_PRESSURE)
    ALLOCATE(AIP_T(0:MT))          ;AIP_T    = 0.0_SP
#   endif
    
#   if defined (SPHERICAL)
    DO I = 1,N
      X_TEMP = (LONC_H(I)-X_CENTER)*DEG2RAD
      Y_TEMP = (LATC(I)-Y_CENTER)*DEG2RAD
      RRR = REARTH*(2.0_SP*ASIN(SQRT(SIN(0.5_SP*Y_TEMP)**2+  &
            COS(Y_CENTER*DEG2RAD)*COS(LATC(I)*DEG2RAD)*   &
	    SIN(0.5_SP*X_TEMP)**2)))
      IF(RRR < 1.0_SP) RRR = 1.0_SP
      
      COR = 2.*7.292e-5_SP * SIN(LATC(I)*DEG2RAD)
      
      WIND_S(I) = SQRT((R_MAX/RRR)**BBB * &
               EXP(1.0_SP-(R_MAX/RRR)**BBB)*V_MAX**2 + &
	       RRR*RRR*COR*COR/4)-RRR*COR/2
      WIND_D(I) = ATAN2(Y_TEMP,X_TEMP)*180.0_SP/3.1415926_SP+90.0_SP
      IF(WIND_D(I) < 0.0_SP) WIND_D(I) = WIND_D(I)+360.0_SP
    END DO

    WNDX = WIND_S*COS(WIND_D*DEG2RAD)
    WNDY = WIND_S*SIN(WIND_D*DEG2RAD)
    
#   if defined (AIR_PRESSURE)
    DO I=1,M
      X_TEMP = (LON_H(I)-X_CENTER)*DEG2RAD
      Y_TEMP = (LAT(I)-Y_CENTER)*DEG2RAD
      RRR = REARTH*(2.0_SP*ASIN(SQRT(SIN(0.5_SP*Y_TEMP)**2+ &
            COS(Y_CENTER*DEG2RAD)*COS(LAT(I)*DEG2RAD)* &
	    SIN(0.5_SP*X_TEMP)**2)))
      IF(RRR < 1.0_SP) RRR = 1.0_SP
      AIP_T(I) = P_C*100._SP+(SLP0-P_C*100._SP)*EXP(-(R_MAX/RRR)**BBB)
    END DO
    AIP = AIP_T
#   endif
#   else
    IF(USE_PROJ)THEN
      DO I = 1,N
        X_TEMP = (LONC_H(I)-X_CENTER)*DEG2RAD
        Y_TEMP = (LATC(I)-Y_CENTER)*DEG2RAD
        RRR = REARTH*(2.0_SP*ASIN(SQRT(SIN(0.5_SP*Y_TEMP)**2+  &
              COS(Y_CENTER*DEG2RAD)*COS(LATC(I)*DEG2RAD)*   &
	      SIN(0.5_SP*X_TEMP)**2)))
        IF(RRR < 1.0_SP) RRR = 1.0_SP

        COR = 2.*7.292e-5_SP * SIN(LATC(I)*DEG2RAD)
      
        WIND_S(I) = SQRT((R_MAX/RRR)**BBB * &
                 EXP(1.0_SP-(R_MAX/RRR)**BBB)*V_MAX**2 + &
	         RRR*RRR*COR*COR/4)-RRR*COR/2
        WIND_D(I) = ATAN2(Y_TEMP,X_TEMP)*180.0_SP/3.1415926_SP+90.0_SP
        IF(WIND_D(I) < 0.0_SP) WIND_D(I) = WIND_D(I)+360.0_SP
      END DO
      WNDX = WIND_S*COS(WIND_D*DEG2RAD)
      WNDY = WIND_S*SIN(WIND_D*DEG2RAD)
    
#     if defined (AIR_PRESSURE)
      DO I=1,M
        X_TEMP = (LON_H(I)-X_CENTER)*DEG2RAD
        Y_TEMP = (LAT(I)-Y_CENTER)*DEG2RAD
        RRR = REARTH*(2.0_SP*ASIN(SQRT(SIN(0.5_SP*Y_TEMP)**2+ &
              COS(Y_CENTER*DEG2RAD)*COS(LAT(I)*DEG2RAD)* &
	      SIN(0.5_SP*X_TEMP)**2)))
        IF(RRR < 1.0_SP) RRR = 1.0_SP
        AIP_T(I) = P_C*100._SP+(SLP0-P_C*100._SP)*EXP(-(R_MAX/RRR)**BBB)
      END DO
      AIP = AIP_T
#     endif
    ELSE
      CALL FATAL_ERROR("USE_PROJ must be true when using cartesian coordinates")
    END IF
#   endif      

    DEALLOCATE(WIND_S,WIND_D)
#   if defined (AIR_PRESSURE)
    DEALLOCATE(AIP_T)
#   endif

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END:UPDATE_HURR "

  END SUBROUTINE UPDATE_HURR

END MODULE MOD_PRO_HURR_MODEL
