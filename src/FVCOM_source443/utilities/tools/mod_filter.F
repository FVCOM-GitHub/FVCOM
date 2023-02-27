Module Mod_Filter
  USE MOD_PREC
  USE MOD_UTILS
  USE CONTROL, only : PI
  implicit none

  ! FUNCTIONS:
  ! PL66_SIZE
  ! PL66_FREQ_CUT
  !
  ! SUBROUTINES:
  ! PL66_WTS
  ! COS_TAPER
  
  contains
  
    INTEGER FUNCTION PL66_SIZE(HAP,DT)
      USE MOD_TIME
      IMPLICIT NONE
      TYPE(TIME),  INTENT(IN)  :: dt  ! Time series interval
      TYPE(TIME),  INTENT(IN)  :: HAP ! Half Amplitude Period

      INTEGER nHalf

      nhalf=ANINT(2.0_DP * SECONDS(HAP) / SECONDS(dt))-1 ! Round to nearest and drop the last
      PL66_SIZE=2*nhalf+1
      
    END FUNCTION PL66_SIZE

    REAL(DP) FUNCTION PL66_FREQ_CUT(HAP,DT)
      USE MOD_TIME
      IMPLICIT NONE
      TYPE(TIME),  INTENT(IN)  :: dt  ! Time series interval
      TYPE(TIME),  INTENT(IN)  :: HAP ! Half Amplitude Period
      REAL(DP) CUTOFF

      CUTOFF = SECONDS(HAP)/SECONDS(DT)
      PL66_FREQ_CUT=1.0_DP/cutoff
      
    END FUNCTION PL66_FREQ_CUT

  SUBROUTINE PL66_WTS(WTS,NSIZE,FREQ_CUT)
    USE CONTROL 

    IMPLICIT NONE
    INTEGER,  INTENT(IN)  :: NSIZE      ! FILTER SIZE
    REAL(SP), INTENT(IN)  :: FREQ_CUT   ! CUT OFF FREQUENCY
    REAL(SP), INTENT(OUT) :: WTS(NSIZE) ! FILTER WTS

    REAL(SP)  :: coef,phi
    REAL(DP)  :: NORM
    INTEGER   :: I,nhalf
    !==============================================================================|
    
    IF(MOD(NSIZE,2)/=0 .OR. NSIZE<=0) CALL FATAL_ERROR("PL66_WTS: invalid size!")

    ! GET THE HALF SIZE
    nhalf = (NSIZE-1)/2
    
    
    DO I=1,nhalf
       phi=PI*I
       coef=(freq_cut**2)*(phi**3)
       WTS(1+NSIZE+I)=(2.0*SIN(2.0*freq_cut*phi)-SIN(freq_cut*phi)-SIN(3.0*freq_cut*phi))/coef
    END DO
    
    
    DO I=1,nhalf
       WTS(I)=WTS(NSIZE-I+1)
    END DO
    
    WTS(nhalf+1)=2*freq_cut
    
    norm=sum(REAL(WTS,DP))
    
    
    WTS = WTS/NORM
    
  END SUBROUTINE PL66_WTS
  !------------------------------------------------------------------------------!
  ! Cosine taper is the extent of the taper at the ends of the data
  ! cosine_taper=0.5 corrisponds to about a 70% cosine filter 
  !------------------------------------------------------------------------------!
  SUBROUTINE COS_TAPER(CS,NSIZE,C_TAP)
    IMPLICIT NONE
    INTEGER,  INTENT(IN)  :: NSIZE      ! FILTER SIZE
    REAL(SP), INTENT(IN)  :: C_TAP      ! TAPER COEFFICIENT
    REAL(SP), INTENT(OUT) :: CS(NSIZE)  ! COSINE TAPER

    integer :: nhalf, I

    IF(MOD(NSIZE,2)/=0 .OR. NSIZE<=0) CALL FATAL_ERROR("COS_TAPER: invalid size!")
    
    ! GET THE HALF SIZE
    nhalf = (NSIZE-1)/2
    
    DO I=1, nhalf
       CS(I)=COS(c_tap*(PI*I)/NSIZE)
    END DO

  END SUBROUTINE COS_TAPER

End Module Mod_Filter
