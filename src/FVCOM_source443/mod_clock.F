!/===========================================================================/
! Copyright (c) 2007, The University of Massachusetts Dartmouth 
! Produced at the School of Marine Science & Technology 
! Marine Ecosystem Dynamics Modeling group
! All rights reserved.
!
! FVCOM has been developed by the joint UMASSD-WHOI research team. For 
! details of authorship and attribution of credit please see the FVCOM
! technical manual or contact the MEDM group.
!
! 
! This file is part of FVCOM. For details, see http://fvcom.smast.umassd.edu 
! The full copyright notice is contained in the file COPYRIGHT located in the 
! root directory of the FVCOM code. This original header must be maintained
! in all distributed versions.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
! AND ANY EXPRESS OR  IMPLIED WARRANTIES, INCLUDING,  BUT NOT  LIMITED TO,
! THE IMPLIED WARRANTIES OF MERCHANTABILITY AND  FITNESS FOR A PARTICULAR
! PURPOSE ARE DISCLAIMED.  
!
!/---------------------------------------------------------------------------/
! CVS VERSION INFORMATION
! $Id$
! $Name$
! $Revision$
!/===========================================================================/

MODULE MOD_CLOCK
   USE MOD_PREC
   USE MOD_TIME
   USE CONTROL, only: USE_REAL_WORLD_TIME,IPT
   USE MOD_UTILS
   IMPLICIT NONE

! SUBROUTINES IN MOD_CLOCK:
!   TIMER ROUTINES
! WATCH_INIT(TYPE::WATCH)
! WATCH_TIME(TYPE::WATCH) RETURN(REAL(DP))
! WATCH_RESET(TYPE::WATCH)
! WATCH_LAP(TYPE::WATCH)
! WATCH_LAP(TYPE::WATCH,INTEGER::UNIT#,CHARACTER(LEN*)::MESSAGE)
! WATCH_REPORT(TYPE::WATCH,INTEGER,CHARACTER(LEN*))
!
!    TIME STRING FORMATTING
! GETTIME(CHARACTER(LEN=13),INTEGER)
! GET_TIMESTAMP(CHARACTER(LEN=*))
!
! REPORT_TIME(INTEGER(ITIME),INTEGER(ITIME),INTEGER(ITIME),TYPE::TIME)


!Michael Dunphy (Michael.Dunphy@dfo-mpo.gc.ca)
!For longer simulations, the SECS/IT printout in the log gets a nonsense value
! due to 32-bit int overflow; the patch changes the relevant ints to 64-bits
!   TYPE WATCH
!      INTEGER  :: COUNT_RATE,COUNT_MAX
!      INTEGER  :: COUNT_0
!      INTEGER  :: COUNT_CURRENT, COUNT_LAST
!      INTEGER  :: Lap_Count
!   END TYPE WATCH

   TYPE WATCH
      INTEGER(ITIME)  :: COUNT_RATE,COUNT_MAX
      INTEGER(ITIME)  :: COUNT_0
      INTEGER(ITIME)  :: COUNT_CURRENT, COUNT_LAST
      INTEGER(ITIME)  :: Lap_Count
   END TYPE WATCH
!Michael Dunphy


   TYPE(WATCH):: Report_watch

   INTERFACE WATCH_LAP
      MODULE PROCEDURE WATCH_LAP_NOPRINT
      MODULE PROCEDURE WATCH_LAP_PRINT
   END INTERFACE

!!$EXAMPLE OF TIMER USE:
!!$
!!$! Declare a watch type
!!$TYPE(WATCH):: atimer
!!$real(sp) :: atime
!!$
!!$! Initialize the Timer
!!$CALL WATCH_INIT(atimer)
!!$
!!$
!!$Do J = 1:KB
!!$   
!!$   Do I = 1:NT
!!$   ! DO SOMETHING
!!$
!!$   END Do
!!$   ! Prints the time take for each internal do-loop
!!$   call Watch_lap(atimer)   ! JUST COUNT FOR LATER
!!$        OR
!!$   call Watch_lap(atimer,ipt,"what am I counting!") ! PRINTS TO SCREEN
!!$END Do
!!$
!!$ ! Prints the average time taken for all laps
!!$ call Watch_report(atimer)
!!$
!!$ ! Return time since watch started
!!$ atime = Watch_time(atimer)
!!$
!!$ ! RESET THE WATCH
!!$ CALL WATCH_RESET(atimer)


   CONTAINS

!==============================================================================!
!  Initialize STOPWATCH                                                        !
!==============================================================================!
SUBROUTINE WATCH_INIT(MYWATCH)
  IMPLICIT NONE
  TYPE(WATCH) :: MYWATCH
  
  MYWATCH%COUNT_RATE    = 0
  MYWATCH%COUNT_MAX     = 0
  MYWATCH%COUNT_CURRENT = 0
  MYWATCH%COUNT_LAST    = 0
  MYWATCH%COUNT_0        = 0
  MYWATCH%Lap_Count     = 0

 CALL SYSTEM_CLOCK (MYWATCH%COUNT_CURRENT, MYWATCH%COUNT_RATE, MYWATCH%COUNT_MAX)
 
 MYWATCH%COUNT_0     = MYWATCH%COUNT_CURRENT
 MYWATCH%COUNT_LAST = MYWATCH%COUNT_CURRENT

     IF(MYWATCH%COUNT_MAX == 0) CALL FATAL_ERROR &
          & ("MOD_CLOCK: FORTRAN DOES NOT RECOGNIZE YOUR SYSTEM CLOCK!",&
          & "INTIRINISIC, SYSTEM_CLOCK returned COUNT_MAX==0")

     IF(MYWATCH%COUNT_RATE .LE. 0) CALL FATAL_ERROR &
          & ("MOD_CLOCK: FORTRAN DOES NOT RECOGNIZE YOUR SYSTEM CLOCK!",&
          & "INTIRINISIC, SYSTEM_CLOCK returned COUNT_RATE .LE. 0")

     IF(MYWATCH%COUNT_0 .LT. 0) CALL FATAL_ERROR &
          & ("MOD_CLOCK: FORTRAN DOES NOT RECOGNIZE YOUR SYSTEM CLOCK!",&
          & "INTIRINISIC, SYSTEM_CLOCK returned COUNT .LT. 0")

END SUBROUTINE WATCH_INIT
!==============================================================================!
!  Retrieve Stopwatch time                                                     !
!  Return a double # of seconds since watch_initialize
!==============================================================================!
   FUNCTION WATCH_TIME(MYWATCH) 
     IMPLICIT NONE  
     TYPE(WATCH) :: MYWATCH
     REAL(DP) :: WATCH_TIME
     INTEGER(ITIME) :: BIG_COUNT
     INTEGER :: MY_COUNT

     CALL SYSTEM_CLOCK(COUNT=MY_COUNT)
     BIG_COUNT = MY_COUNT - MYWATCH%COUNT_0
     IF (BIG_COUNT < 0) BIG_COUNT = BIG_COUNT + MYWATCH%COUNT_MAX

     WATCH_TIME = DBLE(BIG_COUNT) / DBLE(MYWATCH%COUNT_RATE)
     
   END FUNCTION WATCH_TIME

!==============================================================================!
!  Start Watch lapstart                                                             !
!==============================================================================!
   SUBROUTINE WATCH_RESET(MYWATCH) 
     IMPLICIT NONE  
     TYPE(WATCH) :: MYWATCH
     
     CALL SYSTEM_CLOCK(MYWATCH%COUNT_CURRENT)
     MYWATCH%LAP_Count = 0
     MYWATCH%Count_LAST = MYWATCH%COUNT_CURRENT
     MYWATCH%Count_0     = MYWATCH%COUNT_CURRENT

   END SUBROUTINE WATCH_RESET

!==============================================================================!
!  Watch lap                                                             !
!==============================================================================!
   SUBROUTINE WATCH_LAP_NOPRINT(MYWATCH) 
     IMPLICIT NONE  
     TYPE(WATCH) :: MYWATCH

     MYWATCH%COUNT_LAST = MYWATCH%COUNT_CURRENT
     CALL SYSTEM_CLOCK(MYWATCH%COUNT_CURRENT)
     MYWATCH%LAP_Count = MYWATCH%LAP_Count + 1

   END SUBROUTINE WATCH_LAP_NOPRINT

!==============================================================================!
!  Watch Start lap                                                             !
!==============================================================================!
   SUBROUTINE WATCH_START_LAP(MYWATCH) 
     IMPLICIT NONE  
     TYPE(WATCH) :: MYWATCH

     CALL SYSTEM_CLOCK(MYWATCH%COUNT_CURRENT)
     !IGNORE TIME SINCE LAST STOP
     MYWATCH%Count_0 =MYWATCH%Count_0 + (MYWATCH%COUNT_CURRENT - MYWATCH%COUNT_LAST)

   END SUBROUTINE WATCH_START_LAP

!==============================================================================!
!  Watch Stop lap                                                             !
!==============================================================================!
   SUBROUTINE WATCH_STOP_LAP(MYWATCH) 
     IMPLICIT NONE  
     TYPE(WATCH) :: MYWATCH

     CALL SYSTEM_CLOCK(MYWATCH%COUNT_CURRENT)
     MYWATCH%COUNT_LAST = MYWATCH%COUNT_CURRENT
     MYWATCH%LAP_Count = MYWATCH%LAP_Count + 1

   END SUBROUTINE WATCH_STOP_LAP
!==============================================================================!
!  Print Watch Report                                                          !
!==============================================================================!
   SUBROUTINE WATCH_REPORT(MYWATCH,UNIT,MSG) 
     IMPLICIT NONE  
     

     TYPE(WATCH) :: MYWATCH   
     integer, intent(in) :: unit
     Character(len=*), intent(in) :: MSG
     INTEGER(ITIME) :: BIG_COUNT
     Character(len=16)  :: lcnt,ltime
     REAL(DP) :: AVG_TIME

     write(lcnt,'(I8.8)')MYWATCH%LAP_COUNT

     BIG_COUNT = MYWATCH%COUNT_CURRENT - MYWATCH%COUNT_0
     IF (BIG_COUNT < 0) BIG_COUNT = BIG_COUNT + MYWATCH%COUNT_MAX

     AVG_TIME = DBLE(BIG_COUNT) / DBLE(MYWATCH%COUNT_RATE * MYWATCH%LAP_COUNT)

     write(ltime,'(F16.6)') AVG_TIME

     write(unit,*) "! === "//TRIM(msg)//" ==="
     write(unit,*) "! Average Lap Time(seconds):"//TRIM(adjustl(ltime))//"; &
          &Lap Count:"//TRIM(adjustl(lcnt))
     write(unit,*) "! ==============================="
     
  
   END SUBROUTINE WATCH_REPORT
!==============================================================================!
!  Print Watch Print                                                            !
!==============================================================================!
   SUBROUTINE WATCH_LAP_PRINT(MYWATCH,UNIT,MSG) 
     IMPLICIT NONE  
     

     TYPE(WATCH) :: MYWATCH
     integer, intent(in) :: unit
     Character(len=*), intent(in) :: MSG
     INTEGER(ITIME) :: BIG_COUNT
     REAL(DP) :: LAP_TIME
     Character(len=16)  :: ltime
     

     CALL WATCH_LAP_NOPRINT(MYWATCH)
     BIG_COUNT = MYWATCH%COUNT_CURRENT - MYWATCH%COUNT_LAST
     IF (BIG_COUNT < 0) BIG_COUNT = BIG_COUNT + MYWATCH%COUNT_MAX

     LAP_TIME = DBLE(BIG_COUNT) / DBLE(MYWATCH%COUNT_RATE)

     write(ltime,'(F16.6)') LAP_TIME

     write(unit,*) "! === "//TRIM(msg)//" ==="
     write(unit,*) "! Last Lap Time(seconds):"//TRIM(adjustl(ltime))
     write(unit,*) "! ==============================="
  
   END SUBROUTINE WATCH_LAP_PRINT

!==============================================================================!
!   Return a Time String Days:Hours:Minutes:Seconds from Number of Seconds     !
!==============================================================================!

   SUBROUTINE GETTIME(INSTRING,INSECS)

   IMPLICIT NONE
   INTEGER, INTENT(IN) :: INSECS 
   CHARACTER(LEN=13), INTENT(INOUT) :: INSTRING
   CHARACTER(LEN=4)  :: S0
   CHARACTER(LEN=2)  :: S1,S2,S3
   INTEGER :: DTCP,HTCP,MTCP,STCP
 
   DTCP = INSECS/(3600*24)
   HTCP = MOD(INSECS,(3600*24))/3600
   MTCP = MOD(INSECS,(3600))/60
   STCP = INSECS - (DTCP*3600*24 + HTCP*3600 + MTCP*60)

   IF(DTCP >= 10000 ) then
      
      INSTRING = "> 10000 DAYS"
      return
      
   else if (DTCP < 0) THEN
      
      INSTRING = " < 0 DAYS?"
      return
      
   else
      
      WRITE(S0,"(I4.4)")INT(DTCP)
      
   END IF
   
   WRITE(S1,"(I2.2)")INT(HTCP)
   
   WRITE(S2,"(I2.2)")INT(MTCP)
   
   WRITE(S3,"(I2.2)")INT(STCP)
   

   INSTRING = S0//":"//S1//":"//S2//":"//S3
   
   END SUBROUTINE GETTIME
!==============================================================================!
!  Get Current Time As Sting                                                   !
!==============================================================================!
SUBROUTINE GET_TIMESTAMP(TS)
  CHARACTER(LEN=*)  TS
  CHARACTER(LEN=8)  D
  CHARACTER(LEN=10) T
                                                                                      
  CALL DATE_AND_TIME ( DATE=D,TIME=T )
  TS = D(7:8)//'/'//D(5:6)//'/'//D(1:4)//'   '//T(1:2)//':'//T(3:4)
END SUBROUTINE GET_TIMESTAMP
!==============================================================================!
!   Report Calculation Speed and Time to Complete                              !
!==============================================================================!

   SUBROUTINE REPORT_TIME(IINT,ISTART,IEND,STIME)

   IMPLICIT NONE
   INTEGER(itime), INTENT(IN) :: IINT,ISTART,IEND
   Type(TIME),INTENT(IN) :: STIME
   REAL(SP) :: TTCP,TAVE
   CHARACTER(LEN=80)  :: SIMTIME,FINTIME 
   CHARACTER(LEN=22) :: PCOMP
   INTEGER :: I,ICMP,ICMP2
   INTEGER(ITIME) :: BIG_COUNT
   LOGICAL, SAVE :: INITIALIZED

   real(DP):: tmp, seconds
   integer :: hours, minutes
   Character(len=2) :: h, m
   Character(Len=9) :: s
   Character(Len=6) :: d

   IF(.not. Initialized) THEN
      CALL WATCH_INIT(Report_watch)

      TAVE = -1.0_SP
      FINTIME = "   UNKNOWN"
      Initialized = .true.
      
   ELSE
      !
      !  CALCULATE CURRENT TIME
      !
      CALL WATCH_LAP(Report_watch)
      
      !
      !  CALCULATE AVERAGE TIME/ITERATION
      !
      BIG_COUNT = REPORT_WATCH%COUNT_CURRENT - REPORT_WATCH%COUNT_0
      IF (BIG_COUNT < 0) BIG_COUNT = BIG_COUNT + REPORT_WATCH%COUNT_MAX
      
      TAVE = DBLE(BIG_COUNT) / DBLE(REPORT_WATCH%COUNT_RATE * REPORT_WATCH%LAP_COUNT)
      
      IF(TAVE < 0) TAVE = 0.0_SP
      
      !
      !  CALCULATE TIME TO COMPLETION
      !
      TTCP = TAVE*(IEND-IINT)
      
      CALL GETTIME(FINTIME,INT(TTCP))
      
   END IF

!
!  CALCULATE COMPLETION PERCENTAGE GRAPHIC
!
!   IF(MOD(IINT,1) /= 0) RETURN
   ICMP = INT( 100.*FLOAT(IINT-ISTART+1)/FLOAT(IEND-ISTART+1)) 
   ICMP2 = ICMP/5
   PCOMP = " "
   PCOMP(1:1)  = "|"
   PCOMP(22:22) = "|"

   DO I=2,ICMP2+1
     PCOMP(I:I) = "="
   END DO

   if (USE_REAL_WORLD_TIME) then
      !  Get the date and time of the current iterationr
      SIMTIME = Write_DateTime(STIME,6,"UTC")

      
      !
      !  REPORT
      !
      IF(MOD(IINT-1,10_DP) ==0) THEN
         WRITE(IPT,102)
      END IF
      WRITE(IPT,101)IINT,SIMTIME,FINTIME,TAVE,PCOMP
      
   else
      
      tmp = real(STIME%MuSOD,DP) / real(Million,DP)
      
      hours = tmp/3600
      minutes = (tmp-hours*3600)/60
      seconds = mod(tmp,60.0_DP)
      
      write(d,'(i6.6)') STIME%mjd
      write(h,'(i2.2)') hours
      write(m,'(i2.2)') minutes
      write(s,'(F9.6)') seconds

      SIMTIME = "D"//d//"T"//h//":"//m//":"//s
      IF(MOD(IINT-1,10_DP) ==0) THEN
         WRITE(IPT,103)
      END IF

      WRITE(IPT,101)IINT,SIMTIME,FINTIME,TAVE,PCOMP
   end if

   return
101 FORMAT(1X,"!",I9,1X,A26,3X,A13,3X,F8.4,2X,A22)  
102 FORMAT(1X,"!   IINT ",6X,"  SIMTIME(UTC)  ",9X," FINISH IN ",5X," SECS/IT ",1X," PERCENT COMPLETE   ")  
103 FORMAT(1X,"!   IINT ",6X,"  SIMTIME  ",14X," FINISH IN ",5X," SECS/IT ",1X," PERCENT COMPLETE   ")  

   END SUBROUTINE REPORT_TIME

     
END MODULE MOD_CLOCK
