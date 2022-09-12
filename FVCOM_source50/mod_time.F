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

MODULE Mod_Time
  USE MOD_PREC
  IMPLICIT NONE
   
  PUBLIC


! !PUBLIC DATA MEMBERS:
   INTEGER, PARAMETER :: itime = SELECTED_INT_KIND(18)

   INTEGER :: MPI_TIME

!  ! INTEL 9.1.X does not like initialized types
!   TYPE :: TIME
!      INTEGER(itime)  :: MJD     !! MODIFIED JULIAN DAY
!      INTEGER(itime)  :: MuSOD   !! MicroSECOND OF DAY
!   END TYPE TIME

   ! USE THIS FOR INTEL 10.1+ 
   TYPE :: TIME
      INTEGER(itime)  :: MJD   = 0   !! MODIFIED JULIAN DAY
      INTEGER(itime)  :: MuSOD = 0    !! MicroSECOND OF DAY
   END TYPE TIME

   INTERFACE ABS
      MODULE PROCEDURE ABS_TIME
   END INTERFACE


   INTERFACE MOD
      MODULE PROCEDURE modulo_time
   END INTERFACE

   INTERFACE ASSIGNMENT(=)
      MODULE PROCEDURE ASSIGN_TIME
   END INTERFACE

   INTERFACE OPERATOR(*)
      MODULE PROCEDURE TIME_x_int
      MODULE PROCEDURE TIME_x_long
      MODULE PROCEDURE TIME_x_flt
      MODULE PROCEDURE TIME_x_dbl
      MODULE PROCEDURE int_x_TIME
      MODULE PROCEDURE long_x_TIME
      MODULE PROCEDURE flt_x_TIME
      MODULE PROCEDURE dbl_x_TIME
   END INTERFACE

   INTERFACE OPERATOR(/)
      MODULE PROCEDURE TIME_div_int
      MODULE PROCEDURE TIME_div_long
      MODULE PROCEDURE TIME_div_flt
      MODULE PROCEDURE TIME_div_dbl

   END INTERFACE

   INTERFACE OPERATOR(+)
      MODULE PROCEDURE ADD_TIME
      MODULE PROCEDURE ADD_TIME_1
      MODULE PROCEDURE ADD_TIME_1a
      MODULE PROCEDURE ADD_TIME_a1
      MODULE PROCEDURE ADD_TIME_2
      MODULE PROCEDURE ADD_TIME_2a
      MODULE PROCEDURE ADD_TIME_a2
!      MODULE PROCEDURE ADD_TIME_STEP     
   END INTERFACE

   INTERFACE OPERATOR(-)
      MODULE PROCEDURE SUBTRACT_TIME
      MODULE PROCEDURE SUBTRACT_TIME_1
      MODULE PROCEDURE SUBTRACT_TIME_1a
      MODULE PROCEDURE SUBTRACT_TIME_a1
      MODULE PROCEDURE SUBTRACT_TIME_2
      MODULE PROCEDURE SUBTRACT_TIME_2a
      MODULE PROCEDURE SUBTRACT_TIME_a2
!      MODULE PROCEDURE SUBTRACT_TIME_STEP
   END INTERFACE


   INTERFACE OPERATOR(<=)
      MODULE PROCEDURE LE_TIME
   END INTERFACE

   INTERFACE OPERATOR(>=)
      MODULE PROCEDURE GE_TIME
   END INTERFACE

   INTERFACE OPERATOR(==)
      MODULE PROCEDURE EQ_TIME
   END INTERFACE

   INTERFACE OPERATOR(/=)
      MODULE PROCEDURE NE_TIME
   END INTERFACE

   INTERFACE OPERATOR(>)
      MODULE PROCEDURE GT_TIME
   END INTERFACE

   INTERFACE OPERATOR(<)
      MODULE PROCEDURE LT_TIME
   END INTERFACE

   INTERFACE DAYS2TIME
      MODULE PROCEDURE DAYS2TIME_DBL
      MODULE PROCEDURE DAYS2TIME_FLT
      MODULE PROCEDURE DAYS2TIME_INT
      MODULE PROCEDURE DAYS2TIME_LINT
   END INTERFACE
   INTERFACE SECONDS2TIME
      MODULE PROCEDURE SECONDS2TIME_DBL
      MODULE PROCEDURE SECONDS2TIME_FLT
      MODULE PROCEDURE SECONDS2TIME_INT
      MODULE PROCEDURE SECONDS2TIME_LINT
   END INTERFACE

   ! Seconds per day
   integer(itime), parameter :: spd         = 86400   
   integer(itime), parameter :: mspd        = spd * 1000
   integer(itime), parameter :: muspd       = mspd * 1000
   Integer(itime), parameter :: million     = 10**6
 CONTAINS
!======================================================
!======================================================
!======================================================   
   TYPE(TIME) FUNCTION ABS_TIME(A)
     IMPLICIT NONE
     TYPE(TIME), INTENT(IN) ::A

     CALL ADJUST(A)
     
     ABS_TIME%mjd = ABS(A%mjd)
     ABS_TIME%musod = ABS(A%musod)

   END FUNCTION ABS_TIME
!======================================================
   TYPE(TIME) FUNCTION MODULO_TIME(A,B)
! COMMENTS: THIS PROGRAM IS ROBUST BUT SLOW.....
     IMPLICIT NONE
     TYPE(TIME), INTENT(IN) ::A,B
     real(DP)               :: DPA,DPB,DIV,CNT
     TYPE(TIME) ::T, ZT

     ZT%mjd=0
     ZT%MuSOD=0

     DPA = DAYS(A)
     DPB = DAYS(B)
     DIV = DPA/DPB

     CNT = 1.0_DP

     IF(ABS(DIV) .GT. 1E5_DP) CNT = LOG10(ABS(DIV))

     IF( ABS(DIV) .GT. 1E16_DP) THEN
        ! THE DIVISOR EXCEEDS THE NUMERICAL ACCURACY - RETURN A BAD RESULT!
        MODULO_TIME%MJD = -HUGE(ZT%mjd)
        MODULO_TIME%MUSOD = -HUGE(ZT%musod)
        RETURN
     END IF

     

     IF (B .EQ. ZT) THEN
     ! SECOND ARGUMENT IS ZERO - ESCAPE
        MODULO_TIME%MJD = HUGE(ZT%mjd)
        MODULO_TIME%MUSOD = HUGE(ZT%musod)

     ELSE IF(ABS(A) .LT. ABS(B))THEN
     ! THE FIRST IS SMALLER THAN THE SECOND
        MODULO_TIME = A

     ELSE IF(ABS(A) .EQ. ABS(B))THEN
    ! THE FIRST IS SMALLER THAN THE SECOND
        MODULO_TIME = ZT
        
     ELSE IF (A .GT. ZT .and. B .GT. ZT) THEN
        ! BOTH ARE POSITIVE

        T = B * INT((DIV-CNT),ITIME)
        DO WHILE (A .GE. T+B)
           T = T +B
        END DO
        
        MODULO_TIME = A - T
     ELSE IF(A .LT. ZT .and. B .LT. ZT) THEN
        ! BOTH ARE NEGATIVE

        T = B * INT((DIV-CNT),ITIME)
        DO WHILE(A .LE. T+B)
           T = T +B
        END DO
        MODULO_TIME = A - T
     ELSE IF (A .LT. ZT .and. B .GT. ZT) THEN
        ! A IS NEGATIVE AND B IS POSITIVE

        T = B *  INT((DIV+CNT),ITIME)
        DO WHILE(A .LE. T-B)
           T = T -B
        END DO
        MODULO_TIME = A - T

     ELSE IF (A .GT. ZT .and. B .LT. ZT) THEN
        ! A IS POSITIVE AND B IS NEGATIVE

        T = B *  INT((DIV+CNT),ITIME)
        DO WHILE(A .GE. T-B)
           T = T-B
        END DO
        MODULO_TIME = A - T
     ELSE

        ! THIS SHOULD NEVER HAPPEN!
        MODULO_TIME%MJD = -HUGE(ZT%mjd)
        MODULO_TIME%MUSOD = -HUGE(ZT%musod)

     END IF

   END FUNCTION MODULO_TIME
!======================================================
   TYPE(TIME) FUNCTION DAYS2TIME_DBL(DAYS)
     implicit none
     TYPE(TIME) :: MJD
     real(DP), INTENT(IN)    :: DAYS
     real(DP)                :: TEMP

     MJD%mjd   = ANINT(DAYS,itime)
     TEMP = DAYS - ANINT(DAYS,itime)
     MJD%MuSOD = ANINT(temp*real(MUSPD,DP),itime)
     CALL ADJUST(MJD)
     DAYS2TIME_DBL = MJD
   END FUNCTION DAYS2TIME_DBL
!======================================================
   TYPE(TIME) FUNCTION DAYS2TIME_INT(DAYS)
     implicit none
     TYPE(TIME) :: MJD
     INTEGER, INTENT(IN)    :: DAYS

     MJD%mjd   = DAYS
     MJD%MuSOD = 0
     DAYS2TIME_INT = MJD
   END FUNCTION DAYS2TIME_INT
!======================================================
   TYPE(TIME) FUNCTION DAYS2TIME_LINT(DAYS)
     implicit none
     TYPE(TIME) :: MJD
     INTEGER(ITIME), INTENT(IN)    :: DAYS

     MJD%mjd   = DAYS
     MJD%MuSOD = 0
     DAYS2TIME_LINT = MJD
   END FUNCTION DAYS2TIME_LINT
!======================================================
   TYPE(TIME) FUNCTION DAYS2TIME_FLT(DAYS)
     implicit none
     real(SPA), INTENT(IN)    :: DAYS
     real(DP)                :: TEMP
     TEMP = DBLE(DAYS)
     DAYS2TIME_FLT = DAYS2TIME_DBL(TEMP)
   END FUNCTION DAYS2TIME_FLT
!======================================================
   TYPE(TIME) FUNCTION SECONDS2TIME_DBL(SECS)
     implicit none
     TYPE(TIME) :: MJD
     real(DP), INTENT(IN)    :: SECS
     real(DP)                :: TEMP

     MJD%mjd   = INT(SECS/DBLE(SPD),itime)
     TEMP = MOD(SECS,DBLE(SPD))
     MJD%MuSOD = ANINT(temp*real(million,DP),itime)
     CALL ADJUST(MJD)
     SECONDS2TIME_DBL = MJD
   END FUNCTION SECONDS2TIME_DBL
!======================================================
   TYPE(TIME) FUNCTION SECONDS2TIME_INT(SECS)
     implicit none
     TYPE(TIME) :: MJD
     INTEGER, INTENT(IN)    :: SECS

     MJD%mjd   = 0
     MJD%MuSOD = SECS
     CALL ADJUST(MJD)
     SECONDS2TIME_INT = MJD
   END FUNCTION SECONDS2TIME_INT
!======================================================
   TYPE(TIME) FUNCTION SECONDS2TIME_LINT(SECS)
     implicit none
     TYPE(TIME) :: MJD
     INTEGER(ITIME), INTENT(IN)    :: SECS

     MJD%mjd   = 0
     MJD%MuSOD = SECS
     CALL ADJUST(MJD)
     SECONDS2TIME_LINT = MJD
   END FUNCTION SECONDS2TIME_LINT
!======================================================
   TYPE(TIME) FUNCTION SECONDS2TIME_FLT(SECS)
     implicit none
     real(SPA), INTENT(IN)    :: SECS
     
     SECONDS2TIME_FLT = SECONDS2TIME_DBL(DBLE(SECS))

   END FUNCTION SECONDS2TIME_FLT
!======================================================
   FUNCTION TIME2NCITIME(MJD,RJD,D,MS) RESULT(res)
     implicit none
     INTEGER :: RES
     TYPE(TIME), INTENT(IN) :: MJD,RJD
     INTEGER, INTENT(OUT) :: D, MS
     REAL(DP) :: MSEC

     res = -1
     MSEC = dble(MJD%MuSod) / 1000.0_DP
     MS = ANINT(MSEC)

     ! CHECK TO MAKE SURE IT IS NOT TOO LARGE
     IF (ABS(MJD%MJD) .GT. HUGE(D)) THEN
        res =0 
        return
     END IF
     
     D = MJD%MJD - RJD%MJD

   END FUNCTION TIME2NCITIME
!======================================================
   FUNCTION NCITIME(D,MS) RESULT(MJD)
     implicit none
     TYPE(TIME) :: MJD
     INTEGER, INTENT(IN) :: D, MS
     
     MJD%MJD = D

     MJD%MuSod= INT(ms,ITIME)* INT(1000,ITIME)

   END FUNCTION NCITIME
!======================================================
   SUBROUTINE ADJUST(MJD)
     implicit none
     TYPE(TIME)     :: MJD
     integer(itime) :: musec
     integer(itime) :: idays
!     print*,"+++++++++++++begin adjust+++++++++++++++++"
!     write(*,*) "MJD in= ",MJD%MJD
!     write(*,*) "MUSOD in= ",MJD%MuSOD

     musec = mod(MJD%MuSOD, MuSPD)
!     print *,"musec= ", musec

     idays = (MJD%MuSOD - musec)/ MuSPD      
!     print *, "idays= ",idays
     MJD%mjd = MJD%mjd + idays
     MJD%MuSOD = musec 


!     call print_time(mjd,6,"intermediate")

     if (MJD%MuSOD .GT. 0 .AND. MJD%mjd .LT. 0) then
        MJD%mjd=MJD%mjd+1
        MJD%MuSOD=MJD%MuSOD-MuSPD
!        call print_time(mjd,6,"if... out")
        return
     else if (MJD%MuSOD .LT. 0 .AND.MJD%mjd .GT. 0  ) then
        MJD%MuSOD = MuSPD + MJD%MuSOD
        MJD%mjd = MJD%mjd -1
!        call print_time(mjd,6,"else if... out")
        return
     else
!        call print_time(mjd,6,"else... out")
        return
     end if
!     print*,"+++++++++++++END adjust+++++++++++++++++"
     
   END SUBROUTINE ADJUST
!======================================================
 TYPE(TIME) FUNCTION READ_TIME(timestr,status,TZONE)
   implicit none
   include 'fjulian.inc' 
   character(Len=*) :: timestr ! DO NOT USE INTENT ATTRIBUTE
   CHARACTER(LEN=*),OPTIONAL, INTENT(IN) :: TZONE
   integer status
   logical statl
   real(DP) :: SECS
   TYPE(TIME) :: DZONE

   
   if(timestr(1:1)=="-") then
      statl = FJul_ParseTime(timestr(2:),.False.,SECS)
      secs = -secs
   else if(timestr(1:1)=="+") then
      statl = FJul_ParseTime(timestr(2:),.False.,SECS)
   else
      statl = FJul_ParseTime(timestr,.False.,SECS)
   end if
   READ_TIME = seconds2time(SECS)
   status = 0
   if(statl) status = 1

   IF (PRESENT(TZONE)) THEN
      READ_TIME = READ_TIME - TIME_ZONE(TZONE,status)
   END IF

 END FUNCTION READ_TIME
!======================================================
 FUNCTION TIME_ZONE(TZONE,status) RESULT(DZONE)
   TYPE(TIME)                   :: DZONE
   CHARACTER(LEN=*), INTENT(IN) :: TZONE
   integer, intent(out)         :: status
   !! USE STANDARD NAMES FROM www.timeanddate.com
   !! Handle time zones!


   status = 1 
   SELECT CASE (TZONE)
   CASE("UTC")
      DZONE = seconds2time(0.0_DP)
   CASE("NONE")
      DZONE = seconds2time(0.0_DP)
   CASE("none")
      DZONE = seconds2time(0.0_DP)
   CASE("None")
      DZONE = seconds2time(0.0_DP)
   CASE("A")
      DZONE = READ_TIME('01:00:00',status)
   CASE("ACDT")
      DZONE = READ_TIME('10:30:00',status)
   CASE("ACST")
      DZONE = READ_TIME('09:30:00',status)
   CASE("ADT")
      DZONE = READ_TIME('-03:00:00',status)
   CASE("AEDT")
      DZONE = READ_TIME('11:00:00',status)
   CASE("AEST")
      DZONE = READ_TIME('10:00:00',status)
   CASE("AKDT")
      DZONE = READ_TIME('-08:00:00',status)
   CASE("AKST")
      DZONE = READ_TIME('-09:00:00',status)
   CASE("AST")
      DZONE = READ_TIME('-04:00:00',status)
   CASE("AWDT")
      DZONE = READ_TIME('09:00:00',status)
   CASE("AWST")
      DZONE = READ_TIME('08:00:00',status)
   CASE("B")
      DZONE = READ_TIME('02:00:00',status)
   CASE("BST")
      DZONE = READ_TIME('01:00:00',status)
   CASE("CDT")
      DZONE = READ_TIME('-05:00:00',status)
   CASE("CEDT")
      DZONE = READ_TIME('02:00:00',status)
   CASE("CEST")
      DZONE = READ_TIME('02:00:00',status)
   CASE("CET")
      DZONE = READ_TIME('01:00:00',status)
      !! SKIPPING MULTIPLE ABREVIATIONS 'CST'
   CASE("CST")
      DZONE = READ_TIME('-06:00:00',status)
   CASE("CXT")
      DZONE = READ_TIME('07:00:00',status)
   CASE("D")
      DZONE = READ_TIME('04:00:00',status)
   CASE("E")
      DZONE = READ_TIME('05:00:00',status)
   CASE("EDT")
      DZONE = READ_TIME('-04:00:00',status)
   CASE("EEDT")
      DZONE = READ_TIME('03:00:00',status)
   CASE("EEST")
      DZONE = READ_TIME('03:00:00',status)
   CASE("EET")
      DZONE = READ_TIME('02:00:00',status)
      !! AGAIN, SKIPPING MULTIPLE 'EST' defs, sorry australia!
   CASE("EST")
      DZONE = READ_TIME('-05:00:00',status)
   CASE("F")
      DZONE = READ_TIME('06:00:00',status)
   CASE("G")
      DZONE = READ_TIME('07:00:00',status)
   CASE("GMT")
      DZONE = seconds2time(0.0_DP)
   CASE("H")
      DZONE = READ_TIME('08:00:00',status)
   CASE("HAA")
      DZONE = READ_TIME('-03:00:00',status)
   CASE("HAC")
      DZONE = READ_TIME('-05:00:00',status)
   CASE("HADT")
      DZONE = READ_TIME('-09:00:00',status)
   CASE("HAE")
      DZONE = READ_TIME('-04:00:00',status)
   CASE("HAP")
      DZONE = READ_TIME('-07:00:00',status)
   CASE("HAR")
      DZONE = READ_TIME('-06:00:00',status)
   CASE("HAST")
      DZONE = READ_TIME('-10:00:00',status)
   CASE("HAT")
      DZONE = READ_TIME('-02:30:00',status)
   CASE("HAY")
      DZONE = READ_TIME('-08:00:00',status)
   CASE("HNA")
      DZONE = READ_TIME('-04:00:00',status)
   CASE("HNC")
      DZONE = READ_TIME('-06:00:00',status)
   CASE("HNE")
      DZONE = READ_TIME('-05:00:00',status)
   CASE("HNP")
      DZONE = READ_TIME('-08:00:00',status)
   CASE("HNR")
      DZONE = READ_TIME('-07:00:00',status)
   CASE("HNT")
      DZONE = READ_TIME('-3:30:00',status)
   CASE("HNY")
      DZONE = READ_TIME('-09:00:00',status)
   CASE("I")
      DZONE = READ_TIME('09:00:00',status)
   CASE("IST")
      DZONE = READ_TIME('01:00:00',status)
   CASE("K")
      DZONE = READ_TIME('10:00:00',status)
   CASE("L")
      DZONE = READ_TIME('11:00:00',status)
   CASE("M")
      DZONE = READ_TIME('12:00:00',status)
   CASE("MDT")
      DZONE = READ_TIME('-06:00:00',status)
   CASE("MESZ")
      DZONE = READ_TIME('02:00:00',status)
   CASE("MEZ")
      DZONE = READ_TIME('01:00:00',status)
   CASE("MST")
      DZONE = READ_TIME('-07:00:00',status)
   CASE("N")
      DZONE = READ_TIME('-01:00:00',status)
   CASE("NDT")
      DZONE = READ_TIME('-02:30:00',status)
   CASE("NFT")
      DZONE = READ_TIME('11:30:00',status)
   CASE("NST")
      DZONE = READ_TIME('-03:30:00',status)
   CASE("O")
      DZONE = READ_TIME('-02:00:00',status)
   CASE("P")
      DZONE = READ_TIME('-03:00:00',status)
   CASE("PDT")
      DZONE = READ_TIME('-07:00:00',status)
   CASE("PST")
      DZONE = READ_TIME('-08:00:00',status)
   CASE("Q")
      DZONE = READ_TIME('-04:00:00',status)
   CASE("R")
      DZONE = READ_TIME('-05:00:00',status)
   CASE("S")
      DZONE = READ_TIME('-06:00:00',status)
   CASE("T")
      DZONE = READ_TIME('-07:00:00',status)
   CASE("U")
      DZONE = READ_TIME('-08:00:00',status)
   CASE("V")
      DZONE = READ_TIME('-09:00:00',status)
   CASE("W")
      DZONE = READ_TIME('-10:00:00',status)
   CASE("WEDT")
      DZONE = READ_TIME('01:00:00',status)
   CASE("WEST")
      DZONE = READ_TIME('01:00:00',status)
   CASE("WET")
      DZONE = READ_TIME('00:00:00',status)
   CASE("WST")
      write(6,*)"TIMEZONE 'WST' IS AMBIGIOUS! USING WESTERN STANDARD TIME: UTC + 8:00"
      DZONE = READ_TIME('08:00:00',status)
   CASE("X")
      DZONE = READ_TIME('-11:00:00',status)
   CASE("Y")
      DZONE = READ_TIME('-12:00:00',status)
   CASE("Z")
      DZONE = READ_TIME('00:00:00',status)

   CASE DEFAULT
      DZONE = READ_TIME(TZONE,status)
   END SELECT

   if(status==0) write(6,*) "Time_zone failed! :: "//trim(tzone)
   
 END FUNCTION TIME_ZONE
!======================================================
 LOGICAL FUNCTION IS_VALID_TIMEZONE(timezone)
   IMPLICIT NONE
   character(Len=*), intent(in):: timezone
   TYPE(TIME) :: TEST
   integer    :: status

   IS_VALID_TIMEZONE=.false.
   TEST = TIME_ZONE(timezone,status)
   if (status==1) IS_VALID_TIMEZONE=.true.

 END FUNCTION IS_VALID_TIMEZONE
!======================================================
 TYPE(TIME)  FUNCTION READ_DATETIME(timestr,frmt,TZONE,status)
     ! RETURN VALUES FOR STATUS:
     !
     !  STATUS = -1 => SUCCESS
     !
     !  STATUS = 0  => FAILURE
     IMPLICIT NONE
     include 'fjulian.inc' 
     character(Len=*)                       :: timestr ! DO NOT USE INTENT ATTRIBUTE
     character(Len=*), intent(in)           :: frmt
     integer, intent(out)                   :: status
!     CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: TZONE
     CHARACTER(LEN=*), INTENT(IN) :: TZONE

     TYPE(TIME) :: mjd
     TYPE(TIME) :: DZONE
     integer	:: dutc, pos
     real(DP)   :: secs, tai,rmjd
     logical    :: statl
     
     statl = FJul_ParseDT(timestr, frmt, dutc, secs)
     if (.NOT. statl) then
        ! KNOWN WEEKNESS OF PARSE: CAN'T USE '_' AS A DELIMITER
        pos = index(timestr,'_')
        IF(pos == 0 ) return
        ! IF YOU FOUND '_' replace it with a ' '
        timestr(pos:pos) = ' '
        statl = FJul_ParseDT(timestr, frmt, dutc, secs)
        ! IF STILL UNSUCCESSFUL RETURN FAULT
        if(.not. statl) return
     end if
     status = 0 
     if(statl) status = 1
     mjd%MuSOD = ANINT(secs*million,itime)
     tai = FJul_TAIofDUTC(dutc)
     mjd%mjd = ANINT(FJul_MJDofTAI(tai, FJUL_UTC_TYPE),itime)

     READ_DATETIME = MJD - TIME_ZONE(TZONE,status)

     
   END FUNCTION READ_DATETIME
!======================================================
   character(len=80)  FUNCTION WRITE_DATETIME(mjdin,prec,TZONE) 
     IMPLICIT NONE
     include 'fjulian.inc' 
     ! PREC is the number of decimal seconds digits
     integer, intent(IN)                 :: prec
     TYPE(TIME), INTENT(IN)              :: mjdin
     CHARACTER(LEN=*), INTENT(IN) :: TZONE

     TYPE(TIME)  :: mjd
     real(DP)    :: rMJD
     real(DP)    :: tai
     real(DP)    :: secs
     integer :: dutc, status
     real(DP)    :: tmp1,tmp2

     MJD = MJDIN + TIME_ZONE(TZONE,status)


     tmp1 = dble(mjd%MuSOD)
     tmp2 = dble(MuSPD)

     rMJD = mjd%mjd + tmp1/tmp2
     
     secs = tmp1 / dble(million)
     
     tai  = FJul_TAIofMJD(rMJD, FJUL_UTC_TYPE)
     
     dutc = FJul_DUTCofTAI(tai,secs)
 

     call FJul_FormatPDS(dutc, secs, prec, .TRUE., WRITE_DATETIME)
     
   END FUNCTION WRITE_DATETIME
!======================================================
   TYPE(TIME) FUNCTION GET_NOW()
     IMPLICIT NONE
     include 'fjulian.inc' 
     CHARACTER(LEN=8)  D
     CHARACTER(LEN=10) T
     CHARACTER(LEN=5) Z
     CHARACTER(LEN=15) toff
     CHARACTER(LEN=25) TS
     TYPE(TIME) :: tzone
     integer :: dutc
     real(DP) :: secs
     integer status

     CALL DATE_AND_TIME ( DATE=D,TIME=T, ZONE=Z)

     ! GET TIME ZONE
     toff = Z(1:3)//":"//Z(4:5)

     ! GET TIME
     TS = D(1:4)//"/"//D(5:6)//"/"//D(7:8)// &
          & " "//T(1:2)//":"//T(3:4)//":"//T(5:8)

     GET_NOW = READ_DATETIME(TRIM(TS),'ymd',toff,status)

   END FUNCTION GET_NOW
!======================================================
   REAL(DP) FUNCTION SECONDS(MJD)
     IMPLICIT NONE
     TYPE(TIME), INTENT(IN):: MJD

     SECONDS = dble(MJD%mjd * SPD) + dble(MJD%MUSOD)/dble(million)
   END FUNCTION SECONDS
!======================================================
   REAL(DP) FUNCTION DAYS(MJD)
     IMPLICIT NONE
     TYPE(TIME), INTENT(IN):: MJD

     DAYS = dble(MJD%mjd) + dble(MJD%MUSOD)/dble(MUSPD)
   END FUNCTION DAYS
!======================================================
!--MODULE OPERATORS
!======================================================
   TYPE(TIME) FUNCTION int_x_Time(int,MJD)
     IMPLICIT NONE
     TYPE(TIME), INTENT(IN)  :: MJD
     integer, INTENT(IN) :: int
     
     int_x_time%MuSOD=MJD%MuSOD * int
     int_x_time%mjd=MJD%mjd * int


     call adjust(int_x_time)
   END FUNCTION INT_X_TIME
!======================================================
   TYPE(TIME) FUNCTION long_x_Time(long,MJD)
     IMPLICIT NONE
     TYPE(TIME), INTENT(IN)  :: MJD
     integer(itime), INTENT(IN) :: long
     
     long_x_time%MuSOD=MJD%MuSOD * long
     long_x_time%mjd=MJD%mjd * long


     call adjust(long_x_time)
   END FUNCTION LONG_X_TIME
!======================================================
   TYPE(TIME) FUNCTION Time_x_int(MJD,int)
     IMPLICIT NONE
     TYPE(TIME), INTENT(IN)  :: MJD
     integer, INTENT(IN) :: int
     
     Time_x_int%MuSOD=MJD%MuSOD * int
     Time_x_int%mjd=MJD%mjd * int


     call adjust(Time_x_int)
   END FUNCTION TIME_X_INT
!======================================================
   TYPE(TIME) FUNCTION Time_x_long(MJD,long)
     IMPLICIT NONE
     TYPE(TIME), INTENT(IN)  :: MJD
     integer(itime), INTENT(IN) :: long
     
     Time_x_long%MuSOD=MJD%MuSOD * long
     Time_x_long%mjd=MJD%mjd * long

     call adjust(Time_x_long)
   END FUNCTION TIME_X_LONG
!======================================================
   TYPE(TIME) FUNCTION Time_x_flt(MJD,flt)
     IMPLICIT NONE
     TYPE(TIME), INTENT(IN)  :: MJD
     real(SPA), INTENT(in) :: flt
     real(DP) :: DBL

     dbl = flt
     time_x_FLT = MJD * dbl

   END FUNCTION TIME_X_FLT
!======================================================
   TYPE(TIME) RECURSIVE FUNCTION Time_x_dbl(MJD,dbl)
     IMPLICIT NONE
     TYPE(TIME), INTENT(IN)  :: MJD
     real(DP), INTENT(in) :: dbl
     integer(itime) :: int
     real(DP) :: remainder
     real(DP) :: temp

     if (abs(dbl) .gt. 1.0_DP) THEN
        int = anint(dbl,itime)
        remainder = dbl - int
        Time_x_dbl = (MJD * int) + (MJD * remainder)
     else
        temp = real(MJD%MuSOD,DP) * real(dbl,DP)
        Time_x_dbl%MuSOD= anint(temp,itime)

        temp = real(MJD%mjd,DP) * real(dbl,DP)
        int = anint(temp,itime)
        Time_x_dbl%mjd=int
        Time_x_dbl%MuSOD = Time_x_dbl%MuSOD + (temp-int)*MUSPD 
     end if
     call adjust(Time_x_dbl)

   END FUNCTION TIME_X_DBL
!======================================================
   TYPE(TIME) FUNCTION flt_x_time(flt,MJD)
     IMPLICIT NONE
     TYPE(TIME), INTENT(IN)  :: MJD
     real(SPA), INTENT(in) :: flt
     real(DP) :: dbl

     dbl = dble(flt)
!   if(DBG_SET(dbg_log))  print*,'dbl=',dbl
     flt_x_time = MJD * dbl

   END FUNCTION FLT_X_TIME
!======================================================
   TYPE(TIME) FUNCTION dbl_x_time(dbl,MJD)
     IMPLICIT NONE
     TYPE(TIME), INTENT(IN)  :: MJD
     real(DP), INTENT(in) :: dbl

     DBL_X_TIME = MJD * DBL

   END FUNCTION DBL_X_TIME
!======================================================
   TYPE(TIME) FUNCTION Time_div_int(MJD,int)
     IMPLICIT NONE
     TYPE(TIME), INTENT(IN)  :: MJD
     integer, INTENT(IN) :: int
     real(DP) :: dbl

     dbl = dble(int)
     dbl = 1.0_DP / dbl
     TIME_DIV_INT = MJD * dbl

   END FUNCTION TIME_DIV_INT
!======================================================
   TYPE(TIME) FUNCTION Time_div_long(MJD,long)
     IMPLICIT NONE
     TYPE(TIME), INTENT(IN)  :: MJD
     integer(itime), INTENT(IN) :: long
     real(DP) :: dbl

     dbl = dble(long)
     dbl = 1.0_DP / dbl
     TIME_DIV_LONG = MJD * dbl

   END FUNCTION TIME_DIV_LONG
!======================================================
   TYPE(TIME) FUNCTION Time_div_flt(MJD,flt)
     IMPLICIT NONE
     TYPE(TIME), INTENT(IN)  :: MJD
     real(SPA), INTENT(IN) :: flt
     real(DP) :: dbl

     dbl = dble(flt)
     dbl = 1.0_DP / dbl
     TIME_DIV_flt = MJD * dbl

   END FUNCTION TIME_DIV_FLT
!======================================================
   TYPE(TIME) FUNCTION Time_div_dbl(MJD,dbl)
     IMPLICIT NONE
     TYPE(TIME), INTENT(IN)  :: MJD
     real(DP), INTENT(IN) :: dbl

     TIME_DIV_dbl = MJD * (1.0_DP / dbl)

   END FUNCTION TIME_DIV_DBL
!!$!======================================================
!!$   TYPE(TIME) FUNCTION SUBTRACT_TIME_STEP(time1,musec)
!!$     IMPLICIT NONE
!!$     TYPE(TIME), INTENT(IN)  :: time1
!!$     integer(itime), INTENT(IN) :: musec
!!$     
!!$     SUBTRACT_TIME_STEP = ADD_TIME_STEP(time1,-1*musec)
!!$
!!$   END FUNCTION SUBTRACT_TIME_STEP
!======================================================
   TYPE(TIME) FUNCTION ADD_TIME(time1,time2)
     IMPLICIT NONE
     TYPE(TIME), INTENT(IN)  :: time1, time2
     integer(itime) :: musec
     
     ADD_TIME%MuSOD = time1%MuSOD + time2%MuSOD
     ADD_TIME%mjd = time1%mjd + time2%mjd

     call adjust(ADD_TIME)
   END FUNCTION ADD_TIME
!======================================================
   FUNCTION ADD_TIME_1(time1,time2) RESULT(TSUM)
     IMPLICIT NONE
     TYPE(TIME), INTENT(IN)  :: time1(:), time2(:)
     TYPE(TIME), DIMENSION(SIZE(TIME1)) :: TSUM
     INTEGER :: I

     DO I = 1,SIZE(TIME1)
        TSUM(I) = TIME1(I) + TIME2(I)
     END DO
   END FUNCTION ADD_TIME_1
!======================================================
   FUNCTION ADD_TIME_1a(time1,time2) RESULT(TSUM)
     IMPLICIT NONE
     TYPE(TIME), INTENT(IN)  :: time1(:), time2
     TYPE(TIME), DIMENSION(SIZE(TIME1)) :: TSUM
     INTEGER :: I

     DO I = 1,SIZE(TIME1)
        TSUM(I) = TIME1(I) + TIME2
     END DO
   END FUNCTION ADD_TIME_1a
!======================================================
   FUNCTION ADD_TIME_a1(time1,time2) RESULT(TSUM)
     IMPLICIT NONE
     TYPE(TIME), INTENT(IN)  :: time1, time2(:)
     TYPE(TIME), DIMENSION(SIZE(TIME2)) :: TSUM
     INTEGER :: I

     DO I = 1,SIZE(TIME2)
        TSUM(I) = TIME1 + TIME2(I)
     END DO
   END FUNCTION ADD_TIME_a1
!======================================================
   FUNCTION ADD_TIME_2(time1,time2) RESULT(TSUM)
     IMPLICIT NONE
     TYPE(TIME), INTENT(IN)  :: time1(:,:), time2(:,:)
     TYPE(TIME), DIMENSION(SIZE(TIME1,1),size(time1,2)) :: TSUM
     INTEGER :: I,J

     DO I = 1,SIZE(TIME1,1)
        DO J = 1,SIZE(Time1,2)
           TSUM(I,J) = TIME1(I,J) + TIME2(I,J)
        END DO
     END DO
   END FUNCTION ADD_TIME_2
!======================================================
   FUNCTION ADD_TIME_2a(time1,time2) RESULT(TSUM)
     IMPLICIT NONE
     TYPE(TIME), INTENT(IN)  :: time1(:,:), time2
     TYPE(TIME), DIMENSION(SIZE(TIME1,1),size(time1,2)) :: TSUM
     INTEGER :: I,J

     DO I = 1,SIZE(TIME1,1)
        DO J = 1,SIZE(Time1,2)
           TSUM(I,J) = TIME1(I,J) + TIME2
        END DO
     END DO
   END FUNCTION ADD_TIME_2A
!======================================================
   FUNCTION ADD_TIME_a2(time1,time2) RESULT(TSUM)
     IMPLICIT NONE
     TYPE(TIME), INTENT(IN)  :: time1, time2(:,:)
     TYPE(TIME), DIMENSION(SIZE(TIME2,1),size(time2,2)) :: TSUM
     INTEGER :: I,J

     DO I = 1,SIZE(TIME2,1)
        DO J = 1,SIZE(Time2,2)
           TSUM(I,J) = TIME1 + TIME2(I,J)
        END DO
     END DO
   END FUNCTION ADD_TIME_A2
!======================================================
   TYPE(TIME) FUNCTION SUBTRACT_TIME(time1,time2)
     IMPLICIT NONE
     TYPE(TIME), INTENT(IN)  :: time1, time2
     
     SUBTRACT_TIME%MuSOD = time1%MuSOD - time2%MuSOD
     SUBTRACT_TIME%mjd = time1%mjd - time2%mjd

     call adjust(SUBTRACT_TIME)
   END FUNCTION SUBTRACT_TIME
!======================================================
   FUNCTION SUBTRACT_TIME_1(time1,time2) RESULT(TDIFF)
     IMPLICIT NONE
     TYPE(TIME), INTENT(IN)  :: time1(:), time2(:)
     TYPE(TIME), DIMENSION(SIZE(TIME1)) :: TDIFF
     INTEGER :: I

     DO I = 1,SIZE(TIME1)
        TDIFF(I) = TIME1(I) - TIME2(I)
     END DO
   END FUNCTION SUBTRACT_TIME_1
!======================================================
   FUNCTION SUBTRACT_TIME_1a(time1,time2) RESULT(TDIFF)
     IMPLICIT NONE
     TYPE(TIME), INTENT(IN)  :: time1, time2(:)
     TYPE(TIME), DIMENSION(SIZE(TIME2)) :: TDIFF
     INTEGER :: I

     DO I = 1,SIZE(TIME2)
        TDIFF(I) = TIME1 - TIME2(I)
     END DO
   END FUNCTION SUBTRACT_TIME_1A
!======================================================
   FUNCTION SUBTRACT_TIME_a1(time1,time2) RESULT(TDIFF)
     IMPLICIT NONE
     TYPE(TIME), INTENT(IN)  :: time1(:), time2
     TYPE(TIME), DIMENSION(SIZE(TIME1)) :: TDIFF
     INTEGER :: I

     DO I = 1,SIZE(TIME1)
        TDIFF(I) = TIME1(I) - TIME2
     END DO
   END FUNCTION SUBTRACT_TIME_A1
!======================================================
   FUNCTION SUBTRACT_TIME_2(time1,time2) RESULT(TDIFF)
     IMPLICIT NONE
     TYPE(TIME), INTENT(IN)  :: time1(:,:), time2(:,:)
     TYPE(TIME), DIMENSION(SIZE(TIME1,1),SIZE(TIME1,2)) :: TDIFF
     INTEGER :: I,J

     DO I = 1,SIZE(TIME1,1)
        DO J = 1,SIZE(TIME1,2)
           TDIFF(I,J) = TIME1(I,J) - TIME2(I,J)
        END DO
     END DO
   END FUNCTION SUBTRACT_TIME_2
!======================================================
   FUNCTION SUBTRACT_TIME_2a(time1,time2) RESULT(TDIFF)
     IMPLICIT NONE
     TYPE(TIME), INTENT(IN)  :: time1(:,:), time2
     TYPE(TIME), DIMENSION(SIZE(TIME1,1),SIZE(TIME1,2)) :: TDIFF
     INTEGER :: I,J

     DO I = 1,SIZE(TIME1,1)
        DO J = 1,SIZE(TIME1,2)
           TDIFF(I,J) = TIME1(I,J) - TIME2
        END DO
     END DO
   END FUNCTION SUBTRACT_TIME_2A
!======================================================
   FUNCTION SUBTRACT_TIME_a2(time1,time2) RESULT(TDIFF)
     IMPLICIT NONE
     TYPE(TIME), INTENT(IN)  :: time1, time2(:,:)
     TYPE(TIME), DIMENSION(SIZE(TIME2,1),SIZE(TIME2,2)) :: TDIFF
     INTEGER :: I,J

     DO I = 1,SIZE(TIME2,1)
        DO J = 1,SIZE(TIME2,2)
           TDIFF(I,J) = TIME1 - TIME2(I,J)
        END DO
     END DO
   END FUNCTION SUBTRACT_TIME_A2



!!$!======================================================
!!$   TYPE(TIME) FUNCTION  ADD_TIME_STEP(mjd,tstp)
!!$     IMPLICIT NONE
!!$     !INPUT PARAMETERS:
!!$     TYPE(TIME), intent(in)              :: mjd
!!$     integer(itime), intent(in)             :: tstp
!!$     
!!$     ADD_TIME_STEP%MuSOD =  mjd%MuSOD + tstp
!!$     ADD_TIME_STEP%mjd = MJD%mjd
!!$
!!$     call adjust(ADD_TIME_STEP)
!!$   END FUNCTION ADD_TIME_STEP
!======================================================
   SUBROUTINE ASSIGN_TIME(A,B)
     IMPLICIT NONE
     TYPE(TIME), INTENT(OUT) ::A
     TYPE(TIME), INTENT(IN)  ::B

     A%MJD = B%MJD
     A%MUSOD = B%MUSOD

   END SUBROUTINE ASSIGN_TIME
!======================================================   
                       ! time1 <= time2
   LOGICAL FUNCTION LE_TIME(time1,time2)
     IMPLICIT NONE
     TYPE(TIME), INTENT(IN)  :: time1, time2
     TYPE(TIME) dtime     
     LE_TIME = .false.
     dtime = time1 - time2
     if (dtime%MJD .lt. 0 .or. &
          & (dtime%MJD .EQ. 0 .and. dtime%MuSOD .LE. 0) ) LE_TIME = .TRUE.
   END FUNCTION LE_TIME
!======================================================
                       ! time1 < time2
   LOGICAL FUNCTION LT_TIME(time1,time2)
     IMPLICIT NONE
     TYPE(TIME), INTENT(IN)  :: time1, time2
     TYPE(TIME) dtime     
     LT_TIME = .false.
     dtime = time1 - time2
     if (dtime%MJD .lt. 0 .or. dtime%MuSOD .lt. 0) LT_TIME = .TRUE.
   END FUNCTION LT_TIME
!======================================================
                       ! time1 == time2
   LOGICAL FUNCTION EQ_TIME(time1,time2)
     IMPLICIT NONE
     TYPE(TIME), INTENT(IN)  :: time1, time2
     EQ_TIME = .false.
     if (time1%MJD .EQ. time2%MJD .and. &
          & time1%MuSOD .EQ. time2%MuSOD ) EQ_TIME = .TRUE.
   END FUNCTION EQ_TIME
!======================================================
                       ! time1 /= time2
   LOGICAL FUNCTION NE_TIME(time1,time2)
     IMPLICIT NONE
     TYPE(TIME), INTENT(IN)  :: time1, time2
     NE_TIME = .TRUE.
     if (EQ_TIME(TIME1,TIME2) ) NE_TIME = .FALSE.
   END FUNCTION NE_TIME
!======================================================
                       ! time1 >= time2
   LOGICAL FUNCTION GE_TIME(time1,time2)
     IMPLICIT NONE
     TYPE(TIME), INTENT(IN)  :: time1, time2
     TYPE(TIME) dtime     
     GE_TIME = .false.
     dtime = time1 - time2
     if (dtime%MJD .gt. 0 .or. &
          & (dtime%MJD .EQ. 0 .and. dtime%MuSOD .GE. 0) ) GE_TIME = .TRUE.
   END FUNCTION GE_TIME
!======================================================
                       ! time1 > time2
   LOGICAL FUNCTION GT_TIME(time1,time2)
     IMPLICIT NONE
     TYPE(TIME), INTENT(IN)  :: time1, time2
     TYPE(TIME) dtime     
     GT_TIME = .false.
     dtime = time1 - time2
     if (dtime%MJD .gt. 0 .or. dtime%MuSOD .gt. 0 ) GT_TIME = .TRUE.
   END FUNCTION GT_TIME
!======================================================

   SUBROUTINE PRINT_TIME(mjd,IPT,char)
     implicit none
     CHARACTER(Len=*), INTENT(IN) :: char
     INTEGER, INTENT(IN)   :: IPT
     TYPE(TIME),INTENT(IN) :: mjd
     real(DP)              :: tmp, seconds
     integer :: hours, minutes
     Character(len=3) :: h, m
     Character(Len=10) :: s
     Character(len=8) :: d

     tmp = real(mjd%MuSOD,DP) / real(Million,DP)
     
     hours = tmp/3600
     minutes = (tmp-hours*3600)/60
     seconds = mod(tmp,60.0_DP)

     write(d,'(i8.7)') mjd%mjd
     write(h,'(i3.2)') hours
     write(m,'(i3.2)') minutes
     write(s,'(F10.6)') seconds

     d = adjustl(d)
     h = adjustl(h)
     m = adjustl(m)
     s = adjustl(s)

     write(ipt,*)"!========"//trim(char)// "=========="
     write(ipt,*)"!    Day #    :", mjd%mjd
     write(ipt,*)"! MicroSecond #:", mjd%MuSOD
     write(ipt,*)"! (Human Time=d "//TRIM(d)//"::h"//TRIM(h)//":m"//TRIM(m)//":s"//TRIM(s)//")"
     write(ipt,*)"!=========================="
     
   END SUBROUTINE PRINT_TIME

   SUBROUTINE PRINT_REAL_TIME(mjd,IPT,char,TZONE)
     implicit none
     CHARACTER(Len=*), INTENT(IN) :: char
     INTEGER, INTENT(IN)   :: IPT
     TYPE(TIME),INTENT(IN) :: mjd
     CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: TZONE
     CHARACTER(LEN=120) :: string_local, string_utc


     IF(PRESENT(TZONE)) THEN
        ! IF THIS IS AN IDEAL TIME OR CASE THEN USE PRINT_TIME
        IF(TZONE=="none" .or. TZONE == "NONE") THEN
           CALL PRINT_TIME(MJD,IPT,CHAR)
           RETURN
        END IF
        
        string_local =  WRITE_DATETIME(mjd,6,TZONE)
        string_local = "! (Local Time="//trim(string_local)//"); time zone: "//trim(TZONE)
     END IF

     string_utc =  WRITE_DATETIME(mjd,6,"UTC")

     write(ipt,*)"!========"//trim(char)// "=========="
     write(ipt,*)"!    Day #    :", mjd%mjd
     write(ipt,*)"! MicroSecond #:", mjd%MuSOD
     write(ipt,*)"! (Date Time="//trim(string_utc)//")"
     IF(PRESENT(TZONE)) write(ipt,*) string_local
     write(ipt,*)"!=========================="
     
   END SUBROUTINE PRINT_REAL_TIME
!-----------------------------------------------------------------------
!!  Get the the month and total days of present month
!======================================================
   SUBROUTINE Now_2_month_days(TTime,Pyear,Pmonth,Pmdays)
     IMPLICIT NONE
     include 'fjulian.inc'
     TYPE(TIME), INTENT(IN)  ::TTime
     Integer,    INTENT(OUT) :: Pyear,Pmonth,Pmdays
     Integer :: Iyear,Imonth
!======================================================
     CHARACTER(LEN=8)  D
     CHARACTER(LEN=10) T
     CHARACTER(LEN=5) Z
     CHARACTER(LEN=15) toff
     CHARACTER(LEN=25) TS
     TYPE(TIME) :: tzone
     integer :: dutc
     real(DP) :: secs
     integer status
     TYPE(TIME)  ::Time1, Time2
     CHARACTER(LEN=120) :: string_local, string_utc
     CHARACTER(LEN=10) :: string_tt

     !CALL DATE_AND_TIME ( DATE=D,TIME=T, ZONE=Z)

       !string_utc =  WRITE_DATETIME(mjd,6,"UTC")
       string_local =  WRITE_DATETIME(ttime,6,"UTC")

      string_tt=trim(string_local)
      D(1:4)=string_tt(1:4)
      D(5:6)=string_tt(6:7)
      D(7:8)=string_tt(9:10)

     ! GET TIME ZONE
     !toff = Z(1:3)//":"//Z(4:5)
     toff = '0:00' 

     !!  get the month from D(5:6)
     ! TS1
     read(D(1:4),'(I4)') IYear
     read(D(5:6),'(I2)') Pmonth

     Pyear = IYear

     D(7:8)='01'
     T(1:8)='00000000'

     ! First day of Present month
     TS = D(1:4)//"/"//D(5:6)//"/"//D(7:8)// &
          & " "//T(1:2)//":"//T(3:4)//":"//T(5:8)

     Time1 = READ_DATETIME(TRIM(TS),'ymd',toff,status)

     ! First day of next month
       Imonth=Pmonth+1
       if(Imonth>12) then  !! to next year
          Imonth=1
          Iyear =Iyear+1
       endif 
       write(D(1:4),'(I4.4)') Iyear
       write(D(5:6),'(I2.2)') Imonth
       
       
     TS = D(1:4)//"/"//D(5:6)//"/"//D(7:8)// &
          & " "//T(1:2)//":"//T(3:4)//":"//T(5:8)

     Time2 = READ_DATETIME(TRIM(TS),'ymd',toff,status)

     Pmdays = Time2%mjd -Time1%mjd

   END SUBROUTINE Now_2_month_days
!!

   SUBROUTINE Now_2_days_test
     IMPLICIT NONE
     include 'fjulian.inc'
     Integer :: Iyear,Imonth
     Integer :: Pmonth,Pmdays
!======================================================
     CHARACTER(LEN=8)  D
     CHARACTER(LEN=10) T
     CHARACTER(LEN=5) Z
     CHARACTER(LEN=15) toff
     CHARACTER(LEN=25) TS
     TYPE(TIME) :: tzone
     integer :: dutc
     real(DP) :: secs
     integer status
     TYPE(TIME)  ::Time1, Time2
     CHARACTER(LEN=120) :: string_local, string_utc
     CHARACTER(LEN=10) :: string_tt

     !CALL DATE_AND_TIME ( DATE=D,TIME=T, ZONE=Z)

!Michael Dunphy (Michael.Dunphy@dfo-mpo.gc.ca)
!A memory checker found this 10-byte write to an 8-byte string
!     D='1990-01-01'
     D='19900101'
!Michael Dunphy
     toff = '0:00'

     ! TS1
      IYear=01
      Pmonth=01
     write(D(1:4),'(I4.4)') IYear
     write(D(5:6),'(I2.2)') Pmonth

     D(7:8)='01'
     T(1:8)='00000000'

     ! First day of Present month
     TS = D(1:4)//"/"//D(5:6)//"/"//D(7:8)// &
          & " "//T(1:2)//":"//T(3:4)//":"//T(5:8)

     Time1 = READ_DATETIME(TRIM(TS),'ymd',toff,status)

     ! TS1
      IYear=01995
      Pmonth=01
     write(D(1:4),'(I4.4)') IYear
     write(D(5:6),'(I2.2)') Pmonth

     D(7:8)='01'
     T(1:8)='00000000'

     ! First day of Present month
     TS = D(1:4)//"/"//D(5:6)//"/"//D(7:8)// &
          & " "//T(1:2)//":"//T(3:4)//":"//T(5:8)

     Time1 = READ_DATETIME(TRIM(TS),'ymd',toff,status)


     end subroutine Now_2_days_test
 
end module Mod_Time

