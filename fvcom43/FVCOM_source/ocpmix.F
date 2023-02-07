!     Ocean Pack miscellaneous routines
!
!     real function DTTIME
!     subroutine    DTINTI
!     subroutine    DTRETI
!     char function DTTIWR
!     REPARM
!     INAR2D
!     STRACE
!     MSGERR
!     TABHED
!     FOR
!     logical function EQREAL /* Checks whether REAL1 is appr.            
!                                equal to REAL2                 */        
!     LSPLIT            /* splits an input line into data items */        
!     BUGFIX                                                              
!     COPYCH (copied from file OCPDPN)                                    
!
!*******************************************************************
!                                                                  *
      REAL FUNCTION DTTIME (INTTIM)
!                                                                  *
!*******************************************************************
!
      USE OCPCOMM1                                                        
      USE OCPCOMM2                                                        
      USE OCPCOMM3                                                        
      USE OCPCOMM4                                                        
!
      IMPLICIT NONE
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering                              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmers: R.C. Ris, N. Booij,                          |
!     |              IJ.G. Haagsma, A.T.M.M. Kieftenburg,         |
!     |              M. Zijlema, E.E. Kriezi,                     |
!     |              R. Padilla-Hernandez, L.H. Holthuijsen       |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 2004-2005  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. Authors
!
!     30.74: IJsbrand Haagsma (Include version)
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!      9705, May  97: month number is checked
!     30.74, Nov. 97: Prepared for version with INCLUDE statements
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     DTTIME gives time in seconds from a reference day
!            it also initialises the reference day
!
!  3. Method
!
!     every fourth year is a leap-year, but not the century-years, however
!     also leap-years are: year 0, 1000, 2000 etc.
!     1 jan of year 0 is daynumber 1.
!
!  4. Argument variables
!
!     INTTIM(1): year
!           (2): month
!           (3): day
!           (4): hour
!           (5): minute
!           (6): second
!
      INTEGER INTTIM(6)                                                   
!
!  5. PARAMETER VARIABLES
!
!  6. LOCAL VARIABLES
!
!     IDYMON : number of days of each month (February counts as 28 days)
!     IYEAR  : number of years after substacking the centuries
!     IYRM1  : ??
!     IDNOW  : ??
!     I      : ??
!     II     : ??
!
      INTEGER IDYMON(12), IYEAR, IYRM1, IDNOW, I, II
!
!     LEAPYR : Whether year in INTTIM(1) is a leapyear
!     LOGREF : ??
!
      LOGICAL LEAPYR, LOGREF
!
!     REFDAY  day number of the reference day; the reference time is 0:00
!            of the reference day; the first day entered is used as
!             reference day.
!
!
!  8. SUBROUTINE USED
!
!  9. SUBROUTINES CALLING
!
! 10. ERROR MESSAGES
!
! 11. REMARKS
!
! 12. STRUCTURE
!
! 13. SOURCE TEXT
!
      SAVE LOGREF, IDYMON
      DATA IDYMON /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
      DATA LOGREF /.FALSE./
!
      IYEAR = INTTIM(1)
      IYRM1 = IYEAR-1
      LEAPYR=(MOD(IYEAR,4) == 0 .AND. MOD(IYEAR,100) /= 0) .OR.     &
              MOD(IYEAR,1000) == 0
      IDNOW=0
      IF(INTTIM(2) > 12)THEN                                           
        WRITE (PRINTF, 8) INTTIM(2), (INTTIM(II), II=1,6)                 
   8    FORMAT (' erroneous month ', I2, ' in date/time ', 6I4)           
      ELSE IF(INTTIM(2) > 1)THEN                                       
        DO 10 I = 1,INTTIM(2)-1
          IDNOW=IDNOW+IDYMON(I)
  10    CONTINUE
      ENDIF                                                               
      IDNOW=IDNOW+INTTIM(3)
      IF(LEAPYR .AND. INTTIM(2) > 2) IDNOW=IDNOW+1
      IDNOW = IDNOW + IYEAR*365 + IYRM1/4 - IYRM1/100 + IYRM1/1000 + 1
      IF(IYEAR == 0) IDNOW=IDNOW-1
      IF(.NOT.LOGREF)THEN
        REFDAY = IDNOW
        LOGREF = .TRUE.
        DTTIME = 0.
      ELSE
        DTTIME = REAL(IDNOW-REFDAY) * 24.*3600.
      ENDIF
      DTTIME = DTTIME + 3600.*REAL(INTTIM(4)) + 60.*REAL(INTTIM(5)) +     &
               REAL(INTTIM(6))
      RETURN
      END FUNCTION DTTIME
!*****************************************************************
!                                                                *
   SUBROUTINE INAR2D (ARR, MGA,                                  &
                      NDSL,                                      &
		      NDSD, IDFM, RFORM,                         &
                      IDLA, VFAC,                                &
		      NHED, NHEDF)
!                                                                *
!*****************************************************************
!
   USE OCPCOMM1                                                        
   USE OCPCOMM2                                                        
   USE OCPCOMM3                                                        
   USE OCPCOMM4                                                        
!
   IMPLICIT NONE
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering                              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmers: R.C. Ris, N. Booij,                          |
!     |              IJ.G. Haagsma, A.T.M.M. Kieftenburg,         |
!     |              M. Zijlema, E.E. Kriezi,                     |
!     |              R. Padilla-Hernandez, L.H. Holthuijsen       |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 2004-2005  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. Authors
!
!     30.72: IJsbrand Haagsma
!     30.74: IJsbrand Haagsma (Include version)
!     30.82: IJsbrand Haagsma
!     34.01: Jeroen Adema
!     40.00: Nico Booij
!     40.02: IJsbrand Haagsma
!     40.03: Nico Booij
!     40.08: Erick Rogers
!     40.13: Nico Booij
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     01.05, Feb. 90: Before reading values in the array are divided by VFAC,
!                     in order to retain correct values for points where no
!                     value was given
!     01.06, Apr. 91: i/o status is printed if read error occurs
!     30.72, Sept 97: Changed DO-block with one CONTINUE to DO-block with
!                     two CONTINUE's
!     30.72, Sept 97: Corrected reading of heading lines for SERIES of files
!                     in dynamic mode
!     30.74, Nov. 97: Prepared for version with INCLUDE statements
!     40.00, July 98: SWAN specific statements modified
!                     unformatted read: heading lines also read unformatted
!                     distinction between NDSD (data file) and NDSL (file list)
!     30.82, Sep. 98: Added INQUIRE statement to produce correct file name in
!                     case of a read error
!     34.01, Feb. 99: Introducing STPNOW
!     40.02, Sep. 00: Replaced computed GOTO with CASE construct
!     40.02, Sep. 00: Replaced reserved words IOSTAT with IOERR and STATUS with IERR
!     40.03, Jul. 00: END= added to READ statement for correct reading of series
!                     of files
!     40.03, Jul. 00: TRIM used to improve readability of message
!     40.13, Apr. 01: END=930 added in READ statement; corresponding error message added
!     40.08, Mar. 03: Changed an INQUIRE statement so that it does not produce
!                     misleading results.
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Reads a 2d array from dataset
!     is used to read e.g. bathymetry, one component of wind velocity
!
!  3. METHOD
!
!  4. ARGUMENT VARIABLES
!
!     IDFM   : input    format index
!     IDLAM  : input    lay-out indicator
!     MXA    : input    number of points along x-side of grid
!     MYA    : input    number of points along y-side of grid
!     NDSD   : input    unit number of the file from which to read the dataset
!     NDSL   : input    unit number of the file containing the list of filenames
!     NHEDF  : input    number of heading lines in the file (first lines).
!     NHEDL  : input    number of heading lines in the file
!                       before each array
!
   INTEGER   IDFM, IDLA, MGA, NDSD, NDSL, NHED, NHEDF
!
!     ARR    : input    results appear in this array
!     RFORM  : input    format used in reading data (char. string)
!     VFAC   : input    factor by which data must be multiplied.
!
   REAL      ARR(MGA), VFAC
!
   CHARACTER RFORM *(*)
!
!  5. PARAMETER VARIABLES
!
!  6. LOCAL VARIABLES
!
!     IERR   : ??
!     IENT   : number of entries into this subroutine
!     IOERR  : input   0 : Full messages printed
!                      -1: Only error messages printed
!                      -2: No messages printed
!              output  error indicator
!     IH     : ??
!     IX     : ??
!     IY     : ??
!     NUMFIL : ??
!
   INTEGER   IERR, IENT, IOERR, IH, IX, IY, NUMFIL                     
!
!     HEDLIN : Content of a header line
!
   CHARACTER HEDLIN *80
!
!  8. SUBROUTINE USED
!
   LOGICAL STPNOW                                                      
!
!  9. SUBROUTINES CALLING
!
! 10. ERROR MESSAGES
!
! 11. REMARKS
!
! 12. STRUCTURE
!
! 13. SOURCE TEXT
!
   INTEGER IG
   SAVE IENT
   DATA IENT /0/
   CALL STRACE (IENT, 'INAR2D')
!
999  IF(NDSD < 0) RETURN                                               
!     no reading from file due to open error
!
!     *** NUMFIL is the number of that is open in one time step  **
      NUMFIL = 0                                                          
!      IF(ITEST >= 100)THEN
!        WRITE (PRINTF, 12) MXA, MYA, NDSD, IDFM, RFORM,                   40.00
!     &  IDLA, VFAC, NHED
!  12    FORMAT (' * TEST INAR2D *', 4I4, 1X, A16, I3, 1X, E12.4, I3)
!      ENDIF
!
!     Read heading lines, and print the same:
!
  11  IF (NHED.GT.0) THEN
        IF (IDFM.LT.0) THEN                                               
          IF (ITEST.GE.30)                                      &
                  WRITE (PRINTF, '(I3,A)') NHED, ' Heading lines'        
          DO 28 IH=1, NHED
            READ (NDSD, END=910)                                          
  28      CONTINUE
        ELSE
          DO 30 IH=1, NHED
            READ (NDSD, '(A80)', end=910) HEDLIN                          
            IF (IH.EQ.1) WRITE (PRINTF, '(A)') ' **  Heading lines  **'
            WRITE (PRINTF, '(A4,A80)') ' -> ', HEDLIN
  30      CONTINUE
        ENDIF
      ENDIF
!
!     divide existing values in the array by VFAC
!
      DO IG = 1, MGA                                                   
        ARR(IG) = ARR(IG) / VFAC
      END DO	
!
!     start reading of 2D-array
!
      READ(NDSD, END=910, ERR=920, IOSTAT=IERR) (ARR(IG), IG=1,MGA)       
      GOTO 900                                                           
!
!     *** End of data file, in case SERIES next file is opened
!     *** unit = NDSD is closed before the next one is opened
!
 910  CONTINUE
      CLOSE(NDSD)
      NUMFIL = NUMFIL + 1
      IF (NUMFIL .GE. 2) GO TO 911
      IF (NDSL.GT.0) THEN
        READ (NDSL, '(A)', END=930) FILENM                                
        IF (IDFM.NE.-1) THEN
          IOERR = 0
          CALL FOR (NDSD, FILENM, 'OF', IOERR)                            
          IF (STPNOW()) RETURN                                            
        ELSE
          IOERR = 0
          CALL FOR (NDSD, FILENM, 'OU', IOERR)                            
          IF (STPNOW()) RETURN                                            
        ENDIF
!
!       Read heading lines, and print these:
!                                                                         
  2     IF (NHEDF.GT.0) THEN                                              
          IF (IDFM.LT.0) THEN                                             
            IF (ITEST.GE.30) WRITE (PRINTF, '(I3,A,A)') NHEDF,      &
                 ' Heading lines at begin of file ', TRIM(FILENM)        
            DO 828 IH=1, NHEDF                                            
              READ (NDSD)                                                 
 828        CONTINUE
          ELSE                                                            
            WRITE (PRINTF, '(A,A,A)') ' **  Heading lines file ',   &
           TRIM(FILENM), ' **'                                           
            DO 830 IH=1, NHEDF                                            
              READ (NDSD, '(A80)') HEDLIN                                 
              WRITE (PRINTF, '(A4,A80)') ' -> ', HEDLIN                   
 830        CONTINUE                                                      
          ENDIF                                                           
        ENDIF                                                             
        GO TO 11
      ENDIF
!
!     error message when end of file is encountered
!
!     --- initialize FILENM so that previous value is not used            40.08
!         in case unit NDSD does not exist                                40.08
 911  FILENM='DUMMY'
!     --------------------------------------------------------------------40.08
!     THIS INQUIRE STATEMENT IS PROBLEMATIC, SINCE (AT LEAST              40.08
!     SOMETIMES) NDSD HAS ALREADY BEEN CLOSED, SO THE INQUIRE             40.08
!     STATEMENT SHOULD NOT WORK.                                          40.08
!     --------------------------------------------------------------------40.08
      INQUIRE (UNIT=NDSD, NAME=FILENM)
      CALL MSGERR (2, 'Unexpected end of file while reading '//       &
                      TRIM(FILENM))                                      
      NDSD = 0                                                            
      IDLA = -1
!     Value of IDLA=-1 signals end of file to calling program
!
      GOTO 900
!
!     --- initialize FILENM                                               
 920  FILENM='DUMMY'                                                      
      INQUIRE (UNIT=NDSD, NAME=FILENM)                                    
      CALL MSGERR (2, 'Error while reading file '//TRIM(FILENM))          
      WRITE (PRINTF, 922) IERR                                            
 922  FORMAT (' i/o status ', I6)                                         
      IDLA = -2                                                           
!     Value of IDLA=-2 signals read error to calling program
!
!     Multiply all values in the array by VFAC
!
 900  DO IG = 1, MGA                                                  
        ARR(IG) = ARR(IG) * VFAC
      END DO	
!
 990  IF (ITEST.GE.100 .OR. IDLA.LT.0) THEN
!        DO 996 IY=MYA, 1, -1
!          WRITE (PRINTF, 994) (ARR(IX,IY), IX=1,MXA)
! 994      FORMAT ((1X, 10E12.4))
! 996    CONTINUE
      ENDIF
      RETURN

!     No more files in NDSL:
!     --- initialize FILENM                                               
 930  FILENM='DUMMY'                                                      
      INQUIRE (UNIT=NDSL, NAME=FILENM)                                    
      CALL MSGERR (2, 'Series of input files ended in '//TRIM(FILENM))    
      RETURN                                                              

      END SUBROUTINE INAR2D
!
!*****************************************************************
!                                                                *
   SUBROUTINE STRACE (IENT, SUBNAM)
!                                                                *
!*****************************************************************
!
   USE OCPCOMM1                                                        
   USE OCPCOMM2                                                        
   USE OCPCOMM3                                                        
   USE OCPCOMM4                                                        
   USE M_PARALL                                                        
!
   IMPLICIT NONE
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering                              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmers: R.C. Ris, N. Booij,                          |
!     |              IJ.G. Haagsma, A.T.M.M. Kieftenburg,         |
!     |              M. Zijlema, E.E. Kriezi,                     |
!     |              R. Padilla-Hernandez, L.H. Holthuijsen       |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 2004-2005  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. AUTHORS
!
!     40.41: Marcel Zijlema
!
!  1. UPDATES
!
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. PURPOSE
!
!     This subroutine produces depending on the value of 'ITRACE'
!     a message containing the name 'SUBNAM'. the purpose of this
!     action is to detect the entry of a subroutine.
!
!  3. METHOD
!
!     the first executable statement of subroutine 'AAA' has to
!     be : CALL STRACE(IENT,'AAA')
!     further is necessary : DATA IENT/0/
!     IF ITRACE=0, no message
!     IF ITRACE>0, a message is printed up to ITRACE times
!
!  4. ARGUMENT VARIABLES
!
!     IENT   :  i/o    Number of entries into the calling subroutine
!
   INTEGER IENT
!
!     SUBNAM :  inp    name of the calling subroutine.
!
   CHARACTER SUBNAM *(*)
!
!  5. PARAMETER VARIABLES
!
!  6. LOCAL VARIABLES
!                                                                         40.31
!$    LOGICAL,EXTERNAL :: OMP_IN_PARALLEL                                 40.31
!
!  8. SUBROUTINE USED
!
!  9. SUBROUTINES CALLING
!
! 10. ERROR MESSAGES
!
! 11. REMARKS
!
! 12. STRUCTURE
!
! 13. SOURCE TEXT
!
   IF(ITRACE == 0) RETURN
   IF(IENT > ITRACE) RETURN
!$ IF(OMP_IN_PARALLEL())THEN                                         
!$OMP MASTER                                                              
!$   IENT=IENT+1                                                      
!$   WRITE (PRTEST, 10) SUBNAM                                        
!$   IF(SCREEN /= PRINTF) WRITE (SCREEN, 10) SUBNAM                  
!$OMP END MASTER                                                          
!$ ELSE                                                                
     IENT=IENT+1
     WRITE (PRTEST, 10) SUBNAM
     IF(SCREEN /= PRINTF .AND. INODE == MASTER) WRITE (SCREEN, 10) SUBNAM  
!$ ENDIF                                                               
10 FORMAT (' ++ trace subr: ',A)
   RETURN
   END SUBROUTINE STRACE
 
!*****************************************************************
!                                                                *
   SUBROUTINE MSGERR (LEV,STRING)
!                                                                *
!*****************************************************************
!
   USE OCPCOMM1                                                        
   USE OCPCOMM2                                                        
   USE OCPCOMM3                                                        
   USE OCPCOMM4                                                        
   USE M_PARALL                                                        
!
   IMPLICIT NONE
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering                              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmers: R.C. Ris, N. Booij,                          |
!     |              IJ.G. Haagsma, A.T.M.M. Kieftenburg,         |
!     |              M. Zijlema, E.E. Kriezi,                     |
!     |              R. Padilla-Hernandez, L.H. Holthuijsen       |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 2004-2005  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. AUTHORS
!
!     40.02: IJsbrand Haagsma
!     40.03, 40.13: Nico Booij
!     40.30: Marcel Zijlema
!     40.41: Marcel Zijlema
!
!  1. UPDATES
!
!     40.03, Aug. 00: variable ERRFNM introduced in order to get correct
!                     message on UNIX system
!     40.02, Sep. 00: Removed STOP statement
!     40.13, Nov. 01: OPEN statement instead of CALL FOR
!                     to prevent recursive subroutines calling
!     40.30, Jan. 03: introduction distributed-memory approach using MPI
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. PURPOSE
!
!     Error messages are produced by subroutine MSGERR. if necessary
!     the value of LEVERR is increased.
!     In case of a high error level an error message file is opened
!
!  3. METHOD
!
!  4. ARGUMENT VARIABLES
!
!     LEV    : indicates how severe the present error is
!     STRING : contents of the present error message
!
   INTEGER   LEV
!
   CHARACTER STRING*(*)
!
!  5. PARAMETER VARIABLES
!
!  6. LOCAL VARIABLES
!
!     IERR   : if non-zero error message file was already opened unsuccessfully
!     IERRF  : unit reference number of the error message file
!     ILPOS  : actual length of error message filename
!
   INTEGER, SAVE :: IERR=0, IERRF=0                                    
   INTEGER ILPOS                                                       
!
!     ERRM   : error message prefix
!
   CHARACTER (LEN=17) :: ERRM                                          
!
!     ERRFNM : name of error message file
!
   CHARACTER (LEN=LENFNM), SAVE :: ERRFNM = 'Errfile'                  
!
!  8. SUBROUTINE USED
!
!     ---
!
!  9. SUBROUTINES CALLING
!
! 10. ERROR MESSAGES
!
! 11. REMARKS
!
! 12. STRUCTURE
!
! 13. SOURCE TEXT
!
!
   IF(LEV > LEVERR) LEVERR=LEV
   IF(LEV == 0)THEN
     ERRM = 'Message          '
   ELSE IF(LEV == 1)THEN
     ERRM = 'Warning          '
   ELSE IF(LEV == 2)THEN
     ERRM = 'Error            '
   ELSE IF(LEV == 3)THEN
     ERRM = 'Severe error     '
   ELSE
     ERRM = 'Terminating error'
   ENDIF
   WRITE (PRINTF,12) ERRM, STRING
12 FORMAT (' ** ', A, ': ',A)
   IF(LEV > MAXERR)THEN
     IF(IERRF == 0)THEN
       IF(IERR /= 0) RETURN
!
!      append node number to ERRFNM in case of                         
!      parallel computing                                              
!
       IF(PARLL)THEN                                                 
         ILPOS = INDEX ( ERRFNM, ' ' )-1                              
         WRITE(ERRFNM(ILPOS+1:ILPOS+4),13) INODE                      
13       FORMAT('-',I3.3)                                             
       END IF                                                          
!
       IERRF = 17                                                          
       OPEN (UNIT=IERRF, FILE=ERRFNM, FORM='FORMATTED')                    
     ENDIF
     WRITE (IERRF,14) ERRM, STRING
14   FORMAT (A, ': ',A)
   ENDIF
!
   RETURN
!
   END SUBROUTINE MSGERR
 
!
!*****************************************************************
!                                                                *
   LOGICAL FUNCTION STPNOW()                                          
!                                                                *
!*****************************************************************
!
   USE OCPCOMM4                                                      
!
   IMPLICIT NONE
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering                              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmers: R.C. Ris, N. Booij,                          |
!     |              IJ.G. Haagsma, A.T.M.M. Kieftenburg,         |
!     |              M. Zijlema, E.E. Kriezi,                     |
!     |              R. Padilla-Hernandez, L.H. Holthuijsen       |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 2004-2005  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. Authors
!
!     30.82, Feb. 99: IJsbrand Haagsma
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     30.82: New function
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Function determines wheter the SWAN program should be stopped
!     due to a terminating error
!
!  3. Method
!
!     Compares two common variables (the maximum allowable error-level,
!     MAXERR and the actual error-level: LEVERR).
!
!  4. ARGUMENT VARIABLES
!
!  5. PARAMETER VARIABLES
!
!  6. LOCAL VARIABLES
!
!     IENT  : Number of entries into this subroutine
!
   INTEGER IENT
!
!  8. SUBROUTINE USED
!
!  9. SUBROUTINES CALLING
!
! 10. ERROR MESSAGES
!
! 11. REMARKS
!
! 12. STRUCTURE
!
! 13. SOURCE TEXT
!
   SAVE  IENT
   DATA  IENT /0/
   CALL  STRACE (IENT,'STPNOW')
!
   IF(LEVERR >= 4)THEN
     STPNOW = .TRUE.
   ELSE
     STPNOW = .FALSE.
   END IF
   IF(MAXERR == -1) STPNOW = .FALSE.
!
   RETURN
   END FUNCTION STPNOW
 
!*****************************************************************
!                                                                *
   SUBROUTINE FOR (IUNIT, DDNAME, SF, IOSTAT)
!                                                                *
!*****************************************************************
!
   USE OCPCOMM1                                                        
   USE OCPCOMM2                                                        
   USE OCPCOMM3                                                        
   USE OCPCOMM4                                                        
!
   IMPLICIT NONE
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering                              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmers: R.C. Ris, N. Booij,                          |
!     |              IJ.G. Haagsma, A.T.M.M. Kieftenburg,         |
!     |              M. Zijlema, E.E. Kriezi,                     |
!     |              R. Padilla-Hernandez, L.H. Holthuijsen       |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 2004-2005  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. Authors
!
!     30.13: Nico Booij
!     30.70: Nico Booij
!     30.82: IJsbrand Haagsma
!     34.01: IJsbrand Haagsma
!     40.00, 40.03: Nico Booij
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     30.13, Jan. 96: new structure
!     30.70, Feb. 98: terminating error if input file does not exist
!     30.82, Nov. 98: Introduced recordlength of 1000 for new files to
!                     avoid errors on the Cray-J90
!     34.01, Feb. 99: STOP statement removed
!     40.00, Feb. 99: DIRCH2 replaces DIRCH1 in filenames
!     40.03, May  00: modification for Linux: local copy of filename
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  1. PURPOSE
!
!     General open file routine.
!
!  2. METHOD
!
!     FORTRAN 77 OPEN option.
!                INQUIRE
!
!  3. METHOD
!
!  4. ARGUMENT VARIABLES
!
!       IUNIT   int     input   =0 : get free unit number
!                               >0 : fixed unit number
!                       output  allocated unit number
!       DDNAME  char    input   ddname/filename string (empty if IUNIT>0)
!       SF      char*2  input   file qualifiers
!                               1st char: O(ld),N(ew),S(cratch),U(nknown)
!                               2nd char: F(ormatted),U(nformatted)
!       IOSTAT  int     input   0 : Full messages printed
!                               -1: Only error messages printed
!                               -2: No messages printed
!                       output  error indicator
!
   INTEGER   IUNIT, IOSTAT
   CHARACTER DDNAME*(LENFNM), SF*2                                     
!
!  5. PARAMETER VAR. (CONSTANTS)
!
!     Error codes:
!
!       IOSTAT = IESUCC No errors
!       IOSTAT > 0      I/O error
!       IOSTAT = IENUNF No free unit number found
!       IOSTAT = IEUNBD Specified unit number out of bounds
!       IOSTAT = IENODD No filename supplied with IUNIT=0
!       IOSTAT = IEDDNM Incorrect filename supplied with IUNIT>0
!       IOSTAT = IEEXST Specified unit number does not exist
!       IOSTAT = IEOPEN Specified unit number already opened
!       IOSTAT = IESTAT Error in file qualifiers
!       IOSTAT = IENSCR Named scratch file
!       IOSTAT = IENSIO No specified I/O error
!
   INTEGER  IESUCC, IENUNF, IEUNBD, IENODD,                     &
            IEDDNM, IEEXST, IEOPEN, IESTAT, IENSCR
   PARAMETER (IESUCC=  0,IENUNF= -1,IEUNBD= -2,IENODD= -3,      &
              IEDDNM= -4,IEEXST= -5,IEOPEN= -6,IESTAT= -7,      &
	      IENSCR=-12)
!
!  EMPTY    blank string
!
   CHARACTER  EMPTY*(*)
   PARAMETER (EMPTY= '        ')
!
!  6. LOCAL VARIABLES
!
!     IENT      number of entries into this subroutine
!     IFO       format index
!     IFUN      free unit number
!     II        counter
!     IOSTTM    aux. error index
!     IS        file status index
!     IUTTM     aux. unit number
!
   INTEGER   IENT, IFO, IFUN, II, IOSTTM, IS, IUTTM
!
!     EXIST     if true, file exists
!     OPENED    if true, file is opened
!
   LOGICAL   EXIST, OPENED
!
!     S
!     F
!     FILTTM   auxiliary
!     FISTAT   file status, values: OLD, NEW, UNKNOWN
!     FORM     formatting, values: FORMATTED, UNFORMATTED
!     DDNAME_L local copy of DDNAME                                       
!
   CHARACTER S, F, FILTTM *(LENFNM), DDNAME_L *(LENFNM)                
   CHARACTER *11 FISTAT(4),FORM(2)
!
!  4. SUBROUTINES USED
!
!
!  5. ERROR MESSAGES
!
!       and error messages added using MSGERR
!
!
!  6. REMARKS
!
!       Free unit number search interval: FUNLO<=IUNIT<=FUNHI
!       FUNLO, FUNHI, IUNMIN and IUNMAX were initialized by OCPINI,
!       they are transmitted via module OCPCOMM4
!
!  7. STRUCTURE
!
!       ----------------------------------------------------------------
!       Check file qualifiers
!       ----------------------------------------------------------------
!       If IUNIT = 0
!       Then If DDNAME = ' '
!            Then error message
!            Else Inquire to find if file exists and is opened,
!                 and if so, to find correct unit number
!                 If file is not opened
!                 Then get a free unit number, assign value to IUNIT
!                      open the file
!                 Else assign correct unit number to IUNIT
!       Else Inquire to find if file exists and is opened,
!                   and if so, to find correct filename
!            If file with unit nr IUNIT is already open
!            Then If filename does not correspond to DDNAME
!                 Then Close file with old filename and unit IUNIT
!                      Open file with new filename DDNAME and unit IUNIT
!            Else If DDNAME is not empty
!                 Then Open file with new filename DDNAME and unit IUNIT
!                 Else Open file with unit IUNIT
!       ----------------------------------------------------------------
!
!  8. SOURCE TEXT
!
   SAVE      IENT, IFUN
!
   DATA FISTAT(1),FISTAT(2) / 'OLD','NEW'/                            &
        FISTAT(3),FISTAT(4) / 'SCRATCH','UNKNOWN'/                    &
	FORM(1),FORM(2) / 'FORMATTED','UNFORMATTED'/
!
   DATA IENT /0/, IFUN /0/
   CALL STRACE (IENT, 'FOR')
!
   IF(ITEST >= 80) WRITE (PRTEST, 2) IUNIT, DDNAME, SF, IOSTAT
2  FORMAT (' Entry FOR: ', I3, 1X, A36, A2, I7)
   DDNAME_L = DDNAME                                                   
!
!     check file qualifiers
!
   IF((IUNIT /= 0) .AND. ((IUNIT < IUNMIN) .OR. (IUNIT > IUNMAX)))THEN
     IF(IOSTAT > -2) CALL MSGERR (3, 'Unit number out of range')
     IOSTAT= IEUNBD
     RETURN
   END IF
!
   S   = SF(1:1)
   F   = SF(2:2)
   IS  = INDEX('ONSU',S)
   IFO = INDEX('FU',F)
   IF((IS == 0) .OR. (IFO == 0))THEN
     IF(IOSTAT > -2) CALL MSGERR (3,'Error in file qualifiers')
     IOSTAT= IESTAT
     RETURN
   END IF
!
   IF((S == 'S') .AND. (DDNAME /= EMPTY))THEN
     IF(IOSTAT > -2) CALL MSGERR (3, 'Named scratch file')
     IOSTAT= IENSCR
     RETURN
   END IF
!
   IF(DDNAME /= EMPTY)THEN                                           
!       directory separation character is replaced in filenames           
     DO II = 1, LEN(DDNAME)
       IF(DDNAME(II:II) == DIRCH1) DDNAME(II:II) = DIRCH2             
     ENDDO
   ENDIF
!
   IF(IUNIT == 0)THEN
     IF(DDNAME == EMPTY)THEN
       IF(IOSTAT > -1) CALL MSGERR (3, 'No filename given')
       IOSTAT= IENODD
       RETURN
     ELSE
!         Was the file opened already ?
       INQUIRE (FILE=DDNAME, IOSTAT=IOSTTM, EXIST=EXIST,           &
	        OPENED=OPENED, NUMBER=IUTTM)
       IF(IOSTTM /= IESUCC)THEN
         IF(IOSTAT > -1)                                           &
	   CALL MSGERR (2,'Inquire failed, filename: '//DDNAME_L)          
         IOSTAT = IOSTTM
         RETURN
       ENDIF
!         If file does not exist, print term. error
       IF(IS == 1 .AND. .NOT. EXIST)THEN                           
         CALL MSGERR (4,'File cannot be opened/does not exist: '//DDNAME_L)
         IOSTAT = IEEXST
       END IF
       IF(OPENED)THEN
         IF(IOSTAT.GT.-1)                                          &
           CALL MSGERR (2, 'File is already opened: '//DDNAME_L)    
         IOSTAT = IEOPEN
         IUNIT = IUTTM
         RETURN
       ENDIF
!         Assign free unit number
       IF(IFUN == 0)THEN
         IFUN = FUNLO
       ELSE
         IFUN = IFUN + 1
       ENDIF
       IUNIT = IFUN
       IF(IUNIT > FUNHI)THEN
         IF(IOSTAT > -2) CALL MSGERR (3, 'All free units used')
         IOSTAT= IENUNF
       ENDIF
     END IF
     OPEN (UNIT=IUNIT,ERR=999,IOSTAT=IOSTTM,FILE=DDNAME,              &
!/Cray     RECL=1000,                                                 &
!/SGI      RECL=1000,                                                 &
!CVIS      SHARED,                                                    &
           STATUS=FISTAT(IS),ACCESS='SEQUENTIAL',FORM=FORM(IFO))      
   ELSE
     INQUIRE (UNIT=IUNIT, NAME=FILTTM, IOSTAT=IOSTTM,                 &
              EXIST=EXIST, OPENED=OPENED)
     IF(IOSTTM /= IESUCC)THEN
       IF(IOSTAT > -1)                                               &
         CALL MSGERR (2,'Inquire failed, filename: '//FILTTM)
       IOSTAT = IOSTTM
       RETURN
     ENDIF
     IF(OPENED)THEN
       IF(IOSTAT > -1)THEN
         CALL MSGERR (1,'File is already opened, filename: '//FILTTM)
       ENDIF
       IF(FILTTM /= DDNAME .AND. FILTTM /= EMPTY)THEN
         IF(IOSTAT > -2)THEN
           WRITE (PRINTF, '(A, I4, 6A)') ' unit', IUNIT,              &
	          ' filenames: ', FILTTM, ' and: ', DDNAME
           CALL MSGERR (2, 'filename and unit number inconsistent')
         ENDIF
         IOSTAT = IEDDNM
!          close old file and open new one with given filename
         CLOSE (IUNIT)
         OPEN (UNIT=IUNIT,ERR=999,IOSTAT=IOSTTM,STATUS=FISTAT(IS),    &
!/Cray         RECL=1000,                                             &
!/SGI          RECL=1000,                                             &
!CVIS          SHARED,                                                &
               FILE=DDNAME,ACCESS='SEQUENTIAL',FORM=FORM(IFO))
         IF(IOSTTM /= IESUCC) IOSTAT = IOSTTM
         GOTO 80
       ENDIF
       IOSTAT = IEOPEN
       RETURN
     END IF
     IF(DDNAME /= EMPTY)THEN
       OPEN (UNIT=IUNIT,ERR=999,IOSTAT=IOSTTM,STATUS=FISTAT(IS),      &
!/Cray       RECL=1000,                                               &         
!/SGI        RECL=1000,                                               &         
!CVIS        SHARED,                                                  & 
             FILE=DDNAME,ACCESS='SEQUENTIAL',FORM=FORM(IFO))
     ELSE
       OPEN (UNIT=IUNIT,ERR=999,IOSTAT=IOSTTM,STATUS=FISTAT(IS),      &
!/Cray       RECL=1000,                                               &
!/SGI        RECL=1000,                                               &
!CVIS        SHARED,                                                  &
             ACCESS='SEQUENTIAL',FORM=FORM(IFO))
     END IF
   END IF
   HIOPEN = IFUN
80 IF(ITEST >= 30) WRITE (PRINTF, 82) IUNIT, DDNAME, SF
82 FORMAT (' File opened: ', I6, 2X, A36, 2X, A2)
   RETURN
!
!  in case file cannot be opened:
!
999 IF(IOSTAT > -2)THEN
     CALL MSGERR (3, 'File open failed, filename: '//DDNAME_L)         
     WRITE (PRINTF,15) DDNAME, IOSTTM, SF
15   FORMAT (' File -> ', A36, 2X, ' IOSTAT=', I6, 4X, A2)
   ENDIF
   IUNIT = -1
   IOSTAT= IOSTTM

   RETURN
   END SUBROUTINE FOR
 
!***********************************************************************
!                                                                      *
    LOGICAL FUNCTION EQREAL (REAL1, REAL2 )                            
!                                                                      *
!***********************************************************************
!
    USE OCPCOMM1                                                        
    USE OCPCOMM2                                                        
    USE OCPCOMM3                                                        
    USE OCPCOMM4                                                        
!
    IMPLICIT NONE
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering                              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmers: R.C. Ris, N. Booij,                          |
!     |              IJ.G. Haagsma, A.T.M.M. Kieftenburg,         |
!     |              M. Zijlema, E.E. Kriezi,                     |
!     |              R. Padilla-Hernandez, L.H. Holthuijsen       |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 2004-2005  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. Authors
!
!     30.72 IJsbrand Haagsma
!     30.60 Nico Booij
!     40.04 Annette Kieftenburg
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     30.72, Oct. 97: Changed from EXCYES to make floating point point comparisons
!     30.60, July 97: new subroutine (EXCYES)
!     40.04, Aug. 00: introduced EPSILON and TINY
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     to determine whether a value (usually a value read from file)
!     is an exception value or not
!     Later (30.72) used to make comparisons of floating points within reasonable bounds
!
!  3. Method (updated...)
!
!     Checks whether ABS(REAL1-REAL2) .LE. TINY(REAL1) or whether this        40.04
!     difference is .LE. then EPS (= EPSILON(REAL1)*ABS(REAL1-REAL2) )        40.04
!
!  4. Argument variables
!
!     REAL1  : input    value that is to be tested
!     REAL2  : input    given exception value
!
    REAL      REAL1, REAL2
!
!  5. Parameter variables
!
!  6. Local variables
!
!     EPS    : Small number (related to REAL1 and its difference with REAL2)
!     IENT   : Number of entries into this subroutine
!
    REAL      EPS
    INTEGER   IENT
!
!  8. Subroutines used
!
!  9. Subroutines calling
!
!     SWREAD
!     SWDIM
!     SIRAY
!     SWBOUN
!     SWODDC
!     SWOEXD
!     SWOEXA
!     SWOEXF
!     SWPLOT
!     SWSPEC
!     ISOLIN
!     SNYPT2
!     INCTIM
!     INDBLE
!
! 10. Error messages
!
! 11. Remarks
!
! 12. Structure
!
! 13. Source text
!
    SAVE IENT
    DATA IENT/0/
    CALL STRACE(IENT,'EQREAL')
    EQREAL = .FALSE.
!
    EPS = EPSILON(REAL1)*ABS(REAL1-REAL2)                                  
    IF (EPS ==0) EPS = TINY(REAL1)                                         
    IF (ABS(REAL1-REAL2) .GT. TINY(REAL1)) THEN                            
      IF (ABS(REAL1-REAL2) .LT. EPS) EQREAL = .TRUE.                       
    ELSE                                                                   
      EQREAL = .TRUE.                                                      
    ENDIF                                                                  
    RETURN
    END FUNCTION EQREAL
!*****************************************************************
!                                                                *
      SUBROUTINE DTRETI (TSTRNG, IOPT, TIMESC)
!                                                                *
!*****************************************************************
!
      IMPLICIT NONE
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering                              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmers: R.C. Ris, N. Booij,                          |
!     |              IJ.G. Haagsma, A.T.M.M. Kieftenburg,         |
!     |              M. Zijlema, E.E. Kriezi,                     |
!     |              R. Padilla-Hernandez, L.H. Holthuijsen       |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 2004-2005  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. AUTHORS
!
!  1. UPDATES
!
!  2. PURPOSE
!
!  3. METHOD
!
!  4. ARGUMENT VARIABLES
!
!     IOPT   : input    option number
!
      INTEGER IOPT
!
!     TIMESC : output   time in seconds from given reference day REFDAY
!
      REAL    TIMESC
!
!     TSTRNG : input    time string
!
      CHARACTER  TSTRNG *(*)
!
!  5. PARAMETER VARIABLES
!
!  6. LOCAL VARIABLES
!
!     ITIME  : ??
!
      INTEGER ITIME(6)
!
!     DTTIME : Gives time in seconds from a reference day it also initialises the
!              reference day
!
      REAL    DTTIME
!
!  8. SUBROUTINE USED
!
!     DTSTTI   (installation dependent subroutines)
!
!  9. SUBROUTINES CALLING
!
! 10. ERROR MESSAGES
!
! 11. REMARKS
!
! 12. STRUCTURE
!
! 13. SOURCE TEXT
!
      CALL DTSTTI (IOPT, TSTRNG, ITIME)
      TIMESC = DTTIME (ITIME)
      RETURN
      END SUBROUTINE DTRETI
