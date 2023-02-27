!     Last change:  NB   04 Jan 2001   12:00 pm
!
!          OCEAN PACK - Installation dependent subroutines
!
!*****************************************************************
!                                                                *
   SUBROUTINE OCPINI (INIFIL, LREAD, INERR)
!                                                                *
!*****************************************************************
!
   USE OCPCOMM1                                                        
   USE OCPCOMM2                                                        
   USE OCPCOMM3                                                        
   USE OCPCOMM4                                                        
   USE M_PARALL                                                        

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
!     30.82: IJsbrand Haagsma
!     34.01: IJsbrand Haagsma
!     40.00, 40.03: Nico Booij
!     40.30: Marcel Zijlema
!     40.31: Marcel Zijlema
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     10.02, July 94: New argument INIFIL
!                     Check on validity period now uses OCDTIM
!     30.74, Nov. 97: Prepared for version with INCLUDE statements
!     30.82, Nov. 98: Introduced recordlength of 1000 for file PRINT to
!                     avoid error-messages on the Cray-J90 and SGI Origin 200
!     34.01, Feb. 99: Changed STOP statements for MSGERR(4,'message')
!                     calls
!     34.01, Feb. 99: Opens a file 'screen' when unitnr in swaninit<>6
!     40.00, Feb. 99: Directory separation characters included in init file
!                     these characters are used in subr FOR
!     40.03, May  00: backslash replaced by CHAR(92) because of problems on Linux
!     40.30, Jan. 03: introduction distributed-memory approach using MPI
!     40.31, Nov. 03: removing HPGL-functionality
!     40.41, Sep. 04: includes speed processors in initialisation file
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     subroutine initialises a number of common variables
!     opens standard input and output files, if necessary
!
!  4. Argument variables
!
!     INERR : output     Number of the initialisation error
!
   INTEGER INERR
!
!     INIFIL  inp  char  name of initialisation file
!     LREAD   inp  log   if True: command input file must be opened
!                        and command reading must be initialised
!
   LOGICAL    LREAD, FILEXI
   CHARACTER  INPFIL *40, OUTFIL *40, INIFIL *(*), TSTFIL *40,       &
              PLTOPT *4, TIMSTR *24,                                 &
	      OUTFO *40, TSTFO *40, TXT*120
   INTEGER    PRCTIM(6), INIVER, INIVEF
   INTEGER    PFROPT                                                   
   REAL       PLPARM(10)                                               
   INTEGER    NUMM(10)                                                 
   LOGICAL    STPNOW
   DATA PRCTIM /0,0,0,0,0,0/
!
!  version of initialisation file                                      
   INIVER = 4                                                         
   INIVEF = -1                                                         
   INERR  = 0                                                          
!
!  see whether initialisation file exists
!
   INQUIRE (FILE=INIFIL, EXIST=FILEXI)
   IF(FILEXI)THEN
!
!    read initialisation file
!
     OPEN (11, FILE=INIFIL, STATUS='OLD',                         &     
!CVIS             SHARED,                                        &     
                ERR=950)                                              
     READ (11, *, ERR=930, END=930) INIVEF 
     IF(INIVEF > INIVER .OR. INIVEF <= 0) GOTO 935                   
     READ (11, 120, ERR=930, END=930) INST
     READ (11, *,   ERR=930, END=930) INPUTF
     READ (11, 120, ERR=930, END=930) INPFIL
     READ (11, *,   ERR=930, END=930) PRINTF
     READ (11, 120, ERR=930, END=930) OUTFIL
     READ (11, *,   ERR=930, END=930) PRTEST
     READ (11, 120, ERR=930, END=930) TSTFIL
     READ (11, *,   ERR=930, END=930) SCREEN
     READ (11, *,   ERR=930, END=930) IUNMAX
     READ (11, 130, ERR=930, END=930) COMID
     READ (11, 130, ERR=930, END=930) TABC
     IF(INIVEF >= 2)THEN                                             
       READ (11, 130, ERR=930, END=930) DIRCH1
       READ (11, 130, ERR=930, END=930) DIRCH2
     ELSE
!DOS          DIRCH1 = CHAR(47)                                               
!DOS          DIRCH2 = CHAR(92)                                               
!UNIX          DIRCH1 = CHAR(92)                                               
!UNIX          DIRCH2 = CHAR(47)                                               
          DIRCH1 = CHAR(92)                                               
          DIRCH2 = CHAR(47)                                               
     ENDIF
     IF(INIVEF < 3)THEN
       READ (11, 140, ERR=930, END=930) PLTOPT
       READ (11, *,   ERR=930, END=930) NPLP
       READ (11, *,   ERR=930, END=930) (PLPARM(II),II=1,NPLP)
       READ (11, *,   ERR=930, END=930) PFROPT
     END IF
     READ (11, *,   ERR=930, END=930) ITMOPT
!JQI     IF(INIVEF > 3)THEN
!JQI       IF(PARLL)THEN
!JQI         DO JJ = 1, NPROC
!JQI           READ (11, 145, ERR=150, END=150) IW, TXT
!JQI           CALL TXPBLA(TXT,IF,IL)
!JQI           IPOS = 0
!JQI           DO II = IF, IL
!JQI             K = ICHAR(TXT(II:II))
!JQI             IF(K >= 48 .AND. K <= 57)THEN
!JQI               IPOS = IPOS + 1
!JQI               NUMM(IPOS) = K - 48
!JQI             END IF
!JQI           END DO
!JQI           J = 0
!JQI           DO II = 1, IPOS
!JQI             J = J + NUMM(II)*10**(IPOS-II)
!JQI           END DO
!JQI           IWEIG(J) = IW
!JQI         END DO
!JQI       END IF
!JQI     END IF
120  FORMAT (A40)
130  FORMAT (A1)
140  FORMAT (A4)
145  FORMAT (I5,A)
150  CLOSE (11)
   ELSE
!
!    REFERENCE NUMBERS AND NAMES OF STANDARD FILES
!
     INPUTF = 3
     INPFIL = 'INPUT'
     PRINTF = 4
     OUTFIL = 'PRINT'
!    unit ref. numbers for output to screen and to separate
!    test print file:
     PRTEST = PRINTF
     TSTFIL = '    '
     SCREEN = 6
     IUNMAX = 99
!    TABC is the Tab character (interpreted as blank in command reading)
!DOS        TABC   = '	'
!UNIX        TABC   = '	'
!    COMID is the comment identifier (usually $)
     COMID  = '$'
!    DIRCH1 is directory separation character as appears in input file
!    DIRCH2 is directory separation character replacing DIRCH1         
!DOS        DIRCH1 =  CHAR(47)                                                
!DOS        DIRCH2 =  CHAR(92)                                                
!UNIX        DIRCH1 =  CHAR(92)                                                
!UNIX        DIRCH2 =  CHAR(47)                                                
!    INST = name of institute, max. 40 characters
     INST = 'Delft University of Technology'
     ITMOPT = 1
!
   ENDIF
   OUTFO = OUTFIL                                                      
   TSTFO = TSTFIL                                                      

!  --- append node number to OUTFIL and TSTFIL                         
!      in case of parallel computing                                   

!JQI   IF(PARLL)THEN                                                     
!JQI     ILPOS = INDEX ( OUTFIL, ' ' )-1
!JQI     IF(ILPOS > 0)THEN
!JQI       WRITE(OUTFIL(ILPOS+1:ILPOS+4),180) INODE
!JQI     ELSE
!JQI       GOTO 920
!JQI     END IF
!JQI     ILPOS = INDEX ( TSTFIL, ' ' )-1
!JQI     IF(ILPOS > 0)THEN
!JQI       WRITE(TSTFIL(ILPOS+1:ILPOS+4),180) INODE
!JQI     END IF
!JQI180  FORMAT('-',I3.3)
!JQI     CALL SWSYNC
!JQI     IF(STPNOW()) RETURN
!JQI   END IF

   IF(INIVEF < INIVER)THEN
!
!    write initialisation file
!
     OPEN  (12, FILE=INIFIL, STATUS='UNKNOWN', FORM='FORMATTED', ERR=950)
     WRITE (12, 210) INIVER, 'version of initialisation file'          
     WRITE (12, 220) INST,   'name of institute'
     WRITE (12, 210) INPUTF, 'command file ref. number'
     WRITE (12, 220) INPFIL, 'command file name'
     WRITE (12, 210) PRINTF, 'print file ref. number'
     WRITE (12, 220) OUTFO,  'print file name'                         
     WRITE (12, 210) PRTEST, 'test file ref. number'
     WRITE (12, 220) TSTFO,  'test file name'                          
     WRITE (12, 210) SCREEN, 'screen ref. number'
     WRITE (12, 210) IUNMAX, 'highest file ref. number'
     WRITE (12, 230) COMID,  'comment identifier'
     WRITE (12, 230) TABC,   'TAB character'
     WRITE (12, 230) DIRCH1, 'dir sep char in input file'
     WRITE (12, 230) DIRCH2, 'dir sep char replacing previous one'     
     WRITE (12, 210) ITMOPT, 'default time coding option'
!JQI     IF(PARLL)THEN
!JQI       DO II = 1, NPROC
!JQI         WRITE (12, 240) IWEIG(II), 'speed of processor ',II
!JQI       END DO
!JQI     END IF
     CLOSE (12)
210  FORMAT (I5, T41, A)
220  FORMAT (A40, A)
230  FORMAT (A1, T41, A)
240  FORMAT (I5, T41, A19, I3)
   ENDIF
!
   IUNMIN = 0
   FUNLO = 21
   FUNHI = IUNMAX
!
   CALL OCDTIM (PRCTIM)
!
!  initialise command reader
!
   IF(OUTFIL /= '    ')THEN
!    WRITE (*,*) ' Open print file ', PRINTF, OUTFIL
     OPEN (UNIT=PRINTF, FILE=OUTFIL, STATUS='UNKNOWN', FORM='FORMATTED', &
!/Cray         RECL=2000,                                                &
!/SGI          RECL=2000,                                                & 
           ERR=920)                                     
!    WRITE (*,*) ' Print file opened ', PRINTF, OUTFIL
     CALL DTTIST (ITMOPT, TIMSTR, PRCTIM)
     WRITE (PRINTF, 12) TIMSTR
12   FORMAT ('1',//,20X, 'Execution started at ',A, //)
   ENDIF
   IF(PRTEST /= PRINTF) OPEN (UNIT=PRTEST, FILE=TSTFIL, ERR=922)
   IF(SCREEN /= 6) OPEN(UNIT=SCREEN, FILE='screen',ERR=960)           
   IF(LREAD)THEN
     IF (INPFIL /= '    ') OPEN (UNIT=INPUTF, FILE=INPFIL, STATUS='OLD', &
!CVIS           SHARED,                                                  &
                ERR=910)                                               
     CALL RDINIT
     IF (INPFIL /= '    ') CLOSE (UNIT=INPUTF) !conflict with open statement in function SCAN_FILE2, so close it.
   ENDIF
!
   RETURN
!
910 CALL MSGERR(4,'Input file missing')                                 
   RETURN                                                              
!
920 INERR=920                                                           
   IF (INODE == MASTER ) WRITE(*,*) 'Cannot open PRINT file '     
   RETURN                                                             
!
922 CALL MSGERR(4,'Cannot open test file: '//TSTFIL)                    
   RETURN
!
930 INERR=930                                                           
   IF(INODE == MASTER) WRITE(*,*) 'Error reading initialisation file '  
   RETURN                                                              
!
935 INERR=935                                                           
   IF(INODE == MASTER ) WRITE(*,*) 'Incorrect version of initialisation file '
   RETURN                                                              
!
950 INERR=950                                                           
   IF(INODE == MASTER) WRITE(*,*) 'Error opening initialisation file ' 
   RETURN                                                              
!
960 CALL MSGERR(4,'Error opening output file: screen')                  
   RETURN                                                              
!
   END SUBROUTINE OCPINI
 
!*****************************************************************
!                                                                *
   SUBROUTINE OCDTIM (PRCTIM)
!                                                                *
!*****************************************************************
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
!     30.07
!     30.70: Nico Booij
!     30.82: IJsbrand Haagsma
!     40.02: IJsbrand Haagsma
!
!  1. Updates
!
!     30.07, Oct. 95: option DEC added
!     30.70, Sep. 97: adaptation in view of year 2000
!     30.82, Mar. 99: Adapted to Fortran 90 standard
!     40.02, Sep. 00: Removed all platform dependent Fortran 77 statements
!
!  2. PURPOSE
!
!       get time of processing, using processor dependent routines
!
!  3. PARAMETER LIST
!
!       PRCTIM  outp   int   time array: elements: year, month, day,
!                            hour, minute, second
!
!  4. SUBROUTINES USED
!
!       GETDAT, GETTIM or other
!
!  5. ERROR MESSAGES
!
!       ----
!
!  6. REMARKS
!
!       This function uses a processor dependent subroutines (GETDAT,
!       GETTIM), therefore adaptations are necessary when compiled
!       at a different computer system environment.
!
!  7. STRUCTURE
!
!       ---------------------------------------------------------
!       Call DATE and TIME routines (system dependent)
!       decode YEAR, MONTH and DAY
!       assemble DATE string
!       decode HOUR, MINUTE, SECOND
!       assemble TIME string
!       ---------------------------------------------------------
!
!  8. SOURCE TEXT
!
   INTEGER PRCTIM(6)
!
!     Call DATE and TIME routines
!
!     --------Fortran 90 date-time routines --------
!
   CHARACTER TIMSTR *24, CDUMMY *5                                     
   INTEGER   IDUMMY(8)                                                 
!
   CALL DATE_AND_TIME (TIMSTR(1:8), TIMSTR(10:20), CDUMMY, IDUMMY)     
   CALL DTSTTI (1, TIMSTR, PRCTIM)                                     

   RETURN

   END SUBROUTINE OCDTIM
 
!*****************************************************************
!                                                                *
   SUBROUTINE DTSTTI (IOPT, TIMSTR, DTTIME)
!                                                                *
!*****************************************************************
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
!     Updates
!
!       ver 30.70, Sep 1997 by N.Booij: adaptation in view year 2000
!
!     Function:
!
!       transform time string into integer time array
!
!     Argument list:
!
!       IOPT    input  int   option number
!                            1: ISO notation   19870530.153000
!                            2: (HP compiler): 30-May-87 15:30:00
!                            3: (old Lahey)    05/30/87 15:30:00          
!                            4:                         15:30:00
!                            5:                87/05/30 15:30:00
!                            6: WAM            8705301530
!
!       TIMSTR  input  char  time string
!       DTTIME  outp   int   time array: elements: year, month, day,
!                            hour, minute, second
!
!    Remarks:
!     Options can be added by the user
!     existing options should not be changed
!
!     Source:
!
   INTEGER    IOPT, DTTIME(6)
   CHARACTER  TIMSTR *24, MONC(12) *3, MONCI *3
   DATA MONC /'JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL',       &
              'AUG', 'SEP', 'OCT', 'NOV', 'DEC'/
!
   IF(IOPT == 1)THEN
     READ (TIMSTR, '(I4,I2,I2,1X,3I2)', ERR=98) (DTTIME(II), II=1,6)
	
     print*,timstr,dttime,'ppppppppppppppppppppppp'
	
   ELSE IF(IOPT == 2)THEN
     READ (TIMSTR, '(I2,1X,A3,1X,I2,3(1X,I2))', ERR=98)              &
       DTTIME(3), MONCI, DTTIME(1), (DTTIME(II), II=4,6)
     IF(DTTIME(1) < 10)THEN                                         
       DTTIME(1) = 2000 + DTTIME(1)                                    
     ELSE
       DTTIME(1) = 1900 + DTTIME(1)
     ENDIF
     DTTIME(2) = 0
     DO 20 IMM = 1, 12                                                 
       CALL UPCASE (MONCI)
       IF (MONCI.NE.MONC(IMM)) GOTO 20
       DTTIME(2) = IMM
       GOTO 90
20   CONTINUE
     CALL MSGERR (2, 'incorrect month string: '//MONCI)                
   ELSE IF(IOPT == 3)THEN
     READ (TIMSTR, '(I2,5(1X,I2))', ERR=98)                          &
       DTTIME(2), DTTIME(3), DTTIME(1), (DTTIME(II), II=4,6)
     IF(DTTIME(1) < 10)THEN                                         
       DTTIME(1) = 2000 + DTTIME(1)                                    
     ELSE
       DTTIME(1) = 1900 + DTTIME(1)
     ENDIF
   ELSE IF(IOPT == 4)THEN
     READ (TIMSTR, '(I2,2(1X,I2))', ERR=98) (DTTIME(II), II=4,6)
     DO II = 1, 3
       DTTIME(II) = 0
     END DO
   ELSE IF(IOPT == 5)THEN                                            
     READ (TIMSTR, '(I2,5(1X,I2))', ERR=98) (DTTIME(II), II=1,6)       
     IF(DTTIME(1) < 10)THEN                                         
       DTTIME(1) = 2000 + DTTIME(1)                                    
     ELSE
       DTTIME(1) = 1900 + DTTIME(1)
     ENDIF
   ELSE IF(IOPT == 6)THEN                                            
     READ (TIMSTR, '(5I2)', ERR=98) (DTTIME(II), II=1,5)
     DTTIME(6) = 0.
     IF(DTTIME(1) < 10)THEN                                         
       DTTIME(1) = 2000 + DTTIME(1)                                    
     ELSE
       DTTIME(1) = 1900 + DTTIME(1)
     ENDIF
   ELSE
     CALL MSGERR (2, 'wrong time coding option in subroutine DTSTTI')
   ENDIF
90 RETURN
98 CALL MSGERR (2, 'time string unreadable: '//TIMSTR)
   RETURN
   END SUBROUTINE DTSTTI
 
!*****************************************************************
!                                                                *
   SUBROUTINE DTTIST (IOPT, TIMSTR, DTTIME)
!                                                                *
!*****************************************************************
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
!     Updates
!
!       ver 30.70, Sep 1997 by N.Booij: adaptation in view year 2000
!
!     Function:
!
!       transform integer time array into time string
!
!     Argument list:
!
!       IOPT    input  int   option number (see subr. DTSTTI)
!       TIMSTR  outp   char  time string
!       DTTIME  input  int   time array: elements: year, month, day,
!                            hour, minute, second
!
!     Source:
!
   INTEGER    IOPT, DTTIME(6)
   CHARACTER  TIMSTR *24, MONC(12) *3
   DATA MONC /'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul',      &
              'Aug', 'Sep', 'Oct', 'Nov', 'Dec'/
!
   TIMSTR = '    '
   IF(IOPT == 1)THEN
     WRITE (TIMSTR, 12) (DTTIME(II), II=1,6)
12   FORMAT (I4,I2,I2,'.',3I2)
     LTS = 15
   ELSE IF(IOPT == 2)THEN
     IF(DTTIME(1) >= 2000)THEN                                       
       DTTIME(1) = DTTIME(1) - 2000                                    
     ELSE
       DTTIME(1) = DTTIME(1) - 1900                                    
     ENDIF
     WRITE (TIMSTR, 22) DTTIME(3), MONC(DTTIME(2)), DTTIME(1),         &
                       (DTTIME(II), II=4,6)
22   FORMAT (I2,'-',A3,'-',I2,'.',I2,':',I2,':',I2)
     LTS = 18
   ELSE IF(IOPT == 3)THEN
     IF(DTTIME(1) >= 2000)THEN                                       
       DTTIME(1) = DTTIME(1) - 2000                                    
     ELSE
       DTTIME(1) = DTTIME(1) - 1900                                    
     ENDIF
     WRITE (TIMSTR, 32) DTTIME(2), DTTIME(3), DTTIME(1), (DTTIME(II), II=4,6)
32   FORMAT (I2,'/',I2,'/',I2,'.',I2,':',I2,':',I2)
     LTS = 17
   ELSE IF(IOPT == 4)THEN
     WRITE (TIMSTR, 42) (DTTIME(II), II=4,6)
42   FORMAT (I2,':',I2,':',I2)
     LTS = 8
   ELSE IF(IOPT == 5)THEN                                            
     IF(DTTIME(1) >= 2000)THEN                                       
       DTTIME(1) = DTTIME(1) - 2000                                    
     ELSE
       DTTIME(1) = DTTIME(1) - 1900                                    
     ENDIF
     WRITE (TIMSTR, 52) (DTTIME(II), II= 1,6)                          
52   FORMAT (I2,'/',I2,'/',I2,'.',I2,':',I2,':',I2)
     LTS = 17
   ELSE IF(IOPT == 6)THEN                                            
     IF(DTTIME(1) >= 2000)THEN                                       
       DTTIME(1) = DTTIME(1) - 2000                                    
     ELSE
       DTTIME(1) = DTTIME(1) - 1900                                    
     ENDIF
     WRITE (TIMSTR, 62) (DTTIME(II), II=1,5)                           
62   FORMAT (5I2)
     LTS = 10
   ELSE
     CALL MSGERR (2, 'wrong time coding option in subroutine DTTIST')
   ENDIF
!
   DO IC = 1, LTS
     IF (TIMSTR(IC:IC).EQ.' ') TIMSTR(IC:IC) = '0'
   END DO  
!
   RETURN
   END SUBROUTINE DTTIST
 
