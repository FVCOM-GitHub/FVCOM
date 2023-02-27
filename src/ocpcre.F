!               OCEAN PACK  command reading routines
!
!  Contents of this file:
!     RDINIT
!     NWLINE
!     INKEYW
!     INREAL
!     INDBLE
!     ININTG
!     INCSTR
!     INCTIM
!     ININTV
!     LEESEL
!     GETKAR
!     PUTKAR
!     UPCASE
!     KEYWIS
!     WRNKEY
!     IGNORE
!     RDHMS
!
!****************************************************************
!                                                               *
   SUBROUTINE RDINIT
!                                                               *
!****************************************************************
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
!     Initialises the command reading system
!
!  3. METHOD
!
!  4. ARGUMENT VARIABLES
!
!  5. PARAMETER VARIABLES
!
!  6. LOCAL VARIABLES
!
!     IENT   : Number of entries into this subroutine
!
   INTEGER   IENT
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
   SAVE IENT
   DATA IENT/0/
   CALL STRACE (IENT,'RDINIT')
   KAR = ';'
   KARNR = LINELN + 1                                                  
   ELTYPE = 'USED'
   BLANK = '    '
   RETURN
   END SUBROUTINE RDINIT
 
!****************************************************************
!                                                               *
   SUBROUTINE NWLINE
!                                                               *
!****************************************************************
!
   USE OCPCOMM1                                                        
   USE OCPCOMM2                                                        
!   USE OCPCOMM3                                                        
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
!  0. AUTHORS
!
!     34.01: IJsbrand Haagsma
!     40.03: Nico Booij
!     40.41: Marcel Zijlema
!
!  1. UPDATES
!
!     34.01, Feb. 99: Changed STOP statement in a MSGERR(4,'message')
!     40.03, Apr. 99: length of command lines changed from 80 to LINELN (=120)
!                     name of input file included in error message
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. PURPOSE
!
!     Jumps to reading of the next input line,
!     if the end of the previous one is reached.
!
!  3. METHOD
!
!  4. ARGUMENT VARIABLES
!
!  5. PARAMETER VARIABLES
!
!  6. LOCAL VARIABLES
!
!     IENT   : Number of entries into this subroutine
!
   INTEGER   IENT
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
   SAVE IENT
   DATA  IENT/0/
   CALL STRACE (IENT,'NWLINE')
5  IF((ELTYPE == 'USED').OR.(ELTYPE == 'EOR')) CALL LEESEL
!   print*,'after calling LEESEL. ***********',ELTYPE
   IF(ELTYPE == 'EOF') GOTO 90
   IF(ELTYPE == 'KEY' .AND. KEYWRD /= '        ') GOTO 50
   IF(ELTYPE == 'INT') GOTO 50
   IF(ELTYPE == 'REAL') GOTO 50
   IF(ELTYPE == 'CHAR') GOTO 50
   IF(KARNR <= LINELN) GOTO 50                                        
!  The end of the previous line is reached, there are no more
!  unprocessed data items on that line.
!  Jump to new line can take place.
   WRITE (PRINTF,9) '    '
9  FORMAT (A4)
   KARNR=0
   KAR=' '
   ELTYPE='USED'
   GOTO 5
90 IF(ITEST >= 10)THEN
     INQUIRE (UNIT=INPUTF, NAME=FILENM)                                
     WRITE (PRINTF, *) ' end of input file '//FILENM                   
   ENDIF
50 RETURN
   END SUBROUTINE NWLINE
 
!****************************************************************
!                                                               *
   SUBROUTINE INKEYW (KONT, CSTA)
!                                                               *
!****************************************************************
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
!  0. AUTHORS
!
!     40.41: Marcel Zijlema
!
!  1. UPDATES
!
!     ver 30.70, Jan. 1998: data type 'OTHR' is condidered
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. PURPOSE
!
!     this subroutine reads a keyword.
!
!  3. METHOD
!
!  4. ARGUMENT VARIABLES
!
!     KONT   : action to be taken if no keyword is found in input:
!              'REQ' (required) error message
!              'STA' (standard) the value of csta is assigned to keywrd.
!
!     CSTA   : see above.
!
   CHARACTER CSTA *(*), KONT *(*)
!
!  5. PARAMETER VARIABLES
!
!  6. LOCAL VARIABLES
!
!     IENT   : Number of entries into this subroutine
!     LENS   : length of default string (CSTA)
!
   INTEGER   IENT, LENS
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
   SAVE IENT
   DATA  IENT/0/
   CALL  STRACE ( IENT, 'INKEYW')
!
!     if necessary, a new data item is read.
!
   IF(ELTYPE == 'KEY' .AND. KEYWRD == '        ') GOTO 510
   IF(ELTYPE == 'KEY') GOTO 900
   IF(ELTYPE == 'EOR') GOTO 510
   IF(ELTYPE == 'USED') GOTO 510
   GOTO 520
510 CALL LEESEL
520 IF(ELTYPE == 'KEY') GOTO 900
!  KEYWORD IS READ
   IF((KONT == 'STA').OR.(KONT == 'NSKP'))THEN
     LENS = LEN(CSTA)
     IF(LENS >= 8)THEN
       KEYWRD = CSTA(1:8)
     ELSE
       KEYWRD = '        '
       KEYWRD(1:LENS) = CSTA
     ENDIF
     GOTO 900
   ENDIF
!  at the end of the input 'STOP' is generated.
   IF(ELTYPE == 'EOF')THEN
     KEYWRD='STOP'
     CALL MSGERR (2, 'STOP statement is missing')
     GOTO 900
   ENDIF
!  ----------------------------------------------------------
!  Data appear where a keyword is expected.
!  The user must be informed.
!  ----------------------------------------------------------
   IF(ELTYPE == 'EOR')THEN
     KEYWRD = '        '
     GOTO 900
   ENDIF
   IF(ELTYPE == 'INT')THEN
     CALL MSGERR (2, 'Data field skipped:'//ELTEXT)
     GOTO 510
   ENDIF
   IF(ELTYPE == 'REAL')THEN
     CALL MSGERR (2, 'Data field skipped:'//ELTEXT)
     GOTO 510
   ENDIF
   IF(ELTYPE == 'CHAR' .OR. ELTYPE == 'OTHR')THEN                    
     CALL MSGERR (2, 'Data field skipped:'//ELTEXT)
     GOTO 510
   ENDIF
   IF(ELTYPE == 'EMPT')THEN
     CALL MSGERR (2, 'Empty data field skipped')
     GOTO 510
   ENDIF
   CALL MSGERR (3, 'Error subr. INKEYW')
!  ----------------------------------------------------------
900 IF(ITEST >= 10) WRITE (PRINTF,910) KEYWRD
910 FORMAT (' KEYWORD: ',A8)
   RETURN
   END SUBROUTINE INKEYW
 
!****************************************************************
!                                                               *
   SUBROUTINE INREAL (NAAM, R, KONT, RSTA)
!                                                               *
!****************************************************************
!
   USE OCPCOMM1                                                        
   USE OCPCOMM2                                                        
   USE OCPCOMM3                                                        
   USE OCPCOMM4 
   USE MOD_UTILS                                                       
!
   IMPLICIT NONE
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 2004-2005  Delft University of Technology
!     FVCOM-SWAVE; a third generation wave model
!     Copyright (C) 2008-2009  University of Massachusetts Dartmouth
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
!     Jianhua Qi Dec. 20 2006
!
!  1. UPDATES
!
!  2. PURPOSE
!
!     Reads a REAL number in free format.
!
!  3. METHOD
!
!  4. ARGUMENT VARIABLES
!
!     R      : The value of the variable that is to be read.
!     RSTA   : Reference value needed for KONT='STA'or 'RQI'
!
   REAL      R, RSTA
!
!     KONT   : What to do with the varible?
!              ='REQ'; Variable is required
!              ='UNC'; If no variable, then variable will not be changed
!              ='STA'; If no variable, then variable will get value of RSTA
!              ='RQI'; Variable may not have the value of RSTA
!              ='REP' (REPEAT)
!              ='NSKP' (NO SKIP) IF DATA ITEM IS OF DIFFERENT TYPE,
!                VALUE IS LEFT UNCHANGED.
!     NAAM   : Name of the variable according to the user manual.
!
   CHARACTER NAAM *(*), KONT *(*)
!
!  5. PARAMETER VARIABLES
!
!  6. LOCAL VARIABLES
!
!     INPFIL
!     ISCAN
!     FNTMP
!
   CHARACTER(LEN=7) INPFIL
   INTEGER  ISCAN  
   REAL(SP)  FNTMP
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
!   INTEGER, SAVE IENT
!   DATA IENT /0/
!   CALL STRACE ( IENT, 'INREAL')
!
   INPFIL = "./INPUT"

   CHGVAL = .FALSE.
   ISCAN = SCAN_FILE2(INPFIL,NAAM,FSCAL = FNTMP)
   IF(ISCAN == 0)THEN
     IF(ABS(R-FNTMP) > 1.0E-6)CHGVAL = .TRUE.
     R = FNTMP
   ELSE
     IF(KONT == 'STA')THEN
       R = RSTA
     ELSE IF(KONT == 'REQ')THEN  
       WRITE(PRINTF,*)'ERROR READING ',NAAM,': ',ISCAN
       CALL PSTOP
     END IF
   END IF

!JQI   IF(KONT == 'STA')THEN
!JQI     R = RSTA
!JQI   ELSE  
!JQI     ISCAN = SCAN_FILE(INPFIL,NAAM,FSCAL = FNTMP)
!JQI     IF(ISCAN == 0)THEN
!JQI       IF(ABS(R-FNTMP) > 1.0E-6)CHGVAL = .TRUE.
!JQI       R = FNTMP
!JQI     ELSE
!JQI       IF(KONT /= 'UNC')THEN
!JQI         WRITE(PRINTF,*)'ERROR READING ',NAAM,': ',ISCAN
!JQI         CALL PSTOP
!JQI       END IF
!JQI     END IF  
!JQI   END IF        

   RETURN
   END SUBROUTINE INREAL
 
!****************************************************************
!                                                               *
   SUBROUTINE ININTG (NAAM, IV, KONT, ISTA)
!                                                               *
!****************************************************************
!
   USE OCPCOMM1                                                        
   USE OCPCOMM2                                                        
   USE OCPCOMM3                                                        
   USE OCPCOMM4                                                        
   USE MOD_UTILS                                                        
!
   IMPLICIT NONE
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 2004-2005  Delft University of Technology
!     FVCOM-SWAVE; a third generation wave model
!     Copyright (C) 2008-2009  University of Massachusetts Dartmouth
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
!     Jianhua Qi
!
!  1. UPDATES
!
!  2. PURPOSE
!
!     Reads an integer number, in free format
!
!  3. METHOD
!
!  4. ARGUMENT VARIABLES
!
!     IV     :  integer variable which is to be assigned a value
!     ISTA   :  default value
!
   INTEGER   IV, ISTA
!
!     NAAM   :  name of the variable according to the user manual
!     KONT   : What to do with the variable?
!              ='REQ'; error message if no value is found in the input file
!              ='UNC'; If no value, then variable will not be changed
!              ='STA'; If no value, then variable will get default value
!              ='RQI'; Variable may not have the value of RSTA
!              ='REP'  (repeat)
!              ='NSKP' (no skip) if data item is of different type,
!                      value is left unchanged.
!
   CHARACTER NAAM *(*), KONT *(*)
!
!  5. PARAMETER VARIABLES
!
!     PARAMETERS: SEE SUBR. INREAL
!
!  6. LOCAL VARIABLES
!
!     INPFIL
!     ISCAN
!     INTMP
!
   CHARACTER(LEN=7) INPFIL
   INTEGER  ISCAN  
   INTEGER  INTMP
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
!   SAVE IENT
!   DATA  IENT /0/
!   CALL  STRACE ( IENT, 'ININTG')
!
   INPFIL = "./INPUT"

!JQIJQI   IF(KONT == 'STA')THEN
!JQIJQI     IV = ISTA
!JQIJQI   ELSE  
!JQIJQI     ISCAN = SCAN_FILE(INPFIL,NAAM,ISCAL = INTMP)
!JQIJQI     IF(ISCAN == 0)THEN
!JQIJQI       IV = INTMP
!JQIJQI     ELSE
!JQIJQI       IF(KONT /= 'UNC')THEN
!JQIJQI         WRITE(PRINTF,*)'ERROR READING ',NAAM,': ',ISCAN
!JQIJQI         CALL PSTOP
!JQIJQI       END IF
!JQIJQI     END IF  
!JQIJQI   END IF        
      
   ISCAN = SCAN_FILE2(INPFIL,NAAM,ISCAL = INTMP)
   IF(ISCAN == 0)THEN
     IV = INTMP
   ELSE
     IF(KONT == 'STA')THEN
       IV = ISTA
     ELSE IF(KONT == 'REQ')THEN    
       WRITE(PRINTF,*)'ERROR READING ',NAAM,': ',ISCAN
       CALL PSTOP
     END IF
   END IF    
      
   RETURN
   END SUBROUTINE ININTG
 
!****************************************************************
!                                                               *
   SUBROUTINE INLOGC (NAAM, L, KONT, LSTA)
!                                                               *
!****************************************************************
!
   USE OCPCOMM1                                                        
   USE OCPCOMM2                                                        
   USE OCPCOMM3                                                        
   USE OCPCOMM4 
   USE MOD_UTILS                                                       
!
   IMPLICIT NONE
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 2004-2005  Delft University of Technology
!     FVCOM-SWAVE; a third generation wave model
!     Copyright (C) 2008-2009  University of Massachusetts Dartmouth
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
!     Jianhua Qi Dec. 20 2006
!
!  1. UPDATES
!
!  2. PURPOSE
!
!     Reads a REAL number in free format.
!
!  3. METHOD
!
!  4. ARGUMENT VARIABLES
!
!     R      : The value of the variable that is to be read.
!     RSTA   : Reference value needed for KONT='STA'or 'RQI'
!
   LOGICAL      L, LSTA
!
!     KONT   : What to do with the varible?
!              ='REQ'; Variable is required
!              ='UNC'; If no variable, then variable will not be changed
!              ='STA'; If no variable, then variable will get value of RSTA
!              ='RQI'; Variable may not have the value of RSTA
!              ='REP' (REPEAT)
!              ='NSKP' (NO SKIP) IF DATA ITEM IS OF DIFFERENT TYPE,
!                VALUE IS LEFT UNCHANGED.
!     NAAM   : Name of the variable according to the user manual.
!
   CHARACTER NAAM *(*), KONT *(*)
!
!  5. PARAMETER VARIABLES
!
!  6. LOCAL VARIABLES
!
!     INPFIL
!     ISCAN
!     FNTMP
!
   CHARACTER(LEN=7) INPFIL
   INTEGER  ISCAN  
   LOGICAL LTMP
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
!   INTEGER, SAVE IENT
!   DATA IENT /0/
!   CALL STRACE ( IENT, 'INREAL')
!
   INPFIL = "./INPUT"
   ISCAN = SCAN_FILE2(INPFIL,NAAM,LVAL = LTMP)
   IF(ISCAN == 0)THEN
     L = LTMP
   ELSE
     IF(KONT == 'STA')THEN
       L = LSTA
     ELSE IF(KONT == 'REQ')THEN 
       WRITE(PRINTF,*)'ERROR READING ',NAAM,': ',ISCAN
       CALL PSTOP
     END IF
   END IF  

!JQI   IF(KONT == 'STA')THEN
!JQI     L = LSTA
!JQI   ELSE  
!JQI     ISCAN = SCAN_FILE(INPFIL,NAAM,LVAL = LTMP)
!JQI     IF(ISCAN == 0)THEN
!JQI       L = LTMP
!JQI     ELSE
!JQI       IF(KONT /= 'UNC')THEN
!JQI         WRITE(PRINTF,*)'ERROR READING ',NAAM,': ',ISCAN
!JQI         CALL PSTOP
!JQI       END IF
!JQI     END IF  
!JQI   END IF        

   RETURN
   END SUBROUTINE INLOGC
!****************************************************************
!                                                               *
   SUBROUTINE INCSTR (NAAM, C, KONT, CSTA)
!                                                               *
!****************************************************************
!
   USE OCPCOMM1                                                        
   USE OCPCOMM2                                                        
   USE OCPCOMM3                                                        
   USE OCPCOMM4 
   USE MOD_UTILS                                                       
!
   IMPLICIT NONE
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 2004-2005  Delft University of Technology
!     FVCOM-SWAVE; a third generation wave model
!     Copyright (C) 2008-2009  University of Massachusetts Dartmouth
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
!     Jianhua Qi
!
!  1. UPDATES
!
!  2. PURPOSE
!
!     Reads a string in free format
!
!  3. METHOD
!
!  4. ARGUMENT VARIABLES
!
!     NAAM   : name of the variable according to the user manual
!     KONT   : What to do with the variable?
!              ='REQ'; error message if no value is found in the input file
!              ='UNC'; If no value, then variable will not be changed
!              ='STA'; If no value, then variable will get default value
!              ='RQI'; Variable may not have the value of CSTA
!              ='REP'  (repeat)
!              ='NSKP' (no skip) if data item is of different type,
!                      value is left unchanged.
!     C      : string that is to be read from input file
!     CSTA   : default value of the string
!
   CHARACTER NAAM *(*), KONT *(*), C *(*), CSTA *(*)
!
!  5. PARAMETER VARIABLES
!
!  6. LOCAL VARIABLES
!
!     INPFIL
!     ISCAN
!     INTMP
!
   CHARACTER(LEN=7) INPFIL
   INTEGER  ISCAN  
   CHARACTER(LEN=80) CHTMP
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
!   SAVE IENT
!   DATA  IENT /0/
!   CALL  STRACE ( IENT, 'INCSTR')
!
   INPFIL = "./INPUT"
   ISCAN = SCAN_FILE2(INPFIL,NAAM,CVAL = CHTMP)
   IF(ISCAN == 0)THEN
     C = TRIM(CHTMP)
   ELSE
     IF(KONT == 'STA')THEN
       C = CSTA
     ELSE IF(KONT == 'REQ')THEN  
       WRITE(PRINTF,*)'ERROR READING ',NAAM,': ',ISCAN
       CALL PSTOP
     END IF
   END IF  

!JQI   IF(KONT == 'STA')THEN
!JQI     C = CSTA
!JQI   ELSE  
!JQI     ISCAN = SCAN_FILE(INPFIL,NAAM,CVAL = CHTMP)
!JQI     IF(ISCAN == 0)THEN
!JQI       C = TRIM(CHTMP)
!JQI     ELSE
!JQI       IF(KONT /= 'UNC')THEN
!JQI         WRITE(PRINTF,*)'ERROR READING ',NAAM,': ',ISCAN
!JQI         CALL PSTOP
!JQI       END IF
!JQI     END IF  
!JQI   END IF        
      
   RETURN
   END SUBROUTINE INCSTR
 
!****************************************************************
!                                                               *
   SUBROUTINE INCTIM (IOPTIM, NAAM, RV, KONT, RSTA)
!                                                               *
!****************************************************************
!
   USE OCPCOMM1                                                        
   USE OCPCOMM2                                                        
   USE OCPCOMM3                                                        
   USE OCPCOMM4                                                        
   USE MOD_UTILS                                                       
!
   IMPLICIT NONE
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 2004-2005  Delft University of Technology
!     FVCOM-SWAVE; a third generation wave model
!     Copyright (C) 2008-2009  University of Massachusetts Dartmouth
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
!     Jianhua Qi
!
!  1. UPDATES
!
!  2. PURPOSE
!
!     Reads and interprets a time string
!
!  3. METHOD
!
!  4. ARGUMENT VARIABLES
!
!     IOPTIM   int   inp   time reading option (see subr DTSTTI)
!
!
   INTEGER   IOPTIM
!
!     RV     : variable that is to be assigned a value
!     RSTA   : default value
!
   REAL      RV, RSTA
!
!     NAAM   : name of the variable according to the user manual
!     KONT   : What to do with the variable?
!              ='REQ'; error message if no value is found in the input file
!              ='UNC'; If no value, then variable will not be changed
!              ='STA'; If no value, then variable will get default value
!              ='RQI'; Variable may not have the value of RSTA
!              ='REP'  (repeat)
!              ='NSKP' (no skip) if data item is of different type,
!                      value is left unchanged.
!
   CHARACTER NAAM *(*), KONT *(*)
!
!  5. PARAMETER VARIABLES
!
!     PARAMETERS: SEE PROGRAM DOCUMENTATION.
!
!  6. LOCAL VARIABLES
!
!     IENT   : Number of entries into this subroutine
!     LENMN  : length of the string NAAM
!
!     INPFIL
!     ISCAN
!     INTMP
!
   CHARACTER(LEN=7) INPFIL
   INTEGER  ISCAN  
   REAL FNTMP
   CHARACTER(LEN=80) CHTMP
   CHARACTER(LEN=24) C
   REAL R
   INTEGER    IENT, LENNM,INTR
!
!     NAAM_L : local copy of NAAM                                         
!
   CHARACTER (LEN=40) :: NAAM_L                                        
!
!     EQREAL : logical function, True if arguments are equal
!
   LOGICAL    EQREAL                                                   
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
   SAVE IENT
   DATA IENT /0/
   CALL STRACE ( IENT, 'INCTIM')
!
   INPFIL = "./INPUT"
   ISCAN = SCAN_FILE2(INPFIL,NAAM,CVAL = CHTMP)
   IF(ISCAN == 0)THEN
     C = TRIM(CHTMP)
   ELSE
!JQI     IF(KONT == 'STA')THEN
!JQI       RV = RSTA
!JQI     ELSE IF(KONT == 'REQ')THEN  
       WRITE(PRINTF,*)'ERROR READING ',NAAM,': ',ISCAN
       CALL PSTOP
!JQI     END IF
   END IF  

!JQI   IF(KONT == 'STA')THEN
!JQI     RV = RSTA
!JQI   ELSE  
!JQI!     ISCAN = SCAN_FILE(INPFIL,NAAM,FSCAL = FNTMP)
!JQI     ISCAN = SCAN_FILE(INPFIL,NAAM,CVAL = CHTMP)
!JQI     IF(ISCAN == 0)THEN
!JQI!       R = FNTMP
!JQI!       INTR = INT(R)
!JQI!       WRITE(C,'(F24.6)') R     !R
!JQI       C = TRIM(CHTMP)
!JQI!       print*,'INTR=',INTR,R,'C=',c
!JQI!       print*,'C=',c
!JQI     ELSE
!JQI       IF(KONT /= 'UNC')THEN
!JQI         WRITE(PRINTF,*)'ERROR READING ',NAAM,': ',ISCAN
!JQI         CALL PSTOP
!JQI       END IF
!JQI     END IF  
!JQI   END IF        

!   CALL DTRETI (ELTEXT(1:LENCST), IOPTIM, RV)
   CALL DTRETI (C, IOPTIM, RV)
!JQI     IF (.NOT.EQREAL(RV,RSTA)) CHGVAL = .TRUE.                        

   RETURN
   END SUBROUTINE INCTIM
 
!*******************************************************************
!                                                                  *
   SUBROUTINE ININTV (NAME, RVAR, KONT, RSTA)                          
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
!  0. AUTHORS
!
!     40.41: Marcel Zijlema
!
!  1. UPDATES
!
!     Dec 1995, ver 30.09 : new subroutine
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. PURPOSE
!
!     Read a time interval in the form: number  DAY/HR/MIN/SEC
!
!  3. METHOD
!
!  4. ARGUMENT VARIABLES
!
!     NAAM   : name of the variable according to the user manual
!     KONT   : What to do with the variable?
!              ='REQ'; error message if no value is found in the input file
!              ='UNC'; If no value, then variable will not be changed
!              ='STA'; If no value, then variable will get default value
!              ='RQI'; Variable may not have the value of RSTA
!              ='REP'  (repeat)
!              ='NSKP' (no skip) if data item is of different type,
!                      value is left unchanged.
!
   CHARACTER NAME *(*), KONT *(*)
!
!     RSTA   : default value
!     RVAR   : variable that is to be assigned a value
!
   REAL      RSTA, RVAR
!
!  5. PARAMETER VARIABLES
!
!  6. LOCAL VARIABLES
!
!     IENT   : Number of entries into this subroutine
!
   INTEGER   IENT
!
!     FAC    : a factor, value depends on unit of time used
!     RI     : auxiliary variable
!
   REAL      FAC, RI
!
!     KEYWIS : logical function, True if keyword encountered is equal to
!              keyword in user manual
!
   LOGICAL   KEYWIS
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
!     -------------------------------------------------------------
!     Call INREAL to read number of time units
!     If a value was read
!     Then Read time unit
!          Case time unit is
!          DAY: Fac = 24*3600
!          HR:  Fac = 3600
!          MI:  Fac = 60
!          SEC: Fac = 1
!     Else Fac = 1
!     -------------------------------------------------------------
!     Interval in seconds = Fac * number of time units
!     -------------------------------------------------------------
!
! 13. SOURCE TEXT
!
   SAVE IENT
   DATA IENT /0/
   CALL STRACE (IENT, 'ININTV')                                        
!
   CALL INREAL (NAME, RI, KONT, RSTA)
   IF(CHGVAL)THEN
     CALL INKEYW ('STA', 'S')
     IF(KEYWIS('DA'))THEN
       FAC = 24.*3600.
     ELSE IF(KEYWIS('HR'))THEN
       FAC = 3600.
     ELSE IF(KEYWIS('MI'))THEN
       FAC = 60.
     ELSE
       CALL IGNORE ('S')
       FAC = 1.
     ENDIF
   ELSE
     FAC = 1.
   ENDIF
   RVAR = FAC * RI
   RETURN
   END SUBROUTINE ININTV
 
!****************************************************************
!                                                               *
   SUBROUTINE LEESEL
!                                                               *
!****************************************************************
!
   USE OCPCOMM1                                                        
   USE OCPCOMM2                                                        
!   USE OCPCOMM3                                                        
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
!  0. AUTHORS
!
!     40.13: Nico Booij
!     40.41: Marcel Zijlema
!
!  1. UPDATES
!
!     Jan. 1994, mod. 20.05: ELREAL is made double precision
!     40.13, Jan. 01: ! is now added as comment sign
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. PURPOSE
!
!     reads a new data item from the string 'KAART'.
!     type of the item is determined, and the contents appears
!     in ELTEXT, ELINT, or ELREAL, as the case may be.
!     the following types are distinguished:
!     'KEY'   keyword
!     'INT'   integer or real number
!     'REAL'  real number
!     'CHAR'  character string enclosed in quotes
!     'EMPT'  empty data field
!     'OTHR'  non-empty data item not recognized as real, int or char,
!             possibly a time string
!     'EOF'   end of input file
!
!     'EOR'   end of repeat, or end of record
!     'ERR'   error
!     'USED'  used, item last read is processed already.
!
!  3. METHOD
!
!     difference between comment signs $ and !:                           40.13
!     everything on an input line behind a ! is ignored
!     text between two $-signs (on one line) is intepreted as comment
!     text behind two $-signs is intepreted as valid input
!
!  4. ARGUMENT VARIABLES
!
!  5. PARAMETER VARIABLES
!
!  6. LOCAL VARIABLES
!
!     IENT   : Number of entries into this subroutine
!     IRK    : auxiliary value used to detect errors
!     ISIGN1 : sign of mantissa part
!     ISIGN2 : sign of exponent part
!     ISTATE : state of the number reading process
!     J      : counter
!     JJ     : counter
!     JKAR   : counts the number of characters in the data field
!     NREP   : repetition number
!     NUM1   : value of integer part of mantissa
!     NUM2   : exponent value
!
   INTEGER   IENT, IRK, ISIGN1, ISIGN2, ISTATE, J, JJ, JKAR, NREP,     &
             NUM1, NUM2
!
!     RMANT  : real mantissa value
!
   DOUBLE PRECISION RMANT
!
!     QUOTE  : the quote character
!
   CHARACTER QUOTE *1
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
   SAVE  IENT, QUOTE, NREP                                             
   DATA  QUOTE/''''/ , IENT/0/, NREP/1/
   CALL  STRACE ( IENT, 'LEESEL')
!
   IF(NREP > 1)THEN
     NREP = NREP - 1
     GOTO 190
   ENDIF
!
!  initialisations
!
2  NREP = 1
   DO J=1,LINELN,4                                                  
     ELTEXT(J:J+3) = '    '
   END DO  
   JKAR = 1
   ELINT=0
   ELREAL=0.
!
!  start processing data item
!
   IF(KARNR == 0) GOTO 12
!  process a new character
10 IF(KAR == '!' .OR. KARNR > LINELN)THEN                           
!    end of the line is reached, if repetition factor is >1
!    the data item is assumed to be empty
     IF(NREP > 1) GOTO 28
!    end of the line is reached, if no repetition factor appears
!    the data item is assumed to be of type 'EOR'
     ELTYPE='EOR'
     IF(KAR == '!') KARNR = LINELN+1                                  
     GOTO 190
   ENDIF
!  skip leading blanks or Tab characters
11 IF(KAR /= ' ' .AND. KAR /= TABC) GOTO 20
  
!      print*,'before getkar 1'
12 CALL GETKAR
!      print*,'after getkar 1'
!     end of input file was reached
   IF(ELTYPE == 'EOF')THEN
!       generate keyword STOP
     ELTEXT='STOP'
     GOTO 190
   ENDIF
   GOTO 10
!  if character is comma, empty data field
20 IF(KAR /= ',') GOTO 30
!      print*,'before getkar 2'
   CALL GETKAR
!      print*,'after getkar 2'
28 ELTYPE='EMPT'
   GOTO 190
!  Notice: jump to label 28 (empty data field)
!  if after repetition a comment, a keyword, end of record etc. is found.
!  --------------------------------------------------------
!  see whether end of repeat (; or /) is marked
30 IF(INDEX(';/',KAR) > 0)THEN                                      
     IF(NREP > 1) GOTO 28
     ELTYPE='EOR'
!      print*,'before getkar 3'
     CALL GETKAR
!      print*,'after getkar 3'
     GOTO 190
   ENDIF
!  ( marks the beginning of a data item group; is ignored
38 IF(KAR == '(') GOTO 12
!  --------------------------------------------------------
!  comment; data enclosed in comment identifiers is interpreted as comment
!   print*,'KAR=',KAR,'COMID=',COMID,NREP
40 IF(KAR == COMID)THEN
     IF(NREP > 1) GOTO 28
!      print*,'before getkar 4'
41   CALL  GETKAR
!      print*,'after getkar 4'
     IF(KARNR > LINELN) GOTO 10                                      
     IF(KAR /= COMID) GOTO 41
     GOTO 12
   ENDIF
!  -------------------------------------------------------
!  if item is a number, read this integer or real number
!
!  integer number:  SIGN1]NUM1
!  real:            SIGN1]NUM1].]MANT]E]SIGN2]NUM2
!     ISTATE =    10     9    8 7    6 5     4    3
!  SIGN1, SIGN2:     + OR -
!  NUM1, NUM2, MANT: series of digits
!     -------------------------------------------------------
50 IF(INDEX('+-.0123456789',KAR) == 0) GOTO 80
   NUM1=0
   NUM2=0
   ISIGN1=1
   ISIGN2=1
   ISTATE=10
   IRK=0
   RMANT=0.
   ELTYPE='INT'
   IF(INDEX('+-',KAR) == 0) GOTO 52
   ISTATE=9
   IF(KAR == '-') ISIGN1=-1
   CALL PUTKAR (ELTEXT, KAR, JKAR)
!      print*,'before getkar 5'
   CALL GETKAR
!      print*,'after getkar 5'
!     ****  part before decimal point  ****
52 IF(INDEX('0123456789',KAR) == 0) GOTO 54
   IRK=1
   ISTATE=8
   NUM1=10*NUM1+INDEX('123456789',KAR)
   CALL PUTKAR (ELTEXT, KAR, JKAR)
!      print*,'before getkar 6'
   CALL GETKAR
!      print*,'after getkar 6'
   GOTO 52
54 IF(KAR /= '.') GOTO 56
   ISTATE=7
   ELTYPE='REAL'
   CALL PUTKAR (ELTEXT, KAR, JKAR)
!      print*,'before getkar 7'
   CALL GETKAR
!      print*,'after getkar 7'
56 JJ=-1
!     ****  part after decimal point  ****
57 IF(INDEX('0123456789',KAR) == 0) GOTO 58
   IRK=1
   ISTATE=6
   RMANT = RMANT + DBLE(INDEX('123456789',KAR))*1.D1**JJ               
   JJ=JJ-1
   CALL PUTKAR (ELTEXT, KAR, JKAR)
!      print*,'before getkar 8'
   CALL GETKAR
!      print*,'after getkar 8'
   GOTO 57
58 IF(ISTATE >= 9 .OR. IRK == 0) GOTO 120
!     ****  exponent part  ****
   IF(INDEX('DdEe^',KAR) == 0) GOTO 66
   ISTATE=5
   IRK=0
   IF(ELTYPE == 'INT') ELTYPE='REAL'
   CALL PUTKAR (ELTEXT, KAR, JKAR)
!      print*,'before getkar 9'
   CALL GETKAR
!      print*,'after getkar 9'
   IF(INDEX('+-',KAR) == 0) GOTO 62
   IF(KAR == '-') ISIGN2=-1
   ISTATE=4
   CALL PUTKAR (ELTEXT, KAR, JKAR)
!      print*,'before getkar 10'
   CALL GETKAR
!      print*,'after getkar 10'
62 IF(INDEX('0123456789',KAR) == 0) GOTO 66
   IRK=1
   ISTATE=3
   NUM2=10*NUM2+INDEX('123456789',KAR)
   CALL PUTKAR (ELTEXT, KAR, JKAR)
!      print*,'before getkar 11'
   CALL GETKAR
!      print*,'after getkar 11'
   GOTO 62
!     ****  a number is put together  ****
66 IF(IRK == 0) GOTO 120
   IF(INDEX('+-.',KAR) >= 1) ELTYPE='OTHR'
   ISTATE=2
   IF(ITEST >= 330) WRITE (PRINTF,699) ELTYPE, ISIGN1, NUM1,            &
                                       RMANT, ISIGN2, NUM2
699 FORMAT (1X, A4, 2I6, F12.9, 2I6)
   IF(ELTYPE == 'REAL') ELREAL =                                        &
             ISIGN1*(DBLE(NUM1)+RMANT) * 1.D1**(ISIGN2*NUM2)            
   IF(ELTYPE == 'INT') ELINT = ISIGN1*NUM1
   LENCST = JKAR - 1                                                   
!     skip trailing blanks
67 IF(KAR /= ' ' .AND. KAR /= TABC) GOTO 68
   ISTATE=1
!      print*,'before getkar 12'
   CALL GETKAR
!      print*,'after getkar 12'
   GOTO 67
!     If a * is encountered now, it is interpreted as a repetition factor.
68 IF(KAR == '*')THEN
     IF(ELTYPE == 'INT' .AND. ELINT > 0)THEN
       NREP = ELINT
       ELINT = 0
!      print*,'before getkar 13'
       CALL GETKAR
!      print*,'after getkar 13'
       GOTO 10
     ELSE
       CALL MSGERR (2, 'Wrong repetition factor')
!      print*,'before getkar 14'
       CALL GETKAR
!       print*,'after getkar 14'
       GOTO 190
     ENDIF
   ENDIF
69 IF(KAR == ',')THEN
!      print*,'before getkar 15'
     CALL GETKAR
!       print*,'after getkar 15'
     GOTO 190
   ENDIF
   IF(ISTATE == 1) GOTO 190
   IF(INDEX(' ;',KAR) /= 0 .OR. KAR == TABC)THEN
     GOTO 190
   ENDIF
!     number is not followed by , blank or tab; type is made OTHR:
   GOTO 120
!     ----------------------------------------------------------
!     a character string is read; it start and ends with a quote
!     ----------------------------------------------------------
80 IF(KAR == QUOTE)THEN                                              
     ELTYPE='CHAR'
     LENCST = 0                                                          
     JJ=1
!      print*,'before getkar 16'
82   CALL GETKAR
!       print*,'after getkar 16'
!       end of the string: end of record or closing quote
     IF(KARNR > LINELN) GOTO 190                                     
     IF(KAR == QUOTE)THEN
!      print*,'before getkar 17'
       CALL GETKAR
!      print*,'after getkar 17'
!         new character is not a quote; end of the string
       IF(KAR /= QUOTE) GOTO 88
!         double quote is read as a single quote; continue
     ENDIF
!       put the character into ELTEXT
84   ELTEXT(JJ:JJ) = KAR
     LENCST = JJ
     JJ=JJ+1
     GOTO 82
!       process characters behind the string
!      print*,'before getkar 18'
87   CALL GETKAR
!      print*,'after getkar 18'
!       skip trailing blanks
88   IF(KAR == ' ' .OR. KAR == TABC) GOTO 87
     IF(KAR /= ',') GOTO 190
!      print*,'before getkar 19'
     CALL GETKAR
!      print*,'after getkar 19'
     GOTO 190
   ENDIF
!     -------------------------------------------------------
!     a keyword is read
!     a keyword starts with a letter (upper or lower case)
!     -------------------------------------------------------
90 CALL UPCASE (KAR)
!   print*,'KAR= ',KAR,'QUOTE= ',QUOTE
   IF(INDEX('ABCDEFGHIJKLMNOPQRSTUVWXYZ',KAR) > 0)THEN              
     IF(NREP > 1) GOTO 28
     ELTYPE='KEY'
     ISTATE=2
     JJ=1
92   ELTEXT(JJ:JJ) = KAR
     LENCST = JJ                                                       
!      print*,'before getkar 20'
     CALL GETKAR
!      print*,'after getkar 20'
     CALL UPCASE (KAR)
     JJ=JJ+1
!       next characters: letters, digits or - _ .
     IF(INDEX('ABCDEFGHIJKLMNOPQRSTUVWXYZ',KAR) >= 1) GOTO 92
     IF(INDEX('0123456789-_.',KAR) >= 1) GOTO 92
!       keyword is read
     KEYWRD = ELTEXT(1:8)
!     print*,'KEYWRD= ',KEYWRD
!       trailing blanks or tab char are skipped
94   IF(KAR /= ' ' .AND. KAR /= TABC) GOTO  96
!      print*,'before getkar 21'
     CALL GETKAR
!      print*,'after getkar 21'
     GOTO 94
!       closure character  : or = is processed
96   IF(INDEX('=:',KAR) == 0) GOTO 190
!      print*,'before getkar 22'
     CALL GETKAR
!      print*,'after getkar 22'
     GOTO 190
   ENDIF
!     --------------------------------------------------
!     continuation symbol is read
!     --------------------------------------------------
100 IF(INDEX('_&',KAR) == 0) GOTO 120
   IF(NREP > 1) GOTO 28
110 KARNR=0
   GOTO 12
!     --------------------------------------------------
!     other type of data
!     --------------------------------------------------
120 ELTYPE='OTHR'
122 ELTEXT(JKAR:JKAR) = KAR                                             
   LENCST = JKAR                                                       
   JKAR=JKAR+1
!      print*,'before getkar 23'
   CALL GETKAR
!       print*,'after getkar 23'
   IF(INDEX(' ,;', KAR) >= 1 .OR. KAR == TABC) GOTO 126
   GOTO 122
!      print*,'before getkar 24'
126 CALL GETKAR
!       print*,'after getkar 24'
! 127  CALL MSGERR (3, 'Read error in: ')
!      WRITE (PRINTF,129) ELTEXT
! 129  FORMAT (A)                                                         
!      RETURN
!     --------------------------------------------------
!     test output and return to calling program
!     --------------------------------------------------
190 IF(ITEST >= 120) WRITE (PRTEST, 199) KAR, KARNR, ELTYPE, ELREAL,   &
                                         ELINT, NREP, ELTEXT(1:LENCST)
199 FORMAT (' test LEESEL: ', A1, 1X, I4, 1X, A4, D12.4, 2I6, 2X, A)    
   RETURN
   END SUBROUTINE LEESEL

!****************************************************************
!                                                               *
   SUBROUTINE GETKAR
!                                                               *
!****************************************************************
!
   USE OCPCOMM1                                                        
   USE OCPCOMM2                                                        
!   USE OCPCOMM3                                                        
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
!  0. AUTHORS
!
!     40.13: Nico Booij
!     40.41: Marcel Zijlema
!
!  1. UPDATES
!
!     40.13, Jan. 2001: TRIM used to limit output
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. PURPOSE
!
!     This procedure reads a next character (KAR) from the string KAART.
!     The position of this character in KAART is indicated by KARNR.
!     If needed, a new input line is read.
!     At the end of the input file ELTYPE is made 'EOF'.
!
!  3. METHOD
!
!  4. ARGUMENT VARIABLES
!
!  5. PARAMETER VARIABLES
!
!  6. LOCAL VARIABLES
!
!     IENT   : Number of entries into this subroutine
!
   INTEGER   IENT
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
   SAVE IENT
   DATA  IENT /0/
   CALL STRACE (IENT, 'GETKAR')
   
   IF(KARNR == 0)THEN
     READ(INPUTF, '(A)', END=20) KAART
     
!     print*,'KAART=',trim(KAART),lineln
     IF(ITEST >= -10) WRITE (PRINTF, '(1X,A)') TRIM(KAART)                   
     KARNR=1
   ENDIF
   IF(KARNR > LINELN)THEN                                           
     KAR=';'
     GOTO 90
   ENDIF
   KAR = KAART(KARNR:KARNR)
!   print*,'KAR=',trim(KAR)
   KARNR=KARNR+1
   GOTO 90
!     end of file is encountered
20 ELTYPE='EOF'
   KAR='@'
90 IF(ITEST >= 320) WRITE (PRINTF, '(" Test GETKAR", 2X, A4, 2X, A1, I4)')  &
                    ELTYPE, KAR, KARNR

   RETURN
   END SUBROUTINE GETKAR
 
!****************************************************************
!                                                               *
   SUBROUTINE PUTKAR (LTEXT, KARR, JKAR)
!                                                               *
!****************************************************************
!
   USE OCPCOMM1                                                        
   USE OCPCOMM2                                                        
!   USE OCPCOMM3                                                        
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
!     this procedure inserts a character (KARR) usually read by GETKAR
!     into the string LTEXT, usually equal to ELTEXT, in the place
!     JKAR. After this JKAR is increased by 1.
!
!  3. METHOD
!
!  4. ARGUMENT VARIABLES
!
!     JKAR   : counts the number of characters in a data field
!
   INTEGER  JKAR
!
!     LTEXT  : a character string; after a number of calls it should
!              contain the character representation of a data field
!     KARR   : character to be inserted into LTEXT
!
   CHARACTER LTEXT *(*), KARR *1
!
!  5. PARAMETER VARIABLES
!
!  6. LOCAL VARIABLES
!
!     IENT   : Number of entries into this subroutine
!
   INTEGER   IENT
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
   SAVE IENT
   DATA IENT /0/
   CALL STRACE (IENT, 'PUTKAR')
   
   IF(JKAR > LEN(LTEXT)) CALL MSGERR (2, 'PUTKAR, string too long')
   LTEXT(JKAR:JKAR) = KARR
   LENCST = JKAR
   JKAR = JKAR + 1
   RETURN
   END SUBROUTINE PUTKAR

!****************************************************************
!                                                               *
   SUBROUTINE UPCASE (CHARST)
!                                                               *
!****************************************************************
!
   USE OCPCOMM1                                                        
   USE OCPCOMM2                                                        
!   USE OCPCOMM3                                                        
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
!     changes all characters of the string CHARST from lower to
!     upper case
!
!  3. METHOD
!
!  4. ARGUMENT VARIABLES
!
!     CHARST : a character string
!
   CHARACTER*(*) CHARST                                                
!
!  5. PARAMETER VARIABLES
!
!  6. LOCAL VARIABLES
!
!     IC     : sequence number of a character in the string CHARST
!     IENT   : Number of entries into this subroutine
!     KK     : position of a character in a given string
!     LLCC   : length of the given character string
!
   INTEGER   IC, IENT, KK, LLCC
!
!     ABCUP  : A to Z upper case characters
!     ABCLO  : a to z lower case characters
!     CC     : a character
!
   CHARACTER ABCUP *26, ABCLO *26, CC *1
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
   SAVE IENT
   DATA  IENT /0/
   DATA  ABCUP /'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
   DATA  ABCLO /'abcdefghijklmnopqrstuvwxyz'/
   CALL STRACE (IENT, 'UPCASE')
!
   LLCC = LEN(CHARST)
   DO IC = 1, LLCC
     CC = CHARST(IC:IC)
     KK = INDEX(ABCLO, CC)
     IF(KK /= 0) CHARST(IC:IC) = ABCUP(KK:KK)
   END DO  
   RETURN
   END SUBROUTINE UPCASE
 
!****************************************************************
!                                                               *
   LOGICAL FUNCTION KEYWIS (STRING)
!                                                               *
!****************************************************************
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
!  0. AUTHORS
!
!     40.41: Marcel Zijlema
!
!  1. UPDATES
!
!     40.00, July
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. PURPOSE
!
!     This procedure tests whether a keyword given by the user
!     coincides with a keyword known in the program (i.e. string).
!     if so, keywis is made .True., otherwise it is .False.
!     also ELTYPE is made 'USED', so that next element can be read.
!
!  3. METHOD
!
!  4. ARGUMENT VARIABLES
!
!     STRING : a keyword which is compared with a keyword found in the input file
!
   CHARACTER STRING *(*)
!
!  5. PARAMETER VARIABLES
!
!  6. LOCAL VARIABLES
!
!     IENT   : Number of entries into this subroutine
!     J      : counter
!     LENSS  : length of the keyword STRING
!
   INTEGER   IENT, J, LENSS
!
!     KAR1   : a character of the keyword appearing in the input file
!     KAR2   : corresponding character in the STRING
!
   CHARACTER KAR1 *1, KAR2 *1
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
   SAVE IENT
   DATA IENT /0/
   CALL STRACE (IENT, 'KEYWIS')
!
   KEYWIS = .FALSE.
   IF(ELTYPE == 'USED') GOTO 30
!
   KEYWIS=.TRUE.
   LENSS = LEN (STRING)
   DO  J=1, LENSS
     KAR1 = KEYWRD(J:J)
     KAR2 = STRING(J:J)
     IF(KAR1 /= KAR2 .AND. KAR2 /= ' ')THEN                          
       KEYWIS=.FALSE.
       GOTO 30
     ENDIF
   END DO
   IF(ELTYPE == 'KEY') ELTYPE = 'USED'
30 RETURN
   END FUNCTION KEYWIS

!****************************************************************
!                                                               *
   SUBROUTINE  WRNKEY
!                                                               *
!****************************************************************
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
!     THIS PROCEDURE PRODUCES AN ERROR MESSAGE
!     IT IS CALLED IF AN ILLEGAL KEYWORD IS FOUND IN THE
!     USER'S INPUT. IT MAKES ELTYPE = 'USED'
!
!  3. METHOD
!
!  4. ARGUMENT VARIABLES
!
!  5. PARAMETER VARIABLES
!
!  6. LOCAL VARIABLES
!
!     IENT   : Number of entries into this subroutine
!
   INTEGER   IENT
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
   SAVE IENT
   DATA IENT /0/
   CALL STRACE (IENT, 'WRNKEY')
!
   CALL MSGERR (2, 'Illegal keyword: '//KEYWRD)
   ELTYPE = 'USED'
   RETURN
   END SUBROUTINE WRNKEY
 
!****************************************************************
!                                                               *
   SUBROUTINE IGNORE (STRING)
!                                                               *
!****************************************************************
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
!     This procedure calls subroutine INKEYW to read a keyword.
!     if this keyword is equal to string, eltype is made 'USED'.
!     it is used if a keyword can occur in the input which
!     does not lead to any action.
!
!  3. METHOD
!
!  4. ARGUMENT VARIABLES
!
!     STRING : keyword (if appearing in input file) that can be ignored
!
   CHARACTER STRING *(*)
!
!  5. PARAMETER VARIABLES
!
!  6. LOCAL VARIABLES
!
!     IENT   : Number of entries into this subroutine
!
   INTEGER   IENT
!
!     KEYWIS : logical function
!
   LOGICAL   KEYWIS
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
   SAVE IENT
   DATA IENT /0/
   CALL STRACE (IENT, 'IGNORE')
!
   CALL INKEYW ('STA', 'XXXX')
   IF(KEYWIS(STRING)) RETURN
   IF(KEYWIS('XXXX')) RETURN
   IF(ITEST >= 60) WRITE (PRINTF, 5) KEYWRD, ELTYPE
5  FORMAT (' NOT IGNORED: ', A, 2X, A)
   RETURN
   END SUBROUTINE IGNORE
 
