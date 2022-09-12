!
!***********************************************************************
!                                                                      *
#  if defined (WAVE_CURRENT_INTERACTION)
   SUBROUTINE SWREAD !(COMPUT)                                          
!                                                                      *
!***********************************************************************
!
!     Reading and processing of the user commands describing the model
!
!***********************************************************************

   USE TIMECOMM                                                        
   USE OCPCOMM1                                                        
   USE OCPCOMM2                                                        
   USE OCPCOMM4                                                        
   USE SWCOMM1                                                         
   USE SWCOMM2                                                         
   USE SWCOMM3                                                         
   USE SWCOMM4                                                         
   USE OUTP_DATA                                                       
   USE M_SNL4                                                          !
   USE M_GENARR                                                        
   USE M_OBSTA                                                         
!   USE M_PARALL
   USE MOD_UTILS
!   USE MOD_USGRID
#  if defined (EXPLICIT)
   USE MOD_ACTION_EX
#  else      
   USE MOD_ACTION_IM
#  endif   
#  if defined (MULTIPROCESSOR)   
   USE MOD_PAR
#  endif  
   USE ALL_VARS
   USE VARS_WAVE
   USE MOD_OBCS
   USE MOD_NCLL
   USE MOD_NCTOOLS
   USE MOD_FORCE
   USE MOD_NESTING, ONLY : NESTING_ON_WAVE
#  if defined (VEGETATION)
   USE MOD_VEGETATION
#  endif
   IMPLICIT NONE
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 2008  Delft University of Technology
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
!  0. Authors
!     Jianhua Qi
!
!  1. Updates
!
!  2. Purpose
!
!     Reading and processing of the user commands describing the model
!
!  3. Method (updated 30.72)
!
!  4. Argument variables (updated 30.72)
!
!  5. Parameter variables
!
!  6. Local variables
!

   INTEGER, PARAMETER :: MXINCL = 10                                   
#  if defined (VEGETATION)
   REAL              :: HH, CC, DD
   INTEGER           :: NN
#  endif

   INTEGER, SAVE     :: INCNUM(1:MXINCL) = 0                           
   INTEGER, SAVE     :: INCLEV = 1                                     
!
   INTEGER           :: IERR = 0                                       
   INTEGER           :: ICNL4, ILAMBDA                                 
   INTEGER           :: ITMP1, ITMP2, ITMP3, ITMP4                     
!
   LOGICAL           :: MORE                                           
   LOGICAL, SAVE     :: LOBST = .FALSE.                                
!
   INTEGER        :: LREF, LREFDIFF, LRFRD                             
   REAL           :: KDIF                                              
   REAL           :: POWD                                              
   REAL           :: POWS, DUM                                         
   REAL           :: FD1, FD2, FD3, FD4                                
   REAL, ALLOCATABLE :: RLAMBDA(:)                                     

   TYPE(OBSTDAT), POINTER :: OBSTMP                                    
   TYPE(OBSTDAT), SAVE, POINTER :: COBST                               

   TYPE(OPSDAT), POINTER :: OPSTMP                                     

   TYPE XYPT                                                           
     REAL                :: X, Y
     TYPE(XYPT), POINTER :: NEXTXY
   END TYPE XYPT

   TYPE(XYPT), TARGET  :: FRST                                         
   TYPE(XYPT), POINTER :: CURR, TMP                                    

#  if defined (VEGETATION)
   TYPE VEGPT
     INTEGER              :: N
     REAL                 :: H, D, C
     TYPE(VEGPT), POINTER :: NEXTV
   END TYPE VEGPT

   TYPE(VEGPT), TARGET  :: FRSTV
   TYPE(VEGPT), POINTER :: CURRV, TMPV
#  endif
   LOGICAL :: KEYWIS                                                   
   LOGICAL :: EQREAL                                                   

   LOGICAL :: FOUND                                                    
   LOGICAL, SAVE :: RUNMADE = .FALSE.                                  

   LOGICAL, SAVE :: LOGCOM(1:6) = .FALSE.                              
   INTEGER   TIMARR(6)                                                 
!   CHARACTER PSNAME *8, PNAME *8, COMPUT *(*), PTYPE *1, DTTIWR *18    
   CHARACTER PSNAME *8, PNAME *8, PTYPE *1, DTTIWR *18    
   INTEGER, SAVE :: IENT = 0       ! number of entries to this subr
   INTEGER, SAVE :: LWINDR = 0     ! if non-zero, there is wind
   INTEGER, SAVE :: LWINDM = 3     ! type of wind growth formulation
   INTEGER, ALLOCATABLE :: IARR(:) 
   INTEGER :: GEN,IOS,LMAX                                    

   CHARACTER(LEN=7) INPFIL
   CHARACTER(LEN=80) CHTMP
   CHARACTER(LEN=80) PROJJ,SET,MODE,COORDINATES,BOUNDARY,CONSTANT,FRIC_FORM
   REAL FNTMP
   INTEGER INTMP
   LOGICAL LTMP,STAT,CART,CHECK
   INTEGER  INTVEC(150),ISCAN   !,KTEMP
   
   CHARACTER(LEN=80) INP_BOT_NAME,INP_CUR_NAME,INP_WI_NAME,      &
                     INP_FR_NAME,INP_WLEV_NAME,PROP_CHOICE
		     
   CHARACTER(LEN=20) CHARWAVEPERIOD
   CHARACTER(LEN=20) GROWTH
   CHARACTER(LEN=20) CNL4,OBST_TYPE,RDIFF,OFF_TYPE,NUM_CHOICE,STATL, &
                     SIGIM
   CHARACTER(LEN=20) WCAP		     
   CHARACTER(LEN=6)  UNIT
   LOGICAL QUAD,AGROW,BRE,MDIAL,FRICL,VEGE,OBST,TRI,REFL,RFD,SETUP,    &
           DIFFR,OFF,PROP,FLUXLIM,DIR,COM_STAT	
   LOGICAL INP_BOT,INP_CUR,INP_CUR_STAT,INP_CUR_SERI,INP_WI,INP_WI_STAT, &
           INP_WI_SERI,INP_FR,INP_FR_STAT,INP_FR_SERI,INP_WLEV,          &
	   INP_WLEV_STAT,INP_WLEV_SERI
   LOGICAL LAM,LIM,LINE,NUM	   
   	   
	   
!   REAL, ALLOCATABLE :: DEPTH_GL(:),UXB_GL(:,:),UYB_GL(:,:),             &
!           WXI_GL(:,:),WYI_GL(:,:),FRIC_GL(:,:),WLEVL_GL(:,:)
   REAL, ALLOCATABLE :: UXB_GL(:,:),UYB_GL(:,:),FRIC_GL(:,:),WLEVL_GL(:,:)
   REAL :: IFLTMP,NTIME_TMP,TEMPXY	   
   
   REAL :: ALTMP,FRLOW,FRHIG,GAMMA,TMPDIR,TRCF,HGT,OGAM,OBET,REF0,DIFF
   INTEGER :: MSS,IVAL
   INTEGER :: I,IGR2,IGRD,IGR1,J,JJ,NUMCOR,I1,NCNT
   INTEGER :: ITRAS	   	      
   INTEGER, ALLOCATABLE :: TEMP1(:),TEMP2(:)
   REAL    :: DEGCNV
!
!  ***** read command *****
!
!
!  ------------------------------------------------------------------
!
!  PROJECT      reading of project title and description
!
! ===============================================================
!
! PROJect  'NAME'  'NR'
!
!          'title1'
!
!          'title2'
!
!          'title3'
!
! ===============================================================
!

   CALL INCSTR('PROJECT',PROJJ,'REQ',' ')
   
   IF(PROJJ /= 'default')THEN
     CALL INCSTR('NAME',PROJID,'UNC',' ')
     CALL INCSTR('NR',PROJNR,'REQ',' ')
     CALL INCSTR('title1',PROJT1,'UNC',' ')   
     CALL INCSTR('title2',PROJT2,'UNC',' ')    
     CALL INCSTR('title3',PROJT3,'UNC',' ')
   END IF  
    
!
!  ------------------------------------------------------------------
!
!  SET       setting physical parameters and error counters
!
! =============================================================
!
!   SET  [level]  [nor]  [depmin]  [maxmes]        &
!        [maxerr]  [grav]  [rho]  [inrhog]         &
!        [hsrerr]  CARTesian/NAUTical  [pwtail]    &
!        [froudmax]  [printf]  [prtest]
!
! =============================================================
!
   CALL INCSTR('SET',SET,'REQ',' ')
   
   IF(SET /= 'default')THEN
     CALL INREAL('LEVEL', WLEV,  'UNC',0.0)
     CALL INREAL('NOR',   DNORTH,'UNC',0.0)
     CALL INREAL('DEPMIN',DEPMIN,'UNC',0.0)
     CALL ININTG('MAXMES',MAXMES,'UNC',0)
     CALL ININTG('MAXERR',MAXERR,'UNC',0)
     CALL INREAL('GRAV',  GRAV_W,  'UNC',0.0)
     CALL INREAL('RHO',   RHO_W,   'UNC',0.0)
     CALL ININTG('INRHOG',INRHOG,'UNC',0)
     IF(INRHOG == 0)THEN                                            
       OVUNIT(7) = 'm2/s'
       OVUNIT(9) = 'm2/s'
       OVUNIT(19) = 'm3/s'
       OVUNIT(21) = 'm2s'
       OVUNIT(22) = 'm2'
       OVUNIT(29) = 'm2'
     ELSE
       OVUNIT(7) = 'W/m2'
       OVUNIT(9) = 'W/m2'
       OVUNIT(19) = 'W/m'
       OVUNIT(21) = 'Js/m2'
       OVUNIT(22) = 'J/m2'
       OVUNIT(29) = 'J/m2'
     ENDIF                                                            
      
     CALL INREAL('HSRERR',HSRERR,'UNC',0.0)

     CALL INLOGC('NAUTICAL',BNAUT,'UNC','T')
     IF(BNAUT) OUTPAR(4) = 0.   
     IF(ITEST >= 20)THEN                                        
       WRITE(PRTEST,*) ' Set NAUT command: BNAUT=', BNAUT            
     END IF                                                            

     CALL INREAL('PWTAIL',PWTAIL(1),'UNC',0.0)
     IF(CHGVAL)THEN
       IF(PWTAIL(1) <= 1.) CALL MSGERR (3, 'Incorrect PWTAIL')
       PWTAIL(3) = PWTAIL(1) + 1.
     END IF
     
     CALL INREAL('FROUDMAX',PNUMS(18),'UNC',0.0)
     CALL ININTG('PRINTF',  PRINTF,   'UNC',0)
     CALL ININTG('PRTEST',  PRTEST,   'UNC',0)
   END IF 
   
!
!  ------------------------------------------------------------------
!
!  MODE  : Set STATionary, DYNamic (NONSTAtionary) or 1D SWAN model    
!
! ========================================
!
!        | -> STAtionary    |     | -> TWODimensional |                   
! MODE  <                    >   <                     >  (NOUPDATe)      
!        |    NONSTationary |     |    ONEDimensional |                   
!
! ========================================
!
   CALL INCSTR('MODE',MODE,'REQ',' ')

   IF(MODE /= 'default')THEN
     CALL INLOGC('STATIONARY',STAT,'UNC','T')

     IF(.NOT. STAT)THEN
       IF(NSTATM == 0) CALL MSGERR (2, 'Mode Nonst incorrect here')   
       NSTATM = 1                                                      
       NSTATC = 1                                                      
!
!      switch on flag for computation of default initial condition     
       ICOND = 1                                                       
!      no relaxation                                                   
       PNUMS(30) = 0.                                                  
     ELSE  
       IF(NSTATM == 1) CALL MSGERR (2, 'Mode STAT incorrect here')    
       NSTATM = 0                                                      
     ENDIF                                                             

     CALL INLOGC('ONED',ONED,'UNC','T')
   END IF

   CALL INLOGC('ACUPDAT',ACUPDA,'UNC','F')

!
!  ------------------------------------------------------------------
!
!  COORD : spherical or cartesian coordinates
!
! =============================================================
!
!   COORDinates  /  -> CARTesian               \   REPeating
!                \ SPHErical [rearth]   UM/QC  /
!
! =============================================================
!
   CALL INCSTR('COORDINATES',COORDINATES,'REQ',' ')
   
   IF(COORDINATES /= 'default')THEN
     CALL ININTG('KSPHER',KSPHER,'UNC',0)

     IF(KSPHER ==1)THEN
#    if defined (SPHERICAL)
       CALL INREAL ('REARTH', REARTH2, 'UNC', 0.)  !Jianzhong Ge                     
       LENDEG = REARTH2 * PI_W / 180.                                   
!      change properties of output quantities Xp and Yp
       OVUNIT(1) = 'degr'
       OVLLIM(1) = -200.
       OVULIM(1) =  400.
       OVLEXP(1) = -180.
       OVHEXP(1) =  360.
       OVEXCV(1) = -999.
       OVUNIT(2) = 'degr'
       OVLLIM(2) = -100.
       OVULIM(2) =  100.
       OVLEXP(2) = -90.
       OVHEXP(2) =  90.
       OVEXCV(2) = -999.

       CALL ININTG('PROJ_METHOD',PROJ_METHOD,'UNC',0)
#      endif
     ELSE IF(KSPHER /= 0)THEN
       WRITE(PRINTF,*)'PARAMETER KSPHER NOT EQUAL TO 0 OR 1: ',KSPHER
       CALL PSTOP
     END IF
   END IF
   
   CALL ININTG('KREPTX',KREPTX,'STA',0)
       
!
!     ------------------------------------------------------------------
!
!     CGRID     definition of computational grid
!
! ===========================================================================
!
!   CGRID  / REGular [xpc] [ypc] [alpc] [xlenc] [ylenc] [mxc] [myc] \
!          \ CURVilinear [mxc] [myc]  [excval]                      /    &
!          / CIRcle               \
!          \ SECtor [dir1] [dir2] /  [mdc]  [flow]  [fhig]  [msc]
!
! ===========================================================================
!
!   CALL INCSTR('GRIDFILE',GRIDFILE,'REQ',' ')
!   INQUIRE(FILE=TRIM(GRIDFILE),EXIST=CHECK) 
!   IF(CHECK)THEN
!     OPEN(1, FILE=TRIM(GRIDFILE),STATUS='OLD')
!   
!   ELSE
!     WRITE(PRINTF,*) TRIM(GRIDFILE), ' DOES NOT EXIST.'
!     PRINT*, TRIM(GRIDFILE), ' DOES NOT EXIST.'
!     CALL PSTOP
!   END IF
   
   LOGCOM(2) = .TRUE.   
   
   OPTG = 1

   CALL INREAL('ALPC',ALPC,'UNC',0.0)
!
!  ALPC is made to be between -PI and PI
!
   ALTMP = ALPC / 360.
   ALPC = PI2_W * (ALTMP - NINT(ALTMP))
   CVLEFT = .TRUE.                                                 

   CALL INLOGC('FULCIR',FULCIR,'UNC','T')
   IF(.NOT. FULCIR)THEN
     CALL INREAL('DIR1',SPDIR1,'REQ',0.0)
     CALL INREAL('DIR2',SPDIR2,'REQ',0.0)
   END IF  

   CALL ININTG('MDC',MDC,'REQ',0)
   CALL INREAL('FLOW',FRLOW,'REQ',0.0)
   CALL INREAL('FHIGH',FRHIG,'REQ',0.0)
   CALL ININTG('MSC',MSS,'REQ',0)
   IVAL = MAX(ABS(MSS*NINT(FRLOW)),ABS(MSS*NINT(FRHIG)),              &
              ABS(NINT(FRLOW)*NINT(FRHIG)))                          
   IF(IVAL >= 999**2)THEN                                          
     CALL MSGERR(3,'At least, FLOW and FHIGH or MSC must be given!')  
   END IF                                                            
   IVAL = MAX(ABS(NINT(FRLOW)),ABS(NINT(FRHIG)),ABS(MSS))            
   IF(IVAL /= 999)THEN                                             
     GAMMA = EXP(ALOG(FRHIG/FRLOW)/REAL(MSS)) - 1.                  
     WRITE (PRINTF,'(A,(F7.4))')                                      &
                  ' Resolution in sigma-space: df/f = ',GAMMA  
   ELSE IF(MSS == -999)THEN                                        
     MSS = NINT(ALOG(FRHIG/FRLOW)/ALOG(1.1))                        
     WRITE (PRINTF,'(A,I3)')                                          &
                  ' Number of meshes in sigma-space: MSC-1 = ',MSS  
   ELSE IF(FRLOW == -999.)THEN                                     
     FRLOW = FRHIG/EXP(ALOG(1.1)*REAL(MSS))                         
     WRITE (PRINTF,'(A,(F7.4))')                                      &
                  ' Lowest discrete frequency (in Hz): FLOW = ',FRLOW  
   ELSE IF(FRHIG == -999.)THEN                                     
     FRHIG = FRLOW*EXP(ALOG(1.1)*REAL(MSS))                         
     WRITE (PRINTF,'(A,(F7.4))')                                      &
                  ' Highest discrete frequency (in Hz): FHIGH = ',FRHIG  
   END IF                                                            
   SLOW = 2.*PI_W*FRLOW
   SHIG = 2.*PI_W*FRHIG
   MSC  = MSS+1
!
!    ***** MCGRD is the total number of nodes   *****
!    ***** MDC is the number of steps in Theta direction
!              as parts of a circle
!    ***** MSC is the number of frequencies in sigma direction (logaritmic)
     IF(FULCIR)THEN
       DDIR  = PI2_W / MDC                                               
!
!      modification of SPDIR1 first installed with version 30.72, then
!      reversed and reinstalled with 40.13
!      purpose: prevent problem with grids under 45 degrees            
       SPDIR1 = ALPC + 0.5 * DDIR                                      
     ELSE
       IF(BNAUT)THEN                                                 
!        swap values of SPDIR1 and SPDIR2, and transform               
         TMPDIR = SPDIR1                                               
         SPDIR1 = 180. + DNORTH - SPDIR2                               
         SPDIR2 = 180. + DNORTH - TMPDIR                               
       ENDIF                                                           
       SPDIR1 = SPDIR1 * PI_W / 180.                                     
       SPDIR2 = SPDIR2 * PI_W / 180.                                     
       IF(SPDIR2 < SPDIR1) SPDIR2 = SPDIR2 + PI2_W                     
       DDIR = (SPDIR2-SPDIR1) / REAL(MDC)                              
       MDC = MDC + 1                                                   
     ENDIF
!
     IF(.NOT.ALLOCATED(SPCSIG)) ALLOCATE(SPCSIG(MSC))                  
     IF(.NOT.ALLOCATED(SPCDIR)) ALLOCATE(SPCDIR(MDC,6))                
     CALL SSFILL(SPCSIG,SPCDIR)                                        
!
     IF(ITEST >= 20)THEN
       IF(OPTG == 1)WRITE (PRINTF,6048)                              
       IF(OPTG == 3)WRITE (PRINTF,6049)                              
       WRITE (PRINTF,6045) SLOW, SHIG, FRINTF
       WRITE (PRINTF,6046) MCGRD,NGL,MDC,MSC
       WRITE (PRINTF,6047) DDIR
!       WRITE (PRINTF,6047) DX,DY,DDIR
6048   FORMAT ('GRID: REGULAR RECTANGULAR')                            
6049   FORMAT ('GRID: CURVILINEAR')                                    
6045   FORMAT (' S-low: ', F6.3,' S-hig: ', F6.3, ' frintf: ', F6.3)
6046   FORMAT (' MCGRD: ',I6,' MDC: ',I6,' MSC: ',I6)
6047   FORMAT (' DDIR: ', F6.3)
     ENDIF
!
     IF(.NOT.ALLOCATED(XCGRID)) ALLOCATE(XCGRID(MT))              
     IF(.NOT.ALLOCATED(YCGRID)) ALLOCATE(YCGRID(MT))              
!
     IF(OPTG == 1)THEN                                             
!      *** The coordinates of computational points in ***
!      *** regular grid are computed                  ***
       COSPC = COS(ALPC)                                               
       SINPC = SIN(ALPC)                                               
       DO I = 1,MT
         XCGRID(I)=VX(I)                                            
         YCGRID(I)=VY(I)                                            
       END DO   
     ENDIF
!
!    *** the computational grid is included in output data  ***        
     IF(MCGRD > 0)THEN                            
       ALLOCATE(OPSTMP)                                             
       OPSTMP%PSNAME = 'COMPGRID'                                   
       IF(OPTG == 1)THEN                                          
         OPSTMP%PSTYPE = 'F'                                        
         OPSTMP%OPR(1) = VXMIN       !XPC                                        
         OPSTMP%OPR(2) = VYMIN       !YPC                                        
         OPSTMP%OPR(3) = 1.0         !XCLEN                                      
         OPSTMP%OPR(4) = 1.0         !YCLEN                                      
         OPSTMP%OPR(5) = ALPC                                       
       ELSE                                                         
         OPSTMP%PSTYPE = 'H'                                        
         OPSTMP%OPR(1) = FLOAT(MXC-1)                                
         OPSTMP%OPR(2) = FLOAT(MYC-1)                                
         OPSTMP%OPR(3) = 0.                                          
         OPSTMP%OPR(4) = 0.                                          
         OPSTMP%OPR(5) = ALPC                                        
       ENDIF                                                          
       OPSTMP%OPI(1) = MCGRD       !MXC                                            
       ALLOCATE(OPSTMP%XP(0))                                         
       ALLOCATE(OPSTMP%YP(0))                                         
       NULLIFY(OPSTMP%NEXTOPS)                                        
       IF(.NOT.LOPS)THEN                                          
         FOPS = OPSTMP                                               
         COPS => FOPS                                                
         LOPS = .TRUE.                                               
       ELSE                                                           
         COPS%NEXTOPS => OPSTMP                                      
         COPS => OPSTMP                                              
       END IF                                                         
     ENDIF
     
!     ------------------------------------------------------------------
!
!     INPUT    definition of input grids
!
! ============================================================================
!
!   INPgrid                                                                &
!      BOTtom / WLEVel / CURrent / VX / VY / FRiction / WInd / WX / WY     &
!      / REG [xpinp] [ypinp] [alpinp]  [mxinp] [myinp]  [dxinp] [dyinp] \
!      \ CURVilinear [stagrx] [stagry] [mxinp] [myinp]                  /  &
!      (NONSTATionary [tbeginp] [deltinp] SEC/MIN/HR/DAY [tendinp])
!
! ============================================================================
   IGR2 = 0
   INP_BOT = .TRUE.
   IF(INP_BOT)THEN
     IGRD = 1
     IGR1 = 1
     LOGCOM(3) =.TRUE.                                                 
     PSNAME = 'BOTTGRID'
    
     ALLOCATE(DEPTH(MT))      ;  DEPTH    = 0.0

     DEPTH = H
   END IF

# if defined (WAVE_CURRENT_INTERACTION) && !defined (OFFLINE)
# if defined (WAVE_ONLY)
    ICUR = 0
# else
    ICUR = 1
# endif 
# endif 
   CALL INLOGC('INP_CUR',INP_CUR,'REQ','T')
   IF(INP_CUR)THEN
     IGRD = 1
     IGR1  = 2
     IGR2 = 3
     ICUR   = 1
     PSNAME = 'VXGRID  '

     CALL INCSTR('INP_CUR_NAME',INP_CUR_NAME,'REQ',' ')
     INQUIRE(FILE=TRIM(INP_CUR_NAME),EXIST=CHECK)
     IF(CHECK)THEN
       OPEN(1,FILE=TRIM(INP_CUR_NAME),STATUS='OLD',FORM='UNFORMATTED')
     ELSE
       WRITE(PRINTF,*) TRIM(INP_CUR_NAME), ' DOES NOT EXIT.'
       CALL PSTOP
     END IF
    
     CALL INLOGC('INP_CUR_STAT',INP_CUR_STAT,'REQ','T')
     IF(.NOT. INP_CUR_STAT)THEN
       IF(NSTATM == 0) CALL MSGERR(3,                                    &
              'In stationary mode, the current field also have to be stationary')
       NSTATM = 1
       CALL ININTG('CUR_TBEGINP',IFLBEG(IGRD),'REQ',0)
       CALL ININTG('CUR_DELTINP',IFLINT(IGRD),'REQ',0)
       CALL ININTG('CUR_TENDINP',IFLEND(IGRD),'REQ',0)
       IFLDYN(IGRD) = 1
       IFLTIM(IGRD) = IFLBEG(IGRD)
       IF(IGR2 > 0)THEN
         IFLBEG(IGR2) = IFLBEG(IGRD)                                 
         IFLINT(IGR2) = IFLINT(IGRD)                                 
         IFLEND(IGR2) = IFLEND(IGRD)                                 
         IFLDYN(IGR2) = IFLDYN(IGRD)                                 
         IFLTIM(IGR2) = IFLTIM(IGRD)                                 
       ENDIF
       
       CALL INLOGC('INP_CUR_SERI',INP_CUR_SERI,'REQ','T')
       IF(INP_CUR_SERI)THEN
         INP_CUR_NTIME = (IFLEND(IGRD)-IFLBEG(IGRD)+1)/IFLINT(IGRD)
         ALLOCATE(UXB_GL(MGL,INP_CUR_NTIME)); UXB_GL = 0.0
         ALLOCATE(UYB_GL(MGL,INP_CUR_NTIME)); UYB_GL = 0.0

         ALLOCATE(UXB(MT,INP_CUR_NTIME));  UXB = 0.0
         ALLOCATE(UYB(MT,INP_CUR_NTIME));  UYB = 0.0

	 DO J=1,INP_CUR_NTIME
	   READ(1) (UXB_GL(I,J),UYB_GL(I,J),I=1,MGL)
	 END DO
       ELSE	 
         INP_CUR_NTIME = 1
         ALLOCATE(UXB_GL(MGL,1)); UXB_GL = 0.0
         ALLOCATE(UYB_GL(MGL,1)); UYB_GL = 0.0
	 
         ALLOCATE(UXB(MT,1)); UXB = 0.0
         ALLOCATE(UYB(MT,1)); UYB = 0.0

	 READ(1) (UXB_GL(I,1),UYB_GL(I,1),I=1,MGL)
       END IF	 
     ELSE
       INP_CUR_NTIME = 1
       ALLOCATE(UXB_GL(MGL,1)); UXB_GL = 0.0
       ALLOCATE(UYB_GL(MGL,1)); UYB_GL = 0.0

       ALLOCATE(UXB(MT,1)); UXB = 0.0
       ALLOCATE(UYB(MT,1)); UYB = 0.0

       READ(1) (UXB_GL(I,1),UYB_GL(I,1),I=1,MGL)
     END IF  
     CLOSE(1)   
     
     IF(SERIAL)THEN
       UXB = UXB_GL
       UYB = UYB_GL
     END IF
     
#  if defined (MULTIPROCESSOR)
     IF(PAR)THEN
       DO J=1,INP_CUR_NTIME 
         DO I=1,M
	   UXB(I,J) = UXB_GL(NGID(I),J)
	   UYB(I,J) = UYB_GL(NGID(I),J)
	 END DO
	 DO I=1,NHN
	   UXB(I+M,J) = UXB_GL(HN_LST(I),J)
	   UYB(I+M,J) = UYB_GL(HN_LST(I),J)         
       	 END DO
       END DO
     END IF
#  endif       	 	   
   END IF
   IF(.NOT.ALLOCATED(UXB))    ALLOCATE(UXB(MT,1))
   IF(.NOT.ALLOCATED(UYB))    ALLOCATE(UYB(MT,1))
   
   VARWI = .FALSE.
   INP_WI = .FALSE.
   INP_WI_SERI = .FALSE.

   IF(WIND_ON .AND. WIND_KIND == VRBL) INP_WI = .TRUE.
   IF(INP_WI)THEN
     LWINDR = 2                                                        
     IWIND  = LWINDM                                                   
     IGRD = 5
     IGR1   = 5
     IGR2 = 6
     VARWI  = .TRUE.
     PSNAME = 'WVGRID  '
     INP_WI_SERI = .TRUE.
     
     NSTATM = 1
     CALL SURFACE_WIND2WAVE(IGRD)

     IFLDYN(IGRD) = 1
     IFLTIM(IGRD) = IFLBEG(IGRD)
     IF(IGR2 > 0)THEN
       IFLBEG(IGR2) = IFLBEG(IGRD)                                 
       IFLINT(IGR2) = IFLINT(IGRD)                                 
       IFLEND(IGR2) = IFLEND(IGRD)                                 
       IFLDYN(IGR2) = IFLDYN(IGRD)                                 
       IFLTIM(IGR2) = IFLTIM(IGRD)                                 
     ENDIF
   END IF
   
   CALL INLOGC('INP_FR',INP_FR,'REQ','T')
   IF(INP_FR)THEN
     IGRD = 4
     IGR1   = 4
     VARFR  = .TRUE.
     PSNAME = 'FRICGRID'
     
     CALL INCSTR('INP_FR_NAME',INP_FR_NAME,'REQ',' ')
     INQUIRE(FILE=TRIM(INP_FR_NAME),EXIST=CHECK)
     IF(CHECK)THEN
       OPEN(1,FILE=TRIM(INP_FR_NAME),STATUS='OLD',FORM='UNFORMATTED')
     ELSE
       WRITE(PRINTF,*) TRIM(INP_FR_NAME), ' DOES NOT EXIT.'
       CALL PSTOP
     END IF
    
     CALL INLOGC('INP_FR_STAT',INP_FR_STAT,'REQ','T')
     IF(.NOT. INP_FR_STAT)THEN
       IF(NSTATM == 0) CALL MSGERR(3,                                    &
              'In stationary mode, the friction field also have to be stationary')
       NSTATM = 1
       CALL ININTG('FR_TBEGINP',IFLBEG(IGRD),'REQ',0)
       CALL ININTG('FR_DELTINP',IFLINT(IGRD),'REQ',0)
       CALL ININTG('FR_TENDINP',IFLEND(IGRD),'REQ',0)
       IFLDYN(IGRD) = 1
       IFLTIM(IGRD) = IFLBEG(IGRD)
       IF(IGR2 > 0)THEN
         IFLBEG(IGR2) = IFLBEG(IGRD)                                 
         IFLINT(IGR2) = IFLINT(IGRD)                                 
         IFLEND(IGR2) = IFLEND(IGRD)                                 
         IFLDYN(IGR2) = IFLDYN(IGRD)                                 
         IFLTIM(IGR2) = IFLTIM(IGRD)                                 
       ENDIF
       
       CALL INLOGC('INP_FR_SERI',INP_FR_SERI,'REQ','T')
       IF(INP_FR_SERI)THEN
         INP_FR_NTIME = (IFLEND(IGRD)-IFLBEG(IGRD)+1)/IFLINT(IGRD)
         ALLOCATE(FRIC_GL(MGL,INP_FR_NTIME)); FRIC_GL = 0.0
	 ALLOCATE(FRIC(MT,INP_FR_NTIME));     FRIC    = 0.0
	 DO J=1,INP_FR_NTIME
	   READ(1) (FRIC_GL(I,J),I=1,MGL)
	 END DO
       ELSE	 
         INP_FR_NTIME = 1
         ALLOCATE(FRIC_GL(MGL,1));     FRIC_GL = 0.0
	 ALLOCATE(FRIC(MT,1));         FRIC    = 0.0
	 READ(1) (FRIC_GL(I,1),I=1,MGL)
       END IF	 
     ELSE
       INP_FR_NTIME = 1
       ALLOCATE(FRIC_GL(MGL,1));       FRIC_GL = 0.0
       ALLOCATE(FRIC(MT,1));           FRIC    = 0.0
       READ(1) (FRIC_GL(I,1),I=1,MGL)
     END IF  
     CLOSE(1)
     
     IF(SERIAL)THEN
       FRIC = FRIC_GL
     END IF
     
#    if defined (MULTIPROCESSOR)
     IF(PAR)THEN
       DO J=1,INP_FR_NTIME
         DO I=1,M
	   FRIC(I,J) = FRIC_GL(NGID(I),J)
	 END DO
	 DO I=1,NHN
	   FRIC(I+M,J) = FRIC_GL(HN_LST(I),J)
	 END DO
       END DO
     END IF
#    endif       	            
   END IF
   IF(.NOT.ALLOCATED(FRIC))    ALLOCATE(FRIC(MT,1))
   
   CALL INLOGC('INP_WLEV',INP_WLEV,'REQ','T')
   IF(INP_WLEV)THEN
     IGRD = 7
     IGR1   = 7
     VARWLV = .TRUE.                                                   
     PSNAME = 'WLEVGRID'
     
     CALL INCSTR('INP_WLEV_NAME',INP_WLEV_NAME,'REQ',' ')
     INQUIRE(FILE=TRIM(INP_WLEV_NAME),EXIST=CHECK)
     IF(CHECK)THEN
       OPEN(1,FILE=TRIM(INP_WLEV_NAME),STATUS='OLD',FORM='UNFORMATTED')
     ELSE
       WRITE(PRINTF,*) TRIM(INP_WLEV_NAME), ' DOES NOT EXIT.'
       CALL PSTOP
     END IF
    
     CALL INLOGC('INP_WLEV_STAT',INP_WLEV_STAT,'REQ','T')
     IF(.NOT. INP_WLEV_STAT)THEN
       IF(NSTATM == 0) CALL MSGERR(3,                                    &
              'In stationary mode, the water level also have to be stationary')
       NSTATM = 1
       CALL ININTG('WLEV_TBEGINP',IFLBEG(IGRD),'REQ',0)
       CALL ININTG('WLEV_DELTINP',IFLINT(IGRD),'REQ',0)
       CALL ININTG('WLEV_TENDINP',IFLEND(IGRD),'REQ',0)
       IFLDYN(IGRD) = 1
       IFLTIM(IGRD) = IFLBEG(IGRD)
       IF(IGR2 > 0)THEN
         IFLBEG(IGR2) = IFLBEG(IGRD)                                 
         IFLINT(IGR2) = IFLINT(IGRD)                                 
         IFLEND(IGR2) = IFLEND(IGRD)                                 
         IFLDYN(IGR2) = IFLDYN(IGRD)                                 
         IFLTIM(IGR2) = IFLTIM(IGRD)                                 
       ENDIF
       
       CALL INLOGC('INP_WLEV_SERI',INP_WLEV_SERI,'REQ','T')
       IF(INP_WLEV_SERI)THEN
         INP_WLEV_NTIME = (IFLEND(IGRD)-IFLBEG(IGRD)+1)/IFLINT(IGRD)
         ALLOCATE(WLEVL_GL(MGL,INP_WLEV_NTIME)); WLEVL_GL = 0.0
	 ALLOCATE(WLEVL(MT,INP_WLEV_NTIME));     WLEVL    = 0.0
	 DO J=1,INP_WLEV_NTIME
	   READ(1) (WLEVL_GL(I,J),I=1,MGL)
	 END DO
       ELSE	 
         INP_WLEV_NTIME = 1
         ALLOCATE(WLEVL_GL(MGL,1));    WLEVL_GL = 0.0
	 ALLOCATE(WLEVL(MT,1));        WLEVL    = 0.0
	 READ(1) (WLEVL_GL(I,1),I=1,MGL)
       END IF	 
     ELSE
       INP_WLEV_NTIME = 1
       ALLOCATE(WLEVL_GL(MGL,1));    WLEVL_GL = 0.0
       ALLOCATE(WLEVL(MT,1));        WLEVL    = 0.0
       READ(1) (WLEVL_GL(I,1),I=1,MGL)
     END IF  
     CLOSE(1)
     
     IF(SERIAL)THEN
       WLEVL = WLEVL_GL
     END IF
     
#    if defined (MULTIPROCESSOR)
     IF(PAR)THEN
       DO J=1,INP_WLEV_NTIME
         DO I=1,M
	   WLEVL(I,J) = WLEVL_GL(NGID(I),J)
	 END DO
	 DO I=1,NHN
	   WLEVL(I+M,J) = WLEVL_GL(HN_LST(I),J)
	 END DO
       END DO
     END IF
#    endif       	          
   END IF
   IF(.NOT.ALLOCATED(WLEVL)) ALLOCATE(WLEVL(MT,1))
 
!
!     ------------------------------------------------------------------
!
!     WIND      parameters uniform wind field
!
! ==========================================
!
!   WIND  [vel]  [dir]  [astd]                                            
!
! ==========================================
!
!   CALL INLOGC('WIND',VARWI1,'UNC','T')
    IF(WIND_ON .AND. WIND_KIND == CNSTNT) VARWI1 = .TRUE.
    IF(VARWI1)THEN
     VARWI = .NOT. VARWI1
     LWINDR = 1
     IWIND  = LWINDM
     
     U10  = SQRT(WIND_X*WIND_X+WIND_Y*WIND_Y)
     WDIP = ATAN2(WIND_Y,WIND_X)
     WDIP = WDIP*180.0/PI_W
!
!    *** Convert (if necessary) WDIP from nautical degrees ***         
!    *** to cartesian degrees                              ***         
!
!
     WDIP = DEGCNV (WDIP)
!
     IF(IWIND == 0) IWIND = 4                                         
     ALTMP = WDIP / 360.
     WDIP = PI2_W * (ALTMP - NINT(ALTMP))
     IF(JWX2 <= 1)THEN                                               
       MCMVAR = MCMVAR + 2                                             
       JWX2   = MCMVAR - 1
       JWY2   = MCMVAR
     ENDIF
   ENDIF

!
!  -----------------------------------------------------------------------
!
!  BOUNdary  defining boundary conditions                              
!
!==============================================
!
   NESTING = .FALSE.
   CALL INCSTR('BOUNDARY',BOUNDARY,'REQ',' ')
   
   IF(BOUNDARY == 'default' .and. NESTING_ON_WAVE)   &
       CALL FATAL_ERROR("The parameter NESTING_ON_WAVE in ***_run.nml should be .FALSE.",  &
                        "or parameter BOUNDARY in INPUT should be defined.")
   
   IF(BOUNDARY /= 'default')THEN

!
!----------------Determine Number of Nodes on Outer Boundary-------------------!
!
     IOBCN_GL_W = IOBCN_GL

     IOBCN_W = 0

     IF(IOBCN_GL_W > 0)THEN

!------------Read in Open Boundary Nodes and Temperature/Salinity Conditions---!

       ALLOCATE(I_OBC_GL_W(IOBCN_GL_W))

       I_OBC_GL_W = I_OBC_GL

!----------------------Make Sure It Is In Global Domain------------------------!

       DO I=1,IOBCN_GL_W
         IF((I_OBC_GL_W(I) > MGL))THEN
           WRITE(IPT,*)'==================ERROR=================================='
           WRITE(IPT,*)'OPEN BOUNDARY NODE NUMBER',I,'IS NOT IN THE'
           WRITE(IPT,*)'GLOBAL DOMAIN'
           WRITE(IPT,*)'CHECK INPUT FILE AND ENSURE OPEN BOUNDARY NODES <= ',MGL
           WRITE(IPT,*)'========================================================='
           CALL PSTOP
         END IF
       END DO

!----------Shift Open Boundary Node List,Type-----------!

       IF(SERIAL)THEN

         IOBCN_W    = IOBCN_GL_W

         ALLOCATE(I_OBC_N_W(IOBCN_W))
         I_OBC_N_W = I_OBC_GL_W
       END IF

#      if defined (MULTIPROCESSOR)
       IF(PAR)THEN
         ALLOCATE(TEMP1(IOBCN_GL_W))

         NCNT = 0
         !!SET UP LOCAL OPEN BOUNDARY NODES
         DO I=1,IOBCN_GL_W
           I1 = NLID( I_OBC_GL_W(I) )
           IF(I1 /= 0)THEN
             NCNT = NCNT + 1
             TEMP1(NCNT) = I1
           END IF
         END DO
         IOBCN_W = NCNT

         IF(NCNT > 0)THEN
           ALLOCATE(I_OBC_N_W(NCNT))
           I_OBC_N_W  = TEMP1(1:NCNT)
         END IF

         DEALLOCATE(TEMP1)
       END IF
#      endif

!adjust ISONB_W according to FVCOM
     ISONB_W = ISONB

     DO I=1,MT
        IF(ISONB(I)/=0) ISONB_W(I)=1 !separate intetior points and outer boundary
     ENDDO

     DO I=1,IOBCN_W
        I1=I_OBC_N_W(I)
        ISONB_W(I1)=2                ! seprate open boundary from outer boundary
     ENDDO
         
     END IF !!IOBCN_GL_W > 0

     CALL SWBOUN ( )  
   ELSE
     IOBCN_W = 0
     ISONB_W = ISONB
     DO I=1,MT
        IF(ISONB_W(I)==2) ISONB_W(I)=1
        IF(ISONB_W(I)==3) ISONB_W(I)=0
     ENDDO
           
   END IF
   
!
!  ------------------------------------------------------------------
!
! INITial conditions  Definition of initial conditions for MODE DYNAMIC
!
   IF(KEYWIS('INIT'))THEN                                           
!JQI     CALL INITVA( AC2   , SPCSIG, SPCDIR, KGRPNT, XCGRID, YCGRID, XYTST )
!JQI     IF(STPNOW()) RETURN                                              
!JQI     GO TO 100
   ENDIF
!     ------------------------------------------------------------------
!
!     GEN       the generation of mode                          
!
!     ------------------------------------------------------------------

   CALL ININTG('GEN',GEN,'REQ',0)

!
!     GEN=1        deep water with 1st generation                          
!
! ==========================================
!
! GEN=1     [cf10]  [cf20]  [cf30]  [cf40]  [edmlpm]  [cdrag]  [umin]  [cfpm]
!
! ==========================================
!
   IF(GEN == 1)THEN
     LWINDM = 1                                                        
     IF(LWINDR > 0) IWIND = LWINDM                                 
     IQUAD = 0
     IWCAP = 0

     CALL INREAL('GEN1_CF10',  PWIND(1), 'UNC',0.0)
     CALL INREAL('GEN1_CF20',  PWIND(2), 'UNC',0.0)
     CALL INREAL('GEN1_CF30',  PWIND(3), 'UNC',0.0)
     CALL INREAL('GEN1_CF40',  PWIND(4), 'UNC',0.0)
     CALL INREAL('GEN1_EDMLPM',PWIND(10),'UNC',0.0)
     CALL INREAL('GEN1_CDRAG', PWIND(11),'UNC',0.0)
     CALL INREAL('GEN1_UMIN',  PWIND(12),'UNC',0.0)
     CALL INREAL('GEN1_CFPM',  PWIND(13),'UNC',0.0)
      
!    if [pwtail] is not changed, make it 5
     IF(EQREAL(PWTAIL(1),4.))THEN                                    
       PWTAIL(1) = 5.                                                  
       PWTAIL(3) = PWTAIL(1) + 1.                                      
     ENDIF
   END IF     

!
!     ------------------------------------------------------------------
!
!     GEN2        deep water with 2nd generation                          
!
! ==========================================
!
! GEN2 [cf10] [cf20] [cf30] [cf40] [cf50] [cf60] [edmlpm] [cdrag] [umin] [cfpm]
!
! ==========================================
!
   IF(GEN == 2)THEN
     LWINDM = 2                                                        
     IF(LWINDR > 0) IWIND = LWINDM                                 
     IQUAD = 0
     IWCAP = 0

     CALL INREAL('GEN2_CF10',  PWIND(1), 'UNC',0.0)
     CALL INREAL('GEN2_CF20',  PWIND(2), 'UNC',0.0)
     CALL INREAL('GEN2_CF30',  PWIND(3), 'UNC',0.0)
     CALL INREAL('GEN2_CF40',  PWIND(4), 'UNC',0.0)
     CALL INREAL('GEN2_CF50',  PWIND(5), 'UNC',0.0)
     CALL INREAL('GEN2_CF60',  PWIND(6), 'UNC',0.0)
     CALL INREAL('GEN2_EDMLPM',PWIND(10),'UNC',0.0)
     CALL INREAL('GEN2_CDRAG', PWIND(11),'UNC',0.0)
     CALL INREAL('GEN2_UMIN',  PWIND(12),'UNC',0.0)
     CALL INREAL('GEN2_CFPM',  PWIND(13),'UNC',0.0)
      
!    if [pwtail] is not changed, make it 5
     IF(EQREAL(PWTAIL(1),4.))THEN                                    
       PWTAIL(1) = 5.                                                  
       PWTAIL(3) = PWTAIL(1) + 1.                                      
     ENDIF
   END IF     

!
!     ------------------------------------------------------------------
!
!     GEN3        deep water with 3d generation                           
!
! ======================================================================================
!
!       |   JANSsen [cds1] [delta]                         |              
! GEN3 <                                                    >     &
!       | ->KOMen   [cds2] [stpm] [powst] [delta] [powk]   |              
!       |                                                  |
!       |   YAN     (NOT documented)                       |
!       |                                                  |
!       |   WESTHuysen [cds2] [br] [p0] [powst] [powk]     |              
!
!  (QUADrupl [iquad] [limiter] [lambda] [cnl4] [csh1] [csh2] [csh3]) (AGROW [a])
!
! ======================================================================================
!
   IF(GEN == 3)THEN
     GROWTH = 'KOM'
     CALL INCSTR('GROWTH',GROWTH,'STA','KOM')
     
     IF(GROWTH == 'JANS')THEN
       LWINDM = 4
       IF(LWINDR > 0) IWIND = LWINDM
!      *** whitecapping according to Janssen (1989, 1991) ***
       IWCAP = 2
       CALL INREAL('GEN3_JANS_CDS1', PWCAP(3),'UNC',0.0)
       CALL INREAL('GEN3_JANS_DELTA',PWCAP(4),'UNC',0.0)
       
!
!      Recalculate coefficientes that are actually used in the         
!      whitecapping routines                                           
       PWCAP(1)  = PWCAP(3) * (PWCAP(2) ** PWCAP(9))                   
       PWCAP(10) = PWCAP(4)                                            
       JUSTAR = MCMVAR+1                                               
       JZEL   = MCMVAR+2                                               
       JCDRAG = MCMVAR+3                                               
       JTAUW  = MCMVAR+4                                               
       MCMVAR = MCMVAR+4                                               
!      if [pwtail] is not changed, make it 5
       IF(EQREAL(PWTAIL(1),4.))THEN                                  
         PWTAIL(1) = 5.                                                
         PWTAIL(3) = PWTAIL(1) + 1.                                    
       ENDIF
   
     ELSE IF(GROWTH == 'YAN')THEN
!      option not documented in user manual
       LWINDM = 5                                                      
       IF(LWINDR > 0) IWIND = LWINDM                               
     ELSE IF(GROWTH == 'WESTH')THEN                                    
!      *** wind according to (adapted) Yan (1987) ***                  
       LWINDM = 5                                                      
       IF(LWINDR > 0) IWIND = LWINDM                               
!      *** whitecapping according to Alves and Banner (2003) ***       
       IWCAP = 7                                                       

       CALL INREAL('GEN_WESTH_CDS2', PWCAP(1), 'STA',5.0E-5)
       CALL INREAL('GEN_WESTH_BR',   PWCAP(12),'STA',1.75E-3)                
       CALL INREAL('GEN_WESTH_P0',   PWCAP(10),'STA',4.)                     
       CALL INREAL('GEN_WESTH_POWST',PWCAP(9), 'STA',0.)                     
       CALL INREAL('GEN_WESTH_POWK', PWCAP(11),'STA',0.)                     
     ELSE IF(GROWTH == 'KOM')THEN
       LWINDM = 3                                                      
       IF(LWINDR > 0) IWIND = LWINDM                               
!      *** whitecapping according to Komen et al. (1984) ***
       IWCAP = 1

       CALL INREAL('GEN3_KOM_CDS2',PWCAP(1),'UNC',0.0)
       CALL INREAL('GEN3_KOM_STPM',PWCAP(2),'UNC',0.0)
     ELSE
       WRITE(*,*)
       CALL PSTOP      
     ENDIF

!    ***  parameters nonlinear 4 wave interactions ***
     QUAD = .TRUE.
     CALL INLOGC('QUAD',QUAD,'STA','T')
     
     IF(QUAD)THEN
       CALL ININTG('IQUAD',IQUAD,'UNC',0)
      
!      *** if quadruplets are activated then LIMITER = 0.1 ***
!      *** the standard value, however, is 1.e20           ***

       CALL INREAL('LAMBDA',PQUAD(1),'UNC',0.0)
       CALL INREAL('CNL4',  PQUAD(2),'UNC',0.0)
       CALL INREAL('CSH1',  PQUAD(3),'UNC',0.0)
       CALL INREAL('CSH2',  PQUAD(4),'UNC',0.0)
       CALL INREAL('CSH3',  PQUAD(5),'UNC',0.0)
     ELSE
       IQUAD = 0
     END IF 
     
     AGROW = .FALSE.
     CALL INLOGC('AGROW',AGROW,'STA','F')

     IF(AGROW)THEN
       CALL INREAL('A',PWIND(31),'STA',0.0015)
     ELSE
       IF(NSTATM == 1 .AND. ICOND == 0) ICOND = 1
     END IF  
   END IF
!
!     ------------------------------------------------------------------
!
!     WCAP        parameters whitecapping
!
! =============================================================
!
!        | ->KOMen    [cds2] [stpm] [powst] [delta] [powk]  |             
!        |                                                  |
!        |   JANSsen  [cds1]  [delta] [pwtail]              |
! WCAP  <                                                    >    (NOT documented)
!        |   LHIG     [cflhig]                              |
!        |                                                  |
!        |   BJ       [bjstp] [bjalf]                       |
!        |                                                  |
!        |   KBJ      [bjstp] [bjalf] [kconv]               |
!        |                                                  |             
!        |   CSM      [cst] [pow]                           |             
!        |                                                  |
!        |   AB       [cds2] [br] [p0] [powst] [powk]       |             
!
! =============================================================
!
   WCAP = 'KOM '
   CALL INCSTR('WCAP',WCAP,'STA','KOM')

   IF(WCAP == 'KOM')THEN
!    *** whitecapping according to Komen et al. (1984) ***
     IWCAP = 1
     CALL INREAL('WCAP_KOM_CDS2', PWCAP(1), 'UNC',0.0)
     CALL INREAL('WCAP_KOM_STPM', PWCAP(2), 'UNC',0.0)
     CALL INREAL('WCAP_KOM_POWST',PWCAP(9), 'UNC',0.0)
     CALL INREAL('WCAP_KOM_DELTA',PWCAP(10),'UNC',0.0)
     CALL INREAL('WCAP_KOM_POWK', PWCAP(11),'UNC',0.0)
   ELSE IF(WCAP == 'JANS')THEN
!    *** whitecapping according to Janssen (1989, 1991) ***
     IWCAP = 2
     CALL INREAL('WCAP_JANS_CDS1', PWCAP(3),'UNC',0.0)  
     CALL INREAL('WCAP_JANS_DELTA',PWCAP(4),'UNC',0.0)  
!
!    Recalculate coefficients that are actually used in the          
!    whitecapping routines                                           
     PWCAP(1)  = PWCAP(3) * (PWCAP(2) ** PWCAP(9))                   
     PWCAP(10) = PWCAP(4)                                            
     JUSTAR = MCMVAR+1                                               
     JZEL   = MCMVAR+2                                               
     JCDRAG = MCMVAR+3                                               
     JTAUW  = MCMVAR+4                                               
     MCMVAR = MCMVAR+4                                               

     CALL INREAL('WCAP_JANS_PWTAIL',PWTAIL(1),'STA',5.0)
     IF(PWTAIL(1) <= 1.) CALL MSGERR (3, 'Incorrect PWTAIL')        
     PWTAIL(3) = PWTAIL(1) + 1.                                      
   ELSE IF(WCAP == 'LHIG')THEN
!    *** whitecapping according to Longuett Higgins ***
     IWCAP = 3
     PWCAP(5) = 1.0
     CALL INREAL('WCAP_LHIG_CFLHIG',PWCAP(5),'UNC',0.0)
   ELSE IF(WCAP == 'BJ')THEN
!    *** whitecapping according to Battjes/Janssen formulation ***
     IWCAP = 4
     PWCAP(6) = 0.88
     CALL INREAL('WCAP_BJ_BJSTP',PWCAP(6),'UNC',0.0)
     PWCAP(7) = 1.0
     CALL INREAL('WCAP_BJ_BJALF',PWCAP(7),'UNC',0.0)
   ELSE IF(WCAP == 'KBJ')THEN
!    *** whitecapping according to a combination of Komen et al ***
!    *** and Battjes/Janssen                                    ***
     IWCAP = 5
     PWCAP(6) = 0.88
     CALL INREAL('WCAP_KBJ_BJSTP',PWCAP(6),'UNC',0.0)
     PWCAP(7) = 1.0
     CALL INREAL('WCAP_KBJ_BJALF',PWCAP(7),'UNC',0.0)
     PWCAP(8) = 0.75
     CALL INREAL('WCAP_KBJ_KCONV',PWCAP(8),'UNC',0.0)
   ELSE IF(WCAP == 'CSM')THEN                                    
!    *** whitecapping according to cumulative steepness method  ***  
     IWCAP = 6                                                       
     PWCAP(12) = 4.0
     CALL INREAL('WCAP_CSM_CST',PWCAP(12),'UNC',0.0)
     PWCAP(13) = 2.0
     CALL INREAL('WCAP_CSM_POW',PWCAP(13),'UNC',0.0)
   ELSE IF(WCAP == 'AB')THEN                                     
!    *** whitecapping according to Alves and Banner (2003) ***       
     IWCAP = 7                                                       
     CALL INREAL('WCAP_AB_CDS2',PWCAP(1), 'STA',5.0E-5)
     CALL INREAL('WCAP_AB_BR',  PWCAP(12),'STA',1.75E-3)
     CALL INREAL('WCAP_AB_PO',  PWCAP(10),'STA',4.0)
     CALL INREAL('WCAP_AB_POWST',PWCAP(9),'STA',0.0)
     CALL INREAL('WCAP_AB_POWK',PWCAP(11),'STA',0.0)
   ELSE IF(WCAP == 'OFF')THEN
     IWCAP = 0
   ELSE
!JQI          CALL WRNKEY
   ENDIF

!     ------------------------------------------------------------------
!     ***  parameters nonlinear 4 wave interactions ***
     QUAD = .TRUE.
     CALL INLOGC('QUAD',QUAD,'STA','T')
     
     IF(QUAD)THEN
!
! ===================================================================
!
!  QUADrupl [iquad] [limiter] [lambda] [cnl4] [csh1] [csh2] [csh3]        
!
! ===================================================================
!
       CALL ININTG('IQUAD',IQUAD,'UNC',0)
      
!      *** if quadruplets are activated then LIMITER = 0.1 ***
!      *** the standard value, however, is 1.e20           ***
       CALL INREAL('LIMITER',PNUMS(20),'STA',0.1)           
       CALL INREAL('LAMBDA', PQUAD(1), 'UNC',0.0)
       CALL INREAL('CNL4',   PQUAD(2), 'UNC',0.0)
       CALL INREAL('CSH1',   PQUAD(3), 'UNC',0.0)
       CALL INREAL('CSH2',   PQUAD(4), 'UNC',0.0)
       CALL INREAL('CSH3',   PQUAD(5), 'UNC',0.0)
     ELSE
       IQUAD = 0
     ENDIF

! ----------------------------------------------------------------    
!                            | CNL4 < [cnl4] >               |         
!  MDIA LAMbda < [lambda] > <                                 >        
!                            | CNL4_12 < [cnl4_1] [cnl4_2] > |         
! ----------------------------------------------------------------    

   MDIAL = .FALSE.
   CALL INLOGC('MDIA',MDIAL,'UNC','F')

   IF(MDIAL)THEN                                            
     IF(ALLOCATED(LAMBDA))THEN                                       
       DEALLOCATE (LAMBDA,CNL4_1,CNL4_2)                               
     ENDIF                                                             
     ALLOCATE (RLAMBDA(1000))                                          

     LAM = .FALSE.
     CALL INLOGC('MDIA_LAM',LAM,'STA','F')
     
     IF(LAM)THEN                                           
       MORE = .TRUE.                                                   
       ILAMBDA = 0                                                     
       DO WHILE (MORE)                                                 
         CALL INREAL('MDIA_LAMBDA',RLAMBDA(ILAMBDA+1),'STA',-1.)            
         IF(RLAMBDA(ILAMBDA+1) < 0.) THEN                            
           MORE = .FALSE.                                              
         ELSE                                                          
           ILAMBDA = ILAMBDA + 1                                       
         ENDIF                                                         
       ENDDO                                                           
       MDIA = ILAMBDA                                                  
       ALLOCATE (LAMBDA(MDIA))                                         
       LAMBDA(1:MDIA) = RLAMBDA(1:MDIA)                                
       DEALLOCATE (RLAMBDA)                                            
     ELSE                                                              
!JQI   CALL WRNKEY                                                     
     END IF                                                            
     ALLOCATE (CNL4_1(MDIA),CNL4_2(MDIA))                              
     CALL INCSTR('MDIA_CNL4C',CNL4,'UNC',' ')
     IF(CNL4 == 'CNL4_12')THEN                                       
       DO ICNL4=1,MDIA
         CALL INREAL('MDIA_CNL4_1',CNL4_1(ICNL4),'REQ',0.)                  
         CALL INREAL('MDIA_CNL4_2',CNL4_2(ICNL4),'REQ',0.)                  
       END DO                                                          
     ELSE IF(CNL4 == 'CNL4')THEN                                      
       DO ICNL4=1,MDIA                                                 
         CALL INREAL('MDIA_CNL4',CNL4_1(ICNL4),'REQ',0.)                    
       END DO                                                          
       CNL4_2 = CNL4_1                                                 
     ELSE                                                              
!JQI          CALL WRNKEY                                                     
     END IF                                                            
     CNL4_1 = CNL4_1 * ((2.*PI_W)**9)                                    
     CNL4_2 = CNL4_2 * ((2.*PI_W)**9)                                    
   ENDIF                                                               

!
! ------------------------------------------------------------------
!
!     BREAK     parameters surf breaking
!
! ============================================================
!
!           | -> CONstant [alpha] [gamma]                                      |
! BREaking <                                                                    >
!           |    VARiable [alpha] [gammin] [gammax] [gamneg] [coeff1] [coeff2] |
!
! ============================================================
!
  BRE = .FALSE.
  CALL INLOGC('BRE',BRE,'STA','F')
  IF(BRE)THEN
    CALL INCSTR('CONSTANT',CONSTANT,'STA','CON')

    IF(CONSTANT == 'CON')THEN
      ISURF = 1
      CALL INREAL('BRE_CON_ALPHA',PSURF(1),'STA',1.0)
      CALL INREAL('BRE_CON_GAMMA',PSURF(2),'STA',0.73)
    ELSE IF(CONSTANT == 'VAR')THEN                                      
      ISURF = 2
      CALL INREAL('BRE_VAR_ALPHA', PSURF(1),'STA',1.5)
      CALL INREAL('BRE_VAR_GAMMIN',PSURF(4),'STA',0.55)
      CALL INREAL('BRE_VAR_GAMMAX',PSURF(5),'STA',0.81)
      CALL INREAL('BRE_VAR_GAMNEG',PSURF(6),'STA',0.73)
      CALL INREAL('BRE_VAR_COEFF1',PSURF(7),'STA',0.88)
      CALL INREAL('BRE_VAR_COEFF2',PSURF(8),'STA',0.012)
    ENDIF
  ELSE
    ISURF = 0
  ENDIF
!
! ------------------------------------------------------------------
!
!     FRICTION  setting parameters for bottom friction
!
! ===============================================
!
!                   | -> JONSWAP       [cfjon]
!                   |
!   FRICTION       <     COLLINS       [cfw]    &
!                   |
!                   |                  [cfc]    (NOT documented)
!                   |
!                   |    MADSEN        [kn]
!
! ===============================================
!
   FRICL = .FALSE.
   CALL INLOGC('FRICTION',FRICL,'UNC','F')
   
   IF(FRICL)THEN
     IBOT = 1
     CALL INCSTR('FRIC_FORM',FRIC_FORM,'STA','JON')
     IF(FRIC_FORM == 'JON')THEN
       IBOT = 1
       CALL INREAL('CFJON',PBOT(3),'UNC',0.)
     ELSE IF(FRIC_FORM == 'COLL')THEN                                     
       IBOT = 2
       CALL INREAL('CFW',PBOT(2),'UNC',0.)
       CALL INREAL('CFC',PBOT(1),'UNC',0.)
     ELSE IF(FRIC_FORM == 'MAD')THEN
       IBOT = 3
       CALL INREAL('KN',PBOT(5),'UNC',0.)
     ELSE
!JQI       CALL WRNKEY
     ENDIF
   ENDIF

#  if defined (VEGETATION)
!
!   ------------------------------------------------------------------
!   VEGEtation   < [height]  [diamtr]  [nstems]  [drag] >
!   ------------------------------------------------------------------
   VEGE = .FALSE.
   CALL INLOGC('VEGETATION',VEGE,'UNC','F')
   IF(VEGE)THEN
        IVEG  = 1
        ILMAX = NVEG
!!$        ILMAX = 0
!!$        FRSTV%H = 0.
!!$        FRSTV%D = 0.
!!$        FRSTV%N = 0
!!$        FRSTV%C = 0.
!!$        NULLIFY(FRSTV%NEXTV)
!!$        CURRV => FRSTV
!!$ 303    CALL INREAL ('HEIGHT',HH,'REP',-1.)
!!$        IF ( HH.NE.-1. ) THEN
!!$           IF ( HH.LT.0. ) THEN
!!$              CALL MSGERR (2,'height is negative')
!!$              HH = 0.
!!$           END IF
!!$           CALL INREAL ('DIAMTR',DD,'REQ', 0.)
!!$           IF ( DD.LT.0. ) THEN
!!$              CALL MSGERR (2,'stem diameter is negative')
!!$              DD = 0.
!!$           END IF
!!$           CALL ININTG ('NSTEMS',NN,'REQ', 1 )
!!$           IF ( NN.LE.0 ) THEN
!!$              CALL MSGERR (2,'number of stems is negative or zero')
!!$              NN = 1
!!$           END IF
!!$           CALL INREAL ('DRAG'  ,CC,'REQ', 0.)
!!$           IF ( CC.LT.0. ) THEN
!!$              CALL MSGERR (2,'drag coefficient is negative')
!!$              CC = 0.
!!$           END IF
!!$           ILMAX = ILMAX+1
!!$           ALLOCATE(TMPV)
!!$           TMPV%H = HH
!!$           TMPV%D = DD
!!$           TMPV%N = NN
!!$           TMPV%C = CC
!!$           NULLIFY(TMPV%NEXTV)
!!$           CURRV%NEXTV => TMPV
!!$           CURRV => TMPV
!!$           GO TO 303
!!$        END IF
!!$        IF (ILMAX.EQ.0) CALL MSGERR(2,'No vegetation parameters found')
!!$        IF (.NOT.ALLOCATED(LAYH  )) ALLOCATE(LAYH  (ILMAX))
!!$        IF (.NOT.ALLOCATED(VEGDIL)) ALLOCATE(VEGDIL(ILMAX))
!!$        IF (.NOT.ALLOCATED(VEGNSL)) ALLOCATE(VEGNSL(ILMAX))
!!$        IF (.NOT.ALLOCATED(VEGDRL)) ALLOCATE(VEGDRL(ILMAX))
!!$        CURRV => FRSTV%NEXTV
!!$        DO JJ = 1, ILMAX
!!$           LAYH  (JJ) = CURRV%H
!!$           VEGDIL(JJ) = CURRV%D
!!$           VEGNSL(JJ) = REAL(CURRV%N)
!!$           VEGDRL(JJ) = CURRV%C
!!$           CURRV => CURRV%NEXTV
!!$        END DO
!!$        DEALLOCATE(TMPV)

        IF (.NOT.ALLOCATED(LAYH  )) ALLOCATE(LAYH  (0:MT,ILMAX))
        IF (.NOT.ALLOCATED(VEGDIL)) ALLOCATE(VEGDIL(0:MT,ILMAX))
        IF (.NOT.ALLOCATED(VEGNSL)) ALLOCATE(VEGNSL(0:MT,ILMAX))
        IF (.NOT.ALLOCATED(VEGDRL)) ALLOCATE(VEGDRL(0:MT,ILMAX))

        DO I = 1, ILMAX
            LAYH  (:,I) = veg%plant(:,I,phght)
            VEGDIL(:,I) = veg%plant(:,I,pdiam)
            VEGNSL(:,I) = veg%plant(:,I,pdens)
            VEGDRL(:,I) = CD_VEG(I)
        END DO
!!$# if defined (MULTIPROCESSOR)
!!$        IF(PAR)THEN
!!$          CALL AEXCHANGE(NC,MYID,NPROCS,LAYH)  
!!$          CALL AEXCHANGE(NC,MYID,NPROCS,VEGDIL)  
!!$          CALL AEXCHANGE(NC,MYID,NPROCS,VEGNSL)  
!!$          CALL AEXCHANGE(NC,MYID,NPROCS,VEGDRL)  
!!$        END IF
!!$# endif
   ENDIF
#  endif
!
! ------------------------------------------------------------------
!
!     ***  parameters nonlinear 3 wave interactions
!
! ============================================================
!
!  TRIad [trfac] [cutfr] [urcrit] [urslim]                                
!
! ============================================================
!
   TRI = .FALSE.
   CALL INLOGC('TRIAD',TRI,'UNC','F')
   
   IF(TRI)THEN
     ITRIAD = 1                                                        
     CALL INREAL('TRFAC' ,PTRIAD(1),'UNC',0.0)
     CALL INREAL('CUTFR', PTRIAD(2),'STA',2.5)                      
     CALL INREAL('URCRIT',PTRIAD(4),'UNC',0.0)                     
     CALL INREAL('URSLIM',PTRIAD(5),'STA',0.01)                    
!    reserve space in array Compda for Ursell number                   
     IF(JURSEL <= 1)THEN
       MCMVAR = MCMVAR + 1
       JURSEL = MCMVAR                                                 
     ENDIF
   ENDIF

!
! ------------------------------------------------------------------
!
!     LIM     setting parameters in conjunction with action limiter
!
! ============================================================
!
!   LIMiter [ursell] [qb]                                                 
!
! ============================================================
!
   LIM = .FALSE.
   CALL INLOGC('LIMITER',LIM,'UNC','F')
   
   IF(LIM)THEN                                            
     CALL INREAL('URSELL',PTRIAD(3),'UNC',0.)                          
     CALL INREAL('QB'    ,PNUMS(28),'UNC',0.)                          
   ENDIF                                                               

!
! ------------------------------------------------------------------
!
!     *** OBSTACLE   Definition of obstacles in comp grid. ***            
!
! ============================================================
!
!             |  TRANSm [trcoef]              |
! OBSTacle   <                                 >             &
!             |  DAM    [hgt] [alpha] [beta]  |
!                                                                         
!                        | -> RSPEC                    |                  
!      ( REFLec [reflc] <                              > )   &
!                        |    RDIFF [POWS][POWD][Kdif] |                  
!
!             LINe < [xp] [yp] >                                          
!
! ============================================================
!
   OBST = .FALSE.
   CALL INLOGC('OBSTACLE',OBST,'UNC','F')
   
   IF(OBST)THEN
!
     ALLOCATE(OBSTMP)                                                  
     OBSTMP%TRCOEF(1) = 0.                                             
     OBSTMP%TRCOEF(2) = 0.                                             
     OBSTMP%TRCOEF(3) = 0.                                             
     OBSTMP%RFCOEF(1) = 0.                                             
     OBSTMP%RFCOEF(2) = 0.                                             
     OBSTMP%RFCOEF(3) = 0.                                             
     OBSTMP%RFCOEF(4) = 0.                                             
     OBSTMP%RFCOEF(5) = 0.                                             
     OBSTMP%RFCOEF(6) = 0.                                             
     OBSTMP%RFCOEF(7) = 0.                                             
     OBSTMP%RFCOEF(8) = 0.                                             
!
!    data concerning transmission of energy over/through the obstacle
!
     CALL INCSTR('OBST_TYPE',OBST_TYPE,'UDF',' ')

     IF(OBST_TYPE == 'TRANS')THEN
       ITRAS = 0
       CALL INREAL('TRCOEF',TRCF,'UDF', 0.)
       IF((TRCF < 0.) .OR. (TRCF > 1.))THEN                        
         CALL MSGERR(3,'Transmission coeff. [trans] is not allowed ')  
         CALL MSGERR(3,'to be greater than 1 or smaller than 0!  ')    
       ENDIF                                                           
       OBSTMP%TRCOEF(1) = TRCF                                         
     ELSE IF(OBST_TYPE == 'DAM')THEN
       ITRAS = 1
       CALL INREAL('HGT'  , HGT , 'UDF', 0.  )                        
       CALL INREAL('ALPHA', OGAM, 'STA', 2.6 )                        
       CALL INREAL('BETA' , OBET, 'STA', 0.15)                        
       OBSTMP%TRCOEF(1) = HGT                                          
       OBSTMP%TRCOEF(2) = OGAM                                         
       OBSTMP%TRCOEF(3) = OBET                                         
     ELSE
!      if no transmission options are activated, there will be 0 transmission
       ITRAS = 0
       TRCF  = 0.                                                      
       OBSTMP%TRCOEF(1) = TRCF                                         
     ENDIF
     OBSTMP%TRTYPE = ITRAS                                             
!                                                                      
!    data on reflection by the obstacle, activated by keyword REFL     
!
     REFL = .FALSE.
     CALL INLOGC('REFLEC',REFL,'DFL','T')

     IF(REFL)THEN                                         
       LREF = 1                                                     
       IF(.NOT.FULCIR)THEN                                           
         CALL MSGERR(3,'Reflections will only be calculated if     ')  
         CALL MSGERR(3,'the spectral directions cover the full     ')  
         CALL MSGERR(3,'circle.                                    ')  
       ENDIF                                                           
       CALL INREAL('REFLC',REF0,'STA',1.)                          
       DUM = REF0*REF0+TRCF*TRCF                                       
       IF(ITRAS == 0 .AND. ((DUM < 0.) .OR. (DUM > 1.)))THEN         
         CALL MSGERR(3,'Kt^2 + Kr^2 > 1 or Kt^2 + Kr^2 < 0! ')        
       ENDIF                                                          
       OBSTMP%RFCOEF(1) = REF0                                        
!                                                                     
       CALL INCSTR('RDIFF',RDIFF,'UDF',' ')

       IF(RDIFF == 'RDIFF')THEN                                       
         LREFDIFF = 1                                                
         CALL INREAL('POWS',POWS,'UDF',1.)                       
         DUM = MOD(POWS,1.)                                          
         IF(POWS < 0.)THEN                                           
           CALL MSGERR(3,'Power POWS is not a positive number! ')    
         ENDIF                                                       
         IF(DUM /= 0.)THEN                                         
           CALL MSGERR(3,'Power POWS is not an integer number! ')     
         ENDIF                                                         
         OBSTMP%RFCOEF(2) = POWS                                       
!                                                                      
         CALL INREAL('POWD',POWD,'STA',1.)                         
         DUM = MOD(POWD,1.)                                            
         IF(POWD < 0.)THEN                                           
           CALL MSGERR(3,'Power POWD is not a positive number! ')     
         ENDIF                                                        
         IF(DUM /= 0.)THEN                                         
           CALL MSGERR(3,'Power POWD is not an integer number! ')      
         ENDIF                                                         
         OBSTMP%RFCOEF(3) = POWD                                       
         CALL INREAL ('KDIF', KDIF, 'STA', 0.)                         
         IF((KDIF < 0.) .OR. (KDIF > 1.))THEN                        
           CALL MSGERR(3,'KDIF > 1 or KDIF < 0! ')                     
         ENDIF                                                         
         OBSTMP%RFCOEF(4) = KDIF                                       
       ELSE IF(RDIFF == 'RSPEC')THEN                                                           
         LREFDIFF = 0                                                  
       ENDIF                                                           
       OBSTMP%RFTYP2 = LREFDIFF                                        
!
       RFD = .FALSE.
       CALL INLOGC('RFD',RFD,'DFL','T')

       IF(RFD)THEN                                         
         LRFRD = 1                                                     
         CALL INREAL('FD1',FD1,'REQ', 0.565)                        
         CALL INREAL('FD2', FD2, 'REQ', 0.75)                         
         CALL INREAL('FD3', FD3, 'REQ', -21.)                         
         CALL INREAL('FD4', FD4, 'REQ', 0.066)                        
         IF(FD3 >= 0.)THEN                                         
           CALL MSGERR(3,'Positive frequency dependency! ')            
         ENDIF                                                         
         OBSTMP%RFCOEF(5) = FD1                                        
         OBSTMP%RFCOEF(6) = FD2                                        
         OBSTMP%RFCOEF(7) = FD3                                        
         OBSTMP%RFCOEF(8) = FD4                                        
       ELSE                                                            
         LRFRD = 0                                                     
       ENDIF                                                           
       OBSTMP%RFTYP3 = LRFRD                                           

       IF((REF0 < 0.) .OR. (REF0 > 1.))THEN                        
         CALL MSGERR(3,'Reflection coeff. [reflc] is not allowed ')    
         CALL MSGERR(3,'to be greater than 1 or smaller than 0!  ')    
       ENDIF                                                           
     ELSE
!      if there is no keyword REFL, there will be no reflection        
       LREF = 0
     ENDIF                                                             
     OBSTMP%RFTYP1 = LREF                                              
!
!    location of obstacles
!
!    *** NUMCOR : Number of corners ***
     NUMCOR = 0
     CALL INLOGC('LINE',LINE,'REQ','T')
     IF(LINE)THEN
       FRST%X = 0.                                                     
       FRST%Y = 0.                                                     
       NULLIFY(FRST%NEXTXY)                                            
       CURR => FRST                                                    
       DO
!        read coordinates of one corner point of the obstacle
!JQI         CALL READXY ('XP', 'YP', XP, YP, 'REP', -1.E10, -1.E10)
!JQI         IF(XP < -.9E10) GOTO 101
!JQI         NUMCOR = NUMCOR + 1
!JQI         ALLOCATE(TMP)                                                 
!JQI         TMP%X = XP                                                    
!JQI         TMP%Y = YP                                                    
!JQI         NULLIFY(TMP%NEXTXY)                                           
!JQI         CURR%NEXTXY => TMP                                            
!JQI         CURR => TMP                                                   
       ENDDO
     ENDIF
101  OBSTMP%NCRPTS = NUMCOR                                            
!    store coordinates for corner points in array of obstacle data     
     ALLOCATE(OBSTMP%XCRP(NUMCOR))                                     
     ALLOCATE(OBSTMP%YCRP(NUMCOR))                                     
     CURR => FRST%NEXTXY                                               
     DO JJ = 1, NUMCOR                                                 
       OBSTMP%XCRP(JJ) = CURR%X                                      
       OBSTMP%YCRP(JJ) = CURR%Y                                      
       CURR => CURR%NEXTXY                                           
     END DO                                                            
     DEALLOCATE(TMP)                                                   
     NULLIFY(OBSTMP%NEXTOBST)                                          
     IF(NUMCOR <= 1)THEN
!JQI       CALL MSGERR(1,'No corner points for obstacle were found')
     ELSE
!      *** NUMOBS : Number of obstacles ***
       NUMOBS = NUMOBS + 1
       IF(.NOT.LOBST)THEN                                          
         FOBSTAC = OBSTMP                                             
         COBST => FOBSTAC                                             
         LOBST = .TRUE.                                               
       ELSE                                                           
         COBST%NEXTOBST => OBSTMP                                     
         COBST => OBSTMP                                              
       END IF                                                         
     ENDIF
!
   ENDIF

!
!  ------------------------------------------------------------------
!  SETUP     include wave-induced set-up in SWAN calculation
!
! ============================================================
!
!   SETUP [supcor]                                                        
!
! ============================================================
!
   SETUP = .FALSE.
   LSETUP = 0                                                        
   CALL INLOGC('SETUP',SETUP,'UNC','F')
   
   IF(SETUP)THEN  
     IF(KSPHER /= 0)                          &
       CALL MSGERR (2,'setup will not be compute correctly with spherical coordinates')  
     LSETUP = 1                                                        
     CALL INREAL ('SUPCOR',PSETUP(2),'UNC',0.)                      
!    *** set pointers for setup and saved depth in array COMPDA ***    
     JSETUP = MCMVAR + 1                                               
     JDPSAV = MCMVAR + 2                                               
     MCMVAR = MCMVAR + 2                                               
   ENDIF  

!
!  ------------------------------------------------------------------  
!
!  DIFFRac   include diffraction approximation                         
!                                                                         
! =============================================================           
!                                                                         
!   DIFFRac  [idiffr]  [smpar]  [smnum]  [cgmod]                          
!                                                                         
! =============================================================           
!                                                                         
   DIFFR = .FALSE.
   CALL INLOGC('DIFFRAC',DIFFR,'UNC','F')
   
   IF(DIFFR)THEN                                          
     CALL ININTG('IDIFFR', IDIFFR, 'STA', 1)                          
     CALL INREAL('SMPAR', PDIFFR(1), 'STA', 0.0)                      
     CALL INREAL('SMNUM', PDIFFR(2), 'STA', 0.)                       
     CALL INREAL('CGMOD', PDIFFR(3), 'STA', 1.)                       
   ENDIF                                                               

!
!  ------------------------------------------------------------------
!
!     OFF      switching standard options off
!
! ============================================================
!
!        |  BREaking
!        |
!        |  WCAPping
!        |
!        |  REFrac
! OFF   <
!        |  FSHift
!        |
!        |  QUADrupl
!        |
!        |  WINDGrowth                                                    
!        |
!        |  BNDCHK                                                        
!
! ============================================================
!
   OFF = .FALSE.
   CALL INLOGC('OFF',OFF,'UNC','F')
   
   IF(OFF)THEN
     CALL INCSTR('OFF_REF',OFF_TYPE,'STA','YES')
     IF(OFF_TYPE == 'YES') IREFR = 0
     CALL INCSTR('OFF_FSH',OFF_TYPE,'STA','YES')
     IF(OFF_TYPE == 'YES') ITFRE = 0
     CALL INCSTR('OFF_BRE',OFF_TYPE,'STA','YES')
     IF(OFF_TYPE == 'YES') ISURF = 0
     CALL INCSTR('OFF_WCAP',OFF_TYPE,'STA','YES')
     IF(OFF_TYPE == 'YES') IWCAP = 0
     CALL INCSTR('OFF_QUAD',OFF_TYPE,'STA','YES')
     IF(OFF_TYPE == 'YES') IQUAD = 0
     CALL INCSTR('OFF_WINDG',OFF_TYPE,'STA','YES')
     IF(OFF_TYPE == 'YES') IWIND = 0
     CALL INCSTR('OFF_BNDCHK',OFF_TYPE,'STA','YES')
     IF(OFF_TYPE == 'YES') BNDCHK = .FALSE.        !switch off checking of Hs on boundary  
     CALL INCSTR('OFF_RESCALE',OFF_TYPE,'STA','YES')
     IF(OFF_TYPE == 'YES') BRESCL = .FALSE.        !switch off rescaling solution
   END IF
!
!  ------------------------------------------------------------------
! =======================================================================
!
!            |    BSBT                      |
!   PROP    <                     | SEC |   |
!            |    GSE  [waveage] <  MIN  >  |
!            |                    |  HR |   |
!            |                    | DAY |   |
!
!            FLUXLIM
!
! =======================================================================
!
   PROP = .FALSE.
   CALL INLOGC('PROP',PROP,'UNC','F')
   IF(PROP)THEN                                           
     CALL INCSTR('PROP_CHOICE', PROP_CHOICE,'UDF',' ')

     IF(PROP_CHOICE == 'BSBT')THEN                                         
       PROPSN = 1                                                    
       PROPSS = 1                                                    
     ELSE IF(PROP_CHOICE == 'GSE')THEN                                     
       IF(PROPSN /= 3)THEN                                           
         CALL MSGERR(2,'Anti-GSE only allowed for S&L scheme.')    
       ENDIF                                                         
       CALL ININTV('WAVEAGE', WAVAGE, 'STA', 0.)                     
     ENDIF                                                           

     FLUXLIM = .FALSE.
     CALL INLOGC('FLUXLIM',FLUXLIM,'DFL','T')
     
     IF(FLUXLIM)THEN                                      
       PROPFL   = 1                                                  
       PNUMS(6) = 0.                                                 
     END IF                                                          
   ENDIF                                                             

!
!     ------------------------------------------------------------------
!
!     NUM       setting numerical parameters for schemes, solvers, etc.
!
! =====================================================================
!
!             | -> ACCUR [drel] [dhoval] [dtoval] [npnts] |               
!   NUMeric (<                                             >          &   
!             |    STOPC [dabs] [drel] [curvat] [npnts]   |               
!
!                    | -> STAT  [mxitst] [alfa] |                         
!                   <                            >  [limiter]   )     &   
!                    | NONSTat  [mxitns]        |
!
!           ( DIRimpl [cdd] [cdlim]                                )  &
!
!
!           | -> SIGIMpl [css] [eps2] [outp] [niter]               |
!          (<                                                      >) &
!           |    SIGEXpl [cfl]                 (NOT documented)    |
!           |                                                      |
!           |    FIL     [diffc]               (NOT documented)    |
!
!
!           ( SETUP [eps2] [outp] [niter]                          )      
!
! =====================================================================
!
   NUM = .FALSE.
   CALL INLOGC('NUMERIC',NUM,'UNC','F')
   
   IF(NUM)THEN
     CALL INCSTR('NUM_CHOICE',NUM_CHOICE,'UDF',' ')

!    *** accuracy and criterion to terminate the iteration ***
     IF(NUM_CHOICE == 'ACCUR')THEN
       PNUMS(21) = 0.                                                  
       CALL INREAL ('DREL'   , PNUMS(1) , 'UNC', 0.)
       CALL INREAL ('DHOVAL' , PNUMS(15), 'UNC', 0.)                   
       CALL INREAL ('DTOVAL' , PNUMS(16), 'UNC', 0.)                   
       CALL INREAL ('NPNTS'  , PNUMS(4) , 'UNC', 0.)

       CALL INCSTR('STAT',STATL,'UDF',' ')

       IF(STATL == 'STAT')THEN                                       
         CALL ININTG ('MXITST' , MXITST   , 'REQ', 0 )                 
         CALL INREAL ('ALFA'   , PNUMS(30), 'UNC', 0.)                 
       ELSE IF(STATL == 'ITERMX')THEN                                
         CALL ININTG ('MXITST' , MXITST   , 'REQ', 0 )                 
         MXITNS = MXITST                                               
         CALL MSGERR (1, '[itermx] is replaced; see user manual')      
       ELSE IF(STATL == 'NONST')THEN                                 
         CALL ININTG ('MXITNS' , MXITNS   , 'REQ', 0 )                 
       ENDIF
       CALL INREAL ('LIMITER', PNUMS(20), 'UNC', 0.)                   
     ELSE IF(NUM_CHOICE == 'STOPC')THEN                                   
       PNUMS(21) = 1.                                                  
       CALL INREAL ('DABS'   , PNUMS(2) , 'STA', 0.00)                 
       CALL INREAL ('DREL'   , PNUMS(1) , 'STA', 0.01)                 
       CALL INREAL ('CURVAT' , PNUMS(15), 'STA', 0.005)                
       CALL INREAL ('NPNTS'  , PNUMS(4) , 'UNC', 0.)                   

       CALL INCSTR('STAT',STATL,'UDF',' ')

       IF(STATL == 'STAT')THEN                                       
         CALL ININTG ('MXITST' , MXITST   , 'STA', 50)                 
         CALL INREAL ('ALFA'   , PNUMS(30), 'UNC', 0.)                 
       ELSE IF(STATL == 'ITERMX')THEN                                
         CALL ININTG ('MXITST' , MXITST   , 'REQ', 0 )                 
         MXITNS = MXITST                                               
         CALL MSGERR (1, '[itermx] is replaced; see user manual')      
       ELSE IF(STATL == 'NONST')THEN                                 
         CALL ININTG ('MXITNS' , MXITNS   , 'REQ', 0 )                 
       ENDIF                                                           
       CALL INREAL ('LIMITER', PNUMS(20), 'UNC', 0.)                   
     END IF
!
!    *** numerical scheme in directional space (standard  ***
!    *** option implicit scheme)                          ***
!
     DIR = .FALSE.
     CALL INLOGC('DIRIMPL',DIR,'DFL','T')

     IF(DIR)THEN
       CALL INREAL ('CDD'    , PNUMS(6) , 'UNC', 0.)
       CALL INREAL ('CDLIM'  , PNUMS(17), 'UNC', 0.)                   
       IF(PNUMS(17) < 0.) IREFR = 1                                  
       IF(PNUMS(17) > 0.) IREFR = -1                                 
       IF(EQREAL(PNUMS(17),0.))THEN                                  
         IREFR = 0                                                     
         CALL MSGERR(0, 'Refraction deactivated')                      
       ENDIF
     ENDIF
!
!    *** numerical scheme in frequency space :                  ***
!    *** 1) Fully implicit scheme in frequency space (iterative ***
!    ***    SIP solver)                                         ***
!    *** 2) Explicit scheme in frequency space:                 ***
!    ***    Energy is removed from the spectrum based on a CFL  ***
!    ***    criterion                                           ***
!    *** 3) Explicit scheme in frequency space:                 ***
!    ***    No CFL limitation near blocking point --> unstable  ***
!    ***    integration. Therefore a filter is applied with a   ***
!    ***    diffusion coeff.                                    ***
!
     CALL INCSTR('SIGIMPL',SIGIM,'UDF',' ')

     IF(SIGIM == 'SIGIM')THEN                     
!      *** implicit solver ***
!      ***  This is the default option   PNUMS(8) = 1.  ***
       PNUMS(8) = 1.
       CALL INREAL ('CSS'   , PNUMS(7), 'UNC', 0.)
       CALL INREAL ('EPS1'  , PNUMS(11) , 'UNC', 0.)
       CALL INREAL ('EPS2'  , PNUMS(12) , 'UNC', 0.)
       CALL INREAL ('OUTP'  , PNUMS(13) , 'UNC', 0.)
       CALL INREAL ('NITER' , PNUMS(14), 'UNC', 0.)
     ELSE IF(SIGIM == 'SIGEX')THEN                 
!
!      *** 2) explicit scheme ***
       PNUMS(8) = 2.
       CALL INREAL ('CFL'  , PNUMS(19), 'UNC', 0.)
!
     ELSE IF(SIGIM == 'FIL')THEN
!      *** 3) explicit scheme -> filter the spectrum  ***
       PNUMS(8) = 3.
       CALL INREAL ('DIFFC'  , PNUMS(9), 'UNC', 0.)
!
     END IF

     SETUP = .FALSE.
     CALL INLOGC('SETUP',SETUP,'DFL','T')

     IF(SETUP)THEN                                         
!      *** iterative solver        ***                                 
!      *** Settings for the setup  ***                                 
       CALL INREAL ('EPS2' , PNUMS(23), 'UNC', 0.)                     
       CALL INREAL ('OUTP' , PNUMS(24), 'UNC', 0.)                     
       CALL INREAL ('NITER', PNUMS(25), 'UNC', 0.)                     
     ENDIF                                                             
   ENDIF
!
!     ------------------------------------------------------------------
!
!     **  Command COMPUTE   **                                       
!
!       ==============================================================
!
!                  |  STATionary  [time]                      |
!       COMPute ( <                                            > )        
!                  |                    | -> Sec  |           |
!                  |  ([tbegc] [deltc] <     MIn   > [tendc]) |
!                                       |    HR   |
!                                       |    DAy  |
!
!       ==============================================================
!

   CALL INCSTR('COMPUT',COMPUT,'REQ',' ')
   IF(COMPUT == 'COMP')THEN
     RUNMADE = .TRUE.                                                  
!     CALL INLOGC('COM_STAT',COM_STAT,'REQ','T')
     COM_STAT = .FALSE.
     IF(NSTATM <= 0 .OR. COM_STAT)THEN                         
       IF(NSTATM == -1) NSTATM = 0                                    
       IF(NSTATM > 0) CALL INCTIM (ITMOPT,'TIME',TINIC,'REQ',0.)     
       IF(TINIC < TIMCO)THEN                                      
         CALL MSGERR (2, '[time] before current time')                 
         TINIC = TIMCO                                                 
       ENDIF                                                           
       TFINC = TINIC                                                   
       TIMCO = TINIC                                                   
       DTW = 1.E10                                                      
       RDTIM = 0.                                                      
       NSTATC = 0                                                      
       MTC = 1                                                         
     ELSE

       TINIC = IFLBEG(IGRD)
       IF(TINIC < TIMCO)THEN                                      
         CALL MSGERR (2, 'start time [tbegc] before current time')     
         TINIC = TIMCO                                                 
       ENDIF                                                           
       CALL INREAL('NS_DELTC',DTW,'REQ',1)
       
       CALL INCSTR('TIME_UNIT',UNIT,'UNC','SECOND')
       IF(TRIM(UNIT) == 'DAY' .OR. TRIM(UNIT) == 'day')THEN
         DTW = DTW*24.0*3600.0
       ELSE IF(TRIM(UNIT) == 'HOUR' .OR. TRIM(UNIT) == 'hour')THEN
         DTW = DTW*3600.0
       ELSE IF(TRIM(UNIT) == 'MINUTE' .OR. TRIM(UNIT) == 'minute')THEN
         DTW = DTW*60.0
       ELSE IF(TRIM(UNIT) /= 'SECOND' .AND. TRIM(UNIT) /= 'second')THEN
         CALL MSGERR(3, 'the unit of time step must be DAY,HOUR,MINUTE &
	                 or SECOND')	 
       END IF	 	 	 
       
       NSTATC = 1 
       IMDTW = SECONDS2TIME(DTW)
       TFINC = SECONDS(EndTime)                                                     
!
!      *** tfinc must be greater than tinic **
       DIFF = TFINC - TINIC
       
       if(msr) PRINT*,'TINIC=',TINIC,'TFINC=',TFINC,'DIFF=',DIFF,'MYID=',MYID
       IF(DIFF <= 0.) CALL MSGERR (3,                                &
              'start time [tbegc] greater or equal end time [tendc]')
!
!        **The number of computational steps is calculated
       RDTIM = 1./DTW
       MTC = NINT ((TFINC - TINIC)/DTW)                                 
       IF(MOD(TINIC-TFINC,DTW) > 0.01*DTW .AND.                        &
          MOD(TINIC-TFINC,DTW) < 0.99*DTW)                             &
	       CALL MSGERR (1,                                       &
	       'DTW is not a fraction of the computational period')
       TIMCO = TINIC
     ENDIF
!JQI     IF(NSTATM > 0) CHTIME = DTTIWR(ITMOPT, TIMCO)                   
     NCOMPT = NCOMPT + 1                                               
!JQI     IF(NCOMPT > 50) CALL MSGERR (2,                                 &
!JQI                   'No more than 50 COMPUTE commands are allowed')  
     RCOMPT(NCOMPT,1) = REAL(NSTATC)                                   
     RCOMPT(NCOMPT,2) = REAL(MTC)                                      
     RCOMPT(NCOMPT,3) = TFINC                                          
     RCOMPT(NCOMPT,4) = TINIC                                          
     RCOMPT(NCOMPT,5) = DTW                                             
!    set ITERMX equal to MXITST in case of stationary computations
!    and to MXITNS otherwise
     IF(NSTATC == 0)THEN
       ITERMX = MXITST                                                 
     ELSE
       ITERMX = MXITNS                                                 
     ENDIF
!     RETURN                                                            

     CALL INREAL('SOURCE_TERM_DTMAX',SOURCE_DTMAX,'REQ',1)
     CALL INREAL('SOURCE_TERM_DTMIN',SOURCE_DTMIN,'REQ',1)
   ENDIF

   IF(MSR)PRINT*,'AT THE END OF SWREAD, MCGRD = ',MCGRD
!   STOP
   RETURN
   END SUBROUTINE SWREAD
 
!***********************************************************************
!                                                                      *
   SUBROUTINE SSFILL (SPCSIG, SPCDIR)                                  
!                                                                      *
!***********************************************************************
!
!     Discretisation in frequency (sigma) and direction (theta)
!
!  Method
!
!     Create a logarithmic distribution in frequency between lowest and highest
!     frequencies and store them in the array SPCSIG.
!
!     Create a distribution in direction that does not coincide with the orientation
!     of the computational grid, calculate some geometric derivatives and store
!     it in the array SPCDIR
!
!***********************************************************************
!
   USE OCPCOMM1                                                        
   USE OCPCOMM2                                                        
   USE OCPCOMM3                                                        
   USE OCPCOMM4                                                        
   USE SWCOMM1                                                         
   USE SWCOMM2                                                         
   USE SWCOMM3                                                         
   USE SWCOMM4                                                         
   USE M_WCAP  
      
   IMPLICIT NONE                                                        

   REAL    SPCDIR(MDC,6)                                               
   REAL    SPCSIG(MSC)     
   
   REAL    :: SFAC,OLDDIR
   INTEGER :: ID,IS                                            

!
!  Allocate arrays in module M_WCAP                                    
!
   IF (.NOT.ALLOCATED(SIGPOW)) ALLOCATE (SIGPOW(MSC,6))                
!
!
!  distribution of spectral frequencies
!
!  FRINTF is the frequency integration factor (=df/f)                  
   FRINTF = ALOG(SHIG/SLOW) / FLOAT(MSC-1)                             
   SFAC   = EXP(FRINTF)                                                
   FRINTH = SQRT(SFAC)
!  determine spectral frequencies (logarithmic distribution)
   SPCSIG(1) = SLOW                                                    
   DO IS = 2, MSC
     SPCSIG(IS) = SPCSIG(IS-1) * SFAC                                 
   END DO  
!
!  Calculate powers of sigma and store in global array
!
   SIGPOW(:,1) = SPCSIG                                                
   SIGPOW(:,2) = SPCSIG**2                                             
   SIGPOW(:,3) = SPCSIG * SIGPOW(:,2)                                  
   SIGPOW(:,4) = SPCSIG * SIGPOW(:,3)                                  
   SIGPOW(:,5) = SPCSIG * SIGPOW(:,4)                                  
   SIGPOW(:,6) = SPCSIG * SIGPOW(:,5)                                  
!
!  distribution of spectral directions
!
   DO ID = 1, MDC
     SPCDIR(ID,1) = SPDIR1 + FLOAT(ID-1)*DDIR                         
!    if a direction coincides with a direction of the (regular)       
!    computational grid it is slightly changed
     IF(OPTG == 1)THEN                                              
       IF(ABS(MODULO(SPCDIR(ID,1)-ALPC+0.25*PI_W ,0.5*PI_W)-0.25*PI_W)  &
	      < 1.E-6)THEN
         OLDDIR = SPCDIR(ID,1)
         SPCDIR(ID,1) = OLDDIR + 2.E-6                                
       ENDIF
     ENDIF                                                            
     SPCDIR(ID,2) = COS(SPCDIR(ID,1))                                 
     SPCDIR(ID,3) = SIN(SPCDIR(ID,1))                                 
     SPCDIR(ID,4) = SPCDIR(ID,2) **2                                  
     SPCDIR(ID,5) = SPCDIR(ID,2) * SPCDIR(ID,3)                       
     SPCDIR(ID,6) = SPCDIR(ID,3) **2                                  
   END DO  
!
   RETURN
   END SUBROUTINE SSFILL
#  endif
