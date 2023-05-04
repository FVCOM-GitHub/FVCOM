!************************************************************************
!                                                                       *
   SUBROUTINE SWBOUN ( )         
#  if defined (WAVE_CURRENT_INTERACTION)
!                                                                       *
!************************************************************************
!
!     Reading and processing BOUNDARY command
!
!************************************************************************
!
   USE OCPCOMM2                                                        
   USE OCPCOMM4                                                        
   USE SWCOMM1                                                         
   USE SWCOMM2                                                         
   USE SWCOMM3                                                         
   USE M_BNDSPEC                                                       
   USE MOD_UTILS
   USE VARS_WAVE,ONLY : NESTING,OBC_HS,OBC_DIR,OBC_TPEAK
   USE ALL_VARS,ONLY : MSR,PAR,SERIAL
   USE BCS,ONLY : IOBCN_GL_W,IOBCN_W,ELO_TM,OBC_NTIME,I_OBC_GL_W
   USE MOD_PAR, ONLY : NLID
   USE MOD_NESTING, ONLY : NESTING_ON_WAVE

   IMPLICIT NONE

   INTEGER   IENT,KOUNTR,IX1,IY1,IX2,IY2
   INTEGER   MM,IX,IY,ISIDM,ISIDE,KC,KC2,KC1,IX3,IY3,MP
   INTEGER   IP,II,NBSPSS,NFSEQ,IKO,IKO2,IBSPC1,IBSPC2

   REAL      CRDP, CRDM, SOMX, SOMY
   REAL      XP,YP,XC,YC,RR,DIRSI,COSDIR,SINDIR,DIRSID,DIRREF
   REAL      RLEN1,RDIST,RLEN2,XC1,YC1,XC2,YC2,W1

   LOGICAL   KEYWIS, LOCGRI, CCW, BPARF, BOUNPT,DONALL
   LOGICAL   LFRST1, LFRST2, LFRST3                                    

   INTEGER   UPLO, NUMP

   LOGICAL, SAVE :: LBFILS = .FALSE.                                   
   LOGICAL, SAVE :: LBS    = .FALSE.                                   
   LOGICAL, SAVE :: LBGP   = .FALSE.                                   

   TYPE(BSPCDAT), POINTER :: BFLTMP                                    
   TYPE(BSPCDAT), SAVE, POINTER :: CUBFL                               

   TYPE(BSDAT), POINTER :: BSTMP                                       
   TYPE(BSDAT), SAVE, POINTER :: CUBS                                  

   TYPE(BGPDAT), POINTER :: BGPTMP                                     

   TYPE XYPT                                                           
     INTEGER             :: JX, JY
     TYPE(XYPT), POINTER :: NEXTXY
   END TYPE XYPT

   TYPE(XYPT), TARGET  :: FRST                                         
   TYPE(XYPT), POINTER :: CURR, TMP                                    

   LOGICAL STPNOW                                                      
   LOGICAL EQREAL
   LOGICAL CHECK
   CHARACTER(LEN=12) BOUND_CHOICE,SHAPESPEC,CHAR_WAVE_PERIOD,DSPR
   CHARACTER(LEN=120) :: NESTING_FILE,NCFILE
   INTEGER   :: NTIME,NOBC,IERR,J,NCNT
   REAL,  ALLOCATABLE :: HSC1_TMP(:,:),DIR1_TMP(:,:),TPEAK_TMP(:,:)
   REAL,  ALLOCATABLE :: TIMES(:)
   INTEGER :: I
   REAL :: NTMP,F_TMP
!
   CALL INCSTR('BOUND_CHOICE',BOUND_CHOICE,'UNC',' ')

   IF(BOUND_CHOICE == 'UNIFORM' .OR. BOUND_CHOICE == 'CONSTANT')THEN
!
!    specification of the spectral shape
!
! =========================================================================
!
!                      |  JONswap  [gamma]  |
!                      |                    |    | -> PEAK |
!  BOUNdspec  SHAPe   <   PM                 >  <           >   &
!                      |                    |    | MEAN    |
!                      |  GAUSs  [sigfr]    |
!                      |                    |
!                      |  BIN               |
!
!                     | DEGRees   |
!             DSPR   <             >
!                     | -> POWer  |
!
! =========================================================================
!
     IF(NESTING_ON_WAVE)THEN
       IF(MSR)WRITE(IPT,*) "BOUND_CHOICE = ",BOUND_CHOICE,"NESTING_ON_WAVE = ",NESTING_ON_WAVE
       CALL FATAL_ERROR("The parameter NESTING_ON_WAVE in ***_run.nml should be .FALSE.",  &
                        "or BOUND_CHOICE in INPUT should be NESTING.")
     END IF	
     CALL INCSTR('SHAPESPEC', SHAPESPEC,'UNC',' ')
     IF(SHAPESPEC == 'JON')THEN
       FSHAPE = 2
       CALL INREAL('GAMMA', PSHAPE(1), 'STA', 3.3)                    
     ELSE IF(SHAPESPEC == 'BIN')THEN
       FSHAPE = 3
     ELSE IF(SHAPESPEC == 'PM')THEN
       FSHAPE = 1
     ELSE IF(SHAPESPEC == 'GAUS')THEN
       FSHAPE = 4
       CALL INREAL('SIGFR', SIGMAG, 'STA', 0.01)
!      convert from Hz to rad/s:
       PSHAPE(2) = PI2_W * SIGMAG                                        
     ENDIF
!    PEAK or MEAN frequency
     CALL INCSTR('CHAR_WAVE_PERIOD',CHAR_WAVE_PERIOD,'UNC',' ')
     IF(CHAR_WAVE_PERIOD == 'MEAN')THEN
       FSHAPE = -FSHAPE
     ENDIF
!    directional distribution given by DEGR or by POWER
     CALL INCSTR('DSPR',DSPR,'UNC',' ')
     IF(DSPR == 'DEGREES')THEN
       DSHAPE = 1
     ELSE IF(DSPR == 'POW')THEN
       DSHAPE = 2
     ELSE
       STOP  
     ENDIF

!
     CALL INREAL ('HSIG',  SPPARM(1), 'REQ', 0.)
     CALL INREAL ('PER', SPPARM(2), 'REQ', 0.)
     CALL INREAL ('DIR', SPPARM(3), 'REQ', 0.)
     IF (DSHAPE == 1) THEN
       CALL INREAL ('DD',  SPPARM(4), 'UNC', 30.)
       IF (SPPARM(4) > 360. .OR. SPPARM(4) < 0.) THEN               
         CALL MSGERR (2,'Directional spreading is less than '//    &
                        '0 or larger than 360 degrees') 
       END IF                                                      
     ELSE
       CALL INREAL ('DD',  SPPARM(4), 'UNC', 2.)
       IF (SPPARM(4) <= 0.) THEN             
         CALL MSGERR (2,'Power of cosine is less or equal to zero')   
       END IF                                                   
       IF (SPPARM(4)*DDIR**2/2. > 1.) THEN                     
         CALL MSGERR (2,'distribution too narrow to be represented properly')
         WRITE (PRINTF, 142) SQRT(2./SPPARM(4))*180./PI_W         
142      FORMAT (' Advise: choose Dtheta < ', F8.3, ' degr')      
       END IF                                     
     ENDIF
	 
   ELSE IF(BOUND_CHOICE == 'WAMN')THEN
!
!
   ELSE IF(BOUND_CHOICE == 'WW3' .OR. BOUND_CHOICE == 'WWIII')THEN                     
!
   ELSE IF(BOUND_CHOICE == 'NESTING')THEN
!
!    Nesting SWAN model in larger SWAN model
! ==========================================
!                                | -> CLOS |                              
!     BOUNdnest1  NEST 'fname'  <           >
!                                |  OPEN   |                              
! ==========================================

!     NBFILS = NBFILS + 1
!     ALLOCATE(BFLTMP)                                                  
!     CALL INCSTR ('BOUND_NEST_IDX',FILENM,'REQ', ' ')
!     CALL BCFILE (FILENM, 'NEST')                                              
!     NULLIFY(BFLTMP%NEXTBSPC)                                          
!     IF ( .NOT.LBFILS ) THEN                                           
!       FBNDFIL = BFLTMP                                               
!       CUBFL => FBNDFIL                                               
!       LBFILS = .TRUE.                                                
!     ELSE                                                              
!       CUBFL%NEXTBSPC => BFLTMP                                       
!       CUBFL => BFLTMP                                                
!     END IF 
!=====================================================================!
!    reading the boundary forcing from nesting
!=====================================================================!
     NESTING = .TRUE.
     
     IF(.NOT. NESTING_ON_WAVE)THEN
       IF(MSR)WRITE(IPT,*) "BOUND_CHOICE = ",BOUND_CHOICE,"NESTING_ON_WAVE = ",NESTING_ON_WAVE
       CALL FATAL_ERROR("The parameter BOUND_CHOICE in INPUT should be UNIFORM ",&
        "or NESTING_ON_WAVE in ***_run.nml should be .TRUE.")
     END IF	
!JQI     CALL INCSTR('NESTING_FILE',NESTING_FILE,'REQ',' ')
!JQI     INQUIRE(FILE=TRIM(NESTING_FILE),EXIST=CHECK)   
!JQI      IF(CHECK)THEN
!JQI     ELSE
!JQI       WRITE(PRINTF,*) TRIM(NESTING_FILE), ' DOES NOT EXIT.'
!JQI       CALL PSTOP
 !JQI    END IF

!JQI     IF(MSR)THEN
!JQI       NCFILE = "./"//TRIM(NESTING_FILE)
!JQI       OPEN(1,FILE=NCFILE)
!JQI       NTIME = 0
!JQI       READ(1,*)
!JQI       DO I = 1, 1000
!JQI        READ(1,*,END=100) NTMP
!JQI	NTIME = NTIME + 1
!JQI	DO J = 1,IOBCN_GL_W
!JQI	 READ(1,*)
!JQI	END DO 
!JQI       END DO
!JQI100    CONTINUE
!JQI       REWIND(1)
!JQI      END IF
      
!JQI# if defined (MULTIPROCESSOR)
!JQI      IF(PAR)CALL MPI_BCAST(NTIME,1,MPI_INTEGER,0,MPI_FVCOM_GROUP,IERR)
!JQI# endif
!JQI      ALLOCATE(HSC1_TMP(IOBCN_GL_W,NTIME))
!JQI      ALLOCATE(DIR1_TMP(IOBCN_GL_W,NTIME))
!JQI      ALLOCATE(TPEAK_TMP(IOBCN_GL_W,NTIME))
!JQI      ALLOCATE(TIMES(NTIME))

!JQI      IF(MSR)THEN 
!JQI       READ(1,*)
!JQI       DO I=1,NTIME
!JQI        READ(1,*) TIMES(I)
!JQI	DO J = 1,IOBCN_GL_W
!JQI	 READ(1,*) HSC1_TMP(J,I),TPEAK_TMP(J,I),DIR1_TMP(J,I)
!JQI	END DO
 !JQI      END DO	 
!JQI      END IF  
      
!JQI      CLOSE(1)      	

!!# if defined (NETCDF_IO)


!!     IF(MSR)THEN
!!       NCFILE = "./"//TRIM(NESTING_FILE)
!!       CALL NEST_READ_TIME(NCFILE,NTIME,NOBC)
!!       IF(NOBC/=IOBCN_GL_W)THEN
!!          IF(MSR) WRITE(PRINTF,*)'obc number in INP_OBC_NAME i&
!!               &s different from nesting input file :' &
!!               &,trim(NESTING_FILE),nobc ,'/=',IOBCN_GL_W
!!          CALL PSTOP
!!       END IF

!!     END IF
!!# if defined (MULTIPROCESSOR)
!!     IF(PAR)CALL MPI_BCAST(NTIME,1,MPI_INTEGER,0,MPI_FVCOM_GROUP,IERR)
!!# endif
!!     ALLOCATE(HSC1_TMP(IOBCN_GL_W,NTIME))
!!     ALLOCATE(DIR1_TMP(IOBCN_GL_W,NTIME))
!!     ALLOCATE(TPEAK_TMP(IOBCN_GL_W,NTIME))
!!     ALLOCATE(TIMES(NTIME))
!!     IF(MSR)THEN
!!       CALL NCD_READ_NEST(NCFILE,NTIME,IOBCN_GL_W,TIMES,HSC1_TMP,DIR1_TMP,TPEAK_TMP)  
!!     END IF
!!# endif

!JQI#    if defined (MULTIPROCESSOR)
!JQI     IF(PAR)CALL MPI_BCAST(HSC1_TMP,IOBCN_GL_W*NTIME,MPI_F,0,MPI_FVCOM_GROUP,IERR)
!JQI     IF(PAR)CALL MPI_BCAST(DIR1_TMP,IOBCN_GL_W*NTIME,MPI_F,0,MPI_FVCOM_GROUP,IERR)
!JQI     IF(PAR)CALL MPI_BCAST(TPEAK_TMP,IOBCN_GL_W*NTIME,MPI_F,0,MPI_FVCOM_GROUP,IERR)
!JQI     IF(PAR)CALL MPI_BCAST(TIMES,NTIME,MPI_F,0,MPI_FVCOM_GROUP,IERR)
!JQI#    endif

!JQI     ELO_TM%NTIMES = NTIME
!JQI     OBC_NTIME = NTIME


!JQI     IF(IOBCN_W > 0)THEN
!JQI       ALLOCATE(OBC_HS(IOBCN_W,NTIME))
!JQI       ALLOCATE(OBC_DIR(IOBCN_W,NTIME))
!JQI       ALLOCATE(OBC_TPEAK(IOBCN_W,NTIME))
!JQI       ALLOCATE(ELO_TM%TIMES(NTIME))
!JQI     ELSE
!JQI       ALLOCATE(OBC_HS(1,NTIME))
!JQI       ALLOCATE(OBC_DIR(1,NTIME))
!JQI       ALLOCATE(OBC_TPEAK(1,NTIME))
!JQI       ALLOCATE(ELO_TM%TIMES(NTIME))
!JQI     END IF

!JQI     OBC_HS    = 0.0
!JQI     OBC_DIR   = 0.0
!JQI     OBC_TPEAK = 0.0

!JQI     ELO_TM%TIMES(1:NTIME) = TIMES

!
!---Map to Local Domain--------------------------------------------------------!
!

!JQI     IF(SERIAL)THEN
!JQI        OBC_HS    = HSC1_TMP
!JQI        OBC_DIR   = DIR1_TMP
!JQI        OBC_TPEAK = TPEAK_TMP
!JQI     END IF
!JQI     IF(IOBCN_W > 0)THEN
!JQI#    if defined (MULTIPROCESSOR)
!JQI     NCNT = 0
!JQI     IF(PAR)THEN
!JQI       DO  J=1,IOBCN_GL_W
!JQI         IF(NLID(I_OBC_GL_W(J)) /= 0 )THEN
!JQI           NCNT = NCNT +1
!JQI!           ELSBC(NCNT,:) = RTEMP1(J,:)
!JQI           OBC_HS(NCNT,:)    = HSC1_TMP(J,:)
!JQI           OBC_DIR(NCNT,:)   = DIR1_TMP(J,:)
!JQI           OBC_TPEAK(NCNT,:) = TPEAK_TMP(J,:)
!JQI         END IF
!JQI       END DO
!JQI     END IF
!JQI#   endif
!JQI     END IF

!JQI     DEALLOCATE(HSC1_TMP,DIR1_TMP,TPEAK_TMP)


   ELSE
     CALL PSTOP
   ENDIF

   RETURN
#  endif   
   END SUBROUTINE SWBOUN
 
!*********************************************************************
!                                                                    *
   SUBROUTINE BCFILE (FBCNAM, BCTYPE)                     

!  (This subroutine has not been used and tested yet)
!                                                                    *
!*********************************************************************
!
!     Reads file data for boundary condition
!
!*********************************************************************
!
   USE OCPCOMM1                                                        
   USE OCPCOMM2                                                        
   USE OCPCOMM4                                                        
   USE SWCOMM2                                                         
   USE SWCOMM3                                                         
   USE SWCOMM4                                                         
   USE M_BNDSPEC                                                       
!
   IMPLICIT NONE

   CHARACTER FBCNAM *(*), BCTYPE *(*)

   INTEGER :: ISTATF, NDSL, NDSD, IOSTAT, IERR, NBOUNC, NANG, NFRE
   INTEGER :: IBOUNC, DORDER
   INTEGER :: IENT,IOPTT
   INTEGER :: NHEDF, NHEDT, NHEDS, IFRE , IANG
   INTEGER :: NQUANT, IQUANT, IBC, II, NBGRPT_PREV,IIPT2
   REAL    :: XP, YP, XP2, YP2
   REAL    :: FREQHZ, DIRDEG, DIRRD1,DIRRAD, EXCV
   CHARACTER BTYPE *4, HEDLIN *80
   LOGICAL         CCOORD                                            

!
   NDSL = 0
   IIPT2 = 0                                                           
!  open data file
   NDSD = 0
   IOSTAT = 0
   CALL FOR (NDSD, FILENM, 'OF', IOSTAT)
!
!     --- initialize array BFILED of BSPFIL                               
!      BSPFIL%BFILED = 0                                                   
!
!     start reading from the data file
      READ (NDSD, '(A)') HEDLIN
!      IF (EQCSTR(HEDLIN,'TPAR')) THEN                                     
!        BTYPE  = 'TPAR'
!        ISTATF = 1
!        IOPTT  = 1
!        NBOUNC = 1
!        NANG   = 0
!        NFRE   = 0
!        NHEDF  = 0
!        NHEDT  = 0
!        NHEDS  = 0
!        DORDER = 0
!        ALLOCATE(BSPFIL%BSPFRQ(NFRE))                                     
!        ALLOCATE(BSPFIL%BSPDIR(NANG))                                     
!        IF (NSTATM.EQ.0) CALL MSGERR (3,'time information not allowed in stationary mode')
!        NSTATM = 1
!      ELSE IF (EQCSTR(HEDLIN,'SWAN')) THEN                                
!      ELSE
!        CALL MSGERR (3, 'unsupported boundary data file')
!      ENDIF
!
!      ALLOCATE(BSPFIL%BSPLOC(NBOUNC))                                     
!      DO IBC = 1, NBOUNC
!         BSPFIL%BSPLOC(IBC) = NBSPEC + IBC                                
!      ENDDO
      NBSPEC = NBSPEC + NBOUNC
!
!     store file reading parameters in array BFILED
!
!      IF (ITEST.GE.80) WRITE(PRINTF,81) NBFILS, NBSPEC,(BSPFIL%BFILED(II), II=1,16)                                  
  81  FORMAT (' array BFILED: ', 2I4, 2(/,8I10))
!
      RETURN
      END SUBROUTINE BCFILE
