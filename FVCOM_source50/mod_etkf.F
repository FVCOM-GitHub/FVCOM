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

MODULE MOD_ETKF 
#  if defined (ETKF_ASSIM)
   USE CONTROL
   IMPLICIT NONE
   SAVE

   CHARACTER(LEN=80)  ETKF_INIT   

   INTEGER      ETKF_RUNS
   INTEGER      ETKF_NENS
   INTEGER      DELTA_ASS
   INTEGER      ETKF_INT
   INTEGER      ETKF_NOBSMAX
   INTEGER      ETKF_START
   INTEGER      ETKF_END
   REAL(SP) ::  ETKF_CINF
   REAL(SP) ::  OBSERR_EL
   REAL(SP) ::  OBSERR_UV 
   REAL(SP) ::  OBSERR_T
   REAL(SP) ::  OBSERR_S

   LOGICAL  ::  ETKF_CTR
   
   LOGICAL  ::  EL_ASSIM
   LOGICAL  ::  UV_ASSIM
   LOGICAL  ::  T_ASSIM 
   LOGICAL  ::  S_ASSIM
   
   LOGICAL  ::  EL_OBS
   LOGICAL  ::  UV_OBS
   LOGICAL  ::  T_OBS
   LOGICAL  ::  S_OBS
   
   REAL(DP)  :: INFLREF = 0.0_DP   ! 10.0d9

! INTEGER SCALARS
   INTEGER   ::    I_INITIAL
   INTEGER   ::   IENS
   INTEGER   ::   ICYC
   INTEGER   ::   NCYC
   INTEGER   ::   N1CYC
   INTEGER   ::   I_REFRUN
   INTEGER   ::   NLOC
   
   CHARACTER(LEN=4)   FCYC

   INTEGER  ::  INOOB             !! 72 RESERVE FOR I/O INPUT OF OBSERVATION FILE
   INTEGER  ::  INOKF             !! 73 FILE I/O PIPE NUMBER
   INTEGER  ::  IOBCKF             
  
   INTEGER  ::  IDUM
! INTEGER ARRAYS

   REAL(DP),PARAMETER :: ZEROD = 0.0_DP

! REAL ARRAYS
   REAL(DP),ALLOCATABLE  :: OBSDATA(:)
   REAL(DP),ALLOCATABLE  :: OBSDATA1(:)
   REAL(DP),ALLOCATABLE  :: OBSDATA2(:)
   REAL(DP),ALLOCATABLE  :: MODDATA(:)   
   REAL(DP),ALLOCATABLE  :: HBHT(:,:)
   REAL(DP),ALLOCATABLE  :: BINV(:,:)
   REAL(DP),ALLOCATABLE  :: STTEMP1(:)
   REAL(DP),ALLOCATABLE  :: STTEMP2(:) 
   REAL(DP),ALLOCATABLE  :: STTEMP3(:) 
!   REAL(DP),ALLOCATABLE  :: STTEMP4(:) 
!   REAL(DP),ALLOCATABLE  :: STTEMP5(:)       
   REAL(DP),ALLOCATABLE  :: STINIT(:)
   REAL(DP),ALLOCATABLE  :: STFCT(:) 

! LOGICAL SWITCH
   LOGICAL  L_REFRUN

  
! VARIABLES FOR ETKF
   REAL(DP),ALLOCATABLE  :: SF1(:,:)
   REAL(DP),ALLOCATABLE  :: SF(:,:)   
   REAL(DP),ALLOCATABLE  :: SA(:,:)
   REAL(DP),ALLOCATABLE  :: SFSF(:,:)
   REAL(DP),ALLOCATABLE  :: SFD(:)
   REAL(DP),ALLOCATABLE  :: SFU(:,:)
   REAL(DP),ALLOCATABLE  :: SFH(:,:)
   REAL(DP),ALLOCATABLE  :: SAH(:,:)
   REAL(DP),ALLOCATABLE  :: SFH1(:)
   REAL(DP),ALLOCATABLE  :: PFHT(:,:)
   REAL(DP),ALLOCATABLE  :: ERRVEC(:)
   REAL(DP),ALLOCATABLE  :: STMEAN(:)

   REAL(DP),ALLOCATABLE  :: STTR(:)     
   REAL(DP),ALLOCATABLE  :: WKTMP(:)     
   INTEGER,ALLOCATABLE   :: STLOC(:)
   REAL(DP),ALLOCATABLE  :: WK(:,:)
   REAL(DP),ALLOCATABLE  :: WK2(:,:)

   INTEGER   TIMEN
   REAL(DP),ALLOCATABLE  :: EL_SRS(:,:),SRS_TMP(:)
   REAL(DP),ALLOCATABLE  :: TIME_SER(:)
   REAL(SP),ALLOCATABLE  :: EL_INV(:)
   REAL(SP)  BC_AMP_ERR(6)
   REAL(SP)  BC_PHA_ERR(6)
   REAL(SP)  PHAI_IJ,FORCE
!-----------------------------------------------------------------------------|  

   CONTAINS !-----------------------------------------------------------------|
            !SET_ETKF_PARAM                                                   | 
            !                                                                 | 
            !                                                                 |

   SUBROUTINE SET_ETKF_PARAM
   
   USE CONTROL
   USE MOD_INP
   IMPLICIT NONE
   INTEGER  :: I, ISCAN, KTEMP
   CHARACTER(LEN=120) :: FNAME
   REAL(SP) REALVEC(150)      
   
!  initialize iens
   IENS = 0  

   FNAME = TRIM(INPDIR)//"/"//trim(casename)//"_assim_etkf.dat"
   
!----------------------------------------------------------------------------|
!     "ETKF_INIT"   !! 
!----------------------------------------------------------------------------|     
   ISCAN = SCAN_FILE(TRIM(FNAME),"ETKF_INIT",CVAL = ETKF_INIT)
    IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING ETKF_INIT: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE ETKF_INIT NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
     END IF
     CALL PSTOP 
   END IF

!----------------------------------------------------------------------------|
!     "ETKF_NENS"   !! 
!----------------------------------------------------------------------------|  
   ISCAN = SCAN_FILE(FNAME,"ETKF_NENS",ISCAL = ETKF_NENS)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING ETKF_NENS: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE ETKF_NENS NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
     END IF
     CALL PSTOP 
   END IF   
   
!----------------------------------------------------------------------------|
!     "DELTA_ASS "  !! 
!----------------------------------------------------------------------------|  
   ISCAN = SCAN_FILE(FNAME,"DELTA_ASS",ISCAL = DELTA_ASS)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING DELTA_ASS: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE DELTA_ASS NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
     END IF
     CALL PSTOP 
   END IF   
   
!----------------------------------------------------------------------------|
!     "ETKF_INT "   !! 
!----------------------------------------------------------------------------|  
   ISCAN = SCAN_FILE(FNAME,"ETKF_INT",ISCAL = ETKF_INT)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING ETKF_INT: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE ETKF_INT NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
     END IF
     CALL PSTOP 
   END IF    
   
!----------------------------------------------------------------------------|
!     "ETKF_NOBSMAX "   !! 
!----------------------------------------------------------------------------|  
   ISCAN = SCAN_FILE(FNAME,"ETKF_NOBSMAX",ISCAL = ETKF_NOBSMAX)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING ETKF_NOBSMAX: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE ETKF_NOBSMAX NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
     END IF
     CALL PSTOP 
   END IF       
   
!----------------------------------------------------------------------------|
!     "ETKF_START "   !! 
!----------------------------------------------------------------------------|  
   ISCAN = SCAN_FILE(FNAME,"ETKF_START",ISCAL = ETKF_START)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING ETKF_START: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE ETKF_START NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
     END IF
     CALL PSTOP 
   END IF

!----------------------------------------------------------------------------|
!     "ETKF_END "   !! 
!----------------------------------------------------------------------------|  
   ISCAN = SCAN_FILE(FNAME,"ETKF_END",ISCAL = ETKF_END)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING ETKF_END: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE ETKF_END NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
     END IF
     CALL PSTOP 
   END IF

!------------------------------------------------------------------------------|
!     "BC_AMP_ERR"   !! 
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"BC_AMP_ERR",FVEC = REALVEC,NSZE = KTEMP)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING BC_AMP_ERR: ',ISCAN
     CALL PSTOP
   END IF
   IF(MSR)THEN
     IF(KTEMP /= 6)THEN
       WRITE(*,*)'NUMBER OF SPECIFIED TIDAL COMSTITUENTS ERROR ARE NOT EQUAL TO 6' 
     END IF
   END IF
  
   BC_AMP_ERR(1:6)= REALVEC(1:6)
!------------------------------------------------------------------------------|
!     "BC_PHA_ERR"   !! 
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"BC_PHA_ERR",FVEC = REALVEC,NSZE = KTEMP)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING BC_PHA_ERR: ',ISCAN
     CALL PSTOP
   END IF
   IF(MSR)THEN
     IF(KTEMP /= 6)THEN
       WRITE(*,*)'NUMBER OF SPECIFIED TIDAL COMSTITUENTS ERROR ARE NOT EQUAL TO 6' 
     END IF
   END IF
  
   BC_PHA_ERR(1:6)= REALVEC(1:6)
!----------------------------------------------------------------------------|
!   "ETKF_CINF" MAX LONG DISTANCE OF CORRELATIN 
!----------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"ETKF_CINF",FSCAL = ETKF_CINF)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING ETKF_CINF: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE ETKF_CINF NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
     END IF
     CALL PSTOP  
   END IF

!----------------------------------------------------------------------------|
!   "OBSERR_EL" 
!----------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"OBSERR_EL",FSCAL = OBSERR_EL)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING OBSERR_EL: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE OBSERR_EL NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
     END IF
     CALL PSTOP  
   END IF

!----------------------------------------------------------------------------|
!   "OBSERR_UV" 
!----------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"OBSERR_UV",FSCAL = OBSERR_UV)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING OBSERR_UV: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE OBSERR_UV NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
     END IF
     CALL PSTOP  
   END IF

!----------------------------------------------------------------------------|
!   "OBSERR_T" 
!----------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"OBSERR_T",FSCAL = OBSERR_T)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING OBSERR_T: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE OBSERR_T NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
     END IF
     CALL PSTOP  
   END IF

!----------------------------------------------------------------------------|
!   "OBSERR_S" 
!----------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"OBSERR_S",FSCAL = OBSERR_S)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING OBSERR_S: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE OBSERR_S NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
     END IF
     CALL PSTOP  
   END IF

!----------------------------------------------------------------------------|
!   EL_ASSIM  -  OPTION FOR CHOSING ELEVATION AS ASSIMILATION VARIABLES       
!----------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"EL_ASSIM",LVAL = EL_ASSIM) 
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING EL_ASSIM: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE EL_ASSIM NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
     END IF
     CALL PSTOP
   ENDIF   
   
!----------------------------------------------------------------------------|
!   EL_OBS  -  ELEVATION OBSERVATION DATA OPTION                              
!----------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"EL_OBS",LVAL = EL_OBS) 
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING EL_OBS: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE EL_OBS NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
     END IF
     CALL PSTOP
   ENDIF   

!----------------------------------------------------------------------------|
!   UV_ASSIM  -  OPTION FOR CHOSING CURRENT AS ASSIMILATION VARIABLES         
!----------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"UV_ASSIM",LVAL = UV_ASSIM) 
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING UV_ASSIM: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE UV_ASSIM NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
     END IF
     CALL PSTOP
   ENDIF 

!----------------------------------------------------------------------------|
!   UV_OBS  -  CURRENT OBSERVATION DATA OPTION                                
!----------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"UV_OBS",LVAL = UV_OBS) 
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING UV_OBS: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE UV_OBS NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
     END IF
     CALL PSTOP
   ENDIF 

!----------------------------------------------------------------------------|
!   T_ASSIM  -   OPTION FOR CHOSING TEMPERATURE AS ASSIMILATION VARIABLES     
!----------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"T_ASSIM",LVAL = T_ASSIM) 
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING T_ASSIM: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE T_ASSIM NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
     END IF
     CALL PSTOP
   ENDIF 

!----------------------------------------------------------------------------|
!   T_OBS  -  TEMPERATURE OBERVATION DATA OPTION                              
!----------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"T_OBS",LVAL = T_OBS) 
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING T_OBS: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE T_OBS NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
     END IF
     CALL PSTOP
   ENDIF 

!----------------------------------------------------------------------------|
!   S_ASSIM  -  OPTION FOR CHOSING SALINITY AS ASSIMILATION VARIABLES         
!----------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"S_ASSIM",LVAL = S_ASSIM) 
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING S_ASSIM: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE S_ASSIM NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
     END IF
     CALL PSTOP
   ENDIF

!----------------------------------------------------------------------------|
!   S_OBS  -  SALINITY OBSERVATION DATA OPTION                                
!----------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"S_OBS",LVAL = S_OBS) 
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING S_OBS: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE S_OBS NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
     END IF
     CALL PSTOP
   ENDIF

!----------------------------------------------------------------------------|
!   ETKF_CTR  -  CONTROL RUN OR USE ENSEMBLE AS THE TRUE STATE         
!----------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"ETKF_CTR",LVAL = ETKF_CTR) 
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING ETKF_CTR: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE ETKF_CTR NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
     END IF
     CALL PSTOP
   ENDIF
   
!==============================================================================|
!            SCREEN REPORT OF SET ETKF VARIABlES                        !
!==============================================================================|
   IF(MSR) THEN  
     WRITE(IPT,*) '!                                                    !'     
     WRITE(IPT,*) '!------SPECIFY ETKF DATA ASSIMINATION PARAMETERS-----!'     
     WRITE(IPT,*) '!                                                    !'     
     WRITE(IPT,*) '!  #CONTROL RUN OR USE ENSEMBLE AS THE TRUE STATE    :',ETKF_CTR
     WRITE(IPT,*) '!  # SPECIFICY INITIAL PERTUBATION FIELD             :',ETKF_INIT
     WRITE(IPT,*) '!  # GLOBAL NUMBER OF ENSEMBLES                      :',ETKF_NENS
     WRITE(IPT,*) '!  # ASSIMILATON TIME INTERVAL IN SECONDS            :',DELTA_ASS
     WRITE(IPT,*) '!  # ASSIMILATION TIME INTERVAL/FILE OUTPUT INTERVAL :',ETKF_INT
     WRITE(IPT,*) '!  # MAXIMUM NUMBER OF THE OBSERVATION STATIONS      :',ETKF_NOBSMAX
     WRITE(IPT,*) '!  # ASSIMILATION START TIME                         :',ETKF_START
     WRITE(IPT,*) '!  # ASSIMILATION END TIME                           :',ETKF_END
     WRITE(IPT,*) '!  # TIDAL AMPLITUDE ERROR RANGE SPECIFIED           :',(BC_AMP_ERR(I),I=1,6)
     WRITE(IPT,*) '!  # TIDAL PHASE ERROR RANGE SPECIFIED               :',(BC_PHA_ERR(I),I=1,6) 
     WRITE(IPT,*) '!  # MAX DISTANCE OF CORRELATIN                      :',ETKF_CINF
     WRITE(IPT,*) '!  # ELEVATION AS ASSIMILATION VARIABLES             :',EL_ASSIM
     WRITE(IPT,*) '!  # ELEVATION AS OBSERVATION DATA OPTION            :',EL_OBS
     WRITE(IPT,*) '!  # CURRENTS AS ASSIMILATION VARIABLES              :',UV_ASSIM
     WRITE(IPT,*) '!  # CURRENTS OBSERVATION DATA OPTION                :',UV_OBS
     WRITE(IPT,*) '!  # TEMPERATURE AS ASSIMILATION VARIABLES           :',T_ASSIM
     WRITE(IPT,*) '!  # TEMPERATURE OBSERVATION DATA OPTION             :',T_OBS
     WRITE(IPT,*) '!  # SALINITY AS ASSIMILATION VARIABLES              :',S_ASSIM
     WRITE(IPT,*) '!  # SALINITY OBSERVATION DATA OPTION                :',S_OBS
     IF(EL_OBS) WRITE(IPT,*) '!  # EL OBSERVATION ERROR SPECIFIED                  :',OBSERR_EL
     IF(UV_OBS) WRITE(IPT,*) '!  # UV OBSERVATION ERROR SPECIFIED                  :',OBSERR_UV
     IF(T_OBS)  WRITE(IPT,*) '!  # T OBSERVATION ERROR SPECIFIED                   :',OBSERR_T
     IF(S_OBS)  WRITE(IPT,*) '!  # S OBSERVATION ERROR SPECIFIED                   :',OBSERR_S
     WRITE(IPT,*) '!                                                   !'           
  ENDIF
   
   RETURN
   END SUBROUTINE SET_ETKF_PARAM

   SUBROUTINE ALLOC_VARS_ETKF
   
   USE LIMS
   USE CONTROL
   IMPLICIT NONE  
    
   INTEGER         STDIM               !!dimension of state vector,number of elements for all assimilation
   INTEGER         NDB9 
   REAL(DP)  ::    MEMCNT9 
   
   STDIM = 0
   IF(EL_ASSIM) STDIM = STDIM + MGL
   IF(UV_ASSIM) STDIM = STDIM + 2*NGL*KBM1
   IF(T_ASSIM)  STDIM = STDIM + MGL*KBM1
   IF(S_ASSIM)  STDIM = STDIM + MGL*KBM1

   MEMCNT9 = 0.0_DP
   NDB9    = 2
   
   ALLOCATE(OBSDATA(ETKF_NOBSMAX))    ;OBSDATA   = ZEROD
   ALLOCATE(OBSDATA1(ETKF_NOBSMAX))   ;OBSDATA1  = ZEROD
   ALLOCATE(OBSDATA2(ETKF_NOBSMAX))   ;OBSDATA2  = ZEROD
   ALLOCATE(MODDATA(ETKF_NOBSMAX))    ;MODDATA   = ZEROD
   
   ALLOCATE(HBHT(ETKF_NOBSMAX,ETKF_NOBSMAX))    ;HBHT   = ZEROD
   ALLOCATE(BINV(ETKF_NOBSMAX,ETKF_NOBSMAX))    ;BINV   = ZEROD

   ALLOCATE(STTEMP1(STDIM))           ;STTEMP1   = ZEROD
   ALLOCATE(STTEMP2(STDIM))           ;STTEMP2   = ZEROD   
   ALLOCATE(STTEMP3(STDIM))           ;STTEMP3   = ZEROD   
!   ALLOCATE(STTEMP4(STDIM))           ;STTEMP4   = ZEROD
!   ALLOCATE(STTEMP5(STDIM))           ;STTEMP5   = ZEROD
   ALLOCATE(STINIT(STDIM))            ;STINIT    = ZEROD
   ALLOCATE(STFCT(STDIM))             ;STFCT     = ZEROD
     
   ALLOCATE(WK(ETKF_NOBSMAX,ETKF_NOBSMAX))      ;WK     = ZEROD
   ALLOCATE(WK2(ETKF_NOBSMAX,ETKF_NOBSMAX))     ;WK2    = ZEROD
         
   ALLOCATE(SF1(STDIM,ETKF_NENS))     ;SF1       = ZEROD
   ALLOCATE(SF(STDIM,ETKF_NENS))      ;SF        = ZEROD  
   ALLOCATE(SA(STDIM,ETKF_NENS+ETKF_NOBSMAX))   ;SA        = ZEROD
   ALLOCATE(SFSF(ETKF_NENS,ETKF_NENS)) ;SFSF     = ZEROD
   ALLOCATE(SFU(ETKF_NENS,ETKF_NENS)) ;SFU       = ZEROD
   ALLOCATE(SFD(ETKF_NENS))           ;SFD       = ZEROD
   ALLOCATE(SFH(ETKF_NENS,ETKF_NOBSMAX))        ;SFH       = ZEROD
   ALLOCATE(SFH1(ETKF_NENS))          ;SFH1      = ZEROD
   ALLOCATE(PFHT(STDIM,ETKF_NOBSMAX)) ;PFHT      = ZEROD
   
   ALLOCATE(ERRVEC(STDIM))            ;ERRVEC    = ZEROD
   ALLOCATE(STMEAN(STDIM))            ;STMEAN    = ZEROD
   RETURN
   END SUBROUTINE ALLOC_VARS_ETKF
          
   SUBROUTINE DEALLOC_VARS_ETKF
   USE LIMS

   DEALLOCATE(OBSDATA,OBSDATA1,OBSDATA2,MODDATA,STTEMP1,STTEMP2,STTEMP3)
   DEALLOCATE(STINIT,STFCT,WK,WK2,SF1,SF,SA,SFSF,SFU,SFD,SFH, SFH1,PFHT)
   RETURN
   END SUBROUTINE DEALLOC_VARS_ETKF

   SUBROUTINE SET_ASSIM_ETKF

   USE LIMS
   USE CONTROL
   USE ALL_VARS
   USE BCS
#  if defined (MULTIPROCESSOR)
   USE MOD_PAR
#  endif
   IMPLICIT NONE
   
   INTEGER I,J,IERR
   CHARACTER(LEN=100) MKANLDIR, MKFCTDIR, MKERRDIR, MKOUTDIR
   
   IF(MSR) THEN 
#  if !defined (DOS)
      MKANLDIR = "mkdir -p "//TRIM(OUTDIR)//"/anl"
      MKFCTDIR = "mkdir -p "//TRIM(OUTDIR)//"/fct"
      MKERRDIR = "mkdir -p "//TRIM(OUTDIR)//"/out_err"  
      MKOUTDIR = "mkdir -p "//TRIM(OUTDIR)//"/flow"
#     if !defined (CRAY)
         CALL SYSTEM( TRIM(MKANLDIR) )  
         CALL SYSTEM( TRIM(MKFCTDIR) )
         CALL SYSTEM( TRIM(MKERRDIR) )
         CALL SYSTEM( TRIM(MKOUTDIR) )
#     endif
#     if defined (CRAY)
         CALL CRAY_SYSTEM_CALL(TRIM(MKANLDIR))
         CALL CRAY_SYSTEM_CALL(TRIM(MKFCTDIR))
         CALL CRAY_SYSTEM_CALL(TRIM(MKERRDIR))    
         CALL CRAY_SYSTEM_CALL(TRIM(MKOUTDIR))
#     endif             
#  endif         
   ENDIF

   ICYC = 0
   I_INITIAL = IINT            !! BE CAREFULE ABOUT IINT, IT SHOULD BE RELATED TO THE INITIAL FIELD YOU BUILD
   NCYC  = IEND/DELTA_ASS      !! GET THE ASSIMILATON LOOP NUMBERS           
   N1CYC = (I_INITIAL+DELTA_ASS)/DELTA_ASS   

   IF(MSR) THEN
     OPEN(INOOB,FILE=TRIM(INPDIR)//"/"//trim(casename)//"_assim_etkf.dat",FORM='FORMATTED')

     OPEN(74,FILE=TRIM(OUTDIR)//'/out_err/EnSp.dat')
     WRITE(74,*) ' Icyc       aa        infl1        infl2      inflold    inflation' 
   
     OPEN(75,FILE=TRIM(OUTDIR)//'/out_err/ErrOut.dat')
     WRITE(75,*) ' Icyc  fctobserr fctrmserr anlobserr anlrmserr avgrmserr anlrmserr/avgrmserr'
   ENDIF
   
   IF(MSR) THEN
      IF(ETKF_INIT == 'default') CALL SET_INI_ETKF   ! only used for test case!
      CALL PERT_BC 
   ENDIF
   
#  if defined (MULTIPROCESSOR)
      IF(PAR)CALL MPI_BARRIER(MPI_FVCOM_GROUP,IERR)
#  endif 

   RETURN
   END SUBROUTINE SET_ASSIM_ETKF


   SUBROUTINE ETKF_ASS
   
   USE LIMS
   USE CONTROL
   USE ALL_VARS
   IMPLICIT NONE

   INTEGER  I, J, K, II, JJ, KK
   INTEGER  STDIM
   REAL(DP) DISTST, DELT
   CHARACTER(LEN=80) TEXT    
   
   REAL(DP)  RNOBS, GASDEV
   REAL(DP)  RSCALE
   REAL(DP)  ETKF_CINF2
   CHARACTER(LEN=120) FNAM, GNAM
   CHARACTER(LEN=4)  FENS
   CHARACTER(LEN=4)  JEOF
         
   REAL(DP)  ERR2_INN_FCT, ERR2_INN_ANL
   REAL(SP)  ERR1, ERR2
   REAL(DP)  SUM0, SUM9, AVGRMSERR,FCTRMSERR,ANLRMSERR,FCTOBSERR,ANLOBSERR
   REAL(DP)  INFLATION,INFL1,INFL2,INFL1_2,INFLOLD    

   INTEGER   LWORK4, LDVT, LWORK5, RCODE
   REAL(DP)  VT   

   REAL(DP),ALLOCATABLE    ::   WORK4(:)
   REAL(DP),ALLOCATABLE    ::   WORK5(:)
  
   STDIM = 0
   IF(EL_ASSIM) STDIM = STDIM + MGL
   IF(UV_ASSIM) STDIM = STDIM + 2*NGL*KBM1
   IF(T_ASSIM)  STDIM = STDIM + MGL*KBM1
   IF(S_ASSIM)  STDIM = STDIM + MGL*KBM1  
      
   LDVT   = 1   
   LWORK4 = 5*ETKF_NENS   
   LWORK5 = 5*(ETKF_NENS+ETKF_NENS)

   ALLOCATE(WORK4(LWORK4))  ; WORK4 = ZEROD
   ALLOCATE(WORK5(LWORK5))  ; WORK5 = ZEROD
   ALLOCATE(WKTMP(STDIM))   ; WKTMP = ZEROD
   ALLOCATE(STTR(STDIM))    ; STTR  = ZEROD
   ALLOCATE(STLOC(ETKF_NOBSMAX))   ; STLOC = 0   !!STLOC:    
   
   CALL GETOBSLOC
   CALL ALLOC_VARS_ETKF   
!CC----------------------------------------------------------CC
!CC Get the Observation Covariance Matrix: Wk Nobs * Nobs    CC
!CC----------------------------------------------------------CC

   WK = 0.0_DP
   DO I=1, NLOC
      WK(I,I) = WKTMP(STLOC(I))
      WK(I,I) = WK(I,I)**2.0_DP      
   ENDDO   
   
   DO K=1,ETKF_NENS

      WRITE(FENS,'(I4.4)') K
             
      GNAM= TRIM(OUTDIR)//'/fct/restart'//FENS//'.dat'

      OPEN(INOKF,FILE=TRIM(GNAM),FORM='UNFORMATTED',STATUS='OLD')
      CALL GR2ST(INOKF)        ! RETURN TO SF

      DO I=1,STDIM
        SF(I,K)=STTEMP1(I)
      ENDDO

   ENDDO   

   IF(ETKF_CTR) THEN

      WRITE(FENS,'(I4.4)') ETKF_RUNS
             
      GNAM= TRIM(OUTDIR)//'/fct/restart'//FENS//'.dat'

      OPEN(INOKF,FILE=TRIM(GNAM),FORM='UNFORMATTED',STATUS='OLD')
      CALL GR2ST(INOKF)        ! RETURN TO STFCT

      DO I=1,STDIM
        STFCT(I)=STTEMP1(I)
      ENDDO

      DO I=1, STDIM
        SUM0 = 0.0D0
        DO K=1, ETKF_NENS
          SUM0 = SUM0 + SF(I,K)
        ENDDO
        STMEAN(I) = SUM0/DBLE(ETKF_NENS)
        SUM0 = 0.0D0
        DO K=1, ETKF_NENS
          SF(I,K) = (SF(I,K)-STFCT(I))/DSQRT(DBLE(ETKF_NENS))
          SUM0 = SUM0 + SF(I,K)**2.0
          SF1(I,K) = SF(I,K)
        ENDDO
        ERR2 = ERR2 + SUM0
        STTEMP2(I) = SUM0
      ENDDO
   ELSE
      DO I=1, STDIM
        SUM0 = 0.0D0
        DO K=1, ETKF_NENS
          SUM0 = SUM0 + SF(I,K)
        ENDDO
        STMEAN(I) = SUM0/DBLE(ETKF_NENS)
        SUM0 = 0.0D0
        DO K=1, ETKF_NENS
          SF(I,K) = (SF(I,K)-STMEAN(I))/DSQRT(DBLE(ETKF_NENS-1))
          SUM0 = SUM0 + SF(I,K)**2.0
          SF1(I,K) = SF(I,K)
        ENDDO
        ERR2 = ERR2 + SUM0
        STTEMP2(I) = SUM0
      ENDDO
   ENDIF
   ERR2 = DSQRT(ERR2/DBLE(STDIM)) 

! Also size of total error (full state vector)   ! I marked this part cause we don't have whole domain true value in real case 
!   icycref = icyc*itimes + i_refrun + i_trurun
!   call klm_gettruth(icycref,sttemp4)
!   call subst(stmean,sttemp4,errvec)
!   text='Mean(fct)'
!   call print_err(errvec,stdim,text,err2_tot_fct)


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!  Update ensemble perturbations through EnSRF
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!
!  Calculate sum of eigenvalues for inflation factor
!
   DO J = 1, NLOC
      DO I= 1, ETKF_NENS
         SFH(I,J) = SF(STLOC(J),I)/DSQRT(WK(J,J))
      ENDDO
   ENDDO

   CALL DGEMM('n','t',ETKF_NENS,ETKF_NENS,NLOC,1.0D0,SFH,ETKF_NENS,SFH,ETKF_NENS,0.0D0,SFSF,ETKF_NENS)

! SVD of SfSf: SfSf=(psi)*R^-1*(psi)'=SfH*SfH'=SfU*SfD*SfU'
   CALL DGESVD('A','N',ETKF_NENS,ETKF_NENS,SFSF,ETKF_NENS,SFD,SFU,ETKF_NENS,VT,LDVT,WORK4,LWORK4,RCODE)

   OPEN(INOKF,FILE=TRIM(OUTDIR)//'/out_err/Lambda'//fcyc//'.dat')
   INFL1=0.
   DO J=1,ETKF_NENS
      INFL1 = INFL1 + SFD(J)
      WRITE(INOKF,'(I5,E15.7)') ICYC,SFD(J)
   ENDDO
   CLOSE(INOKF)
   
!
!  Serial Processing
!   
   ETKF_CINF2 = DBLE(ETKF_CINF)/2.0_DP
   DO J = 1, NLOC
!calculate  (H Xf)'(H Xf)
      DO I=1,ETKF_NENS  
         SFH1(I) = SF(STLOC(J),I)
      ENDDO
      SUM0 = 0.0D0
      DO I=1,ETKF_NENS
         SUM0 = SUM0 + SFH1(I)**2
      ENDDO
      SUM0 = SUM0 + WK(J,J)
      SUM0 = SUM0 + DSQRT(WK(J,J)*SUM0)
      DO I = 1,ETKF_NENS
         DO K = 1, ETKF_NENS
            SFSF(I,K) = SFH1(I)*SFH1(K)
         ENDDO
      ENDDO
      DO I=1,ETKF_NENS
         DO K = 1, ETKF_NENS
            SFU(I,K)= 1.0_SP/SUM0*SFSF(I,K)
         ENDDO
      ENDDO

      DO I = 1, STDIM
         DELT = DISTST(STLOC(J),I)

        IF(DELT > ETKF_CINF) THEN
           RSCALE = 0.
        ELSE IF(DELT > ETKF_CINF2) THEN
           RSCALE = 1.0_DP/12.0_DP*(DELT/ETKF_CINF2)**5 &
                    -1.0_DP/2.0_DP*(DELT/ETKF_CINF2)**4 &
                    +5.0_DP/8.0_DP*(DELT/ETKF_CINF2)**3 &
                    +5.0_DP/3.0_DP*(DELT/ETKF_CINF2)**2 &
                    -5.0_DP*(DELT/ETKF_CINF2) &
                    +4.0_DP -2.0_DP/3.0_DP*(ETKF_CINF2/DELT)
        ELSE
           RSCALE = -1.0_DP/4.0_DP*(DELT/ETKF_CINF2)**5 &
                    +1.0_DP/2.0_DP*(DELT/ETKF_CINF2)**4 &
                    +5.0_DP/8.0_DP*(DELT/ETKF_CINF2)**3 &
                    -5.0_DP/3.0_DP*(DELT/ETKF_CINF2)**2 &
                    +1
        ENDIF

        DO JJ = 1, ETKF_NENS
           SUM0 = 0.0D0
           DO K = 1, ETKF_NENS
              SUM0 = SUM0 + SF(I,K)*SFU(K,JJ)
           ENDDO
           SA(I,JJ) = SF(I,JJ) - SUM0*RSCALE
        ENDDO
      ENDDO

      DO I = 1, STDIM
        DO JJ = 1, ETKF_NENS  
           SF(I,JJ) = SA(I,JJ)
        ENDDO
      ENDDO
   
   ENDDO
   
!ccccccccccccccccccccccccccccccccccccccccccc
!
! Ensenble mean or Control Analysis 
!
!ccccccccccccccccccccccccccccccccccccccccccc

! Get the observations and model counterparts

! y=H(x_obs) -> obsdata
!         call klm_getobsdata(stloc,Nloc,icycref,obsdata)
   DO I=1,NLOC
      OBSDATA(I) = STTR(STLOC(I))
   ENDDO

! H(x_fct)   -> moddata
   IF (ETKF_CTR) THEN
!         call klm_getmoddata(stloc,Nloc,stfct,moddata)
!         call copyst(stfct,stinit) 
     DO I=1,NLOC
       MODDATA(I) = STFCT(STLOC(I))
     ENDDO
     STINIT = STFCT
   ELSE
!         call klm_getmoddata(stloc,Nloc,stmean,moddata)
!         call copyst(stmean,stinit) 
     DO I=1,NLOC
       MODDATA(I) = STMEAN(STLOC(I))
     ENDDO
     STINIT = STMEAN
   ENDIF

! Perturb the observation with random errors of sqrt(R(i,i))
   IDUM = -31
   DO I=1, NLOC
     RNOBS = GASDEV(IDUM)
     OBSDATA(I) = OBSDATA(I) + DSQRT(WK(I,I))*RNOBS
   ENDDO

!
! Calculate innovation vector y' = H(x_fct)-H(x_obs) ->moddata
!
   INFL2 = 0.
   DO K = 1, NLOC
     MODDATA(K) = OBSDATA(K) - MODDATA(K)
     INFL2 = INFL2 + MODDATA(K)**2/WK(K,K)
   ENDDO
!   WRITE(*,*)'obs,mod:',obsdata(1),moddata(1),idum
   TEXT ='FCT(INNV)'
   CALL PRINT_ERR(MODDATA,NLOC,TEXT,ERR2_INN_FCT)

!
! Compute correction by applying gain matrix to innovation vector
! Use the previous forecast perturbations before serial processing
! Also apply localization in the control analysis update
!

! calculate Xf'H'
   DO I = 1, ETKF_NENS
     DO J = 1,NLOC
       SAH(I,J) = SF1(STLOC(J),I)
     ENDDO
   ENDDO

! calculate Xf Xf'H'
   DO I = 1, STDIM
     DO J = 1,NLOC
       SUM0 = 0.0D0
       DO K = 1, ETKF_NENS
          SUM0 = SUM0 + SF1(I,K)*SAH(K,J)
       ENDDO
       PFHT(I,J) = SUM0
     ENDDO
   ENDDO   

! Schur product of a correlation fn

    DO I = 1, STDIM
       DO J = 1, NLOC
         DELT = DISTST (STLOC(J),I)

         IF(DELT > ETKF_CINF) THEN
           RSCALE = 0.
         ELSE IF(DELT > ETKF_CINF2) THEN
           RSCALE = 1.0_DP/12.0_DP*(DELT/ETKF_CINF2)**5 &
                    -1.0_DP/2.0_DP*(DELT/ETKF_CINF2)**4 &
                    +5.0_DP/8.0_DP*(DELT/ETKF_CINF2)**3 &
                    +5.0_DP/3.0_DP*(DELT/ETKF_CINF2)**2 &
                    -5.0_DP*(DELT/ETKF_CINF2) &
                    +4.0_DP -2.0_DP/3.0_DP*(ETKF_CINF2/DELT)
         ELSE
           RSCALE = -1.0_DP/4.0_DP*(DELT/ETKF_CINF2)**5 &
                    +1.0_DP/2.0_DP*(DELT/ETKF_CINF2)**4 &
                    +5.0_DP/8.0_DP*(DELT/ETKF_CINF2)**3 &
                    -5.0_DP/3.0_DP*(DELT/ETKF_CINF2)**2 &
                    +1
         ENDIF

         PFHT(I,J) = PFHT(I,J)*RSCALE
       ENDDO
    ENDDO

! calculate H Xf Xf'H'
    DO I = 1 , NLOC
       DO J = 1, NLOC
          HBHT(I,J) = PFHT(STLOC(I),J)
       ENDDO
    ENDDO

! calculate H Xf Xf'H' + R
    DO I = 1, NLOC
       HBHT(I,I) = HBHT(I,I) + WK(I,I)
    ENDDO
    
! Invert (H Xf Xf'H' + R)
    CALL DGESVD('A','N',NLOC,NLOC,HBHT,ETKF_NOBSMAX,OBSDATA2,WK2,ETKF_NOBSMAX,VT,LDVT,WORK5,LWORK5,RCODE)

    DO I = 1,NLOC
       DO J = 1,NLOC
          HBHT(I,J) = WK2(I,J)/OBSDATA2(J)
       ENDDO
    ENDDO
    DO I = 1, NLOC
       DO J=1, NLOC
          SUM0 = 0.0D0
            DO K = 1, NLOC
               SUM0 = SUM0 + HBHT(I,K)*WK2(J,K)
            ENDDO
          BINV(I,J) = SUM0
       ENDDO
    ENDDO    

! calculate (H Xf Xf'H' + R)^-1 (y - H Xf)
    CALL DGEMM('n','n',NLOC,1,NLOC,1.0d0,BINV,ETKF_NOBSMAX,MODDATA,ETKF_NOBSMAX,0.0d0,OBSDATA1,ETKF_NOBSMAX)

! calculate Xf Xf'H' (H Xf Xf'H' + R)^-1 (y - H Xf)
    CALL DGEMM('n','n',STDIM,1,NLOC,1.0d0,PFHT,STDIM,OBSDATA1,ETKF_NOBSMAX,0.0d0,STTEMP1,STDIM)
    
! calculate analysis xa= xf + Xf Xf'H' (H Xf Xf'H' + R)^-1 (y - H Xf)
    DO I = 1, STDIM
       STINIT(I) = STINIT(I) + STTEMP1(I)
    ENDDO

!
! Output resulting analysis                              ! Does not output this netcdf file
!                                                        ! can be modified based on personal use 
!        filename='Contfct.cdf'
!        call plotstate_cdf(filename,stfct)
!        filename='EnsMfct.cdf'
!        call plotstate_cdf(filename,stmean)
!        filename='analysis.cdf'
!        call plotstate_cdf(filename,stinit)
!      if(l_control)then
!        call subst(stinit,stfct,sttemp5)
!      else
!        call subst(stinit,stmean,sttemp5)
!      endif
!        filename='incstate.cdf'
!        call plotstate_cdf(filename,sttemp5)
!        call st2gr(stinit,q2)
!        new=(icyc)*delta_ass+ass_start
!        nl_istep=0
!        call nl_wrmean_cdf(1)

!
! Compute diagnostic on analysis - first error in obs space
!
!        call klm_getmoddata(stloc, Nloc,stinit,moddata)
     DO I=1,NLOC
       MODDATA(I) = STINIT(STLOC(I))
     ENDDO
     DO K = 1, NLOC
       ERRVEC(K) = OBSDATA(K) - MODDATA(K)
     ENDDO
     TEXT='ANL(INNV)'
     CALL PRINT_ERR(ERRVEC,NLOC,TEXT,ERR2_INN_ANL)

! Also size of total error                                  ! does not have whole domain true value in real case
!        call klm_gettruth(icycref,sttemp4)
!        call subst(stinit,sttemp4,errvec)
!        text='ANL(full)'
!        call print_err(errvec,stdim,text,err2_tot_anl)


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!  Update ensemble members 
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

! Calculate an inflation factor=[innovation'*innovation-Nloc]/sum_of_lambda

   INFLATION = (INFL2-NLOC)/INFL1

   IF(INFLATION < INFLREF) THEN
      INFLOLD = INFLATION
      INFLATION = 1.0_DP
   ELSE
      INFLOLD = INFLATION
   ENDIF

   INFL1_2 = DSQRT(DBLE(INFLATION))


! output xa = Xa + xa_control or Xa + xa_mean    

   IF(ETKF_CTR) THEN
     DO K = 1, ETKF_NENS
       WRITE(JEOF,'(I4.4)') K
       DO I=1, STDIM
          STTEMP1(I) = SA(I,K)*DSQRT(DBLE(ETKF_NENS))*INFL1_2 + STINIT(I)
       ENDDO
       FNAM = TRIM(OUTDIR)//'/fct/restart'//JEOF//'.dat'            
       GNAM = TRIM(OUTDIR)//'/anl/restart'//JEOF//'.dat'
       CALL FCT2ANL(FNAM,GNAM)
       DO I = 1, STDIM
         SF(I,K) = STTEMP1(I)
       ENDDO
     ENDDO
   ELSE
     DO K = 1, ETKF_NENS
       WRITE(JEOF,'(I4.4)') K
       DO I=1, STDIM
          STTEMP1(I) = SA(I,K)*DSQRT(DBLE(ETKF_NENS-1.))*INFL1_2 + STINIT(I)
       ENDDO
       FNAM = TRIM(OUTDIR)//'/fct/restart'//JEOF//'.dat'            
       GNAM = TRIM(OUTDIR)//'/anl/restart'//JEOF//'.dat'
       CALL FCT2ANL(FNAM,GNAM)
       DO I=1,STDIM
         SF(I,K) = STTEMP1(I)
       ENDDO
    ENDDO
  ENDIF

!!--------Calculate the average rms error of each ensemble member -------CC

!  TEXT = 'EACHRMSERR'
!  SUM0 = 0.0D0
!  DO K = 1, ETKF_NENS
!    DO I=1, STDIM
!       ERRVEC(I) = SF(I,K) - STTEMP4(I)
!    ENDDO
!    CALL PRINT_ERR(ERRVEC,STDIM,TEXT,AVGRMSERR)
!    SUM0 = SUM0 + AVGRMSERR
!  ENDDO
!  AVGRMSERR = SUM0/DBLE(ETKF_NNES)


!!--------Calculate the Ensemble Mean -------CC

   ERR1 = 0.0D0
   IF(ETKF_CTR) THEN
      DO I=1,STDIM
        SUM0 = 0.0D0
        DO K=1, ETKF_NENS
           SUM0 = SUM0 + SF(I,K)
        ENDDO
        STTEMP1(I) = SUM0/DBLE(ETKF_NENS)
        SUM0 = 0.
        DO K=1,ETKF_NENS
           SF(I,K)=(SF(I,K)-STINIT(I))/DSQRT(DBLE(ETKF_NENS))
           SUM0 = SUM0 + SF(I,K)**2
        ENDDO
        ERR1 = ERR1 + SUM0
        STTEMP3(I) = SUM0
      ENDDO
   ELSE
      DO I=1,STDIM
         SUM0 = 0.0D0
         DO K=1, ETKF_NENS
            SUM0 = SUM0 + SF(I,K)
         ENDDO
         STTEMP1(I) = SUM0/DBLE(ETKF_NENS)
         SUM0 = 0.
         DO K=1, ETKF_NENS
            SF(I,K)=(SF(I,K)-STTEMP1(I))/DSQRT(DBLE(ETKF_NENS-1))
            SUM0 = SUM0 + SF(I,K)**2
         ENDDO
         ERR1 = ERR1 + SUM0
         STTEMP3(I) = SUM0
      ENDDO
   ENDIF
   ERR1 = DSQRT(ERR1/DBLE(STDIM))

   WRITE(74,'(I5,6D16.7)') ICYC, ERR1, ERR2, INFL1, INFL2, INFLOLD, INFLATION
   WRITE(75,'(I4,6E12.4)') ICYC, ERR2_INN_FCT, ERR2_INN_ANL
   
   CALL DEALLOC_VARS_ETKF
   DEALLOCATE(WORK4, WORK5, STTR, STLOC, WKTMP)   
   RETURN   
   END SUBROUTINE ETKF_ASS

   SUBROUTINE GR2ST(INF)
   
   USE LIMS
   USE ALL_VARS
#  if defined (WATER_QUALITY)
   USE MOD_WQM
#  endif
#  if defined (EQUI_TIDE)
   USE MOD_EQUITIDE
#  endif
#  if defined (ATMO_TIDE)
   USE MOD_ATMOTIDE
#  endif
   IMPLICIT NONE
   
   INTEGER I,K,N1,IDUMMY 
   INTEGER INF
   INTEGER ITMP
   
   REAL(SP), ALLOCATABLE :: UTMP(:,:)
   REAL(SP), ALLOCATABLE :: VTMP(:,:)
   REAL(SP), ALLOCATABLE :: S1TMP(:,:)
   REAL(SP), ALLOCATABLE :: T1TMP(:,:)
   REAL(SP), ALLOCATABLE :: ELTMP(:)
   REAL(SP), ALLOCATABLE :: UATMP(:)
   REAL(SP), ALLOCATABLE :: VATMP(:)   
   REAL(SP), ALLOCATABLE :: TMP1(:,:)
   REAL(SP), ALLOCATABLE :: TMP2(:,:)
   REAL(SP), ALLOCATABLE :: TMP3(:,:)
   REAL(SP), ALLOCATABLE :: TMP4(:)
   REAL(SP), ALLOCATABLE :: TMP5(:)
   REAL(SP), ALLOCATABLE :: TMP6(:)
#  if defined (WATER_QUALITY)
   REAL(SP), ALLOCATABLE :: TMP7(:,:,:)
#  endif

   ALLOCATE(UTMP(0:NGL,1:KB))      ; UTMP  = 0.0_SP 
   ALLOCATE(VTMP(0:NGL,1:KB))      ; VTMP  = 0.0_SP
   ALLOCATE(S1TMP(1:MGL,1:KB))     ; S1TMP = 0.0_SP
   ALLOCATE(T1TMP(1:MGL,1:KB))     ; T1TMP = 0.0_SP
   ALLOCATE(ELTMP(1:MGL))          ; ELTMP = 0.0_SP
   ALLOCATE(UATMP(0:NGL))          ; UATMP = 0.0_SP 
   ALLOCATE(VATMP(0:NGL))          ; VATMP = 0.0_SP
   ALLOCATE(TMP1(0:NGL,1:KB))      ; TMP1  = 0.0_SP
   ALLOCATE(TMP2(1:MGL,1:KB))      ; TMP2  = 0.0_SP
   ALLOCATE(TMP3(1:NGL,1:KB))      ; TMP3  = 0.0_SP
   ALLOCATE(TMP4(0:NGL))           ; TMP4  = 0.0_SP
   ALLOCATE(TMP5(1:NGL))           ; TMP5  = 0.0_SP
   ALLOCATE(TMP6(1:MGL))           ; TMP6  = 0.0_SP
#  if defined (WATER_QUALITY)
   ALLOCATE(TMP7(1:MGL,1:KB,1:NB)) ; TMP7  = 0.0_SP 
#  endif

   REWIND(INF)
   READ(INF) ITMP
   READ(INF) ((UTMP(I,K),K=1,KB),I=0,NGL)
   READ(INF) ((VTMP(I,K),K=1,KB),I=0,NGL)
   READ(INF) ((TMP1(I,K),K=1,KB),I=0,NGL)
#  if defined (GOTM)
   READ(INF) ((TMP1(I,K),K=1,KB),I=0,NGL)
   READ(INF) ((TMP1(I,K),K=1,KB),I=0,NGL)
#  else
   READ(INF) ((TMP1(I,K),K=1,KB),I=0,NGL)
   READ(INF) ((TMP1(I,K),K=1,KB),I=0,NGL)
   READ(INF) ((TMP1(I,K),K=1,KB),I=0,NGL)
#  endif
   READ(INF) ((TMP1(I,K),K=1,KB),I=0,NGL)
   READ(INF) ((TMP1(I,K),K=1,KB),I=0,NGL)
   READ(INF) ((TMP1(I,K),K=1,KB),I=0,NGL)
   READ(INF) ((TMP1(I,K),K=1,KB),I=0,NGL)
   READ(INF) ((TMP1(I,K),K=1,KB),I=0,NGL)
   READ(INF) ((TMP1(I,K),K=1,KB),I=0,NGL)

   READ(INF) ((S1TMP(I,K),K=1,KB),I=1,MGL)
   READ(INF) ((T1TMP(I,K),K=1,KB),I=1,MGL)
   READ(INF) ((TMP2(I,K),K=1,KB),I=1,MGL)
   READ(INF) ((TMP2(I,K),K=1,KB),I=1,MGL)
   READ(INF) ((TMP2(I,K),K=1,KB),I=1,MGL)
   READ(INF) ((TMP2(I,K),K=1,KB),I=1,MGL)
   READ(INF) ((TMP3(I,K),K=1,KB),I=1,NGL)
   READ(INF) ((TMP3(I,K),K=1,KB),I=1,NGL)
   READ(INF) ((TMP3(I,K),K=1,KB),I=1,NGL)

   READ(INF) (UATMP(I), I=0,NGL)
   READ(INF) (VATMP(I), I=0,NGL)
   READ(INF) (TMP5(I), I=1,NGL)
   READ(INF) (TMP5(I), I=1,NGL)
   READ(INF) (TMP5(I), I=1,NGL)
   READ(INF) (TMP5(I), I=1,NGL)
   READ(INF) (TMP5(I), I=1,NGL)
   READ(INF) (TMP5(I), I=1,NGL)

   READ(INF) (ELTMP(I), I=1,MGL)
   READ(INF) (TMP6(I), I=1,MGL)
   READ(INF) (TMP6(I), I=1,MGL)
   READ(INF) (TMP6(I), I=1,MGL)
   READ(INF) (TMP6(I), I=1,MGL)

#    if defined (EQUI_TIDE)
     READ(INF) (TMP6(I), I=1,MGL)
#    endif
#    if defined (ATMO_TIDE)
     READ(INF) (TMP6(I), I=1,MGL)
#    endif

#  if defined (WATER_QUALITY)
   DO N1=1,NB
     READ(INF) ((TMP7(I,K,N1),K=1,KB),I=1,MGL)
   END DO
#  endif

   CLOSE(INF) 
   
   IDUMMY = 0
   IF(EL_ASSIM) THEN
     DO I=1, MGL
       IDUMMY = IDUMMY + 1
       STTEMP1(IDUMMY) = ELTMP(I) 
     ENDDO
   ENDIF
   IF(UV_ASSIM) THEN
     DO K=1, KBM1
       DO I=1, NGL
         IDUMMY = IDUMMY + 1
         STTEMP1(IDUMMY) = UATMP(I) 
       ENDDO
     ENDDO
     DO K=1, KBM1
       DO I=1, NGL
         IDUMMY = IDUMMY + 1
         STTEMP1(IDUMMY) = VATMP(I)
       ENDDO
     ENDDO
   ENDIF
   IF(T_ASSIM) THEN
     DO K=1, KBM1
       DO I=1, MGL
         IDUMMY = IDUMMY + 1
         STTEMP1(IDUMMY) = T1TMP(I,K) 
       ENDDO
     ENDDO
   ENDIF
   IF(S_ASSIM) THEN
     DO K=1, KBM1
       DO I=1, MGL
         IDUMMY = IDUMMY + 1
         STTEMP1(IDUMMY) = S1TMP(I,K) 
       ENDDO
     ENDDO
   ENDIF
   
   DEALLOCATE(UTMP,VTMP,S1TMP,T1TMP,ELTMP,UATMP,VATMP,TMP1,TMP2,TMP3,TMP4,TMP5,TMP6)
#  if defined (WATER_QUALITY)
   DEALLOCATE(TMP7)
#  endif
   RETURN
   END SUBROUTINE GR2ST 

   SUBROUTINE FCT2ANL(FLNAME,GLNAME)

   USE LIMS
   USE ALL_VARS
#  if defined (WATER_QUALITY)
   USE MOD_WQM
#  endif
#  if defined (EQUI_TIDE)
   USE MOD_EQUITIDE
#  endif
#  if defined (ATMO_TIDE)
   USE MOD_ATMOTIDE
#  endif
   IMPLICIT NONE

   INTEGER I,K,N1,ITMP
   INTEGER IDUMMY
   CHARACTER(LEN=120)    :: GLNAME, FLNAME
        
   REAL(SP), ALLOCATABLE :: UTMP(:,:)
   REAL(SP), ALLOCATABLE :: VTMP(:,:)
   REAL(SP), ALLOCATABLE :: WTMP(:,:)
#  if defined (GOTM)
   REAL(SP), ALLOCATABLE :: TKETMP(:,:)
   REAL(SP), ALLOCATABLE :: TEPSTMP(:,:)
#  else   
   REAL(SP), ALLOCATABLE :: Q2TMP(:,:)
   REAL(SP), ALLOCATABLE :: Q2LTMP(:,:)
   REAL(SP), ALLOCATABLE :: LTMP(:,:)
#  endif   
   REAL(SP), ALLOCATABLE :: STMP(:,:)
   REAL(SP), ALLOCATABLE :: TTMP(:,:)
   REAL(SP), ALLOCATABLE :: RHOTMP(:,:)
   REAL(SP), ALLOCATABLE :: TMEANTMP(:,:)
   REAL(SP), ALLOCATABLE :: SMEANTMP(:,:)
   REAL(SP), ALLOCATABLE :: RMEANTMP(:,:)   
   
   REAL(SP), ALLOCATABLE :: S1TMP(:,:)
   REAL(SP), ALLOCATABLE :: T1TMP(:,:)
   REAL(SP), ALLOCATABLE :: RHO1TMP(:,:)
   REAL(SP), ALLOCATABLE :: TMEAN1TMP(:,:)
   REAL(SP), ALLOCATABLE :: SMEAN1TMP(:,:)
   REAL(SP), ALLOCATABLE :: RMEAN1TMP(:,:)    
   REAL(SP), ALLOCATABLE :: KMTMP(:,:)
   REAL(SP), ALLOCATABLE :: KHTMP(:,:)
   REAL(SP), ALLOCATABLE :: KQTMP(:,:)   
   
   REAL(SP), ALLOCATABLE :: UATMP(:)
   REAL(SP), ALLOCATABLE :: VATMP(:)
   REAL(SP), ALLOCATABLE :: EL1TMP(:)
   REAL(SP), ALLOCATABLE :: ET1TMP(:)
   REAL(SP), ALLOCATABLE :: H1TMP(:)
   REAL(SP), ALLOCATABLE :: D1TMP(:)    
   REAL(SP), ALLOCATABLE :: DT1TMP(:)
   REAL(SP), ALLOCATABLE :: RTPTMP(:)   
   
   REAL(SP), ALLOCATABLE :: ELTMP(:)
   REAL(SP), ALLOCATABLE :: ETTMP(:)
   REAL(SP), ALLOCATABLE :: HTMP(:)    
   REAL(SP), ALLOCATABLE :: DTMP(:)
   REAL(SP), ALLOCATABLE :: DTTMP(:)   

#  if defined (EQUI_TIDE)
   REAL(SP), ALLOCATABLE :: EL_EQITMP(:)
#  endif
#  if defined (ATMO_TIDE)
   REAL(SP), ALLOCATABLE :: EL_ATMOTMP(:)
#  endif

#  if defined (WATER_QUALITY)
   REAL(SP), ALLOCATABLE :: WQMTMP(:,:,:)
#  endif

   ALLOCATE(UTMP(0:NGL,1:KB))               ; UTMP = 0.0_SP
   ALLOCATE(VTMP(0:NGL,1:KB))               ; VTMP = 0.0_SP
   ALLOCATE(WTMP(0:NGL,1:KB))               ; WTMP = 0.0_SP  
#  if defined (GOTM)
   ALLOCATE(TKETMP(0:NGL,1:KB))             ; TKETMP  = 0.0_SP
   ALLOCATE(TEPSTMP(0:NGL,1:KB))            ; TEPSTMP = 0.0_SP
#  else   
   ALLOCATE(Q2TMP(0:NGL,1:KB))              ; Q2TMP  = 0.0_SP
   ALLOCATE(Q2LTMP(0:NGL,1:KB))             ; Q2LTMP = 0.0_SP
   ALLOCATE(LTMP(0:NGL,1:KB))               ; LTMP   = 0.0_SP
#  endif   
   ALLOCATE(STMP(0:NGL,1:KB))               ; STMP      = 0.0_SP
   ALLOCATE(TTMP(0:NGL,1:KB))               ; TTMP      = 0.0_SP
   ALLOCATE(RHOTMP(0:NGL,1:KB))             ; RHOTMP    = 0.0_SP
   ALLOCATE(TMEANTMP(0:NGL,1:KB))           ; TMEANTMP  = 0.0_SP
   ALLOCATE(SMEANTMP(0:NGL,1:KB))           ; SMEANTMP  = 0.0_SP
   ALLOCATE(RMEANTMP(0:NGL,1:KB))           ; RMEANTMP  = 0.0_SP

   ALLOCATE(S1TMP(1:MGL,1:KB))              ; S1TMP      = 0.0_SP
   ALLOCATE(T1TMP(1:MGL,1:KB))              ; T1TMP      = 0.0_SP
   ALLOCATE(RHO1TMP(1:MGL,1:KB))            ; RHO1TMP    = 0.0_SP
   ALLOCATE(TMEAN1TMP(1:MGL,1:KB))          ; TMEAN1TMP  = 0.0_SP
   ALLOCATE(SMEAN1TMP(1:MGL,1:KB))          ; SMEAN1TMP  = 0.0_SP
   ALLOCATE(RMEAN1TMP(1:MGL,1:KB))          ; RMEAN1TMP  = 0.0_SP
   ALLOCATE(KMTMP(1:NGL,1:KB))              ; KMTMP      = 0.0_SP
   ALLOCATE(KHTMP(1:NGL,1:KB))              ; KHTMP      = 0.0_SP
   ALLOCATE(KQTMP(1:NGL,1:KB))              ; KQTMP      = 0.0_SP  

   ALLOCATE(UATMP(0:NGL))                   ; UATMP      = 0.0_SP
   ALLOCATE(VATMP(0:NGL))                   ; VATMP      = 0.0_SP
   ALLOCATE(EL1TMP(1:NGL))                  ; EL1TMP     = 0.0_SP
   ALLOCATE(ET1TMP(1:NGL))                  ; ET1TMP     = 0.0_SP
   ALLOCATE(H1TMP(1:NGL))                   ; H1TMP      = 0.0_SP
   ALLOCATE(D1TMP(1:NGL))                   ; D1TMP      = 0.0_SP
   ALLOCATE(DT1TMP(1:NGL))                  ; DT1TMP     = 0.0_SP
   ALLOCATE(RTPTMP(1:NGL))                  ; RTPTMP     = 0.0_SP

   ALLOCATE(ELTMP(1:MGL))                   ; ELTMP      = 0.0_SP
   ALLOCATE(ETTMP(1:MGL))                   ; ETTMP      = 0.0_SP
   ALLOCATE(HTMP(1:MGL))                    ; HTMP       = 0.0_SP
   ALLOCATE(DTMP(1:MGL))                    ; DTMP       = 0.0_SP
   ALLOCATE(DTTMP(1:MGL))                   ; DTTMP      = 0.0_SP

#  if defined (EQUI_TIDE)
   ALLOCATE(EL_EQITMP(1:MGL))               ; EL_EQITMP  = 0.0_SP
#  endif
#  if defined (ATMO_TIDE)
   ALLOCATE(EL_ATMOTMP(1:MGL))              ; EL_ATMOTMP = 0.0_SP
#  endif

#  if defined (WATER_QUALITY)
   ALLOCATE(WQMTMP(1:MGL,1:KB,1:NB))        ; WQMTMP     = 0.0_SP
#  endif

   OPEN(INOKF,FILE=TRIM(FLNAME), FORM='UNFORMATTED') 

   REWIND(INOKF)
   READ(INOKF) ITMP
   READ(INOKF) ((UTMP(I,K),K=1,KB),I=0,NGL)
   READ(INOKF) ((VTMP(I,K),K=1,KB),I=0,NGL)
   READ(INOKF) ((WTMP(I,K),K=1,KB),I=0,NGL)
#  if defined (GOTM)
   READ(INOKF) ((TKETMP(I,K),K=1,KB),I=0,NGL)
   READ(INOKF) ((TEPSTMP(I,K),K=1,KB),I=0,NGL)
#  else
   READ(INOKF) ((Q2TMP(I,K),K=1,KB),I=0,NGL)
   READ(INOKF) ((Q2LTMP(I,K),K=1,KB),I=0,NGL)
   READ(INOKF) ((LTMP(I,K),K=1,KB),I=0,NGL)
#  endif
   READ(INOKF) ((STMP(I,K),K=1,KB),I=0,NGL)
   READ(INOKF) ((TTMP(I,K),K=1,KB),I=0,NGL)
   READ(INOKF) ((RHOTMP(I,K),K=1,KB),I=0,NGL)
   READ(INOKF) ((TMEANTMP(I,K),K=1,KB),I=0,NGL)
   READ(INOKF) ((SMEANTMP(I,K),K=1,KB),I=0,NGL)
   READ(INOKF) ((RMEANTMP(I,K),K=1,KB),I=0,NGL)

   READ(INOKF) ((S1TMP(I,K),K=1,KB),I=1,MGL)
   READ(INOKF) ((T1TMP(I,K),K=1,KB),I=1,MGL)
   READ(INOKF) ((RHO1TMP(I,K),K=1,KB),I=1,MGL)
   READ(INOKF) ((TMEAN1TMP(I,K),K=1,KB),I=1,MGL)
   READ(INOKF) ((SMEAN1TMP(I,K),K=1,KB),I=1,MGL)
   READ(INOKF) ((RMEAN1TMP(I,K),K=1,KB),I=1,MGL)
   READ(INOKF) ((KMTMP(I,K),K=1,KB),I=1,NGL)
   READ(INOKF) ((KHTMP(I,K),K=1,KB),I=1,NGL)
   READ(INOKF) ((KQTMP(I,K),K=1,KB),I=1,NGL)

   READ(INOKF) (UATMP(I), I=0,NGL)
   READ(INOKF) (VATMP(I), I=0,NGL)
   READ(INOKF) (EL1TMP(I), I=1,NGL)
   READ(INOKF) (ET1TMP(I), I=1,NGL)
   READ(INOKF) (H1TMP(I), I=1,NGL)
   READ(INOKF) (D1TMP(I), I=1,NGL)
   READ(INOKF) (DT1TMP(I), I=1,NGL)
   READ(INOKF) (RTPTMP(I), I=1,NGL)

   READ(INOKF) (ELTMP(I), I=1,MGL)
   READ(INOKF) (ETTMP(I), I=1,MGL)
   READ(INOKF) (HTMP(I), I=1,MGL)
   READ(INOKF) (DTMP(I), I=1,MGL)
   READ(INOKF) (DTTMP(I), I=1,MGL)

#  if defined (EQUI_TIDE)
   READ(INOKF) (EL_EQITMP(I), I=1,MGL)
#  endif
#  if defined (ATMO_TIDE)
   READ(INOKF) (EL_ATMOTMP(I), I=1,MGL)
#  endif

#  if defined (WATER_QUALITY)
   DO N1=1,NB
     READ(INOKF) ((WQMTMP(I,K,N1),K=1,KB),I=1,MGL)
   END DO
#  endif

   CLOSE(INOKF) 

   IDUMMY = 0
   IF(EL_ASSIM) THEN
      DO I=1, MGL
         IDUMMY = IDUMMY + 1
         ELTMP(I) = STTEMP1(IDUMMY)
      ENDDO
      
!      DTMP = HTMP + ELTMP
!      ETTMP = ELTMP
!      DTTMP = DTMP
      
      DO I=1, NGL
         EL1TMP(I)=(ELTMP(NVG(I,1)) + ELTMP(NVG(I,2)) + ELTMP(NVG(I,3)) )/3.0_DP
      ENDDO
   
!      D1TMP = H1TMP + EL1TMP
!      ET1TMP = EL1TMP
!      DT1TMP = D1TMP  
   
   ENDIF

   IF(UV_ASSIM) THEN
      DO K=1, KBM1
        DO I=1, NGL
          IDUMMY = IDUMMY + 1
          UATMP(I) = STTEMP1(IDUMMY)  
        ENDDO
      ENDDO  
!      UA = U(:,1)

      DO K=1, KBM1
        DO I=1, NGL
          IDUMMY = IDUMMY + 1
          VATMP(I) = STTEMP1(IDUMMY)  
        ENDDO
      ENDDO  
!      VA = V(:,1)      
   ENDIF

   IF(T_ASSIM) THEN
      DO K=1, KBM1
        DO I=1, MGL
          IDUMMY = IDUMMY + 1
          T1TMP(I,K) = STTEMP1(IDUMMY)  
        ENDDO
      ENDDO  
   ENDIF

   IF(S_ASSIM) THEN
      DO K=1, KBM1
        DO I=1, MGL
          IDUMMY = IDUMMY + 1
          S1TMP(I,K) = STTEMP1(IDUMMY)  
        ENDDO
      ENDDO  
   ENDIF

   OPEN(INOKF,FILE=TRIM(GLNAME), FORM='UNFORMATTED') 
 
   REWIND(INOKF)

   WRITE(INOKF) IEND
   WRITE(INOKF) ((UTMP(I,K),    K=1,KB),I=0,NGL)
   WRITE(INOKF) ((VTMP(I,K),    K=1,KB),I=0,NGL)
   WRITE(INOKF) ((WTMP(I,K),    K=1,KB),I=0,NGL)
#  if defined (GOTM)
   WRITE(INOKF) ((TKETMP(I,K),   K=1,KB),I=0,NGL)
   WRITE(INOKF) ((TEPSTMP(I,K),  K=1,KB),I=0,NGL)
#  else
   WRITE(INOKF) ((Q2TMP(I,K),   K=1,KB),I=0,NGL)
   WRITE(INOKF) ((Q2LTMP(I,K),  K=1,KB),I=0,NGL)
   WRITE(INOKF) ((LTMP(I,K)  ,  K=1,KB),I=0,NGL)
#  endif
   WRITE(INOKF) ((STMP(I,K),    K=1,KB),I=0,NGL)
   WRITE(INOKF) ((TTMP(I,K),    K=1,KB),I=0,NGL)
   WRITE(INOKF) ((RHOTMP(I,K),  K=1,KB),I=0,NGL)
   WRITE(INOKF) ((TMEANTMP(I,K),K=1,KB),I=0,NGL)
   WRITE(INOKF) ((SMEANTMP(I,K),K=1,KB),I=0,NGL)
   WRITE(INOKF) ((RMEANTMP(I,K),K=1,KB),I=0,NGL)

   WRITE(INOKF) ((S1TMP(I,K),    K=1,KB),I=1,MGL)
   WRITE(INOKF) ((T1TMP(I,K),    K=1,KB),I=1,MGL)
   WRITE(INOKF) ((RHO1TMP(I,K),  K=1,KB),I=1,MGL)
   WRITE(INOKF) ((TMEAN1TMP(I,K),K=1,KB),I=1,MGL)
   WRITE(INOKF) ((SMEAN1TMP(I,K),K=1,KB),I=1,MGL)
   WRITE(INOKF) ((RMEAN1TMP(I,K),K=1,KB),I=1,MGL)

   WRITE(INOKF) ((KMTMP(I,K),K=1,KB),I=1,NGL)
   WRITE(INOKF) ((KHTMP(I,K),K=1,KB),I=1,NGL)
   WRITE(INOKF) ((KQTMP(I,K),K=1,KB),I=1,NGL)

   WRITE(INOKF) (UATMP(I), I=0,NGL)
   WRITE(INOKF) (VATMP(I), I=0,NGL)

   WRITE(INOKF) (EL1TMP(I), I=1,NGL)
   WRITE(INOKF) (ET1TMP(I), I=1,NGL)
   WRITE(INOKF) (H1TMP(I),  I=1,NGL)
   WRITE(INOKF) (D1TMP(I),  I=1,NGL)
   WRITE(INOKF) (DT1TMP(I), I=1,NGL)
   WRITE(INOKF) (RTPTMP(I), I=1,NGL)

   WRITE(INOKF) (ELTMP(I), I=1,MGL)
   WRITE(INOKF) (ETTMP(I), I=1,MGL)
   WRITE(INOKF) (HTMP(I),  I=1,MGL)
   WRITE(INOKF) (DTMP(I),  I=1,MGL)
   WRITE(INOKF) (DTTMP(I), I=1,MGL)

#  if defined (EQUI_TIDE)
   WRITE(INOKF) (EL_EQITMP(I), I=1,MGL)
#  endif
#  if defined (ATMO_TIDE)
   WRITE(INOKF) (EL_ATMOTMP(I), I=1,MGL)
#  endif

#  if defined (WATER_QUALITY)
   DO N1 = 1, NB
     WRITE(INOKF) ((WQMTMP(I,K,N1),K=1,KB),I=1,MGL)
   END DO
#  endif

   CLOSE(INOKF)
   DEALLOCATE(UTMP,VTMP,WTMP,STMP,TTMP,RHOTMP,TMEANTMP,SMEANTMP,RMEANTMP,S1TMP,T1TMP,RHO1TMP)
   DEALLOCATE(TMEAN1TMP,SMEAN1TMP,RMEAN1TMP,KMTMP,KHTMP,KQTMP,UATMP,VATMP,EL1TMP,ET1TMP,H1TMP)
   DEALLOCATE(D1TMP,DT1TMP,RTPTMP,ELTMP,ETTMP,HTMP,DTMP,DTTMP)
#  if defined (EQUI_TIDE)
   DEALLOCATE(EL_EQITMP)
#  endif
#  if defined (ATMO_TIDE)
   DEALLOCATE(EL_ATMOTMP)
#  endif
#  if defined(GOTM)
   DEALLOCATE(TKETMP,TEPSTMP)
#  else
   DEALLOCATE(Q2TMP,Q2LTMP,LTMP)
#  endif    
#  if defined (WATER_QUALITY)
   DEALLOCATE(WQMTMP)
#  endif

   RETURN
   END SUBROUTINE FCT2ANL

   SUBROUTINE GETOBSLOC
   
   USE LIMS
   USE CONTROL
   IMPLICIT NONE

   INTEGER ::  NUM 
   INTEGER  SWITCH
   SAVE     SWITCH
   INTEGER ::  J,K,PNT
   INTEGER ::  IDUMMY
   INTEGER ::  TMP
   CHARACTER(LEN=24) HEADINFO
   INTEGER STLTMP(ETKF_NOBSMAX)
   INTEGER LAY(ETKF_NOBSMAX)

   NUM     = 0 
   IDUMMY  = 0
   NLOC    = 0
   PNT     = 0
   STLOC   = 0
   STLTMP  = 0 
   LAY     = 0
   STTR    = 0.0_DP
   WKTMP   = 0.0_DP

   IF(ICYC == ETKF_START/DELTA_ASS) SWITCH  = 0 

100 READ(INOOB,*,END=200) HEADINFO
    IF(SWITCH/=1) THEN
      IF(HEADINFO=='!===READ') THEN
        SWITCH = 1
        GOTO 100
      ELSE
        GOTO 100
      ENDIF
    ENDIF 
    IF(TRIM(HEADINFO)=='!STEP=') THEN
      PNT = PNT + 1
    ENDIF
    IF(PNT==2) THEN
       BACKSPACE(INOOB)
       GOTO 200
    ELSE
       IF(TRIM(HEADINFO)=='!EL') THEN
         IF(EL_OBS) THEN
           READ(INOOB,*) NUM
           NLOC = NLOC + NUM
           IF(NLOC>ETKF_NOBSMAX) THEN
             WRITE(IPT,*) 'not enough storage for observations:', 'Nloc=', Nloc, 'Nobsmax=', ETKF_NOBSMAX
             CALL PSTOP
           ENDIF
           READ(INOOB,*) (STLTMP(K), K=1,NLOC) 
           READ(INOOB,*) (STTR(STLTMP(K)), K=1,NLOC)
           DO K=1, NLOC
              STLOC(K) = STLTMP(K)     
              WKTMP(STLOC(K)) = OBSERR_EL      ! VALUE SHOULD PUT TO RUN FILE 
           ENDDO 

         ENDIF
  
         IF(EL_ASSIM) THEN
           IDUMMY = IDUMMY + MGL
         ENDIF
       ENDIF
     
       IF(TRIM(HEADINFO)=='!UV') THEN
         IF(UV_OBS) THEN
           READ(INOOB,*) NUM
           NLOC = NLOC + NUM
           IF(NLOC+NUM>ETKF_NOBSMAX) THEN
             WRITE(IPT,*) 'not enough storage for observations:', 'Nloc=', Nloc+num, 'Nobsmax=', ETKF_NOBSMAX
             CALL PSTOP
           ENDIF
           READ(INOOB,*)  (STLTMP(K), K=NLOC-NUM+1,NLOC)
           READ(INOOB,*)  (LAY(K),    K=NLOC-NUM+1,NLOC)
           DO K=NLOC-NUM+1, NLOC
              STLOC(K) = STLTMP(K) + IDUMMY + NGL*(LAY(K)-1)
           ENDDO   
           DO K=NLOC-NUM+1, NLOC
              WKTMP(STLOC(K)) = OBSERR_UV      ! VALUE SHOULD PUT TO RUN FILE
           ENDDO 
  
           NLOC = NLOC + NUM
           DO K=NLOC-NUM+1, NLOC
              STLOC(K) = STLTMP(K-NUM) + IDUMMY + NGL*KBM1 + NGL*(LAY(K-NUM)-1)
           ENDDO
           DO K=NLOC-NUM+1, NLOC
              WKTMP(STLOC(K)) = OBSERR_UV      ! VALUE SHOULD PUT TO RUN FILE
           ENDDO          
           READ(INOOB,*) (STTR(STLOC(K)),STTR(STLOC(K+NUM)), K=NLOC-2*NUM+1,NLOC-NUM) 
  
         ENDIF
         IF(UV_ASSIM) THEN
           IDUMMY = IDUMMY + 2*NGL*KBM1
         ENDIF
       ENDIF
     
       IF(TRIM(HEADINFO)=='!T') THEN
         IF(T_OBS) THEN
           READ(INOOB,*) NUM
           NLOC = NLOC + NUM
           IF(NLOC>ETKF_NOBSMAX) THEN
             WRITE(IPT,*) 'not enough storage for observations:', 'Nloc=', Nloc, 'Nobsmax=', ETKF_NOBSMAX
             CALL PSTOP
           ENDIF
           READ(INOOB,*)  (STLTMP(K), K=NLOC-NUM+1,NLOC)
           READ(INOOB,*)  (LAY(K),    K=NLOC-NUM+1,NLOC)       
           DO K=NLOC-NUM+1, NLOC
              STLOC(K) = STLTMP(K) + IDUMMY + MGL*(LAY(K)-1)
           ENDDO   
           DO K=NLOC-NUM+1, NLOC
              WKTMP(STLOC(K)) = OBSERR_T      ! FIND THE NAME OF VARIABLE
           ENDDO 
           READ(INOOB,*) (STTR(STLOC(K)), K=NLOC-NUM+1,NLOC) 

         ENDIF   
         IF(T_ASSIM) THEN
            IDUMMY = IDUMMY + MGL*KBM1
         ENDIF
       ENDIF
     
       IF(TRIM(HEADINFO)=='!S') THEN
         IF(S_OBS) THEN
           READ(INOOB,*) NUM
           NLOC = NLOC + NUM
           IF(NLOC>ETKF_NOBSMAX) THEN
             WRITE(IPT,*) 'not enough storage for observations:', 'Nloc=', Nloc, 'Nobsmax=', ETKF_NOBSMAX
             CALL PSTOP
           ENDIF
           READ(INOOB,*)  (STLTMP(K), K=NLOC-NUM+1,NLOC) 
           READ(INOOB,*)  (LAY(K),   K=NLOC-NUM+1,NLOC)       
           DO K=NLOC-NUM+1, NLOC
              STLOC(K) = STLTMP(K) + IDUMMY + MGL*(LAY(K)-1)
           ENDDO   
           DO K=NLOC-NUM+1, NLOC
              WKTMP(STLOC(K)) = OBSERR_S      ! FIND THE NAME OF VARIABLE
           ENDDO 
           READ(INOOB,*) (STTR(STLOC(K)), K=NLOC-NUM+1,NLOC) 

         ENDIF
         IF(S_ASSIM) THEN
            IDUMMY = IDUMMY + MGL*KBM1
         ENDIF
       ENDIF    
    
    ENDIF
    
    GOTO 100
200 CONTINUE

    IF(NLOC==0) THEN 
       WRITE(IPT,*) "!WARNING, NOT OBSERVATION DATA FOUND!"
       CALL PSTOP
    ELSE   
       IF(NLOC>1) THEN
         DO J=1, NLOC-1
           DO K=J+1, NLOC
             IF(STLOC(K)<STLOC(J)) THEN
               TMP = STLOC(J)
               STLOC(J) = STLOC(K)
               STLOC(K) = TMP 
             ENDIF
           ENDDO
         ENDDO
       ENDIF
    ENDIF

    RETURN
   END SUBROUTINE GETOBSLOC

   FUNCTION DISTST(STLOC1,STLOC2)

   USE LIMS
   USE CONTROL
   USE ALL_VARS
   IMPLICIT NONE

!   INTEGER NV, IDUMMY, IDUMMY1
   INTEGER IDUMMY, IDUMMY1
   INTEGER RNG(2)
   INTEGER STLOC1,STLOC2
   INTEGER STB,STE,LOC1,LOC2,K
   REAL(DP) BIG_DIST
   REAL*8   DISTST
   
   BIG_DIST = 9000000.0_DP    !1000 Km
   DISTST = BIG_DIST

   IDUMMY  = 0
   IDUMMY1 = 0
   IF(EL_ASSIM) THEN
      RNG(1)  = IDUMMY + 1
      RNG(2)  = IDUMMY + MGL 
      IDUMMY  = MGL
     
      IF(STLOC1>=RNG(1) .AND. STLOC1<=RNG(2)) THEN
         IF(STLOC2>=RNG(1) .AND. STLOC2<=RNG(2)) THEN
            DISTST = SQRT((XG(STLOC1)-XG(STLOC2))**2 + (YG(STLOC1)-YG(STLOC2))**2)
            RETURN
         ELSE
            RETURN
         ENDIF
      ENDIF
      IDUMMY1 = IDUMMY
   ENDIF

   IF(UV_ASSIM) THEN
      RNG(1) = IDUMMY + 1
      RNG(2) = IDUMMY + KBM1*NGL
      IDUMMY = IDUMMY + KBM1*NGL

      IF(STLOC1>=RNG(1) .AND. STLOC1<=RNG(2)) THEN
         IF(STLOC2>=RNG(1) .AND. STLOC2<=RNG(2)) THEN
            IF( INT((STLOC1-IDUMMY1+1)/NGL) == INT((STLOC2-IDUMMY1+1)/NGL) ) THEN           
               LOC1 = MOD(STLOC1-IDUMMY1+1,NGL)
               LOC2 = MOD(STLOC2-IDUMMY1+1,NGL)     
               DISTST = SQRT((XCG(LOC1)-XCG(LOC2))**2 + (YCG(LOC1)-YCG(LOC2))**2)
               RETURN
            ELSE
               RETURN
            ENDIF
         ELSE
            RETURN 
         ENDIF
      ENDIF
      IDUMMY1 = IDUMMY
      
      RNG(1) = IDUMMY + 1
      RNG(2) = IDUMMY + KBM1*NGL
      IDUMMY = IDUMMY + KBM1*NGL 

      IF(STLOC1>=RNG(1) .AND. STLOC1<=RNG(2)) THEN
         IF(STLOC2>=RNG(1) .AND. STLOC2<=RNG(2)) THEN
            IF( INT((STLOC1-IDUMMY1+1)/NGL) == INT((STLOC2-IDUMMY1+1)/NGL) ) THEN           
               LOC1 = MOD(STLOC1-IDUMMY1+1,NGL)
               LOC2 = MOD(STLOC2-IDUMMY1+1,NGL)     
               DISTST = SQRT((XCG(LOC1)-XCG(LOC2))**2 + (YCG(LOC1)-YCG(LOC2))**2)
               RETURN
            ELSE
               RETURN
            ENDIF
         ELSE
            RETURN 
         ENDIF
      ENDIF
      IDUMMY1 = IDUMMY
   ENDIF   

   IF(T_ASSIM) THEN
      RNG(1) = IDUMMY + 1
      RNG(2) = IDUMMY + KBM1*MGL
      IDUMMY = IDUMMY + KBM1*MGL

      IF(STLOC1>=RNG(1) .AND. STLOC1<=RNG(2)) THEN
         IF(STLOC2>=RNG(1) .AND. STLOC2<=RNG(2)) THEN
            IF( INT((STLOC1-IDUMMY1+1)/MGL) == INT((STLOC2-IDUMMY1+1)/MGL) ) THEN           
               LOC1 = MOD(STLOC1-IDUMMY1+1,MGL)
               LOC2 = MOD(STLOC2-IDUMMY1+1,MGL)     
               DISTST = SQRT((XCG(LOC1)-XCG(LOC2))**2 + (YCG(LOC1)-YCG(LOC2))**2)
               RETURN
            ELSE
               RETURN
            ENDIF
         ELSE
            RETURN 
         ENDIF
      ENDIF
      IDUMMY1 = IDUMMY
   ENDIF
      
   IF(S_ASSIM) THEN
      RNG(1) = IDUMMY + 1
      RNG(2) = IDUMMY + KBM1*MGL
      IDUMMY = IDUMMY + KBM1*MGL

      IF(STLOC1>=RNG(1) .AND. STLOC1<=RNG(2)) THEN
         IF(STLOC2>=RNG(1) .AND. STLOC2<=RNG(2)) THEN
            IF( INT((STLOC1-IDUMMY1+1)/NGL) == INT((STLOC2-IDUMMY1+1)/NGL) ) THEN           
               LOC1 = MOD(STLOC1-IDUMMY1+1,MGL)
               LOC2 = MOD(STLOC2-IDUMMY1+1,MGL)     
               DISTST = SQRT((XCG(LOC1)-XCG(LOC2))**2 + (YCG(LOC1)-YCG(LOC2))**2)
               RETURN
            ELSE
               RETURN
            ENDIF
         ELSE
            RETURN 
         ENDIF
      ENDIF
      IDUMMY1 = IDUMMY
   ENDIF  
   
   WRITE(IPT,*) '!ERROR: COULD NOT CALCULATE THE DISTST!'
   CALL PSTOP 
   
   RETURN
   END FUNCTION DISTST

   FUNCTION GASDEV(IDUM)
   IMPLICIT NONE
   
   INTEGER IDUM
   REAL GASDEV
!CU    USES ran2
!      REAL ran2
   INTEGER iset
   REAL fac,gset,rsq,v1,v2,ran2
   SAVE iset,gset
   DATA iset/0/

   if (iset.eq.0) then
1     v1=2.*ran2(idum)-1.
      v2=2.*ran2(idum)-1.
      rsq=v1**2+v2**2
      if(rsq.ge.1..or.rsq.eq.0.)goto 1
         fac=sqrt(-2.*log(rsq)/rsq)
         gset=v1*fac
         gasdev=v2*fac
         iset=1
      else
         gasdev=gset
         iset=0
      endif
   RETURN
   END FUNCTION GASDEV

   FUNCTION RAN2(IDUM)
   IMPLICIT NONE
   
   INTEGER IDUM,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
   REAL ran2,AM,EPS,RNMX
   PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,&
              IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,&
              NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
   INTEGER idum2,j,k,iv(NTAB),iy
   SAVE iv,iy,idum2
   DATA idum2/123456789/, iv/NTAB*0/, iy/0/
   
   if (idum.le.0) then
      idum=max(-idum,1)
      idum2=idum
      do 11 j=NTAB+8,1,-1
         k=idum/IQ1
         idum=IA1*(idum-k*IQ1)-k*IR1
         if (idum.lt.0) idum=idum+IM1
         if (j.le.NTAB) iv(j)=idum
11    continue
      iy=iv(1)
   endif
   k=idum/IQ1
   idum=IA1*(idum-k*IQ1)-k*IR1
   if (idum.lt.0) idum=idum+IM1
   k=idum2/IQ2
   idum2=IA2*(idum2-k*IQ2)-k*IR2
   if (idum2.lt.0) idum2=idum2+IM2
   j=1+iy/NDIV
   iy=iv(j)-idum2
   iv(j)=idum
   if(iy.lt.1)iy=iy+IMM1
   ran2=min(AM*iy,RNMX)
   return
   END FUNCTION RAN2

   SUBROUTINE PRINT_ERR(ERR,NDIM,TEXT,ERR2)
   USE MOD_PREC
   IMPLICIT NONE
      
   INTEGER NDIM,K
   REAL(DP) ERR1,ERR2,ERR3,ERR(NDIM)
   CHARACTER*15 TEXT

   ERR1 = 0.0_DP
   ERR2 = 0.0_DP
   ERR3 = 0.0_DP
   DO K = 1, NDIM
     IF (ABS(ERR(K)) > ERR1) THEN
        ERR1 = ABS(ERR(K))
     ENDIF
     ERR2 = ERR2 + ERR(K) * ERR(K)
     ERR3 = ERR3 + ABS(ERR(K))
   ENDDO
   ERR2 = DSQRT(ERR2/NDIM)
   ERR3 = ERR3 / NDIM
!   WRITE (6,'(a15,3(a10,e12.4))') text,' max err=', err1,' rms err=',err2,' mean err=', err3

   RETURN
   END SUBROUTINE PRINT_ERR

!==============================================================================
! PERTUBATE THE B.C. FOR INVERSE METHOD IN NON-JULIAN TIDAL SIMULATION
!==============================================================================
   SUBROUTINE PERT_BC

   USE LIMS
   USE CONTROL
   USE ALL_VARS
   USE BCS
   USE MOD_OBCS
#  if defined (MULTIPROCESSOR)
   USE MOD_PAR
#  endif
   IMPLICIT NONE
   
   INTEGER I,J,K,KK,ISEED
   REAL(DP)  RNBC(ETKF_NENS),RNBC_AVG,TMP
   CHARACTER(LEN=4) FNUM 
   CHARACTER(LEN=80)   :: ISTR
   CHARACTER(LEN=80)   :: COMT
   INTEGER,ALLOCATABLE :: NODE_SBC(:)
   REAL(DP),ALLOCATABLE :: APTTMP(:,:),PHAITMP(:,:)
   INTEGER   ISBCN1
   
   
   IF(IBCN_GL(1) > 0)THEN

   IF(S_TYPE == 'non-julian') THEN

     ISTR = "./"//TRIM(INPDIR)//"/"//trim(casename)
     OPEN(INOEL,FILE=TRIM(ISTR)//'_el_obc.dat')

     READ(INOEL ,1000) COMT
     READ(INOEL,*) ISBCN1
!
!-------ENSURE SAME NUMBER OF SPECIFIED OPEN BOUNDARY POINTS AS FILE-casename_obc.dat----|
!
     IF(ISBCN1 /= IBCN_GL(1))THEN
       WRITE(IPT,*)'==================ERROR=================================='
       WRITE(IPT,*)'NUMBER OF OPEN BOUNDARY POINTS IN OPEN BOUNDARY SURFACE'
       WRITE(IPT,*)'ELEVATION FILE IS LARGER THAN NUMBER OF OPEN BOUNDARY '
       WRITE(IPT,*)'POINTS OF PRESCRIBED ELEVATION TYPE IN CASENAME_obc.dat'
       WRITE(IPT,*) 'SEE SUBROUTINE BCS_FORCE'
       WRITE(IPT,*)'========================================================='
       CALL PSTOP
     END IF

!
!----READ IN BOUNDARY POINTS, AMPLITUDES, AND PHASES OF TIDE-------------------|
!
     DEALLOCATE(EMEAN,APT,PHAI)
     ALLOCATE(NODE_SBC(IBCN_GL(1)), EMEAN(IBCN_GL(1)))
     ALLOCATE(APT(IBCN_GL(1),6), PHAI(IBCN_GL(1),6))
     ALLOCATE(APTTMP(IBCN_GL(1),6), PHAITMP(IBCN_GL(1),6))
     APT = 0.0_SP ; PHAI = 0.0_SP ; EMEAN = 0.0_SP
     APTTMP = 0.0_SP; PHAITMP = 0.0_SP
     DO I=1,IBCN_GL(1)
       READ(INOEL,*)  NODE_SBC(I), EMEAN(I)
       READ (INOEL,*) (APT(I,J), J=1,6)
       READ (INOEL,*) (PHAI(I,J), J=1,6)
     END DO

     PHAI = MOD(PHAI,360.0_SP)
     APTTMP  = APT
     PHAITMP = PHAI 

     CLOSE(INOEL)

   ELSE
     WRITE(IPT,*) 'INVERSE METHOD CAN ONLY BE USED FOR NON-JULIAN TIDAL SIMULATION RIGHT NOW!'
     CALL PSTOP   
   ENDIF  
   
   ELSE
     WRITE(IPT,*) 'NO TIDAL B.C.s ARE SPECIFIED TO DO INVERSE, PLEASE CHECK AGAIN!'
     CALL PSTOP
   ENDIF 

   TIMEN = (ETKF_END - ETKF_START)/DELTA_ASS+1
   ALLOCATE(EL_SRS(IBCN_GL(1),TIMEN))    ; EL_SRS   = 0.0_DP
   ALLOCATE(TIME_SER(TIMEN))             ; TIME_SER = 0.0_DP
        
   DO I=1,TIMEN 
     TIME_SER(I)=(DTI*DBLE(ETKF_START+(I-1)*DELTA_ASS)) 
   ENDDO

!  CALL RANDOM_SEED
   RNBC_AVG = 0.0_DP
   TMP = 0.0_DP
   ISEED = -7111
   DO K=1, ETKF_NENS
      RNBC(K) = GASDEV(ISEED)
      RNBC_AVG = RNBC_AVG + RNBC(K)
   ENDDO
   RNBC_AVG = RNBC_AVG/DBLE(ETKF_NENS)
   DO K=1, ETKF_NENS 
      RNBC(K) = RNBC(K) - RNBC_AVG
      TMP = TMP + RNBC(K)
   ENDDO
   WRITE(100,*) TMP

   DO K=1, ETKF_NENS 
      WRITE(FNUM,'(I4.4)') K 

      DO I=1, IBCN_GL(1)
        DO J=1, 6
           APT(I,J)  =  APTTMP(I,J) + BC_AMP_ERR(J)*RNBC(K) 
           PHAI(I,J) = PHAITMP(I,J) + BC_PHA_ERR(J)*RNBC(K)
        ENDDO   
      ENDDO

      PHAI = MOD(PHAI,360.0_SP)

      OPEN(IOBCKF,FILE=TRIM(OUTDIR)//'/out_err/'//'el_srs'//FNUM//'.dat',STATUS='REPLACE')
      DO I = 1,IBCN_GL(1)
        DO J = 1, TIMEN 
          FORCE = 0.0_SP
          DO KK = 1, 6
            PHAI_IJ= PHAI(I,KK)*PI2/360.0_SP 
            FORCE  = FORCE + APT(I,KK)/100._DP * COS(PI2/PERIOD(KK)*(DTI*FLOAT(ETKF_START+(J-1)*DELTA_ASS)) - PHAI_IJ)
!           FORCE  = FORCE + APT(I,KK)/100._DP * COS(PI2/PERIOD(KK)*TIME_SER(J) - PHAI_IJ)
          ENDDO
          EL_SRS(I,J) = FORCE + EMEAN(I)
        END DO
      END DO
      DO J=1, TIMEN
        WRITE(IOBCKF,'(1000(F13.5))') (EL_SRS(I,J),I=1,IBCN_GL(1))
      ENDDO
      CLOSE(IOBCKF)

      OPEN(IOBCKF,FILE=TRIM(OUTDIR)//'/out_err/bc_'//FNUM//'.dat')
      DO I=1, IBCN_GL(1)
        WRITE(IOBCKF,'(I10,1000(F13.5))') NODE_SBC(I),EMEAN(I)
        WRITE(IOBCKF,'(1000(F13.5))') ( APT(I,J),J=1,6)
        WRITE(IOBCKF,'(1000(F13.5))') (PHAI(I,J),J=1,6)
      ENDDO
      CLOSE(IOBCKF) 
   ENDDO

   DEALLOCATE(NODE_SBC,APTTMP,PHAITMP)
   1000 FORMAT(A80)
   RETURN   
   END SUBROUTINE PERT_BC

   SUBROUTINE READ_BC
   
   USE LIMS
   USE CONTROL
   USE ALL_VARS
   USE BCS
   USE MOD_OBCS
#  if defined (MULTIPROCESSOR)
   USE MOD_PAR
#  endif
   IMPLICIT NONE
   
   INTEGER I,J,K
   CHARACTER(LEN=4) FNUM 
   CHARACTER(LEN=80)   :: ISTR
   CHARACTER(LEN=80)   :: COMT
   INTEGER,ALLOCATABLE :: NODE_SBC(:)
   INTEGER   ISBCN1,NCNT,JN
   REAL(SP), ALLOCATABLE :: RTEMP(:),RTEMP1(:,:),RTEMP2(:,:)
   INTEGER,  ALLOCATABLE :: TEMP2(:)
      
   IF(IENS>ETKF_NENS) THEN
      ISTR = "./"//TRIM(INPDIR)//"/"//trim(casename)//'_el_obc.dat'
      OPEN(IOBCKF,FILE=TRIM(ISTR))
      READ(IOBCKF ,*) 
      READ(IOBCKF,*) 
   ELSE
      WRITE(FNUM,'(I4.4)') IENS 
      ISTR = TRIM(OUTDIR)//'/out_err/bc_'//FNUM//'.dat'
      OPEN(IOBCKF,FILE=TRIM(ISTR))
   ENDIF
   
!----READ IN BOUNDARY POINTS, AMPLITUDES, AND PHASES OF TIDE-------------------|
!
     DEALLOCATE(EMEAN,APT,PHAI)
     ALLOCATE(NODE_SBC(IBCN_GL(1)), EMEAN(IBCN_GL(1)))
     ALLOCATE(APT(IBCN_GL(1),6), PHAI(IBCN_GL(1),6))
     APT = 0.0_SP ; PHAI = 0.0_SP ; EMEAN = 0.0_SP
     DO I=1,IBCN_GL(1)
       READ(IOBCKF,*)  NODE_SBC(I),EMEAN(I)
       READ (IOBCKF,*) (APT(I,J), J=1,6)
       READ (IOBCKF,*) (PHAI(I,J), J=1,6)
     END DO
   CLOSE(IOBCKF)  

#    if defined (MULTIPROCESSOR)
     IF(PAR)THEN
     ALLOCATE( TEMP2(IBCN_GL(1)) ,RTEMP(IBCN_GL(1)))
     ALLOCATE( RTEMP1(IBCN_GL(1),6) , RTEMP2(IBCN_GL(1),6))
     NCNT = 0
     DO I=1,IBCN_GL(1)
       IF(NLID(NODE_SBC(I)) /= 0)THEN
         NCNT = NCNT + 1
         TEMP2(NCNT)     = NLID(NODE_SBC(I))
         RTEMP(NCNT)     = EMEAN(I)
         RTEMP1(NCNT,1:6) = APT(I,1:6)
         RTEMP2(NCNT,1:6) = PHAI(I,1:6)
       END IF
     END DO

     IF(NCNT /= IBCN(1))THEN
       WRITE(IPT,*)'==================ERROR=================================='
       WRITE(IPT,*)'LOCAL OPEN BOUNDARY NODE COUNTS DIFFER BETWEEN TIDE'
       WRITE(IPT,*)'FORCING AND OPEN BOUNDARY NODE FILES'
       WRITE(IPT,*)'========================================================='
       CALL PSTOP
     END IF

!
!----TRANSFORM TO LOCAL ARRAYS-------------------------------------------------|
!
     DEALLOCATE(NODE_SBC,EMEAN,APT,PHAI)
     IF(IBCN(1) > 0)THEN
       ALLOCATE(NODE_SBC(IBCN(1)),EMEAN(IBCN(1)))
       ALLOCATE(APT(IBCN(1),6),PHAI(IBCN(1),6))
       NODE_SBC = TEMP2(1:NCNT)
       EMEAN    = RTEMP(1:NCNT)
       APT      = RTEMP1(1:NCNT,1:6)
       PHAI     = RTEMP2(1:NCNT,1:6)
     ELSE
       ALLOCATE(NODE_SBC(1),EMEAN(1))
       ALLOCATE(APT(1,6),PHAI(1,6))
       NODE_SBC = 0.0_SP ; EMEAN = 0.0_SP ; APT = 0.0_SP ; PHAI = 0.0_SP
     END IF

     DEALLOCATE(TEMP2,RTEMP,RTEMP1,RTEMP2)

     END IF !!PAR
#    endif

!
!----MAKE SURE LOCAL NODE NUMBERS OF SPECIFIED NODES MATCHES LOCAL NODE--------|
!----NUMBER OF SPECIFIED NODES IN obc.dat FILE---------------------------------|
!
     DO I=1,IBCN(1)
       JN = OBC_LST(1,I)
       IF(NODE_SBC(I) /= I_OBC_N(JN))THEN
         WRITE(IPT,*)'==================ERROR=================================='
         WRITE(IPT,*)'LOCAL OPEN BOUNDARY NODE LIST DIFFERS BETWEEN TIDE'
         WRITE(IPT,*)'FORCING AND OPEN BOUNDARY NODE (TYPE 1 OR 2) FILES'
         WRITE(IPT,*)'========================================================='
         WRITE(IPT,*)NODE_SBC(I),I_OBC_N(JN)
         CALL PSTOP
       END IF
     END DO

     APT = APT/100.0_SP
     PHAI = MOD(PHAI,360.0_SP)
   
     DEALLOCATE(NODE_SBC)
     RETURN
   END SUBROUTINE READ_BC

   SUBROUTINE INI_VALS
   USE ALL_VARS
#  if defined (MULTIPROCESSOR)
   USE MOD_PAR
#  endif
#  if defined (EQUI_TIDE)
   USE MOD_EQUITIDE
#  endif
#  if defined (ATMO_TIDE)
   USE MOD_ATMOTIDE
#  endif
   IMPLICIT NONE   
   
   U   = 0.0_SP
   V   = 0.0_SP
   W   = 0.0_SP
   UA  = 0.0_SP
   VA  = 0.0_SP
   EL1 = 0.0_SP
   ET1 = 0.0_SP
   EL  = 0.0_SP
   ET  = 0.0_SP
#  if defined (EQUI_TIDE)
   EL_EQI  = 0.0_SP
#  endif
#  if defined (ATMO_TIDE)
   EL_ATMO = 0.0_SP
#  endif
   
   
   RETURN
   END SUBROUTINE INI_VALS

   SUBROUTINE SET_INI_ETKF
     
   USE LIMS
   USE ALL_VARS
   IMPLICIT NONE

   INTEGER I,J,K
   INTEGER STDIM
   REAL(DP) SUM9
   REAL(DP),ALLOCATABLE   ::  STREF(:)
   REAL(DP),ALLOCATABLE   ::  RPETMP(:,:)
   REAL(DP),ALLOCATABLE   ::  RPATMP(:)   
   CHARACTER(LEN=120)     ::  FLNAME
   CHARACTER(LEN=4)       ::  IFIL
   
!   DO I=1, ETKF_NENS
!      WRITE(IFIL,'(I4.4)') I
!#     if !defined (DOS)
!        FLNAME = "cp "//TRIM(INPDIR)//"/"//TRIM(CASENAME)//"_restart.dat"//" "//TRIM(OUTDIR)//"/fct/restart"//IFIL//".dat"  
!#       if !defined (CRAY)
!           CALL SYSTEM( TRIM(FLNAME) )  
!#       endif
!#       if defined (CRAY)
!           CALL CRAY_SYSTEM_CALL(TRIM(FLNAME)) 
!#       endif             
!#     endif       
!
!   ENDDO

   STDIM = 0
   IF(EL_ASSIM) STDIM = STDIM + MGL
   IF(UV_ASSIM) STDIM = STDIM + 2*NGL*KBM1
   IF(T_ASSIM)  STDIM = STDIM + MGL*KBM1
   IF(S_ASSIM)  STDIM = STDIM + MGL*KBM1

   ALLOCATE(STFCT(STDIM))            ; STFCT   = ZEROD   !!state vector of one ensemble forecast
   ALLOCATE(STREF(STDIM))            ; STREF   = ZEROD  
   ALLOCATE(RPETMP(ETKF_NENS,STDIM)) ; RPETMP  = ZEROD
   ALLOCATE(RPATMP(STDIM))           ; RPATMP  = ZEROD

   IF(ETKF_CTR) THEN
      WRITE(IFIL,'(I4.4)') ETKF_NENS+1 
#     if !defined (DOS)
        FLNAME = "cp "//TRIM(INPDIR)//"/"//TRIM(CASENAME)//"_restart.dat"//" "//TRIM(OUTDIR)//"/fct/restart"//IFIL//".dat"
#       if !defined (CRAY)
           CALL SYSTEM( TRIM(FLNAME) )
#       endif
#       if defined (CRAY)
           CALL CRAY_SYSTEM_CALL(TRIM(FLNAME))
#       endif
#     endif
   ENDIF

   DO I=1, ETKF_NENS
      WRITE(IFIL,'(I4.4)') (ETKF_START-DELTA_ASS)/DELTA_ASS/ETKF_INT-5-ETKF_NENS+I
      FLNAME = TRIM(OUTDIR)//'/out/restart'//IFIL//'.dat' ! must output restart file at certain time step intervals
      OPEN(INOKF,FILE=FLNAME,FORM='UNFORMATTED')          ! and put them into an "out" folder.
      CALL GR2ST(INOKF)
      DO J=1, STDIM
         RPETMP(I,J) = STTEMP1(J)
      ENDDO
   ENDDO

   DO I=1, STDIM
     SUM9 = 0.0_DP
     DO K=1, ETKF_NENS
        SUM9 = SUM9 + RPETMP(K,I)
     ENDDO
     RPATMP(I) = SUM9/DBLE(ETKF_NENS)
     DO K=1, ETKF_NENS
        RPETMP(K,I) = RPETMP(K,I) - RPATMP(I)
     ENDDO
   ENDDO  
     
   WRITE(IFIL,'(I4.4)') (ETKF_START-DELTA_ASS)/DELTA_ASS/ETKF_INT
   FLNAME = TRIM(OUTDIR)//'/out/restart'//IFIL//'.dat'
   OPEN(INOKF,FILE=FLNAME,FORM='UNFORMATTED')
   CALL GR2ST(INOKF)
   STREF = STTEMP1

   DO K=1, ETKF_NENS
     DO I=1, STDIM
       STFCT(I) = STREF(I)*0.2_DP + RPETMP(K,I)*0.9_DP
     ENDDO
     WRITE(IFIL,'(I4.4)') K
     FLNAME = TRIM(OUTDIR)//'/fct/restart'//IFIL//'.dat'
     CALL ST2GR(FLNAME)
   ENDDO  

   DEALLOCATE(STFCT,STREF,RPETMP,RPATMP)
   RETURN 
   END SUBROUTINE SET_INI_ETKF

   SUBROUTINE ST2GR(FLNAME)

   USE LIMS
   USE ALL_VARS
#  if defined (WATER_QUALITY)
   USE MOD_WQM
#  endif 
#  if defined (EQUI_TIDE)
   USE MOD_EQUITIDE
#  endif
#  if defined (ATMO_TIDE)
   USE MOD_ATMOTIDE
#  endif
   IMPLICIT NONE

   INTEGER I,K,N1
   INTEGER IDUMMY
   REAL(SP), ALLOCATABLE :: RTP(:)
   CHARACTER(LEN=120)    :: FLNAME
    
   ALLOCATE(RTP(NGL)) ; RTP = 0.0_SP
        
   IDUMMY = 0
   IF(EL_ASSIM) THEN
      DO I=1, MGL
         IDUMMY = IDUMMY + 1
         EL(I) = STFCT(IDUMMY)
      ENDDO
      DO I=1, NGL
         EL1(I)=(EL(NVG(I,1)) + EL(NVG(I,2)) + EL(NVG(I,3)) )/3.0_DP
      ENDDO
   ENDIF

   IF(UV_ASSIM) THEN
      DO K=1, KBM1
        DO I=1, NGL
          IDUMMY = IDUMMY + 1
          U(I,K) = STFCT(IDUMMY)  
        ENDDO
      ENDDO  

      DO K=1, KBM1
        DO I=1, NGL
          IDUMMY = IDUMMY + 1
          V(I,K) = STFCT(IDUMMY)  
        ENDDO
      ENDDO  
   ENDIF

   IF(T_ASSIM) THEN
      DO K=1, KBM1
        DO I=1, MGL
          IDUMMY = IDUMMY + 1
          T1(I,K) = STFCT(IDUMMY)  
        ENDDO
      ENDDO  
   ENDIF

   IF(S_ASSIM) THEN
      DO K=1, KBM1
        DO I=1, MGL
          IDUMMY = IDUMMY + 1
          S1(I,K) = STFCT(IDUMMY)  
        ENDDO
      ENDDO  
   ENDIF

   OPEN(INOKF,FILE=TRIM(FLNAME), FORM='UNFORMATTED') 
 
   REWIND(INOKF)

   WRITE(INOKF) I_INITIAL
   WRITE(INOKF) ((U(I,K),    K=1,KB),I=0,NGL)
   WRITE(INOKF) ((V(I,K),    K=1,KB),I=0,NGL)
   WRITE(INOKF) ((W(I,K),    K=1,KB),I=0,NGL)
#  if defined (GOTM)
   WRITE(INOKF) ((TKE(I,K),   K=1,KB),I=0,NGL)
   WRITE(INOKF) ((TEPS(I,K),  K=1,KB),I=0,NGL)
#  else
   WRITE(INOKF) ((Q2(I,K),   K=1,KB),I=0,NGL)
   WRITE(INOKF) ((Q2L(I,K),  K=1,KB),I=0,NGL)
   WRITE(INOKF) ((L(I,K)  ,  K=1,KB),I=0,NGL)
#  endif
   WRITE(INOKF) ((S(I,K),    K=1,KB),I=0,NGL)
   WRITE(INOKF) ((T(I,K),    K=1,KB),I=0,NGL)
   WRITE(INOKF) ((RHO(I,K),  K=1,KB),I=0,NGL)
   WRITE(INOKF) ((TMEAN(I,K),K=1,KB),I=0,NGL)
   WRITE(INOKF) ((SMEAN(I,K),K=1,KB),I=0,NGL)
   WRITE(INOKF) ((RMEAN(I,K),K=1,KB),I=0,NGL)

   WRITE(INOKF) ((S1(I,K),    K=1,KB),I=1,MGL)
   WRITE(INOKF) ((T1(I,K),    K=1,KB),I=1,MGL)
   WRITE(INOKF) ((RHO1(I,K),  K=1,KB),I=1,MGL)
   WRITE(INOKF) ((TMEAN1(I,K),K=1,KB),I=1,MGL)
   WRITE(INOKF) ((SMEAN1(I,K),K=1,KB),I=1,MGL)
   WRITE(INOKF) ((RMEAN1(I,K),K=1,KB),I=1,MGL)

   WRITE(INOKF) ((KM(I,K),K=1,KB),I=1,NGL)
   WRITE(INOKF) ((KH(I,K),K=1,KB),I=1,NGL)
   WRITE(INOKF) ((KQ(I,K),K=1,KB),I=1,NGL)

   WRITE(INOKF) (UA(I), I=0,NGL)
   WRITE(INOKF) (VA(I), I=0,NGL)

   WRITE(INOKF) (EL1(I), I=1,NGL)
   WRITE(INOKF) (ET1(I), I=1,NGL)
   WRITE(INOKF) (H1(I),  I=1,NGL)
   WRITE(INOKF) (D1(I),  I=1,NGL)
   WRITE(INOKF) (DT1(I), I=1,NGL)
   WRITE(INOKF) (RTP(I), I=1,NGL)

   WRITE(INOKF) (EL(I), I=1,MGL)
   WRITE(INOKF) (ET(I), I=1,MGL)
   WRITE(INOKF) (H(I),  I=1,MGL)
   WRITE(INOKF) (D(I),  I=1,MGL)
   WRITE(INOKF) (DT(I), I=1,MGL)

#  if defined (EQUI_TIDE)
   WRITE(INOKF) (EL_EQI(I), I=1,M)
#  endif
#  if defined (ATMO_TIDE)
   WRITE(INOKF) (EL_ATMO(I), I=1,M)
#  endif

#  if defined (WATER_QUALITY)
   DO N1 = 1, NB
     WRITE(INOKF) ((WQM(I,K,N1),K=1,KB),I=1,MGL)
   END DO
#  endif

   CLOSE(INOKF)
   DEALLOCATE(RTP)

   RETURN
   END SUBROUTINE ST2GR

#  endif
END MODULE MOD_ETKF
