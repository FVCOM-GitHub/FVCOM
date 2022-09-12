!-----reference paper: C. Chen, Paola Malanotte-Rizzoli et al.  2009. Application and comparison of Kalman filters for coastal ocean problems: An experiment with FVCOM, J. Geophys. Res., 114, C05011, doi:10.1029/2007JC004548. 
!-----subroutine: analysis and corresponding subroutines were modified from publicly avaialbe EnKF code package of Geir Evensen.(http://enkf.nersc.no/code)
!-----subroutine: assimilate and correspnding subroutines were written based on Hunt et al. (2007): Efficient data assimilation for spatiotemporal chaos: A
!-----local sensemble transform Kalman Filter. Physisa D 
!-----Subroutine: SEIK_analysis and correspoding subroutines was modified from original code source example from Lars, Nerger: A comparison of Error Subspace Kalman Filter
!(2005),Tellus


MODULE MOD_ENKFA 
#  if defined (ENKF)
   use ENKFVAL
   USE MOD_UTILs
   USE CONTROL
   USE MOD_INPUT, ONLY : NC_DAT
   IMPLICIT NONE
   SAVE
      LOGICAL :: ENKF_TEST=.TRUE.
   real(dp),PARAMETER :: one_db=dble(1.0),zero_db=dble(0.0),two_db=dble(2.0) 
   INTEGER  ::  N1CYC =1
   INTEGER :: INNOB =657
   INTEGER  LWORK4, LDVT, RCODE

   INTEGER I_INITIAL
   REAL(DP),ALLOCATABLE :: STSD(:)      !for normalization use 
   INTEGER,ALLOCATABLE  :: STLOC(:)     !counting number in the state vector of observation location. 
   REAL(DP),ALLOCATABLE :: WKTMP(:)     !TEMP ARRAY FOR WK() 

   REAL(SP),allocatable,dimension(:)  :: x_g,y_g
   REAL(SP),allocatable,dimension(:)  :: xc_g,yc_g
   REAL(DP),ALLOCATABLE  :: WK(:,:)     !!observation error covariance matrix
                                        ! Nobs*Nobs
   REAL(DP),ALLOCATABLE  :: STFCT(:)    !!state vector of one ensemble forecast
                                        !!stfct(stdim)
   REAL(DP),ALLOCATABLE  :: AKF(:,:),TEMPAKF(:,:)    !!state vector of ensemble forecasts,
                                        !!and the difference whith it's mean
                                        !!Akf(stdim,Nens)
                                        !!Akf = Akf -Mkf  
   REAL(DP),ALLOCATABLE  :: AKF_TMP(:,:),MKF_TMP(:),MKF_TMP2(:)
   REAL(DP),ALLOCATABLE  :: MKF(:)      !!mean of state Vector forecast over 
                                        !!Nens ensembles, Mkf(stdim)
   REAL(DP),ALLOCATABLE  :: SK(:,:)     !!H*(Af-Mf),   Sk(1:Nobs,1:Nens)
   REAL(DP),ALLOCATABLE  :: MSK(:)     !!H*Mf,   MSk(1:Nobs)
   REAL(DP),ALLOCATABLE  :: HL_p(:,:)     !!H*(Af),   HL_p(1:Nobs,1:Nens)

                                        ! Tk = Tk^-1,              :Tk(1:Nobs,1:Nobs)
   REAL(DP),ALLOCATABLE  :: DUMMY(:)    !!input for subroutine gaussj,dummy(1:Nobs)
   REAL(DP),ALLOCATABLE  :: STTR1(:)
   REAL(DP),ALLOCATABLE  :: OBSDATA(:)  !!true Obs.          :obsdata(Nobs)
   REAL(DP),ALLOCATABLE  :: ERRVEC(:)   !!the difference between the mean of
                                        ! ensemble forecast and the true
                                        ! = H*(Xf-Xt)
                                        ! = Mkf(stloc(i))-sttr(stloc(i))
                                        !  :errvec(1:Nobs)
                                        ! and
                                        ! = Xf-Xt  :errvec(1:stdim)
   REAL(DP),ALLOCATABLE  :: RPERT(:,:)  !!random matrix, which mean equas zero
                                        !Rpert(Nens,Nobs)!   REAL(DP),ALLOCATABLE  :: RPERTOBS(:,:) !pengfei pertubated measurements
   REAL(DP),ALLOCATABLE  :: DDD(:,:) !pengfei pertubated measurements
   REAL(DP),ALLOCATABLE  :: RPAVE(:)
   REAL(DP),ALLOCATABLE  :: MODDATA(:)  !!1.innovation vector y' = H(x_fct)-H(x_obs)
                                        ! 2.model data in observation location 
                                        ! moddata(1:Nloc)
   REAL(DP),ALLOCATABLE ::  INNOV(:)    !pengfei
                                        ! anlysis vector for mean of ensemble
                                        ! stmean(1:stdim)

   REAL(DP),ALLOCATABLE    ::   WORK4(:)
   INTEGER  TIMEN
   REAL(DP),ALLOCATABLE  :: EL_SRS(:,:),SRS_TMP(:)
   REAL(DP),ALLOCATABLE  :: TIME_SER(:)
   REAL(SP),ALLOCATABLE  :: EL_INV(:)
   REAL(DP) AMP(6),PHAS(6),PHAI_IJ,FORCE
   logical,parameter ::  VERBOSE=.true.
   real(dp),parameter :: TRUNCATION=1.0
   logical,parameter ::  UPDATE_RANDROT=.true.
   real(dp),allocatable :: scaling(:),AKF_LOC(:,:)
        INTEGER :: STDIM_LOC,STDIM
	   REAL(DP),ALLOCATABLE  :: MKF_LOC(:)
   REAL(SP), ALLOCATABLE :: UETMP(:,:)
   REAL(SP), ALLOCATABLE :: VETMP(:,:)
   REAL(SP), ALLOCATABLE :: SETMP(:,:)
   REAL(SP), ALLOCATABLE :: TETMP(:,:)
   REAL(SP), ALLOCATABLE :: ELETMP(:)
   INTEGER,  ALLOCATABLE :: ENKFWETN(:),ENKFWETC(:)
   INTEGER, ALLOCATABLE  :: STINDX(:,:)
    INTEGER,ALLOCATABLE  :: RECALC_RHO_MEAN_ENKF(:)
   INTEGER :: MGL_OBS,NGL_OBS !total number of  node  and cell in your localized region, OBS here is confusion but the region is related to your observation data 
  Type  localization_info
      INTEGER :: cell
      INTEGER :: node
      integer :: nlay_local ! how many layers will be included in your localized region for one cell/node
      integer,allocatable :: NLAY_LOCAL_INDEX(:)  ! the index of the layers for one node/cell
  end type localization_info
  TYPE(localization_info),ALLOCATABLE :: MGL_OBS_INDEX(:),NGL_OBS_INDEX(:) ! the index of node and cell  
      
!-------------------------------------------------------------
!-----------------------------------------------------------------------------|  

   CONTAINS !-----------------------------------------------------------------|
   SUBROUTINE ALLOC_VARS_ENKF
   
   USE LIMS
   USE CONTROL
   IMPLICIT NONE  
    
   INTEGER         NDB9 
   REAL(DP)  ::    MEMCNT9 
   

   MEMCNT9 = 0.0_DP
   NDB9    = 2
   IF (ENKF_LOCALIZED) THEN
       ALLOCATE(MKF_LOC(STDIM_LOC))     ;MKF_LOC       = ZEROD
       ALLOCATE(AKF_LOC(STDIM_LOC,ENKF_NENS)) ;AKF_LOC = ZEROD 
   END IF
   ALLOCATE(MKF_TMP2(STDIM)) ; MKF_TMP2 =ZEROD
   ALLOCATE(MKF_TMP(STDIM)) ; MKF_TMP =ZEROD
   ALLOCATE(MKF(STDIM))     ;MKF       = ZEROD   !!mean of state Vector forecast over Nens ensembles
   ALLOCATE(ERRVEC(STDIM))  ;ERRVEC    = ZEROD   !!the difference between the mean of ensemble forecast and the true
                                               ! anlysis vector for mean of ensemble
               MEMCNT9 = MEMCNT9+STDIM*6*NDB9+NLOC
   ALLOCATE(DUMMY(NLOC))    ;DUMMY     = ZEROD   !!input for subroutine gaussj
   ALLOCATE(OBSDATA(NLOC))  ;OBSDATA   = ZEROD   !!true Obs.         
   ALLOCATE(RPAVE(NLOC))    ;RPAVE     = ZEROD   !! mean of perturbations
   ALLOCATE(MODDATA(NLOC))  ;MODDATA   = ZEROD   !!1.innovation vector y' = H(x_fct)-H(x_obs)
                                               ! 2.model data in observation location 
   ALLOCATE(INNOV(NLOC))    ;INNOV     = ZEROD   !pengfei

               MEMCNT9 = MEMCNT9+(ENKF_NENS+NLOC*4+NLOC)*NDB9
   ALLOCATE(STFCT(STDIM))   ; STFCT =ZEROD
   ALLOCATE(WK(NLOC,NLOC))  ;WK        = ZEROD   !!observation error covariance matrix
!   IF (.NOT. ENKF_LOCALIZED)  THEN
     ALLOCATE(AKF(STDIM,ENKF_NENS)) ;AKF = ZEROD   !!1.state vector of ensemble forecasts,
     ALLOCATE(AKF_TMP(STDIM,ENKF_NENS)) ;AKF_TMP = ZEROD 
     ALLOCATE(TEMPAKF(STDIM,ENKF_NENS)) ;TEMPAKF = ZEROD 
     ALLOCATE(STINDX(STDIM,3)) ; STINDX = 0
!   END IF                                          ! 2.and the difference whith it's mean
                                                 !   Akf = Akf -Mkf   
   
               MEMCNT9 = MEMCNT9+(NLOC*NLOC+STDIM*ENKF_NENS+ENKF_NENS*NLOC)*NDB9
   ALLOCATE(MSK(NLOC))        ;MSK    = ZEROD   !!H*Mf
   ALLOCATE(SK(NLOC,ENKF_NENS))        ;SK    = ZEROD   !!H*(Af-Mf)
   ALLOCATE(HL_p(NLOC,ENKF_NENS))        ;HL_p    = ZEROD   !!H*(Af)
               MEMCNT9 = MEMCNT9+(2*ENKF_NENS*ENKF_NENS+NLOC*ENKF_NENS+NLOC*STDIM)*NDB9
   
   ALLOCATE(RPERT(ENKF_NENS,NLOC))     ;RPERT = ZEROD   !!random matrix, which mean equas zero
   ALLOCATE(DDD(NLOC,ENKF_NENS))  ;DDD = ZEROD ! pengfei PERTURBATED MEASUREMENTS 

       MEMCNT9 = MEMCNT9+(2*NLOC*NLOC+STDIM*NLOC+ENKF_NENS*NLOC)*NDB9
  
   ALLOCATE(scaling(NLOC))  
   ALLOCATE(WORK4(LWORK4))  ; WORK4 = ZEROD
   ALLOCATE(STLOC(NLOC))   ; STLOC = ZEROD 
   ALLOCATE(WKTMP(NLOC))   ; WKTMP = ZEROD
   ALLOCATE(STTR1(STDIM))   ; STTR1 = ZEROD
   ALLOCATE(EL_INV(1:MGL))         ; EL_INV = 0.0_SP
   IF(EL_ASSIM) ALLOCATE(ELETMP(MGL))
   IF(UV_ASSIM) THEN
        ALLOCATE(UETMP(NGL,KBM1))
        ALLOCATE(VETMP(NGL,KBM1))
   END IF
   IF(T_ASSIM)  ALLOCATE(TETMP(MGL,KBM1))
   IF(S_ASSIM)  ALLOCATE(SETMP(MGL,KBM1))
 
   RETURN
   END SUBROUTINE ALLOC_VARS_ENKF


   SUBROUTINE DEALLOC_VARS_ENKF
   USE LIMS
   DEALLOCATE(MKF,ERRVEC)
   IF (ALLOCATED(MKF_TMP)) DEALLOCATE(MKF_TMP)
   IF (ALLOCATED(MKF_TMP2)) DEALLOCATE(MKF_TMP2)
   DEALLOCATE(DUMMY,OBSDATA,RPAVE,MODDATA)
   DEALLOCATE(WK,SK,HL_p,MSK)
   IF (ALLOCATED(AKF)) DEALLOCATE(AKF)
   IF (ALLOCATED(AKF_TMP)) DEALLOCATE(AKF_TMP)
   IF (ALLOCATED(TEMPAKF)) DEALLOCATE(TEMPAKF)
   IF (ALLOCATED(STINDX)) DEALLOCATE(STINDX)
   DEALLOCATE(RPERT,DDD,INNOV)
   DEALLOCATE(SCALING)
     IF (ALLOCATED(ELETMP)) DEALLOCATE(ELETMP)
     IF (ALLOCATED(UETMP)) DEALLOCATE(UETMP)
     IF (ALLOCATED(VETMP)) DEALLOCATE(VETMP)
     IF (ALLOCATED(TETMP)) DEALLOCATE(TETMP)
     IF (ALLOCATED(SETMP)) DEALLOCATE(SETMP)

#  if defined(WET_DRY)
   IF (ALLOCATED(ENKFWETN)) DEALLOCATE(ENFKWETN)
   IF (ALLOCATED(ENKFWETC)) DEALLOCATE(ENKFWETC)
#  endif
   IF (ENKF_LOCALIZED) THEN
        DEALLOCATE(MGL_OBS_INDEX,NGL_OBS_INDEX)
   END IF
   IF (ENKF_LOCALIZED) THEN
      DEALLOCATE(AKF_LOC,MKF_LOC)
   END IF
   DEALLOCATE(STFCT)
   DEALLOCATE(WORK4, STTR1, WKTMP, EL_INV,STLOC)
   RETURN
   END SUBROUTINE DEALLOC_VARS_ENKF
   
   SUBROUTINE SET_ASSIM_ENKF_EVE

   USE LIMS
   USE CONTROL
   USE ALL_VARS
   USE BCS
#  if defined (MULTIPROCESSOR)
   USE MOD_PAR
#  endif
   use mod_obcs
   IMPLICIT NONE

   INTEGER I,J,K,IERR,ios,x1,y1
   CHARACTER(LEN=100) MKANLDIR, MKFCTDIR, MKERRDIR, MKOUTDIR
   CHARACTER(LEN=120)  ::  FLNAME
   CHARACTER(LEN=4)    ::  IFIL
   CHARACTER(LEN=5)    ::  IFIL1
   
!   IF(MSR) THEN 
!#  if !defined (DOS)
!      MKANLDIR = "mkdir -p "//TRIM(OUTPUT_DIR)//"/anl"
!      MKFCTDIR = "mkdir -p "//TRIM(OUTPUT_DIR)//"/fct"
!      MKERRDIR = "mkdir -p "//TRIM(OUTPUT_DIR)//"/out_err"  
!      MKOUTDIR = "mkdir -p "//TRIM(OUTPUT_DIR)//"/flow"
!#     if !defined (CRAY)
!         CALL SYSTEM( TRIM(MKANLDIR) )  
!         CALL SYSTEM( TRIM(MKFCTDIR) )
!         CALL SYSTEM( TRIM(MKERRDIR) )
!         CALL SYSTEM( TRIM(MKOUTDIR) )
!#     endif
!#     if defined (CRAY)
!         CALL CRAY_SYSTEM_CALL(TRIM(MKANLDIR))
!         CALL CRAY_SYSTEM_CALL(TRIM(MKFCTDIR))
!         CALL CRAY_SYSTEM_CALL(TRIM(MKERRDIR))    
!         CALL CRAY_SYSTEM_CALL(TRIM(MKOUTDIR))
!#     endif             
!#  endif         
!   ENDIF

# if defined (DATA_ASSIM)
   FVCOM_RUN_MODE = FVCOM_ENKF_WITH_SSA
# else
   FVCOM_RUN_MODE = FVCOM_ENKF_WITHOUT_SSA
# endif
   allocate(RECALC_RHO_MEAN_ENKF(enkf_nens))
   RECALC_RHO_MEAN_ENKF=0 
   IF (MSR) THEN
   allocate(x_g(mgl))
   allocate(y_g(mgl))
   allocate(xc_g(ngl))
   allocate(yc_g(ngl))

      OPEN(innob,file=trim(INPUT_DIR)//trim(GRID_FILE),status='old',action='read')
        rewind(innob) 
        read(innob,*)
        read(innob,*)
        do i=1,ngl
           read(innob,*)
        end do
         I =0
    DO WHILE (.TRUE.)

       READ(innob,*,IOSTAT=IOS)J,X1,Y1
       IF(IOS<0) exit

       I = I + 1
       IF(I > MGL) THEN
          write(ipt,*) "Read ", I 
          CALL FATAL_ERROR('Number of rows of data in the grid file coordinates exceeds the stated number of nodes ?')
       END IF

#  if defined (SPHERICAL)
       IF(X1 < 0.0) X1 = X1 + 360.0_sp
#  endif

       x_g(I) = X1
       y_g(I) = Y1
    END DO

    if( I .NE. MGL) THEN
       write(ipt,*) "Read, ", I, "rows of data but mgl= ",mgl
       CALL FATAL_ERROR('Number of rows of data in the grid file coordinates does not equal the stated number of nodes ?')
    END if

   DO I=1,NGL
     XC_G(I)  = (x_g(NVG(I,1)) + x_g(NVG(I,2)) + x_g(NVG(I,3)))/3.0_SP
     YC_G(I)  = (Y_G(NVG(I,1)) + Y_G(NVG(I,2)) + Y_G(NVG(I,3)))/3.0_SP
   END DO
 
      CLOSE(innob)
      
!     OPEN(INOOB,FILE=TRIM(INPUT_DIR)//"/"//trim(casename)//"_assim_enkf.dat",FORM='FORMATTED')
     print *, 'open ensp.dat'
     OPEN(74,FILE=TRIM(OUTPUT_DIR)//'/out_err/EnSp.dat',status='replace')
     WRITE(74,*) ' Icyc       aa        infl1        infl2      inflold    inflation' 
     print *, 'open err.dat'
     OPEN(75,FILE='./err.dat',status='replace')
     print *, 'call set_in_enkf'
    CALL SET_INI_ENKF  
   END IF 
   

   RETURN
   END SUBROUTINE SET_ASSIM_ENKF_EVE

   SUBROUTINE ENKF_ASS_EVE
#  if defined(WET_DRY)
   USE MOD_WD, ONLY : WET_DRY_ON
#  endif
   USE LIMS
   USE CONTROL
   USE ALL_VARS
   IMPLICIT NONE
   
   INTEGER I, J, K, JJ, IERR, ENKF_WD
   CHARACTER(LEN=120) FNAM, GNAM    ! gnam : directory for get binary data
                                    ! fnam : directory for store binary data 
   CHARACTER(LEN=4)  FENS,cyclenumber
    CHARACTER(LEN=2) fens2   
   

   INTEGER IDOBS, IDOBS2            ! idobs,idobs2 : input number for random number generator
   REAL(DP)  AA, BB, DELT           ! delt =distst :the distance between two location
   REAL*8    DISTST
   REAL(DP)  HBHT                   ! HBHT: for localization 
   REAL(DP)  RNOBS, GASDEV          ! RNobs= gasdev : random number
   REAL(DP)  AAA,SUM0,RSCALE
   REAL(DP)  ENKF_CINF2
   CHARACTER(LEN=80) TEXT 
   REAL(DP)  SUM9, AVGRMSERR,FCTRMSERR,ANLRMSERR,FCTOBSERR,ANLOBSERR,ANLRMSERR_S,ANLRMSERR_UV
                                    ! avgrmserr : averaged rms error for all ensemble member
                                    ! fctrmserr : rms error of the mean of all ensemble memeber forcast
                                    ! anlrmserr : rms error of the mean of all ensemble memeber analysis
                                    ! fctobserr : rms error between forcast and observation
                                    ! anlobserr : rms error between analysis and boservation
   
   REAL(DP) INFLATION,INFL1,INFL2,INFL1_2,INFLOLD
   REAL(DP) VT
   REAL(DP) ELSTD, USTD, TSTD, SSTD
   INTEGER  IDUMMY 
   INTEGER ,ALLOCATABLE    ::   ISWDN(:)
   INTEGER ,ALLOCATABLE    ::   ISWDC(:)
   INTEGER ,ALLOCATABLE    ::   ISWD(:,:)
print *, '00789 all' 
if (msr) then
   IF (ENKF_LOCALIZED) THEN
      CALL GET_LOCALIZATION_REGION
   END IF
!------------------------end pengfei--------------------------------------
!----------------------------------------------------------

   LDVT   = 1
   LWORK4 = 5*ENKF_NENS
   print *, 'before allocate vars enkf'
end if !(msr)
   CALL ALLOC_VARS_ENKF

IF(PAR)CALL MPI_BARRIER(MPI_FVCOM_GROUP,IERR)
print *, 'before call gr2st_true_field'
GNAM=TRIM(OUTPUT_DIR)//"enkftemp/true_field.nc"
 CALL GR2ST_TRUE_FIELD(GNAM) !return to stfct
print *, 'after call gr2st_true_field'
IF(PAR)CALL MPI_BARRIER(MPI_FVCOM_GROUP,IERR)
if (msr) then
   print *, 'before obs2vec'
STTR1=STFCT
   if (enkf_test) then
   call getobsloc1 !! get the wktmp
   DO I=1,NLOC
      OBSDATA(I) = STTR1(STLOC(I)) !! get the obsdata 
   ENDDO
   else
   CALL OBS2VEC ! get the obsdata and wktmp
   end if
   print *, 'print y'
   OPEN(656,FILE='Y.dat')
   write(656,'(f12.6)') obsdata
   close(656)
!CC----------------------------------------------------------CC
!CC Get the Observation Covariance Matrix: Wk Nobs * Nobs    CC
!CC----------------------------------------------------------CC

   WK = 0.0_DP
   DO I=1, NLOC
      WK(I,I) = WKTMP(I)
      scaling(i)=1./wk(i,i)
      WK(I,I) = WK(I,I)**2.0_DP      
   ENDDO

   STFCT = ZEROD
print *, '008'
   MKF=ZEROD
   write(cyclenumber,'(i4.4)') icyc
END IF !(MSR)
IF(PAR)CALL MPI_BARRIER(MPI_FVCOM_GROUP,IERR)
   DO K=1,ENKF_NENS

      WRITE(FENS2,'(I2.2)') K
      GNAM=TRIM(OUTPUT_DIR)//"enkftemp/enkf_forecast_"//FENS2//".nc"
      CALL GR2ST(GNAM)        ! RETURN TO STFCT
      print *, '010'
      AKF(:,K)=STFCT(:)
   ENDDO
IF(PAR)CALL MPI_BARRIER(MPI_FVCOM_GROUP,IERR)
IF (MSR) THEN
   SK = ZEROD
     OPEN(656,FILE='AKF.dat')
     DO I=1,STDIM
        WRITE(656,'(20f12.6)') (AKF(I,K),K=1,ENKF_NENS)
     END DO
     CLOSE(656)

!CC--------Calculate the Ensemble Mean and Anomalies-------CC            
!CC the Notation is consistent with that by                CC
!CC Evensen and Van Leeuwen, MWR, 1996. P85--              CC 
     MKF= SUM(AKF,2) / DBLE(ENKF_NENS)
     DO K= 1, ENKF_NENS
        AKF(:,K) = AKF(:,K) - MKF(:)
     ENDDO
     OPEN(656,FILE='AKF_DEMEAN.dat')
     DO I=1,STDIM
        WRITE(656,'(20f12.6)') (AKF(I,K),K=1,ENKF_NENS)
     END DO
     CLOSE(656)



 write(fcyc,'(i4.4)') icyc
 CALL VEC2VAR(MKF)
   IF(UV_ASSIM) THEN 
      FNAM = TRIM(OUTPUT_DIR)//'/flow/UVm_fc'//fcyc//'.dat'
      OPEN(INOKF, FILE=TRIM(FNAM), STATUS='UNKNOWN')
      DO I=1,NGL
          WRITE(INOKF,'(2F12.4,2F12.6)') XC_G(I),YC_G(I),UETMP(I,1),VETMP(I,1)
      END DO
      CLOSE(INOKF)
   print *, 'finish flow uv'
   end if
   IF(T_ASSIM) THEN
      FNAM = TRIM(OUTPUT_DIR)//'/flow/Tm_fc'//fcyc//'.dat'
      OPEN(INOKF, FILE=TRIM(FNAM), STATUS='UNKNOWN')
      DO I=1,MGL
          WRITE(INOKF,'(2F12.4,2F12.6)') X_G(I),Y_G(I),TETMP(I,1)
      END DO
   print *, 'finish flow t'
   end if
   IF(S_ASSIM) THEN
      FNAM = TRIM(OUTPUT_DIR)//'/flow/Sm_fc'//fcyc//'.dat'
      OPEN(INOKF, FILE=TRIM(FNAM), STATUS='UNKNOWN')
      DO I=1,MGL
          WRITE(INOKF,'(2F12.4,2F12.6)') X_G(I),Y_G(I),SETMP(I,1)
      END DO
   print *, 'finish flow s '
   end if
!---------- Calculate Sk(1:Nobs,1:Nens)-------------------------C
!-----------   Sk=H_{k}*(Af_{k}-Mkf_{k})=Nobs*Nens     ---------C

   SK = ZEROD
IF (ENKF_TEST) THEN
   DO K=1,ENKF_NENS
     DO I=1,NLOC
       SK(I,K) = AKF(STLOC(I),K)
     ENDDO
   ENDDO
   DO I=1,NLOC
       MSK(I) = MKF(STLOC(I))
   ENDDO 
ELSE        
   DO K=1,ENKF_NENS
      call H_MAPPING(AKF(:,K),SK(:,K))
   ENDDO
   call H_MAPPING(MKF,MSK)
END IF
     OPEN(656,FILE='SK.dat')
     DO I=1,NLOC
        WRITE(656,'(20f12.6)') (SK(I,K),K=1,ENKF_NENS)
     END DO
     CLOSE(656)



   DO I=1,NLOC  

      ERRVEC(I) = MSK(I) - OBSDATA(I) 
   ENDDO
   TEXT='FctObserr'
   CALL PRINT_ERR(ERRVEC,NLOC,TEXT,FCTOBSERR)
      ERRVEC = MKF - STTR1
    

   TEXT='Fctrmserr'
   CALL PRINT_ERR(ERRVEC,STDIM,TEXT,FCTRMSERR)
   print *, 'finish print err'



!C-----------------------------------------------------C
!C-------------------------------------------------------C
!C---------        Assimilation  ------------------------C
!C-------------------------------------------------------C

!C--- Initialize the Random Number Generator for Obs. ---C
!C--- Initial value for IDobs must be negative        ---C
      if (ICYC == 1) then

         IDOBS  = -31
         IDOBS2 = -711
      else
        print *, 'ICYC',ICYC
        open(435,file='idobs.dat',status='old')
        read(435,*) idobs,idobs2 
        close(435)
      end if

       
!CC----------- Generate  Observations  from True State-------------CC
!CC-----------   By adding noises to true state       -------------CC

   DO I=1, NLOC
      RNOBS = GASDEV(IDOBS)
      IF (ENKF_TEST) THEN
      OBSDATA(I) = OBSDATA(I) + DSQRT(WK(I,I))*RNOBS ! represents the observation noised, should open only in twin experiments
      END IF
   ENDDO  ! Please check carefully here, it seems no reason to perturb the observation here! zlai     
       
!C---------------------------------------------------------------------C
!C---  Remove the mean of perturbations added to the observations -----C
!C---------------------------------------------------------------------C

   DO K=1, ENKF_NENS
      DO I=1, NLOC
         RNOBS = GASDEV(IDOBS2)
         RPERT(K,I)= DSQRT(WK(I,I))*RNOBS
      ENDDO
   ENDDO
   open(435,file='idobs.dat')
   write(435,*) idobs,idobs2
   close(435)
   DO I=1, NLOC
      AAA = 0.0_DP
      DO K=1, ENKF_NENS
         AAA = AAA + RPERT(K,I)
      ENDDO
      RPAVE(I) = AAA/DBLE(ENKF_NENS)
      DO K=1, ENKF_NENS
         RPERT(K,I) = RPERT(K,I) - RPAVE(I)
      ENDDO
   ENDDO

!C------------------------------------------------------C
!C------------------------------------------------------C
!C------------------------------------------------------C

!c
!c Calculate innovation vector y' = H(x_fct)-H(x_obs) ->moddata
!c
IF (ENKF_TEST) THEN
   MODDATA=MSK
ELSE
   CALL H_MAPPING(MKF,MODDATA)
END IF
     OPEN(656,FILE='MODDATA.dat')
     DO I=1,NLOC
        WRITE(656,'(20f12.6)') MODDATA(I)
     END DO
     CLOSE(656)

   
   INNOV(:)   = OBSDATA(:) - MODDATA(:)
   


   DO K=1, ENKF_NENS

!!CC----------- Perturb the observations  -------------CC
!CC----------------------------------------------------CC 
      AKF(:,K) = AKF(:,K) + MKF(:)    
IF (ENKF_TEST) THEN
            
      DO I=1,NLOC
         MODDATA(I) = AKF(STLOC(I),K)
         DDD(I,K) = OBSDATA(I)+RPERT(K,I) - MODDATA(I) !pengfei
         HL_P(i,k) = AKF(STLOC(I),K) ! H*(Akf)
      ENDDO

ELSE 
      CALL H_MAPPING(AKF(:,K),MODDATA)   
      CALL H_MAPPING(AKF(:,K),HL_P(:,k))
      DO I=1,NLOC
         DDD(I,K) = OBSDATA(I)+ RPERT(K,I)- MODDATA(I)  !RPERT(K,I)  for ensemble perturbation
      ENDDO
END IF
!----------------pengfei-------------------------    
    END DO ! ENSEMBLE
    PRINT *, 'START ANALSYSIS'
    write(cyclenumber,'(i4.4)') icyc
    AKF_TMP=AKF
    MKF_TMP=MKF
!---------------------GET LOCALIZED MATRIXS AND VECTORS-------
    IF (ENKF_LOCALIZED) THEN
       DO I=1,ENKF_NENS
          CALL LOCAL_MAPPING(AKF(:,I),AKF_LOC(:,I))
       END DO 
    END IF   

!-----------------------------------------------------------
!compare the error with true field
end if!(msr)
IF(PAR)CALL MPI_BARRIER(MPI_FVCOM_GROUP,IERR)
!print *, 'before call gr2st_true_field'
!GNAM=TRIM(OUTPUT_DIR)//"enkftemp/true_field.nc"
! CALL GR2ST_TRUE_FIELD(GNAM) !return to stfct
!print *, 'after call gr2st_true_field'
if (msr) then
      ERRVEC = MKF - STTR1

   TEXT='Fctrmserr'
   CALL PRINT_ERR(ERRVEC,STDIM,TEXT,FCTRMSERR)

!------------------------------------------------------
    IF (ENKF_LOCALIZED) THEN
    
!     IF (ENKF_METHOD==1) CALL analysis(AKF_LOC,wk,RPERT,SK,DDD,INNOV,STDIM_LOC,ENKF_NENS,NLOC,VERBOSE,TRUNCATION,MODE,UPDATE_RANDROT,scaling) !pengfei
!     IF (ENKF_METHOD==2) CALL assimilate(AKF_LOC,WK,SK,innov,STDIM_loc,ENKF_NENS,NLOC)
!     IF (ENKF_METHOD==3) CALL SEIK_analysis(AKF,WK,innov,HL_p,enkf_nens,nloc,stdim,5)
    ELSE
     IF (ENKF_METHOD==1)    CALL analysis(AKF,wk,RPERT,SK,DDD,INNOV,STDIM,ENKF_NENS,NLOC,VERBOSE,TRUNCATION,MODE,UPDATE_RANDROT,scaling) !pengfei
!     IF (ENKF_METHOD==2)    call assimilate(AKF,WK,SK,innov,STDIM,ENKF_NENS,NLOC)
!     IF (ENKF_METHOD==3)    call SEIK_analysis(AKF,WK,innov,HL_p,enkf_nens,nloc,stdim,5)
    END IF
    IF (ENKF_LOCALIZED) THEN
    DO I=1,ENKF_NENS
       CALL LOCAL_REPLACE(AKF(:,I),AKF_LOC(:,I))
    END DO 
    END IF
    MKF=ZEROD
    MKF= SUM(AKF,2) / DBLE(ENKF_NENS)
    DO K=1,ENKF_NENS   

       STFCT(:) = AKF(:,K)
       WRITE(FENS2,'(I2.2)') K
      
       GNAM=TRIM(OUTPUT_DIR)//"enkftemp/enkf_forecast_"//FENS2//".nc"        
!--------------------------------------------------------------------------|
!    read fvcom forecast for every ensembles from fct/                      |
!    write EnKF anlysis to anl/                                             |
!---------------------------------------------------------------------------|
        CALL FCT2ANL_NC(GNAM)

    END DO
        
        STFCT(:)=MKF(:)
        if (NC_OUTPUT_STACK .ne. 0) then 
           print *, 'have not built the formal flexible data output, only working for NC_OUTPUT_STACK =0 , stop...'
           call pstop
           stop
        end if
        GNAM=trim(NC_DAT%FNAME)
        print *, 'dump mean state to nc file:',trim(GNAM)
        CALL FCT2ANL_NC(GNAM) ! DUMP MEAN STATE TO NC FILE

   CALL VEC2VAR(MKF)
   IF(UV_ASSIM) THEN 
      FNAM = TRIM(OUTPUT_DIR)//'/flow/UVm_an'//fcyc//'.dat'
      OPEN(INOKF, FILE=TRIM(FNAM), STATUS='UNKNOWN')
      DO I=1,NGL
          WRITE(INOKF,'(2F12.4,2F12.6)') XC_G(I),YC_G(I),UETMP(I,1),VETMP(I,1)
      END DO
      CLOSE(INOKF)
   end if
   IF(T_ASSIM) THEN
      FNAM = TRIM(OUTPUT_DIR)//'/flow/Tm_an'//fcyc//'.dat'
      OPEN(INOKF, FILE=TRIM(FNAM), STATUS='UNKNOWN')
      DO I=1,MGL
          WRITE(INOKF,'(2F12.4,2F12.6)') X_G(I),Y_G(I),TETMP(I,1)
      END DO
   END IF
   IF(S_ASSIM) THEN
      FNAM = TRIM(OUTPUT_DIR)//'/flow/Sm_an'//fcyc//'.dat'
      OPEN(INOKF, FILE=TRIM(FNAM), STATUS='UNKNOWN')
      DO I=1,MGL
          WRITE(INOKF,'(2F12.4,2F12.6)') X_G(I),Y_G(I),SETMP(I,1)
      END DO
   END IF

      CLOSE(INOKF)
      CALL VEC2VAR(STTR1)     
      IF(UV_ASSIM) THEN
      FNAM = TRIM(OUTPUT_DIR)//'/flow/UVm_tr'//fcyc//'.dat'
      OPEN(INOKF, FILE=TRIM(FNAM), STATUS='UNKNOWN')
      DO I=1,NGL
          WRITE(INOKF,'(2F12.4,2F12.6)') XC_G(I),YC_G(I),UETMP(I,1),VETMP(I,1)
      END DO
      CLOSE(INOKF)
     end if
         IF(T_ASSIM) THEN
      FNAM = TRIM(OUTPUT_DIR)//'/flow/Tm_tr'//fcyc//'.dat'
      OPEN(INOKF, FILE=TRIM(FNAM), STATUS='UNKNOWN')
      DO I=1,MGL
          WRITE(INOKF,'(2F12.4,2F12.6)') X_G(I),Y_G(I),TETMP(I,1)
      END DO
      end if
      IF(S_ASSIM) THEN
      FNAM = TRIM(OUTPUT_DIR)//'/flow/Sm_tr'//fcyc//'.dat'
      OPEN(INOKF, FILE=TRIM(FNAM), STATUS='UNKNOWN')
      DO I=1,MGL
          WRITE(INOKF,'(2F12.4,2F12.6)') X_G(I),Y_G(I),SETMP(I,1)
      END DO
      end if
      print *, 'finish write true field'

!---------------COMPARE WITH TRUE FIELD
IF (ENKF_TEST) THEN
   DO I=1,NLOC
       MSK(I) = MKF(STLOC(I))
      OBSDATA(I) = STTR1(STLOC(I)) !! get the obsdata 
   ENDDO 
ELSE        
   call H_MAPPING(MKF,MSK)
END IF


   DO I=1,NLOC  

      ERRVEC(I) = MSK(I) - OBSDATA(I) 
   ENDDO
   TEXT='AnLObserr'
   CALL PRINT_ERR(ERRVEC,NLOC,TEXT,ANLOBSERR)
      ERRVEC = MKF - STTR1
    

   TEXT='Anlrmserr'
   CALL PRINT_ERR(ERRVEC,STDIM,TEXT,ANLRMSERR)
   print *, 'finish print err'






   print *, 'write rms error'

   print *, 'rms error before and after assimilation is:',FCTRMSERR, ANLRMSERR
   WRITE(*,*) FCTOBSERR,FCTRMSERR,ANLOBSERR,ANLRMSERR
   WRITE(75,'(4f18.8)') FCTOBSERR,FCTRMSERR,ANLOBSERR,ANLRMSERR
   print *, 'finish writing'
!--------------------------------------------------
end if !(msr)
   CALL DEALLOC_VARS_ENKF
   RETURN
   END SUBROUTINE ENKF_ASS_EVE
!------------pengfei-----------------------------------
   SUBROUTINE GET_LOCALIZATION_REGION
   IMPLICIT NONE
   INTEGER :: I,J,K,ITMP
   STDIM_LOC=0
   MGL_OBS=0
   IF (EL_ASSIM .OR. T_ASSIM .OR. S_ASSIM) THEN
      OPEN(2000,FILE=trim(INPUT_DIR)//'local_node.dat',STATUS='OLD')
      DO WHILE(.TRUE.)
         READ(2000,*,END=901) 
         MGL_OBS=MGL_OBS+1
      END DO
901   ALLOCATE(MGL_OBS_INDEX(MGL_OBS))
      REWIND(2000)
      DO I=1,MGL_OBS
         READ(2000,*) MGL_OBS_INDEX(I)%node,MGL_OBS_INDEX(I)%NLAY_LOCAL
         ALLOCATE(MGL_OBS_INDEX(I)%NLAY_LOCAL_INDEX(MGL_OBS_INDEX(I)%NLAY_LOCAL))
         BACKSPACE(2000)
         READ(2000,*) ITMP,ITMP,(MGL_OBS_INDEX(I)%NLAY_LOCAL_INDEX(J),J=1,MGL_OBS_INDEX(I)%NLAY_LOCAL)
      END DO
   CLOSE(2000)
   END IF 
   IF (UV_ASSIM) THEN   
   NGL_OBS=0
   OPEN(2000,FILE=trim(INPUT_DIR)//'local_cell.dat',status='old')
      DO WHILE(.TRUE.)
         READ(2000,*,END=902)
         NGL_OBS=NGL_OBS+1
      END DO
902   ALLOCATE(NGL_OBS_INDEX(NGL_OBS))
      REWIND(2000)
      DO I=1,NGL_OBS
         READ(2000,*) NGL_OBS_INDEX(I)%cell,NGL_OBS_INDEX(I)%NLAY_LOCAL
         ALLOCATE(NGL_OBS_INDEX(I)%NLAY_LOCAL_INDEX(NGL_OBS_INDEX(I)%NLAY_LOCAL))
         BACKSPACE(2000)
         READ(2000,*) ITMP,ITMP,(NGL_OBS_INDEX(I)%NLAY_LOCAL_INDEX(J),J=1,NGL_OBS_INDEX(I)%NLAY_LOCAL)         
      END DO
   CLOSE(2000)
   END IF
   IF(EL_ASSIM) STDIM_LOC = STDIM_LOC + MGL_OBS
   IF(UV_ASSIM) THEN
      DO I=1,NGL_OBS
         STDIM_LOC = STDIM_LOC + 2*NGL_OBS_INDEX(I)%NLAY_LOCAL
      END DO
   END IF
   IF(T_ASSIM)  THEN
      DO I=1,MGL_OBS
         STDIM_LOC = STDIM_LOC + MGL_OBS_INDEX(I)%NLAY_LOCAL
      END DO    
   END IF  
   IF(S_ASSIM)  THEN
      DO I=1,MGL_OBS
         STDIM_LOC = STDIM_LOC + MGL_OBS_INDEX(I)%NLAY_LOCAL
      END DO
   END IF
   END SUBROUTINE GET_LOCALIZATION_REGION
!------------------------
   SUBROUTINE OBS2VEC
   use mod_enkf_obs
   IMPLICIT NONE
   REAL(DP) ::H_VEC(NLOC)
   REAL(DP):: VEC(STDIM)
   real    :: var1,var2,w1,w2
   INTEGER :: I,K,NLAY,IDUMMY
   call VEC2VAR(VEC)
   IDUMMY = 0
   IF(EL_ASSIM) THEN
     DO I=1,N_ASSIM_EL
       NLAY      = EL0_OBS(I)%N_LAYERS
       if (nlay .ne. 1) then
           print *, 'impossible nlay for elevation is ',nlay
           call pstop
       end if
       DO K=1,NLAY
         IDUMMY = IDUMMY + 1
         OBSDATA(IDUMMY)= EL0_OBS(I)%EL(I,K)
         WKTMP(IDUMMY)= OBSERR_EL
       END DO
     END DO
    END IF

   IF(UV_ASSIM) THEN
     DO I=1,N_ASSIM_CUR
       NLAY      = CUR_OBS(I)%N_LAYERS
       DO K=1,NLAY
         IDUMMY = IDUMMY + 1
         OBSDATA(IDUMMY)= CUR_OBS(I)%UO(1,K) ! only one time , i.e. current time
         WKTMP(IDUMMY)= OBSERR_UV
       END DO
     END DO
     DO I=1,N_ASSIM_CUR
       NLAY      = CUR_OBS(I)%N_LAYERS
       DO K=1,NLAY
         IDUMMY = IDUMMY + 1
         OBSDATA(IDUMMY)= CUR_OBS(I)%VO(1,K) ! only one time , i.e. current time
         WKTMP(IDUMMY)= OBSERR_UV
       END DO
     END DO
   END IF



   IF(T_ASSIM) THEN

     DO I=1,N_ASSIM_T
       NLAY      = T0_OBS(I)%N_LAYERS
       DO K=1,NLAY
         IDUMMY = IDUMMY + 1
         OBSDATA(IDUMMY)= T0_OBS(I)%TEMP(i,K) !only one time , i.e. current time
         WKTMP(IDUMMY)= OBSERR_T
       END DO
     END DO
   END IF


   IF(S_ASSIM) THEN
     DO I=1,N_ASSIM_S
       NLAY      = S0_OBS(I)%N_LAYERS
       DO K=1,NLAY
         IDUMMY = IDUMMY + 1
         OBSDATA(IDUMMY)= S0_OBS(I)%SAL(1,K) !only one time , i.e. current time
         WKTMP(IDUMMY)= OBSERR_S
       END DO
     END DO
   END IF
   END SUBROUTINE OBS2VEC
!====================================================
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

   RETURN
   END SUBROUTINE PRINT_ERR
!================================================================
   SUBROUTINE LOCAL_REPLACE(VEC,H_VEC) ! GIVEN A VEC, REPLACE ITS CORRESPONDING ELEMENTS TO H_VEC WITH H_VEC
   use mod_enkf_obs
   IMPLICIT NONE
   REAL(DP) ::H_VEC(NLOC)
   REAL(DP):: VEC(STDIM)
   real    :: var1,var2,w1,w2
   INTEGER :: I,K,NLAY,IDUMMY
   call VEC2VAR(VEC)
   IDUMMY = 0
   IF(EL_ASSIM) THEN
     DO I=1,MGL_OBS
         IDUMMY = IDUMMY + 1
         ELETMP(NGL_OBS_INDEX(I)%node)=H_VEC(IDUMMY)
     END DO
   END IF

   IF(UV_ASSIM) THEN
     DO I=1,NGL_OBS
       DO K=1,NGL_OBS_INDEX(I)%NLAY_LOCAL
         IDUMMY = IDUMMY + 1
         UETMP(NGL_OBS_INDEX(I)%cell,NGL_OBS_INDEX(I)%NLAY_LOCAL_INDEX(K))=H_VEC(IDUMMY)
       END DO
     END DO
     DO I=1,NGL_OBS
       DO K=1,NGL_OBS_INDEX(I)%NLAY_LOCAL
         IDUMMY = IDUMMY + 1
         VETMP(NGL_OBS_INDEX(I)%cell,NGL_OBS_INDEX(I)%NLAY_LOCAL_INDEX(K))=H_VEC(IDUMMY)
       END DO
     END DO
   END IF



   IF(T_ASSIM) THEN
     DO I=1,MGL_OBS
       DO K=1,MGL_OBS_INDEX(I)%NLAY_LOCAL
         IDUMMY = IDUMMY + 1
         TETMP(MGL_OBS_INDEX(I)%node,MGL_OBS_INDEX(I)%NLAY_LOCAL_INDEX(K))=H_VEC(IDUMMY)  
       END DO
     END DO
   END IF


   IF(S_ASSIM) THEN
     DO I=1,MGL_OBS
       DO K=1,MGL_OBS_INDEX(I)%NLAY_LOCAL
         IDUMMY = IDUMMY + 1
         TETMP(MGL_OBS_INDEX(I)%node,MGL_OBS_INDEX(I)%NLAY_LOCAL_INDEX(K))=H_VEC(IDUMMY)
       END DO
     END DO
   END IF
   CALL VAR2VEC(VEC)
   END SUBROUTINE LOCAL_REPLACE
!================================================================
   SUBROUTINE LOCAL_MAPPING(VEC,H_VEC)
   use mod_enkf_obs
   IMPLICIT NONE
   REAL(DP) ::H_VEC(NLOC)
   REAL(DP):: VEC(STDIM)
   INTEGER :: I,K,NLAY,IDUMMY
   call VEC2VAR(VEC)
   IDUMMY = 0
   IF(EL_ASSIM) THEN
     DO I=1,MGL_OBS
         IDUMMY = IDUMMY + 1
         H_VEC(IDUMMY)= ELETMP(NGL_OBS_INDEX(I)%node)
     END DO
    END IF

   IF(UV_ASSIM) THEN
     DO I=1,NGL_OBS
       DO K=1,NGL_OBS_INDEX(I)%NLAY_LOCAL
         IDUMMY = IDUMMY + 1
         H_VEC(IDUMMY)=  UETMP(NGL_OBS_INDEX(I)%cell,NGL_OBS_INDEX(I)%NLAY_LOCAL_INDEX(K))
       END DO
     END DO
     DO I=1,NGL_OBS
       DO K=1,NGL_OBS_INDEX(I)%NLAY_LOCAL
         IDUMMY = IDUMMY + 1
         H_VEC(IDUMMY)=  VETMP(NGL_OBS_INDEX(I)%cell,NGL_OBS_INDEX(I)%NLAY_LOCAL_INDEX(K))
       END DO
     END DO
   END IF



   IF(T_ASSIM) THEN
     DO I=1,MGL_OBS
       DO K=1,MGL_OBS_INDEX(I)%NLAY_LOCAL
         IDUMMY = IDUMMY + 1
         H_VEC(IDUMMY)=  TETMP(MGL_OBS_INDEX(I)%node,MGL_OBS_INDEX(I)%NLAY_LOCAL_INDEX(K))
       END DO
     END DO
   END IF


   IF(S_ASSIM) THEN
     DO I=1,MGL_OBS
       DO K=1,MGL_OBS_INDEX(I)%NLAY_LOCAL
         IDUMMY = IDUMMY + 1
         H_VEC(IDUMMY)=  TETMP(MGL_OBS_INDEX(I)%node,MGL_OBS_INDEX(I)%NLAY_LOCAL_INDEX(K))
       END DO
     END DO
   END IF
   END SUBROUTINE LOCAL_MAPPING
!------------pengfei-----------------------------------




   SUBROUTINE H_MAPPING(VEC,H_VEC)
   use mod_enkf_obs
   IMPLICIT NONE
   REAL(DP) ::H_VEC(NLOC)
   REAL(DP):: VEC(STDIM)
   real    :: var1,var2,w1,w2
   INTEGER :: I,K,NLAY,IDUMMY
   H_VEC=ZERO_DB
   call VEC2VAR(VEC)
   IDUMMY = 0
   IF(EL_ASSIM) THEN
     DO I=1,N_ASSIM_EL
       NLAY      = EL0_OBS(I)%N_LAYERS
       if (nlay .ne. 1) then
           print *, 'impossible nlay for elevation is ',nlay
           call pstop
       end if
       DO K=1,NLAY
         IDUMMY = IDUMMY + 1
         var1 = ELETMP(EL0_OBS(I)%INNODE)
!         var2 = ELETMP(EL0_OBS(I)%INNODE)
!         W1 = EL0_OBS(I)%S_WEIGHT(K,1)
!         W2 = EL0_OBS(I)%S_WEIGHT(K,2)
         H_VEC(IDUMMY)= var1 !*W1 + var2*W2
       END DO
     END DO
    END IF

   IF(UV_ASSIM) THEN
     DO I=1,N_ASSIM_CUR
       NLAY      = CUR_OBS(I)%N_LAYERS
       DO K=1,NLAY
         IDUMMY = IDUMMY + 1
         var1 = UETMP(CUR_OBS(I)%N_CELL,CUR_OBS(I)%S_INT(K,1))
         var2 = UETMP(CUR_OBS(I)%N_CELL,CUR_OBS(I)%S_INT(K,2))
         W1 =   CUR_OBS(I)%S_WEIGHT(K,1)
         W2 =   CUR_OBS(I)%S_WEIGHT(K,2)
         H_VEC(IDUMMY)= var1*W1 + var2*W2
       END DO
     END DO
     DO I=1,N_ASSIM_CUR
       NLAY      = CUR_OBS(I)%N_LAYERS
       DO K=1,NLAY
         IDUMMY = IDUMMY + 1
         var1 = VETMP(CUR_OBS(I)%N_CELL,CUR_OBS(I)%S_INT(K,1))
         var2 = VETMP(CUR_OBS(I)%N_CELL,CUR_OBS(I)%S_INT(K,2))
         W1 =   CUR_OBS(I)%S_WEIGHT(K,1)
         W2 =   CUR_OBS(I)%S_WEIGHT(K,2)
         H_VEC(IDUMMY)= var1*W1 + var2*W2
       END DO
     END DO
   END IF



   IF(T_ASSIM) THEN

     DO I=1,N_ASSIM_T
       NLAY      = T0_OBS(I)%N_LAYERS
       DO K=1,NLAY
         IDUMMY = IDUMMY + 1
         var1 = TETMP(T0_OBS(I)%INNODE,T0_OBS(I)%S_INT(K,1))
         var2 = TETMP(T0_OBS(I)%INNODE,T0_OBS(I)%S_INT(K,2))
         W1 =   T0_OBS(I)%S_WEIGHT(K,1)
         W2 =   T0_OBS(I)%S_WEIGHT(K,2)
         H_VEC(IDUMMY)= var1*W1 + var2*W2
       END DO
     END DO
   END IF


   IF(S_ASSIM) THEN
     DO I=1,N_ASSIM_S
       NLAY      = S0_OBS(I)%N_LAYERS
       DO K=1,NLAY
         IDUMMY = IDUMMY + 1
         var1 = SETMP(S0_OBS(I)%INNODE,S0_OBS(I)%S_INT(K,1))
         var2 = SETMP(S0_OBS(I)%INNODE,S0_OBS(I)%S_INT(K,2))
         W1 =   S0_OBS(I)%S_WEIGHT(K,1)
         W2 =   S0_OBS(I)%S_WEIGHT(K,2)
         H_VEC(IDUMMY)= var1*W1 + var2*W2
       END DO
     END DO
   END IF
   END SUBROUTINE H_MAPPING
!------------pengfei-----------------------------------
   SUBROUTINE VAR2VEC(VEC)
   USE LIMS
   USE CONTROL
   IMPLICIT NONE
   INTEGER :: IDUMMY,I,K
   REAL(DP):: VEC(STDIM)
!   ALLOCATE(VEC(STDIM))
   IDUMMY = 0
   IF(EL_ASSIM) THEN
      DO I=1, MGL
         IDUMMY = IDUMMY + 1
         VEC(IDUMMY) = ELETMP(I)
      ENDDO

   ENDIF

   IF(UV_ASSIM) THEN
      DO K=1, KBM1
        DO I=1, NGL
          IDUMMY = IDUMMY + 1
          VEC(IDUMMY) = UETMP(I,K) 
        ENDDO
      ENDDO


      DO K=1, KBM1
        DO I=1, NGL
          IDUMMY = IDUMMY + 1
          VEC(IDUMMY) = VETMP(I,K)
        ENDDO
      ENDDO

   ENDIF

   IF(T_ASSIM) THEN
      DO K=1, KBM1
        DO I=1, MGL
          IDUMMY = IDUMMY + 1
          VEC(IDUMMY) = TETMP(I,K)
        ENDDO
      ENDDO
   ENDIF

   IF(S_ASSIM) THEN
      DO K=1, KBM1
        DO I=1, MGL
          IDUMMY = IDUMMY + 1
          VEC(IDUMMY) = SETMP(I,K)
        ENDDO
      ENDDO
   ENDIF
!   DEALLOCATE(VEC)
   END SUBROUTINE VAR2VEC

!------------------------------------------------------

!------------------------------------------------------
   SUBROUTINE VEC2VAR(VEC)
   USE LIMS
   USE CONTROL
   IMPLICIT NONE
   INTEGER :: IDUMMY,I,K
   REAL(DP):: VEC(STDIM)
!   ALLOCATE(VEC(STDIM))
   IDUMMY = 0
   IF(EL_ASSIM) THEN
      DO I=1, MGL
         IDUMMY = IDUMMY + 1
         ELETMP(I) = VEC(IDUMMY)
      ENDDO

   ENDIF

   IF(UV_ASSIM) THEN
      DO K=1, KBM1
        DO I=1, NGL
          IDUMMY = IDUMMY + 1
          UETMP(I,K) = VEC(IDUMMY)
        ENDDO
      ENDDO


      DO K=1, KBM1
        DO I=1, NGL
          IDUMMY = IDUMMY + 1
          VETMP(I,K) = VEC(IDUMMY)
        ENDDO
      ENDDO

   ENDIF

   IF(T_ASSIM) THEN
      DO K=1, KBM1
        DO I=1, MGL
          IDUMMY = IDUMMY + 1
          TETMP(I,K) = VEC(IDUMMY)
        ENDDO
      ENDDO
   ENDIF

   IF(S_ASSIM) THEN
      DO K=1, KBM1
        DO I=1, MGL
          IDUMMY = IDUMMY + 1
          SETMP(I,K) = VEC(IDUMMY)
        ENDDO
      ENDDO
   ENDIF
!   DEALLOCATE(VEC)
   END SUBROUTINE VEC2VAR

!---------------------------------------------------------
   SUBROUTINE ENKF_READ_FIELD(GNAM,ELLG,ULG,VLG,T1LG,S1LG)
   USE CONTROL
   USE ENKFVAL
   USE MOD_NCTOOLS
   IMPLICIT NONE
   CHARACTER(len=*) :: GNAM
   TYPE(NCFILE), POINTER :: NCF
   TYPE(TIME) :: tTEST
   LOGICAL :: FOUND
   INTEGER :: NTIMES
   

   TYPE(NCDIM), POINTER :: DIM
   TYPE(NCVAR), POINTER :: VAR

   INTEGER :: nrange,IRANGE(2),CNT, STATUS
   integer :: istkcnt
   character(len=120) :: pathnfile 
  REAL(SP), allocatable,DIMENSION(:,:)   :: UL,VL
  REAL(SP), allocatable,DIMENSION(:,:)  :: T1L,S1L
  REAL(SP), allocatable, DIMENSION(:)      :: ELL

  REAL(SP), DIMENSION(0:NGL,KB),INTENT(OUT)   :: ULG,VLG
  REAL(SP), DIMENSION(0:MGL,KB),INTENT(OUT)   :: T1LG,S1LG
  REAL(SP), DIMENSION(0:MGL),INTENT(OUT)      :: ELLG


  ELLG=0.0
  ULG=0.0
  VLG=0.0
  T1LG=0.0
  S1LG=0.0
   



   NCF => NEW_FILE()
   NCF%FNAME=trim(GNAM)
   
   Call NC_OPEN(NCF)
   CALL NC_LOAD(NCF)

   ! GET THE TIME DIMENSION
   DIM => FIND_UNLIMITED(NCF,FOUND)
   IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN ENKF_READ_FIELD FILE",&
            & "FILE NAME: "//TRIM(NCF%FNAME),&
            &"COULD NOT FIND THE UNLIMITED DIMENSION")
   NTIMES = DIM%DIM
   

   ! TEST THE FILE START AND END
   TTEST = get_file_time(NCF,1)
   
   IF(INTTIME < TTEST) CALL FATAL_ERROR &
        &("ENKF_READ_FIELD: INTTIME IS BEFORE FIRST TIME IN FILE.",&
        & "FILE NAME:"//TRIM(NCF%FNAME))
   
   
   TTEST = get_file_time(NCF,NTIMES)
   IF(INTTIME > TTEST) CALL FATAL_ERROR &
        &("ENKD_READ_FIELD: INTTIME IS AFTER LAST TIME IN FILE.",&
        & "FILE NAME:"//TRIM(NCF%FNAME))


   ! GET THE START AND END INDEX

   IRANGE = 0
   istkcnt=0
   CNT = 0
   DO WHILE (product(irange)==0)
      CNT = CNT +1

      IF (CNT > NTIMES) THEN
         CALL PRINT_TIME(INTTIME,IPT,"INTTIME")
         

         CALL FATAL_ERROR&
           &("ENKF_READ_FIELD: CAN NOT FIND INTTIME IN THE THIS FILE!", &
           & "FILE NAME:"//TRIM(NCF%FNAME))
      END IF

    
      TTEST = get_file_time(NCF,CNT)
      if (msr)      print *, 'ttest',ttest
      if (msr) print *, 'inttime',inttime
      IF (INTTIME == TTEST) THEN
         IRANGE(1) = CNT
         istkcnt=cnt
      END IF

      IF (inttime == TTEST) THEN
         IRANGE(2) = CNT
      END IF

   END DO


   nrange = irange(2)-irange(1)+1

   IF(EL_ASSIM) THEN
      
      VAR => FIND_VAR(NCF,"zeta",FOUND)
      IF(.not. FOUND) CALL FATAL_ERROR &
           & ("IN ENKF_READ_FIELD FILE",&
            & "FILE NAME: "//TRIM(NCF%FNAME),&
            &"COULD NOT FIND THE VARIABLE 'zeta'")
      
      allocate(ELL(1:MGL))

      CALL nc_connect_avar(VAR,ELL)

      CALL NC_READ_VAR(var,stkCNT=istkcnt,parallel=.false.)
      ELLG(1:MGL)=ELL
      deallocate(ELL)
   ENDIF
   
   IF(UV_ASSIM) THEN

      VAR => FIND_VAR(NCF,"u",FOUND)
      IF(.not. FOUND) CALL FATAL_ERROR &
           & ("IN ENKF_READ_FIELD FILE",&
            & "FILE NAME: "//TRIM(NCF%FNAME),&
            &"COULD NOT FIND THE VARIABLE 'u'")
      allocate(UL(1:NGL,KBM1))

      CALL nc_connect_avar(VAR,UL)

      CALL NC_READ_VAR(var,stkcnt=istkcnt,parallel=.false.)

      ULG(1:NGL,1:KBM1)=UL
      deallocate(UL)
      VAR => FIND_VAR(NCF,"v",FOUND)
      IF(.not. FOUND) CALL FATAL_ERROR &
           & ("IN ENKF_READ_FIELD FILE",&
            & "FILE NAME: "//TRIM(NCF%FNAME),&
            &"COULD NOT FIND THE VARIABLE 'v'")

      allocate(VL(1:NGL,KBM1))
      CALL nc_connect_avar(VAR,VL)
      CALL NC_READ_VAR(var,stkcnt=istkcnt,parallel=.false.)
      VLG(1:NGL,1:KBM1)=VL
      deallocate(VL)
   ENDIF
   
   IF(T_ASSIM) THEN
      VAR => FIND_VAR(NCF,"temp",FOUND)
      IF(.not. FOUND) CALL FATAL_ERROR &
           & ("IN ENKF_READ_FIELD FILE",&
            & "FILE NAME: "//TRIM(NCF%FNAME),&
            &"COULD NOT FIND THE VARIABLE 'temp'")

      allocate(T1L(1:MGL,KBM1))
      CALL nc_connect_avar(VAR,T1L)

      CALL NC_READ_VAR(var,stkcnt=istkcnt,parallel=.false.)
      T1LG(1:MGL,1:KBM1)=T1L
      deallocate(T1L)
   ENDIF
   
   IF(S_ASSIM) THEN
      VAR => FIND_VAR(NCF,"salinity",FOUND)
      IF(.not. FOUND) CALL FATAL_ERROR &
           & ("IN ENKF_READ_FIELD  FILE",&
            & "FILE NAME: "//TRIM(NCF%FNAME),&
            &"COULD NOT FIND THE VARIABLE 'salinity'")
      allocate(S1L(1:MGL,KBM1))
      CALL nc_connect_avar(VAR,S1L)

      CALL NC_READ_VAR(var,stkcnt=istkcnt,parallel=.false.)
      S1LG(1:MGL,1:KBM1)=S1L
      deallocate(S1L)
   ENDIF
   
   CALL KILL_FILE(NCF)
   
   
 
 END SUBROUTINE ENKF_READ_FIELD
!--------------end pengfei-------------------------------

   SUBROUTINE GR2ST(GNAM)
   
   USE LIMS
   USE ALL_VARS
#  if defined (WATER_QUALITY)
   USE MOD_WQM
#  endif
!#  if defined (EQUI_TIDE)
!   USE MOD_EQUITIDE
!#  endif
!#  if defined (ATMO_TIDE)
!   USE MOD_ATMOTIDE
!#  endif
   IMPLICIT NONE
   CHARACTER(len=120) :: GNAM
!   CHARACTER(len=8) :: chint
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
   REAL(SP), ALLOCATABLE :: TMP8(:,:)
#  if defined (WATER_QUALITY)
   REAL(SP), ALLOCATABLE :: TMP7(:,:,:)
#  endif

   ALLOCATE(UTMP(0:NGL,1:KB))      ; UTMP  = 0.0_SP 
   ALLOCATE(VTMP(0:NGL,1:KB))      ; VTMP  = 0.0_SP
   ALLOCATE(S1TMP(0:MGL,1:KB))     ; S1TMP = 0.0_SP
   ALLOCATE(T1TMP(0:MGL,1:KB))     ; T1TMP = 0.0_SP
   ALLOCATE(ELTMP(0:MGL))          ; ELTMP = 0.0_SP
   ALLOCATE(UATMP(0:NGL))          ; UATMP = 0.0_SP 
   ALLOCATE(VATMP(0:NGL))          ; VATMP = 0.0_SP
   ALLOCATE(TMP1(0:NGL,1:KB))      ; TMP1  = 0.0_SP
   ALLOCATE(TMP2(0:MGL,1:KB))      ; TMP2  = 0.0_SP
   ALLOCATE(TMP8(0:MGL,1:KB))      ; TMP8  = 0.0_SP
   ALLOCATE(TMP3(0:NGL,1:KB))      ; TMP3  = 0.0_SP
   ALLOCATE(TMP4(0:NGL))           ; TMP4  = 0.0_SP
   ALLOCATE(TMP5(0:NGL))           ; TMP5  = 0.0_SP
   ALLOCATE(TMP6(0:MGL))           ; TMP6  = 0.0_SP
#  if defined (WATER_QUALITY)
   ALLOCATE(TMP7(1:MGL,1:KB,1:NB)) ; TMP7  = 0.0_SP 
#  endif
!   call NCD_READ_ENKF(trim(GNAM),UTMP,VTMP,T1TMP,S1TMP,ELTMP,0)
    CALL ENKF_READ_FIELD(trim(GNAM),ELTMP,UTMP,VTMP,T1TMP,S1TMP)
IF (MSR) THEN
   IDUMMY = 0
   IF(EL_ASSIM) THEN
     DO I=1, MGL
       IDUMMY = IDUMMY + 1
       STFCT(IDUMMY) = ELTMP(I)
       STINDX(IDUMMY,1)=1
       STINDX(IDUMMY,2)=I
       STINDX(IDUMMY,3)=1 
     ENDDO
   ENDIF
   
   IF(UV_ASSIM) THEN
     DO K=1, KBM1
       DO I=1, NGL
         IDUMMY = IDUMMY + 1
         STFCT(IDUMMY) = UTMP(I,k)
         STINDX(IDUMMY,1)=2
         STINDX(IDUMMY,2)=I
         STINDX(IDUMMY,3)=K 
       ENDDO
     ENDDO
     DO K=1, KBM1
       DO I=1, NGL
         IDUMMY = IDUMMY + 1
         STFCT(IDUMMY) = VTMP(I,k)
         STINDX(IDUMMY,1)=2
         STINDX(IDUMMY,2)=I
         STINDX(IDUMMY,3)=K
       ENDDO
     ENDDO
   ENDIF
   
   IF(T_ASSIM) THEN
     DO K=1, KBM1
       DO I=1, MGL
         IDUMMY = IDUMMY + 1
         STFCT(IDUMMY) = T1TMP(I,K)
         STINDX(IDUMMY,1)=1
         STINDX(IDUMMY,2)=I
         STINDX(IDUMMY,3)=K 
       ENDDO
     ENDDO
   ENDIF
   
   IF(S_ASSIM) THEN
     DO K=1, KBM1
       DO I=1, MGL
         IDUMMY = IDUMMY + 1
         STFCT(IDUMMY) = S1TMP(I,K)
         STINDX(IDUMMY,1)=1
         STINDX(IDUMMY,2)=I
         STINDX(IDUMMY,3)=K 
       ENDDO
     ENDDO
   ENDIF
END IF !(msr)
   DEALLOCATE(UTMP,VTMP,S1TMP,T1TMP,ELTMP,UATMP,VATMP,TMP1,TMP2,TMP3,TMP4,TMP5,TMP6,TMP8)
#  if defined (WATER_QUALITY)
   DEALLOCATE(TMP7)
#  endif
   RETURN
   END SUBROUTINE GR2ST 
!#####################################
   SUBROUTINE GR2ST_TRUE_FIELD(GNAM)
   
   USE LIMS
   USE ALL_VARS
#  if defined (WATER_QUALITY)
   USE MOD_WQM
#  endif
!#  if defined (EQUI_TIDE)
!   USE MOD_EQUITIDE
!#  endif
!#  if defined (ATMO_TIDE)
!   USE MOD_ATMOTIDE
!#  endif
   IMPLICIT NONE
   CHARACTER(len=120) :: GNAM
!   CHARACTER(len=8) :: chint
   INTEGER I,K,N1,IDUMMY ,HO
   INTEGER INF
   INTEGER ITMP
   REAL(DP) :: TMP_TIME
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
   REAL(SP), ALLOCATABLE :: TMP8(:,:)
#  if defined (WATER_QUALITY)
   REAL(SP), ALLOCATABLE :: TMP7(:,:,:)
#  endif

   ALLOCATE(UTMP(0:NGL,1:KB))      ; UTMP  = 0.0_SP 
   ALLOCATE(VTMP(0:NGL,1:KB))      ; VTMP  = 0.0_SP
   ALLOCATE(S1TMP(0:MGL,1:KB))     ; S1TMP = 0.0_SP
   ALLOCATE(T1TMP(0:MGL,1:KB))     ; T1TMP = 0.0_SP
   ALLOCATE(ELTMP(0:MGL))          ; ELTMP = 0.0_SP
   ALLOCATE(UATMP(0:NGL))          ; UATMP = 0.0_SP 
   ALLOCATE(VATMP(0:NGL))          ; VATMP = 0.0_SP
   ALLOCATE(TMP1(0:NGL,1:KB))      ; TMP1  = 0.0_SP
   ALLOCATE(TMP2(0:MGL,1:KB))      ; TMP2  = 0.0_SP
   ALLOCATE(TMP8(0:MGL,1:KB))      ; TMP8  = 0.0_SP
   ALLOCATE(TMP3(0:NGL,1:KB))      ; TMP3  = 0.0_SP
   ALLOCATE(TMP4(0:NGL))           ; TMP4  = 0.0_SP
   ALLOCATE(TMP5(0:NGL))           ; TMP5  = 0.0_SP
   ALLOCATE(TMP6(0:MGL))           ; TMP6  = 0.0_SP
#  if defined (WATER_QUALITY)
   ALLOCATE(TMP7(1:MGL,1:KB,1:NB)) ; TMP7  = 0.0_SP 
#  endif
   TMP_TIME=DAYS(INTTIME)
!   CALL NCD_FIND_READ_TIME_ENKF(trim(GNAM),TMP_TIME,HO)
!   call NCD_READ_ENKF(trim(GNAM),UTMP,VTMP,T1TMP,S1TMP,ELTMP,HO)
    CALL ENKF_READ_FIELD(trim(GNAM),ELTMP,UTMP,VTMP,T1TMP,S1TMP)
if (msr) then
   IDUMMY = 0
   IF(EL_ASSIM) THEN
     DO I=1, MGL
       IDUMMY = IDUMMY + 1
       STFCT(IDUMMY) = ELTMP(I)
       STINDX(IDUMMY,1)=1
       STINDX(IDUMMY,2)=I
       STINDX(IDUMMY,3)=1 
     ENDDO
   ENDIF
   
   IF(UV_ASSIM) THEN
     DO K=1, KBM1
       DO I=1, NGL
         IDUMMY = IDUMMY + 1
         STFCT(IDUMMY) = UTMP(I,k)
         STINDX(IDUMMY,1)=2
         STINDX(IDUMMY,2)=I
         STINDX(IDUMMY,3)=K 
       ENDDO
     ENDDO
     DO K=1, KBM1
       DO I=1, NGL
         IDUMMY = IDUMMY + 1
         STFCT(IDUMMY) = VTMP(I,k)
         STINDX(IDUMMY,1)=2
         STINDX(IDUMMY,2)=I
         STINDX(IDUMMY,3)=K
       ENDDO
     ENDDO
   ENDIF
   
   IF(T_ASSIM) THEN
     DO K=1, KBM1
       DO I=1, MGL
         IDUMMY = IDUMMY + 1
         STFCT(IDUMMY) = T1TMP(I,K)
         STINDX(IDUMMY,1)=1
         STINDX(IDUMMY,2)=I
         STINDX(IDUMMY,3)=K 
       ENDDO
     ENDDO
   ENDIF
   
   IF(S_ASSIM) THEN
     DO K=1, KBM1
       DO I=1, MGL
         IDUMMY = IDUMMY + 1
         STFCT(IDUMMY) = S1TMP(I,K)
         STINDX(IDUMMY,1)=1
         STINDX(IDUMMY,2)=I
         STINDX(IDUMMY,3)=K 
       ENDDO
     ENDDO
   ENDIF
end if !(msr)
   DEALLOCATE(UTMP,VTMP,S1TMP,T1TMP,ELTMP,UATMP,VATMP,TMP1,TMP2,TMP3,TMP4,TMP5,TMP6,TMP8)
#  if defined (WATER_QUALITY)
   DEALLOCATE(TMP7)
#  endif
   RETURN
   END SUBROUTINE GR2ST_TRUE_FIELD
!#####################################
   SUBROUTINE GETOBSLOC1
   
   USE LIMS
   USE CONTROL
   IMPLICIT NONE

   INTEGER ::  NUM 
   INTEGER  SWITCH
   SAVE     SWITCH
   INTEGER ::  J,K,PNT
   INTEGER ::  IDUMMY
   INTEGER ::  TMP
   CHARACTER(LEN=80) FILENAME
   CHARACTER(LEN=24) HEADINFO
   INTEGER STLTMP(ENKF_NOBSMAX)
   INTEGER LAY(ENKF_NOBSMAX)

   NUM     = 0 
   IDUMMY  = 0
   NLOC    = 0
   PNT     = 0
   STLOC   = 0
   STLTMP  = 0 
   LAY     = 0
!   STTR    = 0.0_DP
   WKTMP   = 0.0_DP

   SWITCH  = 0
   
   FILENAME = TRIM(INPUT_DIR)//"/"//trim(casename)//"_assim_enkf.nml"
    
   OPEN(INOOB,FILE=TRIM(FILENAME),FORM='FORMATTED')
   
 100 READ(INOOB,*,END=200) HEADINFO
     IF(SWITCH/=1) THEN
       IF(HEADINFO=='!===READ') THEN
         SWITCH = 1
         GOTO 100
       ELSE
         GOTO 100
       ENDIF
     ENDIF 

     IF(TRIM(HEADINFO)=='!EL') THEN
       IF(EL_OBS) THEN
         READ(INOOB,*) NUM
	 NLOC = NLOC + NUM
!         IF(NLOC>ENKF_NOBSMAX) THEN
!           WRITE(IPT,*) 'not enough storage for observations:', 'Nloc=', Nloc, 'Nobsmax=', ENKF_NOBSMAX
!           CALL PSTOP
!         ENDIF
	 READ(INOOB,*)  (STLOC(K), K=1,NLOC)	
	 
	   DO K=1, NLOC   
              WKTMP(k) = OBSERR_EL      ! FIND THE NAME OF VARIABLE
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
!         IF(NLOC+NUM>ENKF_NOBSMAX) THEN
!           WRITE(IPT,*) 'not enough storage for observations:', 'Nloc=', Nloc+num, 'Nobsmax=', ENKF_NOBSMAX
!           CALL PSTOP
!         ENDIF
	 READ(INOOB,*)  (STLOC(K), K=NLOC-NUM+1,NLOC)
	 READ(INOOB,*)  (LAY(K),   K=NLOC-NUM+1,NLOC)
         DO K=NLOC-NUM+1, NLOC
	   STLOC(K)=STLOC(K)+IDUMMY+NGL*(LAY(K)-1)
	   WKTMP(K) = OBSERR_UV
	 ENDDO   
	 
	 NLOC = NLOC + NUM
	 DO K=NLOC-NUM+1, NLOC
	   STLOC(K)=STLOC(K-NUM)+IDUMMY+NGL*KBM1
	   WKTMP(K) = OBSERR_UV
	 ENDDO
	 
       ENDIF
       IF(UV_ASSIM) THEN  
          IDUMMY = IDUMMY + 2*NGL*KBM1
       ENDIF
     ENDIF
     
     IF(TRIM(HEADINFO)=='!T') THEN
       IF(T_OBS) THEN
         READ(INOob,*) NUM
	 NLOC = NLOC + NUM
!         IF(NLOC>ENKF_NOBSMAX) THEN
!           WRITE(IPT,*) 'not enough storage for observations:', 'Nloc=', Nloc, 'Nobsmax=', ENKF_NOBSMAX
!           CALL PSTOP
!         ENDIF
	 READ(INOOB,*)  (STLOC(K), K=NLOC-NUM+1,NLOC)       
         READ(INOOB,*)  (LAY(K),   K=NLOC-NUM+1,NLOC)
         DO K=NLOC-NUM+1, NLOC
	   STLOC(K)=STLOC(K)+IDUMMY+MGL*(LAY(K)-1)
	   WKTMP(K) = OBSERR_T
	 ENDDO   
        
       ENDIF   
       IF(T_ASSIM) THEN
         IDUMMY = IDUMMY + MGL*1
       ENDIF
     ENDIF

     IF(TRIM(HEADINFO)=='!S') THEN
       IF(S_OBS) THEN
         READ(INOOB,*) NUM
	 NLOC = NLOC + NUM
!         IF(NLOC>ENKF_NOBSMAX) THEN
!           WRITE(IPT,*) 'not enough storage for observations:', 'Nloc=', Nloc, 'Nobsmax=', ENKF_NOBSMAX
!           CALL PSTOP
!         ENDIF
	 READ(INOOB,*)  (STLOC(K),K=NLOC-NUM+1,NLOC)        
         READ(INOOB,*)  (LAY(K),  K=NLOC-NUM+1,NLOC)
         DO K=NLOC-NUM+1, NLOC
	   STLOC(K)=STLOC(K)+IDUMMY+MGL*(LAY(K)-1)
	   WKTMP(K) = OBSERR_S
	 ENDDO   
         
       ENDIF
       IF(S_ASSIM) THEN
         IDUMMY = IDUMMY + MGL*KBM1 
       ENDIF 
     ENDIF

     GOTO 100
 200 CONTINUE
     REWIND(INOOB)
     
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
   END SUBROUTINE GETOBSLOC1
!####################################
!#####################################
   SUBROUTINE GETNLOC
   
   USE LIMS
   USE CONTROL
   IMPLICIT NONE

   INTEGER ::  NUM 
   INTEGER  SWITCH
   SAVE     SWITCH
   INTEGER ::  J,K,PNT
   INTEGER ::  IDUMMY
   INTEGER ::  TMP
   CHARACTER(LEN=80) FILENAME
   CHARACTER(LEN=24) HEADINFO
   INTEGER STLTMP(ENKF_NOBSMAX)
   INTEGER LAY(ENKF_NOBSMAX)

   NUM     = 0 
   IDUMMY  = 0
   NLOC    = 0
   PNT     = 0
   STLTMP  = 0 
   LAY     = 0
!   STTR    = 0.0_DP

   SWITCH  = 0
   
   FILENAME = TRIM(INPUT_DIR)//"/"//trim(casename)//"_assim_enkf.nml"
    
   OPEN(INOOB,FILE=TRIM(FILENAME),FORM='FORMATTED')
   REWIND(INOOB)
 100 READ(INOOB,*,END=200) HEADINFO
     IF(SWITCH/=1) THEN
       IF(HEADINFO=='!===READ') THEN
         SWITCH = 1
         GOTO 100
       ELSE
         GOTO 100
       ENDIF
     ENDIF 

     IF(TRIM(HEADINFO)=='!EL') THEN
       IF(EL_OBS) THEN
         READ(INOOB,*) NUM
	 NLOC = NLOC + NUM
	 
       ENDIF

     ENDIF
     
     IF(TRIM(HEADINFO)=='!UV') THEN
       IF(UV_OBS) THEN
         READ(INOOB,*) NUM
	 NLOC = NLOC + 2*NUM 
       ENDIF

     ENDIF
     
     IF(TRIM(HEADINFO)=='!T') THEN
       IF(T_OBS) THEN
         READ(INOob,*) NUM
	 NLOC = NLOC + NUM
       ENDIF   
     ENDIF

     IF(TRIM(HEADINFO)=='!S') THEN
       IF(S_OBS) THEN
         READ(INOOB,*) NUM
	 NLOC = NLOC + NUM
         
       ENDIF
 
     ENDIF

     GOTO 100
 200 CONTINUE
     REWIND(INOOB)
     
    IF(NLOC==0) THEN 
       WRITE(IPT,*) "!WARNING, NOT OBSERVATION DATA FOUND!"
    ENDIF

    RETURN
   END SUBROUTINE GETNLOC



!##################################

   FUNCTION DISTST(STLOC1,STLOC2)

   USE LIMS
   USE CONTROL
   USE ALL_VARS
   IMPLICIT NONE

   INTEGER IDUMMY, IDUMMY1
   INTEGER RNG(2)
   INTEGER STLOC1,STLOC2
   INTEGER STB,STE,LOC1,LOC2,K
   REAL(DP) BIG_DIST
   REAL*8   DISTST
   real     x1,y1,x2,y2
   INTEGER :: num1,num2
   num1=stindx(stloc1,2)
   num2=stindx(stloc2,2)
   if (stindx(stloc1,1)==1 ) then
       x1=xg(num1)
       y1=yg(num1)
   else
       x1=xcg(num1)
       y1=ycg(num1)
   end if
   if (stindx(stloc2,1)==1 ) then
       x2=xg(num2)
       y2=yg(num2)
   else
       x2=xcg(num2)
       y2=ycg(num2)
   end if 
   distst=sqrt((x1-x2)**2+(y1-y2)**2)
   RETURN
   END FUNCTION DISTST


   SUBROUTINE gaussj(a,n,np,b,m,mp)

   USE MOD_PREC
   IMPLICIT NONE
     
   INTEGER m,mp,n,np,NMAX
   REAL(DP) a(np,np),b(np,mp)
   PARAMETER (NMAX=1000)
   INTEGER i,icol,irow,j,k,l,ll,indxc(NMAX),indxr(NMAX),ipiv(NMAX)
   REAL(DP) big,dum,pivinv

   do 11 j=1,n
      ipiv(j)=0
11 continue
   do 22 i=1,n
      big=0.
      do 13 j=1,n
        if(ipiv(j).ne.1)then
          do 12 k=1,n
            if (ipiv(k).eq.0) then
              if (abs(a(j,k)).ge.big)then
                big=abs(a(j,k))
                irow=j
                icol=k
              endif

            else if (ipiv(k).gt.1) then
              pause 'singular matrix in gaussj'
            endif
12       continue
       endif
13   continue
     
     ipiv(icol)=ipiv(icol)+1
     if (irow.ne.icol) then
       do 14 l=1,n
         dum=a(irow,l)
         a(irow,l)=a(icol,l)
         a(icol,l)=dum
14     continue
       do 15 l=1,m
         dum=b(irow,l)
         b(irow,l)=b(icol,l)
         b(icol,l)=dum
15     continue
       endif
       indxr(i)=irow
       indxc(i)=icol
       if (a(icol,icol).eq.0.) pause 'singular matrix in gaussj'

       pivinv=1./a(icol,icol)
       a(icol,icol)=1.
       do 16 l=1,n
         a(icol,l)=a(icol,l)*pivinv
16     continue
       do 17 l=1,m
         b(icol,l)=b(icol,l)*pivinv
17     continue
       do 21 ll=1,n
         if(ll.ne.icol)then
           dum=a(ll,icol)
           a(ll,icol)=0.
           do 18 l=1,n
             a(ll,l)=a(ll,l)-a(icol,l)*dum
18         continue
           do 19 l=1,m
             b(ll,l)=b(ll,l)-b(icol,l)*dum
19         continue
         endif
21     continue
22  continue
    do 24 l=n,1,-1
       if(indxr(l).ne.indxc(l))then

         do 23 k=1,n
           dum=a(k,indxr(l))
           a(k,indxr(l))=a(k,indxc(l))
           a(k,indxc(l))=dum
23       continue
       endif
24  continue
    
    RETURN
    END SUBROUTINE GAUSSJ

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
!-----------------------
   SUBROUTINE FCT2ANL_NC(GLNAME)

   USE LIMS
   USE ALL_VARS
#  if defined (WATER_QUALITY)
   USE MOD_WQM
#  endif
!#  if defined (EQUI_TIDE)
!   USE MOD_EQUITIDE
!#  endif
!#  if defined (ATMO_TIDE)
!   USE MOD_ATMOTIDE
!#  endif
   IMPLICIT NONE

   INTEGER I,K,J,N1,ITMP
   INTEGER IDUMMY
   CHARACTER(LEN=120)    :: GLNAME, FLNAME

     real(sp), allocatable ::  enkfel(:)   
     real(sp), allocatable ::  enkfu(:,:), enkfv(:,:)
     real(sp), allocatable ::  enkft(:,:), enkfs(:,:)
        

   allocate(enkfel(mgl))        ; enkfel = 0.0_sp
   allocate(enkfu(ngl,kbm1))    ; enkfu  = 0.0_sp
   allocate(enkfv(ngl,kbm1))    ; enkfv  = 0.0_sp
   allocate(enkft(mgl,kbm1))    ; enkft  = 0.0_sp
   allocate(enkfs(mgl,kbm1))    ; enkfs  = 0.0_sp
   
!---------------------------------------
   IDUMMY = 0
   IF(EL_ASSIM) THEN
      DO I=1, MGL
         IDUMMY = IDUMMY + 1
         ENKFEL(I) = STFCT(IDUMMY)
      ENDDO
      call NCD_WRITE_EL(GLNAME,MGL,1,ENKFEL) 
   ENDIF
   IF(UV_ASSIM) THEN
     print *, 'WRITE ANALYSIS%%%%%%%%%%%%%%%%%%%%%%%'
      DO K=1, KBM1
        DO I=1, NGL
          IDUMMY = IDUMMY + 1
          ENKFU(I,K) = STFCT(IDUMMY)  
        ENDDO
      ENDDO  

      DO K=1, KBM1
        DO I=1, NGL
          IDUMMY = IDUMMY + 1
          ENKFV(I,K) = STFCT(IDUMMY)  
        ENDDO
      ENDDO
    call NCD_WRITE_U(GLNAME,NGL,KBM1,ENKFU)
    call NCD_WRITE_V(GLNAME,NGL,KBM1,ENKFV)
   ENDIF

   IF(T_ASSIM) THEN
      DO K=1, KBM1
        DO I=1, MGL
          IDUMMY = IDUMMY + 1
          ENKFT(I,K) = STFCT(IDUMMY)  
        ENDDO
      ENDDO 
      call NCD_WRITE_T(GLNAME,MGL,KBM1,ENKFT)
   ENDIF

   IF(S_ASSIM) THEN
      DO K=1, KBM1
        DO I=1, MGL
          IDUMMY = IDUMMY + 1
          ENKFS(I,K) = STFCT(IDUMMY)  
        ENDDO
      ENDDO  
      call NCD_WRITE_S(GLNAME,MGL,KBM1,ENKFS)
   ENDIF

!----------------------------   

   RETURN
   END SUBROUTINE FCT2ANL_NC
!----------------------------\\\

!-----------------------
   SUBROUTINE FCT2ANL(GLNAME)

   USE LIMS
   USE ALL_VARS
#  if defined (WATER_QUALITY)
   USE MOD_WQM
#  endif
!#  if defined (EQUI_TIDE)
!   USE MOD_EQUITIDE
!#  endif
!#  if defined (ATMO_TIDE)
!   USE MOD_ATMOTIDE
!#  endif
   IMPLICIT NONE

   INTEGER I,K,J,N1,ITMP
   INTEGER IDUMMY
   CHARACTER(LEN=120)    :: GLNAME, FLNAME

     real(sp), allocatable ::  enkfel(:)   
     real(sp), allocatable ::  enkfu(:,:), enkfv(:,:)
     real(sp), allocatable ::  enkft(:,:), enkfs(:,:)
        

   allocate(enkfel(mgl))        ; enkfel = 0.0_sp
   allocate(enkfu(ngl,kbm1))    ; enkfu  = 0.0_sp
   allocate(enkfv(ngl,kbm1))    ; enkfv  = 0.0_sp
   allocate(enkft(mgl,kbm1))    ; enkft  = 0.0_sp
   allocate(enkfs(mgl,kbm1))    ; enkfs  = 0.0_sp
   
!---------------------------------------
   IDUMMY = 0
   IF(EL_ASSIM) THEN
     OPEN(91,FILE=TRIM(GLNAME)//'EL', FORM='UNFORMATTED')
     REWIND(91)
     IF (ICYC==0) THEN
        WRITE(91) I_INITIAL
     ELSE
     WRITE(91) IEND
     END IF
      DO I=1, MGL
         IDUMMY = IDUMMY + 1
         ENKFEL(I) = STFCT(IDUMMY)
      ENDDO
   DO I=1,MGL
     WRITE(91) ENKFEL(I) 
   END DO
   CLOSE(91)
   ENDIF
   IF(UV_ASSIM) THEN
     OPEN(91,FILE=TRIM(GLNAME)//'UV', FORM='UNFORMATTED')
     REWIND(91)
     print *, ICYC,ENKF_START/DELTA_ASS,'ICYC,ENKF_START/DELTA_ASS'
     print *, I_INITIAL,IEND
     IF (ICYC==0) THEN
        WRITE(91) I_INITIAL
!        print *, 'I_INITIAL',I_INITIAL
     ELSE
     WRITE(91) IEND
     END IF
!     print *, 'WRITE ANALYSIS%%%%%%%%%%%%%%%%%%%%%%%'
      DO K=1, KBM1
        DO I=1, NGL
          IDUMMY = IDUMMY + 1
          ENKFU(I,K) = STFCT(IDUMMY)  
        ENDDO
      ENDDO  

      DO K=1, KBM1
        DO I=1, NGL
          IDUMMY = IDUMMY + 1
          ENKFV(I,K) = STFCT(IDUMMY)  
        ENDDO
      ENDDO
    DO I=1,NGL
     WRITE(91) (ENKFU(I,K),ENKFV(I,K),K=1,KBM1)   !!U,V
   END DO
   CLOSE(91)
   ENDIF

   IF(T_ASSIM) THEN
     OPEN(91,FILE=TRIM(GLNAME)//'T', FORM='UNFORMATTED')
     REWIND(91)
     IF (ICYC==0) THEN
        WRITE(91) I_INITIAL
     ELSE
     WRITE(91) IEND
     END IF
      DO K=1, KBM1
        DO I=1, MGL
          IDUMMY = IDUMMY + 1
          ENKFT(I,K) = STFCT(IDUMMY)  
        ENDDO
      ENDDO 
    DO I=1,MGL
     WRITE(91) (ENKFT(I,K),K=1,KBM1) 
   END DO
   CLOSE(91) 
   ENDIF

   IF(S_ASSIM) THEN
     OPEN(91,FILE=TRIM(GLNAME)//'S', FORM='UNFORMATTED')
     REWIND(91)
     IF (ICYC==0) THEN
        WRITE(91) I_INITIAL
     ELSE
     WRITE(91) IEND
     END IF

      DO K=1, KBM1
        DO I=1, MGL
          IDUMMY = IDUMMY + 1
          ENKFS(I,K) = STFCT(IDUMMY)  
        ENDDO
      ENDDO  
    DO I=1,MGL
     WRITE(91) (ENKFS(I,K),K=1,KBM1)
   END DO
   CLOSE(91)
   ENDIF

!----------------------------   

   RETURN
   END SUBROUTINE FCT2ANL

   SUBROUTINE SET_INI_ENKF
     
   USE LIMS
   IMPLICIT NONE

   INTEGER I,J,K
   INTEGER IERR
   REAL(DP) SUM9
   REAL(DP),ALLOCATABLE   ::  STREF(:)
   REAL(DP),ALLOCATABLE   ::  RPETMP(:,:)
   REAL(DP),ALLOCATABLE   ::  RPATMP(:)   
   CHARACTER(LEN=120)     ::  FLNAME, FLNAME2
   CHARACTER(LEN=4)       ::  IFIL
   
   STDIM = 0
   print *, 'EL_ASSIM' , EL_ASSIM
      print *, 'UV_ASSIM' , UV_ASSIM
   print *, 'T_ASSIM' , T_ASSIM
   print *, 'S_ASSIM' , S_ASSIM 
   IF(EL_ASSIM) STDIM = STDIM + MGL
   IF(UV_ASSIM) STDIM = STDIM + 2*NGL*KBM1
   IF(T_ASSIM)  STDIM = STDIM + MGL*KBM1
   IF(S_ASSIM)  STDIM = STDIM + MGL*KBM1

   RETURN 
   END SUBROUTINE SET_INI_ENKF

!--------------------------------------------------      
      SUBROUTINE GAUSSJ_2(A,N,NP,B,M,MP)
      IMPLICIT NONE
      INTEGER M,MP,N,NP,NMAX
      REAL*8 A(NP,NP),B(NP,MP)
      PARAMETER (NMAX=50)
      INTEGER I,ICOL,IROW,J,K,L,LL,INDXC(NMAX),INDXR(NMAX),&
             IPIV(NMAX)
      REAL*8 BIG,DUM,PIVINV
      DO J=1,N
         IPIV(J)=0
      ENDDO
      DO 22 I=1,N
         BIG=0.
         DO 13 J=1,N
            IF(IPIV(J).NE.1)THEN
              DO 12 K=1,N
                 IF (IPIV(K).EQ.0) THEN
                   IF (ABS(A(J,K)).GE.BIG)THEN
                     BIG=ABS(A(J,K))
                     IROW=J
                     ICOL=K
                   ENDIF

                 ELSE IF (IPIV(K).GT.1) THEN
                   PAUSE 'SINGULAR MATRIX IN GAUSSJ'
                 ENDIF
12            CONTINUE
            ENDIF
13       CONTINUE
         IPIV(ICOL)=IPIV(ICOL)+1
         IF (IROW.NE.ICOL) THEN
           DO 14 L=1,N
              DUM=A(IROW,L)
              A(IROW,L)=A(ICOL,L)
              A(ICOL,L)=DUM
14         CONTINUE
           DO 15 L=1,M
              DUM=B(IROW,L)
              B(IROW,L)=B(ICOL,L)
              B(ICOL,L)=DUM
15         CONTINUE
         ENDIF
         INDXR(I)=IROW
         INDXC(I)=ICOL
         IF (A(ICOL,ICOL).EQ.0.) PAUSE 'SINGULAR MATRIX IN GAUSSJ'

         PIVINV=1./A(ICOL,ICOL)
         A(ICOL,ICOL)=1.
         DO 16 L=1,N
            A(ICOL,L)=A(ICOL,L)*PIVINV
16       CONTINUE
         DO 17 L=1,M
            B(ICOL,L)=B(ICOL,L)*PIVINV
17       CONTINUE
         DO 21 LL=1,N
            IF(LL.NE.ICOL)THEN
              DUM=A(LL,ICOL)
              A(LL,ICOL)=0.
              DO 18 L=1,N
                 A(LL,L)=A(LL,L)-A(ICOL,L)*DUM
18            CONTINUE
              DO 19 L=1,M
                 B(LL,L)=B(LL,L)-B(ICOL,L)*DUM
19            CONTINUE
            ENDIF
21       CONTINUE
22    CONTINUE

      DO 24 L=N,1,-1
         IF(INDXR(L).NE.INDXC(L))THEN

           DO 23 K=1,N
              DUM=A(K,INDXR(L))
              A(K,INDXR(L))=A(K,INDXC(L))
              A(K,INDXC(L))=DUM
23         CONTINUE
         ENDIF
24    CONTINUE

      RETURN
      END SUBROUTINE GAUSSJ_2
      
   SUBROUTINE WD_READ_SINGLE(FNAME)

!------------------------------------------------------------------------------|
   USE ALL_VARS
   IMPLICIT NONE
   INTEGER, ALLOCATABLE,DIMENSION(:) :: NTEMP1,NTEMP2
   INTEGER I,IINT_TMP
   CHARACTER(LEN=120) :: FNAME
   LOGICAL :: FEXIST
!==============================================================================|


!
!--Make Sure Wet/Dry Flag Data Exists------------------------------------------!
!
   INQUIRE(FILE=TRIM(FNAME),EXIST=FEXIST)
   IF(MSR .AND. .NOT.FEXIST)THEN
     WRITE(IPT,*)'WET/DRY RESTART FILE: ',FNAME,' DOES NOT EXIST'
     WRITE(IPT,*)'HALTING.....'
     CALL PSTOP
   END IF

!
!--Ensure File Header is Consistent with Main Flow Variable Restart File-------!
!
   OPEN(1,FILE=FNAME,FORM='FORMATTED',STATUS='UNKNOWN')
   REWIND(1)
   READ(1,*) IINT_TMP
   READ(1,*)
   IF(IINT_TMP /= IINT .AND. MSR)THEN
     WRITE(IPT,*)'IINT IN ',FNAME,' NOT EQUAL TO IINT'
     WRITE(IPT,*)'FROM MAIN RESTART FILE',IINT,IINT_TMP
     CALL PSTOP
   END IF

!
!--Read Variables--------------------------------------------------------------!
!
   ALLOCATE(NTEMP1(NGL),NTEMP2(MGL))
   READ(1,*) (NTEMP1(I), I=1,NGL)
   READ(1,*) (NTEMP2(I), I=1,MGL)

!
!--Transfer Variables to Local Domains-----------------------------------------!

!
!--Transfer Variables to Local Domains-----------------------------------------!
!

     ENKFWETC = NTEMP1(1:NGL)*ENKFWETC
     ENKFWETN = NTEMP2(1:MGL)*ENKFWETN

   DEALLOCATE(NTEMP1,NTEMP2)
   CLOSE(1)
!

   RETURN
   END SUBROUTINE WD_READ_SINGLE

!==============================================================================|

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine etkf(Xb,R,Yb,innov,ndim,nrens,nrobs)
   implicit none
   integer, intent(in) :: ndim             ! dimension of model state
   integer, intent(in) :: nrens            ! number of ensemble members
   integer, intent(in) :: nrobs            ! number of observations
   real(dp), intent(inout) :: Xb(ndim,nrens)    ! ensemble matrix
   real(dp), intent(inout)    :: R(nrobs,nrobs)   ! matrix holding R (only used if mode=?1 or ?2)
   real(dp), intent(inout)    :: Yb(nrobs,nrens)   ! matrix holding HA`
   real(dp), intent(inout)    :: innov(nrobs)     ! vector holding d-H*mean(A)


   real(dp)  :: Ninnov(nrobs)     ! normalized vector holding R^(-1/2)(d-H*mean(A))

   real(dp) :: NYb(nrobs,nrens)   ! matrix holding HA` Normalized by R^(-1/2)
   real(dp) :: NYbT(nrobs,nrens)   ! NYb^T
   real(dp) :: ZHHZ(nrens,nrens)   ! NYb^T*NYb   
   real(dp) :: E(nrobs,nrens) ! E is the NYb*W*eig as in Bishop's Paper.
   real(dp) :: E_t(nrens,nrobs) ! E^T
   real(dp) :: sig_var(ndim) ! signal variance , the trace of signal covariance.
   integer :: i,j
!--------------------temporay array--------------------------
   real(dp),allocatable :: c(:,:),invr(:,:),w(:,:),w_t(:,:),T(:,:),eig(:),eig_1(:),pa(:,:),xbar(:),X1(:,:),E_t_Ninnov(:),W_E_t_Ninnov(:),Xb_w(:,:)
!   real(dp) :: inflate=0.04
!   integer :: info
!------------------------------------------------------------
! step 0: check if there are observation in the local area and calculate the mean of Xb.
   print *, 'step 0'
   if (nrobs==0) return
!cacluate Xb=Xb-mean(Xb)
   allocate(xbar(ndim))
   xbar=sum(xb,2)
   xbar=xbar/dble(nrens)
  do i=1,nrens
     Xb(:,i)=Xb(:,i)-xbar
  end do
! step 1 : calcuate the NYb=R^(-1/2)*Yb , that is the H(hat)Z in Bishop's paper
    do i=1,nrens
       do j=1,nrobs
          NYb(j,i)=Yb(j,i)/sqrt(R(j,j))
          NYbT(i,j)=NYb(j,i)
       end do
    end do
   print *,'step 1'
! step 2:  calculate the ZHHZ= NYb^T*NYb
   print *, 'step2'
   ZHHZ=matmul(NYbT,NYb)
! step 3: calculate the eigen value and eigen vector of ZHHZ
   print *, 'step 3'
   allocate(W(nrens,nrens)) ! W is eigen vector of ZHHZ, as the C in Bishop's paper.
   allocate(eig(nrens))      ! eig is the eigen value of ZHHZ.
   call eigC(ZHHZ,nrens,W,eig)
   allocate(W_t(nrens,nrens))
   do i=1,nrens
      do j=1,nrens
         W_t(i,j)=W(j,i)
      end do
   end do
   eig_1=eig
   do i=1, nrens
      if (eig_1(i)<=0) then
          eig_1(i)=1        ! set all zero eigenvalues to 1
      end if
   end do
! step 4 calculate the signal variance
allocate(Xb_w(ndim,nrens)) 
Xb_w=matmul(Xb,W)
do i=1,ndim
  do j=1,nrens
      Xb_w(i,j)=Xb_w(i,j)*sqrt(eig(j))/sqrt(eig(j)+1) ! Xb*W*eig^(1/2)*(eig+I)^(-1/2)
  end do
end do 
sig_var=0
do i=1,ndim
   do j=1,nrens
      sig_var(i)=sig_var(i)+Xb_w(i,j)*Xb_w(i,j)
   end do
end do

! step 5.1. calculate E=NYb*W*sqrt(eig)^(-1/2)
   print *, 'step 5.1'
          E=matmul(NYb,W)
          do i=1,nrobs
             do j=1,nrens
                 E(i,j)=E(i,j)/eig_1(j)
                 E_t(j,i)=E(i,j)
             end do
          end do
! step 5.2: calculate the normalized innovation
   print *, 'step 5.2'
   do i=1,nrobs
      Ninnov(i)=innov(i)/sqrt(R(i,i))
   end do
! step 6:  calculate the Kalman Gain* innovation
! step 6.1 : calculate the E^T*normalized innovation
   print *, 'step 6.1'
allocate(E_t_Ninnov(nrens)) 
   E_t_Ninnov=matmul(E_t,Ninnov)
! step 6.2 : calculate the (eig+I)^-1*E_t_Ninnov
   print *, 'step 6.2'
   do i=1,nrens
     E_t_Ninnov(i)=E_t_Ninnov(i)/(eig(i)+1)
   end do
! step 6.3: calculate the sqrt(eig)*(eig+I)^-1*E_t_Ninnov
   print *, 'step 6.3'
   do i=1,nrens
      E_t_Ninnov(i)=E_t_Ninnov(i)*sqrt(eig(i))
   end do
! step 6.4 : calculate the W*sqrt(eig)*(eig+I)^-1*E_t_Ninnov
     print *, 'step 6.4'
     allocate(W_E_t_Ninnov(nrens))
     W_E_t_Ninnov=matmul(W,E_t_Ninnov)
! step 6.5 : calculate the  analysis Xbar=A'* W*sqrt(eig)*(eig+I)^-1*E_t_Ninnov+Xbar
     print *, 'step 6.5'
     Xbar=Xbar+matmul(Xb,W_E_t_Ninnov)
! step 7. calculate the transformation matrix T
    print *, 'step 7'
    do i=1,nrens
       do j=1,nrens
         T(i,j)=W(i,j)/sqrt((eig(j)+1))  ! T=W*(eig+I)^-1/2
       end do
    end do
    T=matmul(T,W_t)
    Xb=matmul(Xb,T)
    do i=1,nrens
       Xb(:,i)=Xb(:,i)+Xbar   ! form the analysis ensemble 
    end do

end subroutine




!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine assimilate(Xb,R,Yb,innov,ndim,nrens,nrobs)
   implicit none
   integer, intent(in) :: ndim             ! dimension of model state
   integer, intent(in) :: nrens            ! number of ensemble members
   integer, intent(in) :: nrobs            ! number of observations
   real(dp), intent(inout) :: Xb(ndim,nrens)    ! ensemble matrix
   real(dp), intent(inout)    :: R(nrobs,nrobs)   ! matrix holding R (only used if mode=?1 or ?2)
   real(dp), intent(inout)    :: Yb(nrobs,nrens)   ! matrix holding HA`
   real(dp), intent(inout)    :: innov(nrobs)     ! vector holding d-H*mean(A)
   integer :: i,j
!--------------------temporay array--------------------------
   real(dp),allocatable :: c(:,:),invr(:,:),w(:,:),eig(:),pa(:,:),xbar(:),X1(:,:),cyb(:),ceb(:,:)
   real(dp) :: inflate=0.04
   integer :: info
!------------------------------------------------------------
! step 0: check if there are observation in the local area.
   print *, 'step 0'
   if (nrobs==0) return
! step 1 : Yb is the Yb in Brian's paper
   print *,'step 1'
! step 2:  cacluate Xb=Xb-mean(Xb)
   print *, 'step2'
   allocate(xbar(ndim))
   xbar=sum(xb,2)
   xbar=xbar/dble(nrens)
  do i=1,nrens
     Xb(:,i)=Xb(:,i)-xbar
  end do
! step 3: skip since we have already processed in the localized domain.  
   print *, 'step 3'
! step 4: C=Yb^{t}*R^{-1}
   print *, 'step 4'
!  If R is variable, then it may be more efficient
!  to solve R*transpose(C)=H for C than to compute inverse(R) directly.
   
         print *, mode,'start dgemm'
!-----------------------pengfei

!        Evaluate R^{-1}
!----------------------------------------------------------------------------
!        1.Compute eigenvalue decomposition of R -> W*eig*W`
!         allocate(W(nrobs,nrobs))
!         allocate(eig(nrobs))
!         allocate(invR(nrobs,nrobs))
!         allocate(X1(nrobs,nrobs))
!         print *, 'start eigC'
!         call eigC(R,nrobs,W,eig)
!         print *, 'finish eigc and start eigsign'
!         call eigsign(eig,nrobs,truncation)  !eig has become 1/eig
!         print *, 'eigsign'
!   do i=1,nrobs
!   do j=1,nrobs
!      X1(i,j)=eig(i)*W(j,i) ! eig*W^{t}
!   enddo
!   enddo
!   invR=matmul(W,X1) ! R^{-1}=W*eig^{-1}*W`
!   deallocate(W,X1,eig)
!---------------------------------------------------------------------
!  1. Now I assume the R is the diagonal matrix
!-----------------------------------------------------------------------------------
!  2.C=Yb^{t}*R^{-1}
    allocate(c(nrens,nrobs))
    allocate(cyb(nrens))
!     C=matmul(transpose(Yb),invR)
!   deallocate(invR)
    do j=1,nrobs
       if (R(j,j)==0) then
          print *, 'R(j,j) can not be 0,j is', j
          call pstop
       end if
       c(:,j)=yb(j,:)/R(j,j)
    enddo
    cyb=matmul(c,innov)  ! C*(y0-yb)
!  step 5-7: Pa=[(k-1)*I/pho+C*Yb]^{-1} to form {w^{a(i)}}
   print *, 'step 5- step 7'          
   allocate(ceb(nrens,nrens))
     ceb=matmul(c,Yb) !C*Yb=Yb^{t}*R^{-1}*Yb
    do j=1,nrens
       ceb(j,j)=ceb(j,j)+(nrens-1)/(inflate+1.0) ! cyb^{-1} will be Pa
    enddo
    allocate(W(nrens,nrens))
    allocate(pa(nrens,nrens))
    call compute_w(nrens,ceb,cyb,w,pa,info)
!  step 8 : x^{a(i)}=mean(xb)+Xb*w^{a(i)}
    Xb=matmul(Xb,w)
    do i=1,nrens
       Xb(:,i)=Xb(:,i)+xbar
    end do
    deallocate(W,pa,xbar,ceb,cyb,c)
end subroutine
!---------------------------------------------------------------------------
    subroutine compute_w(nsol,ceb,cyb,w,pa,info)
!  COMPUTE_W - Compute the matrix W that gives the analysis perturbations.
!  Internal work routine to compute steps 5 to 7 of the Physica D paper.
!  The operation count here depends only on the ensemble size.
!  NSOL  : number of ensemble solutions
!  CYB  : C * (YOBS - H(XBAR)); see Brian's writeup.
!  CEB  :=: (k-1)I/rho + C*EB; overwritten by its eigenvectors on successful
!         return.  Contents destroyed on unsuccessful return.
!  W  :=  NSOL x NSOL matrix giving the linear combination of ensemble
!     differences that generate the analysis differences.
!  PA  := estimated analysis covariance matrix
!  INFO := set to zero unless an error occurs.  If an error occurs, then
!     the output variables are undefined.
!
    integer,intent(in)::nsol
    real(dp),intent(inout)::ceb(nsol,nsol)
    real(dp),intent(in)::cyb(nsol)
    real(dp),intent(out)::w(nsol,nsol)  ! the update matrix W
    real(dp),intent(out)::pa(nsol,nsol)  ! estimated covariance mtx Pa~
    integer,intent(out)::info
!
!  Local variables
!
    real(dp)::lambda(nsol)  ! eigenvalues of CEB
    real(dp)::work(nsol,nsol)
    real(dp)::v(nsol)  ! Pa~*[C*(Yobs-Yb)]; step 7 of Physica D paper
    integer::j,effrank,nfound
!-----------------pengfei-------
    real(dp) :: truncation=1.0
!----------end pengfei------
!
!  Get the eigenvector decomposition of CEB and use it to do the
!  inversion required to form Pa~ and to form the matrix square root.
!
         print *, 'start eigC'
         call eigC(ceb,nsol,work,lambda)
         print *, 'finish eigc and start eigsign'
         call eigsign(lambda,nsol,truncation)
!
!  Step 5. At this point we have the eigenvector decomposition of a symmetric
!  positive definite matrix, and ceb holds the orthonormal eigenvectors.
!  Pa~ is the inverse of CEB.
!   
    ceb=work
    do j=1,nsol
       work(:,j)=ceb(:,j)*lambda(j)
    enddo
    pa=matmul(work,transpose(ceb))
!
!  Step 6. Compute sqrt((k-1)*Pa~).
!
    do j=1,nsol
       work(:,j)=ceb(:,j)*sqrt(lambda(j))
    enddo
    w=matmul(work,transpose(ceb))*sqrt(dble(nsol-1))
!
!  Step 7. Complete the calculation of W by forming v=Pa~*C*(Yobs-Yb) and adding
!  it to each column of the matrix square root.
!
    v=matmul(pa,cyb)
    do j=1,nsol
       w(:,j)=w(:,j)+v
    enddo
    info=0
    return
    end subroutine compute_w


!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine analysis(A, R, E, S, D, innov, ndim, nrens, nrobs, verbose, truncation,mode,update_randrot,scaling)
! Computes the analysed ensemble for A using the EnKF or square root schemes.

!   use mod_anafunc !pengfei
!   use m_multa !pengfei
   implicit none
   integer, intent(in) :: ndim             ! dimension of model state
   integer, intent(in) :: nrens            ! number of ensemble members
   integer, intent(in) :: nrobs            ! number of observations
   
   real(dp), intent(inout) :: A(ndim,nrens)    ! ensemble matrix
   real(dp), intent(inout)    :: R(nrobs,nrobs)   ! matrix holding R (only used if mode=?1 or ?2)
   real(dp), intent(inout)    :: D(nrobs,nrens)   ! matrix holding perturbed measurments
   real(dp), intent(inout)    :: E(nrobs,nrens)   ! matrix holding perturbations (only used if mode=?3)
   real(dp), intent(inout)    :: S(nrobs,nrens)   ! matrix holding HA` 
   real(dp), intent(inout)    :: innov(nrobs)     ! vector holding d-H*mean(A)

   logical, intent(in) :: verbose          ! Printing some diagnostic output

   real(dp), intent(in)    :: truncation       ! The ratio of variaince retained in pseudo inversion (0.99)

   integer, intent(in) :: mode             ! first integer means (EnKF=1, SQRT=2)
                                           ! Second integer is pseudo inversion
                                           !  1=eigen value pseudo inversion of SS'+(N-1)R
                                           !  2=SVD subspace pseudo inversion of SS'+(N-1)R
                                           !  3=SVD subspace pseudo inversion of SS'+EE'

   logical, intent(in) :: update_randrot   ! Normally true; false for all but first grid point
                                           ! updates when using local analysis since all grid
                                           ! points need to use the same rotation.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   real(DP) X5(nrens,nrens)
   integer i,nrmin,iblkmax
   logical lreps


   real(dp), allocatable :: eig(:)
   real(dp), allocatable :: W(:,:)
   real(dp), allocatable :: X2(:,:)
   real(dp), allocatable :: X3(:,:)
   real(dp), allocatable :: Reps(:,:)
!------------------pengfei temp----------------------------
!   real(dp) :: A_pert(ndim,nrens)
!   real(dp) :: ph(ndim,nrobs)
   integer :: j
!------------------end pengfei----------------------------------
!-----------------pengfei----------------
real(dp) :: temp1,temp2
!-----------------end pengfei-------------------------
!---------------------------pengfei-------------------------------------------------
   logical,parameter ::  lscale=.false.
   real(dp) :: scaling(nrobs)
  if (lscale) then
     print *,'----------------------------------------------------------'
     print *,'EnKF: Scale matrices'
     do i=1,nrobs
!        scaling(i)=1./sqrt(obs(i)%var)

        S(i,:)=scaling(i)*S(i,:)
        E(i,:)=scaling(i)*E(i,:)
        D(i,:)=scaling(i)*D(i,:)
        innov(i)=scaling(i)*innov(i)
     enddo

     do j=1,nrobs
        do i=1,nrobs
        R(i,j)=scaling(i)*R(i,j)*scaling(j)
     enddo
     enddo

  endif
!----------------------end pengfei----------------------------------------------

!---------------pengfei--------------
!         open(1999,file='D.dat')
!         do i=1,nrobs
!            write(1999,'(21f18.10)') (D(i,j),j=1,nrens)
!         end do
!         close(1999)


!--------------end pengfei---------------
   lreps=.FALSE.
   if (verbose) print * ,'analysis: verbose is on', mode
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Pseudo inversion of C=SS' +(N-1)*R
   print *,'      analysis: Inversion of C:'
   if (nrobs == 1) then
      nrmin=1
      allocate(W(1,1))
      allocate(eig(1))
      eig(1)=dot_product(S(1,:),S(1,:))+real(nrens-1)*R(1,1)
      eig(1)=1.0/eig(1)
      W(1,1)=1.0

   else
      select case (mode)
      case(11,21)
         nrmin=nrobs
	 print *, mode,'start dgemm'
!         open(1999,file='sk.dat')
!         do i=1,nrobs
!            write(1999,'(21f18.10)') (s(i,j),j=1,nrens)
!         end do
!         close(1999)
!-----------------------pengfei
!------------------------original----------------------------
!        Evaluate R= S*S` + (nrens-1)*R
!         call dgemm('n','t',nrobs,nrobs,nrens, &
!                       1.0, S, nrobs, &
!                            S, nrobs, &
!             real(nrens-1), R, nrobs)
!-------------------------------------------------
         temp2=real(nrens-1)
         temp1=1.0
         call dgemm('n','t',nrobs,nrobs,nrens, &
                       temp1, S, nrobs, &
                            S, nrobs, &
             temp2, R, nrobs)
!         open(1999,file='hphr.dat')
!         do i=1,nrobs
!            write(1999,'(21f18.10)') (R(i,j)/real(nrens-1),j=1,nrobs)
!         end do
!         close(1999)

!----------------end pengfei------------------------------
!-----------------pengfei temp------------------------
!         do i=1,nrens
!            A_pert(:,i)=A(:,i)-sum(A,2)/real(nrens)
!         end do
!         open(1999,file='akf.dat')
!         do i=1,ndim
!            write(1999,'(21f18.10)') (a_pert(i,j),j=1,nrens)
 !        end do
!         close(1999)
!
!print *, "temp output PH`=A'(HA')`/(nrens-1)=A'S`/(nrens-1)"
!         temp2=0.0
!         temp1=1.0/(nrens-1)
!         call dgemm('n','t',ndim,nrobs,nrens, &
!                       temp1,A_pert, ndim, &
!                            S, nrobs, &
!                    temp2, ph, ndim)
!          open(1999,file='ph_t.dat')
!         do i=1,ndim
!            write(1999,'(21f18.10)') (ph(i,j),j=1,nrobs)
!         end do
!         close(1999)
           

!-----------------------------------------------------
         print *, "finish dgemm"
!        Compute eigenvalue decomposition of R -> W*eig*W` 
         allocate(W(nrobs,nrobs))
         allocate(eig(nrobs))
         print *, 'start eigC'
         call eigC(R,nrobs,W,eig)
         print *, 'finish eigc and start eigsign'
         call eigsign(eig,nrobs,truncation)
         print *, 'eigsign'

      case(12,22)
         nrmin=min(nrobs,nrens)
         allocate(W(nrobs,nrmin))
         allocate(eig(nrmin))
         call lowrankCinv(S,R,nrobs,nrens,nrmin,W,eig,truncation)

      case(13,23)
         nrmin=min(nrobs,nrens)
         allocate(W(nrobs,nrmin))
         allocate(eig(nrmin))
         call lowrankE(S,E,nrobs,nrens,nrmin,W,eig,truncation)

      case default
         print *,'analysis: Unknown mode: ',mode
         stop
      end select
   endif








!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Generation of X5 (or representers in EnKF case with few measurements)
   print *,'      analysis: Generation of X5:'
   select case (mode)
   case(11,12,13)
      allocate(X3(nrobs,nrens))
      if (nrobs > 1) then
         call genX3(nrens,nrobs,nrmin,eig,W,D,X3)
      else
         X3=D*eig(1)
      endif

      if (2_8*ndim*nrobs < 1_8*nrens*(nrobs+ndim)) then
!        Code for few observations ( m<nN/(2n-N) )
         if (verbose) print * ,'analysis: Representer approach is used'
         lreps=.true.
         allocate (Reps(ndim,nrobs))
!        Reps=matmul(A,transpose(S))
!-------------------------------------pengfei-------------------------------------------
!         call dgemm('n','t',ndim,nrobs,nrens,1.0,A,ndim,S,nrobs,0.0,Reps,ndim) !original    
     temp1=1.0;temp2=0.0 
     call dgemm('n','t',ndim,nrobs,nrens,temp1,A,ndim,S,nrobs,temp2,Reps,ndim)
         open(1999,file='reps.dat')
         do i=1,ndim
            write(1999,'(21f18.10)') (reps(i,j)/real(nrens-1),j=1,nrobs)
         end do
         close(1999)

!------------------------------------end pengfei--------------------------------------------
      else
         if (verbose) print * ,'analysis: X5 approach is used'
!        X5=matmul(transpose(S),X3)
!-----------------------------------pengfei-----------------------------------------------------
!         call dgemm('t','n',nrens,nrens,nrobs,1.0,S,nrobs,X3,nrobs,0.0,X5,nrens) !original
         temp1=1.0;temp2=0.0
         call dgemm('t','n',nrens,nrens,nrobs,temp1,S,nrobs,X3,nrobs,temp2,X5,nrens)
!-----------------------------end pengfei-----------------------------------------------------------
         do i=1,nrens
            X5(i,i)=X5(i,i)+1.0
         enddo
      endif

   case(21,22,23)
! Mean part of X5
      call meanX5(nrens,nrobs,nrmin,S,W,eig,innov,X5)

! Generating X2
      allocate(X2(nrmin,nrens))
      call genX2(nrens,nrobs,nrmin,S,W,eig,X2)

! Generating X5 matrix
      call X5sqrt(X2,nrobs,nrens,nrmin,X5,update_randrot,mode)

   case default
      print *,'analysis: Unknown flag for mode: ',mode
      stop
   end select


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Generation of inflation
!   call inflationTEST(X5,nrens)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Final ensemble update
   print *,'      analysis: Final ensemble update:'
   print *, 'lreps is :' ,lreps
   if (lreps) then
!     A=A+matmul(Reps,X3)
!----------------------------------pengfei----------------------------------------------
       print *, 'starting dgemm'
!      call dgemm('n','n',ndim,nrens,nrobs,1.0,Reps,ndim,X3,nrobs,1.0,A,ndim) !original
      temp1=1.0;temp2=1.0
      call dgemm('n','n',ndim,nrens,nrobs,temp1,Reps,ndim,X3,nrobs,temp2,A,ndim)
      print *, 'finish dgemm and start dumpX3'
!---------------------------------end pengfei-----------------------------------------------------
      call dumpX3(X3,S,nrobs,nrens)
      print *, 'finish dumpX3'
   else
      iblkmax=min(ndim,200)
      call multa(A, X5, ndim, nrens, iblkmax )
      call dumpX5(X5,nrens)
   endif



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   if (allocated(X2))    deallocate(X2)
   if (allocated(X3))    deallocate(X3)
   if (allocated(eig))   deallocate(eig)
   if (allocated(W))     deallocate(W)
   if (allocated(Reps))  deallocate(Reps)
end subroutine
subroutine randrot(Q,nrens)
   implicit none
   integer, intent(in)  :: nrens
   real(dp),    intent(out) :: Q(nrens,nrens)

   real(dp), dimension(nrens,nrens) ::  A, B
   real(dp) sigma(nrens), work(10*nrens)
   real(dp), parameter :: pi=3.14159253589
   integer ierr
   real(dp) meanB

   call random_number(B)
   call random_number(A)
   Q = sqrt(-2.*log(A+tiny(A))) * cos(2.*pi*B)

!$OMP CRITICAL
! QR factorization
   call dgeqrf(nrens, nrens, Q, nrens, sigma, work, 10*nrens, ierr )
   if (ierr /= 0) print *, 'randrot: dgeqrf ierr=',ierr

! Construction of Q
   call dorgqr(nrens, nrens, nrens, Q, nrens, sigma, work, 10*nrens, ierr )
   if (ierr /= 0) print *, 'randrot: dorgqr ierr=',ierr
!$OMP END CRITICAL


end subroutine randrot
subroutine mean_preserving_rotation(Up,nrens)
! Generates the mean preserving random rotation for the EnKF SQRT algorithm
! using the algorithm from Sakov 2006-07.  I.e, generate rotation Up suceh that
! Up*Up^T=I and Up*1=1 (all rows have sum = 1)  see eq 17.
! From eq 18,    Up=B * Upb * B^T 
! B is a random orthonormal basis with the elements in the first column equals 1/sqrt(nrens)
! Upb = | 1  0 |
!       | 0  U |
! where U is an arbitrary orthonormal matrix of dim nrens-1 x nrens-1  (eq. 19)


!use m_randrot !pengfei mark

implicit none

integer, intent(in)    :: nrens
real(dp),    intent(out)   :: Up(nrens,nrens)

real(dp) B(nrens,nrens),Q(nrens,nrens),R(nrens,nrens),U(nrens-1,nrens-1), Upb(nrens,nrens)

integer j,k
!-----------------pengfei----------------
real(dp) :: temp1,temp2
!-----------------end pengfei-------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Generating the B matrix
! Starting with a random matrix with the correct 1st column
   call random_number(B)
   B(:,1)=1.0/sqrt(real(nrens))

! modified_gram_schmidt is used to create the orthonormal basis
!   do k=1,nrens
!      R(k,k)=sqrt(dot_product(B(:,k),B(:,k)))
!      Q(:,k)=B(:,k)/R(k,k)
!      do j=k+1,nrens
!         R(k,j)=dot_product(Q(:,k),B(:,j))
!         B(:,j)=B(:,j)- Q(:,k)*R(k,j)
!      enddo
!   enddo
!   B=Q

! with overwriting of B
   do k=1,nrens
      R(k,k)=sqrt(dot_product(B(:,k),B(:,k)))
      B(:,k)=B(:,k)/R(k,k)
      do j=k+1,nrens
         R(k,j)=dot_product(B(:,k),B(:,j))
         B(:,j)=B(:,j)- B(:,k)*R(k,j)
      enddo
   enddo


! Check on orthonormality of B
!   do k=1,nrens
!      do j=k,min(k+14,nrens)
!         write(*,'(15f10.4)',advance='no')dot_product(B(:,j),B(:,k))
!      enddo
!      write(*,*)' '
!   enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Creating the orthonormal nrens-1 x nrens-1 U matrix
   call randrot(U,nrens-1)
!! Check on orthonormality of U
!   do k=1,nrens-1
!      do j=k,min(k+14,nrens-1)
!         write(*,'(15f10.4)',advance='no')dot_product(U(:,j),U(:,k))
!      enddo
!      write(*,*)' '
!   enddo
!   stop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Creating the orthonormal nrens x nrens Upb matrix
  Upb(2:nrens,2:nrens)=U(1:nrens-1,1:nrens-1)
  Upb(1,1)=1.0
  Upb(2:nrens,1)=0.0
  Upb(1,2:nrens)=0.0

! Check on orthonormality of Upb
!   do k=1,nrens
!      do j=k,min(k+14,nrens)
!         write(*,'(15f10.4)',advance='no')dot_product(Upb(:,j),Upb(:,k))
!      enddo
!      write(*,*)' '
!   enddo
!   stop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Creating the random orthonormal mean preserving nrens x nrens Upb matrix: Up=B^T Upb B
!---------------------pengfei--------------------------------------------------
!   call dgemm('n','n',nrens,nrens,nrens,1.0,B,nrens,Upb,nrens,0.0,Q,nrens)
!   call dgemm('n','t',nrens,nrens,nrens,1.0,Q,nrens,B,nrens,0.0,Up,nrens)
temp1=1.0;temp2=0.0
   call dgemm('n','n',nrens,nrens,nrens,temp1,B,nrens,Upb,nrens,temp2,Q,nrens)
   call dgemm('n','t',nrens,nrens,nrens,temp1,Q,nrens,B,nrens,temp2,Up,nrens)
!---------------------end pengfei----------------------------------------------------


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Checks
!   do k=1,nrens
!      print *,'Up row sum: ',k,sum(Up(k,1:nrens))
!   enddo

! Check on orthonormality of Up
!   do k=1,nrens
!      do j=k,min(k+14,nrens)
!         write(*,'(15f10.4)',advance='no')dot_product(Up(:,j),Up(:,k))
!      enddo
!      write(*,*)' '
!   enddo
!   stop


end subroutine 
subroutine multa(A, X, ndim, nrens, iblkmax)
implicit none
integer, intent(in) :: ndim
integer, intent(in) :: nrens
integer, intent(in) :: iblkmax
real(dp), intent(in)    :: X(nrens,nrens)
real(dp), intent(inout) :: A(ndim,nrens)
real(dp) v(iblkmax,nrens)  ! Automatic work array

integer ia,ib
!--------pengfei-------------
real(dp) :: temp1,temp2
!--------end pengfei-------------------

do ia = 1,ndim,iblkmax
  ib = min(ia+iblkmax-1,ndim)
  v(1:ib-ia+1,1:nrens) = A(ia:ib,1:nrens)
!---------------------pengfei-------------------------
!  call dgemm('n','n', ib-ia+1, nrens, nrens, &
!              1.0, v(1,1), iblkmax, &
!              X(1,1), nrens, &
!              0.0, A(ia,1), ndim) !original
temp1=1.0;temp2=0.0
  call dgemm('n','n', ib-ia+1, nrens, nrens, &
              temp1, v(1,1), iblkmax, &
              X(1,1), nrens, &
              temp2, A(ia,1), ndim) !original
!-----------------------end pengfei----------------------------
enddo
end subroutine multa
!------------------------------

subroutine lowrankE(S,E,nrobs,nrens,nrmin,W,eig,truncation)
   implicit none
   integer, intent(in)  :: nrobs
   integer, intent(in)  :: nrens
   integer, intent(in)  :: nrmin
   real(dp),    intent(in)  :: S(nrobs,nrens)
   real(dp),    intent(in)  :: E(nrobs,nrens)
   real(dp),    intent(out) :: W(nrobs,nrmin)
   real(dp),    intent(out) :: eig(nrmin)
   real(dp),    intent(in)  :: truncation

   real(dp) U0(nrobs,nrmin),sig0(nrmin)
   real(dp) X0(nrmin,nrens)
   integer i,j

   real(dp) U1(nrmin,nrmin),VT1(1,1)
   real(dp), allocatable :: work(:)
   integer lwork
   integer ierr
!--------pengfei-------------
real(dp) :: temp1,temp2
!--------end pengfei-------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute SVD of S=HA`  ->  U0, sig0
   call  svdS(S,nrobs,nrens,nrmin,U0,sig0,truncation)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute X0=sig0^{*T} U0^T E 

! X0= U0^T R
!--------------------------------pengfei--------------------------------------------------
!   call dgemm('t','n',nrmin,nrens,nrobs, 1.0,U0,nrobs, E,nrobs, 0.0,X0,nrmin) !original
   temp1=1.0;temp2=0.0
   call dgemm('t','n',nrmin,nrens,nrobs, temp1,U0,nrobs, E,nrobs, temp2,X0,nrmin)
!-------------------------------end pengfei----------------------------------------------


   do j=1,nrens
   do i=1,nrmin
      X0(i,j)=sig0(i)*X0(i,j)
   enddo
   enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute singular value decomposition  of X0(nrmin,nrens)
   lwork=2*max(3*nrens+nrobs,5*nrens)
   allocate(work(lwork))
   eig=0.0

   call dgesvd('S', 'N', nrmin, nrens, X0, nrmin, eig, U1, nrmin, VT1, 1, work, lwork, ierr)
   deallocate(work)
   if (ierr /= 0) then
      print *,'mod_anafunc (lowrankE): ierr from call dgesvd 1= ',ierr; stop
   endif

   do i=1,nrmin
      eig(i)=1.0/(1.0+eig(i)**2)
   enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! W = U0 * sig0^{-1} * U1
   do j=1,nrmin
   do i=1,nrmin
      U1(i,j)=sig0(i)*U1(i,j)
   enddo
   enddo
!------------------------pengfei-------------------------------------------------
!   call dgemm('n','n',nrobs,nrmin,nrmin, 1.0,U0,nrobs, U1,nrmin, 0.0,W,nrobs) !original
    temp1=1.0;temp2=0.0
   call dgemm('n','n',nrobs,nrmin,nrmin,temp1,U0,nrobs, U1,nrmin,temp2,W,nrobs)
!-----------------------end pengfei-------------------------------------------------


end subroutine


subroutine eigC(R,nrobs,Z,eig)
! Compute eigenvalue decomposition of R -> Z*eig*Z` 
   integer, intent(in) :: nrobs
   real(dp), intent(in)    :: R(nrobs,nrobs)
   real(dp), intent(out) :: Z(nrobs,nrobs)
   real(dp), intent(out)   :: eig(nrobs)

#ifdef IBM
   real(dp), allocatable :: ap(:)
   integer k
#endif
   real(dp)  RR(nrobs,nrobs)

   real(dp) fwork(8*nrobs)  
   integer iwork(5*nrobs)
   integer ifail(nrobs)
   real(dp) abstol,ddum
   integer idum,neig,ierr
   real(dp), external :: DLAMCH

   idum=1

#ifdef IBM
! Upper packed storage as in ESSL manual
   allocate (ap(nrobs*(nrobs+1)/2) )
   k=0
   do j=1,nrobs
   do i=1,j
       k=k+1
       ap(k)=R(i,j)
   enddo
   enddo
   call dspev(21,ap,eig,Z,nrobs,nrobs,fwork,2*nrobs)
   deallocate(ap)
#else
   abstol=2.0*DLAMCH('S')
   RR=R
   call dsyevx('V', 'A', 'U', nrobs, RR, nrobs, ddum, ddum, idum, idum, abstol, &
            neig, eig, Z, nrobs, fwork, 8*nrobs, iwork, ifail, ierr )
!   print *, 'eig',shape(eig)
 !  write(*,*)  eig
!   print *, 'Z',shape(z)
  ! write(*,'(9f13.6)')  z
   
#endif


end subroutine



subroutine eigsign(eig,nrobs,truncation)
! Returns the inverse of the truncated eigenvalue spectrum
implicit none
integer, intent(in)    :: nrobs
real(dp),    intent(inout) :: eig(nrobs)
real(dp),    intent(in)    :: truncation

integer i,nrsigma
real(dp) sigsum,sigsum1
logical ex

   inquire(file='eigenvalues.dat',exist=ex)
   if (ex) then
      open(10,file='eigenvalues.dat',position='append')
         write(10,'(a,i5,a)')' ZONE  F=POINT, I=',nrobs,' J=1 K=1'
         do i=1,nrobs
            write(10,'(i3,g13.5)')i,eig(nrobs-i+1)
         enddo
      close(10)
   else
      open(10,file='eigenvalues.dat')
         write(10,*)'TITLE = "Eigenvalues of C"'
         write(10,*)'VARIABLES = "obs" "eigenvalues"'
         write(10,'(a,i5,a)')' ZONE  F=POINT, I=',nrobs,' J=1 K=1'
         do i=1,nrobs
            write(10,'(i3,g13.5)')i,eig(nrobs-i+1)
         enddo
      close(10)
   endif

! Significant eigenvalues
   sigsum=sum( eig(1:nrobs) )
   sigsum1=0.0
   nrsigma=0
   do i=nrobs,1,-1
!      print '(a,i5,g13.5)','Eigen values: ',i,eig(i)
      if (sigsum1/sigsum < truncation) then
         nrsigma=nrsigma+1
         sigsum1=sigsum1+eig(i)
         eig(i) = 1.0/eig(i)
      else
         eig(1:i)=0.0
         exit
      endif
   enddo
   write(*,'(2(a,i5))')      '       analysis: Number of dominant eigenvalues: ',nrsigma,' of ',nrobs
   write(*,'(2(a,g13.4),a)') '       analysis: Share (and truncation)        : ',sigsum1/sigsum,' (',truncation,')'


end subroutine



subroutine genX2(nrens,nrobs,idim,S,W,eig,X2)
! Generate X2= (I+eig)^{-0.5} * W^T * S
   implicit none
   integer, intent(in) :: nrens
   integer, intent(in) :: nrobs
   integer, intent(in) :: idim ! idim=nrobs for A4 and nrmin for A5
   real(dp), intent(in)    :: W(idim,nrens)
   real(dp), intent(in)    :: S(nrobs,nrens)
   real(dp), intent(in)    :: eig(idim)
   real(dp), intent(out)   :: X2(idim,nrens)
   integer i,j
!--------pengfei-------------
real(dp) :: temp1,temp2
!--------end pengfei-------------------
!-------------------------pengfei------------------------------------------------
!  call dgemm('t','n',idim,nrens,nrobs,1.0,W,nrobs, S,nrobs, 0.0,X2,idim) !orginal
   temp1=1.0;temp2=0.0
   call dgemm('t','n',idim,nrens,nrobs,temp1,W,nrobs, S,nrobs, temp2,X2,idim)

!-------------------------end pengfei--------------------------------------------

   do j=1,nrens
   do i=1,idim
      X2(i,j)=sqrt(eig(i))*X2(i,j)
   enddo
   enddo

end subroutine



subroutine genX3(nrens,nrobs,nrmin,eig,W,D,X3)
   implicit none
   integer, intent(in) :: nrens
   integer, intent(in) :: nrobs
   integer, intent(in) :: nrmin
   real(dp),    intent(in) :: eig(nrmin)
   real(dp),    intent(in) :: W(nrobs,nrmin)
   real(dp),    intent(in) :: D(nrobs,nrens)
   real(dp),    intent(out) :: X3(nrobs,nrmin)

   real(dp) X1(nrmin,nrobs)
   real(dp) X2(nrmin,nrens)
!--------pengfei-------------
real(dp) :: temp1,temp2
!--------end pengfei-------------------
!-----------------pengfei temp use-----------------
!      real(dp) XX(nrobs,nrobs)

!-------------------end pengfei---------------------
   integer i,j
   do i=1,nrmin
   do j=1,nrobs
      X1(i,j)=eig(i)*W(j,i)
   enddo
   enddo
!-----------------------------pengfei------------------------------------
!---------------pengfei-------------------
!print *, 'temp in genX3 for print value, please mark later,does not affect program'
!     XX=matmul(W,X1)
!       temp1=1.0;temp2=0.0
!      call dgemm('n','n',nrobs,nrobs,nrmin,temp1,W ,nrobs,X1,nrmin,temp2,XX,nrobs)
!          open(1999,file='hphr_inverse.dat')
!         do i=1,nrobs
!            write(1999,'(21f18.10)') (XX(i,j)*real(nrens-1),j=1,nrobs)
!         end do
!         close(1999)

!-----------------------------end pengfei----------------------------------
!------------------------------pengfei------------------------------------------------------
!     X2=matmul(X1,D)
!      call dgemm('n','n',nrmin,nrens,nrobs,1.0,X1,nrmin,D ,nrobs,0.0,X2,nrmin) !original
       temp1=1.0;temp2=0.0
      call dgemm('n','n',nrmin,nrens,nrobs,temp1,X1,nrmin,D ,nrobs,temp2,X2,nrmin)
!     X3=matmul(W,X2)
!      call dgemm('n','n',nrobs,nrens,nrmin,1.0,W ,nrobs,X2,nrmin,0.0,X3,nrobs) !original
      call dgemm('n','n',nrobs,nrens,nrmin,temp1,W ,nrobs,X2,nrmin,temp2,X3,nrobs)
!-----------------------------end pengfei--------------------------------------------------
end subroutine



subroutine meanX5(nrens,nrobs,nrmin,S,W,eig,innov,X5)
   implicit none
   integer, intent(in) :: nrens
   integer, intent(in) :: nrobs
   integer, intent(in) :: nrmin
   real(dp), intent(in)    :: W(nrmin,nrmin)
   real(dp), intent(in)    :: S(nrobs,nrens)
   real(dp), intent(in)    :: eig(nrmin)
   real(dp), intent(in)    :: innov(nrobs)
   real(dp), intent(out)   :: X5(nrens,nrens)

   real(dp) y1(nrmin)
   real(dp) y2(nrmin)
   real(dp) y3(nrobs)
   real(dp) y4(nrens) 
   integer i
!--------pengfei-------------
real(dp) :: temp1,temp2
!--------end pengfei-------------------
   if (nrobs==1) then
      y1(1)=W(1,1)*innov(1)
      y2(1)=eig(1)*y1(1)
      y3(1)=W(1,1)*y2(1)
      y4(:)=y3(1)*S(1,:)
   else
!------------------------------pengfei---------------------------------------
!      call dgemv('t',nrobs,nrmin,1.0,W,nrobs,innov,1,0.0,y1 ,1) !original
      temp1=1.0;temp2=0.0
      call dgemv('t',nrobs,nrmin,temp1,W,nrobs,innov,1,temp2,y1 ,1) 
      y2=eig*y1  
!      call dgemv('n',nrobs,nrmin,1.0,W ,nrobs,y2,1,0.0,y3 ,1) !original
      call dgemv('n',nrobs,nrmin,temp1,W ,nrobs,y2,1,temp2,y3 ,1) 
!      call dgemv('t',nrobs,nrens,1.0,S ,nrobs,y3,1,0.0,y4 ,1) !original
     call dgemv('t',nrobs,nrens,temp1,S ,nrobs,y3,1,temp2,y4 ,1)
!-----------------------------end pengfei------------------------------------------
   endif

   do i=1,nrens
      X5(:,i)=y4(:)
   enddo

! X5=enN + (I - enN) X5  = enN + X5
   X5=1.0/real(nrens) + X5

end subroutine



subroutine X5sqrt(X2,nrobs,nrens,nrmin,X5,update_randrot,mode)
!   use m_randrot
!   use m_mean_preserving_rotation
   implicit none
   integer, intent(in) :: nrobs
   integer, intent(in) :: nrens
   integer, intent(inout) :: nrmin ! note that nrmin=nrobs in a4
   real(dp), intent(in)    :: X2(nrmin,nrens)
   real(dp), intent(inout) :: X5(nrens,nrens)
   logical, intent(in) :: update_randrot
   integer, intent(in) :: mode

   real(dp) X3(nrens,nrens)
   real(dp) X33(nrens,nrens)
   real(dp) X4(nrens,nrens)
   real(dp) IenN(nrens,nrens)
   real(dp), save, allocatable :: ROT(:,:)


   real(dp) U(nrmin,1),sig(nrmin),VT(nrens,nrens)
   real(dp), allocatable, dimension(:)   :: work,isigma
   integer i,j,lwork,ierr
!--------pengfei-------------
real(dp) :: temp1,temp2
!--------end pengfei-------------------
   print *,'      analysis (X5sqrt): update_randrot= ',update_randrot
   if (update_randrot) then
      if (allocated(ROT)) deallocate(ROT)
      allocate(ROT(nrens,nrens))
!!      ROT=0.0
!!      do i=1,nrens
!!         ROT(i,i)=1.0
!!      enddo
!      call randrot(ROT,nrens)  ! old non mean preserving random rotation
      call mean_preserving_rotation(ROT,nrens)
   endif

! SVD of X2
   lwork=2*max(3*nrens+nrens,5*nrens); allocate(work(lwork))
   sig=0.0
   call dgesvd('N', 'A', nrmin, nrens, X2, nrmin, sig, U, nrmin, VT, nrens, work, lwork, ierr)
   deallocate(work)
   if (ierr /= 0) then
      print *,'X5sqrt: ierr from call dgesvd = ',ierr
      stop
   endif


   if (mode == 21) nrmin=min(nrens,nrobs)
   allocate(isigma(nrmin))
   isigma=1.0
   do i=1,nrmin
      if ( sig(i) > 1.0 ) print *,'X5sqrt: WARNING (m_X5sqrt): sig > 1',i,sig(i)
      isigma(i)=sqrt(dmax1(1.0-sig(i)**2,0.0_DP) )
   enddo

   do j=1,nrens
      X3(:,j)=VT(j,:)
   enddo


   do j=1,nrmin
      X3(:,j)=X3(:,j)*isigma(j)
   enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Multiply  X3* V' = (V*sqrt(I-sigma*sigma) * V' to ensure symmetric sqrt and 
! mean preserving rotation.   Sakov paper eq 13
!------------------------------------pengfei---------------------------------------
!   call dgemm('n','n',nrens,nrens,nrens,1.0,X3,nrens,VT,nrens,0.0,X33,nrens) !original
      temp1=1.0;temp2=0.0
      call dgemm('n','n',nrens,nrens,nrens,temp1,X3,nrens,VT,nrens,temp2,X33,nrens) 
!-----------------------------------end pengfei-----------------------------------------
! Check mean preservation X33*1_N = a* 1_N (Sakov paper eq 15)
!   do i=1,nrens
!      print *,'sum(X33)= ',i,sum(X33(i,:))
!   enddo
   
!!X33=X3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!   print '(a)','X5sqrt: sig: '
!   print '(5g11.3)',sig(1:min(nrmin,nrens))
!---------------------------pengfei------------------------------------------------
!   call dgemm('n','n',nrens,nrens,nrens,1.0,X33,nrens,ROT,nrens,0.0,X4,nrens) !original
   temp1=1.0;temp2=0.0
   call dgemm('n','n',nrens,nrens,nrens,temp1,X33,nrens,ROT,nrens,temp2,X4,nrens)
!---------------------------end pengfei-----------------------------------------------------

   IenN=-1.0/real(nrens)
   do i=1,nrens
      IenN(i,i)=  IenN(i,i) + 1.0
   enddo
!-----------------------------------pengfei-------------------------------------------
!   call dgemm('n','n',nrens,nrens,nrens,1.0,IenN,nrens,X4,nrens,1.0,X5,nrens) !original
   temp1=1.0;temp2=1.0
   call dgemm('n','n',nrens,nrens,nrens,temp1,IenN,nrens,X4,nrens,temp2,X5,nrens)
!------------------------------------end pengfei-----------------------------------------------

   deallocate(isigma)

end subroutine



subroutine dumpX3(X3,S,nrobs,nrens)
   implicit none
   integer, intent(in) :: nrens
   integer, intent(in) :: nrobs
   real(dp),    intent(in) :: X3(nrens,nrens)
   real(dp),    intent(in) :: S(nrobs,nrens)
   character(len=2) :: tag2

   tag2(1:2)='X3'
   open(10,file='X5.uf',form='unformatted')
      write(10)tag2,nrens,nrobs,X3,S
   close(10)

end subroutine



subroutine dumpX5(X5,nrens)
   implicit none
   integer, intent(in) :: nrens
   real(dp),    intent(in) :: X5(nrens,nrens)
   integer j
   character(len=2) :: tag2

   tag2(1:2)='X5'
   open(10,file='X5.uf',form='unformatted')
      write(10)tag2,nrens,X5
   close(10)

   open(10,file='X5col.dat')
      do j=1,nrens
         write(10,'(i5,f10.4)')j,sum(X5(:,j))
      enddo
   close(10)

   open(10,file='X5row.dat')
      do j=1,nrens
         write(10,'(i5,f10.4)')j,sum(X5(j,:))/real(nrens)
       enddo
   close(10)
end subroutine



subroutine lowrankCinv(S,R,nrobs,nrens,nrmin,W,eig,truncation)
   implicit none
   integer, intent(in)  :: nrobs
   integer, intent(in)  :: nrens
   integer, intent(in)  :: nrmin
   real(dp),    intent(in)  :: S(nrobs,nrens)
   real(dp),    intent(in)  :: R(nrobs,nrobs)
   real(dp),    intent(out) :: W(nrobs,nrmin)
   real(dp),    intent(out) :: eig(nrmin)
   real(dp),    intent(in)  :: truncation

   real(dp) U0(nrobs,nrmin),sig0(nrmin)
   real(dp) B(nrmin,nrmin),Z(nrmin,nrmin)
   integer i,j
!--------pengfei-------------
real(dp) :: temp1,temp2
!--------end pengfei-------------------
! Compute SVD of S=HA`  ->  U0, sig0
   call  svdS(S,nrobs,nrens,nrmin,U0,sig0,truncation)

! Compute B=sig0^{-1} U0^T R U0 sig0^{-1}
   call lowrankCee(B,nrmin,nrobs,nrens,R,U0,sig0)

! Compute eigenvalue decomposition  of B(nrmin,nrmin)
   call eigC(B,nrmin,Z,eig)

!   print *,'eig:',nrmin
!   print '(6g11.3)',eig

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute inverse diagonal of (I+Lamda)
   do i=1,nrmin
      eig(i)=1.0/(1.0+eig(i))
   enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! W = U0 * sig0^{-1} * Z
   do j=1,nrmin
   do i=1,nrmin
      Z(i,j)=sig0(i)*Z(i,j)
   enddo
   enddo
!--------------------------pengfei----------------------------------------------
!   call dgemm('n','n',nrobs,nrmin,nrmin, 1.0,U0,nrobs, Z,nrmin, 0.0,W,nrobs) !original
  temp1=1.0;temp2=0.0
  call dgemm('n','n',nrobs,nrmin,nrmin,temp1,U0,nrobs, Z,nrmin,temp2,W,nrobs)
!-------------------------end pengfei--------------------------------------------------------

end subroutine




subroutine lowrankCee(B,nrmin,nrobs,nrens,R,U0,sig0)
implicit none
integer, intent(in) :: nrmin
integer, intent(in) :: nrobs
integer, intent(in) :: nrens
real(dp), intent(inout) :: B(nrmin,nrmin)
real(dp), intent(in)    :: R(nrobs,nrobs)
real(dp), intent(in)    :: U0(nrobs,nrmin)
real(dp), intent(in)    :: sig0(nrmin)
real(dp) X0(nrmin,nrobs)
integer  i,j
!--------pengfei-------------
real(dp) :: temp1,temp2
!--------end pengfei-------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute B=sig0^{-1} U0^T R U0 sig0^{-1}

! X0= U0^T R
!------------------------------------pengfei-----------------------------------------
!   call dgemm('t','n',nrmin,nrobs,nrobs, 1.0,U0,nrobs, R,nrobs, 0.0,X0,nrmin) !original
   temp1=1.0;temp2=0.0
   call dgemm('t','n',nrmin,nrobs,nrobs, temp1,U0,nrobs, R,nrobs, temp2,X0,nrmin)
! B= X0 U0
!   call dgemm('n','n',nrmin,nrmin,nrobs, 1.0,X0,nrmin, U0,nrobs, 0.0,B,nrmin) !original
    call dgemm('n','n',nrmin,nrmin,nrobs, temp1,X0,nrmin, U0,nrobs,temp2,B,nrmin)
!------------------------------------end pengfei---------------------------------------------------------
   do j=1,nrmin
   do i=1,nrmin
      B(i,j)=sig0(i)*B(i,j)
   enddo
   enddo

   do j=1,nrmin
   do i=1,nrmin
      B(i,j)=sig0(j)*B(i,j)
   enddo
   enddo

   B=real(nrens-1)*B

end subroutine


subroutine svdS(S,nrobs,nrens,nrmin,U0,sig0,truncation)
   integer, intent(in)  :: nrobs
   integer, intent(in)  :: nrens
   integer, intent(in)  :: nrmin
   real(dp),    intent(in)  :: S(nrobs,nrens)
   real(dp),    intent(out) :: sig0(nrmin)
   real(dp),    intent(in)  :: U0(nrobs,nrmin)
   real(dp),    intent(in)  :: truncation

   real(dp) S0(nrobs,nrens)
   real(dp) VT0(1,1)
   integer ierr
   integer lwork
   real(dp), allocatable, dimension(:)   :: work
   integer nrsigma,i

   real(dp) sigsum,sigsum1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute SVD of S=HA`  ->  U0, sig0
   lwork=2*max(3*nrens+nrobs,5*nrens)
   allocate(work(lwork))

   S0=S
   sig0=0.0
   call dgesvd('S', 'N', nrobs, nrens, S0, nrobs, sig0, U0, nrobs, VT0, nrens, work, lwork, ierr)
   deallocate(work)
   if (ierr /= 0) then
      print *,'svdS: ierr from call dgesvd 0= ',ierr; stop
   endif

   sigsum=0.0
   do i=1,nrmin
      sigsum=sigsum+sig0(i)**2
   enddo

   sigsum1=0.0
! Significant eigenvalues.
   nrsigma=0
   do i=1,nrmin                       
      if (sigsum1/sigsum < truncation) then
         nrsigma=nrsigma+1
         sigsum1=sigsum1+sig0(i)**2
      else
         sig0(i:nrmin)=0.0
         exit
      endif
   enddo

   write(*,'(a,i5,g13.5)') '      analysis svdS: dominant sing. values and share ',nrsigma,sigsum1/sigsum
!   write(*,'(5g11.3)')sig0

   do i=1,nrsigma
       sig0(i) = 1.0/sig0(i)
   enddo

end subroutine


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine SEIK_analysis(ens_p,R,resid_p,HL_p,dim_ens,dim_obs,dim_state,subtype)
   use mod_prec
   implicit none
   integer, intent(in) :: dim_state             ! dimension of model state
   integer, intent(in) :: dim_ens            ! number of ensemble members
   integer, intent(in) :: dim_obs            ! number of observations
   real(dp), intent(inout) :: ens_p(dim_state,dim_ens)    ! ensemble matrix
   real(dp), intent(inout)    :: R(dim_obs,dim_obs)   ! matrix holding R (only used if mode=?1 or ?2)
   real(dp), intent(in)    :: resid_p(dim_obs)     ! vector holding d-H*mean(A)
   real(dp), intent(inout)    :: HL_p(dim_obs,dim_ens)   ! matrix holding HA
   real(dp) :: state_p(dim_state)
   real(dp) :: rdim_ens
   real(dp) ::  invdimens,forget,forget_ana,fac
   integer :: i,j,member,row,rank,dummy,maxblksize,blklower,blkupper
   real(dp),allocatable,dimension(:,:) :: RiHL_p,Uinv,Uinv_p,Uinv_inc,temp_Uinv,TRiHLd,tempUinv,Ttrans,OmegaT,OmegaTsave,Omega,TA,ens_blk
   real(dp),allocatable,dimension(:) :: RiHLd,ipiv,state_inc_p
   INTEGER :: dgesv_info,dpotrf_info,dtrtrs_info,col,Nm1Vsn
!   *** Perform prepoststep for SEIK without re-init (subtype=3)   ***
   INTEGER :: subtype
   logical :: firsttime = .true.
   logical :: storeOmega=.false.
   rank=dim_ens-1
   dummy=0
   forget_ana=1.0
   Nm1Vsn=1
! ***********************************
! *** Compute mean forecast state ***
! ***********************************
  state_p = 0.0
  invdimens = 1.0 / REAL(dim_ens)
  DO member = 1, dim_ens
     DO row = 1, dim_state ! change dim_p to dim_state
        state_p(row) = state_p(row) + invdimens * ens_p(row, member)
     END DO
  END DO
! ******************************************
! ***   Compute analized matrix Uinv     ***
! ***                                    ***
! ***  -1                T        T  -1  ***
! *** U  = forget*N T T + (HL)  R  (HL)  ***
! ***  i                      i  i     i ***
! ******************************************
     ! HL = [Hx_1 ... Hx_N] T
     open(1,file='ha.dat')
     do i=1,dim_obs
        write(1,'(24f18.8)') (hl_p(i,j),j=1,dim_ens)
     end do
     close(1)
     CALL PDAF_seik_matrixT(dim_obs, dim_ens, HL_p)
     open(1,file='hat.dat')
     do i=1,dim_obs
        write(1,'(24f18.8)') (hl_p(i,j),j=1,dim_ens)
     end do
     close(1)
     ! ***                RiHL = Rinv HL                 ***
     ! *** this is implemented as a subroutine thus that ***
     ! *** Rinv does not need to be allocated explicitly ***
     ALLOCATE(RiHL_p(dim_obs, rank))
     open(1,file='r_seik.dat')
     do i=1,dim_obs
        write(1,'(90f18.8)') (R(i,j),j=1,dim_obs)
     end do
     close(1)
     CALL prodRinvA(dummy, dim_obs, rank, R, HL_p, RiHL_p)
     open(1,file='rihlp.dat')
     do i=1,dim_obs
        write(1,'(24f18.8)') (Rihl_p(i,j),j=1,rank)
     end do
     close(1)

     ! *** Initialize Uinv = N T^T T ***
     allocate(Uinv(rank,rank))
     CALL PDAF_seik_Uinv(rank, Uinv)
     open(1,file='uinv1.dat')
     do i=1,rank
        write(1,'(24f18.8)') (Uinv(i,j),j=1,rank)
     end do
     close(1)
     ! *** Finish computation of Uinv  ***
     ! ***   -1          -1    T       ***
     ! ***  U  = forget U  + HL RiHL   ***
     ALLOCATE(Uinv_p(rank, rank))
     ALLOCATE(Uinv_inc(rank, rank))
     CALL dgemm('t', 'n', rank, rank, dim_obs, &
          one_db, HL_p, dim_obs, RiHL_p, dim_obs, &
          zero_db, Uinv_p, rank)  ! Uinv_p=HL_p'*RiHL_p
     Uinv = forget_ana * Uinv + Uinv_p
     open(1,file='uinv_p.dat')
     do i=1,rank
        write(1,'(24f18.8)') (Uinv_p(i,j),j=1,rank)
     end do
     close(1)
     open(1,file='uinv2.dat')
     do i=1,rank
        write(1,'(24f18.8)') (Uinv(i,j),j=1,rank)
     end do
     close(1)
  DEALLOCATE(Uinv_p, Uinv_inc)
! ************************************
! ***      update model state      ***
! ***                              ***
! ***  a   f          T         f  ***
! *** x = x + L U RiHV  (y - H x ) ***
! ***                              ***
! ************************************

  ! ************************
  ! *** RiHLd = RiHV^T d ***
  ! ************************
  ALLOCATE(RiHLd(rank))
     open(1,file='Rihl_p.dat')
     do i=1,dim_obs
        write(1,'(23f18.8)') (Rihl_p(i,j),j=1,rank)
     end do
     close(1)
     CALL dgemv('t', dim_obs, rank, one_db, RiHL_p, &
          dim_obs, resid_p, 1, zero_db, RiHLd, 1)

  ! ****************************************
  ! *** Compute  w = U RiHLd  by solving ***
  ! ***           -1                     ***
  ! ***          U  w = RiHLd            ***
  ! *** for w. We use the LAPACK         ***
  ! *** routine DPOSVX                   ***
  ! ****************************************

  ALLOCATE(temp_Uinv(rank, rank))
  ALLOCATE(ipiv(rank)) ! working on here ,pause
  ! save matrix Uinv
  temp_Uinv = Uinv

  ! call solver (DGESV - LU solver)
  CALL dgesv(rank, 1, temp_Uinv, rank, ipiv, &
       RiHLd, rank, dgesv_info)
  DEALLOCATE(temp_Uinv, ipiv)

 ! *** check if solve was successful
  update: IF (dgesv_info /= 0) THEN
     WRITE (*, '(/2x, a/)') '!!! Problem in solve for state analysis !!!,stop...'
     stop                                        
  ELSE

     ! **************************
     ! *** Compute vector T w ***
     ! **************************

     ALLOCATE(TRiHLd(dim_ens, 1))
!     open(1,file='rihld_ana.dat')
!     do i=1,dim_ens
!        write(1,*) rihld(i)
!     end do
!     close(1)
     CALL PDAF_seik_TtimesA(rank, 1, RiHLd, TRiHLd)
     DEALLOCATE(RiHLd)


     ! **************************
     ! *** Update model state ***
     ! ***    a   f           ***
     ! ***   x = x + LT Tw    ***
     ! **************************

     ALLOCATE(state_inc_p(dim_state))
     CALL dgemv('n', dim_state, dim_ens, one_db, ens_p, &
          dim_state, TRiHLd, 1, zero_db, state_inc_p, 1)
!     open(1,file='TRiHLd_ana.dat')
!     do i=1,dim_ens
!        write(1,*) TRiHLd(i,1)
!     end do
!     close(1)
     DEALLOCATE(TRiHLd)
     state_p = state_p + state_inc_p
     open(1,file='state_ana.dat')
     do i=1,dim_state
        write(1,*) state_p(i)
     end do

     DEALLOCATE(state_inc_p)

  END IF update


!***************************start to resample*********************************888
! ************************************
! *** Compute C^(-1) by Cholesky   ***
! *** decomposition U^(-1) = C C^T ***
! ************************************

  ! initialize Uinv for internal use
  ALLOCATE(tempUinv(rank, rank))
  IF (subtype /= 3) THEN
     tempUinv(:,:) = Uinv(:,:)
  ELSE
     rdim_ens = DBLE(dim_ens)

     ! Initialize matrix T^T
     ALLOCATE(Ttrans(rank, dim_ens))
     DO i = 1, rank
        DO j = 1, dim_ens
           Ttrans(i, j) = -1.0 / rdim_ens
        END DO
     END DO
     DO i = 1, rank
        Ttrans(i, i) = Ttrans(i, i) + 1.0
     END DO
     CALL dgemm('n', 't', rank, rank, dim_ens, &
          rdim_ens, Ttrans, rank, Ttrans, rank, &
          zero_db, tempUinv, rank)
     DEALLOCATE(Ttrans)
  END IF
  ! Cholesky decomposition of tempUinv
   print *, 'tempuinv'
   print *, tempuinv ! pengfei
  CALL dpotrf('l', rank, tempUinv, rank, dpotrf_info)



 ! check if Cholesky decomposition was successful
  CholeskyOK: IF (dpotrf_info == 0) THEN
     ! Decomposition OK, continue


! *************************************************
! *** Generate ensemble of interpolating states ***
! *************************************************

     ! allocate fields
     ALLOCATE(omegaT(rank, dim_ens))

     IF (storeOmega ) THEN
        ALLOCATE(omegaTsave(rank, dim_ens))
     END IF
     Omega_store: IF (storeOmega) THEN
        first: IF (firsttime) THEN
           ! *** At first call to SEIK_RESAMPLE initialize   ***
           ! *** the matrix Omega in SEIK_Omega and store it ***

           ALLOCATE(omega(dim_ens, rank))

           ! *** Generate uniform orthogonal matrix OMEGA ***
           CALL PDAF_seik_omega(rank, Omega, 0)

           ! transpose Omega
           omegaT = TRANSPOSE(Omega)

           ! store transposed Omega
           omegaTsave = omegaT
           firsttime  = .FALSE.

           DEALLOCATE(omega)

        ELSE first
!           IF (mype == 0 .AND. screen > 0) &
                WRITE (*, '(8x, a)') '--- use stored Omega'
           omegaT = omegaTsave
        END IF first

     ELSE Omega_store

        ! *** Initialize the matrix Omega in SEIK_Omega ***
        ! *** each time SEIK_RESAMPLE is called         ***

        ALLOCATE(omega(dim_ens, rank))

        ! *** Generate uniform orthogonal matrix OMEGA ***
        CALL PDAF_seik_omega(rank, Omega, 0)

        ! transpose Omega
        omegaT = TRANSPOSE(Omega)

        DEALLOCATE(omega)
     END IF Omega_store
     ! ***     Generate ensemble of states                             ***
     ! *** x_i = x + sqrt(FAC) L (Omega C^(-1))t                       ***
     ! *** Here FAC depends on the use definition of the covariance    ***
     ! *** matrix P using a factor (r+1)^-1 or r^-1.                   ***

     ! A = (Omega C^(-1)) by solving Ct A = OmegaT for A
     CALL dtrtrs('l', 't', 'n', rank, dim_ens, &
          tempUinv, rank, OmegaT, rank, dtrtrs_info)
     ! check if solve was successful
     solveOK: IF (dtrtrs_info == 0) THEN
        ! Solve for A OK, continue

        ! *** T A' (A' stored in OmegaT) ***
        ALLOCATE(TA(dim_ens, dim_ens))
        CALL PDAF_seik_TtimesA(rank, dim_ens, OmegaT, TA)
        ! *** Block formulation for resampling
        maxblksize = 200
             WRITE (*, '(8x, a, i5)') '--- use blocking with size ', maxblksize

        ALLOCATE(ens_blk(maxblksize, dim_ens))
        blocking: DO blklower = 1, dim_state, maxblksize

           blkupper = MIN(blklower + maxblksize - 1, dim_state)

           ! Store old state ensemble
           DO col = 1, dim_ens
              ens_blk(1 : blkupper - blklower + 1, col) &
                   = ens_p(blklower : blkupper, col)
           END DO

           DO col = 1,dim_ens
              ens_p(blklower : blkupper, col) = state_p(blklower : blkupper)
           END DO
           ! *** X = state+ sqrt(FAC) state_ens T A^T (A^T stored in OmegaT) ***
           IF (Nm1vsN == 1) THEN
              ! Use factor (N-1)^-1
              fac = SQRT(REAL(dim_ens - 1))
           ELSE
              ! Use factor N^-1
              fac = SQRT(REAL(dim_ens))
           END IF

           CALL dgemm('n', 'n', blkupper - blklower + 1, dim_ens, dim_ens, &
                fac, ens_blk(1, 1), maxblksize, TA(1, 1), dim_ens, &
                one_db, ens_p(blklower, 1), dim_state)
        END DO blocking
        DEALLOCATE(ens_blk, TA)

     ELSE SolveOK

        ! Solve for A failed
        WRITE (*, '(/2x, a/)') &
             '!!! Problem with solve for A in SEIK_RESAMPLE !!!'

     ENDIF SolveOK

     DEALLOCATE(OmegaT)

  ELSE CholeskyOK

     ! eigendecomposition failed
     WRITE (*, '(/2x, a/)') &
          '!!! Problem with Cholesky decomposition of Uinv !!!'

  ENDIF CholeskyOK


! ****************
! *** clean up ***
! ****************

  DEALLOCATE(tempUinv)

end subroutine SEIK_analysis
!--------------------------subroutine for SEIK-------------------------------
!$Id$
!BOP
!
! !ROUTINE: PDAF_seik_matrixT --- Operate matrix T on A as AT
!
! !INTERFACE:
SUBROUTINE PDAF_seik_matrixT(dim, dim_ens, A)

! !DESCRIPTION:
! Operate matrix T on another matrix as
!         $A = A T$.
!
! T is a dim_ens x (dim_ens-1) matrix with zero column sums.
! There are two proposed forms of T (ensemble size N):\\
! typeT=0: diag(T)=1-1/N; nondiag(T)=-1/N; 
!          last row= -1/N\\
! typeT=1: diag(T)=1; nondiag(T)=0; last row = -1\\
! We typically use TypeT=0, but both variants are implemented.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2002-01 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
!  USE PDAF_memcounting, &
!       ONLY: PDAF_memcount
 use mod_prec
  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: dim               ! dimension of states
  INTEGER, INTENT(in) :: dim_ens           ! Size of ensemble
  REAL(dp), INTENT(inout) :: A(dim, dim_ens)   ! Input/output matrix

! !CALLING SEQUENCE:
! Called by: PDAF_seik_analysis
! Called by: PDAF_seik_analysis_newT
! Calls PDAF_memcount
!EOP
  
! *** local variables ***
  INTEGER :: row, col  ! counters
  INTEGER :: typeT = 0 ! Which type of T
  REAL(dp) :: invdimens    ! Inverse of ensemble size
  INTEGER, SAVE :: allocflag = 0  ! Flag for dynamic allocation
  REAL(dp), ALLOCATABLE :: rowmean(:) ! Mean values of rows of A


! **********************
! *** INITIALIZATION ***
! **********************

  ALLOCATE(rowmean(dim))
  IF (allocflag == 0) THEN
     ! count allocated memory
!     CALL PDAF_memcount(3, 'r', dim)
     allocflag = 1
  END IF
  rowmean   = 0.0
  invdimens = 1.0 / REAL(dim_ens)

  IF (typeT == 0) THEN

     ! *** Compute row means of A ***
     DO col = 1, dim_ens
        DO row = 1, dim
           rowmean(row) = rowmean(row) + A(row, col)
        END DO
     END DO
     rowmean = invdimens * rowmean

  ELSE

     ! *** Get last column of A ***
     DO row = 1, dim
        rowmean(row) = A(row, dim_ens)
     END DO
     
  END IF


! **********************************************
! ***  Operate T on A                        ***
! ***                                        ***
! *** v^TT = (v_1-mean(v), ... ,v_r-mean(v)) ***
! *** with v = (v_1,v_2, ... ,r_(r+1))       ***
! **********************************************

  DO col = 1, dim_ens - 1
     DO row = 1, dim
        A(row, col) = A(row, col) - rowmean(row)
     END DO
  END DO
  
  DO row = 1, dim
     A(row, dim_ens) = 0.0
  END DO


! ********************
! *** FINISHING UP ***
! ********************

  DEALLOCATE(rowmean)

END SUBROUTINE PDAF_seik_matrixT

SUBROUTINE prodRinvA(step, dim_obs_p, rank, obs_p, A_p, C_p)

! !DESCRIPTION:
! User-supplied routine for PDAF (SEEK, SEIK):
!
! The routine is called during the analysis step.
! It has to compute the product of the inverse of 
! the observation error covariance matrix with
! the matrix of observed EOF modes (SEEK) or 
! observed ensemble perturbations (SEIK).
!
! This routine is called by all filter processes.
!
! Implementation for the dummy model with domain
! decomposition. Here, we assume a diagonal observation
! error covariance matrix with constant variances. 
! Thus, the product can be implemented efficiently 
! as a scaling of each element of the input matrix
! by the inverse variance.
!
! !REVISION HISTORY:
! 2004-10 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
!  USE mod_assimilation, &
!       ONLY: rms_obs
  USE MOD_prec
  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step                ! Current time step
  INTEGER, INTENT(in) :: dim_obs_p           ! PE-local dimension of obs. vector
  INTEGER, INTENT(in) :: rank                ! Rank of initial covariance matrix
  REAL(dp), INTENT(in)    :: obs_p(dim_obs_p,dim_obs_p)    ! R
  REAL(dp), INTENT(in)    :: A_p(dim_obs_p,rank) ! Input matrix from SEEK_ANALYSIS
  REAL(dp), INTENT(out)   :: C_p(dim_obs_p,rank) ! Output matrix

! !CALLING SEQUENCE:
! Called by: PDAF_seek_analysis        (as U_prodRinvA)
! Called by: PDAF_seik_analysis        (as U_prodRinvA)
! Called by: PDAF_seik_analysis_newT   (as U_prodRinvA)
!EOP

! *** local variables ***
  INTEGER :: i, j       ! index of observation component
  REAL(dp) :: ivariance_obs ! inverse of variance of the observations


! **********************
! *** INITIALIZATION ***
! **********************
  
  ! *** initialize numbers
!  ivariance_obs = 1.0 / rms_obs ** 2 ! original

! *************************************
! ***                -1             ***
! ***           C = R   A           ***
! ***                               ***
! *** The inverse observation error ***
! *** covariance matrix is not      ***
! *** computed explicitely.         ***
! *************************************

!  DO j = 1, rank
!     DO i = 1, dim_obs_p
!        C_p(i, j) = ivariance_obs * A_p(i, j)
!     END DO
!  END DO
!--------------------end original-------------
     DO i = 1, dim_obs_p
         DO j = 1, rank
        ivariance_obs=1/(obs_p(i,i))
        C_p(i, j) = ivariance_obs * A_p(i, j)
     END DO
  END DO


END SUBROUTINE prodRinvA
!$Id$
!BOP
!
! !ROUTINE: PDAF_seik_Uinv - Initialize matrix Uinv from matrix T
!
! !INTERFACE:
SUBROUTINE PDAF_seik_Uinv(rank, Uinv)

! !DESCRIPTION:
! Initialize matrix Uinv by
! $U^{-1} = FAC\ T^T T$
! where $FAC$ = rank+1 for a covariance matrix with factor
! (rank+1)$^{-1}$ and $FAC$ = rank for a covariance matrix
! with factor rank$^{-1}$.
!
! There are two proposed forms of T (ensemble size N):\\
! typeT=0: diag(T)=1-1/N; nondiag(T)=-1/N; 
!          last row= -1/N\\
! typeT=1: diag(T)=1; nondiag(T)=0; last row = -1\\
! We typically use TypeT=0, but both variants are implemented.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2002-01 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
!  USE PDAF_mod_filter, &
!       ONLY: Nm1vsN
  use mod_prec
  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: rank             ! Rank of initial covariance matrix
  REAL(dp), INTENT(inout) :: Uinv(rank, rank) ! Inverse of matrix U
!EOP
  INTEGER :: Nm1vsN        ! Flag which definition of P ist used in SEIK
  ! (0): Factor N^-1; (1): Factor (N-1)^-1 - Recommended is 1 for
  ! a real ensemble filter, 0 is for compatibility with older PDAF versions
  
! *** local variables ***
  INTEGER :: row, col       ! counters
  INTEGER :: typeT = 0      ! Choose type of T
  REAL(dp) :: rdivrp1, r2divrp1 ! scaling factors for Uinv

Nm1vsN=1
! ***********************
! *** Initialize Uinv ***
! ***********************

  ttype: IF (typeT == 0) THEN

     ! Scaling factors
     IF (Nm1vsN == 1) THEN
        ! For ensemble covariance matrix with factor (N-1)^-1
        rdivrp1 = REAL(rank) / REAL(rank + 1)
        r2divrp1 = REAL(rank)**2 / REAL(rank + 1)
     ELSE
        ! For ensemble covariance matrix with factor N^-1
        rdivrp1 = 1
        r2divrp1 = REAL(rank)
     END IF

     DO col = 1, rank
        ! non-diagonal elements - upper triangle
        DO row = 1, col - 1
           Uinv(row, col) = - rdivrp1
        END DO
        ! diagonal
        Uinv(col, col) = r2divrp1
        ! non-diagonal elements - lower triangle
        DO row = col + 1, rank
           Uinv(row, col) = -rdivrp1
        END DO
     END DO

  ELSE ttype

     ! Scaling factors
     IF (Nm1vsN == 1) THEN
        ! For ensemble covariance matrix with factor (N-1)^-1
        rdivrp1 = REAL(rank) / REAL(rank + 1)
     ELSE
        ! For ensemble covariance matrix with factor N^-1
        rdivrp1 = 1
     END IF

     DO col = 1, rank
        ! non-diagonal elements - upper triangle
        DO row = 1, col - 1
           Uinv(row, col) = rdivrp1
        END DO
        ! diagonal
        Uinv(col, col) = 2.0 * rdivrp1
        ! non-diagonal elements - lower triangle
        DO row = col + 1, rank
           Uinv(row, col) = rdivrp1
        END DO
     END DO
     
  END IF ttype

! ********************
! *** FINISHING UP ***
! ********************

END SUBROUTINE PDAF_seik_Uinv

!$Id$
!BOP
!
! !ROUTINE: PDAF_seik_TtimesA() --- Operate matrix T on some matrix
!
! !INTERFACE:
SUBROUTINE PDAF_seik_TtimesA(rank, dim_col, A, B)

! !DESCRIPTION:
! Operate matrix T on another matrix as
!                 B = T A\\
! \\
! T is a dim_ens x (dim_ens-1) matrix with zero column sums.
! There are two proposed forms of T (ensemble size N):\\
! typeT=0: diag(T)=1-1/N; nondiag(T)=-1/N; 
!          last row= -1/N\\
! typeT=1: diag(T)=1; nondiag(T)=0; last row = -1\\
!
! !  This is a core routine of PDAF and 
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2002-01 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
!  USE PDAF_memcounting, &
!       ONLY: PDAF_memcount
  use mod_prec
  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: rank               ! Rank of initial covariance matrix
  INTEGER, INTENT(in) :: dim_col            ! Number of columns in A and B
  REAL(dp), INTENT(in)    :: A(rank, dim_col)   ! Input matrix
  REAL(dp), INTENT(out)   :: B(rank+1, dim_col) ! Output matrix (TA)

! !CALLING SEQUENCE:
! Called by: User-provided prepoststep routines for SEIK and LSEIK
! Called by: PDAF_seik_analysis_newT
! Called by: PDAF_seik_resample_newT
! Called by: PDAF_lseik_analysis
! Called by: PDAF_lseik_resample
! Calls: PDAF_memcount
!EOP
  
! *** local variables ***
  INTEGER :: row, col  ! counters
  INTEGER :: typeT = 0 ! Which type of T
  REAL(dp) :: invdimens    ! Inversize of ensemble size
  INTEGER, SAVE :: allocflag = 0  ! Flag for dynamic allocation
  REAL(dp), ALLOCATABLE :: colmean(:) ! Mean values of columns of A


! **********************
! *** INITIALIZATION ***
! **********************

  ALLOCATE(colmean(dim_col))
  IF (allocflag == 0) THEN
     ! count allocated memory
 !    CALL PDAF_memcount(3, 'r', dim_col)
     allocflag = 1
  END IF
  colmean = 0.0
  invdimens = -1.0 / REAL(rank + 1)

  whichT: IF (typeT == 0) THEN

! *** Compute column means of A ***
    DO col = 1, dim_col
       DO row = 1, rank
          colmean(col) = colmean(col) + invdimens * A(row, col)
       END DO
    END DO

  END IF whichT


! ****************************************************
! ***  Operate T on A                              ***
! ***                                              ***
! *** Tv_1 = (v_11-mean(v_1), ... ,v_r1-mean(v_1)) ***
! *** with v_1 = (v_11,v_21, ... ,v_N1  )          ***
! ****************************************************

  ! first DIM rows
  DO col = 1, dim_col
     DO row = 1, rank
        B(row, col) = A(row, col) + colmean(col)
     END DO
  END DO

  ! row RANK+1
  DO col = 1, dim_col
     B(row, col) = colmean(col)
  END DO


! ********************
! *** FINISHING UP ***
! ********************

  DEALLOCATE(colmean)

END SUBROUTINE PDAF_seik_TtimesA
!$Id$
!BOP
!
! !ROUTINE: PDAF_seik_omega - Generate random matrix with special properties
!
! !INTERFACE:
SUBROUTINE PDAF_seik_omega(rank, omega, omegatype)

! !DESCRIPTION:
! Generate a transformation matrix OMEGA for
! the generation and transformation of the 
! ensemble in the SEIK and LSEIK filter.
! Generated is a uniform orthogonal matrix OMEGA
! with R columns orthonormal in $R^{r+1}$
! and orthogonal to (1,...,1)' by iteratively 
! applying the Householder matrix onto random 
! vectors distributed uniformly on the unit sphere.
!
! This version initializes at each iteration step
! the whole Householder matrix and subsequently
! computes Omega using DGEMM from BLAS. All fields are 
! allocated once at their maximum required size.
! (On SGI O2K this is about a factor of 2.5 faster
! than the version applying BLAS DDOT, but requires
! more memory.)
!
! For omegatype/=1 a deterministic omega is computed
! where the Housholder matrix of (1,...,1)' is operated
! on an identity matrix.
!
! !  This is a core routine of PDAF and 
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2002-01 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
!  USE PDAF_memcounting, &
!       ONLY: PDAF_memcount
!  USE PDAF_mod_filtermpi, &
!       ONLY: mype
  use mod_prec
  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: rank      ! Approximated rank of covar matrix
  REAL(dp), INTENT(inout) :: omega(rank+1, rank) ! Matrix Omega
  INTEGER, INTENT(in) :: omegatype ! Select type of omega:
                                   !     (1) generated from random vectors
                                   ! (other) generated from deterministic vectors

! !CALLING SEQUENCE:
! Called by: Used-provided ensemble initialization routine for SEIK/LSEIK
! Called by: PDAF_seik_resample
! Called by: PDAF_lseik_resample
! Calls: dlarnv (LAPACK)
! Calls: dgemm  (BLAS)
!EOP

!  *** local variables ***
  INTEGER :: iter, col, row        ! counters
  INTEGER :: iseed(4)              ! seed array for random number routine
  REAL(dp) :: norm                     ! norm of random vector
  INTEGER :: pflag                 ! pointer flag
  INTEGER, SAVE :: first = 1       ! flag for init of random number seed
  REAL(dp) :: rndval                   ! temporary value for init of Householder matrix
  REAL(dp) :: rndnum                   ! Value of randum entry
  INTEGER, SAVE :: allocflag = 0   ! Flag for dynamic allocation
  REAL(dp), ALLOCATABLE :: rndvec(:)   ! vector of random numbers
  REAL(dp), ALLOCATABLE :: house(:,:)  ! Row of the Householder matrix
  REAL(dp), POINTER :: omega_iter(:,:)         ! Pointer to temporary Omega field
  REAL(dp), POINTER :: omega_itermin1(:,:)     ! Pointer to temporary Omega field
  REAL(dp), ALLOCATABLE, TARGET :: temp1(:,:)  ! fields holding temporary Omega
  REAL(dp), ALLOCATABLE, TARGET :: temp2(:,:)  ! fields holding temporary Omega


! **********************
! *** INITIALIZATION ***
! **********************

  randomega: IF (omegatype == 1) THEN
     ! *** Generate omega by random vectors ***

!     if (mype==0) write (*,'(9x,a)') '--- Compute random Omega'

     ! allocate fields
     ALLOCATE(rndvec(rank + 1))
     ALLOCATE(house(rank + 1, rank))
     ALLOCATE(temp1(rank, rank), temp2(rank, rank))
!     IF (allocflag == 0) THEN
!        ! count allocated memory
!        CALL PDAF_memcount(4, 'r', (rank + 1) + (rank + 1) * rank + 2 * rank**2)
!        allocflag = 1
!     END IF

     ! set pointers
     omega_itermin1 => temp1
     omega_iter     => temp2
     pflag = 0

     ! Initialized seed for random number routine
     IF (first == 1) THEN
        iseed(1) = 1000
        iseed(2) = 2034
        iseed(3) = 0
        iseed(4) = 3
        first = 2
     END IF


! ***************************************
! *** First step of iteration         ***  
! *** Determine omega_iter for iter=1 ***
! ***************************************

     ! Get random number [-1,1]
     print *, 'get rand number seik_gen_omega'
     CALL dlarnv(2, iseed, 1, rndvec(1))
     print *, 'rand number is:', rndvec(1)
  
     IF (rndvec(1) >= 0.0) THEN
        omega_itermin1(1, 1) = +1.0
     ELSE
        omega_itermin1(1, 1) = -1.0
     END IF

! *****************
! *** Iteration ***
! *****************

     iteration: DO iter = 2, rank

! *** Initialize new random vector ***
      
        ! Get random vector of dimension DIM (elements in [-1,1])
        CALL dlarnv(2, iseed, iter, rndvec(1:iter))

        ! Normalize random vector
        norm = 0.0
        DO col = 1, iter
           norm = norm + rndvec(col)**2
        END DO
        norm = SQRT(norm)
        
        DO col = 1, iter
           rndvec(col) = rndvec(col) / norm
        END DO

! *** Compute Householder matrix ***

        ! First ITER-1 rows
        rndval = 1.0 / (ABS(rndvec(iter)) + 1.0)
        housecol: DO col = 1, iter - 1
           houserow: DO row = 1,iter - 1
              house(row, col) = - rndvec(row) * rndvec(col) * rndval
           END DO houserow
        END DO housecol
        
        DO col = 1, iter - 1
           house(col, col) = house(col, col) + 1.0
        END DO

        ! Last row
        housecol2: DO col = 1, iter - 1
           house(iter, col) = - (rndvec(iter) + SIGN(1.0_DP, rndvec(iter))) &
                * rndvec(col) * rndval
        END DO housecol2

! *** Compute omega on this iteration stage ***
        
        ! First iter-1 columns
        CALL dgemm ('n', 'n', iter, iter - 1, iter - 1, &
             one_db, house, rank + 1, omega_itermin1, rank, &
             zero_db, omega_iter, rank)

        ! Final column
        DO row = 1, iter
           omega_iter(row, iter) = rndvec(row)
        END DO

! *** Adjust pointers to temporal OMEGA fields ***

        IF (pflag == 0) THEN
           omega_itermin1 => temp2
           omega_iter     => temp1
           pflag = 1
        ELSE IF (pflag == 1) THEN
           omega_itermin1 => temp1
           omega_iter     => temp2
           pflag = 0
        END IF

     END DO iteration


! ********************************************
! ***            Final step                ***
! *** Projecting orthogonal to (1,...,1)^T ***
! ********************************************
    
     rndvec(1) = 1.0 / SQRT(REAL(rank + 1))

! *** Compute Householder matrix ***

     ! First r rows
     rndval = - rndvec(1) * rndvec(1) / (rndvec(1) + 1.0)
     housecolb: DO col = 1, rank
        houserowb: DO row = 1, rank
           house(row, col) = rndval
        END DO houserowb
     END DO housecolb
     
     DO col = 1, rank
        house(col, col) = house(col, col) + 1.0
     END DO
     
     ! Last row
     rndval = - (rndvec(1) + 1.0) * rndvec(1) / (rndvec(1) + 1.0)
     housecolc: DO col = 1, rank
        house(rank + 1, col) = rndval
     END DO housecolc

! *** Compute omega ***

     ! First iter-1 columns
     CALL dgemm ('n', 'n', rank + 1, rank, rank, &
          one_db, house, rank + 1, omega_itermin1, rank, &
          zero_db, omega, rank + 1)

! *** CLEAN UP ***

     NULLIFY(omega_itermin1, omega_iter)
     DEALLOCATE(temp1, temp2)
     DEALLOCATE(rndvec, house)

  ELSE randomega
     ! *** Generate Omega by deterministic vectors ***
     ! *** given by last step of upper recursion   ***
     ! *** with omega_itermin1 being the identity  ***

!     if (mype==0) write (*,'(8x,a)') '--- Compute deterministic Omega'

     rndnum = 1.0 / SQRT(REAL(rank + 1))

! *** Compute Householder matrix ***

     ! First r rows
     rndval = - rndnum * rndnum / (rndnum + 1.0)
     omegacolb: DO col = 1, rank
        omegarowb: DO row = 1, rank
           omega(row, col) = rndval
        END DO omegarowb
     END DO omegacolb
     
     DO col = 1, rank
        omega(col, col) = omega(col, col) + 1.0
     END DO

     ! Last row
     rndval = - (rndnum + 1.0) * rndnum / (rndnum + 1.0)
     omegacolc: DO col = 1, rank
        omega(rank + 1, col) = rndval
     END DO omegacolc

  END IF randomega

END SUBROUTINE PDAF_seik_omega


# endif
END MODULE MOD_ENKFA
