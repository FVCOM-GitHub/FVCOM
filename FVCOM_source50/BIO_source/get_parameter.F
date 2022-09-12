   MODULE MOD_PARAMETER
!** This routine reads in the control parameters for the biological model.  
      USE MOD_1D
      USE MOD_PHYTOPLANKTON
      USE MOD_ZOOPLANKTON
      USE MOD_BACTERIA
      USE MOD_DETRITUS
      USE MOD_DOM
      USE MOD_NUTRIENT
      Use Input_Util2
!      USE LIMS, ONLY: KB
      IMPLICIT NONE
      SAVE
      CHARACTER(LEN=50)              ::  BIONAME,            PARMINDEX
      CHARACTER(LEN=50), ALLOCATABLE ::  PHYTO_NAME(:),      PHYTO_UNIT(:)
      CHARACTER(LEN=50), ALLOCATABLE ::  ZOO_NAME(:),        ZOO_UNIT(:)
      CHARACTER(LEN=50), ALLOCATABLE ::  NUTRIENT_NAME(:),   NUTRIENT_UNIT(:)
      CHARACTER(LEN=50), ALLOCATABLE ::  DETRITUS_NAME(:),   DETRITUS_UNIT(:)
      CHARACTER(LEN=50), ALLOCATABLE ::  DOM_NAME(:),        DOM_UNIT(:)
      CHARACTER(LEN=50), ALLOCATABLE ::  BACTERIA_NAME(:),   BACTERIA_UNIT(:)
      CHARACTER(LEN=50), ALLOCATABLE ::  AUX_NAME(:),        AUX_UNIT(:)
      CHARACTER(LEN=50)              ::  MODEL_STRUCTURE     
      REAL(SPP)                      ::  TOTALL


      CHARACTER(LEN=250)             ::  biofile 

!*****************************  CONTAINS   ********************************
      CONTAINS
      SUBROUTINE GET_PARAMETER(fname)
      character(len=*)    :: fname
      logical             :: fexist
      INTEGER                        ::  I,J,K

!************************    BEGINNING EXECUTABLE    **********************
      KBVM1=KBV-1

      !ensure sediment parameter file exists    
      inquire(file=trim(fname),exist=fexist)
      if(.not.fexist)then
        write(*,*)'biological parameter file: ',trim(fname),' does not exist'
        write(*,*)'stopping'
        stop
      end if

      OPEN (1, FILE=trim(fname),status='old')
      REWIND (1)
      READ(1,*)PARMINDEX
      DO 10 WHILE (PARMINDEX .NE. 'STRUCTURE')             !DETERMING MODEL STRUCTURE
      READ(1,*)PARMINDEX
 10   CONTINUE
      READ(1,*)PARMINDEX
      MODEL_STRUCTURE=PARMINDEX
      SELECT CASE (PARMINDEX)
         CASE ('NPZ')
            NNN=1
            NNP=1
            NNZ=1
            NNM=0
            NNB=0
            NND=0
            NNA=0
            NTT=NNN+NNP+NNZ+NNM+NNB+NND+NNA   !TOTAL VARIABLE NUMBER FOR VERTICAL MIXING
            INN=1
            INP=INN+NNN
            INZ=INP+NNP
            INM=INZ+NNZ
            INB=INM+NNM
            IND=INB+NNB
            INA=IND+NND
            ALLOCATE (NUTRIENT_NAME(NNN))
            ALLOCATE (NUTRIENT_UNIT(NNN))
            ALLOCATE (PHYTO_NAME(NNP))
            ALLOCATE (PHYTO_UNIT(NNP))
            ALLOCATE (ZOO_NAME(NNZ))
            ALLOCATE (ZOO_UNIT(NNZ))
            NUTRIENT_NAME(1)   = 'Nitrogen'
	    NUTRIENT_UNIT(1)   = 'mmole N m-3'
            PHYTO_NAME(1)      = 'Phytoplankton'
	    PHYTO_UNIT(1)      = 'mmole C m-3'
            ZOO_NAME(1)        = 'Zooplankton'
	    ZOO_UNIT(1)        = 'mmole C m-3'

	  CASE ('NPZD')
            NNN=1
            NNP=1
            NNZ=1
            NNM=0
            NNB=0
            NND=1
            NNA=0
            NTT=NNN+NNP+NNZ+NNM+NNB+NND+NNA   !TOTAL VARIABLE NUMBER FOR VERTICAL MIXING
            INN=1
            INP=INN+NNN
            INZ=INP+NNP
            INM=INZ+NNZ
            INB=INM+NNM
            IND=INB+NNB
            INA=IND+NND
            ALLOCATE (NUTRIENT_NAME(NNN))
            ALLOCATE (NUTRIENT_UNIT(NNN))
            ALLOCATE (PHYTO_NAME(NNP))
            ALLOCATE (PHYTO_UNIT(NNP))
            ALLOCATE (ZOO_NAME(NNZ))
            ALLOCATE (ZOO_UNIT(NNZ))
            ALLOCATE (DETRITUS_NAME(NND))
            ALLOCATE (DETRITUS_UNIT(NND))
            NUTRIENT_NAME(1)   = 'Nitrogen'
	    NUTRIENT_UNIT(1)   = 'mmole N m-3'
            PHYTO_NAME(1)      = 'Phytoplankton'
	    PHYTO_UNIT(1)      = 'mmole C m-3'
            ZOO_NAME(1)        = 'Zooplankton'
	    ZOO_UNIT(1)        = 'mmole C m-3'
            DETRITUS_NAME(1)   = 'Detritus'
	    DETRITUS_UNIT(1)   = 'mmole C m-3'

	  CASE ('FASHAM')
            NNN=2
            NNP=1
            NNZ=1
            NNM=1
            NNB=1
            NND=1
            NNA=0
            NTT=NNN+NNP+NNZ+NNM+NNB+NND+NNA   !TOTAL VARIABLE NUMBER FOR VERTICAL MIXING
            INN=1
            INP=INN+NNN
            INZ=INP+NNP
            INM=INZ+NNZ
            INB=INM+NNM
            IND=INB+NNB
            INA=IND+NND
            ALLOCATE (NUTRIENT_NAME(NNN))
            ALLOCATE (NUTRIENT_UNIT(NNN))
            ALLOCATE (PHYTO_NAME(NNP))
            ALLOCATE (PHYTO_UNIT(NNP))
            ALLOCATE (ZOO_NAME(NNZ))
            ALLOCATE (ZOO_UNIT(NNZ))
            ALLOCATE (DETRITUS_NAME(NND))
            ALLOCATE (DETRITUS_UNIT(NND))
            ALLOCATE (DOM_NAME(NNM))
            ALLOCATE (DOM_UNIT(NNM))
            ALLOCATE (BACTERIA_NAME(NNB))
            ALLOCATE (BACTERIA_UNIT(NNB))
            NUTRIENT_NAME(1)   = 'NH4'
	    NUTRIENT_UNIT(1)   = 'mmole N m-3'
            NUTRIENT_NAME(2)   = 'NO3'
	    NUTRIENT_UNIT(2)   = 'mmole N m-3'
            PHYTO_NAME(1)      = 'Phytoplankton'
	    PHYTO_UNIT(1)      = 'mmole C m-3'
            ZOO_NAME(1)        = 'Zooplankton'
	    ZOO_UNIT(1)        = 'mmole C m-3'
            DETRITUS_NAME(1)   = 'Detritus'
	    DETRITUS_UNIT(1)   = 'mmole C m-3'
            DOM_NAME(1)        = 'DOM'
	    DOM_UNIT(1)        = 'mmole C m-3'
            BACTERIA_NAME(1)   = 'Bacteria'
	    BACTERIA_UNIT(1)   = 'mmole C m-3'

          CASE ('N2P2Z2D2')
            NNN=3
            NNP=2
            NNZ=2
            NNM=0
            NNB=0
            NND=2
            NNA=0
            NTT=NNN+NNP+NNZ+NNM+NNB+NND+NNA   !TOTAL VARIABLE NUMBER FOR VERTICAL MIXING
            INN=1
            INP=INN+NNN
            INZ=INP+NNP
            INM=INZ+NNZ
            INB=INM+NNM
            IND=INB+NNB
            INA=IND+NND
            ALLOCATE (NUTRIENT_NAME(NNN))
            ALLOCATE (NUTRIENT_UNIT(NNN))
            ALLOCATE (PHYTO_NAME(NNP))
            ALLOCATE (PHYTO_UNIT(NNP))
            ALLOCATE (ZOO_NAME(NNZ))
            ALLOCATE (ZOO_UNIT(NNZ))
            ALLOCATE (DETRITUS_NAME(NND))
            ALLOCATE (DETRITUS_UNIT(NND))
            NUTRIENT_NAME(1)   = 'NH4'
	    NUTRIENT_UNIT(1)   = 'mmole N m-3'
            NUTRIENT_NAME(2)   = 'NO3'
	    NUTRIENT_UNIT(2)   = 'mmole N m-3'
            NUTRIENT_NAME(3)   = 'SiO3'
	    NUTRIENT_UNIT(3)   = 'mmole Si m-3'
            PHYTO_NAME(1)      = 'Small_phyto'
	    PHYTO_UNIT(1)      = 'mmole C m-3'
            PHYTO_NAME(2)      = 'Large_phyto'
	    PHYTO_UNIT(2)      = 'mmole C m-3'
            ZOO_NAME(1)        = 'Microzooplankton'
	    ZOO_UNIT(1)        = 'mmole C m-3'
            ZOO_NAME(2)        = 'Mesozooplankton'
	    ZOO_UNIT(2)        = 'mmole C m-3'
            DETRITUS_NAME(1)   = 'Detritus_N'
	    DETRITUS_UNIT(1)   = 'mmole N m-3'
            DETRITUS_NAME(2)   = 'Detritus_Si'
	    DETRITUS_UNIT(2)   = 'mmole Si m-3'

          CASE ('N2P2Z2D2B')
            NNN=2
            NNP=2
            NNZ=2
            NNM=0
            NNB=1
            NND=2
            NNA=0
            NTT=NNN+NNP+NNZ+NNM+NNB+NND+NNA   !TOTAL VARIABLE NUMBER FOR VERTICAL MIXING
            INN=1
            INP=INN+NNN
            INZ=INP+NNP
            INM=INZ+NNZ
            INB=INM+NNM
            IND=INB+NNB
            INA=IND+NND
            ALLOCATE (NUTRIENT_NAME(NNN))
            ALLOCATE (NUTRIENT_UNIT(NNN))
            ALLOCATE (PHYTO_NAME(NNP))
            ALLOCATE (PHYTO_UNIT(NNP))
            ALLOCATE (ZOO_NAME(NNZ))
            ALLOCATE (ZOO_UNIT(NNZ))
            ALLOCATE (DETRITUS_NAME(NND))
            ALLOCATE (DETRITUS_UNIT(NND))
            ALLOCATE (BACTERIA_NAME(NNB))
            ALLOCATE (BACTERIA_UNIT(NNB))
            NUTRIENT_NAME(1)   = 'PO4'
	    NUTRIENT_UNIT(1)   = 'mmole P m-3'
            NUTRIENT_NAME(2)   = 'SiO3'
	    NUTRIENT_UNIT(2)   = 'mmole Si m-3'
            PHYTO_NAME(1)      = 'Small_phyto'
	    PHYTO_UNIT(1)      = 'mmole C m-3'
            PHYTO_NAME(2)      = 'Large_phyto'
	    PHYTO_UNIT(2)      = 'mmole C m-3'
            ZOO_NAME(1)        = 'Microzooplankton'
	    ZOO_UNIT(1)        = 'mmole C m-3'
            ZOO_NAME(2)        = 'Mesozooplankton'
	    ZOO_UNIT(2)        = 'mmole C m-3'
            DETRITUS_NAME(1)   = 'Detritus_P'
	    DETRITUS_UNIT(1)   = 'mmole P m-3'
            DETRITUS_NAME(2)   = 'Detritus_Si'
	    DETRITUS_UNIT(2)   = 'mmole Si m-3'
            BACTERIA_NAME(1)   = 'Bacteria'
	    BACTERIA_UNIT(1)   = 'mmole C m-3'
	            
          CASE ('N2P2Z2D2SB')
            NNN=2
            NNP=2
            NNZ=2
            NNM=1
            NNB=2
            NND=2
            NNA=0
            NTT=NNN+NNP+NNZ+NNM+NNB+NND+NNA   !TOTAL VARIABLE NUMBER FOR VERTICAL MIXING
            INN=1
            INP=INN+NNN
            INZ=INP+NNP
            INM=INZ+NNZ
            INB=INM+NNM
            IND=INB+NNB
            INA=IND+NND
            ALLOCATE (NUTRIENT_NAME(NNN))
            ALLOCATE (NUTRIENT_UNIT(NNN))
            ALLOCATE (PHYTO_NAME(NNP))
            ALLOCATE (PHYTO_UNIT(NNP))
            ALLOCATE (ZOO_NAME(NNZ))
            ALLOCATE (ZOO_UNIT(NNZ))
            ALLOCATE (DETRITUS_NAME(NND))
            ALLOCATE (DETRITUS_UNIT(NND))
            ALLOCATE (DOM_NAME(NNM))
            ALLOCATE (DOM_UNIT(NNM))
            ALLOCATE (BACTERIA_NAME(NNB))
            ALLOCATE (BACTERIA_UNIT(NNB))
            NUTRIENT_NAME(1)   = 'NH4'
	    NUTRIENT_UNIT(1)   = 'mmole N m-3'
            NUTRIENT_NAME(2)   = 'NO3'
	    NUTRIENT_UNIT(2)   = 'mmole N m-3'
            PHYTO_NAME(1)      = 'Small_phyto'
	    PHYTO_UNIT(1)      = 'mmole C m-3'
            PHYTO_NAME(2)      = 'Large_phyto'
	    PHYTO_UNIT(2)      = 'mmole C m-3'
            ZOO_NAME(1)        = 'Microzooplankton'
	    ZOO_UNIT(1)        = 'mmole C m-3'
            ZOO_NAME(2)        = 'Mesozooplankton'
	    ZOO_UNIT(2)        = 'mmole C m-3'
            DETRITUS_NAME(1)   = 'Small_detritus'
	    DETRITUS_UNIT(1)   = 'mmole C m-3'
            DETRITUS_NAME(2)   = 'Large_detritus'
	    DETRITUS_UNIT(2)   = 'mmole C m-3'
            DOM_NAME(1)        = 'DOM'
	    DOM_UNIT(1)        = 'mmole C m-3'
            BACTERIA_NAME(1)   = 'Bacteria'
	    BACTERIA_UNIT(1)   = 'mmole C m-3'

         CASE('NEW')
           READ(1,*)PARMINDEX
           DO 20 WHILE (PARMINDEX.NE.'L_FUNCTION')     
             SELECT CASE (PARMINDEX)                        !ASSIGNMENT FOR NEW MODEL
               CASE ('PHYTOPLANKTON')
                 READ(1,*)NNP
                 ALLOCATE (PHYTO_NAME(NNP), PHYTO_UNIT(NNP))
                 DO J=1,NNP
                   READ(1,*)PHYTO_NAME(J),PHYTO_UNIT(J)
                 END DO

               CASE ('ZOOPLANKTON')
                 READ(1,*)NNZ
                 ALLOCATE (ZOO_NAME(NNZ), ZOO_UNIT(NNZ))
                 DO J=1,NNZ
                   READ(1,*)ZOO_NAME(J),ZOO_UNIT(J)
                 END DO

               CASE ('NUTRIENT')
                 READ(1,*)NNN
                 ALLOCATE (NUTRIENT_NAME(NNN), NUTRIENT_UNIT(NNN))
                 DO J=1,NNN
                   READ(1,*)NUTRIENT_NAME(J),NUTRIENT_UNIT(J)
                 END DO

               CASE ('BACTERIA')
                 READ(1,*)NNB
                 ALLOCATE (BACTERIA_NAME(NNB), BACTERIA_UNIT(NNB))
                 DO J=1,NNB
                   READ(1,*)BACTERIA_NAME(J),BACTERIA_UNIT(J)
                 END DO

               CASE ('DETRITUS')
                 READ(1,*)NND
                 ALLOCATE (DETRITUS_NAME(NND), DETRITUS_UNIT(NND))
                 DO J=1,NND
                   READ(1,*)DETRITUS_NAME(J),DETRITUS_UNIT(J)
                 END DO

               CASE ('DOM')
                 READ(1,*)NNM
                 ALLOCATE (DOM_NAME(NNM), DOM_UNIT(NNM))
                 DO J=1,NNM
                   READ(1,*)DOM_NAME(J),DOM_UNIT(J)
                 END DO

               CASE ('AUXILIARY')
                 READ(1,*)NNA
                 ALLOCATE (AUX_NAME(NNA), AUX_UNIT(NNA))
                 DO J=1,NNA
                   READ(1,*)AUX_NAME(J),AUX_UNIT(J)
                 END DO
             END SELECT
             READ(1,*)PARMINDEX
 20       CONTINUE
            NTT=NNN+NNP+NNZ+NNM+NNB+NND+NNA   !TOTAL VARIABLE NUMBER FOR VERTICAL MIXING
            INN=1
            INP=INN+NNN
            INZ=INP+NNP
            INM=INZ+NNZ
            INB=INM+NNM
            IND=INB+NNB
            INA=IND+NND
      END SELECT

!*******************    BIO_NAME ASSIGNMENT FOR NETCDF OUT PUT   ***************
            ALLOCATE (BIO_NAME(NTT,4))
            BIO_NAME(:,3)  = 'Concentration'
            BIO_NAME(:,4)  = 'Array'
            DO J=INN,NNN
              BIO_NAME(J,1)=NUTRIENT_NAME(J)
              BIO_NAME(J,2)=NUTRIENT_UNIT(J)
            END DO
            DO J=INP,INP+NNP-1
              BIO_NAME(J,1)=PHYTO_NAME(J-INP+1)
              BIO_NAME(J,2)=PHYTO_UNIT(J-INP+1)
            END DO
            DO J=INZ,INZ+NNZ-1
              BIO_NAME(J,1)=ZOO_NAME(J-INZ+1)
              BIO_NAME(J,2)=ZOO_UNIT(J-INZ+1)
            END DO
            DO J=INM,INM+NNM-1
              BIO_NAME(J,1)=DOM_NAME(J-INM+1)
              BIO_NAME(J,2)=DOM_UNIT(J-INM+1)
            END DO
            DO J=INB,INB+NNB-1
              BIO_NAME(J,1)=BACTERIA_NAME(J-INB+1)
              BIO_NAME(J,2)=BACTERIA_UNIT(J-INB+1)
            END DO
            DO J=IND,IND+NND-1
              BIO_NAME(J,1)=DETRITUS_NAME(J-IND+1)
              BIO_NAME(J,2)=DETRITUS_UNIT(J-IND+1)
            END DO
            DO J=INA,INA+NNA-1
              BIO_NAME(J,1)=AUX_NAME(J-INA+1)
              BIO_NAME(J,2)=AUX_UNIT(J-INA+1)
            END DO
            NO3_ON = .FALSE.                                       !NH4 INHIBITION
            PO4_ON = .FALSE.                                       !BAC UPTAKE ON PO4
            DOCDON = .FALSE.                                       !BOTH CONSIDERED
            IF(NNN.GT.1) THEN
              IF(NUTRIENT_NAME(2).EQ.'NO3') NO3_ON = .TRUE.
            ENDIF
            IF(NNM.GE.2) THEN
              IF(DOM_NAME(1).EQ.'DOC'.AND.DOM_NAME(2).EQ.'DON') DOCDON = .TRUE.
            END IF
	    IF(NUTRIENT_NAME(1).EQ.'PO4') PO4_ON = .TRUE.

!************* ALLOCATING                   GLOBAL PARAMETERS      AND INITIALIZAITON
      ALLOCATE (BIO_VAR(KBV,NTT))        ;  BIO_VAR     =  0.0_SPP  !ENSEMBLE VARIABLES FOR MIXING
      ALLOCATE (DELTA_D(KBV))            ;  DELTA_D     =  1.0_SPP  !LAYER THICKNESS
      ALLOCATE (DELTA_Z(KBV))            ;  DELTA_Z     =  1.0_SPP  !THICKNESS BETWEEN LAYER CENTRE
      ALLOCATE (DEPTH_Z(KBV))            ;  DEPTH_Z     =  1.0_SPP  !DEPTH AT EACH LAYER CENTER
      ALLOCATE (KM_BIO(KBV))             ;  KM_BIO      =  0.0_SPP  !DIFFUSIVITY 
      ALLOCATE (T_BIO(KBV))              ;  T_BIO       =  0.0_SPP  !TEMPERATURE

!************* ALLOCATING                   PHYTOPLANKTON  AND     INITIALIZAITON *************
      ALLOCATE (ALPHAP(NNP))             ;  ALPHAP      =  0.0_SPP  !LIGHT FUNCTION COEFFICIENT
      ALLOCATE (ALPHA_U(NNP))            ;  ALPHA_U     =  0.5_SPP  !COEFFICIENT BETWEEN MULTIPLICAITON & MINIMUM
      ALLOCATE (A_TP(NNP))               ;  A_TP        =  0.0_SPP  !EXPONENTIAL FOR T FORCING
                                         ;  ATANU_C     =  0.03_SPP !ATTENUATION COEFFICIENT OF CHLOROPHYLL
                                         ;  ATANU_D     =  0.01_SPP !ATTENUATION COEFFICIENT OF DETRITUS
                                         ;  ATANU_W     =  0.04_SPP !ATTENUATION COEFFICIENT OF WATER
      ALLOCATE (BETAP(NNP))              ;  BETAP       =  0.0_SPP  !PHOTOINHIBITION COEFFICIENT
      ALLOCATE (BIO_P(KBV,NNP))          ;  BIO_P       =  0.0_SPP  !PHYTO BIOMASS
      ALLOCATE (CHL2C(NNP))              ;  CHL2C       =  0.24_SPP !CHLOROPHYLL TO CARBON RATIO IN PHYTO
      ALLOCATE (D_DOM(NNP))              ;  D_DOM       =  0.0_SPP  !ACTIVE EXUDATION COEFFICIENT
      ALLOCATE (DPDOM(NNP))              ;  DPDOM       =  0.0_SPP  !COEFFICIENT OF PHYTO EXUDATION
      ALLOCATE (FVP(KBV,NNP))            ;  FVP         =  0.0_DPP  !PHYTO AFTER THRESHOLD SUBSTRACTION
      ALLOCATE (I_OPT(NNP))              ;  I_OPT       =  100._SPP !OPTIMAL LIGHT FOR PHYTO
                                         ;  IRRAD0      =  200._SPP !IRRADIANCE AT THE SEA SURFACE
      ALLOCATE (K_LIGHT(NNP))            ;  K_LIGHT     =  50._SPP  !HALF_SATURATION OF LIGHT FUNCTIONS
      ALLOCATE (MPD(NNP))                ;  MPD         =  0.0_SPP  !PHYTO MORTALITY EOFFICIENT
      ALLOCATE (M_P(NNP))                ;  M_P         =  1.0_SPP  !PHYTO MORTALITY ORDER
      ALLOCATE (N_P(NNP))                ;  N_P         =  1.0_SPP  !POWER OF LIGHT FUNCTIONS
                                         ;  PARFRAC     =  0.43_SPP
      ALLOCATE (P_D(KBV,NNP))            ;  P_D         =  0.0_SPP  !PHYTO MORTALITY
      ALLOCATE (P_DOM(KBV,NNP))          ;  P_DOM       =  0.0_SPP  !PHYTO EXUDATION OF DOM
      ALLOCATE (P_N(KBV,NNP))            ;  P_N         =  0.0_SPP  !PHYTO RESPIRATION
      ALLOCATE (P_0(NNP))                ;  P_0         =  0.0_SPP  !PHYTO THRESHOLD
                                         ;  RP_T        =  0.0_SPP  !EXPONENTIAL FOR T FORCING ON RESPIRATION
      ALLOCATE (R_P(NNP))                ;  R_P         =  1.0E-7_SPP  !PHYTO RESPIRATION COEFFICIENT
      ALLOCATE (T_OPTP(NNP))             ;  T_OPTP      =  20._SPP  !OPTIMAL TEMPERATURE FOR PHYTO
      ALLOCATE (ULR(NNP))                ;  ULR         =  0.0_SPP  !LIGHT GROWHT FACTOR FOR PHYTO
      ALLOCATE (UMAX(NNP))               ;  UMAX        =  1.16E-5_SPP !PHYTO MAXIMUM GROWTH RATE
      ALLOCATE (UNR(NNP,NNN+1))          ;  UNR         =  1.0_SPP  !NUTRIENT GROWTH FACTOR FOR PHYTO
      ALLOCATE (UNRMIN(NNP))             ;  UNRMIN      =  1.0_SPP  !MINIMUM NUTRIENT FACTOR FOR PHYTO
      ALLOCATE (U_P(KBV,NNP))            ;  U_P         =  0.0_SPP  !PHYTO GROWTH
      ALLOCATE (UTR(NNP))                ;  UTR         =  1.0_SPP  !T GORWTH FACTOR FOR PHYTO
      ALLOCATE (W_P(NNP))                ;  W_P         =  1.0_DPP  !SINKING VELOCITY OF PHYTO
      ALLOCATE (WSNK_P(KBV,NNP))         ;  WSNK_P      =  0.0_DPP  !SINKING FLUX OF PHYTO

!**************ALLOCATING                   ZOOPLANKTON    AND     INITIALIZATION *****************
      ALLOCATE (ACTIVE_R(NNZ))           ;  ACTIVE_R    =  0.0_SPP  !ZOO ACTIVE RESPIRATION COEFFICIENT
      ALLOCATE (A_TZ(NNZ))               ;  A_TZ        =  0.0_SPP  !T EXPONENTAIL COEFFICIENT
      ALLOCATE (BIO_Z(KBV,NNZ))          ;  BIO_Z       =  0.0_SPP  !ZOOPLANKTON BIOMASS
      ALLOCATE (EFFIB(NNB,NNZ))          ;  EFFIB       =  0.3_SPP  !GROWTH EFFICIENCY ON BACTERIA
      ALLOCATE (EFFID(NND,NNZ))          ;  EFFID       =  0.3_SPP  !GROWTH EFFICIENCY ON DETRITUS
      ALLOCATE (EFFIP(NNP,NNZ))          ;  EFFIP       =  0.3_SPP  !GROWTH EFFICIENCY ON PHYTO
      ALLOCATE (EFFIZ(NNZ,NNZ))          ;  EFFIZ       =  0.3_SPP  !GROWTH EFFICIENCY ON ZOOPLANKTON
      ALLOCATE (FVZ(KBV,NNP))            ;  FVZ         =  0.0_SPP  !ZOOPLANKTON SUBSTRACTED OF THRESHOLD
      ALLOCATE (G_B(KBV,NNB,NNZ))        ;  G_B         =  0.0_SPP  !PREDATION ON BACTERIA
      ALLOCATE (G_D(KBV,NND,NNZ))        ;  G_D         =  0.0_SPP  !GRAZING ON DETRITUS
      ALLOCATE (G_MAX(NNZ))              ;  G_MAX       =  1.16E-5_SPP !MAXIMUM GRAZING RATE
      ALLOCATE (G_P(KBV,NNP,NNZ))        ;  G_P         =  0.0_SPP  !GRAZING ON PHYTOPLANKTON
      ALLOCATE (G_RATE(NNZ))             ;  G_RATE      =  0.0_SPP  !GRAZING RATE
      ALLOCATE (G_Z(KBV,NNZ,NNZ))        ;  G_Z         =  0.0_SPP  !PREDATION OF ZOOPLANKTON
      ALLOCATE (K_ZG(NNZ))               ;  K_ZG        =  1.0_SPP  !HALF-SATURATION FOR ZOO GRAZING
      ALLOCATE (M_G(NNZ))                ;  M_G         =  1.0_SPP  !ZOOPLANKTON GRAZING POWER
      ALLOCATE (MZD(NNZ))                ;  MZD         =  0.0_SPP  !ZOOPLANKTON MORTALITY COEFFICIENT
      ALLOCATE (M_Z(NNZ))                ;  M_Z         =  1.0_SPP  !ZOOPLANKTON MORTALITY POWER
      ALLOCATE (P_C(NNZ))                ;  P_C         =  0.0_SPP  !GRAZING THRESHOLD
      ALLOCATE (RECRUIT_LOSS(0:NNZ))     ;  RECRUIT_LOSS=  0.0_SPP  !RECRUITMENT LOSS
      ALLOCATE (R_RECRUIT(NNZ))          ;  R_RECRUIT   =  0.0_SPP  !RECRUIT COEFFICIENT
      ALLOCATE (R_Z(NNZ))                ;  R_Z         =  1.E-7_SPP   !ZOOPLANKTON RESPIRATION COEFFICIENT
      ALLOCATE (SIGMA_B(NNB,NNZ))        ;  SIGMA_B     =  1.0_SPP  !PREFERENCE ON BACTERIA
      ALLOCATE (SIGMA_D(NND,NNZ))        ;  SIGMA_D     =  1.0_SPP  !PREFERENCE ON DETRITUS
      ALLOCATE (SIGMA_P(NNP,NNZ))        ;  SIGMA_P     =  1.0_SPP  !PREFERENCE ON PHYTO
      ALLOCATE (SIGMA_Z(NNZ,NNZ))        ;  SIGMA_Z     =  1.0_SPP  !PREFERENCE ON ZOOPLANKTON
      ALLOCATE (T_OPTZ(NNZ))             ;  T_OPTZ      =  20._SPP  !OPTIMAL TEMPERATURE FOR ZOOPLANKTON
      ALLOCATE (Z_D(KBV,NNZ))            ;  Z_D         =  0.0_SPP  !ZOOPANKTON MORTALITY
      ALLOCATE (Z_N(KBV,NNZ))            ;  Z_N         =  0.0_SPP  !ZOOPLANKTON RESPIRATION
      ALLOCATE (Z_0(NNZ))                ;  Z_0         =  0.0_SPP  !ZOOPLANKTO THRESHOLLD
!************  ALLOCATING                   BACTERIA       AND     INITIALIZATION   *************
      ALLOCATE (A_TB(NNB))               ;  A_TB        =  0.0_SPP  !T EXPONENTAIL COEFFICIENT
      ALLOCATE (BIO_B(KBV,NNB))          ;  BIO_B       =  0.0_SPP  !BACTERIA BIOMASS
      ALLOCATE (B_0(NNB))                ;  B_0         =  0.0_SPP  !BACTERIA THRESHOLD
      ALLOCATE (DELTA_B(NNB))            ;  DELTA_B     =  0.4_SPP  !RATIO BETWEEN NH4 AND DON UPTAKE
      ALLOCATE (EFFIBD(NND,NNB))         ;  EFFIBD      =  0.3_SPP  !GROWTH EFFICIENCY ON DETRITUS
      ALLOCATE (EFFIDOM(NNM,NNB))        ;  EFFIDOM     =  0.3_SPP  !GROWTH EFFICIENCY ON DOM
      ALLOCATE (EFFIN(NNN,NNB))          ;  EFFIN       =  1.0_SPP  !GROWTH EFFICIENCY ON NH4
      ALLOCATE (FVB(KBV,NNB))            ;  FVB         =  0.0_SPP  !BACTERIA SUBSTRACTED OF THRESHOLD
      ALLOCATE (R_B(NNB))                ;  R_B         =  1.E-6_SPP   !BACTERIA RESPIRATION
      ALLOCATE (RESPIRATION_B(KBV,NNB))  ;  RESPIRATION_B= 0.0_SPP  !BACTERIA RESPIRATION
      ALLOCATE (SIGMA_BD(NND,NNB))       ;  SIGMA_BD    =  0.2_SPP  !PREFERENCE ON DETRITUS
      ALLOCATE (SIGMA_DOM(NNM,NNB))      ;  SIGMA_DOM   =  0.2_SPP  !PREFERENCE ON DOM
      ALLOCATE (SIGMA_N(NNN,NNB))        ;  SIGMA_N     =  1.0_SPP  !PREFERENCE ON DON
      ALLOCATE (T_OPTB(NNB))             ;  T_OPTB      =  20._SPP  !OPTIMAL TEMPERATURE FOR BACTERIA
      ALLOCATE (U_B(NNB))                ;  U_B         =  0.0_SPP  !BACTERIA GROWTH RATE
      ALLOCATE (UBMAX(NNB))              ;  UBMAX       =  1.E-5_SPP !MAXIMUM GROWTH RATE OF BACTERIA
      ALLOCATE (UPTAKE_BD(KBV,NND,NNB))  ;  UPTAKE_BD   =  0.0_SPP  !DETRITUS UPTAKE BY BACTERI
      ALLOCATE (UPTAKE_BN(KBV,NNN,NNB))  ;  UPTAKE_BN   =  0.0_SPP  !NUTRIENT UPTAKE BY BACTERIA
      ALLOCATE (UPTAKE_DOM(KBV,NNM,NNB)) ;  UPTAKE_DOM  =  0.0_SPP  !DOM UPTAKE BY BACTERIA
      ALLOCATE (UPTAKE_DON(KBV,NNB))     ;  UPTAKE_DON  =  0.0_SPP  !DON UPTAKE BY BACTERIA

!******************  ALLOCATING             DETRITUS       AND     INITIALIZATION   *************************
      ALLOCATE (AGGREGATION(KBV,0:NND)) ;AGGREGATION    =  0.0_SPP  !AGGREGATION
      ALLOCATE (ALPHA_BD(NND,NNB,NNZ))   ;  ALPHA_BD    =  0.0_SPP  !FEEDING LOSS FROM BACTERIA TO DETRITUS
      ALLOCATE (ALPHA_DAG(NND))          ;  ALPHA_DAG   =  0.0_SPP  !AGGREGATION COEFFICIENT
      ALLOCATE (ALPHA_DD(NND,NND,NNZ))   ;  ALPHA_DD    =  0.0_SPP !FEEDING LOSS FROM DETRITUS TO DETRITUS
      ALLOCATE (ALPHA_DDG(NND))          ;  ALPHA_DDG   =  0.0_SPP  !DISAGGREGATION COEFFICIENT
      ALLOCATE (ALPHA_PD(NND,NNP,NNZ))   ;  ALPHA_PD    =  0.0_SPP  !FEEDING LOSS FROM PHYTO TO DETRITUS
      ALLOCATE (ALPHA_ZD(NND,NNZ,NNZ))   ;  ALPHA_ZD    =  0.0_SPP  !FEEDING LOSS FROM ZOO TO DETRITUS
      ALLOCATE (BIO_D(KBV,NND))          ;  BIO_D       =  0.0_SPP  !DETRITUS ABUNDANCE
      ALLOCATE (DETRITE(KBV,NND))        ;  DETRITE     =  0.0_DPP  !DETRITUS CONCENTRATION 
      ALLOCATE (D_D(NND))                ;  D_D         =  0.0_SPP  !DISSOLUTION COEFFICIENT
      ALLOCATE (D_0(NND))                ;  D_0         =  0.0_SPP  !DETRITUS THRESHOLD
      ALLOCATE (D_RN(NND))               ;  D_RN        =  0.0_SPP  !REMINERALIZATION COEFFICIENT
      ALLOCATE (DISAGGREGATION(KBV,NND+1)); DISAGGREGATION =  0.0_SPP !DISAGGREGATION
      ALLOCATE (D_SINK(KBV,NND))         ;  D_SINK      =  0.0_DPP  !DETRITUS SINKING FLUX
      ALLOCATE (D_SINK1(KBV))            ;  D_SINK1     =  0.0_DPP  !SINKING FLUX WITHIN 1 ITERATION
      ALLOCATE (DISSOLUTION(KBV,NND))    ;  DISSOLUTION =  0.0_SPP  !DISSOLUTION
      ALLOCATE (EPSILON_PD(NND,NNP))     ;  EPSILON_PD  =  0.0_SPP  !PHYTO MORTALITY TO DETRITUS
      ALLOCATE (EPSILON_ZD(NND,NNZ))     ;  EPSILON_ZD  =  0.0_SPP  !ZOO MORTALTIY TO DETRITUS
      ALLOCATE (FVD(KBV,NND))            ;  FVD         =  0.0_SPP  !DETRITUS AFTER THRESHOLD SUBSTRACTION
      ALLOCATE (REMINERALIZATION(KBV,NND)); REMINERALIZATION=0.0_SPP !REMINERALIZATION   
      ALLOCATE (W_D(NND))                ;  W_D         =  0.0_SPP  !DETRITUS SINKING VELOCITY
      ALLOCATE (W_D_TOTAL(NND))          ;  W_D_TOTAL   =  0.0_SPP  !SINKING DISTANCE WITIN 1 T_STEP

!******************  ALLOCATING             DOM            AND     INITIALIZATION  *****************
      ALLOCATE (AGE_DOM(0:NNM))          ;  AGE_DOM     =  0.0_SPP  !DOM AGEING
      ALLOCATE (ALPHA_DOM(NNM))          ;  ALPHA_DOM   =  0.0_SPP  !DOM AGEING COEFFICIENT
      ALLOCATE (ALPHA_PDOM(NNM,NNP))     ;  ALPHA_PDOM  =  0.0_SPP  !PHYTO EXUDATION TO DOM
      ALLOCATE (ALPHA_DDOM(NNM,NND))     ;  ALPHA_DDOM  =  0.0_SPP  !DETRITUS DISSOLUTION TO DOM
      ALLOCATE (ALPHA_ZPDOM(NNM,NNP,NNZ));  ALPHA_ZPDOM =  0.0_SPP  !FEEDING LOSS FROM PHYTO TO DOM
      ALLOCATE (ALPHA_ZDDOM(NNM,NNP,NNZ));  ALPHA_ZDDOM =  0.0_SPP  !FEEDING LOSS FROM DETRITUS TO DOM
      ALLOCATE (ALPHA_ZZDOM(NNM,NNZ,NNZ));  ALPHA_ZZDOM =  0.0_SPP  !FEEDING LOSS FROM ZOO TO DOM
      ALLOCATE (ALPHA_ZBDOM(NNM,NNB,NNZ));  ALPHA_ZBDOM =  0.0_SPP  !FEEDING LOSS FROM BACTERIA TO DOM
      ALLOCATE (BIO_DOM(KBV,NNM))        ;  BIO_DOM     =  0.0_SPP  !DOM CONCENTRATION
      ALLOCATE (DOM_0(NNM))              ;  DOM_0       =  0.0_SPP  !DOM THRESHOLD
      ALLOCATE (FVDOM(KBV,NNM))          ;  FVDOM       =  0.0_SPP  !DOM SUBSTRACTED BY THRESHOLD

!*************  ALLOCATING  NUTRIENT        AND                    INITIALIZATION  ******************
      ALLOCATE (BIO_N(KBV,NNN))          ;  BIO_N       =  0.0_SPP  !NUTRIENT CONCENTRATION
      ALLOCATE (FVNN(KBV,NNN))           ;  FVNN        =  0.0_SPP  !NUTRIENT BUSTRACTED BY THRESHOLD
      ALLOCATE (KSN(NNN,NNP))            ;  KSN         =  0.0_SPP  !HALF-SATURAITON OF NUTRIENT
                                         ;  L_NH4N      =  0.0_SPP  !CRITICL LIGHT FOR NITRIFICATION
      ALLOCATE (NH4_NO3(KBV))            ;  NH4_NO3     =  0.0_SPP  !NITRIFICATION
      ALLOCATE (N2CB(NNN,NNB))           ;  N2CB        =  1.0_SPP  !ELEMENT RATIO IN BACTERIA
      ALLOCATE (N2CD(NNN,NND))           ;  N2CD        =  1.0_SPP  !ELEMENT RATIO IN DETRITUS
      ALLOCATE (N2CP(NNN,NNP))           ;  N2CP        =  1.0_SPP  !ELEMENT RATIO IN PHYTO
      ALLOCATE (N2CZ(NNN,NNZ))           ;  N2CZ        =  1.0_SPP  !ELEMENT RATIO IN ZOOPLANKTON
      ALLOCATE (N2CDOM(NNN,NNM))         ;  N2CDOM      =  1.0_SPP  !ELEMENT RATIO IN DOM
      ALLOCATE (N2CZB(NNN,NNB,NNZ))      ;  N2CZB       =  1.0_SPP  !ELEM. RATIO IN ACTIVE RESPIR. FROM GRAZING ON BACTERIA
      ALLOCATE (N2CZD(NNN,NND,NNZ))      ;  N2CZD       =  1.0_SPP  !ELEM. RATIO IN ACTIVE RESPIR. FROM GRAZING ON DETRITUS
      ALLOCATE (N2CZP(NNN,NNP,NNZ))      ;  N2CZP       =  1.0_SPP  !ELEM. RATIO IN ACTIVE RESPIR. FROM GRAZING ON PHYTO
      ALLOCATE (N2CZZ(NNN,NNZ,NNZ))      ;  N2CZZ       =  1.0_SPP  !ELEM. RATIO IN ACTIVE RESPIR. FROM PREDATION ON ZOO
      ALLOCATE (N_0(NNN))                ;  N_0         =  0.0_SPP  !NUTRIENT THRESHOLD
                                         ;  R_AN        =  0.0_SPP  !NITRIFICATION COEFFICIENT
      ALLOCATE (UPTAKE_PN(KBV,NNN,NNP))  ;  UPTAKE_PN   =  0.0_SPP  !NUTRIENT UPTAKE BY PHYTO

!*************************     FUNCTIONS AND PARAMETER VALUES    ***************
      DO 30 WHILE (PARMINDEX .NE. 'END')
        SELECT CASE (PARMINDEX)
          CASE ('L_FUNCTION') 
            READ (1,*) L_FUNCTION
          CASE ('G_FUNCTION') 
            READ (1,*) G_FUNCTION
          CASE ('ALPHAP');                                 !LIGHT FUNCTION COEFFICIENT
            READ(1,*) (ALPHAP(I),I=1,NNP)
          CASE ('ALPHA_U');                                !COEFFICIENT BETWEEN MULTIPLICATION & MINIMUM
            READ(1,*) (ALPHA_U(I),I=1,NNP)
          CASE ('A_TP')                                    !EXPONENTIAL FOR T FORCING
            READ(1,*) (A_TP(I),I=1,NNP)
          CASE ('ATANU_C')                                 !ATTENUATION COEFFICIENT OF CHLOROPHYLL
            READ(1,*)  ATANU_C
          CASE ('ATANU_D')                                 !ATTENUATION COEFFICIENT OF DETRITUS
            READ(1,*)  ATANU_D
          CASE ('ATANU_W')                                 !ATTENUATION COEFFICIENT OF WATER
            READ(1,*)  ATANU_W
          CASE ('BETAP')                                   !PHOTOINHIBITION COEFFICIENT
            READ(1,*) (BETAP(I),I=1,NNP)
          CASE ('CHL2C')                                   !CHLOROPHYLL TO CARBON RATIO IN PHYTO
            READ(1,*) (CHL2C(I),I=1,NNP)
          CASE ('D_DOM')                                   !ACTIVE EXUDATION COEFFICIENT
             READ(1,*) (D_DOM(I),I=1,NNP)
          CASE ('DPDOM')                                   !PASSIVE PHYTO EXUDATION COEFFICIENT
            READ(1,*) (DPDOM(I),I=1,NNP)
          CASE ('I_OPT')                                   !OPTIMAL LIGHT FOR PHYTO
            READ(1,*) (I_OPT(I),I=1,NNP)
          CASE ('K_LIGHT')                                 !HALF_SATURATION OF LIGHT FUNCTIONS
            READ(1,*) (K_LIGHT(I),I=1,NNP)
          CASE ('MPD')                                     !PHYTO MORTALITY EOFFICIENT
            READ(1,*) (MPD(I),I=1,NNP)
          CASE ('M_P')                                     !PHYTO MORTALITY ORDER
            READ(1,*) (M_P(I),I=1,NNP)
          CASE ('N_P')                                     !POWER OF LIGHT FUNCTIONS
            READ(1,*) (N_P(I),I=1,NNP)
          CASE ('PARFRAC')                                 ! PAR FRACTION
            READ(1,*) PARFRAC
          CASE ('P_0')                                     !PHYTOPANKTON THRESHOLDS
            READ(1,*) (P_0(I),I=1,NNP)
          CASE ('RP_T')                                    !EXPONENTIAL FOR T FORCING ON RESPIRATION
            READ(1,*) RP_T
          CASE ('R_P')                                      !PHYTO RESPIRATION COEFFICIENT
            READ(1,*) (R_P(I),I=1,NNP)
          CASE ('T_OPTP')                                  !OPTIMAL TEMPERATURE FOR PHYTO
           READ(1,*) (T_OPTP(I),I=1,NNP)
          CASE ('UMAX')                                    !PHYTO MAXIMUM GROWTH RATE
            READ(1,*) (UMAX(I),I=1,NNP)
          CASE ('W_P')                                     !SINKING VELOCITY OF PHYTO
            READ(1,*) (W_P(I),I=1,NNP)
!********************  ZOOPLANKTON                         PARAMETERS  **********************!
          CASE ('ACTIVE_R')                                !ACTIVE RESPIRATION COEFFICIENT
            READ(1,*) (ACTIVE_R(I),I=1,NNZ)
          CASE ('A_TZ')                                    !EXPONENTIAL FOR T FORCING
            READ(1,*) (A_TZ(I),I=1,NNZ)
          CASE ('EFFIB')                                   !GROWTH EFFICIENCY ON BACTERIA
            READ(1,*) ((EFFIB(I,J),I=1,NNB),J=1,NNZ)
          CASE ('EFFID')                                   !GROWTH EFFICIENCY ON DETRITUS
            READ(1,*) ((EFFID(I,J),I=1,NND),J=1,NNZ)
         CASE ('EFFIP')                                    !GROWTH EFFICIENCY ON PHYTO
            READ(1,*) ((EFFIP(I,J),I=1,NNP),J=1,NNZ)
          CASE ('EFFIZ')                                   !GROWTH EFFICIENCY ON ZOOPLANKTON
            READ(1,*) ((EFFIZ(I,J),I=1,NNZ),J=1,NNZ)
          CASE ('G_MAX')                                   !MAXIMUM GRAZING RATE
            READ(1,*) (G_MAX(I),I=1,NNZ)
          CASE ('K_ZG')                                    !HALF-SATURATION FOR ZOO GRAZING
            READ(1,*) (K_ZG(I),I=1,NNZ)
          CASE ('M_G')                                     !ZOOPLANKTON GRAZING POWER
            READ(1,*) (M_G(I),I=1,NNZ)
          CASE ('MZD')                                     !ZOOPLANKTON MORTALITY COEFFICIENT
            READ(1,*) (MZD(I),I=1,NNZ)
          CASE ('M_Z')                                     !ZOOPLANKTON MORTALITY POWER
            READ(1,*) (M_Z(I),I=1,NNZ)
          CASE ('P_C')                                     !GRAZING THRESHOLD
            READ(1,*) (P_C(I),I=1,NNZ)
          CASE ('R_RECRUIT')                               !RECRUIT COEFFICIENT
            READ(1,*) (R_RECRUIT(I),I=1,NNZ)
          CASE ('R_Z')                                     !ZOO PASSIVE RESPIRATION COEFFICIENT
            READ(1,*) (R_Z(I),I=1,NNZ)
          CASE ('SIGMA_B')                                 !PREFERENCE ON BACTERIA
            READ(1,*) ((SIGMA_B(I,J),I=1,NNB),J=1,NNZ)
          CASE ('SIGMA_D')                                 !PREFERENCE ON DETRITUS
            READ(1,*) ((SIGMA_D(I,J),I=1,NND),J=1,NNZ)
          CASE ('SIGMA_P')                                 !PREFERENCE ON PHYTO
            READ(1,*) ((SIGMA_P(I,J),I=1,NNP),J=1,NNZ)
          CASE ('SIGMA_Z')                                 !PREFERENCE ON ZOOPLANKTON
            READ(1,*) ((SIGMA_Z(I,J),I=1,NNZ),J=1,NNZ)
          CASE ('T_OPTZ')                                  !OPTIMAL TEMPERATURE FOR ZOOPLANKTON
            READ(1,*) (T_OPTZ(I),I=1,NNZ)
          CASE ('Z_0')                                     !ZOOPLANKTON THRESHOLLD
            READ(1,*) (Z_0(I),I=1,NNZ)

!*********************  BACTERIA                           PARAMETERS  **********************!
          CASE ('A_TB')                                    !EXPONENTAIL FOR T FORCING
            READ(1,*) (A_TB(I),I=1,NNB)
          CASE ('B_0')                                     !BACTERIA THRESHOLLD
            READ(1,*) (B_0(I),I=1,NNB)
          CASE ('DELTA_B')                                 !RATIO BETWEEN NH4 AND DON UPTAKE
            READ(1,*) (DELTA_B(I),I=1,NNB)
          CASE ('EFFIBD')                                 !GROWTH EFFICIENCY ON DETRITUS
            READ(1,*) ((EFFIBD(I,J),I=1,NND),J=1,NNB)
          CASE ('EFFIDOM')                                 !GROWTH EFFICIENCY ON DOM
            READ(1,*) ((EFFIDOM(I,J),I=1,NNM),J=1,NNB)
          CASE ('EFFIN')                                   !GROWTH EFFICIENCY ON NH4
            READ(1,*) ((EFFIN(I,J),I=1,NNN),J=1,NNB)
          CASE ('R_B')                                     !BACTERIA RESPIRATION COEFFICIENT
            READ(1,*) (R_B(I),I=1,NNB)
          CASE ('SIGMA_BD')                               !PREFERENCE ON DOM
            READ(1,*) ((SIGMA_BD(I,J),I=1,NND),J=1,NNB)
          CASE ('SIGMA_DOM')                               !PREFERENCE ON DOM
            READ(1,*) ((SIGMA_DOM(I,J),I=1,NNM),J=1,NNB)
          CASE ('SIGMA_N')                                 !PREFERENCE ON DON
            READ(1,*) ((SIGMA_N(I,J),I=1,NNN),J=1,NNB)
          CASE ('T_OPTB')                                  !OPTIMAL T FOR BACTERIA
            READ(1,*) (T_OPTB(I),I=1,NNB)
          CASE ('UBMAX')                                   !MAXIMUM GROWTH RATE OF BACTERIA
            READ(1,*) (UBMAX(I),I=1,NNB)

!***********************  DETRITUS                         PARAMETERS  **********************
          CASE ('ALPHA_BD')                                !FEEDING LOSS FROM BACTERIA TO DETRITUS
            READ(1,*) (((ALPHA_BD(I,J,K),I=1,NND),J=1,NNB),K=1,NNZ)
          CASE ('ALPHA_DAG')                               !AGGREGATION COEFFICIENT
            READ(1,*) (ALPHA_DAG(I),I=1,NND)
          CASE ('ALPHA_DD')                                !FEEDING LOSS FROM DETRITUS TO DETRITUS
            READ(1,*) (((ALPHA_DD(I,J,K),I=1,NND),J=1,NND),K=1,NNZ)
          CASE ('ALPHA_DDG')                               !DISAGGREGATION COEFFICIENT
            READ(1,*) (ALPHA_DDG(I),I=1,NND)
          CASE ('ALPHA_PD')                                !FEEDING LOSS FROM PHYTO TO DETRITUS
            READ(1,*) (((ALPHA_PD(I,J,K),I=1,NND),J=1,NNP),K=1,NNZ)
          CASE ('ALPHA_ZD')                                !FEEDING LOSS FROM ZOO TO DETRITUS
            READ(1,*) (((ALPHA_ZD(I,J,K),I=1,NND),J=1,NNZ),K=1,NNZ)
          CASE ('D_D')                                     !DISSOLUTION COEFFICIENT
            READ(1,*) (D_D(I),I=1,NND) 
          CASE ('D_0')                                     !DETRITUS THRESHOLD
            READ(1,*) (D_0(I),I=1,NND) 
          CASE ('D_RN')                                    !REMINERALIZATION COEFFICIENT
            READ(1,*) (D_RN(I),I=1,NND)
          CASE ('EPSILON_PD')                              !PHYTO MORTALITY TO DETRITUS
            READ(1,*) ((EPSILON_PD(I,J),I=1,NND),J=1,NNP)
          CASE ('EPSILON_ZD')                              !ZOO MORTALTIY TO DETRITU
            READ(1,*) ((EPSILON_ZD(I,J),I=1,NND),J=1,NNZ)
          CASE ('W_D')                                     !DETRITUS SINKING VELOCITY
            READ(1,*) (W_D(I),I=1,NND)
!**********************  DOM                               PARAMETERS  *************************!
          CASE ('ALPHA_DOM')                               !DOM AGEING COEFFICIENT
            READ(1,*) (ALPHA_DOM(I),I=1,NNM)
          CASE ('ALPHA_PDOM')                              !PHYTO EXUDATION TO DOM
            READ(1,*) ((ALPHA_PDOM(I,J),I=1,NNM),J=1,NNP)
          CASE ('ALPHA_DDOM')                              !DETRITUS DISSOLUTION TO DOM
            READ(1,*) ((ALPHA_DDOM(I,J),I=1,NNM),J=1,NND)
          CASE ('ALPHA_ZPDOM')                             !FEEDING LOSS FROM PHYTO TO DOM
            READ(1,*) (((ALPHA_ZPDOM(I,J,K),I=1,NNM),J=1,NNP),K=1,NNZ)
          CASE ('ALPHA_ZDDOM')                             !FEEDING LOSS FROM DETRITUS TO DOM
            READ(1,*) (((ALPHA_ZDDOM(I,J,K),I=1,NNM),J=1,NND),K=1,NNZ)
          CASE ('ALPHA_ZZDOM')                             !FEEDING LOSS FROM ZOO TO DOM
            READ(1,*) (((ALPHA_ZZDOM(I,J,K),I=1,NNM),J=1,NNZ),K=1,NNZ)
          CASE ('ALPHA_ZBDOM')                             !FEEDING LOSS FROM BACTERIA TO DOM
            READ(1,*) (((ALPHA_ZBDOM(I,J,K),I=1,NNM),J=1,NNB),K=1,NNZ)
          CASE ('DOM_0')                                   !DOM THRESHOLLD
            READ(1,*) (DOM_0(I),I=1,NNM)

!*****************  NUTRIENT                               PARAMETERS *******************!
          CASE ('KSN')                                     !HALF-SATURATION OF NUTRIENT
            READ(1,*) ((KSN(I,J),I=1,NNN),J=1,NNP)
          CASE ('N2CB')                                    !ELEMENT RATIO IN BACTERIA
            READ(1,*) ((N2CB(I,J),I=1,NNN),J=1,NNB)
          CASE ('N2CD')                                    !ELEMENT RATIO IN DETRITUS
            READ(1,*) ((N2CD(I,J),I=1,NNN),J=1,NND)
          CASE ('N2CP')                                    !ELEMENT RATIO IN PHYTO
            READ(1,*) ((N2CP(I,J),I=1,NNN),J=1,NNP)
          CASE ('N2CZ')                                    !ELEMENT RATIO IN ZOOPLANKTON
            READ(1,*) ((N2CZ(I,J),I=1,NNN),J=1,NNZ)
          CASE ('N2CDOM')                                  !ELEMENT RATIO IN DOM
            READ(1,*) ((N2CDOM(I,J),I=1,NNN),J=1,NNM)
          CASE ('N_0')                                     !NUTRIENT THRESHOLD
            READ(1,*) (N_0(I),I=1,NNN)
          CASE ('R_AN')                                    !NITRIFICATION COEFFICIENT
            READ(1,*) R_AN
         END SELECT
         READ (1,*)PARMINDEX
 30    CONTINUE
       CLOSE(1)
!       GOTO 99
!*******************************************************************!
! *********   CHECKING GROWTH EFFICIENCY AND FEEDLING LOSSES    ****!
!*******************************************************************!
       DO I=1,NNZ
         DO J=1,NNP                                    !FOR GRAZING ON PHYTO
           TOTALL=ACTIVE_R(I)+EFFIP(J,I)
	   print*,totall,active_r(i),effip(j,i),nnp,nnz,nnd

           IF(NND.GT.1) THEN
              IF (DETRITUS_NAME(2)=='Detritus_Si') THEN
                 DO K=1,1
                   TOTALL=TOTALL+ALPHA_PD(K,J,I)
                 END DO
              ELSE
                 DO K=1,NND
                   TOTALL=TOTALL+ALPHA_PD(K,J,I)
                 END DO
              END IF
           ELSE
               DO K=1,NND
                  TOTALL=TOTALL+ALPHA_PD(K,J,I)
               END DO
           ENDIF
	   print*,totall,alpha_pd(1,j,i),nnm


           DO K=1,NNM
              TOTALL=TOTALL+ALPHA_ZPDOM(K,J,I)
           END DO
           IF(TOTALL .NE. 1.) THEN
             PRINT*,'EFFICIENCY+LOSS OF ZOO ',I,' GRAZING ON PHYTO ',J,' UNEQUAL TO 1',TOTALL
             STOP
           END IF
         END DO

         DO J=1,NND                                    !FOR GRAZING ON DETRITUS
           TOTALL=ACTIVE_R(I)+EFFID(J,I)
           IF(NND.GT.1) THEN
             IF (DETRITUS_NAME(2)=='Detritus_Si') THEN
               DO K=1,1
                  TOTALL=TOTALL+ALPHA_DD(K,J,I)
               END DO
             ELSE
               DO K=1,NND
                  TOTALL=TOTALL+ALPHA_DD(K,J,I)
               END DO
             END IF
           ENDIF
           DO K=1,NNM
              TOTALL=TOTALL+ALPHA_ZDDOM(K,J,I)
           END DO
           IF(TOTALL > 1.) THEN
             PRINT*,'EFFICIENCY+LOSS OF ZOO ',I,' GRAZING ON DETRITUS ',J,' UNEQUAL TO 1',TOTALL
             STOP
           END IF
         END DO

         DO J=1,NNZ                                    !FOR PREDATION ON ZOO
           TOTALL=ACTIVE_R(I)+EFFIZ(J,I)
           IF(NND.GT.1) THEN
             IF (DETRITUS_NAME(2)=='Detritus_Si') THEN
               DO K=1,1
                 TOTALL=TOTALL+ALPHA_ZD(K,J,I)
               END DO
             ELSE
               DO K=1,NND
                 TOTALL=TOTALL+ALPHA_ZD(K,J,I)
               END DO
             END IF
           ENDIF
           DO K=1,NNM
              TOTALL=TOTALL+ALPHA_ZZDOM(K,J,I)
           END DO
           IF(TOTALL > 1.) THEN
             PRINT*,'EFFICIENCY+LOSS OF ZOO ',I,' PREDATION ON ZOO ',J,' UNEQUAL TO 1',TOTALL
             STOP
           END IF
         END DO

         DO J=1,NNB                                    !FOR PREDATION ON BACTERIA
           TOTALL=ACTIVE_R(I)+EFFIB(J,I)
           IF(NND.GT.1) THEN
             IF (DETRITUS_NAME(2)=='Detritus_Si') THEN
               DO K=1,1
                 TOTALL=TOTALL+ALPHA_PD(K,J,I)
               END DO
             ELSE
               DO K=1,NND
                 TOTALL=TOTALL+ALPHA_BD(K,J,I)
               END DO
             END IF
           ENDIF
           DO K=1,NNM
              TOTALL=TOTALL+ALPHA_ZBDOM(K,J,I)
           END DO
           IF(TOTALL > 1.) THEN
             PRINT*,'EFFICIENCY+LOSS OF ZOO ',I,' PREDATION ON BAC ',J,' UNEQUAL TO 1',TOTALL
             STOP
           END IF
         END DO
         TOTALL=0.
         DO J=1,NND                                   !ZOO MORTALITY TO DETRITUS
           TOTALL=TOTALL+EPSILON_ZD(J,I)
         END DO
           IF(TOTALL > 1.) THEN
             PRINT*,'MORTALITY OF ZOO ',I,' TO DETRITUS DOES NOT EQUAL TO 1',TOTALL
             STOP
           END IF
       END DO !I=1,NNZ
       DO I=1,NNP
         TOTALL=0.
         IF(NND.GT.1) THEN
           IF (DETRITUS_NAME(2)=='Detritus_Si') THEN
             DO J=1,1                                   !PHYTO MORTALITY TO DETRITUS
               TOTALL=TOTALL+EPSILON_PD(J,I)
             END DO
           ELSE
             DO J=1,NND                                   !PHYTO MORTALITY TO DETRITUS
               TOTALL=TOTALL+EPSILON_PD(J,I)
             END DO
           END IF
         ENDIF
           IF(TOTALL > 1.) THEN
             PRINT*,'MORTALITY OF PHYTO ',I,' TO DETRITUS DOES NOT EQUAL TO 1',TOTALL
             STOP
           END IF
       END DO !I=1,NNP
! 99    continue
       GOTO 88  !MOVED TO BIO_INITIAL TO PREVENT MULTIPLE PRINT
!*******************************************************************!
! *********     PRINT OUT MODEL SETUP AND PARAMETER VALUES      ****!
!*******************************************************************!
       PRINT*
       PRINT*,'*****************************************************'
       PRINT*,'**  STRUCTURE AND FUNCTION OF THE BIOLOGICAL MODEL **'
       PRINT*,'*****************************************************'

          PRINT*
          PRINT*,'MODEL STRUCTURE        : ',  MODEL_STRUCTURE
         DO I=1,NNN
          PRINT*,'                         ',  NUTRIENT_NAME(I)
         END DO
         DO I=1,NNP
          PRINT*,'                         ',  PHYTO_NAME(I)
         END DO
         DO I=1,NNZ
          PRINT*,'                         ',  ZOO_NAME(I)
         END DO
         DO I=1,NND
          PRINT*,'                         ',  DETRITUS_NAME(I)
         END DO
         DO I=1,NNM
          PRINT*,'                         ', DOM_NAME(I)
         END DO
         DO I=1,NNB
          PRINT*,'                         ', BACTERIA_NAME(I)
         END DO
          WRITE(*,'(A26,A20)')' LIGHT FUNCTION         : ', L_FUNCTION
          WRITE(*,'(A26,A20)')' GRAZING FUNCTION       : ', G_FUNCTION
          PRINT*
      PRINT*,'*********    PHYTOPLANTON PARAMETERS    **************'
        PRINT*
        IF(L_FUNCTION.NE.'EXP_LIGHT'.AND.L_FUNCTION.NE.'SL62_LIGHT') THEN
          PRINT*,'ALPHA                  : '  , (ALPHAP(I),I=1,NNP)
        END IF
          PRINT*,'L_N COMBINE            : ' , (ALPHA_U(I),I=1,NNP)
          PRINT*,'T FORCING              : ' , (A_TP(I),I=1,NNP)
          PRINT*,'CHL ATTANUATION        : ' , ATANU_C
          PRINT*,'D ATTANUATION          : ', ATANU_D
          PRINT*,'WATER ATTANUATION      : ', ATANU_W
        IF(L_FUNCTION.EQ.'PGH80_LIGHT'.OR.L_FUNCTION.EQ.'V65_LIGHT'.OR. &
           L_FUNCTION.EQ.'BWDC9_LIGHT') THEN       
          PRINT*,'BETAP                  : ', (BETAP(I),I=1,NNP)
        END IF
          PRINT*,'CHL:C                  : ',(CHL2C(I),I=1,NNP)
          PRINT*,'ACTIVE DOM EXUD.       : ',(D_DOM(I),I=1,NNP)
          PRINT*,'PASSIVE DOM EXUD.      : ', (DPDOM(I),I=1,NNP)
        IF(L_FUNCTION.EQ.'SL62_LIGHT'.OR.L_FUNCTION.EQ.'V65_LIGHT'.OR. &
           L_FUNCTION.EQ.'PE78_LIGHT') THEN       
          PRINT*,'OPTIMAL LIGHT          : ',(I_OPT(I),I=1,NNP)
        END IF
        IF(L_FUNCTION.EQ.'MM_LIGHT'.OR.L_FUNCTION.EQ.'LB_LIGHT') THEN 
          PRINT*,'LIGHT HALF SATURATTION : ', (K_LIGHT(I),I=1,NNP)
        END IF
          PRINT*,'MORTALITY              : ', (MPD(I),I=1,NNP)
          PRINT*,'MORTALITY POWER        : ', (M_P(I),I=1,NNP)
        IF(L_FUNCTION.EQ.'LB_LIGHT'.OR.L_FUNCTION.EQ.'V65_LIGHT') THEN 
          PRINT*,'POWER OF LIGHT         : ', (N_P(I),I=1,NNP)
        END IF
          PRINT*,'THRESHOLD              : ', (P_0(I),I=1,NNP)
          PRINT*,'T ON RESPIRATION       : ', RP_T
          PRINT*,'RESPIRATION            : ', (R_P(I),I=1,NNP)
          PRINT*,'OPTIMAL T              : ', (T_OPTP(I),I=1,NNP)
          PRINT*,'MAXIMUM GROWTH         : ', (UMAX(I),I=1,NNP)
          PRINT*,'SINKING VELOCITY       : ', (W_P(I),I=1,NNP)
          PRINT*

          PRINT*,'********    ZOOPLANKTON PARAMETERS   ****************'
          PRINT*
          PRINT*,'ACTIVE RESPIRATION     : ', (ACTIVE_R(I),I=1,NNZ)
          PRINT*,'T FORCING EXPONENTIAL  : ', (A_TZ(I),I=1,NNZ)
        IF (NNB.GE.1) THEN
          PRINT*,'EFFICIENCY ON BACTERIA : ',((EFFIB(I,J),I=1,NNB),J=1,NNZ)
        END IF
        IF (NND.GE.1) THEN
          PRINT*,'EFFICIENCY ON DETRITUS : ',((EFFID(I,J),I=1,NND),J=1,NNZ)
        END IF
          PRINT*,'EFFICIENCY ON PHYTO    : ',((EFFIP(I,J),I=1,NNP),J=1,NNZ)
          PRINT*,'EFFICIENCY ON ZOO      : ',((EFFIZ(I,J),I=1,NNZ),J=1,NNZ)
          PRINT*,'MAX GRAIZING RATE      : ',(G_MAX(I),I=1,NNZ)
        IF(G_FUNCTION.EQ.'RECTI_G'.OR.G_FUNCTION.EQ.'MM1_G' .AND. &
           G_FUNCTION.EQ.'MM2_G') THEN 
          PRINT*,'HALF SATURATION        : ', (K_ZG(I),I=1,NNZ)
        END IF
          PRINT*,'GRAZING POWER          : ',(M_G(I),I=1,NNZ)
          PRINT*,'MORALITY               : ',(MZD(I),I=1,NNZ)
          PRINT*,'MORTALITY POWER        : ',(M_Z(I),I=1,NNZ)
        IF(G_FUNCTION.EQ.'MM2_G') THEN 
          PRINT*,'GRAZING THRESHOLD      : ',(P_C(I),I=1,NNZ)
        END IF
        IF (NNZ.GE.1) THEN
          PRINT*,'RECRUITMENT            : ',(R_RECRUIT(I),I=1,NNZ)
        END IF
          PRINT*,'RESPIRATION            : ',(R_Z(I),I=1,NNZ)
        IF (NNB.GE.1) THEN
          PRINT*,'PREFERENCE ON BACTERIA : ',((SIGMA_B(I,J),I=1,NNB),J=1,NNZ)
        END IF
        IF (NND.GE.1) THEN
          PRINT*,'PREFERENCE ON DETRITUS : ',((SIGMA_D(I,J),I=1,NND),J=1,NNZ)
        END IF
          PRINT*,'PREFERENCE ON PHYTO    : ',((SIGMA_P(I,J),I=1,NNP),J=1,NNZ)
          PRINT*,'PREFERENCE ON ZOO      : ',((SIGMA_Z(I,J),I=1,NNZ),J=1,NNZ)
          PRINT*,'OPTIMAL T              : ',(T_OPTZ(I),I=1,NNZ)
          PRINT*,'ZOO THRESHOLD          : ',(Z_0(I),I=1,NNZ)
          PRINT*
          PRINT*,'*********    NUTRIENT PARAMETERS    ***********'
          PRINT*
          PRINT*,'HALF-SATURATION',((KSN(I,J),I=1,NNN),J=1,NNP)
        IF (NNB.GE.1) THEN
          PRINT*,'ELEMENT RATIO IN BAC.  : ',((N2CB(I,J),I=1,NNN),J=1,NNB)
        END IF
        IF (NND.GE.1) THEN
          PRINT*,'ELEMENT RATIO IN D     : ',((N2CD(I,J),I=1,NNN),J=1,NND)
        END IF
          PRINT*,'ELEMENT RATIO IN PHYTO : ',((N2CP(I,J),I=1,NNN),J=1,NNP)
          PRINT*,'ELEMENT RATIO IN ZOO   : ',((N2CZ(I,J),I=1,NNN),J=1,NNZ)
        IF (NNM.GE.1)  THEN
          PRINT*,'ELEMENT RATIO IN DOM   : ', ((N2CDOM(I,J),I=1,NNN),J=1,NNM)
        END IF
          PRINT*,'THRESHOLD              : ', (N_0(I),I=1,NNN)
        IF (NO3_ON)    THEN
          PRINT*,'NITRIFICATION RATE     : ', R_AN
        END IF
          PRINT*
        IF (NNB.GE.1) THEN
          PRINT*,'*********  BACTERIA PARAMETERS  ************'
          PRINT*
          PRINT*,'T_FORCING              : ',(A_TB(I),I=1,NNB)
          PRINT*,'THRESHOLD              : ', (B_0(I),I=1,NNB)
          PRINT*,'RATIO OF NH4 VS DON    : ',(DELTA_B(I),I=1,NNB)
          PRINT*,'EFFICIENCY OF DETRITUS : ',((EFFIBD(I,J),I=1,NND),J=1,NNB)
          PRINT*,'EFFICIENCY OF DOM      : ',((EFFIDOM(I,J),I=1,NNM),J=1,NNB)
          PRINT*,'EFFICIENCY OF NUTRIENT : ', ((EFFIN(I,J),I=1,NNN),J=1,NNB)
          PRINT*,'RESPIRATION            : ',(R_B(I),I=1,NNB)
          PRINT*,'PREFERENCE ON DETRITUS : ',((SIGMA_BD(I,J),I=1,NND),J=1,NNB)
          PRINT*,'PREFERENCE ON DOM      : ',((SIGMA_DOM(I,J),I=1,NNM),J=1,NNB)
          PRINT*,'PREFERENCE ON NUTRIENT : ',((SIGMA_N(I,J),I=1,NNN),J=1,NNB)
          PRINT*,'OPTIMAL TEMPERATURE    : ',(T_OPTB(I),I=1,NNB)
          PRINT*,'MAXIMUM GROWTH RATE    : ',(UBMAX(I),I=1,NNB)
          PRINT*
        END IF
        IF (NND.GE.1) THEN
          PRINT*, '********  DETRITUS PARAMETERS  *************'
          PRINT*
        IF (NNB.GE.1) THEN
          PRINT*,'GRAZING ON B TO D      : ', (((ALPHA_BD(I,J,K),I=1,NND),J=1,NNB),K=1,NNZ)
        END IF
          PRINT*,'AGGREGATION            : ', (ALPHA_DAG(I),I=1,NND)
          PRINT*,'GRAZING LOSS ON D TO D : ', (((ALPHA_DD(I,J,K),I=1,NND),J=1,NND),K=1,NNZ)
          PRINT*,'DISAGGREGATION         : ', (ALPHA_DAG(I),I=1,NND)
          PRINT*,'GRAZING LOSS ON P TO D : ', (((ALPHA_PD(I,J,K),I=1,NND),J=1,NNP),K=1,NNZ)
          PRINT*,'GRAZING LOSS ON Z TO D : ', (((ALPHA_ZD(I,J,K),I=1,NND),J=1,NNZ),K=1,NNZ)
        IF (NNM.GE.1) THEN
          PRINT*,'DISSOLUTION            : ', (D_DOM(I),I=1,NND)
        END IF
          PRINT*,'REMINERALIZATION       : ', (D_RN(I),I=1,NND)
          PRINT*,'THRESHOLD              : ', (D_0(I),I=1,NND)
          PRINT*,'P MORTALITY TO DETRITUS: ', ((EPSILON_PD(I,J),I=1,NND),J=1,NNP)
          PRINT*,'Z MORTALITY TO DETRITUS: ', ((EPSILON_ZD(I,J),I=1,NND),J=1,NNZ)
          PRINT*,'SINING VELOCITY        : ', (W_D(I),I=1,NND)
          PRINT*
        END IF
        IF (NNM.GE.1) THEN
          PRINT*,'*********     DOM PARAMETERS     *************'
          PRINT*
          PRINT*,'DOM AGEING COEFFICIENT : ', (ALPHA_DOM(I),I=1,NNM)
          PRINT*,'PHYTO EXUDATION        : ', ((ALPHA_PDOM(I,J),I=1,NNM),J=1,NNP)
          PRINT*,'DETRITUS DISSOLUTION   : ', ((ALPHA_DDOM(I,J),I=1,NNM),J=1,NND)
        IF (NNB.GE.1) THEN
          PRINT*,'GRAZING LOSS ON B > DOM: ', (((ALPHA_ZBDOM(I,J,K),I=1,NNM),J=1,NNB),K=1,NNZ)
        END IF
        IF (NND.GE.1) THEN
          PRINT*,'GRAZING LOSS ON D > DOM: ', (((ALPHA_ZDDOM(I,J,K),I=1,NNM),J=1,NND),K=1,NNZ)
        END IF
          PRINT*,'GRAZING LOSS ON P > DOM: ', (((ALPHA_ZPDOM(I,J,K),I=1,NNM),J=1,NNP),K=1,NNZ)
          PRINT*,'GRAZING LOSS ON Z > DOM: ', (((ALPHA_ZZDOM(I,J,K),I=1,NNM),J=1,NNZ),K=1,NNZ)
          PRINT*,'THRESHOLD              : ', (DOM_0(I),I=1,NNM)
          PRINT*
        END IF
       PRINT*
       PRINT*,'*****************************************************'
       PRINT*,'*********    END OF BIOLOGICAL MODEL     ************'
       PRINT*,'*****************************************************'
       PRINT*
 88    CONTINUE
       RETURN
END SUBROUTINE GET_PARAMETER





      SUBROUTINE GET_PARAMETER_NEW(fname)
      character(len=*)    :: fname
      logical             :: fexist
      INTEGER  ::  I,J,K
      integer  :: linenum,k1,iscan
      real(spp),allocatable :: tmp(:)
!************************    BEGINNING EXECUTABLE    **********************
      KBVM1=KBV-1

      !ensure sediment parameter file exists    
      inquire(file=trim(fname),exist=fexist)
      if(.not.fexist)then
        write(*,*)'biological parameter file: ',trim(fname),' does not exist'
        write(*,*)'stopping'
        stop
      end if

!      OPEN (1, FILE='bio_model_new.in')
!      biofile='bio_model_new.in'

      OPEN (1, FILE=trim(fname),status='old')
      biofile = trim(fname)

     !read in model struture
     Call Get_Val(MODEL_STRUCTURE,biofile,'STRUCTURE',line=linenum,echo=.true.) 

      SELECT CASE (MODEL_STRUCTURE)
         CASE ('NPZ')
            NNN=1
            NNP=1
            NNZ=1
            NNM=0
            NNB=0
            NND=0
            NNA=0
            NTT=NNN+NNP+NNZ+NNM+NNB+NND+NNA   !TOTAL VARIABLE NUMBER FOR VERTICAL MIXING
            INN=1
            INP=INN+NNN
            INZ=INP+NNP
            INM=INZ+NNZ
            INB=INM+NNM
            IND=INB+NNB
            INA=IND+NND
            ALLOCATE (NUTRIENT_NAME(NNN))
            ALLOCATE (NUTRIENT_UNIT(NNN))
            ALLOCATE (PHYTO_NAME(NNP))
            ALLOCATE (PHYTO_UNIT(NNP))
            ALLOCATE (ZOO_NAME(NNZ))
            ALLOCATE (ZOO_UNIT(NNZ))
            NUTRIENT_NAME(1)   = 'Nitrogen'
	    NUTRIENT_UNIT(1)   = 'mmole N m-3'
            PHYTO_NAME(1)      = 'Phytoplankton'
	    PHYTO_UNIT(1)      = 'mmole C m-3'
            ZOO_NAME(1)        = 'Zooplankton'
	    ZOO_UNIT(1)        = 'mmole C m-3'

	  CASE ('NPZD')
            NNN=1
            NNP=1
            NNZ=1
            NNM=0
            NNB=0
            NND=1
            NNA=0
            NTT=NNN+NNP+NNZ+NNM+NNB+NND+NNA   !TOTAL VARIABLE NUMBER FOR VERTICAL MIXING
            INN=1
            INP=INN+NNN
            INZ=INP+NNP
            INM=INZ+NNZ
            INB=INM+NNM
            IND=INB+NNB
            INA=IND+NND
            ALLOCATE (NUTRIENT_NAME(NNN))
            ALLOCATE (NUTRIENT_UNIT(NNN))
            ALLOCATE (PHYTO_NAME(NNP))
            ALLOCATE (PHYTO_UNIT(NNP))
            ALLOCATE (ZOO_NAME(NNZ))
            ALLOCATE (ZOO_UNIT(NNZ))
            ALLOCATE (DETRITUS_NAME(NND))
            ALLOCATE (DETRITUS_UNIT(NND))
            NUTRIENT_NAME(1)   = 'Nitrogen'
	    NUTRIENT_UNIT(1)   = 'mmole N m-3'
            PHYTO_NAME(1)      = 'Phytoplankton'
	    PHYTO_UNIT(1)      = 'mmole C m-3'
            ZOO_NAME(1)        = 'Zooplankton'
	    ZOO_UNIT(1)        = 'mmole C m-3'
            DETRITUS_NAME(1)   = 'Detritus'
	    DETRITUS_UNIT(1)   = 'mmole C m-3'

	  CASE ('FASHAM')
            NNN=2
            NNP=1
            NNZ=1
            NNM=1
            NNB=1
            NND=1
            NNA=0
            NTT=NNN+NNP+NNZ+NNM+NNB+NND+NNA   !TOTAL VARIABLE NUMBER FOR VERTICAL MIXING
            INN=1
            INP=INN+NNN
            INZ=INP+NNP
            INM=INZ+NNZ
            INB=INM+NNM
            IND=INB+NNB
            INA=IND+NND
            ALLOCATE (NUTRIENT_NAME(NNN))
            ALLOCATE (NUTRIENT_UNIT(NNN))
            ALLOCATE (PHYTO_NAME(NNP))
            ALLOCATE (PHYTO_UNIT(NNP))
            ALLOCATE (ZOO_NAME(NNZ))
            ALLOCATE (ZOO_UNIT(NNZ))
            ALLOCATE (DETRITUS_NAME(NND))
            ALLOCATE (DETRITUS_UNIT(NND))
            ALLOCATE (DOM_NAME(NNM))
            ALLOCATE (DOM_UNIT(NNM))
            ALLOCATE (BACTERIA_NAME(NNB))
            ALLOCATE (BACTERIA_UNIT(NNB))
            NUTRIENT_NAME(1)   = 'NH4'
	    NUTRIENT_UNIT(1)   = 'mmole N m-3'
            NUTRIENT_NAME(2)   = 'NO3'
	    NUTRIENT_UNIT(2)   = 'mmole N m-3'
            PHYTO_NAME(1)      = 'Phytoplankton'
	    PHYTO_UNIT(1)      = 'mmole C m-3'
            ZOO_NAME(1)        = 'Zooplankton'
	    ZOO_UNIT(1)        = 'mmole C m-3'
            DETRITUS_NAME(1)   = 'Detritus'
	    DETRITUS_UNIT(1)   = 'mmole C m-3'
            DOM_NAME(1)        = 'DOM'
	    DOM_UNIT(1)        = 'mmole C m-3'
            BACTERIA_NAME(1)   = 'Bacteria'
	    BACTERIA_UNIT(1)   = 'mmole C m-3'

          CASE ('N2P2Z2D2')
            NNN=3
            NNP=2
            NNZ=2
            NNM=0
            NNB=0
            NND=2
            NNA=0
            NTT=NNN+NNP+NNZ+NNM+NNB+NND+NNA   !TOTAL VARIABLE NUMBER FOR VERTICAL MIXING
            INN=1
            INP=INN+NNN
            INZ=INP+NNP
            INM=INZ+NNZ
            INB=INM+NNM
            IND=INB+NNB
            INA=IND+NND
            ALLOCATE (NUTRIENT_NAME(NNN))
            ALLOCATE (NUTRIENT_UNIT(NNN))
            ALLOCATE (PHYTO_NAME(NNP))
            ALLOCATE (PHYTO_UNIT(NNP))
            ALLOCATE (ZOO_NAME(NNZ))
            ALLOCATE (ZOO_UNIT(NNZ))
            ALLOCATE (DETRITUS_NAME(NND))
            ALLOCATE (DETRITUS_UNIT(NND))
            NUTRIENT_NAME(1)   = 'NH4'
	    NUTRIENT_UNIT(1)   = 'mmole N m-3'
            NUTRIENT_NAME(2)   = 'NO3'
	    NUTRIENT_UNIT(2)   = 'mmole N m-3'
            NUTRIENT_NAME(3)   = 'SiO3'
	    NUTRIENT_UNIT(3)   = 'mmole Si m-3'
            PHYTO_NAME(1)      = 'Small_phyto'
	    PHYTO_UNIT(1)      = 'mmole C m-3'
            PHYTO_NAME(2)      = 'Large_phyto'
	    PHYTO_UNIT(2)      = 'mmole C m-3'
            ZOO_NAME(1)        = 'Microzooplankton'
	    ZOO_UNIT(1)        = 'mmole C m-3'
            ZOO_NAME(2)        = 'Mesozooplankton'
	    ZOO_UNIT(2)        = 'mmole C m-3'
            DETRITUS_NAME(1)   = 'Detritus_N'
	    DETRITUS_UNIT(1)   = 'mmole N m-3'
            DETRITUS_NAME(2)   = 'Detritus_Si'
	    DETRITUS_UNIT(2)   = 'mmole Si m-3'

          CASE ('N2P2Z2D2B')
            NNN=2
            NNP=2
            NNZ=2
            NNM=0
            NNB=1
            NND=2
            NNA=0
            NTT=NNN+NNP+NNZ+NNM+NNB+NND+NNA   !TOTAL VARIABLE NUMBER FOR VERTICAL MIXING
            INN=1
            INP=INN+NNN
            INZ=INP+NNP
            INM=INZ+NNZ
            INB=INM+NNM
            IND=INB+NNB
            INA=IND+NND
            ALLOCATE (NUTRIENT_NAME(NNN))
            ALLOCATE (NUTRIENT_UNIT(NNN))
            ALLOCATE (PHYTO_NAME(NNP))
            ALLOCATE (PHYTO_UNIT(NNP))
            ALLOCATE (ZOO_NAME(NNZ))
            ALLOCATE (ZOO_UNIT(NNZ))
            ALLOCATE (DETRITUS_NAME(NND))
            ALLOCATE (DETRITUS_UNIT(NND))
            ALLOCATE (BACTERIA_NAME(NNB))
            ALLOCATE (BACTERIA_UNIT(NNB))
            NUTRIENT_NAME(1)   = 'PO4'
	    NUTRIENT_UNIT(1)   = 'mmole P m-3'
            NUTRIENT_NAME(2)   = 'SiO3'
	    NUTRIENT_UNIT(2)   = 'mmole Si m-3'
            PHYTO_NAME(1)      = 'Small_phyto'
	    PHYTO_UNIT(1)      = 'mmole C m-3'
            PHYTO_NAME(2)      = 'Large_phyto'
	    PHYTO_UNIT(2)      = 'mmole C m-3'
            ZOO_NAME(1)        = 'Microzooplankton'
	    ZOO_UNIT(1)        = 'mmole C m-3'
            ZOO_NAME(2)        = 'Mesozooplankton'
	    ZOO_UNIT(2)        = 'mmole C m-3'
            DETRITUS_NAME(1)   = 'Detritus_P'
	    DETRITUS_UNIT(1)   = 'mmole P m-3'
            DETRITUS_NAME(2)   = 'Detritus_Si'
	    DETRITUS_UNIT(2)   = 'mmole Si m-3'
            BACTERIA_NAME(1)   = 'Bacteria'
	    BACTERIA_UNIT(1)   = 'mmole C m-3'
	            
          CASE ('N2P2Z2D2SB')
            NNN=2
            NNP=2
            NNZ=2
            NNM=1
            NNB=2
            NND=2
            NNA=0
            NTT=NNN+NNP+NNZ+NNM+NNB+NND+NNA   !TOTAL VARIABLE NUMBER FOR VERTICAL MIXING
            INN=1
            INP=INN+NNN
            INZ=INP+NNP
            INM=INZ+NNZ
            INB=INM+NNM
            IND=INB+NNB
            INA=IND+NND
            ALLOCATE (NUTRIENT_NAME(NNN))
            ALLOCATE (NUTRIENT_UNIT(NNN))
            ALLOCATE (PHYTO_NAME(NNP))
            ALLOCATE (PHYTO_UNIT(NNP))
            ALLOCATE (ZOO_NAME(NNZ))
            ALLOCATE (ZOO_UNIT(NNZ))
            ALLOCATE (DETRITUS_NAME(NND))
            ALLOCATE (DETRITUS_UNIT(NND))
            ALLOCATE (DOM_NAME(NNM))
            ALLOCATE (DOM_UNIT(NNM))
            ALLOCATE (BACTERIA_NAME(NNB))
            ALLOCATE (BACTERIA_UNIT(NNB))
            NUTRIENT_NAME(1)   = 'NH4'
	    NUTRIENT_UNIT(1)   = 'mmole N m-3'
            NUTRIENT_NAME(2)   = 'NO3'
	    NUTRIENT_UNIT(2)   = 'mmole N m-3'
            PHYTO_NAME(1)      = 'Small_phyto'
	    PHYTO_UNIT(1)      = 'mmole C m-3'
            PHYTO_NAME(2)      = 'Large_phyto'
	    PHYTO_UNIT(2)      = 'mmole C m-3'
            ZOO_NAME(1)        = 'Microzooplankton'
	    ZOO_UNIT(1)        = 'mmole C m-3'
            ZOO_NAME(2)        = 'Mesozooplankton'
	    ZOO_UNIT(2)        = 'mmole C m-3'
            DETRITUS_NAME(1)   = 'Small_detritus'
	    DETRITUS_UNIT(1)   = 'mmole C m-3'
            DETRITUS_NAME(2)   = 'Large_detritus'
	    DETRITUS_UNIT(2)   = 'mmole C m-3'
            DOM_NAME(1)        = 'DOM'
	    DOM_UNIT(1)        = 'mmole C m-3'
            BACTERIA_NAME(1)   = 'Bacteria'
	    BACTERIA_UNIT(1)   = 'mmole C m-3'

         CASE('NEW')
            !ASSIGNMENT FOR NEW MODEL

            !read in number of PHYTOPLANKTON
            Call Get_Val(NNP,biofile,'PHYTOPLANKTON_NUM',line=linenum,echo=.true.)
            ALLOCATE (PHYTO_NAME(NNP), PHYTO_UNIT(NNP))
            !read in PHYTOPLANKTON names and units
            Call Get_Val_Array(PHYTO_NAME,biofile,'PHYTOPLANKTON_NAME',NNP,echo=.true.)
            Call Get_Val_Array(PHYTO_UNIT,biofile,'PHYTOPLANKTON_UNIT',NNP,echo=.true.)


            !read in number of ZOOPLANKTON
            Call Get_Val(NNZ,biofile,'ZOOPLANKTON_NUM',line=linenum,echo=.true.)
            ALLOCATE (ZOO_NAME(NNZ), ZOO_UNIT(NNZ))
            !read in ZOOPLANKTON names and units
            Call Get_Val_Array(ZOO_NAME,biofile,'ZOOPLANKTON_NAME',NNZ,echo=.true.)
            Call Get_Val_Array(ZOO_UNIT,biofile,'ZOOPLANKTON_UNIT',NNZ,echo=.true.)


            !read in number of NUTRIENT
            Call Get_Val(NNN,biofile,'NUTRIENT_NUM',line=linenum,echo=.true.)
            ALLOCATE (NUTRIENT_NAME(NNN), NUTRIENT_UNIT(NNN))
            !read in NUTRIENT names and units
            Call Get_Val_Array(NUTRIENT_NAME,biofile,'NUTRIENT_NAME',NNN,echo=.true.)
            Call Get_Val_Array(NUTRIENT_UNIT,biofile,'NUTRIENT_UNIT',NNN,echo=.true.)


            !read in number of BACTERIA
            Call Get_Val(NNB,biofile,'BACTERIA_NUM',line=linenum,echo=.true.)
            ALLOCATE (BACTERIA_NAME(NNB), BACTERIA_UNIT(NNB))
            !read in BACTERIA names and units
            Call Get_Val_Array(BACTERIA_NAME,biofile,'BACTERIA_NAME',NNB,echo=.true.)
            Call Get_Val_Array(BACTERIA_UNIT,biofile,'BACTERIA_UNIT',NNB,echo=.true.)

            !read in number of DETRITUS
            Call Get_Val(NND,biofile,'DETRITUS_NUM',line=linenum,echo=.true.)
            ALLOCATE (DETRITUS_NAME(NND), DETRITUS_UNIT(NND))
            !read in DETRITUS names and units
            Call Get_Val_Array(DETRITUS_NAME,biofile,'DETRITUS_NAME',NND,echo=.true.)
            Call Get_Val_Array(DETRITUS_UNIT,biofile,'DETRITUS_UNIT',NND,echo=.true.)

            !read in number of DOM
            Call Get_Val(NNM,biofile,'DOM_NUM',line=linenum,echo=.true.)
            ALLOCATE (DOM_NAME(NNM),DOM_UNIT(NNM))
            !read in DOM names and units
            Call Get_Val_Array(DOM_NAME,biofile,'DOM_NAME',NNM,echo=.true.)
            Call Get_Val_Array(DOM_UNIT,biofile,'DOM_UNIT',NNM,echo=.true.)

            !read in number of AUXILIARY
            Call Get_Val(NNA,biofile,'AUXILIARY_NUM',line=linenum,echo=.true.)
            ALLOCATE (AUX_NAME(NNA),AUX_UNIT(NNA))
            !read in AUXILIARY names and units
            Call Get_Val_Array(AUX_NAME,biofile,'AUXILIARY_NAME',NNA,echo=.true.)
            Call Get_Val_Array(AUX_UNIT,biofile,'AUXILIARY_UNIT',NNA,echo=.true.)

            NTT=NNN+NNP+NNZ+NNM+NNB+NND+NNA   !TOTAL VARIABLE NUMBER FOR VERTICAL MIXING
            INN=1
            INP=INN+NNN
            INZ=INP+NNP
            INM=INZ+NNZ
            INB=INM+NNM
            IND=INB+NNB
            INA=IND+NND
      END SELECT

!*******************    BIO_NAME ASSIGNMENT FOR NETCDF OUT PUT   ***************
            ALLOCATE (BIO_NAME(NTT,4))
            BIO_NAME(:,3)  = 'Concentration'
            BIO_NAME(:,4)  = 'Array'
            DO J=INN,NNN
              BIO_NAME(J,1)=NUTRIENT_NAME(J)
              BIO_NAME(J,2)=NUTRIENT_UNIT(J)
            END DO
            DO J=INP,INP+NNP-1
              BIO_NAME(J,1)=PHYTO_NAME(J-INP+1)
              BIO_NAME(J,2)=PHYTO_UNIT(J-INP+1)
            END DO
            DO J=INZ,INZ+NNZ-1
              BIO_NAME(J,1)=ZOO_NAME(J-INZ+1)
              BIO_NAME(J,2)=ZOO_UNIT(J-INZ+1)
            END DO
            DO J=INM,INM+NNM-1
              BIO_NAME(J,1)=DOM_NAME(J-INM+1)
              BIO_NAME(J,2)=DOM_UNIT(J-INM+1)
            END DO
            DO J=INB,INB+NNB-1
              BIO_NAME(J,1)=BACTERIA_NAME(J-INB+1)
              BIO_NAME(J,2)=BACTERIA_UNIT(J-INB+1)
            END DO
            DO J=IND,IND+NND-1
              BIO_NAME(J,1)=DETRITUS_NAME(J-IND+1)
              BIO_NAME(J,2)=DETRITUS_UNIT(J-IND+1)
            END DO
            DO J=INA,INA+NNA-1
              BIO_NAME(J,1)=AUX_NAME(J-INA+1)
              BIO_NAME(J,2)=AUX_UNIT(J-INA+1)
            END DO
            NO3_ON = .FALSE.                                       !NH4 INHIBITION
            PO4_ON = .FALSE.                                       !BAC UPTAKE ON PO4
            DOCDON = .FALSE.                                       !BOTH CONSIDERED
            IF(NNN.GT.1) THEN
              IF(NUTRIENT_NAME(2).EQ.'NO3') NO3_ON = .TRUE.
            ENDIF
            IF(NNM.GE.2) THEN
              IF(DOM_NAME(1).EQ.'DOC'.AND.DOM_NAME(2).EQ.'DON') DOCDON = .TRUE.
            END IF
	    IF(NUTRIENT_NAME(1).EQ.'PO4') PO4_ON = .TRUE.

!************* ALLOCATING                   GLOBAL PARAMETERS      AND INITIALIZAITON
      ALLOCATE (BIO_VAR(KBV,NTT))        ;  BIO_VAR     =  0.0_SPP  !ENSEMBLE VARIABLES FOR MIXING
      ALLOCATE (DELTA_D(KBV))            ;  DELTA_D     =  1.0_SPP  !LAYER THICKNESS
      ALLOCATE (DELTA_Z(KBV))            ;  DELTA_Z     =  1.0_SPP  !THICKNESS BETWEEN LAYER CENTRE
      ALLOCATE (DEPTH_Z(KBV))            ;  DEPTH_Z     =  1.0_SPP  !DEPTH AT EACH LAYER CENTER
      ALLOCATE (KM_BIO(KBV))             ;  KM_BIO      =  0.0_SPP  !DIFFUSIVITY 
      ALLOCATE (T_BIO(KBV))              ;  T_BIO       =  0.0_SPP  !TEMPERATURE

!************* ALLOCATING                   PHYTOPLANKTON  AND     INITIALIZAITON *************
      ALLOCATE (ALPHAP(NNP))             ;  ALPHAP      =  0.0_SPP  !LIGHT FUNCTION COEFFICIENT
      ALLOCATE (ALPHA_U(NNP))            ;  ALPHA_U     =  0.5_SPP  !COEFFICIENT BETWEEN MULTIPLICAITON & MINIMUM
      ALLOCATE (A_TP(NNP))               ;  A_TP        =  0.0_SPP  !EXPONENTIAL FOR T FORCING
                                         ;  ATANU_C     =  0.03_SPP !ATTENUATION COEFFICIENT OF CHLOROPHYLL
                                         ;  ATANU_D     =  0.01_SPP !ATTENUATION COEFFICIENT OF DETRITUS
                                         ;  ATANU_W     =  0.04_SPP !ATTENUATION COEFFICIENT OF WATER
      ALLOCATE (BETAP(NNP))              ;  BETAP       =  0.0_SPP  !PHOTOINHIBITION COEFFICIENT
      ALLOCATE (BIO_P(KBV,NNP))          ;  BIO_P       =  0.0_SPP  !PHYTO BIOMASS
      ALLOCATE (CHL2C(NNP))              ;  CHL2C       =  0.24_SPP !CHLOROPHYLL TO CARBON RATIO IN PHYTO
      ALLOCATE (D_DOM(NNP))              ;  D_DOM       =  0.0_SPP  !ACTIVE EXUDATION COEFFICIENT
      ALLOCATE (DPDOM(NNP))              ;  DPDOM       =  0.0_SPP  !COEFFICIENT OF PHYTO EXUDATION
      ALLOCATE (FVP(KBV,NNP))            ;  FVP         =  0.0_DPP  !PHYTO AFTER THRESHOLD SUBSTRACTION
      ALLOCATE (I_OPT(NNP))              ;  I_OPT       =  100._SPP !OPTIMAL LIGHT FOR PHYTO
                                         ;  IRRAD0      =  200._SPP !IRRADIANCE AT THE SEA SURFACE
      ALLOCATE (K_LIGHT(NNP))            ;  K_LIGHT     =  50._SPP  !HALF_SATURATION OF LIGHT FUNCTIONS
      ALLOCATE (MPD(NNP))                ;  MPD         =  0.0_SPP  !PHYTO MORTALITY EOFFICIENT
      ALLOCATE (M_P(NNP))                ;  M_P         =  1.0_SPP  !PHYTO MORTALITY ORDER
      ALLOCATE (N_P(NNP))                ;  N_P         =  1.0_SPP  !POWER OF LIGHT FUNCTIONS
                                         ;  PARFRAC     =  0.43_SPP
      ALLOCATE (P_D(KBV,NNP))            ;  P_D         =  0.0_SPP  !PHYTO MORTALITY
      ALLOCATE (P_DOM(KBV,NNP))          ;  P_DOM       =  0.0_SPP  !PHYTO EXUDATION OF DOM
      ALLOCATE (P_N(KBV,NNP))            ;  P_N         =  0.0_SPP  !PHYTO RESPIRATION
      ALLOCATE (P_0(NNP))                ;  P_0         =  0.0_SPP  !PHYTO THRESHOLD
                                         ;  RP_T        =  0.0_SPP  !EXPONENTIAL FOR T FORCING ON RESPIRATION
      ALLOCATE (R_P(NNP))                ;  R_P         =  1.0E-7_SPP  !PHYTO RESPIRATION COEFFICIENT
      ALLOCATE (T_OPTP(NNP))             ;  T_OPTP      =  20._SPP  !OPTIMAL TEMPERATURE FOR PHYTO
      ALLOCATE (ULR(NNP))                ;  ULR         =  0.0_SPP  !LIGHT GROWHT FACTOR FOR PHYTO
      ALLOCATE (UMAX(NNP))               ;  UMAX        =  1.16E-5_SPP !PHYTO MAXIMUM GROWTH RATE
      ALLOCATE (UNR(NNP,NNN+1))          ;  UNR         =  1.0_SPP  !NUTRIENT GROWTH FACTOR FOR PHYTO
      ALLOCATE (UNRMIN(NNP))             ;  UNRMIN      =  1.0_SPP  !MINIMUM NUTRIENT FACTOR FOR PHYTO
      ALLOCATE (U_P(KBV,NNP))            ;  U_P         =  0.0_SPP  !PHYTO GROWTH
      ALLOCATE (UTR(NNP))                ;  UTR         =  1.0_SPP  !T GORWTH FACTOR FOR PHYTO
      ALLOCATE (W_P(NNP))                ;  W_P         =  1.0_DPP  !SINKING VELOCITY OF PHYTO
      ALLOCATE (WSNK_P(KBV,NNP))         ;  WSNK_P      =  0.0_DPP  !SINKING FLUX OF PHYTO

!**************ALLOCATING                   ZOOPLANKTON    AND     INITIALIZATION *****************
      ALLOCATE (ACTIVE_R(NNZ))           ;  ACTIVE_R    =  0.0_SPP  !ZOO ACTIVE RESPIRATION COEFFICIENT
      ALLOCATE (A_TZ(NNZ))               ;  A_TZ        =  0.0_SPP  !T EXPONENTAIL COEFFICIENT
      ALLOCATE (BIO_Z(KBV,NNZ))          ;  BIO_Z       =  0.0_SPP  !ZOOPLANKTON BIOMASS
      ALLOCATE (EFFIB(NNB,NNZ))          ;  EFFIB       =  0.3_SPP  !GROWTH EFFICIENCY ON BACTERIA
      ALLOCATE (EFFID(NND,NNZ))          ;  EFFID       =  0.3_SPP  !GROWTH EFFICIENCY ON DETRITUS
      ALLOCATE (EFFIP(NNP,NNZ))          ;  EFFIP       =  0.3_SPP  !GROWTH EFFICIENCY ON PHYTO
      ALLOCATE (EFFIZ(NNZ,NNZ))          ;  EFFIZ       =  0.3_SPP  !GROWTH EFFICIENCY ON ZOOPLANKTON
      ALLOCATE (FVZ(KBV,NNP))            ;  FVZ         =  0.0_SPP  !ZOOPLANKTON SUBSTRACTED OF THRESHOLD
      ALLOCATE (G_B(KBV,NNB,NNZ))        ;  G_B         =  0.0_SPP  !PREDATION ON BACTERIA
      ALLOCATE (G_D(KBV,NND,NNZ))        ;  G_D         =  0.0_SPP  !GRAZING ON DETRITUS
      ALLOCATE (G_MAX(NNZ))              ;  G_MAX       =  1.16E-5_SPP !MAXIMUM GRAZING RATE
      ALLOCATE (G_P(KBV,NNP,NNZ))        ;  G_P         =  0.0_SPP  !GRAZING ON PHYTOPLANKTON
      ALLOCATE (G_RATE(NNZ))             ;  G_RATE      =  0.0_SPP  !GRAZING RATE
      ALLOCATE (G_Z(KBV,NNZ,NNZ))        ;  G_Z         =  0.0_SPP  !PREDATION OF ZOOPLANKTON
      ALLOCATE (K_ZG(NNZ))               ;  K_ZG        =  1.0_SPP  !HALF-SATURATION FOR ZOO GRAZING
      ALLOCATE (M_G(NNZ))                ;  M_G         =  1.0_SPP  !ZOOPLANKTON GRAZING POWER
      ALLOCATE (MZD(NNZ))                ;  MZD         =  0.0_SPP  !ZOOPLANKTON MORTALITY COEFFICIENT
      ALLOCATE (M_Z(NNZ))                ;  M_Z         =  1.0_SPP  !ZOOPLANKTON MORTALITY POWER
      ALLOCATE (P_C(NNZ))                ;  P_C         =  0.0_SPP  !GRAZING THRESHOLD
      ALLOCATE (RECRUIT_LOSS(0:NNZ))     ;  RECRUIT_LOSS=  0.0_SPP  !RECRUITMENT LOSS
      ALLOCATE (R_RECRUIT(NNZ))          ;  R_RECRUIT   =  0.0_SPP  !RECRUIT COEFFICIENT
      ALLOCATE (R_Z(NNZ))                ;  R_Z         =  1.E-7_SPP   !ZOOPLANKTON RESPIRATION COEFFICIENT
      ALLOCATE (SIGMA_B(NNB,NNZ))        ;  SIGMA_B     =  1.0_SPP  !PREFERENCE ON BACTERIA
      ALLOCATE (SIGMA_D(NND,NNZ))        ;  SIGMA_D     =  1.0_SPP  !PREFERENCE ON DETRITUS
      ALLOCATE (SIGMA_P(NNP,NNZ))        ;  SIGMA_P     =  1.0_SPP  !PREFERENCE ON PHYTO
      ALLOCATE (SIGMA_Z(NNZ,NNZ))        ;  SIGMA_Z     =  1.0_SPP  !PREFERENCE ON ZOOPLANKTON
      ALLOCATE (T_OPTZ(NNZ))             ;  T_OPTZ      =  20._SPP  !OPTIMAL TEMPERATURE FOR ZOOPLANKTON
      ALLOCATE (Z_D(KBV,NNZ))            ;  Z_D         =  0.0_SPP  !ZOOPANKTON MORTALITY
      ALLOCATE (Z_N(KBV,NNZ))            ;  Z_N         =  0.0_SPP  !ZOOPLANKTON RESPIRATION
      ALLOCATE (Z_0(NNZ))                ;  Z_0         =  0.0_SPP  !ZOOPLANKTO THRESHOLLD
!************  ALLOCATING                   BACTERIA       AND     INITIALIZATION   *************
      ALLOCATE (A_TB(NNB))               ;  A_TB        =  0.0_SPP  !T EXPONENTAIL COEFFICIENT
      ALLOCATE (BIO_B(KBV,NNB))          ;  BIO_B       =  0.0_SPP  !BACTERIA BIOMASS
      ALLOCATE (B_0(NNB))                ;  B_0         =  0.0_SPP  !BACTERIA THRESHOLD
      ALLOCATE (DELTA_B(NNB))            ;  DELTA_B     =  0.4_SPP  !RATIO BETWEEN NH4 AND DON UPTAKE
      ALLOCATE (EFFIBD(NND,NNB))         ;  EFFIBD      =  0.3_SPP  !GROWTH EFFICIENCY ON DETRITUS
      ALLOCATE (EFFIDOM(NNM,NNB))        ;  EFFIDOM     =  0.3_SPP  !GROWTH EFFICIENCY ON DOM
      ALLOCATE (EFFIN(NNN,NNB))          ;  EFFIN       =  1.0_SPP  !GROWTH EFFICIENCY ON NH4
      ALLOCATE (FVB(KBV,NNB))            ;  FVB         =  0.0_SPP  !BACTERIA SUBSTRACTED OF THRESHOLD
      ALLOCATE (R_B(NNB))                ;  R_B         =  1.E-6_SPP   !BACTERIA RESPIRATION
      ALLOCATE (RESPIRATION_B(KBV,NNB))  ;  RESPIRATION_B= 0.0_SPP  !BACTERIA RESPIRATION
      ALLOCATE (SIGMA_BD(NND,NNB))       ;  SIGMA_BD    =  0.2_SPP  !PREFERENCE ON DETRITUS
      ALLOCATE (SIGMA_DOM(NNM,NNB))      ;  SIGMA_DOM   =  0.2_SPP  !PREFERENCE ON DOM
      ALLOCATE (SIGMA_N(NNN,NNB))        ;  SIGMA_N     =  1.0_SPP  !PREFERENCE ON DON
      ALLOCATE (T_OPTB(NNB))             ;  T_OPTB      =  20._SPP  !OPTIMAL TEMPERATURE FOR BACTERIA
      ALLOCATE (U_B(NNB))                ;  U_B         =  0.0_SPP  !BACTERIA GROWTH RATE
      ALLOCATE (UBMAX(NNB))              ;  UBMAX       =  1.E-5_SPP !MAXIMUM GROWTH RATE OF BACTERIA
      ALLOCATE (UPTAKE_BD(KBV,NND,NNB))  ;  UPTAKE_BD   =  0.0_SPP  !DETRITUS UPTAKE BY BACTERI
      ALLOCATE (UPTAKE_BN(KBV,NNN,NNB))  ;  UPTAKE_BN   =  0.0_SPP  !NUTRIENT UPTAKE BY BACTERIA
      ALLOCATE (UPTAKE_DOM(KBV,NNM,NNB)) ;  UPTAKE_DOM  =  0.0_SPP  !DOM UPTAKE BY BACTERIA
      ALLOCATE (UPTAKE_DON(KBV,NNB))     ;  UPTAKE_DON  =  0.0_SPP  !DON UPTAKE BY BACTERIA

!******************  ALLOCATING             DETRITUS       AND     INITIALIZATION   *************************
      ALLOCATE (AGGREGATION(KBV,0:NND)) ;AGGREGATION    =  0.0_SPP  !AGGREGATION
      ALLOCATE (ALPHA_BD(NND,NNB,NNZ))   ;  ALPHA_BD    =  0.0_SPP  !FEEDING LOSS FROM BACTERIA TO DETRITUS
      ALLOCATE (ALPHA_DAG(NND))          ;  ALPHA_DAG   =  0.0_SPP  !AGGREGATION COEFFICIENT
      ALLOCATE (ALPHA_DD(NND,NND,NNZ))   ;  ALPHA_DD    =  0.0_SPP !FEEDING LOSS FROM DETRITUS TO DETRITUS
      ALLOCATE (ALPHA_DDG(NND))          ;  ALPHA_DDG   =  0.0_SPP  !DISAGGREGATION COEFFICIENT
      ALLOCATE (ALPHA_PD(NND,NNP,NNZ))   ;  ALPHA_PD    =  0.0_SPP  !FEEDING LOSS FROM PHYTO TO DETRITUS
      ALLOCATE (ALPHA_ZD(NND,NNZ,NNZ))   ;  ALPHA_ZD    =  0.0_SPP  !FEEDING LOSS FROM ZOO TO DETRITUS
      ALLOCATE (BIO_D(KBV,NND))          ;  BIO_D       =  0.0_SPP  !DETRITUS ABUNDANCE
      ALLOCATE (DETRITE(KBV,NND))        ;  DETRITE     =  0.0_DPP  !DETRITUS CONCENTRATION 
      ALLOCATE (D_D(NND))                ;  D_D         =  0.0_SPP  !DISSOLUTION COEFFICIENT
      ALLOCATE (D_0(NND))                ;  D_0         =  0.0_SPP  !DETRITUS THRESHOLD
      ALLOCATE (D_RN(NND))               ;  D_RN        =  0.0_SPP  !REMINERALIZATION COEFFICIENT
      ALLOCATE (DISAGGREGATION(KBV,NND+1)); DISAGGREGATION =  0.0_SPP !DISAGGREGATION
      ALLOCATE (D_SINK(KBV,NND))         ;  D_SINK      =  0.0_DPP  !DETRITUS SINKING FLUX
      ALLOCATE (D_SINK1(KBV))            ;  D_SINK1     =  0.0_DPP  !SINKING FLUX WITHIN 1 ITERATION
      ALLOCATE (DISSOLUTION(KBV,NND))    ;  DISSOLUTION =  0.0_SPP  !DISSOLUTION
      ALLOCATE (EPSILON_PD(NND,NNP))     ;  EPSILON_PD  =  0.0_SPP  !PHYTO MORTALITY TO DETRITUS
      ALLOCATE (EPSILON_ZD(NND,NNZ))     ;  EPSILON_ZD  =  0.0_SPP  !ZOO MORTALTIY TO DETRITUS
      ALLOCATE (FVD(KBV,NND))            ;  FVD         =  0.0_SPP  !DETRITUS AFTER THRESHOLD SUBSTRACTION
      ALLOCATE (REMINERALIZATION(KBV,NND)); REMINERALIZATION=0.0_SPP !REMINERALIZATION   
      ALLOCATE (W_D(NND))                ;  W_D         =  0.0_SPP  !DETRITUS SINKING VELOCITY
      ALLOCATE (W_D_TOTAL(NND))          ;  W_D_TOTAL   =  0.0_SPP  !SINKING DISTANCE WITIN 1 T_STEP

!******************  ALLOCATING             DOM            AND     INITIALIZATION  *****************
      ALLOCATE (AGE_DOM(0:NNM))          ;  AGE_DOM     =  0.0_SPP  !DOM AGEING
      ALLOCATE (ALPHA_DOM(NNM))          ;  ALPHA_DOM   =  0.0_SPP  !DOM AGEING COEFFICIENT
      ALLOCATE (ALPHA_PDOM(NNM,NNP))     ;  ALPHA_PDOM  =  0.0_SPP  !PHYTO EXUDATION TO DOM
      ALLOCATE (ALPHA_DDOM(NNM,NND))     ;  ALPHA_DDOM  =  0.0_SPP  !DETRITUS DISSOLUTION TO DOM
      ALLOCATE (ALPHA_ZPDOM(NNM,NNP,NNZ));  ALPHA_ZPDOM =  0.0_SPP  !FEEDING LOSS FROM PHYTO TO DOM
      ALLOCATE (ALPHA_ZDDOM(NNM,NNP,NNZ));  ALPHA_ZDDOM =  0.0_SPP  !FEEDING LOSS FROM DETRITUS TO DOM
      ALLOCATE (ALPHA_ZZDOM(NNM,NNZ,NNZ));  ALPHA_ZZDOM =  0.0_SPP  !FEEDING LOSS FROM ZOO TO DOM
      ALLOCATE (ALPHA_ZBDOM(NNM,NNB,NNZ));  ALPHA_ZBDOM =  0.0_SPP  !FEEDING LOSS FROM BACTERIA TO DOM
      ALLOCATE (BIO_DOM(KBV,NNM))        ;  BIO_DOM     =  0.0_SPP  !DOM CONCENTRATION
      ALLOCATE (DOM_0(NNM))              ;  DOM_0       =  0.0_SPP  !DOM THRESHOLD
      ALLOCATE (FVDOM(KBV,NNM))          ;  FVDOM       =  0.0_SPP  !DOM SUBSTRACTED BY THRESHOLD

!*************  ALLOCATING  NUTRIENT        AND                    INITIALIZATION  ******************
      ALLOCATE (BIO_N(KBV,NNN))          ;  BIO_N       =  0.0_SPP  !NUTRIENT CONCENTRATION
      ALLOCATE (FVNN(KBV,NNN))           ;  FVNN        =  0.0_SPP  !NUTRIENT BUSTRACTED BY THRESHOLD
      ALLOCATE (KSN(NNN,NNP))            ;  KSN         =  0.0_SPP  !HALF-SATURAITON OF NUTRIENT
                                         ;  L_NH4N      =  0.0_SPP  !CRITICL LIGHT FOR NITRIFICATION
      ALLOCATE (NH4_NO3(KBV))            ;  NH4_NO3     =  0.0_SPP  !NITRIFICATION
      ALLOCATE (N2CB(NNN,NNB))           ;  N2CB        =  1.0_SPP  !ELEMENT RATIO IN BACTERIA
      ALLOCATE (N2CD(NNN,NND))           ;  N2CD        =  1.0_SPP  !ELEMENT RATIO IN DETRITUS
      ALLOCATE (N2CP(NNN,NNP))           ;  N2CP        =  1.0_SPP  !ELEMENT RATIO IN PHYTO
      ALLOCATE (N2CZ(NNN,NNZ))           ;  N2CZ        =  1.0_SPP  !ELEMENT RATIO IN ZOOPLANKTON
      ALLOCATE (N2CDOM(NNN,NNM))         ;  N2CDOM      =  1.0_SPP  !ELEMENT RATIO IN DOM
      ALLOCATE (N2CZB(NNN,NNB,NNZ))      ;  N2CZB       =  1.0_SPP  !ELEM. RATIO IN ACTIVE RESPIR. FROM GRAZING ON BACTERIA
      ALLOCATE (N2CZD(NNN,NND,NNZ))      ;  N2CZD       =  1.0_SPP  !ELEM. RATIO IN ACTIVE RESPIR. FROM GRAZING ON DETRITUS
      ALLOCATE (N2CZP(NNN,NNP,NNZ))      ;  N2CZP       =  1.0_SPP  !ELEM. RATIO IN ACTIVE RESPIR. FROM GRAZING ON PHYTO
      ALLOCATE (N2CZZ(NNN,NNZ,NNZ))      ;  N2CZZ       =  1.0_SPP  !ELEM. RATIO IN ACTIVE RESPIR. FROM PREDATION ON ZOO
      ALLOCATE (N_0(NNN))                ;  N_0         =  0.0_SPP  !NUTRIENT THRESHOLD
                                         ;  R_AN        =  0.0_SPP  !NITRIFICATION COEFFICIENT
      ALLOCATE (UPTAKE_PN(KBV,NNN,NNP))  ;  UPTAKE_PN   =  0.0_SPP  !NUTRIENT UPTAKE BY PHYTO

!*************************     FUNCTIONS AND PARAMETER VALUES    ***************
 
!*************************  PHYTOPLANKTON     PARAMETERS  ****************************!
      !read in light function
      Call Get_Val(L_FUNCTION,biofile,'L_FUNCTION',line=linenum,echo=.true.) 

      !read in zooplankton grazing function
      Call Get_Val(G_FUNCTION,biofile,'G_FUNCTION',line=linenum,echo=.true.) 
 
      !read in LIGHT FUNCTION COEFFICIENT
      Call Get_Val_Array(ALPHAP,biofile,'ALPHAP',NNP,echo=.true.)

      !read in COEFFICIENT BETWEEN MULTIPLICATION & MINIMUM
      Call Get_Val_Array(ALPHA_U,biofile,'ALPHA_U',NNP,echo=.true.)

      !read in EXPONENTIAL FOR T FORCING
      Call Get_Val_Array(A_TP,biofile,'A_TP',NNP,echo=.true.,option=.true.)

      !read in ATTENUATION COEFFICIENT OF CHLOROPHYLL
      Call Get_Val(ATANU_C,biofile,'ATANU_C',line=linenum,echo=.true.) 

      !read in ATTENUATION COEFFICIENT OF DETRITUS
      Call Get_Val(ATANU_D,biofile,'ATANU_D',line=linenum,echo=.true.) 

      !read in ATTENUATION COEFFICIENT OF WATER
      Call Get_Val(ATANU_W,biofile,'ATANU_W',line=linenum,echo=.true.) 

      !read in PHOTOINHIBITION COEFFICIENT
      Call Get_Val_Array(BETAP,biofile,'BETAP',NNP,echo=.true.)

      !read in CHLOROPHYLL TO CARBON RATIO IN PHYTO
      Call Get_Val_Array(CHL2C,biofile,'CHL2C',NNP,echo=.true.)

      !read in ACTIVE EXUDATION COEFFICIENT
      Call Get_Val_Array(D_DOM,biofile,'D_DOM',NNP,echo=.true.)

      !read in PASSIVE PHYTO EXUDATION COEFFICIENT
      Call Get_Val_Array(DPDOM,biofile,'DPDOM',NNP,echo=.true.)

      !read in OPTIMAL LIGHT FOR PHYTO
      Call Get_Val_Array(I_OPT,biofile,'I_OPT',NNP,echo=.true.)

      !read in HALF_SATURATION OF LIGHT FUNCTIONS
      Call Get_Val_Array(K_LIGHT,biofile,'K_LIGHT',NNP,echo=.true.)

      !read in PHYTO MORTALITY EOFFICIENT
      Call Get_Val_Array(MPD,biofile,'MPD',NNP,echo=.true.)

      !read in PHYTO MORTALITY ORDER
      Call Get_Val_Array(M_P,biofile,'M_P',NNP,echo=.true.)

      !read in POWER OF LIGHT FUNCTIONS
      Call Get_Val_Array(N_P,biofile,'N_P',NNP,echo=.true.)

      !read in PAR FRACTION
      Call Get_Val(PARFRAC,biofile,'PARFRAC',line=linenum,echo=.true.) 

      !read in PHYTOPANKTON THRESHOLDS
      Call Get_Val_Array(P_0,biofile,'P_0',NNP,echo=.true.,option=.true.)

      !read in EXPONENTIAL FOR T FORCING ON RESPIRATION
      Call Get_Val(RP_T,biofile,'RP_T',line=linenum,echo=.true.,option=.true.) 

      !read in PHYTO RESPIRATION COEFFICIENT
      Call Get_Val_Array(R_P,biofile,'R_P',NNP,echo=.true.)

      !read in OPTIMAL TEMPERATURE FOR PHYTO
      Call Get_Val_Array(T_OPTP,biofile,'T_OPTP',NNP,echo=.true.,option=.true.)

      !read in PHYTO MAXIMUM GROWTH RATE
      Call Get_Val_Array(UMAX,biofile,'UMAX',NNP,echo=.true.)

      !read in SINKING VELOCITY OF PHYTO
      Call Get_Val_Array(W_P,biofile,'W_P',NNP,echo=.true.)


!********************  ZOOPLANKTON     PARAMETERS  **********************!

      !read in ACTIVE RESPIRATION COEFFICIENT
      Call Get_Val_Array(ACTIVE_R,biofile,'ACTIVE_R',NNZ,echo=.true.)

      !read in EXPONENTIAL FOR T FORCING
      Call Get_Val_Array(A_TZ,biofile,'A_TZ',NNZ,echo=.true.,option=.true.)


      !read in GROWTH EFFICIENCY ON BACTERIA
      if(NNB*NNZ/=0)then
        allocate(tmp(NNB*NNZ))      ;  tmp = 0.0
        Call Get_Val_Array(tmp,biofile,'EFFIB',NNB*NNZ,echo=.true.)
        do J=1,NNZ
          EFFIB(1:NNB,J)=tmp(1+(J-1)*NNB:NNB+(J-1)*NNB)
        end do
        deallocate(tmp)
      endif 

      !read in GROWTH EFFICIENCY ON DETRITUS
      if(NND*NNZ/=0)then
        allocate(tmp(NND*NNZ))      ;  tmp = 0.0
        Call Get_Val_Array(tmp,biofile,'EFFID',NND*NNZ,echo=.true.)
        do J=1,NNZ
          EFFID(1:NND,J)=tmp(1+(J-1)*NND:NND+(J-1)*NND)
        end do
        deallocate(tmp)
      endif

      !read in GROWTH EFFICIENCY ON PHYTO
      if(NNP*NNZ/=0)then
        allocate(tmp(NNP*NNZ))      ;  tmp = 0.0
        Call Get_Val_Array(tmp,biofile,'EFFIP',NNP*NNZ,echo=.true.)
        do J=1,NNZ
           EFFIP(1:NNP,J)=tmp(1+(J-1)*NNP:NNP+(J-1)*NNP)
        end do
        deallocate(tmp)
      endif

      !read in GROWTH EFFICIENCY ON ZOOPLANKTON
      if(NNZ*NNZ/=0)then
        allocate(tmp(NNZ*NNZ))      ;  tmp = 0.0
        Call Get_Val_Array(tmp,biofile,'EFFIZ',NNZ*NNZ,echo=.true.)
        do J=1,NNZ
          EFFIZ(1:NNZ,J)=tmp(1+(J-1)*NNZ:NNZ+(J-1)*NNZ)
        end do
        deallocate(tmp)
      endif


      !read in MAXIMUM GRAZING RATE
      Call Get_Val_Array(G_MAX,biofile,'G_MAX',NNZ,echo=.true.)

      !read in HALF-SATURATION FOR ZOO GRAZING
      Call Get_Val_Array(K_ZG,biofile,'K_ZG',NNZ,echo=.true.)

      !read in ZOOPLANKTON GRAZING POWER
      Call Get_Val_Array(M_G,biofile,'M_G',NNZ,echo=.true.)

      !read in ZOOPLANKTON MORTALITY COEFFICIENT
      Call Get_Val_Array(MZD,biofile,'MZD',NNZ,echo=.true.)

      !read in ZOOPLANKTON MORTALITY POWER
      Call Get_Val_Array(M_Z,biofile,'M_Z',NNZ,echo=.true.,option=.true.)

      !read in GRAZING THRESHOLD
      Call Get_Val_Array(P_C,biofile,'P_C',NNZ,echo=.true.)

      !read in RECRUIT COEFFICIENT
      Call Get_Val_Array(R_RECRUIT,biofile,'R_RECRUIT',NNZ,echo=.true.)

      !read in ZOO PASSIVE RESPIRATION COEFFICIENT
      Call Get_Val_Array(R_Z,biofile,'R_Z',NNZ,echo=.true.)

      !read in PREFERENCE ON BACTERIA
      if(NNB*NNZ)then
        allocate(tmp(NNB*NNZ)) ;  tmp = 0.0
        Call Get_Val_Array(tmp,biofile,'SIGMA_B',NNB*NNZ,echo=.true.)
        do j=1,nnz
           SIGMA_B(1:nnb,j)=tmp(1+(j-1)*nnb:nnb+(j-1)*nnb)
        end do
        deallocate(tmp)
      endif

      !read in PREFERENCE ON DETRITUS
      if(NNZ*NND/=0)then
        allocate(tmp(NNZ*NND)) ;  tmp = 0.0
        Call Get_Val_Array(tmp,biofile,'SIGMA_D',NNZ*NND,echo=.true.)
        do j=1,nnz
           SIGMA_D(1:nnd,j)=tmp(1+(j-1)*nnd:nnd+(j-1)*nnd)
        end do
        deallocate(tmp)
      endif

      !read in PREFERENCE ON PHYTO
      if(NNZ*NNP/=0)then
        allocate(tmp(NNZ*NNP)) ;  tmp = 0.0
        Call Get_Val_Array(tmp,biofile,'SIGMA_P',NNZ*NNP,echo=.true.)
        do j=1,nnz
           SIGMA_P(1:nnp,j)=tmp(1+(j-1)*nnp:nnp+(j-1)*nnp)
        end do
        deallocate(tmp)
      endif

      !read in PREFERENCE ON ZOOPLANKTON
      if(NNZ*NNZ)then
        allocate(tmp(NNZ*NNZ)) ;  tmp = 0.0
        Call Get_Val_Array(tmp,biofile,'SIGMA_Z',NNZ*NNZ,echo=.true.)
        do j=1,nnz
           SIGMA_Z(1:nnz,j)=tmp(1+(j-1)*nnz:nnz+(j-1)*nnz)
        end do
        deallocate(tmp)
      endif

      !read in OPTIMAL TEMPERATURE FOR ZOOPLANKTON
      Call Get_Val_Array(T_OPTZ,biofile,'T_OPTZ',NNZ,echo=.true.,option=.true.)

      !read in ZOOPLANKTON THRESHOLLD
      Call Get_Val_Array(Z_0,biofile,'Z_0',NNZ,echo=.true.,option=.true.)


!*********************  BACTERIA        PARAMETERS  **********************!

      !read in EXPONENTAIL FOR T FORCING
      Call Get_Val_Array(A_TB,biofile,'A_TB',NNB,echo=.true.)

      !read in BACTERIA THRESHOLLD
      Call Get_Val_Array(B_0,biofile,'B_0',NNB,echo=.true.,option=.true.)

      !read in RATIO BETWEEN NH4 AND DON UPTAKE
      Call Get_Val_Array(DELTA_B,biofile,'DELTA_B',NNB,echo=.true.)

      !read in GROWTH EFFICIENCY ON DETRITUS
      if(NNB*NND/=0)then
        allocate(tmp(NNB*NND)) ;  tmp = 0.0
        Call Get_Val_Array(tmp,biofile,'EFFIBD',NNB*NND,echo=.true.)
        do j=1,nnb
           EFFIBD(1:nnd,j)=tmp(1+(j-1)*nnd:nnd+(j-1)*nnd)
        end do
        deallocate(tmp)
      endif

      !read in GROWTH EFFICIENCY ON DOM
      if(NNB*NNM/=0)then
        allocate(tmp(NNB*NNM)) ;  tmp = 0.0
        Call Get_Val_Array(tmp,biofile,'EFFIDOM',NNB*NNM,echo=.true.)
        do j=1,nnb
           EFFIDOM(1:nnm,j)=tmp(1+(j-1)*nnm:nnm+(j-1)*nnm)
        end do
        deallocate(tmp)
      endif

      !read in GROWTH EFFICIENCY ON NH4
      if(NNB*NNN/=0)then
        allocate(tmp(NNB*NNN)) ;  tmp = 0.0
        Call Get_Val_Array(tmp,biofile,'EFFIN',NNB*NNN,echo=.true.)
        do j=1,nnb
           EFFIN(1:nnn,j)=tmp(1+(j-1)*nnn:nnn+(j-1)*nnn)
        end do
        deallocate(tmp)
      endif

      !read in BACTERIA RESPIRATION COEFFICIENT
      Call Get_Val_Array(R_B,biofile,'R_B',NNB,echo=.true.)

      !read in PREFERENCE ON DOM
      if(NNB*NND/=0)then
        allocate(tmp(NNB*NND)) ;  tmp = 0.0
        Call Get_Val_Array(tmp,biofile,'SIGMA_BD',NNB*NND,echo=.true.)
        do j=1,nnb
           SIGMA_BD(1:nnd,j)=tmp(1+(j-1)*nnd:nnd+(j-1)*nnd)
        end do
        deallocate(tmp)
      endif

      !read in PREFERENCE ON DOM
      if(NNB*NNM/=0)then
        allocate(tmp(NNB*NNM)) ;  tmp = 0.0
        Call Get_Val_Array(tmp,biofile,'SIGMA_DOM',NNB*NNM,echo=.true.)
        do j=1,nnb
           SIGMA_DOM(1:nnm,j)=tmp(1+(j-1)*nnm:nnm+(j-1)*nnm)
        end do
        deallocate(tmp)
      endif

      !read in PREFERENCE ON DON
      if(NNB*NNN/=0)then
        allocate(tmp(NNB*NNN)) ;  tmp = 0.0
        Call Get_Val_Array(tmp,biofile,'SIGMA_N',NNB*NNN,echo=.true.)
        do j=1,nnb
           SIGMA_N(1:nnn,j)=tmp(1+(j-1)*nnn:nnn+(j-1)*nnn)
        end do
        deallocate(tmp)
      endif

      !read in OPTIMAL T FOR BACTERIA
      Call Get_Val_Array(T_OPTB,biofile,'T_OPTB',NNB,echo=.true.)

      !read in MAXIMUM GROWTH RATE OF BACTERIA
      Call Get_Val_Array(UBMAX,biofile,'UBMAX',NNB,echo=.true.)



!***********************  DETRITUS    PARAMETERS  **********************

      !read in FEEDING LOSS FROM BACTERIA TO DETRITUS
      if(NND*NNB*NNZ/=0)then
        allocate(tmp(NND*NNB*NNZ)) ;  tmp = 0.0
        Call Get_Val_Array(tmp,biofile,'ALPHA_BD',NND*NNB*NNZ,echo=.true.)
        do K=1,NNZ
          do J=1,NNB
            ALPHA_BD(1:NND,J,K)=tmp(1+(J-1)*NND*(K-1):NND+(J-1)*NND*(K-1))
          end do
        end do
        deallocate(tmp)
      endif

      !read in AGGREGATION COEFFICIENT
      Call Get_Val_Array(ALPHA_DAG,biofile,'ALPHA_DAG',NND,echo=.true.,option=.true.)

      !read in FEEDING LOSS FROM DETRITUS TO DETRITUS
      if(NND*NND*NNZ/=0)then
        allocate(tmp(NND*NND*NNZ)) ;  tmp = 0.0
        Call Get_Val_Array(tmp,biofile,'ALPHA_DD',NND*NND*NNZ,echo=.true.)
        do K=1,NNZ
          do J=1,NND
            ALPHA_DD(1:NND,J,K)=tmp(1+(J-1)*NND*(K-1):NND+(J-1)*NND*(K-1))
          end do
        end do
        deallocate(tmp)
      endif

      !read in DISAGGREGATION COEFFICIENT
      Call Get_Val_Array(ALPHA_DDG,biofile,'ALPHA_DDG',NND,echo=.true.,option=.true.)

      !read in FEEDING LOSS FROM PHYTO TO DETRITU
      if(NND*NNP*NNZ/=0)then
        allocate(tmp(NND*NNP*NNZ)) ;  tmp = 0.0
        Call Get_Val_Array(tmp,biofile,'ALPHA_PD',NND*NNP*NNZ,echo=.true.)
        do K=1,NNZ
          do J=1,NNP
            ALPHA_PD(1:NND,J,K)=tmp(1+(J-1)*NND*(K-1):NND+(J-1)*NND*(K-1))
          end do
        end do
        deallocate(tmp)
      endif

      !read in FEEDING LOSS FROM ZOO TO DETRITUS
      if(NND*NNZ*NNZ/=0)then
        allocate(tmp(NND*NNZ*NNZ)) ;  tmp = 0.0
        Call Get_Val_Array(tmp,biofile,'ALPHA_ZD',NND*NNZ*NNZ,echo=.true.)
        do K=1,NNZ
          do J=1,NNZ
            ALPHA_ZD(1:NND,J,K)=tmp(1+(J-1)*NND*(K-1):NND+(J-1)*NND*(K-1))
          end do
        end do
        deallocate(tmp)
      endif

      !read in DISSOLUTION COEFFICIENT
      Call Get_Val_Array(D_D,biofile,'D_D',NND,echo=.true.)

      !read in DETRITUS THRESHOLD
      Call Get_Val_Array(D_0,biofile,'D_0',NND,echo=.true.,option=.true.)

      !read in REMINERALIZATION COEFFICIENT
      Call Get_Val_Array(D_RN,biofile,'D_RN',NND,echo=.true.)

      !read in PHYTO MORTALITY TO DETRITUS
      if(NND*NNP/=0)then
        allocate(tmp(NND*NNP)) ;  tmp = 0.0
        Call Get_Val_Array(tmp,biofile,'EPSILON_PD',NND*NNP,echo=.true.)
        do j=1,NNP
           EPSILON_PD(1:NND,j)=tmp(1+(J-1)*NND:NND+(J-1)*NND)
        end do
        deallocate(tmp)
      endif

      !read in ZOO MORTALTIY TO DETRITU
      if(NND*NNZ/=0)then
        allocate(tmp(NND*NNZ)) ;  tmp = 0.0
        Call Get_Val_Array(tmp,biofile,'EPSILON_ZD',NND*NNZ,echo=.true.)
        do j=1,NNZ
           EPSILON_ZD(1:NND,j)=tmp(1+(J-1)*NND:NND+(J-1)*NND)
        end do
        deallocate(tmp)
      endif

      !read in DETRITUS SINKING VELOCITY
      Call Get_Val_Array(W_D,biofile,'W_D',NND,echo=.true.)

!**********************  DOM      PARAMETERS  *************************!

      !read in DOM AGEING COEFFICIENT
      Call Get_Val_Array(ALPHA_DOM,biofile,'ALPHA_DOM',NNM,echo=.true.)

      !read in PHYTO EXUDATION TO DOM
      if(NNM*NNP/=0)then
        allocate(tmp(NNM*NNP)) ;  tmp = 0.0
        Call Get_Val_Array(tmp,biofile,'ALPHA_PDOM',NNM*NNP,echo=.true.)
        do j=1,NNP
           ALPHA_PDOM(1:NNM,j)=tmp(1+(J-1)*NNM:NNM+(J-1)*NNM)
        end do
        deallocate(tmp)
      endif

      !read in DETRITUS DISSOLUTION TO DOM
      if(NNM*NND/=0)then
        allocate(tmp(NNM*NND)) ;  tmp = 0.0
        Call Get_Val_Array(tmp,biofile,'ALPHA_DDOM',NNM*NND,echo=.true.)
        do j=1,NND
           ALPHA_DDOM(1:NNM,j)=tmp(1+(J-1)*NNM:NNM+(J-1)*NNM)
        end do
        deallocate(tmp)
      endif

      !read in FEEDING LOSS FROM PHYTO TO DOM
      if(NNM*NNP*NNZ/=0)then
        allocate(tmp(NNM*NNP*NNZ)) ;  tmp = 0.0
        Call Get_Val_Array(tmp,biofile,'ALPHA_ZPDOM',NNM*NNP*NNZ,echo=.true.)
        do K=1,NNZ
          do J=1,NNP
            ALPHA_ZPDOM(1:NNM,J,K)=tmp(1+(J-1)*NNM*(K-1):NNM+(J-1)*NNM*(K-1))
          end do
        end do
        deallocate(tmp)
      endif

      !read in FEEDING LOSS FROM DETRITUS TO DOM
      if(NNM*NND*NNZ/=0)then
        allocate(tmp(NNM*NND*NNZ)) ;  tmp = 0.0
        Call Get_Val_Array(tmp,biofile,'ALPHA_ZDDOM',NNM*NND*NNZ,echo=.true.)
        do K=1,NNZ
          do J=1,NND
            ALPHA_ZDDOM(1:NNM,J,K)=tmp(1+(J-1)*NNM*(K-1):NNM+(J-1)*NNM*(K-1))
          end do
        end do
        deallocate(tmp)
      endif

      !read in FEEDING LOSS FROM ZOO TO DOM
      if(NNM*NNZ*NNZ/=0)then
        allocate(tmp(NNM*NNZ*NNZ)) ;  tmp = 0.0
        Call Get_Val_Array(tmp,biofile,'ALPHA_ZZDOM',NNM*NNZ*NNZ,echo=.true.)
        do K=1,NNZ
          do J=1,NNZ
             ALPHA_ZZDOM(1:NNM,J,K)=tmp(1+(J-1)*NNM*(K-1):NNM+(J-1)*NNM*(K-1))
          end do
        end do
        deallocate(tmp)
      endif

      !read in FEEDING LOSS FROM BACTERIA TO DOM
      if(NNM*NNB*NNZ/=0)then
        allocate(tmp(NNM*NNB*NNZ)) ;  tmp = 0.0
        Call Get_Val_Array(tmp,biofile,'ALPHA_ZZDOM',NNM*NNB*NNZ,echo=.true.)
        do K=1,NNZ
          do J=1,NNB
             ALPHA_ZZDOM(1:NNM,J,K)=tmp(1+(J-1)*NNM*(K-1):NNM+(J-1)*NNM*(K-1))
          end do
        end do
        deallocate(tmp)
      endif

      !read in DOM THRESHOLLD
      Call Get_Val_Array(DOM_0,biofile,'DOM_0',NNM,echo=.true.,option=.true.)


!*****************  NUTRIENT         PARAMETERS *******************!

      !read in HALF-SATURATION OF NUTRIENT
      if(NNN*NNP/=0)then
        allocate(tmp(NNN*NNP)) ;  tmp = 0.0
        Call Get_Val_Array(tmp,biofile,'KSN',NNN*NNP,echo=.true.)
        do j=1,NNP
           KSN(1:NNN,j)=tmp(1+(J-1)*NNN:NNN+(J-1)*NNN)
        end do
        deallocate(tmp)
      endif

      !read in ELEMENT RATIO IN BACTERIA
      if(NNN*NNB/=0)then
        allocate(tmp(NNN*NNB)) ;  tmp = 0.0
        Call Get_Val_Array(tmp,biofile,'N2CB',NNN*NNB,echo=.true.)
        do j=1,NNB
           N2CB(1:NNN,j)=tmp(1+(J-1)*NNN:NNN+(J-1)*NNN)
        end do
        deallocate(tmp)
      endif

      !read in ELEMENT RATIO IN DETRITUS
      if(NNN*NND/=0)then
        allocate(tmp(NNN*NND)) ;  tmp = 0.0
        Call Get_Val_Array(tmp,biofile,'N2CD',NNN*NND,echo=.true.)
        do j=1,NND
           N2CD(1:NNN,j)=tmp(1+(J-1)*NNN:NNN+(J-1)*NNN)
        end do
        deallocate(tmp)
      endif

      !read in ELEMENT RATIO IN PHYTO
      if(NNN*NNP/=0)then
        allocate(tmp(NNN*NNP)) ;  tmp = 0.0
        Call Get_Val_Array(tmp,biofile,'N2CP',NNN*NNP,echo=.true.)
        do j=1,NNP
           N2CP(1:NNN,j)=tmp(1+(J-1)*NNN:NNN+(J-1)*NNN)
        end do
        deallocate(tmp)
      endif

      !read in ELEMENT RATIO IN ZOOPLANKTON
      if(NNN*NNZ/=0)then
        allocate(tmp(NNN*NNZ)) ;  tmp = 0.0
        Call Get_Val_Array(tmp,biofile,'N2CZ',NNN*NNZ,echo=.true.)
        do j=1,NNZ
           N2CZ(1:NNN,j)=tmp(1+(J-1)*NNN:NNN+(J-1)*NNN)
        end do
        deallocate(tmp)
      endif

      !read in ELEMENT RATIO IN DOM
      if(NNN*NNM/=0)then
        allocate(tmp(NNN*NNM)) ;  tmp = 0.0
        Call Get_Val_Array(tmp,biofile,'N2CDOM',NNN*NNM,echo=.true.)
        do j=1,NNM
           N2CDOM(1:NNN,j)=tmp(1+(J-1)*NNN:NNN+(J-1)*NNN)
        end do
        deallocate(tmp)
      endif

      !read in NUTRIENT THRESHOLD
      Call Get_Val_Array(N_0,biofile,'N_0',NNM,echo=.true.,option=.true.)

      !read in NITRIFICATION COEFFICIENT
      Call Get_Val(R_AN,biofile,'R_AN',line=linenum,echo=.true.) 

!*******************************************************************!
! *********   CHECKING GROWTH EFFICIENCY AND FEEDLING LOSSES    ****!
!*******************************************************************!
       DO I=1,NNZ
         DO J=1,NNP                                    !FOR GRAZING ON PHYTO
           TOTALL=ACTIVE_R(I)+EFFIP(J,I)
	   print*,totall,active_r(i),effip(j,i),nnp,nnz,nnd

           IF(NND.GT.1) THEN
              IF (DETRITUS_NAME(2)=='Detritus_Si') THEN
                 DO K=1,1
                   TOTALL=TOTALL+ALPHA_PD(K,J,I)
                 END DO
              ELSE
                 DO K=1,NND
                   TOTALL=TOTALL+ALPHA_PD(K,J,I)
                 END DO
              END IF
           ELSE
               DO K=1,NND
                  TOTALL=TOTALL+ALPHA_PD(K,J,I)
               END DO
           ENDIF
	   print*,totall,alpha_pd(1,j,i),nnm


           DO K=1,NNM
              TOTALL=TOTALL+ALPHA_ZPDOM(K,J,I)
           END DO
           IF(TOTALL .NE. 1.) THEN
             PRINT*,'EFFICIENCY+LOSS OF ZOO ',I,' GRAZING ON PHYTO ',J,' UNEQUAL TO 1',TOTALL
             STOP
           END IF
         END DO

         DO J=1,NND                                    !FOR GRAZING ON DETRITUS
           TOTALL=ACTIVE_R(I)+EFFID(J,I)
           IF(NND.GT.1) THEN
             IF (DETRITUS_NAME(2)=='Detritus_Si') THEN
               DO K=1,1
                  TOTALL=TOTALL+ALPHA_DD(K,J,I)
               END DO
             ELSE
               DO K=1,NND
                  TOTALL=TOTALL+ALPHA_DD(K,J,I)
               END DO
             END IF
           ENDIF
           DO K=1,NNM
              TOTALL=TOTALL+ALPHA_ZDDOM(K,J,I)
           END DO
           IF(TOTALL > 1.) THEN
             PRINT*,'EFFICIENCY+LOSS OF ZOO ',I,' GRAZING ON DETRITUS ',J,' UNEQUAL TO 1',TOTALL
             STOP
           END IF
         END DO

         DO J=1,NNZ                                    !FOR PREDATION ON ZOO
           TOTALL=ACTIVE_R(I)+EFFIZ(J,I)
           IF(NND.GT.1) THEN
             IF (DETRITUS_NAME(2)=='Detritus_Si') THEN
               DO K=1,1
                 TOTALL=TOTALL+ALPHA_ZD(K,J,I)
               END DO
             ELSE
               DO K=1,NND
                 TOTALL=TOTALL+ALPHA_ZD(K,J,I)
               END DO
             END IF
           ENDIF
           DO K=1,NNM
              TOTALL=TOTALL+ALPHA_ZZDOM(K,J,I)
           END DO
           IF(TOTALL > 1.) THEN
             PRINT*,'EFFICIENCY+LOSS OF ZOO ',I,' PREDATION ON ZOO ',J,' UNEQUAL TO 1',TOTALL
             STOP
           END IF
         END DO

         DO J=1,NNB                                    !FOR PREDATION ON BACTERIA
           TOTALL=ACTIVE_R(I)+EFFIB(J,I)
           IF(NND.GT.1) THEN
             IF (DETRITUS_NAME(2)=='Detritus_Si') THEN
               DO K=1,1
                 TOTALL=TOTALL+ALPHA_PD(K,J,I)
               END DO
             ELSE
               DO K=1,NND
                 TOTALL=TOTALL+ALPHA_BD(K,J,I)
               END DO
             END IF
           ENDIF
           DO K=1,NNM
              TOTALL=TOTALL+ALPHA_ZBDOM(K,J,I)
           END DO
           IF(TOTALL > 1.) THEN
             PRINT*,'EFFICIENCY+LOSS OF ZOO ',I,' PREDATION ON BAC ',J,' UNEQUAL TO 1',TOTALL
             STOP
           END IF
         END DO
         TOTALL=0.
         DO J=1,NND                                   !ZOO MORTALITY TO DETRITUS
           TOTALL=TOTALL+EPSILON_ZD(J,I)
         END DO
           IF(TOTALL > 1.) THEN
             PRINT*,'MORTALITY OF ZOO ',I,' TO DETRITUS DOES NOT EQUAL TO 1',TOTALL
             STOP
           END IF
       END DO !I=1,NNZ
       DO I=1,NNP
         TOTALL=0.
         IF(NND.GT.1) THEN
           IF (DETRITUS_NAME(2)=='Detritus_Si') THEN
             DO J=1,1                                   !PHYTO MORTALITY TO DETRITUS
               TOTALL=TOTALL+EPSILON_PD(J,I)
             END DO
           ELSE
             DO J=1,NND                                   !PHYTO MORTALITY TO DETRITUS
               TOTALL=TOTALL+EPSILON_PD(J,I)
             END DO
           END IF
         ENDIF
           IF(TOTALL > 1.) THEN
             PRINT*,'MORTALITY OF PHYTO ',I,' TO DETRITUS DOES NOT EQUAL TO 1',TOTALL
             STOP
           END IF
       END DO !I=1,NNP
! 99    continue
       GOTO 88  !MOVED TO BIO_INITIAL TO PREVENT MULTIPLE PRINT
!*******************************************************************!
! *********     PRINT OUT MODEL SETUP AND PARAMETER VALUES      ****!
!*******************************************************************!
       PRINT*
       PRINT*,'*****************************************************'
       PRINT*,'**  STRUCTURE AND FUNCTION OF THE BIOLOGICAL MODEL **'
       PRINT*,'*****************************************************'

          PRINT*
          PRINT*,'MODEL STRUCTURE        : ',  MODEL_STRUCTURE
         DO I=1,NNN
          PRINT*,'                         ',  NUTRIENT_NAME(I)
         END DO
         DO I=1,NNP
          PRINT*,'                         ',  PHYTO_NAME(I)
         END DO
         DO I=1,NNZ
          PRINT*,'                         ',  ZOO_NAME(I)
         END DO
         DO I=1,NND
          PRINT*,'                         ',  DETRITUS_NAME(I)
         END DO
         DO I=1,NNM
          PRINT*,'                         ', DOM_NAME(I)
         END DO
         DO I=1,NNB
          PRINT*,'                         ', BACTERIA_NAME(I)
         END DO
          WRITE(*,'(A26,A20)')' LIGHT FUNCTION         : ', L_FUNCTION
          WRITE(*,'(A26,A20)')' GRAZING FUNCTION       : ', G_FUNCTION
          PRINT*
      PRINT*,'*********    PHYTOPLANTON PARAMETERS    **************'
        PRINT*
        IF(L_FUNCTION.NE.'EXP_LIGHT'.AND.L_FUNCTION.NE.'SL62_LIGHT') THEN
          PRINT*,'ALPHA                  : '  , (ALPHAP(I),I=1,NNP)
        END IF
          PRINT*,'L_N COMBINE            : ' , (ALPHA_U(I),I=1,NNP)
          PRINT*,'T FORCING              : ' , (A_TP(I),I=1,NNP)
          PRINT*,'CHL ATTANUATION        : ' , ATANU_C
          PRINT*,'D ATTANUATION          : ', ATANU_D
          PRINT*,'WATER ATTANUATION      : ', ATANU_W
        IF(L_FUNCTION.EQ.'PGH80_LIGHT'.OR.L_FUNCTION.EQ.'V65_LIGHT'.OR. &
           L_FUNCTION.EQ.'BWDC9_LIGHT') THEN       
          PRINT*,'BETAP                  : ', (BETAP(I),I=1,NNP)
        END IF
          PRINT*,'CHL:C                  : ',(CHL2C(I),I=1,NNP)
          PRINT*,'ACTIVE DOM EXUD.       : ',(D_DOM(I),I=1,NNP)
          PRINT*,'PASSIVE DOM EXUD.      : ', (DPDOM(I),I=1,NNP)
        IF(L_FUNCTION.EQ.'SL62_LIGHT'.OR.L_FUNCTION.EQ.'V65_LIGHT'.OR. &
           L_FUNCTION.EQ.'PE78_LIGHT') THEN       
          PRINT*,'OPTIMAL LIGHT          : ',(I_OPT(I),I=1,NNP)
        END IF
        IF(L_FUNCTION.EQ.'MM_LIGHT'.OR.L_FUNCTION.EQ.'LB_LIGHT') THEN 
          PRINT*,'LIGHT HALF SATURATTION : ', (K_LIGHT(I),I=1,NNP)
        END IF
          PRINT*,'MORTALITY              : ', (MPD(I),I=1,NNP)
          PRINT*,'MORTALITY POWER        : ', (M_P(I),I=1,NNP)
        IF(L_FUNCTION.EQ.'LB_LIGHT'.OR.L_FUNCTION.EQ.'V65_LIGHT') THEN 
          PRINT*,'POWER OF LIGHT         : ', (N_P(I),I=1,NNP)
        END IF
          PRINT*,'THRESHOLD              : ', (P_0(I),I=1,NNP)
          PRINT*,'T ON RESPIRATION       : ', RP_T
          PRINT*,'RESPIRATION            : ', (R_P(I),I=1,NNP)
          PRINT*,'OPTIMAL T              : ', (T_OPTP(I),I=1,NNP)
          PRINT*,'MAXIMUM GROWTH         : ', (UMAX(I),I=1,NNP)
          PRINT*,'SINKING VELOCITY       : ', (W_P(I),I=1,NNP)
          PRINT*

          PRINT*,'********    ZOOPLANKTON PARAMETERS   ****************'
          PRINT*
          PRINT*,'ACTIVE RESPIRATION     : ', (ACTIVE_R(I),I=1,NNZ)
          PRINT*,'T FORCING EXPONENTIAL  : ', (A_TZ(I),I=1,NNZ)
        IF (NNB.GE.1) THEN
          PRINT*,'EFFICIENCY ON BACTERIA : ',((EFFIB(I,J),I=1,NNB),J=1,NNZ)
        END IF
        IF (NND.GE.1) THEN
          PRINT*,'EFFICIENCY ON DETRITUS : ',((EFFID(I,J),I=1,NND),J=1,NNZ)
        END IF
          PRINT*,'EFFICIENCY ON PHYTO    : ',((EFFIP(I,J),I=1,NNP),J=1,NNZ)
          PRINT*,'EFFICIENCY ON ZOO      : ',((EFFIZ(I,J),I=1,NNZ),J=1,NNZ)
          PRINT*,'MAX GRAIZING RATE      : ',(G_MAX(I),I=1,NNZ)
        IF(G_FUNCTION.EQ.'RECTI_G'.OR.G_FUNCTION.EQ.'MM1_G' .AND. &
           G_FUNCTION.EQ.'MM2_G') THEN 
          PRINT*,'HALF SATURATION        : ', (K_ZG(I),I=1,NNZ)
        END IF
          PRINT*,'GRAZING POWER          : ',(M_G(I),I=1,NNZ)
          PRINT*,'MORALITY               : ',(MZD(I),I=1,NNZ)
          PRINT*,'MORTALITY POWER        : ',(M_Z(I),I=1,NNZ)
        IF(G_FUNCTION.EQ.'MM2_G') THEN 
          PRINT*,'GRAZING THRESHOLD      : ',(P_C(I),I=1,NNZ)
        END IF
        IF (NNZ.GE.1) THEN
          PRINT*,'RECRUITMENT            : ',(R_RECRUIT(I),I=1,NNZ)
        END IF
          PRINT*,'RESPIRATION            : ',(R_Z(I),I=1,NNZ)
        IF (NNB.GE.1) THEN
          PRINT*,'PREFERENCE ON BACTERIA : ',((SIGMA_B(I,J),I=1,NNB),J=1,NNZ)
        END IF
        IF (NND.GE.1) THEN
          PRINT*,'PREFERENCE ON DETRITUS : ',((SIGMA_D(I,J),I=1,NND),J=1,NNZ)
        END IF
          PRINT*,'PREFERENCE ON PHYTO    : ',((SIGMA_P(I,J),I=1,NNP),J=1,NNZ)
          PRINT*,'PREFERENCE ON ZOO      : ',((SIGMA_Z(I,J),I=1,NNZ),J=1,NNZ)
          PRINT*,'OPTIMAL T              : ',(T_OPTZ(I),I=1,NNZ)
          PRINT*,'ZOO THRESHOLD          : ',(Z_0(I),I=1,NNZ)
          PRINT*
          PRINT*,'*********    NUTRIENT PARAMETERS    ***********'
          PRINT*
          PRINT*,'HALF-SATURATION',((KSN(I,J),I=1,NNN),J=1,NNP)
        IF (NNB.GE.1) THEN
          PRINT*,'ELEMENT RATIO IN BAC.  : ',((N2CB(I,J),I=1,NNN),J=1,NNB)
        END IF
        IF (NND.GE.1) THEN
          PRINT*,'ELEMENT RATIO IN D     : ',((N2CD(I,J),I=1,NNN),J=1,NND)
        END IF
          PRINT*,'ELEMENT RATIO IN PHYTO : ',((N2CP(I,J),I=1,NNN),J=1,NNP)
          PRINT*,'ELEMENT RATIO IN ZOO   : ',((N2CZ(I,J),I=1,NNN),J=1,NNZ)
        IF (NNM.GE.1)  THEN
          PRINT*,'ELEMENT RATIO IN DOM   : ', ((N2CDOM(I,J),I=1,NNN),J=1,NNM)
        END IF
          PRINT*,'THRESHOLD              : ', (N_0(I),I=1,NNN)
        IF (NO3_ON)    THEN
          PRINT*,'NITRIFICATION RATE     : ', R_AN
        END IF
          PRINT*
        IF (NNB.GE.1) THEN
          PRINT*,'*********  BACTERIA PARAMETERS  ************'
          PRINT*
          PRINT*,'T_FORCING              : ',(A_TB(I),I=1,NNB)
          PRINT*,'THRESHOLD              : ', (B_0(I),I=1,NNB)
          PRINT*,'RATIO OF NH4 VS DON    : ',(DELTA_B(I),I=1,NNB)
          PRINT*,'EFFICIENCY OF DETRITUS : ',((EFFIBD(I,J),I=1,NND),J=1,NNB)
          PRINT*,'EFFICIENCY OF DOM      : ',((EFFIDOM(I,J),I=1,NNM),J=1,NNB)
          PRINT*,'EFFICIENCY OF NUTRIENT : ', ((EFFIN(I,J),I=1,NNN),J=1,NNB)
          PRINT*,'RESPIRATION            : ',(R_B(I),I=1,NNB)
          PRINT*,'PREFERENCE ON DETRITUS : ',((SIGMA_BD(I,J),I=1,NND),J=1,NNB)
          PRINT*,'PREFERENCE ON DOM      : ',((SIGMA_DOM(I,J),I=1,NNM),J=1,NNB)
          PRINT*,'PREFERENCE ON NUTRIENT : ',((SIGMA_N(I,J),I=1,NNN),J=1,NNB)
          PRINT*,'OPTIMAL TEMPERATURE    : ',(T_OPTB(I),I=1,NNB)
          PRINT*,'MAXIMUM GROWTH RATE    : ',(UBMAX(I),I=1,NNB)
          PRINT*
        END IF
        IF (NND.GE.1) THEN
          PRINT*, '********  DETRITUS PARAMETERS  *************'
          PRINT*
        IF (NNB.GE.1) THEN
          PRINT*,'GRAZING ON B TO D      : ', (((ALPHA_BD(I,J,K),I=1,NND),J=1,NNB),K=1,NNZ)
        END IF
          PRINT*,'AGGREGATION            : ', (ALPHA_DAG(I),I=1,NND)
          PRINT*,'GRAZING LOSS ON D TO D : ', (((ALPHA_DD(I,J,K),I=1,NND),J=1,NND),K=1,NNZ)
          PRINT*,'DISAGGREGATION         : ', (ALPHA_DAG(I),I=1,NND)
          PRINT*,'GRAZING LOSS ON P TO D : ', (((ALPHA_PD(I,J,K),I=1,NND),J=1,NNP),K=1,NNZ)
          PRINT*,'GRAZING LOSS ON Z TO D : ', (((ALPHA_ZD(I,J,K),I=1,NND),J=1,NNZ),K=1,NNZ)
        IF (NNM.GE.1) THEN
          PRINT*,'DISSOLUTION            : ', (D_DOM(I),I=1,NND)
        END IF
          PRINT*,'REMINERALIZATION       : ', (D_RN(I),I=1,NND)
          PRINT*,'THRESHOLD              : ', (D_0(I),I=1,NND)
          PRINT*,'P MORTALITY TO DETRITUS: ', ((EPSILON_PD(I,J),I=1,NND),J=1,NNP)
          PRINT*,'Z MORTALITY TO DETRITUS: ', ((EPSILON_ZD(I,J),I=1,NND),J=1,NNZ)
          PRINT*,'SINING VELOCITY        : ', (W_D(I),I=1,NND)
          PRINT*
        END IF
        IF (NNM.GE.1) THEN
          PRINT*,'*********     DOM PARAMETERS     *************'
          PRINT*
          PRINT*,'DOM AGEING COEFFICIENT : ', (ALPHA_DOM(I),I=1,NNM)
          PRINT*,'PHYTO EXUDATION        : ', ((ALPHA_PDOM(I,J),I=1,NNM),J=1,NNP)
          PRINT*,'DETRITUS DISSOLUTION   : ', ((ALPHA_DDOM(I,J),I=1,NNM),J=1,NND)
        IF (NNB.GE.1) THEN
          PRINT*,'GRAZING LOSS ON B > DOM: ', (((ALPHA_ZBDOM(I,J,K),I=1,NNM),J=1,NNB),K=1,NNZ)
        END IF
        IF (NND.GE.1) THEN
          PRINT*,'GRAZING LOSS ON D > DOM: ', (((ALPHA_ZDDOM(I,J,K),I=1,NNM),J=1,NND),K=1,NNZ)
        END IF
          PRINT*,'GRAZING LOSS ON P > DOM: ', (((ALPHA_ZPDOM(I,J,K),I=1,NNM),J=1,NNP),K=1,NNZ)
          PRINT*,'GRAZING LOSS ON Z > DOM: ', (((ALPHA_ZZDOM(I,J,K),I=1,NNM),J=1,NNZ),K=1,NNZ)
          PRINT*,'THRESHOLD              : ', (DOM_0(I),I=1,NNM)
          PRINT*
        END IF
       PRINT*
       PRINT*,'*****************************************************'
       PRINT*,'*********    END OF BIOLOGICAL MODEL     ************'
       PRINT*,'*****************************************************'
       PRINT*
 88    CONTINUE
       RETURN
     END SUBROUTINE GET_PARAMETER_NEW



END MODULE MOD_PARAMETER
