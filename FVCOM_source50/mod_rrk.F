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

MODULE RRKVAL
# if defined (RRKF)
   
   USE CONTROL
   IMPLICIT NONE
  
   REAL(DP),ALLOCATABLE :: STTEMP0(:)
   REAL(DP),ALLOCATABLE :: STTEMP1(:)
   REAL(DP),ALLOCATABLE :: STEOF(:)
   REAL(DP),ALLOCATABLE :: SDEOF(:)
   REAL(DP),ALLOCATABLE :: TRANS(:,:)
   REAL(SP) ::  OBSERR_EL         !!EL OBSERVATION ERROR SPECIFIED
   REAL(SP) ::  OBSERR_UV         !!UV OBSERVATION ERROR SPECIFIED
   REAL(SP) ::  OBSERR_T          !!TEMPERATURE OBSERVATION ERROR SPECIFIED
   REAL(SP) ::  OBSERR_S          !!SALINITY OBSERVATION ERROR SPECIFIED   
   REAL(SP),ALLOCATABLE :: RRKEL(:)
   REAL(SP),ALLOCATABLE :: RRKU(:,:)
   REAL(SP),ALLOCATABLE :: RRKV(:,:)
   REAL(SP),ALLOCATABLE :: RRKT(:,:)
   REAL(SP),ALLOCATABLE :: RRKS(:,:)

   REAL(SP),ALLOCATABLE :: RRKEL3(:)
   REAL(SP),ALLOCATABLE :: RRKU3(:,:)
   REAL(SP),ALLOCATABLE :: RRKV3(:,:)
   REAL(SP),ALLOCATABLE :: RRKT3(:,:)
   REAL(SP),ALLOCATABLE :: RRKS3(:,:)

   REAL(SP), ALLOCATABLE :: RRKU_REF(:,:,:)
   REAL(SP), ALLOCATABLE :: RRKV_REF(:,:,:)
   REAL(SP), ALLOCATABLE :: RRKEL_REF(:,:)
   REAL(SP), ALLOCATABLE :: RRKT_REF(:,:,:)
   REAL(SP), ALLOCATABLE :: RRKS_REF(:,:,:)
   
   REAL(DP) ::  ELSD, USD, TSD, SSD
   REAL(DP),ALLOCATABLE :: STSD(:)
   REAL(DP),ALLOCATABLE :: SEOF(:,:)
   REAL(DP),ALLOCATABLE :: STMEAN(:)
   INTEGER  nrange
   INTEGER  ::  NLOC=0              !!Number of observation locations
   LOGICAL  ::  EL_ASSIM          !!OPTION FOR CHOSING ELEVATION AS ASSIMILATION VARIABLES
   LOGICAL  ::  UV_ASSIM          !!OPTION FOR CHOSING CURRENT AS ASSIMILATION VARIABLES
   LOGICAL  ::  T_ASSIM           !!OPTION FOR CHOSING TEMPERATURE AS ASSIMILATION VARIABLES
   LOGICAL  ::  S_ASSIM           !!OPTION FOR CHOSING SALINITY AS ASSIMILATION VARIABLES

   LOGICAL  ::  EL_OBS            !!OPTION FOR ELEVATION OBSERVATION DATA
   LOGICAL  ::  UV_OBS            !!OPTION FOR CURRENT OBSERVATION DATA
   LOGICAL  ::  T_OBS             !!OPTION FOR TEMPERATURE OBSERVATION DATA
   LOGICAL  ::  S_OBS             !!OPTION FOR SALINITY OBSERVATION DATA
# endif
END MODULE RRKVAL

MODULE MOD_RRK
   USE MOD_INPUT, only : NC_START
   USE RRKVAL
   USE CONTROL
   USE MOD_UTILS
   USE MOD_NCTOOLS
   IMPLICIT NONE
   SAVE
   
   INTEGER  ::  RRK_EOFCONTR
   CHARACTER(LEN=80) :: REF_START_DATE
   CHARACTER(LEN=80) :: REF_END_DATE
   CHARACTER(LEN=80) :: RRK_START_DATE
   CHARACTER(LEN=80) :: RRK_END_DATE

   TYPE(TIME) :: REF_START_TIME
   TYPE(TIME) :: REF_END_TIME
   TYPE(TIME) :: RRK_START_TIME
   TYPE(TIME) :: RRK_END_TIME

   TYPE(TIME) :: RRK_CYC

   CHARACTER(LEN=80) :: RRK_ASSIM_INTERVAL
   TYPE(TIME) :: RRK_INTERVAL


   INTEGER      REF_INT           !!GLOBAL NUMBER OF THE READING FILE INTERVALS 
   INTEGER      RRK_NOBSMAX
   INTEGER      RRK_OPTION        !!OPTION 1 FOR BAROTROPIC CASE; OPTION 2 FOR BAROCLINIC CASE
   INTEGER      RRK_NEOF          !!NUMBER OF THE EOF  
   REAL(SP) ::  RRK_PSIZE         !!PERTURBATION SIZE  
   REAL(SP) ::  RRK_PSCALE        !!PSEUDO MODEL ERROR    
   REAL(SP) ::  RRK_RSCALE        !!SCALE FACTOR APPLIED TO ONE STANDARD DEVIATION FOR R

   LOGICAL  :: RRK_ON

# if defined (RRKF)
   INTEGER  ::  ICYC
!   LOGICAL  ::  EL_ASSIM          !!OPTION FOR CHOSING ELEVATION AS ASSIMILATION VARIABLES
!   LOGICAL  ::  UV_ASSIM          !!OPTION FOR CHOSING CURRENT AS ASSIMILATION VARIABLES  
!   LOGICAL  ::  T_ASSIM           !!OPTION FOR CHOSING TEMPERATURE AS ASSIMILATION VARIABLES
!   LOGICAL  ::  S_ASSIM           !!OPTION FOR CHOSING SALINITY AS ASSIMILATION VARIABLES
!
!   LOGICAL  ::  EL_OBS            !!OPTION FOR ELEVATION OBSERVATION DATA
!   LOGICAL  ::  UV_OBS            !!OPTION FOR CURRENT OBSERVATION DATA  
!   LOGICAL  ::  T_OBS             !!OPTION FOR TEMPERATURE OBSERVATION DATA
!   LOGICAL  ::  S_OBS             !!OPTION FOR SALINITY OBSERVATION DATA

   INTEGER      INORRK            !!FILE I/O PIPE NUMBER 


   TYPE(NCFILE), POINTER :: NCF_RRKFDATA

   TYPE(NCFILE), POINTER :: NCF_RRKRE
   TYPE(NCFILE), POINTER :: NCF_RRKFCT
   LOGICAL :: LOCAL_DISK
   CHARACTER(LEN=80), parameter :: Local_file ="/scratch/rrk_restart.nc"
   CHARACTER(LEN=80), parameter :: group_file ="rrktemp/rrk_restart.nc"

   CHARACTER(LEN=80), parameter :: Local_file1 ="/scratch/rrk_forecast.nc"
   CHARACTER(LEN=80), parameter :: group_file1 ="rrktemp/rrk_forecast.nc"

   INTEGER(ITIME) :: RKINT_COUNT, RKINT_START, RKINT_END


   INTERFACE READRESTART
      MODULE PROCEDURE READRESTART_TIME
      MODULE PROCEDURE READRESTART_STATE
   END INTERFACE

   PRIVATE READRESTART_TIME
   PRIVATE READRESTART_STATE


   NAMELIST /NML_RRKF/        &
         & RRK_ON,            &
         & REF_START_DATE,    &
         & REF_END_DATE,      &
         & RRK_START_DATE,    &
         & RRK_END_DATE,      &
         & RRK_ASSIM_INTERVAL,&
         & RRK_NOBSMAX,       &
         & RRK_NEOF,          &
         & RRK_PSIZE,         &
         & RRK_PSCALE,        &
         & RRK_RSCALE,        &
         & EL_ASSIM,          &
         & EL_OBS,            &
         & UV_ASSIM,          &
         & UV_OBS,            &
         & T_ASSIM,           &
         & T_OBS,             &
         & S_ASSIM,           &
         & S_OBS,             &
         & LOCAL_DISK

   CONTAINS
   
     SUBROUTINE NAME_LIST_INITIALIZE_RRKF
       IMPLICIT NONE
       REF_START_DATE = "RRK REFERENCE START AND END TIME FOR EOF CALCULATION"
       REF_END_DATE   = "Date and Time are specified as a string (example '2007-11-05 00:00:00')"  
       RRK_START_DATE = "RRK ASSIMILATION START AND END TIME"
       RRK_END_DATE   = "For an idealized case specify 'seconds=(flt)','days=(flt)', or 'cycles=(int)'" 
       RRK_ASSIM_INTERVAL= "A length of time: 'seconds= ','days= ', or 'cycles= '"   
       RRK_NOBSMAX    = -1  
       RRK_NEOF       = -1  
       RRK_PSIZE      = -1  
       RRK_PSCALE     = -1.0_SP  
       RRK_RSCALE     = -1.0_SP  
       RRK_ON         = .FALSE.
       EL_ASSIM       = .FALSE.  
       EL_OBS         = .FALSE.  
       UV_ASSIM       = .FALSE.  
       UV_OBS         = .FALSE. 
       T_ASSIM        = .FALSE.   
       T_OBS          = .FALSE.  
       S_ASSIM        = .FALSE.  
       S_OBS          = .FALSE.
       LOCAL_DISK     = .FALSE.

     END SUBROUTINE NAME_LIST_INITIALIZE_RRKF
     
     SUBROUTINE NAME_LIST_PRINT_RRKF
       IMPLICIT NONE

       write(UNIT=IPT,NML=NML_RRKF)

     RETURN
     END SUBROUTINE NAME_LIST_PRINT_RRKF
     
     SUBROUTINE NAME_LIST_READ_RRKF
       IMPLICIT NONE
       CHARACTER(LEN=120) :: FNAME
       INTEGER :: IOS

       FNAME = TRIM(INPUT_DIR)//"/"//trim(casename)//"_assim_rrkf.nml"
       
       CALL FOPEN(KFUNIT,trim(FNAME),'cfr')

       ! Read RRKF parameters
       READ(UNIT=KFUNIT, NML=NML_RRKF,IOSTAT=ios)
       if(ios .NE. 0 ) then
          if(DBG_SET(dbg_log)) write(UNIT=IPT,NML=NML_RRKF)
          Call Fatal_Error("Can Not Read NameList NML_RRKF from file: "//trim(FNAME))
       end if
       
       REWIND(KFUNIT)
       
       if(DBG_SET(dbg_scl)) &
            & write(IPT,*) "Read_Name_List: NML_RRKF"
       
       if(DBG_SET(dbg_scl)) &
            & write(UNIT=IPT,NML=NML_RRKF)
       
       CLOSE(NMLUNIT)
       
       RETURN
     END SUBROUTINE NAME_LIST_READ_RRKF
     
     SUBROUTINE RRK_SET_TIME
       USE MOD_SET_TIME
       IMPLICIT NONE
       CHARACTER(LEN=4) :: FLAG,BFLAG,EFLAG
       INTEGER(ITIME):: begins, ends,steps
       integer :: status
# if defined (DATA_ASSIM)
   FVCOM_RUN_MODE = FVCOM_RRKF_WITH_SSA
# else
   FVCOM_RUN_MODE = FVCOM_RRKF_WITHOUT_SSA
# endif

       
       ! GET THE RRK ASSIMILATION INTERVAL PERIOD
       CALL IDEAL_TIME_STRING2TIME(RRK_ASSIM_INTERVAL,FLAG,RRK_INTERVAL,steps)
       IF (FLAG == 'time') THEN ! IF START AND END TIME WERE SPECIFIED

          IF(MOD(RRK_INTERVAL,IMDTI) /= ZeroTime) THEN
             CALL FATAL_ERROR("RRK_ASSIM_INTERVAL must be an integer number of internal time steps!")
          END IF
          
          steps = seconds(RRK_INTERVAL)/seconds(IMDTI)


       ELSEIF(FLAG == 'step') THEN ! IF SPECIFIED AS A NUMBER OF STEPS
          
          RRK_INTERVAL = steps * IMDTI
          
       END IF

       RKInt_start=0
       RKint_end = steps


       SELECT CASE(USE_REAL_WORLD_TIME)
       CASE(.TRUE.)

          REF_START_TIME = READ_DATETIME(REF_START_DATE,DATE_FORMAT,TIMEZONE,status)
!!$          if (.not. status) &
          if (status == 0) &
               & Call Fatal_Error("Could not read the date string REF_START_DATE: "//trim(REF_START_DATE))


          ! GET THE END TIME
          REF_END_Time = READ_DATETIME(REF_END_DATE,DATE_FORMAT,TIMEZONE,status)
!!$          if (.not. status) &
          if (status == 0) &
               & Call Fatal_Error("Could not read the date string REF_END_DATE: "&
               &//trim(REF_END_DATE))

          ! SANITY CHECK
          if(REF_Start_Time .GE. REF_End_Time) &
               & Call Fatal_Error("RRKF REF_Start_Date exceeds or equal to REF_End_Date")
          
          RRK_START_TIME = READ_DATETIME(RRK_START_DATE,DATE_FORMAT,TIMEZONE,status)
!!$          if (.not. status) &
          if (status == 0) &
               & Call Fatal_Error("Could not read the date string RRK_START_DATE: "//trim(RRK_START_DATE))
          
          
          ! GET THE END TIME
          RRK_END_Time = READ_DATETIME(RRK_END_DATE,DATE_FORMAT,TIMEZONE,status)
!!$          if (.not. status) &
          if (status == 0) &
               & Call Fatal_Error("Could not read the date string RRK_END_DATE: "&
               &//trim(RRK_END_DATE))

          ! SANITY CHECK
          if(RRK_Start_Time .GE. RRK_End_Time) &
               & Call Fatal_Error("RRKF RRK_Start_Date exceeds or equal to RRK_End_Date")
          


       CASE (.FALSE.)

          ! GET THE REFERENCE START AND END INFORMATION
          CALL IDEAL_TIME_STRING2TIME(REF_START_DATE,BFLAG,REF_START_TIME,begins)

          CALL IDEAL_TIME_STRING2TIME(REF_END_DATE,EFLAG,REF_End_TIME,ends)

          ! SANITY CHECK
          IF (BFLAG /= EFLAG) CALL FATAL_ERROR&
               ('IDEALIZED MODEL TIME SPECIFICATION IS INCORRENT',&
               &'RRK REF BEGIN AND END CAN BE IN EITHER CYCLES OR TIME BUT NOT MIXED',&
               & trim(REF_start_date),trim(ref_end_date) )
          
          IF (EFLAG == 'time') THEN ! IF START AND END TIME WERE SPECIFIED

             
             ! DO NOTHING

          ELSE IF(EFLAG == 'step') THEN ! IF START AND END IINT WERE SPECIFIED


             REF_START_TIME = StartTime + IMDTI*(begins - ISTART)

             REF_END_TIME = StartTime + IMDTI*(ends - ISTART)

             
          ELSE
             CALL FATAL_ERROR&
                  &('IDEAL_TIME_STRING2TIME returned invalid flag for rrk ref date and time')

          END IF

          ! GET THE RRK START AND END INFORMATION
          CALL IDEAL_TIME_STRING2TIME(RRK_START_DATE,BFLAG,RRK_START_TIME,begins)

          CALL IDEAL_TIME_STRING2TIME(RRK_END_DATE,EFLAG,RRK_End_TIME,ends)

          ! SANITY CHECK
          IF (BFLAG /= EFLAG) CALL FATAL_ERROR&
               ('IDEALIZED MODEL TIME SPECIFICATION IS INCORRENT',&
               &'RRK BEGIN AND END CAN BE IN EITHER CYCLES OR TIME BUT NOT MIXED',&
               & trim(ref_start_date),trim(ref_end_date) )
          
          IF (EFLAG == 'time') THEN ! IF START AND END TIME WERE SPECIFIED

             
             ! DO NOTHING

          ELSE IF(EFLAG == 'step') THEN ! IF START AND END IINT WERE SPECIFIED


             RRK_START_TIME = StartTime + IMDTI*(begins - ISTART)

             RRK_END_TIME = StartTime + IMDTI*(ends - ISTART)

             
          ELSE
             CALL FATAL_ERROR&
                  &('IDEAL_TIME_STRING2TIME returned invalid flag for rrk ref date and time')

          END IF



       END SELECT


       IF(MOD(REF_START_TIME-StartTime,IMDTI) /= ZeroTime) THEN
          CALL FATAL_ERROR("REF_START_TIME must be an integer number of internal time steps from the start time!")
       END IF
       
       IF(MOD(REF_END_TIME-StartTime,IMDTI) /= ZeroTime) THEN
          CALL FATAL_ERROR("REF_END_TIME must be an integer number of internal time steps from the start time!")
       END IF
       
       
       IF(MOD(RRK_START_TIME-StartTime,IMDTI) /= ZeroTime) THEN
          CALL FATAL_ERROR("RRK_START_TIME must be an integer number of internal time steps from the start time!")
       END IF
       
       IF(MOD(RRK_END_TIME-StartTime,IMDTI) /= ZeroTime) THEN
          CALL FATAL_ERROR("RRK_END_TIME must be an integer number of internal time steps from the start time!")
       END IF
       
       
     END SUBROUTINE RRK_SET_TIME


     SUBROUTINE RRK_SET_STARTUP
       USE MOD_SET_TIME
       USE MOD_INPUT
       IMPLICIT NONE
       character(len=160) :: pathnfile
       TYPE(NCFILE), POINTER :: NCF

       OPEN(75,FILE='./err.dat')

       !RESTART FILE
       RESTART_FILE_NAME = trim(casename)//"_restart_0001.nc"
       pathnfile= trim(OUTPUT_DIR)//TRIM(RESTART_FILE_NAME)
       
       CALL SEARCH_FOR_LAST_MATCHING_NAME(PATHNFILE)
       
       ! OPEN THE FILE AND LOAD FOR STARTUP
       NCF => NEW_FILE()
       NCF%FNAME=trim(pathnfile)
       Call NC_OPEN(NCF)
       CALL NC_LOAD(NCF)
       NC_START => NCF
       
       Nullify(NCF)
       CALL SET_STARTUP_FILE_STACK(REF_START_TIME)

!       STARTUP_TYPE = STARTUP_TYPE_CRASHRESTART

!       STARTUP_UV_TYPE = STARTUP_TYPE_SETVALUES

!       STARTUP_TURB_TYPE = STARTUP_TYPE_SETVALUES

!       STARTUP_TS_TYPE = STARTUP_TYPE_SETVALUES

     END SUBROUTINE RRK_SET_STARTUP


!!=====================================================================================
!!
!! READ IN THE SIMULATION DATA AND CALCULATE THE EOF
!!
!!=====================================================================================

   
   SUBROUTINE RRK_REF 
   USE CONTROL
   USE RRKVAL
   USE MOD_NCTOOLS
   IMPLICIT NONE

   TYPE(NCFILE), POINTER :: NCF
   TYPE(TIME) :: tTEST
   LOGICAL :: FOUND
   INTEGER :: NTIMES
   
   INTEGER IDUMMY, I, K, J
   INTEGER SS_DIM
   INTEGER STDIM
!   REAL(DP) ::  ELSD, USD, TSD, SSD
   INTEGER LWORK4, LDVT
   REAL(DP),ALLOCATABLE :: WORK4(:)
   REAL(DP) :: VT
   REAL(DP),ALLOCATABLE :: RKSF(:,:)
   REAL(DP),ALLOCATABLE :: RKSF1(:,:)
!   REAL(DP),ALLOCATABLE :: SEOF(:,:)
   REAL(DP),ALLOCATABLE :: SFSF(:,:)
   REAL(DP),ALLOCATABLE :: SFD(:)
   REAL(DP),ALLOCATABLE :: SFU(:,:)
   REAL(DP),ALLOCATABLE :: STVAR(:)
!   REAL(DP),ALLOCATABLE :: STSD(:)
   REAL(DP) :: SUM0   
   INTEGER RCODE
   INTEGER IDUMMY1
   REAL(DP),ALLOCATABLE :: RRKU2(:,:,:)
   REAL(DP),ALLOCATABLE :: RRKV2(:,:,:)
   REAL(DP),ALLOCATABLE :: RRKEL2(:,:)
   REAL(DP),ALLOCATABLE :: RRKT2(:,:,:)
   REAL(DP),ALLOCATABLE :: RRKS2(:,:,:)
!   REAL(DP),ALLOCATABLE :: STMEAN(:)

   TYPE(NCDIM), POINTER :: DIM
   TYPE(NCVAR), POINTER :: VAR

   INTEGER :: IRANGE(2),CNT, STATUS
   character(len=160) :: pathnfile 
   
   pathnfile = trim(INPUT_DIR)//trim(casename)//"_true_0001.nc"
   
   NCF => NEW_FILE()
   NCF%FNAME=trim(pathnfile)
   
   Call NC_OPEN(NCF)
   CALL NC_LOAD(NCF)

   ! GET THE TIME DIMENSION
   DIM => FIND_UNLIMITED(NCF,FOUND)
   IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN RRK_REF FILE",&
            & "FILE NAME: "//TRIM(NCF%FNAME),&
            &"COULD NOT FIND THE UNLIMITED DIMENSION")
   NTIMES = DIM%DIM
   

   ! TEST THE FILE START AND END
   TTEST = get_file_time(NCF,1)
   
   IF(REF_START_TIME < TTEST) CALL FATAL_ERROR &
        &("RRK_REF: REF_START_TIME IS BEFORE FIRST TIME IN FILE.",&
        & "FILE NAME:"//TRIM(NCF%FNAME))
   
   
   TTEST = get_file_time(NCF,NTIMES)
   IF(REF_END_TIME > TTEST) CALL FATAL_ERROR &
        &("RRK_REF: REF_START_TIME IS BEFORE FIRST TIME IN FILE.",&
        & "FILE NAME:"//TRIM(NCF%FNAME))


   ! GET THE START AND END INDEX

   IRANGE = 0

   CNT = 0
   DO WHILE (product(irange)==0)
      CNT = CNT +1

      IF (CNT > NTIMES) THEN
         CALL PRINT_TIME(REF_START_TIME,IPT,"RRK_REF START TIME")
         CALL PRINT_TIME(REF_END_TIME,IPT,"RRK_REF END TIME")
         

         CALL FATAL_ERROR&
           &("RRK_REF: CAN NOT FIND SPECIFIED START AND END TIME IN THE REFERENCE FILE!", &
           & "FILE NAME:"//TRIM(NCF%FNAME))
      END IF


      TTEST = get_file_time(NCF,CNT)
      IF (REF_START_TIME == TTEST) THEN
         IRANGE(1) = CNT
      END IF

      IF (REF_END_TIME == TTEST) THEN
         IRANGE(2) = CNT
      END IF

   END DO


   nrange = irange(2)-irange(1)+1

   IF(EL_ASSIM) THEN
      
      VAR => FIND_VAR(NCF,"zeta",FOUND)
      IF(.not. FOUND) CALL FATAL_ERROR &
           & ("IN RRK_REF FILE",&
            & "FILE NAME: "//TRIM(NCF%FNAME),&
            &"COULD NOT FIND THE VARIABLE 'zeta'")

      allocate(RRKEL_REF(MGL,nrange),stat=status)
      IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY FOR RRKEL_REF")
      RRKEL_REF = 0.0_SP

      CALL nc_connect_avar(VAR,RRKEL_REF)

      CALL NC_READ_VAR(var,stkrng=irange,parallel=.false.)


   ENDIF
   
   IF(UV_ASSIM) THEN

      VAR => FIND_VAR(NCF,"u",FOUND)
      IF(.not. FOUND) CALL FATAL_ERROR &
           & ("IN RRK_REF FILE",&
            & "FILE NAME: "//TRIM(NCF%FNAME),&
            &"COULD NOT FIND THE VARIABLE 'u'")

      allocate(RRKU_REF(NGL,kbm1,nrange),stat=status)
      IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY FOR RRKU_REF")
      RRKU_REF = 0.0_SP
      CALL nc_connect_avar(VAR,RRKU_REF)

      CALL NC_READ_VAR(var,stkrng=irange,parallel=.false.)

      VAR => FIND_VAR(NCF,"v",FOUND)
      IF(.not. FOUND) CALL FATAL_ERROR &
           & ("IN RRK_REF FILE",&
            & "FILE NAME: "//TRIM(NCF%FNAME),&
            &"COULD NOT FIND THE VARIABLE 'v'")

      allocate(RRKV_REF(NGL,kbm1,nrange),stat=status)
      IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY FOR RRKV_REF")
      RRKV_REF = 0.0_SP
      CALL nc_connect_avar(VAR,RRKV_REF)

      CALL NC_READ_VAR(var,stkrng=irange,parallel=.false.)

      
   ENDIF
   
   IF(T_ASSIM) THEN
      VAR => FIND_VAR(NCF,"temp",FOUND)
      IF(.not. FOUND) CALL FATAL_ERROR &
           & ("IN RRK_REF FILE",&
            & "FILE NAME: "//TRIM(NCF%FNAME),&
            &"COULD NOT FIND THE VARIABLE 'temp'")

      allocate(RRKT_REF(MGL,kbm1,nrange),stat=status)
      IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY FOR RRKT_REF")
      RRKT_REF = 0.0_SP
      CALL nc_connect_avar(VAR,RRKT_REF)

      CALL NC_READ_VAR(var,stkrng=irange,parallel=.false.)
      
   ENDIF
   
   IF(S_ASSIM) THEN
      VAR => FIND_VAR(NCF,"salinity",FOUND)
      IF(.not. FOUND) CALL FATAL_ERROR &
           & ("IN RRK_REF FILE",&
            & "FILE NAME: "//TRIM(NCF%FNAME),&
            &"COULD NOT FIND THE VARIABLE 'salinity'")

      allocate(RRKS_REF(MGL,kbm1,nrange),stat=status)
      IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY FOR RRKS_REF")
      RRKS_REF = 0.0_SP
      CALL nc_connect_avar(VAR,RRKS_REF)

      CALL NC_READ_VAR(var,stkrng=irange,parallel=.false.)
      
   ENDIF
   
   CALL KILL_FILE(NCF)
   
   WRITE(IPT,*) 'Calculate the EOFs from the control run......'
   
   STDIM = 0
   IF(EL_ASSIM) STDIM = STDIM + MGL
   IF(UV_ASSIM) STDIM = STDIM + 2*NGL*KBM1
   IF(T_ASSIM)  STDIM = STDIM + MGL*KBM1
   IF(S_ASSIM)  STDIM = STDIM + MGL*KBM1
    
   SS_DIM = nrange
   LWORK4 = 5*SS_DIM
   LDVT   = 1   
   
   ALLOCATE(WORK4(LWORK4))             ;WORK4   = ZERO
   ALLOCATE(RKSF(STDIM,SS_DIM))        ;RKSF    = ZERO
   ALLOCATE(RKSF1(STDIM,SS_DIM))       ;RKSF1   = ZERO
   ALLOCATE(SEOF(STDIM,SS_DIM))        ;SEOF    = ZERO
   ALLOCATE(SFSF(SS_DIM,SS_DIM))       ;SFSF    = ZERO
   ALLOCATE(SFD(SS_DIM))               ;SFD     = ZERO
   ALLOCATE(SFU(SS_DIM,SS_DIM))        ;SFU     = ZERO
   ALLOCATE(STVAR(STDIM))              ;STVAR   = ZERO
   ALLOCATE(STSD(STDIM))               ;STSD    = ZERO
   ALLOCATE(RRKU(NGL,KBM1))            ;RRKU    = ZERO
   ALLOCATE(RRKV(NGL,KBM1))            ;RRKV    = ZERO
   ALLOCATE(RRKEL(MGL))                ;RRKEL   = ZERO
   ALLOCATE(RRKT(MGL,KBM1))            ;RRKT    = ZERO
   ALLOCATE(RRKS(MGL,KBM1))            ;RRKS    = ZERO
   ALLOCATE(RRKU2(NGL,KBM1,SS_DIM))    ;RRKU2   = ZERO
   ALLOCATE(RRKV2(NGL,KBM1,SS_DIM))    ;RRKV2   = ZERO
   ALLOCATE(RRKEL2(MGL,SS_DIM))        ;RRKEL2  = ZERO
   ALLOCATE(RRKT2(MGL,KBM1,SS_DIM))    ;RRKT2   = ZERO
   ALLOCATE(RRKS2(MGL,KBM1,SS_DIM))    ;RRKS2   = ZERO
   ALLOCATE(STMEAN(STDIM))             ;STMEAN  = ZERO  
   
   
   ! STORE IN RKSF

    
    DO I = 1, SS_DIM
  
      IDUMMY = 0      

      IF(EL_ASSIM) THEN      
        DO J = 1, MGL
          IDUMMY = IDUMMY + 1
	  RKSF(IDUMMY,I) =  RRKEL_REF(J,I)
        ENDDO 
      ENDIF

      IF(UV_ASSIM) THEN
        DO K = 1, KBM1
          DO J = 1, NGL
	     IDUMMY = IDUMMY + 1
	     RKSF(IDUMMY,I) = RRKU_REF(J,K,I)
	  ENDDO
        ENDDO	

        DO K = 1, KBM1
          DO J = 1, NGL
	     IDUMMY = IDUMMY + 1
	     RKSF(IDUMMY,I) = RRKV_REF(J,K,I)
	  ENDDO
        ENDDO     
      ENDIF

      IF(T_ASSIM) THEN
        DO K = 1, KBM1
          DO J = 1, MGL
	     IDUMMY = IDUMMY + 1
	     RKSF(IDUMMY,I) = RRKT_REF(J,K,I)
	  ENDDO
        ENDDO	
      ENDIF

      IF(S_ASSIM) THEN
	DO K = 1, KBM1
          DO J = 1, MGL
	     IDUMMY = IDUMMY + 1
	     RKSF(IDUMMY,I) = RRKS_REF(J,K,I)
	  ENDDO
        ENDDO	
      ENDIF
      
    ENDDO
    
      
! CALCULATE THE MEAN AND THE STANDARD DEVIATION

    DO I=1, STDIM
      SUM0=0.0_DP
      DO J=1, SS_DIM
        SUM0=SUM0+RKSF(I,J)
      ENDDO
     
      STMEAN(I) = SUM0/DBLE(SS_DIM)
      SUM0=0.0_DP
      DO J=1, SS_DIM
        RKSF1(I,J)=(RKSF(I,J)-STMEAN(I))
	SUM0=SUM0+RKSF1(I,J)**2.0_DP
      ENDDO
           
      STVAR(I)=SUM0
      STSD(I)=DSQRT(SUM0/DBLE(SS_DIM-1))
    ENDDO    
    

    IDUMMY  = 0
    IDUMMY1 = 0
    ELSD    = 0.0_DP
    USD     = 0.0_DP
    TSD     = 0.0_DP
    SSD     = 0.0_DP          

    IF(EL_ASSIM) THEN
      SUM0=0.0_DP
      DO I=1, MGL
        IDUMMY = IDUMMY + 1
        SUM0 = SUM0 + STVAR(IDUMMY)
      ENDDO
      
      ELSD=DSQRT(SUM0/DBLE(MGL)/DBLE(SS_DIM-1))

      DO I=1, MGL
        DO J=1, SS_DIM
          RKSF1(I,J)=RKSF1(I,J)/ELSD/DSQRT(DBLE(SS_DIM-1)) 
        ENDDO
      ENDDO
    
      IDUMMY1 = IDUMMY
    ENDIF

    IF(UV_ASSIM) THEN
      SUM0=0.0_DP
      DO K=1, 2*KBM1
        DO I=1, NGL
           IDUMMY = IDUMMY + 1
	   SUM0 = SUM0 +STVAR(IDUMMY)
        ENDDO
      ENDDO
      USD=DSQRT(SUM0/DBLE(NGL*KBM1)/DBLE(SS_DIM-1))
 
      DO I=IDUMMY1+1, IDUMMY1+2*NGL*KBM1 
        DO J=1, SS_DIM	    
	   RKSF1(I,J)=RKSF1(I,J)/USD/DSQRT(DBLE(KBM1))/DSQRT(DBLE(SS_DIM-1))         
        ENDDO 
      ENDDO
    
      IDUMMY1 = IDUMMY
    ENDIF  
      
    IF(T_ASSIM) THEN
      SUM0=0.0_DP
      DO K=1, KBM1
        DO I=1, MGL
           IDUMMY = IDUMMY + 1
	   SUM0 = SUM0 +STVAR(IDUMMY)
        ENDDO
      ENDDO
      TSD=DSQRT(SUM0/DBLE(MGL*KBM1)/DBLE(SS_DIM-1))
 
      DO I=IDUMMY1+1, IDUMMY1+MGL*KBM1 
        DO J=1, SS_DIM
	  RKSF1(I,J)=RKSF1(I,J)/TSD/DSQRT(DBLE(SS_DIM-1))
        ENDDO 
      ENDDO

      IDUMMY1 = IDUMMY
    ENDIF        
      
    IF(S_ASSIM) THEN
      SUM0=0.0_DP
      DO K=1, KBM1
        DO I=1, MGL
           IDUMMY = IDUMMY + 1
	   SUM0 = SUM0 +STVAR(IDUMMY)
        ENDDO
      ENDDO
      SSD=DSQRT(SUM0/DBLE(MGL*KBM1)/DBLE(SS_DIM-1))
 
      DO I=IDUMMY1+1, IDUMMY1+MGL*KBM1 
        DO J=1, SS_DIM
	   RKSF1(I,J)=RKSF1(I,J)/SSD/DSQRT(DBLE(SS_DIM-1))
        ENDDO 
      ENDDO

      IDUMMY1 = IDUMMY
    ENDIF       


! CALCULATE THE EOFs BY SVD OF RKSF1' RKSF1 = C LAMBDA C', RKSF1 RKSF1' = E LAMBDA E', E = RKSF1 C LAMBDA ^-1/2

    DO I=1, SS_DIM
      DO J=1, SS_DIM
        SUM0=0.0_DP
	DO K=1, STDIM
	  SUM0 = SUM0 + RKSF1(K,I)*RKSF1(K,J)
	ENDDO
	SFSF(I,J) =SUM0
      ENDDO
    ENDDO
 
    CALL DGESVD('A','N',SS_DIM,SS_DIM,SFSF,SS_DIM,SFD,SFU,SS_DIM,VT,LDVT,WORK4,LWORK4,RCODE)    
    
    OPEN(72,FILE=TRIM(OUTPUT_DIR)//'/rrktemp/'//'eigenvalue.dat')
    
    DO I=1, SS_DIM
      WRITE(72,'(I5,E15.7)') I, SFD(I)   
    ENDDO
    DO I=1, SS_DIM
      DO J=1, SS_DIM 
        SFU(I,J)=SFU(I,J)/DSQRT(SFD(J))
      ENDDO
    ENDDO

    DO I=1, STDIM
      DO J=1, SS_DIM
         SUM0=0.0_DP 
         DO K=1, SS_DIM
	    SUM0=SUM0+RKSF1(I,K)*SFU(K,J)
	 ENDDO
         SEOF(I,J)=SUM0 
      ENDDO
    ENDDO
    
    write(1000,*) (seof(i,1),i=1,stdim)
    
    DEALLOCATE(WORK4,RKSF,RKSF1,SFSF,SFD,SFU,STVAR)
    DEALLOCATE(RRKU2,RRKV2,RRKEL2,RRKT2,RRKS2)
    DEALLOCATE(RRKU,RRKV,RRKEL,RRKT,RRKS)

 END SUBROUTINE RRK_REF
 
!------------------------------------------------------------------------------|
!  Main program to calculate the reduced rank kalman gain matrix               |
!------------------------------------------------------------------------------|   
   
   SUBROUTINE RRK_RRK(CHOICE)

   USE LIMS
   USE CONTROL
   USE RRKVAL
   IMPLICIT NONE
    
    INTEGER SS_DIM
    INTEGER STDIM
    INTEGER I_EOF
    INTEGER I,J,K
    INTEGER CHOICE
!    INTEGER NLOC
    INTEGER RCODE
    INTEGER IT
    INTEGER ILAST
    INTEGER IDUMMY
    INTEGER SS_TOT
    INTEGER,ALLOCATABLE  :: STLOC(:)
    REAL(DP) ERR1
    REAL(DP),ALLOCATABLE :: KAL(:,:)
    REAL(DP),ALLOCATABLE :: HEOF(:,:)
    REAL(DP),ALLOCATABLE :: R(:,:)
    REAL(DP),ALLOCATABLE :: HTR(:,:)
    REAL(DP),ALLOCATABLE :: R2(:,:)
    REAL(DP),ALLOCATABLE :: PHT(:,:)
    REAL(DP),ALLOCATABLE :: RALPHA(:,:)
    REAL(DP),ALLOCATABLE :: BETA2(:,:)
    REAL(DP),ALLOCATABLE :: GAMMA(:,:)
    REAL(DP),ALLOCATABLE :: GAMMA2(:,:)
    REAL(DP),ALLOCATABLE :: TRANS2(:,:)
    
! WORK ARRAYS FOR LAPACK SUBROUTINES
    INTEGER LWORK, LWORK2, LDVT
    INTEGER,ALLOCATABLE :: IPIV(:)
    INTEGER,ALLOCATABLE :: IPIV2(:)
    REAL(DP),ALLOCATABLE :: WORK(:)
    REAL(DP),ALLOCATABLE :: WORK2(:)
    REAL(DP) :: VT
    

! CHARACTER STRINGS
    CHARACTER(LEN=80) KALNAME
    CHARACTER(LEN=80) FILENAME
    CHARACTER(LEN=4)  STRCYC
    CHARACTER(LEN=8)  CH8
    CHARACTER(LEN=80) INAME    
    CHARACTER(LEN=80) INAME2
    
!    REAL(DP),ALLOCATABLE :: SEOF(:,:)
!    REAL(DP),ALLOCATABLE :: STSD(:)
!    REAL(DP) :: ELSD, USD, TSD, SSD
    REAL(DP),ALLOCATABLE :: HU(:,:)
    REAL(DP),ALLOCATABLE :: HUL(:,:)
    REAL(DP),ALLOCATABLE :: EVAL(:)
    REAL(DP) :: RRKSUM
      
    STDIM = 0
    IF(EL_ASSIM) STDIM = STDIM + MGL
    IF(UV_ASSIM) STDIM = STDIM + 2*NGL*KBM1
    IF(T_ASSIM)  STDIM = STDIM + MGL*KBM1
    IF(S_ASSIM)  STDIM = STDIM + MGL*KBM1

    SS_DIM = RRK_NEOF
    LWORK  = 4*SS_DIM
    LWORK2 = 4*RRK_NOBSMAX
    SS_TOT = nrange
    LDVT   = 1 
    ILAST  = 10  ! THE NUMBER OF ITERATION IN DOUBLING ALGORITHM

! TEMPORARILY ALLOCATE ARRY TO ARRYS


   ALLOCATE(STLOC(RRK_NOBSMAX))            ;STLOC    = 0
   ALLOCATE(IPIV(SS_DIM))                  ;IPIV     = 0
   ALLOCATE(IPIV2(RRK_NOBSMAX))            ;IPIV2    = 0
   ALLOCATE(KAL(SS_DIM,RRK_NOBSMAX))       ;KAL      = ZERO
   ALLOCATE(HEOF(RRK_NOBSMAX,SS_DIM))      ;HEOF     = ZERO
   ALLOCATE(R(RRK_NOBSMAX,RRK_NOBSMAX))    ;R        = ZERO
   ALLOCATE(HTR(SS_DIM,RRK_NOBSMAX))       ;HTR      = ZERO
   ALLOCATE(R2(RRK_NOBSMAX,RRK_NOBSMAX))   ;R2       = ZERO
   ALLOCATE(PHT(SS_DIM,RRK_NOBSMAX))       ;PHT      = ZERO
   ALLOCATE(RALPHA(SS_DIM,SS_DIM))         ;RALPHA   = ZERO
   ALLOCATE(BETA2(SS_DIM,SS_DIM))          ;BETA2    = ZERO
   ALLOCATE(GAMMA(SS_DIM,SS_DIM))          ;GAMMA    = ZERO
   ALLOCATE(GAMMA2(SS_DIM,SS_DIM))         ;GAMMA2   = ZERO
   ALLOCATE(TRANS2(SS_DIM,SS_DIM))         ;TRANS2   = ZERO 
   ALLOCATE(WORK(LWORK))                   ;WORK     = ZERO
   ALLOCATE(WORK2(LWORK2))                 ;WORK2    = ZERO
!   ALLOCATE(SEOF(STDIM,SS_DIM))            ;SEOF     = ZERO
!   ALLOCATE(STSD(STDIM))                   ;STSD     = ZERO
   ALLOCATE(HU(RRK_NOBSMAX,SS_TOT))        ;HU       = ZERO      
   ALLOCATE(HUL(RRK_NOBSMAX,SS_TOT))       ;HUL      = ZERO
   ALLOCATE(EVAL(SS_TOT))                  ;EVAL     = ZERO


! END ALLOCATION

   IF(CHOICE == 1) THEN
     IF(MSR) WRITE(IPT,*) 'Starting perturb the base state in the eigenvector direction......'
     print *, 'call perturb' 
     CALL PERTURB

     DEALLOCATE(STLOC,IPIV,IPIV2,KAL,HEOF,R,HTR,R2,PHT,RALPHA,TRANS2)
     DEALLOCATE(BETA2,GAMMA,GAMMA2,WORK,WORK2,HU,HUL,EVAL) 

     RETURN
   ENDIF

   IF(CHOICE == 2) THEN
     IF(MSR) WRITE(IPT,*) 'Starting calculate the linearized model matrix in the reduced space......'  
! CALCULATE THE LINEARIZED MODEL MATRIX IN THE REDUCED SPACE
     CALL MREDUCED
     
     IF(MSR) THEN
! OUTPUT LINEAR TRANSITION MATRIX   
       FILENAME = TRIM(OUTPUT_DIR)//'/rrktemp/Amatr.dat' 
       OPEN(INORRK,FILE=TRIM(FILENAME),FORM='UNFORMATTED')
       DO J=1, SS_DIM
         WRITE(INORRK) (TRANS(I,J), I=1, SS_DIM)    
       ENDDO 
       CLOSE(INORRK)

 ! ALSO IN ASCII
       FILENAME = TRIM(OUTPUT_DIR)//'/rrktemp/Amatr1.txt'
       OPEN(INORRK, FILE=TRIM(FILENAME),FORM='FORMATTED')
       DO I=1, SS_DIM
         DO J=1, SS_DIM
           WRITE(INORRK,*) TRANS(J,I)
         ENDDO
       ENDDO
       CLOSE(INORRK)
     ENDIF  

     DEALLOCATE(STLOC,IPIV,IPIV2,KAL,HEOF,R,HTR,R2,PHT,RALPHA,TRANS2)
     DEALLOCATE(BETA2,GAMMA,GAMMA2,WORK,WORK2,HU,HUL,EVAL) 
!     DEALLOCATE(BETA2,GAMMA,GAMMA2,WORK,WORK2,SEOF,STSD,HU,HUL,EVAL)     
!     DEALLOCATE(STTEMP0,STTEMP1,STEOF,SDEOF,TRANS,RRKEL,RRKU,RRKV,RRKT,RRKS)
     RETURN
   ENDIF
       
   IF(CHOICE == 4) THEN
     CALL RRK_ALLOC_VAR
     WRITE(IPT,*) 'Calculate the Kalman gain matrix by doubling algorith......'
     KALNAME = TRIM(OUTPUT_DIR)//'/rrktemp/rrK.dat'

! READ THE EIGENVALUES AND SET THE INITIAL ERROR COVARIANCE   
     INAME = TRIM(OUTPUT_DIR)//'/rrktemp/eigenvalue.dat' 
     OPEN(INORRK,FILE=INAME,STATUS='OLD')
     
      WRITE(IPT,*) 'ss_tot:', ss_tot, ss_dim     
     
     DO I=1, SS_TOT
       READ(INORRK,*) J, EVAL(I)
     ENDDO
     
     DO I=1, SS_DIM
       DO J=1, SS_DIM
         GAMMA2(I,J) = 0._SP
       ENDDO
       GAMMA2(I,I) = EVAL(I)*RRK_PSCALE
     ENDDO
     CLOSE(INORRK)
   
! READ THE LINEAR TRANSITION MATRIX     
     FILENAME = TRIM(OUTPUT_DIR)//'/rrktemp/Amatr.dat' 
     OPEN(INORRK,FILE=TRIM(FILENAME),FORM='UNFORMATTED')
     DO J=1, SS_DIM
       READ(INORRK) (TRANS(I,J), I=1, SS_DIM)    
     ENDDO 
     CLOSE(INORRK)
     
! READ THE OBSERVATION NUMBER AND LOCATION   
     CALL READOBS(STLOC,NLOC)
     
     print *, 'obs=', STLOC(1), NLOC
   
! SET RALPHA   
     DO I=1, SS_DIM
       DO J=1, SS_DIM
         TRANS2(I,J)=TRANS(J,I)
       ENDDO
!      WRITE(IPT,*) 'Amatr:', I, TRANS(I,I)     
     ENDDO

     print *, 'debug - 4'
       
! CALCULATE OBSERVATION MATRIX (EOF TO OBS TRANSFORMATION) H*D*EOFs
     DO J=1, SS_DIM
       DO I=1, STDIM  
          STTEMP1(I) = SEOF(I,J)
       ENDDO
       
       DO I=1, NLOC         
         HEOF(I,J)=STTEMP1(STLOC(I))*STSD(STLOC(I))
        WRITE(IPT,*) 'Heof:', I,J,HEOF(I,J)        
       ENDDO
     ENDDO
    

     FILENAME = TRIM(OUTPUT_DIR)//'/rrktemp/ObsOp.dat'
     OPEN(INORRK,FILE=FILENAME, FORM='UNFORMATTED')
     DO J=1, SS_DIM
       WRITE(INORRK) (HEOF(I,J), I=1, NLOC)    
     ENDDO 
     CLOSE(INORRK)     
     FILENAME = TRIM(OUTPUT_DIR)//'/rrktemp/Heof.txt'
     OPEN(INORRK,FILE=FILENAME, FORM='FORMATTED')
     DO J=1, SS_DIM
       DO I=1, NLOC
         WRITE(INORRK,*) I,J,HEOF(I,J)
       ENDDO    
     ENDDO 
     CLOSE(INORRK)

! OUTPUT RALPHA MATRIX, WHICH IS TRANSPOSE OF TRANS
     FILENAME = TRIM(OUTPUT_DIR)//'/rrktemp/Alpha.dat'
     OPEN(INORRK,FILE=FILENAME, FORM='UNFORMATTED')
     DO J=1, SS_DIM
       WRITE(INORRK) (TRANS2(I,J), I=1, SS_DIM)    
     ENDDO 
     CLOSE(INORRK)     

! OUTPUT GAMMA MATRIX
     FILENAME = TRIM(OUTPUT_DIR)//'/rrktemp/Gamma.dat'
     OPEN(INORRK,FILE=FILENAME, FORM='UNFORMATTED')
     DO J=1, SS_DIM
       WRITE(INORRK) (GAMMA2(I,J), I=1, SS_DIM)    
     ENDDO 
     CLOSE(INORRK) 
     
! CALCULATE REPRESENTATIVE OBSERVATION ERROR DIRECTLY USING EOFs IN UNRESOLVED SUBSPACE
     DO I=SS_DIM+1, SS_TOT
!       CALL READEOF(I,2)
       DO J=1, STDIM  
          STTEMP1(J) = SEOF(J,I)
       ENDDO
       
       DO J=1, NLOC
          HU(J,I) = STTEMP1(STLOC(J))*STSD(STLOC(J))       
       ENDDO
     ENDDO
     DO I=1, NLOC
       RRKSUM=0.0_DP
       DO J=SS_DIM+1, SS_TOT
         HUL(I,J)=RRKSUM+HU(I,J)*EVAL(J)
       ENDDO
     ENDDO
     
     DO I=1, NLOC
       DO J=1, NLOC
          RRKSUM=0.0_DP
	  DO K=SS_DIM+1,SS_TOT
             RRKSUM = RRKSUM + HUL(I,K)*HU(J,K)
          ENDDO
          R(I,J) = RRKSUM
       ENDDO
       
! ASSUME THAT THE MEASUREMENT ERROR IS 1% OF THE ONE STANDARD DEVIATION OF THE CONTROL RUN
       IF(RRK_OPTION == 0) THEN
          R(I,I) = R(I,I) + (STSD(STLOC(I))*0.01_DP)**2.0_DP
       ELSE
          R(I,I) = R(I,I) + (STSD(STLOC(I))/DSQRT(DBLE(KBM1))*0.01_DP)**2.0_DP	  
       ENDIF
     ENDDO
   
     DO I=1, NLOC
       DO J=1, I-1
         R(I,J) = 0.5_DP*(R(I,J)+R(J,I))
         R(J,I) = R(I,J)
       ENDDO
     ENDDO
     FILENAME = TRIM(OUTPUT_DIR)//'/rrktemp/Robs.dat'
     OPEN(INORRK,FILE=FILENAME, FORM='UNFORMATTED')
     DO J=1, NLOC
       WRITE(INORRK) (R(I,J), I=1, NLOC)    
     ENDDO 
     CLOSE(INORRK)  
     FILENAME = TRIM(OUTPUT_DIR)//'/rrktemp/Robs.txt'
     OPEN(INORRK,FILE=FILENAME, FORM='FORMATTED')
     DO J=1, NLOC
       DO I=1, NLOC
         WRITE(INORRK,*) I,J, R(I,J)    
       ENDDO	 
     ENDDO 
     CLOSE(INORRK)   
     
! INVERT OBSERVATION ERROR COVARIANCE     
     CALL DGETRF(NLOC,NLOC,R,RRK_NOBSMAX,IPIV2,RCODE)
     IF(RCODE/=0) WRITE(IPT,*) 'error in computing LU factorization'     
     call DGETRI(NLOC,R,RRK_NOBSMAX,IPIV2,WORK2,LWORK2,RCODE)
     IF(RCODE/=0) WRITE(IPT,*) 'error in inverting the matrix'

! FORM BETA MATRIX = H_T*R**(-1)*H, WHERE H IS A NORMALIZED MATRIX H = H_ORIG*S
     CALL DGEMM('t','n',SS_DIM,NLOC,NLOC,1.0d0,HEOF,RRK_NOBSMAX,R,RRK_NOBSMAX,0.0d0,HTR,SS_DIM) 
     CALL DGEMM('n','n',SS_DIM,SS_DIM,NLOC,1.0d0,HTR,SS_DIM,HEOF,RRK_NOBSMAX,0.0d0,BETA2,SS_DIM)    
     FILENAME = TRIM(OUTPUT_DIR)//'/rrktemp/Beta.dat'
     OPEN(INORRK,FILE=FILENAME, FORM='UNFORMATTED')
     DO J=1, SS_DIM
       WRITE(INORRK) (BETA2(I,J), I=1, SS_DIM)    
     ENDDO 
     CLOSE(INORRK)     
     
     WRITE(IPT,*) 'Running the doubling algorithm......'
   
     DO IT =1, ILAST 
       
! READ GAMMA FROM FILE
       FILENAME =  TRIM(OUTPUT_DIR)//'/rrktemp/Gamma.dat'    
       OPEN(INORRK,FILE=TRIM(FILENAME),FORM='UNFORMATTED')
       DO J=1, SS_DIM
         READ(INORRK) (GAMMA(I,J), I=1, SS_DIM)    
       ENDDO 
       CLOSE(INORRK)
                    
! READ BETA FROM FILE
       FILENAME =  TRIM(OUTPUT_DIR)//'/rrktemp/Beta.dat'    
       OPEN(INORRK,FILE=TRIM(FILENAME),FORM='UNFORMATTED')
       DO J=1, SS_DIM
         READ(INORRK) (BETA2(I,J), I=1, SS_DIM)    
       ENDDO 
       CLOSE(INORRK)
      
! COMPUTE EYE + BETA*GAMMA = STORE TEMPORARILY IN RALPHA (RALPHA = EYE + BETA*GAMMA)
         
       DO I=1, SS_DIM
         DO J=1, SS_DIM
           RALPHA(I,J)=0.0_DP
         ENDDO
         RALPHA(I,I)=1.0_DP
       ENDDO     
       CALL DGEMM('n','n',SS_DIM,SS_DIM,SS_DIM,1.0d0,BETA2,SS_DIM,GAMMA,SS_DIM,1.0d0,RALPHA,SS_DIM)
     
! COMPUTE INVERSE OF ABOBE RALPHA (=EYE+BETA*GAMMA) BY LAPACK ROUTINES AND STORE IT IN BETA (BETA=(EYE+BETA*GAMMA)**(-1))
       CALL DGETRF(SS_DIM,SS_DIM,RALPHA,SS_DIM,IPIV,RCODE)
  
       IF(RCODE/=0) WRITE(IPT,*) 'error in computing LU factorization'
       CALL DGETRI(SS_DIM,RALPHA,SS_DIM,IPIV,WORK,LWORK,RCODE)
       IF(RCODE/=0) WRITE(IPT,*) 'error in inverting the matrix'
       DO I=1, SS_DIM
         DO J=1, SS_DIM
           BETA2(I,J) = RALPHA(I,J)
	 ENDDO
       ENDDO
     
! READ RALPHA FROM FILE
       FILENAME =  TRIM(OUTPUT_DIR)//'/rrktemp/Alpha.dat'    
       OPEN(INORRK,FILE=TRIM(FILENAME),FORM='UNFORMATTED')
       DO J=1, SS_DIM
         READ(INORRK) (RALPHA(I,J), I=1, SS_DIM)    
       ENDDO 
       CLOSE(INORRK)

! COMPUTE PRODUCT (EYE+BETA*GAMMA)**(-1)*RALPHA (GAMMA = BETA'*RALPHA)
       CALL DGEMM('n','n',SS_DIM,SS_DIM,SS_DIM,1.0d0,BETA2,SS_DIM,RALPHA,SS_DIM,0.0d0,GAMMA,SS_DIM)  
          
! OUTPUT THIS PRODUCT TO TEMPORARY FILE
       FILENAME = TRIM(OUTPUT_DIR)//'/rrktemp/temp.02'
       OPEN(INORRK,FILE=TRIM(FILENAME),FORM='UNFORMATTED')
       DO J=1, SS_DIM
         WRITE(INORRK) (GAMMA(I,J), I=1, SS_DIM)    
       ENDDO 
       CLOSE(INORRK)     

! COMPUTE PRODUCT RALPHA*(EYE+BETA*GAMMA)**(-1) AND STORE IN GAMMA (GAMMA = RALPHA*BETA')
       CALL DGEMM('n','n',SS_DIM,SS_DIM,SS_DIM,1.0d0,RALPHA,SS_DIM,BETA2,SS_DIM,0.0d0,GAMMA,SS_DIM) 
   
! READ BACK OLD FILE
       FILENAME = TRIM(OUTPUT_DIR)//'/rrktemp/Beta.dat'
       OPEN(INORRK,FILE=TRIM(FILENAME),FORM='UNFORMATTED')
       DO J=1, SS_DIM
         READ(INORRK) (Beta2(I,J), I=1, SS_DIM)    
       ENDDO 
       CLOSE(INORRK)  

! COMPUTE PRODUCT RALPHA*(EYE+BETA*GAMMA)**(-1)*BETA (RALPHA' = GAMMA'*BETA')
       CALL DGEMM('n','n',SS_DIM,SS_DIM,SS_DIM,1.0d0,GAMMA,SS_DIM,BETA2,SS_DIM,0.0d0,RALPHA,SS_DIM)

! READ BACK OLD RALPHA AND PUT IN GAMMA
       FILENAME = TRIM(OUTPUT_DIR)//'/rrktemp/Alpha.dat'
       OPEN(INORRK,FILE=TRIM(FILENAME),FORM='UNFORMATTED')
       DO J=1, SS_DIM
         READ(INORRK) (GAMMA(I,J), I=1, SS_DIM)    
       ENDDO 
       CLOSE(INORRK)  

! COMPUTE BETA+RALPHA*(EYE+BETA*GAMMA)**(-1)*BETA*RALPHA_T (BETA = BETA+RALPHA'*GAMMA_T)
       CALL DGEMM('n','t',SS_DIM,SS_DIM,SS_DIM,1.0d0,RALPHA,SS_DIM,GAMMA,SS_DIM,1.0d0,BETA2,SS_DIM)

! MAKE SURE BETA IS SYMMETRIC
       DO I=1, SS_DIM
          DO J=1, I-1
             BETA2(I,J) = 0.5_DP*(BETA2(I,J)+BETA2(J,I))
             BETA2(J,I) = BETA2(I,J) 
	  ENDDO
       ENDDO

! SAVE THIS NEW BETA TO FILE
       FILENAME = TRIM(OUTPUT_DIR)//'/rrktemp/Beta.dat'
       OPEN(INORRK,FILE=FILENAME, FORM='UNFORMATTED')
         DO J=1, SS_DIM
         WRITE(INORRK) (BETA2(I,J), I=1, SS_DIM)    
       ENDDO 
       CLOSE(INORRK)      

! READ BACK OLD GAMMA
       FILENAME = TRIM(OUTPUT_DIR)//'/rrktemp/Gamma.dat'
       OPEN(INORRK,FILE=FILENAME, FORM='UNFORMATTED')
       DO J=1, SS_DIM
         READ(INORRK) (GAMMA(I,J), I=1, SS_DIM)    
       ENDDO 
       CLOSE(INORRK)     

! READ TEMPORARY FILE INTO RALPHA
       FILENAME = TRIM(OUTPUT_DIR)//'/rrktemp/temp.02'
       OPEN(INORRK,FILE=FILENAME, FORM='UNFORMATTED')
       DO J=1, SS_DIM
         READ(INORRK) (RALPHA(I,J), I=1, SS_DIM)    
       ENDDO 
       CLOSE(INORRK)       

! COMPUTE GAMMA*(EYE+BETA*GAMMA)**(-1)*RALPHA (BETA'=GAMMA*RALPHA')
       CALL DGEMM('n','n',SS_DIM,SS_DIM,SS_DIM,1.0d0,GAMMA,SS_DIM,RALPHA,SS_DIM,0.0d0,BETA2,SS_DIM)

! READ BACK OLD RALPHA
       FILENAME = TRIM(OUTPUT_DIR)//'/rrktemp/Alpha.dat'
       OPEN(INORRK,FILE=FILENAME, FORM='UNFORMATTED')
       DO J=1, SS_DIM
         READ(INORRK) (RALPHA(I,J), I=1, SS_DIM)    
       ENDDO 
       CLOSE(INORRK)     

! COMPUTE GAMMA + RALPHA_T*GAMMA*(EYE+BETA*GAMMA)**(-1)*RALPHA (GAMMA = GAMMA+RALPHA_T*BETA')
       CALL DGEMM('t','n',SS_DIM,SS_DIM,SS_DIM,1.0d0,RALPHA,SS_DIM,BETA2,SS_DIM,1.0d0,GAMMA,SS_DIM)

! ENSURE GAMMA IS SYMMETRIC
       DO I=1, SS_DIM
          DO J=1, I-1
             GAMMA(I,J) = 0.5_DP*(GAMMA(I,J)+GAMMA(J,I))
             GAMMA(J,I) = GAMMA(I,J) 
	  ENDDO
       ENDDO

! SAVE THIS NEW GAMMA TO FILE
       FILENAME = TRIM(OUTPUT_DIR)//'/rrktemp/Gamma.dat'
       OPEN(INORRK,FILE=FILENAME, FORM='UNFORMATTED')
       DO J=1, SS_DIM
         WRITE(INORRK) (GAMMA(I,J), I=1, SS_DIM)    
       ENDDO 
       CLOSE(INORRK)      
       FILENAME = TRIM(OUTPUT_DIR)//'/rrktemp/Pfctr.txt'
       OPEN(INORRK,FILE=FILENAME, FORM='FORMATTED')
         DO J=1, SS_DIM
           DO I=1, SS_DIM
	      WRITE(INORRK,*) I,J,GAMMA(I,J)    
           ENDDO
	 ENDDO 
       CLOSE(INORRK)

! READ TEMPORARY FILE INTO GAMMA   
       FILENAME = TRIM(OUTPUT_DIR)//'/rrktemp/temp.02'
       OPEN(INORRK,FILE=FILENAME, FORM='UNFORMATTED')
       DO J=1, SS_DIM
         READ(INORRK) (GAMMA(I,J), I=1, SS_DIM)    
       ENDDO 
       CLOSE(INORRK)
   
! COMPUTE PRODUCT RALPHA*(EYE+BETA*GAMMA)**(-1)*RALPHA (BETA = RALPHA*GAMMA)  
       CALL DGEMM('n','n',SS_DIM,SS_DIM,SS_DIM,1.0d0,RALPHA,SS_DIM,GAMMA,SS_DIM,0.0d0,BETA2,SS_DIM)
   
! SAVE THIS NEW RALPHA TO FILE
       FILENAME = TRIM(OUTPUT_DIR)//'/rrktemp/Alpha.dat'
       OPEN(INORRK,FILE=FILENAME, FORM='UNFORMATTED')
         DO J=1, SS_DIM
         WRITE(INORRK) (BETA2(I,J), I=1, SS_DIM)    
       ENDDO 
       CLOSE(INORRK)

!===================================================================================================
! END OF ITERATION FOR THE DOUBLING ALGORITHM          
     ENDDO   
!===================================================================================================

     WRITE(IPT,*) 'Setting up the Kalman gain matrix in the reduced space......'   

! READ IN PF FROM DOUBLING ALGORITHM 
     FILENAME = TRIM(OUTPUT_DIR)//'/rrktemp/Gamma.dat'
     OPEN(INORRK,FILE=FILENAME, FORM='UNFORMATTED')
     DO J=1, SS_DIM
       READ(INORRK) (GAMMA(I,J), I=1, SS_DIM)    
     ENDDO 
     CLOSE(INORRK)
   
! READ IN OBSERVATION OPERATOR   
     FILENAME = TRIM(OUTPUT_DIR)//'/rrktemp/ObsOp.dat'
     OPEN(INORRK,FILE=FILENAME, FORM='UNFORMATTED')
     DO J=1, SS_DIM
       READ(INORRK) (HEOF(I,J), I=1, NLOC)    
     ENDDO 
     CLOSE(INORRK)
      
! READ OBSERVATION ERROR COVARIANCE      
     FILENAME = TRIM(OUTPUT_DIR)//'/rrktemp/Robs.dat'
     OPEN(INORRK,FILE=FILENAME, FORM='UNFORMATTED')
     DO J=1, NLOC
       READ(INORRK) (R(I,J), I=1, NLOC)    
     ENDDO 
     CLOSE(INORRK)
   
! COMPUTE P*H^T (PHT = GAMMA * HEOF^T)   
     CALL DGEMM('n','t',SS_DIM,NLOC,SS_DIM,1.0d0,GAMMA,SS_DIM,HEOF,RRK_NOBSMAX,0.0d0,PHT,SS_DIM)

! COMPUTE H*P*H^T (R2 = H*P*HT)  
     CALL DGEMM('n','n',NLOC,NLOC,SS_DIM,1.0d0,HEOF,RRK_NOBSMAX,PHT,SS_DIM,0.0d0,R2,RRK_NOBSMAX)

! OUPUT BOTH FORCAST AND OBSERVATION ERROR STANDARD DEVIATION IN ASCII
     FILENAME =  TRIM(OUTPUT_DIR)//'/rrktemp/Fctobserr.txt'  
     OPEN(INORRK,FILE=FILENAME, FORM='FORMATTED')
     DO I=1, NLOC
       WRITE(INORRK,*) STLOC(I),DSQRT(R2(I,I)),DSQRT(R(I,I))
     ENDDO
     CLOSE(INORRK)
   
! ADD R + (H P H^t) ---> R   
     DO I = 1, NLOC
        DO J=1, NLOC
           R(J,I) = R(J,I) + R2(J,I)
        ENDDO
     ENDDO
   
! INVERT R (H*P*H^T+R)   
     CALL DGETRF(NLOC,NLOC,R,RRK_NOBSMAX,IPIV2,RCODE)
     IF(RCODE/=0) WRITE(IPT,*) 'error in computing LU factorization'
     CALL DGETRI(NLOC,R,RRK_NOBSMAX,IPIV2,WORK2,LWORK2,RCODE)
     IF(RCODE/=0) WRITE(IPT,*) 'error in inverting the matrix'
     
! COMPUTE KALMAN GAIN: K = P*HT*(H*P*H^T+R)^(-1)   
     CALL DGEMM('n','n',SS_DIM,NLOC,NLOC,1.0d0,PHT,SS_DIM,R,RRK_NOBSMAX,0.0d0,KAL,SS_DIM)
   
! OUTPUT RESULT   
     OPEN(INORRK,FILE=KALNAME, FORM='UNFORMATTED') 
     DO J=1, NLOC
       WRITE(INORRK) (KAL(I,J), I=1, SS_DIM)    
     ENDDO 
     CLOSE(INORRK)     
     FILENAME = TRIM(OUTPUT_DIR)//'/rrktemp/rrK.txt'
     OPEN(INORRK,FILE=FILENAME, FORM='FORMATTED') 
     DO J=1, NLOC
       DO I=1, SS_DIM
	 WRITE(INORRK,*) I,J,KAL(I,J)    
       ENDDO
     ENDDO 
     CLOSE(INORRK)
     
! OUTPUT KALMAN GAIN MATRIX IN THE FULL SPACE: D*ER*KAL
!     DO J=1, SS_DIM
!       CALL READEOF(J,2)
!       DO I=1, STDIM
!         SEOF(I,J)=STTEMP1(I)
!       ENDDO
!     ENDDO
     
     DO J=1, NLOC
       WRITE(STRCYC,'(I4.4)') J
       DO I=1, STDIM
          RRKSUM = 0.0_DP
          DO K=1, SS_DIM
             RRKSUM = RRKSUM + STSD(I)*SEOF(I,K)*KAL(K,J) 
          ENDDO 
	  STTEMP1(I) = RRKSUM
       ENDDO
       FILENAME = TRIM(OUTPUT_DIR)//'/rrktemp/'//'BigK'//STRCYC//'.txt'     
       OPEN(INORRK,FILE=FILENAME,FORM='FORMATTED')
       DO I=1, STDIM
          WRITE(INORRK,*) STTEMP1(I)
       ENDDO
       CLOSE(INORRK)
     ENDDO
     
! COMPUTE ANALYSIS ERROR COVARIANCE: P_AN1 = (I-K*H)*P (JUST AS DIAGNOSIS)
     CALL DGEMM('n','t',SS_DIM,SS_DIM,NLOC,-1.0d0,KAL,SS_DIM,PHT,SS_DIM,1.0d0,GAMMA,SS_DIM)
     FILENAME = TRIM(OUTPUT_DIR)//'/rrktemp/Panlr.txt'
     OPEN(INORRK,FILE=FILENAME, FORM='FORMATTED') 
     DO J=1, SS_DIM
       DO I=1, SS_DIM
	 WRITE(INORRK,*) I,J,GAMMA(I,J)    
       ENDDO
     ENDDO 
     CLOSE(INORRK)
     
     CALL DGEMM('n','n',SS_DIM,SS_DIM,NLOC,1.0d0,KAL,SS_DIM,HEOF,RRK_NOBSMAX,0.0d0,RALPHA,SS_DIM)
     
! COMPUTE (I-KH)M AND ITS SINGULAR VALUE     
     DO I=1, SS_DIM
       DO J=1, SS_DIM
         IF(I/=J) THEN
            RALPHA(I,J) = (-1.0_DP)*RALPHA(I,J)
         ELSE
	    RALPHA(I,I) = 1.0_DP - RALPHA(I,I)
         ENDIF
       ENDDO
     ENDDO
     
     FILENAME = TRIM(OUTPUT_DIR)//'/rrktemp/I_KH.txt'
     OPEN(INORRK,FILE=FILENAME, FORM='FORMATTED') 
     DO J=1, SS_DIM
       DO I=1, SS_DIM
	 WRITE(INORRK,*) I,J,RALPHA(I,J)    
       ENDDO
     ENDDO 
     CLOSE(INORRK)     
     
     FILENAME = TRIM(OUTPUT_DIR)//'/rrktemp/Amatr.dat'
     OPEN(INORRK,FILE=FILENAME, FORM='UNFORMATTED')
     DO J=1, SS_DIM
       READ(INORRK) (TRANS(I,J), I=1, SS_DIM)    
     ENDDO 
     CLOSE(INORRK)     
     CALL DGEMM('n','n',SS_DIM,SS_DIM,SS_DIM,1.0d0,RALPHA,SS_DIM,TRANS,SS_DIM,0.0d0,RALPHA,SS_DIM)
    
     FILENAME = TRIM(OUTPUT_DIR)//'/rrktemp/M_KHM.txt'
     OPEN(INORRK,FILE=FILENAME, FORM='FORMATTED') 
     DO J=1, SS_DIM
       DO I=1, SS_DIM
	 WRITE(INORRK,*) I,J,RALPHA(I,J)    
       ENDDO
     ENDDO 
     CLOSE(INORRK) 
     
     DEALLOCATE(STLOC,IPIV,IPIV2,KAL,HEOF,R,HTR,R2,PHT,RALPHA,TRANS2)
!     DEALLOCATE(BETA2,GAMMA,GAMMA2,WORK,WORK2,SEOF,STSD,HU,HUL,EVAL)
     DEALLOCATE(BETA2,GAMMA,GAMMA2,WORK,WORK2,HU,HUL,EVAL)
!     DEALLOCATE(STTEMP0,STTEMP1,STEOF,SDEOF,TRANS,RRKU,RRKV,RRKEL,RRKT,RRKS)
   ENDIF
   
   RETURN
   END SUBROUTINE RRK_RRK


! UTILITIES PROGRAMS
 
  SUBROUTINE RRK_ALLOC_VAR

   USE CONTROL
   USE RRKVAL
   IMPLICIT NONE
   
   INTEGER STDIM
   INTEGER SS_DIM
   
   STDIM = 0
   IF(EL_ASSIM) STDIM = STDIM + MGL
   IF(UV_ASSIM) STDIM = STDIM + 2*NGL*KBM1
   IF(T_ASSIM)  STDIM = STDIM + MGL*KBM1
   IF(S_ASSIM)  STDIM = STDIM + MGL*KBM1
   
   SS_DIM = RRK_NEOF
   
! ALLOCATE ARRYS

   IF(.not.ALLOCATED(STTEMP0))  ALLOCATE(STTEMP0(STDIM))  ;STTEMP0   = ZERO
   IF(.not.ALLOCATED(STTEMP1))  ALLOCATE(STTEMP1(STDIM))  ;STTEMP1   = ZERO
   IF(.not.ALLOCATED(STEOF))    ALLOCATE(STEOF(STDIM))    ;STEOF     = ZERO
   IF(.not.ALLOCATED(SDEOF))    ALLOCATE(SDEOF(STDIM))    ;SDEOF     = ZERO
   IF(.not.ALLOCATED(TRANS))    ALLOCATE(TRANS(SS_DIM,SS_DIM)) ;TRANS     = ZERO
   IF(.not.ALLOCATED(RRKEL))    ALLOCATE(RRKEL(0:MGL))    ;RRKEL     = ZERO
   IF(.not.ALLOCATED(RRKU))     ALLOCATE(RRKU(0:NGL,KB))  ;RRKU      = ZERO
   IF(.not.ALLOCATED(RRKV))     ALLOCATE(RRKV(0:NGL,KB))  ;RRKV      = ZERO
   IF(.not.ALLOCATED(RRKT))     ALLOCATE(RRKT(0:MGL,KB))  ;RRKT      = ZERO
   IF(.not.ALLOCATED(RRKS))     ALLOCATE(RRKS(0:MGL,KB))  ;RRKS      = ZERO
   
! END ALLOCATION
   
   RETURN
  END SUBROUTINE RRK_ALLOC_VAR

!------------------------------------------------------------------------------|
!  PERTURB THE BASE STATE IN THE EIGENVECTOR DIRECTIONS                        |
!------------------------------------------------------------------------------|   
   
  SUBROUTINE PERTURB

#  if defined(WET_DRY)
     USE MOD_WD
#  endif
   USE LIMS
   USE ALL_VARS
   USE RRKVAL
   USE CONTROL
#  if defined (MULTIPROCESSOR)
     USE MOD_PAR
#  endif
   IMPLICIT NONE
    INTEGER IERR 
    INTEGER SS_DIM
    INTEGER STDIM
    INTEGER I_EOF
    INTEGER II,I,J
    INTEGER IDUMMY
    INTEGER I_START
    CHARACTER(LEN=80) IFILE
    CHARACTER(LEN=80) IFILE2
    CHARACTER(LEN=80) TEMP1FILE
    CHARACTER(LEN=80) TEMP2FILE
    CHARACTER(LEN=4)  FEOF
      
    STDIM = 0
    IF(EL_ASSIM) STDIM = STDIM + MGL
    IF(UV_ASSIM) STDIM = STDIM + 2*NGL*KBM1
    IF(T_ASSIM)  STDIM = STDIM + MGL*KBM1
    IF(S_ASSIM)  STDIM = STDIM + MGL*KBM1
    
    SS_DIM = RRK_NEOF

   ALLOCATE(STTEMP0(STDIM))                ;STTEMP0   = ZERO
   ALLOCATE(STTEMP1(STDIM))                ;STTEMP1   = ZERO
   ALLOCATE(STEOF(STDIM))                  ;STEOF     = ZERO
   ALLOCATE(SDEOF(STDIM))                  ;SDEOF     = ZERO
   ALLOCATE(TRANS(SS_DIM,SS_DIM))          ;TRANS     = ZERO        

    IDUMMY  = 0
    IF(EL_ASSIM) THEN
      DO I=1, MGL
        IDUMMY = IDUMMY + 1
        SDEOF(IDUMMY) = ELSD
      ENDDO
    ENDIF
    IF(UV_ASSIM) THEN
      DO I=1, 2*KBM1
        DO J=1, NGL 
          IDUMMY = IDUMMY + 1
          SDEOF(IDUMMY)=USD*DSQRT(DBLE(KBM1))     
        ENDDO
      ENDDO
    ENDIF
    IF(T_ASSIM) THEN
      DO I=1, KBM1
        DO J=1, MGL 
          IDUMMY = IDUMMY + 1
          SDEOF(IDUMMY)=TSD
        ENDDO
      ENDDO    
    ENDIF
    IF(S_ASSIM) THEN
      DO I=1, KBM1
        DO J=1, MGL 
          IDUMMY = IDUMMY + 1
          SDEOF(IDUMMY)=SSD
        ENDDO
      ENDDO    
    ENDIF
      print *, 'before rrk_rrk read restart, file name:',nc_start%fname ! should be marked
 
    print *, 'call readrestart' 
    CALL READRESTART(NC_START,REF_START_TIME)

!   IF(.not.ALLOCATED(RRKEL))    ALLOCATE(RRKEL3(0:MGL))    ;RRKEL3     = ZERO
!   IF(.not.ALLOCATED(RRKU))     ALLOCATE(RRKU3(0:NGL,KB))  ;RRKU3      = ZERO
!   IF(.not.ALLOCATED(RRKV3))     ALLOCATE(RRKV3(0:NGL,KB))  ;RRKV3      = ZERO
!   IF(.not.ALLOCATED(RRKT3))     ALLOCATE(RRKT3(0:MGL,KB))  ;RRKT3      = ZERO
!   IF(.not.ALLOCATED(RRKS3))     ALLOCATE(RRKS3(0:MGL,KB))  ;RRKS3      = ZERO
   ALLOCATE(RRKU(NGL,KB))            ;RRKU    = ZERO
   ALLOCATE(RRKV(NGL,KB))            ;RRKV    = ZERO
   ALLOCATE(RRKEL(MGL))                ;RRKEL   = ZERO
   ALLOCATE(RRKT(MGL,KB))            ;RRKT    = ZERO
   ALLOCATE(RRKS(MGL,KB))            ;RRKS    = ZERO

print *, allocated(rrku)

#  if defined (MULTIPROCESSOR)
if (par) then
CALL ACOLLECT(MYID,MSRID,NPROCS,NMAP,EL, RRKEL)
CALL ACOLLECT(MYID,MSRID,NPROCS,EMAP,U, RRKU)
CALL ACOLLECT(MYID,MSRID,NPROCS,EMAP,V, RRKV)
CALL ACOLLECT(MYID,MSRID,NPROCS,NMAP,T1, RRKT)
CALL ACOLLECT(MYID,MSRID,NPROCS,NMAP,S1, RRKS)
CALL MPI_BCAST(RRKEL,MGL,MPI_F,0,MPI_FVCOM_GROUP,IERR)
CALL MPI_BCAST(RRKU,NGL*KBM1,MPI_F,0,MPI_FVCOM_GROUP,IERR)
CALL MPI_BCAST(RRKV,NGL*KBM1,MPI_F,0,MPI_FVCOM_GROUP,IERR)
CALL MPI_BCAST(RRKT,MGL*KBM1,MPI_F,0,MPI_FVCOM_GROUP,IERR)
CALL MPI_BCAST(RRKS,MGL*KBM1,MPI_F,0,MPI_FVCOM_GROUP,IERR)

else
   RRKEL=EL
   RRKU=U
   RRKV=V
   RRKT=T
   RRKS=S
end if
#  elif
   RRKEL=EL
   RRKU=U
   RRKV=V
   RRKT=T
   RRKS=S
#  endif
    print *, 'call gr2st'
    CALL GR2ST(0)
    print *, 'finish gr2st'
    
! READ THE EOF AND PERTURB THE BASE STATE IN THIS DIRECTION
  if (msr) then
     print *, 'msr id' , myid
     open(343,file='seof',status='replace',form='unformatted')
     write(343) seof
     close(343)
  else
     print *, 'no msr id', myid
!    ALLOCATE(SEOF(STDIM,SS_DIM))   ! because the seof is calculated only the master cpu in rrk_ref
  end if
  
!write(myid+4000,*) myid 
# if defined (MULTIPROCESSOR)
IF(PAR)CALL MPI_BARRIER(MPI_FVCOM_GROUP,IERR)
     open(343,file='seof',status='old',form='unformatted',action='read')
     read(343) seof
     close(343)
IF(PAR)CALL MPI_BARRIER(MPI_FVCOM_GROUP,IERR)
! IF(PAR)CALL MPI_BCAST(SEOF,STDIM*SS_DIM,MPI_DOUBLE_precision,0,MPI_FVCOM_GROUP,IERR)
write(myid+3000,*) myid
# endif 
write(myid+2000,*) seof
   DO I_EOF=1, SS_DIM
!       CALL READEOF(I_EOF,1)
       DO I=1, STDIM  
          STEOF(I) = SEOF(I,I_EOF)
       ENDDO

! PERTURB THE BASE STATE IN THE DIRECTION OF THE I_EOF'TH EOF AND STORE IT IN THE MP1FILE (FOR RESTART)
       DO I=1, STDIM
          STTEMP1(I) = STTEMP0(I) + STEOF(I)*SDEOF(I)*DBLE(RRK_PSIZE)
           write(500+myid,*) STTEMP0(I), STEOF(I),SDEOF(I)
       ENDDO

     print *, 'call st2gr'
       CALL ST2GR
        do i=1,mgl
           write(400+myid,'(30f18.8)') (rrkt(i,j),j=1,kb)
        end do
      
   IF(SERIAL) THEN
     EL =  RRKEL 
     U  =  RRKU
     V  =  RRKV
     T1 =  RRKT
     S1 =  RRKS
   ELSE 
#  if defined (MULTIPROCESSOR)
     DO I=1, M
        EL(I)=RRKEL(NGID(I))
     ENDDO

     DO I=1, N
       DO J=1, KB
         U(I,J)=RRKU(EGID(I),J)
         V(I,J)=RRKV(EGID(I),J)
       ENDDO
     ENDDO

     DO I=1, M
       DO J=1, KB
         T1(I,J)=RRKT(NGID(I),J)
         S1(I,J)=RRKS(NGID(I),J)
       ENDDO
     ENDDO
#  endif
   ENDIF


! WRITE THE PERTURBED STATE IN TEMP1FILE
     CALL WRITERESTART(NCF_RRKRE,I_EOF)

   ENDDO


  RETURN
  END SUBROUTINE PERTURB

!------------------------------------------------------------------------------|
!  CALCULATE THE LINEARIZED MODEL IN THE REDUCED SPACE                         |
!------------------------------------------------------------------------------|   
   
  SUBROUTINE MREDUCED

   USE MOD_WD
   USE LIMS
   USE CONTROL
   USE RRKVAL
   IMPLICIT NONE
    
   INTEGER SS_DIM
   INTEGER STDIM
   INTEGER I_EOF
   INTEGER I,J
   INTEGER IDUMMY
   INTEGER I_START   
   CHARACTER(LEN=80) IFILE
   CHARACTER(LEN=80) IFILE2
   CHARACTER(LEN=80) TEMPFILE
   CHARACTER(LEN=4)  FEOF
   REAL(DP) ::  SUM0

   STDIM = 0
   IF(EL_ASSIM) STDIM = STDIM + MGL
   IF(UV_ASSIM) STDIM = STDIM + 2*NGL*KBM1
   IF(T_ASSIM)  STDIM = STDIM + MGL*KBM1
   IF(S_ASSIM)  STDIM = STDIM + MGL*KBM1

   SS_DIM = RRK_NEOF
   
   IDUMMY = 0
   IF(EL_ASSIM) THEN
     DO I=1, MGL
       IDUMMY = IDUMMY + 1
       SDEOF(IDUMMY) = ELSD
     ENDDO
   ENDIF
   IF(UV_ASSIM) THEN
     DO I =1, 2*KBM1
       DO J=1, NGL
         IDUMMY = IDUMMY + 1
         IF(RRK_OPTION == 0) THEN
            SDEOF(IDUMMY) = USD
         ELSE
            SDEOF(IDUMMY) = USD*DSQRT(DBLE(KBM1))
         ENDIF
       ENDDO
     ENDDO
   ENDIF
   IF(T_ASSIM) THEN
     DO I =1, KBM1
       DO J=1, MGL
         IDUMMY = IDUMMY + 1
         SDEOF(IDUMMY) = TSD
       ENDDO
     ENDDO
   ENDIF
   IF(S_ASSIM) THEN
     DO I =1, KBM1
       DO J=1, MGL
         IDUMMY = IDUMMY + 1
         SDEOF(IDUMMY) = SSD
       ENDDO
     ENDDO
   ENDIF   
   
   CALL READRESTART(NC_START,REF_START_TIME)
   CALL GR2ST(0)  
   
   
! READ THE FORECAST WHICH IS PROPAGATED FROM THE PERTURBED STATE  

   CALL NC_OPEN(NCF_RRKFCT)
 
   DO I_EOF =1, SS_DIM   
! READ THE EACH PERTURBED STATE AT THE END OF LINEARUZATION TIME STEP

   CALL READRESTART(NCF_RRKFCT,I_EOF) 

   CALL GR2ST(1)   
   
! PROJECT  EVOLVED PERTURBATION ONTO EOFs
     DO I=1, SS_DIM
       DO J=1, STDIM  
          STEOF(J) = SEOF(J,I)
       ENDDO
       SUM0 = 0.0_DP
       DO J=1, STDIM
         SUM0 = SUM0 + STEOF(J)/SDEOF(J)*(STTEMP1(J)-STTEMP0(J))
       ENDDO
       TRANS(I,I_EOF) = SUM0/DBLE(RRK_PSIZE)
           
     ENDDO
! EDN OF LOOP OVER EACH EOF   
   ENDDO
   
  END SUBROUTINE MREDUCED 


!=====================================================================================/
!  CONVERT THE STATE FROM THE GRID SPACE TO THE STATE VECTOR SPACE                    /
!=====================================================================================/
  SUBROUTINE GR2ST(OPT)
 
   USE LIMS
   USE ALL_VARS
   USE RRKVAL
   IMPLICIT NONE
   
    INTEGER IDUMMY, I, J
    INTEGER OPT

    IDUMMY=0
    
    IF(EL_ASSIM) THEN
      DO I=1, MGL
        IDUMMY = IDUMMY + 1
        IF (OPT == 0) THEN
          STTEMP0(IDUMMY) = RRKEL(I)
        ELSE
          STTEMP1(IDUMMY) = RRKEL(I)
        ENDIF	 	
      ENDDO
    ENDIF

    IF(UV_ASSIM) THEN
      DO I=1, KBM1
        DO J=1, NGL
          IDUMMY = IDUMMY + 1
          IF (OPT == 0) THEN
	    STTEMP0(IDUMMY) = RRKU(J,I)

	  ELSE
	    STTEMP1(IDUMMY) = RRKU(J,I)
	  ENDIF    
        ENDDO
      ENDDO
      DO I=1, KBM1
        DO J=1, NGL
          IDUMMY = IDUMMY + 1
          IF (OPT == 0) THEN
	    STTEMP0(IDUMMY) = RRKV(J,I)
	  ELSE
	    STTEMP1(IDUMMY) = RRKV(J,I)
	  ENDIF    
        ENDDO
      ENDDO		
    ENDIF  
  
    IF(T_ASSIM) THEN
      DO I=1, KBM1
        DO J=1, MGL
          IDUMMY = IDUMMY + 1
          IF (OPT == 0) THEN
	    STTEMP0(IDUMMY) = RRKT(J,I)
	  ELSE
	    STTEMP1(IDUMMY) = RRKT(J,I)
	  ENDIF    
        ENDDO
      ENDDO
    ENDIF
  
    IF(S_ASSIM) THEN
      DO I=1, KBM1
        DO J=1, MGL
          IDUMMY = IDUMMY + 1
          IF (OPT == 0) THEN
	    STTEMP0(IDUMMY) = RRKS(J,I)
	  ELSE
	    STTEMP1(IDUMMY) = RRKS(J,I)
	  ENDIF    
        ENDDO
      ENDDO
    ENDIF

  RETURN
  END SUBROUTINE GR2ST

!=====================================================================================/
!  CONVERT THE STATE FROM THE STATE VECTOR SPACE TO THE GRID SPACE                    /
!=====================================================================================/
  SUBROUTINE ST2GR
 
   USE LIMS
   USE ALL_VARS
   USE RRKVAL
   IMPLICIT NONE
    
    INTEGER IDUMMY, I, J

    IDUMMY=0
    IF(EL_ASSIM) THEN
      DO I=1, MGL
        IDUMMY = IDUMMY + 1
        RRKEL(I) = STTEMP1(IDUMMY) 
      ENDDO
    ENDIF
    IF(UV_ASSIM) THEN
      DO I=1, KBM1
        DO J=1, NGL
          IDUMMY = IDUMMY + 1
          RRKU(J,I) = STTEMP1(IDUMMY) 
        ENDDO
      ENDDO
      DO I=1, KBM1
        DO J=1, NGL
          IDUMMY = IDUMMY + 1
          RRKV(J,I) = STTEMP1(IDUMMY)
        ENDDO
      ENDDO
    ENDIF
    IF(T_ASSIM) THEN
      DO I=1, KBM1
        DO J=1, MGL
          IDUMMY = IDUMMY + 1
          RRKT(J,I) = STTEMP1(IDUMMY) 
        ENDDO
      ENDDO
    ENDIF
    IF(S_ASSIM) THEN
      DO I=1, KBM1
        DO J=1, MGL
          IDUMMY = IDUMMY + 1
          RRKS(J,I) = STTEMP1(IDUMMY) 
        ENDDO
      ENDDO
    ENDIF
    
  RETURN
  END SUBROUTINE ST2GR


!=====================================================================================/
!  DETERMINE IF NODES/ELEMENTS ARE WET OR DRY                                         /
!=====================================================================================/
  SUBROUTINE WET_JUDGE_EL

   USE MOD_PREC
   USE ALL_VARS
#  if defined (MULTIPROCESSOR)
   USE MOD_PAR
#  endif
#  if defined(WET_DRY)
   USE MOD_WD
#  endif
   IMPLICIT NONE
   REAL(SP) :: DTMP
   INTEGER  :: ITA_TEMP
   INTEGER  :: I,IL,IA,IB,K1,K2,K3,K4,K5,K6

#  if defined(WET_DRY)

!
!--Determine If Node Points Are Wet/Dry Based on Depth Threshold---------------!
!
   ISWETN = 1
   DO I = 1, M
     DTMP = H(I) + EL(I)
     IF((DTMP - MIN_DEPTH) < 1.0E-5_SP) ISWETN(I) = 0
   END DO

!
!--Determine if Cells are Wet/Dry Based on Depth Threshold---------------------!
!
   ISWETC = 1
   DO I = 1, N
     DTMP =  MAX(EL(NV(I,1)),EL(NV(I,2)),EL(NV(I,3)))  + &
             MIN(  H(NV(I,1)),  H(NV(I,2)),  H(NV(I,3)))
     IF((DTMP - MIN_DEPTH) < 1.0E-5_SP) ISWETC(I) = 0
   END DO

!
!--A Secondary Condition for Nodal Dryness-(All Elements Around Node Are Dry)--!
!
   DO I = 1, M
     IF(SUM(ISWETC(NBVE(I,1:NTVE(I)))) == 0)  ISWETN(I) = 0
   END DO

!
!--Adjust Water Surface So It Does Not Go Below Minimum Depth------------------!
!
   EL = MAX(EL,-H + MIN_DEPTH)

!
!--Recompute Element Based Depths----------------------------------------------!
!
   DO I = 1, N
     EL1(I) = ONE_THIRD*(EL(NV(I,1))+EL(NV(I,2))+EL(NV(I,3)))
   END DO

!
!--Extend Element/Node Based Wet/Dry Flags to Domain Halo----------------------!
!
#  if defined (MULTIPROCESSOR)
   IF(PAR)THEN

     CALL AEXCHANGE(EC,MYID,NPROCS,ISWETC)
     CALL AEXCHANGE(NC,MYID,NPROCS,ISWETN)

   END IF

#  endif
#  endif

   RETURN

  END SUBROUTINE WET_JUDGE_EL

!!$!==============================================================================|
!!$!   DUMP WET/DRY FLAG DATA FOR RESTART                                         |
!!$!==============================================================================|
!!$
!!$   SUBROUTINE WD_DUMP_EL(INF,I_START)
!!$
!!$!------------------------------------------------------------------------------|
!!$
!!$   USE ALL_VARS
!!$#  if defined (MULTIPROCESSOR)
!!$   USE MOD_PAR
!!$#  endif
!!$#  if defined(WET_DRY)
!!$   USE MOD_WD
!!$#  endif
!!$
!!$   IMPLICIT NONE
!!$   INTEGER, ALLOCATABLE,DIMENSION(:) :: NTEMP1,NTEMP2
!!$   INTEGER I, INF
!!$   INTEGER I_START
!!$!==============================================================================|
!!$
!!$#  if defined(WET_DRY)
!!$
!!$   IF(MSR)THEN
!!$     REWIND(INF)
!!$     WRITE(INF,*) I_START
!!$     WRITE(INF,*) NGL,MGL
!!$   END IF
!!$
!!$   IF(SERIAL)THEN
!!$     WRITE(INF,*) (ISWETC(I), I=1,N)
!!$     WRITE(INF,*) (ISWETN(I), I=1,M)
!!$   ELSE
!!$   ALLOCATE(NTEMP1(NGL),NTEMP2(MGL))
!!$#  if defined (MULTIPROCESSOR)
!!$   CALL IGATHER(LBOUND(ISWETC,1),UBOUND(ISWETC,1),N,NGL,1,MYID,NPROCS,EMAP,ISWETC,NTEMP1)
!!$   CALL IGATHER(LBOUND(ISWETN,1),UBOUND(ISWETN,1),M,MGL,1,MYID,NPROCS,NMAP,ISWETN,NTEMP2)
!!$   IF(MSR)THEN
!!$     WRITE(INF,*) (NTEMP1(I), I=1,NGL)
!!$     WRITE(INF,*) (NTEMP2(I), I=1,MGL)
!!$   END IF
!!$   DEALLOCATE(NTEMP1,NTEMP2)
!!$#  endif
!!$   END IF
!!$
!!$   CLOSE(INF)
!!$
!!$#  endif
!!$
!!$   RETURN
!!$   END SUBROUTINE WD_DUMP_EL

  SUBROUTINE READOBS(STLOC,NLOC)
    
   USE LIMS
   USE CONTROL
   IMPLICIT NONE
   
     INTEGER ::  NUM  = 0 
     INTEGER ::  NLOC
     INTEGER ::  SWITCH = 0
     INTEGER ::  J,K
     INTEGER ::  IDUMMY
     INTEGER ::  TMP
     CHARACTER(LEN=80) FILENAME
     CHARACTER(LEN=24) HEADINFO
     INTEGER STLOC(RRK_NOBSMAX)
     INTEGER LAY(RRK_NOBSMAX)

     FILENAME = TRIM(INPUT_DIR)//"/"//trim(casename)//"_assim_rrkf.nml"
    
     OPEN(73,FILE=TRIM(FILENAME),FORM='FORMATTED')
     REWIND(73)

     NLOC = 0
     IDUMMY = 0
     STLOC = ZERO
     
 100 READ(73,'(A24)',END=200) HEADINFO
     IF(SWITCH/=1) THEN
       IF(HEADINFO=='!=== READ IN OBSERVATION') THEN
         SWITCH = 1
         GOTO 100
       ELSE
         GOTO 100
       ENDIF
     ENDIF 
     
     IF(TRIM(HEADINFO)=='!EL') THEN
       IF(EL_OBS) THEN
         READ(73,*) NUM
	 NLOC = NLOC + NUM
         IF(NLOC>RRK_NOBSMAX) THEN
           WRITE(IPT,*) 'not enough storage for observations:', 'Nloc=', Nloc, 'Nobsmax=', RRK_NOBSMAX
           CALL PSTOP
         ENDIF
	 READ(73,*)  (STLOC(K), K=1,NLOC)	
       ENDIF  
       
       IF(EL_ASSIM) THEN
	 IDUMMY = IDUMMY + MGL
       ENDIF
            
     ENDIF
     
     IF(TRIM(HEADINFO)=='!UV') THEN
       IF(UV_OBS) THEN       
         READ(73,*) NUM
	 NLOC = NLOC + NUM
         IF(NLOC+NUM>RRK_NOBSMAX) THEN
           WRITE(IPT,*) 'not enough storage for observations:', 'Nloc=', Nloc+num, 'Nobsmax=', RRK_NOBSMAX
           CALL PSTOP
         ENDIF
	 READ(73,*)  (STLOC(K), K=NLOC-NUM+1,NLOC)
	 READ(73,*)  (LAY(K),   K=NLOC-NUM+1,NLOC)
         DO K=NLOC-NUM+1, NLOC
	   STLOC(K)=STLOC(K)+IDUMMY+NGL*(LAY(K)-1)
	 ENDDO   
         IDUMMY = IDUMMY + NGL*KBM1
	 
	 NLOC = NLOC + NUM
	 DO K=NLOC-NUM+1, NLOC
	   STLOC(K)=STLOC(K-NUM)+NGL*KBM1+NGL*(LAY(K-NUM)-1)	   
	 ENDDO
       ENDIF
       
       IF(UV_ASSIM) THEN  
          IDUMMY = IDUMMY + NGL*KBM1
       ENDIF
     ENDIF
     
     IF(TRIM(HEADINFO)=='!T') THEN
       IF(T_OBS) THEN
         READ(73,*) NUM
	 NLOC = NLOC + NUM
         IF(NLOC>RRK_NOBSMAX) THEN
           WRITE(IPT,*) 'not enough storage for observations:', 'Nloc=', Nloc, 'Nobsmax=', RRK_NOBSMAX
           CALL PSTOP
         ENDIF
	 READ(73,*)  (STLOC(K), K=NLOC-NUM+1,NLOC)       
         READ(73,*)  (LAY(K),   K=NLOC-NUM+1,NLOC)      
         DO K=NLOC-NUM+1, NLOC
	   STLOC(K)=STLOC(K)+IDUMMY+MGL*(LAY(K)-1)
	 ENDDO   
         
       ENDIF   
       
       IF(T_ASSIM) THEN
         IDUMMY = IDUMMY + MGL*KBM1
       ENDIF
     ENDIF
         
     IF(TRIM(HEADINFO)=='!S') THEN
       IF(S_OBS) THEN
         READ(73,*) NUM
	 NLOC = NLOC + NUM
         IF(NLOC>RRK_NOBSMAX) THEN
           WRITE(IPT,*) 'not enough storage for observations:', 'Nloc=', Nloc, 'Nobsmax=', RRK_NOBSMAX
           CALL PSTOP
         ENDIF
	 READ(73,*)  (STLOC(K),K=NLOC-NUM+1,NLOC)        
         READ(73,*)  (LAY(K),  K=NLOC-NUM+1,NLOC)    
         DO K=NLOC-NUM+1, NLOC
	   STLOC(K)=STLOC(K)+IDUMMY+MGL*(LAY(K)-1)
	 ENDDO   
 
       ENDIF
       
       IF(S_ASSIM) THEN
         IDUMMY = IDUMMY + MGL*KBM1 
       ENDIF 
     ENDIF
     
     GOTO 100
 200 CONTINUE

     CLOSE(73)
    
  RETURN 
  END SUBROUTINE READOBS  

  SUBROUTINE READRESTART_STATE(NCF,STATE)
    USE MOD_INPUT
    USE MOD_STARTUP
    IMPLICIT NONE
    TYPE(NCFILE), POINTER :: NCF, NCF_tmP
    INTEGER :: STATE
    NCF_tmP=>NC_START
    IF(.not. associated(NCF)) CALL FATAL_ERROR&
         & ("MOD_RRK: READRESTART: THE FILE POINTER PASSED IS NOT ASSOCIATED?")  
    NCF%FTIME%PREV_STKCNT = STATE
    NC_START=>NCF
    
    

    CALL READ_SSH

    IF(WETTING_DRYING_ON) CALL READ_WETDRY

    CALL READ_UV

    CALL READ_TURB

    CALL READ_TS

    NC_START=>NCF_tmP

  END SUBROUTINE READRESTART_STATE

  SUBROUTINE READRESTART_TIME(NCF,NOW)
    USE MOD_INPUT
    USE MOD_STARTUP
    IMPLICIT NONE
    TYPE(NCFILE), POINTER :: NCF
    TYPE(TIME) :: NOW,TTEST
    INTEGER :: STATE,NSTATE
    TYPE(NCDIM), POINTER :: DIM
    LOGICAL FOUND
    
    IF(.not. associated(NCF)) CALL FATAL_ERROR&
         & ("MOD_RRK: READRESTART: THE FILE POINTER PASSED IS NOT ASSOCIATED?")

    ! GET THE TIME DIMENSION
   DIM => FIND_UNLIMITED(NCF,FOUND)
   IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN RESTART FILE",&
            & "FILE NAME: "//TRIM(NCF%FNAME),&
            &"COULD NOT FIND THE UNLIMITED DIMENSION")
   NSTATE = DIM%DIM
   FOUND = .false.
   DO STATE=1,NSTATE
      
      TTEST = get_file_time(NCF,STATE)
      IF(TTEST == NOW) THEN
         FOUND = .true.
         exit
      END IF
      
   END DO
   
   IF(FOUND)THEN
      CALL READRESTART_STATE(NCF,STATE)
   ELSE
      CALL PRINT_TIME(NOW,IPT,"TIME REQUSTED")
      CALL FATAL_ERROR&
           &("COULD NOT FIND CORRECT TIME IN RESTART FILE",&
           &"FILE NAME: "//TRIM(NCF%FNAME))
 
  END IF

  END SUBROUTINE READRESTART_TIME

  SUBROUTINE WRITERESTART(NCF,STATE)
    IMPLICIT NONE
    TYPE(NCFILE), POINTER :: NCF
    INTEGER :: STATE

    IF(.not. associated(NCF)) CALL FATAL_ERROR&
         & ("MOD_RRK: WRITERESTART: THE FILE POINTER PASSED IS NOT ASSOCIATED?")

    NCF%FTIME%NEXT_STKCNT = state
    CALL NC_WRITE_FILE(NCF,LOCAL_DISK)
    
  END SUBROUTINE WRITERESTART


!======================================================|
   FUNCTION SETUP_RESTART(FNAME) RESULT(NCF)
    USE MOD_NCDIO
    IMPLICIT NONE
    CHARACTER(LEN=*) :: FNAME    
    TYPE(GRID), SAVE :: MYGRID
    TYPE(NCFILE), POINTER :: NCF
    TYPE(NCFILE), POINTER :: NCF2
    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START SETUP_RSTFILE"

    IF(DBG_SET(DBG_LOG)) THEN
       
       write(IPT,*)"!--------------------------------------------------"
       write(IPT,*)"! SETTING UP RESTART FILE OUTPUT..."
   END IF
    NCF => NEW_FILE()
    NCF%FNAME=trim(fname)
    print *, 'running setup_restart:',NCF%FNAME ! should be mark
    NCF%FTIME=>NEW_FTIME()

    CALL SET_FVCOM_GRID(MYGRID)

    CALL DEFINE_DIMENSIONS(MYGRID)

    NCF2 => GRID_FILE_OBJECT(MYGRID)
    NCF => ADD(NCF,NCF2)
!!$    NCF => ADD(NCF,GRID_FILE_OBJECT(MYGRID) )

    NCF2 => TIME_FILE_OBJECT()
    NCF => ADD(NCF,NCF2)
!!$    NCF => ADD(NCF,TIME_FILE_OBJECT() )

    NCF2 => ZETA_FILE_OBJECT()
    NCF => ADD(NCF,NCF2)
!!$    NCF => ADD(NCF,ZETA_FILE_OBJECT() )

    NCF2 => FILE_DATE_OBJECT()
    NCF => ADD(NCF,NCF2)
!!$    NCF => ADD(NCF,FILE_DATE_OBJECT() )

    NCF2 => VELOCITY_FILE_OBJECT()
    NCF => ADD(NCF,NCF2)
!!$    NCF => ADD(NCF,VELOCITY_FILE_OBJECT() )

    NCF2 => AVERAGE_VEL_FILE_OBJECT()
    NCF => ADD(NCF,NCF2)
!!$    NCF => ADD(NCF,AVERAGE_VEL_FILE_OBJECT() )

    NCF2 => VERTICAL_VEL_FILE_OBJECT()
    NCF => ADD(NCF,NCF2)
!!$    NCF => ADD(NCF,VERTICAL_VEL_FILE_OBJECT() )

    NCF2 => TURBULENCE_FILE_OBJECT()
    NCF => ADD(NCF,NCF2)
!!$    NCF => ADD(NCF,TURBULENCE_FILE_OBJECT() )
    
    NCF2 => SALT_TEMP_FILE_OBJECT()
    NCF => ADD(NCF,NCF2)
!!$    NCF => ADD(NCF,SALT_TEMP_FILE_OBJECT() )

    NCF2 => DENSITY_FILE_OBJECT()
    NCF => ADD(NCF,NCF2)
!!$    NCF => ADD(NCF,DENSITY_FILE_OBJECT() )

    NCF2 => RESTART_EXTRAS_FILE_OBJECT()
    NCF => ADD(NCF,NCF2)
!!$    NCF => ADD(NCF,RESTART_EXTRAS_FILE_OBJECT() )
    
    IF(WETTING_DRYING_ON) THEN
       NCF2 => WET_DRY_FILE_OBJECT()
       NCF => ADD(NCF, NCF2)
!!$       NCF => ADD(NCF, WET_DRY_FILE_OBJECT() )
    END IF

#   if defined (BioGen)
    NCF2 => BIO_FILE_OBJECT()
    NCF => ADD(NCF,NCF2)
!!$    NCF => ADD(NCF,BIO_FILE_OBJECT() )
#   endif      

#   if defined (WATER_QUALITY)
    NCF2 => WQM_FILE_OBJECT()
    NCF => ADD(NCF,NCF2)
!!$    NCF => ADD(NCF,WQM_FILE_OBJECT() )
#   endif      

#   if defined (SEDIMENT)
    IF(SEDIMENT_MODEL)THEN
      NCF2 => SED_FILE_OBJECT()
      NCF => ADD(NCF,NCF2)
!!$      NCF => ADD(NCF,SED_FILE_OBJECT() )
    ENDIF
#   endif      

#   if defined (WAVE_CURRENT_INTERACTION)      
    NCF2 => WAVE_PARA_FILE_OBJECT()
    NCF => ADD(NCF,NCF2) 
!!$    NCF => ADD(NCF,WAVE_PARA_FILE_OBJECT() ) 

    NCF2 => WAVE_STRESS_FILE_OBJECT()
    NCF => ADD(NCF,NCF2) 
!!$    NCF => ADD(NCF,WAVE_STRESS_FILE_OBJECT() ) 
#   endif  


#   if defined (NH)
    NCF2 => NH_RST_FILE_OBJECT()
    NCF => ADD(NCF,NCF2)
!!$    NCF => ADD(NCF,NH_RST_FILE_OBJECT() )
#   endif

# if defined (ICE)
      !-----------------------------------------------------------------
      ! state variables
    NCF2 => ICE_RESTART_STATE_FILE_OBJECT()
    NCF => ADD(NCF,NCF2)
!!$    NCF => ADD(NCF,ICE_RESTART_STATE_FILE_OBJECT() )
      !-----------------------------------------------------------------
      ! velocity
      !-----------------------------------------------------------------
    NCF2 => ICE_VEL_FILE_OBJECT()
    NCF => ADD(NCF,NCF2)
!!$    NCF => ADD(NCF,ICE_VEL_FILE_OBJECT() )
      !-----------------------------------------------------------------
      ! fresh water, salt, and heat flux
      !-----------------------------------------------------------------
    NCF2 => ICE_EXTRA_FILE_OBJECT()
    NCF => ADD(NCF,NCF2)
!!$    NCF => ADD(NCF,ICE_EXTRA_FILE_OBJECT() )
# endif

    IF (STARTUP_TYPE /= "crashrestart") THEN
       CALL NC_WRITE_FILE(NCF)
       NCF%FTIME%NEXT_STKCNT = 1
    ELSE
       NCF%CONNECTED = .TRUE.
    END IF

    CALL KILL_DIMENSIONS

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END SETUP_RSTFILE"
  END FUNCTION SETUP_RESTART




# endif
END MODULE MOD_RRK

