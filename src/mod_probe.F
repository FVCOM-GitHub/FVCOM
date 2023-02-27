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

MODULE PROBES 
!==============================================================================!
!  PROBES MOD: CONTROLS TIME SERIES OUTPUT OF QUANTITIES AT SELECT LOCATIONS   !
!                                                                              !
!  CONTROL:  IF PROBE_ON = T IN RUNTIME PARAMETER CONTROL FILE TIME SERIES DUMP! 
!            IS ACTIVATED                                                      !
!                                                                              !
!  SETUP: TIME SERIES NML OBJECTS ARE PUT IN THE RUN FILE FOR EACH PROBE       !
!                                                                              !
!                                                                              !
!  SETUP FILE NAMING CONVENTION:                                               !
!     ***_timeseriesXX.dat                                                     !
!     *** = case id tag           (ex: gom/sat/mhb, etc)                       !
!     XX  = digit from 01 to 99, not sequential                                !
!                                                                              !
!  SETUP VARIABLE DESCRIPTION:                                                 !
!     D_LOC = element/node location of time-series data                        !
!     D_TYP = type of location ("element" or "node")                           !
!     D_TIT = title of location, used for time series output filename          !
!     D_DES = description of data (placed in time series output header)        !
!     O_INT = printing interval for data                                       !
!     K_ONE = initial sigma level 
!     K_TWO = final sigma level 
!     VAR   = variable to output in time series                                !
!     VNAME = variable description (placed in time series output header)       !
!                                                                              !
!     CURRENT VARS SET UP:                                                     !
!       u,v,w,ww,t1,s1,ua,va,rho1                                              !
!                                                                              !
!                                                                              !
!  NOTE 1:							               !
!    IF FILE ALREADY EXIST, A NEW FILE WILL BE CREATED NAMED FILENAME-01       !
!    IF FILENAME-01 EXISTS, A NEW FILE WILL BE CREATED NAMED FILENAME-02       !
!    ETC                                                                       !
!==============================================================================!

  USE MOD_UTILS
  USE MOD_PREC
  USE MOD_TIME
  USE CONTROL
  USE LIMS
  USE MOD_PAR
# if defined (SEDIMENT)

#  if defined (ORIG_SED)
  USE mod_sed
#  elif defined (CSTMS_SED)
  USE mod_sed_cstms
#  endif

# endif  
  IMPLICIT NONE
  SAVE
  !
  !--Probe Object Type               
  !
  TYPE PROBE_OBJ          
     LOGICAL   :: MINE
     TYPE(TIME):: O_INT
     TYPE(TIME):: O_NEXT
     INTEGER  :: D_LOC                                        !!LOCAL ELEMENT/NODE LOCATION   
     INTEGER  :: D_LOC_GL                                     !!GLOBAL ELEMENT/NODE LOCATION   
     INTEGER  :: K_ONE                                        !!INITIAL SIGMA LEVEL
     INTEGER  :: K_TWO                                        !!FINAL SIGMA LEVEL
     INTEGER  :: O_NUM                                        !!OBJECT NUMBER    
     REAL(SP) :: XLOC                                         !!X POSITION OF DATA LOCATION
     REAL(SP) :: YLOC                                         !!Y POSITION OF DATA LOCATION
     REAL(SP) :: LONLOC                                       !!LON POSITION OF DATA LOCATION
     REAL(SP) :: LATLOC                                       !!LAT POSITION OF DATA LOCATION
     REAL(SP) :: DPTH                                         !!BATHYMETRIC DEPTH AT DATA LOC
     REAL(SP), POINTER, DIMENSION(:)              :: VEC      !!USED FOR STORING DATA 
     REAL(SP), POINTER                            :: SCL      !!USED FOR STORING DATA 
     CHARACTER(LEN=80)                            :: D_TIT    !!OBJECT TITLE (FOR FILENAMING) 
     CHARACTER(LEN=80)                            :: D_DES    !!OBJECT DESCRIPTION
     CHARACTER(LEN=80)                            :: VAR      !!VARIABLE TO DUMP
     CHARACTER(LEN=80)                            :: VNAME    !!VARIABLE NAMES 
     CHARACTER(LEN=120)                           :: FILENAME !!FILE NAME      
  END TYPE PROBE_OBJ

   INTERFACE ASSIGNMENT(=)
      MODULE PROCEDURE ASSIGN_PROBE
   END INTERFACE

   INTERFACE MYPROBE
      MODULE PROCEDURE MYPROBE_VEC
      MODULE PROCEDURE MYPROBE_ARR
   END INTERFACE


  !
  !--Probe Namelist             
  !
  CHARACTER(LEN=80) :: PROBE_INTERVAL
  INTEGER           :: PROBE_LOCATION                                        !!LOCAL ELEMENT/NODE LOCATION
  INTEGER           :: PROBE_LEVELS(2)                                        !!FINAL SIGMA LEVEL
  CHARACTER(LEN=80) :: PROBE_TITLE    !!OBJECT TITLE (FOR FILENAMING) 
  CHARACTER(LEN=80) :: PROBE_DESCRIPTION    !!OBJECT DESCRIPTION
  CHARACTER(LEN=80) :: PROBE_VARIABLE      !!VARIABLE TO DUMP
  CHARACTER(LEN=80) :: PROBE_VAR_NAME  !!VARIABLE NAMES 

  NAMELIST /NML_PROBE/       &
       & PROBE_INTERVAL,     &
       & PROBE_LOCATION,     &
       & PROBE_LEVELS,       &
       & PROBE_TITLE,        &
       & PROBE_DESCRIPTION,  &
       & PROBE_VARIABLE,     &
       & PROBE_VAR_NAME

  !
  !--Probe Variables                    
  !
  TYPE(PROBE_OBJ), Allocatable  :: GLB_PROBE(:) !!GLOBAL PROBE OBJECTS
  TYPE(PROBE_OBJ), Allocatable  :: LCL_PROBE(:) !!LOCAL PROBE OBJECTS

  LOGICAL :: PROBE_ON
  INTEGER :: GLOBAL_PROBES, LOCAL_PROBES
  CHARACTER :: INPUT_FNAME

CONTAINS !------------------------------------------------------------------!
  ! SET_PROBES          :   READ PROBE OBJECTS FROM INPUT            !
  ! OPEN_PROBES         :   OPEN PROBE OUTPUT FILES                  !
  ! PROBE_STORE         :   ACCUMULATE DATA IN STORAGE ARRAY         !
  ! DUMP_PROBE_DATA     :   WRITE DATA TO PROBE FILES                !
  ! -----------------------------------------------------------------!

  !==============================================================================|
  !==============================================================================|

  SUBROUTINE INIT_NML_PROBE     
    IMPLICIT NONE

    PROBE_INTERVAL = 'none'
    PROBE_LOCATION   = -1                                     !!LOCAL ELEMENT/NODE LOCATION   
    PROBE_LEVELS     = -1                                      !!FINAL SIGMA LEVEL
    PROBE_TITLE      = 'none'  !!OBJECT TITLE (FOR FILENAMING) 
    PROBE_DESCRIPTION = 'none'   !!OBJECT DESCRIPTION
    PROBE_VARIABLE = 'none'     !!VARIABLE TO DUMP
    PROBE_VAR_NAME = 'none' !!VARIABLE NAMES 

  END SUBROUTINE INIT_NML_PROBE

SUBROUTINE ALLOC_PROBE(PROBE,N)     
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: N
    TYPE(PROBE_OBJ),Allocatable  :: PROBE(:)
    INTEGER :: STATUS,I

    ALLOCATE(PROBE(N),STAT=STATUS)
    IF(STATUS/=0) CALL FATAL_ERROR("MOD_PROBE: COULD NOT ALLOCATE PROBE TYPE")

    ! INITIALIZE TYPE DATA
    PROBE%MINE= .FALSE.
    PROBE%O_INT= ZEROTIME
    PROBE%O_NEXT= ZEROTIME
    PROBE%O_NUM=0
    PROBE%D_LOC=0
    PROBE%D_LOC_GL=0
    PROBE%K_ONE=0
    PROBE%K_TWO=0
    PROBE%xloc=0.0_sp
    PROBE%yloc=0.0_sp
    PROBE%lonloc=0.0_sp
    PROBE%latloc=0.0_sp
    PROBE%dpth=0.0_sp
    PROBE%D_TIT    ='none' 
    PROBE%D_DES    ='none' 
    PROBE%VAR      ='none' 
    PROBE%VNAME    ='none' 
    PROBE%FILENAME ='none' 
    DO i=1,N
       NULLIFY(PROBE(i)%SCL)
       NULLIFY(PROBE(i)%VEC)
    END DO

  END SUBROUTINE ALLOC_PROBE

!======================================================
   SUBROUTINE ASSIGN_PROBE(A,B)
     IMPLICIT NONE
     TYPE(PROBE_OBJ), INTENT(OUT) ::A
     TYPE(PROBE_OBJ), INTENT(IN)  ::B

     A%MINE     = B%MINE
     A%O_INT    = B%O_INT
     A%O_NEXT   = B%O_NEXT
     A%O_NUM    = B%O_NUM
     A%D_LOC    = B%D_LOC
     A%D_LOC_GL = B%D_LOC_GL
     A%K_ONE    = B%K_ONE
     A%K_TWO    = B%K_TWO
     A%xloc     = B%xloc
     A%yloc     = B%yloc
     A%lonloc   = B%lonloc
     A%latloc   = B%latloc
     A%dpth     = B%dpth
     A%D_TIT    = B%D_TIT
     A%D_DES    = B%D_DES
     A%VAR      = B%VAR 
     A%VNAME    = B%VNAME
     A%FILENAME = B%FILENAME

     A%VEC      =>B%VEC
     A%SCL      =>B%SCL

   END SUBROUTINE ASSIGN_PROBE

   SUBROUTINE MYPROBE_VEC(PROBE,VEC)
     USE ALL_VARS, only:XMC,YMC,H1,H,XM,YM,H,LON,LAT,LONC,LATC
     IMPLICIT NONE
     TYPE(PROBE_OBJ) :: PROBE
     REAL(SP), ALLOCATABLE, TARGET :: VEC(:)
     CHARACTER(LEN=80):: cstr1,cstr2,cstr3
     INTEGER :: I,J,IBND,PROCMAX
     !============SEE IF DATA POINT IS IN THE GLOBAL DOMAIN=================
     
     if(dbg_set(dbg_sbr)) write(IPT,*) "START: MYPROBE_VEC"
     
     IF (PROBE%K_ONE /= -1 .or. PROBE%K_TWO /= -1) CALL FATAL_ERROR&
          &("ERROR IN PROBE SETUP: PROBE_LEVELS should not be set for vector variables",&
          & "Do not specify it in the PROBE Namelist object for "//TRIM(PROBE_VARIABLE))

     IF(PROBE%D_LOC_GL  <1)THEN 
        write(cstr1,'(i8)') PROBE%D_LOC_GL
        write(cstr2,'(i8)') NGL
        
        CALL FATAL_ERROR('ERROR IN PROBE SETUP: DATA LOCATION'//TRIM(CSTR1)//' FOR TIME SERIES FILE: '//TRIM(PROBE_VARIABLE),&
             & 'IS NOT IN GLOBAL DOMAIN: 1 --> '//trim(cstr2))
     END IF
     
     IF( UBOUND(VEC,1) == NT)THEN

        IF(PROBE%D_LOC_GL > NGL)THEN 
           
           write(cstr1,'(i8)') PROBE%D_LOC_GL
           write(cstr2,'(i8)') NGL
           
           CALL FATAL_ERROR('ERROR IN PROBE SETUP: DATA LOCATION'//TRIM(CSTR1)//' FOR TIME SERIES FILE: '//TRIM(PROBE_VARIABLE),&
                & 'IS NOT IN GLOBAL DOMAIN: 1 --> '//trim(cstr2))
        END IF

        IF(ELID(PROBE%D_LOC_GL) /= 0) THEN
           
           PROBE%D_LOC = ELID(PROBE%D_LOC_GL)
           PROBE%MINE=.TRUE.
           
           ! METERS
           PROBE%XLOC = XMC(PROBE%D_LOC)
           PROBE%YLOC = YMC(PROBE%D_LOC)
           ! SPHERICAL
           PROBE%LONLOC = LONC(PROBE%D_LOC)
           PROBE%LATLOC = LATC(PROBE%D_LOC)
           ! DEPTH 
           PROBE%DPTH = H1(PROBE%D_LOC) 
           
           PROBE%SCL => VEC(PROBE%D_LOC)
        END IF
        
     ELSE IF ( UBOUND(VEC,1) == MT)THEN
        
        IF(PROBE%D_LOC_GL > MGL)THEN 
           write(cstr1,'(i8)') PROBE%D_LOC_GL
           write(cstr2,'(i8)') MGL
           
           CALL FATAL_ERROR('ERROR IN PROBE SETUP: DATA LOCATION'//TRIM(CSTR1)//' FOR TIME SERIES FILE: '//TRIM(PROBE_VARIABLE),&
                & 'IS NOT IN GLOBAL DOMAIN: 1 --> '//trim(cstr2))
           
        END IF

        IF(NLID(PROBE%D_LOC_GL) == 0) RETURN

        IF(NLID(PROBE%D_LOC_GL) > 0) THEN
           
           IF(NDE_ID(NLID(PROBE%D_LOC_GL)) == 1)THEN   !!BOUNDARY NODE
              DO J=1,NBN
                 IF(BN_LST(J) == PROBE%D_LOC_GL) IBND = J
              END DO
              !----Choose Processor of Lowest ID to be responsible for node
              PROCMAX = 10000
              DO J=1,NPROCS
                 IF(BN_NEY(IBND,J)==1) THEN
                    IF(J < PROCMAX) PROCMAX = J
                 END IF
              END DO
              
              IF(PROCMAX == MYID) THEN
                 PROBE%MINE = .TRUE. !!NOT RESPONSIBLE FOR TIME SERIES
                 PROBE%D_LOC = NLID(PROBE%D_LOC_GL)
                 
                 ! METERS
                 PROBE%XLOC = XM(PROBE%D_LOC)
                 PROBE%YLOC = YM(PROBE%D_LOC)
                 ! SPHERICAL
                 PROBE%LONLOC = LON(PROBE%D_LOC)
                 PROBE%LATLOC = LAT(PROBE%D_LOC)
                 ! DEPTH 
                 PROBE%DPTH = H(PROBE%D_LOC) 

                 PROBE%SCL => VEC(PROBE%D_LOC)
              END IF
            ELSE
              PROBE%MINE = .TRUE. !!NOT RESPONSIBLE FOR TIME SERIES
              PROBE%D_LOC = NLID(PROBE%D_LOC_GL)
                 
              ! METERS
              PROBE%XLOC = XM(PROBE%D_LOC)
              PROBE%YLOC = YM(PROBE%D_LOC)
              ! SPHERICAL
              PROBE%LONLOC = LON(PROBE%D_LOC)
              PROBE%LATLOC = LAT(PROBE%D_LOC)
              ! DEPTH 
              PROBE%DPTH = H(PROBE%D_LOC) 

              PROBE%SCL => VEC(PROBE%D_LOC)
   	    
           END IF
        END IF
        

     ELSE
        CALL FATAL_ERROR('MYPROBE: INVALID VARIABLE SIZE (Not equal MT or NT?)'&
             &,'Variable:'//TRIM(PROBE%VAR))
     END IF
     if(dbg_set(dbg_sbr)) write(IPT,*) "END: MYPROBE_VEC"

   END SUBROUTINE MYPROBE_VEC

   SUBROUTINE MYPROBE_ARR(PROBE,ARR)
     USE ALL_VARS, only:XMC,YMC,H1,H,XM,YM,H,LON,LAT,LONC,LATC
     IMPLICIT NONE
     TYPE(PROBE_OBJ) :: PROBE
     REAL(SP), ALLOCATABLE, TARGET :: ARR(:,:)
     CHARACTER(LEN=80):: cstr1,cstr2,cstr3
     INTEGER :: I,J,IBND,PROCMAX
     !============ENSURE SIGMA RANGE IS TENABLE=============================

     if(dbg_set(dbg_sbr)) write(IPT,*) "START: MYPROBE_ARR"

     IF(PROBE%K_ONE > UBOUND(ARR,2) .or. PROBE%K_TWO > UBOUND(ARR,2) )THEN
        CALL FATAL_ERROR('ERROR IN PROBE SETUP: PROBE LEVEL RANGE NOT CORRECT FOR VARIABLE: '//TRIM(PROBE_VARIABLE),&
             & 'MAKE SURE PROBE LEVELS ARE LESS THAN OR EQUAL TO THE NUMBER OF MODEL LEVELS') 
     END IF

     IF(PROBE%K_ONE < 1 .or. PROBE%K_TWO <1)THEN
        CALL FATAL_ERROR('ERROR IN PROBE SETUP: PROBE LEVEL RANGE NOT CORRECT FOR VARIABLE: '//TRIM(PROBE_VARIABLE),&
             & 'MAKE SURE PROBE LEVELS ARE GREATER THAN OR EQUAL TO ONE') 
     END IF
     
     
     IF(PROBE_LEVELS(1) >  PROBE_LEVELS(2) )THEN
        CALL FATAL_ERROR&
             &('ERROR IN PROBE SETUP: PROBE LEVEL RANGE NOT CORRECT FOR VARIABLE: '//TRIM(PROBE_VARIABLE),&
             & 'THE PROBE LEVEL INTERVAL MUST SPECIFY A VALID RANGE a:b') 
     END IF
     
     
     IF(PROBE%D_LOC_GL  <1)THEN 
        write(cstr1,'(i8)') PROBE%D_LOC_GL
        write(cstr2,'(i8)') NGL
        
        CALL FATAL_ERROR('ERROR IN PROBE SETUP: DATA LOCATION'//TRIM(CSTR1)//' FOR TIME SERIES FILE: '//TRIM(PROBE_VARIABLE),&
             & 'IS NOT IN GLOBAL DOMAIN: 1 --> '//trim(cstr2))
     END IF
     
     IF( UBOUND(ARR,1) == NT)THEN

        IF(PROBE%D_LOC_GL > NGL)THEN 
           
           write(cstr1,'(i8)') PROBE%D_LOC_GL
           write(cstr2,'(i8)') NGL
           
           CALL FATAL_ERROR('ERROR IN PROBE SETUP: DATA LOCATION'//TRIM(CSTR1)//' FOR TIME SERIES FILE: '//TRIM(PROBE_VARIABLE),&
                & 'IS NOT IN GLOBAL DOMAIN: 1 --> '//trim(cstr2))
        END IF

        IF(ELID(PROBE%D_LOC_GL) /= 0) THEN
           
           PROBE%D_LOC = ELID(PROBE%D_LOC_GL)
           PROBE%MINE=.TRUE.
           
           ! METERS
           PROBE%XLOC = XMC(PROBE%D_LOC)
           PROBE%YLOC = YMC(PROBE%D_LOC)
           ! SPHERICAL
           PROBE%LONLOC = LONC(PROBE%D_LOC)
           PROBE%LATLOC = LATC(PROBE%D_LOC)
           ! DEPTH 
           PROBE%DPTH = H1(PROBE%D_LOC) 
           
           PROBE%VEC => ARR(PROBE%D_LOC,PROBE%K_ONE:PROBE%K_TWO)
        END IF


     ELSE IF ( UBOUND(ARR,1) == MT)THEN

        IF(PROBE%D_LOC_GL > MGL)THEN 
           write(cstr1,'(i8)') PROBE%D_LOC_GL
           write(cstr2,'(i8)') MGL
           
           CALL FATAL_ERROR('ERROR IN PROBE SETUP: DATA LOCATION'//TRIM(CSTR1)//' FOR TIME SERIES FILE: '//TRIM(PROBE_VARIABLE),&
                & 'IS NOT IN GLOBAL DOMAIN: 1 --> '//trim(cstr2))
           
        END IF

        IF(NLID(PROBE%D_LOC_GL) == 0) RETURN

        IF(NLID(PROBE%D_LOC_GL) > 0) THEN
           
           IF(NDE_ID(NLID(PROBE%D_LOC_GL)) == 1)THEN   !!BOUNDARY NODE
              DO J=1,NBN
                 IF(BN_LST(J) == PROBE%D_LOC_GL) IBND = J
              END DO
              !----Choose Processor of Lowest ID to be responsible for node
              PROCMAX = 10000
              DO J=1,NPROCS
                 IF(BN_NEY(IBND,J)==1) THEN
                    IF(J < PROCMAX) PROCMAX = J
                 END IF
              END DO
              
              IF(PROCMAX == MYID) THEN
                 PROBE%MINE = .TRUE. !!NOT RESPONSIBLE FOR TIME SERIES
                 PROBE%D_LOC = NLID(PROBE%D_LOC_GL)

                 ! METERS
                 PROBE%XLOC = XM(PROBE%D_LOC)
                 PROBE%YLOC = YM(PROBE%D_LOC)
                 ! SPHERICAL
                 PROBE%LONLOC = LON(PROBE%D_LOC)
                 PROBE%LATLOC = LAT(PROBE%D_LOC)
                 ! DEPTH 
                 PROBE%DPTH = H(PROBE%D_LOC) 

                 PROBE%VEC => ARR(PROBE%D_LOC,PROBE%K_ONE:PROBE%K_TWO)
              END IF
           ELSE
             PROBE%MINE = .TRUE. !!NOT RESPONSIBLE FOR TIME SERIES
             PROBE%D_LOC = NLID(PROBE%D_LOC_GL)

             ! METERS
             PROBE%XLOC = XM(PROBE%D_LOC)
             PROBE%YLOC = YM(PROBE%D_LOC)
             ! SPHERICAL
             PROBE%LONLOC = LON(PROBE%D_LOC)
             PROBE%LATLOC = LAT(PROBE%D_LOC)
             ! DEPTH 
             PROBE%DPTH = H(PROBE%D_LOC) 

             PROBE%VEC => ARR(PROBE%D_LOC,PROBE%K_ONE:PROBE%K_TWO)
   	   
           END IF
        END IF
        

     ELSE
        CALL FATAL_ERROR('MYPROBE: INVALID VARIABLE SIZE?'&
             &,'Variable:'//TRIM(PROBE%VAR))
     END IF

     if(dbg_set(dbg_sbr)) write(IPT,*) "END: MYPROBE_ARR"

   END SUBROUTINE MYPROBE_ARR

 

  SUBROUTINE SET_PROBES(P_ON,NP,FNM)     
    !------------------------------------------------------------------------------|
    !  READ IN TIME SERIES OBJECTS FROM INPUT                                      |
    !------------------------------------------------------------------------------|
    USE MOD_PAR
    USE LIMS
    USE MOD_SET_TIME

    IMPLICIT NONE
    LOGICAL, INTENT(IN)  :: P_ON
    INTEGER, INTENT(IN)  :: NP
    CHARACTER(LEN=*), INTENT(IN):: FNM
    
    CHARACTER(LEN=80):: cstr1,cstr2,cstr3
    LOGICAL FEXIST,ISLOCAL
    INTEGER :: I,J,IERR,IOS,STATUS, N_PROBE
    INTEGER :: PROCMAX, IBND, charnum
    CHARACTER(LEN=120) :: pathnfile
    CHARACTER(LEN=4) :: OFLAG
    TYPE(TIME) :: OTIME
    INTEGER(ITIME) :: OSTEP
    if(DBG_SET(dbg_sbr)) &
         & write(IPT,*) "START: SET_PROBES;"

 
    PROBE_ON=P_ON
    GLOBAL_PROBES=NP
    LOCAL_PROBES=0
    INPUT_FNAME=FNM

    IF (.not. Probe_on) THEN
       if(DBG_SET(dbg_log)) &
            & write(IPT,*) "! Time Series Probes are off"
       return
    ELSE
       if(DBG_SET(dbg_log)) then
          write(IPT,*) "! Time Series Probes are on"
          write(IPT,*) "!    Setting up Probes:"
       end if
    END IF



    charnum = index (PROBES_FILE,".nml")
    if (charnum /= len_trim(PROBES_FILE)-3)&
         & CALL WARNING("PROBES FILE does not end in .nml", &
         & trim(PROBES_FILE))
    ! OPEN FILE - try both with appending input dir and without!
    pathnfile = trim(INPUT_DIR)//trim(PROBES_FILE)
    INQUIRE(FILE=PATHNFILE,EXIST=FEXIST)
    IF(FEXIST) THEN
       Call FOPEN(PROBEUNIT,trim(pathnfile),'cfr')
    ELSE
       pathnfile = trim(PROBES_FILE)
       Call FOPEN(PROBEUNIT,trim(pathnfile),'cfr')
       ! LET FOPEN CALL ERROR IF FILES DOES NOT EXIST
    END IF


    CALL ALLOC_PROBE(GLB_PROBE,GLOBAL_PROBES)

    !------------------------------------------------------------------------------|
    !  READ TIME SERIES NAME LIST                                                  |
    !------------------------------------------------------------------------------|
    N_PROBE = 0
    DO 

       CALL INIT_NML_PROBE

       ISLOCAL = .FALSE.

       READ(UNIT=PROBEUNIT, NML=NML_PROBE,IOSTAT=ios)
       if (IOS /= 0 ) exit

       N_PROBE = N_PROBE +1
       if  (N_PROBE > GLOBAL_PROBES) exit ! To prevent sigsev...


       !------------------------------------------------------------------------------|
       ! SET FUNDEMENTAL DATA
       !------------------------------------------------------------------------------|
       GLB_PROBE(N_PROBE)%K_ONE    = PROBE_LEVELS(1)
       GLB_PROBE(N_PROBE)%K_TWO    = PROBE_LEVELS(2)
       GLB_PROBE(N_PROBE)%D_LOC_GL = PROBE_LOCATION

       GLB_PROBE(N_PROBE)%D_TIT    =PROBE_TITLE 
       GLB_PROBE(N_PROBE)%D_DES    =PROBE_DESCRIPTION
       GLB_PROBE(N_PROBE)%VAR      =PROBE_VARIABLE
       GLB_PROBE(N_PROBE)%VNAME    =PROBE_VAR_NAME

       !------------------------------------------------------------------------------|
       ! SET TIME FOR PROBES
       !------------------------------------------------------------------------------|
       GLB_PROBE(N_PROBE)%O_NEXT = StartTime
       CALL IDEAL_TIME_STRING2TIME(PROBE_INTERVAL,Oflag,otime,ostep)
       IF (OFLAG == 'time') THEN ! IF OUTPUT TIME INTERVAL WAS SPECIFIED
          GLB_PROBE(N_PROBE)%O_INT = otime
          
       ELSE IF(OFLAG == 'step') THEN ! IF OUTPUT CYCLE INTERVAL WAS SPECIFIED          
          GLB_PROBE(N_PROBE)%O_INT = IMDTI * ostep
          
       END IF
       
       IF (GLB_PROBE(N_PROBE)%O_INT <= zerotime) CALL FATAL_ERROR&
            &('ERROR IN PROBE SETUP: Time series output interval is less than or equal to zero!')
       
       !------------------------------------------------------------------------------|
       !  POINT STORAGE DATA AND LOCATION TYPE                                        |
       !------------------------------------------------------------------------------|
       CALL PROBE_STORE(GLB_PROBE(N_PROBE))

       IF (GLB_PROBE(N_PROBE)%MINE) THEN
          LOCAL_PROBES=LOCAL_PROBES+1
          
          GLB_PROBE(N_PROBE)%O_NUM    =N_PROBE
          
          
       END IF

    END DO

    IF(GLOBAL_PROBES .ne. N_PROBE) THEN
       
       if(DBG_SET(dbg_log))  then
          write(ipt,*)"Bad NML_PROBE in the Name List!"
          write(ipt,*)"Specified number of PROBES=",GLOBAL_PROBES
          write(ipt,*)"But Found",N_PROBE, "; Valid PROBE name list objects.(Printing Last)"
          write(UNIT=IPT,NML=NML_PROBE)
       end if

       CALL FATAL_ERROR&
         &('THE NUMBER OF PROBES SPECIFIED IN THE RUN FILE CAN',&
         &'NOT BE FOUND IN THE PROBE FILE:'//trim(PROBES_FILE))
    END IF

    CLOSE(PROBEUNIT)

    
    CALL ALLOC_PROBE(LCL_PROBE,LOCAL_PROBES)

    N_PROBE=0
    DO I = 1,GLOBAL_PROBES
       IF (GLB_PROBE(I)%O_NUM .NE. 0) THEN
          N_PROBE=N_PROBE+1
          LCL_PROBE(N_PROBE)=GLB_PROBE(I)
       END IF
    END DO
    IF(N_PROBE .NE. LOCAL_PROBES) CALL FATAL_ERROR("mod_probe: this should not happen?")


    DEALLOCATE(GLB_PROBE)

    !------------------------------------------------------------------------------|
    !  PRINT STATISTICS ON TIME SERIES OBJECTS TO OUTPUT                           |
    !------------------------------------------------------------------------------|
    IF(GLOBAL_PROBES > 0)THEN

       if(dbg_set(dbg_sbrio)) write(ipt,*) "GLobal probes"&
            &,GLOBAL_PROBES,"Local_probes",Local_Probes,size(lcl_probe)

       IF(DBG_SET(DBG_LOG))THEN
          WRITE(IPT,*)
          WRITE(IPT,*)'!         TIME SERIES OBJECT DATA               '
          WRITE(IPT,*)" OBJ#  PROC GLOBAL  LOCAL  VAR  FILENAME"
       END  IF


       ! THIS IS VERY SLOW - DO NOT USE THIS METHOD INSIDE THE MAIN LOOP!
       DO J=1,NPROCS
#  if defined (MULTIPROCESSOR)
          IF(PAR) CALL MPI_BARRIER(MPI_FVCOM_GROUP,IERR)
#  endif

          IF(MYID == J) THEN
             DO I=1,LOCAL_PROBES  
                WRITE(IPT,101)LCL_PROBE(I)%O_NUM,MYID,LCL_PROBE(I)%D_LOC_GL,LCL_PROBE(I)%D_LOC, &
                     & TRIM(LCL_PROBE(I)%VAR),TRIM(LCL_PROBE(I)%D_TIT)
             END DO
          END IF
       END DO

    END IF

    !------------------------------------------------------------------------------|
    !  OPEN UP OUTPUT FILES AND WRITE HEADER INFORMATION                           |
    !------------------------------------------------------------------------------|
    CALL OPEN_PROBES

    CALL DUMP_PROBE_DATA

      if(dbg_set(dbg_sbr)) write(IPT,*) "END: SET_PROBES"

101 FORMAT(I5,I5,I8,I8,3X,A6,1X,A30)
  END SUBROUTINE SET_PROBES

    !==============================================================================|

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|

    !==============================================================================|

    SUBROUTINE OPEN_PROBES   

      !------------------------------------------------------------------------------|
      !  CREATE FILE NAMES AND WRITE HEADER INFORMATION FOR EACH TS OBJECT           |
      !------------------------------------------------------------------------------|

      USE MOD_PREC
      USE ALL_VARS
      IMPLICIT NONE
      INTEGER I,IUNIT,ICNT
      CHARACTER(LEN=120) :: FNAME
      CHARACTER(LEN=2  ) :: NAC 
      CHARACTER(LEN=3  ) :: APPEND
      LOGICAL FEXIST

      if(dbg_set(dbg_sbr)) write(IPT,*) "START: OPEN_PROBES"
      !
      !--Open Up Files -> If File Exists Create Secondary File (-01,-02, etc)
      !

      DO I=1,Local_Probes 

         ICNT  = 0
         FNAME = TRIM(OUTPUT_DIR)//TRIM(LCL_PROBE(I)%D_TIT)
         INQUIRE(FILE=FNAME,EXIST=FEXIST)
         DO WHILE(FEXIST)
            ICNT = ICNT + 1
            IF(ICNT .GE. 100) CALL FATAL_ERROR&
                 &("Please clean old time seris output in your results directory!")
            WRITE(NAC,'(I2.2)')ICNT
            APPEND = "-"//NAC
            FNAME = TRIM(OUTPUT_DIR)//TRIM(LCL_PROBE(I)%D_TIT)//TRIM(APPEND)
            INQUIRE(FILE=FNAME,EXIST=FEXIST)
         END DO

         IUNIT = LCL_PROBE(I)%O_NUM + 100
         CALL FOPEN(IUNIT,FNAME,'ofr')
         WRITE(IUNIT,*)TRIM(LCL_PROBE(I)%D_DES)
         WRITE(IUNIT,*)TRIM(LCL_PROBE(I)%VNAME)
         WRITE(IUNIT,*)
         CALL PRINT_REAL_TIME(GET_NOW(),IUNIT,"MODEL START DATE")
         WRITE(IUNIT,*)
         WRITE(IUNIT,*)'         K1            K2           '
         WRITE(IUNIT,'(2(I12,3X))')LCL_PROBE(I)%K_ONE,LCL_PROBE(I)%K_TWO
         WRITE(IUNIT,*)'     X(M)          Y(M)            DEPTH(M)'
         WRITE(IUNIT,'(3(F12.3,3X))')LCL_PROBE(I)%XLOC,LCL_PROBE(I)%YLOC,LCL_PROBE(I)%DPTH
         WRITE(IUNIT,*)'     LON           LAT               DEPTH(M)'
         WRITE(IUNIT,'(3(F12.3,3X))')LCL_PROBE(I)%LONLOC,LCL_PROBE(I)%LATLOC,LCL_PROBE(I)%DPTH
         WRITE(IUNIT,*)
         WRITE(IUNIT,*)'DATA FOLLOWS:'
         WRITE(IUNIT,*)'Time(days)    Data...'
         CLOSE(IUNIT)
         LCL_PROBE(I)%FILENAME = FNAME
      END DO

      if(dbg_set(dbg_sbr)) write(IPT,*) "END: OPEN_PROBES"

    END SUBROUTINE OPEN_PROBES
!==============================================================================|

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|

!==============================================================================|

    SUBROUTINE DUMP_PROBE_DATA

!------------------------------------------------------------------------------|
!  WRITE TIME SERIES DATA TO TIME SERIES DATA FILES                            |
!------------------------------------------------------------------------------|

      USE MOD_PREC
      USE ALL_VARS
      IMPLICIT NONE
      INTEGER I,K,K1,K2,IUNIT
      TYPE(TIME) :: OTIME

       if(dbg_set(dbg_sbr)) write(IPT,*) "START: DUMP_PROBE_DATA"

!==============================================================================!
!  MAIN LOOP OVER TIME SERIES OUTPUT                                           !
!==============================================================================!
      DO I=1,LOCAL_PROBES

         !----Return if not on Time Series Write Interval-----------------------------
         IF( LCL_PROBE(I)%O_NEXT .GT. IntTime) CYCLE 

         LCL_PROBE(I)%O_NEXT = IntTime + LCL_PROBE(I)%O_INT
         !----Open File For Write-----------------------------------------------------
         IUNIT = LCL_PROBE(I)%O_NUM + 100
         OPEN(UNIT=IUNIT,FILE=LCL_PROBE(I)%FILENAME,FORM='FORMATTED',POSITION='APPEND')

         !----Write Data to File------------------------------------------------------
         IF(ASSOCIATED(LCL_PROBE(I)%VEC))THEN
            K1=LCL_PROBE(I)%K_one
            K2=LCL_PROBE(I)%K_two
!            WRITE(IUNIT,*) DAYS(IntTime),(LCL_PROBE(I)%VEC(K),K=K1,K2)
            WRITE(IUNIT,'(f15.5,50f9.3)') DAYS(IntTime),(LCL_PROBE(I)%VEC)
         ELSE IF(ASSOCIATED(LCL_PROBE(I)%SCL))THEN
            WRITE(IUNIT,'(f15.5,50f9.3)') DAYS(IntTime),(LCL_PROBE(I)%SCL)
         END IF

      END DO

      if(dbg_set(dbg_sbr)) write(IPT,*) "END: DUMP_PROBE_DATA"
      RETURN
    END SUBROUTINE DUMP_PROBE_DATA

!==============================================================================!
!  END MAIN LOOP OVER TIME SERIES OUTPUT                                       !
!==============================================================================!

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|

!==============================================================================|

    SUBROUTINE PROBE_STORE(PROBE)

      !------------------------------------------------------------------------------|
      !  PUT TIME SERIES DATA IN TEMPORARY ARRAY                                     |
      !------------------------------------------------------------------------------|

      USE MOD_PREC
      USE ALL_VARS
# if defined (ICE)
      USE MOD_ICE
      USE MOD_ICE2D
# endif
      IMPLICIT NONE
      TYPE(PROBE_OBJ) :: PROBE

      if(dbg_set(dbg_sbr)) write(IPT,*) "START: PROBE_STORE"

      !--Store Data In Temporary Array-------------------------------------------
      SELECT CASE(TRIM(PROBE%VAR))

      CASE("u")
         CALL MYPROBE(PROBE,U)
      CASE("v")
         CALL MYPROBE(PROBE,V)
      CASE("w")
         CALL MYPROBE(PROBE,W)
      CASE("ww")
         CALL MYPROBE(PROBE,WW)
      CASE("q2")
         CALL MYPROBE(PROBE,Q2)
      CASE("q2l")
         CALL MYPROBE(PROBE,Q2L)
      CASE("l")
         CALL MYPROBE(PROBE,L)
      CASE("km")
         CALL MYPROBE(PROBE,KM)
      CASE("kq")
         CALL MYPROBE(PROBE,KQ)
      CASE("kh")
         CALL MYPROBE(PROBE,KH)
      CASE("t1")
         CALL MYPROBE(PROBE,T1)
      CASE("s1")
         CALL MYPROBE(PROBE,S1)
      CASE("rho1")
         CALL MYPROBE(PROBE,RHO1)
      CASE("ua")
         CALL MYPROBE(PROBE,UA)
      CASE("va")
         CALL MYPROBE(PROBE,VA)
      CASE("el")
         CALL MYPROBE(PROBE,EL)
# if defined (ICE)
      CASE("aice")
         CALL MYPROBE(PROBE,AICE)
      CASE("vice")
         CALL MYPROBE(PROBE,VICE)
      CASE("uice2")
         CALL MYPROBE(PROBE,UICE2)
      CASE("vice2")
         CALL MYPROBE(PROBE,VICE2)
# endif
# if defined (ORIG_SED)
      CASE("csed")
         CALL MYPROBE(PROBE,CSED)
# endif
      CASE DEFAULT
         CALL FATAL_ERROR&
              &('VARIABLE: '//TRIM(PROBE%VAR)//' HAS NOT BEEN SET UP',&
              & 'FOR TIME SERIES OUTPUT (Did you use CAPITALS by mistake?)',&
              & 'MODIFY MOD_PROBE TO ADD IT!')

      END SELECT

      if(dbg_set(dbg_sbr)) write(IPT,*) "END: PROBE_STORE"

    END SUBROUTINE PROBE_STORE
!==============================================================================|

  END MODULE PROBES 
!==============================================================================|
