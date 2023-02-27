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

MODULE ENKFVAL
# if defined (ENKF)
   USE mod_ncdio, only : update_iodata   
   USE CONTROL
   IMPLICIT NONE
   CHARACTER(LEN=2)   ENKF_NUM
   CHARACTER(LEN=80)  ENKF_INIT   !!INITIAL PERTUBATION FIELD OPTION(usr/default)
   CHARACTER(LEN=80)  ENKF_TYPE   !!ENSEMBLE UPDATE OPTION(default/square root)

   INTEGER      ENKF_NENS         !!GLOBAL NUMBER OF ENSEMBLES
   INTEGER      DELTA_ASS         !!ASSIMILATON TIME INTERVAL IN SECONDS
!   INTEGER      ENKF_INT          !!ASSIMILATION TIME INTERVAL/FILE OUTPUT INTERVAL (>=1) 
   INTEGER      ENKF_NOBSMAX      !!MAXIMUM NUMBER OF THE OBSERVATION STATIONS  
   INTEGER      ENKF_START        !!ASSIMILATION START TIME      
   INTEGER      ENKF_END          !!ASSIMILATION END TIME
   REAL(SP) ::  ENKF_CINF         !!MAX DISTANCE OF CORRELATIN  
   REAL(SP) ::  OBSERR_EL         !!EL OBSERVATION ERROR SPECIFIED   
   REAL(SP) ::  OBSERR_UV         !!UV OBSERVATION ERROR SPECIFIED  
   REAL(SP) ::  OBSERR_T          !!TEMPERATURE OBSERVATION ERROR SPECIFIED
   REAL(SP) ::  OBSERR_S          !!SALINITY OBSERVATION ERROR SPECIFIED

   LOGICAL  ::  EL_ASSIM          !!OPTION FOR CHOSING ELEVATION AS ASSIMILATION VARIABLES
   LOGICAL  ::  UV_ASSIM          !!OPTION FOR CHOSING CURRENT AS ASSIMILATION VARIABLES  
   LOGICAL  ::  T_ASSIM           !!OPTION FOR CHOSING TEMPERATURE AS ASSIMILATION VARIABLES
   LOGICAL  ::  S_ASSIM           !!OPTION FOR CHOSING SALINITY AS ASSIMILATION VARIABLES
   LOGICAL :: ENKF_LOCALIZED      !! only do enkf for certain region
   LOGICAL  ::  EL_OBS            !!OPTION FOR ELEVATION OBSERVATION DATA
   LOGICAL  ::  UV_OBS            !!OPTION FOR CURRENT OBSERVATION DATA  
   LOGICAL  ::  T_OBS             !!OPTION FOR TEMPERATURE OBSERVATION DATA
   LOGICAL  ::  S_OBS             !!OPTION FOR SALINITY OBSERVATION DATA
  
   INTEGER  ::  ENKF_METHOD         !! DIFFERENENT ENKF_METHODS
   INTEGER  ::  MODE              !! EVENSEN_KF_OPTIONS
   INTEGER  ::  INOOB             !! 72 RESERVE FOR I/O INPUT OF OBSERVATION FILE
   INTEGER  ::  INOKF             !! 73 FILE I/O PIPE NUMBER
   INTEGER  ::  IOBCKF            !! 76 I/O PIPE ONLY FOR B.C. TREATMENT IN ENKF 

   INTEGER  ::  NLOC=0              !!Number of observation locations
   INTEGER  ::  IENS              !!ensemble number index
   INTEGER  ::  NCYC              !!Number of cycles for assimilation
   INTEGER  ::  ICYC = 0          !!cycle number index
   CHARACTER(LEN=4)   FCYC 

   REAL(DP),PARAMETER :: ZEROD = 0.0_DP
  
   INTEGER :: KF_TIMES =0 
   INTEGER :: ARCHIVE_TIMES=0
   TYPE(TIME) :: enkf_out_Time
# endif
END MODULE ENKFVAL

!=============================================
!MODULE MOD_STARTUP_MIMIC
!# if defined (ENKF)
!  USE MOD_UTILS
!  USE MOD_NCTOOLS
!  USE MOD_INPUT
!  USE ALL_VARS
!  USE EQS_OF_STATE
!  USE MOD_WD
!  USE SINTER
!
!
!#  if defined (ICE)
!  USE MOD_ICE
!  USE MOD_ICE2D
!#  endif
!
!# if defined (WATER_QUALITY)
!  USE MOD_WQM
!# endif
!  
!# if defined (BioGen)
!  USE MOD_BIO_3D
!# endif
!
!# if defined (NH)
!  USE NON_HYDRO, ONLY: W4ZT, NHQDRX, NHQDRY, NHQDRZ, NHQ2DX, NHQ2DY
!# endif
!  
!  IMPLICIT NONE
!  
!  PRIVATE
!
!  PUBLIC :: READ_SSH1
!  PUBLIC :: READ_UV1
!  PUBLIC :: READ_TURB1
!  PUBLIC :: READ_TS1
!  PUBLIC :: READ_WETDRY1
!
!  
!CONTAINS
!!==============================================================================!
!  SUBROUTINE READ_SSH1(NC_START)
!#   if defined (SEDIMENT)
!    USE MOD_SED, only : morpho_model,sed_hot_start 
!#   endif
!    IMPLICIT NONE
!    TYPE(NCFILE),POINTER :: NC_START
!    TYPE(NCVAR),  POINTER :: VAR
!    TYPE(NCDIM),  POINTER :: DIM
!    LOGICAL :: FOUND
!    INTEGER :: STKCNT
!
!    IF(DBG_SET(DBG_SBR)) write(ipt,*) "Start: READ_SSH"
!    
!    STKCNT = NC_START%FTIME%PREV_STKCNT
!
!    ! LOAD EL
!    VAR => FIND_VAR(NC_START,'zeta',FOUND)
!    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'zeta'&
!         & IN THE HOTSTART FILE OBJECT")
!    CALL NC_CONNECT_AVAR(VAR, EL)
!    CALL NC_READ_VAR(VAR,STKCNT)
!
!    ! LOAD ET
!    VAR => FIND_VAR(NC_START,'et',FOUND)
!    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'et'&
!         & IN THE HOTSTART FILE OBJECT")
!    CALL NC_CONNECT_AVAR(VAR, ET)
!    CALL NC_READ_VAR(VAR,STKCNT)
!
!    !----------------------------------------------------------------
!    ! Read the most recent bathymetry if Morphodynamics is Active
!    !----------------------------------------------------------------
!#   if defined(SEDIMENT)
!    IF(MORPHO_MODEL .and. SED_HOT_START)THEN
!    VAR => FIND_VAR(NC_START,'h',FOUND)
!    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'h'&
!         & IN THE HOTSTART FILE OBJECT")
!    CALL NC_CONNECT_AVAR(VAR, h)
!    CALL NC_READ_VAR(VAR,STKCNT)
!    ENDIF
!#   endif
!
!
!    !----------------------------------------------------------------
!    ! Given SSH and Bathy, Update the Bathymetry 
!    !----------------------------------------------------------------
!    D  = H + EL
!    DT = H + ET
!
!    
!    CALL N2E2D(H,H1)
!    CALL N2E2D(EL,EL1)
!    CALL N2E2D(D,D1)
!    CALL N2E2D(DT,DT1)
!
!    IF(DBG_SET(DBG_SBR)) write(ipt,*) "Start: READ_SSH"
!
!  END SUBROUTINE READ_SSH1
!!==============================================================================!
!!==============================================================================!
!  SUBROUTINE READ_WETDRY1(NC_START)
!    USE MOD_WD
!    IMPLICIT NONE
!    TYPE(NCFILE),POINTER :: NC_START
!    TYPE(NCVAR),  POINTER :: VAR
!    TYPE(NCDIM),  POINTER :: DIM
!    LOGICAL :: FOUND
!    INTEGER :: STKCNT
!
!    IF(DBG_SET(DBG_SBR)) write(ipt,*) "Start: READ_WETDRY"
!
!    STKCNT = NC_START%FTIME%PREV_STKCNT
!
!    ! LOAD ISWETN
!    VAR => FIND_VAR(NC_START,'wet_nodes',FOUND)
!    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'wet_nodes'&
!         & IN THE HOTSTART FILE OBJECT")
!    CALL NC_CONNECT_AVAR(VAR, ISWETN)
!    CALL NC_READ_VAR(VAR,STKCNT)
!
!    ! LOAD ISWETC
!    VAR => FIND_VAR(NC_START,'wet_cells',FOUND)
!    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'wet_cells'&
!         & IN THE HOTSTART FILE OBJECT")
!    CALL NC_CONNECT_AVAR(VAR, ISWETC)
!    CALL NC_READ_VAR(VAR,STKCNT)
!
!
!    ! LOAD ISWETNT
!    VAR => FIND_VAR(NC_START,'wet_nodes_prev_int',FOUND)
!    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'wet_nodes_prev_int'&
!         & IN THE HOTSTART FILE OBJECT")
!    CALL NC_CONNECT_AVAR(VAR, ISWETNT)
!    CALL NC_READ_VAR(VAR,STKCNT)
!
!    ! LOAD ISWETCT
!    VAR => FIND_VAR(NC_START,'wet_cells_prev_int',FOUND)
!    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'wet_cells_prev_int'&
!         & IN THE HOTSTART FILE OBJECT")
!    CALL NC_CONNECT_AVAR(VAR, ISWETCT)
!    CALL NC_READ_VAR(VAR,STKCNT)
!
!    ! LOAD ISWETCE
!    VAR => FIND_VAR(NC_START,'wet_cells_prev_ext',FOUND)
!    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'wet_cells_prev_ext'&
!         & IN THE HOTSTART FILE OBJECT")
!    CALL NC_CONNECT_AVAR(VAR, ISWETC)
!    CALL NC_READ_VAR(VAR,STKCNT)
!
!    IF(DBG_SET(DBG_SBR)) write(ipt,*) "End: READ_WETDRY"
!    
!  END SUBROUTINE READ_WETDRY1
!!==============================================================================!
!  SUBROUTINE READ_TS1(NC_START)
!    IMPLICIT NONE
!    TYPE(NCFILE),POINTER :: NC_START
!    TYPE(NCVAR),  POINTER :: VAR
!    TYPE(NCDIM),  POINTER :: DIM
!    LOGICAL :: FOUND
!    INTEGER :: STKCNT, K
!    REAL(SP), DIMENSION(0:MT,KB) :: PRESSURE
!    
!    IF(DBG_SET(DBG_SBR)) write(ipt,*) "Start: READ_TS" 
!   
!    STKCNT = NC_START%FTIME%PREV_STKCNT
!
!
!    ! LOAD TEMPERATURE
!    VAR => FIND_VAR(NC_START,'temp',FOUND)
!    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'temp'&
!         & IN THE HOTSTART FILE OBJECT")
!    CALL NC_CONNECT_AVAR(VAR, T1)
!    CALL NC_READ_VAR(VAR,STKCNT)
!
!
!    ! LOAD MEAN INITIAL TEMPERATURE
!    VAR => FIND_VAR(NC_START,'tmean1',FOUND)
!    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'tmean1'&
!         & IN THE HOTSTART FILE OBJECT")
!    CALL NC_CONNECT_AVAR(VAR, tmean1)
!    CALL NC_READ_VAR(VAR)
!
!    ! LOAD SALINITY
!    VAR => FIND_VAR(NC_START,'salinity',FOUND)
!    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'saltinity'&
!         & IN THE HOTSTART FILE OBJECT")
!    CALL NC_CONNECT_AVAR(VAR, S1)
!    CALL NC_READ_VAR(VAR,STKCNT)
!
!    ! LOAD MEAN INITIAL SALINITY
!    VAR => FIND_VAR(NC_START,'smean1',FOUND)
!    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'smean1'&
!         & IN THE HOTSTART FILE OBJECT")
!    CALL NC_CONNECT_AVAR(VAR, smean1)
!    CALL NC_READ_VAR(VAR)
!
!
!    ! LOAD DENSITY
!    VAR => FIND_VAR(NC_START,'rho1',FOUND)
!    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'rho1'&
!         & IN THE HOTSTART FILE OBJECT")
!    CALL NC_CONNECT_AVAR(VAR, RHO1)
!    CALL NC_READ_VAR(VAR,STKCNT)
!
!    ! LOAD MEAN DENSITY
!    VAR => FIND_VAR(NC_START,'rmean1',FOUND)
!    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'rmean1'&
!         & IN THE HOTSTART FILE OBJECT")
!    CALL NC_CONNECT_AVAR(VAR, rmean1)
!    CALL NC_READ_VAR(VAR)
!
!    ! AVERAGE FROM CELLS TO FACE CENTERS
!
!
!!JQI    SELECT CASE(SEA_WATER_DENSITY_FUNCTION)
!!JQI    CASE(SW_DENS1)
!
!!JQI       ! SET MEAN DENSITY
!!JQI       DO K=1,KBM1
!!JQI          PRESSURE(:,K) = -GRAV_N*1.025_SP*(ZZ(:,K)*D(:))*0.1_SP
!!JQI       END DO
!!JQI       CALL FOFONOFF_MILLARD(SMEAN1,TMEAN1,PRESSURE,0.0_SP,RMEAN1)
!!JQI       RMEAN1(0,:)=0.0_SP
!!JQI       RMEAN1(:,KB)=0.0_SP
!
!!JQI       ! SET REAL DENSITY
!!JQI       CALL DENS1 ! GENERIC CALL TO FOFONOFF_MILLARD FOR S1,T1...
!       
!!JQI    CASE(SW_DENS2)
!!JQI       ! SET MEAN DENSITY
!!JQI       CALL DENS2G(SMEAN1,TMEAN1,RMEAN1)
!!JQI       RMEAN1(0,:)=0.0_SP
!!JQI       RMEAN1(:,KB)=0.0_SP
!
!!JQI       ! SET REAL DENSITY
!!JQI       CALL DENS2 ! GENERIC CALL TO DENS2G FOR S1,T1...
!
!!JQI    CASE(SW_DENS3)
!
!!JQI       ! SET MEAN DENSITY
!!JQI       DO K=1,KBM1
!!JQI          PRESSURE(:,K) = -GRAV_N*1.025_SP*(ZZ(:,K)*D(:))*0.1_SP
!!JQI       END DO
!!JQI       CALL JACKET_MCDOUGALL(SMEAN1,TMEAN1,PRESSURE,RMEAN1)
!!JQI       RMEAN1(0,:)=0.0_SP
!!JQI       RMEAN1(:,KB)=0.0_SP
!
!!JQI       ! SET REAL DENSITY
!!JQI       CALL DENS3 ! GENERIC CALL TO JACKET_MCDOUGALL FOR S1,T1.
!
!!JQI    CASE DEFAULT
!!JQI       CALL FATAL_ERROR("INVALID DENSITY FUNCTION SELECTED:",&
!!JQI            & "   "//TRIM(SEA_WATER_DENSITY_FUNCTION) )
!!JQI    END SELECT
!    
!    CALL N2E3D(T1,T)
!    CALL N2E3D(S1,S)
!    CALL N2E3D(Tmean1,Tmean)
!    CALL N2E3D(Smean1,Smean)
!    CALL N2E3D(Rmean1,Rmean)
!
!   
!   IF(DBG_SET(DBG_SBR)) write(ipt,*) "End: READ_TS" 
!   
! END SUBROUTINE READ_TS1
!!==============================================================================!
!
!  SUBROUTINE READ_UV1(NC_START)
!    IMPLICIT NONE
!    TYPE(NCFILE),POINTER :: NC_START
!    TYPE(NCVAR),  POINTER :: VAR
!    TYPE(NCDIM),  POINTER :: DIM
!    LOGICAL :: FOUND
!    INTEGER :: STKCNT
!
!    IF(DBG_SET(DBG_SBR)) write(ipt,*) "START: READ_UV" 
!
!
!    STKCNT = NC_START%FTIME%PREV_STKCNT
!
!    VAR => FIND_VAR(NC_START,'u',FOUND)
!    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'u'&
!         & IN THE HOTSTART FILE OBJECT")
!    CALL NC_CONNECT_AVAR(VAR, U)
!    CALL NC_READ_VAR(VAR,STKCNT)
!
!
!    VAR => FIND_VAR(NC_START,'v',FOUND)
!    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'v'&
!         & IN THE HOTSTART FILE OBJECT")
!    CALL NC_CONNECT_AVAR(VAR, V)
!    CALL NC_READ_VAR(VAR,STKCNT)
!
!    VAR => FIND_VAR(NC_START,'omega',FOUND)
!    IF(FOUND) THEN
!       CALL NC_CONNECT_AVAR(VAR, WTS)
!       CALL NC_READ_VAR(VAR,STKCNT)
!       
!       CALL N2E3D(WTS,W)
!    ELSE
!       VAR => FIND_VAR(NC_START,'w',FOUND)
!       IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'w' &
!            & or 'omega' IN THE HOTSTART FILE OBJECT")
!       CALL NC_CONNECT_AVAR(VAR, W)
!       CALL NC_READ_VAR(VAR,STKCNT)
!    END IF
!
!    VAR => FIND_VAR(NC_START,'ua',FOUND)
!    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'ua'&
!         & IN THE HOTSTART FILE OBJECT")
!    CALL NC_CONNECT_AVAR(VAR, UA)
!    CALL NC_READ_VAR(VAR,STKCNT)
!
!    VAR => FIND_VAR(NC_START,'va',FOUND)
!    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'va'&
!         & IN THE HOTSTART FILE OBJECT")
!    CALL NC_CONNECT_AVAR(VAR, VA)
!    CALL NC_READ_VAR(VAR,STKCNT)
!
!    IF(DBG_SET(DBG_SBR)) write(ipt,*) "END: READ_UV" 
!
!  END SUBROUTINE READ_UV1
!!==============================================================================!
!  SUBROUTINE READ_TURB1(NC_START)
!    IMPLICIT NONE
!    TYPE(NCFILE),POINTER :: NC_START
!    TYPE(NCVAR),  POINTER :: VAR
!    TYPE(NCDIM),  POINTER :: DIM
!    LOGICAL :: FOUND
!    INTEGER :: STKCNT
!
!    IF(DBG_SET(DBG_SBR)) write(ipt,*) "START: READ_TURB" 
!
!    STKCNT = NC_START%FTIME%PREV_STKCNT
!
!#    if defined (GOTM)
!
!    VAR => FIND_VAR(NC_START,'tke',FOUND)
!    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'tke'&
!         & IN THE HOTSTART FILE OBJECT")
!    CALL NC_CONNECT_AVAR(VAR, TKE)
!    CALL NC_READ_VAR(VAR,STKCNT)
!
!    VAR => FIND_VAR(NC_START,'teps',FOUND)
!    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'teps'&
!         & IN THE HOTSTART FILE OBJECT")
!    CALL NC_CONNECT_AVAR(VAR, TEPS)
!    CALL NC_READ_VAR(VAR,STKCNT)
!
!    L = .001 
!    L(1:MT,2:KBM1) = (.5544**3)*TKE(1:MT,2:KBM1)**1.5/TEPS(1:MT,2:KBM1)
!    
!# else
!
!    VAR => FIND_VAR(NC_START,'q2',FOUND)
!    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'q2'&
!         & IN THE HOTSTART FILE OBJECT")
!    CALL NC_CONNECT_AVAR(VAR, Q2)
!    CALL NC_READ_VAR(VAR,STKCNT)
!
!    VAR => FIND_VAR(NC_START,'q2l',FOUND)
!    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'q2l'&
!         & IN THE HOTSTART FILE OBJECT")
!    CALL NC_CONNECT_AVAR(VAR, Q2L)
!    CALL NC_READ_VAR(VAR,STKCNT)
!
!    VAR => FIND_VAR(NC_START,'l',FOUND)
!    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'l'&
!         & IN THE HOTSTART FILE OBJECT")
!    CALL NC_CONNECT_AVAR(VAR, L)
!    CALL NC_READ_VAR(VAR,STKCNT)
!
!# endif
!
!    VAR => FIND_VAR(NC_START,'km',FOUND)
!    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'km'&
!         & IN THE HOTSTART FILE OBJECT")
!    CALL NC_CONNECT_AVAR(VAR, km)
!    CALL NC_READ_VAR(VAR,STKCNT)
!
!    VAR => FIND_VAR(NC_START,'kq',FOUND)
!    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'kq'&
!         & IN THE HOTSTART FILE OBJECT")
!    CALL NC_CONNECT_AVAR(VAR, KQ)
!    CALL NC_READ_VAR(VAR,STKCNT)
!
!    VAR => FIND_VAR(NC_START,'kh',FOUND)
!    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'kh'&
!         & IN THE HOTSTART FILE OBJECT")
!    CALL NC_CONNECT_AVAR(VAR, KH)
!    CALL NC_READ_VAR(VAR,STKCNT)
!
!
!    CALL N2E3D(KM,KM1)
!
!    IF(DBG_SET(DBG_SBR)) write(ipt,*) "END: READ_TURB" 
!
!  END SUBROUTINE READ_TURB1
!!==============================================================================|
!# endif
!END MODULE MOD_STARTUP_MIMIC

MODULE MOD_ENKF
   USE ENKFVAL
   USE mod_ncdio, only : update_iodata 
   USE MOD_INPUT, only : NC_START
   USE CONTROL
   USE MOD_UTILS
   USE MOD_NCTOOLS
   IMPLICIT NONE
   SAVE
   

   LOGICAL  :: ENKF_ON

# if defined (ENKF)
   INTEGER  ::  ENKF_memberCONTR !, ICYC
   CHARACTER(LEN=80) :: ENKF_START_DATE
   CHARACTER(LEN=80) :: ENKF_END_DATE

   TYPE(TIME) :: ENKF_START_TIME
   TYPE(TIME) :: ENKF_END_Time

   TYPE(TIME) :: ENKF_cyc

   CHARACTER(LEN=80) :: ENKF_ASSIM_INTERVAL
   TYPE(TIME) :: ENKF_INTERVAL



!   TYPE(NCFILE), POINTER :: NCF_RRKFDATA

   TYPE(NCFILE), POINTER :: NCF_ENKFRE
   TYPE(NCFILE), POINTER :: NCF_ENKFfct
   TYPE(NCFILE), POINTER :: enkf_ncfout,enkf_mean_ncfout
   LOGICAL :: LOCAL_DISK
   CHARACTER(LEN=80), parameter :: Local_file ="/scratch/enkf_restart.nc"
   CHARACTER(LEN=80), parameter :: group_file ="rrktemp/enkf_restart.nc"

   CHARACTER(LEN=80), parameter :: Local_file1 ="/scratch/enkf_forecast.nc"
   CHARACTER(LEN=80), parameter :: group_file1 ="rrktemp/enkf_forecast.nc"

   INTEGER(ITIME) :: EKINT_COUNT, EKINT_START, EKINT_END
!   TYPE(NCFILE), POINTER :: enkf_ncfout_enkf(20)

!   INTERFACE READRESTART
!      MODULE PROCEDURE READRESTART_TIME
!      MODULE PROCEDURE READRESTART_STATE
!   END INTERFACE

!   PRIVATE READRESTART_TIME
!   PRIVATE READRESTART_STATE


   NAMELIST /NML_ENKF/        &
         & ENKF_ON,            &
         & ENKF_START_DATE,    &
         & ENKF_END_DATE,      &
         & ENKF_ASSIM_INTERVAL,&
         & ENKF_NOBSMAX,       &
         & ENKF_NENS,          &
         & ENKF_CINF,          &
         & EKINT_START,       &
         & EL_ASSIM,          &
         & EL_OBS,            &
         & UV_ASSIM,          &
         & UV_OBS,            &
         & T_ASSIM,           &
         & T_OBS,             &
         & S_ASSIM,           &
         & S_OBS,             &
         & ENKF_LOCALIZED,    &
         &  ENKF_METHOD,      &
         &  MODE,             &
         &  OBSERR_EL,        & 
         &  OBSERR_UV,        &
         &  OBSERR_T,         &
         &  OBSERR_S,         &
         & LOCAL_DISK          

!     WRITE(IPT,*) '!  # SPECIFICY INITIAL PERTUBATION FIELD             :',ENKF_INIT
!     WRITE(IPT,*) '!  # GLOBAL NUMBER OF ENSEMBLES                      :',ENKF_NENS
!     WRITE(IPT,*) '!  # ASSIMILATON TIME INTERVAL IN SECONDS            :',DELTA_ASS
!     WRITE(IPT,*) '!  # ASSIMILATION TIME INTERVAL/FILE OUTPUT INTERVAL :',ENKF_INT
!     WRITE(IPT,*) '!  # MAXIMUM NUMBER OF THE OBSERVATION STATIONS      :',ENKF_NOBSMAX
!     WRITE(IPT,*) '!  # ASSIMILATION START TIME                         :',ENKF_START
!     WRITE(IPT,*) '!  # ASSIMILATION END TIME                           :',ENKF_END
!     WRITE(IPT,*) '!  # TIDAL AMPLITUDE ERROR RANGE SPECIFIED           :',(BC_AMP_ERR(I),I=1,6)
!     WRITE(IPT,*) '!  # TIDAL PHASE ERROR RANGE SPECIFIED               :',(BC_PHA_ERR(I),I=1,6)
!     WRITE(IPT,*) '!  # MAX DISTANCE OF CORRELATIN                      :',ENKF_CINF
!     WRITE(IPT,*) '!  # ELEVATION AS ASSIMILATION VARIABLES             :',EL_ASSIM
!     WRITE(IPT,*) '!  # ELEVATION AS OBSERVATION DATA OPTION            :',EL_OBS
!     WRITE(IPT,*) '!  # CURRENTS AS ASSIMILATION VARIABLES              :',UV_ASSIM
!     WRITE(IPT,*) '!  # CURRENTS OBSERVATION DATA OPTION                :',UV_OBS
!     WRITE(IPT,*) '!  # TEMPERATURE AS ASSIMILATION VARIABLES           :',T_ASSIM
!     WRITE(IPT,*) '!  # TEMPERATURE OBSERVATION DATA OPTION             :',T_OBS
!     WRITE(IPT,*) '!  # SALINITY AS ASSIMILATION VARIABLES              :',S_ASSIM
!     WRITE(IPT,*) '!  # SALINITY OBSERVATION DATA OPTION                :',S_OBS
!     IF(EL_OBS) WRITE(IPT,*) '!  # EL OBSERVATION ERROR SPECIFIED                  :',OBSERR_EL
!     IF(UV_OBS) WRITE(IPT,*) '!  # UV OBSERVATION ERROR SPECIFIED                  :',OBSERR_UV
!     IF(T_OBS)  WRITE(IPT,*) '!  # T OBSERVATION ERROR SPECIFIED                   :',OBSERR_T
!     IF(S_OBS)  WRITE(IPT,*) '!  # S OBSERVATION ERROR SPECIFIED                   :',OBSERR_S
!     WRITE(IPT,*) '!                                                   !'

   CONTAINS
   
     SUBROUTINE NAME_LIST_INITIALIZE_ENKF
       IMPLICIT NONE
       ENKF_START_DATE = "RRK ASSIMILATION START AND END TIME"
       ENKF_END_DATE   = "For an idealized case specify 'seconds=(flt)','days=(flt)', or 'cycles=(int)'" 
       ENKF_ASSIM_INTERVAL= "A length of time: 'seconds= ','days= ', or 'cycles= '"   
       ENKF_NOBSMAX    = -1  
       ENKF_NENS       = -1 
       ENKF_CINF       =  1E20 
       EKINT_START=0
       ENKF_ON         = .FALSE.
       EL_ASSIM       = .FALSE.  
       EL_OBS         = .FALSE.  
       UV_ASSIM       = .FALSE.  
       UV_OBS         = .FALSE. 
       T_ASSIM        = .FALSE.   
       T_OBS          = .FALSE.  
       S_ASSIM        = .FALSE.  
       S_OBS          = .FALSE.
       ENKF_LOCALIZED = .FALSE.
       ENKF_METHOD      = 1
       MODE           = 11
       LOCAL_DISK     = .FALSE.
       OBSERR_EL    = 0.001
       OBSERR_UV    = 0.001
       OBSERR_T     = 0.01
       OBSERR_S     = 0.01


     END SUBROUTINE NAME_LIST_INITIALIZE_ENKF
     
     SUBROUTINE NAME_LIST_PRINT_ENKF
       IMPLICIT NONE

       write(UNIT=IPT,NML=NML_ENKF)

     RETURN
     END SUBROUTINE NAME_LIST_PRINT_ENKF
     
     SUBROUTINE NAME_LIST_READ_ENKF
       IMPLICIT NONE
       CHARACTER(LEN=120) :: FNAME
       INTEGER :: IOS

       FNAME = TRIM(INPUT_DIR)//"/"//trim(casename)//"_assim_enkf.nml"
       
       CALL FOPEN(KFUNIT,trim(FNAME),'cfr')

       ! Read RRKF parameters
       READ(UNIT=KFUNIT, NML=NML_ENKF,IOSTAT=ios)
       if(ios .NE. 0 ) then
          if(DBG_SET(dbg_log)) write(UNIT=IPT,NML=NML_ENKF)
          Call Fatal_Error("Can Not Read NameList NML_ENKF from file: "//trim(FNAME))
       end if
       
       REWIND(KFUNIT)
       
       if(DBG_SET(dbg_scl)) &
            & write(IPT,*) "Read_Name_List: NML_ENKF"
       
       if(DBG_SET(dbg_scl)) &
            & write(UNIT=IPT,NML=NML_ENKF)
       
       CLOSE(NMLUNIT)
!       write(*,NML=NML_ENKF) ! should be marked 
       if (ENKF_ASSIM_INTERVAL .ne. NC_OUT_INTERVAL)  then
           print *, "NC_OUT_INTERVAL must be same as ENKF_ASSIM_INTERVAL"
           print *, 'you can set the ENKF_ASSIM_INTERVAL smaller than your observation interal to fit your nc_output_interval even though you dont have observation data. '
           print *, 'stop'
           call pstop
       end if
       print *, 'finish reading namelist enkf'
       RETURN
     END SUBROUTINE NAME_LIST_READ_ENKF
     
     SUBROUTINE ENKF_SET_TIME
       USE MOD_SET_TIME
       IMPLICIT NONE
       CHARACTER(LEN=4) :: FLAG,BFLAG,EFLAG
       INTEGER(ITIME):: begins, ends,steps
       integer :: status


       
       ! GET THE RRK ASSIMILATION INTERVAL PERIOD
       CALL IDEAL_TIME_STRING2TIME(ENKF_ASSIM_INTERVAL,FLAG,ENKF_INTERVAL,steps)
       IF (FLAG == 'time') THEN ! IF START AND END TIME WERE SPECIFIED

          IF(MOD(ENKF_INTERVAL,IMDTI) /= ZeroTime) THEN
             CALL FATAL_ERROR("ENKF_ASSIM_INTERVAL must be an integer number of internal time steps!")
          END IF
          
          steps = seconds(ENKF_INTERVAL)/seconds(IMDTI)


       ELSEIF(FLAG == 'step') THEN ! IF SPECIFIED AS A NUMBER OF STEPS
          
          ENKF_INTERVAL = steps * IMDTI
          
       END IF

!       EKInt_start=0
       EKint_end = steps


       SELECT CASE(USE_REAL_WORLD_TIME)
       CASE(.TRUE.)

!          REF_START_TIME = READ_DATETIME(REF_START_DATE,DATE_FORMAT,TIMEZONE,status)
!          if (.not. status) &
!               & Call Fatal_Error("Could not read the date string REF_START_DATE: "//trim(REF_START_DATE))


!          ! GET THE END TIME
!          REF_END_Time = READ_DATETIME(REF_END_DATE,DATE_FORMAT,TIMEZONE,status)
!          if (.not. status) &
!               & Call Fatal_Error("Could not read the date string REF_END_DATE: "&
!               &//trim(REF_END_DATE))

!          ! SANITY CHECK
!          if(REF_Start_Time .GE. REF_End_Time) &
!               & Call Fatal_Error("RRKF REF_Start_Date exceeds or equal to REF_End_Date")
          
          ENKF_START_TIME = READ_DATETIME(ENKF_START_DATE,DATE_FORMAT,TIMEZONE,status)
!!$          if (.not. status) &
          if (status == 0) &
               & Call Fatal_Error("Could not read the date string ENKF_START_DATE: "//trim(ENKF_START_DATE))
          
          
          ! GET THE END TIME
          ENKF_END_Time = READ_DATETIME(ENKF_END_DATE,DATE_FORMAT,TIMEZONE,status)
!!$          if (.not. status) &
          if (status == 0) &
               & Call Fatal_Error("Could not read the date string ENKF_END_DATE: "&
               &//trim(ENKF_END_DATE))

          ! SANITY CHECK
          if(ENKF_START_TIME .GE. ENKF_END_Time) &
               & Call Fatal_Error("ENKF_START_DATE exceeds or equal to ENKF_END_DATE")
          


       CASE (.FALSE.)
!           print *, 'REF_START_DATE',REF_START_DATE
          ! GET THE REFERENCE START AND END INFORMATION
!          CALL IDEAL_TIME_STRING2TIME(REF_START_DATE,BFLAG,REF_START_TIME,begins)

!          CALL IDEAL_TIME_STRING2TIME(REF_END_DATE,EFLAG,REF_End_TIME,ends)

          ! SANITY CHECK
!          IF (BFLAG /= EFLAG) CALL FATAL_ERROR&
!               ('IDEALIZED MODEL TIME SPECIFICATION IS INCORRENT',&
!               &'RRK REF BEGIN AND END CAN BE IN EITHER CYCLES OR TIME BUT NOT MIXED',&
!               & trim(REF_start_date),trim(ref_end_date) )
          
!          IF (EFLAG == 'time') THEN ! IF START AND END TIME WERE SPECIFIED

             
             ! DO NOTHING

!          ELSE IF(EFLAG == 'step') THEN ! IF START AND END IINT WERE SPECIFIED


!             REF_START_TIME = StartTime + IMDTI*(begins - ISTART)

!            REF_END_TIME = StartTime + IMDTI*(ends - ISTART)

             
!          ELSE
!             print *, 'EFLAG:' , eflag
!             CALL FATAL_ERROR&
!                  &('IDEAL_TIME_STRING2TIME returned invalid flag for rrk ref date and time')
             
!          END IF

          ! GET THE RRK START AND END INFORMATION
          CALL IDEAL_TIME_STRING2TIME(ENKF_START_DATE,BFLAG,ENKF_START_TIME,begins)

          CALL IDEAL_TIME_STRING2TIME(ENKF_END_DATE,EFLAG,ENKF_END_Time,ends)

          ! SANITY CHECK
          IF (BFLAG /= EFLAG) CALL FATAL_ERROR&
               ('IDEALIZED MODEL TIME SPECIFICATION IS INCORRENT',&
               &'ENKF BEGIN AND END CAN BE IN EITHER CYCLES OR TIME BUT NOT MIXED',&
               & trim(enkf_start_date),trim(enkf_end_date) )
          
          IF (EFLAG == 'time') THEN ! IF START AND END TIME WERE SPECIFIED

             
             ! DO NOTHING

          ELSE IF(EFLAG == 'step') THEN ! IF START AND END IINT WERE SPECIFIED

!             ENKF_START_TIME = StartTime + IMDTI*(begins - ISTART)
             ENKF_START_TIME = StartTime + IMDTI*(begins - ISTART+1)

!             ENKF_END_Time = StartTime + IMDTI*(ends - ISTART)
             ENKF_END_Time = StartTime + IMDTI*(ends - ISTART+1)

             
          ELSE
             CALL FATAL_ERROR&
                  &('IDEAL_TIME_STRING2TIME returned invalid flag for rrk ref date and time')

          END IF



       END SELECT


!       IF(MOD(REF_START_TIME-StartTime,IMDTI) /= ZeroTime) THEN
!          CALL FATAL_ERROR("REF_START_TIME must be an integer number of internal time steps from the start time!")
!       END IF
       
!       IF(MOD(REF_END_TIME-StartTime,IMDTI) /= ZeroTime) THEN
!          CALL FATAL_ERROR("REF_END_TIME must be an integer number of internal time steps from the start time!")
!       END IF
       
       
       IF(MOD(ENKF_START_TIME-StartTime,IMDTI) /= ZeroTime) THEN
          CALL FATAL_ERROR("ENKF_START_TIME must be an integer number of internal time steps from the start time!")
       END IF
       
       IF(MOD(ENKF_END_Time-StartTime,IMDTI) /= ZeroTime) THEN
          CALL FATAL_ERROR("ENKF_END_Time must be an integer number of internal time steps from the start time!")
       END IF
       
       
     END SUBROUTINE ENKF_SET_TIME


!     SUBROUTINE RRK_SET_STARTUP
!       USE MOD_SET_TIME
!       USE MOD_INPUT
!       IMPLICIT NONE
!       character(len=160) :: pathnfile
!       TYPE(NCFILE), POINTER :: NCF
!
!       OPEN(75,FILE='./err.dat')
!
!       !RESTART FILE
!       RESTART_FILE_NAME = trim(casename)//"_restart_0001.nc"
!       pathnfile= trim(OUTPUT_DIR)//TRIM(RESTART_FILE_NAME)
!       
!       CALL SEARCH_FOR_LAST_MATCHING_NAME(PATHNFILE)
!       
!       ! OPEN THE FILE AND LOAD FOR STARTUP
!       NCF => NEW_FILE()
!       NCF%FNAME=trim(pathnfile)
!       Call NC_OPEN(NCF)
!       CALL NC_LOAD(NCF)
!       NC_START => NCF
!       
!       Nullify(NCF)
!       CALL SET_STARTUP_FILE_STACK(REF_START_TIME)
!
!!       STARTUP_TYPE = STARTUP_TYPE_CRASHRESTART
!
!!       STARTUP_UV_TYPE = STARTUP_TYPE_SETVALUES
!
!!       STARTUP_TURB_TYPE = STARTUP_TYPE_SETVALUES
!
!!       STARTUP_TS_TYPE = STARTUP_TYPE_SETVALUES
!
!     END SUBROUTINE RRK_SET_STARTUP


!!=====================================================================================
!!
!! READ IN THE SIMULATION DATA AND CALCULATE THE EOF
!!
!!=====================================================================================

   
    SUBROUTINE ENKF_REF 
!     USE CONTROL
!     USE ENKFVAL
!     USE MOD_NCTOOLS
     IMPLICIT NONE
  
!     TYPE(NCFILE), POINTER :: NCF
!     TYPE(TIME) :: tTEST
!     LOGICAL :: FOUND
!     INTEGER :: NTIMES
     
!     INTEGER IDUMMY, I, K, J
!     INTEGER SS_DIM
!     INTEGER STDIM
!  !   REAL(DP) ::  ELSD, USD, TSD, SSD
!     INTEGER LWORK4, LDVT
!     REAL(DP),ALLOCATABLE :: WORK4(:)
!     REAL(DP) :: VT
!     REAL(DP),ALLOCATABLE :: RKSF(:,:)
!     REAL(DP),ALLOCATABLE :: RKSF1(:,:)
!  !   REAL(DP),ALLOCATABLE :: SEOF(:,:)
!     REAL(DP),ALLOCATABLE :: SFSF(:,:)
!     REAL(DP),ALLOCATABLE :: SFD(:)
!     REAL(DP),ALLOCATABLE :: SFU(:,:)
!     REAL(DP),ALLOCATABLE :: STVAR(:)
!  !   REAL(DP),ALLOCATABLE :: STSD(:)
!     REAL(DP) :: SUM0   
!     INTEGER RCODE
!     INTEGER IDUMMY1
!     REAL(DP),ALLOCATABLE :: RRKU2(:,:,:)
!     REAL(DP),ALLOCATABLE :: RRKV2(:,:,:)
!     REAL(DP),ALLOCATABLE :: RRKEL2(:,:)
!     REAL(DP),ALLOCATABLE :: RRKT2(:,:,:)
!     REAL(DP),ALLOCATABLE :: RRKS2(:,:,:)
  !   REAL(DP),ALLOCATABLE :: STMEAN(:)
  
!     TYPE(NCDIM), POINTER :: DIM
     TYPE(NCVAR), POINTER :: VAR
  
!     INTEGER :: IRANGE(2),CNT, STATUS
!     character(len=160) :: pathnfile 
    
!    pathnfile = trim(INPUT_DIR)//trim(casename)//"_true_0001.nc"
!    
!    NCF => NEW_FILE()
!    NCF%FNAME=trim(pathnfile)
!    
!    Call NC_OPEN(NCF)
!    CALL NC_LOAD(NCF)
! 
!    ! GET THE TIME DIMENSION
!    DIM => FIND_UNLIMITED(NCF,FOUND)
!    IF(.not. FOUND) CALL FATAL_ERROR &
!             & ("IN RRK_REF FILE",&
!             & "FILE NAME: "//TRIM(NCF%FNAME),&
!             &"COULD NOT FIND THE UNLIMITED DIMENSION")
!    NTIMES = DIM%DIM
!    
! 
!    ! TEST THE FILE START AND END
!    TTEST = get_file_time(NCF,1)
!    
!    IF(REF_START_TIME < TTEST) CALL FATAL_ERROR &
!         &("RRK_REF: REF_START_TIME IS BEFORE FIRST TIME IN FILE.",&
!         & "FILE NAME:"//TRIM(NCF%FNAME))
!    
!    
!    TTEST = get_file_time(NCF,NTIMES)
!    IF(REF_END_TIME > TTEST) CALL FATAL_ERROR &
!         &("RRK_REF: REF_START_TIME IS BEFORE FIRST TIME IN FILE.",&
!         & "FILE NAME:"//TRIM(NCF%FNAME))
! 
! 
!    ! GET THE START AND END INDEX
! 
!    IRANGE = 0
! 
!    CNT = 0
!    DO WHILE (product(irange)==0)
!       CNT = CNT +1
! 
!       IF (CNT > NTIMES) THEN
!          CALL PRINT_TIME(REF_START_TIME,IPT,"RRK_REF START TIME")
!          CALL PRINT_TIME(REF_END_TIME,IPT,"RRK_REF END TIME")
!          
! 
!          CALL FATAL_ERROR&
!            &("RRK_REF: CAN NOT FIND SPECIFIED START AND END TIME IN THE REFERENCE FILE!", &
!            & "FILE NAME:"//TRIM(NCF%FNAME))
!       END IF
! 
! 
!       TTEST = get_file_time(NCF,CNT)
!       IF (REF_START_TIME == TTEST) THEN
!          IRANGE(1) = CNT
!       END IF
! 
!       IF (REF_END_TIME == TTEST) THEN
!          IRANGE(2) = CNT
!       END IF
! 
!    END DO
! 
! 
!    nrange = irange(2)-irange(1)+1
 
!    IF(EL_ASSIM) THEN
!       
!       VAR => FIND_VAR(NCF,"zeta",FOUND)
!       IF(.not. FOUND) CALL FATAL_ERROR &
!            & ("IN RRK_REF FILE",&
!             & "FILE NAME: "//TRIM(NCF%FNAME),&
!             &"COULD NOT FIND THE VARIABLE 'zeta'")
! 
!       allocate(RRKEL_REF(MGL,nrange),stat=status)
!       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY FOR RRKEL_REF")
!       RRKEL_REF = 0.0_SP
! 
!       CALL nc_connect_avar(VAR,RRKEL_REF)
! 
!       CALL NC_READ_VAR(var,stkrng=irange,parallel=.false.)
! 
! 
!    ENDIF
!    
!    IF(UV_ASSIM) THEN
! 
!       VAR => FIND_VAR(NCF,"u",FOUND)
!       IF(.not. FOUND) CALL FATAL_ERROR &
!            & ("IN RRK_REF FILE",&
!             & "FILE NAME: "//TRIM(NCF%FNAME),&
!             &"COULD NOT FIND THE VARIABLE 'u'")
! 
!       allocate(RRKU_REF(NGL,kbm1,nrange),stat=status)
!       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY FOR RRKU_REF")
!       RRKU_REF = 0.0_SP
!       CALL nc_connect_avar(VAR,RRKU_REF)
! 
!       CALL NC_READ_VAR(var,stkrng=irange,parallel=.false.)
! 
!       VAR => FIND_VAR(NCF,"v",FOUND)
!       IF(.not. FOUND) CALL FATAL_ERROR &
!            & ("IN RRK_REF FILE",&
!             & "FILE NAME: "//TRIM(NCF%FNAME),&
!             &"COULD NOT FIND THE VARIABLE 'v'")
! 
!       allocate(RRKV_REF(NGL,kbm1,nrange),stat=status)
!       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY FOR RRKV_REF")
!       RRKV_REF = 0.0_SP
!       CALL nc_connect_avar(VAR,RRKV_REF)
! 
!       CALL NC_READ_VAR(var,stkrng=irange,parallel=.false.)
! 
!       
!    ENDIF
!    
!    IF(T_ASSIM) THEN
!       VAR => FIND_VAR(NCF,"temp",FOUND)
!       IF(.not. FOUND) CALL FATAL_ERROR &
!            & ("IN RRK_REF FILE",&
!             & "FILE NAME: "//TRIM(NCF%FNAME),&
!             &"COULD NOT FIND THE VARIABLE 'temp'")
! 
!       allocate(RRKT_REF(MGL,kbm1,nrange),stat=status)
!       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY FOR RRKT_REF")
!       RRKT_REF = 0.0_SP
!       CALL nc_connect_avar(VAR,RRKT_REF)
! 
!       CALL NC_READ_VAR(var,stkrng=irange,parallel=.false.)
!       
!    ENDIF
!    
!    IF(S_ASSIM) THEN
!       VAR => FIND_VAR(NCF,"salinity",FOUND)
!       IF(.not. FOUND) CALL FATAL_ERROR &
!            & ("IN RRK_REF FILE",&
!             & "FILE NAME: "//TRIM(NCF%FNAME),&
!             &"COULD NOT FIND THE VARIABLE 'salinity'")
! 
!       allocate(RRKS_REF(MGL,kbm1,nrange),stat=status)
!       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY FOR RRKS_REF")
!       RRKS_REF = 0.0_SP
!       CALL nc_connect_avar(VAR,RRKS_REF)
! 
!       CALL NC_READ_VAR(var,stkrng=irange,parallel=.false.)
!       
!    ENDIF
    
!    CALL KILL_FILE(NCF)
!    
!    WRITE(IPT,*) 'Calculate the EOFs from the control run......'
!    
!    STDIM = 0
!    IF(EL_ASSIM) STDIM = STDIM + MGL
!    IF(UV_ASSIM) STDIM = STDIM + 2*NGL*KBM1
!    IF(T_ASSIM)  STDIM = STDIM + MGL*KBM1
!    IF(S_ASSIM)  STDIM = STDIM + MGL*KBM1
     
!    SS_DIM = nrange
!    LWORK4 = 5*SS_DIM
!    LDVT   = 1   
!    
!    ALLOCATE(WORK4(LWORK4))             ;WORK4   = ZERO
!    ALLOCATE(RKSF(STDIM,SS_DIM))        ;RKSF    = ZERO
!    ALLOCATE(RKSF1(STDIM,SS_DIM))       ;RKSF1   = ZERO
!    ALLOCATE(SEOF(STDIM,SS_DIM))        ;SEOF    = ZERO
!    ALLOCATE(SFSF(SS_DIM,SS_DIM))       ;SFSF    = ZERO
!    ALLOCATE(SFD(SS_DIM))               ;SFD     = ZERO
!    ALLOCATE(SFU(SS_DIM,SS_DIM))        ;SFU     = ZERO
!    ALLOCATE(STVAR(STDIM))              ;STVAR   = ZERO
!    ALLOCATE(STSD(STDIM))               ;STSD    = ZERO
!    ALLOCATE(RRKU(NGL,KBM1))            ;RRKU    = ZERO
!    ALLOCATE(RRKV(NGL,KBM1))            ;RRKV    = ZERO
!    ALLOCATE(RRKEL(MGL))                ;RRKEL   = ZERO
!    ALLOCATE(RRKT(MGL,KBM1))            ;RRKT    = ZERO
!    ALLOCATE(RRKS(MGL,KBM1))            ;RRKS    = ZERO
!    ALLOCATE(RRKU2(NGL,KBM1,SS_DIM))    ;RRKU2   = ZERO
!    ALLOCATE(RRKV2(NGL,KBM1,SS_DIM))    ;RRKV2   = ZERO
!    ALLOCATE(RRKEL2(MGL,SS_DIM))        ;RRKEL2  = ZERO
!    ALLOCATE(RRKT2(MGL,KBM1,SS_DIM))    ;RRKT2   = ZERO
!    ALLOCATE(RRKS2(MGL,KBM1,SS_DIM))    ;RRKS2   = ZERO
!    ALLOCATE(STMEAN(STDIM))             ;STMEAN  = ZERO  
    
    
    ! STORE IN RKSF
 
     
!     DO I = 1, SS_DIM
!   
!       IDUMMY = 0      
! 
!       IF(EL_ASSIM) THEN      
!         DO J = 1, MGL
!           IDUMMY = IDUMMY + 1
! 	  RKSF(IDUMMY,I) =  RRKEL_REF(J,I)
!         ENDDO 
!       ENDIF
! 
!       IF(UV_ASSIM) THEN
!         DO K = 1, KBM1
!           DO J = 1, NGL
! 	     IDUMMY = IDUMMY + 1
! 	     RKSF(IDUMMY,I) = RRKU_REF(J,K,I)
! 	  ENDDO
!         ENDDO	
! 
!         DO K = 1, KBM1
!           DO J = 1, NGL
! 	     IDUMMY = IDUMMY + 1
! 	     RKSF(IDUMMY,I) = RRKV_REF(J,K,I)
! 	  ENDDO
!         ENDDO     
!       ENDIF
! 
!       IF(T_ASSIM) THEN
!         DO K = 1, KBM1
!           DO J = 1, MGL
! 	     IDUMMY = IDUMMY + 1
! 	     RKSF(IDUMMY,I) = RRKT_REF(J,K,I)
! 	  ENDDO
!         ENDDO	
!       ENDIF
! 
!       IF(S_ASSIM) THEN
! 	DO K = 1, KBM1
!           DO J = 1, MGL
! 	     IDUMMY = IDUMMY + 1
! 	     RKSF(IDUMMY,I) = RRKS_REF(J,K,I)
! 	  ENDDO
!         ENDDO	
!       ENDIF
!       
!     ENDDO
!     
!       
! ! CALCULATE THE MEAN AND THE STANDARD DEVIATION
! 
!     DO I=1, STDIM
!       SUM0=0.0_DP
!       DO J=1, SS_DIM
!         SUM0=SUM0+RKSF(I,J)
!       ENDDO
!      
!       STMEAN(I) = SUM0/DBLE(SS_DIM)
!       SUM0=0.0_DP
!       DO J=1, SS_DIM
!         RKSF1(I,J)=(RKSF(I,J)-STMEAN(I))
! 	SUM0=SUM0+RKSF1(I,J)**2.0_DP
!       ENDDO
!            
!       STVAR(I)=SUM0
!       STSD(I)=DSQRT(SUM0/DBLE(SS_DIM-1))
!     ENDDO    
!     
! 
!     IDUMMY  = 0
!     IDUMMY1 = 0
!     ELSD    = 0.0_DP
!     USD     = 0.0_DP
!     TSD     = 0.0_DP
!     SSD     = 0.0_DP          
! 
!     IF(EL_ASSIM) THEN
!       SUM0=0.0_DP
!       DO I=1, MGL
!         IDUMMY = IDUMMY + 1
!         SUM0 = SUM0 + STVAR(IDUMMY)
!       ENDDO
!       
!       ELSD=DSQRT(SUM0/DBLE(MGL)/DBLE(SS_DIM-1))
! 
!       DO I=1, MGL
!         DO J=1, SS_DIM
!           RKSF1(I,J)=RKSF1(I,J)/ELSD/DSQRT(DBLE(SS_DIM-1)) 
!         ENDDO
!       ENDDO
!     
!       IDUMMY1 = IDUMMY
!     ENDIF
! 
!     IF(UV_ASSIM) THEN
!       SUM0=0.0_DP
!       DO K=1, 2*KBM1
!         DO I=1, NGL
!            IDUMMY = IDUMMY + 1
! 	   SUM0 = SUM0 +STVAR(IDUMMY)
!         ENDDO
!       ENDDO
!       USD=DSQRT(SUM0/DBLE(NGL*KBM1)/DBLE(SS_DIM-1))
!  
!       DO I=IDUMMY1+1, IDUMMY1+2*NGL*KBM1 
!         DO J=1, SS_DIM	    
! 	   RKSF1(I,J)=RKSF1(I,J)/USD/DSQRT(DBLE(KBM1))/DSQRT(DBLE(SS_DIM-1))         
!         ENDDO 
!       ENDDO
!     
!       IDUMMY1 = IDUMMY
!     ENDIF  
!       
!     IF(T_ASSIM) THEN
!       SUM0=0.0_DP
!       DO K=1, KBM1
!         DO I=1, MGL
!            IDUMMY = IDUMMY + 1
! 	   SUM0 = SUM0 +STVAR(IDUMMY)
!         ENDDO
!       ENDDO
!       TSD=DSQRT(SUM0/DBLE(MGL*KBM1)/DBLE(SS_DIM-1))
!  
!       DO I=IDUMMY1+1, IDUMMY1+MGL*KBM1 
!         DO J=1, SS_DIM
! 	  RKSF1(I,J)=RKSF1(I,J)/TSD/DSQRT(DBLE(SS_DIM-1))
!         ENDDO 
!       ENDDO
! 
!       IDUMMY1 = IDUMMY
!     ENDIF        
!       
!     IF(S_ASSIM) THEN
!       SUM0=0.0_DP
!       DO K=1, KBM1
!         DO I=1, MGL
!            IDUMMY = IDUMMY + 1
! 	   SUM0 = SUM0 +STVAR(IDUMMY)
!         ENDDO
!       ENDDO
!       SSD=DSQRT(SUM0/DBLE(MGL*KBM1)/DBLE(SS_DIM-1))
!  
!       DO I=IDUMMY1+1, IDUMMY1+MGL*KBM1 
!         DO J=1, SS_DIM
! 	   RKSF1(I,J)=RKSF1(I,J)/SSD/DSQRT(DBLE(SS_DIM-1))
!         ENDDO 
!       ENDDO
! 
!       IDUMMY1 = IDUMMY
!     ENDIF       
! 
! 
! ! CALCULATE THE EOFs BY SVD OF RKSF1' RKSF1 = C LAMBDA C', RKSF1 RKSF1' = E LAMBDA E', E = RKSF1 C LAMBDA ^-1/2
! 
!     DO I=1, SS_DIM
!       DO J=1, SS_DIM
!         SUM0=0.0_DP
! 	DO K=1, STDIM
! 	  SUM0 = SUM0 + RKSF1(K,I)*RKSF1(K,J)
! 	ENDDO
! 	SFSF(I,J) =SUM0
!       ENDDO
!     ENDDO
!  
!     CALL DGESVD('A','N',SS_DIM,SS_DIM,SFSF,SS_DIM,SFD,SFU,SS_DIM,VT,LDVT,WORK4,LWORK4,RCODE)    
!     
!     OPEN(72,FILE=TRIM(OUTPUT_DIR)//'/rrktemp/'//'eigenvalue.dat')
!     
!     DO I=1, SS_DIM
!       WRITE(72,'(I5,E15.7)') I, SFD(I)   
!     ENDDO
!     DO I=1, SS_DIM
!       DO J=1, SS_DIM 
!         SFU(I,J)=SFU(I,J)/DSQRT(SFD(J))
!       ENDDO
!     ENDDO
! 
!     DO I=1, STDIM
!       DO J=1, SS_DIM
!          SUM0=0.0_DP 
!          DO K=1, SS_DIM
! 	    SUM0=SUM0+RKSF1(I,K)*SFU(K,J)
! 	 ENDDO
!          SEOF(I,J)=SUM0 
!       ENDDO
!     ENDDO
     
 !    write(1000,'(871f12.4)') (seof(i,1),i=1,stdim)
     
!     DEALLOCATE(WORK4,RKSF,RKSF1,SFSF,SFD,SFU,STVAR)
!     DEALLOCATE(RRKU2,RRKV2,RRKEL2,RRKT2,RRKS2)
 
  END SUBROUTINE ENKF_REF
 
!------------------------------------------------------------------------------|
!  Main program to calculate the reduced rank kalman gain matrix               |
!------------------------------------------------------------------------------|   
   
!   SUBROUTINE RRK_RRK(CHOICE)
!
!   USE LIMS
!   USE CONTROL
!   USE ENKFVAL
!   IMPLICIT NONE
!    
!    INTEGER SS_DIM
!    INTEGER STDIM
!    INTEGER I_EOF
!    INTEGER I,J,K
!    INTEGER CHOICE
!    INTEGER NLOC
!    INTEGER RCODE
!    INTEGER IT
!    INTEGER ILAST
!    INTEGER IDUMMY
!    INTEGER SS_TOT
!    INTEGER,ALLOCATABLE  :: STLOC(:)
!    REAL(DP) ERR1
!    REAL(DP),ALLOCATABLE :: KAL(:,:)
!    REAL(DP),ALLOCATABLE :: HEOF(:,:)
!    REAL(DP),ALLOCATABLE :: R(:,:)
!    REAL(DP),ALLOCATABLE :: HTR(:,:)
!    REAL(DP),ALLOCATABLE :: R2(:,:)
!    REAL(DP),ALLOCATABLE :: PHT(:,:)
!    REAL(DP),ALLOCATABLE :: RALPHA(:,:)
!    REAL(DP),ALLOCATABLE :: BETA2(:,:)
!    REAL(DP),ALLOCATABLE :: GAMMA(:,:)
!    REAL(DP),ALLOCATABLE :: GAMMA2(:,:)
!    REAL(DP),ALLOCATABLE :: TRANS2(:,:)
!    
!! WORK ARRAYS FOR LAPACK SUBROUTINES
!    INTEGER LWORK, LWORK2, LDVT
!    INTEGER,ALLOCATABLE :: IPIV(:)
!    INTEGER,ALLOCATABLE :: IPIV2(:)
!    REAL(DP),ALLOCATABLE :: WORK(:)
!    REAL(DP),ALLOCATABLE :: WORK2(:)
!    REAL(DP) :: VT
!    
!
!! CHARACTER STRINGS
!    CHARACTER(LEN=80) KALNAME
!    CHARACTER(LEN=80) FILENAME
!    CHARACTER(LEN=4)  STRCYC
!    CHARACTER(LEN=8)  CH8
!    CHARACTER(LEN=80) INAME    
!    CHARACTER(LEN=80) INAME2
!    
!!    REAL(DP),ALLOCATABLE :: SEOF(:,:)
!!    REAL(DP),ALLOCATABLE :: STSD(:)
!!    REAL(DP) :: ELSD, USD, TSD, SSD
!    REAL(DP),ALLOCATABLE :: HU(:,:)
!    REAL(DP),ALLOCATABLE :: HUL(:,:)
!    REAL(DP),ALLOCATABLE :: EVAL(:)
!    REAL(DP) :: RRKSUM
!      
!    STDIM = 0
!    IF(EL_ASSIM) STDIM = STDIM + MGL
!    IF(UV_ASSIM) STDIM = STDIM + 2*NGL*KBM1
!    IF(T_ASSIM)  STDIM = STDIM + MGL*KBM1
!    IF(S_ASSIM)  STDIM = STDIM + MGL*KBM1
!
!    SS_DIM = ENKF_NENS
!    LWORK  = 4*SS_DIM
!    LWORK2 = 4*ENKF_NOBSMAX
!    SS_TOT = nrange
!    LDVT   = 1 
!    ILAST  = 10  ! THE NUMBER OF ITERATION IN DOUBLING ALGORITHM
!
!! TEMPORARILY ALLOCATE ARRY TO ARRYS
!
!
!   ALLOCATE(STLOC(ENKF_NOBSMAX))            ;STLOC    = 0
!   ALLOCATE(IPIV(SS_DIM))                  ;IPIV     = 0
!   ALLOCATE(IPIV2(ENKF_NOBSMAX))            ;IPIV2    = 0
!   ALLOCATE(KAL(SS_DIM,ENKF_NOBSMAX))       ;KAL      = ZERO
!   ALLOCATE(HEOF(ENKF_NOBSMAX,SS_DIM))      ;HEOF     = ZERO
!   ALLOCATE(R(ENKF_NOBSMAX,ENKF_NOBSMAX))    ;R        = ZERO
!   ALLOCATE(HTR(SS_DIM,ENKF_NOBSMAX))       ;HTR      = ZERO
!   ALLOCATE(R2(ENKF_NOBSMAX,ENKF_NOBSMAX))   ;R2       = ZERO
!   ALLOCATE(PHT(SS_DIM,ENKF_NOBSMAX))       ;PHT      = ZERO
!   ALLOCATE(RALPHA(SS_DIM,SS_DIM))         ;RALPHA   = ZERO
!   ALLOCATE(BETA2(SS_DIM,SS_DIM))          ;BETA2    = ZERO
!   ALLOCATE(GAMMA(SS_DIM,SS_DIM))          ;GAMMA    = ZERO
!   ALLOCATE(GAMMA2(SS_DIM,SS_DIM))         ;GAMMA2   = ZERO
!   ALLOCATE(TRANS2(SS_DIM,SS_DIM))         ;TRANS2   = ZERO 
!   ALLOCATE(WORK(LWORK))                   ;WORK     = ZERO
!   ALLOCATE(WORK2(LWORK2))                 ;WORK2    = ZERO
!!   ALLOCATE(SEOF(STDIM,SS_DIM))            ;SEOF     = ZERO
!!   ALLOCATE(STSD(STDIM))                   ;STSD     = ZERO
!   ALLOCATE(HU(ENKF_NOBSMAX,SS_TOT))        ;HU       = ZERO      
!   ALLOCATE(HUL(ENKF_NOBSMAX,SS_TOT))       ;HUL      = ZERO
!   ALLOCATE(EVAL(SS_TOT))                  ;EVAL     = ZERO
!
!
!! END ALLOCATION
!
!   IF(CHOICE == 1) THEN
!     IF(MSR) WRITE(IPT,*) 'Starting perturb the base state in the eigenvector direction......'
!     
!     CALL PERTURB
!
!     DEALLOCATE(STLOC,IPIV,IPIV2,KAL,HEOF,R,HTR,R2,PHT,RALPHA,TRANS2)
!     DEALLOCATE(BETA2,GAMMA,GAMMA2,WORK,WORK2,HU,HUL,EVAL) 
!
!     RETURN
!   ENDIF
!
!   IF(CHOICE == 2) THEN
!     IF(MSR) WRITE(IPT,*) 'Starting calculate the linearized model matrix in the reduced space......'  
!! CALCULATE THE LINEARIZED MODEL MATRIX IN THE REDUCED SPACE
!     CALL MREDUCED
!     
!     IF(MSR) THEN
!! OUTPUT LINEAR TRANSITION MATRIX   
!       FILENAME = TRIM(OUTPUT_DIR)//'/rrktemp/Amatr.dat' 
!       OPEN(INORRK,FILE=TRIM(FILENAME),FORM='UNFORMATTED')
!       DO J=1, SS_DIM
!         WRITE(INORRK) (TRANS(I,J), I=1, SS_DIM)    
!       ENDDO 
!       CLOSE(INORRK)
!
! ! ALSO IN ASCII
!       FILENAME = TRIM(OUTPUT_DIR)//'/rrktemp/Amatr1.txt'
!       OPEN(INORRK, FILE=TRIM(FILENAME),FORM='FORMATTED')
!       DO I=1, SS_DIM
!         DO J=1, SS_DIM
!           WRITE(INORRK,*) TRANS(J,I)
!         ENDDO
!       ENDDO
!       CLOSE(INORRK)
!     ENDIF  
!
!     DEALLOCATE(STLOC,IPIV,IPIV2,KAL,HEOF,R,HTR,R2,PHT,RALPHA,TRANS2)
!     DEALLOCATE(BETA2,GAMMA,GAMMA2,WORK,WORK2,HU,HUL,EVAL) 
!!     DEALLOCATE(BETA2,GAMMA,GAMMA2,WORK,WORK2,SEOF,STSD,HU,HUL,EVAL)     
!!     DEALLOCATE(STTEMP0,STTEMP1,STEOF,SDEOF,TRANS,RRKEL,RRKU,RRKV,RRKT,RRKS)
!     RETURN
!   ENDIF
!       
!   IF(CHOICE == 4) THEN
!     CALL RRK_ALLOC_VAR
!     WRITE(IPT,*) 'Calculate the Kalman gain matrix by doubling algorith......'
!     KALNAME = TRIM(OUTPUT_DIR)//'/rrktemp/rrK.dat'
!
!! READ THE EIGENVALUES AND SET THE INITIAL ERROR COVARIANCE   
!     INAME = TRIM(OUTPUT_DIR)//'/rrktemp/eigenvalue.dat' 
!     OPEN(INORRK,FILE=INAME,STATUS='OLD')
!     
!      WRITE(IPT,*) 'ss_tot:', ss_tot, ss_dim     
!     
!     DO I=1, SS_TOT
!       READ(INORRK,*) J, EVAL(I)
!     ENDDO
!     
!     DO I=1, SS_DIM
!       DO J=1, SS_DIM
!         GAMMA2(I,J) = 0._SP
!       ENDDO
!       GAMMA2(I,I) = EVAL(I)*RRK_PSCALE
!     ENDDO
!     CLOSE(INORRK)
!   
!! READ THE LINEAR TRANSITION MATRIX     
!     FILENAME = TRIM(OUTPUT_DIR)//'/rrktemp/Amatr.dat' 
!     OPEN(INORRK,FILE=TRIM(FILENAME),FORM='UNFORMATTED')
!     DO J=1, SS_DIM
!       READ(INORRK) (TRANS(I,J), I=1, SS_DIM)    
!     ENDDO 
!     CLOSE(INORRK)
!     
!! READ THE OBSERVATION NUMBER AND LOCATION   
!     CALL READOBS(STLOC,NLOC)
!     
!     print *, 'obs=', STLOC(1), NLOC
!   
!! SET RALPHA   
!     DO I=1, SS_DIM
!       DO J=1, SS_DIM
!         TRANS2(I,J)=TRANS(J,I)
!       ENDDO
!!      WRITE(IPT,*) 'Amatr:', I, TRANS(I,I)     
!     ENDDO
!
!     print *, 'debug - 4'
!       
!! CALCULATE OBSERVATION MATRIX (EOF TO OBS TRANSFORMATION) H*D*EOFs
!     DO J=1, SS_DIM
!       DO I=1, STDIM  
!          STTEMP1(I) = SEOF(I,J)
!       ENDDO
!       
!       DO I=1, NLOC         
!         HEOF(I,J)=STTEMP1(STLOC(I))*STSD(STLOC(I))
!        WRITE(IPT,*) 'Heof:', I,J,HEOF(I,J)        
!       ENDDO
!     ENDDO
!    
!
!     FILENAME = TRIM(OUTPUT_DIR)//'/rrktemp/ObsOp.dat'
!     OPEN(INORRK,FILE=FILENAME, FORM='UNFORMATTED')
!     DO J=1, SS_DIM
!       WRITE(INORRK) (HEOF(I,J), I=1, NLOC)    
!     ENDDO 
!     CLOSE(INORRK)     
!     FILENAME = TRIM(OUTPUT_DIR)//'/rrktemp/Heof.txt'
!     OPEN(INORRK,FILE=FILENAME, FORM='FORMATTED')
!     DO J=1, SS_DIM
!       DO I=1, NLOC
!         WRITE(INORRK,*) I,J,HEOF(I,J)
!       ENDDO    
!     ENDDO 
!     CLOSE(INORRK)
!
!! OUTPUT RALPHA MATRIX, WHICH IS TRANSPOSE OF TRANS
!     FILENAME = TRIM(OUTPUT_DIR)//'/rrktemp/Alpha.dat'
!     OPEN(INORRK,FILE=FILENAME, FORM='UNFORMATTED')
!     DO J=1, SS_DIM
!       WRITE(INORRK) (TRANS2(I,J), I=1, SS_DIM)    
!     ENDDO 
!     CLOSE(INORRK)     
!
!! OUTPUT GAMMA MATRIX
!     FILENAME = TRIM(OUTPUT_DIR)//'/rrktemp/Gamma.dat'
!     OPEN(INORRK,FILE=FILENAME, FORM='UNFORMATTED')
!     DO J=1, SS_DIM
!       WRITE(INORRK) (GAMMA2(I,J), I=1, SS_DIM)    
!     ENDDO 
!     CLOSE(INORRK) 
!     
!! CALCULATE REPRESENTATIVE OBSERVATION ERROR DIRECTLY USING EOFs IN UNRESOLVED SUBSPACE
!     DO I=SS_DIM+1, SS_TOT
!!       CALL READEOF(I,2)
!       DO J=1, STDIM  
!          STTEMP1(J) = SEOF(J,I)
!       ENDDO
!       
!       DO J=1, NLOC
!          HU(J,I) = STTEMP1(STLOC(J))*STSD(STLOC(J))       
!       ENDDO
!     ENDDO
!     DO I=1, NLOC
!       RRKSUM=0.0_DP
!       DO J=SS_DIM+1, SS_TOT
!         HUL(I,J)=RRKSUM+HU(I,J)*EVAL(J)
!       ENDDO
!     ENDDO
!     
!     DO I=1, NLOC
!       DO J=1, NLOC
!          RRKSUM=0.0_DP
!	  DO K=SS_DIM+1,SS_TOT
!             RRKSUM = RRKSUM + HUL(I,K)*HU(J,K)
!          ENDDO
!          R(I,J) = RRKSUM
!       ENDDO
!       
!! ASSUME THAT THE MEASUREMENT ERROR IS 1% OF THE ONE STANDARD DEVIATION OF THE CONTROL RUN
!       IF(RRK_OPTION == 0) THEN
!          R(I,I) = R(I,I) + (STSD(STLOC(I))*0.01_DP)**2.0_DP
!       ELSE
!          R(I,I) = R(I,I) + (STSD(STLOC(I))/DSQRT(DBLE(KBM1))*0.01_DP)**2.0_DP	  
!       ENDIF
!     ENDDO
!   
!     DO I=1, NLOC
!       DO J=1, I-1
!         R(I,J) = 0.5_DP*(R(I,J)+R(J,I))
!         R(J,I) = R(I,J)
!       ENDDO
!     ENDDO
!     FILENAME = TRIM(OUTPUT_DIR)//'/rrktemp/Robs.dat'
!     OPEN(INORRK,FILE=FILENAME, FORM='UNFORMATTED')
!     DO J=1, NLOC
!       WRITE(INORRK) (R(I,J), I=1, NLOC)    
!     ENDDO 
!     CLOSE(INORRK)  
!     FILENAME = TRIM(OUTPUT_DIR)//'/rrktemp/Robs.txt'
!     OPEN(INORRK,FILE=FILENAME, FORM='FORMATTED')
!     DO J=1, NLOC
!       DO I=1, NLOC
!         WRITE(INORRK,*) I,J, R(I,J)    
!       ENDDO	 
!     ENDDO 
!     CLOSE(INORRK)   
!     
!! INVERT OBSERVATION ERROR COVARIANCE     
!     CALL DGETRF(NLOC,NLOC,R,ENKF_NOBSMAX,IPIV2,RCODE)
!     IF(RCODE/=0) WRITE(IPT,*) 'error in computing LU factorization'     
!     call DGETRI(NLOC,R,ENKF_NOBSMAX,IPIV2,WORK2,LWORK2,RCODE)
!     IF(RCODE/=0) WRITE(IPT,*) 'error in inverting the matrix'
!
!! FORM BETA MATRIX = H_T*R**(-1)*H, WHERE H IS A NORMALIZED MATRIX H = H_ORIG*S
!     CALL DGEMM('t','n',SS_DIM,NLOC,NLOC,1.0d0,HEOF,ENKF_NOBSMAX,R,ENKF_NOBSMAX,0.0d0,HTR,SS_DIM) 
!     CALL DGEMM('n','n',SS_DIM,SS_DIM,NLOC,1.0d0,HTR,SS_DIM,HEOF,ENKF_NOBSMAX,0.0d0,BETA2,SS_DIM)    
!     FILENAME = TRIM(OUTPUT_DIR)//'/rrktemp/Beta.dat'
!     OPEN(INORRK,FILE=FILENAME, FORM='UNFORMATTED')
!     DO J=1, SS_DIM
!       WRITE(INORRK) (BETA2(I,J), I=1, SS_DIM)    
!     ENDDO 
!     CLOSE(INORRK)     
!     
!     WRITE(IPT,*) 'Running the doubling algorithm......'
!   
!     DO IT =1, ILAST 
!       
!! READ GAMMA FROM FILE
!       FILENAME =  TRIM(OUTPUT_DIR)//'/rrktemp/Gamma.dat'    
!       OPEN(INORRK,FILE=TRIM(FILENAME),FORM='UNFORMATTED')
!       DO J=1, SS_DIM
!         READ(INORRK) (GAMMA(I,J), I=1, SS_DIM)    
!       ENDDO 
!       CLOSE(INORRK)
!                    
!! READ BETA FROM FILE
!       FILENAME =  TRIM(OUTPUT_DIR)//'/rrktemp/Beta.dat'    
!       OPEN(INORRK,FILE=TRIM(FILENAME),FORM='UNFORMATTED')
!       DO J=1, SS_DIM
!         READ(INORRK) (BETA2(I,J), I=1, SS_DIM)    
!       ENDDO 
!       CLOSE(INORRK)
!      
!! COMPUTE EYE + BETA*GAMMA = STORE TEMPORARILY IN RALPHA (RALPHA = EYE + BETA*GAMMA)
!         
!       DO I=1, SS_DIM
!         DO J=1, SS_DIM
!           RALPHA(I,J)=0.0_DP
!         ENDDO
!         RALPHA(I,I)=1.0_DP
!       ENDDO     
!       CALL DGEMM('n','n',SS_DIM,SS_DIM,SS_DIM,1.0d0,BETA2,SS_DIM,GAMMA,SS_DIM,1.0d0,RALPHA,SS_DIM)
!     
!! COMPUTE INVERSE OF ABOBE RALPHA (=EYE+BETA*GAMMA) BY LAPACK ROUTINES AND STORE IT IN BETA (BETA=(EYE+BETA*GAMMA)**(-1))
!       CALL DGETRF(SS_DIM,SS_DIM,RALPHA,SS_DIM,IPIV,RCODE)
!  
!       IF(RCODE/=0) WRITE(IPT,*) 'error in computing LU factorization'
!       CALL DGETRI(SS_DIM,RALPHA,SS_DIM,IPIV,WORK,LWORK,RCODE)
!       IF(RCODE/=0) WRITE(IPT,*) 'error in inverting the matrix'
!       DO I=1, SS_DIM
!         DO J=1, SS_DIM
!           BETA2(I,J) = RALPHA(I,J)
!	 ENDDO
!       ENDDO
!     
!! READ RALPHA FROM FILE
!       FILENAME =  TRIM(OUTPUT_DIR)//'/rrktemp/Alpha.dat'    
!       OPEN(INORRK,FILE=TRIM(FILENAME),FORM='UNFORMATTED')
!       DO J=1, SS_DIM
!         READ(INORRK) (RALPHA(I,J), I=1, SS_DIM)    
!       ENDDO 
!       CLOSE(INORRK)
!
!! COMPUTE PRODUCT (EYE+BETA*GAMMA)**(-1)*RALPHA (GAMMA = BETA'*RALPHA)
!       CALL DGEMM('n','n',SS_DIM,SS_DIM,SS_DIM,1.0d0,BETA2,SS_DIM,RALPHA,SS_DIM,0.0d0,GAMMA,SS_DIM)  
!          
!! OUTPUT THIS PRODUCT TO TEMPORARY FILE
!       FILENAME = TRIM(OUTPUT_DIR)//'/rrktemp/temp.02'
!       OPEN(INORRK,FILE=TRIM(FILENAME),FORM='UNFORMATTED')
!       DO J=1, SS_DIM
!         WRITE(INORRK) (GAMMA(I,J), I=1, SS_DIM)    
!       ENDDO 
!       CLOSE(INORRK)     
!
!! COMPUTE PRODUCT RALPHA*(EYE+BETA*GAMMA)**(-1) AND STORE IN GAMMA (GAMMA = RALPHA*BETA')
!       CALL DGEMM('n','n',SS_DIM,SS_DIM,SS_DIM,1.0d0,RALPHA,SS_DIM,BETA2,SS_DIM,0.0d0,GAMMA,SS_DIM) 
!   
!! READ BACK OLD FILE
!       FILENAME = TRIM(OUTPUT_DIR)//'/rrktemp/Beta.dat'
!       OPEN(INORRK,FILE=TRIM(FILENAME),FORM='UNFORMATTED')
!       DO J=1, SS_DIM
!         READ(INORRK) (Beta2(I,J), I=1, SS_DIM)    
!       ENDDO 
!       CLOSE(INORRK)  
!
!! COMPUTE PRODUCT RALPHA*(EYE+BETA*GAMMA)**(-1)*BETA (RALPHA' = GAMMA'*BETA')
!       CALL DGEMM('n','n',SS_DIM,SS_DIM,SS_DIM,1.0d0,GAMMA,SS_DIM,BETA2,SS_DIM,0.0d0,RALPHA,SS_DIM)
!
!! READ BACK OLD RALPHA AND PUT IN GAMMA
!       FILENAME = TRIM(OUTPUT_DIR)//'/rrktemp/Alpha.dat'
!       OPEN(INORRK,FILE=TRIM(FILENAME),FORM='UNFORMATTED')
!       DO J=1, SS_DIM
!         READ(INORRK) (GAMMA(I,J), I=1, SS_DIM)    
!       ENDDO 
!       CLOSE(INORRK)  
!
!! COMPUTE BETA+RALPHA*(EYE+BETA*GAMMA)**(-1)*BETA*RALPHA_T (BETA = BETA+RALPHA'*GAMMA_T)
!       CALL DGEMM('n','t',SS_DIM,SS_DIM,SS_DIM,1.0d0,RALPHA,SS_DIM,GAMMA,SS_DIM,1.0d0,BETA2,SS_DIM)
!
!! MAKE SURE BETA IS SYMMETRIC
!       DO I=1, SS_DIM
!          DO J=1, I-1
!             BETA2(I,J) = 0.5_DP*(BETA2(I,J)+BETA2(J,I))
!             BETA2(J,I) = BETA2(I,J) 
!	  ENDDO
!       ENDDO
!
!! SAVE THIS NEW BETA TO FILE
!       FILENAME = TRIM(OUTPUT_DIR)//'/rrktemp/Beta.dat'
!       OPEN(INORRK,FILE=FILENAME, FORM='UNFORMATTED')
!         DO J=1, SS_DIM
!         WRITE(INORRK) (BETA2(I,J), I=1, SS_DIM)    
!       ENDDO 
!       CLOSE(INORRK)      
!
!! READ BACK OLD GAMMA
!       FILENAME = TRIM(OUTPUT_DIR)//'/rrktemp/Gamma.dat'
!       OPEN(INORRK,FILE=FILENAME, FORM='UNFORMATTED')
!       DO J=1, SS_DIM
!         READ(INORRK) (GAMMA(I,J), I=1, SS_DIM)    
!       ENDDO 
!       CLOSE(INORRK)     
!
!! READ TEMPORARY FILE INTO RALPHA
!       FILENAME = TRIM(OUTPUT_DIR)//'/rrktemp/temp.02'
!       OPEN(INORRK,FILE=FILENAME, FORM='UNFORMATTED')
!       DO J=1, SS_DIM
!         READ(INORRK) (RALPHA(I,J), I=1, SS_DIM)    
!       ENDDO 
!       CLOSE(INORRK)       
!
!! COMPUTE GAMMA*(EYE+BETA*GAMMA)**(-1)*RALPHA (BETA'=GAMMA*RALPHA')
!       CALL DGEMM('n','n',SS_DIM,SS_DIM,SS_DIM,1.0d0,GAMMA,SS_DIM,RALPHA,SS_DIM,0.0d0,BETA2,SS_DIM)
!
!! READ BACK OLD RALPHA
!       FILENAME = TRIM(OUTPUT_DIR)//'/rrktemp/Alpha.dat'
!       OPEN(INORRK,FILE=FILENAME, FORM='UNFORMATTED')
!       DO J=1, SS_DIM
!         READ(INORRK) (RALPHA(I,J), I=1, SS_DIM)    
!       ENDDO 
!       CLOSE(INORRK)     
!
!! COMPUTE GAMMA + RALPHA_T*GAMMA*(EYE+BETA*GAMMA)**(-1)*RALPHA (GAMMA = GAMMA+RALPHA_T*BETA')
!       CALL DGEMM('t','n',SS_DIM,SS_DIM,SS_DIM,1.0d0,RALPHA,SS_DIM,BETA2,SS_DIM,1.0d0,GAMMA,SS_DIM)
!
!! ENSURE GAMMA IS SYMMETRIC
!       DO I=1, SS_DIM
!          DO J=1, I-1
!             GAMMA(I,J) = 0.5_DP*(GAMMA(I,J)+GAMMA(J,I))
!             GAMMA(J,I) = GAMMA(I,J) 
!	  ENDDO
!       ENDDO
!
!! SAVE THIS NEW GAMMA TO FILE
!       FILENAME = TRIM(OUTPUT_DIR)//'/rrktemp/Gamma.dat'
!       OPEN(INORRK,FILE=FILENAME, FORM='UNFORMATTED')
!       DO J=1, SS_DIM
!         WRITE(INORRK) (GAMMA(I,J), I=1, SS_DIM)    
!       ENDDO 
!       CLOSE(INORRK)      
!       FILENAME = TRIM(OUTPUT_DIR)//'/rrktemp/Pfctr.txt'
!       OPEN(INORRK,FILE=FILENAME, FORM='FORMATTED')
!         DO J=1, SS_DIM
!           DO I=1, SS_DIM
!	      WRITE(INORRK,*) I,J,GAMMA(I,J)    
!           ENDDO
!	 ENDDO 
!       CLOSE(INORRK)
!
!! READ TEMPORARY FILE INTO GAMMA   
!       FILENAME = TRIM(OUTPUT_DIR)//'/rrktemp/temp.02'
!       OPEN(INORRK,FILE=FILENAME, FORM='UNFORMATTED')
!       DO J=1, SS_DIM
!         READ(INORRK) (GAMMA(I,J), I=1, SS_DIM)    
!       ENDDO 
!       CLOSE(INORRK)
!   
!! COMPUTE PRODUCT RALPHA*(EYE+BETA*GAMMA)**(-1)*RALPHA (BETA = RALPHA*GAMMA)  
!       CALL DGEMM('n','n',SS_DIM,SS_DIM,SS_DIM,1.0d0,RALPHA,SS_DIM,GAMMA,SS_DIM,0.0d0,BETA2,SS_DIM)
!   
!! SAVE THIS NEW RALPHA TO FILE
!       FILENAME = TRIM(OUTPUT_DIR)//'/rrktemp/Alpha.dat'
!       OPEN(INORRK,FILE=FILENAME, FORM='UNFORMATTED')
!         DO J=1, SS_DIM
!         WRITE(INORRK) (BETA2(I,J), I=1, SS_DIM)    
!       ENDDO 
!       CLOSE(INORRK)
!
!!===================================================================================================
!! END OF ITERATION FOR THE DOUBLING ALGORITHM          
!     ENDDO   
!!===================================================================================================
!
!     WRITE(IPT,*) 'Setting up the Kalman gain matrix in the reduced space......'   
!
!! READ IN PF FROM DOUBLING ALGORITHM 
!     FILENAME = TRIM(OUTPUT_DIR)//'/rrktemp/Gamma.dat'
!     OPEN(INORRK,FILE=FILENAME, FORM='UNFORMATTED')
!     DO J=1, SS_DIM
!       READ(INORRK) (GAMMA(I,J), I=1, SS_DIM)    
!     ENDDO 
!     CLOSE(INORRK)
!   
!! READ IN OBSERVATION OPERATOR   
!     FILENAME = TRIM(OUTPUT_DIR)//'/rrktemp/ObsOp.dat'
!     OPEN(INORRK,FILE=FILENAME, FORM='UNFORMATTED')
!     DO J=1, SS_DIM
!       READ(INORRK) (HEOF(I,J), I=1, NLOC)    
!     ENDDO 
!     CLOSE(INORRK)
!      
!! READ OBSERVATION ERROR COVARIANCE      
!     FILENAME = TRIM(OUTPUT_DIR)//'/rrktemp/Robs.dat'
!     OPEN(INORRK,FILE=FILENAME, FORM='UNFORMATTED')
!     DO J=1, NLOC
!       READ(INORRK) (R(I,J), I=1, NLOC)    
!     ENDDO 
!     CLOSE(INORRK)
!   
!! COMPUTE P*H^T (PHT = GAMMA * HEOF^T)   
!     CALL DGEMM('n','t',SS_DIM,NLOC,SS_DIM,1.0d0,GAMMA,SS_DIM,HEOF,ENKF_NOBSMAX,0.0d0,PHT,SS_DIM)
!
!! COMPUTE H*P*H^T (R2 = H*P*HT)  
!     CALL DGEMM('n','n',NLOC,NLOC,SS_DIM,1.0d0,HEOF,ENKF_NOBSMAX,PHT,SS_DIM,0.0d0,R2,ENKF_NOBSMAX)
!
!! OUPUT BOTH FORCAST AND OBSERVATION ERROR STANDARD DEVIATION IN ASCII
!     FILENAME =  TRIM(OUTPUT_DIR)//'/rrktemp/Fctobserr.txt'  
!     OPEN(INORRK,FILE=FILENAME, FORM='FORMATTED')
!     DO I=1, NLOC
!       WRITE(INORRK,*) STLOC(I),DSQRT(R2(I,I)),DSQRT(R(I,I))
!     ENDDO
!     CLOSE(INORRK)
!   
!! ADD R + (H P H^t) ---> R   
!     DO I = 1, NLOC
!        DO J=1, NLOC
!           R(J,I) = R(J,I) + R2(J,I)
!        ENDDO
!     ENDDO
!   
!! INVERT R (H*P*H^T+R)   
!     CALL DGETRF(NLOC,NLOC,R,ENKF_NOBSMAX,IPIV2,RCODE)
!     IF(RCODE/=0) WRITE(IPT,*) 'error in computing LU factorization'
!     CALL DGETRI(NLOC,R,ENKF_NOBSMAX,IPIV2,WORK2,LWORK2,RCODE)
!     IF(RCODE/=0) WRITE(IPT,*) 'error in inverting the matrix'
!     
!! COMPUTE KALMAN GAIN: K = P*HT*(H*P*H^T+R)^(-1)   
!     CALL DGEMM('n','n',SS_DIM,NLOC,NLOC,1.0d0,PHT,SS_DIM,R,ENKF_NOBSMAX,0.0d0,KAL,SS_DIM)
!   
!! OUTPUT RESULT   
!     OPEN(INORRK,FILE=KALNAME, FORM='UNFORMATTED') 
!     DO J=1, NLOC
!       WRITE(INORRK) (KAL(I,J), I=1, SS_DIM)    
!     ENDDO 
!     CLOSE(INORRK)     
!     FILENAME = TRIM(OUTPUT_DIR)//'/rrktemp/rrK.txt'
!     OPEN(INORRK,FILE=FILENAME, FORM='FORMATTED') 
!     DO J=1, NLOC
!       DO I=1, SS_DIM
!	 WRITE(INORRK,*) I,J,KAL(I,J)    
!       ENDDO
!     ENDDO 
!     CLOSE(INORRK)
!     
!! OUTPUT KALMAN GAIN MATRIX IN THE FULL SPACE: D*ER*KAL
!!     DO J=1, SS_DIM
!!       CALL READEOF(J,2)
!!       DO I=1, STDIM
!!         SEOF(I,J)=STTEMP1(I)
!!       ENDDO
!!     ENDDO
!     
!     DO J=1, NLOC
!       WRITE(STRCYC,'(I4.4)') J
!       DO I=1, STDIM
!          RRKSUM = 0.0_DP
!          DO K=1, SS_DIM
!             RRKSUM = RRKSUM + STSD(I)*SEOF(I,K)*KAL(K,J) 
!          ENDDO 
!	  STTEMP1(I) = RRKSUM
!       ENDDO
!       FILENAME = TRIM(OUTPUT_DIR)//'/rrktemp/'//'BigK'//STRCYC//'.txt'     
!       OPEN(INORRK,FILE=FILENAME,FORM='FORMATTED')
!       DO I=1, STDIM
!          WRITE(INORRK,*) STTEMP1(I)
!       ENDDO
!       CLOSE(INORRK)
!     ENDDO
!     
!! COMPUTE ANALYSIS ERROR COVARIANCE: P_AN1 = (I-K*H)*P (JUST AS DIAGNOSIS)
!     CALL DGEMM('n','t',SS_DIM,SS_DIM,NLOC,-1.0d0,KAL,SS_DIM,PHT,SS_DIM,1.0d0,GAMMA,SS_DIM)
!     FILENAME = TRIM(OUTPUT_DIR)//'/rrktemp/Panlr.txt'
!     OPEN(INORRK,FILE=FILENAME, FORM='FORMATTED') 
!     DO J=1, SS_DIM
!       DO I=1, SS_DIM
!	 WRITE(INORRK,*) I,J,GAMMA(I,J)    
!       ENDDO
!     ENDDO 
!     CLOSE(INORRK)
!     
!     CALL DGEMM('n','n',SS_DIM,SS_DIM,NLOC,1.0d0,KAL,SS_DIM,HEOF,ENKF_NOBSMAX,0.0d0,RALPHA,SS_DIM)
!     
!! COMPUTE (I-KH)M AND ITS SINGULAR VALUE     
!     DO I=1, SS_DIM
!       DO J=1, SS_DIM
!         IF(I/=J) THEN
!            RALPHA(I,J) = (-1.0_DP)*RALPHA(I,J)
!         ELSE
!	    RALPHA(I,I) = 1.0_DP - RALPHA(I,I)
!         ENDIF
!       ENDDO
!     ENDDO
!     
!     FILENAME = TRIM(OUTPUT_DIR)//'/rrktemp/I_KH.txt'
!     OPEN(INORRK,FILE=FILENAME, FORM='FORMATTED') 
!     DO J=1, SS_DIM
!       DO I=1, SS_DIM
!	 WRITE(INORRK,*) I,J,RALPHA(I,J)    
!       ENDDO
!     ENDDO 
!     CLOSE(INORRK)     
!     
!     FILENAME = TRIM(OUTPUT_DIR)//'/rrktemp/Amatr.dat'
!     OPEN(INORRK,FILE=FILENAME, FORM='UNFORMATTED')
!     DO J=1, SS_DIM
!       READ(INORRK) (TRANS(I,J), I=1, SS_DIM)    
!     ENDDO 
!     CLOSE(INORRK)     
!     CALL DGEMM('n','n',SS_DIM,SS_DIM,SS_DIM,1.0d0,RALPHA,SS_DIM,TRANS,SS_DIM,0.0d0,RALPHA,SS_DIM)
!    
!     FILENAME = TRIM(OUTPUT_DIR)//'/rrktemp/M_KHM.txt'
!     OPEN(INORRK,FILE=FILENAME, FORM='FORMATTED') 
!     DO J=1, SS_DIM
!       DO I=1, SS_DIM
!	 WRITE(INORRK,*) I,J,RALPHA(I,J)    
!       ENDDO
!     ENDDO 
!     CLOSE(INORRK) 
!     
!     DEALLOCATE(STLOC,IPIV,IPIV2,KAL,HEOF,R,HTR,R2,PHT,RALPHA,TRANS2)
!!     DEALLOCATE(BETA2,GAMMA,GAMMA2,WORK,WORK2,SEOF,STSD,HU,HUL,EVAL)
!     DEALLOCATE(BETA2,GAMMA,GAMMA2,WORK,WORK2,HU,HUL,EVAL)
!!     DEALLOCATE(STTEMP0,STTEMP1,STEOF,SDEOF,TRANS,RRKU,RRKV,RRKEL,RRKT,RRKS)
!   ENDIF
!   
!   RETURN
!   END SUBROUTINE RRK_RRK


! UTILITIES PROGRAMS
 
!  SUBROUTINE RRK_ALLOC_VAR
!
!   USE CONTROL
!   USE ENKFVAL
!   IMPLICIT NONE
!   
!   INTEGER STDIM
!   INTEGER SS_DIM
!   
!   STDIM = 0
!   IF(EL_ASSIM) STDIM = STDIM + MGL
!   IF(UV_ASSIM) STDIM = STDIM + 2*NGL*KBM1
!   IF(T_ASSIM)  STDIM = STDIM + MGL*KBM1
!   IF(S_ASSIM)  STDIM = STDIM + MGL*KBM1
!   
!   SS_DIM = ENKF_NENS
!   
!! ALLOCATE ARRYS
!
!   IF(.not.ALLOCATED(STTEMP0))  ALLOCATE(STTEMP0(STDIM))  ;STTEMP0   = ZERO
!   IF(.not.ALLOCATED(STTEMP1))  ALLOCATE(STTEMP1(STDIM))  ;STTEMP1   = ZERO
!   IF(.not.ALLOCATED(STEOF))    ALLOCATE(STEOF(STDIM))    ;STEOF     = ZERO
!   IF(.not.ALLOCATED(SDEOF))    ALLOCATE(SDEOF(STDIM))    ;SDEOF     = ZERO
!   IF(.not.ALLOCATED(TRANS))    ALLOCATE(TRANS(SS_DIM,SS_DIM)) ;TRANS     = ZERO
!   IF(.not.ALLOCATED(RRKEL))    ALLOCATE(RRKEL(0:MGL))    ;RRKEL     = ZERO
!   IF(.not.ALLOCATED(RRKU))     ALLOCATE(RRKU(0:NGL,KB))  ;RRKU      = ZERO
!   IF(.not.ALLOCATED(RRKV))     ALLOCATE(RRKV(0:NGL,KB))  ;RRKV      = ZERO
!   IF(.not.ALLOCATED(RRKT))     ALLOCATE(RRKT(0:MGL,KB))  ;RRKT      = ZERO
!   IF(.not.ALLOCATED(RRKS))     ALLOCATE(RRKS(0:MGL,KB))  ;RRKS      = ZERO
!   
!! END ALLOCATION
!   
!   RETURN
!  END SUBROUTINE RRK_ALLOC_VAR

!------------------------------------------------------------------------------|
!  PERTURB THE BASE STATE IN THE EIGENVECTOR DIRECTIONS                        |
!------------------------------------------------------------------------------|   
   
!  SUBROUTINE PERTURB
!
!#  if defined(WET_DRY)
!     USE MOD_WD
!#  endif
!   USE LIMS
!   USE ALL_VARS
!   USE ENKFVAL
!   USE CONTROL
!#  if defined (MULTIPROCESSOR)
!     USE MOD_PAR
!#  endif
!   IMPLICIT NONE
!    
!    INTEGER SS_DIM
!    INTEGER STDIM
!    INTEGER I_EOF
!    INTEGER II,I,J
!    INTEGER IDUMMY
!    INTEGER I_START
!    CHARACTER(LEN=80) IFILE
!    CHARACTER(LEN=80) IFILE2
!    CHARACTER(LEN=80) TEMP1FILE
!    CHARACTER(LEN=80) TEMP2FILE
!    CHARACTER(LEN=4)  FEOF
!      
!    STDIM = 0
!    IF(EL_ASSIM) STDIM = STDIM + MGL
!    IF(UV_ASSIM) STDIM = STDIM + 2*NGL*KBM1
!    IF(T_ASSIM)  STDIM = STDIM + MGL*KBM1
!    IF(S_ASSIM)  STDIM = STDIM + MGL*KBM1
!    
!    SS_DIM = ENKF_NENS
!
!   ALLOCATE(STTEMP0(STDIM))                ;STTEMP0   = ZERO
!   ALLOCATE(STTEMP1(STDIM))                ;STTEMP1   = ZERO
!   ALLOCATE(STEOF(STDIM))                  ;STEOF     = ZERO
!   ALLOCATE(SDEOF(STDIM))                  ;SDEOF     = ZERO
!   ALLOCATE(TRANS(SS_DIM,SS_DIM))          ;TRANS     = ZERO        
!
!    IDUMMY  = 0
!    IF(EL_ASSIM) THEN
!      DO I=1, MGL
!        IDUMMY = IDUMMY + 1
!        SDEOF(IDUMMY) = ELSD
!      ENDDO
!    ENDIF
!    IF(UV_ASSIM) THEN
!      DO I=1, 2*KBM1
!        DO J=1, NGL 
!          IDUMMY = IDUMMY + 1
!          SDEOF(IDUMMY)=USD*DSQRT(DBLE(KBM1))     
!        ENDDO
!      ENDDO
!    ENDIF
!    IF(T_ASSIM) THEN
!      DO I=1, KBM1
!        DO J=1, MGL 
!          IDUMMY = IDUMMY + 1
!          SDEOF(IDUMMY)=TSD
!        ENDDO
!      ENDDO    
!    ENDIF
!    IF(S_ASSIM) THEN
!      DO I=1, KBM1
!        DO J=1, MGL 
!          IDUMMY = IDUMMY + 1
!          SDEOF(IDUMMY)=SSD
!        ENDDO
!      ENDDO    
!    ENDIF
!    
!    print *, 'before rrk_rrk read restart, file name:',nc_start%fname ! should be marked    
!    CALL READRESTART(NC_START,REF_START_TIME)
!    CALL GR2ST(0)
!
!    
!! READ THE EOF AND PERTURB THE BASE STATE IN THIS DIRECTION
!
!   DO I_EOF=1, SS_DIM
!!       CALL READEOF(I_EOF,1)
!       DO I=1, STDIM  
!          STEOF(I) = SEOF(I,I_EOF)
!       ENDDO
!
!! PERTURB THE BASE STATE IN THE DIRECTION OF THE I_EOF'TH EOF AND STORE IT IN THE MP1FILE (FOR RESTART)
!       DO I=1, STDIM
!          STTEMP1(I) = STTEMP0(I) + STEOF(I)*SDEOF(I)*DBLE(RRK_PSIZE)
!       ENDDO
!       CALL ST2GR
!
!
!   IF(SERIAL) THEN
!     EL =  RRKEL 
!     U  =  RRKU
!     V  =  RRKV
!     T1 =  RRKT
!     S1 =  RRKS
!   ELSE 
!#  if defined (MULTIPROCESSOR)
!     DO I=1, M
!        EL(I)=RRKEL(NGID(I))
!     ENDDO
!
!     DO I=1, N
!       DO J=1, KB
!         U(I,J)=RRKU(EGID(I),J)
!         V(I,J)=RRKV(EGID(I),J)
!       ENDDO
!     ENDDO
!
!     DO I=1, M
!       DO J=1, KB
!         T1(I,J)=RRKT(NGID(I),J)
!         S1(I,J)=RRKS(NGID(I),J)
!       ENDDO
!     ENDDO
!#  endif
!   ENDIF
!
!
!! WRITE THE PERTURBED STATE IN TEMP1FILE
!    print *, 'where to read to perturb:', nc_start%fname
!    print *, 'where to write perturb:', NCF_ENKFRE%fname
!
!     CALL WRITERESTART(NCF_ENKFRE,I_EOF)
!
!   ENDDO
!
!
!  RETURN
!  END SUBROUTINE PERTURB

!------------------------------------------------------------------------------|
!  CALCULATE THE LINEARIZED MODEL IN THE REDUCED SPACE                         |
!------------------------------------------------------------------------------|   
   
!  SUBROUTINE MREDUCED
!
!   USE MOD_WD
!   USE LIMS
!   USE CONTROL
!   USE ENKFVAL
!   IMPLICIT NONE
!    
!   INTEGER SS_DIM
!   INTEGER STDIM
!   INTEGER I_EOF
!   INTEGER I,J
!   INTEGER IDUMMY
!   INTEGER I_START   
!   CHARACTER(LEN=80) IFILE
!   CHARACTER(LEN=80) IFILE2
!   CHARACTER(LEN=80) TEMPFILE
!   CHARACTER(LEN=4)  FEOF
!   REAL(DP) ::  SUM0
!
!   STDIM = 0
!   IF(EL_ASSIM) STDIM = STDIM + MGL
!   IF(UV_ASSIM) STDIM = STDIM + 2*NGL*KBM1
!   IF(T_ASSIM)  STDIM = STDIM + MGL*KBM1
!   IF(S_ASSIM)  STDIM = STDIM + MGL*KBM1
!
!   SS_DIM = ENKF_NENS
!   
!   IDUMMY = 0
!   IF(EL_ASSIM) THEN
!     DO I=1, MGL
!       IDUMMY = IDUMMY + 1
!       SDEOF(IDUMMY) = ELSD
!     ENDDO
!   ENDIF
!   IF(UV_ASSIM) THEN
!     DO I =1, 2*KBM1
!       DO J=1, NGL
!         IDUMMY = IDUMMY + 1
!         IF(RRK_OPTION == 0) THEN
!            SDEOF(IDUMMY) = USD
!         ELSE
!            SDEOF(IDUMMY) = USD*DSQRT(DBLE(KBM1))
!         ENDIF
!       ENDDO
!     ENDDO
!   ENDIF
!   IF(T_ASSIM) THEN
!     DO I =1, KBM1
!       DO J=1, MGL
!         IDUMMY = IDUMMY + 1
!         SDEOF(IDUMMY) = TSD
!       ENDDO
!     ENDDO
!   ENDIF
!   IF(S_ASSIM) THEN
!     DO I =1, KBM1
!       DO J=1, MGL
!         IDUMMY = IDUMMY + 1
!         SDEOF(IDUMMY) = SSD
!       ENDDO
!     ENDDO
!   ENDIF   
!   
!   CALL READRESTART(NC_START,REF_START_TIME)
!   CALL GR2ST(0)  
!   
!   
!! READ THE FORECAST WHICH IS PROPAGATED FROM THE PERTURBED STATE  
!
!   CALL NC_OPEN(NCF_ENKFfct)
! 
!   DO I_EOF =1, SS_DIM   
!! READ THE EACH PERTURBED STATE AT THE END OF LINEARUZATION TIME STEP
!
!   CALL READRESTART(NCF_ENKFfct,I_EOF) 
!
!   CALL GR2ST(1)   
!   
!! PROJECT  EVOLVED PERTURBATION ONTO EOFs
!     DO I=1, SS_DIM
!       DO J=1, STDIM  
!          STEOF(J) = SEOF(J,I)
!       ENDDO
!       SUM0 = 0.0_DP
!       DO J=1, STDIM
!         SUM0 = SUM0 + STEOF(J)/SDEOF(J)*(STTEMP1(J)-STTEMP0(J))
!       ENDDO
!       TRANS(I,I_EOF) = SUM0/DBLE(RRK_PSIZE)
!           
!     ENDDO
!! EDN OF LOOP OVER EACH EOF   
!   ENDDO
!   
!  END SUBROUTINE MREDUCED 


!=====================================================================================/
!  CONVERT THE STATE FROM THE GRID SPACE TO THE STATE VECTOR SPACE                    /
!=====================================================================================/
!  SUBROUTINE GR2ST(OPT)
! 
!   USE LIMS
!   USE ALL_VARS
!   USE ENKFVAL
!   IMPLICIT NONE
!   
!    INTEGER IDUMMY, I, J
!    INTEGER OPT
!
!    IDUMMY=0
!    
!    IF(EL_ASSIM) THEN
!      DO I=1, MGL
!        IDUMMY = IDUMMY + 1
!        IF (OPT == 0) THEN
!          STTEMP0(IDUMMY) = EL(I)
!        ELSE
!          STTEMP1(IDUMMY) = EL(I)
!        ENDIF	 	
!      ENDDO
!    ENDIF
!
!    IF(UV_ASSIM) THEN
!      DO I=1, KBM1
!        DO J=1, NGL
!          IDUMMY = IDUMMY + 1
!          IF (OPT == 0) THEN
!	    STTEMP0(IDUMMY) = U(J,I)
!
!	  ELSE
!	    STTEMP1(IDUMMY) = U(J,I)
!	  ENDIF    
!        ENDDO
!      ENDDO
!      DO I=1, KBM1
!        DO J=1, NGL
!          IDUMMY = IDUMMY + 1
!          IF (OPT == 0) THEN
!	    STTEMP0(IDUMMY) = V(J,I)
!	  ELSE
!	    STTEMP1(IDUMMY) = V(J,I)
!	  ENDIF    
!        ENDDO
!      ENDDO		
!    ENDIF  
!  
!    IF(T_ASSIM) THEN
!      DO I=1, KBM1
!        DO J=1, MGL
!          IDUMMY = IDUMMY + 1
!          IF (OPT == 0) THEN
!	    STTEMP0(IDUMMY) = T(J,I)
!	  ELSE
!	    STTEMP1(IDUMMY) = T(J,I)
!	  ENDIF    
!        ENDDO
!      ENDDO
!    ENDIF
!  
!    IF(S_ASSIM) THEN
!      DO I=1, KBM1
!        DO J=1, MGL
!          IDUMMY = IDUMMY + 1
!          IF (OPT == 0) THEN
!	    STTEMP0(IDUMMY) = S(J,I)
!	  ELSE
!	    STTEMP1(IDUMMY) = S(J,I)
!	  ENDIF    
!        ENDDO
!      ENDDO
!    ENDIF
!
!  RETURN
!  END SUBROUTINE GR2ST

!=====================================================================================/
!  CONVERT THE STATE FROM THE STATE VECTOR SPACE TO THE GRID SPACE                    /
!=====================================================================================/
!  SUBROUTINE ST2GR
! 
!   USE LIMS
!   USE ALL_VARS
!   USE ENKFVAL
!   IMPLICIT NONE
!    
!    INTEGER IDUMMY, I, J
!
!    IDUMMY=0
!    IF(EL_ASSIM) THEN
!      DO I=1, MGL
!        IDUMMY = IDUMMY + 1
!        RRKEL(I) = STTEMP1(IDUMMY) 
!      ENDDO
!    ENDIF
!    IF(UV_ASSIM) THEN
!      DO I=1, KBM1
!        DO J=1, NGL
!          IDUMMY = IDUMMY + 1
!          RRKU(J,I) = STTEMP1(IDUMMY) 
!        ENDDO
!      ENDDO
!      DO I=1, KBM1
!        DO J=1, NGL
!          IDUMMY = IDUMMY + 1
!          RRKV(J,I) = STTEMP1(IDUMMY)
!        ENDDO
!      ENDDO
!    ENDIF
!    IF(T_ASSIM) THEN
!      DO I=1, KBM1
!        DO J=1, MGL
!          IDUMMY = IDUMMY + 1
!          RRKT(J,I) = STTEMP1(IDUMMY) 
!        ENDDO
!      ENDDO
!    ENDIF
!    IF(S_ASSIM) THEN
!      DO I=1, KBM1
!        DO J=1, MGL
!          IDUMMY = IDUMMY + 1
!          RRKS(J,I) = STTEMP1(IDUMMY) 
!        ENDDO
!      ENDDO
!    ENDIF
!    
!  RETURN
!  END SUBROUTINE ST2GR


!=====================================================================================/
!  DETERMINE IF NODES/ELEMENTS ARE WET OR DRY                                         /
!=====================================================================================/
!  SUBROUTINE WET_JUDGE_EL
!
!   USE MOD_PREC
!   USE ALL_VARS
!#  if defined (MULTIPROCESSOR)
!   USE MOD_PAR
!#  endif
!#  if defined(WET_DRY)
!   USE MOD_WD
!#  endif
!   IMPLICIT NONE
!   REAL(SP) :: DTMP
!   INTEGER  :: ITA_TEMP
!   INTEGER  :: I,IL,IA,IB,K1,K2,K3,K4,K5,K6
!
!#  if defined(WET_DRY)
!
!!
!!--Determine If Node Points Are Wet/Dry Based on Depth Threshold---------------!
!!
!   ISWETN = 1
!   DO I = 1, M
!     DTMP = H(I) + EL(I)
!     IF((DTMP - MIN_DEPTH) < 1.0E-5_SP) ISWETN(I) = 0
!   END DO
!
!!
!!--Determine if Cells are Wet/Dry Based on Depth Threshold---------------------!
!!
!   ISWETC = 1
!   DO I = 1, N
!     DTMP =  MAX(EL(NV(I,1)),EL(NV(I,2)),EL(NV(I,3)))  + &
!             MIN(  H(NV(I,1)),  H(NV(I,2)),  H(NV(I,3)))
!     IF((DTMP - MIN_DEPTH) < 1.0E-5_SP) ISWETC(I) = 0
!   END DO
!
!!
!!--A Secondary Condition for Nodal Dryness-(All Elements Around Node Are Dry)--!
!!
!   DO I = 1, M
!     IF(SUM(ISWETC(NBVE(I,1:NTVE(I)))) == 0)  ISWETN(I) = 0
!   END DO
!
!!
!!--Adjust Water Surface So It Does Not Go Below Minimum Depth------------------!
!!
!   EL = MAX(EL,-H + MIN_DEPTH)
!
!!
!!--Recompute Element Based Depths----------------------------------------------!
!!
!   DO I = 1, N
!     EL1(I) = ONE_THIRD*(EL(NV(I,1))+EL(NV(I,2))+EL(NV(I,3)))
!   END DO
!
!!
!!--Extend Element/Node Based Wet/Dry Flags to Domain Halo----------------------!
!!
!#  if defined (MULTIPROCESSOR)
!   IF(PAR)THEN
!
!     CALL AEXCHANGE(EC,MYID,NPROCS,ISWETC)
!     CALL AEXCHANGE(NC,MYID,NPROCS,ISWETN)
!
!   END IF
!
!#  endif
!#  endif
!
!   RETURN
!
!  END SUBROUTINE WET_JUDGE_EL

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
      INTEGER STLOC(ENKF_NOBSMAX)
      INTEGER LAY(ENKF_NOBSMAX)
 
      FILENAME = TRIM(INPUT_DIR)//"/"//trim(casename)//"_assim_enkf.nml"
     
      OPEN(73,FILE=TRIM(FILENAME),FORM='FORMATTED')
      REWIND(73)
 
      NLOC = 0
      IDUMMY = 0
      STLOC = ZERO
      
  100 READ(73,'(A24)',END=200) HEADINFO
      IF(SWITCH/=1) THEN
        IF(HEADINFO=='!===READ IN OBSERVATION') THEN
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
          IF(NLOC>ENKF_NOBSMAX) THEN
            WRITE(IPT,*) 'not enough storage for observations:', 'Nloc=', Nloc, 'Nobsmax=', ENKF_NOBSMAX
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
          IF(NLOC+NUM>ENKF_NOBSMAX) THEN
            WRITE(IPT,*) 'not enough storage for observations:', 'Nloc=', Nloc+num, 'Nobsmax=', ENKF_NOBSMAX
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
          IF(NLOC>ENKF_NOBSMAX) THEN
            WRITE(IPT,*) 'not enough storage for observations:', 'Nloc=', Nloc, 'Nobsmax=', ENKF_NOBSMAX
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
          PRINT *, 'SALINITY READING OBSERVATION',NUM
 	 NLOC = NLOC + NUM
          IF(NLOC>ENKF_NOBSMAX) THEN
            WRITE(IPT,*) 'not enough storage for observations:', 'Nloc=', Nloc, 'Nobsmax=', ENKF_NOBSMAX
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
!-------------old--------------------------------
  SUBROUTINE ENKF_READRESTART_first(NCF)
    USE MOD_INPUT
    USE MOD_STARTUP
    IMPLICIT NONE
    TYPE(NCFILE), POINTER :: NCF, NCF_tmP
    INTEGER :: STATE
    NCF_tmP=>NC_START
    IF(.not. associated(NCF)) CALL FATAL_ERROR&
         & ("MOD_ENKF: READRESTART: THE FILE POINTER PASSED IS NOT ASSOCIATED?")
    NCF%FTIME%PREV_STKCNT = 1
    NC_START=>NCF


    print *, 'read ssh enkf'
    CALL READ_SSH
    print *, 'finish read ssh enkf'

    IF(WETTING_DRYING_ON) CALL READ_WETDRY

    CALL READ_UV

    CALL READ_TURB

    CALL READ_TS

    NC_START=>NCF_tmP

  END SUBROUTINE ENKF_READRESTART_first
!-------------------new-------------------------------------------------------
!-------------old--------------------------------
  SUBROUTINE ENKF_READRESTART_last(NCF)
    USE MOD_INPUT
    USE MOD_STARTUP
    IMPLICIT NONE
    TYPE(NCFILE), POINTER :: NCF, NCF_tmP
    INTEGER :: STATE
   TYPE(NCDIM), POINTER :: DIM
   TYPE(NCVAR), POINTER :: VAR
   LOGICAL :: FOUND
   INTEGER :: NTIMES
   ! GET THE TIME DIMENSION
   DIM => FIND_UNLIMITED(NCF,FOUND)
   IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN ENKF_READ_FIELD FILE",&
            & "FILE NAME: "//TRIM(NCF%FNAME),&
            &"COULD NOT FIND THE UNLIMITED DIMENSION")
   NTIMES = DIM%DIM

    NCF_tmP=>NC_START
    IF(.not. associated(NCF)) CALL FATAL_ERROR&
         & ("MOD_ENKF: READRESTART: THE FILE POINTER PASSED IS NOT ASSOCIATED?")
    NCF%FTIME%PREV_STKCNT = NTIMES
    NC_START=>NCF


    print *, 'read ssh enkf'
    CALL READ_SSH
    print *, 'finish read ssh enkf'

    IF(WETTING_DRYING_ON) CALL READ_WETDRY

    CALL READ_UV

    CALL READ_TURB

    CALL READ_TS

    NC_START=>NCF_tmP

  END SUBROUTINE ENKF_READRESTART_last
!-------------------new-------------------------------------------------------
  SUBROUTINE ENKF_READRESTART(NCF)
   USE CONTROL
   USE ENKFVAL
   USE MOD_NCTOOLS
    USE MOD_INPUT
    USE MOD_STARTUP
    IMPLICIT NONE
    TYPE(NCFILE), POINTER :: NCF, NCF_tmP
    INTEGER :: STATE
   TYPE(TIME) :: tTEST
   LOGICAL :: FOUND
   INTEGER :: NTIMES
   

   TYPE(NCDIM), POINTER :: DIM
   TYPE(NCVAR), POINTER :: VAR

   INTEGER :: nrange,IRANGE(2),CNT, STATUS
   integer :: istkcnt
   character(len=120) :: pathnfile 
   
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
      IF (INTTIME == TTEST) THEN
         IRANGE(1) = CNT
         istkcnt=cnt
      END IF

      IF (inttime == TTEST) THEN
         IRANGE(2) = CNT
      END IF

   END DO


   nrange = irange(2)-irange(1)+1
!------------------FINISH FIND CNT---------------------
    NCF_tmP=>NC_START
    IF(.not. associated(NCF)) CALL FATAL_ERROR&
         & ("MOD_ENKF: READRESTART: THE FILE POINTER PASSED IS NOT ASSOCIATED?")
    NCF%FTIME%PREV_STKCNT = CNT
    NC_START=>NCF


    print *, 'read ssh enkf'
    CALL READ_SSH
    print *, 'finish read ssh enkf'

    IF(WETTING_DRYING_ON) CALL READ_WETDRY

    CALL READ_UV

    CALL READ_TURB

    CALL READ_TS

    NC_START=>NCF_tmP

  END SUBROUTINE ENKF_READRESTART











 
!  SUBROUTINE ENKF_READRESTART(NCF)
!    USE MOD_INPUT
!    USE MOD_STARTUP_mimic
!    IMPLICIT NONE
!    TYPE(NCFILE), POINTER :: NCF
!    INTEGER :: STATE

!    IF(.not. associated(NCF)) CALL FATAL_ERROR&
!         & ("MOD_ENKF: READRESTART: THE FILE POINTER PASSED IS NOT ASSOCIATED?")

!    NCF%FTIME%PREV_STKCNT = 1



!    CALL READ_SSH1(NCF)

!    IF(WETTING_DRYING_ON) CALL READ_WETDRY1(NCF)

!    CALL READ_UV1(NCF)

!    CALL READ_TURB1(NCF)

!    CALL READ_TS1(NCF)

!  END SUBROUTINE ENKF_READRESTART

!  SUBROUTINE READRESTART_STATE(NCF,STATE)
!    USE MOD_INPUT
!    USE MOD_STARTUP_mimic
!    IMPLICIT NONE
!    TYPE(NCFILE), POINTER :: NCF
!    INTEGER :: STATE
!    
!    IF(.not. associated(NCF)) CALL FATAL_ERROR&
!         & ("MOD_ENKF: READRESTART: THE FILE POINTER PASSED IS NOT ASSOCIATED?")
!    
!    NCF%FTIME%PREV_STKCNT = STATE
!
!    
!
!    CALL READ_SSH1(NCF)
!
!    IF(WETTING_DRYING_ON) CALL READ_WETDRY1(NCF)
!
!    CALL READ_UV1(NCF)
!
!    CALL READ_TURB1(NCF)
!
!    CALL READ_TS1(NCF)
!
!  END SUBROUTINE READRESTART_STATE
!
!  SUBROUTINE READRESTART_TIME(NCF,NOW)
!    USE MOD_INPUT
!    USE MOD_STARTUP_MIMIC
!    IMPLICIT NONE
!    TYPE(NCFILE), POINTER :: NCF
!    TYPE(TIME) :: NOW,TTEST
!    INTEGER :: STATE,NSTATE
!    TYPE(NCDIM), POINTER :: DIM
!    LOGICAL FOUND
!    
!    IF(.not. associated(NCF)) CALL FATAL_ERROR&
!         & ("MOD_ENKF: READRESTART: THE FILE POINTER PASSED IS NOT ASSOCIATED?")
!
!    ! GET THE TIME DIMENSION
!   DIM => FIND_UNLIMITED(NCF,FOUND)
!   IF(.not. FOUND) CALL FATAL_ERROR &
!            & ("IN RESTART FILE",&
!            & "FILE NAME: "//TRIM(NCF%FNAME),&
!            &"COULD NOT FIND THE UNLIMITED DIMENSION")
!   NSTATE = DIM%DIM
!   FOUND = .false.
!   DO STATE=1,NSTATE
!      
!      TTEST = get_file_time(NCF,STATE)
!      IF(TTEST == NOW) THEN
!         FOUND = .true.
!         exit
!      END IF
!      
!   END DO
!   
!   IF(FOUND)THEN
!      CALL READRESTART_STATE(NCF,STATE)
!   ELSE
!      CALL PRINT_TIME(NOW,IPT,"TIME REQUSTED")
!      CALL FATAL_ERROR&
!           &("COULD NOT FIND CORRECT TIME IN RESTART FILE",&
!           &"FILE NAME: "//TRIM(NCF%FNAME))
! 
!  END IF
!
!  END SUBROUTINE READRESTART_TIME

!  SUBROUTINE WRITERESTART(NCF,STATE)
!    IMPLICIT NONE
!    TYPE(NCFILE), POINTER :: NCF
!    INTEGER :: STATE
!
!    IF(.not. associated(NCF)) CALL FATAL_ERROR&
!         & ("MOD_ENKF: WRITERESTART: THE FILE POINTER PASSED IS NOT ASSOCIATED?")
!
!    NCF%FTIME%NEXT_STKCNT = state
!     CALL UPDATE_IODATA(NCF,IntTime)
!    CALL NC_WRITE_FILE(NCF,LOCAL_DISK)
!    
!  END SUBROUTINE WRITERESTART
!-----------------old-----------------------
!  SUBROUTINE ENKF_WRITERESTART(NCF)
!    IMPLICIT NONE
!    TYPE(NCFILE), POINTER :: NCF
!    INTEGER :: STATE
!
!    IF(.not. associated(NCF)) CALL FATAL_ERROR&
!         & ("MOD_ENKF: WRITERESTART: THE FILE POINTER PASSED IS NOT ASSOCIATED?")
!
!!    NCF%FTIME%NEXT_STKCNT = state
!    NCF%FTIME%NEXT_STKCNT = 1
!     CALL UPDATE_IODATA(NCF,IntTime)
!    CALL NC_WRITE_FILE(NCF,LOCAL_DISK)
!
!  END SUBROUTINE ENKF_WRITERESTART
!-------------new----------------------------------
  SUBROUTINE ENKF_WRITERESTART(NCF)
    IMPLICIT NONE
    TYPE(NCFILE), POINTER :: NCF
    INTEGER :: STATE
        TYPE(NCFTIME), POINTER :: FTM

    IF(.not. associated(NCF)) CALL FATAL_ERROR&
         & ("MOD_ENKF: WRITERESTART: THE FILE POINTER PASSED IS NOT ASSOCIATED?")

!    NCF%FTIME%NEXT_STKCNT = state
!    NCF%FTIME%NEXT_STKCNT = 1
!     CALL UPDATE_IODATA(NCF,IntTime)
!    CALL NC_WRITE_FILE(NCF,LOCAL_DISK)
    FTM => NCF%FTIME
    CALL UPDATE_IODATA(NCF,IntTime)

    CALL NC_WRITE_FILE(NCF)
    IF (ENKF_memberCONTR == ENKF_NENS) then
    ! ONCE THE FILE IS WRITEN INCRIMENT THE FILE OBJECT TIME
    FTM%PREV_IO = IntTime
    FTM%NEXT_IO = FTM%NEXT_IO + ENKF_INTERVAL

    ! INCRIMENT THE STACK COUNT
    FTM%PREV_STKCNT = FTM%NEXT_STKCNT
    FTM%NEXT_STKCNT = FTM%NEXT_STKCNT + 1
    END IF
  END SUBROUTINE ENKF_WRITERESTART
!=======================================================|

  SUBROUTINE ENKF_SETUP_MEAN_NC(NCF)
  IMPLICIT NONE
  INTEGER :: STKCNT
    TYPE(NCFILE), POINTER :: NCF
    TYPE(NCFTIME), POINTER :: FTM 
       IF(.not. associated(NCF)) CALL FATAL_ERROR&
         & ("MOD_ENKF: WRITERESTART: THE FILE POINTER PASSED IS NOT ASSOCIATED?")
    if (ARCHIVE_TIMES==0) then
     enkf_MEAN_ncfout => COPY_FILE(NCF)
     enkf_MEAN_ncfout%FNAME=TRIM(OUTPUT_DIR)//"enkftemp/enkf_forecast_mean.nc"
     enkf_MEAN_ncfout%FTIME%MAX_STKCNT =0 
    FTM => enkf_mean_ncfout%FTIME

    CALL NC_WRITE_FILE(enkf_mean_ncfout,stkcnt=0) ! because this program is called after the call archive , it FTM%NEXT_IO already became 1, then I don't need below
!    IF (ENKF_memberCONTR == ENKF_NENS) then
!    ! ONCE THE FILE IS WRITEN INCRIMENT THE FILE OBJECT TIME
!
!
!    FTM%PREV_IO = IntTime
!    FTM%NEXT_IO = FTM%NEXT_IO + FTM%INTERVAL
!
!    ! INCRIMENT THE STACK COUNT
!    FTM%PREV_STKCNT = FTM%NEXT_STKCNT
!    FTM%NEXT_STKCNT = FTM%NEXT_STKCNT + 1

!    END IF

    WRITE(IPT,*) "END DUMP mean DATA kf_times =0"


     enkf_MEAN_ncfout => COPY_FILE(NCF)
     enkf_MEAN_ncfout%FNAME=TRIM(OUTPUT_DIR)//"enkftemp/enkf_analysis_mean.nc"
     enkf_MEAN_ncfout%FTIME%MAX_STKCNT =0
    FTM => enkf_mean_ncfout%FTIME

    CALL NC_WRITE_FILE(enkf_mean_ncfout,stkcnt=0) ! because this program is called after the call archive , it FTM%NEXT_IO already became 1, then I don't need below
!    IF (ENKF_memberCONTR == ENKF_NENS) then
!    ! ONCE THE FILE IS WRITEN INCRIMENT THE FILE OBJECT TIME
!
!
!    FTM%PREV_IO = IntTime
!    FTM%NEXT_IO = FTM%NEXT_IO + FTM%INTERVAL
!
!    ! INCRIMENT THE STACK COUNT
!    FTM%PREV_STKCNT = FTM%NEXT_STKCNT
!    FTM%NEXT_STKCNT = FTM%NEXT_STKCNT + 1

!    END IF

    WRITE(IPT,*) "END DUMP mean DATA kf_times =0"



! acturall I didn't use call below

!    else

!    FTM => enkf_mean_ncfout%FTIME
!    CALL UPDATE_IODATA(enkf_mean_ncfout,IntTime)
!
!    CALL NC_WRITE_FILE(enkf_mean_ncfout)
!    IF (ENKF_memberCONTR == ENKF_NENS) then
!    ! ONCE THE FILE IS WRITEN INCRIMENT THE FILE OBJECT TIME
!    FTM%PREV_IO = IntTime
!    FTM%NEXT_IO = FTM%NEXT_IO + FTM%INTERVAL
!
!    ! INCRIMENT THE STACK COUNT
!    FTM%PREV_STKCNT = FTM%NEXT_STKCNT
!    FTM%NEXT_STKCNT = FTM%NEXT_STKCNT + 1
!    END IF
!    WRITE(IPT,*) "END DUMP mean DATA kftimes > 0, FTM%NEXT_STKCNT",FTM%NEXT_STKCNT


    end if


  END SUBROUTINE ENKF_SETUP_MEAN_NC
!=======================================================|
!=======================================================|

  SUBROUTINE ENKF_SETUP_NC(NCF)
  IMPLICIT NONE
  INTEGER :: STKCNT
    TYPE(NCFILE), POINTER :: NCF
    TYPE(NCFTIME), POINTER :: FTM 
       IF(.not. associated(NCF)) CALL FATAL_ERROR&
         & ("MOD_ENKF: WRITERESTART: THE FILE POINTER PASSED IS NOT ASSOCIATED?")
    if (ARCHIVE_TIMES==0) then
     enkf_ncfout => COPY_FILE(NCF)
     enkf_ncfout%FNAME=TRIM(OUTPUT_DIR)//"enkftemp/enkf_output_"//enkf_num//".nc"
     enkf_ncfout%FTIME%MAX_STKCNT =0 
    FTM => enkf_ncfout%FTIME
!    CALL UPDATE_IODATA(enkf_ncfout,IntTime)

    CALL NC_WRITE_FILE(enkf_ncfout,stkcnt=0) ! because this program is called after the call archive , it FTM%NEXT_IO already became 1, then I don't need below
    IF (ENKF_memberCONTR == ENKF_NENS) then
!    ! ONCE THE FILE IS WRITEN INCRIMENT THE FILE OBJECT TIME
!
!
!    FTM%PREV_IO = IntTime
!    FTM%NEXT_IO = FTM%NEXT_IO + FTM%INTERVAL
!
!    ! INCRIMENT THE STACK COUNT
!    FTM%PREV_STKCNT = FTM%NEXT_STKCNT
!    FTM%NEXT_STKCNT = FTM%NEXT_STKCNT + 1

    END IF

    WRITE(IPT,*) "END DUMP_DATA kf_times =0"




    else
    enkf_ncfout%FNAME=TRIM(OUTPUT_DIR)//"enkftemp/enkf_output_"//enkf_num//".nc"
    FTM => enkf_ncfout%FTIME
    CALL UPDATE_IODATA(enkf_ncfout,IntTime)

    CALL NC_WRITE_FILE(enkf_ncfout)
    IF (ENKF_memberCONTR == ENKF_NENS) then
    ! ONCE THE FILE IS WRITEN INCRIMENT THE FILE OBJECT TIME
    FTM%PREV_IO = IntTime
    FTM%NEXT_IO = FTM%NEXT_IO + FTM%INTERVAL

    ! INCRIMENT THE STACK COUNT
    FTM%PREV_STKCNT = FTM%NEXT_STKCNT
    FTM%NEXT_STKCNT = FTM%NEXT_STKCNT + 1
    END IF
    WRITE(IPT,*) "END DUMP_DATA kftimes > 0, FTM%NEXT_STKCNT",FTM%NEXT_STKCNT

!    call NC_WRITE_FILE(enkf_ncfout,stkcnt=STKCNT)
    end if


  END SUBROUTINE ENKF_SETUP_NC
!=======================================================|
   SUBROUTINE ENKF_RESTART
   USE MOD_SET_TIME
   IMPLICIT NONE
   if (kf_times>1) then
     write(enkf_num,'(i2.2)')  ENKF_memberCONTR

          NCF_ENKFre => NEW_FILE()
          NCF_ENKFre%FNAME=TRIM(OUTPUT_DIR)//"enkftemp/enkf_forecast_"//enkf_num//".nc"
         Call NC_OPEN(NCF_ENKFre)
         CALL NC_LOAD(NCF_ENKFre)
         CALL ENKF_READRESTART_last(NCF_ENKFre)
         CALL SET_LAST_FILE_TIME(NCF_ENKFre,StartTime, IINT)
         IntTime = starttime
         ExtTime = starttime
         CALL KILL_FILE(NCF_ENKFre)
   else
     write(enkf_num,'(i2.2)')  ENKF_memberCONTR
          NCF_ENKFre => NEW_FILE()
          NCF_ENKFre%FNAME=TRIM(INPUT_DIR)//"enkf_initial_"//enkf_num//".nc"
          Call NC_OPEN(NCF_ENKFre)
          CALL NC_LOAD(NCF_ENKFre)


         CALL ENKF_READRESTART_last(NCF_ENKFre)
         CALL SET_LAST_FILE_TIME(NCF_ENKFre,StartTime, IINT)
         starttime=enkf_start_time
         IntTime = enkf_start_time
         ExtTime = enkf_start_time
         iint=ekint_start
         CALL KILL_FILE(NCF_ENKFre)
   end if
   END SUBROUTINE ENKF_RESTART















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
!  FUNCTION SETUP_RESTART(FNAME) RESULT(NCF)
!    USE MOD_NCDIO
!    IMPLICIT NONE
!    CHARACTER(LEN=*) :: FNAME
!    TYPE(NCFILE), POINTER :: NCF
!
!    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START: SETUP_RESTART"
!
!
!    NCF => NEW_FILE()
!    NCF%FNAME=trim(fname)
!    print *, 'running setup_restart:',NCF%FNAME ! should be mark 
!    NCF%FTIME=>NEW_FTIME()
!
!    NCF => ADD(NCF,GRID_FILE_OBJECT() )
!
!    NCF => ADD(NCF,RESTART_EXTRAS_FILE_OBJECT() )
!
!    NCF => ADD(NCF,VELOCITY_FILE_OBJECT() )
!
!    NCF => ADD(NCF,AVERAGE_VEL_FILE_OBJECT() )
!
!    NCF => ADD(NCF,VERTICAL_VEL_FILE_OBJECT() )
!
!    NCF => ADD(NCF,TURBULENCE_FILE_OBJECT() )
!    
!    NCF => ADD(NCF,SALT_TEMP_FILE_OBJECT() )
!
!    IF(WETTING_DRYING_ON) THEN
!       NCF => ADD(NCF, WET_DRY_FILE_OBJECT() )
!    END IF
!
!    CALL NC_WRITE_FILE(NCF,LOCAL_DISK)
!    NCF%FTIME%NEXT_STKCNT = 1
!    
!
!
!    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END SETUP_RESTART"
!    
!  END FUNCTION SETUP_RESTART
!

# endif
END MODULE MOD_ENKF

