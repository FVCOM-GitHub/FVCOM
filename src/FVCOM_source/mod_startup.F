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

MODULE MOD_STARTUP
  USE MOD_UTILS
  USE MOD_NCTOOLS
  USE MOD_INPUT
  USE ALL_VARS
  USE EQS_OF_STATE
  USE MOD_WD
  USE SINTER


#  if defined (ICE)
  USE MOD_ICE
  USE MOD_ICE2D
#  endif

# if defined (WATER_QUALITY)
  USE MOD_WQM
# endif
  
# if defined (BioGen)
  USE MOD_BIO_3D
# endif

# if defined (NH)
  USE NON_HYDRO, ONLY: W4ZT, NHQDRX, NHQDRY, NHQDRZ, NHQ2DX, NHQ2DY
# endif

# if defined (DYE_RELEASE)
  USE MOD_DYE, ONLY: DYE
# endif
  
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: READ_SSH
  PUBLIC :: READ_UV
  PUBLIC :: READ_TURB
  PUBLIC :: READ_TS
  PUBLIC :: READ_WETDRY

  PUBLIC :: STARTUP
!# if defined (DATA_ASSIM)
!  PUBLIC :: STARTUP_ASSIM
!# endif 
CONTAINS
!==============================================================================!
  SUBROUTINE STARTUP
# if defined(SEDIMENT) 

# if defined(ORIG_SED)
  USE MOD_SED, ONLY: sed_hot_start
# elif defined(CSTMS_SED)
  USE MOD_SED_CSTMS, ONLY: sed_hot_start
# endif

# endif
    IMPLICIT NONE

    IF(DBG_SET(DBG_SBR)) write(ipt,*) "Start: Startup"

    IF(DBG_SET(DBG_LOG)) THEN
       WRITE(IPT,*  )'!'
       WRITE(IPT,*  )'!           SETTING INITIAL CONDITIONS '
       WRITE(IPT,*  )'!'
    END IF

    IF (ASSOCIATED(NC_START))THEN
       IF(DBG_SET(DBG_LOG)) write(ipt,*) "! NC_START FILE INFO:"
       CALL PRINT_FILE(NC_START)
    END IF
       
!    SET THE SEA SURFACE HEIGHT
    if (dbg_set(dbg_log)) write(ipt,*) &
         & "! INITIALIZING SEA SURFACE HEIGHT"
    SELECT CASE(STARTUP_TYPE)
! =================================================
    CASE(STARTUP_TYPE_HOTSTART)
! =================================================

       CALL READ_SSH

       IF(WETTING_DRYING_ON) CALL READ_WETDRY

# if defined (NH)
       CALL READ_NH
# endif

# if defined (ATMO_TIDE)
       CALL READ_ATMO
# endif

# if defined (EQUI_TIDE)
       CALL READ_EQI
# endif


! ------- New: Karsten Lettmann, June 2016 -----------------------
# if defined(AIR_PRESSURE)
  CALL READ_AIR_PRESS_ELE
# endif
! ---------- end  new --------------------------------------------

!    if (dbg_set(dbg_log)) write(ipt,*) &
!         & "! INITIALIZING SEA ICE"


# if defined (ICE)    
!    if (dbg_set(dbg_log)) write(ipt,*) &
!         & "! INITIALIZING SEA ICE1111"

       CALL  READ_ICE
       CALL  AGGREGATE  !!  ggao 07312008
# endif    

!    if (dbg_set(dbg_log)) write(ipt,*) &
!         & "! INITIALIZING SEA ICE222"

# if defined (SEDIMENT)
  if(sed_hot_start)then
    CALL READ_SED
  endif
# endif

# if defined (MEAN_FLOW)
    CALL READ_MEANFLOW_RST
# endif

# if defined (WAVE_CURRENT_INTERACTION)
    CALL READ_WAVE
# endif

# if defined (DYE_RELEASE)
    CALL READ_DYE
# endif

! =================================================
    CASE(STARTUP_TYPE_CRASHRESTART)
! =================================================
       CALL READ_SSH
     
       IF(WETTING_DRYING_ON) CALL READ_WETDRY

# if defined (NH)
       CALL READ_NH
# endif

# if defined (ATMO_TIDE)
       CALL READ_ATMO
# endif

# if defined (EQUI_TIDE)
       CALL READ_EQI
# endif


! ------- New: Karsten Lettmann, June 2016 -----------------------
# if defined(AIR_PRESSURE)
  CALL READ_AIR_PRESS_ELE
# endif
! ---------- end  new --------------------------------------------


# if defined (ICE)    
    CALL  READ_ICE
    CALL  AGGREGATE   
# endif    

# if defined (SEDIMENT)
  if(sed_hot_start)then
    CALL READ_SED
  endif
# endif

# if defined (WAVE_CURRENT_INTERACTION)
    CALL READ_WAVE
# endif

# if defined (DYE_RELEASE)
    CALL READ_DYE
# endif

! =================================================
    CASE(STARTUP_TYPE_COLDSTART)
! =================================================

       CALL SET_WATER_DEPTH
       IF(WETTING_DRYING_ON) CALL SET_WD_DATA

! =================================================
    CASE DEFAULT
! =================================================
       CALL FATAL_ERROR("STARTUP: UNKNOWN STARTUP TYPE: "//TRIM(STARTUP_TYPE))
    END SELECT


    if (dbg_set(dbg_log)) write(ipt,*) &
         & "! INITIALIZING VELOCITY FIELDS"

! SET STARTUP VALUES FOR VELOCITY
    SELECT CASE (STARTUP_UV_TYPE)
    CASE (STARTUP_TYPE_OBSERVED)
       CALL FATAL_ERROR("I DON'T KNOW HOW TO DO THAT KIND OF STARTUP")
    CASE(STARTUP_TYPE_LINEAR) 
       CALL FATAL_ERROR("I DON'T KNOW HOW TO DO THAT KIND OF STARTUP")
    CASE(STARTUP_TYPE_CONSTANT)
       !CALL FATAL_ERROR("I DON'T KNOW HOW TO DO THAT KIND OF STARTUP")
       CALL SET_CONSTANT_UV
    CASE(STARTUP_TYPE_DEFAULT)
       ! OKAY - it is zero
    CASE(STARTUP_TYPE_SETVALUES)

       CALL READ_UV

    CASE DEFAULT
       CALL FATAL_ERROR("UNKNOWN STARTUP_UV_TYPE")
    END SELECT

    if (dbg_set(dbg_log)) write(ipt,*) &
         & "! INITIALIZING TURBULENCE FIELDS"

! SET STARTUP VALUES FOR TURBULENCE
    SELECT CASE(STARTUP_TURB_TYPE)
    CASE(STARTUP_TYPE_OBSERVED) 
       CALL FATAL_ERROR("I DON'T KNOW HOW TO DO THAT KIND OF STARTUP")
    CASE(STARTUP_TYPE_LINEAR)
       CALL FATAL_ERROR("I DON'T KNOW HOW TO DO THAT KIND OF STARTUP")
    CASE(STARTUP_TYPE_CONSTANT)
       CALL FATAL_ERROR("I DON'T KNOW HOW TO DO THAT KIND OF STARTUP")
    CASE(STARTUP_TYPE_DEFAULT)

       CALL SET_DEFAULT_TURB

    CASE(STARTUP_TYPE_SETVALUES)

       CALL READ_TURB

    CASE DEFAULT
       CALL FATAL_ERROR("UNKNOWN STARTUP_TURB_TYPE")
    END SELECT

    if (dbg_set(dbg_log)) write(ipt,*) &
         & "! INITIALIZING TEMPERATURE AND SALINITY"

! SET STARTUP VALUES FOR TEMPERATER AND SALINITY
    SELECT CASE(STARTUP_TS_TYPE)
    CASE(STARTUP_TYPE_OBSERVED) 

       CALL SET_OBSERVED_TS

    CASE(STARTUP_TYPE_LINEAR)

       CALL SET_LINEAR_TS

    CASE(STARTUP_TYPE_CONSTANT)

       CALL SET_CONSTANT_TS

    CASE(STARTUP_TYPE_DEFAULT)
       CALL FATAL_ERROR("There is no default startup for Temperature and Salinity")
    CASE(STARTUP_TYPE_SETVALUES)

       CALL READ_TS

    CASE DEFAULT
       CALL FATAL_ERROR("UNKNOWN STARTUP_TS_TYPE")
    END SELECT

#   if defined (WATER_QUALITY)
    if (dbg_set(dbg_log)) write(ipt,*) &
         & "! INITIALIZING WATER QUALITY VARIABLES"

! SET STARTUP VALUES FOR WATER QUALITY VARIABLES
    SELECT CASE(STARTUP_WQM_TYPE)
    CASE(STARTUP_TYPE_OBSERVED) 

       CALL SET_OBSERVED_WQM

    CASE(STARTUP_TYPE_LINEAR)

       CALL SET_LINEAR_WQM

    CASE(STARTUP_TYPE_CONSTANT)

       CALL SET_CONSTANT_WQM

    CASE(STARTUP_TYPE_DEFAULT)
       CALL FATAL_ERROR("There is no default startup for Water Quality")
    CASE(STARTUP_TYPE_SETVALUES)

       CALL READ_WQM

    CASE DEFAULT
       CALL FATAL_ERROR("UNKNOWN STARTUP_WQM_TYPE")
    END SELECT
#   endif

#   if defined (BioGen)
    IF(BIOLOGICAL_MODEL)THEN

    if (dbg_set(dbg_log)) write(ipt,*) &
         & "! INITIALIZING BIOLOGICAL VARIABLES"

! SET STARTUP VALUES FOR BIOLOGICAL VARIABLES
    SELECT CASE(STARTUP_BIO_TYPE)
    CASE(STARTUP_TYPE_OBSERVED) 

       CALL SET_OBSERVED_BIO

    CASE(STARTUP_TYPE_LINEAR)

       CALL BIO_INITIAL

    CASE(STARTUP_TYPE_CONSTANT)

       CALL BIO_INITIAL

    CASE(STARTUP_TYPE_DEFAULT)
       CALL FATAL_ERROR("There is no default startup for biological model")
    CASE(STARTUP_TYPE_SETVALUES)

       CALL READ_BIO

    CASE DEFAULT
       CALL FATAL_ERROR("UNKNOWN STARTUP_BIO_TYPE")
    END SELECT

    ENDIF

#   endif

    ! ONCE WE HAVE STARTED THE MODEL POINT THE START FILE OBJECT TO
    ! ITS OWN OUTPUT TO RELOAD OLD TIME STATES
    IF(.not. associated(NC_START,NC_RST)) THEN
       
       ! SOME START UPS DO NOT HAVE A START FILE
       IF(Associated(NC_START)) CALL KILL_FILE(NC_START)

       NC_START => NC_RST
    END IF

    IF(DBG_SET(DBG_SBR)) write(ipt,*) "End: Startup"


  END SUBROUTINE STARTUP
!==============================================================================!

!# if defined (DATA_ASSIM)
!!==============================================================================!
!  SUBROUTINE STARTUP_ASSIM
!# if defined(SEDIMENT) 
!  USE MOD_SED, ONLY: sed_hot_start
!# endif
!    IMPLICIT NONE

!    IF(DBG_SET(DBG_SBR)) write(ipt,*) "Start: Startup_ASSIM"
       
!!    SET THE SEA SURFACE HEIGHT
!    if (dbg_set(dbg_log)) write(ipt,*) &
!         & "! INITIALIZING SEA SURFACE HEIGHT"
!       CALL READ_SSH

!       IF(WETTING_DRYING_ON) CALL READ_WETDRY

!# if defined (NH)
!       CALL READ_NH
!# endif

!# if defined (ATMO_TIDE)
!       CALL READ_ATMO
!# endif

!# if defined (EQUI_TIDE)
!       CALL READ_EQI
!# endif

!# if defined (ICE)    
!       CALL  READ_ICE
!       CALL  AGGREGATE  !!  ggao 07312008
!# endif    

!# if defined (SEDIMENT)
!  if(sed_hot_start)then
!    CALL READ_SED
!  endif
!# endif


!       CALL READ_UV

!       CALL READ_TURB

!       CALL READ_TS
!    if (dbg_set(dbg_log)) write(ipt,*) &
!         & "! INITIALIZING TS"


!#   if defined (WATER_QUALITY)

!       CALL READ_WQM

!#   endif

!#   if defined (BioGen)

!       CALL READ_BIO

!#   endif


!    IF(DBG_SET(DBG_SBR)) write(ipt,*) "End: Startup_Assim"



!  END SUBROUTINE STARTUP_ASSIM
!!==============================================================================!
!# endif

  SUBROUTINE READ_SSH
#   if defined (SEDIMENT)

#   if defined (ORIG_SED)
    USE MOD_SED, only : morpho_model,sed_hot_start 
# elif defined(CSTMS_SED)
  USE MOD_SED_CSTMS, ONLY:morpho_model,sed_hot_start
# endif

#   endif
    IMPLICIT NONE
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCDIM),  POINTER :: DIM
    LOGICAL :: FOUND
    INTEGER :: STKCNT

    IF(DBG_SET(DBG_SBR)) write(ipt,*) "Start: READ_SSH"
    
    STKCNT = NC_START%FTIME%PREV_STKCNT

    ! LOAD EL
    VAR => FIND_VAR(NC_START,'zeta',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'zeta'&
         & IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, EL)
    CALL NC_READ_VAR(VAR,STKCNT)

    ! LOAD ET
    VAR => FIND_VAR(NC_START,'et',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'et'&
         & IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, ET)
    CALL NC_READ_VAR(VAR,STKCNT)

    !----------------------------------------------------------------
    ! Read the most recent bathymetry if Morphodynamics is Active
    !----------------------------------------------------------------
#   if defined(SEDIMENT)
    IF(MORPHO_MODEL .and. SED_HOT_START)THEN
    VAR => FIND_VAR(NC_START,'h',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'h'&
         & IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, h)
    CALL NC_READ_VAR(VAR,STKCNT)
    ENDIF
#   endif


    !----------------------------------------------------------------
    ! Given SSH and Bathy, Update the Bathymetry 
    !----------------------------------------------------------------
    D  = H + EL
    DT = H + ET

    
    CALL N2E2D(H,H1)
    CALL N2E2D(EL,EL1)
    CALL N2E2D(D,D1)
    CALL N2E2D(DT,DT1)

    IF(DBG_SET(DBG_SBR)) write(ipt,*) "Start: READ_SSH"

  END SUBROUTINE READ_SSH
!==============================================================================!
!==============================================================================!
  SUBROUTINE READ_WETDRY
    USE MOD_WD
    IMPLICIT NONE
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCDIM),  POINTER :: DIM
    LOGICAL :: FOUND
    INTEGER :: STKCNT

    IF(DBG_SET(DBG_SBR)) write(ipt,*) "Start: READ_WETDRY"

    STKCNT = NC_START%FTIME%PREV_STKCNT

    ! LOAD ISWETN
    VAR => FIND_VAR(NC_START,'wet_nodes',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'wet_nodes'&
         & IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, ISWETN)
    CALL NC_READ_VAR(VAR,STKCNT)

    ! LOAD ISWETC
    VAR => FIND_VAR(NC_START,'wet_cells',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'wet_cells'&
         & IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, ISWETC)
    CALL NC_READ_VAR(VAR,STKCNT)


    ! LOAD ISWETNT
    VAR => FIND_VAR(NC_START,'wet_nodes_prev_int',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'wet_nodes_prev_int'&
         & IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, ISWETNT)
    CALL NC_READ_VAR(VAR,STKCNT)

    ! LOAD ISWETCT
    VAR => FIND_VAR(NC_START,'wet_cells_prev_int',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'wet_cells_prev_int'&
         & IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, ISWETCT)
    CALL NC_READ_VAR(VAR,STKCNT)

    ! LOAD ISWETCE
    VAR => FIND_VAR(NC_START,'wet_cells_prev_ext',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'wet_cells_prev_ext'&
         & IN THE HOTSTART FILE OBJECT")
    ! ------ new: Karsten Lettmann, May 2016 -------------
    !CALL NC_CONNECT_AVAR(VAR, ISWETC)  ! original line
    CALL NC_CONNECT_AVAR(VAR, ISWETCE)
    ! ------------- end new ------------------------------
    CALL NC_READ_VAR(VAR,STKCNT)

    IF(DBG_SET(DBG_SBR)) write(ipt,*) "End: READ_WETDRY"
    
  END SUBROUTINE READ_WETDRY
!==============================================================================!
# if defined (NH)
  SUBROUTINE READ_NH
    IMPLICIT NONE
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCDIM),  POINTER :: DIM
    LOGICAL :: FOUND
    INTEGER :: STKCNT
    
    IF(DBG_SET(DBG_SBR)) write(ipt,*) "START: READ_NH" 

    STKCNT = NC_START%FTIME%PREV_STKCNT

    ! LOAD w4zt
    VAR => FIND_VAR(NC_START,'W4ZT',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR&
         &("COULD NOT FIND VARIABLE 'W4ZT' IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, W4ZT)
    CALL NC_READ_VAR(VAR,STKCNT)

    ! LOAD NHQDRX
    VAR => FIND_VAR(NC_START,'NHQDRX',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR&
         &("COULD NOT FIND VARIABLE 'NHQDRX' IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, NHQDRX)
    CALL NC_READ_VAR(VAR,STKCNT)

    ! LOAD NHQDRY
    VAR => FIND_VAR(NC_START,'NHQDRY',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR&
         &("COULD NOT FIND VARIABLE 'NHQDRY' IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, NHQDRY)
    CALL NC_READ_VAR(VAR,STKCNT)

    ! LOAD NHQDRZ
    VAR => FIND_VAR(NC_START,'NHQDRZ',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR&
         &("COULD NOT FIND VARIABLE 'NHQDRZ' IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, NHQDRZ)
    CALL NC_READ_VAR(VAR,STKCNT)

    ! LOAD NHQ2DX
    VAR => FIND_VAR(NC_START,'NHQ2DX',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR&
         &("COULD NOT FIND VARIABLE 'NHQ2DX' IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, NHQ2DX)
    CALL NC_READ_VAR(VAR,STKCNT)

    ! LOAD NHQ2DY
    VAR => FIND_VAR(NC_START,'NHQ2DY',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR&
         &("COULD NOT FIND VARIABLE 'NHQ2DY' IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, NHQ2DY)
    CALL NC_READ_VAR(VAR,STKCNT)

    IF(DBG_SET(DBG_SBR)) write(ipt,*) "End: READ_NH"

  END SUBROUTINE READ_NH
# endif
!==============================================================================!
  SUBROUTINE READ_ATMO
    IMPLICIT NONE
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCDIM),  POINTER :: DIM
    LOGICAL :: FOUND
    INTEGER :: STKCNT
    
    IF(DBG_SET(DBG_SBR)) write(ipt,*) "START: READ_ATMO" 

    STKCNT = NC_START%FTIME%PREV_STKCNT

    ! LOAD ATMO
    VAR => FIND_VAR(NC_START,'el_atmo',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR&
         &("COULD NOT FIND VARIABLE 'el_atmo' IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, el_atmo)
    CALL NC_READ_VAR(VAR,STKCNT)

    IF(DBG_SET(DBG_SBR)) write(ipt,*) "End: READ_ATMO"

  END SUBROUTINE READ_ATMO
!==============================================================================!
  SUBROUTINE READ_EQI
    IMPLICIT NONE
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCDIM),  POINTER :: DIM
    LOGICAL :: FOUND
    INTEGER :: STKCNT

    IF(DBG_SET(DBG_SBR)) write(ipt,*) "Start: READ_EQI" 
    
    STKCNT = NC_START%FTIME%PREV_STKCNT

    ! LOAD ATMO
    VAR => FIND_VAR(NC_START,'el_eqi',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR&
         &("COULD NOT FIND VARIABLE 'el_eqi' IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, el_eqi)
    CALL NC_READ_VAR(VAR,STKCNT)

    IF(DBG_SET(DBG_SBR)) write(ipt,*) "End: READ_EQI" 
   
  END SUBROUTINE READ_EQI
!==============================================================================!

! ------- New: Karsten Lettmann, June 2016 -----------------------
# if defined(AIR_PRESSURE)
!==============================================================================!
  SUBROUTINE READ_AIR_PRESS_ELE
    IMPLICIT NONE
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCDIM),  POINTER :: DIM
    LOGICAL :: FOUND
    INTEGER :: STKCNT

    IF(DBG_SET(DBG_SBR)) write(ipt,*) "Start: READ_AIR_PRESS_ELE" 
    
    STKCNT = NC_START%FTIME%PREV_STKCNT

    ! LOAD Air pressure induced elevation
    VAR => FIND_VAR(NC_START,'el_press',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR&
         &("COULD NOT FIND VARIABLE 'el_press' IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, EL_AIR)
    CALL NC_READ_VAR(VAR,STKCNT)

    IF(DBG_SET(DBG_SBR)) write(ipt,*) "End: READ_AIR_PRESS_ELE" 
   
  END SUBROUTINE READ_AIR_PRESS_ELE
!==============================================================================!
# endif  
! ------------- end new ---------------------------------------------------------

# if defined (DYE_RELEASE)
  SUBROUTINE READ_DYE
    IMPLICIT NONE
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCDIM),  POINTER :: DIM
    LOGICAL :: FOUND
    INTEGER :: STKCNT, K

    IF(DBG_SET(DBG_SBR)) write(ipt,*) "Start: READ_DYE"

    STKCNT = NC_START%FTIME%PREV_STKCNT


    ! LOAD DYE
    VAR => FIND_VAR(NC_START,'DYE',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'DYE'&
         & IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, DYE)
    CALL NC_READ_VAR(VAR,STKCNT)

    IF(DBG_SET(DBG_SBR)) write(ipt,*) "End: READ_DYE"

  END SUBROUTINE READ_DYE
# endif
!==============================================================================!
  SUBROUTINE READ_TS
    IMPLICIT NONE
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCDIM),  POINTER :: DIM
    LOGICAL :: FOUND
    INTEGER :: STKCNT, K
    REAL(SP), DIMENSION(0:MT,KB) :: PRESSURE
    
    IF(DBG_SET(DBG_SBR)) write(ipt,*) "Start: READ_TS" 
   
    STKCNT = NC_START%FTIME%PREV_STKCNT


    ! LOAD TEMPERATURE
    VAR => FIND_VAR(NC_START,'temp',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'temp'&
         & IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, T1)
    CALL NC_READ_VAR(VAR,STKCNT)


    ! LOAD MEAN INITIAL TEMPERATURE
    VAR => FIND_VAR(NC_START,'tmean1',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'tmean1'&
         & IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, tmean1)
    CALL NC_READ_VAR(VAR)

    ! LOAD SALINITY
    VAR => FIND_VAR(NC_START,'salinity',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'saltinity'&
         & IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, S1)
    CALL NC_READ_VAR(VAR,STKCNT)

    ! LOAD MEAN INITIAL SALINITY
    VAR => FIND_VAR(NC_START,'smean1',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'smean1'&
         & IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, smean1)
    CALL NC_READ_VAR(VAR)

    ! LOAD DENSITY
    VAR => FIND_VAR(NC_START,'rho1',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'rho1'&
         & IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, RHO1)
    CALL NC_READ_VAR(VAR,STKCNT)

    ! LOAD MEAN DENSITY
    VAR => FIND_VAR(NC_START,'rmean1',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'rmean1'&
         & IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, rmean1)
    CALL NC_READ_VAR(VAR)


    ! AVERAGE FROM CELLS TO FACE CENTERS


!JQI    SELECT CASE(SEA_WATER_DENSITY_FUNCTION)
!JQI    CASE(SW_DENS1)

!JQI       ! SET MEAN DENSITY
!JQI       DO K=1,KBM1
!JQI          PRESSURE(:,K) = -GRAV_N*1.025_SP*(ZZ(:,K)*D(:))*0.1_SP
!JQI       END DO
!JQI       CALL FOFONOFF_MILLARD(SMEAN1,TMEAN1,PRESSURE,0.0_SP,RMEAN1)
!JQI       RMEAN1(0,:)=0.0_SP
!JQI       RMEAN1(:,KB)=0.0_SP

!JQI       ! SET REAL DENSITY
!JQI       CALL DENS1 ! GENERIC CALL TO FOFONOFF_MILLARD FOR S1,T1...
       
!JQI    CASE(SW_DENS2)
!JQI       ! SET MEAN DENSITY
!JQI       CALL DENS2G(SMEAN1,TMEAN1,RMEAN1)
!JQI       RMEAN1(0,:)=0.0_SP
!JQI       RMEAN1(:,KB)=0.0_SP

!JQI       ! SET REAL DENSITY
!JQI       CALL DENS2 ! GENERIC CALL TO DENS2G FOR S1,T1...

!JQI    CASE(SW_DENS3)

!JQI       ! SET MEAN DENSITY
!JQI       DO K=1,KBM1
!JQI          PRESSURE(:,K) = -GRAV_N*1.025_SP*(ZZ(:,K)*D(:))*0.01_SP
!JQI       END DO
!JQI       CALL JACKET_MCDOUGALL(SMEAN1,TMEAN1,PRESSURE,RMEAN1)
!JQI       RMEAN1(0,:)=0.0_SP
!JQI       RMEAN1(:,KB)=0.0_SP

!JQI       ! SET REAL DENSITY
!JQI       CALL DENS3 ! GENERIC CALL TO JACKET_MCDOUGALL FOR S1,T1.

!JQI    CASE DEFAULT
!JQI       CALL FATAL_ERROR("INVALID DENSITY FUNCTION SELECTED:",&
!JQI            & "   "//TRIM(SEA_WATER_DENSITY_FUNCTION) )
!JQI    END SELECT

!J. GE
    T0=T1
    S0=S1
!J. GE    
    CALL N2E3D(T1,T)
    CALL N2E3D(S1,S)
    CALL N2E3D(RHO1,RHO)
    CALL N2E3D(Tmean1,Tmean)
    CALL N2E3D(Smean1,Smean)
    CALL N2E3D(Rmean1,Rmean)

   
   IF(DBG_SET(DBG_SBR)) write(ipt,*) "End: READ_TS" 
   
 END SUBROUTINE READ_TS
!==============================================================================!
  SUBROUTINE SET_OBSERVED_TS
    IMPLICIT NONE
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCDIM),  POINTER :: DIM
    LOGICAL :: FOUND
    INTEGER :: STKCNT
    INTEGER KSL                                 !!NUMBER OF STANDARD SEA LEVELS 
    REAL(SP), ALLOCATABLE, TARGET :: DPTHSL(:)          !!DEPTH AT STANDARD SEA LEVEL
    REAL(SP), ALLOCATABLE, TARGET :: TSL(:,:),SSL(:,:)  !!T/S AT STANDARD SEA LEVEL
    REAL(SP),ALLOCATABLE          :: TA(:),SA(:)

    REAL(SP),DIMENSION(KBM1)      :: TI,SI,ZI
    REAL(SP),DIMENSION(0:MT,KB) :: PRESSURE
    INTEGER :: K, I, IB, IERR
    REAL(SP) :: SCOUNT, FAC, RBUF
    
    IF(DBG_SET(DBG_SBR)) write(ipt,*) "START: SET_OBSERVED_TS" 

    TI       = 0.0_SP
    SI       = 0.0_SP
    ZI       = 0.0_SP
    PRESSURE = 0.0_SP

    STKCNT = NC_START%FTIME%PREV_STKCNT

    DIM => FIND_DIM(NC_START,'ksl',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND DIMENSION 'ksl'&
         & IN THE STARTUP FILE OBJECT")

    KSL = DIM%DIM
    ! ALLOCATE SPACE
    ALLOCATE(TSL(MT,KSL))
    ALLOCATE(SSL(MT,KSL))
    ALLOCATE(DPTHSL(KSL))
    TSL   = 0.0_SP
    SSL   = 0.0_SP
    DPTHSL= 0.0_SP

    ! ALLOCATE SPACE FOR THE AVERAGE TEMPERATURE AT EACH DEPTH
    ALLOCATE(TA(KSL),SA(KSL))
    TA = 0.0_SP
    SA = 0.0_SP

    ! READ THE DATA
    VAR => FIND_VAR(NC_START,'tsl',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'tsl'&
         & IN THE STARTUP FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, TSL)
    CALL NC_READ_VAR(VAR,STKCNT)
    CALL NC_DISCONNECT(VAR)

    VAR => FIND_VAR(NC_START,'ssl',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'ssl'&
         & IN THE STARTUP FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, SSL)
    CALL NC_READ_VAR(VAR,STKCNT)
    CALL NC_DISCONNECT(VAR)


    VAR => FIND_VAR(NC_START,'zsl',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'zsl'&
         & IN THE STARTUP FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, DPTHSL)
    CALL NC_READ_VAR(VAR)
    CALL NC_DISCONNECT(VAR)


!==============================================================================|
!    CALCULATE HORIZONTAL AVERAGE VALUES OF TEMPERATURE/SALINITY/DENSITY       !
!==============================================================================|
   
      where(ssl <0.0_SP) ssl=0.0
 
    ! GET THE AVERAGE T/S AT EACH DEPTH ON STANDARD LEVELS
    IF(SERIAL)THEN
       DO K = 1, KSL
          SCOUNT = 0.0_SP
          DO I = 1, MT
             IF(-H(I) <= DPTHSL(K)) THEN
                SCOUNT = SCOUNT + 1.0_SP
                TA(K) = TA(K) + TSL(I,K)
                SA(K) = SA(K) + SSL(I,K)
             END IF
          END DO
          IF(SCOUNT >= 1.0_SP)THEN
             TA(K)=TA(K)/SCOUNT
             SA(K)=SA(K)/SCOUNT
          END IF
       END DO
    END IF
    
    
#    if defined(MULTIPROCESSOR)
    IF(PAR)THEN
       
       DO K = 1, KSL
          SCOUNT = 0.0_SP
          DO I = 1, M
             IF(-H(I) <= DPTHSL(K)) THEN
                IF(NDE_ID(I) == 0)THEN !!INTERNAL NODE
                   SCOUNT = SCOUNT + 1.0_SP
                   TA(K) = TA(K) + TSL(I,K)
                   SA(K) = SA(K) + SSL(I,K)
                ELSE  !!BOUNDARY NODE, ACCUMULATE FRACTION ONLY
                   DO IB = 1,NBN
                      IF(BN_LOC(IB) == I)THEN
                         FAC = 1.0_SP/FLOAT(BN_MLT(IB))
                         SCOUNT = SCOUNT + FAC 
                         TA(K) = TA(K) + FAC*TSL(I,K)
                         SA(K) = SA(K) + FAC*SSL(I,K)
                      END IF
                   END DO
                END IF
             END IF
          END DO
          CALL MPI_ALLREDUCE(TA(K),RBUF,1,MPI_F,MPI_SUM,MPI_FVCOM_GROUP,IERR)
          TA(K) = RBUF
          CALL MPI_ALLREDUCE(SA(K),RBUF,1,MPI_F,MPI_SUM,MPI_FVCOM_GROUP,IERR)
          SA(K) = RBUF
          CALL MPI_ALLREDUCE(SCOUNT,RBUF,1,MPI_F,MPI_SUM,MPI_FVCOM_GROUP,IERR)
          SCOUNT = RBUF

          IF(SCOUNT >= 1.0_SP)THEN
             TA(K)=TA(K)/SCOUNT
             SA(K)=SA(K)/SCOUNT
          END IF
          
       END DO !!LOOP OVER KSL
    END IF  !!PARALLEL
#    endif
    
    ! NOW INTERPOLATE FROM STANDARD LEVELS TO THE VALUE AT EACH NODE
    DO I=1,MT
       DO K=1,KBM1
          ZI(K)=ZZ(I,K)*D(I)+EL(I)
       END DO
       
       ! LEVEL AVERAGE T AND S
       CALL SINTER_EXTRP_UP(DPTHSL,TA,ZI,TI,KSL,KBM1)
       CALL SINTER_EXTRP_UP(DPTHSL,SA,ZI,SI,KSL,KBM1)
       
       TMEAN1(I,1:KBM1) = TI(:)
       SMEAN1(I,1:KBM1) = SI(:)
              
       ! REAL T AND S
       CALL SINTER_EXTRP_UP(DPTHSL,TSL(I,:),ZI,TI,KSL,KBM1)
       CALL SINTER_EXTRP_UP(DPTHSL,SSL(I,:),ZI,SI,KSL,KBM1)
       
       T1(I,1:KBM1) = TI(:)
       S1(I,1:KBM1) = SI(:)
!J. Ge
       T0(I,1:KBM1) = TI(:)
       S0(I,1:KBM1) = SI(:)
!J. Ge
    END DO

       where(S1<0.0_SP) S1=0.0

    IF(.NOT.BAROTROPIC)THEN
     SELECT CASE(SEA_WATER_DENSITY_FUNCTION)
     CASE(SW_DENS1)

       ! SET MEAN DENSITY
       DO K=1,KBM1
          PRESSURE(:,K) = -GRAV_N*1.025_SP*(ZZ(:,K)*D(:))*0.1_SP
       END DO
       CALL FOFONOFF_MILLARD(SMEAN1,TMEAN1,PRESSURE,0.0_SP,RMEAN1)
       RMEAN1(0,:)=0.0_SP
       RMEAN1(:,KB)=0.0_SP

       ! SET REAL DENSITY
       CALL DENS1 ! GENERIC CALL TO FOFONOFF_MILLARD FOR S1,T1...

     CASE(SW_DENS2)

       ! SET MEAN DENSITY
       CALL DENS2G(SMEAN1,TMEAN1,RMEAN1)
       RMEAN1(0,:)=0.0_SP
       RMEAN1(:,KB)=0.0_SP

       ! SET REAL DENSITY
       CALL DENS2 ! GENERIC CALL TO DENS2G FOR S1,T1...

     CASE(SW_DENS3)
     
      ! SET MEAN DENSITY
       DO K=1,KBM1
          PRESSURE(:,K) = -GRAV_N*1.025_SP*(ZZ(:,K)*D(:))*0.01_SP
       END DO
       CALL JACKET_MCDOUGALL(SMEAN1,TMEAN1,PRESSURE,RMEAN1)
       RMEAN1(0,:)=0.0_SP
       RMEAN1(:,KB)=0.0_SP

       ! SET REAL DENSITY
       CALL DENS3 ! GENERIC CALL TO JACKET_MCDOUGALL FOR S1,T1.

     CASE DEFAULT
      CALL FATAL_ERROR("INVALID DENSITY FUNCTION SELECTED:",&
           & "   "//TRIM(SEA_WATER_DENSITY_FUNCTION) )
     END SELECT
    ELSE
     RHO1   = 2.3E-2_SP
     RHO    = 2.3E-2_SP
     RMEAN1 = 2.3E-2_SP
    END IF 
   

   CALL N2E3D(T1,T)
   CALL N2E3D(S1,S)
   CALL N2E3D(Tmean1,Tmean)
   CALL N2E3D(Smean1,Smean)
   CALL N2E3D(Rmean1,Rmean)
   ! DENSITY IS ALREADY INTERPOLATED TO ELEMENTS FOR RHO IN DENSX

   DEALLOCATE(TSL,SSL,DPTHSL)
   DEALLOCATE(TA,SA)


   IF(DBG_SET(DBG_SBR)) write(ipt,*) "END: SET_OBSERVED_TS" 
    
  END SUBROUTINE SET_OBSERVED_TS
!==============================================================================!
  SUBROUTINE SET_LINEAR_TS
    IMPLICIT NONE

    ! BY DEFINITION OF LINEAR...
    INTEGER, PARAMETER :: KSL =2

    REAL(SP), DIMENSION(KBM1)   :: TI,SI,ZI
    REAL(SP), DIMENSION(KSL)   :: TA,SA,ZA
    INTEGER :: K, I

    IF(DBG_SET(DBG_SBR)) write(ipt,*) "START: SET_LINEAR_TS" 

    ZA(1) = 0.0_SP
    ZA(2) = STARTUP_DMAX
    
    TA(1) = STARTUP_T_VALS(1)
    TA(2) = STARTUP_T_VALS(2)
    
    SA(1) = STARTUP_S_VALS(1)
    SA(2) = STARTUP_S_VALS(2)
    
    ! THE HORIZONTAL AVERAGE VALUE IS ALSO THE TRUE VALUE SINCE THE
    ! LINEAR EQUATION DOES NOT VARY WITH LOCATION 

    DO I = 1, MT
       DO K=1,KBM1
          ! CALCULATE ZI relative to z=0
          ZI(K)=ZZ(I,K)*D(I)+EL(I)
       END DO
    
       CALL SINTER_EXTRP_UP(ZA,TA,ZI,TI,KSL,KBM1)
       CALL SINTER_EXTRP_UP(ZA,SA,ZI,SI,KSL,KBM1)
    
       T1(I,1:kbm1) = TI
       S1(I,1:kbm1) = SI
!J. Ge
       T0(I,1:kbm1) = TI
       S0(I,1:kbm1) = SI
!J. Ge

    END DO

    IF(.NOT.BAROTROPIC)THEN
     SELECT CASE(SEA_WATER_DENSITY_FUNCTION)
     CASE(SW_DENS1)
       CALL DENS1 ! USE GENERIC INTERFACE TO DENSX
     CASE(SW_DENS2)
       CALL DENS2
     CASE(SW_DENS3)
       CALL DENS3
     CASE DEFAULT
       CALL FATAL_ERROR("INVALID DENSITY FUNCTION SELECTED:",&
           & "   "//TRIM(SEA_WATER_DENSITY_FUNCTION) )
     END SELECT
    ELSE
     RHO1 = 2.3E-2_SP
     RHO  = 2.3E-2_SP
    END IF 
   
   ! FOR LINEAR AND CONSTANT, THE MEAN IS EQUAL TO THE INITIAL VALUE

   TMEAN1=T1
   SMEAN1=S1
   RMEAN1=RHO1
   RMEAN=RHO
   CALL N2E3D(T1,T)
   CALL N2E3D(S1,S)
   TMEAN=T
   SMEAN=S
   
    IF(DBG_SET(DBG_SBR)) write(ipt,*) "END: SET_LINEAR_TS" 
    
  END SUBROUTINE SET_LINEAR_TS
!==============================================================================!
  SUBROUTINE SET_CONSTANT_TS
    IMPLICIT NONE
    
    IF(DBG_SET(DBG_SBR)) write(ipt,*) "START: SET_CONSTANT_TS" 


    T1(:,1:KBM1) = STARTUP_T_VALS(1)
    S1(:,1:KBM1) = STARTUP_S_VALS(1)
!J. Ge
    T0(:,1:KBM1) = STARTUP_T_VALS(1)
    S0(:,1:KBM1) = STARTUP_S_VALS(1)
!J. Ge

    IF(.NOT.BAROTROPIC)THEN
     SELECT CASE(SEA_WATER_DENSITY_FUNCTION)
     CASE(SW_DENS1)
       CALL DENS1
     CASE(SW_DENS2)
       CALL DENS2
     CASE(SW_DENS3)
       CALL DENS3
     CASE DEFAULT
       CALL FATAL_ERROR("INVALID DENSITY FUNCTION SELECTED:",&
            & "   "//TRIM(SEA_WATER_DENSITY_FUNCTION) )
     END SELECT
    ELSE
     RHO1 = 2.3E-2_SP
     RHO  = 2.3E-2_SP
    END IF 

   ! FOR LINEAR AND CONSTANT, THE MEAN IS EQUAL TO THE INITIAL VALUE

   TMEAN1=T1
   SMEAN1=S1
   RMEAN1=RHO1

   CALL N2E3D(T1,T)
   CALL N2E3D(S1,S)
   
   TMEAN=T
   RMEAN=RHO
   SMEAN=S

   IF(DBG_SET(DBG_SBR)) write(ipt,*) "END: SET_CONSTANT_TS" 

 END SUBROUTINE SET_CONSTANT_TS
!==============================================================================!
!==============================================================================!
  SUBROUTINE SET_CONSTANT_UV
    IMPLICIT NONE

    IF(DBG_SET(DBG_SBR)) write(ipt,*) "START: SET_CONSTANT_UV"


    U(:,1:KBM1) = STARTUP_U_VALS    
    V(:,1:KBM1) = STARTUP_V_VALS
    UA          = STARTUP_U_VALS
    VA          = STARTUP_V_VALS
    WTS         = 0.0          
    W           = 0.0

   IF(DBG_SET(DBG_SBR)) write(ipt,*) "END: SET_CONSTANT_UV"

 END SUBROUTINE SET_CONSTANT_UV
!==============================================================================!

  SUBROUTINE READ_UV
    IMPLICIT NONE
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCDIM),  POINTER :: DIM
    LOGICAL :: FOUND
    INTEGER :: STKCNT

    IF(DBG_SET(DBG_SBR)) write(ipt,*) "START: READ_UV" 


    STKCNT = NC_START%FTIME%PREV_STKCNT

    VAR => FIND_VAR(NC_START,'u',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'u'&
         & IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, U)
    CALL NC_READ_VAR(VAR,STKCNT)


    VAR => FIND_VAR(NC_START,'v',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'v'&
         & IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, V)
    CALL NC_READ_VAR(VAR,STKCNT)

    VAR => FIND_VAR(NC_START,'omega',FOUND)
    IF(FOUND) THEN
       CALL NC_CONNECT_AVAR(VAR, WTS)
       CALL NC_READ_VAR(VAR,STKCNT)
       
       CALL N2E3D(WTS,W)
    ELSE
       VAR => FIND_VAR(NC_START,'w',FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'w' &
            & or 'omega' IN THE HOTSTART FILE OBJECT")
       CALL NC_CONNECT_AVAR(VAR, W)
       CALL NC_READ_VAR(VAR,STKCNT)
    END IF

    VAR => FIND_VAR(NC_START,'ua',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'ua'&
         & IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, UA)
    CALL NC_READ_VAR(VAR,STKCNT)

    VAR => FIND_VAR(NC_START,'va',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'va'&
         & IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, VA)
    CALL NC_READ_VAR(VAR,STKCNT)

    IF(DBG_SET(DBG_SBR)) write(ipt,*) "END: READ_UV" 

  END SUBROUTINE READ_UV
!==============================================================================!
  SUBROUTINE READ_TURB
    IMPLICIT NONE
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCDIM),  POINTER :: DIM
    LOGICAL :: FOUND
    INTEGER :: STKCNT

    IF(DBG_SET(DBG_SBR)) write(ipt,*) "START: READ_TURB" 

    STKCNT = NC_START%FTIME%PREV_STKCNT

#    if defined (GOTM)

    VAR => FIND_VAR(NC_START,'tke',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'tke'&
         & IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, TKE)
    CALL NC_READ_VAR(VAR,STKCNT)

    VAR => FIND_VAR(NC_START,'teps',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'teps'&
         & IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, TEPS)
    CALL NC_READ_VAR(VAR,STKCNT)

    L = .001 
    L(1:MT,2:KBM1) = (.5544**3)*TKE(1:MT,2:KBM1)**1.5/TEPS(1:MT,2:KBM1)
    
# else

    VAR => FIND_VAR(NC_START,'q2',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'q2'&
         & IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, Q2)
    CALL NC_READ_VAR(VAR,STKCNT)

    VAR => FIND_VAR(NC_START,'q2l',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'q2l'&
         & IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, Q2L)
    CALL NC_READ_VAR(VAR,STKCNT)

    VAR => FIND_VAR(NC_START,'l',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'l'&
         & IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, L)
    CALL NC_READ_VAR(VAR,STKCNT)

# endif

    VAR => FIND_VAR(NC_START,'km',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'km'&
         & IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, KM)
    CALL NC_READ_VAR(VAR,STKCNT)

    VAR => FIND_VAR(NC_START,'kq',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'kq'&
         & IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, KQ)
    CALL NC_READ_VAR(VAR,STKCNT)

    VAR => FIND_VAR(NC_START,'kh',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'kh'&
         & IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, KH)
    CALL NC_READ_VAR(VAR,STKCNT)


    CALL N2E3D(KM,KM1)

    IF(DBG_SET(DBG_SBR)) write(ipt,*) "END: READ_TURB" 

  END SUBROUTINE READ_TURB
!==============================================================================|
   SUBROUTINE SET_DEFAULT_TURB         
!==============================================================================|
!   Initialize Turbulent Kinetic Energy and Length Scale                       |
!==============================================================================|
   IMPLICIT NONE
   INTEGER :: I,K
!==============================================================================|
   IF(DBG_SET(DBG_SBR)) write(ipt,*) "START: SET_DEFAULT_TURB" 
   
   DO K=1,KBM1
        AAM(:,K) = NN_HVC(:)
   END DO

!
!------------------------BOUNDARY VALUES---------------------------------------!
!

   DO I = 1, MT
     KM(I,1)   = 0.0_SP
     KM(I,KB)  = 0.0_SP
     KH(I,1)   = 0.0_SP
     KH(I,KB)  = 0.0_SP
     KQ(I,1)   = 0.0_SP
     KQ(I,KB)  = 0.0_SP
     L(I,1)    = 0.0_SP
     L(I,KB)   = 0.0_SP
     Q2(I,1)   = 0.0_SP
     Q2(I,KB)  = 0.0_SP
     Q2L(I,1)  = 0.0_SP
     Q2L(I,KB) = 0.0_SP
#    if defined (GOTM)
     TKE(I,1)  = 0.0_SP
    TEPS(I,1)  = 0.0_SP
     TKE(I,KB) = 0.0_SP
    TEPS(I,KB) = 0.0_SP
#    endif
   END DO


!
!------------------------INTERNAL VALUES---------------------------------------!
!
   DO  K = 2, KBM1
     DO I = 1, MT
       IF (D(I) > 0.0_SP) THEN
         Q2(I,K)  = 1.E-8
         Q2L(I,K) = 1.E-8
         L(I,K)   = 1.
#        if defined (GOTM)
         TKE(I,K) = Q2(I,K)/2.
         TEPS(I,K)= (.1704*TKE(I,K)**1.5)/L(I,K)
#        endif
         KM(I,K)  = 2.*UMOL
         KQ(I,K)  = 2.*UMOL
         KH(I,K)  = 2.*UMOL / VPRNU
       END IF
     END DO
   END DO
   CALL N2E3D(KM,KM1)

   IF(DBG_SET(DBG_SBR)) write(ipt,*) "END: SET_DEFAULT_TURB" 

 END SUBROUTINE SET_DEFAULT_TURB
!==============================================================================|

!==============================================================================|
!   READ IN STATIC WATER DEPTH AND CALCULATE RELATED QUANTITIES                |
!                                                                              |
!   INPUTS: H(NNODE) BATHYMETRIC DEPTH AT NODES				       |
!   INITIALIZES: D(NNODE) DEPTH AT NODES				       |
!   INITIALIZES: DT(NNODE) ???					               |
!   INITIALIZES: H1(NNODE) BATHYMETRIC DEPTH AT ELEMENTS		       |
!   INITIALIZES: D1(NNODE) DEPTH AT NODES                		       | 
!   INITIALIZES: DT1(NNODE) ??                                   	       |
!==============================================================================|

   SUBROUTINE SET_WATER_DEPTH         
!------------------------------------------------------------------------------|
   USE MOD_OBCS
   USE MOD_NESTING
   IMPLICIT NONE
   REAL(SP) :: TEMP
   INTEGER  :: I,K,J1,J2
!------------------------------------------------------------------------------|
   IF(DBG_SET(DBG_SBR)) write(ipt,*) "START: SET_WATER_DEPTH" 


!
!  ADJUST STATIC HEIGHT AND CALCULATE DYNAMIC DEPTHS (D) AND (DT)
!
   H  = H + STATIC_SSH_ADJ
   D  = H + EL
   DT = H + ET

!
!  ADJUST SIGMA VALUES ON OUTER BOUNDARY 
!

#  if !defined (PLBC)
   IF(IOBCN > 0 .AND. OBC_DEPTH_CONTROL_ON) THEN
     IF(.NOT. NESTING_ON)THEN
#  if defined (WAVE_CURRENT_INTERACTION)     
     IF(.NOT. NESTING_ON_WAVE)THEN
#  endif      
     DO I = 1,IOBCN
       J1 = I_OBC_N(I)
       J2 = NEXT_OBC(I)
       H(J1) = H(J2)
       D(J1) = D(J2)
       DT(J1) = DT(J2)
       
       DO K = 1,KB
         Z(J1,K)   = Z(J2,K)
         ZZ(J1,K)  = ZZ(J2,K)
         DZ(J1,K)  = DZ(J2,K)
         DZZ(J1,K) = DZZ(J2,K)
       END DO	 
     END DO
#  if defined (WAVE_CURRENT_INTERACTION)     
     END IF
#  endif      
     END IF
   END IF
#  endif

! 
!  CALCULATE FACE-CENTERED VALUES OF BATHYMETRY AND DEPTH
!
   DO I=1,NT    
     H1(I)  = (H(NV(I,1))+H(NV(I,2))+H(NV(I,3)))/3.0_SP
     D1(I)  = H1(I)+EL1(I)
     DT1(I) = H1(I)+ET1(I)
     
     DO K = 1,KB
       Z1(I,K)=(Z(NV(I,1),K)+Z(NV(I,2),K)+Z(NV(I,3),K))/3.0_SP
     END DO  
   END DO

!-----AFTER MODIFYING BOUNDARY SIGMA VALUES ------------------------------------!
!-----RECOMPUTE SIGMA DERIVATIVES AND INTRA SIGMA LEVELS AGAIN ON CELL----------!

   DO K=1,KB-1
     DO I=1,NT
       DZ1(I,K)  = Z1(I,K)-Z1(I,K+1)
       ZZ1(I,K)  = .5_SP*(Z1(I,K)+Z1(I,K+1))
     END DO
   END DO

   DO I=1,NT
     ZZ1(I,KB) = 2.0_SP*ZZ1(I,KB-1)-ZZ1(I,KB-2)
   END DO

   DO K=1,KBM2
     DO I=1,NT
       DZZ1(I,K) = ZZ1(I,K)-ZZ1(I,K+1)
     END DO
   END DO
   
   DZZ1(:,KBM1) = 0.0_SP
   DZ1(:,KB)    = 0.0_SP

   IF(DBG_SET(DBG_SBR)) write(ipt,*) "END: SET_WATER_DEPTH" 
   RETURN
 END SUBROUTINE SET_WATER_DEPTH
!==============================================================================|
# if defined(ICE)
!==============================================================================!
  SUBROUTINE READ_ICE
    IMPLICIT NONE
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCDIM),  POINTER :: DIM
    LOGICAL :: FOUND
    INTEGER :: STKCNT
 
    STKCNT = NC_START%FTIME%PREV_STKCNT


    ! LOAD ice area
    VAR => FIND_VAR(NC_START,'AICEN',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'AICEN'&
         & IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, AICEN)
    CALL NC_READ_VAR(VAR,STKCNT)


    ! LOAD vicen
    VAR => FIND_VAR(NC_START,'VICEN',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'VICEN'&
         & IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, VICEN)
    CALL NC_READ_VAR(VAR,STKCNT)


    ! LOAD vsnon
    VAR => FIND_VAR(NC_START,'VSNON',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'VSNON'&
         & IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, VSNON)
    CALL NC_READ_VAR(VAR,STKCNT)


    ! LOAD Tsfcn
    VAR => FIND_VAR(NC_START,'TSFCN',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'TSFCN'&
         & IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, TSFCN)
    CALL NC_READ_VAR(VAR,STKCNT)

    ! LOAD esnon
    VAR => FIND_VAR(NC_START,'ESNON',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'ESNON'&
         & IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, ESNON)
    CALL NC_READ_VAR(VAR,STKCNT)

   ! LOAD eicen
    VAR => FIND_VAR(NC_START,'EICEN',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'EICEN'&
         & IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, EICEN)
    CALL NC_READ_VAR(VAR,STKCNT)

    ! LOAD UUICE
    VAR => FIND_VAR(NC_START,'uuice',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'uuice'&
         & IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, UICE2)
    CALL NC_READ_VAR(VAR,STKCNT)

     ! LOAD VVICE
    VAR => FIND_VAR(NC_START,'vvice',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'vvice'&
         & IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, VICE2)
    CALL NC_READ_VAR(VAR,STKCNT)
   
     ! LOAD stress1
    VAR => FIND_VAR(NC_START,'sig1',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'sig1'&
         & IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, sig1)
    CALL NC_READ_VAR(VAR,STKCNT)

     ! LOAD stress2
    VAR => FIND_VAR(NC_START,'sig2',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'sig2'&
         & IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, sig2)
    CALL NC_READ_VAR(VAR,STKCNT)

     ! LOAD stress12
    VAR => FIND_VAR(NC_START,'sig12',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'sig12'&
         & IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, sig12)
    CALL NC_READ_VAR(VAR,STKCNT)

   
 END SUBROUTINE READ_ICE
 
 
# endif

# if defined (SEDIMENT)
!-----------------------------------------------------------------------------
! Hot start the sediment model variables
!-----------------------------------------------------------------------------
  SUBROUTINE READ_SED

# if defined (ORIG_SED)
    USE MOD_SED
# elif defined (CSTMS_SED)
    USE MOD_SED_CSTMS
# endif

    IMPLICIT NONE
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCDIM),  POINTER :: DIM
    REAL(SP),SAVE, POINTER :: TMP(:),TMP2(:,:)
    LOGICAL :: FOUND
    INTEGER :: STKCNT
    INTEGER :: II

    STKCNT = NC_START%FTIME%PREV_STKCNT

    ! SEDIMENT CONCENTRATIONS
    DO II = 1,NSED
      VAR => FIND_VAR(NC_START,TRIM(sed(ii)%sname),FOUND)
      IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE "//    &
         & TRIM(sed(ii)%sname)//" IN THE HOTSTART FILE OBJECT")
      CALL NC_CONNECT_AVAR(VAR, sed(ii)%conc)
      CALL NC_READ_VAR(VAR,STKCNT)
    END DO

   ! BOTTOM VARIABLES
   DO II = 1,N_BOT_VARS
!     TMP => bottom(:,ii)
     TMP => bottom(1:m,ii)
     VAR => FIND_VAR(NC_START,TRIM(BOT_SNAMES(II)),FOUND)
     IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE "//    &
        & TRIM(BOT_SNAMES(II))//" IN THE HOTSTART FILE OBJECT")
     CALL NC_CONNECT_PVAR(VAR, TMP)
     CALL NC_READ_VAR(VAR,STKCNT)
     NULLIFY(TMP)
   END DO

  ! BED VARIABLES
   DO II = 1,N_BED_CHARS
     TMP2 => bed(1:m,1:nbed,ii)
     VAR => FIND_VAR(NC_START,TRIM(BED_SNAMES(II)),FOUND)
     IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE "//    &
        & TRIM(BED_SNAMES(II))//" IN THE HOTSTART FILE OBJECT")
     CALL NC_CONNECT_PVAR(VAR, TMP2)
     CALL NC_READ_VAR(VAR,STKCNT)
     NULLIFY(TMP2)
   END DO

  END SUBROUTINE READ_SED
# endif

# if defined (WATER_QUALITY)
!==============================================================================!
  SUBROUTINE READ_WQM
    USE MOD_WQM
    IMPLICIT NONE
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCDIM),  POINTER :: DIM
    LOGICAL :: FOUND
    INTEGER :: STKCNT, K, II
!!$    REAL(SP), DIMENSION(0:MT,KB) :: PRESSURE
    
    IF(DBG_SET(DBG_SBR)) write(ipt,*) "Start: READ_WQM" 
   
    STKCNT = NC_START%FTIME%PREV_STKCNT


    ! LOAD WATER QUALITY VARIABLES
    DO II = 1,NB
    VAR => FIND_VAR(NC_START,TRIM(WQM_NAME(II,1)),FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE "//    &
         & TRIM(WQM_NAME(II,1))//" IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, WQM)
    CALL NC_READ_VAR(VAR,STKCNT)


    ! LOAD MEAN INITIAL WATER QUALITY VARIABLES
    VAR => FIND_VAR(NC_START,TRIM(WQM_NAME(II,4)),FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE "//    &
         & TRIM(WQM_NAME(II,4))//" IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, WMEAN)
    CALL NC_READ_VAR(VAR)
    END DO

   IF(DBG_SET(DBG_SBR)) write(ipt,*) "End: READ_WQM" 
   
 END SUBROUTINE READ_WQM
!==============================================================================!
  SUBROUTINE SET_OBSERVED_WQM
    USE MOD_WQM
    IMPLICIT NONE
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCDIM),  POINTER :: DIM
    LOGICAL :: FOUND
    INTEGER :: STKCNT
    INTEGER KSL                                 !!NUMBER OF STANDARD SEA LEVELS 
    REAL(SP), ALLOCATABLE, TARGET :: DPTHSL(:)          !!DEPTH AT STANDARD SEA LEVEL
    REAL(SP), ALLOCATABLE, TARGET :: WQMSL(:,:,:) !!WQM AT STANDARD SEA LEVEL
    REAL(SP),ALLOCATABLE          :: WQMA(:,:),WQMSL_TMP(:,:)

    REAL(SP),DIMENSION(KBM1,NB)    :: WQMI
    REAL(SP),DIMENSION(KBM1)      :: ZI
    REAL(SP),DIMENSION(0:MT,KB) :: PRESSURE
    INTEGER :: K, I, IB, IERR, II
    REAL(SP) :: SCOUNT, FAC, RBUF
    
    IF(DBG_SET(DBG_SBR)) write(ipt,*) "START: SET_OBSERVED_WQM" 

    WQMI     = 0.0_SP
    ZI       = 0.0_SP
    PRESSURE = 0.0_SP

    STKCNT = NC_START%FTIME%PREV_STKCNT

    DIM => FIND_DIM(NC_START,'ksl',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND DIMENSION 'ksl'&
         & IN THE STARTUP FILE OBJECT")

    KSL = DIM%DIM
    ! ALLOCATE SPACE
    ALLOCATE(WQMSL(MT,KSL,NB))
    ALLOCATE(DPTHSL(KSL))
    WQMSL   = 0.0_SP
    DPTHSL= 0.0_SP

    ! ALLOCATE SPACE FOR THE AVERAGE WATER QUALITY VARIABLES AT EACH DEPTH
    ALLOCATE(WQMA(KSL,NB))
    WQMA = 0.0_SP

    VAR => FIND_VAR(NC_START,'zsl',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'zsl'&
         & IN THE STARTUP FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, DPTHSL)
    CALL NC_READ_VAR(VAR)
    CALL NC_DISCONNECT(VAR)

    DO II = 1,NB
    ! READ THE DATA
    ALLOCATE(WQMSL_TMP(MT,KSL))
    WQMSL_TMP(:,:) = WQMSL(:,:,II)
    VAR => FIND_VAR(NC_START,TRIM(WQM_NAME(II,1)),FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE "//      &
      & TRIM(WQM_NAME(II,1))//" IN THE STARTUP FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, WQMSL_TMP)
    CALL NC_READ_VAR(VAR,STKCNT)
    CALL NC_DISCONNECT(VAR)
    
    WQMSL(:,:,II) = WQMSL_TMP(:,:)
    DEALLOCATE(WQMSL_TMP)

!==============================================================================|
!    CALCULATE HORIZONTAL AVERAGE VALUES OF WATER QUALITY VARIABLES            !
!==============================================================================|
    
    ! GET THE AVERAGE WQM AT EACH DEPTH ON STANDARD LEVELS
    IF(SERIAL)THEN
       DO K = 1, KSL
          SCOUNT = 0.0_SP
          DO I = 1, MT
             IF(-H(I) <= DPTHSL(K)) THEN
                SCOUNT = SCOUNT + 1.0_SP
                WQMA(K,II) = WQMA(K,II) + WQMSL(I,K,II)
             END IF
          END DO
          IF(SCOUNT >= 1.0_SP)THEN
             WQMA(K,II)=WQMA(K,II)/SCOUNT
          END IF
       END DO
    END IF
    
#    if defined(MULTIPROCESSOR)
    IF(PAR)THEN
       
       DO K = 1, KSL
          SCOUNT = 0.0_SP
          DO I = 1, M
             IF(-H(I) <= DPTHSL(K)) THEN
                IF(NDE_ID(I) == 0)THEN !!INTERNAL NODE
                   SCOUNT = SCOUNT + 1.0_SP
                   WQMA(K,II) = WQMA(K,II) + WQMSL(I,K,II)
                ELSE  !!BOUNDARY NODE, ACCUMULATE FRACTION ONLY
                   DO IB = 1,NBN
                      IF(BN_LOC(IB) == I)THEN
                         FAC = 1.0_SP/FLOAT(BN_MLT(IB))
                         SCOUNT = SCOUNT + FAC 
                         WQMA(K,II) = WQMA(K,II) + FAC*WQMSL(I,K,II)
                      END IF
                   END DO
                END IF
             END IF
          END DO
          CALL MPI_ALLREDUCE(WQMA(K,II),RBUF,1,MPI_F,MPI_SUM,MPI_FVCOM_GROUP,IERR)
          WQMA(K,II) = RBUF
          CALL MPI_ALLREDUCE(SCOUNT,RBUF,1,MPI_F,MPI_SUM,MPI_FVCOM_GROUP,IERR)
          SCOUNT = RBUF

          IF(SCOUNT >= 1.0_SP)THEN
             WQMA(K,II)=WQMA(K,II)/SCOUNT
          END IF
          
       END DO !!LOOP OVER KSL
    END IF  !!PARALLEL
#    endif
    
    ! NOW INTERPOLATE FROM STANDARD LEVELS TO THE VALUE AT EACH NODE
    DO I=1,MT
       DO K=1,KBM1
          ZI(K)=ZZ(I,K)*D(I)+EL(I)
       END DO
       
       ! LEVEL AVERAGE WQM
       CALL SINTER_EXTRP_UP(DPTHSL,WQMA(:,II),ZI,WQMI(:,II),KSL,KBM1)
       
       WMEAN(I,1:KBM1,II) = WQMI(:,II)
              
       ! REAL WQMs
       CALL SINTER_EXTRP_UP(DPTHSL,WQMSL(I,:,II),ZI,WQMI(:,II),KSL,KBM1)
       
       WQM(I,1:KBM1,II) = WQMI(:,II)

    END DO
    END DO


!   CALL N2E3D(T1,T)
!   CALL N2E3D(S1,S)
!   CALL N2E3D(Tmean1,Tmean)
!   CALL N2E3D(Smean1,Smean)
!   CALL N2E3D(Rmean1,Rmean)
   ! DENSITY IS ALREADY INTERPOLATED TO ELEMENTS FOR RHO IN DENSX

   DEALLOCATE(WQMSL,DPTHSL)
   DEALLOCATE(WQMA)


   IF(DBG_SET(DBG_SBR)) write(ipt,*) "END: SET_OBSERVED_WQM" 
    
  END SUBROUTINE SET_OBSERVED_WQM
!==============================================================================!
  SUBROUTINE SET_LINEAR_WQM
    USE MOD_WQM
    IMPLICIT NONE

    ! BY DEFINITION OF LINEAR...
    INTEGER, PARAMETER :: KSL =2

    REAL(SP), DIMENSION(KBM1,NB) :: WQMI
    REAL(SP), DIMENSION(KBM1)    :: ZI
    REAL(SP), DIMENSION(KSL,NB)  :: WQMA
    REAL(SP), DIMENSION(KSL)     :: ZA
    INTEGER :: K, I,II

    IF(DBG_SET(DBG_SBR)) write(ipt,*) "START: SET_LINEAR_WQM" 

    ZA(1) = 0.0_SP
    ZA(2) = STARTUP_DMAX
    
    WQMA(1,:) = STARTUP_WQM_VALS(1,:)
    WQMA(2,:) = STARTUP_WQM_VALS(2,:)
    
    ! THE HORIZONTAL AVERAGE VALUE IS ALSO THE TRUE VALUE SINCE THE
    ! LINEAR EQUATION DOES NOT VARY WITH LOCATION 

    DO I = 1, MT
       DO K=1,KBM1
          ! CALCULATE ZI relative to z=0
          ZI(K)=ZZ(I,K)*D(I)+EL(I)
       END DO

       DO II=1,NB    
       CALL SINTER_EXTRP_UP(ZA,WQMA(:,II),ZI,WQMI(:,II),KSL,KBM1)
    
       WQM(I,1:kbm1,II) = WQMI(1:KBM1,II)
       END DO

    END DO

   ! FOR LINEAR AND CONSTANT, THE MEAN IS EQUAL TO THE INITIAL VALUE

   WMEAN = WQM
!   CALL N2E3D(T1,T)
!   TMEAN=T
   
    IF(DBG_SET(DBG_SBR)) write(ipt,*) "END: SET_LINEAR_WQM" 
    
  END SUBROUTINE SET_LINEAR_WQM
!==============================================================================!
  SUBROUTINE SET_CONSTANT_WQM
    USE MOD_WQM
    IMPLICIT NONE
    INTEGER :: I,K,II
    
    IF(DBG_SET(DBG_SBR)) write(ipt,*) "START: SET_CONSTANT_WQM" 


    DO II = 1,NB
      DO K = 1,KBM1
        DO I = 1,MT
          WQM(I,K,II) = STARTUP_WQM_VALS(1,II)
	END DO
      END DO
    END DO  	  

   ! FOR LINEAR AND CONSTANT, THE MEAN IS EQUAL TO THE INITIAL VALUE

   WMEAN = WQM

!   CALL N2E3D(T1,T)
!   CALL N2E3D(S1,S)
   
!   TMEAN=T

   IF(DBG_SET(DBG_SBR)) write(ipt,*) "END: SET_CONSTANT_WQM" 

 END SUBROUTINE SET_CONSTANT_WQM
!==============================================================================!
# endif

# if defined (BioGen)
!==============================================================================!
  SUBROUTINE READ_BIO
    IMPLICIT NONE
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCDIM),  POINTER :: DIM
    LOGICAL :: FOUND
    INTEGER :: STKCNT, K, II
    
    IF(DBG_SET(DBG_SBR)) write(ipt,*) "Start: READ_BIO" 
   
    STKCNT = NC_START%FTIME%PREV_STKCNT

    ALLOCATE(BIO_ALL(0:MT,KB,NTT))    ; BIO_ALL     =  0.001_SP
    ALLOCATE(BIO_F(0:MT,KB,NTT))      ; BIO_F       =  0.001_SP
    ALLOCATE(BIO_MEAN(00:MT,KB,NTT))  ; BIO_MEAN    =  0.001_SP
    ALLOCATE(XFLUX_OBCB(0:MT,KB,NTT)) ; XFLUX_OBCB  =  0.0_SP
    ALLOCATE(BIO_MEANN(0:NT,KB,NTT))  ; BIO_MEANN   =  0.001_SP
    ALLOCATE(BIO_VAR_MEAN(0:MT,KB,NTT)) ; BIO_VAR_MEAN   =  0.0_SP


    DO II = 1,NTT
    ! LOAD BIOLOGICAL VARIABLES
      VAR => FIND_VAR(NC_START,TRIM(BIO_NAME(II,1)),FOUND)
      IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE "//    &
         & TRIM(BIO_NAME(II,1))//" IN THE HOTSTART FILE OBJECT")
      CALL NC_CONNECT_AVAR(VAR, BIO_ALL)
      CALL NC_READ_VAR(VAR,STKCNT)


    ! LOAD MEAN INITIAL BIOLOGICAL VARIABLES
!      VAR => FIND_VAR(NC_START,TRIM(WQM_NAME(II,4)),FOUND)
!      IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE "//    &
!         & TRIM(WQM_NAME(II,4))//" IN THE HOTSTART FILE OBJECT")
!      CALL NC_CONNECT_AVAR(VAR, WMEAN)
!      CALL NC_READ_VAR(VAR)

      BIO_MEAN=BIO_ALL

    END DO

   IF(DBG_SET(DBG_SBR)) write(ipt,*) "End: READ_BIO" 
   
 END SUBROUTINE READ_BIO
!==============================================================================!

!==============================================================================!
  SUBROUTINE SET_OBSERVED_BIO
    IMPLICIT NONE
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCATT), POINTER :: ATT
    TYPE(NCDIM),  POINTER :: DIM
    LOGICAL :: FOUND
    INTEGER :: STKCNT,ii
    INTEGER KSL                                        !! NUMBER OF STANDARD SEA LEVELS 
    REAL(SP), ALLOCATABLE, TARGET :: DPTHSL(:)         !! DEPTH AT STANDARD SEA LEVEL
    REAL(SP), ALLOCATABLE, TARGET :: BIOSL_ALL(:,:,:)  !! BIO_ALL AT STANDARD SEA LEVEL
    REAL(SP), ALLOCATABLE, TARGET :: BIO_TEMP(:,:) 
    REAL(SP),DIMENSION(KBM1)      :: ZI,BioI
    REAL(SP),DIMENSION(0:MT,KB) :: PRESSURE
    INTEGER :: K, I, IB, IERR
    REAL(SP) :: SCOUNT, FAC, RBUF

    IF(DBG_SET(DBG_SBR)) write(ipt,*) "START: SET_OBSERVED_BIO" 

    ALLOCATE(BIO_ALL(0:MT,KB,NTT))    ; BIO_ALL     =  0.001_SP
    ALLOCATE(BIO_F(0:MT,KB,NTT))      ; BIO_F       =  0.001_SP
    ALLOCATE(BIO_MEAN(00:MT,KB,NTT))  ; BIO_MEAN    =  0.001_SP
    ALLOCATE(XFLUX_OBCB(0:MT,KB,NTT)) ; XFLUX_OBCB  =  0.0_SP
    ALLOCATE(BIO_MEANN(0:NT,KB,NTT))  ; BIO_MEANN   =  0.001_SP
    ALLOCATE(BIO_VAR_MEAN(0:MT,KB,NTT)) ; BIO_VAR_MEAN   =  0.0_SP

    ZI       = 0.0_SP

    STKCNT = NC_START%FTIME%PREV_STKCNT

    DIM => FIND_DIM(NC_START,'ksl',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND DIMENSION 'ksl'&
         & IN THE BIOLOGICAL PARAMETER FILE OBJECT")

    KSL = DIM%DIM
    ! ALLOCATE SPACE

    ALLOCATE(BIOSL_ALL(MT,KSL,NTT)) ;   BIOSL_ALL  =  0.0_SP
    ALLOCATE(BIO_TEMP(MT,KSL))      ;   BIO_TEMP = 0.0_SP
    ALLOCATE(DPTHSL(KSL))           ;   DPTHSL   = 0.0_SP


    ! READ THE DATA
    DO II = 1,NTT
    ! LOAD BIOLOGICAL VARIABLES
    ! LOAD BIOLOGICAL VARIABLES
      VAR => FIND_VAR(NC_START,TRIM(BIO_NAME(II,1)),FOUND)
      IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE "//    &
         & TRIM(BIO_NAME(II,1))//" IN THE BIOLOGICAL PARAMETER FILE OBJECT")

      ATT => FIND_ATT(VAR,'units',FOUND)
      IF(TRIM(ATT%CHR(1)) /= BIO_NAME(II,2))CALL FATAL_ERROR &
         & ("THE BIOLOGICAL PARAMETER FILE:"//TRIM(NC_START%FNAME),"THE VARIABLE: "&
         &//TRIM(VAR%VARNAME), "THE ATTRIBUTE 'units' IS INCORRECT: EXPE&
         &CTING "//BIO_NAME(II,2)//" and is "//TRIM(ATT%CHR(1)))

      CALL NC_CONNECT_AVAR(VAR, BIO_TEMP)
      CALL NC_READ_VAR(VAR,STKCNT)
      CALL NC_DISCONNECT(VAR)

      BIOSL_ALL(1:MT,1:KSL,II)=BIO_TEMP
    END DO


    VAR => FIND_VAR(NC_START,'zsl',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'zsl'&
         & IN THE BIOLOGICAL PARAMETER FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, DPTHSL)
    CALL NC_READ_VAR(VAR)
    CALL NC_DISCONNECT(VAR)
   
    ! NOW INTERPOLATE FROM STANDARD LEVELS TO THE VALUE AT EACH NODE
    DO II=1,NTT
      DO I=1,MT
         DO K=1,KBM1
            ZI(K)=ZZ(I,K)*D(I)+EL(I)
         END DO
       
              
         ! REAL BIO
         CALL SINTER_EXTRP_UP(DPTHSL,BIOSL_ALL(I,:,II),ZI,BIOI,KSL,KBM1)
       
         BIO_ALL(I,1:KBM1,II) = BIOI(:)

      END DO
    END DO


   DEALLOCATE(BIOSL_ALL,DPTHSL)


   IF(DBG_SET(DBG_SBR)) write(ipt,*) "END: SET_OBSERVED_BIO" 
    
  END SUBROUTINE SET_OBSERVED_BIO

# endif



# if defined (MEAN_FLOW)
 SUBROUTINE READ_MEANFLOW_RST
 USE MOD_MEANFLOW
 USE MOD_OBCS2
 USE MOD_OBCS3
 IMPLICIT NONE
 INTEGER :: INRESMF
 INTEGER :: I,K,IINT_TMP

     IF(nmfcell_i>0)THEN
       INRESMF=544
       CALL FOPEN(INRESMF,TRIM(INPUT_DIR)//TRIM(CASENAME)//'_restart_mf.dat',"cur") 
       READ(INRESMF)IINT_TMP
       IF(IINT_TMP /= IINT)THEN
         WRITE(IPT,*)'IINT IN MEANFLOW RESTART FILE NOT EQUAL TO IINT'
         WRITE(IPT,*)'FROM MAIN RESTART FILE',IINT,IINT_TMP
         CALL PSTOP
       END IF
       READ(INRESMF)(ELT   (I),   I=0,ntidenode)
       READ(INRESMF)(ELTF  (I),   I=0,ntidenode)
       READ(INRESMF)(ELRKT (I),   I=0,ntidenode)
       READ(INRESMF)(ELTDT (I),   I=0,ntidenode)
       READ(INRESMF)(ELP   (I),   I=0,ntidenode)
       READ(INRESMF)(ELPF  (I),   I=0,ntidenode)
       READ(INRESMF)(ELRKP (I),   I=0,ntidenode)
       READ(INRESMF)(UAT   (I),   I=0,ntidecell)
       READ(INRESMF)(VAT   (I),   I=0,ntidecell)
       READ(INRESMF)(UATF  (I),   I=0,ntidecell)
       READ(INRESMF)(VATF  (I),   I=0,ntidecell)
       READ(INRESMF)(UAP   (I),   I=0,ntidecell)
       READ(INRESMF)(VAP   (I),   I=0,ntidecell)
       READ(INRESMF)(UANT  (I),   I=0,  nmfcell)
       READ(INRESMF)(VANT  (I),   I=0,  nmfcell)
       READ(INRESMF)(UAN   (I),   I=0,  nmfcell)
       READ(INRESMF)(VAN   (I),   I=0,  nmfcell)
       READ(INRESMF)(UANP  (I),   I=0,  nmfcell)
       READ(INRESMF)(VANP  (I),   I=0,  nmfcell)

       READ(INRESMF)((UT (I,K),I=0,ntidecell),K=1,KBM1)
       READ(INRESMF)((VT (I,K),I=0,ntidecell),K=1,KBM1)
       READ(INRESMF)((UNT(I,K),I=0,  nmfcell),K=1,KBM1)
       READ(INRESMF)((VNT(I,K),I=0,  nmfcell),K=1,KBM1)
       READ(INRESMF)((UN (I,K),I=0,  nmfcell),K=1,KBM1)
       READ(INRESMF)((VN (I,K),I=0,  nmfcell),K=1,KBM1)

       READ(INRESMF)(UAPF  (I),  I=0,ntidecell)
       READ(INRESMF)(VAPF  (I),  I=0,ntidecell)
       READ(INRESMF)(UANTF (I),  I=0,  nmfcell)
       READ(INRESMF)(VANTF (I),  I=0,  nmfcell)
       READ(INRESMF)(UANF  (I),  I=0,  nmfcell)
       READ(INRESMF)(VANF  (I),  I=0,  nmfcell)
       READ(INRESMF)(UANPF (I),  I=0,  nmfcell)
       READ(INRESMF)(VANPF (I),  I=0,  nmfcell)
       READ(INRESMF)(UARKNT(I),  I=0,  nmfcell)
       READ(INRESMF)(VARKNT(I),  I=0,  nmfcell)
       READ(INRESMF)(UARKN (I),  I=0,  nmfcell)
       READ(INRESMF)(VARKN (I),  I=0,  nmfcell)

       READ(INRESMF)((UNTB(I,K), I=0,nmfcell),K=1,KBM1)
       READ(INRESMF)((VNTB(I,K), I=0,nmfcell),K=1,KBM1)
       READ(INRESMF)((UNB (I,K), I=0,nmfcell),K=1,KBM1)
       READ(INRESMF)((VNB (I,K), I=0,nmfcell),K=1,KBM1)
       
       READ(INRESMF)(FLUXOBN2(I),I=0,IBCN(2))
       READ(INRESMF)(CXOBC(I)   ,I=0,NT     )
       READ(INRESMF)(CYOBC(I)   ,I=0,NT     )
       READ(INRESMF)(FLUXOBC2D_X    (I)  , I=0,nmfcell_i) 
       READ(INRESMF)(FLUXOBC2D_Y    (I)  , I=0,nmfcell_i)  
       READ(INRESMF)((FLUXOBC3D_X   (I,K), I=0,nmfcell_i),K=1,KBM1)
       READ(INRESMF)((FLUXOBC3D_Y   (I,K), I=0,nmfcell_i),K=1,KBM1)    
       READ(INRESMF)((FLUXOBC3D_X_2 (I,K), I=0,nmfcell_i),K=1,KBM1)
       READ(INRESMF)((FLUXOBC3D_Y_2 (I,K), I=0,nmfcell_i),K=1,KBM1)
       READ(INRESMF)(OBC2D_X_TIDE   (I)  , I=1,  nmfcell)
       READ(INRESMF)(OBC2D_Y_TIDE   (I)  , I=1,  nmfcell)

       CLOSE(INRESMF)
       WRITE(IPT,*)'!  MEANFLOW RESTART DATA :    READ     '
     END IF

 END SUBROUTINE READ_MEANFLOW_RST
# endif

# if defined (WAVE_CURRENT_INTERACTION)
!==============================================================================!
  SUBROUTINE READ_WAVE
    USE M_GENARR
    USE SWCOMM3
    USE VARS_WAVE, ONLY : AC2,COMPDA
    
    IMPLICIT NONE
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCDIM),  POINTER :: DIM
    LOGICAL :: FOUND
    INTEGER :: STKCNT
    
    REAL(SP), ALLOCATABLE :: HS_TMP(:),TPEAK_TMP(:),WDIR_TMP(:)
    REAL(SP), ALLOCATABLE :: AC2_TMP(:)
    INTEGER :: IG
    INTEGER :: IP, IDC, ISC                                             
    REAL :: SPPARM_TMP(4)
    
    IF(DBG_SET(DBG_SBR)) write(ipt,*) "Start: READ_WAVE" 
    
    ALLOCATE(HS_TMP(0:MT))       ;HS_TMP    = 0.0_SP
    ALLOCATE(TPEAK_TMP(0:MT))    ;TPEAK_TMP = 0.0_SP
    ALLOCATE(WDIR_TMP(0:MT))     ;WDIR_TMP  = 0.0_SP
   
    STKCNT = NC_START%FTIME%PREV_STKCNT


    ! LOAD Significant Wave Height
    VAR => FIND_VAR(NC_START,'hs',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'hs'&
         & IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, HS_TMP)
    CALL NC_READ_VAR(VAR,STKCNT)


    ! LOAD Relative Peak Period
    VAR => FIND_VAR(NC_START,'tpeak',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'tpeak'&
         & IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, TPEAK_TMP)
    CALL NC_READ_VAR(VAR,STKCNT)

    ! LOAD Wave Direction
    VAR => FIND_VAR(NC_START,'wdir',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'wdir'&
         & IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, WDIR_TMP)
    CALL NC_READ_VAR(VAR,STKCNT)

    SPPARM_TMP(1:4) = SPPARM(1:4)
    
    DO IG = 1, M
     SPPARM(1) = HS_TMP(IG)
     SPPARM(2) = TPEAK_TMP(IG)                                             
     SPPARM(3) = WDIR_TMP(IG)
     SPPARM(4) = 20.0_SP
     
     CALL SSHAPE(AC2(1,1,IG),SPCSIG,SPCDIR,FSHAPE,DSHAPE)
     IF(H(IG) <= DEPMIN .OR. ABS(SPPARM(1)) < 1.0E-6_SP) AC2(:,:,IG) = 0.0_SP
    END DO 

#   if defined (MULTIPROCESSOR)
    IF(PAR)THEN
      ALLOCATE(AC2_TMP(0:MT));   AC2_TMP = 0.0
      DO ISC = 1,MSC
        DO IDC = 1,MDC
          DO IP = 1,M
            AC2_TMP(IP)=AC2(IDC,ISC,IP)
          END DO  
          CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,1,MYID,NPROCS,AC2_TMP)
          CALL AEXCHANGE(NC,MYID,NPROCS,AC2_TMP)   
          DO IP = 1,MT
            AC2(IDC,ISC,IP) = AC2_TMP(IP)
          END DO
        END DO
      END DO
      DEALLOCATE(AC2_TMP)
    END IF    	  
#   endif

    SPPARM(1:4) = SPPARM_TMP(1:4)
    DEALLOCATE(HS_TMP,TPEAK_TMP,WDIR_TMP)
    
    CALL SWANOUT(COMPDA(1,JDP2))
    
   IF(DBG_SET(DBG_SBR)) write(ipt,*) "End: READ_WAVE" 
   
 END SUBROUTINE READ_WAVE




# endif 

END MODULE MOD_STARTUP
