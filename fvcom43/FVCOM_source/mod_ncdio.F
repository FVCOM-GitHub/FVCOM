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

MODULE MOD_NCDIO
  !==============================================================================!
  !  NetCDF Io for FVCOM using CF Metadata Convention                            !
  !                                                                              !
  !    see: http://www.cgd.ucar.edu/cms/eaton/cf-metadata/ for info              !
  !                                                                              !
  !    current time dependent variables set up                                   !
  !         el:    surface elevation                                             !
  !          u:    x-velocity. In spherical coordinate,lon-velocity              !                         
  !          v:    y-velocity. In spherical coordinate,lat-velocity              !                        
  !         ww:    z-velocity                                                    !
  !         kh:    turbulent diffusivity                                         !
  !         km:    turbulent viscosity                                           !
  !         t1:    temperature                                                   !
  !         s1:    salinity                                                      !
  !         ua:    vertically-averaged x-velocity                                !
  !                In spherical coordinate,vertically-averaged lon-velocity      !
  !         va:    vertically-averaged y-velocity                                !
  !                In spherical coordinate,vertically-averaged lat-velocity      !
  !          d:    depth at procs                                                !
  !        dye:    dye at procs                                                  !
  !       aice:    ice concentration on procs                                    !
  !       vice:    ice thichness on procs                                        !
  !      uuice:    ice x-velocity                                                !
  !      vvice:    ice y-velocity                                                !
  !     uuwind:    wind speed in x direction                                     !
  !     vvwind:    wind speed in y direction                                     !
  !     pa_air:    sea level atmospheric pressure                                     !
  !                                                                              !
  !       wd:      wet/dry flag (0 or 1)                                         !
  !                                                                              !
  !       vort:    vorticity                                                     !
  !    to add additional variables:                                              !
  !      1.) add to list above                                                   !
  !      2.) add *_vid to variables vid in section "new variable vid"            !
  !      3.) go to definition section "new variable definition"                  !
  !      4.) add io section "new variable io"                                    !
  !==============================================================================!

! For Unstructred CF standard file - beta test!
!# define UCF

  USE ALL_VARS
  USE MOD_PREC
  USE MOD_NCTOOLS
  USE MOD_UTILS
  USE MOD_TIME
  USE MOD_INPUT


# if defined(DYE_RELEASE)
  USE mod_dye , only: dye_on
# endif 
  

#  if defined (SEDIMENT)

#  if defined (ORIG_SED)
  USE mod_sed
#  elif defined (CSTMS_SED)
  USE mod_sed_cstms
#  endif

! J. Ge for fluid mud layer
#  if defined (FLUID_MUD)
  Use Mod_Fluid_Mud
#  endif
! J. Ge for fluid mud layer
#  endif

#  if defined (WAVE_CURRENT_INTERACTION)    
  USE TIMECOMM   !Jianzhong
  USE SWCOMM3
  USE VARS_WAVE
  USE MOD_WAVE_CURRENT_INTERACTION
#  endif

# if defined (NH)
  USE NON_HYDRO, ONLY: W4ZT, RHS, QP, NHQDRX, NHQDRY, NHQDRZ, NHQ2DX, NHQ2DY
# endif
  implicit none

!# if defined(UCF)
!    ATT  => NC_MAKE_ATT(name='grid',values='') 
!    VAR  => ADD(VAR,ATT)
!
!    ATT  => NC_MAKE_ATT(name='grid_location',values='') 
!    VAR  => ADD(VAR,ATT)
!# endif


  LOGICAL :: VISIT_CMD_DUMP
  LOGICAL :: NCNEST_CMD_DUMP
    

  Character(LEN=50) :: CoordVar 

  LOGICAL, private :: FOUND
  logical, private :: NEED_INIT = .TRUE.

  TYPE(NCDIM), POINTER :: DIM_nele
  TYPE(NCDIM), POINTER :: DIM_node
  TYPE(NCDIM), POINTER :: DIM_three
  TYPE(NCDIM), POINTER :: DIM_four

  TYPE(NCDIM), POINTER :: DIM_siglay
  TYPE(NCDIM), POINTER :: DIM_siglev

# if defined (SEDIMENT)
  TYPE(NCDIM), POINTER :: DIM_bedlay
# endif

  TYPE(NCDIM), POINTER :: DIM_time
  TYPE(NCDIM), POINTER :: DIM_DateStrLen

  TYPE(NCDIM), POINTER :: DIM_nobc
  TYPE(NCDIM), POINTER :: DIM_nlsf

  TYPE(NCDIM), POINTER :: DIM_MaxNode
  TYPE(NCDIM), POINTER :: DIM_MaxElem

  TYPE(NCDIM), POINTER :: DIM_GRID
  TYPE(NCDIM), POINTER :: DIM_NCAT
  TYPE(NCDIM), POINTER :: DIM_ntilay

# if defined (WAVE_CURRENT_INTERACTION)    
  TYPE(NCDIM), POINTER :: DIM_MSC
  TYPE(NCDIM), POINTER :: DIM_MDC
# endif

  ! FILE POINTERS FOR NETCDF AVERAGE OUTPUT
  TYPE(NCFILE), POINTER :: NC_AVG_DATA
  TYPE(NCFILE), POINTER :: NC_AVG_SUM


  ! GRID TYPES FOR DATA AND AVERAGE OUTPUT
  TYPE(GRID), POINTER :: NC_DAT_GRIDS(:)
  TYPE(GRID), POINTER :: NC_SF_GRIDS(:)

  TYPE(GRID), POINTER :: NC_AVG_GRIDS(:)


  save
  
CONTAINS
!=============================================================  
  SUBROUTINE ARCHIVE
#   if defined(BALANCE_2D)
    USE MOD_BALANCE_2D
#   endif
    IMPLICIT NONE

    INTEGER :: STATUS
    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START ARCHIVE"


    if(NEED_INIT) then
       IF(USE_MPI_IO_MODE) THEN
          CALL MPI_IO_SYNCHRONIZE(INIT_CODE)
       ELSE
          CALL CALL_FUNC(INIT_CODE,status)
          IF (status/=0) call fatal_error("ARCHIVE:: Bad INIT_CODE",&
               & "Could not retrieve valid function pointer?")
       END IF

       NEED_INIT = .FALSE.
       VISIT_CMD_DUMP = .FALSE.
    end if

    NCNEST_CMD_DUMP = .FALSE.

    IF(NC_ON)THEN
       ! bounds checking
       !IF(NC_DAT%FTIME%NEXT_IO == IntTime .or. FORCE_ARCHIVE) THEN
       IF(abs(NC_DAT%FTIME%NEXT_IO -IntTime)<0.1_SP*IMDTI .or. FORCE_ARCHIVE) THEN
        
          IF(USE_MPI_IO_MODE) THEN
             CALL MPI_IO_SYNCHRONIZE(NC_CODE)
          ELSE
             CALL CALL_FUNC(NC_CODE,status)
             IF (status/=0) call fatal_error("ARCHIVE:: Bad NC_CODE",&
                  & "Could not retrieve valid function pointer?")
          END IF
       END IF
    END IF

    IF(NCSF_ON)THEN
       ! bounds checking
       !IF(NC_SF%FTIME%NEXT_IO == IntTime .or. FORCE_ARCHIVE) THEN
       IF(abs(NC_SF%FTIME%NEXT_IO -IntTime)<0.1_SP*IMDTI .or. FORCE_ARCHIVE) THEN
        
         U_SURFACE(:)   = U(:,1)
         V_SURFACE(:)   = V(:,1)
         T1_SURFACE(:)  = T1(:,1)
         S1_SURFACE(:)  = S1(:,1)
         VISCOFM_SURFACE(:) = VISCOFM(:,1)
         VISCOFH_SURFACE(:) = VISCOFH(:,1)
       
          IF(USE_MPI_IO_MODE) THEN
             CALL MPI_IO_SYNCHRONIZE(NCSF_CODE)
          ELSE
             CALL CALL_FUNC(NCSF_CODE,status)
             IF (status/=0) call fatal_error("ARCHIVE:: Bad NCSF_CODE",&
                  & "Could not retrieve valid function pointer?")
          END IF
       END IF
    END IF


    IF(NCAV_ON)THEN

!       !qxu IF(NC_AVG%FTIME%NEXT_IO < IntTime) CALL ADD_AVERAGE
!       IF(NC_AVG%FTIME%PREV_IO < IntTime) CALL ADD_AVERAGE
!       write(*,*) 'NC_AVG%FTIME%PREV_IO=',NC_AVG%FTIME%PREV_IO
!       write(*,*) '             IntTime=',IntTime
!       !CALL ADD_AVERAGE
!       !IF((NC_AVG%FTIME%NEXT_IO + NC_AVG%FTIME%INTERVAL)  == IntTime)THEN
!       !DAS{ ADD for AVG TIMES
!       !IF(abs(NC_AVG%FTIME%NEXT_IO + NC_AVG%FTIME%INTERVAL-IntTime)<0.1_SP*IMDTI)THEN
!       IF(abs(NC_AVG%FTIME%NEXT_IO -IntTime)<0.1_SP*IMDTI)THEN
!       !DAS}
!          CALL DIVIDE_AVERAGE
!
!          IF(USE_MPI_IO_MODE) THEN
!             CALL MPI_IO_SYNCHRONIZE(NCAV_CODE)
!          ELSE
!             CALL CALL_FUNC(NCAV_CODE,status)
!             IF (status/=0) call fatal_error("ARCHIVE:: Bad NCAV_CODE",&
!                  & "Could not retrieve valid function pointer?")
!          END IF
!
!          CALL ZERO_AVERAGE
!
!         END IF


          ! Start averaging data at Next_IO
          IF(NC_AVG%FTIME%NEXT_IO < IntTime) CALL ADD_AVERAGE
          write(*,*) 'NC_AVG%FTIME%NEXT_IO=',NC_AVG%FTIME%NEXT_IO
          write(*,*) '             IntTime=',IntTime
          
          ! When average interval is complete - write to file
          IF((NC_AVG%FTIME%NEXT_IO + NC_AVG%FTIME%INTERVAL)  == IntTime)THEN

             CALL DIVIDE_AVERAGE

             IF(USE_MPI_IO_MODE) THEN
                CALL MPI_IO_SYNCHRONIZE(NCAV_CODE)
             ELSE
                CALL CALL_FUNC(NCAV_CODE,status)
                IF (status/=0) call fatal_error("ARCHIVE:: Bad NCAV_CODE",&
                & "Could not retrieve valid function pointer?")
             END IF

             CALL ZERO_AVERAGE

       END IF
    END IF

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END AVG ARCHIVE"

    IF(RST_ON)THEN
       !IF(NC_RST%FTIME%NEXT_IO == IntTime .or.VISIT_CMD_DUMP) then
       IF(abs(NC_RST%FTIME%NEXT_IO-IntTime)<0.1_SP*IMDTI .or.VISIT_CMD_DUMP) then

          VISIT_CMD_DUMP = .false.
          
          NCNEST_CMD_DUMP = .TRUE.

          IF(USE_MPI_IO_MODE) THEN
             CALL MPI_IO_SYNCHRONIZE(RESTART_CODE)
          ELSE
             CALL CALL_FUNC(RESTART_CODE,status)
             IF (status/=0) call fatal_error("ARCHIVE:: Bad RESTART_CODE",&
                  & "Could not retrieve valid function pointer?")
# if defined (MEAN_FLOW)
             CALL MEANFLOW_RST
# endif
          END IF
       END IF
    END IF
#  if defined (BALANCE_2D)
   IF(OUT_BALANCE)THEN
     CALL OUT_TIMESERIES_BALANCE
   ENDIF
#  endif
        
    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END ARCHIVE"
  END SUBROUTINE ARCHIVE
!=============================================================  
  SUBROUTINE INIT_NCDIO
    IMPLICIT NONE

    IF(DBG_SET(DBG_LOG)) THEN
       
       write(IPT,*)"!=============================================================="
       write(IPT,*)"! SETTING UP NCDIO: CREATING AND DUMPING OUTPUT FILE META DATA"
       write(IPT,*)"!=============================================================="
    END IF

    ! DEFINE DIMENSIONS HERE - THEY ARE KILLED AT THE END OF THIS
    ! SETUP SCRIPT...

#  if defined (SPHERICAL)
    CoordVar="lat lon"
# else
    CoordVar="x y"
#endif


    IF(NC_ON) then
       
       IF(.not. ASSOCIATED(NC_DAT)) Call Fatal_Error &
            & ("INIT_NCDIO: THE DATA FILE OBJECT IS NOT ASSOCIATED ")
           
       CALL SETUP_DATFILE


    END IF
  
    IF(NCSF_ON) then
       
       IF(.not. ASSOCIATED(NC_SF)) Call Fatal_Error &
            & ("INIT_NCDIO_SURFACE: THE DATA FILE OBJECT IS NOT ASSOCIATED ")
           
       CALL SETUP_SFFILE


    END IF
  
    IF(NCAV_ON) then
       
       IF(.not. ASSOCIATED(NC_AVG)) Call Fatal_Error &
            & ("INIT_NCDIO: THE AVERAGE FILE OBJECT IS NOT ASSOCIATED ")


       CALL SETUP_AVGFILE
    END IF

    IF(RST_ON) then


       IF(.not. ASSOCIATED(NC_RST)) Call Fatal_Error &
            & ("INIT_NCDIO: THE RESTART FILE OBJECT IS NOT ASSOCIATED ")

       CALL SETUP_RSTFILE
    END IF

    IF(DBG_SET(DBG_LOG)) THEN
       

       write(IPT,*)"! FINISHED NCDIO SETUP!"
       write(IPT,*)"!=============================================================="
    END IF


  END SUBROUTINE INIT_NCDIO
!=============================================================  
! NEED AN INTERFACE WITH NO ARGS FOR FUNCTION POINTERS
!=============================================================  
  SUBROUTINE DUMP_NC_DAT
    IMPLICIT NONE
    INTEGER :: I
    TYPE(NCFILE), POINTER :: NCF
    LOGICAL :: FOUND


    IF (ICING_MODEL .AND. .NOT. IOPROC) CALL ICING(IntTime)

    DO I = 1, SIZE(NC_DAT_GRIDS)
       
       NCF => FIND_FILE(FILEHEAD,NC_DAT_GRIDS(I)%NAME,FOUND)
       IF(.NOT.FOUND) CALL FATAL_ERROR&
            ("DUMP_NC_DAT: CAN NOT FILE FILE OBJECT NAME:"//TRIM(NC_DAT_GRIDS(I)%NAME))

       CALL DUMP_DATA(NCF)

       ! INCASE THE NAME WAS CHANGED INSIDE DUMP_DATA
       NC_DAT_GRIDS(I)%NAME = NCF%FNAME

    END DO

  END SUBROUTINE DUMP_NC_DAT
  !=============================================================  
  SUBROUTINE DUMP_NC_SF
    IMPLICIT NONE
    INTEGER :: I
    TYPE(NCFILE), POINTER :: NCF
    LOGICAL :: FOUND


    IF (ICING_MODEL .AND. .NOT. IOPROC) CALL ICING(IntTime)

    DO I = 1, SIZE(NC_SF_GRIDS)
       
       NCF => FIND_FILE(FILEHEAD,NC_SF_GRIDS(I)%NAME,FOUND)
       IF(.NOT.FOUND) CALL FATAL_ERROR&
            ("DUMP_NC_SF: CAN NOT FILE FILE OBJECT NAME:"//TRIM(NC_SF_GRIDS(I)%NAME))

       CALL DUMP_DATA(NCF)

       ! INCASE THE NAME WAS CHANGED INSIDE DUMP_DATA
       NC_SF_GRIDS(I)%NAME = NCF%FNAME

    END DO

  END SUBROUTINE DUMP_NC_SF
  !=============================================================  

  SUBROUTINE DUMP_NC_AVG
    IMPLICIT NONE
    INTEGER :: I
    TYPE(NCFILE), POINTER :: NCF
    LOGICAL :: FOUND
    
    DO I = 1, SIZE(NC_AVG_GRIDS)
       
       NCF => FIND_FILE(FILEHEAD,NC_AVG_GRIDS(I)%NAME,FOUND)
       IF(.NOT.FOUND) CALL FATAL_ERROR&
            ("DUMP_NC_AVG: CAN NOT FILE FILE OBJECT NAME:"//TRIM(NC_AVG_GRIDS(I)%NAME))

       CALL DUMP_DATA(NCF)

       ! INCASE THE NAME WAS CHANGED INSIDE DUMP_DATA
       NC_AVG_GRIDS(I)%NAME = NCF%FNAME

    END DO

  END SUBROUTINE DUMP_NC_AVG
  !=============================================================  
  SUBROUTINE DUMP_NC_RST
    IMPLICIT NONE
    CALL DUMP_DATA(NC_RST)
  END SUBROUTINE DUMP_NC_RST
!=============================================================  
  SUBROUTINE SETUP_DATFILE
    IMPLICIT NONE
    character(len=80) :: tmp,dat_name
    INTEGER :: NUMG, I

    TYPE(NCFILE), POINTER :: NCF,NCF_TMP,NCF2

    LOGICAL :: INCLUDE_MASTER = .false.

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START SETUP_DATAFILE"

    IF(DBG_SET(DBG_LOG)) THEN
       
       write(IPT,*)"!--------------------------------------------------"
       write(IPT,*)"! SETTING UP DATA FILE OUTPUTS..."
    END IF


    ! Get the suffix from the NC_DAT file
    I = len_trim(NC_DAT%FNAME)
    dat_name = NC_DAT%FNAME(I-7:I)
    

    CALL SETUP_SUBDOMAINS(NC_SUBDOMAIN_FILES,NC_DAT_GRIDS)

    NUMG = size(NC_DAT_GRIDS)

    DO I = 1, NUMG

       ! DEFINEN DIMENSIONS FOR THIS GRID
       CALL DEFINE_DIMENSIONS(NC_DAT_GRIDS(I))
       
       IF(NC_DAT_GRIDS(I)%NAME /="FVCOM") THEN
          ! MAKE NEW FILE OBJECT
          tmp = trim(OUTPUT_DIR)//TRIM(NC_DAT_GRIDS(I)%NAME)//TRIM(DAT_NAME)
          NCF => NEW_FILE(TRIM(TMP))
          
          ! SET THE CURRENT NAME IN THE GRID OBJECT
          NC_DAT_GRIDS(I)%NAME=NCF%FNAME

          ! SET THE FTIME OBJECT
          NCF%FTIME => NEW_FTIME()
          NCF%FTIME = NC_DAT%FTIME
          
          ! MAKE A TEMPORARY POINTER TO ADD THE FILE TO THE LIST
          NCF_TMP => NCF
          FILEHEAD => ADD(FILEHEAD,NCF_TMP)

       ELSE
          ! THIS IS THE FVCOM GRID FILE
          NCF => NC_DAT
          NC_DAT_GRIDS(I)%NAME = NC_DAT%FNAME
          INCLUDE_MASTER = .TRUE.
       END IF




       ! ADD THE DATA OBJECTS
       NCF2 => GRID_FILE_OBJECT(NC_DAT_GRIDS(I))
       NCF => ADD(NCF,NCF2)
!!$       NCF => ADD(NCF,GRID_FILE_OBJECT(NC_DAT_GRIDS(I)) )
       
       NCF2 => TIME_FILE_OBJECT()
       NCF => ADD(NCF,NCF2)
!!$       NCF => ADD(NCF,TIME_FILE_OBJECT() )
       
       NCF2 => ZETA_FILE_OBJECT()
       NCF => ADD(NCF,NCF2)
!!$       NCF => ADD(NCF,ZETA_FILE_OBJECT() )
       
       IF(NC_FILE_DATE) THEN
          NCF2 => FILE_DATE_OBJECT()
          NCF => ADD(NCF,NCF2)
!!$          NCF => ADD(NCF,FILE_DATE_OBJECT() )
       END IF
       
       IF(NC_GRID_METRICS) THEN
          NCF2 => GRID_METRICS_FILE_OBJECT(NC_DAT_GRIDS(I))
          NCF => ADD(NCF,NCF2)
!!$          NCF => ADD(NCF,GRID_METRICS_FILE_OBJECT(NC_DAT_GRIDS(I)) )
       END IF
       
       IF(NC_VELOCITY) THEN
          NCF2 => VELOCITY_FILE_OBJECT()
          NCF => ADD(NCF,NCF2)
!!$          NCF => ADD(NCF,VELOCITY_FILE_OBJECT() )
       END IF
       
       IF(NC_VERTICAL_VEL) THEN
          NCF2 => VERTICAL_VEL_FILE_OBJECT()
          NCF => ADD(NCF,NCF2)
!!$          NCF => ADD(NCF,VERTICAL_VEL_FILE_OBJECT() )
       END IF
       
       IF(NC_AVERAGE_VEL) THEN
          NCF2 => AVERAGE_VEL_FILE_OBJECT()
          NCF => ADD(NCF,NCF2)
!!$          NCF => ADD(NCF,AVERAGE_VEL_FILE_OBJECT() )
       END IF
       
       IF(NC_VORTICITY) THEN
          NCF2 => VORTICITY_FILE_OBJECT()
          NCF => ADD(NCF,NCF2)
!!$          NCF => ADD(NCF,VORTICITY_FILE_OBJECT() )
       END IF
       
       IF(NC_SALT_TEMP) THEN
          NCF2 => SALT_TEMP_FILE_OBJECT()
          NCF => ADD(NCF,NCF2)
!!$          NCF => ADD(NCF,SALT_TEMP_FILE_OBJECT() )
       END IF
       
       IF(NC_TURBULENCE) THEN
          NCF2 => TURBULENCE_FILE_OBJECT()
          NCF => ADD(NCF,NCF2)
!!$          NCF => ADD(NCF,TURBULENCE_FILE_OBJECT() )
       END IF
       
       IF (NC_SURFACE_HEAT .and. HEATING_ON) THEN
          NCF2 => SURFACE_HEATING_FILE_OBJECT()
          NCF => ADD(NCF,NCF2)
!!$          NCF => ADD(NCF,SURFACE_HEATING_FILE_OBJECT() )
       END IF
       
       IF (NC_WIND_VEL) THEN
          NCF2 => WIND_VELOCITY_FILE_OBJECT()
          NCF => ADD(NCF,NCF2)
!!$          NCF => ADD(NCF,WIND_VELOCITY_FILE_OBJECT() )
       END IF
       
       IF (NC_WIND_STRESS .and. WIND_ON) THEN
          NCF2 => WIND_STRESS_FILE_OBJECT()
          NCF => ADD(NCF,NCF2)
!!$          NCF => ADD(NCF,WIND_STRESS_FILE_OBJECT() )
       END IF

#      if defined (AIR_PRESSURE) || (HEATING_CALCULATED)
       IF (NC_ATM_PRESS) THEN
          NCF2 => ATMOSPHERIC_PRESSURE_FILE_OBJECT()
          NCF => ADD(NCF,NCF2)
!!$          NCF => ADD(NCF,ATMOSPHERIC_PRESSURE_FILE_OBJECT() )
       END IF
#      endif       
       
       IF (NC_EVAP_PRECIP .and. PRECIPITATION_ON) THEN
          NCF2 => PRECIPITATION_FILE_OBJECT()
          NCF => ADD(NCF,NCF2)
!!$          NCF => ADD(NCF,PRECIPITATION_FILE_OBJECT() )
       END IF
       
       IF(WETTING_DRYING_ON) THEN
          NCF2 => WET_DRY_FILE_OBJECT()
          NCF => ADD(NCF, NCF2)
!!$          NCF => ADD(NCF, WET_DRY_FILE_OBJECT() )
       END IF
       
       IF(ICING_MODEL) THEN
          NCF2 => ICING_FILE_OBJECT()
          NCF => ADD(NCF, NCF2)
!!$          NCF => ADD(NCF, ICING_FILE_OBJECT() )
       END IF
       
       IF (GROUNDWATER_ON .and. NC_GROUNDWATER) THEN
          NCF2 => GROUNDWATER_FILE_OBJECT()
          NCF => ADD(NCF,NCF2)
!!$          NCF => ADD(NCF,GROUNDWATER_FILE_OBJECT() )
       END IF
       
#   if defined (BioGen)
       IF(BIOLOGICAL_MODEL)THEN
       IF(NC_BIO)THEN
          NCF2 => BIO_FILE_OBJECT()
          NCF => ADD(NCF,NCF2)
!!$          NCF => ADD(NCF,BIO_FILE_OBJECT() )
       END IF
       END IF
#   endif      
       
#   if defined (WATER_QUALITY)
       IF(NC_WQM)THEN
          NCF2 => WQM_FILE_OBJECT()
          NCF => ADD(NCF,NCF2)
!!$          NCF => ADD(NCF,WQM_FILE_OBJECT() )
       END IF
#   endif      
       
#   if defined (SEDIMENT) 
    IF(SEDIMENT_MODEL)THEN
      NCF2 => SED_FILE_OBJECT()
      NCF => ADD(NCF,NCF2)
!!$      NCF => ADD(NCF,SED_FILE_OBJECT() )

! J. Ge for fluid mud layer
#   if defined (FLUID_MUD) 
      NCF2 => FLUID_MUD_FILE_OBJECT()
      NCF => ADD(NCF,NCF2)
!!$      NCF => ADD(NCF,SED_FILE_OBJECT() )
#   endif 
! J. Ge for fluid mud layer

    END IF

#   endif      

#   if defined (DYE_RELEASE)
    IF(DYE_ON)THEN
      NCF2 => DYE_FILE_OBJECT()
      NCF => ADD(NCF,NCF2)
!!$      NCF => ADD(NCF,DYE_FILE_OBJECT())
    ENDIF
#   endif      

#   if defined (WAVE_CURRENT_INTERACTION) 
      IF(NC_WAVE_PARA)THEN
        NCF2 => WAVE_PARA_FILE_OBJECT()
        NCF => ADD(NCF,NCF2)   
!!$        NCF => ADD(NCF,WAVE_PARA_FILE_OBJECT() )   

        NCF2 => WAVE_PARA1_FILE_OBJECT()
        NCF => ADD(NCF,NCF2)   
!!$        NCF => ADD(NCF,WAVE_PARA1_FILE_OBJECT() )   
     END IF
      IF(NC_WAVE_STRESS)THEN
        NCF2 => WAVE_STRESS_FILE_OBJECT()
        NCF => ADD(NCF,NCF2)  
!!$        NCF => ADD(NCF,WAVE_STRESS_FILE_OBJECT() )   !Jianzhong
      END IF
#   endif    

#   if defined (NH)
       IF(NC_NH_QP) THEN
          NCF2 => QP_FILE_OBJECT()
          NCF => ADD(NCF,NCF2)
!!$          NCF => ADD(NCF,QP_FILE_OBJECT() )
       END IF
       IF(NC_NH_RHS) THEN
          NCF2 => RHS_FILE_OBJECT()
          NCF => ADD(NCF,NCF2)
!!$          NCF => ADD(NCF,RHS_FILE_OBJECT() )
       END IF
#   endif      
       
#   if defined (ICE)
       IF(NC_ICE) THEN
          !-----------------------------------------------------------------
          ! state variables
          NCF2 => ICE_DATA_STATE_FILE_OBJECT()
          NCF => ADD(NCF,NCF2)
!!$          NCF => ADD(NCF,ICE_DATA_STATE_FILE_OBJECT() )
          !-----------------------------------------------------------------
          ! velocity
          !-----------------------------------------------------------------
          NCF2 => ICE_VEL_FILE_OBJECT()
          NCF => ADD(NCF,NCF2)
!!$          NCF => ADD(NCF,ICE_VEL_FILE_OBJECT() )
          !-----------------------------------------------------------------
       END IF
# endif
       
       IF (STARTUP_TYPE /= "crashrestart") THEN

          ! IF CRASH RESTART - OVER RIDE THE STKCNT
          NCF%FTIME%NEXT_STKCNT = 0
          CALL NC_WRITE_FILE(NCF)
          NCF%FTIME%NEXT_STKCNT = 1
          
       ELSE
          NCF%CONNECTED = .TRUE.
          NCF%WRITABLE = .TRUE.
       END IF


       CALL KILL_DIMENSIONS
    END DO

    IF(.not. INCLUDE_MASTER) THEN

       CALL KILL_FILE(NC_DAT)
       NC_DAT => FIND_FILE(FILEHEAD,NC_DAT_GRIDS(1)%NAME,FOUND)

       IF(.NOT. FOUND) CALL FATAL_ERROR&
            &("LOGICAL ERROR IN SETUP_DATFILE")
    END IF

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END SETUP_DATAFILE"
  END SUBROUTINE SETUP_DATFILE
!=============================================================  
  SUBROUTINE SETUP_SFFILE
    IMPLICIT NONE
    character(len=80) :: tmp,dat_name
    INTEGER :: NUMG, I

    TYPE(NCFILE), POINTER :: NCF,NCF_TMP,NCF2

    LOGICAL :: INCLUDE_MASTER = .false.

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START SETUP_SFFILE"

    IF(DBG_SET(DBG_LOG)) THEN
       
       write(IPT,*)"!--------------------------------------------------"
       write(IPT,*)"! SETTING UP SURFACE FILE OUTPUTS..."
    END IF


    ! Get the suffix from the NC_SF file
    I = len_trim(NC_SF%FNAME)
    dat_name = NC_SF%FNAME(I-7:I)
    

    CALL SETUP_SUBDOMAINS(NCSF_SUBDOMAIN_FILES,NC_SF_GRIDS)

    NUMG = size(NC_SF_GRIDS)

    DO I = 1, NUMG

       ! DEFINEN DIMENSIONS FOR THIS GRID
       CALL DEFINE_DIMENSIONS_SURFACE(NC_SF_GRIDS(I))
       
       IF(NC_SF_GRIDS(I)%NAME /="FVCOM") THEN
          ! MAKE NEW FILE OBJECT
          tmp = trim(OUTPUT_DIR)//TRIM(NC_SF_GRIDS(I)%NAME)//TRIM(DAT_NAME)
          NCF => NEW_FILE(TRIM(TMP))
          
          ! SET THE CURRENT NAME IN THE GRID OBJECT
          NC_SF_GRIDS(I)%NAME=NCF%FNAME

          ! SET THE FTIME OBJECT
          NCF%FTIME => NEW_FTIME()
          NCF%FTIME = NC_SF%FTIME
          
          ! MAKE A TEMPORARY POINTER TO ADD THE FILE TO THE LIST
          NCF_TMP => NCF
          FILEHEAD => ADD(FILEHEAD,NCF_TMP)

       ELSE
          ! THIS IS THE FVCOM GRID FILE
          NCF => NC_SF
          NC_SF_GRIDS(I)%NAME = NC_SF%FNAME
          INCLUDE_MASTER = .TRUE.
       END IF




       ! ADD THE DATA OBJECTS
       NCF2 => GRID_FILE_OBJECT_SURFACE(NC_SF_GRIDS(I))
       NCF => ADD(NCF,NCF2)
       
       NCF2 => TIME_FILE_OBJECT()
       NCF => ADD(NCF,NCF2)
       
       NCF2 => ZETA_FILE_OBJECT()
       NCF => ADD(NCF,NCF2)
       
       IF(NCSF_FILE_DATE) THEN
          NCF2 => FILE_DATE_OBJECT()
          NCF => ADD(NCF,NCF2)
       END IF
       
       IF(NCSF_GRID_METRICS) THEN
          NCF2 => GRID_METRICS_FILE_OBJECT(NC_SF_GRIDS(I))
          NCF => ADD(NCF,NCF2)
       END IF
       
       IF(NCSF_VELOCITY) THEN
          NCF2 => VELOCITY_FILE_OBJECT_SURFACE()
          NCF => ADD(NCF,NCF2)
       END IF
       
       IF(NCSF_SALT_TEMP) THEN
          NCF2 => SALT_TEMP_FILE_OBJECT_SURFACE()
          NCF => ADD(NCF,NCF2)
       END IF
       
       IF(NCSF_TURBULENCE) THEN
          NCF2 => TURBULENCE_FILE_OBJECT_SURFACE()
          NCF => ADD(NCF,NCF2)
       END IF
       
       IF (NCSF_SURFACE_HEAT .and. HEATING_ON) THEN
          NCF2 => SURFACE_HEATING_FILE_OBJECT()
          NCF => ADD(NCF,NCF2)
       END IF
       
       IF (NCSF_WIND_VEL) THEN
          NCF2 => WIND_VELOCITY_FILE_OBJECT()
          NCF => ADD(NCF,NCF2)
       END IF
       
       IF (NCSF_WIND_STRESS .and. WIND_ON) THEN
          NCF2 => WIND_STRESS_FILE_OBJECT()
          NCF => ADD(NCF,NCF2)
       END IF

#      if defined (AIR_PRESSURE) || (HEATING_CALCULATED)
       IF (NCSF_ATM_PRESS) THEN
          NCF2 => ATMOSPHERIC_PRESSURE_FILE_OBJECT()
          NCF => ADD(NCF,NCF2)
       END IF
#      endif       
       
       IF (NCSF_EVAP_PRECIP .and. PRECIPITATION_ON) THEN
          NCF2 => PRECIPITATION_FILE_OBJECT()
          NCF => ADD(NCF,NCF2)
       END IF
       
       IF(WETTING_DRYING_ON) THEN
          NCF2 => WET_DRY_FILE_OBJECT()
          NCF => ADD(NCF, NCF2)
       END IF
       
       IF(ICING_MODEL) THEN
          NCF2 => ICING_FILE_OBJECT()
          NCF => ADD(NCF, NCF2)
       END IF
       
#   if defined (WAVE_CURRENT_INTERACTION) 
      IF(NCSF_WAVE_PARA)THEN
        NCF2 => WAVE_PARA_FILE_OBJECT()
        NCF => ADD(NCF,NCF2)   
      END IF
#   endif    

#   if defined (ICE)
       IF(NCSF_ICE) THEN
          !-----------------------------------------------------------------
          ! state variables
          NCF2 => ICE_DATA_STATE_FILE_OBJECT()
          NCF => ADD(NCF,NCF2)
          !-----------------------------------------------------------------
          ! velocity
          !-----------------------------------------------------------------
          NCF2 => ICE_VEL_FILE_OBJECT()
          NCF => ADD(NCF,NCF2)
          !-----------------------------------------------------------------
       END IF
# endif
       
       IF (STARTUP_TYPE /= "crashrestart") THEN

          ! IF CRASH RESTART - OVER RIDE THE STKCNT
          NCF%FTIME%NEXT_STKCNT = 0
          CALL NC_WRITE_FILE(NCF)
          NCF%FTIME%NEXT_STKCNT = 1
          
       ELSE
          NCF%CONNECTED = .TRUE.
          NCF%WRITABLE = .TRUE.
       END IF


       CALL KILL_DIMENSIONS_SURFACE
    END DO

    IF(.not. INCLUDE_MASTER) THEN

       CALL KILL_FILE(NC_SF)
       NC_SF => FIND_FILE(FILEHEAD,NC_SF_GRIDS(1)%NAME,FOUND)

       IF(.NOT. FOUND) CALL FATAL_ERROR&
            &("LOGICAL ERROR IN SETUP_SFFILE")
    END IF

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END SETUP_SFFILE"
  END SUBROUTINE SETUP_SFFILE
!=============================================================  
  SUBROUTINE SETUP_AVGFILE
    IMPLICIT NONE
    character(len=80) :: tmp,dat_name
    TYPE(NCVAR),POINTER :: VAR
    TYPE(NCATT),POINTER :: ATT
    TYPE(NCDIM),POINTER :: DIM
    LOGICAL :: FOUND
    INTEGER :: NUMG, I
    TYPE(GRID), SAVE :: MYGRID
    TYPE(NCFILE), POINTER :: NCF,NCF_TMP
    LOGICAL :: INCLUDE_MASTER = .false.


    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START SETUP_AVGFILE"

    !===================================================================
    ! MAKE A FILE TO HOLD POINTERS TO ALL VARIABLES THAT WILL BE AVERAGED
    !===================================================================

    IF(DBG_SET(DBG_LOG)) THEN
       
       write(IPT,*)"!--------------------------------------------------"
       write(IPT,*)"! SETTING UP AVERAGE FILE OUTPUTS..."
    END IF

    CALL SET_FVCOM_GRID(MYGRID)

    CALL DEFINE_DIMENSIONS(MYGRID)

    NC_AVG_DATA => NEW_FILE()

    NCF_TMP => ZETA_FILE_OBJECT()
    NC_AVG_DATA => ADD(NC_AVG_DATA,NCF_TMP)
!!$    NC_AVG_DATA => ADD(NC_AVG_DATA,ZETA_FILE_OBJECT() )

    ! NOW ADD THE OTHER DATA VARIABLES TO THE AVG FILE
    IF(NCAV_VELOCITY) THEN
       NCF_TMP => VELOCITY_FILE_OBJECT()
       NC_AVG_DATA => ADD(NC_AVG_DATA,NCF_TMP)
!!$       NC_AVG_DATA => ADD(NC_AVG_DATA,VELOCITY_FILE_OBJECT() )
    END IF

    IF(NCAV_VERTICAL_VEL) THEN
       NCF_TMP => VERTICAL_VEL_FILE_OBJECT()
       NC_AVG_DATA => ADD(NC_AVG_DATA,NCF_TMP)
!!$       NC_AVG_DATA => ADD(NC_AVG_DATA,VERTICAL_VEL_FILE_OBJECT() )
    END IF

    IF(NCAV_AVERAGE_VEL) THEN
       NCF_TMP => AVERAGE_VEL_FILE_OBJECT()
       NC_AVG_DATA => ADD(NC_AVG_DATA,NCF_TMP)
!!$       NC_AVG_DATA => ADD(NC_AVG_DATA,AVERAGE_VEL_FILE_OBJECT() )
    END IF
    
    IF(NCAV_VORTICITY) THEN
       NCF_TMP => VORTICITY_FILE_OBJECT()
       NC_AVG_DATA => ADD(NC_AVG_DATA,NCF_TMP)
!!$       NC_AVG_DATA => ADD(NC_AVG_DATA,VORTICITY_FILE_OBJECT() )
    END IF
    
    IF(NCAV_SALT_TEMP) THEN
       NCF_TMP => SALT_TEMP_FILE_OBJECT()
       NC_AVG_DATA => ADD(NC_AVG_DATA,NCF_TMP)
!!$       NC_AVG_DATA => ADD(NC_AVG_DATA,SALT_TEMP_FILE_OBJECT() )
    END IF

    IF(NCAV_TURBULENCE) THEN
       NCF_TMP => TURBULENCE_FILE_OBJECT()
       NC_AVG_DATA => ADD(NC_AVG_DATA,NCF_TMP)
!!$       NC_AVG_DATA => ADD(NC_AVG_DATA,TURBULENCE_FILE_OBJECT() )
    END IF

    IF (NCAV_SURFACE_HEAT .and. HEATING_ON) THEN
       NCF_TMP => SURFACE_HEATING_FILE_OBJECT()
       NC_AVG_DATA => ADD(NC_AVG_DATA,NCF_TMP)
 !!$      NC_AVG_DATA => ADD(NC_AVG_DATA,SURFACE_HEATING_FILE_OBJECT() )
    END IF

    IF (NCAV_WIND_VEL) THEN
       NCF_TMP => WIND_VELOCITY_FILE_OBJECT()
       NC_AVG_DATA => ADD(NC_AVG_DATA,NCF_TMP)
!!$       NC_AVG_DATA => ADD(NC_AVG_DATA,WIND_VELOCITY_FILE_OBJECT() )
    END IF

    IF (NCAV_WIND_STRESS .and. WIND_ON) THEN
       NCF_TMP => WIND_STRESS_FILE_OBJECT()
       NC_AVG_DATA => ADD(NC_AVG_DATA,NCF_TMP)
!!$       NC_AVG_DATA => ADD(NC_AVG_DATA,WIND_STRESS_FILE_OBJECT() )
    END IF

#   if defined (AIR_PRESSURE) || (HEATING_CALCULATED)
    IF (NCAV_ATM_PRESS) THEN
       NCF_TMP => ATMOSPHERIC_PRESSURE_FILE_OBJECT()
       NC_AVG_DATA => ADD(NC_AVG_DATA,NCF_TMP)
!!$       NC_AVG_DATA => ADD(NC_AVG_DATA,ATMOSPHERIC_PRESSURE_FILE_OBJECT() )
    END IF
#   endif    

    IF (NCAV_EVAP_PRECIP .and. PRECIPITATION_ON) THEN
       NCF_TMP => PRECIPITATION_FILE_OBJECT()
       NC_AVG_DATA => ADD(NC_AVG_DATA,NCF_TMP)
!!$       NC_AVG_DATA => ADD(NC_AVG_DATA,PRECIPITATION_FILE_OBJECT() )
    END IF

!    IF(WETTING_DRYING_ON) THEN
!       NC_AVG_DATA => ADD(NC_AVG_DATA, WET_DRY_FILE_OBJECT() )
!    END IF

!    IF(ICING_MODEL) THEN
!       NC_AVG_DATA => ADD(NC_AVG_DATA, ICING_FILE_OBJECT() )
!    END IF

    IF (GROUNDWATER_ON .and. NCAV_GROUNDWATER) THEN
       NCF_TMP => GROUNDWATER_FILE_OBJECT()
       NC_AVG_DATA => ADD(NC_AVG_DATA,NCF_TMP)
!!$       NC_AVG_DATA => ADD(NC_AVG_DATA,GROUNDWATER_FILE_OBJECT() )
    END IF

#   if defined (BioGen)
    IF(BIOLOGICAL_MODEL)THEN
    IF(NCAV_BIO)THEN
      NCF_TMP => BIO_FILE_OBJECT()
      NC_AVG_DATA => ADD(NC_AVG_DATA,NCF_TMP)
!!$      NC_AVG_DATA => ADD(NC_AVG_DATA,BIO_FILE_OBJECT() )
    END IF
    END IF
#   endif      

#   if defined (WATER_QUALITY)
    IF(NCAV_WQM)THEN
      NCF_TMP => WQM_FILE_OBJECT()
      NC_AVG_DATA => ADD(NC_AVG_DATA,NCF_TMP)
!!$      NC_AVG_DATA => ADD(NC_AVG_DATA,WQM_FILE_OBJECT() )
    END IF
#   endif      

#   if defined (SEDIMENT)
    IF(SEDIMENT_MODEL)THEN
      NCF_TMP => SED_FILE_OBJECT()
      NC_AVG_DATA => ADD(NC_AVG_DATA,NCF_TMP)
!!$      NC_AVG_DATA => ADD(NC_AVG_DATA,SED_FILE_OBJECT() )
! J. Ge for fluid mud layer
#   if defined (FLUID_MUD)
      NCF_TMP => FLUID_MUD_FILE_OBJECT()
      NC_AVG_DATA => ADD(NC_AVG_DATA,NCF_TMP)
!!$      NC_AVG_DATA => ADD(NC_AVG_DATA,SED_FILE_OBJECT() )
#   endif
! J. Ge for fluid mud layer
    ENDIF
#   endif      

#   if defined (WAVE_CURRENT_INTERACTION)
    IF(NCAV_WAVE_PARA)THEN 
      NCF_TMP => WAVE_PARA_FILE_OBJECT()
      NC_AVG_DATA => ADD(NC_AVG_DATA,NCF_TMP)
!!$      NC_AVG_DATA => ADD(NC_AVG_DATA,WAVE_PARA_FILE_OBJECT() )

      NCF_TMP => WAVE_PARA1_FILE_OBJECT()
      NC_AVG_DATA => ADD(NC_AVG_DATA,NCF_TMP)
!!$      NC_AVG_DATA => ADD(NC_AVG_DATA,WAVE_PARA_FILE_OBJECT() )
    END IF
    IF(NCAV_WAVE_STRESS)THEN 
      NCF_TMP => WAVE_STRESS_FILE_OBJECT()
      NC_AVG_DATA => ADD(NC_AVG_DATA,NCF_TMP)
!!$      NC_AVG_DATA => ADD(NC_AVG_DATA,WAVE_STRESS_FILE_OBJECT() )
    END IF
#   endif   

#   if defined (DYE_RELEASE)
    IF(DYE_ON)THEN
      NCF_TMP => DYE_FILE_OBJECT()
      NC_AVG_DATA => ADD(NC_AVG_DATA,NCF_TMP)
!!$      NC_AVG_DATA => ADD(NC_AVG_DATA,DYE_FILE_OBJECT())
    ENDIF
#   endif      

#   if defined (NH)
    IF(NCAV_NH_QP) THEN
      NCF_TMP => QP_FILE_OBJECT()
      NC_AVG_DATA => ADD(NC_AVG_DATA,NCF_TMP)
!!$      NC_AVG_DATA => ADD(NC_AVG_DATA,QP_FILE_OBJECT() )
    END IF
    IF(NCAV_NH_RHS) THEN
      NCF_TMP => RHS_FILE_OBJECT()
      NC_AVG_DATA => ADD(NC_AVG_DATA,NCF_TMP)
!!$      NC_AVG_DATA => ADD(NC_AVG_DATA,RHS_FILE_OBJECT() )
    END IF
#   endif      
    
#   if defined (ICE)
    IF(NCAV_ICE) THEN
       !-----------------------------------------------------------------
       ! state variables
       NCF_TMP => ICE_DATA_STATE_FILE_OBJECT()
       NC_AVG_DATA => ADD(NC_AVG_DATA,NCF_TMP)
!!$       NC_AVG_DATA => ADD(NC_AVG_DATA,ICE_DATA_STATE_FILE_OBJECT() )
       !-----------------------------------------------------------------
       ! velocity
       !-----------------------------------------------------------------
       NCF_TMP => ICE_VEL_FILE_OBJECT()
       NC_AVG_DATA => ADD(NC_AVG_DATA,NCF_TMP)
!!$       NC_AVG_DATA => ADD(NC_AVG_DATA,ICE_VEL_FILE_OBJECT() )
       !-----------------------------------------------------------------
    END IF
# endif


    !===================================================================
    ! MAKE A COPY OF THE DATA POINTER FILE TO IN WHICH TO DO THE SUM
    !===================================================================
    NC_AVG_SUM => COPY_FILE(NC_AVG_DATA)

    !===================================================================
    ! ALLOCATE ALL ITS DATA POINTERS TO MAKE THE EXTRA STORAGE
    !===================================================================
    CALL ALLOCATE_ASSOCIATED_VARS(NC_AVG_SUM)


    CALL KILL_DIMENSIONS

    

    !===================================================================
    ! NOW MAKE THE ACTUAL FILES WHICH ARE WRITTEN TO DISK
    !===================================================================

    ! Get the suffix from the NC_DAT file
    I = len_trim(NC_AVG%FNAME)
    dat_name = NC_AVG%FNAME(I-11:I)
    
    ! SETUP THE SUBDOMAIN FILES AND MAPS
    CALL SETUP_SUBDOMAINS(NCAV_SUBDOMAIN_FILES,NC_AVG_GRIDS)

    NUMG = size(NC_AVG_GRIDS)
    
    ! FOR EACH SUBDOMAIN MAKE THE FILE
    DO I = 1, NUMG

       ! DEFINEN DIMENSIONS FOR THIS GRID
       CALL DEFINE_DIMENSIONS(NC_AVG_GRIDS(I))
       
       IF(NC_AVG_GRIDS(I)%NAME /="FVCOM") THEN
          ! MAKE NEW FILE OBJECT
          tmp = TRIM(NC_AVG_GRIDS(I)%NAME)//TRIM(DAT_NAME)
          NCF => NEW_FILE(TRIM(TMP))
          
          ! SET THE CURRENT NAME IN THE GRID OBJECT
          NC_AVG_GRIDS(I)%NAME=NCF%FNAME

          ! SET THE FTIME OBJECT
          NCF%FTIME => NEW_FTIME()
          NCF%FTIME = NC_AVG%FTIME
          
          ! MAKE A TEMPORARY POINTER TO ADD THE FILE TO THE LIST
          NCF_TMP => NCF
          FILEHEAD => ADD(FILEHEAD,NCF_TMP)

       ELSE
          ! THIS IS THE AVERAGE FVCOM GRID FILE
          NCF => NC_AVG
          NC_AVG_GRIDS(I)%NAME = NC_AVG%FNAME
          
          INCLUDE_MASTER = .TRUE.
       END IF



       ! ADD THE GRID STUFF THAT DOES NOT COME FROM THE AVERAGE DATA
       NCF_TMP => GRID_FILE_OBJECT(NC_AVG_GRIDS(I))
       NCF => ADD(NCF,NCF_TMP)
!!$       NCF => ADD(NCF,GRID_FILE_OBJECT(NC_AVG_GRIDS(I)) )
       
       NCF_TMP => TIME_FILE_OBJECT()
       NCF => ADD(NCF,NCF_TMP)
!!$       NCF => ADD(NCF,TIME_FILE_OBJECT() )
       
       ATT => FIND_ATT(NCF,"title",FOUND)
       ATT%CHR(1) = trim(case_title)//"; Average output file!"
       NULLIFY(ATT)
       
       IF(NCAV_FILE_DATE) THEN
          NCF_TMP => FILE_DATE_OBJECT()
          NCF => ADD(NCF,NCF_TMP)
!!$          NCF => ADD(NCF,FILE_DATE_OBJECT() )
      END IF
    
       IF(NCAV_GRID_METRICS) THEN
          NCF_TMP => GRID_METRICS_FILE_OBJECT(NC_AVG_GRIDS(I))
          NCF => ADD(NCF,NCF_TMP)
!!$          NCF => ADD(NCF,GRID_METRICS_FILE_OBJECT(NC_AVG_GRIDS(I)) )
       END IF


       ! MAKE A COPY OF THE FILE POINTER TO ADD TO THE ACTUAL FILE
       ! THIS COPY POINTS TO THE MEMORY WHICH CONTAINS THE SUM/AVERAGE
       NCF_TMP =>COPY_FILE(NC_AVG_SUM)

       ! ADJUST THE NODE AND NELE DIMENSION OF THE COPIED FILE TO
       ! MATCH THE NEW FILE

       DIM => FIND_DIM(NCF_TMP,'node',FOUND)
       IF(.not.FOUND) CALL FATAL_ERROR&
            &("LOGICAL ERROR IN SETUP_AVGFILE: CAN'T FIND DIM NODE")
       DIM%DIM = NC_AVG_GRIDS(I)%MGL

       DIM => FIND_DIM(NCF_TMP,'nele',FOUND)
       IF(.not.FOUND) CALL FATAL_ERROR&
            &("LOGICAL ERROR IN SETUP_AVGFILE: CAN'T FIND DIM NELE")
       DIM%DIM = NC_AVG_GRIDS(I)%NGL


       NCF => ADD(NCF,NCF_TMP)

       
       IF (STARTUP_TYPE /= "crashrestart") THEN

          ! IF NOT CRASH RESTART - OVER RIDE THE STKCNT
          NCF%FTIME%NEXT_STKCNT = 0
          CALL NC_WRITE_FILE(NCF)
          NCF%FTIME%NEXT_STKCNT = 1
          
       ELSE
          NCF%CONNECTED = .TRUE.
          NCF%WRITABLE = .TRUE.
       END IF
       
       CALL KILL_DIMENSIONS


    END DO


    IF(.not. INCLUDE_MASTER) THEN

       CALL KILL_FILE(NC_AVG)
       NC_AVG => FIND_FILE(FILEHEAD,NC_AVG_GRIDS(1)%NAME,FOUND)

       IF(.NOT. FOUND) CALL FATAL_ERROR&
            &("LOGICAL ERROR IN SETUP_AVGFILE")
    END IF


    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END SETUP_AVGFILE"
  END SUBROUTINE SETUP_AVGFILE
!============================================================
  SUBROUTINE ADD_AVERAGE
    IMPLICIT NONE

    TYPE(NCVARP),POINTER :: CURRENT_DATA,CURRENT_SUM
    TYPE(NCVAR),POINTER :: VAR_DATA, VAR_SUM

    CURRENT_DATA => NC_AVG_DATA%VARS%NEXT
    CURRENT_SUM => NC_AVG_SUM%VARS%NEXT

    DO
       IF(.NOT. ASSOCIATED(CURRENT_DATA)) THEN
          EXIT ! END OF VAR LIST
       END IF

       IF(.NOT. ASSOCIATED(CURRENT_SUM)) THEN
          CALL FATAL_ERROR("ADD_AVERAGE: SUM AND DATA VAR LISTS DO NOT HAVE THE SAME LENGTH?")
       END IF

       IF(.NOT. ASSOCIATED(CURRENT_SUM%VAR)) THEN
          CALL FATAL_ERROR("ALLOCATE_ASSOCIATED_VARS: FOUND NULL VAR POINTER IN THE LIST")
       END IF

       IF(.NOT. ASSOCIATED(CURRENT_DATA%VAR)) THEN
          CALL FATAL_ERROR("ALLOCATE_ASSOCIATED_VARS: FOUND NULL VAR POINTER IN THE LIST")
       END IF

       VAR_SUM => CURRENT_SUM%VAR
       VAR_DATA => CURRENT_DATA%VAR


       IF(Associated(VAR_DATA%SCL_INT)) VAR_SUM%SCL_INT=VAR_DATA%SCL_INT+VAR_SUM%SCL_INT
       IF(Associated(VAR_DATA%VEC_INT)) VAR_SUM%VEC_INT=VAR_DATA%VEC_INT+VAR_SUM%VEC_INT
       IF(Associated(VAR_DATA%ARR_INT)) VAR_SUM%ARR_INT=VAR_DATA%ARR_INT+VAR_SUM%ARR_INT
       IF(Associated(VAR_DATA%CUB_INT)) VAR_SUM%CUB_INT=VAR_DATA%CUB_INT+VAR_SUM%CUB_INT
       
       IF(Associated(VAR_DATA%SCL_FLT)) VAR_SUM%SCL_FLT=VAR_DATA%SCL_FLT+VAR_SUM%SCL_FLT
       IF(Associated(VAR_DATA%VEC_FLT)) VAR_SUM%VEC_FLT=VAR_DATA%VEC_FLT+VAR_SUM%VEC_FLT
       IF(Associated(VAR_DATA%ARR_FLT)) VAR_SUM%ARR_FLT=VAR_DATA%ARR_FLT+VAR_SUM%ARR_FLT
       IF(Associated(VAR_DATA%CUB_FLT)) VAR_SUM%CUB_FLT=VAR_DATA%CUB_FLT+VAR_SUM%CUB_FLT
       
       IF(Associated(VAR_DATA%SCL_DBL)) VAR_SUM%SCL_DBL=VAR_DATA%SCL_DBL+VAR_SUM%SCL_DBL
       IF(Associated(VAR_DATA%VEC_DBL)) VAR_SUM%VEC_DBL=VAR_DATA%VEC_DBL+VAR_SUM%VEC_DBL
       IF(Associated(VAR_DATA%ARR_DBL)) VAR_SUM%ARR_DBL=VAR_DATA%ARR_DBL+VAR_SUM%ARR_DBL
       IF(Associated(VAR_DATA%CUB_DBL)) VAR_SUM%CUB_DBL=VAR_DATA%CUB_DBL+VAR_SUM%CUB_DBL


       CURRENT_DATA => CURRENT_DATA%NEXT
       CURRENT_SUM => CURRENT_SUM%NEXT
    END DO


  END SUBROUTINE ADD_AVERAGE
!============================================================
  SUBROUTINE DIVIDE_AVERAGE
    IMPLICIT NONE      
    TYPE(NCVARP),POINTER :: CURRENT
    TYPE(NCVAR),POINTER :: VAR
    REAL(DP) :: avg_steps

    avg_steps = SECONDS(NC_AVG%FTIME%INTERVAL)/SECONDS(IMDTI)

    CURRENT => NC_AVG_SUM%VARS%NEXT

    DO
       IF(.NOT. ASSOCIATED(CURRENT)) THEN
          EXIT ! END OF VAR LIST
       END IF

       IF(.NOT. ASSOCIATED(CURRENT%VAR)) THEN
          CALL FATAL_ERROR("ALLOCATE_ASSOCIATED_VARS: FOUND NULL VAR POINTER IN THE LIST")
       END IF

       VAR => CURRENT%VAR


       IF(Associated(VAR%SCL_INT)) VAR%SCL_INT=VAR%SCL_INT/avg_steps
       IF(Associated(VAR%VEC_INT)) VAR%VEC_INT=VAR%VEC_INT/avg_steps
       IF(Associated(VAR%ARR_INT)) VAR%ARR_INT=VAR%ARR_INT/avg_steps
       IF(Associated(VAR%CUB_INT)) VAR%CUB_INT=VAR%CUB_INT/avg_steps
       
       IF(Associated(VAR%SCL_FLT)) VAR%SCL_FLT=VAR%SCL_FLT/avg_steps
       IF(Associated(VAR%VEC_FLT)) VAR%VEC_FLT=VAR%VEC_FLT/avg_steps
       IF(Associated(VAR%ARR_FLT)) VAR%ARR_FLT=VAR%ARR_FLT/avg_steps
       IF(Associated(VAR%CUB_FLT)) VAR%CUB_FLT=VAR%CUB_FLT/avg_steps
       
       IF(Associated(VAR%SCL_DBL)) VAR%SCL_DBL=VAR%SCL_DBL/avg_steps
       IF(Associated(VAR%VEC_DBL)) VAR%VEC_DBL=VAR%VEC_DBL/avg_steps
       IF(Associated(VAR%ARR_DBL)) VAR%ARR_DBL=VAR%ARR_DBL/avg_steps
       IF(Associated(VAR%CUB_DBL)) VAR%CUB_DBL=VAR%CUB_DBL/avg_steps


       CURRENT => CURRENT%NEXT
    END DO
 
  END SUBROUTINE DIVIDE_AVERAGE
!============================================================
  SUBROUTINE ZERO_AVERAGE
    IMPLICIT NONE      
    TYPE(NCVARP),POINTER :: CURRENT
    TYPE(NCVAR),POINTER :: VAR
    REAL(DP) :: avg_steps

    avg_steps = SECONDS(NC_AVG%FTIME%INTERVAL)/SECONDS(IMDTI)

    CURRENT => NC_AVG_SUM%VARS%NEXT

    DO
       IF(.NOT. ASSOCIATED(CURRENT)) THEN
          EXIT ! END OF VAR LIST
       END IF

       IF(.NOT. ASSOCIATED(CURRENT%VAR)) THEN
          CALL FATAL_ERROR("ALLOCATE_ASSOCIATED_VARS: FOUND NULL VAR POINTER IN THE LIST")
       END IF

       VAR => CURRENT%VAR


       IF(Associated(VAR%SCL_INT)) VAR%SCL_INT=0
       IF(Associated(VAR%VEC_INT)) VAR%VEC_INT=0
       IF(Associated(VAR%ARR_INT)) VAR%ARR_INT=0
       IF(Associated(VAR%CUB_INT)) VAR%CUB_INT=0
       
       IF(Associated(VAR%SCL_FLT)) VAR%SCL_FLT=0.0_SPA
       IF(Associated(VAR%VEC_FLT)) VAR%VEC_FLT=0.0_SPA
       IF(Associated(VAR%ARR_FLT)) VAR%ARR_FLT=0.0_SPA
       IF(Associated(VAR%CUB_FLT)) VAR%CUB_FLT=0.0_SPA
       
       IF(Associated(VAR%SCL_DBL)) VAR%SCL_DBL=0.0_DP
       IF(Associated(VAR%VEC_DBL)) VAR%VEC_DBL=0.0_DP
       IF(Associated(VAR%ARR_DBL)) VAR%ARR_DBL=0.0_DP
       IF(Associated(VAR%CUB_DBL)) VAR%CUB_DBL=0.0_DP


       CURRENT => CURRENT%NEXT
    END DO
 
  END SUBROUTINE ZERO_AVERAGE
!=============================================================  
  SUBROUTINE SETUP_RSTFILE
    IMPLICIT NONE
    
    TYPE(GRID), SAVE :: MYGRID
    TYPE(NCFILE), POINTER ::NC_RST2

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START SETUP_RSTFILE"

    IF(DBG_SET(DBG_LOG)) THEN
       
       write(IPT,*)"!--------------------------------------------------"
       write(IPT,*)"! SETTING UP RESTART FILE OUTPUT..."
    END IF

    CALL SET_FVCOM_GRID(MYGRID)

    CALL DEFINE_DIMENSIONS(MYGRID)

    NC_RST2 => GRID_FILE_OBJECT(MYGRID)
    NC_RST => ADD(NC_RST,NC_RST2)
!!$    NC_RST => ADD(NC_RST,GRID_FILE_OBJECT(MYGRID) )

    NC_RST2 => TIME_FILE_OBJECT()
    NC_RST => ADD(NC_RST,NC_RST2)
!!$    NC_RST => ADD(NC_RST,TIME_FILE_OBJECT() )

    NC_RST2 => ZETA_FILE_OBJECT()
    NC_RST => ADD(NC_RST,NC_RST2)
!!$    NC_RST => ADD(NC_RST,ZETA_FILE_OBJECT() )

    NC_RST2 => FILE_DATE_OBJECT()
    NC_RST => ADD(NC_RST,NC_RST2)
!!$    NC_RST => ADD(NC_RST,FILE_DATE_OBJECT() )

    NC_RST2 => VELOCITY_FILE_OBJECT()
    NC_RST => ADD(NC_RST,NC_RST2)
!!$    NC_RST => ADD(NC_RST,VELOCITY_FILE_OBJECT() )

    NC_RST2 => AVERAGE_VEL_FILE_OBJECT()
    NC_RST => ADD(NC_RST,NC_RST2)
!!$    NC_RST => ADD(NC_RST,AVERAGE_VEL_FILE_OBJECT() )

    NC_RST2 => VERTICAL_VEL_FILE_OBJECT()
    NC_RST => ADD(NC_RST,NC_RST2)
!!$    NC_RST => ADD(NC_RST,VERTICAL_VEL_FILE_OBJECT() )

    NC_RST2 => TURBULENCE_FILE_OBJECT()
    NC_RST => ADD(NC_RST,NC_RST2)
!!$    NC_RST => ADD(NC_RST,TURBULENCE_FILE_OBJECT() )
    
    NC_RST2 => SALT_TEMP_FILE_OBJECT()
    NC_RST => ADD(NC_RST,NC_RST2)
!!$    NC_RST => ADD(NC_RST,SALT_TEMP_FILE_OBJECT() )

    NC_RST2 => RESTART_EXTRAS_FILE_OBJECT()
    NC_RST => ADD(NC_RST,NC_RST2)
!!$    NC_RST => ADD(NC_RST,RESTART_EXTRAS_FILE_OBJECT() )
    
    NC_RST2 => DENSITY_FILE_OBJECT()
    NC_RST => ADD(NC_RST,NC_RST2)
!!$    NC_RST => ADD(NC_RST,DENSITY_FILE_OBJECT() )

    IF(WETTING_DRYING_ON) THEN
       NC_RST2 => WET_DRY_FILE_OBJECT()
       NC_RST => ADD(NC_RST, NC_RST2)
!!$       NC_RST => ADD(NC_RST, WET_DRY_FILE_OBJECT() )
    END IF

#   if defined (BioGen)
    IF(BIOLOGICAL_MODEL)THEN
    NC_RST2 => BIO_FILE_OBJECT()
    NC_RST => ADD(NC_RST,NC_RST2)
!!$    NC_RST => ADD(NC_RST,BIO_FILE_OBJECT() )
    END IF
#   endif      

#   if defined (WATER_QUALITY)
    NC_RST2 => WQM_FILE_OBJECT()
    NC_RST => ADD(NC_RST,NC_RST2)
!!$    NC_RST => ADD(NC_RST,WQM_FILE_OBJECT() )
#   endif      

#   if defined (SEDIMENT)
    IF(SEDIMENT_MODEL)THEN
      NC_RST2 => SED_FILE_OBJECT()
      NC_RST => ADD(NC_RST,NC_RST2)
!!$      NC_RST => ADD(NC_RST,SED_FILE_OBJECT() )
! J. Ge for fluid mud layer
#   if defined (FLUID_MUD)
      NC_RST2 => FLUID_MUD_FILE_OBJECT()
      NC_RST => ADD(NC_RST,NC_RST2)
!!$      NC_RST => ADD(NC_RST,SED_FILE_OBJECT() )
#   endif   
! J. Ge for fluid mud layer
    ENDIF
#   endif      

#   if defined (WAVE_CURRENT_INTERACTION)      
    NC_RST2 => WAVE_PARA_FILE_OBJECT()
    NC_RST => ADD(NC_RST,NC_RST2) 
!!$    NC_RST => ADD(NC_RST,WAVE_PARA_FILE_OBJECT() ) 
    NC_RST2 => WAVE_STRESS_FILE_OBJECT()
    NC_RST => ADD(NC_RST,NC_RST2) 
!!$    NC_RST => ADD(NC_RST,WAVE_STRESS_FILE_OBJECT() ) 
#   endif  

#   if defined (DYE_RELEASE)
    IF(DYE_ON)THEN
      NC_RST2 => DYE_FILE_OBJECT()
      NC_RST => ADD(NC_RST,NC_RST2)
    ENDIF
#   endif

#   if defined (NH)
    NC_RST2 => NH_RST_FILE_OBJECT()
    NC_RST => ADD(NC_RST,NC_RST2)
!!$    NC_RST => ADD(NC_RST,NH_RST_FILE_OBJECT() )
#   endif

# if defined (ICE)
      !-----------------------------------------------------------------
      ! state variables
    NC_RST2 => ICE_RESTART_STATE_FILE_OBJECT()
    NC_RST => ADD(NC_RST,NC_RST2)
!!$    NC_RST => ADD(NC_RST,ICE_RESTART_STATE_FILE_OBJECT() )
      !-----------------------------------------------------------------
      ! velocity
      !-----------------------------------------------------------------
    NC_RST2 => ICE_VEL_FILE_OBJECT()
    NC_RST => ADD(NC_RST,NC_RST2)
!!$    NC_RST => ADD(NC_RST,ICE_VEL_FILE_OBJECT() )
      !-----------------------------------------------------------------
      ! fresh water, salt, and heat flux
      !-----------------------------------------------------------------
    NC_RST2 => ICE_EXTRA_FILE_OBJECT()
    NC_RST => ADD(NC_RST,NC_RST2)
!!$    NC_RST => ADD(NC_RST,ICE_EXTRA_FILE_OBJECT() )
# endif

    IF (STARTUP_TYPE /= "crashrestart") THEN
       CALL NC_WRITE_FILE(NC_RST)
       NC_RST%FTIME%NEXT_STKCNT = 1
    ELSE
       NC_RST%CONNECTED = .TRUE.
    END IF

    CALL KILL_DIMENSIONS

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END SETUP_RSTFILE"
  END SUBROUTINE SETUP_RSTFILE
!=============================================================  
  SUBROUTINE DUMP_DATA(NCF)
    IMPLICIT NONE  
    TYPE(NCFILE), POINTER ::NCF
    TYPE(NCFTIME), POINTER :: FTM

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START DUMP_DATA"

    IF(DBG_SET(DBG_IO)) CALL PRINT_FILE(NCF)
    IF(DBG_SET(DBG_sbrio)) CALL PRINT_VAR_LIST(NCF)

    FTM => NCF%FTIME


    IF (FTM%MAX_STKCNT .NE. 0 .AND. &
         & FTM%NEXT_STKCNT > FTM%MAX_STKCNT) THEN

       FTM%NEXT_STKCNT=0
       CALL INCRIMENT_FNAME(NCF%FNAME)
       NCF%CONNECTED=.FALSE.
   

       ! WRITE NEW FILE'S CONSTANT DATA (GRID ETC)
       CALL NC_WRITE_FILE(NCF)

       ! INCRIMENT THE STACK COUNT
       FTM%NEXT_STKCNT = 1

    END IF


    ! IF UPDATE IODATA BECOMES SPECIFIC TO DIFFERENT DATA SETS IT
    ! WILL HAVE TO BE MOVED INSIDE OF THE PARTICULAR OUTPUT STATEMENTS
    CALL UPDATE_IODATA(NCF,IntTime)

    CALL NC_WRITE_FILE(NCF)

    ! ONCE THE FILE IS WRITEN INCRIMENT THE FILE OBJECT TIME
    FTM%PREV_IO = IntTime
    FTM%NEXT_IO = FTM%NEXT_IO + FTM%INTERVAL

    ! INCRIMENT THE STACK COUNT
    FTM%PREV_STKCNT = FTM%NEXT_STKCNT
    FTM%NEXT_STKCNT = FTM%NEXT_STKCNT + 1
    
    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END DUMP_DATA"
    
  END SUBROUTINE DUMP_DATA
!=============================================================  
  SUBROUTINE DEFINE_DIMENSIONS(G)
# if defined (ICE)
    USE mod_ICE
    USE mod_ICE2D
# endif

# if defined (SEDIMENT)

#  if defined (ORIG_SED)
  USE mod_sed
#  elif defined (CSTMS_SED)
  USE mod_sed_cstms
#  endif

! J. Ge for fluid mud layer
# if defined (FLUID_MUD)
  USE Mod_Fluid_Mud
# endif
! J. Ge for fluid mud layer

# endif

!# if defined (WAVE_CURRENT_INTERACTION)
!    USE SWCOMM3
!    USE VARS_WAVE
!    USE MOD_WAVE_CURRENT_INTERACTION
!# endif
    USE BCS
    IMPLICIT NONE
    TYPE(GRID) :: G
    INTEGER :: SBUF, SOURCE
    INTEGER :: RBUF, DEST
    

    DIM_nele       => NC_MAKE_DIM(name='nele',len=G%ngl)
    DIM_node       => NC_MAKE_DIM(name='node',len=G%mgl)

    DIM_three      => NC_MAKE_DIM(name='three',len=3)
    DIM_four       => NC_MAKE_DIM(name='four',len=4)

    DIM_siglev     => NC_MAKE_DIM(name='siglev',len=G%kb)
    DIM_siglay     => NC_MAKE_DIM(name='siglay',len=G%kbm1)

# if defined (WAVE_CURRENT_INTERACTION)
    DIM_MSC        => NC_MAKE_DIM(name='msc',len=G%msc)
    DIM_MDC        => NC_MAKE_DIM(name='mdc',len=G%mdc)
# endif
# if defined (SEDIMENT)
    IF(SEDIMENT_MODEL)THEN
      DIM_bedlay    => NC_MAKE_RUNTIME_DIM(name='bedlay',len=nbed)
    ENDIF
# endif
    
    DIM_DateStrLen => NC_MAKE_DIM(name='DateStrLen',len=DateStrLen)
    DIM_time       => NC_MAKE_DIM(name='time',len=NF90_UNLIMITED)

# if defined (ICE)
    DIM_GRID       => NC_MAKE_DIM(name='dimension2',len=JMT_LOCAL)
    DIM_NCAT       => NC_MAKE_DIM(name='category',len=NCAT)
    DIM_Ntilay     => NC_MAKE_DIM(name='totallayer',len=Ntilay)
# else
    NULLIFY(DIM_GRID,DIM_NCAT,DIM_Ntilay)
# endif   
   
    ! ONLY USED IF OBC IS ON
    IF(OBC_ON) THEN
       DIM_nobc       => NC_MAKE_DIM(name='nobc',len=IOBCN_GL)
    ELSE
       nullify(DIM_NOBC)
    END IF
    ! ONLY USED IF LOND SHORE FLOW BOUNDARY IS ON
    IF(OBC_LONGSHORE_FLOW_ON) THEN
       DIM_nlsf       => NC_MAKE_DIM(name='nlsf',len=nobclsf_GL)
    ELSE
       nullify(DIM_NLSF)
    END IF

    DIM_MaxNode    => NC_MAKE_RUNTIME_DIM(name='maxnode',len=MX_NBR_ELEM+3)

    DIM_MaxElem    => NC_MAKE_RUNTIME_DIM(name='maxelem',len=MX_NBR_ELEM+1)


  END SUBROUTINE DEFINE_DIMENSIONS

!=============================================================  
  SUBROUTINE DEFINE_DIMENSIONS_SURFACE(G)
# if defined (ICE)
    USE mod_ICE
    USE mod_ICE2D
# endif

    USE BCS
    IMPLICIT NONE
    TYPE(GRID) :: G
    INTEGER :: SBUF, SOURCE
    INTEGER :: RBUF, DEST
    

    DIM_nele       => NC_MAKE_DIM(name='nele',len=G%ngl)
    DIM_node       => NC_MAKE_DIM(name='node',len=G%mgl)

    DIM_three      => NC_MAKE_DIM(name='three',len=3)
    DIM_four       => NC_MAKE_DIM(name='four',len=4)

# if defined (WAVE_CURRENT_INTERACTION)
    DIM_MSC        => NC_MAKE_DIM(name='msc',len=G%msc)
    DIM_MDC        => NC_MAKE_DIM(name='mdc',len=G%mdc)
# endif
    
    DIM_DateStrLen => NC_MAKE_DIM(name='DateStrLen',len=DateStrLen)
    DIM_time       => NC_MAKE_DIM(name='time',len=NF90_UNLIMITED)

# if defined (ICE)
    DIM_GRID       => NC_MAKE_DIM(name='dimension2',len=JMT_LOCAL)
    DIM_NCAT       => NC_MAKE_DIM(name='category',len=NCAT)
    DIM_Ntilay     => NC_MAKE_DIM(name='totallayer',len=Ntilay)
# else
    NULLIFY(DIM_GRID,DIM_NCAT,DIM_Ntilay)
# endif   
   
    DIM_MaxNode    => NC_MAKE_RUNTIME_DIM(name='maxnode',len=MX_NBR_ELEM+3)

    DIM_MaxElem    => NC_MAKE_RUNTIME_DIM(name='maxelem',len=MX_NBR_ELEM+1)


  END SUBROUTINE DEFINE_DIMENSIONS_SURFACE

!=============================================================  
  SUBROUTINE KILL_DIMENSIONS
    IMPLICIT NONE

    IF (ASSOCIATED(DIM_NELE)) CALL KILL_DIM(DIM_NELE)
    IF (ASSOCIATED(DIM_NODE)) CALL KILL_DIM(DIM_NODE)

    IF (ASSOCIATED(DIM_THREE)) CALL KILL_DIM(DIM_THREE)
    IF (ASSOCIATED(DIM_FOUR)) CALL KILL_DIM(DIM_FOUR)

    IF (ASSOCIATED(DIM_SIGLAY)) CALL KILL_DIM(DIM_SIGLAY)
    IF (ASSOCIATED(DIM_SIGLEV)) CALL KILL_DIM(DIM_SIGLEV)

    IF (ASSOCIATED(DIM_DATESTRLEN)) CALL KILL_DIM(DIM_DATESTRLEN)
    IF (ASSOCIATED(DIM_TIME)) CALL KILL_DIM(DIM_TIME)

    IF (ASSOCIATED(DIM_GRID)) CALL KILL_DIM(DIM_GRID)
    IF (ASSOCIATED(DIM_NCAT)) CALL KILL_DIM(DIM_NCAT)
    IF (ASSOCIATED(DIM_Ntilay)) CALL KILL_DIM(DIM_Ntilay)

    IF (ASSOCIATED(DIM_nobc)) CALL KILL_DIM(DIM_nobc)
    IF (ASSOCIATED(DIM_nlsf)) CALL KILL_DIM(DIM_nlsf)

    IF (ASSOCIATED(DIM_MaxNode)) CALL KILL_DIM(DIM_MaxNode)
    IF (ASSOCIATED(DIM_MaxElem)) CALL KILL_DIM(DIM_MaxElem)    

  END SUBROUTINE KILL_DIMENSIONS
!=============================================================  
  SUBROUTINE KILL_DIMENSIONS_SURFACE
    IMPLICIT NONE

    IF (ASSOCIATED(DIM_NELE)) CALL KILL_DIM(DIM_NELE)
    IF (ASSOCIATED(DIM_NODE)) CALL KILL_DIM(DIM_NODE)

    IF (ASSOCIATED(DIM_THREE)) CALL KILL_DIM(DIM_THREE)
    IF (ASSOCIATED(DIM_FOUR)) CALL KILL_DIM(DIM_FOUR)

    IF (ASSOCIATED(DIM_DATESTRLEN)) CALL KILL_DIM(DIM_DATESTRLEN)
    IF (ASSOCIATED(DIM_TIME)) CALL KILL_DIM(DIM_TIME)

    IF (ASSOCIATED(DIM_GRID)) CALL KILL_DIM(DIM_GRID)
    IF (ASSOCIATED(DIM_NCAT)) CALL KILL_DIM(DIM_NCAT)
    IF (ASSOCIATED(DIM_Ntilay)) CALL KILL_DIM(DIM_Ntilay)

    IF (ASSOCIATED(DIM_MaxNode)) CALL KILL_DIM(DIM_MaxNode)
    IF (ASSOCIATED(DIM_MaxElem)) CALL KILL_DIM(DIM_MaxElem)    

  END SUBROUTINE KILL_DIMENSIONS_SURFACE
!=============================================================  
  FUNCTION GRID_FILE_OBJECT(G) RESULT(NCF)
    USE MOD_CLOCK
    USE MOD_FORCE
    IMPLICIT NONE
    
    TYPE(GRID) :: G
    INTEGER, POINTER :: partition(:)
    INTEGER :: status, I
    LOGICAL, SAVE :: IOPROC_ALLOCATED = .FALSE.
    TYPE(NCFILE), POINTER :: NCF
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCATT),  POINTER :: ATT

    character(len=100)    :: timestamp, temp, netcdf_convention

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START GRID_FILE_OBJECT"
    
    ! IO PROC MUST ALLOCATE SPACE FOR THE ARRAYS 
    ! THESE ARRAYS MUST HAVE THE ATTRIBUTE SAVE AND FOR CLARITY
    ! SHOULD HAVE THE SAME NAME AS THOSE USED ON THE OTHER PROCESSORS!
    IF(IOPROC .AND. .NOT. IOPROC_ALLOCATED) THEN


       ! MGL AND NGL ARE THE BIGEST THING GOING - USE THEM TO MAKE
       ! SURE WE DON'T GET A SIG SEV - MULTIPLE FILE GRIDS MAY USE
       ! THESE POINTERS!
       IOPROC_ALLOCATED = .true.

       allocate(partition(NGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:partition")
       partition = 0

       allocate(xm(MGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:XM")
       xm = 0.0_SP

       allocate(ym(MGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:YM")
       ym = 0.0_SP

       allocate(LON(MGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:LON")
       lon = 0.0_SP

       allocate(LAT(MGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:LAT")
       lat = 0.0_SP

       allocate(xmc(NGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:XMC")
       xmc = 0.0_SP

       allocate(ymc(NGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:YMC")
       ymc = 0.0_SP

       allocate(LONC(NGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:LONC")
       lonc = 0.0_SP

       allocate(LATC(NGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:LATC")
       latc = 0.0_SP

       allocate(zz(MGL,KB),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:ZZ")
       zz = 0.0_SP

       allocate(z(MGL,KB),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:Z")
       z = 0.0_SP

       allocate(h(MGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:H")
       h = 0.0_SP
       
!# if defined(UCF)
       allocate(zz1(NGL,KB),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:ZZ1")
       zz1 = 0.0_SP

       allocate(z1(NGL,KB),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:Z1")
       z1 = 0.0_SP

       allocate(h1(NGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORcY ON IO PROC FOR OUTPUT DATA:H1")
       h1 = 0.0_SP
!# endif

       
    END IF

    IF(IOPROC) THEN

       ! ALLOCATE NV BASED ON THE GRID DIMENSION
       IF(ASSOCIATED(G%NV)) THEN
          
          IF(UBOUND(G%nv,1)/=DIM_NELE%DIM) THEN
             CALL FATAL_ERROR &
                  &("GRID DATA NV HAS ALREADY BEEN ASSOICATED ON THE IOPROC",&
                  & "AND THE DIMENSION DOES NOT MATCH THE CURRENT FILE!",&
                  &"GRID NAME:"//TRIM(G%NAME))
          END IF

       ELSE
          allocate(G%nv(G%NGL,3),stat=status)
          IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:G%NV")
          G%nv = 0
       END IF
    END IF



    ! ALLOCATE THE NEW FILE OBJECT
    NCF => NEW_FILE()


    ! ADD THE FILE ATTRIBUTES
    ATT => NC_MAKE_ATT(name='title',values=trim(case_title)) 
    NCF => ADD(NCF,ATT)


    ATT => NC_MAKE_ATT(name='institution',values=trim(institution)) 
    NCF => ADD(NCF,ATT)

    ATT => NC_MAKE_ATT(name='source',values=trim(fvcom_version)) 
    NCF => ADD(NCF,ATT)

    call get_timestamp(temp)
    timestamp = 'model started at: '//trim(temp)

    ATT => NC_MAKE_ATT(name='history',values=trim(timestamp)) 
    NCF => ADD(NCF,ATT)

    ATT => NC_MAKE_ATT(name='references',values=trim(fvcom_website)) 
    NCF => ADD(NCF,ATT)

    netcdf_convention = 'CF-1.0'
    ATT => NC_MAKE_ATT(name='Conventions',values=trim(netcdf_convention)) 
    NCF => ADD(NCF,ATT)

#  if defined (SPHERICAL)
    ATT => NC_MAKE_ATT(name='CoordinateSystem',values="GeoReferenced" ) 
    NCF => ADD(NCF,ATT)
# else
    ATT => NC_MAKE_ATT(name='CoordinateSystem',values="Cartesian" ) 
    NCF => ADD(NCF,ATT)
#endif

    ATT => NC_MAKE_ATT(name='CoordinateProjection',values=PROJECTION_REFERENCE ) 
    IF(ASSOCIATED(ATT)) NCF => ADD(NCF,ATT)

    IF(TRIM(PRG_NAME) == "FVCOM") THEN
       
       ATT=> NC_Make_Runtime_Att_CHR(name='Tidal_Forcing',values=TIDE_FORCING_COMMENTS)
       IF(ASSOCIATED(ATT)) NCF => ADD(NCF,ATT)
       
       ATT=> NC_Make_Runtime_Att_CHR(name='River_Forcing',values=RIVER_FORCING_COMMENTS)
       IF(ASSOCIATED(ATT)) NCF => ADD(NCF,ATT)
       
       ATT=> NC_Make_Runtime_Att_CHR(name='GroundWater_Forcing',values=GWATER_FORCING_COMMENTS)
       IF(ASSOCIATED(ATT)) NCF => ADD(NCF,ATT)
       
       ATT=> NC_Make_Runtime_Att_CHR(name='Surface_Heat_Forcing',values=HEAT_FORCING_COMMENTS)
       IF(ASSOCIATED(ATT)) NCF => ADD(NCF,ATT)
       
       ATT=> NC_Make_Runtime_Att_CHR(name='Surface_Wind_Forcing',values=WINDS_FORCING_COMMENTS)
       IF(ASSOCIATED(ATT)) NCF => ADD(NCF,ATT)

!------Jadon for offline wave forcing
# if defined (WAVE_CURRENT_INTERACTION)
       ATT=> NC_Make_Runtime_Att_CHR(name='Surface_Wave_Forcing',values=WAVES_FORCING_COMMENTS)
       IF(ASSOCIATED(ATT)) NCF => ADD(NCF,ATT)
# endif

       ATT=> NC_Make_Runtime_Att_CHR(name='Surface_PrecipEvap_Forcing',values=PRECIP_FORCING_COMMENTS)
       IF(ASSOCIATED(ATT)) NCF => ADD(NCF,ATT)
       
       IF (ICING_MODEL) THEN
          ATT=> NC_Make_Runtime_Att_CHR(name='Icing_Model_Forcing',values=ICING_FORCING_COMMENTS)
          IF(ASSOCIATED(ATT)) NCF => ADD(NCF,ATT)
       END IF
       
       IF (ICE_MODEL) THEN
          ATT=> NC_Make_Runtime_Att_CHR(name='Ice_Model_Forcing',values=ICE_FORCING_COMMENTS)
          IF(ASSOCIATED(ATT)) NCF => ADD(NCF,ATT)
       END IF
       
       
       IF(OBC_LONGSHORE_FLOW_ON) THEN
          aTT=> NC_MAKE_ATT(name='Special_Physical_processes',&
               & values='long shore flow adjustment for thermal wind and win&
               &d driven setup')
          NCF => ADD(NCF,ATT)
       END IF

    END IF

    ! ADD THE VARIABLES

    ! NPROCS
    VAR  => NC_MAKE_AVAR(name='nprocs',values=nprocs)

    ATT  => NC_MAKE_ATT(name='long_name',values='number of processors') 
    VAR  => ADD(VAR,ATT)
    NCF  => ADD(NCF,VAR)

#  if defined (MULTIPROCESSOR)
    ! PARTITION
    ALLOCATE(partition(G%NT)); partition=myid
    VAR  => NC_MAKE_PVAR(name='partition',values=partition,DIM1=DIM_nele)
    NULLIFY(PARTITION)
    ! THIS CAN CREATE A MEMORY LEAK - BE CAREFUL!
    

    ATT  => NC_MAKE_ATT(name='long_name',values='partition') 
    VAR  => ADD(VAR,ATT)

# if defined(UCF)
    ATT  => NC_MAKE_ATT(name='grid',values='grid1') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid_location',values='center') 
    VAR  => ADD(VAR,ATT)
# endif

    NCF  => ADD(NCF,VAR)
# endif



    ! X
    VAR  => NC_MAKE_AVAR(name='x',values=xm,DIM1=DIM_node)
   
    ATT  => NC_MAKE_ATT(name='long_name',values='nodal x-coordinate') 
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='units',values='meters') 
    VAR  => ADD(VAR,ATT)
    NCF  => ADD(NCF,VAR)
    
    ! Y
    VAR  => NC_MAKE_AVAR(name='y',values=ym,DIM1=DIM_node)

    ATT  => NC_MAKE_ATT(name='long_name',values='nodal y-coordinate') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='meters') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)


    ! LON
    VAR  => NC_MAKE_AVAR(name='lon',values=lon,DIM1=DIM_node)

    ATT  => NC_MAKE_ATT(name='long_name',values='nodal longitude') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='standard_name',values='longitude') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='degrees_east') 
    VAR  => ADD(VAR,ATT)
    NCF  => ADD(NCF,VAR)

    ! LAT
    VAR  => NC_MAKE_AVAR(name='lat',values=LAT,DIM1=DIM_node)

    ATT  => NC_MAKE_ATT(name='long_name',values='nodal latitude') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='standard_name',values='latitude') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='degrees_north') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    IF (ALLOCATED(XMC)) THEN
       ! XMC
       VAR  => NC_MAKE_AVAR(name='xc',values=xmc,DIM1=DIM_nele)
       
       ATT  => NC_MAKE_ATT(name='long_name',values='zonal x-coordinate') 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='units',values='meters') 
       VAR  => ADD(VAR,ATT)
       NCF  => ADD(NCF,VAR)
    END IF

    IF (ALLOCATED(YMC)) THEN
       ! YMC
       VAR  => NC_MAKE_AVAR(name='yc',values=ymc,DIM1=DIM_nele)
       
       ATT  => NC_MAKE_ATT(name='long_name',values='zonal y-coordinate') 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='units',values='meters') 
       VAR  => ADD(VAR,ATT)
       
       NCF  => ADD(NCF,VAR)
    END IF

    IF (ALLOCATED(LONC)) THEN
       ! LONC
       VAR  => NC_MAKE_AVAR(name='lonc',values=lonc,DIM1=DIM_nele)
       
       ATT  => NC_MAKE_ATT(name='long_name',values='zonal longitude') 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='standard_name',values='longitude') 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='units',values='degrees_east') 
       VAR  => ADD(VAR,ATT)
       NCF  => ADD(NCF,VAR)
    END IF

    IF (ALLOCATED(LATC)) THEN
       ! LATC
       VAR  => NC_MAKE_AVAR(name='latc',values=LATC,DIM1=DIM_nele)
       
       ATT  => NC_MAKE_ATT(name='long_name',values='zonal latitude') 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='standard_name',values='latitude') 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='units',values='degrees_north') 
       VAR  => ADD(VAR,ATT)
       
       NCF  => ADD(NCF,VAR)
    END IF

    IF (ALLOCATED(zz)) THEN
       ! siglay
       VAR  => NC_MAKE_AVAR(name='siglay',&
            & values=zz,&
            & DIM1= DIM_node,&
            & DIM2= DIM_siglay )
       
       ATT  => NC_MAKE_ATT(name='long_name',values='Sigma Layers') 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='standard_name',values='ocean_sigma/general_coordinate') 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='positive',values='up') 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='valid_min',values=-1.0_SPA) 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='valid_max',values=0.0_SPA) 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='formula_terms',values='sigma: siglay eta: zeta depth: h') 
       VAR  => ADD(VAR,ATT)
       
       NCF  => ADD(NCF,VAR)
    END IF

    IF (ALLOCATED(Z)) THEN
       ! siglev
       VAR  => NC_MAKE_AVAR(name='siglev',&
            & values=z, DIM1= DIM_node, DIM2= DIM_siglev )
       
       ATT  => NC_MAKE_ATT(name='long_name',values='Sigma Levels') 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='standard_name',values='ocean_sigma/general_coordinate') 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='positive',values='up') 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='valid_min',values=-1.0_SPA) 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='valid_max',values=0.0_SPA) 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='formula_terms',values='sigma:siglay eta: zeta depth: h') 
       VAR  => ADD(VAR,ATT)
       
       NCF  => ADD(NCF,VAR)
    END IF

!# if defined(UCF)
    ! ADD SIGMA VALUE AT CELL CENTERS!
    IF (ALLOCATED(zz1)) THEN
       ! siglay
       VAR  => NC_MAKE_AVAR(name='siglay_center',&
            & values=zz1,&
            & DIM1= DIM_nele,&
            & DIM2= DIM_siglay )
       
       ATT  => NC_MAKE_ATT(name='long_name',values='Sigma Layers') 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='standard_name',values='ocean_sigma/general_coordinate') 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='positive',values='up') 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='valid_min',values=-1.0_SPA) 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='valid_max',values=0.0_SPA) 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='formula_terms',values='sigma: siglay_center eta: zeta_center depth: h_center') 
       VAR  => ADD(VAR,ATT)
       
       NCF  => ADD(NCF,VAR)
    END IF

    IF (ALLOCATED(Z1)) THEN
       ! siglev
       VAR  => NC_MAKE_AVAR(name='siglev_center',&
            & values=z1, DIM1= DIM_nele, DIM2= DIM_siglev )
       
       ATT  => NC_MAKE_ATT(name='long_name',values='Sigma Levels') 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='standard_name',values='ocean_sigma/general_coordinate') 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='positive',values='up') 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='valid_min',values=-1.0_SPA) 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='valid_max',values=0.0_SPA) 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='formula_terms',values='sigma:siglay_center eta: zeta_center depth: h_center') 
       VAR  => ADD(VAR,ATT)
       
       NCF  => ADD(NCF,VAR)
    END IF

    IF (ALLOCATED(h1)) THEN
       ! h1  (if morphodynamics is on, dump time-dependent bathymetry)
# if defined(SEDIMENT)
       IF(MORPHO_MODEL)THEN
         VAR  => NC_MAKE_AVAR(name='h_center',values=h1, DIM1= DIM_nele, DIM2= DIM_time)
       ELSE
         VAR  => NC_MAKE_AVAR(name='h_center',values=h1, DIM1= DIM_nele)
       ENDIF
# else
       VAR  => NC_MAKE_AVAR(name='h_center',values=h1, DIM1= DIM_nele)
# endif

       ATT  => NC_MAKE_ATT(name='long_name',values='Bathymetry') 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='standard_name',values='sea_floor_depth_below_geoid') 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='units',values='m')
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='positive',values='down') 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='grid',values='grid1 grid3') 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='coordinates',values='latc lonc') 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='grid_location',values='center') 
       VAR  => ADD(VAR,ATT)

       NCF  => ADD(NCF,VAR)
    END IF


!# endif

    IF (ALLOCATED(h)) THEN
       ! h (if morphodynamics is on, dump time-dependent bathymetry)
# if defined(SEDIMENT)
       IF(MORPHO_MODEL)THEN
         VAR  => NC_MAKE_AVAR(name='h_morpho',values=h, DIM1= DIM_node, DIM2= DIM_time)
         ATT  => NC_MAKE_ATT(name='long_name',values='Bathymetry_Morpho') 
         VAR  => ADD(VAR,ATT)
       
         ATT  => NC_MAKE_ATT(name='standard_name',values='sea_floor_depth_below_geoid') 
         VAR  => ADD(VAR,ATT)

         ATT  => NC_MAKE_ATT(name='units',values='m')
         VAR  => ADD(VAR,ATT)

         ATT  => NC_MAKE_ATT(name='positive',values='down') 
         VAR  => ADD(VAR,ATT)

         NCF  => ADD(NCF,VAR)
       END IF
       VAR  => NC_MAKE_AVAR(name='h',values=h, DIM1= DIM_node)

# else
       VAR  => NC_MAKE_AVAR(name='h',values=h, DIM1= DIM_node)
# endif

       ATT  => NC_MAKE_ATT(name='long_name',values='Bathymetry') 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='standard_name',values='sea_floor_depth_below_geoid') 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='units',values='m')
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='positive',values='down') 
       VAR  => ADD(VAR,ATT)

# if defined(UCF)
       ATT  => NC_MAKE_ATT(name='grid',values='grid1 grid3') 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='coordinates',values="lat lon") 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='grid_location',values='node') 
       VAR  => ADD(VAR,ATT)
# else
       ATT  => NC_MAKE_ATT(name='grid',values='Bathymetry_Mesh') 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar) 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='type',values='data') 
       VAR  => ADD(VAR,ATT)
# endif

       NCF  => ADD(NCF,VAR)
    END IF

# if defined(UCF)
    ! THE UCF STANDARD INCLUDES ONE CONNECTIVITY ARRAY FOR EACH MESH
    ! nv
    VAR  => NC_MAKE_PVAR(name='grid1',&
         & values=G%nv, DIM1= DIM_nele, DIM2= DIM_three)

    ATT  => NC_MAKE_ATT(name='cell_type',values='tri_ccw') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='index_start',values='0') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='standard_name',values='connectivity_array') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates_node',values='lon lat') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! nv
    VAR  => NC_MAKE_PVAR(name='grid2',&
         & values=G%nv, DIM1= DIM_nele, DIM2= DIM_three)

    ATT  => NC_MAKE_ATT(name='cell_type',values='tri_ccw') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='index_start',values='0') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='standard_name',values='connectivity_array') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates_node',values='lon lat zeta') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! nv
    VAR  => NC_MAKE_PVAR(name='grid3',&
         & values=G%nv, DIM1= DIM_nele, DIM2= DIM_three)

    ATT  => NC_MAKE_ATT(name='cell_type',values='tri_ccw') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='index_start',values='0') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='standard_name',values='connectivity_array') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates_node',values='lon lat h') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! nv
    VAR  => NC_MAKE_PVAR(name='grid4',&
         & values=G%nv, DIM1= DIM_nele, DIM2= DIM_three)

    ATT  => NC_MAKE_ATT(name='cell_type',values='tri_ccw') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='index_start',values='0') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='standard_name',values='connectivity_array') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates_node',values='lon lat siglay') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates_cell',values='lonc latc siglay_center') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! nv
    VAR  => NC_MAKE_PVAR(name='grid5',&
         & values=G%nv, DIM1= DIM_nele, DIM2= DIM_three)

    ATT  => NC_MAKE_ATT(name='cell_type',values='tri_ccw') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='index_start',values='0') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='standard_name',values='connectivity_array') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates_node',values='lon lat siglev') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates_cell',values='lonc latc siglev_center') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

# else
    ! nv
    VAR  => NC_MAKE_PVAR(name='nv',&
         & values=G%nv, DIM1= DIM_nele, DIM2= DIM_three)

    ATT  => NC_MAKE_ATT(name='long_name',values='nodes surrounding element') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)
# endif

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END GRID_FILE_OBJECT"


  END FUNCTION GRID_FILE_OBJECT

!=============================================================  
  FUNCTION GRID_FILE_OBJECT_SURFACE(G) RESULT(NCF)
    USE MOD_CLOCK
    USE MOD_FORCE
    IMPLICIT NONE
    
    TYPE(GRID) :: G
    INTEGER, POINTER :: partition(:)
    INTEGER :: status, I
    LOGICAL, SAVE :: IOPROC_ALLOCATED = .FALSE.
    TYPE(NCFILE), POINTER :: NCF
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCATT),  POINTER :: ATT

    character(len=100)    :: timestamp, temp, netcdf_convention

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START GRID_FILE_OBJECT_SURFACE"
    
    ! IO PROC MUST ALLOCATE SPACE FOR THE ARRAYS 
    ! THESE ARRAYS MUST HAVE THE ATTRIBUTE SAVE AND FOR CLARITY
    ! SHOULD HAVE THE SAME NAME AS THOSE USED ON THE OTHER PROCESSORS!
    IF(IOPROC .AND. .NOT. IOPROC_ALLOCATED) THEN


       ! MGL AND NGL ARE THE BIGEST THING GOING - USE THEM TO MAKE
       ! SURE WE DON'T GET A SIG SEV - MULTIPLE FILE GRIDS MAY USE
       ! THESE POINTERS!
       IOPROC_ALLOCATED = .true.

       allocate(partition(NGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:partition")
       partition = 0

       allocate(xm(MGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:XM")
       xm = 0.0_SP

       allocate(ym(MGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:YM")
       ym = 0.0_SP

       allocate(LON(MGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:LON")
       lon = 0.0_SP

       allocate(LAT(MGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:LAT")
       lat = 0.0_SP

       allocate(xmc(NGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:XMC")
       xmc = 0.0_SP

       allocate(ymc(NGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:YMC")
       ymc = 0.0_SP

       allocate(LONC(NGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:LONC")
       lonc = 0.0_SP

       allocate(LATC(NGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:LATC")
       latc = 0.0_SP

       allocate(h(MGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:H")
       h = 0.0_SP
       
!# if defined(UCF)
       allocate(h1(NGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORcY ON IO PROC FOR OUTPUT DATA:H1")
       h1 = 0.0_SP
!# endif

       
    END IF

    IF(IOPROC) THEN

       ! ALLOCATE NV BASED ON THE GRID DIMENSION
       IF(ASSOCIATED(G%NV)) THEN
          
          IF(UBOUND(G%nv,1)/=DIM_NELE%DIM) THEN
             CALL FATAL_ERROR &
                  &("GRID DATA NV HAS ALREADY BEEN ASSOICATED ON THE IOPROC",&
                  & "AND THE DIMENSION DOES NOT MATCH THE CURRENT FILE!",&
                  &"GRID NAME:"//TRIM(G%NAME))
          END IF

       ELSE
          allocate(G%nv(G%NGL,3),stat=status)
          IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:G%NV")
          G%nv = 0
       END IF
    END IF



    ! ALLOCATE THE NEW FILE OBJECT
    NCF => NEW_FILE()


    ! ADD THE FILE ATTRIBUTES
    ATT => NC_MAKE_ATT(name='title',values=trim(case_title)) 
    NCF => ADD(NCF,ATT)


    ATT => NC_MAKE_ATT(name='institution',values=trim(institution)) 
    NCF => ADD(NCF,ATT)

    ATT => NC_MAKE_ATT(name='source',values=trim(fvcom_version)) 
    NCF => ADD(NCF,ATT)

    call get_timestamp(temp)
    timestamp = 'model started at: '//trim(temp)

    ATT => NC_MAKE_ATT(name='history',values=trim(timestamp)) 
    NCF => ADD(NCF,ATT)

    ATT => NC_MAKE_ATT(name='references',values=trim(fvcom_website)) 
    NCF => ADD(NCF,ATT)

    netcdf_convention = 'CF-1.0'
    ATT => NC_MAKE_ATT(name='Conventions',values=trim(netcdf_convention)) 
    NCF => ADD(NCF,ATT)

#  if defined (SPHERICAL)
    ATT => NC_MAKE_ATT(name='CoordinateSystem',values="GeoReferenced" ) 
    NCF => ADD(NCF,ATT)
# else
    ATT => NC_MAKE_ATT(name='CoordinateSystem',values="Cartesian" ) 
    NCF => ADD(NCF,ATT)
#endif

    ATT => NC_MAKE_ATT(name='CoordinateProjection',values=PROJECTION_REFERENCE ) 
    IF(ASSOCIATED(ATT)) NCF => ADD(NCF,ATT)

    IF(TRIM(PRG_NAME) == "FVCOM") THEN
       
       ATT=> NC_Make_Runtime_Att_CHR(name='Tidal_Forcing',values=TIDE_FORCING_COMMENTS)
       IF(ASSOCIATED(ATT)) NCF => ADD(NCF,ATT)
       
       ATT=> NC_Make_Runtime_Att_CHR(name='River_Forcing',values=RIVER_FORCING_COMMENTS)
       IF(ASSOCIATED(ATT)) NCF => ADD(NCF,ATT)
       
       ATT=> NC_Make_Runtime_Att_CHR(name='GroundWater_Forcing',values=GWATER_FORCING_COMMENTS)
       IF(ASSOCIATED(ATT)) NCF => ADD(NCF,ATT)
       
       ATT=> NC_Make_Runtime_Att_CHR(name='Surface_Heat_Forcing',values=HEAT_FORCING_COMMENTS)
       IF(ASSOCIATED(ATT)) NCF => ADD(NCF,ATT)
       
       ATT=> NC_Make_Runtime_Att_CHR(name='Surface_Wind_Forcing',values=WINDS_FORCING_COMMENTS)
       IF(ASSOCIATED(ATT)) NCF => ADD(NCF,ATT)

       ATT=> NC_Make_Runtime_Att_CHR(name='Surface_PrecipEvap_Forcing',values=PRECIP_FORCING_COMMENTS)
       IF(ASSOCIATED(ATT)) NCF => ADD(NCF,ATT)
       
       IF (ICING_MODEL) THEN
          ATT=> NC_Make_Runtime_Att_CHR(name='Icing_Model_Forcing',values=ICING_FORCING_COMMENTS)
          IF(ASSOCIATED(ATT)) NCF => ADD(NCF,ATT)
       END IF
       
       IF (ICE_MODEL) THEN
          ATT=> NC_Make_Runtime_Att_CHR(name='Ice_Model_Forcing',values=ICE_FORCING_COMMENTS)
          IF(ASSOCIATED(ATT)) NCF => ADD(NCF,ATT)
       END IF
       
       
       IF(OBC_LONGSHORE_FLOW_ON) THEN
          aTT=> NC_MAKE_ATT(name='Special_Physical_processes',&
               & values='long shore flow adjustment for thermal wind and win&
               &d driven setup')
          NCF => ADD(NCF,ATT)
       END IF

    END IF

    ! ADD THE VARIABLES

    ! NPROCS
    VAR  => NC_MAKE_AVAR(name='nprocs',values=nprocs)

    ATT  => NC_MAKE_ATT(name='long_name',values='number of processors') 
    VAR  => ADD(VAR,ATT)
    NCF  => ADD(NCF,VAR)

#  if defined (MULTIPROCESSOR)
    ! PARTITION
    ALLOCATE(partition(G%NT)); partition=myid
    VAR  => NC_MAKE_PVAR(name='partition',values=partition,DIM1=DIM_nele)
    NULLIFY(PARTITION)
    ! THIS CAN CREATE A MEMORY LEAK - BE CAREFUL!
    

    ATT  => NC_MAKE_ATT(name='long_name',values='partition') 
    VAR  => ADD(VAR,ATT)

# if defined(UCF)
    ATT  => NC_MAKE_ATT(name='grid',values='grid1') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid_location',values='center') 
    VAR  => ADD(VAR,ATT)
# endif

    NCF  => ADD(NCF,VAR)
# endif



    ! X
    VAR  => NC_MAKE_AVAR(name='x',values=xm,DIM1=DIM_node)
   
    ATT  => NC_MAKE_ATT(name='long_name',values='nodal x-coordinate') 
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='units',values='meters') 
    VAR  => ADD(VAR,ATT)
    NCF  => ADD(NCF,VAR)
    
    ! Y
    VAR  => NC_MAKE_AVAR(name='y',values=ym,DIM1=DIM_node)

    ATT  => NC_MAKE_ATT(name='long_name',values='nodal y-coordinate') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='meters') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)


    ! LON
    VAR  => NC_MAKE_AVAR(name='lon',values=lon,DIM1=DIM_node)

    ATT  => NC_MAKE_ATT(name='long_name',values='nodal longitude') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='standard_name',values='longitude') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='degrees_east') 
    VAR  => ADD(VAR,ATT)
    NCF  => ADD(NCF,VAR)

    ! LAT
    VAR  => NC_MAKE_AVAR(name='lat',values=LAT,DIM1=DIM_node)

    ATT  => NC_MAKE_ATT(name='long_name',values='nodal latitude') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='standard_name',values='latitude') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='degrees_north') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    IF (ALLOCATED(XMC)) THEN
       ! XMC
       VAR  => NC_MAKE_AVAR(name='xc',values=xmc,DIM1=DIM_nele)
       
       ATT  => NC_MAKE_ATT(name='long_name',values='zonal x-coordinate') 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='units',values='meters') 
       VAR  => ADD(VAR,ATT)
       NCF  => ADD(NCF,VAR)
    END IF

    IF (ALLOCATED(YMC)) THEN
       ! YMC
       VAR  => NC_MAKE_AVAR(name='yc',values=ymc,DIM1=DIM_nele)
       
       ATT  => NC_MAKE_ATT(name='long_name',values='zonal y-coordinate') 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='units',values='meters') 
       VAR  => ADD(VAR,ATT)
       
       NCF  => ADD(NCF,VAR)
    END IF

    IF (ALLOCATED(LONC)) THEN
       ! LONC
       VAR  => NC_MAKE_AVAR(name='lonc',values=lonc,DIM1=DIM_nele)
       
       ATT  => NC_MAKE_ATT(name='long_name',values='zonal longitude') 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='standard_name',values='longitude') 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='units',values='degrees_east') 
       VAR  => ADD(VAR,ATT)
       NCF  => ADD(NCF,VAR)
    END IF

    IF (ALLOCATED(LATC)) THEN
       ! LATC
       VAR  => NC_MAKE_AVAR(name='latc',values=LATC,DIM1=DIM_nele)
       
       ATT  => NC_MAKE_ATT(name='long_name',values='zonal latitude') 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='standard_name',values='latitude') 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='units',values='degrees_north') 
       VAR  => ADD(VAR,ATT)
       
       NCF  => ADD(NCF,VAR)
    END IF

!# if defined(UCF)
    IF (ALLOCATED(h1)) THEN
       ! h1  (if morphodynamics is on, dump time-dependent bathymetry)
# if defined(SEDIMENT)
       IF(MORPHO_MODEL)THEN
         VAR  => NC_MAKE_AVAR(name='h_center',values=h1, DIM1= DIM_nele, DIM2= DIM_time)
       ELSE
         VAR  => NC_MAKE_AVAR(name='h_center',values=h1, DIM1= DIM_nele)
       ENDIF
# else
       VAR  => NC_MAKE_AVAR(name='h_center',values=h1, DIM1= DIM_nele)
# endif

       ATT  => NC_MAKE_ATT(name='long_name',values='Bathymetry') 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='standard_name',values='sea_floor_depth_below_geoid') 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='units',values='m')
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='positive',values='down') 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='grid',values='grid1 grid3') 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='coordinates',values='latc lonc') 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='grid_location',values='center') 
       VAR  => ADD(VAR,ATT)

       NCF  => ADD(NCF,VAR)
    END IF


!# endif

    IF (ALLOCATED(h)) THEN
       ! h (if morphodynamics is on, dump time-dependent bathymetry)
# if defined(SEDIMENT)
       IF(MORPHO_MODEL)THEN
         VAR  => NC_MAKE_AVAR(name='h_morpho',values=h, DIM1= DIM_node, DIM2= DIM_time)
         ATT  => NC_MAKE_ATT(name='long_name',values='Bathymetry_Morpho') 
         VAR  => ADD(VAR,ATT)
       
         ATT  => NC_MAKE_ATT(name='standard_name',values='sea_floor_depth_below_geoid') 
         VAR  => ADD(VAR,ATT)

         ATT  => NC_MAKE_ATT(name='units',values='m')
         VAR  => ADD(VAR,ATT)

         ATT  => NC_MAKE_ATT(name='positive',values='down') 
         VAR  => ADD(VAR,ATT)

         NCF  => ADD(NCF,VAR)
       END IF
       VAR  => NC_MAKE_AVAR(name='h',values=h, DIM1= DIM_node)

# else
       VAR  => NC_MAKE_AVAR(name='h',values=h, DIM1= DIM_node)
# endif

       ATT  => NC_MAKE_ATT(name='long_name',values='Bathymetry') 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='standard_name',values='sea_floor_depth_below_geoid') 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='units',values='m')
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='positive',values='down') 
       VAR  => ADD(VAR,ATT)

# if defined(UCF)
       ATT  => NC_MAKE_ATT(name='grid',values='grid1 grid3') 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='coordinates',values="lat lon") 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='grid_location',values='node') 
       VAR  => ADD(VAR,ATT)
# else
       ATT  => NC_MAKE_ATT(name='grid',values='Bathymetry_Mesh') 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar) 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='type',values='data') 
       VAR  => ADD(VAR,ATT)
# endif

       NCF  => ADD(NCF,VAR)
    END IF

# if defined(UCF)
    ! THE UCF STANDARD INCLUDES ONE CONNECTIVITY ARRAY FOR EACH MESH
    ! nv
    VAR  => NC_MAKE_PVAR(name='grid1',&
         & values=G%nv, DIM1= DIM_nele, DIM2= DIM_three)

    ATT  => NC_MAKE_ATT(name='cell_type',values='tri_ccw') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='index_start',values='0') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='standard_name',values='connectivity_array') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates_node',values='lon lat') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! nv
    VAR  => NC_MAKE_PVAR(name='grid2',&
         & values=G%nv, DIM1= DIM_nele, DIM2= DIM_three)

    ATT  => NC_MAKE_ATT(name='cell_type',values='tri_ccw') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='index_start',values='0') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='standard_name',values='connectivity_array') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates_node',values='lon lat zeta') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! nv
    VAR  => NC_MAKE_PVAR(name='grid3',&
         & values=G%nv, DIM1= DIM_nele, DIM2= DIM_three)

    ATT  => NC_MAKE_ATT(name='cell_type',values='tri_ccw') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='index_start',values='0') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='standard_name',values='connectivity_array') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates_node',values='lon lat h') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! nv
    VAR  => NC_MAKE_PVAR(name='grid4',&
         & values=G%nv, DIM1= DIM_nele, DIM2= DIM_three)

    ATT  => NC_MAKE_ATT(name='cell_type',values='tri_ccw') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='index_start',values='0') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='standard_name',values='connectivity_array') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates_node',values='lon lat siglay') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates_cell',values='lonc latc siglay_center') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! nv
    VAR  => NC_MAKE_PVAR(name='grid5',&
         & values=G%nv, DIM1= DIM_nele, DIM2= DIM_three)

    ATT  => NC_MAKE_ATT(name='cell_type',values='tri_ccw') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='index_start',values='0') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='standard_name',values='connectivity_array') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates_node',values='lon lat siglev') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates_cell',values='lonc latc siglev_center') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

# else
    ! nv
    VAR  => NC_MAKE_PVAR(name='nv',&
         & values=G%nv, DIM1= DIM_nele, DIM2= DIM_three)

    ATT  => NC_MAKE_ATT(name='long_name',values='nodes surrounding element') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)
# endif

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END GRID_FILE_OBJECT_SURFACE"


  END FUNCTION GRID_FILE_OBJECT_SURFACE

!=============================================================  
  FUNCTION ZETA_FILE_OBJECT() RESULT(NCF)
    IMPLICIT NONE
    
    INTEGER :: status, I
    LOGICAL, SAVE :: IOPROC_ALLOCATED = .FALSE.
    TYPE(NCFILE), POINTER :: NCF
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCATT),  POINTER :: ATT
    
    
    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START ZETA_FILE_OBJECT"
    
    ! IO PROC MUST ALLOCATE SPACE FOR THE ARRAYS 
    ! THESE ARRAYS MUST HAVE THE ATTRIBUTE SAVE AND FOR CLARITY
    ! SHOULD HAVE THE SAME NAME AS THOSE USED ON THE OTHER PROCESSORS!
    IF(IOPROC .AND. .NOT. IOPROC_ALLOCATED) THEN

       IOPROC_ALLOCATED = .true.

       allocate(EL(MGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:EL")
       EL = 0.0_SP

# if defined(UCF)
       allocate(EL1(NGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:EL1")
       EL1 = 0.0_SP
# endif

    END IF

    ! ALLOCATE THE NEW FILE OBJECT
    NCF => NEW_FILE()

  
    IF (ALLOCATED(EL)) THEN
       ! zeta
       VAR  => NC_MAKE_AVAR(name='zeta',&
            & values=EL, DIM1= DIM_node, DIM2= DIM_time)

       ATT  => NC_MAKE_ATT(name='long_name',values='Water Surface Elevation') 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='units',values='meters') 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='positive',values='up') 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='standard_name',values='sea_surface_height_above_geoid') 
       VAR  => ADD(VAR,ATT)
              
# if defined(UCF)
       ATT  => NC_MAKE_ATT(name='grid',values='grid1 grid2') 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='coordinates',values="time lat lon") 
       VAR  => ADD(VAR,ATT)
       
!!$       ATT  => NC_MAKE_ATT(name='grid_location',values='node') 
       ATT  => NC_MAKE_ATT(name='location',values='node') 
       VAR  => ADD(VAR,ATT)
# else
       ATT  => NC_MAKE_ATT(name='grid',values='Bathymetry_Mesh') 
       VAR  => ADD(VAR,ATT)

!!$       ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar) 
       ATT  => NC_MAKE_ATT(name='coordinates',values="time lat lon") 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='type',values='data') 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='location',values='node') 
       VAR  => ADD(VAR,ATT)
# endif
       
       NCF  => ADD(NCF,VAR)
    END IF

# if defined(UCF)
    IF (ALLOCATED(EL1)) THEN
       ! zeta
       VAR  => NC_MAKE_AVAR(name='zeta_center',&
            & values=EL1, DIM1= DIM_nele, DIM2= DIM_time)

       ATT  => NC_MAKE_ATT(name='long_name',values='Water Surface Elevation') 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='units',values='meters') 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='positive',values='up') 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='standard_name',values='sea_surface_height_above_geoid') 
       VAR  => ADD(VAR,ATT)
              
       ATT  => NC_MAKE_ATT(name='grid',values='grid1 grid2') 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='coordinates',values="time lat lon") 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='grid_location',values='center') 
       VAR  => ADD(VAR,ATT)
       
       NCF  => ADD(NCF,VAR)
    END IF
# endif

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END ZETA_FILE_OBJECT"

  END FUNCTION ZETA_FILE_OBJECT
!=============================================================  
  FUNCTION GRID_METRICS_FILE_OBJECT(G) RESULT(NCF)
    IMPLICIT NONE

    TYPE(GRID) :: G

    INTEGER :: status
    LOGICAL, SAVE :: IOPROC_ALLOCATED = .FALSE.
    TYPE(NCFILE), POINTER :: NCF
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCATT),  POINTER :: ATT

    REAL(SP), POINTER :: vec_flt(:),arr_flt(:)
    INTEGER, POINTER :: vec_int(:),arr_int(:)

    character(len=100)    :: timestamp, temp, netcdf_convention

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START GRID_METRICS_FILE_OBJECT"
    
    ! IO PROC MUST ALLOCATE SPACE FOR THE ARRAYS 
    ! THESE ARRAYS MUST HAVE THE ATTRIBUTE SAVE AND FOR CLARITY
    ! SHOULD HAVE THE SAME NAME AS THOSE USED ON THE OTHER PROCESSORS!

    IF(IOPROC) THEN
       ! ALWAYS ALLOCATE SPECIAL SPACE IN THE GRID TYPE FOR INDEX ARRAYS!

       IF(.NOT. ASSOCIATED(G%NBE)) THEN
          allocate(G%NBE(G%NGL,3),stat=status)
          IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:G%NBE")
          G%NBE = 0
       ELSE
          IF(UBOUND(G%NBE,1)/=DIM_NELE%DIM) CALL FATAL_ERROR &
               &("GRID DATA NBE HAS ALREADY BEEN ASSOICATED ON THE IOPROC",&
               & "AND THE DIMENSION DOES NOT MATCH THE CURRENT FILE!",&
               &"GRID NAME:"//TRIM(G%NAME))
       END IF
          
       IF(.NOT. ASSOCIATED(G%NBSN)) THEN
          allocate(G%NBSN(G%MGL,DIM_MaxNode%DIM),stat=status)
          IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:G%NBSN")
          G%NBSN = 0
       ELSE
          IF(UBOUND(G%NBSN,1)/=DIM_NODE%DIM) CALL FATAL_ERROR &
               &("GRID DATA NBSN HAS ALREADY BEEN ASSOICATED ON THE IOPROC",&
               & "AND THE DIMENSION DOES NOT MATCH THE CURRENT FILE!",&
               &"GRID NAME:"//TRIM(G%NAME))
       END IF

       IF(.NOT. ASSOCIATED(G%NBVE)) THEN
          allocate(G%NBVE(G%MGL,DIM_MaxElem%DIM),stat=status)
          IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:G%NBVE")
          G%NBVE = 0  
       ELSE
          IF(UBOUND(G%NBVE,1)/=DIM_NODE%DIM) CALL FATAL_ERROR &
               &("GRID DATA NBVE HAS ALREADY BEEN ASSOICATED ON THE IOPROC",&
               & "AND THE DIMENSION DOES NOT MATCH THE CURRENT FILE!",&
               &"GRID NAME:"//TRIM(G%NAME))
       END IF
       

    END IF


    IF(IOPROC .AND. .NOT. IOPROC_ALLOCATED) THEN

       ! ONLY ALLOCATE THIS SPACE ONCE - IT IS SHARED NO MATTER WHAT
       ! DOMAIN IS BEING SAVED!
       IOPROC_ALLOCATED = .true.

       allocate(NTSN(MGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:NTSN")
       NTSN = 0

       allocate(NTVE(MGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:NTVE")
       NTVE = 0

       allocate(A1U(NGL,4),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:A1U")
       A1U = 0.0_SP

       allocate(A2U(NGL,4),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:A2U")
       A2U = 0.0_SP

       allocate(AWX(NGL,3),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:AWX")
       AWX = 0.0_SP

       allocate(AWY(NGL,3),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:AWY")
       AWY = 0.0_SP

       allocate(AW0(NGL,3),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:AW0")
       AW0 = 0.0_SP

       allocate(ART2(MGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:ART2")
       ART2 = 0.0_SP

       allocate(ART1(MGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:ART1")
       ART1 = 0.0_SP


    END IF


  ! ALLOCATE THE NEW FILE OBJECT
    NCF => NEW_FILE()


    ! NBE
    VAR  => NC_MAKE_PVAR(name='nbe',&
         & values=G%NBE, DIM1= DIM_nele, DIM2= DIM_three)

    ATT  => NC_MAKE_ATT(name='long_name',values='elements surrounding each element') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! NTSN
    VAR  => NC_MAKE_AVAR(name='ntsn',&
         & values=NTSN, DIM1= DIM_node)

    ATT  => NC_MAKE_ATT(name='long_name',values='#nodes surrounding each node') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! NBSN
    VAR  => NC_MAKE_PVAR(name='nbsn',&
         & values=G%NBSN, DIM1= DIM_node, DIM2= DIM_MaxNode)

    ATT  => NC_MAKE_ATT(name='long_name',values='nodes surrounding each node') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! NTVE
    VAR  => NC_MAKE_AVAR(name='ntve',&
         & values=NTVE, DIM1= DIM_node)

    ATT  => NC_MAKE_ATT(name='long_name',values='#elems surrounding each node') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! NBVE
    VAR  => NC_MAKE_PVAR(name='nbve',&
         & values=G%NBVE, DIM1= DIM_node, DIM2= DIM_MaxElem)

    ATT  => NC_MAKE_ATT(name='long_name',values='elems surrounding each node') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! A1U
    VAR  => NC_MAKE_AVAR(name='a1u',&
         & values=a1u, DIM1= DIM_nele, DIM2= DIM_four)

    ATT  => NC_MAKE_ATT(name='long_name',values='a1u') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! A2U
    VAR  => NC_MAKE_AVAR(name='a2u',&
         & values=a2u, DIM1= DIM_nele, DIM2= DIM_four)

    ATT  => NC_MAKE_ATT(name='long_name',values='a2u') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! AW0
    VAR  => NC_MAKE_AVAR(name='aw0',&
         & values=aw0, DIM1= DIM_nele, DIM2= DIM_three)

    ATT  => NC_MAKE_ATT(name='long_name',values='aw0') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! AWX
    VAR  => NC_MAKE_AVAR(name='awx',&
         & values=awx, DIM1= DIM_nele, DIM2= DIM_three)

    ATT  => NC_MAKE_ATT(name='long_name',values='awx') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! AWY
    VAR  => NC_MAKE_AVAR(name='awy',&
         & values=awy, DIM1= DIM_nele, DIM2= DIM_three)

    ATT  => NC_MAKE_ATT(name='long_name',values='awy') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! ART2
    VAR  => NC_MAKE_AVAR(name='art2',&
         & values=art2, DIM1= DIM_node)

    ATT  => NC_MAKE_ATT(name='long_name',values='Area of elements around a node') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! ART1
    VAR  => NC_MAKE_AVAR(name='art1',&
         & values=art1, DIM1= DIM_node)

    ATT  => NC_MAKE_ATT(name='long_name',values='Area of Node-Base Con&
         &trol volume') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)


    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END GRID_METRICS_FILE_OBJECT"

  END FUNCTION GRID_METRICS_FILE_OBJECT

!=============================================================  
  FUNCTION FILE_DATE_OBJECT(SIZE) RESULT(NCF)
    USE MOD_CLOCK
    IMPLICIT NONE

    INTEGER :: status
    TYPE(NCFILE), POINTER :: NCF
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCATT),  POINTER :: ATT
    INTEGER, OPTIONAL :: SIZE
    CHARACTER(LEN=80),pointer :: Data_vec(:)
    
    
    IF(PRESENT(SIZE)) THEN
       ALLOCATE(DATA_vec(SIZE))
    ELSE
       ALLOCATE(DATA_vec(1))
    END IF


    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START FILE_DATE_OBJECT"
    
  ! ALLOCATE THE NEW FILE OBJECT
    NCF => NEW_FILE()

    ! FILE_DATE
    VAR  => NC_MAKE_PVAR(name='file_date', values=Data_vec, DIM1= DIM_DateStrLen, DIM2 = DIM_time)
    VAR%SCL_CHR => VAR%VEC_CHR(1)

    IF(USE_REAL_WORLD_TIME) THEN
       ATT  => NC_MAKE_ATT(name='time_zone',values=TRIM(TIMEZONE))
    ELSE
       ATT  => NC_MAKE_ATT(name='time_zone',values="UTC")
    END IF

    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

  END FUNCTION FILE_DATE_OBJECT


!=============================================================  
  FUNCTION VELOCITY_FILE_OBJECT() RESULT(NCF)
    IMPLICIT NONE

    INTEGER :: status
    LOGICAL, SAVE :: IOPROC_ALLOCATED = .FALSE.
    TYPE(NCFILE), POINTER :: NCF
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCATT),  POINTER :: ATT

    character(len=100)    :: timestamp, temp, netcdf_convention

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START VELOCITY_FILE_OBJECT"
    
    ! IO PROC MUST ALLOCATE SPACE FOR THE ARRAYS 
    ! THESE ARRAYS MUST HAVE THE ATTRIBUTE SAVE AND FOR CLARITY
    ! SHOULD HAVE THE SAME NAME AS THOSE USED ON THE OTHER PROCESSORS!
    IF(IOPROC .AND. .NOT. IOPROC_ALLOCATED) THEN

       IOPROC_ALLOCATED = .true.

       allocate(U(NGL,KB),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:U")
       U = 0.0_sp

       allocate(V(NGL,KB),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:V")
       V = 0.0_sp

       allocate(TAUBM(NGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:WUBOT")
       TAUBM   = 0.0_SP

    END IF


  ! ALLOCATE THE NEW FILE OBJECT
    NCF => NEW_FILE()


    ! U
    VAR  => NC_MAKE_AVAR(name='u',&
         & values=U, DIM1= DIM_nele, DIM2= DIM_siglay, DIM3 = DIM_time)

# if defined(UCF)
    ATT  => NC_MAKE_ATT(name='long_name',values='East Velocity Component') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='standard_name',values='eastward_sea_water_velocity') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='m s-1') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='grid5') 
    VAR  => ADD(VAR,ATT)

!!$    ATT  => NC_MAKE_ATT(name='grid_location',values='center') 
    ATT  => NC_MAKE_ATT(name='location',values='face') 
    VAR  => ADD(VAR,ATT)

!!$    ATT  => NC_MAKE_ATT(name='coordinates',values='latc lonc siglay_center') 
    ATT  => NC_MAKE_ATT(name='coordinates',values='time siglay latc lonc') 
    VAR  => ADD(VAR,ATT)

# else
    ATT  => NC_MAKE_ATT(name='long_name',values='Eastward Water Velocity') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='standard_name',values='eastward_sea_water_velocity') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='meters s-1') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates',values='time siglay latc lonc') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='mesh',values='fvcom_mesh') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='location',values='face') 
    VAR  => ADD(VAR,ATT)
# endif

    NCF  => ADD(NCF,VAR)



    ! V
    VAR  => NC_MAKE_AVAR(name='v',&
         & values=V, DIM1= DIM_nele, DIM2= DIM_siglay, DIM3 = DIM_time)

# if defined(UCF)
    ATT  => NC_MAKE_ATT(name='long_name',values='North Velocity Component') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='standard_name',values='Northward_sea_water_velocity') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='m s-1') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='grid5') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid_location',values='center') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates',values='latc lonc siglay_center') 
    VAR  => ADD(VAR,ATT)

# else
    ATT  => NC_MAKE_ATT(name='long_name',values='Northward Water Velocity') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='standard_name',values='Northward_sea_water_velocity') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='meters s-1') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates',values='time siglay latc lonc') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='mesh',values='fvcom_mesh') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='location',values='face') 
    VAR  => ADD(VAR,ATT)
# endif

    NCF  => ADD(NCF,VAR)

  VAR  => NC_MAKE_AVAR(name='tauc',&
       & values=taubm, DIM1= DIM_NELE, DIM2= DIM_time)

  ATT  => NC_MAKE_ATT(name='long_name',values='bed stress magnitude from currents')
  VAR  => ADD(VAR,ATT)

  ATT  => NC_MAKE_ATT(name='note1',values='this stress is bottom boundary condtion on velocity field')
  VAR  => ADD(VAR,ATT)

  ATT  => NC_MAKE_ATT(name='note2',values='dimensions are stress/rho')
  VAR  => ADD(VAR,ATT)

  ATT  => NC_MAKE_ATT(name='units',values='m^2 s^-2')
  VAR  => ADD(VAR,ATT)

  ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid')
  VAR  => ADD(VAR,ATT)

  ATT  => NC_MAKE_ATT(name='type',values='data')
  VAR  => ADD(VAR,ATT)

  ATT  => NC_MAKE_ATT(name='coordinates',values='time latc lonc')
  VAR  => ADD(VAR,ATT)
  
  ATT  => NC_MAKE_ATT(name='mesh',values='fvcom_mesh')
  VAR  => ADD(VAR,ATT)

  ATT  => NC_MAKE_ATT(name='location',values='face') 
  VAR  => ADD(VAR,ATT)

  NCF  => ADD(NCF,VAR)


    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END VELOCITY_FILE_OBJECT"

  END FUNCTION VELOCITY_FILE_OBJECT
!=============================================================  
  FUNCTION VORTICITY_FILE_OBJECT() RESULT(NCF)
    IMPLICIT NONE

    INTEGER :: status
    LOGICAL, SAVE :: IOPROC_ALLOCATED = .FALSE.
    TYPE(NCFILE), POINTER :: NCF
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCATT),  POINTER :: ATT

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START VORTIICTY_FILE_OBJECT"
    
    ! IO PROC MUST ALLOCATE SPACE FOR THE ARRAYS 
    ! THESE ARRAYS MUST HAVE THE ATTRIBUTE SAVE AND FOR CLARITY
    ! SHOULD HAVE THE SAME NAME AS THOSE USED ON THE OTHER PROCESSORS!
    IF(IOPROC .AND. .NOT. IOPROC_ALLOCATED) THEN

       IOPROC_ALLOCATED = .true.

       allocate(VORT(MGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:VORT")
       VORT = 0.0_SP

    END IF


  ! ALLOCATE THE NEW FILE OBJECT
    NCF => NEW_FILE()


    ! VORTICITY
    VAR  => NC_MAKE_AVAR(name='vorticity',&
         & values=VORT, DIM1= DIM_node, DIM2= DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Ertels 2d potential vorticity') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='s-1') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)



    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END VORTICITY_FILE_OBJECT"

  END FUNCTION VORTICITY_FILE_OBJECT
!=============================================================  
  FUNCTION AVERAGE_VEL_FILE_OBJECT() RESULT(NCF)
    IMPLICIT NONE

    INTEGER :: status
    LOGICAL, SAVE :: IOPROC_ALLOCATED = .FALSE.
    TYPE(NCFILE), POINTER :: NCF
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCATT),  POINTER :: ATT

    character(len=100)    :: timestamp, temp, netcdf_convention

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START AVERAGE_VEL_FILE_OBJECT"
    
    ! IO PROC MUST ALLOCATE SPACE FOR THE ARRAYS 
    ! THESE ARRAYS MUST HAVE THE ATTRIBUTE SAVE AND FOR CLARITY
    ! SHOULD HAVE THE SAME NAME AS THOSE USED ON THE OTHER PROCESSORS!
    IF(IOPROC .AND. .NOT. IOPROC_ALLOCATED) THEN

       IOPROC_ALLOCATED = .true.

       allocate(UA(NGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:UA")
       UA = 0.0_SP

       allocate(VA(NGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:VA")
       VA = 0.0_SP

    END IF


  ! ALLOCATE THE NEW FILE OBJECT
    NCF => NEW_FILE()


    ! UA
    VAR  => NC_MAKE_AVAR(name='ua',&
         & values=UA, DIM1= DIM_nele, DIM2= DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Vertically Averaged x-velocity') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='meters s-1') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)


    ! VA
    VAR  => NC_MAKE_AVAR(name='va',&
         & values=VA, DIM1= DIM_nele, DIM2= DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Vertically Averaged y-velocity') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='meters s-1') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)


    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END AVERAGE_VEL_FILE_OBJECT"

  END FUNCTION AVERAGE_VEL_FILE_OBJECT
  
!=============================================================  
  FUNCTION VERTICAL_VEL_FILE_OBJECT() RESULT(NCF)
    IMPLICIT NONE

    INTEGER :: status
    LOGICAL, SAVE :: IOPROC_ALLOCATED = .FALSE.
    TYPE(NCFILE), POINTER :: NCF
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCATT),  POINTER :: ATT

    character(len=100)    :: timestamp, temp, netcdf_convention

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START VERTICAL_VEL_FILE_OBJECT"
    
    ! IO PROC MUST ALLOCATE SPACE FOR THE ARRAYS 
    ! THESE ARRAYS MUST HAVE THE ATTRIBUTE SAVE AND FOR CLARITY
    ! SHOULD HAVE THE SAME NAME AS THOSE USED ON THE OTHER PROCESSORS!
    IF(IOPROC .AND. .NOT. IOPROC_ALLOCATED) THEN

       IOPROC_ALLOCATED = .true.

       allocate(WW(NGL,KB),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:WW")
       WW = 0.0_SP

       allocate(WTS(MGL,KB),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:WTS")
       WTS = 0.0_SP

    END IF


  ! ALLOCATE THE NEW FILE OBJECT
    NCF => NEW_FILE()


    ! WTS
    VAR  => NC_MAKE_AVAR(name='omega',&
         & values=WTS, DIM1= DIM_node, DIM2= DIM_siglev, DIM3 = DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Vertical Sigma Coordinate Velocity') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='s-1') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)



    ! WW
    VAR  => NC_MAKE_AVAR(name='ww',&
         & values=WW, DIM1= DIM_nele, DIM2= DIM_siglay, DIM3 = DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Upward Water Velocity') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='meters s-1') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END VERTICAL_VEL_FILE_OBJECT"

  END FUNCTION VERTICAL_VEL_FILE_OBJECT
!=============================================================  
  FUNCTION SALT_TEMP_FILE_OBJECT() RESULT(NCF)
    IMPLICIT NONE

    INTEGER :: status
    LOGICAL, SAVE :: IOPROC_ALLOCATED = .FALSE.
    TYPE(NCFILE), POINTER :: NCF
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCATT),  POINTER :: ATT

    character(len=100)    :: timestamp, temp, netcdf_convention

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START SALT_TEMP_FILE_OBJECT"
    
    ! IO PROC MUST ALLOCATE SPACE FOR THE ARRAYS 
    ! THESE ARRAYS MUST HAVE THE ATTRIBUTE SAVE AND FOR CLARITY
    ! SHOULD HAVE THE SAME NAME AS THOSE USED ON THE OTHER PROCESSORS!
    IF(IOPROC .AND. .NOT. IOPROC_ALLOCATED) THEN

       IOPROC_ALLOCATED = .true.

       allocate(T1(MGL,KB),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:T1")
       T1 = 0.0_SP

       allocate(S1(MGL,KB),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:S1")
       S1 = 0.0_SP

    END IF


  ! ALLOCATE THE NEW FILE OBJECT
    NCF => NEW_FILE()


    ! T
    VAR  => NC_MAKE_AVAR(name='temp',&
         & values=T1, DIM1= DIM_node, DIM2= DIM_siglay, DIM3 = DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='temperature') 
    VAR  => ADD(VAR,ATT)

# if defined(UCF)
    ATT  => NC_MAKE_ATT(name='standard_name',values='sea_water_potential_temperature') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='degC') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='grid4') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates',values='lon lat siglay') 
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='grid_location',values='node') 
    VAR  => ADD(VAR,ATT)    

# else
    ATT  => NC_MAKE_ATT(name='standard_name',values='sea_water_temperature') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='degrees_C') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

!!$    ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar) 
    ATT  => NC_MAKE_ATT(name='coordinates',values='time siglay lat lon') 
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)    

    ATT  => NC_MAKE_ATT(name='mesh',values='fvcom_mesh') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='location',values='node') 
    VAR  => ADD(VAR,ATT)
# endif

    NCF  => ADD(NCF,VAR)


    ! S
    VAR  => NC_MAKE_AVAR(name='salinity',&
         & values=S1, DIM1= DIM_node, DIM2= DIM_siglay, DIM3 = DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='salinity') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='standard_name',values='sea_water_salinity') 
    VAR  => ADD(VAR,ATT)

# if defined(UCF)
    ATT  => NC_MAKE_ATT(name='units',values='1e-3') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='grid4') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates',values='lon lat siglay') 
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='grid_location',values='node') 
    VAR  => ADD(VAR,ATT)    

# else

    ATT  => NC_MAKE_ATT(name='units',values='1e-3') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

!!$    ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar) 
    ATT  => NC_MAKE_ATT(name='coordinates',values='time siglay lat lon') 
    VAR  => ADD(VAR,ATT)    

    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='mesh',values='fvcom_mesh') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='location',values='node') 
    VAR  => ADD(VAR,ATT)
# endif

    NCF  => ADD(NCF,VAR)


    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END SALT_TEMP_FILE_OBJECT"

  END FUNCTION SALT_TEMP_FILE_OBJECT
!=============================================================  
  FUNCTION DENSITY_FILE_OBJECT() RESULT(NCF)
    IMPLICIT NONE

    INTEGER :: status
    LOGICAL, SAVE :: IOPROC_ALLOCATED = .FALSE.
    TYPE(NCFILE), POINTER :: NCF
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCATT),  POINTER :: ATT

    character(len=100)    :: timestamp, temp, netcdf_convention

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START DENSITY_FILE_OBJECT"
    
    ! IO PROC MUST ALLOCATE SPACE FOR THE ARRAYS 
    ! THESE ARRAYS MUST HAVE THE ATTRIBUTE SAVE AND FOR CLARITY
    ! SHOULD HAVE THE SAME NAME AS THOSE USED ON THE OTHER PROCESSORS!
    IF(IOPROC .AND. .NOT. IOPROC_ALLOCATED) THEN

       IOPROC_ALLOCATED = .true.

       allocate(RHO1(MGL,KB),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:RHO1")
       RHO1 = 0.0_SP

       allocate(RMEAN1(MGL,KB),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:RMEAN1")
       RMEAN1 = 0.0_SP

    END IF


  ! ALLOCATE THE NEW FILE OBJECT
    NCF => NEW_FILE()


    ! DENSITY
    VAR  => NC_MAKE_AVAR(name='rho1',&
         & values=RHO1, DIM1= DIM_node, DIM2= DIM_siglay, DIM3 = DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='density') 
    VAR  => ADD(VAR,ATT)

# if defined(UCF)
    ATT  => NC_MAKE_ATT(name='standard_name',values='sea_water_potential_density') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='kg/m3') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='grid4') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates',values='lon lat siglay') 
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='grid_location',values='node') 
    VAR  => ADD(VAR,ATT)    

# else
    ATT  => NC_MAKE_ATT(name='standard_name',values='sea_water_density') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='kg/m3') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar) 
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)    
# endif

    NCF  => ADD(NCF,VAR)


    ! RMEAN1
    VAR  => NC_MAKE_AVAR(name='rmean1',&
         & values=RMEAN1, DIM1= DIM_node, DIM2= DIM_siglay)

    ATT  => NC_MAKE_ATT(name='long_name',values='mean density') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='standard_name',values='sea_water_density') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='kg/m3') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='SigmaLayer_Mesh') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)


    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END DENSITY_FILE_OBJECT"

  END FUNCTION DENSITY_FILE_OBJECT
!=============================================================  
  FUNCTION TURBULENCE_FILE_OBJECT() RESULT(NCF)
    IMPLICIT NONE

    INTEGER :: status
    LOGICAL, SAVE :: IOPROC_ALLOCATED = .FALSE.
    TYPE(NCFILE), POINTER :: NCF
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCATT),  POINTER :: ATT
    TYPE(NCDIM),  POINTER :: DIM1
    TYPE(NCDIM),  POINTER :: DIM2
    TYPE(NCDIM),  POINTER :: DIM3

    character(len=100)    :: timestamp, temp, netcdf_convention

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START TURBULENCE_FILE_OBJECT"
    
    ! IO PROC MUST ALLOCATE SPACE FOR THE ARRAYS 
    ! THESE ARRAYS MUST HAVE THE ATTRIBUTE SAVE AND FOR CLARITY
    ! SHOULD HAVE THE SAME NAME AS THOSE USED ON THE OTHER PROCESSORS!
    IF(IOPROC .AND. .NOT. IOPROC_ALLOCATED) THEN

       IOPROC_ALLOCATED = .true.

       allocate(VISCOFM(NGL,KB),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:VISCOFM")
       VISCOFM = 0.0_SP

       allocate(VISCOFH(MGL,KB),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:VISCOFH")
       VISCOFH = 0.0_SP
       
       allocate(KM(MGL,KB),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:KM")
       KM = 0.0_SP

       allocate(KH(MGL,KB),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:KH")
       KH = 0.0_SP

       allocate(KQ(MGL,KB),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:KQ")
       KQ = 0.0_SP
       
# if defined (GOTM)
       allocate(TKE(MGL,KB),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:TKE")
       TKE = 0.0_SP

       allocate(TEPS(MGL,KB),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:TEPS")
       TEPS = 0.0_SP
# else
       allocate(Q2(MGL,KB),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:Q2")
       Q2 = 0.0_SP

       allocate(Q2L(MGL,KB),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:Q2L")
       Q2L = 0.0_SP
       
       allocate(L(MGL,KB),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:L")
       L = 0.0_SP
# endif


    END IF


  ! ALLOCATE THE NEW FILE OBJECT
    NCF => NEW_FILE()


    ! VISCOFM
    VAR  => NC_MAKE_AVAR(name='viscofm',&
         & values=viscofm, DIM1= DIM_nele, DIM2= DIM_siglay, DIM3 = DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Horizontal Turbulent Eddy Viscosity For Momentum') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='m 2 s-1') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar) 
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! VISCOFH
    VAR  => NC_MAKE_AVAR(name='viscofh',&
         & values=viscofh, DIM1= DIM_node, DIM2= DIM_siglay, DIM3 = DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Horizontal Turbulent Eddy Viscosity For Scalars') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='m 2 s-1') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar) 
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! KM
    VAR  => NC_MAKE_AVAR(name='km',&
         & values=km, DIM1= DIM_node, DIM2= DIM_siglev, DIM3 = DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Turbulent Eddy Viscosity For Momentum') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='m 2 s-1') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar) 
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! KH
    VAR  => NC_MAKE_AVAR(name='kh',&
         & values=kh, DIM1= DIM_node, DIM2= DIM_siglev, DIM3 = DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Turbulent Eddy Viscosity For Scalars') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='m 2 s-1') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar) 
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! KQ
    VAR  => NC_MAKE_AVAR(name='kq',&
         & values=kq, DIM1= DIM_node, DIM2= DIM_siglev, DIM3 = DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Turbulent Eddy Viscosity For Q2/Q2L') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='m 2 s-1') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar) 
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

# if defined (GOTM)
    ! TKE
    VAR  => NC_MAKE_AVAR(name='tke',&
         & values=tke, DIM1= DIM_node, DIM2= DIM_siglev, DIM3 = DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Turbulent Kinetic Energy') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='m^2/s^2') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar) 
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! TEPS
    VAR  => NC_MAKE_AVAR(name='teps',&
         & values=teps, DIM1= DIM_node, DIM2= DIM_siglev, DIM3 = DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Turbulent Kinetic Energy Dissipation Rate') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='m^2/s^3') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar) 
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)
# else
    ! Q2
    VAR  => NC_MAKE_AVAR(name='q2',&
         & values=q2, DIM1= DIM_node, DIM2= DIM_siglev, DIM3 = DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Turbulent Kinetic Energy') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='m2 s-2') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar) 
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! Q2L
    VAR  => NC_MAKE_AVAR(name='q2l',&
         & values=q2l, DIM1= DIM_node, DIM2= DIM_siglev, DIM3 = DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Turbulent Kinetic Ene&
         &rgy X Turbulent Macroscale') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='m3 s-2') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar) 
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! L
    VAR  => NC_MAKE_AVAR(name='l',&
         & values=l, DIM1= DIM_node, DIM2= DIM_siglev, DIM3 = DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Turbulent Macroscale') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='m3 s-2') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar) 
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

# endif

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END TURBULENCE_FILE_OBJECT"

  END FUNCTION TURBULENCE_FILE_OBJECT
!=============================================================  
  FUNCTION SURFACE_HEATING_FILE_OBJECT() RESULT(NCF)
    IMPLICIT NONE

    INTEGER :: status
    LOGICAL, SAVE :: IOPROC_ALLOCATED = .FALSE.
    TYPE(NCFILE), POINTER :: NCF
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCATT),  POINTER :: ATT

    character(len=100)    :: timestamp, temp, netcdf_convention

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START SURFACE_HEATING_FILE_OBJECT"
    
    ! IO PROC MUST ALLOCATE SPACE FOR THE ARRAYS 
    ! THESE ARRAYS MUST HAVE THE ATTRIBUTE SAVE AND FOR CLARITY
    ! SHOULD HAVE THE SAME NAME AS THOSE USED ON THE OTHER PROCESSORS!
    IF(IOPROC .AND. .NOT. IOPROC_ALLOCATED) THEN

       IOPROC_ALLOCATED = .true.

       allocate(SWRAD_WATTS(MGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:SWRAD")
       swrad_watts = 0.0_SP

       allocate(WTSURF_WATTS(MGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:WTSURF")
       WTSURF_watts = 0.0_SP

#      if defined (HEATING_CALCULATED) || (HEATING_SOLAR) || (ICE)
       allocate(HSENS_WATTS(MGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:HSENS")
       HSENS_watts = 0.0_SP

       allocate(HLAT_WATTS(MGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:HLAT")
       HLAT_watts = 0.0_SP

       allocate(LWRAD_WATTS(MGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:LWRAD")
       LWRAD_watts = 0.0_SP
#      endif

    END IF


  ! ALLOCATE THE NEW FILE OBJECT
    NCF => NEW_FILE()


    ! SWRAD
    VAR  => NC_MAKE_AVAR(name='short_wave',&
         & values=swrad_WATTS, DIM1= DIM_node, DIM2= DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Short Wave Radiation') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='W m-2') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

!!$    ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar) 
    ATT  => NC_MAKE_ATT(name='coordinates',values='time lat lon') 
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='mesh',values='fvcom_mesh') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='location',values='node') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! WTSURF - NET HEAT FLUX
    VAR  => NC_MAKE_AVAR(name='net_heat_flux',&
         & values=WTSURF_WATTS, DIM1= DIM_node, DIM2= DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Surface Net Heat Flux') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='W m-2') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

!!$    ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar) 
    ATT  => NC_MAKE_ATT(name='coordinates',values='time lat lon') 
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='mesh',values='fvcom_mesh') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='location',values='node') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

#   if defined (HEATING_CALCULATED) || (HEATING_SOLAR) || (ICE)
    !EJA edit to add Sensible, Latent, and Longwave radiation to output - 05/03/2016
    ! HSENS
    VAR  => NC_MAKE_AVAR(name='sensible_heat_flux',&
         & values=HSENS_WATTS, DIM1= DIM_node, DIM2= DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Sensible Heat Flux')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='W m-2')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid')
    VAR  => ADD(VAR,ATT)

!!$    ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar)
    ATT  => NC_MAKE_ATT(name='coordinates',values='time lat lon')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='mesh',values='fvcom_mesh')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='location',values='node')
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! HLAT
    VAR  => NC_MAKE_AVAR(name='latent_heat_flux',&
         & values=HLAT_WATTS, DIM1= DIM_node, DIM2= DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Latent Heat Flux')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='W m-2')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid')
    VAR  => ADD(VAR,ATT)

!!$    ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar)
    ATT  => NC_MAKE_ATT(name='coordinates',values='time lat lon')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='mesh',values='fvcom_mesh')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='location',values='node')
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! LWRAD
    VAR  => NC_MAKE_AVAR(name='long_wave',&
         & values=LWRAD_WATTS, DIM1= DIM_node, DIM2= DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Long Wave Radiation')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='W m-2')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid')
    VAR  => ADD(VAR,ATT)

!!$    ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar)
    ATT  => NC_MAKE_ATT(name='coordinates',values='time lat lon')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='mesh',values='fvcom_mesh')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='location',values='node')
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)
#   endif

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END SURFACE_HEATING_FILE_OBJECT"

  END FUNCTION SURFACE_HEATING_FILE_OBJECT
!=============================================================  
  FUNCTION WIND_VELOCITY_FILE_OBJECT() RESULT(NCF)
    IMPLICIT NONE

    INTEGER :: status
    LOGICAL, SAVE :: IOPROC_ALLOCATED = .FALSE.
    TYPE(NCFILE), POINTER :: NCF
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCATT),  POINTER :: ATT

    character(len=100)    :: timestamp, temp, netcdf_convention

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START WIND_VELOCITY_FILE_OBJECT"
    
    ! IO PROC MUST ALLOCATE SPACE FOR THE ARRAYS 
    ! THESE ARRAYS MUST HAVE THE ATTRIBUTE SAVE AND FOR CLARITY
    ! SHOULD HAVE THE SAME NAME AS THOSE USED ON THE OTHER PROCESSORS!
    IF(IOPROC .AND. .NOT. IOPROC_ALLOCATED) THEN

       IOPROC_ALLOCATED = .true.

       allocate(UUWIND(NGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:UUWIND")
       UUWIND = 0.0_SP

       allocate(VVWIND(NGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:VVWIND")
       VVWIND = 0.0_SP

    END IF


  ! ALLOCATE THE NEW FILE OBJECT
    NCF => NEW_FILE()


    ! UUWIND
    VAR  => NC_MAKE_AVAR(name='uwind_speed',&
         & values=uuwind, DIM1= DIM_nele, DIM2= DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Eastward Wind Velocity') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='standard_name',values='eastward wind') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='meters s-1') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates',values='time latc lonc') 
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='mesh',values='fvcom_mesh') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='location',values='face') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! VVWIND
    VAR  => NC_MAKE_AVAR(name='vwind_speed',&
         & values=vvwind, DIM1= DIM_nele, DIM2= DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Northward Wind Velocity') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='standard_name',values='northward wind') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='meters s-1') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates',values='time latc lonc') 
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='mesh',values='fvcom_mesh') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='location',values='face') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)


    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END WIND_VELOCITY_FILE_OBJECT"
    
  END FUNCTION WIND_VELOCITY_FILE_OBJECT
!=============================================================  
# if defined (AIR_PRESSURE) || (HEATING_CALCULATED)
  FUNCTION ATMOSPHERIC_PRESSURE_FILE_OBJECT() RESULT(NCF)
    IMPLICIT NONE

    INTEGER :: status
    LOGICAL, SAVE :: IOPROC_ALLOCATED = .FALSE.
    TYPE(NCFILE), POINTER :: NCF
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCATT),  POINTER :: ATT

    character(len=100)    :: timestamp, temp, netcdf_convention

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START ATMOSPHERIC_PRESSURE_FILE_OBJECT"
    
    ! IO PROC MUST ALLOCATE SPACE FOR THE ARRAYS 
    ! THESE ARRAYS MUST HAVE THE ATTRIBUTE SAVE AND FOR CLARITY
    ! SHOULD HAVE THE SAME NAME AS THOSE USED ON THE OTHER PROCESSORS!
    IF(IOPROC .AND. .NOT. IOPROC_ALLOCATED) THEN

       IOPROC_ALLOCATED = .true.

       allocate(PA_AIR(MGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:PA_AIR")
       pa_air = 0.0_SP

    END IF


  ! ALLOCATE THE NEW FILE OBJECT
    NCF => NEW_FILE()


    ! PA_AIR
    VAR  => NC_MAKE_AVAR(name='atmos_press',&
         & values=pa_air, DIM1= DIM_node, DIM2= DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Atmospheric Pressure') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='pascals') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar) 
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END ATMOSPHERIC_PRESSURE_FILE_OBJECT"

  END FUNCTION ATMOSPHERIC_PRESSURE_FILE_OBJECT
# endif
!=============================================================  
  FUNCTION WIND_STRESS_FILE_OBJECT() RESULT(NCF)
    IMPLICIT NONE

    INTEGER :: status
    LOGICAL, SAVE :: IOPROC_ALLOCATED = .FALSE.
    TYPE(NCFILE), POINTER :: NCF
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCATT),  POINTER :: ATT

    character(len=100)    :: timestamp, temp, netcdf_convention

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START WIND_STRESS_FILE_OBJECT"
    
    ! IO PROC MUST ALLOCATE SPACE FOR THE ARRAYS 
    ! THESE ARRAYS MUST HAVE THE ATTRIBUTE SAVE AND FOR CLARITY
    ! SHOULD HAVE THE SAME NAME AS THOSE USED ON THE OTHER PROCESSORS!
    IF(IOPROC .AND. .NOT. IOPROC_ALLOCATED) THEN

       IOPROC_ALLOCATED = .true.

       allocate(WUSURF_save(NGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:WUSURF")
       WUSURF_save = 0.0_SP

       allocate(WVSURF_save(NGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:WVSURF")
       WVSURF_save = 0.0_SP

    END IF


  ! ALLOCATE THE NEW FILE OBJECT
    NCF => NEW_FILE()


    ! UUWIND
    VAR  => NC_MAKE_AVAR(name='uwind_stress',&
         & values=WUSURF_save, DIM1= DIM_nele, DIM2= DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Eastward Wind Stress') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='standard_name',values='Wind Stress') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='Pa') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! VVWIND
    VAR  => NC_MAKE_AVAR(name='vwind_stress',&
         & values=WVSURF_save, DIM1= DIM_nele, DIM2= DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Northward Wind Stress') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='standard_name',values='Wind Stress') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='Pa') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)


    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END WIND_STRESS_FILE_OBJECT"
    
  END FUNCTION WIND_STRESS_FILE_OBJECT
!=============================================================  
!=============================================================  
  FUNCTION PRECIPITATION_FILE_OBJECT() RESULT(NCF)
    IMPLICIT NONE

    INTEGER :: status
    LOGICAL, SAVE :: IOPROC_ALLOCATED = .FALSE.
    TYPE(NCFILE), POINTER :: NCF
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCATT),  POINTER :: ATT
    TYPE(NCDIM),  POINTER :: DIM1
    TYPE(NCDIM),  POINTER :: DIM2
    TYPE(NCDIM),  POINTER :: DIM3

    character(len=100)    :: timestamp, temp, netcdf_convention

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START PRECIPITATION_FILE_OBJECT"
    
    ! IO PROC MUST ALLOCATE SPACE FOR THE ARRAYS 
    ! THESE ARRAYS MUST HAVE THE ATTRIBUTE SAVE AND FOR CLARITY
    ! SHOULD HAVE THE SAME NAME AS THOSE USED ON THE OTHER PROCESSORS!
    IF(IOPROC .AND. .NOT. IOPROC_ALLOCATED) THEN

       IOPROC_ALLOCATED = .true.

       allocate(QPREC(MGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:QPREC2")
       QPREC = 0.0_SP

       allocate(QEVAP(MGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:QEVAP2")
       QEVAP = 0.0_SP

    END IF


  ! ALLOCATE THE NEW FILE OBJECT
    NCF => NEW_FILE()


    ! PRECIPITATION
    VAR  => NC_MAKE_AVAR(name='precip',&
         & values=qprec, DIM1= DIM_node, DIM2= DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Precipitation') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='description',values='Precipitation, ocean &
         &lose water is negative') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='m s-1') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar) 
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! EVAPORATION
    VAR  => NC_MAKE_AVAR(name='evap',&
         & values=QEVAP, DIM1= DIM_node, DIM2= DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Evaporation') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='description',values='Evaporation, ocean &
         &lose water is negative') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='m s-1') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar) 
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)


    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END PRECIPITATION_FILE_OBJECT"
    
  END FUNCTION PRECIPITATION_FILE_OBJECT

!=============================================================  
!=============================================================  
!=============================================================  
  FUNCTION RESTART_EXTRAS_FILE_OBJECT() RESULT(NCF)
    USE BCS
    IMPLICIT NONE

    INTEGER :: status
    LOGICAL, SAVE :: IOPROC_ALLOCATED = .FALSE.
    TYPE(NCFILE), POINTER :: NCF
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCATT),  POINTER :: ATT

    character(len=100)    :: timestamp, temp, netcdf_convention

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START RESTART_EXTRAS_FILE_OBJECT"
    
    ! IO PROC MUST ALLOCATE SPACE FOR THE ARRAYS 
    ! THESE ARRAYS MUST HAVE THE ATTRIBUTE SAVE AND FOR CLARITY
    ! SHOULD HAVE THE SAME NAME AS THOSE USED ON THE OTHER PROCESSORS!
    IF(IOPROC .AND. .NOT. IOPROC_ALLOCATED) THEN

       IOPROC_ALLOCATED = .true.

       allocate(COR(NGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:COR")
       COR = 0.0_SP

       allocate(CC_SPONGE(NGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:CC_SPONGE")
       CC_SPONGE = 0.0_SP

       allocate(ET(MGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:EL")
       ET = 0.0_SP


! OPEN BOUNDARY SETTINGS
       IF (OBC_ON) THEN
          ! NEED SPECIAL VARIABLE FOR GLOBAL NODE NUMBER OF LOCAL NODES...
          allocate(I_OBC_N_OUTPUT(IOBCN_GL),stat=status)
          IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:I_OBC_N_OUTPUT")
          I_OBC_N_OUTPUT = 0
          
          
          allocate(type_obc(IOBCN_GL),stat=status)
          IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:type_obc")
          type_OBC = 0
       END IF

! LONG SHORE FLOW BOUNDARY SETTINGS
       IF(OBC_LONGSHORE_FLOW_ON) THEN
          allocate(IBCLSF_OUTPUT(NOBCLSF_GL),stat=status)
          IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:IBCLSF_OUTPUT")
          IBCLSF_OUTPUT = 0
          
          allocate(RBC_GEO(NOBCLSF_GL),stat=status)
          IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:RBC_GEO")
          RBC_GEO = 0.0_SP
          
          allocate(RBC_WDF(NOBCLSF_GL),stat=status)
          IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:RBC_WDF")
          RBC_WDF = 0.0_SP
       END IF
       
          
       allocate(TMEAN1(MGL,KB),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:Tmean1")
       TMEAN1 = 0.0_SP
       
       allocate(SMEAN1(MGL,KB),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:Smean1")
       SMEAN1 = 0.0_SP

# if defined (EQUI_TIDE)
       ALLOCATE(EL_EQI(MGL), stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:EL_EQI")
       EL_EQI = 0.0_SP
# endif

# if defined (ATMO_TIDE)
       ALLOCATE(EL_ATMO(MGL), stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:EL_ATMO")
       EL_ATMO = 0.0_SP
# endif

    ! ------ New: Karsten Lettmann, 2016, June ---------------
    ! add the atmospheric pressure elevation
# if defined(AIR_PRESSURE)
       ALLOCATE(EL_AIR(MGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:EL_AIR")
       EL_AIR = 0.0_SP
# endif      
    ! -------- end new -------------------------------------


    END IF

  ! ALLOCATE THE NEW FILE OBJECT
    NCF => NEW_FILE()


    ! COR
    VAR  => NC_MAKE_AVAR(name='cor',&
         & values=COR, DIM1= DIM_nele)

    ATT  => NC_MAKE_ATT(name='long_name',values='Coriolis Parameter') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='s-1') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! CC_SPONGE
    VAR  => NC_MAKE_AVAR(name='cc_sponge',&
         & values=cc_sponge, DIM1= DIM_nele)

    ATT  => NC_MAKE_ATT(name='long_name',values='Sponge Layer Parameter') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='nd') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! et
    VAR  => NC_MAKE_AVAR(name='et',&
         & values=et, DIM1= DIM_node, DIM2= DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Water Surface Elevation At Last Timestep') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='meters') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='positive',values='up') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='standard_name',values='sea_surface_elevation') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='SSH_Mesh') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! TMEAN1
    VAR  => NC_MAKE_AVAR(name='tmean1',&
         & values=tmean1, DIM1= DIM_node, DIM2= DIM_siglay )

    ATT  => NC_MAKE_ATT(name='long_name',values='mean initial temperature') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='standard_name',values='sea_water_temperature') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='degrees_C') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='SigmaLayer_Mesh') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! SMEAN1
    VAR  => NC_MAKE_AVAR(name='smean1',&
         & values=smean1, DIM1= DIM_node, DIM2= DIM_siglay )

    ATT  => NC_MAKE_ATT(name='long_name',values='mean initial salinity') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='standard_name',values='sea_water_temperature') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='1e-3') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='SigmaLayer_Mesh') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    IF(OBC_ON)THEN
       ! OBC_GRID
       VAR  => NC_MAKE_AVAR(name='obc_nodes',&
            & values=I_OBC_N_OUTPUT, DIM1= DIM_nobc)
       
       ATT  => NC_MAKE_ATT(name='long_name',values='Open Boundary Node Number') 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='grid',values='obc_grid') 
       VAR  => ADD(VAR,ATT)
       
       NCF  => ADD(NCF,VAR)
       
       ! OBC_TYPE
       VAR  => NC_MAKE_AVAR(name='obc_type',&
            & values=type_obc, DIM1= DIM_nobc)
       
       ATT  => NC_MAKE_ATT(name='long_name',values='Open Boundary Type') 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='grid',values='obc_grid') 
       VAR  => ADD(VAR,ATT)
       
       NCF  => ADD(NCF,VAR)
    END IF

    IF(OBC_LONGSHORE_FLOW_ON)THEN
       ! long shore flow grid
       VAR  => NC_MAKE_AVAR(name='lsf_nodes',&
            & values=ibclsf_output, DIM1= DIM_nlsf)
       
       ATT  => NC_MAKE_ATT(name='long_name',values='Longshore Flow Node Number') 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='grid',values='lsf_grid') 
       VAR  => ADD(VAR,ATT)
       
       NCF  => ADD(NCF,VAR)
       
       ! 
       VAR  => NC_MAKE_AVAR(name='wdf',&
            & values=RBC_WDF, DIM1= DIM_nlsf)
       
       ATT  => NC_MAKE_ATT(name='long_name',values='Wind Driven Flow Adjustment Scaling') 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='valid_range',values='[0 1]') 
       VAR  => ADD(VAR,ATT)       

       ATT  => NC_MAKE_ATT(name='grid',values='lsf_grid') 
       VAR  => ADD(VAR,ATT)
       
       NCF  => ADD(NCF,VAR)

       VAR  => NC_MAKE_AVAR(name='geo',&
            & values=RBC_GEO, DIM1= DIM_nlsf)
       
       ATT  => NC_MAKE_ATT(name='long_name',values='Thermal Wind Flow Adjustment Scaling') 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='valid_range',values='[0 1]') 
       VAR  => ADD(VAR,ATT)       

       ATT  => NC_MAKE_ATT(name='grid',values='lsf_grid') 
       VAR  => ADD(VAR,ATT)
       
       NCF  => ADD(NCF,VAR)
    END IF


# if defined (EQUI_TIDE)
    VAR  => NC_MAKE_AVAR(name='el_eqi',&
         & values=el_eqi, DIM1= DIM_node, DIM2= DIM_time)
    
    ATT  => NC_MAKE_ATT(name='long_name',values='Equilibrium tide adjustment height') 
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='units',values='meters') 
    VAR  => ADD(VAR,ATT)
    
    NCF  => ADD(NCF,VAR)
# endif
    
# if defined (ATMO_TIDE)
    VAR  => NC_MAKE_AVAR(name='el_atmo',&
         & values=el_atmo, DIM1= DIM_node, DIM2= DIM_time)
    
    ATT  => NC_MAKE_ATT(name='long_name',values='Atmospheric tide adjustment height') 
    VAR  => ADD(VAR,ATT)
     
    ATT  => NC_MAKE_ATT(name='units',values='meters') 
    VAR  => ADD(VAR,ATT)
       
    NCF  => ADD(NCF,VAR)
# endif


    ! ------ New: Karsten Lettmann, 2016, June ---------------
    ! add the atmospheric pressure elevation
# if defined(AIR_PRESSURE)
    VAR  => NC_MAKE_AVAR(name='el_press',&
         & values=EL_AIR, DIM1= DIM_node, DIM2= DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='AIR-Pressure Induced Water Surface Elevation') 
    VAR  => ADD(VAR,ATT)
       
    ATT  => NC_MAKE_ATT(name='units',values='meters') 
    VAR  => ADD(VAR,ATT)
       
    NCF  => ADD(NCF,VAR)
# endif      
    ! -------- end new -------------------------------------


    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END RESTART_EXTRAS_FILE_OBJECT"
    
  END FUNCTION RESTART_EXTRAS_FILE_OBJECT
!=============================================================  
!=============================================================  
  FUNCTION WET_DRY_FILE_OBJECT() RESULT(NCF)
    USE MOD_WD
    IMPLICIT NONE

    INTEGER :: status
    LOGICAL, SAVE :: IOPROC_ALLOCATED = .FALSE.
    TYPE(NCFILE), POINTER :: NCF
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCATT),  POINTER :: ATT

    character(len=100)    :: timestamp, temp, netcdf_convention

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START WET_DRY_FILE_OBJECT"
    
    ! IO PROC MUST ALLOCATE SPACE FOR THE ARRAYS 
    ! THESE ARRAYS MUST HAVE THE ATTRIBUTE SAVE AND FOR CLARITY
    ! SHOULD HAVE THE SAME NAME AS THOSE USED ON THE OTHER PROCESSORS!
    IF(IOPROC .AND. .NOT. IOPROC_ALLOCATED) THEN

       IOPROC_ALLOCATED = .true.

       allocate(ISWETN(MGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:ISWETN")
       ISWETN = 0

       allocate(ISWETC(NGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:ISWETC")
       ISWETC = 0 

       allocate(ISWETNT(MGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:ISWETNT")
       ISWETN = 0

       allocate(ISWETCT(NGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:ISWETCT")
       ISWETCT = 0 

       allocate(ISWETCE(NGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:ISWETCE")
       ISWETCE = 0


    END IF


  ! ALLOCATE THE NEW FILE OBJECT
    NCF => NEW_FILE()


    ! WET NODES
    VAR  => NC_MAKE_AVAR(name='wet_nodes',&
         & values=ISWETN, DIM1= DIM_node, DIM2= DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Wet_Nodes') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates',values='time lat lon') 
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='mesh',values='fvcom_mesh') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='location',values='node') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! WET CELLS
    VAR  => NC_MAKE_AVAR(name='wet_cells',&
         & values=ISWETC, DIM1= DIM_nele, DIM2= DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Wet_Cells') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates',values='time latc lonc') 
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='mesh',values='fvcom_mesh') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='location',values='face') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! WET NODES AT LAST INT STEP
    VAR  => NC_MAKE_AVAR(name='wet_nodes_prev_int',&
         & values=ISWETNT, DIM1= DIM_node, DIM2= DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Wet_Nodes_At_Previous_Internal_Step') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates',values='time lat lon') 
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='mesh',values='fvcom_mesh') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='location',values='node') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! WET CELLS AT LAST EXT STEP
    VAR  => NC_MAKE_AVAR(name='wet_cells_prev_int',&
         & values=ISWETCT, DIM1= DIM_nele, DIM2= DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Wet_Cells_At_Previous_Internal_Step') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates',values='time latc lonc') 
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='mesh',values='fvcom_mesh') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='location',values='face') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! WET CELLS AT LAST EXT STEP
    VAR  => NC_MAKE_AVAR(name='wet_cells_prev_ext',&
         & values=ISWETCE, DIM1= DIM_nele, DIM2= DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Wet_Cells_At_Previous_External_Step') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)


    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END WET_DRY_FILE_OBJECT"
    
  END FUNCTION WET_DRY_FILE_OBJECT

!=============================================================  
!=============================================================  
  FUNCTION ICING_FILE_OBJECT() RESULT(NCF)
    IMPLICIT NONE

    INTEGER :: status
    LOGICAL, SAVE :: IOPROC_ALLOCATED = .FALSE.
    TYPE(NCFILE), POINTER :: NCF
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCATT),  POINTER :: ATT

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START ICING_FILE_OBJECT"
    
    ! IO PROC MUST ALLOCATE SPACE FOR THE ARRAYS 
    ! THESE ARRAYS MUST HAVE THE ATTRIBUTE SAVE AND FOR CLARITY
    ! SHOULD HAVE THE SAME NAME AS THOSE USED ON THE OTHER PROCESSORS!
    IF(IOPROC .AND. .NOT. IOPROC_ALLOCATED) THEN

       IOPROC_ALLOCATED = .true.

       allocate(ICING_0kts(MGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:ICING_0kts")
       ICING_0kts = 0.0_SP

       allocate(ICING_10kts(MGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:ICING_10kts")
       ICING_10kts = 0.0_SP

       allocate(icing_wndX(MGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:ICING_wndY")
       ICING_wndX = 0.0_SP

       allocate(icing_wndY(MGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:ICING_wndX")
       ICING_wndY = 0.0_SP

       allocate(icing_satmp(MGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:ICING_satmp")
       ICING_satmp = 0.0_SP

    END IF


  ! ALLOCATE THE NEW FILE OBJECT
    NCF => NEW_FILE()


    ! ICING_0KTS
    VAR  => NC_MAKE_AVAR(name='icing_0kts',&
         & values=icing_0kts, DIM1= DIM_node, DIM2= DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Icing Hazard@0knots') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='m C s^-1') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! ICING_10KTS
    VAR  => NC_MAKE_AVAR(name='icing_10kts',&
         & values=icing_10kts, DIM1= DIM_node, DIM2= DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Icing Hazard@10knots') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='m C s^-1') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! ICING_WNDX
    VAR  => NC_MAKE_AVAR(name='icing_wndx',&
         & values=icing_wndx, DIM1= DIM_node, DIM2= DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Icing Wind x-direction') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='m s^-1') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! ICING_WNDY
    VAR  => NC_MAKE_AVAR(name='icing_wndy',&
         & values=icing_wndy, DIM1= DIM_node, DIM2= DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Icing Wind y-direction') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='m s^-1') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! ICING_SATMP
    VAR  => NC_MAKE_AVAR(name='icing_satmp',&
         & values=icing_satmp, DIM1= DIM_node, DIM2= DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Icing Surface Air Temperature') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='degrees_C') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)


    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END ICING_FILE_OBJECT"
    
  END FUNCTION ICING_FILE_OBJECT
!=============================================================  
!=============================================================  
  FUNCTION GROUNDWATER_FILE_OBJECT() RESULT(NCF)
    IMPLICIT NONE

    INTEGER :: status
    LOGICAL, SAVE :: IOPROC_ALLOCATED = .FALSE.
    TYPE(NCFILE), POINTER :: NCF
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCATT),  POINTER :: ATT

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START: GROUNDWATER_FILE_OBJECT"
    
    ! IO PROC MUST ALLOCATE SPACE FOR THE ARRAYS 
    ! THESE ARRAYS MUST HAVE THE ATTRIBUTE SAVE AND FOR CLARITY
    ! SHOULD HAVE THE SAME NAME AS THOSE USED ON THE OTHER PROCESSORS!
    IF(IOPROC .AND. .NOT. IOPROC_ALLOCATED) THEN

       IOPROC_ALLOCATED = .true.

       allocate(BFWDIS(MGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:BFWDIS")
       BFWDIS = 0.0_SP

       IF(GROUNDWATER_TEMP_ON) THEN
          allocate(BFWTMP(MGL),stat=status)
          IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:BFWTMP")
          BFWTMP = 0.0_SP
       END IF

       IF(GROUNDWATER_SALT_ON) THEN
          allocate(BFWSLT(MGL),stat=status)
          IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:BFWSLT")
          BFWSLT = 0.0_SP
       END IF


    END IF


  ! ALLOCATE THE NEW FILE OBJECT
    NCF => NEW_FILE()


    ! GROUNDWATER VOLUME FLUX
    VAR  => NC_MAKE_AVAR(name='groundwater_flux',&
         & values=bfwdis, DIM1= DIM_node, DIM2= DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='groundwater volume flux') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='m3 s-1') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! GROUNDWATER INFLOW TEMPERATURE
    IF(GROUNDWATER_TEMP_ON) THEN
       VAR  => NC_MAKE_AVAR(name='groundwater_temp',&
            & values=bfwdis, DIM1= DIM_node, DIM2= DIM_time)
       
       ATT  => NC_MAKE_ATT(name='long_name',values='groundwater inflow temperature') 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='units',values='degrees_C') 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='type',values='data') 
       VAR  => ADD(VAR,ATT)
       
       NCF  => ADD(NCF,VAR)
    END IF

    ! GROUNDWATER INFLOW SALINITY
    IF(GROUNDWATER_SALT_ON) THEN
       VAR  => NC_MAKE_AVAR(name='groundwater_salt',&
            & values=bfwdis, DIM1= DIM_node, DIM2= DIM_time)
       
       ATT  => NC_MAKE_ATT(name='long_name',values='groundwater inflow salinity') 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='units',values='1e-3') 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='type',values='data') 
       VAR  => ADD(VAR,ATT)
       
       NCF  => ADD(NCF,VAR)
    END IF


    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END GROUND_FILE_OBJECT"
    
  END FUNCTION GROUNDWATER_FILE_OBJECT
!=============================================================  

# if defined (NH)
!=============================================================
  FUNCTION QP_FILE_OBJECT() RESULT(NCF)
    IMPLICIT NONE

    INTEGER :: status
    LOGICAL, SAVE :: IOPROC_ALLOCATED = .FALSE.
    TYPE(NCFILE), POINTER :: NCF
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCATT),  POINTER :: ATT

    character(len=100)    :: timestamp, temp, netcdf_convention

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START QP_FILE_OBJECT"

    ! IO PROC MUST ALLOCATE SPACE FOR THE ARRAYS
    ! THESE ARRAYS MUST HAVE THE ATTRIBUTE SAVE AND FOR CLARITY
    ! SHOULD HAVE THE SAME NAME AS THOSE USED ON THE OTHER PROCESSORS!
    IF(IOPROC .AND. .NOT. IOPROC_ALLOCATED) THEN

       IOPROC_ALLOCATED = .true.

       allocate(QP(MGL,KBM1),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:QP")
       QP = 0.0_sp

    END IF

  ! ALLOCATE THE NEW FILE OBJECT
    NCF => NEW_FILE()

    ! QP
    VAR  => NC_MAKE_AVAR(name='qp',&
         & values=QP, DIM1= DIM_node, DIM2= DIM_siglay, DIM3 = DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Nonhydrostatic Pertubation Pressure')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='N M-2')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data')
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END QP_FILE_OBJECT"

  END FUNCTION QP_FILE_OBJECT
!=============================================================
  FUNCTION RHS_FILE_OBJECT() RESULT(NCF)
    IMPLICIT NONE

    INTEGER :: status
    LOGICAL, SAVE :: IOPROC_ALLOCATED = .FALSE.
    TYPE(NCFILE), POINTER :: NCF
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCATT),  POINTER :: ATT

    character(len=100)    :: timestamp, temp, netcdf_convention

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START RHS_FILE_OBJECT"

    ! IO PROC MUST ALLOCATE SPACE FOR THE ARRAYS
    ! THESE ARRAYS MUST HAVE THE ATTRIBUTE SAVE AND FOR CLARITY
    ! SHOULD HAVE THE SAME NAME AS THOSE USED ON THE OTHER PROCESSORS!
    IF(IOPROC .AND. .NOT. IOPROC_ALLOCATED) THEN

       IOPROC_ALLOCATED = .true.

       allocate(RHS(MGL,KBM1),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:RHS")
       RHS = 0.0_sp

    END IF

  ! ALLOCATE THE NEW FILE OBJECT
    NCF => NEW_FILE()

    ! RHS
    VAR  => NC_MAKE_AVAR(name='rhs',&
         & values=RHS, DIM1= DIM_node, DIM2= DIM_siglay, DIM3 = DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Right Hand Side Term')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='Kg M-3 S-2')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data')
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END RHS_FILE_OBJECT"

  END FUNCTION RHS_FILE_OBJECT
!=============================================================
  FUNCTION NH_RST_FILE_OBJECT() RESULT(NCF)
    IMPLICIT NONE

    INTEGER :: status
    LOGICAL, SAVE :: IOPROC_ALLOCATED = .FALSE.
    TYPE(NCFILE), POINTER :: NCF
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCATT),  POINTER :: ATT

    character(len=100)    :: timestamp, temp, netcdf_convention

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START NH_RST_FILE_OBJECT"
    
    ! IO PROC MUST ALLOCATE SPACE FOR THE ARRAYS 
    ! THESE ARRAYS MUST HAVE THE ATTRIBUTE SAVE AND FOR CLARITY
    ! SHOULD HAVE THE SAME NAME AS THOSE USED ON THE OTHER PROCESSORS!
    IF(IOPROC .AND. .NOT. IOPROC_ALLOCATED) THEN

       IOPROC_ALLOCATED = .true.

       allocate(W4ZT(MGL,KB),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:W4ZT")
       W4ZT = 0.0_SP

       allocate(NHQDRX(NGL,KBM1),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:NHQDRX")
       NHQDRX = 0.0_SP

       allocate(NHQDRY(NGL,KBM1),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:NHQDRY")
       NHQDRY = 0.0_SP

       allocate(NHQDRZ(MGL,KBM1),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:NHQDRZ")
       NHQDRZ = 0.0_SP

       allocate(NHQ2DX(NGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:NHQ2DX")
       NHQ2DX = 0.0_SP

       allocate(NHQ2DY(NGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:NHQ2DY")
       NHQ2DY = 0.0_SP

    END IF

  ! ALLOCATE THE NEW FILE OBJECT
    NCF => NEW_FILE()

    ! W4ZT
    VAR  => NC_MAKE_AVAR(name='w4zt',&
         & values=W4ZT, DIM1= DIM_node, DIM2= DIM_siglev, DIM3 = DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Temp Vertical Velocity in Z')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='meters s-1')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data')
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! NHQDRX
    VAR  => NC_MAKE_AVAR(name='nhqdrx',&
         & values=NHQDRX, DIM1= DIM_nele, DIM2= DIM_siglay, DIM3 = DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Nonhydrostatic Pressure Gradient XCor')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='N M-3')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data')
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! NHQDRY
    VAR  => NC_MAKE_AVAR(name='nhqdry',&
         & values=NHQDRY, DIM1= DIM_nele, DIM2= DIM_siglay, DIM3 = DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Nonhydrostatic Pressure Gradient YCor')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='N M-3')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data')
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! NHQDRZ
    VAR  => NC_MAKE_AVAR(name='nhqdrz',&
         & values=NHQDRZ, DIM1= DIM_node, DIM2= DIM_siglay, DIM3 = DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Nonhydrostatic Pressure Gradient ZCor')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='N M-3')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data')
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! NHQ2DX
    VAR  => NC_MAKE_AVAR(name='nhq2dx',&
         & values=NHQ2DX, DIM1= DIM_nele, DIM2= DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Nonhydrostatic Pressure Gradient XCor 2D')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='N M-3')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data')
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! NHQ2DY
    VAR  => NC_MAKE_AVAR(name='nhq2dy',&
         & values=NHQ2DY, DIM1= DIM_nele, DIM2= DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Nonhydrostatic Pressure Gradient YCor 2D')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='N M-3')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data')
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)


    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END NH_RST_FILE_OBJECT"
    
  END FUNCTION NH_RST_FILE_OBJECT
# endif

!!! ggao/0104/2008  !! restart file for ice model
!!---------------------------------------------------------------------------
# if defined (ICE)
!=============================================================  
  FUNCTION ICE_DATA_STATE_FILE_OBJECT() RESULT(NCF)
    
    USE mod_ICE
    USE mod_ICE2D

    IMPLICIT NONE

    INTEGER :: status
    LOGICAL, SAVE :: IOPROC_ALLOCATED = .FALSE.
    TYPE(NCFILE), POINTER :: NCF
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCATT),  POINTER :: ATT

    REAL(SP),SAVE, POINTER :: AICE_P(:)
    REAL(SP),SAVE, POINTER :: VICE_P(:)

#   if defined (HEATING_CALCULATED) || (HEATING_SOLAR)
!   afm for heat flux output 20160513 & EJA 20160921 ---
    REAL(SP),SAVE, POINTER :: FSH_P(:), FLH_P(:),  &
                              FLWUP_P(:), FLW_P(:), TSFC_P(:)
!-----------------------------------------
#   endif

    character(len=100)    :: timestamp, temp, netcdf_convention

    integer :: NI,K,I,J

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START ICE_STATE_FILE_OBJECT"
    
    ! IO PROC MUST ALLOCATE SPACE FOR THE ARRAYS 
       
  ! ALLOCATE THE NEW FILE OBJECT
    NCF => NEW_FILE()


   ! ICE AREA 

    AICE_P => aice(1:m,1)

!    VAR  => NC_MAKE_AVAR(name='aice',&
!         & values=aice, DIM1= DIM_node, DIM2= Dim_grid, DIM3 = DIM_time)

    VAR  => NC_MAKE_PVAR(name='aice',&
         & values=aice_P, DIM1= DIM_node, DIM2 = DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='ice area') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='standard_name',values='ice_area_category') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar) 
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)    

    NCF  => ADD(NCF,VAR)
    
   ! ICE VOLUME
    VICE_P => vice(1:m,1)

!    VAR  => NC_MAKE_AVAR(name='vice',&
!         & values=vice, DIM1= DIM_node, DIM2= Dim_grid, DIM3 = DIM_time)

    VAR  => NC_MAKE_PVAR(name='vice',&
         & values=vice_P, DIM1= DIM_node, DIM2= DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='ice volume per unit grid area') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='standard_name',values='ice_volume_category') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='m') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar) 
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)    

    NCF  => ADD(NCF,VAR)


#   if defined (HEATING_CALCULATED) || (HEATING_SOLAR)
!   afm 20160513 & EJA 20160921 - For additional outputs of ice-related variables 
!----------------------
   ! ICE SURFACE TEMPERATURE
    TSFC_P => TSFC(1:m,1)

    VAR  => NC_MAKE_PVAR(name='tsfc',&
         & values=TSFC_P, DIM1= DIM_node, DIM2= DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='ice surface temperature')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='standard_name',values='ice_surface_temperature')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='deg C')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar)
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data')
    VAR  => ADD(VAR,ATT)

!    ATT  => NC_MAKE_ATT(name='missing_value',values=0.0_SP)
!    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)
    
    !Edit by EJA - add heat fluxes via ice model if "Surface Heating" NC is not turned on
    IF (.not. NC_SURFACE_HEAT) THEN
   ! SENSIBLE HEAT FLUX AT WATER SURFACE
    FSH_P => FSH(1:m,1)

    VAR  => NC_MAKE_PVAR(name='sensible_heat_flux',&
         & values=FSH_P, DIM1= DIM_node, DIM2= DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Sensible Heat Flux')
    VAR  => ADD(VAR,ATT)

    !ATT  => NC_MAKE_ATT(name='standard_name',values='fsh')
    !VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='W m-2')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar)
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='mesh',values='fvcom_mesh')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='location',values='node')
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)


   ! LATENT HEAT FLUX AT WATER SURFACE
    FLH_P => FLH(1:m,1)

    VAR  => NC_MAKE_PVAR(name='latent_heat_flux',&
         & values=FLH_P, DIM1= DIM_node, DIM2= DIM_time)
    
    ATT  => NC_MAKE_ATT(name='long_name',values='Latent Heat Flux')
    VAR  => ADD(VAR,ATT)

    !ATT  => NC_MAKE_ATT(name='standard_name',values='flh')
    !VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='W m-2')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar)
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='mesh',values='fvcom_mesh')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='location',values='node')
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)


   ! LONGWAVE RADIATION AT WATER SURFACE
    FLW_P => FLWTOT(1:m,1)

    VAR  => NC_MAKE_PVAR(name='long_wave',&
         & values=FLW_P, DIM1= DIM_node, DIM2= DIM_time)
    
    ATT  => NC_MAKE_ATT(name='long_name',values='Long Wave Radiation')
    VAR  => ADD(VAR,ATT)

    !ATT  => NC_MAKE_ATT(name='standard_name',values='flwp')
    !VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='W m-2')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar)
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='mesh',values='fvcom_mesh')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='location',values='node')
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)


   ! INCOMING LONGWAVE RADIATION AT WATER SURFACE
   ! FLW_P => FLW(1:m,1)

   ! VAR  => NC_MAKE_PVAR(name='flw',&
   !      & values=FLW_P, DIM1= DIM_node, DIM2= DIM_time)
    
   ! ATT  => NC_MAKE_ATT(name='long_name',values='incoming longwave radiation at water surface')
   ! VAR  => ADD(VAR,ATT)

   ! ATT  => NC_MAKE_ATT(name='standard_name',values='flw')
   ! VAR  => ADD(VAR,ATT)

   ! ATT  => NC_MAKE_ATT(name='units',values='W/m2')
   ! VAR  => ADD(VAR,ATT)

   ! ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid')
   ! VAR  => ADD(VAR,ATT)

   ! ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar)
   ! VAR  => ADD(VAR,ATT)

   ! ATT  => NC_MAKE_ATT(name='type',values='data')
   ! VAR  => ADD(VAR,ATT)

   ! NCF  => ADD(NCF,VAR)
    END IF

!----------------------
#   endif

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END ICE_STATE0_FILE_OBJECT"

  END FUNCTION ICE_DATA_STATE_FILE_OBJECT
!=============================================================


!=============================================================  
  FUNCTION ICE_RESTART_STATE_FILE_OBJECT() RESULT(NCF)
    USE mod_ICE
    USE mod_ICE2D
    IMPLICIT NONE

    INTEGER :: status
    LOGICAL, SAVE :: IOPROC_ALLOCATED = .FALSE.
    TYPE(NCFILE), POINTER :: NCF
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCATT),  POINTER :: ATT

    character(len=100)    :: timestamp, temp, netcdf_convention

    REAL(SP),SAVE, POINTER :: AICEN_P(:,:)
    REAL(SP),SAVE, POINTER :: VICEN_P(:,:)
    REAL(SP),SAVE, POINTER :: VSNON_P(:,:)
    REAL(SP),SAVE, POINTER :: Esnon_P(:,:)
    REAL(SP),SAVE, POINTER :: Tsfcn_P(:,:)
    REAL(SP),SAVE, POINTER :: Eicen_P(:,:)

    integer :: NI,K,I,J

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START ICE_STATE_FILE_OBJECT"
    
    
  ! ALLOCATE THE NEW FILE OBJECT
    NCF => NEW_FILE()

   ! ICE AREA 
   
    AICEN_P => aicen(1:imt_local,1,1:ncat)
 
    VAR  => NC_MAKE_PVAR(name='AICEN',&
         & values=AICEN_P, DIM1= DIM_node, DIM2= DIM_NCAT, DIM3 = DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='ice area') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='standard_name',values='ice_area_category') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar) 
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)    

    NCF  => ADD(NCF,VAR)
    
   ! ICE VOLUME
      VICEN_P => vicen(1:imt_local,1,1:ncat)
    
    VAR  => NC_MAKE_PVAR(name='VICEN',&
         & values=Vicen_P, DIM1= DIM_node, DIM2= DIM_NCAT, DIM3 = DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='ice volume per unit grid area') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='standard_name',values='ice_volume_category') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='m') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar) 
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)    

    NCF  => ADD(NCF,VAR)
    
   ! SNOW VOLUME
       Vsnon_P =>  vsnon(1:imt_local,1,1:ncat)
    
    VAR  => NC_MAKE_PVAR(name='VSNON',&
         & values=Vsnon_P, DIM1= DIM_node, DIM2= DIM_NCAT, DIM3 = DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='snow volume per unit grid area') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='standard_name',values='snow_volume_category') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='m') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar) 
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)    

    NCF  => ADD(NCF,VAR)
     

   ! ice/snow surface temperature
   
    Tsfcn_P =>  Tsfcn(1:imt_local,1,1:ncat)
    
    
    VAR  => NC_MAKE_PVAR(name='TSFCN',&
         & values=Tsfcn_P, DIM1= DIM_node, DIM2= DIM_NCAT, DIM3 = DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='ice/snow surface temperature') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='standard_name',values='ice_surface_temperature') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='deg C') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar) 
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)    

    NCF  => ADD(NCF,VAR)
     
   !  snow energy 
       Esnon_P => esnon(1:imt_local,1,1:ncat)
    
    VAR  => NC_MAKE_PVAR(name='ESNON',&
         & values=Esnon_P, DIM1= DIM_node, DIM2= DIM_NCAT, DIM3 = DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='snow energy ') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='standard_name',values='snow_energy ') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='J/m$^2$') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar) 
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)    

    NCF  => ADD(NCF,VAR)
    
!!-------------------------------------------------------------------    
   !  ice energy category and layer

    Eicen_P => eicen(1:imt_local,1,1:ntilay)
   
    VAR  => NC_MAKE_PVAR(name='EICEN',&
         & values=Eicen_P, DIM1= DIM_node, DIM2= DIM_ntilay, DIM3 = DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='ice energy ') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='standard_name',values='ice_energy ') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='J/m$^2$') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar) 
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)    

    NCF  => ADD(NCF,VAR)

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END ICE_STATE_FILE_OBJECT"

  END FUNCTION ICE_RESTART_STATE_FILE_OBJECT
!=============================================================

!=============================================================  
  FUNCTION ICE_VEL_FILE_OBJECT() RESULT(NCF)
    USE mod_ICE
    USE mod_ICE2D
    IMPLICIT NONE

    INTEGER :: status
    LOGICAL, SAVE :: IOPROC_ALLOCATED = .FALSE.
    TYPE(NCFILE), POINTER :: NCF
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCATT),  POINTER :: ATT

    character(len=100)    :: timestamp, temp, netcdf_convention

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START ICE_VEL_FILE_OBJECT"
    
    ! IO PROC MUST ALLOCATE SPACE FOR THE ARRAYS 
    ! THESE ARRAYS MUST HAVE THE ATTRIBUTE SAVE AND FOR CLARITY
    ! SHOULD HAVE THE SAME NAME AS THOSE USED ON THE OTHER PROCESSORS!
    IF(IOPROC .AND. .NOT. IOPROC_ALLOCATED) THEN

       IOPROC_ALLOCATED = .true.

       allocate(UICE2(NGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:UICE")
       UICE2 = 0.0_SP

       allocate(VICE2(NGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:VICE")
       VICE2 = 0.0_SP

    END IF


  ! ALLOCATE THE NEW FILE OBJECT
    NCF => NEW_FILE()


    ! UICE2
    VAR  => NC_MAKE_AVAR(name='uuice',&
         & values=Uice2, DIM1= DIM_nele, DIM2= DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Ice x-velocity') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='meters s-1') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)


    ! VICE2
    VAR  => NC_MAKE_AVAR(name='vvice',&
         & values=VICE2, DIM1= DIM_nele, DIM2= DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Ice y-velocity') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='meters s-1') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)


    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END ICE_VEL_FILE_OBJECT"

  END FUNCTION ICE_VEL_FILE_OBJECT

!============================================================= 
  FUNCTION ICE_EXTRA_FILE_OBJECT() RESULT(NCF)
    USE mod_ICE
    USE mod_ICE2D
    IMPLICIT NONE

    INTEGER :: status
    LOGICAL, SAVE :: IOPROC_ALLOCATED = .FALSE.
    TYPE(NCFILE), POINTER :: NCF
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCATT),  POINTER :: ATT

    character(len=100)    :: timestamp, temp, netcdf_convention

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START ICE_EXTRA_FILE_OBJECT"

    ! IO PROC MUST ALLOCATE SPACE FOR THE ARRAYS
    ! THESE ARRAYS MUST HAVE THE ATTRIBUTE SAVE AND FOR CLARITY
    ! SHOULD HAVE THE SAME NAME AS THOSE USED ON THE OTHER PROCESSORS!
    IF(IOPROC .AND. .NOT. IOPROC_ALLOCATED) THEN

       IOPROC_ALLOCATED = .true.

       allocate(SIG1(MGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:SIG1")
       SIG1 = 0.0_SP

       allocate(SIG2(MGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:SIG2")
       SIG2 = 0.0_SP

       allocate(SIG12(MGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:SIG12")
       SIG12 = 0.0_SP

    END IF
    
  ! ALLOCATE THE NEW FILE OBJECT
    NCF => NEW_FILE()

      !-----------------------------------------------------------------
      ! fresh water, salt, and heat flux
      !-----------------------------------------------------------------

   ! fresh water
    
    VAR  => NC_MAKE_AVAR(name='fresh',&
         & values=fresh, DIM1= DIM_node,Dim2 =DIM_GRID, DIM3= DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='flux of water ice to ocean') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='standard_name',values='fresh_water') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='kg/m^2/s') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar) 
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)    

    NCF  => ADD(NCF,VAR)

   ! salt
    
    VAR  => NC_MAKE_AVAR(name='fsalt',&
         & values=fsalt, DIM1= DIM_node,Dim2 =DIM_GRID, DIM3= DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='flux of salt ice to ocean') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='standard_name',values='flux_salt') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='kg/m^2/s') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar) 
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)    

    NCF  => ADD(NCF,VAR)


   ! net heat
    
    VAR  => NC_MAKE_AVAR(name='fhnet',&
         & values=fhnet, DIM1= DIM_node,Dim2 =DIM_GRID, DIM3= DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='heat flux ice to ocean') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='standard_name',values='flux_heat') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='W/m^2') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar) 
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)    

    NCF  => ADD(NCF,VAR)
    
      !-----------------------------------------------------------------
      ! ice strength
      !-----------------------------------------------------------------
   ! ice strength
    
    VAR  => NC_MAKE_AVAR(name='strength',&
         & values=strength, DIM1= DIM_node,Dim2 =DIM_GRID, DIM3= DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='ice strength') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='standard_name',values='ice_strength') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='N/m') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar) 
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)    

    NCF  => ADD(NCF,VAR)

      !-----------------------------------------------------------------
      ! ice mask for dynamics
      !-----------------------------------------------------------------
    ! ice cell mask
    VAR  => NC_MAKE_AVAR(name='isicec',&
         & values=isicec, DIM1= DIM_nele, DIM2= DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='ice cell mask') 
    VAR  => ADD(VAR,ATT)

!    ATT  => NC_MAKE_ATT(name='units',values='') 
!    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

      !-----------------------------------------------------------------
      ! Ice stress 1
      !-----------------------------------------------------------------

    VAR  => NC_MAKE_AVAR(name='sig1',&
         & values=sig1, DIM1= DIM_node,DIM2= DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='ice stress')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='standard_name',values='ice stress divergence')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='N/m')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar)
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data')
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

      !-----------------------------------------------------------------
      ! Ice stress 2
      !-----------------------------------------------------------------

    VAR  => NC_MAKE_AVAR(name='sig2',&
         & values=sig2, DIM1= DIM_node,DIM2= DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='ice stress2')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='standard_name',values='ice stress tension')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='N/m')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar)
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data')
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)


      !-----------------------------------------------------------------
      ! Ice stress 12
      !-----------------------------------------------------------------

    VAR  => NC_MAKE_AVAR(name='sig12',&
         & values=sig12, DIM1= DIM_node,DIM2= DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='ice stress12')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='standard_name',values='ice shear stress')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='N/m')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar)
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data')
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)


    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END ICE_EXTRA_FILE_OBJECT"

  END FUNCTION ICE_EXTRA_FILE_OBJECT
!=============================================================
# endif
!!---------------------------------------------------------------------------
!!---------------------------------------------------------------------------
!! ggao/0104/2008  !! restart file for ice model
!=============================================================  
# if defined (BioGen)
  FUNCTION BIO_FILE_OBJECT() RESULT(NCF)
    USE MOD_BIO_3D
   IMPLICIT NONE

   INTEGER :: status
   LOGICAL, SAVE :: IOPROC_ALLOCATED = .FALSE.
   TYPE(NCFILE), POINTER :: NCF
   TYPE(NCVAR),  POINTER :: VAR
   TYPE(NCATT),  POINTER :: ATT
   INTEGER :: II
   REAL(SP), POINTER :: BIO_ALL_PNT(:,:)

   IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START BIO_FILE_OBJECT"
   
   IF(IOPROC .AND. .NOT. IOPROC_ALLOCATED) THEN

     IOPROC_ALLOCATED = .true.

     ALLOCATE(BIO_ALL(MGL,KB,NTT),STAT=STATUS)
     IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:BIO_ALL")
     BIO_ALL = 0.0_SP

   END IF

  ! ALLOCATE THE NEW FILE OBJECT
   NCF => NEW_FILE()

   DO II=1,NTT
!JQI     BIO_ALL_PNT(:,:) => BIO_ALL(:,:,II)
     BIO_ALL_PNT => BIO_ALL(1:M,1:KB,II)
    
   ! IO PROC MUST ALLOCATE SPACE FOR THE ARRAYS 
   ! THESE ARRAYS MUST HAVE THE ATTRIBUTE SAVE AND FOR CLARITY
   ! SHOULD HAVE THE SAME NAME AS THOSE USED ON THE OTHER PROCESSORS!


!  ! ALLOCATE THE NEW FILE OBJECT
!   NCF => NEW_FILE()

  ! BIO
     VAR  => NC_MAKE_PVAR(name=TRIM(BIO_NAME(II,1)),       &
           & values=BIO_ALL_PNT, DIM1= DIM_node, DIM2= DIM_siglay, DIM3 = DIM_time)

     ATT  => NC_MAKE_ATT(name='long_name',values=TRIM(BIO_NAME(II,3))) 
     VAR  => ADD(VAR,ATT)

!    ATT  => NC_MAKE_ATT(name='standard_name',values='sea_water_temperature') 
!    VAR  => ADD(VAR,ATT)

     ATT  => NC_MAKE_ATT(name='units',values=TRIM(BIO_NAME(II,2)) )
     VAR  => ADD(VAR,ATT)

     ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
     VAR  => ADD(VAR,ATT)

     ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar) 
     VAR  => ADD(VAR,ATT)
    
     ATT  => NC_MAKE_ATT(name='type',values='data') 
     VAR  => ADD(VAR,ATT)    

     NCF  => ADD(NCF,VAR)
  
     NULLIFY(BIO_ALL_PNT)

   END DO

   IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END BIO_FILE_OBJECT"

  END FUNCTION BIO_FILE_OBJECT
# endif
!=============================================================
# if defined (WATER_QUALITY)
  FUNCTION WQM_FILE_OBJECT() RESULT(NCF)
    USE CONTROL, ONLY : NB, WQM_NAME
    USE MOD_WQM, ONLY : WQM  
    IMPLICIT NONE

   INTEGER :: status
   LOGICAL, SAVE :: IOPROC_ALLOCATED = .FALSE.
   TYPE(NCFILE), POINTER :: NCF
   TYPE(NCVAR),  POINTER :: VAR
   TYPE(NCATT),  POINTER :: ATT
   INTEGER :: II
   REAL(SP),ALLOCATABLE :: WQM_TMP(:,:)

   IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START WQM_FILE_OBJECT"
   
   DO II = 1,NB
     ALLOCATE(WQM_TMP(0:MT,KB))
     WQM_TMP(:,:) = WQM(:,:,II)
    
   ! IO PROC MUST ALLOCATE SPACE FOR THE ARRAYS 
   ! THESE ARRAYS MUST HAVE THE ATTRIBUTE SAVE AND FOR CLARITY
   ! SHOULD HAVE THE SAME NAME AS THOSE USED ON THE OTHER PROCESSORS!
   IF(IOPROC .AND. .NOT. IOPROC_ALLOCATED) THEN

     IOPROC_ALLOCATED = .true.

     ALLOCATE(WQM_TMP(MGL,KB),STAT=STATUS)
     IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:WQM")
     WQM_TMP = 0.0_SP

   END IF


  ! ALLOCATE THE NEW FILE OBJECT
   NCF => NEW_FILE()

  ! WQM
    
     VAR  => NC_MAKE_AVAR(name=TRIM(WQM_NAME(II,1)),&
           & values=WQM_TMP, DIM1= DIM_node, DIM2= DIM_siglay, DIM3 = DIM_time)

     ATT  => NC_MAKE_ATT(name='long_name',values=TRIM(WQM_NAME(II,3))) 
     VAR  => ADD(VAR,ATT)

!    ATT  => NC_MAKE_ATT(name='standard_name',values='sea_water_temperature') 
!    VAR  => ADD(VAR,ATT)

     ATT  => NC_MAKE_ATT(name='units',values=TRIM(WQM_NAME(II,2))) 
     VAR  => ADD(VAR,ATT)

     ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
     VAR  => ADD(VAR,ATT)

     ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar) 
     VAR  => ADD(VAR,ATT)
    
     ATT  => NC_MAKE_ATT(name='type',values='data') 
     VAR  => ADD(VAR,ATT)    

     NCF  => ADD(NCF,VAR)
  
     DEALLOCATE(WQM_TMP)
   END DO

   IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END WQM_FILE_OBJECT"

  END FUNCTION WQM_FILE_OBJECT
# endif
!=============================================================  
  FUNCTION TIME_FILE_OBJECT() RESULT(NCF)
   IMPLICIT NONE

   INTEGER :: status
   TYPE(NCFILE), POINTER :: NCF
   TYPE(NCVAR),  POINTER :: VAR
   TYPE(NCATT),  POINTER :: ATT
   INTEGER :: II

   IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START TIME_OBJECT"

  ! ALLOCATE THE NEW FILE OBJECT
   NCF => NEW_FILE()

   ! IINT
   VAR  => IINT_OBJECT(DIM=Dim_time)
   
   NCF  => ADD(NCF,VAR)
   
   ! time
   VAR => FLOAT_TIME_OBJECT &
        &(USE_MJD=use_real_world_time, &
        & DIM=DIM_TIME)
   
   NCF  => ADD(NCF,VAR)
   
   
   ! Itime
   VAR  => ITIME_OBJECT &
        &(Use_MJD=use_real_world_time, &
        & DIM=DIM_TIME)
   
   NCF  => ADD(NCF,VAR)
   
   ! Itime2
   VAR => ITIME2_OBJECT &
        &(Use_MJD=use_real_world_time, &
        & DIM=DIM_TIME)
   
   NCF => ADD(NCF,VAR)
   
   IF (use_real_world_time) THEN
      
      VAR => DATETIME_OBJECT &
           &(DIMSTR=DIM_DateStrLen,&
           & DIMTIME=DIM_TIME,&
           TIMEZONE=TIMEZONE)
      
      NCF  => ADD(NCF,VAR)
   END IF

   IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END TIME_OBJECT"

 END FUNCTION TIME_FILE_OBJECT
 !=============================================================  

# if defined (WAVE_CURRENT_INTERACTION)

  FUNCTION WAVE_PARA_FILE_OBJECT() RESULT(NCF)

   USE SWCOMM3
   USE VARS_WAVE
   USE MOD_WAVE_CURRENT_INTERACTION

   IMPLICIT NONE

    INTEGER :: status
    LOGICAL, SAVE :: IOPROC_ALLOCATED = .FALSE.
    TYPE(NCFILE), POINTER :: NCF
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCATT),  POINTER :: ATT

    character(len=100)    :: timestamp, temp, netcdf_convention

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START WAVE_PARA_FILE_OBJECT"
    
    ! IO PROC MUST ALLOCATE SPACE FOR THE ARRAYS 
    ! THESE ARRAYS MUST HAVE THE ATTRIBUTE SAVE AND FOR CLARITY
    ! SHOULD HAVE THE SAME NAME AS THOSE USED ON THE OTHER PROCESSORS!
    IF(IOPROC .AND. .NOT. IOPROC_ALLOCATED) THEN

       IOPROC_ALLOCATED = .true.

       allocate(HSC1(MGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:HSC1")
       HSC1 = 0.0_sp

       allocate(DIRDEG1(MGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:DIRDEG1")
       DIRDEG1 = 0.0_sp

       allocate(TPEAK(MGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:TPEAK")
       TPEAK = 0.0_sp

       allocate(WLEN(MGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:WLEN")
       WLEN = 0.0_sp

    END IF


  ! ALLOCATE THE NEW FILE OBJECT
    NCF => NEW_FILE()

    ! Wave height
    VAR  => NC_MAKE_AVAR(name='hs',&
         & values=HSC1, DIM1= DIM_node, DIM2 = DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Significant Wave Height') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='meters') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    !  Wave Direction
    VAR  => NC_MAKE_AVAR(name='wdir',&
         & values=DIRDEG1, DIM1= DIM_node, DIM2 = DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Wave Direction') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='degree') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)


    !  Relative Peak Period
    VAR  => NC_MAKE_AVAR(name='tpeak',&
         & values=tpeak, DIM1= DIM_node, DIM2 = DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values=' Relative Peak Period') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='s') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    !  Wavelength
    VAR  => NC_MAKE_AVAR(name='wlen',&
         & values=wlen, DIM1= DIM_node, DIM2 = DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Wavelength') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='meter') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END WAVE_PARA_FILE_OBJECT"


  END  FUNCTION WAVE_PARA_FILE_OBJECT 

   FUNCTION WAVE_PARA1_FILE_OBJECT() RESULT(NCF)

   USE SWCOMM3
   USE VARS_WAVE
   USE MOD_WAVE_CURRENT_INTERACTION

   IMPLICIT NONE

    INTEGER :: status
    LOGICAL, SAVE :: IOPROC_ALLOCATED = .FALSE.
    TYPE(NCFILE), POINTER :: NCF
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCATT),  POINTER :: ATT

    character(len=100)    :: timestamp, temp, netcdf_convention

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START WAVE_PARA1_FILE_OBJECT"
    
    ! IO PROC MUST ALLOCATE SPACE FOR THE ARRAYS 
    ! THESE ARRAYS MUST HAVE THE ATTRIBUTE SAVE AND FOR CLARITY
    ! SHOULD HAVE THE SAME NAME AS THOSE USED ON THE OTHER PROCESSORS!
    IF(IOPROC .AND. .NOT. IOPROC_ALLOCATED) THEN

       IOPROC_ALLOCATED = .true.

       allocate(QB1(MGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:QB1")
       QB1 = 0.0_sp

       allocate(UB_SWAN(MGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:UB_SWAN")
       UB_SWAN = 0.0_sp

       allocate(PWAVE_BOT(MGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:PWAVE_BOT")
       PWAVE_BOT = 0.0_sp

       allocate(DIRBOT(MGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:DIRBOT")
       DIRBOT = 0.0_sp

    END IF


  ! ALLOCATE THE NEW FILE OBJECT
    NCF => NEW_FILE()

    !  Fraction of Breaking Waves
    VAR  => NC_MAKE_AVAR(name='qb',&
         & values=qb1, DIM1= DIM_node, DIM2 = DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Fraction of Breaking Waves') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='meter') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

!   bottom orbital velocity 
    VAR  => NC_MAKE_AVAR(name='ubot',&
         & values=ub_swan, DIM1= DIM_node, DIM2 = DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='bottom wave orbital velocity')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='m/s')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data')
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

!   bottom orbital period 
    VAR  => NC_MAKE_AVAR(name='tmbot',&
         & values=pwave_bot, DIM1= DIM_node, DIM2 = DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='bottom wave orbital period')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='s')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data')
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

!   bottom wave direction 
    VAR  => NC_MAKE_AVAR(name='dirbot',&
         & values=DIRBOT, DIM1= DIM_node, DIM2 = DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='bottom wave direction')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='degree')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data')
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END WAVE_PARA1_FILE_OBJECT"


  END  FUNCTION WAVE_PARA1_FILE_OBJECT 

 FUNCTION WAVE_STRESS_FILE_OBJECT() RESULT(NCF)

   USE SWCOMM3
   USE VARS_WAVE
   USE MOD_WAVE_CURRENT_INTERACTION

   IMPLICIT NONE

    INTEGER :: status
    LOGICAL, SAVE :: IOPROC_ALLOCATED = .FALSE.
    TYPE(NCFILE), POINTER :: NCF
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCATT),  POINTER :: ATT

    character(len=100)    :: timestamp, temp, netcdf_convention

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START WAVE_STRESS_FILE_OBJECT"
    
    ! IO PROC MUST ALLOCATE SPACE FOR THE ARRAYS 
    ! THESE ARRAYS MUST HAVE THE ATTRIBUTE SAVE AND FOR CLARITY
    ! SHOULD HAVE THE SAME NAME AS THOSE USED ON THE OTHER PROCESSORS!
    IF(IOPROC .AND. .NOT. IOPROC_ALLOCATED) THEN

       IOPROC_ALLOCATED = .true.

       allocate(U_STOKES_3D(NGL,KB),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:U_STOKES_3D")
       U_STOKES_3D = 0.0_sp

       allocate(V_STOKES_3D(NGL,KB),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:V_STOKES_3D")
       V_STOKES_3D = 0.0_sp

       allocate(WAVESTRX_3D(NGL,KB),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:WAVESTRX_3D")
       WAVESTRX_3D = 0.0_sp

       allocate(WAVESTRY_3D(NGL,KB),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:WAVESTRY_3D")
       WAVESTRY_3D = 0.0_sp

       allocate(WAVESTRX_2D(NGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:WAVESTRX_2D")
       WAVESTRX_3D = 0.0_sp

       allocate(WAVESTRY_2D(NGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:WAVESTRY_2D")
       WAVESTRY_3D = 0.0_sp


    END IF


  ! ALLOCATE THE NEW FILE OBJECT
    NCF => NEW_FILE()


    ! U_STOKES_3D
    VAR  => NC_MAKE_AVAR(name='u_stokes_3d',&
         & values=U_STOKES_3D, DIM1= DIM_nele, DIM2= DIM_siglay, DIM3 = DIM_time)


    ATT  => NC_MAKE_ATT(name='long_name',values='Eastward Stokes Velocity') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='m/s') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)



    ! V_STOKES_3D
    VAR  => NC_MAKE_AVAR(name='v_stokes_3d',&
         & values=V_STOKES_3D, DIM1= DIM_nele, DIM2= DIM_siglay, DIM3 = DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Northward Stokes Velocity') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='m/s') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! U_STRESS_3D
    VAR  => NC_MAKE_AVAR(name='u_stress_3d',&
         & values=WAVESTRX_3D, DIM1= DIM_nele, DIM2= DIM_siglay, DIM3 = DIM_time)


    ATT  => NC_MAKE_ATT(name='long_name',values='Eastward 3-d stress') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='Pa') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)



    ! V_STRESS_3D
    VAR  => NC_MAKE_AVAR(name='v_stress_3d',&
         & values=WAVESTRY_3D, DIM1= DIM_nele, DIM2= DIM_siglay, DIM3 = DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Northward 3-d stress') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='Pa') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)


    ! U_STRESS_2D
    VAR  => NC_MAKE_AVAR(name='u_stress_2d',&
         & values=WAVESTRX_2D, DIM1= DIM_nele, DIM2 = DIM_time)


    ATT  => NC_MAKE_ATT(name='long_name',values='Depth-averaged Eastward 3-d stress') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='Pa') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)



    ! V_STRESS_2D
    VAR  => NC_MAKE_AVAR(name='v_stress_2d',&
         & values=WAVESTRY_2D, DIM1= DIM_nele, DIM2 = DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Depth-averaged Northward 3-d stress') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='Pa') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END WAVE_STRESS_FILE_OBJECT"


  END  FUNCTION WAVE_STRESS_FILE_OBJECT 


# endif


! J. Ge for fluid mud layer
# if defined (FLUID_MUD)
  FUNCTION FLUID_MUD_FILE_OBJECT() RESULT(NCF)
    USE MOD_SED
    IMPLICIT NONE
   INTEGER :: status
   LOGICAL, SAVE :: IOPROC_ALLOCATED1 = .FALSE.
   LOGICAL, SAVE :: IOPROC_ALLOCATED2 = .FALSE.
   LOGICAL, SAVE :: IOPROC_ALLOCATED3 = .FALSE.
   TYPE(NCFILE), POINTER :: NCF
   TYPE(NCVAR),  POINTER :: VAR
   TYPE(NCATT),  POINTER :: ATT
   INTEGER :: II, IERR, len
   INTEGER :: DEST, SOURCE, TAGI,TAGC
   INTEGER, PARAMETER :: TAG = 70077
   INTEGER, PARAMETER :: TAG_I = 71077
   INTEGER, PARAMETER :: TAG_C = 72077
# if defined(MULTIPROCESSOR)
   INTEGER :: STAT(MPI_STATUS_SIZE)
# endif

   REAL(SP),SAVE, POINTER :: SLD_P(:,:)
   REAL(SP),SAVE, POINTER :: BLD_P(:)
   REAL(SP),SAVE, POINTER :: BED_P(:,:)
   REAL(SP),SAVE, POINTER :: BOT_P(:)

   IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START FLUID_MUD_FILE_OBJECT"  
  ! ALLOCATE THE NEW FILE OBJECT
   NCF => NEW_FILE()


  IF(IOPROC .AND. .NOT. IOPROC_ALLOCATED1) THEN

     IOPROC_ALLOCATED1 = .true.

     allocate(D_FMLcm(MGL),stat=status)
     IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:ELF_FML")
     D_FMLcm = 0.0_SP

     allocate(WETNODES_FML(MGL),stat=status)
     IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:ELF_FML")
     WETNODES_FML = 0.0_SP

     allocate(WETCELLS_FML(NGL),stat=status)
     IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:ELF_FML")
     WETCELLS_FML = 0.0_SP

     allocate(WETNODES_CUR(MGL),stat=status)
     IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:ELF_FML")
     WETNODES_CUR = 0.0_SP

     allocate(WETCELLS_CUR(NGL),stat=status)
     IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:ELF_FML")
     WETCELLS_CUR = 0.0_SP

    allocate(UA_FML(NGL),stat=status)
     IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:UA_FML")
     UA_FML = 0.0_sp

     allocate(VA_FML(NGL),stat=status)
     IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:UA_FML")
     VA_FML = 0.0_sp

  END IF


    ! fluid_mud_thickness
    VAR  => NC_MAKE_AVAR(name='fluid_mud_thickness',&
       & values=D_FMLcm, DIM1= DIM_node, DIM2= DIM_time)
         
    ATT  => NC_MAKE_ATT(name='long_name',values='thickness of fluid mud layer')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='note1',values='this is for near bed process')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='note2',values='dimensions are length')
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='units',values='cm')
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data')
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)


    ! wet nodes FML
    VAR  => NC_MAKE_AVAR(name='wet_nodes_fml',&
       & values=WETNODES_FML, DIM1= DIM_node, DIM2= DIM_time)
         
    ATT  => NC_MAKE_ATT(name='long_name',values='wet node of fml')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='note1',values='this is for near bed process')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='note2',values='dimensions are length')
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='units',values='/')
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data')
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)


    ! wet cells FML
    VAR  => NC_MAKE_AVAR(name='wet_cells_fml',&
       & values=WETCELLS_FML, DIM1= DIM_nele, DIM2= DIM_time)
         
    ATT  => NC_MAKE_ATT(name='long_name',values='wet cell of fml')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='note1',values='this is for near bed process')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='note2',values='dimensions are length')
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='units',values='/')
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data')
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)


    ! wet nodes CUR
    VAR  => NC_MAKE_AVAR(name='wet_nodes_cur',&
       & values=WETNODES_CUR, DIM1= DIM_node, DIM2= DIM_time)
         
    ATT  => NC_MAKE_ATT(name='long_name',values='wet node of current')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='note1',values='this is for near bed process')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='note2',values='dimensions are length')
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='units',values='/')
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data')
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)


    ! wet cells CUR
    VAR  => NC_MAKE_AVAR(name='wet_cells_cur',&
       & values=WETCELLS_CUR, DIM1= DIM_nele, DIM2= DIM_time)
         
    ATT  => NC_MAKE_ATT(name='long_name',values='wet cell of current')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='note1',values='this is for near bed process')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='note2',values='dimensions are length')
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='units',values='/')
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data')
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)


    ! fluid mud source/sink
    VAR  => NC_MAKE_AVAR(name='fluid_mud_source_sink',&
       & values=Source_Fluid_Mud, DIM1= DIM_node, DIM2= DIM_time)
         
    ATT  => NC_MAKE_ATT(name='long_name',values='source and sink term of fluid mud')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='note1',values='this is for near bed process')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='note2',values='dimensions are m/s')
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='units',values='m/s')
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data')
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! ua_fml
    VAR  => NC_MAKE_AVAR(name='ua_fml',&
         & values=ua_fml, DIM1= DIM_nele, DIM2 = DIM_time)


    ATT  => NC_MAKE_ATT(name='long_name',values='x-component flow speed of fluid m&
         &ud layer') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='m/s') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)


    ! va_fml
    VAR  => NC_MAKE_AVAR(name='va_fml',&
         & values=va_fml, DIM1= DIM_nele, DIM2 = DIM_time)


    ATT  => NC_MAKE_ATT(name='long_name',values='y-component flow speed of fluid m&
         &ud layer') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='m/s') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)


   IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END FLUID_MUD_FILE_OBJECT"
  End FUNCTION FLUID_MUD_FILE_OBJECT
# endif
! J. Ge for fluid mud layer



!=============================================================
# if defined (SEDIMENT)
  FUNCTION SED_FILE_OBJECT() RESULT(NCF)

#  if defined (ORIG_SED)
  USE mod_sed
#  elif defined (CSTMS_SED)
  USE mod_sed_cstms
#  endif

    IMPLICIT NONE

   INTEGER :: status
   LOGICAL, SAVE :: IOPROC_ALLOCATED1 = .FALSE.
   LOGICAL, SAVE :: IOPROC_ALLOCATED2 = .FALSE.
   LOGICAL, SAVE :: IOPROC_ALLOCATED3 = .FALSE.
   TYPE(NCFILE), POINTER :: NCF
   TYPE(NCVAR),  POINTER :: VAR
   TYPE(NCATT),  POINTER :: ATT
   INTEGER :: II, IERR, len
   INTEGER :: DEST, SOURCE, TAGI,TAGC
   INTEGER, PARAMETER :: TAG = 70077
   INTEGER, PARAMETER :: TAG_I = 71077
   INTEGER, PARAMETER :: TAG_C = 72077
# if defined(MULTIPROCESSOR)
   INTEGER :: STAT(MPI_STATUS_SIZE)
# endif

   REAL(SP),SAVE, POINTER :: SLD_P(:,:)
   REAL(SP),SAVE, POINTER :: BLD_P(:)
   REAL(SP),SAVE, POINTER :: BED_P(:,:)
   REAL(SP),SAVE, POINTER :: BOT_P(:)

   IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START SED_FILE_OBJECT"
   
  ! ALLOCATE THE NEW FILE OBJECT
   NCF => NEW_FILE()


  !-------------------------------------------------- 
  !bottom stress on the nodes - includes waves
  !-------------------------------------------------- 
  IF(IOPROC .AND. .NOT. IOPROC_ALLOCATED1) THEN

     IOPROC_ALLOCATED1 = .true.

     allocate(TAUBM_N(MGL),stat=status)
     IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:TAUB_N")
     TAUBM_N = 0.0_SP

  END IF

  VAR  => NC_MAKE_AVAR(name='taub_total',&
       & values=taubm_n, DIM1= DIM_node, DIM2= DIM_time)
         
  ATT  => NC_MAKE_ATT(name='long_name',values='total bed stress magnitude')
  VAR  => ADD(VAR,ATT)

  ATT  => NC_MAKE_ATT(name='note1',values='this stress is used to force sediment')
  VAR  => ADD(VAR,ATT)

  ATT  => NC_MAKE_ATT(name='note2',values='dimensions are stress/rho')
  VAR  => ADD(VAR,ATT)
    
  ATT  => NC_MAKE_ATT(name='units',values='m^2 s^-2')
  VAR  => ADD(VAR,ATT)
    
  ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid')
  VAR  => ADD(VAR,ATT)

  ATT  => NC_MAKE_ATT(name='type',values='data')
  VAR  => ADD(VAR,ATT)

  NCF  => ADD(NCF,VAR)

  !-------------------------------------------------- 
  !sediment concentration
  !-------------------------------------------------- 
  ! IO PROC MUST ALLOCATE SPACE FOR THE ARRAYS 
  ! THESE ARRAYS MUST HAVE THE ATTRIBUTE SAVE AND FOR CLARITY
  ! SHOULD HAVE THE SAME NAME AS THOSE USED ON THE OTHER PROCESSORS!
# if defined(MULTIPROCESSOR)
  IF(USE_MPI_IO_MODE) THEN

     IF(MSR) THEN
        DEST = IOPROCID - 1
        CALL MPI_SEND(NSED,1,MPI_INTEGER,DEST,TAG,MPI_COMM_WORLD,IERR)


        DO II = 1,nsed
           
           len = len_trim(sed(ii)%sname)
           CALL MPI_SEND(len,1,MPI_INTEGER,DEST,TAG_I+II,MPI_COMM_WORLD,IERR)
           
           CALL MPI_SEND(sed(ii)%sname(1:len),len,MPI_CHARACTER,DEST,TAG_C+II,MPI_COMM_WORLD,IERR)
           
        END DO

        

     ELSEIF(IOPROC) THEN
        SOURCE = MSRID -1
        CALL MPI_RECV(NSED,1,MPI_INTEGER,SOURCE,TAG,MPI_COMM_WORLD,STAT,IERR)

        DO II = 1,nsed
           CALL MPI_RECV(len,1,MPI_INTEGER,SOURCE,TAG_I+II,MPI_COMM_WORLD,STAT,IERR)
           
           CALL MPI_RECV(sed(ii)%sname,len,MPI_CHARACTER,SOURCE,TAG_C+II,MPI_COMM_WORLD,STAT,IERR)
           
           sed(ii)%sname = sed(ii)%sname(1:len)
           
        END DO


     END IF

  END IF
# endif

# if defined(MULTIPROCESSOR)
  IF(USE_MPI_IO_MODE) THEN

     IF(MSR) THEN
        DEST = IOPROCID - 1
        CALL MPI_SEND(NSED,1,MPI_INTEGER,DEST,TAG,MPI_COMM_WORLD,IERR)


        DO II = 1,nsed

           len = len_trim(sed(ii)%sname2)
           CALL MPI_SEND(len,1,MPI_INTEGER,DEST,TAG_I+II,MPI_COMM_WORLD,IERR)

           CALL MPI_SEND(sed(ii)%sname2(1:len),len,MPI_CHARACTER,DEST,TAG_C+II,MPI_COMM_WORLD,IERR)

        END DO



     ELSEIF(IOPROC) THEN
        SOURCE = MSRID -1
        CALL MPI_RECV(NSED,1,MPI_INTEGER,SOURCE,TAG,MPI_COMM_WORLD,STAT,IERR)

        DO II = 1,nsed
           CALL MPI_RECV(len,1,MPI_INTEGER,SOURCE,TAG_I+II,MPI_COMM_WORLD,STAT,IERR)

           CALL MPI_RECV(sed(ii)%sname2,len,MPI_CHARACTER,SOURCE,TAG_C+II,MPI_COMM_WORLD,STAT,IERR)

           sed(ii)%sname2 = sed(ii)%sname2(1:len)

        END DO


     END IF

  END IF
# endif

  if(bedload)then
  Do II=1,nsed        
    BLD_P => sed(ii)%bflx(1:m) 

    VAR  => NC_MAKE_PVAR(name=TRIM(sed(ii)%sname2),&
        & values=BLD_P, DIM1= DIM_node,  DIM2 = DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values=TRIM('Bedload Mass '//TRIM(sed(ii)%sname))) 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='kg') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    NULLIFY(BLD_P)
  
  END DO
  endif

  if(susload)then
  Do II=1,nsed        
    SLD_P => sed(ii)%conc(1:m,1:kbm1) 

    VAR  => NC_MAKE_PVAR(name=TRIM(sed(ii)%sname),&
        & values=SLD_P, DIM1= DIM_node, DIM2= DIM_siglay, DIM3 = DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values=TRIM('Suspended Sed Conc of '//TRIM(sed(ii)%sname))) 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='g/L') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    NULLIFY(SLD_P)
  
  END DO

    if(nsed>1)then
      VAR  => NC_MAKE_AVAR(name='todal_ssc', &
          & values=csed, DIM1= DIM_node, DIM2= DIM_siglay, DIM3 = DIM_time)

      ATT  => NC_MAKE_ATT(name='long_name',values=TRIM('Total Suspended Sed &
           &Conc of all sediment classes')) 
      VAR  => ADD(VAR,ATT)

      ATT  => NC_MAKE_ATT(name='units',values='g/L') 
      VAR  => ADD(VAR,ATT)

      ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
      VAR  => ADD(VAR,ATT)

      NCF  => ADD(NCF,VAR)
    endif

  endif

  Do II=1,nsed
    BLD_P => sed(ii)%frac(1:m,1)

    VAR  => NC_MAKE_PVAR(name=TRIM(sed(ii)%sname)//'_bedfrac',&
        & values=BLD_P, DIM1= DIM_node, DIM2 = DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values=TRIM('Fraction of '//TRIM(sed(ii)%sname)//' in surface layer'))
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='-')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid')
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    NULLIFY(BLD_P)

  END DO


  !-------------------------------------------------- 
  !Bed Params (thickness, age, and porosity)
  !-------------------------------------------------- 
  if(sed_dumpbed)then
    
    ! IO PROC MUST ALLOCATE SPACE FOR THE ARRAYS 
    ! THESE ARRAYS MUST HAVE THE ATTRIBUTE SAVE AND FOR CLARITY
    ! SHOULD HAVE THE SAME NAME AS THOSE USED ON THE OTHER PROCESSORS!
    IF(IOPROC .AND. .NOT. IOPROC_ALLOCATED2) THEN
 
      IOPROC_ALLOCATED2 = .true.
 
      ALLOCATE(BED(MGL,NBED,N_BED_CHARS),STAT=STATUS)
      IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:SED_BED")
      BED = 0.0_SP
 
    END IF

    Do II=1,n_bed_chars
      BED_P => bed(1:m,1:nbed,ii) 
  
      VAR  => NC_MAKE_PVAR(name=TRIM(BED_SNAMES(II)),&
          & values=BED_P, DIM1= DIM_node, DIM2= DIM_bedlay, DIM3 = DIM_time)
  
      ATT  => NC_MAKE_ATT(name='long_name',values=TRIM(BED_LNAMES(II))) 
      VAR  => ADD(VAR,ATT)
  
      ATT  => NC_MAKE_ATT(name='units',values=TRIM(BED_UNITS(II))) 
      VAR  => ADD(VAR,ATT)
  
      ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
      VAR  => ADD(VAR,ATT)

      NCF  => ADD(NCF,VAR)
  
      NULLIFY(BED_P)
    
    END DO
  endif !sed_dumpbed   

  !-------------------------------------------------- 
  !Bottom Params 
  !-------------------------------------------------- 
  if(sed_dumpbot)then
    
    ! IO PROC MUST ALLOCATE SPACE FOR THE ARRAYS 
    ! THESE ARRAYS MUST HAVE THE ATTRIBUTE SAVE AND FOR CLARITY
    ! SHOULD HAVE THE SAME NAME AS THOSE USED ON THE OTHER PROCESSORS!
    IF(IOPROC .AND. .NOT. IOPROC_ALLOCATED3) THEN
 
      IOPROC_ALLOCATED3 = .true.
 
      ALLOCATE(BOTTOM(MGL,N_BOT_VARS),STAT=STATUS)
      IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:SED_BOT")
      BOTTOM = 0.0_SP
 
    END IF

    Do II=1,n_bot_vars  
      BOT_P => bottom(1:m,ii) 

      VAR  => NC_MAKE_PVAR(name=TRIM(BOT_SNAMES(II)),&
        & values=BOT_P, DIM1= DIM_node,  DIM2 = DIM_time)

      ATT  => NC_MAKE_Att(name='long_name',values=TRIM(BOT_LNAMES(II))) 
      VAR  => ADD(VAR,ATT)

      ATT  => NC_MAKE_Att(name='units',values=TRIM(BOT_UNITS(II))) 
      VAR  => ADD(VAR,ATT)

      ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
      VAR  => ADD(VAR,ATT)

      NCF  => ADD(NCF,VAR)

      NULLIFY(BOT_P)
  
    End Do 
  Endif ! sed_dumpbot 

   IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END SED_FILE_OBJECT"

  END FUNCTION SED_FILE_OBJECT
# endif


!=============================================================  
# if defined (DYE_RELEASE)
  FUNCTION DYE_FILE_OBJECT() RESULT(NCF)
    USE MOD_DYE

    IMPLICIT NONE

    INTEGER :: status, I
    LOGICAL, SAVE :: IOPROC_ALLOCATED = .FALSE.
    TYPE(NCFILE), POINTER :: NCF
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCATT),  POINTER :: ATT


    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START DYE_FILE_OBJECT"

    ! IO PROC MUST ALLOCATE SPACE FOR THE ARRAYS 
    ! THESE ARRAYS MUST HAVE THE ATTRIBUTE SAVE AND FOR CLARITY
    ! SHOULD HAVE THE SAME NAME AS THOSE USED ON THE OTHER PROCESSORS!
    IF(IOPROC .AND. .NOT. IOPROC_ALLOCATED) THEN

       IOPROC_ALLOCATED = .true.

       allocate(DYE(MGL,KB),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:EL")
       DYE = 0.0_SP

    END IF

    ! ALLOCATE THE NEW FILE OBJECT
    NCF => NEW_FILE()


    ! Dye
    VAR  => NC_MAKE_AVAR(name='DYE',&
         & values=DYE, DIM1= DIM_node, DIM2= DIM_siglay, DIM3=DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Water Surface Elevation') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='meters') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='positive',values='up') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='standard_name',values='sea_surface_height_above_geoid') 
    VAR  => ADD(VAR,ATT)

# if defined(UCF)
    ATT  => NC_MAKE_ATT(name='grid',values='grid4') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates',values="lat lon") 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid_location',values='node') 
    VAR  => ADD(VAR,ATT)
# else
    ATT  => NC_MAKE_ATT(name='grid',values='SigmaLayer_Mesh') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar) 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)
# endif

    NCF  => ADD(NCF,VAR)


    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END DYE_FILE_OBJECT"

  END FUNCTION DYE_FILE_OBJECT
# endif

!=============================================================  

SUBROUTINE UPDATE_IODATA(NCF,NOW)
  IMPLICIT NONE
  TYPE(TIME), INTENT(IN) :: NOW
  TYPE(NCVAR), POINTER :: VAR1,VAR2
  TYPE(NCFILE), POINTER :: NCF
  LOGICAL :: FOUND
!====================================================================
! THIS SUBROUTINE IS IN CHARGE OF UPDATING ANY VARIABLES FOR IO WHICH
! ARE NOT ALREADY UPDATED BY FVCOM DURING THE MAIN LOOP:
!====================================================================
  IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START UPDATE_IODATA"
  
  if(.not. Associated(NCF)) CALL FATAL_ERROR&
       &("UPDATE_IODATA: THE FILE OBJECT IS NOT ASSOCIATED!")
  
  VAR1 => FIND_VAR(NCF,"time",FOUND)
  IF(FOUND) CALL UPDATE_FLOAT_TIME(VAR1,NOW)

  VAR1 => FIND_VAR(NCF,"Itime",FOUND)
  IF(FOUND) THEN
     VAR2 => FIND_VAR(NCF,"Itime2",FOUND)
     IF (.NOT.FOUND) THEN
        CALL WARNING&
             & ("FOUND ONLY PART OF INTEGER TIME VARIABLE IN OUT PUT FILE!")
     ELSE
        CALL UPDATE_ITIME(VAR1,VAR2,NOW)
     END IF
  END IF

  VAR1 => FIND_VAR(NCF,"Times",FOUND)
  IF(FOUND) CALL UPDATE_DATETIME(VAR1,NOW)


  ! IINT IS A LONG LONG INTEGER BUT WE CAN'T WRITE LONG INTEGERS...
  VAR1 => FIND_VAR(NCF,"iint",FOUND)
  IF(FOUND) VAR1%scl_int = IINT

  VAR1 => FIND_VAR(NCF,"file_date",FOUND)
  IF(FOUND) CALL UPDATE_DATETIME(VAR1,GET_NOW())

  IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END UPDATE_IODATA"

END SUBROUTINE UPDATE_IODATA
!=============================================================  
!
!-------------------------------------------------------------
!
!=============================================================  
! PROGRAMS FOR THE OUTPUT OF SUBDOMAINS IN FVCOM
!=============================================================  
SUBROUTINE SETUP_SUBDOMAINS(SUB_FILES,GRIDS)
  IMPLICIT NONE
  TYPE(NCFILE), POINTER :: NCF
  TYPE(GRID), POINTER :: GRIDS(:)
  CHARACTER(LEN=*) :: SUB_FILES
  CHARACTER(LEN=80),ALLOCATABLE :: FNAMES(:)
  CHARACTER(LEN=80) :: FILE,PATH,EXT
  INTEGER :: NUMF,I, STATUS

  INTEGER, POINTER:: NID(:),EID(:)

  if(DBG_SET(dbg_sbr))  write(IPT,*) "START Setup_subdomains"



  IF(SERIAL .and. SUB_FILES /= "FVCOM") THEN 
     CALL WARNING&
          &("SETUP_SUBDOMAIN: YOU CAN NOT USE SUBDOMAIN OUTPUT DURING A SINGLE PROCESSOR MODLE RUN!")
     SUB_FILES = "FVCOM"
  ELSEIF(len_trim(SUB_FILES) == 0 ) THEN
     CALL FATAL_ERROR("THE SUBDOMAIN FILE LIST PASSED TO SETUP_SUBDOMAINS IS EMPTY",&
          &"PLEASE CHECK YOUR NAME LIST FILE FOR THE NC_SUBDOMAIN_FILES",&
          &"AND NCAV_SUBDOMAIN_FILES ENTRIES!")
     
  END IF



  CALL SPLIT_STRING(SUB_FILES,",",FNAMES)
  NUMF =SIZE(FNAMES)

  
  if(DBG_SET(dbg_log))  write(IPT,*) "! NUMBER OF DOMAINS TO OUTPUT:",NUMF
  ALLOCATE(GRIDS(NUMF),STAT=status)
  IF(STATUS /=0) CALL FATAL_ERROR("NCDIO: COULD NOT ALLOCATE SUBDOMAIN GRIDS!")

#  if defined (MULTIPROCESSOR)
  DO I = 1, NUMF

     IF(FNAMES(I) == "FVCOM") THEN
        if(DBG_SET(dbg_log))  write(IPT,*) "! SETTING FVCOM DOMAIN OUTPUT"
        if(DBG_SET(dbg_log))  write(IPT,*) "! DIMENSIONS: MGL =",MGL
        if(DBG_SET(dbg_log))  write(IPT,*) "! DIMENSIONS: NGL =",NGL

        CALL SET_FVCOM_GRID(GRIDS(I))
        CYCLE
     END IF

     if(DBG_SET(dbg_log))  write(IPT,*) "! READING SUBDOMAIN FILE:"//TRIM(FNAMES(I))
  
     ! FIND THE LOCAL CELLS IN EACH SUBDOMAIN FILE
     CALL LOAD_SUBDOMAIN_FILE(FNAMES(I),EID,NID)
     ! SPECIFY BY BOX or RADIUS

     ! SET ELID/NLID FOR THE SUBDOMAIN BASED ON THE LOCAL MEMBERS
     CALL GENMAP_SUBDOMAIN(EID,NID,GRIDS(I))

     ! SET THE INDEX ARRAYS FOR THE GRID
     CALL SET_SUBDOMAIN_GRID(GRIDS(I))

     CALL PATH_SPLIT(FNAMES(I),PATH,FILE,EXT)

     GRIDS(I)%NAME=FILE


     if(DBG_SET(dbg_log))  write(IPT,*) "! SET SUBDOMAIN OUTPUT"
     if(DBG_SET(dbg_log))  write(IPT,*) "! DIMENSIONS: MGL =",GRIDS(I)%MGL
     if(DBG_SET(dbg_log))  write(IPT,*) "! DIMENSIONS: NGL =",GRIDS(I)%NGL
     

     IF(GRIDS(I)%MGL <1 .or. GRIDS(I)%NGL <1) CALL FATAL_ERROR &
          &("SUBDOMAIN USER ERROR: NO VALID CELLS OR NODES WERE FOUND",&
          & "IN THE SPECIFIED REGION!")

  END DO
# else

  if(DBG_SET(dbg_log))  write(IPT,*) "! SETTING FVCOM DOMAIN OUTPUT"
  if(DBG_SET(dbg_log))  write(IPT,*) "! DIMENSIONS: MGL =",MGL
  if(DBG_SET(dbg_log))  write(IPT,*) "! DIMENSIONS: NGL =",NGL
  
  CALL SET_FVCOM_GRID(GRIDS(1))

# endif

END SUBROUTINE SETUP_SUBDOMAINS
#  if defined (MULTIPROCESSOR)
SUBROUTINE GENMAP_SUBDOMAIN(EID_L,NID_L,LG)
!==============================================================================|
!
! CREATE A GLOBAL TO LOCAL MAP FOR SUBDOMAIN OUTPUT
! USES DATA READ INTO: 
!                     
! Creates:             MAP LINK LIST ENTRY FOR IO
!                      
!
!==============================================================================|
!  USE MOD_NESTING
  USE MOD_PAR
  IMPLICIT NONE
  TYPE(GRID) :: LG !The Local Grid

  integer :: SENDER,RECVER, ierr, I, NCNT, NSZE, I1, status,J,lb,ub

  INTEGER, POINTER :: EID_L(:),NID_L(:)

  INTEGER, POINTER :: EID(:),NID(:)

  INTEGER, POINTER :: ESD_GL(:), NSD_GL(:)

  INTEGER, POINTER :: TEMP1(:),TEMP2(:)

  TYPE(MAP), POINTER :: E_MAP(:),N_MAP(:)


  if (dbg_set(dbg_SBR)) &
       & write(IPT,*) "START: GENMAP_SUBDOMAIN"


  IF(.not. IOPROC) THEN !ONLY THE FVCOM GROUP PROCS CAN DO THIS
     
     ! COLLECT TO GLOBAL ON ALL PROCS
     ALLOCATE(EID(NGL),STAT=STATUS)
     IF(STATUS/=0) CALL FATAL_ERROR("GENMAP_SUBDOMAIN COULD NOT ALLOCATE EID")
     EID = 0

     CALL PCOLLECT(MYID,MSRID,NPROCS,EMAP,EID_L,EID)

     DEALLOCATE(EID_L)

     SENDER = MSRID -1
     CALL MPI_BCAST(EID,NGL,MPI_INTEGER,SENDER,MPI_FVCOM_GROUP,IERR)
     
     NCNT = SUM(EID)
     LG%NGL = NCNT ! SET DIMENSION

     ! MAKE A LIST OF THE GLOBAL ELEMENTS IN THE SUBDOMAIN
     ALLOCATE(ESD_GL(NCNT))
     I1 = 1
     DO I = 1,NGL
        IF(EID(I) ==1) THEN
           ESD_GL(I1) = I
           I1 = I1 +1
        END IF
     END DO

     DEALLOCATE(EID)

     
     ! COLLECT TO GLOBAL ON ALL PROCS
     ALLOCATE(NID(MGL),STAT=STATUS)
     IF(STATUS/=0) CALL FATAL_ERROR("GENMAP_SUBDOMAIN COULD NOT ALLOCATE NID")
     NID = 0

     CALL PCOLLECT(MYID,MSRID,NPROCS,NMAP,NID_L,NID)

     DEALLOCATE(NID_L)

     SENDER = MSRID -1
     CALL MPI_BCAST(NID,MGL,MPI_INTEGER,SENDER,MPI_FVCOM_GROUP,IERR)
     
     NCNT = SUM(NID)
     LG%MGL = NCNT ! SET DIMENSION


     ! MAKE A LIST OF THE GLOBAL ELEMENTS IN THE SUBDOMAIN
     ALLOCATE(NSD_GL(NCNT))
     I1 = 1
     DO I = 1,MGL
        IF(NID(I) ==1) THEN
           NSD_GL(I1) = I
           I1 = I1 +1
        END IF
     END DO

     DEALLOCATE(NID)

  END IF

  ! SET DIMENSIONS IN THE GRID OBJECT
  LG%KB = KB
  LG%KBM1 = KBM1

  ! SET FOR IOPROC!
  SENDER = MSRID -1
  CALL MPI_BCAST(LG%MGL,1,MPI_INTEGER,SENDER,MPI_COMM_WORLD,IERR)
  CALL MPI_BCAST(LG%NGL,1,MPI_INTEGER,SENDER,MPI_COMM_WORLD,IERR)


  IF(.NOT. IOPROC) THEN
!============================================
! Make a list of the local Nesting nodes
!============================================

     !!SET UP LOCAL NESTING NODES
     ALLOCATE(TEMP1(0:LG%MGL));      TEMP1=0
     ALLOCATE(TEMP2(0:LG%MGL));      TEMP2=0
     
     NCNT = 0
     DO I=1,LG%MGL
        I1 = NLID( NSD_GL(I) )
        IF(I1 /= 0)THEN
           NCNT = NCNT + 1
           TEMP1(NCNT) = I1
           TEMP2(NCNT) = I
        END IF
     END DO
     
     ! SET LOCAL NUMBER OF NESTING NODES
     LG%M = NCNT

     ! SET GLOBAL TO LOCAL MAP FOR THIS DOMAIN
     ALLOCATE(LG%NGID(0:LG%M),stat=status)
     if(status /= 0) call fatal_error("GENMAP_NESTING: can not allocate:LG%NGID")

     ALLOCATE(LG%NLID(0:LG%M),stat=status)
     if(status /= 0) call fatal_error("GENMAP_NESTING: can not allocate:LG%NLID")

     LG%NLID = TEMP1(0:NCNT)

     LG%NGID = TEMP2(0:NCNT)
     
     DEALLOCATE(TEMP1)
     DEALLOCATE(TEMP2)



     !!SET UP LOCAL + HALO NESTING NODES
     ALLOCATE(TEMP1(0:LG%MGL));      TEMP1=0
     ALLOCATE(TEMP2(0:LG%MGL));      TEMP2=0
     
     NCNT = 0
     DO I=1,LG%MGL
        I1 = NLID_X( NSD_GL(I) )
        IF(I1 > M)THEN
           NCNT = NCNT + 1
           TEMP1(NCNT) = I1
           TEMP2(NCNT) = I
        END IF
     END DO
     
     ! SET LOCAL NUMBER OF NESTING NODES
     LG%MT = LG%M + NCNT

     ! SET GLOBAL TO LOCAL MAP FOR THIS DOMAIN
     ALLOCATE(LG%NGID_X(0:LG%MT),stat=status)
     if(status /= 0) call fatal_error("GENMAP_NESTING: can not allocate:LG%NGID")

     ALLOCATE(LG%NLID_X(0:LG%MT),stat=status)
     if(status /= 0) call fatal_error("GENMAP_NESTING: can not allocate:LG%NLID")

     LG%NLID_X(0:LG%M) = LG%NLID

     LG%NGID_X(0:LG%M) = LG%NGID

     lb = LG%M+1
     ub = LG%MT
     LG%NLID_X(lb:ub) = TEMP1(1:NCNT)

     LG%NGID_X(lb:ub) = TEMP2(1:NCNT)

     DEALLOCATE(TEMP1)
     DEALLOCATE(TEMP2)
     DEALLOCATE(NSD_GL)


!============================================
! Make a list of the local Nesting elements
!============================================

     !!SET UP LOCAL OPEN BOUNDARY NODES
     ALLOCATE(TEMP1(0:LG%NGL));      TEMP1=0
     ALLOCATE(TEMP2(0:LG%NGL));      TEMP2=0
     
     NCNT = 0
     DO I=1,LG%NGL
        I1 =  ELID(ESD_GL(I)) 
        IF(I1 /= 0)THEN
           NCNT = NCNT + 1
           TEMP1(NCNT) = I1
           TEMP2(NCNT) = I
        END IF
     END DO
     
     ! SET LOCAL NUMBER OF NESTING ELEMENTS
     LG%N = NCNT

     ! SET GLOBAL TO LOCAL MAP FOR THIS DOMAIN
     ALLOCATE(LG%EGID(0:LG%N),stat=status)
     if(status /= 0) call fatal_error("GENMAP_NESTING: can not allocate:LG%NGID")

     ALLOCATE(LG%ELID(0:LG%N),stat=status)
     if(status /= 0) call fatal_error("GENMAP_NESTING: can not allocate:LG%NLID")

     LG%ELID = TEMP1(0:NCNT)

     LG%EGID = TEMP2(0:NCNT)
     
     DEALLOCATE(TEMP1)
     DEALLOCATE(TEMP2)


     !!SET UP LOCAL+HALO OPEN BOUNDARY NODES
     ALLOCATE(TEMP1(0:LG%NGL));      TEMP1=0
     ALLOCATE(TEMP2(0:LG%NGL));      TEMP2=0
     
     NCNT = 0
     DO I=1,LG%NGL
        I1 =  ELID_X(ESD_GL(I)) 
        IF(I1 > N)THEN
           NCNT = NCNT + 1
           TEMP1(NCNT) = I1
           TEMP2(NCNT) = I
        END IF
     END DO
     
     ! SET LOCAL NUMBER OF NESTING ELEMENTS
     LG%NT = LG%N + NCNT

     ! SET GLOBAL TO LOCAL MAP FOR THIS DOMAIN
     ALLOCATE(LG%EGID_X(0:LG%NT),stat=status)
     if(status /= 0) call fatal_error("GENMAP_NESTING: can not allocate:LG%NGID")

     ALLOCATE(LG%ELID_X(0:LG%NT),stat=status)
     if(status /= 0) call fatal_error("GENMAP_NESTING: can not allocate:LG%NLID")

     LG%ELID_X(0:LG%N) = LG%ELID

     LG%EGID_X(0:LG%N) = LG%EGID

     lb = LG%N+1
     ub = LG%NT
     LG%ELID_X(lb:ub) = TEMP1(1:NCNT)

     LG%EGID_X(lb:ub) = TEMP2(1:NCNT)
     
     DEALLOCATE(TEMP1)
     DEALLOCATE(TEMP2)
     DEALLOCATE(ESD_GL)


  END IF ! END IF IOPROC

  !==============================================================================|
  !   SET UP ELEMENT MAPPING FOR GLOBAL 2 LOCAL TRANSFER OF BC'S                 | 
  !   BOUNDARY MAP :: BCMAP(NPROCS)                                              |
  !     BCMAP(1-->NPROCS)%NSIZE  :: NUMBER OF BOUNDARY NODES IN EACH DOM         |
  !     BCMAP(1-->NPROCS)%LOC_2_GL(NSIZE) :: LOCAL TO GLOBAL MAPPING IN EACH DOM |
  !==============================================================================|

  ! ELEMENTS: HALO
  E_MAP => MAKE_MAP(MYID,NPROCS,LG%NGL,NT,LG%EGID_X,LG%ELID_X)
  CALL ADD_MAP2LIST(INTERNAL_MAPS,E_MAP)

  NULLIFY(E_MAP)

  E_MAP => MAKE_MAP(MYID,NPROCS,LG%NGL,LG%NT,LG%EGID_X)
  CALL ADD_MAP2LIST(INTERNAL_MAPS,E_MAP)

  NULLIFY(E_MAP)

  ! ELEMENTS: INTERNAL
  E_MAP => MAKE_MAP(MYID,NPROCS,LG%NGL,N,LG%EGID,LG%ELID)
  CALL ADD_MAP2LIST(INTERNAL_MAPS,E_MAP)

  NULLIFY(E_MAP)

  E_MAP => MAKE_MAP(MYID,NPROCS,LG%NGL,LG%N,LG%EGID)
  CALL ADD_MAP2LIST(INTERNAL_MAPS,E_MAP)

  NULLIFY(E_MAP)


  

  ! NODES: HALO
  N_MAP => MAKE_MAP(MYID,NPROCS,LG%MGL,MT,LG%NGID_X,LG%NLID_X)
  CALL ADD_MAP2LIST(INTERNAL_MAPS,N_MAP)

  NULLIFY(N_MAP)


  N_MAP => MAKE_MAP(MYID,NPROCS,LG%MGL,LG%MT,LG%NGID_X)
  CALL ADD_MAP2LIST(INTERNAL_MAPS,N_MAP)

  NULLIFY(N_MAP)

  ! NODES: INTERNAL
  N_MAP => MAKE_MAP(MYID,NPROCS,LG%MGL,M,LG%NGID,LG%NLID)
  CALL ADD_MAP2LIST(INTERNAL_MAPS,N_MAP)

  NULLIFY(N_MAP)


  N_MAP => MAKE_MAP(MYID,NPROCS,LG%MGL,LG%M,LG%NGID)
  CALL ADD_MAP2LIST(INTERNAL_MAPS,N_MAP)

  NULLIFY(N_MAP)

! DO NOT DEALLOCATE MAP!!!

  if (dbg_set(dbg_sbr)) &
       & write(IPT,*) "END: GENMAP_SUBDOMAIN"   
  RETURN
END SUBROUTINE GENMAP_SUBDOMAIN
!========================================================================
!--------------------------------------------------------------------
!========================================================================
SUBROUTINE LOAD_SUBDOMAIN_FILE(FNAME,EID,NID)
  IMPLICIT NONE
  CHARACTER(LEN=*) :: FNAME
  INTEGER :: I, STATUS, ISCAN, SENDER,IERR
  INTEGER, POINTER:: EID(:),NID(:)
  REAL(SP), POINTER:: DIST(:)

  CHARACTER(LEN=80) ::  STR1,STR2
  CHARACTER(LEN=160) :: PATHNFILE

  CHARACTER(LEN=10) :: temp

  INTEGER :: UNITS =0
  INTEGER, PARAMETER :: MTRS=1
  INTEGER, PARAMETER :: DGRS=2

  INTEGER :: MODE =0
  INTEGER, PARAMETER :: rds=1
  INTEGER, PARAMETER :: box=2


  REAL(SP) :: BNDS_LR(2) = 0.0_SP
  REAL(SP) :: BNDS_TB(2) = 0.0_SP

  REAL(SP) :: CENTER(2) = 0.0_SP
  REAL(SP) :: RADIUS = 0.0_SP


  if (dbg_set(dbg_sbr)) &
       & write(IPT,*) "START: LOAD_SUBDOMAIN_FILE"   

  IF(IOPROC) RETURN


  IF(MSR) THEN
     PATHNFILE = TRIM(INPUT_DIR)//TRIM(FNAME)
     CALL FOPEN(SUBDUNIT,trim(pathnfile),'cfr')
     
     ! GET SUBDOMAIN MODE
     ISCAN = SCAN_FILE(SUBDUNIT,"SUBDOMAIN MODE",CVAL = STR1)
     IF(ISCAN /= 0) then
        write(temp,'(I2)') ISCAN
        call fatal_error('Improper formatting of SUBDOMAIN FILE:'//TRIM(FNAME),&
             & 'ISCAN ERROR# '//trim(temp),&
             & 'The header must contain: "SUBDOMAIN MODE"', &
             & 'Followed by one of two defined types: "box" or "radius" ')
     END IF
  
     IF(STR1 == 'box') THEN
        MODE = box
     ELSEIF(STR1=='radius') THEN
        MODE = rds
     ELSE
        CALL FATAL_ERROR&
             &('Improper formatting of SUBDOMAIN FILE:'//TRIM(FNAME),&
             & 'The header must contain: "SUBDOMAIN MODE"', &
             & 'Followed by one of two defined types: "box" or "radius" ')
     END IF

     ! GET SUBDOMAIN MODE
     ISCAN = SCAN_FILE(SUBDUNIT,"UNITS",CVAL = STR1)
     IF(ISCAN /= 0) then
        write(temp,'(I2)') ISCAN
        call fatal_error('Improper formatting of SUBDOMAIN FILE:'//TRIM(FNAME),&
             & 'ISCAN ERROR# '//trim(temp),&
             & 'The header must contain: "UNITS"', &
             & 'Followed by one of two defined types: "meters" or "degrees" ')
     END IF
  
     IF(STR1 == 'meters') THEN
        UNITS = mtrs
     ELSEIF(STR1=='degrees') THEN
        UNITS = dgrs
     ELSE
        CALL FATAL_ERROR&
             &('Improper formatting of SUBDOMAIN FILE:'//TRIM(FNAME),&
             & 'The header must contain: "UNITS"', &
             & 'Followed by one of two defined types: "meters" or "degrees" ')
     END IF
     

     SELECT CASE(MODE)
     CASE(BOX)

        ISCAN = SCAN_FILE(SUBDUNIT,"TOP BOTTOM",FVEC = BNDS_TB ,NSZE = I)
        IF(ISCAN /= 0) then
           write(temp,'(I2)') ISCAN
           call fatal_error('Improper formatting of SUBDOMAIN FILE:'//TRIM(FNAME),&
                & 'ISCAN ERROR# '//trim(temp),&
                & 'For SUBDOMAIN MODE "box", The header must conatain "TOP BOTTOM"',&
                & 'Followed by 2 real values')
        END IF
        IF(I /= 2) CALL FATAL_ERROR&
             &('Improper formatting of SUBDOMAIN FILE:'//TRIM(FNAME),&
                & 'For SUBDOMAIN MODE "box", The header must conatain "TOP BOTTOM"',&
                & 'Followed by 2 real values')


        ISCAN = SCAN_FILE(SUBDUNIT,"LEFT RIGHT",FVEC = BNDS_LR ,NSZE = I)
        IF(ISCAN /= 0) then
           write(temp,'(I2)') ISCAN
           call fatal_error('Improper formatting of SUBDOMAIN FILE:'//TRIM(FNAME),&
                & 'ISCAN ERROR# '//trim(temp),&
                & 'For SUBDOMAIN MODE "box", The header must conatain "LEFT RIGHT"',&
                & 'Followed by 2 real values')
        END IF
        IF(I /= 2) CALL FATAL_ERROR&
             &('Improper formatting of SUBDOMAIN FILE:'//TRIM(FNAME),&
                & 'For SUBDOMAIN MODE "box", The header must conatain "LEFT RIGHT"',&
                & 'Followed by 2 real values')



     CASE(RDS)

        ISCAN = SCAN_FILE(SUBDUNIT,"CENTER",FVEC = CENTER ,NSZE = I)
        IF(ISCAN /= 0) then
           write(temp,'(I2)') ISCAN
           call fatal_error('Improper formatting of SUBDOMAIN FILE:'//TRIM(FNAME),&
                & 'ISCAN ERROR# '//trim(temp),&
                & 'For SUBDOMAIN MODE "radius", The header must conatain "CENTER"',&
                & 'Followed by 2 real values')
        END IF
        IF(I /= 2) CALL FATAL_ERROR&
             &('Improper formatting of SUBDOMAIN FILE:'//TRIM(FNAME),&
                & 'For SUBDOMAIN MODE "box", The header must conatain "CENTER"',&
                & 'Followed by 2 real values')


        ISCAN = SCAN_FILE(SUBDUNIT,"RADIUS",FSCAL = RADIUS)
        IF(ISCAN /= 0) then
           write(temp,'(I2)') ISCAN
           call fatal_error('Improper formatting of SUBDOMAIN FILE:'//TRIM(FNAME),&
                & 'ISCAN ERROR# '//trim(temp),&
                & 'For SUBDOMAIN MODE "radius", The header must conatain "RADIUS"',&
                & 'Followed by a real value')
        END IF


     END SELECT


     CLOSE(SUBDUNIT)

  END IF
   
  ! SEND THIS DATA TO PARRALEL PROCESSORS
  SENDER = MSRID-1 ! SEND FROM MASTER
  CALL MPI_BCAST(MODE,1,MPI_INTEGER,SENDER,MPI_FVCOM_GROUP,IERR)

  CALL MPI_BCAST(UNITS,1,MPI_INTEGER,SENDER,MPI_FVCOM_GROUP,IERR)


  SELECT CASE(MODE) 
  CASE(BOX)
     
     CALL MPI_BCAST(BNDS_LR,2,MPI_F,SENDER,MPI_FVCOM_GROUP,IERR)
     CALL MPI_BCAST(BNDS_TB,2,MPI_F,SENDER,MPI_FVCOM_GROUP,IERR)

  CASE(RDS)

     CALL MPI_BCAST(CENTER,2,MPI_F,SENDER,MPI_FVCOM_GROUP,IERR)
     CALL MPI_BCAST(RADIUS,1,MPI_F,SENDER,MPI_FVCOM_GROUP,IERR)
     
  CASE DEFAULT
     CALL FATAL_ERROR&
          &("UNKNOWN VALUE IN SUBDOMAIN MODE - MPI ERROR!")

  END SELECT


  IF(DBG_SET(DBG_IO)) THEN
     WRITE(IPT,*) "! ==========================================="
     IF(MODE == box) THEN
        WRITE(IPT,*) "! SUBDOMAIN SETTINGS: MODE= BOX"
        WRITE(IPT,*) "! BOUNDS - LEFT/RIGHT =", BNDS_LR
        WRITE(IPT,*) "! BOUNDS - TOP/BOTTOM =", BNDS_TB
     ELSE
        WRITE(IPT,*) "! SUBDOMAIN SETTINGS: MODE= RADIUS"
        WRITE(IPT,*) "! CENTER LOCATION =", CENTER
        WRITE(IPT,*) "! RADIUS =", RADIUS
     END IF
     
     IF(UNITS == MTRS) THEN
        WRITE(IPT,*) "! UNITS = METERS"
     ELSE
        WRITE(IPT,*) "! UNITS = DEGREES"
     END IF
     WRITE(IPT,*)"! ==========================================="
  END IF



  ! DETERMINE THE LOCAL ELEMENT ID's WITHIN THE SUBDOMAIN

  ALLOCATE(EID(0:NT)); EID = 0

  SELECT CASE(UNITS)
  CASE(MTRS)

# if defined(SPHERICAL)
     IF(.not. USE_PROJ) CALL FATAL_ERROR&
          &("TO SPECIFY A SUBDOMAIN USING CARTESIAN PROJECTION",&
          & "YOU MUST COMPILE WITH PROJ WHEN USING THE SPHERICAL MODEL")
# endif

     SELECT CASE(MODE) 
     CASE(BOX)
        
        WHERE ( YMC < BNDS_TB(1) .and. YMC > BNDS_TB(2) .and. XMC < BNDS_LR(2) .and. XMC > BNDS_LR(1) )
           EID = 1
        END WHERE


     CASE(RDS)
     
        ALLOCATE(DIST(0:NT))
        DIST = (XMC-CENTER(1))**2 + (YMC - CENTER(2))**2

        WHERE (DIST < RADIUS**2)
           EID = 1
        END WHERE

        DEALLOCATE(DIST)
        
     END SELECT

  CASE(DGRS)

# if !defined(SPHERICAL)
     IF(.not. USE_PROJ) CALL FATAL_ERROR&
          &("TO SPECIFY A SUBDOMAIN USING SPHERICAL COORDINATES",&
          & "YOU MUST COMPILE WITH PROJ WHEN USING THE CARTESIAN MODEL")
# endif

     SELECT CASE(MODE) 
     CASE(BOX)

        ! IF BNDS_LR(1) < BNDS_LR(2) THE SELECTED AREA TAKES ONE ORIENTATION
        ! IF BNDS_LR(1) > BNDS_LR(2) THE SELECTED AREA IS FROM THE OTHER SIDE
        WHERE ( LATC < BNDS_TB(1) .and. LATC > BNDS_TB(2) .and. LONC < BNDS_LR(2) .and. LONC > BNDS_LR(1) )
           EID = 1
        END WHERE

     CASE(RDS)

        ALLOCATE(DIST(0:NT))
        DO I = 1, NT
           CALL ARC(CENTER(1),CENTER(2),LONC(I),LATC(I),DIST(I))
        END DO
        
        WHERE (DIST < RADIUS)
           EID = 1
        END WHERE

        DEALLOCATE(DIST)

     END SELECT

  END SELECT
     


  ALLOCATE(NID(0:MT)); NID = 0
  
  DO I = 1, NT
     
     IF(EID(I) == 1) THEN
        NID(NV(I,1:3)) = 1
     END IF
  END DO


END SUBROUTINE LOAD_SUBDOMAIN_FILE

# endif

!==============================================================================|
!
! THIS SUBOURTINE IS USED TO TRANSFER ALL THE DATA IN THE FVCOM GRID
  ! TO A LOCAL SUBDOMAIN. IT WILL CORRECTLY TAKE THE VALUES/INDEXS
  ! AND SHIFT THEM INTO THE SUBDOMAIN REFERENCE
!
!==============================================================================|
   SUBROUTINE SET_FVCOM_GRID(G)
     IMPLICIT NONE
     TYPE(GRID):: G
     INTEGER :: I,J,STATUS
     
    if(DBG_SET(dbg_sbr)) &
         & write(IPT,*) "START SET_FVCOM_GRID"

     ! SET UP GLOBAL GRID INDEX - USED TO MAKE INDEX ARRAYS FOR OUTPUT
     ! TO FILE

     G%NAME = "FVCOM"
     G%UNITS = "Not used yet"
     
     G%MGL  = MGL
     G%NGL  = NGL
     G%KB   = KB
     G%KBM1 = KBM1
     G%KBM2 = KBM2
#  if defined (WAVE_CURRENT_INTERACTION)
     G%MSC  = MSC
     G%MDC  = MDC
#  endif     
     
     NULLIFY(G%NV)

     NULLIFY(G%XM)
     NULLIFY(G%YM)
     NULLIFY(G%XMC)
     NULLIFY(G%YMC)

     NULLIFY(G%LAT)
     NULLIFY(G%LON)
     NULLIFY(G%LATC)
     NULLIFY(G%LONC)

     NULLIFY(G%NBE)
     NULLIFY(G%NBSN)
     NULLIFY(G%NBVE)
     
     IF(.NOT. IOPROC) THEN
        G%M  = M
        G%MT = MT
        G%N  = N
        G%NT = NT
        
        G%ELID   => ELID
        G%ELID_X => ELID_X
        G%NLID   => NLID
        G%NLID_X => NLID_X
        
        G%EGID   => EGID
        G%EGID_X => EGID_X
        G%NGID   => NGID
        G%NGID_X => NGID_X
        
        
        
        IF(PAR) THEN
           ! MUST ALLOCATE SPACE AND CREATE GLOBAL INDEX FOR OUTPUT

           ALLOCATE(G%NV(G%NT,3), STAT=STATUS)
           IF(STATUS /=0) CALL FATAL_ERROR("FVCOM2GRID: COULD NOT ALLOCATE G%NV")
           G%NV=0
           DO I = 1,G%NT
              G%NV(I,:) = G%NGID_X(NV(I,1:3))
           END DO
        

           ALLOCATE(G%NBE(0:N,3), STAT=STATUS)
           IF(STATUS /=0) CALL FATAL_ERROR("FVCOM2GRID: COULD NOT ALLOCATE G%NBE")
           G%NBE = 0
           DO I = 1,N
              G%NBE(I,:) = G%EGID_X(NBE(I,:))
           END DO


           ALLOCATE(G%NBSN(0:M,MX_NBR_ELEM+3), STAT=STATUS)
           IF(STATUS /=0) CALL FATAL_ERROR("FVCOM2GRID: COULD NOT ALLOCATE G%NBSN")
           G%NBSN = 0

           ALLOCATE(G%NBVE(0:M,MX_NBR_ELEM+1), STAT=STATUS)
           IF(STATUS /=0) CALL FATAL_ERROR("FVCOM2GRID: COULD NOT ALLOCATE G%NBVE")
           G%NBVE = 0

           DO I = 1,M
              G%NBVE(I,:) = G%EGID_X(NBVE(I,:))

              G%NBSN(I,:) = G%NGID_X(NBSN(I,:))
           END DO


        ELSE
           ! IF NOT PARALLEL JUST POINT TO FVCOM DATA FOR OUTPUT
           G%NV => NV
           
           G%NBE  => NBE
           G%NBVE => NBVE
           G%NBSN => NBSN

        END IF

     END IF
     
     if(DBG_SET(dbg_sbr)) &
          & write(IPT,*) "END SET_FVCOM_GRID"
    
   END SUBROUTINE SET_FVCOM_GRID
  

  SUBROUTINE SET_SUBDOMAIN_GRID(G)

    ! THIS CODE SETS THE VALUES FOR NBE, NV, NBVE, NBSN

    IMPLICIT NONE
    TYPE(GRID) :: G
    LOGICAL :: TEST
    INTEGER :: I,J,STATUS
    INTEGER, POINTER :: INV_E(:),INV_N(:)

    if(DBG_SET(dbg_sbr)) &
         & write(IPT,*) "START SET_SUBDOMAIN_GRID"


    IF(.NOT. IOPROC) THEN
       
       ! CREATE THE INVERSE ARRAY FOR NLID_X 
       ALLOCATE(INV_N(0:MT));INV_N=0
       DO I=1,G%MT
          INV_N(G%NLID_X(I))=I
       END DO


       ALLOCATE(INV_E(0:NT));INV_E=0
       DO I=1,G%NT
          INV_E(G%ELID_X(I))=I
       END DO

       
       ! MAKE NV - CONNECTIVITY
       ALLOCATE(G%NV(G%NT,3), STAT=STATUS)
       IF(STATUS /=0) CALL FATAL_ERROR("SET_SUBDOMAIN_GRID: COULD NOT ALLOCATE G%NV")
       
       ! CREAT THE GLOBAL INDEXED CONNECTIVITY ARRAY FOR THIS SUBDOMAIN
       DO I = 1,G%NT
          !       write(ipt,*) "NV(G%ELID_X(I),1:3)", NV(G%ELID_X(I),1:3)
          !       write(ipt,*) "INV(NV(G%ELID_X(I),1:3))" ,INV(NV(G%ELID_X(I),1:3))
          G%NV(I,:) = G%NGID_X(INV_N(NV(G%ELID_X(I),1:3)))
       END DO
       

       ! CREATE THE INVERSE ARRAY FOR NLID AND ELID


       ! MAKE NBE - ELEMENTS SURROUNDING EACH ELEMENT
       ALLOCATE(G%NBE(0:G%N,3), STAT=STATUS)
       IF(STATUS /=0) CALL FATAL_ERROR("FVCOM2GRID: COULD NOT ALLOCATE G%NBE")
       G%NBE = 0
       DO I = 1,G%N
          G%NBE(I,:) = G%EGID_X(INV_E(NBE(G%ELID(I),:)))
       END DO
       

       ! MAKE NBSN AND NBVE
       ALLOCATE(G%NBSN(0:G%M,MX_NBR_ELEM+3), STAT=STATUS)
       IF(STATUS /=0) CALL FATAL_ERROR("FVCOM2GRID: COULD NOT ALLOCATE G%NBSN")
       G%NBSN = 0
       
       ALLOCATE(G%NBVE(0:G%M,MX_NBR_ELEM+1), STAT=STATUS)
       IF(STATUS /=0) CALL FATAL_ERROR("FVCOM2GRID: COULD NOT ALLOCATE G%NBVE")
       G%NBVE = 0
       
       DO I = 1,G%M
          G%NBVE(I,:) = G%EGID_X(INV_E(NBVE(G%NLID(I),:)))
          
          G%NBSN(I,:) = G%NGID_X(INV_N(NBSN(G%NLID(I),:)))
       END DO
       
       
       DEALLOCATE(INV_N)
       DEALLOCATE(INV_E)
    

    END IF


    if(DBG_SET(dbg_sbr)) &
         & write(IPT,*) "END SET_SUBDOMAIN_GRID"

  END SUBROUTINE SET_SUBDOMAIN_GRID

!=============================================================  
!=============================================================  
!=============================================================  
!=============================================================      
  SUBROUTINE SETUP_MPI_IO_MODE(TF,COMMGRP) 
    !===================================================================================|
    !  INITIALIZE MPI ENVIRONMENT                                                       |
    !===================================================================================|
    LOGICAL, INTENT(INOUT) :: TF
    INTEGER, INTENT(OUT)   :: COMMGRP
    INTEGER :: total_group  
    INTEGER :: fvcom_group
    INTEGER :: SBUF,RBUF, trueval, i
    INTEGER, allocatable :: fvcom_subset(:)
    INTEGER IERR


    if(DBG_SET(dbg_sbr)) &
         & write(IPT,*) "STARTING SETUP_MPI_IO_MODE"


    if (NPROCS .LE. 3) THEN
       if (TF) &
            & CALL WARNING("FVCOM CAN NOT USE MPI IO MODE WHEN RUN &
            &ON LESS THAN 4 PROCESSORS", "CONTINUING WITHOUT THIS OPTION!")
       TF = .false.
    end if

    IF (TF) THEN

# if defined (MULTIPROCESSOR)
       ! MODIFY CONTROL VARIABLE EFFECTED BY THE USE OF MPI IO MODE
       NULLIFY(NPROCS_TOTAL)
       ALLOCATE(NPROCS_TOTAL)
       NPROCS_TOTAL = NPROCS

       IOPROCID = NPROCS

       IOPROC=.FALSE.
       IF(MYID==NPROCS) IOPROC=.TRUE.

       ! NOW REDUCE THE NUMBER OF PROCESSORS BY ONE
       NPROCS= NPROCS - 1

       ! RETURN HANDLE TO LIST OF GROUPS
       call mpi_comm_group(mpi_comm_world,total_group,ierr)

       !first: define subset
       !comp procs will have process_id/myid:  0-->(nprocs-1)/1-->nprocs
       !mr_printy  will have process_id/myid:  nprocs/nprocs_tot
       allocate(fvcom_subset(nprocs)) ; fvcom_subset = 0
       do i=1,nprocs
          fvcom_subset(i) = i-1
       end do
       call mpi_group_incl(total_group,nprocs,fvcom_subset,fvcom_group,ierr)
       call mpi_comm_create(mpi_comm_world,fvcom_group,COMMGRP,ierr)
       deallocate(fvcom_subset)


       call mpi_barrier(mpi_comm_world,ierr)

       IF(.not. IOPROC) then
          sbuf = myid
          rbuf = 0
          call mpi_allreduce(sbuf,rbuf,1,mpi_integer,mpi_sum,COMMGRP,ierr)

          trueval = 0
          do i=1,nprocs
             trueval = trueval + i
          end do
          
          if( trueval .NE. rbuf) &
               & CALL FATAL_ERROR("TESTING GROUP COMMUNICATION FOR MPI &
               &IO MODE FAILED")
       END IF
       
       
       sbuf = myid 
       rbuf = 0
       call mpi_allreduce(sbuf,rbuf,1,mpi_integer,mpi_sum,MPI_COMM_WORLD,ierr)

       trueval = 0
       do i=1,nprocs_total
          trueval = trueval + i
       end do

       if( trueval .NE. rbuf ) &
            & CALL FATAL_ERROR("TESTING WORLD COMMUNICATION FOR MPI &
            &IO MODE FAILED")

       IF (DBG_SET(dbg_log)) &
            & write(IPT,*) "!  MPI IO MODE IS ACTIVE"
# endif

    ELSE 

       ! MUST SET INTENT(OUT) COMMGRP
# if defined (MULTIPROCESSOR)
       ! IF MULTIPROCESSOR BUT NOT MPI IO SET IT TO COMM WORLD!
       COMMGRP = MPI_COMM_WORLD
# else
       ! IF NOT MULTIPROCESSOR SET A DUMMY VALUE
       COMMGRP = -999

# endif


       IF (DBG_SET(dbg_log)) &
            & write(IPT,*) "!  MPI IO MODE IS NOT ACTIVE"

    END IF

    if(DBG_SET(dbg_sbr)) &
         & write(IPT,*) "END SETUP_MPI_IO_MODE"
  END SUBROUTINE SETUP_MPI_IO_MODE


# if defined (MULTIPROCESSOR)

  SUBROUTINE MPI_IO_LOOP
    implicit none
    integer RBUF, J, IERR
    integer SOURCE
    integer status
    INTEGER STAT(MPI_STATUS_SIZE)
    character(len=5) :: lpc

    if(DBG_SET(dbg_sbr)) &
         & write(IPT,*) "START MPI_IO_LOOP"

    IN_MPI_IO_LOOP = .TRUE.

    if (MYID .NE. IOPROCID) then
       IF(DBG_SET(DBG_LOG))THEN
          WRITE(IPT,*) "!++++++++++++++++++++++++++++++++++"
          WRITE(IPT,*) "! I AM NOT THE IO PROC:"
          WRITE(IPT,*) "! SETTING IN_MPI_IO_LOOP = T"
          WRITE(IPT,*) "!++++++++++++++++++++++++++++++++++"
          if(DBG_SET(dbg_sbr)) &
               & write(IPT,*) "END MPI_IO_LOOP"
       END IF
       
       RETURN

    else   
       IF(DBG_SET(DBG_LOG))THEN
          WRITE(IPT,*) "!++++++++++++++++++++++++++++++++++++++++++++"
          WRITE(IPT,*) "! I AM THE IO PROC: STARTING IO LOOP"
          WRITE(IPT,*) "!++++++++++++++++++++++++++++++++++++++++++++"
       END IF
    end if


    J = 1
    DO 
       write(lpc,'(I5.5)') J
       IF(DBG_SET(DBG_IO)) write(IPT,*)"========IO LOOP COUNT "//lpc//"==============="

       !      if (J==5) then
       !         call system("sleep 5")
       !         call report_error("testing!")
       !      end if

       IF(DBG_SET(DBG_IO)) WRITE(IPT,*) "WAITING FOR IO CODE:"

       STAT=0
!       SOURCE = MPI_ANY_SOURCE ! ALLOW ANY PROCESSOR TO COMMUNICATE WITH THE IO PROC
       SOURCE = 0 ! ONLY THE MASTER CAN COMMMUNCATE WITH IO PROC
       CALL MPI_RECV(RBUF,1,MPI_INTEGER,SOURCE&
            &,SYNC_TAG,MPI_COMM_WORLD,STAT,IERR)

       select case(RBUF)
       case(EXT_CODE)
          ! JUST CALL STOP...
          IF(DBG_SET(DBG_IO)) WRITE(IPT,*) "IO PROC RECIEVED EXIT CODE:"
          IF(DBG_SET(DBG_IO)) WRITE(IPT,*) "THATS ALL FOLKS!"

          Call mpi_finalize(IERR)
          stop

       case(WAIT_CODE)

          IF(DBG_SET(DBG_IO)) write(IPT,*) "OTHER PROCS WAITING FOR ME:"

          Call mpi_barrier(mpi_comm_world,ierr)

          IF(DBG_SET(DBG_IO)) write(IPT,*) "FINISHED WRITING LAST FILE: GO"

       case default

          IF(DBG_SET(DBG_IO)) write(IPT,*) "IO PROC GOT GO FOR CALL:",RBUF

          Call mpi_barrier(mpi_comm_world,ierr)

          IF(DBG_SET(DBG_IO)) write(IPT,*) "IO PROC SINKED FOR GO"

          call CALL_FUNC(RBUF,status)
          IF (status/=0) call fatal_error("IO PROC Recieved bad go code",&
               & "Could not retrieve valid function pointer?")

          IF(DBG_SET(DBG_IO)) write(IPT,*) "FINISHED MPI_IO CALL!"

       end select

       J = J + 1

    END DO

  END SUBROUTINE MPI_IO_LOOP



  SUBROUTINE MPI_IO_SYNCHRONIZE(CODE)
    implicit none
    INTEGER, INTENT(IN) ::  CODE
    integer :: ierr, STATUS
    INTEGER STAT(MPI_STATUS_SIZE)

    if(DBG_SET(dbg_sbr)) &
         & write(IPT,*) "START MPI_IO_SYNCHRONIZE"



    IF (MYID .EQ. IOPROCID) CALL FATAL_ERROR("IO PROC SHOULD NEVER&
         & BE IN MPI_IO_SYNCHRONIZE")


    IF(DBG_SET(DBG_IO)) write(IPT,*) "SYNCHRONIZE FVCOM_GROUP WITH IOPROC:"

    select CASE(CODE)
    case(WAIT_CODE)

       if (MYID==1) then
          IF(DBG_SET(DBG_IO)) write(IPT,*) "MASTER SENDING: WAIT_CODE"
          CALL MPI_SEND(CODE,1,MPI_INTEGER,IOPROCID-1,SYNC_TAG&
               &,MPI_COMM_WORLD,IERR)
       else
          IF(DBG_SET(DBG_IO)) write(IPT,*) "WAITING FOR IO PROC"
       end if

       call mpi_barrier(mpi_comm_world,ierr)

       IF(DBG_SET(DBG_IO)) write(IPT,*) "IO PROC HAS FINISHED LAS&
            &T WRITE COMMAND"

    case default

       if (MYID==1) then
          IF(DBG_SET(DBG_IO)) write(IPT,*) "MASTER SENDING GO CODE:",CODE
          CALL MPI_SEND(CODE,1,MPI_INTEGER,IOPROCID-1,SYNC_TAG,MPI_COMM_WORLD,IERR)
       else
          IF(DBG_SET(DBG_IO)) write(IPT,*) "WAITING FOR IO PROC"
       end if

       call mpi_barrier(mpi_comm_world,ierr)

          call CALL_FUNC(CODE,status)
          IF (status/=0) call fatal_error("MPI_IO_SYNCHRONIZE: Bad go code",&
               & "Could not retrieve valid function pointer?")

          IF(DBG_SET(DBG_IO)) write(IPT,*) "FINISHED MPI_IO CALL!"

    end select


  END SUBROUTINE MPI_IO_SYNCHRONIZE

# else
    SUBROUTINE MPI_IO_SYNCHRONIZE(CODE)
    implicit none
    integer code

    CALL FATAL_ERROR("THIS IS A DUMMY ROUTINE:","MPI_IO_SYNCHRONIZ&
         &E only exists when FVCOM is compiled with multiprocessor",&
         & "USE_MPI_IO_MODE should not be .true.")
  END SUBROUTINE MPI_IO_SYNCHRONIZE

# endif


# if defined (MEAN_FLOW)
    SUBROUTINE MEANFLOW_RST
    USE MOD_MEANFLOW
    USE MOD_OBCS2
    USE MOD_OBCS3
    IMPLICIT NONE
    INTEGER :: I,K
    CHARACTER(LEN=4) :: FILENUMBER

    IF(nmfcell_i>0)THEN
      MF_RST_STCNT = MF_RST_STCNT + 1
      WRITE(FILENUMBER,'(I4.4)')MF_RST_STCNT
      OPEN(1,FILE=TRIM(OUTPUT_DIR)//'re_mf'//FILENUMBER//'.dat',FORM='UNFORMATTED'&
           &,STATUS='UNKNOWN')      

      WRITE(1)IINT
      WRITE(1)(ELT   (I),   I=0,ntidenode)
      WRITE(1)(ELTF  (I),   I=0,ntidenode)
      WRITE(1)(ELRKT (I),   I=0,ntidenode)
      WRITE(1)(ELTDT (I),   I=0,ntidenode)
      WRITE(1)(ELP   (I),   I=0,ntidenode)
      WRITE(1)(ELPF  (I),   I=0,ntidenode)
      WRITE(1)(ELRKP (I),   I=0,ntidenode)
      WRITE(1)(UAT   (I),   I=0,ntidecell)
      WRITE(1)(VAT   (I),   I=0,ntidecell)
      WRITE(1)(UATF  (I),   I=0,ntidecell)
      WRITE(1)(VATF  (I),   I=0,ntidecell)
      WRITE(1)(UAP   (I),   I=0,ntidecell)
      WRITE(1)(VAP   (I),   I=0,ntidecell)
      WRITE(1)(UANT  (I),   I=0,  nmfcell)
      WRITE(1)(VANT  (I),   I=0,  nmfcell)
      WRITE(1)(UAN   (I),   I=0,  nmfcell)
      WRITE(1)(VAN   (I),   I=0,  nmfcell)
      WRITE(1)(UANP  (I),   I=0,  nmfcell)
      WRITE(1)(VANP  (I),   I=0,  nmfcell)

      WRITE(1)((UT  (I,K), I=0,ntidecell),K=1,KBM1)
      WRITE(1)((VT  (I,K), I=0,ntidecell),K=1,KBM1)
      WRITE(1)((UNT (I,K), I=0,  nmfcell),K=1,KBM1)
      WRITE(1)((VNT (I,K), I=0,  nmfcell),K=1,KBM1)
      WRITE(1)((UN  (I,K), I=0,  nmfcell),K=1,KBM1)
      WRITE(1)((VN  (I,K), I=0,  nmfcell),K=1,KBM1)

      WRITE(1)(UAPF  (I),  I=0,ntidecell)
      WRITE(1)(VAPF  (I),  I=0,ntidecell)
      WRITE(1)(UANTF (I),  I=0,  nmfcell)
      WRITE(1)(VANTF (I),  I=0,  nmfcell)
      WRITE(1)(UANF  (I),  I=0,  nmfcell)
      WRITE(1)(VANF  (I),  I=0,  nmfcell)
      WRITE(1)(UANPF (I),  I=0,  nmfcell)
      WRITE(1)(VANPF (I),  I=0,  nmfcell)
      WRITE(1)(UARKNT(I),  I=0,  nmfcell)
      WRITE(1)(VARKNT(I),  I=0,  nmfcell)
      WRITE(1)(UARKN (I),  I=0,  nmfcell)
      WRITE(1)(VARKN (I),  I=0,  nmfcell)

      WRITE(1)((UNTB (I,K), I=0,nmfcell),K=1,KBM1)
      WRITE(1)((VNTB (I,K), I=0,nmfcell),K=1,KBM1)
      WRITE(1)((UNB  (I,K), I=0,nmfcell),K=1,KBM1)
      WRITE(1)((VNB  (I,K), I=0,nmfcell),K=1,KBM1)
      
      WRITE(1)(FLUXOBN2(I),I=0,IBCN(2))
      WRITE(1)(CXOBC(I),I=0,NT)
      WRITE(1)(CYOBC(I),I=0,NT)
      WRITE(1)(FLUXOBC2D_X    (I)  , I=0,nmfcell_i) 
      WRITE(1)(FLUXOBC2D_Y    (I)  , I=0,nmfcell_i)  
      WRITE(1)((FLUXOBC3D_X   (I,K), I=0,nmfcell_i),K=1,KBM1)
      WRITE(1)((FLUXOBC3D_Y   (I,K), I=0,nmfcell_i),K=1,KBM1)    
      WRITE(1)((FLUXOBC3D_X_2 (I,K), I=0,nmfcell_i),K=1,KBM1)
      WRITE(1)((FLUXOBC3D_Y_2 (I,K), I=0,nmfcell_i),K=1,KBM1)
      WRITE(1)(OBC2D_X_TIDE   (I)  , I=1,  nmfcell)
      WRITE(1)(OBC2D_Y_TIDE   (I)  , I=1,  nmfcell)

      CLOSE(1)
      WRITE(IPT,*)'!  MEANFLOW DATA         :    DONE     '
    ENDIF

    END SUBROUTINE MEANFLOW_RST
# endif

!=============================================================  
  FUNCTION VELOCITY_FILE_OBJECT_SURFACE() RESULT(NCF)
    IMPLICIT NONE

    INTEGER :: status
    LOGICAL, SAVE :: IOPROC_ALLOCATED = .FALSE.
    TYPE(NCFILE), POINTER :: NCF
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCATT),  POINTER :: ATT

    character(len=100)    :: timestamp, temp, netcdf_convention

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START VELOCITY_FILE_OBJECT_SURFACE"
    
    ! IO PROC MUST ALLOCATE SPACE FOR THE ARRAYS 
    ! THESE ARRAYS MUST HAVE THE ATTRIBUTE SAVE AND FOR CLARITY
    ! SHOULD HAVE THE SAME NAME AS THOSE USED ON THE OTHER PROCESSORS!
    IF(IOPROC .AND. .NOT. IOPROC_ALLOCATED) THEN

       IOPROC_ALLOCATED = .true.

       allocate(U_SURFACE(NGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:U")
       U_SURFACE = 0.0_sp

       allocate(V_SURFACE(NGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:V")
       V_SURFACE = 0.0_sp

       allocate(TAUBM(NGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:WUBOT")
       TAUBM   = 0.0_SP

    END IF


  ! ALLOCATE THE NEW FILE OBJECT
    NCF => NEW_FILE()


    ! U_SURFACE
    VAR  => NC_MAKE_AVAR(name='u_surface',&
         & values=U_SURFACE, DIM1= DIM_nele, DIM2 = DIM_time)

# if defined(UCF)
    ATT  => NC_MAKE_ATT(name='long_name',values='East Surface Velocity Component') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='standard_name',values='eastward_seasurface_water_velocity') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='m s-1') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='grid5') 
    VAR  => ADD(VAR,ATT)

!!$    ATT  => NC_MAKE_ATT(name='grid_location',values='center') 
    ATT  => NC_MAKE_ATT(name='location',values='face') 
    VAR  => ADD(VAR,ATT)

!!$    ATT  => NC_MAKE_ATT(name='coordinates',values='latc lonc siglay_center') 
    ATT  => NC_MAKE_ATT(name='coordinates',values='time latc lonc') 
    VAR  => ADD(VAR,ATT)

# else
    ATT  => NC_MAKE_ATT(name='long_name',values='Eastward Surface Water Velocity') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='standard_name',values='eastward_seasurface_water_velocity') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='meters s-1') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates',values='time latc lonc') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='mesh',values='fvcom_mesh') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='location',values='face') 
    VAR  => ADD(VAR,ATT)
# endif

    NCF  => ADD(NCF,VAR)



    ! V_SURFACE
    VAR  => NC_MAKE_AVAR(name='v_surface',&
         & values=V_SURFACE, DIM1= DIM_nele, DIM2 = DIM_time)

# if defined(UCF)
    ATT  => NC_MAKE_ATT(name='long_name',values='North Surface Velocity Component') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='standard_name',values='Northward_seasurface_water_velocity') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='m s-1') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='grid5') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid_location',values='center') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates',values='latc lonc') 
    VAR  => ADD(VAR,ATT)

# else
    ATT  => NC_MAKE_ATT(name='long_name',values='Northward Surface Water Velocity') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='standard_name',values='Northward_seasurface_water_velocity') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='meters s-1') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates',values='time latc lonc') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='mesh',values='fvcom_mesh') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='location',values='face') 
    VAR  => ADD(VAR,ATT)
# endif

    NCF  => ADD(NCF,VAR)

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END VELOCITY_FILE_OBJECT_SURFACE"

  END FUNCTION VELOCITY_FILE_OBJECT_SURFACE
!=============================================================  
  FUNCTION SALT_TEMP_FILE_OBJECT_SURFACE() RESULT(NCF)
    IMPLICIT NONE

    INTEGER :: status
    LOGICAL, SAVE :: IOPROC_ALLOCATED = .FALSE.
    TYPE(NCFILE), POINTER :: NCF
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCATT),  POINTER :: ATT

    character(len=100)    :: timestamp, temp, netcdf_convention

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START SALT_TEMP_FILE_OBJECT_SURFACE"
    
    ! IO PROC MUST ALLOCATE SPACE FOR THE ARRAYS 
    ! THESE ARRAYS MUST HAVE THE ATTRIBUTE SAVE AND FOR CLARITY
    ! SHOULD HAVE THE SAME NAME AS THOSE USED ON THE OTHER PROCESSORS!
    IF(IOPROC .AND. .NOT. IOPROC_ALLOCATED) THEN

       IOPROC_ALLOCATED = .true.

       allocate(T1_SURFACE(MGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:T1_SURFACE")
       T1_SURFACE = 0.0_SP

       allocate(S1_SURFACE(MGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:S1_SURFACE")
       S1_SURFACE = 0.0_SP

    END IF


  ! ALLOCATE THE NEW FILE OBJECT
    NCF => NEW_FILE()


    ! T_SURFACE
    VAR  => NC_MAKE_AVAR(name='temp_surface',&
         & values=T1_SURFACE, DIM1= DIM_node, DIM2 = DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='surface temperature') 
    VAR  => ADD(VAR,ATT)

# if defined(UCF)
    ATT  => NC_MAKE_ATT(name='standard_name',values='seasurface_water_potential_temperature') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='degC') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='grid4') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates',values='lon lat') 
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='grid_location',values='node') 
    VAR  => ADD(VAR,ATT)    

# else
    ATT  => NC_MAKE_ATT(name='standard_name',values='seasurface_water_temperature') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='degrees_C') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

!!$    ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar) 
    ATT  => NC_MAKE_ATT(name='coordinates',values='time lat lon') 
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)    

    ATT  => NC_MAKE_ATT(name='mesh',values='fvcom_mesh') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='location',values='node') 
    VAR  => ADD(VAR,ATT)
# endif

    NCF  => ADD(NCF,VAR)


    ! S_SURFACE
    VAR  => NC_MAKE_AVAR(name='salinity_surface',&
         & values=S1_SURFACE, DIM1= DIM_node, DIM2 = DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='surface salinity') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='standard_name',values='seasurface_water_salinity') 
    VAR  => ADD(VAR,ATT)

# if defined(UCF)
    ATT  => NC_MAKE_ATT(name='units',values='1e-3') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='grid4') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates',values='lon lat') 
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='grid_location',values='node') 
    VAR  => ADD(VAR,ATT)    

# else

    ATT  => NC_MAKE_ATT(name='units',values='1e-3') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

!!$    ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar) 
    ATT  => NC_MAKE_ATT(name='coordinates',values='time lat lon') 
    VAR  => ADD(VAR,ATT)    

    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='mesh',values='fvcom_mesh') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='location',values='node') 
    VAR  => ADD(VAR,ATT)
# endif

    NCF  => ADD(NCF,VAR)


    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END SALT_TEMP_FILE_OBJECT_SURFACE"

  END FUNCTION SALT_TEMP_FILE_OBJECT_SURFACE
!=============================================================  
  FUNCTION TURBULENCE_FILE_OBJECT_SURFACE() RESULT(NCF)
    IMPLICIT NONE

    INTEGER :: status
    LOGICAL, SAVE :: IOPROC_ALLOCATED = .FALSE.
    TYPE(NCFILE), POINTER :: NCF
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCATT),  POINTER :: ATT
    TYPE(NCDIM),  POINTER :: DIM1
    TYPE(NCDIM),  POINTER :: DIM2

    character(len=100)    :: timestamp, temp, netcdf_convention

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START TURBULENCE_FILE_OBJECT_SURFACE"
    
    ! IO PROC MUST ALLOCATE SPACE FOR THE ARRAYS 
    ! THESE ARRAYS MUST HAVE THE ATTRIBUTE SAVE AND FOR CLARITY
    ! SHOULD HAVE THE SAME NAME AS THOSE USED ON THE OTHER PROCESSORS!
    IF(IOPROC .AND. .NOT. IOPROC_ALLOCATED) THEN

       IOPROC_ALLOCATED = .true.

       allocate(VISCOFM_SURFACE(NGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:VISCOFM_SURFACE")
       VISCOFM_SURFACE = 0.0_SP

       allocate(VISCOFH_SURFACE(MGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:VISCOFH_SURFACE")
       VISCOFH_SURFACE = 0.0_SP
       
    END IF


  ! ALLOCATE THE NEW FILE OBJECT
    NCF => NEW_FILE()


    ! VISCOFM_SURFACE
    VAR  => NC_MAKE_AVAR(name='viscofm_surface',&
         & values=viscofm_surface, DIM1= DIM_nele, DIM2 = DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Surface Horizontal Turbulent Eddy Viscosity For Momentum') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='m 2 s-1') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar) 
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! VISCOFH_SURFACE
    VAR  => NC_MAKE_AVAR(name='viscofh_surface',&
         & values=viscofh_surface, DIM1= DIM_node, DIM2 = DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Surface Horizontal Turbulent Eddy Viscosity For Scalars') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='m 2 s-1') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar) 
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END TURBULENCE_FILE_OBJECT_SURFACE"

  END FUNCTION TURBULENCE_FILE_OBJECT_SURFACE
!=============================================================  


END MODULE MOD_NCDIO
