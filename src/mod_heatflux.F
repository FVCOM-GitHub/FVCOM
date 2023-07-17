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

!==============================================================================|
!   CONTROL VARIABLES                                                          |
!==============================================================================|

# if defined (HEATING_CALCULATED)
MODULE MOD_HEATFLUX

   USE ALL_VARS
   USE MOD_UTILS
   USE MOD_PREC

   IMPLICIT NONE
   SAVE

   LOGICAL HEATING_CALCULATE_ON
   LOGICAL HEATING_FRESHWATER
   CHARACTER(LEN=80) HEATING_CALCULATE_TYPE
   CHARACTER(LEN=80) HEATING_CALCULATE_FILE
   CHARACTER(LEN=80) HEATING_CALCULATE_KIND
   CHARACTER(LEN=80) COARE_VERSION
!   REAL(SP) :: ZUU                            !!HEIGHT OF WIND SPEED (M) ! Siqi Li, 2021-01-27
   REAL(SP) :: ZTT                            !!HEIGHT OF AIR TEMPERATURE (M)
   REAL(SP) :: ZQQ                            !!HEIGHT OF RELATIVE HUMIDITY (M)
   REAL(SP) :: AIR_TEMPERATURE
   REAL(SP) :: RELATIVE_HUMIDITY
   REAL(SP) :: SURFACE_PRESSURE
   REAL(SP) :: LONGWAVE_RADIATION
   REAL(SP) :: SHORTWAVE_RADIATION
   REAL(SP) :: HEATING_LONGWAVE_PERCTAGE_IN_HEATFLUX
   REAL(SP) :: HEATING_LONGWAVE_LENGTHSCALE_IN_HEATFLUX
   REAL(SP) :: HEATING_SHORTWAVE_LENGTHSCALE_IN_HEATFLUX
   
   NAMELIST /NML_HEATING_CALCULATED/           &
        & HEATING_CALCULATE_ON,                &
	& HEATING_CALCULATE_TYPE,              &
	& HEATING_CALCULATE_FILE,              &
	& HEATING_CALCULATE_KIND,              &
        & HEATING_FRESHWATER,                  &
	& COARE_VERSION,                       &
!	& ZUU,                                 & ! Siqi Li, 2021-01-27
	& ZTT,                                 &
	& ZQQ,                                 &
	& AIR_TEMPERATURE,                     &
	& RELATIVE_HUMIDITY,                   &
	& SURFACE_PRESSURE,                    &
	& LONGWAVE_RADIATION,                  &
	& SHORTWAVE_RADIATION,                 &
        & HEATING_LONGWAVE_PERCTAGE_IN_HEATFLUX,    &
        & HEATING_LONGWAVE_LENGTHSCALE_IN_HEATFLUX, &
        & HEATING_SHORTWAVE_LENGTHSCALE_IN_HEATFLUX
   
   REAL(SP), ALLOCATABLE :: CORRG(:),CORR(:)           !!LATITUDE OF NODES
   REAL(SP), ALLOCATABLE :: RAIN(:),CP40(:),SIGH(:),ZI40(:)
   REAL(SP), ALLOCATABLE :: USRCOARE(:)              !MDR 3-19-2013 for wind wave mix
   REAL(SP), PARAMETER :: emmiss = 0.97_SP           ! ejw 8/16/2006 based on ROMS
   REAL(SP), PARAMETER :: StefBo = 5.67e-8_SP        ! ejw 8/16/2006 based on ROMS
   REAL(SP), PARAMETER :: fmiss  = -9999._SP

!==============================================================================|


   CONTAINS !---------------------------------------------------------------------|
            ! HEATING_CALCULATE_NAMELIST_INITIALIZE : Initialize the values in    |
	    !                                     namelist NML_HEATING_CALCULATED |
            ! HEATING_CALCULATE_NAMELIST_PRINT      : Print the values of namelist|
	    !                                         NML_HEATING_CALCULATED      |
            ! HEATING_CALCULATE_NAMELIST_READ       : Read the values of namelist |
	    !                                         NML_HEATING_CALCULATED      |
	    !---------------------------------------------------------------------|
	    
!==============================================================================|
!   Input Parameters Which Control the Calculation of Heat Flux                |
!==============================================================================|

   SUBROUTINE HEATING_CALCULATE_NAMELIST_INITIALIZE
   USE control, only : casename
   IMPLICIT NONE
   
   HEATING_CALCULATE_ON   = .FALSE.                
   HEATING_FRESHWATER     = .FALSE.                
   HEATING_CALCULATE_TYPE = "'flux' or 'body'"              
   HEATING_CALCULATE_FILE = trim(casename)//"_hfx.nc"              
   HEATING_CALCULATE_KIND = "Options:"//TRIM(CNSTNT)//","//TRIM(STTC)//","//TRIM(TMDPNDNT)//","//TRIM(PRDC)//","//TRIM(VRBL)              
   COARE_VERSION          = "'COARE26Z' or 'COARE40VN'"
!   ZUU                    = 2  ! Unit = m         ! Siqi Li, 2021-01-27                        
   ZTT                    = 2  ! Unit = m                                 
   ZQQ                    = 2  ! Unit = m                                 
   AIR_TEMPERATURE        = 0.0_SP                     
   RELATIVE_HUMIDITY      = 0.0_SP                   
   SURFACE_PRESSURE       = 0.0_SP                    
   LONGWAVE_RADIATION     = 0.0_SP                  
   SHORTWAVE_RADIATION    = 0.0_SP
   HEATING_LONGWAVE_PERCTAGE_IN_HEATFLUX = 0.78_SP
   HEATING_LONGWAVE_LENGTHSCALE_IN_HEATFLUX = 1.4_SP
   HEATING_SHORTWAVE_LENGTHSCALE_IN_HEATFLUX= 6.3_SP
   
   RETURN
   END SUBROUTINE HEATING_CALCULATE_NAMELIST_INITIALIZE
   
!------------------------------------------------------------------------------|
   SUBROUTINE HEATING_CALCULATE_NAMELIST_PRINT
   USE CONTROL, ONLY : IPT
   
   IMPLICIT NONE
   
   WRITE(UNIT=IPT,NML=NML_HEATING_CALCULATED)
   
   RETURN
   END SUBROUTINE HEATING_CALCULATE_NAMELIST_PRINT  
   
!------------------------------------------------------------------------------|
   SUBROUTINE HEATING_CALCULATE_NAMELIST_READ    
   USE CONTROL, ONLY : casename,NMLUNIT
   
   IMPLICIT NONE
   
   INTEGER :: IOS, I
   CHARACTER(LEN=120) :: FNAME
   CHARACTER(LEN=160) :: PATHNFILE
   
   IF(DBG_SET(DBG_SBR)) &
         & WRITE(IPT,*) "Subroutine Begins: Read_Heating_Calculate_Namelist;"

   IOS = 0

   FNAME = "./"//trim(casename)//"_run.nml"

   CALL FOPEN(NMLUNIT,trim(FNAME),'cfr')

   !READ NAME LIST FILE
    REWIND(NMLUNIT)

   ! Read IO Information
   READ(UNIT=NMLUNIT, NML=NML_HEATING_CALCULATED,IOSTAT=ios)
   if(ios .NE. 0 ) Then
     if(DBG_SET(dbg_log)) write(UNIT=IPT,NML=NML_HEATING_CALCULATED)
     Call Fatal_error("Can Not Read NameList NML_HEATING_CALCULATED from file: "//trim(FNAME))
   end if
   CLOSE(NMLUNIT)
   
   IF(HEATING_CALCULATE_ON .AND. .NOT. HEATING_ON)THEN
     HEATING_ON = .TRUE.
   ELSE IF(HEATING_CALCULATE_ON .AND. HEATING_ON)THEN
     CALL FATAL_ERROR("IN NAMELIST CONTROL FILE, IF HEATING_CALCULATE_ON = TRUE, ", &
                      "THEN SET HEATING_ON = FALSE") 
   END IF 
   IF(HEATING_CALCULATE_TYPE == 'flux') HEATING_TYPE = 'flux'
   IF(HEATING_CALCULATE_TYPE == 'body') HEATING_TYPE = 'body'

   RETURN
   END SUBROUTINE HEATING_CALCULATE_NAMELIST_READ

!==============================================================================|

END MODULE MOD_HEATFLUX
# endif


