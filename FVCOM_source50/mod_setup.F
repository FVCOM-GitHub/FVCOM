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

MODULE MOD_SETUP
  USE ALL_VARS
  USE MOD_PAR
  USE MOD_INPUT
  USE MOD_NCDIO
  USE MOD_NCTOOLS
  USE MOD_UTILS
  USE MOD_OBCS
  USE MOD_FORCE

  IMPLICIT NONE

  SAVE
  ! COORDINATE VARIABLES
  REAL(SP), ALLOCATABLE, TARGET :: X_GBL(:),Y_GBL(:)
  REAL(SP), ALLOCATABLE, TARGET :: X_LCL(:),Y_LCL(:)
  
  ! CORIOLIS VARIABLES
  REAL(SP), ALLOCATABLE, TARGET :: C_LCL(:)     !!CORIOLIS PARAMETER AT NODES

  ! DEPTH VARIABLE
  REAL(SP), ALLOCATABLE, TARGET :: H_LCL(:)     !! DEPTH PARAMETER AT NODES

  ! SPONGE LAYER VARIABLE
  INTEGER, ALLOCATABLE, TARGET  :: N_SPG(:)
  REAL(SP), ALLOCATABLE, TARGET :: R_SPG(:),C_SPG(:),X_SPG(:),Y_SPG(:)
  INTEGER :: NSPONGE

  PRIVATE :: SIGMA_GEOMETRIC, SIGMA_GENERALIZED, SIGMA_TANH

CONTAINS

!==============================================================================!
  SUBROUTINE SETUP_CENTER_COORDS
    USE MOD_SPHERICAL
    IMPLICIT NONE
    INTEGER I,J,IERR,STATUS,SENDID
    REAL(SP) SBUF

    INTEGER K,ITMP
    REAL(DP) VX1,VY1,VX2,VY2,VX3,VY3
    REAL(DP) EVX12,EVX13,EVX23,EVY12,EVY13,EVY23,EVXY
    REAL(DP) VX12,VY12,VX23,VY23,VX31,VY31
    INTEGER :: SENDER

    ! David Added:
    REAL(SP), allocatable :: xc_buf(:), lonc_buf(:)
    REAL(SP), allocatable :: yc_buf(:), latc_buf(:)
 
    IF(DBG_SET(DBG_SBR)) write(IPT,*) "SETUP_CENTER_COORDS: START"
    
    
!==============================================================================|
!   SET UP LOCAL MESH (HORIZONTAL COORDINATES)                                 |
!==============================================================================|
!--------------CALCULATE GLOBAL MINIMUMS AND MAXIMUMS--------------------------!
    
#  if !defined (SPHERICAL) 
! IF THIS IS NOT A SPHERICAL MODEL REMOVE MIN COORDINATE VALUE    
    
    VXMIN = MINVAL(VX(1:MT)) ; VXMAX = MAXVAL(VX(1:MT))
    VYMIN = MINVAL(VY(1:MT)) ; VYMAX = MAXVAL(VY(1:MT))
    
#  if defined (MULTIPROCESSOR) 
    IF(PAR) THEN
       !GLOBAL MIN
       SBUF = VXMIN
       CALL MPI_ALLREDUCE(SBUF,VXMIN,1,MPI_F,MPI_MIN,MPI_FVCOM_GROUP,IERR)
       SBUF = VYMIN
       CALL MPI_ALLREDUCE(SBUF,VYMIN,1,MPI_F,MPI_MIN,MPI_FVCOM_GROUP,IERR)
       !GLOBAL MAX
       SBUF = VXMAX
       CALL MPI_ALLREDUCE(SBUF,VXMAX,1,MPI_F,MPI_MAX,MPI_FVCOM_GROUP,IERR)
       SBUF = VYMAX
       CALL MPI_ALLREDUCE(SBUF,VYMAX,1,MPI_F,MPI_MAX,MPI_FVCOM_GROUP,IERR)
    END IF
#endif
    !--------------SHIFT GRID TO UPPER RIGHT CARTESIAN-----------------------------!
    
    VX = VX - VXMIN
    VY = VY - VYMIN
    VX(0) = 0.0_SP ; VY(0) = 0.0_SP
# if defined (DATA_ASSIM)
!-----------------add by Jadon Ge----------------------    
    IF(.NOT.ALLOCATED(XG))  ALLOCATE(XG(0:MGL))
    IF(.NOT.ALLOCATED(YG))  ALLOCATE(YG(0:MGL))
    IF(.NOT.ALLOCATED(XCG)) ALLOCATE(XCG(0:NGL))
    IF(.NOT.ALLOCATED(YCG)) ALLOCATE(YCG(0:NGL))
# if defined (MULTIPROCESSOR)
    IF(PAR)THEN
      SENDER = MSRID -1
      CALL MPI_BCAST(XG,MGL,MPI_F,SENDER,MPI_FVCOM_GROUP,IERR)
      CALL MPI_BCAST(YG,MGL,MPI_F,SENDER,MPI_FVCOM_GROUP,IERR)
      CALL MPI_BCAST(XCG,NGL,MPI_F,SENDER,MPI_FVCOM_GROUP,IERR)
      CALL MPI_BCAST(YCG,NGL,MPI_F,SENDER,MPI_FVCOM_GROUP,IERR)
    END IF
# endif
    XG = XG - VXMIN
    YG = YG - VYMIN
    XG(0) = 0.0_SP ; YG(0) = 0.0_SP

    XCG = XCG - VXMIN
    YCG = YCG - VYMIN
    XCG(0) = 0.0_SP ; YCG(0) = 0.0_SP
!------------------------------------------------------
# endif
    
    !--------------CALCULATE GLOBAL ELEMENT CENTER GRID COORDINATES----------------!
    CALL N2E2D(VX,XC)
    CALL N2E2D(VY,YC)
    
    XMC = XC + VXMIN
    XMC(0)= 0.0_SP
    YMC = YC + VYMIN
    YMC(0)= 0.0_SP

    IF (USE_PROJ ) THEN
        IF (SERIAL) THEN
           CALL Meters2Degrees(XMC(1:NT),YMC(1:NT)&
                &   ,PROJECTION_REFERENCE,LONC(1:NT),LATC(1:NT),NT)
# if defined (MULTIPROCESSOR)
        ELSE
            IF (MSR) THEN
                allocate(xc_buf(0:NGL),lonc_buf(0:NGL), stat=ierr)
                allocate(yc_buf(0:NGL),latc_buf(0:NGL), stat=ierr)
            END IF
            
            call ACOLLECT(MYID,MSRID,NPROCS,EMAP,XMC,xc_buf)
            call ACOLLECT(MYID,MSRID,NPROCS,EMAP,YMC,yc_buf)
            
            IF (MSR) THEN
                CALL Meters2Degrees(xc_buf(1:NGL),yc_buf(1:NGL)&
                    &   ,PROJECTION_REFERENCE,lonc_buf(1:NGL),latc_buf(1:NGL),NGL)
            END IF
            
            call ADEAL(MYID,MSRID,NPROCS,EXMAP,lonc_buf,lonc)
            call ADEAL(MYID,MSRID,NPROCS,EXMAP,latc_buf,latc)
            
            IF (MSR) THEN
                deallocate(xc_buf, lonc_buf)
                deallocate(yc_buf, latc_buf)
            END IF

# endif
        END IF

    END IF

    
# else
    ! SPHERICAL

    DO I=1,NT
       VX1=VX(NV(I,1))
       VY1=VY(NV(I,1))
       VX2=VX(NV(I,2))
       VY2=VY(NV(I,2))
       VX3=VX(NV(I,3))
       VY3=VY(NV(I,3))
       
       inner: DO K=1,1000000
          
          EVX12=VX2-VX1
          EVX13=VX3-VX1
          EVX23=VX3-VX2
          
          IF(EVX12 >  180.0_SP)THEN
             EVX12 = -360.0_SP+EVX12
          ELSE IF(EVX12 < -180.0_SP)THEN
             EVX12 =  360.0_SP+EVX12
          END IF

          IF(EVX13 >  180.0_SP)THEN
             EVX13 = -360.0_SP+EVX13
          ELSE IF(EVX13 < -180.0_SP)THEN
             EVX13 =  360.0_SP+EVX13
          END IF

          IF(EVX23 >  180.0_SP)THEN
             EVX23 = -360.0_SP+EVX23
          ELSE IF(EVX23 < -180.0_SP)THEN
             EVX23 =  360.0_SP+EVX23
          END IF
          
          EVX12=ABS(EVX12)
          EVX13=ABS(EVX13)
          EVX23=ABS(EVX23)
          
          EVY12=ABS(VY2-VY1)
          EVY13=ABS(VY3-VY1)
          EVY23=ABS(VY3-VY2)
          
          EVXY=1.E-10_SP
          
          IF((EVX12 < EVXY) .AND.(EVX13 < EVXY) .AND. (EVX23 < EVXY) &
               & .AND. (EVY12 < EVXY) .AND. (EVY13 < EVXY)  &
               & .AND. (EVY23 < EVXY)) THEN

             XC(I)=VX1
             YC(I)=VY1
             exit inner
          ELSE
             CALL ARCC(VX1,VY1,VX2,VY2,VX12,VY12)
             CALL ARCC(VX2,VY2,VX3,VY3,VX23,VY23)
             CALL ARCC(VX3,VY3,VX1,VY1,VX31,VY31)
             
             VX1=VX12
             VY1=VY12
             VX2=VX23
             VY2=VY23
             VX3=VX31
             VY3=VY31
          END IF
       END DO inner
       
    END DO

    
    XC(0) = 0.0_SP ; YC(0) = 0.0_SP

# if defined (DATA_ASSIM)
!-----------------add by Jadon Ge----------------------
# if defined (MULTIPROCESSOR)
    IF(.NOT.ALLOCATED(XG))  ALLOCATE(XG(0:MGL))
    IF(.NOT.ALLOCATED(YG))  ALLOCATE(YG(0:MGL))
    IF(.NOT.ALLOCATED(XCG)) ALLOCATE(XCG(0:NGL))
    IF(.NOT.ALLOCATED(YCG)) ALLOCATE(YCG(0:NGL))
    IF(PAR)THEN
      SENDER = MSRID -1
      CALL MPI_BCAST(XG,MGL,MPI_F,SENDER,MPI_FVCOM_GROUP,IERR)
      CALL MPI_BCAST(YG,MGL,MPI_F,SENDER,MPI_FVCOM_GROUP,IERR)
      CALL MPI_BCAST(XCG,NGL,MPI_F,SENDER,MPI_FVCOM_GROUP,IERR)
      CALL MPI_BCAST(YCG,NGL,MPI_F,SENDER,MPI_FVCOM_GROUP,IERR)
    END IF
# endif
    XG(0) = 0.0_SP ; YG(0) = 0.0_SP
    XCG(0) = 0.0_SP ; YCG(0) = 0.0_SP
!------------------------------------------------------
# endif

    LONC=XC
    LATC=YC

    IF (USE_PROJ ) THEN
        IF (SERIAL) THEN
           CALL Degrees2Meters(XC(1:NT),YC(1:NT)&
                &   ,PROJECTION_REFERENCE,XMC(1:NT),YMC(1:NT),NT)
# if defined (MULTIPROCESSOR)
        ELSE
            IF (MSR) THEN
                allocate(xc_buf(0:NGL),lonc_buf(0:NGL), stat=ierr)
                allocate(yc_buf(0:NGL),latc_buf(0:NGL), stat=ierr)
            END IF
            
            call ACOLLECT(MYID,MSRID,NPROCS,EMAP,XC,lonc_buf)
            call ACOLLECT(MYID,MSRID,NPROCS,EMAP,YC,latc_buf)
            
            IF (MSR) THEN
                CALL Degrees2Meters(lonc_buf(1:NGL),latc_buf(1:NGL)&
                    &   ,PROJECTION_REFERENCE,xc_buf(1:NGL),yc_buf(1:NGL),NGL)
            END IF
            
            call ADEAL(MYID,MSRID,NPROCS,EXMAP,xc_buf,xmc)
!JQI            call ADEAL(MYID,MSRID,NPROCS,EXMAP,xc_buf,ymc)
            call ADEAL(MYID,MSRID,NPROCS,EXMAP,yc_buf,ymc)
            
            IF (MSR) THEN
                deallocate(xc_buf, lonc_buf)
                deallocate(yc_buf, latc_buf)
            END IF

# endif
        END IF
        
    END IF


# endif

    IF(DBG_SET(DBG_SBR)) write(IPT,*) "SETUP_CENTER_COORDS: END"

  END SUBROUTINE SETUP_CENTER_COORDS
  !==============================================================================|
  !   SET HORIZONTAL MIXING_COEFFICIENT                                      |
  !==============================================================================|
  SUBROUTINE SETUP_HORIZONTAL_MIXING_COEFFICIENT
    IMPLICIT NONE
  
    if (HORIZONTAL_MIXING_KIND .eq. STTC) THEN

       if(DBG_SET(DBG_LOG)) then
          write(IPT,*) "! "
          write(IPT,*) "! Setting Staticly Variable viscosity"
       end if

       CALL LOAD_HORIZONTAL_MIXING_COEFFICIENT(NN_HVC,CC_HVC)
       
       
    else if(HORIZONTAL_MIXING_KIND .eq. CNSTNT) THEN
       
       CC_HVC=HORIZONTAL_MIXING_COEFFICIENT
       CC_HVC(0)=0.0_SP
       NN_HVC=HORIZONTAL_MIXING_COEFFICIENT
       NN_HVC(0)=0.0_SP
       
    else
       CALL FATAL_ERROR&
            &("HORIZONTAL_MIXING_KIND ERROR",&
            & "This should not happen")
       
    end if
    
    
  END SUBROUTINE SETUP_HORIZONTAL_MIXING_COEFFICIENT
  !==============================================================================|
  !   SET BOTTOM ROUGHNESS                                                       |
  !==============================================================================|
  SUBROUTINE SETUP_BOTTOM_ROUGHNESS
    IMPLICIT NONE
  
    if (BOTTOM_ROUGHNESS_KIND .eq. STTC) THEN
  
       if(DBG_SET(DBG_LOG)) then
          write(IPT,*) "! "
          write(IPT,*) "! Setting Staticly Variable Bottom Roughness"
       end if

!Changed by Siqi Li@2015/08/29
!       CALL LOAD_BOTTOM_ROUGHNESS(CC_Z0B)
       SELECT CASE(BOTTOM_ROUGHNESS_TYPE)
       CASE(BR_UDEF)
         CALL LOAD_DRAG_COEFFICIENT(CBC_UD)
       CASE DEFAULT
         CALL LOAD_BOTTOM_ROUGHNESS(CC_Z0B)
       END SELECT
!Siqi Li
              
    else if(BOTTOM_ROUGHNESS_KIND .eq. CNSTNT) THEN
       
       CC_Z0B    = BOTTOM_ROUGHNESS_LENGTHSCALE
       CC_Z0B(0) = 0.0_SP
       
    else
       CALL FATAL_ERROR&
            &("HORIZONTAL_MIXING_KIND ERROR",&
            & "This should not happen")
       
    end if
    
    
  END SUBROUTINE SETUP_BOTTOM_ROUGHNESS
  !==============================================================================|
  !   SET UP LOCAL MESH (BATHYMETRIC DEPTH)                                      |
  !==============================================================================|
  SUBROUTINE SETUP_DEPTH
    IMPLICIT NONE
    INTEGER :: IERR
    REAL(SP):: SBUF
    
    IF(DBG_SET(DBG_SBR)) write(IPT,*) "SETUP_DEPTH: START"
    HMAX = MAXVAL(H_LCL(1:MT))
    HMIN = MINVAL(H_LCL(1:MT))

# if  defined (MULTIPROCESSOR)
    IF(PAR)THEN
       SBUF = HMAX
       CALL MPI_ALLREDUCE(SBUF,HMAX,1,MPI_F,MPI_MAX,MPI_FVCOM_GROUP,IERR)
       SBUF = HMIN
       CALL MPI_ALLREDUCE(SBUF,HMIN,1,MPI_F,MPI_MIN,MPI_FVCOM_GROUP,IERR)
    END IF
# endif

    H = H_LCL

    IF(DBG_SET(DBG_SBR)) write(IPT,*) "SETUP_DEPTH: END"

  END SUBROUTINE SETUP_DEPTH

!==============================================================================|
!   SET UP LOCAL CORIOLIS FORCE                                                |
!==============================================================================|
  SUBROUTINE SETUP_CORIOLIS 
    IMPLICIT NONE
    integer:: I
    !--------------READ IN CORIOLIS PARAMETER--------------------------------------!

    IF(DBG_SET(DBG_SBR)) write(IPT,*) "SETUP_CORIOLIS: START"

    CALL N2E2D(C_LCL,COR)

#   if defined (SPHERICAL)
    IF(EQUATOR_BETA_PLANE)THEN
      F_ALFA = 1.0_SP-0.8_SP*EXP(-(COR/2.2_SP)**2)
    END IF 
#   endif

    COR = 2.*7.292e-5_SP * SIN(COR * DEG2RAD)
    
    !!  ggao for equatoral min (4deg)
!JQI20201123 -- start comment
!    IF(.NOT. EQUATOR_BETA_PLANE)THEN
!     WHERE(COR <  1.e-5_SP .AND. COR > 0.0_SP) COR =  1.e-5_SP
!     WHERE(COR > -1.e-5_SP .AND. COR < 0.0_SP) COR = -1.e-5_SP
!    END IF 
!JQI20201123 -- end comment

    IF(DBG_SET(DBG_SBR)) write(IPT,*) "SETUP_CORIOLIS: END"

  END SUBROUTINE SETUP_CORIOLIS
    
    
!==============================================================================|
!   COMPUTE GRAVITY VARIED WITH LATITUDE                                       |
!==============================================================================|

  SUBROUTINE SETUP_GRAVITY
    IMPLICIT NONE
    INTEGER I

    IF(DBG_SET(DBG_SBR)) write(IPT,*) "SETUP_GRAVITY: START"

# if defined (SPHERICAL)
    DO I=1,MT
       GRAV_N(I) = 9.780327*(1.0+0.0053024*SIN(VY(I)*DEG2RAD)**2-0.0000058*SIN(2.0*VY(I)*DEG2RAD)**2)
    END DO
    DO I=1,NT
       GRAV_E(I) = 9.780327*(1.0+0.0053024*SIN(YC(I)*DEG2RAD)**2-0.0000058*SIN(2.0*YC(I)*DEG2RAD)**2)
    END DO
# else
    GRAV_N = GRAV
    GRAV_E = GRAV
# endif        

    IF(DBG_SET(DBG_SBR)) write(IPT,*) "SETUP_GRAVITY: END"
  END SUBROUTINE SETUP_GRAVITY
    
!==============================================================================|
!   COMPUTE SPONGE LAYER FOR OPEN BOUNDARY DAMPING                             |
!==============================================================================|

  SUBROUTINE SETUP_SPONGE
    USE MOD_SPHERICAL
    IMPLICIT NONE
    REAL(SP)  TEMP,DTMP,C_SPONGE
    INTEGER :: I1, I, SENDER, IERR
    REAL(DP) X1_DP,Y1_DP,X2_DP,Y2_DP,DTMP_DP

    !--SET SPONGE PARAMETERS-------------------------------------------------------|
    
    IF(DBG_SET(DBG_SBR)) write(IPT,*) "SETUP_SPONGE: START"

    IF (NSPONGE ==0) RETURN

! NOTE: X_SPG/Y_SPG COORDINATES MUST BE AJUSTED FOR VXMIN/VYMIN

# if !defined (SPHERICAL)
    X_SPG = X_SPG - VXMIN
    Y_SPG = Y_SPG - VYMIN

    DO I=1,NT
       DO I1=1,NSPONGE
          DTMP=(XC(I)-X_SPG(I1))**2+(YC(I)-Y_SPG(I1))**2
          DTMP=SQRT(DTMP)/R_SPG(I1)

          IF(DTMP <= 1.0_SP) THEN
             C_SPONGE=C_SPG(I1)*(1.0_SP-DTMP)
             CC_SPONGE(I)=MAX(C_SPONGE,CC_SPONGE(I))
          END IF
       END DO
    END DO

# else
    ! SPHERICAL

    DO I=1,NT
       DO I1=1,NSPONGE
          X1_DP=XC(I)
          Y1_DP=YC(I)
          X2_DP=X_SPG(I1)
          Y2_DP=Y_SPG(I1)
          CALL ARC(X1_DP,Y1_DP,X2_DP,Y2_DP,DTMP_DP)
          DTMP=DTMP_DP/R_SPG(I1)

          IF(DTMP <= 1.0_SP) THEN
             C_SPONGE=C_SPG(I1)*(1.0_SP-DTMP)
             CC_SPONGE(I)=MAX(C_SPONGE,CC_SPONGE(I))
          END IF
       END DO
    END DO

# endif    

    DEALLOCATE(N_SPG,R_SPG,C_SPG,X_SPG,Y_SPG)

    IF(DBG_SET(DBG_SBR)) write(IPT,*) "SETUP_SPONGE: END"
    
    RETURN
  END SUBROUTINE SETUP_SPONGE
!==============================================================================|
  
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|

!==============================================================================|
  SUBROUTINE COORDINATE_UNITS(XL,YL)
    IMPLICIT NONE
    REAL(SP), ALLOCATABLE :: XL(:),YL(:)
    integer status, ierr


    REAL(SP), allocatable :: x_buf(:), lon_buf(:)
    REAL(SP), allocatable :: y_buf(:), lat_buf(:)

    IF(DBG_SET(DBG_SBR)) write(IPT,*) "COODINATE_UNITS: START"

# if defined (SPHERICAL)

    IF (GRID_FILE_UNITS == 'degrees') THEN
       ! COPY DATA INTO VX,VY
       VX = XL
       VY = YL
       
       LON = XL
       LAT = YL
       
       IF ( USE_PROJ ) THEN
          ! USE PROJECTION TOOL BOX TO CONVERT TO METERS
        IF (SERIAL) THEN
           CALL Degrees2meters(XL(1:MT),YL(1:MT), &
               & PROJECTION_REFERENCE,XM(1:MT),YM(1:MT),MT)
# if defined (MULTIPROCESSOR)
        ELSE
            IF (MSR) THEN
                allocate(x_buf(0:MGL),lon_buf(0:MGL), stat=ierr)
                allocate(y_buf(0:MGL),lat_buf(0:MGL), stat=ierr)
            END IF
            
            call ACOLLECT(MYID,MSRID,NPROCS,NMAP,XL,lon_buf)
            call ACOLLECT(MYID,MSRID,NPROCS,NMAP,YL,lat_buf)
            
            IF (MSR) THEN
                CALL Degrees2meters(lon_buf(1:MGL),lat_buf(1:MGL)&
                    &   ,PROJECTION_REFERENCE,x_buf(1:MGL),y_buf(1:MGL),MGL)
            END IF
            
            call ADEAL(MYID,MSRID,NPROCS,NXMAP,x_buf,XM)
            call ADEAL(MYID,MSRID,NPROCS,NXMAP,y_buf,YM)
            
            IF (MSR) THEN
                deallocate(x_buf, lon_buf)
                deallocate(y_buf, lat_buf)
            END IF
          
# endif
          END IF  
        END IF
       
    ELSE IF(GRID_FILE_UNITS == 'meters') THEN

       IF ( USE_PROJ ) THEN

          ! COPY DATA INTO VX,VY
          XM = XL
          YM = YL
          
          ! USE PROJECTION TOOL BOX TO CONVERT TO METERS
        IF (SERIAL) THEN
           CALL meters2degrees(XL(1:MT),YL(1:MT), &
            & PROJECTION_REFERENCE,VX(1:MT),VY(1:MT),MT)
# if defined (MULTIPROCESSOR)
        ELSE
            IF (MSR) THEN
                allocate(x_buf(0:MGL),lon_buf(0:MGL), stat=ierr)
                allocate(y_buf(0:MGL),lat_buf(0:MGL), stat=ierr)
            END IF

            call ACOLLECT(MYID,MSRID,NPROCS,NMAP,XL,x_buf)
            call ACOLLECT(MYID,MSRID,NPROCS,NMAP,YL,y_buf)

            IF (MSR) THEN
                CALL meters2degrees(x_buf(1:MGL),y_buf(1:MGL)&
                    &   ,PROJECTION_REFERENCE,lon_buf(1:MGL),lat_buf(1:MGL),MGL)
            END IF

            call ADEAL(MYID,MSRID,NPROCS,NXMAP,lon_buf,VX)
            call ADEAL(MYID,MSRID,NPROCS,NXMAP,lat_buf,VY)

            IF (MSR) THEN
                deallocate(x_buf, lon_buf)
                deallocate(y_buf, lat_buf)
            END IF

# endif
        END IF
       
          LON = VX
          LAT = VY

       ELSE
           CALL FATAL_ERROR('You must specify a valid projection reference'&
               &,'and compile with PROJ to use files with cartesian coordinates in spherical mode')

       END IF

       
    ELSE 
       CALL FATAL_ERROR('UNRECOGNIZED GRID_FILE_UNITS: '//GRID_FILE_UNITS)
       
    END IF
    
    
# else 
!! IF NOT SPHERICAL CASE
    
    IF (GRID_FILE_UNITS == 'degrees') THEN
       
       IF (USE_PROJ) THEN
          
          LON = XL
          LAT = YL
          
          ! USE PROJECTION TOOL BOX TO CONVERT TO METERS
        IF (SERIAL) THEN
           CALL Degrees2meters(XL(1:MT),YL(1:MT), &
               & PROJECTION_REFERENCE,VX(1:MT),VY(1:MT),MT)
# if defined (MULTIPROCESSOR)
        ELSE
            IF (MSR) THEN
                allocate(x_buf(0:MGL),lon_buf(0:MGL), stat=ierr)
                allocate(y_buf(0:MGL),lat_buf(0:MGL), stat=ierr)
            END IF
            
            call ACOLLECT(MYID,MSRID,NPROCS,NMAP,XL,lon_buf)
            call ACOLLECT(MYID,MSRID,NPROCS,NMAP,YL,lat_buf)
            
            IF (MSR) THEN
                CALL Degrees2meters(lon_buf(1:MGL),lat_buf(1:MGL)&
                    &   ,PROJECTION_REFERENCE,x_buf(1:MGL),y_buf(1:MGL),MGL)
            END IF
            
            call ADEAL(MYID,MSRID,NPROCS,NXMAP,x_buf,VX)
            call ADEAL(MYID,MSRID,NPROCS,NXMAP,y_buf,VY)
            
            IF (MSR) THEN
                deallocate(x_buf, lon_buf)
                deallocate(y_buf, lat_buf)
            END IF
          
# endif
        END IF         
          XM = VX
          YM = VY
          
       ELSE
          
          CALL FATAL_ERROR('You must specify a valid projection reference'&
               &,'and compile with PROJ to use files with latitude and longitue in cartesian mode')
       END IF
       

    ELSE IF(GRID_FILE_UNITS == 'meters') THEN

       VX = XL
       VY = YL
       
       XM = XL
       YM = YL
       
       ! USE PROJECTION TOOL BOX TO CONVERT TO DEGREES
       IF (USE_PROJ ) THEN
        IF (SERIAL) THEN
           CALL meters2degrees(XL(1:MT),YL(1:MT), &
            & PROJECTION_REFERENCE,LON(1:MT),LAT(1:MT),MT)
# if defined (MULTIPROCESSOR)
        ELSE
            IF (MSR) THEN
                allocate(x_buf(0:MGL),lon_buf(0:MGL), stat=ierr)
                allocate(y_buf(0:MGL),lat_buf(0:MGL), stat=ierr)
            END IF

            call ACOLLECT(MYID,MSRID,NPROCS,NMAP,XL,x_buf)
            call ACOLLECT(MYID,MSRID,NPROCS,NMAP,YL,y_buf)

            IF (MSR) THEN
                CALL meters2degrees(x_buf(1:MGL),y_buf(1:MGL)&
                    &   ,PROJECTION_REFERENCE,lon_buf(1:MGL),lat_buf(1:MGL),MGL)
            END IF

            call ADEAL(MYID,MSRID,NPROCS,NXMAP,lon_buf,LON)
            call ADEAL(MYID,MSRID,NPROCS,NXMAP,lat_buf,LAT)

            IF (MSR) THEN
                deallocate(x_buf, lon_buf)
                deallocate(y_buf, lat_buf)
            END IF

# endif
        END IF
          
       END IF

    ELSE 
       CALL FATAL_ERROR('UNRECOGNIZED GRID_FILE_UNITS: '//TRIM(GRID_FILE_UNITS))
       
    END IF

# endif

    IF(DBG_SET(DBG_SBR)) write(IPT,*) "COODINATE_UNITS: END"
    
  END SUBROUTINE COORDINATE_UNITS

!==============================================================================|
! SETUP THE SIGMA COORDINATES FOR THE MODEL                                    |
!==============================================================================|
!==============================================================================|
! This program is used to set up the coordinate in the vertical.               !
!								               !
! case(1): sigma levels                                                        !
! sigma levels are determined by a formula of                                  !
!                      sigma(k)=-[(k-1)/(kb-1)]^k11                            !
!    p_sigma=1: uniform sigma layers                                           !
!    p_sigma=2: layers satisfying a parabolic function with high               !
!               vertical resolution near the surface and bottom.               !
!    p_sigma can be used any real number                                       !
!									       !
! case(2): general vertical level                                              !
! vertical levels are determined by the formula                                !
!                tanh[(dl+du)((kbm1-k)/kbm1)-dl]+tanh(dl)                      !
!        z(k)= ------------------------------------------  - 1                 !
!                      tanh(dl) + tanh(du)                                     !
!                                                                              !
! case(3): constant layer transformation                                       !
! four values need to be specified:                                            !
!  DUU the upper boundaries, up to which the co-ordinates are parallel must be !
!      defined.                                                                !
!  DLL the lower boundaries, up to which the co-ordinates are parallel must be !
!      defined.                                                                !
!  HMIN1 the minimum water depth at which the layers are constant. If H < HMIN1!
!      then sigma co-ordinates are used.                                       !
!                                                                              !
! Reference of case(2), case(3) and case(4):                                   !
! Pietrzak, J.D., Jan B. Jakobson, Hans Burchard, Hans Jacob Vested, Ole       !
! Petersen , 2002. A three-dimensional hydrostatic model for coastal and ocean !
! modelling using a generalised topography following co-ordinate system. Ocean !
! Modelling 4, 173-205                                                         !
!                                                                              !
!  calculates: z(m,kb) vertical levels					       !
!  calculates: dz(m,kb-1) delta between vertical levels		               !
!  calculates: zz(m,kb-1) intra-vertical levels				       !
!  calculates: dzz(m,kb-2) delta between intra-vertical levels		       !
!==============================================================================|

  SUBROUTINE SETUP_SIGMA          
    !==============================================================================|
    IMPLICIT NONE
    INTEGER :: K,KK
    INTEGER :: I
    REAL(SP):: ZTMP(KB)
    REAL(SP):: X1,X2,X3  
    REAL(SP):: DR,RCL,RCU
    !==============================================================================|

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*)"SETUP_SIGMA: START"

    IF(DBG_SET(DBG_SBRIO)) THEN
       WRITE(IPT,*)"==================="
       WRITE(IPT,*)"  SET_SIGMA IO"
       WRITE(IPT,*)"  STYPE = "//TRIM(STYPE)
       WRITE(IPT,*)"  P_SIGMA = ", P_SIGMA
       WRITE(IPT,*)"  KB = ", KB
       WRITE(IPT,*)"==================="
    END IF

    !---------------------------------------------
    !---------------------------------------------
    SELECT CASE(STYPE)
    !---------------------------------------------
    !SIGMA_COORDINATE_TYPE = UNIFORM (DEGENERATE CASE OF GEOMETRIC)
    CASE(STYPE_UNIFORM)
    !---------------------------------------------
       IF(P_SIGMA > 1 .AND. MOD(KB,2) == 0) &
            CALL FATAL_ERROR('SETUP_SIGMA: COORDINATE TYPE:'//trim(STYPE)&
            &,'kb shoude be an odd number for this type of sigma coordinates....' )
       CALL SIGMA_GEOMETRIC
    !---------------------------------------------
    !SIGMA_COORDINATE_TYPE = GEOMETRIC
    CASE(STYPE_GEOMETRIC)
    !---------------------------------------------

       IF(P_SIGMA > 1 .AND. MOD(KB,2) == 0) &
            CALL FATAL_ERROR('SETUP_SIGMA: COORDINATE TYPE:'//trim(STYPE)&
            &,'kb shoude be an odd number for this type of sigma coordinates....' )
       CALL SIGMA_GEOMETRIC
    !---------------------------------------------
    ! SIGMA_COORDINATE_TYPE = TANH
    CASE(STYPE_TANH)
    !---------------------------------------------   
       CALL SIGMA_TANH
    !---------------------------------------------   
    !SIGMA_COORDINATE_TYPE = GENERALIZED
    CASE(STYPE_GENERALIZED)
    !---------------------------------------------   
       CALL SIGMA_GENERALIZED

    ! THIS IS A CURRENTLY UNUSED METHOD
!!$!---------------------------------------------   
!!$!SIGMA_COORDINATE_TYPE = WHAT THE HELL IS THIS?
!!$   CASE("UNKNOWN")
!!$!---------------------------------------------   
!!$
!!$    CALL SIGMA_UNKNOWN
    CASE DEFAULT
       CALL FATAL_ERROR("SET_SIGMA: REACHED DEFAULT CASE FOR SIGMA COOR&
            &DINATE TYPE")
    END SELECT

    !---------COMPUTE SIGMA DERIVATIVES AND INTRA SIGMA LEVELS---------------------!
#    if defined (MULTIPROCESSOR)
    IF(PAR)CALL AEXCHANGE(EC,MYID,NPROCS,Z1)
    IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,Z)
#    endif

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*)"SETUP_SIGMA: END"

    RETURN

  END SUBROUTINE SETUP_SIGMA
!==============================================================================|
  SUBROUTINE SIGMA_GEOMETRIC
      IMPLICIT NONE
      INTEGER :: I,K
      REAL(SP):: ZTMP(KB)
      !orginal formula to set sigma   
      IF(P_SIGMA == 1)THEN
         DO K=1,KB
            ZTMP(K) = -((K-1)/FLOAT(KB-1))**P_SIGMA 
         END DO
      ELSE
         DO K=1,(KB+1)/2
            ZTMP(K) = -((K-1)/FLOAT((KB+1)/2-1))**P_SIGMA/2 
         END DO
         DO K=(KB+1)/2+1,KB
            ZTMP(K) = ((KB-K)/FLOAT((KB+1)/2-1))**P_SIGMA/2-1.0
         END DO
      END IF

      DO I=1,M
         DO K=1,KB
            Z(I,K)=ZTMP(K)
         END DO
      END DO

      DO I=1,N
         DO K=1,KB
            Z1(I,K)=(Z(NV(I,1),K)+Z(NV(I,2),K)+Z(NV(I,3),K))/3.0_SP
         END DO
      END DO
    END SUBROUTINE SIGMA_GEOMETRIC
!--------------------------------------------------------------------    
    SUBROUTINE SIGMA_GENERALIZED
      IMPLICIT NONE
      INTEGER :: I,K, kk
      REAL(SP):: X1,X2,X3  
      REAL(SP):: DR,RCL,RCU
 
      DO I=1,M
         IF(H(I) < HMIN1)THEN
            Z(I,1)=0.0
            DL2=0.001;DU2=0.001
            DO K=1,KBM1
               X1=DL2+DU2
               X1=X1*(KBM1-K)/KBM1
               X1=X1-DL2
               X1=TANH(X1)
               X2=TANH(DL2)
               X3=X2+TANH(DU2)

               Z(I,K+1)=(X1+X2)/X3-1.0_SP
            END DO
         ELSE
            DR=(H(I)-DUU-DLL)/H(I)/(KB-KU-KL-1)

            Z(I,1)=0.0_SP

            DO K=2,KU+1
               Z(I,K)=Z(I,K-1)-ZKU(K-1)/H(I)
            END DO

            DO K=KU+2,KB-KL
               Z(I,K)=Z(I,K-1)-DR
            END DO

            KK=0
            DO K=KB-KL+1,KB
               KK=KK+1
               Z(I,K)=Z(I,K-1)-ZKL(KK)/H(I)
            END DO
         END IF
      END DO

      DO I=1,N
         DO K=1,KB
            Z1(I,K)=(Z(NV(I,1),K)+Z(NV(I,2),K)+Z(NV(I,3),K))/3.0_SP
         END DO
      END DO
    END SUBROUTINE SIGMA_GENERALIZED
!--------------------------------------------------------------------
    SUBROUTINE SIGMA_TANH
      IMPLICIT NONE
      INTEGER :: I,K
      REAL(SP):: X1,X2,X3  
    
      Z=0.0;Z1=0.0
      DO K=1,KBM1
         X1=DL2+DU2
         X1=X1*(KBM1-K)/KBM1
         X1=X1-DL2
         X1=TANH(X1)
         X2=TANH(DL2)
         X3=X2+TANH(DU2)
         DO I=1,M
            Z(I,K+1)=(X1+X2)/X3-1.0_SP
         END DO
         DO I=1,N
            Z1(I,K+1)=(X1+X2)/X3-1.0_SP
         END DO
      END DO
    END SUBROUTINE SIGMA_TANH
!--------------------------------------------------------------------
!!$ SUBROUTINE SIGMA_UNKNOWN
!!$   DO I=1,M
!!$     IF(H(I) < HMIN1)THEN
!!$       RCU=-DUU/HMIN1
!!$       RCL=DLL/HMIN1-1
!!$       DR=(RCL-RCU)/(KB-KU-KL-1)
!!$
!!$       DO K=1,KU
!!$         ZKU(K)=RCU/KU
!!$       END DO
!!$       DO K=1,KL
!!$         ZKL(K)=(-1.0_SP-RCL)/KL
!!$       END DO
!!$
!!$       Z(I,1)=0.0_SP
!!$
!!$       DO K=2,KU+1
!!$         Z(I,K)=Z(I,K-1)+ZKU(K-1)
!!$       END DO
!!$
!!$       DO K=KU+2,KB-KL
!!$         Z(I,K)=Z(I,K-1)+DR
!!$       END DO
!!$
!!$       KK=0
!!$       DO K=KB-KL+1,KB
!!$         KK=KK+1
!!$         Z(I,K)=Z(I,K-1)+ZKL(KK)
!!$       END DO
!!$
!!$     ELSE
!!$       DR=(H(I)-DUU-DLL)/H(I)/(KB-KU-KL-1)
!!$
!!$       Z(I,1)=0.0_SP
!!$
!!$       DO K=2,KU+1
!!$         Z(I,K)=Z(I,K-1)-ZKU(K-1)/H(I)
!!$       END DO
!!$
!!$       DO K=KU+2,KB-KL
!!$         Z(I,K)=Z(I,K-1)-DR
!!$       END DO
!!$
!!$       KK=0
!!$       DO K=KB-KL+1,KB
!!$         KK=KK+1
!!$         Z(I,K)=Z(I,K-1)-ZKL(KK)/H(I)
!!$       END DO
!!$     END IF
!!$   END DO
!!$   
!!$   DO I=1,N
!!$     DO K=1,KB
!!$       Z1(I,K)=(Z(NV(I,1),K)+Z(NV(I,2),K)+Z(NV(I,3),K))/3.0_SP
!!$     END DO
!!$   END DO    
!!$END SUBROUTINE SIGMA_UNKNOWN

!==============================================================================|

  SUBROUTINE SETUP_SIGMA_DERIVATIVES
    USE ALL_VARS
    IMPLICIT NONE
    INTEGER :: K, I


    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*)"SETUP_SIGMA_DERIVATIVES: START"

    DO K=1,KB-1
       DO I=1,MT
          DZ(I,K)  = Z(I,K)-Z(I,K+1)
          ZZ(I,K)  = .5_SP*(Z(I,K)+Z(I,K+1))
       END DO
       DO I=1,NT
          DZ1(I,K)  = Z1(I,K)-Z1(I,K+1)
          ZZ1(I,K)  = .5_SP*(Z1(I,K)+Z1(I,K+1))
       END DO
    END DO

    DO I=1,MT
       ZZ(I,KB) = 2.0_SP*ZZ(I,KB-1)-ZZ(I,KB-2)
    END DO
    DO I=1,NT
       ZZ1(I,KB) = 2.0_SP*ZZ1(I,KB-1)-ZZ1(I,KB-2)
    END DO

    DO K=1,KBM2
       DO I=1,MT
          DZZ(I,K) = ZZ(I,K)-ZZ(I,K+1)
       END DO
       DO I=1,NT
          DZZ1(I,K) = ZZ1(I,K)-ZZ1(I,K+1)
       END DO
    END DO

    DZZ(:,KBM1) = 0.0_SP
    DZ(:,KB)    = 0.0_SP
    DZZ1(:,KBM1) = 0.0_SP
    DZ1(:,KB)    = 0.0_SP


    !----------OUTPUT VALUES-TO INFOFILE-------------------------------------------!

    IF(DBG_SET(DBG_LOG)) THEN
       WRITE(IPT,*  )'!'
       WRITE(IPT,*  )'!'
       WRITE(IPT,*)'!                SIGMA LAYER INFO     '
       WRITE(IPT,*) "SIGMA TYPE:",TRIM(STYPE)
       WRITE(IPT,70)
       SELECT CASE(STYPE)
       CASE(STYPE_UNIFORM)
          DO K=1,KB
             WRITE(IPT,80) K,Z(1,K),ZZ(1,K),DZ(1,K),DZZ(1,K)
          END DO
       CASE(STYPE_RESTART)
          DO K=1,KB
             WRITE(IPT,80) K,Z(1,K),ZZ(1,K),DZ(1,K),DZZ(1,K)
          END DO
       CASE(STYPE_GEOMETRIC)
          DO K=1,KB
             WRITE(IPT,80) K,Z(1,K),ZZ(1,K),DZ(1,K),DZZ(1,K)
          END DO
       CASE(STYPE_TANH)
          DO K=1,KB
             WRITE(IPT,80) K,Z(1,K),ZZ(1,K),DZ(1,K),DZZ(1,K)
          END DO
       CASE(STYPE_GENERALIZED) ! THIS IS CASE SPECIFIC 
          WRITE(IPT,*)"SET CASE SPECIFIC GENERALIZED SIGMA LAYER OUTPUT TO SCREEN &
               &IN mod_setup.F"
       END SELECT
       WRITE(IPT,*  )'!'
    END IF

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*)"END SETUP_SIGMA_DERIVATIVES"

    !----------FORMAT STATEMENTS---------------------------------------------------!

70  FORMAT(2x,'k',13x,'z',11x,'zz',11x,'dz',11x,'dzz')
80  FORMAT(' ',i5,4f13.8)


    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*)"SETUP_SIGMA_DERIVATIVES: END"
  END SUBROUTINE SETUP_SIGMA_DERIVATIVES
  
END MODULE MOD_SETUP

