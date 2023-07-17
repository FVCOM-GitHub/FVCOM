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

MODULE EQS_OF_STATE


  ! THIS MODULE CONTAINS THREE EQUATIONS OF STATE WHICH CAN BE ACCESSED
  ! USING A NUMBER OF DIFFERENT INTERFACES. THE ORIGINAL FVCOM
  ! INTERFACE (DENS1,DENS2,DENS3) STILL EXISTS, WHERE BOTH RHO AND
  ! RHO1 ARE UPDATED, OR THE USER CAN ACCESS THE UNDERLIEING
  ! FUNCTIONS AND PASS AN ARRAY, 1D OR 2D TO THE SUBROUTINE.
  !
  ! SEE DETAILS OF METHOD IN THE HEADER FOR EACH SUBROUTINE
  !
  ! SUBROUTINES PUBLICLY AVAILABLE IN THIS MODULE:
  !=================================================================================
  ! DENS1 - INTERFACE: CALL DENS1
  ! DENS2 - INTERFACE: CALL DENS2
  ! DENS3 - INTERFACE: CALL DENS3
  ! - THESE ROUTINES USE THE VALUES IN S1,T1,ZZ TO UPDATE THE DENSITY IN RHO AND RHO1
  !==================================================================================
  ! FOFONOFF_MILLARD - INTERFACE: CALL FOFONOFF_MILLARD_2D(S,T,Z,PREF,RHO)
  ! REAL(SP), INTENT(IN),DIMENSION(:) :: S,T,P    ! SALINITY,TEMPERATURE, PRESSURE
  ! REAL(SP), INTENT(OUT),DIMENSION(:) :: RHO     ! RESULT - DENSITY
  ! REAL(SP), INTENT(IN) :: PREF                  ! REFERENCE PRESSURE
  !    or
  ! REAL(SP), INTENT(IN),DIMENSION(:,:) :: S,T,P  ! SALINITY,TEMPERATURE, PRESSURE
  ! REAL(SP), INTENT(OUT),DIMENSION(:,:) :: RHO   ! RESULT - DENSITY
  ! REAL(SP), INTENT(IN) :: PREF                  ! REFERENCE PRESSURE
  !=================================================================================
  ! 
  !=================================================================================
  ! DENS2G - INTERFACE: CALL DENS2G(S,T,RHO)
  ! REAL(SP), INTENT(IN),DIMENSION(:) :: S,T      ! SALINITY,TEMPERATURE
  ! REAL(SP), INTENT(OUT),DIMENSION(:) :: RHO     ! DENSITY
  !    or
  ! REAL(SP), INTENT(IN),DIMENSION(:,:) :: S,T,P  ! SALINITY,TEMPERATURE
  ! REAL(SP), INTENT(OUT),DIMENSION(:,:) :: RHO   ! DENSITY
  !=================================================================================
  ! 
  !=================================================================================
  ! JACKET_MCDOUGALL - INTERFACE: CALL JACKET_MCDOUGALL(S,T,Z,RHO)
  ! REAL(SP), INTENT(IN),DIMENSION(:) :: S,T,P    ! SALINITY,TEMPERATURE, PRESSURE
  ! REAL(SP), INTENT(OUT),DIMENSION(:) :: RHO     ! RESULT - DENSITY
  !    or
  ! REAL(SP), INTENT(IN),DIMENSION(:,:) :: S,T,P  ! SALINITY,TEMPERATURE, PRESSURE
  ! REAL(SP), INTENT(OUT),DIMENSION(:,:) :: RHO   ! RESULT - DENSITY
  !=================================================================================
  !=================================================================================
  !=================================================================================

  USE ALL_VARS
  USE MOD_UTILS
  IMPLICIT NONE

  PUBLIC

  INTERFACE FOFONOFF_MILLARD
     MODULE PROCEDURE FOFONOFF_MILLARD_1D
     MODULE PROCEDURE FOFONOFF_MILLARD_2D
  END INTERFACE

  INTERFACE DENS2G
     MODULE PROCEDURE DENS2_1D
     MODULE PROCEDURE DENS2_2D
  END INTERFACE

  INTERFACE JACKET_MCDOUGALL
     MODULE PROCEDURE JACKET_MCDOUGALL_1D
     MODULE PROCEDURE JACKET_MCDOUGALL_2D
  END INTERFACE

  PRIVATE JACKET_MCDOUGALL_1D
  PRIVATE JACKET_MCDOUGALL_2D

  PRIVATE DENS2_1D
  PRIVATE DENS2_2D

  PRIVATE FOFONOFF_MILLARD_1D
  PRIVATE FOFONOFF_MILLARD_2D

  PRIVATE SVAN
  PRIVATE THETA
  PRIVATE ATG

CONTAINS

  !==============================================================================|
  !   Calculate Potential Density Based on Potential Temp and Salinity           |
  !     Pressure effects are incorported (Can Model Fresh Water < 4 Deg C)       |
  !     Ref:  algorithms for computation of fundamental properties of            |
  !         seawater , Fofonoff and Millard.				       |
  !                                                                              |
  !  calculates: rho1(nnode) density at nodes				       |
  !  calculates: rho (ncell) density at elements				       |
  !==============================================================================|
  SUBROUTINE FOFONOFF_MILLARD_2D(MYS,MYT,MYP,PREF, MYRHO)               
    
    !------------------------------------------------------------------------------|
    IMPLICIT NONE
    REAL(SP), INTENT(IN), DIMENSION(:,:) :: MYS,MYT,MYP
    REAL(SP), INTENT(IN) :: PREF
    REAL(SP), INTENT(INOUT), DIMENSION(:,:):: MYRHO
    INTEGER :: I,K,ub1,lb1,ub2,lb2
    REAL(SP)  ::PT,SVA,SIGMA
    !==============================================================================|

    ! SET DIMS USING MY S
    lb1 = lbound(mys,1)
    ub1 = ubound(mys,1)

    lb2 = lbound(mys,2)
    ub2 = ubound(mys,2)

    ! CHECK MY T
    if(lb1 /= lbound(myt,1)) CALL FATAL_ERROR("EQS_OF_STATE:: Dimension mismatch!")
    if(ub1 /= ubound(myt,1)) CALL FATAL_ERROR("EQS_OF_STATE:: Dimension mismatch!")

    if(lb2 /= lbound(myt,2)) CALL FATAL_ERROR("EQS_OF_STATE:: Dimension mismatch!")
    if(ub2 /= ubound(myt,2)) CALL FATAL_ERROR("EQS_OF_STATE:: Dimension mismatch!")

    ! CHECK MY P
    if(lb1 /= lbound(myp,1)) CALL FATAL_ERROR("EQS_OF_STATE:: Dimension mismatch!")
    if(ub1 /= ubound(myp,1)) CALL FATAL_ERROR("EQS_OF_STATE:: Dimension mismatch!")

    if(lb2 /= lbound(myp,2)) CALL FATAL_ERROR("EQS_OF_STATE:: Dimension mismatch!")
    if(ub2 /= ubound(myp,2)) CALL FATAL_ERROR("EQS_OF_STATE:: Dimension mismatch!")

    ! CHECK MY RHO
    if(lb1 /= lbound(myrho,1)) CALL FATAL_ERROR("EQS_OF_STATE:: Dimension mismatch!")
    if(ub1 /= ubound(myrho,1)) CALL FATAL_ERROR("EQS_OF_STATE:: Dimension mismatch!")

    if(lb2 /= lbound(myrho,2)) CALL FATAL_ERROR("EQS_OF_STATE:: Dimension mismatch!")
    if(ub2 /= ubound(myrho,2)) CALL FATAL_ERROR("EQS_OF_STATE:: Dimension mismatch!")

    

    DO I=lb1,ub1
       DO K=lb2,ub2  

          PT  = THETA(MYS(I,K),MYT(I,K),MYP(I,K),PREF)
          SVA = SVAN(MYS(I,K),PT,MYP(I,K),SIGMA)
          MYRHO(I,K) = SIGMA*1.e-3_SP
       END DO
    END DO

    !----------------TRANSFORM TO FACE CENTER--------------------------------------

  END SUBROUTINE FOFONOFF_MILLARD_2D
  !==============================================================================|
  SUBROUTINE FOFONOFF_MILLARD_1D(MYS,MYT,MYP,PREF, MYRHO)               
    
    !------------------------------------------------------------------------------|
    IMPLICIT NONE
    REAL(SP), INTENT(IN), DIMENSION(:) :: MYS,MYT,MYP
    REAL(SP), INTENT(IN) :: PREF
    REAL(SP), INTENT(INOUT), DIMENSION(:):: MYRHO
    INTEGER :: I,K,ub1,lb1
    REAL(SP)  ::PT,SVA,SIGMA
    !==============================================================================|

    ! SET DIMS USING MY S
    lb1 = lbound(mys,1)
    ub1 = ubound(mys,1)

    ! CHECK MY T
    if(lb1 /= lbound(myt,1)) CALL FATAL_ERROR("EQS_OF_STATE:: Dimension mismatch!")
    if(ub1 /= ubound(myt,1)) CALL FATAL_ERROR("EQS_OF_STATE:: Dimension mismatch!")

    ! CHECK MY P
    if(lb1 /= lbound(myp,1)) CALL FATAL_ERROR("EQS_OF_STATE:: Dimension mismatch!")
    if(ub1 /= ubound(myp,1)) CALL FATAL_ERROR("EQS_OF_STATE:: Dimension mismatch!")

    ! CHECK MY RHO
    if(lb1 /= lbound(myrho,1)) CALL FATAL_ERROR("EQS_OF_STATE:: Dimension mismatch!")
    if(ub1 /= ubound(myrho,1)) CALL FATAL_ERROR("EQS_OF_STATE:: Dimension mismatch!")

    DO K=lb1,ub1
       
       PT  = THETA(MYS(K),MYT(K),MYP(K),PREF)
       SVA = SVAN(MYS(K),PT,MYP(K),SIGMA)
       MYRHO(K) = SIGMA*1.e-3_SP
    END DO

    !----------------TRANSFORM TO FACE CENTER--------------------------------------

  END SUBROUTINE FOFONOFF_MILLARD_1D
  !==============================================================================|

  ! GENERIC CALL FOR FVCOM
  !==============================================================================|
  SUBROUTINE DENS1                

    !------------------------------------------------------------------------------|
    IMPLICIT NONE
    INTEGER :: K
    REAL(SP), PARAMETER ::PR = 0.0_SP
    REAL(SP), DIMENSION(0:MT,1:KB) :: RZU
    !==============================================================================|
    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START: DENS1"

    ! The thickness of the water column
    ! Is not the depth realtive to z=0
    DO K=1,KBM1
       RZU(:,K) = -GRAV_N*1.025_SP*(ZZ(:,K)*D(:))*0.1_SP
    END DO

    CALL FOFONOFF_MILLARD_2D(S1,T1,RZU,PR,RHO1)
    RHO1(:,KB)=0.0_SP
    RHO1(0,:)=0.0_SP


    !----------------TRANSFORM TO FACE CENTER--------------------------------------

    CALL N2E3D(RHO1,RHO)

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END: DENS1"
    RETURN
  END SUBROUTINE DENS1
  !==============================================================================|

  !==============================================================================|
  !     COMPUTE DENSITY USING SALINITY AND POTENTIAL TEMP                        |
  !	APPEARS TO BE BASED ON THE ROMS CODE?			                 | 
  !  CALCULATES: RHO1(M) DENSITY AT NODES			                 |
  !  CALCULATES: RHO (N) DENSITY AT ELEMENTS			                 |
  !==============================================================================|

  SUBROUTINE DENS2_2D(MYS,MYT,MYRHO)               

    !==============================================================================|
    IMPLICIT NONE
    REAL(SP), INTENT(IN), DIMENSION(:,:) :: MYS,MYT
    REAL(SP), INTENT(INOUT), DIMENSION(:,:):: MYRHO


    INTEGER :: I,K,ub1,lb1,ub2,lb2
    !==============================================================================|

    ! SET DIMS USING MY S
    lb1 = lbound(MYS,1)
    ub1 = ubound(MYS,1)

    lb2 = lbound(MYS,2)
    ub2 = ubound(MYS,2)

    ! CHECK MY T
    if(lb1 /= lbound(myt,1)) CALL FATAL_ERROR("EQS_OF_STATE:: Dimension mismatch!")
    if(ub1 /= ubound(myt,1)) CALL FATAL_ERROR("EQS_OF_STATE:: Dimension mismatch!")

    if(lb2 /= lbound(myt,2)) CALL FATAL_ERROR("EQS_OF_STATE:: Dimension mismatch!")
    if(ub2 /= ubound(myt,2)) CALL FATAL_ERROR("EQS_OF_STATE:: Dimension mismatch!")

    ! CHECK MY RHO
    if(lb1 /= lbound(myrho,1)) CALL FATAL_ERROR("EQS_OF_STATE:: Dimension mismatch!")
    if(ub1 /= ubound(myrho,1)) CALL FATAL_ERROR("EQS_OF_STATE:: Dimension mismatch!")

    if(lb2 /= lbound(myrho,2)) CALL FATAL_ERROR("EQS_OF_STATE:: Dimension mismatch!")
    if(ub2 /= ubound(myrho,2)) CALL FATAL_ERROR("EQS_OF_STATE:: Dimension mismatch!")

    
    !
    !  CALCULATE DENSITY FROM EQUATION OF STATE
    !
    MYRHO =MYS*MYS*MYS*6.76786136E-6_SP &
         & - MYS*MYS*4.8249614E-4_SP &
         & + MYS*8.14876577E-1_SP &
         & - 0.22584586E0_SP

    MYRHO = MYRHO * &
         & ( MYT*MYT*MYT*1.667E-8_SP &
         & - MYT*MYT*8.164E-7_SP &
         & + MYT*1.803E-5_SP &
         & )
    
    MYRHO = MYRHO+1.0_SP &
         & - MYT*MYT*MYT*1.0843E-6_SP &
         & + MYT*MYT*9.8185E-5_SP &
         & - MYT*4.786E-3_SP
    
    MYRHO = MYRHO* &
         & ( MYS*MYS*MYS*6.76786136E-6_SP &
         & - MYS*MYS*4.8249614E-4_SP &
         & + MYS*8.14876577E-1_SP &
         & +3.895414E-2_SP &
         & )
    
    MYRHO = MYRHO &
         & - (MYT-3.98_SP) *(MYT-3.98_SP) * (MYT+283.0_SP) / (503.57_SP*(MYT+67.26_SP))
    
    !
    !  CALCULATE RHO
    !
    MYRHO =  MYRHO*1.e-3_SP

  END SUBROUTINE DENS2_2D
  !==============================================================================|
  SUBROUTINE DENS2_1D(MYS,MYT,MYRHO)               

    !==============================================================================|
    IMPLICIT NONE
    REAL(SP), INTENT(IN), DIMENSION(:) :: MYS,MYT
    REAL(SP), INTENT(INOUT), DIMENSION(:):: MYRHO

    INTEGER :: I,K,ub1,lb1,ub2,lb2
    !==============================================================================|

    ! SET DIMS USING MY S
    lb1 = lbound(MYS,1)
    ub1 = ubound(MYS,1)

     ! CHECK MY T
    if(lb1 /= lbound(myt,1)) CALL FATAL_ERROR("EQS_OF_STATE:: Dimension mismatch!")
    if(ub1 /= ubound(myt,1)) CALL FATAL_ERROR("EQS_OF_STATE:: Dimension mismatch!")

    ! CHECK MY RHO
    if(lb1 /= lbound(myrho,1)) CALL FATAL_ERROR("EQS_OF_STATE:: Dimension mismatch!")
    if(ub1 /= ubound(myrho,1)) CALL FATAL_ERROR("EQS_OF_STATE:: Dimension mismatch!")
    
    !
    !  CALCULATE DENSITY FROM EQUATION OF STATE
    !
    MYRHO =MYS*MYS*MYS*6.76786136E-6_SP &
         & - MYS*MYS*4.8249614E-4_SP &
         & + MYS*8.14876577E-1_SP &
         & - 0.22584586E0_SP

    MYRHO = MYRHO * &
         & ( MYT*MYT*MYT*1.667E-8_SP &
         & - MYT*MYT*8.164E-7_SP &
         & + MYT*1.803E-5_SP &
         & )
    
    MYRHO = MYRHO+1.0_SP &
         & - MYT*MYT*MYT*1.0843E-6_SP &
         & + MYT*MYT*9.8185E-5_SP &
         & - MYT*4.786E-3_SP
    
    MYRHO = MYRHO* &
         & ( MYS*MYS*MYS*6.76786136E-6_SP &
         & - MYS*MYS*4.8249614E-4_SP &
         & + MYS*8.14876577E-1_SP &
         & +3.895414E-2_SP &
         & )
    
    MYRHO = MYRHO &
         & - (MYT-3.98_SP) *(MYT-3.98_SP) * (MYT+283.0_SP) / (503.57_SP*(MYT+67.26_SP))
    
    !
    !  CALCULATE RHO
    !
    MYRHO =  MYRHO*1.e-3_SP

  END SUBROUTINE DENS2_1D
  !==============================================================================|
  SUBROUTINE DENS2               

    !==============================================================================|
    IMPLICIT NONE
    !==============================================================================|
    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START: DENS2"

    !
    !  CALCULATE DENSITY FROM EQUATION OF STATE
    !
    CALL DENS2_2D(S1,T1,RHO1)
    RHO1(:,KB)=0.0_SP
    RHO1(0,:)=0.0_SP


    !
    !  AVERAGE FROM NODES TO FACE CENTERS
    !
    CALL N2E3D(RHO1,RHO)

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END: DENS2"
    RETURN
  END SUBROUTINE DENS2
  !==============================================================================|

  !==============================================================================|
  !     COMPUTE IN SITU DENSITY - 1000  USING SALINITY, POTENTIAL TEMP,          |
  !     AND PRESSURE FROM A POLYNOMIAL EXPRESSION (JACKETT & MCDOUGALL,          |
  !     1995). IT ASSUMES  NO  PRESSURE  VARIATION  ALONG GEOPOTENTIAL           |
  !     SURFACES, THAT IS, DEPTH (METERS; NEGATIVE) AND PRESSURE (DBAR           |
  !     ASSUMED NEGATIVE HERE) ARE INTERCHANGEABLE.                              |
  !                                                                              |
  !     check Values: (T=3 C, S=35.5 PSU, Z=-5000 m)                             |
  !        RHOF  = 1050.3639165364     (kg/m3)                                   |
  !        DEN1  = 1028.2845117925     (kg/m3)                                   |
  !                                                                              |
  !  Reference:                                                                  |
  !                                                                              |
  !  Jackett, D. R. and T. J. McDougall, 1995, Minimal Adjustment of             |
  !    Hydrostatic Profiles to Achieve Static Stability, J. of Atmos.            |
  !    and Oceanic Techn., vol. 12, pp. 381-389.                                 |
  !									       | 
  !    CALCULATES: RHO1(M) DENSITY AT NODES		                       |
  !    CALCULATES: RHO (N) DENSITY AT ELEMENTS				       |
  !==============================================================================|
  SUBROUTINE JACKET_MCDOUGALL_2D(MYS,MYT,MYP,MYRHO)      

    !==============================================================================|
    IMPLICIT NONE
    REAL(SP), INTENT(IN), DIMENSION(:,:) :: MYS,MYT,MYP
    REAL(SP), INTENT(INOUT), DIMENSION(:,:):: MYRHO

    REAL(SP) :: TF,SF,sqrtSF,PBAR,TEMP(10),BULK,BULK0,BULK1,BULK2
    INTEGER :: I,K,ub1,lb1,ub2,lb2
    !==============================================================================|
    !  Polynomial  expansion  coefficients for the computation of in situ          |
    !  density  via  the  nonlinear  equation of state  for seawater as a          |
    !  function of potential temperature, salinity, and pressure (Jackett          |
    !  and McDougall, 1995).                                                       |
    REAL(SP), PARAMETER :: A00 = +1.965933e+04_SP
    REAL(SP), PARAMETER :: A01 = +1.444304e+02_SP
    REAL(SP), PARAMETER :: A02 = -1.706103e+00_SP
    REAL(SP), PARAMETER :: A03 = +9.648704e-03_SP
    REAL(SP), PARAMETER :: A04 = -4.190253e-05_SP
    REAL(SP), PARAMETER :: B00 = +5.284855e+01_SP
    REAL(SP), PARAMETER :: B01 = -3.101089e-01_SP
    REAL(SP), PARAMETER :: B02 = +6.283263e-03_SP
    REAL(SP), PARAMETER :: B03 = -5.084188e-05_SP
    REAL(SP), PARAMETER :: D00 = +3.886640e-01_SP
    REAL(SP), PARAMETER :: D01 = +9.085835e-03_SP
    REAL(SP), PARAMETER :: D02 = -4.619924e-04_SP
    REAL(SP), PARAMETER :: E00 = +3.186519e+00_SP
    REAL(SP), PARAMETER :: E01 = +2.212276e-02_SP
    REAL(SP), PARAMETER :: E02 = -2.984642e-04_SP
    REAL(SP), PARAMETER :: E03 = +1.956415e-06_SP
    REAL(SP), PARAMETER :: F00 = +6.704388e-03_SP
    REAL(SP), PARAMETER :: F01 = -1.847318e-04_SP
    REAL(SP), PARAMETER :: F02 = +2.059331e-07_SP
    REAL(SP), PARAMETER :: G00 = +1.480266e-04_SP
    REAL(SP), PARAMETER :: G01 = +2.102898e-04_SP
    REAL(SP), PARAMETER :: G02 = -1.202016e-05_SP
    REAL(SP), PARAMETER :: G03 = +1.394680e-07_SP
    REAL(SP), PARAMETER :: H00 = -2.040237e-06_SP
    REAL(SP), PARAMETER :: H01 = +6.128773e-08_SP
    REAL(SP), PARAMETER :: H02 = +6.207323e-10_SP

    REAL(SP), PARAMETER :: Q00 = +9.99842594e+02_SP
    REAL(SP), PARAMETER :: Q01 = +6.793952e-02_SP
    REAL(SP), PARAMETER :: Q02 = -9.095290e-03_SP
    REAL(SP), PARAMETER :: Q03 = +1.001685e-04_SP
    REAL(SP), PARAMETER :: Q04 = -1.120083e-06_SP
    REAL(SP), PARAMETER :: Q05 = +6.536332e-09_SP
    REAL(SP), PARAMETER :: U00 = +8.24493e-01_SP
    REAL(SP), PARAMETER :: U01 = -4.08990e-03_SP
    REAL(SP), PARAMETER :: U02 = +7.64380e-05_SP
    REAL(SP), PARAMETER :: U03 = -8.24670e-07_SP
    REAL(SP), PARAMETER :: U04 = +5.38750e-09_SP
    REAL(SP), PARAMETER :: V00 = -5.72466e-03_SP
    REAL(SP), PARAMETER :: V01 = +1.02270e-04_SP
    REAL(SP), PARAMETER :: V02 = -1.65460e-06_SP
    REAL(SP), PARAMETER :: W00 = +4.8314e-04_SP


    ! SET DIMS USING MY S
    lb1 = lbound(mys,1)
    ub1 = ubound(mys,1)

    lb2 = lbound(mys,2)
    ub2 = ubound(mys,2)

    ! CHECK MY T
    if(lb1 /= lbound(myt,1)) CALL FATAL_ERROR("EQS_OF_STATE:: Dimension mismatch!")
    if(ub1 /= ubound(myt,1)) CALL FATAL_ERROR("EQS_OF_STATE:: Dimension mismatch!")

    if(lb2 /= lbound(myt,2)) CALL FATAL_ERROR("EQS_OF_STATE:: Dimension mismatch!")
    if(ub2 /= ubound(myt,2)) CALL FATAL_ERROR("EQS_OF_STATE:: Dimension mismatch!")

    ! CHECK MY P
    if(lb1 /= lbound(myp,1)) CALL FATAL_ERROR("EQS_OF_STATE:: Dimension mismatch!")
    if(ub1 /= ubound(myp,1)) CALL FATAL_ERROR("EQS_OF_STATE:: Dimension mismatch!")

    if(lb2 /= lbound(myp,2)) CALL FATAL_ERROR("EQS_OF_STATE:: Dimension mismatch!")
    if(ub2 /= ubound(myp,2)) CALL FATAL_ERROR("EQS_OF_STATE:: Dimension mismatch!")

    ! CHECK MY RHO
    if(lb1 /= lbound(myrho,1)) CALL FATAL_ERROR("EQS_OF_STATE:: Dimension mismatch!")
    if(ub1 /= ubound(myrho,1)) CALL FATAL_ERROR("EQS_OF_STATE:: Dimension mismatch!")

    if(lb2 /= lbound(myrho,2)) CALL FATAL_ERROR("EQS_OF_STATE:: Dimension mismatch!")
    if(ub2 /= ubound(myrho,2)) CALL FATAL_ERROR("EQS_OF_STATE:: Dimension mismatch!")
 


    !
    !  CALCULATE DENSITY FROM EQUATION OF STATE
    !
    DO I=lb1,ub1
       DO K=lb2,ub2
          TF = MYT(I,K)
          SF = MYS(I,K)
          sqrtSF = sqrt(SF)

          PBAR = MYP(I,K)

          !  Compute density (kg/m3) at standard one atmosphere pressure
          TEMP(1)=Q00+TF*(Q01+TF*(Q02+TF*(Q03+TF*(Q04+TF*Q05))))
          TEMP(2)=U00+TF*(U01+TF*(U02+TF*(U03+TF*U04)))
          TEMP(3)=V00+TF*(V01+TF*V02)
          MYRHO(I,K)=TEMP(1)+SF*(TEMP(2)+sqrtSF*TEMP(3)+SF*W00)

          !  Compute secant bulk modulus (BULK = BULK0 + BULK1*PBAR + BULK2*PBAR*PBAR)
          TEMP(4)=A00+TF*(A01+TF*(A02+TF*(A03+TF*A04)))
          TEMP(5)=B00+TF*(B01+TF*(B02+TF*B03))
          TEMP(6)=D00+TF*(D01+TF*D02)
          TEMP(7)=E00+TF*(E01+TF*(E02+TF*E03))
          TEMP(8)=F00+TF*(F01+TF*F02)
          TEMP(9)=G01+TF*(G02+TF*G03)
          TEMP(10)=H00+TF*(H01+TF*H02)

          BULK0=TEMP(4)+SF*(TEMP(5)+sqrtSF*TEMP(6))
          BULK1=TEMP(7)+SF*(TEMP(8)+sqrtSF*G00)
          BULK2=TEMP(9)+SF*TEMP(10)
          BULK = BULK0 + PBAR * (BULK1 + PBAR * BULK2)

          !  Compute "in situ" density anomaly (kg/m3)
          MYRHO(I,K)=(MYRHO(I,K)*BULK)/(BULK-PBAR)
          MYRHO(I,K)= MYRHO(I,K)-1000.0_SP
       END DO
    END DO
 
    !
    !  CALCULATE RHO1
    !
    MYRHO =  MYRHO*1.e-3_SP
    

  END SUBROUTINE JACKET_MCDOUGALL_2D
!==============================================================================|
  SUBROUTINE JACKET_MCDOUGALL_1D(MYS,MYT,MYP,MYRHO)      

    !==============================================================================|
    IMPLICIT NONE
    REAL(SP), INTENT(IN), DIMENSION(:) :: MYS,MYT,MYP
    REAL(SP), INTENT(INOUT), DIMENSION(:):: MYRHO

    REAL(SP) :: TF,SF,sqrtSF,PBAR,TEMP(10),BULK,BULK0,BULK1,BULK2
    INTEGER :: K,ub1,lb1
    !==============================================================================|
    !  Polynomial  expansion  coefficients for the computation of in situ          |
    !  density  via  the  nonlinear  equation of state  for seawater as a          |
    !  function of potential temperature, salinity, and pressure (Jackett          |
    !  and McDougall, 1995).                                                       |
    REAL(SP), PARAMETER :: A00 = +1.965933e+04_SP
    REAL(SP), PARAMETER :: A01 = +1.444304e+02_SP
    REAL(SP), PARAMETER :: A02 = -1.706103e+00_SP
    REAL(SP), PARAMETER :: A03 = +9.648704e-03_SP
    REAL(SP), PARAMETER :: A04 = -4.190253e-05_SP
    REAL(SP), PARAMETER :: B00 = +5.284855e+01_SP
    REAL(SP), PARAMETER :: B01 = -3.101089e-01_SP
    REAL(SP), PARAMETER :: B02 = +6.283263e-03_SP
    REAL(SP), PARAMETER :: B03 = -5.084188e-05_SP
    REAL(SP), PARAMETER :: D00 = +3.886640e-01_SP
    REAL(SP), PARAMETER :: D01 = +9.085835e-03_SP
    REAL(SP), PARAMETER :: D02 = -4.619924e-04_SP
    REAL(SP), PARAMETER :: E00 = +3.186519e+00_SP
    REAL(SP), PARAMETER :: E01 = +2.212276e-02_SP
    REAL(SP), PARAMETER :: E02 = -2.984642e-04_SP
    REAL(SP), PARAMETER :: E03 = +1.956415e-06_SP
    REAL(SP), PARAMETER :: F00 = +6.704388e-03_SP
    REAL(SP), PARAMETER :: F01 = -1.847318e-04_SP
    REAL(SP), PARAMETER :: F02 = +2.059331e-07_SP
    REAL(SP), PARAMETER :: G00 = +1.480266e-04_SP
    REAL(SP), PARAMETER :: G01 = +2.102898e-04_SP
    REAL(SP), PARAMETER :: G02 = -1.202016e-05_SP
    REAL(SP), PARAMETER :: G03 = +1.394680e-07_SP
    REAL(SP), PARAMETER :: H00 = -2.040237e-06_SP
    REAL(SP), PARAMETER :: H01 = +6.128773e-08_SP
    REAL(SP), PARAMETER :: H02 = +6.207323e-10_SP

    REAL(SP), PARAMETER :: Q00 = +9.99842594e+02_SP
    REAL(SP), PARAMETER :: Q01 = +6.793952e-02_SP
    REAL(SP), PARAMETER :: Q02 = -9.095290e-03_SP
    REAL(SP), PARAMETER :: Q03 = +1.001685e-04_SP
    REAL(SP), PARAMETER :: Q04 = -1.120083e-06_SP
    REAL(SP), PARAMETER :: Q05 = +6.536332e-09_SP
    REAL(SP), PARAMETER :: U00 = +8.24493e-01_SP
    REAL(SP), PARAMETER :: U01 = -4.08990e-03_SP
    REAL(SP), PARAMETER :: U02 = +7.64380e-05_SP
    REAL(SP), PARAMETER :: U03 = -8.24670e-07_SP
    REAL(SP), PARAMETER :: U04 = +5.38750e-09_SP
    REAL(SP), PARAMETER :: V00 = -5.72466e-03_SP
    REAL(SP), PARAMETER :: V01 = +1.02270e-04_SP
    REAL(SP), PARAMETER :: V02 = -1.65460e-06_SP
    REAL(SP), PARAMETER :: W00 = +4.8314e-04_SP


    ! SET DIMS USING MY S
    lb1 = lbound(mys,1)
    ub1 = ubound(mys,1)

    ! CHECK MY T
    if(lb1 /= lbound(myt,1)) CALL FATAL_ERROR("EQS_OF_STATE:: Dimension mismatch!")
    if(ub1 /= ubound(myt,1)) CALL FATAL_ERROR("EQS_OF_STATE:: Dimension mismatch!")

    ! CHECK MY P
    if(lb1 /= lbound(myp,1)) CALL FATAL_ERROR("EQS_OF_STATE:: Dimension mismatch!")
    if(ub1 /= ubound(myp,1)) CALL FATAL_ERROR("EQS_OF_STATE:: Dimension mismatch!")

    ! CHECK MY RHO
    if(lb1 /= lbound(myrho,1)) CALL FATAL_ERROR("EQS_OF_STATE:: Dimension mismatch!")
    if(ub1 /= ubound(myrho,1)) CALL FATAL_ERROR("EQS_OF_STATE:: Dimension mismatch!")

    !
    !  CALCULATE DENSITY FROM EQUATION OF STATE
    !
    DO K=lb1,ub1
       TF = MYT(K)
       SF = MYS(K)
       sqrtSF = sqrt(SF)
       
       PBAR = MYP(K)
       
       !  Compute density (kg/m3) at standard one atmosphere pressure
       TEMP(1)=Q00+TF*(Q01+TF*(Q02+TF*(Q03+TF*(Q04+TF*Q05))))
       TEMP(2)=U00+TF*(U01+TF*(U02+TF*(U03+TF*U04)))
       TEMP(3)=V00+TF*(V01+TF*V02)
       MYRHO(K)=TEMP(1)+SF*(TEMP(2)+sqrtSF*TEMP(3)+SF*W00)
       
       !  Compute secant bulk modulus (BULK = BULK0 + BULK1*PBAR + BULK2*PBAR*PBAR)
       TEMP(4)=A00+TF*(A01+TF*(A02+TF*(A03+TF*A04)))
       TEMP(5)=B00+TF*(B01+TF*(B02+TF*B03))
       TEMP(6)=D00+TF*(D01+TF*D02)
       TEMP(7)=E00+TF*(E01+TF*(E02+TF*E03))
       TEMP(8)=F00+TF*(F01+TF*F02)
       TEMP(9)=G01+TF*(G02+TF*G03)
       TEMP(10)=H00+TF*(H01+TF*H02)
       
       BULK0=TEMP(4)+SF*(TEMP(5)+sqrtSF*TEMP(6))
       BULK1=TEMP(7)+SF*(TEMP(8)+sqrtSF*G00)
       BULK2=TEMP(9)+SF*TEMP(10)
       BULK = BULK0 + PBAR * (BULK1 + PBAR * BULK2)
       
       !  Compute "in situ" density anomaly (kg/m3)
       MYRHO(K)=(MYRHO(K)*BULK)/(BULK-PBAR)
       MYRHO(K)= MYRHO(K)-1000.0_SP
    END DO
     
    !
    !  CALCULATE RHO1
    !
    MYRHO =  MYRHO*1.e-3_SP
    

  END SUBROUTINE JACKET_MCDOUGALL_1D
!==============================================================================|
  SUBROUTINE DENS3        
!==============================================================================|
    IMPLICIT NONE
    REAL(SP), DIMENSION(0:MT,KB) :: MYP
    INTEGER :: K
    
    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START: DENS3"

    DO K=1,KBM1
       MYP(:,K) = -GRAV_N(:)*1.025_SP*ZZ(:,K)*D(:) *0.01_SP
    END DO

    CALL JACKET_MCDOUGALL_2D(S1,T1,MYP,RHO1)
    RHO1(:,KB)=0.0_SP
    RHO1(0,:)=0.0_SP

    !
    !  AVERAGE FROM NODES TO FACE CENTERS
    !
    CALL N2E3D(RHO1,RHO)

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END: DENS3"
    RETURN
  END SUBROUTINE DENS3
  !==============================================================================!
  FUNCTION SVAN(S4,T4,P04,SIGMA)
    !==============================================================================!
    ! specific volume anomaly (steric anomaly) based on 1980 equation              |
    ! of state for seawater and 1978 practerical salinity scale.                   |
    ! references:                                                                  |
    ! millero, et al (1980) deep-sea res.,27a,255-264                              |
    ! millero and poisson 1981,deep-sea res.,28a pp 625-629.                       |
    ! both above references are also found in unesco report 38 (1981)              |
    !                                                                              |
    ! units:                                                                       |
    !       pressure        p04       decibars                                     |
    !       temperature     t4        deg celsius (ipts-68)                        |
    !       salinity        s4        (ipss-78)                                    |
    !       spec. vol. ana. svan     m**3/kg *1.0e-8                               |
    !       density ana.    sigma    kg/m**3                                       |
    !                                                                              |
    ! check value: svan=981.3021 e-8 m**3/kg. for s = 40 (ipss-78),                |
    ! t = 40 deg c, p0= 10000 decibars.                                            |
    ! check value: sigma = 59.82037  kg/m**3. for s = 40 (ipss-78) ,               |
    ! t = 40 deg c, p0= 10000 decibars.                                            |
    !==============================================================================!


    USE MOD_PREC
    IMPLICIT NONE
    REAL(SP) :: SVAN
    REAL(SP), INTENT(IN)  :: S4,T4,P04
    REAL(SP), INTENT(OUT) :: SIGMA
    REAL(SP) P4,SIG,SR,RR1,RR2,RR3,V350P,DK
    REAL(SP) A4,B4,C4,D4,E4,AA1,BB1,AW,BW,KK,K0,KW,K35,SVA
    REAL(SP) GAM,PK,DVAN,DR35P

    REAL(SP), PARAMETER :: R3500 = 1028.1063_SP
    REAL(SP), PARAMETER :: RR4   = 4.8314E-4_SP
    REAL(SP), PARAMETER :: DR350 = 28.106331_SP


    !   rr4 is refered to as  c  in millero and poisson 1981
    ! convert pressure to bars and take square root salinity.

    P4=P04/10.0_SP
    SR = SQRT(ABS(S4))

    ! pure water density at atmospheric pressure
    !   bigg p.h.,(1967) br. j. applied physics 8 pp 521-537.
    !

    RR1=((((6.536332E-9_SP*T4-1.120083E-6_SP)*T4+1.001685E-4_SP)*T4 &
         -9.095290E-3_SP)*T4+6.793952E-2_SP)*T4-28.263737_SP


    ! seawater density atm press.
    !  coefficients involving salinity
    !  rr2 = a   in notation of millero and poisson 1981

    RR2=(((5.3875E-9_SP*T4-8.2467E-7_SP)*T4+7.6438E-5_SP)*T4-4.0899E-3_SP)*T4 &
         +8.24493E-1_SP

    !  rr3 = b4  in notation of millero and poisson 1981

    RR3=(-1.6546E-6_SP*T4+1.0227E-4_SP)*T4-5.72466E-3_SP

    !  international one-atmosphere equation of state of seawater

    SIG=(RR4*S4+RR3*SR+RR2)*S4+RR1

    ! specific volume at atmospheric pressure

    V350P = 1.0_SP/R3500
    SVA = -SIG*V350P/(R3500+SIG)
    SIGMA=SIG+DR350

    !  scale specific vol. anamoly to normally reported units

    SVAN=SVA*1.0E+8_SP
    IF(P4 == 0.0_SP) RETURN

    !-------------------------------------------------------------|
    !    new high pressure equation of sate for seawater          |
    !                                                             |
    !        millero, el al., 1980 dsr 27a, pp 255-264            |
    !        constant notation follows article                    |
    !-------------------------------------------------------------|
    ! compute compression terms

    E4  = (9.1697E-10*T4+2.0816E-8_SP)*T4-9.9348E-7_SP
    BW  = (5.2787E-8_SP*T4-6.12293E-6_SP)*T4+3.47718E-5_SP
    B4  = BW + E4*S4

    D4  = 1.91075E-4
    C4  = (-1.6078E-6_SP*T4-1.0981E-5_SP)*T4+2.2838E-3_SP
    AW  = ((-5.77905E-7_SP*T4+1.16092E-4_SP)*T4+1.43713E-3_SP)*T4 &
         -0.1194975_SP
    A4  = (D4*SR + C4)*S4 + AW

    BB1 = (-5.3009E-4_SP*T4+1.6483E-2_SP)*T4+7.944E-2_SP
    AA1 = ((-6.1670E-5_SP*T4+1.09987E-2_SP)*T4-0.603459_SP)*T4+54.6746
    KW  = (((-5.155288E-5_SP*T4+1.360477E-2_SP)*T4-2.327105_SP)*T4 &
         +148.4206_SP)*T4-1930.06_SP
    K0  = (BB1*SR + AA1)*S4 + KW

    ! evaluate pressure polynomial
    !-----------------------------------------------------|
    !   k equals the secant bulk modulus of seawater      |
    !   dk=k(s,t,p)-k(35,0,p)                             |
    !   k35=k(35,0,p)                                     |
    !-----------------------------------------------------|

    DK = (B4*P4 + A4)*P4 + K0
    K35  = (5.03217E-5_SP*P4+3.359406_SP)*P4+21582.27_SP
    GAM=P4/K35
    PK = 1.0_SP - GAM
    SVA = SVA*PK + (V350P+SVA)*P4*DK/(K35*(K35+DK))

    !  scale specific vol. anamoly to normally reported units

    SVAN=SVA*1.0E+8_SP
    V350P = V350P*PK

    !----------------------------------------------------------|
    ! compute density anamoly with respect to 1000.0 kg/m**3   |
    !  1) dr350: density anamoly at 35 (ipss-78),              |
    !                               0 deg. c and 0 decibars    |
    !  2) dr35p: density anamoly at 35 (ipss-78),              |
    !                               0 deg. c, pres. variation  |
    !  3) dvan : density anamoly variations involving specific |
    !            volume anamoly                                |
    !                                                          |
    ! check values: sigma = 59.82037 kg/m**3                   |
    ! for s = 40 (ipss-78), t = 40 deg c, p0= 10000 decibars.  |
    !----------------------------------------------------------|

    DR35P=GAM/V350P
    DVAN=SVA/(V350P*(V350P+SVA))
    SIGMA=DR350+DR35P-DVAN

    RETURN
  END FUNCTION SVAN
  !==============================================================================!
  !==============================================================================!

  !==============================================================================|
  FUNCTION THETA(S4,T04,P04,PR)
    !==============================================================================|
    ! to compute local potential temperature at pr using                           |
    ! bryden 1973 polynomial for adiabatic lapse rate and                          |
    ! runge-kutta 4th order integration algorithm.                                 |
    ! ref: bryden,h.,1973,deep-sea res.,20,401-408;                                |
    ! fofonoff,n.,1977,deep-sea res.,24,489-491                                    |
    !                                                                              |
    ! units:                                                                       |
    !       pressure        p04       decibars                                     |
    !       temperature     t04       deg celsius (ipts-68)                        |
    !       salinity         s4        (ipss-78)                                   |
    !       reference prs    pr       decibars                                     |
    !       potential tmp.  theta     deg celsius                                  |
    ! checkvalue:                                                                  |
    !             theta= 36.89073 c,s=40 (ipss-78),                                |
    !             t0=40 deg c,p0=10000 decibars,pr=0 decibars                      |
    !                                                                              |
    !                                                                              |
    ! set up intermediate temperature and pressure variables.                      |
    !==============================================================================|

    USE MOD_PREC
    IMPLICIT NONE
    REAL(SP) :: THETA
    REAL(SP), INTENT(IN) :: S4,T04,P04,PR
    REAL(SP) :: P4,T4,H4,Q4,XK

    
    !==============================================================================|

    P4 = P04
    T4 = T04
    H4 = PR - P4
    XK = H4*ATG(S4,T4,P4)
    T4 = T4 + .5_SP*XK
    Q4 = XK
    P4 = P4 + 0.5_SP*H4
    XK = H4*ATG(S4,T4,P4)
    T4 = T4 + 0.29289322_SP*(XK-Q4)
    Q4 = 0.58578644_SP*XK + 0.121320344_SP*Q4
    XK = H4*ATG(S4,T4,P4)
    T4 = T4 + 1.707106781_SP*(XK-Q4)
    Q4 = 3.414213562_SP*XK - 4.121320344_SP*Q4
    P4 = P4 + 0.5_SP*H4
    XK = H4*ATG(S4,T4,P4)
    THETA = T4 + (XK-2.0_SP*Q4)/6.0_SP

    RETURN
  END FUNCTION THETA
  !==============================================================================|

  !==============================================================================|
  ! adiabatic temperature gradient deg c per decibar    			       |
  ! ref: bryden, h., 1973,deep-sea res.,20,401-408                               |
  !                                                                              |
  ! units:                                                                       |
  !       pressure        P4        decibars                                     |
  !       temperature     T4        deg celsius(ipts-68)                         |
  !       salinity        s4        (ipss-78)                                    |
  !       adiabatic      atg        deg. c/decibar                               |
  ! checkvalue: atg=3.255976e-4 c/dbar for s=40 (ipss-78),                       |
  ! t=40 deg c,p0=10000 decibars                                                 |
  !==============================================================================|

  REAL(SP) FUNCTION ATG(S4,T4,P4)
    USE MOD_PREC

    !------------------------------------------------------------------------------|

    IMPLICIT NONE
    REAL(SP), INTENT(IN) :: S4,T4,P4
    REAL(SP)  :: DS

    !==============================================================================|

    DS  = S4 - 35.0_SP
    ATG = (((-2.1687e-16_SP*T4+1.8676e-14_SP)*T4-4.6206e-13_SP)*P4 &
         +((2.7759e-12_SP*T4-1.1351e-10_SP)*DS+((-5.4481e-14_SP*T4 &
         +8.733e-12_SP)*T4-6.7795e-10_SP)*T4+1.8741e-8_SP))*P4 &
         +(-4.2393e-8_SP*T4+1.8932e-6_SP)*DS &
         +((6.6228e-10_SP*T4-6.836e-8_SP)*T4+8.5258e-6_SP)*T4+3.5803e-5_SP

    RETURN
  END FUNCTION ATG
  !==============================================================================|

END MODULE EQS_OF_STATE
