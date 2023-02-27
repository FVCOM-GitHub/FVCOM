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

MODULE MOD_VECTOR_PROJECTION
  
  USE ALL_VARS
# if defined (MULTIPROCESSOR)
  USE MOD_PAR
# endif
  IMPLICIT NONE
     
CONTAINS
# if defined (PROJ) && !defined (SPHERICAL)

  !======================================================================
  ! Project the current velocity from Cartisian Coordiante to 
  ! Geo-referenced Coordinate.
  ! 
  ! When simulated in Cartisian Coordinate, the current velocities (u,v)
  ! point to the direction of x- and y- axes. By this projection, (u,v)
  ! can be converted/rotated to the velocities in Geo-referenced 
  ! Coordinate (u_NORTH, v_EAST), pointing the east and north, 
  ! respectively.
  ! 
  ! Thank Dr. Vitalii Sheremet for helping with the math.
  !
  ! In order to invoke this projection, the 'Flag_6 = -DPROJ' in make.inc
  ! needs to be turned on and the 'PREJECTION_REFERENCE' in your nml file
  ! needs to be set with a valid paramemters.
  !
  ! There are two parts for the projection.
  !   ---SETUP_UV_PROJECTION : calculate the paramters used in projection  
  !   --- UV_PROJECTION      : do the projection for (u,v) and (ua,va)
  !
  ! Siqi Li, 20221005
  !======================================================================
  SUBROUTINE SETUP_UV_PROJECTION

    INTEGER  :: J, I1, I2, I3
    REAL(SP) :: LON1, LON2, LON3, LAT1, LAT2, LAT3, X1, X2, X3, Y1, Y2, Y3, AREA

    DO J = 1, N
      I1 = NV(J,1);   I2 = NV(J,2);   I3 = NV(J,3)
      X1 = VX(I1);    X2 = VX(I2);    X3 = VX(I3)
      Y1 = VY(I1);    Y2 = VY(I2);    Y3 = VY(I3)
      LON1 = LON(I1); LON2 = LON(I2); LON3 = LON(I3)
      LAT1 = LAT(I1); LAT2 = LAT(I2); LAT3 = LAT(I3)

      AREA = (X2-X1)*(Y3-Y1) - (X3-X1)*(Y2-Y1)

      ACLON_X(J) =  REARTH * COSD(LATC(J)) * ((LON2-LON1)*(Y3-Y1)-(LON3-LON1)*(Y2-Y1)) / AREA*PI/180.
      ACLON_Y(J) = -REARTH * COSD(LATC(J)) * ((LON2-LON1)*(X3-X1)-(LON3-LON1)*(X2-X1)) / AREA*PI/180.
      ALAT_X(J)  =  REARTH                 * ((LAT2-LAT1)*(Y3-Y1)-(LAT3-LAT1)*(Y2-Y1)) / AREA*PI/180.
      ALAT_Y(J)  = -REARTH                 * ((LAT2-LAT1)*(X3-X1)-(LAT3-LAT1)*(X2-X1)) / AREA*PI/180.
    END DO

  END SUBROUTINE SETUP_UV_PROJECTION  

  SUBROUTINE UV_PROJECTION

    INTEGER  :: J, K

    DO J = 1, N
      UA_EAST(J) = ACLON_X(J)*UA(J) + ACLON_Y(J)*VA(J)
      VA_NORTH(J) =  ALAT_X(J)*UA(J) +  ALAT_Y(J)*VA(J)
      DO K = 1, KBM1
        U_EAST(J,K) = ACLON_X(J)*U(J,K) + ACLON_Y(J)*V(J,K)
        V_NORTH(J,K) =  ALAT_X(J)*U(J,K) +  ALAT_Y(J)*V(J,K)
      END DO
    END DO

# if defined (MULTIPROCESSOR)
    CALL AEXCHANGE(EC,MYID,NPROCS,U_EAST,V_NORTH)  
    CALL AEXCHANGE(EC,MYID,NPROCS,UA_EAST,VA_NORTH)
# endif

  END SUBROUTINE UV_PROJECTION
# endif


END MODULE MOD_VECTOR_PROJECTION
