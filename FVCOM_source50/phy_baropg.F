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
!   This subroutine is used to calculate the baroclinic pressure               !
!   gradient in the standarded z-levels. The water column is divided           !
!   into 600 standard levels, and pressure gradient is then determined         !
!   at each level. The resulting pressure gradients are converted back         !
!   to sigma-levels through vertical interpolation approach.                   !
!==============================================================================|

   SUBROUTINE PHY_BAROPG          

!------------------------------------------------------------------------------|

   USE ALL_VARS
   USE SINTER
   IMPLICIT NONE
   INTEGER, PARAMETER :: KBB=601
   INTEGER, PARAMETER :: KBBM1=KBB-1
   REAL(SP) RHOZ(M,KBBM1),RHOZTMP(KBBM1),PHY_Z(KBBM1)
   REAL(SP) RHOS(KBM1),SIGZTMP(KBM1),SIG_Z(M,KBM1)
   REAL(SP) PB1(0:KBBM1),PB2(0:KBBM1),PB3(0:KBBM1)
   REAL(SP) PBXZ(KBBM1),PBYZ(KBBM1)
   REAL(SP) PBXS(KBM1) ,PBYS(KBM1)
   REAL(SP) AREAX1,AREAX2,AREAX3,AREAY1,AREAY2,AREAY3
   REAL(SP) RHOZI1,RHOZI2,RHOZI3,TMP,TEMP,DELTZ
   REAL(SP) GRAV1,GRAV2,GRAV3
   INTEGER  I,K,J1,J2,J3,NTMP

   ! USE RAMP CALCULATED IN 'internal_step.F'

!--CALCULATE Z-LEVELS TO MAX DEPTH---------------------------------------------|
   
   DELTZ=HMAX/FLOAT(KBBM1)

   DO K=1,KBBM1
     PHY_Z(K)=(0.5_SP-FLOAT(K))*DELTZ
   END DO

!--LINEARLY INTERPOLATE TO OBTAIN DENSITY VALUES AT Z LEVELS-------------------|

   DO I=1,M
     DO K=1,KBM1
       SIG_Z(I,K)=ZZ(I,K)*DT(I)+ET(I)
       SIGZTMP(K)=SIG_Z(I,K)
       RHOS(K)=RHO1(I,K)
     END DO

     CALL SINTER_EXTRP_DOWN(SIGZTMP,RHOS,PHY_Z,RHOZTMP,KBM1,KBBM1)

     DO K=1,KBBM1
       RHOZ(I,K)=RHOZTMP(K)
     END DO
   END DO

   DO I=1,N
     J1=NV(I,1)
     J2=NV(I,2)
     J3=NV(I,3)
     NTMP=0
     PB1(0)=0.0_SP
     PB2(0)=0.0_SP
     PB3(0)=0.0_SP
     DO K=1,KBBM1
       TMP=FLOAT(K)*DELTZ
       IF((H(J1) < TMP.OR.H(J2) < TMP.OR.H(J3) < TMP)) THEN
         PB1(K)=0.0_SP
         PB2(K)=0.0_SP
         PB3(K)=0.0_SP
       ELSE
         RHOZI1=0.5_SP*(RHOZ(J2,K)+RHOZ(J3,K))
         RHOZI2=0.5_SP*(RHOZ(J3,K)+RHOZ(J1,K))
         RHOZI3=0.5_SP*(RHOZ(J1,K)+RHOZ(J2,K))
	 GRAV1 =0.5_SP*(GRAV_N(J2)+GRAV_N(J3))
	 GRAV2 =0.5_SP*(GRAV_N(J3)+GRAV_N(J1))
	 GRAV3 =0.5_SP*(GRAV_N(J1)+GRAV_N(J2))
         PB1(K)=PB1(K-1)+GRAV1*RHOZI1*DELTZ
         PB2(K)=PB2(K-1)+GRAV2*RHOZI2*DELTZ
         PB3(K)=PB3(K-1)+GRAV3*RHOZI3*DELTZ
         NTMP=NTMP+1
       END IF
     END DO
     AREAX1=(VY(J3)-VY(J2))*DELTZ
     AREAY1=(VX(J2)-VX(J3))*DELTZ
     AREAX2=(VY(J1)-VY(J3))*DELTZ
     AREAY2=(VX(J3)-VX(J1))*DELTZ
     AREAX3=(VY(J2)-VY(J1))*DELTZ
     AREAY3=(VX(J1)-VX(J2))*DELTZ
     DO K=1,KBBM1
       PBXZ(K)=AREAX1*PB1(K)+AREAX2*PB2(K)+AREAX3*PB3(K)
       PBYZ(K)=AREAY1*PB1(K)+AREAY2*PB2(K)+AREAY3*PB3(K)
       PBXZ(K)=PBXZ(K)/ART(I)/DELTZ
       PBYZ(K)=PBYZ(K)/ART(I)/DELTZ
     END DO

     DO K=1,KBM1
       J1=NV(I,1)
       J2=NV(I,2)
       J3=NV(I,3)
       SIGZTMP(K)=(SIG_Z(J1,K)+SIG_Z(J2,K)+SIG_Z(J3,K))/3.
     END DO

     IF(NTMP == 0) THEN
       DO K=1,KBM1
         PBXS(K)=0.0_SP
         PBYS(K)=0.0_SP
       END DO
     ELSE IF(NTMP == 1) THEN
       DO K=1,KBM1
         PBXS(K)=PBXZ(1)
         PBYS(K)=PBYZ(1)
       END DO
     ELSE
       CALL SINTER_EXTRP_DOWN(PHY_Z,PBXZ,SIGZTMP,PBXS,NTMP,KBM1)
       CALL SINTER_EXTRP_DOWN(PHY_Z,PBYZ,SIGZTMP,PBYS,NTMP,KBM1)
     END IF

     DO K=1,KBM1
       DRHOX(I,K)=-PBXS(K)*DT1(I)*DZ1(I,K)*ART(I)*RAMP
       DRHOY(I,K)=-PBYS(K)*DT1(I)*DZ1(I,K)*ART(I)*RAMP
     END DO
   END DO

   RETURN
   END SUBROUTINE PHY_BAROPG
!==============================================================================|
