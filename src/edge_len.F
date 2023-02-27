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

SUBROUTINE edge_len
  USE all_vars
  use mod_spherical
  use mod_northpole
  use mod_utils

  implicit none
  INTEGER  :: I,IP,I1,I2,IA,IB,J,J1
  REAL(SP) :: XI, YI, X11,X33,Y11,Y33

  INTEGER :: JTMP,J2

#  if defined (SPHERICAL)
  REAL(DP) :: XTMP1,XTMP
  REAL(DP) :: X1_DP,Y1_DP,X2_DP,Y2_DP,XII,YII
  REAL(DP) :: X11_TMP,Y11_TMP,X33_TMP,Y33_TMP
  REAL(DP) :: VX1_TMP,VY1_TMP,VX2_TMP,VY2_TMP
  REAL(DP) :: TXPI_A,TYPI_A
  REAL(DP) :: TXPI_B,TYPI_B
#  endif


  ! Distance between control volue edge and the node

  DO I=1,NCV_I
     IA=NIEC(I,1)
     IB=NIEC(I,2)     
     XI=0.5_SP*(XIJE(I,1)+XIJE(I,2))
     YI=0.5_SP*(YIJE(I,1)+YIJE(I,2))
#      if defined (SPHERICAL)
     X1_DP=XIJE(I,1)
     Y1_DP=YIJE(I,1)
     X2_DP=XIJE(I,2)
     Y2_DP=YIJE(I,2)
     CALL ARCC(X2_DP,Y2_DP,X1_DP,Y1_DP,XII,YII)
     XI=XII		
     XTMP  = XI*TPI-VX(IA)*TPI
     XTMP1 = XI-VX(IA)
     IF(XTMP1 >  180.0_SP)THEN
        XTMP = -360.0_SP*TPI+XTMP
     ELSE IF(XTMP1 < -180.0_SP)THEN
        XTMP =  360.0_SP*TPI+XTMP
     END IF

     DLTXNCVE(I,1)=XTMP*COS(DEG2RAD*VY(IA))    
     DLTYNCVE(I,1)=(YI-VY(IA))*TPI

     XTMP  = XI*TPI-VX(IB)*TPI
     XTMP1 = XI-VX(IB)
     IF(XTMP1 >  180.0_SP)THEN
        XTMP = -360.0_SP*TPI+XTMP
     ELSE IF(XTMP1 < -180.0_SP)THEN
        XTMP =  360.0_SP*TPI+XTMP
     END IF

     DLTXNCVE(I,2)=XTMP*COS(DEG2RAD*VY(IB)) 
     DLTYNCVE(I,2)=(YI-VY(IB))*TPI
#      else
     DLTXNCVE(I,1)=XI-VX(IA)
     DLTYNCVE(I,1)=YI-VY(IA)
     DLTXNCVE(I,2)=XI-VX(IB)
     DLTYNCVE(I,2)=YI-VY(IB)
#      endif

  END DO


  ! Set the distance between Nodes

  IF(MAXVAL(NTSN) > 13) CALL FATAL_ERROR &
       & ("THERE ARE MORE THAN 12 NODES AROUND ONE NODE:",&
       "PLEASE INCREASE THE SIZE OF DLTXPI AND DLTYPI IN MOD_MAIN",&
       "BUT REALLY, WHAT IS WRONG WITH YOUR MESH?")

  DO I=1,M
     DO J=1,NTSN(I)-1
        I1=NBSN(I,J)
        I2=NBSN(I,J+1)
# if defined (SPHERICAL)
        XTMP  = VX(I2)*TPI-VX(I1)*TPI
        XTMP1 = VX(I2)-VX(I1)
        IF(XTMP1 >  180.0_SP)THEN
           XTMP = -360.0_SP*TPI+XTMP
        ELSE IF(XTMP1 < -180.0_SP)THEN
           XTMP =  360.0_SP*TPI+XTMP
        END IF
        DLTXTRIE(i,j) =XTMP*COS(DEG2RAD*VY(I))
        DLTYTRIE(i,j) =(VY(I1)-VY(I2))*TPI

# else

        DLTYTRIE(i,j) = VY(I1)-VY(I2)
        DLTXTRIE(i,j) = VX(I2)-VX(I1)
# endif

     END DO
  END DO

  ! Set the distance between Nodes for the North Pole region
# if defined (SPHERICAL)
  DO IP=1,MP
     I = MP_LST(IP)
     DO J=1,NTSN(I)-1
        I1=NBSN(I,J)
        I2=NBSN(I,J+1)

        VX1_TMP = REARTH * COS(VY(I1)*PI/180.0_SP) * COS(VX(I1)*PI/180.0_SP) &
             * 2._SP /(1._SP+sin(VY(I1)*PI/180.0_SP))
        VY1_TMP = REARTH * COS(VY(I1)*PI/180.0_SP) * SIN(VX(I1)*PI/180.0_SP) &
             * 2._SP /(1._SP+sin(VY(I1)*PI/180.0_SP))

        VX2_TMP = REARTH * COS(VY(I2)*PI/180.0_SP) * COS(VX(I2)*PI/180.0_SP) &
             * 2._SP /(1._SP+sin(VY(I2)*PI/180.0_SP))
        VY2_TMP = REARTH * COS(VY(I2)*PI/180.0_SP) * SIN(VX(I2)*PI/180.0_SP) &
             * 2._SP /(1._SP+sin(VY(I2)*PI/180.0_SP))

        DLTXTRIE(i,j) = (VX2_TMP-VX1_TMP)/(2._SP /(1._SP+sin(VY(I)*PI/180.0_SP)))
        DLTYTRIE(i,j) = (VY1_TMP-VY2_TMP)/(2._SP /(1._SP+sin(VY(I)*PI/180.0_SP)))

        IF(I /= NODE_NORTHPOLE)THEN
           TXPI_A = DLTXTRIE(i,j)
           TYPI_A = DLTYTRIE(i,j)

           TXPI_B = TYPI_A*COS(VX(I)*PI/180._SP)&
                & -TXPI_A*SIN(VX(I)*PI/180._SP)

           TYPI_B = TXPI_A*COS(VX(I)*PI/180._SP)&
                & +TYPI_A*SIN(VX(I)*PI/180._SP)

           DLTXTRIE(i,j) = TXPI_B
           DLTYTRIE(i,j) =-TYPI_B
        END IF

     END DO
  END DO
# endif

  ! Set the distance between triangle edge centers
  IF(MAXVAL(NTVE) > 13) CALL FATAL_ERROR &
       & ("THERE ARE MORE THAN 12 CELLS AROUND ONE NODE:",&
       "PLEASE INCREASE THE SIZE OF DLVISCXPI AND DLVISCYPI IN MOD_MAIN",&
       "BUT REALLY, WHAT IS WRONG WITH YOUR MESH?")


  DO I=1,M
     DO J=1,NTVE(I)
        I1=NBVE(I,J)
        JTMP=NBVT(I,J)
        J1=JTMP+1-(JTMP+1)/4*3
        J2=JTMP+2-(JTMP+2)/4*3
        X11=0.5_SP*(VX(I)+VX(NV(I1,J1)))
        Y11=0.5_SP*(VY(I)+VY(NV(I1,J1)))
!        X22=XC(I1)
!        Y22=YC(I1)
        X33=0.5_SP*(VX(I)+VX(NV(I1,J2)))
        Y33=0.5_SP*(VY(I)+VY(NV(I1,J2)))

#      if defined (SPHERICAL)
        X1_DP=VX(I)
        Y1_DP=VY(I)
        X2_DP=VX(NV(I1,J1))
        Y2_DP=VY(NV(I1,J1))
        CALL ARCC(X2_DP,Y2_DP,X1_DP,Y1_DP,X11_TMP,Y11_TMP)
        X11=X11_TMP
        Y11=Y11_TMP
        X2_DP=VX(NV(I1,J2))
        Y2_DP=VY(NV(I1,J2))
        CALL ARCC(X2_DP,Y2_DP,X1_DP,Y1_DP,X33_TMP,Y33_TMP)
        X33=X33_TMP
        Y33=Y33_TMP

        XTMP  = X33*TPI-X11*TPI
        XTMP1 = X33-X11
        IF(XTMP1 >  180.0_SP)THEN
           XTMP = -360.0_SP*TPI+XTMP
        ELSE IF(XTMP1 < -180.0_SP)THEN
           XTMP =  360.0_SP*TPI+XTMP
        END IF

        DLTYECEC(I,J)=(Y11-Y33)*TPI
        DLTXECEC(I,J)=XTMP*COS(DEG2RAD*VY(I))
#      else

        DLTYECEC(I,J)=(Y11-Y33)
        DLTXECEC(I,J)=(X33-X11)

#      endif


! Set the distance between the node and the edge Center
! NOTE: THE SIGN MATTERS!
#        if defined (SPHERICAL)
         XTMP  = X11*TPI-VX(I)*TPI
         XTMP1 = X11-VX(I)
         IF(XTMP1 >  180.0_SP)THEN
	   XTMP = -360.0_SP*TPI+XTMP
         ELSE IF(XTMP1 < -180.0_SP)THEN
	   XTMP =  360.0_SP*TPI+XTMP
         END IF  
         DLTXNEC(I,J)= XTMP*COS(DEG2RAD*VY(I))
         DLTYNEC(I,J)= (VY(I)-Y11)*TPI

#        else
         DLTYNEC(I,J)=(VY(I)-Y11)
         DLTXNEC(I,J)=(X11-VX(I))
#        endif


     END DO
  END DO


   END SUBROUTINE edge_len
