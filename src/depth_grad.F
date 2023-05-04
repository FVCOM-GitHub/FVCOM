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
!     CALCULATE THE GRADIENT OF THE WATER DEPTH WITH X ANDY                    |
!==============================================================================|

   SUBROUTINE DEPTH_GRADIENT      

!==============================================================================|
   USE ALL_VARS
   USE MOD_PAR

   IMPLICIT NONE
   REAL(SP) :: PHPX,PHPY,F1
   INTEGER  :: I,J,I1,I2
# if defined (SPHERICAL)
   REAL(DP) :: XTMP,XTMP1
# endif
!==============================================================================|

!----------CALCULATE DERIVATIVES OF DEPTH WITH X AND Y AT NODES----------------!

   DO I=1,M
     PHPX = 0.0_SP ; PHPY = 0.0_SP
     DO J=1,NTSN(I)-1
       I1 = NBSN(I,J)
       I2 = NBSN(I,J+1)
       F1 = 0.50_SP*(H(I1)+H(I2)) 
#      if defined (SPHERICAL)
       PHPX=PHPX+F1*(VY(I1)-VY(I2))*TPI
!<for global case when VX change from 360 to 0 or from 0 to 360 (degree)
       XTMP  = VX(I2)*TPI-VX(I1)*TPI  
       XTMP1 = VX(I2)-VX(I1)  
       IF(XTMP1 >  180.0_SP)THEN
         XTMP = -360.0_SP*TPI+XTMP     
       ELSE IF(XTMP1 < -180.0_SP)THEN
         XTMP =  360.0_SP*TPI+XTMP     
       END IF
       PHPY=PHPY+F1*XTMP*COS(DEG2RAD*VY(I))
!       PHPY=PHPY+F1*(VX(I2)-VX(I1))*TPI*COS(DEG2RAD*TY)
!>
#      else
       PHPX = PHPX + F1*(VY(I1)-VY(I2))
       PHPY = PHPY + F1*(VX(I2)-VX(I1))
#      endif
     END DO
     PHPX = PHPX/ART2(I)
     PHPY = PHPY/ART2(I)

     IF(PHPX==0.0_SP .AND.  PHPY ==0.0_SP)THEN
       SITA_GD(I) = 0.0_SP ; PHPN(I) = 0.0_SP
     ELSE
       SITA_GD(I) = ATAN2(PHPY,PHPX)
       PHPN(I)    = SQRT(PHPY*PHPY+PHPX*PHPX)
     END IF

   END DO

#  if defined (NH) && (MULTIPROCESSOR)
   IF(PAR) CALL AEXCHANGE(NC,MYID,NPROCS,SITA_GD,PHPN)
#  endif

   RETURN
   END SUBROUTINE DEPTH_GRADIENT
!==============================================================================|
