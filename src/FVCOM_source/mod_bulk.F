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

MODULE mod_bulk
  USE MOD_PREC
  implicit none

CONTAINS



  SUBROUTINE ASIMPLE_DRAG(spdx,spdy,strx,stry)
    IMPLICIT NONE
    REAL(SP),ALLOCATABLE, TARGET, INTENT(IN)  :: SPDX(:),SPDY(:)
    REAL(SP),ALLOCATABLE, TARGET, INTENT(INOUT) :: STRX(:),STRY(:)

    REAL(SP), POINTER :: SPDXP(:),SPDYP(:)
    REAL(SP), POINTER :: STRXP(:),STRYP(:)
    
    SPDXP => SPDX
    SPDYP => SPDY

    STRXP => STRX
    STRYP => STRY
    
    CALL PSIMPLE_DRAG(SPDXP,SPDYP,STRXP,STRYP)


  END SUBROUTINE ASIMPLE_DRAG



  SUBROUTINE PSIMPLE_DRAG(spdx,spdy,strx,stry)
    IMPLICIT NONE
    REAL(SP), POINTER,INTENT(IN)  :: SPDX(:),SPDY(:)
    REAL(SP), POINTER, INTENT(INOUT) :: STRX(:),STRY(:)
    INTEGER :: I, N
    REAL(SP) :: CD, WDS, TX, TY


    IF(.not.Associated(SPDX)) WRITE(6,*) "SIMPLE DRAG: SPDX is not associated"
    IF(.not.Associated(SPDY)) WRITE(6,*) "SIMPLE DRAG: SPDY is not associated"
    IF(.not.Associated(STRX)) WRITE(6,*) "SIMPLE DRAG: STRX is not associated"
    IF(.not.Associated(STRY)) WRITE(6,*) "SIMPLE DRAG: STRY is not associated"

    N = UBOUND(SPDX,1)


    IF(N /= UBOUND(SPDY,1)) WRITE(6,*) "SIMPLE DRAG: MIS-MATCHED DIMENSIONS"
    IF(N /= UBOUND(STRY,1)) WRITE(6,*) "SIMPLE DRAG: MIS-MATCHED DIMENSIONS"
    IF(N /= UBOUND(STRX,1)) WRITE(6,*) "SIMPLE DRAG: MIS-MATCHED DIMENSIONS"

    DO I=1,N
       TX = SPDX(I)
       TY = SPDY(I)
       WDS=SQRT(TX*TX+TY*TY)
       CD=1.2E-3
       IF (WDS >= 11.0_SP) CD=(0.49_SP+0.065_SP*WDS)*1.E-3_SP
       IF (WDS >= 25.0_SP) CD=(0.49_SP+0.065_SP*25.0_SP)*1.E-3_SP

       STRX(I) = 1.2_SP*CD*TX*WDS
       STRY(I) = 1.2_SP*CD*TY*WDS

    END DO
    

  END SUBROUTINE PSIMPLE_DRAG
  


END MODULE mod_bulk
