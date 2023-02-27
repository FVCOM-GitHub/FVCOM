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
!                                                                	       |
!  this subroutine linearly interpolates and extrapolates an     	       | 
!  array b.                                                      	       |
!                                                                	       |
!  x(m1) must be descending                                                    |
!  a(x) given function                                           	       |
!  b(y) found by linear interpolation and extrapolation          	       |	
!  y(n1) the desired depths                                      	       |	
!  m1   the number of points in x and a                          	       |	
!  n1   the number of points in y and b                         	       |	
!                                                                              |	
!==============================================================================|

MODULE SINTER
  IMPLICIT NONE

  CONTAINS

   SUBROUTINE SINTER_EXTRP_UP(X,A,Y,B,M1,N1)              

!==============================================================================|
   USE MOD_PREC
   IMPLICIT NONE
   INTEGER, INTENT(IN)   :: M1,N1
   REAL(SP), INTENT(IN)  :: X(M1),A(M1),Y(N1)
   REAL(SP), INTENT(OUT) :: B(N1) 
   INTEGER   I,J,NM
!==============================================================================|

!
!  EXTRAPOLATION
!
   DO I=1,N1
       IF (Y(I) > X(1 )) B(I) = A(1) + ((A(1)-A(2))/(X(1)-X(2))) * (Y(I)-X(1))
       IF (Y(I) < X(M1)) B(I) = A(M1)
   END DO

!
!  INTERPOLATION
!
   NM = M1 - 1
   DO I=1,N1
     DO J=1,NM
        IF (Y(I) <=  X(J) .AND. Y(I) >= X(J+1)) &
              B(I) = A(J) - (A(J)- A(J+1)) * (X(J)-Y(I)) / (X(J)-X(J+1))
     END DO
   END DO

   RETURN
 END SUBROUTINE SINTER_EXTRP_UP
!==============================================================================|

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!==============================================================================|
   SUBROUTINE SINTER_EXTRP_DOWN(X,A,Y,B,M1,N1)              
!  for baroclinic interpolation                                                |
!==============================================================================|


!==============================================================================|
   USE MOD_PREC
   IMPLICIT NONE
   INTEGER, INTENT(IN)  :: M1,N1
   REAL(SP),  INTENT(IN)  :: X(M1),A(M1),Y(N1)
   REAL(SP),  INTENT(OUT) :: B(N1) 
   INTEGER :: I,J,NM
!==============================================================================|

!
!  EXTRAPOLATION 
!
   DO I=1,N1
     IF(Y(I) > X(1 )) B(I) = A(1)
     IF(Y(I) < X(M1)) B(I)=A(M1)+(A(M1-1)-A(M1))*(Y(I)-X(M1))/(X(M1-1)-X(M1))
   END DO

!
!  INTERPOLATION
!
   NM = M1 - 1
   DO I=1,N1
     DO J=1,NM
       IF (Y(I)<=X(J).AND.Y(I)>=X(J+1)) &
          B(I) = A(J) - (A(J)- A(J+1)) *(X(J)-Y(I)) / (X(J)-X(J+1))
     END DO
   END DO

   RETURN
 END SUBROUTINE SINTER_EXTRP_DOWN
!==============================================================================|


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!==============================================================================|
   SUBROUTINE SINTER_EXTRP_NONE(X,A,Y,B,M1,N1)              
!  for t&s obc interpolation                                                   |
!==============================================================================|


   USE MOD_PREC
   IMPLICIT NONE
   INTEGER, INTENT(IN)  :: M1,N1
   REAL(SP),  INTENT(IN)  :: X(M1),A(M1),Y(N1)
   REAL(SP),  INTENT(OUT) :: B(N1) 
   INTEGER :: I,J,NM
!==============================================================================|


!
!   EXTRAPOLATION
!
   DO I=1,N1
     IF (Y(I) > X(1 )) B(I) = A(1)
     IF (Y(I) < X(M1)) B(I) = A(M1)
   END DO

!
!  INTERPOLATION
!
   NM = M1 - 1
   DO I=1,N1
     DO J=1,NM
       IF (Y(I) <= X(J).AND.Y(I) >= X(J+1))  &
       B(I) = A(J) - (A(J)- A(J+1)) * (X(J)-Y(I)) / (X(J)-X(J+1))
     END DO
   END DO

   RETURN
 END SUBROUTINE SINTER_EXTRP_NONE
!==============================================================================|

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!==============================================================================|
   SUBROUTINE SINTER_EXTRP_BOTH(X,A,Y,B,M1,N1)              
!==============================================================================|


!==============================================================================|
   USE MOD_PREC
   IMPLICIT NONE
   INTEGER, INTENT(IN)  :: M1,N1
   REAL(SP),  INTENT(IN)  :: X(M1),A(M1),Y(N1)
   REAL(SP),  INTENT(OUT) :: B(N1) 
   INTEGER :: I,J,NM
!==============================================================================|

!
!  EXTRAPOLATION 
!
   DO I=1,N1
     IF(Y(I) > X(1 )) B(I) = A(1) + ((A(1)-A(2))/(X(1)-X(2))) * (Y(I)-X(1))
     IF(Y(I) < X(M1)) B(I)=A(M1)+(A(M1-1)-A(M1))*(Y(I)-X(M1))/(X(M1-1)-X(M1))
   END DO

!
!  INTERPOLATION
!
   NM = M1 - 1
   DO I=1,N1
     DO J=1,NM
       IF (Y(I)<=X(J).AND.Y(I)>=X(J+1)) &
          B(I) = A(J) - (A(J)- A(J+1)) *(X(J)-Y(I)) / (X(J)-X(J+1))
     END DO
   END DO

   RETURN
 END SUBROUTINE SINTER_EXTRP_BOTH
!==============================================================================|

END MODULE SINTER
