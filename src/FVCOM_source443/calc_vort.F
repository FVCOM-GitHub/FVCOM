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
!  CALCULATE VORTICITY FIELD AT NODES                                          |
!    OPTIONAL INPUT: VTYPE                                                     | 
!              1    CURL OF VERTICALLY-AVERAGED VELOCITY (DEFAULT)             |
!              2    CURL OF FLUX (For Streamfunctions)                         |
!    RETURNS:  VORT(1:M)                                                       |
!                                                                              |
!    CHECK, INITIALIZE WITH ANALYTICAL VELOCITY FIELD (UA = -.5YC, VA = .5XC)  !
!==============================================================================|
   SUBROUTINE CALC_VORT(VTYPE) 
!==============================================================================|
   USE ALL_VARS
   IMPLICIT NONE
   INTEGER, OPTIONAL, INTENT(IN) :: VTYPE
   REAL(SP) :: XFLUX(0:M)
   REAL(SP) :: UIJ,VIJ,EXFLUX,AVE
   INTEGER  :: I,J,IA,IB,I1,CNT,JNODE,MY_VTYPE
   INTEGER, ALLOCATABLE  :: LIST(:)
!==============================================================================|

  !---------------------------------------------
  !Process input arguments
  !---------------------------------------------
  MY_VTYPE = 1
  IF(present(VTYPE))Then
    MY_VTYPE = VTYPE
  ENDIF

  !---------------------------------------------
  !Initialize Fields
  !---------------------------------------------
  ALLOCATE(LIST(M)) ; LIST = 0
  VORT = 0.0

  !===================================================================
  ! Calculate vorticity at all nodes using c.v. half edges
  !   Note minus sign (nodes are order clockwise [against convention])
  !===================================================================

  SELECT CASE(MY_VTYPE)

  !----------------------------------------------
  CASE(1) !Curl of Vert-Avged Velocity
  !----------------------------------------------
  DO I=1,NCV
    I1  = NTRG(I)
    IA  = NIEC(I,1)
    IB  = NIEC(I,2)
    UIJ = UA(I1)
    VIJ = VA(I1) 
    EXFLUX = -(UIJ*DLTXE(I) + VIJ*DLTYE(I))  
    VORT(IA) = VORT(IA)-EXFLUX
    VORT(IB) = VORT(IB)+EXFLUX
  END DO
  !----------------------------------------------
  CASE(2) !Curl of Flux
  !----------------------------------------------
  DO I=1,NCV
    I1  = NTRG(I)
    IA  = NIEC(I,1)
    IB  = NIEC(I,2)
    UIJ = UA(I1)
    VIJ = VA(I1) 
    EXFLUX = -D1(I1)*(UIJ*DLTXE(I) + VIJ*DLTYE(I))  
    VORT(IA) = VORT(IA)-EXFLUX
    VORT(IB) = VORT(IB)+EXFLUX
  END DO
  !----------------------------------------------
  CASE DEFAULT
    !ERROR AND HALT 
  END SELECT

  VORT = VORT/ART1

  !===================================================================
  !Correction at Boundaries (May no longer be necessary for 2.5+)
  !===================================================================

   do i=1,m
     if(isonb(i) > 0)then
     ave = 0.
     cnt = 0
     do j=1,ntsn(i)
       jnode = nbsn(i,j)
       if(jnode /= i .and. isonb(jnode) == 0)then
         ave = ave + vort(jnode)
         cnt = cnt + 1
       end if
     end do
     vort(i) = ave/float(cnt)
     if(cnt == 0) list(i) = 1
     end if
   end do

   do i=1,m
     if(list(i) > 0)then
     ave = 0.
     cnt = 0
     do j=1,ntsn(i)
       jnode = nbsn(i,j)
       if(jnode /= i .and. list(jnode) == 0)then
         ave = ave + vort(jnode)
         cnt = cnt + 1
       end if
     end do
     vort(i) = ave/float(cnt)
     end if
   end do

       
   RETURN
   END SUBROUTINE CALC_VORT  
!==============================================================================|

