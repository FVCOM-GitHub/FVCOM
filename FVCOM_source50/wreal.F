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
!  Compute Cartesian Vertical Velocity                                         |
!==============================================================================|
   SUBROUTINE WREAL               
!==============================================================================|
   USE ALL_VARS
   USE MOD_UTILS
#  if defined (WET_DRY)
   USE MOD_WD
#  endif
   USE MOD_NESTING
   IMPLICIT NONE
   REAL(SP) :: DDDX,DDDY,DEDX,DEDY,ETF1AA,WW1,WW2 
#   if defined (ICE_EMBEDDING)
   REAL(SP) :: ICEF1AA       !!yding 
   REAL(SP) :: DHIDX,DHIDY   !!yding
#   endif
   INTEGER  :: I,K,J1,J2,J3
   INTEGER  :: J, II
   REAL(SP) :: U_TMP, V_TMP
!==============================================================================|

   IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "Start: wreal"

!------------------------------------------------------------------------------!
!  SAVE OMEGA VELOCITY FROM PREVIOUS TIMESTEP (USED FOR LAGRANGIAN TRACKING)   !
!------------------------------------------------------------------------------!

   WTTS = WTS

!------------------------------------------------------------------------------!
!  CALCULATE A NEW OMEGA VELOCITY                                              !
!------------------------------------------------------------------------------!

!!===========yding==================
#   if defined (ICE_EMBEDDING)
       DO I=1,N
#  if defined (WET_DRY)
    IF(ISWETC(I) == 1)THEN
#  endif
     J1=NV(I,1)
     J2=NV(I,2)
     J3=NV(I,3)
     DDDX=AWX(I,1) * D(J1)+AWX(I,2) * D(J2)+AWX(I,3)*D(J3)
     DDDY=AWY(I,1) * D(J1)+AWY(I,2) * D(J2)+AWY(I,3)*D(J3)
     DEDX=AWX(I,1)*ELF(J1)+AWX(I,2)*ELF(J2)+AWX(I,3)*ELF(J3)
     DEDY=AWY(I,1)*ELF(J1)+AWY(I,2)*ELF(J2)+AWY(I,3)*ELF(J3)

     DHIDX=AWX(I,1)*QTHICERK(J1)+AWX(I,2)*QTHICERK(J2)+AWX(I,3)*QTHICERK(J3)
     DHIDY=AWY(I,1)*QTHICERK(J1)+AWY(I,2)*QTHICERK(J2)+AWY(I,3)*QTHICERK(J3)

     ETF1AA=ONE_THIRD*(EL(NV(I,1))+EL(NV(I,2))+EL(NV(I,3)))

     ICEF1AA=ONE_THIRD*(QTHICE_incr(NV(I,1))+QTHICE_incr(NV(I,2))+QTHICE_incr(NV(I,3)))  !!--yding
    

      DO K=1,KBM1
      WW1=0.5_SP*(W(I,K)+W(I,K+1))+U(I,K)*(ZZ1(I,K)*DDDX+DEDX-DHIDX)+ &
                                V(I,K)*(ZZ1(I,K)*DDDY+DEDY-DHIDY)
      WW2=(ZZ1(I,K)+1.)*((ETF1AA-ET1(I))/DTI - ICEF1AA/DTI)
      WW(I,K)=WW1+WW2
     END DO
#  if defined (WET_DRY)
    ELSE
     DO K=1,KBM1
      WW(I,K)=0.0_SP
     END DO
    END IF
#  endif
   END DO


#   else


   DO I=1,N
#  if defined (WET_DRY)
    IF(ISWETC(I) == 1)THEN
#  endif
     J1=NV(I,1)
     J2=NV(I,2)
     J3=NV(I,3)
     DDDX=AWX(I,1) * D(J1)+AWX(I,2) * D(J2)+AWX(I,3)*D(J3)
     DDDY=AWY(I,1) * D(J1)+AWY(I,2) * D(J2)+AWY(I,3)*D(J3)
     DEDX=AWX(I,1)*ELF(J1)+AWX(I,2)*ELF(J2)+AWX(I,3)*ELF(J3)
     DEDY=AWY(I,1)*ELF(J1)+AWY(I,2)*ELF(J2)+AWY(I,3)*ELF(J3)
     ETF1AA=ONE_THIRD*(EL(NV(I,1))+EL(NV(I,2))+EL(NV(I,3)))
     DO K=1,KBM1
      WW1=0.5_SP*(W(I,K)+W(I,K+1))+U(I,K)*(ZZ1(I,K)*DDDX+DEDX)+ &
                                V(I,K)*(ZZ1(I,K)*DDDY+DEDY)
      WW2=(ZZ1(I,K)+1.)*(ETF1AA-ET1(I))/DTI
      WW(I,K)=WW1+WW2
     END DO
#  if defined (WET_DRY)
    ELSE
     DO K=1,KBM1
      WW(I,K)=0.0_SP
     END DO
    END IF
#  endif
   END DO

#  endif
!!ice_embedding  yding

# if defined(MULTIPROCESSOR)
  IF(PAR)CALL AEXCHANGE(EC,MYID,NPROCS,WW)
# endif

   IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "End: wreal"

   END SUBROUTINE WREAL
!==============================================================================|
