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

# if !defined (SEMI_IMPLICIT)
!==============================================================================|
!  CALCULATE n+1 time step perturbation water surface elevation                |
!==============================================================================|
   SUBROUTINE EXTELPF_EDGE(KKT)       
!==============================================================================|
   USE ALL_VARS
   USE BCS

   USE MOD_OBCS2
   USE MOD_OBCS3

   IMPLICIT NONE
   INTEGER, INTENT(IN) :: KKT
   REAL(SP) :: XFLUX(0:MT)
   REAL(SP) :: DIJ,UIJ,VIJ,DTK,EXFLUX,FXFLUX,FFLUX(1:IOBCN),TP,CC,CP
   INTEGER  :: I,J,I1,IA,IB,JJ,J1,J2,J3,II,JN

!---------ACCUMULATE FLUX BY LOOPING OVER CONTROL VOLUME HALF EDGES------------!

   XFLUX = 0.0_SP
   DO II=1,NOBCV
     I = NOBEDGE_LST(II)
     I1  = NTRG(I)
     IA  = NIEC(I,1)
     IB  = NIEC(I,2)
     J1  = I_OBC_NODE(nv(I1,1))
     J2  = I_OBC_NODE(nv(I1,2))
     J3  = I_OBC_NODE(nv(I1,3))
     DIJ = H1(I1) + one_third*(ELT(J1)+ELT(J2)+ELT(J3))
     UIJ = UAT(I_OBC_CELL(I1))
     VIJ = VAT(I_OBC_CELL(I1))
     EXFLUX = DIJ*(-UIJ*DLTYE(I) + VIJ*DLTXE(I))  

#    if defined (PLBC)
     EXFLUX = DIJ*(-UIJ*DLTYE(I) + 0.0_SP*DLTXE(I))  
#    endif

     XFLUX(IA) = XFLUX(IA)-EXFLUX
     XFLUX(IB) = XFLUX(IB)+EXFLUX
   END DO

   DO I = 1, IOBCN
      J = I_OBC_N(I)
      J1= I_OBC_NODE(J)
      FXFLUX = -(ELTF(J1)-ELRKT(J1))*ART1(J)/(ALPHA_RK(KKT)*DTE) - XFLUX(J)
      FFLUX(I)= FXFLUX / (H(J)+ELT(J1)) * ELP(J1)
   END DO
       
   XFLUX = 0.0_SP
   DO II=1,NOBCV
     I = NOBEDGE_LST(II)
     I1  = NTRG(I)
     IA  = NIEC(I,1)
     IB  = NIEC(I,2)
     J1  = I_OBC_NODE(nv(I1,1))
     J2  = I_OBC_NODE(nv(I1,2))
     J3  = I_OBC_NODE(nv(I1,3))
     DIJ = one_third * (ELP(J1)+ELP(J2)+ELP(J3))
     UIJ = UAT(I_OBC_CELL(I1))
     VIJ = VAT(I_OBC_CELL(I1))
     EXFLUX = DIJ*(-UIJ*DLTYE(I) + VIJ*DLTXE(I))

#    if defined (PLBC)
     EXFLUX = DIJ*(-UIJ*DLTYE(I) + 0.0_SP*DLTXE(I))
#    endif
  
     XFLUX(IA) = XFLUX(IA)-EXFLUX
     XFLUX(IB) = XFLUX(IB)+EXFLUX
   END DO

   DO I = 1, IOBCN
      J = I_OBC_N(I)
      FFLUX(I) = FFLUX(I) + XFLUX(J)
   END DO

   XFLUX = 0.0_SP
   DO II=1,NOBCV
     I = NOBEDGE_LST(II)
     I1  = NTRG(I)
     IA  = NIEC(I,1)
     IB  = NIEC(I,2)
     J1  = I_OBC_NODE(nv(I1,1))
     J2  = I_OBC_NODE(nv(I1,2))
     J3  = I_OBC_NODE(nv(I1,3))
     DIJ = H1(I1) + one_third*(ELT(J1)+ELT(J2)+ELT(J3)+ &
           ELP(J1)+ ELP(J2)+ELP(J3))
     UIJ = UAP(I_OBC_CELL(I1))
     VIJ = VAP(I_OBC_CELL(I1))
     EXFLUX = DIJ*(-UIJ*DLTYE(I) + VIJ*DLTXE(I))  

#    if defined (PLBC)
     EXFLUX = DIJ*(-UIJ*DLTYE(I) + 0.0_SP*DLTXE(I))
#   endif

     XFLUX(IA) = XFLUX(IA)-EXFLUX
     XFLUX(IB) = XFLUX(IB)+EXFLUX
   END DO

   DO I = 1, IOBCN
      J = I_OBC_N(I)
      FFLUX(I) = FFLUX(I) + XFLUX(J) + FLUXOBN2(I)
   END DO
   FLUXOBN2 = 0.0_SP

   IF (IBCN(1) > 0) THEN
   DO I  = 1, IBCN(1)
      JN = OBC_LST(1,I)
      J  = I_OBC_N(JN)
      J1 = I_OBC_NODE(J)
      ELPF(J1) = ELRKP(J1) - ALPHA_RK(KKT)*DTE*FFLUX(JN)/ART1(J)
   END DO
   END IF

   TP = 10800.0_SP       ! 3 hours
   IF (IBCN(4) > 0) THEN
   DO I  = 1, IBCN(4)
      JN = OBC_LST(4,I)
      J  = I_OBC_N(JN)
      J1 = I_OBC_NODE(J)

      I1 = I_OBC_NODE(NEXT_OBC(JN))
      CC = SQRT(GRAV_N(J)*D(J))*ALPHA_RK(KKT)*DTE/DLTN_OBC(JN)
      CP = CC + 1.0_SP

      ELPF(J1) = (CC*ELPF(I1) + ELRKP(J1)*(1.0_SP-ALPHA_RK(KKT)*DTE/TP))/CP
   END DO
   END IF
 
   IF((IBCN(2)>0).or.(IBCN(3)>0).or.(IBCN(5)>0)) print *,"error in EXTELPF_EDGE.F"

   RETURN
   END SUBROUTINE EXTELPF_EDGE
!==============================================================================|
# endif
