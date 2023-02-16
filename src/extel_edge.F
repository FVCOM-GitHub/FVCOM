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
!  CALCULATE FLUXES OF FREE SURFACE ELEVATION (CONTINUITY) EQUATION            |
!==============================================================================|
   SUBROUTINE EXTEL_EDGE(K)       
# if !defined (SEMI_IMPLICIT)
!==============================================================================|
   USE ALL_VARS
   USE BCS
   USE MOD_OBCS

#  if defined (SPHERICAL)
   USE MOD_NORTHPOLE
#  endif        

#  if defined (BALANCE_2D)
   USE MOD_BALANCE_2D
#  endif

!  ggao 0903/2007
#  if defined (ICE)
   use mod_ice,only : isicen  !,fresh,RHOW
#  endif

#  if defined (THIN_DAM)
   use mod_dam
#  endif


   IMPLICIT NONE
   INTEGER, INTENT(IN) :: K
   REAL(SP) :: XFLUX(0:MT)
   REAL(SP) :: DIJ,UIJ,VIJ,DTK,UN,EXFLUX
   INTEGER  :: I,J,I1,IA,IB,JJ,J1,J2

#  if defined (BALANCE_2D)
   REAL(SP), DIMENSION(0:MT) :: XFLUXU,XFLUXV
   REAL(SP), DIMENSION(0:NT) :: XFLUXU1,XFLUXV1
   REAL(SP) :: EXFLUX_U,EXFLUX_V
#  endif

#  if defined (THIN_DAM)
   INTEGER :: JN
#  endif
!==============================================================================|
   IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "Start: extel_edge.F",K
!----------INITIALIZE FLUX ARRAY ----------------------------------------------!

   XFLUX = 0.0_SP

#  if defined (BALANCE_2D)
   XFLUXU= 0.0_SP
   XFLUXV= 0.0_SP
   XFLUXU1= 0.0_SP
   XFLUXV1= 0.0_SP
#  endif
!---------ACCUMULATE FLUX BY LOOPING OVER CONTROL VOLUME HALF EDGES------------!

   DO I=1,NCV
     I1  = NTRG(I)
     IA  = NIEC(I,1)
     IB  = NIEC(I,2)
     DIJ = D1(I1)

     UIJ = UA(I1)
     VIJ = VA(I1)
     EXFLUX = DIJ*(-UIJ*DLTYE(I) + VIJ*DLTXE(I))  

#    if defined (PLBC)
     EXFLUX = DIJ*(-UIJ*DLTYE(I) + 0.0_SP*DLTXE(I)) 
#    endif

     XFLUX(IA) = XFLUX(IA)-EXFLUX
     XFLUX(IB) = XFLUX(IB)+EXFLUX

#  if defined (BALANCE_2D)
     EXFLUX_U = -DIJ*UIJ*DLTYE(I)
     EXFLUX_V =  DIJ*VIJ*DLTXE(I)
     XFLUXU(IA) = XFLUXU(IA)-EXFLUX_U
     XFLUXU(IB) = XFLUXU(IB)+EXFLUX_U
     XFLUXV(IA) = XFLUXV(IA)-EXFLUX_V
     XFLUXV(IB) = XFLUXV(IB)+EXFLUX_V
#  endif

   END DO

!   write(ipt,*) "after control volume flux",maxval(xflux),minval(xflux)
!   write(ipt,*) "after control volume flux",maxloc(xflux),minloc(xflux)

#  if defined (SPHERICAL)
   CALL EXTEL_EDGE_XY(K,XFLUX)
#  endif

!--ADD EVAPORATION AND PRECIPITATION TERMS-------------------------------------!
   IF (PRECIPITATION_ON) THEN
#  if defined (ICE)
      WHERE(ISICEN==1)
         QPREC2=QPREC
         QEVAP2=QEVAP
      END WHERE

!!$       DO I=1,M
!!$        IF(ISICEN(i)==1) THEN
!!$          QPREC2(I)= QPREC3(I)
!!$          QEVAP2(I)= QEVAP3(I)
!!$        END IF
!!$       END DO
#  endif

!qxu   XFLUX = XFLUX+(QEVAP2-QPREC2)*ROFVROS*ART1 
!qxu---the evap is negative for evaporating in ocean
      XFLUX = XFLUX-(QEVAP2+QPREC2)*ROFVROS*ART1 
   END IF
    
!--ADD GROUND WATER TERM-------------------------------------------------------!

   IF(GROUNDWATER_ON) THEN
      XFLUX = XFLUX - BFWDIS2
   END IF

!--SAVE ACCUMULATED FLUX ON OPEN BOUNDARY NODES AND ZERO OUT OPEN BOUNDARY FLUX!

   IF(IOBCN > 0) THEN  
     DO I=1,IOBCN
       XFLUX_OBCN(I)=XFLUX(I_OBC_N(I))
       XFLUX(I_OBC_N(I)) = 0.0_SP
     END DO
   END IF


!---------ADJUST FLUX FOR FRESH WATER DISCHARGE--------------------------------!

   IF(NUMQBC >= 1) THEN   
     IF(RIVER_INFLOW_LOCATION == 'node') THEN
       DO J=1,NUMQBC
         JJ=INODEQ(J)
         XFLUX(JJ)=XFLUX(JJ)-QDIS(J)
       END DO
     ELSE IF(RIVER_INFLOW_LOCATION == 'edge') THEN
       DO J=1,NUMQBC
         J1=N_ICELLQ(J,1)
         J2=N_ICELLQ(J,2)
         XFLUX(J1)=XFLUX(J1)-QDIS(J)*RDISQ(J,1)
         XFLUX(J2)=XFLUX(J2)-QDIS(J)*RDISQ(J,2)
       END DO
     END IF
   END IF


!----------PERFORM UPDATE ON ELF-----------------------------------------------!

   DTK = ALPHA_RK(K)*DTE
   ELF = ELRK - DTK*XFLUX/ART1
!!#  if defined (THIN_DAM)
!!   DO I=1,NODE_DAM1_N
!!     JN=I_NODE_DAM1_N(I,1)
!!     ELF(JN)=ELRK(JN)-DTK*(XFLUX(JN)+XFLUX(I_NODE_DAM1_N(I,2)))&
!!          &/(ART1(JN)+ART1(I_NODE_DAM1_N(I,2)))
!!     ELF(I_NODE_DAM1_N(I,2))=ELF(JN)
!!   END DO
!!   DO I=1,NODE_DAM2_N
!!     JN=I_NODE_DAM2_N(I,1)
!!     ELF(JN)=ELRK(JN)-DTK*(XFLUX(JN)+XFLUX(I_NODE_DAM2_N(I,2))&
!!          &+XFLUX(I_NODE_DAM2_N(I,3)) )&
!!          &/(ART1(JN)+ART1(I_NODE_DAM2_N(I,2))+ART1(I_NODE_DAM2_N(I,3)))
!!     ELF(I_NODE_DAM2_N(I,2))=ELF(JN)
!!     ELF(I_NODE_DAM2_N(I,3))=ELF(JN)
!!   END DO
!!   DO I=1,NODE_DAM3_N
!!     JN=I_NODE_DAM3_N(I,1)
!!     ELF(JN)=ELRK(JN)-DTK*(XFLUX(JN)+XFLUX(I_NODE_DAM3_N(I,2))&
!!          &+XFLUX(I_NODE_DAM3_N(I,3))+XFLUX(I_NODE_DAM3_N(I,4)) )&
!!          &/(ART1(JN)+ART1(I_NODE_DAM3_N(I,2))+ART1(I_NODE_DAM3_N(I&
!!          &,3))+ART1(I_NODE_DAM3_N(I,4)))
!!     ELF(I_NODE_DAM3_N(I,2))=ELF(JN)
!!     ELF(I_NODE_DAM3_N(I,3))=ELF(JN)
!!     ELF(I_NODE_DAM3_N(I,4))=ELF(JN)
!!   END DO
!!#  endif

!
!--STORE VARIABLES FOR MOMENTUM BALANCE CHECK----------------------------------|
!
#  if defined (BALANCE_2D)

   DO I=1,N
     XFLUXU1(I)=ONE_THIRD*(XFLUXU(NV(I,1))+ XFLUXU(NV(I,2))+ XFLUXU(NV(I,3)))
     XFLUXV1(I)=ONE_THIRD*(XFLUXV(NV(I,1))+ XFLUXV(NV(I,2))+ XFLUXV(NV(I,3)))
   END DO
   
   IF(K == 4) THEN
     DIVX2D2 = DIVX2D2 + XFLUXU1/ART/FLOAT(ISPLIT)            !dUD/dx
     DIVY2D2 = DIVY2D2 + XFLUXV1/ART/FLOAT(ISPLIT)            !dVD/dy
     DEDT2   = DEDT2 + (ELF1-ELRK1*(H1+ELRK1)/(H1+ELF1))/DTE/FLOAT(ISPLIT)    
   END IF     
#  endif
  
   IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END: extel_edge.F"

#  endif
   RETURN
   END SUBROUTINE EXTEL_EDGE
!==============================================================================|
