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
!    This subroutine is used to calculate the area of individual               !
!    triangle based on the three vertex coordinates and also calculate         !
!    the sigma-surface area of individual control volume consisted of          !
!    triangles with a common node point                                        !
!									       !
! calculates: art(ntri)   = area of element (triangle) 			       !
! calculates: art1(nnode) = area of interior cv (for node value integration)   !
! calculates: art2(nnode) = sum area of all cells around node		       !
!==============================================================================|

   SUBROUTINE CELL_AREA 

!==============================================================================!
   USE ALL_VARS
   USE MOD_UTILS
   USE MOD_PAR
   USE MOD_SPHERICAL
   IMPLICIT NONE
   REAL(SP), ALLOCATABLE :: XX(:),YY(:) 
   REAL(SP) :: ARTMAX,ARTTOT,ARTMIN
   REAL(SP) :: ART1MAX,ART1TOT,ART1MIN
   INTEGER  :: I,J,II,J1,J2,MAX_NBRE
   
   REAL(SP) :: SBUF
   INTEGER  :: IERR
   CHARACTER(LEN=10) :: TSTR
   CHARACTER(LEN=80) ::MSG

#  if defined (SPHERICAL)   
   REAL(DP) :: VX1,VX2,VX3, VY1,VY2,VY3,side1,side2,side3,area1               
   REAL(DP) :: xx1,yy1,xx2,yy2,XXC,YYC,&
               X1_DP, Y1_DP, X2_DP, Y2_DP, X3_DP, Y3_DP
#  endif   
!==============================================================================!

    IF (DBG_SET(DBG_SBR)) WRITE(IPT,*) "START: CELL_AREA"
    
!---------------INITIALIZE ARRAYS----------------------------------------------!

   ART  = 0.0_SP ; ART1 = 0.0_SP ; ART2 = 0.0_SP
   MAX_NBRE = MAXVAL(NTVE)+1
   ALLOCATE(XX(2*MAX_NBRE+1),YY(2*MAX_NBRE+1))
   XX = 0.0_SP ; YY = 0.0_SP

!---------------COMPUTE AREA OF TRIANGLES USING CROSS PRODUCT------------------!
  
#  if defined (SPHERICAL)
   DO I=1,NT 
     VX1=VX(NV(I,1))
     VX2=VX(NV(I,2))
     VX3=VX(NV(I,3))
     VY1=VY(NV(I,1))
     VY2=VY(NV(I,2))
     VY3=VY(NV(I,3))
     CALL ARC(VX2,VY2,VX3,VY3,side1)
     CALL ARC(VX3,VY3,VX1,VY1,side2)
     CALL ARC(VX1,VY1,VX2,VY2,side3)
     CALL AREA(side1,side2,side3,area1)
     ART(I)=AREA1       
   ENDDO
#  else
    DO I=1,NT 
    ART(I) = (VX(NV(I,2)) - VX(NV(I,1))) * (VY(NV(I,3)) - VY(NV(I,1))) - &
             (VX(NV(I,3)) - VX(NV(I,1))) * (VY(NV(I,2)) - VY(NV(I,1)))
    END DO
   ART    = ABS(.5_SP*ART)
#  endif    

!---------------COMPUTE MESH STATISTICS----------------------------------------!

   ARTMIN = MINVAL(ART(1:N))
   ARTMAX = MAXVAL(ART(1:N))
   ARTTOT =    SUM(ART(1:N))


   IF(artmin .LT. 1.0E-6_SP) THEN
      MSG = ""
      IF (PAR) THEN
         MSG = "Proc#"
         WRITE(tstr,'(I3)') MYID
         MSG=TRIM(MSG)//trim(TSTR)//"; "
      END IF

      MSG = TRIM(MSG)//"Min Triangle Area="
      WRITE(tstr,'(F9.6)') artmin
      MSG=TRIM(MSG)//trim(TSTR)
      
      I = MINLOC(ART(1:N),DIM=1)

      MSG = TRIM(MSG)//"; EGID="
      WRITE(tstr,'(I7)') EGID(I)
      MSG=TRIM(MSG)//trim(TSTR)

      WRITE(IPT,*) "*****************************"
      WRITE(IPT,*) TRIM(MSG)
      WRITE(IPT,*) "*****************************"
   END IF

#  if defined (MULTIPROCESSOR)
   IERR=0
   SBUF=ARTMIN
   IF(PAR)CALL MPI_ALLREDUCE(SBUF,ARTMIN,1,MPI_F,MPI_MIN,MPI_FVCOM_GROUP,IERR)
   SBUF=ARTMAX
   IF(PAR)CALL MPI_ALLREDUCE(SBUF,ARTMAX,1,MPI_F,MPI_MAX,MPI_FVCOM_GROUP,IERR)
   SBUF=ARTTOT
   IF(PAR)CALL MPI_ALLREDUCE(SBUF,ARTTOT,1,MPI_F,MPI_SUM,MPI_FVCOM_GROUP,IERR)
#  endif


   IF (DBG_SET(DBG_SCL)) THEN
      WRITE(IPT,*) "!     Minimum Triangle Area: ", ARTMIN
      WRITE(IPT,*) "!     Maximum Triangle Area: ", ARTMAX
      WRITE(IPT,*) "!     Total Triangle Area  : ", ARTTOT
   END IF
   IF(artmin.LT. 1.0E-6_SP) CALL WARNING("CELL_AREA: TRIANGLE AREA IS SMALL (LT 1e-6)")

!-------COMPUTE CONTROL VOLUME ART1: CV FOR FLUXES OF NODAL BASED VALUES-------!

   DO I=1,M
     IF(ISONB(I) == 0) THEN
       DO J=1,NTVE(I)
         II=NBVE(I,J)
         J1=NBVT(I,J)
         J2=J1+1-INT((J1+1)/4)*3
#        if defined (SPHERICAL)
         xx1=VX(NV(II,J1))
         yy1=VY(NV(II,J1))
         xx2=VX(NV(II,J2))
         yy2=VY(NV(II,J2))
         CALL ARCC(xx1,yy1,xx2,yy2,XXC,YYC) 
         XX(2*J-1)=XXC
         YY(2*J-1)=YYC           
         XX(2*J)=XC(II)
         YY(2*J)=YC(II)            
#        else
         XX(2*J-1)=(VX(NV(II,J1))+VX(NV(II,J2)))*0.5_SP-VX(I)
         YY(2*J-1)=(VY(NV(II,J1))+VY(NV(II,J2)))*0.5_SP-VY(I)
         XX(2*J)=XC(II)-VX(I)
         YY(2*J)=YC(II)-VY(I)
#        endif
       END DO
       XX(2*NTVE(I)+1)=XX(1)
       YY(2*NTVE(I)+1)=YY(1)

       DO J=1,2*NTVE(I)
#         if defined (SPHERICAL)
          X1_DP=XX(J)
          Y1_DP=YY(J)
          X2_DP=XX(J+1)
          Y2_DP=YY(J+1)
          X3_DP=VX(I)
          Y3_DP=VY(I)
          CALL ARC(X1_DP,Y1_DP,X2_DP,Y2_DP,side1)
          CALL ARC(X2_DP,Y2_DP,X3_DP,Y3_DP,side2)
          CALL ARC(X3_DP,Y3_DP,X1_DP,Y1_DP,side3)
          CALL AREA(side1,side2,side3,area1)        
          ART1(I)=ART1(I)+area1            
#         else
          ART1(I)=ART1(I)+0.5_SP*(XX(J+1)*YY(J)-XX(J)*YY(J+1))
#         endif
       END DO
       ART1(I)=ABS(ART1(I))
     ELSE
       DO J=1,NTVE(I)
         II=NBVE(I,J)
         J1=NBVT(I,J)
         J2=J1+1-INT((J1+1)/4)*3
#        if defined(SPHERICAL)
         xx1=VX(NV(II,J1))
         yy1=VY(NV(II,J1))
         xx2=VX(NV(II,J2))
         yy2=VY(NV(II,J2))
         CALL ARCC(xx1,yy1,xx2,yy2,XXC,YYC) 
                  
         XX(2*J-1)=XXC
         YY(2*J-1)=YYC             
             
         XX(2*J)=XC(II)
         YY(2*J)=YC(II) 
#        else
         XX(2*J-1)=(VX(NV(II,J1))+VX(NV(II,J2)))*0.5_SP-VX(I)
         YY(2*J-1)=(VY(NV(II,J1))+VY(NV(II,J2)))*0.5_SP-VY(I)
         XX(2*J)=XC(II)-VX(I)
         YY(2*J)=YC(II)-VY(I)
#        endif	                
       END DO
       J=NTVE(I)+1
       II=NBVE(I,J-1)
       J1=NBVT(I,NTVE(I))
       J2=J1+2-INT((J1+2)/4)*3
#      if defined (SPHERICAL)
       xx1 =  VX(NV(II,J1))
       yy1 =  VY(NV(II,J1))
       xx2 =  VX(NV(II,J2))
       yy2 =  VY(NV(II,J2))
       CALL ARCC(xx1,yy1,xx2,yy2,XXC,YYC) 
     
       XX(2*J-1)=XXC
       YY(2*J-1)=YYC 
          
       XX(2*J)=VX(I)
       YY(2*J)=VY(I)          

       XX(2*J+1)=XX(1)
       YY(2*J+1)=YY(1)

#      else
       XX(2*J-1)=(VX(NV(II,J1))+VX(NV(II,J2)))*0.5_SP-VX(I)
       YY(2*J-1)=(VY(NV(II,J1))+VY(NV(II,J2)))*0.5_SP-VY(I)

       XX(2*J)=VX(I)-VX(I)
       YY(2*J)=VY(I)-VY(I)

       XX(2*J+1)=XX(1)
       YY(2*J+1)=YY(1)
#      endif       

       DO J=1,2*NTVE(I)+2
#       if defined (SPHERICAL)
        X1_DP= XX(J)
        Y1_DP= YY(J)
        X2_DP= XX(J+1)
        Y2_DP= YY(J+1)
        X3_DP= VX(I)
        Y3_DP= VY(I)
        CALL ARC(X1_DP, Y1_DP, X2_DP, Y2_DP, side1)
        CALL ARC(X2_DP, Y2_DP, X3_DP, Y3_DP, side2)
        CALL ARC(X3_DP, Y3_DP, X1_DP, Y1_DP, side3)

        CALL AREA(side1,side2,side3,area1)
        
        ART1(I)=ART1(I)+area1             
#       else
        ART1(I)=ART1(I)+0.5_SP*(XX(J+1)*YY(J)-XX(J)*YY(J+1))
#       endif
       END DO
       ART1(I)=ABS(ART1(I))
     END IF
   ENDDO

!---------------COMPUTE MESH STATISTICS----------------------------------------!

   ART1MIN = MINVAL(ART1(1:M))
   ART1MAX = MAXVAL(ART1(1:M))
   ART1TOT =    SUM(ART1(1:M))

   IF(art1min .LT. 1.0E-6_SP) THEN
      MSG = ""
      IF (PAR) THEN
         MSG = "Proc#"
         WRITE(tstr,'(I3)') MYID
         MSG=TRIM(MSG)//trim(TSTR)//"; "
      END IF

      MSG = TRIM(MSG)//"Min Control Volume Area="
      WRITE(tstr,'(F9.6)') art1min
      MSG=TRIM(MSG)//trim(TSTR)
      
      I = MINLOC(ART1(1:M),DIM=1)

      MSG = TRIM(MSG)//"; NGID="
      WRITE(tstr,'(I7)') NGID(I)
      MSG=TRIM(MSG)//trim(TSTR)

      IF(ISONB(I)==0) THEN
         MSG = TRIM(MSG)//"; Node is interior"
      ELSEIF(ISONB(I)==1) THEN
         MSG = TRIM(MSG)//"; Node is on solid boundary"
      ELSEIF(ISONB(I)==2) THEN
         MSG = TRIM(MSG)//"; Node is on open boundary"
      ELSE
         MSG = TRIM(MSG)//"; ISONB has bad value!"
      END IF


      WRITE(IPT,*) "*****************************"
      WRITE(IPT,*) TRIM(MSG)
      WRITE(IPT,*) "*****************************"
   END IF

#  if defined (MULTIPROCESSOR)
   IERR=0

   SBUF=ART1MIN
   IF(PAR)CALL MPI_ALLREDUCE(SBUF,ART1MIN,1,MPI_F,MPI_MIN,MPI_FVCOM_GROUP,IERR)
   SBUF=ART1MAX
   IF(PAR)CALL MPI_ALLREDUCE(SBUF,ART1MAX,1,MPI_F,MPI_MAX,MPI_FVCOM_GROUP,IERR)
   SBUF=ART1TOT
   IF(PAR)CALL MPI_ALLREDUCE(SBUF,ART1TOT,1,MPI_F,MPI_SUM,MPI_FVCOM_GROUP,IERR)
#  endif


   IF (DBG_SET(DBG_SCL)) THEN
      WRITE(IPT,*) "!     Minimum Node Control Volume Area: ", ART1MIN
      WRITE(IPT,*) "!     Maximum Node Control Volume Area: ", ART1MAX
      WRITE(IPT,*) "!     Total Node Control Volume Area  : ", ART1TOT
   END IF
   IF(art1min.LT. 1.0E-6_SP) CALL WARNING(" CELL_AREA: NODAL CONTROL VOLUME IS SMALL (LT 1e-6)")




!---COMPUTE AREA OF CONTROL VOLUME ART2(I) = SUM(ALL TRIS SURROUNDING NODE I)--!

   DO I=1,M
     ART2(I) = SUM(ART(NBVE(I,1:NTVE(I))))
   END DO
   
   ART(0) = ART(1) 
   ART1(0) = ART1(1) 
   !   IF(NT > N)ART(N+1:NT) = ART(N)
   IF(MT > M)ART2(M+1:MT) = ART2(M)
   IF(MT > M)ART1(M+1:MT) = ART1(M)
   DEALLOCATE(XX,YY)
   
   ! NOTES: SHOULD MAKE AN ARRAY TO STORE 1/ART, 1/ART2 and 1/ART2
   ! IT is faster and safer

   IF (DBG_SET(DBG_LOG)) WRITE(IPT,*) "!  CELL AREA             :    COMPLETE"

   IF (DBG_SET(DBG_SBR)) WRITE(IPT,*) "END: CELL_AREA"
   RETURN
   END SUBROUTINE CELL_AREA
!==============================================================================|
