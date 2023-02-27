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
!   Calculate the Turbulent Kinetic Energy and Mixing Length Based  on         |
!   The Mellor-Yamada Level 2.5 Turbulent Closure Model                        |
!==============================================================================|

   SUBROUTINE ADV_Q_RK(Q,QB,QF)               

!------------------------------------------------------------------------------|
   USE MOD_UTILS
   USE ALL_VARS
   USE MOD_PAR
   USE MOD_WD
   USE MOD_SPHERICAL
   USE MOD_NORTHPOLE

   IMPLICIT NONE
   REAL(SP), DIMENSION(0:MT,KB)     :: Q,QB,QF,XFLUX
   REAL(SP), DIMENSION(0:MT)        :: PUPX,PUPY,PVPX,PVPY  
   REAL(SP), DIMENSION(0:MT)        :: PQPX,PQPY,PQPXD,PQPYD,VISCOFF
   REAL(SP), DIMENSION(3*(NT),KBM1) :: DTIJ 
   REAL(SP), DIMENSION(3*(NT),KBM1) :: UVN
   REAL(SP) :: UTMP,VTMP,SITAI,FFD,FF1 !,X11,Y11,X22,Y22,X33,Y33,TMP1,TMP2,XI,YI
   REAL(SP) :: DXA,DYA,DXB,DYB,FIJ1,FIJ2,UN
   REAL(SP) :: TXX,TYY,FXX,FYY,VISCOF,EXFLUX,TEMP,STPOINT
   REAL(SP) :: FACT,FM1
   INTEGER  :: I,I1,I2,IA,IB,J,J1,J2,K,JTMP,JJ,II
   REAL(SP) :: Q1MIN, Q1MAX, Q2MIN, Q2MAX

   REAL(SP) :: QMEAN1
   REAL(SP), DIMENSION(0:NT,KB)    :: UQ,VQ

   REAL(SP), ALLOCATABLE :: UQ1(:,:),VQ1(:,:)

   IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "Start: adv_q"

!------------------------------------------------------------------------------!

   QMEAN1 = 1.E-8

   ALLOCATE(UQ1(0,0))
   ALLOCATE(VQ1(0,0))   

   SELECT CASE(HORIZONTAL_MIXING_TYPE)
   CASE ('closure')
      FACT = 1.0_SP
      FM1  = 0.0_SP
   CASE('constant')
      FACT = 0.0_SP
      FM1  = 1.0_SP
   CASE DEFAULT
      CALL FATAL_ERROR("UNKNOW HORIZONTAL MIXING TYPE:",&
           & TRIM(HORIZONTAL_MIXING_TYPE) )
   END SELECT
   
!
!--Initialize Fluxes-----------------------------------------------------------!
!
   QF    = 0.0_SP
   XFLUX = 0.0_SP
   
   UQ = 0.0_SP
   VQ = 0.0_SP
   UVN = 0.0_SP
   
   DO K=2,KBM1
     DO I=1,NT
       UQ(I,K) = (U(I,K)*DZ1(I,K-1)+U(I,K-1)*DZ1(I,K))/(DZ1(I,K)+DZ1(I,K-1))
       VQ(I,K) = (V(I,K)*DZ1(I,K-1)+V(I,K-1)*DZ1(I,K))/(DZ1(I,K)+DZ1(I,K-1))
     END DO
   END DO     

!
!--Loop Over Control Volume Sub-Edges And Calculate Normal Velocity------------!
!
   DO I=1,NCV
     I1=NTRG(I)
     DO K=2,KBM1
       DTIJ(I,K)=DT1(I1)*DZZ1(I1,K-1)
       UVN(I,K) = VQ(I1,K)*DLTXE(I) - UQ(I1,K)*DLTYE(I)
     END DO
   END DO

!
!--Calculate the Advection and Horizontal Diffusion Terms----------------------!
!

   DO K=2,KBM1
     PQPX  = 0.0_SP 
     PQPY  = 0.0_SP 
     PQPXD = 0.0_SP 
     PQPYD = 0.0_SP
     DO I=1,M
       DO J=1,NTSN(I)-1
         I1=NBSN(I,J)
         I2=NBSN(I,J+1)

         FFD=0.5_SP*(Q(I1,K)+Q(I2,K)-QMEAN1-QMEAN1)
         FF1=0.5_SP*(Q(I1,K)+Q(I2,K))

         PQPX(I)=PQPX(I)+FF1*DLTYTRIE(i,j)
         PQPY(I)=PQPY(I)+FF1*DLTXTRIE(i,j)
         PQPXD(I)=PQPXD(I)+FFD*DLTYTRIE(i,j)
         PQPYD(I)=PQPYD(I)+FFD*DltXTRIE(i,j)

       END DO
       PQPX(I)=PQPX(I)/ART2(I)
       PQPY(I)=PQPY(I)/ART2(I)
       PQPXD(I)=PQPXD(I)/ART2(I)
       PQPYD(I)=PQPYD(I)/ART2(I)
     END DO

     DO I=1,M 
       VISCOFF(I) = (VISCOFH(I,K)*DZ(I,K-1)+VISCOFH(I,K-1)*DZ(I,K))/  &
                    (DZ(I,K)+DZ(I,K-1))  
     END DO

     DO I=1,NCV_I
       IA=NIEC(I,1)
       IB=NIEC(I,2)

        FIJ1=Q(IA,K)+DLTXNCVE(I,1)*PQPX(IA)+DLTYNCVE(I,1)*PQPY(IA)
        FIJ2=Q(IB,K)+DLTXNCVE(I,2)*PQPX(IB)+DLTYNCVE(I,2)*PQPY(IB)


       Q1MIN=MINVAL(Q(NBSN(IA,1:NTSN(IA)-1),K))
       Q1MIN=MIN(Q1MIN, Q(IA,K))
       Q1MAX=MAXVAL(Q(NBSN(IA,1:NTSN(IA)-1),K))
       Q1MAX=MAX(Q1MAX, Q(IA,K))
       Q2MIN=MINVAL(Q(NBSN(IB,1:NTSN(IB)-1),K))
       Q2MIN=MIN(Q2MIN, Q(IB,K))
       Q2MAX=MAXVAL(Q(NBSN(IB,1:NTSN(IB)-1),K))
       Q2MAX=MAX(Q2MAX, Q(IB,K))
       IF(FIJ1 < Q1MIN) FIJ1=Q1MIN
       IF(FIJ1 > Q1MAX) FIJ1=Q1MAX
       IF(FIJ2 < Q2MIN) FIJ2=Q2MIN
       IF(FIJ2 > Q2MAX) FIJ2=Q2MAX
    
       UN=UVN(I,K)

       ! David moved HPRNU and added HVC
       VISCOF=(FACT*0.5_SP*(VISCOFF(IA)*NN_HVC(IA)+VISCOFF(IB)*NN_HVC(IB)) + FM1*0.5_SP*(NN_HVC(IA)+NN_HVC(IB)))/HPRNU

       TXX=0.5_SP*(PQPXD(IA)+PQPXD(IB))*VISCOF
       TYY=0.5_SP*(PQPYD(IA)+PQPYD(IB))*VISCOF

       FXX=-DTIJ(I,K)*TXX*DLTYE(I)
       FYY= DTIJ(I,K)*TYY*DLTXE(I)

       EXFLUX=-UN*DTIJ(I,K)* &
          ((1.0_SP+SIGN(1.0_SP,UN))*FIJ2+(1.0_SP-SIGN(1.0_SP,UN))*FIJ1)*0.5_SP+FXX+FYY

       XFLUX(IA,K)=XFLUX(IA,K)+EXFLUX
       XFLUX(IB,K)=XFLUX(IB,K)-EXFLUX

     END DO

#    if defined (SPHERICAL)
     CALL ADV_Q_XY(XFLUX,PQPX,PQPY,PQPXD,PQPYD,VISCOFF,Q,UQ,VQ,K,UQ1,VQ1,0.0_SP)
#    endif  

   END DO !!SIGMA LOOP

!
!-Accumulate Fluxes at Boundary Nodes
!
# if defined (MULTIPROCESSOR)
   IF(PAR)CALL NODE_MATCH(0,NBN,BN_MLT,BN_LOC,BNC,MT,KB,MYID,NPROCS,XFLUX)
# endif

 
!--------------------------------------------------------------------
!   The central difference scheme in vertical advection
!--------------------------------------------------------------------
   DO K=2,KBM1
     DO I=1,M
#    if defined (WET_DRY)
       IF(ISWETN(I)*ISWETNT(I) == 1) THEN
#    endif
         TEMP=WTS(I,K-1)*Q(I,K-1)-WTS(I,K+1)*Q(I,K+1)
         XFLUX(I,K)=XFLUX(I,K)+TEMP*ART1(I)*DZZ(I,K-1)/(DZ(I,K-1)+DZ(I,K))
#    if defined (WET_DRY)
       END IF
#    endif
     END DO
   END DO  !! SIGMA LOOP

!
!--Update Q or QL-------------------------------------------------------------!
!

#  if defined (WET_DRY)
   DO I=1,M
     IF(ISWETN(I)*ISWETNT(I) == 1 )THEN
        DO K=2,KBM1
           QF(I,K)=(QB(I,K)-XFLUX(I,K)/ART1(I)*(DTI/(DT(I)*DZZ(I,k-1))))*(DT(I)/D(I))
        END DO
     ELSE
        DO K=2,KBM1
           QF(I,K)=QB(I,K)
        END DO
     END IF
   END DO
# else

   DO I=1,M
      DO K=2,KBM1
         QF(I,K)=(QB(I,K)-XFLUX(I,K)/ART1(I)*(DTI/(DT(I)*DZZ(I,k-1))))*(DT(I)/D(I))
      END DO
   END DO
#  endif

   IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "End: adv_q_rk"

   END SUBROUTINE ADV_Q_RK
!==============================================================================|
