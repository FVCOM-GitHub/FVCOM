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
!   Calculate Advection and Horizontal Diffusion Terms for Temperature         |
!==============================================================================|

SUBROUTINE ADV_T_RK(TB1)               

  !------------------------------------------------------------------------------|

  USE ALL_VARS
  USE MOD_UTILS
  USE MOD_PAR
  USE BCS
  USE MOD_OBCS
  USE MOD_WD 
  USE MOD_SPHERICAL
  USE MOD_NORTHPOLE
# if defined (THIN_DAM)
  USE MOD_DAM,only : kdam,N_DAM_MATCH,IS_DAM
# endif

# if defined (TVD)
  USE MOD_TVD
# endif

  IMPLICIT NONE
  REAL(SP), DIMENSION(0:MT,KB)      :: XFLUX,XFLUX_ADV,RF
  REAL(SP), DIMENSION(0:MT)         :: PUPX,PUPY,PVPX,PVPY  
  REAL(SP), DIMENSION(0:MT)         :: PTPX,PTPY,PTPXD,PTPYD,VISCOFF
  REAL(SP), DIMENSION(3*(NT),KBM1)  :: DTIJ 
  REAL(SP), DIMENSION(3*(NT),KBM1)  :: UVN
  REAL(SP) :: UTMP,VTMP,SITAI,FFD,FF1 !,X11,Y11,X22,Y22,X33,Y33,TMP1,TMP2,XI,YI
  REAL(SP) :: DXA,DYA,DXB,DYB,FIJ1,FIJ2,UN,TTIME,ZDEP
  REAL(SP) :: TXX,TYY,FXX,FYY,VISCOF,EXFLUX,TEMP,STPOINT,STPOINT1,STPOINT2
  REAL(SP) :: FACT,FM1
  INTEGER  :: I,I1,I2,IA,IB,J,J1,J2,K,JTMP,JJ,II
  REAL(SP) :: T1MIN, T1MAX, T2MIN, T2MAX

# if defined (TVD)
  REAL(SP), DIMENSION(1:M)             :: dstA,dstB
  REAL(SP), DIMENSION(1:NCV)           :: XUA,XUB,YUA,YUB
  REAL(SP) :: LX,LY,minA,minB,T_upstreamA,T_upstreamB,A_limiter
  REAL(SP) :: B_limiter,smoothrA,smoothrB
# endif

# if defined (LIMITER_VER_ADV)
   REAL(SP)  CONV_T(1:KB), DISS_T(1:KB)
   REAL(SP)  SL_H(0:KB), T_TEMP(0:KB)
   REAL(SP)  SL_U, SL_F
# endif

# if defined (THIN_DAM)
  INTEGER :: NX
# endif

  REAL(SP), DIMENSION(0:MT,KB)     :: TB1      !! temporary temperature in RK

  !------------------------------------------------------------------------------!

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
  XFLUX     = 0.0_SP
  XFLUX_ADV = 0.0_SP

  !
  !--Loop Over Control Volume Sub-Edges And Calculate Normal Velocity------------!
  !
  DO I=1,NCV
     I1=NTRG(I)
     DO K=1,KBM1
        DTIJ(I,K) = DT1(I1)*DZ1(I1,K)
        ! USE U,V
        UVN(I,K)  = V(I1,K)*DLTXE(I) - U(I1,K)*DLTYE(I) 
     END DO
  END DO

  !
  !--Add the Shortwave Radiation Body Force--------------------------------------!
  !
  RF = 0.0_sp
  IF(HEATING_ON) THEN
     SELECT CASE(HEATING_TYPE)
     CASE('body')

        DO  K=1,KBM1
           DO  I=1,M
              ZDEP=0.5_SP*(Z(I,K)+Z(I,K+1))*DT(I)
              RF(I,K)=-SWRAD(I)*((RHEAT/ZETA1)*EXP(ZDEP/ZETA1) &
                   +((1-RHEAT)/ZETA2)*EXP(ZDEP/ZETA2))*DT(I)*DZ(I,K)
           END DO
        END DO

     CASE('flux')
        RF = 0.0_SP
     CASE('surface')
        RF = 0.0_SP
     CASE DEFAULT
        CALL FATAL_ERROR('The surface heating type is set incorrectly:',&
             & TRIM(HEATING_TYPE))
     END SELECT
  END IF

!--ADJUST VOLUME'S HEAT CONTENT FOR EVAPORATION AND PRECIPITATION ------------------!
  IF (PRECIPITATION_ON) THEN
     RF(:,1)=RF(:,1)+ROFVROS*(QEVAP+QPREC)*T1(:,1)
  END IF

  !
  !--Calculate the Advection and Horizontal Diffusion Terms----------------------!
  !

  DO K=1,KBM1
     PTPX  = 0.0_SP
     PTPY  = 0.0_SP
     PTPXD = 0.0_SP
     PTPYD = 0.0_SP
     DO I=1,M
        DO J=1,NTSN(I)-1
           I1=NBSN(I,J)
           I2=NBSN(I,J+1)
	 
#    if defined (WET_DRY)
           IF(ISWETN(I1) == 0 .AND. ISWETN(I2) == 1)THEN
            FFD=0.5_SP*(T1(I,K)+T1(I2,K)-TMEAN1(I,K)-TMEAN1(I2,K))
            FF1=0.5_SP*(T1(I,K)+T1(I2,K))
       	   ELSE IF(ISWETN(I1) == 1 .AND. ISWETN(I2) == 0)THEN
            FFD=0.5_SP*(T1(I1,K)+T1(I,K)-TMEAN1(I1,K)-TMEAN1(I,K))
            FF1=0.5_SP*(T1(I1,K)+T1(I,K))
	   ELSE IF(ISWETN(I1) == 0 .AND. ISWETN(I2) == 0)THEN
            FFD=0.5_SP*(T1(I,K)+T1(I,K)-TMEAN1(I,K)-TMEAN1(I,K))
            FF1=0.5_SP*(T1(I,K)+T1(I,K))
	   ELSE
            FFD=0.5_SP*(T1(I1,K)+T1(I2,K)-TMEAN1(I1,K)-TMEAN1(I2,K))
            FF1=0.5_SP*(T1(I1,K)+T1(I2,K))
	   END IF 
#    else	 
         FFD=0.5_SP*(T1(I1,K)+T1(I2,K)-TMEAN1(I1,K)-TMEAN1(I2,K))
         FF1=0.5_SP*(T1(I1,K)+T1(I2,K))
#    endif	 
	 
           PTPX(I)=PTPX(I)+FF1*DLTYTRIE(i,j)
           PTPY(I)=PTPY(I)+FF1*DLTXTRIE(i,j)
           PTPXD(I)=PTPXD(I)+FFD*DLTYTRIE(i,j)
           PTPYD(I)=PTPYD(I)+FFD*DLTXTRIE(i,j)

        END DO
! gather all neighboring control volumes connecting at dam node 
# if defined (THIN_DAM)
       IF(IS_DAM(I)==1.AND.K<=KDAM(I))THEN
         DO NX=1,N_DAM_MATCH(I,1)
           DO J=1,NTSN(N_DAM_MATCH(I,NX+1))-1
             I1=NBSN(N_DAM_MATCH(I,NX+1),J)
             I2=NBSN(N_DAM_MATCH(I,NX+1),J+1)
               FFD=0.5_SP*(T1(I1,K)+T1(I2,K)-TMEAN1(I1,K)-TMEAN1(I2,K))
               FF1=0.5_SP*(T1(I1,K)+T1(I2,K))
             PTPX(I)=PTPX(I)+FF1*DLTYTRIE(N_DAM_MATCH(I,NX+1),j)
             PTPY(I)=PTPY(I)+FF1*DLTXTRIE(N_DAM_MATCH(I,NX+1),j)
             PTPXD(I)=PTPXD(I)+FFD*DLTYTRIE(N_DAM_MATCH(I,NX+1),j)
             PTPYD(I)=PTPYD(I)+FFD*DLTXTRIE(N_DAM_MATCH(I,NX+1),j)
           END DO
         END DO
       END IF
# endif
# if !defined (THIN_DAM)
        PTPX(I)=PTPX(I)/ART2(I)
        PTPY(I)=PTPY(I)/ART2(I)
        PTPXD(I)=PTPXD(I)/ART2(I)
        PTPYD(I)=PTPYD(I)/ART2(I)
# else
        IF(IS_DAM(I)==1.AND.K<=KDAM(I))THEN
          PTPX(I)=PTPX(I)/(ART2(I)+SUM(ART2(N_DAM_MATCH(I,2:1+N_DAM_MATCH(I,1)))))
          PTPY(I)=PTPY(I)/(ART2(I)+SUM(ART2(N_DAM_MATCH(I,2:1+N_DAM_MATCH(I,1)))))
          PTPXD(I)=PTPXD(I)/(ART2(I)+SUM(ART2(N_DAM_MATCH(I,2:1+N_DAM_MATCH(I,1)))))
          PTPYD(I)=PTPYD(I)/(ART2(I)+SUM(ART2(N_DAM_MATCH(I,2:1+N_DAM_MATCH(I,1)))))
        ELSE
          PTPX(I)=PTPX(I)/ART2(I)
          PTPY(I)=PTPY(I)/ART2(I)
          PTPXD(I)=PTPXD(I)/ART2(I)
          PTPYD(I)=PTPYD(I)/ART2(I)
        END IF
        
# endif

     END DO

     IF(K == KBM1)THEN
        DO I=1,M
           PFPXB(I) = PTPX(I)
           PFPYB(I) = PTPY(I)
        END DO
     END IF

     DO I=1,M

        VISCOFF(I) = VISCOFH(I,K)

     END DO
     IF(K == KBM1) THEN
        AH_BOTTOM(1:M) = (FACT*VISCOFF(1:M) + FM1) * NN_HVC(1:M)
     END IF


     DO I=1,NCV_I
        IA=NIEC(I,1)
        IB=NIEC(I,2)

# if !defined (TVD)
          FIJ1=T1(IA,K)+DLTXNCVE(I,1)*PTPX(IA)+DLTYNCVE(I,1)*PTPY(IA)
          FIJ2=T1(IB,K)+DLTXNCVE(I,2)*PTPX(IB)+DLTYNCVE(I,2)*PTPY(IB)

        T1MIN=MINVAL(T1(NBSN(IA,1:NTSN(IA)-1),K))
        T1MIN=MIN(T1MIN, T1(IA,K))
        T1MAX=MAXVAL(T1(NBSN(IA,1:NTSN(IA)-1),K))
        T1MAX=MAX(T1MAX, T1(IA,K))
        T2MIN=MINVAL(T1(NBSN(IB,1:NTSN(IB)-1),K))
        T2MIN=MIN(T2MIN, T1(IB,K))
        T2MAX=MAXVAL(T1(NBSN(IB,1:NTSN(IB)-1),K))
        T2MAX=MAX(T2MAX, T1(IB,K))
        IF(FIJ1 < T1MIN) FIJ1=T1MIN
        IF(FIJ1 > T1MAX) FIJ1=T1MAX
        IF(FIJ2 < T2MIN) FIJ2=T2MIN
        IF(FIJ2 > T2MAX) FIJ2=T2MAX
# else
! ------------------------------------------------------------------------------------------ !
!                                 Drawing the TVD scheme
! ------------------------------------------------------------------------------------------ !
!     If A is upstream of B
       T_upstreamA=T1(Anear_node(I),K)+PTPX(Anear_node(I))*XUAdist(I)+PTPY(Anear_node(I))*YUAdist(I)
       T_upstreamB=T1(Bnear_node(I),K)+PTPX(Bnear_node(I))*XUBdist(I)+PTPY(Bnear_node(I))*YUBdist(I)                

      IF (T1(IA,K).EQ.T1(IB,K)) THEN
          smoothrA    = 0.0_SP
          smoothrB    = 0.0_SP
         ELSE IF (T1(IA,K).LT.1.E-10.AND.T1(IB,K).LT.1.E-10) THEN
          smoothrA    = 0.0_SP
          smoothrB    = 0.0_SP
         ELSE
!         T_upstreamA is the "far upstream node", T1(IA,K) is upstream of T1(IB,K), vica versa for "B"
          smoothrA    = (T1(IA,K)-T_upstreamA)/(T1(IB,K)-T1(IA,K))
          smoothrB    = (T1(IB,K)-T_upstreamB)/(T1(IA,K)-T1(IB,K))
       END IF

!     The Flux-limiter 
!      (superbee)
       A_limiter = MAX(0.0_SP, MIN(1.0_SP, 2.0_SP*smoothrA), MIN(2.0_SP, smoothrA))
       B_limiter = MAX(0.0_SP, MIN(1.0_SP, 2.0_SP*smoothrB), MIN(2.0_SP, smoothrB))

!     (minmod)
!      A_limiter = MAX(0.0_SP,MIN(1.0_SP,smoothrA))
!      B_limiter = MAX(0.0_SP,MIN(1.0_SP,smoothrB))

!      Calculating the salinity in-between IA and IB
!       T12_A     = T1(IA,K)+0.5_SP*A_LIMITER*(T1(IB,K)-T1(IA,K))
!       T12_B     = T1(IB,K)+0.5_SP*B_LIMITER*(T1(IA,K)-T1(IB,K))
       FIJ1     = T1(IA,K)+0.5_SP*A_LIMITER*(T1(IB,K)-T1(IA,K))
       FIJ2     = T1(IB,K)+0.5_SP*B_LIMITER*(T1(IA,K)-T1(IB,K))
# endif

        UN=UVN(I,K)

        ! David moved HPRNU and added HVC
        VISCOF=(FACT*0.5_SP*(VISCOFF(IA)*NN_HVC(IA)+VISCOFF(IB)*NN_HVC(IB)) + FM1*0.5_SP*(NN_HVC(IA)+NN_HVC(IB)))


        TXX=0.5_SP*(PTPXD(IA)+PTPXD(IB))*VISCOF
        TYY=0.5_SP*(PTPYD(IA)+PTPYD(IB))*VISCOF

        FXX=-DTIJ(I,K)*TXX*DLTYE(I)
        FYY= DTIJ(I,K)*TYY*DLTXE(I)

        EXFLUX=-UN*DTIJ(I,K)* &
             ((1.0_SP+SIGN(1.0_SP,UN))*FIJ2+(1.0_SP-SIGN(1.0_SP,UN))*FIJ1)*0.5_SP+FXX+FYY

        XFLUX(IA,K)=XFLUX(IA,K)+EXFLUX
        XFLUX(IB,K)=XFLUX(IB,K)-EXFLUX

        XFLUX_ADV(IA,K)=XFLUX_ADV(IA,K)+(EXFLUX-FXX-FYY)
        XFLUX_ADV(IB,K)=XFLUX_ADV(IB,K)-(EXFLUX-FXX-FYY)

#      if defined (THIN_DAM)
       IF(K<=KDAM(IA).AND.IS_DAM(IA)==1)THEN
         IF(N_DAM_MATCH(IA,1)==1)THEN
           XFLUX(N_DAM_MATCH(IA,2),K) = XFLUX(N_DAM_MATCH(IA,2),K) + EXFLUX
           XFLUX_ADV(N_DAM_MATCH(IA,2),K) = XFLUX_ADV(N_DAM_MATCH(IA,2),K) +(EXFLUX-FXX-FYY)
         END IF
         IF(N_DAM_MATCH(IA,1)==2)THEN
           XFLUX(N_DAM_MATCH(IA,2),K) = XFLUX(N_DAM_MATCH(IA,2),K) + EXFLUX
           XFLUX(N_DAM_MATCH(IA,3),K) = XFLUX(N_DAM_MATCH(IA,3),K) + EXFLUX
           XFLUX_ADV(N_DAM_MATCH(IA,2),K) = XFLUX_ADV(N_DAM_MATCH(IA,2),K) +(EXFLUX-FXX-FYY)
           XFLUX_ADV(N_DAM_MATCH(IA,3),K) = XFLUX_ADV(N_DAM_MATCH(IA,3),K) +(EXFLUX-FXX-FYY)
         END IF
         IF(N_DAM_MATCH(IA,1)==3)THEN
           XFLUX(N_DAM_MATCH(IA,2),K) = XFLUX(N_DAM_MATCH(IA,2),K) + EXFLUX
           XFLUX(N_DAM_MATCH(IA,3),K) = XFLUX(N_DAM_MATCH(IA,3),K) + EXFLUX
           XFLUX(N_DAM_MATCH(IA,4),K) = XFLUX(N_DAM_MATCH(IA,4),K) + EXFLUX
           XFLUX_ADV(N_DAM_MATCH(IA,2),K) = XFLUX_ADV(N_DAM_MATCH(IA,2),K) +(EXFLUX-FXX-FYY)
           XFLUX_ADV(N_DAM_MATCH(IA,3),K) = XFLUX_ADV(N_DAM_MATCH(IA,3),K) +(EXFLUX-FXX-FYY)
           XFLUX_ADV(N_DAM_MATCH(IA,4),K) = XFLUX_ADV(N_DAM_MATCH(IA,4),K) +(EXFLUX-FXX-FYY)
         END IF
       END IF
       IF(K<=KDAM(IB).AND.IS_DAM(IB)==1)THEN
         IF(N_DAM_MATCH(IB,1)==1)THEN
           XFLUX(N_DAM_MATCH(IB,2),K) = XFLUX(N_DAM_MATCH(IB,2),K) - EXFLUX
           XFLUX_ADV(N_DAM_MATCH(IB,2),K) = XFLUX_ADV(N_DAM_MATCH(IB,2),K) - (EXFLUX-FXX-FYY)
         END IF
         IF(N_DAM_MATCH(IB,1)==2)THEN
           XFLUX(N_DAM_MATCH(IB,2),K) = XFLUX(N_DAM_MATCH(IB,2),K) - EXFLUX
           XFLUX(N_DAM_MATCH(IB,3),K) = XFLUX(N_DAM_MATCH(IB,3),K) - EXFLUX
           XFLUX_ADV(N_DAM_MATCH(IB,2),K) = XFLUX_ADV(N_DAM_MATCH(IB,2),K) - (EXFLUX-FXX-FYY)
           XFLUX_ADV(N_DAM_MATCH(IB,3),K) = XFLUX_ADV(N_DAM_MATCH(IB,3),K) - (EXFLUX-FXX-FYY)
         END IF
         IF(N_DAM_MATCH(IB,1)==3)THEN
           XFLUX(N_DAM_MATCH(IB,2),K) = XFLUX(N_DAM_MATCH(IB,2),K) - EXFLUX
           XFLUX(N_DAM_MATCH(IB,3),K) = XFLUX(N_DAM_MATCH(IB,3),K) - EXFLUX
           XFLUX(N_DAM_MATCH(IB,4),K) = XFLUX(N_DAM_MATCH(IB,4),K) - EXFLUX
           XFLUX_ADV(N_DAM_MATCH(IB,2),K) = XFLUX_ADV(N_DAM_MATCH(IB,2),K) - (EXFLUX-FXX-FYY)
           XFLUX_ADV(N_DAM_MATCH(IB,3),K) = XFLUX_ADV(N_DAM_MATCH(IB,3),K) - (EXFLUX-FXX-FYY)
           XFLUX_ADV(N_DAM_MATCH(IB,4),K) = XFLUX_ADV(N_DAM_MATCH(IB,4),K) - (EXFLUX-FXX-FYY)
         END IF
       END IF
#      endif

     END DO

#    if defined (SPHERICAL)
       CALL ADV_T_XY(XFLUX,XFLUX_ADV,PTPX,PTPY,PTPXD,PTPYD,VISCOFF,K,0.0_SP)
#    endif  

  END DO !! K LOOP

# if defined (MULTIPROCESSOR)
  IF(PAR)CALL NODE_MATCH(0,NBN,BN_MLT,BN_LOC,BNC,MT,KB,MYID,NPROCS,XFLUX,XFLUX_ADV)
# endif

  DO K=1,KBM1
     IF(IOBCN > 0) THEN
        DO I=1,IOBCN
           I1=I_OBC_N(I)
           XFLUX_OBC(I,K)=XFLUX_ADV(I1,K)
        END DO
     END IF
  END DO

  !--------------------------------------------------------------------
  !   The central difference scheme in vertical advection
  !--------------------------------------------------------------------
  DO I=1, M
#    if defined (WET_DRY)
     IF(ISWETN(I)*ISWETNT(I) == 1) THEN
#    endif

#       if defined (LIMITER_VER_ADV)
        T_TEMP(0)  = -T1(I,1)
        T_TEMP(KB) = -T1(I,KBM1)
        SL_H(0)    = DZ(I,1)
        SL_H(KB)   = DZ(I,KBM1)
        DO K=1, KBM1
          T_TEMP(K) = T1(I,K)
          SL_H(K)   = DZ(I,K)
        ENDDO

        DO K=2, KBM1
          CONV_T(K) = WTS(I,K)*(T_TEMP(K)+T_TEMP(K-1))*0.5_SP
          SL_U = 2.0_SP*(T_TEMP(K)-T_TEMP(K+1))/(SL_H(K)+SL_H(K+1))
          SL_F = 2.0_SP*(T_TEMP(K-2)-T_TEMP(K-1))/(SL_H(K-2)+SL_H(K-1))
          DISS_T(K) = 0.5_SP*ABS(WTS(I,K))*(T_TEMP(K-1)-T_TEMP(K)-0.5_SP*LIMLED2(SL_U,SL_F,2.0_SP)*(SL_H(K-1)+SL_H(K)))
        ENDDO
        CONV_T(1)  = 0.0_SP
        DISS_T(1)  = 0.0_SP
        CONV_T(KB) = 0.0_SP
        DISS_T(KB) = 0.0_SP
#       endif

        DO K=1, KBM1

#          if defined (LIMITER_VER_ADV)
           TEMP = CONV_T(K)-CONV_T(K+1)+DISS_T(K+1)-DISS_T(K)
#          endif

#          if !defined (LIMITER_VER_ADV)

             IF(K == 1) THEN
              TEMP=-WTS(I,K+1)*(T1(I,K)*DZ(I,K+1)+T1(I,K+1)*DZ(I,K))/   &
                   (DZ(I,K)+DZ(I,K+1))
             ELSE IF(K == KBM1) THEN
              TEMP= WTS(I,K)*(T1(I,K)*DZ(I,K-1)+T1(I,K-1)*DZ(I,K))/(DZ(I,K)+DZ(I,K-1))
             ELSE
              TEMP= WTS(I,K)*(T1(I,K)*DZ(I,K-1)+T1(I,K-1)*DZ(I,K))/(DZ(I,K)+DZ(I,K-1))-&
                   WTS(I,K+1)*(T1(I,K)*DZ(I,K+1)+T1(I,K+1)*DZ(I,K))/(DZ(I,K)+DZ(I,K+1))
             END IF

#          endif

           IF(ISONB(I) == 2) THEN
              !         XFLUX(I,K)=TEMP*ART1(I)/DZ(I,K)
              XFLUX(I,K)=TEMP*ART1(I)
           ELSE
              !         XFLUX(I,K)=XFLUX(I,K)+TEMP*ART1(I)/DZ(I,K)
              XFLUX(I,K)=XFLUX(I,K)+TEMP*ART1(I)
#          if defined (THIN_DAM)
              IF(IS_DAM(I)==1.AND.K<=KDAM(I))THEN
                IF(N_DAM_MATCH(I,1)==1)THEN
                  XFLUX(N_DAM_MATCH(I,2),K) = XFLUX(N_DAM_MATCH(I,2),K)+TEMP*ART1(I)
                END IF
                IF(N_DAM_MATCH(I,1)==2)THEN
                  XFLUX(N_DAM_MATCH(I,2),K) = XFLUX(N_DAM_MATCH(I,2),K)+TEMP*ART1(I) 
                  XFLUX(N_DAM_MATCH(I,3),K) = XFLUX(N_DAM_MATCH(I,3),K)+TEMP*ART1(I)
                END IF
                IF(N_DAM_MATCH(I,1)==3)THEN
                  XFLUX(N_DAM_MATCH(I,2),K) = XFLUX(N_DAM_MATCH(I,2),K)+TEMP*ART1(I)
                  XFLUX(N_DAM_MATCH(I,3),K) = XFLUX(N_DAM_MATCH(I,3),K)+TEMP*ART1(I)
                  XFLUX(N_DAM_MATCH(I,4),K) = XFLUX(N_DAM_MATCH(I,4),K)+TEMP*ART1(I)
                END IF                
              END IF
#          endif

           END IF
        ENDDO
#    if defined (WET_DRY)
     END IF
#    endif
  END DO  !!K LOOP
  !
  !--Set Boundary Conditions-For Fresh Water Flux--------------------------------!
  !
  IF(RIVER_TS_SETTING == 'calculated') THEN
     IF(RIVER_INFLOW_LOCATION == 'node') THEN
        IF(NUMQBC > 0) THEN
           DO J=1,NUMQBC
              JJ=INODEQ(J)
              STPOINT=TDIS(J)  
              DO K=1,KBM1
              !   STPOINT    = T1(JJ,K)
                 !             XFLUX(JJ,K)= XFLUX(JJ,K)-QDIS(J)*VQDIST(J,K)*STPOINT/DZ(JJ,K)
                 XFLUX(JJ,K)= XFLUX(JJ,K)-QDIS(J)*VQDIST(J,K)*STPOINT
              END DO
           END DO
        END IF
     ELSE IF(RIVER_INFLOW_LOCATION == 'edge') THEN
        IF(NUMQBC > 0) THEN
           DO J=1,NUMQBC
              J1=N_ICELLQ(J,1)
              J2=N_ICELLQ(J,2)
              STPOINT=TDIS(J) 
              DO K=1,KBM1
                 !STPOINT1 = T1(J1,K)
                 !STPOINT2 = T1(J2,K)
                 XFLUX(J1,K)=XFLUX(J1,K)-QDIS(J)*RDISQ(J,1)*VQDIST(J,K)*STPOINT     !1/DZ1(J1,K)
                 XFLUX(J2,K)=XFLUX(J2,K)-QDIS(J)*RDISQ(J,2)*VQDIST(J,K)*STPOINT     !2/DZ1(J2,K)
              END DO
           END DO
        END IF
     END IF
  END IF
  !---------------------------------------------------------------------

  ! APPLY GROUND WATER TEMPERATURE FORCING
  IF(GROUNDWATER_ON .and. GROUNDWATER_TEMP_ON)THEN
     DO I=1,M
        XFLUX(I,KBM1)=XFLUX(I,KBM1)-BFWDIS(I)*BFWTMP(I)
     END DO
  ELSEIF(GROUNDWATER_ON) THEN
     DO I=1,M
          XFLUX(I,KBM1)=XFLUX(I,KBM1)-BFWDIS(I)*T1(I,KBM1)
     END DO
  END IF


  !
  !--Update Temperature----------------------------------------------------------!
  !


#  if defined (WET_DRY)
  DO I=1,M
     IF(ISWETN(I)*ISWETNT(I) == 1 )THEN
        DO K=1,KBM1
           XFLUX(I,K) = XFLUX(I,K) - RF(I,K)*ART1(I)    !/DZ(I,K)
#      if !defined (THIN_DAM)
           TF1(I,K)=(TB1(I,K)-XFLUX(I,K)/ART1(I)*(DTI/(DT(I)*DZ(I,K))))*(DT(I)/DTFA(I))
#      else
           IF(IS_DAM(I)==1.AND.K<=KDAM(I))THEN
             TF1(I,K)=(TB1(I,K)-XFLUX(I,K)/(ART1(I)&
                  &+SUM(ART1(N_DAM_MATCH(I,2:1+N_DAM_MATCH(I,1)))))*(DTI/(DT(I)*DZ(I,K))))*(DT(I)/DTFA(I))
           ELSE
             TF1(I,K)=(TB1(I,K)-XFLUX(I,K)/ART1(I)*(DTI/(DT(I)*DZ(I,K))))*(DT(I)/DTFA(I))
           END IF
#      endif
        END DO
     ELSE
        DO K=1,KBM1
           TF1(I,K)=TB1(I,K)
        END DO
     END IF
  END DO
# else
  DO I=1,M
     DO K=1,KBM1
        XFLUX(I,K) = XFLUX(I,K) - RF(I,K)*ART1(I)    !/DZ(I,K)
#      if !defined (THIN_DAM)
        TF1(I,K)=(TB1(I,K)-XFLUX(I,K)/ART1(I)*(DTI/(DT(I)*DZ(I,K))))*(DT(I)/DTFA(I))
#      else
        IF(IS_DAM(I)==1.AND.K<=KDAM(I))THEN
          TF1(I,K)=(TB1(I,K)-XFLUX(I,K)/(ART1(I)&
               &+SUM(ART1(N_DAM_MATCH(I,2:1+N_DAM_MATCH(I,1)))))*(DTI/(DT(I)*DZ(I,K))))*(DT(I)/DTFA(I))
        ELSE
          TF1(I,K)=(TB1(I,K)-XFLUX(I,K)/ART1(I)*(DTI/(DT(I)*DZ(I,K))))*(DT(I)/DTFA(I))
        END IF
#      endif
     END DO
  END DO
# endif



  IF(DBG_SET(DBG_SBR)) WRITE(IPT,*)"End: adv_t_rk"

END SUBROUTINE ADV_T_RK
!==============================================================================|
