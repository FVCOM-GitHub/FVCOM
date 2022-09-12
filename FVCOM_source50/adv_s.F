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
!   Calculate Advection and Horizontal Diffusion Terms for Salinity            |
!==============================================================================|

SUBROUTINE ADV_S               

  !------------------------------------------------------------------------------|

  USE ALL_VARS
  USE MOD_UTILS
  USE BCS
  USE MOD_OBCS
  USE MOD_PAR
  USE MOD_WD
  USE MOD_SPHERICAL
  USE MOD_NORTHPOLE
# if defined (SEMI_IMPLICIT)
  USE MOD_SEMI_IMPLICIT
# endif

# if defined (THIN_DAM)
  USE MOD_DAM, only : kdam,N_DAM_MATCH,IS_DAM
# endif

# if defined (TVD)
  USE MOD_TVD
# endif

  IMPLICIT NONE
  REAL(SP), DIMENSION(0:MT,KB)     :: XFLUX,XFLUX_ADV
  REAL(SP), DIMENSION(0:MT)        :: PUPX,PUPY,PVPX,PVPY  
  REAL(SP), DIMENSION(0:MT)        :: PSPX,PSPY,PSPXD,PSPYD,VISCOFF
  REAL(SP), DIMENSION(3*(NT),KBM1) :: DTIJ 
  REAL(SP), DIMENSION(3*(NT),KBM1) :: UVN
  REAL(SP) :: UTMP,VTMP,SITAI,FFD,FF1 !,X11,Y11,X22,Y22,X33,Y33 !,TMP1,TMP2 !,XI,YI
  REAL(SP) :: DXA,DYA,DXB,DYB,FIJ1,FIJ2,UN
  REAL(SP) :: TXX,TYY,FXX,FYY,VISCOF,EXFLUX,TEMP,STPOINT
  REAL(SP) :: FACT,FM1
  INTEGER  :: I,I1,I2,IA,IB,J,J1,J2,K,JTMP,JJ,II
# if defined (TVD)
  REAL(SP), DIMENSION(1:M)             :: dstA,dstB
  REAL(SP), DIMENSION(1:NCV)           :: XUA,XUB,YUA,YUB
  REAL(SP) :: LX,LY,minA,minB,S_upstreamA,S_upstreamB,A_limiter
  REAL(SP) :: B_limiter,smoothrA,smoothrB
# else
  REAL(SP) :: S1MIN, S1MAX, S2MIN, S2MAX
# endif


!!$#  if defined (SPHERICAL)
!!$  REAL(DP) :: TY,TXPI,TYPI
!!$  REAL(DP) :: XTMP1,XTMP
!!$  REAL(DP) :: X1_DP,Y1_DP,X2_DP,Y2_DP,XII,YII
!!$  REAL(DP) :: X11_TMP,Y11_TMP,X33_TMP,Y33_TMP
!!$  REAL(DP) :: VX1_TMP,VY1_TMP,VX2_TMP,VY2_TMP
!!$  REAL(DP) :: TXPI_TMP,TYPI_TMP
!!$#  endif

#  if defined (SEMI_IMPLICIT)
   REAL(SP) :: UN1
   REAL(SP), DIMENSION(3*(NT),KBM1) :: UVN1
   REAL(SP), DIMENSION(3*(NT),KBM1) :: DTIJ1
#  endif
#  if defined (MPDATA)
  REAL(SP) :: SMIN,SMAX,XXXX,Vadj
  REAL(SP), DIMENSION(0:MT,KB)     :: S1_S, S1m, S1_FRESHm    !! temporary salinity in modified upwind
  REAL(SP), DIMENSION(0:MT,KB)     :: S1_SF   !! temporary salinity in modified upwind
  REAL(SP), DIMENSION(0:MT,KB)     :: WWWS     
  REAL(SP), DIMENSION(0:MT,KB)     :: WWWSF   
  REAL(SP), DIMENSION(0:MT)        :: DTWWWS  
  REAL(SP), DIMENSION(0:MT,KB)     :: ZZZFLUX !! temporary total flux in corrected part
  REAL(SP), DIMENSION(0:MT,KB)     :: BETA    !! temporary beta coefficient in corrected part
  REAL(SP), DIMENSION(0:MT,KB)     :: BETAIN  !! temporary beta coefficient in corrected part
  REAL(SP), DIMENSION(0:MT,KB)     :: BETAOUT !! temporary beta coefficient in corrected part
  REAL(SP), DIMENSION(0:MT,KB)     :: S1_FRESH    !! for source term which also bring mass volume
  REAL(SP), DIMENSION(0:MT,KB)     :: OFFS     !! Offset to run the simulation in single precision
  INTEGER ITERA, NTERA
#  endif

# if defined (NH) || defined (LIMITER_VER_ADV)
   REAL(SP)  CONV_S(1:KB), DISS_S(1:KB)
   REAL(SP)  SL_H(0:KB), S_TEMP(0:KB)
   REAL(SP)  SL_U, SL_F
# endif

# if defined (THIN_DAM)
  INTEGER :: NX
# endif


  IF(DBG_SET(DBG_SBR)) WRITE(IPT,*)"Start: adv_s"
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

# if defined (MPDATA)
  ! Adding offset to force FVCOM not to crash when salinities approach zero
      OFFS = 1.0_SP
      S1  = S1+OFFS
      SMEAN1 = SMEAN1 + OFFS
# endif

  !
  !--Initialize Fluxes-----------------------------------------------------------!
  !
  XFLUX     = 0.0_SP
  XFLUX_ADV = 0.0_SP

  !
  !--Loop Over Control Volume Sub-Edges And Calculate Normal Velocity------------!
  !
!!#  if !defined (WET_DRY)
  DO I=1,NCV
     I1=NTRG(I)
     !     DTIJ(I)=DT1(I1)
     DO K=1,KBM1
       DTIJ(I,K) = DT1(I1)*DZ1(I1,K)
       ! USE U,V
       UVN(I,K)  = V(I1,K)*DLTXE(I) - U(I1,K)*DLTYE(I)
#      if defined (SEMI_IMPLICIT)
       DTIJ1(I,K) = D1(I1)*DZ1(I1,K)
       UVN1(I,K) = VF(I1,K)*DLTXE(I) - UF(I1,K)*DLTYE(I)
#      endif
     END DO
  END DO
!!#  else
!!  DO I=1,NCV
!!     I1=NTRG(I)
!!     !     DTIJ(I)=DT1(I1)
!!     DO K=1,KBM1
!!       DTIJ(I,K) = DT1(I1)*DZ1(I1,K)
!!       ! USE US,VS
!!       UVN(I,K) = VS(I1,K)*DLTXE(I) - US(I1,K)*DLTYE(I)
!!#      if defined (SEMI_IMPLICIT)
!!       DTIJ1(I,K) = D1(I1)*DZ1(I1,K)
!!       UVN1(I,K) = VF(I1,K)*DLTXE(I) - UF(I1,K)*DLTYE(I)
!!#      endif
!!     END DO
!!  END DO
!!#  endif

  !
  !--Calculate the Advection and Horizontal Diffusion Terms----------------------!
  !

  DO K=1,KBM1
     PSPX  = 0.0_SP 
     PSPY  = 0.0_SP 
     PSPXD = 0.0_SP 
     PSPYD = 0.0_SP
     DO I=1,M
        DO J=1,NTSN(I)-1
           I1=NBSN(I,J)
           I2=NBSN(I,J+1)

#    if defined (WET_DRY)
!J. Ge for tracer advection
           IF(BACKWARD_ADVECTION.eqv..FALSE.)THEN

            IF(ISWETN(I1) == 0 .AND. ISWETN(I2) == 1)THEN
             FFD=0.5_SP*(S1(I,K)+S1(I2,K)-SMEAN1(I,K)-SMEAN1(I2,K))
             FF1=0.5_SP*(S1(I,K)+S1(I2,K))
            ELSE IF(ISWETN(I1) == 1 .AND. ISWETN(I2) == 0)THEN
             FFD=0.5_SP*(S1(I1,K)+S1(I,K)-SMEAN1(I1,K)-SMEAN1(I,K))
             FF1=0.5_SP*(S1(I1,K)+S1(I,K))
            ELSE IF(ISWETN(I1) == 0 .AND. ISWETN(I2) == 0)THEN
             FFD=0.5_SP*(S1(I,K)+S1(I,K)-SMEAN1(I,K)-SMEAN1(I,K))
             FF1=0.5_SP*(S1(I,K)+S1(I,K))
            ELSE
             FFD=0.5_SP*(S1(I1,K)+S1(I2,K)-SMEAN1(I1,K)-SMEAN1(I2,K))
             FF1=0.5_SP*(S1(I1,K)+S1(I2,K))
            END IF

           ELSE
         
            IF(BACKWARD_STEP==1)THEN
              IF(ISWETN(I1) == 0 .AND. ISWETN(I2) == 1)THEN
               FFD=0.5_SP*((S0(I,K)+S1(I,K))*0.5+(S0(I2,K)+S1(I2,K))*0.5-SMEAN1(I,K)-SMEAN1(I2,K))
               FF1=0.5_SP*((S0(I,K)+S1(I,K))*0.5+(S0(I2,K)+S1(I2,K))*0.5)
              ELSE IF(ISWETN(I1) == 1 .AND. ISWETN(I2) == 0)THEN
               FFD=0.5_SP*((S0(I1,K)+S1(I1,K))*0.5+(S0(I,K)+S1(I,K))*0.5-SMEAN1(I1,K)-SMEAN1(I,K))
               FF1=0.5_SP*((S0(I1,K)+S1(I1,K))*0.5+(S0(I,K)+S1(I,K))*0.5)
	      ELSE IF(ISWETN(I1) == 0 .AND. ISWETN(I2) == 0)THEN
               FFD=0.5_SP*((S0(I,K)+S1(I,K))*0.5+(S0(I,K)+S1(I,K))*0.5-SMEAN1(I,K)-SMEAN1(I,K))
               FF1=0.5_SP*((S0(I,K)+S1(I,K))*0.5+(S0(I,K)+S1(I,K))*0.5)
	      ELSE
               FFD=0.5_SP*((S0(I1,K)+S1(I1,K))*0.5+(S0(I2,K)+S1(I2,K))*0.5-SMEAN1(I1,K)-SMEAN1(I2,K))
               FF1=0.5_SP*((S0(I1,K)+S1(I1,K))*0.5+(S0(I2,K)+S1(I2,K))*0.5)
	      END IF 
            ELSEIF(BACKWARD_STEP==2)THEN
              IF(ISWETN(I1) == 0 .AND. ISWETN(I2) == 1)THEN
               FFD=0.5_SP*((S2(I,K)+S0(I,K)+S1(I,K))/3.0_SP+(S2(I2,K)+S0(I2,K)+S1(I2,K))/3.0_SP-SMEAN1(I,K)-SMEAN1(I2,K))
               FF1=0.5_SP*((S2(I,K)+S0(I,K)+S1(I,K))/3.0_SP+(S2(I2,K)+S0(I2,K)+S1(I2,K))/3.0_SP)
              ELSE IF(ISWETN(I1) == 1 .AND. ISWETN(I2) == 0)THEN
               FFD=0.5_SP*((S2(I1,K)+S0(I1,K)+S1(I1,K))/3.0_SP+(S2(I,K)+S0(I,K)+S1(I,K))/3.0_SP-SMEAN1(I1,K)-SMEAN1(I,K))
               FF1=0.5_SP*((S2(I1,K)+S0(I1,K)+S1(I1,K))/3.0_SP+(S2(I,K)+S0(I,K)+S1(I,K))/3.0_SP)
	      ELSE IF(ISWETN(I1) == 0 .AND. ISWETN(I2) == 0)THEN
               FFD=0.5_SP*((S2(I,K)+S0(I,K)+S1(I,K))/3.0_SP+(S2(I,K)+S0(I,K)+S1(I,K))/3.0_SP-SMEAN1(I,K)-SMEAN1(I,K))
               FF1=0.5_SP*((S2(I,K)+S0(I,K)+S1(I,K))/3.0_SP+(S2(I,K)+S0(I,K)+S1(I,K))/3.0_SP)
	      ELSE
               FFD=0.5_SP*((S2(I1,K)+S0(I1,K)+S1(I1,K))/3.0_SP+(S2(I2,K)+S0(I2,K)+S1(I2,K))/3.0_SP-SMEAN1(I1,K)-SMEAN1(I2,K))
               FF1=0.5_SP*((S2(I1,K)+S0(I1,K)+S1(I1,K))/3.0_SP+(S2(I2,K)+S0(I2,K)+S1(I2,K))/3.0_SP)
	      END IF 
            ENDIF

           END IF
!J. Ge for tracer advection
#    else
!J. Ge for tracer advection	 
           IF(BACKWARD_ADVECTION.eqv..FALSE.)THEN
             FFD=0.5_SP*(S1(I1,K)+S1(I2,K)-SMEAN1(I1,K)-SMEAN1(I2,K))
             FF1=0.5_SP*(S1(I1,K)+S1(I2,K))
           ELSE
             IF(BACKWARD_STEP==1)THEN
               FFD=0.5_SP*((S0(I1,K)+S1(I1,K))*0.5+(S0(I2,K)+S1(I2,K))*0.5-SMEAN1(I1,K)-SMEAN1(I2,K))
               FF1=0.5_SP*((S0(I1,K)+S1(I1,K))*0.5+(S0(I2,K)+S1(I2,K))*0.5)
             ELSEIF(BACKWARD_STEP==2)THEN
               FFD=0.5_SP*((S2(I1,K)+S0(I1,K)+S1(I1,K))/3.0_SP+(S2(I2,K)+S0(I2,K)+S1(I2,K))/3.0_SP-SMEAN1(I1,K)-SMEAN1(I2,K))
               FF1=0.5_SP*((S2(I1,K)+S0(I1,K)+S1(I1,K))/3.0_SP+(S2(I2,K)+S0(I2,K)+S1(I2,K))/3.0_SP)
             ENDIF
           ENDIF
!J. Ge for tracer advection
#    endif	 
	 

!!$#        if defined (SPHERICAL)
!!$           XTMP  = VX(I2)*TPI-VX(I1)*TPI
!!$           XTMP1 = VX(I2)-VX(I1)
!!$           IF(XTMP1 >  180.0_SP)THEN
!!$              XTMP = -360.0_SP*TPI+XTMP
!!$           ELSE IF(XTMP1 < -180.0_SP)THEN
!!$              XTMP =  360.0_SP*TPI+XTMP
!!$           END IF
!!$           TXPI=XTMP*COS(DEG2RAD*VY(I))
!!$           TYPI=(VY(I1)-VY(I2))*TPI
!!$           ! ERROR HERE
!!$!#    if defined (NORTHPOLE)
!!$           IF(NODE_NORTHAREA(I) == 1)THEN
!!$              VX1_TMP = REARTH * COS(VY(I1)*DEG2RAD) * COS(VX(I1)*DEG2RAD) &
!!$                   * 2._SP /(1._SP+sin(VY(I1)*DEG2RAD))
!!$              VY1_TMP = REARTH * COS(VY(I1)*DEG2RAD) * SIN(VX(I1)*DEG2RAD) &
!!$                   * 2._SP /(1._SP+sin(VY(I1)*DEG2RAD))
!!$
!!$              VX2_TMP = REARTH * COS(VY(I2)*DEG2RAD) * COS(VX(I2)*DEG2RAD) &
!!$                   * 2._SP /(1._SP+sin(VY(I2)*DEG2RAD))
!!$              VY2_TMP = REARTH * COS(VY(I2)*DEG2RAD) * SIN(VX(I2)*DEG2RAD) &
!!$                   * 2._SP /(1._SP+sin(VY(I2)*DEG2RAD))
!!$
!!$              TXPI = (VX2_TMP-VX1_TMP)/(2._SP /(1._SP+sin(VY(I)*DEG2RAD)))
!!$              TYPI = (VY1_TMP-VY2_TMP)/(2._SP /(1._SP+sin(VY(I)*DEG2RAD)))
!!$              IF(I /= NODE_NORTHPOLE)THEN
!!$                 TXPI_TMP = TYPI*COS(VX(I)*DEG2RAD)-TXPI*SIN(VX(I)*DEG2RAD)
!!$                 TYPI_TMP = TXPI*COS(VX(I)*DEG2RAD)+TYPI*SIN(VX(I)*DEG2RAD)
!!$                 TYPI_TMP = -TYPI_TMP
!!$
!!$                 TXPI = TXPI_TMP
!!$                 TYPI = TYPI_TMP
!!$              END IF
!!$           END IF
!!$# endif
!!$           ! END ERROR
!!$#        else
!!$           PSPX(I)=PSPX(I)+FF1*(VY(I1)-VY(I2))
!!$           PSPY(I)=PSPY(I)+FF1*(VX(I2)-VX(I1))
!!$           PSPXD(I)=PSPXD(I)+FFD*(VY(I1)-VY(I2))
!!$           PSPYD(I)=PSPYD(I)+FFD*(VX(I2)-VX(I1))
!!$#        endif

           
           PSPX(I)=PSPX(I)+FF1*DLTYTRIE(i,j)
           PSPY(I)=PSPY(I)+FF1*DLTXTRIE(i,j)
           PSPXD(I)=PSPXD(I)+FFD*DLTYTRIE(i,j)
           PSPYD(I)=PSPYD(I)+FFD*DLTXTRIE(i,j)
           

        END DO
! gather all neighboring control volumes connecting at dam node 
# if defined (THIN_DAM)
       IF(IS_DAM(I)==1.AND.K<=KDAM(I))THEN
         DO NX=1,N_DAM_MATCH(I,1)
           DO J=1,NTSN(N_DAM_MATCH(I,NX+1))-1
             I1=NBSN(N_DAM_MATCH(I,NX+1),J)
             I2=NBSN(N_DAM_MATCH(I,NX+1),J+1)
!J. Ge for tracer advection
             IF(BACKWARD_ADVECTION.eqv..FALSE.)THEN
               FFD=0.5_SP*(S1(I1,K)+S1(I2,K)-SMEAN1(I1,K)-SMEAN1(I2,K))
               FF1=0.5_SP*(S1(I1,K)+S1(I2,K))
             ELSE
               IF(BACKWARD_STEP==1)THEN
                FFD=0.5_SP*((S0(I1,K)+S1(I1,K))*0.5+(S0(I2,K)+S1(I2,K))*0.5-SMEAN1(I1,K)-SMEAN1(I2,K))
                FF1=0.5_SP*((S0(I1,K)+S1(I1,K))*0.5+(S0(I2,K)+S1(I2,K))*0.5)
               ELSEIF(BACKWARD_STEP==2)THEN
                FFD=0.5_SP*((S2(I1,K)+S0(I1,K)+S1(I1,K))/3.0_SP+(S2(I2,K)+S0(I2,K)+S1(I2,K))/3.0_SP-SMEAN1(I1,K)-SMEAN1(I2,K))
                FF1=0.5_SP*((S2(I1,K)+S0(I1,K)+S1(I1,K))/3.0_SP+(S2(I2,K)+S0(I2,K)+S1(I2,K))/3.0_SP)
               ENDIF
             ENDIF
!J. Ge for tracer advection
             PSPX(I)=PSPX(I)+FF1*DLTYTRIE(N_DAM_MATCH(I,NX+1),j)
             PSPY(I)=PSPY(I)+FF1*DLTXTRIE(N_DAM_MATCH(I,NX+1),j)
             PSPXD(I)=PSPXD(I)+FFD*DLTYTRIE(N_DAM_MATCH(I,NX+1),j)
             PSPYD(I)=PSPYD(I)+FFD*DLTXTRIE(N_DAM_MATCH(I,NX+1),j)
           END DO
         END DO
       END IF
# endif

# if !defined (THIN_DAM)
        PSPX(I)=PSPX(I)/ART2(I)
        PSPY(I)=PSPY(I)/ART2(I)
        PSPXD(I)=PSPXD(I)/ART2(I)
        PSPYD(I)=PSPYD(I)/ART2(I)
# else
        IF(IS_DAM(I)==1.AND.K<=KDAM(I))THEN
          PSPX(I)=PSPX(I)/(ART2(I)+SUM(ART2(N_DAM_MATCH(I,2:1+N_DAM_MATCH(I,1)))))
          PSPY(I)=PSPY(I)/(ART2(I)+SUM(ART2(N_DAM_MATCH(I,2:1+N_DAM_MATCH(I,1)))))
          PSPXD(I)=PSPXD(I)/(ART2(I)+SUM(ART2(N_DAM_MATCH(I,2:1+N_DAM_MATCH(I,1)))))
          PSPYD(I)=PSPYD(I)/(ART2(I)+SUM(ART2(N_DAM_MATCH(I,2:1+N_DAM_MATCH(I,1)))))
        ELSE
          PSPX(I)=PSPX(I)/ART2(I)
          PSPY(I)=PSPY(I)/ART2(I)
          PSPXD(I)=PSPXD(I)/ART2(I)
          PSPYD(I)=PSPYD(I)/ART2(I)
        END IF
        
# endif

     END DO

     IF(K == KBM1)THEN
        DO I=1,M
           PFPXB(I) = PSPX(I)
           PFPYB(I) = PSPY(I)
        END DO
     END IF

     DO I = 1,M
        VISCOFF(I)=VISCOFH(I,K)
     END DO

     IF(K == KBM1) THEN
        AH_BOTTOM(1:M) = (FACT*VISCOFF(1:M) + FM1)*NN_HVC(1:M)
     END IF


     DO I=1,NCV_I
        IA=NIEC(I,1)
        IB=NIEC(I,2)


!!$        XI=0.5_SP*(XIJE(I,1)+XIJE(I,2))
!!$        YI=0.5_SP*(YIJE(I,1)+YIJE(I,2))
!!$#      if defined (SPHERICAL)
!!$        X1_DP=XIJE(I,1)
!!$        Y1_DP=YIJE(I,1)
!!$        X2_DP=XIJE(I,2)
!!$        Y2_DP=YIJE(I,2)
!!$        CALL ARCC(X2_DP,Y2_DP,X1_DP,Y1_DP,XII,YII)
!!$        XI=XII		
!!$        XTMP  = XI*TPI-VX(IA)*TPI
!!$        XTMP1 = XI-VX(IA)
!!$        IF(XTMP1 >  180.0_SP)THEN
!!$           XTMP = -360.0_SP*TPI+XTMP
!!$        ELSE IF(XTMP1 < -180.0_SP)THEN
!!$           XTMP =  360.0_SP*TPI+XTMP
!!$        END IF
!!$
!!$        DXA=XTMP*COS(DEG2RAD*VY(IA))    
!!$        DYA=(YI-VY(IA))*TPI
!!$        XTMP  = XI*TPI-VX(IB)*TPI
!!$        XTMP1 = XI-VX(IB)
!!$        IF(XTMP1 >  180.0_SP)THEN
!!$           XTMP = -360.0_SP*TPI+XTMP
!!$        ELSE IF(XTMP1 < -180.0_SP)THEN
!!$           XTMP =  360.0_SP*TPI+XTMP
!!$        END IF
!!$
!!$        DXB=XTMP*COS(DEG2RAD*VY(IB)) 
!!$        DYB=(YI-VY(IB))*TPI
!!$#      else
!!$        DXA=XI-VX(IA)
!!$        DYA=YI-VY(IA)
!!$        DXB=XI-VX(IB)
!!$        DYB=YI-VY(IB)
!!$#      endif

# if !defined (TVD)
!J. Ge for tracer advection
        IF(BACKWARD_ADVECTION.eqv..FALSE.)THEN
          FIJ1=S1(IA,K)+DLTXNCVE(I,1)*PSPX(IA)+DLTYNCVE(I,1)*PSPY(IA)
          FIJ2=S1(IB,K)+DLTXNCVE(I,2)*PSPX(IB)+DLTYNCVE(I,2)*PSPY(IB)
        ELSE
          IF(BACKWARD_STEP==1)THEN
            FIJ1=(S0(IA,K)+S1(IA,K))*0.5+DLTXNCVE(I,1)*PSPX(IA)+DLTYNCVE(I,1)*PSPY(IA)
            FIJ2=(S0(IB,K)+S1(IB,K))*0.5+DLTXNCVE(I,2)*PSPX(IB)+DLTYNCVE(I,2)*PSPY(IB)
          ELSEIF(BACKWARD_STEP==2)THEN
            FIJ1=(S2(IA,K)+S0(IA,K)+S1(IA,K))/3.0_SP+DLTXNCVE(I,1)*PSPX(IA)+DLTYNCVE(I,1)*PSPY(IA)
            FIJ2=(S2(IB,K)+S0(IB,K)+S1(IB,K))/3.0_SP+DLTXNCVE(I,2)*PSPX(IB)+DLTYNCVE(I,2)*PSPY(IB)
          ENDIF
        ENDIF
!J. Ge for tracer advection
        S1MIN=MINVAL(S1(NBSN(IA,1:NTSN(IA)-1),K))
        S1MIN=MIN(S1MIN, S1(IA,K))
        S1MAX=MAXVAL(S1(NBSN(IA,1:NTSN(IA)-1),K))
        S1MAX=MAX(S1MAX, S1(IA,K))
        S2MIN=MINVAL(S1(NBSN(IB,1:NTSN(IB)-1),K))
        S2MIN=MIN(S2MIN, S1(IB,K))
        S2MAX=MAXVAL(S1(NBSN(IB,1:NTSN(IB)-1),K))
        S2MAX=MAX(S2MAX, S1(IB,K))
        IF(FIJ1 < S1MIN) FIJ1=S1MIN
        IF(FIJ1 > S1MAX) FIJ1=S1MAX
        IF(FIJ2 < S2MIN) FIJ2=S2MIN
        IF(FIJ2 > S2MAX) FIJ2=S2MAX
# else
! ------------------------------------------------------------------------------------------ !
!                                     Drawing the TVD scheme
! ------------------------------------------------------------------------------------------ !
!     If A is upstream of B
       S_upstreamA=S1(Anear_node(I),K)+PSPX(Anear_node(I))*XUAdist(I)+PSPY(Anear_node(I))*YUAdist(I)
!     If B is upstream of A
       S_upstreamB=S1(Bnear_node(I),K)+PSPX(Bnear_node(I))*XUBdist(I)+PSPY(Bnear_node(I))*YUBdist(I)                

!     Smoothness-parameter
!       IF (S1(IA,K).EQ.S1(IB,K)) THEN
!          smoothrA    = 0.0_SP
!          smoothrB    = 0.0_SP
!         ELSE
!     S_upstreamA if far upstream, and S1(IA,K) is upstream of S1(IB,K)
!          smoothrA    = (S1(IA,K)-S_upstreamA)/(S1(IB,K)-S1(IA,K))
!          smoothrB    = (S1(IB,K)-S_upstreamB)/(S1(IA,K)-S1(IB,K))
!       END IF

       IF (S1(IA,K).EQ.S1(IB,K)) THEN
          smoothrA    = 0.0_SP
          smoothrB    = 0.0_SP
         ELSE IF (S1(IA,K).LT.1.E-10.AND.S1(IB,K).LT.1.E-10) THEN
          smoothrA    = 0.0_SP
          smoothrB    = 0.0_SP
         ELSE
!         S1_upstreamA is the "far upstream node", S1(IA,K) is upstream of S1(IB,K)
          smoothrA    = (S1(IA,K)-S_upstreamA)/(S1(IB,K)-S1(IA,K))
          smoothrB    = (S1(IB,K)-S_upstreamB)/(S1(IA,K)-S1(IB,K))
       END IF

!     Flux-limiter
!     (superbee) - Feel free to use other limiters!
       A_limiter = MAX(0.0_SP, MIN(1.0_SP, 2.0_SP*smoothrA), MIN(2.0_SP, smoothrA))
       B_limiter = MAX(0.0_SP, MIN(1.0_SP, 2.0_SP*smoothrB), MIN(2.0_SP, smoothrB))

!     Estimating the salinity at the cellwall in-betweem IA and IB
!       S12_A     = S1(IA,K)+0.5_SP*A_LIMITER*(S1(IB,K)-S1(IA,K))
!       S12_B     = S1(IB,K)+0.5_SP*B_LIMITER*(S1(IA,K)-S1(IB,K))
       FIJ1     = S1(IA,K)+0.5_SP*A_LIMITER*(S1(IB,K)-S1(IA,K))
       FIJ2     = S1(IB,K)+0.5_SP*B_LIMITER*(S1(IA,K)-S1(IB,K))

# endif
        UN=UVN(I,K)
#      if defined (SEMI_IMPLICIT)
       UN1=UVN1(I,K)
#      endif

        !VISCOF=HORCON*(FACT*(VISCOFF(IA)+VISCOFF(IB))*0.5_SP + FM1)
        
        ! David moved HPRNU and added HVC
        VISCOF=(FACT*0.5_SP*(VISCOFF(IA)*NN_HVC(IA)+VISCOFF(IB)*NN_HVC(IB)) + FM1*0.5_SP*(NN_HVC(IA)+NN_HVC(IB)))
        
        TXX=0.5_SP*(PSPXD(IA)+PSPXD(IB))*VISCOF
        TYY=0.5_SP*(PSPYD(IA)+PSPYD(IB))*VISCOF

        FXX=-DTIJ(I,K)*TXX*DLTYE(I)
        FYY= DTIJ(I,K)*TYY*DLTXE(I)

#       if !defined (SEMI_IMPLICIT)
        EXFLUX=-UN*DTIJ(I,K)* &
             ((1.0_SP+SIGN(1.0_SP,UN))*FIJ2+(1.0_SP-SIGN(1.0_SP,UN))*FIJ1)*0.5_SP+FXX+FYY
#       else
        EXFLUX=-UN*DTIJ(I,K)* &
           ((1.0_SP+SIGN(1.0_SP,UN))*FIJ2+(1.0_SP-SIGN(1.0_SP,UN))*FIJ1)*0.5_SP
        EXFLUX=(1.0_SP-IFCETA)*EXFLUX+IFCETA*(-UN1*DTIJ1(I,K)*((1.0_SP+SIGN(1.0_SP,UN1))*FIJ2+(1.0_SP-SIGN(1.0_SP,UN1))*FIJ1)*0.5_SP)+FXX+FYY
#       endif

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
#    if !defined (SEMI_IMPLICIT)
       CALL ADV_S_XY(XFLUX,XFLUX_ADV,PSPX,PSPY,PSPXD,PSPYD,VISCOFF,K,0.0_SP)
#    else
       CALL ADV_S_XY(XFLUX,XFLUX_ADV,PSPX,PSPY,PSPXD,PSPYD,VISCOFF,K,IFCETA)
#    endif
#    endif  

  END DO !!SIGMA LOOP

  !
  !-Accumulate Fluxes at Boundary Nodes
  !
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

  !--Set Boundary Conditions-For Fresh Water Flux--------------------------------!
  !
# if defined (MPDATA)

  !   S. HU
  !   Using smolarkiewicz, P. K; A fully multidimensional positive definite advection
  !   transport algorithm with small implicit diffusion, Journal of Computational
  !   Physics, 54, 325-362, 1984
  !-----------------------------------------------------------------

  !-----combine all the horizontal flux first-----------------------------------

  !-------fresh water part--------------

  S1_FRESH=S1


  IF(RIVER_TS_SETTING == 'calculated') THEN
     IF(RIVER_INFLOW_LOCATION == 'node') THEN
        IF(NUMQBC > 0) THEN
           DO J=1,NUMQBC
              JJ=INODEQ(J)
              STPOINT=SDIS(J)+OFFS(1,1)
              DO K=1,KBM1
                 S1_FRESH(JJ,K)=SDIS(J)+OFFS(1,1)
                 XFLUX(JJ,K)=XFLUX(JJ,K) - QDIS(J)*VQDIST(J,K)*STPOINT !/DZ(JJ,K)
              END DO
           END DO
        END IF
     ELSE IF(RIVER_INFLOW_LOCATION == 'edge') THEN
        IF(NUMQBC > 0) THEN
           DO J=1,NUMQBC
              J1=N_ICELLQ(J,1)
              J2=N_ICELLQ(J,2)
              STPOINT=SDIS(J)+OFFS(1,1)
              DO K=1,KBM1
                 S1_FRESH(J1,K)=SDIS(J)+OFFS(1,1)
                 S1_FRESH(J2,K)=SDIS(J)+OFFS(1,1)
                 XFLUX(J1,K)=XFLUX(J1,K)-  &
                      QDIS(J)*RDISQ(J,1)*VQDIST(J,K)*STPOINT    !/DZ1(J1,K)
                 XFLUX(J2,K)=XFLUX(J2,K)-  &
                      QDIS(J)*RDISQ(J,2)*VQDIST(J,K)*STPOINT    !/DZ1(J2,K)
              END DO
           END DO
        END IF
     END IF
  END IF
  !

  ! -------------------------------------------------------------------
  !        The horizontal term of advection is neglected here
  ! -------------------------------------------------------------------
  DO K=1,KBM1
     DO I=1,M
        IF(ISONB(I) == 2) THEN
           XFLUX(I,K)=0._SP
        ENDIF
     END DO
  END DO

  ! Initialize variables of MPDATA
  S1_S=0._SP
  S1_SF=0._SP
  WWWS=0._SP
  WWWSF=0._SP
  DTWWWS=0._SP
  ZZZFLUX=0._SP
  BETA=0._SP
  BETAIN=0._SP
  BETAOUT=0._SP

! --------------------------------------------------------------------
!              Accounting for precipitation
! --------------------------------------------------------------------
! Added by Akvaplan, October 2018.
!
! Precipitation influences the surface salinity upon impact with the 
! ocean surface (total impact integrated as a flux over a timestep). 
! We therefore adjust the thresholds (SMIN/SMAX) at the nodes to 
! make sure that MPDATA allows changes in salinity due to precipitatio
! and evaporation.
! --------------------------------------------------------------------
  S1m = S1
  S1_FRESHm = S1_FRESH

! Adjusting the tresholds to include precipitation at cell I
   DO I=1,M
      S1m(I,1)       = (S1(I,1)-OFFS(1,1))*(D(I)*DZ(I,1)- &
                       (QPREC(I)+QEVAP(I))*DTI*ROFVROS)/(D(I)*DZ(I,1))
      S1_FRESHm(I,1) = (S1_FRESH(I,1)-OFFS(1,1))*(D(I)*DZ(I,1)- &
                       (QPREC(I)+QEVAP(I))*DTI*ROFVROS)/(D(I)*DZ(I,1))
   END DO

! Including the offset
  S1_FRESHm(:,1) = S1_FRESHm(:,1)+OFFS(:,1)
  S1m(:,1)       = S1m(:,1)+OFFS(:,1)

! Including the offset-flux of rainwater (This is zero in absence of
! the offset)
  DO I=1,M
     XFLUX(I,1)=XFLUX(I,1)-(QPREC(I)+QEVAP(I))*ART1(I)*OFFS(1,1)*ROFVROS
  END DO


! -----------------------------------------------------------------
!                  Back to the original MPDATA
! -----------------------------------------------------------------

  !!   first loop for vertical upwind
  !!   flux including horizontal and vertical upwind
  DO K=1,KBM1
     DO I=1,M
#    if defined (WET_DRY)
        IF(ISWETN(I)*ISWETNT(I) == 1) THEN
#    endif
           IF(K == 1) THEN
              TEMP = -(WTS(I,K+1)-ABS(WTS(I,K+1)))*S1(I,K)   &
                   -(WTS(I,K+1)+ABS(WTS(I,K+1)))*S1(I,K+1) &
                   +(WTS(I,K)+ABS(WTS(I,K)))*S1(I,K)    
           ELSE IF(K == KBM1) THEN
              TEMP = +(WTS(I,K)-ABS(WTS(I,K)))*S1(I,K-1)     &
                   +(WTS(I,K)+ABS(WTS(I,K)))*S1(I,K)
           ELSE
              TEMP = -(WTS(I,K+1)-ABS(WTS(I,K+1)))*S1(I,K)   &
                   -(WTS(I,K+1)+ABS(WTS(I,K+1)))*S1(I,K+1) &
                   +(WTS(I,K)-ABS(WTS(I,K)))*S1(I,K-1)     &
                   +(WTS(I,K)+ABS(WTS(I,K)))*S1(I,K)
           END IF
           TEMP = 0.5_SP*TEMP 
           
         ! The minimum value is the one where the entire neighbouring cell is transported into cell I
         ! and the precipitation to cell I further diminish the salt concentration.
          
           IF(K == 1)THEN
              SMAX = MAXVAL(S1(NBSN(I,1:NTSN(I)),K))
              SMIN = MINVAL(S1(NBSN(I,1:NTSN(I)),K))
              SMAX = MAX(SMAX,S1(I,K+1),S1(I,K),S1_FRESH(I,K))
              SMIN = MIN(SMIN,S1(I,K+1),S1(I,K),S1_FRESH(I,K))
           ELSEIF(K == KBM1) THEN
              SMAX = MAXVAL(S1(NBSN(I,1:NTSN(I)),K))
              SMIN = MINVAL(S1(NBSN(I,1:NTSN(I)),K))
              SMAX = MAX(SMAX,S1(I,K-1),S1(I,K),S1_FRESH(I,K))
              SMIN = MIN(SMIN,S1(I,K-1),S1(I,K),S1_FRESH(I,K))
           ELSE
              SMAX = MAXVAL(S1(NBSN(I,1:NTSN(I)),K))
              SMIN = MINVAL(S1(NBSN(I,1:NTSN(I)),K))
              SMAX = MAX(SMAX,S1(I,K+1),S1(I,K-1),S1(I,K),S1_FRESH(I,K))
              SMIN = MIN(SMIN,S1(I,K+1),S1(I,K-1),S1(I,K),S1_FRESH(I,K))
           END IF

         ! Total (horizontal + vertical) flux to the cell
           ZZZFLUX(I,K) = TEMP*(DTI/DT(I))/DZ(I,K) + XFLUX(I,K)/ART1(I)*(DTI/DT(I))/DZ(I,K)

         ! Updated salinity without limiter minus current salinity
           XXXX = ZZZFLUX(I,K)*DT(I)/DTFA(I)+S1(I,K)-S1(I,K)*DT(I)/DTFA(I)
           
         ! For clarity: The line above can be rewritten;
         ! XXXX = S1(I,K) - (S1(I,K) - ZZZFLUX(I,K))*DT(I)/DTFA(I)

         ! --------------------------------------------------------------------------------
         ! An interpretation of the next lines:
         ! --------------------------------------------------------------------------------
         ! If the advected salinity is bigger or smaller than the largest/smallest salinities
         ! in- and around the node the previous timestep, then MPDATA force the solution not 
         ! to over/undershoot by multiplying it with BETA. This is a dubious procedure...
         ! --------------------------------------------------------------------------------
         IF((ABS(XXXX).LT.ABS(SMAX-S1(I,K))).AND.(XXXX.LT.0.0_SP)) THEN
            S1_SF(I,K) = S1(I,K)-XXXX
         ELSE IF((ABS(XXXX).LT.ABS(SMIN-S1(I,K))).AND.(XXXX.GT.0.0_SP)) THEN
            S1_SF(I,K) = S1(I,K)-XXXX
         ELSE IF(XXXX.EQ.0.0_SP) THEN
            S1_SF(I,K) = S1(I,K)-XXXX
         ELSE
             BETA(I,K)=0.5*(1.-SIGN(1._SP,XXXX)) * (SMAX-S1(I,K))/(ABS(XXXX)+1.E-10) &
                 +0.5*(1.-SIGN(1._SP,-XXXX)) * (S1(I,K)-SMIN)/(ABS(XXXX)+1.E-10)

         ! Why do they store BETA - that is just using memory for no reason...
             S1_SF(I,K)=S1(I,K)-MIN(1._SP,BETA(I,K))*XXXX
         END IF
#    if defined (WET_DRY)
        END IF
#    endif
     END DO
  END DO  !! SIGMA LOOP

  !----------------------------------------------------------------------------------------
  NTERA = 4
  DO ITERA=1,NTERA   !! Smolaricizw Loop 
     IF(ITERA == 1)THEN
        WWWSF  = WTS
        S1_S   = S1_SF
        DTWWWS = DT
     ELSE
        WWWSF  = WWWS
        S1_S   = S1_SF
        DTWWWS = DTFA
     END IF
     DO K=2,KBM1
        DO I=1,M
           TEMP=ABS(WWWSF(I,K))-DTI*(WWWSF(I,K))*(WWWSF(I,K))/DZ(I,K)/DTWWWS(I)
           WWWS(I,K)=TEMP*(S1_S(I,K-1)-S1_S(I,K))/(ABS(S1_S(I,K-1))+ABS(S1_S(I,K))+1.E-14)

           IF(TEMP < 0.0_SP .OR. S1_S(I,K) == 0.0_SP)THEN 
              WWWS(I,K)=0. 
           END IF
        END DO
     END DO
     DO I=1,M
        WWWS(I,1)=0.0_SP
     END DO

     DO I=1,M
        SMAX = MAXVAL(S1(NBSN(I,1:NTSN(I)),1))
        SMIN = MINVAL(S1(NBSN(I,1:NTSN(I)),1))
        SMAX = MAX(SMAX,S1(I,2),S1(I,1),S1_FRESH(I,1))
        SMIN = MIN(SMIN,S1(I,2),S1(I,1),S1_FRESH(I,1))

        TEMP=0.5*((WWWS(I,2)+ABS(WWWS(I,2)))*S1_S(I,2))*(DTI/DTFA(I))/DZ(I,1)
        BETAIN(I,1)=(SMAX-S1_S(I,1))/(TEMP+1.E-10)

        TEMP=0.5*((WWWS(I,1)+ABS(WWWS(I,1)))*S1_S(I,1)-        &
             (WWWS(I,2)-ABS(WWWS(I,2)))*S1_S(I,1))*(DTI/DTFA(I))/DZ(I,1)
        BETAOUT(I,1)=(S1_S(I,1)-SMIN)/(TEMP+1.E-10)

        WWWSF(I,1)=0.5*MIN(1.,BETAOUT(I,1))*(WWWS(I,1)+ABS(WWWS(I,1))) + &
             0.5*MIN(1.,BETAIN(I,1))*(WWWS(I,1)-ABS(WWWS(I,1)))
     END DO

     DO K=2,KBM1-1
        DO I=1,M
           SMAX = MAXVAL(S1(NBSN(I,1:NTSN(I)),K))
           SMIN = MINVAL(S1(NBSN(I,1:NTSN(I)),K))
           SMAX = MAX(SMAX,S1(I,K+1),S1(I,K-1),S1(I,K),S1_FRESH(I,K))
           SMIN = MIN(SMIN,S1(I,K+1),S1(I,K-1),S1(I,K),S1_FRESH(I,K))

           TEMP=0.5*((WWWS(I,K+1)+ABS(WWWS(I,K+1)))*S1_S(I,K+1)-  &
                (WWWS(I,K)-ABS(WWWS(I,K)))*S1_S(I,K-1))*(DTI/DTFA(I))/DZ(I,K)
           BETAIN(I,K)=(SMAX-S1_S(I,K))/(TEMP+1.E-10)

           TEMP=0.5*((WWWS(I,K)+ABS(WWWS(I,K)))*S1_S(I,K)-        &
                (WWWS(I,K+1)-ABS(WWWS(I,K+1)))*S1_S(I,K))*(DTI/DTFA(I))/DZ(I,K)
           BETAOUT(I,K)=(S1_S(I,K)-SMIN)/(TEMP+1.E-10)

           WWWSF(I,K)=0.5*MIN(1.,BETAIN(I,K-1),BETAOUT(I,K))*(WWWS(I,K)+ABS(WWWS(I,K))) + &
                0.5*MIN(1.,BETAIN(I,K),BETAOUT(I,K-1))*(WWWS(I,K)-ABS(WWWS(I,K)))
        END DO
     END DO


     K=KBM1
     DO I=1,M
        SMAX = MAXVAL(S1(NBSN(I,1:NTSN(I)),K))
        SMIN = MINVAL(S1(NBSN(I,1:NTSN(I)),K))
        SMAX = MAX(SMAX,S1(I,K-1),S1(I,K),S1_FRESH(I,K))
        SMIN = MIN(SMIN,S1(I,K-1),S1(I,K),S1_FRESH(I,K))

        TEMP=0.5*((WWWS(I,K+1)+ABS(WWWS(I,K+1)))*S1_S(I,K+1)-  &
             (WWWS(I,K)-ABS(WWWS(I,K)))*S1_S(I,K-1))*(DTI/DTFA(I))/DZ(I,K)
        BETAIN(I,K)=(SMAX-S1_S(I,K))/(TEMP+1.E-10)

        TEMP=0.5*((WWWS(I,K)+ABS(WWWS(I,K)))*S1_S(I,K)-        &
             (WWWS(I,K+1)-ABS(WWWS(I,K+1)))*S1_S(I,K))*(DTI/DTFA(I))/DZ(I,K)
        BETAOUT(I,K)=(S1_S(I,K)-SMIN)/(TEMP+1.E-10)

        WWWSF(I,K)=0.5*MIN(1.,BETAIN(I,K-1),BETAOUT(I,K))*(WWWS(I,K)+ABS(WWWS(I,K))) + &
             0.5*MIN(1.,BETAIN(I,K),BETAOUT(I,K-1))*(WWWS(I,K)-ABS(WWWS(I,K)))
     END DO


     WWWS=WWWSF 

     DO K=1,KBM1
        DO I=1,M
#      if defined (WET_DRY)
           IF(ISWETN(I)*ISWETNT(I) == 1) THEN
#      endif
              IF(K == 1) THEN
                 TEMP = -(WWWS(I,K+1)-ABS(WWWS(I,K+1)))*S1_S(I,K)   &
                      -(WWWS(I,K+1)+ABS(WWWS(I,K+1)))*S1_S(I,K+1) &
                      +(WWWS(I,K)+ABS(WWWS(I,K)))*S1_S(I,K)
              ELSE IF(K == KBM1) THEN
                 TEMP = +(WWWS(I,K)-ABS(WWWS(I,K)))*S1_S(I,K-1)     &
                      +(WWWS(I,K)+ABS(WWWS(I,K)))*S1_S(I,K)
              ELSE
                 TEMP = -(WWWS(I,K+1)-ABS(WWWS(I,K+1)))*S1_S(I,K)   &
                      -(WWWS(I,K+1)+ABS(WWWS(I,K+1)))*S1_S(I,K+1) &
                      +(WWWS(I,K)-ABS(WWWS(I,K)))*S1_S(I,K-1)     &
                      +(WWWS(I,K)+ABS(WWWS(I,K)))*S1_S(I,K)
              END IF
              TEMP = 0.5_SP*TEMP
              S1_SF(I,K)=(S1_S(I,K)-TEMP*(DTI/DTFA(I))/DZ(I,K)) 
#      if defined (WET_DRY)
           END IF
#      endif
        END DO
     END DO  !! SIGMA LOOP
  END DO  !! Smolarvizw Loop
  !--------------------------------------------------------------------------
  ! End of smolarkiewicz upwind loop
  !--------------------------------------------------------------------------

      S1 = S1 - OFFS
      S1_SF = S1_SF - OFFS
      SMEAN1 = SMEAN1 - OFFS
#  endif


# if ! defined(MPDATA)
# if defined (ONE_D_MODEL)
    XFLUX = 0.0_SP
# endif


  !--------------------------------------------------------------------
  !   The central difference scheme in vertical advection
  !--------------------------------------------------------------------
  DO I=1,M
#    if defined (WET_DRY)
     IF(ISWETN(I)*ISWETNT(I) == 1) THEN
#    endif

#       if defined (NH) || defined (LIMITER_VER_ADV)
        S_TEMP(0)  = -S1(I,1)
        S_TEMP(KB) = -S1(I,KBM1)
        SL_H(0)    = DZ(I,1)
        SL_H(KB)   = DZ(I,KBM1)
        DO K=1, KBM1
          S_TEMP(K) = S1(I,K)
          SL_H(K)   = DZ(I,K)
        ENDDO

        DO K=2, KBM1
          CONV_S(K) = WTS(I,K)*(S_TEMP(K)+S_TEMP(K-1))*0.5_SP
          SL_U = 2.0_SP*(S_TEMP(K)-S_TEMP(K+1))/(SL_H(K)+SL_H(K+1))
          SL_F = 2.0_SP*(S_TEMP(K-2)-S_TEMP(K-1))/(SL_H(K-2)+SL_H(K-1))
          DISS_S(K) = 0.5_SP*ABS(WTS(I,K))*(S_TEMP(K-1)-S_TEMP(K)-0.5_SP*LIMLED2(SL_U,SL_F,2.0_SP)*(SL_H(K-1)+SL_H(K)))
        ENDDO
        CONV_S(1)  = 0.0_SP
        DISS_S(1)  = 0.0_SP
        CONV_S(KB) = 0.0_SP
        DISS_S(KB) = 0.0_SP
#       endif

        DO K=1, KBM1

#          if defined (NH) || defined (LIMITER_VER_ADV)
           TEMP = CONV_S(K)-CONV_S(K+1)+DISS_S(K+1)-DISS_S(K)
#          endif

#          if !defined (NH) && !defined (LIMITER_VER_ADV)

!J. Ge for tracer advection
           IF(BACKWARD_ADVECTION.eqv..FALSE.)THEN
             IF(K == 1) THEN
              TEMP=-WTS(I,K+1)*(S1(I,K)*DZ(I,K+1)+S1(I,K+1)*DZ(I,K))/   &
                   (DZ(I,K)+DZ(I,K+1))
             ELSE IF(K == KBM1) THEN
              TEMP= WTS(I,K)*(S1(I,K)*DZ(I,K-1)+S1(I,K-1)*DZ(I,K))/(DZ(I,K)+DZ(I,K-1))
             ELSE
              TEMP= WTS(I,K)*(S1(I,K)*DZ(I,K-1)+S1(I,K-1)*DZ(I,K))/(DZ(I,K)+DZ(I,K-1))-&
                   WTS(I,K+1)*(S1(I,K)*DZ(I,K+1)+S1(I,K+1)*DZ(I,K))/(DZ(I,K)+DZ(I,K+1))
             END IF
           ELSE
             IF(BACKWARD_STEP==1)THEN
               IF(K == 1) THEN
                TEMP=-WTS(I,K+1)*((S0(I,K)+S1(I,K))*0.5*DZ(I,K+1)+(S0(I,K+1)+S1(I,K+1))*0.5*DZ(I,K))/   &
                   (DZ(I,K)+DZ(I,K+1))
               ELSE IF(K == KBM1) THEN
                TEMP= WTS(I,K)*((S0(I,K)+S1(I,K))*0.5*DZ(I,K-1)+(S0(I,K-1)+S1(I,K-1))*0.5*DZ(I,K))/(DZ(I,K)+DZ(I,K-1))
               ELSE
                 TEMP= WTS(I,K)*((S0(I,K)+S1(I,K))*0.5*DZ(I,K-1)+(S0(I,K-1)+S1(I,K-1))*0.5*DZ(I,K))/(DZ(I,K)+DZ(I,K-1))-&
                   WTS(I,K+1)*((S0(I,K)+S1(I,K))*0.5*DZ(I,K+1)+(S0(I,K+1)+S1(I,K+1))*0.5*DZ(I,K))/(DZ(I,K)+DZ(I,K+1))
               END IF
             ELSEIF(BACKWARD_STEP==2)THEN
               IF(K == 1) THEN
                TEMP=-WTS(I,K+1)*((S2(I,K)+S0(I,K)+S1(I,K))/3.0_SP*DZ(I,K+1)+(S2(I,K+1)+S0(I,K+1)+S1(I,K+1))/3.0_SP*DZ(I,K))/   &
                   (DZ(I,K)+DZ(I,K+1))
               ELSE IF(K == KBM1) THEN
                TEMP= WTS(I,K)*((S2(I,K)+S0(I,K)+S1(I,K))/3.0_SP*DZ(I,K-1)+(S2(I,K-1)+S0(I,K-1)+S1(I,K-1))/3.0_SP*DZ(I,K))/(DZ(I,K)+DZ(I,K-1))
               ELSE
                 TEMP= WTS(I,K)*((S2(I,K)+S0(I,K)+S1(I,K))/3.0_SP*DZ(I,K-1)+(S2(I,K-1)+S0(I,K-1)+S1(I,K-1))/3.0_SP*DZ(I,K))/(DZ(I,K)+DZ(I,K-1))-&
                   WTS(I,K+1)*((S2(I,K)+S0(I,K)+S1(I,K))/3.0_SP*DZ(I,K+1)+(S2(I,K+1)+S0(I,K+1)+S1(I,K+1))/3.0_SP*DZ(I,K))/(DZ(I,K)+DZ(I,K+1))
               END IF
             ENDIF
           ENDIF
!J. Ge for tracer advection

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
  END DO  !! SIGMA LOOP

  !
  !--Set Boundary Conditions-For Fresh Water Flux--------------------------------!
  !
  IF(RIVER_TS_SETTING == 'calculated') THEN
     IF(RIVER_INFLOW_LOCATION == 'node') THEN
        IF(NUMQBC > 0) THEN
           DO J=1,NUMQBC
              JJ=INODEQ(J)
              STPOINT=SDIS(J)
              DO K=1,KBM1
                 !             XFLUX(JJ,K)=XFLUX(JJ,K) - QDIS(J)*VQDIST(J,K)*STPOINT/DZ(JJ,K)
                 XFLUX(JJ,K)=XFLUX(JJ,K) - QDIS(J)*VQDIST(J,K)*STPOINT
              END DO
           END DO
        END IF
     ELSE IF(RIVER_INFLOW_LOCATION == 'edge') THEN
        IF(NUMQBC > 0) THEN
           DO J=1,NUMQBC
              J1=N_ICELLQ(J,1)
              J2=N_ICELLQ(J,2)
              STPOINT=SDIS(J)
              DO K=1,KBM1
                 XFLUX(J1,K)=XFLUX(J1,K)-QDIS(J)*RDISQ(J,1)*VQDIST(J,K)*STPOINT     !/DZ1(J1,K)
                 XFLUX(J2,K)=XFLUX(J2,K)-QDIS(J)*RDISQ(J,2)*VQDIST(J,K)*STPOINT     !/DZ1(J2,K)
              END DO
           END DO
        END IF
     END IF
  END IF
!---------------------------------------------------------------------
# endif


  ! APPLY GROUND WATER SALINITY FORCING
  IF(GROUNDWATER_ON .and. GROUNDWATER_SALT_ON)THEN
     DO I=1,M
        XFLUX(I,KBM1)=XFLUX(I,KBM1)-BFWDIS(I)*BFWSLT(I)
     END DO
  ELSEIF(GROUNDWATER_ON) THEN
     DO I=1,M
!        XFLUX(I,KBM1)=XFLUX(I,KBM1)-BFWDIS(I)*S1(I,KBM1)/DZ(I,KBM1)
!J. Ge for tracer advection
        IF(BACKWARD_ADVECTION.eqv..FALSE.)THEN
          XFLUX(I,KBM1)=XFLUX(I,KBM1)-BFWDIS(I)*S1(I,KBM1)
        ELSE
          IF(BACKWARD_STEP==1)THEN
            XFLUX(I,KBM1)=XFLUX(I,KBM1)-BFWDIS(I)*(S0(I,KBM1)+S1(I,KBM1))*0.5
          ELSEIF(BACKWARD_STEP==2)THEN
            XFLUX(I,KBM1)=XFLUX(I,KBM1)-BFWDIS(I)*(S2(I,KBM1)+S0(I,KBM1)+S1(I,KBM1))/3.0_SP
          ENDIF
        ENDIF
!J. Ge for tracer advection
     END DO
  END IF


  !--Update Salinity-------------------------------------------------------------!
  !
#  if defined (WET_DRY)
  DO I=1,M
     IF(ISWETN(I)*ISWETNT(I) == 1 )THEN
        DO K=1,KBM1
#    if !defined (MPDATA)     
           !       SF1(I,K)=(S1(I,K)-XFLUX(I,K)/ART1(I)*(DTI/DT(I)))*(DT(I)/DTFA(I))
#      if !defined (SEMI_IMPLICIT)

#      if !defined (THIN_DAM)
           SF1(I,K)=(S1(I,K)-XFLUX(I,K)/ART1(I)*(DTI/(DT(I)*DZ(I,K))))*(DT(I)/DTFA(I))
#      else
           IF(IS_DAM(I)==1.AND.K<=KDAM(I))THEN
              SF1(I,K)=(S1(I,K)-XFLUX(I,K)/(ART1(I)&
                   &+SUM(ART1(N_DAM_MATCH(I,2:1+N_DAM_MATCH(I,1)))))*(DTI/(DT(I)*DZ(I,K))))*(DT(I)/DTFA(I))
           ELSE
              SF1(I,K)=(S1(I,K)-XFLUX(I,K)/ART1(I)*(DTI/(DT(I)*DZ(I,K))))*(DT(I)/DTFA(I))
           END IF
#      endif

#      else
       SF1(I,K)=(TSC(I,K)-RK_TS(STAGE)*XFLUX(I,K)/ART1(I)*(DTI/(DT(I)*DZ(I,K))))*(DT(I)/DTFA(I))
#      endif
#    else
           SF1(I,K)=S1_SF(I,K)
#    endif              
        END DO

     ELSE ! IF NOT WET S STAYS THE SAME
        DO K=1,KBM1
#      if !defined (SEMI_IMPLICIT)
           SF1(I,K)=S1(I,K)
#      else
           SF1(I,K)=TSC(I,K)
#      endif
        END DO
     END IF
  END DO

# else

  DO I=1,M
     DO K=1,KBM1
#    if !defined (MPDATA)     
        !       SF1(I,K)=(S1(I,K)-XFLUX(I,K)/ART1(I)*(DTI/DT(I)))*(DT(I)/DTFA(I))
#      if !defined (SEMI_IMPLICIT)

#      if !defined (THIN_DAM)
        SF1(I,K)=(S1(I,K)-XFLUX(I,K)/ART1(I)*(DTI/(DT(I)*DZ(I,K))))*(DT(I)/DTFA(I))
#      else
        IF(IS_DAM(I)==1.AND.K<=KDAM(I))THEN
          SF1(I,K)=(S1(I,K)-XFLUX(I,K)/(ART1(I)&
               &+SUM(ART1(N_DAM_MATCH(I,2:1+N_DAM_MATCH(I,1)))))*(DTI/(DT(I)*DZ(I,K))))*(DT(I)/DTFA(I))
        ELSE
          SF1(I,K)=(S1(I,K)-XFLUX(I,K)/ART1(I)*(DTI/(DT(I)*DZ(I,K))))*(DT(I)/DTFA(I))
        END IF
#      endif

#      else
       SF1(I,K)=(TSC(I,K)-RK_TS(STAGE)*XFLUX(I,K)/ART1(I)*(DTI/(DT(I)*DZ(I,K))))*(DT(I)/DTFA(I))
#      endif
#    else
        SF1(I,K)=S1_SF(I,K)
#    endif              
     END DO

  END DO

# endif

  IF(DBG_SET(DBG_SBR)) WRITE(IPT,*)"End: adv_s"

END SUBROUTINE ADV_S
!==============================================================================|
