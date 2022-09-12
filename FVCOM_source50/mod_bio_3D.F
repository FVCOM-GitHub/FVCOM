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

!   THIS PROGRAM LINKS 1D BIOLOGICAL CALCULATION TO FVCOM 3D COMPUTAITION
MODULE MOD_BIO_3D 
#  if defined (BioGen)
   USE ALL_VARS
!   USE MOD_NCDIO
   USE BCS
   USE MOD_OBCS
# if defined (MULTIPROCESSOR)
   USE MOD_PAR   
# endif
#  if defined (WET_DRY) 
   USE MOD_WD
#  endif
#  if defined (SPHERICAL)
   USE MOD_SPHERICAL
#  endif

!   use netcdf

   USE MOD_1D
   USE MOD_PHYTOPLANKTON !,   ONLY: BIO_P,NNP,INP,IRRAD0,PARFRAC
   USE MOD_ZOOPLANKTON !,     ONLY: BIO_Z,NNZ,INZ
   USE MOD_BACTERIA !,        ONLY: BIO_B,NNB,INB
   USE MOD_DETRITUS !,        ONLY: BIO_D,NND,IND
   USE MOD_DOM !,             ONLY: BIO_DOM,NNM,INM
   USE MOD_NUTRIENT !,        ONLY: BIO_N,NNN,INN
   USE MOD_PARAMETER

   IMPLICIT NONE
   SAVE
   REAL(SP), ALLOCATABLE, TARGET ::  BIO_ALL(:,:,:)        !3D BIO_VARIABLES
   REAL(SP), ALLOCATABLE ::  BIO_F(:,:,:)          !FORECASTED VARIABLES
   REAL(SP), ALLOCATABLE ::  BIO_MEAN(:,:,:)       !MEAN VARIABLES
   REAL(SP), ALLOCATABLE ::  XFLUX_OBCB(:,:,:)     !OPEN BOUNDARY FLUX
   REAL(SP), ALLOCATABLE ::  BIO_MEANN(:,:,:)      !MEAN IN CELLS
!************   FOR NETCDF OUTPUT
   INTEGER, ALLOCATABLE :: trcsid(:)
   REAL(SP),ALLOCATABLE :: BIO_VAR_MEAN(:,:,:)

!!$   CHARACTER(LEN=80) STARTUP_BIO_TYPE  !!TYPE OF BIO START

!!$   NAMELIST /NML_BIOLOGICAL/     &
!!$       & STARTUP_BIO_TYPE       ! constant, linear, observed, set values

 
   CONTAINS !------------------------------------------------------------------!
            ! BIO_3D1D             :ADVANCE TMODEL USING GOTM LIBRARIES        !
            ! BIO_ADV              :ADVECTION OF BIOLOGICAL STATE VARIABLES    !
            ! BIO_BCOND            :BOUNDARY CONDITION                         !
            ! BIO_EXCHANGE         :INTERPROCESSOR EXCHANGE                    !
            ! BIO_INITIAL          :INITIALIZATION                             !
            ! BIO_HOT_START        :HOT_START FROM BIO_RESTART.NC              !
            ! NAME_LIST_INITIALIZE_BIO                                         !
            ! NAME_LIST_PRINT_BIO                                              !
            !------------------------------------------------------------------!
SUBROUTINE BIO_3D1D
     IMPLICIT NONE
     SAVE
     INTEGER  :: I,J,K,L
     REAL(SP) :: SPCP,ROSEA,SPRO !,BIO_VAR_MEAN(M,KBM1,NTT)

  if(dbg_set(dbg_sbr)) write(ipt,*)&
       & "Start: BIO_3D1D"

     SPCP  = 4.2174E3_SP                        !HEAT SPECIFIC CAPACITY
     ROSEA = 1.023E3_SP                         !RHO OF SEA WATER
     SPRO=SPCP*ROSEA
!     BIO_VAR_MEAN = 0.0_SP
!---------------------------                    !MAIN LOOP OVER ELEMENTS

#    if defined (ONE_D_MODEL)
     DO I=M,M
#    else
     DO I=1,M
#    endif
       DO K=1,KBM1                              !3D TO 1D FIELD
         DO J=1,NNN
            BIO_N(K,J)=BIO_ALL(I,K,J+INN-1)
         END DO
         DO J=1,NNP
            BIO_P(K,J)=BIO_ALL(I,K,J+INP-1)
         END DO
         DO J=1,NNZ
            BIO_Z(K,J)=BIO_ALL(I,K,J+INZ-1)
         END DO
         DO J=1,NNM
            BIO_DOM(K,J)=BIO_ALL(I,K,J+INM-1)
         END DO
         DO J=1,NNB
            BIO_B(K,J)=BIO_ALL(I,K,J+INB-1)
         END DO
         DO J=1,NND
            BIO_D(K,J)=BIO_ALL(I,K,J+IND-1)
         END DO
         DELTA_D(K)=DZ(I,K)*D(I)                   !LAYER THICKNESS
         DELTA_Z(K)=DZZ(I,K)*D(I)                  !DISTANCE BETWEEN LAYERS
         DEPTH_Z(K)=Z(I,K)*D(I)                    !LAYER CENTER DEPTH
         IRRAD0=-SWRAD(I)*PARFRAC*SPRO/RAMP      !PAR FRACTION
         L_NH4N=30._SPP                          !NITRIFICATION USE	 
         T_BIO(K)=T1(I,K)
       END DO                                    !K=1,KB
       T_STEP=DTI
       CALL ZOOPLANKTON
       CALL PHYTOPLANKTON
       CALL BACTERIA
       CALL DETRITUS
       CALL DOM
       CALL NUTRIENT
        DO K=1,KBM1                                !1D TO 3D FIELD
          DO J=1,NNN
             BIO_ALL(I,K,J+INN-1)=BIO_N(K,J)
          END DO
          DO J=1,NNP
             BIO_ALL(I,K,J+INP-1)=BIO_P(K,J)
          END DO
          DO J=1,NNZ
             BIO_ALL(I,K,J+INZ-1)=BIO_Z(K,J)
          END DO
          DO J=1,NNM
             BIO_ALL(I,K,J+INM-1)=BIO_DOM(K,J)
          END DO
          DO J=1,NNB
             BIO_ALL(I,K,J+INB-1)=BIO_B(K,J)
          END DO
          DO J=1,NND
             BIO_ALL(I,K,J+IND-1)=BIO_D(K,J)
          END DO
        END DO 
      KM_BIO(:)=KH(I,:)
      BIO_VAR(1:KBV,1:NTT)=BIO_ALL(I,1:KBV,1:NTT)
      CALL BIO_MIXING 
      BIO_ALL(I,1:KBV,1:NTT)=BIO_VAR(1:KBV,1:NTT)
      END DO !I=1,M 
#    if !defined (ONE_D_MODEL)
#    if defined (MULTIPROCESSOR)
       IF(PAR)CALL BIO_EXCHANGE
#    endif
      CALL BIO_ADV
#  if defined(MULTIPROCESSOR)
      IF(PAR)THEN
        DO J=1,NTT
         CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,KB,MYID,NPROCS,BIO_F(:,:,J))
        END DO
        CALL BIO_EXCHANGE
      END IF  
#    endif
      CALL BIO_BCOND
      BIO_ALL=BIO_F                                   !UPDATE
#    endif
!    end if defined 1D
      WHERE (BIO_ALL < 0.001) BIO_ALL=0.001

  if(dbg_set(dbg_sbr)) write(ipt,*)&
       & "End: BIO_3D1D"
!      IF(MOD(IINT-1,CDF_INT)==0) CALL BIO_OUT_NETCDF


END SUBROUTINE BIO_3D1D

!=============================================================================!
   SUBROUTINE BIO_ADV  
!=============================================================================!
!                                                                             !
!  This subroutine is used to calculate the advection and horizontal diffusion!
!  terms for the state variables of the adjustable biomodel                   !
!=============================================================================!

   USE ALL_VARS
   USE BCS
   USE MOD_OBCS
   USE MOD_PAR   
   USE MOD_WD
   USE MOD_SPHERICAL
#  if defined (SEMI_IMPLICIT)
   USE MOD_SEMI_IMPLICIT, ONLY : IFCETA 
#  endif
#  if defined (TVD)
   USE MOD_TVD
#  endif

   IMPLICIT NONE

   REAL(SP), DIMENSION(0:MT,KB,ntt)  :: XFLUX,RF,XFLUX_ADV
   REAL(SP), DIMENSION(M)           :: PUPX,PUPY,PVPX,PVPY
   REAL(SP), DIMENSION(M)           :: PFPX,PFPY,PFPXD,PFPYD,VISCOFF
   REAL(SP), DIMENSION(3*(NT),KBM1) :: DTIJ
   REAL(SP), DIMENSION(3*(NT),KBM1) :: UVN
   REAL(SP) :: FFD,FF1   !,X11,Y11,X22,Y22,X33,Y33,TMP1,TMP2,XI,YI
   REAL(SP) :: DXA,DYA,DXB,DYB,FIJ1,FIJ2,UN
   REAL(SP) :: TXX,TYY,FXX,FYY,VISCOF,EXFLUX,TEMP
   REAL(SP) :: FACT,FM1
   REAL(SP) :: TT,TTIME,STPOINT
   INTEGER  :: I,I1,I2,IA,IB,J,J1,J2,JTMP,K,JJ,N1
   REAL(SP) :: WQM1MIN, WQM1MAX, WQM2MIN, WQM2MAX

   REAL(SP), ALLOCATABLE :: DBIODIS(:,:,:)  !!BIOLOGICAL VARIABLES DISCHARGE DATA
!   REAL(SP), ALLOCATABLE :: BIODIS(:,:)     !!BIOLOGICAL VARIABLES AT CURRENT TIME

!!$# if defined (SPHERICAL)
!!$   REAL(SP) :: ty,txpi,typi
!!$   REAL(DP) :: XTMP,XTMP1
!!$   REAL(DP) :: X1_DP,Y1_DP,X2_DP,Y2_DP,XII,YII
!!$   REAL(DP) :: X11_TMP,Y11_TMP,X33_TMP,Y33_TMP
!!$# endif

#  if defined (SEMI_IMPLICIT)
   REAL(SP) :: UN1
   REAL(SP), DIMENSION(3*(NT),KBM1) :: UVN1
   REAL(SP), DIMENSION(3*(NT),KBM1) :: DTIJ1
#  endif
# if defined (TVD)
   REAL(SP) :: S_upstreamA,S_upstreamB
   REAL(SP) :: A_limiter,B_limiter,smoothrA,smoothrB
#  endif
#  if defined (MPDATA)
   REAL(SP) :: WQMMIN,WQMMAX,XXXX
   REAL(SP), DIMENSION(0:MT,KB)     :: WQM_S    !! temporary BIO in modified upwind
   REAL(SP), DIMENSION(0:MT,KB)     :: WQM_SF   !! temporary BIO in modified upwind
   REAL(SP), DIMENSION(0:MT,KB)     :: WWWS     
   REAL(SP), DIMENSION(0:MT,KB)     :: WWWSF   
   REAL(SP), DIMENSION(0:MT)        :: DTWWWS  
   REAL(SP), DIMENSION(0:MT,KB)     :: ZZZFLUX !! temporary total flux in corrected part
   REAL(SP), DIMENSION(0:MT,KB)     :: BETA    !! temporary beta coefficient in corrected part
   REAL(SP), DIMENSION(0:MT,KB)     :: BETAIN  !! temporary beta coefficient in corrected part
   REAL(SP), DIMENSION(0:MT,KB)     :: BETAOUT !! temporary beta coefficient in corrected part
   REAL(SP), DIMENSION(0:MT,KB)     :: BIO_FRESH    !! for source term
   REAL(SP), DIMENSION(0:MT,KB,NTT) :: OFFS     !! Offset to avoid values less than zero.
   INTEGER ITERA, NTERA
#  endif

!------------------------------------------------------------------------------
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
  ! Adding offset to force FVCOM not to crash when variable BIO_ALL approach zero
      OFFS = 1.0_SP
      BIO_ALL = BIO_ALL+OFFS
      BIO_MEAN = BIO_MEAN+OFFS
# endif

!
!--Initialize Fluxes-----------------------------------------------------------
!
   XFLUX = 0.0_SP
   XFLUX_ADV = 0.0_SP
!
!--Loop Over Control Volume Sub-Edges And Calculate Normal Velocity------------
!
   DO I=1,NCV
     I1=NTRG(I)
     DO K=1,KBM1
       DTIJ(I,K) = DT1(I1)*DZ1(I1,K)
       UVN(I,K)=V(I1,K)*DLTXE(I) - U(I1,K)*DLTYE(I) 
#      if defined (SEMI_IMPLICIT)
       DTIJ1(I,K) = D1(I1)*DZ1(I1,K)
       UVN1(I,K) = VF(I1,K)*DLTXE(I) - UF(I1,K)*DLTYE(I)
#      endif
     END DO
   END DO

!!$   TTIME=THOUR

   RF = 0.0_SP

!--Calculate the Advection and Horizontal Diffusion Terms----------------------

   DO N1=1,NTT
     DO K=1,KBM1
       PFPX  = 0.0_SP
       PFPY  = 0.0_SP
       PFPXD = 0.0_SP
       PFPYD = 0.0_SP

       DO I=1,M
         DO J=1,NTSN(I)-1
           I1=NBSN(I,J)
           I2=NBSN(I,J+1)

#    if defined (WET_DRY)
         IF(ISWETN(I1) == 0 .AND. ISWETN(I2) == 1)THEN
          FFD=0.5_SP*(BIO_ALL(I,K,N1)+BIO_ALL(I2,K,N1)           &
	      -BIO_MEAN(I,K,N1)-BIO_MEAN(I2,K,N1))
          FF1=0.5_SP*(BIO_ALL(I,K,N1)+BIO_ALL(I2,K,N1))
	 ELSE IF(ISWETN(I1) == 1 .AND. ISWETN(I2) == 0)THEN
          FFD=0.5_SP*(BIO_ALL(I1,K,N1)+BIO_ALL(I,K,N1)           &
	      -BIO_MEAN(I1,K,N1)-BIO_MEAN(I,K,N1))
          FF1=0.5_SP*(BIO_ALL(I1,K,N1)+BIO_ALL(I,K,N1))
	 ELSE IF(ISWETN(I1) == 0 .AND. ISWETN(I2) == 0)THEN
          FFD=BIO_ALL(I,K,N1)-BIO_MEAN(I,K,N1)
          FF1=BIO_ALL(I,K,N1)
	 ELSE
          FFD=0.5_SP*(BIO_ALL(I1,K,N1)+BIO_ALL(I2,K,N1)          &
	      -BIO_MEAN(I1,K,N1)-BIO_MEAN(I2,K,N1))
          FF1=0.5_SP*(BIO_ALL(I1,K,N1)+BIO_ALL(I2,K,N1))
	 END IF 
#    else	 
           FFD=0.5_SP*(BIO_ALL(I1,K,N1)+BIO_ALL(I2,K,N1)          &
               -BIO_MEAN(I1,K,N1)-BIO_MEAN(I2,K,N1))
           FF1=0.5_SP*(BIO_ALL(I1,K,N1)+BIO_ALL(I2,K,N1))
#    endif

           PFPX(I)=PFPX(I)+FF1*DLTYTRIE(I,J)
           PFPY(I)=PFPY(I)+FF1*DLTXTRIE(I,J)
           PFPXD(I)=PFPXD(I)+FFD*DLTYTRIE(I,J)
           PFPYD(I)=PFPYD(I)+FFD*DLTXTRIE(I,J)
         END DO
         PFPX(I)=PFPX(I)/ART2(I)
         PFPY(I)=PFPY(I)/ART2(I)
         PFPXD(I)=PFPXD(I)/ART2(I)
         PFPYD(I)=PFPYD(I)/ART2(I)
       END DO

      IF(K == KBM1)THEN
        DO I=1,M
          PFPXB(I) = PFPX(I)
          PFPYB(I) = PFPY(I)
        END DO
      END IF


       DO I=1,M
         VISCOFF(I)=VISCOFH(I,K)  !CALCULATED IN viscofh.F
       END DO

       DO I=1,NCV_I
         IA=NIEC(I,1)
         IB=NIEC(I,2)

# if !defined (TVD)
!         FIJ1=BIO_ALL(IA,K,N1)+DXA*PFPX(IA)+DYA*PFPY(IA)
!         FIJ2=BIO_ALL(IB,K,N1)+DXB*PFPX(IB)+DYB*PFPY(IB)
         FIJ1=BIO_ALL(IA,K,N1)+DLTXNCVE(I,1)*PFPX(IA)+DLTYNCVE(I,1)*PFPY(IA)
         FIJ2=BIO_ALL(IB,K,N1)+DLTXNCVE(I,2)*PFPX(IB)+DLTYNCVE(I,2)*PFPY(IB)

         WQM1MIN=MINVAL(BIO_ALL(NBSN(IA,1:NTSN(IA)-1),K,N1))
         WQM1MIN=MIN(WQM1MIN, BIO_ALL(IA,K,N1))
         WQM1MAX=MAXVAL(BIO_ALL(NBSN(IA,1:NTSN(IA)-1),K,N1))
         WQM1MAX=MAX(WQM1MAX, BIO_ALL(IA,K,N1))
         WQM2MIN=MINVAL(BIO_ALL(NBSN(IB,1:NTSN(IB)-1),K,N1))
         WQM2MIN=MIN(WQM2MIN, BIO_ALL(IB,K,N1))
         WQM2MAX=MAXVAL(BIO_ALL(NBSN(IB,1:NTSN(IB)-1),K,N1))
         WQM2MAX=MAX(WQM2MAX, BIO_ALL(IB,K,N1))
         IF(FIJ1 < WQM1MIN) FIJ1=WQM1MIN
         IF(FIJ1 > WQM1MAX) FIJ1=WQM1MAX
         IF(FIJ2 < WQM2MIN) FIJ2=WQM2MIN
         IF(FIJ2 > WQM2MAX) FIJ2=WQM2MAX

# else
! ------------------------------------------------------------------------------------------ !
!                                     Drawing the TVD scheme
! ------------------------------------------------------------------------------------------ !
!     If A is upstream of B
       S_upstreamA=BIO_ALL(Anear_node(I),K,N1)+PFPX(Anear_node(I))*XUAdist(I)+PFPY(Anear_node(I))*YUAdist(I)
!     If B is upstream of A
       S_upstreamB=BIO_ALL(Bnear_node(I),K,N1)+PFPX(Bnear_node(I))*XUBdist(I)+PFPY(Bnear_node(I))*YUBdist(I)                

!     Smoothness-parameter
       IF (BIO_ALL(IA,K,N1).EQ.BIO_ALL(IB,K,N1)) THEN
          smoothrA    = 0.0_SP
          smoothrB    = 0.0_SP
       ELSE IF (BIO_ALL(IA,K,N1).LT.1.E-10.AND.BIO_ALL(IB,K,N1).LT.1.E-10) THEN
          smoothrA    = 0.0_SP
          smoothrB    = 0.0_SP
       ELSE
          smoothrA    = (BIO_ALL(IA,K,N1)-S_upstreamA)/(BIO_ALL(IB,K,N1)-BIO_ALL(IA,K,N1))
          smoothrB    = (BIO_ALL(IB,K,N1)-S_upstreamB)/(BIO_ALL(IA,K,N1)-BIO_ALL(IB,K,N1))
       END IF

!     Flux-limiter
!     (superbee) - Feel free to use other limiters!
       A_limiter = MAX(0.0_SP, MIN(1.0_SP, 2.0_SP*smoothrA), MIN(2.0_SP, smoothrA))
       B_limiter = MAX(0.0_SP, MIN(1.0_SP, 2.0_SP*smoothrB), MIN(2.0_SP, smoothrB))

!     Estimating the variable at the cellwall in-betweem IA and IB
       FIJ1     = BIO_ALL(IA,K,N1)+0.5_SP*A_LIMITER*(BIO_ALL(IB,K,N1)-BIO_ALL(IA,K,N1))
       FIJ2     = BIO_ALL(IB,K,N1)+0.5_SP*B_LIMITER*(BIO_ALL(IA,K,N1)-BIO_ALL(IB,K,N1))

# endif
         UN=UVN(I,K)
#        if defined (SEMI_IMPLICIT)
         UN1=UVN1(I,K)
#        endif
  
!         VISCOF=HORCON*(FACT*(VISCOFF(IA)+VISCOFF(IB))*0.5_SP + FM1)
        ! David moved HPRNU and added HVC
        VISCOF=(FACT*0.5_SP*(VISCOFF(IA)*NN_HVC(IA)+VISCOFF(IB)*NN_HVC(IB)) + FM1*0.5_SP*(NN_HVC(IA)+NN_HVC(IB)))

         TXX=0.5_SP*(PFPXD(IA)+PFPXD(IB))*VISCOF
         TYY=0.5_SP*(PFPYD(IA)+PFPYD(IB))*VISCOF

         FXX=-DTIJ(I,K)*TXX*DLTYE(I)
         FYY= DTIJ(I,K)*TYY*DLTXE(I)

#        if !defined (SEMI_IMPLICIT)
         EXFLUX=-UN*DTIJ(I,K)*                           &
                ((1.0_SP+SIGN(1.0_SP,UN))*FIJ2+          &
                 (1.0_SP-SIGN(1.0_SP,UN))*FIJ1)*0.5_SP   &
                +FXX+FYY
#        else
         EXFLUX=-UN*DTIJ(I,K)* &
            ((1.0_SP+SIGN(1.0_SP,UN))*FIJ2+(1.0_SP-SIGN(1.0_SP,UN))*FIJ1)*0.5_SP
         EXFLUX=(1.0_SP-IFCETA)*EXFLUX+IFCETA*(-UN1*DTIJ1(I,K)*((1.0_SP+SIGN(1.0_SP,UN1))*FIJ2+  &
            (1.0_SP-SIGN(1.0_SP,UN1))*FIJ1)*0.5_SP)+FXX+FYY
#        endif
 
         XFLUX(IA,K,N1)=XFLUX(IA,K,N1)+EXFLUX
         XFLUX(IB,K,N1)=XFLUX(IB,K,N1)-EXFLUX

       XFLUX_ADV(IA,K,N1)=XFLUX_ADV(IA,K,N1)+(EXFLUX-FXX-FYY)
       XFLUX_ADV(IB,K,N1)=XFLUX_ADV(IB,K,N1)-(EXFLUX-FXX-FYY)

       END DO !to M

#    if defined (SPHERICAL)
!     CALL ADV_T_XY(XFLUX(:,:,N1),XFLUX_ADV(:,:,N1),PTPX,PTPY,PTPXD,PTPYD,VISCOFF,K)
#    endif

     END DO !to KBM1
   END DO !to nnt

!
!-Accumulate Fluxes at Boundary Nodes
!
# if defined (MULTIPROCESSOR)
      DO N1=1,NTT
!      IF(PAR)CALL NODE_MATCH(0,NBN,BN_MLT,BN_LOC,BNC,MT,KB,MYID,NPROCS,       &
!                            XFLUX(:,:,I))    
   IF(PAR)CALL NODE_MATCH(0,NBN,BN_MLT,BN_LOC,BNC,MT,KB,MYID,NPROCS,       &
                            XFLUX(:,:,N1),XFLUX_ADV(:,:,N1))
      END DO
# endif

!#  if !defined (MPDATA)
   DO N1=1,NTT
     DO K=1,KBM1
        IF(IOBCN > 0) THEN
          DO I=1,IOBCN
            I1=I_OBC_N(I)
            XFLUX_OBCB(I,K,N1)=XFLUX_ADV(I1,K,N1)
          END DO
        END IF
      END DO
    END DO
!# endif


   DO N1=1,ntt
!#  if !defined (MPDATA)
!
!--Calculate the Vertical Terms------------------------------------------------
!
!     DO K=1,KBM1
!       DO I=1,M
!#      if defined (WET_DRY)
!       IF(ISWETN(I)*ISWETNT(I) == 1) THEN
!#      endif
!         IF(K == 1) THEN  !Is there any violation ?
!           TEMP=-WTS(I,K+1)*(BIO_ALL(I,K,N1)*DZ(I,K+1)+BIO_ALL(I,K+1,N1)*DZ(I,K))/   &
!	         (DZ(I,K)+DZ(I,K+1))
!         ELSE IF(K == KBM1) THEN
!           TEMP=WTS(I,K)*(BIO_ALL(I,K,N1)*DZ(I,K-1)+BIO_ALL(I,K-1,N1)*DZ(I,K))/      &
!	         (DZ(I,K)+DZ(I,K-1))
!         ELSE
!           TEMP=WTS(I,K)*(BIO_ALL(I,K,N1)*DZ(I,K-1)+BIO_ALL(I,K-1,N1)*DZ(I,K))/      &
!	         (DZ(I,K)+DZ(I,K-1))-  &
!                WTS(I,K+1)*(BIO_ALL(I,K,N1)*DZ(I,K+1)+BIO_ALL(I,K+1,N1)*DZ(I,K))/    &
!		 (DZ(I,K)+DZ(I,K+1))
!         END IF

!
!--Total Fluxes ---------------------------------------------------------------
!
!         IF(ISONB(I) == 2) THEN
!           XFLUX(I,K,N1)=TEMP*ART1(I)
!         ELSE
!           XFLUX(I,K,N1)=XFLUX(I,K,N1)+TEMP*ART1(I)
!         END IF
!#    if defined (WET_DRY)
!       END IF
!#    endif
!       END DO  !i=1,M
!     END DO    !k=1,kbm1

!--Set Boundary Conditions-For Fresh Water Flux--------------------------------!
!
!     IF(RIVER_TS_SETTING == 'calculated') THEN
!       IF(RIVER_INFLOW_LOCATION == 'node') THEN
!         IF(NUMQBC > 0) THEN
!           DO J=1,NUMQBC
!             JJ=INODEQ(J)
!             STPOINT=BIODIS(J,N1)
!             DO K=1,KBM1
!               XFLUX(JJ,K,N1)=XFLUX(JJ,K,N1) - QDIS(J)*VQDIST(J,K)*STPOINT
!             END DO
!           END DO
!         END IF
!       ELSE IF(RIVER_INFLOW_LOCATION == 'edge') THEN
!         IF(NUMQBC > 0) THEN
!           DO J=1,NUMQBC
!             J1=N_ICELLQ(J,1)
!             J2=N_ICELLQ(J,2)
!             STPOINT=BIODIS(J,N1) !!ASK LIU SHOULD THIS BE STPOINT1(J1)/STPOINT2(J2)
!             DO K=1,KBM1
!               XFLUX(J1,K,N1)=XFLUX(J1,K,N1)-QDIS(J)*RDISQ(J,1)*VQDIST(J,K)*STPOINT
!               XFLUX(J2,K,N1)=XFLUX(J2,K,N1)-QDIS(J)*RDISQ(J,2)*VQDIST(J,K)*STPOINT
!             END DO
!           END DO
!         END IF
!       END IF
!     END IF

!#  else
#   if defined (MPDATA)

!--------------------------------------------------------------------------------
!   Using smolarkiewicz, P. K; A fully multidimensional positive definite advection
!   TEMPport algorithm with small implicit diffusion, Journal of Computational
!   Physics, 54, 325-362, 1984
!-----------------------------------------------------------------        

   BIO_FRESH=BIO_ALL(:,:,N1)

   IF(RIVER_TS_SETTING == 'calculated') THEN
     IF(RIVER_INFLOW_LOCATION == 'node') THEN
       IF(NUMQBC > 0) THEN
         DO J=1,NUMQBC
           JJ=INODEQ(J)
           STPOINT=BIODIS(J,N1)+OFFS(1,1,1)
           DO K=1,KBM1
	    BIO_FRESH(JJ,K)=BIODIS(J,N1)+ OFFS(1,1,1)  
            XFLUX(JJ,K,N1)=XFLUX(JJ,K,N1) - QDIS(J)*VQDIST(J,K)*STPOINT
           END DO
         END DO
       END IF
     ELSE IF(RIVER_INFLOW_LOCATION == 'edge') THEN
       IF(NUMQBC > 0) THEN
         DO J=1,NUMQBC
           J1=N_ICELLQ(J,1)
           J2=N_ICELLQ(J,2)
           STPOINT=BIODIS(J,N1)+OFFS(1,1,1) 
           DO K=1,KBM1
             BIO_FRESH(J1,K)=BIODIS(J,N1)+OFFS(1,1,1)
             BIO_FRESH(J2,K)=BIODIS(J,N1)+OFFS(1,1,1)
             XFLUX(J1,K,N1)=XFLUX(J1,K,N1)-QDIS(J)*RDISQ(J,1)*VQDIST(J,K)*STPOINT
             XFLUX(J2,K,N1)=XFLUX(J2,K,N1)-QDIS(J)*RDISQ(J,2)*VQDIST(J,K)*STPOINT
           END DO
         END DO
       END IF
     END IF
   END IF
!

! The horizontal term of advection is neglected here
   DO K=1,KBM1
     DO I=1,M
       IF(ISONB(I) == 2) THEN
         XFLUX(I,K,N1)=0._SP
       ENDIF
     END DO
   END DO

! Initialize variables of MPDATA
   WQM_S=0._SP
   WQM_SF=0._SP
   WWWS=0._SP
   WWWSF=0._SP
   DTWWWS=0._SP
   ZZZFLUX=0._SP
   BETA=0._SP
   BETAIN=0._SP
   BETAOUT=0._SP

!!   first loop for vertical upwind
!!   flux including horizontal and vertical upwind
   DO K=1,KBM1
     DO I=1,M
#    if defined (WET_DRY)
       IF(ISWETN(I)*ISWETNT(I) == 1) THEN
#    endif
         IF(K == 1) THEN
           TEMP = -(WTS(I,K+1)-ABS(WTS(I,K+1)))*BIO_ALL(I,K,N1)   &
                  -(WTS(I,K+1)+ABS(WTS(I,K+1)))*BIO_ALL(I,K+1,N1) &
                  +(WTS(I,K)+ABS(WTS(I,K)))*BIO_ALL(I,K,N1)    
         ELSE IF(K == KBM1) THEN
           TEMP = +(WTS(I,K)-ABS(WTS(I,K)))*BIO_ALL(I,K-1,N1)     &
                  +(WTS(I,K)+ABS(WTS(I,K)))*BIO_ALL(I,K,N1)
         ELSE
           TEMP = -(WTS(I,K+1)-ABS(WTS(I,K+1)))*BIO_ALL(I,K,N1)   &
                  -(WTS(I,K+1)+ABS(WTS(I,K+1)))*BIO_ALL(I,K+1,N1) &
                  +(WTS(I,K)-ABS(WTS(I,K)))*BIO_ALL(I,K-1,N1)     &
                  +(WTS(I,K)+ABS(WTS(I,K)))*BIO_ALL(I,K,N1)
         END IF
         TEMP = 0.5_SP*TEMP 

         IF(K == 1)THEN
           WQMMAX = MAXVAL(BIO_ALL(NBSN(I,1:NTSN(I)),K,N1))
           WQMMIN = MINVAL(BIO_ALL(NBSN(I,1:NTSN(I)),K,N1))
           WQMMAX = MAX(WQMMAX,BIO_ALL(I,K+1,N1),BIO_ALL(I,K,N1),BIO_FRESH(I,K))
           WQMMIN = MIN(WQMMIN,BIO_ALL(I,K+1,N1),BIO_ALL(I,K,N1),BIO_FRESH(I,K))
         ELSEIF(K == KBM1) THEN
           WQMMAX = MAXVAL(BIO_ALL(NBSN(I,1:NTSN(I)),K,N1))
           WQMMIN = MINVAL(BIO_ALL(NBSN(I,1:NTSN(I)),K,N1))
           WQMMAX = MAX(WQMMAX,BIO_ALL(I,K-1,N1),BIO_ALL(I,K,N1),BIO_FRESH(I,K))
           WQMMIN = MIN(WQMMIN,BIO_ALL(I,K-1,N1),BIO_ALL(I,K,N1),BIO_FRESH(I,K))
         ELSE
           WQMMAX = MAXVAL(BIO_ALL(NBSN(I,1:NTSN(I)),K,N1))
           WQMMIN = MINVAL(BIO_ALL(NBSN(I,1:NTSN(I)),K,N1))
           WQMMAX = MAX(WQMMAX,BIO_ALL(I,K+1,N1),BIO_ALL(I,K-1,N1),BIO_ALL(I,K,N1),BIO_FRESH(I,K))
           WQMMIN = MIN(WQMMIN,BIO_ALL(I,K+1,N1),BIO_ALL(I,K-1,N1),BIO_ALL(I,K,N1),BIO_FRESH(I,K))
         END IF

         ! Total (horizontal + vertical) flux to the cell
         ZZZFLUX(I,K) = TEMP*(DTI/DT(I))/DZ(I,K) + XFLUX(I,K,N1)/ART1(I)*(DTI/DT(I))/DZ(I,K) 

         ! Updated variable without limiter minus current variable
         XXXX = ZZZFLUX(I,K)*DT(I)/DTFA(I)+BIO_ALL(I,K,N1)-BIO_ALL(I,K,N1)*DT(I)/DTFA(I) 

           ! Added by Akvaplan 2018 to avoid single precision-crash
         IF((ABS(XXXX).LT.ABS(WQMMAX-BIO_ALL(I,K,N1))).AND.(XXXX.LT.0.0_SP)) THEN
            WQM_SF(I,K) = BIO_ALL(I,K,N1)-XXXX
         ELSE IF((ABS(XXXX).LT.ABS(WQMMIN-BIO_ALL(I,K,N1))).AND.(XXXX.GT.0.0_SP)) THEN 
            WQM_SF(I,K) = BIO_ALL(I,K,N1)-XXXX
         ELSE IF(XXXX.EQ.0) THEN
            WQM_SF(I,K) = BIO_ALL(I,K,N1)
         ELSE
            BETA(I,K)=0.5*(1.-SIGN(1._SP,XXXX)) * (WQMMAX-BIO_ALL(I,K,N1))/(ABS(XXXX)+1.E-10) &                               
                 +0.5*(1.-SIGN(1._SP,-XXXX)) * (BIO_ALL(I,K,N1)-WQMMIN)/(ABS(XXXX)+1.E-10)                                     
            WQM_SF(I,K)=BIO_ALL(I,K,N1)-MIN(1._SP,BETA(I,K))*XXXX
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
       WQM_S   = WQM_SF
       DTWWWS = DT
     ELSE
       WWWSF  = WWWS
       WQM_S   = WQM_SF
       DTWWWS = DTFA
     END IF
     DO K=2,KBM1
       DO I=1,M
         TEMP=ABS(WWWSF(I,K))-DTI*(WWWSF(I,K))*(WWWSF(I,K))/DZ(I,K)/DTWWWS(I)
         WWWS(I,K)=TEMP*(WQM_S(I,K-1)-WQM_S(I,K))/(ABS(WQM_S(I,K-1))+ABS(WQM_S(I,K))+1.E-14)

         IF(TEMP < 0.0_SP .OR. WQM_S(I,K) == 0.0_SP)THEN 
           WWWS(I,K)=0.0_SP 
         END IF
       END DO 
     END DO

     DO I=1,M
       WWWS(I,1)=0.0_SP
     END DO

     DO I=1,M
       WQMMAX = MAXVAL(BIO_ALL(NBSN(I,1:NTSN(I)),1,N1))
       WQMMIN = MINVAL(BIO_ALL(NBSN(I,1:NTSN(I)),1,N1))
       WQMMAX = MAX(WQMMAX,BIO_ALL(I,2,N1),BIO_ALL(I,1,N1),BIO_FRESH(I,1))
       WQMMIN = MIN(WQMMIN,BIO_ALL(I,2,N1),BIO_ALL(I,1,N1),BIO_FRESH(I,1))
 
       TEMP=0.5*((WWWS(I,2)+ABS(WWWS(I,2)))*WQM_S(I,2))*(DTI/DTFA(I))/DZ(I,1)
       BETAIN(I,1)=(WQMMAX-WQM_S(I,1))/(TEMP+1.E-10)

       TEMP=0.5*((WWWS(I,1)+ABS(WWWS(I,1)))*WQM_S(I,1)-        &
	         (WWWS(I,2)-ABS(WWWS(I,2)))*WQM_S(I,1))*(DTI/DTFA(I))/DZ(I,1)
       BETAOUT(I,1)=(WQM_S(I,1)-WQMMIN)/(TEMP+1.E-10)

       WWWSF(I,1)=0.5*MIN(1.,BETAOUT(I,1))*(WWWS(I,1)+ABS(WWWS(I,1))) + &
                  0.5*MIN(1.,BETAIN(I,1))*(WWWS(I,1)-ABS(WWWS(I,1)))
     END DO

     DO K=2,KBM1-1
       DO I=1,M
         WQMMAX = MAXVAL(BIO_ALL(NBSN(I,1:NTSN(I)),K,N1))
         WQMMIN = MINVAL(BIO_ALL(NBSN(I,1:NTSN(I)),K,N1))
         WQMMAX = MAX(WQMMAX,BIO_ALL(I,K+1,N1),BIO_ALL(I,K-1,N1),BIO_ALL(I,K,N1),BIO_FRESH(I,K))
         WQMMIN = MIN(WQMMIN,BIO_ALL(I,K+1,N1),BIO_ALL(I,K-1,N1),BIO_ALL(I,K,N1),BIO_FRESH(I,K))

         TEMP=0.5*((WWWS(I,K+1)+ABS(WWWS(I,K+1)))*WQM_S(I,K+1)-  &
	           (WWWS(I,K)-ABS(WWWS(I,K)))*WQM_S(I,K-1))*(DTI/DTFA(I))/DZ(I,K)
         BETAIN(I,K)=(WQMMAX-WQM_S(I,K))/(TEMP+1.E-10)

         TEMP=0.5*((WWWS(I,K)+ABS(WWWS(I,K)))*WQM_S(I,K)-        &
	           (WWWS(I,K+1)-ABS(WWWS(I,K+1)))*WQM_S(I,K))*(DTI/DTFA(I))/DZ(I,K)
         BETAOUT(I,K)=(WQM_S(I,K)-WQMMIN)/(TEMP+1.E-10)

         WWWSF(I,K)=0.5*MIN(1.,BETAIN(I,K-1),BETAOUT(I,K))*(WWWS(I,K)+ABS(WWWS(I,K))) + &
                    0.5*MIN(1.,BETAIN(I,K),BETAOUT(I,K-1))*(WWWS(I,K)-ABS(WWWS(I,K)))
       END DO
     END DO

     K=KBM1
     DO I=1,M
       WQMMAX = MAXVAL(BIO_ALL(NBSN(I,1:NTSN(I)),K,N1))
       WQMMIN = MINVAL(BIO_ALL(NBSN(I,1:NTSN(I)),K,N1))
       WQMMAX = MAX(WQMMAX,BIO_ALL(I,K-1,N1),BIO_ALL(I,K,N1),BIO_FRESH(I,K))
       WQMMIN = MIN(WQMMIN,BIO_ALL(I,K-1,N1),BIO_ALL(I,K,N1),BIO_FRESH(I,K))

       TEMP=0.5*((WWWS(I,K+1)+ABS(WWWS(I,K+1)))*WQM_S(I,K+1)-  &
	         (WWWS(I,K)-ABS(WWWS(I,K)))*WQM_S(I,K-1))*(DTI/DTFA(I))/DZ(I,K)
       BETAIN(I,K)=(WQMMAX-WQM_S(I,K))/(TEMP+1.E-10)

       TEMP=0.5*((WWWS(I,K)+ABS(WWWS(I,K)))*WQM_S(I,K)-        &
	         (WWWS(I,K+1)-ABS(WWWS(I,K+1)))*WQM_S(I,K))*(DTI/DTFA(I))/DZ(I,K)
       BETAOUT(I,K)=(WQM_S(I,K)-WQMMIN)/(TEMP+1.E-10)

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
             TEMP = -(WWWS(I,K+1)-ABS(WWWS(I,K+1)))*WQM_S(I,K)   &
                    -(WWWS(I,K+1)+ABS(WWWS(I,K+1)))*WQM_S(I,K+1) &
                    +(WWWS(I,K)+ABS(WWWS(I,K)))*WQM_S(I,K)
           ELSE IF(K == KBM1) THEN
             TEMP = +(WWWS(I,K)-ABS(WWWS(I,K)))*WQM_S(I,K-1)     &
                    +(WWWS(I,K)+ABS(WWWS(I,K)))*WQM_S(I,K)
           ELSE
             TEMP = -(WWWS(I,K+1)-ABS(WWWS(I,K+1)))*WQM_S(I,K)   &
                    -(WWWS(I,K+1)+ABS(WWWS(I,K+1)))*WQM_S(I,K+1) &
                    +(WWWS(I,K)-ABS(WWWS(I,K)))*WQM_S(I,K-1)     &
                    +(WWWS(I,K)+ABS(WWWS(I,K)))*WQM_S(I,K)
           END IF
           TEMP = 0.5_SP*TEMP
           WQM_SF(I,K)=(WQM_S(I,K)-TEMP*(DTI/DTFA(I))/DZ(I,K)) 
#      if defined (WET_DRY)
         END IF
#      endif
       END DO
     END DO  !! SIGMA LOOP
   END DO  !! Smolarvizw Loop
!--------------------------------------------------------------------------
! End of smolarkiewicz upwind loop
!--------------------------------------------------------------------------
      BIO_ALL(:,:,N1) = BIO_ALL(:,:,N1)-OFFS(1,1,1)
      WQM_SF = WQM_SF-OFFS(1,1,1)
      BIO_MEAN(:,:,N1) = BIO_MEAN(:,:,N1)-OFFS(1,1,1)
     
#  endif

#  if !defined (MPDATA)
!
!--Calculate the Vertical Terms------------------------------------------------
!
     DO K=1,KBM1
       DO I=1,M
#      if defined (WET_DRY)
       IF(ISWETN(I)*ISWETNT(I) == 1) THEN
#      endif
         IF(K == 1) THEN  
           TEMP=-WTS(I,K+1)*(BIO_ALL(I,K,N1)*DZ(I,K+1)+BIO_ALL(I,K+1,N1)*DZ(I,K))/   &
	         (DZ(I,K)+DZ(I,K+1))
         ELSE IF(K == KBM1) THEN
           TEMP=WTS(I,K)*(BIO_ALL(I,K,N1)*DZ(I,K-1)+BIO_ALL(I,K-1,N1)*DZ(I,K))/      &
	         (DZ(I,K)+DZ(I,K-1))
         ELSE
           TEMP=WTS(I,K)*(BIO_ALL(I,K,N1)*DZ(I,K-1)+BIO_ALL(I,K-1,N1)*DZ(I,K))/      &
	         (DZ(I,K)+DZ(I,K-1))-  &
                WTS(I,K+1)*(BIO_ALL(I,K,N1)*DZ(I,K+1)+BIO_ALL(I,K+1,N1)*DZ(I,K))/    &
		 (DZ(I,K)+DZ(I,K+1))
         END IF

!
!--Total Fluxes ---------------------------------------------------------------
!
         IF(ISONB(I) == 2) THEN
           XFLUX(I,K,N1)=TEMP*ART1(I)
         ELSE
           XFLUX(I,K,N1)=XFLUX(I,K,N1)+TEMP*ART1(I)
         END IF
#    if defined (WET_DRY)
       END IF
#    endif
       END DO  !i=1,M
     END DO    !k=1,kbm1

!--Set Boundary Conditions-For Fresh Water Flux--------------------------------!
!
     IF(RIVER_TS_SETTING == 'calculated') THEN
       IF(RIVER_INFLOW_LOCATION == 'node') THEN
         IF(NUMQBC > 0) THEN
           DO J=1,NUMQBC
             JJ=INODEQ(J)
             STPOINT=BIODIS(J,N1)
             DO K=1,KBM1
               XFLUX(JJ,K,N1)=XFLUX(JJ,K,N1) - QDIS(J)*VQDIST(J,K)*STPOINT
             END DO
           END DO
         END IF
       ELSE IF(RIVER_INFLOW_LOCATION == 'edge') THEN
         IF(NUMQBC > 0) THEN
           DO J=1,NUMQBC
             J1=N_ICELLQ(J,1)
             J2=N_ICELLQ(J,2)
             STPOINT=BIODIS(J,N1) 
             DO K=1,KBM1
               XFLUX(J1,K,N1)=XFLUX(J1,K,N1)-QDIS(J)*RDISQ(J,1)*VQDIST(J,K)*STPOINT
               XFLUX(J2,K,N1)=XFLUX(J2,K,N1)-QDIS(J)*RDISQ(J,2)*VQDIST(J,K)*STPOINT
             END DO
           END DO
         END IF
       END IF
     END IF
# endif

!--Update Variables--------------------------------
!

#    if defined (WET_DRY)
     DO I = 1,M
       IF(ISWETN(I)*ISWETNT(I) == 1 )THEN
         DO K = 1, KBM1
       XFLUX(I,K,N1) = XFLUX(I,K,N1) - RF(I,K,N1)*ART1(I)     !/DZ(K)
#    if !defined (MPDATA) 
           BIO_F(I,K,N1)=(BIO_ALL(I,K,N1)-XFLUX(I,K,N1)/ART1(I)*(DTI/(DT(I)*DZ(I,K))))*   &
                         (DT(I)/DTFA(I)) 
!JQI 04/29/2019 D(I) -> DTFA   (DT(I)/D(I)) 
#    else
           BIO_F(I,K,N1)=WQM_SF(I,K)
#    endif  		 			 
         END DO 
       ELSE
         DO K=1,KBM1
           BIO_F(I,K,N1)=BIO_ALL(I,K,N1)
         END DO
       END IF
     END DO
#    else
     DO I = 1,M
       DO K = 1, KBM1
         XFLUX(I,K,N1) = XFLUX(I,K,N1) - RF(I,K,N1)*ART1(I)     !/DZ(K)
#    if !defined (MPDATA) 
         BIO_F(I,K,N1)=(BIO_ALL(I,K,N1)-XFLUX(I,K,N1)/ART1(I)*(DTI/(DT(I)*DZ(I,K))))*   &
                       (DT(I)/DTFA(I)) 
!JQI 04/29/2019 D(I) -> DTFA    (DT(I)/D(I)) 
#    else
         BIO_F(I,K,N1)=WQM_SF(I,K)
#    endif  		 			 
       END DO 
     END DO	 
#    endif
   END DO !do N1=1,ntt
   RETURN
   END SUBROUTINE BIO_ADV



   SUBROUTINE BIO_BCOND
!==============================================================================|
!   Set Boundary Conditions for BioGen                                         |
!==============================================================================|

!------------------------------------------------------------------------------|
   USE ALL_VARS
   USE BCS
   USE MOD_OBCS
   IMPLICIT NONE
   REAL(SP) :: T2D,T2D_NEXT,T2D_OBC,XFLUX2D,TMP,RAMP_BIO
   INTEGER  :: I,J,K,J1,J11,J22,NCON2,N1
!   REAL(SP), ALLOCATABLE :: WDIS(:,:)     !!FRESH WATER QUALITY AT CURRENT TIME
   REAL(SP) ::WQMMAX,WQMMIN
   CHARACTER(LEN=80) RIVER_TS_SETTING_TMP

!   ALLOCATE(WDIS(NUMQBC,ntt))     ;WDIS      = ZERO
!------------------------------------------------------------------------------|


!
!--SET CONDITIONS FOR FRESH WATER INFLOW---------------------------------------|
!
   IF(RIVER_TS_SETTING == 'specified') THEN
     IF(NUMQBC > 0) THEN
       IF(RIVER_INFLOW_LOCATION == 'node') THEN
         DO I=1,NUMQBC
           J11=INODEQ(I)
           DO K=1,KBM1
             DO N1=1,NTT
               BIO_F(J11,K,N1) = BIODIS(I,N1)
             END DO
           END DO
         END DO
       ELSE IF(RIVER_INFLOW_LOCATION == 'edge') THEN
         DO I=1,NUMQBC
           J11=N_ICELLQ(I,1)
           J22=N_ICELLQ(I,2)
           DO K=1,KBM1
             DO N1=1,NTT
               BIO_F(J11,K,N1)=BIODIS(I,N1)
               BIO_F(J22,K,N1)=BIODIS(I,N1)
             END DO
           END DO
         END DO
       END IF
     END IF
   END IF

       
   IF(IOBCN > 0) THEN
!
!  SET CONDITIONS ON OUTER BOUNDARY
!
     IF(OBC_BIO_NUDGING) CALL UPDATE_OBC_BIO(IntTime,BIO_OBC)
     RAMP_BIO = TANH(FLOAT(IINT)/FLOAT(IRAMP+1))
     DO N1=1,NTT
       DO I=1,IOBCN
         J=I_OBC_N(I)
         J1=NEXT_OBC(I)
         T2D=0.0_SP
         T2D_NEXT=0.0_SP
         XFLUX2D=0.0_SP
         DO K=1,KBM1
           T2D=T2D+BIO_ALL(J,K,N1)*DZ(J,K)
           T2D_NEXT=T2D_NEXT+BIO_F(J1,K,N1)*DZ(J1,K)
           XFLUX2D=XFLUX2D+XFLUX_OBCB(I,K,N1)           !*DZ(K)
         END DO
         IF(UARD_OBCN(I) > 0.0_SP) THEN
           TMP=XFLUX2D+T2D*UARD_OBCN(I)
           T2D_OBC=(T2D*DT(J)-TMP*DTI/ART1(J))/D(J)
           DO K=1,KBM1
!            BIO_ALL(J,K,N1)=T2D_OBC+(BIO_ALL(J1,K,N1)-T2D_NEXT)
            BIO_F(J,K,N1)=T2D_OBC+(BIO_F(J1,K,N1)-T2D_NEXT)
           END DO

           DO K=1,KBM1
             WQMMAX = MAXVAL(BIO_ALL(NBSN(J,1:NTSN(J)),K,N1))
             WQMMIN = MINVAL(BIO_ALL(NBSN(J,1:NTSN(J)),K,N1))
         
             IF(K == 1)THEN
               WQMMAX = MAX(WQMMAX,(BIO_ALL(J,K,N1)*DZ(J,K+1)+BIO_ALL(J,K+1,N1)*DZ(J,K))/  &
	               (DZ(J,K)+DZ(J,K+1)))
               WQMMIN = MIN(WQMMIN,(BIO_ALL(J,K,N1)*DZ(J,K+1)+BIO_ALL(J,K+1,N1)*DZ(J,K))/  &
	               (DZ(J,K)+DZ(J,K+1)))
             ELSE IF(K == KBM1)THEN
               WQMMAX = MAX(WQMMAX,(BIO_ALL(J,K,N1)*DZ(J,K-1)+BIO_ALL(J,K-1,N1)*DZ(J,K))/  &
	               (DZ(J,K)+DZ(J,K-1)))
               WQMMIN = MIN(WQMMIN,(BIO_ALL(J,K,N1)*DZ(J,K-1)+BIO_ALL(J,K-1,N1)*DZ(J,K))/  &
	               (DZ(J,K)+DZ(J,K-1)))
             ELSE
               WQMMAX = MAX(WQMMAX,(BIO_ALL(J,K,N1)*DZ(J,K-1)+BIO_ALL(J,K-1,N1)*DZ(J,K))/  &
	               (DZ(J,K)+DZ(J,K-1)), &
                       (BIO_ALL(J,K,N1)*DZ(J,K+1)+BIO_ALL(J,K+1,N1)*DZ(J,K))/  &
		       (DZ(J,K)+DZ(J,K+1)))
               WQMMIN = MIN(WQMMIN,(BIO_ALL(J,K,N1)*DZ(J,K-1)+BIO_ALL(J,K-1,N1)*DZ(J,K))/  &
	               (DZ(J,K)+DZ(J,K-1)), &
                       (BIO_ALL(J,K,N1)*DZ(J,K+1)+BIO_ALL(J,K+1,N1)*DZ(J,K))/  &
		       (DZ(J,K)+DZ(J,K+1)))
             END IF
 
             IF(WQMMIN-BIO_F(J,K,N1) > 0.0_SP)BIO_F(J,K,N1) = WQMMIN
             IF(BIO_F(J,K,N1)-WQMMAX > 0.0_SP)BIO_F(J,K,N1) = WQMMAX

           END DO

         ELSE

           IF(OBC_BIO_NUDGING) THEN
             DO K=1,KBM1
               BIO_F(J,K,N1) = BIO_ALL(J,K,N1) - OBC_BIO_NUDGING_TIMESCALE*RAMP_BIO*(BIO_ALL(J,K,N1)&
                      &-BIO_OBC(I,K,N1))
             END DO
           ELSE
             DO K=1,KBM1
               BIO_F(J,K,N1)=BIO_ALL(J,K,N1)
             END DO
           END IF
         END IF
       END DO
     END DO !!OUTER LOOP OVER BIO-VARIABLES

   END IF

!
!--SET BOUNDARY CONDITIONS-----------------------------------------------------|
!
       BIO_ALL(0,:,:)=ZERO
   RETURN
   END SUBROUTINE BIO_BCOND
!==============================================================================!

SUBROUTINE BIO_EXCHANGE
!==============================================================================!
!     PERFORM DATA EXCHANGE FOR the Generalized biological model               |
!==============================================================================!
#if defined (MULTIPROCESSOR)
!     USE ALL_VARS
     USE MOD_PAR
     USE LIMS
     USE CONTROL
     IMPLICIT NONE
     INTEGER :: I3
     REAL(SP),ALLOCATABLE :: BIO_ALL_T(:,:),BIO_MEAN_T(:,:),BIO_F_T(:,:)
     DO I3=1,NTT
      ALLOCATE(BIO_ALL_T(0:MT,KB))
      ALLOCATE(BIO_MEAN_T(0:MT,KB))
      ALLOCATE(BIO_F_T(0:MT,KB))
      
      BIO_ALL_T(:,:)  = BIO_ALL(:,:,I3)
      BIO_MEAN_T(:,:) = BIO_MEAN(:,:,I3)
      BIO_F_T(:,:)    = BIO_F(:,:,I3)
      
      IF(PAR) CALL AEXCHANGE(NC,MYID,NPROCS,BIO_ALL_T)
      IF(PAR) CALL AEXCHANGE(NC,MYID,NPROCS,BIO_MEAN_T)
      IF(PAR) CALL AEXCHANGE(NC,MYID,NPROCS,BIO_F_T)

      BIO_ALL(:,:,I3)  = BIO_ALL_T(:,:)  
      BIO_MEAN(:,:,I3) = BIO_MEAN_T(:,:) 
      BIO_F(:,:,I3)    = BIO_F_T(:,:)    
      
      DEALLOCATE(BIO_ALL_T,BIO_MEAN_T,BIO_F_T)
     END DO 
   RETURN
#endif
  END SUBROUTINE BIO_EXCHANGE



  SUBROUTINE BIO_PARAMETER_REPORT
    IMPLICIT NONE

    INTEGER :: I,J,K

    !*******************************************************************!
    ! *********     PRINT OUT MODEL SETUP AND PARAMETER VALUES      ****!
    !*******************************************************************!
    IF (MSR) THEN
       PRINT*
       PRINT*,'*****************************************************'
       PRINT*,'**  STRUCTURE AND FUNCTION OF THE BIOLOGICAL MODEL **'
       PRINT*,'*****************************************************'

       PRINT*
       PRINT*,'MODEL STRUCTURE        : ',  MODEL_STRUCTURE
       DO I=1,NNN
          PRINT*,'                         ',  NUTRIENT_NAME(I)
       END DO
       DO I=1,NNP
          PRINT*,'                         ',  PHYTO_NAME(I)
       END DO
       DO I=1,NNZ
          PRINT*,'                         ',  ZOO_NAME(I)
       END DO
       DO I=1,NND
          PRINT*,'                         ',  DETRITUS_NAME(I)
       END DO
       DO I=1,NNM
          PRINT*,'                         ', DOM_NAME(I)
       END DO
       DO I=1,NNB
          PRINT*,'                         ', BACTERIA_NAME(I)
       END DO
       WRITE(*,'(A26,A20)')' LIGHT FUNCTION         : ', L_FUNCTION
       WRITE(*,'(A26,A20)')' GRAZING FUNCTION       : ', G_FUNCTION
       PRINT*
       PRINT*,'*********    PHYTOPLANTON PARAMETERS    **************'
       PRINT*
       IF(L_FUNCTION.NE.'EXP_LIGHT'.AND.L_FUNCTION.NE.'SL62_LIGHT') THEN
          PRINT*,'ALPHA                  : '  , (ALPHAP(I),I=1,NNP)
       END IF
       PRINT*,'L_N COMBINE            : ' , (ALPHA_U(I),I=1,NNP)
       PRINT*,'T FORCING              : ' , (A_TP(I),I=1,NNP)
       PRINT*,'CHL ATTANUATION        : ' , ATANU_C
       PRINT*,'D ATTANUATION          : ', ATANU_D
       PRINT*,'WATER ATTANUATION      : ', ATANU_W
       IF(L_FUNCTION.EQ.'PGH80_LIGHT'.OR.L_FUNCTION.EQ.'V65_LIGHT'.OR. &
            L_FUNCTION.EQ.'BWDC9_LIGHT') THEN       
          PRINT*,'BETAP                  : ', (BETAP(I),I=1,NNP)
       END IF
       PRINT*,'CHL:C                  : ',(CHL2C(I),I=1,NNP)
       PRINT*,'ACTIVE DOM EXUD.       : ',(D_DOM(I),I=1,NNP)
       PRINT*,'PASSIVE DOM EXUD.      : ', (DPDOM(I),I=1,NNP)
       IF(L_FUNCTION.EQ.'SL62_LIGHT'.OR.L_FUNCTION.EQ.'V65_LIGHT'.OR. &
            L_FUNCTION.EQ.'PE78_LIGHT') THEN       
          PRINT*,'OPTIMAL LIGHT          : ',(I_OPT(I),I=1,NNP)
       END IF
       IF(L_FUNCTION.EQ.'MM_LIGHT'.OR.L_FUNCTION.EQ.'LB_LIGHT') THEN 
          PRINT*,'LIGHT HALF SATURATTION : ', (K_LIGHT(I),I=1,NNP)
       END IF
       PRINT*,'MORTALITY              : ', (MPD(I),I=1,NNP)
       PRINT*,'MORTALITY POWER        : ', (M_P(I),I=1,NNP)
       IF(L_FUNCTION.EQ.'LB_LIGHT'.OR.L_FUNCTION.EQ.'V65_LIGHT') THEN 
          PRINT*,'POWER OF LIGHT         : ', (N_P(I),I=1,NNP)
       END IF
       PRINT*,'THRESHOLD              : ', (P_0(I),I=1,NNP)
       PRINT*,'T ON RESPIRATION       : ', RP_T
       PRINT*,'RESPIRATION            : ', (R_P(I),I=1,NNP)
       PRINT*,'OPTIMAL T              : ', (T_OPTP(I),I=1,NNP)
       PRINT*,'MAXIMUM GROWTH         : ', (UMAX(I),I=1,NNP)
       PRINT*,'SINKING VELOCITY       : ', (W_P(I),I=1,NNP)
       PRINT*

       PRINT*,'********    ZOOPLANKTON PARAMETERS   ****************'
       PRINT*
       PRINT*,'ACTIVE RESPIRATION     : ', (ACTIVE_R(I),I=1,NNZ)
       PRINT*,'T FORCING EXPONENTIAL  : ', (A_TZ(I),I=1,NNZ)
       IF (NNB.GE.1) THEN
          PRINT*,'EFFICIENCY ON BACTERIA : ',((EFFIB(I,J),I=1,NNB),J=1,NNZ)
       END IF
       IF (NND.GE.1) THEN
          PRINT*,'EFFICIENCY ON DETRITUS : ',((EFFID(I,J),I=1,NND),J=1,NNZ)
       END IF
       PRINT*,'EFFICIENCY ON PHYTO    : ',((EFFIP(I,J),I=1,NNP),J=1,NNZ)
       PRINT*,'EFFICIENCY ON ZOO      : ',((EFFIZ(I,J),I=1,NNZ),J=1,NNZ)
       PRINT*,'MAX GRAIZING RATE      : ',(G_MAX(I),I=1,NNZ)
       IF(G_FUNCTION.EQ.'RECTI_G'.OR.G_FUNCTION.EQ.'MM1_G' .AND. &
            G_FUNCTION.EQ.'MM2_G') THEN 
          PRINT*,'HALF SATURATION        : ', (K_ZG(I),I=1,NNZ)
       END IF
       PRINT*,'GRAZING POWER          : ',(M_G(I),I=1,NNZ)
       PRINT*,'MORALITY               : ',(MZD(I),I=1,NNZ)
       PRINT*,'MORTALITY POWER        : ',(M_Z(I),I=1,NNZ)
       IF(G_FUNCTION.EQ.'MM2_G') THEN 
          PRINT*,'GRAZING THRESHOLD      : ',(P_C(I),I=1,NNZ)
       END IF
       IF (NNZ.GE.1) THEN
          PRINT*,'RECRUITMENT            : ',(R_RECRUIT(I),I=1,NNZ)
       END IF
       PRINT*,'RESPIRATION            : ',(R_Z(I),I=1,NNZ)
       IF (NNB.GE.1) THEN
          PRINT*,'PREFERENCE ON BACTERIA : ',((SIGMA_B(I,J),I=1,NNB),J=1,NNZ)
       END IF
       IF (NND.GE.1) THEN
          PRINT*,'PREFERENCE ON DETRITUS : ',((SIGMA_D(I,J),I=1,NND),J=1,NNZ)
       END IF
       PRINT*,'PREFERENCE ON PHYTO    : ',((SIGMA_P(I,J),I=1,NNP),J=1,NNZ)
       PRINT*,'PREFERENCE ON ZOO      : ',((SIGMA_Z(I,J),I=1,NNZ),J=1,NNZ)
       PRINT*,'OPTIMAL T              : ',(T_OPTZ(I),I=1,NNZ)
       PRINT*,'ZOO THRESHOLD          : ',(Z_0(I),I=1,NNZ)
       PRINT*
       PRINT*,'*********    NUTRIENT PARAMETERS    ***********'
       PRINT*
       PRINT*,'HALF-SATURATION',((KSN(I,J),I=1,NNN),J=1,NNP)
       IF (NNB.GE.1) THEN
          PRINT*,'ELEMENT RATIO IN BAC.  : ',((N2CB(I,J),I=1,NNN),J=1,NNB)
       END IF
       IF (NND.GE.1) THEN
          PRINT*,'ELEMENT RATIO IN D     : ',((N2CD(I,J),I=1,NNN),J=1,NND)
       END IF
       PRINT*,'ELEMENT RATIO IN PHYTO : ',((N2CP(I,J),I=1,NNN),J=1,NNP)
       PRINT*,'ELEMENT RATIO IN ZOO   : ',((N2CZ(I,J),I=1,NNN),J=1,NNZ)
       IF (NNM.GE.1)  THEN
          PRINT*,'ELEMENT RATIO IN DOM   : ', ((N2CDOM(I,J),I=1,NNN),J=1,NNM)
       END IF
       PRINT*,'THRESHOLD              : ', (N_0(I),I=1,NNN)
       IF (NO3_ON)    THEN
          PRINT*,'NITRIFICATION RATE     : ', R_AN
       END IF
       PRINT*
       IF (NNB.GE.1) THEN
          PRINT*,'*********  BACTERIA PARAMETERS  ************'
          PRINT*
          PRINT*,'T_FORCING              : ',(A_TB(I),I=1,NNB)
          PRINT*,'THRESHOLD              : ', (B_0(I),I=1,NNB)
          PRINT*,'RATIO OF NH4 VS DON    : ',(DELTA_B(I),I=1,NNB)
          PRINT*,'EFFICIENCY OF DETRITUS : ',((EFFIBD(I,J),I=1,NND),J=1,NNB)
          PRINT*,'EFFICIENCY OF DOM      : ',((EFFIDOM(I,J),I=1,NNM),J=1,NNB)
          PRINT*,'EFFICIENCY OF NUTRIENT : ', ((EFFIN(I,J),I=1,NNN),J=1,NNB)
          PRINT*,'RESPIRATION            : ',(R_B(I),I=1,NNB)
          PRINT*,'PREFERENCE ON DETRITUS : ',((SIGMA_BD(I,J),I=1,NND),J=1,NNB)
          PRINT*,'PREFERENCE ON DOM      : ',((SIGMA_DOM(I,J),I=1,NNM),J=1,NNB)
          PRINT*,'PREFERENCE ON NUTRIENT : ',((SIGMA_N(I,J),I=1,NNN),J=1,NNB)
          PRINT*,'OPTIMAL TEMPERATURE    : ',(T_OPTB(I),I=1,NNB)
          PRINT*,'MAXIMUM GROWTH RATE    : ',(UBMAX(I),I=1,NNB)
          PRINT*
       END IF
       IF (NND.GE.1) THEN
          PRINT*, '********  DETRITUS PARAMETERS  *************'
          PRINT*
          IF (NNB.GE.1) THEN
             PRINT*,'GRAZING ON B TO D      : ', (((ALPHA_BD(I,J,K),I=1,NND),J=1,NNB),K=1,NNZ)
          END IF
          PRINT*,'AGGREGATION            : ', (ALPHA_DAG(I),I=1,NND)
          PRINT*,'GRAZING LOSS ON D TO D : ', (((ALPHA_DD(I,J,K),I=1,NND),J=1,NND),K=1,NNZ)
          PRINT*,'DISAGGREGATION         : ', (ALPHA_DAG(I),I=1,NND)
          PRINT*,'GRAZING LOSS ON P TO D : ', (((ALPHA_PD(I,J,K),I=1,NND),J=1,NNP),K=1,NNZ)
          PRINT*,'GRAZING LOSS ON Z TO D : ', (((ALPHA_ZD(I,J,K),I=1,NND),J=1,NNZ),K=1,NNZ)
          IF (NNM.GE.1) THEN
             PRINT*,'DISSOLUTION            : ', (D_DOM(I),I=1,NND)
          END IF
          PRINT*,'REMINERALIZATION       : ', (D_RN(I),I=1,NND)
          PRINT*,'THRESHOLD              : ', (D_0(I),I=1,NND)
          PRINT*,'P MORTALITY TO DETRITUS: ', ((EPSILON_PD(I,J),I=1,NND),J=1,NNP)
          PRINT*,'Z MORTALITY TO DETRITUS: ', ((EPSILON_ZD(I,J),I=1,NND),J=1,NNZ)
          PRINT*,'SINING VELOCITY        : ', (W_D(I),I=1,NND)
          PRINT*
       END IF
       IF (NNM.GE.1) THEN
          PRINT*,'*********     DOM PARAMETERS     *************'
          PRINT*
          PRINT*,'DOM AGEING COEFFICIENT : ', (ALPHA_DOM(I),I=1,NNM)
          PRINT*,'PHYTO EXUDATION        : ', ((ALPHA_PDOM(I,J),I=1,NNM),J=1,NNP)
          PRINT*,'DETRITUS DISSOLUTION   : ', ((ALPHA_DDOM(I,J),I=1,NNM),J=1,NND)
          IF (NNB.GE.1) THEN
             PRINT*,'GRAZING LOSS ON B > DOM: ', (((ALPHA_ZBDOM(I,J,K),I=1,NNM),J=1,NNB),K=1,NNZ)
          END IF
          IF (NND.GE.1) THEN
             PRINT*,'GRAZING LOSS ON D > DOM: ', (((ALPHA_ZDDOM(I,J,K),I=1,NNM),J=1,NND),K=1,NNZ)
          END IF
          PRINT*,'GRAZING LOSS ON P > DOM: ', (((ALPHA_ZPDOM(I,J,K),I=1,NNM),J=1,NNP),K=1,NNZ)
          PRINT*,'GRAZING LOSS ON Z > DOM: ', (((ALPHA_ZZDOM(I,J,K),I=1,NNM),J=1,NNZ),K=1,NNZ)
          PRINT*,'THRESHOLD              : ', (DOM_0(I),I=1,NNM)
          PRINT*
       END IF
       PRINT*
       PRINT*,'*****************************************************'
       PRINT*,'*********    END OF BIOLOGICAL MODEL     ************'
       PRINT*,'*****************************************************'
    END IF !(MSR)

  END SUBROUTINE BIO_PARAMETER_REPORT



   SUBROUTINE BIO_INITIAL
!=============================================================================!
! THS PROGRAM INITIALIZES THE 3D BIOLOGICAL FIELD FOR THE GENERALIZED         !
! BIOLOGICAL MODEL: BIO_ALL(I,K,Nl),N1=1,NTT),AND MEAN VALUES BIO_MEAN        !
! EACH BIOLOGICAL STATE VARIABLE HAS AN INDEPENDENT INITIAL CONDITION FILE    !
! PLACED IN INPDIR. THEY SHOULD BE NAME AS "NUTRIENT_INI_1", "NUTRIENT_INI_2",!
! "PHYTOPLANKTON_INI_1", "ZOOPLANKTON_INI_1", "BACTERIA_INI_1", 'DETRITUS_    !
! INI_1", "DOM_INI_1" AND SO FORTH. THREE TYPES OF INITIAL CONDITIONS WERE    !
! IMPLEMENTED: (1) 'CONSTANT': A SINGLE VALUE; (2) 'LINEAR':WITH AT LEAST TWO !
! PAIRS OF VALUES WITH DEPTH. VARIABLE VALUES WILL BE LINEARLY INTERPOLATED   !
! BETWEEN THE VALUES GIVEN), (3) "DATA": OBSERVATION DATA SHOULD BE INTER-    !
! POLATED ONTO THE GRID POINTS AT STANDARD LEVELS. VARIABLE VALUES WILL BE    !
! INTERPOLATED AT EACH GRID POINT FROM THE DATA. THE TYPE OF INITIAL CONDI-   !
! TION SHOULD BE PUT ON THE FIRST LINE OF EACH INITIAL FILE                   !
!=============================================================================!
  USE MOD_NCTOOLS
      IMPLICIT NONE
      INTEGER :: I,J,K,LL,N_DATA
      CHARACTER(LEN=80) :: ISTR
      CHARACTER(LEN=1)  :: BIO_NUMBER
      CHARACTER(LEN=10) :: INI_TYPE
      REAL(SP), DIMENSION(KBM1)    :: ZM           !GRID DEPTH
      REAL(SP), DIMENSION(500)     :: DEPTH_STD    !STANDARD DEPTH OF DATA
      REAL(SP), DIMENSION(500)     :: DATA_BIO     !STANDARD DATA FOR LINEAR INTERPOLATION
      REAL(SP), DIMENSION(KB)      :: DATA_INT     !INTERPOLDATED VALUES

      INTEGER KSL      
      REAL(SP), ALLOCATABLE :: DPTHSL(:)
      REAL(SP), ALLOCATABLE,DIMENSION(:,:) :: TEMPB        !TEMPERAL FOR DATA INPUT
      REAL(SP), ALLOCATABLE,DIMENSION(:,:) :: DATA_3D      !3D OBSERVATION DATA

      ALLOCATE(BIO_ALL(0:MT,KB,NTT))    ; BIO_ALL     =  0.001_SP
      ALLOCATE(BIO_F(0:MT,KB,NTT))      ; BIO_F       =  0.001_SP
      ALLOCATE(BIO_MEAN(0:MT,KB,NTT))   ; BIO_MEAN    =  0.001_SP
      ALLOCATE(XFLUX_OBCB(0:MT,KB,NTT)) ; XFLUX_OBCB  =  0.0_SP
      ALLOCATE(BIO_MEANN(0:NT,KB,NTT))  ; BIO_MEANN   =  0.001_SP
      ALLOCATE(BIO_VAR_MEAN(0:MT,KB,NTT)) ; BIO_VAR_MEAN   =  0.0_SP

!*********      NUTRIENT INITIAL CONDITIONS   *************
      DO LL=1,NNN
        WRITE(BIO_NUMBER,'(I1.1)')LL
        OPEN(1,FILE=TRIM(INPUT_DIR)//'/NUTRIENT_INI_'//TRIM(BIO_NUMBER)//'.dat',STATUS='old')
        READ(1,*)INI_TYPE
        IF (TRIM(INI_TYPE).EQ.'CONSTANT') THEN
          READ(1,*)DATA_BIO(1)
          DO I=1,M
            DO K=1,KB
              BIO_ALL(I,K,LL+INN-1)=DATA_BIO(1)
            END DO
          END DO
        END IF
        IF (TRIM(INI_TYPE).EQ.'LINEAR') THEN
          N_DATA=1
 11       READ(1,*,END=12)DEPTH_STD(N_DATA),DATA_BIO(N_DATA)
          N_DATA=N_DATA+1
          GO TO 11
 12       N_DATA=N_DATA-1
          DO I=1,M
            DO K= 1,KBM1
              ZM(K) =ZZ(I,K)*(D(I)+EL(I))
            END DO
            CALL SINTER_BIO(DEPTH_STD,DATA_BIO,ZM,DATA_INT,N_DATA,KBM1)
            DO K =1,KBM1
              BIO_ALL(I,K,LL+INN-1) = DATA_INT(K)
            END DO
          END DO
        END IF !LINEAR

        IF (TRIM(INI_TYPE).EQ.'DATA') THEN
	  ALLOCATE(DPTHSL(KSL))
	  ALLOCATE(TEMPB(MGL,KSL))
	  ALLOCATE(DATA_3D(M,KSL))
	  READ(1,*) KSL
          READ(1,*) (DPTHSL(K),K=1,KSL)

          DO I=1,MGL 
            READ(1,*) (TEMPB(I,K), K=1,KSL)
          END DO

          IF(SERIAL)THEN
            DATA_3D = TEMPB
          END IF

#   if defined (MULTIPROCESSOR)
          IF(PAR)THEN
            DO I=1,M
              DATA_3D(I,1:KSL) = TEMPB(NGID(I),1:KSL)
            END DO
          END IF
#    endif
          DO I=1,M
            DO K=1,KSL
              DATA_BIO(K)=DATA_3D(I,K)
            END DO
            DO K= 1,KBM1
              ZM(K) =ZZ(I,K)*D(I)+EL(I)
            END DO
            CALL SINTER_BIO(DPTHSL,DATA_BIO,ZM,DATA_INT,KSL,KBM1)
            DO K =1,KBM1
              BIO_ALL(I,K,LL+INN-1) = DATA_INT(K)
            END DO
          END DO
	  DEALLOCATE(DPTHSL,TEMPB,DATA_3D)
        END IF !DATA
      END DO !L=1,NNN; NUTRIENT INITIALIZATION
!*********      PHYTOPLANKTON INITIAL CONDITIONS   *************
      DO LL=1,NNP
        WRITE(BIO_NUMBER,'(I1.1)')LL
        OPEN(1,FILE=TRIM(INPUT_DIR)//'/PHYTOPLANKTON_INI_'//TRIM(BIO_NUMBER)//'.dat',STATUS='old')
        READ(1,*)INI_TYPE
        IF (TRIM(INI_TYPE).EQ.'CONSTANT') THEN
          READ(1,*)DATA_BIO(1)
          DO I=1,M
            DO K=1,KB
              BIO_ALL(I,K,LL+INP-1)=DATA_BIO(1)
            END DO
          END DO
        END IF
        IF (TRIM(INI_TYPE).EQ.'LINEAR') THEN
          N_DATA=1
 21       READ(1,*,END=22)DEPTH_STD(N_DATA),DATA_BIO(N_DATA)
          N_DATA=N_DATA+1
          GO TO 21
 22       N_DATA=N_DATA-1
          DO I=1,M
            DO K= 1,KBM1
              ZM(K) =ZZ(I,K)*D(I)+EL(I)
            END DO
            CALL SINTER_BIO(DEPTH_STD,DATA_BIO,ZM,DATA_INT,N_DATA,KBM1)
            DO K =1,KBM1
              BIO_ALL(I,K,LL+INP-1) = DATA_INT(K)
            END DO
          END DO
        END IF !LINEAR

        IF (TRIM(INI_TYPE).EQ.'DATA') THEN
	  ALLOCATE(DPTHSL(KSL))
	  ALLOCATE(TEMPB(MGL,KSL))
	  ALLOCATE(DATA_3D(M,KSL))
          DO I=1,MGL 
            READ(1,*) (TEMPB(I,K), K=1,KSL)
          END DO

          IF(SERIAL)THEN
            DATA_3D = TEMPB
          END IF

#   if defined (MULTIPROCESSOR)
          IF(PAR)THEN
            DO I=1,M
              DATA_3D(I,1:KSL) = TEMPB(NGID(I),1:KSL)
            END DO
          END IF
#    endif
          DO I=1,M
            DO K=1,KSL
              DATA_BIO(K)=DATA_3D(I,K)
            END DO
            DO K= 1,KBM1
              ZM(K) =ZZ(I,K)*D(I)+EL(I)
            END DO
            CALL SINTER_BIO(DPTHSL,DATA_BIO,ZM,DATA_INT,KSL,KBM1)
            DO K =1,KBM1
              BIO_ALL(I,K,LL+INP-1) = DATA_INT(K)
            END DO
          END DO
	  DEALLOCATE(DPTHSL,TEMPB,DATA_3D)
        END IF !DATA
      END DO !L=1,NNP; PHYTOPLANKTON INITIALIZATION

!*********      ZOOPLANKTON INITIAL CONDITIONS   *************
      DO LL=1,NNZ
        WRITE(BIO_NUMBER,'(I1.1)')LL
        OPEN(1,FILE=TRIM(INPUT_DIR)//'/ZOOPLANKTON_INI_'//TRIM(BIO_NUMBER)//'.dat',STATUS='old')
        READ(1,*)INI_TYPE
        IF (TRIM(INI_TYPE).EQ.'CONSTANT') THEN
          READ(1,*)DATA_BIO(1)
          DO I=1,M
            DO K=1,KB
              BIO_ALL(I,K,LL+INZ-1)=DATA_BIO(1)
            END DO
          END DO
        END IF
        IF (TRIM(INI_TYPE).EQ.'LINEAR') THEN
          N_DATA=1
 31       READ(1,*,END=32)DEPTH_STD(N_DATA),DATA_BIO(N_DATA)
          N_DATA=N_DATA+1
          GO TO 31
 32       N_DATA=N_DATA-1
          DO I=1,M
            DO K= 1,KBM1
              ZM(K) =ZZ(I,K)*D(I)+EL(I)
            END DO
            CALL SINTER_BIO(DEPTH_STD,DATA_BIO,ZM,DATA_INT,N_DATA,KBM1)
            DO K =1,KBM1
              BIO_ALL(I,K,LL+INZ-1) = DATA_INT(K)
            END DO
          END DO
        END IF !LINEAR

        IF (TRIM(INI_TYPE).EQ.'DATA') THEN
	  ALLOCATE(DPTHSL(KSL))
	  ALLOCATE(TEMPB(MGL,KSL))
	  ALLOCATE(DATA_3D(M,KSL))
          DO I=1,MGL 
            READ(1,*) (TEMPB(I,K), K=1,KSL)
          END DO

          IF(SERIAL)THEN
            DATA_3D = TEMPB
          END IF

#   if defined (MULTIPROCESSOR)
          IF(PAR)THEN
            DO I=1,M
              DATA_3D(I,1:KSL) = TEMPB(NGID(I),1:KSL)
            END DO
          END IF
#    endif
          DO I=1,M
            DO K=1,KSL
              DATA_BIO(K)=DATA_3D(I,K)
            END DO
            DO K= 1,KBM1
              ZM(K) =ZZ(I,K)*D(I)+EL(I)
            END DO
            CALL SINTER_BIO(DPTHSL,DATA_BIO,ZM,DATA_INT,KSL,KBM1)
            DO K =1,KBM1
              BIO_ALL(I,K,LL+INZ-1) = DATA_INT(K)
            END DO
          END DO
	  DEALLOCATE(DPTHSL,TEMPB,DATA_3D)
        END IF !DATA
      END DO !L=1,NNZ; ZOOPLANKTON INITIALIZATION


!*********      BACTERIA INITIAL CONDITIONS   *************
      DO LL=1,NNB
        WRITE(BIO_NUMBER,'(I1.1)')LL
        OPEN(1,FILE=TRIM(INPUT_DIR)//'/BACTERIA_INI_'//TRIM(BIO_NUMBER)//'.dat',STATUS='old')
        READ(1,*)INI_TYPE
        IF (TRIM(INI_TYPE).EQ.'CONSTANT') THEN
          READ(1,*)DATA_BIO(1)
          DO I=1,M
            DO K=1,KB
              BIO_ALL(I,K,LL+INB-1)=DATA_BIO(1)
            END DO
          END DO
        END IF
        IF (TRIM(INI_TYPE).EQ.'LINEAR') THEN
          N_DATA=1
 41       READ(1,*,END=42)DEPTH_STD(N_DATA),DATA_BIO(N_DATA)
          N_DATA=N_DATA+1
          GO TO 41
 42       N_DATA=N_DATA-1
          DO I=1,M
            DO K= 1,KBM1
              ZM(K) =ZZ(I,K)*D(I)+EL(I)
            END DO
            CALL SINTER_BIO(DEPTH_STD,DATA_BIO,ZM,DATA_INT,N_DATA,KBM1)
            DO K =1,KBM1
              BIO_ALL(I,K,LL+INB-1) = DATA_INT(K)
            END DO
          END DO
        END IF !LINEAR

        IF (TRIM(INI_TYPE).EQ.'DATA') THEN
	  ALLOCATE(DPTHSL(KSL))
	  ALLOCATE(TEMPB(MGL,KSL))
	  ALLOCATE(DATA_3D(M,KSL))
          DO I=1,MGL 
            READ(1,*) (TEMPB(I,K), K=1,KSL)
          END DO

          IF(SERIAL)THEN
            DATA_3D = TEMPB
          END IF

#   if defined (MULTIPROCESSOR)
          IF(PAR)THEN
            DO I=1,M
              DATA_3D(I,1:KSL) = TEMPB(NGID(I),1:KSL)
            END DO
          END IF
#    endif
          DO I=1,M
            DO K=1,KSL
              DATA_BIO(K)=DATA_3D(I,K)
            END DO
            DO K= 1,KBM1
              ZM(K) =ZZ(I,K)*D(I)+EL(I)
            END DO
            CALL SINTER_BIO(DPTHSL,DATA_BIO,ZM,DATA_INT,KSL,KBM1)
            DO K =1,KBM1
              BIO_ALL(I,K,LL+INB-1) = DATA_INT(K)
            END DO
          END DO
	  DEALLOCATE(DPTHSL,TEMPB,DATA_3D)
        END IF !DATA
      END DO !L=1,NNB; BACTERIA INITIALIZATION


!*********      DOM INITIAL CONDITIONS   *************
      DO LL=1,NNM
        WRITE(BIO_NUMBER,'(I1.1)')LL
        OPEN(1,FILE=TRIM(INPUT_DIR)//'/DOM_INI_'//TRIM(BIO_NUMBER)//'.dat',STATUS='old')
        READ(1,*)INI_TYPE
        IF (TRIM(INI_TYPE).EQ.'CONSTANT') THEN
          READ(1,*)DATA_BIO(1)
          DO I=1,M
            DO K=1,KB
              BIO_ALL(I,K,LL+INM-1)=DATA_BIO(1)
            END DO
          END DO
        END IF
        IF (TRIM(INI_TYPE).EQ.'LINEAR') THEN
          N_DATA=1
 51       READ(1,*,END=52)DEPTH_STD(N_DATA),DATA_BIO(N_DATA)
          N_DATA=N_DATA+1
          GO TO 51
 52       N_DATA=N_DATA-1
          DO I=1,M
            DO K= 1,KBM1
              ZM(K) =ZZ(I,K)*D(I)+EL(I)
            END DO
            CALL SINTER_BIO(DEPTH_STD,DATA_BIO,ZM,DATA_INT,N_DATA,KBM1)
            DO K =1,KBM1
              BIO_ALL(I,K,LL+INM-1) = DATA_INT(K)
            END DO
          END DO
        END IF !LINEAR

        IF (TRIM(INI_TYPE).EQ.'DATA') THEN
	  ALLOCATE(DPTHSL(KSL))
	  ALLOCATE(TEMPB(MGL,KSL))
	  ALLOCATE(DATA_3D(M,KSL))
          DO I=1,MGL 
            READ(1,*) (TEMPB(I,K), K=1,KSL)
          END DO

          IF(SERIAL)THEN
            DATA_3D = TEMPB
          END IF

#   if defined (MULTIPROCESSOR)
          IF(PAR)THEN
            DO I=1,M
              DATA_3D(I,1:KSL) = TEMPB(NGID(I),1:KSL)
            END DO
          END IF
#    endif
          DO I=1,M
            DO K=1,KSL
              DATA_BIO(K)=DATA_3D(I,K)
            END DO
            DO K= 1,KBM1
              ZM(K) =ZZ(I,K)*D(I)+EL(I)
            END DO
            CALL SINTER_BIO(DPTHSL,DATA_BIO,ZM,DATA_INT,KSL,KBM1)
            DO K =1,KBM1
              BIO_ALL(I,K,LL+INM-1) = DATA_INT(K)
            END DO
          END DO
	  DEALLOCATE(DPTHSL,TEMPB,DATA_3D)
        END IF !DATA
      END DO !L=1,NNM; DOM INITIALIZATION

!*********      DETRITUS INITIAL CONDITIONS   *************
      DO LL=1,NND
        WRITE(BIO_NUMBER,'(I1.1)')LL
        OPEN(1,FILE=TRIM(INPUT_DIR)//'/DETRITUS_INI_'//TRIM(BIO_NUMBER)//'.dat',STATUS='old')
        READ(1,*)INI_TYPE
        IF (TRIM(INI_TYPE).EQ.'CONSTANT') THEN
          READ(1,*)DATA_BIO(1)
          DO I=1,M
            DO K=1,KB
              BIO_ALL(I,K,LL+IND-1)=DATA_BIO(1)
            END DO
          END DO
        END IF
        IF (TRIM(INI_TYPE).EQ.'LINEAR') THEN
          N_DATA=1
 61       READ(1,*,END=62)DEPTH_STD(N_DATA),DATA_BIO(N_DATA)
          N_DATA=N_DATA+1
          GO TO 61
 62       N_DATA=N_DATA-1
          DO I=1,M
            DO K= 1,KBM1
              ZM(K) =ZZ(I,K)*D(I)+EL(I)
            END DO
            CALL SINTER_BIO(DEPTH_STD,DATA_BIO,ZM,DATA_INT,N_DATA,KBM1)
            DO K =1,KBM1
              BIO_ALL(I,K,LL+IND-1) = DATA_INT(K)
            END DO
          END DO
        END IF !LINEAR

        IF (TRIM(INI_TYPE).EQ.'DATA') THEN
	  ALLOCATE(DPTHSL(KSL))
	  ALLOCATE(TEMPB(MGL,KSL))
	  ALLOCATE(DATA_3D(M,KSL))
          DO I=1,MGL 
            READ(1,*) (TEMPB(I,K), K=1,KSL)
          END DO

          IF(SERIAL)THEN
            DATA_3D = TEMPB
          END IF

#   if defined (MULTIPROCESSOR)
          IF(PAR)THEN
            DO I=1,M
              DATA_3D(I,1:KSL) = TEMPB(NGID(I),1:KSL)
            END DO
          END IF
#    endif
          DO I=1,M
            DO K=1,KSL
              DATA_BIO(K)=DATA_3D(I,K)
            END DO
            DO K= 1,KBM1
              ZM(K) =ZZ(I,K)*D(I)+EL(I)
            END DO
            CALL SINTER_BIO(DPTHSL,DATA_BIO,ZM,DATA_INT,KSL,KBM1)
            DO K =1,KBM1
              BIO_ALL(I,K,LL+IND-1) = DATA_INT(K)
            END DO
          END DO
	  DEALLOCATE(DPTHSL,TEMPB,DATA_3D)
        END IF !DATA
      END DO !L=1,NND; DETRITUS INITIALIZATION
      WHERE (BIO_ALL < 0.001) BIO_ALL=0.001
      BIO_MEAN=BIO_ALL   !3D ASSIGNMENT
#    if defined (MULTIPROCESSOR)
      IF (PAR) CALL BIO_EXCHANGE
#    endif
!!$      CALL BIO_NETCDF_HEADER
   RETURN
END SUBROUTINE BIO_INITIAL


SUBROUTINE BIO_HOT_START
! THIS SUBROUTINE READS IN RESTART BIOLOGICAL DATA FROM THE NETCDF FILE restart_bio.nc  !
! IT QUERIES THE DIMENSION AND TAKES THE LAST TIME STEP AS THE RESTART DATA.            !

      USE NETCDF  

      IMPLICIT NONE
      SAVE
      INTEGER  ::  I,J,K,IERR,NC_FID,N_START,VARID,DIMID,DIMS(3)
      REAL(SP) ::  TEMPB(MGL,KBM1)
      CHARACTER(LEN=20) ::  TEMPNAME,time
      ALLOCATE(BIO_ALL(0:MT,KB,NTT))    ; BIO_ALL     =  0.001_SP
      ALLOCATE(BIO_F(0:MT,KB,NTT))      ; BIO_F       =  0.001_SP
      ALLOCATE(BIO_MEAN(00:MT,KB,NTT))  ; BIO_MEAN    =  0.001_SP
      ALLOCATE(XFLUX_OBCB(0:MT,KB,NTT)) ; XFLUX_OBCB  =  0.0_SP
      ALLOCATE(BIO_MEANN(0:NT,KB,NTT))  ; BIO_MEANN   =  0.001_SP
      ALLOCATE(BIO_VAR_MEAN(0:MT,KB,NTT)) ; BIO_VAR_MEAN   =  0.0_SP

!*******************   START EXECUTABLE      *******************!

!      IF(MSR) THEN
        IERR = NF90_OPEN('restart_bio.nc',NF90_NOWRITE,NC_FID)
        IF(IERR /=NF90_NOERR)THEN
           WRITE(*,*)' ERROR IN OPENNING restart_bio.nc '
           STOP
        END IF
!      END IF
      IERR = NF90_INQ_DIMID(NC_FID,'time',DIMID)   !GET time DIMENSION ID
      IF(IERR /=NF90_NOERR)THEN
         WRITE(*,*)' ERROR GETTING time ID: '
         STOP
      END IF
      IERR = NF90_INQUIRE_DIMENSION(NC_FID,DIMID,TEMPNAME,N_START)
      IF(IERR /=NF90_NOERR)THEN
         WRITE(*,*)' ERROR GETTING TIME STEPS '
         STOP
      END IF
!**************   DETERMINE START POINT
      DIMS(1) = 1
      DIMS(2) = 1
      DIMS(3) = N_START
      DO I=1,NTT
         IERR = nf90_inq_varid(NC_FID,TRIM(BIO_NAME(I,1)),VARID)
         IF(IERR /=NF90_NOERR)THEN
           WRITE(*,*)'ERROR GETTING THE ID OF ',TRIM(BIO_NAME(I,1))
           STOP
         END IF

         IERR = nf90_get_var(NC_FID,VARID,TEMPB,START=DIMS)
         IF(IERR /=NF90_NOERR)THEN
            WRITE(*,*)'ERROR GETTING RESTART DATA OF ',TRIM(BIO_NAME(I,1))
            STOP
         END IF

      IF(SERIAL)THEN
        BIO_ALL(:,:,I)=TEMPB(:,:)
      END IF
#   if defined (MULTIPROCESSOR)
          IF(PAR)THEN
            DO J=1,M
              BIO_ALL(J,1:KBM1,I)=TEMPB(NGID(J),1:KBM1)
            END DO
          END IF
#    endif
      END DO  !DO I=1,NNT
      BIO_MEAN=BIO_ALL
#    if defined (MULTIPROCESSOR)
      IF (PAR) CALL BIO_EXCHANGE
#    endif
!      CALL BIO_NETCDF_HEADER
      IERR=NF90_CLOSE(NC_FID)

END SUBROUTINE BIO_HOT_START

!==============================================================================|

   SUBROUTINE SINTER_BIO(X,A,Y,B,M1,N1)              

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
   END SUBROUTINE SINTER_BIO
!==============================================================================|
!
!=============================================================================!
!!$   SUBROUTINE NAME_LIST_INITIALIZE_BIO
!!$   USE CONTROL
   
!!$   IMPLICIT NONE
   
   !--Parameters in NameList NML_BIOLOGICAL
!!$   STARTUP_BIO_TYPE = "'constant' 'linear' 'observed' or 'set values'"
   
!!$   RETURN
!!$   END SUBROUTINE NAME_LIST_INITIALIZE_BIO
!=============================================================================!
!
!=============================================================================!   
!!$   SUBROUTINE NAME_LIST_PRINT_BIO
!!$   USE CONTROL
   
!!$   IMPLICIT NONE
   
!!$   write(UNIT=IPT,NML=NML_BIOLOGICAL)
   
!!$   RETURN
!!$   END SUBROUTINE NAME_LIST_PRINT_BIO
!=============================================================================!   
   

#    endif
!    END IF DEFINED BioGen AT THE BEGINNING
END MODULE MOD_BIO_3D
