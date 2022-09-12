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

!! NEED TO ADD NORTH POLE ROUTINE FOR DYE!!!
#  if defined (DYE_RELEASE)
MODULE MOD_DYE

   USE MOD_PREC
!   USE MOD_INP
   USE CONTROL
   USE MOD_PAR
   USE MOD_SET_TIME
   IMPLICIT NONE
   SAVE
!
!--VARIABLES for SPECIFY DYE RELEASE                 
     INTEGER, PARAMETER :: KSPE_DYE_MAX = 100
     INTEGER, PARAMETER :: MSPE_DYE_MAX = 200
     LOGICAL  :: DYE_ON                        !!RELEASE DYE ACTIVE
     CHARACTER(LEN=80) :: DYE_RELEASE_START
     CHARACTER(LEN=80) :: DYE_RELEASE_STOP

     TYPE(TIME) :: DYESTART
     TYPE(TIME) :: DYESTOP


     INTEGER  :: KSPE_DYE                      !!NUMBER OF SIGMA LAYER FOR SPECIFY DYE RELEASE 
     INTEGER  :: MSPE_DYE                      !!NUMBER OF NODE FOR SPECIFY DYE  RELEASE
     INTEGER  :: K_SPECIFY(KSPE_DYE_MAX)       !!NO of sigma layer for specify dye release
     INTEGER  :: M_SPECIFY(MSPE_DYE_MAX)       !!NO of node for specify dye release
     REAL(SP) :: DYE_SOURCE_TERM               !!Specify source term value of dye releasing
!--VARIABLES OF DYE     
     REAL(SP), ALLOCATABLE, TARGET :: DYE(:,:)       !!DYE CONCENTRATION AT NODE
     REAL(SP), ALLOCATABLE :: DYEF(:,:)      !!DYE CONCENTRATION FROM PREVIOUS TIME
!     REAL(SP), ALLOCATABLE :: DYEMEAN(:,:)   !!MEAN INITIAL DYE - not used!

#  if defined (MPDATA)
     REAL(SP), ALLOCATABLE :: DYE_S(:,:)    !! temporary dye in modified upwind
     REAL(SP), ALLOCATABLE :: DYE_SF(:,:)    !! temporary dye in modified upwind

     REAL(SP), ALLOCATABLE :: WWWS(:,:)    !! temporary dye in modified upwind
     REAL(SP), ALLOCATABLE :: WWWSF(:,:)    !! temporary dye in modified upwind
     REAL(SP), ALLOCATABLE :: DTWWWS(:)    !! temporary dye in modified upwind
     REAL(SP), ALLOCATABLE :: zzzflux(:,:)    !! temporary total flux in corrected part
     REAL(SP), ALLOCATABLE :: beta(:,:)    !! temporary beta coefficient in corrected part

     REAL(SP), ALLOCATABLE :: betain(:,:)    !! temporary beta coefficient in corrected part
     REAL(SP), ALLOCATABLE :: betaout(:,:)    !! temporary beta coefficient in corrected part
# endif

     NAMELIST /NML_DYE_RELEASE/     &
          & DYE_ON,                 &
          & DYE_RELEASE_START,      &
          & DYE_RELEASE_STOP,       &
	  & KSPE_DYE,               &
	  & MSPE_DYE,               &
	  & K_SPECIFY,              &
	  & M_SPECIFY,              &
	  & DYE_SOURCE_TERM

    integer itera,ntera
     real ssss
   CONTAINS !------------------------------------------------------------------!
            ! ALLOC_VARS_DYE  : Allocate and Initialize Arrays of dye          !
            ! NAME_LIST_READ_DYE : specify dye soerce term and parameter      !
            ! ADV_DYE       : Horizontal Advection/Diffusion of dye Variables  !
            ! VDIF_DYE      : Vertical Diffusion of dye Variables              !
            ! INITIAL_DYE   : Initialize for dye Variables - not used!         !
            ! BCOND_DYE     : Boundary Conditions (River Flux) of DYE Variables!

!==============================================================================|
!    Allocate and Initialize Arrays of Dye                                     !
!==============================================================================|

   SUBROUTINE ALLOC_VARS_DYE

!==============================================================================!
   USE ALL_VARS
   IMPLICIT NONE
   INTEGER NCT
   INTEGER NDB, status
   INTEGER(ITIME) :: tstep
   CHARACTER(LEN=4) :: BFLAG,EFLAG
# if !defined (DOUBLE_PRECISION)
   NDB = 1       !!GWC BASE THIS ON KIND
# else
   NDB = 2
# endif
   
   NCT = NT*3
!==============================================================================!
!  ALLOCATE:                                                                   !
!==============================================================================!
   ALLOCATE(DYE(0:MT,KB))           ;DYE     = ZERO
   ALLOCATE(DYEF(0:MT,KB))          ;DYEF    = ZERO
!   ALLOCATE(DYEMEAN(0:MT,KB))       ;DYEMEAN    = ZERO
#  if defined (MPDATA)
   ALLOCATE(DYE_S(0:MT,KB))         ;DYE_S   = ZERO
   ALLOCATE(DYE_SF(0:MT,KB))        ;DYE_SF  = ZERO
   ALLOCATE(WWWS(0:MT,KB))          ;WWWS    = ZERO
   ALLOCATE(WWWSF(0:MT,KB))         ;WWWSF   = ZERO
   ALLOCATE(DTWWWS(0:MT))           ;DTWWWS  = ZERO

   ALLOCATE(zzzflux(0:MT,KB))       ;zzzflux = ZERO 
   ALLOCATE(beta(0:MT,KB))          ;beta    = ZERO  
   ALLOCATE(betain(0:MT,KB))        ;betain  = ZERO
   ALLOCATE(betaout(0:MT,KB))       ;betaout = ZERO
# endif

!   add convert time
    SELECT CASE(USE_REAL_WORLD_TIME)
      CASE(.TRUE.)
          
        DYESTART = READ_DATETIME(DYE_RELEASE_START,DATE_FORMAT,TIMEZONE,status)
        if (status == 0 ) &
           & Call Fatal_Error("Could not read the date string START_DATE: "//trim(DYE_RELEASE_START))

        DYESTOP = READ_DATETIME(DYE_RELEASE_STOP,DATE_FORMAT,TIMEZONE,status)
        if (status == 0) &
           & Call Fatal_Error("Could not read the date string END_DATE: "//trim(DYE_RELEASE_STOP))

        if(DYESTART .GT. DYESTOP) &
           & Call Fatal_Error("Runfile Start_Date exceeds or equal to Stop_Date")
          
      CASE(.FALSE.) ! THIS MODEL IS USING IDEALIZED TIME

        ! GET THE START AND END INFORMATION
        CALL IDEAL_TIME_STRING2TIME(DYE_RELEASE_START,BFLAG,DYESTART,tstep)
        CALL IDEAL_TIME_STRING2TIME(DYE_RELEASE_STOP,EFLAG,DYESTOP,tstep)

        ! SANITY CHECK
        IF (BFLAG /= EFLAG) CALL FATAL_ERROR&
           ('IDEALIZED MODEL TIME SPECIFICATION IS INCORRENT',&
           &'BEGIN AND END CAN BE IN EITHER CYCLES OR TIME BUT NOT MIXED',&
           & trim(dye_release_start),trim(dye_release_stop) )

    END SELECT
!   end add convert time


   RETURN
   END SUBROUTINE ALLOC_VARS_DYE


!==============================================================================|
!   Specify source term                |
!==============================================================================|
   SUBROUTINE NAME_LIST_READ_DYE
   USE CONTROL

   IMPLICIT NONE
   integer :: ios, i
   Character(Len=120):: FNAME


   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "Subroutine Begins: set_dye_param;"

    ios = 0

    FNAME = "./"//trim(casename)//"_run.nml"

    if(DBG_SET(dbg_io)) &
         & write(IPT,*) "Set_dye_param: File: ",trim(FNAME)

    CALL FOPEN(NMLUNIT,trim(FNAME),'cfr')

    !READ NAME LIST FILE

    ! Read Model Dye Release Settings
    READ(UNIT=NMLUNIT, NML=NML_DYE_RELEASE,IOSTAT=ios)
    if(ios .NE. 0 ) then
       if(DBG_SET(dbg_log)) write(UNIT=IPT,NML=NML_DYE_RELEASE)
       Call Fatal_Error('Can Not Read NameList NML_DYE_RELEASE from file:' //trim(FNAME))
    end if


    REWIND(NMLUNIT)

    if(DBG_SET(dbg_scl)) &
         & write(IPT,*) "Read_Name_List:"

    if(DBG_SET(dbg_scl)) &
         & write(UNIT=IPT,NML=NML_DYE_RELEASE)


    IF(KSPE_DYE > KSPE_DYE_MAX)THEN
      CALL FATAL_ERROR("KSPE_DYE > KSPE_DYE_MAX")
    END IF
    IF(MSPE_DYE > MSPE_DYE_MAX) THEN
      CALL FATAL_ERROR("MSPE_DYE > MSPE_DYE_MAX")
    END IF
      
    CLOSE(NMLUNIT)


!==============================================================================|
!            SCREEN REPORT OF SET DYE RELEASE VARIABlES                        !
!==============================================================================|
!!$   IF(MSR) THEN  
!!$     WRITE(IPT,*) '!                                                   !'     
!!$     WRITE(IPT,*) '!------SPECIFY DYE RELEASE VARIABlES----------------!'     
!!$     WRITE(IPT,*) '!                                                   !'     
!!$     WRITE(IPT,*) '!  # DYE_ON              :',DYE_ON
!!$     WRITE(IPT,*) '!  # DYE_RELEASE_START   :',DYE_RELEASE_START
!!$     WRITE(IPT,*) '!  # DYE_RELEASE_STOP    :',DYE_RELEASE_STOP
!!$     WRITE(IPT,*) '!  # KSPE_DYE            :',KSPE_DYE
!!$     WRITE(IPT,*) '!  # K_SPECIFY           :',K_SPECIFY
!!$     WRITE(IPT,*) '!  # MSPE_DYE            :',MSPE_DYE
!!$     WRITE(IPT,*) '!  # M_SPECIFY           :',M_SPECIFY
!!$   END IF
!!$   RETURN
   END SUBROUTINE NAME_LIST_READ_DYE


!==============================================================================|
!   Calculate Advection and Horizontal Diffusion Terms for DYE                 |
!==============================================================================|

   SUBROUTINE ADV_DYE               

!------------------------------------------------------------------------------|

   USE ALL_VARS
   USE MOD_UTILS
   USE BCS
   USE MOD_OBCS
   USE MOD_PAR
   USE MOD_WD
   USE MOD_SPHERICAL
   USE MOD_NORTHPOLE
#  if defined (SEMI_IMPLICIT)
   USE MOD_SEMI_IMPLICIT
#  endif
#  if defined (TVD)
   USE MOD_TVD
#  endif

   IMPLICIT NONE
   REAL(SP), DIMENSION(0:MT,KB)     :: XFLUX,XFLUX_ADV
   REAL(SP), DIMENSION(M)           :: PUPX,PUPY,PVPX,PVPY  
   REAL(SP), DIMENSION(M)           :: PDYEPX,PDYEPY,PDYEPXD,PDYEPYD,VISCOFF
   REAL(SP), DIMENSION(3*(NT),KBM1) :: DTIJ 
   REAL(SP), DIMENSION(3*(NT),KBM1) :: UVN
   REAL(SP) :: UTMP,VTMP,SITAI,FFD,FF1 !,X11,Y11,X22,Y22,X33,Y33,TMP1,TMP2,XI,YI
   REAL(SP) :: DXA,DYA,DXB,DYB,FIJ1,FIJ2,UN
   REAL(SP) :: TXX,TYY,FXX,FYY,VISCOF,EXFLUX,TEMP,STPOINT
   REAL(SP) :: FACT,FM1
   REAL(SP) :: s1min, s1max, s2min, s2max,SMIN,SMAX,xxxx
   INTEGER  :: I,I1,I2,IA,IB,J,J1,J2,K,JTMP,JJ,KK,IP,II
!!$# if defined (SPHERICAL)
!!$   REAL(SP) :: ty,txpi,typi
!!$   REAL(DP) :: XTMP1,XTMP
!!$   REAL(DP) :: X1_DP,Y1_DP,X2_DP,Y2_DP,XII,YII
!!$# endif
#  if defined (WET_DRY)
   INTEGER :: N_NTVE
#  endif
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
   REAL(SP), DIMENSION(0:MT,KB)     :: OFFS     !! Offset to run the simulation in single precision
#  endif

  IF(DBG_SET(DBG_SBR)) WRITE(IPT,*)"Start: adv_dye"
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
  ! Adding offset to force FVCOM not to crash when dyes approach zero
      OFFS = 1.0_SP
      DYE  = DYE+OFFS
# endif

!
!--Initialize Fluxes and dyef----------------------------------------------------------!
!
   XFLUX     = 0.0_SP
   XFLUX_ADV = 0.0_SP
   DYEF      = 0.0_SP
!
!--Loop Over Control Volume Sub-Edges And Calculate Normal Velocity------------!
!
!!#  if !defined (WET_DRY)
   DO I=1,NCV
     I1=NTRG(I)
!     DTIJ(I)=DT1(I1)
     DO K=1,KBM1
       DTIJ(I,K)=DT1(I1)*DZ1(I1,K)
       UVN(I,K) = V(I1,K)*DLTXE(I) - U(I1,K)*DLTYE(I)
#      if defined (SEMI_IMPLICIT)
       DTIJ1(I,K) = D1(I1)*DZ1(I1,K)
       UVN1(I,K) = VF(I1,K)*DLTXE(I) - UF(I1,K)*DLTYE(I)
#      endif
     END DO
   END DO
!!#  else
!!   DO I=1,NCV
!!     I1=NTRG(I)
!!!     DTIJ(I)=DT1(I1)
!!     DO K=1,KBM1
!!       DTIJ(I,K)=DT1(I1)*DZ1(I1,K)
!!       UVN(I,K) = VS(I1,K)*DLTXE(I) - US(I1,K)*DLTYE(I)
!!#      if defined (SEMI_IMPLICIT)
!!       DTIJ1(I,K) = D1(I1)*DZ1(I1,K)
!!       UVN1(I,K) = VF(I1,K)*DLTXE(I) - UF(I1,K)*DLTYE(I)
!!#      endif
!!     END DO
!!   END DO
!!#  endif

!
!--Calculate the Advection and Horizontal Diffusion Terms----------------------!
!

   DO K=1,KBM1
      PDYEPX  = 0.0_SP 
      PDYEPY  = 0.0_SP 
      PDYEPXD = 0.0_SP 
      PDYEPYD = 0.0_SP
     DO I=1,M
       DO J=1,NTSN(I)-1
         I1=NBSN(I,J)
         I2=NBSN(I,J+1)

#    if defined (WET_DRY)
         IF(ISWETN(I1) == 0 .AND. ISWETN(I2) == 1)THEN
          FFD=0.5_SP*(DYE(I,K)+DYE(I2,K))!-DYEMEAN(I,K)-DYEMEAN(I2,K))
          FF1=0.5_SP*(DYE(I,K)+DYE(I2,K))
	 ELSE IF(ISWETN(I1) == 1 .AND. ISWETN(I2) == 0)THEN
          FFD=0.5_SP*(DYE(I1,K)+DYE(I,K))!-DYEMEAN(I1,K)-DYEMEAN(I,K))
          FF1=0.5_SP*(DYE(I1,K)+DYE(I,K))
	 ELSE IF(ISWETN(I1) == 0 .AND. ISWETN(I2) == 0)THEN
          FFD=0.5_SP*(DYE(I,K)+DYE(I,K))!-DYEMEAN(I,K)-DYEMEAN(I,K))
          FF1=0.5_SP*(DYE(I,K)+DYE(I,K))
	 ELSE
          FFD=0.5_SP*(DYE(I1,K)+DYE(I2,K))!-DYEMEAN(I1,K)-DYEMEAN(I2,K))
          FF1=0.5_SP*(DYE(I1,K)+DYE(I2,K))
	 END IF 
#    else	 
         FFD=0.5_SP*(DYE(I1,K)+DYE(I2,K))!-DYEMEAN(I1,K)-DYEMEAN(I2,K))
         FF1=0.5_SP*(DYE(I1,K)+DYE(I2,K))
#    endif

!!$#        if defined (SPHERICAL)
!!$         XTMP  = VX(I2)*TPI-VX(I1)*TPI
!!$	 XTMP1 = VX(I2)-VX(I1)
!!$	 IF(XTMP1 >  180.0_SP)THEN
!!$	   XTMP = -360.0_SP*TPI+XTMP
!!$	 ELSE IF(XTMP1 < -180.0_SP)THEN
!!$	   XTMP =  360.0_SP*TPI+XTMP
!!$	 END IF  
!!$         TXPI=XTMP*COS(DEG2RAD*VY(I))
!!$         TYPI=(VY(I1)-VY(I2))*TPI
!!$
!!$         PDYEPX(I)=PDYEPX(I)+FF1*TYPI
!!$         PDYEPY(I)=PDYEPY(I)+FF1*TXPI
!!$         PDYEPXD(I)=PDYEPXD(I)+FFD*TYPI
!!$         PDYEPYD(I)=PDYEPYD(I)+FFD*TXPI
!!$#        else
!!$         PDYEPX(I)=PDYEPX(I)+FF1*(VY(I1)-VY(I2))
!!$         PDYEPY(I)=PDYEPY(I)+FF1*(VX(I2)-VX(I1))
!!$         PDYEPXD(I)=PDYEPXD(I)+FFD*(VY(I1)-VY(I2))
!!$         PDYEPYD(I)=PDYEPYD(I)+FFD*(VX(I2)-VX(I1))
!!$#        endif

         PDYEPX(I)=PDYEPX(I)+FF1*DLTYTRIE(i,j)
         PDYEPY(I)=PDYEPY(I)+FF1*DLTXTRIE(i,j)
         PDYEPXD(I)=PDYEPXD(I)+FFD*DLTYTRIE(i,j)
         PDYEPYD(I)=PDYEPYD(I)+FFD*DLTXTRIE(i,j)

       END DO
       PDYEPX(I)=PDYEPX(I)/ART2(I)
       PDYEPY(I)=PDYEPY(I)/ART2(I)
       PDYEPXD(I)=PDYEPXD(I)/ART2(I)
       PDYEPYD(I)=PDYEPYD(I)/ART2(I)
     END DO
          
     IF(K == KBM1)THEN
       DO I=1,M
         PFPXB(I) = PDYEPX(I)
         PFPYB(I) = PDYEPY(I)
       END DO
     END IF

     DO I=1,M
       VISCOFF(I)=VISCOFH(I,K)
     END DO

     IF(K == KBM1) THEN
!       AH_BOTTOM(1:M) = HORCON*(FACT*VISCOFF(1:M) + FM1)
       AH_BOTTOM(1:M) = (FACT*VISCOFF(1:M) + FM1)*NN_HVC(1:M)
    END IF


     DO I=1,NCV_I
       IA=NIEC(I,1)
       IB=NIEC(I,2)


!!$       XI=0.5_SP*(XIJE(I,1)+XIJE(I,2))
!!$       YI=0.5_SP*(YIJE(I,1)+YIJE(I,2))
!!$#      if defined (SPHERICAL)
!!$       X1_DP=XIJE(I,1)
!!$       Y1_DP=YIJE(I,1)
!!$       X2_DP=XIJE(I,2)
!!$       Y2_DP=YIJE(I,2)
!!$       CALL ARCC(X2_DP,Y2_DP,X1_DP,Y1_DP,XII,YII)
!!$       XI=XII		
!!$       XTMP  = XI*TPI-VX(IA)*TPI
!!$       XTMP1 = XI-VX(IA)
!!$       IF(XTMP1 >  180.0_SP)THEN
!!$         XTMP = -360.0_SP*TPI+XTMP
!!$       ELSE IF(XTMP1 < -180.0_SP)THEN
!!$         XTMP =  360.0_SP*TPI+XTMP
!!$       END IF	 
!!$
!!$       DXA=XTMP*COS(DEG2RAD*VY(IA))    
!!$       DYA=(YI-VY(IA))*TPI
!!$       XTMP  = XI*TPI-VX(IB)*TPI
!!$       XTMP1 = XI-VX(IB)
!!$       IF(XTMP1 >  180.0_SP)THEN
!!$         XTMP = -360.0_SP*TPI+XTMP
!!$       ELSE IF(XTMP1 < -180.0_SP)THEN
!!$         XTMP =  360.0_SP*TPI+XTMP
!!$       END IF	 
!!$
!!$       DXB=XTMP*COS(DEG2RAD*VY(IB)) 
!!$       DYB=(YI-VY(IB))*TPI
!!$#      else
!!$       DXA=XI-VX(IA)
!!$       DYA=YI-VY(IA)
!!$       DXB=XI-VX(IB)
!!$       DYB=YI-VY(IB)
!!$#      endif

# if !defined (TVD)
       FIJ1=DYE(IA,K)+DLTXNCVE(I,1)*PDYEPX(IA)+DLTYNCVE(I,1)*PDYEPY(IA)
       FIJ2=DYE(IB,K)+DLTXNCVE(I,2)*PDYEPX(IB)+DLTYNCVE(I,2)*PDYEPY(IB)

       S1MIN=MINVAL(DYE(NBSN(IA,1:NTSN(IA)-1),K))
       S1MIN=MIN(S1MIN, DYE(IA,K))
       S1MAX=MAXVAL(DYE(NBSN(IA,1:NTSN(IA)-1),K))
       S1MAX=MAX(S1MAX, DYE(IA,K))
       S2MIN=MINVAL(DYE(NBSN(IB,1:NTSN(IB)-1),K))
       S2MIN=MIN(S2MIN, DYE(IB,K))
       S2MAX=MAXVAL(DYE(NBSN(IB,1:NTSN(IB)-1),K))
       S2MAX=MAX(S2MAX, DYE(IB,K))
       IF(FIJ1 < S1MIN) FIJ1=S1MIN
       IF(FIJ1 > S1MAX) FIJ1=S1MAX
       IF(FIJ2 < S2MIN) FIJ2=S2MIN
       IF(FIJ2 > S2MAX) FIJ2=S2MAX

# else
! ------------------------------------------------------------------------------------------ !
!                                     Drawing the TVD scheme
! ------------------------------------------------------------------------------------------ !
!     If A is upstream of B
       S_upstreamA=DYE(Anear_node(I),K)+PDYEPX(Anear_node(I))*XUAdist(I)+PDYEPY(Anear_node(I))*YUAdist(I)
!     If B is upstream of A
       S_upstreamB=DYE(Bnear_node(I),K)+PDYEPX(Bnear_node(I))*XUBdist(I)+PDYEPY(Bnear_node(I))*YUBdist(I)                

!     Smoothness-parameter
       IF (DYE(IA,K).EQ.DYE(IB,K)) THEN
          smoothrA    = 0.0_SP
          smoothrB    = 0.0_SP
       ELSE IF (DYE(IA,K).LT.1.E-10.AND.DYE(IB,K).LT.1.E-10) THEN
          smoothrA    = 0.0_SP
          smoothrB    = 0.0_SP
       ELSE
          smoothrA    = (DYE(IA,K)-S_upstreamA)/(DYE(IB,K)-DYE(IA,K))
          smoothrB    = (DYE(IB,K)-S_upstreamB)/(DYE(IA,K)-DYE(IB,K))
       END IF

!     Flux-limiter
!     (superbee) - Feel free to use other limiters!
       A_limiter = MAX(0.0_SP, MIN(1.0_SP, 2.0_SP*smoothrA), MIN(2.0_SP, smoothrA))
       B_limiter = MAX(0.0_SP, MIN(1.0_SP, 2.0_SP*smoothrB), MIN(2.0_SP, smoothrB))

!     Estimating the dye at the cellwall in-betweem IA and IB
       FIJ1     = DYE(IA,K)+0.5_SP*A_LIMITER*(DYE(IB,K)-DYE(IA,K))
       FIJ2     = DYE(IB,K)+0.5_SP*B_LIMITER*(DYE(IA,K)-DYE(IB,K))

# endif
       UN=UVN(I,K)
#      if defined (SEMI_IMPLICIT)
       UN1=UVN1(I,K)
#      endif

      ! VISCOF=HORCON*(FACT*(VISCOFF(IA)+VISCOFF(IB))*0.5_SP + FM1)

       ! David moved HPRNU and added VHC
        VISCOF=(FACT*0.5_SP*(VISCOFF(IA)*NN_HVC(IA)+VISCOFF(IB)*NN_HVC(IB)) +  &
	        FM1*0.5_SP*(NN_HVC(IA)+NN_HVC(IB)))

       TXX=0.5_SP*(PDYEPXD(IA)+PDYEPXD(IB))*VISCOF
       TYY=0.5_SP*(PDYEPYD(IA)+PDYEPYD(IB))*VISCOF

       FXX=-DTIJ(I,K)*TXX*DLTYE(I)
       FYY= DTIJ(I,K)*TYY*DLTXE(I)

#      if !defined (SEMI_IMPLICIT)
       EXFLUX=-UN*DTIJ(I,K)* &
          ((1.0_SP+SIGN(1.0_SP,UN))*FIJ2+(1.0_SP-SIGN(1.0_SP,UN))*FIJ1)*0.5_SP+FXX+FYY
#       else
        EXFLUX=-UN*DTIJ(I,K)* &
           ((1.0_SP+SIGN(1.0_SP,UN))*FIJ2+(1.0_SP-SIGN(1.0_SP,UN))*FIJ1)*0.5_SP
        EXFLUX=(1.0_SP-IFCETA)*EXFLUX+IFCETA*(-UN1*DTIJ1(I,K)*  &
	       ((1.0_SP+SIGN(1.0_SP,UN1))*FIJ2+(1.0_SP-SIGN(1.0_SP,UN1))*FIJ1)*0.5_SP)+FXX+FYY
#       endif

       XFLUX(IA,K)=XFLUX(IA,K)+EXFLUX
       XFLUX(IB,K)=XFLUX(IB,K)-EXFLUX

       XFLUX_ADV(IA,K)=XFLUX_ADV(IA,K)+(EXFLUX-FXX-FYY)
       XFLUX_ADV(IB,K)=XFLUX_ADV(IB,K)-(EXFLUX-FXX-FYY)
     END DO

#    if defined (SPHERICAL)
#    if !defined (SEMI_IMPLICIT)
!       CALL ADV_DYE_XY(XFLUX,XFLUX_ADV,PDYEPX,PDYEPY,PDYEPXD,PDYEPYD,VISCOFF,K,0.0_SP)
#    else
!       CALL ADV_DYE_XY(XFLUX,XFLUX_ADV,PDYEPX,PDYEPY,PDYEPXD,PDYEPYD,VISCOFF,K,IFCETA)
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
  ENDDO
 

!   --------------------------------------------------------------------
!   The central difference scheme in vertical advection
!   --------------------------------------------------------------------
!   DO K=1,KBM1
!     DO I=1,M
!#    if defined (WET_DRY)
!       IF(ISWETN(I)*ISWETNT(I) == 1) THEN
!#    endif
!       IF(K == 1) THEN
!         TEMP=-WTS(I,K+1)*(DYE(I,K)*DZ(K+1)+DYE(I,K+1)*DZ(K))/(DZ(K)+DZ(K+1))
!       ELSE IF(K == KBM1) THEN
!         TEMP= WTS(I,K)*(DYE(I,K)*DZ(K-1)+DYE(I,K-1)*DZ(K))/(DZ(K)+DZ(K-1))
!       ELSE
!         TEMP= WTS(I,K)*(DYE(I,K)*DZ(K-1)+DYE(I,K-1)*DZ(K))/(DZ(K)+DZ(K-1))-&
!               WTS(I,K+1)*(DYE(I,K)*DZ(K+1)+DYE(I,K+1)*DZ(K))/(DZ(K)+DZ(K+1))
!       END IF

!       IF(ISONB(I) == 2) THEN
!         XFLUX(I,K)=TEMP*ART1(I)/DZ(K)
!       ELSE
!         XFLUX(I,K)=XFLUX(I,K)+TEMP*ART1(I)/DZ(K)
!       END IF
!#    if defined (WET_DRY)
!       END IF
!#    endif
!     END DO
!   END DO  !! SIGMA LOOP

!   -------------------------------------------------------------------------------
!   -------------------------------------------------------------------------------

#   if defined (MPDATA)

!   --------------------------------------------------------------------------------
!   ---------------------------------------------------------------------------------
!   Songhu {
!   Using smolarkiewicz, P. K; A fully multidimensional positive definite advection
!   transport algorithm with small implicit diffusion, Journal of Computational
!   Physics, 54, 325-362, 1984
        
! The horizontal term of advection is neglected here

   DO K=1,KBM1
     DO I=1,M
       IF(ISONB(I) == 2) THEN
         XFLUX(I,K)=0._SP
       ENDIF
     ENDDO
   ENDDO

! Initialize variables of MPDATA
   DYE_S=0._SP
   DYE_SF=0._SP
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
           TEMP = -(WTS(I,K+1)-ABS(WTS(I,K+1)))*DYE(I,K)   &
                  -(WTS(I,K+1)+ABS(WTS(I,K+1)))*DYE(I,K+1) &
                  +(WTS(I,K)+ABS(WTS(I,K)))*DYE(I,K)    
         ELSE IF(K == KBM1) THEN
           TEMP =  +(WTS(I,K)-ABS(WTS(I,K)))*DYE(I,K-1)     &
                   +(WTS(I,K)+ABS(WTS(I,K)))*DYE(I,K)
         ELSE
           TEMP = -(WTS(I,K+1)-ABS(WTS(I,K+1)))*DYE(I,K)   &
                  -(WTS(I,K+1)+ABS(WTS(I,K+1)))*DYE(I,K+1) &
                  +(WTS(I,K)-ABS(WTS(I,K)))*DYE(I,K-1)     &
                  +(WTS(I,K)+ABS(WTS(I,K)))*DYE(I,K)
         END IF
         TEMP = 0.5_SP*TEMP 

         IF(K == 1)THEN
           SMAX = MAXVAL(DYE(NBSN(I,1:NTSN(I)),K))
           SMIN = MINVAL(DYE(NBSN(I,1:NTSN(I)),K))
           SMAX = MAX(SMAX,DYE(I,K+1),DYE(I,K))
           SMIN = MIN(SMIN,DYE(I,K+1),DYE(I,K))
         ELSE IF(K == KBM1)THEN
           SMAX = MAXVAL(DYE(NBSN(I,1:NTSN(I)),K))
           SMIN = MINVAL(DYE(NBSN(I,1:NTSN(I)),K))
           SMAX = MAX(SMAX,DYE(I,K-1),DYE(I,K))
           SMIN = MIN(SMIN,DYE(I,K-1),DYE(I,K))
         ELSE
           SMAX = MAXVAL(DYE(NBSN(I,1:NTSN(I)),K))
           SMIN = MINVAL(DYE(NBSN(I,1:NTSN(I)),K))
           SMAX = MAX(SMAX,DYE(I,K+1),DYE(I,K-1),DYE(I,K))
           SMIN = MIN(SMIN,DYE(I,K+1),DYE(I,K-1),DYE(I,K))
         END IF

        ZZZFLUX(I,K)= TEMP*(DTI/DT(I))/DZ(I,K) + XFLUX(I,K)/ART1(I)*(DTI/DT(I))/DZ(I,K) 
        XXXX= ZZZFLUX(I,K)*DT(I)/DTFA(I)+DYE(I,K)-DYE(I,K)*DT(I)/DTFA(I) 

        ! Added by Akvaplan 2018 to avoid single precision-crash
        IF((ABS(XXXX).LT.ABS(SMAX-DYE(I,K))).AND.(XXXX.LT.0.0_SP)) THEN
          DYE_SF(I,K) = DYE(I,K)-XXXX
        ELSE IF((ABS(XXXX).LT.ABS(SMIN-DYE(I,K))).AND.(XXXX.GT.0.0_SP)) THEN 
          DYE_SF(I,K) = DYE(I,K)-XXXX
        ELSE IF(XXXX.EQ.0) THEN
          DYE_SF(I,K) = DYE(I,K)
        ELSE
          BETA(I,K)=0.5*(1.-SIGN(1.,XXXX)) * (SMAX-DYE(I,K))/(ABS(XXXX)+1.E-10) &
                   +0.5*(1.-SIGN(1.,-XXXX)) * (DYE(I,K)-SMIN)/(ABS(XXXX)+1.E-10)

          DYE_SF(I,K)=DYE(I,K)-MIN(1.,BETA(I,K))*XXXX
        END IF

#    if defined (WET_DRY)
       END IF
#    endif
     END DO
   END DO  !! SIGMA LOOP
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   NTERA = 4
   DO ITERA=1,NTERA   !! Smolaricizw Loop 
     IF(ITERA == 1)THEN
       WWWSF = WTS
       DYE_S =DYE_SF
       DTWWWS = DT
     ELSE
       WWWSF = WWWS
       DYE_S =DYE_SF
       DTWWWS = DTFA
     END IF
     DO K=2,KBM1
       DO I=1,M
         TEMP=ABS(WWWSF(I,K))-DTI*(WWWSF(I,K))*(WWWSF(I,K))/DZ(I,K)/DTWWWS(I)
         WWWS(I,K)=TEMP*(DYE_S(I,K-1)-DYE_S(I,K))/(ABS(DYE_S(I,K-1))+     &
	           ABS(DYE_S(I,K))+1.E-14)

         IF(TEMP < 0.0_SP .OR. DYE_S(I,K) == 0.0_SP)THEN 
           WWWS(I,K)=0.0_SP 
         END IF
       END DO
     END DO
     DO I=1,M
       WWWS(I,1)=0.0_SP
     END DO


     DO I=1,M
       SMAX = MAXVAL(DYE(NBSN(I,1:NTSN(I)),1))
       SMIN = MINVAL(DYE(NBSN(I,1:NTSN(I)),1))
       SMAX = MAX(SMAX,DYE(I,2),DYE(I,1))
       SMIN = MIN(SMIN,DYE(I,2),DYE(I,1))
 
       TEMP=0.5*((WWWS(I,2)+ABS(WWWS(I,2)))*DYE_S(I,2))*(DTI/DTFA(I))/DZ(I,1)
       BETAIN(I,1)=(SMAX-DYE_S(I,1))/(TEMP+1.E-10)

       TEMP=0.5*((WWWS(I,1)+ABS(WWWS(I,1)))*DYE_S(I,1)-    &
            (WWWS(I,2)-ABS(WWWS(I,2)))*DYE_S(I,1))*(DTI/DTFA(I))/DZ(I,1)
       BETAOUT(I,1)=(DYE_S(I,1)-SMIN)/(TEMP+1.E-10)

       WWWSF(I,1)=0.5*MIN(1.,BETAOUT(I,1))*(WWWS(I,1)+ABS(WWWS(I,1))) + &
                  0.5*MIN(1.,BETAIN(I,1))*(WWWS(I,1)-ABS(WWWS(I,1)))
     END DO

     DO K=2,KBM1-1
       DO I=1,M
         SMAX = MAXVAL(DYE(NBSN(I,1:NTSN(I)),K))
         SMIN = MINVAL(DYE(NBSN(I,1:NTSN(I)),K))
         SMAX = MAX(SMAX,DYE(I,K+1),DYE(I,K-1),DYE(I,K))
         SMIN = MIN(SMIN,DYE(I,K+1),DYE(I,K-1),DYE(I,K))
 
         TEMP=0.5*((WWWS(I,K+1)+ABS(WWWS(I,K+1)))*DYE_S(I,K+1)-   &
	      (WWWS(I,K)-ABS(WWWS(I,K)))*DYE_S(I,K-1))*(DTI/DTFA(I))/DZ(I,K)
         BETAIN(I,K)=(SMAX-DYE_S(I,K))/(TEMP+1.E-10)

         TEMP=0.5*((WWWS(I,K)+ABS(WWWS(I,K)))*DYE_S(I,K)-         &
	      (WWWS(I,K+1)-ABS(WWWS(I,K+1)))*DYE_S(I,K))*(DTI/DTFA(I))/DZ(I,K)
         BETAOUT(I,K)=(DYE_S(I,K)-SMIN)/(TEMP+1.E-10)

         WWWSF(I,K)=0.5*MIN(1.,BETAIN(I,K-1),BETAOUT(I,K))*(WWWS(I,K)+ABS(WWWS(I,K))) + &
                    0.5*MIN(1.,BETAIN(I,K),BETAOUT(I,K-1))*(WWWS(I,K)-ABS(WWWS(I,K)))
       END DO
     END DO
     
     K=KBM1
     DO I=1,M
       SMAX = MAXVAL(DYE(NBSN(I,1:NTSN(I)),K))
       SMIN = MINVAL(DYE(NBSN(I,1:NTSN(I)),K))
       SMAX = MAX(SMAX,DYE(I,K-1),DYE(I,K))
       SMIN = MIN(SMIN,DYE(I,K-1),DYE(I,K))
 
       TEMP=0.5*((WWWS(I,K+1)+ABS(WWWS(I,K+1)))*DYE_S(I,K+1)-   &
            (WWWS(I,K)-ABS(WWWS(I,K)))*DYE_S(I,K-1))*(DTI/DTFA(I))/DZ(I,K)
       BETAIN(I,K)=(SMAX-DYE_S(I,K))/(TEMP+1.E-10)

       TEMP=0.5*((WWWS(I,K)+ABS(WWWS(I,K)))*DYE_S(I,K)-         &
	    (WWWS(I,K+1)-ABS(WWWS(I,K+1)))*DYE_S(I,K))*(DTI/DTFA(I))/DZ(I,K)
       BETAOUT(I,K)=(DYE_S(I,K)-SMIN)/(TEMP+1.E-10)

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
             TEMP = -(WWWS(I,K+1)-ABS(WWWS(I,K+1)))*DYE_S(I,K)   &
                    -(WWWS(I,K+1)+ABS(WWWS(I,K+1)))*DYE_S(I,K+1) &
                    +(WWWS(I,K)+ABS(WWWS(I,K)))*DYE_S(I,K)
           ELSE IF(K == KBM1) THEN
             TEMP =  +(WWWS(I,K)-ABS(WWWS(I,K)))*DYE_S(I,K-1)     &
                     +(WWWS(I,K)+ABS(WWWS(I,K)))*DYE_S(I,K)
           ELSE
             TEMP = -(WWWS(I,K+1)-ABS(WWWS(I,K+1)))*DYE_S(I,K)   &
                    -(WWWS(I,K+1)+ABS(WWWS(I,K+1)))*DYE_S(I,K+1) &
                    +(WWWS(I,K)-ABS(WWWS(I,K)))*DYE_S(I,K-1)     &
                    +(WWWS(I,K)+ABS(WWWS(I,K)))*DYE_S(I,K)
           END IF
           TEMP = 0.5_SP*TEMP
           DYE_SF(I,K)=(   DYE_S(I,K)-TEMP*(DTI/DTFA(I))/DZ(I,K)   ) 
#      if defined (WET_DRY)
         END IF
#      endif
       END DO
     END DO  !! SIGMA LOOP
   END DO  !! Smolarvizw Loop


!--------------------------------------------------------------------------
! End of smolarkiewicz upwind loop
!--------------------------------------------------------------------------

      DYE = DYE - OFFS
      DYE_SF = DYE_SF - OFFS
!--------------------------------------------------------------------------
#  endif

# if ! defined(MPDATA)
!--------------------------------------------------------------------
!   The central difference scheme in vertical advection
!--------------------------------------------------------------------

!!$  ORIGINAL - before switch to match adv_t
!!$   DO K=1,KBM1
!!$     DO I=1,M
!!$#    if defined (WET_DRY)
!!$       IF(ISWETN(I)*ISWETNT(I) == 1) THEN
!!$#    endif
!!$       IF(K == 1) THEN
!!$         TEMP=-WTS(I,K+1)*(DYE(I,K)*DZ(I,K+1)+DYE(I,K+1)*DZ(I,K))/  &
!!$	      (DZ(I,K)+DZ(I,K+1))
!!$       ELSE IF(K == KBM1) THEN
!!$         TEMP= WTS(I,K)*(DYE(I,K)*DZ(I,K-1)+DYE(I,K-1)*DZ(I,K))/    &
!!$	      (DZ(I,K)+DZ(I,K-1))
!!$       ELSE
!!$         TEMP= WTS(I,K)*(DYE(I,K)*DZ(I,K-1)+DYE(I,K-1)*DZ(I,K))/    &
!!$	      (DZ(I,K)+DZ(I,K-1))-&
!!$               WTS(I,K+1)*(DYE(I,K)*DZ(I,K+1)+DYE(I,K+1)*DZ(I,K))/  &
!!$	      (DZ(I,K)+DZ(I,K+1))
!!$       END IF
!!$
!!$       IF(ISONB(I) == 2) THEN
!!$         XFLUX(I,K)=TEMP*ART1(I)    !/DZ(K)
!!$       ELSE
!!$         XFLUX(I,K)=XFLUX(I,K)+TEMP*ART1(I)    !/DZ(K)
!!$       END IF
!!$#    if defined (WET_DRY)
!!$       END IF
!!$#    endif
!!$     END DO
!!$   END DO  !! SIGMA LOOP


   
   DO I=1,M
#    if defined (WET_DRY)
      IF(ISWETN(I)*ISWETNT(I) == 1) THEN
#    endif
         
      DO K=1, KBM1
         
         IF(K == 1) THEN
            TEMP=-WTS(I,K+1)*(DYE(I,K)*DZ(I,K+1)+DYE(I,K+1)*DZ(I,K))/  &
                 (DZ(I,K)+DZ(I,K+1))
         ELSE IF(K == KBM1) THEN
            TEMP= WTS(I,K)*(DYE(I,K)*DZ(I,K-1)+DYE(I,K-1)*DZ(I,K))/    &
                 (DZ(I,K)+DZ(I,K-1))
         ELSE
            TEMP= WTS(I,K)*(DYE(I,K)*DZ(I,K-1)+DYE(I,K-1)*DZ(I,K))/    &
                 (DZ(I,K)+DZ(I,K-1))-&
                 WTS(I,K+1)*(DYE(I,K)*DZ(I,K+1)+DYE(I,K+1)*DZ(I,K))/  &
                 (DZ(I,K)+DZ(I,K+1))
         END IF
         
         IF(ISONB(I) == 2) THEN
            XFLUX(I,K)=TEMP*ART1(I)    !/DZ(K)
         ELSE
            XFLUX(I,K)=XFLUX(I,K)+TEMP*ART1(I)    !/DZ(K)
         END IF
      END DO

#    if defined (WET_DRY)
      END IF
#    endif
   END DO  !! SIGMA LOOP



!
!--Set Boundary Conditions-For Fresh Water Flux--------------------------------!
!
!   IF(POINT_ST_TYPE == 'calculated') THEN
!     IF(INFLOW_TYPE == 'node') THEN
!       IF(NUMQBC > 0) THEN
!         DO J=1,NUMQBC
!           JJ=INODEQ(J)
!           STPOINT=SDIS(J)
!           DO K=1,KBM1
!!             XFLUX(JJ,K)=XFLUX(JJ,K) - QDIS(J)*VQDIST(J,K)*STPOINT/DZ(K)
!             XFLUX(JJ,K)=XFLUX(JJ,K) - QDIS(J)*VQDIST(J,K)*STPOINT
!           END DO
!         END DO
!       END IF
!     ELSE IF(INFLOW_TYPE == 'edge') THEN
!       IF(NUMQBC > 1) THEN
!         DO J=1,NUMQBC
!           J1=N_ICELLQ(J,1)
!           J2=N_ICELLQ(J,2)
!           STPOINT=SDIS(J) !!ASK LIU SHOULD THIS BE STPOINT1(J1)/STPOINT2(J2)
!           DO K=1,KBM1
!!             XFLUX(J1,K)=XFLUX(J1,K)-QDIS(J)*RDISQ(J,1)*VQDIST(J,K)*STPOINT/DZ(K)
!!             XFLUX(J2,K)=XFLUX(J2,K)-QDIS(J)*RDISQ(J,2)*VQDIST(J,K)*STPOINT/DZ(K)
!             XFLUX(J1,K)=XFLUX(J1,K)-QDIS(J)*RDISQ(J,1)*VQDIST(J,K)*STPOINT
!             XFLUX(J2,K)=XFLUX(J2,K)-QDIS(J)*RDISQ(J,2)*VQDIST(J,K)*STPOINT
!           END DO
!         END DO
!       END IF
!     END IF
!   END IF

#  endif
! NOT MPDATA


!
!--Update Dye-------------------------------------------------------------!
!

#  if defined (WET_DRY)
   DO I=1,M
     IF(ISWETN(I)*ISWETNT(I) == 1 )THEN
       DO K=1,KBM1

#  if defined (MPDATA)
          DYEF(I,K)=DYE_SF(I,K)
#  else     
          

#      if !defined (SEMI_IMPLICIT)
            DYEF(I,K)=(DYE(I,K)-XFLUX(I,K)/ART1(I)*(DTI/(DT(I)*DZ(I,K))))*(DT(I)/DTFA(I))
#      else
!            DYEF(I,K)=(DYE(I,K)-XFLUX(I,K)/ART1(I)*(DTI/(DT(I)*DZ(I,K))))*(DT(I)/DTFA(I))
            DYEF(I,K)=(TSC(I,K)-RK_TS(STAGE)*XFLUX(I,K)/ART1(I)*(DTI/(DT(I)*DZ(I,K))))*(DT(I)/DTFA(I))
#      endif
        
          
#  endif 
!(MPDATA)
         END DO
      ELSE  

         DO K=1,KBM1
#      if !defined (SEMI_IMPLICIT)
            DYEF(I,K)=DYE(I,K)
#      else
!            DYEF(I,K)=DYE(I,K)
            DYEF(I,K)=TSC(I,K)
#      endif
         END DO
      END IF
 END DO
#  else  
!NOT WET DRY
 DO I=1,M
    DO K=1,KBM1

#    if !defined (MPDATA)     
#      if !defined (SEMI_IMPLICIT)
       DYEF(I,K)=(DYE(I,K)-XFLUX(I,K)/ART1(I)*(DTI/(DT(I)*DZ(I,K))))*(DT(I)/DTFA(I))
#      else
!       DYEF(I,K)=(DYE(I,K)-XFLUX(I,K)/ART1(I)*(DTI/(DT(I)*DZ(I,K))))*(DT(I)/DTFA(I))
       DYEF(I,K)=(TSC(I,K)-RK_TS(STAGE)*XFLUX(I,K)/ART1(I)*(DTI/(DT(I)*DZ(I,K))))*(DT(I)/DTFA(I))
#      endif
#    else
       DYEF(I,K)=DYE_SF(I,K)
#    endif              

    END DO
 END DO
#  endif


!   
!---- specify the source term------------------------------------------|
!
!   IF(IINT.GE.IINT_SPE_DYE_B.AND.IINT.LE.IINT_SPE_DYE_E) THEN
!     DO KK=1, KSPE_DYE
!        K=K_SPECIFY(KK)
!        DO J=1,MSPE_DYE
!           I = M_SPECIFY(J)
!           DYEF(I,K)= DYE_SOURCE_TERM
!        ENDDO
!     ENDDO
!   ENDIF

  IF(DBG_SET(DBG_SBR)) WRITE(IPT,*)"End: adv_dye"

   RETURN
   END SUBROUTINE ADV_DYE
!==============================================================================|

!==============================================================================|
!   Calculate Advection and Horizontal Diffusion Terms for DYE                 |
!==============================================================================|

   SUBROUTINE ADV_DYE_RK(DYEB)               

!------------------------------------------------------------------------------|

   USE ALL_VARS
   USE MOD_UTILS
   USE BCS
   USE MOD_OBCS
   USE MOD_PAR
   USE MOD_WD
   USE MOD_SPHERICAL
   USE MOD_NORTHPOLE
#  if defined (TVD)
   USE MOD_TVD
#  endif

   IMPLICIT NONE
   REAL(SP), DIMENSION(0:MT,KB)     :: XFLUX,XFLUX_ADV,DYEB
   REAL(SP), DIMENSION(M)           :: PUPX,PUPY,PVPX,PVPY  
   REAL(SP), DIMENSION(M)           :: PDYEPX,PDYEPY,PDYEPXD,PDYEPYD,VISCOFF
   REAL(SP), DIMENSION(3*(NT),KBM1) :: DTIJ 
   REAL(SP), DIMENSION(3*(NT),KBM1) :: UVN
   REAL(SP) :: UTMP,VTMP,SITAI,FFD,FF1 
   REAL(SP) :: DXA,DYA,DXB,DYB,FIJ1,FIJ2,UN
   REAL(SP) :: TXX,TYY,FXX,FYY,VISCOF,EXFLUX,TEMP,STPOINT
   REAL(SP) :: FACT,FM1
   REAL(SP) :: s1min, s1max, s2min, s2max,SMIN,SMAX,xxxx
   INTEGER  :: I,I1,I2,IA,IB,J,J1,J2,K,JTMP,JJ,KK,IP,II
#  if defined (WET_DRY)
   INTEGER :: N_NTVE
#  endif
# if defined (TVD)
   REAL(SP) :: S_upstreamA,S_upstreamB
   REAL(SP) :: A_limiter,B_limiter,smoothrA,smoothrB
#  endif

  IF(DBG_SET(DBG_SBR)) WRITE(IPT,*)"Start: adv_dye_rk"
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
!--Initialize Fluxes and dyef----------------------------------------------------------!
!
   XFLUX     = 0.0_SP
   XFLUX_ADV = 0.0_SP
   DYEF      = 0.0_SP
!
!--Loop Over Control Volume Sub-Edges And Calculate Normal Velocity------------!
!
   DO I=1,NCV
     I1=NTRG(I)
     DO K=1,KBM1
       DTIJ(I,K)=DT1(I1)*DZ1(I1,K)
       UVN(I,K) = V(I1,K)*DLTXE(I) - U(I1,K)*DLTYE(I)
     END DO
   END DO
!
!--Calculate the Advection and Horizontal Diffusion Terms----------------------!
!

   DO K=1,KBM1
      PDYEPX  = 0.0_SP 
      PDYEPY  = 0.0_SP 
      PDYEPXD = 0.0_SP 
      PDYEPYD = 0.0_SP
     DO I=1,M
       DO J=1,NTSN(I)-1
         I1=NBSN(I,J)
         I2=NBSN(I,J+1)

#    if defined (WET_DRY)
         IF(ISWETN(I1) == 0 .AND. ISWETN(I2) == 1)THEN
          FFD=0.5_SP*(DYE(I,K)+DYE(I2,K))!-DYEMEAN(I,K)-DYEMEAN(I2,K))
          FF1=0.5_SP*(DYE(I,K)+DYE(I2,K))
	 ELSE IF(ISWETN(I1) == 1 .AND. ISWETN(I2) == 0)THEN
          FFD=0.5_SP*(DYE(I1,K)+DYE(I,K))!-DYEMEAN(I1,K)-DYEMEAN(I,K))
          FF1=0.5_SP*(DYE(I1,K)+DYE(I,K))
	 ELSE IF(ISWETN(I1) == 0 .AND. ISWETN(I2) == 0)THEN
          FFD=0.5_SP*(DYE(I,K)+DYE(I,K))!-DYEMEAN(I,K)-DYEMEAN(I,K))
          FF1=0.5_SP*(DYE(I,K)+DYE(I,K))
	 ELSE
          FFD=0.5_SP*(DYE(I1,K)+DYE(I2,K))!-DYEMEAN(I1,K)-DYEMEAN(I2,K))
          FF1=0.5_SP*(DYE(I1,K)+DYE(I2,K))
	 END IF 
#    else	 
         FFD=0.5_SP*(DYE(I1,K)+DYE(I2,K))!-DYEMEAN(I1,K)-DYEMEAN(I2,K))
         FF1=0.5_SP*(DYE(I1,K)+DYE(I2,K))
#    endif

         PDYEPX(I)=PDYEPX(I)+FF1*DLTYTRIE(i,j)
         PDYEPY(I)=PDYEPY(I)+FF1*DLTXTRIE(i,j)
         PDYEPXD(I)=PDYEPXD(I)+FFD*DLTYTRIE(i,j)
         PDYEPYD(I)=PDYEPYD(I)+FFD*DLTXTRIE(i,j)

       END DO
       PDYEPX(I)=PDYEPX(I)/ART2(I)
       PDYEPY(I)=PDYEPY(I)/ART2(I)
       PDYEPXD(I)=PDYEPXD(I)/ART2(I)
       PDYEPYD(I)=PDYEPYD(I)/ART2(I)
     END DO
          
     IF(K == KBM1)THEN
       DO I=1,M
         PFPXB(I) = PDYEPX(I)
         PFPYB(I) = PDYEPY(I)
       END DO
     END IF

     DO I=1,M
       VISCOFF(I)=VISCOFH(I,K)
     END DO

     IF(K == KBM1) THEN
       AH_BOTTOM(1:M) = (FACT*VISCOFF(1:M) + FM1)*NN_HVC(1:M)
    END IF


     DO I=1,NCV_I
       IA=NIEC(I,1)
       IB=NIEC(I,2)

# if !defined (TVD)
       FIJ1=DYE(IA,K)+DLTXNCVE(I,1)*PDYEPX(IA)+DLTYNCVE(I,1)*PDYEPY(IA)
       FIJ2=DYE(IB,K)+DLTXNCVE(I,2)*PDYEPX(IB)+DLTYNCVE(I,2)*PDYEPY(IB)

       S1MIN=MINVAL(DYE(NBSN(IA,1:NTSN(IA)-1),K))
       S1MIN=MIN(S1MIN, DYE(IA,K))
       S1MAX=MAXVAL(DYE(NBSN(IA,1:NTSN(IA)-1),K))
       S1MAX=MAX(S1MAX, DYE(IA,K))
       S2MIN=MINVAL(DYE(NBSN(IB,1:NTSN(IB)-1),K))
       S2MIN=MIN(S2MIN, DYE(IB,K))
       S2MAX=MAXVAL(DYE(NBSN(IB,1:NTSN(IB)-1),K))
       S2MAX=MAX(S2MAX, DYE(IB,K))
       IF(FIJ1 < S1MIN) FIJ1=S1MIN
       IF(FIJ1 > S1MAX) FIJ1=S1MAX
       IF(FIJ2 < S2MIN) FIJ2=S2MIN
       IF(FIJ2 > S2MAX) FIJ2=S2MAX

# else
! ------------------------------------------------------------------------------------------ !
!                                     Drawing the TVD scheme
! ------------------------------------------------------------------------------------------ !
!     If A is upstream of B
       S_upstreamA=DYE(Anear_node(I),K)+PDYEPX(Anear_node(I))*XUAdist(I)+PDYEPY(Anear_node(I))*YUAdist(I)
!     If B is upstream of A
       S_upstreamB=DYE(Bnear_node(I),K)+PDYEPX(Bnear_node(I))*XUBdist(I)+PDYEPY(Bnear_node(I))*YUBdist(I)                

!     Smoothness-parameter
       IF (DYE(IA,K).EQ.DYE(IB,K)) THEN
          smoothrA    = 0.0_SP
          smoothrB    = 0.0_SP
       ELSE IF (DYE(IA,K).LT.1.E-10.AND.DYE(IB,K).LT.1.E-10) THEN
          smoothrA    = 0.0_SP
          smoothrB    = 0.0_SP
       ELSE
          smoothrA    = (DYE(IA,K)-S_upstreamA)/(DYE(IB,K)-DYE(IA,K))
          smoothrB    = (DYE(IB,K)-S_upstreamB)/(DYE(IA,K)-DYE(IB,K))
       END IF

!     Flux-limiter
!     (superbee) - Feel free to use other limiters!
       A_limiter = MAX(0.0_SP, MIN(1.0_SP, 2.0_SP*smoothrA), MIN(2.0_SP, smoothrA))
       B_limiter = MAX(0.0_SP, MIN(1.0_SP, 2.0_SP*smoothrB), MIN(2.0_SP, smoothrB))

!     Estimating the dye at the cellwall in-betweem IA and IB
       FIJ1     = DYE(IA,K)+0.5_SP*A_LIMITER*(DYE(IB,K)-DYE(IA,K))
       FIJ2     = DYE(IB,K)+0.5_SP*B_LIMITER*(DYE(IA,K)-DYE(IB,K))

# endif
       UN=UVN(I,K)

       ! David moved HPRNU and added VHC
        VISCOF=(FACT*0.5_SP*(VISCOFF(IA)*NN_HVC(IA)+VISCOFF(IB)*NN_HVC(IB)) +  &
	        FM1*0.5_SP*(NN_HVC(IA)+NN_HVC(IB)))

       TXX=0.5_SP*(PDYEPXD(IA)+PDYEPXD(IB))*VISCOF
       TYY=0.5_SP*(PDYEPYD(IA)+PDYEPYD(IB))*VISCOF

       FXX=-DTIJ(I,K)*TXX*DLTYE(I)
       FYY= DTIJ(I,K)*TYY*DLTXE(I)

       EXFLUX=-UN*DTIJ(I,K)* &
          ((1.0_SP+SIGN(1.0_SP,UN))*FIJ2+(1.0_SP-SIGN(1.0_SP,UN))*FIJ1)*0.5_SP+FXX+FYY

       XFLUX(IA,K)=XFLUX(IA,K)+EXFLUX
       XFLUX(IB,K)=XFLUX(IB,K)-EXFLUX

       XFLUX_ADV(IA,K)=XFLUX_ADV(IA,K)+(EXFLUX-FXX-FYY)
       XFLUX_ADV(IB,K)=XFLUX_ADV(IB,K)-(EXFLUX-FXX-FYY)
     END DO

#    if defined (SPHERICAL)
!       CALL ADV_DYE_XY(XFLUX,XFLUX_ADV,PDYEPX,PDYEPY,PDYEPXD,PDYEPYD,VISCOFF,K,0.0_SP)
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
  ENDDO

!--------------------------------------------------------------------
!   The central difference scheme in vertical advection
!--------------------------------------------------------------------

   DO I=1,M
#    if defined (WET_DRY)
      IF(ISWETN(I)*ISWETNT(I) == 1) THEN
#    endif
         
      DO K=1, KBM1
         
         IF(K == 1) THEN
            TEMP=-WTS(I,K+1)*(DYE(I,K)*DZ(I,K+1)+DYE(I,K+1)*DZ(I,K))/  &
                 (DZ(I,K)+DZ(I,K+1))
         ELSE IF(K == KBM1) THEN
            TEMP= WTS(I,K)*(DYE(I,K)*DZ(I,K-1)+DYE(I,K-1)*DZ(I,K))/    &
                 (DZ(I,K)+DZ(I,K-1))
         ELSE
            TEMP= WTS(I,K)*(DYE(I,K)*DZ(I,K-1)+DYE(I,K-1)*DZ(I,K))/    &
                 (DZ(I,K)+DZ(I,K-1))-&
                 WTS(I,K+1)*(DYE(I,K)*DZ(I,K+1)+DYE(I,K+1)*DZ(I,K))/  &
                 (DZ(I,K)+DZ(I,K+1))
         END IF
         
         IF(ISONB(I) == 2) THEN
            XFLUX(I,K)=TEMP*ART1(I)    !/DZ(K)
         ELSE
            XFLUX(I,K)=XFLUX(I,K)+TEMP*ART1(I)    !/DZ(K)
         END IF
      END DO

#    if defined (WET_DRY)
      END IF
#    endif
   END DO  !! SIGMA LOOP

!
!--Update Dye-------------------------------------------------------------!
!

#  if defined (WET_DRY)
   DO I=1,M
     IF(ISWETN(I)*ISWETNT(I) == 1 )THEN
       DO K=1,KBM1
         DYEF(I,K)=(DYEB(I,K)-XFLUX(I,K)/ART1(I)*(DTI/(DT(I)*DZ(I,K))))*(DT(I)/DTFA(I))
       END DO
     ELSE  
       DO K=1,KBM1
         DYEF(I,K)=DYEB(I,K)
       END DO
     END IF
   END DO
#  else  
!NOT WET DRY
   DO I=1,M
     DO K=1,KBM1
       DYEF(I,K)=(DYEB(I,K)-XFLUX(I,K)/ART1(I)*(DTI/(DT(I)*DZ(I,K))))*(DT(I)/DTFA(I))
     END DO
   END DO
#  endif

  IF(DBG_SET(DBG_SBR)) WRITE(IPT,*)"End: adv_dye_rk"

   RETURN
   END SUBROUTINE ADV_DYE_RK
!==============================================================================|



!==============================================================================|
!     this subroutine is used to calculate the dye                             !
!     by including vertical diffusion implicitly.                              !
!==============================================================================|

   SUBROUTINE VDIF_DYE(F)                

!------------------------------------------------------------------------------|

   USE ALL_VARS
   USE MOD_UTILS
   USE BCS
#  if defined (WET_DRY)
   USE MOD_WD
#  endif
   IMPLICIT NONE
   INTEGER :: I,K,J,KI,KK
   REAL(DP) :: TMP,TMP1,TMP2,TMP3,QTMP,GW,ZDEP,FKH,UMOLPR
   REAL(SP), DIMENSION(0:MT,KB)  :: F
   REAL(DP), DIMENSION(M,KB)     :: FF,AF,CF,VHF,VHPF,RAD
   REAL(DP), DIMENSION(M)        :: KHBOTTOM,WFSURF,SWRADF

   IF(DBG_SET(DBG_SBR)) WRITE(IPT,*)"Start: vdif_dye :"

   UMOLPR = UMOL*1.E0_SP

!------------------------------------------------------------------------------!
!                                                                              !
!        the following section solves the equation                             !
!         dti*(kh*f')'-f=-fb                                                   !
!                                                                              !
!------------------------------------------------------------------------------!


   DO K = 2, KBM1
     DO I = 1, M
#  if !defined (WET_DRY)
!       IF (D(I) > 0.0_SP) THEN
#  else
       IF(ISWETN(I) == 1)THEN
#  endif
!         FKH=0.0_SP
!         DO J=1,NTVE(I)
!           FKH=FKH+KH(NBVE(I,J),K)
!         END DO
!         FKH=FKH/FLOAT(NTVE(I))
         FKH = KH(I,K)

         IF(K == KBM1) THEN
           KHBOTTOM(I)=FKH
         END IF

!         AF(I,K-1)=-DTI*(FKH+UMOLPR)/(DZ(K-1)*DZZ(K-1)*D(I)*D(I))
!         CF(I,K)=-DTI*(FKH+UMOLPR)/(DZ(K)*DZZ(K-1)*D(I)*D(I))
         AF(I,K-1)=-DTI*(FKH+UMOLPR)/(DZ(I,K-1)*DZZ(I,K-1)*D(I)*D(I))
         CF(I,K)=-DTI*(FKH+UMOLPR)/(DZ(I,K)*DZZ(I,K-1)*D(I)*D(I))
# if defined(WET_DRY)
       END IF
# endif
     END DO
   END DO


!------------------------------------------------------------------------------!
!     the net heat flux input.                                                 !
!     the method shown below can be used when we turn off the                  !
!     body force in subroutine advt. be sure this method could                 !
!     cause the surface overheated if the vertical resolution                  !
!     is not high enough.                                                      !
!------------------------------------------------------------------------------!

     SWRADF = 0.0_SP
     WFSURF = 0.0_SP
     RAD    = 0.0_SP


!------------------------------------------------------------------------------!
!   surface bcs; wfsurf                                                        !
!------------------------------------------------------------------------------!

   DO I = 1, M
#  if !defined (WET_DRY)
!     IF (D(I) > 0.0_SP) THEN
#  else
     IF(ISWETN(I) == 1)THEN
#  endif
       VHF(I,1) = AF(I,1) / (AF(I,1)-1.)
       VHPF(I,1) = -DTI *(WFSURF(I)-SWRADF(I) &
                   +RAD(I,1)-RAD(I,2)) / (-DZ(I,1)*D(I)) - F(I,1)
       VHPF(I,1) = VHPF(I,1) / (AF(I,1)-1.)
#  if defined (WET_DRY)
     END IF
# endif
   END DO

   DO K = 2, KBM2
     DO I = 1, M
#  if !defined (WET_DRY)
!       IF (D(I) > 0.0_SP) THEN
#  else
       IF(ISWETN(I) == 1)THEN
#  endif
         VHPF(I,K)=1./ (AF(I,K)+CF(I,K)*(1.-VHF(I,K-1))-1.)
         VHF(I,K) = AF(I,K) * VHPF(I,K)
         VHPF(I,K) = (CF(I,K)*VHPF(I,K-1)-DBLE(F(I,K)) &
                     +DTI*(RAD(I,K)-RAD(I,K+1))/(D(I)*DZ(I,K)))*VHPF(I,K)
#  if defined (WET_DRY)
       END IF
# endif
     END DO
   END DO


!!$   DO  K = 1, KBM1
!!$     DO  I = 1, M
!!$#  if !defined (WET_DRY)
!!$!       IF (D(I) > 0.0_SP) THEN
!!$#  else
!!$       IF(ISWETN(I) == 1)THEN
!!$#  endif
!!$         FF(I,K) = F(I,K)
!!$#  if defined (WET_DRY)
!!$       END IF
!!$# endif
!!$     END DO
!!$   END DO


#  if !defined (WET_DRY)
  FF(1:M,1:KBM1) = F(1:M,1:KBM1)
#  else
  DO  K = 1, KBM1
     DO  I = 1, M
        IF(ISWETN(I) == 1)THEN
           FF(I,K) = F(I,K)
        END IF
     END DO
  END DO
#  endif



   DO I = 1, M
#  if !defined (WET_DRY)
     IF (ISONB(I) /= 2) THEN
#  else
     IF(ISWETN(I) == 1 .AND.ISONB(I) /= 2)THEN
#  endif
       TMP1=PFPXB(I)*COS(SITA_GD(I))+PFPYB(I)*SIN(SITA_GD(I))
       TMP2=AH_BOTTOM(I)*PHPN(I)
       TMP3=KHBOTTOM(I)+UMOLPR+AH_BOTTOM(I)*PHPN(I)*PHPN(I)
       TMP=TMP1*TMP2/TMP3*(KHBOTTOM(I)+UMOLPR)
! 
!       IF (TMP1 > 0.0_SP) TMP=0.0_SP
       TMP=0.0_SP
! 
       GW=0.0_SP
!!$       IF(IBFW > 0) THEN
!!$!         DO J=1,IBFW
!!$!           IF(I == NODE_BFW(J)) THEN
!!$!!             QTMP=-(F(I,KBM1)*D1(I)*DZ(KBM1)*BFWDIS(J))/ &
!!$!!                   (D1(I)*DZ(KBM1)*ART1(I)+BFWDIS(J))
!!$!!             GW=DTI/D1(I)/DZ(KBM1)*QTMP
!             QTMP=-(F(I,KBM1)*D(I)*DZ(I,KBM1)*BFWDIS(J))/ &
!!$!                   (D(I)*DZ(I,KBM1)*ART1(I)+BFWDIS(J))
!!$!             GW=DTI/D(I)/DZ(I,KBM1)*QTMP
!!$!             TMP=0.0_SP
!!$!           END IF
!!$!         END DO
!!$!       END IF

       FF(I,KBM1) = ((CF(I,KBM1)*VHPF(I,KBM2)-FF(I,KBM1)-GW &
               +DTI*(RAD(I,KBM1)-RAD(I,KB)-TMP)/(D(I)*DZ(I,KBM1))) &
                /(CF(I,KBM1)*(1.-VHF(I,KBM2))-1.))
     END IF
   END DO

   DO  K = 2, KBM1
     KI = KB - K
     DO  I = 1, M
#  if !defined (WET_DRY)
       IF (D(I) > 0.0_SP .AND.ISONB(I) /= 2) THEN
#  else
       IF(ISWETN(I) == 1 .AND.ISONB(I) /= 2)THEN
#  endif
         FF(I,KI) = (VHF(I,KI)*FF(I,KI+1)+VHPF(I,KI))
       END IF
     END DO
   END DO

   DO I = 1, M
#  if defined (WET_DRY)
     IF(ISWETN(I)*ISWETNT(I) == 1 )then
#  endif
       DO K = 1, KBM1
         F(I,K) = FF(I,K)
       END DO
#  if defined (WET_DRY)
     END IF
#  endif
   END DO

!   
!---- specify the source term------------------------------------------|
!
   IF(INTTIME.GE.DYESTART.AND.INTTIME.LE.DYESTOP) THEN

!    IF(SERIAL)THEN
!     DO KK=1, KSPE_DYE
!       K=K_SPECIFY(KK)
!       DO J=1,MSPE_DYE
!         I = M_SPECIFY(J)
!         F(I,K)= DYE_SOURCE_TERM
!       END DO
!     END DO
!    END IF 

!# if defined (MULTIPROCESSOR)
!    IF(PAR)THEN

     DO KK=1, KSPE_DYE
        K=K_SPECIFY(KK)
        DO J=1,MSPE_DYE
           I = NLID(M_SPECIFY(J))
           F(I,K)= DYE_SOURCE_TERM
        ENDDO
     ENDDO
!    END IF 
!#  endif
   ENDIF

   IF(DBG_SET(DBG_SBR)) WRITE(IPT,*)"End: vdif_dye"
   RETURN
   END SUBROUTINE VDIF_DYE
!==============================================================================|
!==============================================================================|
!  AVERAGE THE dye                                                     |
!==============================================================================|

   SUBROUTINE FCT_DYE

!==============================================================================|
   USE ALL_VARS
   USE MOD_UTILS
   USE BCS
   USE MOD_OBCS
   IMPLICIT NONE
   REAL(SP):: SMAX,SMIN
   INTEGER :: I,J,K
!==============================================================================|

  IF(HEATING_TYPE == 'body') RETURN

  IF(DBG_SET(DBG_SBR)) WRITE(IPT,*)"Start: fct_dye"

  nodes: DO I=1,M
     IF(IOBCN > 0)THEN
       DO J=1,IOBCN
         IF(I == I_OBC_N(J)) CYCLE nodes
       END DO
     END IF  	 

     IF(NUMQBC > 0)THEN
       DO J=1,NUMQBC
         IF(RIVER_INFLOW_LOCATION == 'node')THEN
	   IF(I == INODEQ(J)) CYCLE nodes
	 END IF  
         IF(RIVER_INFLOW_LOCATION == 'edge')THEN
	   IF(I == N_ICELLQ(J,1) .OR. I == N_ICELLQ(J,2)) CYCLE nodes
	 END IF  
       END DO
     END IF

     DO K=1,KBM1
       SMAX = MAXVAL(DYE(NBSN(I,1:NTSN(I)),K))
       SMIN = MINVAL(DYE(NBSN(I,1:NTSN(I)),K))

       IF(K == 1)THEN
         SMAX = MAX(SMAX,(DYE(I,K)*DZ(I,K+1)+DYE(I,K+1)*DZ(I,K))/  &
	        (DZ(I,K)+DZ(I,K+1)))
         SMIN = MIN(SMIN,(DYE(I,K)*DZ(I,K+1)+DYE(I,K+1)*DZ(I,K))/  &
	        (DZ(I,K)+DZ(I,K+1)))
       ELSE IF(K == KBM1)THEN
         SMAX = MAX(SMAX,(DYE(I,K)*DZ(I,K-1)+DYE(I,K-1)*DZ(I,K))/  &
	        (DZ(I,K)+DZ(I,K-1)))
         SMIN = MIN(SMIN,(DYE(I,K)*DZ(I,K-1)+DYE(I,K-1)*DZ(I,K))/  &
	        (DZ(I,K)+DZ(I,K-1)))
       ELSE
         SMAX = MAX(SMAX,(DYE(I,K)*DZ(I,K-1)+DYE(I,K-1)*DZ(I,K))/  &
	        (DZ(I,K)+DZ(I,K-1)), &
                 (DYE(I,K)*DZ(I,K+1)+DYE(I,K+1)*DZ(I,K))/   &
		 (DZ(I,K)+DZ(I,K+1)))
         SMIN = MIN(SMIN,(DYE(I,K)*DZ(I,K-1)+DYE(I,K-1)*DZ(I,K))/  &
	        (DZ(I,K)+DZ(I,K-1)), &
                 (DYE(I,K)*DZ(I,K+1)+DYE(I,K+1)*DZ(I,K))/   &
		 (DZ(I,K)+DZ(I,K+1)))
       END IF

       IF(SMIN-DYEF(I,K) > 0.0_SP)DYEF(I,K) = SMIN
       IF(DYEF(I,K)-SMAX > 0.0_SP)DYEF(I,K) = SMAX

     END DO
   END DO nodes

   IF(DBG_SET(DBG_SBR)) WRITE(IPT,*)"End: fct_dye"
   RETURN
   END SUBROUTINE FCT_DYE
!==============================================================================|


!==============================================================================|
!    Initialize dye fields (dye)                                               !
!    Calculate Mean Fields (DYEMEAN)                                           !
!==============================================================================|

!!$   SUBROUTINE INITIAL_DYE
!!$
!!$!==============================================================================!
!!$   USE ALL_VARS
!!$#  if defined (MULTIPROCESSOR)
!!$   USE MOD_PAR
!!$#  endif
!!$   IMPLICIT NONE
!!$   INTEGER :: I,K,ierr
!!$!   real(sp), allocatable :: temp_dye(:,:),temp_dyemean(:,:)
!!$
!!$!==============================================================================!
!!$
!!$   DO I=1,M
!!$     DO K=1,KB
!!$        DYE(I,K)=0.0_SP
!!$!        DYEMEAN(I,K)=0.0_SP
!!$     END DO
!!$   END DO
!!$
!!$
!!$      
!!$
!!$
!!$!! read global dye data into temporary arrays
!!$!   allocate(temp_dye(0:mgl,kb))
!!$!   allocate(temp_dyemean(0:mgl,kb))
!!$!   open(78,file='gom_bioini.dat')
!!$!   do i=1,MGL
!!$!      read(78,*)(temp_dye(i,k),k=1,kb)
!!$!   enddo
!!$!!   do i=1,MGL
!!$!!      read(78,*)(temp_dyemean(i,k),k=1,kb)
!!$!!   enddo
!!$!   temp_dyemean = temp_dye
!!$!   close(78)
!!$
!!$!! broadcast to all processors
!!$!   call mpi_bcast(temp_dye,kb*(mgl+1), mpi_f,0,MPI_FVCOM_GROUP,ierr)
!!$!   call mpi_bcast(temp_dyemean,kb*(mgl+1), mpi_f,0,MPI_FVCOM_GROUP,ierr)
!!$
!!$!! transform to local arrays 
!!$!   if(serial)then
!!$!     dye = temp_dye
!!$!     dyemean = temp_dyemean
!!$!   end if
!!$
!!$!# if defined (MULTIPROCESSOR)
!!$!   if(par)then
!!$!     do i=1,m
!!$!        dye(i,:) = temp_dye(ngid(i),:)
!!$!        dyemean(i,:) = temp_dyemean(ngid(i),:)
!!$!     end do
!!$!     if(par)call exchange(nc,mt,kb,myid,nprocs,dye,dyemean)
!!$!   end if
!!$
!!$!   deallocate(temp_dye)
!!$!   deallocate(temp_dyemean)
!!$!# endif
!!$    
!!$   RETURN
!!$   END SUBROUTINE INITIAL_DYE
!==============================================================================|
!==============================================================================|
!   Set Boundary Conditions on DYE                                             |
!==============================================================================|

   SUBROUTINE BCOND_DYE     

!------------------------------------------------------------------------------|
   USE ALL_VARS
   USE BCS
   USE MOD_OBCS
   IMPLICIT NONE
   REAL(SP) :: S2D,S2D_NEXT,S2D_OBC,T2D,T2D_NEXT,T2D_OBC,XFLUX2D,TMP
   INTEGER  :: I,J,K,J1,J11,J22
   REAL(SP) :: SMAX,SMIN
!------------------------------------------------------------------------------|


       
   IF(IOBCN > 0) THEN


!
!  SET dye CONDITIONS ON OUTER BOUNDARY
!
     DO I=1,IOBCN
       J=I_OBC_N(I)
       J1=NEXT_OBC(I)
       S2D=0.0_SP
       S2D_NEXT=0.0_SP
       XFLUX2D=0.0_SP
       DO K=1,KBM1
         S2D=S2D+DYE(J,K)*DZ(J,K)
         S2D_NEXT=S2D_NEXT+DYEF(J1,K)*DZ(J1,K)
         XFLUX2D=XFLUX2D+XFLUX_OBC(I,K)                 !*DZ(K)
       END DO
 
       IF(UARD_OBCN(I) > 0.0_SP) THEN
         TMP=XFLUX2D+S2D*UARD_OBCN(I)
         S2D_OBC=(S2D*DT(J)-TMP*DTI/ART1(J))/D(J)
         DO K=1,KBM1
!           DYEF(J,K)=S2D_OBC+(DYEF(J1,K)-S2D_NEXT)  !!bug 2 
           DYEF(J,K)=DYEF(J1,K)
          END DO

         DO K=1,KBM1
           SMAX = MAXVAL(DYE(NBSN(J,1:NTSN(J)),K))
           SMIN = MINVAL(DYE(NBSN(J,1:NTSN(J)),K))

           IF(K == 1)THEN
            SMAX = MAX(SMAX,(DYE(J,K)*DZ(J,K+1)+DYE(J,K+1)*DZ(J,K))/  &
                   (DZ(J,K)+DZ(J,K+1)))
            SMIN = MIN(SMIN,(DYE(J,K)*DZ(J,K+1)+DYE(J,K+1)*DZ(J,K))/  &
                   (DZ(J,K)+DZ(J,K+1)))
           ELSE IF(K == KBM1)THEN
            SMAX = MAX(SMAX,(DYE(J,K)*DZ(J,K-1)+DYE(J,K-1)*DZ(J,K))/  &
                   (DZ(J,K)+DZ(J,K-1)))
            SMIN = MIN(SMIN,(DYE(J,K)*DZ(J,K-1)+DYE(J,K-1)*DZ(J,K))/  &
                   (DZ(J,K)+DZ(J,K-1)))
           ELSE
            SMAX = MAX(SMAX,(DYE(J,K)*DZ(J,K-1)+DYE(J,K-1)*DZ(J,K))/  &
                   (DZ(J,K)+DZ(J,K-1)),                             &
                   (DYE(J,K)*DZ(J,K+1)+DYE(J,K+1)*DZ(J,K))/           &
                   (DZ(J,K)+DZ(J,K+1)))
            SMIN = MIN(SMIN,(DYE(J,K)*DZ(J,K-1)+DYE(J,K-1)*DZ(J,K))/  &
                   (DZ(J,K)+DZ(J,K-1)),                             &
                   (DYE(J,K)*DZ(J,K+1)+DYE(J,K+1)*DZ(J,K))/           &
                   (DZ(J,K)+DZ(J,K+1)))
           END IF

           IF(SMIN-DYEF(J,K) > 0.0_SP) DYEF(J,K) = SMIN
           IF(DYEF(J,K)-SMAX > 0.0_SP) DYEF(J,K) = SMAX

         END DO

        ELSE
         DO K=1,KBM1
           DYEF(J,K)=DYE(J,K)
         END DO
       END IF
     END DO
   END IF

!
!--SET BOUNDARY CONDITIONS-----------------------------------------------------|
!
!   DO K=1,KBM1
!     DYE(0,K)=0.0_SP
!   END DO

   RETURN
   END SUBROUTINE BCOND_DYE
!==============================================================================|
!
!==============================================================================|     
   SUBROUTINE NAME_LIST_INITIALIZE_DYE
   USE CONTROL
   
   IMPLICIT NONE
   
   !--Parameters in NameList NML_DYE_RELEASE
   DYE_ON          = .FALSE.
   DYE_RELEASE_START = 'Date or time to start dye release: Format the same as START_DATE'
   DYE_RELEASE_STOP  = 'Date or time to stop dye release: Format the same as START_DATE'
   KSPE_DYE        = 0
   MSPE_DYE        = 0
   K_SPECIFY       = 0
   M_SPECIFY       = 0
   DYE_SOURCE_TERM = 1.0

   RETURN
   END SUBROUTINE NAME_LIST_INITIALIZE_DYE  
!==============================================================================!
!
!==============================================================================!   
   
   SUBROUTINE SETUP_DYE
     implicit none
     integer status
     CHARACTER(LEN=4) :: BFLAG,EFLAG
     integer(itime) :: idstart, idstop

     SELECT CASE(USE_REAL_WORLD_TIME)
       CASE(.TRUE.)

          ! GET THE START TIME
          DYEstart = READ_DATETIME(DYE_RELEASE_START,DATE_FORMAT,TIMEZONE,status)
          if (status == 0) &
               & Call Fatal_Error("Could not read the date string DYE_RELEASE_START: "//trim(DYE_RELEASE_START))
     
     ! GET THE START TIME
          DYEstop = READ_DATETIME(DYE_RELEASE_STOP,DATE_FORMAT,TIMEZONE,status)
          if (status == 0) &
               & Call Fatal_Error("Could not read the date string DYE_RELEASE_STOP: "//trim(DYE_RELEASE_STOP))
     
       CASE (.FALSE.)

          
          ! GET THE START AND END INFORMATION
          CALL IDEAL_TIME_STRING2TIME(DYE_RELEASE_START,BFLAG,DYEstart,Idstart)
          CALL IDEAL_TIME_STRING2TIME(DYE_RELEASE_STOP,EFLAG,DYEstop,Idstop)
          
          ! SANITY CHECK
          IF (BFLAG /= EFLAG) CALL FATAL_ERROR&
               ('IDEALIZED DYE MODEL RELEASE TIME SPECIFICATION IS INCORRENT',&
               &'START AND STOP CAN BE IN EITHER CYCLES OR TIME BUT NOT MIXED',&
               & trim(DYE_RELEASE_START),trim(DYE_RELEASE_STOP) )
          
          IF (BFLAG == 'time') THEN ! IF START AND END TIME WERE SPECIFIED

             ! NOTHING TO BE DONE
             
          ELSE IF(BFLAG == 'step') THEN ! IF START AND END IINT WERE SPECIFIED

             ! Get the start time for the dye from the cycle number
             DYEstart = StartTime + IMDTI *(Idstart-Istart)

             ! Get the start time for the dye from the cycle number
             DYEstart = StartTime + IMDTI *(Idstop-Istart)

          ELSE
             CALL FATAL_ERROR('IDEAL_TIME_STRING2TIME returned invalid flag')

          END IF


          
          if(DBG_SET(dbg_log)) then 
             call print_time(DYESTART,ipt,"dyestart")
             call print_time(DYESTOP,ipt,"dyestop")
          end if



          if (DYESTART > DYESTOP)&
               & Call Fatal_Error("The dye release start time must be &
               &less than or equal to the dye release stop time")


          if (STARTTIME > DYESTART)&
               & Call Fatal_Error("The dye release start time must be &
               &greater than or equal to the model start time")
          
          if (ENDTIME < DYESTOP)&
               & Call Fatal_Error("The dye release stop time must be l&
               &ess than or equal to the model end time")


       END SELECT


   END SUBROUTINE SETUP_DYE
!==============================================================================!
!
!==============================================================================!     
   SUBROUTINE NAME_LIST_PRINT_DYE
   USE CONTROL
   
   IMPLICIT NONE
   
   write(UNIT=IPT,NML=NML_DYE_RELEASE)
     
   RETURN
   END SUBROUTINE NAME_LIST_PRINT_DYE
!==============================================================================!   
END MODULE MOD_DYE
# else
   SUBROUTINE SETUP_DYE
     implicit none
     return
   END SUBROUTINE SETUP_DYE
# endif

