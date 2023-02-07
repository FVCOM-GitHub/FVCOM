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

!==============================================================================!
! TYPE_OBC = 1: Surface Elevation Specified (Tidal Forcing) (ASL)              !
! TYPE_OBC = 2: AS TYPE_OBC=1 AND NON-LINEAR FLUX FOR CURRENT AT OPEN BOUNDARY !
! TYPE_OBC = 3: Zero Surface Elevation Boundary Condition (ASL_CLP)            !
! TYPE_OBC = 4: AS TYPE_OBC=3 AND NON-LINEAR FLUX FOR CURRENT AT OPEN BOUNDARY !
! TYPE_OBC = 5: GRAVITY-WAVE RADIATION IMPLICIT OPEN BOUNDARY CONDITION (GWI)  !
! TYPE_OBC = 6: AS TYPE_OBC=5 AND NON-LINEAR FLUX FOR CURRENT AT OPEN BOUNDARY !
! TYPE_OBC = 7: BLUMBERG AND KHANTA IMPLICIT OPEN BOUNDARY CONDITION (BKI)     !
! TYPE_OBC = 8: AS TYPE_OBC=7 AND NON-LINEAR FLUX FOR CURRENT AT OPEN BOUNDARY !
! TYPE_OBC = 9: ORLANSKI RADIATION EXPLICIT OPEN BOUNDARY CONDITION (ORE)      !
! TYPE_OBC =10: AS TYPE_OBC=9 AND NON-LINEAR FLUX FOR CURRENT AT OPEN BOUNDARY !
!                                                                              !
! TYPE_TSOBC = 1: THE PERTURBATION OF TEMPERATURE AND SALINITY AT OPEN BOUNDARY!
!                 ARE EQUAL TO THAT AT NEXT_OBC                                !
! TYPE_TSOBC = 2: GRAVITY-WAVE RADIATION IMPLICIT OPEN BOUNDARY CONDITION FOR  !
!                 THE PERTURBATION OF TEMPERATURE AND SALINITY                 !
! TYPE_TSOBC = 3: BLUMBERG AND KHANTA IMPLICIT OPEN BOUNDARY CONDITION FOR     !
!                 THE PERTURBATION OF TEMPERATURE AND SALINITY                 !
! TYPE_TSOBC = 4: ORLANSKI RADIATION EXPLICIT OPEN BOUNDARY CONDITION FOR      !
!                 THE PERTURBATION OF TEMPERATURE AND SALINITY                 !
!==============================================================================!

MODULE MOD_OBCS
  USE MOD_PREC
  USE MOD_UTILS
  USE CONTROL,  type_tsobc => obc_ts_type
  USE MOD_TIME
  USE MOD_FORCE
  USE BCS
  IMPLICIT NONE
  SAVE

  !--Open Boundary Types, Lists, Pointers
# if defined (TIDE_OUTPUT)
  INTEGER               :: IOBCNODE_GL      !! = ntidenode_GL
  INTEGER,  ALLOCATABLE :: I_OBCNODE_GL(:)  !! = I_TIDENODE_GL
  INTEGER               :: IOBCELL_GL       !! = ntidecell_GL
  INTEGER,  ALLOCATABLE :: I_OBCELL_GL(:)   !! = I_TIDECELL_GL
# endif

  !   INTEGER               :: IOBCN_GL         !!GLOBAL NUMBER OF OPEN BOUNDARY NODES
  !   INTEGER               :: IOBCN            !!LOCAL NUMBER OF OPEN BOUNDARY NODES
  !   INTEGER,  ALLOCATABLE :: I_OBC_GL(:)      !!GLOBAL ID OF OPEN BOUNDARY NODES
  !   INTEGER,  ALLOCATABLE :: I_OBC_N(:)       !!OPEN BOUNDARY NODE LIST
  INTEGER,  ALLOCATABLE :: NEXT_OBC(:)      !!INTERIOR NEIGHBOR OF OPEN BOUNDARY NODE
  INTEGER,  ALLOCATABLE :: NEXT_OBC2(:)     !!INTERIOR NEIGHBOR OF NEXT_OBC
  !   INTEGER,  ALLOCATABLE :: TYPE_OBC(:)      !!OUTER BOUNDARY NODE TYPE (FOR SURFACE ELEVATION)
  !   INTEGER,  ALLOCATABLE :: TYPE_OBC_GL(:)   !!OUTER BOUNDARY NODE TYPE (FOR SURFACE ELEVATION)
  INTEGER               :: IBCN(5)          !!NUMBER OF EACH TYPE OF OBN IN LOCAL  DOM
  INTEGER               :: IBCN_GL(5)       !!NUMBER OF EACH TYPE OF OBN IN GLOBAL DOM
  INTEGER,  ALLOCATABLE :: OBC_LST(:,:)     !!MAPPING OF OPEN BOUNDARY ARRAYS TO EACH TYPE
  INTEGER,  ALLOCATABLE :: NADJN_OBC(:)     !!NUMBER OF ADJACENT OPEN BOUNDARY NODES TO OBN
  INTEGER,  ALLOCATABLE :: ADJN_OBC(:,:)    !!ADJACENT OBN's of OBN
  INTEGER,  ALLOCATABLE :: NADJC_OBC(:)     !!NUMBER OF ADJACENT OPEN BOUNDARY CELLS TO OBN
  INTEGER,  ALLOCATABLE :: ADJC_OBC(:,:)    !!ADJACENT OPEN BOUNDARY CELLS

  !-- Open Boundary TEMPERATURE AND SALINITY FOR NUDGING
  REAL(SP), ALLOCATABLE :: TEMP_OBC(:,:)      !!OPEN BOUNDARY TEMPERATURE (see mod_force.F)     
  REAL(SP), ALLOCATABLE :: SALT_OBC(:,:)      !!OPEN BOUNDARY SALT (see mod_force.F)   
# if defined (WATER_QUALITY)
  REAL(SP), ALLOCATABLE :: WQM_OBC(:,:,:)      !!OPEN BOUNDARY WQM (see mod_force.F)   
# endif
# if defined (BioGen)
  REAL(SP), ALLOCATABLE :: BIO_OBC(:,:,:)      !!OPEN BOUNDARY BIO (see mod_force.F)   
# endif
  !--Open Boundary Metrics
  INTEGER,  ALLOCATABLE :: NFLUXF_OBC(:)    !!NUMBER OF FLUX SEGMENTS TO OBN
  REAL(SP), ALLOCATABLE :: FLUXF_OBC(:,:)   !!FLUX FRACTION ON EACH SIDE OF OBN
  REAL(SP), ALLOCATABLE :: NXOBC(:)         !!INWARD POINTING X-NORMAL OF OBN
  REAL(SP), ALLOCATABLE :: NYOBC(:)         !!INWARD POINTING Y-NORMAL OF OBN
  REAL(SP), ALLOCATABLE :: DLTN_OBC(:)      !!DISTANCE BETWEEN NEXT_OBC AND OBN NORMAL TO BOUNDARY

  !--Previous Time Level Free Surface Fields
  REAL(SP), ALLOCATABLE :: ELM1(:)          !!SURFACE ELEV FROM PREVIOUS TIME LEVEL (ORLANSKI COND)
  REAL(SP), ALLOCATABLE :: ELM2(:)          !!SURFACE ELEV FROM PREVIOUS TIME LEVEL (ORLANSKI COND)
  REAL(SP), ALLOCATABLE :: T1M1(:,:)          !!TEMPERATURE FROM PREVIOUS TIME LEVEL (ORLANSKI COND)
  REAL(SP), ALLOCATABLE :: T1M2(:,:)          !!TEMPERATURE FROM PREVIOUS TIME LEVEL (ORLANSKI COND)
  REAL(SP), ALLOCATABLE :: S1M1(:,:)          !!SALINITY FROM PREVIOUS TIME LEVEL (ORLANSKI COND)
  REAL(SP), ALLOCATABLE :: S1M2(:,:)          !!SALINITY FROM PREVIOUS TIME LEVEL (ORLANSKI COND)

  !--Nonlinear Velocity Open Boundary Condition Arrays
  REAL(SP), ALLOCATABLE :: FLUXOBN(:)
  REAL(SP), ALLOCATABLE :: IUCP(:)
  REAL(SP), ALLOCATABLE :: XFLUX_OBCN(:)
  REAL(SP), ALLOCATABLE :: UARD_OBCN(:)
  REAL(SP), ALLOCATABLE :: XFLUX_OBC(:,:)
# if defined (ICE)
! - begin afm 20171113 for ice open boundary -------------
  REAL(SP), ALLOCATABLE :: UIARD_OBCN(:)
  REAL(SP), ALLOCATABLE :: XFLUX_ICE_OBC(:)
! - end afm 20171113 for ice open boundary -------------
# endif

CONTAINS
  !==========================================================================|
  SUBROUTINE ALLOC_OBC_DATA

    !--------------------------------------------------------------------------|
    !  ALLOCATE AND INITIALIZE SURFACE ELEVATION ARRAYS FOR                    |
    !  TIME STEPS N-1 AND N-2                                                  |
    !--------------------------------------------------------------------------|

    USE ALL_VARS
#   if defined (BioGen)
    USE MOD_1D, ONLY : NTT
#   endif    
    IMPLICIT NONE

    ALLOCATE(ELM1(0:MT))             ;ELM1       = ZERO 
    ALLOCATE(ELM2(0:MT))             ;ELM2       = ZERO 
    ALLOCATE(NEXT_OBC(IOBCN))        ;NEXT_OBC   = 0
    ALLOCATE(NEXT_OBC2(IOBCN))       ;NEXT_OBC2  = 0
    ALLOCATE(FLUXOBN(1:NT))          ;FLUXOBN    = ZERO
    ALLOCATE(IUCP(0:NT))             ;IUCP       = 1
    ALLOCATE(XFLUX_OBCN(IOBCN+1))    ;XFLUX_OBCN = ZERO
    ALLOCATE(UARD_OBCN(IOBCN+1))     ;UARD_OBCN  = ZERO
    ALLOCATE(XFLUX_OBC(IOBCN+1,KBM1));XFLUX_OBC  = ZERO
#   if defined (ICE)
! - begin afm 20171113 for ice open boundary -------------
    ALLOCATE(UIARD_OBCN(IOBCN+1))     ;UIARD_OBCN  = ZERO
    ALLOCATE(XFLUX_ICE_OBC(IOBCN+1))     ;XFLUX_ICE_OBC  = ZERO
! - end afm 20171113 for ice open boundary -------------
#   endif

    ALLOCATE(TEMP_OBC(IOBCN,KBM1))   ;TEMP_OBC   = ZERO
    ALLOCATE(SALT_OBC(IOBCN,KBM1))   ;SALT_OBC   = ZERO
#   if defined (WATER_QUALITY)
    ALLOCATE(WQM_OBC(IOBCN,KBM1,NB))   ;WQM_OBC   = ZERO
#   endif
#   if defined (BioGen)
    IF(BIOLOGICAL_MODEL)THEN
      ALLOCATE(BIO_OBC(IOBCN,KBM1,NTT));BIO_OBC   = ZERO
    END IF
#   endif

    IF(TYPE_TSOBC == 4)THEN
       ALLOCATE(T1M1(0:MT,KBM1))        ;T1M1       = ZERO 
       ALLOCATE(T1M2(0:MT,KBM1))        ;T1M2       = ZERO 
       ALLOCATE(S1M1(0:MT,KBM1))        ;S1M1       = ZERO 
       ALLOCATE(S1M2(0:MT,KBM1))        ;S1M2       = ZERO 
    END IF
    RETURN
  END SUBROUTINE ALLOC_OBC_DATA
  !==========================================================================|

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|

  !==========================================================================|
  SUBROUTINE ASSIGN_ELM1_TO_ELM2

    !--------------------------------------------------------------------------|
    !  Assign ELM1 to ELM2 and EL to ELM1                                      |
    !--------------------------------------------------------------------------|

    USE ALL_VARS
    IMPLICIT NONE

    ELM2 = ELM1
    ELM1 = EL

    IF(TYPE_TSOBC == 4)THEN
       T1M2 = T1M1
       T1M1 = T1

       S1M2 = S1M1
       S1M1 = S1
    END IF

    RETURN
  END SUBROUTINE ASSIGN_ELM1_TO_ELM2
  !==========================================================================|

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|

# if defined (EQUI_TIDE)
  !==========================================================================|
  SUBROUTINE ELEVATION_EQUI ! ONLY USED if defined (EQUI_TIDE)
    !--------------------------------------------------------------------------|
    !  Surface Elevation of EQUILIBRIUM TIDE                                   |
    !--------------------------------------------------------------------------|
    USE ALL_VARS
    IMPLICIT NONE

    TYPE(TIME):: TIME_ELAPSED
    INTEGER :: I,J,K
!    REAL(SP), ALLOCATABLE :: FORCE(:)
    REAL(DP):: TIME_FLT,PHAI_IJ

    Integer :: year_cur,month_cur,mdays_cur,ic1,kdm,kd0,intdays,equi_iv,equi_l,equi_iuu
    real*8  :: equi_d1,equi_h,equi_pp,equi_s,equi_p,equi_enp,equi_dh,equi_dpp,equi_ds,equi_dp,equi_dnp
    real*8  :: hourm,dtau,equi_hh,equi_tau,vdbl,slat
    real*8  :: equi_rr,equi_uudbl,equi_uu,sumc,sums
    type(time) :: time_adjust

# if !defined(SPHERICAL)
    IF(.not.USE_PROJ) CALL FATAL_ERROR &
         & ("ELEVATION_EQI: ILLEGAL ATTEMPT TO USE EQUILIBRIUM TIDE IN CARTESIAN MODE",&
         & "YOU MUST BE USING PROJ4 TO DO THIS",&
         & "THE POLICE HAVE BEEN NOTIFIED AND WILL BE THERE SHORTLY")
# endif


# if !defined (SEMI_IMPLICIT)

    ! THIS CODE IS CALLED DURING EACH RK STAGE
    ! DECIDE WHICH TIME YOU WANT TO USE!
    TIME_ELAPSED = RKTime - SpecTime
# else
   TIME_ELAPSED = IntTime - SpecTime
# endif

!! IN MY JUDGEMENT, A TIMESERIES FORCING DOES NOT MAKE SENSE FOR ATMO/EQUITIDE


!    SELECT CASE(TIDE_FORCING_TYPE)
!    CASE(TIDE_FORCING_TIMESERIES)
       !
       !-Julian: Set Elevation Based on Linear Interpolation Between Two Data Times-|
       !
!       CALL FATAL_ERROR("ELEVATION_EQUI: Not set up for Julian Timeseries forcing yet?")

!    CASE(TIDE_FORCING_SPECTRAL)
       !
       !-Non-Julian: Set Elevation of Equilibrium Tide -----------------------------|
       !

       ELF_EQI=0.0_SP

       IF(USE_REAL_WORLD_TIME) THEN  

         time_adjust%mjd = 678942
         time_adjust%musod = 0

         CALL Now_2_month_days(IntTime,year_cur,month_cur,mdays_cur)
         IF(month_equi/=month_cur) THEN

           ic1 = year_cur/100
           year_cur = year_cur - ic1*100
           CALL GDAY1(16,month_cur,year_cur,ic1,kdm)

           TIME_16DAY%MJD = kdm
           TIME_16DAY%MuSOD = 0
           month_equi = month_cur

           hourm = dfloat(kdm)*24.0d0
           equi_d1 = hourm/24.0d0
           CALL GDAY1(31,12,99,18,kd0)
           equi_d1 = equi_d1-dfloat(kd0)-0.5d0
           call astro(equi_d1,equi_h,equi_pp,equi_s,equi_p,equi_enp,equi_dh,equi_dpp,equi_ds,equi_dp,equi_dnp)
           intdays=idint(hourm/24.0d0)
           equi_hh=hourm-dfloat(intdays*24)
           equi_tau=equi_hh/24.0d0+equi_h-equi_s
           dtau=365.0d0+equi_dh-equi_ds

           do j=1, nTideComps
             vdbl=equi_astr(equi_id(j))%ii*equi_tau+equi_astr(equi_id(j))%jj*equi_s+equi_astr(equi_id(j))%kk*equi_h  &
                 +equi_astr(equi_id(j))%ll*equi_p+equi_astr(equi_id(j))%mm*equi_enp+equi_astr(equi_id(j))%nn*equi_pp &
                 +equi_astr(equi_id(j))%semi
             equi_iv = int(vdbl)
             equi_iv = (equi_iv/2)*2
             EQUI_VN(J) = vdbl - equi_iv

             do i=1, mt

               slat = vy(i)
               if(abs(slat).lt.5.0) slat=sign(5.0,slat) 
               slat=dsin(pi*slat/180.0d0)

               sumc = 1.0D0
               sums = 0.0D0
               do k=1, equi_astr(equi_id(j))%nj
                 equi_rr = equi_astr(equi_id(j))%ee(k) 
                 equi_l  = equi_astr(equi_id(j))%ir(k)+1                    
                 goto (901, 902, 903) equi_l
  902            equi_rr = equi_astr(equi_id(j))%ee(k)*0.36309D0*(1.0D0-5.0D0*slat*slat)/slat
                 goto 901
  903            equi_rr = equi_astr(equi_id(j))%ee(k)*2.59808D0*slat
  901            continue
                 equi_uudbl = equi_astr(equi_id(j))%ldel(k)*equi_p+equi_astr(equi_id(j))%mdel(k)*equi_enp+equi_astr(equi_id(j))%ndel(k)*equi_pp+equi_astr(equi_id(j))%ph(k)
                 equi_iuu = int(equi_uudbl)
                 equi_uu  = equi_uudbl - equi_iuu
                 
                 sumc = sumc + equi_rr*dcos(equi_uu*pi2)
                 sums = sums + equi_rr*dsin(equi_uu*pi2)
               enddo

               EQUI_F(I,J)=dsqrt(sumc*sumc+sums*sums)
               EQUI_U(I,J)=datan2(sums,sumc)/pi2 
             enddo

           enddo

         ENDIF

         DO J = 1,nTideComps

            TIME_FLT = SECONDS( (TIME_ELAPSED-(TIME_16DAY-SpecTime-time_adjust)) * PI2/PERIOD(J))
            SELECT CASE(TIDE_TYPE(J))
            CASE(SEMIDIURNAL)

               DO I = 1, MT
                  PHAI_IJ = LON(I)*PI2/360.0_SP+(EQUI_VN(J)+EQUI_U(I,J))*PI2/2.0_DP
                  ELF_EQI(I) = ELF_EQI(I) + EQUI_F(I,J)*BETA_EQI(J)*APT_EQI(J)*COS(LAT(I)*PI2/360.0_SP)**2         &
                       *DCOS(TIME_FLT+2.0_DP*PHAI_IJ)
               END DO
            CASE(DIURNAL)

               DO I = 1, MT
                  PHAI_IJ = LON(I)*PI2/360.0_SP+(EQUI_VN(J)+EQUI_U(I,J))*PI2
                  ELF_EQI(I) = ELF_EQI(I) + EQUI_F(I,J)*BETA_EQI(J)*APT_EQI(J)*SIN(2.0_SP*LAT(I)*PI2/360.0_SP)     &
                       *DCOS(TIME_FLT+PHAI_IJ)
               END DO

            CASE DEFAULT
               CALL FATAL_ERROR("ELEVATION_EQUI: INVALID TIDE_TYPE:",&
                    &"Choices:"//TRIM(DIURNAL)//" or "//TRIM(SEMIDIURNAL) )
            END SELECT
         END DO         

       ELSE

         DO J = 1,nTideComps

            TIME_FLT = SECONDS(TIME_ELAPSED * PI2/PERIOD(J))             
            SELECT CASE(TIDE_TYPE(J))
            CASE(SEMIDIURNAL)

               DO I = 1, MT
                  PHAI_IJ = LON(I)*PI2/360.0_SP
                  ELF_EQI(I) = ELF_EQI(I) + BETA_EQI(J)*APT_EQI(J)*COS(LAT(I)*PI2/360.0_SP)**2         &
                       *DCOS(TIME_FLT+2.0_DP*PHAI_IJ)
               END DO
            CASE(DIURNAL)
             
               DO I = 1, MT   
                  PHAI_IJ = LON(I)*PI2/360.0_SP
                  ELF_EQI(I) = ELF_EQI(I) + BETA_EQI(J)*APT_EQI(J)*SIN(2.0_SP*LAT(I)*PI2/360.0_SP)     &
                       *DCOS(TIME_FLT+PHAI_IJ)
               END DO
             
            CASE DEFAULT
               CALL FATAL_ERROR("ELEVATION_EQUI: INVALID TIDE_TYPE:",&
                    &"Choices:"//TRIM(DIURNAL)//" or "//TRIM(SEMIDIURNAL) )
            END SELECT
         END DO

       ENDIF

       ELF_EQI = ELF_EQI * RAMP
       
!    CASE DEFAULT
!       CALL FATAL_ERROR("ELEVATION_EQUI: INVALID TIDE FORCING TYPE ?")
!    END SELECT

    RETURN
  END SUBROUTINE ELEVATION_EQUI
  !==========================================================================|
# endif

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|

  !==========================================================================|
  SUBROUTINE ELEVATION_ATMO ! ONLY USED if defined (ATMO_TIDE)
    !--------------------------------------------------------------------------|
    !  Surface Elevation of ATMOSPHERIC TIDE                                   |
    !--------------------------------------------------------------------------|
    USE ALL_VARS
    IMPLICIT NONE

    REAL(SP),PARAMETER :: APT_ATMO = 0.0113_SP       
    REAL(DP),PARAMETER :: ALFA_ATMO = 112.0_SP*PI2/360.0_SP
    REAL(SP),PARAMETER :: S2PERIOD = 43200.0_SP

    INTEGER :: I,J
    TYPE(TIME):: TIME_ELAPSED
    REAL(DP):: TIME_FLT,PHAI_IJ
    REAL(SP):: FORCE

# if !defined(SPHERICAL)
    IF(.not.USE_PROJ) CALL FATAL_ERROR &
         & ("ELEVATION_ATMO: ILLEGAL ATTEMPT TO USE ATMOSPHERIC TIDE IN CARTESIAN MODE",&
         & "YOU MUST BE USING PROJ4 TO DO THIS",&
         & "THE POLICE HAVE BEEN NOTIFIED AND WILL BE THERE SHORTLY")
# endif



# if !defined (SEMI_IMPLICIT)

    ! THIS CODE IS CALLED DURING EACH RK STAGE
    ! DECIDE WHICH TIME YOU WANT TO USE!
    TIME_ELAPSED = RKTime - SpecTime
# else
   TIME_ELAPSED = IntTime - SpecTime
# endif

!! IN MY JUDGEMENT, A TIMESERIES FORCING DOES NOT MAKE SENSE FOR ATMO/EQUITIDE

!    SELECT CASE(TIDE_FORCING_TYPE)
!    CASE(TIDE_FORCING_TIMESERIES)
       !
       !-Julian: Set Elevation Based on Linear Interpolation Between Two Data Times-|
       !
!       CALL FATAL_ERROR("ELEVATION_EQUI: Not set up for Julian Timeseries forcing yet?")

!    CASE(TIDE_FORCING_SPECTRAL)
       !
       !-Non-Julian: Set Elevation of Equilibrium Tide -----------------------------|
       !

       TIME_FLT = SECONDS(TIME_ELAPSED * PI2/S2PERIOD)

       DO I = 1, MT
          PHAI_IJ = LON(I)*PI2/360.0_SP
          FORCE = APT_ATMO*DCOS(TIME_FLT+2.0_DP*PHAI_IJ-ALFA_ATMO)
          ELF_ATMO(I) = FORCE * RAMP
       END DO

!    CASE DEFAULT
!       CALL FATAL_ERROR("ELEVATION_ATMO: INVALID TIDE FORCING TYPE ?")
!    END SELECT



    RETURN
  END SUBROUTINE ELEVATION_ATMO
  !==========================================================================|

!==============================================================================|
!     SET BOUNDARY CONDITIONS FOR CALCULATING atmospher pressure external model|
!                                                                              |
!==============================================================================|

   SUBROUTINE BCOND_PA_AIR
#  if defined (AIR_PRESSURE)
!==============================================================================|
   USE ALL_VARS
#  if defined (MULTIPROCESSOR)
   USE MOD_PAR, ONLY : NC,AEXCHANGE
#  endif
   IMPLICIT NONE

   TYPE(TIME) :: NOW
   REAL(SP) :: FACT,UFACT,TX,TY
!    for ice model
   REAL(SP) :: pa_air2
   INTEGER  I,J,L1,L2,IERR

# if defined (SEMI_IMPLICIT)
   NOW = IntTime
# else
   NOW = ExtTime
#endif

!==============================================================================|
!   SURFACE AIR PRESSURE FORCING FOR INTERNAL MODE                                          !
!==============================================================================|
   CALL UPDATE_AIRPRESSURE(NOW,PA_AIR)

#  if defined (MULTIPROCESSOR)
   IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,PA_AIR)         
#  endif

!!$   ELF_AIR = -(PA_AIR-SLP0)/(1025._SP*GRAV_N)
   ELF_AIR(1:MT) = -(PA_AIR(1:MT)-SLP0)/(1025._SP*GRAV_N(1:MT))

 
   RETURN
#  endif
   END SUBROUTINE BCOND_PA_AIR
   
   
!==============================================================================|
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|

  !==========================================================================|
  SUBROUTINE BCOND_ASL

    !--------------------------------------------------------------------------|
    !  Surface Elevation Boundary Condition (Tidal Forcing)                    |
    !--------------------------------------------------------------------------|

    USE ALL_VARS, only : ELF
    IMPLICIT NONE

    INTEGER :: I,II,J,JN
    TYPE(TIME):: TIME_ELAPSED
    REAL(DP):: TIME_FLT,PHAI_IJ
    REAL(SP):: FORCE
    REAL(SP), ALLOCATABLE :: BND_ELV(:)


    ALLOCATE(BND_ELV(IOBCN))
    !
    !-Julian: Set Elevation Based on Linear Interpolation Between Two Data Times-|
    !
    SELECT CASE(TIDE_FORCING_TYPE)
    CASE(TIDE_FORCING_TIMESERIES)

# if !defined (SEMI_IMPLICIT)

       ! THIS CODE IS CALLED DURING EACH RK STAGE
       ! DECIDE WHICH TIME YOU WANT TO USE!
       CALL UPDATE_TIDE(RKTime,BND_ELV)
       !      CALL UPDATE_TIDE(ExtTime,BND_ELV)
# else
       CALL UPDATE_TIDE(IntTime,BND_ELV)
# endif


       DO J = 1, IBCN(1)
          JN = OBC_LST(1,J)
          II = I_OBC_N(JN)

          ELF(II) = BND_ELV(J)*RAMP
       END DO

       !
       !-Non-Julian: Set Elevation Based on Input Amplitude and Phase of Tidal Comps-|
       !

    CASE(TIDE_FORCING_SPECTRAL)



# if !defined (SEMI_IMPLICIT)

       ! THIS CODE IS CALLED DURING EACH RK STAGE
       ! DECIDE WHICH TIME YOU WANT TO USE!
       TIME_ELAPSED = RKTime - SpecTime
       !    TIME_ELAPSED = ExtTime - SpecTime
# else
       TIME_ELAPSED = IntTime - SpecTime
# endif

       !      write(ipt,*) "-----------------------------------------"
       DO I = 1, IBCN(1)
          JN = OBC_LST(1,I)
          II = I_OBC_N(JN)
          FORCE = 0.0_SP
          DO J = 1,nTideComps
             PHAI_IJ = PHAI(I,J)*PI2/360.0_SP
             TIME_FLT = SECONDS(TIME_ELAPSED * PI2/PERIOD(J))
             ! Take cosine of a double precision number to preserve accuracy!
             FORCE = APT(I,J)*DCOS(TIME_FLT -PHAI_IJ) + FORCE
          END DO
          FORCE = FORCE + EMEAN(I)
          ELF(II) = FORCE * RAMP
          !         write(ipt,*) "ELF=",ELF(II)

       END DO

    CASE DEFAULT
       CALL FATAL_ERROR("BCOND_ASL: INVALID TIDE FORCING TYPE ?")
    END SELECT

    DEALLOCATE(BND_ELV)

    RETURN
  END SUBROUTINE BCOND_ASL
!!$!==========================================================================|

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|

  !==========================================================================|
  SUBROUTINE BCOND_ASL_CLP

    !--------------------------------------------------------------------------|
    !  Zero Surface Elevation Boundary Condition                               |
    !--------------------------------------------------------------------------|

    USE ALL_VARS
    IMPLICIT NONE

    INTEGER :: I,II,JN

    DO I = 1, IBCN(2)
       JN = OBC_LST(2,I)
       II = I_OBC_N(JN)  
       ELF(II) = 0.0_SP
    END DO

    RETURN
  END SUBROUTINE BCOND_ASL_CLP
  !==========================================================================|

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|

  !==========================================================================|
  SUBROUTINE BCOND_GWI(K_RK)

    !--------------------------------------------------------------------------|
    !  GRAVITY-WAVE RADIATION IMPLICIT OPEN BOUNDARY CONDITION (GWI)           |
    !--------------------------------------------------------------------------|

    USE ALL_VARS
    IMPLICIT NONE

    INTEGER :: I1,I2,J,JN,K_RK
    REAL(SP):: CC,CP,ALPHA_RK_TMP

    IF(K_RK == 0)THEN
      ALPHA_RK_TMP = 1.0
    ELSE
      ALPHA_RK_TMP = ALPHA_RK(K_RK)
    END IF

    DO J = 1,IBCN(3)
      JN = OBC_LST(3,J)
      I1 =  I_OBC_N(JN) 
      I2 = NEXT_OBC(JN) 
#     if defined (SEMI_IMPLICIT)
      CC = SQRT(GRAV_N(I1)*H(I1))*DTE/DLTN_OBC(JN)
      CP = CC + 1.0_SP
      ELF(I1) = (CC*ELF(I2) + EL(I1))/CP
#     else
      CC = SQRT(GRAV_N(I1)*H(I1))*DTE/DLTN_OBC(JN)*ALPHA_RK_TMP
      CP = CC + 1.0_SP
      ELF(I1) = (CC*ELF(I2) + ELRK(I1))/CP
#     endif     
    END DO

    RETURN
  END SUBROUTINE BCOND_GWI
  !==========================================================================|

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|

  !==========================================================================|
  SUBROUTINE BCOND_BKI(K_RK)

    !--------------------------------------------------------------------------|
    !  BLUMBERG AND KHANTA IMPLICIT OPEN BOUNDARY CONDITION                    |
    !--------------------------------------------------------------------------|

    USE ALL_VARS
    IMPLICIT NONE

    INTEGER :: I1,I2,J,JN,K_RK
    REAL(SP):: CC,CP,ALPHA_RK_TMP

    IF(K_RK == 0)THEN
      ALPHA_RK_TMP = 1.0
    ELSE
      ALPHA_RK_TMP = ALPHA_RK(K_RK)
    END IF

    DO J = 1,IBCN(4)
      JN = OBC_LST(4,J)
      I1 = I_OBC_N(JN)
      I2 = NEXT_OBC(JN)
#     if defined (SEMI_IMPLICIT)
      CC = SQRT(GRAV_N(I1)*H(I1))*DTE/DLTN_OBC(JN)
      CP = CC + 1.0_SP
      ELF(I1) = (CC*ELF(I2) + EL(I1)*(1.0_SP-DTE/10800.0_SP))/CP
#     else
      CC = SQRT(GRAV_N(I1)*H(I1))*DTE/DLTN_OBC(JN)*ALPHA_RK_TMP
      CP = CC + 1.0_SP
      ELF(I1) = (CC*ELF(I2) + ELRK(I1)*(1.0_SP-DTE/10800.0_SP))/CP
#     endif     
    END DO

    RETURN
  END SUBROUTINE BCOND_BKI
  !==============================================================================|

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|

  !==============================================================================|
  SUBROUTINE BCOND_ORE

    !------------------------------------------------------------------------------|
    !  ORLANSKI RADIATION EXPLICIT OPEN BOUNDARY CONDITION (ORE)                   |
    !------------------------------------------------------------------------------|

    USE ALL_VARS
    IMPLICIT NONE

    INTEGER  :: I1,I2,I3,J,JN
    REAL(SP) :: CL, MU

    DO J = 1, IBCN(5)
       JN = OBC_LST(5,J)
       I1 = I_OBC_N(JN)
       I2 = NEXT_OBC(JN)
       I3 = NEXT_OBC2(JN)

       CL = (ELM2(I2)-EL(I2))/(EL(I2)+ELM2(I2)-2.0*ELM1(I3))
       IF(CL >= 1.0)THEN
          MU = 1.0
       ELSE IF(CL > 0.0 .AND. CL < 1.0)THEN
          MU = CL
       ELSE
          MU = 0.0
       END IF

       ELF(I1)=(ELM1(I1)*(1.0-MU)+2.0*MU*EL(I2))/(1.0+MU)
    END DO

    RETURN
  END SUBROUTINE BCOND_ORE
  !==============================================================================|

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|

  !========================================================================

  SUBROUTINE SETUP_OBCTYPES
    USE LIMS
    IMPLICIT NONE

    INTEGER              :: I,I1,I2,NCNT,IERR,J
    INTEGER, ALLOCATABLE :: TEMP1(:),TEMP2(:),TEMP3(:),TEMP4(:),TEMP5(:),TEMP6(:),TEMP7(:),ITEMP(:)

    ! CREATE THE BC TYPE INTEGERS USED IN THE MODEL
    CALL SEPARATE_OBC ! (MOD_OBCS.F)

    !==============================================================================|
    !   REPORT AND CHECK OBC RESULTS                                               |
    !==============================================================================|
    ALLOCATE(TEMP1(NPROCS),TEMP2(NPROCS),TEMP3(NPROCS),TEMP4(NPROCS))
    ALLOCATE(TEMP5(NPROCS),TEMP6(NPROCS),TEMP7(NPROCS))
    TEMP1(1)  = IOBCN
    TEMP2(1) = IBCN(1)
    TEMP3(1) = IBCN(2)
    TEMP4(1) = IBCN(3)
    TEMP5(1) = IBCN(4)
    TEMP6(1) = IBCN(5)
    !   TEMP7(1) = NOBCGEO

# if defined (MULTIPROCESSOR)
    CALL MPI_ALLGATHER(IOBCN,  1,MPI_INTEGER,TEMP1,1,MPI_INTEGER,MPI_FVCOM_GROUP,IERR)
    CALL MPI_ALLGATHER(IBCN(1),1,MPI_INTEGER,TEMP2,1,MPI_INTEGER,MPI_FVCOM_GROUP,IERR)
    CALL MPI_ALLGATHER(IBCN(2),1,MPI_INTEGER,TEMP3,1,MPI_INTEGER,MPI_FVCOM_GROUP,IERR)
    CALL MPI_ALLGATHER(IBCN(3),1,MPI_INTEGER,TEMP4,1,MPI_INTEGER,MPI_FVCOM_GROUP,IERR)
    CALL MPI_ALLGATHER(IBCN(4),1,MPI_INTEGER,TEMP5,1,MPI_INTEGER,MPI_FVCOM_GROUP,IERR)
    CALL MPI_ALLGATHER(IBCN(5),1,MPI_INTEGER,TEMP6,1,MPI_INTEGER,MPI_FVCOM_GROUP,IERR)
    !   CALL MPI_GATHER(NOBCGEO,1,MPI_INTEGER,TEMP7,1,MPI_INTEGER,0,MPI_FVCOM_GROUP,IERR)
# endif


    IF(DBG_SET(DBG_LOG)) THEN
       WRITE(IPT,*)  '! BCOND TYPE             : TOTAL BC NODES; BC NODES IN EACH DOMAIN'
       WRITE(IPT,100)'!  IOBCN                 :',IOBCN_GL,   (TEMP1(I),I=1,NPROCS)
       WRITE(IPT,100)'!  IBCN(1)               :',IBCN_GL(1), (TEMP2(I),I=1,NPROCS)
       WRITE(IPT,100)'!  IBCN(2)               :',IBCN_GL(2), (TEMP3(I),I=1,NPROCS)
       WRITE(IPT,100)'!  IBCN(3)               :',IBCN_GL(3), (TEMP4(I),I=1,NPROCS)
       WRITE(IPT,100)'!  IBCN(4)               :',IBCN_GL(4), (TEMP5(I),I=1,NPROCS)
       WRITE(IPT,100)'!  IBCN(5)               :',IBCN_GL(5), (TEMP6(I),I=1,NPROCS)
       !   IF(MSR)WRITE(IPT,100)'!  NOBCGEO               :'
       !,NOBCGEO_GL, (TEMP7(I),I=1,NPROCS)
    END IF

    DEALLOCATE(TEMP1,TEMP2,TEMP3,TEMP4,TEMP5,TEMP6,TEMP7)

    RETURN
100 FORMAT(1X,A26,I6," =>",2X,4(I5,1H,))

  END SUBROUTINE SETUP_OBCTYPES
  !==============================================================================|

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|

  !==============================================================================|
  SUBROUTINE SEPARATE_OBC

    !------------------------------------------------------------------------------|
    ! Accumulate separately the amounts of nodes for 11 types of open boundaries   |
    !------------------------------------------------------------------------------|
    USE ALL_VARS
    IMPLICIT NONE

    INTEGER :: I,I1,I2,I3,I4,I5,II,J

    IBCN = 0
    IBCN_GL = 0

    DO I = 1, IOBCN_GL
       IF(TYPE_OBC_GL(I) == 1 .OR. TYPE_OBC_GL(I) == 2)  IBCN_GL(1) = IBCN_GL(1) + 1
       IF(TYPE_OBC_GL(I) == 3 .OR. TYPE_OBC_GL(I) == 4)  IBCN_GL(2) = IBCN_GL(2) + 1
       IF(TYPE_OBC_GL(I) == 5 .OR. TYPE_OBC_GL(I) == 6)  IBCN_GL(3) = IBCN_GL(3) + 1
       IF(TYPE_OBC_GL(I) == 7 .OR. TYPE_OBC_GL(I) == 8)  IBCN_GL(4) = IBCN_GL(4) + 1
       IF(TYPE_OBC_GL(I) == 9 .OR. TYPE_OBC_GL(I) == 10) IBCN_GL(5) = IBCN_GL(5) + 1
    END DO

    DO I = 1, IOBCN
       IF(TYPE_OBC(I) == 1 .OR. TYPE_OBC(I) == 2)  IBCN(1) = IBCN(1) + 1
       IF(TYPE_OBC(I) == 3 .OR. TYPE_OBC(I) == 4)  IBCN(2) = IBCN(2) + 1
       IF(TYPE_OBC(I) == 5 .OR. TYPE_OBC(I) == 6)  IBCN(3) = IBCN(3) + 1
       IF(TYPE_OBC(I) == 7 .OR. TYPE_OBC(I) == 8)  IBCN(4) = IBCN(4) + 1
       IF(TYPE_OBC(I) == 9 .OR. TYPE_OBC(I) == 10) IBCN(5) = IBCN(5) + 1
    END DO

    I1 = 0
    I2 = 0
    I3 = 0
    I4 = 0
    I5 = 0


    ALLOCATE(OBC_LST(5,MAXVAL(IBCN))) ; OBC_LST = 0

    DO I=1,IOBCN
       IF(TYPE_OBC(I) == 1 .OR. TYPE_OBC(I) == 2)THEN
          I1 = I1 + 1
          OBC_LST(1,I1) = I
       ELSE IF(TYPE_OBC(I) == 3 .OR. TYPE_OBC(I) == 4)THEN
          I2 = I2 + 1
          OBC_LST(2,I2) = I
       ELSE IF(TYPE_OBC(I) == 5 .OR. TYPE_OBC(I) == 6)THEN
          I3 = I3 + 1
          OBC_LST(3,I3) = I
       ELSE IF(TYPE_OBC(I) == 7 .OR. TYPE_OBC(I) == 8)THEN
          I4 = I4 + 1
          OBC_LST(4,I4) = I
       ELSE IF(TYPE_OBC(I) == 9 .OR. TYPE_OBC(I) == 10)THEN
          I5 = I5 + 1
          OBC_LST(5,I5) = I
       END IF
    END DO

    RETURN
  END SUBROUTINE SEPARATE_OBC
  !==============================================================================|

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|

  !==============================================================================|
  SUBROUTINE SETUP_OBC

    !------------------------------------------------------------------------------!
    USE ALL_VARS
#  if defined (SPHERICAL)
    USE MOD_SPHERICAL
#  endif
# if defined (MULTIPROCESSOR)
    USE MOD_PAR  
# endif
    IMPLICIT NONE

    REAL(SP) :: DXC,DYC,DXN,DYN,CROSS,E1,E2,DOTMAX,DOT,DX,DY,DXN_TMP,DYN_TMP
    INTEGER  :: I,J,JJ,INODE,JNODE,I1,I2,IC,N1,N2,N3
    LOGICAL  :: DEBUG

    REAL(SP), ALLOCATABLE :: NXOBC_TMP(:),NYOBC_TMP(:)

# if defined (SPHERICAL)
    REAL(DP) x1_dp,y1_dp,x2_dp,y2_dp,side
# endif
    !------------------------------------------------------------------------------!


    IF(DBG_SET(DBG_LOG))WRITE(IPT,*) "!"
    IF(DBG_SET(DBG_LOG))WRITE(IPT,*) "!           SETTING UP OPEN BOUNDARY CONDITIONS"
    IF(DBG_SET(DBG_LOG))WRITE(IPT,*) "!"

    !--Determine Adjacent Open Boundary Points-------------------------------------!

    ALLOCATE(NADJN_OBC(IOBCN))  ; NADJN_OBC = 0
    ALLOCATE(ADJN_OBC(IOBCN,2)) ; ADJN_OBC = 0

    DO I=1,IOBCN
       INODE = I_OBC_N(I)
       DO J=1,NTSN(INODE)-1
          JNODE = NBSN(INODE,J)
          IF(ISONB(JNODE) == 2 .AND. INODE /= JNODE)THEN
             NADJN_OBC(I) = NADJN_OBC(I) + 1
             ADJN_OBC(I,NADJN_OBC(I)) = JNODE
          END IF
       END DO
    END DO

    DO I=1,IOBCN
       IF(NADJN_OBC(I) == 0)THEN
          WRITE(*,*)'NO ADJACENT NODE FOUND FOR BOUNDARY NODE',I
          WRITE(*,*)'IN PROCESSOR',MYID
          CALL PSTOP
       END IF
    END DO

    !--Determine Adjacent Cells-(Nonlinear Only)-----------------------------------!
    !--Simultaneously Determine INWARD Pointing Normal NXOBC,NYOBC                 !

    ALLOCATE(NADJC_OBC(IOBCN))  ; NADJC_OBC = 0
    ALLOCATE(ADJC_OBC(IOBCN,2)) ; ADJC_OBC = 0
    ALLOCATE(NXOBC(IOBCN)) ; NXOBC = 0
    ALLOCATE(NYOBC(IOBCN)) ; NYOBC = 0
    ALLOCATE(NXOBC_TMP(IOBCN)) ; NXOBC_TMP = 0
    ALLOCATE(NYOBC_TMP(IOBCN)) ; NYOBC_TMP = 0

    DO I=1,IOBCN
       I1 = I_OBC_N(I)

       !!Mark First Cell on Boundary Edge Adjacent to Node I
       I2 = ADJN_OBC(I,1)
       DO J = 1, NTVE(I1)
          IC = NBVE(I1,J)
          N1 = NV(IC,1) ; N2 = NV(IC,2) ; N3 = NV(IC,3)
          IF( N1-I2 == 0 .OR. N2-I2 == 0 .OR. N3-I2 == 0)THEN
# if defined (SPHERICAL)
             x1_dp=VX(I1); y1_dp=VY(I1)
             x2_dp=VX(I2); y2_dp=VY(I2)
             CALL ARCX(x1_dp,y1_dp,x2_dp,y2_dp,side)
             DXN = side; DYN = (VY(I2)-VY(I1))*TPI

             x2_dp=XC(IC); y2_dp=YC(IC)
             CALL ARCX(x1_dp,y1_dp,x2_dp,y2_dp,side)
             DXC = side; DYC = (YC(IC)-VY(I1))*TPI
# else
             DXN = VX(I2)-VX(I1) ; DYN = VY(I2)-VY(I1)
             DXC = XC(IC)-VX(I1) ; DYC = YC(IC)-VY(I1)
# endif
             CROSS = SIGN(1.0_SP,DXC*DYN - DYC*DXN)
             NXOBC_TMP(I) =  CROSS*DYN/SQRT(DXN**2 +DYN**2)
             NYOBC_TMP(I) = -CROSS*DXN/SQRT(DXN**2 +DYN**2)
             NXOBC(I) = NXOBC_TMP(I)
             NYOBC(I) = NYOBC_TMP(I)
             NADJC_OBC(I) = NADJC_OBC(I) + 1
             ADJC_OBC(I,NADJC_OBC(I)) = IC
             IF(MOD(TYPE_OBC(I),2) == 1)THEN
                !!Node is Linear, Mark Cell as Linear for Flux Update
                IUCP(IC) = 0
             END IF
          END IF
       END DO

       IF(NADJN_OBC(I) > 1)THEN
          I2 = ADJN_OBC(I,2)
          DO J = 1, NTVE(I1)
             IC = NBVE(I1,J)
             N1 = NV(IC,1) ; N2 = NV(IC,2) ; N3 = NV(IC,3)
             IF( N1-I2 == 0 .OR. N2-I2 == 0 .OR. N3-I2 == 0)THEN
# if defined (SPHERICAL)
                x1_dp=VX(I1); y1_dp=VY(I1)
                x2_dp=VX(I2); y2_dp=VY(I2)
                CALL ARCX(x1_dp,y1_dp,x2_dp,y2_dp,side)
                DXN = side; DYN = (VY(I2)-VY(I1))*TPI

                x2_dp=XC(IC); y2_dp=YC(IC)
                CALL ARCX(x1_dp,y1_dp,x2_dp,y2_dp,side)
                DXC = side; DYC = (YC(IC)-VY(I1))*TPI
# else
                DXN = VX(I2)-VX(I1) ; DYN = VY(I2)-VY(I1)
                DXC = XC(IC)-VX(I1) ; DYC = YC(IC)-VY(I1)
# endif
                CROSS = SIGN(1.0_SP,DXC*DYN - DYC*DXN)
                NXOBC_TMP(I) = NXOBC_TMP(I) + CROSS*DYN/SQRT(DXN**2 +DYN**2)
                NYOBC_TMP(I) = NYOBC_TMP(I) - CROSS*DXN/SQRT(DXN**2 +DYN**2)
                NXOBC(I) = NXOBC_TMP(I)/SQRT(NXOBC_TMP(I)**2 + NYOBC_TMP(I)**2)
                NYOBC(I) = NYOBC_TMP(I)/SQRT(NXOBC_TMP(I)**2 + NYOBC_TMP(I)**2)

                NADJC_OBC(I) = NADJC_OBC(I) + 1
                ADJC_OBC(I,NADJC_OBC(I)) = IC
                IF(MOD(TYPE_OBC(I),2) == 1)THEN
                   !!Node is Linear, Mark Cell as Linear for Flux Update
                   IUCP(IC) = 0
                END IF
             END IF
          END DO
       END IF
    END DO

    DEALLOCATE(NXOBC_TMP,NYOBC_TMP)

    !--Determine Adjacent FluxFractions--------------------------------------------!

    ALLOCATE(NFLUXF_OBC(IOBCN))  ; NFLUXF_OBC = 0
    ALLOCATE(FLUXF_OBC(IOBCN,2)) ; FLUXF_OBC  = 0
    DO I=1,IOBCN
       IF(NADJN_OBC(I) == 1)THEN
          NFLUXF_OBC(I) = 1
          FLUXF_OBC(I,1) = 1.
          FLUXF_OBC(I,2) = 0.
       ELSE
          NFLUXF_OBC(I) = 2
          N1   = I_OBC_N(I)
          N2   = ADJN_OBC(I,1)
          N3   = ADJN_OBC(I,2)
# if defined (SPHERICAL)
          x1_dp=VX(N2)
          y1_dp=VY(N2)
          x2_dp=VX(N1)
          y2_dp=VY(N1)
          CALL ARCX(x1_dp,y1_dp,x2_dp,y2_dp,side)
          E1 = SQRT( side**2 + ((VY(N1)-VY(N2))*TPI)**2)

          x1_dp=VX(N3)
          y1_dp=VY(N3)
          CALL ARCX(x1_dp,y1_dp,x2_dp,y2_dp,side)
          E2 = SQRT( side**2 + ((VY(N1)-VY(N3))*TPI)**2)
# else
          E1 = SQRT( (VX(N1)-VX(N2))**2 + (VY(N1)-VY(N2))**2)
          E2 = SQRT( (VX(N1)-VX(N3))**2 + (VY(N1)-VY(N3))**2)
# endif
          FLUXF_OBC(I,1) =  E1/(E1+E2)
          FLUXF_OBC(I,2) =  E2/(E1+E2)
       END IF
    END DO
    !--Determine 1st Layer Neighbor for Open Boundary Points-----------------------!
    !--Node Chosen is Node That is Connected to OBC Node and is Oriented           !
    !--Most Normal to the Boundary.  It is not Necessarily the Closest Node        !
    !--Determine also DLTN_OBC, the normal component of the distance between       !
    !--Next_obc and the open boundary node                                         !

    DO I=1,IOBCN
       DOTMAX =  -1.0
       INODE = I_OBC_N(I)
       DO J=1,NTSN(INODE)-1
          JNODE = NBSN(INODE,J)
          IF(ISONB(JNODE) /= 2 .AND. INODE /= JNODE)THEN
# if defined (SPHERICAL)
             x1_dp=VX(INODE)
             y1_dp=VY(INODE)
             x2_dp=VX(JNODE)
             y2_dp=VY(JNODE)
             CALL ARCX(x1_dp,y1_dp,x2_dp,y2_dp,side)
             DXN_TMP = side
             DYN_TMP = (VY(JNODE)-VY(INODE))*TPI
# else
             DXN_TMP = VX(JNODE)-VX(INODE)
             DYN_TMP = VY(JNODE)-VY(INODE)
# endif
             DXN = DXN_TMP/SQRT(DXN_TMP**2 + DYN_TMP**2)
             DYN = DYN_TMP/SQRT(DXN_TMP**2 + DYN_TMP**2)
             DOT = DXN*NXOBC(I) + DYN*NYOBC(I)
             IF(DOT > DOTMAX)THEN
                DOTMAX = DOT
                NEXT_OBC(I) = JNODE
             END IF
          END IF
       END DO
    END DO

    !--Determine 2nd Layer Neighbor for Open Boundary Points-----------------------!

    DO I=1,IOBCN
       DOTMAX =  -1.0
       INODE = NEXT_OBC(I)
       IF(INODE > M) cycle
       DO J=1,NTSN(INODE)-1
          JNODE = NBSN(INODE,J)
          IF(ISONB(JNODE) /= 2)THEN
# if defined (SPHERICAL)
             x1_dp=VX(INODE)
             y1_dp=VY(INODE)
             x2_dp=VX(JNODE)
             y2_dp=VY(JNODE)
             CALL ARCX(x1_dp,y1_dp,x2_dp,y2_dp,side)
             DXN_TMP = side
             DYN_TMP = (VY(JNODE)-VY(INODE))*TPI
# else
             DXN_TMP = VX(JNODE)-VX(INODE)
             DYN_TMP = VY(JNODE)-VY(INODE)
# endif
             DXN = DXN_TMP/(SQRT(DXN_TMP**2 + DYN_TMP**2) + 1.0e-9_SP)
             DYN = DYN_TMP/(SQRT(DXN_TMP**2 + DYN_TMP**2) + 1.0e-9_SP)
             DOT = DXN*NXOBC(I) + DYN*NYOBC(I)
             IF(DOT > DOTMAX)THEN
                DOTMAX = DOT
                NEXT_OBC2(I) = JNODE
             END IF
          END IF
       END DO
    END DO

    !--Determine DLTN_OBC----------------------------------------------------------!
    ALLOCATE(DLTN_OBC(IOBCN))
    DO I=1,IOBCN
       I1 = I_OBC_N(I)
       I2 = NEXT_OBC(I)

# if defined (SPHERICAL)
       x1_dp=VX(I1)
       y1_dp=VY(I1)
       x2_dp=VX(I2)
       y2_dp=VY(I2)
       call ARCX(x1_dp,y1_dp,x2_dp,y2_dp,side)
       DX = side
       DY = (VY(I2)-VY(I1))*TPI
# else
       DX = VX(I2)-VX(I1)
       DY = VY(I2)-VY(I1)
# endif
       DLTN_OBC(I) = ABS(DX*NXOBC(I) + DY*NYOBC(I))
    END DO

!!$   RETURN
!!$!--Dump Information to Matlab Files for Checking-------------------------------!
!!$
!!$   OPEN(UNIT=81,FILE='mesh.scatter',FORM='formatted')
!!$   DO I=1,M
!!$     WRITE(81,*)vx(i),vy(i)
!!$   END DO
!!$   CLOSE(81)
!!$   OPEN(UNIT=81,FILE='nextobc.scatter',FORM='formatted')
!!$   DO I=1,IOBCN
!!$     I1 = NEXT_OBC(I)
!!$     WRITE(81,*)VX(I1),VY(I1)
!!$   END DO
!!$   CLOSE(81)
!!$   OPEN(UNIT=81,FILE='nextobc2.scatter',FORM='formatted')
!!$   DO I=1,IOBCN
!!$     I1 = NEXT_OBC2(I)
!!$     WRITE(81,*)VX(I1),VY(I1)
!!$   END DO
!!$   CLOSE(81)
!!$   OPEN(UNIT=81,FILE='iobcn.scatter',FORM='formatted')
!!$   DO I=1,IOBCN
!!$     I1 = I_OBC_N(I)
!!$     WRITE(81,*)VX(I1),VY(I1) 
!!$   END DO
!!$   CLOSE(81)
!!$   OPEN(UNIT=81,FILE='obcnorm.scatter',FORM='formatted')
!!$   DO I=1,IOBCN
!!$     I1 = I_OBC_N(I)
!!$     WRITE(81,*)NXOBC(I1),NYOBC(I1) 
!!$   END DO
!!$   CLOSE(81)
!!$   OPEN(UNIT=81,FILE='nonlinear.scatter',FORM='formatted')
!!$   DO I=1,IOBCN
!!$     IF(NADJC_OBC(I) > 0) WRITE(81,*)XC(ADJC_OBC(I,1)),YC(ADJC_OBC(I,1)) 
!!$     IF(NADJC_OBC(I) > 1) WRITE(81,*)XC(ADJC_OBC(I,2)),YC(ADJC_OBC(I,2)) 
!!$   END DO
!!$   CLOSE(81)
!!$   OPEN(UNIT=81,FILE='linear.scatter',FORM='formatted')
!!$   DO I=1,N
!!$     IF(IUCP(I)==0)THEN
!!$       WRITE(81,*)XC(I),YC(I) 
!!$     END IF
!!$   END DO
!!$   CLOSE(81)

    RETURN
  END SUBROUTINE SETUP_OBC
  !==============================================================================|

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|

  !==============================================================================|
# if !defined (SEMI_IMPLICIT)
  SUBROUTINE FLUX_OBN(K)

    USE ALL_VARS
    IMPLICIT NONE

    INTEGER, INTENT(IN)  :: K
    INTEGER              :: I,J,C1,C2
    REAL(SP)             :: FF,FLUX 

    FLUXOBN = 0.0_SP

    DO I = 1, IOBCN

       J = I_OBC_N(I)
       !Compute Boundary Flux From Continuity Flux Defect
       FLUX = -(ELF(J)-ELRK(J))*ART1(J)/(ALPHA_RK(K)*DTE)-XFLUX_OBCN(I)

       !Set Flux In Adjacent Nonlinear BC Element 1 (If Exists)
       IF(NADJC_OBC(I) > 0) THEN
          C1 = ADJC_OBC(I,1) 
          FF = FLUXF_OBC(I,1)
          FLUXOBN(C1) = FLUXOBN(C1) + FF*FLUX 
       END IF

       !Set Flux In Adjacent Nonlinear BC Element 2 (If Exists)
       IF(NADJC_OBC(I) > 1) THEN
          C2 =  ADJC_OBC(I,2) 
          FF = FLUXF_OBC(I,2)
          FLUXOBN(C2) = FLUXOBN(C2) + FF*FLUX 
       END IF

    END DO

    RETURN
  END SUBROUTINE FLUX_OBN
# endif
  !========================================================================

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|

  !==============================================================================|
  SUBROUTINE BCOND_T_PERTURBATION(T2D_NEXT,T2D,TTMP,I,J,J1)
    !==============================================================================|
    ! Calculate the OBC for temperature perturbation                               |
    !==============================================================================|
    USE ALL_VARS
    IMPLICIT NONE

    !   INTEGER :: I1,I2,J,JN
    INTEGER :: I,J,J1,J2,K
    REAL(SP):: CC,CP,MU,CL
    REAL(SP):: PERT_NEXT,PERT,T2D_NEXT,T2D
    REAL(SP):: T2D_NEXT1,TM12D_NEXT2,TM12D_NEXT1,TM22D_NEXT1
    REAL(SP):: TTMP(IOBCN,KBM1)

    SELECT CASE(TYPE_TSOBC)

    CASE(1)
       DO K=1,KBM1
          TTMP(I,K) = TF1(J1,K) - T2D_NEXT
       END DO
    CASE(2)
       CC = SQRT(GRAV_N(J)*H(J))*DTI/DLTN_OBC(I)
       CP = CC + 1.0_SP
       DO K=1,KBM1
          PERT_NEXT = TF1(J1,K) - T2D_NEXT
          PERT      = T1(J,K) - T2D
          TTMP(I,K) = (CC*PERT_NEXT + PERT)/CP
       END DO
    CASE(3)
       CC = SQRT(GRAV_N(J)*H(J))*DTI/DLTN_OBC(I)
       CP = CC + 1.0_SP
       DO K=1,KBM1
          PERT_NEXT = TF1(J1,K) - T2D_NEXT
          PERT      = T1(J,K) - T2D
          TTMP(I,K) = (CC*PERT_NEXT + PERT*(1.0_SP - DTI/10800.0_SP))/CP
       END DO
    CASE(4)
       J2 = NEXT_OBC2(I)
       T2D_NEXT1  =0.0_SP
       TM12D_NEXT2=0.0_SP
       TM12D_NEXT1=0.0_SP
       TM22D_NEXT1=0.0_SP
       DO K=1,KBM1
          T2D_NEXT1  =T2D_NEXT1  +T1(J1,K)*DZ(J1,K)
          TM12D_NEXT2=TM12D_NEXT2+T1M1(J2,K)*DZ(J2,K)
          TM12D_NEXT1=TM12D_NEXT1+T1M1(J,K)*DZ(J,K)
          TM22D_NEXT1=TM22D_NEXT1+T1M2(J1,K)*DZ(J1,K)
       END DO

       DO K=1,KBM1
          CL = ((T1M2(J1,K)-TM22D_NEXT1)-(T1(J1,K)-T2D_NEXT1))/   &
               ((T1(J1,K)-T2D_NEXT1)+(T1M2(J1,K)-TM22D_NEXT1)     &
               -2.0*(T1M1(J2,K)-TM12D_NEXT2))
          IF(CL >= 1.0)THEN
             MU = 1.0
          ELSE IF(CL > 0.0 .AND. CL < 1.0)THEN
             MU = CL
          ELSE
             MU = 0.0
          END IF

          TTMP(I,K)=((T1M1(J,K)-TM12D_NEXT1)*(1.0-MU)     &
               +2.0*MU*(T1(J1,K)-T2D_NEXT1))/(1.0+MU)
       END DO

    CASE DEFAULT
       CALL FATAL_ERROR("INVALID OBC_TS_TYPE IN NML_OPEN_BOUNDARY_CONTROL"&
            &, "See mod_obcs.F")

    END SELECT

    RETURN
  END SUBROUTINE BCOND_T_PERTURBATION
  !========================================================================

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|

  !==============================================================================|
  SUBROUTINE BCOND_S_PERTURBATION(S2D_NEXT,S2D,STMP,I,J,J1)
    !==============================================================================|
    ! Calculate the OBC for salinity perturbation                                  |
    !==============================================================================|
    USE ALL_VARS
    IMPLICIT NONE

    INTEGER :: I,J,J1,J2,K
    REAL(SP):: CC,CP,MU,CL
    REAL(SP):: PERT_NEXT,PERT,S2D_NEXT,S2D
    REAL(SP):: S2D_NEXT1,SM12D_NEXT2,SM12D_NEXT1,SM22D_NEXT1
    REAL(SP):: STMP(IOBCN,KBM1)

    SELECT CASE(TYPE_TSOBC)

    CASE(1)
       DO K=1,KBM1
          STMP(I,K) = SF1(J1,K) - S2D_NEXT
       END DO
    CASE(2)
       CC = SQRT(GRAV_N(J)*H(J))*DTI/DLTN_OBC(I)
       CP = CC + 1.0_SP
       DO K=1,KBM1
          PERT_NEXT = SF1(J1,K) - S2D_NEXT
          PERT      = S1(J,K) - S2D
          STMP(I,K) = (CC*PERT_NEXT + PERT)/CP
       END DO
    CASE(3)
       CC = SQRT(GRAV_N(J)*H(J))*DTI/DLTN_OBC(I)
       CP = CC + 1.0_SP
       DO K=1,KBM1
          PERT_NEXT = SF1(J1,K) - S2D_NEXT
          PERT      = S1(J,K) - S2D
          STMP(I,K) = (CC*PERT_NEXT + PERT*(1.0_SP - DTI/10800.0_SP))/CP
       END DO
    CASE(4)
       J2 = NEXT_OBC2(I)
       S2D_NEXT1  =0.0_SP
       SM12D_NEXT2=0.0_SP
       SM12D_NEXT1=0.0_SP
       SM22D_NEXT1=0.0_SP
       DO K=1,KBM1
          S2D_NEXT1  =S2D_NEXT1  +S1(J1,K)*DZ(J1,K)
          SM12D_NEXT2=SM12D_NEXT2+S1M1(J2,K)*DZ(J2,K)
          SM12D_NEXT1=SM12D_NEXT1+S1M1(J,K)*DZ(J,K)
          SM22D_NEXT1=SM22D_NEXT1+S1M2(J1,K)*DZ(J1,K)
       END DO

       DO K=1,KBM1
          CL = ((S1M2(J1,K)-SM22D_NEXT1)-(S1(J1,K)-S2D_NEXT1))/   &
               ((S1(J1,K)-S2D_NEXT1)+(S1M2(J1,K)-SM22D_NEXT1)     &
               -2.0*(S1M1(J2,K)-SM12D_NEXT2))
          IF(CL >= 1.0)THEN
             MU = 1.0
          ELSE IF(CL > 0.0 .AND. CL < 1.0)THEN
             MU = CL
          ELSE
             MU = 0.0
          END IF

          STMP(I,K)=((S1M1(J,K)-SM12D_NEXT1)*(1.0-MU)     &
               +2.0*MU*(S1(J1,K)-S2D_NEXT1))/(1.0+MU)
       END DO

    CASE DEFAULT
       CALL FATAL_ERROR("INVALID OBC_TS_TYPE IN NML_OPEN_BOUNDARY_CONTROL"&
            &, "See mod_obcs.F")

    END SELECT

    RETURN
  END SUBROUTINE BCOND_S_PERTURBATION
  !========================================================================

  SUBROUTINE GDAY1(IDD,IMM,IYY,ICC,KD)
!
!  given day,month,year and century(each 2 digits), gday returns
!  the day#, kd based on the gregorian calendar.
!  the gregorian calendar, currently 'universally' in use was
!  initiated in europe in the sixteenth century. note that gday
!  is valid only for gregorian calendar dates.
!
!  kd=1 corresponds to january 1, 0000
!
!  note that the gregorian reform of the julian calendar 
!  omitted 10 days in 1582 in order to restore the date
!  of the vernal equinox to march 21 (the day after
!  oct 4, 1582 became oct 15, 1582), and revised the leap 
!  year rule so that centurial years not divisible by 400
!  were not leap years.
!
!  this routine was written by eugene neufeld, at ios, in june 1990.
!
    integer idd, imm, iyy, icc, kd
    integer ndp(13)
    integer ndm(12)
    data ndp/0,31,59,90,120,151,181,212,243,273,304,334,365/
    data ndm/31,28,31,30,31,30,31,31,30,31,30,31/
!
!  test for invalid input:
    if(icc.lt.0)then
!     write(11,5000)icc
      call pstop
    endif
    if(iyy.lt.0.or.iyy.gt.99)then
!     write(11,5010)iyy
      call pstop
    endif
    if(imm.le.0.or.imm.gt.12)then
!     write(11,5020)imm
      call pstop
      endif
    if(idd.le.0)then
!     write(11,5030)idd
      call pstop
    endif
    if(imm.ne.2.and.idd.gt.ndm(imm))then
!     write(11,5030)idd
      call pstop
    endif
    if(imm.eq.2.and.idd.gt.29)then
!     write(11,5030)idd
      call pstop
    endif
    if(imm.eq.2.and.idd.gt.28.and.((iyy/4)*4-iyy.ne.0.or.(iyy.eq.0.and.(icc/4)*4-icc.ne.0)))then
!     write(11,5030)idd
      call pstop
    endif
5000  format(' input error. icc = ',i7)
5010  format(' input error. iyy = ',i7)
5020  format(' input error. imm = ',i7)
5030  format(' input error. idd = ',i7)
!
!  calculate day# of last day of last century:
    kd = icc*36524 + (icc+3)/4
!
!  calculate day# of last day of last year:
    kd = kd + iyy*365 + (iyy+3)/4
!
!  adjust for century rule:
!  (viz. no leap-years on centurys except when the 2-digit
!  century is divisible by 4.)
    if(iyy.gt.0.and.(icc-(icc/4)*4).ne.0) kd=kd-1
!  kd now truly represents the day# of the last day of last year.
!
!  calculate day# of last day of last month:
    kd = kd + ndp(imm)
!
!  adjust for leap years:
    if(imm.gt.2.and.((iyy/4)*4-iyy).eq.0.and.((iyy.ne.0).or.(((icc/4)*4-icc).eq.0))) kd=kd+1
!  kd now truly represents the day# of the last day of the last
!  month.
!
!  calculate the current day#:
    kd = kd + idd

  RETURN
  END SUBROUTINE GDAY1 

  subroutine astro(d1,h,pp,s,p,np,dh,dpp,ds,dp,dnp)
! this subroutine calculates the following five ephermides
! of the sun and moon
! h = mean longitude of the sum
! pp = mean longitude of the solar perigee
! s = mean longitude of the moon
! p = mean longitude of the lunar perigee
! np = negative of the longitude of the mean ascending node
! and their rates of change.
! units for the ephermides are cycles and for their derivatives
! are cycles/365 days
! the formulae for calculating this ephermides were taken from
! pages 98 and 107 of the explanatory supplement to the
! astronomical ephermeris and the american ephermis and
! nautical almanac (1961)
!
  implicit real*8(a-h,o-z)
  real*8 np
  d2=d1*1.d-4
  f=360.d0
  f2=f/365.d0
  h=279.696678d0+.9856473354d0*d1+.00002267d0*d2*d2
  pp=281.220833d0+.0000470684d0*d1+.0000339d0*d2*d2+.00000007d0*d2**3
  s=270.434164d0+13.1763965268d0*d1-.000085d0*d2*d2+.000000039d0*d2**3
  p=334.329556d0+.1114040803d0*d1-.0007739d0*d2*d2-.00000026d0*d2**3
  np=-259.183275d0+.0529539222d0*d1-.0001557d0*d2*d2-.00000005d0*d2**3
  h=h/f
  pp=pp/f
  s=s/f
  p=p/f
  np=np/f
  h=h-dint(h)
  pp=pp-dint(pp)
  s=s-dint(s)
  p=p-dint(p)
  np=np-dint(np)
  dh=.9856473354d0+2.d-8*.00002267d0*d1
  dpp=.0000470684d0+2.d-8*.0000339d0*d1+3.d-12*.00000007d0*d1**2
  ds=13.1763965268d0-2.d-8*.000085d0*d1+3.d-12*.000000039d0*d1**2
  dp=.1114040803d0-2.d-8*.0007739d0*d1-3.d-12*.00000026d0*d1**2
  dnp=+.0529539222d0-2.d-8*.0001557d0*d1-3.d-12*.00000005d0*d1**2
  dh=dh/f2
  dpp=dpp/f2
  ds=ds/f2
  dp=dp/f2
  dnp=dnp/f2
  return
  end subroutine astro

!!$

END MODULE MOD_OBCS
