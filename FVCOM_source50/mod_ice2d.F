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

MODULE MOD_ICE2D
#  if defined (ICE)
  USE ALL_VARS
  USE MOD_PREC
#  if defined (SPHERICAL)
  USE MOD_SPHERICAL
  !#  if defined (NORTHPOLE)
  USE MOD_NORTHPOLE
  !#  endif
#  endif
  
# if defined (MULTIPROCESSOR)
  USE MOD_PAR
# endif
  
  !===============================================================================|
  !   use CICE      
  !===============================================================================|
  USE ICE_KINDS_MOD
  USE ICE_STATE
  USE ICE_MECHRED
  USE ICE_CONSTANTS
  USE ICE_CALENDAR, only: dtice,dyn_dt,ndyn_dt,dtei
  
  IMPLICIT NONE
  SAVE
  
  Logical :: EVP_ON   
  REAL(SP), PARAMETER   :: ECC = 4.0_SP
  REAL(SP), PARAMETER   :: ECCI = C1I/ECC
  real (kind=dbl_kind) :: &
       Tdte             &! subcycling timestep for EVP dynamics, s
       ,  tdamp2           ! 2(wave damping time scale T)
  
  REAL(kind=dbl_kind)   ::  &
       dte2T      &  ! dte/2T
       ,denom1     &  ! constants for stress equation
       ,  denom2        !
  
  real(SP), ALLOCATABLE :: strocnxE(:),strocnyE(:)! ice-ocean stress
  real(SP), ALLOCATABLE :: STRAIRXE(:),STRAIRYE(:)
  
  REAL(SP), ALLOCATABLE :: AICETMP(:),AIU(:)
  REAL(SP), ALLOCATABLE :: UICE2(:),VICE2(:),UICE2F(:),VICE2F(:)

  REAL(SP), ALLOCATABLE :: UICE2N(:), VICE2N(:)

  REAL(SP), ALLOCATABLE :: UICE2RK(:),VICE2RK(:)
  REAL(SP), ALLOCATABLE :: PSTXICE(:),PSTYICE(:)
  REAL(SP), ALLOCATABLE :: ADVSTRX(:),ADVSTRY(:)
  REAL(SP), ALLOCATABLE :: ADVUICE(:),ADVVICE(:)
  REAL(SP), ALLOCATABLE :: IMASS1(:),IMASS(:),PICE(:)
  REAL(SP), ALLOCATABLE :: TAUXWI(:),TAUYWI(:),VREL(:)
  INTEGER,  ALLOCATABLE :: ISICEN(:),ISICEC(:),ISICEN_OLD(:),ISICEC_OLD(:)

  INTEGER,  ALLOCATABLE :: NN_STRESS(:),NN_STRESS_OLD(:)




#  if defined (ICE_EMBEDDING)  
!----J. Ge for ice baroclinic pressure gradient-----!
  REAL(SP) :: NET_DRHOX ! net baroclinic gradient force for whole ice
  REAL(SP) :: NET_DRHOY ! net baroclinic gradient force for whole ice
  integer,  allocatable :: iwbnd(:) ! cell boundary between ice and water 
  integer,  allocatable :: Kzice(:)  ! bottom sigma index of icepack
  REAL(SP), ALLOCATABLE :: DRHOX_ICE (:,:) !baroclinic gradient force
  REAL(SP), ALLOCATABLE :: DRHOY_ICE (:,:) !baroclinic gradient force
!----J. Ge for ice baroclinic pressure gradient-----!
#  endif

  real (kind=dbl_kind), dimension (:,:) ,allocatable,save ::&     
!     &   shear       ! strain rate II component (1/s)
       divu       &   ! strain rate I component, velocity divergence (1/s)
       ,Delta         ! function of strain rates (1/s)
  
  real (kind=SP), dimension (:) ,allocatable,save ::&
       DivuN        & ! strain rate I component, velocity divergence (1/s)
       ,DeltaN      & ! function of strain rates (1/s)
       ,DivuE       & ! strain rate I component, velocity divergence (1/s)
       ,DeltaE        ! function of strain rates (1/s)
  
  integer (kind=int_kind) ::  &
        kdyn        & ! type of dynamics ( 1 = evp )
      , ndte          ! number of subcycles:  ndte=dyn_dt/dte
  
  REAL(DP), PARAMETER :: eyc =0.36
  REAL(DP) :: rcon
  
  CHARACTER(LEN=80) :: RHEOLOGY

  REAL(SP), PARAMETER :: DRAGW = 0.0055_SP*RHOW
  ! drag coefficient for water on ice *rhow (kg/m^3)
# if defined (ICE_FRESHWATER)
  ! afm 20151112 & EJA 20160921 - turning angle 0
  REAL(SP), PARAMETER :: COSW = C1I  !cos(ocean turning angle) !turning angle = 0
  REAL(SP), PARAMETER :: SINW = C0I  !sin(ocean turning angle) !turning angle = 0
# else
  REAL(SP), PARAMETER :: COSW = 0.9063_SP  !cos(ocean turning angle) !turning angle = 25.
  REAL(SP), PARAMETER :: SINW = 0.4226_SP  !sin(ocean turning angle) !turning angle = 25.
# endif

# if defined (ICE_FRESHWATER)
  ! afm 20151112 & EJA 20160921 - turning angle 0
  REAL(SP), PARAMETER :: COSA = c1i  !cos(wind turning angle) !turning angle = 0.
  REAL(SP), PARAMETER :: SINA = c0i  !sin(wind turning angle) !turning angle = 0.
# else
  REAL(SP), PARAMETER :: COSA = 0.9063_SP  !cos(wind turning angle) !turning angle = 25.
  REAL(SP), PARAMETER :: SINA = 0.4226_SP  !sin(wind turning angle) !turning angle = 25.
# endif

  ! afm 20151112 & EJA 20160921 - not used for freshwater
# if !defined (ICE_FRESHWATER) 
  REAL(SP), PARAMETER :: COSA30 = 0.86602_SP  !cos(wind turning angle) !turning angle = 30.
  REAL(SP), PARAMETER :: SINA30 = 0.50000_SP  !sin(wind turning angle) !turning angle = 30.
 
  REAL(SP), PARAMETER :: COSA20 = 0.93969_SP  !cos(wind turning angle) !turning angle = 20.
  REAL(SP), PARAMETER :: SINA20 = 0.34202_SP  !sin(wind turning angle) !turning angle = 20.
# endif

  !!---------------------------------------------------------------------------------   
  real (kind=dbl_kind), dimension (:) ,allocatable,save ::&
       prs_sig       ! replacement pressure, for stress calc
  REAL(SP), ALLOCATABLE :: SIG1(:), SIG2(:) 
  REAL(SP), ALLOCATABLE :: SIG11(:),SIG22(:),SIG12(:)
  REAL(SP), ALLOCATABLE :: PSIG1(:), PSIG2(:) !!principal_stress
  
  CONTAINS

!===============================================================================|
!===============================================================================|
   SUBROUTINE ALLOC_UVICE
   
   IMPLICIT NONE
  
   CHARACTER(LEN=120) :: FNAME
   INTEGER  :: ISCAN
   
   ALLOCATE(UICE2(0:NT))     ;    UICE2 = 0.0_SP
   ALLOCATE(VICE2(0:NT))     ;    VICE2 = 0.0_SP
   ALLOCATE(AIU(0:NT))       ;      AIU = 0.0_SP
 
   ALLOCATE(strocnxE(0:NT))  ; strocnxE =0.0_SP
   ALLOCATE(strocnyE(0:NT))  ; strocnyE =0.0_SP

   ALLOCATE(ISICEC(0:NT))    ; ISICEC     = 0
   ALLOCATE(ISICEC_OLD(0:NT)); ISICEC_OLD = 0
   ALLOCATE(ISICEN(0:MT))    ; ISICEN     = 0
   ALLOCATE(ISICEN_OLD(0:MT)); ISICEN_OLD = 0

   ALLOCATE(NN_STRESS(0:MT));     NN_STRESS     = 0
   ALLOCATE(NN_STRESS_OLD(0:MT)); NN_STRESS_OLD = 0

#  if defined (ICE_EMBEDDING)  
!----J. Ge for ice baroclinic pressure gradient-----!
   allocate(iwbnd(0:nt))     ;  iwbnd = 0
   allocate(kzice(0:nt))     ;  kzice  = 0
   allocate(DRHOX_ICE(0:nt,kb)) ; DRHOX_ICE = 0.0_sp
   allocate(DRHOY_ICE(0:nt,kb)) ; DRHOY_ICE = 0.0_sp
!----J. Ge for ice baroclinic pressure gradient-----!
#  endif

   !  for ridge redistribution
   allocate(Divu(M,1))       ;  Divu  =0.0
   allocate(Delta(M,1))      ;  Delta =0.0

   allocate(DivuN(0:MT))     ; DivuN  =0.0
   allocate(DeltaN(0:MT))    ; DeltaN =0.0

   allocate(DivuE(0:NT))     ; DivuE  =0.0
   allocate(DeltaE(0:NT))    ; DeltaE =0.0

   ALLOCATE(SIG1(0:MT))      ;  SIG1 =0.0_SP
   ALLOCATE(SIG2(0:MT))      ;  SIG2 =0.0_SP
   ALLOCATE(SIG11(0:MT))     ;  SIG11=0.0_SP
   ALLOCATE(SIG22(0:MT))     ;  SIG22=0.0_SP
   ALLOCATE(SIG12(0:MT))     ;  SIG12=0.0_SP

   ALLOCATE(PRS_SIG(0:MT))   ;  PRS_SIG=0.0_SP
   ALLOCATE(PSIG1(0:MT))     ;  PSIG1  =0.0_SP
   ALLOCATE(PSIG2(0:MT))     ;  PSIG2  =0.0_SP

      dyn_dt=dtice/real(ndyn_dt,kind=dbl_kind)

# if defined (ICE_FRESHWATER)
! afm 20151113 & EJA 20160921
! ndte=120-->30, based on time stepping analysis, it should be allowable
! and saves some computational time
     ndte = 30
# else
     ndte =120 !!  
# endif
 
      Tdte = dyn_dt/real(ndte,kind=dbl_kind)        ! s
      dtei = c1i/Tdte              ! 1/s
      tdamp2 = c2I*eyc*dyn_dt         ! s

      ! constants for stress equation
      dte2T = Tdte/tdamp2                    ! unitless
      denom1 = c1i/(c1i+dte2T)
      denom2 = c1i/(c1i+dte2T*ecc)

      rcon  = 1230._dbl_kind*eyc*dyn_dt*dtei**2  ! kg/s

      IF(MSR) write(ipt,*)' ICE dynamice time step  dyn_dt =' , dyn_dt
   
   RETURN
   END SUBROUTINE ALLOC_UVICE
!===============================================================================|
!
!===============================================================================|
   SUBROUTINE ICE_RUNGE_KUTTA

   IMPLICIT NONE
!   integer (kind=int_kind), intent(in) :: kstrngth
   integer (kind=int_kind) :: kstrngth
   INTEGER :: K   ,KID
   integer :: I,J,I1
   real, allocatable :: rsub(:,:)
   real(SP) :: aice_tmp
   REAL(SP)  ::  UI_TMP,VI_TMP
   INTEGER   ::  N_TMP
   INTEGER   :: IUVN
   REAL(SP)  :: ART_TMP
!===============================================================================|
   !  temp use for calculate the air -ice stress
   real :: vmag_tmp      &! surface wind magnitude   (m/s)
      ,      rd_tmp      &! sqrt of exchange coefficient (momentum)
      ,     rdn_tmp      &! sqrt of exchange coefficient (momentum)
      ,     tau_tmp      &! stress at zlvl
      ,   ustar_tmp      &! ustar (m/s)
      ,    rhoa_tmp      &! air density (kg/m^3)
      ,windspeed_tmp

!===============================================================================|
   ! major/minor axis length ratio, squared
      
   ALLOCATE(AICETMP(0:MT));   AICETMP = 0.0_SP
   ALLOCATE(UICE2F(0:NT)) ;   UICE2F  = 0.0_SP
   ALLOCATE(VICE2F(0:NT)) ;   VICE2F  = 0.0_SP
   ALLOCATE(UICE2RK(0:NT));   UICE2RK = 0.0_SP
   ALLOCATE(VICE2RK(0:NT));   VICE2RK = 0.0_SP
   ALLOCATE(PSTXICE(0:NT));   PSTXICE = 0.0_SP
   ALLOCATE(PSTYICE(0:NT));   PSTYICE = 0.0_SP
   ALLOCATE(ADVSTRX(0:NT));   ADVSTRX = 0.0_SP
   ALLOCATE(ADVSTRY(0:NT));   ADVSTRY = 0.0_SP
   ALLOCATE(ADVUICE(0:NT));   ADVUICE = 0.0_SP
   ALLOCATE(ADVVICE(0:NT));   ADVVICE = 0.0_SP
   ALLOCATE(IMASS1(0:NT)) ;   IMASS1  = 0.0_SP
   ALLOCATE(PICE(0:MT))   ;   PICE    = 0.0_SP
   ALLOCATE(IMASS(0:MT))  ;   IMASS   = 0.0_SP
   ALLOCATE(TAUXWI(0:NT)) ;   TAUXWI  = 0.0_SP
   ALLOCATE(TAUYWI(0:NT)) ;   TAUYWI  = 0.0_SP
   ALLOCATE(VREL(0:NT))   ;   VREL  = 0.0_SP

   ALLOCATE(UICE2N(0:MT)); UICE2N =0.0_SP
   ALLOCATE(VICE2N(0:MT)); VICE2N =0.0_SP

   !-----------------------------------------------------------------
   ! pressure and forcing terms; set sigma=0 for no ice;
   ! initialize ice velocity in cells previously empty to ocn current
   !-----------------------------------------------------------------

   CALL ICE_STRENGTH(KSTRENGTH)
   DO I=1,M !T
      PICE(I) = STRENGTH(I,1)
   END DO
   PICE =PICE*RAMP
# if defined (MULTIPROCESSOR)
   IF(PAR) CALL AEXCHANGE(NC,MYID,NPROCS,PICE)
# endif

   !-----------------------------------------------------------------
   ! total mass of ice and snow, centered in T-cell
   ! NOTE: vice and vsno must be up to date in all grid cells,
   !       including ghost cells
   !-----------------------------------------------------------------

    ISICEC_OLD = ISICEC  !  ice in last step
    ISICEC = 0
    ISICEN = 0
    IMASS  =0.0_SP
    IMASS1 =0.0_SP

   DO I=1,M 
      IF (AICE(I,1) >p001) then
      IMASS(I) = (RHOI*VICE(I,1) + RHOS*VSNO(I,1)) ! kg/m^2
      IF(IMASS(I) >P01)  ISICEN(I) =1
      END IF
   END DO
# if defined (MULTIPROCESSOR)
      if(PAR) CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,1,MYID,NPROCS,IMASS)
      IF(PAR) CALL AEXCHANGE(NC,MYID,NPROCS,IMASS)
# endif

   CALL  N2E2D(IMASS,IMASS1)

   DO I=1,M !T 
     AICETMP(I) = AICE(I,1)
   END DO
   AIU =0.0_SP

# if defined (MULTIPROCESSOR)
       IF(PAR) CALL AEXCHANGE(NC,MYID,NPROCS,AICETMP)
# endif

   CALL N2E2D(AICETMP,AIU)

     !EXCHANGE ELEMENT-BASED VALUES ACROSS THE INTERFACE
#    if defined (MULTIPROCESSOR)
     IF(PAR)CALL AEXCHANGE(EC,MYID,NPROCS,IMASS1,AIU)
#    endif
   
   DO I=1,N !T
      IF(sum(ISICEN(NV(I,1:3))) > 1 .and.IMASS1(I)>p01  &
         .and. AIU(I)>P001) &
      ISICEC(i)=1
   END DO

#    if defined (MULTIPROCESSOR)
     IF(PAR)CALL AEXCHANGE(EC,MYID,NPROCS,ISICEC)
#    endif

   DO I=1,N 
     UICE2(I) =UICE2(I)*ISICEC(I)
     VICE2(I) =VICE2(I)*ISICEC(I)
   END DO
#    if defined (MULTIPROCESSOR)
     IF(PAR)CALL AEXCHANGE(EC,MYID,NPROCS,UICE2,VICE2)
#    endif

   !-----------------------------------------------------------   
   !  calculate the wind stress on ice
   !-----------------------------------------------------------
     allocate(STRAIRXE(0:NT),STRAIRYE(0:NT))
     STRAIRXE = 0.0_SP
     STRAIRYE = 0.0_SP

      !------------------------------------------------------------
      ! Define functions
      !------------------------------------------------------------
      !------------------------------------------------------------
      ! Compute turbulent flux coefficients, wind stress, and 
      !------------------------------------------------------------
      !   rdn_tmp  = vonkar/log(zref/iceruf) ! neutral coefficient
!        rdn_tmp = 0.4/log(10./0.0005)
         DO I=1,N
          windspeed_tmp=SQRT(uuwind(i)*uuwind(i)+vvwind(i)*vvwind(i))
!           if ( aicen(i,j,ni) > puny) then
      !------------------------------------------------------------
      ! define some needed variables
      !------------------------------------------------------------
         vmag_tmp = max(1.0, windspeed_tmp)

# if defined (ICE_FRESHWATER)
! afm turning angle 0 ------------- EJA 20160921
         STRAIRXE(I) =1.3E-3*(1.1+0.04*vmag_tmp)*vmag_tmp*(UUWIND(I)*COSA - VVWIND(I)*SINA)*AIU(I)
         STRAIRYE(I) =1.3E-3*(1.1+0.04*vmag_tmp)*vmag_tmp*(UUWIND(I)*SINA + VVWIND(I)*COSA)*AIU(I)
# else
         IF(vmag_tmp<15.0_SP) THEN
         STRAIRXE(I) =1.3E-3*(1.1+0.04*vmag_tmp)*vmag_tmp*(UUWIND(I)*COSA30 - VVWIND(I)*SINA30)*AIU(I) 
         STRAIRYE(I) =1.3E-3*(1.1+0.04*vmag_tmp)*vmag_tmp*(UUWIND(I)*SINA30 + VVWIND(I)*COSA30)*AIU(I) 
         ELSE
         STRAIRXE(I) =1.3E-3*(1.1+0.04*vmag_tmp)*vmag_tmp*(UUWIND(I)*COSA20 - VVWIND(I)*SINA20)*AIU(I)
         STRAIRYE(I) =1.3E-3*(1.1+0.04*vmag_tmp)*vmag_tmp*(UUWIND(I)*SINA20 + VVWIND(I)*COSA20)*AIU(I)
         ENDIF
# endif

      END DO

        STRAIRXE=STRAIRXE*RAMP
        STRAIRYE=STRAIRYE*RAMP     


if (.true.) then !dwang
!------------J. Ge for baroclinic pressure gradient force----------------------!
#  if defined (ICE_EMBEDDING)
  
   call find_ice_water_boundary
   call get_kzice
   call baropg_ice

#  endif 
end if !dwang
!------------------------------------------------------------------------------!



!-----------------------------------------------------------
!
!------BEGIN MAIN LOOP OVER EXTERNAL MODEL 4 STAGE RUNGE-KUTTA INTEGRATION-----!
!
      EVP_ON= .TRUE.
      
      IF(EVP_ON) THEN

         !CALL ICEOCEAN_STRESS
         DO KID=1,ndte              ! subcycling
            
            UICE2N=0.0_SP;VICE2N=0.0_SP
            DO I=1,M
               IF(ISICEN(I)==1)THEN
                  IUVN=0
                  ART_TMP=0.0_SP
                  DO I1=1, NTVE(I)
                     UICE2N(I) =UICE2N(I)+UICE2(NBVE(I,I1))*ISICEC(NBVE(I,I1))
                     VICE2N(I) =VICE2N(I)+VICE2(NBVE(I,I1))*ISICEC(NBVE(I,I1))
                     IUVN   =IUVN+ISICEC(NBVE(I,I1))
                  END DO
                  IF(IUVN>0) THEN
                     UICE2N(I) = UICE2N(I) /FLOAT(IUVN)
                     VICE2N(I) = VICE2N(I) /FLOAT(IUVN)
                  ENDIF
               ENDIF
            END DO
# if defined(MULTIPROCESSOR)
            IF(PAR)CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,1,MYID,NPROCS,UICE2N,VICE2N)
            IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,UICE2N,VICE2N)
# endif
            
            CALL ICEOCEAN_STRESS
            UICE2RK  = UICE2
            VICE2RK  = VICE2
            CALL ICE_EVP(ADVSTRX,ADVSTRY,KID)
            DO K=1,4
               !CALL ICE_EVP(ADVSTRX,ADVSTRY)
               CALL ICE_UV(K)
               UICE2  = UICE2F
               VICE2  = VICE2F

!dwang
               !for global case: limit ice velocity of the new ice near the open water
               IF(K==4) then
                  DO I=1,N
                     UICE2(I)=MAX(MIN(UICE2(i),1.00_SP),-1.00_SP)
                     VICE2(I)=MAX(MIN(VICE2(i),1.00_SP),-1.00_SP)
                  END DO
               ENDIF

!dwang

               !EXCHANGE ELEMENT-BASED VALUES ACROSS THE INTERFACE
#     if defined (MULTIPROCESSOR)
               IF(PAR)CALL AEXCHANGE(EC,MYID,NPROCS,UICE2,VICE2)
#     endif
            END DO     !! END RUNGE-KUTTA LOOP
         END DO     !! END ice subcycle LOOP
         
         DO I=1,M
            SIG1(I) =SIG1(I) *ISICEN(I)
            SIG2(I) =SIG2(I) *ISICEN(I)
            SIG12(I) =SIG12(I) *ISICEN(I)
            ISICEN_OLD(I) =ISICEN(I)
         END DO
# if defined(MULTIPROCESSOR)
         IF(PAR)CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,1,MYID,NPROCS,SIG1,SIG2,SIG12)
         IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,SIG1,SIG2,SIG12)
# endif
         
      ELSE
         !==========================================================================
         Divu=0.0_SP;     Delta=0.0_SP
         DivuN=0.0_SP;    DeltaN=0.0_SP
         
         DO I=1,M
            IF(ISICEN(I) ==1) THEN
               DO I1=1, NTVE(I)
                  divuN(I) =divuN(I)+divuE(NBVE(I,I1))
                  DeltaN(I) =DeltaN(I)+DeltaE(NBVE(I,I1))
               END DO
               DivuN(I)  =DivuN(I)/float(NTVE(I))
               DeltaN(I) =DeltaN(I)/float(NTVE(I))
            ENDIF
         END DO
         
      ENDIF !! EVP selection

      Divu(1:M,1) =DivuN(1:M)
      Delta(1:M,1)=DeltaN(1:M)

      !gao for global case
      DO I=1,M
       !IF(I==NODE_NORTHPOLE) THEN
       IF(VY(I)>=89.5_SP) THEN
          Divu (I,1)=0.0_SP
          Delta(I,1)=0.0_SP
       END IF
      END DO


!==========================================================================
      DEALLOCATE(UICE2N,VICE2N)
      
      DEALLOCATE(AICETMP)
      DEALLOCATE(UICE2F,VICE2F)
      DEALLOCATE(UICE2RK,VICE2RK)
      DEALLOCATE(PSTXICE,PSTYICE)
      DEALLOCATE(ADVSTRX,ADVSTRY)
      DEALLOCATE(ADVUICE,ADVVICE)
      DEALLOCATE(IMASS1,IMASS,PICE)
      DEALLOCATE(TAUXWI,TAUYWI,vrel)
      DEALLOCATE(STRAIRXE,STRAIRYE)
   
      RETURN
    END SUBROUTINE ICE_RUNGE_KUTTA
!==========================================================================
!
!==============================================================================|
    SUBROUTINE ICEOCEAN_STRESS
      
      IMPLICIT NONE

      REAL(SP) :: WATERX,WATERY
      INTEGER  :: I
      
      DO I=1,N !T
         IF(ISICEC(I) == 1)THEN
            WATERX = U(I,1)*COSW - V(I,1)*SINW
            WATERY = V(I,1)*COSW + U(I,1)*SINW
            ! (magnitude of relative ocean current)*rhow*drag*aice
            VREL(I) = AIU(i)*DRAGW*SQRT((U(I,1)-UICE2(I))**2+(V(I,1)-VICE2(I))**2)  ! m/s
            ! ice/ocean stress
            TAUXWI(I) = VREL(I)*WATERX ! NOTE this is not the entire
            TAUYWI(I) = VREL(I)*WATERY ! ocn stress term
         END IF
      END DO
      
      RETURN
    END SUBROUTINE ICEOCEAN_STRESS
    !==============================================================================|
    !
    !==============================================================================|
    !     ACCUMLATE FLUXES FOR ICE                                                 |
    !==============================================================================|
    
    SUBROUTINE ICE_UV(K)       
      
      !==============================================================================|

#  if defined (SPHERICAL)
      !#  if defined (NORTHPOLE)
      USE MOD_NORTHPOLE
      !#  endif   
#  endif   

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: K
      REAL(SP), DIMENSION(0:NT) :: RESXICE,RESYICE
      REAL(SP) :: UICE2FT,VICE2FT
      REAL(SP) :: UI,VI
      INTEGER  :: I,J
      real(SP) :: WIND_SP
      real(SP) :: cca,ccb,ab2,cc1,cc2
      real(SP) :: TAUX, TAUY
      
      !==============================================================================|
      !
      !--ACCUMULATE RESIDUALS FOR ICE EQUATIONS--------------------------------------|
      !
      
      UICE2FT = UICE2F(0)
      VICE2FT = VICE2F(0)

      DO I=1,N
         IF(ISICEC(I) == 1)THEN
            !--UPDATE----------------------------------------------------------------------|
            ! alpha, beta are defined in Hunke and Dukowicz (1997), section 3.2
            cca = imass1(i)/Tdte + vrel(I) * cosw*ALPHA_RK(K) ! alpha, kg/m^2 s
            ccb = (cor(i)*imass1(i) + vrel(I) * sinw)*ALPHA_RK(K)  ! beta,  kg/m^2 s
            ab2 = cca**2 + ccb**2
            
            ! divergence of the internal stress tensor
            !  ADVSTRX(I),ADVSTRY(I)
            
            ! finally, the velocity components

!----------------J. Ge for baroclinic pressure gradient force----------------!
# if !defined (ICE_EMBEDDING)
            cc1 = (ADVSTRX(I) - IMASS1(I)*PSTXICE(I) + TAUXWI(I) +STRAIRXE(I))&
                 &       *ALPHA_RK(K) + imass1(I)/Tdte*UICE2RK(I)
            cc2 = (ADVSTRY(I) - IMASS1(I)*PSTYICE(I) + TAUYWI(I) +STRAIRYE(I))&
                 &        *ALPHA_RK(K)+ imass1(I)/Tdte*VICE2RK(I)

# else
            cc1 = (ADVSTRX(I) - IMASS1(I)*PSTXICE(I) + net_drhox + TAUXWI(I) +STRAIRXE(I))&
                 &       *ALPHA_RK(K) + imass1(I)/Tdte*UICE2RK(I)
            cc2 = (ADVSTRY(I) - IMASS1(I)*PSTYICE(I) + net_drhoy + TAUYWI(I) +STRAIRYE(I))&
                 &        *ALPHA_RK(K)+ imass1(I)/Tdte*VICE2RK(I)
# endif
!----------------J. Ge for baroclinic pressure gradient force----------------! 
           
            UICE2F(I) = (cca*cc1 + ccb*cc2)/ab2              ! m/s
            VICE2F(I) = (cca*cc2 - ccb*cc1)/ab2
            
         END IF  !  end isicec   
      END DO
      
#  if defined (SPHERICAL)
      CALL ICE_UV_XY(K)
#  endif
      
      DO I=1,N
         IF(ISBCE(I) == 1) THEN
            UI= UICE2F(I)*(SIN(ALPHA(I)))**2-VICE2F(I)*SIN(ALPHA(I))*COS(ALPHA(I))
            VI=-UICE2F(I)*SIN(ALPHA(I))*COS(ALPHA(I))+VICE2F(I)*(COS(ALPHA(I)))**2
!            UICE2F(I)=UI
!            VICE2F(I)=VI
         END IF
         UICE2F(I) =UICE2F(I)*ISICEC(I)
         VICE2F(I) =VICE2F(I)*ISICEC(I)
      END DO
      
      VICE2F(0) = VICE2FT
      UICE2F(0) = UICE2FT
      
      RETURN
    END SUBROUTINE ICE_UV
!==============================================================================|
!
!==============================================================================|

#  if defined (SPHERICAL)
    !#  if defined (NORTHPOLE)
    !
    !==============================================================================|
    SUBROUTINE ICE_EVP_XY(PSIG11PX,PSIG12PY,PSIG12PX,PSIG22PY,KI) 
      !------------------------------------------------------------------------------|
      !   USE BCS
      !   USE MOD_OBCS
      
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: KI
      REAL(SP), DIMENSION(0:MT) :: PUPX,PVPX,PUPY,PVPY
      REAL(SP) :: X11,Y11,X22,Y22,X33,Y33,TMP1,TMP2
      
      REAL(SP) :: E11,E12,E22
      Real(SP) :: DD_divu,DT_tension,DS_shear
      Real(SP) :: DELT,ZETA,ETA
      Real(SP) :: PDelta
      REAL(SP) :: SIG11IJ,SIG12IJ,SIG22IJ
      
      REAL(DP) :: ELIJ,UIJ,VIJ,EXFLUX
      INTEGER  :: I,J,II,K,I1,IA,IB,J1,J2,JTMP
      
      REAL(DP) :: UIJ_TMP,VIJ_TMP
      REAL(DP) :: VX1_TMP,VY1_TMP,VX2_TMP,VY2_TMP
      REAL(DP) :: XC_TMP,YC_TMP
      REAL(DP) :: DLTXE_TMP,DLTYE_TMP,DLTXC_TMP,DLTYC_TMP
      REAL(DP) :: STXJ1,STYJ1,STXJ2,STYJ2
      
      REAL(SP),DIMENSION(0:NT) :: PSIG11PX,PSIG12PY
      REAL(SP),DIMENSION(0:NT) :: PSIG12PX,PSIG22PY
      
      !================================================================================|
      !----------INITIALIZE STRESS ARRAY ----------------------------------------------!
      ! ALL OTHER PROCESSORS ESCAPE HERE
      IF (NODE_NORTHPOLE .EQ. 0) RETURN
      
      DO II=1,NPCV
         I = NCEDGE_LST(II)
         IA  = NIEC(I,1)
         IB  = NIEC(I,2)
         IF(IA == NODE_NORTHPOLE)THEN
            PUPX(IA) = 0.0_SP
            PVPX(IA) = 0.0_SP
            PUPY(IA) = 0.0_SP
            PVPY(IA) = 0.0_SP
         END IF
         IF(IB == NODE_NORTHPOLE)THEN  
            PUPX(IB) = 0.0_SP
            PVPX(IB) = 0.0_SP
            PUPY(IB) = 0.0_SP
            PVPY(IB) = 0.0_SP
         END IF
      END DO
      
      !---------ACCUMULATE STRESS BY LOOPING OVER CONTROL VOLUME HALF EDGES----------!
      
      DO II=1,NPCV
         I = NCEDGE_LST(II)
         I1  = NTRG(I)
         IA  = NIEC(I,1)
         IB  = NIEC(I,2)
         
         IF(ISICEN(IA) == 1 .OR. ISICEN(IB) == 1)THEN
            UIJ = UICE2(I1)
            VIJ = VICE2(I1)
            
            IF(IA == NODE_NORTHPOLE .OR. IB == NODE_NORTHPOLE)THEN
               UIJ_TMP = -VIJ*COS(XC(I1)*DEG2RAD)-UIJ*SIN(XC(I1)*DEG2RAD)
               VIJ_TMP = -VIJ*SIN(XC(I1)*DEG2RAD)+UIJ*COS(XC(I1)*DEG2RAD)
               
               VX1_TMP = REARTH * DCOS(YIJE(I,1)*DEG2RAD) * DCOS(XIJE(I,1)*DEG2RAD)&
                    * 2._DP /(1._DP+Dsin(YIJE(I,1)*DEG2RAD))
               VY1_TMP = REARTH * DCOS(YIJE(I,1)*DEG2RAD) * DSIN(XIJE(I,1)*DEG2RAD)&
                    * 2._DP /(1._DP+Dsin(YIJE(I,1)*DEG2RAD))
               
               VX2_TMP = REARTH * DCOS(YIJE(I,2)*DEG2RAD) * DCOS(XIJE(I,2)*DEG2RAD)&
                    * 2._DP /(1._DP+Dsin(YIJE(I,2)*DEG2RAD))
               VY2_TMP = REARTH * DCOS(YIJE(I,2)*DEG2RAD) * DSIN(XIJE(I,2)*DEG2RAD)&
                    * 2._DP /(1._DP+Dsin(YIJE(I,2)*DEG2RAD))
               
               DLTXE_TMP = VX2_TMP-VX1_TMP
               DLTYE_TMP = VY2_TMP-VY1_TMP
               
               EXFLUX = -UIJ_TMP*DLTYE_TMP                  
               IF(IA == NODE_NORTHPOLE) THEN    
                  PUPX(IA) = PUPX(IA)-EXFLUX/ART1(IA)
               ELSE IF(IB == NODE_NORTHPOLE)THEN
                  PUPX(IB) = PUPX(IB)+EXFLUX/ART1(IB)
               END IF

               EXFLUX = -VIJ_TMP*DLTYE_TMP                  
               IF(IA == NODE_NORTHPOLE) THEN    
                  PVPX(IA) = PVPX(IA)-EXFLUX/ART1(IA)
               ELSE IF(IB == NODE_NORTHPOLE)THEN
                  PVPX(IB) = PVPX(IB)+EXFLUX/ART1(IB)
               END IF
               
               EXFLUX =  UIJ_TMP*DLTXE_TMP  
               IF(IA == NODE_NORTHPOLE) THEN    
                  PUPY(IA) = PUPY(IA)-EXFLUX/ART1(IA)
               ELSE IF(IB == NODE_NORTHPOLE)THEN
                  PUPY(IB) = PUPY(IB)+EXFLUX/ART1(IB)
               END IF
               
               EXFLUX =  VIJ_TMP*DLTXE_TMP  
               IF(IA == NODE_NORTHPOLE) THEN    
                  PVPY(IA) = PVPY(IA)-EXFLUX/ART1(IA)
               ELSE IF(IB == NODE_NORTHPOLE)THEN
                  PVPY(IB) = PVPY(IB)+EXFLUX/ART1(IB)
               END IF
            END IF
         END IF
      END DO
!!!==================================================================================================
      I = NODE_NORTHPOLE

      IF(.false.) then
         PUPX(I)=0.0_SP
         PUPY(I)=0.0_SP
         PVPX(I)=0.0_SP
         PVPY(I)=0.0_SP
         
         DO J=1,NTVE(I)
            I1=NBVE(I,J)
            JTMP=NBVT(I,J)
            J1=JTMP+1-(JTMP+1)/4*3
            J2=JTMP+2-(JTMP+2)/4*3
            !       X11=0.5_SP*(VX(I)+VX(NV(I1,J1)))
            !       Y11=0.5_SP*(VY(I)+VY(NV(I1,J1)))
            VX1_TMP = REARTH * DCOS(VY(I)*DEG2RAD) * DCOS(VX(I)*DEG2RAD)&
                 * 2._DP /(1._DP+Dsin(VY(I)*DEG2RAD))
            VY1_TMP = REARTH * DCOS(VY(I)*DEG2RAD) * DSIN(VX(I)*DEG2RAD)&
                 * 2._DP /(1._DP+Dsin(VY(I)*DEG2RAD))
            
            VX2_TMP = REARTH * DCOS(VY(NV(I1,J1))*DEG2RAD) * DCOS(VX(NV(I1,J1))*DEG2RAD)&
                 * 2._DP /(1._DP+Dsin(VY(NV(I1,J1))*DEG2RAD))
            VY2_TMP = REARTH * DCOS(VY(NV(I1,J1))*DEG2RAD) * DSIN(VX(NV(I1,J1))*DEG2RAD)&
                 * 2._DP /(1._DP+Dsin(VY(NV(I1,J1))*DEG2RAD))
            
            X11=0.5_SP*(VX1_TMP+VX2_TMP)
            Y11=0.5_SP*(VY1_TMP+VY2_TMP)
            
            !       X22=XC(I1)
            !       Y22=YC(I1)
            
            XC_TMP = REARTH * DCOS(YC(I1)*DEG2RAD) * DCOS(XC(I1)*DEG2RAD)&
                 * 2._DP /(1._DP+Dsin(YC(I)*DEG2RAD))
            YC_TMP = REARTH * DCOS(YC(I)*DEG2RAD) * DSIN(XC(I1)*DEG2RAD)&
                 * 2._DP /(1._DP+Dsin(YC(I)*DEG2RAD))
            
            X22=XC_TMP
            Y22=YC_TMP
            
            !       X33=0.5_SP*(VX(I)+VX(NV(I1,J2)))
            !       Y33=0.5_SP*(VY(I)+VY(NV(I1,J2)))
            
            VX1_TMP = REARTH * DCOS(VY(I)*DEG2RAD) * DCOS(VX(I)*DEG2RAD)&
                 * 2._DP /(1._DP+Dsin(VY(I)*DEG2RAD))
            VY1_TMP = REARTH * DCOS(VY(I)*DEG2RAD) * DSIN(VX(I)*DEG2RAD)&
                 * 2._DP /(1._DP+Dsin(VY(I)*DEG2RAD))
            
            VX2_TMP = REARTH * DCOS(VY(NV(I1,J2))*DEG2RAD) * DCOS(VX(NV(I1,J2))*DEG2RAD)&
                 * 2._DP /(1._DP+Dsin(VY(NV(I1,J2))*DEG2RAD))
            VY2_TMP = REARTH * DCOS(VY(NV(I1,J2))*DEG2RAD) * DSIN(VX(NV(I1,J2))*DEG2RAD)&
                 * 2._DP /(1._DP+Dsin(VY(NV(I1,J2))*DEG2RAD))

            X33=0.5_SP*(VX1_TMP+VX2_TMP)
            Y33=0.5_SP*(VY1_TMP+VY2_TMP)
            
            UIJ = UICE2(I1)
            VIJ = VICE2(I1)
            UIJ_TMP = -VIJ*COS(XC(I1)*DEG2RAD)-UIJ*SIN(XC(I1)*DEG2RAD)
            VIJ_TMP = -VIJ*SIN(XC(I1)*DEG2RAD)+UIJ*COS(XC(I1)*DEG2RAD)
            
            PUPX(I)=PUPX(I)+UIJ_TMP*(Y11-Y33)
            PUPY(I)=PUPY(I)+UIJ_TMP*(X33-X11)
            PVPX(I)=PVPX(I)+VIJ_TMP*(Y11-Y33)
            PVPY(I)=PVPY(I)+VIJ_TMP*(X33-X11)
            
         END DO
         
         PUPX(I) = PUPX(I)/ART1(I)
         PVPX(I) = PVPX(I)/ART1(I)
         PUPY(I) = PUPY(I)/ART1(I)
         PVPY(I) = PVPY(I)/ART1(I)
      END IF !! 
!!!==================================================================================================
      
      IF(ISICEN(I) == 0)RETURN
      
      E11 = PUPX(I)
      E12 = 0.5_SP*(PUPY(I)+PVPX(I))
      E22 = PVPY(I)
      
      ! divergence  =  e_11 + e_22
      DD_divu = E11+E22
      ! tension strain rate  =  e_11 - e_22
      DT_tension = E11-E22
      ! shearing strain rate  =  e_12
      DS_shear = 2.0_SP*E12  
      
      DELT = DD_divu*DD_divu+ECCI*(DT_tension*DT_tension+DS_shear*DS_shear)
      DELT = SQRT(DELT)
      !  DELT= min(max(DELT,2.E-9),PICE(I)/(2.0*4.E+8))
      !  DELT= min(max(DELT,2.E-8),PICE(I)/(2.0*4.E+8))
      !  For 1D case
      !--------------------------------------------------  
      !  calculate stress for ridging
      
      DELT= max(DELT,2.E-9)
      IF (KI == ndte) THEN 
         DeltaN (I)= DELT 
         divuN  (I)= DD_divu 
      ENDIF
      !DELT= min(max(DELT,2.E-8),PICE(I)/(2.0*4.E+8))
      !--------------------------------------------------  
      !       Pressure/Delta
      PDelta = PICE(I)*1.0E-6/DELT
      !Pdelta =min(Pice(I)*1.E-6/Delt,Rcon*Art1(I)*1.0E-6)
      !!=========================================================================================| 
      !     first step stress equations
      !     computing sigma1 and sigma2   
      !     SIG1 (I) = ( SIG1 (I) + dte2T*(PDelta*DD_divu-Pice(I))*1.0E+6)  * denom1 
      SIG1 (I) = ( SIG1 (I) + dte2T*(PDelta*(DD_divu-DELT))*1.0E+6)*denom1 
      SIG2 (I) = ( SIG2 (I) + dte2T*PDelta*DT_tension      *1.0E+6)*denom2 
      SIG12(I) = ( SIG12(I) + dte2T*Pdelta*0.5_SP*DS_shear *1.0E+6)*denom2
      
      !     recover sigma11 and sigma22
      SIG11(I)=(SIG1(I)+SIG2(I))*0.5_SP
      SIG22(I)=(SIG1(I)-SIG2(I))*0.5_SP
      
      !!=========================================================================================|
      prs_sig(I)=PDelta*DELT*1.E+6
      
      DO II=1,NPE
         I = NPEDGE_LST(II)
         IA = IEC(I,1)
         IB = IEC(I,2)
         
         IF(CELL_NORTHAREA(IA) == 1)THEN
            PSIG11PX(IA) = 0.0_SP
            PSIG12PY(IA) = 0.0_SP
            
            PSIG12PX(IA) = 0.0_SP
            PSIG22PY(IA) = 0.0_SP 
            
            PSTXICE(IA) = 0.0_SP
            PSTYICE(IA) = 0.0_SP
         END IF
         IF(CELL_NORTHAREA(IB) == 1)THEN  
            PSIG11PX(IB) = 0.0_SP
            PSIG12PY(IB) = 0.0_SP
            
            PSIG12PX(IB) = 0.0_SP
            PSIG22PY(IB) = 0.0_SP 
            
            PSTXICE(IB) = 0.0_SP
            PSTYICE(IB) = 0.0_SP
         END IF
      END DO
      
      DO II=1,NPE
         I=NPEDGE_LST(II)
         IA=IEC(I,1)
         IB=IEC(I,2)
         J1=IENODE(I,1)
         J2=IENODE(I,2)
         
         !     IF(ISICEC(IA) == 1 .OR. ISICEC(IB) == 1)THEN
         ELIJ=0.5_SP*(EL(J1)+EL(J2))
         ! no surface tilt contribution ggao 0516
         ELIJ=0.0
         
         SIG11IJ=0.5_SP*(SIG11(J1)+SIG11(J2))
         SIG22IJ=0.5_SP*(SIG22(J1)+SIG22(J2))
         SIG12IJ=0.5_SP*(SIG12(J1)+SIG12(J2))
         
         VX1_TMP = REARTH * COS(VY(IENODE(I,1))*DEG2RAD) * COS(VX(IENODE(I,1))*DEG2RAD) &
              * 2._DP /(1._DP+sin(VY(IENODE(I,1))*DEG2RAD))
         VY1_TMP = REARTH * COS(VY(IENODE(I,1))*DEG2RAD) * SIN(VX(IENODE(I,1))*DEG2RAD) &
              * 2._DP /(1._DP+sin(VY(IENODE(I,1))*DEG2RAD))
         
         VX2_TMP = REARTH * COS(VY(IENODE(I,2))*DEG2RAD) * COS(VX(IENODE(I,2))*DEG2RAD) &
              * 2._DP /(1._DP+sin(VY(IENODE(I,2))*DEG2RAD))
         VY2_TMP = REARTH * COS(VY(IENODE(I,2))*DEG2RAD) * SIN(VX(IENODE(I,2))*DEG2RAD) &
              * 2._DP /(1._DP+sin(VY(IENODE(I,2))*DEG2RAD))

         DLTXC_TMP = VX2_TMP-VX1_TMP
         DLTYC_TMP = VY2_TMP-VY1_TMP
         
         EXFLUX    = SIG11IJ*DLTYC_TMP
         IF(CELL_NORTHAREA(IA) == 1 .AND. CELL_NORTHAREA(IB) == 1)THEN
            PSIG11PX(IA) = PSIG11PX(IA) - EXFLUX/ART(IA)
            PSIG11PX(IB) = PSIG11PX(IB) + EXFLUX/ART(IB)
         ELSE IF(CELL_NORTHAREA(IA) == 1 .AND. CELL_NORTHAREA(IB) /= 1)THEN
            PSIG11PX(IA) = PSIG11PX(IA) - EXFLUX/ART(IA)
         ELSE IF(CELL_NORTHAREA(IB) == 1 .AND. CELL_NORTHAREA(IA) /= 1)THEN  
            PSIG11PX(IB) = PSIG11PX(IB) + EXFLUX/ART(IB)
         END IF
         
         EXFLUX    = SIG12IJ*DLTXC_TMP
         IF(CELL_NORTHAREA(IA) == 1 .AND. CELL_NORTHAREA(IB) == 1)THEN
            PSIG12PY(IA) = PSIG12PY(IA) + EXFLUX/ART(IA)
            PSIG12PY(IB) = PSIG12PY(IB) - EXFLUX/ART(IB)
         ELSE IF(CELL_NORTHAREA(IA) == 1 .AND. CELL_NORTHAREA(IB) /= 1)THEN
            PSIG12PY(IA) = PSIG12PY(IA) + EXFLUX/ART(IA)
         ELSE IF(CELL_NORTHAREA(IB) == 1 .AND. CELL_NORTHAREA(IA) /= 1)THEN  
            PSIG12PY(IB) = PSIG12PY(IB) - EXFLUX/ART(IB)
         END IF

         EXFLUX    = SIG12IJ*DLTYC_TMP
         IF(CELL_NORTHAREA(IA) == 1 .AND. CELL_NORTHAREA(IB) == 1)THEN
            PSIG12PX(IA) = PSIG12PX(IA) - EXFLUX/ART(IA)
            PSIG12PX(IB) = PSIG12PX(IB) + EXFLUX/ART(IB)
         ELSE IF(CELL_NORTHAREA(IA) == 1 .AND. CELL_NORTHAREA(IB) /= 1)THEN
            PSIG12PX(IA) = PSIG12PX(IA) - EXFLUX/ART(IA)
         ELSE IF(CELL_NORTHAREA(IB) == 1 .AND. CELL_NORTHAREA(IA) /= 1)THEN  
            PSIG12PX(IB) = PSIG12PX(IB) + EXFLUX/ART(IB)
         END IF

         EXFLUX    = SIG22IJ*DLTXC_TMP
         IF(CELL_NORTHAREA(IA) == 1 .AND. CELL_NORTHAREA(IB) == 1)THEN
            PSIG22PY(IA) = PSIG22PY(IA) + EXFLUX/ART(IA)
            PSIG22PY(IB) = PSIG22PY(IB) - EXFLUX/ART(IB)
         ELSE IF(CELL_NORTHAREA(IA) == 1 .AND. CELL_NORTHAREA(IB) /= 1)THEN
            PSIG22PY(IA) = PSIG22PY(IA) + EXFLUX/ART(IA)
         ELSE IF(CELL_NORTHAREA(IB) == 1 .AND. CELL_NORTHAREA(IA) /= 1)THEN  
            PSIG22PY(IB) = PSIG22PY(IB) - EXFLUX/ART(IB)
         END IF
         
         !      ACCUMULATE BAROTROPIC FLUX
         IF(CELL_NORTHAREA(IA) == 1 .AND. CELL_NORTHAREA(IB) == 1)THEN
            PSTXICE(IA)=PSTXICE(IA)-GRAV_E(IA)*ELIJ*(DLTYC_TMP/ART(IA))/          &
                 (2._SP /(1._SP+sin(YC(IA)*DEG2RAD)))
            PSTYICE(IA)=PSTYICE(IA)+GRAV_E(IA)*ELIJ*(DLTXC_TMP/ART(IA))/           &
                 (2._SP /(1._SP+sin(YC(IA)*DEG2RAD)))
            PSTXICE(IB)=PSTXICE(IB)+GRAV_E(IB)*ELIJ*(DLTYC_TMP/ART(IB))/           &
                 (2._SP /(1._SP+sin(YC(IB)*DEG2RAD)))
            PSTYICE(IB)=PSTYICE(IB)-GRAV_E(IB)*ELIJ*(DLTXC_TMP/ART(IB))/           &
                 (2._SP /(1._SP+sin(YC(IB)*DEG2RAD)))
         ELSE IF(CELL_NORTHAREA(IA) == 1 .AND. CELL_NORTHAREA(IB) /= 1)THEN
            PSTXICE(IA)=PSTXICE(IA)-GRAV_E(IA)*ELIJ*(DLTYC_TMP/ART(IA))/           &
                 (2._SP /(1._SP+sin(YC(IA)*DEG2RAD)))
            PSTYICE(IA)=PSTYICE(IA)+GRAV_E(IA)*ELIJ*(DLTXC_TMP/ART(IA))/           &
                 (2._SP /(1._SP+sin(YC(IA)*DEG2RAD)))
         ELSE IF(CELL_NORTHAREA(IB) == 1 .AND. CELL_NORTHAREA(IA) /= 1)THEN  
            PSTXICE(IB)=PSTXICE(IB)+GRAV_E(IB)*ELIJ*(DLTYC_TMP/ART(IB))/           &
                 (2._SP /(1._SP+sin(YC(IB)*DEG2RAD)))
            PSTYICE(IB)=PSTYICE(IB)-GRAV_E(IB)*ELIJ*(DLTXC_TMP/ART(IB))/           &
                 (2._SP /(1._SP+sin(YC(IB)*DEG2RAD)))
         END IF
         !     END IF  !  is ice
      END DO
      
      RETURN
    END SUBROUTINE ICE_EVP_XY
    !==============================================================================|
    !
    !==============================================================================|
    !     ACCUMLATE FLUXES FOR ICE                                                 |
    !==============================================================================|
    
    SUBROUTINE ICE_UV_XY(K)       
      
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: K
      REAL(SP), DIMENSION(0:NT) :: RESXICE,RESYICE
      REAL(SP) :: UICE2FT,VICE2FT
      INTEGER  :: I,II
      
      REAL(SP) :: UICE_TMP,VICE_TMP,STRAIRX_TMP,STRAIRY_TMP
      REAL(SP) :: PSTXICE_TMP,PSTYICE_TMP
      REAL(SP) :: UICE2RK_TMP,VICE2RK_TMP,UICE2F_TMP,VICE2F_TMP
      REAL(SP) :: DLTYC_TMP ,DLTXC_TMP

      real(SP) :: cca,ccb,ab2,cc1,cc2
      real(SP) :: UW_TMP,VW_TMP
      real(SP) :: TAUXWI_tmp,TAUYWI_tmp
      !==============================================================================|
      !
      !--ACCUMULATE RESIDUALS FOR ICE EQUATIONS--------------------------------------|
      !
      ! ALL OTHER PROCESSORS ESCAPE HERE
      IF (NODE_NORTHPOLE .EQ. 0) RETURN
      
      DO II=1,NP
         I=NP_LST(II)
         IF(ISICEC(I) == 1)THEN
            UICE_TMP = -VICE2(I)*COS(XC(I)*DEG2RAD)-UICE2(I)*SIN(XC(I)*DEG2RAD)
            VICE_TMP = -VICE2(I)*SIN(XC(I)*DEG2RAD)+UICE2(I)*COS(XC(I)*DEG2RAD)
            
            UW_TMP = -V(I,1)*COS(XC(I)*DEG2RAD)-U(I,1)*SIN(XC(I)*DEG2RAD)
            VW_TMP = -V(I,1)*SIN(XC(I)*DEG2RAD)+U(I,1)*COS(XC(I)*DEG2RAD)
            
            TAUXWI_TMP = -TAUYWI(I)*COS(XC(I)*DEG2RAD)-TAUXWI(I)*SIN(XC(I)*DEG2RAD)
            TAUYWI_TMP = -TAUYWI(I)*SIN(XC(I)*DEG2RAD)+TAUXWI(I)*COS(XC(I)*DEG2RAD)
            
            STRAIRX_TMP = -STRAIRYE(I)*COS(XC(I)*DEG2RAD)-STRAIRXE(I)*SIN(XC(I)*DEG2RAD)
            STRAIRY_TMP = -STRAIRYE(I)*SIN(XC(I)*DEG2RAD)+STRAIRXE(I)*COS(XC(I)*DEG2RAD)
            
            PSTXICE_TMP = -PSTYICE(I)*COS(XC(I)*DEG2RAD)-PSTXICE(I)*SIN(XC(I)*DEG2RAD)
            PSTYICE_TMP = -PSTYICE(I)*SIN(XC(I)*DEG2RAD)+PSTXICE(I)*COS(XC(I)*DEG2RAD)
            !
            !--UPDATE----------------------------------------------------------------------|
            !
            ! alpha, beta are defined in Hunke and Dukowicz (1997), section 3.2
            !        cca = imass1(i)/dtice + vrel * cosw         ! alpha, kg/m^2 s
            cca = imass1(i)/Tdte + vrel(I) * cosw*ALPHA_RK(K) ! alpha, kg/m^2 s
            !        ccb = fm(i,j)        + vrel * sinw        ! beta,  kg/m^2 s
            ccb = (cor(i)*imass1(i) + vrel(I) * sinw)*ALPHA_RK(K)  ! beta,  kg/m^2 s
            
            ab2 = cca**2 + ccb**2
            
            ! divergence of the internal stress tensor
            !  ADVSTRX(I),ADVSTRY(I)

            UICE2RK_TMP = -VICE2RK(I)*COS(XC(I)*DEG2RAD)-UICE2RK(I)*SIN(XC(I)*DEG2RAD)
            VICE2RK_TMP = -VICE2RK(I)*SIN(XC(I)*DEG2RAD)+UICE2RK(I)*COS(XC(I)*DEG2RAD)
            
            ! finally, the velocity components
            cc1 = (ADVSTRX(I) - IMASS1(I)*PSTXICE_TMP + TAUXWI_TMP+STRAIRX_TMP)&
                 &       *ALPHA_RK(K) + imass1(I)/Tdte*UICE2RK_TMP
            cc2 = (ADVSTRY(I) - IMASS1(I)*PSTYICE_TMP + TAUYWI_TMP+STRAIRY_TMP)&
                 &        *ALPHA_RK(K)+ imass1(I)/Tdte*VICE2RK_TMP
            
            UICE2F_TMP = (cca*cc1 + ccb*cc2)/ab2              ! m/s
            VICE2F_TMP = (cca*cc2 - ccb*cc1)/ab2
            
            UICE2F(I)  = VICE2F_TMP*COS(XC(I)*DEG2RAD)-UICE2F_TMP*SIN(XC(I)*DEG2RAD)
            VICE2F(I)  = UICE2F_TMP*COS(XC(I)*DEG2RAD)+VICE2F_TMP*SIN(XC(I)*DEG2RAD)
            VICE2F(I)  = -VICE2F(I)

         else 
            UICE2F(I)=0.0_SP
            VICE2F(I)=0.0_SP
         endif
         
      END DO
      
      RETURN
    END SUBROUTINE ICE_UV_XY
    !==============================================================================|
    !==============================================================================|
    !#  endif
    !end if defined (NORTHPOLE)
#  endif
    !end if defined (SPHERICAL)
    
    
    !==============================================================================|
    !   Calculate Advection  for ICE variables                                     |
    !==============================================================================|
    
    SUBROUTINE ICE_ADV_fvcom(var_ice,cuice,cvice)       
      
      !------------------------------------------------------------------------------|
      IMPLICIT NONE
      
      REAL(SP), DIMENSION(0:MT)    :: XFLUX
      REAL(SP), DIMENSION(0:MT)    :: PUPX,PUPY,PVPX,PVPY  
      REAL(DP), DIMENSION(0:MT)    :: PIPX,PIPY
      REAL(SP), DIMENSION(3*(NT))  :: DTIJ 
      REAL(SP), DIMENSION(3*(NT))  :: UVN
      
      REAL(SP), DIMENSION(0:NT)  :: cuice,cvice
      REAL(SP), DIMENSION(0:MT)  :: var_ice,var_ice1
      
      REAL(SP) :: UTMP,VTMP,SITAI,FFD,FF1,X11,Y11,X22,Y22,X33,Y33,TMP1,TMP2,XI,YI
      REAL(SP) :: DXA,DYA,DXB,DYB,FIJ1,FIJ2,UN,TTIME,ZDEP
      REAL(SP) :: TXX,TYY,FXX,FYY,VISCOF,EXFLUX,TEMP,STPOINT,STPOINT1,STPOINT2
      REAL(SP) :: FACT,FM1
      INTEGER  :: I,I1,I2,IA,IB,J,J1,J2,K,JTMP,JJ,II
#  if defined (SPHERICAL)
      REAL(DP) :: TY,TXPI,TYPI
      REAL(DP) :: XTMP1,XTMP
      REAL(DP) :: X1_DP,Y1_DP,X2_DP,Y2_DP,XII,YII
      REAL(DP) :: X11_TMP,Y11_TMP,X33_TMP,Y33_TMP
      !# if defined (NORTHPOLE)
      REAL(DP) :: VX1_TMP,VY1_TMP,VX2_TMP,VY2_TMP
      REAL(DP) :: TXPI_TMP,TYPI_TMP
      !# endif
#  endif
      
      REAL(SP) :: VI1MIN, VI1MAX, VI2MIN, VI2MAX
      REAL(SP) :: XMIN, XMAX
      REAL(SP) :: ART_TMP 
      
      
      !------------------------------------------------------------------------------!
      !
      !--Initialize Fluxes-----------------------------------------------------------!
      !
      XFLUX     = 0.0_SP
      UVN       =0.0_SP
      !
      !--Loop Over Control Volume Sub-Edges And Calculate Normal Velocity------------!
      !
      DO I=1,NCV
         I1=NTRG(I)
         UVN(I) = cvice(I1)*DLTXE(I) - cuice(I1)*DLTYE(I)
      END DO
      
      !
      !--Calculate the Advection and Horizontal Diffusion Terms----------------------!
      !
      
      PIPX  = 0.0_SP
      PIPY  = 0.0_SP
      DO I=1,M
         DO J=1,NTSN(I)-1
            I1=NBSN(I,J)
            I2=NBSN(I,J+1)
            FFD=0.5_SP*(var_ice(I1)+var_ice(I2))
            FF1=0.5_SP*(var_ice(I1)+var_ice(I2))
            
            PIPX(I)=PIPX(I)+FF1*DLTYTRIE(i,j)
            PIPY(I)=PIPY(I)+FF1*DLTXTRIE(i,j)
            
         END DO
         PIPX(I)=PIPX(I)/ART2(I)
         PIPY(I)=PIPY(I)/ART2(I)
         
      END DO
      
      DO I=1,NCV_I
         
         I1=NTRG(I)
         
         IF(ISICEC(I1)==1) THEN
         UVN(I) = cvice(I1)*DLTXE(I) - cuice(I1)*DLTYE(I)
         
         IA=NIEC(I,1)
         IB=NIEC(I,2)

         FIJ1=var_ice(IA)+DLTXNCVE(I,1)*PIPX(IA)+DLTYNCVE(I,1)*PIPY(IA)
         FIJ2=var_ice(IB)+DLTXNCVE(I,2)*PIPX(IB)+DLTYNCVE(I,2)*PIPY(IB)
         
         VI1MIN=MINVAL(var_ice(NBSN(IA,1:NTSN(IA)-1)))
         VI1MIN=MIN(VI1MIN, var_ice(IA))
         VI1MAX=MAXVAL(var_ice(NBSN(IA,1:NTSN(IA)-1)))
         VI1MAX=MAX(VI1MAX, var_ice(IA))
         VI2MIN=MINVAL(var_ice(NBSN(IB,1:NTSN(IB)-1)))
         VI2MIN=MIN(VI2MIN, var_ice(IB))
         VI2MAX=MAXVAL(var_ice(NBSN(IB,1:NTSN(IB)-1)))
         VI2MAX=MAX(VI2MAX, var_ice(IB))
         IF(FIJ1 < VI1MIN) FIJ1=VI1MIN
         IF(FIJ1 > VI1MAX) FIJ1=VI1MAX
         IF(FIJ2 < VI2MIN) FIJ2=VI2MIN
         IF(FIJ2 > VI2MAX) FIJ2=VI2MAX
         
         UN=UVN(I)
         
         EXFLUX=-UN* &
              ((1.0_SP+SIGN(1.0_SP,UN))*FIJ2+(1.0_SP-SIGN(1.0_SP,UN))*FIJ1)*0.5_SP 

         IF(ISICEN(IA)==0 .and. EXFLUX >0.0_SP) EXFLUX=0.0
         IF(ISICEN(IB)==0 .and. EXFLUX <0.0_SP) EXFLUX=0.0
         
         XFLUX(IA)=XFLUX(IA)+EXFLUX
         XFLUX(IB)=XFLUX(IB)-EXFLUX
         
         END IF  ! just calculate the ice
      END DO
      
#    if defined (SPHERICAL)
      !#    if defined (NORTHPOLE)
      CALL ICE_ADV_XY(XFLUX,PIPX,PIPY,var_ice,cuice,cvice)
      !#    endif
#    endif  
      
# if defined (MULTIPROCESSOR)
      IF(PAR)CALL NODE_MATCH(0,NBN,BN_MLT,BN_LOC,BNC,MT,1,MYID,NPROCS,XFLUX)
# endif
      
      !--Update ice variables----------------------------------------------------------!
      !
     
    var_ICE1=0.0_SP
# if defined (ICE_FRESHWATER)    
    DO I=1,M
! afm 20151207 & EJA 20160921 - Fix bug.
! 1. ART_TMP: Underestimates a control volume (It only counts the ice-covered
!             cell), resulting in overestimating flux.
! 2. "IF(ART_TMP>0.0)&...": Ice only advects to a cell where ice already exists.
! Put them into a classical advection form: var_ICE(I)=var_ice(I)-XFLUX(I)*dyn_dt/ART1(i)
!
!!       IF(ISICEN(I)==1)THEN
!       ART_TMP=0.0
!       DO J=1,NTVE(I)
!!          ART_TMP=ART_TMP+ISICEC(NBVE(I,J))*ART(NBVE(I,J))
!          ART_TMP=ART_TMP+real(ISICEC(NBVE(I,J)),SP)*ART(NBVE(I,J))
!       END DO
!        IF(ART_TMP>0.0)&
!          var_ICE1(I)=var_ice(I)-XFLUX(I)*dyn_dt/(ART_TMP/3.0_SP)
!!       ENDIF
!       var_ICE(I)=var_ICE1(I)
!       if (VAR_ICE(I) < 0.0)VAR_ICE(I) = c0i
        var_ICE(I)=var_ice(I)-XFLUX(I)*dyn_dt/ART1(i)
    END DO
# else
    DO I=1,M
       IF(ISICEN(I)==1)THEN
       ART_TMP=0.0
       DO J=1,NTVE(I)
          ART_TMP=ART_TMP+ISICEC(NBVE(I,J))*ART(NBVE(I,J))
       END DO
        IF(ART_TMP>0.0)&
          var_ICE1(I)=var_ice(I)-XFLUX(I)*dyn_dt/(ART_TMP/3.0_SP)
       ENDIF
       var_ICE(I)=var_ICE1(I)
       if (VAR_ICE(I) < 0.0)VAR_ICE(I) = c0i
    END DO
# endif      
      
! - begin afm 20171113 for ice open boundary -------------
      XFLUX_ICE = XFLUX
! - end afm 20171113 for ice open boundary -------------
      RETURN
    END SUBROUTINE ICE_ADV_fvcom
    !==============================================================================|
    
    !#  if defined (NORTHPOLE)
#  if defined (SPHERICAL)
    
    !==============================================================================|
    SUBROUTINE ICE_ADV_XY(XFLUX,PIPX,PIPY,var_ice,cuice,cvice)               
      !------------------------------------------------------------------------------|
      
      IMPLICIT NONE
      REAL(SP), DIMENSION(0:MT)     :: XFLUX,XFLUX_ADV
      REAL(DP), DIMENSION(0:MT)     :: PIPX,PIPY
      REAL(SP), DIMENSION(3*(NT))   :: DTIJ 
      REAL(SP) :: XI,YI
      REAL(SP) :: DXA,DYA,DXB,DYB,FIJ1,FIJ2 
      REAL(SP) :: TXX,TYY,FXX,FYY,VISCOF   
      REAL(SP) :: FACT,FM1
      INTEGER  :: I,I1,I2,IA,IB,J,J1,J2,K,JTMP,JJ,II
      REAL(SP) :: TXPI,TYPI
      
      REAL(SP) :: VX_TMP,VY_TMP,VX1_TMP,VY1_TMP,VX2_TMP,VY2_TMP,VX3_TMP,VY3_TMP
      REAL(SP) :: XI_TMP,YI_TMP,VXA_TMP,VYA_TMP,VXB_TMP,VYB_TMP
      REAL(SP) :: UIJ_TMP,VIJ_TMP,DLTXE_TMP,DLTYE_TMP,UVN_TMP,EXFLUX_TMP
      REAL(SP) :: PUPX_TMP,PUPY_TMP,PVPX_TMP,PVPY_TMP
      REAL(SP) :: PIPX_TMP,PIPY_TMP
      REAL(SP) :: U_TMP,V_TMP
      REAL(SP) :: X11,Y11,X22,Y22,X33,Y33,TMP1,TMP2
      integer :: KK
      
      REAL(SP), DIMENSION(0:NT)  :: cuice,cvice
      REAL(SP), DIMENSION(0:MT)  :: var_ice,var_ice1
      !!  ggao edge calculation
      REAL(SP) :: XIJE1_TMP,YIJE1_TMP,XIJE2_TMP,YIJE2_TMP
      REAL(SP) :: VI1MIN, VI1MAX, VI2MIN, VI2MAX
      
      !------------------------------------------------------------------------------!
      !
      !--Initialize Fluxes-----------------------------------------------------------!
      !
      ! ALL OTHER PROCESSORS ESCAPE HERE
      IF (NODE_NORTHPOLE .EQ. 0) RETURN
      
      DO II=1,NPCV
         I = NCEDGE_LST(II)
         IA = NIEC(I,1)
         IB = NIEC(I,2)
         IF(IA == NODE_NORTHPOLE)THEN
            XFLUX(IA) = 0.0_SP
            XFLUX_ADV(IA) = 0.0_SP
         ELSE IF(IB == NODE_NORTHPOLE)THEN  
            XFLUX(IB) = 0.0_SP
            XFLUX_ADV(IB) = 0.0_SP
         END IF
      END DO
      
      !
      !--Loop Over Control Volume Sub-Edges And Calculate Normal Velocity------------!
      !
      DO II=1,NPCV
         I = NCEDGE_LST(II)
         I1=NTRG(I)
      END DO

      !
      !--Calculate the Advection and Horizontal Diffusion Terms----------------------!
      !
      I = NODE_NORTHPOLE
      
      IF(I==0)  RETURN
      
      DO II=1,NPCV
         I = NCEDGE_LST(II)
         I1=NTRG(I)
         IA=NIEC(I,1)
         IB=NIEC(I,2)
         
         IF((IA <= M .AND. IB <= M) .AND. I1 <= N)THEN

            !       XI=0.5_SP*(XIJE(I,1)+XIJE(I,2))
            !       YI=0.5_SP*(YIJE(I,1)+YIJE(I,2))
            !!  ggao edge calculation
            XIJE1_TMP = REARTH * COS(YIJE(I,1)*DEG2RAD) * COS(XIJE(I,1)*DEG2RAD) &
                 * 2._SP /(1._SP+sin(YIJE(I,1)*DEG2RAD))
            YIJE1_TMP = REARTH * COS(YIJE(I,1)*DEG2RAD) * SIN(XIJE(I,1)*DEG2RAD) &
                 * 2._SP /(1._SP+sin(YIJE(I,1)*DEG2RAD))
            
            XIJE2_TMP = REARTH * COS(YIJE(I,2)*DEG2RAD) * COS(XIJE(I,2)*DEG2RAD) &
                 * 2._SP /(1._SP+sin(YIJE(I,2)*DEG2RAD))
            YIJE2_TMP = REARTH * COS(YIJE(I,2)*DEG2RAD) * SIN(XIJE(I,2)*DEG2RAD) &
                 * 2._SP /(1._SP+sin(YIJE(I,2)*DEG2RAD))
            XI_TMP =0.5_SP*(XIJE1_TMP+XIJE2_TMP)
            YI_TMP =0.5_SP*(YIJE1_TMP+YIJE2_TMP)
            
            !!  change end 0220/2008

            IF(IA == NODE_NORTHPOLE .OR. IB == NODE_NORTHPOLE)THEN
               VXA_TMP = REARTH * COS(VY(IA)*DEG2RAD) * COS(VX(IA)*DEG2RAD) &
                    * 2._SP /(1._SP+sin(VY(IA)*DEG2RAD))
               VYA_TMP = REARTH * COS(VY(IA)*DEG2RAD) * SIN(VX(IA)*DEG2RAD) &
                    * 2._SP /(1._SP+sin(VY(IA)*DEG2RAD))
               
               VXB_TMP = REARTH * COS(VY(IB)*DEG2RAD) * COS(VX(IB)*DEG2RAD) &
                    * 2._SP /(1._SP+sin(VY(IB)*DEG2RAD))
               VYB_TMP = REARTH * COS(VY(IB)*DEG2RAD) * SIN(VX(IB)*DEG2RAD) &
                    * 2._SP /(1._SP+sin(VY(IB)*DEG2RAD))
               
               !         IF(IA == NODE_NORTHPOLE)THEN
               DXA=XI_TMP-VXA_TMP
               DYA=YI_TMP-VYA_TMP
               !         ELSE IF(IB == NODE_NORTHPOLE)THEN
               DXB=XI_TMP-VXB_TMP
               DYB=YI_TMP-VYB_TMP
               !	 END IF
               !       END IF
               
               IF(IA == NODE_NORTHPOLE)THEN
                  PIPX_TMP=-PIPY(IB)*COS(VX(IB)*DEG2RAD)-PIPX(IB)*SIN(VX(IB)*DEG2RAD)
                  PIPY_TMP=-PIPY(IB)*SIN(VX(IB)*DEG2RAD)+PIPX(IB)*COS(VX(IB)*DEG2RAD)       
                  
                  FIJ1=var_ice(IA)+DXA*PIPX(IA)+DYA*PIPY(IA)
                  FIJ2=var_ice(IB)+DXB*PIPX_TMP+DYB*PIPY_TMP
                  
               ELSE IF(IB == NODE_NORTHPOLE)THEN
                  PIPX_TMP=-PIPY(IA)*COS(VX(IA)*DEG2RAD)-PIPX(IA)*SIN(VX(IA)*DEG2RAD)
                  PIPY_TMP=-PIPY(IA)*SIN(VX(IA)*DEG2RAD)+PIPX(IA)*COS(VX(IA)*DEG2RAD)
                  
                  FIJ1=var_ice(IA)+DXA*PIPX_TMP+DYA*PIPY_TMP
                  FIJ2=var_ice(IB)+DXB*PIPX(IB)+DYB*PIPY(IB)
               END IF
               
               VI1MIN=MINVAL(var_ice(NBSN(IA,1:NTSN(IA)-1)))
               VI1MIN=MIN(VI1MIN, var_ice(IA))
               VI1MAX=MAXVAL(var_ice(NBSN(IA,1:NTSN(IA)-1)))
               VI1MAX=MAX(VI1MAX, var_ice(IA))
               VI2MIN=MINVAL(var_ice(NBSN(IB,1:NTSN(IB)-1)))
               VI2MIN=MIN(VI2MIN, var_ice(IB))
               VI2MAX=MAXVAL(var_ice(NBSN(IB,1:NTSN(IB)-1)))
               VI2MAX=MAX(VI2MAX, var_ice(IB))
               IF(FIJ1 < VI1MIN) FIJ1=VI1MIN
               IF(FIJ1 > VI1MAX) FIJ1=VI1MAX
               IF(FIJ2 < VI2MIN) FIJ2=VI2MIN
               IF(FIJ2 > VI2MAX) FIJ2=VI2MAX
               
               !      IF(IA == NODE_NORTHPOLE .OR. IB == NODE_NORTHPOLE)THEN
               UIJ_TMP = -cvice(I1)*COS(XC(I1)*DEG2RAD)-cuice(I1)*SIN(XC(I1)*DEG2RAD)
               VIJ_TMP = -cvice(I1)*SIN(XC(I1)*DEG2RAD)+cuice(I1)*COS(XC(I1)*DEG2RAD)
               
               VX1_TMP = REARTH * COS(YIJE(I,1)*DEG2RAD) * COS(XIJE(I,1)*DEG2RAD)
               VY1_TMP = REARTH * COS(YIJE(I,1)*DEG2RAD) * SIN(XIJE(I,1)*DEG2RAD)
               
               VX2_TMP = REARTH * COS(YIJE(I,2)*DEG2RAD) * COS(XIJE(I,2)*DEG2RAD)
               VY2_TMP = REARTH * COS(YIJE(I,2)*DEG2RAD) * SIN(XIJE(I,2)*DEG2RAD)
               
               DLTXE_TMP = VX2_TMP-VX1_TMP
               DLTYE_TMP = VY2_TMP-VY1_TMP
               
               
               UVN_TMP = VIJ_TMP*DLTXE_TMP - UIJ_TMP*DLTYE_TMP
               EXFLUX_TMP = -UVN_TMP*((1.0_SP+SIGN(1.0_SP,UVN_TMP))*FIJ2+   &
                    (1.0_SP-SIGN(1.0_SP,UVN_TMP))*FIJ1)*0.5_SP
               
               IF(IA == NODE_NORTHPOLE)THEN
                  XFLUX(IA)=XFLUX(IA)+EXFLUX_TMP 
               ELSE IF(IB == NODE_NORTHPOLE)THEN
                  XFLUX(IB)=XFLUX(IB)-EXFLUX_TMP 
               END IF
            END IF
         END IF
      END DO
# if defined (MULTIPROCESSOR)
      !IF(PAR)CALL NODE_MATCH(0,NBN,BN_MLT,BN_LOC,BNC,MT,1,MYID,NPROCS,XFLUX)
# endif
            
      RETURN
    END SUBROUTINE ICE_ADV_XY
    !==============================================================================|
    !#  endif
#  endif
    !end if defined (NORTHPOLE)
    
    !==============================================================================|
    SUBROUTINE ICE_EVP(XFLUX,YFLUX,KI)     
      !------------------------------------------------------------------------------|
      !   USE BCS
      !   USE MOD_OBCS
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: KI
      REAL(SP), ALLOCATABLE :: PUPX(:),PVPX(:),PUPY(:),PVPY(:)
      
      REAL(SP) :: X11,Y11,X22,Y22,X33,Y33,TMP1,TMP2
      REAL(SP) :: ELIJ,UIJ,VIJ,EXFLUX
      INTEGER  :: I,J,K,I1,IA,IB,J1,J2,JTMP
   
      REAL(SP) :: E11,E12,E22
      Real(SP) :: DD_divu,DT_tension,DS_shear
      Real(SP) :: DELT,ZETA,ETA
      Real(SP) :: PDelta
      REAL(SP) :: SIG11IJ,SIG12IJ,SIG22IJ
      REAL(SP) :: XFLUX(0:NT),YFLUX(0:NT)

      !   REAL(SP), ALLOCATABLE :: UICE2N(:), VICE2N(:)
      REAL(SP) :: XTMP,XTMP1
      INTEGER, ALLOCATABLE :: IUVN(:)
      
      REAL(DP) :: DLTXE_TMP,DLTYE_TMP,DLTXC_TMP,DLTYC_TMP
      
      REAL(SP), DIMENSION(0:NT) :: PSIG11PX,PSIG12PY
      REAL(SP), DIMENSION(0:NT) :: PSIG12PX,PSIG22PY
      
      REAL(SP)  :: SIG_TMP
      INTEGER   :: N_TMP
      
      REAL(SP) :: UI,VI
      !==============================================================================|
      ! First calculate the velocity gradient on the nodes
      !----------INITIALIZE FLUX ARRAY ----------------------------------------------!
      
      ALLOCATE(PUPX(0:MT))     ;PUPX= 0.0_SP
      ALLOCATE(PVPX(0:MT))     ;PVPX= 0.0_SP
      ALLOCATE(PUPY(0:MT))     ;PUPY= 0.0_SP
      ALLOCATE(PVPY(0:MT))     ;PVPY= 0.0_SP
      
      Delta= 0.0_SP
      divu = 0.0_SP
      
!!!!================================================================================================
      DO I=1,NCV_I
         I1  = NTRG(I)
         IA  = NIEC(I,1)
         IB  = NIEC(I,2)
         
         IF(ISICEN(IA) == 1 .OR. ISICEN(IB) == 1)THEN
            IF(ISICEC(I1)==1) then
               UIJ = UICE2(I1)
               VIJ = VICE2(I1)
            ELSE
               UIJ = UICE2N(IA)*ISICEN(IA)+UICE2N(IB)*ISICEN(IB)
               VIJ = VICE2N(IA)*ISICEN(IA)+VICE2N(IB)*ISICEN(IB)
            ENDIF
            
            !         UIJ = UICE2(I1)
            !         VIJ = VICE2(I1)
            
            EXFLUX = -UIJ*DLTYE(I)
            PUPX(IA) = PUPX(IA)-EXFLUX
            PUPX(IB) = PUPX(IB)+EXFLUX
            
            EXFLUX = -VIJ*DLTYE(I)
            PVPX(IA) = PVPX(IA)-EXFLUX
            PVPX(IB) = PVPX(IB)+EXFLUX

            EXFLUX =  UIJ*DLTXE(I)
            PUPY(IA) = PUPY(IA)-EXFLUX
            PUPY(IB) = PUPY(IB)+EXFLUX
            
            EXFLUX =  VIJ*DLTXE(I)
            PVPY(IA) = PVPY(IA)-EXFLUX
            PVPY(IB) = PVPY(IB)+EXFLUX
            !!  slip boundary
            !! for the solid boundary edge
            !!     isonb(i)=1:  node on the solid boundary  
            !!   normal velocity == 0

            if(.false.) then
            IF(ISONB(IA)==1 .and. ISONB(IB) == 1) THEN
               UI= UICE2(I1)*(SIN(ALPHA(I1)))**2-VICE2(I1)*SIN(ALPHA(I1))*COS(ALPHA(I1))
               VI=-UICE2(I1)*SIN(ALPHA(I1))*COS(ALPHA(I1))+VICE2(I1)*(COS(ALPHA(I1)))**2
               
#           if defined (SPHERICAL)
               XTMP1 = VX(IB)-VX(IA)
               XTMP = XTMP1*TPI
               IF(XTMP1 > 180.0)THEN
                  XTMP = -360.0*TPI+XTMP
               ELSE IF(XTMP1 < -180.0)THEN
                  XTMP =  360.0*TPI+XTMP
               END IF
               DLTXE_TMP = XTMP*COS(DEG2RAD*YC(I1))
               DLTYE_TMP = (VY(IB)-VY(IA))*TPI
#          else
               DLTXE_TMP = VX(IB)-VX(IA)
               DLTYE_TMP = VY(IB)-VY(IA)
#          endif
               
               EXFLUX =  -UI*DLTYE_TMP*0.5_SP
               PUPX(IA) = PUPX(IA)-EXFLUX
               PUPX(IB) = PUPX(IB)+EXFLUX
               
               EXFLUX = -VI*DLTYE_TMP*0.5_SP
               PVPX(IA) = PVPX(IA)-EXFLUX
               PVPX(IB) = PVPX(IB)+EXFLUX
               
               EXFLUX =  UI*DLTXE_TMP*0.5_SP
               PUPY(IA) = PUPY(IA)-EXFLUX
               PUPY(IB) = PUPY(IB)+EXFLUX
               
               EXFLUX =  VI*DLTXE_TMP*0.5_SP
               PVPY(IA) = PVPY(IA)-EXFLUX
               PVPY(IB) = PVPY(IB)+EXFLUX
            END IF
            !!   end solid boundary---non-slip
           endif  !! slip solid boundary condition
            
         END IF
      END DO
!!!-----------------------------------------------------------------------------
# if defined(MULTIPROCESSOR)
      IF(PAR)CALL NODE_MATCH(0,NBN,BN_MLT,BN_LOC,BNC,MT,1,MYID,NPROCS,PUPX,PVPX)
      IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,PUPX,PVPX)
      IF(PAR)CALL NODE_MATCH(0,NBN,BN_MLT,BN_LOC,BNC,MT,1,MYID,NPROCS,PUPY,PVPY)
      IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,PUPY,PVPY)
# endif

      PUPX = PUPX/ART1
      PVPX = PVPX/ART1
      PUPY = PUPY/ART1
      PVPY = PVPY/ART1
      
      DO I=1,M
         IF(ISICEN(I) == 1)THEN
            E11 = PUPX(I)
            E12 = 0.5_SP*(PUPY(I)+PVPX(I))
            E22 = PVPY(I)
            
            ! divergence  =  e_11 + e_22
            DD_divu = E11+E22
            ! tension strain rate  =  e_11 - e_22
            DT_tension = E11-E22
            ! shearing strain rate  =  e_12
            DS_shear = 2.0_SP*E12  
            
            DELT = DD_divu*DD_divu+ECCI*(DT_tension*DT_tension+DS_shear*DS_shear)
            DELT = SQRT(DELT)
            !
            !--------------------------------------------------  
            !       calculate stress for ridging
            DELT= max(DELT,2.E-9)
            IF (KI == ndte) THEN 
               DeltaN (I)= DELT 
               DivuN  (I)= DD_divu 
            ENDIF

            !--------------------------------------------------  
            PDelta = PICE(I)*1.0E-6/DELT
            
            !!====================================================================================|
            !     first step stress equations
            !     computing sigma1 and sigma2
            !     SIG1 (I) = ( SIG1 (I) + dte2T*(PDelta*DD_divu-Pice(I))*1.0E+6)  * denom1
            SIG1 (I) = ( SIG1 (I) + dte2T*(PDelta*(DD_divu-DELT))*1.0E+6)*denom1
            SIG2 (I) = ( SIG2 (I) + dte2T*PDelta*DT_tension      *1.0E+6)*denom2
            SIG12(I) = ( SIG12(I) + dte2T*Pdelta*0.5_SP*DS_shear *1.0E+6)*denom2
            
            !     recover sigma11 and sigma22
            SIG11(I)=(SIG1(I)+SIG2(I))*0.5_SP
            SIG22(I)=(SIG1(I)-SIG2(I))*0.5_SP
            
            !!===================================================================================|
            prs_sig(I)=PDelta*DELT*1.E+6
            
         ENDIF       
      END DO
      
# if defined(MULTIPROCESSOR)
      IF(PAR)CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,1,MYID,NPROCS,SIG1,SIG2)
      IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,SIG1,SIG2)
      IF(PAR)CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,1,MYID,NPROCS,SIG11,SIG12,SIG22)
      IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,SIG11,SIG12,SIG22)
# endif
      
      PSIG11PX = 0.0_SP
      PSIG12PY = 0.0_SP
      PSIG12PX = 0.0_SP
      PSIG22PY = 0.0_SP 
      
      DO I=1,NE
         IA=IEC(I,1)
         IB=IEC(I,2)
         J1=IENODE(I,1)
         J2=IENODE(I,2)
         
         IF(ISICEC(IA) == 1 .OR. ISICEC(IB) == 1)THEN
            ELIJ=0.5_SP*(EL(J1)+EL(J2))
            ! no surface tilt contribution ggao 0516 
            ELIJ=0.0    
            
            SIG11IJ=0.5_SP*(SIG11(J1)+SIG11(J2))
            SIG22IJ=0.5_SP*(SIG22(J1)+SIG22(J2))
            SIG12IJ=0.5_SP*(SIG12(J1)+SIG12(J2))
            
#      if defined (SPHERICAL)
            !for spherical coordinator and domain across 360^o longitude         
            XTMP  = VX(J2)*TPI-VX(J1)*TPI
            XTMP1 = VX(J2)-VX(J1)
            IF(XTMP1 >  180.0_SP)THEN
               XTMP = -360.0_SP*TPI+XTMP
            ELSE IF(XTMP1 < -180.0_SP)THEN
               XTMP =  360.0_SP*TPI+XTMP
            END IF
            !-------------------------------------------------------------
            EXFLUX    = SIG11IJ*DLTYC(I)
            PSIG11PX(IA) = PSIG11PX(IA) - EXFLUX/ART(IA)
            PSIG11PX(IB) = PSIG11PX(IB) + EXFLUX/ART(IB)
            
            EXFLUX    = SIG12IJ*XTMP*COS(DEG2RAD*YC(IA))
            PSIG12PY(IA) = PSIG12PY(IA) + EXFLUX/ART(IA)
            EXFLUX    = SIG12IJ*XTMP*COS(DEG2RAD*YC(IB))
            PSIG12PY(IB) = PSIG12PY(IB) - EXFLUX/ART(IB)
            
            EXFLUX    = SIG12IJ*DLTYC(I)
            PSIG12PX(IA) = PSIG12PX(IA) - EXFLUX/ART(IA)
            PSIG12PX(IB) = PSIG12PX(IB) + EXFLUX/ART(IB)
            
            EXFLUX    = SIG22IJ*XTMP*COS(DEG2RAD*YC(IA))
            PSIG22PY(IA) = PSIG22PY(IA) + EXFLUX/ART(IA)
            EXFLUX    = SIG22IJ*XTMP*COS(DEG2RAD*YC(IB))       
            PSIG22PY(IB) = PSIG22PY(IB) - EXFLUX/ART(IB)
            
            !      ACCUMULATE BAROTROPIC FLUX
            PSTXICE(IA)=PSTXICE(IA)-GRAV_E(IA)*ELIJ*DLTYC(I)/ART(IA)
            PSTYICE(IA)=PSTYICE(IA)+GRAV_E(IA)*ELIJ*XTMP*COS(DEG2RAD*YC(IA))/ART(IA)
            
            PSTXICE(IB)=PSTXICE(IB)+GRAV_E(IB)*ELIJ*DLTYC(I)/ART(IB)
            PSTYICE(IB)=PSTYICE(IB)-GRAV_E(IB)*ELIJ*XTMP*COS(DEG2RAD*YC(IB))/ART(IB)
            
#      else
            EXFLUX    = SIG11IJ*DLTYC(I)
            PSIG11PX(IA) = PSIG11PX(IA) - EXFLUX/ART(IA)
            PSIG11PX(IB) = PSIG11PX(IB) + EXFLUX/ART(IB)
            
            EXFLUX    = SIG12IJ*DLTXC(I)
            PSIG12PY(IA) = PSIG12PY(IA) + EXFLUX/ART(IA)
            PSIG12PY(IB) = PSIG12PY(IB) - EXFLUX/ART(IB)
            
            EXFLUX    = SIG12IJ*DLTYC(I)
            PSIG12PX(IA) = PSIG12PX(IA) - EXFLUX/ART(IA)
            PSIG12PX(IB) = PSIG12PX(IB) + EXFLUX/ART(IB)
            
            EXFLUX    = SIG22IJ*DLTXC(I)
            PSIG22PY(IA) = PSIG22PY(IA) + EXFLUX/ART(IA)
            PSIG22PY(IB) = PSIG22PY(IB) - EXFLUX/ART(IB)
            
            !      ACCUMULATE BAROTROPIC FLUX
            PSTXICE(IA)=PSTXICE(IA)-GRAV_E(IA)*ELIJ*DLTYC(I)/ART(IA)
            PSTYICE(IA)=PSTYICE(IA)+GRAV_E(IA)*ELIJ*DLTXC(I)/ART(IA)
            
            PSTXICE(IB)=PSTXICE(IB)+GRAV_E(IB)*ELIJ*DLTYC(I)/ART(IB)
            PSTYICE(IB)=PSTYICE(IB)-GRAV_E(IB)*ELIJ*DLTXC(I)/ART(IB)
            
#    endif     
         END IF
      END DO
   
#  if defined (SPHERICAL)
      CALL ICE_EVP_XY(PSIG11PX,PSIG12PY,PSIG12PX,PSIG22PY,KI)
#  endif
      
      XFLUX=0.0_SP
      YFLUX=0.0_SP

      DO I=1,N
         IF(sum(ISICEN(NV(I,1:3)))==3) THEN
            XFLUX(I) =  PSIG11PX(I)+PSIG12PY(I)
            YFLUX(I) =  PSIG12PX(I)+PSIG22PY(I)
         END IF
      END DO
      
      DEALLOCATE(PUPX,PUPY,PVPX,PVPY)
      
      RETURN
    END SUBROUTINE ICE_EVP
    !==============================================================================|
    !=======================================================================
    !  principal_stress - computes principal stress for yield curve
    !
    !
    subroutine principal_stress
      !
      ! Computes principal stresses for comparison with the theoretical
      ! yield curve; northeast values
      !
      integer (kind=int_kind) :: i
      do i=1,M
         !if(prs_sig(i) > puny) then
         if(prs_sig(i) > 1.E-8) then
            Psig1(i)=(p5*(sig1(i)+sqrt(sig2(i)**2+c4I*sig12(i)**2)))/prs_sig(i)
            Psig2(i)=(p5*(sig1(i)-sqrt(sig2(i)**2+c4I*sig12(i)**2)))/prs_sig(i)
         else
            Psig1(i)=1.0E+2
            Psig2(i)=1.0E+2
         endif
      enddo
      
    end subroutine principal_stress
        




#  if defined (ICE_EMBEDDING)
!==============================================================================|
    ! subroutine to find ice-water boundary
    ! yding
    !
    ! J. Ge modified for parallel running and 
    !    ice baroclinic pressure gradient
    !
    subroutine FIND_ICE_WATER_BOUNDARY

    implicit none
    !
    INTEGER :: I,K,IOUTTMP,ierr 
    !
    !

    do i = 1, N
      if(aice(nv(i,1),1)>1.0e-5 .and. aice(nv(i,2),1)<1.0e-5 .and. aice(nv(i,3),1)<1.0e-5) then
        iwbnd(i) = 1
      else if (aice(nv(i,2),1)>1.0e-5 .and. aice(nv(i,1),1)<1.0e-5 .and. aice(nv(i,3),1)<1.0e-5) then
        iwbnd(i) = 1
      else if (aice(nv(i,3),1)>1.0e-5 .and. aice(nv(i,1),1)<1.0e-5 .and. aice(nv(i,2),1)<1.0e-5) then
        iwbnd(i) = 1
      else if (aice(nv(i,1),1)>1.0e-5 .and. aice(nv(i,2),1)>1.0e-5 .and. aice(nv(i,3),1)<1.0e-5) then
        iwbnd(i) = 1
      else if (aice(nv(i,1),1)>1.0e-5 .and. aice(nv(i,3),1)>1.0e-5 .and. aice(nv(i,2),1)<1.0e-5) then
        iwbnd(i) = 1
      else if (aice(nv(i,2),1)>1.0e-5 .and. aice(nv(i,3),1)>1.0e-5 .and. aice(nv(i,1),1)<1.0e-5) then
        iwbnd(i) = 1
      else
        iwbnd(i) = 0
      end if
    end do 
   end subroutine FIND_ICE_WATER_BOUNDARY 





!=============================================================================|
!   calculate the bottom sigma index of embedding icepack
!=============================================================================|
   SUBROUTINE GET_Kzice
   
   IMPLICIT NONE
   INTEGER  :: I,I1,I2,I3,I4,J,K
   REAL(SP) :: DTMP,ice_depth

   if(dbg_set(dbg_sbr)) write(ipt,*) "Start: get_kzice "    

   DO I=1,N

     if(iwbnd(i)==0)cycle

     D1(I) = H1(I)+ELF1(I)
     ice_depth = max(QTHICE(NV(I,1)),QTHICE(NV(I,2)),QTHICE(NV(I,3)))
     Kzice(I) = 0
     DTMP = 0.0_SP

     DO K=1,KB
       DTMP = DTMP + DZ1(I,K)*D1(I)
       IF(DTMP <= ice_depth)THEN
         Kzice(I) = Kzice(I) + 1
       END IF
     END DO  	 
     
   END DO

   if(dbg_set(dbg_sbr)) write(ipt,*) "End: get_kzice " 

   Return

 End SUBROUTINE GET_Kzice


!==============================================================================|
!     CALCULATE THE BAROCLINIC PRESSURE GRADIENT IN SIGMA COORDINATES          |
!==============================================================================|

   SUBROUTINE BAROPG_ICE

!==============================================================================|
   USE ALL_VARS
   USE MOD_SPHERICAL
   USE MOD_NORTHPOLE
   USE MOD_WD

#  if defined (MULTIPROCESSOR)
   USE MOD_PAR
#  endif

   IMPLICIT NONE
   REAL(SP) :: RIJK(0:N,3,KBM1), DRIJK1(0:N,3,KBM1), DRIJK2(0:N,KBM1)
   REAL(SP) :: TEMP,DIJ,DRHO1,DRHO2
   INTEGER  :: I,K,J,J1,J2,IJK
#  if defined (SPHERICAL)
   REAL(SP) :: XTMP,XTMP1
#  endif

   integer :: ierr
   integer,allocatable, dimension(:)   :: iwbndtmp
   real(sp),allocatable,dimension(:)   :: art1tmp
   real(sp),allocatable,dimension(:,:) :: aicetmp,vicetmp,drhoxtmp,drhoytmp
   real(dp)::total_mass,total_drhox,total_drhoy


!==============================================================================|

   IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "Start: baropg.F"
   
   ! USE RAMP CALCULATED IN 'internal_step.F'

!----------SUBTRACT MEAN DENSITY TO MINIMIZE ROUNDOFF ERROR--------------------!

   RHO1(:,1:KBM1) = RHO1(:,1:KBM1) - RMEAN1(:,1:KBM1)
   RHO = RHO - RMEAN 

!----------INITIALIZE ARRAYS---------------------------------------------------!

   DRHOX_ice  = 0.0_SP
   DRHOY_ice  = 0.0_SP
   RMEAN(0,:) = 0.0_SP
   RHO(0,:)   = 0.0_SP
   RIJK       = 0.0_SP
   DRIJK1     = 0.0_SP
   DRIJK2     = 0.0_SP

!----------CALCULATE AVERAGE DENSITY ON EACH EDGE------------------------------!

   DO K=1,KBM1
     DO I=1,N
!------J. Ge-------
       if(iwbnd(i)==0)cycle
       if(k>kzice(i))cycle
!------J. Ge-------
       DO J=1,3
         J1=J+1-INT((J+1)/4)*3
         J2=J+2-INT((J+2)/4)*3
         RIJK(I,J,K)  = 0.5_SP*(RHO1(NV(I,J1),K)+RHO1(NV(I,J2),K))
       END DO
     END DO
   END DO

   DO I=1,N
!------J. Ge-------
     if(iwbnd(i)==0)cycle
!------J. Ge-------
     DO J=1,3
       DRIJK1(I,J,1)=RIJK(I,J,1)*(-ZZ1(I,1))
       DO K=2,KBM1
!------J. Ge-------
         if(k>kzice(i))cycle
!------J. Ge-------
         DRIJK1(I,J,K)=0.5_SP*(RIJK(I,J,K-1)+RIJK(I,J,K))*(ZZ1(I,K-1)-ZZ1(I,K))
         DRIJK1(I,J,K)=DRIJK1(I,J,K)+DRIJK1(I,J,K-1)
       END DO
     END DO
   END DO

   DO I=1,N
!------J. Ge-------
     if(iwbnd(i)==0)cycle
!------J. Ge-------
     DRIJK2(I,1)=0.0_SP
     DO K=2,KBM1
!------J. Ge-------
         if(k>kzice(i))cycle
!------J. Ge-------
       DRIJK2(I,K)=0.5_SP*(ZZ1(I,K-1)+ZZ1(I,K))*(RHO(I,K)-RHO(I,K-1)) 
       DRIJK2(I,K)=DRIJK2(I,K-1)+DRIJK2(I,K)
     END DO
   END DO

   DO I = 1, N
!------J. Ge-------
     if(iwbnd(i)==0)cycle
!------J. Ge-------
#  if defined (WET_DRY)
    IF(ISWETCT(I)*ISWETC(I) == 1 .AND. &
      (H(NV(I,1)) > STATIC_SSH_ADJ .OR. H(NV(I,2)) > STATIC_SSH_ADJ .OR. H(NV(I,3)) > STATIC_SSH_ADJ))THEN
#  endif
     DO K=1,KBM1
!------J. Ge-------
        if(k>kzice(i))cycle
!------J. Ge-------
        DO J = 1, 3
          J1=J+1-INT((J+1)/4)*3
          J2=J+2-INT((J+2)/4)*3
          IJK=NBE(I,J)
          DIJ=0.5_SP*(DT(NV(I,J1))+DT(NV(I,J2)))

#         if defined (SPHERICAL)
          DRHO1=-DELTUY(I,J)*DRIJK1(I,J,K)*DT1(I)
          DRHO2=-DELTUY(I,J)*DIJ*DRIJK2(I,K)
#         else
          DRHO1=(VY(NV(I,J1))-VY(NV(I,J2)))*DRIJK1(I,J,K)*DT1(I)
          DRHO2=(VY(NV(I,J1))-VY(NV(I,J2)))*DIJ*DRIJK2(I,K)
#         endif
          DRHOX_ice(I,K)=DRHOX_ice(I,K)+DRHO1+DRHO2

#         if defined (SPHERICAL)
          XTMP  = VX(NV(I,J2))*TPI-VX(NV(I,J1))*TPI
          XTMP1 = VX(NV(I,J2))-VX(NV(I,J1))
          IF(XTMP1 >  180.0_SP)THEN
            XTMP = -360.0_SP*TPI+XTMP
          ELSE IF(XTMP1 < -180.0_SP)THEN
            XTMP =  360.0_SP*TPI+XTMP
          END IF  

          DRHO1=XTMP*COS(DEG2RAD*YC(I))*DRIJK1(I,J,K)*DT1(I)
          DRHO2=XTMP*COS(DEG2RAD*YC(I))*DIJ*DRIJK2(I,K)
!          DRHO1=DELTUX(I,J)*DRIJK1(I,J,K)*DT1(I)
!          DRHO2=DELTUX(I,J)*DIJ*DRIJK2(I,K)
#         else
	  DRHO1=(VX(NV(I,J2))-VX(NV(I,J1)))*DRIJK1(I,J,K)*DT1(I)
          DRHO2=(VX(NV(I,J2))-VX(NV(I,J1)))*DIJ*DRIJK2(I,K)
#         endif
          DRHOY_ice(I,K)=DRHOY_ice(I,K)+DRHO1+DRHO2

       END DO
     END DO
#  if defined (WET_DRY)
    END IF
#  endif
   END DO

#  if defined (SPHERICAL)
   CALL BAROPG_XY(DRIJK1,DRIJK2)
#  endif  


!----------MULTIPLY BY GRAVITY AND ELEMENT DEPTH-------------------------------!

   DO K=1,KBM1
     DRHOX_ice(:,K)=DRHOX_ice(:,K)*DT1(:)*DZ1(:,K)*GRAV_E(:)*RAMP
     DRHOY_ice(:,K)=DRHOY_ice(:,K)*DT1(:)*DZ1(:,K)*GRAV_E(:)*RAMP
   END DO

!----------ADD MEAN DENSITY BACK ON--------------------------------------------!

   RHO1 = RHO1 + RMEAN1
   RHO  = RHO  + RMEAN





   IF(SERIAL)THEN
!----------calculate the whole ice mass----------------------------------------!
    total_mass = 0.0
    do i=1,m
      total_mass = total_mass + art1(i)*vice(i,1)*aice(i,1)*rhoi
    end do

!----------accumulate all gradient force---------------------------------------!
    total_drhox = 0.0
    total_drhoy = 0.0
    do i=1,n
       if(iwbnd(i)==0)cycle
       total_drhox = total_drhox + sum(DRHOX_ice (i,:))
       total_drhoy = total_drhoy + sum(DRHOy_ice (i,:))
    end do

   ENDIF


#  if defined (MULTIPROCESSOR) 
   IF(PAR)THEN

     CALL MPI_BARRIER(MPI_FVCOM_GROUP,IERR)

     ALLOCATE(aicetmp(0:mGL,1))  
     allocate(art1tmp(0:mgl)) 
     allocate(vicetmp(0:mgl,1)) 
     allocate(drhoxtmp(0:ngl,kb))
     allocate(drhoytmp(0:ngl,kb))
     allocate(iwbndtmp(0:ngl)) 

     CALL ACOLLECT(MYID,MSRID,NPROCS,NMAP,art1,art1tmp)
     CALL ACOLLECT(MYID,MSRID,NPROCS,NMAP,vice,vicetmp)
     CALL ACOLLECT(MYID,MSRID,NPROCS,NMAP,aice,aicetmp)
     CALL ACOLLECT(MYID,MSRID,NPROCS,EMAP,DRHOX_ice,drhoxtmp)
     CALL ACOLLECT(MYID,MSRID,NPROCS,EMAP,DRHOY_ice,drhoytmp)
     CALL ACOLLECT(MYID,MSRID,NPROCS,EMAP,iwbnd,iwbndtmp)

!----------calculate the whole ice mass----------------------------------------!
     total_mass = 0.0
     do i=1,mgl
       total_mass = total_mass + art1tmp(i)*vicetmp(i,1)*aicetmp(i,1)*rhoi
     end do

!----------accumulate all gradient force---------------------------------------!
     total_drhox = 0.0
     total_drhoy = 0.0
     do i=1,ngl
       if(iwbndtmp(i)==0)cycle
       total_drhox = total_drhox + sum(DRHOXtmp(i,:))
       total_drhoy = total_drhoy + sum(DRHOytmp(i,:))
     end do     

   END IF

#  endif

!----------put the total gradient force on whole icepack-----------------------!
   net_DRHOX = total_drhox / total_mass
   net_DRHOY = total_drhoy / total_mass


   IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "End: baropg.F"

   RETURN
   END SUBROUTINE BAROPG_ICE
!==============================================================================|





!==============================================================================|

#  endif


#  endif
  END MODULE MOD_ICE2D
 
