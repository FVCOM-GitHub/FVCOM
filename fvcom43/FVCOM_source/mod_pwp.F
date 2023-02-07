  MODULE MOD_PWP
# if defined (PWP)
    USE MOD_PREC, ONLY : SP, DP

    IMPLICIT NONE

    PRIVATE
    PUBLIC          :: SETUP_PWP
    PUBLIC          :: PWPGO 
    PUBLIC          :: NAME_LIST_INITIALIZE_PWP
    PUBLIC          :: NAME_LIST_PRINT_PWP
    PUBLIC          :: NAME_LIST_READ_PWP

    REAL(SP),PUBLIC :: UPPER_DEPTH_LIMIT
    REAL(SP),PUBLIC :: LOWER_DEPTH_LIMIT
    REAL(SP),PUBLIC :: VERTICAL_RESOLUTION
    REAL(SP),PUBLIC :: BULK_RICHARDSON
    REAL(SP),PUBLIC :: GRADIENT_RICHARDSON
    INTEGER,PUBLIC,ALLOCATABLE :: MLD(:)   

    INTEGER  :: N_ORDER
    INTEGER  :: K_MAX
    INTEGER  :: ML
    REAL(SP) :: D_UPPER
    REAL(SP) :: D_LOWER
    REAL(SP) :: DDZ
    REAL(SP) :: DDDZ
    REAL(SP) :: TR
    REAL(SP) :: SR
    REAL(SP) :: DR
    REAL(SP) :: D_ALPHA
    REAL(SP) :: D_BETA

    INTEGER , ALLOCATABLE   ::  K_IDX(:)
    INTEGER , ALLOCATABLE   ::  NZ(:)
    REAL(SP), ALLOCATABLE   ::  TTT(:)
    REAL(SP), ALLOCATABLE   ::  SSS(:)
    REAL(SP), ALLOCATABLE   ::  UUU(:)
    REAL(SP), ALLOCATABLE   ::  VVV(:)
    REAL(SP), ALLOCATABLE   ::  DENSITY(:)

    REAL(SP), ALLOCATABLE   ::  AVG_U(:,:)
    REAL(SP), ALLOCATABLE   ::  AVG_V(:,:)

!    logical check

    NAMELIST /NML_PWP/                 &
            & UPPER_DEPTH_LIMIT,       &
            & LOWER_DEPTH_LIMIT,       &
            & VERTICAL_RESOLUTION,     &
            & BULK_RICHARDSON,         &
            & GRADIENT_RICHARDSON

    CONTAINS

!-----------------------------------------------------
!
!-----------------------------------------------------
    SUBROUTINE NAME_LIST_INITIALIZE_PWP

    IMPLICIT NONE

    UPPER_DEPTH_LIMIT   = 20.0_SP
    LOWER_DEPTH_LIMIT   = 200.0_SP
    VERTICAL_RESOLUTION = 1.0_SP 
    BULK_RICHARDSON     = 0.65_SP
    GRADIENT_RICHARDSON = 0.25_SP

    END SUBROUTINE NAME_LIST_INITIALIZE_PWP

!-----------------------------------------------------
!
!-----------------------------------------------------
    SUBROUTINE NAME_LIST_PRINT_PWP
    USE CONTROL, ONLY : IPT

    IMPLICIT NONE

    WRITE (UNIT=IPT,NML=NML_PWP)

    RETURN
    END SUBROUTINE NAME_LIST_PRINT_PWP

!-----------------------------------------------------
!
!-----------------------------------------------------
    SUBROUTINE NAME_LIST_READ_PWP

    USE CONTROL
    USE MOD_UTILS
    IMPLICIT NONE
    INTEGER  :: ISCAN
    CHARACTER(LEN=120) :: FNAME
    INTEGER :: ios

    ios = 0

    FNAME = "./"//trim(casename)//"_run.nml"

    if(DBG_SET(dbg_io)) &
           & write(IPT,*) "Read_Name_List: File: ",trim(FNAME)

    CALL FOPEN(NMLUNIT,trim(FNAME),'cfr')

    ! Read NH Settings
    READ(UNIT=NMLUNIT, NML=NML_PWP,IOSTAT=ios)
    if(ios .NE. 0 ) then
      if(DBG_SET(dbg_log)) write(UNIT=IPT,NML=NML_PWP)
      Call Fatal_Error("Can Not Read NameList NML_PWP from file: "//trim(FNAME))
    endif

    REWIND(NMLUNIT)
    CLOSE(NMLUNIT)

    D_UPPER = UPPER_DEPTH_LIMIT
    D_LOWER = LOWER_DEPTH_LIMIT
    DDZ     = VERTICAL_RESOLUTION

    RETURN
    END SUBROUTINE NAME_LIST_READ_PWP

!-----------------------------------------------------
!
!-----------------------------------------------------
    SUBROUTINE SETUP_PWP

    USE MOD_PREC, ONLY : SP, DP
    USE ALL_VARS, ONLY : GRAV, M, N, MT, NT, MGL, KBM1
    IMPLICIT NONE 

    K_MAX = ANINT(D_LOWER/DDZ) + 1

    ALLOCATE(K_IDX(0:MT))             ; K_IDX  = -99
    ALLOCATE(NZ(0:MT))                ; NZ     = -99
    ALLOCATE(MLD(0:MT))               ; MLD    = -99

    ALLOCATE(AVG_U(0:MT,KBM1))        ; AVG_U  = 0.0_SP
    ALLOCATE(AVG_V(0:MT,KBM1))        ; AVG_V  = 0.0_SP

    CALL FIND_NZ

    END SUBROUTINE SETUP_PWP

!-----------------------------------------------------
!
!-----------------------------------------------------
    SUBROUTINE FIND_NZ

    USE ALL_VARS
    IMPLICIT NONE

    INTEGER :: I, K, KK

    DO I=1, M
      IF(H(I)>D_UPPER) THEN
        DO K=1, KB
          IF(-H(I)*Z(I,K)>D_LOWER) THEN
            K_IDX(I) = K
            EXIT
          ENDIF
          K_IDX(I) = K
        ENDDO

        IF( (-H(I)*Z(I,K_IDX(I)))>D_LOWER ) THEN
          NZ(I) = K_MAX
        ELSE
          NZ(I) = ANINT(-H(I)*Z(I,K_IDX(I))/DDZ) + 1
        ENDIF

      ENDIF
    ENDDO

    END SUBROUTINE FIND_NZ

!-----------------------------------------------------
!
!-----------------------------------------------------
    SUBROUTINE PWPGO

    USE ALL_VARS
#   if defined(MULTIPROCESSOR)
    USE MOD_PAR, only: EC, AEXCHANGE
#   endif
!    use mod_par, only: ngid
!    use mod_utils
    IMPLICIT NONE

    INTEGER :: I, J, K
!    integer :: iunit

!    real(sp) :: rv
!    real(sp), allocatable :: rc(:)

#   if defined(MULTIPROCESSOR)
    IF(PAR) CALL AEXCHANGE(EC,MYID,NPROCS,UF,VF)
#   endif
    AVG_U  = 0.0_SP
    AVG_V  = 0.0_SP
    DO K=1, KBM1
      DO I=1, M
        DO J=1, NTVE(I)
          AVG_U(I,K) = AVG_U(I,K) + UF(NBVE(I,J),K)*ART(NBVE(I,J))
          AVG_V(I,K) = AVG_V(I,K) + VF(NBVE(I,J),K)*ART(NBVE(I,J))
        ENDDO
        AVG_U(I,K) = AVG_U(I,K)/ART2(I)
        AVG_V(I,K) = AVG_V(I,K)/ART2(I)
      ENDDO
    ENDDO

!    if(iint==210386_DP) call pshutdown

    MLD = -99
    DO N_ORDER = 1, M 

!      if(ngid(n_order)==9810) then
!        iunit = 800
!        check = .true.
!      elseif(ngid(n_order)==9801) then
!        iunit = 810
!        check = .true.
!      elseif(ngid(n_order)==9825) then
!        iunit = 820
!        check = .true.
!      elseif(ngid(n_order)==12004) then
!        iunit = 830
!        check = .true.
!      elseif(ngid(n_order)==7748) then
!        iunit = 840
!        check = .true.
!      else
!        check = .false.
!      endif

      IF(NZ(N_ORDER)>-90) THEN

        CALL GEN_PROFILE

!        if(check) then
!          write(iunit+1,'(i10,500e20.6)') iint, (-float(k-1)*dddz,k=1,nz(n_order))
!          write(iunit+1,'(i10,500e20.6)') iint, (uuu(k),k=1,nz(n_order))
!          write(iunit+1,'(i10,500e20.6)') iint, (vvv(k),k=1,nz(n_order))
!          write(iunit+1,'(i10,500e20.6)') iint, (ttt(k),k=1,nz(n_order))
!          write(iunit+1,'(i10,500e20.6)') iint, (sss(k),k=1,nz(n_order))
!          write(iunit+1,'(i10,500e20.6)') iint, (density(k),k=1,nz(n_order))
!        endif

        CALL MLDEP

!        if(check) then
!          write(iunit+2,'(i10,6e20.6)') iint, FLOAT(ML)*DDDZ 
!        endif

        CALL FREE_CONVECTION
        CALL MLDEP

!        if(check) then
!          write(iunit+3,'(i10,6e20.6)') iint, FLOAT(ML)*DDDZ
!        endif

        CALL BULK_MIXING  !(check,iunit)

!        allocate(rc(nz(n_order)))
        CALL GRADIENT_MIXING    !(rc)
        CALL MLDEP

!        if(check) then
!          write(iunit+5,'(i10,500e20.6)') iint, FLOAT(ML)*DDDZ, (rc(k),k=1,nz(n_order))
!        endif
!        deallocate(rc) 

        DO K=1, KBM1
          IF( -D(N_ORDER)*ZZ(N_ORDER,K)>FLOAT(ML)*DDDZ ) EXIT
          MLD(N_ORDER) = K
        ENDDO

!        if(check) then
!          write(iunit+6,'(2i10)') iint, MLD(N_ORDER)
!        endif

        DEALLOCATE(TTT,SSS,UUU,VVV,DENSITY)
      ENDIF

    ENDDO


    END SUBROUTINE PWPGO

!-----------------------------------------------------
!
!-----------------------------------------------------
    SUBROUTINE GEN_PROFILE

    USE ALL_VARS
    IMPLICIT NONE

    INTEGER :: I, J, K
    REAL(SP)  ::  Z_, GG

    I = N_ORDER

    ALLOCATE(UUU(NZ(I)))       ; UUU     = -999.0_SP
    ALLOCATE(VVV(NZ(I)))       ; VVV     = -999.0_SP
    ALLOCATE(TTT(NZ(I)))       ; TTT     = -999.0_SP
    ALLOCATE(SSS(NZ(I)))       ; SSS     = -999.0_SP
    ALLOCATE(DENSITY(NZ(I)))   ; DENSITY = -999.0_SP

    DDDZ = -D(I)*Z(I,K_IDX(I))/REAL(NZ(I)-1)
    DO J=1, NZ(I)
      Z_ = (J-1)*DDDZ
      IF(Z_<=-D(I)*ZZ(I,1)) THEN
        UUU(J) = AVG_U(I,1)
        VVV(J) = AVG_V(I,1)
        TTT(J) = T1(I,1)
        SSS(J) = S1(I,1)
      ELSE IF(Z_>=-D(I)*ZZ(I,KBM1)) THEN
        UUU(J) = AVG_U(I,KBM1)
        VVV(J) = AVG_V(I,KBM1)
        TTT(J) = T1(I,KBM1)
        SSS(J) = S1(I,KBM1)
      ELSE
        IF( K_IDX(I)>KBM1 ) THEN
          DO K=2,K_IDX(I)-1
            IF( Z_>-D(I)*ZZ(I,K-1) .AND. Z_<=-D(I)*ZZ(I,K) ) THEN
              UUU(J) = ( AVG_U(I,K-1)*(-D(I)*ZZ(I,K)-Z_) + AVG_U(I,K)*(Z_+D(I)*ZZ(I,K-1)) )/(-D(I)*ZZ(I,K)+D(I)*ZZ(I,K-1))
              VVV(J) = ( AVG_V(I,K-1)*(-D(I)*ZZ(I,K)-Z_) + AVG_V(I,K)*(Z_+D(I)*ZZ(I,K-1)) )/(-D(I)*ZZ(I,K)+D(I)*ZZ(I,K-1))
              TTT(J) = ( T1(I,K-1)*(-D(I)*ZZ(I,K)-Z_) + T1(I,K)*(Z_+D(I)*ZZ(I,K-1)) )/(-D(I)*ZZ(I,K)+D(I)*ZZ(I,K-1))
              SSS(J) = ( S1(I,K-1)*(-D(I)*ZZ(I,K)-Z_) + S1(I,K)*(Z_+D(I)*ZZ(I,K-1)) )/(-D(I)*ZZ(I,K)+D(I)*ZZ(I,K-1))
              EXIT          
            ENDIF
          ENDDO
        ELSE
          DO K=2,K_IDX(I)
            IF( Z_>-D(I)*ZZ(I,K-1) .AND. Z_<=-D(I)*ZZ(I,K) ) THEN
              UUU(J) = ( AVG_U(I,K-1)*(-D(I)*ZZ(I,K)-Z_) + AVG_U(I,K)*(Z_+D(I)*ZZ(I,K-1)) )/(-D(I)*ZZ(I,K)+D(I)*ZZ(I,K-1))
              VVV(J) = ( AVG_V(I,K-1)*(-D(I)*ZZ(I,K)-Z_) + AVG_V(I,K)*(Z_+D(I)*ZZ(I,K-1)) )/(-D(I)*ZZ(I,K)+D(I)*ZZ(I,K-1))
              TTT(J) = ( T1(I,K-1)*(-D(I)*ZZ(I,K)-Z_) + T1(I,K)*(Z_+D(I)*ZZ(I,K-1)) )/(-D(I)*ZZ(I,K)+D(I)*ZZ(I,K-1))
              SSS(J) = ( S1(I,K-1)*(-D(I)*ZZ(I,K)-Z_) + S1(I,K)*(Z_+D(I)*ZZ(I,K-1)) )/(-D(I)*ZZ(I,K)+D(I)*ZZ(I,K-1))
              EXIT
            ENDIF
          ENDDO
        ENDIF
      ENDIF
    ENDDO

    CALL CAL_DENSITY

!    TR = TTT(1)
!    SR = SSS(1)
!    DR = SGT(TR,SR,GG)
!    D_ALPHA = SGT(TR+0.5_SP,SR,GG) - SGT(TR-0.5_SP,SR,GG)
!    D_BETA  = SGT(TR,SR+0.5_SP,GG) - SGT(TR,SR-0.5_SP,GG)

!    DO J=1, NZ(I)
!      DENSITY(J) = DR + (TTT(J)-TR)*D_ALPHA + (SSS(J)-SR)*D_BETA 
!    ENDDO

    END SUBROUTINE GEN_PROFILE

!-----------------------------------------------------
!
!-----------------------------------------------------
    SUBROUTINE CAL_DENSITY

    USE ALL_VARS
    USE EQS_OF_STATE
    IMPLICIT NONE

    INTEGER :: I, K
    REAL(SP), PARAMETER ::PR = 0.0_SP
    REAL(SP), ALLOCATABLE, DIMENSION(:) :: RZU

    I = N_ORDER 

    SELECT CASE(SEA_WATER_DENSITY_FUNCTION)
    CASE(SW_DENS1)

      ALLOCATE(RZU(NZ(I)))       ; RZU = -999.0_SP
      DO K=1,NZ(I)
        RZU(K) = GRAV_N(I)*1.025_SP*((K-1)*DDDZ)*0.1_SP
      END DO
      CALL FOFONOFF_MILLARD(SSS,TTT,RZU,PR,DENSITY)
      DEALLOCATE(RZU)

    CASE(SW_DENS2)

      CALL DENS2G(SSS,TTT,DENSITY)

    CASE(SW_DENS3)

       ALLOCATE(RZU(NZ(I)))       ; RZU = -999.0_SP
       DO K=1,NZ(I)
         RZU(K) = GRAV_N(I)*1.025_SP*((K-1)*DDDZ)*0.01_SP
       END DO
       CALL JACKET_MCDOUGALL(SSS,TTT,RZU,DENSITY)
       DEALLOCATE(RZU)

    CASE DEFAULT
      CALL FATAL_ERROR("INVALID DENSITY FUNCTION SELECTED:",&
           & "   "//TRIM(SEA_WATER_DENSITY_FUNCTION) )
    END SELECT

    DENSITY = DENSITY*1000.0_SP

    END SUBROUTINE CAL_DENSITY

!-----------------------------------------------------
!
!-----------------------------------------------------
    SUBROUTINE MLDEP
    
    IMPLICIT NONE

    INTEGER  :: I, J
    REAL(SP) :: DEPS, DDENSITY

    I = N_ORDER

    DEPS = 1.e-4 

    DO J=1, NZ(I)-1
      DDENSITY = ABS(DENSITY(J+1)-DENSITY(J))
      IF(DDENSITY>DEPS) EXIT
    ENDDO
    ML = J

    END SUBROUTINE MLDEP

!-----------------------------------------------------
!
!-----------------------------------------------------
    SUBROUTINE FREE_CONVECTION

    IMPLICIT NONE

    INTEGER :: I, J

    I = N_ORDER

    DO J=2, NZ(I)
      IF(DENSITY(J)>DENSITY(J-1)) EXIT
      IF(DENSITY(J)<DENSITY(J-1)) CALL MIX(J)
    ENDDO

!    CALL CAL_DENSITY

    END SUBROUTINE FREE_CONVECTION

!-----------------------------------------------------
!
!-----------------------------------------------------
    SUBROUTINE BULK_MIXING   !(check,iunit)

!    use all_vars, only : iint
    IMPLICIT NONE

    INTEGER :: I, J, MML
    REAL(SP) :: DELR, DS2, RO, G, RV

!    logical check
!    integer iunit

    G  = 9.8_SP
    RO = 1.024E3_SP

    I = N_ORDER

    DO J= ML, NZ(I)-1
      DELR = (DENSITY(J+1)-DENSITY(1))/RO
      DS2  = (UUU(J+1)-UUU(1))**2 + (VVV(J+1)-VVV(1))**2 + 1.0E-8

      RV   = G*DELR*FLOAT(J)*DDDZ/DS2

!      if(check) write(iunit+4,'(2i10,4f20.12)') iint, j, rv, DELR, FLOAT(J)*DDDZ, DS2

      IF(RV>BULK_RICHARDSON) RETURN

      CALL MIX(J+1)
!      CALL CAL_DENSITY
    ENDDO

    END SUBROUTINE BULK_MIXING

!-----------------------------------------------------
!
!-----------------------------------------------------
    SUBROUTINE GRADIENT_MIXING    !(R)

    IMPLICIT NONE

    INTEGER  :: I, J, JS, KCD, J1, J2
    REAL(SP) :: G, RO, DV, RS, DDENSITY
    REAL(SP) :: R(NZ(N_ORDER))

    KCD = 0
    G   = 9.8_SP
    RO  = 1.024E3_SP

    I = N_ORDER

    J1 = 1
    J2 = NZ(I)-1
 10 CONTINUE

    DO J=J1, J2
      DDENSITY = DENSITY(J+1) - DENSITY(J)
      IF(DDENSITY<1.0E-3) DDENSITY = 1.0E-3
      DV = (UUU(J+1)-UUU(J))**2 + (VVV(J+1)-VVV(J))**2
      IF(DV<1.0E-6) DV = 1.0E-6
      R(J) = G*DDDZ*DDENSITY/(DV*RO)
    ENDDO  

    RS = R(1)
    JS = 1
    DO J=2,NZ(I)-1
      IF(R(J)<RS) THEN
        RS = R(J)
        JS = J
      ENDIF 
    ENDDO

    IF(RS>GRADIENT_RICHARDSON) RETURN

    IF(JS>=KCD) KCD = JS + 1
    CALL STIR(GRADIENT_RICHARDSON,RS,JS)
!    CALL CAL_DENSITY

    J1 = JS - 2
    IF(J1<1) J1 = 1
    J2 = JS + 2
    IF(J2>NZ(I)-1) J2 = NZ(I)-1
    GOTO 10

    END SUBROUTINE GRADIENT_MIXING

!-----------------------------------------------------
!
!-----------------------------------------------------
    SUBROUTINE MIX(K_DEP)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: K_DEP
    INTEGER             :: K
    REAL(SP)            :: BT, BS, BU, BV, BD
    
    BT = 0.0_SP
    BS = 0.0_SP
    BU = 0.0_SP
    BV = 0.0_SP
    BD = 0.0_SP
    DO K=1, K_DEP
      BT = BT + TTT(K)
      BS = BS + SSS(K)
      BU = BU + UUU(K)
      BV = BV + VVV(K)
      BD = BD + DENSITY(K)
    ENDDO
    BT = BT / FLOAT(K_DEP)
    BS = BS / FLOAT(K_DEP)
    BU = BU / FLOAT(K_DEP)
    BV = BV / FLOAT(K_DEP)
    BD = BD / FLOAT(K_DEP)
    DO K=1, K_DEP
      TTT(K) = BT
      SSS(K) = BS
      UUU(K) = BU
      VVV(K) = BV
      DENSITY(K) = BD
    ENDDO

    END SUBROUTINE MIX

!-----------------------------------------------------
!
!-----------------------------------------------------
    SUBROUTINE STIR(RC,R,J)

    IMPLICIT NONE

    INTEGER  :: J
    REAL(SP) :: RC, R, RCON, RNEW, F, DA

    RCON = 0.02_SP + (RC-R)/2.0_SP
    RNEW = RC + RCON/5.0_SP

    F  = 1.0_SP - R/RNEW

    DA = (UUU(J+1)-UUU(J))*F/2.0_SP
    UUU(J+1) = UUU(J+1) - DA
    UUU(J)   = UUU(J) + DA

    DA = (VVV(J+1)-VVV(J))*F/2.0_SP
    VVV(J+1) = VVV(J+1) - DA
    VVV(J)   = VVV(J) + DA

    DA = (TTT(J+1)-TTT(J))*F/2.0_SP
    TTT(J+1) = TTT(J+1) - DA
    TTT(J)   = TTT(J) + DA

    DA = (SSS(J+1)-SSS(J))*F/2.0_SP
    SSS(J+1) = SSS(J+1) - DA
    SSS(J)   = SSS(J) + DA

    DA = (DENSITY(J+1)-DENSITY(J))*F/2.0_SP
    DENSITY(J+1) = DENSITY(J+1) - DA
    DENSITY(J)   = DENSITY(J) + DA

    END SUBROUTINE STIR

    FUNCTION SG0(S)
    IMPLICIT NONE
    REAL(SP) :: S, SG0

!   A sigma-0 subroutine neede by the sigma-t subroutine;
!   taken from seaprop.

!   sigma-0 knudsen

    SG0 = ((6.76786136E-6_SP*S-4.8249614E-4_SP)*S+0.814876577_SP)*S-0.0934458632_SP
    
    RETURN
    END FUNCTION SG0

    FUNCTION SGT(T,S,SG)
    IMPLICIT NONE
    REAL(SP)  :: T, S, SG, SGT

!   A sigma-t subroutine taken from seaprop;
!   sigma-t knudsen

    SG = SG0(S)
    SGT = ((((-1.43803061e-7_SP*T-1.98248399e-3_SP)*T-0.545939111_SP)*T     &
       +4.53168426_SP)*T)/(T+67.26_SP)+((((1.667e-8_SP*T-8.164e-7_SP)*T     &
       +1.803e-5_SP)*T)*SG+((-1.0843e-6_SP*T+9.8185e-5_SP)*T-4.7867e-3_SP)*T   &
       +1.0_SP)*SG
    
    RETURN
    END FUNCTION SGT

# endif
  END MODULE MOD_PWP
