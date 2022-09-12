# if defined (WAVE_CURRENT_INTERACTION)

!******************************************************************
!
!   SUBROUTINE SWCOMP(AC1,IT,CROSS)                                  
   SUBROUTINE SWCOMP(IT)                                  
!
!******************************************************************
!
!     The aim of this model is to simulate the wave energy in
!     shallow water areas. In the subroutine SWCOMP the main processes
!     taking place in the shallow water zone are determined in
!     several subroutines.
!     The input for this subroutine comes from SWANPRE1 and SWANPRE2.
!     The output is send to the subroutines SWANOUT1 and SWANOUT2.
!     The output consist of some characteristic
!     wave parameters and the wave action density. The equations are
!     all based on the action density N which is a function of the
!     spatial position (x,y), the relative frequency (s) and the
!     spectral direction (d).
!
!
!     Keywords:
!     Action density, propagation terms, refraction, reflection,
!     white capping, wave breaking, bottom friction, nonlinear
!     and nonhomogeneous wind- and current-fields, wave blocking,
!     fully spectral description, nonlinear wave-wave interaction,
!     higher order upwind schemes, flux limiting, SIP solver
!
!     Coefficients for the arrays:
!     -----------------------------
!                         default
!                         value:
!
!     PBOT(1)   = CFC      0.005    (Collins equation)
!     PBOT(2)   = CFW      0.01     (Collins equation)
!     PBOT(3)   = GAMJNS   0.0038   (Jonswap equation)
!     PBOT(4)   = MF      -0.08     (Madsen equation)
!     PBOT(5)   = KN       0.05     (bottom roughness)
!
!     ISURF                1        (Constant breaking coefficient)
!                          2        (variable breaking coefficient
!                                    according to Nelson (1994))
!     PSURF(1)  = ALFA    SWCOMP 1.0      (Battjes Janssen)
!     PSURF(2)  = GAMMA    0.8      (Breaking criterium)
!
!     PWCAP(1)  = ALFAWC   2.36e-5  (Empirical coefficient)
!     PWCAP(2)  = ALFAPM   3.02E-3  (Alpha of Pierson Moskowitz frequency)
!
!     PWIND(1)  = CF10     188.0    (second generation wind growth model)
!     PWIND(2)  = CF20     0.59     (second generation wind growth model)
!     PWIND(3)  = CF30     0.12     (second generation wind growth model)
!     PWIND(4)  = CF40     250.0    (second generation wind growth model)
!     PWIND(5)  = CF50     0.0023   (second generation wind growth model)
!     PWIND(6)  = CF60    -0.2233   (second generation wind growth model)
!     PWIND(7)  = CF70     0.       (second generation wind growth model)
!     PWIND(8)  = CF80    -0.56     (second generation wind growth model)
!     PWIND(9)  = RHOAW    0.00125  (density air / density water)
!     PWIND(10) = EDMLPM   0.0036   (limit energy Pierson Moskowitz)
!     PWIND(11) = CDRAG    0.0012   (drag coefficient)
!     PWIND(12) = UMIN     1.0      (minimum wind velocity)
!     PWIND(13) = PMLM     0.13     (  )
!
!     PNUMS(1)  = DREL     relative error in Hs and Tm
!     PNUMS(2)  = DHABS    absolute error in Hs
!     PNUMS(3)  = DTABS    absolute error in Tm
!     PNUMS(4)  = NPNTS    number of points were accuracy is reached
!
!     PNUMS(5)  = NOT USED
!     PNUMS(6)  = CDD      blending parameter for finite differences
!                          in theta space
!     PNUMS(7)  = CSS      blending parameter for finite differences
!                          in sigma space
!     PNUMS(8)  = NUMFRE  SWCOMP numerical scheme in frequency space :
!                          1) implicit scheme
!                          2) explicit scheme CFL limited
!                          3) explicit scheme filter after iteration
!     PNUMS(9)  = DIFFC    if explicit scheme is used, then numerical
!                          diffusion coefficient can be chosen
!     PNUMS(12) = EPS2     termination criterion in relative sense for a
!                          penta-diagonal solver
!     PNUMS(13) = OUTP     request for output for a penta-diagonal solver
!     PNUMS(14) = NITER    maximum number of iterations for a penta-diagonal
!                          solver
!     PNUMS(15) = DHOVAL   global error in Hs
!               = CURVAT   curvature of Hs meant for convergence check
!     PNUMS(16) = DTOVAL   global error in Tm01
!     PNUMS(17) = CDLIM    coefficient of limitation of Ctheta
!     PNUMS(18) = FROUDMAX maximum Froude number for reduction of currents
!     PNUMS(19) = CFL      CFL criterion for option explicit scheme
!                          in frequency space (see PNUMS(8))
!     PNUMS(20) = GRWMX    maximum growth in spectral bin
!
!     PNUMS(21) = STOPC    type of stopping criterion:
!                          0: standard SWAN based on relative and global
!                             errors of Hs and Tm01,
!                          1: based on absolute, relative and curvature
!                             errors of Hs
!
!     PNUMS(30) = ALFA     relaxation parameter for under-relaxation method
!
!     arrays for the 4-wave interactions:
!
!     WWINT ( 1 = IDP    WWAWG ( = AGW1    WWSWG ( = SWG1
!             2 = IDP1           = AWG2            = SWG2
!             3 = IDM            = AWG3            = SWG3
!             4 = IDM1           = AWG4            = SWG4
!             5 = ISP            = AWG5            = SWG5
!             6 = ISP1           = AWG6            = SWG6
!             7 = ISM            = AWG7            = SWG7
!             8 = ISM1           = AWG8 )          = SWG8  )
!             9 = ISLOW
!             10= ISHGH
!             11= ISCLW
!             12= ISCHG
!             13= IDLOW
!             14= IDHGH
!             15= MSC4MI
!             16= MSC4MA
!             17= MDC4MI
!             18= MDC4MA
!             19= MSCMAX
!             20= MDCMAX
!             21= IDPP
!             22= IDMM
!             23= ISPP
!             24= ISMM )
!
!******************************************************************

   USE TIMECOMM                                                        
   USE OCPCOMM1                                                        
   USE OCPCOMM2                                                        
   USE OCPCOMM3                                                        
   USE OCPCOMM4                                                        
   USE SWCOMM1                                                         
   USE SWCOMM2                                                         
   USE SWCOMM3                                                         
   USE SWCOMM4 
   USE M_GENARR                                                        
   USE m_constants                                                     
   USE m_xnldata                                                       
   USE m_fileio  

#  if defined (EXPLICIT)
   USE MOD_ACTION_EX
#  else   
   USE MOD_ACTION_IM
#  endif   
#  if defined (MULTIPROCESSOR)
   USE MOD_PAR
   USE MOD_NESTING, ONLY : NESTING_ON_WAVE, NESTING_TYPE_WAVE, SET_VAR_WAVE, SET_VAR_WAVE_AC2
#  endif    
!   USE ALL_VARS    !, ONLY : M,MT,NGL,SERIAL,PAR,MYID,NPROCS,MSR,     &
!                   !     IOBCN,I_OBC_N,AC2,COMPDA,ISBCE 
   USE VARS_WAVE
#  if defined (WAVE_SETUP)
   USE MOD_WAVESETUP
#  endif   
 
#  if defined (PLBC)
   USE MOD_PERIODIC_LBC
#  endif   

   IMPLICIT NONE                                                      

   REAL    SIGLOW
!
!     ************************************************************************
!     *                                                                      *
!     *                  MAIN SUBROUTINE OF COMPUTATIONAL PART               *
!     *                                                                      *
!     *                               -- SWCOMP --                           *
!     *                                                                      *
!     ************************************************************************
!
   INTEGER :: ITER  ,IX    ,IY    ,IS    ,IT                           
   INTEGER :: IP, IDC, ISC                                             
   INTEGER :: INOCNV                                                   
   INTEGER :: INOCNT                                                   
!
   REAL ::  DDX   ,DDY   ,ACCUR ,XIS   ,SNLC1 ,DAL1  ,DAL2  ,DAL3      
!
   LOGICAL :: PRECOR
!
   INTEGER :: MNISL, MXNFL, MXNFR, NPFL, NPFR, NWETP, IDEBUG           
   REAL    :: FRAC                                                     
   PARAMETER (IDEBUG=0)                                                
!
   INTEGER   ISTAT, IF1, IL1                                           
   CHARACTER*20 NUMSTR, CHARS(1)                                       
   CHARACTER*80 MSGSTR                                                 
!
   INTEGER IARR(10)                                                    
   REAL     ARR(10)                                                    

   INTEGER JDUM, JS, JE, JWFRS, JWFRE, INCJ, JJ, JNODE,             &
           JSD, JED, IE, INCI, II, III, LSTCP                       

   INTEGER IX1, IX2, IY1, IY2                                        
!
!  Add variables for the XNL interface (quadruplet interaction)        
   INTEGER :: IXGRID, IXQUAD, IQERR                                    
!
   INTEGER :: MSTPDA                                          
!JQI   INTEGER :: CROSS(2,M)
!
!JQI   REAL AC1(MDC,MSC,0:MT)                                         

   REAL WWAWG(8), WWSWG(8)                                             

   INTEGER, DIMENSION(:), ALLOCATABLE :: IDCMIN, IDCMAX,            &
                                         ISCMIN, ISCMAX              
   INTEGER WWINT(24)                                                  

   REAL, DIMENSION(:,:,:), ALLOCATABLE :: CAX,CAY,CAX1,CAY1,CAS,CAD
   REAL, DIMENSION(:,:), ALLOCATABLE :: CGO,KWAVE                      
   REAL, DIMENSION(:,:), ALLOCATABLE :: ALIMW                          
   REAL, DIMENSION(:,:), ALLOCATABLE :: UE,SA1,SA2,SFNL               
   REAL, DIMENSION(:,:), ALLOCATABLE :: DA1C,DA1P,DA1M,             &
                                        DA2C,DA2P,DA2M,DSNL           
   REAL, DIMENSION(:,:,:), ALLOCATABLE :: MEMNL4                       
   REAL, DIMENSION(:,:,:), ALLOCATABLE :: OBREDF                       
   REAL, DIMENSION(:,:), ALLOCATABLE :: REFLSO                         
   REAL, DIMENSION(:), ALLOCATABLE :: HSAC1,HSAC2,SACC1,SACC2          
   REAL, DIMENSION(:), ALLOCATABLE :: HSAC0,HSDIFC                     
   REAL, DIMENSION(:,:), ALLOCATABLE :: SETPDA                         
   LOGICAL, DIMENSION(:), ALLOCATABLE :: GROWW                         
   REAL, ALLOCATABLE :: SWTSDA(:,:,:,:)                                
   REAL, DIMENSION(:,:,:), ALLOCATABLE :: SWMATR                       
   INTEGER, ALLOCATABLE :: ISLMIN(:), NFLIM(:), NRSCAL(:)    
   
   INTEGER :: J1,J2,KWIND,KWCAP,KQUAD,IDTOT,IOB,IOB_NODE
   REAL    :: GRWOLD,ALFAT         

   REAL(SP), ALLOCATABLE :: N32_TMP(:),AC2_TMP(:)
   REAL(SP), ALLOCATABLE :: N32_TTMP(:,:,:),AC2_TTMP(:,:,:)
   INTEGER, ALLOCATABLE :: PASS(:)
   REAL,    ALLOCATABLE :: PASS_TMP(:)
   REAL,    ALLOCATABLE :: COMPDA_TMP1(:)
   
   REAL :: AWX_W,AWY_W,RWX,RWY,WIND10,THETAW
!   LOGICAL, DIMENSION(:), ALLOCATABLE :: PASS(:)
   REAL :: TY,XTMP,XTMP1,DELTX,DELTY,ALPHA1,ALPHA2,ALPHA3
!     Add variables for OMP thread parameters.                            
      INTEGER I1GRD,I2GRD,I1MYC,I2MYC
      
   INTEGER :: IERR 
   REAL :: DEP_TMP                                       
!
!-----------------------------------------------------------------------
!                      End of variable definition
!-----------------------------------------------------------------------
!
   IF(IT == 1 .AND. ITEST >= 1)THEN                                
     WRITE(PRINTF,333) 'SWAN'                                          
333  FORMAT(/,                                                          &
     '----------------------------------------------------------------' &
     ,/,                                                                &
     '                  COMPUTATIONAL PART OF ', A                      &
     ,/,                                                                &
     '----------------------------------------------------------------' &
     ,/)
!
     IF(ONED)THEN                                                  
       WRITE(PRINTF,*) 'One-dimensional mode of SWAN is activated'     
     ENDIF                                                             
!
     WRITE(PRINTF,7001) NGL
7001 FORMAT(' The Numbers of Cell: NGC   ',I12)
     WRITE(PRINTF,7002) MCGRD                                          
7002 FORMAT(' The Numbers of Node: MGC   ',I12)                     
     WRITE(PRINTF,7101) MSC,MDC
7101 FORMAT('                      : MSC    ',I12  ,' MDC   ',I12)
     WRITE(PRINTF,7201) MTC                                            
7201 FORMAT('                      : MTC    ',I12)                     
     WRITE(PRINTF,7301) NSTATC, ITERMX
7301 FORMAT('                      : NSTATC ',I12  ,' ITERMX',I12)
     WRITE(PRINTF,7013) ITFRE,IREFR
7013 FORMAT(' Propagation flags    : ITFRE  ',I12  ,' IREFR ',I12)
     WRITE(PRINTF,7014) IBOT,ISURF
7014 FORMAT(' Source term flags    : IBOT   ',I12  ,' ISURF ',I12)
     WRITE(PRINTF,7114) IWCAP,IWIND
7114 FORMAT('                      : IWCAP  ',I12  ,' IWIND ',I12)
     WRITE(PRINTF,7015) ITRIAD,IQUAD
7015 FORMAT('                      : ITRIAD ',I12  ,' IQUAD ',I12)
# if defined (VEGETATION)
     WRITE(PRINTF,7016) IVEG
7016 FORMAT('                      : IVEG   ',I12)
# endif
     WRITE(PRINTF,7104) EXP(ALOG(SHIG/SLOW)/REAL(MSC-1))-1.,DDIR/DEGRAD
7104 FORMAT(' Spectral bin         : df/f   ',E12.4,' DDIR  ',E12.4)   
     WRITE(PRINTF,7003) GRAV_W, RHO_W                                      
7003 FORMAT(' Physical constants   : GRAV   ',E12.4,' RHO   ',E12.4)   
     WRITE(PRINTF,7027) U10 , WDIC/DEGRAD
7027 FORMAT(' Wind input           : WSPEED ',E12.4,' DIR   ',E12.4)
     WRITE(PRINTF,7123) PWTAIL(1),PWTAIL(2)
7123 FORMAT(' Tail parameters      : E(f)   ',E12.4,' E(k)  ',E12.4)
     WRITE(PRINTF,7133) PWTAIL(3),PWTAIL(4)
7133 FORMAT('                      : A(f)   ',E12.4,' A(k)  ',E12.4)
     WRITE(PRINTF,8013) PNUMS(1), PNUMS(4)                             
8013 FORMAT(' Accuracy parameters  : DREL   ',E12.4,' NPNTS ',E12.4)   
     IF(PNUMS(21) == 0.)THEN                                         
       WRITE(PRINTF,8213) PNUMS(15),PNUMS(16)                         
     ELSE IF(PNUMS(21) == 1.)THEN                                    
       WRITE(PRINTF,8214) PNUMS(2),PNUMS(15)                          
     END IF                                                            
8213 FORMAT('                      : DHOVAL ',E12.4,' DTOVAL',E12.4)   
8214 FORMAT('                      : DHABS  ',E12.4,' CURVAT',E12.4)   
     WRITE(PRINTF,3613) PNUMS(20)
3613 FORMAT('                      : GRWMX  ',E12.4)
     WRITE(PRINTF,8508) WLEV, DEPMIN                                   
8508 FORMAT(' Drying/flooding      : LEVEL  ',E12.4,' DEPMIN',E12.4)   

     IF(BNAUT)THEN                                                   
       WRITE (PRINTF,8510) 'nautical '                                
     ELSE                                                              
       WRITE (PRINTF,8510) 'Cartesian'                                
     END IF                                                            
8510 FORMAT(' The ',A9,' convention for wind and wave directions is used')        

     WRITE(PRINTF,8513) PNUMS(7), PNUMS(6)                             
8513 FORMAT(' Scheme freq. space   : CSS    ',E12.4,' CDD   ',E12.4)   
!
     IF((DYNDEP .OR. ICUR == 1) .AND. INT(PNUMS(8)) == 1)THEN      
       WRITE(PRINTF,*) 'Solver is SIP'                                
       WRITE(PRINTF,2213) PNUMS(12), INT(PNUMS(13))
2213   FORMAT('                      : EPS2   ',E12.4,' OUTPUT',I12)
       WRITE(PRINTF,8223) INT(PNUMS(14))
8223   FORMAT('                      : NITER  ',I12)
     END IF
!
     IF(ICUR > 0)THEN
       WRITE (PRINTF,7115) 'on'
     ELSE
       WRITE (PRINTF,7115) 'off'
     ENDIF
7115 FORMAT(' Current is ', A3)
!
     IF(IQUAD > 0)THEN
       WRITE(PRINTF,8413) IQUAD
8413   FORMAT(' Quadruplets          : IQUAD  ',I12)
       IF(IQUAD <= 3 .OR. IQUAD == 8)THEN                           
         WRITE(PRINTF,8414) PQUAD(1), PQUAD(2)                       
8414     FORMAT('                      : LAMBDA ',E12.4,' CNL4  ',E12.4)
         WRITE(PRINTF,8415) PQUAD(3), PQUAD(4)                       
8415     FORMAT('                      : CSH1   ',E12.4,' CSH2  ',E12.4)
         WRITE(PRINTF,8416) PQUAD(5)                                 
8416     FORMAT('                      : CSH3   ',E12.4)             
       END IF                                                         
       WRITE(PRINTF,8417) PTRIAD(3)                                   
8417   FORMAT(' Maximum Ursell nr for Snl4 :  ',E12.4)                
     ELSE
       WRITE (PRINTF, *) 'Quadruplets is off'
     END IF                                                            

     IF(ITRIAD > 0)THEN                                             
       WRITE(PRINTF,9413) ITRIAD, PTRIAD(1)
9413   FORMAT(' Triads               : ITRIAD ',I12  ,' TRFAC ',E12.4)
       WRITE(PRINTF,9414) PTRIAD(2), PTRIAD(4)                        
9414   FORMAT('                      : CUTFR  ',E12.4,' URCRI ',E12.4)
       WRITE(PRINTF,9415) PTRIAD(5)                                   
9415   FORMAT(' Minimum Ursell nr for Snl3 :  ',E12.4)                
     ELSE
       WRITE (PRINTF, *) 'Triads is off'
     END IF                                                            

     IF(IBOT == 2)THEN
       WRITE(PRINTF,7005) PBOT(2),PBOT(1)
7005   FORMAT(' Collins (`72)        : CFW    ',E12.4,' CFC   ',E12.4)
     ELSE IF(IBOT == 3)THEN
       WRITE(PRINTF,7335) PBOT(4),PBOT(5)
7335   FORMAT(' Madsen et al. (`84)  : MF     ',E12.4,' KN    ',E12.4)
     ELSE IF(IBOT == 1)THEN
       WRITE(PRINTF,7325) PBOT(3)
7325   FORMAT(' JONSWAP (`73)        : GAMMA  ',E12.4)
     ELSE
       WRITE (PRINTF, *) 'Bottom friction is off'
     ENDIF
!
# if defined (VEGETATION)
!
      IF (IVEG.EQ.1) THEN
         WRITE(PRINTF,7006)
 7006    FORMAT(' Vegetation due to Dalrymple (1984)')
      ELSE
         WRITE (PRINTF, *) 'Vegetation is off'
      ENDIF
!
# endif
     IF(IWCAP == 1)THEN
       WRITE(PRINTF,6005) PWCAP(1),PWCAP(2)
6005   FORMAT(' W-cap Komen (`84)    : EMPCOF ',E12.4,' APM   ',E12.4)
     ELSE IF(IWCAP == 2)THEN
       WRITE(PRINTF,6335) PWCAP(3),PWCAP(4)
6335   FORMAT(' W-cap Janssen (`90)  : CFJANS ',E12.4,' DELTA ',E12.4)
     ELSE IF(IWCAP == 3)THEN
       WRITE(PRINTF,6135) PWCAP(5)
6135   FORMAT(' W-cap Longuet-Higgins: CFLHIG ',E12.4)
     ELSE IF(IWCAP == 4)THEN
       WRITE(PRINTF,6136) PWCAP(6), PWCAP(7)
6136   FORMAT(' W-cap Battjes/Janssen: BJSTP  ',E12.4,' BJALF ',E12.4)
     ELSE IF(IWCAP == 5)THEN
       WRITE(PRINTF,6136) PWCAP(6), PWCAP(7)
       WRITE(PRINTF,6137) PWCAP(8)
6137   FORMAT('                      : KCONV  ',E12.4)
     ELSE IF(IWCAP == 6)THEN                                         
       WRITE(PRINTF,6138) PWCAP(12), PWCAP(13)                         
6138   FORMAT(' W-cap cum. steepness : CST    ',E12.4,' POW   ',E12.4)
     ELSE IF(IWCAP == 7)THEN       
       WRITE(PRINTF,6139) PWCAP(1), PWCAP(12)                          
       WRITE(PRINTF,6140) PWCAP(9), PWCAP(11)                          
6139   FORMAT(' W-cap Alves-Banner   : CDS2   ',E12.4,' BR    ',E12.4)
6140   FORMAT('                      : POWST  ',E12.4,' POWK  ',E12.4)
     ELSE
       WRITE (PRINTF, *) 'Whitecapping is off'
     ENDIF
!
     IF(ISURF == 1)THEN
       WRITE(PRINTF,7012) PSURF(1),PSURF(2)
7012   FORMAT(' Battjes&Janssen (`78): ALPHA  ',E12.4,' GAMMA ',E12.4)
     ELSE IF(ISURF == 2)THEN                                        
       WRITE(PRINTF,7212) PSURF(1), PSURF(4), PSURF(5)
7212   FORMAT(' Nelson (`94): ALPHA  ',E12.4,' GAMmin ',E12.4, ' GAMmax ',E12.4)
     ELSE
       WRITE (PRINTF, *) 'Surf breaking is off'
     ENDIF
!
     IF(LSETUP > 0)THEN
       WRITE(PRINTF,7109) PSETUP(2)
7109   FORMAT(' Set-up               : SUPCOR ',E12.4)
     ELSE
       WRITE (PRINTF, *) 'Set-up is off'
     END IF
!
     IF(IDIFFR == 1)THEN                                             
       WRITE(PRINTF,7110) PDIFFR(1), NINT(PDIFFR(2))
7110   FORMAT(' Diffraction          : SMPAR  ',E12.4,' SMNUM ',I12)
     ELSE
       WRITE (PRINTF, *) 'Diffraction is off'
     ENDIF
!
     WRITE(PRINTF,7126) PWIND(14), PWIND(15)
7126 FORMAT(' Janssen (`89,`90)    : ALPHA  ',E12.4,' KAPPA ',E12.4)
     WRITE(PRINTF,7136) PWIND(16), PWIND(17)
7136 FORMAT(' Janssen (`89,`90)    : RHOA   ',E12.4,' RHOW  ',E12.4)
     WRITE(PRINTF,*)
     WRITE(PRINTF,1012) PWIND(1), PWIND(2)
1012 FORMAT(' 1st and 2nd gen. wind: CF10   ',E12.4,' CF20  ',E12.4)
     WRITE(PRINTF,1013) PWIND(3), PWIND(4)
1013 FORMAT('                      : CF30   ',E12.4,' CF40  ',E12.4)
     WRITE(PRINTF,1014) PWIND(5), PWIND(6)
1014 FORMAT('             SWCOMP         : CF50   ',E12.4,' CF60  ',E12.4)
     WRITE(PRINTF,1015) PWIND(7), PWIND(8)
1015 FORMAT('                      : CF70   ',E12.4,' CF80  ',E12.4)
     WRITE(PRINTF,1016) PWIND(9), PWIND(10)
1016 FORMAT('                      : RHOAW  ',E12.4,' EDMLPM',E12.4)
     WRITE(PRINTF,1017) PWIND(11), PWIND(12)
1017 FORMAT('                      : CDRAG  ',E12.4,' UMIN  ',E12.4)
     WRITE(PRINTF,1018) PWIND(13)
1018 FORMAT('                      : LIM_PM ',E12.4)
!
     WRITE(PRINTF,*)
     IF(ITEST > 2)THEN
       DO IS = 1, MSC
         WRITE(PRINTF,*)' IS and SPCSIG(IS)    :',IS,SPCSIG(IS)        
       ENDDO
       WRITE(PRINTF,*)
     ENDIF
   END IF
!
!  *** prepare ranges of spectral space, constants and      ***        
!  *** weight factors for nonlinear 4 wave interactions     ***        
!                                                                         
   IF(IQUAD >=1)                                              &
      CALL FAC4WW (XIS   ,SNLC1 ,DAL1  ,DAL2  ,DAL3  ,SPCSIG, &
		   WWINT ,WWAWG ,WWSWG )                               
!                                                                         
! *** Indexing and bounds for SWMAT arrays                                
!                                                                         
   JMATD = 1                                                           
   JMATR = 2                                                           
   JMATL = 3                                                           
   JMATU = 4                                                           
   JMAT5 = 5                                                           
   JMAT6 = 6                                                           
   JDIS0 = 7                                                           
   JDIS1 = 8                                                           
   JLEK1 = 9                                                           
   JAOLD = 10                                                          
   MSWMATR = 10                                                        
   JABIN = 1                                                           
   JABLK = 2                                                           

!----------------------------------------------------------------------   
!     Begin allocate shared arrays.                                       
!----------------------------------------------------------------------   
!
   ALLOCATE(HSAC1(M))                                              
   ALLOCATE(HSAC2(M))                                              
   ALLOCATE(SACC1(M))                                              
   ALLOCATE(SACC2(M))                                              
   ALLOCATE(HSAC0(M))                                              
   ALLOCATE(HSDIFC(M))                                             
!
   ALLOCATE(ISLMIN(M))                                             
   ALLOCATE(NFLIM(M))                                              
   ALLOCATE(NRSCAL(M))                                             
!
   IF(IQUAD >= 1)THEN
!    *** quadruplets ***
     IF(IQUAD >= 3)THEN                                          
!      *** prior to every iteration full directional domain ***
       ALLOCATE(MEMNL4(MDC,MSC,M),STAT=ISTAT)                      
       IF(ISTAT /= 0)THEN                                          
         CHARS(1) = NUMSTR(ISTAT,RNAN,'(I6)')                         
         CALL TXPBLA(CHARS(1),IF1,IL1)                                
         MSGSTR = 'Allocation problem: array MEMNL4 and return code is '//   &
	          CHARS(1)(IF1:IL1)                                          
         CALL MSGERR ( 4, MSGSTR )                                    
         RETURN                                                       
       END IF                                                          
     ELSE                                                              
!    *** iquad < 3 ***                                                 
       ALLOCATE(MEMNL4(0,0,0))                                         
     END IF
   ELSE
!    *** no quadruplets ***
     ALLOCATE(MEMNL4(0,0,0))                                           
   ENDIF

   ALLOCATE(SWTSDA(MDC,MSC,NPTSTA,MTSVAR))                             
   SWTSDA = 0.                                                         
!
!  *** In case of SETUP expand array for setup data ***                
!
   IF(LSETUP > 0)THEN                                               
     MSTPDA = 23                                                       
     ALLOCATE(SETPDA(M,MSTPDA)) 
   ELSE                                                                
     ALLOCATE(SETPDA(0,0))                                             
   END IF

!----------------------------------------------------------------------   
!     End allocate shared arrays.                                         
!----------------------------------------------------------------------   

!----------------------------------------------------------------------   
!     Begin initialization of shared arrays.                              
!----------------------------------------------------------------------   

   HSAC1 = 0.                                                          
   HSAC2 = 0.                                                          
   SACC1 = 0.                                                          
   SACC2 = 0.                                                          
   HSAC0 = 0.                                                          
   HSDIFC= 0.                                                          
   IF(IQUAD >= 3) MEMNL4 = 0.                                       
   IF(LSETUP > 0) SETPDA = 0.                                      

!----------------------------------------------------------------------   
!     End initialization shared arrays.                                   
!----------------------------------------------------------------------   

!----------------------------------------------------------------------   
!     Begin allocate private arrays.                                      
!----------------------------------------------------------------------   

   ALLOCATE(CAX(MDC,MSC,MICMAX))                                       
   ALLOCATE(CAY(MDC,MSC,MICMAX))                                       
   ALLOCATE(CAX1(MDC,MSC,MICMAX))                                      
   ALLOCATE(CAY1(MDC,MSC,MICMAX))                                      
   ALLOCATE(CAS(MDC,MSC,MICMAX))                                       
   ALLOCATE(CAD(MDC,MSC,MICMAX))                                       
   ALLOCATE(CGO(MSC,MICMAX))                                           
   ALLOCATE(KWAVE(MSC,MICMAX))                                         
   ALLOCATE(SWMATR(MDC,MSC,MSWMATR))                                   
   ALLOCATE(ALIMW(MDC,MSC))                                            
   ALLOCATE(GROWW(MDC*MSC))                                            
   ALLOCATE(IDCMIN(MSC))                                               
   ALLOCATE(IDCMAX(MSC))                                               
   ALLOCATE(ISCMIN(MDC))                                               
   ALLOCATE(ISCMAX(MDC))                                               
!  *** quadruplets ***
   IF(IQUAD >= 1)THEN
     ALLOCATE(UE(MSC4MI:MSC4MA,MDC4MI:MDC4MA))                         
     ALLOCATE(SA1(MSC4MI:MSC4MA,MDC4MI:MDC4MA))                        
     ALLOCATE(SA2(MSC4MI:MSC4MA,MDC4MI:MDC4MA))                        
     ALLOCATE(SFNL(MSC4MI:MSC4MA,MDC4MI:MDC4MA))                       
     IF(IQUAD == 1)THEN
!      *** semi-implicit calculation ***
       ALLOCATE(DA1C(MSC4MI:MSC4MA,MDC4MI:MDC4MA))                     
       ALLOCATE(DA1P(MSC4MI:MSC4MA,MDC4MI:MDC4MA))                     
       ALLOCATE(DA1M(MSC4MI:MSC4MA,MDC4MI:MDC4MA))                     
       ALLOCATE(DA2C(MSC4MI:MSC4MA,MDC4MI:MDC4MA))                     
       ALLOCATE(DA2P(MSC4MI:MSC4MA,MDC4MI:MDC4MA))                     
       ALLOCATE(DA2M(MSC4MI:MSC4MA,MDC4MI:MDC4MA))                     
       ALLOCATE(DSNL(MSC4MI:MSC4MA,MDC4MI:MDC4MA))                     
     ELSE                                                              
!    *** iquad > 1 ***                                                 
       ALLOCATE(DA1C(0,0))                                             
       ALLOCATE(DA1P(0,0))                                             
       ALLOCATE(DA1M(0,0))                                             
       ALLOCATE(DA2C(0,0))                                             
       ALLOCATE(DA2P(0,0))                                             
       ALLOCATE(DA2M(0,0))                                             
       ALLOCATE(DSNL(0,0))                                             
     END IF                                                            
   ELSE
!    *** no quadruplets ***
     ALLOCATE(UE(0,0))                                                 
     ALLOCATE(SA1(0,0))                                                
     ALLOCATE(SA2(0,0))                                                
     ALLOCATE(SFNL(0,0))                                               
     ALLOCATE(DA1C(0,0))                                               
     ALLOCATE(DA1P(0,0))                                               
     ALLOCATE(DA1M(0,0))                                               
     ALLOCATE(DA2C(0,0))                                               
     ALLOCATE(DA2P(0,0))                                               
     ALLOCATE(DA2M(0,0))                                               
     ALLOCATE(DSNL(0,0))                                               
   END IF
!
!  *** for obstacles, to store the transmission coefficients ***
!  *** and contribution to the source terms                  ***
!
   ALLOCATE(OBREDF(MDC,MSC,2))                                         
   ALLOCATE(REFLSO(MDC,MSC))                                           

!----------------------------------------------------------------------   
!     End allocate private arrays.                                        
!----------------------------------------------------------------------   

!----------------------------------------------------------------------   
!     Begin initialization of private arrays.                             
!----------------------------------------------------------------------   

   CAX    = 0.                                                         
   CAY    = 0.                                                         
   CAX1   = 0.                                                         
   CAY1   = 0.                                                         
   CAS    = 0.                                                         
   CAD    = 0.                                                         
   CGO    = 0.                                                         
   KWAVE  = 0.                                                         
   SWMATR = 0.                                                         
   ALIMW  = 0.                                                         
   IF(IQUAD >= 1)THEN                                              
     UE   = 0.                                                         
     SA1  = 0.                                                         
     SA2  = 0.                                                         
     SFNL = 0.                                                         
     IF(IQUAD == 1)THEN                                            
       DA1C = 0.                                                       
       DA1P = 0.                                                       
       DA1M = 0.                                                       
       DA2C = 0.                                                       
       DA2P = 0.                                                       
       DA2M = 0.                                                       
       DSNL = 0.                                                       
     END IF                                                            
   END IF                                                              

!----------------------------------------------------------------------   
!     End initialization private arrays.                                  
!----------------------------------------------------------------------   

!
!  *** initialise values for determining the accuracy that ***
!  *** has been reached                                    ***
!  *** This is done in parallel within OpenMP environment  ***         
!
   CALL INSAC (SPCSIG,COMPDA(1,JDP2),HSAC2,SACC2)   
!
!  *** To obtain a first estimate of energy density in a    ***
!  *** gridpoint considered we run the SWAN model (in case  ***
!  *** of active wind) in a second generation mode first.   ***
!  *** After 1st iteration, the options, as defined by the  ***
!  *** user, are re-activated.                              ***
!  *** This first guess is not used in nonstationary        ***
!  *** computations (NSTATC>0), or if a restart file was    ***
!  *** used (ICOND=4)                                       ***
!
   IF(IWIND >= 3 .AND. NSTATC == 0 .AND. ICOND /= 4)THEN           
!     --- first guess will be used
      PRECOR = .TRUE.
   ELSE
      PRECOR = .FALSE.
   END IF
!
!     *** call initialization procedure of XNL to create *.BQF ***        
!     *** interaction files                                    ***        
!                                                                         
   IF(IQUAD == 51 .OR. IQUAD == 52 .OR. IQUAD == 53)THEN 
     CALL init_constants                                              
     IXQUAD = IQUAD - 50                                              
     IF(MSR)THEN                                        
       WRITE(SCREEN,*) 'GurboQuad initialization'                    
       WRITE(SCREEN,*) 'gravity               :', GRAV_W               
       WRITE(SCREEN,*) 'pftail                :', PWTAIL(1)          
       WRITE(SCREEN,*) 'number of sigma values:', MSC                
       WRITE(SCREEN,*) 'number of directions  :', MDC                
       WRITE(SCREEN,*) 'IQ_QUAD               :', IXQUAD             
     END IF                                                           
!                                                                         
     IXGRID = 3                                                       
     J1 = I1GRD                                                       
     IF(J1 == 1) J1 = 2                                              
     J2 = I2GRD                
     CALL xnl_init(SPCSIG  , SPCDIR(:,1)*180./PI_W , MSC , MDC ,        &
                   -PWTAIL(1), GRAV_W   , COMPDA(J1:J2,JDP2) ,          &
!                   J2-J1+1   , IXQUAD , IXGRID   ,INODE,IQERR )       
                   J2-J1+1   , IXQUAD , IXGRID   ,IQERR )       
   END IF                                                              

   DO 450 ITER = 1, ITERMX                                             

!    initialise local (thread private) counter for SIP solver          
     INOCNT = 0                                                        
!
!    initialise Dissipation and Leak to 0 for each iteration
!    this is done in parallel within OpenMP environment                
!
     DO IP = 1,M                                               
       COMPDA(IP,JDISS) = 0.
       COMPDA(IP,JLEAK) = 0.
     ENDDO
!
!    initialise Ursell number to 0 for each iteration                  
!    this is done in parallel within OpenMP environment                
!
     IF(ITRIAD > 0)THEN
       DO IP = 1,M                                            
         COMPDA(IP,JURSEL) = 0.
       ENDDO
     ENDIF
!                                                                         
!    *** IQUAD = 3: the nonlinear wave interactions are     ***        
!    *** calculated just once for an iteration. First,      ***        
!    *** set the auxiliary array equal zero before a        ***        
!    *** new iteration                                      ***        
!    *** This is done in parallel within OpenMP environment ***        
!                                                                         
     IF(IQUAD >= 3)THEN                                          
       DO IP = 1,M                                             
         DO ISC = 1,MSC                                                
           DO IDC = 1,MDC                                              
             MEMNL4(IDC,ISC,IP)=0.                                     
           END DO                                                      
         END DO                                                        
       END DO                                                          
     END IF                                                            
!                                                                         
!    *** If a current is present and a penta-diagonal solver    ***    
!    *** is employed, it is possible that the solver does       ***    
!    *** not converged. For this, the counter INOCNV represents ***    
!    *** the number of geographical points in which the solver  ***    
!    *** did not converged                                      ***    
!                                                                         
     INOCNV = 0                                                        
!
     IF(ITEST >= 30 .OR. IDEBUG == 1)THEN
       ISLMIN(:) = 9999                                               
       NFLIM (:) = 0                                                  
       NRSCAL(:) = 0                                                  
       IARR = 0                                                       
       ARR  = 0.                                                      
     END IF
!
     IF(PRECOR)THEN                                                
!       *** third generation wave input ***
       IF(ITER == 1)THEN
!        *** save settings of 3rd generation model             ***   
!        *** bottom friction, surf breaking and triads may     ***
!        *** still active                                      ***
         KWIND  = IWIND
         KWCAP  = IWCAP
         KQUAD  = IQUAD
!        ***  save maximum change per bin and under-relaxation ***   
         GRWOLD = PNUMS(20)
         ALFAT  = PNUMS(30)                                          

!        first guess settings                                        

!        if 1st generation is to be used as first guess, replace     
!        the next statement by IWIND = 1                             
         IWIND  = 2
         IWCAP  = 0
         IQUAD  = 0
         PNUMS(20) = 1.E22
!        ***  under-relaxation parameter is PNUMS(30) and            
!             temporarily set to zero                                
         PNUMS(30) = 0.                                              
       ELSE IF ( ITER == 2 ) THEN
         IWIND  = KWIND
         IWCAP  = KWCAP
         IQUAD  = KQUAD
         PNUMS(20) = GRWOLD
         PNUMS(30) = ALFAT                                           
       ENDIF
!
     ENDIF

     IF(PRECOR .AND. ITER <= 2)THEN                              
       WRITE(PRINTF,*)' -----------------------------------------',     &
                      '----------------------'
       IF(ITER == 1)THEN
         WRITE(PRINTF,*) ' First guess by 2nd generation model',        &
	                 ' flags for first iteration:'
       ELSE IF(ITER == 2)THEN
         WRITE(PRINTF,*) ' Options given by user are activated',        &
	                 ' for proceeding calculation:'
       ENDIF
       WRITE(PRINTF,2001) ITER, PNUMS(20), PNUMS(30)
2001   FORMAT('  ITER    ',I4,' GRWMX    ',E12.4,' ALFA     ',E12.4)
       WRITE(PRINTF,2002) IWIND, IWCAP, IQUAD
2002   FORMAT('  IWIND   ',I4,' IWCAP   ',I4,' IQUAD   ',I4)
       WRITE(PRINTF,2003) ITRIAD, IBOT , ISURF
2003   FORMAT('  ITRIAD  ',I4,' IBOT    ',I4,' ISURF   ',I4)
       WRITE(PRINTF,*)' -----------------------------------------',     &
                      '----------------------'
     ENDIF

!    --- calculate diffraction parameter and its derivatives           
!JQI     IF(IDIFFR > 0)                                                    &
!JQI          CALL DIFPAR( AC2   , SPCSIG, KGRPNT, COMPDA(1,JDP2),         &
!JQI	               CROSS , XCGRID, YCGRID, XYTST  )                 

!
!            *** START ITERATION PROCESS WITH 4 SWEEPS ***
!
!    *** loop over sweep directions ***
!
       CALL ACTION_ALLO

       !if(msr)print*,'begin ***************'
       DO JDUM = 1,M
         DEP_TMP = COMPDA(JDUM,JDP2)
	 IF(DEP_TMP > DEPMIN)THEN      
         CALL SWOMPU1(DTW      ,SNLC1    ,DAL1     ,DAL2     ,       &
	              DAL3     ,XIS      ,SWTSDA   ,INOCNT   ,       &
		      ITER     ,CGO      ,CAX      ,CAY      ,       &
		      CAS      ,CAD      ,SWMATR   ,KWAVE    ,       &
		      ALIMW    ,GROWW    ,UE       ,SA1      ,       &
		      SA2      ,DA1C     ,DA1P     ,DA1M     ,       &
		      DA2C     ,DA2P     ,DA2M     ,SFNL     ,       &
		      DSNL     ,MEMNL4   ,IDCMIN   ,IDCMAX   ,       &
		      WWINT    ,WWAWG    ,WWSWG    ,ISCMIN   ,       &
		      ISCMAX   ,AC1      ,IT       ,CROSS    ,       &
		      OBREDF   ,REFLSO   ,ISLMIN   ,NFLIM    ,       &
		      NRSCAL   ,CAX1     ,CAY1     ,JDUM     ,       &
		      IDTOT    )
	 END IF	      
       END DO  ! JDUM

!       print*,"fshape=",fshape,"dshape=",dshape,"pshape(2)=",pshape(2), pi2*0.01
!       stop
       DO IOB = 1, IOBCN_W
         IOB_NODE = I_OBC_N_W(IOB)
         CALL SSHAPE(N32(1,1,IOB_NODE),SPCSIG,SPCDIR,FSHAPE,DSHAPE)
       END DO
 
#      if defined (MULTIPROCESSOR)
       IF(NESTING_ON_WAVE)THEN
        IF(NESTING_TYPE_WAVE == 'wave parameters')THEN
  
         CALL SET_VAR_WAVE(WaveTime,HSC1=HSC1)
         CALL SET_VAR_WAVE(WaveTime,TPEAK=TPEAK)
         CALL SET_VAR_WAVE(WaveTime,DIRDEG1=DIRDEG1)
       
         DO IOB = 1, IOBCN_W
          IOB_NODE = I_OBC_N_W(IOB)
	  SPPARM(1) = HSC1(IOB_NODE)
	  SPPARM(2) = TPEAK(IOB_NODE)
	  SPPARM(3) = DIRDEG1(IOB_NODE)
	  SPPARM(4) = 20.0_SP 
          CALL SSHAPE(N32(1,1,IOB_NODE),SPCSIG,SPCDIR,FSHAPE,DSHAPE)
         END DO
 
        ELSE IF(NESTING_TYPE_WAVE == 'spectral density')THEN
  
         ALLOCATE(N32_TTMP(MDC,MSC,0:MT));  N32_TTMP = 0.0_SP
         CALL SET_VAR_WAVE_AC2(WaveTime,AC2=N32_TTMP)
	 
	 DO IOB = 1, IOBCN_W
	   IOB_NODE = I_OBC_N_W(IOB)
	   N32(:,:,IOB_NODE) = N32_TTMP(:,:,IOB_NODE)
	 END DO
	 DEALLOCATE(N32_TTMP)  
  
        ELSE
         CALL FATAL_ERROR("THE PARAMETER NESTING_TYPE_WAVE SHOULD BE 'wave parameters' or 'spectral density'")
        END IF
       END IF
#      endif


#      if defined (MULTIPROCESSOR)
       IF(PAR)THEN
         ALLOCATE(N32_TMP(0:MT));  N32_TMP = 0.0
         DO ISC = 1,MSC
           DO IDC = 1,MDC
!             DO IP = 1,MT
             DO IP = 1,M
               N32_TMP(IP)=N32(IDC,ISC,IP)
	     END DO  
             CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,1,MYID,NPROCS,N32_TMP)
!             CALL EXCHANGE(NC,MT,1,MYID,NPROCS,N32_TMP) 
             CALL AEXCHANGE(NC,MYID,NPROCS,N32_TMP)  
	     DO IP = 1,MT
	       N32(IDC,ISC,IP) = N32_TMP(IP)
	     END DO  
           END DO
         END DO
       DEALLOCATE(N32_TMP)
       END IF    

!JQI       IF(PAR)CALL MPI_BARRIER(MPI_FVCOM_GROUP,IERR)
#      endif

!       print*,"U10=",U10,"WDIC=",WDIC
     !  if(msr)print*,'between ********'
       CALL ADV_N(DTW)

#  if defined (MULTIPROCESSOR)
#  if defined (ICE)
      CALL ICETT2
      DEALLOCATE(Sice)
#  endif
#  endif

!JQI#      if defined (MULTIPROCESSOR)
!JQI       IF(PAR)THEN
!JQI         ALLOCATE(AC2_TMP(0:MT));   AC2_TMP = 0.0
!JQI         DO ISC = 1,MSC
!JQI           DO IDC = 1,MDC
!JQI!             DO IP = 1,MT
!JQI             DO IP = 1,M
!JQI	       AC2_TMP(IP)=AC2(IDC,ISC,IP)
!JQI	     END DO  
!JQI             CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,1,MYID,NPROCS,AC2_TMP)
!JQI!             CALL EXCHANGE(NC,MT,1,MYID,NPROCS,AC2_TMP)
!JQI            CALL AEXCHANGE(NC,MYID,NPROCS,AC2_TMP)   
!JQI             DO IP = 1,MT
!JQI	       AC2(IDC,ISC,IP) = AC2_TMP(IP)
!JQI	     END DO
!JQI          END DO
!JQI         END DO
!JQI         DEALLOCATE(AC2_TMP)
!JQI       END IF    	  

!JQI       IF(PAR)CALL MPI_BARRIER(MPI_FVCOM_GROUP,IERR)
!JQI#      endif
      ! if(msr)print*,'after ********'

# if defined (MULTIPROCESSOR)
  IF(PAR)THEN
    ALLOCATE(COMPDA_TMP1(0:MT));   COMPDA_TMP1 = 0.0
    DO JDUM = 1,M
      COMPDA_TMP1(JDUM)=COMPDA(JDUM,JDP2)
    END DO  
!    CALL EXCHANGE(NC,MT,1,MYID,NPROCS,COMPDA_TMP1)
    CALL AEXCHANGE(NC,MYID,NPROCS,COMPDA_TMP1)   
    DO JDUM = 1,MT
      COMPDA(JDUM,JDP2) = COMPDA_TMP1(JDUM)
    END DO
    DEALLOCATE(COMPDA_TMP1)
  END IF    	  
# endif
       
       ALLOCATE(PASS_TMP(0:MT)); PASS_TMP = 1.0   !.TRUE.      !1
       DO JDUM = 1, M
         IF(ISONB_W(JDUM) == 1) THEN
!          compute absolute wind velocity                                    
           IF(VARWI)THEN
             AWX_W = COMPDA(JDUM,JWX2)    !WX2(IG)
             AWY_W = COMPDA(JDUM,JWY2)    !WY2(IG)
!	     write(100+myid,*) awx,awy
           ELSE
             AWX_W = U10 * COS(WDIC)
             AWY_W = U10 * SIN(WDIC)
           END IF
!          compute relative wind velocity                                      
           IF(ICUR == 0)THEN
             RWX = AWX_W
             RWY = AWY_W
           ELSE
             RWX = AWX_W - COMPDA(1,JVX2)  !UX2(IG)
             RWY = AWY_W - COMPDA(1,JVY2)  !UY2(IG)
           END IF
!          compute absolute value of relative wind velocity                   
           WIND10 = SQRT(RWX**2+RWY**2)
           IF(WIND10 > 0.)THEN
             THETAW = ATAN2 (RWY,RWX)
             THETAW = MOD ( (THETAW + PI2_W) , PI2_W )
           ELSE
             THETAW = 0.
           END IF
	 ALPHA3 = THETAW    !3.1415926*3.0/2.0    !0.0      !3.1415926/2.0 
	 ALPHA3 = ALPHA3*180.0/3.1415926
!	 write(100+myid,*) alpha3
!         if(jdum == 1004)then
!	   print*,"ALPHA ",JDUM,ALPHA(JDUM),ALPHA3
!	 end if  
!LWU
!JQI  	   IF(ABS(ALPHA(JDUM)-ALPHA3) <= 90.0 .OR.    &
!JQI	      ABS(ALPHA(JDUM)-360.0-ALPHA3) <= 90.0)THEN
!JQI	     PASS_TMP(JDUM) = 0.0     !.FALSE.    !0
	     
!	     write(100+myid,*) "PASS=",PASS(JDUM),NGID(JDUM)
!JQI  	   END IF  
         END IF
       END DO
!       write(100+myid,*)
#      if defined (MULTIPROCESSOR)
!       IF(PAR)CALL EXCHANGE(NC,MT,1,MYID,NPROCS,PASS_TMP)
       IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,PASS_TMP)     
#      endif       
       ALLOCATE(PASS(MT))
       PASS(1:MT) = INT(PASS_TMP(1:MT))
       DEALLOCATE(PASS_TMP)
       
       DO JDUM = 1,M   
         DEP_TMP = COMPDA(JDUM,JDP2)
	 IF(DEP_TMP > DEPMIN .AND. PASS(JDUM) == 1)THEN
           CALL SWOMPU2(DTW      ,SNLC1    ,DAL1     ,DAL2     ,        &
	                DAL3     ,XIS      ,SWTSDA   ,INOCNT   ,        &
		        ITER     ,CGO      ,CAX      ,CAY      ,        &
		        CAS      ,CAD      ,SWMATR   ,KWAVE    ,        &
		        ALIMW    ,GROWW    ,UE       ,SA1      ,        &
		        SA2      ,DA1C     ,DA1P     ,DA1M     ,        &
		        DA2C     ,DA2P     ,DA2M     ,SFNL     ,        &
		        DSNL     ,MEMNL4   ,IDCMIN   ,IDCMAX   ,        &
		        WWINT    ,WWAWG    ,WWSWG    ,ISCMIN   ,        &
		        ISCMAX   ,AC1      ,IT       ,CROSS    ,        &
		        OBREDF   ,REFLSO   ,ISLMIN   ,NFLIM    ,        &
		        NRSCAL   ,CAX1     ,CAY1     ,JDUM     ,        &
		        IDTOT    )

         ENDIF		      
       END DO  ! JDUM
       DEALLOCATE(PASS)
       
       DO IOB = 1, IOBCN_W
         IOB_NODE = I_OBC_N_W(IOB)
         CALL SSHAPE(AC2(1,1,IOB_NODE),SPCSIG,SPCDIR,FSHAPE,DSHAPE)
       END DO

#      if defined (MULTIPROCESSOR)
       IF(NESTING_ON_WAVE)THEN
        IF(NESTING_TYPE_WAVE == 'wave parameters')THEN
  
         CALL SET_VAR_WAVE(WaveTime,HSC1=HSC1)
         CALL SET_VAR_WAVE(WaveTime,TPEAK=TPEAK)
         CALL SET_VAR_WAVE(WaveTime,DIRDEG1=DIRDEG1)
       
         DO IOB = 1, IOBCN_W
          IOB_NODE = I_OBC_N_W(IOB)
	  SPPARM(1) = HSC1(IOB_NODE)
	  SPPARM(2) = TPEAK(IOB_NODE)
	  SPPARM(3) = DIRDEG1(IOB_NODE)
	  SPPARM(4) = 20.0_SP 

          CALL SSHAPE(AC2(1,1,IOB_NODE),SPCSIG,SPCDIR,FSHAPE,DSHAPE)
         END DO

        ELSE IF(NESTING_TYPE_WAVE == 'spectral density')THEN
  
         ALLOCATE(AC2_TTMP(MDC,MSC,0:MT));  AC2_TTMP = 0.0_SP
         CALL SET_VAR_WAVE_AC2(WaveTime,AC2=AC2_TTMP)
	 
         DO IOB = 1, IOBCN_W
           IOB_NODE = I_OBC_N_W(IOB)
           AC2(:,:,IOB_NODE) = AC2_TTMP(:,:,IOB_NODE)
         END DO
         DEALLOCATE(AC2_TTMP)  
  
        ELSE
         CALL FATAL_ERROR("THE PARAMETER NESTING_TYPE_WAVE SHOULD BE 'wave parameters' or 'spectral density'")
        END IF
       END IF
#      endif       

# if defined(PLBC)
       DO ISC = 1,MSC
           DO IDC = 1,MDC
           CALL replace_ac2(IDC,ISC)
          ENDDO
       ENDDO
# endif

#      if defined (MULTIPROCESSOR)
       IF(PAR)THEN
         ALLOCATE(AC2_TMP(0:MT));   AC2_TMP = 0.0
         DO ISC = 1,MSC
           DO IDC = 1,MDC
!             DO IP = 1,MT
             DO IP = 1,M
	       AC2_TMP(IP)=AC2(IDC,ISC,IP)
	     END DO  
             CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,1,MYID,NPROCS,AC2_TMP)
 !            CALL EXCHANGE(NC,MT,1,MYID,NPROCS,AC2_TMP)
             CALL AEXCHANGE(NC,MYID,NPROCS,AC2_TMP)  
             DO IP = 1,MT
	       AC2(IDC,ISC,IP) = AC2_TMP(IP)
	     END DO
           END DO
         END DO
         DEALLOCATE(AC2_TMP)
       END IF    	  

!JQI       IF(PAR)CALL MPI_BARRIER(MPI_FVCOM_GROUP,IERR)
#      endif

       CALL ACTION_DEALLO
       
!----------------------------------------------------------------------   
!     Each thread sum contributions to the global INOCNV counter          
!     which counts the number of grid points over the four sweeps         
!     in which the SIP solver did not converge.                           
!----------------------------------------------------------------------   
     INOCNV = INOCNV + INOCNT                                            

!
!    *** compute wave-induced setup ***                                
!
#    if defined(WAVE_SETUP)
     IF(LSETUP > 0)THEN 
       print*,"LSETUP= *",LSETUP
       CALL WAVE_INDUCED_FORCE(COMPDA(1,JDP2),COMPDA(1,JDPSAV))
       print*,"after calling WAVE_INDUCED_FORCE"
       CALL WAVE_INDUCED_SETUP(COMPDA(1,JDP2),COMPDA(1,JDPSAV))  
       print*,"after calling WAVE_INDUCED_SETUP"
     ENDIF  
#    endif                                                                
!
!----------------------------------------------------------------------   
!     End master thread region.                                           
!----------------------------------------------------------------------   
!
!    *** check if numerical accuracy has been reached       ***
!    *** this is done in parallel within OpenMP environment ***        
!
     IF(PNUMS(21) == 0.)THEN                                         
!       CALL SACCUR (COMPDA(1,JDP2),KGRPNT        ,XYTST         ,       &
!                    AC2           ,SPCSIG        ,ACCUR         ,       &
!		     HSAC1         ,HSAC2         ,SACC1         ,       &
!		     SACC2         ,COMPDA(1,JDHS),COMPDA(1,JDTM),       &
!		     I1MYC         ,I2MYC         )    
     ELSE IF(PNUMS(21) == 1.)THEN                                    
!       CALL SWSTPC (HSAC0         ,HSAC1         ,HSAC2         ,       &
!                    SACC1         ,SACC2         ,HSDIFC        ,       &
!		     COMPDA(1,JDHS),COMPDA(1,JDTM),COMPDA(1,JDP2),       &
!                    ACCUR         ,I1MYC         ,I2MYC         )                           
     END IF                                                           

450 CONTINUE                                                            
470 CONTINUE

!---------------------------------------------------------------------- 
!  Begin deallocate private arrays.                                    
!----------------------------------------------------------------------
   DEALLOCATE(IDCMIN)                                                  
   DEALLOCATE(IDCMAX)                                                  
   DEALLOCATE(ISCMIN)                                                  
   DEALLOCATE(ISCMAX)                                                  
   DEALLOCATE(CGO)                                                     
   DEALLOCATE(KWAVE)                                                   
   DEALLOCATE(CAX)                                                     
   DEALLOCATE(CAY)                                                     
   DEALLOCATE(CAS)                                                     
   DEALLOCATE(CAD)                                                     
   DEALLOCATE(CAX1)                                                    
   DEALLOCATE(CAY1)                                                    
   DEALLOCATE(ALIMW)                                                   
   DEALLOCATE(UE)                                                      
   DEALLOCATE(SA1)                                                     
   DEALLOCATE(SA2)                                                     
   DEALLOCATE(SFNL)                                                    
   DEALLOCATE(DA1C)                                                    
   DEALLOCATE(DA1P)                                                    
   DEALLOCATE(DA1M)                                                    
   DEALLOCATE(DA2C)                                                    
   DEALLOCATE(DA2P)                                                    
   DEALLOCATE(DA2M)                                                    
   DEALLOCATE(DSNL)                                                    
   DEALLOCATE(OBREDF)                                                  
   DEALLOCATE(REFLSO)                                                  
   DEALLOCATE(GROWW)                                                   
   DEALLOCATE(SWMATR)                                                  
!----------------------------------------------------------------------   
!  End deallocate private arrays.                                      
!----------------------------------------------------------------------

!----------------------------------------------------------------------   
!  Begin deallocate shared arrays.                                     
!----------------------------------------------------------------------   
   DEALLOCATE(SETPDA)                                                  
   DEALLOCATE(HSAC1)                                                   
   DEALLOCATE(HSAC2)                                                   
   DEALLOCATE(SACC1)                                                   
   DEALLOCATE(SACC2)                                                   
   DEALLOCATE(HSAC0)                                                   
   DEALLOCATE(HSDIFC)                                                  
   DEALLOCATE(ISLMIN)                                                  
   DEALLOCATE(NFLIM)                                                   
   DEALLOCATE(NRSCAL)                                                  
   DEALLOCATE(MEMNL4)                                                  
   DEALLOCATE(SWTSDA)                                                  
!----------------------------------------------------------------------
!  End deallocate shared arrays.                                       
!----------------------------------------------------------------------
!
   RETURN
   END SUBROUTINE SWCOMP
!
!************************************************************************
   SUBROUTINE SWOMPU1 (DTW      ,SNLC1    ,DAL1     ,DAL2     ,         &
                       DAL3     ,XIS      ,SWTSDA   ,INOCNV   ,         &
		       ITER     ,CGO      ,CAX      ,CAY      ,         &
		       CAS      ,CAD      ,SWMATR   ,KWAVE    ,         &
		       ALIMW    ,GROWW    ,UE       ,SA1      ,         &
		       SA2      ,DA1C     ,DA1P     ,DA1M     ,         &
		       DA2C     ,DA2P     ,DA2M     ,SFNL     ,         &
		       DSNL     ,MEMNL4   ,IDCMIN   ,IDCMAX   ,         &
		       WWINT    ,WWAWG    ,WWSWG    ,ISCMIN   ,         &
		       ISCMAX   ,AC1      ,IT       ,CROSS    ,         &
		       OBREDF   ,REFLSO   ,ISLMIN   ,NFLIM    ,         &
		       NRSCAL   ,CAX1     ,CAY1     ,IG       ,         &
		       IDTOT    )
!
!************************************************************************
!
!     This subroutine computes the wave spectrum for one sweep
!     direction, and is called four times per iteration.
!
!************************************************************************
!
   USE OCPCOMM1                                                        
   USE OCPCOMM2                                                        
   USE OCPCOMM3                                                        
   USE OCPCOMM4 
   USE M_GENARR                                                       
   USE SWCOMM1                                                         
   USE SWCOMM2                                                         
   USE SWCOMM3                                                         
   USE SWCOMM4                                                         
#  if defined (EXPLICIT)
   USE MOD_ACTION_EX
#  else   
   USE MOD_ACTION_IM
#  endif   
!   USE ALL_VARS, ONLY : M,MT,NTVE,AC2,COMPDA                                                   
   USE VARS_WAVE, ONLY : M,MT,NTVE,AC2,COMPDA                                                   
#  if defined (SPHERICAL)
   USE MOD_SPHERICAL
#  endif      

   IMPLICIT NONE
   
   INTEGER ITER, IT, IDC, ISC, IG

   INTEGER  IC    ,IS    ,                    &
	    IDWMIN,IDWMAX,IDTOT ,ISTOT ,IDDLOW,IDDTOP,                &
	    ISSTOP,INOCNV                                        
   INTEGER  LINK(MICMAX)                                           
!
   REAL     DTW    ,                                     &
            ETOT  ,AC2TOT,ABRBOT,HM    ,HS    ,QBLOC ,                &
	    SMESPC,KMESPC,ETOTW ,WIND10,FPM   ,                       &
	    THETAW,SNLC1 ,DAL1  ,DAL2  ,DAL3  ,XIS   ,                &
	    UFRIC ,SMEBRK,SZEROC,EPS2WC,                              &
	    DISWCP,WCPSME,WCPKME,WCPQB ,WCPHM
!
   LOGICAL  INSIDE                                                     
!
   INTEGER :: IDCMIN(MSC)                              
   INTEGER :: IDCMAX(MSC)    ,ISCMIN(MDC)    ,ISCMAX(MDC)
   INTEGER :: WWINT(*)
   INTEGER :: CROSS(2,M) 
!
!  *** number of arrays for SWAN ***
!
   REAL  :: AC1(MDC,MSC,0:MT)                                         
   REAL  :: CGO(MSC,MICMAX)          ,                                 &
            CAX(MDC,MSC,MICMAX)      ,                                 &
	    CAY(MDC,MSC,MICMAX)      ,                                 &
	    CAX1(MDC,MSC,MICMAX)     ,                                 &
	    CAY1(MDC,MSC,MICMAX)     ,                                 &
	    CAS(MDC,MSC,MICMAX)      ,                                 &
	    CAD(MDC,MSC,MICMAX)                                        
   REAL  :: ALIMW(MDC,MSC)
   REAL  :: SWMATR(MDC,MSC,MSWMATR)                                    
   REAL  :: KWAVE(MSC,MICMAX),                                         &
            UE(MSC4MI:MSC4MA , MDC4MI:MDC4MA )    ,                    &
	    SA1(MSC4MI:MSC4MA , MDC4MI:MDC4MA )   ,                    &
	    SA2(MSC4MI:MSC4MA , MDC4MI:MDC4MA )   ,                    &
	    DA1C(MSC4MI:MSC4MA , MDC4MI:MDC4MA )  ,                    &
	    DA1P(MSC4MI:MSC4MA , MDC4MI:MDC4MA )  ,                    &
	    DA1M(MSC4MI:MSC4MA , MDC4MI:MDC4MA )  ,                    &
	    DA2C(MSC4MI:MSC4MA , MDC4MI:MDC4MA )  ,                    &
	    DA2P(MSC4MI:MSC4MA , MDC4MI:MDC4MA )  ,                    &
	    DA2M(MSC4MI:MSC4MA , MDC4MI:MDC4MA )  ,                    &
	    SFNL(MSC4MI:MSC4MA , MDC4MI:MDC4MA )  ,                    &
	    DSNL(MSC4MI:MSC4MA , MDC4MI:MDC4MA )  ,                    &
	    MEMNL4(MDC,MSC,M)               ,                      &
	    SWTSDA(MDC,MSC,NPTSTA,MTSVAR)         ,                    &
	    WWAWG(*)                              ,                    &
	    WWSWG(*)                              ,                    &
	    RDX(10)       ,RDY(10)               ,                     &
	    OBREDF(MDC,MSC,2)                    ,                     &
	    REFLSO(MDC,MSC)                                       
!
   LOGICAL  GROWW(MDC,MSC)
!
   INTEGER :: ISLMIN(M), NFLIM(M), NRSCAL(M)               
!
!  *** Get grid point numbers for points in computational stencil ***

   IF(IG >= 1)THEN                                             
     ICMAX = NTVE(IG)+1
   ENDIF

!  *** Compute wavenumber KWAVE and group velocity CGO   ***
!  *** in the gridpoints of the stencil                  ***

   CALL SWAPAR (IG, ICMAX, COMPDA(1,JDP2), KWAVE, CGO) 
!
!  *** compute minimum and maximum counter (IDCMIN and ***
!  *** IDCMAX) ***
!
   CALL SWPSEL (IDCMIN         ,IDCMAX         ,CAX            ,  &
                CAY            ,ISCMIN         ,ISCMAX         ,  &
		IDTOT          ,ISTOT          ,IDDLOW         ,  &
		IDDTOP         ,ISSTOP         ,COMPDA(1,JDP2) ,  &
		COMPDA(1,JVX2) ,COMPDA(1,JVY2) ,SPCDIR         )

!  *** compute the propagation velocities CAS and CAD   ***
!  *** for the central gridpoint only and only for the  ***
!  *** directional domain : IDCMIN-1 until IDCMAX+1     ***

   CALL SPROSD (KWAVE          ,CAS            ,CAD            ,    &
                CGO            ,COMPDA(1,JDP2) ,COMPDA(1,JDP1) ,    &
		SPCDIR(1,2)    ,SPCDIR(1,3)    ,COMPDA(1,JVX2) ,    &
		COMPDA(1,JVY2) ,IDCMIN         ,IDCMAX         ,    &
		SPCDIR(1,4)    ,SPCDIR(1,6)    ,SPCDIR(1,5)    ,    &
		RDX            ,RDY            ,CAX            ,    &
		CAY            ,IG             )

#  if defined (SPHERICAL)
   IF(KSPHER > 0 .AND. IREFR /= 0)THEN                              
!
!     *** compute the change of propagation velocity CAD   ***
!     *** due to the use of spherical coordinates          ***
!
     CALL DSPHER (CAD, CAX, CAY, IG, SPCDIR(1,2), SPCDIR(1,3))             
   ENDIF
#  endif   
!
   IF(IDTOT > 0)THEN                                              
!
!    *** initialize friction velocity and Fpm frequency even ***
!    *** when there is no wind input                         ***
!
     UFRIC = 1.E-15
     FPM   = 1.E-15

!
!    *** If there are obstacles crossing the points in the stencil ***
!    *** then the transmission and reflection coeff. are computed  ***
!    *** and also the contribution to the source term              ***
!
     IF(NUMOBS /= 0)THEN
!
!      *** OBREDF(:,:,2) are the transmission coeff for the two links ***
!      *** in the stencil (between the three point on the stencil)    ***
!      *** REFLSO(:,:) contains the contribution to the source term   ***
!
       DO IDC = 1, MDC
         DO ISC = 1, MSC
           OBREDF(IDC,ISC,1) = 1.                                      
           OBREDF(IDC,ISC,2) = 1.                                      
           REFLSO(IDC,ISC  ) = 0.
         ENDDO
       ENDDO
!
         LINK(1) = CROSS(1,KCGRD(1))
         LINK(2) = CROSS(2,KCGRD(1))

       IF(LINK(1) /= 0 .OR. LINK(2) /= 0)THEN                    
!
!JQI         CALL SWTRCF (COMPDA(1,JWLV2),                                 &
!JQI	              COMPDA(1,JHS), LINK, OBREDF,                     &
!JQI		      AC2, REFLSO, KGRPNT, XCGRID,                     &
!JQI		      YCGRID, CAX, CAY, RDX, RDY, LSWMAT(1,1,JABIN),   &
!JQI		      SPCSIG)                                         
       ENDIF
!
     ENDIF
!
     IF((DYNDEP .OR. ICUR == 1) .AND. ITFRE /= 0)THEN               
       CALL ALGORITHM_FCT(CAS,IG,DTW,IDCMIN,IDCMAX)
     ELSE
       DO ISC = 1,MSC
	 DO IDC = 1,MDC
           N31(IDC,ISC) = AC2(IDC,ISC,IG)
	 END DO
       END DO
     END IF


     IF(IREFR /= 0)THEN                                              
       IF(PROPFL == 0)THEN    
         CALL ALGORITHM_CRANK_NICOLSON(CAD,IG,DTW,IDCMIN,IDCMAX,DDIR)
       END IF
     END IF
    ENDIF

   RETURN
   END SUBROUTINE SWOMPU1
 
!
!************************************************************************
!
   SUBROUTINE SWOMPU2 (DTW      ,SNLC1    ,DAL1     ,DAL2     ,         &
                       DAL3     ,XIS      ,SWTSDA   ,INOCNV   ,         &
		       ITER     ,CGO      ,CAX      ,CAY      ,         &
		       CAS      ,CAD      ,SWMATR   ,KWAVE    ,         &
		       ALIMW    ,GROWW    ,UE       ,SA1      ,         &
		       SA2      ,DA1C     ,DA1P     ,DA1M     ,         &
		       DA2C     ,DA2P     ,DA2M     ,SFNL     ,         &
		       DSNL     ,MEMNL4   ,IDCMIN   ,IDCMAX   ,         &
		       WWINT    ,WWAWG    ,WWSWG    ,ISCMIN   ,         &
		       ISCMAX   ,AC1      ,IT       ,CROSS    ,         &
		       OBREDF   ,REFLSO   ,ISLMIN   ,NFLIM    ,         &
		       NRSCAL   ,CAX1     ,CAY1     ,IG       ,         &
		       IDTOT    )
!
!************************************************************************

   USE OCPCOMM1                                                        
   USE OCPCOMM2                                                        
   USE OCPCOMM3                                                        
   USE OCPCOMM4 
   USE M_GENARR                                                       
   USE SWCOMM1                                                         
   USE SWCOMM2                                                         
   USE SWCOMM3                                                         
   USE SWCOMM4                                                         
#  if defined (EXPLICIT)
   USE MOD_ACTION_EX
#  else    
   USE MOD_ACTION_IM
#  endif   
#  if defined (MULTIPROCESSOR)
   USE MOD_PAR
#  endif
!   USE ALL_VARS, ONLY : M,MT,NTVE,SERIAL,PAR,MYID,IOBCN,I_OBC_N,AC2,COMPDA
   USE VARS_WAVE, ONLY : M,MT,NTVE,SERIAL,PAR,MYID,IOBCN_W,I_OBC_N_W,AC2,COMPDA,   &
                         SOURCE_DTMAX,SOURCE_DTMIN
   
   IMPLICIT NONE                                                

   INTEGER ITER,   IT
   LOGICAL STPNOW                                                      

   INTEGER  IC    ,IS    ,                       &
            IG    ,IDC   ,ISC   ,                       &
	    IDWMIN,IDWMAX,IDTOT ,ISTOT ,IDDLOW,IDDTOP,                &
	    ISSTOP,INOCNV                                        
   INTEGER  LINK(MICMAX)                                           
!
   REAL     DTW    ,                                     &
            ETOT  ,AC2TOT,ABRBOT,HM    ,HS    ,QBLOC ,                &
	    SMESPC,KMESPC,ETOTW ,WIND10,FPM   ,                       &
	    THETAW,SNLC1 ,DAL1  ,DAL2  ,DAL3  ,XIS   ,                &
	    UFRIC ,SMEBRK,SZEROC,EPS2WC,                              &
	    DISWCP,WCPSME,WCPKME,WCPQB ,WCPHM

   REAL     DTDYN,DTTOT,NSTEPS_WAVE
   
   
   REAL :: XP,TPI2,FACP,DAMAX,XR,AFILT,DTS,DTMAX,XFILT,AMAX,AFAC,OFFSET
   REAL :: DTMIN,HDT
   INTEGER :: IDT,ISMAX,ISP1
   LOGICAL :: SHAVE
   REAL :: DAM(MDC,MSC)   
!
   LOGICAL  INSIDE                                                     
!
   INTEGER :: IDCMIN(MSC)                              
   INTEGER :: IDCMAX(MSC)    ,ISCMIN(MDC)    ,ISCMAX(MDC)
   INTEGER :: WWINT(*)
   INTEGER :: CROSS(2,M)
!
!  *** number of arrays for SWAN ***
!
   REAL  :: AC1(MDC,MSC,0:MT)                                         
   REAL  :: CGO(MSC,MICMAX)          ,                                 &
            CAX(MDC,MSC,MICMAX)      ,                                 &
	    CAY(MDC,MSC,MICMAX)      ,                                 &
	    CAX1(MDC,MSC,MICMAX)     ,                                 &
	    CAY1(MDC,MSC,MICMAX)     ,                                 &
	    CAS(MDC,MSC,MICMAX)      ,                                 &
	    CAD(MDC,MSC,MICMAX)                                        
   REAL  :: ALIMW(MDC,MSC)
   REAL  :: SWMATR(MDC,MSC,MSWMATR)                                    
   REAL  :: KWAVE(MSC,MICMAX),                                         &
            UE(MSC4MI:MSC4MA , MDC4MI:MDC4MA )    ,                    &
	    SA1(MSC4MI:MSC4MA , MDC4MI:MDC4MA )   ,                    &
	    SA2(MSC4MI:MSC4MA , MDC4MI:MDC4MA )   ,                    &
	    DA1C(MSC4MI:MSC4MA , MDC4MI:MDC4MA )  ,                    &
	    DA1P(MSC4MI:MSC4MA , MDC4MI:MDC4MA )  ,                    &
	    DA1M(MSC4MI:MSC4MA , MDC4MI:MDC4MA )  ,                    &
	    DA2C(MSC4MI:MSC4MA , MDC4MI:MDC4MA )  ,                    &
	    DA2P(MSC4MI:MSC4MA , MDC4MI:MDC4MA )  ,                    &
	    DA2M(MSC4MI:MSC4MA , MDC4MI:MDC4MA )  ,                    &
	    SFNL(MSC4MI:MSC4MA , MDC4MI:MDC4MA )  ,                    &
	    DSNL(MSC4MI:MSC4MA , MDC4MI:MDC4MA )  ,                    &
	    MEMNL4(MDC,MSC,MT)               ,                      &
	    SWTSDA(MDC,MSC,NPTSTA,MTSVAR)         ,                    &
	    WWAWG(*)                              ,                    &
	    WWSWG(*)                              ,                    &
	    RDX(10)       ,RDY(10)               ,                     &
	    OBREDF(MDC,MSC,2)                    ,                     &
	    REFLSO(MDC,MSC)                                       
!
   LOGICAL  GROWW(MDC,MSC) 
!
   INTEGER :: ISLMIN(M), NFLIM(M), NRSCAL(M) 
   INTEGER :: IDDUM   
   INTEGER :: IOB,IOB_NODE           
!
!  *** Get grid point numbers for points in computational stencil ***

   IF(IG >= 1)THEN                                             
     ICMAX = NTVE(IG)+1
   ENDIF
!
!  *** If there are obstacles crossing the points in the stencil ***
!  *** then fall back to first order scheme                      ***
!

!
!     *** determine whether the point is a test point ***
!
   IPTST  = 0
   TESTFL = .FALSE.

!  *** Compute wavenumber KWAVE and group velocity CGO   ***
!  *** in the gridpoints of the stencil                  ***

   CALL SWAPAR(IG, ICMAX, COMPDA(1,JDP2), KWAVE, CGO) 
!
!  *** compute minimum and maximum counter (IDCMIN and ***
!  *** IDCMAX)   ***
!
   CALL SWPSEL (IDCMIN         ,IDCMAX         ,CAX            ,  &
                CAY            ,ISCMIN         ,ISCMAX         ,  &
		IDTOT          ,ISTOT          ,IDDLOW         ,  &
		IDDTOP         ,ISSTOP         ,COMPDA(1,JDP2) ,  &
		COMPDA(1,JVX2) ,COMPDA(1,JVY2) ,SPCDIR         )

   IF(IDTOT > 0)THEN                                              
!
!    *** initialize friction velocity and Fpm frequency even ***
!    *** when there is no wind input                         ***
!
     UFRIC = 1.E-15
     FPM   = 1.E-15
     IF(IWIND >= 1)THEN
!
!      *** compute the wind speed, mean wind direction, the    ***
!      *** PM frequency, wind friction velocity U*  and the    ***
!      *** minimum and maximum counters for active wind input  ***
!
       CALL WINDP1 (WIND10         ,THETAW         ,IDWMIN         , &
                    IDWMAX         ,FPM            ,UFRIC          , &
		    COMPDA(1,JWX2) ,COMPDA(1,JWY2) ,SPCDIR         , &
		    COMPDA(1,JVX2) ,COMPDA(1,JVY2) ,SPCSIG         , &
		    IG   )  
     END IF
!
!    *** estimate action density in case of first iteration ***
!    *** in stationary mode (since it is zero in first      ***
!        stationary run)                                    ***
!
     IF(ICOND /= 4 .AND. ITER == 1 .AND. NSTATC == 0)THEN            
       CALL SPREDT (AC2     ,CAX     ,CAY     ,IDCMIN  ,IDCMAX  ,   &
		    ISSTOP  ,RDX     ,RDY     )
     END IF

!print*, 'after SPREDT'
     DTDYN  = 0.0
     DTTOT  = 0.0
     NSTEPS_WAVE = 0.0
     
     XP = 0.15
     TPI2 = 2.0*PI_W
     FACP   = XP / PI_W * 0.62E-3 * TPI2**4 / GRAV_W**2
     DO ISC = 1,MSC
       DO IDC = 1,MDC
         DAM(IDC,ISC) = FACP/(SPCSIG(ISC)*KWAVE(ISC,1)**3)
       END DO
     END DO 	 
     
     DO
       NSTEPS_WAVE = NSTEPS_WAVE + 1

!
!       Calculate various integral parameters for use in the source terms
!
     CALL SINTGRL (SPCDIR          ,KWAVE   ,COMPDA(1,JDP2)   ,QBLOC   , &
                   COMPDA(1,JURSEL),RDX     ,RDY              ,AC2TOT  , &
		   ETOT            ,ABRBOT  ,COMPDA(1,JUBOT)  ,HS      , &
		   COMPDA(1,JQB)   ,HM      ,KMESPC           ,SMEBRK  , &
		   COMPDA(1,JPBOT) ,IG      )    
    !print*,'after SINTGRL'                  
!
     COMPDA(IG,JHS) = HS                                         
!
!    *** If there are obstacles crossing the points in the stencil ***
!    *** then the transmission and reflection coeff. are computed  ***
!    *** and also the contribution to the source term              ***
!
     IF(NUMOBS /= 0)THEN
!
!      *** OBREDF(:,:,2) are the transmission coeff for the two links ***
!      *** in the stencil (between the three point on the stencil)    ***
!      *** REFLSO(:,:) contains the contribution to the source term   ***
!
       DO IDC = 1, MDC
         DO ISC = 1, MSC
           OBREDF(IDC,ISC,1) = 1.                                      
           OBREDF(IDC,ISC,2) = 1.                                      
           REFLSO(IDC,ISC  ) = 0.
         ENDDO
       ENDDO
!
         LINK(1) = CROSS(1,KCGRD(1))
         LINK(2) = CROSS(2,KCGRD(1))

       IF(LINK(1) /= 0 .OR. LINK(2) /= 0)THEN                    
!
!JQI         CALL SWTRCF (COMPDA(1,JWLV2),                                 &
!JQI	              COMPDA(1,JHS), LINK, OBREDF,                     &
!JQI		      AC2, REFLSO, KGRPNT, XCGRID,                     &
!JQI		      YCGRID, CAX, CAY, RDX, RDY, LSWMAT(1,1,JABIN),   &
!JQI		      SPCSIG)                                         
       ENDIF
!
     ENDIF
!
!    *** compute source terms and fill the matrix ***
!
     CALL SOURCE (ITER   ,KWAVE               ,SPCSIG              ,   &
     SPCDIR(1,2)         ,SPCDIR(1,3)         ,COMPDA(1,JDP2)      ,   &
     SWMATR(1,1,JMATD)   ,SWMATR(1,1,JMATR)   ,ABRBOT              ,   &
     KMESPC              ,SMESPC              ,COMPDA(1,JUBOT)     ,   &
     UFRIC               ,COMPDA(1,JVX2)      ,COMPDA(1,JVY2)      ,   &
     IDCMIN              ,IDCMAX              ,IDDLOW              ,   &
     IDDTOP              ,IDWMIN              ,IDWMAX              ,   &
     ISSTOP              ,SWTSDA(1,1,1,JPWNDA),SWTSDA(1,1,1,JPWNDB),   &
     SWTSDA(1,1,1,JPWCAP),SWTSDA(1,1,1,JPBTFR),SWTSDA(1,1,1,JPWBRK),   &
     SWTSDA(1,1,1,JP4S)  ,SWTSDA(1,1,1,JP4D)  ,SWTSDA(1,1,1,JPTRI) ,   &
#    if defined (VEGETATION)
     SWTSDA(1,1,1,JPVEGT),COMPDA(1,JNPLA2)                         ,   &
#    endif
     HS                  ,ETOT                ,QBLOC               ,   &
     THETAW              ,HM                  ,FPM                 ,   &
     WIND10              ,ETOTW               ,GROWW               ,   &
     ALIMW               ,SMEBRK              ,SNLC1               ,   &
     DAL1                ,DAL2                ,DAL3                ,   &
     UE                  ,SA1                 ,SA2                 ,   &
     DA1C                ,DA1P                ,DA1M                ,   &
     DA2C                ,DA2P                ,DA2M                ,   &
     SFNL                ,DSNL                ,MEMNL4              ,   &
     WWINT               ,WWAWG               ,WWSWG               ,   &
     CGO                 ,COMPDA(1,JUSTAR)    ,COMPDA(1,JZEL)      ,   &
     SPCDIR              ,SWMATR(1,1,JDIS0)   ,SWMATR(1,1,JDIS1)   ,   &
     SZEROC              ,EPS2WC              ,DISWCP              ,   &
     WCPSME              ,WCPKME              ,WCPQB               ,   &
     WCPHM               ,XIS                 ,COMPDA(1,JFRC2)     ,   &
     IT                  ,COMPDA(1,JURSEL)    ,REFLSO              ,   &
     IG                  )   
     
     OFFSET = 1.0
     XR = 0.10
     XFILT = 0.05
     AMAX = MAXVAL(AC2(:,:,IG))
!!     DTMAX = 300.0
!!     DTMIN = 10.0      !30.0    !5.0    !0.5    !1.0   !0.1
     DTMAX = SOURCE_DTMAX
     DTMIN = SOURCE_DTMIN
     DTS = MIN(DTW-DTTOT,DTMAX)
     AFILT = MAXVAL(DAM)
     AFILT = MAX(AFILT,XFILT*AMAX)
     
     ISMAX = 1
     ISP1   = WWINT(6)
     DO IS = 1, MSC
       IF(SPCSIG(IS) < (PTRIAD(2) * SMEBRK))THEN
         ISMAX = IS
       END IF
     END DO
     ISMAX = MAX ( ISMAX , ISP1 )

!print*,'after source 1'
     DO ISC = 1,ISMAX 
       DO IDDUM = 1,MDC
         IDC = MOD ( IDDUM - 1 + MDC , MDC ) + 1
         DAMAX = MIN(DAM(IDC,ISC),MAX(XR*AC2(IDC,ISC,IG),AFILT))
	 AFAC = 1.0/MAX(1.0E-10,ABS(SWMATR(IDC,ISC,JMATR)/DAMAX))
	 DTS = MIN(DTS,AFAC/(MAX(1.0E-10,1.0+OFFSET*AFAC*             &
	       MIN(0.0,SWMATR(IDC,ISC,JMATD)))))
       END DO
     END DO  	       
     
!print*,'after source 2'
     DTS = MAX(0.5,DTS)
     DTDYN = DTDYN + DTS
     
     IDT = 1 + INT(0.99*(DTW-DTTOT)/DTS)
     DTS = (DTW-DTTOT)/REAL(IDT)
     SHAVE = DTS < DTMIN .AND. DTS < DTW-DTTOT
     DTS = MAX(DTS,MIN(DTMIN,DTW-DTTOT))
     HDT = OFFSET*DTS
     DTTOT = DTTOT + DTS
     
     IF(SHAVE)THEN
       DO ISC = 1,ISMAX    !MSC
         DO IDDUM = 1,MDC
           IDC = MOD ( IDDUM - 1 + MDC , MDC ) + 1
	   SWMATR(IDC,ISC,JMATR)=SWMATR(IDC,ISC,JMATR)*DTS/          &
	        MAX(1.0,(1.0-HDT*SWMATR(IDC,ISC,JMATD)))
           SWMATR(IDC,ISC,JMATR)=SIGN(MIN(DAM(IDC,ISC),ABS(SWMATR(IDC,ISC,JMATR))), &
	                        SWMATR(IDC,ISC,JMATR))
           AC2(IDC,ISC,IG) = MAX(0.0,AC2(IDC,ISC,IG)+SWMATR(IDC,ISC,JMATR))
	 END DO
       END DO	
     ELSE
       DO ISC = 1,ISMAX    !MSC
         DO IDDUM = 1,MDC
           IDC = MOD ( IDDUM - 1 + MDC , MDC ) + 1
	   SWMATR(IDC,ISC,JMATR)=SWMATR(IDC,ISC,JMATR)*DTS/          &
	        MAX(1.0,(1.0-HDT*SWMATR(IDC,ISC,JMATD)))   
	   AC2(IDC,ISC,IG)=MAX(0.0, AC2(IDC,ISC,IG)+SWMATR(IDC,ISC,JMATR))
	 END DO
       END DO
     END IF  	   	  			 		
!print*,'after source 3'
     IF(DTTOT >= 0.9999*DTW) EXIT
!
     IF(.NOT.ACUPDA)THEN                                             
       IF(TESTFL .AND. ITEST >= 30) WRITE (PRINTF, *) ' No update'    
       GOTO 700                                                        
     ENDIF                                                             
!
!    preparatory steps before solution of linear system
!
!     CALL SOLPRE(AC2                ,SWMATR(1,1,JAOLD)  ,             &
!                 SWMATR(1,1,JMATR)  ,SWMATR(1,1,JMATL)  ,             &
!		 SWMATR(1,1,JMATD)  ,SWMATR(1,1,JMATU)  ,             &
!		 SWMATR(1,1,JMAT5)  ,SWMATR(1,1,JMAT6)  ,             &
!		 IDCMIN             ,IDCMAX             ,             &
!		 LSWMAT(1,1,JABIN)  ,                                 &
!		 IDTOT              ,ISTOT              ,             &
!		 IDDLOW             ,IDDTOP             ,             &
!		 ISSTOP             ,                                 &
!		 SPCSIG                                 )            
!
!
!       *** test output ***
!    *** if negative action density occur rescale with a factor ***
!    *** only the sector computed is rescaled !!                ***
!
!!!     IF(BRESCL) CALL RESCALE(AC2, ISSTOP, IDCMIN, IDCMAX, NRSCAL)     
!
!    limit the change of the spectrum
!
!     IF(PNUMS(20) < 100.) CALL PHILIM (AC2, SWMATR(1,1,JAOLD),         &
!                                       CGO, KWAVE,                     &
!			               SPCSIG, LSWMAT(1,1,JABIN),      &
!				       ISLMIN, NFLIM,                  &
!				       QBLOC)                   
!
!    *** reduce the computed energy density if the value is  ***
!    *** larger then the limit value as computed in SWIND    ***
!    *** in case of first or second generation mode          ***
!
!     IF(IWIND == 1 .OR. IWIND == 2)                                    &
!       CALL WINDP3 (ISSTOP, ALIMW, AC2, GROWW, IDCMIN, IDCMAX )      
!
!
!    calculate Dissipation and Leak in all points
!
700  CONTINUE                                                          
!     CALL ADDDIS (COMPDA(1,JDISS)    ,COMPDA(1,JLEAK)    ,              &
!                  AC2                ,LSWMAT(1,1,JABIN)  ,              &
!		  SWMATR(1,1,JDIS0)  ,SWMATR(1,1,JDIS1)  ,              &
!		  SWMATR(1,1,JLEK1)  ,SPCSIG                         )
!
   END DO

!   DO IOB = 1, IOBCN
!     IOB_NODE = I_OBC_N(IOB)
!     CALL SSHAPE(AC2(1,1,IOB_NODE),SPCSIG,SPCDIR,FSHAPE,DSHAPE)
!   END DO    
 END IF

 RETURN
 END SUBROUTINE SWOMPU2
 
!****************************************************************
!
      SUBROUTINE SACCUR (DEP2       ,KGRPNT     ,XYTST      ,        &
			 AC2        ,SPCSIG     ,ACCUR      ,        &
			 HSACC1     ,HSACC2     ,SACC1      ,        &
			 SACC2      ,DELHS      ,DELTM      ,        &
			 I1MYC      ,I2MYC      )   
!     (This subroutine has not been tested yet)
!
!****************************************************************
!
      USE OCPCOMM1                                                        
      USE OCPCOMM2                                                        
      USE OCPCOMM3                                                        
      USE OCPCOMM4                                                        
      USE SWCOMM1                                                         
      USE SWCOMM2                                                         
      USE SWCOMM3                                                         
      USE SWCOMM4                                                         
      USE M_PARALL   
      USE ALL_VARS, ONLY : MT                                                     
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering                              |
!     | Fluid Mechanics Section                                   |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer(s): R.C. Ris                                   |
!     |                Modified by R. Padilla and N. Booij        |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 2004-2005  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. Authors
!
!     30.72: IJsbrand Haagsma
!     30.82: IJsbrand Haagsma
!     40.03: Nico Booij
!     40.22: John Cazes and Tim Campbell
!     40.30: Marcel Zijlema
!     40.31: Tim Campbell and John Cazes
!     40.41: Marcel Zijlema
!
!  1. Update
!
!     30.72, Nov. 97: Declaration of DDIR, PI and PI2 removed because
!                     they are common and already declared in the
!                     INCLUDE file
!     30.72, Feb. 98: Introduced generic names XCGRID, YCGRID and SPCSIG for SWAN
!     30.82, Aug. 99: Introduced a new overall measure for checking accuracy
!     30.82, Aug. 99: Changed all variables INDEX to INDX, since INDEX is reserved
!     40.03, Feb. 00: test level of message changed
!     40.22, Sep. 01: Added initialization of SACC1 and HSACC1 elements   40.22
!                     that are not wet points.                            40.22
!     40.30, Mar. 03: introduction distributed-memory approach using MPI
!     40.31, Jul. 03: some improvements and corrections w.r.t. OpenMP
!     40.41, Aug. 04: add some test output for checking accuracy
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!     40.55, Dec. 05: introducing vegetation model
!
!  2. Purpose
!
!     To check the accuracy of the final computation. If a certain
!     accuracy has been reached then terminate the iteration process
!
!  3. Method
!
!     ---
!
!  4. Argument variables
!
!     AC2         action density
!     ACCUR       indicates percentage of grid points in
!                 which accuracy is reached
!     DELHS       difference in Hs between last 2 iterations
!     DELTM       difference in Tm01 between last 2 iterations
!     DEP2        depth
!     HSACC1      significant wave height at iter-1
!     HSACC2      significant wave height at iter
!     I1MYC       lower index for thread loop over y-grid row
!     I2MYC       upper index for thread loop over y-grid row
!     KGRPNT      indirect addressing
!     SACC1       mean wave frequency at iter-1
!     SACC2       mean wave frequency at iter
!     SPCSIG      relative frequencies in computational domain
!                 in sigma-space                                          30.72
!     XYTST       coordinates of test points
!
      INTEGER I1MYC, I2MYC                                                
      REAL    ACCUR
      INTEGER KGRPNT(MXC,MYC)
      INTEGER XYTST(2*NPTST)                                              
      REAL    AC2(MDC,MSC,0:MT) ,               &
              DEP2(MT)          ,               &
	      HSACC1(MT)        ,               &
	      HSACC2(MT)        ,               &
	      SACC1(MT)         ,               &
	      SACC2(MT)         ,               &
	      DELHS(MT)         ,               &
	      DELTM(MT)
      REAL    SPCSIG(MSC)                                                
!
!  6. Local variables
!
      INTEGER  IS    ,ID    ,WETGRD,IACCUR,IX,IY,IX1,IX2,IY1,IY2
      INTEGER  WETGRDt, IACCURt                                          
!
!     INDX  : counter                                                     30.82
!     NINDX : number of gridpoints to average over                        30.82
!
      INTEGER INDX, NINDX                                                 
      INTEGER NINDXt                                                      
!
      REAL    SME_T ,SME_B ,                                     &
              TMREL ,HSREL ,TMABS ,HSABS                                  
      REAL    ARR(3)                                                      
!
!     HSMN2 : mean Hs over space at current iteration level               30.82
!     HSOVAL: Overall accuracy measure for Hs                             30.82
!     SMN2  : mean Tm over space at current iteration level               30.82
!     TMOVAL: Overall accuracy measure for Tm                             30.82
!
      REAL    HSMN2, HSOVAL, SMN2, TMOVAL                                 
      REAL    HSMN2t, SMN2t                                               
!
!     LHEAD : logical indicating to write header                          
!     TSTFL : indicates whether grid point is a test point                
!
      LOGICAL LHEAD, TSTFL                                                
!
!  7. Common blocks used
!
!
!     Place local summed variables in common block so they will           40.31
!     be scoped as shared                                                 40.31
      COMMON/SACCUR_MT_COM/HSMN2,SMN2,NINDX,WETGRD,IACCUR                 
!
!  8. Subroutines used
!
!     EQREAL           Boolean function which compares two REAL values
!     STRACE           Tracing routine for debugging
!     STPNOW           Logical indicating whether program must
!                      terminated or not
!     SWREDUCE         Performs a global reduction
!
      LOGICAL EQREAL, STPNOW
!
!  9. Subroutines calling
!
!     SWCOMP (in SWANCOM1)
!
! 12. Structure
!
!   ---------------------------------------------------------------
!   If not the first iteration, the do
!     Set old values in dummy array
!   ---------------------------------------------------------------
!   Do for every x and y
!     Compute the mean action density frequency SACC1 and the
!     and the significant waveheight HSACC1
!   ---------------------------------------------------------------
!   If relative error for mean frequency or significant wave height
!      > certain given value then increase variable with one and
!      compute the relative number of gridpoints in where the accuracy
!      has not been reached
!   ---------------------------------------------------------------
!   End of the subroutine SACCUR
!   ----------------------------------------------------------------
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SACCUR')
!
!$OMP MASTER                                                              40.31
!     Master thread initialize the shared variables                       40.31
      HSMN2  = 0.                                                         
       SMN2  = 0.                                                         
      NINDX  = 0                                                          
      WETGRD = 0                                                          
      IACCUR = 0                                                          
!$OMP END MASTER                                                          
!$OMP BARRIER                                                             
!
      IF ( LMXF ) THEN                                                    
         IX1 = 1
      ELSE
         IX1 = 1+IHALOX
      END IF
      IF ( LMXL ) THEN
         IX2 = MXC
      ELSE
         IX2 = MXC-IHALOX
      END IF
      IF ( LMYF ) THEN
         IY1 = I1MYC
      ELSE
         IY1 = 1+IHALOY
      END IF
      IF ( LMYL ) THEN
         IY2 = I2MYC
      ELSE
         IY2 = MYC-IHALOY
      END IF
!
!     *** If the computation is non steady : check the gridpoints ***
!     *** at the different "timesteps" if they are still the same ***
!     *** then WETGRD is the same. If not: change this subroutine ***
!     ***                                                         ***
!     ***   +++++++++++++++++++     +++++++++++++++++++++         ***
!     ***   +++++++++++++++++++     +++++++++++++++++++++         ***
!     ***   ++++++++     ++++++     +++++++++      ++++++         ***
!     ***   +                 +     ++++                +         ***
!     ***   +                 +     +++                 +         ***
!     ***   +                 +     ++                  +         ***
!     ***   +                 +     +                   +         ***
!     ***   +       t=0       +     +      t=t+1        +         ***
!     ***   +                 +     +                   +         ***
!     ***   +++++++++++++++++++     +++++++++++++++++++++         ***
!     ***                                                         ***
!
      WETGRDt = 0                                                         
      DO 101 IX = IX1, IX2
      DO 100 IY = IY1, IY2
        INDX = KGRPNT(IX,IY)
        IF (DEP2(INDX) .GT. DEPMIN) THEN                                  
          HSACC1(INDX) = MAX( 1.E-20 , HSACC2(INDX) )                     
          SACC1(INDX)  = MAX( 1.E-20 , SACC2(INDX)  )                     
          WETGRDt = WETGRDt + 1                                           
!       Added to initialize HSACC1 and SACC1 values.                      
        ELSE                                                              
          HSACC1(INDX) = 0.                                               
          SACC1(INDX)  = 0.                                               
        END IF                                                            
 100  CONTINUE
 101  CONTINUE
!
!     *** first criterion to terminate the iteration process ***
!     ***                                                    ***
!     *** RELATIVE error :                                   ***
!     ***               Hs2  - Hs1      Tm2 - Tm1            ***
!     ***      DREL  =  ----------  and ----------           ***
!     ***                   Hs1            Tm1               ***
!     ***                                                    ***
!     *** ABSOLUTE error :                                   ***
!     ***                                                    ***
!     ***      DHABS  =  Hs2  - Hs1 < PNUMS(2)               ***
!     ***      DTABS  =  Tm2  - Tm1 < PNUMS(3)               ***
!     ***                                                    ***
!
      DO 201 IX = IX1, IX2
      DO 200 IY = IY1, IY2
        INDX = KGRPNT(IX,IY)
!
!       *** Compute the mean ENERGY DENSITY frequency and    ***
!       *** significant waveheight over the full spectrum    ***
!       *** per gridpoint                                    ***
!
        IF (DEP2(INDX) .GT. DEPMIN) THEN                                  
          SME_T  = 0.
          SME_B  = 0.
          DO 180 IS = 1, MSC
            DO 170 ID = 1, MDC
              ACS2  = SPCSIG(IS)**2 * AC2(ID,IS,INDX)                     
              ACS3  = SPCSIG(IS) * ACS2                                   
              SME_B = SME_B + ACS2
              SME_T = SME_T + ACS3
 170        CONTINUE
 180      CONTINUE
          SME_B = SME_B * FRINTF * DDIR                                   
          SME_T = SME_T * FRINTF * DDIR                                   
!
!         *** mean frequency and significant wave height per gridpoint ***
!
          IF ( SME_B .LE. 0. ) THEN
            SME_B = 1.E-20
            SACC2(INDX) = 1.E-20                                          
            HSACC2(INDX) = 1.E-20                                         
          ELSE
            SACC2(INDX) = MAX ( 1.E-20 , (SME_T / SME_B) )                
            HSACC2(INDX) = MAX ( 1.E-20 , (4. * SQRT(SME_B)) )            
          END IF
        END IF
 200  CONTINUE
 201  CONTINUE
!
!     *** the mean significant waveheight and the mean  ***
!     *** relative frequency over the gridpoints which  ***
!     *** depth is larger than 0 m.                     ***
!     *** The amount of gridpoints is denoted with the  ***
!     *** variable : NINDX                              ***
!     *** These values are used to compute the SRELF    ***
!     *** and the HSRELF instead of SACC2 and HSACC2    ***
!
!     Note that initialization of ACCUR is not needed since it is         
!     not being summed upon                                               
      IACCURt = 0                                                         
      PI2_W     = 2.*PI_W
!
!     Calculate the mean Hs and Tm over all wet gridpoints. These means
!     are then used as an overall accuracy measure.
!
      HSMN2t= 0.                                                          
       SMN2t= 0.                                                          
      NINDXt= 0                                                           
!
      DO 301 IX = IX1, IX2
      DO 300 IY = IY1, IY2
        INDX = KGRPNT(IX,IY)
        IF (DEP2(INDX).GT.DEPMIN) THEN                                    
          HSMN2t = HSMN2t + HSACC2(INDX)                                  
          SMN2t  =  SMN2t +  SACC2(INDX)                                  
          NINDXt = NINDXt + 1                                             
        END IF                                                            
 300  CONTINUE                                                            
 301  CONTINUE
!
!     Global sum of NINDX                                                 
!$OMP ATOMIC                                                              
      NINDX = NINDX + NINDXt                                              
!$OMP ATOMIC                                                              
      HSMN2 = HSMN2 + HSMN2t                                              
!$OMP ATOMIC                                                              
      SMN2 = SMN2 + SMN2t                                                 
!
!$OMP BARRIER                                                             
!$OMP MASTER                                                              
      ARR(1) = HSMN2                                                      
      ARR(2) = SMN2                                                       
      ARR(3) = REAL(NINDX)                                                
!      CALL SWREDUCE( ARR, 3, SWREAL, SWSUM )                              
!DOS      CALL SWREDUCE( NINDX, 1, SWINT, SWSUM )                         
!DOS      ARR(3) = REAL(NINDX)                                            
!OMP#ifndef _OPENMP                                                       
      IF (STPNOW()) RETURN                                                
!OMP#endif                                                                
      HSMN2 = ARR(1) / ARR(3)                                             
       SMN2 = ARR(2) / ARR(3)                                             
!$OMP END MASTER                                                          
!$OMP BARRIER                                                             
!
!     Calculate a set of accuracy parameters based on relative, absolute
!     and overall accuracy measures for Hs and Tm
!
      LHEAD=.TRUE.                                                        
      DO 401 IX = IX1, IX2
      DO 400 IY = IY1, IY2
        INDX = KGRPNT(IX,IY)

!       --- determine whether the point is a test point                   40.41

        TSTFL = .FALSE.                                                   
        IF (NPTST.GT.0) THEN                                              
           DO 20 II = 1, NPTST                                            
              IF (IX.NE.XYTST(2*II-1)) GOTO 20                            
              IF (IY.NE.XYTST(2*II  )) GOTO 20                            
              TSTFL = .TRUE.                                              
  20       CONTINUE                                                       
        END IF                                                            

        IF ( DEP2(INDX) .GT. DEPMIN ) THEN                                
          TMREL  = ABS ( SACC2(INDX) - SACC1(INDX) ) /                &
	                 SACC1(INDX)                                      
          TMABS  = ABS ( ( PI2_W/SACC2(INDX)) - (PI2_W/SACC1(INDX)) )         
          TMOVAL = ABS ( SACC2(INDX) - SACC1(INDX) ) / SMN2               
!
          HSREL  = ABS ( HSACC2(INDX) - HSACC1(INDX) ) /              &
	                 HSACC1(INDX)                                     
          HSABS  = ABS ( HSACC2(INDX) - HSACC1(INDX) )                    
          HSOVAL = ABS ( HSACC2(INDX) - HSACC1(INDX) ) / HSMN2            
!
          IF (EQREAL(SACC1(INDX),1.E-20) .OR.                         &
	      EQREAL(SACC2(INDX),1.E-20) ) THEN                           
             DELTM(INDX) = 0.                                             
          ELSE                                                            
             DELTM(INDX) = TMABS
          END IF                                                          
          DELHS(INDX) = HSABS
!
!         *** gridpoint in which mean period and wave height ***
!         *** have reached required accuracy                 ***
!
          IF ( ITEST .GE. 30 .AND. TESTFL) THEN
            WRITE(PRINTF,3002) SACC2(INDX), SACC1(INDX),             &
	                       HSACC2(INDX), HSACC1(INDX)                 
 3002       FORMAT(' SACCUR: SA2 SA1 HSA2 HSA1       :',4E12.4)
            WRITE(PRINTF,2002) TMREL, HSREL, TMABS, HSABS
 2002       FORMAT(' SACCUR: TMREL HSREL TMABS HSABS :',4E12.4)
          ENDIF
!
          IF ( (TMREL .LE. PNUMS(1) .OR. TMOVAL .LE. PNUMS(16)) .AND.    &
	       (HSREL .LE. PNUMS(1) .OR. HSOVAL .LE. PNUMS(15)) ) THEN    
            IACCURt = IACCURt + 1                                         
          END IF
!
          IF (TSTFL) THEN                                                 
             IF (LHEAD) WRITE(PRINTF,501)                                 
             WRITE(PRINTF,502) IX+MXF-2, IY+MYF-2, HSREL, HSOVAL,       &
	                       TMREL, TMOVAL                              
 501         FORMAT(25X,'dHrel          ','dHoval         ',            &
                   'dTm01rel       ','dTm01oval ')                   
 502         FORMAT(1X,SS,'(IX,IY)=(',I5,',',I5,')','  ',               &
                    1PE13.6E2,'  ',1PE13.6E2,'  ',1PE13.6E2,'  ',       &
		    1PE13.6E2)                                            
             LHEAD=.FALSE.                                                
          END IF                                                          
        ELSE
!         *** otherwise set arrays equal 0 ***
          DELTM(INDX) = 0.0
          DELHS(INDX) = 0.0
        END IF
!
!       Test output at test points
!
        IF ( ITEST .GE. 30 .AND. TESTFL) THEN                             
           WRITE(PRINTF,1003) TMREL, TMABS, TMOVAL                        
 1003      FORMAT(' SACCUR: TMREL, TMABS, TMOVAL  :',3E12.4)              
           WRITE(PRINTF,1004) HSREL, HSABS, HSOVAL                        
 1004      FORMAT(' SACCUR: HSREL, HSABS, HSOVAL  :',3E12.4)              
        END IF                                                            
!
 400  CONTINUE
 401  CONTINUE
!
!     Global sum of IACCUR and WETGRD                                     
!$OMP ATOMIC                                                              
      IACCUR = IACCUR + IACCURt                                           
!$OMP ATOMIC                                                              
      WETGRD = WETGRD + WETGRDt                                           
!
!$OMP BARRIER                                                             
!$OMP MASTER                                                              
!
!      CALL SWREDUCE ( IACCUR, 1, SWINT, SWSUM )                           
!OMP#ifndef _OPENMP                                                       
      IF (STPNOW()) RETURN                                                
!OMP#endif                                                                
      ACCUR  = REAL(IACCUR) * 100. / ARR(3)
!$OMP END MASTER                                                          
!$OMP BARRIER                                                             
!
!     *** test output ***
!
!$OMP MASTER                                                              
      IF ( ITEST .GE. 30 ) THEN                                           
        WRITE(PRINTF,1002) PNUMS(1), PNUMS(2), PNUMS(3)
 1002   FORMAT(' SACCUR: PNUMS(1) DHABS DTABS  :',3E12.4)
        WRITE(PRINTF,1008) NINT(ARR(3)),IACCUR,ACCUR
 1008   FORMAT(' SACCUR: WETGRD IACCUR ACCUR   :',2I8,E12.4)
      END IF
!$OMP END MASTER                                                          
!
      RETURN
      END SUBROUTINE SACCUR
 
!
!****************************************************************
!
   SUBROUTINE INSAC (SPCSIG,DEP2,HSACC2,SACC2)                
!
!****************************************************************
!
!     To check the accuracy of the final computation. If a certain
!     accuracy has been reached then quit the iteration process
!
!****************************************************************
!
   USE OCPCOMM1                                                        
   USE OCPCOMM2                                                        
   USE OCPCOMM3                                                        
   USE OCPCOMM4                                                        
   USE SWCOMM1                                                         
   USE SWCOMM2                                                         
   USE SWCOMM3                                                         
   USE SWCOMM4 
!   USE ALL_VARS
   USE VARS_WAVE
   
   IMPLICIT NONE                                                        

   REAL    :: SPCSIG(MSC)                                                 
   INTEGER :: IS,ID,IG      
   INTEGER :: I1MYC, I2MYC                                               
   INTEGER :: KGRPNT(MXC,MYC)
   REAL    :: SME_T ,SME_B                                               
   REAL    :: DEP2(M) ,HSACC2(M) ,SACC2(M) 
   REAL    :: ACS2,ACS3	                                                 
!
   DO IG = 1,M
!
!    *** Compute the mean ENERGY DENSITY frequency SACC2  ***
!    *** and the wavenumber HSACC2 average over the full  ***
!    *** spectrum per gridpoint                           ***
!
     IF(DEP2(IG) > DEPMIN)THEN                                  
       SME_T  = 0.
       SME_B  = 0.
       DO IS = 1, MSC
         DO ID = 1, MDC
           ACS2 = SPCSIG(IS)**2 * AC2(ID,IS,IG)                       
           ACS3 = SPCSIG(IS) * ACS2                                    
           SME_B = SME_B + ACS2
           SME_T = SME_T + ACS3
	 END DO
       END DO	   
       SME_B = SME_B * FRINTF * DDIR                                   
       SME_T = SME_T * FRINTF * DDIR                                   
!
!      *** mean frequency and significant wave height ***
!      *** per gridpoint                               ***
!
       IF(SME_B <= 0.)THEN
         SACC2(IG)  = 1.E-20
         HSACC2(IG) = 1.E-20
       ELSE
         SACC2(IG)  = MAX ( 1.E-20 , (SME_T / SME_B) )
         HSACC2(IG) = MAX ( 1.E-20 , (4. * SQRT(SME_B)) )
       END IF
     ELSE
       SACC2(IG)  = 0.
       HSACC2(IG) = 0.
     END IF
   END DO  
!
   RETURN
   END SUBROUTINE INSAC
 
!****************************************************************
!
   SUBROUTINE SINTGRL(SPCDIR  ,KWAVE   ,DEP2    ,QB_LOC  ,            &
                      URSELL  ,RDX     ,RDY     ,AC2TOT  ,            &
		      ETOT    ,ABRBOT  ,UBOT    ,HS      ,            &
		      QB      ,HM      ,KMESPC  ,SMEBRK  ,            &
		      TMBOT   ,IG      )      
!
!****************************************************************
!
!     To compute several integrals used in SWAN and some general parameters
!
!  Method
!
!     The total energy ETOT is calculate as the following integral
!
!     ETOT = Integrate [ AC2(theta,sigma) sigma dsigma dtheta ]
!
!     To avoid too high dissipation by breaking, ETOT is maximised by a maximum
!     total energy EMAX based on the maximum wave heigth HM:
!
!     HM    = PSURF(2) depth
!
!                      2
!     EMAX  = 0.25 * HM
!
!     When EMAX > ETOT, then the action density AC2 is reduced by EMAX/ETOT
!
!     In the physicaly unrealistic case that ETOT <= 0, the integrals and other parameters
!     get values that represent a steady sea-state with wind close to zero.
!
!     The following integrals are calculated:
!
!                                                       2
!     AB2   = Integrate [ (AC2(theta,sigma) sigma / Sinh [ K(sigma) depth ]) dsigma dtheta ]
!
!     ACTOT = Integrate [ AC2(theta,sigma) dsigma dtheta ]
!
!     EDRKTOT=Integrate [ (AC2(theta,sigma) sigma / Sqrt [ K(sigma) ]) dsigma dtheta ]
!
!     EKTOT = Integrate [ AC2(theta,sigma) K(sigma) sigma dsigma dtheta ]
!
!                                               2
!     ETOT1 = Integrate [ AC2(theta,sigma) sigma dsigma dtheta ]
!
!                                               3
!     ETOT2 = Integrate [ AC2(theta,sigma) sigma dsigma dtheta ]
!
!                                               5
!     ETOT4 = Integrate [ AC2(theta,sigma) sigma dsigma dtheta ]
!
!                                                3      2
!     UB2   = Integrate [ (AC2(theta,sigma) sigma / Sinh [ K(sigma) depth ]) dsigma dtheta ]
!
!     For reasons of ??, in the calculation of UB2, AB2, ETOTM2, ETOTM4, the high frequency
!     tail is ignored.
!
!     Based on these integrals the following parameters are calculated:
!
!     ABRBOT  = Sqrt [ 2 AB2 ]
!     HS      = 4 Sqrt [ ETOT ]
!                               2
!     KM_WAM  = (ETOT / EDRKTOT)
!     QB      : computed in the subroutine FRABRE
!     SIGM01  = ETOT1 / ETOT
!     SIGM_10 = ETOT  / ACTOT
!     UBOT    = Sqrt [ UB2 ]
!     TMBOT   = 2*PI AB2 / UB2
!
!---------------------------------------------------------
!   Determine ETOT
!   If energy level is too large compared with depth
!   Then reduce action densities
!   ------------------------------------------------------
!   For all spectral frequencies do
!       determine wavenumber K and K*depth
!       For every spectral direction do
!           add AC2 to sum of action densities
!       --------------------------------------------------
!       add contributions to various moments of energy density
!   ---------------------------------------------------------
!   add tail contributions
!   determine average frequency and wavenumber
!   ----------------------------------------------------------
!   If B&J surf breaking is used
!   Then call FRABRE to compute fraction of breaking waves
!   ----------------------------------------------------------
!   determine orbital motion near the bottom
!   ----------------------------------------------------------
!****************************************************************
!
   USE OCPCOMM4                                                        
   USE SWCOMM1                                                         
   USE SWCOMM3                                                         
   USE SWCOMM4                                                         
   USE M_WCAP
!   USE ALL_VARS, ONLY : MT ,AC2
   USE VARS_WAVE, ONLY : MT ,AC2
!
   IMPLICIT NONE

   REAL, INTENT(IN)     :: DEP2(MT)
   REAL, INTENT(IN)     :: KWAVE(MSC,MICMAX)                           
   REAL, INTENT(IN)     :: RDX(10), RDY(10)                            
   REAL, INTENT(IN)     :: SPCDIR(MDC,6)
   REAL, INTENT(IN OUT) :: QB(MT)
   REAL, INTENT(IN OUT) :: UBOT(MT)
   REAL, INTENT(IN OUT) :: URSELL(MT)
   REAL, INTENT(IN OUT) :: TMBOT(MT)                                
   REAL, INTENT(OUT)    :: ABRBOT, ETOT, HM, HS, QB_LOC
   REAL, INTENT(OUT)    :: AC2TOT, KMESPC, SMEBRK
   INTEGER, SAVE :: IENT = 0
   REAL              :: AB2, EMAX, UB2
   REAL              :: FRINTF_X_DDIR
   REAL              :: BRCOEF   ! variable breaking coefficient (calc. in BRKPAR)  40.22
   REAL              :: ETOT_DSIG(MSC)                                 
   REAL              :: ACTOT_DSIG(MSC), ETOT_DSHKD2_DSIG(MSC)         
   REAL              :: ETOT_DRK_DSIG(MSC), ETOT_SIG2_DSIG(MSC)        
   REAL              :: ETOT_K_DSIG(MSC), ETOT_SIG_DSIG(MSC)           
   REAL              :: ETOT_SIG2_DSHKD2_DSIG(MSC)                     
   REAL              :: ETOT_SIG4_DSIG(MSC)                            
   REAL              :: SINH_K_X_DEP_2(MSC)  
   INTEGER           :: IG                          
!
!  Initialisation
!
   KM_WAM  = 10.
   KM01    = 10.
!
   SIGM01  = 10.
   SIGM_10 = 10.
!
   HS      = 0.
   HM      = 0.1
!
   QB(IG)   = 0.
   ABRBOT   = 0.001
   UBOT(IG) = 0.
   TMBOT(IG)= 0.
!
!  Calculate total spectral energy:
!
!
   FRINTF_X_DDIR = FRINTF * DDIR
   ETOT_DSIG(:)  = SUM(AC2(:,:,IG),DIM=1) * SIGPOW(:,2) *     &
                   FRINTF_X_DDIR
   ETOT          = SUM(ETOT_DSIG)
!
!  *** add high frequency tail ***
!
   ETOT = ETOT + ETOT_DSIG(MSC) * PWTAIL(6) / FRINTF
!
   IF(ISURF == 1)THEN
     HM   = PSURF(2) * DEP2(IG)                                  
   ELSE IF(ISURF >= 2)THEN
!    Calulate the correct breaking coefficient BRCOEF
     CALL BRKPAR (BRCOEF  ,SPCDIR(1,2), SPCDIR(1,3), AC2     ,      &
                  SIGPOW(:,1), DEP2, IG )            
!!$                  SIGPOW(:,1), DEP2, RDX, RDY, IG )            
     HM   = BRCOEF * DEP2(IG)                                   
   ELSE
!    breaking disabled, assign very high value to HM
     HM   = 100.                                                       
   ENDIF
!
   EMAX = 0.25 * HM**2
!
!  --- reduce action density if necessary
!
   IF(ACUPDA .AND. ETOT > EMAX .AND. ISURF >= 1)THEN           
     AC2(:,:,IG) = MAX(0.,(EMAX/ETOT)*AC2(:,:,IG))
!
     IF(TESTFL.AND.ITEST >= 80) WRITE (PRTEST,7) DEP2(IG), EMAX, ETOT
7    FORMAT (' energy is reduced in SINTGRL', 4(1x, e12.4))
!
!    Correct value for ETOT
!
     ETOT = EMAX
!
   ENDIF
!
   IF(ETOT > 0.)THEN
!
!  Calculate all other integrals
!
!
     SINH_K_X_DEP_2(:)   = SINH(MIN(30.,KWAVE(:,1)*DEP2(IG)))**2
     ACTOT_DSIG(:)       = SUM(AC2(:,:,IG),DIM=1) * SIGPOW(:,1) * FRINTF_X_DDIR
     ETOT_DSIG(:)        = ACTOT_DSIG(:) * SIGPOW(:,1)
     ETOT_SIG_DSIG(:)    = ACTOT_DSIG(:) * SIGPOW(:,2)
     ETOT_SIG2_DSIG(:)   = ACTOT_DSIG(:) * SIGPOW(:,3)
     ETOT_SIG4_DSIG(:)   = ACTOT_DSIG(:) * SIGPOW(:,5)
     ETOT_DRK_DSIG(:)    = ETOT_DSIG(:) / SQRT(KWAVE(:,1))
     ETOT_K_DSIG(:)      = ETOT_DSIG(:) * KWAVE(:,1)
     ETOT_DSHKD2_DSIG(:) = ETOT_DSIG(:) / SINH_K_X_DEP_2(:)
     ETOT_SIG2_DSHKD2_DSIG(:) = ETOT_DSHKD2_DSIG(:) * SIGPOW(:,2)
!
     ACTOT         = SUM(ACTOT_DSIG)
     ETOT1         = SUM(ETOT_SIG_DSIG)
     ETOT2         = SUM(ETOT_SIG2_DSIG)
     ETOT4         = SUM(ETOT_SIG4_DSIG)
     EDRKTOT       = SUM(ETOT_DRK_DSIG)
     EKTOT         = SUM(ETOT_K_DSIG)
     UB2           = SUM(ETOT_SIG2_DSHKD2_DSIG)
     AB2           = SUM(ETOT_DSHKD2_DSIG)
!
!  Add high frequency tails
!
     ACTOT       = ACTOT + PWTAIL(5) * ACTOT_DSIG(MSC) / FRINTF
     ETOT1       = ETOT1 + PWTAIL(7) * ETOT_DSIG(MSC) * SIGPOW(MSC,1) / FRINTF
     EDRKTOT     = EDRKTOT + PWTAIL(5) * ETOT_DSIG(MSC) / (SQRT(KWAVE(MSC,1)) * FRINTF)
     EKTOT       = EKTOT + PWTAIL(8) * ETOT_DSIG(MSC) * KWAVE(MSC,1) / FRINTF
!
!  Calculate the mean frequencies SIGM01 and SIGM_10,
!  mean wavenumbers KM_WAM, KM01 and significant waveheight HS
!
     IF(ETOT1 > 0.) SIGM01  = ETOT1 / ETOT
     IF(EKTOT > 0.) KM01    = EKTOT / ETOT
     IF(ACTOT > 0.) SIGM_10 = ETOT / ACTOT
     IF(EDRKTOT > 0.)THEN
       KM_WAM  = ( ETOT / EDRKTOT )**2.
     ENDIF
     IF(ETOT > 1.E-20)THEN
       HS       = 4. * SQRT (ETOT)
     END IF
!
!  Calculate QB, when breaking is activated
!
     IF(ISURF > 0) CALL FRABRE (HM, ETOT, QB(IG))
!
!  Calculate the orbital velocity UBOT, orbital excursion ABRBOT and
!  bottom wave period TMBOT                                            
!
     IF(UB2 > 0.) UBOT(IG) = SQRT ( UB2 )
     IF(AB2 > 0.) ABRBOT = SQRT (2. *  AB2)
     IF(UB2 > 0. .AND. AB2 > 0.) TMBOT(IG) = PI2_W*SQRT(AB2/UB2)  
!
   ENDIF
!
!  *** calculate Ursell number ***
!  *** update only for first encounter in a sweep                      
!
   IF(ITRIAD > 0)THEN                                             
     URSELL(IG) = (GRAV_W*HS) /                                    &
                          (2.*SQRT(2.)*SIGM01**2*DEP2(IG)**2)
   ELSE                                                              
     URSELL(IG) = 0.                                             
   END IF                                                            
!
!  *** test output ***
!
!   IF (TESTFL .AND. ITEST.GE.60) THEN
!     WRITE(PRTEST, 901) ETOT, HS, SIGM_10, KM_WAM, ABRBOT
!901  FORMAT (' SINTGRL: ETOT Hs Sigma K Aorb', 5(1X, E11.4))
!   END IF
!
!
!  Set variables used outside the whitecapping scope
!
   AC2TOT = ACTOT
   KMESPC = KM_WAM
   SMEBRK = SIGM01
   QB_LOC = QB(IG)
!
   RETURN
!
   END SUBROUTINE SINTGRL
!
!****************************************************************
!
      SUBROUTINE SOLPRE (AC2         ,AC2OLD      ,             &
                         IMATRA      ,IMATLA      ,             &
			 IMATDA      ,IMATUA      ,             &
			 IMAT5L      ,IMAT6U      ,             &
			 IDCMIN      ,IDCMAX      ,             &
			 ANYBIN      ,                          &
			 IDTOT       ,ISTOT       ,             &
			 IDDLOW      ,IDDTOP      ,             &
			 ISSTOP      ,                          &
			 SPCSIG                   )                      

!    (This subroutine has not been tested yet)
!
!****************************************************************
!
      USE SWCOMM2                                                        
      USE SWCOMM3                                                        
      USE SWCOMM4                                                        
      USE OCPCOMM4                                                       
      USE M_PARALL
      USE ALL_VARS, ONLY : MT                                                       
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering                              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmers: R.C. Ris, N. Booij,                          |
!     |              IJ.G. Haagsma, A.T.M.M. Kieftenburg,         |
!     |              M. Zijlema, E.E. Kriezi,                     |
!     |              R. Padilla-Hernandez, L.H. Holthuijsen       |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 2004-2005  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. Authors
!
!     40.00: Nico Booij
!     40.23: Marcel Zijlema
!     40.30: Marcel Zijlema
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     40.00, Feb. 99: New subroutine common tasks before solution of linear
!                     system (software moved from SOLBAND, SOLMAT and SOLMT1)
!     40.23, Aug. 02: implementation of under-relaxation technique
!     40.30, Mar. 03: correcting indices of test point with offsets MXF, MYF
!     40.41, Aug. 04: code optimized
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Copy local spectrum to array AC2OLD before solving system,
!     apply under-relaxation approach in active bins, fill matrix
!     arrays for non-active bins and write test output
!
!  3. Method
!
!     ---
!
!  4. Argument variables
!
!        one and more dimensional arrays:
!        ---------------------------------
!        AC2       4D    Action density as function of D,S,X,Y and T
!        AC2OLD    2D    Values of action density at previous iteration
!        IMATDA    2D    Coefficients of diagonal of matrix
!        IMATLA    2D    Coefficients of lower diagonal of matrix
!        IMATUA    2D    Coefficients of upper diagonal of matrix
!        IMATRA    2D    Coefficients of right hand side of matrix
!        IMAT5L    2D    Coefficients for implicit calculation in
!                        frequency space (lower diagonal)
!        IMAT6U    2D    Coefficients for implicit calculation in
!                        frequency space (upper diagonal)
!        SPCSIG    1D    Relative frequencies in sigma-space
!
      REAL     AC2(MDC,MSC,0:MT)            ,                    &
               IMATRA(MDC,MSC)              ,                    &
	       IMATLA(MDC,MSC)              ,                    &
	       IMATDA(MDC,MSC)              ,                    &
	       IMATUA(MDC,MSC)              ,                    &
	       IMAT5L(MDC,MSC)              ,                    &
	       IMAT6U(MDC,MSC)              ,                    &
	       AC2OLD(MDC,MSC)              ,                    &
	       SPCSIG(MSC)
!
!        IDCMIN    1D    Integer array containing minimum counter
!        IDCMAX    1D    Integer array containing maximum counter
!
      INTEGER  IDCMIN(MSC)                  ,                    &
               IDCMAX(MSC)
!
!        ANYBIN    2D    Logical array. if a certain bin is enclosed
!                        in a sweep then ANYBIN is TRUE . array is
!                        used to determine whether some coefficients
!                        in the array have to be changed
!
      LOGICAL  ANYBIN(MDC,MSC)
!
!  7. Common blocks used
!
!
!     5. SUBROUTINES CALLING
!
!        SWOMPU
!
!     6. SUBROUTINES USED
!
!        NONE
!
!     7. ERROR MESSAGES
!
!        ---
!
!     8. REMARKS
!
!
!     9. STRUCTURE
!
!
!     10. SOURCE
!
!************************************************************************
!
      INTEGER  IS      ,ID      ,IDDUM   ,                        &
               IDDLOW  ,                                          &
	       IDDTOP  ,IDTOT   ,ISTOT   ,ISSTOP
      REAL     ALFA                                                      
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SOLPRE')
!
!     --- apply under-relaxation approach, if requested
!
      ALFA = PNUMS(30)                                                    
      IF (ALFA.GT.0.) THEN                                                
         DO IS = 1, ISSTOP
            DO IDDUM = IDCMIN(IS), IDCMAX(IS)
               ID = MOD(IDDUM-1 + MDC, MDC) + 1
               IMATDA(ID,IS) = IMATDA(ID,IS) + ALFA*SPCSIG(IS)            
               IMATRA(ID,IS) = IMATRA(ID,IS) + ALFA*SPCSIG(IS)*       &
	                       AC2(ID,IS,KCGRD(1))        
            END DO
         END DO
      END IF
!
!     --- when ambient currents are involved or when the spectral space
!         is not a circular one (use SECTOR instead of CIRCLE in command
!         CGRID), some bins do not fall within the current sweep (in
!         particular, when SECTOR = 0 or 4 i.c. ICUR=1 or SECTOR = 4 i.c.
!         FULCIR=.FALSE., see routine SWPSEL for meaning of SECTOR). For
!         such bins, the corresponding rows in the matrix are reset such
!         that the solution AC2 does not change: the main diagonal is set
!         to 1, the off-diagonals are set to 0 and the righ-hand side is
!         set to AC2.
!
      IF ( ICUR.EQ.1 .OR. .NOT.FULCIR ) THEN                              
         DO IS = 1, MSC
           DO ID = 1, MDC
             IF ( .NOT. ANYBIN(ID,IS) ) THEN
               IMATLA(ID,IS) = 0.
               IMATDA(ID,IS) = 1.
               IMATUA(ID,IS) = 0.
               IMATRA(ID,IS) = AC2(ID,IS,KCGRD(1))                        
               IMAT5L(ID,IS) = 0.
               IMAT6U(ID,IS) = 0.
             END IF
           ENDDO
         ENDDO
      END IF
!
!     *** the action density is stored in an auxiliary array AC2OLD ***
!
      DO IS = 1, MSC
        DO ID = 1, MDC
          AC2OLD(ID,IS) = AC2(ID,IS,KCGRD(1))
        ENDDO
      ENDDO
!
!     *** test output ***
!
      IF ( TESTFL .AND. ITEST .GE. 70 ) THEN
        WRITE (PRINTF,6120) IXCGRD(1)+MXF-2, IYCGRD(1)+MYF-2              
 6120   FORMAT (' SOLPRE: Matrix values for point:', 2I5)
        WRITE (PRINTF,6121)
 6121   FORMAT ('  bin     diagonal     r.h.s.        ID-1        ID+1',  &
                '         IS-1         IS+1')
!
        DO IS = 1, ISSTOP
          ID_MIN = IDCMIN(IS)
          ID_MAX = IDCMAX(IS)
          DO IDDUM = ID_MIN, ID_MAX
            ID = MOD(IDDUM-1 + MDC, MDC) + 1
            IF ( DYNDEP .OR. ICUR .EQ. 1 ) THEN                           
              WRITE(PRINTF,6620) ID, IS, IMATDA(ID,IS), IMATRA(ID,IS),   &
	      IMATLA(ID,IS), IMATUA(ID,IS), IMAT5L(ID,IS), IMAT6U(ID,IS)
 6620         FORMAT(2I3,6(1X,E12.4))
            ELSE
              WRITE(PRINTF,6620) ID, IS, IMATDA(ID,IS), IMATRA(ID,IS),   &
	      IMATLA(ID,IS), IMATUA(ID,IS)
            ENDIF
          ENDDO
        ENDDO
      END IF
!
      RETURN
      END
!
!****************************************************************
!
   SUBROUTINE SOURCE (ITER       ,KWAVE      ,SPCSIG     ,      &
		      ECOS       ,ESIN       ,DEP2       ,      &
		      IMATDA     ,IMATRA     ,ABRBOT     ,      &
		      KMESPC     ,SMESPC     ,UBOT       ,      &
		      UFRIC      ,UX2        ,UY2        ,      &
		      IDCMIN     ,IDCMAX     ,IDDLOW     ,      &
		      IDDTOP     ,IDWMIN     ,IDWMAX     ,      &
		      ISSTOP     ,PLWNDA     ,PLWNDB     ,      &
		      PLWCAP     ,PLBTFR     ,PLWBRK     ,      &
		      PLNL4S     ,PLNL4D     ,PLTRI      ,      &
# if defined (VEGETATION)
                      PLVEGT     ,NPLA2      ,                  &
# endif
		      HS         ,ETOT       ,QBLOC      ,      &
		      THETAW     ,HM         ,FPM        ,      &
		      WIND10     ,ETOTW      ,GROWW      ,      &
		      ALIMW      ,SMEBRK     ,SNLC1      ,      &
		      DAL1       ,DAL2       ,DAL3       ,      &
		      UE         ,SA1        ,SA2        ,      &
		      DA1C       ,DA1P       ,DA1M       ,      &
		      DA2C       ,DA2P       ,DA2M       ,      &
		      SFNL       ,DSNL       ,MEMNL4     ,      &
		      WWINT      ,WWAWG      ,WWSWG      ,      &
		      CGO        ,USTAR      ,ZELEN      ,      &
		      SPCDIR     ,DISSC0     ,DISSC1     ,      &
		      SZEROC     ,EPS2WC     ,DISWCP     ,      &
		      WCPSME     ,WCPKME     ,WCPQB      ,      &
		      WCPHM      ,XIS        ,FRCOEF     ,      &
		      IT         ,URSELL     ,REFLSO     ,      &
		      IG         )
!     (This subroutine has not been fully tested. Only tested for 
!      case IQUAD = 2)
!
!****************************************************************
!
!     to compute the source terms, i.e., bottom friction,
!     wave breaking, wind input, white capping and non linear
!     wave wave interactions
!
!     Coefficients for the arrays:
!     -----------------------------
!                         default
!                         value:
!
!     PBOT(1)   = CFC      0.005    (Putnam and Collins equation)
!     PBOT(2)   = CFW      0.01     (Putnam and Collins equation)
!     PBOT(3)   = GAMJNS   0.0038   (Jonswap equation)
!     PBOT(4)   = MF      -0.08     (Madsen et al. equation)
!     PBOT(5)   = KN       0.05     (Madsen et al. bottom roughness)
!
!     PNUMS(*)  =
!
!     PSURF(1)  = ALFA     1.0      (Battjes & Janssen, 1978)
!     PSURF(2)  = GAMMA    0.8      (Breaking criterium)
!
!     PWCAP(1)  = ALFAWC   2.36e-5  (Emperical coefficient)
!     PWCAP(2)  = ALFAPM   3.02E-3  (Alpha of Pierson Moskowitz frequency)
!     PWCAP(3)  = CFJANS   4.5
!     PWCAP(4)  = DELTA    0.5
!     PWCAP(5)  = CFLHIG   1.
!     PWCAP(6)  = GAMBTJ   0.88     (Steepness limited wave breaking )
!
!     PWIND(1)  = CF10     188.0    (second generation wind growth model)
!     PWIND(2)  = CF20     0.59     (second generation wind growth model)
!     PWIND(3)  = CF30     0.12     (second generation wind growth model)
!     PWIND(4)  = CF40     250.0    (second generation wind growth model)
!     PWIND(5)  = CF50     0.0023   (second generation wind growth model)
!     PWIND(6)  = CF60    -0.2233   (second generation wind growth model)
!     PWIND(7)  = CF70     0.       (second generation wind growth model)
!     PWIND(8)  = CF80    -0.56     (second generation wind growth model)
!     PWIND(9)  = RHOAW    0.00125  (density air / density water)
!     PWIND(10) = EDMLPM   0.0036   (limit energy Pierson Moskowitz)
!     PWIND(11) = CDRAG    0.0012   (drag coefficient)
!     PWIND(12) = UMIN     1.0      (minimum wind velocity)
!     PWIND(13) = PMLM     0.13     (  )
!
!     arrays for Janssen (`89)
!     -----------
!     PWIND(14) 1D    alfa (which is tuned at 0.01)
!     PWIND(15) 1D    Kappa ( 0.41)
!     PWIND(16) 1D    Rho air (1.32)
!     PWIND(17) 1D    Rho water (1025)
!
!   ------------------------------------------------------------
!   If SBOT is on (IBOT > 0 ) then,
!     Call SBOT  to compute the source term due to bottom friction
!      according to Hasselmann et al. (1974), Putnam and Jonsson (1949)
!      or Madsen et al. (1991)
!   ------------------------------------------------------------
!   If SVEG is on (IVEG > 0 ) then,
!     Call SVEG  to compute the source term due to vegetation
!      dissipation according to Dalrymple (1984)
!   ------------------------------------------------------------
!   If SSURF is on (ISURF > 0 ) then,
!     Call SSURF to compute the source term due to wave breaking
!       according to Battjes and Janssen (1978)
!   ------------------------------------------------------------
!   IF IWIND =1 OR IWIND =2 THEN
!     Call WNDPAR (first or second generation mode of source terms
!                  using the DOLPHIN-B formulations)
!
!   else if IWIND = 3 then
!     input source term according to Snyder (1981)
!     Call SWIND3
!   else if IWIND = 4 then
!     input source term according to Janssen (1989,1991)
!     Call SWIND4
!   else if IWIND = 5 then
!     input source term according to Yan (1989) [reduces to Snyder form
!     for low frequencies and to Plant's (1982) form for high freq.
!     Call SWIND5
!   ------------------------------------------------------------
!   If IWCAP > 1 then
!     Call SWCAP to compute the source term for white capping             40.02
!   ---------------------------------------------------------------------
!   If ITRIAD > 0 then                                                    30.80
!      Call SWLTA to compute the nonlinear 3 wave-wave interactions
!      based on the LTA technique (Eldeberky, 1996)
!   ---------------------------------------------------------------------
!   If Ursell < Urmax
!   Then If IQUAD = 1
!        Then Call SWSNL1 to compute the nonlinear 4-wave interactions
!             semi-implicit per sweep direction
!        Else if IQUAD = 2
!        Then Call SWSNL2 to compute the nonlinear 4-wave interactions
!             fully explicit per sweep direction
!        Else if IQUAD = 3
!        Then Call SWSNL3 to compute the nonlinear 4-wave interactions
!             fully explicit per iteration
!             Call FILNL3 to get values for interactions from array
!             for full circle
!   ---------------------------------------------------------------------
!
!****************************************************************

   USE OCPCOMM1                                                        
   USE OCPCOMM2                                                        
   USE OCPCOMM3                                                        
   USE OCPCOMM4                                                        
   USE SWCOMM1                                                         
   USE SWCOMM2                                                         
   USE SWCOMM3                                                         
   USE SWCOMM4                                                         
   USE M_SNL4  
!   USE ALL_VARS, ONLY : MT,AC2                                                        
   USE VARS_WAVE, ONLY : MT,AC2                                                        

   IMPLICIT NONE                                                      

   REAL    ECOS(MDC)
   REAL    ESIN(MDC)
   REAL    SPCDIR(MDC,6)                                               
   REAL    SPCSIG(MSC)                                                 

   INTEGER IQERR                                                    
!
   INTEGER           :: ID, IDIA, IDC, IERR, IENT, IS, ISC, IT         
   INTEGER           :: N2, LMAX  
   INTEGER           :: IG                                     

   REAL              :: DQ, DQ2, DT2                                   
   REAL, ALLOCATABLE :: ACRIAM(:,:), SNLRIAM(:,:), SIGRIAM(:)          

   INTEGER  ITER    ,IDWMIN  ,IDWMAX  ,ISSTOP  ,             &
            IDDTOP  ,IDDLOW  ,IX      ,IY
!
   REAL     ABRBOT  ,ETOT    ,HM      ,QBLOC   ,ETOTW   ,             &
            FPM     ,WIND10  ,THETAW  ,SMESPC  ,KMESPC  ,             &
	    SNLC1   ,FACHFR  ,DAL1    ,DAL2    ,DAL3    ,             &
	    UFRIC   ,SMEBRK  ,HS      ,SZEROC  ,EPS2WC  ,             &
	    DISWCP  ,WCPQB   ,WCPHM   ,WCPSME  ,WCPKME  ,             &
	    XIS
!
   REAL  :: DEP2(MT)
   REAL  :: ALIMW(MDC,MSC)
   REAL  :: IMATDA(MDC,MSC)
   REAL  :: IMATRA(MDC,MSC)
   REAL  :: KWAVE(MSC,MICMAX)                                          
   REAL  :: UBOT(MT)
   REAL  :: UX2(MT)
   REAL  :: UY2(MT)
   REAL  :: UE(MSC4MI:MSC4MA , MDC4MI:MDC4MA )
   REAL  :: SA1(MSC4MI:MSC4MA , MDC4MI:MDC4MA )
   REAL  :: SA2(MSC4MI:MSC4MA , MDC4MI:MDC4MA )
   REAL  :: DA1C(MSC4MI:MSC4MA , MDC4MI:MDC4MA )
   REAL  :: DA1P(MSC4MI:MSC4MA , MDC4MI:MDC4MA )
   REAL  :: DA1M(MSC4MI:MSC4MA , MDC4MI:MDC4MA )
   REAL  :: DA2C(MSC4MI:MSC4MA , MDC4MI:MDC4MA )
   REAL  :: DA2P(MSC4MI:MSC4MA , MDC4MI:MDC4MA )
   REAL  :: DA2M(MSC4MI:MSC4MA , MDC4MI:MDC4MA )
   REAL  :: SFNL(MSC4MI:MSC4MA , MDC4MI:MDC4MA )
   REAL  :: DSNL(MSC4MI:MSC4MA , MDC4MI:MDC4MA )
   REAL  :: MEMNL4(MDC,MSC,MT)
   REAL  :: PLWNDA(MDC,MSC,NPTST)
   REAL  :: PLWNDB(MDC,MSC,NPTST)
   REAL  :: PLWCAP(MDC,MSC,NPTST)
   REAL  :: PLBTFR(MDC,MSC,NPTST)
   REAL  :: PLWBRK(MDC,MSC,NPTST)
   REAL  :: PLNL4S(MDC,MSC,NPTST)
   REAL  :: PLNL4D(MDC,MSC,NPTST)
   REAL  :: PLTRI(MDC,MSC,NPTST)
#  if defined (VEGETATION)
   REAL  :: PLVEGT(MDC,MSC,NPTST)
   REAL  :: NPLA2(MT)
#  endif
   REAL  :: WWAWG(*)
   REAL  :: WWSWG(*)
   REAL  :: CGO(MSC,MICMAX)                                            
   REAL  :: USTAR(MT)
   REAL  :: ZELEN(MT)
   REAL  :: DISSC0(MDC,MSC)
   REAL  :: DISSC1(MDC,MSC)
   REAL  :: URSELL(MT)                                              
   REAL  :: FRCOEF(MT)                                              
   REAL  :: REFLSO(MDC,MSC)                                            
!
   INTEGER  IDCMIN(MSC)    ,IDCMAX(MSC)    ,WWINT(*)
!
   LOGICAL  GROWW(MDC,MSC)

!
!  *** set coefficients in the arrays 0 ***
!
   DO IS = 1, MSC                                                      
     DO ID = 1, MDC
       IMATRA(ID,IS) = 0.
       IMATDA(ID,IS) = 0.
     ENDDO
   ENDDO
!
!  *** set all dissipation coeff at 0 ***
!
   DO ISC = 1, MSC
     DO IDC = 1, MDC
       DISSC0(IDC,ISC) = 0.
       DISSC1(IDC,ISC) = 0.
     ENDDO
   ENDDO
!
!   PRINT*,'IBOT=',IBOT
   IF(IBOT >= 1)THEN
!
!    *** wave-bottom interactions ***
!
     CALL SBOT (ABRBOT   ,DEP2     ,ECOS     ,ESIN     ,KWAVE    ,     &
                SPCSIG   ,UBOT     ,UX2      ,UY2      ,IDCMIN   ,     &
		IDCMAX   ,ISSTOP   ,VARFR    ,FRCOEF   ,IG       ,     &
		IMATRA   )
   END IF
!
#  if defined (VEGETATION)
   IF ( IVEG.GE.1 ) THEN
!
!     *** wave-vegetation interactions ***       
!
!        *** energy dissipation according to Dalrymple (1984)  
!
         CALL SVEG (DEP2   ,IMATDA   ,ETOT   ,SMEBRK    ,  KMESPC &
              &,PLVEGT   ,  IDCMIN , IDCMAX   ,ISSTOP ,DISSC1    ,   &
              & NPLA2    ,      IG , IMATRA  )
   END IF
#  endif
!
!   PRINT*,'ISURF=',ISURF
   IF(ISURF >= 1)THEN
!
!    *** wave breaking with Kirby type formulation (f/fm)^2 ***
!
!    *** wave breaking according to Battjes and Janssen (1978) ***
!
!     print*,'before SSURF, IG=',IG
     CALL SSURF (ETOT    ,HM      ,QBLOC   ,SMEBRK  ,                  &
                 IMATRA  ,IMATDA  ,IDCMIN  ,IDCMAX  ,                  &
		 PLWBRK  ,ISSTOP  ,IG      )
!
   END IF

!   PRINT*,'IWIND=',IWIND
   IF(IWIND >= 3)THEN
!
!    *** linear wind input according to Cavaleri and Malanotte ***
!    *** Rizolli (1981) for a third generation mode of SWAN    ***
!
!      PRINT*,'PWIND(31)=',PWIND(31)
     IF (PWIND(31) > 1.E-20)                                           &
       CALL SWIND0 (IDCMIN  ,IDCMAX  ,ISSTOP  ,SPCSIG  ,THETAW  ,      &
		    UFRIC   ,FPM     ,PLWNDA  ,IMATRA  ,SPCDIR  )
   ENDIF
!
!   print*,'IWIND=',IWIND
!   stop
   IF(IWIND == 1 .OR. IWIND == 2 )THEN    
!print*, 'before WNDPAR'                      
!
     CALL WNDPAR (ISSTOP,IDWMIN,IDWMAX,IDCMIN,IDCMAX,                  &
                  DEP2  ,WIND10,THETAW,AC2   ,KWAVE ,                  &
		  IMATRA,IMATDA,SPCSIG,CGO   ,ALIMW ,                  &
		  GROWW ,ETOTW ,PLWNDA,PLWNDB,SPCDIR,                  &
		  ITER  ) 
!print*,'after WNDPAR'         
!
!
   ELSE IF(IWIND == 3)THEN
!
!    *** Wind input according to Snyder et al (1981) ***
!
!     print*,'before SWIND3'
     CALL SWIND3 (SPCSIG  ,THETAW  ,IMATDA  ,KWAVE   ,IMATRA  ,        &
		  IDCMIN  ,IDCMAX  ,UFRIC   ,FPM     ,PLWNDB  ,        &
		  ISSTOP  ,SPCDIR  ,IG      )
!     print*,'after SWIND3'
!
   ELSE IF(IWIND == 4)THEN
!
!    *** Wind input according to Janssen (1989,1991) ***
!
     CALL SWIND4  (IDWMIN  ,IDWMAX  ,SPCSIG  ,WIND10  ,THETAW  ,       &
                   XIS     ,DDIR    ,KWAVE   ,IMATRA  ,IMATDA  ,       &
		   IDCMIN  ,IDCMAX  ,UFRIC   ,PLWNDB  ,ISSTOP  ,       &
		   ITER    ,USTAR   ,ZELEN   ,SPCDIR  ,IT      ,       &
		   IG      )
!
   ELSE IF(IWIND == 5)THEN
!
!    *** Wind input according to Yan (1989) ***
!
     CALL SWIND5 (SPCSIG  ,THETAW  ,ISSTOP  ,UFRIC   ,KWAVE   ,        &
                  IMATRA  ,IMATDA  ,IDCMIN  ,IDCMAX  ,PLWNDB  ,        &
		  SPCDIR  ,IG      )

   END IF
!
!  Calculate whitecapping source term (six formulations)               
!
!   print*,'IWCAP=',IWCAP
   IF(IWCAP >= 1) CALL SWCAP (SPCDIR  ,SPCSIG  ,KWAVE   ,IDCMIN  ,     &
                              IDCMAX  ,ISSTOP  ,ETOT    ,IMATDA  ,     &
			      IMATRA  ,PLWCAP  ,CGO     ,UFRIC   ,     &
			      DEP2    ,DISSC1  ,DISSC0  ,IG      )    
!
!  compute nonlinear interactions, starting with triads                 
!
!   print*,'ITRIAD=',ITRIAD,ICUR,ITER
   IF(ITRIAD > 0)THEN
!
!    *** compute the 3 wave-wave interactions if in each ***            
!    *** geographical gridpoint a continuous spectrum    ***
!    *** is present, i.e., after first iteration         ***
!
     IF(ICUR == 0 .AND. ITER >= 1)THEN
!
       CALL SWLTA ( IG     ,DEP2   ,CGO    ,SPCSIG ,                   &
	            KWAVE  ,IMATRA ,IMATDA ,IDDLOW ,                   &
		    IDDTOP ,ISSTOP ,IDCMIN ,IDCMAX ,                   &
		    HS     ,SMEBRK ,PLTRI  ,URSELL )
!
     ELSE IF(ICUR == 1 .AND. ITER > 1)THEN
!
       CALL SWLTA ( IG     ,DEP2   ,CGO    ,SPCSIG ,                   &
	            KWAVE  ,IMATRA ,IMATDA ,IDDLOW ,                   &
		    IDDTOP ,ISSTOP ,IDCMIN ,IDCMAX ,                   &
		    HS     ,SMEBRK ,PLTRI  ,URSELL )
!
     ENDIF
!
!print*,'after SWLTA'
   ENDIF
!
!  --- compute quadruplet interactions if Ursell number < Urmax        
!      (usually, Urmax = 0.1, but here Urmax = 10)                     
!
   IF(URSELL(IG) < PTRIAD(3))THEN                             
!
!    *** compute the counters for the nonlinear four ***
!    *** wave-wave interactions in spectral space    ***
!    *** and high frequency factor                   ***
!
!     print*,'IQUAD=',IQUAD
     IF(IQUAD >= 1)THEN
       CALL RANGE4(WWINT,IDDLOW,IDDTOP)                            
       FACHFR = 1. / XIS ** PWTAIL(1)                                
     ENDIF

!
     IF(IQUAD == 1)THEN
!
!    *** semi-implicit calculation for al the bins that fall ***
!    *** within a sweep. No additional array is required     ***
!
       CALL SWSNL1 (WWINT   ,WWAWG   ,WWSWG   ,IDCMIN  ,IDCMAX  ,      &
                    UE      ,SA1     ,SA2     ,DA1C    ,DA1P    ,      &
		    DA1M    ,DA2C    ,DA2P    ,DA2M    ,SPCSIG  ,      &
		    SNLC1   ,KMESPC  ,FACHFR  ,ISSTOP  ,DAL1    ,      &
		    DAL2    ,DAL3    ,SFNL    ,DSNL    ,DEP2    ,      &
		    AC2     ,IMATDA  ,IMATRA  ,PLNL4S  ,PLNL4D  ,      &
		    IDDLOW  ,IDDTOP  )   
!print*,'after SWSNL1'  
!
     ELSE IF(IQUAD == 2)THEN
!
!      *** fully explicit calculation for al the bins that fall ***
!      *** within a sweep. No additional array is required      ***
!
!       print*,'before calling SWSNL2'
       CALL SWSNL2 (IDDLOW  ,IDDTOP  ,WWINT   ,WWAWG   ,UE      ,      &
                    SA1     ,ISSTOP  ,SA2     ,SPCSIG  ,SNLC1   ,      &
		    DAL1    ,DAL2    ,DAL3    ,SFNL    ,DEP2    ,      &
		    KMESPC  ,IMATDA  ,IMATRA  ,FACHFR  ,PLNL4S  ,      &
		    IDCMIN  ,IDCMAX  ,IG      ,CGO     ,WWSWG   ) 
!       print*,'after calling SWSNL2'
!
     ELSE IF(IQUAD == 3)THEN
!
!      *** fully explicit calculation of the 4 wave-wave inter-  ***
!      *** actions for the full circle (1 -> MDC). An additional ***
!      *** array is required in which the values are stored prior***
!      *** to every iteration                                    ***
!
       IF(ITER == 1)THEN
!
!        *** calculate the interactions every sweep in each grid ***
!        *** point for the first iteration to ensure stable      ***
!        *** behaviour of the model                              ***
!
         CALL SWSNL3 (WWINT   ,WWAWG   ,UE      ,SA1     ,SA2     ,    &
	              SPCSIG  ,SNLC1   ,DAL1    ,DAL2    ,DAL3    ,    &
		      SFNL    ,DEP2    ,KMESPC  ,MEMNL4  ,FACHFR  )  
!
       ELSE IF(ITER > 1)THEN
!
         CALL SWSNL3 (WWINT   ,WWAWG   ,UE      ,SA1     ,SA2     ,    &
	              SPCSIG  ,SNLC1   ,DAL1    ,DAL2    ,DAL3    ,    &
		      SFNL    ,DEP2    ,KMESPC  ,MEMNL4  ,FACHFR  ) 
!
       ENDIF
!print*,'after SWSNL3'
!
!      *** Get source term value of additional array for the bin   ***
!      *** that fall within a sweep and store in right hand vector ***
!
       CALL FILNL3 (IDCMIN  ,IDCMAX  ,IMATRA  ,IMATDA  ,           &
                    MEMNL4  ,PLNL4S  ,ISSTOP  ,IG                  ) 

     ELSE IF(IQUAD == 4)THEN                                       

!      Multiple DIA according to Hashimoto (1999)                      

       IF(ITER >= 1)THEN
         DO IDIA=1,MDIA                                                
           PQUAD(1) = LAMBDA(IDIA)                                     
           CALL FAC4WW (XIS   ,SNLC1 ,DAL1  ,DAL2  ,DAL3  ,SPCSIG,     &
	                WWINT ,WWAWG ,WWSWG )           
           FACHFR = 1. / XIS ** PWTAIL(1)                              
           CALL RANGE4 (WWINT ,IDDLOW,IDDTOP )                         
           CALL SWSNL4 (WWINT   ,WWAWG   ,SPCSIG  ,SNLC1   ,           &
		        DAL1    ,DAL2    ,DAL3    ,DEP2    ,           &
			KMESPC  ,MEMNL4  ,FACHFR  ,IDIA    ,           &
			ITER )                                
         END DO
       ENDIF

!      Fill the matrix per sweep even though the quadruplets are calculated
!      only once per iteration

       CALL FILNL3 (IDCMIN  ,IDCMAX  ,IMATRA  ,IMATDA  ,           &
	            MEMNL4  ,PLNL4S  ,ISSTOP  ,IG                  )     
!
     ELSE IF(IQUAD == 6)THEN                                      

!      Calculate the quadruplets using the shallow water expression    
!      of Hashimoto, Komatsu and Masuda                                

!      Avoid calculation of the quadruplets in more than one sweep     

       IF((ITER >= 1) .AND.(.NOT.ONED))THEN               

         N2   = INT(MDC/2)+1                                           
         DT2  = DDIR/2.                                                
         DQ   = LOG((SHIG/SLOW)**(1./MSC))                             
         DQ2  = DQ/2.                                                  

!        The action density array should be extended to LMAX using     
!        a diagnostic tail.                                            

!        RIAM and SWAN use swapped indices for frequencies and         
!        directions                                                    

         LMAX = NINT(LOG(2.*SHIG/SLOW)/LOG(SHIG/SLOW)*MSC)             

         ALLOCATE(ACRIAM(LMAX,MDC),SNLRIAM(LMAX,MDC),SIGRIAM(LMAX))    
	 DO IS = 1, MSC                                                
           SIGRIAM(IS)=SPCSIG(IS)                                      
           DO ID = 1, MDC                                              
              ACRIAM(IS,ID) = AC2(ID,IS,KCGRD(1))                      
           ENDDO                                                       
         ENDDO                                                         

!        We need to multiply the action density with rho*grav in order 
!        to get action based on true energy and with Cg / k to get     
!        an action density spectrum in terms of wavenumber vector in
!        stead of frequency and direction. Wavenumber vector is used   
!        for the action density in the exact computations of Hashimoto.

         DO IS = 1, MSC                                                
           DO ID = 1, MDC                                              
             ACRIAM(IS,ID) = ACRIAM(IS,ID) * RHO_W * GRAV_W *            &
	                     CGO(IS,1) / KWAVE(IS,1)                  
           ENDDO                                                       
         ENDDO                                                         

!        extend the range from MSC to LMAX                             
         DO IS = MSC+1, LMAX                                           
           SIGRIAM(IS)=SIGRIAM(IS-1)*(SPCSIG(MSC)/SPCSIG(MSC-1))      
           DO ID = 1, MDC                                             
             ACRIAM(IS,ID) = ACRIAM(MSC,ID) *                        &
                             (SIGRIAM(MSC)/SIGRIAM(IS))**(PWTAIL(3)+3)      
           ENDDO                                                      
         ENDDO                                                         
         CALL RIAM_SLW(LMAX     ,MDC      ,N2       ,GRAV_W     ,      &
	               DEP2(IG) ,DQ       ,DQ2      ,DDIR     ,      &
		       DT2      ,SIGRIAM  ,RHO_W      ,ACRIAM   ,      &
		       SNLRIAM  ,1        )   

!        divide SNL4 by rho and grav to get energy source term         
!        formulated in variance density                                
!        divide by sigma to get action density source term             

         DO ID = 1, MDC
           DO IS = 1, MSC                                              
             MEMNL4(ID,IS,KCGRD(1)) = SNLRIAM(IS,ID) /                    &
	                         (SPCSIG(IS) * RHO_W * GRAV_W)       
           ENDDO                                                       
         ENDDO                                                         

         DEALLOCATE(ACRIAM,SNLRIAM,SIGRIAM)                            

       ENDIF                                                           

!       IF(ITEST >= 30)THEN                                           
!         WRITE (PRTEST,*) '+SOURCE: IX, IY: ',               &
!                          IX, IY                               
!       ENDIF                                                           
!       IF(TESTFL .AND. ITEST >= 100)THEN                               
!         DO IS=1, MSC                                                  
!           DO ID = 1, MDC                                             
!             WRITE(PRINTF,9200) IS,ID,MEMNL4(ID,IS,KCGRD(1))         
! 9200        FORMAT(' SOURCE: IS ID MEMNL(): ',2I6,E12.4)            
!           ENDDO                                                      
!         ENDDO                                                         
!       ENDIF                                                           

       CALL FILNL3 (IDCMIN  ,IDCMAX  ,IMATRA  ,IMATDA  ,           &
                    MEMNL4  ,PLNL4S  ,ISSTOP  ,IG      )      

!       IF (ITEST >= 30)THEN                                           
!         WRITE (PRTEST,*) '+SOURCE: ITER, IQUAD: ',          & 
!                          ITER, IQUAD                          
!       ENDIF                                                           
!
     ELSE IF(IQUAD == 8)THEN                                      
!
!      --- fully explicit calculation of the 4 wave-wave
!          interactions for the full circle. The interactions
!          in neighbouring bins are interpolated in piecewise
!          constant manner. An additional array is required in
!          which the values are stored prior to every iteration
!
       IF(ITER == 1)THEN
!
!        *** calculate the interactions every sweep in each grid ***
!        *** point for the first iteration to ensure stable      ***
!        *** behaviour of the model                              ***
!
         CALL SWSNL8 (WWINT   ,UE      ,SA1     ,SA2     ,SPCSIG  ,   &
	              SNLC1   ,DAL1    ,DAL2    ,DAL3    ,SFNL    ,   &
	              DEP2    ,KMESPC  ,MEMNL4  ,FACHFR  )
!
       ELSE IF(ITER > 1)THEN
!
         CALL SWSNL8 (WWINT   ,UE      ,SA1     ,SA2     ,SPCSIG  ,   &
	              SNLC1   ,DAL1    ,DAL2    ,DAL3    ,SFNL    ,   &
	              DEP2    ,KMESPC  ,MEMNL4  ,FACHFR  )
!
       ENDIF
!
!      *** get source term value of additional array for the bin   ***
!      *** that fall within a sweep and store in right hand vector ***
!
       CALL FILNL3 (IDCMIN  ,IDCMAX  ,IMATRA  ,IMATDA  ,          &
	            MEMNL4  ,PLNL4S  ,ISSTOP  ,IG      )  
!
     ELSE IF((IQUAD == 51).OR.(IQUAD == 52).OR.(IQUAD == 53))THEN     
!
!      Calculate the quadruplets using the XNL interface of            
!      G. van Vledder                                                  
!
!      Avoid calculation of the quadruplets in more than one sweep     
!
       IF((ITER >= 1) .AND. (.NOT.ONED))THEN   
!
!         IF(ITEST >= 30) WRITE(PRINTF,'(A,3I6,F12.2)')               &
!	        'SOURCE XNL: iter, kcgrd(1), iquad, depth:',  &
!	        ITER, KCGRD(1), IQUAD, DEP2(KCGRD(1))            

!JQI         CALL SWINTFXNL(AC2,SPCSIG,SPCDIR,MDC,MSC,MCGRD,              &
!JQI	                DEP2,IQUAD,MEMNL4,KCGRD,ICMAX,IQERR)        
!
       ENDIF                                                           
!
!       IF(ITEST >= 30)THEN                                           
!         WRITE (PRTEST,*) '+SOURCE: IX, IY: ',                &
!                          IX, IY                             
!       ENDIF                                                           
!       IF(TESTFL .AND. ITEST >= 100)THEN                               
!         DO IS=1, MSC                                                 
!           DO ID = 1, MDC                                            
!             WRITE(PRINTF,9300) IS,ID,MEMNL4(ID,IS,KCGRD(1))         
! 9300        FORMAT(' SOURCE: IS ID MEMNL(): ',2I6,E12.4)            
!           ENDDO                                                     
!         ENDDO                                                        
!       ENDIF                                                           
!
       CALL FILNL3 (IDCMIN  ,IDCMAX  ,IMATRA  ,IMATDA  ,           &
                    MEMNL4  ,PLNL4S  ,ISSTOP  ,IG      )   
!
!       IF(ITEST >= 30)THEN                                           
!         WRITE (PRTEST,*) '+SOURCE: ITER, IQUAD, IQERR: ',     &
!	                   ITER, IQUAD, IQERR                  
!       ENDIF                                                           
!
     ENDIF
   ENDIF
!
!  --- add contribution due to reflection of obstacles
!
   IF(NUMOBS /= 0)THEN                                               
!JQI     DO IS = 1, MSC                                                   
!JQI       DO ID = 1, MDC                                                
!JQI         IMATRA(ID,IS) = IMATRA(ID,IS) + REFLSO(ID,IS)          
!JQI       END DO                                                        
!JQI     END DO                                                           
   END IF                                                             

   RETURN
   END SUBROUTINE SOURCE
 
!
!************************************************************************
!                                                                       *
      SUBROUTINE PHILIM(AC2,AC2OLD,CGO,KWAVE,SPCSIG,ANYBIN,ISLMIN,NFLIM, &
                        QB_LOC)

!     (This subroutine has not been used and tested yet)
!                                                                       *
!************************************************************************
!
      USE SWCOMM3 
      USE ALL_VARS, ONLY : MT                                                        
!
      IMPLICIT NONE
!
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering                              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmers: R.C. Ris, N. Booij,                          |
!     |              IJ.G. Haagsma, A.T.M.M. Kieftenburg,         |
!     |              M. Zijlema, E.E. Kriezi,                     |
!     |              R. Padilla-Hernandez, L.H. Holthuijsen       |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 2004-2005  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. Authors
!
!     30.82: IJsbrand Haagsma
!     40.16: IJsbrand Haagsma
!     40.23: Marcel Zijlema
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     30.82, Feb. 99: New subroutine
!     40.16, Dec. 01: Implemented limiter switch
!     40.23, Aug. 02: Store number of frequency use of limiter
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Limits the change in action density between two iterations to a
!     certain percentage of the (directionally independent) Phillips
!     equilibrium level
!
!  3. Method
!
!     The maximum change of energy density per bin is related to
!     the (directionally independent) Phillips equilibrium level.
!     This change is estimated as:
!
!     |D E(s)| = factor * alpha_PM * g^2 / (s^5)
!
!     in which the Phillips' constant for a Pierson-Moskowitz
!     spectrum (alpha_PM) is taken to be 0.0081. Note that this
!     is a measure for a 1D spectrum. In SWAN, factor = 0.1
!     (stored in PNUMS(20)).
!
!     In terms of action density, we have
!
!     |D N(s)| = factor * alpha_PM * g^2 / (s^6)
!
!     Expressing in wave number k, this becomes with a deep
!     water approach of s^2 = gk:
!
!     |D N(s)| = factor * alpha_PM / (s^2 k^2)
!
!     Furthermore, with s = 2*k*c_g (deep water), we finally
!     have (Ris, 1997, p.36):
!
!                          alpha_PM
!     |D N(s,t)| = factor -----------
!                         2 s k^3 c_g
!
!     In cases where waves are breaking the dissipation of energy
!     is not limited. This is assumed to be the case when the
!     fraction of breaking waves Qb is more than 1.e-5.
!
!  4. Argument variables
!
      LOGICAL ANYBIN(MDC,MSC)
!
!     ISLMIN: Lowest sigma-index occured in applying limiter
!     NFLIM : Number of frequency use of limiter
!     QB_LOC: Local value of Qb (fraction of breaking waves)
!
      INTEGER ISLMIN(MT), NFLIM(MT)
      REAL    QB_LOC
      REAL    AC2(MDC,MSC,0:MT)
      REAL    AC2OLD(MDC,MSC)
      REAL    CGO(MSC,ICMAX)
      REAL    KWAVE(MSC,ICMAX)
      REAL    SPCSIG(MSC)
!
!  6. Local variables
!
!     ID    : Counter for directional (theta) space
!     IS    : Counter for frequency (sigma) space
!
      INTEGER ID,IS
!
!     DAC2MX: Maximum deviation of action density AC2 between iterations
!
      REAL    DAC2MX
!
! 13. Source text
!
      IF (MSC.GT.3) THEN
         DO IS=1,MSC
            DAC2MX=ABS((PNUMS(20)*0.0081)/                         &
	               (2.*SPCSIG(IS)*(KWAVE(IS,1)**3)*CGO(IS,1)))
            DO ID=1,MDC
               IF (ANYBIN(ID,IS) .AND.                             &
	           AC2(ID,IS,KCGRD(1)).GT.AC2OLD(ID,IS)+DAC2MX) THEN
                  AC2(ID,IS,KCGRD(1))=AC2OLD(ID,IS)+DAC2MX
                  NFLIM(KCGRD(1)) = NFLIM(KCGRD(1)) + 1
                  ISLMIN(KCGRD(1)) = MIN(IS,ISLMIN(KCGRD(1)))
               END IF
            END DO
            IF (QB_LOC.LT.PNUMS(28)) THEN                                
               DO ID=1,MDC
                  IF (ANYBIN(ID,IS) .AND.                          &
		     AC2(ID,IS,KCGRD(1)).LT.AC2OLD(ID,IS)-DAC2MX) THEN
                     AC2(ID,IS,KCGRD(1))=AC2OLD(ID,IS)-DAC2MX
                     NFLIM(KCGRD(1)) = NFLIM(KCGRD(1)) + 1
                     ISLMIN(KCGRD(1)) = MIN(IS,ISLMIN(KCGRD(1)))
                  END IF
               END DO
            END IF
         END DO
      END IF
      RETURN
      END SUBROUTINE PHILIM
!****************************************************************
!
      SUBROUTINE RESCALE (AC2, ISSTOP, IDCMIN, IDCMAX, NRSCAL)

!     (This subroutine has not been used and tested yet)
!
!****************************************************************

      USE SWCOMM3                                                         
      USE SWCOMM4                                                         
      USE OCPCOMM4                                                        
      USE M_PARALL 
      USE ALL_VARS, ONLY : MT                                                       

      IMPLICIT NONE
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering                              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmers: R.C. Ris, N. Booij,                          |
!     |              IJ.G. Haagsma, A.T.M.M. Kieftenburg,         |
!     |              M. Zijlema, E.E. Kriezi,                     |
!     |              R. Padilla-Hernandez, L.H. Holthuijsen       |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 2004-2005  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!  0. Authors
!
!     40.00: Nico Booij
!     40.23: Marcel Zijlema
!     40.30: Marcel Zijlema
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     40.00, Feb. 99: New subroutine (software moved from subroutines
!                     SOLBAND, SOLMT1 and SOLMAT
!     40.23, Aug. 02: Store number of frequency use of rescaling
!     40.30, Mar. 03: correcting indices of test point with offsets MXF, MYF
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Remove negative values from a computed action density spectrum
!
!  3. Method
!
!     Make negative action densities 0 at the expense of other action densities
!     for the frequency
!
!  4. Argument variables
!
!     AC2         action densities
!
      REAL        AC2(MDC,MSC,0:MT)
!
!     ISSTOP      maximum frequency counter in this sweep
!
      INTEGER     ISSTOP
!
!     IDCMIN      Integer array containing minimum counter of directions
!     IDCMAX      Integer array containing maximum counter
!     NRSCAL      Number of frequency use of rescaling
!
      INTEGER     IDCMIN(MSC), IDCMAX(MSC)
      INTEGER     NRSCAL(MT)
!
!  7. Common blocks used
!
!
!     5. SUBROUTINES CALLING
!
!        SWOMPU
!
!     6. SUBROUTINES USED
!
!        NONE
!
!     7. ERROR MESSAGES
!
!        ---
!
!     8. REMARKS
!
!        ---
!
!     9. STRUCTURE
!
!   -------------------------------------------------------------
!   For all frequencies do
!       Make ATOT equal to integral of action density over direction
!       Make ATOTP equal to integral of positive action density
!       Determine FACTOR
!       If negative values do occur
!       Then for all directions do
!            If action density is negative
!            Then make action density =0
!            Else multiply action density by FACTOR
!   ------------------------------------------------------------
!
!     10. SOURCE
!
!************************************************************************
!
!         local variables
!
!         IS         counter of frequency
!         ID         counter of direction
!         IDDUM      uncorrected counter of direction
!
      INTEGER  IS      ,ID      ,IDDUM,   IENT
!
!         ATOT       integral of action density for one frequency
!         ATOTP      integral of positive action density for one frequency
!         FACTOR
!
      REAL     ATOT    ,ATOTP   ,FACTOR
!
!         NEGVAL      if True, there are negative values in the spectrum
!
      LOGICAL  NEGVAL
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'RESCALE')
!
!     *** if negative action density occur rescale with a factor ***
!     *** only the sector computed is rescaled !!                ***
!
      DO 180 IS = 1 , ISSTOP
        ATOT   = 0.
        ATOTP  = 0.
        FACTOR = 0.
        NEGVAL = .FALSE.
        DO 160 IDDUM = IDCMIN(IS), IDCMAX(IS)
          ID = MOD ( IDDUM - 1 + MDC , MDC ) + 1
          ATOT = ATOT + AC2(ID,IS,KCGRD(1))
          IF ( AC2(ID,IS,KCGRD(1)) .LT. 0. ) THEN
            NRSCAL(KCGRD(1)) = NRSCAL(KCGRD(1)) + 1
            NEGVAL = .TRUE.
          ELSE
            ATOTP = ATOTP + AC2(ID,IS,KCGRD(1))
          END IF
 160    CONTINUE
        IF (NEGVAL) THEN
          IF ( ATOTP .LT. 1.E-15 ) ATOTP = 1.E-15
          FACTOR = ATOT / ATOTP
!
!         *** rescale ***
!
          DO 170 IDDUM = IDCMIN(IS), IDCMAX(IS)
            ID = MOD ( IDDUM - 1 + MDC , MDC ) + 1
            IF ( AC2(ID,IS,KCGRD(1)) .LT. 0.) THEN
              AC2(ID,IS,KCGRD(1)) = 0.
            END IF
            IF ( FACTOR .GE. 0. ) THEN
              AC2(ID,IS,KCGRD(1)) = FACTOR * AC2(ID,IS,KCGRD(1))
            ENDIF
 170      CONTINUE
!
          IF ( ITEST .GE. 120 .AND. TESTFL )                       &
	    WRITE (PRINTF, 171) IXCGRD(1)+MXF-2, IYCGRD(1)+MYF-2, IS,   &
	        FACTOR , ATOT, ATOTP
 171      FORMAT(' Rescale in Point, Isig, Factor, ATOT, ATOTP:',     &
                 3I4, 3(1X,E11.4))
        ENDIF
 180  CONTINUE
      RETURN
      END SUBROUTINE RESCALE
!****************************************************************
!
      SUBROUTINE SWSIP ( AC2   , IMATDA, IMATRA, IMATLA, IMATUA,     &
                         IMAT5L, IMAT6U, AC2OLD, REPS  , MAXIT ,     &
			 IAMOUT, INOCNV, IDDLOW, IDDTOP, ISSTOP,     &
			 IDCMIN, IDCMAX )                        

!     (This subroutine has not been used and tested yet)
!
!****************************************************************
!
      USE SWCOMM1                                                         
      USE SWCOMM3                                                         
      USE SWCOMM4                                                         
      USE OCPCOMM4                                                        
      USE M_PARALL  
      USE ALL_VARS, ONLY : MT                                                      

      IMPLICIT NONE
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: M. Zijlema                                    |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 2004-2005  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. Authors
!
!     40.23: Marcel Zijlema
!     40.30: Marcel Zijlema
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     40.23, Oct. 02: New subroutine
!     40.30, Mar. 03: introduction distributed-memory approach using MPI
!     40.41, Mar. 04: parameter ALFA set to 0.0, extra test output
!                     and some corrections if SECTOR=0
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Solves penta-diagonal system of equations in
!     spectral space by means of Stone's SIP solver
!
!  4. Argument variables
!
!     AC2         action density
!     AC2OLD      action density at previous iteration
!     IAMOUT      control parameter indicating the amount of
!                 output required
!                 0: no output
!                 1: only fatal errors will be printed
!                 2: gives output concerning the iteration process
!                 3: additional information about the iteration
!                    is printed
!     IDCMAX      maximum counter in directional space
!     IDCMIN      minimum counter in directional space
!     IDDLOW      minimum direction that is propagated within a sweep
!     IDDTOP      maximum direction that is propagated within a sweep
!     IMAT5L      coefficients of lower diagonal in sigma-space
!     IMAT6U      coefficients of upper diagonal in sigma-space
!     IMATDA      coefficients of main diagonal
!     IMATLA      coefficients of lower diagonal in theta-space
!     IMATUA      coefficients of upper diagonal in theta-space
!     IMATRA      right-hand side
!     INOCNV      integer indicating number of grid points in which
!                 solver does not converged
!     ISSTOP      maximum frequency counter in a sweep
!     MAXIT       the maximum number of iterations to be performed in
!                 the linear solver
!     REPS        accuracy with respect to the right-hand side used
!                 in the following termination criterion:
!
!                 ||b-Ax || < reps*||b||
!                       k
!
      INTEGER IAMOUT, INOCNV, IDDLOW, IDDTOP, ISSTOP, MAXIT
      INTEGER IDCMIN(MSC), IDCMAX(MSC)
      REAL    REPS
      REAL    AC2(MDC,MSC,0:MT),                         &
              IMATDA(MDC,MSC), IMATRA(MDC,MSC),           &
	      IMAT5L(MDC,MSC), IMAT6U(MDC,MSC),           &
	      IMATLA(MDC,MSC), IMATUA(MDC,MSC),           &
	      AC2OLD(MDC,MSC)
!
!  5. Parameter variables
!
!     ALFA        relaxation parameter used in the SIP solver
!     SMALL :     a small number
!
      REAL    ALFA, SMALL
      PARAMETER (ALFA=0.0,SMALL=1.E-15)                                   
!
!  6. Local variables
!
!     BNORM :     2-norm of right-hand side vector
!     CMAT5L:     coefficients of lower diagonal in sigma-space
!                 obtained by an incomplete lower-upper factorization
!     CMAT6U:     coefficients of upper diagonal in sigma-space
!                 obtained by an incomplete lower-upper factorization
!     CMATDA:     coefficients of main diagonal obtained by an
!                 incomplete lower-upper factorization
!     CMATLA:     coefficients of lower diagonal in theta-space
!                 obtained by an incomplete lower-upper factorization
!     CMATUA:     coefficients of upper diagonal in theta-space
!                 obtained by an incomplete lower-upper factorization
!     EPSLIN:     required accuracy in the linear solver
!     ICONV :     indicator for convergence (1=yes, 0=no)
!     ID    :     loop counter
!     IDDL  :     minimum counter in theta-space of modulo MDC
!     IDDT  :     maximum counter in theta-space of modulo MDC
!     IDDUM :     loop counter
!     IDM   :     index of point ID-1
!     IDMAX :     local array of maximum counter in theta-space
!     IDMIN :     local array of minimum counter in theta-space
!     IDP   :     index of point ID+1
!     IENT  :     number of entries
!     IS    :     loop counter
!     ISM   :     index of point IS-1
!     ISP   :     index of point IS+1
!     IT    :     iteration count
!     LOPERI:     auxiliary vector meant for computation in
!                 periodic theta-space
!     P1    :     auxiliary factor
!     P2    :     auxiliary factor
!     P3    :     auxiliary factor
!     RES   :     the residual vector
!     RNORM :     2-norm of residual vector
!     RNRM0 :     2-norm of initial residual vector
!     UEPS  :     minimal accuracy based on machine precision
!     UPPERI:     auxiliary vector meant for computation in
!                 periodic theta-space
!
      INTEGER ICONV, ID, IDDL, IDDT, IDDUM, IDM, IDP, IENT,     &
              IS, ISM, ISP, IT
      INTEGER IDMIN(MSC), IDMAX(MSC)
      REAL    BNORM, EPSLIN, P1, P2, P3, RNORM, RNRM0, UEPS
      REAL    RES(MDC,MSC)   , CMATDA(MDC,MSC),                 &
              CMAT5L(MDC,MSC), CMAT6U(MDC,MSC),                 &
	      CMATLA(MDC,MSC), CMATUA(MDC,MSC),                 &
	      LOPERI(MSC)    , UPPERI(MSC)
!
!  7. Common blocks used
!
!
!  8. Subroutines used
!
!     STRACE           Tracing routine for debugging
!
!  9. Subroutines calling
!
!     SWOMPU  (in SWANCOM1)
!
! 12. Structure
!
!     The system of equations is solved using an incomplete
!     factorization technique called Strongly Implicit Procedure
!     (SIP) as described in
!
!     H.L. Stone
!     Iterative solution of implicit approximations of
!     multidimensional partial differential equations
!     SIAM J. of Numer. Anal., vol. 5, 530-558, 1968
!
!     This method constructs an incomplete lower-upper factorization
!     that has the same sparsity as the original matrix. Hereby, a
!     parameter alfa is used, which should be 0.0 in case of SWAN         40.41
!     (when alfa > 0.95, the method may diverge).
!
!     Afterward, the resulting system is solved in an iterative manner
!     by forward and backward substitutions.
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWSIP')

!     --- initialize arrays

      RES    = 0.
      CMATDA = 0.
      CMATLA = 0.
      CMATUA = 0.
      CMAT5L = 0.
      CMAT6U = 0.
      LOPERI = 0.
      UPPERI = 0.

!     --- in case of periodicity in theta-space, store values
!         of matrix coefficients corresponding to left bottom and
!         right top

      DO IS = 1, ISSTOP
         IF ( IDCMIN(IS).EQ.1 .AND. IDCMAX(IS).EQ.MDC ) THEN
           UPPERI(IS) = IMATLA(  1,IS)
           LOPERI(IS) = IMATUA(MDC,IS)
         END IF
      END DO

!     --- when no bins fall within the sweep, i.e. SECTOR = 0,
!         reset the bounds of sector as 1..MDC (routine SOLPRE
!         has clear the rows in the matrix that do not belong
!         to the sweep)

      DO IS = 1, ISSTOP
         IF ( IDCMIN(IS).LE.IDCMAX(IS) ) THEN
            IDMIN(IS) = IDCMIN(IS)
            IDMAX(IS) = IDCMAX(IS)
         ELSE
            IDMIN(IS) = 1
            IDMAX(IS) = MDC
         END IF
      END DO

      IT    = 0
      ICONV = 0

!     --- construct L and U matrices (stored in CMAT[xx])

      BNORM = 0.

      IS    = 1
      IDDUM = IDMIN(IS)
      ID    = MOD ( IDDUM - 1 + MDC , MDC ) + 1

      CMAT5L(ID,IS) = IMAT5L(ID,IS)
      CMATLA(ID,IS) = IMATLA(ID,IS)
      CMATDA(ID,IS) = 1./(IMATDA(ID,IS)+SMALL)
      CMAT6U(ID,IS) = IMAT6U(ID,IS)*CMATDA(ID,IS)
      CMATUA(ID,IS) = IMATUA(ID,IS)*CMATDA(ID,IS)
      BNORM = BNORM + IMATRA(ID,IS)*IMATRA(ID,IS)

      DO IDDUM = IDMIN(IS)+1, IDMAX(IS)
         ID  = MOD ( IDDUM - 1 + MDC , MDC ) + 1
         IDM = MOD ( IDDUM - 2 + MDC , MDC ) + 1

         P2 = ALFA*CMAT6U(IDM,IS)
         CMAT5L(ID,IS) = IMAT5L(ID,IS)
         CMATLA(ID,IS) = IMATLA(ID,IS)/(1.+P2)

         P2 = P2*CMATLA(ID,IS)
         P3 = IMATDA(ID,IS) + P2                        &
	     -CMATLA(ID,IS)*CMATUA(IDM,IS )             &
	     +SMALL
         CMATDA(ID,IS) = 1./P3
         CMAT6U(ID,IS) = (IMAT6U(ID,IS)-P2)*CMATDA(ID,IS)
         CMATUA(ID,IS) =      IMATUA(ID,IS)*CMATDA(ID,IS)
         BNORM = BNORM + IMATRA(ID,IS)*IMATRA(ID,IS)
      END DO

      DO IS = 2, ISSTOP
         ISM = IS - 1

         IDDUM = IDMIN(IS)
         ID  = MOD ( IDDUM - 1 + MDC , MDC ) + 1

         P1 = ALFA*CMATUA(ID,ISM)
         CMAT5L(ID,IS) = IMAT5L(ID,IS)/(1.+P1)
         CMATLA(ID,IS) = IMATLA(ID,IS)
         P1 = P1*CMAT5L(ID,IS)
         P3 = IMATDA(ID,IS) + P1                        &
	     -CMAT5L(ID,IS)*CMAT6U(ID,ISM)              &
	     +SMALL
         CMATDA(ID,IS) = 1./P3
         CMAT6U(ID,IS) =      IMAT6U(ID,IS)*CMATDA(ID,IS)
         CMATUA(ID,IS) = (IMATUA(ID,IS)-P1)*CMATDA(ID,IS)
         BNORM = BNORM + IMATRA(ID,IS)*IMATRA(ID,IS)

         DO IDDUM = IDMIN(IS)+1, IDMAX(IS)
            ID  = MOD ( IDDUM - 1 + MDC , MDC ) + 1
            IDM = MOD ( IDDUM - 2 + MDC , MDC ) + 1

            P1 = ALFA*CMATUA(ID ,ISM)
            P2 = ALFA*CMAT6U(IDM,IS )
            CMAT5L(ID,IS) = IMAT5L(ID,IS)/(1.+P1)
            CMATLA(ID,IS) = IMATLA(ID,IS)/(1.+P2)
            P1 = P1*CMAT5L(ID,IS)
            P2 = P2*CMATLA(ID,IS)
            P3 = IMATDA(ID,IS) + P1 + P2                &
	        -CMAT5L(ID,IS)*CMAT6U(ID ,ISM)          &
		-CMATLA(ID,IS)*CMATUA(IDM,IS )          &
		+SMALL
            CMATDA(ID,IS) = 1./P3
            CMAT6U(ID,IS) = (IMAT6U(ID,IS)-P2)*CMATDA(ID,IS)
            CMATUA(ID,IS) = (IMATUA(ID,IS)-P1)*CMATDA(ID,IS)
            BNORM = BNORM + IMATRA(ID,IS)*IMATRA(ID,IS)
         END DO
      END DO
      BNORM = SQRT(BNORM)

      EPSLIN = REPS*BNORM
      UEPS   = 1000.*UNDFLW*BNORM
      IF ( EPSLIN.LT.UEPS .AND. BNORM.GT.0. ) THEN
         IF ( IAMOUT.GE.1 ) THEN
            WRITE (PRINTF,'(A)')                            &
	          ' ++ SWSIP: the required accuracy is too small'
            WRITE (PRINTF,*)                                &
	          '           required accuracy    = ',EPSLIN
            WRITE (PRINTF,*)                                &
	          '           appropriate accuracy = ',UEPS
         END IF
         EPSLIN = UEPS
      END IF

!     --- solve the system by forward and backward substitutions
!         in an iterative manner

 10   IF ( ICONV.EQ.0 .AND. IT.LT.MAXIT ) THEN

         IT    = IT + 1
         ICONV = 1
         RNORM = 0.

         IS  = 1
         ISP = IS + 1

         IDDUM = IDMIN(IS)
         ID    = MOD ( IDDUM - 1 + MDC , MDC ) + 1
         IDP   = MOD ( IDDUM     + MDC , MDC ) + 1
         IDDT  = MOD ( IDMAX(IS) - 1 + MDC , MDC ) + 1

         RES(ID,IS) = IMATRA(ID,IS)                           &
	              -IMATDA(ID,IS)*AC2(ID ,IS ,KCGRD(1))    &
		      -IMAT6U(ID,IS)*AC2(ID ,ISP,KCGRD(1))    &
		      -IMATUA(ID,IS)*AC2(IDP,IS ,KCGRD(1))    &
		      -UPPERI(IS)*AC2(IDDT,IS,KCGRD(1))
         RNORM = RNORM + RES(ID,IS)*RES(ID,IS)
         RES(ID,IS) = RES(ID,IS)*CMATDA(ID,IS)

         DO IDDUM = IDMIN(IS)+1, IDMAX(IS)-1
            ID  = MOD ( IDDUM - 1 + MDC , MDC ) + 1
            IDM = MOD ( IDDUM - 2 + MDC , MDC ) + 1
            IDP = MOD ( IDDUM     + MDC , MDC ) + 1

            RES(ID,IS) = IMATRA(ID,IS)                        &
	                 -IMATDA(ID,IS)*AC2(ID ,IS ,KCGRD(1)) &
			 -IMAT6U(ID,IS)*AC2(ID ,ISP,KCGRD(1)) &
			 -IMATLA(ID,IS)*AC2(IDM,IS ,KCGRD(1)) &
			 -IMATUA(ID,IS)*AC2(IDP,IS ,KCGRD(1))
            RNORM = RNORM + RES(ID,IS)*RES(ID,IS)
            RES(ID,IS) = (RES(ID,IS) - CMATLA(ID,IS)*RES(IDM,IS))*   &
	                  CMATDA(ID,IS)
         END DO

         IDDUM = IDMAX(IS)
         ID    = MOD ( IDDUM - 1 + MDC , MDC ) + 1
         IDM   = MOD ( IDDUM - 2 + MDC , MDC ) + 1
         IDDL  = MOD ( IDMIN(IS) - 1 + MDC , MDC ) + 1

         RES(ID,IS) = IMATRA(ID,IS)                           &
	              -IMATDA(ID,IS)*AC2(ID ,IS ,KCGRD(1))    &
		      -IMAT6U(ID,IS)*AC2(ID ,ISP,KCGRD(1))    &
		      -IMATLA(ID,IS)*AC2(IDM,IS ,KCGRD(1))    &
		      -LOPERI(IS)*AC2(IDDL,IS,KCGRD(1))
         RNORM = RNORM + RES(ID,IS)*RES(ID,IS)
         RES(ID,IS) = (RES(ID,IS) - CMATLA(ID,IS)*RES(IDM,IS))*  &
	               CMATDA(ID,IS)

         DO IS = 2, ISSTOP-1
            ISM = IS - 1
            ISP = IS + 1

            IDDUM = IDMIN(IS)
            ID    = MOD ( IDDUM - 1 + MDC , MDC ) + 1
            IDP   = MOD ( IDDUM     + MDC , MDC ) + 1
            IDDT  = MOD ( IDMAX(IS) - 1 + MDC , MDC ) + 1

            RES(ID,IS) = IMATRA(ID,IS)                        &
	                -IMATDA(ID,IS)*AC2(ID ,IS ,KCGRD(1))  &
			-IMAT5L(ID,IS)*AC2(ID ,ISM,KCGRD(1))  &
			-IMAT6U(ID,IS)*AC2(ID ,ISP,KCGRD(1))  &
			-IMATUA(ID,IS)*AC2(IDP,IS ,KCGRD(1))  &
			-UPPERI(IS)*AC2(IDDT,IS,KCGRD(1))
            RNORM = RNORM + RES(ID,IS)*RES(ID,IS)
            RES(ID,IS) = (RES(ID,IS) - CMAT5L(ID,IS)*RES(ID,ISM))* &
	                  CMATDA(ID,IS)

            DO IDDUM = IDMIN(IS)+1, IDMAX(IS)-1
               ID  = MOD ( IDDUM - 1 + MDC , MDC ) + 1
               IDM = MOD ( IDDUM - 2 + MDC , MDC ) + 1
               IDP = MOD ( IDDUM     + MDC , MDC ) + 1

               RES(ID,IS) = IMATRA(ID,IS)                         &
	                   -IMATDA(ID,IS)*AC2(ID ,IS ,KCGRD(1))   &
			   -IMAT5L(ID,IS)*AC2(ID ,ISM,KCGRD(1))   &
			   -IMAT6U(ID,IS)*AC2(ID ,ISP,KCGRD(1))   &
			   -IMATLA(ID,IS)*AC2(IDM,IS ,KCGRD(1))   &
			   -IMATUA(ID,IS)*AC2(IDP,IS ,KCGRD(1))
               RNORM = RNORM + RES(ID,IS)*RES(ID,IS)
               RES(ID,IS) = (RES(ID,IS) - CMAT5L(ID,IS)*RES(ID ,ISM)  &
	                   - CMATLA(ID,IS)*RES(IDM,IS ))*             &
			   CMATDA(ID,IS)
            END DO

            IDDUM = IDMAX(IS)
            ID    = MOD ( IDDUM - 1 + MDC , MDC ) + 1
            IDM   = MOD ( IDDUM - 2 + MDC , MDC ) + 1
            IDDL  = MOD ( IDMIN(IS) - 1 + MDC , MDC ) + 1

            RES(ID,IS) = IMATRA(ID,IS)                            &
	                -IMATDA(ID,IS)*AC2(ID ,IS ,KCGRD(1))      &
			-IMAT5L(ID,IS)*AC2(ID ,ISM,KCGRD(1))      &
			-IMAT6U(ID,IS)*AC2(ID ,ISP,KCGRD(1))      &
			-IMATLA(ID,IS)*AC2(IDM,IS ,KCGRD(1))      &
			-LOPERI(IS)*AC2(IDDL,IS,KCGRD(1))
            RNORM = RNORM + RES(ID,IS)*RES(ID,IS)
            RES(ID,IS) = (RES(ID,IS) - CMAT5L(ID,IS)*RES(ID ,ISM)   &
	                - CMATLA(ID,IS)*RES(IDM,IS ))*              &
			  CMATDA(ID,IS)

         END DO

         IS  = ISSTOP
         ISM = IS - 1

         IDDUM = IDMIN(IS)
         ID    = MOD ( IDDUM - 1 + MDC , MDC ) + 1
         IDP   = MOD ( IDDUM     + MDC , MDC ) + 1
         IDDT  = MOD ( IDMAX(IS) - 1 + MDC , MDC ) + 1

         RES(ID,IS) = IMATRA(ID,IS)                         &
	             -IMATDA(ID,IS)*AC2(ID ,IS ,KCGRD(1))   &
		     -IMAT5L(ID,IS)*AC2(ID ,ISM,KCGRD(1))   &
		     -IMATUA(ID,IS)*AC2(IDP,IS ,KCGRD(1))   &
		     -UPPERI(IS)*AC2(IDDT,IS,KCGRD(1))
         RNORM = RNORM + RES(ID,IS)*RES(ID,IS)
         RES(ID,IS) = (RES(ID,IS) - CMAT5L(ID,IS)*RES(ID,ISM))*   &
	               CMATDA(ID,IS)

         DO IDDUM = IDMIN(IS)+1, IDMAX(IS)-1
            ID  = MOD ( IDDUM - 1 + MDC , MDC ) + 1
            IDM = MOD ( IDDUM - 2 + MDC , MDC ) + 1
            IDP = MOD ( IDDUM     + MDC , MDC ) + 1

            RES(ID,IS) = IMATRA(ID,IS)                       &
	                -IMATDA(ID,IS)*AC2(ID ,IS ,KCGRD(1)) &
			-IMAT5L(ID,IS)*AC2(ID ,ISM,KCGRD(1)) &
			-IMATLA(ID,IS)*AC2(IDM,IS ,KCGRD(1)) &
			-IMATUA(ID,IS)*AC2(IDP,IS ,KCGRD(1))
            RNORM = RNORM + RES(ID,IS)*RES(ID,IS)
            RES(ID,IS) = (RES(ID,IS) - CMAT5L(ID,IS)*RES(ID ,ISM)  &
	                - CMATLA(ID,IS)*RES(IDM,IS ))*             &
			  CMATDA(ID,IS)
         END DO

         IDDUM = IDMAX(IS)
         ID    = MOD ( IDDUM - 1 + MDC , MDC ) + 1
         IDM   = MOD ( IDDUM - 2 + MDC , MDC ) + 1
         IDDL  = MOD ( IDMIN(IS) - 1 + MDC , MDC ) + 1

         RES(ID,IS) = IMATRA(ID,IS)                         &
	             -IMATDA(ID,IS)*AC2(ID ,IS ,KCGRD(1))   &
		     -IMAT5L(ID,IS)*AC2(ID ,ISM,KCGRD(1))   &
		     -IMATLA(ID,IS)*AC2(IDM,IS ,KCGRD(1))   &
		     -LOPERI(IS)*AC2(IDDL,IS,KCGRD(1))
         RNORM = RNORM + RES(ID,IS)*RES(ID,IS)
         RES(ID,IS) = (RES(ID,IS) - CMAT5L(ID,IS)*RES(ID ,ISM)     &
	             - CMATLA(ID,IS)*RES(IDM,IS ))*                &
		       CMATDA(ID,IS)

         IF ( RNORM.GT.1.E8 ) THEN
            IT = MAXIT + 1
            ICONV = 0
            GOTO 10
         END IF
         RNORM=SQRT(RNORM)
         IF ( IAMOUT.EQ.3 .AND. IT.EQ.1 ) RNRM0 = RNORM

         IF ( IAMOUT.EQ.2 ) THEN
            WRITE (PRINTF,'(A,I3,A,E13.6)')               &
	            ' ++ SWSIP: iter = ',IT,'    res = ',RNORM
         END IF

         IS    = ISSTOP
         IDDUM = IDMAX(IS)
         ID    = MOD ( IDDUM - 1 + MDC , MDC ) + 1

         AC2(ID,IS,KCGRD(1)) = AC2(ID,IS,KCGRD(1)) + RES(ID,IS)

         DO IDDUM = IDMAX(IS)-1, IDMIN(IS), -1
            ID  = MOD ( IDDUM - 1 + MDC , MDC ) + 1
            IDP = MOD ( IDDUM     + MDC , MDC ) + 1

            RES(ID,IS) = RES(ID,IS) - CMATUA(ID,IS)*RES(IDP,IS)
            AC2(ID,IS,KCGRD(1)) = AC2(ID,IS,KCGRD(1)) + RES(ID,IS)
         END DO

         DO IS = ISSTOP-1, 1, -1
            ISP = IS + 1

            IDDUM = IDMAX(IS)
            ID  = MOD ( IDDUM - 1 + MDC , MDC ) + 1

            RES(ID,IS) = RES(ID,IS) - CMAT6U(ID,IS)*RES(ID,ISP)
            AC2(ID,IS,KCGRD(1)) = AC2(ID,IS,KCGRD(1)) + RES(ID,IS)

            DO IDDUM = IDMAX(IS)-1, IDMIN(IS), -1
               ID  = MOD ( IDDUM - 1 + MDC , MDC ) + 1
               IDP = MOD ( IDDUM     + MDC , MDC ) + 1

               RES(ID,IS) = RES(ID,IS) - CMAT6U(ID,IS)*RES(ID ,ISP)   &
	                               - CMATUA(ID,IS)*RES(IDP,IS )
               AC2(ID,IS,KCGRD(1)) = AC2(ID,IS,KCGRD(1)) + RES(ID,IS)
            END DO
         END DO

         IF ( RNORM.GT.UNDFLW**2 .AND. RNORM.GT.EPSLIN ) ICONV = 0
         GOTO 10

      END IF

!     --- investigate the reason to stop

      IF ( ICONV.EQ.0 ) THEN
         AC2(:,:,KCGRD(1)) = AC2OLD(:,:)
         INOCNV = INOCNV + 1
      END IF
      IF ( ICONV.EQ.0 .AND. IAMOUT.GE.1 ) THEN
!         IF (ERRPTS.GT.0.AND.INODE.EQ.MASTER) THEN
!            WRITE(ERRPTS,100) IXCGRD(1)+MXF-1, IYCGRD(1)+MYF-1, 2
!         END IF
 100     FORMAT (I4,1X,I4,1X,I2)
         WRITE (PRINTF,'(A,I5,A,I5,A)')                     &
	      ' ++ SWSIP: no convergence in grid point (',  &
	       IXCGRD(1)+MXF-1,',',IYCGRD(1)+MYF-1,')'
         WRITE (PRINTF,'(A,I3)')                            &
	      '           total number of iterations     = ',IT
         WRITE (PRINTF,'(A,E13.6)')                         &
	      '           2-norm of the residual         = ',RNORM
         WRITE (PRINTF,'(A,E13.6)')                         &
	      '           required accuracy              = ',EPSLIN
      ELSE IF ( IAMOUT.EQ.3 ) THEN
         WRITE (PRINTF,'(A,E13.6)')                         &
	      ' ++ SWSIP: 2-norm of the initial residual = ',RNRM0
         WRITE (PRINTF,'(A,I3)')                            &
	      '           total number of iterations     = ',IT
         WRITE (PRINTF,'(A,E13.6)')                         &
	      '           2-norm of the residual         = ',RNORM
      END IF

!     --- test output

      IF ( TESTFL .AND. ITEST.GE.120 ) THEN
         WRITE(PRTEST,*)
         WRITE(PRTEST,*) '  Subroutine SWSIP'
         WRITE(PRTEST,*)
         WRITE(PRTEST,200) KCGRD(1), MDC, MSC
 200     FORMAT(' SWSIP : POINT MDC MSC              :',3I5)
         WRITE(PRTEST,250) IDDLOW, IDDTOP, ISSTOP
 250     FORMAT(' SWSIP : IDDLOW IDDTOP ISSTOP       :',3I4)
         WRITE(PRTEST,*)
         WRITE(PRTEST,*) ' coefficients of matrix and rhs  '
         WRITE(PRTEST,*)
         WRITE(PRTEST,'(A111)')                             &
	        ' IS ID         IMATLA         IMATDA'//    &
		'         IMATUA         IMATRA         IMAT5L'//   &
		'         IMAT6U            AC2'
         DO IDDUM = IDDLOW, IDDTOP
            ID = MOD ( IDDUM - 1 + MDC , MDC ) + 1
            DO IS = 1, ISSTOP
               WRITE(PRTEST,300) IS, ID,                            &
	                         IMATLA(ID,IS), IMATDA(ID,IS),      &
				 IMATUA(ID,IS), IMATRA(ID,IS),      &
				 IMAT5L(ID,IS), IMAT6U(ID,IS),      &
				 AC2(ID,IS,KCGRD(1))
 300           FORMAT(2I3,7E15.7)
            END DO
         END DO
         WRITE(PRTEST,*)
         WRITE(PRTEST,*)'IS ID      LPER          UPER '
         DO IDDUM = IDDLOW, IDDTOP
            ID = MOD ( IDDUM - 1 + MDC , MDC ) + 1
            IF ( ID.EQ.1 .OR. ID.EQ.MDC ) THEN
               DO IS = 1, ISSTOP
                  WRITE(PRTEST,350) IS, ID, LOPERI(IS), UPPERI(IS)
 350              FORMAT(2I3,2E15.7)
               END DO
            END IF
         END DO
      END IF

!     --- set matrix coefficients to zero

      IMATDA = 0.
      IMATRA = 0.
      IMATLA = 0.
      IMATUA = 0.
      IMAT5L = 0.
      IMAT6U = 0.

      RETURN
      END SUBROUTINE SWSIP
!****************************************************************
!
      SUBROUTINE SWSOR ( AC2   , IMATDA, IMATRA, IMATLA, IMATUA, &
                         IMAT5L, IMAT6U, AC2OLD, REPS  , MAXIT , &
			 IAMOUT, INOCNV, IDDLOW, IDDTOP, ISSTOP, &
			 IDCMIN, IDCMAX )
!  (This subroutine has not been used and tested yet) 
!
!****************************************************************
!
      USE SWCOMM1
      USE SWCOMM3
      USE SWCOMM4
      USE OCPCOMM4
      USE M_PARALL
      USE ALL_VARS, ONLY : MT

      IMPLICIT NONE
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: M. Zijlema                                    |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 2004-2005  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. Authors
!
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     40.41, Nov. 04: New subroutine
!
!  2. Purpose
!
!     Solves penta-diagonal system of equations in
!     spectral space with point SOR method
!
!  4. Argument variables
!
!     AC2         action density
!     AC2OLD      action density at previous iteration
!     IAMOUT      control parameter indicating the amount of
!                 output required
!                 0: no output
!                 1: only fatal errors will be printed
!                 2: gives output concerning the iteration process
!                 3: additional information about the iteration
!                    is printed
!     IDCMAX      maximum counter in directional space
!     IDCMIN      minimum counter in directional space
!     IDDLOW      minimum direction that is propagated within a sweep
!     IDDTOP      maximum direction that is propagated within a sweep
!     IMAT5L      coefficients of lower diagonal in sigma-space
!     IMAT6U      coefficients of upper diagonal in sigma-space
!     IMATDA      coefficients of main diagonal
!     IMATLA      coefficients of lower diagonal in theta-space
!     IMATUA      coefficients of upper diagonal in theta-space
!     IMATRA      right-hand side
!     INOCNV      integer indicating number of grid points in which
!                 solver does not converged
!     ISSTOP      maximum frequency counter in a sweep
!     MAXIT       the maximum number of iterations to be performed in
!                 the linear solver
!     REPS        relative accuracy of the final approximation
!
      INTEGER IAMOUT, INOCNV, IDDLOW, IDDTOP, ISSTOP, MAXIT
      INTEGER IDCMIN(MSC), IDCMAX(MSC)
      REAL    REPS
      REAL    AC2(MDC,MSC,0:MT),                          &
              IMATDA(MDC,MSC), IMATRA(MDC,MSC),           &
	      IMAT5L(MDC,MSC), IMAT6U(MDC,MSC),           &
	      IMATLA(MDC,MSC), IMATUA(MDC,MSC),           &
	      AC2OLD(MDC,MSC)
!
!  5. Parameter variables
!
!     OMEG  :     relaxation parameter
!
      REAL    OMEG
      PARAMETER (OMEG=0.8)
!
!  6. Local variables
!
!     AC2I  :     intermediate action density
!     ICONV :     indicator for convergence (1=yes, 0=no)
!     ID    :     loop counter in theta-space
!     IDDL  :     minimum counter in theta-space of modulo MDC
!     IDDT  :     maximum counter in theta-space of modulo MDC
!     IDDUM :     loop counter
!     IDINF :     index of point ID with largest error in solution
!     IDM   :     index of point ID-1
!     IDMAX :     local array of maximum counter in theta-space
!     IDMIN :     local array of minimum counter in theta-space
!     IDP   :     index of point ID+1
!     IENT  :     number of entries
!     INVMDA:     inverse of main diagonal
!     IS    :     loop counter in sigma-space
!     ISINF :     index of point IS with largest error in solution
!     ISM   :     index of point IS-1
!     ISP   :     index of point IS+1
!     IT    :     iteration count
!     LOPERI:     auxiliary vector meant for computation in
!                 periodic theta-space
!     RES   :     residual
!     RESM  :     inf-norm of residual vector
!     RESM0 :     inf-norm of initial residual vector
!     UPPERI:     auxiliary vector meant for computation in
!                 periodic theta-space
!
      INTEGER ICONV, ID, IDINF, IDDL, IDDT, IDDUM, IDM, IDP, IENT,  &
              IS, ISINF, ISM, ISP, IT
      INTEGER IDMIN(MSC), IDMAX(MSC)
      REAL    AC2I, RES, RESM, RESM0
      REAL    LOPERI(MSC), UPPERI(MSC), INVMDA(MDC,MSC)
!
!  7. Common blocks used
!
!
!  8. Subroutines used
!
!     MSGERR           Writes error message
!     STRACE           Tracing routine for debugging
!
!  9. Subroutines calling
!
!     SWOMPU  (in SWANCOM1)
!
! 12. Structure
!
!     The system of equations is solved using the SOR technique in
!     pointwise manner
!     Note that with omeg=1, the Gauss-Seidel method is recovered
!
!     Convergence is reached, if the difference between two consecutive
!     iteration levels measured w.r.t. the maximum norm is smaller than
!     given tolerance
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWSOR')

!     --- initialize arrays

      LOPERI = 0.
      UPPERI = 0.
      INVMDA = 0.

!     --- in case of periodicity in theta-space, store values
!         of matrix coefficients corresponding to left bottom and
!         right top

      DO IS = 1, ISSTOP
         IF ( IDCMIN(IS).EQ.1 .AND. IDCMAX(IS).EQ.MDC ) THEN
           UPPERI(IS) = IMATLA(  1,IS)
           LOPERI(IS) = IMATUA(MDC,IS)
         END IF
      END DO

!     --- when no bins fall within the sweep, i.e. SECTOR = 0,
!         reset the bounds of sector as 1..MDC (routine SOLPRE
!         has clear the rows in the matrix that do not belong
!         to the sweep)

      DO IS = 1, ISSTOP
         IF ( IDCMIN(IS).LE.IDCMAX(IS) ) THEN
            IDMIN(IS) = IDCMIN(IS)
            IDMAX(IS) = IDCMAX(IS)
         ELSE
            IDMIN(IS) = 1
            IDMAX(IS) = MDC
         END IF
      END DO

!     --- store inverse of main diagonal

      DO IS = 1, ISSTOP
         DO IDDUM = IDMIN(IS), IDMAX(IS)
            ID  = MOD ( IDDUM - 1 + MDC , MDC ) + 1
            IF ( IMATDA(ID,IS).NE.0. ) THEN
               INVMDA(ID,IS) = 1./IMATDA(ID,IS)
            ELSE
               CALL MSGERR ( 4,'Main diagonal of spectral matrix is zero!' )
               RETURN
            END IF
         END DO
      END DO

      IT    = 0
      ICONV = 0

!     --- start iteration process

 10   IF ( ICONV.EQ.0 .AND. IT.LT.MAXIT ) THEN

         IT    = IT + 1
         ICONV = 1
         RESM  = 0.
         IDINF = 0
         ISINF = 0

         IS  = 1
         ISP = IS + 1

         IDDUM = IDMIN(IS)
         ID    = MOD ( IDDUM - 1 + MDC , MDC ) + 1
         IDP   = MOD ( IDDUM     + MDC , MDC ) + 1
         IDDT  = MOD ( IDMAX(IS) - 1 + MDC , MDC ) + 1

         AC2I = IMATRA(ID,IS)                               &
	       -IMAT6U(ID,IS)*AC2(ID ,ISP,KCGRD(1))         &
	       -IMATUA(ID,IS)*AC2(IDP,IS ,KCGRD(1))         &
	       -UPPERI(IS)*AC2(IDDT,IS,KCGRD(1))
         AC2I = AC2I*OMEG*INVMDA(ID,IS)+(1.-OMEG)*AC2(ID,IS,KCGRD(1))

         RES = ABS(AC2(ID,IS,KCGRD(1)) - AC2I)
         IF ( RES.GT.RESM ) THEN
            RESM  = RES
            IDINF = ID
            ISINF = IS
         END IF
         AC2(ID,IS,KCGRD(1)) = AC2I

         DO IDDUM = IDMIN(IS)+1, IDMAX(IS)-1
            ID  = MOD ( IDDUM - 1 + MDC , MDC ) + 1
            IDM = MOD ( IDDUM - 2 + MDC , MDC ) + 1
            IDP = MOD ( IDDUM     + MDC , MDC ) + 1

            AC2I = IMATRA(ID,IS)                            &
	          -IMAT6U(ID,IS)*AC2(ID ,ISP,KCGRD(1))      &
		  -IMATLA(ID,IS)*AC2(IDM,IS ,KCGRD(1))      &
		  -IMATUA(ID,IS)*AC2(IDP,IS ,KCGRD(1))
            AC2I = AC2I*OMEG*INVMDA(ID,IS)+(1.-OMEG)*AC2(ID,IS,KCGRD(1))

            RES = ABS(AC2(ID,IS,KCGRD(1)) - AC2I)
            IF ( RES.GT.RESM ) THEN
               RESM  = RES
               IDINF = ID
               ISINF = IS
            END IF
            AC2(ID,IS,KCGRD(1)) = AC2I

         END DO

         IDDUM = IDMAX(IS)
         ID    = MOD ( IDDUM - 1 + MDC , MDC ) + 1
         IDM   = MOD ( IDDUM - 2 + MDC , MDC ) + 1
         IDDL  = MOD ( IDMIN(IS) - 1 + MDC , MDC ) + 1

         AC2I = IMATRA(ID,IS)                                &
	       -IMAT6U(ID,IS)*AC2(ID ,ISP,KCGRD(1))          &
	       -IMATLA(ID,IS)*AC2(IDM,IS ,KCGRD(1))          &
	       -LOPERI(IS)*AC2(IDDL,IS,KCGRD(1))
         AC2I = AC2I*OMEG*INVMDA(ID,IS)+(1.-OMEG)*AC2(ID,IS,KCGRD(1))

         RES = ABS(AC2(ID,IS,KCGRD(1)) - AC2I)
         IF ( RES.GT.RESM ) THEN
            RESM  = RES
            IDINF = ID
            ISINF = IS
         END IF
         AC2(ID,IS,KCGRD(1)) = AC2I

         DO IS = 2, ISSTOP-1
            ISM = IS - 1
            ISP = IS + 1

            IDDUM = IDMIN(IS)
            ID    = MOD ( IDDUM - 1 + MDC , MDC ) + 1
            IDP   = MOD ( IDDUM     + MDC , MDC ) + 1
            IDDT  = MOD ( IDMAX(IS) - 1 + MDC , MDC ) + 1

            AC2I = IMATRA(ID,IS)                          &
	          -IMAT5L(ID,IS)*AC2(ID ,ISM,KCGRD(1))    &
		  -IMAT6U(ID,IS)*AC2(ID ,ISP,KCGRD(1))    &
		  -IMATUA(ID,IS)*AC2(IDP,IS ,KCGRD(1))    &
		  -UPPERI(IS)*AC2(IDDT,IS,KCGRD(1))
            AC2I = AC2I*OMEG*INVMDA(ID,IS)+(1.-OMEG)*AC2(ID,IS,KCGRD(1))

            RES = ABS(AC2(ID,IS,KCGRD(1)) - AC2I)
            IF ( RES.GT.RESM ) THEN
               RESM  = RES
               IDINF = ID
               ISINF = IS
            END IF
            AC2(ID,IS,KCGRD(1)) = AC2I

            DO IDDUM = IDMIN(IS)+1, IDMAX(IS)-1
               ID  = MOD ( IDDUM - 1 + MDC , MDC ) + 1
               IDM = MOD ( IDDUM - 2 + MDC , MDC ) + 1
               IDP = MOD ( IDDUM     + MDC , MDC ) + 1

               AC2I = IMATRA(ID,IS)                        &
	             -IMAT5L(ID,IS)*AC2(ID ,ISM,KCGRD(1))  &
		     -IMAT6U(ID,IS)*AC2(ID ,ISP,KCGRD(1))  &
		     -IMATLA(ID,IS)*AC2(IDM,IS ,KCGRD(1))  &
		     -IMATUA(ID,IS)*AC2(IDP,IS ,KCGRD(1))
               AC2I = AC2I*OMEG*INVMDA(ID,IS)+             &
	              (1.-OMEG)*AC2(ID,IS,KCGRD(1))

               RES = ABS(AC2(ID,IS,KCGRD(1)) - AC2I)
               IF ( RES.GT.RESM ) THEN
                  RESM  = RES
                  IDINF = ID
                  ISINF = IS
               END IF
               AC2(ID,IS,KCGRD(1)) = AC2I

            END DO

            IDDUM = IDMAX(IS)
            ID    = MOD ( IDDUM - 1 + MDC , MDC ) + 1
            IDM   = MOD ( IDDUM - 2 + MDC , MDC ) + 1
            IDDL  = MOD ( IDMIN(IS) - 1 + MDC , MDC ) + 1

            AC2I = IMATRA(ID,IS)                            &
	          -IMAT5L(ID,IS)*AC2(ID ,ISM,KCGRD(1))      &
		  -IMAT6U(ID,IS)*AC2(ID ,ISP,KCGRD(1))      &
		  -IMATLA(ID,IS)*AC2(IDM,IS ,KCGRD(1))      &
		  -LOPERI(IS)*AC2(IDDL,IS,KCGRD(1))
            AC2I = AC2I*OMEG*INVMDA(ID,IS)+(1.-OMEG)*AC2(ID,IS,KCGRD(1))

            RES = ABS(AC2(ID,IS,KCGRD(1)) - AC2I)
            IF ( RES.GT.RESM ) THEN
               RESM  = RES
               IDINF = ID
               ISINF = IS
            END IF
            AC2(ID,IS,KCGRD(1)) = AC2I

         END DO

         IS  = ISSTOP
         ISM = IS - 1

         IDDUM = IDMIN(IS)
         ID    = MOD ( IDDUM - 1 + MDC , MDC ) + 1
         IDP   = MOD ( IDDUM     + MDC , MDC ) + 1
         IDDT  = MOD ( IDMAX(IS) - 1 + MDC , MDC ) + 1

         AC2I = IMATRA(ID,IS)                             &
	       -IMAT5L(ID,IS)*AC2(ID ,ISM,KCGRD(1))       &
	       -IMATUA(ID,IS)*AC2(IDP,IS ,KCGRD(1))       &
	       -UPPERI(IS)*AC2(IDDT,IS,KCGRD(1))
         AC2I = AC2I*OMEG*INVMDA(ID,IS)+(1.-OMEG)*AC2(ID,IS,KCGRD(1))

         RES = ABS(AC2(ID,IS,KCGRD(1)) - AC2I)
         IF ( RES.GT.RESM ) THEN
            RESM  = RES
            IDINF = ID
            ISINF = IS
         END IF
         AC2(ID,IS,KCGRD(1)) = AC2I

         DO IDDUM = IDMIN(IS)+1, IDMAX(IS)-1
            ID  = MOD ( IDDUM - 1 + MDC , MDC ) + 1
            IDM = MOD ( IDDUM - 2 + MDC , MDC ) + 1
            IDP = MOD ( IDDUM     + MDC , MDC ) + 1

            AC2I = IMATRA(ID,IS)                         &
	          -IMAT5L(ID,IS)*AC2(ID ,ISM,KCGRD(1))   &
		  -IMATLA(ID,IS)*AC2(IDM,IS ,KCGRD(1))   &
		  -IMATUA(ID,IS)*AC2(IDP,IS ,KCGRD(1))
            AC2I = AC2I*OMEG*INVMDA(ID,IS)+(1.-OMEG)*AC2(ID,IS,KCGRD(1))

            RES = ABS(AC2(ID,IS,KCGRD(1)) - AC2I)
            IF ( RES.GT.RESM ) THEN
               RESM  = RES
               IDINF = ID
               ISINF = IS
            END IF
            AC2(ID,IS,KCGRD(1)) = AC2I

         END DO

         IDDUM = IDMAX(IS)
         ID    = MOD ( IDDUM - 1 + MDC , MDC ) + 1
         IDM   = MOD ( IDDUM - 2 + MDC , MDC ) + 1
         IDDL  = MOD ( IDMIN(IS) - 1 + MDC , MDC ) + 1

         AC2I = IMATRA(ID,IS)                             &
	       -IMAT5L(ID,IS)*AC2(ID ,ISM,KCGRD(1))       &
	       -IMATLA(ID,IS)*AC2(IDM,IS ,KCGRD(1))       &
	       -LOPERI(IS)*AC2(IDDL,IS,KCGRD(1))
         AC2I = AC2I*OMEG*INVMDA(ID,IS)+(1.-OMEG)*AC2(ID,IS,KCGRD(1))

         RES = ABS(AC2(ID,IS,KCGRD(1)) - AC2I)
         IF ( RES.GT.RESM ) THEN
            RESM  = RES
            IDINF = ID
            ISINF = IS
         END IF
         AC2(ID,IS,KCGRD(1)) = AC2I

         IF ( RESM.GT.1.E8 ) THEN
            IT = MAXIT + 1
            ICONV = 0
            GOTO 10
         END IF
         IF ( IAMOUT.EQ.2 ) THEN
            WRITE (PRINTF,'(A,I3,A,E13.6,A,I3,A,I3,A)')     &
	          ' ++ SWSOR: iter = ',IT,'    res = ',RESM, &
		  ' in (ID,IS) = (',IDINF,',',ISINF,')'
         END IF
         IF ( IAMOUT.EQ.3 .AND. IT.EQ.1 ) RESM0 = RESM

         IF ( RESM.GT.REPS ) ICONV = 0
         GOTO 10

      END IF

!     --- investigate the reason to stop

      IF ( ICONV.EQ.0 ) THEN
         AC2(:,:,KCGRD(1)) = AC2OLD(:,:)
         INOCNV = INOCNV + 1
      END IF
      IF ( ICONV.EQ.0 .AND. IAMOUT.GE.1 ) THEN
!         IF (ERRPTS.GT.0.AND.INODE.EQ.MASTER) THEN
!            WRITE(ERRPTS,100) IXCGRD(1)+MXF-1, IYCGRD(1)+MYF-1, 2
!         END IF
 100     FORMAT (I4,1X,I4,1X,I2)
         WRITE (PRINTF,'(A,I5,A,I5,A)')                       &
	      ' ++ SWSOR: no convergence in grid point (',    &
	      IXCGRD(1)+MXF-1,',',IYCGRD(1)+MYF-1,')'
         WRITE (PRINTF,'(A,I3)')                              &
	      '           total number of iterations       = ',IT
         WRITE (PRINTF,'(A,E13.6)')                           &
	      '           inf-norm of the residual         = ',RESM
         WRITE (PRINTF,'(A,E13.6)')                           &
	      '           required accuracy                = ',REPS
      ELSE IF ( IAMOUT.EQ.3 ) THEN
         WRITE (PRINTF,'(A,E13.6)')                           &
	      ' ++ SWSOR: inf-norm of the initial residual = ',RESM0
         WRITE (PRINTF,'(A,I3)')                              &
	      '           total number of iterations       = ',IT
         WRITE (PRINTF,'(A,E13.6)')                           &
	      '           inf-norm of the residual         = ',RESM
      END IF

!     --- test output

      IF ( TESTFL .AND. ITEST.GE.120 ) THEN
         WRITE(PRTEST,*)
         WRITE(PRTEST,*) '  Subroutine SWSOR'
         WRITE(PRTEST,*)
         WRITE(PRTEST,200) KCGRD(1), MDC, MSC
 200     FORMAT(' SWSOR : POINT MDC MSC              :',3I5)
         WRITE(PRTEST,250) IDDLOW, IDDTOP, ISSTOP
 250     FORMAT(' SWSOR : IDDLOW IDDTOP ISSTOP       :',3I4)
         WRITE(PRTEST,*)
         WRITE(PRTEST,*) ' coefficients of matrix and rhs  '
         WRITE(PRTEST,*)
         WRITE(PRTEST,'(A111)')                                    &
	        ' IS ID         IMATLA         IMATDA'//           &
		'         IMATUA         IMATRA         IMAT5L'//  &
		'         IMAT6U            AC2'
         DO IDDUM = IDDLOW, IDDTOP
            ID = MOD ( IDDUM - 1 + MDC , MDC ) + 1
            DO IS = 1, ISSTOP
               WRITE(PRTEST,300) IS, ID,                           &
	                         IMATLA(ID,IS), IMATDA(ID,IS),     &
				 IMATUA(ID,IS), IMATRA(ID,IS),     &
				 IMAT5L(ID,IS), IMAT6U(ID,IS),     &
				 AC2(ID,IS,KCGRD(1))
 300           FORMAT(2I3,7E15.7)
            END DO
         END DO
         WRITE(PRTEST,*)
         WRITE(PRTEST,*)'IS ID      LPER          UPER '
         DO IDDUM = IDDLOW, IDDTOP
            ID = MOD ( IDDUM - 1 + MDC , MDC ) + 1
            IF ( ID.EQ.1 .OR. ID.EQ.MDC ) THEN
               DO IS = 1, ISSTOP
                  WRITE(PRTEST,350) IS, ID, LOPERI(IS), UPPERI(IS)
 350              FORMAT(2I3,2E15.7)
               END DO
            END IF
         END DO
      END IF

!     --- set matrix coefficients to zero

      IMATDA = 0.
      IMATRA = 0.
      IMATLA = 0.
      IMATUA = 0.
      IMAT5L = 0.
      IMAT6U = 0.

      RETURN
      END SUBROUTINE SWSOR
!****************************************************************
!
   SUBROUTINE SWMTLB ( N1, N2, M1, M2 )

!  (This subroutine has not been used and tested yet) 
!
!****************************************************************
!
   USE OCPCOMM4                                                        
!
   IMPLICIT NONE
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering                              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmers: R.C. Ris, N. Booij,                          |
!     |              IJ.G. Haagsma, A.T.M.M. Kieftenburg,         |
!     |              M. Zijlema, E.E. Kriezi,                     |
!     |              R. Padilla-Hernandez, L.H. Holthuijsen       |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 2004-2005  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. Authors
!
!     40.31: Tim Campbell and John Cazes
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     40.31, Jul. 03: New subroutine
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Given global loop bounds N1 and N2, compute loop bounds
!     M1 and M2 for calling thread
!
!  4. Argument variables
!
!     M1          lower index of thread loop
!     M2          upper index of thread loop
!     N1          lower index of global loop
!     N2          upper index of global loop
!
   INTEGER N1, N2, M1, M2
!
!  6. Local variables
!
!     ID    :     thread number
!     IENT  :     number of entries
!     NCH   :     auxiliary integer
!     NTH   :     number of threads
!
   INTEGER ID, IENT, NTH, NCH
!$     INTEGER OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM
!$     EXTERNAL OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM
!
!  8. Subroutines used
!
!     ---
!
!  9. Subroutines calling
!
!     ---
!
! 10. Error messages
!
!     ---
!
! 11. Remarks
!
!     ---
!
! 12. Structure
!
!     Description of the pseudo code
!
! 13. Source text
!
!   SAVE IENT
!   DATA IENT/0/
!   IF (LTRACE) CALL STRACE (IENT,'SWMTLB')

   NTH = 1
   ID  = 0
!$  NTH = OMP_GET_NUM_THREADS()
!$  ID  = OMP_GET_THREAD_NUM()
   NCH = (N2-N1+1)/NTH
   M1  = ID*NCH+N1
   M2  = (ID+1)*NCH+N1-1
   IF(ID == NTH-1) M2 = N2

   RETURN
   END SUBROUTINE SWMTLB
 
!****************************************************************
!
      SUBROUTINE SWSTPC ( HSACC0    ,HSACC1    ,HSACC2    ,     &
                          SACC1     ,SACC2     ,HSDIFC    ,     &
			  DELHS     ,DELTM     ,DEP2      ,     &
			  ACCUR     ,I1MYC     ,I2MYC     )

!     (This subroutine has not been used and tested yet)
!
!****************************************************************
!
      USE OCPCOMM4                                                        
      USE SWCOMM3                                                         
      USE SWCOMM4                                                         
      USE M_GENARR
      USE M_PARALL
!      USE ALL_VARS, ONLY : MT,AC2
      USE VARS_WAVE, ONLY : MT,AC2

      IMPLICIT NONE
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering                              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmers: R.C. Ris, N. Booij,                          |
!     |              IJ.G. Haagsma, A.T.M.M. Kieftenburg,         |
!     |              M. Zijlema, E.E. Kriezi,                     |
!     |              R. Padilla-Hernandez, L.H. Holthuijsen       |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 2004-2005  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. Authors
!
!     40.41: Andre van der Westhuysen
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     40.41, Jun. 04: New subroutine
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Check convergence based on the relative, absolute
!     and curvature values of significant wave height
!
!  4. Argument variables
!
!     ACCUR       indicates percentage of grid points in
!                 which accuracy is reached
!     DELHS       difference in Hs between last 2 iterations
!     DELTM       difference in Tm01 between last 2 iterations
!     DEP2        depth
!     HSACC0      significant wave height at iter-2
!     HSACC1      significant wave height at iter-1
!     HSACC2      significant wave height at iter
!     HSDIFC      difference of Hs(i) - Hs(i-2) meant for
!                 computation of curvature of Hs
!     I1MYC       lower index for thread loop over y-grid row
!     I2MYC       upper index for thread loop over y-grid row
!     SACC1       mean wave frequency at iter-1
!     SACC2       mean wave frequency at iter
!
      INTEGER I1MYC, I2MYC
      REAL    ACCUR
      REAL    DEP2(MT)          ,                  &
              HSACC0(MT)        ,                  &
	      HSACC1(MT)        ,                  &
	      HSACC2(MT)        ,                  &
	      SACC1(MT)         ,                  &
	      SACC2(MT)         ,                  &
	      DELHS(MT)         ,                  &
	      DELTM(MT)         ,                  &
	      HSDIFC(MT)
!
!  6. Local variables
!
!     ACS2  :     auxiliary variable
!     ACS3  :     auxiliary variable
!     HSABS :     absolute value of Hs
!     HSCURV:     curvature value of Hs
!     HSDIFO:     previous value of HSDIFC
!     HSREL :     relative value of Hs
!     IACCUR:     indicates number of grid points in which
!                 accuracy is reached
!     IARR  :     auxiliary array meant for global reduction
!     ID    :     counter of direction
!     IENT  :     number of entries
!     II    :     loop variable
!     INDX  :     index for indirect address
!     IS    :     counter of frequency
!     IX    :     loop counter
!     IX1   :     lower index in x-direction
!     IX2   :     upper index in x-direction
!     IY    :     loop counter
!     IY1   :     lower index in y-direction
!     IY2   :     upper index in y-direction
!     LHEAD :     logical indicating to write header
!     TMABS :     absolute value of Tm01
!     TSTFL :     indicates whether grid point is a test point
!     WETGRD:     number of wet grid points
!     XMOM0 :     zeroth moment
!     XMOM1 :     first moment
!
      INTEGER ID, IS, IENT, II, INDX, IX, IY, IX1, IX2, IY1, IY2
      INTEGER :: IG
      INTEGER IACCUR, WETGRD, IACCURt, WETGRDt, IARR(2)
      REAL    ACS2, ACS3, HSREL ,HSABS, HSCURV, HSDIFO, TMABS,      &
              XMOM0, XMOM1
      LOGICAL LHEAD, TSTFL
!
!  7. Common blocks used
!
      COMMON/SWSTPC_MT_COM/WETGRD,IACCUR
!
!     SWSTPC_MT_COM    place local summed variables WETGRD and IACCUR
!                      in common block so they will be scoped as shared
!
!  8. Subroutines used
!
!     EQREAL           Boolean function which compares two REAL values
!     STRACE           Tracing routine for debugging
!     STPNOW           Logical indicating whether program must
!                      terminated or not
!     SWREDUCE         Performs a global reduction
!
      LOGICAL EQREAL, STPNOW
!
!  9. Subroutines calling
!
!     SWCOMP (in SWANCOM1)
!
! 12. Structure
!
!     master thread initialize the shared variables
!     store Hs and Tm as old values and count number of wet grid points
!     compute new values of Hs and Tm
!     calculate a set of accuracy parameters based on relative,
!         absolute and curvature values of Hs and check accuracy
!     global sum of IACCUR and WETGRD
!     carry out reductions across all nodes
!
! 13. Source text

      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWSTPC')

!     --- master thread initialize the shared variables
!$OMP MASTER
      WETGRD = 0
      IACCUR = 0
!$OMP END MASTER
!$OMP BARRIER

      IF ( LMXF ) THEN
         IX1 = 1
      ELSE
         IX1 = 1+IHALOX
      END IF
      IF ( LMXL ) THEN
         IX2 = MXC
      ELSE
         IX2 = MXC-IHALOX
      END IF
      IF ( LMYF ) THEN
         IY1 = I1MYC
      ELSE
         IY1 = 1+IHALOY
      END IF
      IF ( LMYL ) THEN
         IY2 = I2MYC
      ELSE
         IY2 = MYC-IHALOY
      END IF

!     --- store Hs and Tm as old values and count number of wet grid points

      WETGRDt = 0
!      DO IX = IX1, IX2
         DO IG = 1, MT
            INDX = KGRPNT(IG)
            IF ( DEP2(INDX).GT.DEPMIN ) THEN
               HSACC0(INDX) = MAX( 1.E-20 , HSACC1(INDX) )
               HSACC1(INDX) = MAX( 1.E-20 , HSACC2(INDX) )
               SACC1 (INDX) = MAX( 1.E-20 , SACC2 (INDX) )
               WETGRDt = WETGRDt + 1
            ELSE
               HSACC0(INDX) = 0.
               HSACC1(INDX) = 0.
               SACC1 (INDX) = 0.
            END IF
         END DO
!      END DO

!     --- compute new values of Hs and Tm

!      DO IX = IX1, IX2
         DO IG = 1, MT
            INDX = KGRPNT(IG)

            IF ( DEP2(INDX).GT.DEPMIN ) THEN

               XMOM0 = 0.
               XMOM1 = 0.
               DO IS = 1, MSC
                  DO ID = 1, MDC
                     ACS2  = SPCSIG(IS)**2 * AC2(ID,IS,INDX)
                     ACS3  = SPCSIG(IS) * ACS2
                     XMOM0 = XMOM0 + ACS2
                     XMOM1 = XMOM1 + ACS3
                  END DO
               END DO
               XMOM0 = XMOM0 * FRINTF * DDIR
               XMOM1 = XMOM1 * FRINTF * DDIR

               IF ( XMOM0.GT.0. ) THEN
                  HSACC2(INDX) = MAX ( 1.E-20 , 4.*SQRT(XMOM0) )
                  SACC2 (INDX) = MAX ( 1.E-20 , (XMOM1/XMOM0) )
               ELSE
                  HSACC2(INDX) = 1.E-20
                  SACC2 (INDX) = 1.E-20
               END IF

            END IF

         END DO
!      END DO

      IACCURt = 0

!     --- calculate a set of accuracy parameters based on relative,
!         absolute and curvature values of Hs and check accuracy

      LHEAD=.TRUE.
!      DO IX = IX1, IX2
         DO IG = 1, MT
            INDX = KGRPNT(IG)

!           --- determine whether the point is a test point

            TSTFL = .FALSE.
            IF (NPTST.GT.0) THEN
               DO 20 II = 1, NPTST
                  IF (IX.NE.XYTST(2*II-1)) GOTO 20
                  IF (IY.NE.XYTST(2*II  )) GOTO 20
                  TSTFL = .TRUE.
  20           CONTINUE
            END IF

            DELHS(INDX) = 0.0
            DELTM(INDX) = 0.0
            IF ( DEP2(INDX).GT.DEPMIN ) THEN

               HSABS = ABS ( HSACC2(INDX) - HSACC1(INDX) )
               HSREL = HSABS / HSACC2(INDX)
               TMABS = ABS ( (PI2_W/SACC2(INDX)) - (PI2_W/SACC1(INDX)) )

               HSDIFO       = HSDIFC(INDX)
               HSDIFC(INDX) = 0.5*( HSACC2(INDX) - HSACC0(INDX) )
               HSCURV       = ABS(HSDIFC(INDX) - HSDIFO)/HSACC2(INDX)

               DELHS(INDX) = HSABS
               IF (EQREAL(SACC1(INDX),1.E-20) .OR.               &
	           EQREAL(SACC2(INDX),1.E-20) ) THEN
                  DELTM(INDX) = 0.
               ELSE
                  DELTM(INDX) = TMABS
               END IF

!              --- add gridpoint in which wave height has reached
!                  required accuracy

               IF ( HSCURV.LE.PNUMS(15) .AND.                    &
	           (HSREL.LE.PNUMS(1) .OR. HSABS.LE.PNUMS(2)) ) THEN
                  IACCURt = IACCURt + 1
               END IF

!               IF (TSTFL) THEN
!                  IF (LHEAD) WRITE(PRINTF,501)
!                  WRITE(PRINTF,502) IX+MXF-2, IY+MYF-2, HSABS, HSREL,  &
!		                    HSCURV
! 501              FORMAT(25X,'dHabs          ','dHrel          ',      &
!                        'Curvature ')
! 502              FORMAT(1X,SS,'(IX,IY)=(',I5,',',I5,')','  ',         &
!                    1PE13.6E2,'  ',1PE13.6E2,'  ',1PE13.6E2)
!                  LHEAD=.FALSE.
!               END IF

            END IF

         END DO
!      END DO
!
!     --- global sum of IACCUR and WETGRD
!$OMP ATOMIC
      IACCUR = IACCUR + IACCURt
!$OMP ATOMIC
      WETGRD = WETGRD + WETGRDt

!     --- carry out reductions across all nodes

!!$OMP BARRIER
!!$OMP MASTER
!      IARR(1) = IACCUR
!      IARR(2) = WETGRD
!      CALL SWREDUCE ( IARR, 2, SWINT, SWSUM )
!!OMP#ifndef _OPENMP
!      IF (STPNOW()) RETURN
!!OMP#endif
      IACCUR = IARR(1)
      WETGRD = IARR(2)
      ACCUR  = REAL(IACCUR) * 100. / REAL(WETGRD)
!$OMP END MASTER
!$OMP BARRIER
!
!     --- test output
!
!$OMP MASTER
      IF ( ITEST.GE.30 ) THEN
        WRITE(PRINTF,1002) PNUMS(1), PNUMS(2), PNUMS(15)
 1002   FORMAT(' SWSTPC: DHREL DHABS CURV      :',3E12.4)
        WRITE(PRINTF,1008) WETGRD,IACCUR,ACCUR
 1008   FORMAT(' SWSTPC: WETGRD IACCUR ACCUR   :',2I8,E12.4)
      END IF
!$OMP END MASTER

      RETURN
      END SUBROUTINE SWSTPC

# endif


!****************************************************************
!
  SUBROUTINE ICETT2
!
!****************************************************************
!
  USE M_GENARR
  USE SWCOMM2
  USE SWCOMM3
  USE SWCOMM4
  USE OCPCOMM4
  USE M_PARALL
  USE ALL_VARS
  USE VARS_WAVE

  INTEGER :: I,ISS,ID

  DO I = 1,MT
   DO ISS = 1,MSC
    DO ID = 1,MDC
     if (Sice(ID,ISS,I)>0) then
       AC2(ID,ISS,I)=AC2(ID,ISS,I)-Sice(ID,ISS,I)
       AC2(ID,ISS,I)=MAX(0.0,AC2(ID,ISS,I))
     endif
    END DO
   END DO
  END DO

  RETURN
  END SUBROUTINE ICETT2

