!
!******************************************************************
!
   SUBROUTINE FAC4WW (XIS   ,SNLC1 ,DAL1  ,DAL2  ,DAL3  ,SPCSIG,  &
		      WWINT ,WWAWG ,WWSWG )        
!(This subroutine has not been tested yet)
!
!******************************************************************

!
   USE SWCOMM3                                                         
   USE SWCOMM4                                                         
   USE OCPCOMM4                                                        
   USE M_SNL4                                                          
!
!  Purpose :
!
!     Calculate interpolation constants for Snl.
!
   REAL    SPCSIG(MSC)                                                 
   INTEGER     MSC2  ,MSC1  ,IS    ,IDP   ,IDP1  ,                &
               IDM   ,IDM1  ,ISP   ,ISP1  ,ISM   ,ISM1  ,         &
	       IDPP  ,IDMM  ,ISPP  ,ISMM  ,                       &
	       ISLOW ,ISHGH ,ISCLW ,ISCHG ,IDLOW ,IDHGH ,         &
               MSCMAX,MDCMAX                     

   REAL        SNLC1 ,LAMM2 ,LAMP2 ,DELTH3,                       &
               AUX1  ,DELTH4,DAL1  ,DAL2  ,DAL3  ,CIDP  ,WIDP  ,  &
	       WIDP1 ,CIDM  ,WIDM  ,WIDM1 ,XIS   ,XISLN ,WISP  ,  &
	       WISP1 ,WISM  ,WISM1 ,AWG1  ,AWG2  ,AWG3  ,AWG4  ,  &
	       AWG5  ,AWG6  ,AWG7  ,AWG8  ,SWG1  ,SWG2  ,SWG3  ,  &
	       SWG4  ,SWG5  ,SWG6  ,SWG7  ,SWG8  ,FREQ  ,         &
	       RADE                                            
!
   REAL       WWAWG(*),WWSWG(*)
!
   INTEGER    WWINT(*)
!

   IF(ALLOCATED(AF11)) DEALLOCATE(AF11)                               

!  *** Compute frequency indices                               ***
!  *** XIS is the relative increment of the relative frequency ***
!
   MSC2   = INT ( FLOAT(MSC) / 2.0 )
   MSC1   = MSC2 - 1

   XIS    = SPCSIG(MSC2) / SPCSIG(MSC1)  

!
!  *** set values for the nonlinear four-wave interactions ***
!
   SNLC1  = 1. / GRAV_W**4                                               
!
   LAMM2  = (1.-PQUAD(1))**2                                           
   LAMP2  = (1.+PQUAD(1))**2                                           
   DELTH3 = ACOS( (LAMM2**2+4.-LAMP2**2) / (4.*LAMM2) )
   AUX1   = SIN(DELTH3)
   DELTH4 = ASIN(-AUX1*LAMM2/LAMP2)
!
   DAL1   = 1. / (1.+PQUAD(1))**4                                      
   DAL2   = 1. / (1.-PQUAD(1))**4                                      
   DAL3   = 2. * DAL1 * DAL2
!
!  *** Compute directional indices in sigma and theta space ***
!
   CIDP   = ABS(DELTH4/DDIR)                                           
   IDP   = INT(CIDP)
   IDP1  = IDP + 1
   WIDP   = CIDP - REAL(IDP)
   WIDP1  = 1.- WIDP
!
   CIDM   = ABS(DELTH3/DDIR)                                           
   IDM   = INT(CIDM)
   IDM1  = IDM + 1
   WIDM   = CIDM - REAL(IDM)
   WIDM1  = 1.- WIDM
   XISLN  = LOG( XIS )
!
   ISP    = INT( LOG(1.+PQUAD(1)) / XISLN )                            
   ISP1   = ISP + 1
   WISP   = (1.+PQUAD(1) - XIS**ISP) / (XIS**ISP1 - XIS**ISP)          
   WISP1  = 1. - WISP
!
   ISM    = INT( LOG(1.-PQUAD(1)) / XISLN )                            
   ISM1   = ISM - 1
   WISM   = (XIS**ISM -(1.-PQUAD(1))) / (XIS**ISM - XIS**ISM1)         
   WISM1  = 1. - WISM
!
!  *** Range of calculations ***
!
   ISLOW =  1  + ISM1
   ISHGH = MSC + ISP1 - ISM1
   ISCLW =  1
   ISCHG = MSC - ISM1
   IDLOW = 1 - MDC - MAX(IDM1,IDP1)
   IDHGH = MDC + MDC + MAX(IDM1,IDP1)
!
   MSC4MI = ISLOW
   MSC4MA = ISHGH
   MDC4MI = IDLOW
   MDC4MA = IDHGH
   MSCMAX = MSC4MA - MSC4MI + 1
   MDCMAX = MDC4MA - MDC4MI + 1
!
!  *** Interpolation weights ***
!
   AWG1   = WIDP  * WISP
   AWG2   = WIDP1 * WISP
   AWG3   = WIDP  * WISP1
   AWG4   = WIDP1 * WISP1
!
   AWG5   = WIDM  * WISM
   AWG6   = WIDM1 * WISM
   AWG7   = WIDM  * WISM1
   AWG8   = WIDM1 * WISM1
!
!  *** quadratic interpolation ***
!
   SWG1   = AWG1**2
   SWG2   = AWG2**2
   SWG3   = AWG3**2
   SWG4   = AWG4**2
!
   SWG5   = AWG5**2
   SWG6   = AWG6**2
   SWG7   = AWG7**2
   SWG8   = AWG8**2
!
!  --- determine discrete counters for piecewise                       
!      constant interpolation                                          
!
   IF(AWG1 < AWG2)THEN
     IF(AWG2 < AWG3)THEN
       IF(AWG3 < AWG4)THEN
         ISPP=ISP
         IDPP=IDP
       ELSE
         ISPP=ISP
         IDPP=IDP1
       END IF
     ELSE IF(AWG2 < AWG4)THEN
       ISPP=ISP
       IDPP=IDP
     ELSE
       ISPP=ISP1
       IDPP=IDP
     END IF
   ELSE IF(AWG1 < AWG3)THEN
     IF(AWG3 < AWG4)THEN
       ISPP=ISP
       IDPP=IDP
     ELSE
       ISPP=ISP
       IDPP=IDP1
     END IF
   ELSE IF(AWG1 < AWG4)THEN
     ISPP=ISP
     IDPP=IDP
   ELSE
     ISPP=ISP1
     IDPP=IDP1
   END IF
   IF(AWG5 < AWG6)THEN
     IF(AWG6 < AWG7)THEN
       IF(AWG7 < AWG8)THEN
         ISMM=ISM
         IDMM=IDM
       ELSE
         ISMM=ISM
         IDMM=IDM1
       END IF
     ELSE IF(AWG6 < AWG8)THEN
       ISMM=ISM
       IDMM=IDM
     ELSE
       ISMM=ISM1
       IDMM=IDM
     END IF
   ELSE IF(AWG5 < AWG7)THEN
     IF(AWG7 < AWG8)THEN
       ISMM=ISM
       IDMM=IDM
     ELSE
       ISMM=ISM
       IDMM=IDM1
     END IF
   ELSE IF(AWG5 < AWG8)THEN
     ISMM=ISM
     IDMM=IDM
   ELSE
     ISMM=ISM1
     IDMM=IDM1
   END IF
!
!  *** fill the arrays *
!
   WWINT(1) = IDP
   WWINT(2) = IDP1
   WWINT(3) = IDM
   WWINT(4) = IDM1
   WWINT(5) = ISP
   WWINT(6) = ISP1
   WWINT(7) = ISM
   WWINT(8) = ISM1
   WWINT(9) = ISLOW
   WWINT(10)= ISHGH
   WWINT(11)= ISCLW
   WWINT(12)= ISCHG
   WWINT(13)= IDLOW
   WWINT(14)= IDHGH
   WWINT(15)= MSC4MI
   WWINT(16)= MSC4MA
   WWINT(17)= MDC4MI
   WWINT(18)= MDC4MA
   WWINT(19)= MSCMAX
   WWINT(20)= MDCMAX
   WWINT(21)= IDPP                                                     
   WWINT(22)= IDMM                                                     
   WWINT(23)= ISPP                                                     
   WWINT(24)= ISMM                                                     
!
   WWAWG(1) = AWG1
   WWAWG(2) = AWG2
   WWAWG(3) = AWG3
   WWAWG(4) = AWG4
   WWAWG(5) = AWG5
   WWAWG(6) = AWG6
   WWAWG(7) = AWG7
   WWAWG(8) = AWG8
!
   WWSWG(1) = SWG1
   WWSWG(2) = SWG2
   WWSWG(3) = SWG3
   WWSWG(4) = SWG4
   WWSWG(5) = SWG5
   WWSWG(6) = SWG6
   WWSWG(7) = SWG7
   WWSWG(8) = SWG8

   ALLOCATE (AF11(MSC4MI:MSC4MA))                                      

!  *** Fill scaling array (f**11)                     ***
!  *** compute the radian frequency**11 for IS=1, MSC ***
!
   DO IS=1, MSC
     AF11(IS) = ( SPCSIG(IS) / ( 2. * PI_W ) )**11                       
   END DO  
!
!  *** compute the radian frequency for the IS = MSC+1, ISHGH ***
!
   FREQ   = SPCSIG(MSC) / ( 2. * PI_W )                                  
   DO IS = MSC+1, ISHGH
     FREQ   = FREQ * XIS
     AF11(IS) = FREQ**11
   END DO  
!
!  *** compute the radian frequency for IS = 0, ISLOW ***
!
   FREQ   = SPCSIG(1) / ( 2. * PI_W )                                    
   DO IS = 0, ISLOW, -1
     FREQ   = FREQ / XIS
     AF11(IS) = FREQ**11
   END DO  

   RETURN
   END SUBROUTINE FAC4WW
 
!
!****************************************************************
!
   SUBROUTINE SWLTA ( IG    ,DEP2   ,CGO    ,SPCSIG ,          &
                      KWAVE  ,IMATRA ,IMATDA ,IDDLOW ,          &
		      IDDTOP ,ISSTOP ,IDCMIN ,IDCMAX ,          &
		      HS     ,SMEBRK ,PLTRI  ,URSELL )
! (This subroutine has not been tested yet)		      
!
!****************************************************************
!
   USE OCPCOMM4
   USE SWCOMM3
   USE SWCOMM4
   USE VARS_WAVE, ONLY : MT, AC2
!
   IMPLICIT NONE
!
!  2. Purpose
!
!     In this subroutine the triad-wave interactions are calculated
!     with the Lumped Triad Approximation of Eldeberky (1996). His
!     expression is based on a parametrization of the biphase (as
!     function of the Ursell number), is directionally uncoupled and
!     takes into account for the self-self interactions only.
!
!     For a full description of the equations reference is made
!     to PhD thesis of Eldeberky (1996). Here only the main expressions
!     are given.
!
!  3. Method
!
!     The parametrized biphase is given by (see eq. 3.19):
!
!                                  0.2
!     beta = - pi/2 + pi/2 tanh ( ----- )
!                                   Ur
!
!     The Ursell number is calculated in routine SINTGRL
!
!     The source term as function of frequency p is (see eq. 7.25):
!
!             +      -
!     S(p) = S(p) + S(p)
!
!     in which
!
!      +
!     S(p) = alpha Cp Cg,p (R(p/2,p/2))**2 sin (|beta|) ( E(p/2)**2 -2 E(p) E(p/2) )
!
!      -          +
!     S(p) = - 2 S(2p)
!
!     with alpha a tunable coefficient and R(p/2,p/2) is the interaction
!     coefficient of which the expression can be found in Eldeberky (1996);
!     see eq. 7.26.
!
!     Note that a slightly adapted formulation of the LTA is used in
!     in the SWAN model:
!
!     - Only positive contributions to higher harmonics are considered
!       here (no energy is transferred to lower harmonics).
!
!     - The mean frequency in the expression of the Ursell number
!       is calculated according to the first order moment over the
!       zeroth order moment (personal communication, Y.Eldeberky, 1997).
!
!     - The interactions are calculated up to 2.5 times the mean
!       frequency only.
!
!     - Since the spectral grid is logarithmically distributed in frequency
!       space, the interactions between central bin and interacting bin
!       are interpolated such that the distance between these bins is
!       factor 2 (nearly).
!
!     - The interactions are calculated in terms of energy density
!       instead of action density. So the action density spectrum
!       is firstly converted to the energy density grid, then the
!       interactions are calculated and then the spectrum is converted
!       to the action density spectrum back.
!
!     - To ensure numerical stability the Patankar rule is used.
!
   INTEGER IDDLOW, IDDTOP, ISSTOP
   INTEGER IDCMIN(MSC), IDCMAX(MSC)

   REAL :: HS, SMEBRK
!  REAL :: AC2(MDC,MSC,0:MT)

   REAL :: CGO(MSC,MICMAX)
   REAL :: DEP2(MT)
   REAL :: IMATDA(MDC,MSC), IMATRA(MDC,MSC)
   REAL :: SPCSIG(MSC)
   REAL :: KWAVE(MSC,MICMAX)
   REAL :: PLTRI(MDC,MSC,NPTST)
   REAL :: URSELL(MT)

   INTEGER I1, I2, ID, IG, IDDUM, IENT, II, IS, ISM, ISM1, ISMAX, ISP, ISP1
   REAL    AUX1, AUX2, BIPH, C0, CM, DEP, DEP_2, DEP_3, E0, EM,     &
           FT, RINT, SIGPI, SINBPH, STRI, WISM, WISM1, WISP, WISP1, &
	   W0, WM, WN0, WNM, XIS, XISLN
   REAL, ALLOCATABLE :: E(:), SA(:,:)

!   DEP   = DEP2(IGC)
   DEP   = DEP2(IG)
   DEP_2 = DEP**2
   DEP_3 = DEP**3
!
!  --- compute some indices in sigma space
!
   I2     = INT (FLOAT(MSC) / 2.)
   I1     = I2 - 1
   XIS    = SPCSIG(I2) / SPCSIG(I1)
   XISLN  = LOG( XIS )

   ISP    = INT( LOG(2.) / XISLN )
   ISP1   = ISP + 1
   WISP   = (2. - XIS**ISP) / (XIS**ISP1 - XIS**ISP)
   WISP1  = 1. - WISP

   ISM    = INT( LOG(0.5) / XISLN )
   ISM1   = ISM - 1
   WISM   = (XIS**ISM -0.5) / (XIS**ISM - XIS**ISM1)
   WISM1  = 1. - WISM

   ALLOCATE (E (1:MSC))
   ALLOCATE (SA(1:MDC,1:MSC+ISP1))
   E  = 0.
   SA = 0.
!
!  --- compute maximum frequency for which interactions are calculated
!
   ISMAX = 1
   DO IS = 1, MSC
     IF(SPCSIG(IS) < ( PTRIAD(2) * SMEBRK) )THEN
       ISMAX = IS
     ENDIF
   ENDDO
   ISMAX = MAX ( ISMAX , ISP1 )
!
!  --- compute 3 wave-wave interactions
!
!   IF( URSELL(IGC) >= PTRIAD(5) )THEN
   IF( URSELL(IG) >= PTRIAD(5) )THEN
!
!    --- calculate biphase
!
!     BIPH   = (0.5*PI_W)*(TANH(PTRIAD(4)/URSELL(IGC))-1.)
     BIPH   = (0.5*PI_W)*(TANH(PTRIAD(4)/URSELL(IG))-1.)
     SINBPH = ABS( SIN(BIPH) )
!
     DO II = IDDLOW, IDDTOP
       ID = MOD ( II - 1 + MDC , MDC ) + 1
!
!      --- initialize array with E(f) for the direction considered
!
       DO IS = 1, MSC
!         E(IS) = AC2(ID,IS,IGC) * 2. * PI_W * SPCSIG(IS)
         E(IS) = AC2(ID,IS,IG) * 2. * PI_W * SPCSIG(IS)
       END DO
!
       DO IS = 1, ISMAX

         E0  = E(IS)
         W0  = SPCSIG(IS)
         WN0 = KWAVE(IS,1)
         C0  = W0 / WN0

         IF( IS > -ISM1 )THEN
           EM  = WISM * E(IS+ISM1)       + WISM1 * E(IS+ISM)
           WM  = WISM * SPCSIG(IS+ISM1)  + WISM1 * SPCSIG(IS+ISM)
           WNM = WISM * KWAVE(IS+ISM1,1) + WISM1 * KWAVE(IS+ISM,1)
           CM  = WM / WNM
         ELSE
           EM  = 0.
           WM  = 0.
           WNM = 0.
           CM  = 0.
         END IF

         AUX1 = WNM**2 * ( GRAV_W * DEP + 2.*CM**2 )
         AUX2 = WN0 * DEP * ( GRAV_W * DEP +                            &
	        (2./15.) * GRAV_W * DEP_3 * WN0**2 -                    &
		(2./ 5.) * W0**2 * DEP_2 )
         RINT = AUX1 / AUX2
         FT = PTRIAD(1) * C0 * CGO(IS,1) * RINT**2 * SINBPH

         SA(ID,IS) = MAX(0., FT * ( EM * EM - 2. * EM * E0 ))

       END DO
     END DO
!
!    ---  put source term together
!
     DO IS = 1, ISSTOP
       SIGPI = SPCSIG(IS) * 2. * PI_W
       DO IDDUM = IDCMIN(IS), IDCMAX(IS)
         ID = MOD ( IDDUM - 1 + MDC , MDC ) + 1
!
         STRI = SA(ID,IS) - 2.*(WISP  * SA(ID,IS+ISP1) +     &
	        WISP1 * SA(ID,IS+ISP ))
!
!        --- store results in rhs and main diagonal according
!            to Patankar-rules
!
         IF(TESTFL) PLTRI(ID,IS,IPTST) = STRI / SIGPI
         IF(STRI > 0.)THEN
           IMATRA(ID,IS) = IMATRA(ID,IS) + STRI / SIGPI
         ELSE
           IMATDA(ID,IS) = IMATDA(ID,IS) - STRI /     &
!	                   MAX(1.E-18,AC2(ID,IS,IGC)*SIGPI)
	                   MAX(1.E-18,AC2(ID,IS,IG)*SIGPI)
         END IF
       END DO
     END DO

   END IF

   DEALLOCATE(E,SA)

   RETURN
   END SUBROUTINE SWLTA

!
!******************************************************************
!
   SUBROUTINE RANGE4 (WWINT,IDDLOW,IDDTOP)                            
!
!******************************************************************
!
!     calculate the minimum and maximum counters in frequency and
!     directional space which fall with the calculation for the
!     nonlinear wave-wave interactions.
!
!     Method :  review for the counters :
!
!                            Frequencies -->
!                 +---+---------------------+---------+- IDHGH
!              d  | 3 :          2          :    2    |
!              i  + - + - - - - - - - - - - + - - - - +- MDC
!              r  |   :                     :         |
!              e  | 3 :  original spectrum  :    1    |
!              c  |   :                     :         |
!              t. + - + - - - - - - - - - - + - - - - +- 1
!                 | 3 :          2          :    2    |
!                 +---+---------------------+---------+- IDLOW
!                 |   |                     |    ^    |
!             ISLOW   1                     MSC  |  ISHGH
!                     ^                          |
!                     |                          |
!                    ISCLW                     ISCHG
!              lowest discrete               highest discrete
!                central bin                   central bin
!
!
!       The directional counters depend on the numerical method that
!       is used.
!****************************************************************
!
   USE SWCOMM3                                                         
   USE SWCOMM4                                                         
   USE OCPCOMM4 
   IMPLICIT NONE

   INTEGER     IDDLOW,IDDTOP                                           
!
   INTEGER     WWINT(*)
!
!   *** Range in directional domain ***
!
   IF(IQUAD < 3 .AND. IQUAD > 0)THEN                         
!   *** counters based on bins which fall within a sweep ***
     WWINT(13) = IDDLOW - MAX( WWINT(4), WWINT(2) )
     WWINT(14) = IDDTOP + MAX( WWINT(4), WWINT(2) )
   ELSE
!   *** counters initially based on full circle ***
     WWINT(13) = 1   - MAX( WWINT(4), WWINT(2) )
     WWINT(14) = MDC + MAX( WWINT(4), WWINT(2) )
   END IF
!
!  *** error message ***
!
!   IF(WWINT(9) < WWINT(15) .OR. WWINT(10) > WWINT(16) .OR.          &
!      WWINT(13) < WWINT(17) .OR. WWINT(14) > WWINT(18))THEN
!     WRITE (PRINTF,900) IXCGRD(1), IYCGRD(1),                       &
!                        WWINT(9) ,WWINT(10) ,WWINT(13) ,WWINT(14),  &
!			WWINT(15),WWINT(16) ,WWINT(17) ,WWINT(18)
!900  FORMAT ( ' ** Error : array bounds and maxima in subr RANGE4, ', &
!              ' point ', 2I5,                                         &
!	      /,'            ISL,ISH : ',2I4, '   IDL,IDH : ',2I4,    &
!	      /,'            SMI,SMA : ',2I4, '   DMI,DMA : ',2I4)
!     IF(ITEST > 50) WRITE (PRTEST, 901) MSC, MDC, IDDLOW, IDDTOP
!901  FORMAT (' MSC, MDC, IDDLOW, IDDTOP: ', 4I5)
!   ENDIF
!
!  test output
!
!   IF(TESTFL .AND. ITEST >= 60)THEN
!     WRITE(PRTEST,911) WWINT(4), WWINT(2), WWINT(8), WWINT(6)
!911  FORMAT (' RANGE4: IDM1 IDP1 ISM1 ISP1    :',4I5)
!     WRITE(PRTEST,916) WWINT(11), WWINT(12), IQUAD
!916  FORMAT (' RANGE4: ISCLW ISCHG IQUAD      :',3I5)
!     WRITE (PRTEST,917) WWINT(9), WWINT(10), WWINT(13), WWINT(14)
!917  FORMAT (' RANGE4: ISLOW ISHGH IDLOW IDHGH:',4I5)
!     WRITE (PRTEST,919) WWINT(15), WWINT(16), WWINT(17), WWINT(18)
!919  FORMAT (' RANGE4: MS4MI MS4MA MD4MI MD4MA:',4I5)
!     WRITE(PRINTF,*)
!   END IF
!
   RETURN
   END SUBROUTINE RANGE4

!
!*******************************************************************
!
   SUBROUTINE FILNL3 (IDCMIN  ,IDCMAX  ,IMATRA  ,IMATDA  ,      &
                      MEMNL4  ,PLNL4S  ,ISSTOP  ,IGC            )    
!(This subroutine has not been tested yet)
!
!*******************************************************************
!
   USE SWCOMM3                                                         
   USE SWCOMM4                                                         
   USE OCPCOMM4 
!   USE ALL_VARS, ONLY : MT,AC2                                                       
   USE VARS_WAVE, ONLY : MT,AC2                                                       
   
   IMPLICIT NONE
!
!  2. Purpose
!
!     Fill the IMATRA/IMATDA arrays with the nonlinear wave-wave interaction
!     source term for a gridpoint ix,iy per sweep direction
!
!*******************************************************************
!
   INTEGER   IS,ID,ISSTOP,IDDUM,IGC
!
   REAL      IMATRA(MDC,MSC)           ,                          &
             IMATDA(MDC,MSC)           ,                          &
!	     AC2(MDC,MSC,0:MT)        ,                          &
	     PLNL4S(MDC,MSC,NPTST)     ,                          &
	     MEMNL4(MDC,MSC,MT)                                
!
   INTEGER   IDCMIN(MSC),IDCMAX(MSC)
!
!   SAVE IENT
!   DATA IENT/0/
!   IF (LTRACE) CALL STRACE (IENT,'FILNL3')
!
   DO IS=1, ISSTOP
     DO IDDUM = IDCMIN(IS), IDCMAX(IS)
       ID = MOD ( IDDUM - 1 + MDC , MDC ) + 1
       IF(TESTFL) PLNL4S(ID,IS,IPTST) = MEMNL4(ID,IS,IGC)         
       IF(MEMNL4(ID,IS,IGC) > 0.)THEN                          
          IMATRA(ID,IS) = IMATRA(ID,IS) + MEMNL4(ID,IS,IGC)
       ELSE
          IMATDA(ID,IS) = IMATDA(ID,IS) - MEMNL4(ID,IS,IGC) /     &
	                  MAX(1.E-18,AC2(ID,IS,IGC))
       END IF
     END DO
   END DO    
!
!   IF(TESTFL .AND. ITEST >= 50)THEN
!     WRITE(PRINTF,9000) IDCMIN(1),IDCMAX(1),MSC,ISSTOP
!9000 FORMAT(' FILNL3: ID_MIN ID_MAX MSC ISTOP :',4I6)
!     IF(ITEST >= 100)THEN
!       DO IS=1, ISSTOP
!         DO IDDUM = IDCMIN(IS), IDCMAX(IS)
!           ID = MOD ( IDDUM - 1 + MDC , MDC ) + 1
!           WRITE(PRINTF,6001) IS,ID,MEMNL4(ID,IS,KCGRD(1))             
!6001       FORMAT(' FILNL3: IS ID MEMNL()          :',2I6,E12.4)
!         ENDDO
!       ENDDO
!     ENDIF
!   ENDIF
!
   RETURN
   END SUBROUTINE FILNL3

!
!*********************************************************************
   SUBROUTINE SWSNL8 (WWINT   ,UE      ,SA1     ,SA2     ,SPCSIG  ,     &
                      SNLC1   ,DAL1    ,DAL2    ,DAL3    ,SFNL    ,     &
		      DEP2    ,KMESPC  ,MEMNL4  ,FACHFR  )
		      
!  (This subroutine has not been tested yet)		      
!*********************************************************************
!
   USE SWCOMM3                                                         
   USE SWCOMM4                                                         
   USE OCPCOMM4                                                        
   USE M_SNL4
!   USE ALL_VARS, ONLY : MT,AC2
   USE VARS_WAVE, ONLY : MT,AC2
!
   IMPLICIT NONE
!
!  2. Purpose
!
!     Calculate non-linear interaction using the discrete interaction
!     approximation (Hasselmann and Hasselmann 1985; WAMDI group 1988)
!     for the full circle.
!
!  3. Method
!
!     Discrete interaction approximation. To make interpolation simple,
!     the interactions are calculated in a "folded" space.
!
!                            Frequencies -->
!                 +---+---------------------+---------+- IDHGH
!              d  | 3 :          2          :    2    |
!              i  + - + - - - - - - - - - - + - - - - +- MDC
!              r  |   :                     :         |
!              e  | 3 :  original spectrum  :    1    |
!              c  |   :                     :         |
!              t. + - + - - - - - - - - - - + - - - - +- 1
!                 | 3 :          2          :    2    |
!                 +---+---------------------+---------+- IDLOW
!                 |   |                     |     ^   |
!              ISLOW  1                    MSC    |   ISHGH
!                     |                           |
!                   ISCLW                        ISCHG
!              lowest discrete               highest discrete
!                central bin                   central bin
!
!                            1 : Extra tail added beyond MSC
!                            2 : Spectrum copied outside ID range
!                            3 : Empty bins at low frequencies
!
!     ISLOW =  1  + ISM1
!     ISHGH = MSC + ISP1 - ISM1
!     ISCLW =  1
!     ISCHG = MSC - ISM1
!     IDLOW =  1  - MAX(IDM1,IDP1)
!     IDHGH = MDC + MAX(IDM1,IDP1)
!
!     Note: using this subroutine requires an additional array
!           with size MXC*MYC*MDC*MSC.
!
   INTEGER WWINT(*)
!
   REAL    DAL1, DAL2, DAL3, FACHFR, KMESPC, SNLC1
!   REAL    AC2(MDC,MSC,0:MT)
   REAL    DEP2(MT)
   REAL    MEMNL4(MDC,MSC,MT)
   REAL    SA1(MSC4MI:MSC4MA,MDC4MI:MDC4MA)
   REAL    SA2(MSC4MI:MSC4MA,MDC4MI:MDC4MA)
   REAL    SFNL(MSC4MI:MSC4MA,MDC4MI:MDC4MA)
   REAL    SPCSIG(MSC)
   REAL    UE(MSC4MI:MSC4MA,MDC4MI:MDC4MA)
!
!*******************************************************************
!
   INTEGER   IS      ,ID      ,ID0     ,I       ,J       ,         &
             ISLOW   ,ISHGH   ,IDLOW   ,IDHGH   ,                  &
	     IDPP    ,IDMM    ,ISPP    ,ISMM    ,                  &
	     ISCLW   ,ISCHG   ,IDDUM   ,IENT
!
   REAL      X       ,X2      ,CONS    ,FACTOR  ,SNLCS1  ,SNLCS2  ,  &
             SNLCS3  ,E00     ,EP1     ,EM1     ,EP2     ,EM2     ,  &
	     SA1A    ,SA1B    ,SA2A    ,SA2B    ,                    &
	     JACOBI  ,SIGPI
!
!   SAVE IENT
!   DATA IENT/0/
!   IF (LTRACE) CALL STRACE (IENT,'SWSNL8')
!
   ISLOW  = WWINT(9)
   ISHGH  = WWINT(10)
   ISCLW  = WWINT(11)
   ISCHG  = WWINT(12)
   IDLOW  = WWINT(13)
   IDHGH  = WWINT(14)
   IDPP   = WWINT(21)
   IDMM   = WWINT(22)
   ISPP   = WWINT(23)
   ISMM   = WWINT(24)
!
!  *** Calculate prop. constant.                           ***
!  *** Calculate factor R(X) to calculate the NL wave-wave ***
!  *** interaction for shallow water                       ***
!  *** SNLC1 = 1/GRAV**4                                   ***
!
   SNLCS1 = PQUAD(3)
   SNLCS2 = PQUAD(4)
   SNLCS3 = PQUAD(5)
!   X      = MAX ( 0.75 * DEP2(IGC) * KMESPC , 0.5 )
   X      = MAX ( 0.75 * DEP2(kcgrd(1)) * KMESPC , 0.5 )
   X2     = MAX ( -1.E15, SNLCS3*X)
   CONS   = SNLC1 * ( 1. + SNLCS1/X * (1.-SNLCS2*X) * EXP(X2))
   JACOBI = 2. * PI_W
!
!  *** extend the area with action density at periodic boundaries ***
!
   DO IDDUM = IDLOW, IDHGH
     ID = MOD ( IDDUM - 1 + MDC , MDC ) + 1
     DO IS=1, MSC
!       UE (IS,IDDUM) = AC2(ID,IS,IGC) * SPCSIG(IS) * JACOBI
       UE (IS,IDDUM) = AC2(ID,IS,kcgrd(1)) * SPCSIG(IS) * JACOBI
     ENDDO
   ENDDO
!
   DO ID = IDLOW, IDHGH
     DO IS = MSC+1, ISHGH
       UE(IS,ID) = UE(IS-1,ID) * FACHFR
     ENDDO
   ENDDO
!
!  *** Calculate (unfolded) interactions ***
!  *** Energy at interacting bins        ***
!
   DO ID = 1, MDC
     DO IS = ISCLW, ISCHG
       E00 = UE(IS     ,ID     )
       EP1 = UE(IS+ISPP,ID+IDPP)
       EM1 = UE(IS+ISMM,ID-IDMM)
       EP2 = UE(IS+ISPP,ID-IDPP)
       EM2 = UE(IS+ISMM,ID+IDMM)
!
!      Contribution to interactions
!
       FACTOR = CONS * AF11(IS) * PQUAD(2) * E00
!
       SA1A   = E00 * ( EP1*DAL1 + EM1*DAL2 )
       SA1B   = SA1A - EP1*EM1*DAL3
       SA2A   = E00 * ( EP2*DAL1 + EM2*DAL2 )
       SA2B   = SA2A - EP2*EM2*DAL3
!
       SA1 (IS,ID) = FACTOR * SA1B
       SA2 (IS,ID) = FACTOR * SA2B
!
     ENDDO
   ENDDO
!
!  *** Fold interactions to side angles -> domain in theta is ***
!  *** periodic                                               ***
!
   DO ID = 1, IDHGH - MDC
     ID0   = 1 - ID
     DO IS = ISCLW, ISCHG
       SA1 (IS,MDC+ID) = SA1 (IS,  ID   )
       SA2 (IS,MDC+ID) = SA2 (IS,  ID   )
       SA1 (IS,  ID0 ) = SA1 (IS,MDC+ID0)
       SA2 (IS,  ID0 ) = SA2 (IS,MDC+ID0)
     ENDDO
   ENDDO
!
!  *** Put source term together (To save space I=IS and ***
!  *** J=MDC is used)  ----                             ***
!
   DO I = 1, MSC
     SIGPI = SPCSIG(I) * JACOBI
     DO J = 1, MDC
       SFNL(I,J) =   - 2. * ( SA1(I,J) + SA2(I,J) )                 &
                   + ( SA1(I-ISPP,J-IDPP) + SA2(I-ISPP,J+IDPP) )    &
		   + ( SA1(I-ISMM,J+IDMM) + SA2(I-ISMM,J-IDMM) )
!
!      *** store value in auxiliary array and use values in ***
!      *** next four sweeps (see subroutine FILNL3)         ***
!
!       MEMNL4(J,I,IGC) = SFNL(I,J) / SIGPI
       MEMNL4(J,I,kcgrd(1)) = SFNL(I,J) / SIGPI
     ENDDO
   ENDDO
!
!  *** value source term in every bin ***
!
!   IF(ITEST >= 150 .AND. TESTFL)THEN
!     DO I=1, MSC
!       DO J=1, MDC
!         WRITE(PRINTF,2006) I,J,MEMNL4(J,I,KCGRD(1)),SFNL(I,J),SPCSIG(I)
!2006     FORMAT (' I J MEMNL() SFNL() SPCSIG:',2I4,3E12.4)
!       ENDDO
!     ENDDO
!   END IF
!
   RETURN
!
   END SUBROUTINE SWSNL8

!
!********************************************************************
!
!     <<  Numerical Computations of the Nonlinear Energy Transfer
!           of Gravity Wave Spectra in Finite Water Depth  >>
!
!                  << developed by Noriaki Hashimoto >>
!
!        References:  N. Hashimoto et al. (1998)
!                     Komatsu and Masuda (2000) (in Japanese)
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
      SUBROUTINE RIAM_SLW(LMAX     ,N        ,N2       ,G        ,  &
                          H        ,DQ       ,DQ2      ,DT       ,  &
			  DT2      ,W        ,P        ,ACT      ,  &
			  SNL      ,MINT     )   

!     (This subroutine has not been tested yet)

      USE M_SNL4                                                          

!     LMAX
!     N     : number of directional bins                                  
!     N2
!     G     : gravitational acceleration                                  
!     H     : depth                                                       
!     DQ
!     DQ2
!     DT    : size of the directional bins (Delta Theta)                  
!     DT2
!     W     : discretised frequency array                                 
!     P     : density of water                                            
!     ACT   : action density                                              
!     SNL   : Quadruplet source term                                      
!     MINT

!     IW4   : counter for the 4th frequency                               
!     W4    : frequency of the fourth quadruplet component                
!     AK4   : wavenumber of the fourth quadruplet component               
!     DNA   : coefficient in eq. 17 of Hashimoto (98)                     
!             = 2 w4 k4 / Cg(k4)                                          
!     IW3L
!     II
!     JJ
!     DI
!     DJ
!     CGK4  : group velocity for the fourth quadruplet component          

      REAL :: W(LMAX), ACT(LMAX,N), SNL(LMAX,N)                           

!     Initialisation of the quadruplet source term                        

      SNL(:,:) = 0.                                                       
!
!     =================
      DO 130 IW4=1,LMAX
!     =================
!
      W4=W(IW4)

!     WAVE converts nondimensional Sqr(w)d/g to nondimensional kd

      AK4=WAVE(W4**2*H/G)/H

!     CGCMP computes group velocity

      CALL CGCMP(G,AK4,H,CGK4)

!     Calculates the coefficient in equation (17) of Hashimoto (1998)

      DNA=2.*W4*AK4/CGK4

      IW3L=MAX0(1,NINT(IW4-ALOG(3.)/DQ))

!     ------------------------------------------------------------
      CALL PRESET(LMAX,N,IW3L,IW4,N2,G,H,W4,AK4,W,DQ,DQ2,DT,DT2,P)        
!     ------------------------------------------------------------
!
!     ==============
      DO 130 IT4=1,N
!     ==============
!
!     ================
      CURRIAM => FRIAM                                                    
      DO                !    (W1 - W3 - W4 - W2)                          
!     ================
!
      K1=CURRIAM%II(1)+IW4                                                
      K2=CURRIAM%II(2)+IW4                                                
      K3=CURRIAM%II(3)+IW4                                                
!
      IF(K1.LT.1   ) GO TO 120                                            
      IF(K2.GT.LMAX) GO TO 120                                            
!
      GG=CURRIAM%SSS*DNA                                                  
!
      M1P= CURRIAM%JJ(1)+IT4                                              
      M1N=-CURRIAM%JJ(1)+IT4                                              
      M2P= CURRIAM%JJ(2)+IT4                                              
      M2N=-CURRIAM%JJ(2)+IT4                                              
      M3P= CURRIAM%JJ(3)+IT4                                              
      M3N=-CURRIAM%JJ(3)+IT4                                              
!
      M1P=MOD(M1P,N)
      M1N=MOD(M1N,N)
      M2P=MOD(M2P,N)
      M2N=MOD(M2N,N)
      M3P=MOD(M3P,N)
      M3N=MOD(M3N,N)
!
      IF(M1P.LT.1) M1P=M1P+N
      IF(M1N.LT.1) M1N=M1N+N
      IF(M2P.LT.1) M2P=M2P+N
      IF(M2N.LT.1) M2N=M2N+N
      IF(M3P.LT.1) M3P=M3P+N
      IF(M3N.LT.1) M3N=M3N+N
!
      A4=ACT(IW4,IT4)
!
      IF(MINT.EQ.1) THEN
!
        K01=0
        K02=0
        K03=0
        M01=0
        M02=0
        M03=0
        IF(CURRIAM%DI(1).LE.1.) K01= 1                                    
        IF(CURRIAM%DI(2).LE.1.) K02= 1                                    
        IF(CURRIAM%DI(3).LE.1.) K03= 1                                    
        IF(CURRIAM%DJ(1).LE.1.) M01=-1                                    
        IF(CURRIAM%DJ(2).LE.1.) M02=-1                                    
        IF(CURRIAM%DJ(3).LE.1.) M03=-1                                    
!
        K11=K1+K01
        K21=K2+K02
        K31=K3+K03
        IF(K11.GT.LMAX) K11=LMAX
        IF(K21.GT.LMAX) K21=LMAX
        IF(K31.GT.LMAX) K31=LMAX
!
        K111=K1+K01-1
        K211=K2+K02-1
        K311=K3+K03-1
        IF(K111.LT.1) K111=1
        IF(K211.LT.1) K211=1
        IF(K311.LT.1) K311=1
!
        M1P1=M1P+M01
        M2P1=M2P+M02
        M3P1=M3P+M03
        IF(M1P1.LT.1) M1P1=N
        IF(M2P1.LT.1) M2P1=N
        IF(M3P1.LT.1) M3P1=N
!
        M1N1=M1N-M01
        M2N1=M2N-M02
        M3N1=M3N-M03
        IF(M1N1.GT.N) M1N1=1
        IF(M2N1.GT.N) M2N1=1
        IF(M3N1.GT.N) M3N1=1
!
        M1P11=M1P+M01+1
        M2P11=M2P+M02+1
        M3P11=M3P+M03+1
        IF(M1P11.GT.N) M1P11=1
        IF(M2P11.GT.N) M2P11=1
        IF(M3P11.GT.N) M3P11=1
!
        M1N11=M1N-M01-1
        M2N11=M2N-M02-1
        M3N11=M3N-M03-1
        IF(M1N11.LT.1) M1N11=N
        IF(M2N11.LT.1) M2N11=N
        IF(M3N11.LT.1) M3N11=N
!
        DI1=CURRIAM%DI(1)                                                 
        DI2=CURRIAM%DI(2)                                                 
        DI3=CURRIAM%DI(3)                                                 
        DJ1=CURRIAM%DJ(1)                                                 
        DJ2=CURRIAM%DJ(2)                                                 
        DJ3=CURRIAM%DJ(3)                                                 
!
        A1P=(DI1*(DJ1*ACT(K11, M1P1)+ACT(K11, M1P11))                 &
	        + DJ1*ACT(K111,M1P1)+ACT(K111,M1P11))                 &
		/((1.+DI1)*(1.+DJ1))                                      
!
        A1N=(DI1*(DJ1*ACT(K11, M1N1)+ACT(K11, M1N11))                 &
	        + DJ1*ACT(K111,M1N1)+ACT(K111,M1N11))                 &
		/((1.+DI1)*(1.+DJ1))                                  
!
        A2P=(DI2*(DJ2*ACT(K21, M2P1)+ACT(K21, M2P11))                 &
	        + DJ2*ACT(K211,M2P1)+ACT(K211,M2P11))                 &
		/((1.+DI2)*(1.+DJ2))                                  
!
        A2N=(DI2*(DJ2*ACT(K21, M2N1)+ACT(K21, M2N11))                 &
	        + DJ2*ACT(K211,M2N1)+ACT(K211,M2N11))                 &
		/((1.+DI2)*(1.+DJ2))                                
!
        A3P=(DI3*(DJ3*ACT(K31, M3P1)+ACT(K31, M3P11))                 &
	        + DJ3*ACT(K311,M3P1)+ACT(K311,M3P11))                 &
		/((1.+DI3)*(1.+DJ3))                                 
!
        A3N=(DI3*(DJ3*ACT(K31, M3N1)+ACT(K31, M3N11))                 &
	        + DJ3*ACT(K311,M3N1)+ACT(K311,M3N11))                 &
		/((1.+DI3)*(1.+DJ3))                                 
!
      ELSE
!
        IF(K1.LT.1) K1=1
        IF(K2.LT.1) K2=1
        IF(K3.LT.1) K3=1
!
        A1P=ACT(K1,M1P)
        A1N=ACT(K1,M1N)
        A2P=ACT(K2,M2P)
        A2N=ACT(K2,M2N)
        A3P=ACT(K3,M3P)
        A3N=ACT(K3,M3N)
!
      ENDIF
!
      W1P2P=A1P+A2P
      W1N2N=A1N+A2N
      S1P2P=A1P*A2P
      S1N2N=A1N*A2N
      W3P4 =A3P+A4
      W3N4 =A3N+A4
      S3P4 =A3P*A4
      S3N4 =A3N*A4
!
      XP=S1P2P*W3P4-S3P4*W1P2P
      XN=S1N2N*W3N4-S3N4*W1N2N
!
      SNL( K1,M1P)=SNL( K1,M1P)-XP*GG
      SNL( K1,M1N)=SNL( K1,M1N)-XN*GG
      SNL( K2,M2P)=SNL( K2,M2P)-XP*GG
      SNL( K2,M2N)=SNL( K2,M2N)-XN*GG
      SNL( K3,M3P)=SNL( K3,M3P)+XP*GG
      SNL( K3,M3N)=SNL( K3,M3N)+XN*GG
      SNL(IW4,IT4)=SNL(IW4,IT4)+(XP+XN)*GG
!
! ============
  120 CONTINUE
      IF (.NOT.ASSOCIATED(CURRIAM%NEXTRIAM)) EXIT                         
      CURRIAM => CURRIAM%NEXTRIAM                                         
      END DO                                                              
! ============
!
! ============
  130 CONTINUE
! ============
!
      RETURN
      END SUBROUTINE RIAM_SLW

!
!*******************************************************************
!
   SUBROUTINE SWSNL4 (WWINT   ,WWAWG   ,SPCSIG  ,SNLC1   ,         &
		      DAL1    ,DAL2    ,DAL3    ,DEP2    ,         &
		      KMESPC  ,MEMNL4  ,FACHFR  ,IDIA    ,         &
		      ITER    )

!  (This subroutine has not been tested yet)		      
!
!*******************************************************************
!
   USE SWCOMM3                                                         
   USE SWCOMM4                                                         
   USE OCPCOMM4                                                        
   USE M_SNL4
!   USE ALL_VARS, ONLY : MT,AC2
   USE VARS_WAVE, ONLY : MT,AC2
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering                              |
!     | Fluid Mechanics Section                                   |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmers: H.L. Tolman, R.C. Ris                        |
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
!     40.17: IJsbrand Haagsma
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     40.17, Dec. 01: New Subroutine based on SWSNL3
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Calculate non-linear interaction using the discrete interaction
!     approximation (Hasselmann and Hasselmann 1985; WAMDI group 1988)
!     for the full circle (option if a current is present). Note: using
!     this subroutine requires an additional array with size
!     (MXC*MYC*MDC*MSC). This requires more internal memory but can
!     speed up the computations sigificantly if a current is present.
!
!  3. Method
!
!     Discrete interaction approximation. To make interpolation simple,
!     the interactions are calculated in a "folded" space.
!
!                            Frequencies -->
!                 +---+---------------------+---------+- IDHGH
!              d  | 3 :          2          :    2    |
!              i  + - + - - - - - - - - - - + - - - - +- MDC
!              r  |   :                     :         |
!              e  | 3 :  original spectrum  :    1    |
!              c  |   :                     :         |
!              t. + - + - - - - - - - - - - + - - - - +- 1
!                 | 3 :          2          :    2    |
!                 +---+---------------------+---------+- IDLOW
!                 |   |                     |     ^   |
!              ISLOW  1                    MSC    |   ISHGH
!                     |                           |
!                   ISCLW                        ISCHG
!              lowest discrete               highest discrete
!                central bin                   central bin
!
!                            1 : Extra tail added beyond MSC
!                            2 : Spectrum copied outside ID range
!                            3 : Empty bins at low frequencies
!
!     ISLOW =  1  + ISM1
!     ISHGH = MSC + ISP1 - ISM1
!     ISCLW =  1
!     ISCHG = MSC - ISM1
!     IDLOW =  1  - MAX(IDM1,IDP1)
!     IDHGH = MDC + MAX(IDM1,IDP1)
!
!       Relative offsets of interpolation points around central bin
!       "#" and corresponding numbers of AWGn :
!
!               ISM1  ISM
!                5        7    T |
!          IDM1   +------+     H +
!                 |      |     E |      ISP      ISP1
!                 |   \  |     T |       3           1
!           IDM   +------+     A +        +---------+  IDP1
!                6       \8      |        |         |
!                                |        |  /      |
!                           \    +        +---------+  IDP
!                                |      /4           2
!                              \ |  /
!          -+-----+------+-------#--------+---------+----------+
!                                |           FREQ.
!
!
!  4. Argument variables
!
!     MCGRD : number of wet grid points of the computational grid
!     MDC   : grid points in theta-direction of computational grid
!     MDC4MA: highest array counter in directional space (Snl4)
!     MDC4MI: lowest array counter in directional space (Snl4)
!     MSC   : grid points in sigma-direction of computational grid
!     MSC4MA: highest array counter in frequency space (Snl4)
!     MSC4MI: lowest array counter in frequency space (Snl4)
!     WWINT : counters for quadruplet interactions
!
   INTEGER WWINT(*)
   INTEGER IDIA
!
!     AC2   : action density
!     AF11  : scaling frequency
!     DAL1  : coefficient for the quadruplet interactions
!     DAL2  : coefficient for the quadruplet interactions
!     DAL3  : coefficient for the quadruplet interactions
!     DEP2  : depth
!     FACHFR
!     KMESPC: mean average wavenumber over full spectrum
!     MEMNL4
!     PI    : circular constant
!     SA1   : interaction contribution of first quadruplet (unfolded space)
!     SA2   : interaction contribution of second quadruplet (unfolded space)
!     SFNL
!     SNLC1
!     SPCSIG: relative frequencies in computational domain in sigma-space
!     UE    : "unfolded" spectrum
!     WWAWG : weight coefficients for the quadruplet interactions
!
   REAL    DAL1, DAL2, DAL3, FACHFR, KMESPC, SNLC1
!   REAL    AC2(MDC,MSC,0:MT)
   REAL    DEP2(MT)
   REAL    MEMNL4(MDC,MSC,MT)
   REAL    SPCSIG(MSC)
   REAL    WWAWG(*)
!
!  6. Local variables
!
   REAL, ALLOCATABLE :: SA1(:,:), SA2(:,:), SFNL(:,:), UE(:,:)
!
!  7. Common blocks used
!
!
!  8. Subroutines used
!
!     ---
!
!  9. Subroutines calling
!
!     SOURCE (in SWANCOM1)
!
! 12. Structure
!
!     -------------------------------------------
!       Initialisations.
!       Calculate proportionality constant.
!       Prepare auxiliary spectrum.
!       Calculate (unfolded) interactions :
!       -----------------------------------------
!         Energy at interacting bins
!         Contribution to interactions
!         Fold interactions to side angles
!       -----------------------------------------
!       Put source term together
!     -------------------------------------------
!
! 13. Source text
!
!*******************************************************************
!
   INTEGER   IS      ,ID      ,ID0     ,I       ,J       ,          &
             ISHGH   ,IDLOW   ,IDHGH   ,ISP     ,ISP1    ,          &
	     IDP     ,IDP1    ,ISM     ,ISM1    ,IDM     ,IDM1    , &
	     ISCLW   ,ISCHG
!
   REAL      X       ,X2      ,CONS    ,FACTOR  ,SNLCS2  ,          &
             SNLCS3  ,E00     ,EP1     ,EM1     ,EP2     ,EM2     , &
	     SA1A    ,SA1B    ,SA2A    ,SA2B    ,                   &
	     AWG1    ,AWG2    ,AWG3    ,AWG4    ,AWG5    ,AWG6    , &
	     AWG7    ,AWG8    ,                                     &
	     JACOBI  ,SIGPI
!
!   SAVE IENT
!   DATA IENT/0/
!   IF (LTRACE) CALL STRACE (IENT,'SWSNL4')

   ALLOCATE(SA1(MSC4MI:MSC4MA,MDC4MI:MDC4MA))
   ALLOCATE(SA2(MSC4MI:MSC4MA,MDC4MI:MDC4MA))
   ALLOCATE(SFNL(MSC4MI:MSC4MA,MDC4MI:MDC4MA))
   ALLOCATE(UE(MSC4MI:MSC4MA,MDC4MI:MDC4MA))
!
   IDP    = WWINT(1)
   IDP1   = WWINT(2)
   IDM    = WWINT(3)
   IDM1   = WWINT(4)
   ISP    = WWINT(5)
   ISP1   = WWINT(6)
   ISM    = WWINT(7)
   ISM1   = WWINT(8)
   ISLOW  = WWINT(9)
   ISHGH  = WWINT(10)
   ISCLW  = WWINT(11)
   ISCHG  = WWINT(12)
   IDLOW  = WWINT(13)
   IDHGH  = WWINT(14)
!
   AWG1 = WWAWG(1)
   AWG2 = WWAWG(2)
   AWG3 = WWAWG(3)
   AWG4 = WWAWG(4)
   AWG5 = WWAWG(5)
   AWG6 = WWAWG(6)
   AWG7 = WWAWG(7)
   AWG8 = WWAWG(8)
!
!  *** Initialize auxiliary arrays per gridpoint ***
!
   DO ID = MDC4MI, MDC4MA
     DO IS = MSC4MI, MSC4MA
       UE(IS,ID)   = 0.
       SA1(IS,ID)  = 0.
       SA2(IS,ID)  = 0.
       SFNL(IS,ID) = 0.
     ENDDO
   ENDDO
!
!  *** Calculate prop. constant.                           ***
!  *** Calculate factor R(X) to calculate the NL wave-wave ***
!  *** interaction for shallow water                       ***
!  *** SNLC1 = 1/GRAV**4                                   ***
!
   SNLCS1 = PQUAD(3)                                                   
   SNLCS2 = PQUAD(4)                                                   
   SNLCS3 = PQUAD(5)                                                   
!   X      = MAX ( 0.75 * DEP2(IGC) * KMESPC , 0.5 )
   X      = MAX ( 0.75 * DEP2(kcgrd(1)) * KMESPC , 0.5 )
   X2     = MAX ( -1.E15, SNLCS3*X)
   CONS   = SNLC1 * ( 1. + SNLCS1/X * (1.-SNLCS2*X) * EXP(X2))
   JACOBI = 2. * PI_W
!
!  *** extend the area with action density at periodic boundaries ***
!
   DO IDDUM = IDLOW, IDHGH
     ID = MOD ( IDDUM - 1 + MDC , MDC ) + 1
     DO IS=1, MSC
!       UE (IS,IDDUM) = AC2(ID,IS,IGC) * SPCSIG(IS) * JACOBI
       UE (IS,IDDUM) = AC2(ID,IS,kcgrd(1)) * SPCSIG(IS) * JACOBI
     ENDDO
   ENDDO
!
   DO IS = MSC+1, ISHGH
     DO ID = IDLOW, IDHGH
       UE(IS,ID) = UE(IS-1,ID) * FACHFR
     ENDDO
   ENDDO
!
!  *** Calculate (unfolded) interactions ***
!  *** Energy at interacting bins        ***
!
   DO IS = ISCLW, ISCHG
     DO ID = 1, MDC
       E00    =        UE(IS      ,ID      )
       EP1    = AWG1 * UE(IS+ISP1,ID+IDP1) +                           &
                AWG2 * UE(IS+ISP1,ID+IDP ) +                           &
		AWG3 * UE(IS+ISP ,ID+IDP1) +                           &
		AWG4 * UE(IS+ISP ,ID+IDP )
       EM1    = AWG5 * UE(IS+ISM1,ID-IDM1) +                           &
                AWG6 * UE(IS+ISM1,ID-IDM ) +                           &
		AWG7 * UE(IS+ISM ,ID-IDM1) +                           &
		AWG8 * UE(IS+ISM ,ID-IDM )
       EP2    = AWG1 * UE(IS+ISP1,ID-IDP1) +                           &
                AWG2 * UE(IS+ISP1,ID-IDP ) +                           &
		AWG3 * UE(IS+ISP ,ID-IDP1) +                           &
		AWG4 * UE(IS+ISP ,ID-IDP )
       EM2    = AWG5 * UE(IS+ISM1,ID+IDM1) +                           &
                AWG6 * UE(IS+ISM1,ID+IDM ) +                           &
		AWG7 * UE(IS+ISM ,ID+IDM1) +                           &
		AWG8 * UE(IS+ISM ,ID+IDM )
!
!      Contribution to interactions
!
       FACTOR = CONS * AF11(IS) * E00
!
       SA1A   = E00 * ( EP1*DAL1 + EM1*DAL2 ) * CNL4_1(IDIA)
       SA1B   = SA1A - EP1*EM1*DAL3 * CNL4_2(IDIA)
       SA2A   = E00 * ( EP2*DAL1 + EM2*DAL2 ) * CNL4_1(IDIA)
       SA2B   = SA2A - EP2*EM2*DAL3 * CNL4_2(IDIA)
!

       SA1 (IS,ID) = FACTOR * SA1B
       SA2 (IS,ID) = FACTOR * SA2B
!
!       IF(ITEST >= 100 .AND. TESTFL)THEN
!         WRITE(PRINTF,9002) E00,EP1,EM1,EP2,EM2
!9002     FORMAT (' E00 EP1 EM1 EP2 EM2  :',5E11.4)
!         WRITE(PRINTF,9003) SA1A,SA1B,SA2A,SA2B
!9003     FORMAT (' SA1A SA1B SA2A SA2B  :',4E11.4)
!         WRITE(PRINTF,9004) IS,ID,SA1(IS,ID),SA2(IS,ID)
!9004     FORMAT (' IS ID SA1() SA2()    :',2I4,2E12.4)
!         WRITE(PRINTF,9005) FACTOR,JACOBI
!9005     FORMAT (' FACTOR JACOBI        : ',2E12.4)
!       END IF
!
     ENDDO
   ENDDO
!
!  *** Fold interactions to side angles -> domain in theta is ***
!  *** periodic                                               ***
!
   DO ID = 1, IDHGH - MDC
     ID0   = 1 - ID
     DO IS = ISCLW, ISCHG
       SA1 (IS,MDC+ID) = SA1 (IS,  ID   )
       SA2 (IS,MDC+ID) = SA2 (IS,  ID   )
       SA1 (IS,  ID0 ) = SA1 (IS,MDC+ID0)
       SA2 (IS,  ID0 ) = SA2 (IS,MDC+ID0)
     ENDDO
   ENDDO
!
!  *** Put source term together (To save space I=IS and ***
!  *** J=MDC is used)                                   ***
!
   FAC = 1.                                                            

   DO I = 1, MSC
     SIGPI = SPCSIG(I) * JACOBI                                        
     DO J = 1, MDC
       SFNL(I,J) = - 2. * ( SA1(I,J) + SA2(I,J) )                        &
                   + AWG1 * ( SA1(I-ISP1,J-IDP1) + SA2(I-ISP1,J+IDP1) )  &
		   + AWG2 * ( SA1(I-ISP1,J-IDP ) + SA2(I-ISP1,J+IDP ) )  &
		   + AWG3 * ( SA1(I-ISP ,J-IDP1) + SA2(I-ISP ,J+IDP1) )  &
		   + AWG4 * ( SA1(I-ISP ,J-IDP ) + SA2(I-ISP ,J+IDP ) )  &
		   + AWG5 * ( SA1(I-ISM1,J+IDM1) + SA2(I-ISM1,J-IDM1) )  &
		   + AWG6 * ( SA1(I-ISM1,J+IDM ) + SA2(I-ISM1,J-IDM ) )  &
		   + AWG7 * ( SA1(I-ISM ,J+IDM1) + SA2(I-ISM ,J-IDM1) )  &
		   + AWG8 * ( SA1(I-ISM ,J+IDM ) + SA2(I-ISM ,J-IDM ) )
!
!      *** store value in auxiliary array and use values in ***
!      *** next four sweeps (see subroutine FILNL3)         ***
!
       IF(IDIA == 1)THEN
!         MEMNL4(J,I,IGC) = FAC * SFNL(I,J) / SIGPI
         MEMNL4(J,I,kcgrd(1)) = FAC * SFNL(I,J) / SIGPI
       ELSE
!         MEMNL4(J,I,IGC) = MEMNL4(J,I,IGC) + FAC * SFNL(I,J) / SIGPI
         MEMNL4(J,I,kcgrd(1)) = MEMNL4(J,I,kcgrd(1)) + FAC * SFNL(I,J) / SIGPI
       END IF
     ENDDO
   ENDDO
!
!  *** test output ***
!
!   IF(ITEST >= 50 .AND. TESTFL)THEN
!     WRITE(PRINTF,*)
!     WRITE(PRINTF,*) ' SWSNL4 subroutine '
!     WRITE(PRINTF,9011) IDP, IDP1, IDM, IDM1
!9011 FORMAT (' IDP IDP1 IDM IDM1     :',4I5)
!     WRITE (PRINTF,9013) ISP, ISP1, ISM, ISM1
!9013 FORMAT (' ISP ISP1 ISM ISM1     :',4I5)
!     WRITE (PRINTF,9015) ISLOW, ISHGH, IDLOW, IDHGH
!9015 FORMAT (' ISLOW ISHG IDLOW IDHG :',4I5)
!     WRITE(PRINTF,9016) ISCLW, ISCHG, JACOBI
!9016 FORMAT (' ICLW ICHG JACOBI      :',2I5,E12.4)
!     WRITE (PRINTF,9017) AWG1, AWG2, AWG3, AWG4
!9017 FORMAT (' AWG1 AWG2 AWG3 AWG4   :',4E12.4)
!     WRITE (PRINTF,9018) AWG5, AWG6, AWG7, AWG8
!9018 FORMAT (' AWG5 AWG6 AWG7 AWG8   :',4E12.4)
!     WRITE (PRINTF,9019) MSC4MI, MSC4MA, MDC4MI, MDC4MA
!9019 FORMAT (' S4MI S4MA D4MI D4MA   :',4I6)
!     WRITE(PRINTF,9020) SNLC1,X,X2,CONS
!9020 FORMAT (' SNLC1  X  X2  CONS    :',4E12.4)
!     WRITE(PRINTF,9021) DEP2(KCGRD(1)),KMESPC,FACHFR,PI
!9021 FORMAT (' DEPTH KMESPC FACHFR PI:',4E12.4)
!     WRITE(PRINTF,*)
!
!    *** value source term in every bin ***
!
!     IF(ITEST >= 150)THEN
!       DO I=1, MSC
!         DO J=1, MDC
!           WRITE(PRINTF,2006) I,J,MEMNL4(J,I,KCGRD(1)),SFNL(I,J),     &
!	                      SPCSIG(I)
!2006       FORMAT (' I J MEMNL() SFNL() SPCSIG:',2I4,3E12.4)
!         ENDDO
!       ENDDO
!     END IF
!   END IF
!
   DEALLOCATE (SA1, SA2, SFNL, UE)

   RETURN
!
   END SUBROUTINE SWSNL4

!
!************************************************************
!

   SUBROUTINE SWSNL3 (WWINT   ,WWAWG   ,UE      ,SA1     ,SA2     ,    &
                      SPCSIG  ,SNLC1   ,DAL1    ,DAL2    ,DAL3    ,    &
		      SFNL    ,DEP2    ,KMESPC  ,MEMNL4  ,FACHFR  )    
!  (This subroutine has not been tested yet)		      
!
!*******************************************************************
!
   USE SWCOMM3                                                         
   USE SWCOMM4                                                         
   USE OCPCOMM4                                                        
   USE M_SNL4 
!   USE ALL_VARS, ONLY : MT,AC2                
   USE VARS_WAVE, ONLY : MT,AC2                
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering                              |
!     | Fluid Mechanics Section                                   |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmers: H.L. Tolman, R.C. Ris                        |
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
!     40.17: IJsbrand Haagsma
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     30.72, Feb. 98: Introduced generic names XCGRID, YCGRID and SPCSIG for SWAN
!     40.17, Dec. 01: Implemented Multiple DIA
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Calculate non-linear interaction using the discrete interaction
!     approximation (Hasselmann and Hasselmann 1985; WAMDI group 1988)
!     for the full circle (option if a current is present). Note: using
!     this subroutine requires an additional array with size
!     (MXC*MYC*MDC*MSC). This requires more internal memory but can
!     speed up the computations sigificantly if a current is present.
!
!  3. Method
!
!     Discrete interaction approximation. To make interpolation simple,
!     the interactions are calculated in a "folded" space.
!
!                            Frequencies -->
!                 +---+---------------------+---------+- IDHGH
!              d  | 3 :          2          :    2    |
!              i  + - + - - - - - - - - - - + - - - - +- MDC
!              r  |   :                     :         |
!              e  | 3 :  original spectrum  :    1    |
!              c  |   :                     :         |
!              t. + - + - - - - - - - - - - + - - - - +- 1
!                 | 3 :          2          :    2    |
!                 +---+---------------------+---------+- IDLOW
!                 |   |                     |     ^   |
!              ISLOW  1                    MSC    |   ISHGH
!                     |                           |
!                   ISCLW                        ISCHG
!              lowest discrete               highest discrete
!                central bin                   central bin
!
!                            1 : Extra tail added beyond MSC
!                            2 : Spectrum copied outside ID range
!                            3 : Empty bins at low frequencies
!
!     ISLOW =  1  + ISM1
!     ISHGH = MSC + ISP1 - ISM1
!     ISCLW =  1
!     ISCHG = MSC - ISM1
!     IDLOW =  1  - MAX(IDM1,IDP1)
!     IDHGH = MDC + MAX(IDM1,IDP1)
!
!       Relative offsets of interpolation points around central bin
!       "#" and corresponding numbers of AWGn :
!
!               ISM1  ISM
!                5        7    T |
!          IDM1   +------+     H +
!                 |      |     E |      ISP      ISP1
!                 |   \  |     T |       3           1
!           IDM   +------+     A +        +---------+  IDP1
!                6       \8      |        |         |
!                                |        |  /      |
!                           \    +        +---------+  IDP
!                                |      /4           2
!                              \ |  /
!          -+-----+------+-------#--------+---------+----------+
!                                |           FREQ.
!
!
!  4. Argument variables
!
!     MCGRD : number of wet grid points of the computational grid
!     MDC   : grid points in theta-direction of computational grid
!     MDC4MA: highest array counter in directional space (Snl4)
!     MDC4MI: lowest array counter in directional space (Snl4)
!     MSC   : grid points in sigma-direction of computational grid
!     MSC4MA: highest array counter in frequency space (Snl4)
!     MSC4MI: lowest array counter in frequency space (Snl4)
!     WWINT : counters for quadruplet interactions
!
   INTEGER WWINT(*)
!
!     AC2   : action density
!     AF11  : scaling frequency
!     DAL1  : coefficient for the quadruplet interactions
!     DAL2  : coefficient for the quadruplet interactions
!     DAL3  : coefficient for the quadruplet interactions
!     DEP2  : depth
!     FACHFR
!     KMESPC: mean average wavenumber over full spectrum
!     MEMNL4
!     PI    : circular constant
!     SA1   : interaction contribution of first quadruplet (unfolded space)
!     SA2   : interaction contribution of second quadruplet (unfolded space)
!     SFNL
!     SNLC1
!     SPCSIG: relative frequencies in computational domain in sigma-space
!     UE    : "unfolded" spectrum
!     WWAWG : weight coefficients for the quadruplet interactions
!
   REAL    DAL1, DAL2, DAL3, FACHFR, KMESPC, SNLC1                     
!   REAL    AC2(MDC,MSC,0:MT)
   REAL    DEP2(MT)
   REAL    MEMNL4(MDC,MSC,MT)
   REAL    SA1(MSC4MI:MSC4MA,MDC4MI:MDC4MA)
   REAL    SA2(MSC4MI:MSC4MA,MDC4MI:MDC4MA)
   REAL    SFNL(MSC4MI:MSC4MA,MDC4MI:MDC4MA)
   REAL    SPCSIG(MSC)                                                 
   REAL    UE(MSC4MI:MSC4MA,MDC4MI:MDC4MA)
   REAL    WWAWG(*)
!
!  7. Common blocks used
!
!
!  8. Subroutines used
!
!     ---
!
!  9. Subroutines calling
!
!     SOURCE (in SWANCOM1)
!
! 12. Structure
!
!     -------------------------------------------
!       Initialisations.
!       Calculate proportionality constant.
!       Prepare auxiliary spectrum.
!       Calculate (unfolded) interactions :
!       -----------------------------------------
!         Energy at interacting bins
!         Contribution to interactions
!         Fold interactions to side angles
!       -----------------------------------------
!       Put source term together
!     -------------------------------------------
!
! 13. Source text
!
!*******************************************************************
!
   INTEGER   IS      ,ID      ,ID0     ,I       ,J       ,          &
             ISHGH   ,IDLOW   ,IDHGH   ,ISP     ,ISP1    ,          &
	     IDP     ,IDP1    ,ISM     ,ISM1    ,IDM     ,IDM1    , &
	     ISCLW   ,ISCHG
!
   REAL      X       ,X2      ,CONS    ,FACTOR  ,SNLCS2  ,          &
             SNLCS3  ,E00     ,EP1     ,EM1     ,EP2     ,EM2     , &
	     SA1A    ,SA1B    ,SA2A    ,SA2B    ,                   &
	     AWG1    ,AWG2    ,AWG3    ,AWG4    ,AWG5    ,AWG6    , &
	     AWG7    ,AWG8    ,                                     &
	     JACOBI  ,SIGPI
!
!   SAVE IENT
!   DATA IENT/0/
!   IF (LTRACE) CALL STRACE (IENT,'SWSNL3')
!
   IDP    = WWINT(1)
   IDP1   = WWINT(2)
   IDM    = WWINT(3)
   IDM1   = WWINT(4)
   ISP    = WWINT(5)
   ISP1   = WWINT(6)
   ISM    = WWINT(7)
   ISM1   = WWINT(8)
   ISLOW  = WWINT(9)
   ISHGH  = WWINT(10)
   ISCLW  = WWINT(11)
   ISCHG  = WWINT(12)
   IDLOW  = WWINT(13)
   IDHGH  = WWINT(14)
!
   AWG1 = WWAWG(1)
   AWG2 = WWAWG(2)
   AWG3 = WWAWG(3)
   AWG4 = WWAWG(4)
   AWG5 = WWAWG(5)
   AWG6 = WWAWG(6)
   AWG7 = WWAWG(7)
   AWG8 = WWAWG(8)
!
!  *** Initialize auxiliary arrays per gridpoint ***
!
   DO ID = MDC4MI, MDC4MA
     DO IS = MSC4MI, MSC4MA
       UE(IS,ID)   = 0.
       SA1(IS,ID)  = 0.
       SA2(IS,ID)  = 0.
       SFNL(IS,ID) = 0.
     ENDDO
   ENDDO
!
!  *** Calculate prop. constant.                           ***
!  *** Calculate factor R(X) to calculate the NL wave-wave ***
!  *** interaction for shallow water                       ***
!  *** SNLC1 = 1/GRAV**4                                   ***         
!
   SNLCS1 = PQUAD(3)                                                   
   SNLCS2 = PQUAD(4)                                                   
   SNLCS3 = PQUAD(5)                                                   
!   X      = MAX ( 0.75 * DEP2(IGC) * KMESPC , 0.5 )               
   X      = MAX ( 0.75 * DEP2(kcgrd(1)) * KMESPC , 0.5 )               
   X2     = MAX ( -1.E15, SNLCS3*X)
   CONS   = SNLC1 * ( 1. + SNLCS1/X * (1.-SNLCS2*X) * EXP(X2))
   JACOBI = 2. * PI_W
!
!  *** extend the area with action density at periodic boundaries ***
!
   DO IDDUM = IDLOW, IDHGH
     ID = MOD ( IDDUM - 1 + MDC , MDC ) + 1
     DO IS=1, MSC
!       UE (IS,IDDUM) = AC2(ID,IS,IGC) * SPCSIG(IS) * JACOBI       
       UE (IS,IDDUM) = AC2(ID,IS,kcgrd(1)) * SPCSIG(IS) * JACOBI       
     ENDDO
   ENDDO
!
   DO IS = MSC+1, ISHGH
     DO ID = IDLOW, IDHGH
       UE(IS,ID) = UE(IS-1,ID) * FACHFR
     ENDDO
   ENDDO
!
!  *** Calculate (unfolded) interactions ***
!  *** Energy at interacting bins        ***
!
   DO IS = ISCLW, ISCHG
     DO ID = 1, MDC
       E00    =        UE(IS      ,ID      )
       EP1    = AWG1 * UE(IS+ISP1,ID+IDP1) +                       &
                AWG2 * UE(IS+ISP1,ID+IDP ) +                       &
		AWG3 * UE(IS+ISP ,ID+IDP1) +                       &
		AWG4 * UE(IS+ISP ,ID+IDP )
       EM1    = AWG5 * UE(IS+ISM1,ID-IDM1) +                       &
                AWG6 * UE(IS+ISM1,ID-IDM ) +                       &
		AWG7 * UE(IS+ISM ,ID-IDM1) +                       &
		AWG8 * UE(IS+ISM ,ID-IDM )
       EP2    = AWG1 * UE(IS+ISP1,ID-IDP1) +                       &
                AWG2 * UE(IS+ISP1,ID-IDP ) +                       &
		AWG3 * UE(IS+ISP ,ID-IDP1) +                       &
		AWG4 * UE(IS+ISP ,ID-IDP )
       EM2    = AWG5 * UE(IS+ISM1,ID+IDM1) +                       &
                AWG6 * UE(IS+ISM1,ID+IDM ) +                       &
		AWG7 * UE(IS+ISM ,ID+IDM1) +                       &
		AWG8 * UE(IS+ISM ,ID+IDM )
!
!      Contribution to interactions
!
       FACTOR = CONS * AF11(IS) * E00
!
       SA1A   = E00 * ( EP1*DAL1 + EM1*DAL2 ) * PQUAD(2)               
       SA1B   = SA1A - EP1*EM1*DAL3 * PQUAD(2)                         
       SA2A   = E00 * ( EP2*DAL1 + EM2*DAL2 ) * PQUAD(2)               
       SA2B   = SA2A - EP2*EM2*DAL3 * PQUAD(2)                         
!
       SA1 (IS,ID) = FACTOR * SA1B
       SA2 (IS,ID) = FACTOR * SA2B
!
!       IF(ITEST >= 100 .AND. TESTFL)THEN
!         WRITE(PRINTF,9002) E00,EP1,EM1,EP2,EM2
!9002     FORMAT (' E00 EP1 EM1 EP2 EM2  :',5E11.4)
!         WRITE(PRINTF,9003) SA1A,SA1B,SA2A,SA2B
!9003     FORMAT (' SA1A SA1B SA2A SA2B  :',4E11.4)
!         WRITE(PRINTF,9004) IS,ID,SA1(IS,ID),SA2(IS,ID)
!9004     FORMAT (' IS ID SA1() SA2()    :',2I4,2E12.4)
!         WRITE(PRINTF,9005) FACTOR,JACOBI
!9005     FORMAT (' FACTOR JACOBI        : ',2E12.4)
!       END IF
!
     ENDDO
   ENDDO
!
!  *** Fold interactions to side angles -> domain in theta is ***
!  *** periodic                                               ***
!
   DO ID = 1, IDHGH - MDC
     ID0   = 1 - ID
     DO IS = ISCLW, ISCHG
       SA1 (IS,MDC+ID) = SA1 (IS,  ID   )
       SA2 (IS,MDC+ID) = SA2 (IS,  ID   )
       SA1 (IS,  ID0 ) = SA1 (IS,MDC+ID0)
       SA2 (IS,  ID0 ) = SA2 (IS,MDC+ID0)
     ENDDO
   ENDDO
!
!  *** Put source term together (To save space I=IS and ***
!  *** J=MDC is used)  ----                             ***
!
   DO I = 1, MSC
     SIGPI = SPCSIG(I) * JACOBI                                        
     DO J = 1, MDC
       SFNL(I,J) = - 2. * ( SA1(I,J) + SA2(I,J) )                       &
                   + AWG1 * ( SA1(I-ISP1,J-IDP1) + SA2(I-ISP1,J+IDP1) ) &
		   + AWG2 * ( SA1(I-ISP1,J-IDP ) + SA2(I-ISP1,J+IDP ) ) &
		   + AWG3 * ( SA1(I-ISP ,J-IDP1) + SA2(I-ISP ,J+IDP1) ) &
		   + AWG4 * ( SA1(I-ISP ,J-IDP ) + SA2(I-ISP ,J+IDP ) ) &
		   + AWG5 * ( SA1(I-ISM1,J+IDM1) + SA2(I-ISM1,J-IDM1) ) &
		   + AWG6 * ( SA1(I-ISM1,J+IDM ) + SA2(I-ISM1,J-IDM ) ) &
		   + AWG7 * ( SA1(I-ISM ,J+IDM1) + SA2(I-ISM ,J-IDM1) ) &
		   + AWG8 * ( SA1(I-ISM ,J+IDM ) + SA2(I-ISM ,J-IDM ) )
!
!      *** store value in auxiliary array and use values in ***
!      *** next four sweeps (see subroutine FILNL3)         ***
!
!       MEMNL4(J,I,IGC) = SFNL(I,J) / SIGPI                        
       MEMNL4(J,I,kcgrd(1)) = SFNL(I,J) / SIGPI                        
     ENDDO
   ENDDO
!
!  *** test output ***
!
!   IF(ITEST >= 50 .AND. TESTFL)THEN
!     WRITE(PRINTF,*)
!     WRITE(PRINTF,*) ' SWSNL3 subroutine '
!     WRITE(PRINTF,9011) IDP, IDP1, IDM, IDM1
!9011 FORMAT (' IDP IDP1 IDM IDM1     :',4I5)
!     WRITE (PRINTF,9013) ISP, ISP1, ISM, ISM1
!9013 FORMAT (' ISP ISP1 ISM ISM1     :',4I5)
!     WRITE (PRINTF,9015) ISLOW, ISHGH, IDLOW, IDHGH
!9015 FORMAT (' ISLOW ISHG IDLOW IDHG :',4I5)
!     WRITE(PRINTF,9016) ISCLW, ISCHG, JACOBI
!9016 FORMAT (' ICLW ICHG JACOBI      :',2I5,E12.4)
!     WRITE (PRINTF,9017) AWG1, AWG2, AWG3, AWG4
!9017 FORMAT (' AWG1 AWG2 AWG3 AWG4   :',4E12.4)
!     WRITE (PRINTF,9018) AWG5, AWG6, AWG7, AWG8
!9018 FORMAT (' AWG5 AWG6 AWG7 AWG8   :',4E12.4)
!     WRITE (PRINTF,9019) MSC4MI, MSC4MA, MDC4MI, MDC4MA
!9019 FORMAT (' S4MI S4MA D4MI D4MA   :',4I6)
!     WRITE(PRINTF,9020) SNLC1,X,X2,CONS
!9020 FORMAT (' SNLC1  X  X2  CONS    :',4E12.4)
!    WRITE(PRINTF,9021) DEP2(KCGRD(1)),KMESPC,FACHFR,PI
!9021 FORMAT (' DEPTH KMESPC FACHFR PI:',4E12.4)
!     WRITE(PRINTF,*)
!
!    *** value source term in every bin ***
!
!     IF(ITEST >= 150)THEN
!       DO I=1, MSC
!         DO J=1, MDC
!           WRITE(PRINTF,2006) I,J,MEMNL4(J,I,KCGRD(1)),SFNL(I,J),      &
!	                      SPCSIG(I)                                
!2006       FORMAT (' I J MEMNL() SFNL() SPCSIG:',2I4,3E12.4)           
!         ENDDO
!       ENDDO
!     END IF
!   END IF
!
   RETURN
!
   END SUBROUTINE SWSNL3

!
!*******************************************************************
!
   SUBROUTINE SWSNL2 (IDDLOW  ,IDDTOP  ,WWINT   ,WWAWG   ,UE      , &
                      SA1     ,ISSTOP  ,SA2     ,SPCSIG  ,SNLC1   , &
		      DAL1    ,DAL2    ,DAL3    ,SFNL    ,DEP2    , &
		      KMESPC  ,IMATDA  ,IMATRA  ,FACHFR  ,PLNL4S  , &
		      IDCMIN  ,IDCMAX  ,IG      ,CGO     ,WWSWG   )  
!
!*******************************************************************
!
!     Calculate non-linear interaction using the discrete interaction
!     approximation (Hasselmann and Hasselmann 1985; WAMDI group 1988)
!
!*******************************************************************
   USE SWCOMM3                                                         
   USE SWCOMM4                                                         
   USE OCPCOMM4                                                        
   USE M_SNL4  
!   USE ALL_VARS, ONLY : MT,AC2                   
   USE VARS_WAVE, ONLY : MT,AC2                   
   
   IMPLICIT NONE                                     
!
   REAL    :: SPCSIG(MSC)                                                 
   INTEGER :: IS     ,ID     ,I      ,J      ,IG     ,ISHGH  ,  &
              ISSTOP ,ISP    ,ISP1   ,IDP    ,IDP1   ,ISM    ,ISM1   ,  &
	      IDM    ,IDM1   ,ISCLW  ,ISCHG  ,                          &
	      IDLOW  ,IDHGH  ,IDDLOW ,IDDTOP ,IDCLOW ,IDCHGH    
!
   REAL    :: X      ,X2     ,CONS   ,FACTOR ,SNLCS1 ,SNLCS2 ,SNLCS3 ,  &
              E00    ,EP1    ,EM1    ,EP2    ,EM2    ,SA1A   ,SA1B   ,  &
	      SA2A   ,SA2B   ,KMESPC ,FACHFR ,AWG1   ,AWG2   ,AWG3   ,  &
	      AWG4   ,AWG5   ,AWG6   ,AWG7   ,AWG8   ,DAL1   ,DAL2   ,  &
	      DAL3           ,JACOBI ,SIGPI                            
!
    REAL   :: DEP2(MT)                              ,                   &
	      UE(MSC4MI:MSC4MA , MDC4MI:MDC4MA )    ,                   &
	      SA1(MSC4MI:MSC4MA , MDC4MI:MDC4MA )   ,                   &
	      SA2(MSC4MI:MSC4MA , MDC4MI:MDC4MA )   ,                   &
	      DA1C(MSC4MI:MSC4MA , MDC4MI:MDC4MA )  ,                   &
	      DA1P(MSC4MI:MSC4MA , MDC4MI:MDC4MA )  ,                   &
	      DA1M(MSC4MI:MSC4MA , MDC4MI:MDC4MA )  ,                   &
	      DA2C(MSC4MI:MSC4MA , MDC4MI:MDC4MA )  ,                   &
	      DA2P(MSC4MI:MSC4MA , MDC4MI:MDC4MA )  ,                   &
	      DA2M(MSC4MI:MSC4MA , MDC4MI:MDC4MA )  ,                   &
	      SFNL(MSC4MI:MSC4MA , MDC4MI:MDC4MA)   ,                   &
	      DFNL(MSC4MI:MSC4MA , MDC4MI:MDC4MA)   ,                   &
	      IMATRA(MDC,MSC)                       ,                   &
	      IMATDA(MDC,MSC)                       ,                   &
	      PLNL4S(MDC,MSC,NPTST)                 ,                   &
	      WWAWG(*)                              ,                   &
	      CGO(MSC,MICMAX)
!
   INTEGER :: WWINT(*)         ,IDCMIN(MSC)      ,IDCMAX(MSC)
   INTEGER :: ISLOW,IIID,IDDUM,ID0
   REAL    :: SNLC1,SWG1,SWG2,SWG3,SWG4,SWG5,SWG6,SWG7,SWG8
   REAL    :: WWSWG(*),PI3
   
!
   LOGICAL   PERCIR
!
   IDP    = WWINT(1)
   IDP1   = WWINT(2)
   IDM    = WWINT(3)
   IDM1   = WWINT(4)
   ISP    = WWINT(5)
   ISP1   = WWINT(6)
   ISM    = WWINT(7)
   ISM1   = WWINT(8)
   ISLOW  = WWINT(9)
   ISHGH  = WWINT(10)
   ISCLW  = WWINT(11)
   ISCHG  = WWINT(12)
   IDLOW  = WWINT(13)
   IDHGH  = WWINT(14)

   AWG1 = WWAWG(1)
   AWG2 = WWAWG(2)
   AWG3 = WWAWG(3)
   AWG4 = WWAWG(4)
   AWG5 = WWAWG(5)
   AWG6 = WWAWG(6)
   AWG7 = WWAWG(7)
   AWG8 = WWAWG(8)
!
   SWG1 = WWSWG(1)
   SWG2 = WWSWG(2)
   SWG3 = WWSWG(3)
   SWG4 = WWSWG(4)
   SWG5 = WWSWG(5)
   SWG6 = WWSWG(6)
   SWG7 = WWSWG(7)
   SWG8 = WWSWG(8)
!
!  *** Initialize auxiliary arrays per gridpoint ***
!
   DO ID = MDC4MI, MDC4MA
     DO IS = MSC4MI, MSC4MA
       UE(IS,ID)   = 0.
       SA1(IS,ID)  = 0.
       SA2(IS,ID)  = 0.
       
       DA1C(IS,ID) = 0.
       DA1P(IS,ID) = 0.
       DA1M(IS,ID) = 0.
       DA2C(IS,ID) = 0.
       DA2P(IS,ID) = 0.
       DA2M(IS,ID) = 0.
       
       SFNL(IS,ID) = 0.
       DFNL(IS,ID) = 0.
       
     ENDDO
   ENDDO
!
!  *** Calculate prop. constant.                           ***
!  *** Calculate factor R(X) to calculate the NL wave-wave ***
!  *** interaction for shallow water                       ***
!  *** SNLC1 = 1/GRAV**4                                   ***         
!
   SNLCS1 = PQUAD(3)                                                   
   SNLCS2 = PQUAD(4)                                                   
   SNLCS3 = PQUAD(5)    

   X      = MAX ( 0.75 * DEP2(IG) * KMESPC , 0.5 )
   X2     = MAX ( -1.E15, SNLCS3*X)
   CONS   = SNLC1 * ( 1. + SNLCS1/X * (1.-SNLCS2*X) * EXP(X2))
   JACOBI = 2. * PI_W
!
!  *** check whether the spectral domain is periodic in ***
!  *** direction space and if so modify boundaries      ***
!
   PERCIR = .FALSE.
   IF(IDDLOW == 1 .AND. IDDTOP == MDC)THEN
!    *** periodic in theta -> spectrum can be folded  ***
!    *** (can only occur in presence of a current)    ***
     IDCLOW = 1
     IDCHGH = MDC
     IIID   = 0
     PERCIR = .TRUE.
   ELSE
!    *** different sectors per sweep -> extend range with IIID ***
     IIID   = MAX ( IDM1 , IDP1 )
     IDCLOW = IDLOW
     IDCHGH = IDHGH
   ENDIF
!
!  *** Prepare auxiliary spectrum               ***
!  *** set action original spectrum in array UE ***
!
!print*,'swsnl2-1'
!print*,'IDLOW,IDHGH,IIID',IDLOW,IDHGH,IIID
!print*,'low high MDC=',IDLOW - IIID, IDHGH + IIID,MDC

   DO IDDUM = IDLOW - IIID , IDHGH + IIID
     ID = MOD ( IDDUM - 1 + MDC , MDC ) + 1
     DO IS = 1, MSC
!      print*,'ID,IS,IG=',ID,IS,IG
       UE(IS,IDDUM) = AC2(ID,IS,IG) * SPCSIG(IS) * JACOBI      
     ENDDO
   ENDDO
!print*,'swsnl2-2'
!
!  *** set values in the areas 2 for IS > MSC+1 ***
!
   DO IS = MSC+1, ISHGH
     DO ID = IDLOW - IIID , IDHGH + IIID
       UE (IS,ID) = UE(IS-1,ID) * FACHFR
     ENDDO
   ENDDO
!
!  *** Calculate interactions      ***
!  *** Energy at interacting bins  ***
!
   DO IS = ISCLW, ISCHG
     DO ID = IDCLOW , IDCHGH
       E00 =        UE(IS     ,ID      )
       EP1 = AWG1 * UE(IS+ISP1,ID+IDP1) +                   &
             AWG2 * UE(IS+ISP1,ID+IDP ) +                   &
	     AWG3 * UE(IS+ISP ,ID+IDP1) +                   &
	     AWG4 * UE(IS+ISP ,ID+IDP )
       EM1 = AWG5 * UE(IS+ISM1,ID-IDM1) +                   &
             AWG6 * UE(IS+ISM1,ID-IDM ) +                   &
	     AWG7 * UE(IS+ISM ,ID-IDM1) +                   &
	     AWG8 * UE(IS+ISM ,ID-IDM )
!
       EP2 = AWG1 * UE(IS+ISP1,ID-IDP1) +                   &
             AWG2 * UE(IS+ISP1,ID-IDP ) +                   &
	     AWG3 * UE(IS+ISP ,ID-IDP1) +                   &
	     AWG4 * UE(IS+ISP ,ID-IDP )
       EM2 = AWG5 * UE(IS+ISM1,ID+IDM1) +                   &
             AWG6 * UE(IS+ISM1,ID+IDM ) +                   &
	     AWG7 * UE(IS+ISM ,ID+IDM1) +                   &
	     AWG8 * UE(IS+ISM ,ID+IDM )
!
!      *** Contribution to interactions                          ***
!      *** CONS is the shallow water factor for the NL interact. ***
!
       FACTOR = CONS * AF11(IS) * E00
!
       SA1A = E00 * ( EP1*DAL1 + EM1*DAL2 ) * PQUAD(2)               
       SA1B = SA1A - EP1*EM1*DAL3 * PQUAD(2)                         
       SA2A = E00 * ( EP2*DAL1 + EM2*DAL2 ) * PQUAD(2)               
       SA2B = SA2A - EP2*EM2*DAL3 * PQUAD(2)                         
!
       SA1 (IS,ID) = FACTOR * SA1B
       SA2 (IS,ID) = FACTOR * SA2B
       
       DA1C(IS,ID) = CONS*AF11(IS)*(SA1A+SA1B)
       DA1P(IS,ID) = FACTOR*(DAL1*E00-DAL3*EM1)*PQUAD(2)
       DA1M(IS,ID) = FACTOR*(DAL2*E00-DAL3*EP1)*PQUAD(2)
       
       DA2C(IS,ID) = CONS*AF11(IS)*(SA2A+SA2B)
       DA2P(IS,ID) = FACTOR*(DAL1*E00-DAL3*EM2)*PQUAD(2)
       DA2M(IS,ID) = FACTOR*(DAL2*E00-DAL3*EP2)*PQUAD(2)
       
     ENDDO
   ENDDO
!
!  *** Fold interactions to side angles if spectral domain ***
!  *** is periodic in directional space                    ***
!
   IF(PERCIR)THEN
     DO ID = 1, IDHGH - MDC
       ID0   = 1 - ID
       DO IS = ISCLW, ISCHG
         SA1 (IS,MDC+ID) = SA1 (IS ,  ID    )
         SA2 (IS,MDC+ID) = SA2 (IS ,  ID    )
         SA1 (IS,  ID0 ) = SA1 (IS , MDC+ID0)
         SA2 (IS,  ID0 ) = SA2 (IS , MDC+ID0)
       
       DA1C(IS,MDC+ID) = DA1C(IS,ID)
       DA1P(IS,MDC+ID) = DA1P(IS,ID)
       DA1M(IS,MDC+ID) = DA1M(IS,ID)
       DA1C(IS,ID0) = DA1C(IS,MDC+ID0)
       DA1P(IS,ID0) = DA1P(IS,MDC+ID0)
       DA1M(IS,ID0) = DA1M(IS,MDC+ID0)
       
       DA2C(IS,MDC+ID) = DA2C(IS,ID)
       DA2P(IS,MDC+ID) = DA2P(IS,ID)
       DA2M(IS,MDC+ID) = DA2M(IS,ID)
       DA2C(IS,ID0) = DA2C(IS,MDC+ID0)
       DA2P(IS,ID0) = DA2P(IS,MDC+ID0)
       DA2M(IS,ID0) = DA2M(IS,MDC+ID0)
       ENDDO
     ENDDO
   ENDIF
!
!  ***  Put source term together (To save space I=IS and J=ID ***
!  ***  is used)                                              ***
!
   PI3 = (2.0*PI_W)**3
   DO I = 1, ISSTOP
     SIGPI = SPCSIG(I) * JACOBI 
     
     DO J = IDCMIN(I), IDCMAX(I)
       ID = MOD ( J - 1 + MDC , MDC ) + 1
       SFNL(I,ID) = - 2. * ( SA1(I,J) + SA2(I,J) )                       &
                    + AWG1 * ( SA1(I-ISP1,J-IDP1) + SA2(I-ISP1,J+IDP1) ) &
		    + AWG2 * ( SA1(I-ISP1,J-IDP ) + SA2(I-ISP1,J+IDP ) ) &
		    + AWG3 * ( SA1(I-ISP ,J-IDP1) + SA2(I-ISP ,J+IDP1) ) &
		    + AWG4 * ( SA1(I-ISP ,J-IDP ) + SA2(I-ISP ,J+IDP ) ) &
		    + AWG5 * ( SA1(I-ISM1,J+IDM1) + SA2(I-ISM1,J-IDM1) ) &
		    + AWG6 * ( SA1(I-ISM1,J+IDM ) + SA2(I-ISM1,J-IDM ) ) &
		    + AWG7 * ( SA1(I-ISM ,J+IDM1) + SA2(I-ISM ,J-IDM1) ) &
		    + AWG8 * ( SA1(I-ISM ,J+IDM ) + SA2(I-ISM ,J-IDM ) )

       DFNL(I,ID) =   - 2. * ( DA1C(I,J) + DA2C(I,J) )                      &
                    + SWG1 * ( DA1P(I-ISP1,J-IDP1) + DA2P(I-ISP1,J+IDP1) )  &
		    + SWG2 * ( DA1P(I-ISP1,J-IDP ) + DA2P(I-ISP1,J+IDP ) )  &
		    + SWG3 * ( DA1P(I-ISP ,J-IDP1) + DA2P(I-ISP ,J+IDP1) )  &
		    + SWG4 * ( DA1P(I-ISP ,J-IDP ) + DA2P(I-ISP ,J+IDP ) )  &
		    + SWG5 * ( DA1M(I-ISM1,J+IDM1) + DA2M(I-ISM1,J-IDM1) )  &
		    + SWG6 * ( DA1M(I-ISM1,J+IDM ) + DA2M(I-ISM1,J-IDM ) )  &
		    + SWG7 * ( DA1M(I-ISM ,J+IDM1) + DA2M(I-ISM ,J-IDM1) )  &
		    + SWG8 * ( DA1M(I-ISM ,J+IDM ) + DA2M(I-ISM ,J-IDM ) )
		    
!
!      *** store results in rhv ***
!      *** store results in rhs and main diagonal according ***        
!      *** to Patankar-rules                                ***        
!
       IF(TESTFL) PLNL4S(ID,I,IPTST) =  SFNL(I,ID) / SIGPI  
       IMATRA(ID,I) = IMATRA(ID,I) + SFNL(I,ID) / SIGPI
       IMATDA(ID,I) = IMATDA(ID,I) + DFNL(I,ID)     
	 
     ENDDO
   ENDDO

   RETURN
   END SUBROUTINE SWSNL2

!
!********************************************************************
!
   SUBROUTINE SWSNL1 (WWINT   ,WWAWG   ,WWSWG   ,IDCMIN  ,IDCMAX  , &
                      UE      ,SA1     ,SA2     ,DA1C    ,DA1P    , &
		      DA1M    ,DA2C    ,DA2P    ,DA2M    ,SPCSIG  , &
		      SNLC1   ,KMESPC  ,FACHFR  ,ISSTOP  ,DAL1    , &
		      DAL2    ,DAL3    ,SFNL    ,DSNL    ,DEP2    , &
		      AC2     ,IMATDA  ,IMATRA  ,PLNL4S  ,PLNL4D  , &
		      IDDLOW  ,IDDTOP  )                           

!(This subroutine has not been tested yet and not useful any more)
!
!********************************************************************
!
   USE SWCOMM3                                                         
   USE SWCOMM4                                                         
   USE OCPCOMM4                                                        
   USE M_SNL4   
   USE ALL_VARS, ONLY : MT                                                       
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering                              |
!     | Fluid Mechanics Section                                   |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmers: H.L. Tolman, R.C. Ris                        |
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
!     40.13: Nico Booij
!     40.17: IJsbrand Haagsma
!     40.23: Marcel Zijlema
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     30.72, Feb. 98: Introduced generic names XCGRID, YCGRID and SPCSIG for SWAN
!     40.17, Dec. 01: Implentation of Multiple DIA
!     40.23, Aug. 02: some corrections
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Calculate non-linear interaction using the discrete interaction
!     approximation (Hasselmann and Hasselmann 1985; WAMDI group 1988),
!     including the diagonal term for the implicit integration.
!
!     The interactions are calculated for all bin's that fall
!     within a sweep. No additional auxiliary array is required (see
!     SWSNL3)
!
!  3. Method
!
!     Discrete interaction approximation.
!
!     Since the domain in directional domain is by definition not
!     periodic, the spectral space can not beforehand
!     folded to the side angles. This can only be done if the
!     full circle has to be calculated
!
!
!                            Frequencies -->
!                 +---+---------------------+---------+- IDHGH
!              d  | 3 :          2          :    2    |
!              i  + - + - - - - - - - - - - + - - - - +- MDC
!              r  |   :                     :         |
!              e  | 3 :  original spectrum  :    1    |
!              c  |   :                     :         |
!              t. + - + - - - - - - - - - - + - - - - +- 1
!                 | 3 :          2          :    2    |
!                 +---+---------------------+---------+- IDLOW
!                 |   |                     |    ^    |
!             ISLOW   1                     MSC  |    ISHGH
!                     ^                          |
!                     |                          |
!                    ISCLW                     ISCHG
!              lowest discrete               highest discrete
!                central bin                   central bin
!
!                            1 : Extra tail added beyond MSC
!                            2 : Spectrum copied outside ID range
!                            3 : Empty bins at low frequencies
!
!     ISLOW =  1  + ISM1
!     ISHGH = MSC + ISP1 - ISM1
!     ISCLW =  1
!     ISCHG = MSC - ISM1
!     IDLOW =  IDDLOW - MAX(IDM1,IDP1)
!     IDHGH =  IDDTOP + MAX(IDM1,IDP1)
!
!     For the meaning of the counters on the right hand side of the
!     above equations see section 4.
!
!  4. Argument variables
!
!     SPCSIG: Relative frequencies in computational domain in sigma-space 
!
   REAL    SPCSIG(MSC)                                                 
!
!     Data in PARAMETER statements :
!     ----------------------------------------------------------------
!       DAL1    Real  LAMBDA dependend weight factors (see FAC4WW)
!       DAL2    Real
!       DAL3    Real
!       ITHP, ITHP1, ITHM, ITHM1, IFRP, IFRP1, IFRM, IFRM1
!               Int.  Counters of interpolation point relative to
!                     central bin, see figure below (set in FAC4WW).
!       NFRLOW, NFRHGH, NFRCHG, NTHLOW, NTHHGH
!               Int.  Range of calculations, see section 2.
!       AF11    R.A.  Scaling array (Freq**11).
!       AWGn    Real  Interpolation weights, see numbers in fig.
!       SWGn    Real  Id. squared.
!       UE      R.A.  "Unfolded" spectrum.
!       SA1     R.A.  Interaction constribution of first and second
!       SA2     R.A.    quadr. respectively (unfolded space).
!       DA1C, DA1P, DA1M, DA2C, DA2P, DA2M
!               R.A.  Idem for diagonal matrix.
!       PERCIR        full circle or sector
!     ----------------------------------------------------------------
!
!       Relative offsets of interpolation points around central bin
!       "#" and corresponding numbers of AWGn :
!
!               ISM1  ISM
!                5        7    T |
!          IDM1   +------+     H +
!                 |      |     E |      ISP      ISP1
!                 |   \  |     T |       3           1
!           IDM   +------+     A +        +---------+  IDP1
!                6       \8      |        |         |
!                                |        |  /      |
!                           \    +        +---------+  IDP
!                                |      /4           2
!                              \ |  /
!          -+-----+------+-------#--------+---------+----------+
!                                |           FREQ.
!
!  7. Common blocks used
!
!
!  8. Subroutines used
!
!     ---
!
!  9. Subroutines calling
!
!     SOURCE (in SWANCOM1)
!
! 12. Structure
!
!     -------------------------------------------
!       Initialisations.
!       Calculate proportionality constant.
!       Prepare auxiliary spectrum.
!       Calculate interactions :
!       -----------------------------------------
!         Energy at interacting bins
!         Contribution to interactions
!         Fold interactions to side angles
!       -----------------------------------------
!       Put source term together
!     -------------------------------------------
!
! 13. Source text
!
!*************************************************************
!
   INTEGER   IS     ,ID     ,I      ,J      ,                       &
             ISHGH  ,IDLOW  ,ISP    ,ISP1   ,IDP    ,IDP1   ,       &
	     ISM    ,ISM1   ,IDHGH  ,IDM    ,IDM1   ,ISCLW  ,       &
	     ISCHG  ,IDDLOW ,IDDTOP                                    
!
   REAL      X      ,X2     ,CONS   ,FACTOR ,SNLCS1 ,SNLCS2 ,SNLCS3, &
             E00    ,EP1    ,EM1    ,EP2    ,EM2    ,SA1A   ,SA1B  , &
	     SA2A   ,SA2B   ,KMESPC ,FACHFR ,AWG1   ,AWG2   ,AWG3  , &
	     AWG4   ,AWG5   ,AWG6   ,AWG7   ,AWG8   ,DAL1   ,DAL2  , &
	     DAL3   ,SNLC1  ,SWG1   ,SWG2   ,SWG3   ,SWG4   ,SWG5  , &
	     SWG6   ,SWG7   ,SWG8           ,JACOBI ,SIGPI            
!
   REAL      AC2(MDC,MSC,0:MT)                    ,                 &
             DEP2(MT)                           ,                 &
	     UE(MSC4MI:MSC4MA , MDC4MI:MDC4MA )    ,                 &
	     SA1(MSC4MI:MSC4MA , MDC4MI:MDC4MA )   ,                 &
	     SA2(MSC4MI:MSC4MA , MDC4MI:MDC4MA )   ,                 &
	     DA1C(MSC4MI:MSC4MA , MDC4MI:MDC4MA )  ,                 &
	     DA1P(MSC4MI:MSC4MA , MDC4MI:MDC4MA )  ,                 &
	     DA1M(MSC4MI:MSC4MA , MDC4MI:MDC4MA )  ,                 &
	     DA2C(MSC4MI:MSC4MA , MDC4MI:MDC4MA )  ,                 &
	     DA2P(MSC4MI:MSC4MA , MDC4MI:MDC4MA )  ,                 &
	     DA2M(MSC4MI:MSC4MA , MDC4MI:MDC4MA )  ,                 &
	     SFNL(MSC4MI:MSC4MA , MDC4MI:MDC4MA )  ,                 &
	     DSNL(MSC4MI:MSC4MA , MDC4MI:MDC4MA )  ,                 &
	     IMATDA(MDC,MSC)                       ,                 &
	     IMATRA(MDC,MSC)                       ,                 &
	     PLNL4S(MDC,MSC,NPTST)                 ,                 &
	     PLNL4D(MDC,MSC,NPTST)                 ,                 &
	     WWAWG(*)                              ,                 &
	     WWSWG(*)
!
   INTEGER   IDCMIN(MSC),IDCMAX(MSC),WWINT(*)
!
   LOGICAL   PERCIR
!
!   SAVE IENT
!   DATA IENT/0/
!   IF (LTRACE) CALL STRACE (IENT,'SWSNL1')
!
   IDP    = WWINT(1)
   IDP1   = WWINT(2)
   IDM    = WWINT(3)
   IDM1   = WWINT(4)
   ISP    = WWINT(5)
   ISP1   = WWINT(6)
   ISM    = WWINT(7)
   ISM1   = WWINT(8)
   ISLOW  = WWINT(9)
   ISHGH  = WWINT(10)
   ISCLW  = WWINT(11)
   ISCHG  = WWINT(12)
   IDLOW  = WWINT(13)
   IDHGH  = WWINT(14)
!
   AWG1 = WWAWG(1)
   AWG2 = WWAWG(2)
   AWG3 = WWAWG(3)
   AWG4 = WWAWG(4)
   AWG5 = WWAWG(5)
   AWG6 = WWAWG(6)
   AWG7 = WWAWG(7)
   AWG8 = WWAWG(8)
!
   SWG1 = WWSWG(1)
   SWG2 = WWSWG(2)
   SWG3 = WWSWG(3)
   SWG4 = WWSWG(4)
   SWG5 = WWSWG(5)
   SWG6 = WWSWG(6)
   SWG7 = WWSWG(7)
   SWG8 = WWSWG(8)
!
!  *** Initialize auxiliary arrays per gridpoint ***
!
   DO ID = MDC4MI, MDC4MA
     DO IS = MSC4MI, MSC4MA
       UE(IS,ID)   = 0.
       SA1(IS,ID)  = 0.
       SA2(IS,ID)  = 0.
       SFNL(IS,ID) = 0.
       DA1C(IS,ID) = 0.
       DA1P(IS,ID) = 0.
       DA1M(IS,ID) = 0.
       DA2C(IS,ID) = 0.
       DA2P(IS,ID) = 0.
       DA2M(IS,ID) = 0.
       DSNL(IS,ID) = 0.
     ENDDO
   ENDDO
!
!  *** Calculate factor R(X) to calculate the NL wave-wave ***
!  *** interaction for shallow water                       ***
!  *** SNLC1 = 1/GRAV**4                                   ***         
!
   SNLCS1 = PQUAD(3)                                                   
   SNLCS2 = PQUAD(4)                                                   
   SNLCS3 = PQUAD(5)                                                   
!   X      = MAX ( 0.75 * DEP2(IGC) * KMESPC , 0.5 )
   X      = MAX ( 0.75 * DEP2(kcgrd(1)) * KMESPC , 0.5 )
   X2     = MAX ( -1.E15, SNLCS3*X)
   CONS   = SNLC1 * ( 1. + SNLCS1/X * (1.-SNLCS2*X) * EXP(X2))
   JACOBI = 2. * PI_W
!
!  *** check whether the spectral domain is periodic in ***
!  *** directional space and if so, modify boundaries   ***
!
   PERCIR = .FALSE.
   IF(IDDLOW == 1 .AND. IDDTOP == MDC)THEN
!    *** periodic in theta -> spectrum can be folded    ***
!    *** (can only be present in presence of a current) ***
     IDCLOW = 1
     IDCHGH = MDC
     IIID   = 0
     PERCIR = .TRUE.
   ELSE
!    *** different sectors per sweep -> extend range with IIID ***
     IIID   = MAX ( IDM1 , IDP1 )
     IDCLOW = IDLOW
     IDCHGH = IDHGH
   ENDIF
!
!  *** Prepare auxiliary spectrum               ***
!  *** set action original spectrum in array UE ***
!
   DO IDDUM = IDLOW - IIID, IDHGH + IIID
     ID = MOD ( IDDUM - 1 + MDC , MDC ) + 1
     DO IS = 1, MSC
!       UE(IS,IDDUM) = AC2(ID,IS,IGC) * SPCSIG(IS) * JACOBI        
       UE(IS,IDDUM) = AC2(ID,IS,kcgrd(1)) * SPCSIG(IS) * JACOBI        
     ENDDO
   ENDDO
!
!  *** set values in area 2 for IS > MSC+1  ***
!
   DO IS = MSC+1, ISHGH
     DO ID = IDLOW - IIID , IDHGH + IIID
       UE (IS,ID) = UE(IS-1,ID) * FACHFR
     ENDDO
   ENDDO
!
!  *** Calculate interactions      ***
!  *** Energy at interacting bins  ***
!
   DO IS = ISCLW, ISCHG
     DO ID = IDCLOW, IDCHGH
       E00    =        UE(IS      ,ID      )
       EP1    = AWG1 * UE(IS+ISP1,ID+IDP1) +                  &
                AWG2 * UE(IS+ISP1,ID+IDP ) +                  &
		AWG3 * UE(IS+ISP ,ID+IDP1) +                  &
		AWG4 * UE(IS+ISP ,ID+IDP )
       EM1    = AWG5 * UE(IS+ISM1,ID-IDM1) +                  &
                AWG6 * UE(IS+ISM1,ID-IDM ) +                  &
		AWG7 * UE(IS+ISM ,ID-IDM1) +                  &
		AWG8 * UE(IS+ISM ,ID-IDM )
!
       EP2    = AWG1 * UE(IS+ISP1,ID-IDP1) +                  &
                AWG2 * UE(IS+ISP1,ID-IDP ) +                  &
		AWG3 * UE(IS+ISP ,ID-IDP1) +                  &
		AWG4 * UE(IS+ISP ,ID-IDP )
       EM2    = AWG5 * UE(IS+ISM1,ID+IDM1) +                  &
                AWG6 * UE(IS+ISM1,ID+IDM ) +                  &
		AWG7 * UE(IS+ISM ,ID+IDM1) +                  &
		AWG8 * UE(IS+ISM ,ID+IDM )
!
!      *** Contribution to interactions                          ***
!      *** CONS is the shallow water factor for the NL interact. ***
!
       FACTOR = CONS * AF11(IS) * E00
!
       SA1A   = E00 * ( EP1*DAL1 + EM1*DAL2 ) * PQUAD(2)               
       SA1B   = SA1A - EP1*EM1*DAL3 * PQUAD(2)                         
       SA2A   = E00 * ( EP2*DAL1 + EM2*DAL2 ) * PQUAD(2)               
       SA2B   = SA2A - EP2*EM2*DAL3 * PQUAD(2)                         
!
       SA1 (IS,ID) = FACTOR * SA1B
       SA2 (IS,ID) = FACTOR * SA2B
!
!       IF(ITEST >= 100 .AND. TESTFL)THEN
!         WRITE(PRINTF,9002) E00,EP1,EM1,EP2,EM2
!9002     FORMAT (' E00 EP1 EM1 EP2 EM2  :',5E11.4)
!         WRITE(PRINTF,9003) SA1A,SA1B,SA2A,SA2B
!9003     FORMAT (' SA1A SA1B SA2A SA2B  :',4E11.4)
!         WRITE(PRINTF,9004) IS,ID,SA1(IS,ID),SA2(IS,ID)
!9004     FORMAT (' IS ID SA1() SA2()    :',2I4,2E12.4)
!         WRITE(PRINTF,9005) FACTOR
!9005     FORMAT (' FACTOR               : ',E12.4)
!       END IF
!
       DA1C(IS,ID) = CONS * AF11(IS) * ( SA1A + SA1B )
       DA1P(IS,ID) = FACTOR * ( DAL1*E00 - DAL3*EM1 ) * PQUAD(2)       
       DA1M(IS,ID) = FACTOR * ( DAL2*E00 - DAL3*EP1 ) * PQUAD(2)       
!
       DA2C(IS,ID) = CONS * AF11(IS) * ( SA2A + SA2B )
       DA2P(IS,ID) = FACTOR * ( DAL1*E00 - DAL3*EM2 ) * PQUAD(2)       
       DA2M(IS,ID) = FACTOR * ( DAL2*E00 - DAL3*EP2 ) * PQUAD(2)       
     ENDDO
   ENDDO
!
!  *** Fold interactions to side angles if spectral domain ***
!  *** is periodic in directional space                    ***
!
   IF(PERCIR)THEN
     DO ID = 1, IDHGH - MDC
       ID0   = 1 - ID
       DO IS = ISCLW, ISCHG
         SA1 (IS,MDC+ID) = SA1 (IS,  ID   )
         SA2 (IS,MDC+ID) = SA2 (IS,  ID   )
         DA1C(IS,MDC+ID) = DA1C(IS,  ID   )
         DA1P(IS,MDC+ID) = DA1P(IS,  ID   )
         DA1M(IS,MDC+ID) = DA1M(IS,  ID   )
         DA2C(IS,MDC+ID) = DA2C(IS,  ID   )
         DA2P(IS,MDC+ID) = DA2P(IS,  ID   )
         DA2M(IS,MDC+ID) = DA2M(IS,  ID   )
!
         SA1 (IS,  ID0 ) = SA1 (IS, MDC+ID0)
         SA2 (IS,  ID0 ) = SA2 (IS, MDC+ID0)
         DA1C(IS,  ID0 ) = DA1C(IS, MDC+ID0)
         DA1P(IS,  ID0 ) = DA1P(IS, MDC+ID0)
         DA1M(IS,  ID0 ) = DA1M(IS, MDC+ID0)
         DA2C(IS,  ID0 ) = DA2C(IS, MDC+ID0)
         DA2P(IS,  ID0 ) = DA2P(IS, MDC+ID0)
         DA2M(IS,  ID0 ) = DA2M(IS, MDC+ID0)
       ENDDO
     ENDDO
   ENDIF
!
!  *** Put source term together (To save space I=IS and J=ID ***
!  *** is used)                                              ***
!
   PI3   = (2. * PI_W)**3
   DO I = 1, ISSTOP
     SIGPI = SPCSIG(I) * JACOBI                                        
     DO J = IDCMIN(I), IDCMAX(I)
       ID = MOD ( J - 1 + MDC , MDC ) + 1
       SFNL(I,ID) = - 2. * ( SA1(I,J) + SA2(I,J) )                       &
                    + AWG1 * ( SA1(I-ISP1,J-IDP1) + SA2(I-ISP1,J+IDP1) ) &
		    + AWG2 * ( SA1(I-ISP1,J-IDP ) + SA2(I-ISP1,J+IDP ) ) &
		    + AWG3 * ( SA1(I-ISP ,J-IDP1) + SA2(I-ISP ,J+IDP1) ) &
		    + AWG4 * ( SA1(I-ISP ,J-IDP ) + SA2(I-ISP ,J+IDP ) ) &
		    + AWG5 * ( SA1(I-ISM1,J+IDM1) + SA2(I-ISM1,J-IDM1) ) &
		    + AWG6 * ( SA1(I-ISM1,J+IDM ) + SA2(I-ISM1,J-IDM ) ) &
		    + AWG7 * ( SA1(I-ISM ,J+IDM1) + SA2(I-ISM ,J-IDM1) ) &
		    + AWG8 * ( SA1(I-ISM ,J+IDM ) + SA2(I-ISM ,J-IDM ) )
!
       DSNL(I,ID) = - 2. * ( DA1C(I,J) + DA2C(I,J) )                       &
                    + SWG1 * ( DA1P(I-ISP1,J-IDP1) + DA2P(I-ISP1,J+IDP1) ) &
		    + SWG2 * ( DA1P(I-ISP1,J-IDP ) + DA2P(I-ISP1,J+IDP ) ) &
		    + SWG3 * ( DA1P(I-ISP ,J-IDP1) + DA2P(I-ISP ,J+IDP1) ) &
		    + SWG4 * ( DA1P(I-ISP ,J-IDP ) + DA2P(I-ISP ,J+IDP ) ) &
		    + SWG5 * ( DA1M(I-ISM1,J+IDM1) + DA2M(I-ISM1,J-IDM1) ) &
		    + SWG6 * ( DA1M(I-ISM1,J+IDM ) + DA2M(I-ISM1,J-IDM ) ) &
		    + SWG7 * ( DA1M(I-ISM ,J+IDM1) + DA2M(I-ISM ,J-IDM1) ) &
		    + SWG8 * ( DA1M(I-ISM ,J+IDM ) + DA2M(I-ISM ,J-IDM ) )
!
!      *** store results in IMATDA and IMATRA ***
!
       IF(TESTFL) THEN
         PLNL4S(ID,I,IPTST) = SFNL(I,ID) / SIGPI                       
         PLNL4D(ID,I,IPTST) = -1. * DSNL(I,ID) / PI3                   
       END IF
!
       IMATRA(ID,I) = IMATRA(ID,I) + SFNL(I,ID) / SIGPI
       IMATDA(ID,I) = IMATDA(ID,I) - DSNL(I,ID) / PI3
!
!       IF(ITEST >= 90 .AND. TESTFL) THEN
!         WRITE(PRINTF,9006) I,J,SFNL(I,ID),DSNL(I,ID),SPCSIG(I)                                                    30.72
!9006     FORMAT (' IS ID SFNL DSNL SPCSIG:',2I4,3E12.4)                
!       END IF
!
     ENDDO
   ENDDO
!

   RETURN
   END SUBROUTINE SWSNL1

!=====================================================================
!  (The functions and subroutines below have not been tested yet)
!=====================================================================

   REAL FUNCTION WAVE(D)

!  Transforms nondimensional Sqr(w)d/g into nondimensional kd using the
!  dispersion relation and an iterative method

   IF(D-10.) 2,2,1
 1 XX=D
   GO TO 6
 2 IF(D-1.) 3,4,4
 3 X=SQRT(D)
   GO TO 5
 4 X=D
 5 COTHX=1.0/TANH(X)
   XX=X-(X-D*COTHX)/(1.+D*(COTHX**2-1.))
   E=1.-XX/X
   X=XX
   IF(ABS(E)-0.0005) 6,5,5
 6 WAVE=XX
   RETURN
   END FUNCTION WAVE


   SUBROUTINE CGCMP(G,AK,H,CG)

!  Calculates group velocity Cg based on depth and wavenumber
!  Includes a deep water limit for kd > 10.

!  G     : gravitational acceleration
!  AK    : wave number
!  H     : depth
!  CG    : group velocity
!  AKH   : depth x Wave number (kd)
!  RGK   : square root of gk
!  SECH2 : the square of the secant Hyperbolic of kd
!          = 1 - TANH(kd)**2

!  Calculation of group velocity Cg

   AKH=AK*H
   IF(AKH <= 10.)THEN                                                 

!    Shallow water:

!    Cg = (1/2) (1 + (2 kd) / Sinh (2 kd)) (w / k)
!       = (1/2) (1 + (kd SECH2) / Tanh (kd)) (w / k)
!       = (1/2) (1 + (kd SECH2) / Tanh (kd)) (Sqrt (gk Tanh(kd)) / k)
!       = (1/2) (Sqrt (gk) (Sqrt(Tanh (kd)) + kd SECH2 / Sqrt (Tanh (kd)) / k
!       = (1/2) (Sqrt (k) g (Tanh (kd) + kd SECH2) / (k Sqrt (g Tanh (kd)))
!       = (g Tanh(kd) + gkd SECH2) / (2 Sqrt(gk Tanh(kd))

     SECH2=1.-TANH(AKH)**2
     CG=(G*AKH*SECH2+G*TANH(AKH))/(2.*SQRT(G*AK*TANH(AKH)))

   ELSE                                                                

!    Deep water:

!    Cg = w / (2 k)
!       = g / (2 Sqrt (gk))

     RGK=SQRT(G*AK)
     CG=G/(2.*RGK)

   ENDIF                                                               
!
   RETURN

   END SUBROUTINE CGCMP

   SUBROUTINE PRESET(LMAX,N,IW3L,IW4,N2,G,H,W4,AK4,W,DQ,DQ2,DT,DT2,P)  

!(This subroutine has not been tested yet)

!  T1A   : Angle between theta_1 and theta_a
!  T2A   : Angle between theta_2 and theta_a
!  T34   : Angle between theta_3 and theta_4
!  WA    : wa = w1 + w2 = w3 + w4 (resonance conditions)
!  AKA   : absolute value of ka = k1 + k2 = k3 + k4 (resonance conditions)
!  TA    : Angle theta_a representing vector ka

   REAL :: W(LMAX)                                                     
!
!  ================
   DO 120 IT34=1,N2
!  ================
!
     T34=DT*(IT34-1)
     DT3=DT
     IF(IT34 == 1 .OR. IT34 == N2) DT3=DT2
!
!    ===================
     DO 110 IW3=IW3L,IW4
!    ===================
!
       W3=W(IW3)
       AK3=WAVE(W3**2*H/G)/H
!
       WA=W3+W4
       AKA=SQRT(AK3*AK3+AK4*AK4+2.*AK3*AK4*COS(T34))
       TA=ATAN2(AK3*SIN(T34),AK3*COS(T34)+AK4)
       R=SQRT(G*AKA*TANH(AKA*H/2.))/WA-1./SQRT(2.)
!
       TL=0.
       ACS=AKA/(2.*WAVE(WA**2*H/(4.*G))/H)
       ACS=SIGN(1.,ACS)*MIN(1.,ABS(ACS))                                   
       IF(R < 0.) TL=ACOS(ACS)
       IT1S=NINT(TL/DT)+1
!
!     ===================
       DO 100 IT1A=IT1S,N2
!     ===================
!
         T1A=DT*(IT1A-1)
!
         DT1=DT
         IF(IT1A == 1 .OR. IT1A == N2) DT1=DT2
!
!     -----------------------------------------------------------------
         CALL FINDW1W2(LMAX,IW3,W,G,H,W1,W2,W3,WA,AK1,AK2,AKA,T1A,T2A,IND)   
!     -----------------------------------------------------------------
         IF(IND < 0) GO TO 100
!
         IF(W1 <= W3 .AND. W3 <= W4 .AND. W4 <= W2)THEN
!JQI        CALL KERNEL(N,G,H,W1,W2,W3,W4,AK1,AK2,AK3,AK4,AKA,           &
!JQI	            T1A,T2A,T34,TA,DQ,DT,DT1,DT3,P)                       
         ENDIF
!
! ============
  100 CONTINUE
  110 CONTINUE
  120 CONTINUE
! ============
!
   RETURN
   END SUBROUTINE PRESET

   SUBROUTINE FINDW1W2(LMAX,IW3,W,G,H,W1,W2,W3,WA,AK1,AK2,AKA,T1A,T2A,IND)  
   REAL  :: W(LMAX)                                                    
!
   IND=0
   EPS=0.0005
!
   X1=0.005
   X2=W3
!
!     ---------------------------------------
110 CALL FDW1(G,H,AKA,T1A,WA,X1,X2,X,EPS,M)
!     ---------------------------------------
!
   IF(M.EQ.0) GO TO 999
!
   W1=X
!
   AK1=WAVE(W1**2*H/G)/H
   AK2=SQRT(AKA*AKA+AK1*AK1-2.*AKA*AK1*COS(T1A))
   W2=SQRT(G*AK2*TANH(AK2*H))
   IF(W1.GT.W2) GO TO 999
   T2A=ATAN2(-AK1*SIN(T1A),AKA-AK1*COS(T1A))
   RETURN
!
999 IND=-999
   RETURN
   END SUBROUTINE FINDW1W2

   SUBROUTINE FDW1(G,H,AKA,T1A,WA,X1,X2,X,EPS,M)
!
   M=0
   F1=FUNCW1(G,H,X1,AKA,T1A,WA)
   F2=FUNCW1(G,H,X2,AKA,T1A,WA)
20 IF(F1*F2) 18,19,21
18 M=M+1
   X=X2-((X2-X1)/(F2-F1)*F2)
   IF((ABS(X-X2)/ABS(X))-EPS) 22,23,23
23 IF((ABS(X-X1)/ABS(X))-EPS) 22,24,24
22 RETURN
24 F=FUNCW1(G,H,X,AKA,T1A,WA)
   FM=F*F1
   IF(FM) 25,22,26
25 X2=X
   F2=F
   GO TO 20
26 X1=X
   F1=F
   GO TO 20
19 M=M+1
   IF(F1) 17,16,17
16 X=X1
   GO TO 22
17 X=X2
   GO TO 22
21 M=0
   RETURN
   END SUBROUTINE FDW1

   FUNCTION FUNCW1(G,H,X,AKA,T1A,WA)
!
   AK1=WAVE(X**2*H/G)/H
   AK2=SQRT(AKA*AKA+AK1*AK1-2.*AKA*AK1*COS(T1A))
   FUNCW1=WA-SQRT(G*AK1*TANH(AK1*H))-SQRT(G*AK2*TANH(AK2*H))
!
   RETURN
   END FUNCTION FUNCW1
